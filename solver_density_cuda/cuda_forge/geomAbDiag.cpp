// V0 の評価の経路 (plans/active/architecture-float-state-double-geometry.md §4.2c、段 ①)。説明は geomAbDiag.hpp。
// 既定 off: 環境変数 FORGE_DIAG_GEOMAB_DUMP / FORGE_DIAG_GEOMAB_REF が無ければ、フックは呼ばれず何も確保しない。

#include "cuda_forge/geomAbDiag.hpp"
#include "cuda_forge/scalarTransport_d.cuh"   // ScalarTransportDesc の定義だけ (起動関数はフックの呼び出し側から受け取る)
#include "cuda_forge/cudaWrapper.cuh"

#include <highfive/H5File.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace geomAbDiag {
namespace {

// 粘性の面の流束の成分 (viscousFlux_d の faceFlux の並び)
constexpr int kViscComp = 6;
const char* const kViscCompNames = "0:momentum_x 1:momentum_y 2:momentum_z 3:energy (as added to res_roe) "
                                   "4:heat (conduction after Taw/Qw replacement) 5:work (tau.U_f)";
const char* const kSign = "face flux F is added to res[ic0] and subtracted from res[ic1] (ic0/ic1 = plane_cells[2*ip+0/1]); "
                          "NaN = face not evaluated by the kernel";

// 粘性カーネルが読むセル配列 (面の幾何と座標は /mesh に別に書く)。存在するものだけ記録・上書きする。
// 使わない構成で未初期化のまま渡らない配列 (Tau_Wall 等) も記録されるが、カーネルは同じ条件で nullptr を受けるので読まない。
const char* const kViscInputs[] = {
    "cp", "thermCond", "vis_lam", "vis_turb",
    "ro", "roUx", "roUy", "roUz", "roe", "Ux", "Uy", "Uz", "P", "Ht", "sonic", "T",
    "dUxdx", "dUxdy", "dUxdz", "dUydx", "dUydy", "dUydz", "dUzdx", "dUzdy", "dUzdz",
    "dTdx", "dTdy", "dTdz", "axisym_divU",
    "Tau_Wall", "Qw_Wall", "Taw_Prim_Overlay", "Taw_HTnx", "Taw_HTny", "Taw_HTnz", "Taw_diag",
    "k", "dKdx", "dKdy", "dKdz", "sstF1"};
const char* const kViscRes[] = {"res_ro", "res_roUx", "res_roUy", "res_roUz", "res_roe"};
// 両カーネルが読む面の幾何 (デバイスの最終値 = 軸対称 method 0 では半径の重みを掛けた後)
const char* const kPlaneGeom[] = {"fx", "sx", "sy", "sz", "ss"};

struct Cmp { long long n = 0; double maxAbs = 0.0; };

struct HookRec {
    bool reached = false;
    int nComp = 0;                                         // 面の流束の成分数
    std::vector<std::string> inputNames;                   // 記録 (DUMP) / 上書き (REF) したセル配列
    std::map<std::string, std::vector<flow_float>> inputs;
    std::vector<std::string> resNames;                     // 照合した残差 (と輸送対角)
    std::vector<flow_float*> resPtr;                       // 本番の書き先 (resNames と同順)
    std::map<std::string, std::vector<flow_float>> resBefore, resAfter, resOld, resOldRepeat;
    std::vector<flow_float> fluxOld, fluxOldRepeat, fluxNew, fluxRef;
    Cmp resVsProd, resRepeat, fluxRepeat;
    long long nResEntries = 0, nFacesEvaluated = 0;
    std::string line;
    // スカラーの記述子 (属性に書く)
    std::vector<std::string> phiNames, F1Names;
    std::vector<double> sigma, sigma2, sigmaLam;
    std::vector<int> diffusion;
};

struct State {
    int mode = -1;
    std::string path;       // DUMP: 出力 / REF: 入力 (DUMP の出力)
    std::string outPath;    // REF の出力
    bool armed = false;
    bool geomOverwritten = false;
    long long eMismatch = -1, eNonFinite = -1;
    std::vector<int> planeCells, normalHalo;
    HookRec scalar, viscous;
};
State g;

std::vector<flow_float> d2h(const flow_float* d, size_t n)
{
    std::vector<flow_float> h(n);
    gpuErrchk( cudaMemcpy(h.data(), d, n*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    return h;
}
void h2d(const std::vector<flow_float>& h, flow_float* d)
{
    gpuErrchk( cudaMemcpy(d, h.data(), h.size()*sizeof(flow_float), cudaMemcpyHostToDevice) );
}
flow_float* dalloc(size_t n)
{
    flow_float* p = nullptr;
    gpuErrchk( cudaMalloc((void**)&p, n*sizeof(flow_float)) );
    return p;
}
flow_float* dcopy(const flow_float* src, size_t n)
{
    flow_float* p = dalloc(n);
    gpuErrchk( cudaMemcpy(p, src, n*sizeof(flow_float), cudaMemcpyDeviceToDevice) );
    return p;
}
flow_float* dfill(size_t n, int byte)   // byte = 0 → 0、0xFF → NaN (全ビット 1)
{
    flow_float* p = dalloc(n);
    gpuErrchk( cudaMemset(p, byte, n*sizeof(flow_float)) );
    return p;
}

// ビット単位の比較 (NaN 同士でもビットが違えば数える)。maxAbs は両方有限の要素だけ。
Cmp cmpBits(const std::vector<flow_float>& a, const std::vector<flow_float>& b)
{
    Cmp c;
    if (a.size() != b.size()) { c.n = -1; return c; }
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::memcmp(&a[i], &b[i], sizeof(flow_float)) == 0) continue;
        ++c.n;
        const double d = std::fabs((double)a[i] - (double)b[i]);
        if (std::isfinite(d) && d > c.maxAbs) c.maxAbs = d;
    }
    return c;
}
void addCmp(Cmp& acc, const Cmp& c)
{
    if (c.n < 0 || acc.n < 0) { acc.n = -1; return; }
    acc.n += c.n;
    if (c.maxAbs > acc.maxAbs) acc.maxAbs = c.maxAbs;
}

// var.c_d の中で同じデバイスポインタを持つ名前 (無ければ空)
std::string nameOf(variables& var, const flow_float* p)
{
    if (p == nullptr) return std::string();
    for (const auto& kv : var.c_d) if (kv.second == p) return kv.first;
    return std::string();
}

std::string refOutPath(const std::string& in)
{
    const std::string ext = ".h5";
    if (in.size() > ext.size() && in.compare(in.size() - ext.size(), ext.size(), ext) == 0)
        return in.substr(0, in.size() - ext.size()) + "_ref.h5";
    return in + "_ref.h5";
}

[[noreturn]] void refuse(const std::string& why)
{
    fprintf(stderr, "[geomab] refused: %s\n", why.c_str());
    fflush(stdout);
    std::exit(EXIT_FAILURE);
}

// REF: DUMP の出力から 1 本読み、大きさを確かめてデバイスへ上書きする
void overwriteFromDump(HighFive::File& f, const std::string& dsName, flow_float* d, size_t n)
{
    if (!f.exist(dsName)) refuse("dump has no dataset " + dsName);
    std::vector<flow_float> h;
    f.getDataSet(dsName).read(h);   // float32 → double は丸めなしで広がる
    if (h.size() != n) refuse("dataset " + dsName + " has " + std::to_string(h.size()) + " entries, expected " + std::to_string(n));
    h2d(h, d);
}

// REF: 面の幾何 (fx・sx..ss) を DUMP の値 (float を double に広げたもの) で 1 回だけ上書きする。
// FP64 のビルド自身の面ベクトルとは入れ替えない (§4.2c 2.)。
void overwriteGeomOnce(HighFive::File& f, mesh& msh, variables& var)
{
    if (g.geomOverwritten) return;
    for (const char* k : kPlaneGeom) overwriteFromDump(f, std::string("/mesh/") + k, var.p_d.at(k), (size_t)msh.nPlanes);
    g.geomOverwritten = true;
    printf("[geomab] REF: plane geometry fx, sx, sy, sz, ss overwritten with the dumped float values (widened)\n");
}

long long countFinite(const std::vector<flow_float>& v, size_t n)
{
    long long c = 0;
    for (size_t i = 0; i < n && i < v.size(); ++i) if (std::isfinite((double)v[i])) ++c;
    return c;
}

void writeVec(HighFive::File& h5, const std::string& name, const std::vector<flow_float>& v)
{
    h5.createDataSet(name, v);
}
void writeMat(HighFive::File& h5, const std::string& name, const std::vector<flow_float>& v, size_t rows, size_t cols)
{
    auto ds = h5.createDataSet<flow_float>(name, HighFive::DataSpace(std::vector<size_t>{rows, cols}));
    ds.write_raw(v.data());
}

void writeHook(HighFive::File& h5, const std::string& grp, const HookRec& r, size_t nP, bool dump)
{
    if (!r.reached) return;
    if (dump) {
        for (const auto& nm : r.inputNames) writeVec(h5, "/" + grp + "/inputs/" + nm, r.inputs.at(nm));
        writeMat(h5, "/" + grp + "/old/face_flux", r.fluxOld, (size_t)r.nComp, nP);
        writeMat(h5, "/" + grp + "/new/face_flux", r.fluxNew, (size_t)r.nComp, nP);
        for (const auto& nm : r.resNames) {
            writeVec(h5, "/" + grp + "/check/" + nm + "/before_production", r.resBefore.at(nm));
            writeVec(h5, "/" + grp + "/check/" + nm + "/after_production", r.resAfter.at(nm));
            writeVec(h5, "/" + grp + "/check/" + nm + "/old_arm_from_before", r.resOld.at(nm));
        }
        h5.createAttribute(grp + "_input_names", r.inputNames);
        h5.createAttribute(grp + "_check_names", r.resNames);
        h5.createAttribute(grp + "_old_vs_production_mismatch", r.resVsProd.n);
        h5.createAttribute(grp + "_old_vs_production_maxabs", r.resVsProd.maxAbs);
        h5.createAttribute(grp + "_old_repeat_residual_mismatch", r.resRepeat.n);
        h5.createAttribute(grp + "_old_repeat_face_flux_mismatch", r.fluxRepeat.n);
        h5.createAttribute(grp + "_check_entries", r.nResEntries);
        h5.createAttribute(grp + "_check_line", r.line);
    } else {
        writeMat(h5, "/" + grp + "/ref/face_flux", r.fluxRef, (size_t)r.nComp, nP);
        h5.createAttribute(grp + "_overwritten_inputs", r.inputNames);
    }
    h5.createAttribute(grp + "_faces_evaluated", r.nFacesEvaluated);
    h5.createAttribute(grp + "_n_components", r.nComp);
}

}  // namespace

int mode()
{
    if (g.mode >= 0) return g.mode;
    const char* d = std::getenv("FORGE_DIAG_GEOMAB_DUMP");
    const char* r = std::getenv("FORGE_DIAG_GEOMAB_REF");
    const bool hd = (d != nullptr && *d != '\0');
    const bool hr = (r != nullptr && *r != '\0');
    if (hd && hr) refuse("FORGE_DIAG_GEOMAB_DUMP and FORGE_DIAG_GEOMAB_REF are both set; set only one");
    g.mode = hd ? 1 : (hr ? 2 : 0);
    if (hd) g.path = d;
    if (hr) { g.path = r; g.outPath = refOutPath(g.path); }
    return g.mode;
}

bool armed() { return g.armed; }
void arm(bool on) { g.armed = on; }

bool begin(solverConfig& cfg, mesh& msh, variables& var)
{
    const int m = mode();
    const size_t nP = (size_t)msh.nPlanes, nCa = (size_t)msh.nCells_all;
    if (m == 1) {
        printf("[geomab] FORGE_DIAG_GEOMAB_DUMP=%s: one assembly (= step 1), viscous + k/omega diffusion evaluated with the old "
               "difference (old arm) and with e32 (new arm) into scratch buffers, then exit without update (flow_float %zu bytes)\n",
               g.path.c_str(), sizeof(flow_float));
        if (sizeof(flow_float) != 4) printf("[geomab] WARNING: DUMP is meant for the float build (this build has flow_float = %zu bytes)\n", sizeof(flow_float));
    } else if (m == 2) {
        printf("[geomab] FORGE_DIAG_GEOMAB_REF=%s: one assembly (= step 1), inputs overwritten from the dump, kernels evaluated once "
               "with this build's own e, written to %s, then exit without update (flow_float %zu bytes)\n",
               g.path.c_str(), g.outPath.c_str(), sizeof(flow_float));
        if (sizeof(flow_float) != 8) printf("[geomab] WARNING: REF is meant for the FP64 build (this build has flow_float = %zu bytes)\n", sizeof(flow_float));
    } else {
        return false;
    }

    // 接続 (並びは plane_cells[2*ip+k] = planes[ip].iCells[k]、setMeshMap_d と同じ)
    g.planeCells.assign(2*nP, -1);
    for (size_t ip = 0; ip < nP; ++ip) {
        const auto& pc = msh.planes[ip].iCells;
        for (size_t k = 0; k < 2 && k < pc.size(); ++k) g.planeCells[2*ip + k] = pc[k];
    }
    g.normalHalo.assign((size_t)msh.nNormal_halo_Planes, -1);
    if (msh.nNormal_halo_Planes > 0) {
        gpuErrchk( cudaMemcpy(g.normalHalo.data(), msh.normal_halo_planes_d, g.normalHalo.size()*sizeof(geom_int), cudaMemcpyDeviceToHost) );
    }

    // e32 と、今のカーネルが作る差 fl(ccx[ic1] − ccx[ic0]) (デバイスの ccx を flow_float のまま引く) の照合。
    // FP64 のビルドでは座標が double のままなのでビット一致が期待値 (plan §5.1 #3 の合格条件)。
    {
        const auto cx = d2h(var.c_d.at("ccx"), nCa), cy = d2h(var.c_d.at("ccy"), nCa), cz = d2h(var.c_d.at("ccz"), nCa);
        const auto& ex = var.p.at("ge_x"); const auto& ey = var.p.at("ge_y"); const auto& ez = var.p.at("ge_z");
        long long nMis = 0, nNonFin = 0;
        for (size_t ip = 0; ip < nP; ++ip) {
            const int i0 = g.planeCells[2*ip + 0], i1 = g.planeCells[2*ip + 1];
            if (i0 < 0 || i1 < 0) continue;
            const flow_float dx = cx[i1] - cx[i0], dy = cy[i1] - cy[i0], dz = cz[i1] - cz[i0];
            if (!std::isfinite((double)ex[ip]) || !std::isfinite((double)ey[ip]) || !std::isfinite((double)ez[ip])) ++nNonFin;
            if (std::memcmp(&dx, &ex[ip], sizeof(flow_float)) != 0 || std::memcmp(&dy, &ey[ip], sizeof(flow_float)) != 0
                || std::memcmp(&dz, &ez[ip], sizeof(flow_float)) != 0) ++nMis;
        }
        g.eMismatch = nMis; g.eNonFinite = nNonFin;
        printf("[geomab] e: %zu planes; e32 differs bitwise from the current difference fl(cc[ic1]) - fl(cc[ic0]) in %lld planes "
               "(FP64 build: expected 0); non-finite e32 in %lld planes\n", nP, nMis, nNonFin);
    }

    if (m == 2) {
        // DUMP の出力と大きさ・接続が一致することを確かめる (違えば面の番号で比べられない)
        HighFive::File f(g.path, HighFive::File::ReadOnly);
        auto chk = [&](const char* a, long long want) {
            if (!f.hasAttribute(a)) refuse(std::string("dump has no attribute ") + a);
            const long long got = f.getAttribute(a).read<long long>();
            if (got != want) refuse(std::string(a) + " differs: dump " + std::to_string(got) + ", this run " + std::to_string(want));
        };
        chk("nCells", msh.nCells); chk("nCells_all", msh.nCells_all); chk("nPlanes", msh.nPlanes);
        chk("nNormalPlanes", msh.nNormalPlanes); chk("nNormal_halo_Planes", msh.nNormal_halo_Planes);
        std::vector<int> pc, nh;
        f.getDataSet("/mesh/plane_cells").read(pc);
        f.getDataSet("/mesh/normal_halo_planes").read(nh);
        if (pc != g.planeCells) refuse("plane_cells differ from the dump (connectivity must be identical)");
        if (nh != g.normalHalo) refuse("normal_halo_planes differ from the dump");
        printf("[geomab] REF: sizes and connectivity match the dump\n");
    }
    (void)cfg;
    return true;
}

void scalarBefore(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                  const ScalarTransportDesc* descs, int n, ScalarArmFn armFn)
{
    HookRec& r = g.scalar;
    if (r.reached) return;   // 組立 1 回の中で 1 度だけ
    const size_t nCa = (size_t)msh.nCells_all, nP = (size_t)msh.nPlanes;
    r.nComp = n;

    // カーネルが読むセル配列 (φ・F1・ρ・μ・μ_t) と、書き先 (残差・輸送対角)
    std::vector<const flow_float*> inPtr;
    auto addIn = [&](const flow_float* p, const std::string& fallback) {
        if (p == nullptr) return;
        for (auto q : inPtr) if (q == p) return;
        std::string nm = nameOf(var, p);
        if (nm.empty()) nm = fallback;
        inPtr.push_back(p); r.inputNames.push_back(nm);
    };
    for (int s = 0; s < n; ++s) {
        addIn(descs[s].phi, "phi" + std::to_string(s));
        addIn(descs[s].F1, "F1_" + std::to_string(s));
        r.phiNames.push_back(nameOf(var, descs[s].phi));
        r.F1Names.push_back(nameOf(var, descs[s].F1));
        r.sigma.push_back(descs[s].sigma); r.sigma2.push_back(descs[s].sigma2); r.sigmaLam.push_back(descs[s].sigma_lam);
        r.diffusion.push_back((cfg.scalarDiffusion == 1 && descs[s].diffusion == 1) ? 1 : 0);
    }
    addIn(var.c_d.at("ro"), "ro"); addIn(var.c_d.at("vis_lam"), "vis_lam"); addIn(var.c_d.at("vis_turb"), "vis_turb");
    for (int s = 0; s < n; ++s) {
        r.resNames.push_back(nameOf(var, descs[s].res_rho_phi)); r.resPtr.push_back(descs[s].res_rho_phi);
        r.resNames.push_back(nameOf(var, descs[s].transport_diag)); r.resPtr.push_back(descs[s].transport_diag);
    }

    std::vector<flow_float*> resAlt(n), diagAlt(n);
    auto runArm = [&](bool fromProduction, const flow_float* gex, const flow_float* gey, const flow_float* gez,
                      std::vector<flow_float>& fluxOut, std::map<std::string, std::vector<flow_float>>* resOut) {
        for (int s = 0; s < n; ++s) {
            resAlt[s]  = fromProduction ? dcopy(descs[s].res_rho_phi, nCa)   : dfill(nCa, 0);
            diagAlt[s] = fromProduction ? dcopy(descs[s].transport_diag, nCa) : dfill(nCa, 0);
        }
        flow_float* ff = dfill((size_t)n*nP, 0xFF);
        armFn(cfg, cuda_cfg, msh, var, descs, n, resAlt.data(), diagAlt.data(), gex, gey, gez, ff);
        fluxOut = d2h(ff, (size_t)n*nP);
        if (resOut != nullptr) {
            for (int s = 0; s < n; ++s) {
                (*resOut)[r.resNames[2*s + 0]] = d2h(resAlt[s], nCa);
                (*resOut)[r.resNames[2*s + 1]] = d2h(diagAlt[s], nCa);
            }
        }
        for (int s = 0; s < n; ++s) { cudaFree(resAlt[s]); cudaFree(diagAlt[s]); }
        cudaFree(ff);
    };

    if (g.mode == 1) {
        for (size_t i = 0; i < inPtr.size(); ++i) r.inputs[r.inputNames[i]] = d2h(inPtr[i], nCa);
        for (size_t i = 0; i < r.resPtr.size(); ++i) r.resBefore[r.resNames[i]] = d2h(r.resPtr[i], nCa);
        // 旧腕 (今の座標の差)、同じ条件の旧腕 (atomicAdd の順序の目安)、新腕 (e32)
        runArm(true, nullptr, nullptr, nullptr, r.fluxOld, &r.resOld);
        runArm(true, nullptr, nullptr, nullptr, r.fluxOldRepeat, &r.resOldRepeat);
        runArm(false, var.p_d.at("ge_x"), var.p_d.at("ge_y"), var.p_d.at("ge_z"), r.fluxNew, nullptr);
        r.nFacesEvaluated = countFinite(r.fluxOld, nP);
        printf("[geomab] scalar diffusion (%d scalars): inputs recorded, old/new arms evaluated on %lld faces\n", n, r.nFacesEvaluated);
    } else {
        HighFive::File f(g.path, HighFive::File::ReadOnly);
        for (size_t i = 0; i < inPtr.size(); ++i)
            overwriteFromDump(f, "/scalar/inputs/" + r.inputNames[i], const_cast<flow_float*>(inPtr[i]), nCa);
        overwriteGeomOnce(f, msh, var);
        runArm(false, var.p_d.at("ge_x"), var.p_d.at("ge_y"), var.p_d.at("ge_z"), r.fluxRef, nullptr);
        r.nFacesEvaluated = countFinite(r.fluxRef, nP);
        printf("[geomab] scalar diffusion (%d scalars): %zu input arrays overwritten from the dump, reference arm evaluated on %lld faces\n",
               n, inPtr.size(), r.nFacesEvaluated);
    }
    r.reached = true;
}

void scalarAfter(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                 const ScalarTransportDesc* descs, int n)
{
    (void)cfg; (void)cuda_cfg; (void)var; (void)descs; (void)n;
    HookRec& r = g.scalar;
    if (g.mode != 1 || !r.reached || !r.resAfter.empty()) return;
    const size_t nCa = (size_t)msh.nCells_all;
    for (size_t i = 0; i < r.resPtr.size(); ++i) {
        const std::string& nm = r.resNames[i];
        r.resAfter[nm] = d2h(r.resPtr[i], nCa);
        addCmp(r.resVsProd, cmpBits(r.resOld.at(nm), r.resAfter.at(nm)));
        addCmp(r.resRepeat, cmpBits(r.resOld.at(nm), r.resOldRepeat.at(nm)));
        r.nResEntries += (long long)nCa;
    }
    r.fluxRepeat = cmpBits(r.fluxOld, r.fluxOldRepeat);
    char buf[512];
    snprintf(buf, sizeof(buf), "[geomab] scalar diffusion check: old arm (from the pre-call residual) vs production, bitwise mismatch %lld of %lld "
             "entries (residual + transport_diag, max|d| %.3e); old-arm repeat: residual mismatch %lld, face-flux mismatch %lld",
             r.resVsProd.n, r.nResEntries, r.resVsProd.maxAbs, r.resRepeat.n, r.fluxRepeat.n);
    r.line = buf;
    printf("%s\n", r.line.c_str());
}

void viscousBefore(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, ViscousArmFn armFn)
{
    HookRec& r = g.viscous;
    if (r.reached) return;
    const size_t nCa = (size_t)msh.nCells_all, nP = (size_t)msh.nPlanes;
    r.nComp = kViscComp;

    std::vector<flow_float*> inPtr;
    for (const char* k : kViscInputs) {
        auto it = var.c_d.find(k);
        if (it == var.c_d.end() || it->second == nullptr) continue;
        inPtr.push_back(it->second); r.inputNames.push_back(k);
    }
    for (const char* k : kViscRes) { r.resNames.push_back(k); r.resPtr.push_back(var.c_d.at(k)); }

    auto runArm = [&](bool fromProduction, const flow_float* gex, const flow_float* gey, const flow_float* gez,
                      std::vector<flow_float>& fluxOut, std::map<std::string, std::vector<flow_float>>* resOut) {
        flow_float* res[5];
        for (int c = 0; c < 5; ++c) res[c] = fromProduction ? dcopy(r.resPtr[c], nCa) : dfill(nCa, 0);
        flow_float* ff = dfill((size_t)kViscComp*nP, 0xFF);
        armFn(cfg, cuda_cfg, msh, var, res, gex, gey, gez, ff);
        fluxOut = d2h(ff, (size_t)kViscComp*nP);
        if (resOut != nullptr) for (int c = 0; c < 5; ++c) (*resOut)[r.resNames[c]] = d2h(res[c], nCa);
        for (int c = 0; c < 5; ++c) cudaFree(res[c]);
        cudaFree(ff);
    };

    if (g.mode == 1) {
        for (size_t i = 0; i < inPtr.size(); ++i) r.inputs[r.inputNames[i]] = d2h(inPtr[i], nCa);
        for (size_t i = 0; i < r.resPtr.size(); ++i) r.resBefore[r.resNames[i]] = d2h(r.resPtr[i], nCa);
        runArm(true, nullptr, nullptr, nullptr, r.fluxOld, &r.resOld);
        runArm(true, nullptr, nullptr, nullptr, r.fluxOldRepeat, &r.resOldRepeat);
        runArm(false, var.p_d.at("ge_x"), var.p_d.at("ge_y"), var.p_d.at("ge_z"), r.fluxNew, nullptr);
        r.nFacesEvaluated = countFinite(r.fluxOld, nP);
        printf("[geomab] viscous flux: %zu input arrays recorded, old/new arms evaluated on %lld faces\n", inPtr.size(), r.nFacesEvaluated);
    } else {
        HighFive::File f(g.path, HighFive::File::ReadOnly);
        for (size_t i = 0; i < inPtr.size(); ++i)
            overwriteFromDump(f, "/viscous/inputs/" + r.inputNames[i], inPtr[i], nCa);
        overwriteGeomOnce(f, msh, var);
        runArm(false, var.p_d.at("ge_x"), var.p_d.at("ge_y"), var.p_d.at("ge_z"), r.fluxRef, nullptr);
        r.nFacesEvaluated = countFinite(r.fluxRef, nP);
        printf("[geomab] viscous flux: %zu input arrays overwritten from the dump, reference arm evaluated on %lld faces\n",
               inPtr.size(), r.nFacesEvaluated);
    }
    r.reached = true;
}

void viscousAfter(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg; (void)var;
    HookRec& r = g.viscous;
    if (g.mode != 1 || !r.reached || !r.resAfter.empty()) return;
    const size_t nCa = (size_t)msh.nCells_all;
    for (size_t i = 0; i < r.resPtr.size(); ++i) {
        const std::string& nm = r.resNames[i];
        r.resAfter[nm] = d2h(r.resPtr[i], nCa);
        addCmp(r.resVsProd, cmpBits(r.resOld.at(nm), r.resAfter.at(nm)));
        addCmp(r.resRepeat, cmpBits(r.resOld.at(nm), r.resOldRepeat.at(nm)));
        r.nResEntries += (long long)nCa;
    }
    r.fluxRepeat = cmpBits(r.fluxOld, r.fluxOldRepeat);
    char buf[512];
    snprintf(buf, sizeof(buf), "[geomab] viscous check: old arm (from the pre-call residual) vs production, bitwise mismatch %lld of %lld "
             "entries (res_ro..res_roe, max|d| %.3e); old-arm repeat: residual mismatch %lld, face-flux mismatch %lld",
             r.resVsProd.n, r.nResEntries, r.resVsProd.maxAbs, r.resRepeat.n, r.fluxRepeat.n);
    r.line = buf;
    printf("%s\n", r.line.c_str());
}

int finish(solverConfig& cfg, mesh& msh, variables& var)
{
    const size_t nP = (size_t)msh.nPlanes, nCa = (size_t)msh.nCells_all;
    if (!g.viscous.reached) fprintf(stderr, "[geomab] WARNING: the viscous hook was not reached in the assembly\n");
    if (!g.scalar.reached)  printf("[geomab] note: the k/omega diffusion hook was not reached (not SST, or scalarDiffusion off)\n");
    if (!g.viscous.reached && !g.scalar.reached) refuse("neither hook was reached; nothing to write");

    const bool dump = (g.mode == 1);
    const std::string out = dump ? g.path : g.outPath;
    {
        HighFive::File h5(out, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
        h5.createAttribute("geomab_mode", std::string(dump ? "dump" : "ref"));
        h5.createAttribute("plan", std::string("plans/active/architecture-float-state-double-geometry.md §4.2c (V0)"));
        h5.createAttribute("flow_float_bytes", (long long)sizeof(flow_float));
        h5.createAttribute("nCells", (long long)msh.nCells);
        h5.createAttribute("nCells_all", (long long)msh.nCells_all);
        h5.createAttribute("nPlanes", (long long)msh.nPlanes);
        h5.createAttribute("nNormalPlanes", (long long)msh.nNormalPlanes);
        h5.createAttribute("nNormal_halo_Planes", (long long)msh.nNormal_halo_Planes);
        h5.createAttribute("discretization", cfg.discretization);
        h5.createAttribute("isAxisymmetric", (long long)cfg.isAxisymmetric);
        h5.createAttribute("axisymMethod", (long long)cfg.axisymMethod);
        h5.createAttribute("mesh_file", cfg.meshFileName);
        h5.createAttribute("input_value_file", cfg.valueFileName);
        h5.createAttribute("e_mismatch_vs_current_difference", g.eMismatch);
        h5.createAttribute("e_nonfinite", g.eNonFinite);
        h5.createAttribute("sign_convention", std::string(kSign));
        h5.createAttribute("viscous_components", std::string(kViscCompNames));
        h5.createAttribute("layout", std::string("face_flux [n_components, nPlanes]; plane_cells [2*ip+k]; cell arrays [nCells_all] incl. ghosts"));
        if (!dump) h5.createAttribute("dump_path", g.path);

        // 接続・座標・e (DUMP は float の座標と e32 に加えて double の e64、REF はこのビルド自身の座標と e)
        h5.createDataSet("/mesh/plane_cells", g.planeCells);
        h5.createDataSet("/mesh/normal_halo_planes", g.normalHalo);
        for (const char* k : {"ccx", "ccy", "ccz"}) writeVec(h5, std::string("/mesh/") + k, d2h(var.c_d.at(k), nCa));
        for (const char* k : {"ge_x", "ge_y", "ge_z"}) writeVec(h5, std::string("/mesh/") + k, var.p.at(k));
        if (dump) {
            // 面の幾何はカーネルが読んだデバイスの値 (REF はこれで上書きする)
            for (const char* k : kPlaneGeom) writeVec(h5, std::string("/mesh/") + k, d2h(var.p_d.at(k), nP));
            if (msh.cc64.size() == 3*nCa) {
                std::vector<double> e64[3];
                for (int k = 0; k < 3; ++k) e64[k].assign(nP, std::numeric_limits<double>::quiet_NaN());
                for (size_t ip = 0; ip < nP; ++ip) {
                    const int i0 = g.planeCells[2*ip + 0], i1 = g.planeCells[2*ip + 1];
                    if (i0 < 0 || i1 < 0) continue;
                    for (int k = 0; k < 3; ++k) e64[k][ip] = msh.cc64[3*(size_t)i1 + k] - msh.cc64[3*(size_t)i0 + k];
                }
                h5.createDataSet("/mesh/ge64_x", e64[0]);
                h5.createDataSet("/mesh/ge64_y", e64[1]);
                h5.createDataSet("/mesh/ge64_z", e64[2]);
            }
        }

        if (g.scalar.reached) {
            h5.createAttribute("scalar_phi_names", g.scalar.phiNames);
            h5.createAttribute("scalar_F1_names", g.scalar.F1Names);
            h5.createAttribute("scalar_sigma", g.scalar.sigma);
            h5.createAttribute("scalar_sigma2", g.scalar.sigma2);
            h5.createAttribute("scalar_sigma_lam", g.scalar.sigmaLam);
            h5.createAttribute("scalar_diffusion", g.scalar.diffusion);
        }
        writeHook(h5, "scalar", g.scalar, nP, dump);
        writeHook(h5, "viscous", g.viscous, nP, dump);
    }
    printf("[geomab] wrote %s; exiting without any update\n", out.c_str());
    fflush(stdout);
    return 0;
}

}  // namespace geomAbDiag
