#include <cstdlib>
#include <limits>
#include <set>
#include <map>
#include <cstring>
#include <algorithm>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include "../condensationTransport_d.cuh"   // cond_prop_opts() (凝縮物性オプション → CondArgs)
#include "../condensationEOS_d.cuh"         // cond_face_h_cpg (SLAU CPG 二相の面エンタルピー)
#include "../condensationSourceF_d.cuh"     // cond_face_h_cpg_f / cond_tab_latent_f (float 面潜熱, plan condensation-float-speedup)
#include "convectiveFlux_d.cuh"
#include "lowMachPrecond_d.cuh"
#include "speciesTransport_d.cuh"  // species_roY_device_ptr()
#include "passiveTransport_d.cuh"  // 受動種 S3 (passiveScalarScheme 1)
#include "condensationProperties_d.cuh"  // n2_latent (二相エネルギー流束の潜熱補正)
#include "convectiveFlux_common_d.cuh"

#include <stdexcept>
#include "convectiveFlux_slau_d.inc.cuh"
#include "farfieldFlux_d.inc.cuh"          // 遠方境界 farfield (node、plan boundary-node-farfield-characteristic)

#include "legacy/convectiveFlux_ausm_keep_d.inc.cuh"



#include "convectiveFlux_hlle_d.inc.cuh"



#include "convectiveFlux_roe_d.inc.cuh"

#include "convectiveFlux_keep_d.inc.cuh"

#include "legacy/convectiveFlux_keepslau_d.inc.cuh"

#include "convectiveFlux_boundary_d.inc.cuh"




// ---- 遠方境界 farfield の面の値配列 (plan boundary-node-farfield-characteristic §4.3) ----
// 化学種・k・ω の輸送カーネルが、farfield の境界半割面で流入 (質量流束 < 0) のときに運ぶ外側状態の値。
// 面 (plane) ごとの配列で、farfield 以外の面は NaN (= 使わない。既存の境界はビット不変)。
namespace {
struct FfFaceArrays {
    bool init = false;
    std::vector<flow_float*> Y;       // 化学種ごと
    flow_float** Ydev = nullptr;      // デバイス上のポインタ配列 (S3 カーネル用)
    flow_float* k = nullptr;
    flow_float* om = nullptr;
    geom_int nPlanes = 0;
};
FfFaceArrays& ffFace() { static FfFaceArrays a; return a; }

void ffFaceAlloc(solverConfig& cfg, mesh& msh)
{
    FfFaceArrays& A = ffFace();
    if (A.init) return;
    A.init = true;
    bool any = false;
    for (auto& bc : msh.bconds) any = any || (bc.bcondKind == "farfield");
    if (!any) return;
    // 周期境界と節点を共有する farfield は初版の対象外 (plan §2)。起動時に拒否する
    {
        std::vector<char> per((size_t)msh.nCells + 1, 0);
        for (auto& bc : msh.bconds) if (bc.bcondKind == "periodic") for (auto c : bc.iCells) if (c >= 0 && (size_t)c < per.size()) per[c] = 1;
        for (auto& bc : msh.bconds) {
            if (bc.bcondKind != "farfield") continue;
            for (auto c : bc.iCells) if (c >= 0 && (size_t)c < per.size() && per[c]) {
                std::cerr << "Error: farfield (physID " << bc.physID << ") の節点 " << c << " が周期境界と共有されている (plan boundary-node-farfield-characteristic §2: 初版は非対応)" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }
    A.nPlanes = msh.nPlanes;
    std::vector<flow_float> nanv((size_t)msh.nPlanes, std::numeric_limits<flow_float>::quiet_NaN());
    auto mk = [&]() {
        flow_float* d = nullptr;
        CHECK_CUDA_ERROR(cudaMalloc(&d, sizeof(flow_float) * (size_t)msh.nPlanes));
        CHECK_CUDA_ERROR(cudaMemcpy(d, nanv.data(), sizeof(flow_float) * (size_t)msh.nPlanes, cudaMemcpyHostToDevice));
        return d;
    };
    const int nY = (cfg.thermalMethod == 2 && cfg.nSpecies >= 2) ? cfg.nSpecies : 0;
    for (int sp = 0; sp < nY; ++sp) A.Y.push_back(mk());
    if (nY > 0) {
        CHECK_CUDA_ERROR(cudaMalloc(&A.Ydev, sizeof(flow_float*) * (size_t)nY));
        CHECK_CUDA_ERROR(cudaMemcpy(A.Ydev, A.Y.data(), sizeof(flow_float*) * (size_t)nY, cudaMemcpyHostToDevice));
    }
    if (cfg.LESorRANS == 2 && cfg.RANSmodel == 1) { A.k = mk(); A.om = mk(); }
    std::cout << "[farfield] 面の値配列を確保 (化学種 " << nY << "、k/ω " << (A.k ? "あり" : "なし") << "、" << msh.nPlanes << " 面)" << std::endl;
}
}  // namespace

flow_float* farfieldFaceScalar(const std::string& name)
{
    FfFaceArrays& A = ffFace();
    if (!A.init) return nullptr;
    if (name == "k") return A.k;
    if (name == "omega") return A.om;
    if (name.size() > 1 && name[0] == 'Y') {
        const int sp = std::atoi(name.c_str() + 1);
        if (sp >= 0 && sp < (int)A.Y.size()) return A.Y[sp];
    }
    return nullptr;
}

flow_float** farfieldFaceYDevice() { return ffFace().Ydev; }

static void farfieldFlux_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, bcond& bc, mesh& msh, variables& var, int sstEnergyK)
{
    ffFaceAlloc(cfg, msh);
    FfFaceArrays& A = ffFace();
    const int nY = (cfg.thermalMethod == 2 && cfg.nSpecies >= 2) ? cfg.nSpecies : 0;
    FfGas gas{cfg.thermalMethod, (double)cfg.gamma, (double)cfg.cp, thermo_species_device_ptr(), cfg.nSpecies};
    FfInf inf{};
    auto hv = [&](const char* nm, double dflt) -> double {
        auto it = bc.bvar.find(nm);
        return (it != bc.bvar.end() && !it->second.empty()) ? (double)it->second[0] : dflt;
    };
    inf.r = hv("ro", 0.0); inf.u[0] = hv("Ux", 0.0); inf.u[1] = hv("Uy", 0.0); inf.u[2] = hv("Uz", 0.0); inf.p = hv("Ps", 0.0);
    inf.k = hv("k", 0.0); inf.om = hv("omega", 0.0);
    for (int sp = 0; sp < nY; ++sp) inf.Y[sp] = hv(("Y" + std::to_string(sp)).c_str(), sp == 0 ? 1.0 : 0.0);
    static bool s_logged = false;
    if (!s_logged) {
        s_logged = true;
        std::cout << "[farfield] physID " << bc.physID << ": 自由流 ρ " << inf.r << ", u (" << inf.u[0] << ", " << inf.u[1] << ", " << inf.u[2]
                  << "), P " << inf.p << ", k " << inf.k << ", ω " << inf.om << std::endl;
    }
    const bool rans = (cfg.LESorRANS == 2 && cfg.RANSmodel == 1);
    // 診断ダンプ (env FORGE_DUMP_FARFIELD=<path>、既定 off。出力専用): 最初の呼び出し (FORGE_DUMP_FARFIELD_CALLS=n なら最初の n 回)
    // で面ごとの流束・外側状態を書く。1 回目は <path>.<physID>.csv、2 回目以降は <path>.<physID>.<回>.csv
    static std::map<int, int> s_calls;   // physID ごとの呼び出し回数 (= assembleResidual の回数)
    const int call = ++s_calls[bc.physID];
    static const int s_maxCalls = [] { const char* c = std::getenv("FORGE_DUMP_FARFIELD_CALLS"); return c ? std::max(1, std::atoi(c)) : 1; }();
    float* dumpBuf = nullptr;
    const char* dumpPath = std::getenv("FORGE_DUMP_FARFIELD");
    const size_t nb = bc.iPlanes.size();
    if (dumpPath && *dumpPath && call <= s_maxCalls && nb > 0) {
        CHECK_CUDA_ERROR(cudaMalloc(&dumpBuf, sizeof(float) * nb * FF_DUMP_NF));
    }
    farfield_flux_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(
        (geom_int)nb, bc.map_bplane_plane_d, bc.map_bplane_cell_d,
        var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
        var.c_d["ro"], var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"], var.c_d["P"],
        (nY > 0) ? species_roY_device_ptr() : nullptr,
        rans ? var.c_d["k"] : nullptr, rans ? var.c_d["omega"] : nullptr,
        gas, inf, (flow_float)cfg.pRef, sstEnergyK,
        var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"],
        var.p_d["massflux"], A.Ydev, A.k, A.om, dumpBuf);
    gpuErrchk( cudaPeekAtLastError() );
    if (dumpBuf != nullptr) {
        gpuErrchkKernelSync();
        std::vector<float> h(nb * FF_DUMP_NF);
        CHECK_CUDA_ERROR(cudaMemcpy(h.data(), dumpBuf, sizeof(float) * h.size(), cudaMemcpyDeviceToHost));
        CHECK_CUDA_ERROR(cudaFree(dumpBuf));
        std::ofstream o(std::string(dumpPath) + "." + std::to_string(bc.physID) + (call > 1 ? "." + std::to_string(call) : std::string()) + ".csv");
        o << "ip,ic,nx,ny,nz,S,F_ro,F_roUx,F_roUy,F_roUz,F_roe,R_ro,R_Ux,R_Uy,R_Uz,R_P,R_k,R_om,R_Y0,vacuum,hll,pRef,c_i,c_R\n";
        o.precision(9);
        for (size_t b = 0; b < nb; ++b) {
            for (int q = 0; q < FF_DUMP_NF; ++q) o << (q ? "," : "") << (double)h[b * FF_DUMP_NF + q];
            o << "\n";
        }
        if (call == 1) std::cout << "[FORGE_DUMP_FARFIELD] physID " << bc.physID << ": " << nb << " 面を書いた (最初の " << s_maxCalls << " 回)" << std::endl;
    }
    // 退避・置換の計数 (累積。増えたときだけ表示)
    static unsigned long long s_hll = 0, s_vac = 0;
    unsigned long long hll = 0, vac = 0;
    CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&hll, g_ffHll, sizeof(unsigned long long)));
    CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&vac, g_ffVac, sizeof(unsigned long long)));
    if (hll != s_hll || vac != s_vac) {
        std::cout << "[farfield] 累積: HLL 退避 " << hll << " 面・回、真空/非物理の置換 " << vac << " 面・回 (評価区間では 0 が合格条件)" << std::endl;
        s_hll = hll; s_vac = vac;
    }
}

void convectiveFlux_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_ns)
{
    // free-stream 保存: 基準静圧 pRef と参照一様流を device 定数へ転送 (既定 0.0 でビット不変)。
    // cfg 値は run 中不変なので値が変わったときだけ転送する (毎ステップの同期 H2D ×5 を回避, 2026-09-12)。
    static double s_pRef = -1.0e300, s_roRef = -1.0e300, s_uRefX = -1.0e300, s_uRefY = -1.0e300, s_uRefZ = -1.0e300;
    if (s_pRef != cfg.pRef) {
        flow_float pRef_h = static_cast<flow_float>(cfg.pRef);
        CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_pRef, &pRef_h, sizeof(flow_float)));
        s_pRef = cfg.pRef;
    }
    if (s_roRef != cfg.roRef || s_uRefX != cfg.uRefX || s_uRefY != cfg.uRefY || s_uRefZ != cfg.uRefZ) {
        s_roRef = cfg.roRef; s_uRefX = cfg.uRefX; s_uRefY = cfg.uRefY; s_uRefZ = cfg.uRefZ;
        flow_float roRef_h = static_cast<flow_float>(cfg.roRef);
        flow_float uRefX_h = static_cast<flow_float>(cfg.uRefX);
        flow_float uRefY_h = static_cast<flow_float>(cfg.uRefY);
        flow_float uRefZ_h = static_cast<flow_float>(cfg.uRefZ);
        CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_roRef, &roRef_h, sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_uRefX, &uRefX_h, sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_uRefY, &uRefY_h, sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_uRefZ, &uRefZ_h, sizeof(flow_float)));
    }

    // D2a/D4 診断フラグを env から 1 度だけ device へ設定 (既定 off = ビット不変)。
    {
        static bool s_init = false;
        if (!s_init) {
            int c1 = 0, clog = 0; flow_float th = 0.05f, lth = 0.3f, blend = 0.0f;
            if (const char* e = getenv("FORGE_CONTACT_1ST"))       c1   = atoi(e);
            if (const char* e = getenv("FORGE_CONTACT_THRESH"))    th   = (flow_float)atof(e);
            if (const char* e = getenv("FORGE_CONTACT_LOG"))       clog = atoi(e);
            if (const char* e = getenv("FORGE_CONTACT_LOG_THRESH"))lth  = (flow_float)atof(e);
            if (const char* e = getenv("FORGE_CONTACT_BLEND"))     blend= (flow_float)atof(e);
            int fthy = 0;
            if (const char* e = getenv("FORGE_FACE_THERMOY"))      fthy = atoi(e);
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_faceThermoY,      &fthy, sizeof(int)));
            int fhd = 0;
            if (const char* e = getenv("FORGE_DIAG_FACE_H_DOUBLE")) fhd = atoi(e);
            if (fhd) printf("[DIAG] FORGE_DIAG_FACE_H_DOUBLE: SLAU の TP 多成分の面エンタルピーを double で評価する (診断)\n");
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_faceHDouble,      &fhd,  sizeof(int)));
            const int rem = (cfg.discretization == "node") ? 1 : 0;   // node は常にエッジ中点再構成
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_reconEdgeMid,     &rem,  sizeof(int)));
            const int sfr = cfg.speciesFaceReconstruction;   // config 由来 (env でない)
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_speciesFaceRecon, &sfr,  sizeof(int)));
            const int rycl = cfg.multispeciesRhoYCommonLimiter;   // config 由来 (opt-in 診断)
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_rhoYCommonLim,    &rycl, sizeof(int)));
            const unsigned long long zero = 0ULL;
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_speciesOvershoot, &zero, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contact1st,       &c1,   sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactThresh,    &th,   sizeof(flow_float)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactLog,       &clog, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactLogThresh, &lth,  sizeof(flow_float)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactBlend,     &blend,sizeof(flow_float)));
            if (const char* e = getenv("FORGE_DIAG_FACE_VEL_CELL")) {   // 診断介入 (数値を変える。既定 off)
                std::string v(e);
                const size_t c = v.find(':');
                if (c != std::string::npos) {
                    const long long fid = std::atoll(v.substr(0, c).c_str()), nid = std::atoll(v.substr(c + 1).c_str());
                    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_diagVelCellFace, &fid, sizeof(long long)));
                    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_diagVelCellNode, &nid, sizeof(long long)));
                    std::cout << "[FORGE_DIAG_FACE_VEL_CELL] **警告: 数値を変える診断介入** 面 " << fid << " の節点 " << nid
                              << " 側の再構成速度をセル値に戻す (生産で使わない)" << std::endl;
                }
            }
            const int brd = ((cfg.badReconDiag > 0) ? 1 : 0);
            const flow_float brro = (flow_float)cfg.roMin, brp = (flow_float)cfg.pMin;
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_badReconDiag,  &brd,  sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_badReconRoMin, &brro, sizeof(flow_float)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_badReconPMin,  &brp,  sizeof(flow_float)));
            // W2 フォールバックの面カウンタ (§4.23)。opt-in のときだけ確保する。
            const int brf = cfg.badReconFallback;
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_badReconHyst, &brf, sizeof(int)));
            if (brf > 0) {
                signed char* cnt_d = nullptr;
                CHECK_CUDA_ERROR(cudaMalloc(&cnt_d, sizeof(signed char)*msh.nPlanes));
                CHECK_CUDA_ERROR(cudaMemset(cnt_d, 0, sizeof(signed char)*msh.nPlanes));
                CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_badReconCnt, &cnt_d, sizeof(signed char*)));
                std::cout << "[convectiveFlux] badReconFallback: " << brf
                          << " 回の訪問だけ 1 次化する面カウンタを " << msh.nPlanes << " 面ぶん確保した" << std::endl;
            }
            s_init = true;
        }
    }

    // W2 V1 (plan convection-node-wall-reconstruction §6.4): 非物理な再構成の発火を一定間隔で 1 行印字しリセット。
    if (cfg.badReconDiag > 0) {
        static int s_br_call = 0;
        const int interval = cfg.badReconDiag;
        if ((s_br_call % interval) == 0) {
            unsigned long long faces=0, nro=0, np=0, tot=0;
            gpuErrchkKernelSync();
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&faces, g_badReconFaces, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&nro,   g_badReconRo,    sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&np,    g_badReconP,     sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&tot,   g_badReconTotal, sizeof(unsigned long long)));
            unsigned long long act=0;
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&act, g_badReconActive, sizeof(unsigned long long)));
            printf("BADRECON call=%d faces=%llu/%llu [ro=%llu P=%llu] active1st=%llu\n",
                   s_br_call, faces, tot, nro, np, act);
            const unsigned long long z = 0ULL;
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_badReconFaces, &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_badReconRo,    &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_badReconP,     &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_badReconTotal, &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_badReconActive, &z, sizeof(unsigned long long)));
        }
        s_br_call++;
    }

    // rho-Y 共通リミタ診断: 一定間隔で device カウンタを読み出して 1 行印字しリセット (opt-in 時のみ)。
    if (cfg.multispeciesRhoYCommonLimiter == 1) {
        static int s_rycl_call = 0;
        const int interval = 200;
        if ((s_rycl_call % interval) == 0) {
            int psimin=0, ymin=0, ymax=0;
            unsigned long long lt001=0, lt01=0, byrho=0, bysp=0, fb=0, ovs=0;
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&psimin, g_psiRhoY_min_scaled, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&lt001,  g_psiRhoY_lt001, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&lt01,   g_psiRhoY_lt01,  sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&byrho,  g_rhoYMinByRho,  sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&bysp,   g_rhoYMinBySpecies, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&fb,     g_rhoYFallback,  sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&ovs,    g_speciesOvershoot, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&ymin,   g_Yface_min_scaled, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&ymax,   g_Yface_max_scaled, sizeof(int)));
            printf("RHOYLIM call=%d minPsiRhoY=%.4f lt0.01=%llu lt0.1=%llu minBy[rho=%llu sp=%llu] overshoot=%llu fallback=%llu Yface=[%.5f,%.5f]\n",
                   s_rycl_call, (double)psimin*1e-6, lt001, lt01, byrho, bysp, ovs, fb,
                   (double)ymin*1e-6, (double)ymax*1e-6);
            // 次区間用にリセット (min/max は両端へ)。
            const int big=2000000, nbig=-2000000; const unsigned long long z=0ULL;
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_psiRhoY_min_scaled, &big, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_Yface_min_scaled,   &big, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_Yface_max_scaled,   &nbig, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_psiRhoY_lt001, &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_psiRhoY_lt01,  &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_rhoYMinByRho,  &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_rhoYMinBySpecies, &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_rhoYFallback,  &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_speciesOvershoot, &z, sizeof(unsigned long long)));
        }
        s_rycl_call++;
    }

    // initialize
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_ro"]  , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roUx"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roUy"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roUz"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roe"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.p_d["massflux"], 0.0, msh.nPlanes*sizeof(flow_float)));

    dim3 dimGrid_normal_halo = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));

    // node-centered 弱形式 (Phase 2): 主対流ループは「内部双対面のみ」を処理し、末尾に並ぶ全非 periodic 境界 plane
    // (nBoundaryHaloPlanes = wall+inlet+outlet+slip...) を除外する。全境界は別途 convectiveFlux_boundary_d が bvar を
    // R 状態とする弱形式で担う。これは境界ノードが物理境界上に乗る node-centered で ghost を主ループの右状態に食わせる
    // と退化幾何 (d_along_n=0) により出口列/コーナーで近壁 BL が崩壊するため (case/26 で実証)。
    // **periodic も除外する** (median-dual M4, §4.5): node の周期境界は DOF 同一視 (periodicNodeGather + 合併体積)
    // で扱い、継ぎ目に双対面を作らない。周期半割面を主ループに入れると外向き境界パッチを CV 間面と誤用し発散する。
    // よって node の主ループ境界 = 純内部双対面 (nNormalPlanes)。cell モードは全 plane を主ループでゴースト処理 (従来どおり)。
    const geom_int convPlaneBound = (cfg.discretization == "node")
                                  ? msh.nNormalPlanes : msh.nNormal_halo_Planes;

    // (検証A) node 弱形式の plane 振り分けを 1 度だけログ: 主ループ面 = 内部(+periodic)、弱形式面 = 全境界半割面。
    if (cfg.discretization == "node") {
        static bool s_planeLogDone = false;
        if (!s_planeLogDone) {
            s_planeLogDone = true;
            const geom_int mainPlanes = convPlaneBound;
            const geom_int periodicInMain = mainPlanes - msh.nNormalPlanes;  // 主ループのうち内部以外 = periodic
            printf("[node weak-form] conv main planes = %ld (internal %ld + periodic %ld), "
                   "boundary weak planes = %ld (wall %ld + non-wall %ld)\n",
                   (long)mainPlanes, (long)msh.nNormalPlanes, (long)periodicInMain,
                   (long)msh.nBoundaryHaloPlanes, (long)msh.nWallHaloPlanes,
                   (long)(msh.nBoundaryHaloPlanes - msh.nWallHaloPlanes));
        }
    }

    // -----------------------
    // *** sum over planes ***
    // -----------------------
    // 非平衡凝縮: 二相エネルギー流束補正用の総液相分率 g。off / 未登録なら nullptr (補正なし)。
    // 現状 nCondSpecies=1 (g_0 が総 g)。多成分は総和 g 配列を別途用意する (TODO)。
    flow_float* cond_g = nullptr;
    if (cfg.condensation == 1 && var.nCondSpeciesRegistered >= 1) cond_g = var.c_d["g_0"];

    // 全スキーム共通の引数バンドルを 1 度だけ構築 (SLAU/HLLE/ROE で使い回す)。
    FaceGeom geom {
        msh.nCells, msh.nPlanes, msh.nNormalPlanes, msh.map_plane_cells_d,
        convPlaneBound, msh.normal_halo_planes_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],
        var.p_d["massflux"],
        msh.wall_flag_d };
    PrimState st {
        var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"],
        var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"], var.c_d["P"], var.c_d["Ht"], var.c_d["sonic"],
        var.c_d["T"] };
    ResidualOut reso {
        var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"] };
    LimiterFields lim {
        var.c_d["limiter_ro"], var.c_d["limiter_Ux"], var.c_d["limiter_Uy"], var.c_d["limiter_Uz"], var.c_d["limiter_P"],
        var.c_d["ducros"], (cfg.reconT == 1 && var.c_d.count("limiter_T")) ? var.c_d["limiter_T"] : nullptr };
    GradFields grd {
        var.c_d["drodx"], var.c_d["drody"], var.c_d["drodz"],
        var.c_d["dUxdx"], var.c_d["dUxdy"], var.c_d["dUxdz"],
        var.c_d["dUydx"], var.c_d["dUydy"], var.c_d["dUydz"],
        var.c_d["dUzdx"], var.c_d["dUzdy"], var.c_d["dUzdz"],
        var.c_d["dPdx"] , var.c_d["dPdy"] , var.c_d["dPdz"],
        var.c_d["dTdx"] , var.c_d["dTdy"] , var.c_d["dTdz"] };
    // SST 全エネルギー E_t = E_m + ρk (sstEnergyIncludesK): 面エンタルピー +(5/3)k, 圧力流束 p* = p + (2/3)ρk。
    const bool sstEnergyK = (cfg.sstEnergyIncludesK != 0 && cfg.LESorRANS == 2 && cfg.RANSmodel == 1);
    if (sstEnergyK && !(cfg.solver == "SLAU" || cfg.solver == "SLAU2")) {
        throw std::runtime_error("turbulence.sstEnergyIncludesK=1 is implemented for solver SLAU/SLAU2 only");
    }
    CondArgs cnd { cfg.cp, cond_g, var.c_d["T"], cfg.condModel, sstEnergyK ? var.c_d["k"] : nullptr, sstEnergyK ? 1 : 0 };
    cnd.cprops = condProps_make(cfg.condModel, cond_prop_opts(cfg));   // σ 倍率・N2 低温物性を面エンタルピーの潜熱にも反映
    cnd.Yw     = cfg.condVaporMassFraction;
    cnd.tables    = cond_tables_device();
    cnd.condFloat = (cfg.condFloat != 0 && cnd.tables.valid) ? 1 : 0;

    if (cfg.solver == "SLAU" || cfg.solver == "SLAU2") {
        int slauVariant = (cfg.solver == "SLAU2") ? 2 : 1;
        SpeciesArgs spA {
            cfg.thermalMethod, thermo_species_device_ptr(), thermo_species_device_ptr_f(), cfg.nSpecies,
            species_roY_device_ptr(),
            species_Y_device_ptr(), species_dYdx_device_ptr(), species_dYdy_device_ptr(), species_dYdz_device_ptr(),
            species_limiterY_device_ptr(),
            (cfg.speciesFaceReconstruction >= 2) ? species_Yface_alloc(msh.nPlanes) : nullptr,
            var.c_d["Rmix"] };
        // 受動種 S3 (passiveScalarScheme 1 かつ speciesFaceReconstruction>=2): ψ_P で再構成した upwind 面値を Pface へ。
        if (passiveSchemeEnabled(cfg) && cfg.speciesFaceReconstruction >= 2) {
            spA.nPassive       = passive_count();
            spA.nPassiveUnit   = (passive_tracer_index() == 0) ? 1 : 0;
            spA.P_recon        = passive_P_device_ptr();
            spA.dPdx_recon     = passive_dPdx_device_ptr();
            spA.dPdy_recon     = passive_dPdy_device_ptr();
            spA.dPdz_recon     = passive_dPdz_device_ptr();
            spA.limiterP_recon = passive_limiter_device_ptr();
            spA.Pface_out      = passive_Pface_alloc(msh.nPlanes);
        }

        SLAU_d<<<dimGrid_normal_halo , cuda_cfg.dimBlock>>> (
            cfg.convMethod, cfg.limiter, slauVariant, cfg.reconT,
            (cfg.slauWallNormalChi > 0 ? 1 : 0),   // 未解決の auto (-1) を有効扱いにしない
            cfg.slauContactFloor,
            cfg.lowMachPrecond, cfg.precondEps,
            cfg.lowMachThornber,
            cfg.gamma,
            spA, cnd, geom, st, reso, lim, grd
        ) ;

    } else if (cfg.solver == "HLLE") {
        HLLE_d<<<dimGrid_normal_halo , cuda_cfg.dimBlock>>> (
            cfg.convMethod, cfg.limiter,
            cfg.gamma,
            cnd, geom, st, reso, lim, grd
        );

    } else if (cfg.solver == "ROE") {
        //ROE_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>> ( 
        ROE_d<<<dimGrid_normal_halo , cuda_cfg.dimBlock>>> (
            cfg.convMethod, cfg.limiter,
            cfg.gamma,
            cfg.roeEntropyFixCoeff,
            cfg.thermalMethod, thermo_species_device_ptr(), cfg.nSpecies,
            cnd, geom, st, reso, lim, grd
        );

    } else if (cfg.solver == "KEEP") {
        // 純粋 KEEP 中心流束 + opt-in ES 散逸レイヤ (keepDissType, 既定 0=散逸なし・ビット不変)。
        // LES/ILES 向け低散逸対流。SGS 散逸は WALE が担う。
        // TP (thermalMethod==2): 単成分 (Step 3)・多成分 (Step 4) とも scalar/matrix 対応。
        // 多成分 matrix の市松プラトーバグは根治済 (真因: ΣρY≠ρ 共通モードノイズの w 増幅。
        // カーネル内 Y 正規化で除去。plan 変更ログ 2026-07-19 深掘り参照)。
        KEEP_d<<<dimGrid_normal_halo , cuda_cfg.dimBlock>>> (
            cfg.gamma,
            cfg.thermalMethod, thermo_species_device_ptr(),
            cfg.nSpecies, species_roY_device_ptr(),
            cfg.keepDissType, cfg.keepDissCoeff,
            cfg.keepDissCprime, cfg.precondEps,
            cfg.keepDissJump, cfg.keepDissPrecond,
            // f_d 駆動 σ ブレンド (keepDissFdBlend=1 かつ DES のときのみ fd_shield を渡す。
            // config 検証で DESmode>0 は保証済み → fd_shield は毎サブ反復 turbulent_viscosity が更新)
            (cfg.keepDissFdBlend == 1) ? var.c_d["fd_shield"] : nullptr,
            (cfg.DESmode == 1) ? 1 : 0,   // fd_shield の向き (DDES: 1=LES / IDDES: 1=RANS)
            cfg.keepDissCoeffMax,
            cfg.keepDissCbCoeff, cfg.keepDissCbEps,
            cfg.keepDissOpBlendRaw,
            grd,
            geom, st, reso
        );

    } else {
        std::cerr << "Error: unsupported solver name " << cfg.solver
                  << " (enabled: SLAU, HLLE, ROE, KEEP)" << std::endl;
        exit(EXIT_FAILURE);
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();


    // cell モード: SLAU/ROE/HLLE は主ループが normal_halo_planes_d 経由で全境界 plane をゴースト処理するため、
    // この dedicated 境界カーネルは二重計上になるのでスキップする。
    // node モード (弱形式・Phase 2): 主ループは全非 periodic 境界 plane を除外した (convPlaneBound=内部+periodic) ので、
    // 全境界 (wall/inlet/outlet/slip...) の境界寄与を本 convectiveFlux_boundary_d が bvar を R 状態とする弱形式で担う。
    // (壁: bvar u=0 → mdot=0, pressure-only。出口: bvar Uxb/Psb。slip: Uxb=0+Psb=P[ic]。入口: 固定状態。)
    const bool nodeMode = (cfg.discretization == "node");
    bool skipBoundaryFluxKernel = (!nodeMode)
                               && (cfg.solver == "SLAU"
                                || cfg.solver == "SLAU2"
                                || cfg.solver == "ROE"
                                || cfg.solver == "HLLE"
                                || cfg.solver == "KEEP");

    for (auto& bc : msh.bconds)
    {
        if (bc.bcondKind == "periodic") {
            continue;
        }
        if (skipBoundaryFluxKernel) {
            continue;
        }
        if (bc.bcondKind == "farfield") {   // 遠方境界: 専用の HLLC 流束 (既存の境界流束は通らない)
            farfieldFlux_d_wrapper(cfg, cuda_cfg, bc, msh, var, sstEnergyK ? 1 : 0);
            continue;
        }
        convectiveFlux_boundary_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> (
            cfg.gamma,
            // 移流基準差分 (plan §8.5): 主ループ KEEP_d の advGauge (d_roRef>0 && CPG) と厳密同条件。
            // CV の全面に載らないと telescoping が破れるため、境界だけの on/off は不可。
            ((cfg.solver == "KEEP" && cfg.roRef > 0.0 && cfg.thermalMethod != 2) ? 1 : 0),
            // node × outlet_statPress: p_tilde を内部値 Ps[ic] で評価 (§2.11)。cell 経路は不変。
            ((nodeMode && bc.bcondKind == "outlet_statPress") ? 1 : 0),
            // mesh structure
            bc.iPlanes.size(),
            bc.map_bplane_plane_d,  
            bc.map_bplane_cell_d,  
            bc.map_bplane_cell_ghst_d,

            var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
            var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
            var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],
            (nodeMode ? var.p_d["massflux"] : nullptr),   // node 弱形式境界のみ massflux 書き戻し (cell は主ループが書く)

            // basic variables
            var.c_d["ro"] ,
            var.c_d["roUx"] ,
            var.c_d["roUy"] ,
            var.c_d["roUz"] ,
            var.c_d["roe"] ,
            var.c_d["Ux"]  ,
            var.c_d["Uy"]  ,
            var.c_d["Uz"]  ,
            var.c_d["P"]  ,
            var.c_d["Ht"]  ,
            var.c_d["sonic"]  ,
            var.c_d["T"]  ,

            bc.bvar_d["ro"],
            bc.bvar_d["roUx"],
            bc.bvar_d["roUy"],
            bc.bvar_d["roUz"],
            bc.bvar_d["roe"],
            bc.bvar_d["Ux"],
            bc.bvar_d["Uy"],
            bc.bvar_d["Uz"],
            bc.bvar_d["Tt"],
            bc.bvar_d["Pt"],
            bc.bvar_d["Ts"],
            bc.bvar_d["Ps"],
 
            var.c_d["res_ro"] ,
            var.c_d["res_roUx"] ,
            var.c_d["res_roUy"] ,
            var.c_d["res_roUz"] ,
            var.c_d["res_roe"]  ,
            // sstEnergyIncludesK: 内部側 k と境界側 k (bvar kb / 入口 k / 無ければ内部値)
            sstEnergyK ? var.c_d["k"] : nullptr,
            // 境界側 k: 入口 (Dirichlet, bvar "k" = 指定値) のみ bvar を使う。壁/出口/slip の bvar kb は node では未充填になり得るため内部値 (Neumann)
            (sstEnergyK && bc.bcondKind.rfind("inlet", 0) == 0 && bc.bvar_d.count("k")) ? bc.bvar_d["k"] : nullptr,
            sstEnergyK ? 1 : 0
        ) ;
    }


    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

    // 診断ダンプ (env `FORGE_DUMP_MASSFLUX=<path>`、既定 off。**数値の振る舞いは変えない**)。
    // 面流束 massflux[nPlanes] を **最初の呼び出しだけ** ホストへ写して raw float32 で書く。
    // plan convection-slau-wall-normal-chi §6 V5 (#10d): 場 (res_*.h5) は残差の atomicAdd で
    // 1 step でもビット再現しないが、massflux[ip] は 1 面 = 1 スレッドが非 atomic に書くので
    // 面レベルのビット比較ができる。**第 2 評価以降は flag 0/1 で状態が違う**ので 1 回目に限る。
    {
        static bool s_mfDumped = false;
        static int  s_cfCall = 0;
        ++s_cfCall;                      // この wrapper の呼び出し回数 (1 = 第 1 評価)
        if (!s_mfDumped) {
            const char* mfPath = std::getenv("FORGE_DUMP_MASSFLUX");
            if (mfPath && *mfPath) {
                s_mfDumped = true;
                std::vector<flow_float> mf(msh.nPlanes);
                CHECK_CUDA_ERROR(cudaMemcpy(mf.data(), var.p_d["massflux"],
                                            msh.nPlanes*sizeof(flow_float), cudaMemcpyDeviceToHost));
                std::ofstream ofs(mfPath, std::ios::binary);
                if (ofs) {
                    ofs.write(reinterpret_cast<const char*>(mf.data()),
                              (std::streamsize)(msh.nPlanes*sizeof(flow_float)));
                    std::cout << "[FORGE_DUMP_MASSFLUX] wrote " << msh.nPlanes
                              << " faces to " << mfPath << " (call " << s_cfCall << ")\n";
                } else {
                    std::cout << "[FORGE_DUMP_MASSFLUX] cannot open " << mfPath << '\n';
                }
                // **カーネルが実際に読む状態**も同じ呼び出しで書く。res_*.h5 の値を代理に使うと
                // 1 ulp ずれることがある (§6 V5 の P2 初版が 3 面で FAIL した原因)。
                // 並び: ro, Ux, Uy, Uz, Ps, sonic を各 nCells 個、この順に連結。
                {
                    const std::string sPath = std::string(mfPath) + ".state";
                    std::ofstream sfs(sPath, std::ios::binary);
                    if (sfs) {
                        std::vector<flow_float> buf(msh.nCells);
                        for (const char* nm : {"ro", "Ux", "Uy", "Uz", "P", "sonic"}) {
                            CHECK_CUDA_ERROR(cudaMemcpy(buf.data(), var.c_d[nm],
                                                        msh.nCells*sizeof(flow_float), cudaMemcpyDeviceToHost));
                            sfs.write(reinterpret_cast<const char*>(buf.data()),
                                      (std::streamsize)(msh.nCells*sizeof(flow_float)));
                        }
                        std::cout << "[FORGE_DUMP_MASSFLUX] wrote state (ro,Ux,Uy,Uz,P,sonic x "
                                  << msh.nCells << ") to " << sPath << '\n';
                    }
                }
            }
        }
    }

}
// =============================================================================
// 局所帳簿ダンプ (plan tooling-nozzle-sern-3d §5.1 R5h、codex diagnose 2026-09-27)。
//   FORGE_DUMP_LEDGER=<path>        : 出力先 (CSV)。未設定なら全関数 no-op (解はビット同一)。
//   FORGE_DUMP_LEDGER_NODES=<a,b,..>: 印を付ける節点 ID (0 始まり)。all で全節点。
//   FORGE_DUMP_LEDGER_CALLS=<n>     : 記録する assembleResidual の呼び出し数 (既定 2)。
// 出力専用: 状態・残差をホストへ写して書くだけで、どの配列も書き換えない。
//   <path>          行 = call,tag,node,field,value (段ごとの状態・残差)
//   <path>.faces    行 = call + SLAU_d が記録した面の LEDGER_FACE_NF 列 (列名は 1 行目)
// =============================================================================
namespace {
struct LedgerState {
    bool init = false, on = false;
    std::string path;
    std::vector<long long> nodes;
    unsigned char* flag_d = nullptr;
    float* buf_d = nullptr;
    unsigned int cap = 0;
    int call = 0, maxCalls = 2;
};
LedgerState& ledger() { static LedgerState s; return s; }
}

static void ledgerInitOnce(mesh& msh)
{
    LedgerState& L = ledger();
    if (L.init) return;
    L.init = true;
    const char* p = std::getenv("FORGE_DUMP_LEDGER");
    const char* n = std::getenv("FORGE_DUMP_LEDGER_NODES");
    if (!p || !*p || !n || !*n) return;
    L.path = p;
    if (const char* c = std::getenv("FORGE_DUMP_LEDGER_CALLS")) L.maxCalls = std::max(1, std::atoi(c));
    std::string s(n);
    size_t pos = 0;
    if (s == "all") {   // 全節点 (環境変数に全 ID を並べると ARG_MAX を超える)
        for (long long id = 0; id < (long long)msh.nCells; ++id) L.nodes.push_back(id);
        pos = s.size();
    }
    while (pos < s.size()) {
        size_t q = s.find(',', pos);
        if (q == std::string::npos) q = s.size();
        if (q > pos) {
            const long long id = std::atoll(s.substr(pos, q - pos).c_str());
            if (id >= 0 && id < (long long)msh.nCells) L.nodes.push_back(id);
        }
        pos = q + 1;
    }
    if (L.nodes.empty()) { std::cout << "[FORGE_DUMP_LEDGER] 有効な節点が無いので無効\n"; return; }
    std::vector<unsigned char> flag(msh.nCells, 0);
    for (long long id : L.nodes) flag[id] = 1;
    CHECK_CUDA_ERROR(cudaMalloc(&L.flag_d, msh.nCells));
    CHECK_CUDA_ERROR(cudaMemcpy(L.flag_d, flag.data(), msh.nCells, cudaMemcpyHostToDevice));
    L.cap = (unsigned int)std::min<size_t>(200000, 64 * L.nodes.size() + 1024);
    CHECK_CUDA_ERROR(cudaMalloc(&L.buf_d, sizeof(float) * (size_t)L.cap * LEDGER_FACE_NF));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_ledgerFaceCap, &L.cap, sizeof(unsigned int)));
    const unsigned int zero = 0;
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_ledgerFaceCount, &zero, sizeof(unsigned int)));
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_ledgerFaceBuf, &L.buf_d, sizeof(float*)));
    { std::ofstream o(L.path); o << "call,tag,node,field,value\n"; }
    {
        std::ofstream o(L.path + ".faces");
        o << "call,ip,ic0,ic1,sx,sy,sz,ss,ro_L,ro_R,P_L,P_R,Pf_L,Pf_R,Ux_L,Uy_L,Uz_L,Ux_R,Uy_R,Uz_R,h_p,h_m,c_hat,M_hat,chi,chi_mass,"
             "Vn_p,Vn_m,p_tilde_r,mdot,F_ro,F_roUx,F_roUy,F_roUz,F_roe,conv_scheme,P_del,fx,limiter_ro_0,limiter_ro_1,limiter_P_0\n";
    }
    L.on = true;
    std::cout << "[FORGE_DUMP_LEDGER] " << L.nodes.size() << " 節点に印、最初の " << L.maxCalls
              << " 回の assembleResidual を " << L.path << " に記録する\n";
}

void ledgerBeginAssemble(mesh& msh)
{
    ledgerInitOnce(msh);
    LedgerState& L = ledger();
    if (!L.on) return;
    ++L.call;
    // 記録する回だけ SLAU_d の面記録を有効にする (それ以外は nullptr で完全に no-op)
    const unsigned char* f = (L.call <= L.maxCalls) ? L.flag_d : nullptr;
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_ledgerFlag, &f, sizeof(const unsigned char*)));
}

void ledgerCapture(mesh& msh, variables& var, const char* tag, bool residual)
{
    LedgerState& L = ledger();
    if (!L.on || L.call > L.maxCalls) return;
    static const char* stateF[] = {"ro","roUx","roUy","roUz","roe","roK","roOmega","roY0","roY1","T","P","Ux","Uy","Uz","Y0","Y1","Ht","sonic","Rmix"};
    static const char* resF[]   = {"res_ro","res_roUx","res_roUy","res_roUz","res_roe","res_roK","res_roOmega","res_roY0","res_roY1"};
    std::ofstream o(L.path, std::ios::app);
    o.precision(9);
    const size_t nf = residual ? sizeof(resF)/sizeof(resF[0]) : sizeof(stateF)/sizeof(stateF[0]);
    for (size_t i = 0; i < nf; ++i) {
        const char* nm = residual ? resF[i] : stateF[i];
        auto it = var.c_d.find(nm);
        if (it == var.c_d.end() || it->second == nullptr) continue;
        for (long long id : L.nodes) {
            flow_float v;
            CHECK_CUDA_ERROR(cudaMemcpy(&v, it->second + id, sizeof(flow_float), cudaMemcpyDeviceToHost));
            o << L.call << ',' << tag << ',' << id << ',' << nm << ',' << (double)v << '\n';
        }
    }
    (void)msh;
}

void ledgerFlushFaces()
{
    LedgerState& L = ledger();
    if (!L.on || L.call > L.maxCalls) return;
    unsigned int cnt = 0;
    CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&cnt, g_ledgerFaceCount, sizeof(unsigned int)));
    const unsigned int n = std::min(cnt, L.cap);
    std::vector<float> buf((size_t)n * LEDGER_FACE_NF);
    if (n > 0) CHECK_CUDA_ERROR(cudaMemcpy(buf.data(), L.buf_d, sizeof(float) * buf.size(), cudaMemcpyDeviceToHost));
    std::ofstream o(L.path + ".faces", std::ios::app);
    o.precision(9);
    auto asInt = [](float x) { int i; std::memcpy(&i, &x, sizeof(int)); return i; };
    for (unsigned int k = 0; k < n; ++k) {
        const float* r = buf.data() + (size_t)k * LEDGER_FACE_NF;
        o << L.call;
        for (int j = 0; j < LEDGER_FACE_NF; ++j) {
            o << ',';
            if (j == 0 || j == 1 || j == 2 || j == 34) o << asInt(r[j]); else o << (double)r[j];
        }
        o << '\n';
    }
    if (cnt > L.cap) std::cout << "[FORGE_DUMP_LEDGER] 面バッファ不足: " << cnt << " > " << L.cap << "\n";
    const unsigned int zero = 0;
    CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_ledgerFaceCount, &zero, sizeof(unsigned int)));
    std::cout << "[FORGE_DUMP_LEDGER] call " << L.call << ": " << n << " 面を記録\n";
}
