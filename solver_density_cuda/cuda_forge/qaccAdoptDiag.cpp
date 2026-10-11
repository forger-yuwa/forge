// 診断 FORGE_DIAG_QACC_ADOPT (plans/active/architecture-float-state-double-geometry.md §5.1 #18・§6.34、既定 off)。
// 記録する量と時点は qaccAdoptDiag.hpp。ここはホスト側の読み取り (D2H) と mask の集計・h5 の書き出し。

#include "cuda_forge/qaccAdoptDiag.hpp"
#include "cuda_forge/cudaWrapper.cuh"

#include <highfive/H5File.hpp>

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <set>
#include <string>
#include <vector>

namespace qaccAdoptDiag {

namespace detail {
bool g_enabled = false;
}

namespace {

constexpr int kNQ = 5;
const char* const kQ[kNQ] = {"ro", "roUx", "roUy", "roUz", "roe"};

struct State {
    std::string prefix;          // 出力の接頭辞 (<接頭辞>_<step>.h5)
    std::set<int> steps;         // 対象の step (iStep + 1)
    std::set<int> written;       // 書いた step
};
State g;

template <typename T>
std::vector<T> d2hFlag(const T* p, size_t n)
{
    std::vector<T> h;
    if (p == nullptr || n == 0) return h;
    h.resize(n);
    gpuErrchk( cudaMemcpy(h.data(), p, n*sizeof(T), cudaMemcpyDeviceToHost) );
    return h;
}

template <typename T>
void putIfAny(HighFive::File& h5, const std::string& name, const std::vector<T>& v)
{
    if (!v.empty()) h5.createDataSet(name, v);
}

// flag が無い (空) ときは条件を満たさないとみなす。
inline bool flagOn(const std::vector<geom_int>& f, size_t i) { return !f.empty() && f[i] == 1; }

// "1,30000" を読む。不正なら止める (黙って違う step を書かない)。
std::set<int> parseSteps(const char* s)
{
    std::set<int> out;
    std::string str(s);
    size_t pos = 0;
    while (true) {
        const size_t comma = str.find(',', pos);
        const std::string tok = str.substr(pos, (comma == std::string::npos) ? std::string::npos : comma - pos);
        errno = 0;
        char* end = nullptr;
        const long n = std::strtol(tok.c_str(), &end, 10);
        if (tok.empty() || errno != 0 || end == tok.c_str() || *end != '\0' || n < 1 || n > 1000000000L) {
            fprintf(stderr, "[qaccAdopt] FORGE_DIAG_QACC_ADOPT_STEPS='%s' の '%s' は 1 以上の整数でない。止める\n", s, tok.c_str());
            std::exit(EXIT_FAILURE);
        }
        out.insert((int)n);
        if (comma == std::string::npos) break;
        pos = comma + 1;
    }
    return out;
}

std::string stepsText(const std::set<int>& st)
{
    std::string t;
    for (int v : st) t += (t.empty() ? "" : ",") + std::to_string(v);
    return t;
}

}  // namespace

void init(solverConfig& cfg, mesh& msh, variables& var)
{
    const char* prefix = std::getenv("FORGE_DIAG_QACC_ADOPT");
    const char* stepsEnv = std::getenv("FORGE_DIAG_QACC_ADOPT_STEPS");
    if (prefix == nullptr || *prefix == '\0') {
        if (stepsEnv != nullptr && *stepsEnv != '\0')
            printf("[qaccAdopt] FORGE_DIAG_QACC_ADOPT_STEPS=%s は FORGE_DIAG_QACC_ADOPT が無いので無視する\n", stepsEnv);
        return;
    }
    g = State{};
    g.prefix = prefix;
    g.steps = (stepsEnv != nullptr && *stepsEnv != '\0') ? parseSteps(stepsEnv) : std::set<int>{1};

    // qAccumulatorFP64 0 なら正本が無いので何もしない (計算は続ける)。
    if (cfg.qAccumulatorFP64 != 1 || var.qacc_d[0] == nullptr) {
        printf("[qaccAdopt] qAccumulatorFP64 が 0 (FP64 の正本が無い) なので何もしない。FORGE_DIAG_QACC_ADOPT=%s には書かない\n",
               prefix);
        fflush(stdout);
        return;
    }
    detail::g_enabled = true;
    printf("[qaccAdopt] 有効: step %s (= iStep + 1) の reconcile の直前に、量ごとの採用の mask (Q != (float)Qacc) と Q・Qacc・"
           "座標・flag を %s_<step>.h5 に書く (内点 %ld CV)\n", stepsText(g.steps).c_str(), prefix, (long)msh.nCells);
    if (cfg.mainLoopCount() < *g.steps.rbegin())
        printf("[qaccAdopt] 注意: nStepOuter (%d) が対象の最大の step %d より小さい。その step は書かない\n",
               cfg.mainLoopCount(), *g.steps.rbegin());
    fflush(stdout);
}

void beforeReconcile(solverConfig& cfg, mesh& msh, variables& var, int iStep)
{
    if (!detail::g_enabled) return;
    const int step = iStep + 1;
    if (g.steps.count(step) == 0 || g.written.count(step) != 0) return;
    // 起動時に確かめているが、黙って違うものを書かないようここでも確かめる。
    if (var.qacc_d[0] == nullptr) {
        printf("[qaccAdopt] 未対応: step %d: FP64 の正本が確保されていない。何も書かない (計算は続ける)\n", step);
        fflush(stdout);
        g.written.insert(step);
        return;
    }

    const size_t nC = (size_t)msh.nCells;
    std::vector<std::vector<double>> Q(kNQ), Qacc(kNQ);
    std::vector<std::vector<uint8_t>> mask(kNQ);
    for (int k = 0; k < kNQ; ++k) {
        std::vector<flow_float> q(nC);
        Qacc[k].resize(nC);
        if (nC > 0) {
            gpuErrchk( cudaMemcpy(q.data(), var.c_d.at(kQ[k]), nC*sizeof(flow_float), cudaMemcpyDeviceToHost) );
            gpuErrchk( cudaMemcpy(Qacc[k].data(), var.qacc_d[k], nC*sizeof(double), cudaMemcpyDeviceToHost) );
        }
        mask[k].resize(nC);
        Q[k].resize(nC);
        for (size_t i = 0; i < nC; ++i) {
            // reconcile (qaccReconcileCell) と同じ判定: q != (flow_float)qacc (NaN も採用に数える)。
            mask[k][i] = (q[i] != static_cast<flow_float>(Qacc[k][i])) ? 1 : 0;
            Q[k][i] = static_cast<double>(q[i]);
        }
    }
    const std::vector<geom_int> axisF  = d2hFlag(msh.axis_flag_d, nC);
    const std::vector<geom_int> wallF  = d2hFlag(msh.wall_flag_d, nC);
    const std::vector<geom_int> isoF   = d2hFlag(msh.iso_wall_flag_d, nC);
    const std::vector<geom_int> bnodeF = d2hFlag(msh.bnode_flag_d, nC);

    // 格子の節点の数 (flag ごと)。
    long long nAxis = 0, nWall = 0, nIso = 0, nBnode = 0, nAxisWall = 0, nAxisIso = 0;
    for (size_t i = 0; i < nC; ++i) {
        const bool a = flagOn(axisF, i), w = flagOn(wallF, i);
        nAxis += a; nWall += w; nIso += flagOn(isoF, i); nBnode += flagOn(bnodeF, i);
        nAxisWall += (a && w); nAxisIso += (a && flagOn(isoF, i));
    }
    // 量ごとの採用数と内訳。内点 = bnode_flag 0 (bnode_flag が無ければ -1)。
    std::vector<long long> nAdopt(kNQ, 0), nAxisA(kNQ, 0), nWallA(kNQ, 0), nIsoA(kNQ, 0), nAxisWallA(kNQ, 0),
                           nInteriorA(kNQ, bnodeF.empty() ? -1 : 0);
    for (int k = 0; k < kNQ; ++k) {
        for (size_t i = 0; i < nC; ++i) {
            if (!mask[k][i]) continue;
            ++nAdopt[k];
            const bool a = flagOn(axisF, i), w = flagOn(wallF, i);
            nAxisA[k] += a; nWallA[k] += w; nIsoA[k] += flagOn(isoF, i); nAxisWallA[k] += (a && w);
            if (!bnodeF.empty() && bnodeF[i] == 0) ++nInteriorA[k];
        }
    }

    const std::string path = g.prefix + "_" + std::to_string(step) + ".h5";
    try {
        HighFive::File h5(path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
        h5.createAttribute("plan", std::string("plans/active/architecture-float-state-double-geometry.md §5.1 #18・§6.34"));
        h5.createAttribute("step", (long long)step);
        h5.createAttribute("iStep", (long long)iStep);
        h5.createAttribute("when", std::string("right before updateVariablesOuter (the reconcile) in advanceImplicitSteady: "
                                               "after the commit, SST, species and the post-commit isothermal-wall pin"));
        h5.createAttribute("mask_rule", std::string("mask = (Q != (flow_float)Qacc), same predicate as qaccReconcileCell "
                                                    "(1 = the reconcile adopts the FP32 value)"));
        h5.createAttribute("flow_float_bytes", (long long)sizeof(flow_float));
        h5.createAttribute("nCells", (long long)msh.nCells);
        h5.createAttribute("quantities", std::string("ro,roUx,roUy,roUz,roe"));
        h5.createAttribute("discretization", cfg.discretization);
        h5.createAttribute("mesh_file", cfg.meshFileName);
        h5.createAttribute("isAxisymmetric", (long long)cfg.isAxisymmetric);
        h5.createAttribute("qAccumulatorFP64", (long long)cfg.qAccumulatorFP64);
        const char* axDiag = std::getenv("FORGE_QACC_AXISYM_DIAG");
        h5.createAttribute("FORGE_QACC_AXISYM_DIAG", std::string((axDiag && *axDiag) ? axDiag : "(unset)"));
        h5.createAttribute("n_adopt_per_qty", nAdopt);
        h5.createAttribute("n_adopt_axis_per_qty", nAxisA);
        h5.createAttribute("n_adopt_wall_per_qty", nWallA);
        h5.createAttribute("n_adopt_iso_wall_per_qty", nIsoA);
        h5.createAttribute("n_adopt_axis_and_wall_per_qty", nAxisWallA);
        h5.createAttribute("n_adopt_interior_per_qty", nInteriorA);
        h5.createAttribute("interior_rule", std::string("interior = bnode_flag 0 (-1 if bnode_flag is absent)"));
        h5.createAttribute("n_nodes_axis", nAxis);
        h5.createAttribute("n_nodes_wall", nWall);
        h5.createAttribute("n_nodes_iso_wall", nIso);
        h5.createAttribute("n_nodes_bnode", nBnode);
        h5.createAttribute("n_nodes_axis_and_wall", nAxisWall);
        h5.createAttribute("n_nodes_axis_and_iso_wall", nAxisIso);
        for (int k = 0; k < kNQ; ++k) {
            h5.createDataSet(std::string("/") + kQ[k] + "/mask", mask[k]);
            h5.createDataSet(std::string("/") + kQ[k] + "/Q", Q[k]);
            h5.createDataSet(std::string("/") + kQ[k] + "/Qacc", Qacc[k]);
        }
        if (msh.cc64.size() >= 3*nC) putIfAny(h5, "/cells/cc64", std::vector<double>(msh.cc64.begin(), msh.cc64.begin() + 3*nC));
        putIfAny(h5, "/cells/axis_flag",     axisF);
        putIfAny(h5, "/cells/wall_flag",     wallF);
        putIfAny(h5, "/cells/iso_wall_flag", isoF);
        putIfAny(h5, "/cells/bnode_flag",    bnodeF);
    } catch (const std::exception& ex) {
        fprintf(stderr, "[qaccAdopt] %s を書けない: %s。止める\n", path.c_str(), ex.what());
        std::exit(EXIT_FAILURE);
    }
    g.written.insert(step);

    printf("[qaccAdopt] step %d: %s に書いた。採用 (Q != (float)Qacc) ro %lld roUx %lld roUy %lld roUz %lld roe %lld; "
           "うち軸 %lld/%lld/%lld/%lld/%lld、壁 %lld/%lld/%lld/%lld/%lld、等温壁 %lld/%lld/%lld/%lld/%lld、"
           "内点 %lld/%lld/%lld/%lld/%lld。格子: 軸 %lld・壁 %lld・等温壁 %lld・境界 %lld・軸と壁の重なり %lld 節点\n",
           step, path.c_str(), nAdopt[0], nAdopt[1], nAdopt[2], nAdopt[3], nAdopt[4],
           nAxisA[0], nAxisA[1], nAxisA[2], nAxisA[3], nAxisA[4],
           nWallA[0], nWallA[1], nWallA[2], nWallA[3], nWallA[4],
           nIsoA[0], nIsoA[1], nIsoA[2], nIsoA[3], nIsoA[4],
           nInteriorA[0], nInteriorA[1], nInteriorA[2], nInteriorA[3], nInteriorA[4],
           nAxis, nWall, nIso, nBnode, nAxisWall);
    fflush(stdout);
}

void finalize()
{
    if (detail::g_enabled) {
        std::set<int> missing;
        for (int s : g.steps) if (g.written.count(s) == 0) missing.insert(s);
        if (!missing.empty())
            printf("[qaccAdopt] 対象の step %s に達しなかった。その分は書いていない\n", stepsText(missing).c_str());
    }
    detail::g_enabled = false;
}

}  // namespace qaccAdoptDiag
