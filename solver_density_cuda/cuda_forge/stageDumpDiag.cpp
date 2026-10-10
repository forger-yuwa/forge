// 診断 FORGE_DIAG_STAGE_DUMP (plans/active/architecture-float-state-double-geometry.md §5.1 #17・§6.30、既定 off)。
// 記録する量と時点は stageDumpDiag.hpp。ここはホスト側の読み取り (D2H) と照合・h5 の書き出し。

#include "cuda_forge/stageDumpDiag.hpp"
#include "cuda_forge/cudaWrapper.cuh"

#include <highfive/H5File.hpp>

#include <array>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <string>
#include <vector>

namespace stageDumpDiag {

namespace detail {
bool g_enabled = false;
bool g_active  = false;
}

namespace {

constexpr int kNQ = 5;
const char* const kQ[kNQ]   = {"ro", "roUx", "roUy", "roUz", "roe"};
const char* const kRes[kNQ] = {"res_ro", "res_roUx", "res_roUy", "res_roUz", "res_roe"};
const char* const kN[kNQ]   = {"roN", "roUxN", "roUyN", "roUzN", "roeN"};

using Set = std::array<std::vector<double>, kNQ>;

struct Stage {
    bool got = false;
    Set v;
};

struct State {
    std::string path;          // 出力の h5
    int step = 1;              // 対象の step (iStep + 1)
    bool done = false;         // 書いた / 未対応と出した
    std::string pathLabel;     // lineImplicit / blockDPLUR / scalarDPLUR
    bool freezeTurb = false, freezeSpecies = false;
    std::string reason;        // 対象の step の中で見つけた未対応の理由 (空なら対応)
    int nCommit = 0;           // 対象の step の中の commit の回数 (1 を期待)
    std::string commitKernel;
    std::string dSource;
    Stage Q0, Qasm, R, b, d, q, qfinal;
};
State g;

bool envFlag(const char* name)
{
    const char* e = std::getenv(name);
    return e != nullptr && std::atoi(e) != 0;
}

// デバイスの flow_float の配列 [n] を読み、double に広げる (値はそのまま)。
std::vector<double> d2hWide(const flow_float* p, size_t n)
{
    std::vector<flow_float> h(n);
    if (n > 0) gpuErrchk( cudaMemcpy(h.data(), p, n*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    return std::vector<double>(h.begin(), h.end());
}

void readSet(variables& var, const char* const names[kNQ], size_t n, Stage& st)
{
    for (int k = 0; k < kNQ; ++k) st.v[k] = d2hWide(var.c_d.at(names[k]), n);
    st.got = true;
}

template <typename T>
std::vector<T> d2hFlag(const T* p, size_t n)
{
    std::vector<T> h;
    if (p == nullptr || n == 0) return h;
    h.resize(n);
    gpuErrchk( cudaMemcpy(h.data(), p, n*sizeof(T), cudaMemcpyDeviceToHost) );
    return h;
}

bool sameBits(double a, double b)
{
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

// commit の式 q == fl(b + d) のビット不一致の数 (5 量の合計)。
// float のビルドでも b + d を double で足してから float に丸めた値は float の和と同じ (double の仮数は 2·24 + 2 ビット以上)。
long long commitMismatch(const Stage& b, const Stage& d, const Stage& q, std::vector<long long>& perQ)
{
    perQ.assign(kNQ, 0);
    long long n = 0;
    for (int k = 0; k < kNQ; ++k) {
        const size_t m = q.v[k].size();
        for (size_t i = 0; i < m; ++i) {
            const flow_float s  = static_cast<flow_float>(b.v[k][i] + d.v[k][i]);
            const flow_float qq = static_cast<flow_float>(q.v[k][i]);
            if (std::memcmp(&s, &qq, sizeof(flow_float)) != 0) { ++perQ[k]; ++n; }
        }
    }
    return n;
}

long long countDiff(const Stage& a, const Stage& c, std::vector<long long>& perQ)
{
    perQ.assign(kNQ, 0);
    long long n = 0;
    for (int k = 0; k < kNQ; ++k) {
        const size_t m = a.v[k].size();
        for (size_t i = 0; i < m; ++i)
            if (!sameBits(a.v[k][i], c.v[k][i])) { ++perQ[k]; ++n; }
    }
    return n;
}

void putSet(HighFive::File& h5, const std::string& grp, const Stage& st, const std::string& when)
{
    for (int k = 0; k < kNQ; ++k) h5.createDataSet("/" + grp + "/" + kQ[k], st.v[k]);
    h5.getGroup("/" + grp).createAttribute("when", when);
}

template <typename T>
void putIfAny(HighFive::File& h5, const std::string& name, const std::vector<T>& v)
{
    if (!v.empty()) h5.createDataSet(name, v);
}

void clearStages()
{
    for (Stage* s : {&g.Q0, &g.Qasm, &g.R, &g.b, &g.d, &g.q, &g.qfinal}) {
        s->got = false;
        for (auto& x : s->v) std::vector<double>().swap(x);
    }
}

void giveUp(const std::string& why)
{
    printf("[stageDump] 未対応: step %d: %s。何も書かない (計算は続ける)\n", g.step, why.c_str());
    fflush(stdout);
    g.done = true;
    detail::g_active = false;
    detail::g_enabled = false;
    clearStages();
}

}  // namespace

void init(solverConfig& cfg, mesh& msh)
{
    const char* path = std::getenv("FORGE_DIAG_STAGE_DUMP");
    const char* stepEnv = std::getenv("FORGE_DIAG_STAGE_STEP");
    if (path == nullptr || *path == '\0') {
        if (stepEnv != nullptr && *stepEnv != '\0')
            printf("[stageDump] FORGE_DIAG_STAGE_STEP=%s は FORGE_DIAG_STAGE_DUMP が無いので無視する\n", stepEnv);
        return;
    }
    g = State{};
    g.path = path;
    if (stepEnv != nullptr && *stepEnv != '\0') {
        errno = 0;
        char* end = nullptr;
        const long n = std::strtol(stepEnv, &end, 10);
        if (errno != 0 || end == stepEnv || *end != '\0' || n < 1 || n > 1000000000L) {
            fprintf(stderr, "[stageDump] FORGE_DIAG_STAGE_STEP='%s' は 1 以上の整数でない。止める\n", stepEnv);
            std::exit(EXIT_FAILURE);
        }
        g.step = (int)n;
    }
    g.freezeTurb = envFlag("FORGE_FREEZE_TURB");
    g.freezeSpecies = envFlag("FORGE_FREEZE_SPECIES");

    // 対象は定常の陰解法の commit (advanceImplicitSteady → implicitNonlinearUpdate → applyBlock/ScalarImplicitCorrection)。
    std::string why;
    if (cfg.gpu != 1)
        why = "gpu 0 (CPU の commit にはフックが無い)";
    else if (cfg.isImplicit != 1)
        why = "陽解法 (定常の陰解法の commit だけを記録する)";
    else if (cfg.unsteady != 0)
        why = "dual-time (unsteady 1。commit は in-place の別の経路)";
    else if (cfg.qAccumulatorFP64 == 1)
        why = "qAccumulatorFP64 1 (commit の基準が roN でなく FP64 の正本)";
    else if (cfg.blockDPLUR == 1 && cfg.updateGuardAlpha > (flow_float)0.0)
        why = "updateGuardAlpha > 0 (commit が掛ける縮小率 s を写していないので d が決まらない)";
    if (!why.empty()) {
        g.done = true;
        printf("[stageDump] 未対応: %s。FORGE_DIAG_STAGE_DUMP=%s には何も書かない (計算は続ける)\n", why.c_str(), path);
        fflush(stdout);
        return;
    }
    g.pathLabel = (cfg.blockDPLUR == 1) ? ((cfg.lineImplicit == 1) ? "lineImplicit" : "blockDPLUR") : "scalarDPLUR";
    detail::g_enabled = true;
    printf("[stageDump] 有効: step %d (= iStep + 1) の 1 回の非線形更新で、ro..roe の Q0・Q_asm・R・b・d・q・q_final を "
           "double で %s に書く (経路 %s、nImplicitLines %ld、FORGE_FREEZE_TURB %d)\n",
           g.step, path, g.pathLabel.c_str(), (long)msh.nImplicitLines, g.freezeTurb ? 1 : 0);
    if (cfg.mainLoopCount() < g.step)
        printf("[stageDump] 注意: nStepOuter (%d) が対象の step %d より小さい。このままでは書かない\n", cfg.mainLoopCount(), g.step);
    fflush(stdout);
}

void beginStep(solverConfig& cfg, mesh& msh, variables& var, int iStep)
{
    (void)cfg;
    if (!detail::g_enabled || g.done || iStep + 1 != g.step) return;
    clearStages();
    g.reason.clear();
    g.nCommit = 0;
    g.commitKernel.clear();
    g.dSource.clear();
    readSet(var, kQ, (size_t)msh.nCells, g.Q0);
    detail::g_active = true;
}

void beforeSolve(solverConfig& cfg, mesh& msh, variables& var)
{
    (void)cfg;
    if (!detail::g_active) return;
    if (g.R.got) { g.reason = "1 step に線形 solve が 2 回ある"; return; }
    readSet(var, kQ, (size_t)msh.nCells, g.Qasm);
    readSet(var, kRes, (size_t)msh.nCells, g.R);
}

void afterCommit(solverConfig& cfg, mesh& msh, variables& var, const char* kernel,
                 const char* const dqNames[5], flow_float guardAlpha)
{
    (void)cfg;
    if (!detail::g_active) return;
    ++g.nCommit;
    if (g.nCommit > 1) { g.reason = "1 step に平均流の commit が 2 回ある"; return; }
    // 起動時に拒否しているが、構成の外から変わった場合に黙って違うものを書かないよう、ここでも確かめる (commitLossDiag と同じ)。
    if (var.qacc_d[0] != nullptr) { g.reason = "qAccumulatorFP64 の commit (基準が FP64 の正本)"; return; }
    if (guardAlpha > (flow_float)0.0) { g.reason = "updateGuardAlpha > 0 の commit (縮小率 s を写していない)"; return; }
    const size_t n = (size_t)msh.nCells;
    readSet(var, kN, n, g.b);
    readSet(var, dqNames, n, g.d);
    readSet(var, kQ, n, g.q);
    g.commitKernel = kernel;
    g.dSource.clear();
    for (int k = 0; k < kNQ; ++k) g.dSource += (k ? "," : "") + std::string(dqNames[k]);
}

void endStep(solverConfig& cfg, mesh& msh, variables& var)
{
    if (!detail::g_active) return;
    detail::g_active = false;
    const size_t nC = (size_t)msh.nCells;
    readSet(var, kQ, nC, g.qfinal);

    if (g.reason.empty()) {
        if (!g.R.got)      g.reason = "線形 solve の直前のフック (implicitNonlinearUpdate) を通らなかった";
        else if (!g.b.got) g.reason = "平均流の commit のフック (applyBlock/ScalarImplicitCorrection_d_wrapper) を通らなかった";
    }
    if (!g.reason.empty()) { giveUp(g.reason); return; }

    std::vector<long long> mmQ, bQ0, qfQ;
    const long long nMis = commitMismatch(g.b, g.d, g.q, mmQ);
    const long long nB = countDiff(g.b, g.Q0, bQ0);
    const long long nQf = countDiff(g.qfinal, g.q, qfQ);

    try {
        HighFive::File h5(g.path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
        h5.createAttribute("plan", std::string("plans/active/architecture-float-state-double-geometry.md §5.1 #17・§6.30"));
        h5.createAttribute("step", (long long)g.step);
        h5.createAttribute("iStep", (long long)(g.step - 1));
        h5.createAttribute("flow_float_bytes", (long long)sizeof(flow_float));
        h5.createAttribute("geom_float_bytes", (long long)sizeof(geom_float));
        h5.createAttribute("path", g.pathLabel);
        h5.createAttribute("path_note", std::string("lineImplicit = block DPLUR + line Thomas (lineImplicit 1); "
                                                    "blockDPLUR = block DPLUR without lines (the 'point' path in the code); "
                                                    "scalarDPLUR = blockDPLUR 0"));
        h5.createAttribute("commit_kernel", g.commitKernel);
        h5.createAttribute("d_source", g.dSource);
        h5.createAttribute("quantities", std::string("ro,roUx,roUy,roUz,roe"));
        h5.createAttribute("nCells", (long long)msh.nCells);
        h5.createAttribute("nCells_all", (long long)msh.nCells_all);
        h5.createAttribute("nImplicitLines", (long long)msh.nImplicitLines);
        h5.createAttribute("discretization", cfg.discretization);
        h5.createAttribute("mesh_file", cfg.meshFileName);
        h5.createAttribute("isAxisymmetric", (long long)cfg.isAxisymmetric);
        h5.createAttribute("blockDPLUR", (long long)cfg.blockDPLUR);
        h5.createAttribute("lineImplicit", (long long)cfg.lineImplicit);
        h5.createAttribute("lineDtDirectional", (long long)cfg.lineDtDirectional);
        h5.createAttribute("lineDtDirectionalCap", (double)cfg.lineDtDirectionalCap);
        h5.createAttribute("implicitThermalJacobian", (long long)cfg.implicitThermalJacobian);
        h5.createAttribute("implicitSolvePrecision", (long long)cfg.implicitSolvePrecision);
        h5.createAttribute("lowMachPrecond", (long long)cfg.lowMachPrecond);
        h5.createAttribute("nStepInner", (long long)cfg.nStepInner);
        h5.createAttribute("implicitRelax", (double)cfg.implicitRelax);
        h5.createAttribute("updateGuardAlpha", (double)cfg.updateGuardAlpha);
        h5.createAttribute("qAccumulatorFP64", (long long)cfg.qAccumulatorFP64);
        h5.createAttribute("speciesImplicitCoupling", (long long)cfg.speciesImplicitCoupling);
        h5.createAttribute("FORGE_FREEZE_TURB", (long long)(g.freezeTurb ? 1 : 0));
        h5.createAttribute("FORGE_FREEZE_SPECIES", (long long)(g.freezeSpecies ? 1 : 0));
        h5.createAttribute("n_commit_mismatch_q_vs_fl_b_plus_d", nMis);
        h5.createAttribute("n_commit_mismatch_per_qty", mmQ);
        h5.createAttribute("n_bits_diff_b_vs_Q0_per_qty", bQ0);
        h5.createAttribute("n_bits_diff_qfinal_vs_q_per_qty", qfQ);

        putSet(h5, "Q0", g.Q0, "start of the outer step (advanceImplicitSteady entry, before assembly)");
        putSet(h5, "Q_asm", g.Qasm, "ro..roe after assembly, right before the linear solve "
                                    "(state the residual/Jacobian were built from; includes in-place pins, EOS floor, BCs; not committed)");
        putSet(h5, "R", g.R, "res_ro..res_roe right before the linear solve (RHS arrays read by block DPLUR / line Thomas)");
        putSet(h5, "b", g.b, "roN..roeN read right after the commit kernel (base of ro = roN + d)");
        putSet(h5, "d", g.d, "update added by the commit kernel (" + g.dSource + "; implicitRelax already applied)");
        putSet(h5, "q", g.q, "ro..roe right after the commit kernel (before SST, species, post-commit pin, updateVariablesOuter)");
        putSet(h5, "q_final", g.qfinal, "ro..roe at the end of the outer step (after post-commit pin, updateVariablesOuter, outputs)");

        // メタ (領域の切り分け用)
        if (msh.cc64.size() >= 3*nC) putIfAny(h5, "/cells/cc64", std::vector<double>(msh.cc64.begin(), msh.cc64.begin() + 3*nC));
        putIfAny(h5, "/cells/bnode_flag",    d2hFlag(msh.bnode_flag_d, nC));
        putIfAny(h5, "/cells/axis_flag",     d2hFlag(msh.axis_flag_d, nC));
        putIfAny(h5, "/cells/wall_flag",     d2hFlag(msh.wall_flag_d, nC));
        putIfAny(h5, "/cells/iso_wall_flag", d2hFlag(msh.iso_wall_flag_d, nC));
        if (cfg.lineImplicit == 1) {
            putIfAny(h5, "/cells/line_prev", d2hFlag(msh.line_prev_d, nC));
            putIfAny(h5, "/cells/line_next", d2hFlag(msh.line_next_d, nC));
        }
    } catch (const std::exception& ex) {
        fprintf(stderr, "[stageDump] %s を書けない: %s。止める\n", g.path.c_str(), ex.what());
        std::exit(EXIT_FAILURE);
    }

    printf("[stageDump] step %d: %s に書いた (経路 %s、%s、d = %s、flow_float %zu B)。commit の照合 q == fl(b + d) の不一致 %lld、"
           "b ≠ Q0 %lld 値、q_final ≠ q %lld 値\n",
           g.step, g.path.c_str(), g.pathLabel.c_str(), g.commitKernel.c_str(), g.dSource.c_str(), sizeof(flow_float),
           nMis, nB, nQf);
    if (nMis != 0)
        printf("[stageDump] 注意: commit の照合が合わない (b・d が commit の実際の入力でない可能性)。量ごと ro %lld roUx %lld "
               "roUy %lld roUz %lld roe %lld\n", mmQ[0], mmQ[1], mmQ[2], mmQ[3], mmQ[4]);
    fflush(stdout);
    g.done = true;
    detail::g_enabled = false;
    clearStages();
}

void finalize()
{
    if (detail::g_enabled && !g.done)
        printf("[stageDump] 対象の step %d に達しなかった。%s には何も書いていない\n", g.step, g.path.c_str());
    detail::g_enabled = false;
    detail::g_active = false;
    clearStages();
}

}  // namespace stageDumpDiag
