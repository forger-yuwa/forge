// 診断 FORGE_DIAG_QUANTIZE_STATE (plans/active/architecture-float-state-double-geometry.md §5.1 #19・#20・§6.36・§6.40、既定 off)。
// 対象・時点・拒否・自己点検の書式は quantizeStateDiag.hpp。ここはホスト側の段取り (環境変数・拒否・カーネルの呼び出し) と
// 自己点検の読み取り (D2H) と集計。1 値の処理と集計・書式・名前の検査の関数は quantizeStateDiag_d.cuh (単体試験と共用)。

#include "cuda_forge/quantizeStateDiag.hpp"
#include "cuda_forge/quantizeStateDiag_d.cuh"
#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <string>
#include <vector>

namespace quantizeStateDiag {

namespace detail {
bool g_enabled = false;
bool g_speciesPending = false;
}

namespace {

struct State {
    std::vector<std::string> targets;     // 丸める保存量 (Q)。FORGE_DIAG_QUANTIZE_STATE_ONLY があれば選んだ量だけ
    std::vector<std::string> histN;       // 点検する履歴 Q_N (丸める量の分だけ)
    std::vector<std::string> histM;       // 点検する履歴 Q_M (丸める量の分だけ)
    std::vector<std::string> speciesN;    // (b) roY{s}N (丸める化学種の分だけ)
    std::vector<std::string> speciesM;    // (b) roY{s}M (丸める化学種の分だけ)
    bool only = false;                    // FORGE_DIAG_QUANTIZE_STATE_ONLY で丸める量を選んだ
    std::vector<std::string> unselected;  // 選ばなかった対象 (丸めない。自己点検の step で丸めの前後のビットを比べる)
    std::set<int> checkSteps;             // 自己点検の step (iStep + 1)
    std::set<int> checked;                // (a) を出した step
    std::set<int> unselChecked;           // 選ばなかった量の丸め前後の比較を出した step
    std::set<int> pendingSpecies;         // (a) を出し、(b) を待っている step
};
State g;

std::string joinNames(const std::vector<std::string>& v)
{
    std::string t;
    for (const auto& s : v) t += (t.empty() ? "" : ", ") + s;
    return t;
}

std::string stepsText(const std::set<int>& st)
{
    std::string t;
    for (int v : st) t += (t.empty() ? "" : ",") + std::to_string(v);
    return t;
}

// "1,30000" を読む。不正なら止める (黙って違う step を点検しない)。
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
            fprintf(stderr, "[quantizeState] FORGE_DIAG_QUANTIZE_STATE_CHECK='%s' の '%s' は 1 以上の整数でない。止める\n",
                    s, tok.c_str());
            std::exit(EXIT_FAILURE);
        }
        out.insert((int)n);
        if (comma == std::string::npos) break;
        pos = comma + 1;
    }
    return out;
}

[[noreturn]] void reject(const std::string& why)
{
    fprintf(stderr, "[quantizeState] 拒否: %s\n", why.c_str());
    fprintf(stderr, "[quantizeState] FORGE_DIAG_QUANTIZE_STATE=1 の対応: FP64 のビルド・GPU・node・定常の陰解法 "
                    "(timeIntegration 11・unsteady 0)・qAccumulatorFP64 0\n");
    std::exit(EXIT_FAILURE);
}

// デバイスの配列 [n] を host に写す (値はそのまま)。
std::vector<flow_float> d2h(variables& var, const std::string& name, size_t n)
{
    std::vector<flow_float> h(n);
    if (n > 0) gpuErrchk( cudaMemcpy(h.data(), var.c_d.at(name), n*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    return h;
}

// 群ごとに全 CV と内点 (ic < nCells) の表せない値を数え、1 行出す。0 でなければ内点の行も出し、node の入口のピン
// (scalarDirichletPin、組立の中で ρY・ρk・ρω を Dirichlet 値に書く。speciesTransport_d.cu の species_dirichlet_boundary_d)
// があれば、それを除いた内点の個数も続ける。空の群 (FORGE_DIAG_QUANTIZE_STATE_ONLY で化学種だけを選んだときの Q_N 等) は
// 行に出さない (未設定では空の群は無いので、行は #19 のまま)。
void countAndPrint(int step, const std::vector<std::pair<std::string, const std::vector<std::string>*>>& spec,
                   mesh& msh, variables& var, const std::string& when)
{
    const size_t nAll = (size_t)msh.nCells_all, nInt = (size_t)msh.nCells;
    std::vector<flow_float> pin;
    {
        auto it = var.c_d.find("scalarDirichletPin");
        if (it != var.c_d.end() && it->second != nullptr) pin = d2h(var, "scalarDirichletPin", nInt);
    }
    long long nPin = -1;
    if (!pin.empty()) { nPin = 0; for (flow_float p : pin) nPin += (p == (flow_float)1.0); }
    std::vector<QuantizeStateCountGroup> all, interior, noPin;
    for (const auto& sp : spec) {
        if (sp.second->empty()) continue;
        QuantizeStateCountGroup ga{sp.first, *sp.second, {}}, gi{sp.first, *sp.second, {}}, gp{sp.first, *sp.second, {}};
        for (const auto& name : *sp.second) {
            const std::vector<flow_float> h = d2h(var, name, nAll);
            ga.counts.push_back(quantizeStateCountNotRepresentable(h, 0, nAll));
            gi.counts.push_back(quantizeStateCountNotRepresentable(h, 0, nInt));
            gp.counts.push_back(pin.empty() ? 0 : quantizeStateCountNotRepresentable(h, 0, nInt, &pin));
        }
        all.push_back(ga);
        interior.push_back(gi);
        noPin.push_back(gp);
    }
    long long total = 0, totalInt = 0, totalNoPin = 0;
    const std::string body = quantizeStateFormat(all, &total);
    const std::string bodyInt = quantizeStateFormat(interior, &totalInt);
    const std::string bodyNoPin = quantizeStateFormat(noPin, &totalNoPin);
    printf("%s\n", quantizeStateCheckLine(step, body, (long long)msh.nCells_all, when).c_str());
    if (total != 0) printf("%s\n", quantizeStateInteriorLine(step, (long long)msh.nCells, bodyInt, nPin, bodyNoPin).c_str());
    fflush(stdout);
}

}  // namespace

void init(solverConfig& cfg, mesh& msh, variables& var, bool sst)
{
    const char* env = std::getenv("FORGE_DIAG_QUANTIZE_STATE");
    const char* checkEnv = std::getenv("FORGE_DIAG_QUANTIZE_STATE_CHECK");
    const char* onlyEnv = std::getenv("FORGE_DIAG_QUANTIZE_STATE_ONLY");
    const bool on = (env != nullptr && std::strcmp(env, "1") == 0);
    if (env != nullptr && *env != '\0' && !on && std::strcmp(env, "0") != 0) {
        fprintf(stderr, "[quantizeState] FORGE_DIAG_QUANTIZE_STATE='%s' は 1 (有効) か 0 (無効) でない。止める\n", env);
        std::exit(EXIT_FAILURE);
    }
    if (!on) {
        if (checkEnv != nullptr && *checkEnv != '\0')
            printf("[quantizeState] FORGE_DIAG_QUANTIZE_STATE_CHECK=%s は FORGE_DIAG_QUANTIZE_STATE=1 が無いので無視する\n",
                   checkEnv);
        if (onlyEnv != nullptr && *onlyEnv != '\0')
            printf("[quantizeState] FORGE_DIAG_QUANTIZE_STATE_ONLY=%s は FORGE_DIAG_QUANTIZE_STATE=1 が無いので無視する\n",
                   onlyEnv);
        return;
    }

    // 拒否。どれも「黙って違う試験にしない」ため。
    if (sizeof(flow_float) != 8)
        reject("float のビルド (sizeof(flow_float) = " + std::to_string(sizeof(flow_float)) +
               ")。状態がもともと float32 なので丸めても何も変わらない (FP64 のビルド専用)");
    if (cfg.gpu != 1)
        reject("gpu 0 (丸めはデバイスの配列に当てるカーネル。CPU の経路にはフックが無い)");
    if (cfg.discretization != "node")
        reject("discretization '" + cfg.discretization + "' (node のみ。cell は使わない)");
    if (cfg.timeIntegration != 11 || cfg.unsteady != 0)
        reject("timeIntegration " + std::to_string(cfg.timeIntegration) + "・unsteady " + std::to_string(cfg.unsteady) +
               " (定常の陰解法 timeIntegration 11・unsteady 0 のみ。フックは advanceImplicitSteady にだけある)");
    if (cfg.qAccumulatorFP64 == 1)
        reject("qAccumulatorFP64 1 との併用 (FP64 の正本と reconcile が状態の積み方を変えるので、末尾の量子化だけの試験にならない)");

    g = State{};
    g.targets = {"ro", "roUx", "roUy", "roUz", "roe"};
    g.histN   = {"roN", "roUxN", "roUyN", "roUzN", "roeN"};
    g.histM   = {"roM", "roUxM", "roUyM", "roUzM", "roeM"};
    if (sst) {
        for (const char* n : {"roK", "roOmega"}) g.targets.emplace_back(n);
        for (const char* n : {"roKN", "roOmegaN"}) g.histN.emplace_back(n);
        for (const char* n : {"roKM", "roOmegaM"}) g.histM.emplace_back(n);
    }
    if (var.nSpeciesRegistered >= 2) {
        for (const auto& n : var.speciesVarNames) {
            g.targets.push_back(n);
            g.speciesN.push_back(n + "N");
            g.speciesM.push_back(n + "M");
        }
    }
    // 名前の配列がデバイスにあることを確かめる (無ければ内部の不整合なので止める)。選ぶ前の全量で見る (選ばなかった対象も
    // 丸めの前後の比較で読む)。
    for (const auto* list : {&g.targets, &g.histN, &g.histM, &g.speciesN, &g.speciesM}) {
        for (const auto& n : *list) {
            auto it = var.c_d.find(n);
            if (it == var.c_d.end() || it->second == nullptr) {
                fprintf(stderr, "[quantizeState] 内部エラー: デバイスの配列 %s が無い。止める\n", n.c_str());
                std::exit(EXIT_FAILURE);
            }
        }
    }
    g.checkSteps = (checkEnv != nullptr && *checkEnv != '\0') ? parseSteps(checkEnv) : std::set<int>{1};

    // FORGE_DIAG_QUANTIZE_STATE_ONLY (§5.1 #20): 丸める量を選ぶ。設定されていれば (空の値も) 検査し、対象に無い名前・空の要素・
    // 重複は止める。選んだ量に合わせて、自己点検 (a) の履歴と (b) の化学種の履歴も選んだ量の分だけにする。
    std::string what = joinNames(g.targets);   // 起動時の行の「丸める量」の欄 (未設定なら #19 のまま)
    if (onlyEnv != nullptr) {
        std::vector<bool> sel;
        std::string err;
        if (!quantizeStateParseOnly(onlyEnv, g.targets, &sel, &err)) {
            fprintf(stderr, "[quantizeState] FORGE_DIAG_QUANTIZE_STATE_ONLY='%s': %s。この構成の丸めの対象は %s (カンマ区切り、"
                            "空白なし)。止める\n", onlyEnv, err.c_str(), joinNames(g.targets).c_str());
            std::exit(EXIT_FAILURE);
        }
        what = quantizeStateOnlyText(g.targets, sel);
        // targets の並びは [ro..roe (、roK・roOmega)] が histN・histM と、[roY0..] が speciesN・speciesM と同じ添字。
        const size_t nFlow = g.histN.size();
        if (g.histM.size() != nFlow || g.speciesN.size() != g.speciesM.size() ||
            g.targets.size() != nFlow + g.speciesN.size() || sel.size() != g.targets.size()) {
            fprintf(stderr, "[quantizeState] 内部エラー: 対象 %zu・履歴 %zu/%zu・化学種の履歴 %zu/%zu の数が合わない。止める\n",
                    g.targets.size(), g.histN.size(), g.histM.size(), g.speciesN.size(), g.speciesM.size());
            std::exit(EXIT_FAILURE);
        }
        State f;
        for (size_t k = 0; k < g.targets.size(); ++k) {
            if (!sel[k]) { f.unselected.push_back(g.targets[k]); continue; }
            f.targets.push_back(g.targets[k]);
            if (k < nFlow) {
                f.histN.push_back(g.histN[k]);
                f.histM.push_back(g.histM[k]);
            } else {
                f.speciesN.push_back(g.speciesN[k - nFlow]);
                f.speciesM.push_back(g.speciesM[k - nFlow]);
            }
        }
        g.targets.swap(f.targets);
        g.histN.swap(f.histN);
        g.histM.swap(f.histM);
        g.speciesN.swap(f.speciesN);
        g.speciesM.swap(f.speciesM);
        g.unselected.swap(f.unselected);
        g.only = true;
    }

    detail::g_enabled = true;
    detail::g_speciesPending = false;
    printf("[quantizeState] 診断用: step の末尾で保存量 %zu 本を float32 に丸める (%s)。全 CV (nCells_all %ld)、"
           "時点は advanceImplicitSteady の commit 後の等温壁のピンの後・updateVariablesOuter の前。自己点検 step %s\n",
           g.targets.size(), what.c_str(), (long)msh.nCells_all, stepsText(g.checkSteps).c_str());
    if (g.only)
        printf("[quantizeState] 自己点検の step では、表せない値は選んだ量とその履歴だけを数え、選ばなかった量は丸めのカーネルの"
               "直前と直後で生のビットを比べる (変わった値 0 が期待)\n");
    // 丸めない保存量を黙らせない (拒否はしない。§5.1 #19 の対象は上の量だけ)。
    std::vector<std::string> notQ;
    if (var.transitionRegistered != 0) notQ.push_back("遷移モデル (roGamma, roReth)");
    if (var.tracerRegistered != 0)     notQ.push_back("受動トレーサ (roXi)");
    if (!var.condMomentConsNames.empty()) notQ.push_back("凝縮モーメント (" + joinNames(var.condMomentConsNames) + ")");
    if (!notQ.empty()) {
        printf("[quantizeState] 注意: 次の保存量は丸めない (double のまま次の step へ持ち越す):");
        for (size_t i = 0; i < notQ.size(); ++i) printf("%s %s", i ? " /" : "", notQ[i].c_str());
        printf("\n");
    }
    if (cfg.mainLoopCount() < *g.checkSteps.rbegin())
        printf("[quantizeState] 注意: nStepOuter (%d) が自己点検の最大の step %d より小さい。その step は点検しない\n",
               cfg.mainLoopCount(), *g.checkSteps.rbegin());
    fflush(stdout);
}

void quantize(cudaConfig& cuda_cfg, mesh& msh, variables& var, int iStep)
{
    if (!detail::g_enabled) return;
    // 量を選んだとき (FORGE_DIAG_QUANTIZE_STATE_ONLY) の自己点検の step: 選ばなかった対象を丸めのカーネルの直前と直後に
    // host に写し、生のビットが変わった値の個数を量ごとに出す (カーネルは選んだ配列にしか書かないので 0 が期待。
    // 配列の取り違え・別名を見る)。未設定では写さない (#19 のまま)。
    const int step = iStep + 1;
    const bool cmp = g.only && g.checkSteps.count(step) != 0 && g.unselChecked.count(step) == 0;
    const size_t nAll = (size_t)msh.nCells_all;
    std::vector<std::vector<flow_float>> before;
    if (cmp) for (const auto& n : g.unselected) before.push_back(d2h(var, n, nAll));
    for (const auto& n : g.targets) quantizeStateToFloat32_d_wrapper(cuda_cfg, msh.nCells_all, var.c_d.at(n));
    if (cmp) {
        std::vector<long long> changed;
        for (size_t k = 0; k < g.unselected.size(); ++k)
            changed.push_back(quantizeStateCountChanged(before[k], d2h(var, g.unselected[k], nAll)));
        printf("%s\n", quantizeStateUnselectedLine(step, g.unselected, changed, (long long)msh.nCells_all, nullptr).c_str());
        fflush(stdout);
        g.unselChecked.insert(step);
    }
}

void afterUpdateOuter(mesh& msh, variables& var, int iStep)
{
    if (!detail::g_enabled) return;
    const int step = iStep + 1;
    if (g.checkSteps.count(step) == 0 || g.checked.count(step) != 0) return;
    countAndPrint(step, {{"Q", &g.targets}, {"Q_N", &g.histN}, {"Q_M", &g.histM}}, msh, var, quantizeStateWhenOuter());
    g.checked.insert(step);
    if (!g.speciesN.empty()) {
        g.pendingSpecies.insert(step);
        detail::g_speciesPending = true;
    }
}

void afterSpeciesUpdateOuter(mesh& msh, variables& var, int iStep)
{
    if (!detail::g_enabled || !detail::g_speciesPending) return;
    const int stepNow = iStep + 1;
    // 対象の step の次の step の化学種の更新で数える。それより前の step のものが残っていたら (化学種の更新を通らなかった)、
    // 数えずにその旨を出す (後の step の roY{s}N は対象の step の末尾の ρY から取ったものでない)。
    std::set<int> rest;
    for (int st : g.pendingSpecies) {
        if (st == stepNow - 1) {
            countAndPrint(st, {{"roY_N", &g.speciesN}, {"roY_M", &g.speciesM}}, msh, var, quantizeStateWhenSpecies(stepNow));
        } else if (st < stepNow - 1) {
            printf("[quantizeState] step %d: 未対応: step %d で化学種の更新を通らなかった。roY_N・roY_M は数えない\n", st, st + 1);
            fflush(stdout);
        } else {
            rest.insert(st);
        }
    }
    g.pendingSpecies.swap(rest);
    detail::g_speciesPending = !g.pendingSpecies.empty();
}

void finalize()
{
    if (detail::g_enabled) {
        std::set<int> missing;
        for (int s : g.checkSteps) if (g.checked.count(s) == 0) missing.insert(s);
        if (!missing.empty())
            printf("[quantizeState] 自己点検の step %s に達しなかった。その分は数えていない\n", stepsText(missing).c_str());
        if (!g.pendingSpecies.empty())
            printf("[quantizeState] step %s の roY_N・roY_M は、次の step の化学種の更新に達しなかったので数えていない\n",
                   stepsText(g.pendingSpecies).c_str());
    }
    detail::g_enabled = false;
    detail::g_speciesPending = false;
}

}  // namespace quantizeStateDiag
