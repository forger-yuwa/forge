#include "output/outputFieldNames.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

// 出力する場の量を config output.level で絞る (procedures/solver-settings.md「output」)。
//   level 2: output_cellValNames 全部 (従来)。level 0/1: 下の基本集合 + extraFields を output_cellValNames の順で。
//   h0 (全エンタルピー) は level>=1 で合成出力 (Ht [+k]) し、属性 h0_includes_k を付ける。
// extraFields は登録済みの cell 変数なら何でも末尾に足す (extraOnly_cellValNames の wall_y_eff・dY{s}d* もこれで出る)。
// (旧 output.cpp の effectiveOutputNames。確保側 hostCellSet と共用するためここへ移した。plan architecture-solver-host-memory §4.3)
static std::list<std::string> effectiveOutputNames(const solverConfig& cfg, const variables& var)
{
    // level 2 は output_cellValNames 全部、level 0/1 は基本集合。
    // **extraFields はどの level でも効く** (2026-09-24, codex result m1): 以前は level>=2 で
    // 即 return していたため、下の「確保済み変数を出力する」処理へ到達せず、`level: 2` の run では
    // `res_ro` などを指定しても黙って出なかった。
    std::vector<std::string> base;
    if (cfg.outputLevel >= 2) {
        for (const auto& n : var.output_cellValNames) base.push_back(n);
    } else {
        for (const char* n : {"ro","roUx","roUy","roUz","roe","roK","roOmega"}) base.push_back(n);
        for (const auto& n : var.speciesVarNames) base.push_back(n);            // roY{s}
        for (const auto& n : var.condMomentConsNames) base.push_back(n);        // 凝縮モーメント保存量
        if (var.tracerRegistered != 0) base.push_back("roXi");                  // 受動トレーサ保存量 (restart 用)
        if (var.transitionRegistered != 0) { base.push_back("roGamma"); base.push_back("roReth"); }   // 遷移モデル保存量 (restart 用)
        if (cfg.outputLevel >= 1) {
            if (var.tracerRegistered != 0) base.push_back("Xi");
            if (var.transitionRegistered != 0) { for (const char* n : {"gammaTr","reTheta","gammaEff"}) base.push_back(n); }
            for (const char* n : {"P","T","Ux","Uy","Uz","k","omega","sonic","vis_lam","vis_turb","wall_dist"}) base.push_back(n);
            for (const auto& n : var.speciesVarNames) base.push_back(n.substr(2));   // Y{s}
            for (const auto& n : var.condMomentConsNames) base.push_back(n.substr(2));
        }
    }
    for (const auto& n : cfg.outputExtraFields) base.push_back(n);

    std::list<std::string> out;
    for (const auto& n : var.output_cellValNames) {
        if (std::find(base.begin(), base.end(), n) != base.end()) out.push_back(n);
    }
    // **extraFields は登録済みの cell 変数なら何でも出せる** (2026-09-24)。
    // 以前は `output_cellValNames` に入っているものしか受け付けず、`res_ro` のように
    // `cellValNames` には在って出力候補に入っていない診断量を**指定しても黙って無視**していた。
    // 丸めの内訳を場で測るのに残差そのものが要る (plan time_integration-fp64-accumulator §5.1 S6)。
    // 受付の判定は登録 (cellValNames) で行う (旧: c.count || c_d.count。c_d は実行中に operator[] でキーが増えうるので、
    // 確保の時点と出力の時点で判定が食い違わないよう登録で見る。監査 §2 (iv))。
    for (const auto& n : cfg.outputExtraFields) {
        if (std::find(out.begin(), out.end(), n) != out.end()) continue;
        if (std::find(var.output_cellValNames.begin(), var.output_cellValNames.end(), n) != var.output_cellValNames.end()) continue;
        if (std::find(var.cellValNames.begin(), var.cellValNames.end(), n) != var.cellValNames.end()) {
            out.push_back(n);   // 登録済みの診断量 (res_* など)
            continue;
        }
        static bool warned = false;
        if (!warned) { std::cerr << "[output] extraFields: '" << n << "' は確保されていない変数なので無視する\n"; warned = true; }
    }
    return out;
}

std::list<std::string> dualTimeHistoryNames(const solverConfig& cfg, const variables& var)
{
    if (!(cfg.unsteady == 1 && cfg.dualTime == 1)) return {};
    std::list<std::string> hist = {"roN", "roUxN", "roUyN", "roUzN", "roeN", "roKN", "roOmegaN"};
    for (const auto& nm : var.speciesVarNames) hist.push_back(nm + "P");
    // 受動種の履歴は scheme 1 (BDF あり) のときだけ (scheme 0 は物理時間項を持たずシフトしない = 無効な履歴; codex result M3)。
    if (cfg.passiveScalarScheme == 1) {
        if (var.tracerRegistered != 0) hist.push_back("roXiP");
        for (const auto& nm : var.condMomentConsNames) hist.push_back(nm + "P");
    }
    return hist;
}

OutputFieldPlan outputFieldPlan(const solverConfig& cfg, const variables& var)
{
    OutputFieldPlan plan;
    plan.solution = effectiveOutputNames(cfg, var);
    plan.h0 = (cfg.outputLevel >= 1) && var.c.count("Ht") && var.c.count("k");
    if (plan.h0) plan.h0Deps = {"Ht", "k"};
    plan.checkpoint = dualTimeHistoryNames(cfg, var);
    return plan;
}

std::list<std::string> initialValueNames(const variables& var)
{
    // A: 流れの保存量・壁距離・k/ω (常に)
    std::list<std::string> n = {"ro", "roUx", "roUy", "roUz", "roe", "wall_dist", "roK", "roOmega"};
    // B: 化学種 roY{s}, Y{s} (nSpecies >= 2)
    if (var.nSpeciesRegistered >= 2)
        for (int s = 0; s < var.nSpeciesRegistered; s++) { n.push_back("roY" + std::to_string(s)); n.push_back("Y" + std::to_string(s)); }
    // C: 遷移モデルの保存量 (入力に無いときは書かないが、名前は常に入れる)
    if (var.transitionRegistered != 0) { n.push_back("roGamma"); n.push_back("roReth"); }
    // D: 受動トレーサ
    if (var.tracerRegistered != 0) { n.push_back("roXi"); n.push_back("Xi"); }
    // E: 凝縮モーメントの保存量と原始量
    for (const auto& cons : var.condMomentConsNames) { n.push_back(cons); n.push_back(cons.substr(2)); }
    return n;
}

std::list<std::string> implicitDiagCellNames()
{
    return {"ro", "Ux", "Uy", "Uz", "sonic", "vis_turb", "dt_local"};
}

std::list<std::string> pinDiagResidualNames(const variables& var)
{
    std::list<std::string> names = {"scalarDirichletPin"};
    for (int k = 0; k < var.nSpeciesRegistered; ++k) names.push_back("res_roY" + std::to_string(k));
    if (var.tracerRegistered != 0) names.push_back("res_roXi");
    for (const auto& nm : var.condMomentConsNames) names.push_back("res_" + nm);
    return names;
}

std::list<std::string> pinDiagStateNames(const variables& var)
{
    std::list<std::string> names;
    for (int k = 0; k < var.nSpeciesRegistered; ++k) names.push_back("Y" + std::to_string(k));
    if (var.tracerRegistered != 0) names.push_back("Xi");
    for (const auto& nm : var.condMomentConsNames) names.push_back(nm.substr(2));
    return names;
}

std::set<std::string> hostCellSet(const solverConfig& cfg, const variables& var)
{
    std::set<std::string> H;
    // 出力 (level・extraFields・h0 の依存・/CHECKPOINT)。NaN ダンプ (dumpSolutionH5_force) も同じ関数を通る。
    const OutputFieldPlan plan = outputFieldPlan(cfg, var);
    H.insert(plan.solution.begin(), plan.solution.end());
    H.insert(plan.h0Deps.begin(), plan.h0Deps.end());
    H.insert(plan.checkpoint.begin(), plan.checkpoint.end());
    // 初期場の読込 (readValueHDF5)。dual-time の復元 (initDualTimeHistory) は plan.checkpoint と同じ名前。
    for (const auto& n : initialValueNames(var)) H.insert(n);
    // FORGE_IMPLICIT_DIAG_CSV (空でない) かつ timeIntegration 11: main.cpp の gatherImplicitDiagSnapshot
    {
        const char* path = std::getenv("FORGE_IMPLICIT_DIAG_CSV");
        if (path != nullptr && path[0] != '\0' && cfg.timeIntegration == 11)
            for (const auto& n : implicitDiagCellNames()) H.insert(n);
    }
    // FORGE_PIN_DIAG != 0 かつ node かつ dual-time: main.cpp の pinRowDiagnosticResidual / pinRowDiagnosticState
    {
        const char* e = std::getenv("FORGE_PIN_DIAG");
        if (e && std::atoi(e) != 0 && cfg.discretization == "node" && cfg.unsteady == 1 && cfg.dualTime == 1) {
            for (const auto& n : pinDiagResidualNames(var)) H.insert(n);
            for (const auto& n : pinDiagStateNames(var)) H.insert(n);
        }
    }
    // lineImplicit: main.cpp の buildImplicitLines がホストの ccx/ccy/ccz を渡す (GPU 経路では書かれない 0 のまま。監査 §5 危険 5)
    if (cfg.lineImplicit == 1) { H.insert("ccx"); H.insert("ccy"); H.insert("ccz"); }
    // 登録されていない名前は確保できないので落とす (参照されたら hostCell が名前つきで止める)
    for (auto it = H.begin(); it != H.end(); ) {
        if (std::find(var.cellValNames.begin(), var.cellValNames.end(), *it) == var.cellValNames.end()) it = H.erase(it);
        else ++it;
    }
    return H;
}

void registerOutputDiagnostics(variables& var)
{
    // 診断 (FORGE_OUT_RESIDUALS=1): 流れ残差場と陰的補正 dq を h5 出力へ追加する
    // (サブ反復収縮の空間局在の測定用。既定 off = 出力不変)。書かれる値は「最終サブ反復・
    // 最終 sweep 時点」の res_* (BDF 項込み R*) と dq_block_new_* (implicitRelax 適用後)。
    // 注意: blockDPLURSolve は sweep 毎に new/old を swap するため、最終補正は dq_block_old_* に
    // 残る (dq_block_new_* は 1 sweep 前 — 2026-09-03 Codex 指摘で修正)。
    // (旧 main.cpp の main() にあった登録。確保より前へ移した。条件は同じ。level < 2 では出力に効かない・level 2 では
    //  output_cellValNames と重複する既存の挙動もそのまま。監査 §5 危険 8)
    if (const char* e = std::getenv("FORGE_OUT_RESIDUALS"); e && std::atoi(e) != 0) {
        for (const char* n : {"res_ro","res_roUx","res_roUy","res_roUz","res_roe",
                              "dq_block_old_0","dq_block_old_1","dq_block_old_2",
                              "dq_block_old_3","dq_block_old_4"})
            var.output_cellValNames.push_back(n);
        std::printf("[FORGE_OUT_RESIDUALS] residual/dq fields added to h5 outputs\n");
    }
    // FORGE_RESID_SNAP=1: dual-time の subiter 0 直後の res/dq を未使用スロット (res_*_m / dq_*_new
    // スカラー枠) へ退避して出力に含める → 局所収縮率 g=|dq_final|/|dq_sub0| を場で測れる。
    if (const char* e = std::getenv("FORGE_RESID_SNAP"); e != nullptr) {  // 値は退避 subiter 番号 ("0" も有効)
        for (const char* n : {"res_ro_m","res_roUx_m","res_roUy_m","res_roUz_m","res_roe_m",
                              "dq_ro_new","dq_roUx_new","dq_roUy_new","dq_roUz_new","dq_roe_new"})
            var.output_cellValNames.push_back(n);
        std::printf("[FORGE_RESID_SNAP] subiter-0 residual/dq snapshots added to h5 outputs\n");
    }
}

void applyEnvGatedRemovals(variables& var)
{
    // (旧 variables::allocVariables の先頭。確保集合 H を作る前に外すため切り出した。条件・順序は同じ)
    auto removeAll = [&](const char* n) {
        var.cellValNames.remove(n);
        var.output_cellValNames.remove(n);
        var.c.erase(n);
        var.c_d.erase(n);
    };
    // W-I 実力診断 (§4.2) は FORGE_WI_FORCE_DIAG=1 のときだけ有効。OFF なら
    // 確保も出力もしない (通常経路のメモリ・D2H・HDF5 を増やさない)。
    {
        const char* e = std::getenv("FORGE_WI_FORCE_DIAG");
        if (!(e && std::atoi(e) != 0)) {
            for (const char* n : {"wi_ftan", "wi_fnrm", "wi_fnrm_abs", "wi_ftan_res", "wi_eheat", "wi_ework"}) removeAll(n);
        }
    }
    // E3 (§5) の wf_sprod も env ゲート (OFF ならメモリ・D2H・HDF5 を増やさない)
    {
        const char* e = std::getenv("FORGE_WF_OMEGA_SOURCE");
        if (!(e && std::atoi(e) != 0)) removeAll("wf_sprod");
    }
    // 閉包則診断 (§5.2 ③) の wf_g も env ゲート
    {
        const char* e = std::getenv("FORGE_WF_CLOSURE_DIAG");
        if (!(e && std::atoi(e) != 0)) removeAll("wf_g");
    }
    // 代表点幾何診断 (§3.1) も env ゲート
    {
        const char* e = std::getenv("FORGE_WF_REP_DIAG");
        if (!(e && std::atoi(e) != 0)) {
            for (const char* n : {"rep_id", "rep_y", "rep_dist", "rep_cos", "rep_toff", "rep_wdratio",
                                  "rep_nx", "rep_ny", "rep_nz"}) removeAll(n);
        }
    }
    // omega 項別収支 (§4.1) も env ゲート
    {
        const char* e = std::getenv("FORGE_OMEGA_BUDGET");
        if (!(e && std::atoi(e) != 0)) {
            for (const char* n : {"omg_prod", "omg_dest", "omg_cross", "omg_trans", "omg_axisym"}) removeAll(n);
        }
    }
}
