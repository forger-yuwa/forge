#pragma once

// 出力・checkpoint・初期場の読込・環境変数で有効になる診断が「ホストで読み書きするセル変数の名前」を 1 か所で作る
// (plan architecture-solver-host-memory §4.3、監査 notes/investigations/2026-10-07-host-memory-audit.md §3)。
// 出力側 (output.cpp)・確保側 (allocVariables に渡すホスト確保集合 H)・dual-time の復元 (main.cpp の
// initDualTimeHistory)・診断 (main.cpp) が同じ関数を使うので、出力を変えても確保と食い違わない
// (食い違ったら variables::hostCell が名前つきで止める)。

#include <list>
#include <set>
#include <string>

#include "input/solverConfig.hpp"
#include "variables.hpp"

// writeSolutionH5_XDMF が書く名前の計画。
struct OutputFieldPlan {
    std::list<std::string> solution;     // /VALUE と XDMF の Attribute の順 (重複も今のまま残す。監査 §5 危険 8)
    bool h0 = false;                      // /VALUE/h0 (全エンタルピー) を合成して書くか (level >= 1)
    std::list<std::string> h0Deps;        // h0 を組むために D2H する名前 (Ht, k)
    std::list<std::string> checkpoint;    // dual-time の /CHECKPOINT 履歴 (= dualTimeHistoryNames)
};

// 出力名の計画。入力は cfg (output.level / extraFields / dual-time) と登録済みの変数 (registerOutputDiagnostics・
// applyEnvGatedRemovals を済ませた後)。extraFields の受付は「登録済み (cellValNames)」で判定する。
OutputFieldPlan outputFieldPlan(const solverConfig& cfg, const variables& var);

// dual-time の履歴 (前物理レベル): 流れ *N、化学種 roY{s}P、受動種 <cons>P (passiveScalarScheme 1 のとき)。
// unsteady == 1 かつ dualTime == 1 でなければ空。/CHECKPOINT の書出し (output.cpp) と復元 (main.cpp) が共用する。
std::list<std::string> dualTimeHistoryNames(const solverConfig& cfg, const variables& var);

// readValueHDF5 がホストに書いて (または読んで) H2D する名前 (監査 §1 の A〜E)。
std::list<std::string> initialValueNames(const variables& var);

// FORGE_IMPLICIT_DIAG_CSV の診断 (main.cpp の gatherImplicitDiagSnapshot) が D2H して読む名前。
std::list<std::string> implicitDiagCellNames();
// FORGE_PIN_DIAG の診断 (main.cpp の pinRowDiagnosticResidual / pinRowDiagnosticState) が D2H して読む名前。
std::list<std::string> pinDiagResidualNames(const variables& var);   // 先頭は scalarDirichletPin
std::list<std::string> pinDiagStateNames(const variables& var);

// ホスト確保集合 H: GPU 経路でホスト側のセル配列を nCells_all 長で確保する名前
// (上の出力・checkpoint・初期場の和 ∪ 有効な診断の名前 ∪ lineImplicit の ccx/ccy/ccz)。登録されていない名前は含めない。
std::set<std::string> hostCellSet(const solverConfig& cfg, const variables& var);

// FORGE_OUT_RESIDUALS / FORGE_RESID_SNAP で出力名 (output_cellValNames) を足す。確保 (hostCellSet) より前に呼ぶ。
void registerOutputDiagnostics(variables& var);

// 環境変数で有効になる診断変数 (wi_* / wf_sprod / wf_g / rep_* / omg_*) を、無効なら登録ごと外す
// (cellValNames / output_cellValNames / c / c_d)。確保 (hostCellSet・allocVariables) より前に呼ぶ。
void applyEnvGatedRemovals(variables& var);
