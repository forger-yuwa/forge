#pragma once

// 本段の監視用の出力 (plans/active/tooling-nozzle-card-local-mesh.md §4.4 1′(b)(c)(d))。
// どれも**出力だけ**で、流束・勾配・リミッタ・境界条件・陰解法・ソースには触らない。キーが無ければ何も書かない
// (協調停止は常に有効だが、STOP_REQUEST が無ければ 1 step に stat を 1 回するだけ)。
// 仕様は procedures/solver-settings.md「output」。
//
//   (b) 協調停止: 作業ディレクトリ (= run dir) の STOP_REQUEST を毎 step stat 1 回で確かめ、見つけた step で
//       全場・境界の h5・節点の出力・質量流束の行を間隔に関係なく書き、main ループを抜けて正常終了する (終了コード 0)。
//       残差履歴を書き終えた最後に、書いたものを確かめてから STOPPED_AT (JSON) を置く (= 完了の通知)。確認・通知に
//       失敗したら STOPPED_AT を置かず終了コード 3 で終える。
//   (b') 随時の書き出し: 同じく毎 step SNAPSHOT_REQUEST を stat で確かめ、見つけた step で同じ出力を**止めずに**書く。
//       書けたことを確かめて SNAPSHOT_DONE (失敗なら SNAPSHOT_FAILED) に JSON を 1 行書き足してから SNAPSHOT_REQUEST を消す。
//       失敗しても計算は続ける。
//   (c) 節点の状態: output.nodeProbeFile の節点の ρ・Ux・Uy・Uz・P・T・sonic を res_nodes_<step>.h5 に。
//   (d) 境界の質量流束: output.boundaryMassflux: 1 で boundary_massflux.csv に境界ごとの Σ massflux を 1 行。

#include <string>
#include <vector>

#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

namespace monitorOutputs {

// 起動時の検査と準備 (格子・境界条件を読んだ後、main ループの前に 1 回)。不正な入力は std::runtime_error を投げる:
//   - 作業ディレクトリに STOP_REQUEST / STOPPED_AT / SNAPSHOT_REQUEST が既にある (前の run の残り; そのまま回すと 1 step で止まる・
//     見張りが完了と取り違える)
//   - nodeProbeFile: node 以外の離散化・ファイルが開けない・数でない行・範囲外・重複・節点が 0 個
//   - boundaryMassflux: node 以外の離散化・gpu 0・軸対称の axisymMethod 1 (未対応)
void init(const solverConfig& cfg, const mesh& msh);

// 協調停止: STOP_REQUEST を stat 1 回で確かめる。見つけたら (初回だけ) 中身を読み、log に 1 行出して true。
// iStep は writeStepOutputs に渡る番号 (= 完了した step 数 = 全場のファイル名の番号)。
bool pollStopRequest(int iStep);
bool stopRequested();
int  stoppedStep();
// 停止の step の出力 (全場・境界・節点・流量) の書き込みで出た例外を記録する (writeStepOutputs が捕まえて渡す)。
void noteStopWriteError(const std::string& what);
// 協調停止の締め (main の最後、残差履歴を flush した後)。停止の step の出力の例外が無いこと、残差履歴のストリームが健全で
// 最後の行の step 列 = 停止の step − 1 であること、流量の CSV の最後の行の step = 停止の step であること、書くべきファイル
// (全場・outputHDFflg 1 の境界の h5/xmf・節点の h5・CSV) がディスクにあって空でないことを確かめ、**すべて通ったときだけ**
// STOPPED_AT (JSON: step・理由・ファイルの一覧と大きさ) を一時ファイルに書いて rename する。戻り値は main の終了コード:
// 0 = 完了、3 = 確認または通知の書き込みに失敗 (STOPPED_AT は置かない)。
int finishStop(const solverConfig& cfg, const mesh& msh, bool residualStreamOk);

// 随時の書き出し: SNAPSHOT_REQUEST を stat 1 回で確かめ、あれば中身を読んで true (この step の出力を force で書く)。
// 要求を消せなかった run では以後 false。
bool pollSnapshotRequest(int iStep);
// 随時の書き出しの step の出力で出た例外を記録する (writeStepOutputs が捕まえて渡す; 計算は止めない)。
void noteSnapshotWriteError(const std::string& what);
// 随時の書き出しの締め (同じ step の出力の直後)。出力を確かめ、SNAPSHOT_DONE (成功: step・pid・中身・ファイルの一覧と大きさ) か
// SNAPSHOT_FAILED (失敗: step・pid・中身・理由) に JSON を 1 行書き足してから SNAPSHOT_REQUEST を消す。
void finishSnapshot(const mesh& msh);
// 随時の書き出しを step k で処理したか (成功・失敗を問わず)。main がその step の残差の行を flush するのに使う
// (plans/active/tooling-settle-convergence-and-monitor.md §5.1 #1b)。
bool snapshotHandledAt(int k);
// 協調停止・随時の書き出しの step k に、全場・境界・節点・流量のほかに書いたファイル (ソルバ内 CHT の固体の状態など)。
// 締めがそれらも存在と大きさを確かめ、STOPPED_AT / SNAPSHOT_DONE の一覧に載せる (同じ step で呼ぶたびに足す)。
void noteStepFiles(int k, const std::vector<std::string>& files);

// (c) 節点の状態 res_nodes_<step>.h5。outStepInterval ごと・最後の step・force (協調停止の step)。
void writeNodeProbes(const solverConfig& cfg, const mesh& msh, variables& var, int iStep, bool force);
// (d) boundary_massflux.csv に 1 行。間隔は (c) と同じ。
void writeBoundaryMassflux(const solverConfig& cfg, const mesh& msh, variables& var, int iStep, bool force);

// 診断 FORGE_DIAG_OUTPUT_INVARIANCE (plans/active/tooling-settle-convergence-and-monitor.md §5.1 #1c) 専用。通常の経路からは呼ばない。
// 起動時に有効にした (c)・(d) を一時的に止める (suspend = true) / 起動時の状態に戻す (false)。
void diagSuspendMonitorKeys(bool suspend);
// 起動時に (c)・(d) のどちらかを有効にしたか。
bool monitorKeysEnabled();

}  // namespace monitorOutputs
