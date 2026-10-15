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

#include <sys/resource.h>

#include <string>
#include <vector>

#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

namespace monitorOutputs {

// 試験用の口 (既定 off; plans/active/tooling-settle-convergence-and-monitor.md §5.1 #7c、result 段レビュー M2)。
// 環境変数 FORGE_TEST_OUTPUT_FSIZE="<ファイル名の正規表現 (ECMAScript, 名前全体に一致)>:<bytes>" のとき、協調停止・随時の書き出しの
// step (force) に書く出力のうち名前が一致するファイルを書く間だけ、プロセスのファイルの大きさの上限 (RLIMIT_FSIZE) を <bytes> にして
// SIGXFSZ を無視する → 上限を超える write が EFBIG で失敗し、途中まで書かれた**非空の**ファイルが残る (ディスクの満杯の模擬)。
// 例: FORGE_TEST_OUTPUT_FSIZE='res_[0-9]+\.xmf:200'。未設定・force でない・名前が一致しないときは何もしない (環境変数を 1 度読む
// だけで、出力のビットも挙動も変えない)。上限をかけている間は log に書かない (書けば落ちる; 解除のときに標準出力・標準エラーの状態を戻す)。
class TestFsizeLimit {
public:
    TestFsizeLimit(const std::string& fname, bool force);
    ~TestFsizeLimit();
    TestFsizeLimit(const TestFsizeLimit&) = delete;
    TestFsizeLimit& operator=(const TestFsizeLimit&) = delete;
private:
    bool active_ = false;
    struct rlimit old_{};
    void (*oldHandler_)(int) = nullptr;
};

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
// (全場・outputHDFflg 1 の境界の h5/xmf・節点の h5・CSV) がディスクにあって空でなく、途中で切れていないこと (h5 は HDF5 で
// 読み取りで開ける = 記録された大きさより短くない、xmf は最後の行が </Xdmf>) を確かめ、**すべて通ったときだけ**
// STOPPED_AT (JSON: step・理由・ファイルの一覧と大きさ) を一時ファイルに書いて rename する。戻り値は main の終了コード:
// 0 = 完了、3 = 確認または通知の書き込みに失敗 (STOPPED_AT は置かない)。
int finishStop(const solverConfig& cfg, const mesh& msh, bool residualStreamOk);

// 終了コードの保全 (plans/active/tooling-settle-convergence-and-monitor.md §5.1 #7c)。協調停止・随時の書き出しの step の出力の失敗を
// 記録したとき (noteStopWriteError・noteSnapshotWriteError・締めの確認の失敗) だけ、プロセスの終了時に HDF5 のライブラリの後始末を
// 飛ばして setExitStatus の値で終える (標準出力・標準エラーは flush する)。HDF5 1.10 は close に失敗した dataset の ID を残し、
// 終了時の後始末 (H5_term_library → H5D_close) がそれを閉じ直して SIGSEGV になる (2026-10-16 に確認) — 停止の失敗の終了コード 3 が
// 139 になり、随時の書き出しの失敗が run の最後の失敗になる。失敗が無い run では何もしない (後始末も今どおり)。
// main は最後の return の直前に呼ぶ (呼ぶ前に exit したら 1)。
void setExitStatus(int rc);

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
// 書き込みの失敗は例外 (HighFive の close はデストラクタで失敗を握りつぶすので、最後に flush を明示して失敗をここで出す)。
void writeNodeProbes(const solverConfig& cfg, const mesh& msh, variables& var, int iStep, bool force);
// (d) boundary_massflux.csv に 1 行。間隔は (c) と同じ。
void writeBoundaryMassflux(const solverConfig& cfg, const mesh& msh, variables& var, int iStep, bool force);

// 診断 FORGE_DIAG_OUTPUT_INVARIANCE (plans/active/tooling-settle-convergence-and-monitor.md §5.1 #1c) 専用。通常の経路からは呼ばない。
// 起動時に有効にした (c)・(d) を一時的に止める (suspend = true) / 起動時の状態に戻す (false)。
void diagSuspendMonitorKeys(bool suspend);
// 起動時に (c)・(d) のどちらかを有効にしたか。
bool monitorKeysEnabled();

}  // namespace monitorOutputs
