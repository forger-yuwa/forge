# 諮問: V0 後半の切り分け — 「加算順を固定した 1 更新 A/B」の最小の実装設計 (2026-10-03)

関連 plan: `plans/active/boundary-node-farfield-characteristic.md` §5.1 #5a (V0 後半 FAIL、ユーザ決定 2026-10-03「切り分ける」)。前回諮問: `notes/reviews/2026-10-03-farfield-v0b-result-diagnose.md` (推奨 = 旧・新の同一入力で加算順を固定した完全な 1 outer step の A/B、同一版内のビット一致が前提ゲート)。

## 事実
- 旧 = `~/sglsq/forge_2fa3826c` (commit 2fa3826c)、新 = build-ff (commit a6ceee0b 系、sha b0240cd7)。`git diff --stat 2fa3826c a6ceee0b -- solver_density_cuda` は 25 ファイル。farfield 以外のソルバ変更: `boundaryCond.{cpp,hpp}`、`convectiveFlux_common_d.cuh` (+11、FORGE_DIAG_FACE_VEL_CELL)、`convectiveFlux_d.cu` (+283、farfield 振り分け・ダンプ・帳簿)、`convectiveFlux_slau_d.inc.cuh` (+23、FACE_VEL_CELL 介入)、`passiveKernels_d.cuh` (9)、`ransBoundary_d.cu` (6)、`ransTransport_d.cu` (2)、`scalarTransport_d.{cu,cuh}` (16+4、farfield 面値の読み分け)、`speciesTransport_d.cu` (7)、`solverConfig.hpp` (2)、`main.cpp` (11)、`probe/point_probes.cu` (1)。
- `cuda_forge` の `atomicAdd` は 36 ファイル 478 か所、うち `res_*` への加算は 162 か所。node 残差・輸送対角・勾配の多くが atomicAdd。
- 既存の診断: `FORGE_DUMP_MASSFLUX` (初回評価の massflux と ro,Ux,Uy,Uz,P,sonic)、`FORGE_DUMP_LEDGER` (段別残差の帳簿、節点・変数ごと)。決定的な加算の仕組みは無い (前回諮問で確認)。
- 実行は AWS g5 (A10G)。SERN g3 は約 100 万節点級? (run_0971)。

## 設計候補 (どれも未実装)
A. **固定小数点の影アキュムレータ**: 診断ビルド (`-DFORGE_DET_ATOMIC`) で `atomicAdd(float*, float)` をラッパーに置き換え、値を int64 固定小数点に変換して影配列へ整数 atomicAdd (結合則が成り立つので順序非依存)、段の終わりに float へ戻す。宛先ポインタ → 影配列の対応表が要る。両版 (旧・新のソース) に同じパッチを当ててビルド。
B. **倍精度アキュムレータ + 最後に float へ丸め**: 影配列を double にし atomicAdd(double)。順序依存は残るが、float への最終丸めでほぼ消える (丸め境界に当たる値だけ非決定的)。ビット一致の保証はないが確率的にほぼ決定的。
C. **静的な切り分け**: 旧→新の非 farfield 差分を行単位で確認し、farfield の無い構成で実行経路に入る変更が無いこと (FACE_VEL_CELL は env 未設定で -1、面値ポインタは nullptr 経路) を証明し、加えて新バイナリで farfield コードを `#if 0` で除いたビルドとの 1 step 比較 (これも非決定的なので不十分?)。
D. **部分決定化**: 差分のあるカーネル (scalarTransport・speciesTransport・passive・rans) と流れの残差だけを A/B で決定化し、他は対照で非決定のまま。
E. **同一ソースの対照**: 新ソースから farfield 追加分だけを除いたビルド (= 旧と同じ作用素のはず) を作り、多 step で旧と同じ統計か比べる (統計比較なので推奨外?)。

## 問い
1. 前回の推奨を満たす最小の設計はどれか (A〜E 以外も可)。特に「同一版内のビット一致」を 1 outer step (内反復 5) で得る現実的な方法。
2. 旧ソースにも同じ診断パッチを当てる必要があるか (旧バイナリは改変できないので再ビルドになる。再ビルド旧と配布済み旧の同一性はどう担保するか)。
3. 合否 (事前閾値) と、A (一致) / B (不一致) それぞれの次の一手。V0 合格に昇格できる条件はあるか (前回は「1 更新一致は V0 合格に昇格させない」)。
