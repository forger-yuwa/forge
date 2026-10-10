# 熱伝導の近傍 K を入れたライン粘性 Jacobian の破綻は、面エンタルピーの float の評価が要るか (粘性ヤコビアン plan の再開)

## メタ

- **area**: `time_integration`
- **status**: `in_progress`
- **related_docs**:
  - `methods/time_integration/implementation.md` の `lineViscCoupling: 2`・`3` の節 (診断用、本線不採用)
- **related_plans**:
  - [`../accepted/time_integration-line-viscous-jacobian.md`](../accepted/time_integration-line-viscous-jacobian.md) (親。§6.16 の面エンタルピーの A/B、§6.17 の再開の条件)
  - [`time_integration-implicit-thermal-jacobian.md`](time_integration-implicit-thermal-jacobian.md) §6.2・§6.3 (point 仕上げでの面エンタルピーの A/B、NOT_SUPPORT)
  - [`time_integration-line-implicit-speed.md`](time_integration-line-implicit-speed.md) (本線 B0 = 値 0 + 方向別 + キー 5 + 上限 50)
- **created**: `2026-10-10`
- **owner**: Claude (ユーザ 2026-10-10「速度上がる上がらないでなく、筋のいい手法がいいと思ってる」、マスク 7 + 面エンタルピーの double の A/B の提案に「よい」)

## 1. 目的

ライン面に薄層の粘性・熱伝導の Jacobian を入れる値 2・3 は、既往の探索では、熱伝導の近傍 K を入れる (マスク 7、設計どおりの全部入り) と case/45 の方向別 dt で発散し、
熱伝導の近傍 K を外す (マスク 5) と 2000 step 回った。ただしこの 7/5 の比較は事前確認が許容外で無効 (探索的な観測) なので、熱伝導の K が破綻に要るという帰属は確定していない。
マスク 5 は温度をライン内で結合しないので、親 plan の動機 (壁際のエントロピーのモードをライン内で結合する) を満たさない回避策である。
親 plan §6.16 では、面エンタルピーを float で評価すると、エネルギーの行の実残差の方向微分が double の約 28 倍に見えた (‖J_t p‖ 148.9 → 5.28、S0・p7・ライン 2183・ε = 1e-6 に限る)。
そこで「値 3・マスク 7 の破綻 (非有限化) が、面エンタルピーの精度の切替だけで避けられるか」を、診断の切替だけを変える A/B で確かめる。
支持されれば、筋のいい形 (熱伝導も結合した全部入り) を残差の評価の精度とセットで成り立たせる道を調べる意味がある。棄却されれば「この切替だけでは回避できない」まで分かる。
この A/B が測るのは、切替による 2000 step 以内の非有限化の回避までである。精度依存と近似 Jacobian の問題は排他的でなく併存しうるので、両側が発散しても熱伝導の K に原因を絞れない (codex 諮問 2026-10-10)。

## 2. スコープ

- **やる**: 既存のバイナリと診断の切替 (`lineViscCoupling` 2・3、`FORGE_LVC_TERMS`、`FORGE_DIAG_FACE_H_DOUBLE`) だけで組む短い A/B (§6)。コードは変えない。
- **やらない**: 面エンタルピーの double の既定化、マスク 5 の設定キー化、総時間・速度の評価、float ビルドでの試験 (いずれも本 A/B の結果を見て別に決める)。

## 3. 関連 docs と前提

- 親 plan の観測 (同じ出発点 = `run_0183` の res_100000、キー 5・方向別・上限なし・cfl 4・緩和 0.7):
  - 値 2 (マスク 7): `run_0261`・`run_0300` とも 29 step で非有限。
  - 値 3・マスク 7: `run_0301` 122 step、`run_0305` 150 step、`run_0311` 122 step で非有限 (同じ入力でも再実行で破綻の step が変わる)。
  - 値 3・マスク 5: `run_0306` は 2000 step 有限 (ただし事前確認が許容外で比較は無効・探索的な観測、親 §6.9)。値 3・マスク 15 (密度の列だけ外す): `run_0312` 566 step で非有限。
  - これらは別のバイナリ (lineE〜lineG、従来の配列の並び)。今回の lineM_fp64 は LAYOUT2 が既定で、値 3・マスク 5 の経路は従来の並びとビット一致を確かめてある (速度 plan §6.21、`run_0380`)。
- `FORGE_DIAG_FACE_H_DOUBLE=1` は SLAU の TP 多成分の面エンタルピー (`convectiveFlux_slau_d.inc.cuh` の `thermo_h_mix_f`) を、同じ NASA 係数・datum の double 版で評価する。既定の経路は FP64 ビルドでも T・Y を float にして float を返す。
- point 仕上げの残差の床には、この切替は 10 % 以上効かなかった (thermal-jacobian plan §6.3、NOT_SUPPORT)。今回は別の問い (ライン陰解法の破綻) である。

## 4. 設計方針

- 変えるのは環境変数 `FORGE_DIAG_FACE_H_DOUBLE` の有無だけ。値 (`lineViscCoupling` 2 / 3) とマスク (既定 7) は腕ごとに固定し、他の設定は親 plan の E1b (`run_0311`) と同じにする。
- 値 3 (スカラー対角を重ねた版) を主の対にする。既知の破綻が 122〜150 step で出るので判定が速く、親 plan の作用素の確認もこの値で行った。
  値 2 (スカラー対角なし = 筋のいい形に近い) の対は、主の 4 本を判定した後に回すかを決める (探索扱い、codex 諮問 2026-10-10)。
- 判定は「有限か非有限か」の二値と、非有限になった step。再実行で破綻の step が揺れるので、主の対は各 2 本にする。

## 5. 実装ステップ

1. 台本 `case/45.isobutane_m6_d155/lvcfh.sh` と判定 `lvcfh_judge.py` を書き、run の前に commit する。
2. AWS (`~/forge-wallfit/case/45.isobutane_m6_d155/`) で逐次に回す。
3. 判定し、結果を §6.2 に書く (解釈は上位に諮ってから)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | §4・§6 の諮問と codex plan 段 | 事前登録 (§6) を上位に諮り、codex plan 段の採否を §6.1 に書く。run の前 | F |
| 2 | 台本・判定の実装と主の 4 本の run | `lvcfh.sh`・`lvcfh_judge.py`、case/45 の `run_0540`〜`run_0543`。合格条件: 判定器のゲートが全部通り INVALID の run が無い | O |
| 3 | 結果の解釈と次の一手 | §6 の分岐で判定し、諮問の後に §6.2 へ。次の確認の候補 (諮問 2026-10-10): 支持なら、保存した増幅前の共通の状態・共通の方向で double の経路の実残差の応答と近似作用素を照合 (差分の再現と再評価のノイズを先に検査)。棄却なら、未完了の U-J の列ごとの照合を先に。どちらも長い run より先。値 2 の対を回すかもここで決める | F |

## 6. 検証 (事前登録、2026-10-10、run の前。codex 諮問 [記録](../../notes/reviews/2026-10-10-lvc-faceh-design-diagnose.md) の採否を反映)

- **共通の設定**: 出発 `case/45.isobutane_m6_d155/run_0183_ns_coldmesh_tw300_ext` の res_100000 (sha256 207d39f0…)。`cold_cfl.py prep run_0183_ns_coldmesh_tw300_ext <run> --steps 2000 --out 100 --cfl 4
  --limiter-ref-from run_0183_ns_coldmesh_tw300_ext --line dir --extra res_ro,volume --itj 5 --lvc 3` (親 plan の `run_0311` と同じ。出力の間隔だけ 100 にして、既知の破綻 122〜150 step より前の場を残す)。
  上限なし・緩和 0.7・sweep 5・`implicitSolvePrecision` は設定しない (既定 0)。マスクは既定 (7) で `FORGE_LVC_TERMS` を設定しない。
  序盤 200 step は E1b と同じ 112 節点の帳簿 (`FORGE_DUMP_LEDGER`、両側同じ)。
  バイナリは全 run で `lineM_fp64` (`~/forge-linespeed-fp64`、sha256 05ad8bdf…、FP64、LAYOUT2 既定)。FP64 ビルドでも LHS 全体が倍精度という意味ではない。
- **腕** (逐次、この順): `run_0540_lvcfh_v3_a1` (値 3・既定)、`run_0541_lvcfh_v3_b1` (値 3・`FORGE_DIAG_FACE_H_DOUBLE=1`)、`run_0542_lvcfh_v3_a2`、`run_0543_lvcfh_v3_b2`。
  最大 2000 step、場は 100 step ごと、全残差は毎 step。値 2 の対 (`run_0544` 以降) は主の判定の後に決める。
- **各 run の分類** (`lvcfh_judge.py`):
  - **DIVERGED**: forge の `[detectNaN] Non-finite value detected … at step N` があるか、残差の CSV に非有限がある。破綻の step N を記録する。
  - **FINITE**: 必須の残差の列 (ρ・ρu・ρv・ρE・k・ω・Y0・Y1) がそろい、step 0〜1999 の `outer_begin` の行が一意・連続、全行が有限・非負、`RUN_RC` 0、
    100 step ごとの場の保存量・P・T が有限で ρ・P・T が正。
  - それ以外 (列の欠け・読めない行・負の値・読めない場・途中で止まった) は **INVALID**。INVALID の run が 1 本でもあれば判定しない (判定器は終了コード 1)。
- **分岐** (A = a1・a2、B = b1・b2):
  1. A が 2 本とも DIVERGED、B が 2 本とも FINITE → **この条件・期間の非有限化の回避を支持**。H2 (近似 Jacobian の問題) の否定、H3 (遅らせるだけ) の除外、量子化による増幅の機構の確定には使わない。
     「安定」「収束」とは呼ばない。
  2. A・B とも 2 本とも DIVERGED → **この切替だけでは非有限化を回避できない**。精度依存そのものや元の破綻の機構への寄与は否定しない。熱伝導の K の式・符号の誤りとも決めない。
     B の破綻の step が 2 本とも A の 2 本の最大の 2 倍を超えれば「遅らせる」と記録だけする (同じ機構かは別に確かめる)。
  3. A のどちらかが FINITE → この新しいバイナリで既知の破綻を再現しないので帰属不能 (旧 run を対照に置き換えない)。
  4. A が 2 本とも DIVERGED で B が分かれる → 判別不能。
- **入力のゲート** (1 つでも外れたら判定しない): 台本は起動前に、`run_0183` の最後の res が res_100000.h5 で sha256 が 207d39f0… であることを確かめる。
  判定器は、出発の場の sha256、全 run の `COLD_PAIR.json` の出発の場、格子の実体 (各 run の `nozzle.h5` の MESH を読み直したハッシュが run_0183 と同じ)、
  入力ファイル (境界条件・化学種・壁・`MESH_QUALITY.txt`) が run_0183 と sha256 で一致、`MESH_QUALITY.txt` の VERDICT が PASS、
  `solverConfig.yaml` が §6 の期待どおり (cfl・cfl_pseudo 4、緩和 0.7、キー 5、値 3、方向別 1、ライン 1、blockDPLUR 1、detectNaN 1、timeIntegration 11、nStepInner 5、
  nStepOuter 2000、出力の間隔 100、extraFields [res_ro, volume]、上限と `implicitSolvePrecision` は書かない) かつ全 run で一致、
  起動ログで `lineViscCoupling` = 3、切替の表示が腕と合う、`FORGE_LVC_TERMS` の表示が無い、実効の並びが LAYOUT2、forge の sha256 が lineM_fp64、を確かめる。
- **証拠のゲート** (欠けた run は INVALID): 全 run で初期場 `res_0` が読めて有限・正、`nozzle.h5` がある。DIVERGED の run は、破綻より前に保存されるはずの場 (`res_100` など) が読めて有限・正、
  `res_nan_<N>.h5` が読める、帳簿のヘッダー・112 節点・呼び出し 1〜min(200, N) がそろう。FINITE の run は、100 step ごとの全部の場、帳簿の呼び出し 1〜200、
  区間付きの `check_convergence --segment` の VERDICT (PASS は要求しない) がある。CSV は必須の列・ヘッダーの重複・行の列数を検査し、読めない入力は例外にせず理由付きの INVALID にする。
  破綻の step は、ログの `detectNaN` と CSV の最初の非有限を別に記録し、小さいほうを代表にする。
- **記録だけにするもの**: 各 run の全残差の最大/開始、3 倍を超えた最初の step、最初に非有限になった変数 (`detectNaN` の行)。FINITE の run は末尾 500 step の log10(rms) の傾き × 500 と
  `check_convergence --segment` の VERDICT (区間 0〜1999)。帳簿は選んだ節点だけなので、全域の原因の確定には使わない。
- **やらないこと**: 自動の延長、結果を見てから期間・閾値を変えること、マスク・CFL・緩和・LHS の精度を同時に変えること、この A/B だけで double の評価を既定にすること・診断の機能を消すこと。
- **証拠の保存**: 初期場 (`res_0`)・破綻前の場 (`res_100` など)・`res_nan_*.h5`・格子 (`nozzle.h5`)・帳簿は、result 段のレビューが済むまで消さない。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| 諮問 (§4・§6 の設計) | 2026-10-10 | [2026-10-10-lvc-faceh-design-diagnose.md](../../notes/reviews/2026-10-10-lvc-faceh-design-diagnose.md) | 値 3 を主に各 2 本・2000 step を先に、Major 3・Minor 1 | 全件採用: (M) 分岐の結論を「非有限化の回避」に限定し H1・H2 は排他的でないと §1・§6 に明記、(M) 判定器の FINITE に必須列・全行の有限・非負・読めないファイルの検査を足し INVALID なら終了コード 1 (人工入力で確認)、(M) 出力を 100 step ごとにし、序盤の帳簿を両側に足し、初期場・破綻前の場・res_nan・格子をレビューまで残す、(m) 起動前に出発の場の名前と sha256 を照合し、LAYOUT2 と ISP 0 をゲートに。値 2 の対は主の判定の後に決める |
| plan | 2026-10-10 | [2026-10-10-time_integration-line-viscous-jacobian-faceh-plan.md](../../notes/reviews/2026-10-10-time_integration-line-viscous-jacobian-faceh-plan.md) | GO-with-changes, C0/M3/m2 | 全件採用 (判定器を書き直し、人工入力で main() まで確認: 支持・破綻前の場の欠け・共通の CFL の誤り・step 列の欠け・CSV の非有限がログより早い): M1 分類と独立の証拠のゲート (格子の実体のハッシュ・初期場・破綻前の場・res_nan・帳簿の節点と呼び出し・FINITE 側の区間付き VERDICT)、M2 §6 の期待の設定と照合 (a1 との一致だけにしない)・メッシュ品質の VERDICT、M3 CSV の必須ヘッダー・重複・列数を検査し例外を理由付き INVALID に、異常でも JSON を書き終了コード 1、m4 破綻の step をログと CSV で別に記録し小さいほうを代表に、m5 §1 の 7/5 の断定を探索的な観測に弱めた |

## 7. 影響範囲

- コードは変えない。case/45 に台本・判定と run を足すだけ。

## 変更ログ

- 2026-10-10: 起票 (draft)。親 plan (accepted) §6.17 の再開の条件 2 の延長として、ユーザの合意で立てた。
- 2026-10-10: 設計の諮問と codex plan 段を全件採用して in_progress (§1・§4・§6 を改訂、判定器を書き直し)。主の 4 本 (`run_0540`〜`run_0543`) を投入する。
