forge (自作の圧縮性 FVM ソルバ。CUDA/float32、cell 中心と node 中心 median-dual の 2 離散化、現在は node 主体。
SLAU/Roe/KEEP、block-DPLUR 陰解法、SST、多成分 TP、凝縮、軸対称、ノズル設計ツール design/forge_design を含む) の
リポジトリに対する**外部レビュー**を依頼する。忖度なしで、主張はコードと実測 (run の数値) で検証すること。
結論が「この計画/結果は誤り」でも構わない。両論併記で逃げず、推奨は 1 つに絞ること。

ルール:
- **ファイルを変更しない** (read-only サンドボックスで動いている。読む・実行して確認するのは可)。
- 出力は日本語。識別子・ファイル名は原語のまま。
- 指摘は **Critical / Major / Minor** の重大度付きで、必ず根拠 (`ファイル:行` または `run_*` の数値) と対案をセットで書く。
- リポジトリのルールは `AGENTS.md`、現在仕様は `methods/`、運用手順は `procedures/`、設計判断は `plans/`。
  用語や設定の意味は推測せず `procedures/solver-settings.md` / `procedures/recommended-settings.md` を読むこと。
- 収束の判定は `solver_density_cuda/tools/check_convergence.py <run_dir>` (各 run の `CONVERGENCE_VERDICT.txt`)、
  派生量の定常性は `check_quasisteady.py` の VERDICT を根拠にする。`rms_ro` 単独やスナップショット 1 枚で判断しない。

## 依頼: 診断・設計判断の諮問 (stage = diagnose)

あなたは forge の**診断・設計判断係**である。呼び出し側は実装と run を進めている別のモデル (Claude) で、
**もっともらしい真因に飛びつく前に**あなたに諮っている。仕事は手を動かすことではなく、**次の一手を 1 つに絞ること**。

### 前提
- あなたは呼び出し側の会話を見ていない。下のブリーフと、自分で読んだファイルだけが根拠になる。
  足りなければ推測で埋めずに「何が足りないか」を返す。
- ブリーフは「観測事実 / 期待値と出典 / 再現条件 / 実施済みの操作と結果 / 仮説」に分かれて渡される約束である。
  **観測事実と呼び出し側の解釈が混ざっていたら、まず分け直す**。呼び出し側の要約より、run の数値・コード・
  設定ファイルを自分で確かめた内容を優先する。
- forge を起動しない。`python3` による `residual_history.csv` / `res_*.h5` の読み取りは**統計量だけ**を出す
  (全量ダンプ・長いログ全文をコンテキストに流さない。`*.log`・`*.vtu`・`plans/README.md` は読まない)。

### 診断の作法
1. **「除外済み」というラベルを信用せず、潰した証拠を確認する** (run パス・設定差分・判定区間・VERDICT)。
   証拠が足りない・判定期間が短い・変えた設定が実際には効いていない (YAML の階層違い等) なら**候補へ戻す**。
   証拠が十分な候補は出し直さない。
2. **症状と原因を分ける**。`detectNaN` が指す変数は結果であって原因ではない (EOS 床 → 負密度 → 圧力暴走 → ω の実績)。
   後処理のアーチファクト (2 列混在の抽出、`centCoords` の置換、ソルバ `ypls` の退化) を先に疑う。
3. **このリポジトリで繰り返された真因**を照合する: 投入設定の不整合 (IC と BC、亜音速に超音速 BC)、
   押し出し 2 ノード spanwise、float32 桁落ち (双対幾何・r 重み)、stale build、cross-mesh IC の基底不一致、
   絶対値のゼロ割ガード、境界ノードの凍結、YAML キーの階層違いで黙って無視される設定。
4. 仮説は**確度順に最大 3 つ**。第 1 仮説には根拠を `ファイル:行` か run の数値で付ける。示せないものは「未確認」と明記。
5. **判別する A/B を 1 つだけ**提案する。安く短く回せて、結果がどちらに出ても仮説が 1 つ消えるもの。
   「A なら仮説 1、B なら仮説 2」を先に書く (結果を見てから解釈を作らない)。
6. 少数点の一致・短い窓の値・未収束のトランジェント同士の比較を根拠にしない。

### 設計判断 (plan §4・§6、codex 指摘の採否、result 段の解釈) を諮られたとき
- 採否は指摘ごとに「採用 / 却下 / 要再検証」と理由。根拠が示されていない指摘は自分で該当箇所を読んでから判定する。
- 検証計画は「何が出たら方針が誤りと言えるか」が定量的に書かれているかを見る。
- 既定値の変更・opt-in 機能の削除は、plan の処置欄とユーザ決定の履歴を確認してから判断する
  (「opt-in 残置」は削除対象でない)。
- result 段の解釈は、主張ごとに根拠 run・判定ツールの VERDICT・判定区間が揃っているかを確かめる
  (過渡ピークを定常値と、抽出アーチファクトを物理と誤認した実績は「予想どおり」に見える場面で起きた)。

あなたの結論は**仮説**であって確定ではない。呼び出し側はこの A/B を回して確かめ、plan への反映も呼び出し側が行う。

## ブリーフ (`notes/reviews/briefs/2026-10-10-lvc-faceh-result.md`)

# 諮問ブリーフ: 値 3・マスク 7 の破綻と面エンタルピーの精度の A/B の結果の解釈 (2026-10-10)

日付 2026-10-10。諮問先 codex (diagnose)。AGENTS.md の条件 7 (result の解釈を確定する前)。
plan: `plans/active/time_integration-line-viscous-jacobian-faceh.md` (§6 が事前登録、§6.1 に設計の諮問と plan 段の採否)。親: `plans/accepted/time_integration-line-viscous-jacobian.md`。
判定の出力: `case/45.isobutane_m6_d155/_band_ab/cold_pair/lvcfh_judge.json`、判定器 `lvcfh_judge.py`、台本 `lvcfh.sh` (commit 5e2ac3c5、run の前に commit)。

## 観測事実

- 入力のゲート 45 項目はすべて合格。4 本とも証拠のゲート (初期場・res_100・res_nan・帳簿 1〜min(200, N)・格子) を満たして DIVERGED に分類された。
- 判定器の出力: `主 (値 3): 棄却 — この切替だけでは非有限化を回避できない; 切替でも全部 DIVERGED。切替の破綻 step [121, 122]、既定の最大 122`。
- 破綻の step (ログの detectNaN / CSV の最初の非有限): A1 123/122、A2 123/122、B1 122/121、B2 123/122。最初の非有限はすべて `ro`。
- 序盤の残差 (`rms_ro` / `rms_roUy`、outer_begin):

  | step | A1 | A2 | B1 | B2 | 旧 run_0311 (lineG、従来の並び) |
  | --- | --- | --- | --- | --- | --- |
  | 1 | 1.089e-5 / 2.708e-3 | 同 | 同 | 同 | 同 |
  | 5 | 3.044e-5 / 1.135e-2 | 3.044e-5 / 1.135e-2 | 3.045e-5 / 1.135e-2 | 3.045e-5 / 1.135e-2 | 3.044e-5 / 1.135e-2 |
  | 50 | 3.900e-3 / 1.124 | 3.900e-3 / 1.124 | 3.898e-3 / 1.123 | 3.899e-3 / 1.123 | 3.900e-3 / 1.124 |
  | 100 | 1.297e-2 / 3.380 | 1.296e-2 / 3.379 | 1.317e-2 / 3.431 | 1.318e-2 / 3.427 | 1.296e-2 / 3.380 |
  | 120 | 1.156e-2 / 3.036 | 1.147e-2 / 3.004 | 1.164e-2 / 3.067 | 1.151e-2 / 2.945 | 1.159e-2 / 3.039 |

- 全残差の最大/開始と 3 倍を超えた最初の step は 4 本でほぼ同じ (ρ 約 1670 倍・step 3、ρv 約 3240〜3290 倍・step 2、ρE 約 221 倍・step 22、k 約 19 倍・step 48、ω 約 2100 倍・step 1)。

## 期待値と出典

- 事前登録の分岐 2 (§6): A・B とも 2 本とも DIVERGED → 「この切替だけでは非有限化を回避できない」。精度依存そのもの・元の機構への寄与は否定しない。熱伝導の K の式・符号の誤りとも決めない。
- 設計の諮問 (`notes/reviews/2026-10-10-lvc-faceh-design-diagnose.md`) は、棄却なら次に未完了の U-J の列ごとの照合 (同じ薄層モデルに対する D/K・壁拘束の整合) を先に行うとした。

## 再現条件

- 出発 `run_0183` の res_100000 (sha256 207d39f0…)、値 3・マスク 7・キー 5・方向別・上限なし・cfl 4・緩和 0.7・sweep 5、lineM_fp64 (05ad8bdf…、LAYOUT2)。
  A = 既定、B = `FORGE_DIAG_FACE_H_DOUBLE=1`。run: `case/45.isobutane_m6_d155/run_0540_lvcfh_v3_a1`〜`run_0543_lvcfh_v3_b2`。

## 仮説 (呼び出し側の読み。確定していない)

- H1 (面エンタルピーの float の量子化が値 3・マスク 7 の破綻を起こす) は、この期間・条件では支持されない。A と B の経過は step 50 まで 4 桁一致し、破綻の step も同じ。
  切替は step 100 付近で 1〜2 % の差を生むだけで、増幅の速さ (e 倍の step 数) を変えていない。
- 親 §6.16 の「エネルギーの行の方向微分が double で 1/28」は、ε = 1e-6 の微小な摂動での局所の応答で、反復の増幅を支配する有限振幅の応答とは別だった、と読める (推定)。
- 破綻は ρv が step 2、ω が step 1 で 3 倍を超える早い増幅で始まり、面エンタルピーの経路と無関係な成分が先に育つ。
- 旧バイナリ (lineG、従来の並び) の run_0311 とも経過が一致するので、LAYOUT2 やバイナリの違いは効いていない。

## 問い

1. 登録の判定 (棄却) をこの範囲で書いてよいか。H1 の否定を「この期間・条件で切替だけでは回避できない」より強く書ける根拠があるか (例: 経過が 4 桁一致することを根拠にしてよいか)。
2. 親 §6.16 の局所の応答 (1/28) との関係の書き方。
3. 次の一手は、設計の諮問どおり U-J の列ごとの照合でよいか。ユーザの方針は「速度でなく筋のいい手法」。U-J の後に、熱伝導も結合した全部入りを成り立たせるために何を確かめるべきか (長い run は要らない)。
4. 証拠 (場・帳簿・res_nan) は result 段のレビューまで残す登録だが、この結果で残す範囲を変えるべきか。

## 関連 plan 全文 (`plans/active/time_integration-line-viscous-jacobian-faceh.md`)

```markdown
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
```

## 参考: `case/45.isobutane_m6_d155/_band_ab/cold_pair/lvcfh_judge.json`

```
{
 "plan": "time_integration-line-viscous-jacobian-faceh §6",
 "gates": [
  {
   "check": "出発の場の sha256 が事前に固定した値",
   "ok": true,
   "value": "207d39f0e7f4aa03"
  },
  {
   "check": "run_0540_lvcfh_v3_a1: forge の sha256 が lineM_fp64",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0540_lvcfh_v3_a1: 出発の場が run_0183_ns_coldmesh_tw300_ext/res_100000.h5",
   "ok": true,
   "value": [
    "run_0183_ns_coldmesh_tw300_ext",
    "res_100000.h5",
    "207d39f0e7f4aa03"
   ]
  },
  {
   "check": "run_0540_lvcfh_v3_a1: 格子の実体が出発 run と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0540_lvcfh_v3_a1: 入力ファイルが出発 run と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0540_lvcfh_v3_a1: メッシュ品質の VERDICT が PASS",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0540_lvcfh_v3_a1: 設定が §6 の期待どおり",
   "ok": true,
   "value": {}
  },
  {
   "check": "run_0540_lvcfh_v3_a1: 設定が a1 と一致",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0540_lvcfh_v3_a1: 起動ログの lineViscCoupling = 3",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0540_lvcfh_v3_a1: 切替の表示 = 0",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0540_lvcfh_v3_a1: マスクの表示が無い (= 既定 7)",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0540_lvcfh_v3_a1: 実効の並びが LAYOUT2",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0541_lvcfh_v3_b1: forge の sha256 が lineM_fp64",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0541_lvcfh_v3_b1: 出発の場が run_0183_ns_coldmesh_tw300_ext/res_100000.h5",
   "ok": true,
   "value": [
    "run_0183_ns_coldmesh_tw300_ext",
    "res_100000.h5",
    "207d39f0e7f4aa03"
   ]
  },
  {
   "check": "run_0541_lvcfh_v3_b1: 格子の実体が出発 run と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0541_lvcfh_v3_b1: 入力ファイルが出発 run と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0541_lvcfh_v3_b1: メッシュ品質の VERDICT が PASS",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0541_lvcfh_v3_b1: 設定が §6 の期待どおり",
   "ok": true,
   "value": {}
  },
  {
   "check": "run_0541_lvcfh_v3_b1: 設定が a1 と一致",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0541_lvcfh_v3_b1: 起動ログの lineViscCoupling = 3",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0541_lvcfh_v3_b1: 切替の表示 = 1",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0541_lvcfh_v3_b1: マスクの表示が無い (= 既定 7)",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0541_lvcfh_v3_b1: 実効の並びが LAYOUT2",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0542_lvcfh_v3_a2: forge の sha256 が lineM_fp64",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0542_lvcfh_v3_a2: 出発の場が run_0183_ns_coldmesh_tw300_ext/res_100000.h5",
   "ok": true,
   "value": [
    "run_0183_ns_coldmesh_tw300_ext",
    "res_100000.h5",
    "207d39f0e7f4aa03"
   ]
  },
  {
   "check": "run_0542_lvcfh_v3_a2: 格子の実体が出発 run と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0542_lvcfh_v3_a2: 入力ファイルが出発 run と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0542_lvcfh_v3_a2: メッシュ品質の VERDICT が PASS",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0542_lvcfh_v3_a2: 設定が §6 の期待どおり",
   "ok": true,
   "value": {}
  },
  {
   "check": "run_0542_lvcfh_v3_a2: 設定が a1 と一致",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0542_lvcfh_v3_a2: 起動ログの lineViscCoupling = 3",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0542_lvcfh_v3_a2: 切替の表示 = 0",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0542_lvcfh_v3_a2: マスクの表示が無い (= 既定 7)",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0542_lvcfh_v3_a2: 実効の並びが LAYOUT2",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0543_lvcfh_v3_b2: forge の sha256 が lineM_fp64",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0543_lvcfh_v3_b2: 出発の場が run_0183_ns_coldmesh_tw300_ext/res_100000.h5",
   "ok": true,
   "value": [
    "run_0183_ns_coldmesh_tw300_ext",
    "res_100000.h5",
    "207d39f0e7f4aa03"
   ]
  },
  {
   "check": "run_0543_lvcfh_v3_b2: 格子の実体が出発 run と同じ",
   "ok": true,
   "value": "3fe163acbbbc26df"
  },
  {
   "check": "run_0543_lvcfh_v3_b2: 入力ファイルが出発 run と同じ",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0543_lvcfh_v3_b2: メッシュ品質の VERDICT が PASS",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0543_lvcfh_v3_b2: 設定が §6 の期待どおり",
   "ok": true,
   "value": {}
  },
  {
   "check": "run_0543_lvcfh_v3_b2: 設定が a1 と一致",
   "ok": true,
   "value": []
  },
  {
   "check": "run_0543_lvcfh_v3_b2: 起動ログの lineViscCoupling = 3",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0543_lvcfh_v3_b2: 切替の表示 = 1",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0543_lvcfh_v3_b2: マスクの表示が無い (= 既定 7)",
   "ok": true,
   "value": null
  },
  {
   "check": "run_0543_lvcfh_v3_b2: 実効の並びが LAYOUT2",
   "ok": true,
   "value": null
  }
 ],
 "gates_ok": true,
 "runs": {
  "a1": {
   "class": "DIVERGED",
   "run": "run_0540_lvcfh_v3_a1",
   "run_rc": "1",
   "detectNaN": [
    "ro",
    123
   ],
   "last_step": 122,
   "csv_first_nonfinite": 122,
   "diverged_step": 122,
   "steps_log_csv": [
    123,
    122
   ],
   "evidence_problems": [],
   "records": {
    "rms_ro": {
     "max_over_start": 1667.0653504466914,
     "first_step_over_3x": 3
    },
    "rms_roUx": {
     "max_over_start": 290.45536912058475,
     "first_step_over_3x": 29
    },
    "rms_roUy": {
     "max_over_start": 3239.7684697243167,
     "first_step_over_3x": 2
    },
    "rms_roe": {
     "max_over_start": 220.53172452283633,
     "first_step_over_3x": 22
    },
    "rms_roK": {
     "max_over_start": 19.124863582548905,
     "first_step_over_3x": 48
    },
    "rms_roOmega": {
     "max_over_start": 2093.252516677221,
     "first_step_over_3x": 1
    },
    "rms_roY0": {
     "max_over_start": 1667.0711429739308,
     "first_step_over_3x": 3
    },
    "rms_roY1": {
     "max_over_start": 1667.0629946000245,
     "first_step_over_3x": 3
    }
   }
  },
  "b1": {
   "class": "DIVERGED",
   "run": "run_0541_lvcfh_v3_b1",
   "run_rc": "1",
   "detectNaN": [
    "ro",
    122
   ],
   "last_step": 121,
   "csv_first_nonfinite": 121,
   "diverged_step": 121,
   "steps_log_csv": [
    122,
    121
   ],
   "evidence_problems": [],
   "records": {
    "rms_ro": {
     "max_over_start": 1669.8168248205088,
     "first_step_over_3x": 3
    },
    "rms_roUx": {
     "max_over_start": 290.0195358911928,
     "first_step_over_3x": 29
    },
    "rms_roUy": {
     "max_over_start": 3288.195560766119,
     "first_step_over_3x": 2
    },
    "rms_roe": {
     "max_over_start": 221.2331877392665,
     "first_step_over_3x": 22
    },
    "rms_roK": {
     "max_over_start": 18.81673121516417,
     "first_step_over_3x": 48
    },
    "rms_roOmega": {
     "max_over_start": 2101.8341130704407,
     "first_step_over_3x": 1
    },
    "rms_roY0": {
     "max_over_start": 1669.8273996372782,
     "first_step_over_3x": 3
    },
    "rms_roY1": {
     "max_over_start": 1669.818049857188,
     "first_step_over_3x": 3
    }
   }
  },
  "a2": {
   "class": "DIVERGED",
   "run": "run_0542_lvcfh_v3_a2",
   "run_rc": "1",
   "detectNaN": [
    "ro",
    123
   ],
   "last_step": 122,
   "csv_first_nonfinite": 122,
   "diverged_step": 122,
   "steps_log_csv": [
    123,
    122
   ],
   "evidence_problems": [],
   "records": {
    "rms_ro": {
     "max_over_start": 1668.7103693541837,
     "first_step_over_3x": 3
    },
    "rms_roUx": {
     "max_over_start": 290.40037747104185,
     "first_step_over_3x": 29
    },
    "rms_roUy": {
     "max_over_start": 3238.729949580921,
     "first_step_over_3x": 2
    },
    "rms_roe": {
     "max_over_start": 220.38501533170393,
     "first_step_over_3x": 22
    },
    "rms_roK": {
     "max_over_start": 19.203158609218796,
     "first_step_over_3x": 48
    },
    "rms_roOmega": {
     "max_over_start": 2095.5968133858814,
     "first_step_over_3x": 1
    },
    "rms_roY0": {
     "max_over_start": 1668.7162754217245,
     "first_step_over_3x": 3
    },
    "rms_roY1": {
     "max_over_start": 1668.70818862442,
     "first_step_over_3x": 3
    }
   }
  },
  "b2": {
   "class": "DIVERGED",
   "run": "run_0543_lvcfh_v3_b2",
   "run_rc": "1",
   "detectNaN": [
    "ro",
    123
   ],
   "last_step": 122,
   "csv_first_nonfinite": 122,
   "diverged_step": 122,
   "steps_log_csv": [
    123,
    122
   ],
   "evidence_problems": [],
   "records": {
    "rms_ro": {
     "max_over_start": 1671.4053214248606,
     "first_step_over_3x": 3
    },
    "rms_roUx": {
     "max_over_start": 290.64146605804444,
     "first_step_over_3x": 29
    },
    "rms_roUy": {
     "max_over_start": 3284.60766007267,
     "first_step_over_3x": 2
    },
    "rms_roe": {
     "max_over_start": 221.60403553174555,
     "first_step_over_3x": 22
    },
    "rms_roK": {
     "max_over_start": 18.83263613723407,
     "first_step_over_3x": 48
    },
    "rms_roOmega": {
     "max_over_start": 2111.8970497838063,
     "first_step_over_3x": 1
    },
    "rms_roY0": {
     "max_over_start": 1671.4161308322439,
     "first_step_over_3x": 3
    },
    "rms_roY1": {
     "max_over_start": 1671.4065672062472,
     "first_step_over_3x": 3
    }
   }
  }
 },
 "main_v3": {
  "verdict": "棄却",
  "detail": "切替でも全部 DIVERGED。切替の破綻 step [121, 122]、既定の最大 122"
 },
 "VERDICT": "主 (値 3): 棄却 — この切替だけでは非有限化を回避できない (精度依存そのもの・元の機構への寄与は否定しない); 切替でも全部 DIVERGED。切替の破綻 step [121, 122]、既定の最大 122"
}
```

## 参考: `notes/reviews/2026-10-10-lvc-faceh-design-diagnose.md`

```
# codex 諮問 (diagnose): lvc-faceh-design

- **brief**: [`notes/reviews/briefs/2026-10-10-lvc-faceh-design.md`](../../notes/reviews/briefs/2026-10-10-lvc-faceh-design.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **date**: 2026-10-10
- **commit**: `0be79325` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 3.4 min, rc=0
- **結論**: **判定器・入力ゲート・証拠保存を修正し、値3・マスク7で面エンタルピー精度だけを変える各2本・最大2000 stepのA/Bを先に行う。**
- **extra**: `plans/accepted/time_integration-line-viscous-jacobian.md`, `notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 3 / Minor 1）

| 論点 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| このA/BでH1とH2を判別する | **要再検証・Major** | [plan:72](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:72)の分岐で測れるのは、切替による**2000 step以内の非有限化回避**まで。Bが有限でも、近似Jacobianの問題と精度依存は併存し得る。両側が発散しても、[plan:22](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:22)の「熱伝導K側に原因が残る」には絞れない。支持は「この期間の回避を支持」、棄却は「この切替だけでは回避できない」と記録する。 |
| `FINITE`判定の実装 | **要再検証・Major** | [判定器:61](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcfh_judge.py:61)は存在する残差列だけを検査する。メモリ内の模擬入力で、**残差列なし／負の`rms_ro`だけ**のCSVが、ともに`FINITE`になった。[判定器:87](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcfh_judge.py:87)。必須の平均流・乱流・化学種残差列、全行の有限性・非負性、欠損・破損ファイルを検査し、異常は理由付き`INVALID`にする。さらに[191行](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcfh_judge.py:191)はrun分類が`INVALID`でも入力ゲートが通れば終了コード0になるため、これも直す。 |
| 200 stepごとの場と判定後の削除 | **却下・Major** | 既往の破綻は29、122〜150 stepなので、[現在の出力間隔](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:66)では初期場と非有限化後の場しか残らない可能性が高い。[後片付け規定](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:85)はそれらも消す。既存の[序盤ledger採取](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/e1b.sh:24)を両側同条件で追加し、初期場・床到達前の状態・`res_nan_*`・対応する格子をresultレビューまで保持する。選択節点のledgerだけで全域の原因を確定しない。 |
| 初期場・実効設定のゲート | **採用、補強・Minor** | [準備コード:111](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/cold_cfl.py:111)は親runの最新出力を選ぶ。現判定器には[指定の`res_100000`との事後照合](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcfh_judge.py:155)があるが、実行前にもファイル名・事前固定したハッシュを照合する。またLAYOUT2には[自動フォールバック](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2965)があるため、実効の並びと`implicitSolvePrecision=0`を確認する。FP64ビルドはLHS全体の倍精度を意味しない。 |
| 値3を主、各2本、最大2000 step | **採用** | 値3は[スカラー対角を保持](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:969)し、既存の方向微分診断ともつながる。2000 stepは早期破綻の回避を調べる期間として妥当。ただし各2本は探索的な再現確認で、破綻確率や長期安定性は保証しない。副の値2・1/2本は探索扱いに限り、今回は主の4本を先に判定する。 |

結論: **判定器・入力ゲート・証拠保存を修正し、値3・マスク7で面エンタルピー精度だけを変える各2本・最大2000 stepのA/Bを先に行う。**

第1仮説: 面エンタルピー評価の混合精度が、値3・マスク7の早期増幅に寄与している。確度: **中**。非有限化を回避できるかは未確認。

- 根拠: [切替実装:401](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:401)はTP多成分の面エンタルピー評価を変更する。親planの[§6.16記録](/home/sano/work/forge-faceh/plans/accepted/time_integration-line-viscous-jacobian.md:398)では、S0・p7・ライン2183の自由エネルギー行で‖Jₜp‖が148.9→5.28、近似作用素との相対差が0.975→0.389。ただし、これは局所応答の記録であり、破綻回避の測定ではない。
- 反証条件: ゲートを満たしたBが2本とも2000 step以内に非有限化すれば、**「この切替だけで非有限化を回避できる」という強い予測**を棄却する。精度依存そのものや、元の破綻機構への寄与までは否定しない。

第2仮説: 薄層近似・壁拘束・ライン外結合・SST更新などとの不整合が、精度変更後も増幅を残す。確度: **中**。熱伝導K単独への帰属は未確認。実装自身も[厳密Jacobianではなく前処理近似](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:94)としている。

第3仮説: 精度変更は非有限化を2000 stepより後へ遅らせるだけである。確度: **低・未確認**。Bが2000 step有限でも除外できない。

判別A/B:

- 対象は`case/45.isobutane_m6_d155/run_0540_lvcfh_v3_a1/`〜`run_0543_lvcfh_v3_b2/`の4本。Aは切替なし、Bは`FORGE_DIAG_FACE_H_DOUBLE=1`。他の実効設定・入力・バイナリを固定する。
- 全残差を毎step検査し、最初の増幅・非有限化のstep、保存場の有限性・正値性、序盤の局所状態を記録する。`FINITE`側には区間0〜1999の`check_convergence --segment`のVERDICTを添える。
- **Aが2本とも非有限、Bが2本とも有限** → この条件・期間の非有限化回避を支持。H2の否定、H3の除外、量子化による増幅機構の確定には使わない。
- **両側が2本とも非有限** → 切替単独による回避を棄却。破綻stepの延長は観測として記録し、同じ機構かは別途確認する。
- Aに有限が混じれば帰属不能、Bが分かれれば判別不能、入力・証拠に不備があれば`INVALID`。欠けた反復を旧runで補わない。

支持された場合の次の確認は、**保存した増幅前の共通状態・共通方向で、double経路の実残差応答と近似作用素を照合すること**。差分幅の再現性と再評価ノイズを先に検査し、過去に無効だった有限振幅評価も各評価経路内で閉じる。

棄却された場合の次の確認は、**未完了のU-J列ごとの照合**。[既存の合格は行列全体の最大値による正規化](/home/sano/work/forge-faceh/methods/time_integration/implementation.md:543)に限られる。まず同じ薄層モデルに対するD/K・壁拘束の整合を確認し、それが通ってから実残差との差を調べる。どちらも長期runより先に行う。

やらない方がよいこと:

- Bの2000 step完走を「安定」「収束」「H1の機構が証明された」とする。
- 両側の発散から、熱伝導Kの式・符号が誤りだと決める。
- マスク・CFL・緩和・LHS精度まで同時に変える。
- 判定後すぐに破綻の証拠を削除する、結果を見て期間や閾値を変更する。
- この結果だけでdouble評価を既定化する、診断機能を削除する。

呼び出し側の前提への異議: **H1とH2は排他的ではない。** また、マスク5の比較は親planで探索的観測に格下げされており、「熱伝導Kが必要」は確定していない。point仕上げの保存判定`NOT_SUPPORT`は確認できたが、今回のライン破綻について面エンタルピー精度を除外する証拠にはならない。「LHSは変えない」は固定状態での式についてであり、反復で状態が変わった後のLHS数値まで同一という意味ではない。

不足情報: 対象の旧run・初期HDF5・方向微分の元配列・AWSバイナリはローカルにない。旧run数値はplanの記録として扱い、独立再計算はしていない。コード、保存されたpoint判定集計、判定器の模擬入力を確認した。`forge`は起動せず、ファイルも変更していない。**plan未反映。反映先は対象planの§4・§6・§5.1。**
```

## 参考: `plans/accepted/time_integration-line-viscous-jacobian.md`

```
# ライン陰解法の K と対角に薄層の粘性・熱伝導の Jacobian を入れる (`lineViscCoupling: 2`)

## メタ

- **area**: `time_integration`
- **status**: `done` (2026-10-10。結論は「本線不採用・診断として残置」、§6.17。2026-10-09 に §6.7 で保留し、同日ユーザ指示「つづけて」で原因の切り分けを再開した: §6.8〜§6.16)
- **related_docs**:
  - `methods/time_integration/implementation.md` の「line-implicit」「v2」(`lineViscCoupling`) と「エネルギー行の熱伝導 Jacobian」
- **related_plans**:
  - [`time_integration-implicit-thermal-jacobian.md`](../active/time_integration-implicit-thermal-jacobian.md) (キー 5 と上限 R、§6.0 の run_0223/0224)
  - [`time_integration-line-implicit-speed.md`](../active/time_integration-line-implicit-speed.md) (同じカーネルの速度の改善)
  - [`time_integration-line-implicit-viscous-v2.md`](time_integration-line-implicit-viscous-v2.md) (`lineViscCoupling: 1` = スカラーの結合)
  - [`tooling-nozzle-isothermal-wall-chain.md`](../active/tooling-nozzle-isothermal-wall-chain.md) §5.1 #27 (冷却壁の遅い過渡)
- **created**: `2026-10-09`
- **owner**: Claude (ユーザ指示 2026-10-09「K_ij に粘性の寄与は入れたほうがいいんじゃない」「プラン書いて諮っていいよ」「まずは粘性込みのヤコビアンを試してほしい」)

## 1. 目的

ライン陰解法の近傍行列 K は対流流束の固有値分割の負の側 (−A⁻) だけで、固有値は音波の V ± c とせん断・エントロピーの V (V = 面の法線速度) である。
壁法線のラインでは壁際で V ≈ 0 なので、K が結合するのは音波 (圧力と法線速度) だけで、**せん断 (接線の速度) とエントロピー (等圧での T・ρ) のモードはライン内で結合していない**。
これらを壁法線につなぐのは物理的には粘性と熱伝導だけで、LHS では対角のスカラー 2ν_f δ_f/dcc_f しか入っていない。
冷却ノズル (case/45) の遅い過渡は近壁の T・ρ と境界層の速度なので、薄層の粘性・熱伝導の Jacobian をライン面の K と対角に入れて、
(1) 方向別の擬似 dt の上限 R を外しても安定か、(2) 遅い過渡の step 数がさらに縮むかを確かめる。

## 2. スコープ

- **やる**: ライン面 (Thomas が厳密に解く面) の粘性・熱伝導の薄層 Jacobian (`lineViscCoupling: 2`)、純伝導の試験 (U2)、case/45 の A/B。
- **やらない**: ライン外の面 (従来のスカラーのまま)、SST の k・ω のライン化 (別議題)、Thomas の速度 ([time_integration-line-implicit-speed](../active/time_integration-line-implicit-speed.md))。

## 3. 関連 docs と前提

- 残差 (`viscousFlux_d.cu` の内部面、node・既定 `heatCorrSU2 0`): 面 f (ic0 → ic1) で
  τ = μ_f δ/dcc (u₁ − u₀) + (面平均の勾配による転置・発散・非直交の項)、q = k_f δ/dcc (T₁ − T₀) + (勾配の項)、
  `res_ρu[ic0] += τ`、`res_ρE[ic0] += τ·u_f + q` (u_f = f u₀ + (1 − f) u₁)、ic1 は符号を反転。μ_f = f μ₀ + (1−f) μ₁ (層流 + μ_t)、k_f = `tc_face` (層流 k + c_p μ_t/Pr_t の面の値)。
- LHS の形 (`implicit_defect_correction_block_d` と `lineThomas*`): 節点 i の行は D_i ΔQ_i − Σ_{line j} K_ij ΔQ_j = rhs_i、D = V/Δτ + Σ A⁺ + (粘性の対角)、K = −A⁻ (+ `lineViscCoupling 1` の α I)。
  拡散の残差 R_i = α (Q_j − Q_i) なら D += α、K += α の符号系。
- `lineViscCoupling 1` (スカラー α I を 5 行すべて、対角 2α → α) は case/45 の directional で 20 step で発散した (run_0211)。連続の行にも ρ の拡散を入れ、エネルギーの行を ρE で結合する点が物理の Jacobian と違う (仮説、確かめていない)。

## 4. 設計方針 (案。codex の plan 段と諮問の前)

### 4.1 式 (ライン面 f、自節点 i・ライン上の隣 j、`isLineFace` のとき)

β = μ_f δ/dcc、κ = k_f δ/dcc (残差と同じ面の値)、n̂ = 面の単位法線、P = I + ⅓ n̂ n̂ᵀ (薄層の応力の法線成分 4/3)、f_i = 自節点の補間の重み (ic0 なら f、ic1 なら 1 − f)、
Δu = u_j − u_i、ū = f_i u_i + (1 − f_i) u_j。保存量 Q = (ρ, ρu, ρv, ρw, ρE) について

  ∂u/∂Q = (1/ρ) [−u, I₃, 0]、 ∂T/∂Q = (γ/c_p)(1/ρ) [−(e − ½|u|²), −uᵀ, 1] (e = ρE/ρ − ½|u|²、キー 1 と同じ)。

- 連続の行: 0 (粘性の流束がない)。
- 運動量の行: D_i += β P ∂u_i/∂Q_i、 K_ij += β P ∂u_j/∂Q_j。
- エネルギーの行: D_i += κ ∂T_i/∂Q_i + β (P ū − f_i P Δu)ᵀ ∂u_i/∂Q_i、 K_ij += κ ∂T_j/∂Q_j + β (P ū + (1 − f_i) P Δu)ᵀ ∂u_j/∂Q_j
  (粘性の仕事 τ·ū、τ = β P Δu の微分)。
- 隣 j が壁の節点 (速度の Dirichlet) なら運動量の項と仕事の項の K を 0、等温壁の節点なら熱伝導の項の K を 0 (境界条件で値が固定され Δu_j = 0・ΔT_j = 0)。
  自節点が壁の節点のときは従来どおり `rowDec` の行が単位行になる。
- ライン面では従来のスカラー 2ν δ/dcc (全行) を入れない (上の式で置き換える)。ライン外の面は従来どおり (キー `implicitThermalJacobian` のビットもライン外の面にだけ効く)。
- 物性 (μ、k、c_p、γ) は凍結する。勾配の項 (転置・発散・非直交) は入れない (薄層近似)。

**codex の反映 (2026-10-09、§6.1)**: P = I + ⅓ n̂ n̂ᵀ は現在の離散残差の「勾配を凍結した厳密な微分」ではない (残差の直接の速度差の項は βI で、
転置・発散は面の勾配から入る)。**薄層近似の前処理行列**として採用する。粘性の仕事は τ = β P Δu の薄層の流束モデルの微分で組む
(残差の τ を使う積の微分 (βPδΔu)·ū + τ_res·δū は採らない — 面の応力を LHS に渡す経路がないため。選んだ流束モデルに対して単体照合する)。
等温壁: **値 2 では強制の等温壁の節点の行 4 を Δ(ρE)_w − e_w Δρ_w = 0 (ΔT_w = 0、キーのビット 2 と同じ行) にする**。この拘束の下で、隣が等温壁のときの熱伝導の K を消す。
キー 5 の壁の行 (Δ(ρE)_w = 0) のままでは、静止壁で ΔT_w = −e_w Δρ_w/(ρ_w c_v) (CPG 300 K・Δρ/ρ 1 % で −3 K) が線形系の中に残る。
速度の Dirichlet の壁では Δ(ρu)_w = 0 かつ u_w = 0 なので、運動量・仕事の K の寄与は消さなくても 0 になる (消しても同じ)。

### 4.2 予想される危険

- 連続の行はライン面の粘性の減衰を失う (物理として正しいが、従来は 2α の人工の減衰があった)。V/Δτ が小さい方向別では ρ の行は対流 (音波) の結合と A⁺ だけになる。
  提案の粘性の行列は各節点の速度・温度を保つ密度の補正 δQ = δρ(1, u, v, w, E) を全行で消す (codex が独立のモデルで確認) ので、キー 1 の「行ごとの減衰の基準の食い違い」とは違う。
  ただし人工の減衰を除いた結果、既存の対流・境界・ライン外の結合の誤差が育つ可能性は残る。
- k・ω は従来どおり分離した点陰解法で同じ Δτ を使う。SST の更新を止めた試行 (run_0222) で増幅が弱まったが、凍結の実効性の確認 (壁・入口のピンを除いた roK・roOmega・μ_t) は未完了で、寄与の分離はできていない。
- 勾配の項を入れないので、非直交が強い面や軸対称の τ_θθ では LHS と残差の差が残る。

### 4.3 実装

- 共通関数 `block_dplur::accumulate_thinlayer_visc_jacobian` (`block_dplur_jacobian_d.cuh`、host/device) が面 1 つ分の D_i と K_ij を組む (単体試験と同じ関数)。
- `implicit_defect_correction_block_d` に節点ごとの層流粘性 `vis_lam` を渡す引数を足し、値 2 のときは `thermCond`・`cp`・`fx` も (キーのビット 1 によらず) 渡す (wrapper の `FORGE_BDPLUR_ARGS`)。
- ライン面 (`isLineFace`) で値 2 なら従来のスカラー 2ν δ/dcc の代わりに共通関数の D を足し、`storeLU` の sweep で K を `Kprev`/`Knext` に足す (`rowDec` の行は 0)。
- 強制の等温壁の節点の行 4 は値 2 のとき [−e_w, 0, 0, 0, 1] (§4.1)。
- 設定の検査 (`solverConfig.cpp`): 値 0〜2。2 は node・`lineImplicit 1`・`blockDPLUR 1`・`timeIntegration 11`・`lowMachPrecond < 2`・`heatCorrSU2 0`・
  `nodeIsothermalEnergyBC 0`・(SST のとき) `wallTreatmentSST 0` を要求し、それ以外は起動時に拒否する。`implicitThermalJacobian` との併用は許す (値 1 との併用拒否は残す)。

## 5. 実装ステップ

1. 共通関数と host の単体試験 (`solver_density_cuda/tools/test_line_visc_jacobian.cpp`)。
2. カーネル・wrapper・設定の検査・`methods/time_integration/implementation.md`。FP64 (AWS) と FP32 のビルド。
3. §6 の検証。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | codex plan 段と諮問の採否 | 判断: 2026-10-09 plan 段 GO-with-changes (C0/M5/m2) と諮問 (diagnose) を全件採用 (§6.1)。等温壁の拘束、入力経路と対応範囲、単体試験、判定の分離、同じ水準までの壁時計 | F |
| 2 | 共通関数と単体試験 (U-J) | 合格は §6 U-J。**2026-10-10 訂正 (codex result M1)**: 実装した試験は差を行列全体の最大値で割る判定で PASS。登録した列ごとの相対誤差 ≤ 1e-6 は未確認 (再開の課題、§6.17) | O |
| 3 | 実装 (§4.3) とビルド | 合格: FP64・FP32 のビルドが通る、値 2 + キー 0 で起動する | O |
| 4 | U0・U2・U3・FP32 (§6) | U2 不合格 (950 step)・U3 は Couette が判定不能 (node で壁を動かせない) で Poiseuille に組み直して不合格 (収束が遅い)、§6.0。U0・FP32 は未実施。**処置 (2026-10-10、§6.17)**: 不合格・判定不能は閾値を変えずに残す。U0 は未実施のまま閉じる (値 0 の一部の条件だけ、速度 plan の B1 の短期の一致が代わりの証拠になる)。FP32 は未実施のまま閉じ、値 2・3 の動作確認は FP64 ビルドに限ると適用範囲を書く | O |
| 5 | V-n1 (上限なし 2000 step の A/B)・V-n2 (上限 50 の同じ水準までの壁時計) | V-n1 は B 不合格 (29 step で発散)、上限 50 でも成長 (§6.0)。判断: 2026-10-09 codex 諮問 — 値 2 は本線に使わない、V-n2 へ進む条件を満たさない | F |
| 5b | 実ライン行列の host 再解 A/B (§6.2) | 完了 (§6.3): CUDA の解は書き出した行列の解と一致 (≤ 4e-12; 行列が残差の Jacobian として正しいことは示していない、§6.14)、密度が主の仮説は棄却、弱いモードは壁の近くの ρE が主、σ = 1 で全成分が約 1/25 | O |
| 5c | 次の一手 (値 2 の扱い) | 判断: 2026-10-09 codex 諮問 — 全行のスカラー対角を戻す診断の A/B (§6.4)。エネルギーの行だけの正則化・熱伝導 K の除去を先にしない、値 2 の系統はまだ閉じない | F |
| 5d | 値 3 (診断) の実装と §6.4 の A/B | 完了 (§6.5): B 不合格 → 十分性を棄却 | O |
| 5e | 次の一手 | 判断: 2026-10-09 codex 諮問 — 値 0 でキー 5/7 だけを変える A/B を一組 (§6.6) 行い、その後は値 2・3 を保留して本線の評価に戻る | F |
| 5f | §6.6 の A/B | 完了 (§6.7): 分岐 2 (壁の拘束の変更だけでは壊れない)。値 2・3 は保留 | O |
| 5g | E1 (§6.8 改訂): マスク 7/5 の A/B と事前の行列の確認 | 完了 (§6.9): 熱伝導の近傍 K を除くと 2000 step 有限。観測は分岐 1 の形だが、マスク 0 の事前確認が float の約 2 ulp で許容外 (許容を double の前提で置いた設定の誤り) なので、**登録上は比較無効・探索的な観測** (2026-10-10 に要約を訂正。旧要約「(分岐 1)」は誤り) | O |
| 5h | 次の一手 (熱伝導の K の何が効くか) | 判断: 2026-10-09 codex 諮問 — §6.9 は探索的な観測に格下げ、密度の列だけの A/B (§6.10) を先に、ゲートを直す | F |
| 5i | §6.10 の A/B | 完了 (§6.11): 分岐 2 = 密度の列を除いても 566 step で非有限。同じ破綻の機構にとって必要かは未確定 (2026-10-10 に要約を訂正、codex result M3) | O |
| 5j | 次の一手 | 判断: 2026-10-09 codex 諮問 — 列の切り分けをやめ作用素の切り分け (方向微分) の一組で区切る、またはマスク 5 の加速の評価へ (§6.12)。ユーザの選択 (2026-10-09): 方向微分で作用素を確かめてから区切る | F |
| 5k | 方向微分による作用素の確認 (§6.13) | 完了 (§6.14): S0・S1 とも分岐 (d) 判別不能。判断: 2026-10-09 codex 諮問 — (d) で閉じ、追加は面エンタルピーの float/double の A/B 一組だけ、その後は本線へ | O |
| 5l | 面エンタルピーの精度の A/B (§6.15) | 完了 (§6.16): 精度依存を支持・H-c の説明としても支持 (S0・p7・ライン 2183・エネルギーの行・ε = 1e-6 に限る)。探索はここで終える | O |
| 5m | 本線へ戻る | 判断: 2026-10-09 codex 諮問 — 面エンタルピーは既定の精度のまま、本線へ戻る。値 2・3 は保留のまま (既定 0)。マスク 5 の長期評価・熱伝導 K の列の削除は続けない。速度は plan time_integration-line-implicit-speed の本線 (方向別 dt + 上限・point 仕上げを含む総壁時計) で評価する。2026-10-09 追記: ユーザの指摘を受け、マスク 5 は本 plan では追わず、速度 plan §5.1 #9 の比較の腕としてだけ総壁時計で評価する | F |
| 5n | 残差評価の精度の監査 (移管) | point 仕上げで面エンタルピーの float/double が残差の停滞に効くかの A/B は [time_integration-implicit-thermal-jacobian](../active/time_integration-implicit-thermal-jacobian.md) §5.1 #5 へ移した (本 plan では回さない) | — |
| 6 | result 段のレビューと採否 | 閉じ方 (codex 2026-10-09): 本線不採用・値 2・3 と診断の切替 (`FORGE_LVC_TERMS`・`FORGE_LINE_DUMP_*`・`FORGE_DIAG_FACE_H_DOUBLE`) は診断として残置の判断文書として accepted へ。U0・FP32 は未実施のまま (完了扱いにしない)、§6.13 の (d)・各試験の不合格・判定不能を残す。判断: 2026-10-10 codex result 段 GO-with-changes (C0/M3/m1) を全件採用 (§6.1) — U-J の合格範囲、U2 の「収束」、E1・密度の列の要約を訂正し、総時間の引用に留保を足して accepted へ移した | F |

## 6. 検証 (事前登録、2026-10-09、codex 反映後)

- **U-J (Jacobian の単体照合、host)**: 共通関数の D = −∂R_i/∂Q_i・K = ∂R_i/∂Q_j を、同じ薄層の流束モデル R (τ = βPΔu、q = κΔT、仕事 τ·ū) の中心差分と照合する。
  状態は CPG・TP 相当 (節点ごとの γ・c_p)、高速の接線流 (|u| 1700 m/s)、異なる密度、f_i = f と 1 − f の両方向、壁の拘束 (速度・温度) を含む 200 組。
  合格: 列ごとの相対誤差 ≤ 1e-6 (double)、零空間 δQ = δρ(1, u, v, w, E) で D・K の作用 ≤ 1e-12 (相対)、等温壁の行の拘束で ΔT_w = 0 (≤ 1e-12)。
  短いライン (8 節点) を組んで block-Thomas の前進消去・後退代入 (host の写し) と密行列の解の差 ≤ 1e-10。実装した CUDA の K の配置 (Kprev/Knext の符号) は U2・U3 で確かめる。
- **U0 (既定と値 1 が不変)**: `lineViscCoupling 0` と 1 で、旧 (4d394a71 の linespeed) と新のバイナリの 1 step・20 step の場の差が、同じバイナリの再実行の最大の 3 倍以内 (速度の plan §6 と同じ形)。
- **U2 (純伝導)**: case/52 の流体の板 (U1 と同じ、下 300 K・上 350 K、静止)、ライン + 方向別。値 0 と 2 を cfl 5・50 で 5000 step・25 ごと
  (`case/52.conjugate_slab/run_0008_u2_lvc0_cfl5`・`run_0009_u2_lvc2_cfl5`・`run_0010_u2_lvc0_cfl50`・`run_0011_u2_lvc2_cfl50`)。
  合格: 値 2 が L∞ < 0.05 K に 250 step 以内 (point キー 0 cfl 20 の 2500 step の 1/10)、NaN・非物理値なし。
- **U3 (Couette、せん断と粘性の仕事)**: 同じ板の上壁を**紙面に垂直な z 方向**に U = 300 m/s で動かす (側面は slip 壁なので x 方向には流せない。2D の格子で u_z(y) だけが立つ)。
  両壁 300 K、空気 CPG、1013 Pa、定数 μ = 1.81e-5・k = 0.0241。解析解 u_z = U y/H、T = T_w + (μU²/2k)(y/H)(1 − y/H) (中央で +8.45 K)。
  値 0 と 2 を cfl 5・50 で 5000 step・25 ごと (`run_0012_u3_lvc0_cfl5`・`run_0013_u3_lvc2_cfl5`・`run_0014_u3_lvc0_cfl50`・`run_0015_u3_lvc2_cfl50`)。
  合格: 値 2 で u_z の L∞ < 0.1 % U・T の L∞ < 1 % (8.45 K の) に、同じ cfl の値 0 より少ない step で入る。NaN なし。
- **U3 の組み直し (2026-10-09、事後。Couette は node で壁を動かせず成立しなかった、§6.0)**: Poiseuille。同じ板の側面を周期、両壁 300 K で止めたまま、
  体積力 (単位体積あたり、`bodyForce`) f = 8μU/H² で押す (中央 U = 300 m/s、f = 434.4 N/m³)。解析解 u = (f/2μ) y (H − y)、T = T_w + (f²/192μk)(H⁴ − (H − 2y)⁴) (中央で +22.5 K、
  体積力の仕事と粘性の仕事が打ち消して k T'' + μ(u')² = 0)。ライン + 方向別、値 0 と 2 × cfl 5・50、5000 step・25 ごと (`case/52.conjugate_slab/run_0020〜0023_u3p_*`、`u3p.sh`)。
  合格: 値 2 で u の L∞ < 0.1 % (中央の速度) かつ T の L∞ < 5 % (22.5 K の。温度の 4 次式の離散化の誤差は 16 セルで約 1.6 % の見積もり) に、同じ cfl の値 0 より少ない step で入る。NaN なし。
  最終の L∞ も記録する (離散化の誤差の水準)。
- **FP32**: U2 と U3 を FP32 のビルドで 1 本ずつ。合格: NaN なし (FP32 のビルドは FP64 の結果の後に作る)。
- **V-n1 (十分性の試験、診断の A/B)**: run_0183 の res_100000 から、同じ新バイナリで `lineViscCoupling` 0 と 2 だけを変える (キー 5・方向別・上限なし・cfl 4、2000 step)。
  全残差を毎 step、場を 200 step ごと。短期合格 (B = 値 2): 非有限・非物理値なし、線形解の失敗 0、全残差の最大が開始の 10 倍以内、
  末尾 500 step (1500〜1999) の log10(rms) の最小二乗の傾き × 500 < 0.1 桁 (全列)。線形解 (Thomas) の失敗は現在の forge では数えて出力していない
  (`line_fail_d` は内部のフラグだけ) ので、この項目は判定に使えないと記録する。run 名: `case/45.isobutane_m6_d155/run_0260_vn1_lvc0` (A)・`run_0261_vn1_lvc2` (B)、
  バイナリは 5ab83056 + typedef double (`~/forge-linevisc-fp64`)、`lineDtDirectionalCap` は書かない (0)、場は 200 step ごと。
  A (値 0) で成長が再現し B で止まれば第 1 仮説を支持。B でも同じモードが育てば「この変更だけで十分」を棄却。A が再現しなければ判別不能。2000 step の合格は短期安定の判定に限る。
- **V-n2 (速さ)**: 上限 50 で値 0 と 2 を同じ新バイナリ・同じ IC (run_0183 の res_100000) で回し、**事前登録の同じ状態の水準 (欠損 |Σ| ≤ 0.1 kg/s と |res_ro| の局所ノルムの減少、
  θ_r の `check_quasisteady` STEADY (末尾 2 万 step、4 点以上、drift 0.0005))** に達するまでの step 数と、専有 GPU・profiler なしの ms/step の積 (壁時計) で比べる。
  固定 step の比較は途中の診断。ms/step の増分の上限 10 %。V-n2 は値 2 の上限なしの安定 (V-n1) の結果を見てから上限なしでも行う。
- **判定の分離**: 短期安定 (上記)、収束 (`check_convergence --segment`、区間を応答に書く)、採用 (状態の水準と切り戻し、thermal-jacobian plan §6.0) を別に書く。
  Q_w は `cold_series.py` (CPG の定数で再構成した共通の後処理の指標) なので、MW の値や熱収支の根拠に使わない。

### 6.0 結果 (2026-10-09)

- **U-J**: PASS (中心差分との差 2.3e-8、零空間 2.5e-16、等温壁の拘束 1.5e-16、8 節点ラインの解 1.6e-15、float と double 1.3e-5)。
  **訂正 (2026-10-10、codex result M1)**: 実装した試験 (`test_line_visc_jacobian.cpp`) は中心差分との差を D・K それぞれの行列全体の最大値で割って判定している。事前登録した「列ごとの相対誤差 ≤ 1e-6」は確かめていない (小さい列の相対誤差はこの判定では保証されない)。PASS は全行列の最大値で正規化した判定に限る。
- **U2** (`case/52.conjugate_slab/run_0008〜0011`): 値 0 cfl 5 は 5000 step で L∞ 0.198 K、値 2 cfl 5 は 0.094 K、値 0 cfl 50 は非有限、**値 2 cfl 50 は 950 step で 0.05 K、5000 step で 1.4e-6 K**。
  値 0 が壊れる cfl 50 で、値 2 は温度の解析解誤差がこの水準まで下がった (熱伝導の項のライン内の結合は効いている、という観測) が、事前登録の「250 step 以内」は**不合格**。残差の収束は確かめていない (`check_convergence` の VERDICT なし。2026-10-10 に「収束」の表現を訂正、codex result M2)。静止流なので運動量・仕事の項は試していない。
- **U3**: 判定不能。z 方向に壁を動かす版 (`run_0012〜0015`) も、側面を周期にして x 方向に動かす版 (`run_0016〜0019`、事後に組み直した) も、壁の節点の速度が 0 のまま。
  node の壁は `nodeWallDirichlet` で no-slip に固定する作りで、bcond の Ux/Uz は node では効かない。せん断・仕事の項は**未検証**。
- **U3 (Poiseuille、`case/52.conjugate_slab/run_0020〜0023_u3p_*`)**: 5000 step の速度の L∞ (中央の速度に対する比) は値 0 cfl 5 で 0.276、値 2 cfl 5 で 0.259、値 0 cfl 50 は非有限、
  値 2 cfl 50 は 0.121 → どれも許容差 (0.1 %) に未到達で**不合格**。値 2 は値 0 が壊れる cfl 50 で安定 (運動量の結合の符号の誤りなら壊れるはず、という間接の証拠)。
  cfl を 10 倍にしても 1 step あたりの減衰は 2.7e-4 → 4.3e-4 (1.6 倍) で、U2 の熱 (5.6e-3) より 13 倍遅い。この板は縦横比 1.6 なので、ライン外の面に残るスカラー 2α_x が
  ライン方向の滑らかなモードの応答を抑えている可能性 (仮説、未確認。ノズルの近壁では α_x/α_y = 1/AR² で小さい)。
- **V-n1** (`case/45.isobutane_m6_d155/run_0260_vn1_lvc0`・`run_0261_vn1_lvc2`): A (値 0) は残差の最大/開始 ρ 572・ρv 942 まで振れて回復 (末尾 500 step の傾き × 500 は −0.13〜−0.34 桁、k +0.04)。
  **B (値 2) は 1 step 目から爆発して 29 step で非有限** (rms_ro 8.0e-6 → 1.8e-4 → 3.9e-3 → 2.6e-2)。ρ の非有限は列 60〜70 の全層、T は列 65 の層 52/53 で 50 K / 4327 K と交互。
  → 事前登録の判定で B は不合格 (「この変更だけで十分」を棄却)。即座の発散は登録時に想定していなかった。
- **切り分けの計測** (原因は未確定、codex 諮問中): `run_0264_diag_lineonly_lvc2` (Δτ は point) は `run_0201` (値 0) とほぼ同じ推移 (step 99 の ρ 0.585 vs 0.586)。
  `run_0265_diag_linedir_cap50_lvc2` (上限 50) は step 99 で ρ 98 倍・k 7000 倍に成長 (値 0 の run_0223 は 0.53 倍)。
  `run_0266/0267_diag_s3` (3 step・毎 step): 1 step 目の ρ の変化の最大は値 0 で 9.3e-3 (列 37・1 層)、値 2 で 0.29 kg/m³ (列 12・5 層)、2 step 目の P の変化は 1.5e3 vs 9.1e4 Pa。
  Δτ が小さいと無害で、Δτ を伸ばすほど速く壊れる。

- **codex 諮問 (2026-10-09、[記録](../../notes/reviews/2026-10-09-line-viscous-jacobian-divergence-diagnose.md)) の採否 (全件採用)**:
  第 1 仮説 (確度 中) = 全行のスカラー対角を外したので、方向別の Δτ で密度・圧力が主の補正を抑えきれず初回から過大な補正を返す。
  「ライン方向に一様な圧力のモードでほぼ特異」という限定した機構は未確認 (実行列には質量項・ライン外の A⁺・境界の半割面・軸対称の項が残る)。
  等温壁の行も値 0 → 2 で同時に変えているので、次の A/B では両側で同じにする。V-n1・U2・U3 の不合格と判定不能は閾値を変えずに記録する。
  値 2 は本線に使わない。連続の行だけの正則化 (案 b・c) を先に実装しない。U3 の遅さをライン外のスカラーだけで説明しない (ライン外には音響の対角 cΔy もあり、静止一様の概算でスカラーの約 56 倍)。
  **次の一手 (§6.2 に事前登録)**: 1 step 目の実際のライン行列を取り出し、従来のスカラー対角を戻さない/戻すだけの host の再解 A/B。

### 6.2 事前登録: 実ライン行列の host 再解 A/B (2026-10-09)

- **採取**: `run_0183` の res_100000 から値 2 (キー 5・方向別・上限なし・cfl 4) を 1 step。デバッグの環境変数で、ライン上の全節点の D (25)・Kprev・Knext (各 25)・
  各 sweep の rhs (5) と Thomas 後の補正 (5)・ライン面のスカラー 2ν_eff δ/dcc の和・拘束の行のフラグ (壁・等温壁・軸)・V/Δτ を、列 12・33・65・1640・2183 を含むラインについて書き出す。
- **host 再解**: 同じ D・K・rhs で (A) σ = 0 (そのまま) と (B) σ = 1 (スカラーの和を拘束でない行の対角へ戻す) を密行列で解く。(A) は CUDA の補正と照合する (相対差 ≤ 1e-8 でなければ実装経路を疑う)。
  補正は保存量から δu = (δ(ρu) − u δρ)/ρ、δT = {δ(ρE) − u·δ(ρu) + (½|u|² − e) δρ}/(ρ c_v) を作って評価する。固定の物理尺度 (ρ・u・T の代表値) で行・列を無次元化し、
  拘束の自由度を消去して最小特異値とその右特異ベクトル (弱いモード)、rhs のそのモードへの射影を出す。
- **分岐**: (A) の大きな補正が密度が主の弱いモードに集中し (無次元の補正の二乗ノルムの 90 % 以上が密度の成分)、(B) でその密度の補正のノルムが (A) の 1/10 以下 → 第 1 仮説を支持
  (この場合も仕事の正しさや長期の安定は判定しない)。その構造が見えない → 狭い準特異の仮説を棄却。host 解が CUDA と合わない → 実装経路。合っていて後の sweep だけで増幅 → 第 2 仮説 (lag の増幅)。
  (B) が全補正を一様に小さくしただけ → 正則化への感度を示しただけとし、圧力のモードの説の支持にはしない。

### 6.3 §6.2 の結果 (2026-10-09)

採取: `case/45.isobutane_m6_d155/run_0288_dump_lvc2` (lineB2 のバイナリ、値 2・キー 5・方向別・上限なし・1 step、`run_0288_dump_lvc2_linedump/`)、
列 12・33・65・1640・2183 の壁の節点を含むライン 5 本 × 121 節点、5 sweep。解析 `linedump_analyze.py` → `_band_ab/cold_pair/linedump_run_0288_analysis.json`。
- **実装の経路**: host の再解 (σ = 0) × implicitRelax と CUDA の dq の相対差は 2e-13〜4e-12 → 書き出した D・K・rhs の系を CUDA の Thomas が正しく解いている
  (2026-10-10 訂正: 当初「実装経路の誤りは無い」と書いたが、示したのは解法の整合だけ。D/K の組立・物性の入力・実残差の Jacobian との整合は示していない。末尾と §6.14)。
- **大きな補正は 1 回目の sweep から** (入口寄りの 3 本で |δρ/ρ| 0.97〜1.3 %、|δu| 0.8〜1.0 m/s、|δT| 7.8〜9.9 K が sweep 0〜4 でほぼ一定) → 後の sweep の lag の増幅ではない。
- **弱いモードは密度が主ではない**: 無次元化した行列の最小の右特異ベクトルの成分は ρE 94 %・ρ 5 %・ρv 1 % (入口寄りの 3 本)、位置は壁から 6〜12 節点目。
  1 回目の補正の成分も ρE 94 %。→ 事前登録の分岐で**「密度が主の準特異」の狭い仮説は棄却**。
- **σ = 1 (スカラーを戻す)**: 入口寄りの 3 本で最小特異値 5.8e-5〜8.3e-5 → 1.4e-3〜1.8e-3 (20〜23 倍)、補正は全成分がそろって約 1/25 (|δT| 8.6 → 0.35 K)。
  → 事前登録の分岐で「正則化への感度を示しただけ」に当たる。下流の 2 本 (1640・2183) は V/Δτ がスカラーの 0.1〜0.15 倍で、補正は 1/10 か不変。
- 観測のまとめ: 値 2 では入口寄りのラインの壁から 6〜12 節点目に、尺度化した保存量で ρE が主 (94 %) の弱いモードができ、1 step 目から大きな補正を返す。
  **訂正 (2026-10-09、codex 諮問)**: 「温度が主のモード」とは言えない。保存量を ρ・ρc・ρc² で尺度化すると、静止・δT = 0 の密度の変化でも δ(ρE) = e δρ のため ρE の割合が 95 % になる
  (代数的な反例)。また SVD は拘束の自由度を消去していない。原始量 (δρ/ρ・δu/c・δT/T) に直し、拘束を消去した系で判定し直す (§6.4)。
  |δT| 7.8〜9.9 K は緩和前の値で、緩和 0.7 を掛けると 5.99 / 6.92 / 5.47 K (列 12・33・65)。列 2183 で「不変」なのは最大の密度・温度の補正で、保存量の補正のノルムは 0.448 倍、最大の速度の補正は 0.137 倍。
  host と CUDA の一致は「採取した系の解法との整合」の確認であって、D/K の組立・物性の入力・実残差との整合の確認ではない。

### 6.4 事前登録: 全行のスカラー対角を戻す A/B (2026-10-09、codex 諮問 [記録](../../notes/reviews/2026-10-09-line-viscous-thermal-weak-mode-diagnose.md))

- **診断用の値 3** (`lineViscCoupling: 3`): 値 2 と同じ薄層の D/K と同じ壁の拘束に、ライン面のスカラー 2ν_eff δ/dcc を全行の対角へ足す (拘束の行は後で単位行に上書きされるので「全自由行」)。
  値 2 と同じ設定の検査。恒久の採用はしない (診断)。
- **A/B**: 同じ新バイナリ・run_0183 の res_100000・キー 5・方向別・上限なし・cfl 4・緩和 0.7・sweep 5 で、A = 値 2 (`run_0300_vn1b_lvc2`)、B = 値 3 (`run_0301_vn1b_lvc3`)、最大 2000 step、
  場は 200 step ごと、全残差は毎 step。上限 50 の組は同時に足さない。B の 1 step 目のライン行列を書き出し (`run_0302_dump_lvc3`)、host の再解 (σ = 0 = そのまま) と CUDA の差を先に確かめる。
- **分岐**: A で既知の早期の増幅が再現し、B が非有限・非物理値なし、全残差の最大/開始 ≤ 10、末尾 500 step の log10(rms) の傾き × 500 < 0.1 桁 (全列) を満たす
  → 「全行の対角の復元がこの条件の短期安定に十分」を支持 (熱伝導が真因・Jacobian が正しいとは判定しない)。B がどれかを満たさない → 十分性を棄却。A が再現しない → 判別不能。
  `check_convergence` の VERDICT は両側で別に記録する。V-n2 へ進むなら目的量の `check_quasisteady` と共通の終了条件までの総壁時計を要求する。
- 解析の改訂 (`linedump_analyze.py`): 弱いモードを原始量に直した割合 (δρ/ρ・δu/c_ref・δT/T_ref、T_ref 300 K) と、拘束の自由度を消去した系 (壁の運動量 0、等温壁 δ(ρE) = e_w δρ、軸の半径運動量 0) の SVD を併記する。

### 6.5 §6.4 の結果 (2026-10-09、バイナリ 6d49738e + typedef double、sha256 8e989dc6…)

- **改訂した解析で値 2 の書き出し (run_0288) を見直し**: 拘束の自由度を消去し原始量 (δρ/ρ・δu/300 m/s・δT/300 K) に直すと、入口寄りの 3 本の弱いモードは δT 81〜83 %・δρ/ρ 14〜16 %・δu 3 %、
  1 回目の補正は δT 83〜86 %。この尺度では温度が主 (尺度の選び方に依存する)。下流の 2 本は混ざり、列 2183 の弱いモードは δu 81 %。最小特異値は消去の前後でほぼ同じ。
- **A/B**: A = 値 2 (`case/45.isobutane_m6_d155/run_0300_vn1b_lvc2`) は 29 step で非有限 (run_0261 と同じ step で再現、step 2 で ρ 487 倍・ρv 1270 倍)。
  B = 値 3 (`run_0301_vn1b_lvc3`) は 122 step で非有限。ρv が先に育つ (step 1: 2.6 倍、5: 11 倍、20: 34 倍、50: 1080 倍)。
  → **事前登録の分岐で B 不合格 = 「全行の対角の復元がこの条件の短期安定に十分」は棄却**。
  B の 1 step 目の書き出し (`run_0302_dump_lvc3`) は host の再解と CUDA の差 ≤ 5e-13 (採取した系の解法とは整合)。
- **観測 (別のバイナリの比較なので参考)**: 値 0 (スカラーだけ、薄層の結合なし) は同じ条件の `run_0260_vn1_lvc0` で 2000 step 生き残る (序盤の残差/開始は 0.5〜2 倍、最大 ρ 572 倍で回復)。
  値 3 = 値 0 のスカラー + 薄層の D/K + 等温壁の拘束の行の変更、なので「薄層の D/K を足すこと」か「等温壁の行の変更」のどちらか (または両方) が不安定を強めている。まだ分けていない。

### 6.6 事前登録: 等温壁の拘束の行だけの A/B (2026-10-09、codex 諮問 [記録](../../notes/reviews/2026-10-09-line-viscous-lvc3-result-diagnose.md))

- 同じバイナリ (lineE、6d49738e + typedef double、sha256 8e989dc6…)・run_0183 の res_100000 から、`lineViscCoupling: 0` に固定し、
  A = `implicitThermalJacobian: 5` (`case/45.isobutane_m6_d155/run_0303_wallA_tj5`)、B = `7` (= 5 + 等温壁の拘束の行、`run_0304_wallB_tj7`)。
  方向別・上限なし・cfl 4・緩和 0.7・sweep 5、最大 2000 step。全残差は毎 step、場は 200 step ごと (`extraFields: [res_ro, volume]`)、
  序盤 200 step は列 12・33・40・50・65・1640・2183 の壁から 0〜15 層目 (112 節点) の状態と残差を `FORGE_DUMP_LEDGER` で毎 step 記録する。実効の設定の差がキーだけであることを COLD_PAIR.json の config_diff で確かめる。
- **分岐** (V-n1 の「最大/開始 ≤ 10」はこの因果の判別には使わない。対照の run_0260 自体が ρ 572 倍まで振れるため):
  1. A が 2000 step 有限、B が途中で非有限か非物理値 → 壁の拘束の変更だけで不安定になることを支持し、「薄層の変更が破綻に必要」を棄却 (同じ原因とはまだ断定しない)。
  2. 両側とも 2000 step 有限・非物理値なし → この条件・期間で「壁の拘束の変更だけによる破綻」を棄却。薄層の変更とその壁の拘束との相互作用は残す (過渡の場の一致や収束は意味しない)。
  3. A も破綻 → 基準が再現しないので帰属不能 (run_0260 を対照として補わない。バイナリ・実効設定・restart を確かめる)。
  両側で全残差の最大/開始、最初に増えた step、末尾 500 step の傾き、`check_convergence` の VERDICT と判定区間を記録する。
- **この A/B の後**: 値 2・3 は本線に使わないまま保留 (コードは消さず既定も変えない)。本線は「方向別 + キー 5 + 上限 50 で加速し point で仕上げる」の、仕上げ込みの共通の終了条件までの壁時計の評価に戻る。
- **訂正 (codex)**: 値 0 + キー 5 は「スカラーだけ」ではなく、熱伝導の自己側の Jacobian (キー 5 の対角) を持つ。値 0 → 3 は主に運動量・仕事の D と近傍の K、壁の拘束の変更を足す比較である。
  §6.5 の原始量の割合 (δT 81〜83 %) は保存量の尺度の縮約行列の最小の右特異ベクトルを原始量で表示したもので、外側の反復の不安定なモードや熱伝導が原因であることの証明には使わない。

### 6.7 §6.6 の結果と保留 (2026-10-09)

- A = キー 5 (`case/45.isobutane_m6_d155/run_0303_wallA_tj5`)、B = キー 7 (`run_0304_wallB_tj7`)、ともに値 0・同じバイナリ (lineE)・2000 step。config_diff はキーの値だけが違う。
  **両側とも 2000 step 有限、200 step ごとの場に非有限・非正の ρ・P・T なし → 分岐 2: この条件・期間で「壁の拘束の変更だけによる破綻」は棄却**。
  値 2・3 の早期の発散には薄層の D/K (またはそれと壁の拘束の相互作用) が要る (帰属の範囲はここまで)。
  残差の最大/開始は A で ρ 550・ρv 873 倍、B で 625・1020 倍、3 倍を超えた最初の step は A で ρv 5・ρ 210、B で ρv 120・ρ 141。
  末尾 500 step の傾き × 500 は A で −2.92 (ρ)、B で −0.26 (ρ)。check_convergence (0〜1999): A は NOT CONVERGED (stalled/plateau、低下 2.7 桁)、B は NOT CONVERGED (still converging、低下 0.8 桁)。
  → 拘束の行の変更は戻りを遅くするが、単独では壊れない。A は旧バイナリの run_0260 の振れ方 (最大 ρ 572 倍) を再現した。
  帳簿 (112 節点 × 200 step、`ledger.csv`) は採取済みで、未解析 (保留中の解析の材料として残す)。
- **保留 (codex 諮問の結論どおり)**: 値 2・3 は本線に使わず、ここで保留する。コードは消さず、既定 (値 0) は変えない。再開するときの候補は「薄層の D/K のどの項 (運動量・仕事・熱伝導の K) が要るか」の単因子の切り分け。
  本線は「方向別 + キー 5 + 上限 50 で加速し point で仕上げる」の、仕上げ込みの共通の終了条件までの壁時計の評価に戻る (thermal-jacobian plan §6.0)。

### 6.8 事前登録: 原因の切り分けの再開 (2026-10-09、ユーザ指示「つづけて」)

目的は値 2・3 を本線に使うことではなく、**ノズルの方向別の Δτ でなぜ壊れるか**を切り分けること。
- **E1 (項の切り分け、ノズル)**: 値 3 に、薄層の項を選ぶ診断のマスク `FORGE_LVC_TERMS` (ビット 1 = 運動量の D/K、2 = 熱伝導の D/K、4 = 粘性の仕事の D/K、既定 7) を足す。
  熱伝導の項を外したときはキー 5 の対角 (従来のスカラー + 熱伝導の自己側) をライン面にも入れて、基準を「値 0 + キー 5」とそろえる (熱伝導の D は キー 5 の対角と同じ式なので、
  ビット 2 だけの版は「値 0 + キー 5 + 熱伝導の近傍 K」になる)。値 3 の等温壁の行 ([−e_w,0,0,0,1]) は全版で同じ。
  同じ新バイナリ・run_0183 の res_100000・キー 5・方向別・上限なし・cfl 4・2000 step で、マスク 1 (`run_0305_e1_m1`)・2 (`run_0306_e1_m2`)・4 (`run_0307_e1_m4`)。
  対照は §6.7 の B (`run_0304_wallB_tj7` = 値 0 + キー 7 = 値 0 + キー 5 + 等温壁の行、2000 step 有限) と §6.5 の値 3 (`run_0301`、マスク 7、122 step で非有限)。
  **分岐**: 単独の項で 2000 step 以内に非有限・非物理値が出たら、その項は「値 0 + キー 7 に足すと破綻を起こすのに十分」。どれも出なければ「単独では起こさず、組み合わせが要る」。
  複数出たら、非有限までの step 数の短い順に記録する (原因の順位ではない)。
- **E2 (対流のある平板、case/48)**: `case/48.flat_plate_cooled_m4` の B (300 K、y₁ 3 µm、M 4.19、低 Re SST) の収束場 (`run_0005_B_tw300_y3/res_48000.h5`、単精度) を、
  FP64 の変換器で変換し直した同じ格子に載せ、同じ新バイナリで cfl 2 (この case の生産の値)・2000 step:
  P = point (ラインなし、対照)、A = ライン + 方向別 + キー 5 + 値 0、B = 同 + 値 2 (いずれも上限なし)。
  **分岐**: P・A が有限で B だけ破綻 → 平板 (対流・乱流・冷却壁) で再現する (以後の切り分けを安い平板でできる)。P・A・B とも有限 → 平板では再現せず、ノズルに固有の条件
  (軸対称・縮流部の低マッハ・縦横比 約 4000 など) が要る。A も破綻 → 平板の方向別そのものが成り立たないので、上限 50 で取り直す (B も同じ条件で)。
- どちらも全残差を毎 step、場は 200 step ごと。`check_convergence` の VERDICT と判定区間を記録する。

**§6.8 の改訂 (2026-10-09、codex 点検 [記録](../../notes/reviews/2026-10-09-line-viscous-resume-diagnose.md) を全件採用)**:
- **E1 を作り直す**: マスク `FORGE_LVC_TERMS` は ビット 1 = 運動量の D/K、ビット 2 = **熱伝導の近傍 K だけ**、ビット 4 = 仕事の D/K。熱伝導の D は共通の経路で常に 1 回だけ入れる
  (キー 5 の処理を足すとスカラーが二重になる、という指摘)。最初の一組は **A = マスク 7 (`run_0305_e1_m7`)、B = マスク 5 (`run_0306_e1_m5`、熱伝導の K だけを抜く)**、
  値 3・キー 5・方向別・上限なし・cfl 4・2000 step、同じ新バイナリ、逐次。序盤 200 step は §6.6 と同じ 112 節点の帳簿。
  **事前の確認**: 1 step 目の行列を書き出し (マスク 7・5・0 と、値 0 + キー 7)、D・rhs・拘束の行が 7 と 5 で一致し、K の差が熱伝導の成分 (行 4) だけであること、
  マスク 0 の D が値 0 + キー 7 の D と相対 1e-12 以内で一致すること (丸めの順序は違う)。どれかが外れればこの比較は無効。
  **分岐**: A だけ破綻し B が 2000 step 有限・非物理値なし → この条件・期間で熱伝導の K を除くことが破綻の回避に十分 (熱伝導単独の誤りや長期の安定は証明しない)。
  A・B とも破綻 → 「熱伝導の K が破綻に必要」を棄却し運動量・仕事の側を候補に残す。A が破綻を再現しない → 帰属不能 (run_0301 を代わりの対照にしない)。
  全域の最初の非有限・非正の ρ/P/T、EOS 床、δρ/ρ・δu・δT の増える場所を記録し、両側の `check_convergence` の VERDICT と判定区間を保存する。実効の `implicitSolvePrecision` も記録する (FP64 ビルドでも既定は float)。
- **E2 は後回し**: 準備のスクリプトは restart の失敗を確かめずに起動できる穴があるので、restart の終了コード・保存量の移送の検査・メッシュの品質 (座標の一致) を起動の条件にする。
  値 2 の側だけ壁の拘束が変わるので、行うなら両側キー 7 にそろえる。再現しなかったときの結論は「この平板の条件では 2000 step 以内に再現しなかった」まで (ノズル固有とは言わない)。
  平板の起点は「既存の発達場」と呼ぶ (残差の収束と目的量の定常性は別に示す)。縦横比ではなく実測の Δτ と V/Δτ で比べる。

### 6.9 E1 (改訂) の結果 (2026-10-09、バイナリ lineF = 4fdc2c0e + typedef double、sha256 89ea9438…、実効の `implicitSolvePrecision` 0 = LHS は float で組む)

- **事前の確認** (`e1_compare.py`、1 step 目の書き出し `run_0307〜0310_e1dump_*`): マスク 7 と 5 は D が完全一致、rhs (sweep 0) の差 1.2e-12、K の行 0〜3 は完全一致、差は行 4 (熱伝導) だけ
  (最大 5.6e3) → 7/5 の比較の前提は成り立つ。マスク 0 と値 0 + キー 7 は K が完全一致、D は行 4 だけ最大 4.9e-4 (値 −4050.89、相対 1.2e-7。4050 付近の float32 の ulp は 2.44e-4 なので**約 2 ulp**。当初「1 ulp」と書いたのは誤り、codex 指摘) で、
  **事前に置いた許容 1e-12 には外れた**。原因は LHS を float で組むこと (ISP 0) と、熱伝導の D の計算順序が 2 つの経路で違うことで、許容を double の前提で置いた設定の誤り。
  基準どおりなら「無効」だが、差は float の丸めだけと記録する (許容を緩めるのは事後の変更)。
- **A/B**: A = マスク 7 (`case/45.isobutane_m6_d155/run_0305_e1_m7`) は **150 step で非有限** (DIVERGED)。序盤の推移は値 3 (run_0301) と同じ (step 1: ρ 1.37・ρv 2.6 倍、50: ρ 489 倍) → 破綻を再現。
  B = マスク 5 = 熱伝導の近傍 K だけを抜く (`run_0306_e1_m5`) は **2000 step 有限**、最大/開始 ρ 650・ρv 1060 倍、3 倍を超えた最初は ρ 141・ρv 120 step (値 0 + キー 7 の run_0304 と同じ)、
  末尾 500 step の傾き × 500 は ρ −1.71 桁、check_convergence (0〜1999) NOT CONVERGED (still converging、低下 3.0 桁)。200 step ごとの場の非有限・非正は後で確かめる。
  → 観測は事前の分岐 1 の形 (A 破綻・B 有限) だが、**事前確認 (マスク 0 と値 0 + キー 7) が許容に外れたので、登録上は「比較は無効」のまま**。
  「熱伝導の近傍 K を除くと破綻を避けられる」は、その解釈を支持する**探索的な観測**として扱う (codex 諮問 [記録](../../notes/reviews/2026-10-09-line-viscous-heatK-diagnose.md))。
- **監査 (2026-10-09)**: 両 run とも forge の sha256 89ea9438…、restart は種とビット一致 (9 量)、マスクはログで確認 (マスク 5 は起動時に印字、7 は既定)。
  マスク 5 の最終の場 (res_2000) に非有限・非正の ρ・P・T は無い。**中間のスナップショットは片付けで消していたので途中の場は監査できない** (私の手順の誤り。以後は監査が済むまで消さない)。
  forge_run.log に EOS の床の記録は無いが、forge が床の到達をログに出すかは未確認。マスク 5 の判定は区間 0〜1999 で NOT CONVERGED (still converging、低下 3.0 桁) で、収束とは扱わない。

### 6.10 事前登録: 熱伝導の K の密度の列だけの A/B (2026-10-09、codex 諮問 [記録](../../notes/reviews/2026-10-09-line-viscous-heatK-diagnose.md))

- マスクのビット 8 (診断) = 熱伝導の近傍 K のうち密度の列 K[4][0] (−κ(γ_j/c_p,j)(e_j − ½|u_j|²)/ρ_j) だけを外す。運動量・エネルギーの列、粘性の仕事の K[4][0]、D、壁の拘束は残す。
  これは全体の K が持つ「速度と温度を保つ補正 δρ(1,u,v,w,E) を消す」性質を壊すので、診断専用 (安定になってもその列の微分が誤りという証拠にはしない)。
- **A = マスク 7、B = マスク 15 (= 7 + 8)**。値 3・キー 5・方向別・上限なし・cfl 4・緩和 0.7・sweep 5・ISP 0、run_0183 の res_100000 から、同じ新バイナリで逐次に最大 2000 step、
  場は 200 step ごと (**監査が済むまで消さない**)、全残差は毎 step、序盤 200 step は 112 節点の帳簿。run 名 `run_0311_e1b_m7`・`run_0312_e1b_m15`。
- **起動の条件 (ゲート)**: 1 step 目の書き出し (マスク 7・15・5) で、(i) D が 7・15 で完全一致、(ii) 拘束のフラグが一致、(iii) rhs (sweep 0) の差 ≤ 1e-10 (再実行の差の桁)、
  (iv) K は行 0〜3 と行 4 の列 1〜4 が 7・15 で完全一致、行 4 の列 0 は 15 と 5 で完全一致 (どちらも熱伝導の密度の項を含まない、同じ演算の順序)、7 と 15 では違う。
  どれかが外れたら比較の run を起動しない (`e1_compare.py` の終了コードで止める)。
- **分岐**: A が既知の破綻を再現し、B が 2000 step 有限・非物理値なし → 密度の列の除去がこの条件・期間の破綻の回避に十分 (H1 の限定版を支持。式の誤り・TP の基準の問題・恒久の修正の妥当性は証明しない)。
  両側で同じ前駆の破綻 → 「密度の列が必要」を棄却 (H2/H3 を確定したとは言わない)。A が再現しない、または B だけ別の破綻 → 判別不能 (旧 run を対照に置き換えない)。
  1 step 目の書き出しで、熱伝導の K の作用を密度・運動量・エネルギーの列に分けた量 (K_heat ΔQ の各項と和) を A について出す (過渡の診断量)。
- マスク 5 (熱伝導の K を全部外す) を加速の候補として評価するのは、この A/B と監査の後。評価するなら、上限なし・上限 50・既存の方式について、point の仕上げを含む同じ終了条件までの総壁時計で比べる。

### 6.11 §6.10 の結果 (2026-10-09、バイナリ lineG = 846727be + typedef double、sha256 c0b64905…、ISP 0)

- **ゲート通過** (`run_0313〜0315_e1bdump_*`、`e1b_gate.py`): D・拘束のフラグが 7・15 で完全一致、rhs の差 1.9e-11、K の行 0〜3 と行 4 の列 1〜4 は 7・15 で完全一致、
  行 4 の列 0 は 15 と 5 で完全一致・7 と 15 で違う。A の 1 step 目 (sweep 0) の熱伝導の K の作用の二乗平均は密度の列 4.48・運動量の列 0.010・エネルギーの列 3.58・和 7.62 (打ち消さず足し合う)。
- **A/B**: A = マスク 7 (`case/45.isobutane_m6_d155/run_0311_e1b_m7`) は 122 step で非有限 (前回の run_0305 は 150 step、同じ入力でも再実行で破綻の step が変わる)。
  B = マスク 15 = 熱伝導の K の密度の列だけを外す (`run_0312_e1b_m15`) は **566 step で非有限**。前兆はどちらも ρv が先 (step 2 で A 5.3 倍・B 3.3 倍)、B は最後に k の残差が 1.3e6 倍。
  保存した場 (A は res_0 のみ、B は res_0・200・400) に非有限・非正の ρ/P/T は無い。
  → **事前の分岐 2: 「熱伝導の K の密度の列が破綻に必要」は棄却** (密度の列は破綻を早めるが、無くても壊れる)。H2/H3 を確定したとは言わない。
  **訂正 (codex 諮問 [記録](../../notes/reviews/2026-10-09-line-viscous-heatK-energycol-diagnose.md))**: 「残る候補はエネルギーの列」は強すぎる。単一の列の除去は人工的な結合を作り、
  熱伝導の K の作用は全列の和で κδT になる。「密度の列なしでも壊れた」は記録どおりだが、「元と同じ破綻の機構に密度の列は不要」は、発生の場所・流束の項・床の到達の順序の照合が足りない。

### 6.12 次の一手の候補と調査の区切り (2026-10-09、codex 諮問。ユーザは作用素の切り分けを選択 → §6.13)

- **作用素の切り分け (推奨の一組)**: 同じバイナリでマスク 7 / 5 (熱伝導の K の有無だけ) を新しい run で逐次に最大 2000 step。初期の状態と、A で最初の増幅が出た床到達前の状態を保存し、
  その同じ状態から両方式の補正 p を計算して、実残差の方向微分 J_true p ≈ {R(Q+εp) − R(Q−εp)}/(2ε) (ε と ε/2 で安定、差分の信号が再評価のノイズを上回る) と
  近似 Jacobian の作用 J_approx p を同じ p で比べる (熱伝導・応力と仕事・対流・ソースの内訳)。緩和後の実補正 p で R(Q+p) − R(Q) と J_true p も比べる。
  診断の閾値 (事前登録): 差分の再現誤差 1 % 以内、10 % 以上の乖離を調査の対象とする。分岐: 微小補正から不整合があり A の前兆の場所・成分に強く出る → 第 1 仮説 (近似 Jacobian と実残差の食い違い) を支持、
  微小補正では整合し実補正でだけ乖離 → 第 2 仮説 (有限振幅・物性の変化・床) を優先、両方整合 → H2 を支持しない、A が再現しない・差分がノイズに支配される → 判別不能。
  摂動した場で実残差を評価する道具 (forge で場を 0 step から評価し全保存量の残差を書く経路) が新たに要る。
- **区切り**: 上の一組で得た範囲を記録して原因の調査は一区切りにする (codex 推奨)。真因の完全な同定を速度の評価の必須条件にはしない。
- **マスク 5 を加速の候補として評価するなら**: 新しい比較で 2000 step の全残差・全域の非有限/非物理値・床の到達を監査し、実効の設定・restart・行列の差を確かめてから、
  上限なし・上限 50・既存の方式について、point の仕上げを含む同じ終了条件までの総壁時計で比べる。共通の条件に届かない、または仕上げ込みで速くなければマスク 5 の評価は終える (θ の調整を際限なく続けない)。
- run の監査の注意: e1b.sh は計算後の echo が関数の終了コードになり計算の失敗を伝えない。次の比較では新しいディレクトリを必須にし、計算の終了コード・書き出しの時刻と実行の識別子を確かめる。

### 6.13 事前登録: 方向微分による作用素の確認 (2026-10-09、ユーザの選択「方向微分で作用素を確かめる」、§6.12 の一組)

- **目的**: 値 3 (薄層の粘性・熱伝導の Jacobian + スカラー対角) のライン行列が、実残差の Jacobian の作用をどれだけ再現しているかを、
  同じ状態・同じ方向で直接測る。マスク 7 (熱伝導の近傍 K あり) と 5 (なし) の違いが実残差に近づく向きか遠ざかる向きかも見る。
- **状態**: S0 = `run_0183` の `res_100000` (= `run_0313`/`run_0315` の 1 step 目の状態。res_0 と ρ はビット一致、ρE は forge が起動時に掛ける等温壁のピンで ≤ 7.5e-9 の差、毎回同じ)。
  S1 = 新しいマスク 7 の run (`run_0316_jp_m7_s1`、同じ設定・バイナリ lineG・20 step・10 ごとに出力) の 20 step 目。
  run_0311 では 20 step で rms ρv が初期の約 35 倍、破綻 (122 step) の約 1/6 の位置 = 増幅が進んでいて床に届く前の状態。
  S1 で非有限・非正があれば S1 は判別不能として打ち切る。S1 の 1 step 目の書き出し (`run_0318_jpdump_s1_m7`・`run_0319_jpdump_s1_m5`) は
  S0 の 5 本のラインに、S1 の出力で |res_ρv| が最大の節点を含むラインを足す。7 と 5 の書き出しが同じ状態・同じ rhs (相対 ≤ 1e-9)・D の行 0〜3 が一致することを確かめてから進む。
- **方向**: p7 / p5 = マスク 7 / 5 の書き出しの最後の sweep の補正 ÷ 緩和 0.7 (= 緩和前のライン解。書き出しの行列で解き直すと残差 1e-14)。書き出したライン上だけ非零。
  ρ・ρu・ρE に s·p を足し、ρk・ρω は保存量のまま固定 (ソルバの 1 step でも乱流は流れの更新と別)、ρY_k は Y_k を固定 (ライン行列は組成を凍結した 5×5)。
- **差分**: 中心差分 J_t p = −{R(Q+εp) − R(Q−εp)}/(2ε)、ε = 1e-2 (S0 で max |δρ/ρ| ≈ 7e-6、|δu| ≈ 6e-4 m/s)。ε/2 = 5e-3 と比べた相対差 (ライン・行ごと) ≤ 1 % を「差分が有効」とする。
  同じ Q を 2 回評価した差 (atomicAdd の順序の揺れ) を方向微分の尺度に換算して記録する。
  1 % を外れた行は ε = 1e-1 と 5e-2 で 1 回だけ再試行し、なお外れればその行は判別不能 (ノイズまたは微分できない切替)。
- **比較**: J_a p = (D − V/Δτ I) p − K_prev p_{k−1} − K_next p_{k+1} (書き出しの D・K・dt_local・体積; D の時間項は float で組むので相対 6e-8 の丸めは残る)。
  ライン・行ごと (拘束の行を除く) に相対差 ‖J_a p − J_t p‖/‖J_t p‖、大きさの比、cos。作用素 7・5 をそれぞれ方向 p7・p5 に掛ける (2×2)。
  有限振幅: 実際に掛けた補正 0.7p について −{R(Q+0.7p) − R(Q)} と 0.7 J_t p の相対差。
  エネルギーの行の実の内訳: `FORGE_WI_FORCE_DIAG` の `wi_eheat` (内部面の熱伝導) と `wi_ework` (内部面の粘性の仕事) の方向微分、残り = 対流・境界面・ソース。運動量の行の内訳は取らない。
  出力の残差 (`res_*`) が書き出しの rhs_s0 と同じ R(Q) であることを q0 で確かめる (rhs は float で持つので相対 1e-6 程度は丸め)。
- **分岐** (§6.12 のとおり。「強く出る」= 10 % 以上の乖離が ρv の行 (A の前兆の成分) か、7 と 5 で違うエネルギーの行に、少なくとも 1 本のラインで出る):
  (a) 微小補正 (ε) で 10 % 以上の乖離が強く出る → 第 1 仮説 (近似 Jacobian と実残差の食い違い) を支持。7 と 5 のどちらが実残差に近いかも書く。
  (b) 微小補正では全行 10 % 未満で、有限振幅 (0.7p) でだけ 10 % 以上 → 第 2 仮説 (有限振幅・物性の変化・床) を優先。
  (c) 両方 10 % 未満 → H2 (近似 Jacobian の誤り) を支持しない。破綻はライン上の作用素の不一致でなく、反復の側 (ライン外との Jacobi の結合・擬似時間刻み) に残る。
  (d) 差分が有効でない行が判定対象の半数以上、または S1 が得られない → 判別不能。
  S0 と S1 は別々に判定して両方書く。
- **対象外**: ライン外の隣接節点への作用 (書き出していない K_offline) は比べない (応答の大きさだけ記録)。番号が隣り合うラインが含まれると互いの摂動が混ざりうるので、そのときは明記する。
- **run**: `run_0316_jp_m7_s1` (S1 の生成)、`run_0317_jp_s0_*` (S0 の 1 step の評価 12 本: q0・q0b・p7/p5 × ±ε・±ε/2・0.7p)、`run_0318`/`run_0319_jpdump_s1_*` (S1 の書き出し)、
  `run_0320_jp_s1_*` (S1 の評価 12 本)。`jprobe.sh` は各 forge の終了コードと出力の時刻を確かめ、外れたら全体を止める。評価 run の全節点の `res_*`・`wi_*` は
  `_jprobe/npz/` に抜いてから h5 を消す (q0 の res_1 は残す、摂動した入力の場は sha256 を記録して消す)。
- 結果の解釈は codex に諮ってから §6.14 に書く。どの分岐でも原因の調査はこの一組で区切る (§6.12)。
- **事後の追加 (2026-10-09、結果を見た後)**: ε = 1e-2 では差分の再現が S0 のほぼ全行・S1 の大半の行で 1 % を外れた (同じ入力の再評価の差は 1e-8 以下)。事前登録どおり ε = 1e-1 (ε/2 = 5e-2) で 1 回だけ再試行する (`run_0321_jp_*r_*`)。あわせて、判定には使わない妥当性の診断として S0・方向 p7 で ε = 1e-6 の応答を測る (`run_0322_jp_s0p7n_*`; 真の応答は ε = 1e-2 の 1e-4 倍なので、それより大きく出る分は入力に依存する丸めの揺れ)。残差の経路には double のビルドでも float の箇所がある (輸送物性の表引き `transport_mix_Y_tab` は T を float にして μ・λ を float で返す)。**訂正 (2026-10-09)**: ここで挙げた「リミッタの `fabsf`」(`limiter_d.cu:432`) は診断 `g_limDiag` の中だけで本計算の経路ではない。代わりに見落としていた経路として、TP の SLAU の面エンタルピー (`convectiveFlux_slau_d.inc.cuh:401-402` `thermo_h_mix_f`、double のビルドでも T・Y を float にして float を返す、`thermoFloat` と独立) と、`thermoFloat` 1 (既定) の γ・音速 (`dependentVariables_d.cu:186,225`) がある (codex 諮問 2026-10-09)。

### 6.14 §6.13 の結果 (2026-10-09、バイナリ lineG、codex 諮問 [記録](../../notes/reviews/2026-10-09-line-viscous-jprobe-diagnose.md))

- **正式判定: S0・S1 とも分岐 (d) 判別不能**。差分の再現 ≤ 1 % を満たした行 (ε = 1e-2 で合格した行はそのまま、外れた行だけ ε = 1e-1 で再試行) は
  S0 で 20 行中 3 行 (ライン 2183 の mass・momx、ライン 1640 の mass)、S1 で 24 行中 4 行 (ライン 2183 の mass・momx・momy、ライン 1640 の mass)、方向 p7・p5 とも同じ。
  同じ入力を 2 回評価した差 (atomicAdd の順序) は ε = 1e-2 の応答の 1e-13〜1e-8 で、差分の失敗の原因ではない。出力の残差は書き出しの rhs_s0 と相対 2e-8 で一致 (rhs は float)。
  合格した少数の行での近似と実の相対差は 0.20〜0.78 (作用素 7 と 5 はこれらの行で同一)。
- **測定が成り立たなかった原因は未同定**。残差の応答 ΔR(s) = R(Q+sp) − R(Q) を s = ±5e-7〜±1e-1 の 13 点で見ると (S0・p7、`_jprobe/curve.py`・`jumps.py`)、
  大半の節点で s ≤ 1e-6 では比例する一方、s = 5e-3〜1e-1 では線形から外れる量が s に比例せず (mass の行で 90 % 点 5e-7 → 1.4e-6、s を 20 倍にして 2.7 倍)、
  s = 1e-6 でも一部の節点で大きく外れる。リミッタの min・頭打ちは入力が連続なら値の不連続を作らないので、これだけでは説明にならない (codex)。
  候補として残すのは float の評価経路 (TP の面エンタルピー `thermo_h_mix_f`、`thermoFloat` の γ・音速、float の輸送物性の表引き) と、再構成・境界などの分岐。
  `jumps.py` が数えたのは「最小幅の中心差分から 10 % 以上外れた節点」で数学的な不連続点ではない。
- **探索 (事前登録の外、ε = 1e-6、判定には使わない)**: 差分が再現したライン (S0 の 2183、S1 の 65・1640・2183・12 の p5 の一部の行) では、
  **安定化項を含む近似作用素** (値 3 は D にライン面のスカラー 2ν·δ/dcc も足している; `jprobe.py` が引くのは時間項だけ) と評価した残差の方向微分の相対差は
  mass・運動量の行で 0.22〜1.5 (熱伝導の K と無関係な行も含む)、ライン 2183 のエネルギーの行は大きさの比 |J_a p|/|J_t p| が 0.028 (S0)・0.024 (S1)。
  エネルギーの行の実の応答に対する内部面の熱伝導・粘性の仕事の方向微分のノルム比は熱伝導 0〜0.044・仕事 ≤ 0.001 (足して 1 になる寄与率ではない)。
  これは「近似 Jacobian の欠陥」や破綻の原因の証拠ではない (安定化項・1 次の流束分離・凍結した物性・再構成の差を含む)。
  0.7p の試験は選んだライン上だけの補正の有限振幅の応答で、全域の 1 step の予測誤差とは別物。
- **主張しないこと**: ε = 1e-6 の結果で (d) を (a) に置き換えない。跳ぶ節点を事後に除いて「微分が成り立った」と扱わない。マスク 7 の破綻を熱伝導の K の誤りに帰さない。
- codex の採否 (Major 4・Minor 2、全件採用): H-a の原因の説明は却下し「測定不成立」と「原因未同定」を分けて記録 (上)、H-b は呼び方を限定して記録 (上)、
  H-c は面エンタルピーの float の経路の A/B 一組を §6.15 に事前登録、H-d の「float では説明できない」は面エンタルピーの経路を見落としていたので撤回、
  集計 (S1・p5・ライン 12 の momx の再現差は 2.5e-3 で 1e-3 を超える; 内訳はノルム比) を訂正。

### 6.15 事前登録: TP の面エンタルピーの精度だけを変える A/B (2026-10-09、codex 諮問の判別 A/B)

- **変える一点**: SLAU の TP 多成分の面エンタルピー (`convectiveFlux_slau_d.inc.cuh` の `thermo_h_mix_f(spf, ...)`) を、診断の環境変数
  `FORGE_DIAG_FACE_H_DOUBLE=1` のときだけ同じ NASA 係数・datum の double 版 `thermo_h_mix(sp, ...)` (温度・組成を double のまま) で評価する。
  既定 (変数なし) の経路は変えない。単成分の経路 (`thermo_h_mass_f`) は対象外。A = 切替なし、B = 切替あり、同じ新バイナリ。
- **固定**: S0 (`run_0183` の `res_100000`)、方向 p7 (`run_0313` の書き出し、作り直さない)、摂動の支持 (書き出した 5 本のライン)、近似作用素 (`run_0313` の D・K)、設定は §6.13 と同じ。
- **評価**: 各側で q0・q0 の再評価・±1e-6・±5e-7 の 6 本 (1 step)、計 12 本。run 名 `run_0323_jph_{a,b}_{q0,q0b,pe,me,ph,mh}`。対象は事前に **ライン 2183 のエネルギーの行 (拘束の行を除く)** に固定。
- **有効の条件**: 両側でその行の差分の再現 (ε/ε/2) ≤ 1 %、同じ入力の再評価の差 ÷ 信号 ≤ 1e-3。
- **分岐**: B で J_t p が A から 10 % 以上変わる → 面の熱力学の精度に依存することを支持。さらに近似作用素との相対差が A の半分以下 → H-c (エネルギーの行の大差) の説明としても支持。
  変化が 10 % 未満 → この条件での主要因説を棄却。どちらかの側で有効の条件を満たさない → 判別不能のまま終える。
  **どの分岐でも探索はここで終え** (§5.1 #5m)、面エンタルピーを double にする変更を本線に入れるかは別の判断にする (速度・回帰への影響を測っていない)。
  面ごとの Δh は測らない (計測の経路がない)。mass の行の不成立は面エンタルピーでは説明しない (対象外)。
- ビルド: 新しい commit + typedef double を AWS で別の作業ツリーにビルドし、`cold_cfl.py` の `ALT_BINARIES` に `lineH_fp64` として sha256 を登録する。

### 6.16 §6.15 の結果 (2026-10-09、バイナリ lineH = 5197e00e + typedef double、sha256 561813b3…、codex 諮問 [記録](../../notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md))

- `run_0323_jph_{a,b}_*` (A = 既定、B = `FORGE_DIAG_FACE_H_DOUBLE=1`)、S0・方向 p7・ε = 1e-6、ライン 2183 の自由なエネルギーの行 (`_jprobe/fh_judge.json`):

  | | A (float の面エンタルピー) | B (double) |
  | --- | --- | --- |
  | 差分の再現 ε/ε/2 | 9.5e-6 | 7.1e-4 |
  | 同じ入力の再評価 ÷ 信号 | 7.9e-10 | 2.0e-5 |
  | ‖J_t p‖ | 148.9 | 5.28 |
  | 近似作用素 7 との相対差 (大きさの比、cos) | 0.975 (0.028、0.89) | 0.389 (0.80、0.93) |

  両側とも有効の条件を満たし、‖J_t^B − J_t^A‖/‖J_t^A‖ = 0.966 → **事前登録の分岐: 面の熱力学の精度に依存することを支持、近似作用素との相対差が半分以下 (0.399 倍) なので H-c の説明としても支持**。
  主張は S0・p7・ライン 2183・自由なエネルギーの行・ε = 1e-6 に限る。§6.13 の正式判定 (S0・S1 とも (d)) は上書きしない。
- 記録の仕方 (codex): 「約 36 倍の大きさの隔たりには面エンタルピーの評価の精度が強く寄与した」。切替は温度だけでなく組成・係数・演算・戻り値の精度も変えるので、
  「h が凍結して ṁ·∂h が抜けた」という機構は面ごとの Δh を測っていないため未確定。B でも相対差 38.9 % が残る (ライン面のスカラーを D から除いても 0.389 のまま = この対象ではスカラーが主因ではない)。
- 切替なし (A) の q0 の残差は旧バイナリ lineG と同じ入力の再評価の差の桁で一致 (res_roe 最大 2.3e-10)。B の q0 は A の q0 から ‖R_B − R_A‖/‖R_A‖ = 5.5e-3 (エネルギー、他の行は ≤ 6e-11)。
  これは「同じ状態で残差の評価が変わる」ことまでで、残差のノルムが 0.55 % 増減した意味ではない。**S0 (`run_0183`) は収束済みではない** (`check_convergence` NOT CONVERGED、stalled/plateau、rms_roOmega RISING; `_band_ab/cold_pair/gates_aws.json`)。
  ブリーフの「収束した S0」は誤り。残差の床・ライン陰解法の破綻との因果は測っていない → point 仕上げでの A/B を [time_integration-implicit-thermal-jacobian](../active/time_integration-implicit-thermal-jacobian.md) §5.1 #5 に登録 (本 plan では回さない)。
- **無効にした数字**: B 側の「有限振幅」「実補正の予測差」は pp に旧 float 経路の `run_0317_jp_s0p7_pp` を使ったので A と B の評価の経路が混ざっている。掲載しない (正式判定は参照していない)。
- **やらないこと** (codex): 面エンタルピーの double を既定にする、値 2・3 を本線に戻す、熱伝導の K の列の削除を続ける、datum の大きさを測らずに原因とする。

### 6.17 閉じるときの処置 (2026-10-10、諮問 [記録](../../notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md) の「plan を閉じる際」の要求に沿って書く)

**結論**: 値 2・3 は本線に採らない。既定 (`lineViscCoupling: 0`) を維持する。値 2・3 と診断の切替は、診断用としてコードに残す。

**総時間の結果** (速度 plan [time_integration-line-implicit-speed](../active/time_integration-line-implicit-speed.md) の測定の引用):

- §6.15 (終わりの基準は水準だけ、事後の集計): 粘性入りのライン L5 (値 3・`FORGE_LVC_TERMS=5` = 熱伝導の近傍 K なし・上限なし、`case/45.isobutane_m6_d155/run_0353_m9_L5`) は 115000 step・約 1.18 h。
  値 0・キー 5・上限 50 のライン L0 (`run_0223` → `run_0224` → `run_0252`、古いバイナリ) は 120000 step・約 1.16 h。差は約 70 s で、出力の間隔 1 回ぶん (約 180 s) より小さいので**判別不能**。
- §6.14 (その後に外した E2 の基準): L5 1.38 ± 0.08 h、L0 1.49 ± 0.08 h で、こちらも判別不能。
- §6.20: 値 0 の全長ライン B0 (LAYOUT2 の単価) は 125000 step・1.119 ± 0.049 h。
- 速度 plan §6.15 は「本番は値 0。粘性入りは診断用の環境変数に頼る経路なので本番に置かない」とした。本 plan の結論もこれに従う。
- 留保 (codex result m4): 上の時間は専有の単価で換算した到達時間の推定で、収束解どうしの比較ではない。L5 のライン段・point 段とも `check_convergence` は NOT CONVERGED、到達の窓の θ_r は TRANSIENT-UNSETTLED だった (速度 plan §6.14)。§6.20 の B0 は「2 出力連続」の到達条件で数えたので、§6.15 の L0・L5 と到達条件が違う。

**試験ごとの処置** (不合格・判定不能・未実施を消さない):

| 試験 | 結果 | 処置 |
| --- | --- | --- |
| U-J | PASS (§6.0)。ただし行列全体の最大値で正規化した判定 | 薄層の流束モデルに対する単体照合で、差を D・K の全体の最大値で割った判定の合格にとどまる。登録した列ごとの相対誤差 ≤ 1e-6 は未確認 (再開の課題)。実残差の Jacobian との整合も示していない (§6.14 は判別不能) |
| U0 | 未実施 | 未実施のまま閉じる。**代わりの証拠** (事後に見つけたもの。事前登録の条件とは違う): 速度 plan §6.0 の B1 の短期の一致 (`run_0276`〜`run_0281_abB_*`、`_band_ab/cold_pair/AB_abB.json`) は、5ab83056 より前の 4d394a71 (forge の sha256 d8b06ebc…) と 5ab83056 を含む 64ed2cd6 (7a7e9eb9…) を、値 0・キー 5・方向別・上限 50 で比べている。20 step の 8 量の相対 RMS・最大絶対差は、同じバイナリの再実行の最大の 3 倍以内だった (再実行は各側 1 組だけ)。**範囲**: 値 0 のこの設定に限る。値 1・point・ISP 1 は確かめていない。**移管先**: 値 1 は [`procedures/solver-settings.md`](../../procedures/solver-settings.md) に「5ab83056 以降の不変を確かめていない」と書き、使う前に旧バイナリとの短期の一致を取る。値 0 の他の経路は main へ統合するときの回帰で見る |
| U2 | 不合格 (0.05 K まで 950 step、基準は 250 step 以内) | 記録のまま。値 0 が壊れる cfl 50 で値 2 は温度の解析解誤差が 5000 step で 1.4e-6 K まで下がった、という観測も残す (残差の収束は未確認) |
| U3 | Couette は判定不能、Poiseuille は不合格 | 記録のまま。せん断・粘性の仕事の項は未検証 |
| FP32 | 未実施 | 未実施のまま閉じる。**適用範囲**: 値 2・3 の動作は FP64 ビルド (typedef double) でだけ確かめた。FP32 ビルドでは起動も確かめていない。再開するときの最初の試験に入れる (下の再開の条件) |
| V-n1 | B 不合格 (29 step で非有限) | 記録のまま |
| V-n2 | 未実施 | 事前の条件 (V-n1 の合格) を満たさないので回さない。代わりに値 3・マスク 5 を速度 plan の総時間で評価し、判別不能だった (上) |
| §6.2・§6.3 (host の再解) | 解法の整合を確認、密度が主の準特異の仮説は棄却 | 記録のまま。「実装経路の誤りは無い」の表現は 2026-10-10 に訂正した |
| §6.4・§6.5 (値 3) | B 不合格 (十分性を棄却) | 記録のまま |
| §6.6・§6.7 (壁の拘束) | 分岐 2 | 記録のまま |
| §6.8・§6.9 (E1) | 事前確認が許容外のため比較無効 (探索的な観測) | 記録のまま。§5.1 #5g の旧要約「(分岐 1)」を訂正した |
| E2 (平板) | 未実施 (後回し) | 未実施のまま閉じる。準備の穴 (restart の終了コードを確かめずに起動できる) も直していない |
| §6.10・§6.11 (密度の列) | 分岐 2 | 記録のまま。除いても 566 step で非有限だったが、同じ破綻の機構にとって必要かは未確定 |
| §6.13・§6.14 (方向微分) | S0・S1 とも (d) 判別不能 | 記録のまま。§6.16 の結果で上書きしない |
| §6.15・§6.16 (面エンタルピー) | 限定付きの支持 (S0・p7・ライン 2183・自由なエネルギーの行・ε = 1e-6) | 記録のまま。機構 (ṁ·∂h の欠落) は未確定。B 側の有限振幅の指標は無効 |
| 残差評価の精度の監査 | — | [time_integration-implicit-thermal-jacobian](../active/time_integration-implicit-thermal-jacobian.md) §5.1 #5 へ移した (本 plan では回さない) |

**診断コードの残置の範囲** (コードは変えない):

- `lineViscCoupling` 2・3 と、その設定の検査 (`solverConfig.cpp`)。値 3 は値 2 にライン面のスカラー対角を足した診断用の値。
- `FORGE_LVC_TERMS`: 値 2・3 の項のマスク。ビット 1 = 運動量の D/K、2 = 熱伝導の近傍 K、4 = 粘性の仕事の D/K、8 = 熱伝導の K の密度の列を外す。既定は 7。
- `FORGE_LINE_DUMP_DIR`・`FORGE_LINE_DUMP_CALL`・`FORGE_LINE_DUMP_NODES`: ライン行列の書き出し。
- `FORGE_DIAG_FACE_H_DOUBLE`: TP の SLAU の面エンタルピーを double で評価する。
- 単体試験 `solver_density_cuda/tools/test_line_visc_jacobian.cpp`。case/45 の解析スクリプト `linedump_analyze.py`・`e1_compare.py`・`e1b_gate.py`・`jprobe.py`・`jprobe.sh`・`jprobe_fh.sh`。
- 既定の経路 (値 0、診断の環境変数なし) について確かめた範囲: 値 0 は上の B1 の短期の一致、`FORGE_DIAG_FACE_H_DOUBLE` なしは §6.16 の q0 (旧バイナリ lineG と、同じ入力の再評価の差の桁で一致)。
  ビット一致は確かめていない (forge は同じバイナリの再実行でもビット一致しない。thermal-jacobian plan §6.0 の V0)。
- 本番の設定に値 2・3 と診断の環境変数を書かない。速度 plan の L5 の腕は `FORGE_LVC_TERMS=5` に頼っていた。
- 残る不整合 (今回は直さない): `solverConfig.cpp` の `implicitThermalJacobian` と値 1 の併用を拒否するメッセージは「use lineViscCoupling 2」と値 2 を勧めている。値 2 は本線に採らないので、次にこのファイルを触るときに直す。

**再開の条件** (どれかに当たったら新しい plan を起こして再開する。本 plan は判断文書として accepted に残す):

1. 値 0 のライン + 方向別 dt の総時間が、ライン内で結合していないせん断・エントロピーのモード (§1) の遅さで律速されていると測定で示されたとき。
2. 残差の評価の精度を上げた経路 (面エンタルピーの double など) が既定になったとき。§6.16 の限定付きの支持があるので、§6.13 の作用素の確認をその経路で取り直す意味がある。
3. 再開したら、次を先に行う: 単因子の切り分け (運動量・仕事・熱伝導の K のどれが破綻に要るか、§6.7 の候補)、U-J の列ごとの照合、FP32 ビルドでの起動と U2・U3、値 1 の U0、E2 の準備の穴の修正。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| result | 2026-10-10 | [2026-10-10-time_integration-line-viscous-jacobian-result.md](../../notes/reviews/2026-10-10-time_integration-line-viscous-jacobian-result.md) | GO-with-changes, C0/M3/m1 | 全件採用 (根拠の箇所を確かめた): M1 U-J の試験は行列全体の最大値で正規化しており、登録した列ごとの判定は未確認 → §5.1 #2・§6.0・§6.17・methods を訂正し再開の課題に。M2 U2 の「収束」は温度の解析解誤差だけが根拠 → §6.0・§6.17・case/52 README を「誤差が下がった・残差の収束は未確認」に。M3 E1 は比較無効、密度の列は「除いても壊れる・同じ機構に必要かは未確定」→ §5.1 #5i・§6.17・case/45 README を訂正。m4 総時間の引用に「専有単価の推定・NOT CONVERGED・TRANSIENT-UNSETTLED・B0 は到達条件が違う」の留保を §6.17 に追加。B1 を U0 の補助証拠とする扱いは codex も可 (U0 未実施は維持) |
| plan | 2026-10-09 | [2026-10-09-time_integration-line-viscous-jacobian-plan.md](../../notes/reviews/2026-10-09-time_integration-line-viscous-jacobian-plan.md) | GO-with-changes, C0/M5/m2 | 全件採用: M1 等温壁の行を ΔT_w = 0 の拘束に (§4.1)、M2 入力経路と対応範囲の限定 (§4.3)、M3 単体照合 U-J・U3・FP32 (§6)、M4 判定の分離と `check_quasisteady` (§6)、M5 同じ水準までの壁時計 (V-n2)、m6 P は薄層の前処理行列と明記 (§4.1)、m7 SST の寄与を未分離と修正 (§4.2) |
| 諮問 (面エンタルピーの A/B) | 2026-10-09 | [2026-10-09-line-viscous-faceh-diagnose.md](../../notes/reviews/2026-10-09-line-viscous-faceh-diagnose.md) | 限定付きの支持として記録、機構 (ṁ·∂h の欠落) は未確定、S0 は収束済みでない、B の有限振幅の指標は無効、本線へ戻る、Major 3 | 全件採用: §6.16、§5.1 #5l〜#5n・#6、残差評価の精度の監査を implicit-thermal-jacobian §5.1 #5 へ |
| 諮問 (方向微分の結果) | 2026-10-09 | [2026-10-09-line-viscous-jprobe-diagnose.md](../../notes/reviews/2026-10-09-line-viscous-jprobe-diagnose.md) | (d) で閉じる、H-a の原因の説明は却下、H-b は呼び方を限定、面エンタルピーの float の経路の A/B 一組だけ、その後は本線、Major 4・Minor 2 | 全件採用: §6.14 に結果と訂正、§6.15 に A/B を事前登録、§5.1 #5l・#5m |
| 諮問 (密度の列の結果) | 2026-10-09 | [2026-10-09-line-viscous-heatK-energycol-diagnose.md](../../notes/reviews/2026-10-09-line-viscous-heatK-energycol-diagnose.md) | 「残るのはエネルギーの列」は強すぎ、作用素の切り分け (方向微分) を推奨、その後は区切る、スクリプトの終了コードの不備 | 全件採用: §6.11 を訂正、§6.12 に候補と区切りの条件を記録 |
| 諮問 (熱伝導の K) | 2026-10-09 | [2026-10-09-line-viscous-heatK-diagnose.md](../../notes/reviews/2026-10-09-line-viscous-heatK-diagnose.md) | 分岐 1 の確定は却下 (事前確認の外れ・監査未了)、ulp の訂正、ゲートの不備、H1 の「TP の負の e が原因」は却下、密度の列だけの A/B を推奨 | 全件採用: §6.9 を探索的な観測に、監査を記録、§6.10 を事前登録、e1_compare をゲートに |
| 諮問 (再開の点検) | 2026-10-09 | [2026-10-09-line-viscous-resume-diagnose.md](../../notes/reviews/2026-10-09-line-viscous-resume-diagnose.md) | E1 の組み方 (スカラーの二重計上)・E2 の準備の穴・平板の結論の範囲、Major 5・Minor 1 | 全件採用: マスクは熱伝導の K だけを切り替える、7/5 を先に、行列の事前確認、E2 は後回し |
| 諮問 (値 3 の結果) | 2026-10-09 | [2026-10-09-line-viscous-lvc3-result-diagnose.md](../../notes/reviews/2026-10-09-line-viscous-lvc3-result-diagnose.md) | 十分性の棄却は採用、「薄層の D/K が原因」の断定は却下 (壁の拘束と未分離)、P→I を先にしない | 全件採用: 値 0 でキー 5/7 の A/B を §6.6 に事前登録、その後は値 2・3 を保留、値 0 + キー 5 の記述と SVD の解釈を訂正 |
| 諮問 (弱いモード) | 2026-10-09 | [2026-10-09-line-viscous-thermal-weak-mode-diagnose.md](../../notes/reviews/2026-10-09-line-viscous-thermal-weak-mode-diagnose.md) | 第 1 仮説 (中) = 全行のスカラー対角を外したことで減衰が不足、Major 6・Minor 1 | 全件採用: 「温度が主」の解釈を訂正、拘束を消去した原始量の SVD、値 3 の A/B を §6.4 に事前登録、δT を緩和後の値に訂正 |
| 諮問 (発散) | 2026-10-09 | [2026-10-09-line-viscous-jacobian-divergence-diagnose.md](../../notes/reviews/2026-10-09-line-viscous-jacobian-divergence-diagnose.md) | 第 1 仮説 (中) = スカラー対角の除去で密度・圧力の補正が過大、Major 多数 | 全件採用: 実ライン行列の host 再解 A/B を §6.2 に事前登録、値 2 は本線に使わない、連続の行だけの正則化を先にしない、Q_w の転記を訂正 |
| plan (諮問) | 2026-10-09 | [2026-10-09-line-viscous-jacobian-and-v0-diagnose.md](../../notes/reviews/2026-10-09-line-viscous-jacobian-and-v0-diagnose.md) | 式・符号は採用 (独立モデルで数値微分と 3.6e-15)、壁拘束は条件付き、検証は修正 | 全件採用: 仕事の微分のモデルを固定 (§4.1)、V-n1 を同じ新バイナリの値 0/2 の上限なし 2000 step の A/B に (§6)、run_0211 の発散を ρ の人工拡散だけに帰さない。実残差との差の測定は未実施 (CFD の挙動で間接に見るだけ) |

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/timeIntegration_d.cu`、`solver_density_cuda/input/solverConfig.cpp`。`lineViscCoupling 0/1` の経路は変えない。

## 変更ログ

- 2026-10-09: 起票 (draft)。
- 2026-10-09: codex plan 段と諮問を全件採用して in_progress (§4.1・§4.3・§6 を改訂)。
- 2026-10-09: 実装 (5ab83056)、U2・U3・V-n1 (§6.0)、実ライン行列の解析 (§6.3)、値 3 の A/B (§6.5)、壁の拘束の A/B (§6.7)。値 2・3 は発散し、薄層の結合が要因と分かったところで保留。
- 2026-10-09: E1 (§6.9、比較無効・探索的な観測)、密度の列の A/B (§6.11、分岐 2)、方向微分による作用素の確認 (§6.14、(d) 判別不能)、面エンタルピーの精度の A/B (§6.16、限定付きの支持)。変更ログの前行の「薄層の結合が要因と分かった」は、§6.7 の帰属の範囲 (値 2・3 の早期の発散には薄層の D/K、またはそれと壁の拘束の相互作用が要る) に限る。
- 2026-10-10: 閉じるときの処置を §6.17 に記録 (試験ごとの処置・診断コードの残置の範囲・再開の条件・U0 と FP32 の扱い・総時間の引用)。§5.1 #4・#5g・§6.3 の古い要約を訂正。codex result 段 (GO-with-changes、C0/M3/m1) を全件採用して訂正し、status done (本線不採用・診断として残置) で `plans/accepted/` へ移した。コードは変えていない。methods/time_integration/implementation.md と procedures/solver-settings.md の `lineViscCoupling` の記述を同期。
```

## 出力形式 (この形のまま)

```
結論: <次にやる一手を 1 文で>
第 1 仮説: <内容>  確度: <高/中/低>
  根拠: <ファイル:行 / run パスと数値>
  反証条件: <何が観測されたらこの仮説は誤りか>
第 2・第 3 仮説: <あれば 1 行ずつ>
判別 A/B: <変える設定 1 点、回す長さ、見る量>  → A なら … / B なら …
やらない方がよいこと: <呼び出し側が取りそうな誤った一手>
呼び出し側の前提への異議: <ブリーフの枠組み・除外判断・指標の定義で受け入れなかったものと理由。無ければ「無し」>
不足情報: <あれば>
```
設計判断・採否を諮られた場合は、上の前に「採否表 (指摘ごとに 採用/却下/要再検証 と理由)」を置いてよい。
