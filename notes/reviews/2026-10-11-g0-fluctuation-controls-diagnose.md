# codex 諮問 (diagnose): g0-fluctuation-controls

- **brief**: [`notes/reviews/briefs/2026-10-11-g0-fluctuation-controls.md`](../../notes/reviews/briefs/2026-10-11-g0-fluctuation-controls.md)
- **plan**: [`plans/active/axisymmetric-freestream-hoop-gauge.md`](../../plans/active/axisymmetric-freestream-hoop-gauge.md)
- **date**: 2026-10-11
- **commit**: `e979a5c9` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.9 min, rc=0
- **結論**: **G0 の新格子・同一 FP64 バイナリ・同じ後期保存場から、`mesh.axisSegmentRWeight` だけを 1/0 に変えた固定 60,000 step の A/B を行う。**
- **extra**: `plans/active/architecture-float-state-double-geometry.md`, `case/45.isobutane_m6_d155/hp7.sh`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 4 / Minor 1）

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| Major | **3 要因を区別する必要性は採用。直ちに `K0`・`OLD` を長時間投入する案は却下** | 比較の配置自体は有用だが、`OLD` 対 `run_0452` にはビルド以外に実行環境・再実行差・診断の有無が混ざる。旧腕は `FORGE_DIAG_COMMIT_LOSS=2500` を付ける一方、新腕は付けない（[v4.sh:23](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/v4.sh:23)、[hp7.sh:9](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/hp7.sh:9)）。**まず同じ保存場・新格子・同一バイナリからキーだけを変える A/B に絞る。** |
| Major | **停止規則で止めた各末尾の RMS 比較は却下** | 停止自体が θ_r の変動に依存するため、比較窓が結果に選ばれる。提案の 10〜14 万 step も、ブリーフでは `run_0452` が 12.5 万で停止しており共通窓にならない。`m9_watch` は 2 万 step の両端差を使う（[m9_watch.py:116](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/m9_watch.py:116)）。**揺れの診断は固定長・固定窓とし、到達で止めない。** |
| Major | **「RMS 比 > 2 なら原因が効いた」は要再検証** | 4 万 step・間隔 2,500 では差分は 16 個程度。今回の FP64 構成の再実行差も未提示であり、2 倍に統計的な裏付けはない。**2 倍は暫定の実用効果量として登録し、再現性が確立した因果判定とは呼ばない。** 生系列の振幅・ドリフトも併記する。 |
| Major | **§6.20 の停止比較は採用、因果解釈は却下** | §6.11 は初めから「停止規則のもとでの比較」であり、非定常を許している。一方、§6.20(a) の「差が消えれば hoop と float の組み合わせが原因」、(b) の「差が半分以上残れば hoop は主因でない」は、格子・ビルド・停止時点も変わった比較からは導けない（[float plan:629](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:629)）。**フラグの数値判定は保存し、原因の断定を外す。** |
| Minor | **差分の負相関を雑音の証拠とする解釈は却下** | ブリーフの相関 −0.1〜−0.6 だけでは機構を識別できない。独立な変動 εₙ でも、差分 εₙ−εₙ₋₁ の隣接相関は −1/2 になる。**元系列の自己相関と振幅を確認し、「出力間差分が負相関」と記述する。** |

結論: **G0 の新格子・同一 FP64 バイナリ・同じ後期保存場から、`mesh.axisSegmentRWeight` だけを 1/0 に変えた固定 60,000 step の A/B を行う。**

第 1 仮説: **G0 では、キー 1 による面ベクトルと `ss` の変更が、後期反復の大きな Q_w の揺れを維持している。** 確度: **低**
  
  根拠: キー 1 は `sx..sz = W_f`、`ss = ‖W_f‖` に置き換え、キー 0 は従来の半径重みを使う（[variables.cpp:782](/home/sano/work/forge-integ-1005/solver_density_cuda/variables.cpp:782)）。したがって自由流の圧力釣り合い以外にも影響しうる。ブリーフの `case/45.isobutane_m6_d155/run_0483_hp7_k1/` の Q_w 差分 RMS 0.077 % は旧腕の 0.006 % より大きいが、**この数値は未再計算で、キーへの帰属は未確認**。G1/G1x のキー 0 にも揺れがあるという報告から、全格子共通の「キー 1 だけが原因」という説にはしない。

  反証条件: 下記の固定窓で、キー 1 の対照に揺れが再現しているにもかかわらず、キー 0 にしても Q_w の揺れが半減しない場合、「この状態・期間でキー 1 が揺れを 2 倍以上に維持する」という仮説は支持されない。**過去の状態形成への関与まで除外しない。**

第 2 仮説: **新しい境界半割面の幾何など、キー 0/1 に共通する変更が揺れを維持している。** 未確認。半割面の変更はキー 0 にも効き、勾配・再構成の幾何にも及びうる（[hoop plan:154](/home/sano/work/forge-integ-1005/plans/active/axisymmetric-freestream-hoop-gauge.md:154)）。

第 3 仮説: **長い過渡と停止時点・標本化の違いが、新旧の揺れの差を増幅して見せている。** 未確認。ビルド由来の差も除外済みではない。

判別 A/B:

- **起点:** `case/45.isobutane_m6_d155/run_0483_hp7_k1/` の取得時に保存されている最新の完全な出力を、step・ハッシュを固定して使う。両腕を新しい run に作り、`restart_field.py` で同じ 9 保存量を移す。`run_0483` 本体の停止規則は変えない。
- **変更点:** A = キー 1 の継続、B = キー 0。新格子・FP64 バイナリ・BC・B0 の実効設定・環境変数・抽出器を共通にする。キーの実効 ON/OFF と初期保存量の一致を確認する。
- **長さ:** 60,000 step 固定、2,500 step ごと。最初の 20,000 step を切替応答として分け、**20,000〜60,000 の 17 点**を主窓にする。全出力を残す。到達による停止や結果を見た窓の延長はしない。
- **物差し:** 主量は符号つき Q_w、副量は θ_r(40/70/94)。各量について、`J = RMS(qᵢ − qᵢ₋₁) / |窓平均|`、線形トレンド除去後の RMS、最大−最小、ドリフトを記録する。分母が非有限・ゼロなら判定不能。前後の 2 万 step 窓も同じ式で集計する。
- **事前の分岐:** A の揺れが維持され、B/A の `J` とトレンド除去後 RMS が**全窓・前後半とも ≤ 0.5**なら、第 1 仮説を暫定支持する。半減が持続しなければ、この窓での「キーによる 2 倍以上の維持」は支持せず、共通変更・履歴の候補を残す。A 自体が静まる、窓ごとに結論が逆転する場合は判別不能。これは再実行による信頼区間ではない。
- 全残差・保存場の非有限を確認し、同じ区間の `check_convergence` と、対象系列の `check_quasisteady` の VERDICT を残す。**非定常でも反復の揺れは比較できるが、定常解の比較とは呼ばない。**

やらない方がよいこと: **`run_0483` が未到達だから停止閾値を緩めたり、20 万 step を越えて延長したりしない。** §6.20 は登録どおり到達または上限まで待つ。FP64 が上限でも未到達なら、停止時の精度比較は「判定不能」とする。判定器もその扱いになっている（[fl1_judge.py:102](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/fl1_judge.py:102)）。揺れの増加だけで hoop 修正を撤回することも勧めない。

呼び出し側の前提への異議:

- **「到達」は揺れの小ささを保証しない。** 実際の `m9_watch.py` の判定関数に、平均の ±0.1 % を交互に取る合成系列を渡したところ、差分 RMS **0.2000003 %**、両端ドリフト **0 %**で、**`REACH` at 25,000**となった。これは実 run の結果ではなく、停止規則の性質の確認である。
- **case/48 の抽出アーチファクトを今回の主因にはしない。** 現行 G0 抽出器は固定の `common_yb` を使い、座標は結果ファイルではなく入力 `nozzle.h5` から読む（[cold_xcheck.py:49](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_xcheck.py:49)、[cold_xcheck.py:107](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_xcheck.py:107)）。動的な `je` 切替という既知の機構は、このコードにはない。

不足情報: 対象 run の `m9_watch.json`、実効設定・バイナリの記録、残差、保存場、正式 VERDICT はこの作業ツリーに無く、AWS 上の RMS・非有限・定常性は独立検証できていない。したがって原因の確定はしない。**ファイル変更・forge 起動は行っていない。plan 未反映であり、採用する場合の反映先は hoop plan §4 と float plan §6.20。**
