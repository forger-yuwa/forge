# codex レビュー: tooling-nozzle-throat-monotone-r2 (result)

- **plan**: [`plans/active/tooling-nozzle-throat-monotone-r2.md`](../../plans/active/tooling-nozzle-throat-monotone-r2.md)
- **stage**: `result` (diff base `c821b71e`)
- **date**: 2026-10-07
- **commit**: `6970876c` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 7.2 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M2/m2
- **extra**: `case/45.isobutane_m6_d155/README.md`, `notes/reviews/2026-10-07-throat-mono-result-interpretation-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

単調拘束の実装と形状ゲート S1〜S8 は支持する。NS の登録窓の準定常判定も再現した。
ただし、生産準備経路の確認と実務判定器の前提検査を補うまで、`accepted` への移動は保留する。

1. **Major — 生産入口の A/B が実際の準備処理を通していない。**  
   [throat_mono_entry_ab.py:26](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/throat_mono_entry_ab.py:26) は、共通の設計結果から `integral_delta_r` と `PhysicalNozzleWall` を直接呼ぶ比較である。`prepare_ns`・メッシュ生成・変換後座標と接続は比較していない。これは較正値の解決を確認する証拠にはなるが、生産入口全体の再現確認には足りない。

   実際、[prep_c2pin.py:11](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/prep_c2pin.py:11) は `cfl_main=5.0`・12000 step を固定しており、検証時は [run_mono_ns_chain.sh:26](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/run_mono_ns_chain.sh:26) で CFL 1・60000 step に書き換えている。YAML の較正値追加だけでは、この検証レシピまで再現しない。

   **対案:** 生産で使うコマンドと引数を確定し、実入口を通す 0-step 準備 A/B を行う。実効 config、設計 spline、δ_r、物理壁、変換後メッシュの座標・接続を保存して比較する。既存の CFD 結果は保持する。

2. **Major — E′ の判定器が証拠欠損でも「採用」を返す。**  
   [throat_mono_practical_eval.py:47](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/throat_mono_practical_eval.py:47) は CSV の枚数・有限性だけを確認し、残差判定、実壁、IC の証拠を読まずに採否を決める。メモリ上の負例検査で、**CSV だけ存在し、それらの証拠が一切ない入力でも「単調壁を候補形状として採用」**を返した。

   [plan §6 E4:220](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-throat-monotone-r2.md:220) の「DIVERGED・欠損は保留」は、E′ で置き換えると明記されていない。既存の前提検査テストが通っても、最終的な採用判定には接続されていない。

   **対案:** E′ の集計前に、同一判定区間の残差、A/B の壁設定・実壁証拠、IC 記録を照合する。ユーザ承認済みの IC 関門上書きと旧形式記録の例外は明示的に扱い、欠損・発散・壁取り違えは保留にする。今回の結果が発散していたという指摘ではない。

3. **Minor — IC 比較の差の向きが逆。**  
   [計算コード:60](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/throat_mono_practical_eval.py:60) と保存 JSON の `ic_beta_minus_alpha` は **β−α＝番号写像−最近傍**で、オーバーシュート η0.1 は **+0.0014334544 %pt**。一方、[plan:244](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-throat-monotone-r2.md:244) は α−β と記載し、[報告生成コード:119](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/throat_mono_summary_pptx.py:119) も「最近傍−番号写像」に同じ正値を表示する。

   **対案:** 表記を β−α に統一するか、α−β を使う箇所では値を負にする。IC 依存を除外できないという結論は変わらない。

4. **Minor — 残作業 #9 は壁形状だけの切り分けにならない。**  
   [plan §5.1 #9:162](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-throat-monotone-r2.md:162) は凍結元を `run_0062` から `run_0143` に替える比較を「壁形状だけが変わる」としている。しかし、[run_0062 の記録:121](/home/sano/work/forge/case/45.isobutane_m6_d155/run_0062_euler_wallfit_fit_r1_ext6k/prepare_info.json:121) は **1100×65**、[run_0143 用 YAML:63](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/problem_d155_euler_pin_G1_recal_mono.yaml:63) は **2000×97**で、出口較正も異なる。この比較から壁形状の寄与を単独で帰属できない。

   **対案:** 同条件の腕 A `run_0140` と腕 B `run_0143` から同じ方法で初期線を抽出する比較に直す。将来課題として移す場合も、この交絡を訂正して移管先を残す。

確認できた裏付けは以下のとおり。

- 単調拘束・メッシュ設定・物理壁・判定関数・排除厚さの5テストは合格。形状ゲートもファイルを書かずに再計算し、**`VERDICT: PASS`**。
- `case/45.isobutane_m6_d155/run_0149_ns_mono_final_ext/`：連結窓 60000〜80000 の4量は **`OVERALL: ALL STEADY`**。
- `case/45.isobutane_m6_d155/run_0148_ns_mono_final_cond/`：窓 14000〜18000 の4量は **`OVERALL: ALL STEADY`**。
- 残差の保存判定は **`NOT CONVERGED (stalled/plateau)`**、品質記録は **`VERDICT: PASS (AR<=5000, skew<=0.90)`**。収束したとは判断しない。
- Euler の保存集計から6量の計算を照合し、最大 `(D＋2SE)/Δq = 0.417795` を確認。ただし、Euler 原系列と新規 NS の残差 CSV・保存場はローカルになく、原系列からの再集計と NaN の再検査はできていない。[run 索引](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:21)。

**推奨は、上記1→2→3→4の順に修正し、限定した生産採用として閉じること。** その後に生産 YAML・`methods/design/overview.md`・報告を同期し、残作業の完了／移管を §5.1、レビュー採否を §6.1 に記録してから `accepted` へ移す。ファイル変更・forge 起動は行っておらず、本レビューの提案は **plan 未反映**。

指摘数: Critical 0 / Major 2 / Minor 2
