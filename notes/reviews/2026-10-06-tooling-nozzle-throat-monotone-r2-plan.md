# codex レビュー: tooling-nozzle-throat-monotone-r2 (plan)

- **plan**: [`plans/active/tooling-nozzle-throat-monotone-r2.md`](../../plans/active/tooling-nozzle-throat-monotone-r2.md)
- **stage**: `plan`
- **date**: 2026-10-06
- **commit**: `4b23209a` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.8 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2
- **extra**: `case/45.isobutane_m6_d155/README.md`, `notes/reviews/2026-10-06-throat-monotone-r2-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

単調拘束の方針は妥当です。再計算でも r″ 最大値は 0.510164 → 0.500000、積分方法への感度は D₃=2.40%、D₄=5.13% でした。  
ただし、評価器・物理壁ゲート・生産化条件に穴があり、現状の計画では「検証合格」を保証できません。

1. **Major — 生産化条件が統一されておらず、凝縮の合否条件も欠けています。**

   **根拠:** [plan:79](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-throat-monotone-r2.md:79) は「Euler 合格なら生産キー化」、§5.1 は NS 合格後、[§8:176](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-throat-monotone-r2.md:176) は S・E の保留を諮問しただけでも完了できる記述です。[N:154](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-throat-monotone-r2.md:154) も、準定常判定の実行は要求しますが、必要な VERDICT と凝縮の合格条件を明記していません。

   実際、参照 `case/45.isobutane_m6_d155/run_0117_ns_recal_final_ext/quantities_series.csv` を再判定すると、40000–60000 step の `wave01` は **DRIFTING**、出口 M・オーバーシュート・δ_E/δ_C は **STEADY**。保存済み `CONVERGENCE_VERDICT.txt` は **NOT CONVERGED** です。「run_0117 と同じ」だけでは、数値合格と過去の個別免除を区別できません。

   **対案:** 全節を **S → E → dry NS → 凝縮評価 → result レビュー → 生産化** に統一してください。dry の対象量・判定窓・必要 VERDICT、凝縮4量の STEADY と残差 RISING なし、未達時の扱いを明記する。旧壁へのユーザー判断は新壁に継承せず、「機能実装完了」「生産採用」「性能認定」を分けます。

2. **Major — 評価器の改修範囲に、E3・E4 の判定実装が抜けています。**

   **根拠:** [plan:101](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-throat-monotone-r2.md:101) の一般化は主に座標・run 名・出口標本です。しかし [eval_wallfit_euler.py:170](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/eval_wallfit_euler.py:170) は絶対差による対称判定で、計算した `noninferior` を [総合判定:188](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/eval_wallfit_euler.py:188) に使いません。例えば Δ=0.003、U=0.0003、B−A=−0.010 なら、E3 は採用ですが既存総合判定は「差あり」です。

   また、準定常ツール呼び出しは `--tail` 未指定。既定は40%なので、18枚なら末尾8枚となり、E4 の「末尾5枚」と一致しません。

   **対案:** E2 の判定対象、η0 の除外規則、片側の非劣化判定、保留条件を独立した判定関数として実装対象に追加してください。準定常と T は同じ末尾5枚を使い、改善・悪化・境界・未定常・欠損データを入力した判定テストを設けます。

3. **Major — 出口 M の評価誤差が、出口規格化オーバーシュートへ伝播していません。**

   **根拠:** [plan:143](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-throat-monotone-r2.md:143) は `E_exit` を出口 M のみに適用します。一方、[eval_wallfit_euler.py:74](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/eval_wallfit_euler.py:74) の規格化オーバーシュートは `100·(M/M_exit−1)` です。

   M≈M_exit≈6 では、出口評価が 1e−4 動くと、この量は約 **0.00167 %pt** 動きます。許容幅 0.003 %pt の約56%であり、出口 M だけに不確かさを入れる扱いは不十分です。

   **対案:** E2 の出口標本・定義の各変種について、`|M_exit−6|` と規格化オーバーシュートを両腕とも再計算し、**B−A の感度を各量の U に入れる**ことを明記してください。

4. **Major — `wall_design.csv` は、単調壁がメッシュに使われた証拠になりません。**

   **根拠:** [plan:134](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-throat-monotone-r2.md:134) は同 CSV との照合を要求します。しかし [runner_axismach.py:568](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:568) が保存するのは **当てはめ前の `wall_inv`** です。今回の変更では MOC 点群は変わらず、この CSV も変わりません。

   run_0114 の保存点群から再計算したところ、単調壁の MOC 点上最大偏差は約 **1.055 µm**、現行壁との最大差は約 **0.501 µm**。参照対象の取り違えが変更量を上回ります。

   **対案:** 当てはめ後の spline のノット・係数・設定を保存し、変換後の実壁節点位置で解析壁を評価して照合してください。A/B の実座標差と float32 丸め後の差も記録する。`prepare_ns` の情報には現在 `wall_fit` がないため、[runner_axismach.py:1093](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1093) への明示的な追加も実装項目にします。

5. **Major — S6 の例外条件は、物理壁の山を実質的に制限していません。**

   **根拠:** [plan:119](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-throat-monotone-r2.md:119) は、増加を δ_r″ の寄与で説明できれば合格です。[PhysicalNozzleWall.r:497](/home/sano/work/forge-integ-1005/design/forge_design/geometry/wall_axismach.py:497) より、対象区間では常に  
   `r_phys‴ = r_design‴ + δ_r‴`。  
   設計壁が単調なら、物理壁に残る増加は必然的に δ_r 側で説明できます。極値の個数も、山の高さを制限しません。

   **対案:** 生産経路で再計算する方針は維持し、物理壁の最大増加量 `maxₓ<ᵧ[r_phys″(y)−r_phys″(x)]` に許容値を設定してください。δ_r による説明は診断情報とし、合格条件から外す。残存増加の許容値を決められなければ S6 は保留です。

6. **Minor — 正規化後の許容誤差を、導関数の単調性保証と同一視しています。**

   **根拠:** [plan:68](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-throat-monotone-r2.md:68)・S1 は行正規化後の 1e−10 を使います。今回のノットから再構成した3階微分行列の最大行スケールは **5.3752e7**。該当行では、正規化後 1e−10 は未正規化で約 **0.00538** に相当します。r″≤0.5 の検査だけでは、その上限より下での増加を検出できません。また63行目には旧停止条件も残っています。

   **対案:** 解法用の正規化許容差と、形状用の導関数許容差を分けてください。未正規化係数・区間多項式の極値でも単調性を検査し、停止条件の記述を一本化します。

7. **Minor — §9 の「厳密評価値」に別腕の値と格子評価が混在しています。**

   **根拠:** [throat_mono_integration_ab.py:82](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/throat_mono_integration_ab.py:82) 以降の最大値は均等点列上の評価です。採用予定の `grid` 腕を区間端点で厳密評価すると、`max|r⁗|` は **58.83977**。plan の **58.3** は `gauss` 腕の値です。乗数最小値も `grid` は **4.10e−10**、記載の約2e−9は `gauss` 側です。

   **対案:** 採用腕を明示して数値を訂正し、S8 の極値評価を区間多項式で実装してください。今回の大小関係と積分方法の採否は変わりません。

**推奨は、既定経路を保った opt-in の単調拘束を採用し、上記を修正してから実装へ進むことです。** 優先順は **①生産化ゲート → ②③判定器と不確かさ → ④⑤実形状の証拠と物理壁ゲート → ⑥⑦数値仕様の整理**です。

目的は未解決の形状要求として適切で、既存 accepted plan の重複実装には当たりません。変更は設計側の spline に閉じるため、現時点で CUDA・block-DPLUR・周期境界の変更は不要です。case/45 の node 検証を中心にする範囲も現行検証方針と整合します。

ファイル変更なし、**plan 未反映**。参照 run は主作業ツリー `/home/sano/work/forge/` で確認しました。準定常は CSV から再判定、収束は保存済み VERDICT を照合しています。

指摘数: Critical 0 / Major 5 / Minor 2
