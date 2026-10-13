# codex レビュー: tooling-nozzle-deltastar-band-root-zone (plan)

- **plan**: [`plans/active/tooling-nozzle-deltastar-band-root-zone.md`](../../plans/active/tooling-nozzle-deltastar-band-root-zone.md)
- **stage**: `plan`
- **date**: 2026-10-13
- **commit**: `c13be1fa` (feature/nozzle-deliverables)
- **codex**: effort `high`, 4.8 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2
- **extra**: `case/70.m5_air_wt_demo/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
探索を積分域に限る修正と、実寸を保ってカードから再設計する方針は妥当です。  
ただし、異常値の拒否、端点の丸め、準定常の受入れ条件に実装前の修正が必要です。

目的は未解決の実装不具合を正しく捉えています。`plans/README.md` と既存の accepted plan を確認しましたが、今回の探索不具合は解決されていません。変更を設計・後処理層に置く構成も適切です。CUDA の流束・block-DPLUR・周期処理を変更する必要はなく、構造 node 格子を対象とした case/45・70 の検証は、現行の node 限定方針と整合します。

以下はコード、既存成果物、書き込みを伴わない再現試験に基づく指摘です。新しい CFD 計算は実行していません。

1. **Major — `nonmono_zone` を hard に追加するだけでは、不合格点が壁更新に使われます。**

   **根拠:** [deltastar_loop.py:85](/home/sano/work/forge-deliv/design/forge_design/feedback/deltastar_loop.py:85) は hard 点の平滑化重みをゼロにしますが、平滑化結果 `d_ext` を全点で評価して `d_next` に入れます。`held` は記録されるだけです。さらに [同:335](/home/sano/work/forge-deliv/design/forge_design/feedback/deltastar_loop.py:335) の寸法解きは、有限な `delta_E` を読み、`held` を検査しません。

   メモリ上の再現試験では、不合格点の前回値を `0.5`、周囲の抽出値を `0.1`、`omega=1` とすると、出力は **`held=1, delta_target=0.1, delta_E=0.1`** でした。§4.1 の「前回値を保持する」は現実装と異なります。

   **対案:** 今回の NS 後の寸法解きでは、使用域 `[0,x_e]` に hard／held 点があれば停止する契約にしてください。V1 に「断面で NaN」だけでなく、平滑化・CSV・寸法解きまで通して不合格が伝播する試験を追加してください。

2. **Major — 単調性検査を D≥0 の枝に限定すると、壁側の負の参照流束を取り逃がします。**

   **根拠:** [deltastar.py:270](/home/sano/work/forge-deliv/design/forge_design/metrics/deltastar.py:270) の D<0 は `_negative_delta_r` に直行し、[同:126](/home/sano/work/forge-deliv/design/forge_design/metrics/deltastar.py:126) は非正の参照流束を `1e-30` に置換します。

   合成断面 `r_w=1, q_E=1, q_NS=max(−0.1+4y,0), y_b=0.1` を実行すると、入力流束は全点非負でも、帯から壁への外挿で `q_ref,wall≈−0.1` となり、**D=−0.00778864、δ_r=−4.979×10¹³、理由は soft の `negative_deficit` のみ**でした。現計画の修正対象外の枝なので、この異常は残ります。

   **対案:** 参照流束の有限性・正値性と探索域の妥当性を、D の符号分岐より前に検査してください。正の参照に対する正常な負の δ_r は維持します。V1 にこの負例を、V2 に新しい拒否条件を含めてください。

3. **Major — #18 の片側許容は、正常な float32 の丸めを拒否します。**

   **根拠:** §4.2 は出口を超えた表を一律例外にします。しかし、実在する `case/70.m5_air_wt_demo/run_0009_ns_c2_final/` では、

   - [wall_repr.json:14](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/run_0009_ns_c2_final/wall_repr.json:14): 出口 `48.03176227798793`
   - `delta_r_next.csv` 最終行: `48.03176498413086`
   - 差: **出口を 2.70614×10⁻⁶ r_t 超過**

   となっています。この表の終点は、同じ出口を float32 の実寸座標に丸め、既存の抽出と同じ換算をした値で再現できました。丸めは内外両方向に起きます。

   **対案:** 出典と端点対応を照合したうえで、入口・出口とも **絶対差 ≤ ε_f32** を認めて正規化してください。範囲不足や別設計の表は引き続き拒否し、`applicability` は変更しません。V4 に内側・外側への丸めと、負座標の入口を追加してください。

4. **Major — 再設計の入力となる δ_E 分布自体の準定常確認が抜けています。**

   **根拠:** §5.1 #5 の延長 A/B は Euler 参照の影響を検査しますが、NS は `run_0010` の最終場に固定しています。V3 にも NS の時系列による δ_E の受入れ条件がありません。移管元の [deliverables plan:263](/home/sano/work/forge-deliv/plans/accepted/tooling-nozzle-wt-deliverables.md:263) は「使う δ_r」の準定常確認も明示しています。

   再実行したツールでは `run_0010_ns_c2b_final` は **`PASS (converged)`**、出口コア M・流量系列は **`OVERALL: ALL STEADY`** でした。一方、同 run の [MANIFEST.json:865](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/run_0010_ns_c2b_final/deliverables/MANIFEST.json:865) では壁熱量が **`TRANSIENT-UNSETTLED`** です。これは δ_E の非定常を証明しませんが、出口 M と流量だけで冷却境界層全体の整定を代用できないことを示します。

   **対案:** 受入れ済み Euler を固定し、NS の複数の末尾スナップショットから同じ抽出・平滑化で δ_E(x) を作ってください。`check_quasisteady` の VERDICT に加え、使用域全体の最大変化を、例えば既登録の感度閾値と同じ **1e-4 r_t 以下**で判定してください。V5b の固定点差を解釈する際にも同じ確認が必要です。

5. **Major — 既定の `STEADY` は、出口 M の ±0.05% を保証するには緩すぎます。**

   **根拠:** [qs_series.py:43](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/qs_series.py:43) は `--drift`・`--osc` を指定していません。[check_quasisteady.py:543](/home/sano/work/forge-deliv/solver_density_cuda/tools/check_quasisteady.py:543) の既定値は、それぞれ **5%・10%** です。

   同ツールの判定関数に合成系列 `Mₙ=5.01−0.032768·0.8ⁿ`、n=0…5 を与えると、**`STEADY`、漸近値 5.01** が返りました。最終値 **4.99926258** は V5 の範囲内ですが、漸近値は範囲外です。つまり現条件では過渡の通過点を合格にできます。

   **対案:** M 専用の準定常許容を明記し、末尾変化と推定残変化を **2.5e-4 以下**など、合否半幅 0.0025 より十分小さくしてください。単調系列では漸近値も合否範囲内であることを要求します。流量・δ_E は別の許容で評価してください。

6. **Minor — V5 の「同じ IC」は事実と異なります。**

   **根拠:** §4.3 の新 NS は `run_0010` を IC にしますが、[run_0010 の prepare_info.json:636](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/run_0010_ns_c2b_final/prepare_info.json:636) の IC は `run_0009_ns_c2_final` です。したがって既存 `run_0010` との比較は、壁と r_t だけを変えた同一起点 A/B ではありません。

   **対案:** V5 は各設計の受入れ済み定常解どうしの比較と記述してください。上記の準定常ゲートを満たせば、同一 IC を合格の前提にする必要はありません。

7. **Minor — 実装手順と仕様文書の更新時期を揃える必要があります。**

   **根拠:** [plan:111](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-deltastar-band-root-zone.md:111) の実装ステップは §4.3 の「1〜4」で終わり、最終 NS・V5 を含みません。§7 の影響範囲にも、新たに変更する `wt_card` と CLI が抜けています。また、[AGENTS.md:166](/home/sano/work/forge-deliv/AGENTS.md:166) は設計変更の実装前に `methods` 更新を求めていますが、計画では完了条件にしかありません。

   **対案:** `methods` 更新を実装前に置き、カード→補正表→`prepare_ns`→寸法照合の一貫性試験、最終 NS、result レビューまでを §5・§7 に揃えてください。

推奨は、**探索域限定＋実寸固定のカード再設計という現在の経路を維持し、指摘1〜5をこの順で計画に反映してから実装すること**です。特に hard 判定の伝播と符号分岐前の検査を先に固め、その後に端点処理、入力 δ_E と最終 M の受入れ条件を確定してください。帯選択法の変更や Mach の係数調整を追加する必要は、現時点ではありません。

読み取り専用の依頼に従い、ファイルは変更していません。上記提案は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 2
