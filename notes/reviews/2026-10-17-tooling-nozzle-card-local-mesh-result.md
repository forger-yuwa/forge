**未成立 (incomplete)** — incomplete: 証拠不足: V8 の AWS B `~/m55_p5/cmp_float60.log` と、４ run の step 60,000 に対応する比較入力または抽出記録。手元の写しでは、特に δ_E の差 0.015〜0.039 mm を独立に照合できませ …

# codex レビュー: tooling-nozzle-card-local-mesh (result)

- **plan**: [`plans/accepted/tooling-nozzle-card-local-mesh.md`](../../plans/accepted/tooling-nozzle-card-local-mesh.md)
- **stage**: `result` (diff base `48cab832`)
- **date**: 2026-10-17
- **commit**: `fbbb9c99` (feature/card-local-mesh)
- **分類**: `incomplete`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 7.6 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-17-tooling-nozzle-card-local-mesh-result.codex.log` (codex: incomplete)
- **判定**: **NO-GO**, 指摘 C0/M2/m2
- **focus**: §8 の完了の範囲 (2026-10-17 ユーザ決定: M5.5 系は hand_mesh.py、生成器への組み込みと係数の試行錯誤のループは tooling-design-orchestrator §5.1 #10 へ) で、§6 の V2〜V9 の判定が事前登録どおりか、case/71 run_0014 の accepted (公差 0.1 %、ユーザ決定) と手順書の試走 A の記録が主張を支えるかを見る。証拠: case/71.m55_combustion_wt/_runs/ (run_0014_chain_v2/accepted_2026-10-16・qs_trial_a_evidence・lcscan_e6・lcscan_forge7d) と手元の run の写し

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: NO-GO
`run_0014` の合格値は再現できましたが、判定ツールに誤って `accepted` に進む経路が２件あります。
V7 の解釈と現在仕様の文書にも修正が必要です。V8 の一次証拠が不足しているため、レビュー全体は未成立です。

1. **Major — 未監視の線を、準定常確認なしで合格にできる。**  
   [wt_qualify.py:535](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:535) は、評価線と停止判定の監視線が異なっても警告だけで処理を続けます。`run_0014_chain_v2` に `cli_tols={"line_eta": 0.06}` を渡すと、`matches_settle_line=false` のまま `accepted` を再現しました。最大偏差 0.061881 %・片振幅 0.061542 %という単一時刻の値だけでは、その線の定常性を保証できません。  
   **対案:** 監視対象との不一致・監視量の欠落・停止 step との不一致は `incomplete` にする。別の線や時刻を認める場合は、その対象について `check_quasisteady.py` の合格を要求する。

2. **Major — 指定 step・全壁の出力が揃っていなくても壁解像が PASS になる。**  
   新しい呼出側 [wt_qualify.py:309](/home/sano/work/forge-cardmesh/design/forge_design/feedback/wt_qualify.py:309) は、ツールの `VERDICT` を受理するだけです。[check_wall_resolution.py:335](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/check_wall_resolution.py:335) は指定 step が無い場合に最新出力へ代替し、存在しない指定壁も検出しません。実際に `run_0014_chain_v2 --step 46170 --groups wall,missing_wall` は、**step 46171 の wall だけで `VERDICT: PASS`、終了コード 0**でした。  
   **対案:** 全対象壁について指定 step のダンプが存在することを先に検査し、欠落時は判定不能で停止する。別 step への代替は禁止し、呼出側でも評価済みの壁・step を照合する。

3. **Minor — V7 に、事前登録で支持されない因果断定が残っている。**  
   [plan:157](/home/sano/work/forge-cardmesh/plans/accepted/tooling-nozzle-card-local-mesh.md:157) の「平均の過大は途中の厚すぎが主因だった」は、同文書の修正方針と矛盾します。片振幅は 0.0788 → 0.0761 %で半減せず、保存された収束判定も不合格です。準定常判定の `ALL STEADY` は、V7 の主因説の支持を意味しません。  
   **対案:** 当時の事前登録に沿って「判別不能」と明記し、観測事実を「平均の偏りは改善、傾きの原因は未確定」に統一する。

4. **Minor — 現在仕様と完了範囲の文書が同期していない。**  
   [methods/design/overview.md:232](/home/sano/work/forge-cardmesh/methods/design/overview.md:232) は今回の終端 C⁺ の補間を説明せず、[同:622](/home/sano/work/forge-cardmesh/methods/design/overview.md:622) にも３係数補正の仕様がありません。また [plan:250](/home/sano/work/forge-cardmesh/plans/accepted/tooling-nozzle-card-local-mesh.md:250) は「既存カードは新しい既定で格子が変わる」としており、生成器への組込みを延期した §8 と矛盾します。  
   **対案:** 実装済みの補間・３係数・qualification を現在仕様へ反映し、§7 と計画索引を今回の完了範囲に合わせる。後続計画 #10・#12 への移管自体は確認できています。

確認できた証拠は次のとおりです。`check_settled.py --no-write` を再実行し、設定された量・区間について以下を確認しました。

| run（リポジトリ相対パス） | 判定区間／停止 step | VERDICT |
|---|---|---|
| `case/71.m55_combustion_wt/run_0014_chain_v2/` | 36,000–46,000／46,171 | `SETTLED` |
| `case/90.qs_trial_a/run_0003_ns_p0/` | 32,000–42,000／42,060 | `SETTLED` |
| `case/90.qs_trial_a/run_0005_ns_p1/` | 20,000–30,000／30,088 | `SETTLED` |

３本とも最終 `res_<停止step>.h5` の全 `VALUE/*` と全残差履歴に NaN/Inf はありません。`run_0014` の出口コア M = 5.502841913、線の最大偏差 0.062305 %・片振幅 0.057285 %、公差 0.1 %での `accepted` を再現しました。壁解像も正しい停止 step では `PASS`、最大 y₁⁺ = 0.999、超過面積 0.0 %です。

V2・V3 の保存記録は品質 `PASS`、第一層の要求超過なし。V9 の保存 JSON を再集計すると、2400／4800 の出口位置差は最大 1.99383e−6 rₜ、直線当てはめ残差は 1.34e−6／3.67e−7 rₜで、記載値を裏付けます。run 索引は [case/71 README](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:23)、試走の証拠は [試走 A README](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_runs/qs_trial_a_evidence/README.md)です。

**推奨:** `active` のまま２件の誤受理経路と文書を修正し、V8 の一次記録を揃えてから結果レビューを完了させてください。ファイル変更はしていません。

証拠不足: V8 の AWS B `~/m55_p5/cmp_float60.log` と、４ run の step 60,000 に対応する比較入力または抽出記録。手元の写しでは、特に δ_E の差 0.015〜0.039 mm を独立に照合できませんでした。

指摘数: Critical 0 / Major 2 / Minor 2
