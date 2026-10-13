# codex レビュー: tooling-nozzle-deltastar-band-root-zone (result)

- **plan**: [`plans/active/tooling-nozzle-deltastar-band-root-zone.md`](../../plans/active/tooling-nozzle-deltastar-band-root-zone.md)
- **stage**: `result` (diff base `2dda16c5`)
- **date**: 2026-10-13
- **commit**: `84503ceb` (feature/nozzle-deliverables)
- **codex**: effort `high`, 6.3 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M0/m3
- **extra**: `case/70.m5_air_wt_demo/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
コード修正の受入れを推奨します。追加の Critical／Major 不具合は確認できませんでした。  
M5 デモは未達です。移管前に、仕様・完了条件・索引の記述を以下の順で揃えてください。

1. **Minor — §4.2 の端点許容が、実装と採用済みレビューに反しています。**

   **根拠:** [plan:87](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-deltastar-band-root-zone.md:87) は片側の不足だけを許し、超過を例外としています。一方、[deltastar_loop.py:249](/home/sano/work/forge-deliv/design/forge_design/feedback/deltastar_loop.py:249) は両側の絶対差を判定します。これは §6.1 で採用した方針どおりです。

   実データでは入口が設計端より外側に **1.039×10⁻⁷ r_t**、出口が内側に **1.542×10⁻⁶ r_t** ずれています。再確認でも、両側とも許容の 0.5 倍は受理、2 倍は拒否されました。

   **対案:** §4.2 を「入口・出口とも |Δ| ≤ ε_f32」に訂正してください。実装の変更は不要です。

2. **Minor — コード修正に絞った完了条件と、未達の検証条件の関係が曖昧です。**

   **根拠:** [plan:196](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-deltastar-band-root-zone.md:196) は「V1〜V4 が合格」としていますが、[V3e:160](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-deltastar-band-root-zone.md:160) は未達です。§4.4 に逸脱の承認はありますが、完了条件から明示的に外したのは V5 だけです。また、§5.1 #5 は Euler 参照の受入れ作業を未整理のまま残しています。

   再実行した収束判定も、`case/70.m5_air_wt_demo/run_0012_euler_c2b/`・`run_0012b_euler_c2b_cont12k/`・`run_0013_euler_ns2/` は **NOT CONVERGED** でした。

   **対案:** §8 に受入れ対象を列挙し、V3e・Euler 参照の厳密受入れ・V5 を「未達、今回のコード受入れ対象外」と明記してください。#7・#10・#11 は残し、後続の active plan への引継ぎ先を示してください。追加 CFD を今回の移管条件にする必要はありません。

3. **Minor — 索引と影響範囲に古い記述が残っています。**

   **根拠:** [plans/README.md:62](/home/sano/work/forge-deliv/plans/README.md:62) は「result 段は未実施 (#12)」ですが、現在の該当項目は #8 で、レビュー記録も存在します。[plan:191](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-deltastar-band-root-zone.md:191) の影響範囲には、変更した `wt_card.py`、設計壁ゲート、MOC 出口切断が抜けています。

   **対案:** 索引・§7・残作業の状態を現在の実装に合わせ、移動時に `active/` への参照リンクも更新してください。#9 の `de35cb0b` の連絡も、実施済みか引継ぎ済みかを記録してください。

検証では、次を確認しました。

- 指定 diff に CUDA ソルバの変更はなく、境界半割面・周期 seam・軸の処理を直接変更していません。
- V1 は旧実装の δ_r = r_w を再現し、新実装の真値に対する相対誤差は **1.29×10⁻⁷**。正常な合成断面 13 個はビット同一でした。既存 `run_deltastar_tests.py` と単調性ゲート 42 件も合格しました。
- case/45 の 4 対は、保存された新旧 `delta_r_equiv.csv`・`delta_r_next.csv` が実際にバイト同一でした。
- `_ns2_extract/` から補正表を再生成すると保存表と全値一致し、カードの SHA-256 照合も通りました。
- 出口切断 A/B の記録は、c2b **+0.845760 mm**、case/45 生産 **+5.353145 mm** を裏付けています。この既定挙動の変更は手順書にも記載済みです。

M5 の根拠は次のとおりです。出口 M は過渡値です。

| `case/70.m5_air_wt_demo/` 配下 | 収束 VERDICT | 出口 M 専用の準定常 VERDICT |
|---|---|---|
| `run_0014_ns_ns2/` | **PASS (converged)**、main 単独区間 | **DRIFTING**、最終 M = 5.002198 |
| `run_0014b_ns_ns2_cont60k/` | **NOT CONVERGED** | **DRIFTING** |
| `run_0014c_ns_ns2_cont120k/` | **NOT CONVERGED** | **DRIFTING**、今回の再判定でも再現 |
| `run_0014d_ns_ns2_cfl5/` | **NOT CONVERGED** | **DRIFTING** |

根拠は各 run の `CONVERGENCE_VERDICT.txt`、`quantities_series.csv`、[専用判定記録](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/_v5_judgment/v5_judgment.txt)、[run 一覧](/home/sano/work/forge-deliv/case/70.m5_air_wt_demo/README.md:10)です。したがって、60k 時点の値を定常性能の改善率として扱うことはできません。今回の「V5 未達」への訂正は妥当です。

**推奨は、上記3点を直したうえで、コード修正の完了として accepted に移すことです。** 全43本の回帰試験は書き込みを伴うため再実行しておらず、「41本合格・2本は main でも失敗」は本レビューで独立再現した範囲には含めません。

ファイルは変更していません。指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 0 / Minor 3
