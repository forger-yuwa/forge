**未成立 (incomplete)** — incomplete: 証拠不足: `run_0002_ns` の全残差履歴、`stage_manifest.json`・起動時設定、`CONVERGENCE_VERDICT.txt`、`stop_rule.json` 全体（帯・判定内訳・δ系列を含む）、終了状態の記録。 / AWS A で実際に使用した `solverConfig.yaml`・`bcondConfig.yaml`・`prepare_info.json` とバイナリの精度・鮮度の記録。 / 保持済み保存場と対応メッシュ、固定 Euler 参照、報告対象量の `check_quasisteady.py` の VERDICT、局所壁解像の判定。

# codex 諮問 (diagnose): m55-ns-unsettled

- **brief**: [`notes/reviews/briefs/2026-10-15-m55-ns-unsettled.md`](../../notes/reviews/briefs/2026-10-15-m55-ns-unsettled.md)
- **plan**: [`plans/active/tooling-nozzle-card-local-mesh.md`](../../plans/active/tooling-nozzle-card-local-mesh.md)
- **date**: 2026-10-15
- **commit**: `f31c017c` (feature/card-local-mesh)
- **分類**: `incomplete`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 2.2 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-m55-ns-unsettled-diagnose.codex.log` (codex: incomplete)
- **結論**: 長時間の継続投入を保留し、AWS A の本段記録を回収して、既存保存場による θ_r の測定帯の A/B 再評価を先に行う。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 対象 | 判断 | 重大度・根拠・対案 |
|---|---|---|
| 同一設定で延長する | **要再検証** | **Major**：緩慢な過渡とは整合するが、入口直管への原因の局在は未確認。関連 plan も「温度境界層と連動、原因の一意な特定ではない」としている（`plans/active/tooling-nozzle-cooled-wall-ns-settling.md:44`）。まず本段の実効設定・全残差・停止規則の記録を確認する。 |
| 出力間隔を 5,000 step にする | **却下** | **Major**：窓 20,000 step・最低 9 点なので間隔上限は 2,500。`resolve_settings` を読み取りのみで実行し、5,000 は `use=False`、2,000 は `use=True` を確認した（`design/forge_design/evaluate/stop_rule.py:242`）。対案は **2,000 step、`evaluate.stop_rule: true` を明示**。既定の暗黙指定では、使用不能時に停止規則なしへ移る（`runner_axismach.py:2171`）。 |
| 追加 12 万、上限 15 万 step | **却下** | **Major**：提示表の 150k→180k から τ≈62,139 step。ゼロへ単一指数で減衰すると仮定しても、182k から瞬時 r の閾値まで約117k、窓平均の遅れ約10k、保持40kで、**質量条件だけで追加約168k**。θ_r の D も同じ τ で減るという未検証の仮定では、保持込み約169k〜208k。対案は予測の後続窓での検証であり、この外挿を新しい上限の根拠にはしない。保持は実装上必須（`stop_rule.py:613`）。 |
| 継続の判定区間・`--from-floor` | **条件付き採用／却下** | **Major**：実効設定が同じなら、元の本段から継続への残差連結は妥当。ただし `--segment` は指定 run の manifest 内しか連結しない（`solver_density_cuda/tools/check_convergence.py:219`）。183k の場へ戻るなら、それ以後の旧履歴を重ねず、区間と step 対応を記録する。**未収束の元 run を `--from-floor` に使う案は却下**。参照自身の通常判定 PASS が必要（同`:359`）。 |
| 現在の M・θ_r を δ* 補正へ渡す | **却下** | **Major**：提示された停止規則は未到達で、Euler も plan 上で参照として未適格（`plans/active/tooling-nozzle-card-local-mesh.md:71`）。対案は停止規則、全保存量の収束、派生量の準定常、δ_E・δ_target の2窓、局所壁解像を確認してから渡すこと（`design/forge_design/feedback/deltastar_loop.py:905`）。M の +0.4 %だけでは適否を判定できない。 |

結論: **長時間の継続投入を保留し、AWS A の本段記録を回収して、既存保存場による θ_r の測定帯の A/B 再評価を先に行う。**

第 1 仮説: 実際の流れ場に緩慢な緩和が残っているが、「入口直管の冷却境界層が真因」とまでは特定できない。 確度: 中（提示表を条件付きで信用した場合）
  根拠: ブリーフの `case/71.m55_combustion_wt/run_0002_ns` では、150k→180k に r が 0.00551→0.00340、θ_r が3断面とも増加。複数量の変化は残存過渡と整合する。ただし原記録は未確認で、τ は物理時間ではなく反復 step の尺度。
  反証条件: 同じ保存場を適切な固定帯で再評価すると θ_r の変化が停止閾値以下になる場合、θ_r の増加を境界層の未発達の根拠とした部分を撤回する。入口以外に変化が集中すれば、入口直管への局在仮説を撤回する。

第 2・第 3 仮説:
- 第 2：初期の固定帯が後期境界層の縁を覆わず、θ_r の変化を歪めている可能性。未確認。後期の被覆は plan 自身が未確認としている（`tooling-nozzle-cooled-wall-ns-settling.md:122`）。
- 第 3：BC・実効設定・実行バイナリの不整合が緩和を妨げている可能性。未確認。問題 YAML と commit 名だけでは除外できない。

判別 A/B: **追加 forge は 0 step**。保持済みの100k〜183kの同じ保存場列・同じ Euler スナップショット・同じ断面を使い、変えるのは積分帯 y_b だけ。A は元の `stop_rule.json` の固定帯、B は後期場から抽出して全比較時点に固定した帯とする。既存の `theta_r_sections`（`stop_rule.py:279`）で再評価する。両腕とも同じ長い診断窓で Dθ を20k step当たりに換算する。
  → **両腕とも Dθ > 0.05 %なら**「帯選択だけで未到達を説明できる」という第2仮説を棄却する。**Aだけ超過しBが閾値以下なら**測定帯への依存が判別され、θ_r を根拠とした延長予測を撤回する。これは測定感度の診断であり、間引いた系列による正式な停止判定ではない。

やらない方がよいこと: 5,000 step 間隔で延長する、r の瞬時値だけで終了時期を決める、閾値を緩める、未収束 run を参照床にする、M が目標に近いことを理由に δ* 補正へ進むこと。

呼び出し側の前提への異議:
- **Major：継続では「窓だけ」がリセットされるわけではない。** 現実装は起動ごとに `band=None` から始め、最初の出力で帯を再抽出する（`stop_rule.py:848`, `:891`）。残差を連結できても、θ_r は測定帯の同一性を別途確認しないと連結できない。対案は帯を明示して別系列として判定すること。
- 「NaN・Inf は全出力で0」は提示された監査報告であり、全残差履歴の有限性や収束の代わりにはならない。
- θ_r は実装上 r_t 単位（`stop_rule.py:333`）。値そのものの妥当性には、測定帯・縁の状態・断面分布が必要。Euler の軸上 M と NS の出口コア平均 M も異なる指標である。

不足情報: 延長時の容量計画には、空き容量、実測の1枚サイズ、必要な壁出力・後処理成果物の容量も必要。

証拠不足: `run_0002_ns` の全残差履歴、`stage_manifest.json`・起動時設定、`CONVERGENCE_VERDICT.txt`、`stop_rule.json` 全体（帯・判定内訳・δ系列を含む）、終了状態の記録。
証拠不足: AWS A で実際に使用した `solverConfig.yaml`・`bcondConfig.yaml`・`prepare_info.json` とバイナリの精度・鮮度の記録。
証拠不足: 保持済み保存場と対応メッシュ、固定 Euler 参照、報告対象量の `check_quasisteady.py` の VERDICT、局所壁解像の判定。
