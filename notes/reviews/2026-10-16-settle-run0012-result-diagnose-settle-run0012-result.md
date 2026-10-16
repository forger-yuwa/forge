# codex 諮問 (diagnose): settle-run0012-result

- **brief**: [`notes/reviews/briefs/2026-10-16-settle-run0012-result.md`](../../notes/reviews/briefs/2026-10-16-settle-run0012-result.md)
- **plan**: [`plans/active/tooling-settle-convergence-and-monitor.md`](../../plans/active/tooling-settle-convergence-and-monitor.md)
- **date**: 2026-10-16
- **commit**: `90913054` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 2.9 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-16-settle-run0012-result-diagnose-settle-run0012-result.codex.log` (codex: ok)
- **結論**: (i) を選び、停止規則と閾値を固定したまま、本段上限だけ260,000へ延ばす新しい自由継続runを1本、事前登録して実施する。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 4 / Minor 1）

| 重大度 | 判断対象 | 採否・根拠・対案 |
|---|---|---|
| Major | (i) 新しい run で自由継続を長くする | **採用**。`run_0012_ns_free_judge/settle_virtual_stop/VERDICT.txt` は **未成立**。136,000〜200,000 の33格子点には反例がなく、必要な216,000までの8格子点が欠けている。元の結果は変更せず、上限260,000・t_s上限180,000・検証長80,000を次の試験として事前登録する。 |
| Major | 残差の条件は不要／(iii) 保持の規則を直す | **現時点では却下**。同ディレクトリの `grid_verdicts.csv` では、不合格8点はすべて残差で、流量は全80点で合格。t_s′からの検証も超過なしなので、**この軌道では残差による54,000 stepの追加待機の利益は観測されなかった**。ただし、他条件でも不要という証拠ではない。固定したAND条件を維持して検証を完結させる。 |
| Major | 今回の結果でV4へ進む | **却下**。[plan:129](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:129) は結果Aを前提としている。今回の未成立をA相当に扱わず、次の自由継続試験の結果Aを待つ。 |
| Major | #1cをV1(f′)・(g′)合格にする | **要再検証**。診断はhostコピーを広く比較し、host差があればFAILにする実装（[main.cpp:2665](/home/sano/work/forge-cardmesh/solver_density_cuda/main.cpp:2665)、[main.cpp:2763](/home/sano/work/forge-cardmesh/solver_density_cuda/main.cpp:2763)）。一方、契約は「次の数値更新が読む状態」の不変性である。**全配列差0の不合格は残し、契約上の合否は保留**する。対案は、差が出た配列ごとに「読まれない／読む前に必ず上書き／上書き前に読む」をコード位置付きで網羅すること。間接的なH2D転送も含め、最後の分類では値の不変性が必要。`farfield` は実際にhostを読む（[convectiveFlux_d.cu:119](/home/sano/work/forge-cardmesh/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:119)）ので、未試験のまま除外しない。NaNが出ない継続試験だけではこの証明を代替できない。 |
| Minor | t_s′以後の最大比0.071 | **訂正を採用**。`span_deviation_ts_prime.csv` の82,000〜162,000に限った最大比は **0.0665224461**（146,000、壁熱流）。対案は区間を固定してこの値を記載すること。結論は変わらない。 |

結論: **(i) を選び、停止規則と閾値を固定したまま、本段上限だけ260,000へ延ばす新しい自由継続runを1本、事前登録して実施する。**

第 1 仮説: 今回は停止規則の危険性が現れたのではなく、保持のリセットで検証期間が足りなくなったもので、完全な80,000 stepの追跡では結果Aになる。 確度: 中  
  根拠: `case/71.m55_combustion_wt/run_0012_ns_free_judge/settle_virtual_stop/` のCSVから既存の `emulate_v2` を再実行し、t_s＝136,000、t_s′＝82,000を再現した。94,000での超過は0.0228282848 > 0.0211892991で、保持リセットは[stop_rule.py:629](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:629)の仕様どおり。読めた136,000〜200,000では全条件がSETTLED、目的量の累積差／許容の最大は0.055201。判定コード2本のSHA256も保存記録と一致した。  
  反証条件: 次のrunでt_s ≤ 180,000を得た後、t_s〜t_s＋80,000の目的量の累積差が許容を超える、または全条件がNOT_SETTLED／DIVERGEDになる。

第 2 仮説: #1cのhost差は出力用コピーの更新であり、試験した数値経路には影響しない。確度: 中。ただし全読み手の網羅とFP64への適用は未確認で、一般的な無影響とは確定できない。

判別 A/B: **変更する計算設定は本段上限200,000→260,000の1点だけ**。初期場・バイナリ・格子・閾値・W＝20,000・保持40,000・出力間隔2,000を維持し、途中restartを挟まない。判定側のt_s上限は260,000−80,000＝180,000として事前登録する。比較は既存ツールによる同一軌道の「仮想停止時点」と「その後80,000 step」で行い、対照用の別runは不要。  
→ **A**：全区間で累積差／許容 ≤ 1、全条件SETTLEDなら、この条件・検証長での早すぎる停止という仮説を退け、V4へ進む。  
→ **B**：どちらかに反例が出れば「観測期間不足だけ」という第1仮説を退け、停止規則の設計へ戻る。  
→ t_s > 180,000・候補なし・必要データ不足なら再び未成立。合格するまで上限を延ばす運用にはしない。

やらない方がよいこと: 94,000の1点を見て閾値や連続超過回数を調整すること、run_0012の末尾にrestart後の系列を継ぎ足して「同じ自由軌道の追跡」とすること、host差を消すためだけに出力コピーを復元すること。

呼び出し側の前提への異議: 「1点の僅かな超過」は**判定格子上の1点**であり、生の残差の単発異常とは限らない。[settle.py:1049](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/settle.py:1049)では窓内のlog残差から傾きを求めている。また、「残差が独立の情報を足した」と「その情報が早すぎる停止を防いだ」は別で、今回確認できたのは判定の差だけである。

不足情報: #1cの配列別原報告・全読み手と上書き順序の対応表・FP64での適用根拠。これらはV1の合格確定に必要だが、今回の未成立を保持して次の検証を選ぶ判断は妨げない。

証拠: 足りる
