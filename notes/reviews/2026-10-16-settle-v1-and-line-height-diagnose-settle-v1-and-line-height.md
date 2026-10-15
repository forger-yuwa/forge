# codex 諮問 (diagnose): settle-v1-and-line-height

- **brief**: [`notes/reviews/briefs/2026-10-16-settle-v1-and-line-height.md`](../../notes/reviews/briefs/2026-10-16-settle-v1-and-line-height.md)
- **plan**: [`plans/active/tooling-settle-convergence-and-monitor.md`](../../plans/active/tooling-settle-convergence-and-monitor.md)
- **date**: 2026-10-16
- **commit**: `9f9d029b` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 6.5 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-16-settle-v1-and-line-height-diagnose-settle-v1-and-line-height.codex.log` (codex: ok)
- **結論**: V1 の合格確定を保留し、同じ数値状態で随時出力の有無だけを変える、出力処理前後の状態不変 A/B を次に行う。
- **extra**: `plans/active/tooling-nozzle-card-local-mesh.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表です。**V1 の代替検査は補助証拠として採用しますが、「数値を変えない」の確証には不足します。η = 0.03 は近軸の評価線として採用できます。**

| 重大度 | 論点・採否 | 根拠と対案 |
|---|---|---|
| **Major** | 「600 step の差が同じ桁だから無影響」：**却下** | 旧バイナリの再実行差は実在します。ただし、その大きさに隠れた変更は排除できません。追加出力は device→host コピーと host 上の集計です（[monitorOutputs.cpp:214](/home/sano/work/forge-cardmesh/solver_density_cuda/output/monitorOutputs.cpp:214)、[同:515](/home/sano/work/forge-cardmesh/solver_density_cuda/output/monitorOutputs.cpp:515)）。対案は、下記の**出力処理前後で数値状態が変わらない検査**です。 |
| **Major** | V1 (a)・(f)・(g) の変更：**採用。ただし事後改訂と明記** | [plan:111](/home/sano/work/forge-cardmesh/plans/active/tooling-settle-convergence-and-monitor.md:111) は独立した run のビット一致を要求しています。旧条件を「判定不能」として残し、改訂版の条件・理由・適用日を追記してください。既存結果を元の条件の PASS に書き換えることは認めません。 |
| **Major** | η = 0.03 なら軸の影響を避けられる：**要再検証** | `eta_line` は径方向を線形補間します（[nozzle_report.py:138](/home/sano/work/forge-cardmesh/design/forge_design/report/nozzle_report.py:138)）。今の NS 格子では第1層に一致しますが、将来 η₁ > 0.03 なら軸節点が混ざります。**丸め誤差を除いて全対象列で η₁ ≤ 0.03**を確認し、補間節点と軸の重みを記録してください。 |
| **Major** | η = 0.03 を試験部全体の一様性と扱う：**却下** | `run_0010_ns_b0/res_86000.h5` を再抽出すると、η = 0.03 と 0.1 の点ごとの差は最大 **0.042737 percentage points、x = 4.0 m**。測る位置を変えれば拾う分布も変わります。新指標を「η = 0.03 の近軸線の一様性」と明記し、η = 0.1 は比較用に残してください。 |
| **Minor** | 最大偏差 0.1 %・片振幅 0.05 % の継承：**採用** | これは測定値から推定した定数ではなく、ユーザ了承済みの要求公差です（[card plan:154](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:154)）。新しい評価線にも同じ要求を課して構いません。ただし、旧評価と同じ厳しさ・同じ評価範囲を保証するものではありません。 |
| **Minor** | 次の run から変更し、`run_0010` は旧判定を保持：**採用** | [card plan:158](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:158) の扱いは妥当です。η・x 範囲・補間・標本点数・公差を一緒に版管理してください。`settle_presets` は既に 0.03、`nozzle_report.metrics` はまだ 0.1 なので、報告と監視の定義を揃える必要があります。 |

既存データから確認できたことを補足します。

一時 run の親ディレクトリを  
`/tmp/claude-1000/-home-sano-work-forge/e254baf2-8857-4243-99d2-53804f553ad9/scratchpad/runs/`  
として、既存の `check_field_regress.py` を **step 600・factor 2** で実行しました。保存量・化学種保存量・P・T・h0 が対象です。

| 比較 | 対照の反復 run | L2／L∞ のノイズ床比の最大 | VERDICT |
|---|---|---:|---|
| キー無し `run_t1_new` | `run_t1_base_a/b` | 1.48 | `PASS` |
| キーあり `run_t1_newkeys` | `run_t1_base_a/b` | 1.66 | `PASS` |
| 随時出力 `run_g1_snap` | `run_g0_noreq_a/b` | 1.14 | `PASS` |
| 出力失敗後の継続 `run_g3_snapfail` | `run_g0_noreq_a/b` | 1.01 | `PASS` |

これは**既存ツールの許容内だったという事後評価**です。境界出力・残差の数値を含めた V1 全体の合格でも、統計的同等性の証明でもありません。候補は各1本です。また、最初の面流束と入力状態のダンプは、旧・新・キーあり・最終バイナリで SHA-256 が一致しました。

`case/71.m55_combustion_wt/run_0010_ns_b0/res_86000.h5` では、既存の `eta_line` と `line_quantities` を使って次を再現しました。以下は**保存時点の診断値**です。

| η | 最大絶対偏差 | 片振幅 |
|---:|---:|---:|
| 0.03 | 0.062304 % | 0.057286 % |
| 0.06 | 0.061878 % | 0.061527 % |
| 0.1 | 0.094020 % | 0.076146 % |

η = 0.03 でも片振幅は公差を超えます。第1・第2層の平均が近くても、両線の局所差は最大 **0.009893 percentage points**あり、「分布が同じ」とは言えません。

既存の η = 0.1 の系列に対する再実行結果は `OVERALL: ALL STEADY`（2,000〜88,000 step）。これは η = 0.03 の定常性を保証しません。`run_0010` の保存済み収束判定も `NOT CONVERGED` のままです。

V1 の改訂文は、次の内容を推奨します。

- **(a)** キー無しでは数値計算経路・既存出力の仕様を維持する。再実行で決定的な項目はビット一致を要求し、非決定的な場は `check_field_regress.py` の固定した比較量・step・許容で判定する。
- **(f)** 出力処理は、次の数値更新が読む状態を変更しないことを直接検査する。別 run の全場ビット一致は要求しない。
- **(g)** 成功・失敗時の既存の動作条件を維持し、数値への影響は (f) と同じ状態不変の検査で確認する。継続場の比較は非退行の補助検査とする。

現在の600 step の結果は改訂の参考資料とし、新条件の検証結果とは区別してください。

結論: **V1 の合格確定を保留し、同じ数値状態で随時出力の有無だけを変える、出力処理前後の状態不変 A/B を次に行う。**

第 1 仮説: 観測された場の差の主成分は既存の実行非決定性であり、追加出力による数値状態の直接変更ではない。確度: **中**。  
根拠: 旧バイナリ反復でも差があり、上記4比較はノイズ床2倍以内。追加処理は [monitorOutputs.cpp:509](/home/sano/work/forge-cardmesh/solver_density_cuda/output/monitorOutputs.cpp:509)・[同:521](/home/sano/work/forge-cardmesh/solver_density_cuda/output/monitorOutputs.cpp:521) で読み取りコピーを行う。ただし、原因が `atomicAdd` だけであることは未確認。  
反証条件: 数値更新を挟まない出力処理だけで、次の更新が読む状態に1 bitでも変更が出る。

第 2 仮説: 強制全場出力で更新される host 側の写し等が、後続計算に影響している。確度: **低・未確認**。全場出力は実際に host 配列を更新するため、device 配列だけの検査では不十分（[output.cpp:119](/home/sano/work/forge-cardmesh/solver_density_cuda/output/output.cpp:119)）。

判別 A/B: 同じ node・GPU・陰解法の状態を、通常出力日でも最終 step でもない位置で同期・退避する。**A は要求なし、B は `SNAPSHOT_REQUEST` あり**とし、出力処理だけを実行して、次の更新が読む device／host 状態を比較する。数値更新は両腕の間で **0 step**。既存の [pdeSnapshot／pdeRestoreMismatch](/home/sano/work/forge-cardmesh/solver_density_cuda/main.cpp:2247) を再利用するが、現在は大きさ不明の面配列を省くため、対象配列の被覆と必要な host 状態の検査を補う。→ **差が0なら、試験した経路の直接変更説を除外**／**差があれば、非決定性だけで説明する説を棄却**。差0でも将来の全軌道のビット一致までは主張しない。

やらない方がよいこと: cell・陽解法へ変更して代替合格にすること、観測差に合わせて許容を広げること、η = 0.03 の単一保存場から定常性や格子独立性を主張すること。

呼び出し側の前提への異議: 「NS の軸と第1層の差が Euler より大きい → NS の軸節点の癖」は原因の確定になっていません。方程式・壁・格子が異なります。また、η 固定で保てるのは正規化した位置であり、壁形状を変えた際の実寸半径ではありません。

不足情報: 正式な V1 完了には、上記の状態不変検査、境界・残差を含む改訂条件の結果が必要です。新しい線の定常性には η = 0.03 自身の時系列と VERDICT が必要です。いずれも今回の「次の一手」の選択は妨げません。

証拠: 足りる
