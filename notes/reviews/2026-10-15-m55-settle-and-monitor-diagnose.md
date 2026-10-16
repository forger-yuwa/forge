# codex 諮問 (diagnose): m55-settle-and-monitor

- **brief**: [`notes/reviews/briefs/2026-10-15-m55-settle-and-monitor.md`](../../notes/reviews/briefs/2026-10-15-m55-settle-and-monitor.md)
- **plan**: [`plans/active/tooling-nozzle-card-local-mesh.md`](../../plans/active/tooling-nozzle-card-local-mesh.md)
- **date**: 2026-10-15
- **commit**: `c0215c6e` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 16.3 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-m55-settle-and-monitor-diagnose.codex.log` (codex: ok)
- **結論**: 4 係数化や追加 NS に進む前に、`run_0010_ns_b0/res_86000.h5` を固定し、δ_E の抽出帯位置だけを変える A/B を行う。** §4.4 の方針は採用できるが、現案のまま実装すると準定常の誤判定と停止時の全場欠落が起こる。
- **extra**: `case/71.m55_combustion_wt/README.md`, `solver_density_cuda/tools/check_convergence.py`, `design/forge_design/evaluate/stop_rule.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **4 係数化や追加 NS に進む前に、`run_0010_ns_b0/res_86000.h5` を固定し、δ_E の抽出帯位置だけを変える A/B を行う。** §4.4 の方針は採用できるが、現案のまま実装すると準定常の誤判定と停止時の全場欠落が起こる。

指定プロンプト全文と関連コード・保存済み判定・CSV・最終 HDF5 を確認した。ファイル変更、forge 起動は行っていない。

| 方針 | 採否 | 理由 |
|---|---|---|
| NS 本段の B0 既定化 | 採用 | ユーザ決定の範囲どおり。前段と Euler は別設定を維持する |
| 低下桁数から「残差・収支・目的量の準定常」へ変更 | 採用 | ただし log10 残差への現行判定の直接適用は却下 |
| 境界を密、全場を疎に出力 | 採用 | 出力量の定義と、停止時の最終全場の保証が必要 |
| θ_r・δ の時系列を停止条件から外す | 採用 | ユーザ了承済み。残す観測量で遅い変化を監視する |
| 保持を 4 万から 2 万 step に短縮 | 要再検証 | 旧停止条件の実績を、新停止条件の保持時間へ直接移せない |
| 初回本段の limiter 基準値を後続 pass に継承 | 採用 | 基準長・保存精度・継承範囲も明示する |
| 直接 δ_E 経路を既定から外す | 採用 | opt-in は残す。削除の根拠にはしない |

**確認した観測事実**

根拠 run は `case/71.m55_combustion_wt/run_0008_ns_b0/` と `case/71.m55_combustion_wt/run_0010_ns_b0/`。恒久索引は [case README の計算 run 一覧](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/README.md:21)。

- `run_0008`：保存済み本段判定は `PASS`。停止は 90,000 step。
- `run_0010`：保存済み本段判定は **`NOT CONVERGED`**。`--from-floor` も同判定。新方式で合格したとはまだ言えない。
- V5c 系列は、既存 `check_quasisteady` の判定関数で停止時点までに切り直しても、両 run とも **`VERDICT: STEADY`** を再現した。
- `run_0010/res_86000.h5` の η = r/r_w = 0.1、x = 4.0〜7.95 m は、最大絶対偏差 **0.094020%：PASS**、半振幅 **0.076146%：FAIL**。出口コア M は **5.5028409：V5b FAIL**。
- 最終場の非有限は 0。保存済み壁解像判定は **PASS**、メッシュ品質判定も **PASS**（AR 最大 4848.6、緩和上限 5000、skew 最大 0.421）。

主要根拠は [本段の収束判定](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_runs/run_0010_ns_b0/judge/check_convergence_main.txt)、[V5c の準定常判定](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_runs/run_0010_ns_b0/judge/check_quasisteady_v5c.txt)、最終場 `res_86000.h5`。準定常と設計精度の合否は別である。

**Major 1 — log10 残差に現行 `check_quasisteady` をそのまま使うと誤判定する。**

[classify:280](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/check_quasisteady.py:280) は、傾き・振幅を「入力系列の平均の絶対値」で割る。対数残差には不適切である。既存関数に、2 万 step で残差が 1e−8 → 1.2e−8 と **20% 増える系列**を渡した結果：

```text
log10(rms)     → STEADY  drift=1.0%
log10(rms) + 8 → DRIFTING drift=200.0%
```

残差の単位を変えただけで判定が変わる。また、既定の物理量の drift 許容 5% は、M = 5.5 なら ΔM = 0.275 に相当し、今回の設計公差より大幅に緩い。

対案は、**対数残差には絶対的な decade 幅、物理量には設計公差から決めた固定尺度**を使うこと。初期検証用の候補は次のとおり。ただし、既定として検証済みの値ではない。

- 窓 W = 20,000 step。2,000 step 出力なら端点込み 11 点以上を要求し、欠測を黙って詰めない。
- 残差の傾きは |傾き|W ≤ log10(1.05) ≈ 0.0212 decade。振幅は例えば p95−p5 ≤ 0.1 decade とし、別に生系列の突発的増大を検査する。
- 2W の長い窓でも上昇を検査し、4 万 step で 10% を超える持続的増大を通さない。全履歴の NaN/Inf 検査は残す。
- M の時間変化許容は、例えば設計公差の 1/10 から検証する。平均がゼロに近い符号付き偏差を、その平均で割らない。
- `OSCILLATING` は定常計算の停止条件では不合格。終了コードだけでなく判定名を確認する。

全有効保存量・乱流・成分の列を対象にし、恒等ゼロの非活動成分と、欠測・異常なゼロを区別する必要がある。

**Major 2 — 現在の境界 HDF5 は、全場の境界ノード値の抜粋ではない。**

[output.cpp:425](/home/sano/work/forge-cardmesh/solver_density_cuda/output/output.cpp:425) が保存するのは `bc.bvar`。特に出口の `Ps` は **規定背圧の入力値**であり、超音速流出時の実際の圧力へ更新されない。[boundaryCond_d.cu:570](/home/sano/work/forge-cardmesh/solver_density_cuda/cuda_forge/boundaryCond_d.cu:570)

したがって、境界 HDF5 の `Ps`・`ro` から音速を再構成すると、超音速出口でも不整合な M を作り得る。`Pt`・`Tt` も測定された全圧・全温とは限らない。現在の出力には音速と局所組成も不足する。

対案は用途を分けること。

- **M・流れ角**：比較対象と同じ内点／境界ノードの状態から、速度・`sonic` を保存する。
- **離散保存則の収支**：ソルバが実際に使った境界 `massflux` を積算する。[convectiveFlux_boundary_d.inc.cuh:188](/home/sano/work/forge-cardmesh/solver_density_cuda/cuda_forge/convection/convectiveFlux_boundary_d.inc.cuh:188)
- 現行の折れ線積分は「断面流量の診断値」として残せるが、数値流束の総和と同一とは扱わない。

非平面出口では −ρv dx を残し、接続順序・向き・重複面を処理する。軸ノードの r = 0 だけを見て、そのノードに属する有限の双対面流束までゼロにしてはいけない。軸対称の r 重みと 2π の二重適用も避ける。

最終全場から既存関数で得た入口・出口流量は 121.359227 / 121.354974 kg/s、相対差約 3.50e−5。ただし、これは新しい境界出力方式の同値性を証明しない。

**Major 3 — 「全場 20,000 ごと＋最後」は、現在の停止方法では保証されない。**

[output.cpp:320](/home/sano/work/forge-cardmesh/solver_density_cuda/output/output.cpp:320) の「最後」は `nStepOuter` 到達を意味する。一方、監視側は [stop_rule.py:784](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/stop_rule.py:784) で SIGTERM、必要なら SIGKILL を送る。

例えば 86,000 step で停止条件を満たしても、全場間隔が 20,000 なら最後の全場は 80,000 のままになり得る。そこで δ_E を抽出すると、停止判定と異なる状態を報告する。

対案は **協調停止要求 → 同一 step の全場・境界出力完了 → 完了通知 → 終了**。簡易実装なら、全場出力が完了した時点で停止条件を再評価し、その保存 step を判定対象として固定する。

現在も、`run_0010` の選定場は 86,000、V5c 保存判定は 88,000、本段残差判定は 88,519 までと終端が異なる。新仕様では残差の step 番号対応も含めてそろえる必要がある。

**Major 4 — 軸 M と総熱量だけでは、残すべき空間変化を見落とす。**

固定された軸の偏り自体は時間差で消える。しかし、**軸が定常なら η = 0.1 も定常**とは言えない。最大偏差と半振幅だけが一定でも、分布が移動する場合はある。

V5c と合わせて η = 0.1 の線を監視し、位置ごとの変化も判定する案を採用する。ただし現在の [point_probes.cu:152](/home/sano/work/forge-cardmesh/solver_density_cuda/probe/point_probes.cu:152) は `centCoords` に対する最近傍選択で、補間された η = 0.1 線ではない。既存の `eta_line` と対応する節点・補間重みを固定し、音速も出す必要がある。軸は補助監視として残せる。

壁熱量も同様で、総和だけでは局所的な増減が相殺する。現在仕様は [overview.md:1377](/home/sano/work/forge-cardmesh/methods/design/overview.md:1377) のとおり、保存形の `iface_Qf_eff` と壁節点ごとの準定常を正本とする。`qwall` は符号も意味も異なる診断値で、同一視できない。

θ_r・δ の時系列は復活させなくてよい。その代わり、残す壁量・コア分布でどこまで遅いモードを観測できるかを明示する。新規則が、既知の未整定例 `run_0002_ns` を途中で合格させないことは必須の回帰条件になる。

**Major 5 — 保持 2 万 step の根拠は、旧停止条件に限定される。**

保存系列を既存 `first_stop` で再生すると、保持だけを 20,000 に変えた停止は：

| run | 保持 20,000 | 保持 40,000 |
|---|---:|---:|
| `run_0008_ns_b0` | 70,000 | 90,000 |
| `run_0010_ns_b0` | 66,000 | 86,000 |

これは短縮候補の根拠にはなる。しかし旧条件には θ_r があり、新条件では最初の成立がもっと早まる可能性がある。

対案は **当面 40,000 を維持**。新条件で「最初に成立する時点」から再生して検証した後に、B0 の対象範囲を限定して 20,000 を採用する。今回の 2 run だけから全ノズルの既定にはしない。

**Major 6 — limiter 継承には `L_ref` と実効値の記録が不足している。**

[main.cpp:1768](/home/sano/work/forge-cardmesh/solver_density_cuda/main.cpp:1768) は ρ・P・a の体積平均に加え、未指定の基準長をメッシュ境界箱から決める。壁変更時に 3 値だけ固定すると、基準長が変わって作用素も変わり得る。また、貼り付け用ログは有効数字 10 桁で、FP64 の厳密な再現には不足する。

対案：

- 同一設計系列の restart・小さな壁変更では、ρ・P・a・`L_ref` を継承する。
- 新しい運転条件・物理問題では、初回本段で自動決定し直す。
- 実効値を機械可読形式・往復可能な精度で保存し、出所と `stage_manifest` に反映する。
- 現在必須の `main_line.limiter_ref` の検証も更新する。[runner_axismach.py:1533](/home/sano/work/forge-cardmesh/design/forge_design/evaluate/runner_axismach.py:1533)

**Major 7 — V7 の結果から「δ の過大補正が主因だった」と確定できない。**

[事前登録 V7](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:143) は、δ 誤差だけでなく **V5c の両量の半減**を要求している。

実際には最大偏差は 0.2247 → 0.0940% と改善したが、半振幅は **0.07879 → 0.07615%**。ほぼ残っている。さらに旧収束判定は `NOT CONVERGED` なので、登録条件を厳密に守れば判別不能。準定常を根拠に比較を認めても、結果 A にはならない。

対案は「平均偏差は改善したが、傾きの原因は未確定」と記録すること。同じ 3 係数の再較正で最大誤差が 0.624 → 0.612 mm にしかならない事実は、**その近似族の改善余地が小さい**証拠であり、4 係数化が Mach 分布を改善する証拠ではない。

また、[最終場の物理量検査](/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/_runs/run_0010_ns_b0/judge/final_field_check.json) は非有限 0 だが、総合表示は `NG`。P₀/Pt 最大 1.004546、入口の 10 点で T₀ > Tt + 1 K を記録している。これだけで傾きの原因とは断定しないが、「有限だから物理的にも合格」とする新ゲートにはできない。

**第 1 仮説：残りの傾きには、3 係数では表現できない有効壁形状の残差が寄与している。ただし、現行 δ_E の符号分布を信頼できることが前提。確度：中。**

根拠：δ_in−δ_E は x = 4〜6 m で負、7〜8 m で正。有効壁を `r_wall,NS−δ_E` と解釈すれば、設計壁に対して前半が細く後半が太い分布になり、M が下流へ増える観測と定性的に整合する。

ただし、x ≈ 1 m の +0.62 mm から試験部への影響は局所の面積則だけでは決まらない。V7 で半振幅が改善していないことも、主因と断定できない理由である。

反証条件：抽出帯位置だけで、判断に使っている δ 誤差の大きさ・符号分布が有意に変わること。その場合、現行 δ_E を根拠に形状自由度を追加する判断を撤回する。

**第 2 仮説：残りの δ 誤差には抽出帯位置への依存が含まれる。確度：中〜低。**

過去の `refine=4/32` は固定帯内の積分標本化を検証したもの。帯の位置への感度までは除外していない。[plan §5.1 #3k](/home/sano/work/forge-cardmesh/plans/active/tooling-nozzle-card-local-mesh.md:106)

**第 3 仮説：NS コア格子または出口境界の影響が残る。確度：低、未除外。**

Euler の平坦さだけでは、異なる NS 格子・境界層を持つ計算の離散化誤差を除外できない。最終出口には壁を含め 14/92 点の M < 1 があり、「出口は全点超音速なので背圧は無関係」という除外も成立しない。B0 が原因だと示す証拠はない。

**判別 A/B：同一保存場の抽出帯位置だけを変更する。forge は 0 step。**

- 固定：`run_0010_ns_b0/res_86000.h5`、参照 `run_0009_euler/res_12000.h5` と SHA、一定 x への補間、`refine`、平滑化 λ = 1・ノット 2 r_t・重み。
- A：現在選択された帯位置 y_b。
- B：各断面の帯位置を **1.25 y_b** に固定して再抽出する。帯選択の梯子は再探索しない。
- 既存の `deltastar_from_core_matched_euler(..., yb_sens_factors=[1.25])` と `smooth_delta_quintic` を再利用する。[deltastar.py:679](/home/sano/work/forge-cardmesh/design/forge_design/metrics/deltastar.py:679)
- A が保存済み δ_E を再現することを前提とし、x = 0.5〜7.9 m の平滑化後の差を比較する。事前の許容は **max|δ_E,B−δ_E,A| ≤ 0.05 mm**、かつ x = 4〜6 m の負・7〜7.9 m の正が維持されること。

**許容内なら**、今回の残差分布を説明する主要因としての帯位置依存を退け、形状残差仮説を残す。**許容外なら**、第 1 仮説の「現行 δ_E を形状誤差の根拠にできる」という前提を退け、4 係数化を止める。欠測・根なしなら比較不成立とし、合格扱いしない。どちらの結果も、Mach 傾きの因果まで確定するものではない。

やらない方がよいこと：同じ 3 係数の再較正を繰り返すこと、δ の最大誤差だけで 4 係数目を追加すること、B0 採用を理由に保持を直ちに短縮すること、軸 M の定常だけで V5c の定常を代用すること。

呼び出し側の前提への異議：**「境界 HDF5 は全場の代用品」「初期残差が小さいから今回の未収束判定は無視できる」「平均の改善で過大補正の主因説が確定した」には同意しない。** 新しい停止基準への変更自体には賛成するが、過去の判定を書き換えることと、新基準を検証することは分けるべきである。

証拠: 足りる（設計上の欠陥の指摘と、次の A/B の選定について。傾きの真因確定には至らない）。
