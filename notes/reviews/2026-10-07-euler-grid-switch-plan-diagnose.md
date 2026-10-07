# codex 諮問 (diagnose): euler-grid-switch-plan

- **brief**: [`notes/reviews/briefs/2026-10-07-euler-grid-switch-plan.md`](../../notes/reviews/briefs/2026-10-07-euler-grid-switch-plan.md)
- **plan**: [`plans/active/verification-case45-euler-total-enthalpy.md`](../../plans/active/verification-case45-euler-total-enthalpy.md)
- **date**: 2026-10-07
- **commit**: `670bb0fc` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.8 min, rc=0
- **結論**: **Euler 専用設定と上記の前提を登録し、legacy MOC・新配点で E4 の出口較正を先に行う。**
- **extra**: `plans/active/discretization-moc-axis-limit-and-corrector.md`, `plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0／Major 6／Minor 1）

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| Major | **Euler 専用設定を設ける案を採用。YAML の値だけ直す案は却下** | 現在は両経路とも `p.mesh` を読む。[runner_axismach.py:647](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:647)。`mesh_euler` を独立した設定とし、Euler の既定を全断面 `wall_first_frac=0.005`、スロート別指定・軸側 cap なしにする。既存 `mesh` は NS 用として維持し、Euler への暗黙継承を禁止する。 |
| Major | **E2 を「新配点なら設計評価量も準定常」と読むことは却下** | 保存 JSON の B は全温偏差の絶対幅を満たすが、13 枚窓の正式判定は `DRIFTING` 25 列・`OSCILLATING` 2 列。残差も `NOT CONVERGED`。出口 M・波・傾きの準定常は別に確認する。 |
| Major | **E4 は legacy MOC を固定して実施する案を採用** | E2 の B は既に `analytic + converge`。[euler_t0_e2.py:65](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/euler_t0_e2.py:65)。これを legacy の較正基準にすると MOC の変更が混ざる。単調壁・`legacy + fixed2` を E4 の基準にする。 |
| Major | **出口コア M の共通標本化を採用。Euler の較正成功を NS に移せるとの前提は要再検証** | 出口指標は最終断面の η∈[0.05, 0.7] の節点単純平均で、標本位置にも依存する。[accepted plan:108](/home/sano/work/forge-integ-1005/plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md:108)。固定した共通 η 標本で較正し、自格子平均も併記する。NS の合否は NS で判定する。 |
| Major | **V5d は共通の較正値・共通の初期化手順で比較する案を採用** | `Md_moc_offset` は MOC と軸則に渡され、壁を変える。[runner_axismach.py:336](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:336)。腕ごとに較正してから比較すると、MOC と較正の複合効果になる。両腕に E4 の同じ値を使う。 |
| Major | **旧 NS をそのまま U4 の対照にする案は、較正値が変わる場合は却下** | U4 は `poly` だけを変え、比較元を `run_0147+0149` としている。[upstream plan:166](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md:166)。新較正値で legacy・ramp の対照を更新してから、poly、MOC の順に進む。 |
| Minor | **実効設定の記録と手順書の更新を採用** | `prepare_info.json` の mesh 欄は現在 3 項目だけ。[runner_axismach.py:751](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:751)。全 `Mesh2DParams`・採用元・座標／接続ハッシュを残す。[手順書:116](/home/sano/work/forge-integ-1005/procedures/nozzle-design-workflow.md:116)の「NS と同じ格子で較正」も更新する。出口半径の表記は **0.775 m ± 0.1 mm** に直す。 |

以下を §4・§6 の登録案とする。**plan 未反映**。依頼どおりファイルは変更していない。

**実装方法**

`prepare` は `mesh_euler`、`prepare_ns` は既存の `mesh` を読む。名前の違う二つのブロックを混ぜて補完しない。

既存の Euler YAML に `mesh` だけがある場合は、移行先を示して停止する。黙って NS の配点を使うことも、既存指定を無視して既定へ落とすことも避ける。両ブロックを持つ問題では、使用したブロックを明記する。ローダーから runner までの受け渡しも検証対象に含める。

case/45 の `mesh_euler` には **2000×97、現行の軸方向配置、全域 0.005** を明示する。2000×97 は今回の評価条件であり、全ケースの汎用既定にする根拠はない。移行の合格条件は、Euler が E2-B の実効座標を再現し、NS は移行前と座標・接続がビット一致すること。

**E4：出口較正**

- 単調壁、`legacy + fixed2`、凍結初期線 `run_0062`、r_t、BC、熱物性、バイナリを固定する。開始値は δ₀＝`Md_moc_offset`＝+3.770e−4。
- 両段階とも、新メッシュ上に同じ手順で作った等エントロピー IC から開始する。異常を含む G1 の Euler 場を移送しない。起動前の全節点で |T₀−1600|≤1 K を確認する。
- E2 と同じ **soft 3000＋本段 54000、1000 step ごと出力**。判定窓は本段 42000〜54000 の 13 枚と末尾 5 枚。54000 は予算であり、収束の予測ではない。
- 出口は各設計の実際の最終断面。η∈[0.05, 0.7] について、既存 NS 基準格子から採った固定 η 列を登録し、同じ線形補間・単純平均で M_common を評価する。自格子平均は別列に残す。軸上 M や x_E の値に置換しない。
- 残差は本段区間の `check_convergence` で確認する。全不合格列が停滞だけの場合は、今回の実務較正に限って許容するが、`NOT CONVERGED` を保持する。RISING、still converging、DIVERGED、入力不備は保留。
- 全温には E2 の健全性条件を課す。出口 M は両窓で `check_quasisteady` が `STEADY`、さらに **13 枚の幅と、末尾 5 枚・その直前 5 枚の平均差が各 ≤5e−5**。T₀ の絶対幅 0.1 K を Mach 数の判定に転用しない。
- 合格目標は **13 枚すべてで |M_common−6|≤1e−4**。開始値で満たせば据え置く。外れる場合にだけ、#11f と同じ係数 1 の式  
  **δ₁＝δ₀−(平均 M_common−6)**  
  で一度更新し、作り直した壁の独立した run で同じ条件を検証する。失敗時は追加補正を重ねず保留する。

合格した場を新しい Euler 参照に指定し、δ_E 抽出・C2 較正・報告の参照先を同期する。**初期線の凍結源 `run_0062` は変更しない。**

**V5d：新格子での MOC 比較**

腕 B＝単調壁・`legacy + fixed2`、腕 M＝単調壁・`analytic + converge`。**両腕とも E4 の較正値を使い、等エントロピー IC から各 3 本、計 6 本**を推奨する。条件と評価器を事前登録しておけば、E4 の最終確認 run を B の 1 本として共用できる。

旧第 7 run は番号写像 IC の影響を確認するものだった。今回は番号写像を使わないため、その検査を機械的に継承しない。ただし結論は「指定した等エントロピー初期化での比較」に限定し、IC 非依存とは主張しない。

長さ・出力間隔・窓は E4 と同じ。V5 の 6 量、Δq、腕平均・SE、D±2SE の判定は維持する。準定常は両窓で確認し、`exit_M_dev` だけの既存絶対許容例外も維持する。加えて、各量の13枚窓の幅が Δq/10 以下、符号付き出口 M の幅が 1e−5 以下であることを登録する。SE は自己相関を補正しない実務指標のままである。

全 run に E4 と同じ全温・残差の前提を課す。評価座標 X_E・X_F と出口の共通 η 列も固定する。旧 G1 の run は対照に混ぜない。

出口較正の据え置き判定と、M＝6 の達成は分ける。**腕間差が 1e−4 以下でも、腕 M 自身が E4 の絶対目標に入るとは限らない。** 腕 M が再較正を要する場合は NS を開始せず、新しい共通較正値で両腕を再確認する。共通値で両者を通せない場合、MOC 単独の比較と各設計の最適較正を同じ試験として扱わない。

**NS の順序と本数**

順序は **設定分離 → E4 → V5d → 較正値の確定 → NS 対照更新 → U4 → V5′**。

較正値が変わった場合、変更を分離する最少の構成は次の 3 条件になる。

| 条件 | 上流物理壁 | MOC | 比較の役割 |
|---|---|---|---|
| N0 | ramp | legacy | 新較正値での対照 |
| N1 | poly | legacy | N0 との差で U4 |
| N2 | poly | analytic | N1 との差で V5′ |

3 条件で較正値・k_f・r_t・NS 格子・実効設定を固定する。途中で C2 により k_f・r_t を変えた場合、その比較は単独変更ではなくなるため、対照も更新する。

U4 の本段 80000、窓 60000〜80000、未達時の追加 20000 を 1 回という枠と、性能・壁解像・量別準定常のゲートを維持する。cross-mesh 移送と段階起動も維持する。凝縮についても U4・V5′ の既存の評価義務は残る。**3 条件は最少の比較構成であり、C2 修正や凝縮を含む総 run 数が 3 本で済むという意味ではない。**

結論: **Euler 専用設定と上記の前提を登録し、legacy MOC・新配点で E4 の出口較正を先に行う。**

第 1 仮説: **新配点と G1 由来 IC の排除により、出口較正を判定できる状態が得られる。** 確度: **中**  
根拠: `case/45.isobutane_m6_d155/run_0162_euler_t0cluster_u5em3/` の保存評価では、42000〜54000 の最大全温偏差 0.103227 K、27 列の最大時間幅 0.052117 K。ただし記録された残差判定は **`NOT CONVERGED (stalled/plateau)`**、E2 全体は **`判別不能`**。これは出口 M の準定常の証拠ではない。  
反証条件: legacy の新規場が固定窓で全温・出口 M の前提を満たさなければ、「この手順と予算で較正可能」を退ける。

第 2 仮説: **Euler の新較正値だけでは G1 の NS の出口目標を満たせない可能性がある。** 未確認。配点も粘性も異なるため、NS ゲートで確認する。

判別 A/B: **E4 の補正前 A と補正後 B で、変える入力は `Md_moc_offset` だけ。** 各 soft 3000＋本段 54000、上記の出口 M・全温・残差を評価する。A が絶対目標内なら再補正不要。A が目標外で B が全条件を満たせば係数 1 の一回補正を支持する。両者の前提が成立しても B が目標外なら、その補正則を棄却する。前提未達は判別不能とする。

やらない方がよいこと: G1 の異常場の再利用、腕別較正後の差を MOC 単独効果と呼ぶこと、旧 NS と新較正・poly の NS を直接比較して U4 とすること、合格するまで窓や補正回数を変えること。

呼び出し側の前提への異議: **E2-B は「全温偏差が指定の絶対幅内」であり、「流れや設計評価量が準定常」とは未確認。** また、2000×97 をそろえても Euler と NS の離散化誤差はそろわない。

不足情報: ローカルには E2 の run 本体がなく、確認できたのは保存 JSON 内の数値・判定記録とコードである。HDF5・残差 CSV の再検査は実施していない。実測結果の台帳は [case README](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:106)。参照時の HEAD は `670bb0fc` で、ブリーフ記載の `7fb7041e` と異なるため、実装・投入時には対象 commit とバイナリを固定すること。
