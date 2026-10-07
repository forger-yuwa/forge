# codex 諮問 (diagnose): ns-n012-exitM-deficit

- **brief**: [`notes/reviews/briefs/2026-10-07-ns-n012-exitM-deficit.md`](../../notes/reviews/briefs/2026-10-07-ns-n012-exitM-deficit.md)
- **plan**: [`plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md`](../../plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md)
- **date**: 2026-10-07
- **commit**: `34fbc2ec` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 6.0 min, rc=0
- **結論**: **較正値・k_f・r_t・判定条件を固定し、N0・N1・N2 を登録どおり各 20000 step、1 回だけ延長する。**
- **extra**: `plans/active/discretization-moc-axis-limit-and-corrector.md`, `plans/active/verification-case45-euler-total-enthalpy.md`, `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 登録どおり 3 本を各 20000 step 延長 | **採用** | [upstream plan:173](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md:173) の規定どおり。目的は残る準定常未達の確認。出口 M の回復は保証しない。 |
| H1：「かさ上げが消えた」、応答係数 0.93 | **要再検証／Major** | ΔM/Δoffset = 0.92953 は再計算できた。しかし旧生産→N0 では `r_t` も 0.0766539→0.0766653655 m に変わり、初期化・段階起動も異なる。[問題記録:20](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/ns_n012_problems.json:20)、[run 索引:105](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:105)。**較正変更との整合はあるが、単独の感度係数ではない**。全温異常だけの寄与とも断定しない。 |
| H2(a)：NS から offset を約 +1.6e−3 補正 | **現方針では却下／Major** | 現在仕様は「NS では較正し直さない」。過去の中心合わせ案も撤回済み。[現在仕様:361](/home/sano/work/forge-integ-1005/methods/design/overview.md:361)、[採用済み plan:227](/home/sano/work/forge-integ-1005/plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md:227)。これは既存較正の継続ではなく設計方針の変更。交絡した 0.93 を使った補正は行わず、延長後も未達なら候補 (c) に戻す。 |
| H2(b)：直ちに C2 をやり直す | **要再検証／Major** | δ_E/δ_C ≈ 1.0017 は現在の数値ゲート内。旧値からの +0.00180 には Euler 参照変更も含まれ、境界層の物理的変化と分離されていない。E5 は参照時刻への感度を調べる試験であり、旧参照との差の原因まで識別しない。[Euler plan:117](/home/sano/work/forge-integ-1005/plans/active/verification-case45-euler-total-enthalpy.md:117)。まず E5 の結果を確認し、C2 が出口 M を回復させるとは先取りしない。 |
| H3：近零量への相対判定が未達を生む | **機構は採用、救済への使用は却下／Major** | 判定器は `abs(mean)` で正規化するため、この説明はコードと整合する。[check_quasisteady.py:280](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_quasisteady.py:280)。ただし、値が上限より小さいことは準定常の証明ではない。E2 の絶対幅条件を U4 に事後移植せず、今回の VERDICT を維持する。 |
| 差が時間幅以下なので U4・V5′ を採用 | **却下／Major** | 評価器自身が時間幅を「差の信頼区間ではない」と定義している。[ns_n012_eval.py:446](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/ns_n012_eval.py:446)。共通の絶対未達だけで変更を不採用ともできない。**候補は維持、生産採用は保留**。N2 のオーバーシュート低下も未整定系列の参考差として残す。 |

結論: **較正値・k_f・r_t・判定条件を固定し、N0・N1・N2 を登録どおり各 20000 step、1 回だけ延長する。**

第 1 仮説: 出口 M の未達は、上流多項式化・MOC 変更に固有ではない共通の不足で、登録済みの追加 20000 step だけでは解消しない。確度: **中**。

根拠: [評価 JSON](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/ns_n012_eval.json:140) に記録された本段 60000〜80000 の結果は次のとおり。

| 根拠 run（リポジトリ相対パス） | 出口 M の窓平均 | 下限との差 | 出口 M の VERDICT |
|---|---:|---:|---|
| `case/45.isobutane_m6_d155/run_0165_ns_n012_N0/` | 5.99852654 | −2.73460e−4 | `STEADY`、単調減少・漸近推定 5.99851 |
| `case/45.isobutane_m6_d155/run_0166_ns_n012_N1/` | 5.99851269 | −2.87312e−4 | `STEADY` |
| `case/45.isobutane_m6_d155/run_0167_ns_n012_N2/` | 5.99853980 | −2.60200e−4 | `STEADY` |

各窓の最大値でも下限を 2.50〜2.79e−4 下回る。3 本とも、本段区間の残差判定は **`NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)`**。これは収束解の認定ではない。

出口 M は最終断面の η∈[0.05, 0.7] の節点平均で、Euler 参照はその算出式に入らない。[exitM_sampling_ab.py:24](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/exitM_sampling_ab.py:24)。したがって、**Euler 参照ファイルの差し替えだけによる後処理上の M 低下**ではない。一方、共通不足の内訳が境界層補正・離散化・初期化履歴のどれかは未確定。

反証条件: 設定を変えない延長後、出口 M が登録帯 5.9988〜6.0012 に入り、対応する準定常条件も満たす。その条件について「追加 20000 step では解消しない」を棄却する。

第 2 仮説: 旧較正から新較正への変更が、旧生産→N0 の低下の主要因である。確度: **中**。符号・大きさは整合するが、寸法変更との寄与分離は未確認。0.93 は補正式に使える係数ではない。

第 3 仮説: オーバーシュートの準定常未達には、近零の平均による正規化が強く効いている。確度: **高（判定機構について）**。実際の時間変化が無視できるかは未確認。

判別 A/B: **変える要因は本段の計算長だけ**。A＝既存 80000 step、B＝同一設定で 100000 step。各 `res_80000.h5` から新しい run へ `restart_field.py` で継ぎ、soft 段を挟まず 20000 step、5000 step ごとに保存する。実装はこの設定差だけを検査している。[ns_n012.py:727](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/ns_n012.py:727)。

- B の判定窓は登録済みの通算 **80000〜100000、5 枚**。出口 M・オーバーシュート・波・δ_E/δ_C と全残差を再判定する。
- **不足が持続する結果なら**、「登録の追加 20000 step で出口ゲートを回復できる」を棄却し、候補 (c) として未達のままユーザ判断へ戻す。
- **出口ゲートと準定常条件を満たす結果なら**、第 1 仮説を該当条件について棄却する。
- 準定常も満たさない結果は判別不能。再延長や窓変更で救済しない。この A/B は「較正が真因」を証明する試験ではない。

やらない方がよいこと: NS の欠損を使った `Md_moc_offset` の即時補正、C2 と較正の同時変更、相対差だけでの生産採用、未達を合格と読み替えた凝縮への移行。旧壁についてのユーザ決定を、今回の新しい未達への包括的免除にもしない。

呼び出し側の前提への異議:

- **Minor**：「他のゲートは合格」は広すぎる。N0 の波は値の上限を満たすが、準定常は **`DRIFTING`**。オーバーシュートは N0/N2 が **`DRIFTING`**、N1 が **`OSCILLATING`（0.002292 ± 0.00017 %）**。[評価 JSON:140](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/ns_n012_eval.json:140)。値の合否と準定常の合否を分けて記録する。
- 漸近値 5.99851 は増分減衰に基づく外挿であり、延長結果の保証ではない。「回復しそうにない」は支持できるが、「延長で動かない」とは断言しない。

不足情報: 当該 run の HDF5・残差 CSV・全系列 CSV・元の判定ファイルはこの作業ツリーに無く、独立した再判定はできない。評価器の SHA256 と問題 YAML 6 本のハッシュは保存記録との一致を確認した。**E5 の結果は未提供**。また、ローカルの `e4_recal_eval.json` は d1 が「未実施」の旧版で、E4V 合格は plan の記録までの確認に留まる。

run の恒久索引は [case README:113](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:113)。ファイル変更・forge 起動は行っていない。**本諮問の提案は plan 未反映**で、呼び出し側が upstream plan §6 U4・§9 に採否を記録する。
