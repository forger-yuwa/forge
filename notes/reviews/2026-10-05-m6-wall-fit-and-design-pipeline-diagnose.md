# codex 諮問 (diagnose): m6-wall-fit-and-design-pipeline

- **brief**: [`notes/reviews/briefs/2026-10-05-m6-wall-fit-and-design-pipeline.md`](../../notes/reviews/briefs/2026-10-05-m6-wall-fit-and-design-pipeline.md)
- **plan**: [`plans/active/verification-m6-axis-wave-mesh-su2.md`](../../plans/active/verification-m6-axis-wave-mesh-su2.md)
- **date**: 2026-10-05
- **commit**: `ed1fab3e` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 8.6 min, rc=0
- **結論**: ?
- **extra**: `plans/active/tooling-design-problem-campaign-recipe.md`, `plans/accepted/tooling-nozzle-axismach-physical-throat.md`, `case/45.isobutane_m6_d155/moc_wall_fit_probe.py`, `case/45.isobutane_m6_d155/_band_ab/moc_wall_fit_probe.json`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：**Critical なし。現案の生産採用は保留**。ファイル変更・forge 起動は行っていない。

| 対象 | 重大度・採否 | 根拠と対案 |
|---|---|---|
| A1：`2400 + joint_fit(λ=1e-7)` を生産化 | **Major／要再検証** | 全区間では位置・角度ゲートに不合格。既存の集計は x<0.5 を除いている。[`moc_wall_fit_probe.py:38`](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/moc_wall_fit_probe.py:38)。下記のスロート局所診断を先に行う。位置と壁角を同時に合わせる方向自体は採用。 |
| A1：始点・出口・物理壁の拘束 | **Major／修正して採用** | 現ケースの壁始点は **x₀=0**。旧構成の放物線区間は存在しない。試作には出口半径のハード拘束がなく、`PhysicalNozzleWall` は補正後の点群から出口微分を再推定する。[試作:31](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/moc_wall_fit_probe.py:31)、[物理壁:325](/home/sano/work/forge-integ-1005/design/forge_design/geometry/wall_axismach.py:325)。設計壁の出口位置を拘束し、補正後は δ_r′・δ_r″を含む端条件と全接合を検査する。物理壁の出口角を無条件にゼロへ固定してはいけない。 |
| A2・B3：`run_0051` を全面的な合格基準にする | **Major／却下** | 再実行した判定は **残差 `NOT CONVERGED`、壁解像 `FAIL`**。波・オーバーシュートの時系列も後述の判定では `DRIFTING`。歴史的な回帰参照として保持し、設計資格の合格とは分ける。 |
| A3：現在の r″高周波指標 | **Minor／改訂を採用** | 偶数窓の移動平均と端のゼロ化が入り、実曲率変化と数値振動を分離できない。[試作:46](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/moc_wall_fit_probe.py:46)。さらに報告ツールは壁 CSV を別のスプラインで再補間している。[`nozzle_report.py:410`](/home/sano/work/forge-integ-1005/design/forge_design/report/nozzle_report.py:410)。実際の壁関数の解析微分を評価する。 |
| B1：C2 を新 recipe にする | **Major／条件付き採用** | `axismach.contur_c2/v1` を追加し、r_t 解き・k_f 較正・再測定・停止判定は recipe 内部に置く。campaign は候補探索と評価順を持つ。既存の[責務分担 §4.1–4.3](/home/sano/work/forge-integ-1005/plans/active/tooling-design-problem-campaign-recipe.md:66)と整合する。ただし固定寸法の連成と入力成果物の来歴を明文化する必要がある。 |
| B2：壁変更を先に生産へ入れてから再現化 | **Major／順序を変更** | 評価量の定義と旧基準を固定する前に壁・物性・実行経路を変えると差を帰属できない。TP の ω 床初期化にも CPG 逆算が残る。[`runner_axismach.py:886`](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:886)。旧基準の保存→評価量統一・熱力学修正→最小 C2 campaign→検証済み壁表現→寸法固定・探索、の順を推奨する。 |

**結論:** 次は、同じ `n_axis_inv=2400` の MOC 点群と細分したスロート側ノットを固定し、正則化だけを変える形状 A/B で、全区間の忠実度・滑らかさ・物理壁への引き継ぎを検証する。

**第 1 仮説:** 現案の形状ゲート不合格には、スロート付近の表現自由度不足と正則化が寄与しており、`n_axis_inv` 増加だけでは解消しない。 **確度: 高**

根拠：現作業ツリーの `design_chain` と試作の `joint_fit` を、書き込みを伴わず呼び出して再測定した。

| n_axis_inv | h₀ | λ | 全 MOC 点上 max｜Δr｜[r_t] | 全 MOC 点上 max｜Δθ｜ |
|---:|---:|---:|---:|---:|
| 1200 | 0.05 | 1e-7 | 4.915e-5 | 0.05211° |
| 2400 | 0.05 | 1e-7 | 4.737e-5 | 0.05091° |
| 2400 | 0.05 | 0 | 1.852e-5 | 0.03872° |
| 2400 | 0.0125 | 0 | **3.905e-7** | **0.001757°** |
| 提案ゲート | — | — | ≤5e-6 | ≤0.005° |

`2400・λ=1e-7` の位置最大誤差は x=0.0520、角度最大誤差は x=0.0258。ブリーフの帯別集計では両方とも対象外だった。なお、同じ点群の**現行補間壁も全点角度誤差は 0.05039°**なので、近喉部の問題すべてを新しい当てはめのせいにはできない。

h₀だけを細かくした比較では、制御点数 222→258 により位置・角度ゲートを満たした。ただし、これは**忠実度の確認であり、滑らかさ・単調性・CFD 性能の合格ではない**。

反証条件：修正版の全区間評価で上記の局所細分効果が再現しない、または点間の曲率振動・非単調性が増えて必要な滑らかさを満たせない場合、この細分案を対策として棄却する。

**第 2 仮説:** 下流の r″の毛羽には全点通過補間による増幅が寄与する。 **確度: 中**。提示 JSON は支持するが、現在の重み `w=ds/ds.mean()` では点数とともにデータ項が増え、固定 λ の相対的な強さが変わるため、解像度比較には正則化強度の変化も混ざる。[試作:27](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/moc_wall_fit_probe.py:27)。

**第 3 仮説:** C2 の出口較正は他の Mach・寸法にも有用だが、M6 の一回更新だけでは一般化できない。 **確度: 低・未確認**。r_t 比の −0.2 乗換算は初期予測として扱い、最終 NS の再測定で閉じる必要がある。

**判別 A/B:**

- 共通：同一 MOC 点群、`n_axis_inv=2400`、h₀=0.0125、h₁=0.5、同じ重み・端拘束・δ_r 関数。
- 変更は **λだけ**。A＝0、B＝1e-9。
- 長さ：**CFD 0 step**。設計壁と `PhysicalNozzleWall` 再構成まで。
- 判定：全 MOC 点で |Δr|≤5e-6 r_t・|Δθ|≤0.005°。加えて全ノット区間の点間極値・単調性、始点と出口、接合の C⁰/C¹/C²、物理壁の解析微分を検査する。端部をゼロ埋めして合格させない。
- **A 合格・B 不合格なら** λ=1e-9 でも正則化が忠実度を壊している。**両者合格で B が曲率振動を減らすなら** B を CFD 検証候補にする。**A も点間・物理壁ゲートで落ちるなら**、点上誤差だけで局所細分を採用できるという仮説を棄却する。

**やらない方がよいこと:**

- 現状の帯別 JSON を根拠に生産壁を切り替える。
- 1200→2400、壁表現、物性更新、壁メッシュ変更を一度に入れて差を評価する。
- 高解像度の当てはめ壁を真値とみなす。同じ近似空間に由来する共通誤差が残る。
- `run_0051` に数値が近いことだけで新 recipe を `accepted` にする。
- 再実行差が大きいほど許容差を広げ、物性変更や未定常まで「ノイズ」として吸収する。

**呼び出し側の前提への異議:**

1. **「最終設計は全項目合格」は範囲を限定する必要がある。**

   対象は `case/45.isobutane_m6_d155/run_0051_ns_final_c2/`。作業ツリーには run がないため、元の `/home/sano/work/forge/` 配下の同名 run を読んだ。保存来歴は commit `bc1c84d1`、バイナリ SHA256 `6d5933e1…`。

   - `check_convergence.py`、step 0–11999：**`NOT CONVERGED (stalled/plateau)`**。
   - 保存済みメッシュ判定：**`SOFT-PASS`**、最大 AR 1208.7。
   - `check_wall_resolution.py --groups wall`、step 12000：**`FAIL`**。y₁⁺最大 14.881、1 超の面積割合 **31.0%**。最大位置 x≈−0.02495 m。
   - step 1000–12000 の12出力について `VALUE/*` の NaN/Inf は 0。
   - 報告ツールと同じ量を抽出し、`check_quasisteady.classify_series` の既定閾値・末尾40%で判定：オーバーシュート **`DRIFTING`**、波指標 **`DRIFTING`**、出口コア M **`STEADY`**。出口コア M は単調減少で、漸近推定 **5.99992**。

   末尾5枚の変動が登録済みの絶対許容幅に入ることと、量自身が定常化したことは別である。旧用途判定を遡って書き換える必要はないが、新 recipe の全面的な資格証明にはできない。参照台帳は [case README:63](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:63)、根拠成果物は同 run の `residual_history.csv`、`res_1000.h5`～`res_12000.h5`、`res_wall_3_12000.h5`。

2. **A2 の評価量を先に固定する。**

   現在の「波」は圧力そのものではなく、Mach 偏差から P-spline トレンドを引いた量である。ノットも評価窓端に依存する。[`nozzle_report.py:206`](/home/sano/work/forge-integ-1005/design/forge_design/report/nozzle_report.py:206)、[`deltastar.py:296`](/home/sano/work/forge-integ-1005/design/forge_design/metrics/deltastar.py:296)。

   共通の評価座標・窓・ノット原点・重み・正規化・版を固定し、実圧力の指標も追加する。出口コア M の現報告値は節点算術平均なので、質量流束平均へ変更するなら旧結果も同じ関数で再抽出する。

   CFD 比較は形状ゲート通過後、同一 n・同一ビルド・同一実効設定で行う。再現許容差内でも、波≤0.01%、オーバーシュート≤+0.035%、出口コア M＝6±0.02%などの**用途上の絶対ゲートは独立して維持**する。

3. **A3 はスケール依存の局所トレンド除去を推奨する。**

   実際の壁関数の r″に対し、対称な局所多項式フィルタで複数の固定窓幅を評価する。端部は別判定にし、自然な曲率変化まで「ノイズ」と断定しない。高解像度との差は離散化感度の補助指標にする。

   新指標に旧閾値 1e-5 をそのまま移植せず、滑らかな基準曲線と既知振幅・波長の擾乱で応答を確認してから事前登録する。

4. **B1 の寸法固定は r_t 解きだけでは完結しない。**

   `from_length` の `L_total` は r_t 単位で、設計スロート起点。[`runner_axismach.py:386`](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:386)。実寸長・入口径を固定するなら、r_t 更新ごとに `L_total` と `r_inlet` も更新し、長さの起点が設計スロートか物理スロートかを定義する。

   実際、現在の YAML の r_t と `r_inlet` から得る入口半径は **0.4980966 m**で、コメントの0.5 mとは一致しない。[problem:29](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/problem_d155_ns_c2final.yaml:29)。

   recipe は寸法残差と出口 δ 残差をともに検査し、上限反復で未達なら `rejected`、必須成果物欠落なら `incomplete`。凝縮評価は確定した物理壁を参照する後段として campaign に配置する。

5. **B2・B3 は「旧結果の再現」と「修正後の設計資格」を分ける。**

   推奨順は、旧入力・物性・壁・メッシュ・初期 snapshot・バイナリのハッシュ保存と M6 の再現幅取得 → #4 評価量統一 → #3 熱力学修正 → #6・#7 の最小契約 → search なしの C2 recipe → 壁変更 → 寸法固定・Pareto 探索。

   同一条件の再現許容差は量ごとに **max（3反復の最大差×3、事前登録した絶対下限）**とする。既存 plan の M＝1e-5、流量＝相対1e-4、T＝0.01 K は候補にできるが、別 case のノイズ実測を M6 に流用しない。波・オーバーシュートの下限は別途必要。許容差が用途上の余裕を食い潰す場合は判定不能とする。

   main の物性変更や TP 初期化修正による差は再現ノイズに含めず、変更影響として報告する。これらを直した recipe に、旧 `run_0051` との数値一致を必須条件として課してはいけない。

**不足情報:** M6 の同一条件3反復による評価量別の再現幅、固定寸法の許容差と長さの起点、出口 δ の較正が一点値か帯平均かという統一契約、C2 最終 r_t・k_f 解きの実行可能な手順が不足している。提示された `calibrate_contur.py` は全域3係数較正、`prep_contur_cal.py` は受け取った係数で準備するスクリプトであり、C2 全手順の再現器にはなっていない。

**plan 未反映**。依頼の read-only 条件に従った。呼び出し側で `verification-m6-axis-wave-mesh-su2.md` §5.1 #12・§9、および `tooling-design-problem-campaign-recipe.md` §4・§5.1・§6へ採否と検証条件を反映すること。
