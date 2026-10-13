# 帯局所抽出の等価半径の探索を積分域に限る (δ_r = r_w の誤抽出) と δ_r 表の端の float32 (#18)、M5 デモの出口 M を閉じる

## メタ

- **area**: `tooling / design (δ* 抽出・設計チェーン)`
- **status**: `draft`  <!-- 2026-10-13 起票 (ユーザ決定「ある程度できたら、2 に進んでもらっていい」= 抽出の頑健化と #18 を別 plan で起こし、codex に諮ってから実装 → NS の分布全体で壁と r_t を作り直して出口 M を 4.9975〜5.0025 に) -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) 「排除厚さ (δ*) の帯局所抽出」の節
  - [`procedures/nozzle-design-workflow.md`](../../procedures/nozzle-design-workflow.md) (NS 後の寸法 `solve_rt_throat --prev`)
- **related_plans**:
  - [`../accepted/tooling-nozzle-deltastar-core-matched-euler.md`](../accepted/tooling-nozzle-deltastar-core-matched-euler.md) (帯局所抽出 `band_local_deficit` の起源)
  - [`verification-m6-axis-wave-mesh-su2.md`](verification-m6-axis-wave-mesh-su2.md) §5.1 #7a (帯の位置の選び方 `band_select="edge"`; 別セッションの plan、本計画からは書き込まない)
  - [`../accepted/tooling-nozzle-wt-deliverables.md`](../accepted/tooling-nozzle-wt-deliverables.md) §5.1 #18・#19 (本計画へ移す)
  - [`tooling-nozzle-isothermal-wall-chain.md`](tooling-nozzle-isothermal-wall-chain.md) (冷却壁の符号付き δ_r `_negative_delta_r`)
- **created**: `2026-10-13`
- **owner**: `sano`

## 1. 目的

case/70 の M5 デモ (等温壁 350 K、断熱壁温の約 1150 K に対して強く冷やす) で、出口コア M が 4.9907 (−0.19 %) に留まった。
NS 後の δ* の分布で寸法と壁を作り直す経路 (`solve_rt_throat --prev`) が次の 2 つの欠陥で使えなかった。

1. 帯局所抽出 `band_local_deficit` が、スロート直後 (x ≈ 1〜3 r_t) の断面で δ_r = r_w (壁の半径そのもの) を返す (§3)。
2. (#18) 抽出した δ_r の表の終点が、設計壁の出口より float32 の丸めの分だけ手前で止まり、1 本の B-spline の物理壁が作れない。

完了時の状態:

1. 等価半径の探索が積分域 (壁〜帯の下端) に限られ、その区間で累積流量が単調でなければ hard の不合格になる。
2. δ_r の表の端が、同じ出口から来たことを照合したうえで設計の出口に揃う (許容を緩めない)。
3. case/70 で NS 後の δ_E 全分布から r_t と壁を作り直した最終 NS の出口コア M が 4.9975〜5.0025 (5.000 ± 0.05 %、M6 生産の判定基準) に入る。

## 2. スコープ

- **やる**: `band_local_deficit` の等価半径の探索、hard 不合格の理由の追加、`solve_rt_throat --prev` の表の端の照合と揃え、
  case/70 の NS 後の作り直しと最終 NS、試験。
- **やらない**: 帯の位置の選び方 (`edge_band_positions` の正値だけの梯子) の変更、soft 不合格の扱いの変更 (§4.1 の後に再評価)、
  `core_matched_deficit` の変更 (参照 q_ref = α q_E は正なので累積流量は単調で、同じ欠陥は無い)、CONTUR の熱閉包、M6 生産 (別セッション) の再抽出。

## 3. 関連 docs と前提 (観測事実, 2026-10-13)

- **欠陥の機構** (`design/forge_design/metrics/deltastar.py:255-283`): 参照 $q_{ref}(y)=(a+b\,(y-y_b))\,q_E$ は帯 $[y_b,(1+w)y_b]$ の線形フィットを
  **全半径へ外挿**する。累積流量 $F(y)=2\pi\int q_{ref}\,r\,dr$ (壁から) を全半径で作り、`np.interp(D, F_from_wall, r_desc)` で等価半径を引く。
  帯の傾き $b$ が大きいとコアで $q_{ref}<0$ になり、$F$ が途中から減って末尾が $D$ を下回る。`np.interp` は `x ≥ xp[-1]` で `fp[-1]`
  (軸の半径 ≈ 0) を返すので **δ_r = r_w** になる。直前の容量の検査 `D > F_zone` (帯の下端での $F$) は通っている (根は帯の中にある) ので検出されない。
  `band_select="edge"` は帯を固定 (`max_retry = 0`) するので、帯の傾きの再試行 (`slope0 > slope_max`) も効かない。
- **case/70 の実測** (`case/70.m5_air_wt_demo/ab_root_zone.py`、Euler 参照 `run_0003_euler_36k`、forge 0 step): 腕 orig = 現行、腕 zone = 探索を壁〜$y_b$ の区間の $F$ に限る。

| NS run | 腕 | δ_r = r_w の断面数 | max \|raw\| (0 ≤ x ≤ 3) [r_t] | δ_E (平滑化) x = 0 / 1 / 2 / 3 / 10 / 48 [r_t] | 流量からの実効値 δ_eff,ṁ [r_t] |
| --- | --- | --- | --- | --- | --- |
| run_0009 | orig | 99 (x 1.086〜2.831) | 1.48 | +0.2403 / +0.5821 / +0.7737 / +0.6993 / 0.0411 / 0.2857 | −0.003421 |
| run_0009 | zone | 0 | 0.0034 | **−0.003601** / −0.002972 / −0.001804 / +0.000092 / 0.0308 / 0.2857 | (同上) |
| run_0010 | orig | 38 | 1.49 | +0.0525 / +0.1959 / +0.3002 / +0.2995 / 0.0352 / 0.2867 | −0.002978 |
| run_0010 | zone | 0 | 0.0029 | **−0.003022** / −0.002336 / −0.001150 / +0.000832 / 0.0320 / 0.2867 | (同上) |

  - zone の腕のスロートの δ_E は、独立な流量からの実効値 (`massflow_ratio`: $\delta_{eff}=r_{throat,phys}-\sqrt{\dot m_{NS}/\dot m_E}$、
    `deltastar.py:93-117`) と run_0010 で 1.5 %、run_0009 で 5 % で一致する (どちらも負: 冷却壁で壁近傍の ρu が参照より多い)。
  - zone の腕で、探索区間の中の $F$ が単調でない呼び出しは 0 / 114,323、全半径の $F$ の末尾が $D$ を下回る呼び出しは 5,521 (4.8 %)。
  - ~~x ≥ 10 の値は両腕でほぼ同じ~~ **[訂正 2026-10-13、codex 諮問 Major]**: x = 10 の δ_E は run_0010 で 0.03524 → 0.03198 (**−9.3 %**)、run_0009 で 0.04107 → 0.03076。
    x = 30・48 は相対 1e-3 以内。x = 0 や x = 10 が変わるのは、帯の位置 (edge) が梯子 (帯を固定した δ_r の列) の値から決まり、
    梯子に誤った値が入っていた (κ の移動中央値 半幅 3 r_t・P-spline ノット 6 r_t で下流へ広がる) ため。
  - **欠陥は最終の raw に出なくても梯子に入る**: run_0006 (orig で δ_r = r_w の断面 0) でも zone の腕で 598 呼び出しが「全半径の F の末尾 < D」に当たり、
    縁 y* が x 0.875〜1.06 の 3 断面で変わり (NaN の有無は 16 断面で違う)、raw が 340/1600 断面で変わった (最大 3.8e-4 r_t、x ≈ 6 で 2.5 %)。
    δ_E は x = 0 で −0.000887 → −0.000884、x = 48 で同値。**上の 33 個の走査 (最終の raw だけ) は過小評価** (梯子は記録に残らない)。
  - 手元の抽出の記録 33 個を走査すると、同じ症状 (δ_r = r_w) は case/70 の 3 個と **case/45 `run_0089_ns_c2pin_pass1` の縮流部 8 断面** (x −5.08〜−4.45) にあった
    (断熱壁でも起きる; 冷却壁に固有ではない)。M6 生産の抽出は AWS にあり未走査。
- **#18 の実測**: `problem_m5_air_c2b.yaml` の設計の出口 x_e = 48.031377427 r_t、run_0010 の抽出の表の終点 48.031375885 r_t (短い分 1.54e-6 r_t = 3.9e-8 m)。
  格子の座標は float32 [m] で、x_e·S = 1.1995145 m の 1 ulp は 1.19e-7 m (4.8e-6 r_t)。抽出は x を r_t 単位の float32 で持つ (1 ulp 3.8e-6 r_t)。
  短い分は 2 回の丸めの範囲に入る。`SingleBSplinePhysicalWall.applicability` (`wall_axismach.py:1452-1455`) は表が [0, x_e] を覆うことを要求して止める (正しい検査)。
- `solve_rt_throat --prev` (`deltastar_loop.py:288-357`) は同じ問題 YAML (r_t 単位の設計が同じ) で r_t だけを解く。設計の x_e (r_t 単位) は反復で変わらない。

## 4. 設計方針

### 4.1 等価半径の探索を積分域に限る

- `band_local_deficit` の $D \ge 0$ の枝で、壁から帯の下端 $y_b$ までの区間の累積流量 $F$ (区間の端は $y_b$ で線形補間した点を足す) だけを使って $D$ を引く。
  容量の検査 `D > F_zone` はこの区間の末尾で行う (今と同じ値)。
- 区間の中で $F$ が単調増加でない (= 区間の中で $q_{ref} \le 0$) なら、理由 `nonmono_zone` を **hard** (`HARD` に追加) にして δ_r を NaN にする
  (呼び出し側は前回値を保持する; 今の hard の扱いと同じ)。今回の測定では 0 件だが、帯の外挿が壁側でも負になるような場合の保険。
- $F$ が全半径で単調な (= 今まで正しく動いていた) 断面では、`np.interp` が同じ区間の同じ 2 点を使うので**結果はビット同一**。
- `core_matched_deficit` は変えない ($q_{ref}=\alpha q_E > 0$ で単調、かつ `D > F_from_wall[-1]` を検査している)。

### 4.2 δ_r の表の端 (#18) と新しい設計への写し

**[2026-10-13 改定]** §5.1 #2 の寸法の A/B で、r_t 単位の設計を据え置く経路は棄却 (長さ +4.4 mm)。新しい設計は**実寸の長さを保って**カードから作り直すので、
r_t 単位の x_e も変わる (約 −0.4 %)。表は前の設計の x から新しい設計の x へ写す必要があり、#18 もこの写しの中で扱う。

- **前の設計の終点に揃える** (#18): 表の終点 $x_{last}$ と**前の設計** (prev run の `wall_repr.json` の `domain[1]`) の出口 $x_{e,prev}$ を比べ、
  $0 \le x_{e,prev} - x_{last} \le \varepsilon_{f32}$ なら終点の x を $x_{e,prev}$ に置き換える (値はそのまま; 新しい点は足さない)。
  $\varepsilon_{f32} = \mathrm{ulp}_{f32}(x_{e,prev} S_{prev})/S_{prev} + \mathrm{ulp}_{f32}(x_{e,prev})$ (格子の座標 [m] と抽出の r_t 単位の 2 回の丸め)。
  それより大きい差・超過は例外。始点も同じ規則で前の設計の入口端に揃える。`SingleBSplinePhysicalWall.applicability` は緩めない。
- **区間で正規化した座標で写す**: 壁温分布 (plan tooling-nozzle-wall-temperature-profile §4.1) と同じ $s$ (直管 [−2,−1]・縮流部 [−1,0]・拡大部 [0,1]、
  区切りの点はそれぞれの設計の壁から) で、前の設計の x → s → 新しい設計の x。端は端に写るので覆い漏れが起きない。
  実寸の長さを保つ作り直しでは、s の写しと実寸の x の写しの差は長さの公差 (1 mm) の程度 (δ_E の傾きで出口 δ に 1e-3 r_t 程度)。
- 大きさは今の `solve_rt_throat --prev` と同じ $\delta(x; r_t) = \delta_E\,(r_t/S_{prev})^{-0.2}$。

### 4.3 M5 デモを閉じる経路 (事前登録)

1. **Euler 参照 (前の設計)**: `problem_m5_air_c2b_euler.yaml` (r_t 0.024973561 m、run_0010 と同じ設計) で Euler (`run_0012_euler_c2b`、run_0003 と同じ段取り・本段 36000 step)。
2. **抽出**: run_0010 × run_0012 を修正後の抽出 (§4.1) に通す。合否は §6 V3。
3. **カードから作り直す**: 要求カード (c2b と同じ実寸の値) の r_t の固定点で、内側の寸法を CFD 前の積分法の代わりに 2. の δ_E (§4.2 の写しと換算) で解く
   (`wt_card` に NS の補正の入力を足す: 前の run・`delta_r_next.csv`)。長さは実寸で保ち、寸法ゲートで確かめる。最終の r_t での補正の表を YAML の横に書き、
   YAML の記録 (`_card`) に出典 (前の run・表・sha256) を残す。NS の準備は `prepare_ns(delta_r_csv = その表, offset = radial)` (CLI を足す)。
4. **Euler 参照 (新しい設計)**: 新しい設計の Euler (`run_0013_euler_<slug>`)。出口コア M の評価と流量の比、次の抽出の参照に使う。
5. **最終 NS**: 新しい設計 (`run_0014_ns_<slug>`)、IC = run_0010 の場の cross-mesh、`--stages none`・cfl 1・60000 step・`interfaceDiag: 1` (run_0010 と同じ)。
6. 合否は §6 V5。加えて 1 回の作り直しで固定点に近いか (§6 V5b) を記録する。FAIL なら再実行や係数の調整に進まず諮問する (AGENTS.md 条件 3)。

## 5. 実装ステップ

1. §4.1 (抽出) と §6 V1・V2。
2. §4.2 (#18) と §6 V4。
3. §4.3 の 1〜4 (case/70)。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | plan 段の諮問とレビュー | 諮問は済み (§6.1: 機構は採用、Major 2 件を採用して §3・§6 V2/V3・§4.3 を改定)。plan 段のレビューは #2 の結果で §4.3 を確定してから | F |
| 2 | ~~寸法の A/B~~ **完了 2026-10-13: r_t 単位の据え置きを棄却** (B ×1.004091・拡大部長 +4.42 mm) (codex 諮問の推奨) | forge 0 step。同じ問題 YAML (`problem_m5_air_c2b.yaml`)・同じ物理スロート径 (25 mm)・同じ壁の構築で、補正だけを A = run_0010 の壁の補正 (CONTUR k_f 0.8246、`solve_rt_throat` の CFD 前の経路)・B = run_0010 の zone の腕の δ_E 全分布 (`--prev` の経路; #18 は計測のため表の終点を x_e に揃えた写しで回避) にして寸法を解く。見る量: r_t、物理スロートの位置、実寸の長さ (直管・縮流部・スロート→出口) とカードとの差、Euler 参照 (run_0003) に対する尺度の 2 乗比。A は run_0010 の r_t (0.024973561 m) を再現すること (しなければ比較不成立)。判定 (事前): **A が公差内・B が公差外なら r_t 単位の設計を据え置く経路を棄却** (カードから作り直す)、**B も公差内なら寸法の公差破綻の仮説を棄却** | O |
| 3 | 等価半径の探索 | §4.1。触る: `design/forge_design/metrics/deltastar.py`・`design/tests/run_deltastar_tests.py`。合格: §6 V1・V2 | O |
| 4 | 表の端と写し、カードの NS 補正 | §4.2・§4.3 3。触る: `design/forge_design/feedback/deltastar_loop.py`・`intake/wt_card.py`・試験。合格: §6 V4 (写しの試験を含む) | O |
| 5 | Euler 参照 (前の設計) と受入れ | §4.3 1。`case/70/run_0012_euler_c2b`。**受入れ (codex 諮問 2026-10-13 ns-rebuild-route の採用)**: `check_convergence --segment` PASS かつ `qs_series.py` (出口コア M・無次元の流量の系列を `check_quasisteady --series-csv` で末尾 5 枚) が STEADY。**判別 A/B** (変えるのは Euler の本段の反復数だけ): A = run_0012 (36000 step)、B = run_0012 を複製して同じ設定で 12000 step 継続 (`run_0012b_*`)。両方から run_0010 を同じ修正済みの抽出に通し、事前の識別閾値 |ΔM_exit core| ≤ 2.5e-4・流量の相対差 ≤ 1e-4・max|Δδ_E| ≤ 1e-4 r_t。A・B とも PASS/STEADY で差が閾値内 → 参照の過渡が主因の仮説を棄却し A を使う。A が不合格か延長で閾値を超える → 36000 step の参照をそのまま使う方針を棄却 (諮問)。B も不合格なら判別未了 (過渡どうしの差を補正の効果と解釈しない) | O |
| 5b | ~~帯位置の固定の A/B~~ **完了 2026-10-13: 帯選択の寄与を支持** (x = 5 の 50k→55k: 再選択 5.1e-4・固定 2.4e-5 r_t ≤ 1e-4)。帯を固定しても δ_E は通算 120k まで直線的に動く (スロート +2e-6 r_t / 1k step; `run_0010b_ns_c2b_cont60k`)。Euler A/B の差も帯位置の差を伴う (y_b 最大 27 %)。**帯位置の固定の A/B** (諮問 2026-10-13 reference-acceptance) | forge 0 step。Euler は run_0012 に固定、run_0010 の 40k・45k・50k・55k・60k を再抽出。腕 A = 各時刻で edge を再選択 (`_deltaE_series_run_0010/`)、腕 B = 50k の最終の帯位置 y_b(x) を全時刻で固定 (流束のフィット・積分・平滑化は各時刻)。見る量: x = 5 の 50k→55k の差、使用域の最大差、y_star・y_b・hard/held の変化。判定: A で段差が再現し B で x = 5 の差 ≤ 1e-4 r_t → 帯選択の寄与を支持; B でも > 1e-4 r_t → 帯選択だけで説明する仮説を棄却 (流れの過渡とは断定しない)。B は診断用で、その表を設計入力に使わない | O |
| 6 | 抽出と作り直し | §4.3 2〜3。合格: §6 V3、寸法ゲート PASS (長さ 1 mm・径 0.02 mm) | O |
| 7 | Euler (新しい設計) と最終 NS | §4.3 4〜5。`case/70/run_0013_euler_<slug>`・`run_0014_ns_<slug>`。合格: run_0013 は #5 と同じ受入れ (収束 PASS・出口コア M と流量が STEADY)、最終 NS は §6 V5 (V5b は記録) | O |
| 8 | result 段 | 諮問とレビュー | F |
| 9 | (別セッションへ) M6 生産の抽出の走査 | 同じ症状 (δ_r = r_w) が M6 生産の抽出 (AWS) にあるかを M6 のセッションに伝える (調整板と SendMessage) | O |

## 6. 検証 (事前登録)

- **V1 単体** (合成の断面): 帯の傾きを大きくして $q_{ref}$ がコアで負になり、全半径の $F$ の末尾が $D$ を下回る断面を作る。
  旧実装は δ_r = r_w を返し (この試験が欠陥を検出できることの確認)、新実装は解析的に作った真の δ_r と相対 1e-6。区間の中で $F$ が単調でない断面は
  `nonmono_zone` (hard) で δ_r が NaN。既存の `run_deltastar_tests.py` が全件合格。
- **V2 後方互換**: 欠陥に当たる呼び出し (全半径の F の末尾 < D、梯子と感度の呼び出しを含む) が **0 件**の抽出で、新旧の
  `delta_r_equiv.csv`・`delta_r_next.csv` が**バイト同一**。対象は手元の case/45 の NS/Euler の対から当たり件数 0 のものを 1 つ (件数を数える計測フックで選ぶ)。
  当たり件数が 0 でない抽出 (case/70 run_0006: 598 件) は差を記録する (合否なし)。**[2026-10-13 改定]** 当初の「δ_r = r_w の断面が無い抽出」は梯子の汚染を見落とす。
- **V3 抽出の妥当性** (case/70 run_0010 × 新しい Euler 参照): (a) δ_r = r_w の断面が 0、(b) $0 \le x \le 3$ の $|\delta_E|$ の最大 ≤ 0.01 r_t、
  (c) スロートの δ_E と流量からの実効値 δ_eff,ṁ の差 ≤ 0.2·|δ_eff,ṁ| + 3e-4 r_t、(d) x ≥ 20 の δ_E が orig の腕と相対 1 % 以内
  (**[2026-10-13 改定、codex 諮問 Major]** 当初の x ≥ 10 は梯子の汚染の届く範囲 [x ≈ 3 + 半幅 3 + ノット 2 つ 12 ≈ 18] に入っていた。x = 10 では −9.3 % で、旧の腕が汚染側)。
  **既存の A/B (`_ab_root_zone/`) は誤抽出が消えることの証拠であって V3 の合格ではない** (V3 は新しい Euler 参照で回し直す)。
  **注記**: (c) の許容は run_0009/run_0010 の結果 (差 5 %・1.5 %) を見た後に書いた。根拠は事前の見積もり: 2 つの量は定義が違い
  (帯局所の欠損 と 断面全体の流量の比)、コアの流量係数の差 (`cd_rel` 1.006) と流量の測定の広がり (`mdot_ns_spread` 2e-5 相対) を含む。
  0.2 は「符号と桁が合う」ことの確認として置く。
- **V3e 入力 δ_E の準定常** (codex plan レビュー M4 と諮問 reference-acceptance): 受入れた Euler を固定し、NS の事前固定の 5 枚・20k step の窓 (5k 間隔) から
  同じ抽出・平滑化で δ_E(x) を作る。`check_quasisteady --series-csv` (窓全体) の VERDICT を残し、使用域 [0, x_e] で max_x(max_t δ_E − min_t δ_E) ≤ 1e-4 r_t・hard/held なし。
  単調な量は推定の残りの変化も同じ許容、推定不能は保留。合格した窓の**点ごとの平均**を入力 δ_E にする (不合格を平均で消さない)。
  **2026-10-13 の run_0010 の 40k〜60k の窓は不合格** (x = 2〜12 で 50k→55k に 5e-4 r_t の段差、x = 48 が単調に 2.2e-4 r_t 減る)。
- **V4 #18 と写し**: (a) 終点が $\varepsilon_{f32}$ 以内で前の設計の出口より短い表 → 揃えて物理壁が作れる、(b) $\varepsilon_{f32}$ を超えて短い → 例外、
  (c) 前の設計の記録 (`wall_repr.json`) が無い・表の範囲と合わない → 例外、(d) 同じ設計への写しは恒等 (相対 1e-15)、
  (e) 長さを変えた設計への写しで区切りの点 (−2・−1・0・1) が新しい設計の区切りの点に写る (相対 1e-12)、(f) 既存の `solve_rt_throat --prev` の経路 (同じ設計) の出力が今と同じ。
- **V5b 固定点の近さ** (記録、合否なし): run_0014 × run_0013 の修正後の抽出の δ_E と、作り直しに使った補正の表の差 (スロート・x = 10・出口)。
  1 回で固定点に近いか (M6 の実績では比 0.995〜1.010) の確認で、2 回目に進むかは諮問で決める。
- **V5 M5 デモ (最終 NS)**: 出口コア M (`nozzle_report` の `exit_core.M_mean`、η 0.05〜0.7) が **4.9975 ≤ M ≤ 5.0025**。
  出口コア M の準定常は**専用の許容** (codex plan レビュー M5): 末尾の変化と推定の残りの変化がともに ≤ 2.5e-4、単調なら漸近値も合否範囲内、推定不能は合格にしない
  (既定の `check_quasisteady` の drift 5 %・osc 10 % は合否の半幅 0.0025 に対して緩すぎる)。比較の対照は「各設計の受入れた定常解どうし」(run_0010 の IC は run_0009 で、同じ IC の A/B ではない; m6)。
  加えて `check_convergence --segment` の VERDICT、`check_quasisteady` で出口コア M・ṁ_NS/ṁ_Euler・最大 M・最大 p が STEADY、
  壁解像 (`check_wall_resolution.py`) PASS、メッシュ品質 PASS、成果物の寸法ゲート PASS。収束が PASS でなければ M の合否を出さない。
  比較の対照は run_0010 (同じ IC・段取り・step 数で、壁と r_t だけが違う)。出口コア M は新しい設計の Euler (run_0013) と同じ評価で出す。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | 2026-10-13 | [2026-10-13-tooling-nozzle-deltastar-band-root-zone-plan.md](../../notes/reviews/2026-10-13-tooling-nozzle-deltastar-band-root-zone-plan.md) | GO-with-changes, C0/M5/m2 | **全件採用** (採否は諮問 [2026-10-13-deltastar-reference-acceptance-diagnose.md](../../notes/reviews/2026-10-13-deltastar-reference-acceptance-diagnose.md) で確認): M1 使用域 [0, x_e] の hard/held で NS 後の寸法を止める + 抽出→平滑化→CSV→寸法の伝播試験 (§4.1・§6 V1) / M2 符号分岐の前に q_ref の有限・正と探索域を検査 (§4.1・V1 の負例) / M3 表の端は両側 \|Δ\| ≤ ε_f32 (§4.2・V4) / M4 入力 δ_E の準定常 (§6 V3e) / M5 出口 M 専用の準定常の許容 (§6 V5) / m6 V5 は各設計の受入れた定常解どうし / m7 methods を実装前に (済) と §5・§7 の整備 |
| diagnose (参照と δ_E の受入れ) | 2026-10-13 | [2026-10-13-deltastar-reference-acceptance-diagnose.md](../../notes/reviews/2026-10-13-deltastar-reference-acceptance-diagnose.md) | Euler A (run_0012)・B (run_0012b) とも受入れない (収束 PASS なし; 使用域でも x = 10 で 1.47e-4 r_t > 1e-4)。入力 δ_E の 40k〜60k の窓は不合格 (x = 48 が単調に 2.2e-4 r_t 減る、x = 3 DRIFTING)。Major 3 (呼び出し側の前提への異議): 使用域への限定で合格にしない・x < 0 を壁が使わないことから設計への無影響は導けない・plateau = float32 の床は未立証 | **採用**: 閾値と受入れ条件は事前登録のまま維持 (§5.1 #5・§6 V3e)。最終 NS を保留し、帯位置を毎回選ぶ / 50k に固定する抽出の A/B (§5.1 #5b) を先に行う。やらないこと: 使用域への限定で合格にする・末尾 2 枚の平均・原因の確認前に CFL と緩和を変える |
| diagnose (経路の改定) | 2026-10-13 | [2026-10-13-deltastar-ns-rebuild-route-diagnose.md](../../notes/reviews/2026-10-13-deltastar-ns-rebuild-route-diagnose.md) | カードから実寸を保って作り直す判断は採用。Major 1: Euler 参照の受入れが「完走・NaN なし」で弱い (未収束の参照の過渡が δ_E に入り得る; 確度 中) | **採用**: §5.1 #5・#7 に受入れ (収束 PASS + 出口コア M・流量が STEADY) と反復数だけを変える判別 A/B (閾値 2.5e-4・1e-4・1e-4 r_t) を入れた。不足情報の run_0010 の対象量の準定常は `qs_series.py` で補った: 出口コア M・無次元の流量とも STEADY (末尾 5 枚、M 4.990617〜4.990682) |
| diagnose | 2026-10-13 | [2026-10-13-deltastar-band-root-zone-diagnose.md](../../notes/reviews/2026-10-13-deltastar-band-root-zone-diagnose.md) | 機構 (全半径の F の `np.interp`) は採用・確度高。Major 2: (1) x = 10 で δ_E −9.27 %、V3(d) は満たさない (2) 寸法の変化は ≈1e-4 でなく約 0.4 % の見積もり | **全件採用**: (1) §3 を訂正、V3(d) を x ≥ 20 に改定し既存の A/B を V3 の合格扱いしない、V2 を「当たり件数 0」に改定 (2) §4.3 の「≈1e-4」を撤回、§5.1 #2 に寸法の A/B を入れて経路を決める。不足情報の指摘 (run_0010 の出口コア M・流量比そのものの準定常の VERDICT が無い、Euler run_0003 は未収束) は §6 V5 の判定に出口コア M と流量比の `check_quasisteady` を含めることで対応 |

## 7. 影響範囲

- `metrics/deltastar.py` の `band_local_deficit` (全ての帯局所抽出; 欠陥の無い断面はビット同一)、`feedback/deltastar_loop.py` の `solve_rt_throat --prev`。
- M6 生産 (case/45) の抽出も同じ関数を通る。過去の抽出の結果は変えない (再抽出は M6 のセッションの判断)。

## 8. 完了条件

- §6 V1〜V5 が合格、plan 段・result 段の codex レビューの行、methods の同期。

## 9. 変更ログ

- 2026-10-13: 起票。症状 (x ≈ 1〜3 の δ_r ≈ r_w) の機構を `np.interp` の全半径の探索と特定し、探索を積分域に限る腕で消えることと、
  スロートの δ_E が流量からの実効値と 1.5〜5 % で一致することを測った (§3、`case/70.m5_air_wt_demo/ab_root_zone.py`)。
  deliverables plan §5.1 #18・#19 をここへ移す。
- 2026-10-13: codex 諮問 (diagnose) を全件採用。§3 の「x ≥ 10 はほぼ同じ」を訂正 (x = 10 で −9.3 %)、梯子の汚染 (run_0006) を追記、V2・V3(d) を改定、§4.3 の寸法の見積もりを撤回し §5.1 #2 に寸法の A/B を入れた。
- 2026-10-13: §5.1 #2 の寸法の A/B (`case/70.m5_air_wt_demo/sizing_ab.py`、`_sizing_ab/`): A は run_0010 の r_t を再現 (0.024973561 m)、B (修正後の δ_E) は r_t 0.025075725 m (×1.004091)・拡大部長 +4.422 mm・直管 +0.818 mm・縮流部 +1.227 mm で公差 (1 mm) の外 → **事前の判定どおり r_t 単位の設計を据え置く経路を棄却**。§4.2 (写し) と §4.3 (カードから作り直す) を改定。codex の近似の予測 (×1.00403・+4.8 mm) と一致。
- 2026-10-13: 経路の改定の諮問を採用 (Euler 参照の受入れゲートと反復数の A/B)。`case/70.m5_air_wt_demo/qs_series.py` を追加し、run_0010 の出口コア M・流量が STEADY を確認。
- 2026-10-13: plan 段のレビュー (GO-with-changes, M5/m2) と参照の受入れの諮問を採用。Euler A/B は受入れ不可 (両腕とも収束 PASS なし、使用域でも 1.47e-4 r_t)、入力 δ_E の窓も不合格。最終 NS を保留し §5.1 #5b (帯位置の固定の A/B) を追加。methods を実装前に更新した (§4.1 の決定分)。
- 2026-10-13: §5.1 #5b — 段差は帯位置の切り替わり (固定で 2.4e-5 r_t)。帯を固定した δ_E は run_0010 を 120k まで延ばしても直線的に動く (冷却境界層が未整定; 総熱量 TRANSIENT-UNSETTLED と整合)。出口コア M は 4.99068 → 4.99077。
