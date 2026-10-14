# 諮問ブリーフ: 370 K・H の最初の固定壁の pass (V-h370) の事前登録と、#36a〜#36c で出た判断待ちの 5 点 (2026-10-14)

- 依頼者: 主セッション `M6ノズルの設計` (エスカレーション条件 1: §6 V-h370 を新規に書く、条件 3: 事前に決めた許容との比較が FAIL、条件 4: 新しい run の前)
- plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` §4.9・§5.1 #36・#36a〜#36d・§6
- 前の諮問: `notes/reviews/2026-10-14-half-inlet-370k-pass-inputs-diagnose.md` (NS 1 本 + H の Euler 参照 1 本 + 抽出の参照の A/B、停止規則と受け渡しの制御を先に、更新に進む条件と設計を終える条件を分ける)
- ユーザ決定 (2026-10-14): H を使う、370 K は B0 で進める、H の格子は Gref から。
- 実装: commit `48951bf5` (feature/nozzle-wall-fit-and-pipeline、push 済み)。forge・AWS はまだ回していない。

## 観測事実

### #36a 停止規則と受け渡しの関門 (実装済み、試験 ALL PASS)

- `design/forge_design/evaluate/stop_rule.py` は origin/feature/nozzle-cooled-wall-settling の 60de8125 (θ_r の断面の基準を設計壁の出口 `physical_wall.domain[1]` に直した版) を写した。forge-25 の 300 K の H の結果 (719dfbf4 = 66364fff と同一) は旧い基準 (x_E、case/45 で 15.9・27.9・37.8 r_t) で出たもので、断面が違う。
- 有効化は `evaluate.stop_rule: true`。断面は `stop_rule_thresholds.x_fracs` (x_w = 95.2974 に対する比; 40・70・94 r_t → 0.419739・0.734543・0.986386)。継続 `stop_rule_hold_steps` (最初の到達から H step のあいだ全出力で条件を保てば REACHED。途中で外れたら `hold.breaks` に記録して次の到達から数え直す。上限までに保ち切れなければ CENSORED)。
- 状態は `STOP_RULE_STATE.json` (REACHED / CENSORED / DIVERGED / EXEC_ERROR / DATA_ERROR / UNAVAILABLE / NOT_USED)。要求したのに使えない設定は run dir を作る前か S1 の前に例外。
- 受け渡しの関門 (`deltastar_loop.handoff_gate`): **状態が REACHED かつ `check_convergence.py --segment` が PASS のときだけ**抽出へ進む。上書きは `--allow-not-reached` (NOT_USED・CENSORED・UNAVAILABLE) と `--allow-unconverged` (NOT_CONVERGED)。DIVERGED・EXEC_ERROR・DATA_ERROR・記録なし・判定不能は上書きしても通さない。抽出は `band_select="edge"`。
- **未決 (実装者の指摘)**: 収束の判定は既定の `--drop 3`。B0 の本段の全域の残差は横ばいになる (H・300 K で 8.6e−8、core-grid plan §4.101) ので、REACHED でも NOT CONVERGED になって止まる見込み。
- 関門が見るのは停止の状態と本段の全残差だけ。§4.9 の他の条件 (保存則・局所 y₁⁺・抽出の被覆・δ_E と δ_target の隣り合う 2 窓・held) は親が別に確かめる。
- UNAVAILABLE (見張りの途中で帯の位置が決まらない) なら本段は fallback_steps (= nStepOuter = 200000) まで回る。

### #36b H の Euler 参照 (経路は実装、2 か所で止まった)

- 経路: `_physical_upstream` の検査 → CFD ピンの照合 (凍結の 16 項目) → `euler_reference_wall`。壁は x ≥ 0 が設計壁 S (+0、ビット同一)、x < 0 は直管 3.2425 と (−9, 3.2425, 0, 0) から x = 0 の S の (r, r′, r″) へつなぐ 5 次 Hermite。格子は `mesh_euler` (2000 × 97・全域 0.005・throat_refine 4・width 3)。NS の壁 − Euler の壁は x ≥ 0 で δ_r (差 9e−16)、x < 0 では x = 0 の δ_r の (値, 1 階, 2 階) に C² でつながる 5 次 Hermite (旧い参照の −248 mm の形の差は入らない)。
- **穴 1 (δ_r = 0 の壁が物理壁のゲートで落ちる)**: 凍結した S は出口端で r′ がわずかに負 (最小 −4.0e−8 at x 95.02、負の区間 x 94.90〜95.13、累積の減少 6.1e−9 r_t = 0.47 nm)。設計壁のゲート `monotone_gate` (`design/forge_design/geometry/wall_axismach.py:55`、閾値 :51-52) は合格。物理壁の `PhysicalNozzleWall._init_poly` → `wall_global_min` のゲート「スロートの後で r′ ≥ −1e−12」(`THROAT_R1_TOL`、同 :965) で不合格。δ_r > 0 の NS の壁では δ_r′ > 0 が隠すので出ない。
  - 案 A: Euler の参照を「H の上流を付けた設計壁」として作り、設計壁のゲート (`AxisMachCFDWall.validate` と同じ) で判定する。共有のゲートは変えない。壁ファイルは別の記録。
  - 案 B: `_init_poly` で δ_r が恒等的に 0 のときだけ「スロートの後」を S の `monotone_gate` で判定する。単一 B-spline のゲートにも同じ扱いが要るので表現は legacy に固定。
  - 案 C: S を直す → 凍結した設計を変えるので不可。
- **穴 2 (δ_in の照合が事前の許容を超えた)**: 許容 `DELTA_IN_TOL_MM = 1e−3` は実装者が測る前に決めた値 (前の諮問の時間変動の予算 0.01 mm の 1/10)。CFD 0 step・生成の座標で、NS (Gref の H370 の格子) と Euler 参照 (上の壁・`mesh_euler`) を照合:
  - 被覆: 全列 OK (旧い参照は終端 95.2276 で、NS の最後の 2 列 x 95.236〜95.297 が外れる)。
  - x ≥ 0 の max|δ_in − δ_r| = **1.137e−3 mm** (x 0.0107)。許容を超えるのは x 0.009〜0.151 の 24 列。区間ごとの最大は x 1〜40 で 3.6e−4 mm、試験部 1.2e−4 mm、94〜出口 1.2e−7 mm。
  - 内訳: 誤差はすべて抽出器が Euler の壁を列の間で線形補間する分 (S″·h²/8、スロートの列の間隔 0.0155 r_t)。NS の節点と真の間隔の差は 8e−13 mm 以下。
  - δ* ループへの影響 (式からの見積もり、ω = 0.5): `d_next = (1−ω)·d_in + ω·d_ext` なので固定点の δ_r は δ_E + e、残差 |δ_E − δ_in| は 2e (約 2.3 µm) から下がらない (設計を終える条件 0.05 mm に対して小さい)。
  - 案: 許容を決め直す (事後)、抽出器で Euler の壁を厳密に評価する (`metrics/deltastar.py` の変更、全チェーンに効く)、Euler のスロートの配置を細かくする (見積もり: ni 2200 で 0.95 µm、throat_refine 5 で 0.80 µm、ni 3000 で 0.51 µm)。

### #36c 370 K・H・B0 の問題 YAML (完了、CFD 0 step)

- `case/45.isobutane_m6_d155/problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml`: #34 の 300 K の H から Tw 300 → 370 と evaluate (implicit_relax 0.7・main_line B0・nStepOuter 200000・outStepInterval 2500、V-g0c-chain・C50 と同じ) だけが違う。**`evaluate.stop_rule` 系のキーはまだ書いていない**。
- δ_r は `_band_ab/half_inlet/delta_r_initial_C.csv` (sha 24251cb6、断熱の積分法で作り直した表) を prepare_ns の引数で渡す。
- 準備の再現 (変換の前まで): CFD ピンの照合 16 項目、格子 4719 × 160 (cap 0.015) は #34 の H とビット同一、生成の座標で AR 4574・skew 0.320 (上限 5000)、第一層 26.7 nm〜11.8 µm、壁 `wall_isothermal` Ts 370.0。S1 (cfl 0.5・1 次・3,000)・S2 (1.0・1 次・3,000) に B0 のキーなし、本段 cfl 4・2 次・B0・200,000・2,500、implicitRelax 0.7 は 3 段とも。
- **リミッタの基準値** (`main_line.limiter_ref`: ρ 0.8739869154、P 387491.0573、a 359.7712185。本段の `space` の先頭に入る基準値) は C の値のまま。準一次元の面積加重の見積もりで、C は ρ 0.872・P 387,957 (登録値と 0.3 % 以内)、H では ρ −35 %・P −37 %・a −3 %。
- Euler の YAML (`problem_d155_euler_h_ref.yaml`) の evaluate は v5d_M のまま (12000 step、cfl_main 2)。

### 参考: forge-25 の H・300 K (core-grid plan §4.104・§4.105、こちらでは未確認)

H (Gref、300 K) を B0 (`run_0917_hs_b0`) と point (`run_0918_hs_point`) で落ち着くまで回し、8 量の差はすべて 0.002 % 以下。B0 は本段の通算 12.25 万 step 以内に停止規則 (旧い断面) に到達、M6 の断面 (40・70・94、|欠損| ≤ 0.1 kg/s) でも同じ時点で到達。手順は「到達の後も 4 万 step 続けて規則を保ち、量の資格 (d ≤ τ/8・U ≤ τ/4・STEADY) を確かめる」。V-tw370 (今の形 C・G0c、370 K) は R0・C1 が 200,000 step で CENSORED (結果 B)。

## 期待値と出典

- 更新に進む条件と設計を終える条件 (前の諮問、plan §4.9): 更新へは REACHED・NaN/Inf なし・同一設定区間の全残差判定・保存則・局所 y₁⁺・抽出の全域被覆・δ_E と δ_target の隣り合う 2 つの 2 万 step 窓で全分布の変動 ≤ 0.01 mm・準定常・held/欠測なし。設計を終えるには加えて未緩和 max|δ_E − δ_in| ≤ 0.05 mm、物理壁の pass 間変化 ≤ 0.05 mm、出口 M・流量・寸法、真の最小断面積。0.05 mm・0.01 mm は検証済みの精度でなく暫定の予算。

## 問い

1. **穴 1**: 案 A・B のどちらか (または別案)。
2. **穴 2**: 事前の許容 1e−3 mm で FAIL した。許容を決め直すのは事後の緩和に当たるので、どう扱うか (抽出器を厳密化する・Euler を細かくする・この誤差を固定点のずれ e として予算に入れる、など)。最初の固定壁の pass (壁は更新しない) の前に直す必要があるか。
3. **リミッタの基準値**: H で C の値のままでよいか、H の値に直すか。直すなら何で決めるか (前の pass の場の平均、準一次元の見積もり、など)。B0 の設定の一部なので、変えると B0 の既往の結果 (forge-25 の H・300 K を含む) との対応が崩れる。
4. **受け渡しの収束判定**: B0 の本段の全残差は横ばいになる。関門の `check_convergence.py --segment` の判定の基準を何にするか (`--from-floor` の参照、落ち着いた区間の傾きの有意性、`--allow-unconverged` を理由つきで常用、など)。事前登録に書く形で。
5. **継続の意味**: 「外れたら数え直す」でよいか、「外れたら不合格」にするか。
6. **point の腕**: 370 K でも point を 1 本足すか (前の諮問の「自己の落ち着きだけでは共通の偏りを検出できない」への補い。forge-25 の 300 K の H では B0 と point の差 ≤ 0.002 %)。足すなら何を比べるか。
7. **V-h370 の事前登録の形**: run の本数・step の上限・保存間隔・停止規則の設定 (x_fracs・hold 4 万)・判定の量と区間・Euler 参照の本段の上限 (前の諮問の提案は 6 万)・抽出の参照の A/B (腕 A 旧 run_0174 の共通区間、腕 B H の Euler 参照) の合否。最小に。

## 読んでよいファイル

- `plans/active/tooling-nozzle-isothermal-wall-chain.md` §4.9・§5.1 #33〜#36d・§6 (V-tw370 の結果)
- `design/forge_design/evaluate/stop_rule.py`、`design/forge_design/evaluate/runner_axismach.py` (`_stop_rule_settings`・`run_staged_ns`・`euler_reference_wall`・`delta_r_zero`・`_main_line_spec`)、`design/forge_design/feedback/deltastar_loop.py` (`handoff_gate`・`run_pass`)、`design/forge_design/metrics/delta_in_check.py`、`design/forge_design/metrics/deltastar.py`、`design/forge_design/geometry/wall_axismach.py` (:51-55・:965 付近)
- `case/45.isobutane_m6_d155/problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml`・`problem_d155_euler_h_ref.yaml`・`h370_problems.py`・`_band_ab/half_inlet/h370_prep_check.json`・`h370_problems.json`
- `notes/reviews/2026-10-14-half-inlet-370k-pass-inputs-diagnose.md`・`2026-10-14-half-inlet-370k-design-diagnose.md`
