# 諮問ブリーフ: NS の 3 条件 (N0・N1・N2) が出口コア M のゲートを 3 本そろって約 3e-4 下回った — 次の手

- 日付: 2026-10-07
- 関連 plan: `plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md` (§6 U4・§9 末尾)、`plans/active/discretization-moc-axis-limit-and-corrector.md` (§6 V5′・§9 末尾)、`plans/active/verification-case45-euler-total-enthalpy.md` (E4V・E5)
- エスカレーション条件: 3 (事前登録した合否の FAIL)・4 (登録に無い修正・再実行へ進む前)

## 1. 観測事実

評価器 `case/45.isobutane_m6_d155/ns_n012_eval.py` (sha256 f1129721…)、出力 `case/45.isobutane_m6_d155/_band_ab/ns_n012_eval.json`・`ns_n012_eval.log`。コード a6ce9390 (AWS、バイナリは run の RUN_PROVENANCE)。3 本とも RUN_RC 0、NaN なし、本段区間 `check_convergence --segment` は NOT CONVERGED (stalled/plateau、全不合格列が STALLED、RISING なし)。判定窓 60000〜80000 の 5 枚。

| 量 (窓の平均) | N0 run_0165 | N1 run_0166 | N2 run_0167 | ゲート | 旧生産 run_0147+0149 |
| --- | --- | --- | --- | --- | --- |
| 出口コア M | 5.998527 (未達) | 5.998513 (未達) | 5.998540 (未達) | 5.9988〜6.0012 | 5.998871 |
| 準定常 出口コア M | STEADY (単調減少、漸近値 5.99851) | STEADY | STEADY | STEADY | |
| オーバーシュート η0.1 [%] | 0.00238 | 0.00229 | −0.00136 | ≤ +0.035 | 0.00787 |
| 準定常 オーバーシュート | DRIFTING (drift 34 %、tail-end に極値) | OSCILLATING (0.00229 ± 0.00017) | DRIFTING (drift 17 %) | STEADY | |
| 波 η0.1 [%] | 0.00651 | 0.00664 | 0.00659 | ≤ 0.01 | 0.00642 |
| 準定常 波 | DRIFTING (drift 7.6 %) | STEADY | STEADY | STEADY | |
| δ_E/δ_C (x_F) | 1.00167 | 1.00166 | 1.00167 | 1 ± 0.005 | 0.99987 |
| 出口半径 [m] | 0.7749987 | 0.7749987 | 0.7749977 | 0.775 ± 0.0001 | |
| 壁解像 y1+ > 1 の面積 | 3.6 % | 3.6 % | 3.6 % | ≤ 5 % | |

- 差 (判定なし、窓の平均の差、括弧は両窓の幅の大きい方):
  - N1 − N0 (U4): 出口コア M −1.4e-5 (3.3e-5)、オーバーシュート −8.4e-5 % (9.2e-4)、波 +1.3e-4 % (4.5e-4)、δ_E/δ_C −9.1e-6 (2.8e-5)。
  - N2 − N1 (V5′): 出口コア M +2.7e-5 (2.1e-5)、オーバーシュート −3.66e-3 % (4.9e-4)、波 −5.2e-5 % (4.8e-4)、δ_E/δ_C +1.7e-5 (3.2e-5)。
  - N0 − 旧生産: 出口コア M −3.44e-4 (3.3e-5)、δ_E/δ_C +1.80e-3 (2.8e-5)、オーバーシュート −5.5e-3 %。
- 凝縮 (run_0168〜0170) は「dry が全ゲート合格のときだけ」の登録どおり作っていない。

## 2. 期待値と出典

- ゲート: 出口コア M 6.000 ± 0.02 % (verification-m6 §5.1 #9 の登録、`plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md:75`)。他は monotone plan の N の条件の転記 (upstream plan §6 U4)。
- 延長: 「未達の量があれば延長 1 回 (20000 step)。なお未達なら未達のままユーザに判断を仰ぐ」(upstream plan §6 U4)。
- Euler の出口較正 (E4V、run_0164): 新しい格子 (mesh_euler、2000 × 97・全域 0.005) で M_common 13 枚の平均 5.99999968。較正値 6.8825162455159465e-06 (旧値 +3.770e-4 は全温が異常な G1 の Euler 由来)。plan には「NS の出口の合否は NS で判定する (Euler の較正の合格を NS に移せるとは限らない)」と書いてある (euler plan §6 E4)。

## 3. 再現条件

- 3 条件で固定: 較正値 6.8825162455159465e-06、r_t 0.07666536551630307 m (N0 の問題で CFD 前の solve_rt で出口半径 0.775 m に解き直し、3 条件で共通)、k_f 1.054129117086371 (旧生産と同じ)、NS 格子 G1 (2000 × 97、壁際細分化)、数値設定は生産の投入スクリプトと同じ。IC は run_0149 の res_20000 から cross-mesh + 段階起動 (soft 3000・mid 3000・本段 80000)。
- N0 = ramp・legacy MOC、N1 = poly・legacy、N2 = poly・analytic + converge。Euler 参照 (δ_E の抽出): N0・N1 は run_0164、N2 は run_0174 (V5d の腕 M)。

## 4. 実施済みの操作と結果

- 出口較正のやり直し (E4 段 1 run_0163 は判別不能、E4V run_0164 合格)。
- r_t の解き直し (出口半径ゲートは 3 本とも合格)。
- 評価器のゲート判定 (上の表)。延長はまだ回していない。
- E5 (δ_E の参照場の感度: NS の最終場を固定し Euler 参照の窓 13 枚を差し替えた δ_E/δ_C の幅、合格 ≤ 5e-4) を AWS で実行中。結果が出たら追記する。

## 5. 仮説 (呼び出し側の読み、確定していない)

- H1: 出口コア M の未達は、Euler の出口較正を新しい格子でやり直したことで、旧生産が持っていたかさ上げ (G1 の Euler の誤差に合わせた +3.77e-4) が無くなった結果。NS の出口 M と Euler の出口 M の差は、旧生産で −1.13e-3 (5.99887 − 5.99999)、今回 −1.48e-3 (5.99852 − 5.99999968)。N0 − 旧生産 −3.44e-4 は較正値の変化 −3.70e-4 とほぼ同じ (応答係数 ≈ 0.93)。
- H2: NS の出口 M を 6 にするには、NS の結果から較正値を 1 回補正する (+1.48e-3 / 0.93 ≈ +1.6e-3) か、C2 (k_f・r_t を NS の δ_E から解く) をやり直す必要がある。較正値を +1.6e-3 動かすと出口半径が約 +0.5 mm 動く (+1.4e-5 で +4.7 µm の比例) ので、r_t の解き直しも要る。
- H3: オーバーシュートの準定常の未達は、平均が零に近い量に相対の drift の判定を当てた結果で、値そのものは上限の 1/15。E2 で同じ問題に絶対の幅の条件を足した前例がある (登録の後で条件を足すのは避けたい)。

## 6. 諮りたいこと

1. 登録どおり延長 (3 本 × 20000 step、GPU 約 15 分) を回すべきか。出口コア M は STEADY (単調減少の漸近値 5.99851) で、延長で下限に届かない見込み。オーバーシュートと波の準定常の未達には延長が効くかもしれない。
2. 出口コア M の未達をどう扱うか。候補: (a) NS の結果から較正値を 1 回補正し、r_t を解き直した N0′ を回す (登録の変更が要る)、(b) C2 で k_f・r_t を NS の δ_E から解き直す、(c) 未達のままユーザに判断を仰ぐ。設計の目標は「出口の M ≈ 6 厳密 (物理壁、粘性込み)」(case/45 の問題の冒頭のユーザ条件)。
3. U4 (上流の多項式化) と V5′ (MOC) の判断は、3 条件が同じ絶対ゲートで落ちている状態でどう進めるか。N1 − N0・N2 − N1 の差は小さい (N2 のオーバーシュートの低下を除き、時間の幅と同じ桁以下)。
4. H1〜H3 の読みに誤りや見落としがあれば指摘してほしい (特に H1 の「かさ上げ」の解釈)。

## 読んでよいファイル

- `plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md`、`plans/active/discretization-moc-axis-limit-and-corrector.md`、`plans/active/verification-case45-euler-total-enthalpy.md`、`plans/accepted/tooling-nozzle-cfd-pinned-initial-line.md`、`plans/accepted/tooling-nozzle-throat-monotone-r2.md`
- `case/45.isobutane_m6_d155/_band_ab/ns_n012_eval.json`・`ns_n012_eval.log`、`case/45.isobutane_m6_d155/ns_n012_eval.py`・`ns_n012.py`・`run_ns_n012.sh`、`case/45.isobutane_m6_d155/README.md`
