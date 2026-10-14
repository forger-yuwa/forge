# 諮問: 熱伝導 Jacobian のキー 1 がノズルを不安定化させた理由の見立てと、対策 (b) キー 5 の妥当性

日付 2026-10-09。諮問先 codex (diagnose)。エスカレーション条件 4 (原因を書く前) と 6 (cuda_forge の数値の振る舞いの変更)。
plan: `plans/active/time_integration-implicit-thermal-jacobian.md` (§4.1・§4.3・§6・§6.0 を全文読むこと)。作業ツリー `/home/sano/work/forge-integ-1005` (commit 06d1b149)。
コード: `solver_density_cuda/cuda_forge/timeIntegration_d.cu` の `implicit_defect_correction_block_d` の粘性の対角 (ビット 1・4 の分岐、`thermalJac`)。
ユーザ: 「じゃあ b」(行 4 に従来のスカラーも残す案)。

## 観測事実 (§6.0 にも記載)

- キー 1 (行 4 を温度の熱伝導の Jacobian に置き換え、行 0〜3 は従来のスカラー A = Σ 2ν_eff δ/dcc):
  - point cfl 4 (run_0216、case/45 の 300 K 冷却ノズル、run_0183 の res_100000 から): 最初の 200 step は従来 (run_0191) より速く残差が下がる (rms_roe 10.4 → 4.5) が、
    300 step 前後から運動量の残差が育ち 922 step で NaN。非有限は x_w 69.9〜70.5 の壁から 79〜89 層目。キー 0 (run_0191) は同じ条件で 40000 step 安定。
  - その場所 (run_0208 の res_200000、キー 0 の収束に近い場、列 4307・壁から 88 層目): μ_t/μ 474 (列内の最大)、ρ 0.031 kg/m³、T 261 K、|u| 1710 m/s、c 324 m/s、
    e = ρE/ρ − ½|u|² = −1.16e5 J/kg (TP の顕熱の基準)、½|u|² = 1.46e6 J/kg、γ 1.377、c_v 775。
    概算で A/(V/Δτ) = 0.71、B/(V/Δτ) = 0.54 (B = Σ k_f δ/dcc /(ρ c_v)、k_f は乱流の分を含む面の値)。
  - directional (cfl 4/2/1) + キー 1 は 1 層目の熱伝導の段の成長を止めたが、数百 step で対流の段が育って発散 (run_0212・0214・0215)。
- 純伝導 (U1、静止、一様格子、D = Δτκ/Δy² ≤ 0.5) ではキー 0 と 1 は同じ解に収束し、キー 1 は同じ cfl で 12〜40 % 遅い。cfl 50 は両方とも負の P・ρ で壊れる。

## 呼び出し側の見立て (仮説、確かめていない)

行 0〜3 は物理の微分でない人工的な減衰 A (スペクトル半径の近似) のまま、行 4 だけを温度の微分 B·ρc_v·ΔT に替えたので、行ごとに歩幅の縮め方の基準が食い違う。
温度を変えない残差 R = R₀(1, u, v, w, E) を対角ブロックだけで解くと (V/Δτ = 1 に規格化、対流 A⁺ は無視):
Δρ = R₀/(1 + A)、ΔρE は行 4 から ΔρE(1 + B) = E R₀ + B(e + ½|u|²)Δρ、よって ΔT = (E/c_v)(Δρ/ρ)·A/(1 + B)。
x = 70 の値で E/c_v ≈ 1740 K、A/(1 + B) ≈ 0.46 → Δρ/ρ 1 % ごとに約 8 K の偽の ΔT。壁際 (u ≈ 0) では E/c_v ≈ T で小さい。
NaN の場所 (μ_t 最大・高速) と一致する。

## 対策 (b) = キー 5 (ビット 1 + 4、実装済み 06d1b149)

行 4 = A·e₅ᵀ (従来のスカラー) + B·[−(e−½|u|²), −u, −v, −w, 1] (温度の項)。行 0〜3 は従来どおり。
全行に同じ A がかかるので、温度を変えない補正は従来どおり向きを保って縮み、温度が変わる補正にだけ熱伝導の見込みが加わる。
事前登録 (plan §6.0): Vb-point (point cfl 4 + キー 5、20000 step、早い判定 2000 step)、Vb-dir (directional cfl 4 + キー 5、200 step、帳簿)。

## 問い

1. 上の見立て (行ごとの減衰の基準の食い違いが、高速・高 μ_t の所で偽の ΔT を作る) は、観測 (point の NaN の場所・時期、directional の対流の段の成長) と整合するか。対流 A⁺ を入れると結論は変わるか。別の機構の候補は。
2. (b) は見立てに対する対策として妥当か。A と B の両方が行 4 にかかることで、温度のモードの減衰が二重になる影響 (過剰な減衰で遅くなるだけか、別の不整合を生むか)。
3. 事前登録の Vb-point・Vb-dir の判定に穴はないか。(b) が point で安定でも directional で対流の段が育つ場合、次に何を見るべきか。

## 読んでよいもの

- 上記 plan、`methods/time_integration/implementation.md` の「エネルギー行の熱伝導 Jacobian」節、`solver_density_cuda/cuda_forge/timeIntegration_d.cu`・`viscousFlux_d.cu`
- `plans/active/tooling-nozzle-isothermal-wall-chain.md` §5.1 #27 (経緯)、`case/45.isobutane_m6_d155/README.md` (run_0191・0203〜0216 の行)
- 過去の諮問: `notes/reviews/2026-10-09-linedir-divergence-2-diagnose.md`、`notes/reviews/2026-10-09-time_integration-implicit-thermal-jacobian-plan.md`
