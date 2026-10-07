# 諮問: 壁温が変わる前提で、積分法 (CONTUR) の物性をどの温度・どのモデルで評価するか

日付 2026-10-08。諮問先 codex (diagnose、`~/.config/forge/diagnose-backend` = codex)。エスカレーション条件 1 (plan §4 の設計方針を新規に書く)。
作業ツリー `/home/sano/work/forge-integ-1005` (ブランチ `feature/nozzle-wall-fit-and-pipeline`)。
関連 plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` (等温壁チェーン、§4.6 壁温影響の評価方針)、
`plans/active/verification-m6-axis-wave-mesh-su2.md` §5.1 #8d と §9 2026-10-08 (CONTUR 較正の 4 係数化の前段)。

## ユーザの発言 (2026-10-08)

- 「これから壁温設定はガンガン変わり得るからね。気を付けてね」
- 「等温壁にするかもしれないし、壁温分布を与えるかも。CONTUR で物性をどの温度で評価するか？は任意性がある話かと思うので、どうするべきかは考えてよ」

## 観測事実

対象は case/45 (イソブタン燃焼ガス M6、semi-perfect、組成 Y = CO2 0.1677 / H2O 0.0858 / O2 0.0220 / N2 0.7245、Tt 1600 K、Pt 5.5 MPa)。
生産の NS = `case/45.isobutane_m6_d155/run_0167_ns_n012_N2` + `run_0179_ns_n012_N2_ext` (断熱壁、SST 低 Re、viscMethod 2 = 種ごとの CEA 輸送物性)。
CONTUR の実装 = `design/forge_design/feedback/deltastar_integral.py` (全文を読むこと)。測定スクリプト = `case/45.isobutane_m6_d155/delta_contur_compare.py` (`extract`/`plot`/`tw`/`taw`)。

1. **CONTUR の中で温度・物性が入る場所** (実装を読んで洗い出したもの。漏れがあれば指摘してほしい):
   - (a) 温度分布 Eq. 69 を**温度**で書いている: T(u) = T_w + a(T_aw − T_w)u + (T_e − a(T_aw − T_w) − T_w)u²、ρ/ρ_e = T_e/T。
   - (b) 断熱壁温 T_aw = T_e(1 + r(γ_e − 1)/2 M²)、r = Pr^{1/3}、Pr = 0.72 固定、γ_e は縁の局所値。
   - (c) 粘性 μ(T) は**空気の Sutherland** (1.716e-5, S 110.4) を μ_e・μ_w の両方に使う (`feedback/deltastar._sutherland`)。
   - (d) 圧縮性変換: F_c = [∫₀¹(ρ/ρ_e)^{1/2}du]^{-2} (分布から)、F_Rδ = μ_e/μ_w (壁温で評価)、C_f = k_f · C_fi(R_θi)/F_c (van Driest II 型)。
   - (e) N(R_δ) の表と R_θc は縁の μ_e。
   - (f) 入口の θ0 は平板の式 (縁の物性)。
   - (g) 縁の状態 (T_e, p_e, ρ_e, u_e, γ_e) は semi-perfect の気体 (NASA-9 の c_p(T)) から作っていて、ここは NS と同じ気体。
2. **T_aw の式が全温を超える**: 試験部 x = 40〜90 r_t で CONTUR の T_aw = 1654〜1663 K (Tt 1600 K)。NS の断熱壁温は 1470〜1482 K。
   - エンタルピーで回復させた h_aw = h_e + r(h_0 − h_e) (r = 0.72^{1/3} = 0.896、同じ semi-perfect の h(T)) は 1473〜1484 K で、NS との差は +1.5〜4 K。
   - NS の実効の回復係数 (エンタルピー基準) は 0.893〜0.895。
   - 混合気の Pr (下の 4) の 1/3 乗 r = 0.909 だと 1488〜1498 K で、NS から +16 K 離れる。
3. **T_aw の式の違いが δ_r を動かす量** (生産の k_f = 1.0541 のまま、`delta_contur_compare.py taw`、出力 `_band_ab/delta_contur/taw_sensitivity.json`):
   - 比較の式は全温基準 T_aw = T_e + r(T_t − T_e) (エンタルピー形の近似)。
   - 等温壁でも T_aw は Eq. 69 の (T_aw − T_w) の項に入るので差は消えない。
   - x = 40 / 70 / x_F (95.2) の δ_r の相対変化:
     - 断熱: −3.2 / −3.0 / −3.4 %
     - T_w 1000 K: −3.0 / −3.0 / −3.1 %
     - T_w 600 K: −4.2 / −4.1 / −4.2 %
     - T_w 300 K: −5.5 / −5.3 / −5.4 %
   - 出口の δ_r の壁温依存:
     - 現行の式: 断熱 56.17 mm → 300 K で 56.70 mm (+0.95 %、冷やすと厚くなる)。
     - 全温基準: 54.25 → 53.67 mm (−1.1 %)。
4. **粘性・Pr の食い違い** (NS と同じ CEA 輸送物性の独立参照実装 `solver_density_cuda/tests/unit/transport_reference.py` で計算):
   - μ_mix/μ_Sutherland(空気): 250 K で 0.956、350 K で 0.973、600 K で 1.007、1000 K で 1.050、1470 K で 1.089、1600 K で 1.099。
   - Pr_mix = μc_p/λ は 0.745〜0.758 (CONTUR は 0.72)。
   - 断熱 (T_w ≈ 1470 K) では F_Rδ = μ_e/μ_w が、空気の式だと混合気より約 11 % 大きい。T_w = 300 K では食い違いの向きが変わる。
5. **NS の壁温を CONTUR に与えた実験** (`delta_contur_compare.py tw`、NS の断熱壁温 `ns_wall_T.csv` を `Tw_table` で与える):
   - δ_E/δ_C (k_f = 1) は、出口で 1.047 → 1.056、試験部の振れは 2.93 → 3.06 %。δ_C はさらに薄くなる。
   - 出口に合わせた k_f は 1.0565 → 1.0682 (+1.1 %)。
   - 試験部の傾き (振れ 2.8〜3.1 %) は変わらない。δ の不足と傾きは T_aw の誤差では説明できない。
6. **NS と CONTUR の摩擦**: 試験部の平均で C_f(NS)/C_f(CONTUR, k_f = 1) = 0.979。δ_E/δ_C (k_f = 1) は 1.047〜1.08 (NS が厚い)。
7. **既存の冷却壁の NS データ**:
   - case/48 平板: M4.19、空気 CPG、Sutherland、Pr 0.72。`case/48.flat_plate_cooled_m4/README.md` と isothermal plan §5.1 #2 を参照。
     - C_f/VD-II: T_w 300 K で 0.97〜0.99、700 K で 0.91〜0.95、断熱で 0.87〜0.92 (断熱は Re_θ 1400〜3500 と低い)。
     - δ* の NS/CONTUR: 300 K で +3.5〜6 %、700 K で +8〜16 %。θ は ±0 / ±3 %。
   - case/44 ノズル: va3 vitiated air の TP、Tt 1060 K。`run_0108` (断熱) と `run_0109`/`run_0114` (300 K) は同じ壁・同じメッシュ。
     - 冷却で出口コア M は +0.35〜0.48 %。
     - 300 K の積分法初期壁 (現行の式) は冷却の効果を過小評価したが、NS の δ* 反復 1 回 (`run_0115`) でゲートに入った。
8. **配管の現状**:
   - `spec.wall_thermal` は `adiabatic | isothermal (一様 Tw)` だけ。これが単一ソースとして NS の bcond と CONTUR の `thermal_bc` の両方を作る (`design/forge_design/probdef.py`)。
   - CONTUR は `Tw_table` を受けられるが、問題 YAML からは渡せない。
   - ソルバは `wall_isothermal` + `ints: {wallProfile: 1}` で T_w(x) を受けられる (`methods/boundary.md`「壁温分布の入力」、node で実測済み)。
9. **較正の記録**:
   - 生産の k_f (`deltastar_initializer.cf_scale` = 1.054129117086371) は断熱壁の NS で出口に合わせた値。
   - YAML には、どの壁温条件で較正したかの記録が無い。prepare も壁温条件の違いを検査しない。

## 期待値と出典

- 断熱壁の回復は全温を超えない (r < 1 ならエンタルピーで h_aw < h_0)。超えるのは式の誤り。
- Crocco–Busemann / Walz の関係はエンタルピーの関係 (c_p 一定のとき温度で書ける)。
- van Driest II は冷却壁の C_f で最も良いとされる (Hopkins & Inouye 1971)。Eckert の参照エンタルピーも広く使われる。
  CONTUR (Sivells AEDC-TR-78-63) の F_c・F_Rδ は van Driest II 型。

## 実施済みの操作

上の 2〜6 は CFD 0 step の手元計算 (`design/.venv-opt`)。コードの既定は何も変えていない。生産の壁・run は変えていない。

## 仮説 (私の案)

物性の評価を 2 層に分ける。

- **第 1 層 (任意性の無い整合)**: NS と同じ気体で書き直す。CPG + Sutherland + Pr 0.72 の問題では、今とビット同一 (または丸めの差) になる。
  - (a) Eq. 69 をエンタルピーで書き、T(u) は気体の h(T) の逆で戻す。
  - (b) h_aw = h_e + r(h_0 − h_e)。r は 0.72^{1/3} のままにして、NS の回復 (0.893〜0.895) と照合する。混合気の Pr は使わない (NS とずれるため)。
  - (c) μ(T) を NS と同じ輸送モデルにする (gas.transport があれば混合気、無ければ Sutherland)。
  - 合格の目安: CONTUR の T_aw が NS の断熱壁温と ±0.5 % で合う (case/45・case/44)。
- **第 2 層 (本当に任意な閉包)**: 圧縮性変換をどの参照温度で取るか。候補は今の van Driest II 型 (F_c・F_Rδ = μ_e/μ_w) と Eckert の参照エンタルピー。
  - 好みで選ばず、壁温を跨いだ交差検証で選ぶ。
  - 較正 (k_f、必要なら k_N) を 1 つの壁温条件で決め、別の壁温条件の NS (case/48 の 3 条件、case/44 の断熱/300 K、case/45 の断熱) を予測させる。
  - δ* と C_f の最大誤差が小さい方を採る。基準は結果を見る前に登録する。
- **較正の記録とガード**:
  - `deltastar_initializer` に較正時の `wall_thermal` (と、第 1 層の版) を記録する。
  - prepare は、問題の `spec.wall_thermal` と違えば止める (`require_moc_gate` と同じ形)。
  - 壁温を変えたら「NS の δ* 反復 1 回で較正し直す」を既定の手順にする。第 2 層で壁温を跨いで当たることが示せたら緩める。
- **壁温分布**:
  - `spec.wall_thermal` に `mode: profile` (表) を足し、NS の `wallProfile` CSV と CONTUR の `Tw_table` の両方をそこから作る。
  - 座標は物理長 [m] (スロート原点) にする。r_t を解き直しても壁温の位置が動かないようにするため。

## 問い

1. 温度・物性が入る場所の洗い出し (1 の a〜g) に漏れはあるか。2 層に分けることは妥当か。第 1 層に入れたもののうち、実は任意性があるもの (r に何を使うか、Pr、Walz の a、ρ/ρ_e = T_e/T の組成一定の仮定) はあるか。
2. 第 1 層の r: NS と合う 0.72^{1/3} を採るのは、NS (SST の Pr_lam 0.72・Pr_t 0.9) に合わせる意味で正しいか。物理の正しさ (混合気の Pr) と NS との整合のどちらを優先すべきか。
3. 第 2 層の選び方 (交差検証、使えるデータ、事前登録する基準) は妥当か。case/48 は低 Re_θ・空気 CPG、case/44 は vitiated air の TP、case/45 は燃焼ガスで断熱だけ。この組で判定できるか。足りないなら、どの NS を 1〜2 本足すのが最小か (例: case/45 の生産形状で 300 K を 1 本)。
4. 較正のガード (壁温条件が違えば止める) は強すぎ / 弱すぎないか。壁温分布のときの「同じ条件」の判定はどうするか (表の一致? 代表値の許容差?)。
5. 第 1 層を入れると生産の壁 (今の CONTUR + k_f) と食い違う。生産は今のまま記録として残し、新しい設計から新しい式を使う (k_f は新しい式で較正し直す)、という切り替えで問題ないか。
6. #8d (CONTUR の 4 係数化) との順序。第 1 層を先に入れてから #8d の較正をするべきか (第 1 層で δ_C が 3〜5 % 動くので、先に較正すると無駄になる)。

## 読んでよいもの

- `design/forge_design/feedback/deltastar_integral.py`、`design/forge_design/feedback/deltastar.py` (`_sutherland`)、`design/forge_design/probdef.py` (`wall_thermal`)、`design/forge_design/evaluate/runner_axismach.py` (`prepare_ns` の `thermal_bc`、`require_moc_gate`)
- `case/45.isobutane_m6_d155/delta_contur_compare.py` と `_band_ab/delta_contur/*.json` (`summary.json`・`tw_experiment.json`・`taw_sensitivity.json`)
- `plans/active/tooling-nozzle-isothermal-wall-chain.md`、`plans/active/verification-m6-axis-wave-mesh-su2.md` (§5.1 #8d、§9 の 2026-10-08)
- `case/48.flat_plate_cooled_m4/README.md`、`case/44.*/README.md` の run_0108〜0117 の行
- `methods/design/overview.md` の「壁の熱境界条件」節、`methods/boundary.md` の「壁温分布の入力」
- `solver_density_cuda/tests/unit/transport_reference.py`
