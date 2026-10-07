# 諮問: Euler の全温の超過は保存状態にある (V5c) — MOC の判断と、この異常の扱いの次の手

日付 2026-10-07。諮問先 codex (diagnose)。エスカレーション条件 4 (plan に無い次の手へ進む前) と 7 (解釈を確定する前)。V5c の登録に「保存状態側の診断へ進む (次の手は改めて諮る)」とある。
plan: `plans/active/discretization-moc-axis-limit-and-corrector.md` (§6 V5・V5b・V5c と位置づけ、§9 の 2026-10-07 の記録)。作業ツリー `/home/sano/work/forge-integ-1005` (HEAD bddb2932)。

## 観測事実

- V5c の出力 `case/45.isobutane_m6_d155/_band_ab/moc_v5c_thermo_ab.json` (評価器 6a91da03…、登録 4a484270)。保存された P・T・h0 から求めた全温 (A) と、保存量から独立に復元した全温 (B) は、全 56 標本で閾値内 (P 傾き ≤ 5.9e-5 %pt、全温 ≤ 0.65 K)。Tt (1600 K) を 100 K 以上超える全温が 55 標本で両経路に再現 (最大 +373 K)。層ごとの記録 (壁からの j_w) は JSON の `samples[*].layers`。
- 超過の時間変化: 番号写像の IC の元 run_0114 (旧壁の Euler、収束場) で +241 K → 腕 B run_0143 +334 K、腕 M run_0150 +373 K (step 18000)。ISEN は step 18000 で +75 K → 延長で +343 K (通算 54000)。超過の節点は数千〜1 万個、スロート付近の壁際 (y ≈ スロート半径、x ≈ 0〜0.5 r_t が最大)、ISEN では x −0.96〜7.3 m に広がる。全温の中央値は 1599.998 K。
- 等エントロピー IC の初期場 (run_0153 step 0) にも Tt を超える節点が 95,983 個あり、最大 1915 K。
- 設定 (run_0150 の `solverConfig.yaml`): Euler (非粘性、`wall_dist` 1e30)、node、2 次 (convMethod 1)、limiter 2、SLAU、block-DPLUR (timeIntegration 11)、cfl 2・implicitRelax 0.7、nStepInner 5、lowMachPrecond 0、TP (化学種 2 つ: MIXDRY・H2O)、壁はすべり壁。soft 段は 1 次・cfl 0.5。
- 残差は腕 B・腕 M・ISEN とも本段・延長で stalled/plateau (0.1 dec)。
- MOC の変更の位置づけ (前回の諮問で採用済み): 候補のまま、生産採用と既定の切り替えは保留。V5・V5b は保留。
- 設計チェーンでの Euler の使われ方: CFD でピン止めする初期線は Euler run_0062 の場から取る (`geometry.initial_line_run`)。出口較正 (`Md_moc_offset`) は細分格子の Euler で決めた (run_0113+0114)。単調壁の採用 (E′) も Euler の比較。
- 生産の評価は NS (run_0147・0148)。NS の場の全温の超過は今回調べていない。

## 問い

1. この全温の超過は、Euler の比較 (V5) と設計チェーン (初期線のピン止め・出口較正) にどの程度効く可能性があるか。MOC の判断を進める前に、この異常の切り分けを先にすべきか。
2. 切り分けの最初の一手 (安い順に)。例: 既存の場で超過の分布を壁からの層・x で詳しく見る (0 step)、全温の輸送の保存性を局所に見る (0 step)、NS の場 (run_0147 など) にも同じものがあるか見る (0 step、AWS の保存物)、設定を変えた短い run (1 次・リミッタ・SLAU の壁の扱い) の A/B。どれを、何を判別するために。
3. MOC の判断は、この切り分けを待つべきか、並行して NS の評価 (上流の多項式化の NS の後の V5′) へ進めてよいか。
4. この異常は別の plan に切り出すべきか (対象の範囲と担当)。

## 読んでよいもの

- 上記 plan、`plans/accepted/tooling-nozzle-throat-monotone-r2.md`
- `case/45.isobutane_m6_d155/_band_ab/{moc_v5c_thermo_ab.json, moc_v5c_thermo_ab.log, moc_v5b_ext_eval.json, moc_v5_euler_eval.json}`
- `case/45.isobutane_m6_d155/{moc_v5c_thermo_ab.py, moc_v5_euler.py}`、`case/45.isobutane_m6_d155/README.md` (run 一覧)
- `solver_density_cuda/tools/total_quantities.py`、`procedures/solver-settings.md`、`methods/` の境界・node の節 (すべり壁)

編集は禁止。
