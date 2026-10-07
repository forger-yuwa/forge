# 諮問: V5b (Euler の延長の診断) も保留 — 解釈と、MOC の軸処理を生産に入れるかの判断材料

日付 2026-10-07。諮問先 codex (diagnose)。エスカレーション条件 7 (result 段の解釈を確定する前) と 4 (plan に無い次の手へ進む前)。
plan: `plans/active/discretization-moc-axis-limit-and-corrector.md` (§6 V4・V5・V5b、§9 の 2026-10-07 の記録を全文読むこと)。作業ツリー `/home/sano/work/forge-integ-1005` (HEAD dda97dbb)。

## 観測事実

- V5b の評価器 (sha256 0ea2fe9a…、登録 4dd94af2、実装時の定義は §9) の出力: `case/45.isobutane_m6_d155/_band_ab/moc_v5b_ext_eval.json`・`.log`。親 (0〜18000) の時系列は `_band_ab/v5_series/<run>/wallfit_series_v5.csv`、子 (延長、local step + 18000) は `_band_ab/v5b_series/<run>/wallfit_series_v5b.csv`。
- **総合: 保留**。7 本すべてが減衰側ではない。run_0154 (腕 B r1) と run_0157 (腕 M r1) は exit_M_dev が持続振動側、他の量は判別不能。残り 5 本は全量で判別不能 (exit_core_M だけは 6 本で減衰側)。
- 窓 A (通算 24000〜36000) → 窓 B (42000〜54000) の平均の移動は、許容幅の 1/10 を超える量が多い (腕 B・腕 M とも P 傾き +1.5〜6.0e-3、M 波 +0.7〜3.1e-4)。延長しても一方向に動き続けている量がある。
- **窓 B の参考値** (判定ではない):
  - D = 腕 M − 腕 B: M 波 +1.8e-4 (2SE 1.5e-4)、P 波 +1.4e-3 (8.3e-4)、オーバーシュート −2.3e-3 (1.3e-3)、出口規格化オーバーシュート −3.0e-3 (1.2e-3)、P 傾き −4.4e-2 (2.4e-3)、exit_M_dev +1.5e-5 (5.0e-6)。出口コア M 腕 B 5.999987・腕 M 6.000029。
  - ISEN − 腕 M (同じ壁、IC だけ違う): P 傾き +4.9e-2 (2SE 8.0e-3)、出口規格化オーバーシュート +3.4e-3 (1.8e-3)、出口コア M −7.2e-5 (1.5e-5)。
  - 窓 B の P 傾きの平均: 腕 B 0.207・腕 M 0.164・ISEN 0.213 (%pt)。**同じ壁で IC だけを変えると 0.05 違い、その差は腕 B と腕 M の差 (−0.044) と同じ大きさ**。
- 腕 B・腕 M の IC は run_0114 (旧壁の Euler の収束場) の保存量を番号写像で写したもの。ISEN は等エントロピー IC から段階起動。
- V5 (窓 6000〜18000) の参考値でも、P 傾きの D は −4.3e-2 で、ISEN − 腕 M は +2.7e-2 だった。
- 単調壁の採用 (2026-10-06、`plans/accepted/tooling-nozzle-throat-monotone-r2.md` §6 E′) は、腕 A (今の壁) 対 腕 B (単調壁) を同じ種類の Euler の比較で行い、P 傾きの差を検出した (0.1926 → 0.2008 %pt、Δq の 37 %) と記録している。両腕とも IC は run_0114 の写像。
- MOC の変更そのものの証拠 (CFD を使わない): 放射源流の厳密解で壁の誤差が 2 次に (V0、次数 1.99)、case/45 で始点の角度差 0.0219° → 0.0054°、r″ の山 1.0e-2 → 3.5e-5、設計壁の変化最大 6.2 µm (x = 0.40 r_t)、出口 7e-8 µm (V4)。
- Euler の最終場では、どの腕でもスロート付近の壁際の数千節点で全温が Tt を数百 K 超える (§9 の観測)。原因は調べていない。
- run_0154〜0160 の post-check (NaN・残差判定) は実行中。親の run_0150〜0153 は NaN なし、本段の残差は腕 B・腕 M が stalled/plateau、ISEN が still converging だった。

## 期待値と出典

- V5 の目的: 新しい壁が生産の Euler 格子で今の単調壁より悪くないこと (E′ の許容幅)、出口較正の要否、IC 依存の有無。
- 許容幅 Δq: M 波 0.001、P 波 0.010、オーバーシュート 0.003、出口規格化オーバーシュート 0.003、\|P 傾き\| 0.03 (%pt)、\|出口コア M − 6\| 1.8e-4。

## 問い

1. P 傾きの「同じ壁で IC により 0.05 違う」は何を意味するか (Euler の解が IC に依存する別の状態に落ちている、ゆっくりした過渡がまだ続いている、抽出の問題、など)。確かめ方 (CFD 0 step でできる切り分けを優先)。これは単調壁の採用 (E′) の根拠にも影響するか。
2. この Euler の設定で、壁の µm 級の違いを E′ の許容幅で比べることに意味があるか。無いなら、MOC の変更を生産に入れるかの判断を何で行うべきか (設計側の証拠だけで進める、NS で見る、Euler の比較方法を作り直す、など)。
3. 次の手として推奨するもの。候補: (a) MOC の変更は設計側の証拠 (V0・V4) で生産候補とし、Euler の V5 は「この設定では判別できない」と記録して閉じ、NS (V5′、上流の多項式化の NS の後) へ進む。(b) Euler の IC 依存と非定常を先に切り分ける。(c) 生産に入れない。どれを、どの条件で勧めるか。ユーザに判断を仰ぐ材料として、何を並べるべきか。
4. スロート付近の全温の超過は、上の判断に関係するか。

## 読んでよいもの

- 上記 plan、`plans/accepted/tooling-nozzle-throat-monotone-r2.md` (§6 E′・§9)
- `case/45.isobutane_m6_d155/_band_ab/{moc_v5_euler_eval.json, moc_v5b_ext_eval.json, moc_v5b_ext_eval.log, wallfit_series_v5.json, wallfit_series_v5b.json, moc_v5_ic_inspection.json}`、`_band_ab/v5_series/`・`_band_ab/v5b_series/`
- `case/45.isobutane_m6_d155/{moc_v5b_ext_eval.py, moc_v5_euler_eval.py, eval_wallfit_euler.py, moc_v5_euler.py}`、`solver_density_cuda/tools/check_quasisteady.py`
- `methods/design/overview.md` (軸上の解析極限の節)

編集は禁止。
