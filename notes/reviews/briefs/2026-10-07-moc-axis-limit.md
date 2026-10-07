# 諮問: 逆 MOC の軸上ソース項を解析極限 θ_r にし予測修正を収束させる方針 (§4) と検証計画 (§6)

日付 2026-10-07。諮問先 codex (diagnose)。エスカレーション条件 1 (plan §4・§6 を新規に書く)。
plan: `plans/active/discretization-moc-axis-limit-and-corrector.md` (全文を読むこと)。作業ツリー `/home/sano/work/forge-integ-1005`。

## 観測事実

- **現行の軸上の扱い** (`design/forge_design/geometry/moc_kernel.py:406-412`): 軸上の端点 (r ≤ 1e-9、または r < 0.05 × 相手の r) では、sinθ/r を相手の点の値で代用している。相手も軸上なら 0。
- **予測修正** (`moc_kernel.py:428`、`n_corr=2`): 予測 1 回と修正 2 回。軸の 1 段目で修正子が振動する (比 ≈ −1/2)。
- **再出発試験** (`case/45.isobutane_m6_d155/throat_moc_restart_test.py` → `_band_ab/throat_moc_restart_test.json`):
  - 生産の手順では、網の列を初期線に戻して作り直すと、列の 1 段目で −0.05〜−0.11° ずれる。
  - 修正子を収束させ、軸上のソース項を常に 0 にした試験用の手順 (K2c) では、≤ 1e-7°。
- **事前試算** (`throat_moc_fix_probe.py`・`throat_moc_fix_probe_n4800.py` → `_band_ab/throat_moc_fix_probe.json`、case/45 の単調壁の生産問題): 第 1 点 (n_start 41、x = 0.025) の角度差 Δθ = θ − atan(x/R)。
  - 生産 (n_axis 2400・dx0 0.03): 0.022°。K2c: 0.018°。
  - 修正子だけを収束させた版 (conv) と K2c は、n_start 321 でもほぼ同じ (0.019° 対 0.019°)。
  - n_axis 4800 で dx0 を 0.03 / 0.015 / 0.0075 にすると、生産 0.020 / 0.012 / 0.009°、K2c 0.017 / 0.010 / 0.008°。
  - 壁の変化は最大 5〜6 µm (x ≈ 4。n_axis 2400 → 4800 による)。出口 ≤ 0.005 µm。
- **所要時間** (design_chain 1 回): 生産 2.3 s、n_axis 4800 で 6 s、K2c (修正子 20 回固定) で 8.6 s (n_axis 2400)、28 s (n_axis 4800・dx0 0.0075)。
- **θ_r の大きさ**: x_A (M 1.215、M′ 0.635、γ ≈ 1.27) で θ_r ≈ 0.103 rad/r_t。

## 期待値と出典

- A11 (`plans/accepted/discretization-moc-axisymmetric-source-term.md`): 源項を点自身で評価するよう直し、放射源流の厳密解で誤差が 10〜40 倍下がり、収束次数が 1.02 → 1.7 になった。
- methods の「今後の課題」に、軸上の解析極限 ∂θ/∂r|r=0 = −½ d ln F/dx が挙がっている。

## 方針 (plan §4)

- 軸上の端点の sinθ/r に θ_r = ½ √(M² − 1) · ν_M(M) · M′(x) を使う。ν_M は MOC 本体と同じ気体モデルの ν(M) から取る (semi-perfect は表の 3 次スプライン微分)。M・M′ は軸則から、x_A ではアンカーから取る。
- 修正子は、θ・ν の変化が 1e-12 を下回るまで回す (上限 50 回)。
- 新しい YAML キーで選び、既定はビット同一。

## 問い

1. θ_r の式 (semi-perfect を含む) と、軸上の端点への入れ方は正しいか。
   - 予測段と修正段で軸上の点の θ_r をどう扱うか。
   - 軸上の点の M・M′ を軸則から取ることと、網の中で計算した値との整合。
   - 初期線の軸端の扱い。
2. 修正子の収束判定 (1e-12、上限 50) は妥当か。軸の 1 段目の振動が収束しない場合の扱いはどうするか。
3. `AXIS_LIMIT_FRAC` の分岐 (軸に近い軸外の点) を変えないことは妥当か。新しい極限と混在して不整合を作らないか。
4. §6 の V1〜V6 の合格条件は妥当か。V2 (放射源流) で軸上の極限の効果が見えるか (放射源流では軸上の M′ が既知)。V4 の「case/45 に入れる条件」は結果を見る前の登録として十分か。
5. 生産に入れる (case/45) 場合、Euler 1 本 (V5) で足りるか。NS をやり直す条件は何か (壁の変化が µm 級の前提)。

## 読んでよいもの

- 上記 plan
- `plans/accepted/discretization-moc-axisymmetric-source-term.md`
- `plans/accepted/tooling-nozzle-throat-monotone-r2.md` (§3 仮説 H、§9 2026-10-06〜07)
- `design/forge_design/geometry/{moc_kernel.py, moc_inverse.py}`、`design/forge_design/evaluate/runner_axismach.py` (`design_chain`)
- `case/45.isobutane_m6_d155/{throat_moc_restart_test.py, throat_moc_fix_probe.py, throat_moc_fix_probe_n4800.py}`
- `methods/design/overview.md` (A11 の節と「今後の課題」)

編集は禁止。
