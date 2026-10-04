# 諮問ブリーフ: D-7275 試験 26 の残差プラトーが出口の角に局在し、2 回の対処で解けない (AGENTS 条件 2)

- plan: `plans/active/case-hypersonic-gap-heating-validation.md` §4.12、§6 G15、§5.1 #61 の末尾
- 台帳・事前登録: `case/60.flatplate_d7275_m7/README.md`、`case/60.flatplate_d7275_m7/acceptance.json` (T4-0b-WD・T4-0b-OUT の outcome)
- 生成器: `case/60.flatplate_d7275_m7/gen_runs.py`、メッシュ `case/60.flatplate_d7275_m7/mesh/fp_d7275_y3.geo` (と `_A`/`_B`)
- 前回の諮問: `notes/reviews/2026-09-27-t4-0-results-and-gap-handoff-diagnose.md`
- run は AWS。数値は下に全部書く。読まないこと: plan の §4.12・§6・#61 以外。

## 1. 観測事実

- 設定: 2D node、SLAU 2 次・limiter 2、SST (dilatation 2・Kato-Launder)、燃焼ガス 5 成分 NASA-9、全域 FP64、陰解法 block-DPLUR cfl 1.5 relax 1.0、
  領域 x ∈ [−0.1, 2.8] m・y ∈ [0, 0.5] m、平板 (等温 300 K) は x ∈ [0, 2.6] m、上流 x < 0 と下流 x > 2.6 の底面は slip、上端 y = 0.5 は slip、出口 x = 2.8。
- 残差 (run_0001 の判定区間 + run_0002): 2.2–3.9 dec で下げ止まり (rms_ro ≈ 1e-6)。FP64 なので丸め床ではない。
- **残差場** (`FORGE_OUT_RESIDUALS=1` + extraFields、run_0007): res_ro・roUx・roUy・roe の二乗和の **99.99–100 % が x ≥ 2.6 m (節点の 6.3 %)**。
  最大は**出口∩上端 slip の角** (x ≈ 2.79 m、y ≈ 0.47–0.50 m) と**出口∩下流 slip の角** (y = 0)。比較域 1.07–2.46 m と境界層 (wd < 3 cm) は 0.00 %。
- **対処 1 (壁距離)**: 下流 slip を別 physID にして壁距離に含める A/B → 末尾残差 0.53–0.73 倍 (1/10 に届かず)、q の差 ≤ 2e-6。
- **対処 2 (出口の種類)**: outlet_statPress (Ps = p∞) → outflow → 同じ長さの対照と末尾残差 1.00 倍、残差は依然 100 % が x ≥ 2.6 m、q の差 ≤ 5e-6。
- 比較量 St (比較 10 点) は全 run で ALL STEADY、R 平均 1.313 は 3 run で同じ。

## 2. 期待値と出典

- G15 の前提ゲート: check_convergence の同一設定区間 PASS。
- 現行レシピ (procedures/recommended-settings.md) は超音速 node 出口を outflow とする。

## 5. 仮説 (未確認)

- H1: 上端 slip∩出口の角ノード (2 つの境界の共有ノード) の扱いが、定常解を持たない (境界条件が矛盾する) — 前縁の Mach 波が上端か出口を通る位置と関係するかもしれない。
- H2: 下流 slip∩出口の角 (後流の中心線) の問題。
- 比較量には効いていない (2 回の A/B で 5e-6 以下)。

## 6. 聞きたいこと

1. **G15 の収束ゲートをどう扱うか**: (a) 角の原因をさらに追う (次の一手を 1 つ)、(b) 比較量に効かないことを 2 回の A/B と残差場で示したので、「残差は出口の角に局在・比較域の残差は 0.00 %」と明記して判定区間を比較域に限る (事前登録の変更になる)、(c) 領域を平板端 (x = 2.6 m) で切って角を比較域から遠ざける/無くす別メッシュで取り直す、のどれか。
2. (c) を選ぶなら、何を事前に決めておくべきか。
