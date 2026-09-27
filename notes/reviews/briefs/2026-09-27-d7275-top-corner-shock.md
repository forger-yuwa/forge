# 諮問ブリーフ: D-7275 試験 26 — 平板を出口まで延ばしたら残差は上端 slip∩出口の角 (前縁衝撃波の到達点) へ移った

- plan: `plans/active/case-hypersonic-gap-heating-validation.md` §4.12、§6 G15、§5.1 #61 の末尾
- 台帳・事前登録: `case/60.flatplate_d7275_m7/README.md`、`case/60.flatplate_d7275_m7/acceptance.json` (T4-0b-WD / OUT / RELAX / EXT / WALL の outcome)
- 前回の諮問: `notes/reviews/2026-09-27-d7275-outlet-corner-plateau-diagnose.md`、同 brief `notes/reviews/briefs/2026-09-27-d7275-outlet-corner-plateau.md`
- メッシュ生成器: `case/56.gap_tp1187/gen_mesh.py` (`--wall-dn` を追加)。読まないこと: plan の §4.12・§6・#61 以外。

## 1. 観測事実

- 設定は前回ブリーフと同じ (2D node、SLAU 2 次、SST、燃焼ガス 5 成分、全域 FP64、block-DPLUR cfl 1.5 relax 1.0、領域 x ∈ [−0.1, 2.8]、y ∈ [0, 0.5]、上端 slip、出口 outlet_statPress Ps = p∞)。
- **T4-0b-WALL** (事前登録 commit 679eb81e、ユーザ承認): メッシュ `fp_d7275_y3_wall` は元と形状・分割が同一 (.geo の Physical 行だけ差)、下流 2.6–2.8 m の底辺を slip → 等温壁 300 K。run `case/60.flatplate_d7275_m7/run_0012_t26_wall` (run_0002 最終場から interp_field、12 量、40k step)。
  - check_convergence: `NOT CONVERGED (stalled/plateau)`。流れ・化学種 2.0–2.7 dec で横ばい、rms_roK 6.5 dec・rms_roOmega 7.4 dec で下降中。
  - step 16–20k の平均残差 / 対照 run_0009 (元メッシュ、同設定、run_0002 から 20k) の同区間: ro 0.27、roUx 0.26、roUy 0.94、roe 0.29、roY 0.27、roK 0.006、roOmega 2.8。step 36–40k もほぼ同じ (roOmega 1.4)。事前登録の読み = **not_resolved** (1/10 に届かず)。
  - **残差場 (step 40000)**: res_ro・roUx・roUy・roe の二乗和の 99.7–99.9 % が **上端 slip 沿い y 0.45–0.50 m、x 2.70–2.81 m** (最大 (2.758, 0.500))。x ∈ [2.75, 2.78) が約 53 %、[2.78, 2.81) が約 22 %。壁近傍・比較域は ≤ 0.1 %。延長前は二乗和の大半が slip 後流にあった (T4-0b-EXT で 96–98 % が 2.6–3.5 m の slip 区間) — それは消えた。
  - res_roOmega の 75 % は延長した壁の第一層 (x ≥ 2.6、y < 1 mm)、最大 (2.692, 0)。ただし roOmega 残差は下降中。
  - **前縁衝撃波の位置** (P > 1.01 p∞ の最上点): x = 0.5/1.0/1.5/2.0/2.5/2.7/2.8 m で y = 0.088/0.166/0.250/0.333/0.420/0.445/0.472 m (約 9.5°、P 最大 ≈ 1.08 p∞)。上端 (y = 0.5) の P は x = 2.7 で 1.001、2.8 で 1.005 p∞ → **衝撃波の裾が上端∩出口の角に掛かっている**。
- 比較量: St 比較 10 点は run_0002 と表示桁で一致 (R_A 平均 1.313)、パネル平均 q_w 89.02 kW/m²。check_quasisteady (St 10 点) ALL STEADY。
- 補足: 事前登録注記の「比較 10 点は x ≤ 2.46 m」は誤りで、最下流点は x = 2.546 m (平板端の 5.4 cm 上流)。

## 2. 期待値と出典

- G15 前提ゲート: check_convergence の同一設定区間 PASS。
- Cary (case/59) は平板端 = 出口で全 9 系列 PASS。Cary の上端・高さと衝撃波の関係は未確認 (比較していない)。

## 3. 実施済みの対処 (AGENTS 条件 2 は既に成立)

壁距離 (run_0004/0005)・出口の種類 (run_0008)・緩和 (run_0009/0010)・slip 延長 (run_0011)・壁を出口まで (run_0012)。比較量の変化はすべて 5e-6 以下。

## 5. 仮説 (未確認)

- H1': 前縁衝撃波が上端 slip∩出口の角に当たり、角ノード (2 境界の共有) で反射/外挿の扱いが定常解を持たない。対処案: 高さ H を 0.5 → 0.8 m 程度に上げ、衝撃波を出口の上端から十分下で抜けさせる (壁は出口まで、形状・分割は x 方向同一、y は同じ y1・等比で層を足す)。
- H1'': 上端を slip でなく超音速の流入/外挿に変える案もあるが、BC の変更は影響が広いので H を上げる方を先に考えている。
- roOmega の壁第一層は延長壁の入口遷移 (初期場が slip 後流の場) の過渡で、下降中なので待てば下がると見ている (未確認)。

## 6. 聞きたいこと

1. 次の一手を 1 つ: (a) H を上げる (いくつ・何を事前に固定するか)、(b) 上端 BC を変える、(c) その他。
2. その run の事前登録の読み (PASS 以外の合格条件を置くなら何か) と、run_0012 を対照にするか。
3. 比較量が 6 回の変更で 5e-6 以下しか動かないことを、G15 の判定にどこまで使ってよいか (前回は「比較域だけで判定しない」だった)。
