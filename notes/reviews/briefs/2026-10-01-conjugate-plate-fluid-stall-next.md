# 諮問: 共役平板の流体停滞は連成でなく流体側 — 次の切り分けをどう設計するか

AGENTS.md 条件 7 (結果の解釈)・4 (plan に無い次の run)・1 (§4 の追加)。plan `plans/active/boundary-cht-conjugate-flat-plate.md` §4.1・§5.1 #1b2 (commit は本 brief と同時)。判定出力 `case/65.conjugate_flat_plate/AB_JUDGE_run0011_0012.txt`、判定器 `ab_judge.py`。

## 観測事実 (AWS FP64 `86115cb1`)
- B1 本段: `run_0011_ab_c1_n64_coupled` (warmup 5000) / `run_0012_ab_c1_n64_fixed_tw` (warmup 900000)。`run_0007_c1_n64` step 600000 から各 60000 step。
- 事前登録の判定: 「動的な連成だけが停滞を維持する」を棄却 (B/A 残差比 0.996〜0.999、両者横ばい max/min 1.002、A は元 run を比 1.00 で再現、界面残差最大位置は元・A とも x/L 0.0082)。
- 流体の局所変動の振幅 B/A 0.95〜1.03。前縁帯の界面熱流束 (iface_q_eff) の時間変動 (標準偏差/評価窓平均 q 1379 W/m²) は A 0.1146・B 0.1148 → **壁温を固定しても前縁の界面熱流束は 11 % 揺れる**。これが G-if ① (res_abs 678 W/m²) の主成分と見える (未検証)。
- `check_convergence.py`: 両 run NOT CONVERGED (stalled, 0.1 桁)。元 run は 600000 step で rms_ro 1.8 桁・rms_roUy 2.0 桁・rms_roe 1.8 桁で横ばい、rms_roUx だけ 4.6 桁。
- 変動の位置 (B、2000 step ごとの全場 30 枚、時間標準偏差/尺度): 最大は後縁の直後 x/L 1.006〜1.03 の y=0 (slip_down) で P 1.3e-3 q∞ (≈0.9 Pa)・T 8e-5×10 K・Uy 4.4e-5 U∞。次が前縁の直後 x/L 0.008 (板の壁上 P 8.1e-4)。上流 slip (x<0) は P 3.3e-4。
- スペクトル (B、毎 step 60000 点): 上流 slip の P は周期 7.2 step にパワーの 92 % (単一の鋭い振動)。前縁帯・下流は広帯域 (最大ピークでも 0.4〜1.7 %)。
- 設定: node、SLAU、limiter 2 (Venkatakrishnan、limiterScaled 既定 1・venkatK 既定 0.05)、陰解法 blockDPLUR cfl_pseudo 2・implicitRelax 0.7・nStepInner 4・timeIntegration 11、定常 (unsteady 0)、M0.1・Re_L 1e4・定数物性、低マッハ前処理なし。前縁 (slip→no-slip) と後縁 (no-slip→slip) は特異点。
- 既知の未修正欠陥: node slip + 接線方向の密度勾配で市松状の接線速度 (`notes/investigations/node-slip-tangential-density-spurious-flow.md`)。今回の最大変動は後縁直後の slip 上 (後流の温度=密度勾配がある) で、形は未確認。
- 登録の完了条件 (§6): C1/C2 × n16/32/64 の 6 本が流体収束・G-if を含めて合格。

## 仮説
H1: リミッタのチャタリング (特異点近傍でリミッタ値が step ごとに切り替わる) による有界振動。
H2: 陰解法の内部反復 (implicitRelax 0.7、nStepInner 4) の擬似時間の振動モード。上流 slip の周期 7 step はこれか H1。
H3: node slip の既知欠陥 (後縁直後の slip で後流の密度勾配)。
H4: M0.1 で前処理なしの低マッハ硬さ。

## 私の案
- 次の A/B は同じ起点 (run_0012 相当、壁温固定で連成を除いた流体単独) から 1 キーだけ変える。
  - (a) リミッタを凍結 (forge に凍結キーがあれば) または limiter 0 (1 次)。1 次は解が変わる。
  - (b) implicitRelax 1.0 / cfl_pseudo 1。
- 観測は B1 と同じ `ab_extract.py`・`ab_judge.py` の仕組み。

## 質問
1. 結果の解釈 (連成は原因でない、前縁の界面熱流束の揺れは流体由来) は妥当か。書いてよい範囲は。
2. 次の判別の設計 (どの 1 キーを変えるか、判別条件)。H3 を他と分けられる A/B は。
3. C の完了条件 (6 本全合格) に向け、停滞が解けない場合の扱いをどう先に決めておくべきか。
