# 諮問: farfield #4f (新輸送での格子差・幅) と SERN R9/R9b (種 DB 段 3 + #14) の結果解釈 (2026-10-02)

関連 plan: `plans/active/boundary-node-farfield-characteristic.md` (§5.1 #4e・#4f)、`plans/active/tooling-nozzle-sern-chain.md` (R9 行)。
run 一覧: `case/46.sern_design/README.md` (run_1014–1027)。run はすべて AWS `~/forge-r8/case/46.sern_design/` (手元には無い)。
評価スクリプト: `case/46.sern_design/v3_farfield_eval.py` (窓条件: 末尾 10000 step 平均 / a = max|値−平均| / 前窓差 ≤ 0.1ε / D = |Δ平均| + a1 + a2)。
過去の諮問: `notes/reviews/2026-10-01-farfield-r8-status-review-diagnose.md`。

## 1. 観測事実

### 1.1 R9 (段 3、CEA 熱物性) — コード bfd8c3c8 (+ランナー修正 dba4b582)、バイナリ sha256 2ac853cf
- 2D m6_on、`problem_moo_frozen_tp_cycle3op_transport.yaml`、ランナー既定の段階起動、3 反復ずつ。A = run_1008–1010 (旧 dc560910)、B = run_1014–1016。
- 全 6 本 GATES PASS・4 量 STEADY・残差は全列プラトー (0.2–0.4 桁、上昇なし、A と同じ)。
- |B−A| / N (N = max(3×群内最大差, 1e-6)): C_T 5.5e-7/5.0e-5、C_T_with_shear 1.4e-6/4.8e-5、C_L 2.0e-6/2.0e-4、C_M 4.3e-5/4.5e-3 → 事前登録で PASS。
- B の群内最大差は A の 4–5 倍 (C_M 1.5e-3 vs 3.6e-4)。原因は未調査。
- 初回投入はランナーの段間継承 (`_species_signature`) が段 3 の記録形式 (`Tbounds`/`nasa9_intervals`) を読めず KeyError で停止 → 修正・試験追加して同名で再投入。

### 1.2 #4f (段 3 バイナリ、生産 YAML = side_far farfield + 種ごとの輸送、リミッタ基準値固定)
- run_1017 = g3 2.50 H を run_1012 (旧バイナリ、同設定) から +20000: run_1012 との D = C_T 8.4e-6・C_T_with_shear 8.3e-6・C_L 2.9e-5・C_M 6.4e-4 (τ = 0.2ε 内)。
- run_1018 (g4 2.50)・1019 (g3 3.42)・1020 (g3 4.35) は窓条件未達 (C_T_with_shear・C_M) → +20000 の run_1021・1022・1023 で全量 OK、全 run GATES PASS。延長の restart は種 DB #3d どおり許可なしで通った。
- (a) 幅 (g3): 2.50 vs 3.42 の D = C_T 7.5e-6・C_T_with_shear 8.4e-6・C_L 8.7e-5・C_M 2.31e-3。2.50 vs 4.35 = 6.8e-6・8.2e-6・8.0e-5・2.17e-3。3.42 vs 4.35 = 7.8e-6・7.6e-6・3.4e-5・7.3e-4。ε = 5e-4 (C_T 系・C_L)、5e-3 (C_M) → 全量 D ≤ ε。
- (b) G = g4 − g3 (2.50 H、run_1021 − run_1017): C_T −1.10e-4、C_T_with_shear −1.94e-4、C_L −5.50e-4、C_M +1.67e-2。|G| + D(2.50 vs 4.35): 1.16e-4・2.03e-4・6.30e-4・1.88e-2。§8 許容 C_T・C_L 2e-3、C_M 5e-2。
- 旧輸送 (#4c) の G: C_L −0.00052、C_M +0.0159。

### 1.3 R9b (#14、LJ 出典の既定 [gri30, svehla1962]) — コード f28ca2fa、バイナリ sha256 bec8db5f
- 段 3 以降のソルバ差分は種データ・DB/輸送 DB 読み込みだけ (カーネル変更なし)。
- 解決済み LJ: EXH σ 3.36827→3.36691、ε/k 212.555→212.484 (構成種の質量分率平均)、AMB 不変。`speciesDiffusionMethod` は既定 1 (kinetic 混合平均)。
- 2D: B = run_1024–1026 (3 本 GATES PASS・STEADY)。B − A / N: C_T −6.1e-7/5.0e-5、C_T_with_shear −6.4e-8/4.8e-5、C_T_friction +5.5e-7/2.2e-6、C_L −5.0e-6/2.0e-4、C_M +1.05e-4/4.5e-3 (事前登録で合否なし、記録)。
- 3D: run_1027 = run_1017 から新バイナリで +20000 (ハッシュ変更のため許可なしの restart は拒否 → `--force-species`)。run_1017 との D: 6.4e-6・6.4e-6・2.6e-5・5.8e-4、全量 ≤ τ → 事前規則で #4f の結論を持ち越し。

## 2. 期待値と出典
- 判定規則・ε・τ・§8 許容はすべて run 前に plan に書いた (farfield plan #4e/#4f、SERN plan R9 の事前登録)。
- 依頼元の予測: 段 3 は H2O 分子量の分 (1e-6 程度)、#14 は LJ 由来の拡散係数変化 (最大 +34 %、OH–H; EXH は lump なので薄まる)。

## 3. 再現条件
- 2D: `python3 -m forge_design.evaluate.runner_sern case/46.sern_design/problem_moo_frozen_tp_cycle3op_transport.yaml <run> --op m6_on`、`FORGE_CUDA_BLOCKSIZE=128`。
- 3D: `case/46.sern_design/prod_restart_setup.py` (初回) と `batch_4f*.sh` / `batch_r9b_3d.sh` (scratchpad、AWS に配置) の手順。

## 4. 仮説・解釈 (確定したい内容)
1. 生産 3D (g3、2.50 H farfield、種ごとの輸送、#14 後) の必要幅は 2.50 H、格子差込みで §8 内 — 限定: m6_on・g3/g4 の 2 水準・固定基準値。
2. R9・R9b は生産の力係数をノイズ床以下でしか動かさない (2D 3 反復、3D 1 本)。
3. 「新輸送で G がほぼ不変」は、輸送切替 (C_L +0.23 %・C_M −0.22 %) が g3・g4 で同程度に効いたため、と読んでよいか (g4 で輸送 A/B はしていない)。

## 5. 問い
1. 4 章の解釈に誤り・言い過ぎはあるか。特に (a) 2 水準 (g3/g4) で「格子差込みで §8 内」と書く範囲、(b) 3D の持ち越し判定を 1 本で行ったこと、(c) 1.1 の B 群のばらつき増大を放置してよいか。
2. farfield plan を accepted にするための残り (V2c 判定不能・V2d-2 時間精度・独立参照未収束、result 段 codex レビュー) について、前回諮問 (2026-10-01) の道筋から変えるべき点はあるか。
3. 生産設定の「限定運用」(精度認定しない) を解除してよい条件は満たされたか。
