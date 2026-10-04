# 諮問ブリーフ: G1 の解釈の確定と G2 への移行

エスカレーション条件 7 (result の解釈を確定する前)。plan: `plans/active/condensation-two-phase-default.md` §5.1 #4・#4r・#4s・#4sr・#4rr、§6 G1・G2。
前回諮問: `notes/reviews/2026-10-04-twophase-g1-result-diagnose.md` (指数スケール A/B と改訂尺度の条件)。
成果物: `notes/investigations/2026-10-04-twophase-g1/evidence/` (G1_VERDICT.txt 旧基準、G1_VERDICT_revised.txt 改訂、SCALE_AB_VERDICT.txt)、
判定 `notes/investigations/2026-10-04-twophase-g1/g1_judge.py` (`--revised-4sr`)、A/B `scale_ab.py`・`solver_density_cuda/tests/unit/tp_jl_scale_ab.cu`。

## 観測 (要約、数値は plan #4sr・#4rr)
- A/B: 保存入力 A は旧尺度 574 面超過、液 × 2⁶⁴ の B は全面旧尺度内 (最大比 0.19)、A の Jl_f は保存値とビット一致。
- 改訂尺度 E_abs = (3·abs(ctg)+1)·2⁻¹⁵⁰ (導出は #4sr) で (ii') PASS、液の最大比 0.995。(i) PASS (比 0.000)。(iii) 規則上棄却 (記録)。

## 問い
1. G1 を「改訂尺度のもとで成立」と確定してよいか。最大比 0.995 は導出の上界がほぼ実現したことを示す — 導出に見落とし (例: ctg = ct·geo 自体の丸め、ir の丸め、FMA の有無) があれば上界を超えていたはずで、余裕が無いこと自体をどう評価するか。導出をもう一段厳密にすべきか。
2. G1 の (i) の結果 (OFF の総水分拡散流束が実質 0) と (iii) を、既定化の根拠としてどう書くか (plan §3 の「CFD レベルの外部参照なし」の限定のもとで)。
3. G2 (生産レシピ: cfl_pseudo 4・implicitRelax 0.7・nStepInner 4・Solver 1・nn0) の事前登録は §6 G2 のとおりで足りるか。ノイズ σ を測る同条件反復の本数・長さ (48000 step を 2 本?)、IC (run_0482 res_48000 から ON で開始し OFF の場から始めない?)、判定区間、7 量の抽出マスクの数値照合の方法を具体化してほしい。比較相手 ON@cfl2 は run_0521 (通算 32000 step、cfl 2・relax 1・nStepInner 5)。
4. G3 (作用素収支・更新写像収支) は G2 と並行してよいか、G2 の後か。
読んでよいファイル: 上記、`plans/active/condensation-two-phase-transport.md` (#4j)、`procedures/recommended-settings.md`、`case/16.nozzle_wys/README.md`、`solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh`。
