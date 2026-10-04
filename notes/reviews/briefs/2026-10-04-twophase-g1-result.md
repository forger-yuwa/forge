# 諮問ブリーフ: G1 診断の結果の扱い (事前登録の float/double 判定が非正規化数の範囲で FAIL)

エスカレーション条件 3 (事前登録の比較が FAIL・判定不能) と 7。plan: `plans/active/condensation-two-phase-default.md` §5.1 #4 (事前登録の詳細と実装上の選択)、#4r (結果、数値はすべてここ)、§6 G1。
成果物: `notes/investigations/2026-10-04-twophase-g1/evidence/` (G1_VERDICT.txt・forge_run.log・IC_FROM.txt)、判定 `notes/investigations/2026-10-04-twophase-g1/g1_judge.py`、診断実装 `solver_density_cuda/cuda_forge/twoPhaseFaceDiag_d.cuh`・`speciesTransport_d.cu` 診断節・`main.cpp runTpFacesDiag` (commit 3d6e853f)。

## 観測事実
- (ii) 液: 574 面で abs(Jl_f − Jl_d) > 8ε₃₂·A_l (A_l = abs(ct·geo)·(abs(ρg₀/ρ₀)+abs(ρg₁/ρ₁)))。574 面はすべて abs(Jl_d) < FLT_MIN (1.18e-38)、絶対誤差最大 5.4e-45。例: rg0 4.4e-39, rg1 0 → Jl_f −7.006e-45, Jl_d −6.311e-45。正規化数の範囲の面 (abs(Jl_d) ≥ FLT_MIN) では最大比 0.19・超過 0。分子蒸気は全面で最大比 0.068。
- (i) 生の比較: 乱流域 (μt,f/μ_f ≥ 1) 50446 面で max abs(J_w^A) 1.42e-10 vs 0.05·max abs(J_l^B) 7.18e-7。
- (iii) 帯 0.1〜0.4・0.4〜1.6 mm は上下とも支持、帯 0〜0.1 mm は棄却 (分子蒸気 Σ −8.9e-9 = 壁向き、液 Σ ~−2.4e-23 = 液がほぼ無い)。事前登録は「面のある全ての (帯 × 壁) が支持のときだけ支持」。
- P0 (OFF 写し vs 本番) ok (比 0.28〜0.30)。

## 問い
1. (ii) の FAIL をどう扱うか。事前登録の尺度が非正規化数 (相対精度が無い範囲) を想定していなかった欠陥として、絶対の床 (例: 8ε₃₂·A_l + c·FLT_TRUE_MIN、または真値 < FLT_MIN の面を「float で表現不能 = 判定対象外」と記録) を足して判定し直してよいか。結果を見た後の改訂になる — 許されるならその条件と記録の仕方を。許されないなら G1 をどう再実施するか (新しい入力で改訂した尺度を事前登録して取り直す、など)。
2. 物理的に、液の非正規化数 (ρg ~1e-39) の面がある (液がほぼ無い領域の境界) ことは、ON の作用素や本番で問題になるか (float の FTZ、非正規化数の流束が残差に入ること)。
3. (iii) の壁直近帯の棄却をどう読むか。液が無い帯では「液は壁向き」は検出不能で、蒸気が壁向きなのは何を意味するか (壁で総水分が減る観測 #4j と整合か)。事前登録の「全帯で支持」は判別として妥当だったか。(iii) は既定化のゲートではなく記録。
4. G1 として「既定化を止める/進める」の判断に必要なものは揃ったか。G2 (生産レシピ) に進んでよいか。
読んでよいファイル: 上記、`plans/active/condensation-two-phase-transport.md` (#4j・#4k)、`methods/condensation.md` §7c、`solver_density_cuda/cuda_forge/twoPhaseDiffusion_d.cuh`。
