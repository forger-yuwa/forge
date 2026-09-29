# 諮問ブリーフ: 共役ベンチマークの再評価後の結論と閉じ方 (2026-09-30)

AGENTS.md 条件 **7** (result 段の解釈を確定する前)・**3**。plan `plans/active/boundary-cht-conjugate-benchmarks.md` §5.1 #6・#6b。
前回の諮問 `notes/reviews/2026-09-30-conjugate-benchmarks-results-diagnose.md` は全件採用・再評価済み (#6b)。

## 読んでよいファイル
plan §5.1 #6・#6b の行、case/64 と case/65 の run_0005〜0011 の `EVAL_CONJ.txt`・`SERIES_CONJ.txt`、`case/64.conjugate_pipe_wall/ab_levels3/EVAL_CONJ_a1_r16_levels3.txt`

## 再評価の結果 (要約。詳細は #6b)
- A: r32 (A1・A2) と r64 (A1) は全項目 PASS。A2 r64 は Q_up だけ判定不能 (領域切断の U)。r16 は A1・A2 とも FAIL (一部判定不能)。A/B で r16 の FAIL は forge の粗格子誤差に帰属。準定常は全節点 PASS。
- C: 主判定 6 本 PASS、準定常 PASS、G-cons PASS。前提ゲートの流体収束と G-if は FAIL (界面残差は板の端に集中、窓内 0.14 W/m²)。前回の推奨どおり C は判定不能として保留。

## 私の閉じ方の案
- A: 「固体の軸方向伝導・熱抵抗が効く軸対称の共役伝熱 (上流へ 15 % / 35 % の熱が回り込む) で、forge の界面温度・熱流束は、forge の流れ場を固定した独立参照解と N_r 32・64 で登録許容内 (A2 r64 の Q_up のみ判定不能)。最も粗い N_r 16 は不合格 (forge の空間離散化誤差、格子倍増で約 1/4)」— 代表格子の選定は事後と明記。
- C: 判定不能として保留し、前縁・後縁の slip 境界と共役壁の角での揺れ (既知の node slip 欠陥が疑わしい) を別の残作業に切り出す。
- plan は A を成果として result レビューに回し、C は未完了のまま active に残すか、plan を A と C に分けるか。

## 諮りたいこと (推奨を 1 つに絞る)
1. A の結論文として上の案は妥当か (過大な主張・欠けた限定)。
2. C の扱い (保留のまま plan を閉じる / plan を分ける / C を直してから閉じる) の推奨。
3. result レビューに進んでよいか、進む前に足りないもの。
