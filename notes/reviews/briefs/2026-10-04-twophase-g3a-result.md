# 諮問ブリーフ: G3-a の結果の解釈と、射影と戦う機構への対処の方向

エスカレーション条件 7 (result の解釈を確定する前)、3 (§6 の事前登録との比較が FAIL)、6 (対処が `cuda_forge/` の数値の振る舞いを変える)。
plan: `plans/active/condensation-two-phase-default.md` §5.1 #4g3 (設計)・#4g3a (ΔF の扱い)・#4g3r (結果)・#4g3n (対処の候補)・#4pj/#4pjr/#4pjg (経緯)。
前回諮問: `notes/reviews/2026-10-04-twophase-g3-design-diagnose.md`、`notes/reviews/2026-10-04-twophase-projection-fight-diagnose.md` (射影の判別 #4pj)。

## 観測事実
- コード版: build 8a9b673b (AWS)。回帰 (診断分岐を入れた既定経路、ON/OFF、200 step、旧 6663cbc3 ×3 vs 新 ×2): 両方 PASS (`notes/investigations/2026-10-04-twophase-g3/evidence/g3areg_{on,off}.txt`)。
- 診断 run (0 step、更新なし): ON `case/16.nozzle_wys/run_0604`–`0606` (入力 = G2 の ON 生産レシピ最終場 run_0561–0563 res_48000、7 量 ALL STEADY)、OFF `run_0607`–`0609` (run_0567–0569)。判定 `evidence/G3A_pair{0,1,2}.txt` (`g3a_judge.py`)。
- 判定器の不具合 (結果を見る前の不具合として修正): 出口マーカーを種別 `outlet*` で探していたが case/16 は `outflow` → 出口面 0 で ΔF = 0。修正後 120/120 面。
- ΔF (3 本平均の ON−OFF、run 間の標本標準偏差から誤差 2·√(s_on²/3+s_off²/3)、時間平均ではない): g 2.213e-4 ± 3.3e-8、Q2 2.12 ± 4.6e-4、Q1 2.84e7 ± 3.1e4、Q0 5.49e14 ± 2.5e12、総水分 w 4.3e-9 ± 1.9e-8 (判定不能)。
- G3-a (1) 組立の閉鎖: 全 6 run・全 Ω・全成分 PASS (最悪 abs(E)/B 7.3e-3)。
- G3-a (2) 物理の収支 abs(ΣR_final)+B < 0.1ΔF: OFF は全 Ω PASS (w 除く)。ON: g PASS、**Q2・Q1・Q0 が全域 FAIL** (0.365 対 0.212 / 4.9e7 対 2.8e6 / 1.1e14 対 5.5e13)、Q1 は x35–70 wd<0.4 と x[70,出口] でも FAIL。3 本とも同じ。
- 射影集合 (run_0592 で射影が働いた 4953 節点、`extract_projection_mask.py`、座標チェックサム一致) の段分解 (`boundary_push.py`、`evidence/BPUSH_*.txt`、全成分正の 4755〜4831 節点):
  - 境界への張り付き (1e-3 以内、ρ_l = 1000 近似): y = √x (Q2² = Q1·Q3) が 89〜95 %、y = x² 0.3 %、x = 1 0.1 %。中央 x ≈ 0.27〜0.32、y ≈ 0.52〜0.55。
  - 余裕 H_up = ln Q1 + ln g − 2 ln Q2 (≥ 定数 が実現可能) の段ごとの dH/dt = Σ a_c R_c/(V q_c): ソース段 (cond_src − tp_diff) が 97〜99 % の節点で負、中央 −2.2e5〜−3.2e5 /s。二相拡散段は中央 +9.0e3 /s (負は 13 %)、移流段は中央 +1.7 /s、ソース後の段は 0。合計が負の節点 96〜99.7 %。
  - ソース分岐: 蒸発 97 % (4734〜4810)、成長 21、なし 0 (全成分正の中で)。飽和度 S の中央 0.27 (S ≤ 1 が 98.8 %)。成長分岐の Sg < 0 クリップ作動 0。
- 段ごとの残差の符号 (pair0、射影集合 4953): Q2 は二相拡散で増 (+、98 %) とソースで減、合計は増 (96.6 %)。Q1 はソースで減 (99.6 %)、合計は減 (93 %)。
- 先行の観測 (#4pjr/#4pjg): 射影節点は開始状態がすでに境界上、毎更新 Q2 +9.8e-4・Q1 −1.5e-3 (相対) の増分を射影が 70〜97 % 戻す。g・Q0 はほぼ不変。液の上限クランプ作動 0、θ は全成分共通、Q の非負カットは Q1 16・Q2 1 節点。

## 期待値と出典
- 事前登録 (plan §6 G3、#4g3): (1) abs(E) ≤ B、(2) abs(ΣR_final)+B < 0.1ΔF (液・総水分・Q_n)、差が誤差以下なら判定不能。
- 蒸発ソースの形 (`solver_density_cuda/cuda_forge/condensationSource_d.cuh` `cond_evap_source_rate`、condLimiterMode 1): 一様 ṙ で S_Q1 = q0ṙ、S_Q2 = 2q1ṙ、S_g = 4πρ_l q2ṙ、**S_Q0 = 0** (「数密度は消滅まで保存。消滅 r30 < 2r_min と Q0 = 0 の不整合は実現可能性クランプが確定する」と注記)。q1e = min(q1, q0 r30)、q2e = min(q2, q0 r30²) で上限整合。
- 射影 (`condensationRealizability_d.cuh` / float 実体 `cond_realizability_clamp_f_d` in `condensationTransport_d.cu`): 無次元 (x, y) を可行域 0≤x≤1, x² ≤ y ≤ √x へ最近点射影、Q0 と g を保ち Q1・Q2 を書き戻す。

## 仮説 (主セッション)
H1: 射影と戦っているのは蒸発の一様 ṙ 形が数を除かない (S_Q0 = 0) こと。多分散分布を一様に左へずらすと最小の液滴が先に半径 0 に達するが、モーメント式はそれを除かないので、分布は「半径 0 への集積」側の境界 y = √x (2 点求積の一方の節点が r = 0) の外へ出る。二相拡散は高温で不飽和の壁近傍へ液を運び込みこの蒸発域を作る役 (OFF は液を混ぜないので起きにくい)。輸送そのものは境界の内向き。
H1 からの対処の候補 (#4g3n): 境界 y = √x 上の分布は w₀δ(0) + w₁δ(r₁) (r₁ = Q2/Q1、w₁ = Q1²/Q2、w₁r₁³ = Q2²/Q1 = Q3)。**半径 0 の数 w₀ を除く (Q0 ← Q1²/Q2、Q1・Q2・Q3 不変)** のが物理的に整合した補正。入れ方 (a) 上側境界の外では射影の代わりに Q0 を下げる、(b) 蒸発ソースに半径 0 での数の流出 S_Q0 < 0 を入れる (分布の再構成 — 2 点求積など — が要る)。

## 問い
1. H1 の解釈は観測から支持されるか。段分解の H_up の作り方 (ρ が斉次 0 次で消える、ρ_l の近似は張り付き判定にしか使っていない) に誤りはないか。ソース段の符号が「数を除かない」ことに起因すると言うには、他に何を確かめるべきか (例: 同じ節点で S_Q0 に半径 0 の数流出を仮に入れたら H_up の符号が変わるかの 0 step 試算)。
2. G3-a (2) の ON FAIL の扱い: 「射影が残差の外で釣り合っているので流束では閉じない」という読みでよいか。事前登録の合否はどう記録すべきか (FAIL のまま、既定化は保留継続)。ΔF を 3 本の最終場の平均差で代用したこと、w の判定不能の扱いは妥当か。
3. 対処の方向: (a) 射影を「上側境界の外では Q0 を下げる」補正に替える、(b) 蒸発ソースに数の流出を入れる、(c) その他 (文献の標準: Massot/Laurent/Kah らの EMSM・EQMOM の蒸発時の r = 0 数流束、Fox の QMOM 蒸発 等)。どれを採るべきか、既定経路 (OFF) の結果も変わる (OFF でも射影は少しは働く) ことをどう扱うか、この plan の範囲か別 plan か。
4. 次の一手として事前登録すべき A/B (0 step 試算・1 更新・延長 run) と合否。

読んでよいファイル: plan 上記、`notes/investigations/2026-10-04-twophase-g3/` (`g3a_judge.py`・`boundary_push.py`・`extract_projection_mask.py`・`evidence/`)、`solver_density_cuda/cuda_forge/condensationSource_d.cuh`・`condensationSourceKernels_d.cuh`・`condensationRealizability_d.cuh`・`condensationTransport_d.cu`、`methods/condensation.md`、`plans/accepted/condensation-evaporation.md`、`plans/accepted/condensation-source-limiter-steady.md`。
