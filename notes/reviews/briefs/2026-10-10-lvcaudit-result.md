# 諮問ブリーフ: 製品の経路の照合 (plan time_integration-line-viscous-jacobian-faceh §6.7) の結果の解釈と次の一手 (2026-10-10)

日付 2026-10-10。諮問先 codex (diagnose)。AGENTS.md の条件 3 (事前登録の比較が FAIL)・7 (result の解釈を確定する前)・1 (次の §6 を書く)。
plan: `plans/active/time_integration-line-viscous-jacobian-faceh.md` §6.7 (事前登録)・§6.8 (結果)。コード: `solver_density_cuda/cuda_forge/timeIntegration_d.cu` (監査の塊)、
`solver_density_cuda/tools/line_audit_helper.cpp`、`case/45.isobutane_m6_d155/lvcaudit_judge.py` (commit 732270b4・fe878c4b で採取の前に commit)。

## 観測事実

§6.8 の表のとおり。要点:
- VERDICT: `FAIL (最初に外れる段 C); 対照とのビット一致が成り立たないので、通常の経路への結論は保留`。S・G・B は不一致 0 (B は 47645 項目)。D は不一致 0・判別不能 3。
- C の不一致 125 はすべて対流の K の `K[3][3]` (2D の z 方向)。登録外の調査で、double に対して CUDA も host の float も数百 ulp ずれる (大きな項の相殺)。
- T: D が監査用と通常のビルドで最大 9.5e-7 (float の 1 ulp) 違う。Kprev・Knext・状態・dt・フラグ・スカラーの粘性の和は再実行どうしも監査用もビット一致。rhs は再実行でも揺れ、監査用との差はその範囲。dq は D の差の結果として違う。
- A: 壁法線のライン面 294 面で β・κ が double の幾何の係数と 1e-3 以上違う。最大 10.8 %、壁から 0〜15 番目の節点で 3〜11 %、52 番目まで 1e-3 超。
  原因の経路は ISP 0 で dcc を float にした座標の差で作ること (`timeIntegration_d.cu` の粘性の幾何)。例: 壁の隣の dcc が double 2.679e-8・float 3.003e-8 (y の ulp 7.5e-9)。
  B が合格しているので、製品は float の式どおりに計算していて、ずれは float の座標から来る (式の誤りではない)。
- この dcc は、値 0 (本線 B0) のスカラーの粘性の対角 2ν·δ/dcc にも同じ幾何で使われている。残差 (FP64 ビルドの `viscousFlux_d.cu`) は double の座標の差を使う。
- 既往: 値 3・マスク 7 は case/45 の方向別 dt で約 122 step で非有限 (§6.2)。共通関数は多倍長の参照で全列合格 (§6.6)。

## 期待値と出典

- §6.7 の分岐: B・C・D に不一致 → FAIL。T が成り立たない → 通常の経路への結論は保留。A は別に書き、再現の PASS と両立、破綻の原因の断定にも係数の生成の問題の除外にも使わない。
- 元の float 化の plan (`architecture-float-state-double-geometry`、元のセッションの担当) は、残差側の幾何を double の差のベクトルにする段を進めている (stage 1: `ge_x` 等)。LHS (陰解法の係数) の幾何はその対象か、私は確かめていない。

## 仮説 (呼び出し側の読み)

- H1: 登録の FAIL (段 C) は、対流の K の float の相殺に対して許容の尺度 (結果の列の大きさ) を誤って選んだことによるもので、製品の組立の誤りではない。
- H2: A の壁際の β・κ の最大 11 % のずれ (面ごとにばらつく) は、薄層の熱伝導・粘性の K と D を壁際で実残差と食い違わせる。値 3・マスク 7 の破綻に寄与しているかもしれない (未確認)。
- H3: T の D の 1 ulp の違いは、監査の書き込みがコンパイラの命令の組み方 (FMA の縮約など) を変えたもので、数値の性質は変えない。

## 問い

1. 登録の FAIL と、登録外の調査の記録の書き方。H1・H3 はどこまで言えるか。
2. A の発見 (壁際の LHS の係数が float の座標で最大 11 % ずれる) の意味と、次の一手。候補:
   (a) LHS の幾何の差だけを double の座標で取る診断の切替 (例: `dcc_x = ST(ccx[o] − ccx[ic])` を double で引いてから ST へ、既定はビット不変) を足し、値 3・マスク 7 で破綻が止まるかの A/B (cuda_forge の数値の変更、別に登録)。
   (b) `implicitSolvePrecision 1` (LHS 全体を double) で値 3・マスク 7 を回す A/B (コードは変えない、精度をまとめて変える)。
   (c) 先に元のセッションの float 化の plan に申し送り、LHS の幾何もその対象に入れてもらう。
   長い run を使わずに、筋のいい手法 (全部入り) に向けて判断に要る最小の組を示してほしい。
3. 値 0 (本線) のスカラーの対角も同じ幾何を使う。本線への影響をどう扱うべきか (元のセッションへの申し送りの文面)。
