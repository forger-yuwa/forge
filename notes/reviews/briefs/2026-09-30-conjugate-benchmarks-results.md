# 諮問ブリーフ: 共役ベンチマーク (A 厚肉管・C 共役平板) の結果の解釈 (2026-09-30)

AGENTS.md エスカレーション条件 **3** (登録基準の FAIL・判定不能)・**7** (結果の解釈を確定する前)。
plan `plans/active/boundary-cht-conjugate-benchmarks.md` §5.1 #6 に結果、§4.5〜§4.7・§6 に登録文。

## 読んでよいファイル
- plan の §4.5〜§4.7・§6・§5.1 #6 の行
- 各 run の `EVAL_CONJ.txt`・`CONVERGENCE_CHECK.txt`・`SERIES_CONJ.txt`・`CHT_INTERFACE_VERDICT.txt`・`CHT_BALANCE_VERDICT.txt`・`eval_conj_600000.csv`
  (case/64.conjugate_pipe_wall/run_0005〜0011_*、case/65.conjugate_flat_plate/run_0005〜0010_*)
- 評価器 `case/64.conjugate_pipe_wall/eval_conj.py`・`series_conj.py`

## 観測事実 (§5.1 #6 の行に詳細)
- A: 前提ゲート (収束・準定常全節点・G-if・G-cons) は 6 本とも PASS。主判定は r32・r64 が PASS (A1 温度 0.30/0.075 %、A2 0.09/0.02 %)、**r16 は A1・A2 とも判定不能** (A1: 温度 1.19 % FAIL、熱流束 U 1.1 % > 上限 0.67 %。A2: q 1.78 %+U 0.53 % FAIL、Q_up の U 超過)。温度差は格子倍増で約 1/4。
- A2 r16 の元の run (`run_0008`) は step 32300 で連成の安全停止 (max|dTw| 10 更新連続増 6.1e-4→1.3e-3 K)。発散手順の対処 1 として `Df_scale` 20 の `run_0011` を回し完走 (他の格子は 5)。
- C: 主判定は 6 本すべて PASS。**前提ゲートの流体収束 (rms_roe 1.4 桁で横ばい) と G-if (① res_abs 5.8 W/m² > 1.0、④ res_solid 1.9e-9 > 1e-9) が 6 本とも NOT CONVERGED**。準定常 (評価窓の全節点・積分量) は PASS。揺れは後縁のすぐ下流 (x 10.1〜10.7 mm) の slip 境界に集中、前縁上流の slip 節点に法線速度 0.72 m/s。既知の未修正欠陥 (node slip + 接線方向の密度勾配で偽の流れ) と同じ形。

## 登録の曖昧さ (私の誤り)
§6 は「forge の格子は各 3 水準」とだけ書き、**どの格子で主判定の合否を決めるかを書いていない**。

## 仮説
- H1: A の r16 は格子が粗いだけで、誤差は 2 次で減っている (1.19 → 0.30 → 0.075 %)。
- H2: C の前提ゲート FAIL は評価窓の外 (slip 境界) の揺れで、窓内の伝熱比較を損なっていない (窓内の界面温度の揺れ 2e-6 K)。

## 諮りたいこと (推奨を 1 つに絞る)
1. A の合否をどの格子で決めるべきか、r16 の判定不能をどう記録するか (登録の曖昧さを事後にどう扱うか)。
2. C の前提ゲート FAIL の扱い: (a) 判定不能として主張しない、(b) slip を使わない境界 (後縁の下流も断熱の no-slip 壁にする、上流 slip を無くす) で登録し直して回し直す、(c) 例外として受け入れ主張を限定、のどれか。
3. A2 r16 だけ Df_scale が違うことの扱い。
4. result レビューに進んでよいか、足りないものは何か。
