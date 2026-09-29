# 諮問ブリーフ: A の result レビュー (NO-GO C0/M5/m1) の採否と、領域切断 (d) の適用 (2026-09-30)

AGENTS.md 条件 **5** (Major の採否)・**1** (§4.6 の設計判断)。レビュー `notes/reviews/2026-09-30-boundary-cht-conjugate-benchmarks-result.md`、plan `plans/active/boundary-cht-conjugate-benchmarks.md`。

## 採否案 (実装済み: commit 直前の HEAD。評価器 `case/64.conjugate_pipe_wall/eval_conj.py`・`series_conj.py`、ツール `check_quasisteady.py --abs-scale`)
- M1 採用: Q_up/Q_tot を forge は forge の総入熱 (固体ダンプの q_hole)、参照は各変種の Robin 入熱で割る。U は各変種の比そのものの差。
- M2 採用 → **ただし問題あり (下記)**: 登録どおり上流を 1.5 倍 (−120R) に延ばす変種を実装 (延長部は forge の入口列の ρ・u・v・T、圧力は入口勾配で線形)。
- M3 採用: series を `check_quasisteady.py --abs-scale 1` (新オプション: drift/fluct を登録尺度で割る。`test_gate_bad_input.py` PASS) に渡す。NaN・欠落・重複・尺度 0 を REFUSED (`test_series_conj.py` 5/5)。eval も期待節点集合と照合。
- M4 採用: 一次データを回収して台帳、評価条件別のファイル名 (`eval_conj_<st>_L<levels>.csv`)、ログと CSV を同時に再生成 (これから)。
- M5 採用: 固体の効果 (ΔT_s・Q_ax/Q_tot) の U を参照の変種から出し「効果 ≥ 5U かつ効果 > 許容幅」を判定、時系列を series に含めた。私の「許容より 1〜2 桁大きい」は A2 の ΔT_s について誤り (3.3 倍) — 訂正する。
- m6 採用: README・methods・plans/README を再評価後に同期。

## (d) の問題 — 実測
合成の正常ケース (参照解そのものを forge として書いた run、`test_eval_conj.py`) で、上流延長の変種を U に入れると**全長の温度の U が 2.46 %** (上限 0.33 %) になり判定不能。原因: 上流の断熱壁区間では**粘性散逸が熱を足し続ける** (壁温は x = −60R で入口から 0.14 K 上昇) ので、入口を −120R に移すと −80R に届くガスの温度が 0.1〜0.17 K 変わる。**これは参照解の数値誤差ではなく問題定義 (入口位置) の感度**で、forge と参照解は同じ入口位置・同じ BC の同じ問題を解いている。
(前回の plan レビュー M4 が (d) を要求した趣旨は「参照解の不確かさ」、result レビュー M2 は「登録どおり延長せよ」。延長は問題を変えてしまう。)

## 諮りたいこと (推奨を 1 つに絞る)
1. A の (d) をどうするか: (i) 延長部の散逸を止めた変種 (問題を変えずに領域だけ延ばす) で測る、(ii) (d) は「forge と参照が同一領域・同一 BC を解く」ので比較の不確かさに入れず、問題定義の感度として別に報告する (事後の設計変更として記録)、(iii) 登録どおり延長を U に入れて判定不能を受け入れる、のどれか。
2. 他の採否 (M1・M3〜M5・m6) に問題はないか。
