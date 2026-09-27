# 諮問ブリーフ: 軸対称 fem2d の V-ax0〜4 の解釈を確定してよいか (2026-09-27)

AGENTS.md 条件 **7** (result 段の解釈を確定する前、codex result レビューの前)。

## 読んでよいファイル
- plan `plans/active/boundary-cht-axisymmetric-fem2d.md` の §6 (V-ax0〜4 の登録と、そこに追記した結果・再登録)、§5.1 #1〜#12 (特に #12 総括)
- 流体側 plan `plans/active/axisymmetric-graded-grid-static-gas.md` の §5.1 #2〜#4
- case README: `case/61.conjugate_annulus/README.md`・`case/62.conjugate_disk/README.md`・`case/58.conjugate_slot/README.md` (run 一覧)
- 実装差分: `git diff 1a389fef..HEAD -- solver_density_cuda/conjugateWall.cpp solver_density_cuda/conjugate/ solver_density_cuda/tools/check_cht_balance.py solver_density_cuda/tools/cht_loop.py`
- 評価器: `case/61.conjugate_annulus/eval_vax2.py`・`sens_vax2.py`、`case/62.conjugate_disk/eval_vax2b.py`・`series_vax2b.py`、`solver_density_cuda/tools/test_solid_fem2d_axisym.py`
- run の成果物はローカルにある (大きな h5 は必要な配列だけ numpy で)

## 総括 (plan §5.1 #12 の写し)
- V-ax0 PASS (ローカル負例 5 件・受理 1 件・method 1 の NaN)
- V-ax1: (a)(b)(c‴)(d) PASS。(e) 収支の (c‴) 格子で 2.58e-12 > 1e-12 → **FAIL を例外受け入れ (ユーザ決定)** (零和性の丸め × 絶対温度の FP64 演算床)。旧 (c) の FAIL (最初の格子・同方向対角) も履歴に残す
- V-ax2 PASS: 本体 `run_0005` 壁温 0.0375 % of 降下・G-if・残差 (一様 IC 単一区間)・準定常。連成格子感度 0.00655 K、交差の和 0.00656 K、流体半径 0.00176 K (許容 0.023581)
- V-ax2b PASS: 一様 N_r 8/16/32 (一様 IC)、比 0.317/0.298、**非一様 N_r=32 は純伝導 IC に再登録**して 0.0858 % (一様 IC からは流体側の偽流れで失敗 — 流体側 plan で未解決)
- V-ax3 PASS (新旧差 1.18e-3 K ≤ 2ε₀ 3.96e-3 K、平面 case/58 10000 step)
- V-ax4 PASS (再開 1.9e-8 ≤ 2ε₀ 2.97e-8 K、起動時の状態照合、評価器の REFUSED)
- 観測: 最粗の感度 run (軸 2) で `iface_q_eff` 0.76 % (本体基準 0.5 % を超えるが感度 run の合否項目でない)

## 諮りたいこと (推奨を 1 つに)
1. この総括で**拒否を解除し plan を `done` にしてよい**と書けるか (PASS / 条件つき / 判定不能)。書くなら**主張の範囲の 1 文**。
2. 見落とし (例外受け入れ・再登録・保証範囲の書き方、methods/boundary.md の「計画中」の更新要否、流体側 plan の未解決の扱い)。
3. 次の codex result レビューに回してよいか。
