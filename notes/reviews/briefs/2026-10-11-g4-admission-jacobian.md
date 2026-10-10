# 諮問: g4 の投入条件で「スケール済みヤコビアンの最小 ≥ 0.65」だけ FAIL — 基準の 0 格子も同値 (2026-10-11)

関連 plan: `plans/active/tooling-sern-te-wake-grid.md` §5.1 #4 (g4 の投入条件は g4 の 0 格子を基準に `--admission`)、起点 plan `convection-zero-thickness-edge-reconstruction.md` §6.0 の投入条件の表 (codex diagnose 2026-10-08 te-wake-blend-length で事前登録。0.65 は g3 の再生成値 0.67245 を見て決めた)。
事前の規則: 1 つでも外れたら投入しない (判定不能、L_b を変えて再挑戦しない — 諮り直す)。エスカレーション条件 3。記録の写し `notes/investigations/2026-10-11-g4-admission/`。

## 観測事実 (AWS-B、HEAD d734392d のメッシャ、現行の変換器 36b061b9)
- g4 の 0 格子 `case/46.sern_design/run_1080_tewake_g4_A0_m10` と 1.0 格子 `run_1081_tewake_g4_B10_m10` (YAML `problem_3d_prod_3op_wallres_lswx08_tewake_g4_{A0,B10}.yaml`、g4 の格子仕様は `problem_3d_prod_m6on_g4.yaml` と同じ: nj_ext_top 67・first_top_frac 1.5e-4・nz_out 32・nj_top 83・nj_bot 65・first_wall_frac 4e-5・first_wall_frac_far 4e-3・wall_frac_blend_len 3.0)。節点 3,221,994、ヘキサ 3,122,688。
- `--admission`: 前提・最初の下流の辺 (2.1994°)・後縁の折れ (下 2.1994°/上 1.6047°)・復帰区間の折れ (**5.9999°** ≤ 6°、位置 i=110・j=8・y = −0.33 m、A は 2.59°)・固定すべき量・層の Δy (最大 1.65 %)・非正のヤコビアン 0・skew 最大 0.70117 (A と同値、増加 0)・`check_mesh_quality` B SOFT-PASS = A SOFT-PASS・変形領域の新しい AR 超過 0・`check_dual_closure` PASS は**すべて PASS**。**FAIL は「スケール済みヤコビアンの最小 (B) 0.45235 ≥ 0.65」の 1 行だけで、A (0 格子) も 0.45235**。
- 位置: 0.65 未満のヘキサは 40,700 個 (全体の 1.3 %)、すべて x/H 1.5〜5、最小は x/H 3.1・y/H −4.65 (プルームの下の外部領域の遠方、後縁 x/H = 1.2 から遠く、変形の範囲 x/H ∈ (1.2, 2.2] の外)。g3 (run_1079) は最小 0.67245 (x/H 1.88・y/H −0.63)、0.65 未満 0 個。
- g4 の双子の割り当て (0 step の A/B、`twin_ic_check.py`): **B_OK** (領域別の初期化の誤割当 0、双子 11,620 組 [側壁 11,480・角 140]、壁節点 16,612 も 0、prepare が書いた初期場とも一致)、現行の補間 (A) なら 5,810 件の誤割当。
- 過去の g4 (run_0970 → run_0972、GATES PASS・床 0) も同じ g4 の格子仕様 (te_wake なし) で回っている。

## 問い
1. 0.65 の絶対閾値は g3 の値を見て決めたもので、g4 の基準の格子自体が満たさない。局所変形の投入条件としては「B が A より悪化しない」(最小値が A 以上、変形領域で 0.65 未満を新しく作らない) と読み替えてよいか。読み替えは結果を見た後の規則の変更になるので、認めるなら条件を書く。
2. 認めないなら、g4 の基準の格子の品質の扱い (g4 の格子仕様の見直し・別の水準) をどうするか。
3. 復帰区間の折れ 5.9999° が 6° の閾値の直下である点の扱い (A は 2.59°、g3 の B は 5.08°)。
