# 諮問: g4 の初期場の作り方 (cross-mesh の最近傍補間はカウルの双子節点を取り違える) (2026-10-11)

関連 plan: `plans/active/tooling-sern-te-wake-grid.md` §5.1 #4・§6 #4 (g4 は「run_1079 の最終場から interp_field」と書いた。codex diagnose `2026-10-08-te-wake-ab-result-diagnose.md` の提案)。エスカレーション条件 4 (承認済みの手順からの変更)。

## 観測事実
- `solver_density_cuda/tools/interp_field.py` は KD-tree の最近傍 (3D は x,y,z) で原始量を移し、座標が一致する 2 節点は区別しない (双子の対策の引数は無い、grep で確認)。
- SERN の旧メッシャはカウル (厚さ 0、後縁の上流 i < i_te) の上下の面を**座標が同じ別 ID の節点 (dup1)**で持つ (`mesh_sern3d.py:271`)。g3 (run_1079 の格子) も g4 も同じ作り。
- 過去の事故 (メモリ interp-field-coincident-nodes-trap): case/46 run_0009 で、最近傍補間により排気側の壁節点が外部流の 2 kPa を持ち、2 次化した step 7 で発散。合成場で 134/134 station の誤写像を確認。
- g4 の格子仕様 (`problem_3d_prod_m6on_g4.yaml` との差): nj_ext_top 57→67、first_top_frac 6e-4→1.5e-4、nz_out 15→32、nj_top 73→83、nj_bot 55→65、first_wall_frac 1.6e-4→4e-5、`first_wall_frac_far: 4e-3`・`wall_frac_blend_len: 3.0` を追加 (x の station 数は同じ)。**g4 は first_wall_frac_far の帯の補間があり、te_wake の局所変形 (中間線の Hermite) と組み合わせるのは初めて** (投入条件の判定で確かめる)。
- 過去の g4 (run_0970_3d_g4_chidef → run_0972 cont40k、GATES PASS・床 0) は MOC の初期場から段階起動 (暖機 ramp → soft → mid → 本段) で作った。段階起動のレシピは R-b で検証済み (skill sern-eval)。
- B の格子の初期場の道具 `restart_field_deformed.py` は双子を検出して拒否する (動いた節点のみ対象、同じ接続の格子どうし用で g3→g4 には使えない)。

## 案
- (a) **g4 は MOC の初期場から生産の段階起動** (run_0970 と同じ流れ) → 本段で判定区間 20000 step (+20000 を 1 回まで)。g3 側は run_1079 の最終場からの継続。履歴は違うが、両側とも判定区間で STEADY を要求するので定常解の比較になる。
- (b) **双子を区別する補間**: interp_field に照会点を双対重心 (`CELLS/centCoords`、双子は上下それぞれの側に重心を持つ) にする選択肢を足し、src・dst とも重心で最近傍を取る (値は節点の値)。ツールの変更 + 試験 (双子の対応を合成格子で確認)。
- (c) interp_field のまま回し、双子の節点の値を事後に検査 (片側の値になっていたら作り直す)。

## 問い
1. どれにするか。(a) なら §6 #4 の比較の読み方 (初期場の履歴の違い) に注記が要るか。(b) なら合格条件。
2. g4 の 0 格子 (投入条件の基準用) は流れを回さない、のままでよいか。
