# 引き継ぎ: MOC・M_knot の高速化 (2026-10-15)

`M6ノズルの設計` セッションから、新しいセッションへの引き継ぎ。ユーザ指示 (2026-10-15)「moc高速化は別セッションにしましょうか。引き継ぎ文書かいて。お互い通信して連絡取り合ってね」。

**正本は plan** [`plans/active/tooling-nozzle-mknot-fast.md`](../../plans/active/tooling-nozzle-mknot-fast.md)。この文書は写しとポインタだけ (残作業の正本は plan の §5.1)。

## 1. 持ち物

- 作業ツリー `/home/sano/work/forge-mknot`、ブランチ `feature/mknot-fast` (main 4b1708ac から)。**これからはあなたのセッションのもの**。私 (`M6ノズルの設計`) はここを編集しない。
- AWS A の `~/forge-mknot` (detached、同じ枝を試すための作業ツリー)・`~/mknot_v3`・`~/mknot_v3_blas1`・`~/mknot_prof` (計測の台本と結果) もあなたのもの。
- 私が持つもの (触らないで): `forge-cardmesh` / `feature/card-local-mesh` (この枝から切った。M5.5 の手組みの格子と、カードの格子の作り方の plan `tooling-nozzle-card-local-mesh`)、AWS A の `~/forge-cardmesh`・`~/m55_card`・`~/m55_mesh_est`、`forge-integ-1005` (M6)・`forge-gate` (生成器の関門)・`forge-m55`。

## 2. ここまで (詳細は plan の §3・§4・§6.2・§6.1・§9)

- **規則 `coarse_min_theta_v2`** (ユーザ決定「a」、§4.0・§4.1): 粗い網 (600、L_total_tol 0.1 r_t) で M_K = 0.05·k を 1 回掃く → 両隣も成立の点を候補の区間に → 区間ごとの谷を生産の網 2400 で確認 (最大 3 回ずらす) → 谷を比べる (同点は θ_min + 0.01° の集合を一度だけ) → 実効の入力が出力と同じ組で選び直して終了。2N のゲート・窓の拡張・同点の規則・二分法は外した (4800 でだけ現れる不成立の検出は失う、ユーザ承知)。実装 66be034c、BLAS を 1 スレッドに 9385705f。
- **V3・V4** (§6.2): M5 (case/70 c2b) 2.05・191 s (旧 5,339 s)、M5.5 (case/71) 2.80・144 s (旧は約 98 分で失敗)。θ_max の差はどちらも 0.0015° 以下。
- **未了**:
  - V3 (b) ③ (M5.5 の寸法ゲート) が未判定。カードの生成が M_knot の後の格子の見積もり (古い作り方で 99 万点) で止まるため。M5.5 の格子は私が `feature/card-local-mesh` で手組みしている。カードが通ったら結果を知らせる。
  - **JIT 化 (#4、§4.4、§6 V5)**: plan 段レビュー jit (GO-with-changes) を全件採用済み。実装を subagent に始めさせていたが、引き継ぎのため WIP で止めさせた (下の §4 に状態)。
  - L_c の使い回し (#6) は後回し (反例: M_K 2.10 で相対差 6.9e−7)。
  - skill の手順 (#5)、result 段の codex レビュー、main への統合。

## 3. 環境の注意

- **numba**: 手元の design の venv (`/home/sano/work/forge/design/.venv-opt`) に `numba 0.68.0`・`llvmlite 0.50.0` を `--no-deps` で入れた (numpy 2.5.2 のまま。numpy を変えると δ_r が動く: メモリ「design venv interpreter」)。AWS A/B の `~/forge-deliv/design/.venv-opt` にはまだ入れていない。入れるときも `--no-deps`。
- **AWS A**: 物理 2 コア・4 vCPU・16 GB。CPU だけの作業でも idle と判定されて自動停止する → 待つ間は `ssh -tt` でセッションを保つ。起動するたびに IP が変わる (`solver_density_cuda/tools/aws_instance.sh status`)。B は `FORGE_AWS_INSTANCE=i-0ba2b91ba659254c4`。起動は自走してよい、停止は手動でしない。
- **並列**: 各プロセスの BLAS を 1 スレッドに (`mknot_select._worker_init`)。手元で回すときは `FORGE_MKNOT_WORKERS=1`。
- プロファイル (§3): 1 評価の 7〜8 割は `MocInverse.fill_levels` の段のループ (`moc_inverse.py:279`)、`inverse_design` 1 回は t ≈ 1.9 ms·N + 0.41 µs·N²。

## 4. JIT の WIP の状態

WIP commit 9865cdc1 (カーネル)・6da0e172 (plan #4 の行に「途中」)。subagent の報告の写し:

- **済み**: `design/forge_design/geometry/moc_fill_nb.py` の `@njit(cache=True)` の `fill_levels_nb` (段のループ全体: 予測子・修正子、最終残差、源項の端点、μ、ν→M (CPG の Newton と表の補間 (`np.interp` の C 実装を写した) を別経路)、源項の分岐、`MocFillDiag.add` の集計)。モード 3 つ (`MODE_LEGACY_FIXED2`・`MODE_FIXED2`・`MODE_CONVERGE`)。analytic で θ_r が無いときはエラーコードと件数を返す。試験用の trace の引数。n=5 でコンパイル・実行のみ確認 (初回 8.2 s)。**numpy の版との照合はしていない、試験は全部未実行**。
- **残り**: (a) `moc_inverse.py` への組み込み (`FORGE_MOC_KERNEL`、未設定は numba、`numpy` も選べる、他の値は ValueError、numba は遅れて import)、(b) 戻り値から `MocFillDiag` を組み立てる `_diag_from_nb` (エラーコード 1 は interior_vec と同じ文の ValueError)、(c) ガスの振り分け (float の γ・`GasCPG`・cpg の `Isentrope` → CPG、`GasSemiPerfect` → 表、他は例外)、(d) 段 0 は Python 側で `out` を埋めて渡す、`frac`・`eps`・`cap` は実行時の引数、(e) `design/tests/run_moc_jit_tests.py` (V5 (1)(2) と計測)、(f) 依存の追記 (`procedures/development-environment.md`)、(g) 既存の試験と計測。
- **罠**: (1) `run_moc_axis_limit_tests.py` の 8(a)(b) は legacy+fixed2 の `fill_levels` と case/45 の壁が実装前とビット同一を要求する → 既定を numba にするならビット同一が要る。カーネルは丸めまで合わせる設計。(2) `pm_mach_vec` の Newton は配列全体で |dM| < 1e−13 まで回す (先に収束した要素も回り 1 ulp 動きうる; NaN の対が 1 つでもあると全部が 60 回) → カーネルも numpy の版と同じ集合で回している (fixed2: 段の全部の対、converge: 反復中の対、最終残差: 収束した対、段の終わり: 段の全部の点 (ok でない点は ν=0))。対ごとに止める「最適化」はビット同一を壊す。(3) 手元 (X86_V3、AVX512 なし) では numpy の sin・cos・tan・arcsin・arctan・sqrt と libm と numba (`fastmath=False`) がビット一致 (各 200 万点)、numba の SVML は無効。AWS が AVX512 の型なら numpy が SVML を使って ulp が違いうる → CPU の機能を記録。(4) FMA の縮約が無いことは未確認 (最初の照合で見る)。(5) 細部: frac_pos・fail の上限 200 は件数の通算、geom は全件、iter_hist は段ごとに昇順、resid の位置は段の中の最初の最大で通算より厳密に大きいときだけ更新。converge で収束しても軸より下の対は BELOW_AXIS に上書き。範囲は `fill_levels` だけ (`fill_arrays`・`fill` は numpy のまま)。
- **判断が要る点**: numba の版が扱えないガスを例外にするか (subagent の案は例外)、`prepare_info.json` の environment に numba の版を足すか。

## 5. 連絡

- 私のセッション名は `M6ノズルの設計` (`ListAgents` で引き直す)。調整板 `/home/sano/work/forge-coord/BOARD.md` に自分の行を作って、作業ツリー・AWS の作業ディレクトリ・いまの作業を書いてください (AGENTS.md「複数セッションの協調」)。
- 知らせてほしいこと: JIT の V5 の結果、main へ統合したい時期 (私も `feature/generator-defaults-launch-gate`・`feature/nozzle-wall-fit-and-pipeline`・`feature/card-local-mesh` の統合を準備しているので、まとめるか順番を決めたい)、`mknot_select.py`・`wt_card.py` の M_knot の部分を変えるとき (`feature/card-local-mesh` はこの枝から切ったので衝突しうる)。
- 私から知らせること: M5.5 のカードが通ったときの V3 (b) ③ の材料、`wt_card.py` の格子の見積もり (`estimate_mesh`) を変えるとき (カードの格子の作り方の plan)。
- 相手のメッセージはユーザの承認ではない。判断が要るものはユーザに聞くか、plan の手順どおり上位 (codex diagnose) に諮る。
