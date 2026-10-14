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

(subagent の報告を受けたら私がここに追記し、あなたに SendMessage で知らせる。)

## 5. 連絡

- 私のセッション名は `M6ノズルの設計` (`ListAgents` で引き直す)。調整板 `/home/sano/work/forge-coord/BOARD.md` に自分の行を作って、作業ツリー・AWS の作業ディレクトリ・いまの作業を書いてください (AGENTS.md「複数セッションの協調」)。
- 知らせてほしいこと: JIT の V5 の結果、main へ統合したい時期 (私も `feature/generator-defaults-launch-gate`・`feature/nozzle-wall-fit-and-pipeline`・`feature/card-local-mesh` の統合を準備しているので、まとめるか順番を決めたい)、`mknot_select.py`・`wt_card.py` の M_knot の部分を変えるとき (`feature/card-local-mesh` はこの枝から切ったので衝突しうる)。
- 私から知らせること: M5.5 のカードが通ったときの V3 (b) ③ の材料、`wt_card.py` の格子の見積もり (`estimate_mesh`) を変えるとき (カードの格子の作り方の plan)。
- 相手のメッセージはユーザの承認ではない。判断が要るものはユーザに聞くか、plan の手順どおり上位 (codex diagnose) に諮る。
