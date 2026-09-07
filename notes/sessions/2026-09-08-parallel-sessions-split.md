# 並行セッションの切り分け (2026-09-08)

同じ作業ツリー `/home/sano/work/forge` ・同じブランチ `feature/sern-design` で 2 つの Claude セッションが同時に作業している
ことが判明した (A の commit b12d74a4 の後に B の 1bd2f69d / c612b09e が同ブランチに入っていた)。衝突を避けるための取り決め。

| | セッション A (`session_01Qc4xBiaFFrVj1ZSKJmBnfE`, jsonl 59c62107) | セッション B (`session_01J2GHA3y8BzKnuMY5BnvVz2`, jsonl 4def020c) |
|---|---|---|
| 担当 | **case/16 Wyslouzil**: 2D/3D node NS・凝縮・AWS run・後処理・壁圧比較、変換器/ツール (`tools/`)、case README、Artifact ページ | **SST モデル整合**: `turbulence.sst*` の既定・整理 (F1 前処理, 等方項の順序)、plan `turbulence-sst-consistency-options` / `turbulence-sst-energy-includes-k` (R3: E_t = ρ(e+u²/2+k))、case/26 回帰 |
| 触るコード | `case/16.nozzle_wys/*`, `solver_density_cuda/tools/*`, `solver_density_cuda/mesh/*` (変換器), `input/calcWallDistance*` | `solver_density_cuda/cuda_forge/{ransSource,ransBoundary,ransTransport,scalarTransport,viscousFlux,turbulent_viscosity,dependentVariables}*`, `input/solverConfig.*` (turbulence キー), `main.cpp` |
| case/16 の run 番号 | `run_02xx` (次は run_0234) | `run_03xx` (run_0303〜 使用中) |
| バイナリ | `/home/sano/work/forge-bin-sessionA/{forge,convertGmshToForge}` を **FORGE_BIN で固定** (B のリビルドの影響を受けない; 更新は A が明示的にコピー) | `solver_density_cuda/build/` を自由にリビルド |
| AWS g5 (`~/forge`, `case/16.nozzle_wys`) | A が使用 (3D run)。B は使わない (使うなら本ノートに追記) | — |
| main へのマージ | **`/home/sano/work/forge-main` (worktree) で `git merge --no-ff feature/sern-design` → push**。共有ツリーで `git checkout main` しない | 同左 |

共通ルール:
1. commit は自分の担当ファイルだけを `git add` する (`git add -A` 禁止)。commit 前に `git pull --ff-only origin feature/sern-design`。
2. 共有文書 (`plans/README.md`, `procedures/solver-settings.md`, `methods/index.md`, `case/16 README` の run 表) は追記のみ・pull 直後に編集・すぐ commit。
3. 相手の plan (`plans/active/turbulence-sst-*` は B) は編集しない。A の知見は自分の plan/README/memory に書き、必要なら本ノートで相手に伝える。
4. 履歴書き換え (rebase/amend/force push) 禁止。作業ツリーの `git checkout <branch>` 禁止。
5. 検証 A/B で「旧バイナリ」が要るときは `FORGE_BIN=<path>` (run_case.sh 対応済, B 追加)。

状態 (2026-09-08 01:00): B が `sstOmegaProdFromPk` / `sstSigmaBlend` を既定 ON 化 (1bd2f69d, 回帰 §3.6 済) → A はこれを前提に case/16 の
今後の run を回す (case/16 の既存 run_02xx は既定 OFF 時代の結果)。A の未完: 3D SST 凝縮 (run_0229) の壁圧比較は済、`sstIsotropicStress` /
`sstEnergyKSource` は B の R3 plan で置換予定 (A は触らない)。
