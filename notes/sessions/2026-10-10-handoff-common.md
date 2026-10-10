# 引き継ぎの共通ルール (2026-10-10、ライン陰解法の速度・float 化のセッションから)

以下の個別メモ (`2026-10-10-handoff-*.md`) を受け取ったセッションが守ること。残作業の正本は各 plan の §5.1 で、このメモは写しとポインタ。

- **作業ツリーとブランチを分ける。**
  - 元のセッションは `/home/sano/work/forge-integ-1005` (ブランチ `feature/nozzle-wall-fit-and-pipeline`) で float 化を実装中。ここで編集・`git add -A`・`commit -a` をしない。
  - 自分の作業ツリーを作る: `git -C /home/sano/work/forge-integ-1005 worktree add /home/sano/work/forge-<slug> -b feature/<slug> origin/feature/nozzle-wall-fit-and-pipeline`。
  - commit はパス指定 (`git commit -- <paths>`)。自分のブランチに push する。元のブランチへのマージは区切りで相談する。
- **触らない plan**: `plans/active/time_integration-line-implicit-speed.md` と `plans/active/architecture-float-state-double-geometry.md` (元のセッションが編集中)。結果は自分の担当の plan か notes に書き、元のセッションが速度 plan へ反映する。
- **case/45 の run 番号**: 元のセッションは 03xx を使う。引き継いだ側は個別メモで割り当てた帯を使う。case README の表は、自分の run の行だけを編集の直前に読み直して足す。
- **AWS** (`solver_density_cuda/tools/aws_instance.sh status`、IP は起動のたびに変わる。鍵 `~/.ssh/test.pem`)。
  - インスタンスは共有。起動は許可済み (この作業の範囲)。手動で止めない (30 分の空きで自動停止)。
  - 元のセッションは、速さを測るときに GPU を専有したい。自分の run の前に `nvidia-smi` と `pgrep -x forge` で他の run を確かめる。
  - 判定が済んだ run は、途中の全場 `res_<n>.h5` と `nozzle.h5` をその日のうちに消す (2026-10-10 にディスクが 99 % になった)。
  - 手順と罠は skill `forge-aws-run`。
- **AWS のバイナリ**: `~/forge-linespeed-fp64/solver_density_cuda/build/forge` (= `cold_cfl.py` の `lineM_fp64`、LAYOUT2 既定、FP64) を `FORGE_BIN` で固定して使ってよい。
  - このツリー、`~/forge-linespeed-f32`、`~/forge-fgeom-*` は元のセッションの比較の基準なので、**再ビルドしない**。別のバイナリが要るときは、別のディレクトリにビルドする。
- 判定・報告は AGENTS.md のとおり (VERDICT の貼付、事前登録、codex の plan/result 段、上位への諮問の条件)。
