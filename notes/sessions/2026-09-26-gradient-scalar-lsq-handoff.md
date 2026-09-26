# 引き継ぎ: node スカラー勾配の LSQ 統一 (2026-09-26)

残作業の正本は plan [`plans/active/gradient-scalar-lsq-unification.md`](../../plans/active/gradient-scalar-lsq-unification.md) §5.1。ここは写しとポインタ。

## 場所

- ワークツリー `/home/sano/work/forge-sern-design`、ブランチ `feature/sern-design` (push 済み)。
- **未 commit**: `solver_density_cuda/` の hook 8 ファイル (#4a、既定 off・出力専用)。同じ diff を `notes/sessions/gradient-scalar-lsq-unification-4a-hook.patch` に保存。`solver_density_cuda/build` は未追跡なので `git add solver_density_cuda` は使わない (個別 add)。
- AWS: 共有インスタンス `i-0b1a5e0b8dc152f00` (IP は `aws_instance.sh status`、再起動で変わる)。認証 `AWS_SHARED_CREDENTIALS_FILE=/home/sano/.aws-wsl/credentials AWS_CONFIG_FILE=/home/sano/.aws-wsl/config AWS_PROFILE=forge`、鍵 `~/.ssh/test.pem`。手動 stop 禁止。`~/forge-pgrad-new` は 1a334275 + hook patch でビルド済み (`~/sglsq/forge_hook4a`)。旧 = `~/sglsq/forge_base_f67fe877`。スクラッチ run は `~/sglsq/`。
- 検証ハーネス: `case/09.Taylor-Green/_g0_lsq_seam/` (g_suite.py `s0`、s1_nointerference.py、pregather_check.py、s0y_species.py)。
- run は AWS で回す (ユーザ指示、2D も)。block は `FORGE_CUDA_BLOCKSIZE=128 FORGE_CUDA_BLOCKSIZE_SMALL=128`。

## ここまで

- Phase 1 opt-in 実装 (`mesh.scalarGradient: gg|lsq`、既定 gg) は commit 7ebf6d7d (レビュー済み)。
- S0 (a–e) PASS、S1(1) gg 旧 vs 新 PASS、#4b = **A** (tgv の NS 配列の不一致は周期 gather の atomicAdd 順序差)、#4c (5 種 Y) PASS。

## 未決 (次セッションの最初にやること)

上位判断 (Fable 上限時は codex、ユーザのルール変更) で決める:

1. **#4a hook の diff レビュー → commit 可否**。implementer の解釈: lsq の dY/dξ は wrapper 内で合併済みなので、main の pre-gather ダンプから外し wrapper 内で別 tag (`species_lsq.loop1`・`passive_lsq.loop1`) で書く。
2. **「ダンプ有効で res 不変」の判定規則** (事前に無かった)。S1 規則の流用では tgv (dUxdz 32 点・dUxdy 4 点、最大 6e-8) と case48 (res_1 roe 2 点、最大 0.0156) が FAIL。親の案: ダンプ有りでも「gather 前の局所配列と `FORGE_DUMP_MASSFLUX` の面流束がダンプ無しとビット一致」を決定的に示し、場は §6 S1 の訂正版で判定 (事後規則と明記)。
3. **§6 S1 (2) の訂正の適用** (#4b = A で前提は成立) と、S0/S1 を PASS として閉じて **S2 (物理 A/B、§6 表の 4 ケース + case/16 の 2 双子 + FCT smoke) に進むか**。

その後: S2/S3 → codex result 1 回目 → Phase 2 (既定化。前提は別 plan `tooling-stage-manifest-launch-binding` #2g — 未起票)。

## 関連

- 前提 plan (accepted): `plans/accepted/boundary-node-periodic-gradient-fix.md` (継続課題 §5.1 #7)。
- 回転周期 plan の #0a に node + type 1 周期の起動エラー化と NS gather の並進判定を登録済み。
- 後続項目 (§7): `periodicGather1ToRoot_d` の決定的 gather 化。
