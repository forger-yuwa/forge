# main 統合 (2026-10-04): feature/sern-design + feature/species-transport

統合ブランチ `integrate/main-2026-10-04` = origin/main (bcc68797) ← feature/sern-design (4703a727) ← feature/species-transport (8cc24a5b)。
main は 2026-09-18 以降 sern-design の祖先で独自の内容を持たない (`git diff sern...main` 空)。species-transport は 2026-10-03 (653dfd7e) に
sern-design へ一度取り込まれており、差は species 側の 36 commit (二相拡散既定化 plan の S0・G1〜G3 診断。いずれも既定 OFF / 環境変数で有効)。

## 衝突の解消 (3 件、すべて両側の追加)
- `cuda_forge/passiveKernels_d.cuh` `species_advection_faceY_d`: 引数を `ffY` (sern: farfield 面組成) と `tpo` (species: 診断 G3-a) の両方に。
- `cuda_forge/speciesTransport_d.cu`: 化学種の呼び出しは両方を渡す。受動種の S3 呼び出し (自動マージ側) は `ffY = nullptr` を補う (sern でも受動種には渡していなかった)。
- `tools/stage_manifest.py`: hard キーに `conjugate.solid_sha1` (sern) と `twophase_diffusion_effective` (species) の両方。署名は sern の `stage_key(cfg, bc, run_dir=None)`。

## 事前登録 (run 前に commit)
Python: `test_stage_manifest_scalar_gradient.py`・`test_gate_bad_input.py`・`tests/unit/test_twophase_state.py` PASS (済)。`check_plans.py` の FAIL が親 2 本 (28 / 29 件) から増えない (済: 28)。
`tests/unit/*.py` は親ブランチで PASS のものが統合版でも PASS。

AWS (g5.xlarge): 統合版をクリーンビルド (ビルド成功、`cuobjdump` でレジスタ数を記録)。CUDA 単体試験 `test_passive_fct`・`test_passive_scalar`・
`test_farfield_flux`・`test_renorm_gate`・`test_twophase_kernel` がビルドでき、親での結果と同じ判定。

回帰 (`check_field_regress.py`、参照 ×3 vs 統合版 ×2、200 step、`FORGE_CUDA_BLOCKSIZE=256`): すべて VERDICT PASS を合格とする。

| 腕 | 入力 | 参照バイナリ | 設定 | 何を見るか |
|---|---|---|---|---|
| A1 | case/16 湿り凝縮 二相拡散 ON (G2 生産 config run_0561, IC run_0482 res_48000) | species 8a9b673b | 両方に `space.slauWallNormalChi: 0`・`mesh.scalarGradient: gg` を明記 (sern の既定変更 2 件を外す) | species の既定 OFF 経路・二相拡散が統合で変わらない |
| A2 | 同 OFF (run_0567) | species 8a9b673b | 同上 | 同上 |
| A3 | case/16 乾き (run_0524 config + res_48000) | sern-design 4703a727 (クリーンビルド) | 既定のまま | sern の既定経路が species の取り込みで変わらない |
| A4 | case/44 dual-time Euler 凝縮 (G0 c44 入力) | sern-design 4703a727 | 既定のまま | 受動種 S3 移流 (衝突解消箇所) と FCT |

FAIL のときは main へ入れない (原因を切り分ける)。
