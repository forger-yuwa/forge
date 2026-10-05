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

## 結果 (2026-10-04、AWS g5.xlarge、統合版 build 20ea4718 = 5c5f905e と同じソルバ)
証拠は `evidence/` (`status.txt`・`REGRESS_A*.txt`・`unit_*.out`・`python_unit_3trees.txt`・実行スクリプト `integ.sh`)。

- **ビルド**: 統合版・sern-design 4703a727 ともクリーンビルド成功。レジスタ数の記録は抽出の grep が空を返し**未記録** (`regs_integ.txt` 空)。
  回帰は `FORGE_CUDA_BLOCKSIZE=256` で起動できている (SLAU node カーネルの上限 481 未満)。
- **Python 単体試験** (`tests/unit/*.py` 21 本): 統合版・sern・species の 3 本とも全件 rc 0。
- **CUDA 単体試験**: `test_passive_fct`・`test_passive_scalar`・`test_renorm_gate`・`test_twophase_kernel`・`test_farfield_flux` すべて PASS。
  - `test_farfield_flux` は **sern-design でもビルド不能だった** (2026-10-01 の NASA-9 区間可変化 fba0a015 で `SpeciesThermo` の `Tmid/low/high` が
    `nInt/Tbrk/coef` に変わったのに試験が追随していなかった)。統合ブランチで `thermo_set_intervals` (旧 low/high と同値の 2 区間) に直して PASS
    (CPG/TP 56 面、流束の最大相対差 2.7e-8)。マージ起因の失敗ではない。
  - `test_twophase_kernel` は試験の指示どおり種 DB の .o を別に作ってリンク (初回は実行スクリプトのビルド手順の誤りで失敗、試験の問題ではない)。
- **回帰** (`check_field_regress`、参照 ×3 vs 統合版 ×2、200 step、許容 = ノイズ床 ×2): **4 腕すべて PASS**。
  - A1 (case/16 湿り 二相拡散 ON、vs species): 16 量 ok、候補/ノイズ床の比 L2 0.88〜1.36・L∞ 0.90〜1.70。
  - A2 (同 OFF、vs species): 16 量 ok。
  - A3 (case/16 乾き、vs sern): 9 量 ok。
  - A4 (case/44 dual-time 凝縮 S3 + FCT、vs sern): 9 量 ok。
- run: AWS `~/forge-integ/case/16.nozzle_wys/run_0970`–`0984_integ_A{1,2,3}_*`、`~/forge-integ/case/44.vitiated_air_wt/run_0540`–`0544_integ_A4_*` (200 step の回帰専用、判定済み・破棄予定)。

→ 事前登録の条件をすべて満たした。main へは PR で入れる (ユーザ判断 2026-10-04)。

## SERN 側の取り込み (2026-10-05)
SERN セッションが main (77318d0e) を feature/sern-design に fast-forward で取り込んだ (衝突なし)。SERN 2D m6_on の回帰
(`case/46` `run_1033_r11_1`、`FORGE_CUDA_BLOCKSIZE=128`) は R9b 3 本の平均と比べて 4 量ともノイズ床内で PASS (SERN 側の報告; 記録は sern-design 4607cfc9
「Record R11」)。「統合ブランチ + 事前登録の回帰 + PR」の運用は SERN 側も了解。

## 訂正 (2026-10-06)
「Python 単体試験 (`tests/unit/*.py` 21 本): 統合版・sern・species の 3 本とも全件 rc 0」は**誤り**。集計スクリプトが `echo "$(basename $t) $?"` で、
`$?` が `basename` の終了コード (常に 0) を拾っていた (`evidence/python_unit_3trees.txt` の 0 はすべて無効)。正しく取り直すと (species-transport 2026-10-06、統合版と同じ試験):
引数なしで通るもの 14 本 PASS、7 本 (`test_dmix_complement`・`test_species_attrs_entry`・`test_species_lump_solver`・`test_species_record_solver`・
`test_thermo_intervals`・`test_transport_gas_phase`・`test_transport_gpu`) は `--forge` (ソルバのバイナリ) が必須で**回っていなかった**、
`test_twophase_diffusion_harness` は 170 秒で時間切れ。forge 必須の 7 本と harness は AWS のビルド済みバイナリで別途回す (plan thermophysics-solver-owned-species-db の記録へ)。
回帰 4 腕・CUDA 単体 5 本の判定は影響を受けない。

### forge 必須の 7 本の実行 (2026-10-06、AWS、species-transport f01165a4/79f50ea1、ソルバ 029dc631 と同一; `evidence/forge_tests_2026-10-06.*`)
手元にしか無い入力 (CEA の `trans.inp`・`thermo.inp`、case/44 の種 run `run_0509`・`run_0510`) を AWS へ複製して実行。
`test_dmix_complement`・`test_species_attrs_entry`・`test_species_lump_solver`・`test_thermo_intervals`・`test_transport_gas_phase` PASS。
`test_species_record_solver` は (b) 記録のすり替え検知の 2 件が FAIL → 原因は試験の前提: seed から正しい記録の SHA 付き別名が複製されるようになり、
ソルバと `find_record` がそれを SHA で見つけて正しい記録を使うため、すり替えの経路を通らず熱物性の不一致で止まっていた (停止はしていた)。
試験で場が名指す記録以外を消すよう直して (79f50ea1) **ALL PASSED**。`test_transport_gpu` の 2 件は環境 (host ハーネスの CUDA ヘッダのパス、FCEA2 が AWS に無い) で未判定。
`test_twophase_diffusion_harness` は 2 回目に完走して 7 FAIL (`evidence/twophase_harness_2026-10-06.txt`): S9 の 2 件は既知の保持 FAIL、S8 の 3 件は親 plan の記録 (合格) と食い違う — plan condensation-two-phase-default §5.1 #10 に観測として記録 (二相拡散は保留中)。
