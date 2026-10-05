# lump を含む化学種拡散の縮約 (#7) — 検証記録

plan: [`plans/active/thermophysics-solver-owned-species-db.md`](../../../plans/active/thermophysics-solver-owned-species-db.md) §4.4 確定版・§5.1 #7・§6 V4 (事前登録 2026-10-05)。
設計の諮問: [`notes/reviews/2026-10-05-lump-blanc-diffusion-diagnose.md`](../../reviews/2026-10-05-lump-blanc-diffusion-diagnose.md)。

## 7a 参照実装 (double、ソルバ非依存) — `evidence/reference_7a.txt`
`solver_density_cuda/tests/unit/test_lump_diffusion_reduction.py`: V4a (合成例・va3) PASS、V4d (SERN 15 状態で ε_新 ≤ ε_旧) PASS。

## 7b 実装の検証 (AWS g5.xlarge、build 91a79925、旧 = 20ea4718 (main 相当))
実行スクリプト `evidence/v7.sh`、ログ `evidence/status.txt`。

| 判定 | 結果 |
|---|---|
| host 単体 `test_lump_diffusion_host.cpp` (ローカル) | ALL PASS (V4b va3 の H2O が full config の現行関数と 16/16 ビット一致、V4e float vs double ≤2.6e-7、V4a 恒等式 2.4e-16)。AWS では `vector_types.h` が無くビルド不可 (ローカルのみ) |
| V4f GPU の D vs 独立 double 参照 (SERN 2D の場、62067 セル、11 実種) | **PASS** 最大相対差 4.8e-7、非有限 0、縮約フラグと記録の lump 有無が一致 (`evidence/V4f_sern.txt`) |
| V4b lump 無しの不変 (case/16、旧 ×3 vs 新 ×2、200 step) | 乾き・湿り OFF・湿り ON の 3 腕とも `check_field_regress` **PASS** (`evidence/REGRESS_*.txt`) |
| 性能 (SERN 2D、10000 step 継続 ×2 本ずつ) | 旧 1.93/1.94 ms/step → 新 2.25/2.22 ms/step (+16 %)。**後の GPU 占有測定で +9〜10 % と訂正** (下の 7e; この組は GPU 共有の混入を含む) |
| V4g 記録: SERN 2D の場の変化 (同じ起点から 10000 step、旧 ×2 vs 新 ×2) | 流れの量 (ρ・ρU・ρe・P・T・k・ω・h0・μt) は run 間のノイズ床内 (比 ≤1.38)。**化学種の質量分率は変わる**: Y の差は L2 でノイズ床の 9.8 倍・L∞ 31.8 倍 (L∞ 5.8e-3) (`evidence/SERN_change*.txt`)。SERN の目的量 (力係数) への影響は SERN 側の回帰で確認する |

レジスタ数: AWS に `cuobjdump` が無く記録できていない (これまでの `regs_*.txt` が空だった原因)。ビルド時の `-Xptxas -v` で取る必要がある。回帰は `FORGE_CUDA_BLOCKSIZE=128`/`256` で起動できている。

run (AWS `~/forge-integ`): case/46 `run_1900_lumpdiff_probe_new`・`run_1901`–`1904_lumpdiff_cont_{old,new}{1,2}`、case/16 `run_0985`–`0999_lumpdiff_{dry,wetoff,weton}_*`。
初回投入は AWS の作業ツリーの未 commit の変更で checkout が中断し旧コミットをビルドしたため無効 (削除し同名で再実行)、2 回目は `__constant__` の定義位置でビルドエラー (91a79925 で修正)。

## 7c 表引き (2026-10-05) — 撤回
二元係数 D_rq·P を Neufeld のクランプ点で分けた ln T の 3 次 Hermite で引く版 (36cadb3b)。host の V4e は二元 1.9e-7・縮約 2.9e-7 で合格
(初版は ln T を float で作り ln f を float で持ったため 2.13e-6 で基準 2e-6 を超えた; 基準は変えず精度を直した)。AWS (`evidence/v7c_*`):
V4f PASS、表 vs 式の場 (10000 step、2 本ずつ) はノイズ床内。**1 step の時間 (SERN 2D)**: 旧 1.93、式 2.27、表 2.30 ms
(1 回目の組 8.36/6.32/2.33 ms は他セッションと GPU を共有した時間帯で使わない)。**表引きでは速くならない** → 外して式のまま (ソルバは 91a79925 と同一)。
run: case/46 `run_1905_lumpdiff_tab_probe`・`run_1906`–`1911_lumpdiff_tabperf_{old,form,tab}{1,2}` (AWS `~/forge-integ`、破棄予定)。

## 7e 局所配列なし・GPU 占有時の性能 (2026-10-06)
実種ごとにスカラーで積む形 (029dc631) は GPU の probe 出力が配列版とバイト一致 (`run_1912_lumpdiff_noarr_probe`)。
性能は GPU に他プロセスがいないときだけ測り、測定中に出た回は捨てて測り直した (`evidence/v7p.sh`、`evidence/v7e_perf_clean.txt`、`run_1919_lumpdiff_perf_clean`):
旧 1.96/1.94、配列版 2.13/2.13、配列なし 2.14/2.12 ms/step → **+9〜10 %**。`run_1913`–`1918` は他セッションと GPU を共有していて無効 (`evidence/v7e_status_contaminated.txt`)。
