# node 定常 run の残差床が起点 run の床に戻らない (2026-09-26 調査メモ)

起点: plan [`gradient-scalar-lsq-unification.md`](../../plans/active/gradient-scalar-lsq-unification.md) §5.1 #5a・#5g・#5l。
**このメモへの移管は、スカラー LSQ の Phase 1 GO を意味しない** (S2 の収束ゲート未達は plan 側に残る)。

## 観測

S2 の双子 (現行 `f99f236d`、gg・lsq とも) を、2026-09-12〜17 に作った収束 run の最終場から `restart_field.py` (保存量ビット一致) で
継続すると、次の case で**起点の残差床 (末尾 20 % 平均) の 1.5 倍以内に戻らない** (`check_floor_ratio`、修正後ツールでも同じ。
`case/09.Taylor-Green/_g0_lsq_seam/FLOOR_RECHECK.txt`):

| case | 起点 run | 現行 gg の末尾平均 / 起点床 (代表列) | 物理量 |
| --- | --- | --- | --- |
| case/48 冷却平板 | `run_0005_B_tw300_y3` (09-12) | rms_ro 1.85×・roUy 9.1×・roe 1.88× | ゲート内・STEADY |
| case/40 軸対称ノズル | `run_0045_node_yp1_outletfix_cont` (08-11) | rms_roOmega 56×・ro 1.6× | η_CF 0.9754 → 0.9779 (設定差: 廃止キー `nodeAxisDirichlet`) |
| case/16 0476 | `run_0476_tp5_node_euler_profile_base` (09-17) | 全列 21–64× | Y0 5.5e-5・Xi 4.8e-4 (旧 vs 現 gg) |
| case/16 0482 | `run_0482_passive_wys_s1_sfr2_c1` (09-17) | 通常判定 3.8 桁 PASS → 現行 plateau 2.5 桁 | 壁 p/p0 L∞ 0.16 % (旧 vs 現 gg) |
| case/39 周期丘 (**戻る**) | `run_0039_r1_gradfix_new_ext` (09-26、同日のバイナリ系列、`slauWallNormalChi: 0` 明示) | 0.99–1.00× | — |

## 切り分け済み

- **chi 単独ではない**: case/48 で `slauWallNormalChi: 0` だけ変えても 1.82×・9.08× (plan #5a = B)。
- **本 plan の commit 群 (7ebf6d7d..f99f236d) ではない (case/48 gg に限る)**: 実装直前版 `36d8ba03` でも同じ 24000 + restart + 24000 で
  1.83×・9.10×、現行 gg との床比 0.99–1.09 (plan #5g = A)。他 case・lsq には一般化しない。

## 残る候補 (未検証)

- 起点作成 (09-12〜09-17) 以降の他の commit (09-12 以降 `cuda_forge`・`input` に 96 commit)。旧起点のバイナリで同じ restart をすれば床に戻るかで判別できる。
- restart そのもの (SST の restart 非忠実性の既知の事例がある)。旧起点バイナリで restart して戻らなければこちら。
- 実施する場合は、判定規則を先に決めてから 1 case (case/48) で。
