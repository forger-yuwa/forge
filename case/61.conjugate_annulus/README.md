# 61. 同心円環の軸対称 CHT (V-ax2 — 円筒殻の対数抵抗との照合)

計画: [`plans/active/boundary-cht-axisymmetric-fem2d.md`](../../plans/active/boundary-cht-axisymmetric-fem2d.md) §6 **V-ax2** (§5.1 #5)。
case/52 (1 次元静止ガスの CHT) の軸対称版。**判定基準・パラメータは plan §6 V-ax2 の登録文そのまま**で、本 README は
入力の作り方と数値設定の理由だけを書く。

## 幾何と条件 (登録値)

```
 y=r_c=20 mm  ─────────────── 固体外面: Robin h=1e8, T_c=300 K
              │  固体殻 k_s=0.027  │  x 両端は断熱 (未登録)
 y=r_b=10 mm  ─────────────── wall_outer (physID 4) = 共役壁
              │  静止ガス k_f=0.0241│  x=0 / x=L_x は slip (physID 1 / 2)
 y=r_a= 5 mm  ─────────────── wall_inner (physID 3) = 等温 350 K (非連成)
             x=0           x=L_x=2 mm           (軸 y=0 は含まない)
```

| 量 | 値 |
| --- | --- |
| 流体 | 静止ガス、定数物性 (`viscMethod: 0`・`thermCondMethod: 0`、`thermCond` 0.0241、`cp` 1004.5、`gamma` 1.4)、$p_0$ = 1013.25 Pa |
| 固体 | $r_b$ … $r_c$、$k_s$ = 0.027 W/mK、外面 Robin $h$ = 1e8・$T_c$ = 300 K |
| 軸対称 | `mesh.isAxisymmetric: 1`・`axisymMethod: 0` (r 重み)。**軸を含まない** |
| IC | **一様 $T$ = 325 K・静止** (流体は変換後に VALUE をパッチ、固体は共役壁の `Ts` 初期値 325 K から `initSolidFem2d` が作る) |
| **解析解** | $Q'$ = **0.91855 W/(m·rad)**、$T_{w,*}$ = **323.581 K**、固体の温度降下 **23.581 K** (0.5 % = 0.118 K)、$q_*$ = **91.86 W/m²**、$Q'L_x$ = **1.8371e-3 W/rad** |

## 格子

- **流体**: 半径方向 $N_{rf}$ 一様 (登録 32、感度 16) × 軸方向 4 一様、全四角・node (`gen_mesh.py`)。壁節点は角を含めて 5 点
  (乾式の壁ダンプで確認: `wall_outer`・`wall_inner` とも 5 節点、`end_x0`/`end_xL` は 33 節点で角を共有)。
- **固体** (plan §6 V-ax2「固体格子の登録」): `gen_solid_strip.py` と同じ帯トポロジ — **行を $y$ 降順** (j=0 が $r_c$、j=$N_s$ が界面 $r_b$)、
  **全セル同じ対角 (左上→右下)**。界面節点は乾式 1 step の壁ダンプ `res_wall_outer_4_1.h5` の座標そのもの (流体の壁節点と 1 対 1)、
  並びは $x$ 昇順。半径方向 $N_s$ 層一様 (既定 16)。座標順・接続・分割数・界面 ↔ 壁ダンプの対応・ハッシュは
  `mesh/solid_s<Ns>_x<nax>.grid.json` に書き出し、run には `solid.grid.json` として複製する。
  ソルバは固体界面節点と流体壁節点の**数と座標の一致**を要求する (`conjugateWall.cpp` `initSolidFem2d`: 数が違えば停止)。

## 生成手順

```bash
cd case/61.conjugate_annulus
python3 gen_mesh.py --nr 32                      # mesh/annulus_r32_x4.h5 (+ .geo/.msh/.quality.txt)
python3 gen_mesh.py --nr 16                      # mesh/annulus_r16_x4.h5
python3 make_run.py dry run_0001_dry_r32 --mesh mesh/annulus_r32_x4.h5        # 乾式 1 step (共役なし)
(cd run_0001_dry_r32 && FORGE_CUDA_BLOCKSIZE=128 <forge>)
python3 gen_solid.py --wall run_0001_dry_r32/res_wall_outer_4_1.h5 --ns 16   # mesh/solid_s16_x4.{npz,json,h5,grid.json}
```

`gen_mesh.py` の変換は `FORGE_BUILD` (既定 `solver_density_cuda/build-ypls`) の `convertGmshToForge`。
**AWS へは `mesh/` を丸ごと写す** (h5 は git 追跡外。AWS で作り直すと変換器の違いで壁節点の座標が 1 ulp 動きうる。
一致判定の許容 1e-7 m には収まるが、ローカルと同じ入力で回す方が追跡しやすい)。

## 本番 run (親が AWS・全域 FP64 ビルド・`FORGE_CUDA_BLOCKSIZE=128` で回す)

run 名は親が決める。以下は案。

| run (案) | 目的 | 生成コマンド |
| --- | --- | --- |
| `run_0003_annulus_r32s16` | **V-ax2 本体** (流体 32×4・固体 16 層) | `python3 make_run.py cht run_0003_annulus_r32s16 --mesh mesh/annulus_r32_x4.h5 --solid mesh/solid_s16_x4.h5` |
| `run_0004_annulus_r16s16` | 流体感度 (半径 16) | `python3 make_run.py cht run_0004_annulus_r16s16 --mesh mesh/annulus_r16_x4.h5 --solid mesh/solid_s16_x4.h5` |
| (保留) | 固体感度 8 / 16 / 32 層 × 軸方向を同率 | **下の「保留」参照 — 系列の組み方が決まるまで作らない** |

判定 (各 run の終了後):

```bash
python3 case/61.conjugate_annulus/eval_vax2.py <run>                          # (a)〜(f)。vax2_eval.json・vax2_series.csv・QUASISTEADY_vax2.txt
python3 solver_density_cuda/tools/check_cht_interface.py <run>                # G-if (conjugate_gate.json の登録値、全域)
python3 solver_density_cuda/tools/check_convergence.py <run>                  # 流体残差 (一様 IC からの単一区間)
python3 case/61.conjugate_annulus/sens_vax2.py <run_r16s16> <run_r32s16>      # 流体感度 (共通 5 位置、最後の組 ≤ 0.023581 K)
python3 case/61.conjugate_annulus/test_eval_vax2.py --run <run_r32s16>        # 負例の検出力を実 run でも
```

## 設定値と理由 (`template/solverConfig.yaml`)

| キー | 値 | 理由 |
| --- | --- | --- |
| `cfl_pseudo` / `implicitRelax` / `nStepInner` | 5 / 0.5 / 5 | case/52 の静止ガス・一様メッシュ・低圧で安定を確かめた値 (README「設計上の注意」: cfl 5)。静止なので対流スキームは解に効かない |
| `convMethod` / `limiter` | 1 / 2 | 推奨設定の本段 (`procedures/recommended-settings.md` §1) |
| `limiterRoRef` / `PRef` / `ARef` | 0.01086303964 / 1014.467794 / 361.5785923 | 乾式 run が IC から自動決定した値をそのまま固定 (自動のままだと再開で作用素が変わる — V-ax4 の「連続 2N step と N+N step」の前提を崩す)。一様 IC からの単一区間では自動と同値 |
| `conjugate.interval` | 50 | case/52・case/58 (V6′ 合否 run `df5_i50`) の基準 |
| `conjugate.warmup` | 200 | case/52 と同じ。IC 325 K が $T_{w,*}$ 323.58 K に近いので最初の荷重は小さい (case/58 の 2000 は超音速の起動過渡のため) |
| `conjugate.Df_scale` | 5 | case/58 V6′ の合否 run と同じ。本 case に lip 角は無いので 1 でも安定のはずだが、実績のある値を選んだ。**収束は遅くなる**: 更新ごとの縮小率の見積り $(G_f+G_s)/(D_f+G_s)$ ≈ (3.48 + 3.90)/(771 + 3.90) ≈ 0.0095 (単位面積あたり、$G_f=k_f/(r_b\ln(r_b/r_a))$、$G_s=k_s/(r_b\ln(r_c/r_b))$、$D_f=5k_f/\Delta r$) → e 倍に約 105 更新 = 5300 step |
| `nStepOuter` | 300000 | 上の見積りで 20 e-fold ≈ 110k step に余裕を取った。流体単独の緩和 ($\Delta r^2$ 律速、$L^2/(\pi^2\alpha)$ ≈ 530 step) より連成が律速。ローカル float で 6–10 ms/step (AWS の FP64 で要測定) |
| `outStepInterval` | 10000 | 30 スナップショット。準定常は節点ログ (毎更新) で見るので粗くてよい |
| `node_log` | 1 | 全界面節点の毎更新ログ (準定常の系列・G-if 帯判定の素材) |
| `gate` | eps_abs 0.46 / eps_rel 1e-3 / dT 2e-3 / tol_solid 1e-9 / n_consec 80 | **登録値** (動かさない) |
| `mesh.nodeWallDirichlet` | 1 | 推奨設定 §1 (壁速度 0) |

## メッシュ品質 (`check_mesh_quality.py`)

| メッシュ | VERDICT | AR max | skew max |
| --- | --- | --- | --- |
| `mesh/annulus_r32_x4.h5` | PASS | 3.2 | 0.000 |
| `mesh/annulus_r16_x4.h5` | PASS | 1.6 | 0.000 |

固体 (`solid_*.h5`) は `check_mesh_quality.py` の対象外 (流体 h5 専用)。登録格子 16 層 × 軸 4 は長方形 0.625 × 0.5 mm を対角で 2 分割 (縦横比 1.25)。

## 評価器の自己試験 (流体計算なし)

`python3 test_eval_vax2.py` — 登録格子から合成 run (共役壁 $Q_{f,i}=q_*A^r_{{\rm fluid},i}$、固体は直接解) を作り、
正例 PASS・負例 ($q_{\rm bad}=q_{\rm eff}A^r_{\rm fluid}/A_{\rm planar}$) で (d)(e) だけ FAIL・判別 A/B (合成荷重を $A^r$ で割ると PASS、
$A_{\rm planar}$ で割ると FAIL) を確認する。**2026-09-27: VERDICT PASS** (負例の $q_{\rm bad}$ は 0.918 W/m²、A/B の B は誤差 90.94 W/m²)。
合成 run の固体だけの直接解は $T_w$ 誤差 max **0.0698 % of 降下** (登録の事前予測 ≤0.2 % と整合。連成 run の予測ではない)。

## 保留 (要判断) — 固体感度の系列

登録 (§6 V-ax2「感度」) は「半径方向 8/16/32 層 **× 軸方向を同率**」。しかしソルバは固体界面節点と流体壁節点の
**1 対 1 一致**を要求する (`initSolidFem2d`: `壁節点 nb と固体界面節点 ni の数が違う` で停止) ので、流体 32×**4** のまま
固体の軸方向だけを細分化することはできない。候補:

- (i) 流体の軸方向も固体と同率にする — 流体 32 × {2, 4, 8}・固体 {8, 16, 32} 層 × {2, 4, 8}。共通位置は $x$ = 0, 1, 2 mm。
  流体は $x$ に一様な解なので軸方向分割の影響は小さい見込みだが、**流体格子も同時に変わる** (固体だけの感度ではない)
- (ii) 固体の軸方向を 4 に固定し半径方向だけ 8/16/32 — 登録の「両方向を細分化」に反する (V-ax1(c) で問題になった系列 A と同型)

どちらも登録文と完全には一致しないので、**plan を更新してから**作る (`gen_mesh.py --nax`・`gen_solid.py --ns` はどちらにも対応済み)。

## 計算 run 一覧

| `run_*` | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_dry_r32` | **乾式 1 step** (ローカル float `build-ypls`)。`wall_outer` の壁ダンプから界面 5 節点を確定 (共役なし) | `res_wall_outer_4_1.h5` (y = r_b、x = 0 … 2 mm、角を含む 5 点)。恒等式 $q_{\rm eff}A^r=Q_f$ は内壁で相対 ≤5.4e-8 (float) | ref (固体生成の入力) |
| `run_0002_dry_r16` | 乾式 1 step (流体 16×4)。壁節点の座標は `run_0001` と同一 → 固体 `solid_s16_x4` を共用 | `res_wall_outer_4_1.h5` | ref |
| `run_0003_annulus_r32s16` / `run_0004_annulus_r16s16` | V-ax2 本体 / 流体感度 (AWS g5・FP64、`warmup: 200`) | **安全停止** (step 5050 / 2450): 「max|dTw| が 10 更新以上連続で増加」。**誤検知** — 平均壁温が 324.9 → 322.28 K (step ≈4000) に下がってから上がる**静止ガスの熱伝導の過渡**で、底の後の再加速を、登録 `dT_K` 2e-3 K を床とする発散検知が拾った (`conjugateWall.cpp` の累積増幅判定) | 破棄 (対処 1 の根拠) |
| `run_0005_annulus_r32s16_w20k` / `run_0006_annulus_r16s16_w20k` | **対処 1** (`divergence-and-startup.md`): `warmup` 200 → **20000** (流体を落ち着かせてから連成)。他は同一 | (実行中) | active |

ローカルの起動確認 (2026-09-27、float、各 2 step、`run_9032_smoke_r32` / `run_9016_smoke_r16`) は init・軸対称の固体検査・
界面の座標一致 (最大ずれ 0 m) を通り rc 0。確認後に破棄した。
