# 62. 同軸円板の軸対称 CHT (V-ax2b — 半径が変わる界面の局所整合)

計画: [`plans/active/boundary-cht-axisymmetric-fem2d.md`](../../plans/active/boundary-cht-axisymmetric-fem2d.md) §6 **V-ax2b** (§5.1 #6)。
界面・Robin 辺の半径が 4 倍変わるので、plan §3 の「流体の半辺積分と固体の集中量 $\int N_ir\,ds$ の違い」が局所の壁温の歪みとして
出るならここで出る。**判定基準・パラメータは登録文そのまま**。**FAIL したら荷重の再配分を入れずに止めて諮る** (plan)。
共通部品は [`../61.conjugate_annulus/axcht.py`](../61.conjugate_annulus/axcht.py) (case/61 と同じ作法)。

## 幾何と条件 (登録値)

```
 r=r_2=20 mm ┌──────────┬──────────┐   side_r2 (physID 2, slip) / 固体の r 端は断熱
             │ 静止ガス  │  固体     │
 wall_hot    │ k_f      │  k_s      │  Robin h=1e8, T_c=300 K (x = H+t)
 350 K (3)   │ 0.0241   │  0.027    │
 r=r_1= 5 mm └──────────┴──────────┘   side_r1 (physID 1, slip)     軸は含まない
            x=0     x=H=5 mm     x=H+t=10 mm
                    wall_cj (physID 4) = 共役壁
```

| 量 | 値 |
| --- | --- |
| 流体・IC・数値設定 | case/61 と同じ (静止ガス、定数物性、$p_0$ 1013.25 Pa、一様 325 K・静止、固体も 325 K) |
| **解析解** ($x$ の 1 次関数で $r$ に依らない) | $q=(350-300)/(H/k_f+t/k_s+1/h)$ = **127.339 W/m²**、$T_{w,*}=T_c+q(t/k_s+1/h)$ = **323.581 K**、降下 **23.581 K**、$Q_{\rm tot}=q(r_2^2-r_1^2)/2$ = **2.3876e-2 W/rad** |

(登録文の式は $1/h$ を省いているが、$h$ = 1e8 の寄与は $T_{w,*}$ で 1e-6 K 程度。評価器は含めている)

## 格子 (事前固定、plan レビュー M5)

- **流体**: $x$ 方向 16 一様 × $r$ 方向 $N_r$ = 8 / 16 / 32 一様、および $N_r$ = 32 の等比 1.1 (**間隔は $r$ が大きいほど広い**: 0.0746 → 1.431 mm、`gen_mesh.py` が向きを検査)
- **固体**: $x$ 方向 16 層一様 × $r$ 方向 $N_r$ (界面節点 = 乾式の壁ダンプ `res_wall_cj_4_1.h5` の座標、流体の壁節点と 1 対 1)。
  **四角形は左下→右上の対角**。`gen_solid_strip.py` と同じ帯トポロジ (行 j=0 が背面 $x=H+t$、j=16 が界面、界面の並びは **$r$ 降順**) で
  組むと対角 (j,i)–(j+1,i+1) は $x$ も $r$ も減る向き = 左下↔右上になる (case/58 は固体が −x 側なので見かけが鏡映)。`gen_solid.py` が向きを検査する
- 格子の記述は `mesh/solid_disk_<tag>.grid.json` (run には `solid.grid.json`)

## 生成手順

```bash
cd case/62.conjugate_disk
python3 gen_mesh.py --nr 8;  python3 gen_mesh.py --nr 16;  python3 gen_mesh.py --nr 32
python3 gen_mesh.py --nr 32 --ratio 1.1                         # mesh/disk_r32_g1p1.h5
python3 gen_mesh.py --nr 16 --axis-probe                        # 軸負例 mesh/disk_axisprobe_r16.h5
python3 make_run.py dry run_0001_dry_r8u --mesh mesh/disk_r8_u.h5      # 以下 r16u / r32u / r32g も同様
(cd run_0001_dry_r8u && FORGE_CUDA_BLOCKSIZE=128 <forge>)
python3 gen_solid.py --wall run_0001_dry_r8u/res_wall_cj_4_1.h5        # mesh/solid_disk_r8_u.{npz,json,h5,grid.json}
```

**AWS へは `mesh/` を丸ごと写す** (case/61 と同じ理由)。

## 本番 run (親が AWS・全域 FP64・`FORGE_CUDA_BLOCKSIZE=128`。run 名は案)

| run (案) | 目的 | 生成コマンド |
| --- | --- | --- |
| `run_0006_disk_r8u` | 系列 $N_r$=8 | `python3 make_run.py cht run_0006_disk_r8u --mesh mesh/disk_r8_u.h5 --solid mesh/solid_disk_r8_u.h5` |
| `run_0007_disk_r16u` | 系列 $N_r$=16 | `python3 make_run.py cht run_0007_disk_r16u --mesh mesh/disk_r16_u.h5 --solid mesh/solid_disk_r16_u.h5` |
| `run_0008_disk_r32u` | **合否** $N_r$=32 一様 | `python3 make_run.py cht run_0008_disk_r32u --mesh mesh/disk_r32_u.h5 --solid mesh/solid_disk_r32_u.h5` |
| `run_0009_disk_r32g` | **合否** $N_r$=32 等比 1.1 | `python3 make_run.py cht run_0009_disk_r32g --mesh mesh/disk_r32_g1p1.h5 --solid mesh/solid_disk_r32_g1p1.h5` |
| `run_0010_axisprobe_fp64` | 軸負例 1 step (恒等式が `axisRFloor` の床で成り立つか。FP64 で取り直す) | `python3 make_run.py axisprobe run_0010_axisprobe_fp64 --mesh mesh/disk_axisprobe_r16.h5` |

判定:

```bash
python3 case/62.conjugate_disk/eval_vax2b.py <run>                 # (a) T_w・(b) G-cons・(c) 連成保存・(d) 恒等式・(e) 準定常
python3 solver_density_cuda/tools/check_cht_interface.py <run>     # G-if (V-ax2 と同じ登録値)
python3 solver_density_cuda/tools/check_convergence.py <run>
python3 case/62.conjugate_disk/series_vax2b.py <r8u> <r16u> <r32u> --nonuniform <r32g>   # N_r=32 の 0.5 % と誤差比 ≤ 1/2.5
python3 case/62.conjugate_disk/eval_vax2b.py <axisprobe_run> --identity-only
```

## 設定値

case/61 と同じ (`template/solverConfig.yaml`、理由は [case/61 README](../61.conjugate_annulus/README.md)「設定値と理由」)。差:

- リミッタ基準値は `run_0001_dry_r8u` の自動値 (0.01086303964 / 1015.685644 / 361.7917356) を系列 4 本と軸負例で共通に固定
- 連成の収束の見積り: $D_f=5k_f/\Delta x$ = 386、$G_f=k_f/H$ = 4.82、$G_s=k_s/t$ = 5.40 [W/m²K] → 更新ごと約 0.026、e 倍に約 38 更新 = 1900 step。`nStepOuter` 300000 は case/61 と揃えた (過大だが同じ扱いにする)
- G-if は V-ax2 と**同じ数値** (eps_abs 0.46 / eps_rel 1e-3 / dT 2e-3 / tol_solid 1e-9 / n_consec 80) を使う。登録文「V-ax2 と同じ G-if」を数値の一致と読んだ。
  円板の $q$ = 127.3 W/m² に 0.5 % を当て直すと eps_abs は 0.64 なので、**0.46 はそれより厳しい**
- 軸負例: `mesh.axisRFloor: 5.0e-4`。軸節点の壁半割面の重心は $r=\Delta r/4$ = 0.3125 mm で、既定の床 (1e-20 m) では床に当たる面が無いので上げた (次の面は 1.25 mm)。
  両壁とも非連成の等温 (350 / 300 K) にして 1 step 目から $Q_f\ne0$ にした (界面の軸上節点は CHT ではソルバが拒否するので、恒等式の負例は `interfaceDiag` 単独で取る)

## メッシュ品質 (`check_mesh_quality.py`)

| メッシュ | VERDICT | AR max | skew max |
| --- | --- | --- | --- |
| `mesh/disk_r8_u.h5` | PASS | 6.0 | 0.000 |
| `mesh/disk_r16_u.h5` | PASS | 3.0 | 0.000 |
| `mesh/disk_r32_u.h5` | PASS | 1.5 | 0.000 |
| `mesh/disk_r32_g1p1.h5` | PASS | 4.6 | 0.000 |
| `mesh/disk_axisprobe_r16.h5` | PASS | 4.0 | 0.000 |

固体は対象外 (流体 h5 専用)。固体の長方形は $\Delta x$ 0.3125 mm × $\Delta r$ (流体と同じ)。

## 評価器の自己試験 (流体計算なし)

`python3 test_eval_vax2b.py` — 合成 run (流体半辺の厳密な面積分 $Q_{f,i}=qA^r_{{\rm fluid},i}$ を固体円板に与えて直接解) で
**2026-09-27: VERDICT PASS**。界面温度の最大誤差は $N_r$ = 8 / 16 / 32 一様で **0.94206 / 0.29831 / 0.08796 % of 降下** —
plan の参考予測 (codex の独立計算) 0.9421 / 0.2983 / 0.08796 % と比 0.99995–1.00002 で一致。非一様は **0.16168 %** ($r_2$ 端)。
**これは荷重移送だけの試験で、連成計算の予測ではない** (登録文どおり)。判別 A/B (合成荷重 $q_*A_i^r$ を $A^r$ で割る / $A_{\rm planar}$ で割る) は
固体側・流体側とも A PASS・B FAIL。

## 計算 run 一覧

| `run_*` | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_dry_r8u` | 乾式 1 step (ローカル float)、流体 16×8 | `res_wall_cj_4_1.h5` (x = H、9 節点) | ref (固体生成の入力) |
| `run_0002_dry_r16u` | 同、16×16 | 17 節点 | ref |
| `run_0003_dry_r32u` | 同、16×32 | 33 節点 | ref |
| `run_0004_dry_r32g` | 同、16×32 等比 1.1 | 33 節点 | ref |
| `run_0005_axisprobe_r16` | 軸負例 1 step (**ローカル float = 恒等式は判定不能**、経路確認のみ) | 床に当たった面 2 (両壁の軸節点)。恒等式 相対 ≤6.1e-8 (float の丸め)、床を外すと 0.375 (検出力あり)。`eval_vax2b.py --identity-only` は REFUSED (FP32) | active (経路確認) |
| `run_0006_disk_r8u` 〜 `run_0009_disk_r32g` | V-ax2b 系列 (AWS FP64、`warmup: 200`) | r8/r16/r32 一様は step 2050–2150 で**安全停止** (case/61 と同じ誤検知の型: 壁温が下がってから上がる過渡)。**r32 非一様は step 200 の最初の更新で固体が 454.6 K** (範囲外) — 流体が一様 IC のまま連成を始めた起動の罠 (親 plan #99 ⑤ と同型) | 破棄 (対処 1 の根拠) |
| `run_0010_axisprobe_fp64` | 軸負例 1 step (FP64、`axisRFloor` 5e-4) | **`eval_vax2b.py --identity-only` PASS**: 床に当たった面を含め q_eff·A^r = Q_f 相対 0、raw 1.5e-16、床を外すと 0.375 (検出力あり) | ref |
| `run_0011_disk_r8u_w20k` 〜 `run_0013_disk_r32u_w20k` | **対処 1**: `warmup` 200 → **20000**。他は同一 | **一様系列**: 全界面節点 max 誤差 0.504 / 0.160 / 0.0477 % (比 1/3.15, 1/3.36)、N_r=32 PASS、恒等式 ≤1.7e-16、G-cons・連成保存・準定常・G-if・流体残差 PASS。**非一様が未達なので V-ax2b は未完了** | active |
| `run_0014_disk_r32g_w20k` | 同上 (非一様 r 方向 等比 1.1) | **step 20700 で安全停止 (固体 276 K)** — 2 回目の失敗。**原因は流体側**: warmup 終了時 (連成前) に静止ガスで max|U| **0.231 m/s** (一様 r32 は 4.3e-4)、**P 728–742 Pa** (一様は 1027、加熱なので上がるはず)、共役壁の壁熱流束が r≈6 mm で隣接節点ごとに符号反転 (±3000–4000 W/m²、物理 ≈120)。**軸対称 node の流体が非一様格子で静止を保てない**疑い → codex (diagnose) に諮り中 | 保留 (流体側の異常) |
| `run_0015_hoop_A_hoop0` / `run_0015_hoop_B_hoop1` | **判別 A/B** (plan §5.1 #6): `run_0014` の入力で `mesh.hoopAreaFromClosure` 0 / 1 だけを変え 500 step (連成前) | step 250–500 の max: A 21.86 m/s・19.67 Pa・市松 4.53e5 W/m²、B 21.88・19.63・4.51e5 → **B/A ≈ 1 で閉包補正は効かない** (第 1 仮説を下げる) | ref |

ローカルの起動確認 (2026-09-27、float、各 2 step) は 4 本とも init・固体検査・界面の座標一致 (ずれ 0 m) を通り rc 0。
評価器の経路確認 (`warmup 0`・`interval 1` の 3 step) も通した。いずれも確認後に破棄した。
