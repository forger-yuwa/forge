# 58. 深いスロット片壁の CHT (V6′ — 伝導漸近との照合)

計画: [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md) §6 **V6′**

**他セッションの資産に依存しない自己完結の CHT 検証**。深いスロットの深部は $Pe\ll1$ で伝導支配になるので、
**幅方向 1 次元伝導 + 固体の直列抵抗**の解析解に漸近しなければならない。**一様壁温だと漸近解が
$T\equiv T_w$・$q\equiv0$ の自明解 (null 検査) にしかならない**ので、両壁に温度差をつけるのが要点。

## 幾何と条件

```
 y=H  ┌──────────────────────────────┐ top (slip)
 y=0  ┤──────┬──┬────────────────────┤ plate (断熱)
             │  │ ← スロット W=1 mm
             │  │   前壁 slot_front = 唯一の共役壁 (固体帯、背面 Robin h=1e8 の T_c=300 K)
             │  │   後壁 slot_back  = 固定等温 500 K (熱側・非結合)
 y=-D        └──┘ slot_bottom (断熱)      D=20 mm  (D/W = 20)
```

| 量 | 値 |
| --- | --- |
| 自由流 | $M$=2、$T_\infty$=300 K、$p_\infty$=10 kPa、$U_\infty$=694.38 m/s、$\rho_\infty$=0.116144 |
| 乱流 | **なし (層流)** — v1a の制約を素直に満たすため |
| 物性 | **定数** (`viscMethod: 0`, `visc` 1.8469e-5, `thermCond` 0.0445, `cp` 1004.5, `gamma` 1.4)。**$Pr$=0.417** (空気の 0.72 ではない。深部の伝導漸近には効かないが、自由流の境界層を空気と読まないこと) |
| 固体 | $t$=1 mm、$k_s$=**0.05** W/mK ($R_s/R_f$=0.890)、16 層、5457 節点 / 10240 三角形、界面 321 点が壁節点と 1 対 1 |
| メッシュ | 32641 節点・全四角。幅方向 41 (Δx 25 µm)、深さ方向 321、チャネル法線 81。**品質 `VERDICT: PASS`** (AR 最大 94.4、skew 0.000) |
| **漸近解** | $q_*$ = **4709.0 W/m²**、$T_{w1,*}$ = **394.180 K**、固体の温度上昇 **94.18 K** (0.5 % ゲート = **0.471 K**) |

**向きに注意**: 背面を加熱してガスより熱い固体にする構成は**ソルバが拒否する** — 安全停止が固体温度を
[min($T_c$)−20, 流体の最大全温+20] に制限している (`conjugateWall.cpp`:774-800)。**冷却剤が冷側・ガスが熱側**にする。

## 生成物

```bash
python3 case/58.conjugate_slot/gen_mesh.py                       # .geo -> .msh
# 乾式 1 step で slot_front の壁ダンプを得る (角節点の physID 所属を推測しないため)
python3 case/58.conjugate_slot/gen_solid_strip.py --wall <run>/res_slot_front_5_1.h5 --nl 16 \
        --out case/58.conjugate_slot/mesh/solid_front_tc300_nl16
python3 solver_density_cuda/tools/solid_mesh_to_h5.py \
        --solid .../solid_front_tc300_nl16.json --out .../solid_front_tc300_nl16.h5
```

## 計算 run 一覧

| `run_*` | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| `run_0001_dry` | **乾式 1 step** — `slot_front` の壁ダンプから**界面節点 321 点を確定**させるためだけ (`conjugate` なし)。`plate` が等温だと共有角ガードで弾かれるので断熱にした | `res_slot_front_5_1.h5` (321 節点、$x$=0 一定・$y$ 0 … −20 mm)。メッシュ品質 `VERDICT: PASS` | ref (固体生成の入力) |
| `run_0002_v6p_spinup` / `run_0002_v6p_base` | **破棄**。IC が `uniform_p101325_u10` (P=101325 / u=10) で入口と桁違いだったため $p$ が **2.18e7 Pa** に達し `rms_roUy` RISING | 発散の記録。`initial: "slot_m2"` 追加の根拠 | 破棄 |
| `run_0003_v6p_spinup` | `initial: "slot_m2"` で soft 2000 → mid 4000。**壁温は `slot_front` 700 K** (本段と不一致だったのが後の問題の原因) | `rms_ro` −3.9 桁・全列 `falling`、$p$ 9915–14360 Pa | 中継 (破棄可) |
| `run_0004_v6p_base` | 破棄。`run_0003` (壁 700 K) から本段へ入り、**step 0 で安全停止** (固体が 712.7 K) | 「restart 元と本段で壁温が違うと最初の結合更新が過大」の記録 | 破棄 |
| `run_0005_v6p_base` | 破棄。向きを冷却剤 300 K / 後壁 500 K に入れ替えたが、restart 元がまだ壁 700 K で同じ停止 | 同上 | 破棄 |
| `run_0006_v6p_spinup` | **本段と壁温を揃えた spinup** (`slot_front` 394 K / `slot_back` 500 K)、soft 2000 → mid 6000 | `rms_ro` **−4.3 桁**。**以降の restart 元** | active (**IC の正本**) |
| `run_0007_v6p_base` | `Df_scale` **1.0** で本段。**step 2200 で安全停止** — `dTw_max` の比 1.195/1.404/1.312 = **増幅** | 振れているのは **lip 角 1 ノード**。`conjugate_history.csv` | active (発散の記録・§5.1 #98) |
| `run_0008_v6p_df5` | **判別 A/B — 変えたのは `Df_scale` 1.0 → 5.0 だけ** (3000 step = 20 更新、壁ダンプ 50 step 刻み) | **完走**。`dTw_max` 28.4 → 0.87 K 単調減衰、`res_abs` 3.5e5 → 7.8e3。**lip 角の $\Delta Q_f/\Delta T_w$ = 0.0694 W/mK** (完全緩和予測 0.062 に一致) → $\lambda$=+0.36。**事前登録の A 判定「比 ≤0.7」は未達** (0.181→0.755→0.97) | active (**機構の根拠**) |
| `run_0009_v6p_prod` | **float32 の報告 run** (`Df_scale` 5、`interval` 50、100k step)。**登録は「合否は FP64 で取る」なので合否そのものが無い** | **判定不能**。旧登録の帯で (a) 0.6353 / (b) 0.1328 / (c) 0.7429 / (d) 0.6101 % (許容 0.5、**帯平均は 0.0126–0.0734 %**)。**超過は帯の最上端 4–7 点だけ**。深部は $y$=−9.445 mm で $T_w$ **394.172 K** (解析 394.180)・$q$ **4709** (解析 4709.0)。**連成の伝達差は帯内 max 0.0097 %**・固体–流体の $T$ 差 **1.5e-5 K**。G-if は登録どおり `NOT CONVERGED` (① max 33.9 > 23.5、② 2.41e-3 > 1e-3) で lip 行のリップルに支配される | active (**報告項目**・§5.1 #97) |

| `run_0010_v6p_fp64` | **全域 FP64 ビルド** (`build-fp64`) で `run_0009` の最終場から 100k step。config は `run_0009` と同一 (`Df_scale` 5 / `interval` 50)。**登録上これが合否を取る run** | `eval_v6p.py` は **6 条件 PASS** ((a) 0.1528 / (b) **0.0088** / (c) 0.2047 / (d) 0.1959 % / (e) 0.0041 % / (f) 0.0000 K、帯 15.6 W)。帯平均は解析解と $T_{w1}$ **+0.009 K**・$q_{w1}$ **−0.025 W/m²**。深部 $Pe$ **6.5e-5 〜 2.0e-4** (float32 の 1.4–2.1e-3 から 1 桁下降 = §5.1 #87 の予測を裏付け)。**ただし総合は「判定不能」** (§5.1 #100): **(e) が保存性を測っていなかった** (`q_iface` は $Q_f$ のコピー。正しい物理残差は 0.0025 %) / **帯内 G-if の 80 更新連続の証拠が無い**。`check_convergence` `NOT CONVERGED`、全域 G-if `NOT CONVERGED` (① 28.4 > 23.5) | active (**本命・判定不能**) |
| `run_0011_v6p_df20_i50` / `run_0011_v6p_df5_i200` / `run_0011_v6p_df20_i200` | **$K$/$D_f$ 感度 2×2** (`Df_scale` 5/20 × `interval` 50/200)。すべて FP64・`run_0009` の最終場から 100k step・同一 IC | `eval_v6p.py` は 3 本とも 6 条件 PASS。**共通帯 211 行での `5/50` との最大差**: 温度 `20/50` 0.0231 / `5/200` 0.0211 / **`20/200` 0.2186 %** (許容 0.5)。**ただし正しい物理残差では `20/200` が 0.4607 %** で許容 0.1 % を超過。**局所量は `20/50` と `5/200` が各 13 点 `DRIFTING`/`TRANSIENT-UNSETTLED`**、`20/200` は温度 7 点・熱流束 6 点 → **感度も合格としない**。所要 224–703 s | active (感度・判定保留) |
| `run_0012_v6p_df5_i50_nlog` / `run_0012_v6p_df20_i50_nlog` / `run_0012_v6p_df5_i200_nlog` / `run_0012_v6p_df20_i200_nlog` | **`run_0010`/`run_0011` の 4 条件の再実行 (AWS g5・FP64)**。入力 (`mesh.h5` = `run_0009` の最終場・`solid.h5`・config) はバイト同一で、**差は `conjugate.node_log: 1` だけ** (界面節点ごとの $r_i,Q_{f,i},\Delta T_i,T_{w,i}$ を毎更新 `conjugate_iface_log_5.csv` に書く。数値は変えない)。**帯内 G-if (§5.1 #101 ②) と局所量の準定常 (③) を取るための run** | 全 run NaN なし・流体残差 `NOT CONVERGED (stalled/plateau)`。`df5_i50` は旧 `run_0010` を温度 0.0001 % で再現。**`df5_i50` (合否 run)**: `eval_v6p.py` 6 条件 PASS ((e) **0.0011 %**)、**帯内 G-if PASS** (① 0.240 W/m² / ② 2.26e-5 / ③ 2.66e-5 K、`CHT_INTERFACE_BAND_VERDICT.txt`)、**準定常 帯平均+局所 211 節点 ALL STEADY** (`QUASISTEADY_band.txt`)。全域 G-if は NOT CONVERGED (① 32.8、報告項目)。**比較 run**: `df20_i50`・`df5_i200` は 6 条件・帯内 G-if PASS だが帯上端 (lip 側) 15 節点の温度が DRIFTING/UNSETTLED (+0.05…+0.12 K、基準へ漸近中)。**`df20_i200` は (e) 0.5050 % FAIL・帯内 G-if NOT CONVERGED** (① 23.9)。感度 (`sens_v6p.py`、共通帯 211 節点): 温度 0.0231 / 0.0212 / 0.2186 %。成果物: `v6p_band.json`・`v6p_band_series.csv`・`conjugate_iface_log_5.csv`・`residual_history.png`。**解釈は codex (diagnose) 諮問中** | active |
| `run_0013_v6p_df20_i50_ext394k` / `run_0013_v6p_df5_i200_ext394k` / `run_0013_v6p_df20_i200_ext394k` | **延長判別** (plan §5.1 #102、事前登録)。`run_0012_*_nlog` の最終状態 (流体 `restart_field.py --keep-src-dtype`・壁温 `wall_profile_5.csv`+`wallProfile: 1`・固体 `conjugate_state_5.h5`) から **+294k step (累積 394k)**、他は同一・`outStepInterval` 10000。AWS g5・FP64・`FORGE_CUDA_BLOCKSIZE=128` | 共通帯 211 節点固定。**`20/50`・`5/200` は A** (局所 ALL STEADY・帯内 G-if PASS・(e) 0.0025 / 0.0015 %・基準との温度差 0.0000 %)。**`20/200` は判別未完** (G-if PASS・(e) 0.0278 %・温度差 0.0149 % だが帯上端 5 節点が TRANSIENT-UNSETTLED、基準へ単調に漸近中 −0.014 K)。NaN なし・流体残差 NOT CONVERGED | active |
| `run_0014_v6p_df20_i200_ext688k` | **`20/200` の再延長** (plan §5.1 #103、事前登録)。`run_0013_v6p_df20_i200_ext394k` の最終状態から計算長だけ +294k (累積 688k) | **登録条件 A**: 固定帯 211 節点適格、(a)〜(f) PASS ((e) 0.0018 %)、局所 211 節点 ALL STEADY、帯内 G-if PASS (① 0.196)、基準との差 温度 0.0004 % / $q$ 0.0017 %。流体残差・全域 G-if は NOT CONVERGED (報告項目) | active |

**次** (2026-09-26 の codex diagnose の結論。正本は
[plan §5.1 #101](../../plans/active/boundary-conjugate-heat-transfer.md)):
**① `eval_v6p.py` の (e) を直す** — 固体の物理作用素から $Q_{\rm sol}=(K_su-b_s)_{\rm iface}$ を組む
(現行は `q_iface` = 流体荷重のコピーを比べているだけ)。**② 帯内 G-if を毎更新で測る診断出力**を足す
(節点ごとの $r_i,Q_{f,i},A_i,\Delta T_i$ を節点 ID + 更新番号つき。固体ダンプに足すだけでは不足)。
**③ `check_quasisteady --series-csv`** で帯平均と**局所量**を判定。**④ 感度を共通帯・節点ごとに**取り直す。
