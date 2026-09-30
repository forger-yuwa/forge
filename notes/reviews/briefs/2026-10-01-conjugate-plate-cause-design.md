# 諮問: case C (共役平板) 前提ゲート FAIL の原因切り分け A/B の設計

AGENTS.md 条件 1 (plan §4・§6 の新規)・2・4。plan `plans/active/boundary-cht-conjugate-flat-plate.md` (§5.1 #1)。発注元 `plans/accepted/boundary-cht-conjugate-benchmarks.md` §4.3・§4.5・§4.7・§6。

## 観測事実 (run は case/65.conjugate_flat_plate/、AWS FP64 `86115cb1`、600000 step 単一区間、node・SLAU・陰解法 blockDPLUR cfl_pseudo 2・implicitRelax 0.7・limiter 2)
- 構成: 一様流 M0.1 Re_L 1e4、板 0≤x≤L (L 10 mm) の上面が共役壁 (physID 5, wall_isothermal conjugate:1)、前縁上流 L/2 (physID 4)・後縁下流 L/2 (physID 6)・上端 6δ (physID 3) は slip。入口 inlet_uniformVelocity、出口 outlet_statPress。固体下面 Robin h 1e8・T_h 310 K、端面断熱。conjugate: fem2d, q_eff, warmup 5000, interval 50, Df_scale 5, flux_avg 1。
- 流体収束 (`CONVERGENCE_CHECK.txt`、全 6 本): NOT CONVERGED (stalled)。n64 C1: rms_ro 1.8 桁・rms_roUy 2.0 桁・rms_roe 1.8 桁で横ばい、rms_roUx は 4.6 桁。n16/n32 も同型 (1.3〜1.8 桁)。
- G-if (`CHT_INTERFACE_VERDICT.txt`、末尾 80 更新):
  - n16/n32: ① res_abs 6.4 / 5.8 W/m² (> 1.0)、②③ OK、④ res_solid 3.7e-9 / 1.9e-9 (> 1e-9)。
  - n64 (C1・C2): ① 678 W/m²、② 0.074、③ dTw 4.1e-3 K、④ 1.9e-9 すべて NG。
- **④ res_solid は 2^-31 (4.657e-10) の整数倍の値しか取らない** (全 6 本、後半の distinct 値: n64 {4.66e-10, 6.99e-10, 9.31e-10, 1.40e-9, 1.86e-9}、n16 {9.31e-10 … 7.45e-9})。Robin h 1e8 × T_h 310 × 節点幅 ~6e-5 m ≈ 2e6 W/m の項の和の丸め床と見える。A (h_o 570) では ④ は PASS。
- n64 C1 の res_abs の時系列 (`conjugate_history.csv`): 5000 step 8314 → 以降 14〜742 W/m² で減衰せず振れ続ける (step 54550: 198、104100: 30、252750: 655、500500: 658、599600: 316)。末尾 400 更新で min 13.7 / max 742 / mean 183、dTw_max 1e-4〜4.4e-3 K。q_total は −14.94〜−14.96 W/m で安定。
- 界面残差の位置 (前回の評価): n64 は前縁直後 x/L 0.004〜0.015 に集中、n32 は後縁 x/L 1.0。評価窓 x/L 0.2〜0.9 内は 0.14〜0.15 W/m²。
- 場の時間変動が最大なのは後縁のすぐ下流 (x 10.1〜10.7 mm) の slip 境界 (最後の 2 枚で P 3 Pa、T 3e-3 K)。前縁上流の slip 節点に法線速度 0.72 m/s (U∞ 34.7 m/s)。
- 主判定 (独立参照解との比較) は 6 本 PASS。準定常 (評価窓内全節点) PASS。
- 既知: node slip + 接線方向の密度勾配で偽の流れが定在する未修正欠陥 (`notes/investigations/node-slip-tangential-density-spurious-flow.md`)。
- 中間スナップショットは AWS 上 (AWS は現在停止中。起動にはユーザ確認が要る)。手元には最終場なし (設定・ログ・評価出力のみ)。

## 仮説 (未検証)
H1: ④ は丸め床で、登録の tol_solid 1e-9 が h 1e8 では達成不能 (ゲート側の問題)。
H2: 流体の停滞は連成と無関係 (slip 境界の既知欠陥 / 前縁・後縁の特異点 / 出口)。
H3: n64 前縁の界面残差の振れは連成反復 (Df_scale 5、前縁の h の特異性) の不安定、または H2 の流体の揺れに駆動されたもの。

## 私の A/B 案 (n64 C1 の最終場から restart_field で同一メッシュ再開)
- B1: 連成を止め、壁を最終の界面温度分布で固定した等温壁 → 流体収束が回復するか (H2 vs 連成由来)。
- B2: Df_scale 20 (連成を強く減衰) → 前縁の振れが消えるか (H3 の連成側)。
- B3: slip_up/slip_down を別境界 (例: 板の上流・下流も no-slip 断熱壁、あるいは領域を短く) → 流体停滞・後縁下流の揺れが消えるか。
- H1: 丸め床を解析的に見積もり、tol_solid を h・T_h・節点幅の積で相対化する案 (登録の変更になる)。

## 質問
1. 判別として上の B1〜B3 は妥当か。優先順と合否条件 (何が出たらどの仮説を棄却するか) を事前登録したい。
2. H1 の ④ のゲート閾値の扱い (事後変更になる) をどう書くべきか。
3. 他に先に確かめるべき候補はあるか。
