# forge_design capability registry (対応状況の正本)

設計キャンペーンの立ち上げ時 (インテーク対話) に「何が選べるか / 何が要実装か」を
判断するための対応表。**メニューの品揃えは runner が YAML から solverConfig に翻訳
できる範囲で決まる**ため、対応状況は問題タイプ (runner) ごとに記す。

- **更新規律 (必須)**: アダプタ・メニュー・runner キーを追加/検証/廃止したら本表を
  同期する (case README の run 一覧表と同じ精神)。状態を上げるには検証ケースを 1 つ
  伴うこと。
- 凡例: ✅ 検証済み / 🔧 実装あり (検証限定・条件付き) / 🧩 リポジトリに資産あり・
  設計パイプラインへの配管が未 / 📋 未実装 (plans/ に切り出してから)

## 1. 問題タイプ ([`probdef.py`](forge_design/probdef.py) `KNOWN_TYPES`)

| type | runner | 状態 | 備考 |
| --- | --- | --- | --- |
| `thruster_bell` | `evaluate/runner.py` | ✅ | MOO ドライバ ([`opt/driver.py`](forge_design/opt/driver.py)) の対象。case/40 run_0074 で 69/70 PASS |
| `thruster_bell` (TP 排気) | `evaluate/runner.py` (`gas.model: semiperfect` + `gas.exhaust`) | ✅ | ロケット排気の凍結 TP (CEA2 の燃焼室平衡組成を lump `EXH` に、`gas/rocket_exhaust.py`)。1D 理想値は CEA 凍結と 2.1e-5 以内。推力の CV 会計 (`thrust_cv`)・ゲート (`tp_gates`、plan [tooling-nozzle-bell-two-mode](../plans/active/tooling-nozzle-bell-two-mode.md) §4.1)・同一メッシュの継続。TOP 基準ケース (case/68 run_0008) と (D) の 1 点 (run_0009) が全ゲート PASS |
| `thruster_bell` (D) 直接モード | `geometry.bell_mode: direct_bspline` (`geometry/wall_bell_c2.py`) | ✅ | C2 の 5 次 B-spline、dv $(L,\theta_e,\theta_{max})$、$x_m$ は内部決定。前進 MOC ゲート (`geometry/moc_forward.py`) は記録だけ。MOO は case/68 run_0010 |
| `thruster_bell` (M) 特性曲線モード | `geometry/rao_axisym.py`・`geometry.bell_mode: moc_rao` (`geometry/wall.py` の `RaoBellWall`) | 🔧 | TP の軸対称 Rao (γ-free の制御面条件、円弧 kernel・B で C1、受理は det_D > 0)。S1 を通過 (plan §5.1 #5)。runner は生成器の壁ファイル (`geometry.rao_wall_file`) を読み、(D) と同じ収縮部につなぐ (スロートで C2)。L・ε・R_t・排気ファイルの不一致で止める。生成器の推力を N・η に換算して `WALL_DESIGN.json` に記録 (V2 の参照値)。Euler (`evaluate.model: euler`、slip 壁・壁際を細分しない格子): prepare は変換器の手前まで試験済み、推力の CV 会計は `thrust_cv(wall_model="slip")`、ゲートは y1+ を対象外 (理由を GATES.json に)。forge での Euler run は未実施 (plan §5.1 #6 (W2)) |
| `wind_tunnel_axisym` | `evaluate/runner_wt.py` | 🔧 | 旧生産系 (B8)。axismach に置換済み、回帰対照として保持 |
| `wind_tunnel_axisym_walldriven` | `evaluate/runner_walldriven.py` | 🔧 | 壁駆動系。axismach に置換済み |
| `wind_tunnel_axisym_axismach` | `evaluate/runner_axismach.py` | ✅ | 現行生産 (case/42 isobutane, case/44 vitiated air)。ガス・凝縮メニューが最も厚い |
| `sern_2d` | `evaluate/runner_sern.py` | ✅ | ⑤ SERN: 燃焼器出口 starting line + 2D 平面最大推力 key point 逆設計 (`geometry/moc_sern.py`) + メッシュ/評価 (`meshing/mesh_sern.py`, `evaluate/runner_sern.py`, `metrics/sern_forces.py`)。**cell Euler で MOC・NASA TM X-71972 傾向と照合済み** (case/46 run_0002–0007)。RANS は node/cell とも SST で完走 (node の発散は stage 間 interp 移植が真因、index コピーで解決)、δ* 一発補正 (`metrics/sern_deltastar.py`)、多作動点 MOO (`opt/driver_sern.py`; case/46 run_0010 Euler 17 PASS, run_0017 node SST 10 PASS + C_M 制約) まで検証済み。plan: [`plans/active/tooling-nozzle-sern-chain.md`](../plans/active/tooling-nozzle-sern-chain.md) (2026-09-04 起票) |

## 2. ガスモデル・物理メニュー (`evaluate.*` / `gas.*` キー)

| メニュー | キー | 対応 runner | 状態 | 備考・罠 |
| --- | --- | --- | --- | --- |
| CPG (定比熱) | `gas.model: cpg` (既定) | 全 | ✅ | |
| semi-perfect TP (NASA-9) | `gas.model: semiperfect`, `evaluate.cfd_gas` | axismach | ✅ | CEA 凍結流と 0.04% 一致 (case/44)。**`thermo_href_temp: 298.15` 必須** (絶対基準 h だと χ_eos 桁違い→発散) |
| TP 擬似種の分割 | `evaluate.tp_species: {mode: full\|lumped, lumps, keep}` (旧 `pseudo` / `split_h2o` / SERN `[EXH, AIR]` は別名) | axismach / sern | ✅ | `full` = 各種を CEA NASA-9 の独立種、`lumped` = ユーザ指定の lump に畳む (plan thermophysics-cea-mole-fraction-species) |
| 組成のモル分率入力 | `gas.composition_basis: mole\|mass`, `gas.condensing_species`, `gas.species_db` (CEA 直読み DB) | axismach / sern | ✅ | 換算は `gas/composition.py` に一元化、`species_db.yaml` に由来コメント、`species_meta.yaml` に機械可読メタ |
| 凝縮 (非平衡/平衡) | `evaluate.condensation: {…}` → forge `condensation` ブロック素通し | axismach | ✅ | `condEquilibrium: 2` (EOS 拘束形) を既定推奨。蒸発は既定 ON。h0 保存を確認する |
| 乱流 SST | `evaluate.turbulence: sst\|none` | runner.py (bell) | ✅ | axismach は Euler/NS 系が主 |
| 化学反応 — 凍結流 | (= TP semi-perfect) | axismach | ✅ | 組成固定 |
| 化学反応 — 平衡ブラケット | CEA2 で凍結↔平衡を挟む | — | 🧩 | CEA2 native ビルド済み (`.venv-cea/nasa_cea/FCEA2`)。設計側ガスモデルへの組込みが未 |
| 化学反応 — 有限速度 | — | — | 📋 | **大型** (剛性ソース + point-implicit + 検証ケース)。多成分 implicit の生成エンタルピー増幅 (既知) が本丸化する。選ばれたら専用 plan |
| 高度・背圧 (P_amb/T_amb) | — | — | 📋 | プリュームチャンバー配管 (§4) とセットで実装 |
| 軸対称オプション | `evaluate.axisym_method`, `axis_r_floor` | axismach | ✅ | node は軸ノード DOF 整合セット常時 ON |

## 3. メッシュ経路

| 経路 | 実体 | 状態 | 備考・罠 |
| --- | --- | --- | --- |
| プロシージャル構造 quad | [`meshing/mesh2d.py`](forge_design/meshing/mesh2d.py) (msh4.1 直書き) | ✅ | 設計ループの既定。dv→壁形状→メッシュを毎評価生成。品質ゲート込み |
| Gmsh `.geo` | `procedures/calculation-workflow.md` | ✅ | 固定形状評価向き。プリューム付きテンプレ = `case/23.axi_nozzle/mesh_axisym_m4/axi_nozzle_2d_publicrao_plume.geo` |
| Fluent CFF HDF5 | `tools/fluent_h5_to_forge.py` | ✅ | tet+prism/非直交 hex とも検証済み (fan 完全収束, StaticMixer 完走)。**レシピ厳守**: 1 次で投入 / IC を入口流速に整合 / 陰解法 / `space.pRef`=動作静圧 / `outlet_statPress`+逆流 Pt,Tt (calculation-workflow「Fluent メッシュの取り込み」が正本)。多面体は可視化不可・単位 `--scale`・`mapping.yaml` 整備。設計ループ (dv→再メッシュ) には未接続 = 固定形状評価用 |
| Netgen/FreeCAD | case/37 ブリッジ (msh4.1)、`.venv-mesh` | 🔧 | 3D 複雑形状 (ピントル)。boolean の罠は case/37 メモ参照 |

共通の罠: node 計算のメッシュは **node config で変換必須** / RANS メッシュは no-slip
壁 bcond で変換しないと wall_dist=0 / 変換後は `check_mesh_quality.py` VERDICT 必須。

## 4. 評価領域・目的関数

| 項目 | 状態 | 備考 |
| --- | --- | --- |
| ノズル内部評価 (η_CF, Cd, 軸 M, ε_M, 推力) | ✅ | `metrics/`, `evaluate/health.py` |
| プリューム/チャンバー (高度評価) | 🧩 | case/23 に実証資産一式: plume `.geo` 6 ブロックトポロジ / `setInitial_plume_outer.py` (outer 大気静止 IC) / 段階起動レシピ (`run_ext_sweep.sh`) / `chamber_metric*.py`。要実装: `mesh2d` へのブロック移植 + `plume:` YAML 節 + P_amb 配管 + 推力の (P_e−P_amb)A_e 化。低高度は剥離・非定常で OSCILLATING 前提・評価コスト増 |
| 熱構造 (片方向: q_w → FEM 熱伝導/熱応力) | 📋 | forge 側の壁熱流束は `sstEnergyWallFunction` (Kader q_w) あり。OSS 候補: CalculiX (軸対称第一候補) / Code_Aster。ファイル受け渡しで §5 の VERDICT 契約に載せる |
| 共役熱伝達 (双方向 CHT) | 📋 | 片方向とは別次元の難度 (界面反復・緩和・安定性)。選ばれたら専用 plan |

## 5. ゲート・台帳 (どの段も VERDICT をファイルで残すのが契約)

| 項目 | 実体 | 状態 |
| --- | --- | --- |
| 残差収束 / 発散検出 | `tools/check_convergence.py` (`CONVERGENCE_VERDICT.txt`) | ✅ |
| 派生量の準定常 | `tools/check_quasisteady.py` (`QUASISTEADY_VERDICT.txt`) | ✅ |
| 目的量の定常 (軸 M トレンド+SE) | `evaluate/health.py::axis_M_series_steady` | ✅ |
| 物理ゲート (η, Cd; Cd>1 は SUSPECT) | `opt/driver.py` | ✅ |
| 失敗分類 `fail_class` (RETRYABLE/INFEASIBLE/SUSPECT) | `opt/driver.py` ledger (commit a9e923cb) | ✅ |
| RETRYABLE の自動再投入 | — | 📋 |
| feasibility classifier (FAIL 点のサロゲート活用) | — | 📋 |

## 6. 最適化・帰還

| 項目 | 実体 | 状態 |
| --- | --- | --- |
| LHS DOE / Kriging / NSGA-II + EHVI / ledger / warm seed | `opt/` | ✅ (case/40) |
| Euler 帰還ループ (特性線 C⁻ マップで壁修正) | `feedback/euler_loop.py` | ✅ (case/41) |
| NS δ* 補正 | `feedback/deltastar.py` | ✅ (相関で十分・無帰還 RANS 0.533%) |

## 7. 外部参照ツール

| ツール | 実体 | 状態 | 用途 |
| --- | --- | --- | --- |
| NASA CEA2 | `.venv-cea/nasa_cea/FCEA2` (native) | ✅ | 平衡/凍結の理論参照。pi/p は 1 problem 16 個制限 |
| SU2 | `procedures/su2-cross-check.md` | ✅ | 同一メッシュ・同一 BC の切り分け比較 |
