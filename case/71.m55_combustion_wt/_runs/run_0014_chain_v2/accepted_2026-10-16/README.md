# ノズルの成果物一式 — `run_0014_chain_v2`

- **寸法の出荷ゲート**: **PASS** (`DIMENSIONS.json`: STEP とメッシュから測り直した寸法と要求カードの公差)
- **壁面の熱量を熱設計に使えるか (`thermal_use`)**: **PASS**
- **設計壁と物理壁の比較** (`wall_compare/`): Δr = r_物理 − r_設計 は設計スロート (x = 0) で -0.24341 mm・出口で 41.4569 mm。prepare の δ_r の記録との照合: **PASS** (`wall_compare/wall_design_vs_physical.json` の `delta_record_check`)
- **Fluent での読み込み**: 未確認 (手元に Fluent が無い。`fluent/nozzle_2d_axi.msh` を読んで `Mesh > Check` の結果を教えてほしい)
- 総合: **PASS**

## 単位・座標

- 長さ: CAD (STEP・形状表) は mm、メッシュ (msh・Fluent・CGNS) と壁面量表は m。
- x = 流れ方向 (回転軸)、y (= r) = 半径。原点 = 設計スロート (非粘性の設計壁の x = 0)。物理スロート (境界層込みの最小径) の位置は `DIMENSIONS.json`。
- 2D 軸対称の上半分 (r ≥ 0)。Fluent は 2D Axisymmetric、x 軸が対称軸。

## フォルダ

| フォルダ | 中身 |
| --- | --- |
| `cad/` | `wall_physical.step` (2D の壁曲線、5 次 B-spline 1 本)・`wall_surface.step` (3D の内面)・`wall_surface.stl`・検査の JSON |
| `mesh/` | `nozzle.msh` (解析に使った msh4.1)・`nozzle_v22.msh` (msh2.2、番号は gmsh が振り直す)・`mesh_structure.json` (構造ファイル: 半径方向の格子線の格子は v1、壁法線の近壁層の格子は v2 = 節点の規則 `wall_normal_blend_v1` と重み β・壁の傾き)・`rebuild_from_step.py` (STEP から同じ構造で作り直す: `python rebuild_from_step.py STEP mesh_structure.json OUT [--compare nozzle.msh]`、gmsh が要る)・その実行例 `rebuilt_from_step*` |
| | 作り直しの出力: `OUT.msh` = forge の変換器に入れる形 (元と同じ番号・書式、座標 17 桁)、`OUT_cad.msh` = 壁の節点が STEP の曲線に分類された gmsh の msh4.1 (**forge の変換器には入れない**: ファイル内の節点の順が番号順でない)、`OUT_v22.msh`・`OUT.cgns` (番号・物理タグは gmsh が振り直す。物理名で対応付ける) |
| `fluent/` | `nozzle_2d_axi.msh` (Fluent 形式、境界の名前と型つき)・`nozzle_2d_axi.cgns`・`boundary_conditions.md` (動作圧 0 Pa で絶対圧を入れる) |
| `wall_compare/` | `wall_design_vs_physical_mm.csv` (同じ x [mm] で設計壁 = Euler で回す壁 (境界層なし) と物理壁 = 作る壁 (設計壁 + 境界層の厚み) の半径・差 Δr = r_物理 − r_設計・壁角。x の点は形状表と同じ。xlsx のシート `wall_compare` にも入れた)・図 (`*_full.png` 全長の壁の形・`*_dr.png` Δr の分布・`*_throat.png` スロート付近の拡大)・`wall_design_vs_physical.json` (列の意味・区間・δ_r の記録との照合) |
| `excel/` | `wall_shape_mm.csv` (形状表)・`wall_cfd_nodes.csv` (壁節点ごとの圧力・温度・せん断・熱量)・`conditions.csv`・`nozzle_tables.xlsx` |

## 注意

- forge は節点中心、Fluent はセル中心。同じメッシュでも壁の第 1 未知数までの距離は Fluent が約半分になる。
- 熱量の列は `.status.json` の状態 (`valid` / `not_applicable` / `missing` / `diagnostic` / `nonfinite`) と `thermal_use` を見てから使う。`diagnostic` (再構成勾配の `qwall`) は正本でない。
- 軸対称の壁への熱量は `iface_Qf_eff` (1 rad あたりの積分値) × 2π から作った (`iface_q_eff` は軸対称では r·q で熱流束ではない)。
- 壁の両端の節点 (入口・出口の角) の熱量は信頼できない (別の壁との境目で再構成勾配との比が 2〜4 倍になった例がある。原因は未確認)。`excel/wall_cfd_nodes.status.json` の `end_nodes` を見る。
- STEP を CAD で編集したら `rebuild_from_step.py` で作り直し、`check_mesh_quality.py` で格子品質を確かめる。
- 壁法線の近壁層の格子 (`mesh_structure.json` が v2) の作り直しは、STEP の曲線の傾きで近壁の節点を壁の法線上に置く。入口・出口の面は直線 (x 一定) のままなので、端の列の節点が直線から外れることがある (外れる量は `mesh/rebuilt_from_step_compare.json` の `end_columns`)。構造ファイルを作るときに、格子の生成の記録 (run の `prepare_info.json` の `mesh.params`) の規則のパラメータで節点を置き直して元と照合した (`MANIFEST.json` の `structure_rule_check`。記録が無ければ「未実施」)。
- STL (`cad/wall_surface.stl`) は表示・形の確認用の近似: 弦誤差の設定 線の偏差 0.05 mm・角度の偏差 0.05 rad (三角形 822660 枚)。寸法の正本は STEP。
- CGNS (`fluent/nozzle_2d_axi.cgns`) を Fluent で読むと**境界の型は付かない** (入口・出口・壁・軸を Fluent 側で付け直す)。gmsh で読み直すと物理タグの番号も振り直されるので、境界は**名前**で対応付ける。型つきの Fluent 形式は `fluent/nozzle_2d_axi.msh`。
- 設計壁と物理壁の比較 (`wall_compare/`) の Δr は区間で意味が違う:
  - x ≥ 0 (設計スロートから出口): Δr = δ_r = 境界層の補正 (物理壁 = 設計壁 S + δ_r)。冷却壁ではスロート付近で負になりうる
  - 縮流部 (x < 0): 物理壁は多項式 Q、設計壁は Hermite で別の曲線。Δr は境界層の厚みではなく 2 つの上流の形の違い (x = 0 で δ_r(0) に連続)
  - 直管: Δr = 物理壁の r_U − 設計壁の r_U (この壁では 0.0 mm)

## 出典

- run: `/home/sano/work/forge-cardmesh/case/71.m55_combustion_wt/run_0014_chain_v2`
- 作成: 2026-10-16 15:15:15、forge_design の git commit `4e8949bd`
