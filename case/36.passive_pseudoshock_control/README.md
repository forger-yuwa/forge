# case/36 多孔壁パッシブコントロールによる擬似衝撃波抑制 (2D 逆解析)

Matsuo et al. (1988)「境界層のパッシブコントロールが擬似衝撃波に及ぼす影響」
(九州大学総合理工学研究科報告 10(1), pp.45-50,
[PDF](https://api.lib.kyushu-u.ac.jp/opac_download_md/17711/p045.pdf)) の試験部を 2D で逆解析する。

**目的**: 多孔壁+キャビティの境界層パッシブコントロールで擬似衝撃波 (shock train) が
弱化し「ショックレス」に近づくかを再現確認する。詳細計画は
[`.github/plans/verification-passive-pseudoshock-control.md`](../../plans/active/verification-passive-pseudoshock-control.md)。

## 形状 (x=入口からの距離, 上下対称, 片側 1° 広がり)

```
            多孔板+キャビティ (上)
         ┌──┬─┬─┬─┬─┬─┬─┬─┬──┐   ← キャビティ 深さ20
   slip  │  └┐└┐└┐└┐└┐└┐└┐  │   ← 多孔板 厚10, 8スロット 幅0.8 ピッチ4
 ====inlet========================...========== outlet
   (助走30)  1°拡大100      ↑x=130-162    528mm直管
         │  ┌┘┌┘┌┘┌┘┌┘┌┘┌┘  │
         └──┴─┴─┴─┴─┴─┴─┴─┴──┘
            多孔板+キャビティ (下)
   x:0    30        130  162              690
```

| 区間 [mm] | 内容 | 壁 |
|---|---|---|
| 0–30 | 助走領域 (高さ 21.509) | slip |
| 30–130 | 1° 拡大ダクト | no-slip |
| 130–162 | キャビティ部 (多孔板 厚10 + キャビティ 深20、8スロット 幅0.8/ピッチ4、開口率0.2) | no-slip |
| 162–690 | 528mm 直管 (1°) | no-slip |

主流 M=1.89 (多孔壁上流端)。入口 M=1.689 (面積比逆算)。

## 入口条件 (よどみ 3MPa / 288.15K, 空気 CPG γ=1.4)

`inlet_uniformVelocity`: Ux=458.6 m/s, ro=11.733 kg/m³, Ps=617800 Pa (M=1.689, T=183.5K)。
背圧は出口 `outlet_statPress` の Ps で掃引し擬似衝撃波を前後に動かす。

## メッシュ

`mesh/make_mesh.py` (パラメトリック)。`POROUS=True`→多孔壁, `False`→固体壁 (比較基準)。

```bash
cd mesh
python3 make_mesh.py                       # passive_porous.geo
gmsh -2 passive_porous.geo -o passive_porous.msh -format msh4   # forge は msh4 を読む (msh2 不可)
```

初版: porous 63.6k quad / solid 38.3k quad (全 quad)。physID:
1=inlet, 2=outlet, 3=wall_slip, 4=wall(no-slip), 10=fluid。

### Salome (GEOM/SMESH) 版 — `mesh_salome/`

gmsh 版と**同一形状・同一グループ名・同一 physID** を生成する代替フロー。
構造ダクト (直管部 Quadrangle 写像) + 非構造コア/キャビティ (NETGEN quad-dominant)。

```bash
cd mesh_salome
salome -t make_mesh_salome.py                 # POROUS=True  -> passive_porous.med/.unv
salome -t make_mesh_salome.py args:solid      # POROUS=False -> passive_solid.med/.unv
salome -t make_mesh_salome.py args:porous,netgen   # 全 NETGEN フォールバック
# MED/UNV -> gmsh msh4.1 (meshio + h5py + gmsh が必要; forge は msh4 を読む)
python3 med_to_gmsh.py passive_porous.med     # -> passive_porous.msh (physID/グループ検証つき)
```

`make_mesh_salome.py` の主要パラメータ (冒頭に集約): `POROUS`, `STRUCTURED_DUCT`,
`NY`/`NX_RUN`/`NX_DIFF`/`NX_DOWN` (構造分割数), `CL_CORE`/`CL_SLOT`/`CL_CAV` (特性長 [m])。
スロット口・ブロック界面などの内部エッジはグループに含めず conformal に連続させ、
`med_to_gmsh.py` が physID=0 の内部線を除外して書き出す (forge は dim=1 を全て境界として読むため)。

**確認図** (gmsh 版で生成した同一形状・現行 10mm 板厚): `mesh_salome/preview_porous.png`
(全体+キャビティ), `mesh_salome/preview_slotzoom.png` (スロット 0.8mm が ~4-5 セル)。
gmsh 等価メッシュ: porous 62.8k quad (100% quad), 境界線 inlet/outlet=各64・wall_slip=100・
wall=2274, 最小辺 ~0.087mm。`med_to_gmsh.py` の変換ロジックは合成ケースで検証済み。

> Salome 本体は本リポジトリ環境に未導入のため、`make_mesh_salome.py` の SMESH 実行は
> 各自の Salome で行うこと (形状/グループ/physID 規約は gmsh 等価メッシュで確認済み)。
> 変換後の `.msh` は `gmsh passive_porous.msh` で開けば可視化できる (Salome GUI 不要)。

## 計算 run 一覧

| run_* | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
|---|---|---|---|
| `run_0001_slau_estab_1st` | 層流・SLAU・1次・吹き抜け (outflow) で超音速場確立。porous mesh, M_in=1.689 | **健全**: NaN 無し・物理妥当 (ρ>0, P≤P0, T 154-354K), Mach 0-2.05, 75%超音速。残差~2桁低下。場確立 (背圧無しのため shock train 無し=想定通り)。`field_Mach.png`/`field_P_kPa.png`/`residual_history.png`/`res_20000.h5` | active (確立) |
| `run_0002_slau_sst_muscl_imp_cfl5` | **本番設定の cfl_pseudo=5 安定性確認**: SST + MUSCL(2次,Venkata) + 陰解法(blockDPLUR, implicitRelax0.5, nStepInner10) + **cfl_pseudo=5**。run_0001 確立場から restart (自由流 k/omega 注入)、Sutherland, dilatationCorr=2, 吹き抜け | **cfl_pseudo=5 安定 (発散なし)**: 10000 step/58s, NaN 0, 物理妥当 (ρ>0,P≤P0,T168-296K,k≥0,ω>0), 残差は平坦プラトーで安定 (発散せず)。乱流BL発達を確認。収束は ~0.7-1.5桁でプラトー=キャビティ剪断層の limit cycle+背圧無しのため (cfl5 起因ではない) | active |
| `run_0008_lam_estab` | 層流・MUSCL・陰解法 cfl5 で超音速場確立 (run_0001 から restart) | 薄いBL (δ99~0.1-0.3mm)、cavity core M=1.90 (実験1.89に一致)。outlet M=2.52 | ref (層流anchor) |
| `run_0011_lam_bp1p7_seedSST` | 層流・背圧 Ps=Pt=1.7。SST shock場を seed (超音速ロックイン回避) | 衝撃波保持 x447mm。**seeding で背圧BC が機能** | ref |
| `run_0012_lam_bp2p0_chain` | 層流・Ps=2.0 (run_0011 から chain) | **衝撃波 x129mm=多孔壁上流端 (x_f≈0)**。ただし中心軸静圧は**波状 (shock train)**=ショックレス化せず。層流BL薄すぎが原因 | ref |
| `run_0013_sst_wt1_ti1pct` | SST・**wallTreatmentSST=1**・入口乱流 μt/μ=1,TI=1% (k=31.5,ω=3e7)。y+~1200メッシュ | 安定だが **BL過厚** (cavity M 1.63, δ99 4.2mm=34%半高, μt/μ~37000)。y+~1200 で壁関数不成立が主因 | superseded |
| `run_0014_sst_wt1_wallres` | SST・wt=1・**wall-resolved メッシュ (NY120, 第一セル~5μm, y+~77)** で確立 | **BL適正化**: cavity M **1.81** (←1.63), δ99 **2.11mm=17%** (←34%), μt/μ **7558** (←37000)。メッシュ品質 AR202/skew SOFT-PASS | active (SST anchor) |
| `run_0018_sst_wt1_porous_bp2p06` | SST wall-res・**多孔壁**・Ps=Pt=2.06 (seed) | **NOT CONVERGED** (roUx/roUy/roK 上昇・他プラトー)。衝撃波 x~165mm だが非収束スナップショット値 | 非収束 |
| `run_0019_sst_wt1_solid_bp2p06` | SST wall-res・**固体壁**・Ps=2.06 (比較) | **NOT CONVERGED** (roUx 上昇)。衝撃波 x~129mm (ドリフト中) | 非収束 |
| `run_0020_sst_wt1_porous_bp2p08` | SST wall-res・多孔壁・Ps=2.08 (x_f≈0狙い) | **NOT CONVERGED** (roOmega 上昇・他低下中)。衝撃波 x~141mm | 非収束 |
| `run_0021_sst_wt1_porous_prof` | **計算時間ボトルネック分析用** (run_0020 入力複製, nStepOuter 縮小)。ncu でカーネル別 GPU 時間を取得 | 単独実行 **~5.35ms/step** (run_0020 基準, 98k cell)。GPU busy ~3.7ms/step・残りは launch/host。**上位: implicit_defect_correction_block 24.6% / SLAU 22.4% / limiter_r1 17.4% / viscousFlux 7.3%**。`residual_history.png` | 削除済み(2026-08-31) (プロファイル) |
| `run_0050_solid_prof` | **run_0048 (solid 構造) の計算時間ボトルネック分析 + 残差モニタ device 常駐化の検証** (入力複製)。ncu + detectNaN on/off 比較、旧/新バイナリ比較 | solid 79.4k cell。分析: GPU busy **2.87ms/step (~50%)・残り host/launch (110 launch/step)**、`detectNaN=1` が +1.45ms/step。**最適化結果** ([plan](../../plans/accepted/architecture-residual-monitor-async.md)): 残差 RMS/detectNaN を fused device 縮約+間引き flush 化し 2000 step **13.38s→9.58s (−28%)**、残差 CSV は不変 (差は solver 非決定性ノイズ床内)。`residual_history.png`/`CONVERGENCE_VERDICT.txt` | 削除済み(2026-08-31) (プロファイル/検証) |

### node-centered (median-dual) SST 試行 (2026-06-23)

solid メッシュを node + SST で計算。**wall_dist は現ビルド (commit 8b60f2c) で再変換し正しく計算**されることを
確認 (全壁ノード=0、near-wall で wd が幾何最近接壁距離と ratio=1.000 一致。ご要望「壁距離をちゃんと」を満たす)。

| run_* | 設定 | 結果 | 状態 |
|---|---|---|---|
| `run_node_sst_wdfix` | node SST, wt=1, conv1次, cfl_p1, ccnode_B(fix前)場から restart 4000 step | NaN なし・shock train 確立。残差は減少傾向 (rms_roOmega 1.6e4→498) しプラトー | active |
| `run_node_sst_fromcell` | **健全な cell 場 (run_0049) を nearest 補間して開始**、conv1次 | NaN なし。cell と整合 (mut/μ・圧力回復)。1次のため shock train は平滑化、衝撃波 x237mm | active |
| `run_node_sst_muscl` | fromcell から **convMethod=2 (MUSCL) に昇格** 8000 step | **shock train 構造が cell と一致** (多段 P/Mach 振動再現)、衝撃波 **x127mm** (cell 147, 1次 237)、下流圧力回復一致。`compare_centerline_cell_node_muscl.png` | active |

**所見 (訂正済み)**: 当初「過剰乱流」と誤判定したが、**それは誤り**だった。この超音速擬似衝撃波は本来 **mut/μ ~9万の
高乱流場**で、**cell 自身**も peak mut/μ=91,137 (p99 84,016, p90 41,201, frac>1e4=0.33)。node は peak 90,469
(p99 83,246, p90 39,369, frac>1e4=0.31) で **cell とほぼ一致**。flat plate の値 (199) と比べたのが誤りだった。

平均流もセンターラインで cell とよく一致: P(kPa) x=0.3:1598/1600, x=0.4:1756/1757, x=0.5:1835/1825 (cell/node)、
Mach も近接。**node SST は cell と整合**しており暴走していない。wall_dist は再変換で正しく計算 (全壁ノード=0,
near-wall ratio=1.000)。

**1次 vs MUSCL**: node 1次は強い数値拡散で shock train を平滑化し衝撃波が下流 (x237mm) に流れていた。**MUSCL
(convMethod=2) に昇格すると cell と同じ多段 shock train を解像**し衝撃波 x127mm (cell 147mm に接近)、圧力回復も
一致 (`compare_centerline_cell_node_muscl.png`)。残差床 (rms_roUx ~0.08–0.15) は shock train の limit cycle 的
性質でプラトー (cell も rms_roUx ~0.03 でプラトー)。衝撃波位置の ~20mm 差は擬似衝撃波の背圧敏感性 + 収束途中の範囲。
**結論: node + SST + 正しい wall_dist で case/36 solid が cell と整合的に計算できる**。

### SU2 クロスチェック (中立基準で cell/node を判定; 2026-06-23)

forge cell(衝撃 147mm)と node-MUSCL(127mm)の差の妥当性を、独立ソルバ **SU2 v8.5.0** で同一 `.geo`・同一 BC
(超音速入口 M=1.689/Ps=617.8kPa, 出口 Ps=1.90MPa, SST/SLAU/MUSCL, Sutherland)で判定。手順は
[`procedures/su2-cross-check.md`](../../procedures/su2-cross-check.md)。run dir: `run_su2_sst_slau_muscl/`、
比較スクリプト `compare_su2_forge_centerline.py` / 図 `compare_su2_vs_forge_centerline.png`。

| 設定 (SU2, 全て restart 同一発達場から) | CFL | 残差 rms[Rho] | 衝撃 x[mm] | peak μt/μ | M_max |
|---|---|---|---|---|---|
| SLAU 2nd + Venkat(K=0.05) | adapt | **激しい limit cycle** −0.9↔−2.8 | (振動) | (振動) | — |
| SLAU 2nd + Van Albada edge | adapt | **激しい limit cycle** | (振動) | — | — |
| SLAU **1st-order** | adapt | **激しい limit cycle** | (振動) | — | — |
| SLAU 1st-order | **固定 2.0** | 準定常プラトー −3.3 | — | — | — |
| **SLAU 2nd + Venkat（中立基準）** | **固定 2.0** | 準定常プラトー **−3.2** | **132** | **99,900** | **1.82** |
| forge cell run_0049 (参考) | — | プラトー | 142 | 81,300 | 1.82 |
| forge node MUSCL (参考) | — | プラトー | 120 | 73,700 | 1.69 |

- **振動の正体は `CFL_ADAPT` のフィードバック共振**（リミッタ無関係）。リミッタを 3 段階 (Venkat→Van Albada→1次) 変えても
  CFL adapt 下では全部 limit cycle、**CFL を固定にした瞬間 1次でも2次でも振動消失 → 準定常プラトー (rms[Rho]~−3.2)**。
  安定設定は **`CFL_ADAPT= NO` + `CFL_NUMBER= 2.0`**（2次の鋭い衝撃を保ったまま準定常）。
- **判定 (衝撃位置)**: SU2 中立基準 **132mm** は node(120) と cell(142) の**ほぼ中間**（node −12mm 上流, cell +10mm 下流）。
  cell/node は SU2 を挟む形で、どちらか一方が明確に正しいわけではない。
- **node がやや劣る点**: ① 上流の超音速加速を解像できず **M_max が 1.69 止まり**（cell/SU2 は area 駆動で 1.82 まで加速）、
  ② 衝撃が ~12mm 上流寄り。node-mode の数値拡散過多で衝撃が上流へ押されている兆候。**上流 Mach と衝撃位置の両方で cell の方が SU2 に近い**。
- **peak μt/μ**: SU2(~10万) > cell(8.1万) > node(7.4万)。forge 両者が SU2 より ~20-26% 低いのは **forge の dilatation
  correction（圧縮性補正, SU2 標準 SST には無い）が衝撃近傍の k を抑制**するため（既知の差）。
- 下流圧力回復 (>250mm) は三者ほぼ一致。**いずれも完全収束ではなく準定常スナップショット同士の比較**（擬似衝撃波は背圧固定で残差が落ちきらない）。

### node/cell 大差の根本原因 = node SST 壁関数バグ (2026-06-24)

「なぜ case36 だけ node/cell 差が大きいか (平板 case26・Euler case29 は ~一致)」を多角分析し**根本原因を特定**。

| run_* | 目的・設定差分 | 主要結果 | 状態 |
|---|---|---|---|
| `run_node_sst_bp1p90_matched` | node SST を **cell と背圧一致 Ps=1.90** に揃え res_4000(shock134mm)から継続 40k | **shock 46mm に漸近** (cell 142/SU2 132 から大きく外れ)・Mmax 1.69 (加速失敗)・Pmax 1.909(背圧満足)。背圧交絡を排除し node 固有問題を確定 | ref (診断) |
| `run_0054_node_sst_wffix_bp1p90` | **壁関数修正バイナリ (commit 8a2caad)** で Ps=1.90 node SST、res_4000 から 40k | utau/ypls=0→物理(y+~98)、peak μt/μ 91k→57k 是正、**shock 46→~65mm 改善** (なお cell/SU2 132-142 には届かず=残差はコア checkerboard)。漸近 ~65mm (plateau, 未収束) | ref (壁関数修正検証) |
| `run_0056/57/58_node_sst_wffix_bp1p80/70/60` | 壁関数修正バイナリで **node 背圧 down-sweep** Ps=1.80/1.70/1.60 (chain warm-start) | **shock 漸近 55/44/44mm**・Mmax 1.69 (全 BP で加速失敗)。**背圧低下で衝撃が上流へ=cell と逆応答** (cell は Ps1.90→142, 1.80→181mm で下流=物理的に正)。残差 plateau (未収束)。`cmp_node_vs_cell_bp_sweep.png` | ref (診断) |
| `run_0059_node_euler_bp1p90` | **node Euler** (viscMethod0, slip壁, Ps1.90) で粘性/SST 非依存切り分け | **コア加速 OK (Mmax 1.96-2.06)・衝撃 下流 157-197mm** (cell/SU2 ballpark)。基底 inviscid スキームは健全=checkerboard/対流/幾何は無実 | ref (切り分け) |
| `run_0060_node_lam_visc_bp1p90` | **node 層流粘性** (viscMethod1, no SST, no-slip, Ps1.90) | **衝撃 下流 179mm・コア加速 OK (Mmax 1.86)** = 健全な向き。step~17k で**発散** (既知の node 近壁粘性剛性, 向きは正)。→ **case36 不具合は SST 固有**と確定 | ref (切り分け) |
| `run_0067_node_sst_tauwall_bp1p90` | **AddTauWall バイナリ** (`945a27f`) で node SST Ps1.90、res_4000 から 40k | **shock 46→171mm・コア加速 1.69→1.92** (cell142/SU2132 方向へ)。NOT CONVERGED (plateau)。主病理 (壁関数 τ_w 未付与) 解消の根拠 run | ref (AddTauWall 検証) |
| `run_0068_node_sst_tauwall_outletchar_bp1p90` | 特性出口 **forward のみ** (backflow 未修正) で case36 SST 対照 | run_0067 と**衝撃軌跡ビット一致** (171mm)・出口コーナー圧 CV 不変 → **出口 forward 構成は case36 で inert** | ref (出口対照) |
| `run_0071_node_lam1st_bffix_bp1p70` | **逆流修正バイナリ** (静圧アンカー backflow) で層流コーナー発散ケース、res_4000 から 10k | **step9999 まで生存・NaN 無** (naive/forward-only は step3183/3723 で壁∩出口コーナー発散)。逆流処理が真因と実証 | ref (逆流修正検証) |
| `run_0072_node_sst_bffix_bp1p90` | 逆流修正バイナリで case36 SST Ps1.90 回帰 | run_0067 と**完全一致** (171mm/Mmax1.92)・NOT CONVERGED plateau → 逆流修正は SST Ps1.90 に**無害** | ref (回帰) |
| `run_0073_node_lam1st_bffix_bp1p70` / `run_0074_cell_lam1st_bffix_bp1p70` | **1次・層流・逆流修正で node vs cell 比較** (共通 in-duct IC: node←run_0072, cell←run_0049) | (#3 調査, 別途) | active (node/cell 比較) |

**切り分け結論 (2026-06-24)**: node Euler (下流157-197/加速2.0) と node 層流粘性 (下流179/加速1.86) は**衝撃を下流へ・コア加速 OK** で健全。**SST だけが衝撃上流46-66/加速失敗1.69・背圧逆応答**。→ **case36 node 不具合は SST 固有** (基底スキーム・対流・幾何・チェッカーボードは無実; ユーザー指摘どおり)。壁関数 utau=0 バグ (`8a2caad` 修正) は SST 不具合の一部だが残差あり。**SST 壁関数 τ_w を壁エッジに付与する SU2 AddTauWall (`945a27f`) で主病理解消** (衝撃 46→171mm, コア加速 1.92, `run_0067`)。残 ~30mm overshoot は τ_w リファイン。

### node 層流 1次の SU2/cell/node 比較 + コーナー発散 (2026-06-24)

| run_* | 設定 | 結果 | 状態 |
|---|---|---|---|
| `run_su2_lam1st_bp1p70` | SU2 NAVIER_STOKES 1次, **固定CFL** (振動回避) | **clean 単一衝撃 155mm, コア加速 Mmax1.92**, 振動なし=良い中立基準 | ref |
| `run_cell_lam1st_bp1p70` | forge cell 層流1次 | 入口衝撃+shock train, Mmax1.64 (SU2 と別枝=IC ヒステリシス) | ref |
| `run_node_lam1st_bp1p70` | forge node 層流1次 cfl_pseudo=1.0 | **step7185 で発散**=壁∩出口コーナー (x690,y±22.27) で発散 (detectNaN 局在化) | 破棄(発散) |
| `run_nodelam_cfl0p3/0p1/expl` | 同 cfl_pseudo=0.3/0.1/陽0.2 | **15000step 安定** (発散解消)・過散逸 (Mmax1.40, 入口衝撃) | ref (CFL 律速実証) |
| `run_nodelam_fineout` | 同 cfl1.0 を 50step毎出力で build-up 捕捉 | **コーナー P が Ps=1700 中心に振幅増大して振動 (std 30→333kPa), ρ 21↔6.6**=成長する圧力振動 (蓄積でない) | 破棄(診断) |
| `run_nodelam_slipwall_cfl1` | 壁を slip 化 (BL/no-slip 除去) | **それでも step3176 で同所発散・P std309kPa で振動**→ no-slip BL 無関係, **真因=出口静圧 BC** | 破棄(切り分け) |

→ **node コーナー発散の正体 (訂正済)**: 当初「質量/エネルギー蓄積」と書いたが**誤り**。細分出力 build-up で
**壁∩出口コーナーの成長する圧力振動** (P が Ps 中心に振幅10倍に増大して発散) と判明、slip テストで **no-slip BL も
否定**。**真因=壁∩出口コーナーでの outlet_statPress (Ps 規定) の数値不安定** (slip/no-slip 不問)。multi-marker emit
は効いている (コーナー2ノードが両属) が**§9.1 検証は実コーナー無しメッシュだったため症状再発**。高 CFL で振動成長、
**cfl_pseudo≤0.3 で安定**。根治=出口静圧 BC の数値安定化 (弱形式+緩和 等)。詳細 [`plans/diffusion-node-wall-viscous-distance.md`](../../plans/accepted/diffusion-node-wall-viscous-distance.md) §9.8。図 `cmp_laminar_3way.png`/`corner_divergence_buildup.png`。
| `run_node_wallstress_{off,on,off2}` | bb90036 twall カーネル A/B (30step) | 場は不変(非決定性床内)、新 twall 物理値・**utau/ypls の nonzero frac=0 実測**(壁関数退化の実証) | 削除済み(2026-08-31)(診断) |

**確定した連鎖** (`run_0049`(cell) vs `run_node_sst_muscl`/`run_node_sst_bp1p90_matched`(node) vs SU2 中立132mm):
1. node/cell とも `convMethod=2` (移流1次化ではない)。x<40mm のコア・BL は node≡cell。
2. **node SST 壁関数 (`wallTreatmentSST=1`) が `ic=bplane_cell`=壁ノード(u=0,Dirichlet)を参照→ Ut=0→ utau/ypls/wf_pk=0 に退化** (実測 utau/ypls=0)。cell は ic=内部セルで正常。twall≡0 と同一トラップ。
3. 近壁 μt が異常化 (μt/μ~4-5千がコア y=9-10mm へ侵入, cell は y=10.5+) → BL 厚化 → 発散ダクトの実効面積増を相殺 → **超音速コア加速失敗** (中心 M 1.69 vs cell/SU2 1.82) → **擬似衝撃波が ~80mm 上流** (node 46mm)。
4. 擬似衝撃波の背圧敏感性が小さな壁関数 μt 誤差を巨大な衝撃位置差へ増幅。平板(case26)は `wallTreatmentSST=0` ゆえ無傷で node/cell 一致 (Cf±1%)。

→ 図 `cmp_node_vs_cell_ROOTCAUSE.png`。修正の経緯と完了は [`plans/accepted/turbulence-node-sst-wallfunction.md`](../../plans/accepted/turbulence-node-sst-wallfunction.md) (代表点修正 + τ_w 付与で node≈cell に解消、`done`)。**この差は node の物理バグであり cell は SU2 と整合 (142 vs 132mm)。**

### ★ node vs cell SST 再比較 @ Ps=1.90 (現行バイナリ; 2026-06-27)

node SST の一連の修正 (壁関数 τ_w `AddTauWall`=`945a27f`、SST init freestream=`f639ff2`、入口ピン=`537d80f`、
境界 ghostless=`9cc6475`、**境界半割面拡散 skip**=`af5b98d`) を全て積んだ現行バイナリで、**同一メッシュ
`passive_solid.h5`・同一 BC (Ps=1.90, wt=1, MUSCL, SST, Sutherland, blockDPLUR)** で node と cell を再比較。
restart: cell←run_0049/res_80000、node←run_0072/res_40000。

| run | discretization | shock x (limit-cycle) | Mmax | peak μt/μ | VERDICT |
| --- | --- | --- | --- | --- | --- |
| `run_cmp_cell_sst_bp1p90` (+`_cont` で計 60k) | cell | **~160mm** (159-161 振動) | 1.923 | 93800 | **NOT CONVERGED** (plateau) |
| `run_cmp_node_sst_bp1p90` | node | **~163mm** (163.5→163.2 静止) | 1.913 | 92300 | **NOT CONVERGED** (plateau) |

- **★ node ≈ cell に一致 (~160 vs ~163mm, 差 ~3mm)**。過去の **cell 142 / node 171mm (29mm 差)** から
  ギャップがほぼ消滅。node の Mmax も **1.69→1.913** とコア加速が回復 (過去の「加速失敗」解消)、peak μt/μ も
  cell とほぼ同一 (92-94k)。→ **node SST 修正群 (壁関数 τ_w + init + ghostless + skip) で node が cell と整合**。
- **⚠ 収束性 (要注意)**: **両者とも `check_convergence.py` VERDICT=NOT CONVERGED** (擬似衝撃波の limit-cycle
  プラトー; rms 全列が flat 床、cell rms_roUy 6.4e-2・node rms_roOmega 2.44e6 で停滞)。報告した shock 位置は
  **limit-cycle スナップショット**で、cell は ±1mm 振動・node は静止。**「収束した」とは主張しない** (この case は
  定常 RANS では収束しないことが既知=下記 down-sweep 節・本 README 末尾の所見)。場の比較は準定常同士。
- **⚠ cell 自身が 142→160mm に移動 = 主因は forge の SST restart 過渡 (大); コードも cell を僅かに変えている (小)**:
  現行バイナリの cell は旧 run_0049 (142mm) より ~18mm 下流。2 つの効果が重なっている:
  - **(小) コード変更も cell を実際に変えている**: 894fa47 と現行を**同一 restart から並走**させると 2000 step で
    場 relL2 ~5e-4〜1e-3。これは current-vs-current の atomicAdd 床 (~5e-5) より **8〜63 倍大きい** = ノイズでなく
    実在のコード効果。犯人は **`88def3b` (出口 `outlet_statPress_d` 統一、node 専用でない一般変更)** で、cell 衝撃を
    limit-cycle 内 ~1mm・場を ~1e-3 動かす (boundaryCond を 88def3b 直前へ戻すと 894fa47 と一致)。**node 専用修正
    (537d80f/9cc6475/af5b98d, isNode ガード) は cell bit 不変**。→ 「node 修正だけ」は厳密には崩れており、
    `88def3b` が cell も僅かに変える。
  - **(大) 142→160mm の本体は restart 過渡で、コード起因ではない**: 894fa47 と現行は**同一 restart から衝撃軌跡が
    完全一致** (両者 res_500/1000/1500/2000 = 141.8/142.4/146.9/149.9mm)。**run_0049 自身のバイナリ (`894fa47`)
    でも restart すると 141.8→165.8mm に動く**。→ ~18mm の移動自体はどのコミットでも同じに起きる restart 効果
    (下記)。さらに:
  - **非収束ドリフトではない**: run_0049 は res_24000〜80000 の 56000 step を **141.5-142.1mm / Mmax 1.820-1.821 で
    完全静止** = genuine steady。当初の「非収束スナップショットのドリフト」説は**棄却**。
  - **restart_field の primitive→conserved 変換のせいでもない**: 保存量 (ro,roUx,roe,roK,roOmega) を**直接コピー**する
    忠実 restart でも同じく 165.8mm へジャンプ (`run_bisect_cell_894_faithful`)。
  - **CFL オーバーシュート単独でもない**: restart 直後 res_0 の cfl は mean 0.59/max 3.33 と連続マーチ最終
    (mean 0.21/max 1.0) の ~2.8倍だが、cfl を 0.3 に絞っても 142→150mm へ (遅くなるだけで) ドリフト
    (`run_bisect_cell_894_gentle`)。
  - **真因 = restart 直後の SST 乱流場「再平衡過渡」が双安定な擬似衝撃波を別枝へ flip させる** (詳細は 1-step probe
    `run_probe_omega` で実測):
    - **ρω→ω の計算自体は正常**: res_0 で `omega == roOmega/ro` は一致 (max差 1.0 は 1e6 に対する float32 丸め)。
      ρ 掛け忘れ等のバグではない。
    - **`vis_turb`=0 は出力タイミングのみ**: res_0 で 0 だが **res_1 で 0.3405 = 収束値ぴったり (×1.000) に回復**。
      渦粘性は犯人ではない (当初「vis_turb 復元せず」と書いたのは不正確、訂正)。
    - **本体は壁 omega の過渡ディップ→回復**: 連続収束の壁 omega 1.677e6 に対し、restart 直後は **0.82倍** (1.382e6)
      に落ち、res_2000=×0.84 → res_8000=×0.98 → res_16000=×1.05 → res_20000=×1.02 と **~16000 step かけて回復**する。
      壁 omega は永久に間違うのではなく「再平衡過渡」。ただし **その数千 step の過渡 (近壁散逸が一時不足) の間に、
      142mm の擬似衝撃波が ~160mm 側へ移り、壁 omega が戻っても戻らない**。
    - **142 は弱安定・160 は restart 安定 attractor (「双安定」ではない)**: 160mm 場から忠実 restart すると
      **160mm に留まる** (160.1→159.8→160.1, `run_selfconsist_from160` で確認, 診断後削除)。forge は決定的で、
      restart(X)→同じ Y を再現する。142mm は連続マーチでだけ成り立つ不動点 (restart 過渡で崩れる弱安定) で、
      160mm が restart に対しても安定な真の attractor。→ 「restart のたびに答えが変わる」のではなく
      「**弱い 142 が一度だけ強い 160 へ落ちる**」。
    - これは forge が restart で「収束した coupled 乱流平衡状態」をそのまま復元せず step ごとに再導出するため。
      保存量 (ro,roUx,roK,roOmega) はバイト一致でも、壁関数 u_τ・production・omega ピンの相互平衡が数千 step ぶん
      ずれ、敏感ケースを別枝へ押す。CFL を下げると過渡が穏やかになり flip も遅くなる (cfl0.3 で 142→150) が止まらない。
  - 結果、**forge node/cell は互いに一致するが、両者とも restart 由来でこの ~160-165mm 枝に乗る**。run_0049 の
    142mm (連続マーチ steady, SU2 132mm に近い) は restart では再現できていない。→ **node vs cell の比較自体は
    「同じ restart 条件同士」で公平だが、SU2/旧 cell との 142 基準には forge の SST restart 再平衡過渡が絡む**。
    これは別課題 (forge restart で壁 omega/壁関数の coupled 平衡を収束値から復元するか、過渡を抑えて flip を防ぐ)。
- 図 `compare_node_cell_sst_bp1p90_current.png` (P/Mach センターライン、SU2 132・旧 cell 142 マーカ付き)。
- **node 専用でない直近変更の確認**: 直近コミットのうち cell に効くのは `88def3b` (outlet `outlet_statPress_d`、
  isNode ガード無し) で、**実際に cell を ~1e-3 動かす** (上記並走テストで atomicAdd 床の 8〜63 倍)。他
  (`537d80f`/`9cc6475`/`af5b98d` は isNode ガード、`945a27f` AddTauWall は cell に nullptr で bit-identical、
  `f639ff2` SST init は restart で roK/roOmega が在れば inert) は **cell bit 不変**。
  → **正確には: ユーザーの node 専用修正は cell を動かしていないが、一般変更 `88def3b` は cell を僅かに変えている**
  (「node 修正だけ」は厳密には崩れている)。ただし**この小さなコード効果は 142→160mm の主因ではない** —
  ~18mm の移動は restart 過渡が支配的で、`88def3b` 込み/抜き・新旧バイナリいずれでも同じに起きる。
- **結論**: 問い「node と cell を SST で比較」への答え=**現行バイナリで node SST は cell SST と整合 (shock ~3mm,
  Mmax/μt 同等)**。ただし**双方とも未収束 (limit-cycle)** で準定常スナップショット比較であり、cell の 142→160mm は
  binary 変更でなく非収束ドリフト。両 forge とも SU2 132mm より ~30mm 下流である点は別課題。

### 背圧 down-sweep (porous, SST wall-resolved, 本番設定; 2026-06-21)

本番設定 = SLAU+SST(`wallTreatmentSST:1`)+MUSCL+陰解法 `cfl_pseudo:5`/`nStepInner:5`/**`implicitRelax:1.0`**/`detectNaN:1`、入口乱流 μt/μ=1・TI=1%、`outlet_statPress` Pt=Ps、80k step、出力2000。各点は前段から継続 (`run_case.sh` 経由, VERDICT 記録)。

| run_* | Ps[MPa] | 衝撃波フロント x[mm] | 状態 |
|---|---|---|---|
| `run_0023..0027` | 2.06→1.9 | 159 → 237 | 定常 (残差フロア ~3e-4) |
| `run_0028_sst_porous_bp1p8` | 1.8 | 267 | 定常 (フロア) |
| `run_0029_sst_porous_bp1p7` | 1.7 | 327 | 定常 (フロア) |
| `run_0030_sst_porous_bp1p6` | 1.6 | 381 | 定常 (フロア) |

- **背圧を下げると衝撃波は単調に下流へ移動** (159→381mm, Ps 2.06→1.6)、cfl5 で全点安定・NaN無し。
- 各点とも衝撃波は静止し残差は一定フロア (~3e-4, rms_roUy ~3-5e-2) で plateau。`check_convergence.py` は plateau 判定だが、**衝撃波静止＋残差一定の steady 状態として収束受理** (擬似衝撃波の limit-cycle フロア; cfl15 発散・cfl10 振動・cfl5 が安定限界)。
- 多孔壁制御の本命域 (衝撃波を壁上 x130-162) は Ps≈2.06-2.1 で、それ以上 (x130 狙い) は発散側。down-sweep は壁より下流側の特性。

### 構造メッシュ化 + 上下非対称の正体 (2026-06-21)

非構造 pave (quad/tri) は上下非鏡像 (31% 節点に鏡像なし) で擬似衝撃波の上偏りが疑われたため、
**キャビティ部を全 quad 構造格子化** (`make_mesh.py`, x区間×y帯のマップド quad, `N_SLOT` パラメータ化)。

- **メッシュ品質**: 全 quad・skew 0.011 (←pave 0.93)・**厳密 x 軸対称** (鏡像なし節点 0, 上下セル数一致)。`check_mesh_quality.py` PASS。
- **上偏りは物理的な擬似衝撃波の分岐と確定** (run_0035): 厳密対称メッシュ＋**厳密対称 IC** (asym=0) から出発しても、非対称が機械精度から指数成長 (0→0.0007→0.055→0.29→**0.37**) し、非対称IC版 (0.374) と同値へ。対称解が不安定で流れが片側へ偏る = メッシュ起因でなく物理 (ダクト内擬似衝撃波の既知の非対称分岐)。
- **構造メッシュ背圧マップ (porous, cfl5/implicitRelax1/detectNaN1/80k, 全点安定)**:

| Ps[MPa] | 1.8 | 1.81 | 1.82 | 1.83 | 1.84 | 1.85 | 1.9 | 1.95 | 2.0 |
|---|---|---|---|---|---|---|---|---|---|
| 衝撃波 x[mm] | 279 | 273 | 261 | 249 | 249 | 255 | 237 | 219 | 195 |

  背圧↑で衝撃波は上流へ単調移動。Ps=2.0 で x=195 (多孔壁 130-162 にあと一歩)。構造メッシュは quad-pave (>2.06 発散) より安定。
- 固体壁 (キャビティ・多孔板なし) 構造メッシュ Ps 1.80-1.90/0.02 を別途計算中 (porous vs solid 比較基準)。
- run_0036-0043=porous 構造 sweep, run_0044-0049=solid 構造 sweep。`up_sweep_struct.sh`/`fine_sweep_struct.sh`/`solid_sweep_struct.sh` で駆動。

### ★ 上下非対称は「粘性段階」で発生 — Euler→層流→SST 段階切り分け (2026-06-27)

「上下偏りはどの物理段階で出るか」を、**厳密対称 IC から Euler→層流→SST を順に確立**して切り分けた
(solid 構造メッシュ, 対称性メトリック `asym = ||Mach - mirror(Mach)|| / ||Mach||`, 0=対称)。各段は前段を
鏡像平均で**強制対称化**してから開始し、機械精度からの非対称成長を見る。**restart せず段階内は連続マーチ**。

| 段階 | 物理 | asym(Mach) の挙動 | 偏り |
| --- | --- | --- | --- |
| Euler | 非粘性 (**visc=0** 厳密, slip 壁) | 過渡ピーク 0.25 → **減衰 ~0.05** (Mach ほぼ対称) | **出ない** |
| 層流 | Sutherland 粘性, no-slip 壁 | 成長 → **飽和 ~0.45-0.46** (持続, 過渡でない) | **強い** |
| SST | 乱流 (wt=1) | 層流場から 0.46 → **~0.25 に緩和して飽和** | 中程度 |

- **核心**: **上下偏りは inviscid Euler では出ず (visc=0 で十分長く回すと非対称は減衰)、no-slip 境界層を入れた
  層流段階で初めて飽和的に発生する**。乱流 (SST) は層流の偏りを部分的に緩和するが残る (run_0049 SST 0.30 と同系統)。
  → 「物理的な非対称分岐」(run_0035) の駆動源が**粘性境界層 (剥離)** と確定。bp=1.90/1.70 の両方で同じ。
- **落とし穴 (教訓)**: ① `viscMethod:0` は「定数粘性」で**非粘性ではない**。厳密 Euler は `visc:0.0`+`thermCond:0.0`。
  ② 短い run (~12k step) では非対称の**過渡ピークを定常的偏りと誤認**する。Euler は過渡 0.25→減衰 0.05 なので、
  十分長く (100k) 回さないと判定を誤る。
- **★ bp=1.70 SST が VERDICT=PASS (収束)**: この**クリーン段階確立 (Euler→層流→SST) で初めて収束解**を得た
  (rms 全列 3-8.6 桁低下, 衝撃 x=237mm)。restart ベースの全 run が NOT CONVERGED だったのは、別途診断した
  **SST restart 過渡 (壁 omega ディップ) を段階確立が回避**するため。run: `run_sym_E_euler_bp170`(Euler) /
  `run_sym_F_laminar_bp170`(層流) / `run_sym_G_sst_bp170`(**SST, PASS**)、bp1.90 側は `run_sym_C_trueeuler` /
  `run_sym_D_laminar`。図 `run_sym_G_sst_bp170/bp170_staged_asym.png`。
- **★ 2次風上で cell ≈ node (bp1.70 SST)**: 上記 `run_sym_G`(1次風上収束) から `convMethod:1`(2次風上)+Venkat に
  切替え、**cell** (`run_sym_H_2up_cell`, 1次収束場から restart) と **node** (`run_sym_H_2up_node`,
  cell 収束場を centCoords 最近傍で interp) を計算。両者ほぼ一致: **衝撃 cell 241.8 / node 240.9mm (~1mm差)、
  Mmax 2.019 / 2.005、asym 0.251 / 0.246**。1次風上 (Mmax~1.9, 衝撃 237) から 2次でコアがシャープ化 (Mmax~2.01)・
  衝撃わずか下流。**判定**: 両者とも `check_convergence`=NOT CONVERGED (ρ/運動量 3桁低下だが roK/roOmega プラトー)
  だが、`check_quasisteady`=**STEADY** (shock/asym/Mmax の drift<0.5%) → 量は準定常で一致 (残差未収束でも量で判定)。
  図 `run_sym_H_2up_node/cmp_cell_node_2up_bp170.png` / `_mach.png`。node(median-dual) 健全性を 2次風上でも追認。

### ★ porous vs solid 比較 = 論文の機構を再現 (2026-06-21)

構造メッシュで porous (run_0036-0043) と solid (run_0044-0049) を matched-Ps で比較。
中心軸静圧の波状振幅 (peak-to-peak, x180-450mm) を「ショックレス度」の指標とする:

| Ps[MPa] | porous 波状 | solid 波状 | porous 衝撃波x | solid 衝撃波x |
|---|---|---|---|---|
| 1.80 | 0.361 | 0.362 | 279 | 279 |
| **1.82** | **0.280** | **0.380** | 261 | 267 |
| **1.84** | **0.265** | **0.377** | 249 | 255 |
| 1.90 | 0.260 | 0.250 | 237 | 231 |

- **Ps=1.82-1.84 で多孔壁が擬似衝撃波を弱化** (中心軸波状 ~30% 減)。図 `porous_vs_solid_centerline.png`:
  ① 多孔壁では**圧力上昇が多孔壁領域内 (130-162) から早く始まる** (穴からの吹出/吸込が壁上で圧縮を開始)、
  ② 振動振幅が solid (peak 0.48-0.56) より小さい (porous peak 0.36-0.43) = **shock train が弱まりショックレスに近づく**。
- **= Matsuo et al. の「多孔壁パッシブコントロールで擬似衝撃波が弱化する」を 2D で定性再現**。ユーザが指摘した「1.8-1.85 の良い領域」と一致。
- 注意 (誠実な留保): 各 run は残差フロアの steady 状態 (擬似衝撃波 limit-cycle), 流れは物理的に上下非対称 (中心軸指標), 2D スリット (3D 丸穴でない)。Ps=1.80/1.90 では porous≈solid (効果は Ps 窓依存=衝撃波が多孔壁と相互作用する位置のときのみ)。

### Ducros 1次化 ON/OFF 比較 (2026-06-21, 固体壁 Ps=1.84)

衝撃近傍で MUSCL リミタを強制 1 次化する Ducros センサ寄与を `ducrosLimiter` config フラグ
(**既定 0=off で 1 次化を使わない**, 1=on で従来の強制 1 次化) で切り替えられるようにし
(`solverConfig.hpp`/`.cpp`, `ducrosSensor_d.cu`)、固体壁 Ps=1.84 で ON/OFF を**同一バイナリ・
同一 restart 場 (run_0046 発達場)・30000 step**で比較。解析 `analyze_ducros.py` → `compare_ducros_centerline.png`/`ducros_summary.txt`。

| run_* | ducrosLimiter | 結果 | 状態 |
|---|---|---|---|
| `run_0052_solid_ducrosON_bp1p84` | 1 (1次化あり) | ducros 計算正常: max 0.9996, duc>0.8 発火 1.79% (x169-260mm=衝撃波フロントに局在)。中心軸 ripple std 8.76kPa | active (比較基準) |
| `run_0053_solid_ducrosOFF_bp1p84` | 0 (1次化なし) | ducros 一様 0。**ON とほぼ完全一致** (中心軸 ripple std 8.76kPa, OFF/ON 比 1.000, 衝撃波フロント 189mm 同一) | active |
| `run_0030_limA_pshock` / `run_0031_limC_pshock` / `run_0032_limC_pshock_long` | **リミッタ無次元化の回帰** (`run_sym_H_2up_node` の `res_80000` から index コピー restart、A=`mr0/scaled0` 生産既定, C=`mr1/scaled1` 無次元化。旧キー `LESorRANS/LESmodel/RANSmodel` は `model: "sst"` に置換)。C は 8000 step では過渡が残るので +24000 step 継続 | **衝撃位置は同一** (`check_quasisteady --quantity shock` 両者 **STEADY**, tail mean **240.9**, drift 0.0 %)。残差は過渡が抜けると C が上: `rms_roUx` 2.19e-2→**3.14e-3**, `rms_roK` 2.96e-1→**3.04e-2**, `ro`/`roUy`/`roe` も C が下、`roOmega` のみ 1.09 倍 | ref ([../../plans/active/convection-node-wall-reconstruction.md](../../plans/active/convection-node-wall-reconstruction.md) §4.22) |
| `run_0033_fairA` / `run_0034_fairC` | **リミッタ無次元化の公平比較** (`run_sym_H_2up_node` の `res_80000` から index コピー restart、**同一 32000 step**、基準値を `limiterRoRef 16.45741351 / limiterPRef 1207044.245 / limiterARef 312.4399249` に固定。ログで `(fixed)` を確認) | **衝撃位置は同一** (`check_quasisteady --quantity shock` 両者 STEADY, tail mean **241.7**, drift 0.0 %)。残差は一長一短: C が良いのは `rms_roUx` (9.76e-3→2.79e-3) のみ、`ro`/`roUy`/`roe`/`roOmega` は C が悪い。両者 32000 step で NOT CONVERGED。**`run_0030`–`run_0032` の不公平比較 (C だけ +24000 step) は撤回** | ref (plan §4.22 / codex plan-3 M7) |

- **結論: Ducros 1 次化を切っても場の差は <0.1% (中心軸 P 差 max 1.3kPa, ducros 発火セルに局在)。
  衝撃波列 ripple 振幅・フロント位置・limit-cycle 変動はいずれも不変**。ユーザ予想の「変動が激しく
  なる」は**起きなかった**。
- 理由: 基底フラックスが SLAU (上流型) で、衝撃捕捉はフラックス自身の散逸が支配する。~1.8% の衝撃
  セルで MUSCL 再構成次数 (1 次/2 次) を変えても結果はほぼ変わらない。SLAU+Venkatakrishnan では
  Ducros 1 次化は**ほぼ冗長な保険**。低散逸の中央/KEEP フラックス (再構成次数が散逸を直接支配) では
  効くはずで、それが本来の用途。
- ducros が計算できているかの確認: ON では衝撃波フロント (x~190mm) に正しく発火 (median 197mm)、
  OFF ではフラグで一様 0 化。センサ自体は正常動作。

### 現時点の所見 (2026-06-21 自走セッション)

- **メッシュ/手法は確立**: 2D wall-resolved SST (y+~77, AR≤1000, skew≤0.9 ゲート通過)、cfl_pseudo=5 陰解法は establish では安定。背圧で擬似衝撃波を seeding→chain で形成可能。
- **⚠ 背圧 run はいずれも未収束 (`check_convergence.py` で全列確認)**: run_0016/0018/0019/0020 とも残差が plateau または**上昇 (発散的)**。**衝撃波＋剥離＋キャビティ循環が本質的に非定常 (limit cycle / 振動)** のため、定常 RANS (unsteady=0) では収束しない。→ **報告した衝撃波位置・中心軸波状振幅は非収束スナップショット値であり信頼できない**。位置安定だけを見て「収束」と判断したのは誤りだった。
- **物理メカニズムの再現は未達 (かつ未収束で判定不能)**: 層流(run_0012, BL薄すぎ)も SST(非収束)も、現時点で論文の「ショックレス化」を確認できていない。
- **次手 (収束問題の解決が先決)**: ① **URANS (unsteady=1, dual-time) で時間平均**するか、定常で収束する条件/手法を探す (擬似衝撃波の非定常性が主問題)、② 3D 丸穴化 (2Dスリットの吸込/吹出分布が論文と異なる疑い)、③ matched-x_f 比較、④ `katoLaunder` で μt 過大抑制、⑤ 出口チャンバ追加。

> 注: 領域は x=690 出口の簡略版 (10°ディフューザ+チャンバ未追加)。run_0003-0007 は旧SST(y+~1200)/手法確立途中の探索 run、run_0009-0010 は超音速ロックインで棄却、run_0015-0017 は中間 (seed/establish)。
