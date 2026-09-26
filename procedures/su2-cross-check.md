# forge ↔ SU2 クロスチェック手順

forge の結果(特に軸対称・乱流・近軸挙動)を独立した検証済みソルバ SU2 と突き合わせるための標準手順。
**目的**: forge 固有のスキーム/離散化/数値精度の問題を、同一メッシュ・同一境界条件で SU2 と比較して切り分ける。

> このファイルは [`AGENTS.md`](../AGENTS.md) から参照される。軸対称や乱流の「forge だけ変な値が出る」事象を調べるときは、
> 推測で結論づけず、まず本手順で SU2 と比較すること。

## SU2 の所在と実行

- バイナリ: `.external/su2/bin/SU2_CFD`(v8.5.0 "Harrier", linux64-omp ビルド)。`SU2_SOL` / `SU2_DEF` も同梱。
- 直接実行: `cd <run_dir> && OMP_NUM_THREADS=8 /home/sano/work/forge/.external/su2/bin/SU2_CFD case.cfg`
  - `SU2_RUN` 環境変数は直接 CFD 実行には不要(python ラッパ使用時のみ)。
  - スレッド数は `OMP_NUM_THREADS` で指定。複数 run を並行させるときはコア数(`nproc`)を超えないよう分配する。

## メッシュ生成(forge と同一形状を共有する)

**最重要**: forge と SU2 を**同一の `.geo` から**生成し、「メッシュ違い」を交絡させない。

```bash
# 同一 .geo から両方を出力(gmsh 4.x)
gmsh -2 conical.geo -o conical.msh -format msh2   # forge 用(convertGmshToForge は Gmsh 4.1 形式が必要)
gmsh -2 conical.geo -o conical.su2 -format su2    # SU2 用
```

- **単位スケール**: `.geo` 内の `Mesh.ScalingFactor = 0.001;`(mm→m)は **gmsh が .msh / .su2 両方に適用**するため、
  両者とも SI(メートル)で一致する。SU2 側は `SYSTEM_MEASUREMENTS= SI` のままでよい。
  → 生成後に `.su2` の `NPOIN` 直後の座標が forge メッシュ(`/CELLS/centCoords`)と同レンジか必ず確認する。
- **Physical 名 → SU2 マーカー**: `.geo` の `Physical Curve("inlet"/"outlet"/"wall"/"axis")` がそのまま SU2 のマーカー名になる。
  forge の `bcondConfig.yaml` の physID と対応付けて、下記の cfg マーカーに割り当てる。

## 軸対称 cfg テンプレート(case 29 conical, Pt=4MPa/Tt=1500K/背圧20kPa)

3 本(Euler / 粘性層流 / RANS-SST)を同じメッシュ・BC で回す。差分のみ示す。

共通:
```
AXISYMMETRIC= YES
SYSTEM_MEASUREMENTS= SI
FLUID_MODEL= STANDARD_AIR
GAMMA_VALUE= 1.4
GAS_CONSTANT= 287.058
INIT_OPTION= TD_CONDITIONS
MACH_NUMBER= 0.10
FREESTREAM_PRESSURE= 20000.0
FREESTREAM_TEMPERATURE= 300.0
CONV_NUM_METHOD_FLOW= ROE
MUSCL_FLOW= YES
SLOPE_LIMITER_FLOW= VENKATAKRISHNAN
TIME_DISCRE_FLOW= EULER_IMPLICIT
CFL_NUMBER= 10.0
CFL_ADAPT= YES
CFL_ADAPT_PARAM= ( 0.5, 1.5, 1.0, 1.0E6, 0.001, 50 )
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
MARKER_SYM= ( axis )                                   # 軸 = 対称面 (forge: kind: slip/axis)
INLET_TYPE= TOTAL_CONDITIONS
MARKER_INLET= ( inlet, 1500.0, 4000000.0, 1.0, 0.0, 0.0 )   # Tt, Pt, 流れ方向 (x,y,z)
MARKER_OUTLET= ( outlet, 20000.0 )                     # 背圧 Ps (超音速流出時は無視され外挿)
MESH_FORMAT= SU2
MESH_FILENAME= conical.su2
TABULAR_FORMAT= CSV
OUTPUT_FILES= ( RESTART, PARAVIEW, SURFACE_CSV )
VOLUME_OUTPUT= ( COORDINATES, SOLUTION, PRIMITIVE )
```

- **Euler**: `SOLVER= EULER` / `MARKER_EULER= ( wall )`(壁=スリップ)。
- **粘性層流**: `SOLVER= NAVIER_STOKES` / `MARKER_HEATFLUX= ( wall, 0.0 )`(断熱 no-slip)。
  `VISCOSITY_MODEL= SUTHERLAND` / `MU_REF= 1.716e-5` / `MU_T_REF= 273.15` / `SUTHERLAND_CONSTANT= 110.4` / `PRANDTL_LAM= 0.72`。
  (forge `viscMethod: 1` = Sutherland に対応)
- **RANS-SST**: 粘性層流に加えて `SOLVER= RANS` / `KIND_TURB_MODEL= SST` /
  `FREESTREAM_TURBULENCEINTENSITY= 0.01` / `FREESTREAM_TURB2LAMVISCRATIO= 10.0`。
  - 入口乱流: forge の `k`/`omega` 直接指定と、SU2 の(乱流強度 + 粘性比)は厳密一致しない。
    目標 μt/μlam を合わせる(本ケースは ~10)。絶対量は多少ずれても、**軸中心 k の「スパイク有無(形状)」の比較**には十分。

## 収束確認(必須・[AGENTS.md] の収束ルールに従う)

SU2 のノズル流は **リミットサイクル**で残差が下げ止まることが多い。`rms[Rho]` だけでなく:
- `history.csv` の `rms[Rho]`, `rms[Momentum-Y]`, `rms[TKE]`, `rms[Dissipation]` の**トレンド**を見る(単調減 or 振動プラトーか)。
- 積分量(出口 massflux: `MARKER_ANALYZE_AVERAGE= MASSFLUX`)が定常化しているか。
- 早期の `rms[Rho]` だけ低くても、運動量・乱流残差が下げ止まり/上昇していれば**未収束**と判断する。
- **`rms[RhoE]`(エネルギー)が最も遅い**。`rms[Rho]` が数桁落ちても `rms[RhoE]` が 1 桁未満・`rms[w]` が正のままなら未収束。
- **`su2.log` の `Exit Success` を収束と誤認しない**。セッション中断時の SIGTERM でも SU2 は graceful 終了し
  `Exit Success` + restart を書く。`history.csv` の **iter 数が `ITER` 上限/収束基準に達したか**を必ず確認すること
  (実例 2026-06: nozzle 粘性/SST が 975/752 iter で中断され `rms[RhoE]≈-0.2` のまま「完走」扱い → forge と
  下流で偽の 20-25% 差。`RESTART_SOL=YES`+`READ_BINARY_RESTART=NO`+`SOLUTION_FILENAME` で途中解から継続し
  `rms[RhoE]≈-1.4`・出口積分量ドリフト <0.2% まで発達させると **forge と全域 ≤0.8% 一致**した)。

### `CFL_ADAPT` 共振による偽の limit cycle に注意(衝撃列/擬似衝撃波)

衝撃列 (shock train) や擬似衝撃波 (pseudo-shock) で **残差が 2 桁オーダーで激しく振動して全く収束しない**ときは、
**物理的非定常と決めつける前に `CFL_ADAPT= NO` + 固定 `CFL_NUMBER`(例 2.0)を試す**こと。`CFL_ADAPT` の ramp
(`CFL_ADAPT_PARAM=(0.5,1.2,…)` 等)が残差の増減と共振し、**スキームと無関係に**限界振動を作ることがある。

- 切り分け手順 (実例 2026-06 case 36 擬似衝撃波): リミッタを 3 段階 **Venkatakrishnan → Van Albada edge → 1次(MUSCL OFF)** に
  変えても `CFL_ADAPT= YES` 下では全部 `rms[Rho]` が −0.9↔−2.8 で振動。**`CFL_ADAPT= NO`+`CFL_NUMBER=2.0` にした瞬間、
  1次でも2次でも振動が消え準定常プラトー (`rms[Rho]~−3.2`) に落ちた** → 振動は CFL adapt 由来 (物理でもリミッタでもない)。
- したがって振動を見たら **(1) 固定 CFL で再開 → 準定常に落ちるか** を最優先で確認する。落ちれば CFL 共振、
  固定 CFL でも振動継続なら物理的非定常 (URANS or limit-cycle 時間平均で評価)。
- 背圧固定の擬似衝撃波は固定 CFL でも `rms[Rho]` が −3 付近のプラトー止まりが普通。**衝撃位置が静止し massflux が
  定常**なら準定常スナップショットとして比較してよい (forge 側も同様にプラトー)。

## ライン比較プロトコル(node / cell / SU2 三者・必須)

**離散化・対流・拡散・境界など forge の数値挙動を検証/比較するときは、forge を node と cell の両系統で回し、
さらに同一形状・同一 BC の SU2 とあわせて「サンプリングライン上の物理量」で三者比較する。** 残差や単一スナップショットの
目視ではなく、**ライン上の P / T / 速度**を SU2 基準で突き合わせることを比較の正本とする。

手順:

1. **同一形状で 3 ケースを用意**する。forge は同一メッシュを `discretization: "cell"` と `"node"` で、SU2 は同じ `.geo` から
   生成した真の 2D (または軸対称) `.su2` で回す。**メッシュ違いを交絡させない**(上記「メッシュ生成」)。3 ケースとも
   `check_convergence.py` で `PASS`(または妥当なプラトー)を確認してから比較する(未収束同士を「一致」と呼ばない)。
2. **サンプリングライン**は、壁・よどみ・形状擾乱を避けて流れが素直な位置に取る。目安は主流方向に沿った
   **下端から高さ ~25% 程度の水平ライン** (壁・はく離・衝撃を避け、全長で流体内に収まる高さ)。具体的なライン位置は
   形状依存なので各ケースのドキュメント側で定める。壁面量を見たい場合のみ壁沿いラインを別途取る。
3. **比較する物理量**: 静圧 `P`、温度 `T`、速度 (大きさ `|U|` または成分 `Ux,Uy`)。ライン上で x に対してプロット重ね描きし、
   SU2 を基準に node / cell の相対差を出す。
4. **判定**: 三者が許容内で一致すれば妥当。node/cell 間で差が出る場合はその量・位置を明示し、どちらが SU2 に近いかを述べる
   (cell の atomicAdd 非決定性で残差が落ちきらないケースなど、収束度の差も併記する)。
5. **成果物**: 比較図 (`*_line_compare.png`) と数値差サマリを run に残し、結論にラインと run パスを明示する。

実装の足がかり: SU2 VTU の読み取りは [`case/29.bell_vs_conical/compare_forge_su2.py`](../case/29.bell_vs_conical/compare_forge_su2.py) の
`load_su2()`、forge `res_*.h5` のライン抽出は同 case の比較スクリプトを流用する。

## VTU の読み取り(注意)

SU2 v8 の `flow.vtu` は `NumberOfComponents= "3"`(= 後の空白)など属性が非標準で、**VTK の `vtkXMLUnstructuredGridReader` が失敗する**ことがある。
その場合は appended-raw バイナリを手動パースする(`header_type="UInt64"`、各 DataArray は 8 byte 長さ接頭辞 + Float32 ペイロード、offset は `<AppendedData>` の `_` 直後からの相対)。
実装例: [`case/29.bell_vs_conical/compare_forge_su2.py`](../case/29.bell_vs_conical/compare_forge_su2.py) の `load_su2()`。

## 既知の比較結果(case 29 conical, 2026-06)

軸対称 SST の「軸中心 k スパイク」を SU2 と比較し、**forge 固有・近軸の数値問題**であることを特定:
- SU2(頂点中心・軸上に節点・倍精度)は軸中心の半径速度 `u_r` が 0 から滑らかに立ち上がり、k は軸で最小(スパイク無し)。
- forge(セル中心)は **float32 の陰解法(block-DPLUR)が近軸第一セルの `u_r` を収束させきれず**(陽解法・倍精度では正しい)、
  偽の `∂u_r/∂r` → 偽ひずみ → SST 生産で k がスパイク。フープ項・Kato–Launder は無関係(下流の対症療法)。
- 詳細: [`.github/plans/architecture-axisym-axis-singularity.md`](../plans/accepted/architecture-axisym-axis-singularity.md)。

## CHT (multizone) で突き合わせるとき (2026-09-24 追加)

plan [`boundary-conjugate-heat-transfer`](../plans/accepted/boundary-conjugate-heat-transfer.md) §6 V2 で実施した手順。
実例は [`case/52.conjugate_slab/su2/`](../case/52.conjugate_slab/su2/) (1 次元純伝導の共役スラブ)。

**構成**: 親 config が `SOLVER= MULTIPHYSICS` + `CONFIG_LIST= (fluid.cfg, solid.cfg)`、
`MARKER_ZONE_INTERFACE` と `MARKER_CHT_INTERFACE` に**両側のマーカ名**を並べる。
**`MULTIZONE_MESH= NO` にしてゾーンごとに別メッシュ**を与えるのが楽 (1 ファイルに `IZONE=` で
詰める形式もあるが、マーカ名の衝突を気にせずに済む)。固体は `SOLVER= HEAT_EQUATION` +
`THERMAL_CONDUCTIVITY_CONSTANT` + `MARKER_ISOTHERMAL` で背面を固定する。

**界面の節点は両側で一致させる**。SU2 は内挿できるが、内挿誤差を比較の差に混ぜない。

**罠: 静止流体では既定の収束判定が即座に成立して、固体が未収束のまま終わる。**
純伝導のスラブでは流体の BGS 残差が最初から **−23** で、`CONV_FIELD= AVG_BGS_RES[0]` /
`CONV_RESIDUAL_MINVAL= -12` だと **11 外反復**で「収束」して終了する。そのとき固体側は **−0.62** で、
界面温度は解析解から **1.8 K (11 % of rise)** ずれていた。
→ **判定は固体側 (`AVG_BGS_RES[1]`) に付けるか、固定回数まで回して両ゾーンの残差を記録する**。
3000 外反復で固体が −8.58 まで落ち、界面温度が解析解と **1e-7 K** で一致した。

**流体の温度は restart に `Temperature` 列として出ないことがある** (`INC_DENSITY_MODEL= CONSTANT`
では `Enthalpy` が出る)。$h=c_p(T-T_{\rm ref})$ で、$T_{\rm ref}$ は `FREESTREAM_TEMPERATURE` の既定
298.15 K。固定温度マーカの値 (既知) から $T_{\rm ref}$ を逆算して検算すること。

**熱流束 `HF[*]` は `MARKER_MONITORING` を設定しないと 0 のまま出る**。

### 乱流ケースで CHT を回すとき (2026-09-24、case/48 で実測)

1 次元スラブは通ったが、**乱流平板 (case/48, M4.19, SST) は 4 回設定を変えても定常に達しなかった**。
踏んだ順に記録する (同じ順で踏む可能性が高い)。

- **既存 cfg を流用するとキー重複で即エラー**になる (`VOLUME_OUTPUT: option appears twice`)。
- **multizone では `ITER=` が使えない**。`INNER_ITER` (ゾーン内) と `OUTER_ITER` (連成) に分ける。
  `ITER` を残すと `ITER must not be used when running multizone` で落ちる。
- **再開ファイルはゾーン番号が付く**: `SOLUTION_FILENAME= solution_flow` に対し、実際に読むのは
  `solution_flow_0.csv`。
- **出発点を揃えないと「合わない」ではなく「解けていない」を見る**。自由流から始めた回は
  界面温度が forge と **117 K** 違ったが、残差を分解すると `bgs[w]` が **+6.7 桁**で、単に未収束だった。
  **固定壁の収束解から再開し、さらにその固定壁温を連成の答えの近く (本件では 604 K) に取って暖機する**こと。
  暖機は `rms[Rho]` −12 まで 3677 反復で届いた。
- **CFL 適応の既定 `(0.5, 1.2)` は積が 0.6 < 1** で、残差が振れるたびに CFL が最小へ落ちて戻らない
  ([[su2-cfl-adapt-collapses]])。`(0.8, 1.25)` など積 ≥ 1 にする。
- **収束したかは界面の熱収支で見る**。両ゾーンの界面温度が一致していても (実測: 流体 664.314 /
  固体 664.44 K で 0.13 K 一致) 収束しているとは限らない。**界面温度の一致は連成が効いている証拠であって、
  収束の証拠ではない**。
- **固体側の熱流束を $k_s(T_w-T_{back})/t$ で逆算しない** (2026-09-25 訂正)。**2D 固体ではこの式が面内伝導を
  落とす**ので、収支の食い違いを実際より大きく見せる。実例: 逆算 79 kW/m² vs 流体 `HF[0]` 49.8 kW/m² で
  「**60 % 食い違い**」と報告したが、**SU2 自身が出す `HF[1]` は 46.7 kW/m² で差は 6.21 %**、しかも
  外反復とともに縮んでいた (`HF[0]` −11.3 / `HF[1]` **+127.5** /outer)。
  ~~収支は `HISTORY_OUTPUT` の `HEAT[0]`/`HEAT[1]` (履歴 CSV の `HF[0]`/`HF[1]`) で比べる。~~
  **訂正 (2026-09-26、CHT plan §5.1 #91/#93、codex result 8 巡目)**: `HF[0]` と `HF[1]` は**定義が違う量**なので、両者の差を連成の不釣合いとして読まない。**`HF[1]` = 固体ゾーンにかかる Robin 荷重の和 (適用荷重)**、**`HF[0]` = 流体壁節点の浮動 DOF に依存する WLS 再構成の診断流束**。forge と比べる相手は **`HF[1]` ↔ Σ`iface_Qf_eff`** (どちらも固体に実際にかけた積分済み荷重)。実例 `case/48` `su2_cht_cont` 最終値: `HF[0]` 59725.64 / `HF[1]` 60684.57 W/m で**1.58 % 違うが、これは連成不釣合いではない** (適用荷重どうしは保存則で 0.1 W/m 級に閉じる)。
- **SU2 の固体ゾーンは既定 (`INNER_ITER`=1 × `CFL_NUMBER`=1.25) だと擬似時間がほとんど進まない** (1 外反復 ≈2.8 µs)。固体の冷却過渡を這うだけで「収束しない」「発散する」ように見える。**`solid.cfg` に `CFL_NUMBER= 1000`** を入れる (`case/48` の判別 A/B で確定、CHT plan §5.1 #91)。
- **界面温度は固体の `restart_1.csv` の $y$=0 行から取る**。`AvgTemp[1]` は `CHT_WALL_INTERFACE` では常に 0 (SU2 の仕様、`CHeatSolver.cpp`)、流体壁節点の $T$ は浮動 DOF で前縁では 28 K ずれる (CHT plan §5.1 #93)。
- **「N 反復回した」の N を取り違えない** (同上)。`OUTER_ITER= 200` × `INNER_ITER= 200` は
  **流体反復 40000 だが連成更新は 200 回**である。**CHT で効くのは外反復の数**で、
  界面量の交換は外反復ごとにしか起きない。1 次元スラブでさえ固体残差を下げるのに **3000 外反復**を要した。
  履歴 CSV の `bgs[T][1]` の傾きから所要外反復を見積もってから回すこと
  (case/48 は −0.617・−0.00385 /outer で、判定値 −14 まで約 3500 外反復)。
- **固体ゾーンの再開にも `READ_BINARY_RESTART= NO` が要る**。流体側だけに書いていると
  `Unable to open SU2 restart file solution_1.dat` で落ちる。固体の再開に要るのは 3 行:
  `RESTART_SOL= YES` / `SOLUTION_FILENAME= solution` / `READ_BINARY_RESTART= NO`
  (ゾーン 1 の既定の出力名が `restart_1.csv` なので、それを `solution_1.csv` に置いて読ませる)。
