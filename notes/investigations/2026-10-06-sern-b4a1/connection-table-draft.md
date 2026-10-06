> 2026-10-06、Plan エージェントの下書き (未確定・上位の判断前)。plan `tooling-sern-mesh-blocking.md` §5.1 B4a(1)。

# B4a(1) 全体接続表 — 下書き (⑤ SERN 3D 全ヘキサ・ブロッキング)

対象リポジトリは `/home/sano/work/forge-sern-design`(読んだだけで、編集していません)。以下、行番号の参照は次の略記を使います。
- plan = `plans/active/tooling-sern-mesh-blocking.md`
- M3 = `design/forge_design/meshing/mesh_sern3d.py`
- JM = `case/46.sern_design/cad/hex_junction_model.py`
- RN = `design/forge_design/evaluate/runner_sern3d.py`
- YAML = `case/46.sern_design/problem_3d_prod_m6on_wallres_lswx08.yaml`

**生産形状の数値は自分で計算しました**: `runner_sern.design_from_problem` を `PYTHONDONTWRITEBYTECODE=1` で実行しています。純関数で、ファイル出力はありません。
- L_ramp 10.0916、y_e 2.6834、L_cowl 1.2、y_te −0.10499
- カウル後縁の傾き −0.08749 (−5.00°)。cowl_xy は直線
- θ_b(MOC 自由境界の終端角)−40.24°、θ_e −0.758°
- y_veh = max y_r + 0.06 = 2.7465、テーパ始点 x0 = 0.65 L_ramp = 6.560
- x_out = 12.092

**ここからの重要な発見**: 旧メッシャのカウル後流線は接線方向ではなく θ_b 方向に引かれています(RN:148 `interface_angle = θ_b`、M3:184-187)。その結果、旧の y_bot は y_mid(x_out) − bot_depth = −9.32 − 3 = **−12.32 H** になります(M3:187)。つまり旧の下側境界は、接続模型の YBOT = −1 H(JM:58)よりずっと深い位置にあります。

---

## 1. 固体と流体の定義

記号: z_w = W/2 = 1.0、z_o = z_w + t_sw = 1.005、t_sw = t_c = 0.005 H(plan:769-774)。z_far は **絶対値 2.50 H**(plan:790。接続模型は 2.0、JM:58)。b(x) は機体上面線(M3:403-443)。

### 1a. 固体(推奨案 V1: 機体側面を z = z_o に置き、側壁外面と同じ平面にする)

| 固体 | x | y | z | 旧メッシャとの差 |
|---|---|---|---|---|
| ランプ = 機体下面 (`ramp`) | [−L_up, L_ramp] | y_r(x)。x<0 では 1、[−t_f, x_f2] はフィレット円弧(export_contours.py:26-38) | [0, z_w] | 旧も z≤W/2 は ramp(M3:364)。新では z∈[z_w,z_o] の帯が増える(§4) |
| 機体本体 | [−L_up, L_ramp]。上流端は入口面で切れて穴になる | [y_r, b(x)]。b は x<x0 で y_veh、その後エルミートテーパで y_e+t_base | [0, **z_o**] | 旧は z≤W/2(M3:528-531)。**幅が 2 t_sw = 0.01 H 広がる** |
| 機体上面 (`vehicle_top`) | 同上 | b(x) | [0, z_o] | 旧は k<k_sw のみ(M3:490) |
| 機体側面 (`vehicle_side`) | [−L_up, L_ramp] | [y_r, b] | **z_o** 平面 | 旧は z = W/2(M3:575) |
| 機体ベース (`vehicle_base`) | L_ramp 平面 | [y_e, y_e+t_base] | [0, z_o] | 旧は z≤W/2、t_base 0.02(RN:145 既定) |
| 側壁 | [−L_up, **L_sw=0.8**]。後端は平面 | [y_c, y_r] | [z_w, z_o](+z 側に厚み) | 旧は**厚さ 0**。z=W/2 に D2 双子節点を置いていた(M3:9,276-278) |
| カウル板 | [−L_up, **L_cowl**]。後端は平面 x=L_cowl | 内壁 = y_c(x)。0≤x≤L_cowl では cowl_xy、**x<0 では 0(要定義、§6-8)**。外壁 = 内壁法線方向の外向き(下)オフセット t_c(一定) | [0, **z_o**](側壁の下まで延びる) | 旧は中間線 = 弦 −x·tan_c、厚さは t_c から 0.8 L_cowl 以降 0 へ絞る(M3:212)。z も W/2 手前 2 セルで 0 に絞り(M3:282-288)、z>W/2 には板が無い。**旧の実壁は +t/2 ずれている**(plan:767) |

z>z_o の y=y_r 面は、接続模型では `vehicle` 壁で閉じています(JM:411)。V1 ではここが流体どうしの内部面になります。これは旧の R4e と同じ扱いです(M3:366-367)。

### 1b. 流体領域

| 領域 | 範囲 | 初期場の区分 |
|---|---|---|
| ノズル内 | z∈[0,z_w]、y∈[y_c,y_r]、x∈[−L_up, x_out](N ブロック) | 排気 |
| 側壁の後流 | z∈[z_w,z_o]、y∈[y_c,y_r]、x≥L_sw(SW) | 外気(旧 RN:121 も z≤W/2 だけを排気にしている。境界上の節点の帰属は B4d で決める) |
| カウルの後流 | y∈[y_o,y_c]、z∈[0,z_o]、x≥L_cowl(CW1/CW2) | 外気(推測: 旧の dup1 に相当する扱い) |
| カウル下の外部流 | y∈[YBOT, y_o]、全 z | 外気 |
| 側方の外部流 | z∈[z_o, z_far]、y∈[y_o, b] | 外気 |
| 機体上の外部流 | y∈[b, b+top_depth]、全 z | 外気 |
| ベースの後流 | x>L_ramp、y∈[y_r,ext, y_r,ext+t_base]、z∈[0,z_o] | 外気(旧と同じ) |

外部境界:
- x=−L_up: N ブロックは `inlet_nozzle`、それ以外は `inlet_ext`(旧 M3:351、M3:498、M3:585)
- x_out: `outlet`
- y=YBOT: `bottom`
- y=b+top_depth: `top_out`(**旧は b(x) に沿うオフセット面**、M3:472)
- z=0: `sym`
- z=2.50: `side_far`

## 2. x 方向の区間分割

物理 station は次のとおり(丸めなし)。
- −L_up、−t_f、x_f2(ランプフィレット接点 ±0.013 H)
- L_sw−ZONE、L_sw、L_sw+LWAKE、(L_sw+L_cowl)/2
- L_cowl、L_cowl+LWAKE、L_cowl+ZONE
- x0(テーパ始点)、L_ramp、L_ramp+LWAKE_b、x_out

| 区間 | 流体ブロック数 | 固体のままのブロック | 接続模型との関係 |
|---|---|---|---|
| U0 [−L_up, 0) | 26 | CW1, CW2, SW, V1, V2 | **新規**。A と同じ位相。カウル内壁の折れ (x=0) と フィレットを含む |
| A [0, L_sw) | 26 | 同上 | 模型 A の 16 + 上部 10 |
| B [L_sw, L_cowl) | 27 (+SW) | CW1, CW2, V1, V2 | 模型 B の 17 + 10 |
| C [L_cowl, L_ramp) | 29 (+CW1, CW2) | V1, V2 | 模型 C の 19 + 10 |
| D [L_ramp, x_out] | 31 (+V1, V2) | なし | **新規** |

追加が要るブロックは 12 個です。
- 行 v: V1, V2, Sb_v, Sc_v
- 行 tb: T1b, T2b, Sb_tb, Sc_tb
- 行 tc: T1c, T2c, Sb_tc, Sc_tc

位相は全区間で同じで、区間ごとに違うのは「どのブロックが流体か」だけです(plan:87、JM:403-404 と同じ考え方)。

## 3. 断面 (y–z) のブロック構成

格子点 G(j,k) を接続模型(JM:365)から拡張します。
- ys = [YBOT, y_o−DR, y_o, y_c, y_r, b, b+DR, b+top_depth](j = 0..7)
- zs = [0, z_w, z_o, z_o+DR, z_far](k = 0..4)

```
行 j (y)      列0 [0,z_w]  列1 [z_w,z_o]  列2 [z_o,z_o+DR]  列3 [z_o+DR,z_far]
6 tc [b+DR,top]   T1c        T2c          Sb_tc            Sc_tc     ← 新規
5 tb [b,b+DR]     T1b        T2b          Sb_tb            Sc_tb     ← 新規 (機体上面の壁帯)
4 v  [y_r,b]      V1*        V2*          Sb_v             Sc_v      ← 新規 (*x<L_ramp は機体固体)
3 de [y_c,y_r]    N×4        SW*          Sb_de            Sc_de     (*x<L_sw は側壁固体)
2 c  [y_o,y_c]    CW1*       CW2*         Sb_c             Sc_c      (*x<L_cowl はカウル固体)
1 b  [y_o−DR,y_o] U1b        U2b          Sb_b             Sc_b
0 a  [YBOT,y_o−DR] U1a       U2a          Sb_a             Sc_a
```

### 3a. 各ブロックの面

x 面は、−L_up が入口、x_out が出口、それ以外は隣の区間の同じブロックと共有します(端面は §3b)。
表記: 壁(tag)、境界(tag)、共有(相手)。「/」の左が上流側、右が下流側です。

| ブロック | −y | +y | −z | +z |
|---|---|---|---|---|
| U1a / U2a | 境界 bottom | U1b / U2b | sym / U1a | U2a / Sb_a |
| U1b / U2b | U1a / U2a | 壁 cowl_out (x<L_cowl) / CW1・CW2 | sym / U1b | U2b / Sb_b |
| CW1 (x≥L_cowl) | U1b | N_BOT (e8) | sym | CW2 |
| CW2 (x≥L_cowl) | U2b | SW | CW1 | Sb_c |
| N_TOP | (N 内部) | 壁 ramp e1 (x<L_ramp) / **V1** | sym e4 | (N 内部) |
| N_SIDE | — | — | (N 内部) | 壁 sidewall_in e5 (x<L_sw) / SW |
| N_BOT | 壁 cowl_in e8 (x<L_cowl) / CW1 | (N 内部) | sym e9 | (N 内部) |
| N_CORE | — | — | sym e11 | (N 内部) |
| SW (x≥L_sw) | 壁 cowl_in (x<L_cowl) / CW2 | **壁 vehicle (帯) (x<L_ramp) / V2** | N_SIDE | Sb_de |
| Sb_a, Sb_b | bottom / Sb_a | Sb_b / Sb_c | U2a / U2b | Sc_a / Sc_b |
| Sb_c | Sb_b | Sb_de | **壁 cowl_side (x<L_cowl) / CW2** | Sc_c |
| Sb_de | Sb_c | **Sb_v**(模型では vehicle 壁) | 壁 sidewall_out (x<L_sw) / SW | Sc_de |
| Sb_v | Sb_de | Sb_tb | **壁 vehicle_side (x<L_ramp) / V2** | Sc_v |
| Sb_tb / Sb_tc | Sb_v / Sb_tb | Sb_tc / 境界 top_out | T2b / T2c | Sc_tb / Sc_tc |
| Sc_* | 行の連鎖(Sc_a の −y は bottom) | 行の連鎖(Sc_tc の +y は top_out) | Sb_* | 境界 side_far |
| V1 (x≥L_ramp) | N_TOP | T1b | sym | V2 |
| V2 (x≥L_ramp) | SW | T2b | V1 | Sb_v |
| T1b / T2b | 壁 vehicle_top (x<L_ramp) / V1・V2 | T1c / T2c | sym / T1b | T2b / Sb_tb |
| T1c / T2c | T1b / T2b | 境界 top_out | sym / T1c | T2c / Sb_tc |

N の内部共有面(e2, e3, e6, e7, e10)は全区間で共有のままです(plan:314)。

### 3b. 壁が内部面に変わる位置

| x | 壁 → 共有面 | 新しく出る端面 |
|---|---|---|
| L_sw | e5 sidewall_in → N_SIDE\|SW。z_o 面の sidewall_out → SW\|Sb_de。SW の下面と上面は固体どうしだった面が露出して壁になる(cowl_in / vehicle 帯) | SW の断面 = `sidewall_end` |
| L_cowl | e8 cowl_in → N_BOT\|CW1。SW 下面 → SW\|CW2。U1b/U2b 上面の cowl_out → U\|CW。Sb_c 左面の cowl_side → CW2\|Sb_c | CW1+CW2 の断面 = `cowl_base` |
| L_ramp | e1 ramp → N_TOP\|V1。SW 上面の vehicle 帯 → SW\|V2。vehicle_side → V2\|Sb_v。vehicle_top → V\|T | V1+V2 の断面 = `vehicle_base` |
| −L_up | (壁はそのまま入口面まで続く) | 固体 5 ブロックの断面は入口面の「穴」で、面を作らない(旧と同じ) |

### 3c. 節点数の同値類(新規・変更分)

- NFY(行 a): 全区間で最大の帯高から決める
- NV(行 v の y): 両端細分。下端は e5 の端 h_end、上端は h1
- 行 tb: NL
- NT(行 tc)
- NFZ: z_far を 2.5 H にした値で再計算
- 列 0 の z 分布(NZ、z=z_w 側への片側細分)は、行 v/tb/tc にもそのまま伝わる(内部面なので無駄はあるが、位相上避けられない)

## 4. タグと力の帰属

| ID | タグ | 新しい面 | 帳簿 | 旧との対応・差 |
|---|---|---|---|---|
| 1 | inlet_nozzle | N の x=−L_up 断面 | 入口項 | 旧の高さは H−t_c/2(上面 ym+t_c/2、M3:212、M3:254)。新は y∈[0,1] の**面積 H·W/2 ちょうど**で、RN:311 の公称値と一致する |
| 2 | inlet_ext | x=−L_up のその他の流体断面 | — | 同じ。接続模型の名前は `ext_in` |
| 3 | outlet | x_out の全 31 ブロック | — | 同じ |
| 4 | ramp | e1×x、x∈[−L_up,L_ramp]、z≤z_w | ノズル | 同じ範囲。RN:295 の z_split は空振りになる(検査用に残す) |
| 5 | cowl_in | e8×x (x<L_cowl) と SW 下面 (L_sw≤x<L_cowl) | ノズル | 側壁下の板上面 (z∈[z_w,z_o]) が新たに露出する |
| 6 | cowl_out | U1b/U2b の上面、x<L_cowl、z≤z_o | ノズル | 旧は z≤W/2 |
| 7 | bottom | 行 a の −y | — | 模型の名前は `far_bottom` |
| 8 | top_out | 行 tc の +y | — | 同じ(ext_top 版) |
| 9 | sym | z=0 | — | 模型の名前は `symmetry` |
| 10 | side_far | z=2.50 | — | 模型の名前は `far_side` |
| 11 | sidewall_in | e5×x、x<L_sw | ノズル | 同じ |
| 12 | sidewall_out | Sb_de 左面 (z=z_o, y∈[y_c,y_r])、x<L_sw | ノズル | 旧は z=W/2 で y∈[ym,yt]。新では板の側面がここから外れて ID 20 になる |
| 14 | vehicle | **SW 上面の帯 (y=y_r, z∈[z_w,z_o], x∈[L_sw,L_ramp])** | 機体 | 旧は R4e 以降 0 面(RN:89) |
| 15 | vehicle_top | T1b/T2b の −y、x<L_ramp | 機体 | 幅が W/2 → z_o |
| 16 | underside_far | なし | — | W_vehicle が None なので 0 面のまま |
| 17 | vehicle_side | Sb_v の −z、x<L_ramp | 機体 | z = W/2 → z_o |
| 18 | vehicle_base | V1+V2 の −x、x=L_ramp | 機体 | 幅が W/2 → z_o |
| **19** | sidewall_end | SW の −x、x=L_sw | ノズル | 新規。expect_dir (−1,0,0) |
| **20** | cowl_side | Sb_c の −z (z=z_o, y∈[y_o,y_c])、**x∈[−L_up, L_cowl]** | ノズル | 新規。(0,0,−1)。L_sw〜L_cowl に限定しない(plan:783) |
| **21** | cowl_base | CW1+CW2 の −x、x=L_cowl | ノズル | 新規。(−1,0,0)。旧は後縁厚 0 |

**側壁厚みの外側に残るランプ高さの帯の分類**: 帯は x∈[L_sw, L_ramp] の y=y_r、z∈[z_w, z_o] です。

- `ramp` のままにすると不整合が出ます。圧力は RN:300-301 で z>W/2 として機体側に入る一方、摩擦は RN:331 で全面がノズル側に入ります。
- そこで `vehicle` (14) を推奨します。これなら圧力と摩擦の両方が機体帳簿に入ります。
- 法線の向きは既存の `vehicle` の expect_dir (0,1,0)(RN:287)がそのまま使えます。

**B5 への申し送り**:
- RN:285-290 の `spec` に 19〜21 が無いので、今のままでは面が黙って欠落します。追加したうえで、欠落はエラーにする。
- bcond(RN:78-100)にも 19〜21 の壁行が要る。

## 5. 延長則と下側境界

**事実**:
- 旧メッシャ:
  - 後流線 y_mid = y_te + (x−L_cowl) tan θ_b、θ_b = −40.2°(M3:186、RN:148)
  - 下端は水平面 y_bot = −12.32 H(M3:187)
  - 下バンドは nj_bot 55 点で後流線側に集中(M3:256-258)。高さは上流で 12.3 H、x_out で 3 H(4 倍縮む)
  - 上バンドは nj_top 73 点の両側 tanh で、x_out では約 12 H
- 接続模型:
  - YBOT は −1 H 固定、後流線は後縁接線 (−5°)(JM:330)
  - NFY は最大帯高で 1 回だけ決めている(JM:508)ので、帯高が縮むと破綻する(plan:1014-1018)

**下側境界とカウル跡線の延長則の候補**:

| 案 | 内容 | 利点 | 欠点 |
|---|---|---|---|
| **L1(推奨)** | 跡線は後縁接線の直線、厚さ t_c 一定のまま x_out まで。下端は**旧と同じ絶対位置の水平面** y_bot = y_te + (x_out−L_cowl) tan θ_b − bot_depth = −12.32 H | 行 a の帯高は 12.2→11.2 H でほぼ一定なので縮みの破綻が出ない。C¹ で折れが無い。領域が旧と同じになり B6 の旧→新比較ができる。bottom は水平の outflow のまま | 行 a が 12 H と深いので、HFAR 0.10 H のままだと NFY が約 135(推測)。HFAR を緩める(AR = HFAR/Δx_min ≤ 1000 なので HFAR ≤ 約 1 H)必要がある。設計点のせん断層 (−40°) が行 a〜c を斜めに横切り、旧のように格子線に沿わない |
| L2 | 跡線を後縁接線から物理長 L_bl(例 0.5〜1 H、要決定)で θ_b へ滑らかに曲げる(エルミート)。下端は L1 と同じ | 設計点で旧と同じ格子の整列。行 a の縮みは 12.2→3 H に収まる(NFY·f ≲ 3 H を事前判定) | 行 de(N・SW・S)の高さが約 12 H に増え、NY が 105 から約 180 へ(推測)。N 内の両端細分の無駄が大きい。θ_b は設計点 m6_on だけの角度で、m4_off / m10_on では外れる。曲げ区間の薄い CW 帯で skew の危険 |
| L3(非推奨) | YBOT(x) = 跡線 − bot_depth(帯高一定) | 分布が全 x で同じ | 下端が下流に向かって下がるので、自由流が outflow 境界から流入する(BC の種別と合わない)。領域も旧と変わる |

**分割数の決め方(全案共通)**: 各帯の点数は全 station と区間中点の最大長から決めます。そのうえで「点数 × 第一間隔 ≤ 最小長 / 1.2」を事前判定し、満たさなければ生成を止めます(§6.5 の破綻を規則として防ぐ)。後流の壁帯 (U の b、リング) は、壁が無くなった下流で物理長のブレンドにより粗くします(df/dx ≤ 4e-3、plan:370、M3:225-250)。

**ランプ後端・機体ベース・機体上の外部流の接続**:

- **E1(推奨、旧と同型)**:
  - ランプは L_ramp より下流を θ_e の直線で延長する(M3:177。接続模型の B-spline 端接線ではなく、明示した θ_e を使う)。
  - b(x) は M3:435-442 と同じ。下流では y_r,ext + t_base。
  - 行 v は全 z で [y_r, b] とし、x<L_ramp で z≤z_o の部分を固体にする。
  - top_out は b + top_depth に沿うオフセット面(旧と同じ領域)。
  - 行 v の高さは 1.75 H(入口)→ 0.02 H(L_ramp)で **87 倍縮む**。NV は最大高から決め(両端 h1 で約 96 点、推測)、**t_base ≥ NV·h_end** を事前判定する(0.0154 ≤ 0.02 で余裕が小さい)。
- E2(行 v を 2 帯に分ける)は、L_ramp で片方の帯が厚さ 0 に潰れるので不可。
- E3(上端を平面 y_veh+2 にする)は、領域が旧と変わる。

## 6. 未決事項(上位の判断が要る点)

1. **機体の幅**:
   - V1: z_o まで。側壁・カウル・機体の外側面が同一平面になる。推奨。
   - V2: 旧どおり W/2。側壁の厚み分が機体側面から出っ張り、上面が段差の壁になる。列 1 の行 v に流体の細い柱が残る。
   - V3: W_vehicle > W+2t_sw。行 de の上に vehicle 下面が残る。列が増える。
2. 帯の分類を `vehicle` にするか `ramp` にするか(§4)。
3. 下側の延長則: L1 / L2 / L3 と、行 a の HFAR。
4. **x<0 のカウル内壁**: 旧は y=0 で、x=0 で −5° 折れる(M3:186)。接続模型の Contour は x<0 を接線で外挿する(JM:14)ので、旧と食い違う。明示的に平らにし、外壁は折れ点でオフセット線どうしの交点(マイタ)で閉じる、という定義が要る。
5. `vehicle_top`・`vehicle_side`・`vehicle_base` の第一層:
   - 一般壁の h1 にするか、旧生産値(top 6e-4 H、z 法線 4e-3 H、ベース後流 1.7e-4 H、YAML:79,121,122)にするか。
   - vehicle_base の h1e と後流の Δx 上限(t_base/5 = 0.004 H か、1e-3 H か)。
6. t_base 0.02 H と NV の整合(§5)。t_base は暫定モデル値(M3:66-68)。
7. 薄い帯(CW の跡 t_c、SW の跡 t_sw)を x_out まで同じ厚さで延ばすか、下流で厚さを広げるか。
8. ランプフィレット(弧長約 0.026 H)の x 解像。旧は両側 6 区間(M3:137-141)。接続模型の HX 0.05 H では 1 セルに入ってしまう。
9. 入口の排気/外気の節点帰属。旧は節点番号算術(RN:116-122)で、ブロックの所有から作り直す(B4d)。

### 格子規模の事前見積もり(推測、±50 %)

接続模型 B(scale 1)は 434.9 万節点 / 154 station = **断面あたり 2.82 万**です(`notes/investigations/2026-10-06-sern-b4a2/jm_B_report.json`。NL 30、NR 62、NY 85、NZ 50、NSW 19、NTC 19、NFY 20、NFZ 20)。

全体に広げたときの断面あたりの内訳:

| 部分 | 見積もり |
|---|---|
| N | NZ×NY + 2·NZ×NR + NY×NR、NY は全長で 105(secx 実測)。約 1.8 万 |
| 列 0 の他の行 | 50 × (NFY + 30 + 19 + NV 約 90 + 30 + NT 約 35)。約 1.3〜1.7 万 |
| 列 1〜3 | (19+30+NFZ 約 25) × 全行 (約 360〜445)。約 2.7〜3.3 万 |
| **合計** | **約 5.7〜6.8 万 / 断面**。旧生産の 127×39 + 57×39 ≈ 0.72 万の約 8〜9 倍 |

station 数は 全長 296(secx)+ 上流約 15 + プルーム約 50 ≈ 360 なので、**約 2000〜2500 万節点**になります。

- 旧生産 g3 は約 192 万です。
- B4-5 後の変換器は 1579 B/節点で、g5 の 16 GB だと**上限は約 900 万**です(plan:831)。

**規模に効く要素**(旧メッシャに無いもの):

| 要素 | 内容 |
|---|---|
| NZ 50 | 側壁内面を h1 で解像。旧は first_z_frac 4e-3 で nz_in 25 |
| リング NR 62 | 両端細分 |
| e5 系 NY 105 | 両端 h1 |
| 薄い帯 | NSW / NTC 各 19 が全行・全列を貫く |
| 端面 x 解像の波及 | h1e が断面全体に効き、HFAR で遠方が細かくなる |
| 行 v の NV 約 96 | — |
| 行 a の NFY | 深い下側 |

**削り方の候補**(どれも上位の判断。§4.12 の予算 200〜300 万は接続模型に対する決定で、全体には未設定):
- 遠方の HFAR を 0.1 → 0.5 H
- C/D 区間の HX_FAR を 0.08 → 0.2 H(station 数を約 2 割削る)
- 側壁の z 解像を旧並みにする(§6.2 の sidewall_in h1 規格と衝突する)

### Critical Files for Implementation
- /home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py
- /home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py
- /home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py
- /home/sano/work/forge-sern-design/plans/active/tooling-sern-mesh-blocking.md
- /home/sano/work/forge-sern-design/case/46.sern_design/cad/export_contours.py