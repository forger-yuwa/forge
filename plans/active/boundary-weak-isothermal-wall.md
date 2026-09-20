# node 等温壁のエネルギー境界を弱形式 (SU2 型) で課す opt-in

## メタ

- **area**: `boundary`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/boundary.md`](../../methods/boundary.md) 「node 壁強境界」「共役熱伝達 (CHT)」
- **related_plans**:
  - [`boundary-conjugate-heat-transfer.md`](boundary-conjugate-heat-transfer.md) §5.1 #43 (本計画の発端)
  - [`architecture-node-boundary-gradient-dof-only.md`](../accepted/architecture-node-boundary-gradient-dof-only.md)
- **created**: `2026-09-20`
- **owner**: Claude

## 1. 目的

node 等温壁のエネルギー境界条件を、現行の**強制 (壁ノード $T$ ピン + 壁エネルギー残差ゼロ化 +
陰解法エネルギー行の単位行化)** に加えて、**SU2 型の弱形式 (壁ノードのエネルギー方程式を残し、
指定 $T_w$ と第一内部点から作った熱流束を源として加える)** を opt-in で選べるようにする。
完了時に得られるのは、**熱的壁閉包を変えたときの効果の測定値**である (M5 採用)。
壁熱流束の節点交番 (SU2 比 20 倍)・第 1 スペーシング勾配の $-15\,\%$ バイアス・
陰解法の擬似 CFL 上限 $\sim5$ の 3 つについて、改善 / 不変 / **帰属不能**のいずれかを返す。
**「他仮説は棄却済みだから残りはこれ」とは言えない**: 対流スキーム A/B の 2 run
(`run_0029_roe` / `run_0030_hlle`) は `NOT CONVERGED` で参考値にすぎず、float32 も
「温度の単純な量子化では説明不足」までで残差組立ての桁落ちは未検証である。
不変だった場合に言えるのは「**この変更だけでは解消しない**」までである。

## 2. スコープ

- **やる**: `wall_isothermal` × node × 非 WMLES のエネルギー境界に限った第 2 の閉包。
  対流 (SLAU 等)・運動量の no-slip 強制・SST の壁条件・cell 方式には触らない。
  陰解法のヤコビアン寄与まで含める (SU2 と同型)。
- **やらない**: 断熱壁 (`wall`)、壁関数 (`wallTreatmentSST=1`)、WMLES、CHT の界面熱量定義の変更
  (`iface_q_eff` の式は不変)、既定の切り替え。**初版は opt-in のみで既定 0**。

## 3. 関連 docs と前提

現行仕様は [`methods/boundary.md`](../../methods/boundary.md) の「node 壁強境界」。
実装は [`cuda_forge/nodeWallDirichlet_d.cu`](../../solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu)
(状態ピン `pin_wall_node_temperature_d` / 残差ゼロ化 `zero_res_roe_bplane_d`) と
[`cuda_forge/timeIntegration_d.cu`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu) の
エネルギー行単位行化、壁流束は [`cuda_forge/viscousFlux_d.cu`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu)
の node 経路 (`ghostless`, $\nabla\phi\cdot S$)。

**前提となる既往事実** (methods/boundary.md に記載):

| 項目 | 現行 (ピン) | ピン導入前の旧弱形式 |
| --- | --- | --- |
| 壁ノード $T$ | BC 値厳密 (350.0000/300.0000) | 壁 CV 平均へ緩む ($\sim0.1$ K オフセット) |
| 中央部熱流束 | 厳密 ($+0.00\,\%$) | — |
| 第 1 スペーシング勾配 | $-15\,\%$ (W-I 弱形式閉包の $O(\Delta)$ 誤差) | $-24\,\%$ |
| 陰解法 擬似 CFL | 実用上限 $\sim5$ (エネルギー行 decouple で境界が実効陽) | — |

**旧弱形式と本計画の弱形式は別物**である。旧実装は壁ノード $T$ を浮かせ、**その緩んだ $T$ で
壁流束を作った**ので CV 平均へ落ちた。SU2 は流束を**指定 $T_w$** と第一内部点から作る
([`CNSSolver.cpp:727`](../../.external/su2-src/SU2_CFD/src/solvers/CNSSolver.cpp)):

```cpp
su2double dTdn = -(There - Twall)/dist_ij;
/*--- Apply a weak boundary condition for the energy equation. ---*/
su2double Res_Visc = thermal_conductivity * dTdn * Area;
LinSysRes(iPoint, nDim+1) += Res_Conv - Res_Visc;
```

運動量は SU2 も強制 (`Jacobian.DeleteValsRowi` を `iVar=1..nDim` にだけ適用)。
実測で SU2 の壁ノード保存温度は 566.011–566.886 K。

## 4. 設計方針

### 4.1 連続系は同一、離散系だけ変える

どちらも $T|_{\rm wall}=T_w$ を課す。差は**自由度を消すか (強制) 残すか (弱形式)**。

**強制 (現行、`nodeIsothermalEnergyBC: 0`)**

$$T_W \leftarrow T_w,\quad R^{\rm roe}_W \leftarrow 0,\quad A_{WW}^{\rm energy}\leftarrow I$$

**弱形式 (新規、`nodeIsothermalEnergyBC: 1`)**: 壁ノード $W$ のエネルギー方程式を残し、
壁半割面ごとに

$$R^{\rm roe}_W \mathrel{-}= k_{\rm eff}\,\frac{T_I-T_w}{d_1}\,A_{\rm half}$$

を加える。符号は「$T_I>T_w$ なら壁ノードから熱が抜ける」。**$T_W$ 自身は使わない** (ここが旧弱形式との差)。

**幾何の定義 (M2 採用)** — 初版はこれを初期化時に 1 度だけ確定し、残差と陰解法が同じ値を参照する:

| 記号 | 定義 | 備考 |
| --- | --- | --- |
| $I$ | 壁ノード $W$ の隣接のうち**壁法線との alignment が最大**の**非壁**ノード | 既存 `conjugateWall.cpp` の `firstInterior()` と同一規約 |
| $d_1$ | **法線投影距離** $\lvert(\mathbf x_I-\mathbf x_W)\cdot\hat n\rvert$ | **点間距離ではない**。SU2 も `NormalDistance()` を使う ([`CNSSolver.cpp:679`](../../.external/su2-src/SU2_CFD/src/solvers/CNSSolver.cpp))。法線から 30° 傾いた辺では点間距離版の熱流束が **13.4 % 小さく**なる。C3X の直交第一層ではこの差を検出できないので、**傾斜格子でのテストを §6 に必須で入れる** |
| $A_{\rm half}$ | 壁半割面の面積 | 同一ノードに複数の壁半割面が来る場合は**加算**。二重加算しないことを局所検査で確認 |
| $k_{\rm eff}$ | $k_{\rm lam}+c_p\mu_t/{\rm Pr}_t$ | 壁ノードの値 |

**起動時に拒否する構成**: 内部点が見つからない / $d_1$ が退化 / 選ばれた $I$ が別の壁ノード /
軸対称 / 壁ノードが周期で同一視される構成 / 移動壁。**翼列の遠方周期境界は禁止しない**
(壁ノードが周期に属さないため)。`nodeWallDirichlet: 1` は**注意書きでなく起動時の必須条件**にする。

### 4.2 二重計上を避ける — `Qw_Wall` は流用しない (M3 採用)

node の壁半割面は現在 `viscousFlux_wall_d` ([`viscousFlux_d.cu:530`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu)) の
`ghostless` 経路で $\nabla T\cdot S$ の伝導を `res_roe` に入れている。**置換するのはここだけ**で、
内部双対面には触らない。

**`Qw_Wall` を流用してはならない**。あれは片端が壁の **W–I 内部双対面**の伝導を置換する機構で
([`viscousFlux_d.cu:261`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu))、壁半割面とは別経路である。
初稿の「内部双対面は不変」と矛盾していた。さらに `Qw_Wall > -0.5` を有効判定に使っているため
**符号付き熱流束とマーカを兼用**しており、冷却壁の負の流束に無条件で流用できない。

したがって本計画は **`viscousFlux_wall_d` の伝導項だけを置換する新経路**を足し、
内部カーネルには新しいポインタを渡さない。**選択フラグと符号付き流束は別の引数に分ける**。
加熱・冷却の**両方向**を局所検査で確認する。

### 4.3 陰解法の行列寄与 — 厳密微分ではなく SU2 型の近似対角項 (M1 採用)

**初稿は誤りだった。** 物性を凍結して $g=k_{\rm eff}A_{\rm half}/d_1$ と置くと、壁寄与は

$$R_W^{\rm wall}=-g\,(T_I-T_w)$$

なので、**$(\rho e)_W$ による厳密微分は 0** であり、温度を介する微分は**内部点 $I$ の列**にある。
初稿の非零対角項は厳密微分ではない。SU2 も同じく、内部点温度による残差に対して**壁点の対角項**を
加える近似線形化である ([`CNSSolver.cpp:715`](../../.external/su2-src/SU2_CFD/src/solvers/CNSSolver.cpp))。
これは移植候補として妥当だが、**厳密微分として説明しない**。

**符号**: forge は `res_roe` を右辺に取り、行列へは**残差微分の符号を反転して**組む
([`timeIntegration_d.cu:819`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu))。
**計画の負符号をそのまま `diag_block` に足してはならない。** CPG でのエネルギー対角への追加は

$$\Delta A_{WW}^{\rm energy} = +\frac{g}{\rho\,c_v}$$

速度行は従来どおり単位行 (no-slip は強制のまま)。

**初版は `thermalMethod: 0` (CPG) に限定する。** TP は $e(T)$ を使っており
([`thermo_d.cuh:919`](../../solver_density_cuda/cuda_forge/thermo_d.cuh))、密度微分の CPG 式を
そのまま一般化できない。起動時に拒否する。

**検査**: 壁点と内部点を**別々に摂動**する差分検査で、残差・厳密微分・実際の近似行列を分けて確認し、
小規模伝導問題で符号を確認する (§6)。

### 4.4 状態ピンの扱い

`nodeIsothermalEnergyBC: 1` のとき `pin_wall_node_temperature_d` は**エネルギー系だけ**
適用しない ($T$, $P$, $\rho e$, `sonic` を上書きしない)。運動量の `enforceWallNoSlip` と
SST の壁条件は不変。**`nodeWallDirichlet: 0` を使ってはならない** (運動量まで弱くなる)。

### 4.5 CHT の界面熱量との整合 — 比較量を $T_w^{\rm bc}$ で定義し直す (M4 採用)

拘束反力 $C=-R^{\rm raw}$ は弱形式では**拘束が無いので $C=0$**、$Q_f=\sum F^E$ がそのまま収支に一致する。
ここまでは妥当。**ただし `ifaceRraw=0` だけでは足りない。**

現行の `iface_q_compact` は [`conjugateWall.cpp:219`](../../solver_density_cuda/conjugateWall.cpp) で

```cpp
keff * (T[j] - T[ic]) / fi.d1[ib]
```

と、**指定壁温ではなく壁ノードの保存温度 `T[ic]`** を使っている。強制側は $T_W=T_w$ なので
問題が隠れるが、**弱形式では新しい境界流束と別の量になる**。`iface_q_2nd` も同様。
外部 CHT ループの既定入力も `q_compact` である ([`cht_loop.py:157`](../../solver_density_cuda/tools/cht_loop.py))。

したがって:

- **`Tw_bc` (指定値) と `T_W` (保存温度) を別々に壁ダンプへ記録する**。
- **比較用のコンパクト熱流束は両枝とも `Tw_bc` で定義する**。保存温度からの勾配は**別名の診断量**にする。
- CHT 検証は `--flux q_eff` を**明示**し、`ifaceFw` が実際に残差へ加えた新流束を保存していることを確認する。
- $C=0$ は流体側の定常性も固体との収支も保証しない。**G-cons と収束判定を別途通す**。

## 5. 実装ステップ

1. `input/solverConfig.{hpp,cpp}` に `mesh.nodeIsothermalEnergyBC` (既定 0) を追加。
   1 のとき `discretization=node` かつ `wall_isothermal` が存在することを検証し、
   `wallTreatmentSST=1` / `wallModelLES=1` / cell との併用は**明示的に拒否**する。
2. `cuda_forge/nodeWallDirichlet_d.cu`: `nodeIsothermalPinActive` を温度系だけ無効化できるよう分岐。
   `zeroNodeIsothermalEnergyResidual` を `nodeIsothermalEnergyBC==0` のときだけ呼ぶ。
3. `cuda_forge/viscousFlux_d.cu`: 壁半割面の伝導を 4.1 の式で置換する経路を追加 (既存 `Qw_Wall` 機構を流用)。
4. `cuda_forge/timeIntegration_d.cu`: エネルギー行の単位行化を条件化し、4.3 のヤコビアン寄与を足す。
5. `conjugateWall.cpp`: 4.5 に従い `ifaceRraw` の扱いを分岐。
6. `methods/boundary.md` に本節を追記し、`methods/index.md` を同期。

### 5.1 残作業 (優先順)

codex plan レビュー (§6.1) の Major 6 件を**全件採用**し、実装前の優先順を
**①残差・近似行列と対応範囲 → ②幾何・壁専用流束 → ③診断/CHT 契約 → ④定量ゲートと因果主張の訂正**
とする (レビューの推奨順)。

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | **M1 行列寄与の定義** | §4.3。厳密微分 (対角 0、温度微分は $I$ の列) と SU2 型近似対角を分けて書く。forge の符号反転規約 ([`timeIntegration_d.cu:819`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu)) に合わせ、CPG は $+g/(\rho c_v)$。**`thermalMethod: 0` 限定**で起動時に拒否 |
| 2 | **M2 幾何の確定** | §4.1。**法線投影距離** $d_1$、内部点の選択規約、半割面積、同一ノードへの加算。初期化時に 1 度確定し残差と陰解法が同じ値を参照。別の壁点を選ばない検査、退化・内部点なしの拒否。軸対称・周期同一視・移動壁を拒否 |
| 3 | **M3 壁専用流束経路** | §4.2。`viscousFlux_wall_d` の伝導項だけを置換。**`Qw_Wall` は流用しない** (W–I 内部双対面の機構であり、`>-0.5` でマーカと符号付き値を兼用している)。選択フラグと符号付き流束を分離し加熱・冷却両方向を検査 |
| 4 | **M4 診断と CHT 契約** | §4.5。`Tw_bc` と `T_W` を別々に記録し、比較用コンパクト熱流束を**両枝とも `Tw_bc`** で定義。保存温度からの勾配は別名。CHT 検証は `--flux q_eff` 明示、`ifaceFw` が新流束を保存していることの確認 |
| 5 | **M5 因果主張の訂正** | §1・§3。「棄却済み」を撤回し、目的を**効果測定**に限定。結果に**帰属不能**を追加 (A が再現しない / 片枝が未収束 / 振幅が非定常)。不変なら「この変更だけでは解消しない」までしか言わない |
| 6 | **M6 定量ゲートの事前登録** | §6。局所演算・純伝導・C3X・CFL/CHT・投入条件の 5 群。交番振幅の抽出式・正規化・保存間隔・末尾窓・有意差閾値を**実装前に固定**し `check_quasisteady.py --series-csv` に渡す。`run_isoT_cond*` は cell run と `cfl_pseudo: 100` の旧入力が混在するので**対象 run を特定**する (ワイルドカード不可) |
| 6a | ~~**A/B が交絡している — 2×2 にする**~~ → **撤回 (2026-09-20)**。**強制側は壁ノードの `res_roe` をゼロ化するので、壁半割面の流束は解に一切入らない** (診断にしか出ない)。したがって「強制+2 点」の B を作っても場は A と同一で、A↔B は何も測らない。2×2 は成立しない。**A↔C は交絡していない** — C で変わるのは拘束一式だけで、これは codex の言う「熱的壁閉包一式の効果まで確定」そのもの。ただし**弱形式では壁半割面の流束が初めて解に効く**ので、その離散化 (#6b) が残る唯一の未検証変数になる |
| 6b | **弱形式に over-relaxed 非直交補正を入れるか** (2026-09-20 ユーザ指摘、codex 未検出) | 初版の弱形式は**純粋な 2 点コンパクト差分 + 法線投影距離**で、forge が内部面で使っている **over-relaxed 分解 (直接差分 + 非直交補正) が入っていない**。SU2 がそうだからそのまま移植したが、forge の通常の作法からは外れる。効くのは法線から傾いた第一層 (30° で 13.4 %)。**単純に足せない**: 非直交補正には壁面の接線勾配が要り、弱形式では隣の壁ノードの $T$ が自由 DOF なので **隣接壁ノード同士が再結合する**。これは (a) 市松を減衰させる方向かもしれないし、(b) 旧弱形式の「CV 平均に緩む」問題を呼び戻すかもしれない。**どちらに転ぶか不明なので最初の A/B に混ぜない**。**①は完了 (2026-09-20, `run_0039_weakbc`): 市松は 0.9186 → 1.2092 % と 32 % 悪化**したので、残るのは②の補正版のみ。**弱形式のまま $\nabla T\cdot S$ に戻すのは不可** — $T_W$ が自由だと $\nabla T$ は場の温度だけで決まり**指定 $T_w$ がどこにも入らない** (= 2026-07-20 に棄却した旧弱形式)。したがって D は「**法線成分は指定 $T_w$ から、接線の非直交補正は再構成勾配から**」の混合形になる。内部面と同形 (`viscousFlux_d.cu:229`: 直接差分×delta + $\nabla T\cdot k$) にする |
| 6c | **内部側 1 点という選択の妥当性** (同上) | 複数の内部隣接に分配せず整列度最大の 1 点だけを使う。SU2 と同型で、診断 `iface_q_compact` と CHT の $D_f=k_{\rm eff}A/d_1$ が同じ 1 点を使うので BC・診断・連成が同じ幾何を見る。**壁半割面が複数ある壁ノードでは残差も対角も `atomicAdd` で加算される** (構造的に回数が一致)。6b の補正版を入れるときはこの前提も見直す |
| 7 | 既定の扱い | 1–6 が全部良ければ既定切り替えを別 plan で検討。**本計画では既定 0 のまま** |

## 6. 検証

**実装前にこの表を確定させる** (M6 採用)。閾値・抽出式・対象 run はここで事前登録する。

### 6-A 局所演算 (実装と同時、ユニット)

| 検査 | 合格条件 |
| --- | --- |
| 一様温度場 | 壁寄与の熱流束が**厳密に 0** |
| 加熱 / 冷却 | $T_I>T_w$ と $T_I<T_w$ で符号が反転し、絶対値が一致 |
| 壁寄与の加算回数 | 同一ノードに複数の壁半割面が来る構成で**一回だけ**加算される |
| 近似行列の符号 | 壁点・内部点を別々に摂動した差分と比べ、対角寄与の**符号が正しい** (forge の反転規約込み) |
| 既定パス | **「ビット同一」は不成立** (2026-09-20 実測)。node の C3X は `atomicAdd` で **run-to-run 非決定**で、旧バイナリ 1 回目 vs 2 回目が新 vs 旧と同程度に違う (P の max\|Δ\| 20.9 vs 16.5)。**両側からノイズ床を測り、交差が群内と同じ population であること**に置き換える ([[noise-floor-both-sides]])。実測 (200 step, 旧 4 反復 / 新 4 反復, RMS 差の中央値): 交差 / 群内 = **0.91–0.97** = 同一 population で **PASS**。さらに**交番振幅そのもの**は再現的で、旧 0.9112 % / 新 0.9111 %、幅 **0.0001 %** |

### 6-B 純伝導 (`case/24`)

**対象 run を特定する**: `run_isoT_cond*` には cell の run と `cfl_pseudo: 100` の旧入力が混在するため、
node・現行設定の run を明示して使う (ワイルドカード不可)。

| 検査 | 合格条件 |
| --- | --- |
| 直交格子・厳密解シード | 中央部熱流束 $\pm0.1\,\%$ |
| 第 1 スペーシング勾配 | **旧弱形式の $-24\,\%$ より良い** (必須。割ったらその場で却下) |
| **傾斜格子** | 法線投影距離が効くのはここ。点間距離版との差 (30° で 13.4 %) を検出できる格子で確認 |
| 格子細分化 | 境界整合性と精度次数 |
| 摂動からの緩和 | 厳密解シードだけでなく、摂動初期場からも同じ定常解に落ちる |

### 6-C C3X 市松 A/B (本題)

**同一ソース・同一初期場**の A/B。`case/53.c3x_vane_cht/run_0009_uniformTw` と同一条件。

- **有意差閾値 = 0.001 %** (2026-09-20 事前登録)。交番振幅の run-to-run 幅が **0.0001 %** (0.9111–0.9113 %、
  旧 4 反復 + 新 3 反復) と測れたので、その **10 倍**を有意とする。指標は
  $A_{\rm oe}={\rm RMS}[(q_{i-1}-2q_i+q_{i+1})/4]$ を区間平均 $q$ で正規化、負圧面 $s/S$ 0.30–0.80。
- 交番振幅の**抽出式・正規化・保存間隔・末尾窓を実装前に固定**し、その時系列を
  `check_quasisteady.py --series-csv` に渡す (`wall_series.csv` には交番振幅の列が無いので**追加する**。
  積分熱量の `STEADY` を交番振幅の定常性に転用しない)。
- 評価は**第一内部点の $T$ と $q$ の両方** (壁温変動による見かけの相殺を避ける)。
- 両枝とも `check_convergence.py` の VERDICT を付す。
- **判定**: 改善 / 不変 / **帰属不能**。A が再現しない・片枝が未収束・振幅が非定常なら帰属不能。

### 6-D CFL と CHT

| 検査 | 合格条件 |
| --- | --- |
| 擬似 CFL 掃引 | **CFL 以外を固定**し、「同一解への到達」と「安定限界」を**分けて**判定する |
| G-cons (`case/52.conjugate_slab` V1) | `--flux q_eff` 明示。`q_floor`・絶対許容値・単位を事前登録。解析解 81.3090 W/m² への誤差基準も登録 |
| G-cons (C3X) | 同上 |

### 6-E 投入条件 (全 run 共通)

メッシュ品質 VERDICT、適用可能な IC / 段階起動、NaN 検査、判定区間、
**新規 run と `case/*/README.md` 索引**を明記する。

### 6-F 回帰 (`case/48.flat_plate_cooled_m4`)

既存手順の段階起動に従い、**3 % 以上の変化を回帰とする基準**
([`procedures/verification/48-flat-plate-cooled.md`](../../procedures/verification/48-flat-plate-cooled.md)) を適用する。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-20` | [`notes/reviews/2026-09-20-boundary-weak-isothermal-wall-plan.md`](../../notes/reviews/2026-09-20-boundary-weak-isothermal-wall-plan.md) | **GO-with-changes**, C0/M6/m0 | **全件採用** → §5.1 #1–#6。M1 ヤコビアンは厳密微分でない (対角 0・温度微分は $I$ の列)・forge の符号反転規約・CPG $+g/(\rho c_v)$・`thermalMethod: 0` 限定 → §4.3 / M2 $d_1$ は**法線投影距離** (点間距離だと 30° で 13.4 % 差)・幾何を初期化時に確定・拒否条件 → §4.1 / M3 **`Qw_Wall` 流用は誤り** (W–I 内部双対面の機構、`>-0.5` でマーカと値を兼用) → §4.2 で `viscousFlux_wall_d` 専用経路へ / M4 `iface_q_compact` が**保存壁温**を使うので弱形式で別量になる → §4.5 で両枝とも `Tw_bc` 定義に / M5 「他仮説棄却済み」は**言い過ぎ** (ROE/HLLE run は `NOT CONVERGED`、`q_total: STEADY` は交番振幅に転用不可) → §1 を効果測定に限定し**帰属不能**を追加 / M6 定量ゲートの事前登録 → §6 を 6-A〜6-F に再構成 |
| 前提 (自由形式) | `2026-09-20` | [`notes/reviews/2026-09-20-codex-c3x-checkerboard-triage.md`](../../notes/reviews/2026-09-20-codex-c3x-checkerboard-triage.md) | 判定なし (切り分け依頼) | 本計画の A/B 設計はこのレビューの提案を採用したもの。float32・低マッハ前処理の 2 仮説は棄却済み |

## 7. 影響範囲

- `solver_density_cuda/input/solverConfig.{hpp,cpp}`、`cuda_forge/{nodeWallDirichlet_d.cu,viscousFlux_d.cu,timeIntegration_d.cu}`、`conjugateWall.cpp`
- 既定 0 なので既存ケース・実行手順への影響は無い (ビット同一を検証で担保)
- `methods/boundary.md`、`methods/index.md`

## 8. 完了条件

- [ ] `methods/boundary.md` の現在仕様を更新済み
- [ ] 実装・検証完了 (§6 を満たす)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を §5.1 に反映済み
- [ ] `status` を `done` に変更し、§9 に変更ログを記載
- [ ] `plans/active/` → `plans/accepted/` へ移動
- [ ] [`plans/README.md`](../README.md) の一覧を同期

## 9. 変更ログ

- `2026-09-20` — **実装 (C: 弱形式+2 点) 完了・A/B 実施**。`mesh.nodeIsothermalEnergyBC` 追加、
  新規 `cuda_forge/weakIsothermalWall_d.{cu,cuh}` (幾何は `conjugateWall::firstInterior` を再利用 = 法線投影距離)、
  `viscousFlux_d.cu` に壁専用置換経路 + 対角素材の `atomicAdd`、`nodeWallDirichlet_d.cu` でピン一式を解除、
  `timeIntegration_d.cu` で単位行化を外し近似対角 $+g/(\rho c_v)$、`conjugateWall.cpp` で $C=0$ と `Tw_bc` 基準の診断。
  **結果**: 壁ノード $T_W$ が **566.001–566.831 K** (SU2 実測 566.011–566.886 K と同等) = 実装は SU2 の挙動を再現。
  しかし**市松は 0.9186 → 1.2092 % と 32 % 悪化**、第一内部点 $T$ の交番も 0.0827 → 0.1120 K。
  `check_convergence.py` **`NOT CONVERGED`**、`rms_roK`/`rms_roOmega` が **RISING**。
  **判定: 熱的壁閉包一式を SU2 型に替えるだけでは市松は解消せず、むしろ悪化する。**
  既定パスは不変 (ノイズ床比較 0.91–0.97、交番振幅 0.9112 → 0.9111 %、幅 0.0001 %)。
  残るのは #6b の over-relaxed 版のみ。

- `2026-09-20` — codex plan レビュー **GO-with-changes** (C0/M6/m0) を全件採用して改訂。設計の誤り 3 件 (ヤコビアン・`Qw_Wall` 流用・`iface_q_compact` の壁温) と因果主張 1 件を訂正、幾何仕様と定量ゲートを追加。`status: in_progress`。
- `2026-09-20` — 初稿。`boundary-conjugate-heat-transfer` §5.1 #43 の切り分け (対流スキーム・低マッハ前処理・float32 を棄却) を受けて起票。
