# node 周期の継ぎ目勾配の修正 (LSQ の 2 重計上、SST k/ω 勾配の未合算)

## メタ

- **area**: `boundary`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/gradient.md`](../../methods/gradient.md) (「境界寄与」node × 周期の継ぎ目)
  - [`methods/discretization.md`](../../methods/discretization.md) §7.3 (LSQ)、§4.5 (周期の DOF 同一視)
- **related_plans**:
  - [`boundary-node-rotational-periodic.md`](boundary-node-rotational-periodic.md) (回転周期。本 plan の完了を前提にする。本 plan は同 plan の §4.0 / §5.1 #0 を切り出したもの)
  - [`discretization-median-dual-3d.md`](discretization-median-dual-3d.md) §4.5 (node 周期の骨格)
- **created**: `2026-09-26`
- **owner**: `CFD Dev`

## 1. 目的

node の周期境界の継ぎ目にある既存の欠陥 2 件を直す。どちらも回転と無関係に、**今の並進周期の node 計算 (非軸対称) すべてに効く**。本 plan の適用範囲は **node ∧ 並進周期 ∧ 非軸対称** (codex plan M2: 軸対称では gather が `periodicNode_d.cu:174` で return し、GG は `A_planar` を使うので前提が異なる)。

1. **LSQ 勾配の 2 重計上**: node の勾配は `gradLSQ: 2` 固定で、継ぎ目で割れた各部分 CV が片側の隣接だけで完全な LSQ 勾配を解き、
   `periodicGradientGather` が和を取る。**線形場で継ぎ目の勾配が正確に 2 倍** (2026-09-26 実測、case/09 TGV 32³、`Ux = 10+y` で x 継ぎ目 1458 点の
   $\partial U_x/\partial y$ = 2.000000、内部 1.000000。`case/09.Taylor-Green/_g0_lsq_seam/`)。
2. **SST の $k,\omega$ 勾配が合算されていない**: `periodicGradientGather` (`main.cpp:1459`) の後で `ransGradient` (`main.cpp:1481`) がゼロから作り直し
   (`ransTransport_d.cu:211`)、その後に gather が無い。F1 (`ransSource_d.cu:344`) は片側の勾配を読む。`ransGradient` は周期半割面も除外していない。

**スカラー全体の点検 (2026-09-26、ユーザ質問「その他スカラーは大丈夫?」)** — main で勾配を作る 4 か所 (`calcGradient`・`speciesGradient`・
`passiveGradient`・`ransGradient`) と継ぎ目の gather の順序を確認した:

| 量 | 勾配 | 継ぎ目 | 判定 |
| --- | --- | --- | --- |
| 速度・密度・圧力・温度 | LSQ | gather は和 | **欠陥 1 (2 倍)** |
| $k,\omega$ | Green–Gauss (`ransGradient`) | gather の**後**に作り直され、周期半割面も除外していない | **欠陥 2 (片側)** |
| 化学種 | Green–Gauss (`speciesFaceReconstruction ≥ 1` のみ) | 周期半割面を除外 (`excludePeriodic`) し、合併体積で割った部分寄与の和 | 正しい |
| 受動種・凝縮モーメント | Green–Gauss (`passiveGradient`) | 同上 | 正しい |
| $\gamma,\ Re_{\theta t}$ | 拡散は 2 点差分 (勾配配列なし)。生成項は速度勾配を読む | 状態は root からミラー | 拡散は影響なし。**生成項は欠陥 1 の影響を受け、欠陥 1 の修正で直る** |

**ユーザ決定 (2026-09-26): スカラーの勾配も LSQ に揃える**。node の勾配は NS の原始変数だけ LSQ (`calcGradient`、`gradLSQ: 2` 固定) で、
$k,\omega$ (`ransGradient`)・化学種 (`species_gradient_d`)・受動種・凝縮モーメント (`passiveGradient`) は Green–Gauss。なお既存の「フォールバック」は
LSQ の退化方向を 0 にするスペクトル打ち切りで、Green–Gauss への差し替えではない (`methods/discretization.md` §7.3.1)。
**進め方**: 本 plan では合併 stencil の LSQ 係数を**どの変数にも使える形**で作る。スカラー勾配の LSQ 統一は、全 node 計算のスカラー勾配を
(継ぎ目以外も) 変える変更なので、本 plan の完了後に別 plan で行う ($k,\omega$ → 化学種 → 受動種・凝縮の順、冷却平板・SERN・凝縮ノズルで回帰。
壁近くの $\omega\sim1/y^2$ と化学種の有界性に注意)。統一すれば欠陥 2 と同じ型 (Green–Gauss の継ぎ目の別扱い) は構造的に無くなる。

## 2. スコープ

- **やる**: 合併 stencil の LSQ 係数の事前計算、$k,\omega$ 勾配の合算位置の修正、並進周期での検証 (G0 拡張、二次場、case/39、case/09)。
- **やらない**: 回転周期 (`boundary-node-rotational-periodic` で本 plan の後に)。Green–Gauss の勾配 (化学種・受動種・凝縮モーメント。現行の和が正しい)。
  陰解法の行縮約 (別 plan)。継ぎ目の部分双対面の合併 (並進の押し出しでは部分面が同一平面で差が無い)。

## 3. 関連 docs と前提

- LSQ の事前計算: `calcGradient_d.cu:577-920` (重み $w=1/|\Delta\mathbf x|^2$ は :668、スペクトル打ち切りは :671)。
- 継ぎ目の gather: `periodicNode_d.cu:166-212`。
- 双対面は primal edge を一意化して生成し、**継ぎ目の接線方向のエッジは両側に存在する** (`gmshReader.hpp:1841`)。

## 4. 設計方針 (2026-09-26 `diagnostician`、codex plan (回転 plan の 2 回目) を採用)

### 4.1 合併 stencil の LSQ (codex plan M1 で精度仕様を固定)

- **適用条件**: node ∧ 並進周期 ∧ 非軸対称。係数の合併と §4.2 の gather は**同じ条件関数**を使う。軸対称の既存経路は変えない。
- **周期像の識別と幾何の照合を分ける**: 同一物理隣接の候補は、隣接の `periodicRoot` が一致すること。そのうえで幾何を
  $\lvert\Delta\mathbf x_{mj}-\Delta\mathbf x_{m'j'}\rvert \le 10^{-4}\,h_{min}$ ($h_{min}$ = その節点の最短内部エッジ長。float32 座標差の丸め ~1e-7 相対より十分大きく、
  隣接間隔より十分小さい) で照合する。**`periodicRoot` が同じでも $\Delta\mathbf x$ が違う (異なる周期像) ものは統合しない**。
- ~~同値類 $E$ ごとに $d_E, w_E$ を最初の incidence の値で 1 つ決める~~ → **決着 (2026-09-26、実装レビュー M1)**: 行列・係数とも**各 incidence の実変位 $d_{mj}$** ($w_{mj}=1/\lvert d_{mj}\rvert^2$) を使い、
  同値類は配分係数 $\alpha$ の決定にだけ使う。代表値で置き換えると、実行時に各部分 CV が読む $\phi_j-\phi_m$ と係数の変位が float32 丸め分ずれ、
  線形場で誤差が出る (codex 反例: $d_0=f32(0.03)-f32(0)$、$d_1=f32(100.03)-f32(100)$)。
- 配分係数 $\alpha = 1/\text{重複数}$。$M_r=\sum_{m,j}\alpha_{mj} w_{mj} d_{mj} d_{mj}^{\mathsf T}$。**スペクトル打ち切りは $M_r$ に 1 回だけ**。
- 各部分 CV の incidence の係数 $c_{mj}=M_{r,\tau}^{+}\,\alpha_{mj}w_{mj} d_{mj}$ を焼き込む。**毎 step の gather は現行の和のまま**。
- 壁∩継ぎ目では実在する内部隣接だけを合併する (壁の疑似点は足さない)。
- 係数は**変数に依らない**形で作る (後続のスカラー LSQ 統一で流用)。
- **不採用**: root の値での上書き、一律 0.5 倍 (4・8 member の角、非対称 stencil、部分的な rank 欠損を扱えない)。
- codex の CPU 数値確認: 重複を持つ非対称 stencil を 2・4・8 member に配分した合併は、一意 stencil の解と最大 3.4e-16 で一致 (式の確認)。

### 4.2 SST の $k,\omega$ 勾配 (本 plan では Green–Gauss のまま合算を直す)

- `ransGradient` の直後 (`ransBlendF1` の前) に $k,\omega$ 専用の gather (Green–Gauss、**周期半割面を除外して積算**、合併体積で割った部分寄与の和)。
  早い方の gather (`main.cpp:1459`) からは $k,\omega$ を外す。順序は `ransGradient → 専用 gather → ransBlendF1 → ransTransport`。
- **LSQ 化は後続 plan** (`diagnostician`: 本 plan で LSQ 化まで入れると、R1 の変化が「継ぎ目の欠陥の修正」と「スキームの変更」の 2 要因になり帰属できない。
  合併 GG の $k,\omega$ は後続の LSQ 化の回帰参照にもなる)。
- **F1 の初回上書き (codex m5)**: `buildScalarDescs` (`ransTransport_d.cu:104`) が初回に `sstF1` を 1 で埋め、直前に計算した F1 を使わない。初期充填を変数初期化時へ移す (**実装済 2026-09-26**: `variables.cpp` の `allocVariables` で 1 を入れ、`buildScalarDescs` は副作用なし)。
  初回 step が変わるので、非周期 run も 1 step 目からビット差が出る (回帰は 2 step 目以降のノイズ床比較で判定)。

### 4.3 影響

継ぎ目の節点の勾配が変わる (意図した修正)。粘性応力・2 次再構成・リミタ・SST 生成・F1 に効く。並進の他の演算経路は不変。

## 5. 実装ステップ

1. codex plan 段。
2. G0 の拡張 (試験を先に: 2/4/8 member、非対称 stencil、壁∩継ぎ目、root 交換) と再現物の保存。
3. 合併 LSQ の実装 (`calcGradient_d.cu`、`mesh.cpp`)。
4. $k,\omega$ 勾配の gather (`main.cpp`、`periodicNode_d.cu`、`ransTransport_d.cu`)。
5. 二次場・case/39・case/09 の検証。
6. docs、codex result。

### 5.1 残作業 (優先順)

**計算資源**: 小規模はローカル、重いものは AWS (ユーザ指示: ローカルで大きな計算をかけない。AWS は他セッションと共有、手動 stop しない)。
**ユーザ指示: 実装前・区切りごとに codex に諮る**。

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | codex plan 段 | Critical/Major は `diagnostician` に諮る | F |
| 2 | G0 拡張と再現物 | 線形場を生成スクリプトで焼く (`initial` ではなく)。入力 h5 生成・集計スクリプト・修正前の勾配出力・revision を `case/09.Taylor-Green/_g0_lsq_seam/` に保存し README に登録。2/4/8 member (面・辺・角)、非対称 stencil、壁∩継ぎ目、root 交換のメッシュ | O |
| 3 (**実装済 2026-09-26**、正式試験は #5) | 合併 LSQ | `calcGradient_d.cu` の事前計算、`mesh.cpp` の group 情報。合格: §6 G0 | O (F レビュー、条件 6) |
| 4 (**実装済 2026-09-26**、正式試験は #5) | $k,\omega$ 勾配の gather | `main.cpp`、`periodicNode_d.cu`、`ransTransport_d.cu`。合格: §6 G1 | O (F レビュー) |
| 3a (**実装済 2026-09-26**) | 実装レビュー M1・m3 | M1: 各 incidence の実変位で組む (§4.1)。m3: F1 初期値を `allocVariables` へ (§4.2) | O (判断: 2026-09-26 `diagnostician`・全件採用) |
| 5a | G0/G1/G2/G2′ ハーネス | #2 の再現物に加え、float32 反例 (期待 1.000000、root 両順序)、ジッタ格子 (±0.2h、決定論的、周期像は同じ量)、CPU double 参照、F1≠1 入力で 1・2 回目に輸送が読む値。合格: §6 G0/G1/G2/G2′ | O |
| 5b | R3 | case/48・case/16・case/44 の勾配配列が旧新ビット同一 (F1 初期化で 1 step 目が変わる run は 2 step 目以降のノイズ床比較) | O |
| 5 | 検証 R1・R2 | §6 R1・R2 (5a・5b の後)。**区切りで codex** | O (結論 F) |
| 6 | docs + codex result | `methods/gradient.md` の「修正中」を外す | F |

## 6. 検証 (測る前に固定)

| # | 試験 | 合格 (測る前に固定) |
| --- | --- | --- |
| G0 | LSQ の局所作用素試験: 各 group の**展開した局所座標**で線形場を作り、BC・ミラー・時間更新の**前**に作用素だけ比較 (三重周期の角で大域線形場は周期条件を満たさないため)。2/4/8 member (面・辺・角)、非対称 stencil、壁∩継ぎ目、root 交換、原点移動、斜め並進、部分で rank 欠損 → 合併で回復する例 | float32 反例 ($d_0=f32(0.03)-f32(0)$、$d_1=f32(100.03)-f32(100)$) は 1.000000 (root 両順序)。非退化方向の最大誤差 ≤ 1e-5 (相対)。ゼロ成分と**非零定数場の勾配**は絶対誤差 ≤ 1e-6×\|φ\|/h。退化方向は同じ打ち切りの参照解と比較 |
| G1 | $k,\omega$ 勾配: 一様直交格子では解析解、非対称・壁∩継ぎ目・2/4/8 member では **CPU double の合併 GG** (同じ面値規則) を参照。F1 は初回と 2 回目の残差組立で、輸送が実際に読む値を検査 | 参照との最大誤差 ≤ 1e-5 (相対) |
| G2 | 作用素の精度: 参照 = 格納済み float32 座標・場を double で評価した合併 stencil LSQ | $\max_i\lvert\nabla\phi_{gpu}-\nabla\phi_{ref}\rvert \le 10^{-5}\,S$、$S=\max_i\lvert\nabla\phi_{ref}\rvert$ |
| G2' | 細分 3 水準 (二次場)、**ジッタ格子** (内部節点を ±0.2h の決定論的擬似乱数で動かす。周期像は同じ量で動かす。対称 stencil では二次場の勾配が厳密になるため。一様格子は「床以下 = 厳密」を別行で記録) | 誤差床 $e_{floor}=10\,\varepsilon_{f32}S\approx1.2\times10^{-6}S$。継ぎ目・内部の誤差がともに床を超える水準だけで次数を計算し、使える水準が 2 未満なら曲率を 10 倍にして再試験。継ぎ目の節点の誤差 (対 解析勾配) の収束次数 ≥ 0.9 (勾配は O(h)。面再構成の O(h²) と混同しない)、各水準で継ぎ目誤差 / 内部誤差 ≤ 2。ゼロ成分は絶対誤差 ≤ $e_{floor}$ |
| R1 | case/39 周期丘: 整備設定 (`wallTreatmentSST: 0`、実際の `kInit/omegaInit`、段階起動、`check_mesh_quality` PASS、$y_1^+$ 報告) で**旧/新バイナリを同一メッシュ・IC から再生成** | 両 run 同一区間で `check_convergence` **PASS**、$C_f$ 3 点・$x_r$ が `--drift 0.002 --osc 0.005` で STEADY。$C_f=\tau_{w,t}/(\tfrac12\rho_bU_b^2)$ ($\tau_{w,t}$ = `twall` を下壁 $+x$ 接線へ射影、符号は `sern_forces.py` の `twall_on_fluid` と同じ。$\rho_b,U_b$ = 丘頂断面 $y\in[h,3.035h]$ のバルク。z は一意 DOF 平均で継ぎ目重複は重み 1/2) を $x/h=0.5,2,6$ で壁ノード線形補間。$x_r$ = 下壁 $C_f$ の負→正の最初のゼロ交差 ($x/h\in[1,8]$、線形補間、交差なしは「未再付着」で判定不能として R1 は落とさない)。CSV `step,Cf_x05,Cf_x2,Cf_x6,xr_h,r_gradu,r_gradk,r_gradw,dF1_inf` を `check_quasisteady.py --series-csv --drift 0.002 --osc 0.005 --tail 0.4` (OSCILLATING は平均±振幅)。継ぎ目指標 (z 継ぎ目の列 vs 隣接内部列、同じ $x$ 集合、末尾平均): $\lvert\nabla\mathbf u\rvert$・$\lvert\nabla k\rvert$・$\lvert\nabla\omega\rvert$ の L2 比が**新で [0.9, 1.1]**、$F_1$ の差の $L^\infty\le0.05$ (旧は記録のみ、≈2 の想定)。$C_f$ 相対 L2 差・$x_r$ 差は記録 (変わるのが正)。新の低下桁数 ≥ 旧 − 0.5 |
| R2 | case/09 TGV: `procedures/verification/09-taylor-green.md` の非粘性 KEEP 基準 (保存誤差の閾値そのまま、初期総量で正規化) + **粘性 SLAU 2 次** 1 本 (勾配を読む経路)。共通固定 dt、終了時刻 $t = 10\,t_c$、一意 DOF | KEEP は既存閾値 ($\lvert K/K_0-1\rvert\lesssim1\%$、$\lvert\Delta S/S_0\rvert\lesssim10^{-4}$、運動量 $\lesssim10^{-6}$)。**SLAU (Re=1600、定数粘性、unsteady・dual-time なし)**: 新 run で質量 $\lvert M-M_0\rvert/M_0\le10^{-6}$、全運動量 $\lvert P_i\rvert/(\rho_0U_0V)\le10^{-6}$、全エネルギー $\lvert E-E_0\rvert/E_0\le10^{-5}$。KE・$S$ 履歴の旧新差は記録のみ (32³ で継ぎ目節点 ≈9 % なので $10^{-3}$–$10^{-2}$ 級の差が出てよい)。定常 PASS は要求しない |
| R3 | 非周期・軸対称の回帰: case/48 (node、非周期)、case/16 (化学種 GG)、case/44 (軸対称) の固定状態 | 勾配配列が旧新で**ビット同一** (継ぎ目が無ければ係数不変、軸対称は経路不変) |

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-26` | [2026-09-26-boundary-node-periodic-gradient-fix-plan.md](../../notes/reviews/2026-09-26-boundary-node-periodic-gradient-fix-plan.md) | **GO-with-changes**, C0/M4/m1 (合併 LSQ と ransGradient 直後の gather を支持) | **全件採用** (2026-09-26 `diagnostician` 判断)。M1 → 周期像の識別と幾何照合を分離、照合許容 1e-4 h_min、同値類で共通の d_E,w_E (§4.1)。M2 → 適用範囲を node ∧ 並進 ∧ 非軸対称に、同じ条件関数 (§1、§4.1)。M3 → G0/G2 は局所座標で作用素比較、G1 は CPU double の合併 GG を参照 (§6)。M4 → 閾値・区間・再生成条件を具体値で固定 (§6)。m5 → F1 の初回充填を初期化時へ (§4.2)。**k/ω の LSQ 化は後続 plan** (本 plan で入れると R1 の変化の帰属ができない) |
| plan (実装レビュー) | `2026-09-26` | [2026-09-26-boundary-node-periodic-gradient-fix-plan-2.md](../../notes/reviews/2026-09-26-boundary-node-periodic-gradient-fix-plan-2.md) | **GO-with-changes**, C0/M2/m2 | **全件採用** (2026-09-26 `diagnostician` 判断)。M1 → 実変位で組む (§4.1、実装済)。M2 → G2′ の誤差床とジッタ格子、R1 の $C_f$/$x_r$ 抽出規則、R2 の保存上限を具体値で固定 (§6)。m3 → F1 初期値を確保時へ (§4.2、実装済)。m4 → G0 再現物の保存 (§5.1 #2・#5a)。順序 G0/G1/G2′/R3 → R1/R2 |

### 6.2 結果

- **予備確認 (修正後、2026-09-26、正式な G0/G1 ではない)**: 同じ TGV・線形場で 1 step。x 継ぎ目 1458 点 (y・z 継ぎ目から 2.5 格子以上) で
  $\partial U_x/\partial y$ = 1.000000 [0.999995, 1.000006] (修正前 2.000000)、接線成分 0。SST を有効にして $k=1+0.1y$、$\omega=100+2y$ を入れると
  $\partial k/\partial y$ = 0.100000、$\partial\omega/\partial y$ = 2.000000 (内部と一致)。合併し直した group は 2977 (打ち切り 0)。
- **G0 (修正前、2026-09-26)**: 並進、case/09 TGV 32³、線形場 `Ux = 10+y` で x 継ぎ目 1458 点の $\partial U_x/\partial y$ = 2.000000 [1.999990, 2.000012]、
  内部 19683 点 1.000000。`case/09.Taylor-Green/_g0_lsq_seam/G0_translational.txt`。

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/calcGradient_d.cu`、`periodicNode_d.cu`、`ransTransport_d.cu`、`main.cpp`、`mesh/mesh.cpp`。
- **既存の node 周期 run (case/39、case/09、周期翼列など) の結果は継ぎ目付近で変わる** (修正)。再現には修正前の commit のバイナリが要る。
- docs: `methods/gradient.md`。

## 8. 完了条件

- [ ] `methods/gradient.md` の「修正中」を外す
- [ ] 実装・検証完了 (§6)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done` に変更し、§9 に変更ログを記載
- [ ] `plans/active/` → `plans/accepted/` へ移動、[`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-26` — 実装レビュー (codex、GO-with-changes C0/M2/m2) を全件採用。M1 (合併 LSQ を各 incidence の実変位で組む) と m3 (F1 初期値 1 を `allocVariables` へ) を実装、ビルド成功。§6 の G2′・R1・R2 を具体値で固定。

- `2026-09-26` — codex plan 段 GO-with-changes (C0/M4/m1) を全件採用して §4/§6 を改訂。回転周期の書きかけ (mesh.cpp の角の割当) は `notes/sessions/boundary-node-rotational-periodic-wip.patch` に退避。
- `2026-09-26` — ユーザ決定: スカラーの勾配も LSQ に揃える (本 plan の後に別 plan)。本 plan の合併 LSQ 係数は変数に依らない形で作る。
- `2026-09-26` — 初稿。回転周期 plan の G0 で見つかった LSQ の 2 重計上と、codex (同 plan の 2 回目) が見つけた SST $k,\omega$ 勾配の未合算を、
  回転と独立に先に直すため切り出した (`diagnostician` 判断)。
