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
| 化学種 | Green–Gauss (`speciesFaceReconstruction ≥ 1` のみ) | 周期半割面を除外 (`excludePeriodic`) し、合併体積で割った部分寄与の和 | ~~正しい~~ → **撤回 (2026-09-26)**: 除外条件 `ic1 < nCells` はゴースト付与 (`mesh.cpp:443-456`) のため一度も成立しない死にコードで、合併は半割面 (φ[ic0]) 込み。継ぎ目に $\phi(S_a+S_b)/V$ の誤差 (tgv で ≈20–27 ε·φ/h)。本 plan の修正対象 (§4.2) |
| 受動種・凝縮モーメント | Green–Gauss (`passiveGradient`、同じ `species_gradient_d`) | 同上 | ~~正しい~~ → **撤回** (同上、本 plan の修正対象) |
| $\gamma,\ Re_{\theta t}$ | 拡散は 2 点差分 (勾配配列なし)。生成項は速度勾配を読む | 状態は root からミラー | 拡散は影響なし。**生成項は欠陥 1 の影響を受け、欠陥 1 の修正で直る** |

**ユーザ決定 (2026-09-26): スカラーの勾配も LSQ に揃える**。node の勾配は NS の原始変数だけ LSQ (`calcGradient`、`gradLSQ: 2` 固定) で、
$k,\omega$ (`ransGradient`)・化学種 (`species_gradient_d`)・受動種・凝縮モーメント (`passiveGradient`) は Green–Gauss。なお既存の「フォールバック」は
LSQ の退化方向を 0 にするスペクトル打ち切りで、Green–Gauss への差し替えではない (`methods/discretization.md` §7.3.1)。
**進め方**: 本 plan では合併 stencil の LSQ 係数を**どの変数にも使える形**で作る。スカラー勾配の LSQ 統一は、全 node 計算のスカラー勾配を
(継ぎ目以外も) 変える変更なので、本 plan の完了後に別 plan で行う ($k,\omega$ → 化学種 → 受動種・凝縮の順、冷却平板・SERN・凝縮ノズルで回帰。
壁近くの $\omega\sim1/y^2$ と化学種の有界性に注意)。統一すれば欠陥 2 と同じ型 (Green–Gauss の継ぎ目の別扱い) は構造的に無くなる。

## 2. スコープ

- **やる**: 合併 stencil の LSQ 係数の事前計算、$k,\omega$ 勾配の合算位置の修正、並進周期での検証 (G0 拡張、二次場、case/39、case/09)。
- **やる (2026-09-26 追加)**: Green–Gauss のスカラー勾配 ($k,\omega$・化学種・受動種・凝縮モーメント) の周期半割面の除外 (§4.2a。除外条件が死んでいたため)。適用は `periodicSeamMergeActive` のときだけ。
- **やらない**: 回転周期 (`boundary-node-rotational-periodic` で本 plan の後に)。~~Green–Gauss の勾配 (化学種・受動種・凝縮モーメント。現行の和が正しい)~~ (撤回、§4.2a)。軸対称×周期のスカラー勾配 (片側 GG + 半割面込みのまま、既存の未修正挙動)。
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

### 4.2a 周期半割面の除外 (2026-09-26 `diagnostician`、G1-a で発見)

- **欠陥**: `calc_scalar_gradient_face_d` (`ransTransport_d.cu:43`) と `species_gradient_d` (`speciesTransport_d.cu:649`、化学種・tracer・凝縮モーメント・受動種の全呼び出し) の除外条件
  `excludePeriodic && ip >= nNormalPlanes && ic1 < nCells` は、`mesh.cpp:443-456` が周期 bcond にもゴースト (`nCells+nGhost`) を付けるので**一度も成立しない**。
  半割面は境界面として φ[ic0] で積算され、合併後に $\phi(S_a+S_b)/V$ が残る (対の半割面の stored float32 `surfVect` が $2.4\times10^{-6}h^2$ 食い違う。内部面だけの合併閉包は厳密 0)。
  NS の GG は `calcGradient_d.cu:1091` で `bcondKind=="periodic"` をホストで skip しており無事。
- **修正**: `mesh` にホストで作る面フラグ `planePeriodic` (nPlanes byte、全 bcond を走査し `bcondKind=="periodic"` の `iPlanes` に 1、node/cell 問わず) を追加し device へ。
  3 カーネルの条件を `excludePeriodic != 0 && planePeriodic[ip] != 0` に置換し、死に条件は削除。`excludePeriodic` は現行どおり `periodicSeamMergeActive` のときだけ 1
  → 非周期 run と cell 周期はビット不変 (R3 維持)。
- **不採用**: bcond 単位のカーネル分割 (NS 流) — スカラー勾配は全 nPlanes を 1 カーネルで回す構造で変更が大きい。ゴーストの有無・bcond 順での判定 (今回の死に条件と同型)。
- フラグ配列は回転周期 plan でも流用する。

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
| 5a | G0/G1/G2/G2′ ハーネス | #2 の再現物に加え、G1-a の CPU float32 再現、線形 $Y$ (化学種) を 1 本焼いて床を 1 行記録、float32 反例 (期待 1.000000、root 両順序)、ジッタ格子 (±0.2h、決定論的、周期像は同じ量)、CPU double 参照、F1≠1 入力で 1・2 回目に輸送が読む値。合格: §6 G0/G1/G2/G2′ | O |
| 5b (**完了 2026-09-26**、PASS) | R3 | case/48・case/16・case/44 の勾配配列が旧新ビット同一 (F1 初期化で 1 step 目が変わる run は 2 step 目以降のノイズ床比較) | O |
| 5c (**実装済 2026-09-26**、G1-a/G1-b PASS、R3 は #5b) | 周期半割面の除外 (§4.2a) | `mesh` に `planePeriodic`、3 カーネルの条件置換。合格: §6 G1-a/G1-b/定数場を 7 run 再実行、R3 ビット同一。`plans/accepted/species-passive-scalar-unification.md` §4.1-5-1 に訂正 1 行 | O (判断: 2026-09-26 `diagnostician`・本 plan で修正) |
| 5d (**完了 2026-09-26**: R2 収支 成立、R1 延長 限定合格) | R1 延長・R2 収支 (判断: 2026-09-26 `diagnostician`) | R1: 800k 場から `restart_field.py` で +800k step (同一設定・同一閾値、10k–1.6M 全系列 `--tail 0.4`)。なお DRIFTING なら「限定合格 (単調収束中・漸近値・末尾全点で ≤ 0.05)」と書きラベルは変えない。F1 差の L2 と x 継ぎ目の比を参考列に。case/39 の `slauWallNormalChi` 実効値を旧新で控える。R2: 既存 snapshot から $r(t)=(dK/dt+\varepsilon)/(K_0/t_c)$ (測る前に固定: 新は $|r|\le0.05$、旧は $t\lesssim7$ で $r>0$ が持続すれば「継ぎ目由来の注入」、符号不定なら「収支が閉じない」に留める) | O |
| 6a (**完了 2026-09-26**、§6.2) | codex result M1 | `speciesTransport_d.cu:690,1232` の `excludePeriodic` を `periodicSeamMergeActive` に統一。合格: 軸対称×周期の小メッシュで旧新の 1 step 後 roY/roXi がビット同一。§7 に「軸対称×周期・回転周期のスカラー勾配は片側 GG + 半割面込みのまま (既存の未修正挙動)」 | O |
| 6b (**完了 2026-09-26**、§6.2) | codex result M2 | `g2p_jitter.py` を訂正後の基準に揃え、ジッタ 128³ を追加 (AWS)。**測る前に固定**: 最細対 64→128 の継ぎ目の次数 (ゼロ成分含む全成分) ≥ 0.9 かつ 128³ で継ぎ目/内部 ≤ 2。32→64 は補助、16→32 は漸近域外として記録のみ。64→128 が 0.9 未満なら判定不能とし、内部の次数を並べて継ぎ目固有でないことだけ書く | O |
| 6c (**完了 2026-09-26**、§6.2) | codex result M3 | `g_suite.py` の G1-b から n_member 倍の自動緩和を削除 (上限 $2N_{max}\varepsilon\phi/h$ かつ継ぎ目/内部 ≤ 2 を機械判定)、GPU 定数場 (k・ω・ξ = const、継ぎ目勾配 ≤ 4ε\|φ\|/h) を追加、7 run 再実行 (AWS)。化学種の直接確認は tracer ξ で代理 (根拠: `passiveGradient_d_wrapper` は同じ `species_gradient_d` を同じ excludePeriodic で呼ぶ、`speciesTransport_d.cu:1220-1232`) | O |
| 6d (**完了 2026-09-26**、§6.2) | codex result M5 | R1 の記述訂正。旧は元 run + 延長を `stage_manifest` で同一実効設定の 1 区間として連結し `--segment` (スパイク込み)。F1 仮説: 同じ restart 入力で旧バイナリ 1 step を `sstSigmaBlend` 0/1 で比較 (スパイクが 1 側だけなら仮説支持、そうでなければ「原因未特定」)。新 dF1_inf の定常性は継続課題 (#7) | O |
| 6e (**完了 2026-09-26**: 判定行を置換、`r2_ke_budget.txt` 更新、診断 ALL HOLD) | codex result M4 | R2 の結論を観測と帰属に限定し、`r2_ke_budget.py` の VERDICT 行を符号条件 (a)(b)(c) に置換、旧 \|r\| ≤ 0.05 行は履歴 | O |
| 6f (**完了 2026-09-26**: `methods/gradient.md`・§2・`plans/README.md`) | codex result m1 | `methods/gradient.md` を実装後の記述に (欠陥は履歴節)、§2 スコープ、`plans/README.md`、§5.1 #5d、case/09・case/39 の run 表 | O |
| 6g | codex result 2 回目 | 6a–6f の後。`--focus` 「M1–M5 の閉じ方、R1 の限定合格 + 継続課題の扱い、R2 の観測限定の文言」 | O (結論 F) |
| 7 | 継続課題 (accepted 後も残す) | 新 R1 の dF1_inf の定常性 (DRIFTING 0.3 %/tail、上限内)。R2 の機構分解 (同一状態の旧新離散残差から KE 仕事を分解) | O |
| 5 (**完了 2026-09-26**、§6.2。codex result で #6a–#6g を追加) | 検証 R1・R2 | §6 R1・R2 (5a・5b の後)。**区切りで codex** | O (結論 F) |
| 6 | docs + codex result | `methods/gradient.md` の「修正中」を外す | F |

## 6. 検証 (測る前に固定)

| # | 試験 | 合格 (測る前に固定) |
| --- | --- | --- |
| G0 | LSQ の局所作用素試験: 各 group の**展開した局所座標**で線形場を作り、BC・ミラー・時間更新の**前**に作用素だけ比較 (三重周期の角で大域線形場は周期条件を満たさないため)。2/4/8 member (面・辺・角)、非対称 stencil、壁∩継ぎ目、root 交換、原点移動、斜め並進、部分で rank 欠損 → 合併で回復する例 | float32 反例 ($d_0=f32(0.03)-f32(0)$、$d_1=f32(100.03)-f32(100)$) は 1.000000 (root 両順序)。非退化方向の最大誤差 ≤ 1e-5 (相対)。ゼロ成分と**非零定数場の勾配**は絶対誤差 ≤ 1e-6×\|φ\|/h。退化方向は同じ打ち切りの参照解と比較 |
| G1 | $k,\omega$ 勾配と化学種・tracer の GG。**訂正の履歴**: 当初の「CPU double 参照に相対 1e-5」を 2026-09-26 に「float32 GG の積算床」として書き直したが、**これは誤り** (`diagnostician` 判断を含む、機序の取り違えでスケールが偶然近かった)。G1-a の CPU float32 再現により、周期半割面の除外条件が死んでおり継ぎ目誤差の主因は対の半割面の閉包差 $\phi(S_a+S_b)/V$ だと特定した (§4.2a、#5c)。**G1-a**: CPU で同じ面値規則・同じ float32 演算 (半割面除外、部分和を float32 で積算 → 合併体積で除算 → 和・broadcast) を再現して GPU と比較。**G1-b**: CPU double の合併 GG と比較。**定数場**: φ = const で継ぎ目の勾配を見る。一様直交・非対称・壁∩継ぎ目・2/4/8 member・mirror の 7 run。F1 は初回と 2 回目の残差組立で輸送が実際に読む値を検査 | G1-a: 差 ≤ $4\varepsilon_{f32}\max|\phi|/h$ (継ぎ目・内部とも)。通らなければ実装のバグであり床と呼ばない。G1-b: $2N_{max}\varepsilon_{f32}\max|\phi|/h$ (面数から導く上限) 以内**かつ**継ぎ目/内部の誤差比 ≤ 2 (内部が厳密 0 の一様格子成分は継ぎ目 ≤ 4 ε)。定数場: 継ぎ目 ≤ $4\varepsilon|\phi|/h$。修正前の床の表は「修正前」として結果欄に残す |
| G2 | 作用素の精度: 参照 = 格納済み float32 座標・場を double で評価した合併 stencil LSQ | $\max_i\lvert\nabla\phi_{gpu}-\nabla\phi_{ref}\rvert \le 10^{-5}\,S$、$S=\max_i\lvert\nabla\phi_{ref}\rvert$ |
| G2' | 細分 3 水準 (二次場)、**ジッタ格子** (内部節点を ±0.2h の決定論的擬似乱数で動かす。周期像は同じ量で動かす。対称 stencil では二次場の勾配が厳密になるため。一様格子は「床以下 = 厳密」を別行で記録) | 誤差床 $e_{floor}=10\,\varepsilon_{f32}S\approx1.2\times10^{-6}S$。継ぎ目・内部の誤差がともに床を超える水準だけで次数を計算し、使える水準が 2 未満なら曲率を 10 倍にして再試験。継ぎ目の節点の誤差 (対 解析勾配) の収束次数 ≥ 0.9 (勾配は O(h)。面再構成の O(h²) と混同しない)、各水準で継ぎ目誤差 / 内部誤差 ≤ 2。**ジッタ格子ではゼロ成分も同じ判定** (非対称 stencil では O(h))。「ゼロ成分は絶対誤差 ≤ $e_{floor}$」は**一様格子の行のみ** (2026-09-26 訂正: 当初の書き方は一様格子前提の誤指定) |
| R1 | case/39 周期丘: 整備設定 (`wallTreatmentSST: 0`、実際の `kInit/omegaInit`、段階起動、`check_mesh_quality` PASS、$y_1^+$ 報告) で**旧/新バイナリを同一メッシュ・IC から再生成** | 両 run 同一区間で `check_convergence` **PASS**、$C_f$ 3 点・$x_r$ が `--drift 0.002 --osc 0.005` で STEADY。$C_f=\tau_{w,t}/(\tfrac12\rho_bU_b^2)$ ($\tau_{w,t}$ = `twall` を下壁 $+x$ 接線へ射影、符号は `sern_forces.py` の `twall_on_fluid` と同じ。$\rho_b,U_b$ = 丘頂断面 $y\in[h,3.035h]$ のバルク。z は一意 DOF 平均で継ぎ目重複は重み 1/2) を $x/h=0.5,2,6$ で壁ノード線形補間。$x_r$ = 下壁 $C_f$ の負→正の最初のゼロ交差 ($x/h\in[1,8]$、線形補間、交差なしは「未再付着」で判定不能として R1 は落とさない)。CSV `step,Cf_x05,Cf_x2,Cf_x6,xr_h,r_gradu,r_gradk,r_gradw,dF1_inf` を `check_quasisteady.py --series-csv --drift 0.002 --osc 0.005 --tail 0.4` (OSCILLATING は平均±振幅)。継ぎ目指標 (z 継ぎ目の列 vs 隣接内部列、同じ $x$ 集合、末尾平均): $\lvert\nabla\mathbf u\rvert$・$\lvert\nabla k\rvert$・$\lvert\nabla\omega\rvert$ の L2 比が**新で [0.9, 1.1]**、$F_1$ の差の $L^\infty\le0.05$ (旧は記録のみ、≈2 の想定)。$C_f$ 相対 L2 差・$x_r$ 差は記録 (変わるのが正)。新の低下桁数 ≥ 旧 − 0.5 |
| R2 | case/09 TGV: `procedures/verification/09-taylor-green.md` の非粘性 KEEP 基準 (保存誤差の閾値そのまま、初期総量で正規化) + **粘性 SLAU 2 次** 1 本 (勾配を読む経路)。共通固定 dt、終了時刻 $t = 10\,t_c$、一意 DOF | KEEP は既存閾値 ($\lvert K/K_0-1\rvert\lesssim1\%$、$\lvert\Delta S/S_0\rvert\lesssim10^{-4}$、運動量 $\lesssim10^{-6}$)。**SLAU (Re=1600、定数粘性、unsteady・dual-time なし)**: 新 run で質量 $\lvert M-M_0\rvert/M_0\le10^{-6}$、全運動量 $\lvert P_i\rvert/(\rho_0U_0V)\le10^{-6}$、全エネルギー $\lvert E-E_0\rvert/E_0\le10^{-5}$。KE・$S$ 履歴の旧新差は記録のみ (32³ で継ぎ目節点 ≈9 % なので $10^{-3}$–$10^{-2}$ 級の差が出てよい)。定常 PASS は要求しない |
| R3 | 非周期・軸対称の回帰: case/48 (node、非周期)、case/16 (化学種 GG)、case/44 (軸対称) の固定状態 | 勾配配列が旧新で**ビット同一** (継ぎ目が無ければ係数不変、軸対称は経路不変)。**判定の具体化 (2026-09-26、測定前)**: 同一入力状態から 1 step、`output.level: 2` の `res_1.h5` の全勾配・リミタ配列を旧 2 回・新 2 回で比較。旧同士がビット一致する配列は旧新もビット一致を要求。旧同士でも一致しない配列 (GG の atomicAdd 順序、`case/48/run_0030/0031` で 1 step の場が非再現と実測済み) は、旧新の最大差 ≤ 旧同士の最大差 × 2 かつ不一致数の桁が同じ。case/44 は入力メッシュが削除済みのため、同じ軸対称 node 経路を通る小さな軸対称メッシュで代替 (理由を結果欄に書く) |

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-26` | [2026-09-26-boundary-node-periodic-gradient-fix-plan.md](../../notes/reviews/2026-09-26-boundary-node-periodic-gradient-fix-plan.md) | **GO-with-changes**, C0/M4/m1 (合併 LSQ と ransGradient 直後の gather を支持) | **全件採用** (2026-09-26 `diagnostician` 判断)。M1 → 周期像の識別と幾何照合を分離、照合許容 1e-4 h_min、同値類で共通の d_E,w_E (§4.1)。M2 → 適用範囲を node ∧ 並進 ∧ 非軸対称に、同じ条件関数 (§1、§4.1)。M3 → G0/G2 は局所座標で作用素比較、G1 は CPU double の合併 GG を参照 (§6)。M4 → 閾値・区間・再生成条件を具体値で固定 (§6)。m5 → F1 の初回充填を初期化時へ (§4.2)。**k/ω の LSQ 化は後続 plan** (本 plan で入れると R1 の変化の帰属ができない) |
| plan (実装レビュー) | `2026-09-26` | [2026-09-26-boundary-node-periodic-gradient-fix-plan-2.md](../../notes/reviews/2026-09-26-boundary-node-periodic-gradient-fix-plan-2.md) | **GO-with-changes**, C0/M2/m2 | **全件採用** (2026-09-26 `diagnostician` 判断)。M1 → 実変位で組む (§4.1、実装済)。M2 → G2′ の誤差床とジッタ格子、R1 の $C_f$/$x_r$ 抽出規則、R2 の保存上限を具体値で固定 (§6)。m3 → F1 初期値を確保時へ (§4.2、実装済)。m4 → G0 再現物の保存 (§5.1 #2・#5a)。順序 G0/G1/G2′/R3 → R1/R2 |
| G1 閾値の判断 | `2026-09-26` | (本 plan §6 G1、`_g0_lsq_seam/G0_translational_m1.txt`) | `diagnostician`: G1 の 1e-5 は float32 GG の床の見落とし | 採用。GG カーネルは変えない。G1 を G1-a (float32 再現でバグ検出) と G1-b (面数から導いた床) に分割。§1 の「正しい」を合併の意味に限定 |
| G1-a FAIL の判断 | `2026-09-26` | (本 plan §4.2a、`_g0_lsq_seam/G_tgv.txt` ほか) | `diagnostician`: 前行の「床」判断は誤り。真因は周期半割面除外の死に条件 | 採用。§4.2a の面フラグで修正 (#5c)、G1 を修正後基準に、G2′ ゼロ成分はジッタ格子で次数判定、root 交換は mirror (点反転) で代替 (root は union-find の最小 index で bcond 順に依らない、`mesh.cpp:704-712`) |
| R1/R2 解釈 | `2026-09-26` | (本 plan §6.2) | `diagnostician`: dF1_inf DRIFTING は 1 回延長、R2 は KE 収支で帰属を確定、抽出の [解釈] 4 点は注記つきで承認 | 採用 (#5d)。R1 の Cf・$x_r$ 変化は「本 plan の修正 (旧新差は本 plan の 7 commit のみ) の効果」と書き「改善」とは書かない。$x_r$ を LES との距離で評価しない。codex result は #5d の後 |
| R2 収支基準 | `2026-09-26` | (本 plan §6.2 R2) | `diagnostician`: 新 $|r|\le0.05$ は数値散逸を見落とした誤指定 | 測定後の訂正として記録し、符号判定 (a)(b)(c) に置換。旧は「注入」と書く (根拠: $K/K_0>1$ と $r>0$ の連続、Π 差し引き後)。0.05 を緩めて通すことはしない。任意で 64³ の新 run で $|r|_{max}$ の格子依存を記録 (§7) |
| result | `2026-09-26` | [2026-09-26-boundary-node-periodic-gradient-fix-result.md](../../notes/reviews/2026-09-26-boundary-node-periodic-gradient-fix-result.md) | **NO-GO**, C0/M5/m1 (中核実装は支持、accepted 不可) | **全件採用** (2026-09-26 `diagnostician`)。M1 → 化学種・受動種の除外条件を `periodicSeamMergeActive` に統一 (#6a)。M2 → G2′ ゼロ成分 16→32 0.849 未達を戻し、ジッタ 128³ を追加 (#6b)。M3 → G1-b の n_member 自動緩和を削除、比 ≤ 2 と GPU 定数場を機械判定、7 run 再実行 (#6c)。M4 → R2 は観測と帰属に限定し機構語を削除、符号条件は「測定後に追加した診断」(#6e、前回の「注入と書いてよい」を一部撤回)。M5 → 「単調収束中・漸近値」を撤回、旧延長は連結区間で判定、F1 仮説を切り分け (#6d)。m1 → docs・スコープ・README 同期 (#6f)。新の dF1_inf DRIFTING は、(i) 旧連結区間 PASS (ii) F1 切り分け実施 (iii) 上限成立 を満たせば accepted の障害にしない (継続課題として §5.1 に残す) |
| #6a–#6d の判断 | `2026-09-26` | (本 plan §6.2「codex result 対応」) | `diagnostician`: #6a は規則上 FAIL を記録し決定論的試験で閉じる、#6c は比を場ごと・定数場の壁節点を除外 (定義の補正)、#6b は判定不能で閉じ障害にしない、#6d は継続課題で可 | 採用。#6a の面流束・状態の旧新ビット一致を確認 (`R6a_massflux_bitcheck.txt`)。規則を分布ベースに直す・閾値 4ε を上げる・#6b を粗さで PASS にすることはしない |

### 6.2 結果

- **予備確認 (修正後、2026-09-26、正式な G0/G1 ではない)**: 同じ TGV・線形場で 1 step。x 継ぎ目 1458 点 (y・z 継ぎ目から 2.5 格子以上) で
  $\partial U_x/\partial y$ = 1.000000 [0.999995, 1.000006] (修正前 2.000000)、接線成分 0。SST を有効にして $k=1+0.1y$、$\omega=100+2y$ を入れると
  $\partial k/\partial y$ = 0.100000、$\partial\omega/\partial y$ = 2.000000 (内部と一致)。合併し直した group は 2977 (打ち切り 0)。
- **G0 (修正前、2026-09-26)**: 並進、case/09 TGV 32³、線形場 `Ux = 10+y` で x 継ぎ目 1458 点の $\partial U_x/\partial y$ = 2.000000 [1.999990, 2.000012]、
  内部 19683 点 1.000000。`case/09.Taylor-Green/_g0_lsq_seam/G0_translational.txt`。

- **ハーネス #5a (2026-09-26、HEAD `1f03c036`、半割面除外の修正前)**: G0 (2/4/8 member・mirror による root 交換・原点移動・ジッタ・壁∩継ぎ目・float32 反例) PASS、G2 ≤ 1.4e-7·S PASS、
  G2′ 次数 0.95/0.93・継ぎ目/内部 ≤ 1.00 PASS、F1 読み取り PASS。**G1-a FAIL** (継ぎ目 4.7–27 ε·max|φ|/h) → §4.2a の欠陥を特定。G2′ ゼロ成分は基準の誤指定 (§6 訂正)。
  修正前の床 (tgv、GPU − CPU double、ε·max|φ|/h): 内部 dK 1.1 / dΩ 1.3 / dξ 0.9、**継ぎ目 dK 20.8 / dΩ 26.5 / dξ 22.7**、LSQ dUx 0.4/0.4。
- **§4.2a 修正後 (2026-09-26、バイナリ sha256 `88e949fa…`)**: 7 run (tgv / tgv_mirror / tgv_shift / tgv_repeat / tgv_bcswap / jitter32 / channel) で G0・G2・G1-a・G1-b すべて PASS。
  継ぎ目の床は tgv で dK 1.7 / dΩ 2.4 / dξ 1.5 (修正前 20.8 / 26.5 / 22.7)。**継ぎ目/内部比は全 run 0.83–1.85 (≤ 2)**。G1-a の最大差は ≤ 1.7 ε (≤ 4 ε)。
  結果 `case/09.Taylor-Green/_g0_lsq_seam/G_*.txt`。**未実施**: GPU の定数場単独試験 (§6 G1、#5a に残す)。

- **R3 (2026-09-26、旧 `1266aba1` sha256 `9b45d018…` / 新 sha256 `88e949fa…`、各 1 step・旧 2・新 2 run、`case/09.Taylor-Green/_g0_lsq_seam/R3.txt`)**: **PASS**。
  case/48 (node 非周期 SST): 勾配・リミタ 26 配列ビット一致、dK/dΩ 4 配列は旧同士の atomicAdd 床と同水準。case/16 (化学種 GG 有効): 30 配列ビット一致。
  軸対称 (case/44 の入力メッシュ削除済みのため 101×41 の軸対称 node メッシュで代替、品質 PASS): 30 配列ビット一致。NaN 0。
  **既定 (`sstSigmaBlend` 1) の case/48 は step 1 の roK/roOmega が旧新で約 5 万点違う** (旧同士 0/2 点) → m3 (F1 初期値 1) の予告どおりの変化。
  帰属確認: 同じ入力で `sstSigmaBlend: 0` にすると旧新差は roK 1 点・roOmega 3 点 (旧同士 1/0 点) に消える (スクラッチ `r3b_case48_*`)。
  **未測定**: ∇Y は出力変数に無く直接比較していない (roY の旧新差が旧同士と同水準であることのみ確認)。case/48 は一様 IC のため dP/dT/dρ とリミタは検出力なし。

- **R2 (2026-09-26、AWS g5、旧 `1266aba1` sha256 `f0505fe8…` / 新 `565959c7` sha256 `0bcff10d…`、block 128、`case/09.Taylor-Green/run_0173..0176`、`r2_conservation.txt`)**:
  32³ node (品質 PASS)、RK4 dt 0.007 × 3572 step (t = 10 t_c)。**4 本とも §6 R2 の保存上限 PASS**、NaN 0。
  KEEP 非粘性 旧/新: |K/K0−1| ≤ 6.8e-3、|ΔS/S0| ≤ 8.4e-6、運動量 ≤ 1.6e-8、旧新差 max|ΔK/K0| 6.0e-8。
  SLAU 粘性 (Re 1600、κ は Pr 0.71 から 1.408e-4、`slauWallNormalChi: 0`) 新: 質量 2.3e-8、運動量 1.2e-8、全エネルギー 8.3e-9。
  **旧新差 (記録のみ)**: max|ΔK/K0| 1.0e-1 (t = 10.85)、最終 4.5e-2、max|Δ(ΔS/S0)| 3.4e-3 — §6 の見込み 1e-3〜1e-2 を超える。
  観測: 旧 SLAU は t ≲ 7 で K/K0 が 1.02 まで増え ΔS < 0、新は単調減少 (`r2_ke_entropy.png`)。**解釈は未確定** (R1 と合わせて `diagnostician`、条件 7)。
  **KE 収支 (#5d、`case/09.Taylor-Green/r2_ke_budget.{py,csv,txt,png}`)**: $r=(dK/dt-\Pi+\varepsilon)/(K_0/t_c)$ ($\varepsilon$ は解像粘性散逸、$\Pi=\int p\nabla\cdot u$)。
  新: $r\le-0.0074$ (全時刻負、$t\le7$ は −0.0074〜−0.0447、$|r|_{max}$ 0.0885 @ t≈13 = 散逸ピーク)。旧: $t\le4.5$ で +0.0096〜+0.0118 が連続、$K/K_0$ 最大 1.0187。
  ~~判定: 新 $|r|\le0.05$~~ → **訂正 (2026-09-26、測定後の訂正であることを明記)**: 当初基準は風上スキーム (SLAU・MUSCL) の数値散逸が $\varepsilon$ に入らないことを見落とした誤指定 (`diagnostician`)。
  符号で判定し直す: (a) 新 $r\le+0.005$ 全時刻 → 最大 −0.0074 **成立**、(b) $r_{old}-r_{new}>0$ が $t\le7$ の全 snapshot → 最小 +0.0176・平均 +0.0257 **成立**、
  (c) 旧の $K/K_0>1$ 区間で $r>0$ が連続 ≥ 10 snapshot → 15 連続 **成立**。$|r|$ の大きさはゲートにしない (数値散逸の量は本 plan の対象外)。
  ~~**解釈**: 旧は継ぎ目由来の非物理なエネルギー注入を持ち、修正後は注入が消えて数値散逸のみ~~ (撤回、codex result M4: $r$ は独自の中心差分と snapshot 差分による後処理で、作用素の不一致と時間微分誤差を含むので機構は言えない)。
  **結論 (観測と帰属に限定)**: 修正後に $K/K_0>1$ の超過が消え ($\max K/K_0$ 1.0187 → 1.0000)、この後処理での収支残差 $r$ が全時刻で負になった (旧は $K/K_0>1$ 区間で $r>0$ が 15 snapshot 連続)。差は本 plan の commit 範囲に帰属。符号条件 (a)(b)(c) は測定後に追加した診断であり保存性ゲートの代替ではない。機構の分解は #7。
- **R1 最終 (2026-09-26、S2 800k + 延長 800k = 連結 10k–1.6M、`case/39.periodic_hills/run_0036..0039`、AWS g5・block 128)**:
  収束: S1→S2 区間 `--segment` 旧新とも **PASS** (roK 旧 4.6 / 新 4.3 桁 ≥ 旧 − 0.5)。延長 (`--from-floor` 各 800k run) 新 **PASS**、旧 NOT CONVERGED — 延長開始直後 step 0 の rms_roK のみ床の 408 倍 (1.34e-6、新は 2.9e-9 で床)、末尾は床の 0.94 倍。
  **仮説 (未検証)**: 旧バイナリは `buildScalarDescs` の初回呼び出しで計算済みの sstF1 を 1 で上書きする (m3 で修正した挙動) ので、restart 直後の 1 step が乱れる。
  準定常 (`--drift 0.002 --osc 0.005 --tail 0.4`、160 snap): 旧 全列 STEADY。新 Cf・$x_r$・継ぎ目比は STEADY、**dF1_inf のみ DRIFTING 0.3 %/tail** (800k 時点 1.6 % から減少)。
  → **限定合格 (差の大きさに対して)**: ~~dF1_inf は単調収束中、漸近値 ≈0.0335~~ (撤回、codex result M5: 末尾 64 点は増加 33・減少 30 で単調でない) → 微小な増加傾向と変動が残るが、観測窓内で上限 0.05 は成立 (0.03334–0.03345)。STEADY は未達でラベルは変えない。定常性は継続課題 (#7)。
  1.6M の値 (旧 → 新): Cf(0.5) −7.398e-3 → −7.396e-3、Cf(2) −4.578e-3 → −4.539e-3、Cf(6) 5.29e-4 → 3.94e-4、$x_r/h$ 4.811 → 5.034、U_b 26.20 → 26.49 m/s。
  z 継ぎ目比 (∇u, ∇k, ∇ω) 2.24 / 1.58 / 0.50 → 0.997 / 0.984 / 1.000 (基準 [0.9, 1.1] 成立)、dF1_inf 0.553 → 0.0334 (≤ 0.05 成立)。
  参考列 (判定外): F1 差 L2 0.058 → 0.0019、x 継ぎ目 (丘頂、両隣列平均との比) 1.87 / 0.63 / 0.50 → 0.94 / 0.99 / 1.00。
  壁解像 PASS (y1+ 最大 0.75 / 0.74、超過 0 %)。`slauWallNormalChi` 実効値 旧新とも 0 (explicit)。
  **解釈 (`diagnostician` 2026-09-26)**: 旧の継ぎ目比は LSQ 2 重計上と k/ω 片側 GG の指紋で、新で消えた。Cf(6) −26 %・$x_r$ 4.81→5.03 は本 plan の修正 (旧新差は 7 commit のみ) による継ぎ目の偽勾配の除去の効果。「改善」とは書かず、LES 参照との距離で評価しない。抽出の [解釈] 4 点 (F1 の再計算と 1 step ずれ、隣接 = 第 1 内部層、x 継ぎ目の除外、L2 比) は `r1_extract.py` の docstring どおり承認。
- ~~**R1 (実行中)**: `case/39.periodic_hills/run_0036_r1_gradfix_old` / `run_0037_r1_gradfix_new` (AWS `~/forge-pgrad-new/`)、メッシュ 80×50×30 y1 1.5e-3 h (品質 PASS、AR 149)、
  1 次 SST 定常 (S1 静止スピンアップ 2000 → S2 成形 IC 800k step)。先行 300k の暫定: 収束 PASS 旧新とも、新の Cf_x6・dF1・r_gradk は DRIFTING (単調減衰)。継ぎ目比 旧 r_gradu 2.24 / r_gradk 1.58 / r_gradw 0.50 / dF1 0.55 → 新 1.00 / 0.99 / 1.00 / 0.03。壁解像 局所 y1+ 最大 0.75、超過 0 %。最終判定は 800k の結果で行う。~~ (上の R1 最終で置き換え)

- **codex result 対応 (#6a–#6d、2026-09-26、新 bd22376d sha256 `f0ad00c8…` / 旧 1266aba1、AWS・block 128、commit 9411df0e)**:
  - **#6a 軸対称×並進周期** (`_g0_lsq_seam/R6a_axi_periodic.txt`): res_0 は 88 配列が 4 run でビット一致。res_1 は R3 規則で **FAIL** (roUx 不一致数 旧同士 85 / 旧新 132、roY1 6 / 11。最大差は旧同士と同値)。
    規則は書き換えず、**決定論的試験で閉じた** (`R6a_massflux_bitcheck.txt`): 初回面流束 5104 面・カーネルが読んだ状態 15006 値が旧新で**ビット一致**、`speciesTransport_d.cu` の差分は除外条件の引数とコメントのみ。
    → 面作用素は同一、場の不一致数の差は atomicAdd 集積順序の統計差。
  - **#6c G1** (7 変種、`G_*.txt`): G0・G2・G1-a・G1-b (上限 $2N_{max}\varepsilon\phi/h$ のみ、自動緩和削除) すべて PASS、G1-b 最大差/閾値 ≤ 0.161。
    **定義の補正 (測定後、理由つき)**: (i) 継ぎ目/内部比は**場ごと** (全成分の最大誤差の比)。一様構造格子では内部の一部成分の誤差が対称性で丸め床になり、成分比は不良条件。成分比 (tgv 系 ω の x 成分 2.13–2.15) は記録のみ。場ごとの比は最大 1.85。
    (ii) GPU 定数場の「継ぎ目」区分から**壁節点を除外** (k/ω と同じ扱い)。channel の ξ は壁∩継ぎ目 96 点・非継ぎ目壁 210 点とも最大 5.31 ε|φ|/h で同値 = 壁半割面込み GG の既存の閉包床で継ぎ目由来でない。壁でない継ぎ目は 0.82 → PASS。
    ハーネスの `G_channel.txt` の VERDICT 行は壁節点込みの FAIL のまま (判定は本節の区分による)。定数場は他 6 変種 PASS (最大 3.03)。
  - **#6b G2′** (`G2p_jitter.txt`、ジッタ 16/32/64/128³): 継ぎ目の次数 (全成分最小) 16→32 0.499 (記録のみ)、32→64 0.876、**64→128 0.855** → 登録どおり**判定不能**。
    内部の次数 0.825–1.000、128³ の継ぎ目/内部 0.77–0.87 → 継ぎ目は内部より悪くない (次数が 0.9 を跨ぐのはジッタ格子での LSQ 自体の性質)。G2 は 128³ で 9.9e-6·S (閾値 1e-5、余裕 1 %)。accepted の障害にしない。
  - **#6d R1**: 連結区間 (S1→S2→S3_ext を 1 区間、`CONVERGENCE_VERDICT_concat.txt`) 旧 PASS (roK 4.6 桁)・新 PASS (4.3 桁)。
    F1 切り分け (`_f1split/F1_SPLIT.txt`): 旧バイナリの restart step 0 の rms_roK は sstSigmaBlend 1 / 0 とも床の 408.7 倍 → **F1 仮説は不支持、原因未特定** (旧バイナリ固有の restart スパイク。新バイナリの自場 restart ではスパイクなし)。
    新の dF1_inf の DRIFTING は継続課題 (#7) とし accepted の障害にしない (条件 (i)(ii)(iii) 成立)。

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/calcGradient_d.cu`、`periodicNode_d.cu`、`ransTransport_d.cu`、`main.cpp`、`mesh/mesh.cpp`。
- **既存の node 周期 run (case/39、case/09、周期翼列など) の結果は継ぎ目付近で変わる** (修正)。再現には修正前の commit のバイナリが要る。
- docs: `methods/gradient.md`。
- 壁半割面込み GG の定数場閉包は float32 で約 5 ε|φ|/h (channel 実測、継ぎ目に依らない) で未対応。
- ジッタ格子での LSQ 勾配の次数は 0.83–1.0 (内部も同じ) で、0.9 を跨ぐ。後続の LSQ 統一 plan で扱う。
- 軸対称×周期・回転周期のスカラー勾配は片側 GG + 半割面込みのまま (既存の未修正挙動、本修正で不変)。
- GG を残す変数 (化学種・受動種・凝縮モーメント) で $N\varepsilon\phi/h$ 床が問題になるなら、$(\phi_f-\phi_i)S_f$ 形を検討する (本 plan では変えない: 全域の GG ビットが変わり R3 の帰属が濁る。k/ω は後続の LSQ 統一 plan で床ごと消え、その plan の §1 にこの床を根拠として引用する)。

## 8. 完了条件

- [ ] `methods/gradient.md` の「修正中」を外す
- [ ] 実装・検証完了 (§6)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done` に変更し、§9 に変更ログを記載
- [ ] `plans/active/` → `plans/accepted/` へ移動、[`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-26` — codex result 対応 #6a–#6f 完了。#6a は決定論的試験 (面流束ビット一致) で閉じ、#6c は比・定数場区分の定義を補正、#6b は判定不能で閉じた (`diagnostician`)。codex result 2 回目へ。

- `2026-09-26` — codex result NO-GO (C0/M5/m1) を全件採用 (#6a–#6g)。R2 の機構の断定と R1 の「単調収束中」を撤回。
- `2026-09-26` — R1 延長 (+800k) 完了: 新は dF1_inf のみ DRIFTING 0.3 %/tail で限定合格 ~~、他の R1 基準はすべて成立~~ (旧延長の収束が NOT CONVERGED、#6d)。R2 は KE 収支の符号判定で成立。検証 (§6) は出そろい、codex result 待ち。

- `2026-09-26` — R3 PASS (3 ケース、NS 勾配・リミタはビット一致、GG k/ω は atomicAdd 床と同水準)。既定 SST で step 1 の roK/roOmega が変わるのは m3 (F1 初期値) によることを `sstSigmaBlend: 0` の A/B で確認。

- `2026-09-26` — §4.2a 実装 (`mesh` の面フラグ `planePeriodic_d`、`calc_scalar_gradient_face_d`・`species_gradient_d` の条件置換)。7 run で G1-a/G1-b PASS、継ぎ目/内部比 ≤ 1.85。

- `2026-09-26` — ハーネス #5a: G0・G2・G2′ 次数・F1 読み取り PASS。**G1-a FAIL** から、周期半割面の除外条件 (`ic1 < nCells`) が死にコードだったと特定。直前の行の「float32 床」は**誤りと訂正** (記録は残す)。§4.2a で修正方針を決定 (`diagnostician`)。

- `2026-09-26` — G0 ハーネス (M1 後): LSQ は継ぎ目・内部とも ε 級で PASS。GG の k/ω は内部でも G1 当初閾値を割る float32 床 ($3$–$61\,\varepsilon\phi/h$) を確認し、`diagnostician` 判断で G1 を書き直した (カーネル不変)。

- `2026-09-26` — 実装レビュー (codex、GO-with-changes C0/M2/m2) を全件採用。M1 (合併 LSQ を各 incidence の実変位で組む) と m3 (F1 初期値 1 を `allocVariables` へ) を実装、ビルド成功。§6 の G2′・R1・R2 を具体値で固定。

- `2026-09-26` — codex plan 段 GO-with-changes (C0/M4/m1) を全件採用して §4/§6 を改訂。回転周期の書きかけ (mesh.cpp の角の割当) は `notes/sessions/boundary-node-rotational-periodic-wip.patch` に退避。
- `2026-09-26` — ユーザ決定: スカラーの勾配も LSQ に揃える (本 plan の後に別 plan)。本 plan の合併 LSQ 係数は変数に依らない形で作る。
- `2026-09-26` — 初稿。回転周期 plan の G0 で見つかった LSQ の 2 重計上と、codex (同 plan の 2 回目) が見つけた SST $k,\omega$ 勾配の未合算を、
  回転と独立に先に直すため切り出した (`diagnostician` 判断)。
