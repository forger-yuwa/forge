forge (自作の圧縮性 FVM ソルバ。CUDA/float32、cell 中心と node 中心 median-dual の 2 離散化、現在は node 主体。
SLAU/Roe/KEEP、block-DPLUR 陰解法、SST、多成分 TP、凝縮、軸対称、ノズル設計ツール design/forge_design を含む) の
リポジトリに対する**外部レビュー**を依頼する。忖度なしで、主張はコードと実測 (run の数値) で検証すること。
結論が「この計画/結果は誤り」でも構わない。両論併記で逃げず、推奨は 1 つに絞ること。

ルール:
- **ファイルを変更しない** (read-only サンドボックスで動いている。読む・実行して確認するのは可)。
- 出力は日本語。識別子・ファイル名は原語のまま。
- 指摘は **Critical / Major / Minor** の重大度付きで、必ず根拠 (`ファイル:行` または `run_*` の数値) と対案をセットで書く。
- リポジトリのルールは `AGENTS.md`、現在仕様は `methods/`、運用手順は `procedures/`、設計判断は `plans/`。
  用語や設定の意味は推測せず `procedures/solver-settings.md` / `procedures/recommended-settings.md` を読むこと。
- 収束の判定は `solver_density_cuda/tools/check_convergence.py <run_dir>` (各 run の `CONVERGENCE_VERDICT.txt`)、
  派生量の定常性は `check_quasisteady.py` の VERDICT を根拠にする。`rms_ro` 単独やスナップショット 1 枚で判断しない。

## 依頼: 検証結果レビュー (stage = result)

対象の plan は下に全文を貼る (`plans/active/boundary-node-periodic-gradient-fix.md`)。実装と検証が終わり、`status: done` にして `plans/accepted/` へ移す直前の段階である。
次を順に評価せよ。

1. **実装 diff の検証**: `git diff 1266aba1...HEAD -- solver_density_cuda design methods procedures` (必要なら `git log 1266aba1..HEAD --oneline`)
   を自分で取り、plan §4 の設計方針どおりに実装されているか、符号・単位・境界 (node の境界半割面、周期 seam、軸)・
   float32 桁落ち・ゼロ割ガードの絶対閾値などの誤りがないかを見る。
2. **検証結果の裏付け**: plan の §6 / 変更ログに書かれた数値・主張 (収束、一致、改善率) を、挙げられている `run_*`
   ディレクトリの `CONVERGENCE_VERDICT.txt` / `residual_history.csv` / `README.md` の run 一覧で確認する。
   主張と実測が食い違う箇所、VERDICT が NOT CONVERGED / DRIFTING のまま「一致」と書いている箇所を挙げよ。
3. **回帰**: 既存機能 (既定値、他ケース) を壊していないか。既定挙動が変わった場合にそれが plan に明記されているか。
4. **文書整合**: `methods/<area>/` の現在仕様、`procedures/` の手順、`methods/index.md`、`plans/README.md` が
   実装と一致しているか。
5. **残作業表**: 未解決事項が §5.1 の残作業表に残っているか (対話で決めて書いていない、が無いか)。

最後に「accepted に移してよいか」を **GO / GO-with-changes / NO-GO** の 1 語で判定し、
GO-with-changes なら移す前に直すべき点を優先順で列挙すること。

## 重点

2 回目 (1 回目 NO-GO C0/M5/m1 への対応確認): (1) M1 の条件統一と #6a を決定論的試験 (面流束・状態のビット一致) で閉じたことの妥当性 (2) #6c の比 (場ごと) と定数場の壁節点除外が事後の緩和でないか (3) #6b 判定不能を accepted の障害にしない扱い (4) R1 の限定合格 + dF1_inf 定常性を継続課題とする扱い、旧延長の連結区間 PASS と F1 切り分け (5) R2 の結論が観測と帰属に限定されているか (6) accepted にしてよいか

## plan 全文 (`plans/active/boundary-node-periodic-gradient-fix.md`)

```markdown
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
```

## 参考: `case/09.Taylor-Green/_g0_lsq_seam/README.md`

```
# G0/G1/G2/G2′: node 並進周期の継ぎ目勾配 (plan boundary-node-periodic-gradient-fix §6)

32³ の三重周期 TGV メッシュ (`case/09.Taylor-Green/mesh`、node 変換) などに線形・二次場を焼き込み、1 step の勾配を
解析値・CPU 参照と比べる。run はすべてスクラッチ (`/tmp/claude-1000/.../scratchpad/g_harness_*`、破棄可)。

## ファイル

| ファイル | 内容 |
| --- | --- |
| `bake_linear.py` / `aggregate.py` | 初版の G0 (Ux=10+y、x 継ぎ目の接線成分のみ) |
| `solverConfig.yaml` / `bcondConfig.yaml` / `solverConfig_sst.yaml` | 初版の設定 (TGV の周期 bcond は以後も共通で使う) |
| `G0_translational.txt` / `G0_translational_m1.txt` | 初版の結果 (修正前 `3de26cba` / 修正後 M1) |
| `gharness.py` | 共通部: h5 読み込み、周期 group の再構成 (面重心 + 並進で相手面を照合、root = 最小 index)、合併 stencil LSQ の double 参照、node GG の float32 再現と double、場の焼き込み、run 作成 (`run_case.sh` を block 128 で起動) |
| `mkmesh.py` | gmsh 構造格子 → 節点ジッタ (±0.2h、index を周期で折り返したハッシュ = 周期像は同じ量) → node 変換 → `check_mesh_quality.py` |
| `g_suite.py <variant>` | G0 (解析値) / G2 (CPU double 参照) / G1-a・G1-b (GG、2026-09-26 #6c から誤差比 ≤ 2 も判定) / 4 量の床 を 1 run で測り、GPU 定数場 (k・ω・ξ = const、別 run `g_harness_<variant>_const`) を判定する。変種は下表 |
| `g0_cpu_cases.py` | G0 の CPU 作用素試験: float32 反例 (root 両順序)、部分 CV の rank 欠損 → 合併で回復、合併後も退化 → 打ち切り参照 |
| `g2p_jitter.py` | G2′: ジッタ格子 16³/32³/64³/128³ の二次場で成分ごと (ゼロ成分含む) の次数・継ぎ目/内部比、一様 32³ の別行 (2026-09-26 #6b で訂正後の基準に揃えた) |
| `r6a_axi_periodic.py` | #6a: 軸対称 × 並進周期の小メッシュ (61×41) で 2 成分 + 受動トレーサ ξ の 1 step を旧新で比較 (§6 R3 と同じ揺れの規則、`r3_compare.py` の関数を流用) |
| `R6a_axi_periodic.txt` | #6a の結果 (末尾に反復 12 本の総当たりのノイズ分布、参考) |
| `f1_read_test.py` | F1 の読み出し (codex m3): 1・2 回目の残差組立で輸送が計算値 F1 を読むか |
| `G_<variant>.txt` / `G0_cpu_cases.txt` / `G2p_jitter.txt` / `G_f1_read.txt` | 結果 (各ファイル末尾に VERDICT) |

場は ρ=1 で焼く (channel は ρ も線形)。w は継ぎ目中心の局所座標 (周期軸を中央で折り返し、group 内は root の値)。
折返し不連続に触れる節点は解析値比較から除く。状態は `res_0.h5`、勾配は `res_1.h5` (k/ω の勾配は res_0 では未計算)。

## 結果 (2026-09-26、バイナリ sha256 `bb55e1bf…`、HEAD `1f03c036`)

| 試験 | 変種 / メッシュ | VERDICT | 主要数値 |
| --- | --- | --- | --- |
| G0 | tgv (一様、2/4/8 member) | PASS | 非零成分 最大誤差/\|a\| ≤ 1.34e-6、ゼロ成分 0、ρ=1 の勾配 0 |
| G0 | tgv_mirror (点対称 = root を継ぎ目の反対側へ) | PASS | 同上。(group, 継ぎ目軸) 3072 組で root が下側にあるのは tgv 2080 / mirror 992 = 全組で側が入れ替わった |
| G0 | tgv_shift (原点 (100,100,100)) | PASS | ≤ 8.1e-7 |
| G0 | jitter32 (非対称 stencil) | PASS | 非零成分 ≤ 1.62e-6、ゼロ成分 ≤ 6.2e-8 |
| G0 | channel (x 等比・y 両壁、壁∩継ぎ目) | PASS | ρ 線形で壁∩継ぎ目 (面・辺) も ≤ 1.5e-6 |
| G0 | CPU: float32 反例 | PASS | 実変位 (現行) は root 両順序で 1.000000。代表変位 (旧) は 0.999980 / 1.000020 |
| G0 | CPU: rank 欠損 → 合併で回復 | PASS | 合併 4.7e-8、部分ごとの和 (旧) は y を 2 重計上 (誤差 0.85) |
| G0 | CPU: 合併後も退化 | PASS | 打ち切り参照との差 4.8e-8 |
| root 順序 | tgv_bcswap / tgv_repeat | PASS | bcond の記述順では root は変わらない (最小 index)。run 間差は atomicAdd 順序の床 ≤ 1.75 ε·max\|φ\|/h |
| G2 | 全 7 run (tgv・mirror・shift・bcswap・repeat・jitter32・channel) | PASS | 最大差/S ≤ 1.4e-7 (折返し不連続を除いた S_clean で ≤ 2.4e-7)。G2′ の 3 水準でも ≤ 2.8e-6 |
| G1-a | 全 7 run | **FAIL** | 継ぎ目で 4.7–27 ε·max\|φ\|/h (内部は ≤ 1.7)。下の「G1-a の不合格」 |
| G1-b | 全 7 run | PASS | 継ぎ目 ≤ 26.5 ε (閾値 2·N_max·n_member = 40)。jitter32 の ω 面は 19.75 / 閾値 20.0 |
| G2′ 次数・比 | ジッタ 16/32/64 | PASS | 継ぎ目の次数 0.95, 0.93、継ぎ目/内部 1.00 / 0.85 / 0.83 |
| G2′ ゼロ成分 | ジッタ 16/32/64 | **FAIL** | ∂ψ/∂z (ψ は z に依らない二次場) が 1.5e-2 → 8.4e-3 → 4.3e-3 (O(h) で収束、床 1.7e-6 を超える)。一様 32³ は 0 |
| F1 読み出し | TGV 32³、F1 計算値 0 / 1 | PASS | 1・2 回目とも Δres_ω/Δres_k = 151.893 / 151.895 (期待 151.893) |
| R3 (2026-09-26、旧 `1266aba1` sha `9b45d018…` / 新 sha `88e949fa…` = §4.2a 修正後) | case48 (fp_y1_12um、SST、一様 IC) / case16 (run_0507 の発達場、speciesFaceReconstruction 1) / axi (101×41 軸対称 node、MIXDRY/H2O、焼き込み場。case/44 のメッシュ削除済みの代替) を旧 2・新 2 で 1 step。`r3_prepare.py` → `r3_compare.py` → `R3.txt`、run はスクラッチ `r3_<case>_<old\|new>_<a\|b>` | PASS (3 ケースとも) | res_1 の勾配・リミタ 30 配列: 旧同士ビット一致の配列は旧新もビット一致。case48 の dK/dΩ (GG) は旧同士でも 62–544 点不一致、旧新は 49–570 点・最大差は旧同士と同値 (ノイズ床)。NaN/Inf なし。化学種 ∇Y は出力に無く直接は比べていない (res_1 の roY は旧同士と同程度の不一致)。case48 の res_1 roK/roOmega は旧新で 5.1 万点不一致 (旧同士 0/2 点) |

### 4 量の床 (GPU − CPU double、ε·max|φ|/h 単位、tgv)

| 量 | 内部 最大 | 継ぎ目 最大 |
| --- | --- | --- |
| dUx (LSQ) | 0.4 | 0.4 |
| dK (GG) | 1.1 | 20.8 |
| dΩ (GG) | 1.3 | 26.5 |
| dξ (GG、化学種 Y の代理) | 0.9 | 22.7 |

化学種 Y は `dY{s}d*` が出力変数に無い (`variables.hpp` の output_cellValNames) ので測れない。受動トレーサ ξ は同じカーネル
(`species_gradient_d`、`speciesTransport_d.cu:1229`) と同じ合併 (`periodicGradientGather`) を通るので代理にした。

### 修正後 (2026-09-26、plan §4.2a、バイナリ sha256 `88e949fa…`)

`G_*.txt` 7 本は修正後の再実行で上書きした (修正前は git の `633997bf`)。G0・G2・G1-a・G1-b すべて PASS。
継ぎ目の床 (tgv、ε·max|φ|/h) dK 1.7 / dΩ 2.4 / dξ 1.5 (修正前 20.8 / 26.5 / 22.7)、継ぎ目/内部比 0.83–1.85。上の表の G1-a 行は修正前。

### G1-a の不合格: 周期半割面が GG から除外されていない (修正前、上記で解消)

`ransTransport_d.cu:43` と `speciesTransport_d.cu:649` の除外条件は `ip >= nNormalPlanes && ic1 < nCells` だが、実行時の
`plane_cells` では周期半割面にも ghost が割り当てられ (`mesh.cpp:443-456`、全 bcond の iBPlanes に ghost を作る)、
`ic1 >= nCells` になる。したがって周期半割面は除外されず、node の境界面規則 (面値 φ[ic0]) で積算されている。

- CPU の float32 再現に「周期半割面も φ[ic0] で積算」を足すと、GPU との差は継ぎ目でも ≤ 1.75 ε·max|φ|/h (全変種、`G_*.txt` の診断列)。
- 合併後に残るのは対の半割面の面ベクトル不一致 (S_a + S_b)·φ/V。tgv で max|S_a+S_b|/|S_a| = 3.7e-6 → 定数場の GG 勾配が継ぎ目で
  3.7e-6·φ/h (≈ 31 ε·φ/h)。**φ に比例** (Δφ ではない) し、float32 演算の丸めではない。初版で「継ぎ目法線成分の ~20 ulp·φ/h の丸め」
  と書いたもの (`G0_translational_m1.txt`) はこれ。
- NS の GG 経路 (`calcGradient_b_d`) は host 側で `bcondKind == "periodic"` を飛ばしている (`calcGradient_d.cu:1091`)。

### G2′ ゼロ成分の不合格について (事実のみ)

ψ = w·(a_x, a_y, 0) + ½κ(w_x² + ½w_y²) は z に依らないが、ジッタ格子の非対称 stencil では LSQ の z 成分に二次項の打ち切り誤差が
乗り、他の成分と同じく O(h) で収束する (1.5e-2 → 8.4e-3 → 4.3e-3、継ぎ目・内部とも)。一様格子 (対称 stencil) では 0。
§6 G2′ の「ゼロ成分は絶対誤差 ≤ e_floor」は一様格子・線形場でしか成り立たない条件で、ジッタ格子の二次場には当てはまらない。
次数・継ぎ目/内部比の判定 (PASS) とは別の行にしてある。

## run (スクラッチ、破棄可)

`/tmp/claude-1000/-home-sano-work-forge/4b0c8643-66fb-4fde-8878-7c8a5104c060/scratchpad/` の下:
`g_harness_{tgv,tgv_mirror,tgv_shift,tgv_repeat,tgv_bcswap,jitter32,channel}`、
`g_harness_g2p_k0.2_N{16,32,64}_j0.2`、`g_harness_g2p_k0.2_N32_j0`、`g_harness_f1_{inf_b0,inf_b1,inf_b1_repeat,tiny_b0,tiny_b1}`、
メッシュは `g_harness_mesh/{box16_j0.2,box32_j0.2,box64_j0.2,box32_j0,box32_j (jitter32 用),box16_j (動作確認のみ),channel}` (各 `quality.txt` は PASS)。

## §5.1 #6a–#6c (2026-09-26、AWS g5・block 128、新 `bd22376d` sha256 `f0ad00c8…` / 旧 `1266aba1` sha256 `f0505fe8…`)

codex result (M1–M3) を受けた再実行。run はすべて AWS のスクラッチ `~/pgrad6_scratch/` (破棄可):
`g_harness_{tgv,tgv_mirror,tgv_shift,tgv_repeat,tgv_bcswap,jitter32,channel}` と各 `_const`、`g_harness_g2p_k0.2_N{16,32,64,128}_j0.2`・`_N32_j0`、
`r6a_prep`・`r6a_axiper_{old,new}_{a..f}`。tgv 系の入力は `lsqseam_m1/Taylor-Green.h5`、jitter32/channel/G2′ 16–64 のメッシュはローカルで作った
`g_harness_mesh/` をコピーして使った (jitter32・channel は AWS で 1 回変換し直された。品質は再判定して PASS)。
変換器は AWS g5 で終了時に GPUassert (exit≠0) を出すが出力 h5 は完全 (既知の罠)。`mkmesh.py`・`r6a_axi_periodic.py` はこの終わり方だけを許す。

### #6c: G1-b の自動緩和を削除、誤差比と GPU 定数場を機械判定 (`G_<variant>.txt`)

- G1-b の閾値は $2N_{max}\varepsilon\max|\phi|/h$ のみ (旧 `g_suite.py` の `max(1e-5·S, …)` と n_member 倍を削除)。
- 誤差比は G1-b 表と同じ量 (全成分の max 誤差、節点ごとの $h$ で正規化) の継ぎ目/内部を場ごとに判定 (≤ 2)。成分ごとの比は参考列。
- **化学種の代理**: 化学種 ∇Y は出力変数に無いので受動トレーサ ξ で代理する。根拠: `passiveGradient_d_wrapper` (`speciesTransport_d.cu:1220-1235`) は
  化学種の `speciesGradient_d_wrapper` (`:675-692`) と同じ `species_gradient_d` カーネルを同じ引数規則 (軸対称なら `A_planar`/`s*_planar`、
  `excludePeriodic = periodicSeamMergeActive(cfg, msh)`、面フラグ `planePeriodic_d`) で呼び、同じ `species_gradient_normalize_d` で割る。
  継ぎ目の合併は `periodicGradientGather` (`periodicNode_d.cu:184-195`) が `dY{s}d*` と受動種勾配の両方を対象にする。

| 変種 | G0 | G2 | G1-a (最大 [ε·max\|φ\|/h]) | G1-b (最大差/閾値) | G1-b 誤差比 (場ごと、最大) | GPU 定数場 継ぎ目最大 [ε\|φ\|/h] k / ω / ξ | root 順序 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| tgv | PASS | PASS | PASS (1.69) | PASS (0.118) | PASS (1.81、ω) | PASS 1.57 / 1.01 / 1.57 | - |
| tgv_mirror | PASS | PASS | PASS (1.69) | PASS (0.118) | PASS (1.85、ω) | PASS 1.57 / 1.01 / 1.57 | - |
| tgv_shift | PASS | PASS | PASS (1.69) | PASS (0.098) | PASS (1.66、ω) | PASS 1.57 / 1.01 / 1.57 | - |
| tgv_repeat | PASS | PASS | PASS (1.69) | PASS (0.118) | PASS (1.81、ω) | PASS 1.57 / 1.01 / 1.57 | PASS (1.69) |
| tgv_bcswap | PASS | PASS | PASS (1.69) | PASS (0.118) | PASS (1.81、ω) | PASS 1.57 / 1.01 / 1.57 | PASS (1.69) |
| jitter32 | PASS | PASS | PASS (1.67) | PASS (0.161) | PASS (1.16、ω) | PASS 2.33 / 3.03 / 2.36 | - |
| channel | PASS | PASS | PASS (0.86) | PASS (0.151) | PASS (1.03、k) | **FAIL** 1.42 / 1.47 / **5.31** | - |

- 成分ごとの比 (参考、判定外) は tgv 系 5 本の ω の x 成分だけが 2 を超える (2.13–2.15)。他は ≤ 2。
- channel の定数場 ξ の 5.31 は**壁の節点**で出ている (壁 ∩ 継ぎ目 96 節点 max 5.31、壁で継ぎ目でない 210 節点も max 5.31、
  壁でない継ぎ目 528 節点 max 0.82、壁でない内部 1155 節点 max 0.99。最大は節点 21 (x 0.778, y 0, z 0、member 2) の y 成分)。
  k・ω は壁の上書き (ピン) に stencil が触れる節点を除外しているので、壁の節点は判定に入っていない (除外 192 / 598)。

### #6a: 軸対称 × 並進周期 (`R6a_axi_periodic.txt`)

61×41 節点 (品質 PASS)、x 並進周期 (41 面、41 節点を合併)、軸・上面 slip、MIXDRY/H2O + ξ、`speciesFaceReconstruction: 1`、乱流なし、1 step。
res_0 は 4 run とも全 88 配列ビット一致。res_1 は §6 R3 の規則 (旧同士一致なら旧新も一致 / 旧同士が揺れるなら最大差 ≤ 2 倍かつ不一致数の桁が同じ) で
**VERDICT FAIL**: roUx (旧同士 85 / 旧新 132、最大差 1.526e-05 で同値) と roY1 (旧同士 6 / 旧新 11、最大差 3.725e-09 で同値) の不一致数の桁が違う。
他の配列 (roY0・roXi・Y0・Y1・Xi・全勾配・リミタを含む) は規則どおり PASS。最大差はすべての配列で旧同士と同値。
参考として旧 6・新 6 本の総当たり (旧同士 15 組・新同士 15 組・旧新 36 組) を同ファイル末尾に置いた (判定なし)。
例: roUx 不一致数 旧同士 77/109/151・新同士 99/126/159・旧新 112/134/181 (min/中央値/max)、roY1 4/11/16・5/11/14・7/12/18、
dXidx・dXidy の旧新 min は 0 / 1。


### #6b: G2′ を訂正後の基準で再判定、ジッタ 128³ を追加 (`G2p_jitter.txt`)

16/32/64 は前回と同じメッシュ (ローカル生成をコピー)、128³ (129³ = 2.1M 節点、品質 PASS) は AWS で生成・変換。κ = 0.2、全水準で全成分が床 (1.7e-6) を超える。
判定は §5.1 #6b の固定どおり (64→128 の継ぎ目の次数が全成分 ≥ 0.9 かつ 128³ で継ぎ目/内部 ≤ 2、32→64 は補助、16→32 は記録のみ)。

| 対 | φ_x | φ_y | φ_z | ψ_x | ψ_y | ψ_z (ゼロ) | 継ぎ目 最小 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 16→32 (記録) 継ぎ目 / 内部 | 0.952 / 0.727 | 0.499 / 0.811 | 0.897 / 0.864 | 0.641 / 0.778 | 0.577 / 0.879 | 0.849 / 0.730 | 0.499 |
| 32→64 (補助) | 0.927 / 0.889 | 1.129 / 0.885 | 0.876 / 0.881 | 0.979 / 0.945 | 1.067 / 0.861 | 0.985 / 0.982 | 0.876 |
| 64→128 (本判定) | 0.936 / 0.885 | **0.855** / 0.976 | **0.869** / 0.939 | 0.969 / 0.825 | 0.911 / 1.000 | 1.009 / 0.961 | **0.855** |

128³ の継ぎ目/内部: 0.77–0.87 (全成分 ≤ 2)。各水準の継ぎ目/内部の最大: 1.00 / 0.95 / 0.85 / 0.87。
**VERDICT: 判定不能** (64→128 の継ぎ目の次数 最小 0.855 < 0.9、φ_y と φ_z)。内部の次数 (64→128) は 0.825–1.000。
G2 (GPU − CPU double)/S は 128³ で 9.9e-6 (16/32/64 は 6.6e-7 / 1.8e-6 / 2.8e-6)。一様格子 N=32 の行 (記録) は床以下、ゼロ成分 0。
```

## 参考: `case/09.Taylor-Green/_g0_lsq_seam/R6a_axi_periodic.txt`

```
#6a: 軸対称 × 並進周期で 1 step 後の VALUE が旧新で不変 (plan boundary-node-periodic-gradient-fix §5.1 #6a、codex result M1)
run root: /home/ubuntu/pgrad6_scratch/r6a_axiper_<old|new>_<a|b>
old_a: forge_sha256: f0505fe8385898a3a5a8438eafb508da9fb67b5fde1045cfcb85950e58c23845 git_head    : bd22376da521335470208b856a4d3ae341adf96b
new_a: forge_sha256: f0ad00c8b3530576f598f81064cb3aefca98560262536a7c7051ae0d64740a39 git_head    : bd22376da521335470208b856a4d3ae341adf96b
mesh quality: VERDICT: PASS (AR<= 1000, skew<= 0.90)
new_a log: [setPeriodicPartner] physID 1 <-> 2: 41 faces, max residual 0 m (face size 0.01581137441 m)
new_a log: [buildPeriodicNodeGroups] node periodic DOF identification: 41 slave CVs merged into masters (union-find).
new_a log: [passive] initial total roXi     1.171355910523e-03 (root-only, before the first physical step)
列: 不一致数 (ビット比較) / 最大絶対差。判定は oo=old_a/old_b と on=old_a/new_a に §6 R3 の規則を当てる。
NaN/Inf: なし (全 run・res_0/res_1 の全 VALUE)

-- res_1.h5 全 VALUE (88 配列)  [VERDICT 対象]
array                 old_a/old_b        new_a/new_b        old_a/new_a  判定
P                    0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
Pk_diag              0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
T                    0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
Taw_diag             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
Ux                   0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
Uy                   0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
Uz                   0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
Xi                   2  5.960e-08       2  5.960e-08       3  5.960e-08  PASS (ノイズ床)
Y0                   1  1.192e-07       0  0.000e+00       2  1.192e-07  PASS (ノイズ床)
Y1                   5  3.725e-09       6  3.725e-09       5  3.725e-09  PASS (ノイズ床)
axisym_divU          0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
cfl                  0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dKdx                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dKdy                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dKdz                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dOmegadx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dOmegady             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dOmegadz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dPdx                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dPdy                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dPdz                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dTdx                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dTdy                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dTdz                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUxdx                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUxdy                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUxdz                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUydx                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUydy                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUydz                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUzdx                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUzdy                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUzdz                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dXidx                9  3.338e-05       9  3.338e-05       9  3.338e-05  PASS (ノイズ床)
dXidy                9  6.866e-05      10  6.056e-05       6  6.866e-05  PASS (ノイズ床)
dXidz                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
delta_les            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
drodx                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
drody                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
drodz                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dt_local             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
ducros               0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
fd_shield            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
fe_iddes             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
gamma                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
h0                   0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
k                    0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
l_des                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
limiter_P            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
limiter_T            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
limiter_Ux           0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
limiter_Uy           0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
limiter_Uz           0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
limiter_Xi          21  5.364e-07      22  5.364e-07      18  5.364e-07  PASS (ノイズ床)
limiter_ro           0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
omega                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
passiveFloorCorr_Xi       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
passiveLimCorr_Xi       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
rd_des               0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
res_ro             460  2.910e-11     411  2.910e-11     510  2.910e-11  PASS (ノイズ床)
res_roK              0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
res_roOmega          0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
res_roUx           500  1.192e-07     463  1.192e-07     607  1.192e-07  PASS (ノイズ床)
res_roUy           639  9.537e-07     406  9.537e-07     540  9.537e-07  PASS (ノイズ床)
res_roe            523  3.052e-05     398  1.526e-05     512  1.526e-05  PASS (ノイズ床)
ro                   3  1.192e-07       7  1.192e-07       8  1.192e-07  PASS (ノイズ床)
roK                  0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
roK_wf               0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
roOmega              0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
roUx                85  1.526e-05     109  1.526e-05     132  1.526e-05  FAIL (不一致数の桁 132 vs 85)
roUy              1154  1.150e-05     836  1.019e-05    1000  9.686e-06  PASS (ノイズ床)
roUz                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
roXi                 5  1.192e-07       9  5.960e-08       9  1.192e-07  PASS (ノイズ床)  <- roY/roXi
roY0                 2  1.192e-07       7  1.192e-07       7  1.192e-07  PASS (ノイズ床)  <- roY/roXi
roY1                 6  3.725e-09      13  3.725e-09      11  3.725e-09  FAIL (不一致数の桁 11 vs 6)  <- roY/roXi
roe                 23  3.125e-02      12  3.125e-02      15  3.125e-02  PASS (ノイズ床)
sonic                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
src_jac_k            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
src_jac_omega        0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
thermCond            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
transport_diag_k       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
transport_diag_omega       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
vis_lam              0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
vis_turb             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
volume               0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
wall_dist            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
wf_irep_flag         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
wf_pk                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS

-- res_0.h5 全 VALUE (88 配列)  [参考: step 前]
array                 old_a/old_b        new_a/new_b        old_a/new_a  判定
P                    0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
Pk_diag              0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
T                    0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
Taw_diag             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
Ux                   0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
Uy                   0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
Uz                   0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
Xi                   0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
Y0                   0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
Y1                   0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
axisym_divU          0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
cfl                  0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dKdx                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dKdy                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dKdz                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dOmegadx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dOmegady             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dOmegadz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dPdx                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dPdy                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dPdz                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dTdx                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dTdy                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dTdz                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUxdx                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUxdy                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUxdz                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUydx                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUydy                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUydz                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUzdx                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUzdy                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dUzdz                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dXidx                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dXidy                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dXidz                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
delta_les            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
drodx                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
drody                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
drodz                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
dt_local             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
ducros               0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
fd_shield            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
fe_iddes             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
gamma                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
h0                   0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
k                    0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
l_des                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
limiter_P            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
limiter_T            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
limiter_Ux           0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
limiter_Uy           0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
limiter_Uz           0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
limiter_Xi           0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
limiter_ro           0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
omega                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
passiveFloorCorr_Xi       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
passiveLimCorr_Xi       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
rd_des               0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
res_ro               0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
res_roK              0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
res_roOmega          0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
res_roUx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
res_roUy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
res_roe              0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
ro                   0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
roK                  0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
roK_wf               0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
roOmega              0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
roUx                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
roUy                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
roUz                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
roXi                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS  <- roY/roXi
roY0                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS  <- roY/roXi
roY1                 0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS  <- roY/roXi
roe                  0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
sonic                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
src_jac_k            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
src_jac_omega        0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
thermCond            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
transport_diag_k       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
transport_diag_omega       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
vis_lam              0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
vis_turb             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
volume               0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
wall_dist            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
wf_irep_flag         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS
wf_pk                0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS

VERDICT #6a (軸対称 × 並進周期、res_1 全 VALUE の旧新): FAIL

## 参考: 反復 run 総当たりの res_1 不一致 (判定なし)。旧 6 本 (old_a,old_b,old_c,old_d,old_e,old_f)・新 6 本 (new_a,new_b,new_c,new_d,new_e,new_f)
列: 不一致数 min/中央値/max、最大絶対差の max。旧同士でも一致しない配列だけ出す (他は全組ビット一致)
array                  旧同士 不一致数 (min/med/max)        最大差         新同士 不一致数 (min/med/max)        最大差          旧新 不一致数 (min/med/max)        最大差
Xi                                      0/4/6  5.960e-08                          2/3/6  5.960e-08                          1/4/6  5.960e-08
Y0                                      0/2/4  1.192e-07                          0/2/4  1.192e-07                          0/2/4  1.192e-07
Y1                                      3/5/7  3.725e-09                          0/4/9  3.725e-09                          1/5/8  3.725e-09
dXidx                                   0/9/9  3.338e-05                          0/0/9  3.338e-05                          0/0/9  3.338e-05
dXidy                                  0/8/10  6.866e-05                         0/4/10  6.866e-05                         1/6/11  6.866e-05
limiter_Xi                            0/12/22  5.364e-07                         0/4/22  5.364e-07                         1/8/23  5.364e-07
res_ro                            330/507/637  2.910e-11                    411/507/660  2.910e-11                    489/620/776  2.910e-11
res_roUx                          400/628/743  2.384e-07                    463/583/746  2.384e-07                    555/691/932  2.384e-07
res_roUy                          464/607/681  9.537e-07                    406/538/758  9.537e-07                    515/667/793  9.537e-07
res_roe                           375/485/589  3.052e-05                    392/472/624  3.052e-05                    492/598/758  3.052e-05
ro                                     2/6/11  1.192e-07                         4/8/12  1.192e-07                         3/8/12  1.192e-07
roUx                               77/109/151  1.526e-05                     99/126/159  1.526e-05                    112/134/181  1.526e-05
roUy                            890/1092/1207  1.150e-05                  834/1054/1268  1.132e-05                  956/1202/1323  1.130e-05
roXi                                   3/9/14  1.192e-07                         6/9/15  1.192e-07                        6/11/15  1.192e-07
roY0                                   2/7/10  1.192e-07                         4/8/12  1.192e-07                         3/8/13  1.192e-07
roY1                                  4/11/16  3.725e-09                        5/11/14  3.725e-09                        7/12/18  3.725e-09
roe                                  15/18/24  3.125e-02                       11/17/23  3.125e-02                       14/20/28  3.125e-02
```

## 参考: `case/09.Taylor-Green/_g0_lsq_seam/R6a_massflux_bitcheck.txt`

```
# #6a 決定論的試験 (plan boundary-node-periodic-gradient-fix §5.1 #6a、2026-09-26 diagnostician 判断)
# 軸対称 × 並進周期 (61×41、x 周期 41 節点合併、MIXDRY/H2O + ξ、speciesFaceReconstruction 1)、同一 res_0 入力から 1 step。
# FORGE_DUMP_MASSFLUX: 初回の面流束 (1 面 1 スレッド・非 atomic) と、そのときカーネルが読んだ状態 (ro,Ux,Uy,Uz,P,sonic)。
# AWS ~/pgrad6_scratch/r6a_axiper_{mfold,mfnew}、mf_{old,new}.bin[.state]
old forge_sha256: f0505fe8385898a3a5a8438eafb508da9fb67b5fde1045cfcb85950e58c23845 (1266aba1)
new forge_sha256: f0ad00c8b3530576f598f81064cb3aefca98560262536a7c7051ae0d64740a39 (bd22376d)
massflux      : 5104 faces   bit-identical  (NaN 0)
state (.state): 15006 values bit-identical  (NaN 0)
diff 1266aba1..bd22376d -- speciesTransport_d.cu: 除外条件の引数 (periodicNodeActive -> periodicSeamMergeActive、nNormalPlanes -> planePeriodic_d) とコメントのみ。
  軸対称では 旧 = 除外条件が死んでいて除外しない / 新 = periodicSeamMergeActive=false で除外しない -> 同じ経路。
VERDICT: PASS (面作用素は同一。R3 規則上の FAIL (roUx・roY1 の不一致数の桁) は場の atomicAdd 集積順序の統計差として、この試験で閉じた)
```

## 参考: `case/09.Taylor-Green/_g0_lsq_seam/G_channel.txt`

```
# channel: plan boundary-node-periodic-gradient-fix §6 G0/G1/G2
run: /home/ubuntu/pgrad6_scratch/g_harness_channel
forge_mtime : 2026-09-26T05:44:54+00:00
forge_sha256: f0ad00c8b3530576f598f81064cb3aefca98560262536a7c7051ae0d64740a39
git_head    : bd22376da521335470208b856a4d3ae341adf96b
harness revision: bd22376 (dirty)
CV 1989, 周期 group 299, member 数の内訳 1:1365, 2:572, 4:52, 周期対の最大残差 6.66e-08, 同値類 (重複 >1) 1176
  注: res_0 の Ux は 306 節点で焼いた値と異なる (境界条件の上書き)。G0 はその stencil を除外
  注: res_0 の Uy は 306 節点で焼いた値と異なる (境界条件の上書き)。G0 はその stencil を除外
  注: res_0 の Uz は 306 節点で焼いた値と異なる (境界条件の上書き)。G0 はその stencil を除外
  注: res_0 の k は 306 節点で焼いた値と異なる (境界条件の上書き)。G0 はその stencil を除外
  注: res_0 の omega は 306 節点で焼いた値と異なる (境界条件の上書き)。G0 はその stencil を除外

## G0: 合併 stencil LSQ の線形場 (解析勾配との差)
閾値: 非零成分 |g_c−a_c| ≤ 1e-5·|a|、ゼロ成分 |g_c| ≤ 1e-6·max|φ|/h (h = 合併体積^(1/3))。対象 = 折返し不連続にも境界条件の上書きにも触れない節点
| 場 | a | 区分 | 節点数 | 非零成分 最大誤差/|a| | ゼロ成分 最大 |g|·h/max|φ| | 判定 |
| --- | --- | --- | --- | --- | --- | --- |
| Ux | (np.float64(0.0), np.float64(1.0), np.float64(0.0)) | 内部 | 585 | 1.43e-06 | 1.38e-07 | ok |
| Ux | (np.float64(0.0), np.float64(1.0), np.float64(0.0)) | 面 (2 member) | 324 | 1.31e-06 | 1.21e-07 | ok |
| Ux | (np.float64(0.0), np.float64(1.0), np.float64(0.0)) | 辺 (4 member) | 36 | 4.77e-07 | 1.07e-07 | ok |
| Uy | (np.float64(1.0), np.float64(0.0), np.float64(0.0)) | 内部 | 585 | 7.15e-07 | 1.11e-07 | ok |
| Uy | (np.float64(1.0), np.float64(0.0), np.float64(0.0)) | 面 (2 member) | 324 | 7.15e-07 | 5.98e-08 | ok |
| Uy | (np.float64(1.0), np.float64(0.0), np.float64(0.0)) | 辺 (4 member) | 36 | 4.77e-07 | 2.67e-08 | ok |
| Uz | (np.float64(0.3), np.float64(-0.7), np.float64(0.5)) | 内部 | 585 | 8.31e-07 | 0.00e+00 | ok |
| Uz | (np.float64(0.3), np.float64(-0.7), np.float64(0.5)) | 面 (2 member) | 324 | 9.03e-07 | 0.00e+00 | ok |
| Uz | (np.float64(0.3), np.float64(-0.7), np.float64(0.5)) | 辺 (4 member) | 36 | 4.38e-07 | 0.00e+00 | ok |
| ρ (線形) | (np.float64(0.8), np.float64(-0.4), np.float64(0.6)) | 内部 | 715 | 1.65e-06 | - | ok (場の丸め ε·max|ρ|/(h·|a|) = 3.8e-06) |
| ρ (線形) | (np.float64(0.8), np.float64(-0.4), np.float64(0.6)) | 内部∩壁 | 130 | 1.47e-06 | - | ok (場の丸め ε·max|ρ|/(h·|a|) = 5.0e-06) |
| ρ (線形) | (np.float64(0.8), np.float64(-0.4), np.float64(0.6)) | 面 (2 member) | 396 | 1.45e-06 | - | ok (場の丸め ε·max|ρ|/(h·|a|) = 3.8e-06) |
| ρ (線形) | (np.float64(0.8), np.float64(-0.4), np.float64(0.6)) | 面 (2 member)∩壁 | 72 | 1.47e-06 | - | ok (場の丸め ε·max|ρ|/(h·|a|) = 5.0e-06) |
| ρ (線形) | (np.float64(0.8), np.float64(-0.4), np.float64(0.6)) | 辺 (4 member) | 44 | 1.34e-06 | - | ok (場の丸め ε·max|ρ|/(h·|a|) = 2.8e-06) |
| ρ (線形) | (np.float64(0.8), np.float64(-0.4), np.float64(0.6)) | 辺 (4 member)∩壁 | 8 | 1.03e-06 | - | ok (場の丸め ε·max|ρ|/(h·|a|) = 3.7e-06) |
VERDICT G0 (channel): PASS

## G2: GPU LSQ と CPU double の合併 stencil LSQ (格納 float32 座標・res_0 の場) の差
閾値: max_i |∇φ_gpu − ∇φ_ref| ≤ 1e-5·S、S = max_i |∇φ_ref| (全節点)
S は折返し不連続の節点で大きくなるので、不連続に触れない節点だけの S_clean と差も併記する (判定は両方)
| 場 | S | 最大差 | 最大差/S | S_clean | 最大差 (clean) / S_clean | 最大差 (ε·max|φ|/h 単位) | 判定 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Ux | 28.62 | 1.70e-06 | 5.94e-08 | 28.62 | 5.94e-08 | 0.6 | ok |
| Uy | 15.49 | 9.20e-07 | 5.94e-08 | 13.92 | 5.28e-08 | 0.6 | ok |
| Uz | 13.89 | 8.51e-07 | 6.13e-08 | 12.81 | 6.01e-08 | 0.6 | ok |
| P | 0 | 0.00e+00 | 0.00e+00 | 0 | 0.00e+00 | 0.0 | ok (参考: 定数場) |
  (打ち切り group 数 0)
VERDICT G2 (channel): PASS

## G1: node GG (k, ω は ransGradient、ξ は受動種 = 化学種と同じ species_gradient_d)
G1-a: GPU − CPU float32 再現 ≤ 4ε·max|φ|/h。G1-b (plan §6 G1、2026-09-26 codex result M3 で自動緩和を削除): GPU − CPU double ≤ 2·N_max·ε·max|φ|/h **かつ** 継ぎ目/内部の誤差比 ≤ 2 (下の「G1-b 誤差比」表。内部が厳密 0 の成分は継ぎ目 ≤ 4ε·max|φ|/h)。h = 合併体積^(1/3) (区分内の最小)、N_max = 区分内の節点の面数の最大。CPU 再現は plan §4.2 どおり周期半割面を除外
診断列 (判定外): 「半割面込み」= 周期半割面も面値 φ[ic0] で積算する CPU float32 再現との差 (GPU の plane_cells では周期半割面の ic1 が ghost なので `excludePeriodic` が効かない仮説の検査)
| 場 | 区分 | 節点数 | N_max | G1-a 最大差 [ε·max|φ|/h] | G1-b 最大差 [ε·max|φ|/h] | G1-b 閾値 [同] | 判定 a / b | 半割面込み再現との差 [同] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| k | 内部 | 1365 | 6 | 0.51 | 1.45 | 12.0 | ok / ok | 0.51 |
| k | 面 (2 member) | 572 | 10 | 0.51 | 1.30 | 20.0 | ok / ok | 4.72 |
| k | 辺 (4 member) | 52 | 16 | 0.34 | 0.76 | 32.0 | ok / ok | 4.85 |
| omega | 内部 | 1365 | 6 | 0.72 | 1.49 | 12.0 | ok / ok | 0.72 |
| omega | 面 (2 member) | 572 | 10 | 0.72 | 1.08 | 20.0 | ok / ok | 5.33 |
| omega | 辺 (4 member) | 52 | 16 | 0.24 | 1.22 | 32.0 | ok / ok | 4.36 |
| Xi | 内部 | 1365 | 6 | 0.86 | 1.81 | 12.0 | ok / ok | 0.86 |
| Xi | 面 (2 member) | 572 | 10 | 0.66 | 1.57 | 20.0 | ok / ok | 4.89 |
| Xi | 辺 (4 member) | 52 | 16 | 0.28 | 0.79 | 32.0 | ok / ok | 4.31 |
  対の周期半割面の面ベクトル不一致 max|S_a+S_b|/|S_a| = 1.08e-06。定数場 φ≡1 の GG 勾配 (半割面込み、double): 継ぎ目 max 6.19e-07·φ/h、内部 max 4.19e-07·φ/h

### G1-b 誤差比 (継ぎ目/内部。判定は場ごと (全成分の max) ≤ 2、内部が厳密 0 なら継ぎ目 ≤ 4ε·max|φ|/h。成分ごとは参考)
| 場 | 成分 | 内部 最大 [ε·max|φ|/h] | 継ぎ目 最大 [同] | 継ぎ目/内部 | 判定 |
| --- | --- | --- | --- | --- | --- |
| k | 全成分 (判定) | 2.33 | 2.40 | 1.03 | ok |
| k | x (参考) | 2.33 | 2.05 | 0.88 | - |
| k | y (参考) | 1.61 | 1.37 | 0.85 | - |
| k | z (参考) | 1.96 | 2.40 | 1.22 | - |
| omega | 全成分 (判定) | 2.34 | 1.76 | 0.75 | ok |
| omega | x (参考) | 2.34 | 1.74 | 0.74 | - |
| omega | y (参考) | 1.32 | 1.16 | 0.88 | - |
| omega | z (参考) | 1.30 | 1.76 | 1.35 | - |
| Xi | 全成分 (判定) | 2.42 | 2.29 | 0.95 | ok |
| Xi | x (参考) | 2.10 | 2.29 | 1.09 | - |
| Xi | y (参考) | 2.42 | 1.70 | 0.70 | - |
| Xi | z (参考) | 1.82 | 1.76 | 0.97 | - |
G1-b 誤差比: ok
VERDICT G1-a (channel): PASS
VERDICT G1-b (channel): PASS (上限 2·N_max・誤差比 ≤ 2 の両方)

## 4 量の床 (GPU − CPU double、ε·max|φ|/h 単位、h = 各節点の合併体積^(1/3))
| 量 | 内部 最大 | 内部 中央値 | 継ぎ目 (member≥2) 最大 | 継ぎ目 中央値 |
| --- | --- | --- | --- | --- |
| dUx (LSQ) | 0.6 | 0.02 | 0.5 | 0.02 |
| dK (GG) | 2.3 | 0.65 | 2.4 | 0.59 |
| dΩ (GG) | 2.3 | 0.53 | 1.8 | 0.55 |
| dξ (GG、Y の代理) | 2.4 | 0.55 | 2.3 | 0.53 |

## G1 定数場 (GPU、k・ω・ξ = const を焼いて 1 step、plan §6 G1・§5.1 #6c)
run: /home/ubuntu/pgrad6_scratch/g_harness_channel_const
forge_mtime : 2026-09-26T05:44:54+00:00
forge_sha256: f0ad00c8b3530576f598f81064cb3aefca98560262536a7c7051ae0d64740a39
git_head    : bd22376da521335470208b856a4d3ae341adf96b
判定: 継ぎ目 (member≥2) で max_c |∂φ/∂x_c| ≤ 4ε|φ|/h。単位 ε|φ|/h (節点ごとの h)。対象 = 境界条件の上書きに stencil が触れない節点 (除外数を併記)
| 場 | φ | 継ぎ目 節点 (除外) | 継ぎ目 最大 [ε|φ|/h] | 内部 節点 (除外) | 内部 最大 [同] (参考) | 判定 |
| --- | --- | --- | --- | --- | --- | --- |
| k | 2 | 432 (192) | 1.42 | 945 (420) | 0.91 | ok |
| omega | 100 | 26 (598) | 1.47 | 62 (1303) | 0.81 | ok |
| Xi | 0.5 | 624 (0) | 5.31 | 1365 (0) | 5.31 | NG |
VERDICT G1 定数場 (channel): FAIL
```

## 参考: `case/09.Taylor-Green/_g0_lsq_seam/G2p_jitter.txt`

```
# G2': 合併 stencil LSQ の二次場の収束 (plan boundary-node-periodic-gradient-fix §6 G2'、§5.1 #6b)
harness revision: bd22376 (dirty)

ジッタ ±0.2h、κ = 0.2、a = (0.3, -0.5, 0.2)。誤差 = 成分ごとの max |∇_gpu − ∇_解析|、e_floor = 10 ε S

## 水準ごと (継ぎ目 / 内部 / 比)
| N | h | 継ぎ目節点 | 内部節点 | S | e_floor | φ_x 継ぎ目 / 内部 (比) | φ_y 継ぎ目 / 内部 (比) | φ_z 継ぎ目 / 内部 (比) | ψ_x 継ぎ目 / 内部 (比) | ψ_y 継ぎ目 / 内部 (比) | ψ_z (ゼロ) 継ぎ目 / 内部 (比) | G2 (GPU−CPU double)/S | run |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 16 | 0.3927 | 1176 | 2196 | 1.334 | 1.6e-06 | 2.612e-02 / 2.620e-02 (1.00) | 1.683e-02 / 2.205e-02 (0.76) | 1.916e-02 / 2.314e-02 (0.83) | 1.869e-02 / 2.372e-02 (0.79) | 1.588e-02 / 2.146e-02 (0.74) | 1.518e-02 / 1.644e-02 (0.92) | 6.6e-07 | g_harness_g2p_k0.2_N16_j0.2 |
| 32 | 0.1963 | 5402 | 24379 | 1.385 | 1.7e-06 | 1.350e-02 / 1.582e-02 (0.85) | 1.191e-02 / 1.256e-02 (0.95) | 1.029e-02 / 1.271e-02 (0.81) | 1.199e-02 / 1.383e-02 (0.87) | 1.065e-02 / 1.167e-02 (0.91) | 8.428e-03 / 9.915e-03 (0.85) | 1.8e-06 | g_harness_g2p_k0.2_N32_j0.2 |
| 64 | 0.0982 | 23066 | 226982 | 1.438 | 1.7e-06 | 7.097e-03 / 8.543e-03 (0.83) | 5.443e-03 / 6.804e-03 (0.80) | 5.605e-03 / 6.904e-03 (0.81) | 6.083e-03 / 7.184e-03 (0.85) | 5.081e-03 / 6.422e-03 (0.79) | 4.257e-03 / 5.021e-03 (0.85) | 2.8e-06 | g_harness_g2p_k0.2_N64_j0.2 |
| 128 | 0.0491 | 95248 | 1953138 | 1.446 | 1.7e-06 | 3.710e-03 / 4.626e-03 (0.80) | 3.010e-03 / 3.458e-03 (0.87) | 3.069e-03 / 3.601e-03 (0.85) | 3.107e-03 / 4.054e-03 (0.77) | 2.703e-03 / 3.211e-03 (0.84) | 2.116e-03 / 2.580e-03 (0.82) | 9.9e-06 | g_harness_g2p_k0.2_N128_j0.2 |
メッシュ品質: N=16: VERDICT: PASS (AR<= 1000, skew<= 0.90); N=32: VERDICT: PASS (AR<= 1000, skew<= 0.90); N=64: VERDICT: PASS (AR<= 1000, skew<= 0.90); N=128: VERDICT: PASS (AR<= 1000, skew<= 0.90)

## 次数 (継ぎ目 / 内部)。16→32 は記録のみ、32→64 は補助、64→128 が本判定
| 対 | 扱い | φ_x | φ_y | φ_z | ψ_x | ψ_y | ψ_z (ゼロ) | 継ぎ目 最小 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 16→32 | 記録のみ | 0.952 / 0.727 | 0.499 / 0.811 | 0.897 / 0.864 | 0.641 / 0.778 | 0.577 / 0.879 | 0.849 / 0.730 | 0.499 |
| 32→64 | 補助 | 0.927 / 0.889 | 1.129 / 0.885 | 0.876 / 0.881 | 0.979 / 0.945 | 1.067 / 0.861 | 0.985 / 0.982 | 0.876 |
| 64→128 | 本判定 | 0.936 / 0.885 | 0.855 / 0.976 | 0.869 / 0.939 | 0.969 / 0.825 | 0.911 / 1.000 | 1.009 / 0.961 | 0.855 |

補助 32→64: 継ぎ目の次数 最小 0.876 (< 0.9)。各水準の継ぎ目/内部 最大: N=16 1.00, N=32 0.95, N=64 0.85, N=128 0.87

一様格子 (ジッタ 0、N=32、記録): 継ぎ目 e = 5.36e-07、内部 e = 6.22e-07、e_floor = 1.6e-06 → 床以下 (対称 stencil で二次場の勾配は厳密)。ゼロ成分 0.0e+00 (≤ e_floor)。run g_harness_g2p_k0.2_N32_j0

VERDICT G2' (ジッタ、64→128 の継ぎ目の次数 全成分 ≥ 0.9 かつ 128³ の継ぎ目/内部 ≤ 2): 判定不能 (64→128 の継ぎ目の次数 最小 0.855 < 0.9)。内部の次数: φ_x 0.885, φ_y 0.976, φ_z 0.939, ψ_x 0.825, ψ_y 1.000, ψ_z (ゼロ) 0.961。128³ の継ぎ目/内部: φ_x 0.80, φ_y 0.87, φ_z 0.85, ψ_x 0.77, ψ_y 0.84, ψ_z (ゼロ) 0.82
VERDICT G2' (一様格子 N=32、ゼロ成分 ≤ e_floor、記録): PASS
```

## 参考: `case/09.Taylor-Green/r2_ke_budget.txt`

```
old: |r| max 0.0990 (t=14.35)、|r_eps| max 0.1019、max K/K0 1.01872、t≤7 の r: 平均 +0.0048 / 正の割合 0.75 / min -0.0170 max +0.0118
new: |r| max 0.0885 (t=13.30)、|r_eps| max 0.0943、max K/K0 1.00000、t≤7 の r: 平均 -0.0209 / 正の割合 0.00 / min -0.0447 max -0.0074
診断 (a) 新 max r = -0.0074 ≤ +0.005: 成立
診断 (b) t≤7 の min(r_old − r_new) = +0.0176 > 0: 成立
診断 (c) 旧 K/K0>1 区間の r>0 最長連続 = 15 snapshot ≥ 10: 成立
DIAGNOSTIC (測定後に追加、保存性ゲートの代替ではない): ALL HOLD
(撤回した旧基準の値、記録のみ: 新 max|r| = 0.0885)
```

## 参考: `case/39.periodic_hills/run_0038_r1_gradfix_old_ext/CONVERGENCE_VERDICT_concat.txt`

```
[segment] 判定区間 = S1_spin -> S2_dev -> S3_ext  (4806000 行) -> run_0038_r1_gradfix_old_ext/residual_history_segment.csv

=== run_0038_r1_gradfix_old_ext  [last step 1601999]  -> PASS (converged) ===
  rms_ro      : init=4.32e-11 peak=1.42e-05 fin=3.83e-10 drop= 4.6dec flat   
  rms_roUx    : init=1.38e-05 peak=2.20e-03 fin=4.07e-08 drop= 4.7dec flat   
  rms_roUy    : init=2.05e-08 peak=2.19e-03 fin=4.59e-08 drop= 4.7dec flat   
  rms_roUz    : init=7.32e-16 peak=1.53e-04 fin=2.84e-08 drop= 3.7dec flat   
  rms_roe     : init=1.12e-05 peak=4.32e+00 fin=1.03e-04 drop= 4.6dec flat   
  rms_roK     : init=5.10e-05 peak=1.30e-04 fin=3.06e-09 drop= 4.6dec flat   
  rms_roOmega : init=3.23e+03 fin=1.29e-03 drop= 6.4dec flat   

OVERALL: ALL PASS
```

## 参考: `case/39.periodic_hills/run_0039_r1_gradfix_new_ext/CONVERGENCE_VERDICT_concat.txt`

```
[segment] 判定区間 = S1_spin -> S2_dev -> S3_ext  (4806000 行) -> run_0039_r1_gradfix_new_ext/residual_history_segment.csv

=== run_0039_r1_gradfix_new_ext  [last step 1601999]  -> PASS (converged) ===
  rms_ro      : init=4.32e-11 peak=1.42e-05 fin=3.84e-10 drop= 4.6dec flat   
  rms_roUx    : init=1.38e-05 peak=2.20e-03 fin=4.09e-08 drop= 4.7dec flat   
  rms_roUy    : init=2.05e-08 peak=2.19e-03 fin=4.62e-08 drop= 4.7dec flat   
  rms_roUz    : init=7.12e-16 peak=1.53e-04 fin=2.89e-08 drop= 3.7dec flat   
  rms_roe     : init=1.12e-05 peak=4.32e+00 fin=1.03e-04 drop= 4.6dec flat   
  rms_roK     : init=5.10e-05 peak=5.79e-05 fin=2.67e-09 drop= 4.3dec flat   
  rms_roOmega : init=3.23e+03 fin=1.22e-03 drop= 6.4dec flat   

OVERALL: ALL PASS
```

## 参考: `case/39.periodic_hills/_f1split/F1_SPLIT.txt`

```
#6d F1 仮説の切り分け (plan boundary-node-periodic-gradient-fix §5.1 #6d、codex result M5)
入力: run_0038_r1_gradfix_old_ext/ic_S3_ext.h5 (= run_0036_r1_gradfix_old/res_800000.h5 の保存量 7 配列とビット一致、input_check.txt。sstF1 は入力に無い)
設定: run_0038 の solverConfig.yaml に nStepOuter 1・outStepInterval 1・turbulence.sstSigmaBlend を明記しただけ。block 128。
AWS run: ~/forge-pgrad-new/case/39.periodic_hills/_f1split/{old_b1,old_b0,new_b1} (スクラッチ)
参照: run_0038 延長の step 0 rms_roK = 1.344498855e-06、--from-floor の参照床 (run_0036 末尾 20 % 平均) 3.29e-09

| run | バイナリ・設定 | step 0 rms_roK | ÷ 参照床 | step 0 rms_roOmega | step 0 rms_ro |
| --- | --- | --- | --- | --- | --- |
| old_b1 | 旧 1266aba1 (sha256 f0505fe8…)、sstSigmaBlend 1 (既定) | 1.344497832e-06 | 408.7 | 1.441476e-03 | 3.874670e-10 |
| old_b0 | 旧 1266aba1、sstSigmaBlend 0 | 1.344499765e-06 | 408.7 | 1.438797e-03 | 3.874867e-10 |
| new_b1 | 新 bd22376d (sha256 f0ad00c8…)、sstSigmaBlend 1 (参考) | 4.996487769e-05 | 15186.9 | 3.076557e-01 | 3.874364e-10 |

判定規則 (§5.1 #6d、測る前に固定): スパイクが sstSigmaBlend 1 側だけ → 仮説支持、両方またはどちらも無し → 原因未特定
スパイク (参照床の 10 倍超を目安に表示): blend 1 あり / blend 0 あり
VERDICT F1 仮説: 原因未特定
```

## 参考: `notes/reviews/2026-09-26-boundary-node-periodic-gradient-fix-result.md`

```
# codex レビュー: boundary-node-periodic-gradient-fix (result)

- **plan**: [`plans/active/boundary-node-periodic-gradient-fix.md`](../../plans/active/boundary-node-periodic-gradient-fix.md)
- **stage**: `result` (diff base `1266aba1`)
- **date**: 2026-09-26
- **commit**: `0977ca72` (feature/sern-design)
- **codex**: effort `high`, 5.6 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m1
- **focus**: (1) dF1_inf の限定合格 (DRIFTING 0.3%/tail、漸近 0.0335、末尾全点 ≤0.05) の妥当性 (2) R2 の KE 収支基準を測定後に符号判定へ訂正したことと、それによる『旧は継ぎ目由来の注入』という帰属の妥当性 (3) R1 抽出の [解釈] 4 点と x 継ぎ目の参考列 (4) §4.2a 面フラグ planePeriodic_d の実装と G1-a/G1-b・R3 の裏付け (5) 旧延長 run の step0 rms_roK スパイクを F1 上書きの仮説とする記述 (6) accepted にしてよいか
- **extra**: `case/09.Taylor-Green/_g0_lsq_seam/README.md`, `case/09.Taylor-Green/_g0_lsq_seam/R3.txt`, `case/09.Taylor-Green/_g0_lsq_seam/G_tgv.txt`, `case/09.Taylor-Green/r2_conservation.txt`, `case/09.Taylor-Green/r2_ke_budget.txt`, `case/39.periodic_hills/run_0039_r1_gradfix_new_ext/QUASISTEADY_1p6M.txt`, `case/39.periodic_hills/run_0038_r1_gradfix_old_ext/CONVERGENCE_VERDICT_fromfloor.txt`, `case/39.periodic_hills/r1_extract.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

合併 LSQ と SST 勾配の gather は支持します。ただし、適用範囲外の軸対称に影響する実装漏れと、訂正後も未達の検証があります。  
R2 の保存性は裏付けられましたが、「注入が消え、数値散逸のみになった」という断定は現在の収支評価では証明できません。

1. **Major — 化学種・受動種の半割面除外が、設計した適用範囲を越えています。**

   **根拠:** [speciesTransport_d.cu:690](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/speciesTransport_d.cu:690) と同ファイル1232行は、`excludePeriodic` に `periodicNodeActive` を渡しています。plan §4.2a が要求する `periodicSeamMergeActive` ではありません。一方、[periodicNode_d.cu:174](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/periodicNode_d.cu:174) は軸対称の勾配 gather を実行しません。

   したがって、**node・軸対称・周期・スカラー再構成有効**では、今回から周期半割面だけが除外され、反対側の寄与は補充されません。定数場でも、除去した面の寄与に対応する `−φ S_periodic/A_planar` が従来勾配との差として生じます。「軸対称の既存経路は変えない」に反します。

   R3 の軸対称ケースは [r3_prepare.py:107](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/r3_prepare.py:107) のとおり入口・出口・壁・軸で、周期境界がありません。既存3ケースの `res_1.h5` を再比較すると **R3: PASS** でしたが、この条件漏れは検出できません。

   **対案:** 両呼出しを `periodicSeamMergeActive` に統一し、軸対称×周期の定数・線形スカラーで旧経路不変を確認してください。回転周期も今回の除外変更から外します。

2. **Major — G2′ のゼロ成分は、訂正後の基準でも一部不合格です。**

   **根拠:** [plan:134](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:134) は、ジッタ格子のゼロ成分にも次数判定を適用しています。保存済み `g_harness_g2p_k0.2_N{16,32,64}_j0.2/res_1.h5` を、ハーネスと同じ抽出条件で再集計しました。

   | 水準 | ゼロ成分の継ぎ目最大誤差 | 収束次数 |
   |---|---:|---:|
   | 16³ | 0.0151796062 | — |
   | 32³ | 0.0084277764 | **0.84891** |
   | 64³ | 0.0042566396 | 0.98544 |

   16³→32³ は要求の **0.9未満**です。誤差は十分に丸め床を超えています。[g2p_jitter.py:99](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/g2p_jitter.py:99) の次数計算は主場だけを対象にし、ゼロ成分には現在も旧基準の「床以下」を適用しています。

   **対案:** ハーネスを訂正後の基準と一致させ、この結果を未達として §5.1 に戻してください。漸近領域への到達を追加水準で確認する検証計画を先に定め、現行結果を保存したうえで再判定すべきです。これは直ちに LSQ 実装の誤りを意味しませんが、「検証完了」は成立しません。

3. **Major — G1 の PASS は、現行仕様の全条件を検査した PASS ではありません。**

   **根拠:** [g_suite.py:265](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/g_suite.py:265) は、G1-b の閾値を `max(1e-5*S, 2*Nmax*εφ/h)` とし、超過時にはさらに `n_member` 倍へ広げます。plan の現行基準には、この自動緩和はありません。また、継ぎ目／内部の誤差比は表示だけで、PASS 判定に入りません。

   保存 HDF5 から7ケースの誤差比を再計算すると、いずれも2以内でした。したがって半割面修正の効果自体は支持できます。しかし、[plan:163](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:163) が明記する **GPU 定数スカラー場の未実施**は残っています。ハーネスの定数場計算は CPU double・半割面込みの診断であり、修正後 GPU の定数保存試験ではありません。

   **対案:** 自動緩和を除き、誤差比と定数場を機械判定に組み込んでください。§5.1 #5a/#5c にある直接の化学種勾配確認も、実施または代理試験で代替する根拠を明示して閉じる必要があります。

4. **Major — R2 の符号判定は再現できますが、「収支の残りは数値散逸のみ」という帰属は過剰です。**

   **根拠:** [r2_ke_budget.py:50](/home/sano/work/forge-sern-design/case/09.Taylor-Green/r2_ke_budget.py:50) は空間微分を独自の周期中心差分で評価し、76行では snapshot 間差分から `dK/dt` を作っています。これはソルバが実際に用いた面応力・圧力仕事・対流流束の離散収支ではありません。

   CSV から次は確認できました。

   - 新の `max r = −0.0074228083`
   - `t≤7` の `min(r_old−r_new) = +0.0175604664`
   - 旧の `K/K0>1` かつ `r>0` は15 snapshot連続
   - 旧の `max K/K0 = 1.0187202`

   しかし `r` には、数値散逸に加えて**空間作用素の不一致、snapshot 時間微分誤差**が含まれます。負の総和から、正の寄与が消失したとは断定できません。旧新差を今回の修正に帰属することと、その差を「継ぎ目でのエネルギー注入」と機構まで確定することは別です。

   さらに [r2_ke_budget.py:105](/home/sano/work/forge-sern-design/case/09.Taylor-Green/r2_ke_budget.py:105) と保存された結果は、今も **`VERDICT (新 |r| ≤ 0.05): FAIL`** です。訂正後の符号条件は判定ツールに反映されていません。

   **対案:** 今回は「修正後に KE 超過が消え、この後処理での収支残差が負になった」という観測までに結論を限定してください。符号条件を正式な保存性ゲートの代替にせず、測定後に追加した診断として実装・記録を同期します。機構の確定は、同一状態に対する旧新の離散残差から KE 仕事を分解する後続検証へ残すのが妥当です。

5. **Major — R1 の限定合格には根拠がありますが、「他の基準はすべて成立」は誤りです。**

   **根拠:** 再実行した `check_quasisteady.py` は、保存記録と同じ結果でした。

   - `case/39.periodic_hills/run_0038_r1_gradfix_old_ext/`：**ALL STEADY**
   - `case/39.periodic_hills/run_0039_r1_gradfix_new_ext/`：**dF1_inf = DRIFTING、NOT ALL STEADY**

   新の末尾64点の `dF1_inf` は **0.03334079～0.03344835**。全点が0.05未満なので、**差の大きさに対する限定合格は支持します**。ただし増加33回・減少30回で、単調ではありません。漸近値0.0335を再現する推定手順も提示されていません。

   また旧延長の [CONVERGENCE_VERDICT_fromfloor.txt:3](/home/sano/work/forge-sern-design/case/39.periodic_hills/run_0038_r1_gradfix_old_ext/CONVERGENCE_VERDICT_fromfloor.txt:3) は **NOT CONVERGED**。`rms_roK` のピーク408.52倍を、末尾0.94倍だけで合格扱いにはできません。[plan:215](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:215) の完了記述と矛盾します。

   F1 初回上書きは旧コードに実在し、スパイク原因の**仮説としては妥当**です。ただし、この restart での原因確認や PASS の代わりにはなりません。

   **対案:** 「単調収束中・漸近値」を「微小な増加傾向と変動が残る、観測窓内で上限成立」に訂正し、未達を §5.1 に残してください。旧側は元 run と延長の同一実効設定区間を正式に連結して判定し、スパイクを含む全履歴を保存します。F1 仮説は同一 restart 入力の短い切り分けで確認してください。新の定常性確認も継続課題として明示します。

6. **Minor — 現在仕様・スコープ・残作業表の同期が不足しています。**

   **根拠:** [methods/gradient.md:156](/home/sano/work/forge-sern-design/methods/gradient.md:156) は現在も「LSQ が2倍」「SST未合算」を現状として記述しています。[plan:47](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:47) は化学種・受動種 GG を対象外としていますが、実装では変更済みです。[plans/README.md:27](/home/sano/work/forge-sern-design/plans/README.md:27) も `draft` のままです。

   §5.1 #5d は古い R2 判定式・閾値を保持しつつ「完了」。G2′ の新基準での未達、旧延長の収束未達、F1 仮説の検証が独立した残作業になっていません。

   **対案:** 現在仕様を実装後の記述へ更新し、過去の欠陥・撤回した判断は履歴として分離してください。上記未解決事項を §5.1 に登録し、case の run 一覧も最終状態へ同期します。

補足すると、R1 抽出の第一内部層との比較、x継ぎ目除外、勾配ノルムの L2 比は、z継ぎ目の診断として妥当です。F1 再計算も今回のモデルの式・定数と対応しています。ただし実際に輸送が読んだ F1 の直接測定ではありません。x継ぎ目列も、今回の格子での参考指標という扱いを支持します。これらを局所勾配全成分の一致証明へ広げるべきではありません。

R2 の保存 CSV は全4本で既定の保存上限を満たしました。`check_convergence.py` の **NOT CONVERGED** は非定常 TGV では想定どおりです。一方、R1 最終 run の残差 CSV と本番 HDF5 はローカルに無いため、収束については保存 VERDICT の照合までで、元履歴からの再判定はできていません。

**推奨は、`active` に留め、1→2・3→5→4→6の順に修正・再判定してから再レビューすることです。** 中核実装を撤回する必要はありませんが、現状での `accepted` 移行は認められません。ファイルは変更しておらず、指摘は plan 未反映です。

指摘数: Critical 0 / Major 5 / Minor 1
```

## 出力形式

1. 冒頭に **判定 (GO / GO-with-changes / NO-GO)** と 3 行以内の要約。
2. 指摘一覧 (Critical → Major → Minor の順、番号付き。各項目に根拠と対案)。
3. 推奨 (1 つに絞る)。
4. 末尾に `指摘数: Critical N / Major N / Minor N` の 1 行。
