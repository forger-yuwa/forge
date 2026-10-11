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

3 回目 (2 回目 NO-GO C0/M4/m1 への対応確認)。§8 に固定した accepted 条件 1–6 に対する採点表 (§6.1『§8 採点』行、§6.2 #8a–#8d) の妥当性を採点してほしい。特に (1) #8a 面寄与ダンプの 3 バイナリ比較と、並進周期の基準を測定前に訂正したこと (2) #8b の定義補正 2 点 (角 group 3 格子×1、丸め床) が事後緩和でないか、H' 追試 (3) #8c 定数場 3 区分 (4) #8d の帰結 (旧の restart スパイク = F1 初回上書き) の書き方 (5) 継続課題 (dF1_inf・壁閉包 5ε・軸対称×周期スカラー未合併・ジッタ格子の LSQ 次数・R2 機構分解) の追跡先。§6.1 の運用規則: §8 の外側の新しい要求で中核の符号・適用条件・単位の誤りを示さないものは後続 plan 項目になる。accepted にしてよいか

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
| 5d (**完了 2026-09-26**: R2 収支 成立、R1 延長 限定合格) | R1 延長・R2 収支 (判断: 2026-09-26 `diagnostician`) | R1: 800k 場から `restart_field.py` で +800k step (同一設定・同一閾値、10k–1.6M 全系列 `--tail 0.4`)。なお DRIFTING なら「限定合格 (~~単調収束中・漸近値・~~末尾全点で ≤ 0.05)」と書きラベルは変えない (取り消し線は codex result M5 で撤回した表現)。F1 差の L2 と x 継ぎ目の比を参考列に。case/39 の `slauWallNormalChi` 実効値を旧新で控える。R2: 既存 snapshot から $r(t)=(dK/dt+\varepsilon)/(K_0/t_c)$ (測る前に固定したのは ~~新は $|r|\le0.05$、旧は $t\lesssim7$ で $r>0$ が持続すれば「継ぎ目由来の注入」~~ — 撤回: 前者は数値散逸の見落とし、後者は機構の断定 (codex result M4)。現行の扱いは §6.2 R2: 符号診断 (a)(b)(c) を測定後の診断として記録し、結論は観測と帰属に限定) | O |
| 6a (2026-09-26 実施、massflux のビット一致は検出力不足で #8a へ) | codex result M1 | `speciesTransport_d.cu:690,1232` の `excludePeriodic` を `periodicSeamMergeActive` に統一。合格: 軸対称×周期の小メッシュで旧新の 1 step 後 roY/roXi がビット同一。§7 に「軸対称×周期・回転周期のスカラー勾配は片側 GG + 半割面込みのまま (既存の未修正挙動)」 | O |
| 6b (**完了 2026-09-26**、§6.2) | codex result M2 | `g2p_jitter.py` を訂正後の基準に揃え、ジッタ 128³ を追加 (AWS)。**測る前に固定**: 最細対 64→128 の継ぎ目の次数 (ゼロ成分含む全成分) ≥ 0.9 かつ 128³ で継ぎ目/内部 ≤ 2。32→64 は補助、16→32 は漸近域外として記録のみ。64→128 が 0.9 未満なら判定不能とし、内部の次数を並べて継ぎ目固有でないことだけ書く | O |
| 6c (**完了 2026-09-26**、§6.2) | codex result M3 | `g_suite.py` の G1-b から n_member 倍の自動緩和を削除 (上限 $2N_{max}\varepsilon\phi/h$ かつ継ぎ目/内部 ≤ 2 を機械判定)、GPU 定数場 (k・ω・ξ = const、継ぎ目勾配 ≤ 4ε\|φ\|/h) を追加、7 run 再実行 (AWS)。化学種の直接確認は tracer ξ で代理 (根拠: `passiveGradient_d_wrapper` は同じ `species_gradient_d` を同じ excludePeriodic で呼ぶ、`speciesTransport_d.cu:1220-1232`) | O |
| 6d (連結区間は完了、F1 切り分けは #8d へ) | codex result M5 | R1 の記述訂正。旧は元 run + 延長を `stage_manifest` で同一実効設定の 1 区間として連結し `--segment` (スパイク込み)。F1 仮説: 同じ restart 入力で旧バイナリ 1 step を `sstSigmaBlend` 0/1 で比較 (スパイクが 1 側だけなら仮説支持、そうでなければ「原因未特定」)。新 dF1_inf の定常性は継続課題 (#7) | O |
| 6e (**完了 2026-09-26**: 判定行を置換、`r2_ke_budget.txt` 更新、診断 ALL HOLD) | codex result M4 | R2 の結論を観測と帰属に限定し、`r2_ke_budget.py` の VERDICT 行を符号条件 (a)(b)(c) に置換、旧 \|r\| ≤ 0.05 行は履歴 | O |
| 6f (**完了 2026-09-26**: `methods/gradient.md`・§2・`plans/README.md`) | codex result m1 | `methods/gradient.md` を実装後の記述に (欠陥は履歴節)、§2 スコープ、`plans/README.md`、§5.1 #5d、case/09・case/39 の run 表 | O |
| 8a (**完了 2026-09-26**、§6.2) | result-2 M1 | 診断ダンプ `FORGE_DUMP_SCALARGRAD=<path>` (既定 off): `species_gradient_d` の初回呼び出しで atomicAdd に渡す同じ値を面ごとの別配列へ非 atomic に書く (ポインタ null で分岐、化学種・受動種の両呼び出し)。R3 (非周期・既定 off) で全場ビット同一を再確認。**3 バイナリ**: 新+ダンプ / 旧 1266aba1+ダンプ / 条件を `periodicNodeActive` に戻した誤り版+ダンプ。**合格 (測る前に固定)**: ~~軸対称×周期と並進周期の両方で 旧+ダンプ = 新+ダンプ が全面ビット同一~~ → **測定前の訂正 (2026-09-26、親)**: 並進周期では旧は除外条件が死んでいて半割面を積算し、新は除外するので、全面一致は誤った期待。**軸対称×周期**: 旧+ダンプ = 新+ダンプ が全面ビット同一。**並進周期**: 周期半割面以外の全面で旧新ビット同一、周期半割面は新で 0・旧で非 0。かつ誤り版は軸対称×周期で周期面の寄与が 0 になり新と差が出る (検出力) | O (判断: 2026-09-26 `diagnostician`、診断ダンプは条件 6 の先例 `FORGE_DUMP_MASSFLUX` と同じ扱い) |
| 8b (**完了 2026-09-26**、§6.2) | result-2 M2 | G2 参照 (`lsq_seam_ref`、GPU と ≤ 1e-5·S 一致済み) をそのまま使い (新規再実装しない)、32³ ジッタ格子の実在の継ぎ目 group (面 2・辺 4・角 8 member、各 ≥ 2、計 ≥ 6) を無次元形状固定で $h,h/2,h/4$ に縮小、二次場で**全成分の次数 ≥ 0.9**。誤差/h の定数を同形状の内部 stencil と比較して表に。注記: 形状固定では次数はほぼ 1 になり、示すのは合併作用素の 1 次整合。ジッタ系列 (0.855) は判定不能 (診断) として保存 | O |
| 8c (**完了 2026-09-26**、§6.2) | result-2 M3 | 定数場を 3 区分で機械判定: 壁なし継ぎ目 ≤ 4ε / 壁∩継ぎ目は格納 float32 面ベクトルの閉包 $\Sigma S_f$ (double 和) から作る期待値 $\phi\Sigma S_f/V$ との差 ≤ 4ε / 非継ぎ目壁は同じ参照で報告。壁の絶対定数保存 (5.3ε) は既知制約として §7 と `methods/gradient.md` に | O |
| 8d (**完了 2026-09-26**、§6.2) | result-2 M4 | `1266aba1` に m3 だけを当てたバイナリ (AWS)、同じ restart 入力 `ic_S3_ext.h5` で 1 step。**規則 (測る前に固定)**: m3 のみ版で rms_roK/床 < 10 → 「F1 初回上書きが原因」、≥ 10 → 「原因未特定 (F1 でない)」 | O |
| 8e (**完了 2026-09-26**: #5d・#7・case/39 README・`methods/gradient.md`) | result-2 m1 | §5.1 #5d の旧記述を履歴化、case/39 README の 800k 行と延長行 (連結区間 PASS)、#7 に追跡先 plan と完了条件 | O |
| 6g (2 回目は 2026-09-26 実施、NO-GO → #8a–#8e) | codex result 2 回目 | 6a–6f の後。`--focus` 「M1–M5 の閉じ方、R1 の限定合格 + 継続課題の扱い、R2 の観測限定の文言」 | O (結論 F) |
| 7 | 継続課題 (accepted 後も残す。追跡先と完了条件) | (1) 新 R1 の dF1_inf の定常性 (DRIFTING 0.3 %/tail、上限内) → 追跡先: 本 plan の §7 (accepted 後は後続の LSQ 統一 plan で R1 を再実行するときに同じ抽出で再判定)。完了条件: 同じ閾値で STEADY。(2) R2 の機構分解 → 追跡先: 後続の LSQ 統一 plan §6。完了条件: 同一状態の旧新の離散残差から KE 仕事を分解し、差の所在 (粘性応力/対流再構成) を示す。(3) 壁半割面込み GG の定数場閉包 (≈5ε) → 追跡先: 後続の LSQ 統一 plan (LSQ 化で定数場は厳密)。完了条件: 壁節点の定数場勾配 ≤ 4ε。(4) ジッタ格子での LSQ 勾配次数 (0.83–1.0、内部も同じ) → 追跡先: 後続の LSQ 統一 plan。(5) 軸対称×周期・回転周期のスカラー勾配の合併 → 追跡先: `boundary-node-rotational-periodic` | O |
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
| result (2 回目) | `2026-09-26` | [2026-09-26-boundary-node-periodic-gradient-fix-result-2.md](../../notes/reviews/2026-09-26-boundary-node-periodic-gradient-fix-result-2.md) | **NO-GO**, C0/M4/m1 (中核実装・条件統一・R1/R2 の扱いは支持) | **全件採用** (2026-09-26 `diagnostician`)。M1 → 面寄与ダンプで 3 バイナリ比較 (#8a)。M2 → G2 参照による固定形状縮小試験 (#8b)。M3 → 定数場 3 区分の機械判定 (#8c)、前回の「壁節点を k/ω と同じく除外」は根拠が誤りで撤回。M4 → m3 のみ版で F1 仮説を識別 (#8d)、前回の「F1 仮説は不支持」は識別能力の無い A/B に基づいていたので撤回。m1 → 同期 (#8e)。**accepted の条件を §8 に固定**。**3 回目の運用規則**: Major が §8 の 1–6 の外側 (新しい要求) にあり、中核の符号・適用条件・単位の誤りを示さないなら後続 plan の項目として登録して accepted へ進む。中核の誤りを示す指摘は採用して戻る |
| §8 採点 | `2026-09-26` | (本 plan §8、§6.2 #8a–#8d) | `diagnostician`: §8 の 1–6 はすべて成立 (#8b は定義補正 2 点つき、$H'$ 追試で角 z 成分も次数 1.0) | 採用。codex result 3 回目へ |

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
    ~~(ii) GPU 定数場の「継ぎ目」区分から**壁節点を除外** (k/ω と同じ扱い)~~ (撤回、codex result-2 M3: k/ω の除外は BC で壊れた stencil で、壁だからではない。#8c で 3 区分判定に置換)。channel の ξ は壁∩継ぎ目 96 点・非継ぎ目壁 210 点とも最大 5.31 ε|φ|/h で同値 = 壁半割面込み GG の既存の閉包床で継ぎ目由来でない。壁でない継ぎ目は 0.82 → PASS。
    ハーネスの `G_channel.txt` の VERDICT 行は壁節点込みの FAIL のまま (判定は本節の区分による)。定数場は他 6 変種 PASS (最大 3.03)。
  - **#6b G2′** (`G2p_jitter.txt`、ジッタ 16/32/64/128³): 継ぎ目の次数 (全成分最小) 16→32 0.499 (記録のみ)、32→64 0.876、**64→128 0.855** → 登録どおり**判定不能**。
    内部の次数 0.825–1.000、128³ の継ぎ目/内部 0.77–0.87 → 継ぎ目は内部より悪くない (次数が 0.9 を跨ぐのはジッタ格子での LSQ 自体の性質)。G2 は 128³ で 9.9e-6·S (閾値 1e-5、余裕 1 %)。accepted の障害にしない。
  - **#6d R1**: 連結区間 (S1→S2→S3_ext を 1 区間、`CONVERGENCE_VERDICT_concat.txt`) 旧 PASS (roK 4.6 桁)・新 PASS (4.3 桁)。
    F1 切り分け (`_f1split/F1_SPLIT.txt`): 旧バイナリの restart step 0 の rms_roK は sstSigmaBlend 1 / 0 とも床の 408.7 倍。~~→ F1 仮説は不支持~~ (撤回、codex result-2 M4: 旧版では初回の σ が両条件とも 0.85/0.5 になるので、この A/B は識別能力が無い) → **このA/B では識別不能**、#8d で識別する。
    新の dF1_inf の DRIFTING は継続課題 (#7) とし accepted の障害にしない (条件 (i)(iii) 成立、(ii) は #8d で完了)。

- **codex result-2 対応 (#8a–#8d、2026-09-26、新 f67fe877 sha256 `712141ec…`、AWS・block 128、commit c222ce62)** — 観測 → 帰属 → 限定つき解釈:
  - **#8a 面寄与ダンプ** (`_g0_lsq_seam/R8a_facedump.txt`、比較用ビルドのパッチ `notes/sessions/boundary-node-periodic-gradient-fix-8a-{olddump,wrong}.patch`):
    軸対称×周期 (5104 面、周期半割面 82) は旧+ダンプと新が species 10208・passive 5104 の全 (変数, 面) でビット同一。
    誤り版 (`periodicNodeActive`) は周期半割面の寄与が全部 0 (164/164、82/82) で新と差 → 検出力あり。
    並進周期 (TGV 32³ 三重周期 + 化学種・ξ、111078 面、周期半割面 6534) は周期半割面以外の全面で旧新ビット同一、周期半割面は新で 0・旧で非 0 (= 意図した修正そのもの)。
    R3 再確認 (新・ダンプ off、`R8a_R3_dumpoff.txt`) 30 配列 PASS。→ **PASS**。#6a の massflux 比較は検出力不足として置き換え (記録は残す)。
  - **#8c 定数場 3 区分** (`G_*.txt` 7 本): G0・G2・G1-a (0.86–2.32 ε)・G1-b (最大/閾値 ≤ 0.176、場ごと比 ≤ 1.81) PASS。壁なし継ぎ目の定数場 最大 3.03 ε|φ|/h。
    channel の壁∩継ぎ目 (ξ、96 点) は格納面ベクトルの閉包参照 $\phi\Sigma S_f/V$ (5.31 ε) との差 0.91 ε。k/ω の壁∩継ぎ目は壁ピンで対象 0。→ **PASS**。
  - **#8b G2′ 固定形状** (`G2p_fixedshape.txt`、`G2p_fixedshape_Halt.txt`): G2 参照 `lsq_merged_ref` をそのまま使い $h,h/2,h/4$。面 3・辺 3 group と内部は全成分 1.0000。
    **測定後の定義補正 (2 点、緩和ではない)**: (i) 「角 ≥ 2 group」は三重周期立方体では 8 隅が 1 group になるので構造的に満たせない指定ミス → box16/32/64 の角 group (ジッタ実現が別で形状も別) の 3 個で代替。
    (ii) 丸め床 (誤差 ≤ $10^3\varepsilon_{f64}S$) 以下の成分は次数を定義できないので判定対象外 — ただしこれに頼らず、全成分非零の二次形式 $H'$ ($H'_{zz}=0.3$、$H'_{xy}=0.2$) で追試:
    角 group の z 誤差は $H$ では 2.8e-17 (H の zz = 0 で消えていた) → $H'$ で 4.9e-4 になり、3 格子とも全成分の次数 1.0000。→ **PASS (定義補正つき)**。
    注記: 形状固定では誤差は $h$ に厳密比例し次数 ≈ 1 は構成上の帰結。示したのは合併作用素の 1 次整合。情報量は誤差定数 $C=e/h$ (中央値 面 0.030・辺 0.0087・角 0.0049 に対し内部 0.031) で、継ぎ目は内部より悪くない。
    ジッタ系列 (64→128 で 0.855) は診断・判定不能のまま保存 (#7)。スクリプトの VERDICT 行は当初の個数規則のままの「判定保留」を出す。
  - **#8d F1 切り分け** (`case/39.periodic_hills/_f1split/F1_SPLIT.txt`、パッチ `…-8d-m3only.patch`): 1266aba1 + m3 のみ版で同じ restart 入力の step 0 rms_roK/床 = **0.987** (旧版 408.5)。
    固定規則 (< 10) により**旧 restart スパイクの原因は F1 の初回上書き**と確定 (旧延長の `--from-floor` NOT CONVERGED の説明)。「m3 が他の restart 問題も直す」とは書かない。

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/calcGradient_d.cu`、`periodicNode_d.cu`、`ransTransport_d.cu`、`main.cpp`、`mesh/mesh.cpp`。
- **既存の node 周期 run (case/39、case/09、周期翼列など) の結果は継ぎ目付近で変わる** (修正)。再現には修正前の commit のバイナリが要る。
- docs: `methods/gradient.md`。
- 壁半割面込み GG の定数場閉包は float32 で約 5 ε|φ|/h (channel 実測、継ぎ目に依らない) で未対応。
- ジッタ格子での LSQ 勾配の次数は 0.83–1.0 (内部も同じ) で、0.9 を跨ぐ。後続の LSQ 統一 plan で扱う。
- 軸対称×周期・回転周期のスカラー勾配は片側 GG + 半割面込みのまま (既存の未修正挙動、本修正で不変)。
- GG を残す変数 (化学種・受動種・凝縮モーメント) で $N\varepsilon\phi/h$ 床が問題になるなら、$(\phi_f-\phi_i)S_f$ 形を検討する (本 plan では変えない: 全域の GG ビットが変わり R3 の帰属が濁る。k/ω は後続の LSQ 統一 plan で床ごと消え、その plan の §1 にこの床を根拠として引用する)。

## 8. 完了条件

**accepted の条件 (2026-09-26 `diagnostician` 固定、codex result 3 回目はこの表で採点する)**:

1. 中核: G0 全変種 PASS、G1-a ≤ 4ε、G1-b (上限 $2N_{max}\varepsilon\phi/h$ + 場ごと比 ≤ 2)、定数場 3 区分 (#8c)、G2 ≤ 1e-5·S 全水準。
2. G2′: 固定形状縮小試験で全成分 ≥ 0.9 (≥ 6 group、#8b)。ジッタ系列は診断として保存 (判定不能ラベルのまま)。
3. 同一性: R3 (非周期ビット同一)、#8a 面寄与ダンプで 旧+ダンプ = 新+ダンプ (軸対称×周期・並進周期)、誤り版で差が出る。
4. R1: 旧新とも連結区間 PASS、dF1_inf 以外 STEADY、dF1_inf は上限 0.05 内で継続課題、#8d の結果 (どちらでも) を記録。
5. R2: 保存上限 4 本 PASS、符号診断は「測定後の診断」と明記、結論は観測と commit 範囲への帰属に限定。
6. 文書同期 + 既知制約の登録 (壁閉包 5ε、軸対称×周期のスカラー勾配は未合併、ジッタ格子の LSQ 次数は後続 LSQ 化 plan)。


- [ ] `methods/gradient.md` の「修正中」を外す
- [ ] 実装・検証完了 (§6)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] `status` を `done` に変更し、§9 に変更ログを記載
- [ ] `plans/active/` → `plans/accepted/` へ移動、[`plans/README.md`](../README.md) を同期

## 9. 変更ログ

- `2026-09-26` — #8a–#8e 完了。§8 の accepted 条件 1–6 がすべて成立 (`diagnostician` 採点)。旧 restart スパイクは F1 初回上書きと確定。codex result 3 回目へ。

- `2026-09-26` — codex result 2 回目 NO-GO (C0/M4/m1) を全件採用 (#8a–#8e)。`diagnostician` の前回判断 2 件 (壁節点除外の根拠、F1 仮説の不支持) を撤回。accepted の条件を §8 に固定。

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

## 参考: `case/09.Taylor-Green/_g0_lsq_seam/R8a_facedump.txt`

```
#8a: スカラー GG の面寄与ダンプ (FORGE_DUMP_SCALARGRAD) を 3 バイナリで比較 (plan boundary-node-periodic-gradient-fix §5.1 #8a、codex result-2 M1)
run root: /home/ubuntu/pgrad6_scratch/r8a_<axi|trans>_<new|olddump|wrong>
比較は float32 のビット比較 (uint32 view)。面 = PLANES の index、周期半割面 = BCONDS/<physID>/iPlanes (kind periodic) の和集合
注: RUN_PROVENANCE の git_head は run_case.sh を置いたリポジトリ (forge-pgrad-new = f67fe877) の値で、バイナリの版ではない。バイナリは forge_sha256 で識別する (new = f67fe877、olddump = 1266aba1 + olddump.patch、wrong = f67fe877 + wrong.patch)

## axi
new: forge_sha256: 712141ec6289165b2b4cba1c11111d410f042ff9a45976da2ae581714f24db0c git_head    : f67fe877ab2fa5fcaa20dd2193fa7a21065b4846
  new log: [buildPeriodicNodeGroups] node periodic DOF identification: 41 slave CVs merged into masters (union-find).
  new log: [FORGE_DUMP_SCALARGRAD] wrote 2 x 5104 x 3 face contributions to /home/ubuntu/pgrad6_scratch/r8a_axi_new/sgdump.species
  new log: [FORGE_DUMP_SCALARGRAD] wrote 1 x 5104 x 3 face contributions to /home/ubuntu/pgrad6_scratch/r8a_axi_new/sgdump.passive
olddump: forge_sha256: 0ef63e4fb18e074601405713f5a49e06dc17fdfa2a11b1f4fc10c180e5585c69 git_head    : f67fe877ab2fa5fcaa20dd2193fa7a21065b4846
  olddump log: [buildPeriodicNodeGroups] node periodic DOF identification: 41 slave CVs merged into masters (union-find).
  olddump log: [FORGE_DUMP_SCALARGRAD] wrote 2 x 5104 x 3 face contributions to /home/ubuntu/pgrad6_scratch/r8a_axi_olddump/sgdump.species
  olddump log: [FORGE_DUMP_SCALARGRAD] wrote 1 x 5104 x 3 face contributions to /home/ubuntu/pgrad6_scratch/r8a_axi_olddump/sgdump.passive
wrong: forge_sha256: d7d3c56fbd91f6158d9d816d079618451ff612edc7a20bc5d137a9ab85722c90 git_head    : f67fe877ab2fa5fcaa20dd2193fa7a21065b4846
  wrong log: [buildPeriodicNodeGroups] node periodic DOF identification: 41 slave CVs merged into masters (union-find).
  wrong log: [FORGE_DUMP_SCALARGRAD] wrote 2 x 5104 x 3 face contributions to /home/ubuntu/pgrad6_scratch/r8a_axi_wrong/sgdump.species
  wrong log: [FORGE_DUMP_SCALARGRAD] wrote 1 x 5104 x 3 face contributions to /home/ubuntu/pgrad6_scratch/r8a_axi_wrong/sgdump.passive
  new res_1.h5 の VALUE 非有限: 0
  olddump res_1.h5 の VALUE 非有限: 0
  wrong res_1.h5 の VALUE 非有限: 0
nPlanes 5104 (内部 4900)、周期半割面 82、それ以外 5022
- species: nVar 2。非有限 new 0, olddump 0, wrong 0
  | 比較 | 周期半割面で差のある (変数,面) | それ以外で差のある (変数,面) |
  | --- | --- | --- |
  | olddump vs new | 0 / 164 | 0 / 10044 |
  | wrong vs new | 164 / 164 | 0 / 10044 |
  | olddump vs wrong | 164 / 164 | 0 / 10044 |
  | バイナリ | 周期半割面で 3 成分とも 0 の (変数,面) | それ以外で 3 成分とも 0 の (変数,面) |
  | --- | --- | --- |
  | new | 0 / 164 | 0 / 10044 |
  | olddump | 0 / 164 | 0 / 10044 |
  | wrong | 164 / 164 | 0 / 10044 |
  判定 (axi: olddump = new 全面ビット同一): ok
  検出力 (wrong は周期半割面が全て 0・その全てで new と差、他の面は new と同一): ok
- passive: nVar 1。非有限 new 0, olddump 0, wrong 0
  | 比較 | 周期半割面で差のある (変数,面) | それ以外で差のある (変数,面) |
  | --- | --- | --- |
  | olddump vs new | 0 / 82 | 0 / 5022 |
  | wrong vs new | 82 / 82 | 0 / 5022 |
  | olddump vs wrong | 82 / 82 | 0 / 5022 |
  | バイナリ | 周期半割面で 3 成分とも 0 の (変数,面) | それ以外で 3 成分とも 0 の (変数,面) |
  | --- | --- | --- |
  | new | 0 / 82 | 0 / 5022 |
  | olddump | 0 / 82 | 0 / 5022 |
  | wrong | 82 / 82 | 0 / 5022 |
  判定 (axi: olddump = new 全面ビット同一): ok
  検出力 (wrong は周期半割面が全て 0・その全てで new と差、他の面は new と同一): ok

## trans
new: forge_sha256: 712141ec6289165b2b4cba1c11111d410f042ff9a45976da2ae581714f24db0c git_head    : f67fe877ab2fa5fcaa20dd2193fa7a21065b4846
  new log: [buildPeriodicNodeGroups] node periodic DOF identification: 3169 slave CVs merged into masters (union-find).
  new log: gradLSQ=2 periodic seam: 2977 groups re-fitted on the merged stencil (0 truncated) (plan boundary-node-periodic-gradient-fix)
  new log: [FORGE_DUMP_SCALARGRAD] wrote 2 x 111078 x 3 face contributions to /home/ubuntu/pgrad6_scratch/r8a_trans_new/sgdump.species
  new log: [FORGE_DUMP_SCALARGRAD] wrote 1 x 111078 x 3 face contributions to /home/ubuntu/pgrad6_scratch/r8a_trans_new/sgdump.passive
olddump: forge_sha256: 0ef63e4fb18e074601405713f5a49e06dc17fdfa2a11b1f4fc10c180e5585c69 git_head    : f67fe877ab2fa5fcaa20dd2193fa7a21065b4846
  olddump log: [buildPeriodicNodeGroups] node periodic DOF identification: 3169 slave CVs merged into masters (union-find).
  olddump log: [FORGE_DUMP_SCALARGRAD] wrote 2 x 111078 x 3 face contributions to /home/ubuntu/pgrad6_scratch/r8a_trans_olddump/sgdump.species
  olddump log: [FORGE_DUMP_SCALARGRAD] wrote 1 x 111078 x 3 face contributions to /home/ubuntu/pgrad6_scratch/r8a_trans_olddump/sgdump.passive
wrong: forge_sha256: d7d3c56fbd91f6158d9d816d079618451ff612edc7a20bc5d137a9ab85722c90 git_head    : f67fe877ab2fa5fcaa20dd2193fa7a21065b4846
  wrong log: [buildPeriodicNodeGroups] node periodic DOF identification: 3169 slave CVs merged into masters (union-find).
  wrong log: gradLSQ=2 periodic seam: 2977 groups re-fitted on the merged stencil (0 truncated) (plan boundary-node-periodic-gradient-fix)
  wrong log: [FORGE_DUMP_SCALARGRAD] wrote 2 x 111078 x 3 face contributions to /home/ubuntu/pgrad6_scratch/r8a_trans_wrong/sgdump.species
  wrong log: [FORGE_DUMP_SCALARGRAD] wrote 1 x 111078 x 3 face contributions to /home/ubuntu/pgrad6_scratch/r8a_trans_wrong/sgdump.passive
  new res_1.h5 の VALUE 非有限: 0
  olddump res_1.h5 の VALUE 非有限: 0
  wrong res_1.h5 の VALUE 非有限: 0
nPlanes 111078 (内部 104544)、周期半割面 6534、それ以外 104544
- species: nVar 2。非有限 new 0, olddump 0, wrong 0
  | 比較 | 周期半割面で差のある (変数,面) | それ以外で差のある (変数,面) |
  | --- | --- | --- |
  | olddump vs new | 13068 / 13068 | 0 / 209088 |
  | wrong vs new | 0 / 13068 | 0 / 209088 |
  | olddump vs wrong | 13068 / 13068 | 0 / 209088 |
  | バイナリ | 周期半割面で 3 成分とも 0 の (変数,面) | それ以外で 3 成分とも 0 の (変数,面) |
  | --- | --- | --- |
  | new | 13068 / 13068 | 0 / 209088 |
  | olddump | 0 / 13068 | 0 / 209088 |
  | wrong | 13068 / 13068 | 0 / 209088 |
  判定 (trans: 周期半割面以外で olddump = new ビット同一 ok、周期半割面 new = 0 ok、olddump ≠ 0 ok): ok
- passive: nVar 1。非有限 new 0, olddump 0, wrong 0
  | 比較 | 周期半割面で差のある (変数,面) | それ以外で差のある (変数,面) |
  | --- | --- | --- |
  | olddump vs new | 6534 / 6534 | 0 / 104544 |
  | wrong vs new | 0 / 6534 | 0 / 104544 |
  | olddump vs wrong | 6534 / 6534 | 0 / 104544 |
  | バイナリ | 周期半割面で 3 成分とも 0 の (変数,面) | それ以外で 3 成分とも 0 の (変数,面) |
  | --- | --- | --- |
  | new | 6534 / 6534 | 0 / 104544 |
  | olddump | 0 / 6534 | 0 / 104544 |
  | wrong | 6534 / 6534 | 0 / 104544 |
  判定 (trans: 周期半割面以外で olddump = new ビット同一 ok、周期半割面 new = 0 ok、olddump ≠ 0 ok): ok

VERDICT #8a axi (軸対称 × 周期、olddump = new 全面ビット同一): PASS
VERDICT #8a trans (並進周期、半割面以外ビット同一・半割面 new 0 / olddump 非 0): PASS
VERDICT #8a 検出力 (wrong は axi で周期半割面の寄与 0 → new と差): PASS
VERDICT #8a (3 条件): PASS
```

## 参考: `case/09.Taylor-Green/_g0_lsq_seam/R8a_R3_dumpoff.txt`

```
# #8a の R3 再確認 (plan §5.1 #8a「R3 (非周期・既定 off) で全場ビット同一を再確認」、2026-09-26)
# この表の「旧」= bd22376d (HEAD 前。f67fe877 との solver 差分は FORGE_DUMP_SCALARGRAD のダンプのみ)、AWS で再ビルド sha256 90f3dc94…
#   (前回 #6a で使った bd22376d のビルド f0ad00c8… は forge-pgrad-new を f67fe877 に更新した際に消えたので作り直した。sha256 は一致しない)
# 「新」= f67fe877 (ダンプ off = FORGE_DUMP_SCALARGRAD 未設定)、sha256 712141ec…。入力は前回 R3 の r3_prep (r3_prepare.py 製) をそのまま転送。
# 判定規則は §6 R3 (r3_compare.py) のまま。

R3: 非周期・軸対称の勾配・リミタ配列 旧新比較 (plan boundary-node-periodic-gradient-fix §6 R3、2026-09-26)
run root: /home/ubuntu/pgrad8/r3_<case>_<old|new>_<a|b>
forge sha256: old 90f3dc94ddab…  new 712141ec6289…
列: 不一致数 (ビット比較) / 最大絶対差。判定は oo=old_a/old_b と on=old_a/new_a に §6 R3 の規則を当てる。

======================================================================================================================
case case48
NaN/Inf: なし (全 run・res_0/res_1 の全 VALUE)

-- res_1.h5 勾配・リミタ (30 配列)  [VERDICT 対象]
array             old_a/old_b        new_a/new_b        old_a/new_a  判定
dKdx           440  1.214e-02     459  1.173e-02     431  2.067e-02  PASS (ノイズ床)   (max|値| 1.167e+05)
dKdy            58  4.779e-01     103  9.366e-01      92  8.234e-04  PASS (ノイズ床)   (max|値| 6.250e+06)
dKdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegadx       562  1.251e+04     540  1.181e+04     486  1.432e+04  PASS (ノイズ床)   (max|値| 1.356e+11)
dOmegady       612  2.097e+06     574  1.573e+06     538  2.097e+06  PASS (ノイズ床)   (max|値| 7.264e+12)
dOmegadz         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dPdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dPdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dPdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dTdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dTdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dTdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUxdx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.195e+04)
dUxdy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 8.333e+05)
dUxdz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUydx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUydy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUydz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
drodx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
drody            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
drodz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_P        0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.000e+00)
limiter_T        0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_Ux       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.000e+00)
limiter_Uy       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.000e+00)
limiter_Uz       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.000e+00)
limiter_ro       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.000e+00)

-- res_0.h5 勾配・リミタ (30 配列)  [参考: step 前]
array             old_a/old_b        new_a/new_b        old_a/new_a  判定
dKdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dKdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dKdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegadx         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegady         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegadz         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dPdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dPdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dPdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dTdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dTdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dTdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUxdx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.195e+04)
dUxdy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 8.333e+05)
dUxdz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUydx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUydy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUydz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
drodx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
drody            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
drodz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_P        0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_T        0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_Ux       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_Uy       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_Uz       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_ro       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)

-- res_1.h5 保存量 (参考、判定なし)
array             old_a/old_b        new_a/new_b        old_a/new_a
ro               8  1.943e-05      20  4.804e-05      19  4.804e-05
roUx            62  2.384e-05      77  4.482e-05      60  4.482e-05
roUy          9429  3.356e-07   10262  8.299e-07    8899  8.299e-07
roUz             0  0.000e+00       0  0.000e+00       0  0.000e+00
roe              1  1.562e-02       2  1.562e-02       1  1.562e-02
roK              0  0.000e+00       0  0.000e+00       0  0.000e+00
roOmega          1  7.812e-03       1  7.812e-03       1  7.812e-03

VERDICT case48: PASS

======================================================================================================================
case case16
NaN/Inf: なし (全 run・res_0/res_1 の全 VALUE)

-- res_1.h5 勾配・リミタ (30 配列)  [VERDICT 対象]
array             old_a/old_b        new_a/new_b        old_a/new_a  判定
dKdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dKdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dKdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegadx         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegady         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegadz         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dPdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 2.844e+06)
dPdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 3.267e+06)
dPdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dTdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 3.525e+04)
dTdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 5.416e+05)
dTdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUxdx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 8.360e+05)
dUxdy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 8.025e+06)
dUxdz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUydx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.751e+05)
dUydy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 8.366e+05)
dUydz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
drodx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 9.396e+01)
drody            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 8.359e+02)
drodz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_P        0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.000e+00)
limiter_T        0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_Ux       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.000e+00)
limiter_Uy       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.000e+00)
limiter_Uz       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.000e+00)
limiter_ro       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.000e+00)

-- res_0.h5 勾配・リミタ (30 配列)  [参考: step 前]
array             old_a/old_b        new_a/new_b        old_a/new_a  判定
dKdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dKdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dKdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegadx         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegady         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegadz         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dPdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 2.844e+06)
dPdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 3.267e+06)
dPdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dTdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 3.525e+04)
dTdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 5.416e+05)
dTdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUxdx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 8.360e+05)
dUxdy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 8.025e+06)
dUxdz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUydx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.751e+05)
dUydy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 8.366e+05)
dUydz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
drodx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 9.396e+01)
drody            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 8.359e+02)
drodz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_P        0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_T        0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_Ux       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_Uy       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_Uz       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_ro       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)

-- res_1.h5 保存量 (参考、判定なし)
array             old_a/old_b        new_a/new_b        old_a/new_a
ro            1832  5.960e-08    1752  5.960e-08    1629  5.960e-08
roUx          6342  2.289e-05    5314  2.289e-05    5180  2.289e-05
roUy         25758  1.144e-05   24186  9.537e-06   23948  1.276e-05
roUz             0  0.000e+00       0  0.000e+00       0  0.000e+00
roe           1360  1.172e-02    1310  1.172e-02    1219  1.172e-02
roK              0  0.000e+00       0  0.000e+00       0  0.000e+00
roOmega          0  0.000e+00       0  0.000e+00       0  0.000e+00
roY0          1821  5.960e-08    1749  5.960e-08    1624  5.960e-08
roY1          2803  1.397e-09    2857  1.397e-09    2731  1.397e-09

VERDICT case16: PASS

======================================================================================================================
case axi
NaN/Inf: なし (全 run・res_0/res_1 の全 VALUE)

-- res_1.h5 勾配・リミタ (30 配列)  [VERDICT 対象]
array             old_a/old_b        new_a/new_b        old_a/new_a  判定
dKdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dKdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dKdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegadx         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegady         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegadz         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dPdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 8.408e+05)
dPdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 8.081e+05)
dPdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dTdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 3.584e+03)
dTdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 3.578e+02)
dTdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUxdx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 5.000e+02)
dUxdy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 5.925e+03)
dUxdz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUydx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 2.513e+02)
dUydy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 4.000e+02)
dUydz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
drodx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 3.626e+00)
drody            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 3.420e+00)
drodz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_P        0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.000e+00)
limiter_T        0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_Ux       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.000e+00)
limiter_Uy       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.000e+00)
limiter_Uz       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.000e+00)
limiter_ro       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.000e+00)

-- res_0.h5 勾配・リミタ (30 配列)  [参考: step 前]
array             old_a/old_b        new_a/new_b        old_a/new_a  判定
dKdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dKdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dKdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegadx         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegady         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dOmegadz         0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dPdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 8.408e+05)
dPdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 8.085e+05)
dPdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dTdx             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 1.137e+02)
dTdy             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 3.578e+02)
dTdz             0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUxdx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 5.000e+02)
dUxdy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 5.925e+03)
dUxdz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUydx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 2.513e+02)
dUydy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 4.000e+02)
dUydz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdy            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
dUzdz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
drodx            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 3.626e+00)
drody            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 3.420e+00)
drodz            0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_P        0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_T        0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_Ux       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_Uy       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_Uz       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)
limiter_ro       0  0.000e+00       0  0.000e+00       0  0.000e+00  PASS   (max|値| 0.000e+00)

-- res_1.h5 保存量 (参考、判定なし)
array             old_a/old_b        new_a/new_b        old_a/new_a
ro              12  1.192e-07      19  1.192e-07      17  1.192e-07
roUx           181  1.526e-05     191  1.526e-05     180  1.526e-05
roUy          1567  1.001e-05    1632  1.526e-05    1716  1.574e-05
roUz             0  0.000e+00       0  0.000e+00       0  0.000e+00
roe             24  3.125e-02      30  3.125e-02      29  3.125e-02
roK              0  0.000e+00       0  0.000e+00       0  0.000e+00
roOmega          0  0.000e+00       0  0.000e+00       0  0.000e+00
roY0            11  1.192e-07      17  1.192e-07      17  1.192e-07
roY1            25  3.725e-09      28  3.725e-09      23  3.725e-09

VERDICT axi: PASS

======================================================================================================================
SUMMARY: case48=PASS  case16=PASS  axi=PASS
VERDICT: PASS
```

## 参考: `case/09.Taylor-Green/_g0_lsq_seam/G_channel.txt`

```
# channel: plan boundary-node-periodic-gradient-fix §6 G0/G1/G2
run: /home/ubuntu/pgrad6_scratch/g_harness_channel
forge_mtime : 2026-09-26T06:22:46+00:00
forge_sha256: 712141ec6289165b2b4cba1c11111d410f042ff9a45976da2ae581714f24db0c
git_head    : f67fe877ab2fa5fcaa20dd2193fa7a21065b4846
harness revision: f67fe87 (dirty)
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
| k | 全成分 (判定) | 2.33 | 2.19 | 0.94 | ok |
| k | x (参考) | 2.33 | 2.05 | 0.88 | - |
| k | y (参考) | 1.61 | 1.37 | 0.85 | - |
| k | z (参考) | 1.96 | 2.19 | 1.12 | - |
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
| dK (GG) | 2.3 | 0.65 | 2.2 | 0.59 |
| dΩ (GG) | 2.3 | 0.52 | 1.8 | 0.54 |
| dξ (GG、Y の代理) | 2.4 | 0.55 | 2.3 | 0.53 |

## G1 定数場 (GPU、k・ω・ξ = const を焼いて 1 step、plan §6 G1・§5.1 #8c の 3 区分)
run: /home/ubuntu/pgrad6_scratch/g_harness_channel_const
forge_mtime : 2026-09-26T06:22:46+00:00
forge_sha256: 712141ec6289165b2b4cba1c11111d410f042ff9a45976da2ae581714f24db0c
git_head    : f67fe877ab2fa5fcaa20dd2193fa7a21065b4846
判定: 壁なし継ぎ目 max_c |∂φ| ≤ 4ε|φ|/h、壁∩継ぎ目 max_c |∂φ − φ·ΣS_f/V| ≤ 4ε|φ|/h (ΣS_f = 格納 float32 面ベクトルの double 和、周期半割面を除く・非周期境界面を含む、V = 合併体積)。非継ぎ目壁は同じ参照で報告、内部は参考。単位 ε|φ|/h (h = 節点の合併体積^(1/3))
除外 (現行どおり): 境界条件の上書きで res_0 が焼いた定数と違う節点と、それに stencil が触れる group (k/ω の壁ピン: 定数場が BC で壊れるので試験にならない)。ξ は上書きされないので除外 0
区分の節点数 (除外前): 壁なし継ぎ目 528、壁∩継ぎ目 96、非継ぎ目壁 210、内部 1155
| 場 | φ | 区分 | 節点 (除外) | 比較 | 最大 [ε|φ|/h] | 参照 |φ·ΣS_f/V| の最大 [同] | 判定 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| k | 2 | 壁なし継ぎ目 | 432 (96) | |∂φ| | 1.41 | 0.55 | ok |
| k | 2 | 壁∩継ぎ目 | 0 (96) | - | - | - | 対象 0 (除外) |
| k | 2 | 非継ぎ目壁 | 0 (210) | - | - | - | - |
| k | 2 | 内部 | 945 (210) | |∂φ| | 0.91 | 0.55 | 報告 |
| omega | 100 | 壁なし継ぎ目 | 26 (502) | |∂φ| | 1.47 | 0.52 | ok |
| omega | 100 | 壁∩継ぎ目 | 0 (96) | - | - | - | 対象 0 (除外) |
| omega | 100 | 非継ぎ目壁 | 0 (210) | - | - | - | - |
| omega | 100 | 内部 | 62 (1093) | |∂φ| | 0.81 | 0.55 | 報告 |
| Xi | 0.5 | 壁なし継ぎ目 | 528 (0) | |∂φ| | 0.91 | 0.57 | ok |
| Xi | 0.5 | 壁∩継ぎ目 | 96 (0) | |∂φ − φΣS/V| | 0.91 | 5.31 | ok |
| Xi | 0.5 | 非継ぎ目壁 | 210 (0) | |∂φ − φΣS/V| | 0.91 | 5.31 | 報告 |
| Xi | 0.5 | 内部 | 1155 (0) | |∂φ| | 0.99 | 0.57 | 報告 |
VERDICT G1 定数場 (channel, 3 区分: 壁なし継ぎ目 ≤ 4ε、壁∩継ぎ目 |∂φ − φΣS_f/V| ≤ 4ε): PASS
```

## 参考: `case/09.Taylor-Green/_g0_lsq_seam/G2p_fixedshape.txt`

```
# G2′ 固定形状縮小試験 (plan boundary-node-periodic-gradient-fix §5.1 #8b、codex result-2 M2)
harness revision: 58ea7151
メッシュ: /tmp/claude-1000/-home-sano-work-forge/4b0c8643-66fb-4fde-8878-7c8a5104c060/scratchpad/g_harness_mesh/box32_j0.2/mesh.h5 (32³、ジッタ ±0.2h、h = 2π/32 = 0.19635)
参照: gharness.lsq_merged_ref (G2 の合併 stencil LSQ の double 参照、GPU と ≤ 1e-5·S で一致済み) をそのまま呼ぶ
場: φ = a·x + ½ xᵀH x、原点 = group の root、a = (0.3, -0.5, 0.2)、H = 0.2·[[1,0,−0.35],[0,0.5,0.15],[−0.35,0.15,0]]。解析勾配 (原点) = a。誤差 = 成分ごとの |∇_ref − a|
縮小: 変位 d → s·d (s = 1, 1/2, 1/4、無次元形状固定)。h_s = s·h

注記: 形状を固定した縮小では LSQ の二次場の誤差は s に厳密に比例する (M = Σ w d dᵀ は s に不変、右辺は s に比例) ので、次数はほぼ 1 になる。示すのは合併作用素の 1 次整合であり、ジッタ系列 (G2p_jitter.txt、64→128 の継ぎ目次数 0.855) の「判定不能」は診断としてそのまま残す (この試験で覆さない)。

実在の継ぎ目 group (member 数: 個数): {2: 2883, 4: 93, 8: 1}

## group ごと (誤差 e_x/e_y/e_z、次数は h→h/2 / h/2→h/4、C = e/h (h = h_s、最細))
| 種別 | root | member | 自己検査 max|Δ| (s=1、全メッシュとの差) | e (s=1) x / y / z | e (s=1/4) x / y / z | 次数 h→h/2 x / y / z | 次数 h/2→h/4 x / y / z | 最小次数 | C = e/h (最細) x / y / z |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 面 (2) | 380 | 2 | 0.0e+00 | 2.106e-04 / 1.197e-03 / 3.372e-03 | 5.265e-05 / 2.991e-04 / 8.431e-04 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0011 / 0.0061 / 0.0172 |
| 面 (2) | 1838 | 2 | 0.0e+00 | 7.839e-03 / 1.768e-03 / 3.645e-03 | 1.960e-03 / 4.419e-04 / 9.112e-04 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0399 / 0.0090 / 0.0186 |
| 面 (2) | 3262 | 2 | 0.0e+00 | 5.132e-03 / 3.675e-03 / 5.918e-03 | 1.283e-03 / 9.187e-04 / 1.479e-03 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0261 / 0.0187 / 0.0301 |
| 辺 (4) | 8 | 4 | 0.0e+00 | 1.373e-03 / 2.973e-04 / 6.608e-04 | 3.433e-04 / 7.433e-05 / 1.652e-04 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0070 / 0.0015 / 0.0034 |
| 辺 (4) | 53 | 4 | 0.0e+00 | 3.234e-03 / 2.530e-03 / 1.141e-03 | 8.086e-04 / 6.326e-04 / 2.853e-04 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0165 / 0.0129 / 0.0058 |
| 辺 (4) | 286 | 4 | 0.0e+00 | 1.197e-03 / 3.329e-05 / 1.717e-03 | 2.993e-04 / 8.323e-06 / 4.292e-04 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0061 / 0.0002 / 0.0087 |
| 角 (8) | 0 | 8 | 0.0e+00 | 4.444e-04 / 9.576e-04 / 2.776e-17 | 1.111e-04 / 2.394e-04 / 2.776e-17 | 1.0000 / 1.0000 / 0.0000 | 1.0000 / 1.0000 / 0.0000 | 0.0000 | 0.0023 / 0.0049 / 0.0000 |
| 内部 (1) | 6146 | 1 | 0.0e+00 | 2.679e-03 / 2.243e-03 / 2.839e-03 | 6.699e-04 / 5.607e-04 / 7.098e-04 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0136 / 0.0114 / 0.0145 |
| 内部 (1) | 11729 | 1 | 0.0e+00 | 4.082e-04 / 7.600e-03 / 9.032e-03 | 1.021e-04 / 1.900e-03 / 2.258e-03 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0021 / 0.0387 / 0.0460 |
| 内部 (1) | 17313 | 1 | 0.0e+00 | 7.299e-03 / 6.119e-03 / 3.743e-03 | 1.825e-03 / 1.530e-03 / 9.356e-04 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0372 / 0.0312 / 0.0191 |
| 内部 (1) | 24768 | 1 | 0.0e+00 | 2.490e-04 / 1.905e-03 / 2.425e-05 | 6.225e-05 / 4.761e-04 / 6.063e-06 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0013 / 0.0097 / 0.0001 |
| 内部 (1) | 30352 | 1 | 0.0e+00 | 9.069e-03 / 2.373e-03 / 5.191e-03 | 2.267e-03 / 5.932e-04 / 1.298e-03 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0462 / 0.0121 / 0.0264 |
| 内部 (1) | 35936 | 1 | 0.0e+00 | 4.755e-03 / 2.811e-03 / 4.717e-03 | 1.189e-03 / 7.028e-04 / 1.179e-03 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0242 / 0.0143 / 0.0240 |

## 誤差定数 C = e/h (成分の最大) の比較
| 区分 | 個数 | C 最小 | C 中央値 | C 最大 |
| --- | --- | --- | --- | --- |
| 面 (2) | 3 | 0.0172 | 0.0301 | 0.0399 |
| 辺 (4) | 3 | 0.0070 | 0.0087 | 0.0165 |
| 角 (8) | 1 | 0.0049 | 0.0049 | 0.0049 |
| 内部 (1) | 6 | 0.0097 | 0.0307 | 0.0462 |

自己検査 (部分メッシュ抽出 = 全メッシュの lsq_merged_ref、s=1): 最大差 0.0e+00 (一致、判定 ≤ 1e-12)

## 参考 (判定外): 32³ 以外のジッタ格子の角 group (同じ生成器、別の格子)
- box16_j0.2 角 group root 0 (member 8): 最小次数 0.0000、e (s=1) = 2.462e-03 / 6.045e-04 / 2.776e-17
- box64_j0.2 角 group root 0 (member 8): 最小次数 0.0000、e (s=1) = 2.607e-04 / 7.168e-04 / 2.776e-17

継ぎ目 group の個数: 面 (2) 3、辺 (4) 3、角 (8) 1、計 7
継ぎ目の最小次数 (全成分・両対): 0.0000、内部の最小次数: 1.0000
誤差が全縮小率で倍精度丸め (≤ 1e-14) 以下の成分 (次数を定義できない。§5.1 #8b の基準に床の扱いは無いので判定には使わず事実のみ): 角 (8) root 0 の z 成分 (e = 2.8e-17 / 2.8e-17 / 2.8e-17)
それらを除いた継ぎ目の最小次数: 1.0000

VERDICT G2′ 固定形状 (面・辺・角 各 ≥ 2、計 ≥ 6、全成分の次数 ≥ 0.9): 判定保留 (基準の個数条件を満たせない: 角 (8) の実在 group が 32³ ジッタ格子に 1 個しか無い (三重周期立方体の 8 隅は 1 つの group)。取った group の最小次数 0.0000 は記録のみ)
```

## 参考: `case/09.Taylor-Green/_g0_lsq_seam/G2p_fixedshape_Halt.txt`

```
# G2′ 固定形状縮小試験 (plan boundary-node-periodic-gradient-fix §5.1 #8b、codex result-2 M2)
harness revision: c222ce62 (dirty)
メッシュ: /tmp/claude-1000/-home-sano-work-forge/4b0c8643-66fb-4fde-8878-7c8a5104c060/scratchpad/g_harness_mesh/box32_j0.2/mesh.h5 (32³、ジッタ ±0.2h、h = 2π/32 = 0.19635)
参照: gharness.lsq_merged_ref (G2 の合併 stencil LSQ の double 参照、GPU と ≤ 1e-5·S で一致済み) をそのまま呼ぶ
場: φ = a·x + ½ xᵀH x、原点 = group の root、a = (0.3, -0.5, 0.2)、H = [[ 0.2 , 0.04,-0.07],
 [ 0.04, 0.1 , 0.03],
 [-0.07, 0.03, 0.06]] (H'、--hess-alt)。解析勾配 (原点) = a。誤差 = 成分ごとの |∇_ref − a|
縮小: 変位 d → s·d (s = 1, 1/2, 1/4、無次元形状固定)。h_s = s·h

注記: 形状を固定した縮小では LSQ の二次場の誤差は s に厳密に比例する (M = Σ w d dᵀ は s に不変、右辺は s に比例) ので、次数はほぼ 1 になる。示すのは合併作用素の 1 次整合であり、ジッタ系列 (G2p_jitter.txt、64→128 の継ぎ目次数 0.855) の「判定不能」は診断としてそのまま残す (この試験で覆さない)。

実在の継ぎ目 group (member 数: 個数): {2: 2883, 4: 93, 8: 1}

## group ごと (誤差 e_x/e_y/e_z、次数は h→h/2 / h/2→h/4、C = e/h (h = h_s、最細))
| 種別 | root | member | 自己検査 max|Δ| (s=1、全メッシュとの差) | e (s=1) x / y / z | e (s=1/4) x / y / z | 次数 h→h/2 x / y / z | 次数 h/2→h/4 x / y / z | 最小次数 | C = e/h (最細) x / y / z |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 面 (2) | 380 | 2 | 0.0e+00 | 7.540e-04 / 1.889e-03 / 3.474e-03 | 1.885e-04 / 4.723e-04 / 8.685e-04 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0038 / 0.0096 / 0.0177 |
| 面 (2) | 1838 | 2 | 0.0e+00 | 8.401e-03 / 3.556e-03 / 3.759e-03 | 2.100e-03 / 8.891e-04 / 9.398e-04 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0428 / 0.0181 / 0.0191 |
| 面 (2) | 3262 | 2 | 0.0e+00 | 5.951e-03 / 4.198e-03 / 6.580e-03 | 1.488e-03 / 1.049e-03 / 1.645e-03 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0303 / 0.0214 / 0.0335 |
| 辺 (4) | 8 | 4 | 0.0e+00 | 1.197e-03 / 9.121e-04 / 7.811e-04 | 2.992e-04 / 2.280e-04 / 1.953e-04 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0061 / 0.0046 / 0.0040 |
| 辺 (4) | 53 | 4 | 0.0e+00 | 4.066e-03 / 3.540e-03 / 6.813e-04 | 1.017e-03 / 8.851e-04 / 1.703e-04 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0207 / 0.0180 / 0.0035 |
| 辺 (4) | 286 | 4 | 0.0e+00 | 1.216e-03 / 5.703e-05 / 1.350e-03 | 3.040e-04 / 1.426e-05 / 3.375e-04 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0062 / 0.0003 / 0.0069 |
| 角 (8) | 0 | 8 | 0.0e+00 | 4.444e-04 / 9.576e-04 / 4.938e-04 | 1.111e-04 / 2.394e-04 / 1.234e-04 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0023 / 0.0049 / 0.0025 |
| 内部 (1) | 6146 | 1 | 0.0e+00 | 3.273e-03 / 4.518e-03 / 4.610e-03 | 8.184e-04 / 1.129e-03 / 1.153e-03 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0167 / 0.0230 / 0.0235 |
| 内部 (1) | 11729 | 1 | 0.0e+00 | 5.394e-04 / 7.929e-03 / 9.933e-03 | 1.349e-04 / 1.982e-03 / 2.483e-03 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0027 / 0.0404 / 0.0506 |
| 内部 (1) | 17313 | 1 | 0.0e+00 | 6.613e-03 / 7.076e-03 / 2.701e-03 | 1.653e-03 / 1.769e-03 / 6.753e-04 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0337 / 0.0360 / 0.0138 |
| 内部 (1) | 24768 | 1 | 0.0e+00 | 3.826e-05 / 2.828e-03 / 9.474e-04 | 9.566e-06 / 7.069e-04 / 2.369e-04 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0002 / 0.0144 / 0.0048 |
| 内部 (1) | 30352 | 1 | 0.0e+00 | 1.088e-02 / 4.212e-03 / 6.492e-03 | 2.721e-03 / 1.053e-03 / 1.623e-03 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0554 / 0.0215 / 0.0331 |
| 内部 (1) | 35936 | 1 | 0.0e+00 | 6.593e-03 / 4.106e-03 / 4.969e-03 | 1.648e-03 / 1.026e-03 / 1.242e-03 | 1.0000 / 1.0000 / 1.0000 | 1.0000 / 1.0000 / 1.0000 | 1.0000 | 0.0336 / 0.0209 / 0.0253 |

## 誤差定数 C = e/h (成分の最大) の比較
| 区分 | 個数 | C 最小 | C 中央値 | C 最大 |
| --- | --- | --- | --- | --- |
| 面 (2) | 3 | 0.0177 | 0.0335 | 0.0428 |
| 辺 (4) | 3 | 0.0061 | 0.0069 | 0.0207 |
| 角 (8) | 1 | 0.0049 | 0.0049 | 0.0049 |
| 内部 (1) | 6 | 0.0144 | 0.0348 | 0.0554 |

自己検査 (部分メッシュ抽出 = 全メッシュの lsq_merged_ref、s=1): 最大差 0.0e+00 (一致、判定 ≤ 1e-12)

## 参考 (判定外): 32³ 以外のジッタ格子の角 group (同じ生成器、別の格子)
- box16_j0.2 角 group root 0 (member 8): 最小次数 1.0000、e (s=1) = 2.462e-03 / 6.045e-04 / 1.924e-03
- box64_j0.2 角 group root 0 (member 8): 最小次数 1.0000、e (s=1) = 2.607e-04 / 7.168e-04 / 1.031e-04

継ぎ目 group の個数: 面 (2) 3、辺 (4) 3、角 (8) 1、計 7
継ぎ目の最小次数 (全成分・両対): 1.0000、内部の最小次数: 1.0000

VERDICT G2′ 固定形状 (面・辺・角 各 ≥ 2、計 ≥ 6、全成分の次数 ≥ 0.9): 判定保留 (基準の個数条件を満たせない: 角 (8) の実在 group が 32³ ジッタ格子に 1 個しか無い (三重周期立方体の 8 隅は 1 つの group)。取った group の最小次数 1.0000 は記録のみ)
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

## #8d (plan §5.1 #8d、codex result-2 M4): m3 だけを当てた旧版で識別 (2026-09-26)
上の #6d の A/B (旧版の sstSigmaBlend 1 / 0) は、旧版では初回の σ が両条件とも 0.85 / 0.5 になるので識別能力が無い (codex result-2 M4)。
m3 のみ版: 1266aba1 + notes/sessions/boundary-node-periodic-gradient-fix-8d-m3only.patch (sstF1 の初期値 1 を allocVariables へ、
buildScalarDescs の初回 1 充填を削除。2ee5f0a7 の F1 部分と同じ。LSQ・gather・半割面除外などの他の変更は入れない)、
AWS g5 でクリーンビルド、sha256 ef6a91c046d8c34e3dc8c365782ec4ca58fc47432d45f60979d283dc101a35ad。block 128。
入力・設定は上と同じ (run_0038 の ic_S3_ext.h5 = run_0036 res_800000 の保存量 7 配列とビット一致を再確認、sstF1 は入力に無い。
run_0038 の solverConfig.yaml に nStepOuter 1・outStepInterval 1・turbulence.sstSigmaBlend 1 を明記しただけ)。
AWS run: ~/forge-pgrad-new/case/39.periodic_hills/_f1split/{m3_b1,old_b1_rerun} (スクラッチ)
参照床: run_0036 の rms_roK の末尾 20 % 平均 (check_convergence --from-floor の定義、outer_end 行 800000 のうち末尾 160000) = 3.291146e-09
(rms_roOmega は 1.253727e-03)

| run | バイナリ・設定 | step 0 rms_roK | ÷ 参照床 | step 0 rms_roOmega | ÷ 同床 | step 0 rms_ro |
| --- | --- | --- | --- | --- | --- | --- |
| m3_b1 | 1266aba1 + m3 のみ (sha256 ef6a91c0…)、sstSigmaBlend 1 | 3.247565328e-09 | 0.987 | 1.277059e-03 | 1.019 | 3.874474e-10 |
| old_b1_rerun | 旧 1266aba1 (sha256 f0505fe8…)、sstSigmaBlend 1 (同セッションの対照、判定外) | 1.344498742e-06 | 408.5 | 1.439659e-03 | 1.148 | 3.874605e-10 |
| (old_b1、上表) | 旧 1266aba1、sstSigmaBlend 1 | 1.344497832e-06 | 408.5 | 1.441476e-03 | 1.150 | 3.874670e-10 |

(上表 #6d の「÷ 参照床 408.7」は 3.29e-09 で割った値。ここでは 3.291146e-09 で割っている)
規則 (§5.1 #8d、測る前に固定): m3 のみ版で rms_roK/床 < 10 → 「F1 初回上書きが原因」、≥ 10 → 「原因未特定 (F1 でない)」
VERDICT #8d: F1 初回上書きが原因 (m3 のみ版の step 0 rms_roK/床 = 0.987 < 10)
```

## 参考: `notes/reviews/2026-09-26-boundary-node-periodic-gradient-fix-result-2.md`

```
# codex レビュー: boundary-node-periodic-gradient-fix (result)

- **plan**: [`plans/active/boundary-node-periodic-gradient-fix.md`](../../plans/active/boundary-node-periodic-gradient-fix.md)
- **stage**: `result` (diff base `1266aba1`)
- **date**: 2026-09-26
- **commit**: `80a659b7` (feature/sern-design)
- **codex**: effort `high`, 4.6 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M4/m1
- **focus**: 2 回目 (1 回目 NO-GO C0/M5/m1 への対応確認): (1) M1 の条件統一と #6a を決定論的試験 (面流束・状態のビット一致) で閉じたことの妥当性 (2) #6c の比 (場ごと) と定数場の壁節点除外が事後の緩和でないか (3) #6b 判定不能を accepted の障害にしない扱い (4) R1 の限定合格 + dF1_inf 定常性を継続課題とする扱い、旧延長の連結区間 PASS と F1 切り分け (5) R2 の結論が観測と帰属に限定されているか (6) accepted にしてよいか
- **extra**: `case/09.Taylor-Green/_g0_lsq_seam/README.md`, `case/09.Taylor-Green/_g0_lsq_seam/R6a_axi_periodic.txt`, `case/09.Taylor-Green/_g0_lsq_seam/R6a_massflux_bitcheck.txt`, `case/09.Taylor-Green/_g0_lsq_seam/G_channel.txt`, `case/09.Taylor-Green/_g0_lsq_seam/G2p_jitter.txt`, `case/09.Taylor-Green/r2_ke_budget.txt`, `case/39.periodic_hills/run_0038_r1_gradfix_old_ext/CONVERGENCE_VERDICT_concat.txt`, `case/39.periodic_hills/run_0039_r1_gradfix_new_ext/CONVERGENCE_VERDICT_concat.txt`, `case/39.periodic_hills/_f1split/F1_SPLIT.txt`, `notes/reviews/2026-09-26-boundary-node-periodic-gradient-fix-result.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

合併 LSQ、SST 勾配の合算位置、`periodicSeamMergeActive` への条件統一は支持します。  
ただし、追加試験2件に検出力の問題があり、G1・G2′の未達を完了扱いにする根拠も不足しています。現状での `accepted` 移行は認められません。

1. **Major — #6a の質量流束ビット一致では、修正対象のスカラー作用素を検証できません。**

   **根拠:** [R6a_massflux_bitcheck.txt:7](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/R6a_massflux_bitcheck.txt:7) が比較したのは `massflux` 5104面と `ro,Ux,Uy,Uz,P,sonic` です。`Y`・`Xi`・その勾配・スカラー面流束は含まれません。

   この試験の `speciesFaceReconstruction: 1` では、組成勾配は面エンタルピーに使われますが、通常の質量流束式には入りません。[convectiveFlux_slau_d.inc.cuh:579](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:579) の `mdot` が一致しても、組成再構成の正しさは保証されません。受動種についても、[speciesTransport_d.cu:1275](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1275) は再構成面値を使う条件を `>=2` としています。

   条件統一自体は正しく、コードから軸対称で旧新とも半割面を除外しないことは確認できます。しかし、このダンプから「スカラー面作用素も同一」「場の差はすべて集積順序」と結論するのは過剰です。

   **対案:** 同一入力で `species_gradient_d` が生成する**集積前の面寄与**を比較する決定論的試験に置き換えてください。化学種・受動種の両呼出しを対象にし、誤った `periodicNodeActive` 条件へ戻すと失敗することも確認します。既存の `R3` 規則上の `FAIL` は履歴として保持します。

2. **Major — G2′の次数未達は、内部との誤差比だけでは閉じられません。**

   **根拠:** [G2p_jitter.txt:20](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/G2p_jitter.txt:20) の本判定は次のとおりです。

   | 成分 | 継ぎ目の次数 | 内部の次数 |
   |---|---:|---:|
   | `φ_y` | **0.855** | 0.976 |
   | `φ_z` | **0.869** | 0.939 |

   どちらも継ぎ目で要求の0.9未満です。128³で誤差が内部より小さいことは支持できますが、**誤差の大小と収束次数は別の条件**です。内部の別成分に0.9未満があることから、この2成分の未達を「継ぎ目固有でない」と確定することもできません。

   [plan:218](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:218) は「判定不能」を記録しながら完了扱いにしていますが、§5.1 #7には次数確認の残作業がありません。

   **対案:** G2′を未達として戻してください。追加試験は、非対称な局所 stencil の無次元形状を固定し、`h,h/2,h/4` に縮小して継ぎ目全成分の次数を測る方式を推奨します。現在のジッタ系列は保存し、格子ごとに評価点集合が変わる最大誤差の診断として併記します。直ちに256³の本計算を追加する必要はありません。

3. **Major — 定数場から壁節点を除く変更は、検証範囲の事後縮小です。**

   **根拠:** [G_channel.txt:98](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/G_channel.txt:98) では、`k/ω` は境界条件で定数場が壊れた stencil を除外しています。一方、`Xi` は**除外0点**で定数場を保ち、最大誤差 `5.31 ε|φ|/h`、正式な判定は **FAIL** です。

   したがって、`Xi` の壁節点除外を「k/ωと同じ扱い」とする [plan:215](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:215) の説明は成立しません。壁∩継ぎ目と非継ぎ目壁で最大値が同じことは既存の壁閉包誤差を示唆しますが、当初対象だった壁∩継ぎ目の検証を不要にはしません。ハーネスも変更後の区分を判定していません。

   **対案:** 壁なし継ぎ目・壁∩継ぎ目・非継ぎ目壁を分けて機械判定してください。壁なしには現行の `4ε` を維持し、壁∩継ぎ目には、格納面ベクトルによる閉包誤差を含む参照との比較を追加します。絶対定数保存の未達は、根拠付きの既知制約として別途登録します。

   なお、**G1-bを場ごとの最大誤差ノルムで評価すること自体は妥当**です。全成分に絶対上限も課しており、成分別の不安定な比を参考値として残す整理は支持します。

4. **Major — #6d の `sstSigmaBlend=0/1` 比較には、F1初回上書き仮説を識別する能力がありません。**

   **根拠:** 旧版 `1266aba1` の `ransTransport_d.cu:105–118` は、`sstSigmaBlend=1` の初回に `F1=1` を書きます。一方、`=0` は固定係数を使います。[scalarTransport_d.cu:213](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/scalarTransport_d.cu:213) の式では、両者とも初回の係数は

   `σ_k=0.85、σ_ω=0.5`

   になります。したがって、[F1_SPLIT.txt:9](/home/sano/work/forge-sern-design/case/39.periodic_hills/_f1split/F1_SPLIT.txt:9) の両条件でスパイクが408.7倍になることは、上書き仮説と矛盾しません。「F1仮説は不支持」とする判断は撤回すべきです。新バイナリとの比較は勾配修正も含むため、この原因だけを分離できません。

   **対案:** 旧版を基準に、**F1初回上書きだけを除いた版**との同一restart入力・初回残差比較を行ってください。それまでは「このA/Bでは識別不能」と記録し、#6dを完了扱いにしないでください。

5. **Minor — 残作業表とrun索引に、撤回済み・更新前の記述が残っています。**

   **根拠:** [plan:124](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:124) の #5d は、現在も「単調収束中・漸近値」「新 `|r|≤0.05`」「継ぎ目由来の注入」を現行手順として保持しています。いずれも結果欄では撤回済みです。

   [case/39 README:78](/home/sano/work/forge-sern-design/case/39.periodic_hills/README.md:78) は800kへの回し直し中のままで、延長行にも連結区間の `PASS` が未反映です。#7には壁閉包・G2′の継続課題もありません。

   **対案:** 現行判定を§5.1・§6・run索引で統一し、旧判断は履歴と明示してください。後続へ移す課題には、追跡先planと完了条件を付けます。

指定の `git diff 1266aba1...HEAD -- solver_density_cuda design methods procedures` を取得し、HEAD `80a659b7` まで確認しました。中核実装に新たな符号・単位・適用条件の誤りは見つけていません。`methods/gradient.md` の現在仕様化と、化学種の代理として同じGGカーネルを通るトレーサを用いる説明も支持します。

R1のCSVを `check_quasisteady.py` で再判定した結果は、旧 `run_0038_r1_gradfix_old_ext/` が **ALL STEADY**、新 `run_0039_r1_gradfix_new_ext/` が **dF1_inf: DRIFTING / NOT ALL STEADY** でした。新の末尾64点は `0.03334079–0.03344835`。**観測窓内の上限成立に限った限定合格**と、定常性を継続課題にする扱いは支持します。連結区間の収束は保存VERDICTで旧新とも **PASS** ですが、残差元CSVがローカルにないため再判定はできていません。

R2は `run_0173–0176` の保存CSVから全4本の保存上限成立を再確認しました。符号診断も `max r_new=-0.0074228083`、早期の旧新差最小 `+0.0175604664`、旧の正値15点連続を再現しました。現在の「観測とcommit範囲への帰属」に限定した結論は妥当です。

**推奨は、`active` に留め、上記4件の検証を閉じて文書を同期した後に再レビューすることです。** ファイルは変更しておらず、指摘はplan未反映です。

指摘数: Critical 0 / Major 4 / Minor 1
```

## 出力形式

1. 冒頭に **判定 (GO / GO-with-changes / NO-GO)** と 3 行以内の要約。
2. 指摘一覧 (Critical → Major → Minor の順、番号付き。各項目に根拠と対案)。
3. 推奨 (1 つに絞る)。
4. 末尾に `指摘数: Critical N / Major N / Minor N` の 1 行。
