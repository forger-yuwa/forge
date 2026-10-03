# 非構造 2 次 FV のリミッタ — 主要コードの既定・比較研究の結論・2018–2026 の新しい研究

## メタ

- **area**: `limiter`
- **status**: `draft` (**調査専用 / コード変更なし / 推奨は決めていない** — 候補の順位付けまで)
- **前提となる調査**: [`limiter-unstructured-convergence-survey.md`](limiter-unstructured-convergence-survey.md)
  (以下「前回調査」)。ε のスケーリング・Wang の範囲 ε・SU2/FUN3D/Wind-US の ψ 凍結・境界の扱い・Berger らの chatter は
  **前回調査に書いたので繰り返さない**。本ノートは「いま何が既定・推奨か」「比較研究は何を『最良』としたか」
  「2018 年以降に何が出たか」「forge にどう当てはまるか」に絞る。
- **related_docs**: [`methods/limiter.md`](../../methods/limiter.md) (無次元化 Venkatakrishnan `limiterScaled: 1` の現行仕様)
- **related_plans**: [`limiter-inlet-column-oscillation.md`](../../plans/accepted/limiter-inlet-column-oscillation.md)、
  [`limiter-config-simplify.md`](../../plans/active/limiter-config-simplify.md) §4.3
- **created**: `2026-10-03`
- **一次資料の置き場**: SU2 ソースは `/home/sano/work/forge/.external/su2-src` (commit `12eb826f`)。
  ダウンロードした PDF (FUN3D 講習資料・CFL3D v5 マニュアル・Ahmad ら 2023・Kitamura–Hashimoto 2016・May–Berger) は
  作業用スクラッチにのみ置いた (リポジトリには入れていない)。

**確認度の表記** (前回調査と同じ):

- **[一次・ソース]** — ソースコードを読んで確認した。
- **[一次・文書]** — 公式マニュアル・講習資料・論文本文 (または公式アブストラクト) を読んで確認した。
- **[二次]** — 他論文の記述・検索結果の要約など。一次資料を直接読めていない。
- **[未確認]** — 記憶・推測、または探したが見つからなかった。採否の根拠に使わないこと。
- 本ノートで一次資料の式から自分で式変形した結果は「**[一次・文書から導出]**」と書く (出典の主張ではない)。

---

## 1. 問い

forge は node 中心 median dual、エッジ中点再構成の 2 次 MUSCL、SLAU 流束、block-DPLUR の定常擬似時間、float32 で動く。
既定リミッタは無次元化 Venkatakrishnan ($\hat\delta=\delta/q_\mathrm{ref}$、$\hat\epsilon^2=(K h_i/L_\mathrm{ref})^3$、$K=0.05$)。
旧 Barth–Jespersen も選べ、Nishikawa R1 は実装だけ残してコメントアウトしてある ([`methods/limiter.md`](../../methods/limiter.md))。

既知の問題として、定常残差の床が旧経路より約 20 倍高く、床は ψ が 0 付近と 1 付近を行き来する少数の節点に集中している。
ε を領域一定 (2e-6、許容オーバーシュートで上限を決めた値) にしても床は約 0.65 倍にしか下がらなかった。
報告量への影響は無い。

この状況で次を問う。

1. 主要コードは何を既定にし、何を推奨しているか (SU2・FUN3D・CFL3D/OVERFLOW・USM3D・Kestrel・TAU・elsA・OpenFOAM・Fluent・
   STAR-CCM+・CONVERGE・Code_Saturne・Cart3D・Loci/CHEM)。
2. 比較研究は、精度 (滑らかな極値を削らないか)・衝撃での頑健性・定常収束のそれぞれで、何を「最良」としているか。
3. 2018〜2026 年の新しい研究 (MLP 系・Nishikawa・微分可能リミッタ・平滑度センサ・a posteriori・機械学習・
   正値性/エントロピー・低マッハ対応) は何をしていて、どこまで実用化され、node 中心エッジベースのコードに入れるといくら掛かるか。
4. それぞれが forge (node 中心 median dual・エッジ中点再構成・陰的定常・float32・強い衝撃と境界層の両立) に合うか。

---

## 2. 主要コードの既定と推奨

### 2.1 一覧

「既定」はオプションを書かないときに使われる値、「推奨」は公式文書・講習資料の推奨文言。両者を分けて書く。

| コード | 離散化 | 既定 | 公式の推奨・実際の使われ方 | 選択肢 (抜粋) | 確認度 |
| --- | --- | --- | --- | --- | --- |
| **SU2** | node 中心 median dual、エッジ | `SLOPE_LIMITER_FLOW = VENKATAKRISHNAN`、`VENKAT_LIMITER_COEFF = 0.05`、`LIMITER_ITER = 999999` (乱流も `VENKATAKRISHNAN`、化学種・熱は `NONE`) | 文書: Venkatakrishnan は BJ に近く「収束性が大きく改善」、`VAN_ALBADA_EDGE` は「BJ よりやや拡散的だが収束は頑健」。Wang 版は「メッシュ寸法に直接依存しないので無次元化なしで使える」 | `NONE` / `VENKATAKRISHNAN` / `VENKATAKRISHNAN_WANG` / `BARTH_JESPERSEN` / `VAN_ALBADA_EDGE` / `NISHIKAWA_R3/R4/R5` (+ 随伴用 `SHARP_EDGES` / `WALL_DISTANCE`) | 一次・ソース (`Common/src/CConfig.cpp:2006–2094`、`config_template.cfg:1549`)、一次・文書 (SU2 docs) |
| **FUN3D** | node 中心 median dual、エッジ | `flux_limiter = "none"` (2014 マニュアルの名前リスト既定) | 2014 マニュアル: 完全気体の極超音速では `hvanleer` か `hvanalbada`。2023 講習 (v14.0): 「`hvanalbada` が推奨リミッタ」、鈍頭物体は `hlle++` + `hvanalbada`、「リミッタの buzzing はよくある、凍結できる」。2010 講習: 随伴を解くなら `hvanleer`/`hvanalbada` (stencil 型が必須) | `none`/`barth`/`venkat`/`minmod`/`vanleer`/`vanalbada`/`smooth`/`h*` 系、v14.0 から `nishikawa`、v14.2 (2025-04) から単独の `hpressure` | 一次・文書 (2014 マニュアル第 6 章、講習資料 2010・2023、リリースノート) |
| **CFL3D** (構造、比較用) | セル中心、κ 法 | `iflim = 0` (無制限) | v5 マニュアル: 「多くの純亜音速流はリミッタ不要」。κ≠1/3 なら minmod (`iflim=2`)、κ=1/3 なら「小さな勾配域で制限を切るカットオフ付きの滑らかリミッタ」(`iflim=3`)。v6: 「κ=1/3 でリミッタが要るときは常に `iflim=4`」— カットオフを総セル数から決める**等方・全域の値** | 0 / 1 (smooth) / 2 (minmod) / 3 / 4 | 一次・文書 (v5 マニュアル、v6 New Features) |
| **OVERFLOW** (構造、比較用) | 有限差分、重合格子 | 既定は未確認 (講習資料は中心差分 `IRHS=0` を最初に挙げる) | 上流化 (`IRHS=3–6`) で `ILIMIT` 1 (Koren) / 2 (minmod) / 3 (van Albada)、強い衝撃用の `DELTA` 修正 | 同左 + WENO | 一次・文書 (2010 講習資料)、既定値は未確認 |
| **USM3D** | セル中心 | 未確認 | 未確認 (前回調査で HANIS 論文に凍結の記述が無いことのみ確認) | 未確認 | 未確認 |
| **Kestrel (KCFD)** | セル中心、最小二乗再構成 | 未確認 | 未確認 | 未確認 | 未確認 (検索で見つけた資料は他コードの記述が混在) |
| **DLR TAU** | node 中心 dual (forge と同系統) | 未確認 | SBPW2 の DLR 発表では AUSMDV 2 次 + **Venkatakrishnan** + Green–Gauss | 2 次上流のリミッタ関数として Barth–Jespersen・Venkatakrishnan・SRR | 一次・文書 (SBPW2 発表資料)、既定は未確認 |
| **elsA** | 構造 (+ elsA-Hybrid) | 未確認 | 未確認 | Roe + MUSCL で van Albada・van Leer・minmod・Superbee | 二次 (CERFACS 技術報告の検索要約) |
| **OpenFOAM** | セル中心 | **既定なし** (`fvSchemes` に必ず書く) | `cellLimited<Venkatakrishnan>` のソースに「1 で切っているので微分不可能、`cubic` を推奨」の注記 | `cellLimited` / `cellMDLimited` / `faceLimited` / `cellLimited<cubic>` / `<Venkatakrishnan>` | 一次・ソース (前回調査 §2.2) |
| **Ansys Fluent** | セル中心 | **standard limiter** (minmod 関数、BJ 由来)、制限方向は **cell-to-face** | 理論ガイド: 微分不可能リミッタは「数桁で見かけの残差が停滞しがち」、それを避けるため differentiable limiter (Venkatakrishnan の修正形)。multidimensional limiter は面法線成分だけを切って散逸を減らす (Makarov を参照) | standard / multidimensional / differentiable、方向 cell-to-face / cell-to-cell / cell-to-node | 一次・文書 (Theory Guide v24.2) |
| **Simcenter STAR-CCM+** | セル中心 | 未確認 | HLPW3 では Siemens が「LSQ/Green–Gauss ハイブリッド + Venkatakrishnan」、Oxford が「Green–Gauss + min-mod」 | Venkatakrishnan / min-mod (他は未確認) | 二次 (HLPW3 要約論文の検索要約) |
| **CONVERGE** | セル中心 | 未確認 | `MUSCL_CVG` は 3D 勾配型リミッタとして BJ 型 minmod か Venkatakrishnan | 未確認 | 二次 (CFD-Wiki FAQ) |
| **Code_Saturne** | セル中心 (主に非圧縮) | 2 次 (中心/SOLU) に **slope test** (振動検出時に局所 1 次へ切替)。`isstpc` の既定は版で記述が異なる | — | slope test / NVD・TVD (`ischcv=4`) | 二次 (検索要約。版ごとの既定は未確認) |
| **NASA Cart3D** | セル中心 cut-cell | `Limiter 1` (= Barth–Jespersen) | 選択肢は「散逸の小さい順」に並ぶ | 0 なし / 1 BJ / 2 van Leer / 3 sin / 4 van Albada / 5 minmod | 一次・文書 (`input.cntl` 解説) |
| **Loci/CHEM** | セル中心 | 未確認 | 未確認 (前回調査: Luke の次元整合 ε は引用経由のみ) | 未確認 | 未確認 |

### 2.2 一覧から読めること

- **非構造の汎用・研究系コードの既定は Venkatakrishnan が多数派** (SU2・TAU の実使用・STAR-CCM+ の実使用・JAX-FVM などの新しいコード
  [一次・文書: arXiv 2607.07385 のアブストラクト])。**商用汎用コード (Fluent) の既定は微分不可能な minmod/BJ 型**で、
  収束が問題になる場合に differentiable を選ばせる構成。Cart3D も BJ が既定 [各行の出典]。
- **NASA FUN3D は既定が「リミッタなし」**で、衝撃がある問題では **h 系 (圧力ベースの補助リミッタ付き)** を推奨し、
  2023 年の講習でも `hvanalbada` を推奨にしている。Nishikawa 型 (`nishikawa`) は v14.0 から使えるが、講習資料で
  推奨に格上げされた記述は今回見つけていない [一次・文書]。
- **構造格子の老舗 (CFL3D) も「亜音速なら無制限、要るときだけ制限」**で、要るときのリミッタは「小さな勾配域で切るカットオフ付き」、
  そのカットオフは v6 で**格子の分割に依らない全域の値 (総セル数基準)** に改められている [一次・文書]。
  これは SU2 既定の「ε は全域定数」(前回調査 §2.2) と同じ向きの設計。
- **SU2 の `BARTH_JESPERSEN` は名前だけで、中身は「ε = 計算機イプシロンの Venkatakrishnan 関数」**
  [一次・ソース: `CLimiterDetails.hpp:145–168` の `limiterFunction` が `venkatFunction(proj, delta, eps2)`、`eps2 = epsilon()`]。
  2026-09 に「正しい $\min(1,\Delta_+/\Delta_-)$ を実装する」PR #2911 が出たが、保守者 (pcarruscag) が
  「**このリミッタは削除しよう、ひどい、何をやってもうまく収束しない**」とコメントし、PR は削除案に書き換えられたうえで
  **未マージのまま閉じられた** (2026-09-18) [一次・文書: GitHub PR #2911 と API の `merged: false`]。
  「本物の BJ は収束しない」は保守者の発言で、定量データは PR に無い。
- **SU2 の `VAN_ALBADA_EDGE` は近傍の最大/最小を使わないエッジ型**で、エッジ両端の差 $\delta_{ij}=q_j-q_i$ と各端の
  U-MUSCL 射影 $\mathrm{proj}$ だけから $\psi=(\delta_{ij}^2+10^{-6}+\mathrm{proj}\,\delta_{ij})/(\mathrm{proj}^2+\delta_{ij}^2+10^{-6})$
  を作る。**$10^{-6}$ は次元付きの固定値** (SU2 既定は SI のまま解くので、変数ごとに効き方が変わる) [一次・ソース:
  `numerics_simd/flow/convection/common.hpp:147–176`]。乱流・化学種へ広げる PR #2908 は未マージで閉じられた [一次・文書: GitHub API]。
  SU2 には 2025-10 に MUSCL の κ 法 (`MUSCL_KAPPA_FLOW`) が入っている (PR #2591、マージ済み) [一次・文書 / 一次・ソース: `config_template.cfg:1547`]。
- **FUN3D はリミッタを「edge 型」と「stencil 型」に分けている** (2010 講習資料):
  edge 型 (`minmod`/`vanleer`/`vanalbada`/`smooth`) は「散逸が小さく hex 格子ではよく働くが、混合要素・四面体では頑健でない」
  「**凍結できず、リミッタの循環で収束が止まることがある**」「随伴に使えない」。stencil 型 (`barth`/`venkat`/`h*`) は
  「頑健だが散逸的、全格子で動く」「**凍結でき**、凍結で収束が進むことがある」「随伴にはこちらが必須」
  [一次・文書: FUN3D 講習 2010 Session 7 p.5]。SU2 文書の「van Albada エッジ型は収束が頑健」とは評価が逆向きで、
  どちらも定量データを示していない。

---

## 3. 比較研究の結論 — 何が、どの基準で「最良」か

結論を先に書くと、**3 つの基準すべてで一位になるリミッタは文献上に無い**。各基準の勝者と、その根拠の強さは次の通り。

### 3.1 基準ごとの整理

| 基準 | 文献上の勝者 | 根拠 | 根拠の強さ |
| --- | --- | --- | --- |
| **滑らかな極値・滑らか域の精度** | Nishikawa $R_p$ (5 次まで保持) > Michalak–Ollivier-Gooch (4 次まで) > Venkatakrishnan (2 次まで) > BJ・minmod 系 | Nishikawa 2022 の理論 + 1D/2D 数値、Ahmad ら 2023 で「一貫して最小散逸」 | 一次・文書 (アブストラクト / 本文)。node 中心エッジベースで検証されている点は forge と同じ |
| **衝撃での頑健性・単調性** | 厳密な局所最大値原理を守るもの (BJ、MLP 系、h 系の圧力補助リミッタ)。Venkatakrishnan と R_p は単調性を厳密には保証しない | Venkatakrishnan は「厳密に単調ではない」(MLP-pw 論文の整理)、FUN3D は極超音速で h 系を推奨、MLP-pw は極超音速での安定性改善を主張 | 一次・文書 (講習資料・アブストラクト)。直接比較の定量はコード・問題依存 |
| **定常収束** | 微分可能かつ 1 を超えない関数 (R_p) ≥ Venkatakrishnan ≫ BJ/minmod。ただし強い衝撃では**どれも停滞しうる** | Venkatakrishnan 1993/1995、Ahmad ら 2023 (超音速バンプで R5 が 1e-12、Venkat も妥当、minmod/van Leer/van Albada は ~3000 反復で停滞)、Oliveira–Azevedo 2026 (V/W/R3 とも強衝撃ケースで停滞)、SU2 保守者の BJ 評価 | 一次・文書。「微分可能性が効く」は一致、「何で停滞が完全に消えるか」は未解決 |

### 3.2 個別の比較研究

- **Ahmad, Park, Nishikawa, Wang, Elmiligui (NASA, AIAA 2023-4405)** — FUN3D (node 中心、エッジ中点再構成、LSQ 勾配) で
  Nishikawa R5・Venkatakrishnan・van Albada・van Leer・minmod を比較 [一次・文書: NTRS 20230004071 本文]。
  - **実装の枠組み**: $\phi_j=\min_i\{\min(1,\phi_{ij})\}$、$\Delta_-=u_{ij}-u_j$ (エッジ中点への射影)、$\Delta_+=(u_\max-u_j)/2$。
    「$\tfrac12$ は古典的リミッタ関数 (minmod 等) をこの枠組みで使えるようにするため」で、そうした古典関数は
    「この枠組み用に設計されておらず、Venkatakrishnan より**かなり散逸的**」。ε は $\{K(6V_j/\pi)^{1/3}\}^3$、$K=1$ (バンプのみ 0.5)。
  - **R5 は上から 1 で抑えられるので外側の $\min(1,\cdot)$ を省ける** → その分の微分不可能点が消える。
  - **超音速バンプ (Roe、凍結なし)**: R5 は 1e-12 まで収束、Venkatakrishnan も「妥当な収束」、minmod/van Leer/van Albada は約 3000 反復で停滞。
    Venkatakrishnan の密度残差場は R5 より「ノイジー」。
  - **C25F (4.7M 節点、SA-neg、1500 反復で凍結)**: 収束が最も良かったのは**van Albada**、R5 と Venkat も良好。R5 が最も強い衝撃を出す。
  - **X-59 C612a (173.5M 節点、U-MUSCL κ=0.5、凍結なし)**: R5 と Venkat が良好、次いで van Albada。
  - **結論**: 「R5 は一貫して最小散逸、定常の反復収束は他と同等か良い」。**最良の収束リミッタは問題で入れ替わる** (C25F では van Albada)。
  - 背景として、ソニックブーム予測ワークショップ第 3 回の参加者は BJ と Venkatakrishnan で衝撃強さ・位置を妥当に解けたが
    「一貫して正確かつ頑健なリミッタは特定できなかった」と書いている [一次・文書: 同論文の緒言]。
- **Oliveira & Azevedo (arXiv 2601.16291, 2026)** — 遷音速 NACA0012、RANS-SA-neg、陰的。Venkatakrishnan・Wang・Nishikawa R3 の 3 者は
  「制御定数を適切な範囲に保てば工学的に同等、散逸特性は違う」、強衝撃ケースでは 3 者とも停滞 [一次・文書: アブストラクト。
  本文の詳細は前回調査 §2.2]。
- **Michalak & Ollivier-Gooch (2008/2009)** — 高次 (4 次) で精度を保つ微分可能リミッタと滑らか域の判定。「2 次の散逸低減にも
  収束性をほとんど犠牲にせず使える」[一次・文書: アブストラクト。本文未読 — 前回調査 §2.3]。
- **Venkatakrishnan (1993/1995)** — BJ 型は「数桁で停滞」、微分可能な修正で収束。May & Berger (SISC 2013) はこれを
  「単調性と定常収束の衝突」と要約し、自分たちの LP リミッタ (勾配の x・y 成分を別々に制限するベクトルリミッタ) でも
  **chatter は残り、スカラーより遅れて始まるだけ**、**$10^{-3}$ 程度の相対オーバーシュートを許すと収束が改善するが新しい極値を許す**と報告
  [一次・文書: May & Berger 論文 §5.1]。また同論文の表 5.1 では、近傍重心への再構成の方がエッジ中点への再構成より誤差が小さい
  (三角形の滑らかな定常ケース) [一次・文書]。
- **MLP (Park–Yoon–Kim JCP 2010、Park–Kim C&F 2012)** — 「MLP 条件は高次精度と**収束性の改善**を、偽振動なしで保証できる」と主張
  [一次・文書: JCP 2010 アブストラクトの検索要約。原文アブストラクトは出版社サイトが拒否、二次寄り]。
  後続の Zhang ら (C&F 2017) は「MLP は精度・頑健性・収束で良い結果を示してきたが、極超音速では安定性と収束に改善の余地がある」と整理
  [一次・文書: arXiv 1710.07187 緒言]。**MLP と R_p/Venkat を同じコードで比べた定量比較は今回見つけていない** [未確認]。
- **Tsoutsanis (JCP 2018)** — MOG リミッタの「上下限を何の近傍から取るか」だけを変え、再構成ステンシル全体から取る
  extended bounds で、セル近傍・頂点近傍の上下限より精度とメッシュ感度が良い [一次・文書: アブストラクト]。収束への言及はアブストラクトに無い。
- **Kitamura–Hashimoto の post limiter (JSFM 2016、JCP 2017、AIAA J 2018)** — 制限なし/制限ありの 2 候補を混ぜる a posteriori 型で、
  1D で約 4 倍・2D で約 16 倍相当の解像度、**NACA0012 (M 0.8) で密度残差の収束向上**を図示 [一次・文書: JSFM 2016 本文 Fig. 2 と緒言]。
  同論文の緒言は「滑らかな格子の遷音速 CFD ではリミッタなしでも破綻しない場合がある (Nishikawa 私信)」とも書く。
- **Diskin・Thomas らの NASA 比較** — node 中心と cell 中心の流束離散化を比べた一連の論文はあるが、**リミッタ同士を比べたものは見つからなかった**
  [未確認]。Jawahar & Kamath (JCP 2000) の重み付き勾配型リミッタ、Li & Ren 系の WBAP (成分ごとの重み付き偏り平均) は
  「精度・頑健性・収束が良い」と二次資料にあるが原典未読 [二次: arXiv 1710.07187 緒言]。

### 3.3 まとめ

- **精度だけなら R_p (特に R5) が現時点の最上位**で、NASA の実機級計算と FUN3D・SU2・scFLOW (2025) への実装がある。
- **強い衝撃の頑健性は「厳密な最大値原理 + 圧力センサ」側** (h 系、MLP-pw) が推奨されている。FUN3D が極超音速で推すのは R_p でなく `hvanalbada`。
- **定常収束は「微分可能 + 1 を超えない + (必要なら) 凍結」が実務の答え**で、それでも強衝撃では停滞が残るというのが 2023〜2026 の報告の一致点。
  停滞の原因について、Oliveira–Azevedo は「近似ヤコビアンがリミッタの非線形性を捉えない」を仮説にとどめ、Nishikawa ら自身は
  「ε を無くして滑らかな極値と一様域を別の方法で検出したい」「凍結を頑健に発火させたい」を未解決課題に挙げる (前回調査 §2.4)。

---

## 4. 最新の研究 (2018–2026)

各項目に「アイデア / 式の骨子 / 主張 / 成熟度 / node 中心エッジベースへの実装コスト」を書く。
コストは forge の構造 (ノード並列で近傍 min/max → 面ごと ψ → ノード min、を既に持つ) を前提にした見立てで、未検証。

### 4.1 Nishikawa の $R_p$ (AIAA 2022-1374) とその普及

- **アイデア**: Venkatakrishnan と同じ「ノードの近傍 min/max + エッジ中点の射影」の枠組みのまま、関数だけを差し替える。
- **式** (R5、FUN3D 実装の記法、$a=|2\Delta_+|$、$b=|\Delta_-|$):
  $\phi=\dfrac{a^5+\epsilon_5+aS_5}{a^5+\epsilon_5+b(a^4+S_5)}$ ($a<2b$)、$\phi=1$ (それ以外)、
  $S_5=8b^2[a^2-2b(a-b)]$、$\epsilon_5=\{K(6V_j/\pi)^{1/3}\}^6$ [一次・文書: Ahmad ら 2023 式 (12)–(14)]。
  SU2 の R3/R4/R5 は $\epsilon_p=(K\cdot\mathrm{RefElemLength})^{p+1}$ (全域定数) [一次・ソース / 一次・文書 (SU2 docs)]。
- **主張**: 5 次まで保持、1D でも完全に微分可能、ε 修正と両立、Venkat より低散逸、**無制限の再構成が有界ならそのまま保つ** (Venkat は保たない)
  [一次・文書: アブストラクト]。
- **成熟度**: FUN3D v14.0 (2023-02)、SU2 (`NISHIKAWA_R3/4/5`)、scFLOW (MSC/Cradle、2025 年の論文で「改良したリミッタ関数」を挙げる。
  どの関数かはアブストラクトで未確認) [一次・文書]。SU2 Conference 2023 で SU2 実装の評価発表 (Sachdeva) があるが本文未読 [二次]。
- **コスト**: forge では `limiterFunction` の差し替え + ε の次数 ($p+1$ 乗) のみ。**小**。
- **極値節点での振る舞い** (forge の問題に直結): $a=0$ (節点が極値、$\Delta_+=0$) のとき $S_5=16b^4$ なので
  $\phi=\epsilon_5/(\epsilon_5+16b^5)$ [一次・文書から導出]。Venkatakrishnan の $\epsilon^2/(\epsilon^2+2\Delta_-^2)$ と同じ型で、
  **$b\gg\epsilon_5^{1/5}$ なら ψ≈0、$b\ll$ なら ψ≈1 の二値化は R_p でも残る**。R_p が消すのは外側の $\min(1,\cdot)$ の折れ目であって、
  「節点が極値か否か」で ψ が跳ぶ機構ではない。

### 4.2 MLP 系 (Park–Kim 2010/2012、Zhang ら 2017・2023)

- **アイデア**: セル中心で、各セル**頂点**での再構成値を、その頂点を共有する全セルの min/max に収める (MLP 条件
  $q^{\min}_{V(l)}\le q_l\le q^{\max}_{V(l)}$) [一次・文書: arXiv 1710.07187 式 (18)–(19)]。面中心でなく頂点で縛るので
  BJ (面中心評価) より厳しいが、上下限を頂点共有の広いステンシルから取るので「多次元性」を持つ。
  MLP-u1/u2 は制限の掛け方の違い (u2 は Venkat 型の微分可能関数を使う、と二次資料) [二次]。
- **MLP-pw (Zhang ら C&F 2017)**: 緩い weak-MLP (滑らか域で散逸減) と厳しい strict-MLP (衝撃で安定) を、**微分可能な圧力重み関数**で混ぜる。
  「極超音速で安定性と収束が改善、圧力変化の無い所では低散逸で接触・膨張を正確に捉える」[一次・文書: アブストラクト]。
- **Zhang・Yuan・Liu (arXiv 2311.17240, 2023)**: 衝撃安定性の議論。リミッタが衝撃近傍で大きな数値誤差を持ち込みうること、
  圧力ベースの衝撃指標で散逸を調整することを提案 [一次・文書: アブストラクト]。
- **成熟度**: 韓国・中国の研究コードで使用。商用・NASA コードでの採用は確認できない [未確認]。
  (Fluent の "multidimensional limiter" は Makarov を参照する別物で、MLP ではない [一次・文書: Fluent 理論ガイド]。)
- **node 中心への写し方**: median dual の「頂点」は辺中点・面重心・要素重心で、それを縛る上下限は「その要素の節点」の min/max になる。
  **node 中心版 MLP の定式化は文献で見つけていない** [未確認]。forge に入れるなら設計から始める必要がある (コスト**中〜大**)。

### 4.3 上下限の取り方を広げる (extended bounds、頂点共有ステンシル)

- **Tsoutsanis (JCP 2018)**: 上下限を再構成ステンシル全体から取ると精度・メッシュ感度が改善 [一次・文書]。MLP の頂点共有ステンシルも、
  Wind-US の `DEBUG 64` (face 近傍 → 頂点共有セル) も同じ方向 (前回調査 §2.6)。
- **forge への含意**: forge の近傍は「内部双対面を共有する実ノード」= 辺でつながる節点だけ。四角形・六面体では、要素を共有するが辺でつながらない
  対角節点が上下限に入っていない。対角を入れると「節点が近傍の極値」になる頻度が下がる [推測]。コスト**中** (節点→要素の接続が要る)。
  単調性は弱まる (Wind-US の注記: TVD は保証しない)。

### 4.4 圧力センサ・第 2 リミッタ (滑らか域でリミッタを切る)

前回調査 §2.7 で FUN3D h 系と Kitamura–Shima を挙げた。ここでは 2018 年以降の動きだけ書く。

- **FUN3D v14.2 (2025-04) の単独 `hpressure`** — h 系の圧力補助リミッタを単独のリミッタとして使える [一次・文書: リリースノート]。
  「圧力ベースだけで制限する」選択肢が公式に追加されたことになる。中身の式は今回の資料に無い [未確認]。
- **Gnoffo の補助リミッタの使われ方**: Kitamura–Hashimoto の post limiter は Gnoffo の補助リミッタ $G$ を
  $(\alpha_\max,\alpha_\min)=(3,2)$ で使い、**$p_\max/p_\min\ge3$ なら $G=0$、$\le2$ なら $G=1$**、その間を滑らかにつなぐ
  [一次・文書: JSFM 2016 式 (4)。補間関数の形は抽出で式が欠けており未確認]。FUN3D の h 系は同じ圧力リミッタを全変数の ψ に掛ける (前回調査)。
- **Kitamura–Shima (AIAA J 2012) の second limiter**: tanh で淀み・亜音速域のリミッタを切る、パラメタなし [二次: 1710.07187 緒言]。式は未確認。
- **成熟度**: FUN3D (h 系は推奨設定)、JAXA FaSTAR (Kitamura 系) [一次・文書 / 二次]。
- **コスト**: forge は近傍 min/max のループを既に持つので、$p_\max/p_\min$ からノードごとの係数を作り
  $\psi\leftarrow 1-G_p(1-\psi)$ 等で緩める形なら**小**。ただし「どう混ぜるか」(h 系は掛け算、post limiter は混合) は出典ごとに違う。
- **限界**: 接触不連続・せん断層・温度/化学種の段差は圧力比で検出できない (Nishikawa–White–O'Connell の指摘、前回調査 §2.7)。

### 4.5 a posteriori 型 (MOOD、post limiter)

- **MOOD (Clain–Diot–Loubère JCP 2011〜)**: 無制限の候補解を作り、正値性・離散最大値原理 (DMP) を破るセルだけ低次に落として再計算。
  高次 (WENO 等) 向けで、セルごとに次数が変わり反復が要るので「並列・陰解法に不向き」と Kitamura–Hashimoto が指摘
  [一次・文書: JSFM 2016 §3]。2018 年以降は Tsoutsanis らの relaxed MOOD (CMAME 2020、AMC 2022) など**高次 FV** 側で発展 [一次・文書: 書誌]。
- **Kitamura–Hashimoto の post limiter** (2 次専用に簡略化した a posteriori):
  1. 同時に作れる 2 候補だけを使う — 無制限の面値 $q_\mathrm{unlim}$ と制限付き $q_\mathrm{lim}$。
  2. 面ごとに、$\rho_\mathrm{unlim}>0$ かつ $p_\mathrm{unlim}>0$ (正値性)、かつ $\rho$・$p$ の面値が**エッジ両端値の間** (面の DMP) なら $G'=1$、でなければ $G'=0$。
  3. $G=G'\cdot$ (Gnoffo の圧力比リミッタ)、$q=G\,q_\mathrm{unlim}+(1-G)\,q_\mathrm{lim}$ で混ぜる。
  4. 非構造三角形では、セル中心→面中心→隣接セル中心の折れ角 $\theta_\mathrm{face}$ による係数 $\psi_\mathrm{face}=0.25(\cos\theta+|\cos\theta|)^2$ で
     $G=\min(G',\psi_\mathrm{face})$ と弱める (一直線上なら 1、90° 超で 0) [一次・文書: JSFM 2016 式 (3)–(6)]。
  5. 3D 版 (AIAA J 2018) は「汚い (形の悪い) セル」を検出してそこではリミッタを残す [一次・文書: アブストラクト]。
  - **主張**: 2D で約 16 倍相当の解像度、NACA0012 で残差収束の向上、FaSTAR に実装 [一次・文書]。
  - **成熟度**: JAXA FaSTAR (研究〜実務)。
  - **コスト**: 判定はエッジ局所 (両端値と面値の比較) なので、**エッジ中点再構成の forge とは相性が良い**。実装は**小〜中**。
    ただし判定 2 は if による二値で、**判定の境界にいる面では反復ごとに $G'$ が 0/1 を行き来しうる** [推測]。
    定常収束への効果は 2D の 1 例の図のみで、陰解法・RANS での報告は見つけていない [未確認]。

### 4.6 エッジ型の滑らかリミッタの再評価

- 近傍の min/max を使わず、エッジ両端の差と各端の射影だけで ψ を決める型 (van Albada エッジ型)。SU2 の `VAN_ALBADA_EDGE` は U-MUSCL (κ) と組み合わせて
  2020 年代に追加・保守されている (§2.2) [一次・ソース]。FUN3D の `vanalbada` も講習資料上は edge 型 [一次・文書]。
- **forge への含意**: **ノードの近傍 min/max も、面ごとの ψ の min も無い**ので、「節点が極値か否か」「どの面が ψ を決めるか」で ψ が跳ぶ機構そのものが無い
  [一次・ソースから導出: SU2 の式は $\delta_{ij}$ と $\mathrm{proj}$ の滑らかな関数]。
  代わりに局所最大値原理の保証が無く、Nishikawa 型の枠組みに当てはめた古典関数は「かなり散逸的」(Ahmad ら)、FUN3D は「四面体・混合要素で頑健でない」
  「凍結できず循環で収束が止まることがある」とする [一次・文書]。SU2 文書は逆に「収束は頑健」とする。**評価が割れている**。
- **コスト**: forge のエッジ中点再構成にそのまま乗る。**小**。ε (SU2 は次元付き $10^{-6}$) は forge では $q_\mathrm{ref}$ で無次元化した値にする必要がある。

### 4.7 微分可能性・随伴・Newton 向け

- 随伴・Newton 系で微分可能リミッタが要るという要求は変わっていない (FUN3D: 随伴には stencil 型が必須、SU2: 離散随伴用の凍結
  `FROZEN_LIMITER_DISC`) [一次・文書 / 一次・ソース]。
- **JAX-FVM (arXiv 2607.07385, 2026)**: 非構造三角形の 2D 圧縮性 FV を全体 AD 可能に書き、リミッタは Venkatakrishnan [一次・文書: アブストラクト]。
  新しいコードでも「微分可能な既定」として Venkatakrishnan が選ばれている例。
- 「リミッタを Newton で線形化して停滞を消す」系統の 2018 年以降の決定的な結果は見つけていない。Mor-Yossef (JCP 2024) は SA 方程式の
  正値性保持の準 Newton で、リミッタの話ではない [一次・文書: 書誌]。

### 4.8 機械学習リミッタ

- **Nguyen(-Fotiadis)・McKerns・Sornborger (PoF 2022)**: 粗視化した 1D Burgers で、高解像度データとの誤差を最小化する区分線形の flux limiter を学習。
  「$\phi(1)=1$ の規則に縛られないのに標準リミッタより誤差が小さい」[一次・文書: アブストラクト]。
- **Probabilistic flux limiters (PoF 2025、arXiv 2405.08185)**: 複数のリミッタ関数を確率付きで選ぶ。同じく 1D Burgers [一次・文書]。
- **成熟度**: 1D・構造格子・非定常の研究段階。**非構造・定常 RANS・陰解法での使用例は見つけていない** [未確認]。
  学習した関数が微分可能・有界である保証も無いので、forge の問題 (定常収束) には現時点で当てはまらない。

### 4.9 正値性・エントロピー

- 2 次 FV で正値性を保つ枠組み (invariant-domain preserving、MUSCL の正値性条件) は理論側で進んでいるが、
  実用コードの対処は依然「非物理な再構成の面だけ 1 次に落とす」(SU2 `bad_recon`、FUN3D「外挿失敗点だけ再計算」) [前回調査 §2.6・§2.5]。
  forge には `badReconFallback` がある ([`methods/limiter.md`](../../methods/limiter.md))。Kitamura の post limiter も正値性判定を含む (§4.5)。
- エントロピー安定 FV (JAX-FVM 等) はリミッタの代わりではなく流束側の話 [一次・文書]。

### 4.10 低マッハ対応

- 低マッハ・淀み域で不要な制限を外す、という意味では §4.4 の Kitamura–Shima (亜音速・淀みで切る) と MLP-pw (圧力変化の無い所で緩める) が該当。
- 前処理 (低マッハ) と組み合わせたリミッタの専用設計は今回見つけていない [未確認]。forge の入口列 (M 0.11 程度) はこの「圧力がほぼ一様な亜音速域」に当たる。

---

## 5. forge への適合性と候補の比較

### 5.1 forge 側の事実の整理 (候補を評価する軸)

1. **ψ の跳びは「節点が極値か否か」と「どの面が min を取るか」の 2 つの離散的な切替から来ている**
   (前回調査 §3.1、plan [`limiter-inlet-column-oscillation`](../../plans/accepted/limiter-inlet-column-oscillation.md) #4r)。
2. **ε を領域一定 2e-6 にしても床は約 0.65 倍** → 増分 $\hat\Delta_-$ ($10^{-4}$〜$10^{-3}$) が許容オーバーシュートで縛った ε より
   まだ桁で大きく、**ε の大きさだけで切替を消すには、単調性を許容値以上に緩める必要がある**と読める [推測]。
   May & Berger の「$10^{-3}$ のオーバーシュートを許すと収束は改善するが新しい極値を許す」と同じ二律背反。
3. forge は**エッジ中点再構成**・node 中心・float32・陰的定常で、FUN3D・SU2・TAU と同じ系統。§4 のうち FUN3D/SU2 で動いている手法は移植の前例がある。
4. 床は少数節点 (入口列・壁の折れ点・出口中心線) に集中し、そこは**圧力がほぼ一様** (入口列) か**幾何の角**。
5. 衝撃 (超音速ノズル・噴流) と境界層の両方が要る。境界層内では、法線方向に単調な速度分布で ψ が下がると壁摩擦・熱流束に効く
   (FUN3D の `--limit_near_walls_less` の注記 — 前回調査 §2.6)。

### 5.2 候補の比較表

効果は §5.1 の機構からの見込みで、いずれも forge では未検証。

| 候補 | 何が変わるか | 実装コスト | 床への期待 (機構) | 精度・衝撃のリスク | 実用の前例 |
| --- | --- | --- | --- | --- | --- |
| **A. エッジ型の滑らかリミッタ** (van Albada エッジ型、SU2 `VAN_ALBADA_EDGE` 相当) | 近傍 min/max と面 min を使わず、エッジ両端差と射影だけで ψ | 小 (forge のエッジ中点再構成に直接乗る。ε の無次元化が要る) | **切替機構 1 の両方 (極値判定・面 min) が無くなる**。ψ は $\delta_{ij}$ と $\mathrm{proj}$ の滑らかな関数 | 局所最大値原理なし。強い衝撃・四面体での頑健性に FUN3D は否定的、SU2 は肯定的。滑らかな極値では $\delta_{ij}\approx0$ で ψ が下がる (切れ目なく)。境界層への影響は未知 | SU2 (流れのみ)、FUN3D `vanalbada`。C25F で最良の収束 (Ahmad ら) |
| **B. 圧力センサで緩める** (h 系 / `hpressure` / Kitamura–Shima 型) | $p_\max/p_\min$ が小さい節点で ψ を 1 側へ寄せる | 小 (近傍ループで $p$ の min/max を取れば済む) | 入口列のような圧力一様域では ψ が 1 に固定され、切替が残っても流束に効かない | 接触・せん断層・温度/化学種段差を見逃す。混ぜ方 (積・混合) と閾値の選択が要る | FUN3D の推奨 (`hvanalbada`)、FaSTAR |
| **C. Nishikawa $R_3$/$R_5$ に差し替え** | 関数のみ。外側 $\min(1,\cdot)$ の折れ目が消え、有界な再構成は保持 | 小 (既存の R1 は別物。$R_p$ を新規に) | **極値節点の二値化は残る** (§4.1 の導出)。面 min の切替も残る。床への効果は小さい見込み | 低散逸側 (衝撃が鋭くなる)。Sod・近傍逸脱ゲートの再測定が要る | FUN3D v14、SU2、scFLOW、NASA 実機級計算 |
| **D. ψ の凍結** (前回調査) | N 反復後に ψ 固定、再構成失敗点だけ再計算 | 小 | 切替そのものを止める | 凍結時点への依存。`limiter-config-simplify` §4.3 の決定の見直しが要る | SU2・FUN3D・Wind-US、NASA の常用 |
| **E. 上下限を要素共有 (対角) 節点まで広げる** (extended bounds / MLP の頂点ステンシル相当) | 極値判定の母集団が増える | 中 (節点→要素接続) | 「節点が近傍の極値」の頻度が下がる。面 min の切替は残る | 単調性が弱まる (Wind-US 注記) | Tsoutsanis 2018 (高次 FV)、Wind-US `DEBUG 64`、MLP |
| **F. post limiter** (面ごとの正値性 + 面 DMP + 圧力比で無制限と制限付きを混合) | 制限が要らない面では無制限値を使う | 小〜中 (エッジ局所判定) | 入口列の面で DMP が満たされれば ψ の値に関係なく無制限側になる。ただし判定が二値で、境界上の面は 0/1 を行き来しうる | 解像度は上がる (2D で約 16 倍相当)。極端な格子では 3D 版の汚いセル判定が要る | FaSTAR (JAXA) |
| **G. node 中心版 MLP** | 双対 CV の頂点で再構成値を要素節点の min/max に収める | 中〜大 (定式化から) | 未知 | 単調性は強い側 | node 中心の前例なし |
| **H. 機械学習リミッタ** | 学習した関数 | 大 | 未知 (微分可能性・有界性の保証なし) | 未知 | 1D 研究段階のみ |

### 5.3 順位 (理由つき。決定ではない)

1. **A. エッジ型の滑らかリミッタ** — forge の観測 (ψ の切替が「極値判定」と「面 min」という離散選択から来る) に対して、
   **その離散選択を構造ごと持たない唯一の候補**。forge のエッジ中点再構成に最も素直に乗り、SU2 の実装 (U-MUSCL と一体) が手元のソースにある
   (`common.hpp:147–176`)。弱点は衝撃での単調性保証の欠如と、FUN3D と SU2 で評価が割れている点で、超音速ノズル・噴流で
   近傍逸脱 (`limiterDiag`) と Sod の再測定が必須。ε の無次元化 (SU2 の $10^{-6}$ は SI の次元付き) を forge 流に直す必要がある。
2. **B. 圧力センサで緩める** — 床が集中する入口列は M 0.11 程度の圧力ほぼ一様域なので、機構上は「切替は残っても流束に効かない」状態にできる。
   NASA FUN3D が衝撃問題の**推奨既定**として長年使っている (h 系) ので、強い衝撃側のリスクが最も小さい候補。
   実装も小さい。弱点は接触・せん断層・温度/化学種段差を見逃すことと、幾何の角・出口中心線の節点が圧力一様域とは限らないこと
   (そこには効かない可能性)。現行 Venkatakrishnan の上に重ねる形で A/B ができる。
3. **C. Nishikawa $R_p$** — 文献上は「精度と収束の両方で Venkat 以上」の最も強い証拠があり、FUN3D・SU2・scFLOW と実装も揃う。
   ただし本ノートの導出 (§4.1) では、**forge の床の機構 (極値節点の二値化と面 min) は R_p でもそのまま残る**ので、
   床を下げる期待は小さい。精度 (衝撃の鋭さ・滑らかな極値) を上げる目的なら最有力。

D (凍結) は前回調査で既に整理済みで、A〜C のどれとも併用できる補完策という位置づけ。E と F は機構上は効きうるが、
E は単調性を弱め、F は判定の二値性が新しい切替を作りうるので A・B より後ろに置いた。G・H は前例・成熟度が足りない。

---

## 6. 未確認事項

1. **FUN3D 現行 (v14.x) マニュアルの `flux_limiter` 既定値**と `smooth` リミッタの中身 (2014 版マニュアルで既定 `none` を確認したのみ)。`hpressure` の式。
2. **USM3D・Kestrel・Loci/CHEM・STAR-CCM+・CONVERGE・elsA・OVERFLOW の既定リミッタ** (公式文書に当たれなかった)。
3. **Code_Saturne の `isstpc` 既定値** (版で記述が割れている)。
4. **Kitamura–Shima (2012) second limiter の式** と、Gnoffo 圧力リミッタの補間関数の正確な形 (JSFM 2016 の式 (4) は抽出で欠落)。
5. **MLP-u1/u2 の定義の違い**と、MLP と R_p/Venkat を同じコードで比べた定量比較の有無。
6. **node 中心 median dual 版 MLP** の定式化の有無。
7. **edge 型 van Albada の収束評価が FUN3D (否定的) と SU2 (肯定的) で割れている理由** — 格子種別 (hex か tet か)、
   κ、凍結の有無のどれが効いているか。
8. **scFLOW 2025 論文の「改良リミッタ関数」が何か** (アブストラクトのみ)。SU2 Conference 2023 の R_p 評価の本文。
9. **SU2 の BJ 削除の議論の続き** (PR #2911 は閉じられた。Issue #1841 の決着、文書の更新の有無)。
10. forge 側: 本ノート §4.1 の「R_p でも極値節点の二値化は残る」は式からの導出で、forge の入口列で確かめていない。
    A (エッジ型) で入口列の ψ が実際に滑らかになるか、Sod・近傍逸脱・超音速ノズルでどれだけ悪化するかも未測定。

## 出典

- SU2 ソース (`/home/sano/work/forge/.external/su2-src`, commit 12eb826f): `Common/src/CConfig.cpp` (2004–2103 行のリミッタ既定)、
  `config_template.cfg` (1546–1583 行)、`SU2_CFD/include/limiters/CLimiterDetails.hpp` (`BARTH_JESPERSEN` 特殊化)、
  `SU2_CFD/include/numerics_simd/flow/convection/common.hpp` (`musclEdgeLimited`)。
  SU2 文書 [Slope Limiters and Shock Resolution](https://su2code.github.io/docs_v7/Slope-Limiters-and-Shock-Resolution/)。
  GitHub [PR #2911](https://github.com/su2code/SU2/pull/2911) (BJ 削除案、未マージで閉鎖)、PR #2908 (`VAN_ALBADA_EDGE` の乱流・化学種、未マージ)、
  PR #2591 (MUSCL κ 法、2025-10 マージ)。
- FUN3D: [Manual Chapter 6 (2014-06-05 版)](https://fun3d.larc.nasa.gov/papers/Website_June2014_Chapter06_Analysis.pdf)、
  [Chapter 1 リリースノート](https://fun3d.larc.nasa.gov/chapter-1.html)、
  [講習 2010 Session 7 (Supersonic/Hypersonic)](https://fun3d.larc.nasa.gov/session7.pdf)、
  [講習 2023 Session 6 (Thermochemical Nonequilibrium)](https://fun3d.larc.nasa.gov/session6_2023.pdf)。
- CFL3D: Krist ら, *CFL3D User's Manual (Version 5.0)*, NASA/TM-1998-208444 ([NTRS 19980218172](https://ntrs.nasa.gov/citations/19980218172))、
  [CFL3D v6 New Features](https://nasa.github.io/CFL3D/Cfl3dv6/cfl3dv6_new.html)。
- OVERFLOW: Nichols, *OVERFLOW 2 Training Class Morning Session*, 2010 ([NASA](https://www.nasa.gov/wp-content/uploads/2025/09/overflow-training2010-morning.pdf))。
- DLR TAU: Kirz & Rudnik, *DLR TAU Simulations for the Second AIAA Sonic Boom Prediction Workshop* ([発表資料](https://lbpw-ftp.larc.nasa.gov/sbpw2/workshop/sbpw2-talks-nearfield/02-sbpw2-dlr-kirz-rudnik.pdf))。
- Ansys Fluent Theory Guide, Gradient Limiters ([v24.2](https://ansyshelp.ansys.com/public////Views/Secured/corp/v242/en/flu_th/flu_th_sec_grad_limit.html))。
- Cart3D: [input.cntl 解説](https://www.nas.nasa.gov/publications/software/docs/cart3d/pages/input_cntl.html)。
- Ashton ら, *3rd High-Lift Workshop Summary Paper — OpenFOAM, STAR-CCM+ & LAVA* ([ResearchGate](https://www.researchgate.net/publication/322312686)) — 二次 (検索要約)。
- CONVERGE: [CFD-Wiki CONVERGE FAQ](https://www.cfd-online.com/Wiki/CONVERGE_FAQ) — 二次。Code_Saturne: [Theory Guide 7.1](https://www.code-saturne.org/documentation/7.1/theory.pdf) — 二次 (検索要約)。
- elsA: CERFACS TR/CFD/11/64 *Development of a new hybrid compressible solver inside the CFD elsA software* — 二次 (検索要約)。
- Nishikawa, H., "New Unstructured-Grid Limiter Functions," AIAA 2022-1374. doi:10.2514/6.2022-1374
- Ahmad, N., Park, M.A., Nishikawa, H., Wang, L., Elmiligui, A., "Evaluation of Limiter Functions for Supersonic Applications," AIAA 2023-4405. [NTRS 20230004071](https://ntrs.nasa.gov/api/citations/20230004071/downloads/aiaa-2023-cst-Ver9.pdf)
- Nakashima, Y., Higo, Y., Nishikawa, H., "Recent Algorithmic Advances in the scFLOW Unstructured-Polyhedral-Grid CFD Solver," AIAA 2025-0075. doi:10.2514/6.2025-0075
- Oliveira, F.B., Azevedo, J.L.F., "A Study of Improved Limiter Formulations for Second-Order Finite Volume Schemes Applied to Unstructured Grids," [arXiv:2601.16291](https://arxiv.org/abs/2601.16291), 2026.
- Michalak, C., Ollivier-Gooch, C., JCP 228, 8693–8711, 2009. doi:10.1016/j.jcp.2009.08.021
- Park, J.S., Yoon, S.-H., Kim, C., "Multi-dimensional limiting process for hyperbolic conservation laws on unstructured grids," JCP 229, 788–812, 2010. doi:10.1016/j.jcp.2009.10.011
- Park, J.S., Kim, C., "Multi-dimensional limiting process for finite volume methods on unstructured grids," Comput. Fluids, 2012. doi:10.1016/j.compfluid.2012.04.015
- Zhang, F., Liu, J., Chen, B., "Modified multi-dimensional limiting process with enhanced shock stability on unstructured grids," Comput. Fluids, 2017. doi:10.1016/j.compfluid.2017.11.019 ([arXiv:1710.07187](https://arxiv.org/abs/1710.07187))
- Zhang, F., Yuan, Z., Liu, J., "A discussion on numerical shock stability of unstructured finite volume method: Riemann solvers and limiters," [arXiv:2311.17240](https://arxiv.org/abs/2311.17240), 2023.
- Tsoutsanis, P., "Extended bounds limiter for high-order finite-volume schemes on unstructured meshes," JCP, 2018. doi:10.1016/j.jcp.2018.02.009
- Kitamura, K., Hashimoto, A., 「高解像度流体計算に向けた a posteriori 制限関数 (第 3 報): 非構造格子への拡張」, 第 30 回数値流体力学シンポジウム B02-2, 2016 ([PDF](https://www2.nagare.or.jp/cfd/cfd30/web_paper/B02-2.pdf))。
- Kitamura, K., Hashimoto, A., "Simple a posteriori slope limiter (Post Limiter) for high resolution and efficient flow computations," JCP 341, 313–340, 2017. doi:10.1016/j.jcp.2017.04.002
- Kitamura, K., Aogaki, T., Inatomi, A., Fukumoto, K., Takahama, T., "Postlimiters and Simple Dirty-Cell Detection for Three-Dimensional Unstructured (Unlimited) Aerodynamic Simulations," AIAA J, 2018. doi:10.2514/1.J056683
- Kitamura, K., Shima, E., "Simple and Parameter-Free Second Slope Limiter for Unstructured Grid Aerodynamic Simulations," AIAA J 50(6), 1415–1426, 2012. doi:10.2514/1.J051269
- May, S., Berger, M., "Two-Dimensional Slope Limiters for Finite Volume Schemes on Non-Coordinate-Aligned Meshes," SIAM J. Sci. Comput., 2013 ([著者版 PDF](https://cs.nyu.edu/~berger/lpLimiter_revised.pdf))。
- Nguyen, N.T.T., McKerns, M.M., Sornborger, A.T., "Machine learning changes the rules for flux limiters," Phys. Fluids 34, 085136, 2022. doi:10.1063/5.0102939
- Nguyen-Fotiadis, N.T.T. ら, "Probabilistic flux limiters," Phys. Fluids 37, 046112, 2025 ([arXiv:2405.08185](https://arxiv.org/abs/2405.08185))。
- de Romémont, G., "JAX-FVM: A differentiable, entropy-stable finite volume solver on unstructured meshes for compressible flows," [arXiv:2607.07385](https://arxiv.org/abs/2607.07385), 2026.
- Farmakis, P.S., Tsoutsanis, P., Nogueira, X., "WENO schemes on unstructured meshes using a relaxed a posteriori MOOD limiting approach," CMAME, 2020. doi:10.1016/j.cma.2020.112921
- Mor-Yossef, Y., "Quasi-Newton positivity-preserving scheme for the Spalart-Allmaras turbulence model using unstructured grids," JCP, 2024. doi:10.1016/j.jcp.2024.113166
