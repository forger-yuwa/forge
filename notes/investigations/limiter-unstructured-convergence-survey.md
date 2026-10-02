# 非構造格子のリミッタによる定常収束の停滞 — 文献・実用コードの対処の調査

## メタ

- **area**: `limiter`
- **status**: `draft` (**調査専用 / コード変更なし / 推奨は決めていない**)
- **related_docs**:
  - [`methods/limiter.md`](../../methods/limiter.md) (無次元化 Venkatakrishnan `limiterScaled: 1` の現行仕様)
- **related_plans**:
  - [`limiter-inlet-column-oscillation.md`](../../plans/active/limiter-inlet-column-oscillation.md) §3.1・§5.1 #2r・#4r・#5er (本調査の動機: case/16 の入口列の床)
  - [`limiter-config-simplify.md`](../../plans/active/limiter-config-simplify.md) §4.3 (「ψ 凍結は今回は入れない」の決定。証拠しだいで見直し可)
- **created**: `2026-10-03`
- **一次資料の置き場**: SU2 ソースは `/home/sano/work/forge/.external/su2-src` (commit `12eb826f`, 2026-04-27)。
  ダウンロードした論文 PDF は作業用スクラッチにのみ置いた (リポジトリには入れていない)。

**確認度の表記** (各項目の末尾に付ける):

- **[一次・ソース]** — ソースコードを読んで確認した。
- **[一次・文書]** — 公式マニュアル・論文本文 (または公式アブストラクト) を読んで確認した。
- **[二次]** — 他論文・検索結果の要約など、一次資料を直接読めていない。
- **[未確認]** — 記憶・推測。採否の根拠に使わないこと。

---

## 1. 問い

forge (node 中心・median dual、2 次 MUSCL、SLAU、block-DPLUR 定常擬似時間、float32) の既定リミッタ
(無次元化 Venkatakrishnan: $\hat\delta=\delta/q_\mathrm{ref}$、$\hat\epsilon^2=(K h_i/L_\mathrm{ref})^3$、$K=0.05$、
ψ はエッジ中点評価の面ごと値の最小) で、定常残差が旧リミッタより約 20 倍高い床で止まる。
床は少数の節点 (no-slip 壁に隣接する入口境界の列、幾何の折れ点の壁節点、出口中心線) に集中し、そこでは

- 節点自身が近傍の極値 ($\delta^\pm=0$ または float32 の数 ulp)、
- $\hat\epsilon\approx$ $q/q_\mathrm{ref}$ の 1 ulp、
- ψ が反復ごとに 0 付近と 1 付近を行き来し、ψ を決める面が 0.85〜1.0/step で入れ替わる、
- CFL 半減は効かない (plan #2r の判定 B)、

ことが観測されている (plan §3.1・#4r・#5er)。**コミュニティはこの種の「リミッタ起因の定常収束停滞」を
どう扱っているか**を、文献と実用コード (SU2 / FUN3D / OpenFOAM / Fluent / Wind-US / CFL3D 等) について整理し、
forge の観測 (極値節点 + ulp 程度の ε̂ + 境界の角) に効きうる候補を並べる。**採否は決めない。**

---

## 2. 手法ごとの整理

### 2.0 一覧

| # | 手法 | 何を変えるか | 狙い | 収束への効果の証拠 | 確認度 |
| --- | --- | --- | --- | --- | --- |
| 1 | Venkatakrishnan (1993/1995) | Barth–Jespersen の min を微分可能関数に、$\epsilon^2=(K\Delta x)^3$ で一様域を不活性化 | 微分可能化 + 一様域で ψ→1 | 原論文: BJ 型は数桁で停滞、修正版で収束 | 一次・文書 (アブストラクト) |
| 2 | Wang (2000) / SU2 `VENKATAKRISHNAN_WANG` | $\epsilon=K'(q_\max^{gl}-q_\min^{gl})$ (全域レンジ) | 格子尺度・次元への依存を消す | 直接の収束比較は確認できず。2026 年の比較では V/W/R3 とも強衝撃ケースで停滞 | 一次・ソース (SU2) / 一次・文書 (比較論文) |
| 2' | SU2 `VENKATAKRISHNAN` | $\epsilon^2=(K\cdot\mathrm{RefElemLength})^3$、**RefElemLength=1 の定数** (セル尺度ではない) | — | — | 一次・ソース |
| 2'' | OpenFOAM `cellLimited<Venkatakrishnan> k` | ε 無しの比の形 + 係数 k で近傍レンジを $(1/k-1)(\max-\min)$ だけ広げる | 許容幅の相対緩和 | 記載なし | 一次・ソース |
| 3 | Michalak & Ollivier-Gooch (2008/2009) | 3 次多項式の微分可能関数 ($y_t$ で 1 に到達)、滑らか域の判定 | 高次精度の保持と収束の両立 | 2008 アブストラクト: 収束良好、2 次にも転用可 | 一次・文書 (アブストラクト)、関数形は OpenFOAM `cubic` で一次・ソース |
| 4 | Nishikawa R_p (2022) | $p=3,4,5$ の微分可能関数、1 を超えない、ε 修正と両立 | 低散逸・高次保持・微分可能 | NASA 2023: 超音速バンプで R5 が 1e-12 まで収束、minmod/vanLeer/vanAlbada は 3000 反復で停滞 | 一次・文書、関数は SU2 ソースで一次・ソース |
| 5a | SU2 `LIMITER_ITER` | 指定反復以降 ψ を再計算しない (流れ側は凍結値を使い続ける) | 停滞の解消 | 文書に効果の記述なし | 一次・ソース + 文書 |
| 5b | FUN3D `--freeze_limiter` / `freeze_limiter_iteration(s)` | 指定反復以降 ψ を凍結、再構成失敗点だけ再計算 | 「stall / ring」する収束の改善 | NASA の実計算で 400〜2000 反復凍結が常用、K=5 + 凍結で機械零 | 一次・文書 |
| 5c | Wind-US `DEBUG 57` | 非構造の slope limiter を mode 反復後に凍結 | 同上 | 記載なし | 一次・文書 |
| 6 | 境界の扱い | 境界値を極値に含める (OpenFOAM)、境界近傍で別係数 (Wind-US 構造格子 `BOUNDARY TVD`)、壁近傍で弱める (FUN3D `--limit_near_walls_less`)、SU2 は特別扱い無し | — | 収束への効果の記述は見つからず | 一次・ソース/文書 |
| 7 | 滑らか域でリミッタを切るセンサ | FUN3D h 系 (Gnoffo の圧力比リミッタ)、Kitamura–Shima の second limiter、Nejat–Ollivier-Gooch の tanh | 亜音速・淀み域の不要な制限を消す | Kitamura–Shima の詳細は未読 | 一次・文書 (FUN3D)、二次 (他) |
| 8 | 丸め・浮動小数点 | 一様域の機械零レベルの雑音でリミッタが入る、min の切替で勾配が回転して chatter | — | Berger ら: scalar min の実装で 3 桁停滞、recentering で解消 | 一次・文書 (Berger)、二次 (Nishikawa の記述) |

以下、項目ごとに式・解く問題・収束の証拠・出典を書く。

### 2.1 Venkatakrishnan (AIAA 93-0880 / JCP 118, 1995) と ε の役割

- **式** (forge・SU2 と同形。$\Delta_-$ = 再構成の増分 (エッジ中点への射影)、$\Delta_+$ = 近傍の最大/最小と自身の差):

  $$\psi=\frac{1}{\Delta_-}\,\frac{(\Delta_+^2+\epsilon^2)\Delta_-+2\Delta_-^2\Delta_+}{\Delta_+^2+2\Delta_-^2+\Delta_+\Delta_-+\epsilon^2},\qquad \epsilon^2=(K\Delta x)^3 .$$

  SU2 の実装は同値な形 `y = δ(δ+proj)+ε²; ψ = (y + δ·proj)/(y + 2 proj²)` [一次・ソース: `SU2_CFD/include/limiters/CLimiterDetails.hpp` `LimiterHelpers::venkatFunction`]。
- **極値節点での振る舞い** (forge の観測に直結): $\Delta_+=0$ なら $\psi=\epsilon^2/(\epsilon^2+2\Delta_-^2)$。
  したがって **$|\Delta_-|\gg\epsilon$ なら ψ≈0、$|\Delta_-|\ll\epsilon$ なら ψ≈1**。ε はまさに「節点が極値だが変動が小さい」
  場合に ψ を 1 側へ寄せるための項で、ε が増分より何桁も小さいと、この保護は働かない。
- **解く問題**: (i) Barth–Jespersen の $\min(1,\cdot)$ の微分不可能性による収束停滞、(ii) 滑らか域での精度低下。
- **収束の証拠**: 1993 論文アブストラクト「多次元では単調解を与えるが、**数桁の残差低下で収束が停滞**する。
  修正したリミッタでは収束した定常解が得られる」[一次・文書: AIAA 93-0880 アブストラクト (OpenAlex 経由)]。
  1995 JCP 版のアブストラクト冒頭も同趣旨 [二次: 検索結果の要約]。
- **K の値**: 実用コードの既定はまちまち。FUN3D 系の論文は $K=1$、$\Delta x=(6V/\pi)^{1/3}$ (median dual 体積 $V$) を
  既定として使用 [一次・文書: Ahmad ら 2023 式 (10)]、Wind-US は `KVENKAT` 既定 20 (`KVenk*Dcell`、Dcell は同体積球の直径)
  [一次・文書: Wind-US `TVD` キーワード]、Nishikawa–White–O'Connell は $K=5$ で機械零まで収束 [一次・文書]。
  Bolsoni–Azevedo (2026) は推奨範囲を $\epsilon_V\in[0.01,10]$ と書く (PDF 抽出で小数点が落ちており範囲の読みは要確認)
  [一次・文書 (抽出不完全)]。
- **ε の次元**: 原形は次元が合わない (長さ³ と 変数²)。FUN3D・原論文は無次元化した変数で使う前提。
  Nishikawa–White–O'Connell は Luke (AIAA 2007-3956) にならい、ε に**局所値** ($\rho_j$、速度は $\max(|u_j|,a_j)$) を掛けて
  次元を合わせている [一次・文書: NTRS 20230004018 式 (4.3) 付近]。Luke の原典は未読 [二次]。

### 2.2 Wang の修正 (IJNMF 33, 2000) と範囲ベースの ε

- **式**: $\epsilon^2=[\epsilon_W(q^{gl}_{\max}-q^{gl}_{\min})]^2$、$\epsilon_W\in[0.01,0.20]$ 推奨、$q^{gl}$ は**全域**の最大・最小
  [一次・文書: Bolsoni–Azevedo arXiv:2601.16291 式 (10)]。Wang 2000 のアブストラクトは「Venkatakrishnan の微分可能リミッタを
  修正した形で用いる」とだけ書き、式は本文 [一次・文書: アブストラクト / 本文は未読]。
- **SU2 `VENKATAKRISHNAN_WANG`**: 変数ごとに全域 (MPI 全体) の min/max を毎回の limiter 計算で集計し、
  `eps2 = max((K·range)², machine_eps)`、K は `VENKAT_LIMITER_COEFF` [一次・ソース: `CLimiterDetails.hpp`]。
  SU2 文書: 「メッシュ尺度に直接依存しないので無次元化なしで使える」[一次・文書: SU2 docs *Slope Limiters and Shock Resolution*]。
- **SU2 `VENKATAKRISHNAN` (既定) の ε**: `eps2 = max((RefElemLength·K)³, machine_eps)`、`RefElemLength = 1.0`
  (US 単位系なら 1/0.3048) で**全域一定**。文書は「$\bar\Delta$ は平均格子寸法」と書くが、ソース上は定数 1 [一次・ソース:
  `Common/src/CConfig.cpp:5021`、`CLimiterDetails.hpp`]。既定 K=0.05 で $\epsilon^2=1.25\times10^{-4}$、$\epsilon\approx0.011$。
  - **SU2 の既定 `REF_DIMENSIONALIZATION` は `DIMENSIONAL`** で、そのとき `Pressure_Ref = Density_Ref = Temperature_Ref = 1`
    (SI のまま解く) [一次・ソース: `CConfig.cpp:1511`、`CEulerSolver.cpp:1011`]。SU2 文書も「無次元モード外で使うと、値の大きい
    場 (圧力・温度) ほど強く制限される」と注意している [一次・文書]。
  - **forge の `methods/limiter.md` は「SU2 は解を無次元化して解いているので δ が O(1)」と書いているが、既定設定については
    当たらない** (無次元化はオプション)。forge の「`venkatK` 0.05 は SU2 既定と同値」も、**K の数値が同じだけで ε の式は別物**
    (SU2 は全域一定、forge は $h_i/L_\mathrm{ref}$ の局所比)。本ノートでは指摘にとどめ、`methods/` は編集していない。
- **OpenFOAM の範囲型緩和**: `cellLimited<Venkatakrishnan> Gauss linear k` の Venkatakrishnan 関数は ε を持たない比の形
  $\phi(r)=(r^2+2r)/(r^2+r+2)$ ($r=\Delta_+/\Delta_-$)。代わりに k<1 のとき近傍の上下限を
  $\Delta=(1/k-1)(q_\max-q_\min)$ だけ広げる (局所レンジに比例) [一次・ソース: OpenFOAM-dev `cellLimitedGrad.C`、
  `VenkatakrishnanGradientLimiter.H`]。同ファイルは「Venkatakrishnan 関数は 1 を超えうるので 1 で切っており微分可能でない、
  `cubic` を推奨」と注記 [一次・ソース]。**局所レンジが 0 (完全に一様) なら広がらない**点に注意。
- **収束の証拠**: Bolsoni–Azevedo (2026, 遷音速 NACA0012, RANS-SA, 陰解法 GMRES) は V・W・R3 の 3 つとも、強い衝撃のケース 3 で
  **6 桁にも届かず停滞**し、リミッタ 0 (1 次) なら機械零まで、2 次でも陽解法並みの小 CFL なら機械零まで行くと報告。
  原因は「近似ヤコビアンが制限付きスキームの非線形性を捉えない可能性」と仮説にとどめている。結論は「制御定数を推奨範囲に保てば
  3 者の解の差は工学的許容内」「Wang 版は格子尺度の差が大きい場合に理論上より頑健」[一次・文書: arXiv:2601.16291 §4.1・§5]。

### 2.3 Michalak & Ollivier-Gooch (AIAA 2008-776 / JCP 228, 2009)

- **内容** (2008 アブストラクト): 2 次用リミッタの高次への拡張、「精度と**効率的な収束**の要件」を論じ新しい制限手順を提案。
  遷音速・超音速の 4 次解で良好な収束、「新リミッタの一部は 2 次スキームの散逸低減にも、収束性をほとんど犠牲にせず使える」
  [一次・文書: AIAA 2008-776 アブストラクト (OpenAlex)]。
- **関数形**: OpenFOAM の `cubic` リミッタが MOG 2008 を出典に実装している形は
  $\phi(r)=r+b r^2+a r^3$ ($r<r_t$)、$\phi=1$ ($r\ge r_t$)、$a=1/r_t^2-2/r_t^3$、$b=-\tfrac32 a r_t-\tfrac{1}{2r_t}$、$r_t\in(1,2)$。
  値と傾きの両方で 1 に滑らかにつながり 1 を超えない [一次・ソース: OpenFOAM-dev `cubicGradientLimiter.H`]。
- **滑らか域の判定 (σ スイッチ)**: MOG 2009 は、近傍レンジ $(q_\max-q_\min)^2$ を $(K\Delta x)^3$ と比べて滑らかに
  リミッタを切るスイッチを導入した、と記憶しているが**本文を読めていない** (UBC リポジトリ・ResearchGate がアクセス拒否) [未確認]。
  Nejat & Ollivier-Gooch (JCP 2008) の tanh によるリミッタの局所化を MOG が改良した、という記述は二次資料にある
  [二次: arXiv:1710.07187 の緒言]。また「Nejat–Ollivier-Gooch は K=10 で凍結不要だった」という要約も検索結果にあるが原典未確認 [二次]。

### 2.4 Nishikawa の R_p リミッタ (AIAA 2022-1374) と Barth–Jespersen の微分不可能性

- **主張** (アブストラクト): (1) 滑らか域で最大 5 次精度を保つ (Venkatakrishnan は 2 次まで、MOG は 4 次まで)、(2) 1 次元でも
  完全に微分可能、(3) Venkatakrishnan の ε 修正 (ほぼ一様な領域でリミッタを滑らかに切る) と両立、(4) Venkatakrishnan より低散逸、
  (5) 無制限の再構成がすでに有界ならそのまま保つ (Venkatakrishnan は保たない)。node 中心エッジベースで検証
  [一次・文書: AIAA 2022-1374 アブストラクト (OpenAlex)]。
- **式** (SU2 実装、$D_p=|\delta|$、$D_m=|\mathrm{proj}|$):
  $D_p>2D_m$ なら 1、それ以外は $\psi=\dfrac{D_p^p+\epsilon_p+D_pS_p}{D_p^p+\epsilon_p+D_m(\delta^{p-1}+S_p)}$、
  $S_3=4D_m^2$ ほか。ε は $\epsilon_p=(K\cdot\mathrm{RefElemLength})^{p+1}$ [一次・ソース: `LimiterHelpers::r3Function` 等]。
  FUN3D は `flux_limiter = 'nishikawa'` (v13.1 以降) [一次・文書: FUN3D リリースノート]。
- **収束の証拠**: Ahmad, Park, Nishikawa, Wang, Elmiligui (NASA, AIAA 2023) — 2D 超音速バンプ (Roe、凍結なし) で R5 は 1e-12 まで収束、
  Venkatakrishnan も妥当、minmod / van Leer / van Albada は約 3000 反復で停滞。R5 は 1 で上から抑えられるので
  $\min(1,\phi)$ が不要になり、その分の微分不可能点が消える。3D の実機形状 (SA-neg) では「凍結 1500 反復」で回している
  [一次・文書: NTRS 20230004071]。
- **Barth–Jespersen**: 微分不可能 ($\min$) で収束が劣る、が通説 [二次: arXiv:1710.07187 緒言、Fluent 理論ガイド
  「非微分可能リミッタは数桁で見かけの残差が停滞する」 一次・文書]。
- **Nishikawa ら自身の課題認識**: 「ε を無くして滑らかな極値と一様域を別の方法で検出したい」「より深い反復収束のために
  リミッタ凍結を頑健に発火させる方法を探している」「K を問題ごとに事前に決めるのは難しく、元の ε は次元が合わない」
  [一次・文書: Nishikawa–White–O'Connell, NTRS 20230004018 §4 末・§5]。

### 2.5 実用コードのリミッタ凍結

| コード | オプション | 意味 | 注意点 | 確認度 |
| --- | --- | --- | --- | --- |
| SU2 | `LIMITER_ITER` (既定 999999) | `InnerIter <= LimiterIter` の間だけ `SetPrimitive_Limiter` を呼ぶ。以後は最後に計算した ψ を流れの再構成で使い続ける | **乱流・化学種 (`CScalarSolver`) は凍結ではなく、`LIMITER_ITER` 以降リミッタを掛けない (ψ=1 相当)** — 本ソース版で確認。これを「凍結に揃える」PR #2910 は 2026-09-17 にマージされず閉じられた (保守者「PR は説明どおりのことをしていない」)。離散随伴用に別途 `FROZEN_LIMITER_DISC` | 一次・ソース (`CEulerSolver.cpp:1684,1826`、`CScalarSolver.inl:100,142`)、一次・文書 (PR) |
| FUN3D | コマンドライン `--freeze_limiter xx` (旧)、名前リスト `&inviscid_flux_method` の `freeze_limiter_iteration(s)` (新) | xx 時間ステップ後に全域の ψ を凍結。「リミッタ使用時に **stall または ring** する収束の改善に有用」。2014 マニュアルは barth / venkat / h 系での使用を挙げる。「stencil 型は freezable、edge 型 (minmod 等) は違う」という区別は FUN3D 講習資料の検索要約でのみ確認 [二次] | 凍結中も再構成は毎ステップ行い、**面への外挿が失敗した点だけ ψ を再計算**する。restart 後も凍結を続けるには `--freeze_limiter 0` が必要。v13.6 で非定常では副反復番号に紐づけ。v14.2 で HLLE++ のショックスイッチも凍結対象に。名前リストの変数名は資料により `freeze_limiter_iteration` / `freeze_limiter_iterations` の両表記があり、現行マニュアルで未確定 | 一次・文書 (2014 マニュアル第 6 章、リリースノート) |
| Wind-US | `DEBUG 57 mode` | 非構造格子の slope limiter を mode 反復後に凍結 | 記載は 1 行のみ | 一次・文書 (`debug.html`) |
| OpenFOAM / Fluent | — | 凍結オプションは見つからなかった。Fluent は「limiter filter」(Standard limiter) があるが中身は非公開 (Ansys 掲示板で説明拒否) | — | 一次・文書 (未発見であることの確認) |
| CFL3D (構造) | `iflim` 3/4 | 凍結ではなく、滑らかリミッタのカットオフを格子寸法 (iflim=4 は総セル数基準の等方カットオフ) で決める | — | 一次・文書 (v6 新機能) |
| USM3D | — | HANIS (非線形制御・GCR) で収束を改善した論文に、リミッタ凍結の記述は見つからなかった | — | 一次・文書 (Pandya ら 2015 を grep) |

**実計算での使われ方** (NASA): 低ブーム機 (4.7M 節点、SA-neg) で「収束を速めるため 1500 反復で凍結」、X-59 模型で
「HANIM 600 反復・旧ソルバ 2000 反復で凍結」[一次・文書: NTRS 20230004071]。Nishikawa–White–O'Connell は 2 次・3 次とも
「K=5、400 反復で凍結、10 桁以上・機械零まで収束」(同じ論文の別ケースは凍結なし) [一次・文書]。
第 1 回 Sonic Boom Prediction Workshop の要約では「力は収束しても圧力シグネチャは変化し、RMS 残差が停滞しているときは
抽出位置の残差を監視する必要があった」と、リミッタによる停滞が誤差源として挙げられている [一次・文書 (Ahmad ら 2023 の引用経由)]。

### 2.6 境界での扱い

- **SU2**: limiter の極値探索は `geometry.nodes->GetPoints(iPoint)` の全エッジ近傍 (境界節点も同じ)、境界・角の特別扱いは無い
  (周期境界だけ通信で min/max を合わせる)。境界条件の流束 (`BC_Inlet` ほか) は節点値で作り、**再構成も ψ も使わない**
  [一次・ソース: `computeLimiters_impl.hpp`、`CEulerSolver.cpp` の `BC_*` に `MUSCL_Reconstruction` の呼び出しが無い]。
  → forge の「入口半割面の流束は ψ を使わず境界状態から作る」「近傍は内部双対面を共有する実ノードのみ」と同型。
  再構成が非物理 (負の P・ρ・音速) になったエッジだけ 1 次に落とす `bad_recon` はある (起動時用と注記)。
- **OpenFOAM**: `cellLimitedGrad` は**境界面の値 (非結合パッチも) を極値の集合に含める** [一次・ソース: `cellLimitedGrad.C`]。
  入口の指定値や壁の値が近傍の上下限を広げる形になる。
- **FUN3D**: `--limit_near_walls_less` (h 系リミッタ用) は壁に近づくとリミッタを「切る」(制限を弱める)。混合要素格子で
  壁熱流束・摩擦が改善、四面体格子では悪化、頑健性が落ちることがある、と注意書き [一次・文書: 2014 マニュアル]。
  境界節点を 1 次にする・角を特別扱いする記述は見つからなかった。
- **Wind-US**: 構造格子の `BOUNDARY TVD` で境界 (結合境界) にだけ別の圧縮係数 (0 で 1 次) を与えられる。非構造は
  `DEBUG 64` で極値探索を face 近傍から node 近傍 (頂点共有セル) に広げられ、「線形関数を任意格子で制限しない」が TVD は
  保証しない [一次・文書]。
- **SU2 `SHARP_EDGES` / `WALL_DISTANCE`**: 鋭い縁・壁からの距離 $d$ で ψ に raised-sine の幾何係数
  ($\tfrac12(1+s+\sin(\pi s)/\pi)$、$s=d/(c\,\epsilon_1)-1$) を掛けて**縁・壁の近くで制限を強める**。設定テンプレート上は随伴流
  (`SLOPE_LIMITER_ADJFLOW`) 用 [一次・ソース + `config_template.cfg`]。
- **入口∩壁の角 (一様入口速度と 0 速度の壁節点が隣接する特異角) を名指しで扱った文献・コードは見つからなかった** [未発見]。

### 2.7 滑らか域でリミッタを切るセンサ (定常 RANS)

- **FUN3D h 系** (`hminmod`/`hvanleer`/`hvanalbada`/`hsmooth`/`hvenkat`): 選んだリミッタに**圧力ベースの経験的リミッタ**を自動で
  掛け合わせる。v14.2 で単独の `hpressure` も追加 [一次・文書: マニュアル・リリースノート]。中身は Gnoffo (AIAA 2010-1271) の
  圧力リミッタで、近傍の最大/最小圧力比で 0〜1 に滑らかに切り替え、全変数の ψ に掛ける。「比較的大きな K を使っても衝撃で
  単調性を保てるが、接触不連続は検出しない」[一次・文書: Nishikawa–White–O'Connell 式 (4.6) と §5]。
- **Kitamura & Shima (AIAA J 50(6), 2012) の second limiter**: tanh 型の関数で亜音速・淀み域のリミッタを切る、パラメタ無し、
  と二次資料は述べる [二次: arXiv:1710.07187]。式・変数 (Mach か圧力か) は原典未読 [未確認]。
- **Nejat & Ollivier-Gooch (JCP 227, 2008)**: tanh でリミッタの活性域を限定 [二次]。
- **Ducros 型**: forge は Roe 経路で既にある (`methods/limiter.md`)。SLAU 経路には掛からない。定常 RANS でリミッタのスイッチに
  使った一次資料は今回見つけていない。
- **CENO の smoothness indicator** (Ivan & Groth, JCP 257, 2014) は高次再構成の切替に使う例 [一次・文書 (存在のみ)]。

### 2.8 浮動小数点・丸めとリミッタの chatter

- **Berger, Aftosmis, Murman (AIAA 2005-490)**: 「非滑らかなリミッタは**ほぼ一様な領域でさえ chatter** を起こし、時間進行でも
  Newton でも定常に達しにくくする」(緒言)。数値例 (ONERA M6, M 0.5) で、境界面 (2:1 界面) のセルに scalar BJ を使うと
  **密度残差が約 3 桁で停滞**、recentering (制限後の勾配を単調性多角形の角に固定する) にすると収束。cut cell の scalar min でも
  同じ停滞が出て、「単調で線形保存なのに chatter するのは scalar min の実装が原因。**勾配ベクトルが反復ごとに過度に回転する**
  のが仮説」と書く [一次・文書: NTRS 20050182793 §V]。
  forge の「ψ を決める面が 0.85〜1.0/step で入れ替わる」観測と同じ型の現象が報告されている唯一の一次資料。
- **Nishikawa (2022)**: 検索結果の要約に「非構造リミッタは機械零レベルの雑音に敏感で、ほぼ一様な領域 (自由流など) を含む問題で
  非線形ソルバを停滞させる」という趣旨の記述が出るが、原文を確認できていない [二次]。
- **単精度 (float32) とリミッタの chatter を直接扱った一次資料は見つからなかった** [未発見]。SU2 は ε の下限を倍精度の
  machine epsilon (`numeric_limits<passivedouble>::epsilon()`) に置いている (= 「ε=0 の割り算回避」だけで、丸め雑音の吸収は
  意図していない) [一次・ソース]。

---

## 3. forge の観測との対応

forge の観測を 3 つに分け、どの手法がどこに作用するかを対応させる (「効く」とは書かない — 機構が噛み合う候補を挙げる)。

### 3.1 「節点自身が近傍の極値 ($\Delta_+\approx0$)」

$\Delta_+=0$ での $\psi=\epsilon^2/(\epsilon^2+2\Delta_-^2)$ により、ψ は **$|\Delta_-|/\epsilon$ の大小だけで 0 と 1 に分かれる**。
forge の入口境界列では、面の増分 $\hat\Delta_-$ (plan #4r: ψδ_m の変動 5e-4 程度、勾配の射影は $10^{-4}$〜$10^{-3}$) が
$\hat\epsilon\approx5\times10^{-8}$ より 3〜4 桁大きいので、極値になった瞬間に ψ≈0、極値でなくなると ψ≈1 になる。
極値か否かは近傍値との ulp 単位の大小で決まる (plan #4r) ので、**ψ の切替は ε が丸め幅程度しかないことの帰結として説明できる**
(観測との整合であって、原因の確定ではない。plan #5er の判定は保留)。

この型に作用する手法:

- **ε を大きくする / 局所 h に依存させない** (2.2 の SU2 定数 ε、Wang の範囲 ε): 増分が ε より小さい限り ψ→1 になるので、
  「極値か否か」の ulp 判定が ψ に伝わらなくなる。目安として、SU2 を無次元化して使った場合の $\epsilon\approx0.011$ (q_ref 比) は、
  forge の入口列の $\hat\epsilon\approx5\times10^{-8}$ の約 $2\times10^5$ 倍で、観測された増分 ($10^{-4}$〜$10^{-3}$) より上にある。
  forge の $\hat\epsilon=(Kh_i/L_\mathrm{ref})^{3/2}$ は近壁の細かい節点ほど小さくなる (近壁 $h\sim40\,\mu$m で $5\times10^{-8}$) 。
- **MOG 型の滑らか域スイッチ**: 近傍レンジが小さい所でリミッタを切る [未確認の詳細]。
- **OpenFOAM の範囲緩和 (k)**: 上下限を局所レンジに比例して広げると $\Delta_+>0$ になる。ただし入口列のように近傍がほぼ同値
  (一様入口) ならレンジ自体が小さく、効き方は不明。
- **極値探索の近傍を広げる** (Wind-US `DEBUG 64`): 極値になりにくくなる。

### 3.2 「$\hat\epsilon\approx$ ulp」(丸め)

ε が float32 の丸め幅と同程度という状況は、ε を「割り算回避」以上の意味で機能させていない状態。文献で丸めと chatter を直接
扱ったのは Berger らの「min の切替による勾配の回転」(2.8) で、対処は recentering (勾配を単調性多角形の角に固定) だった。
float32 そのものを論じた資料は無い。作用する候補は 3.1 と同じ (ε の尺度の変更) に加え、

- **ψ の凍結** (2.5): 切替そのものを止める。FUN3D は凍結後も「外挿が失敗した点だけ再計算」する安全弁を持つ。
- **R_p 系** (2.4): $\min(1,\cdot)$ の折れ目は消えるが、面ごとの min と極値判定は残るので、極値節点の切替は残る見込み [推測]。

### 3.3 「境界の角 (一様入口 × no-slip 壁)」

入口境界の節点は一方の近傍しか持たず、壁節点 (速度 0) が近傍の極値に入る。文献・コードで**この角を名指しで扱った例は見つからなかった**。
近い扱いは:

- **境界値を極値集合に含める** (OpenFOAM): 入口の指定状態が上下限に入る。forge は境界半割面・ghost を除外している (plan §4.2)。
- **境界近傍だけ別の制限** (Wind-US 構造格子 `BOUNDARY TVD`、FUN3D `--limit_near_walls_less`、SU2 `WALL_DISTANCE` の幾何係数):
  向きは逆のもの (制限を強める/弱める) が混在しており、「角の節点を 1 次にする」「角の ψ を固定する」を一般的な処方として
  書いた資料は無い。
- SU2 は境界を特別扱いせず、境界流束も節点値 (1 次) で作る — forge と同型なので、**同じメッシュ・同じ BC で SU2 にも同じ角の
  ψ 振動が出るか**は、forge 固有か否かの切り分けに使える (未実施)。

---

## 4. 候補の比較

採否は決めていない。複雑さは forge への実装量の見立て、効果は §3 の機構からの見込みで、いずれも未検証。

| 候補 | 実装の複雑さ | 入口列の床への期待 (機構) | 精度・衝撃捕獲へのリスク | 文献上の根拠 |
| --- | --- | --- | --- | --- |
| A. ε を局所 h に依存させない定数にする (SU2 型、$\hat\epsilon^2=(K\cdot c)^3$) | 小 (式 1 行 + キー) | 3.1・3.2 に直接作用。増分 < ε の節点で ψ→1 | ε 以下の振動 (q_ref の数 %) を許す。全域で制限が弱まる。既定値の再検証 (Sod 厳密解・近傍逸脱、`limiter-config-simplify` §4) が要る | SU2 既定、Venkatakrishnan 原論文 (無次元変数前提) |
| B. Wang の範囲 ε (変数ごとの全域レンジ) | 小〜中 (全域 reduction を毎回、または凍結) | A と同じ機構。レンジが反復で変わる分、収束までは ε も動く | A と同様。全域レンジが局所の弱い衝撃に比べて大きいと弱い衝撃を見逃す [推測] | SU2 `VENKATAKRISHNAN_WANG`、Bolsoni–Azevedo 2026 |
| C. 滑らか域スイッチ (MOG σ 型、または Kitamura–Shima の低 Mach・淀み域で切る型) | 中 | 入口列は M 0.11 の亜音速・ほぼ一様なので、スイッチで ψ=1 に固定される可能性 | スイッチ自体の閾値と微分可能性。接触不連続の扱い | MOG 2009 [未確認の詳細]、Kitamura–Shima 2012 [二次]、FUN3D h 系 (圧力比) |
| D. ψ の全域凍結 (N 反復後) | 小 (forge には診断用の `FORGE_DIAG_PSI_DUALEVAL` の退避・差し替えの仕組みがある) | 3.2 の切替を止める。plan #5er は D(res_ro)=0.91 で「入口 ψ の差し替えが入口残差に直接大きく寄与する」は保留 | 凍結時点の ψ に依存した別の固定点になる。衝撃が凍結後に動くと過小制限。restart で凍結状態を保存しないと再計算される (FUN3D は明示フラグ)。`limiter-config-simplify` §4.3 の決定の見直しが必要 | SU2 `LIMITER_ITER`、FUN3D `--freeze_limiter`、Wind-US `DEBUG 57`、NASA の実計算で常用 |
| D'. 局所凍結 (切替を繰り返す節点だけ凍結、または FUN3D 同様に再構成失敗点だけ再計算) | 中 | D より介入範囲が狭い | 判定条件の設計が必要。文献に「切替節点だけ凍結」の例は見つからなかった | FUN3D の「失敗点だけ再計算」が部分的に近い |
| E. 境界値 (入口指定状態・壁値) を極値集合に含める | 中 (node の境界半割面の状態を limiter に渡す) | 3.3 に作用。ただし一様入口では入口状態も近傍とほぼ同値で、極値判定は変わらない可能性 | 小さい見込み (境界値は物理的な上下限) | OpenFOAM `cellLimitedGrad` |
| F. 入口列 (境界節点) の内部エッジを 1 次にする / ψ を固定する | 小 | 床の 70 % を占める列の切替を消す | 入口近傍の精度が 1 次。入口境界層の発達に影響しうる。`badReconFallback`・`FORGE_CONTACT_1ST` は代用にならない (plan §4.2) | 構造格子の Wind-US `BOUNDARY TVD` のみ (非構造の一般処方は未発見) |
| G. R_p リミッタ (Nishikawa) に替える | 中 (関数差し替え + ε の次数) | $\min(1,\cdot)$ の折れ目は消える。極値節点の切替そのものは ε 次第で残る見込み | 低散逸側。forge の Sod・近傍逸脱ゲートでの再検証が要る | Nishikawa 2022、Ahmad ら 2023、SU2/FUN3D 実装 |
| H. 極値探索の近傍を広げる / 範囲緩和 (OpenFOAM k) | 中 | 極値になりにくくする | TVD 性が弱まる (Wind-US 注記) | Wind-US `DEBUG 64`、OpenFOAM k |
| I. 陰解法でリミッタを線形化する (Newton 型) | 大 | 「近似ヤコビアンが非線形性を捉えない」仮説への対処 | 実装・コスト大。ψ が微分不可能な点では効かない | Bolsoni–Azevedo の仮説のみ (未検証)、HANIM/HANIS は非線形制御側 |

補足:

- A・B は「ε の尺度」の話で、`limiter-config-simplify` が既定化の根拠にした「K=0.05 は SU2 既定と同値」という前提 (2.2 で指摘した
  式の違い) にも触れる。A/B を試す場合は、既定化のゲートだった厳密解比較と近傍逸脱の再測定がセットになる。
- D は plan #5er の結果 (枝分かれ二重評価で D(res_ro)=0.91, Q=0.70、判定保留) と、`limiter-config-simplify` §4.3 の
  「入れない」の決定の両方に関わる。文献上は最も普及した実用的対処 (SU2・FUN3D・Wind-US が全部持つ) だが、凍結後の解が
  凍結時点に依存する点はどのコードの文書も定量的に論じていない。

---

## 5. 未確認事項

1. **MOG 2009 の滑らか域スイッチ σ の式と閾値**、および ε の扱い (本文未読。UBC・ResearchGate がアクセス拒否)。
2. **Venkatakrishnan 1995 本文の推奨 K と、凍結への言及の有無** (アブストラクトのみ確認)。
3. **Wang 2000 本文の ε の式** (Bolsoni–Azevedo の引用式でのみ確認)。
4. **Nishikawa 2022 本文** (R_p の ε の推奨、機械零雑音への言及の原文、node 中心エッジベースでの境界の扱い)。
5. **Kitamura & Shima 2012 の second limiter の式** (Mach か圧力か、tanh の引数)。
6. **FUN3D の現行マニュアルでの名前リスト変数名** (`freeze_limiter_iteration` か `freeze_limiter_iterations` か) と既定値。
   2014 版マニュアルはコマンドライン `--freeze_limiter` のみ記載。現行の chapter-6 HTML は 404。
7. **Luke (AIAA 2007-3956) の次元整合な ε の式** (Nishikawa–White–O'Connell の引用経由のみ)。
8. **SU2 で restart 時に凍結 ψ が保存されるか** (ソース未確認)。
9. **Fluent の "limiter filter" の中身** (非公開)。
10. **float32 とリミッタ chatter を直接論じた文献** — 見つからなかった。
11. **入口∩壁の角を名指しで扱った文献・コード** — 見つからなかった。
12. forge 側: 同一メッシュ・同一 BC の SU2 (無次元化あり/なし) で入口列に同じ ψ の切替が出るか (未実施。
    [`procedures/su2-cross-check.md`](../../procedures/su2-cross-check.md) の手順で切り分けに使える)。

## 出典

- Venkatakrishnan, V., "On the accuracy of limiters and convergence to steady state solutions," AIAA 93-0880, 1993. doi:10.2514/6.1993-880 (アブストラクトは OpenAlex)。[NTRS 19930040944](https://ntrs.nasa.gov/citations/19930040944)
- Venkatakrishnan, V., "Convergence to steady state solutions of the Euler equations on unstructured grids with limiters," JCP 118, 120–130, 1995. doi:10.1006/jcph.1995.1084
- Wang, Z.J., "A fast nested multi-grid viscous flow solver for adaptive Cartesian/Quad grids," IJNMF 33, 657–680, 2000.
- Michalak, C., Ollivier-Gooch, C., "Limiters for unstructured higher-order accurate solutions of the Euler equations," AIAA 2008-776. / "Accuracy preserving limiter for the high-order accurate solution of the Euler equations," JCP 228, 8693–8711, 2009.
- Nishikawa, H., "New unstructured-grid limiter functions," AIAA 2022-1374. doi:10.2514/6.2022-1374
- Ahmad, N., Park, M.A., Nishikawa, H., Wang, L., Elmiligui, A., "Evaluation of Limiter Functions for Supersonic Applications," AIAA 2023. [NTRS 20230004071](https://ntrs.nasa.gov/api/citations/20230004071/downloads/aiaa-2023-cst-Ver9.pdf)
- Nishikawa, H., White, J.A., O'Connell, M.D., "Toward a Third-Order Accurate, Second-Derivative-Free, Shock-Capturing Finite-Volume Method for Hypersonic Flows on Tetrahedral Grids." [NTRS 20230004018](https://ntrs.nasa.gov/api/citations/20230004018/downloads/NishikawaWhiteOConnell_v7.pdf)
- Oliveira, F.B., Azevedo, J.L.F., "A Study of Improved Limiter Formulations for Second-Order Finite Volume Schemes Applied to Unstructured Grids," [arXiv:2601.16291](https://arxiv.org/abs/2601.16291), 2026.
- Berger, M., Aftosmis, M.J., Murman, S.M., "Analysis of Slope Limiters on Irregular Grids," AIAA 2005-490. [NTRS 20050182793](https://ntrs.nasa.gov/api/citations/20050182793/downloads/20050182793.pdf)
- Gnoffo, P.A., "Updates to Multi-Dimensional Flux Reconstruction for Hypersonic Simulations on Tetrahedral Grids," AIAA 2010-1271.
- Kitamura, K., Shima, E., "Simple and Parameter-Free Second Slope Limiter for Unstructured Grid Aerodynamic Simulations," AIAA J 50(6), 2012. doi:10.2514/1.J051269
- SU2 ソース (`.external/su2-src`, commit 12eb826f): `SU2_CFD/include/limiters/CLimiterDetails.hpp`, `computeLimiters_impl.hpp`, `SU2_CFD/src/solvers/CEulerSolver.cpp`, `SU2_CFD/include/solvers/CScalarSolver.inl`, `Common/src/CConfig.cpp`, `config_template.cfg`。SU2 文書 [Slope Limiters and Shock Resolution](https://su2code.github.io/docs_v7/Slope-Limiters-and-Shock-Resolution/)。[SU2 PR #2910](https://github.com/su2code/SU2/pull/2910)
- FUN3D Manual Chapter 6 (2014-06-05 版 PDF) [fun3d.larc.nasa.gov](https://fun3d.larc.nasa.gov/papers/Website_June2014_Chapter06_Analysis.pdf)、[FUN3D Chapter 1 (リリースノート)](https://fun3d.larc.nasa.gov/chapter-1.html)
- OpenFOAM-dev ソース: `src/finiteVolume/finiteVolume/gradSchemes/limitedGradSchemes/cellLimitedGrad/` (`cellLimitedGrad.C/.H`, `gradientLimiters/VenkatakrishnanGradientLimiter.H`, `cubicGradientLimiter.H`)
- Wind-US User's Guide: [`TVD`](https://www.grc.nasa.gov/www/winddocs/user/keywords/tvd.html), [`BOUNDARY TVD`](https://www.grc.nasa.gov/www/winddocs/user/keywords/bndrytvd.html), [`DEBUG`](https://www.grc.nasa.gov/www/winddocs/user/keywords/debug.html) (57, 64)
- Ansys Fluent Theory Guide, Gradient Limiters ([v24.2](https://ansyshelp.ansys.com/public////Views/Secured/corp/v242/en/flu_th/flu_th_sec_grad_limit.html))、[CFL3D v6 New Features](https://nasa.github.io/CFL3D/Cfl3dv6/cfl3dv6_new.html)
- Pandya, M.J., Diskin, B., Thomas, J.L., Frink, N.T., "Improved Convergence and Robustness of USM3D Solutions on Mixed Element Grids," [NTRS 20150005707](https://ntrs.nasa.gov/citations/20150005707)
