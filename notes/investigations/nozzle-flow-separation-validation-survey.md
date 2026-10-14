# ノズル内剥離 (過膨張) の検証データ調査 — forge 剥離評価機能の第一検証ケース選定 (2026-10-11)

- **目的**: 外気圧を与えて過膨張剥離を評価する機能 (ノズル外の外部領域・プルーム領域つき計算) を forge に作るにあたり、剥離位置・壁圧分布・遷移 NPR を実験と比べる**検証データセット**を選ぶ。
- **対象機種**: ③ ベル (真空用と海面作動の両方)・④ デュアルベル・⑤ SERN。
- **forge 側の前提**: 先に検証するのは冷気 (空気 / N2) の軸対称または 2D 平面の定常 RANS (k-ω SST)。非定常 (URANS・横荷重) は後回し。形状を座標で完全に再現できることを必須条件とする。
- **位置づけ**: 調査記録であり、検証ケースの採否は起票する plan で決める (**plan 未反映**。対応 plan はまだ無い)。親計画は [`plans/active/tooling-nozzle-design-tool.md`](../../plans/active/tooling-nozzle-design-tool.md)。既存サーベイ [nozzle-optimization-tool-survey.md](nozzle-optimization-tool-survey.md) の A9 と付録 10 (デュアルベル未調査) の続き。
- **保存した PDF**: `/home/sano/work/forge/papers/nozzle_separation/` (git 管理外)。一覧は §8。

**凡例 (引用の検証状態)**

| 記号 | 意味 |
| --- | --- |
| ✔ | 一次ソース (本文 PDF、NTRS、出版社 DOI) で書誌と中身を確認した |
| △ | 書誌は DOI (Crossref) か NTRS で確認したが、中身は抄録か二次的な言及だけで確認した |
| ✖ | 確認できなかった |

**輪郭の公開度**: 「全点」= 座標表がある / 「寸法」= 直線・円弧の寸法表で形状が閉じる / 「MOC 再構成」= 設計条件 (設計 Mach・面積比・切断位置など) が書かれていて特性曲線法で作り直せる / 「図のみ」= 図から読み取るしかない / 「非公開」。

---

## 0. 結論 — 第一検証ケースに推す 5 件

| 順 | データセット | 次元・型・作動ガス | 推す理由 | 入手性 |
| --- | --- | --- | --- | --- |
| 1 | **DLR TIC — FSCD/ATAC workshop test case 1A** [S07] | 軸対称・TIC ($M_D$=5.15, $\varepsilon\approx14.0$)・N2 冷気 ($p_0$=2.5 MPa, $T_0$=283 K)・NPR 25.25 | 輪郭 ($X/R^*$ = −4.5〜9.0)、壁圧 $p_w/p_0$、マッハディスク・斜め衝撃・滑り線の座標が**すべて表で**公開されている。同じ条件に 19 本の RANS 解 (SST・SA・k-ω・Baldwin-Lomax) と剥離位置の誤差表が付いている。SST は $X_{sep}$ を −9〜−17 % (上流側) に外し、SA と k-ω は −2.4〜+4.1 % に収まった | DLR elib でオープン ✔ |
| 2 | **NASA MSFC TIC / PAR (TOP 放物線)** [R10][R09] | 軸対称・$\varepsilon$=30.5・$D_t$=38.1 mm・加熱乾燥空気 (66 °C, $p_c$=10.2 atm)・NPR 5.8〜380 | 0° と 180° の 2 方位の壁圧が付録 A に**表**で載る (NPR ごと、$p_w/p_c$ と $p_w/p_a$)。TIC は FSS のみ。PAR は増圧で FSS→RSS (NPR 23.6)、減圧で RSS→FSS (NPR 12.7) のヒステリシスがあり、cap-shock と RSS の検証に使える。TIC 輪郭は設計法 (理想ノズルを等面積比 15° 円錐の 79 % 長で切断) と $x_e/r^*$=13.812、初期膨張 22.68°、スロート下流半径 0.2 $r^*$ から MOC で再構成できる。PAR 輪郭は図のみ | NTRS でオープン ✔ |
| 3 | **NASA LaRC 2D-CD: TP-1704 の B1/B2 [MPR80] と Hunter のノズル [H98]** | 2D 平面 CD・$A_e/A_t$=1.80 / 1.797・空気 (300 K / 約 294 K) | TP-1704 は直線とスロート円弧の**寸法表**で形状が閉じ、上下フラップ 3 列と側壁 1 列の壁圧が NPR 1.96〜8.73 で**表**になっている。Hunter は座標表と半径で形状が閉じ、NPR 1.25〜9.5 の中心線・側壁寄り壁圧 (図)、シュリーレン、推力効率がある。Hunter のノズルには PAB3D・NPARC・Fluent (SST ほか) の RANS 比較が複数ある。NPR < 2 は 3 次元の部分再付着剥離で、2D RANS が外れる既知の領域 | NTRS でオープン ✔ |
| 4 | **NASA TM-86270 SERN 43 形状** [RL84] (補完: TP-1999-209138 [DA99]、NUAA [YXYM15]) | 2D 平面 SERN・空気 300 K・NPR 2〜10 | 上下フラップの座標と中心線壁圧がともに**表**。設計 NPR 3.3〜4.1 より低い側でランプ上に衝撃と剥離 (圧力平坦部) が出る。ただし面積比が小さい航空機用 (内部 1.05〜1.50、外部 1.15〜2.39) で、スクラムジェット SERN の強い過膨張ではない。強い過膨張は [DA99] (設計 NPR 102 の形態を NPR ≤ 13 で試験、3D 剥離を観測) と NUAA [YXYM15] (RSS のランプ側・フラップ側) で補う | NTRS でオープン ✔ |
| 5 | **DLR P6.2 軸対称デュアルベル** [NGS07][NGS10] (CFD: [SchG16]) | 軸対称・TIC ベース + 等圧 (CP) 延長・N2 冷気 | **座標表が公開されたデュアルベルは見つからなかった**。最も条件がそろうのがこのモデルで、$R^*$=10 mm、$\varepsilon_b$=11.3、$\varepsilon_e$=27.1、$L_b/R^*$=6.2、$L_e/R^*$=8.3、変曲角 $\alpha$=9° から MOC で再構成する。遷移 NPR ≈ 48〜51、ヒステリシス約 20 %。P6.2 の別モデル (正勾配延長) では軸対称 RANS (SA・RSM) が遷移 NPR を 1 % で再現した [SchG16] | EUCASS 2007 版は elib でオープン ✔、JPP 版は有償 |

推薦の根拠を三つ補足する。

- 1〜3 は形状が表か寸法で閉じており、壁圧も表で取れる (1 は本文テキストとして抽出できる。2・3・4 の表はスキャン画像で、転記と読み取り検査が要る)。
- 1 と 2 は同じ「TIC・FSS」で、ノズル規模 ($D_t$ 20 / 38 mm)・作動ガス (N2 / 空気)・NPR 範囲が違う。forge の SST が 1 で文献の SST と同じ方向 (上流側) に外れるかを見てから、2 で NPR 依存 ($X_{sep}$ 対 NPR) を確かめる順になる。
- ④ デュアルベルと ⑤ SERN には「形状が全点公開・壁圧が表・RANS 比較あり」の三条件をそろえたデータが無い。④ は 5 の MOC 再構成に頼り、⑤ は 4 で 2D SERN の剥離を確かめたうえで、スクラムジェット SERN 相当は [DA99] と NUAA の定性比較にとどめる。

---

## 1. 調査範囲と選定条件

選定条件は次の 5 つ。

| 条件 | 内容 |
| --- | --- |
| (a) 形状 | 輪郭を座標で再現できる (全点 > 寸法 > MOC 再構成 > 図のみ) |
| (b) 流れ | 冷気・軸対称または 2D 平面・定常で比較できる量がある (時間平均壁圧、剥離位置) |
| (c) 計測 | 壁圧分布を複数 NPR で、剥離位置、可能ならシュリーレン・遷移 NPR・ヒステリシス |
| (d) 先行 CFD | 同じデータに RANS の結果があり、乱流モデルごとのずれ方が分かる |
| (e) 入手性 | NTRS・機関リポジトリでオープン |

検索は NTRS API、Crossref、OpenAlex、Semantic Scholar、DLR elib と Web 検索で行った。本文を読んだのは §8 の PDF と、手元の `papers/nozzle_design/dual_bell/` の 7 本。

---

## 2. データセット

### 2.1 ベル・円錐ノズル (軸対称)

| 文献 | 機関・試験台・条件 | ノズル型・形状公開度 | 計測量・不確かさ | 先行 CFD | forge 適性 |
| --- | --- | --- | --- | --- | --- |
| **[S07]** Stark & Hagemann, EUCASS 2007 (FSCD/ATAC test case 1A)。実験の詳細は [SW06][SW09] ✔ | DLR Lampoldshausen P6.2 水平リグ。乾燥 N2、$p_0$=2.5×10⁶ Pa、$T_0$=283 K、海面雰囲気。test case は NPR=25.25 の 1 点。スロート Re は約 8×10⁶ (本ノートで $p_0$・$T_0$・$D_t$ から等エントロピー関係と Sutherland 則で概算) | TIC、$M_D$=5.15、$D_t$=20 mm、発散部長 90 mm、アクリル製で設計輪郭との差 < 5 µm。**全点** (表 3、収縮部から出口まで約 150 点)。出口 $R/R^*$=3.744 → $\varepsilon$≈14.0 (本ノートで算出) | 壁圧 (ポート間隔 2.5 mm、Kulite、1 kHz、160 Hz LPF) を $p_w/p_0$ で表に。剥離位置 $X_{sep}/R^*$=6.6 (最大勾配の接線と全流れ MOC 壁圧の交点で定義)。マッハディスク $X_{sh}/R^*$=9.85、斜め衝撃・反射衝撃・滑り線・マッハディスクの座標表 (シュリーレン)。不確かさの数値は記載なし | 19 解 (表 1・2)。$X_{sep}$ の誤差: 標準 SST −9.1〜−17.4 %、SST の改良版 −9.7〜−34.4 %、Östlund の SST 改良 +3.3 % (As=2.0。As=3.23 では −10 %)、SA −2.4 / −1.1 / +4.1 % (SA 改良版は −21.8 %)、k-ω Wilcox +2.7 %、k-ω Sarkar (NASA Ruf) −2.4 %、Baldwin-Lomax −13.2 %。マッハディスク位置は標準 SST で −6.8〜−8.7 %。剥離後の逆流域の最大壁圧を 1〜4 % 高く出す解ほど剥離が上流に寄る | **◎** 表で全部そろう。単一 NPR なので NPR 依存は 2 番 [R10] で補う。表 3 は $X/R^*$=4.900 の点が 5.800 の後に重複して載り、反射衝撃と滑り線の座標が同一値 (原典の転記ミスとみられる) |
| **[R10]** Ruf, McDaniels, Brown, AIAA 2010-6813 (NTRS 20100033270) ✔ / **[R09]** AIAA 2009-4856 (NTRS 20100017649) ✔ | NASA MSFC Nozzle Test Facility (真空槽、2 段エジェクタ)。加熱乾燥空気、$p_c$=10.2 atm (150 psia)、66 °C。$p_c$ 一定で槽圧を下げて NPR を上げる。スロート Re は約 5×10⁶ (本ノートで概算) | TIC と PAR (TOP)。$\varepsilon$=30.5、$D_t$=38.1 mm、アルミ製。TIC: 理想輪郭を等面積比 15° 円錐の 79 % 長で切断、$x_e/r^*$=13.812、初期膨張 22.68°、スロート下流半径 0.2 $r^*$ (3.8 mm)、収縮 30°・半径 0.75 in の円弧 → **MOC 再構成**。PAR: TDK の skewed parabola、初期角 40°、出口 $x/r^*$≈12 → **図のみ** (出口角の記載なし) | 静圧ポート各方位約 30 点 (間隔 8.6 mm、4 Hz)。付録 A 表 A1〜A10 に NPR ごとの $p_w/p_c$・$p_w/p_a$ (TIC は NPR 5.8〜380.6)。$p_{sep}/p_a$ 対 $x$ (図 16)、ノズル底面圧、横荷重モーメント (相対値)。PAR: 増圧で NPR 23.1 から不安定化し 23.6 で RSS、NPR 約 55 から部分再付着、62〜70 で全流れ。減圧で RSS→FSS は約 12.7 ([R09] は 23.8 / 13.4 と記載)。注意: スロート直後に 0.025〜0.08 mm の段差 (加工公差) があり、全流れ壁圧が前半で設計値からずれる。ディフューザ入口を出口に近づけると剥離位置が下流へずれる | 設計段階の CFD が RSS 遷移を NPR 約 24 と予測 (モデル名なし)。同じ著者は TC1A に k-ω Sarkar で参加 (−2.4 %) | TIC **◎/○** (MOC 再構成の確度を設計壁圧曲線と照合して確かめる必要)。PAR **○** (RSS・ヒステリシスの唯一の表データだが輪郭が図のみ) |
| [Re01] Reijasse et al., AIAA 2001-3556 △ | ONERA R2Ch 吹き出し風洞 | TIC と TOC (フィルム噴射可)。輪郭: 未確認 | シュリーレン (TOC の cap-shock)、中心線ピトー、壁圧と rms (TIC の干渉域で局所圧の最大 40 %、剥離域で 5〜8 %) | — | △ 形状未確認 |
| [F00] Frey et al., AIAA 2000-3777 △ / [HFK02] Hagemann, Frey, Koschel, JPP 18(3) 2002 △ | DLR P6.2 | 同性能の TOP と TIC。輪郭: 未確認 | 壁圧、カラーシュリーレン、横荷重 (TOP の FSS⇔RSS 遷移時に最大、TIC は約 1/3) | — | △ |
| [V09] Verma, Shock Waves 19(3) 2009 △ / [VSH06] Verma, Stark, Haidn, AST 10(6) 2006 △ | DLR P6.2 | TOP 放物線。輪郭: 未確認 | 衝撃の非定常性と横荷重 | — | △ |
| [O04] Östlund, KTH 博士論文 2004 ✔ / [ODF04] Östlund, Damgaard, Frey, JPP 20(4) 2004 △ | Volvo Aero (VAC) のサブスケール試験 (S1 ほか) | TOC。輪郭: 非公開とみられる (公開資料を見つけられなかった) | 壁圧、横荷重。S1 で $n=p_a/p_c$ < 0.14 で FSS、それ以上で RSS | [OJ99] 標準 2 方程式モデルは剥離を外す。衝撃部の $P_k$ リミタと実現性補正で改善 (論文 3 の要旨 ✔) | ✖ 形状なし |
| [N03] Nguyen et al., FTC 71 2003 △ / [DN04] Deck & Nguyen, AIAA J 42(9) 2004 △ | LEA Poitiers | TOC、ヒステリシス域の横荷重。輪郭: 未確認 | 壁圧変動 | DES (Deck) | △ |
| Martelli, Jaunet et al. arXiv:2512.12317 (2025) △ / JFM 895 A29 (2020) △ | Institut PPRIME S150 | TIC、$D_t$=38 mm、出口径 97 mm ($\varepsilon$≈6.5 は本ノートで算出)、設計 M 3.5、空気約 260 K、NPR 6〜12。輪郭: 本文に記載なし | 壁圧変動 18 点 | DDES (SA ベース) は剥離位置を上流側に外す | △ 形状なし |
| **[CF60]** Campbell & Farley, NASA TN D-467 (1960) ✔ | NASA Lewis。直燃式加熱器による空気 1200 °F (一部非加熱)、NPR 4〜120 | 円錐、半角 15.3 / 25.0 / 29.0°、$D_t$=6.49 in、$\varepsilon$≈10 / 25 / 40 (同じノズルを切断)。寸法は図 1 → **寸法 (図)** | 推力比、壁静圧 (軸方向に一定間隔)、剥離圧の相関。加熱試験の大半で凝縮衝撃があり補正している | — | ○ 形状は単純。作動ガスが燃焼生成物を含み凝縮もあるので冷気検証の第一候補にはしない |
| [MG78] Morrisette & Goldberg, NASA TP-1207 (1978) ✔ | 文献データの集成 (円錐・輪郭・2D くさび) | — | 表 I に出典ごとの条件 | — | 経験式の根拠 (§3) |
| [NC73] Nave & Coffey, AIAA 73-1284 (NTRS 19740030190) △ | Rocketdyne J-2S 実機 ($\varepsilon$=40) と模型 | 輪郭: 非公開 | 剥離パターンの観察、横荷重モデル。RSS の命名元 [O04] | — | ✖ |
| [T09] Tomita et al., Shock Waves 19(3) 2009 △ | JAXA LE-7A サブスケール燃焼試験 | CTP (compressed truncated perfect) ほか。輪郭: 未確認 | 横荷重 | 数値解析は Shimizu et al., AIAA J 44(3) 2006 △ | ✖ 燃焼ガス・形状未確認 |

ユーザが挙げた「Ruf らの TIC / TOC / CTP」について。MSFC の公開資料で確認できたのは **TIC と PAR (TOP)** の 2 種で、CTP は JAXA の LE-7A 系 [T09] と Hoffman の設計法 (AIAA 85-1172 △) に出てくる。Ruf らが J-2X・SSME の模型で行った試験は JANNAF 論文 (2010) として [R10] に引用されているだけで、中身は確認できなかった。

### 2.2 2D 平面 CD ノズル

| 文献 | 機関・試験台・条件 | ノズル型・形状公開度 | 計測量・不確かさ | 先行 CFD | forge 適性 |
| --- | --- | --- | --- | --- | --- |
| **[MPR80]** Mason, Putnam, Re, NASA TP-1704 (1980, NTRS 19800020749) ✔ | NASA LaRC 16-Foot Transonic Tunnel の静的試験設備。空気 300 K、最大 1013 kPa。静止大気へ排気 | 2D-CD 5 種。A1・A2 ($A_e/A_t$=1.09) と B1・B2・B3 ($A_e/A_t$=1.80、$M_d$=2.08、NPR$_d$=8.81)。収縮・発散とも直線、スロートは円弧 ($r_c$=0.68 / 2.74 cm)、収縮角 $\theta$=20.84 / 22.33 / 20.42°、発散角 $\varepsilon$=10.85 / 11.24 / 10.85°。長さ・高さを図 2 の表に記載 → **寸法** | 推力比、流量係数。上下フラップの中心線・$y/(w_t/2)$=0.45・0.875 の 3 列と側壁中心線の静圧を表 I〜V に (B2 は NPR 1.96〜8.73 の 15 点)。B 系は NPR 約 2 でスロート直後から剥離。不確かさの数値: 未確認 | 本報告は 2D 非粘性理論との比較のみ。この B 系を使った RANS 比較は未確認 | **◎** 形状・壁圧とも表。剥離が出るのは NPR ≈ 2〜3 で、側壁列で 2 次元性を確かめられる。表はスキャン画像 |
| **[H98]** Hunter, AIAA 98-3107 (NTRS 20040090513) ✔ / [H04] JPP 20(3) 2004 △ | NASA LaRC 16-Foot 静的試験設備。乾燥空気 $T_0$≈530 °R、$p_a$=14.85 psi。NPR 1.25〜9.5 で 20 点 (うち 15 点が NPR < 6)。設計スロート Re 3.2×10⁶ | 2D-CD、$A_t$=4.317 in²、$A_e/A_t$=1.797、幅 3.990 in、NPR$_d$=8.78、$M_e$=2.07。点 A〜I の座標表 + 収縮円弧 2.000 in R・スロート円弧 0.625 in R・収縮 27.29°・発散 11.01° → **全点 (座標+半径)**。ガラス側壁 | 中心線と側壁から 0.40 in の列の静圧 (図のみ、ポート径 0.020 in)、推力効率、焦点シュリーレン。不確かさ: NPR と $p/p_{0j}$ ±0.28 %、$F/F_i$ ±0.4 %。流れの区分: NPR ≤ 1.8 は 3D・非定常の部分再付着剥離、NPR ≥ 2.0 は全離脱 (2D 的) で、1.8→2.0 で急変 | PAB3D (k-ε + 代数応力 SZL) の 2D 計算は NPR ≥ 2.4 で良好、NPR < 2.0 では全離脱を予測して外れた。ピーク推力効率 0.993 (実験 0.986、差は側壁摩擦)。[HV97] NPARC で 5 モデル比較 (Hunter の修士論文データ)、[HV98] 3D 計算で過膨張時にフラップ側から外気が流入し衝撃背後に側壁の二次流れ。[Ba11] SST と v2-f が最良 (二次的言及のみ △)、[TM19] Fluent で SST が最良 | **◎** 形状は表。壁圧は図の読み取り。低 NPR 側は 2D で合わないことが分かっているので合否から外す |
| [PZJ09] Papamoschou, Zill, Johnson, Shock Waves 19(3) 2009 △ / [PZ04] AIAA 2004-1111 △ | UC Irvine | 平面 CD、平板壁、$A_e/A_t$=1.0〜1.6。[XTP07] によると $A_e/A_t$=1.5 の発散角 3.89° | 非対称な λ 衝撃 (片側の足が大きい) を NPR > 1.4 で観測、中心線圧 | [XTP07] ✔ 2D 定常 RANS (Wilcox k-ω) で衝撃位置・壁圧が実験と良く合い、1.5 < NPR < 2.4 で非対称解 | ○ 非対称剥離の確認用。低 NPR に限られる |
| [BR05] Bourgoing & Reijasse, Shock Waves 14(4) 2005 △ | ONERA S8Ch | 平面ノズル、下流の第 2 スロートで背圧を制御。輪郭: 未確認 | 非対称衝撃、壁圧変動 (60 Hz 帯、200 Hz) ([XTP07] による) | Deck et al. AIAA 2002-406 (SA で初期値により非対称 2 解・対称 1 解) △ | △ |

TP-1704 の B 系と Hunter のノズルは別物である。面積比はほぼ同じだが、スロート半径 (2.74 cm vs 0.625 in = 1.59 cm) と収縮角 (22.33° vs 27.29°)、発散角 (11.24° vs 11.01°) が違う。

### 2.3 デュアルベル

手元の 7 本の評価は §2.3.1、それ以外は §2.3.2。

#### 2.3.1 手元の 7 本 (`papers/nozzle_design/dual_bell/`) の評価

| 資料 | 内容 | 形状公開度 | データ公開度 | 先行 CFD | forge 適性 |
| --- | --- | --- | --- | --- | --- |
| **[GS11a]** Génin & Stark, "Experimental investigation of the inflection geometry on dual bell nozzle flow behavior", AIAA 2011-5612 ✔ | 平面デュアルベル、変曲部 5 形状 (尖り・R15・R30・組合せ)。冷気 (P6.2、乾燥 N2、最大 6 MPa) と温風 (M11.1、空気 95 % + H2O 5 %、$T_0$ 600〜1000 K、最大 3 MPa) | 表 2: $R_{th}$=8 mm、$L_b/R_{th}$=16.25、$L_e/R_{th}$=15、$\varepsilon_b$=3.95、$\varepsilon_e$=7.97、$\alpha$=15°。ベースは $M_D$=2.8 の全長理想ノズル (平面)、延長は等圧 → **MOC 再構成** | 剥離位置対 NPR (図)、変曲近傍の壁圧 (図)、sneak transition の ΔNPR (表 3)、温風の遷移マップ (図)。冷気では延長部の壁圧勾配がわずかに負になり (3D 効果)、本来の遷移が起きなかった。600 K で H2O が凝縮し遷移 NPR が下がった。輪郭 1 は製作不良 | — (先行研究として Martelli et al. の数値解析を引用し、変曲域長さ $L_i/R_{th}$ を数値 1.8〜4、実験約 1 と対比) | ○ 平面・冷気で遷移しない点が検証上の弱み。変曲部の剥離 (sneak transition) の比較には使える |
| **[GGS11]** Génin, Gernoth, Stark, "Experimental and numerical study of heat flux in dual bell nozzles", AIAA 2011-5613 ✔ (JPP 29(1):21–26, 2013 △) | 上と同じ平面モデルの M11 温風試験 ($T_0$=1000 K、$p_0$ 1.5〜3 MPa) と熱流束。冷気シュリーレン | 図 2 (寸法は上と同じ系)。軸対称モデル (サーモグラフィ) は図 12 | NPR 20 / 24.5 / 29.8 の壁圧分布 (図)、剥離衝撃角対 NPR (図)、壁内温度 (1 mm・3 mm) と熱流束 | Fluent v12.1 密度ベース、SA、ハイブリッド 164,675 要素、y⁺ 最大 1.1。海面モードの壁圧は良く一致。NPR を上げると計算の剥離位置が実験より下流に出て、比較には NPR のずらしが要った | ○ SA が下流に外す実例 |
| **[SGWK12]** Stark, Génin, Wagner, Koschel, "The altitude adaptive dual bell nozzle" (2012、掲載会議名は未確認) ✔ | DLR の冷気試験 (DB1〜DB3) の総括、TAU による Vulcain 型デュアルベルの飛行模擬 | 形状なし (Vulcain 型の設計パラメータ表のみ) | 遷移予測式 (§3 の式)、延長長さ比 $L'$ と遷移 NPR・ヒステリシス・遷移時間の関係 (図)。ヒステリシス 20 % を確保するには $L' \ge 0.5$ | TAU (SA、凍結 H2/H2O) で海面モードの延長部壁圧が外気の最大 75 % まで下がる (外部流の吸い出し)。遷移は設計より低高度で起きた | △ 予測式と傾向の確認用 |
| **[SG16]** Stark & Génin, "Sea-level transitioning dual bell nozzles", Space Propulsion 2016 ✔ (CEAS Space J 9(3) 2017 △) | Ariane 5 ECA を想定したシステム検討 | 表 3 に $L_{inf}$・$\varepsilon_{inf}$・$\alpha_{inf}$・$L_{tot}$・$\varepsilon_e$ (MOC 設計の結果) | 流れの計測なし | — | ✖ 検証データではない。④ の設計例として参照 |
| **[HF94]** Horn & Fisher, "Dual-bell altitude compensating nozzles" (Penn State PERC annual report 1993, pp.140–147; NTRS 19940018584) ✔ | Rockwell (El Segundo) の高空模擬槽、空気 (12 lb/s、300 psia)。16:1 Rao ベース + 40:1 延長 4 種 (円錐・Rao・等圧・overturned) | 図 5・6 のみ (Rao 16:1 / 40:1 は再生成できるが延長部の形は図から読むしかない) | 推力係数対 NPR (天秤)、壁圧分布 (図 8)、シュリーレン。等圧・overturned 延長は 30 ms 未満で再現性よく遷移、円錐・Rao 延長は付着が不安定 | — | △ 図のみ |
| [J13] Jones, Bui, Ruf, NASA/TM-2013-216590 ✔ | F-15 による飛行試験の提案 | なし | MSFC NTF でデュアルベルと通常ベルを壁圧ポート付きで静的試験したと記述 (写真のみ) | MSFC の CFD (マッハ数コンタのみ) | ✖ データなし |
| [R16] NASA 20160008031 (MSFC 技術報告 2014、Ruf / Jones) ✔ | 上記飛行試験の進捗 | なし | 冷気試験データの出典として Ruf, McDaniels, Nguyen (JANNAF 2005) を引用 | — | ✖ データなし (JANNAF 論文は入手不可) |

7 本の中で形状を作り直せるのは [GS11a][GGS11] の平面モデルだけで、それも MOC 再構成が要る。軸対称で使えるものは手元に無かったので §2.3.2 で探した。

#### 2.3.2 それ以外のデュアルベル

| 文献 | 機関・条件 | 形状公開度 | 計測量 | 先行 CFD | forge 適性 |
| --- | --- | --- | --- | --- | --- |
| **[NGS07]** Nürnberger-Génin & Stark, EUCASS 2007 ✔ / **[NGS10]** JPP 26(3):497–502, 2010 △ / [NGS09] Shock Waves 19(3) 2009 △ | DLR P6.2 水平リグ、N2、$p_a$=1 bar、供給圧を 0.7〜2.2 bar/s で上下 | 軸対称。$R^*$=10 mm、$\varepsilon_1$=11.3、$\varepsilon_2$=27.1、$L_{base}/R^*$=6.2、$L_{ext}/R^*$=8.3、$\alpha$=9°。ベース TIC、延長は等圧 (MOC) → **MOC 再構成** (ベース TIC の設計 Mach は未記載。$\varepsilon_b$ と $L_b$ の 2 条件から逆算する) | 壁圧 (LF / HF Kulite)、シュリーレン 2 kHz。NPR 30 で変曲部剥離、48 で剥離域の勾配が減り、51 で全流れ。延長部の剥離域圧力は $p_a$ の 92〜98 %。ヒステリシス平均約 20 %。$T_0$ < 260 K で模型表面に着氷し遷移が早まる。[NGS10] は 3 形状と延長の段階的切断で、遷移 NPR が延長長さに線形に依存すると報告 (抄録) | — | **○** 軸対称・冷気で一番条件がそろう。形状再構成の検証が先に要る |
| [SchG16] Schneider & Génin, JPP 32(5):1212–1219, 2016 △ (AIAA 2015-4219 △) | P6.2 冷気試験 (正の壁圧勾配の延長をもつモデル) | 未確認 | — | 2D 軸対称の RANS・URANS。SA と RSM が良好。遷移 NPR を 1 %、ヒステリシスを約 10 % の精度で予測、遷移時間は同じ桁 (抄録) | 先行 CFD の目安 |
| [VSNH10] Verma et al., Shock Waves 20(3) 2010 △ / [VSGH11] Shock Waves 21(2) 2011 △ / [VSH12] JPP 28(6) 2012 △ / [VSH13] JPP 29(3) 2013 △ | DLR P6.2 (海面と高空模擬槽) | 未確認 | 遷移時の剥離特性、ガス密度の影響、Re 数の影響 (Re が下がると変曲域が広がり遷移 NPR が上がる。高空模擬槽の試験は寸法効果が強い) | — | △ (Re 依存の注意として参照) |
| [HTHF02] Hagemann, Terhardt, Haeseler, Frey, JPP 18(1) 2002 △ | TEHORA / TEKAN の冷気・燃焼ガス試験 | 未確認 | 延長形状による遷移の速さ、熱伝達 | — | △ |
| [GW04] Gross & Weiland, JPP 20(3) 2004 △ | (CFD 論文) | — | — | TIC 1・TOP 2・デュアルベル 3 のサブスケール冷気を計算し実験と比較、手法の欠点も指摘 (抄録) | 先行 CFD の目安 |
| [NOM05] Nasuti, Onofri, Martelli, JPP 21(2) 2005 △ / Martelli et al., AIAA J 45(3) 2007 △ | (CFD 論文) | MOC で設計した形状 | — | 延長形状のわずかな違いで遷移挙動が大きく変わる (抄録) | 設計側の参考 |
| [Re11] Reijasse et al., Progress in Propulsion Physics 2 (2011), doi 10.1051/eucass/201102655 △ | ONERA の冷気設備 | 未確認 (EUCASS proceedings サイトが自動取得を拒否) | 壁圧、推力 | — | △ |
| [To06] Tomita et al., AIAA 2006-4375 △ / [Ki09] Kimura et al., AIAA 2009-5149 △ / Proschanka et al., JPP 28 (2012) △ | JAXA 角田・大阪大 | 未確認 | 冷気試験、遷移・横荷重、フィルム冷却による遷移制御 | — | △ |

### 2.4 SERN

| 文献 | 機関・条件 | 形状公開度 | 計測量 | 先行 CFD | forge 適性 |
| --- | --- | --- | --- | --- | --- |
| **[RL84]** Re & Leavitt, NASA TM-86270 (1984, NTRS 19850006488) ✔ | NASA LaRC 16-Foot 静的試験設備、空気約 300 K、最大約 10 atm、NPR 2〜10、静止大気 | 43 形状 (ランプ長・ランプ初期角・ランプ弦角・下フラップ長・下フラップ角の組合せ)。流路幅 10.16 cm、スロート高さ 2.54 cm (公称)、スロート縦横比 4.0、上フラップのスロート半径 1.65 cm、ランプ終端角 2.5°。上フラップ 15 種・下フラップ 9 種の無次元座標を表 II・III、側壁 7 種を図 4 → **全点** | 推力比・合力・推力偏向角・流量係数・ピッチングモーメント (図 5)、上下フラップ中心線の $p/p_{t,j}$ (表 V)。NPR 4 でランプ上に衝撃が立ち、後方で圧力が平坦になる (剥離)。高面積比の DQ3・HQ4 は低 NPR で内部衝撃 | 同設備の SERN を PAB3D で計算した Carlson & Abdol-Hamid, AIAA 93-2571 △ (中身未確認) | **◎/○** 形状と壁圧は表 (スキャン)。剥離の強さは中程度で、側壁による 3D 効果に注意 |
| **[DA99]** Deere & Asbury, NASA/TP-1999-209138 (NTRS 19990052644) ✔ | NASA LaRC Jet Exit Test Facility。NPR ≤ 13 | 平行移動スロートの高速 SERN。スロート位置 3 種 (設計 NPR ≈ 9・42・102) × 凹・凸ランプ。表 1 に主要寸法 ($x_t$, $l_c$, $l_r$, $h_t$, $h_e$, 角度)、ランプ形状は図 → **寸法 (一部図)** | 推力・偏向角、ランプ上の静圧 29 点、paint-oil、焦点シュリーレン。中・高 Mach 形態は全試験点で強い過膨張、3D 剥離を観測 | PAB3D (k-ε) 2D。低 Mach 形態は推力比を 1.5 % 以内。高 Mach 形態の NPR 10・13 はランプ上の大きな剥離で収束せず定性比較のみ。2D では 3D 剥離を表せず推力効率を過大評価 | ○ スクラムジェット SERN の強い過膨張の唯一の公開例。形状の一部が図 |
| **[YXYM15]** Yu Y., Xu J., Yu K., Mo J., Shock Waves 25(6):623–633, 2015 △。実験の記述は **[YT23]** Yu T. et al., JAFM 16(11):2249–2262, 2023 ✔ | NUAA の吹き出し / 吸い込み風洞 (高圧 0.8〜0.9 MPa、真空 3〜5 kPa、最大 0.8 kg/s)。空気 300 K。NPR 2.3〜3.0 付近 | 2D SERN、膨張角 25° のランプ + フラップ、$A_t$=1200 mm²、$A_e/A_t$=2.896、幅 60 mm、光学ガラス側壁、上流 30 mm の円→矩形遷移部。形状は図 2(a) → **図 (寸法の充足度は未確認)** | Z 型シュリーレン (200 mm)、高速度カメラ、ランプ壁圧 (FS の 0.05 %)。NPR 3.0 で RSS (ランプ側)、NPR 2.3 で RSS (フラップ側)、パターン遷移 | Fluent (RNG k-ε + 壁関数)、OpenFOAM rhoCentralFoam / sonicFoam。RSS の 2 パターンを再現。sonicFoam は遷移を 0.0773 NPR 早く起こし、剥離泡後の圧力上昇を外した | ○/△ RSS のランプ側・フラップ側の再現確認用。形状は図から |
| [Tan09] Tan, Jin, Du, Ji, Xu, Journal of Propulsion Technology (推進技術) 30(3):292–296, 2009 △ | 北京航空航天大。高圧空気、NPR 4〜10、風なし | 2D SERN。未確認 | 壁静圧。FSS と衝撃誘起剥離、軸推力係数 0.78〜0.87 (抄録) | 計算あり | △ 中国語・形状未確認 |
| [CRB92] Capone, Re, Bare, NASA TP-3240 (1992, NTRS 19920024949) △ | NASA LaRC 16-Foot 遷音速風洞、M 0.60〜1.20 | 未確認 | SERN の外部流つき性能 | — | 外部流の影響を見るときの候補 (未読) |
| [GHK00] Gruhn, Henckels, Kirschstein, AST 4(8) 2000 △ | DLR Köln 極超音速風洞 | 未確認 | フラップ輪郭と SERN 性能 | — | △ |

ユーザが挙げた「Lv et al. AST 2019 (低 NPR の RSS / FSS、外部流の影響)」は特定できなかった。著者に Lv Z. を含む AST 2019 の論文は Yu K., Xu J., Lv Z., Song G., "Inverse design methodology on a single expansion ramp nozzle for scramjets", AST 92:9–19 [YXLS19] △ で、MOC による逆設計の論文 (剥離実験ではない)。NUAA 系で外部流と剥離を扱うのは Yu Y. の AST 98:105674 (2020) △ と Aerospace 10(11):958 (2023) △ で、どちらも数値計算が主体である (2023 は実験で手法を検証したと抄録に書かれている)。

---

## 3. 剥離判定の経験式

記号: $p_i$ = 剥離直前の壁圧 (文献により $p_s$, $p_{sep}$)、$p_a$ = 外気圧、$p_p$ = 剥離後の平坦部圧、$M_i$ = 剥離直前の壁面 Mach 数、$p_0$ = 全圧。

| 名称 | 式 | 根拠データ・適用範囲 | 出典と検証状態 |
| --- | --- | --- | --- |
| Summerfield (1954) | $p_i/p_a \approx 0.4$ | 1940 年代の円錐ノズル | 原典 (Jet Propulsion 24, 1954) は未読 △。[O04] 式 (19) と [MG78] が引用 ✔ |
| Schilling (1962) | $p_i/p_a = k_1 (p_0/p_a)^{k_2}$。輪郭: $k_1$=0.582, $k_2$=−0.195。円錐: $k_1$=0.541, $k_2$=−0.136 | 円錐・TIC | 原典 (修士論文) 未読。[O04] 式 (20) ✔ |
| Kalt & Badal (1965) | 上式で $k_1$=2/3, $k_2$=−0.2 | 円錐 | 原典 (J. Spacecraft 2(3)) △。[O04] と [Sch73] 式 (10) が同じ係数 ✔ |
| Schmucker (1973) | $p_i/p_a = (1.88\,M_i - 1)^{-0.64}$ | 化学ロケットの剥離データの平均 (冷気・燃焼ガス) | [Sch73] 式 (11) ✔。欧州で最も広く使われた式と [O04] |
| Reshotko–Tucker (1955) | 剥離衝撃前後の Mach 比 0.762 (零圧力勾配) から $p_i/p_a$ を算出 | 半角 ≳ 9〜10° の円錐で良く合う。低角円錐・輪郭ノズルの高 $M_s$ 側では外れる | [MG78] ✔ |
| Morrisette & Goldberg (1978)、輪郭ノズル | $p_s/p_a = 1.082 - 0.363\,M_s + 0.0386\,M_s^2$ | $2.4 \le M_s \le 4.5$ (Roschke & Massier を除く 2 次近似)。層流・遷移域のデータを除外したうえでの相関 | [MG78] 式 (10) ✔ |
| 同、出口近傍の適用限界 | $A_s/A^* = 0.8\,(A_e/A^* - 1) + 1$ より下流では $p_s/p_a$ が増加側にずれる (相関の対象外) | 出口 Mach 2 以下では式 (7) $A_s/A^*=0.8\,A_e/A^*$ より合う | [MG78] 式 (7)(8) ✔ |
| 同、剥離しない上限 | $p_t/p_a = 1 + 0.39\,p_t/p_e$ を超えると剥離しない (出口に斜め衝撃が立つだけ) | — | [MG78] 式 (9) ✔ |
| Zukoski (1967) | $p_p/p_i = 1 + 0.5\,M_i$ (剥離点の圧力上昇のみ) | $M_i$=1.4〜6.0、$Re_{\delta} \ge 10^5$、空気 ($\gamma$=1.4) の段差流 | 原典 AIAA J 5(10) △。[O04] 式 (22) ✔ |
| Frey (2001) | $p_i/p_p = \left[1 + \tfrac{2\gamma}{\gamma+1}\left(M_i^2 \sin^2\beta - 1\right)\right]^{-1}$、$\beta = -4.7\,M_i + 44.5$ [°] | $2.5 \le M_i \le 4.5$。$M_i \approx 4.8$ で極小を持ち、それ以上で傾向が逆転 | [O04] 式 (24) ✔ (原典は Frey の博士論文 2001、未読) |
| Östlund | 斜め衝撃関係で $p_i/p_p$ を与え、$\beta = -3.764\,M_i + 42.878$、$\theta = 1.678\,M_i + 9.347$ [°] | $2.5 \le M_i \le 4.5$、VAC サブスケール試験 | [O04] 式 (23) ✔ (式全体は原典で確認すること) |
| Stark | $p_{sep}/p_a = 1/M_{sep}$ | DLR の冷気・燃焼ガス試験と文献の剥離データベース [S05] | [S07] 図 4 の凡例で式を確認 ✔。原典 [S05] は抄録のみ △ |
| Stark のデュアルベル遷移式 | $\mathrm{NPR}_{tr} = \dfrac{1}{M_e}\left(1 + \dfrac{\gamma-1}{2} M_e^2\right)^{\gamma/(\gamma-1)}$ ($M_e$ = 等圧延長の壁面 Mach 数) | 上の剥離式を延長部に当てはめたもの。延長を全長で使う場合だけ成り立ち、延長を切り詰めると遷移 NPR は線形に上がる | [SGWK12] 式 (1) ✔ |

Frey と Östlund の式は $p_i/p_p$ (剥離点の圧力上昇) を与えるので、$p_i/p_a$ にするには平坦部圧 $p_p/p_a$ が別に要る ([JG15] は Frey の式を $p_p/p_a$=0.9 で使った)。Lawrence と Arens & Spiegler が 1960 年代に指摘したとおり、$p_i/p_a$ を直接与える式は「剥離点の圧力上昇」と「剥離域で外気へ回復する圧力上昇」を一つにまとめており、ばらつきが大きい ([O04] §5.1.2)。NASA は 1976 年に Crocco–Probstein の式に 20 % の余裕を付けて使うよう推奨した ([O04]、原典 NASA SP-8120 は未読)。

[JG15] は DLR の楕円化 TIC (幾何 4: TIC-53 を変形、NPR 20・40) の 1 断面で 4 式を比べ、Summerfield と Schmucker は剥離位置を上流に、Stark (低 NPR) と Frey は実験に近く出たと報告している。

---

## 4. 先行 CFD 検証の知見

### 4.1 軸対称 TIC (FSS)

- **FSCD/ATAC test case 1A** [S07] ✔: SST 系の多くは剥離位置を上流側に予測し (標準 SST で −9〜−17 %)、SA と標準 k-ω (Wilcox、Sarkar 補正) は −2.4〜+4.1 %。SST の改良版はほとんど改善せず、Östlund の改良 (VAC 冷気データで調整) だけが +3.3 %。著者の解釈では、乱流モデルは境界層だけでなく噴流せん断層の巻き込み能力を通じて剥離域の逆流量を決め、逆流域の最大壁圧を 1〜4 % 高く出す解ほど剥離が早まる。斜め衝撃角も大きめに出る。マッハディスク背後の閉じ込め渦と強く湾曲したマッハディスクは格子を細かくすると消えたので、格子依存の人工物と結論している。
- **DLR 楕円化 TIC** [JG15] ✔: TAU の SA は剥離位置を実験より**下流**に、剥離域の壁圧を低く出した。RSM で両方とも改善した。著者は未変形形状と先行研究でも同じ傾向だったと書いている。
- **Östlund & Jaran** [OJ99] (要旨 ✔): 補正なしの標準 2 方程式モデルは VAC ノズルの剥離特性を外した。原因は衝撃での乱流エネルギー生成 $P_k$ の過大で、$P_k$ リミタと弱非線形の実現性補正で改善した。
- **PPRIME TIC の DDES** (arXiv 2025 △): SA ベースの DDES でも剥離位置を上流側に外した (格子起因の可能性を著者が指摘)。

同じノズルで SST と SA を並べたのは TC1A だけだが、「SST は上流寄り、SA は下流寄り」の傾向が複数の報告で見える。forge は SST を使うので、TC1A で forge の SST が文献の SST と同程度 (−10 % 前後) に出るかをまず確かめ、乱流モデルの偏りと forge 固有の誤差を分ける。

### 4.2 TOP / TOC の FSS→RSS

- MSFC の PAR [R09][R10] ✔: 設計段階の CFD が RSS 遷移を NPR 約 24 と予測し、実験は 23.6〜23.8 で遷移した (乱流モデルは記載なし)。RSS は NPR 一定でも強く非定常で、静圧の時間平均を代表値として載せている。
- [FH00][HFK02] (抄録 △): RSS は cap-shock パターン (内部衝撃の軸上での逆マッハ反射) から生じる。
- 定常 RANS で RSS を扱う場合は、増圧と減圧の両系列 (前の NPR の解からの再始動) を回さないとヒステリシスを再現できない。

### 4.3 2D 平面 CD

- Hunter [H98] ✔: 2D の PAB3D は NPR ≥ 2.4 (全離脱・2D 的) で壁圧と推力効率が実験と良く合い、NPR < 2 (3D の部分再付着) では全離脱を予測して外れた。NPR 3.0 では SZL 代数応力モデルが衝撃位置と剥離を最もよく再現した。
- [HV98] (抄録 △): 3D 計算で、過膨張時の外気流入は主にフラップ側から起き、衝撃背後で側壁を横切る強い二次流れがある。2D 計算と 3D 計算の使い分けを推奨している。
- [XTP07] ✔: 2D 定常 RANS (Wilcox k-ω) で $A_e/A_t$=1.5 の衝撃位置・壁圧が合い、1.5 < NPR < 2.4 で非対称解が出た。外部領域は 10 スロート高さ、外周に M=0.1 の流れを与えた。
- [Ba11] (二次言及 △)・[TM19] (抄録 △): Hunter のデータで SST が最良とする報告が複数ある。

### 4.4 デュアルベル

- [GGS11] ✔ (Fluent SA): 海面モード (変曲部で剥離) の壁圧は一致。NPR を上げて剥離点が延長部に入ると、計算の剥離位置が下流に出た。
- [SchG16] (抄録 △): 軸対称 RANS / URANS で SA と RSM が良く、遷移 NPR を 1 %、ヒステリシスを約 10 % で予測。遷移 NPR は供給圧の上げ速度で変わる。
- [SGWK12] ✔: 定常 CFD は遷移の過渡を扱いにくく、設計 (MOC) より低い高度で遷移した。外部流があると延長部の剥離域圧力が外気の 75 % まで下がる (通常のノズルの剥離域は約 90 %)。

### 4.5 SERN

- [DA99] ✔: 2D の k-ε 計算は設計点近傍で推力比を 1.5 % 以内で再現したが、強い過膨張 (高 Mach 形態、NPR 10・13) ではランプ上の大きな剥離で収束しなかった。実験の剥離は 3D で、2D 計算は推力効率を過大評価すると著者が明記している。
- [YT23] ✔: RANS (RNG k-ε、OpenFOAM の 2 ソルバ) で NPR 3.0 の RSS (ランプ側) と NPR 2.3 の RSS (フラップ側) の両方を再現。ソルバの違いで遷移 NPR が 0.08 ずれた。

### 4.6 外部領域・境界条件の慣行 (読んだ文献から)

| 文献 | 外部領域の大きさ | 外部の境界条件 |
| --- | --- | --- |
| [H98] | 出口から下流 30、上流 25、法線方向 25 スロート高さ | 左面に $M_a$=0.025 の亜音速流入、上面は特性 BC、右面は局所 Mach で圧力一定 / 外挿を切替 |
| [XTP07] | 10 スロート高さ | 外周に $p_a$ と M=0.1 |
| [GGS11] | 未記載 | farfield、M=0.01 |
| [DA99] | 遠方境界 106 スロート高さ (低 Mach 形態) | M=0.05 (高 Mach 形態は安定化のため M=0.1) |

静止大気を小さい並行流 (M 0.01〜0.1) で近似するのが一般的で、安定化のためと明記した文献もある。forge の [case/27](../../case/27.axi_nozzle_plume/README.md) はプルーム領域つきの軸対称計算で、プルームの衝撃セルが非定常なため残差が約 2×10⁻³ で頭打ちになったと記録している (ノズル内部は定常)。

---

## 5. forge の検証へ持ち込むときの注意 (候補段階)

採否は plan で決める。ここでは文献から読み取れる注意だけを書く。

1. **剥離位置の定義を実験に合わせる**。TC1A は「最大壁圧勾配の接線と全流れ MOC 壁圧の交点」、Ruf は「ポートごとの $p_{sep}/p_a$」で定義している。forge で $\tau_w$=0 の点を使うと定義が違うので、壁圧ベースの同じ定義で比べる。
2. **2D で合わない領域を合否から外す**。Hunter の NPR < 2、[DA99] の強い過膨張、[GS11a] の平面延長部の負勾配は 3D 効果である。TP-1704 は側壁列の壁圧で 2 次元性を確かめられる。
3. **定常性の判定**。剥離位置・マッハディスク位置はリミットサイクルになりうる (Hunter NPR ≤ 1.8 は実験でも非定常、RSS は強く非定常)。AGENTS.md の準定常確認 (`check_quasisteady.py`) を剥離位置に当て、振動なら平均と振幅で比べる。
4. **ヒステリシスのある量** (RSS 遷移、デュアルベル遷移) は、増圧系列と減圧系列を前の NPR の収束解から再始動して回す。同一メッシュの再始動は `restart_field.py` を使う。
5. **乱流モデルの偏り**。TC1A では標準 SST が −9〜−17 % に外れるのが文献の基準になる。forge の SST の選択肢 (`dilatationCorrection`、`katoLaunder` など。[`procedures/recommended-settings.md`](../../procedures/recommended-settings.md)、[`turbulence-sst-consistency-options.md`](../../plans/active/turbulence-sst-consistency-options.md)) が剥離位置に効くかは TC1A で感度を見る対象になる。[OJ99] が原因に挙げた「衝撃での $P_k$ 過大」に対応するのは生成項側の選択肢 (`katoLaunder` など) で、効くかどうかは未確認。SA など別モデルとの差は [procedures/su2-cross-check.md](../../procedures/su2-cross-check.md) の手順で SU2 (SST と SA) を同じメッシュ・BC で回して切り分けられる。
6. **作動ガス**。TC1A と P6.2 系は N2、MSFC と LaRC は空気。M11 の温風 (H2O 5 %、600 K で凝縮) と [CF60] の直燃加熱空気 (凝縮衝撃あり) は冷気検証から外す。
7. **表データの転記**。TC1A の表 3 は本文テキストで取り出せる (重複点と同一値の 2 列は §2.1 の注記どおり扱う)。TP-1704・TM-86270・Ruf 付録 A はスキャン画像なので、転記値を図と突き合わせる検査が要る。

---

## 6. 未確認事項

| 項目 | 状態 |
| --- | --- |
| Stark 2013 overview [S13] の本文 (剥離データベースの中身) | 抄録のみ。有償 |
| Hadjadj & Onofri 2009 [HO09] の本文 | 書誌のみ。有償 |
| Östlund & Muhammad-Klingmann 2005 [OMK05] の本文 | 同名の博士論文 [O04] で代用。本文の差は未確認 |
| Balabel 2011 [Ba11]、Yaravintelimath 2016 [Y16]、Gross & Weiland 2004 [GW04] の結果の数値 | 本文未読 |
| FSCD/ATAC workshop (2006, Noordwijk) の他の test case (1C: DLR TIC の横荷重など) と、TC1A の NPR 25.25 以外のデータ | proceedings が公開されていない |
| DLR の PAR/TOP、DB1〜DB3、VAC S1、ONERA の TIC/TOC、LEA の TOC、JAXA の CTP の輪郭座標 | 公開資料を見つけられなかった |
| MSFC PAR の出口角・全長 (図から約 $x/r^*$=12 と読める) | 図のみ |
| MSFC TIC の MOC 再構成が実機輪郭を再現するか | 設計壁圧曲線 (図 8) との照合が要る。実機はスロート直後に段差 |
| MSFC のデュアルベル冷気試験データ (Ruf, McDaniels, Nguyen, JANNAF 2005) | 入手不可 |
| NUAA SERN の形状図に全寸法が入っているか ([YXYM15] 図・[YT23] 図 2) | 図の精読が未了 |
| 「Lv et al. AST 2019」 | 該当論文を特定できず (§2.4) |
| TP-1704 のノズル幅と計測不確かさ、TM-86270 の不確かさ | 本文の該当箇所を未確認 |
| Summerfield 1954、Schilling 1962、Frey 2001 (博士論文)、NASA SP-8120 の原典 | 未読 (二次資料の式を使った) |
| Nave & Coffey 1973 の本文 | 抄録のみ |
| [SGWK12] の掲載会議 | PDF にロシア語テンプレートの痕跡があり、会議名を特定できず |

---

## 7. 先行サーベイとの関係

[nozzle-optimization-tool-survey.md](nozzle-optimization-tool-survey.md) の B4.5 は剥離制約を「$\tau_w$ 符号反転 / Summerfield $p_w/p_a \lesssim 0.35$〜0.4」としている。§3 のとおり Summerfield 型の定数は円錐の古いデータ由来で、輪郭ノズル・高 $M_i$ では外れる ([MG78] 式 (10)、[O04] 図 34)。剥離制約の定義を見直すときは本ノートの §3 と §5-1 (定義の違い) を参照する。同サーベイ付録 10 (a) の「デュアルベルの一次確認」は §2.3 で一部進めた (遷移ヒステリシス幅は約 20 % [NGS07]、横荷重の定量は [GS11b] △ で未読)。

---

## 8. 保存した PDF (`/home/sano/work/forge/papers/nozzle_separation/`、git 管理外)

| ファイル | 文献 |
| --- | --- |
| `Stark_Hagemann_2007_EUCASS_numerical_prediction_separated_nozzle_flows_TIC_testcase.pdf` | [S07] (DLR elib 49262) |
| `Ruf_McDaniels_Brown_2010_TIC_and_parabolic_side_load_details_NTRS20100033270.pdf` | [R10] |
| `Ruf_McDaniels_Brown_2009_AIAA-2009-4856_side_load_testing_MSFC_NTRS20100017649.pdf` | [R09] |
| `Mason_Putnam_Re_1980_NASA-TP-1704_throat_contouring_2DCD_NTRS19800020749.pdf` | [MPR80] |
| `Hunter_1998_AIAA-98-3107_separated_nozzle_flows_NTRS20040090513.pdf` | [H98] |
| `Xiao_Tsai_Papamoschou_2007_AIAAJ45-3_numerical_nozzle_separation.pdf` | [XTP07] (UCI 著者サイト) |
| `Re_Leavitt_1984_SERN_static_internal_performance_NTRS19850006488.pdf` | [RL84] |
| `Deere_Asbury_1999_translating_throat_SERN_exp_CFD_NTRS19990052644.pdf` | [DA99] |
| `Yu_T_et_al_2023_JAFM16-11_OpenFOAM_SERN_separation_NUAA_exp.pdf` | [YT23] |
| `Nurnberger-Genin_Stark_2007_EUCASS_experimental_study_dual_bell.pdf` | [NGS07] (DLR elib 63529) |
| `Jack_Genin_2015_AIAA-2015-4154_ovalized_TIC_numerical.pdf` | [JG15] (DLR elib 98160) |
| `Ostlund_2004_PhD_KTH_supersonic_flow_separation_rocket_nozzles.pdf` | [O04] (21 MB) |
| `Morrisette_Goldberg_1978_NASA-TP-1207_separation_criteria_NTRS19780022611.pdf` | [MG78] |
| `Schmucker_1973transl1984_NASA-TM-77396_overexpanded_nozzles_part1_separation_NTRS19840011402.pdf` | [Sch73] |
| `Campbell_Farley_1960_NASA-TN-D-467_conical_CD_nozzles_NTRS19980227347.pdf` | [CF60] |

デュアルベルの 7 本は従来どおり `papers/nozzle_design/dual_bell/`。

---

## 9. 文献 (検証状態つき)

### ベル・円錐 (軸対称)

- [S07] ✔ R. Stark, G. Hagemann, "Current status of numerical flow prediction for separated nozzle flows," 2nd EUCASS, Brussels, 2007. DLR elib 49262 (`https://elib.dlr.de/49262/`).
- [SW06] △ R. Stark, B. Wagner, "Experimental Flow Investigation of a Truncated Ideal Contour Nozzle," AIAA 2006-5208, doi:10.2514/6.2006-5208 (抄録)。
- [SW09] △ R. Stark, B. Wagner, "Experimental study of boundary layer separation in truncated ideal contour nozzles," Shock Waves 19(3):185–191, 2009, doi:10.1007/s00193-008-0174-6 (書誌)。
- [S05] △ R. Stark, "Flow Separation in Rocket Nozzles, a Simple Criteria," AIAA 2005-3940, doi:10.2514/6.2005-3940 (抄録)。
- [S13] △ R. Stark, "Flow Separation in Rocket Nozzles – An Overview," AIAA 2013-3840, doi:10.2514/6.2013-3840 (抄録)。
- [R09] ✔ J. H. Ruf, D. M. McDaniels, A. M. Brown, "Nozzle Side Load Testing and Analysis at Marshall Space Flight Center," AIAA 2009-4856, NTRS 20100017649.
- [R10] ✔ J. H. Ruf, D. M. McDaniels, A. M. Brown, "Details of Side Load Test Data and Analysis for a Truncated Ideal Contour Nozzle and a Parabolic Contour Nozzle," AIAA 2010-6813, doi:10.2514/6.2010-6813, NTRS 20100033270.
- [Re01] △ P. Reijasse et al., "Flow separation experimental analysis in overexpanded subscale rocket-nozzles," AIAA 2001-3556, doi:10.2514/6.2001-3556 (抄録)。
- [F00] △ M. Frey, R. Stark, H. K. Ciezki, F. Quessard, W. Kwan, "Subscale nozzle testing at the P6.2 test stand," AIAA 2000-3777, doi:10.2514/6.2000-3777 (抄録)。
- [FH00] △ M. Frey, G. Hagemann, "Restricted Shock Separation in Rocket Nozzles," JPP 16(3):478–484, 2000, doi:10.2514/2.5593 (抄録)。
- [HFK02] △ G. Hagemann, M. Frey, W. Koschel, "Appearance of Restricted Shock Separation in Rocket Nozzles," JPP 18(3):577–584, 2002, doi:10.2514/2.5971 (抄録)。
- [V09] △ S. B. Verma, "Shock unsteadiness in a thrust optimized parabolic nozzle," Shock Waves 19(3):193–212, 2009, doi:10.1007/s00193-008-0180-8 (書誌)。
- [VSH06] △ S. B. Verma, R. Stark, O. Haidn, "Relation between shock unsteadiness and the origin of side-loads inside a thrust optimized parabolic rocket nozzle," AST 10(6):474–483, 2006, doi:10.1016/j.ast.2006.06.004 (書誌)。
- [O04] ✔ J. Östlund, "Supersonic flow separation with application to rocket engine nozzles," PhD thesis, KTH, TRITA-MEK 2004:07.
- [OMK05] △ J. Östlund, B. Muhammad-Klingmann, "Supersonic Flow Separation with Application to Rocket Engine Nozzles," Appl. Mech. Rev. 58(3):143–177, 2005, doi:10.1115/1.1894402 (抄録)。
- [ODF04] △ J. Östlund, T. Damgaard, M. Frey, "Side-Load Phenomena in Highly Overexpanded Rocket Nozzles," JPP 20(4):695–704, 2004, doi:10.2514/1.3059 (書誌)。
- [OJ99] △ J. Östlund, M. Jaran, "Assessment of turbulence models in overexpanded rocket nozzle flow simulations," AIAA 99-2583, doi:10.2514/6.1999-2583 (要旨は [O04] ✔)。
- [HO09] △ A. Hadjadj, M. Onofri, "Nozzle flow separation," Shock Waves 19(3):163–169, 2009, doi:10.1007/s00193-009-0209-7 (書誌)。
- [NC73] △ L. H. Nave, G. A. Coffey, "Sea level side loads in high-area-ratio rocket engines," AIAA 73-1284, doi:10.2514/6.1973-1284, NTRS 19740030190 (抄録)。
- [N03] △ A. T. Nguyen, H. Deniau, S. Girard, T. Alziary de Roquefort, "Unsteadiness of Flow Separation and End-Effects Regime in a Thrust-Optimized Contour Rocket Nozzle," Flow Turb. Combust. 71:161–181, 2003, doi:10.1023/B:APPL.0000014927.61427.ad (書誌)。
- [DN04] △ S. Deck, A. T. Nguyen, "Unsteady Side Loads in a Thrust-Optimized Contour Nozzle at Hysteresis Regime," AIAA J 42(9):1878–1888, 2004, doi:10.2514/1.2425 (書誌)。
- [DG02] △ S. Deck, P. Guillen, "Numerical Simulation of Side Loads in an Ideal Truncated Nozzle," JPP 18(2):261–269, 2002, doi:10.2514/2.5965 (抄録)。
- [JG15] ✔ S. Jack, C. Génin, "Flow Separation Study in Stiff Ovalized Rocket Nozzles, Part II: Numerical Approach," AIAA 2015-4154, doi:10.2514/6.2015-4154, DLR elib 98160.
- [T09] △ T. Tomita et al., "Experimental evaluation of side-loads in LE-7A prototype engine nozzle," Shock Waves 19(3):213–228, 2009, doi:10.1007/s00193-009-0191-0 (書誌)。
- △ E. Martelli, V. Jaunet, G. Della Posta, M. Bernardini, "Experimental and Numerical Analysis of the Intermittentency in a Nozzle Overexpanded-Flow," arXiv:2512.12317, 2025 (HTML 版の要約)。E. Martelli et al., JFM 895 A29, 2020 (書誌)。
- [CF60] ✔ C. E. Campbell, J. M. Farley, "Performance of Several Conical Convergent-Divergent Rocket-Type Exhaust Nozzles," NASA TN D-467, 1960, NTRS 19980227347.

### 2D 平面 CD

- [MPR80] ✔ M. L. Mason, L. E. Putnam, R. J. Re, "The Effect of Throat Contouring on Two-Dimensional Converging-Diverging Nozzles at Static Conditions," NASA TP-1704, 1980, NTRS 19800020749.
- [H98] ✔ C. A. Hunter, "Experimental, Theoretical, and Computational Investigation of Separated Nozzle Flows," AIAA 98-3107, doi:10.2514/6.1998-3107, NTRS 20040090513.
- [H04] △ C. A. Hunter, "Experimental Investigation of Separated Nozzle Flows," JPP 20(3):527–532, 2004, doi:10.2514/1.4612 (抄録)。
- [HV97] △ A. Hamed, C. Vogiatzis, "Overexpanded Two-Dimensional-Convergent-Divergent Nozzle Flow Simulations, Assessment of Turbulence Models," JPP 13(3):444–445, 1997, doi:10.2514/2.5183 (冒頭部のみ)。
- [HV98] △ A. Hamed, C. Vogiatzis, "Overexpanded Two-Dimensional Convergent-Divergent Nozzle Performance, Effects of Three-Dimensional Flow Interactions," JPP 14(2):234–240, 1998, doi:10.2514/2.5272 (抄録)。
- [Ba11] △ A. Balabel et al., "Assessment of turbulence modeling for gas flow in two-dimensional convergent–divergent rocket nozzle," Appl. Math. Model. 35(7):3408–3422, 2011, doi:10.1016/j.apm.2011.01.013 (書誌。結果は二次的言及のみ)。
- [TM19] △ S. L. B. Tolentino Masgo, "Evaluation of Turbulence Models for the Air Flow in a Planar Nozzle," Ingenius 22:25–37, 2019, doi:10.17163/ings.n22.2019.03 (抄録)。
- [Y16] △ A. Yaravintelimath, B. N. Raghunandan, J. A. Moríñigo, "Numerical prediction of nozzle flow separation: Issue of turbulence modeling," AST 50:31–43, 2016, doi:10.1016/j.ast.2015.12.016 (書誌のみ)。
- [PZ04] △ D. Papamoschou, A. Zill, "Fundamental Investigation of Supersonic Nozzle Flow Separation," AIAA 2004-1111, doi:10.2514/6.2004-1111 (書誌)。
- [PZJ09] △ D. Papamoschou, A. Zill, A. Johnson, "Supersonic flow separation in planar nozzles," Shock Waves 19(3):171–183, 2009, doi:10.1007/s00193-008-0160-z (書誌)。
- [XTP07] ✔ Q. Xiao, H. M. Tsai, D. Papamoschou, "Numerical Investigation of Supersonic Nozzle Flow Separation," AIAA J 45(3):532–541, 2007, doi:10.2514/1.20073.
- [BR05] △ A. Bourgoing, P. Reijasse, "Experimental analysis of unsteady separated flows in a supersonic planar nozzle," Shock Waves 14(4):251–258, 2005, doi:10.1007/s00193-005-0269-2 (書誌)。

### デュアルベル

- [GS11a] ✔ C. Génin, R. Stark, "Experimental investigation of the inflection geometry on dual bell nozzle flow behavior," AIAA 2011-5612, doi:10.2514/6.2011-5612.
- [GGS11] ✔ C. Génin, A. Gernoth, R. Stark, "Experimental and numerical study of heat flux in dual bell nozzles," AIAA 2011-5613, doi:10.2514/6.2011-5613 (雑誌版 JPP 29(1):21–26, 2013, doi:10.2514/1.B34479 △)。
- [SGWK12] ✔ R. Stark, C. Génin, B. Wagner, W. Koschel, "The altitude adaptive dual bell nozzle," 2012 (掲載会議は未確認)。
- [SG16] ✔ R. Stark, C. Génin, "Sea-level transitioning dual bell nozzles," Space Propulsion 2016 (SP2016_3124642)。雑誌版 CEAS Space J 9(3):279–287, 2017, doi:10.1007/s12567-017-0154-8 △。
- [HF94] ✔ M. Horn, S. Fisher, "Dual-bell altitude compensating nozzles," Penn State PERC Annual Report Vol. II, pp.140–147, 1993, NTRS 19940018584 (DLR 文献では NASA-CR-194719 として引用)。
- [J13] ✔ D. S. Jones, T. T. Bui, J. H. Ruf, "Proposed Flight Research of a Dual-Bell Rocket Nozzle Using the NASA F-15 Airplane," NASA/TM-2013-216590.
- [R16] ✔ NASA MSFC, "Altitude Compensating Nozzle," NTRS 20160008031.
- [NGS07] ✔ C. Nürnberger-Génin, R. Stark, "Experimental Study of Dual Bell Nozzles," 2nd EUCASS, 2007, DLR elib 63529.
- [NGS10] △ C. Nürnberger-Génin, R. Stark, "Experimental Study on Flow Transition in Dual Bell Nozzles," JPP 26(3):497–502, 2010, doi:10.2514/1.47282 (抄録)。
- [NGS09] △ C. Nürnberger-Génin, R. Stark, "Flow transition in dual bell nozzles," Shock Waves 19(3):265–270, 2009, doi:10.1007/s00193-008-0176-4 (書誌)。
- [GS11b] △ C. Génin, R. Stark, "Side Loads in Subscale Dual Bell Nozzles," JPP 27(4):828–837, 2011, doi:10.2514/1.B34170 (冒頭部のみ)。
- [VSNH10] △ S. B. Verma, R. Stark, C. Nürnberger-Génin, O. Haidn, "Cold-gas experiments to study the flow separation characteristics of a dual-bell nozzle during its transition modes," Shock Waves 20(3):191–203, 2010, doi:10.1007/s00193-010-0259-x (書誌)。
- [VSGH11] △ S. B. Verma, R. Stark, C. Génin, O. Haidn, "Cold gas dual-bell tests in high-altitude simulation chamber," Shock Waves 21(2):131–140, 2011, doi:10.1007/s00193-011-0302-6 (書誌)。
- [VSH12] △ S. B. Verma, R. Stark, O. Haidn, "Gas-Density Effects on Dual-Bell Transition Behavior," JPP 28(6):1315–1323, 2012, doi:10.2514/1.B34451 (書誌)。
- [VSH13] △ S. B. Verma, R. Stark, O. Haidn, "Reynolds Number Influence on Dual-Bell Transition Phenomena," JPP 29(3):602–609, 2013, doi:10.2514/1.B34734 (抄録)。
- [SchG16] △ D. Schneider, C. Génin, "Numerical Investigation of Flow Transition Behavior in Cold Flow Dual-Bell Rocket Nozzles," JPP 32(5):1212–1219, 2016, doi:10.2514/1.B36010 (抄録)。
- △ C. Génin, R. Stark, O. Haidn, K. Quering, M. Frey, "Experimental and numerical study of dual bell nozzle flow," Progress in Flight Physics 5:363–376, 2013, doi:10.1051/eucass/201305363 (書誌。PDF は取得を拒否された)。
- [HTHF02] △ G. Hagemann, M. Terhardt, D. Haeseler, M. Frey, "Experimental and Analytical Design Verification of the Dual-Bell Concept," JPP 18(1):116–122, 2002, doi:10.2514/2.5905 (抄録)。
- [GW04] △ A. Gross, C. Weiland, "Numerical Simulation of Separated Cold Gas Nozzle Flows," JPP 20(3):509–519, 2004, doi:10.2514/1.2714 (抄録)。
- [NOM05] △ F. Nasuti, M. Onofri, E. Martelli, "Role of Wall Shape on the Transition in Axisymmetric Dual-Bell Nozzles," JPP 21(2):243–250, 2005, doi:10.2514/1.6524 (抄録)。
- [Re11] △ P. Reijasse et al., "Wall pressure and thrust of a dual bell nozzle in a cold gas facility," Progress in Propulsion Physics 2:655–674, 2011, doi:10.1051/eucass/201102655 (書誌)。
- [To06] △ T. Tomita, M. Takahashi, M. Sasaki, H. Tamura, AIAA 2006-4375, doi:10.2514/6.2006-4375 (書誌)。
- [Ki09] △ T. Kimura et al., "Experimental and Analytical Study for Design of Dual-Bell Nozzles," AIAA 2009-5149, doi:10.2514/6.2009-5149 (書誌)。

### SERN

- [RL84] ✔ R. J. Re, L. D. Leavitt, "Static Internal Performance of Single-Expansion-Ramp Nozzles With Various Combinations of Internal Geometric Parameters," NASA TM-86270, 1984, NTRS 19850006488.
- [DA99] ✔ K. A. Deere, S. C. Asbury, "Experimental and Computational Investigation of a Translating-Throat, Single-Expansion-Ramp Nozzle," NASA/TP-1999-209138, NTRS 19990052644.
- [YXYM15] △ Y. Yu, J. Xu, K. Yu, J. Mo, "Unsteady transitions of separation patterns in single expansion ramp nozzle," Shock Waves 25(6):623–633, 2015, doi:10.1007/s00193-015-0595-y (書誌。実験の記述は [YT23] で確認)。
- [YT23] ✔ T. Yu, Y. Yu, Y. P. Mao, Y. L. Yang, S. L. Xu, "Comparative Study of OpenFOAM Solvers on Separation Pattern and Separation Pattern Transition in Overexpanded Single Expansion Ramp Nozzle," JAFM 16(11):2249–2262, 2023, doi:10.47176/jafm.16.11.1751.
- △ Y. Yu, "Over-expanded separation transitions of single expansion ramp nozzle in the accelerating and decelerating processes," AST 98:105674, 2020, doi:10.1016/j.ast.2019.105674 (抄録)。
- △ Y. Yu, Y. Mao, T. Yu, Y. Yang, S. Xu, "The Influence of External Flow Field on the Flow Separation of Overexpanded Single-Expansion Ramp Nozzle," Aerospace 10(11):958, 2023, doi:10.3390/aerospace10110958 (抄録)。
- [YXLS19] △ K. Yu, J. Xu, Z. Lv, G. Song, "Inverse design methodology on a single expansion ramp nozzle for scramjets," AST 92:9–19, 2019, doi:10.1016/j.ast.2019.05.054 (抄録)。
- [Tan09] △ J. Tan, J. Jin, G. Du, H. M. Ji, J. L. Xu, "Experimental and computational investigation of a over-expanded single-expansion-ramp-nozzle," 推進技術 (J. Propulsion Technology) 30(3):292–296, 2009 (北京航空航天大リポジトリの抄録)。
- [CRB92] △ F. J. Capone, R. J. Re, E. A. Bare, "Parametric Investigation of Single-Expansion-Ramp Nozzles at Mach Numbers From 0.60 to 1.20," NASA TP-3240, 1992, NTRS 19920024949 (書誌)。
- △ J. R. Carlson, K. S. Abdol-Hamid, "Prediction of static performance for single expansion ramp nozzles," AIAA 93-2571, doi:10.2514/6.1993-2571 (書誌)。
- [GHK00] △ P. Gruhn, A. Henckels, S. Kirschstein, "Flap contour optimization for highly integrated SERN nozzles," AST 4(8):555–565, 2000, doi:10.1016/S1270-9638(00)01080-4 (書誌)。

### 経験式

- [MG78] ✔ E. L. Morrisette, T. J. Goldberg, "Turbulent-Flow Separation Criteria for Overexpanded Supersonic Nozzles," NASA TP-1207, 1978, NTRS 19780022611.
- [Sch73] ✔ R. H. Schmucker, "Flow Processes in Overexpanded Chemical Rocket Nozzles. Part 1: Flow Separation," NASA TM-77396 (TUM-LRT-TB-7, 1973 の英訳, 1984), NTRS 19840011402.
- △ M. Arens, E. Spiegler, "Shock-induced boundary layer separation in over-expanded conical exhaust nozzles," AIAA J 1(3):578–581, 1963, doi:10.2514/3.1598 (書誌)。
- △ S. Kalt, D. Badal, "Conical rocket nozzle performance under flow-separated conditions," J. Spacecraft Rockets 2(3):447–449, 1965, doi:10.2514/3.28200 (書誌)。
- △ E. E. Zukoski, "Turbulent boundary-layer separation in front of a forward-facing step," AIAA J 5(10):1746–1753, 1967, doi:10.2514/3.4299 (書誌)。
- △ G. L. Romine, "Nozzle Flow Separation," AIAA J 36(9):1618–1625, 1998, doi:10.2514/2.588 (書誌)。
- △ M. Summerfield, C. Foster, W. Swan, "Flow separation in overexpanded supersonic exhaust nozzles," Jet Propulsion 24, 1954 ([O04] の文献表による)。
