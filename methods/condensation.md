# 非平衡凝縮 (4 モーメント方程式)

極超音速ノズルの強膨張で試験ガス (窒素 N$_2$ 等) が**非平衡凝縮**し、潜熱放出で静圧・静温が
上昇する現象を、気相 Euler/NS に**液相のモーメント輸送方程式**を結合して記述する。

主参照: L. Lin, W. Cheng, X. Luo, F. Qin, *On nitrogen condensation in hypersonic nozzle flows*,
Shock Waves 24:179–189 (2014)。物理係数・導出の詳細は
[`papers/.../_summary.md`](../papers/on%20nitrogen%20condensation%20in%20hypersonic%20nozzle%20flows_summary.md)
にまとめてある。実装対応は 本ドキュメントの「実装」節 を参照。

---

本ドキュメントは理論(係数・方程式)と実装(ソース対応)をまとめる。

## 理論

### 1. 支配方程式 — Euler/NS + method of moments

液滴雲はサイズ分布関数 $f(r;x,t)$ の有限個の**モーメント** $Q_n=\int_0^\infty r^n f\,dr$ で表す
(Hill 1966)。凝縮種ごとに次の保存ベクトルを気相に追加する:

$$
U = (\rho,\ \rho u,\ \rho v,\ \rho w,\ \rho E,\ \underbrace{\rho g,\ \rho Q_2,\ \rho Q_1,\ \rho Q_0}_{\text{凝縮種ごと 4 本}})^T
$$

- $\rho Q_0$ : 単位質量あたり液滴数密度
- $\rho Q_1,\ \rho Q_2$ : サイズ分布のモーメント
- $\rho g$ : 液相質量分率 $g = \tfrac{4}{3}\pi\rho_l Q_3$ (実質 $Q_3$ モーメント、$\rho_l$ は液密度)

各モーメント方程式は**移流項 + 凝縮ソース $S$ (核生成 + 成長)** を持つ:

$$
\frac{\partial(\rho Q_n)}{\partial t} + \nabla\!\cdot(\rho Q_n \mathbf u) = S_{Q_n},\qquad
\frac{\partial(\rho g)}{\partial t} + \nabla\!\cdot(\rho g\, \mathbf u) = S_g
$$

気相の質量・運動量・エネルギー式のソースは**ゼロ** (混合気体全体は保存)。気相と液相の結合は
**熱力学 (状態方程式・温度関係) 経由**で、潜熱が静温を上げる。モーメントは気相速度 $\mathbf u$ で
運ばれる受動スカラー (拡散なし) として移流される。

ソース項は核生成項 ($T$ 依存) と成長項 ($T_d$ 依存) の加算分離:

$$
S_{Q_n} = \underbrace{J\, r_*^{\,n}}_{\text{核生成}} + \underbrace{n\,\rho Q_{n-1}\frac{dr}{dt}}_{\text{成長}},\qquad
S_g = \frac{4}{3}\pi\rho_l\Big( J\, r_*^{\,3} + 3\rho Q_2\frac{dr}{dt}\Big)
$$

---

### 2. 核生成モデル (nucleation) — **気相温度 $T$ で評価**

古典核生成理論 (CNT, Becker–Döring, capillarity 近似):

$$
J_{CNT} = K_{CNT}\exp\!\left(-\frac{\Delta G^*}{k_B T}\right),\quad
K_{CNT}=\Big(\frac{2\sigma}{\pi M\upsilon_l}\Big)\Big(\frac{p_v}{k_B T}\Big)^2,\quad
\Delta G^*=\frac{4}{3}\pi r_*^2\sigma
$$

臨界半径 (Kelvin–Thomson):

$$
r_* = \frac{2\sigma}{\rho_l R T\,\ln(p_v/p_{sat})}
$$

過飽和 $S=p_v/p_{sat}(T)$。**核生成は気相温度 $T$ で評価**する (臨界核は sub-nm で自己加熱なしとみなす)。

#### 種ごとの補正 (切替可能)

| モデル enum | 補正 | 用途 |
| --- | --- | --- |
| `CNT` | 補正なし | 基準 |
| `CNT_Iland` | $J=J_{CNT}\exp(A+B/T)$, $A=-55,\ B=4270$ K | **N2** (Iland et al. 2009) |
| `CNT_Kantrowitz` | $J_{noniso}=J_{iso}/(1+b\,q^2)$ 系の非等温補正 | **H2O** (生成中クラスタ自己加熱) |

- **Iland 補正** (N2): 膨張チェンバー実験から決めた経験補正。本論文の圧力域 (2.8–6.0 kPa) で検証済み。
- **Kantrowitz 非等温補正** (H2O): 生成中の臨界核が放出潜熱を捨てきれず自己加熱し $J$ が下がる効果を
  **核生成項そのもの**に織り込む。N2 では使わず (自己加熱は成長側 $T_d$ で扱う、下記)。

  $$
  J_{noniso}=\frac{J_{iso}}{1+\theta},\qquad
  \theta=\frac{2(\gamma_v-1)}{\gamma_v+1}\,b\Big(b-\tfrac12\Big),\qquad b=\frac{L(T)}{R_vT}
  $$

  ここで $\gamma_v=c_{p,v}/c_{v,v}$ は**凝縮種 (蒸気) 自身の比熱比** (Kantrowitz 1951 / Feder et al. 1966 の純蒸気形:
  H2O 蒸気 1.331 [NASA-9, 200–300 K]、N2 1.4)。$\theta$ はクラスタが 1 分子を取り込むたびに受け取る潜熱
  $L$ を、蒸気分子との衝突で持ち去れるエネルギー揺らぎ ($\propto c_{v,v}+R_v/2$) で割った比なので、
  比熱比はキャリア (N2) や気相混合のものではなく蒸気のものを使う ($\frac{2(\gamma_v-1)}{\gamma_v+1}=\frac{R_v}{c_{v,v}+R_v/2}$)。
  **2026-09-10 以前の実装はセルの気相混合 $\gamma=c_{p,gas}/c_{v,gas}$ (H2O–N2 では N2 支配 ≈1.40) を使っていた**
  ($\theta$ が ~17 % 過大 → 純蒸気形の $J$ が ~15 % 過小; Wyslouzil 2D では onset が 0.45 mm 上流へ移り onset 帯の壁圧偏差が −4.8 → −4.5 % に縮む, 2026-09-10 検証)。
  これは**純蒸気形近似の中での係数修正**であり、
  旧挙動は `condKantrowitzGammaMode: 1` で A/B 用に残す。キャリア気体分子との衝突も揺らぎに数える Feder の
  carrier 拡張 ($b^2$ に $p_c/p_v\sqrt{m_v/m_c}$ 重みの項が加わり $\theta$ が大幅に小さくなる; Wedekind et al.) は未実装。
  $J_{pure}\le J_{carrier}\le J_{iso}$ は「同一の核生成障壁・前因子を固定し、衝突による熱除去だけを追加したモデル間」の関係で、
  carrier 中の核生成率の真値の保証範囲ではない (Wedekind et al. は carrier の $pV$ 仕事など逆向きの寄与も区別する)。
  純蒸気形は onset の実験一致を保証するものではない。
  計画: [plans/active/condensation-kantrowitz-gamma-twophase-sonic.md](../plans/active/condensation-kantrowitz-gamma-twophase-sonic.md)。

---

### 3. 成長モデル (droplet growth) — **液滴温度 $T_d$ で評価**

#### 修正 Gyarmathy / Goodheart (N2)

純窒素蒸気で広い Knudsen 数域に適用可能 (Goodheart 2004):

$$
\frac{dr}{dt} = \frac{1}{\rho_l}\,\frac{f_{F\text{-}S}(Kn)\,k\,R T^2}{L^2}\,\ln\!\left(\frac{p}{p_{sat}(T_d)}\right),\qquad
f_{F\text{-}S}(Kn)=\frac{1+2Kn}{r(1+3.42Kn+5.32Kn^2)}\Big(1-\frac{r_*}{r}\Big)
$$

$k$: 熱伝導率、$L$: 潜熱、$Kn=\lambda/2r$ (核生成の nm から成長後の µm までカバー)。
**成長は液滴温度 $T_d$ で評価**する (有限サイズ液滴は蓄積潜熱で $T_d>T_v$ になる)。

#### Hertz–Knudsen (H2O、切替候補)

質量フラックスベースの成長則。H2O 検証 (Wyslouzil) 用に同じ枠組みで切替可能とする。

#### 平衡凝縮 (`condEquilibrium=1`, 2026-08-18)

非平衡モデルの上限 (熱力学的限界) として、核生成・成長を経ず**各セルで局所瞬時平衡**を仮定するモデル。
セル状態 $(\rho, e, Y_w, g)$ で、液相分率を $g\to g+\Delta$ だけ変えたときの温度は二相 EOS
($e=e_{gas}(Y,T)+g(R_wT-L)$、$e$ 一定) の線形化から

$$
T(\Delta)=T+\Delta\,\frac{L-R_wT}{c_{v,g}+gR_w},\qquad
p_v(\Delta)=\rho\,(Y_w-g-\Delta)\,R_w\,T(\Delta)
$$

であり ($Y_w$ は蒸気+液の総水分率、pure-condensible では $Y_w\to1$)、平衡条件
$F(\Delta)\equiv p_v(\Delta)-p_{sat}\!\big(T(\Delta)\big)=0$ を $\Delta\in[-g,\;Y_w-g]$ で二分法+Newton で解いて
$g_{eq}=g+\Delta$ を得る。$F(0)\le0$ かつ $g=0$ (未飽和・液なし) なら $\Delta=0$、$F(-g)<0$ (全蒸発しても未飽和) なら $\Delta=-g$。
ソースは緩和形 $S_g=\alpha\,\rho\,\Delta/\Delta t_{loc}$ (`condEqRelax` $\alpha$、非平衡と同じ θ 律速) で、
陰的ヤコビアンは $\partial S_g/\partial(\rho g)=-\alpha\theta/\Delta t$ (負帰還)。定常では $\Delta\to0$ で線形化誤差も消え、
場は厳密に $S=1$ (飽和線上) か $g=0$ (未飽和) のどちらかになる。蒸発は $\Delta<0$ として同式に含まれる。
モーメント $Q_0$–$Q_2$ は輸送するだけ (ソース 0、液滴径の情報なし)。飽和線は非平衡と同じ **過冷却液** $p_{sat}$
(Murphy–Koop) を使う (氷ではない — CEA 平衡の H₂O(cr) とは基準が違うので比較時に注意)。

**用途**: 非平衡 (Wilson 点で遅れて一気に凝縮) に対し「凝縮し得る最大量と、それによる $T$・$p$・$M$ の変化」の上限評価。
(θ 遅れの無い代数拘束版は次節の `condEquilibrium=2`。)
実装は `condensation_source_d` の `eq` 分岐 (`cond_equilibrium_delta`)。診断 `condTsat_<s>` = $T_{sat}(p_v)$
(過冷却度 $T_{sat}-T$ の評価用、非平衡でも出力) を併せて追加した ($p_{sat}$ の Newton 反転、$p_v\to0$ では 0)。

#### 平衡凝縮 — EOS 拘束形 (`condEquilibrium=2`, 2026-08-19)

上の緩和形 (`condEquilibrium=1`) は「各セルで過冷却ゼロになる $g_{eq}$」をソース項で擬似時間緩和させるため、
定常では飽和線に乗るが、θ 律速 (`condEqDTmax`/`condEqDgMax`) の分だけ **onset 直後 ~2 $r_t$ で $S$ が 1.2 程度残る**
(case/44 `run_0097`)。EOS 拘束形は $g$ を輸送量ではなく**状態量**として扱い、従属変数計算 (`dependentVariables_d`) の中で
$(T, g)$ を保存量 $(\rho, e, Y_w)$ から代数的に決める (湿り度 β の内部反復に相当):

$$
\text{未飽和: } g=0,\ e_{gas}(Y,T)=e_{in},\ p_v=\rho Y_wR_wT\le p_{sat}(T)
\qquad
\text{飽和: } \begin{cases} e_{gas}(Y,T)+g\,(R_wT-L(T))=e_{in}\\ \rho\,(Y_w-g)\,R_w\,T=p_{sat}(T)\end{cases}\quad 0<g\le Y_w
$$

解法: まず $g=0$ で $T_0$ を反転し $p_v(T_0)\le p_{sat}(T_0)$ なら乾き。飽和なら $G(g)\equiv\rho(Y_w-g)R_wT(g)-p_{sat}(T(g))$
($T(g)$ は内側 Newton `cond_T_from_e_carrier`) を $g\in[0,Y_w]$ で解く。$G(0)>0$, $G(Y_w)=-p_{sat}<0$ で $G$ は単調減少
($g$↑ で蒸気↓・潜熱で $T$↑・$p_{sat}$↑) なので**唯一解**。前ステップの $g$ を初期値に括弧付き Newton
($dG/dg=-\rho R_wT+[\rho(Y_w-g)R_w-p_{sat}']\,dT/dg$, $dT/dg=(L-R_wT)/(c_v+g(R_w-L'))$) を回し、括弧外に出れば二分法へ退避する。
pure-condensible では $p_v=\rho(1-g)R_vT$、$Y_w\to1$。

得られた $g_{eq}$ を **`rog` に書き戻す (射影)** ので、SLAU 面温度・潜熱補正・出力 `g_<s>` は全て $g_{eq}$ を見る。
`rog` の輸送方程式は残差 0・ソース 0 にして凍結 (`rms_rog` は恒等 0)、モーメント $Q_0$–$Q_2$ は緩和形と同じく輸送のみ (意味なし)。
診断 `condS_<s>`/`condTsat_<s>` は緩和形と同じく source kernel が書く (飽和セルでは厳密に $S=1$、$T_{sat}=T$)。

緩和形との関係: 定常解の固定点は同一 (どちらも「飽和線上か乾き」)。差は (i) onset 直後の θ 遅れが無い、
(ii) 過渡・非定常でも各ステップで厳密平衡、(iii) $g$ が対流でなく熱力学で決まるので流線に沿った履歴を持たない
(平衡なので本来そうあるべき)。飽和線は緩和形と同じ過冷却液 (Murphy–Koop)。陰解法との結合は従来どおり loose
(NS ブロックの音速と $\kappa$ は全蒸気 frozen のまま。`condSonicModel` の自動解決は `condEquilibrium 2` を対象外にする (§5「二相 frozen 音速」、2026-09-10)。平衡音速への置換は未実装)。
`nCondSpecies>=2` では未対応 (エラー)。実装は `cond_eq_solve_g` / `cond_equilibrium_Tg_{carrier,pure_tp,pure_cpg}` (`condensationEOS_d.cuh`)、
単体テスト `tests/unit/test_cond_equilibrium_eos.cpp`。case/44 `run_0098` (旧条件×va3 形状) で onset 以降 S=1.0000・過冷却 0.000 K、下流は緩和形と一致。設定・検証は
[plans/accepted/condensation-equilibrium-eos.md](../plans/accepted/condensation-equilibrium-eos.md)。

#### 統一駆動力と蒸発 (負成長, `condEvaporation=1`)

上の成長則はいずれも「周囲蒸気圧 $p_v$ と、半径 $r$ の液滴が平衡になる Kelvin 蒸気圧
$p_d(r)=p_{sat}(T)\exp(K_e/r)$, $K_e\equiv 2\sigma/(\rho_l R T)$ の差」を駆動力にしている。
臨界半径は $p_v=p_d(r_*)$ の解 $r_*=K_e/\ln S$ なので

$$
\ln S\left(1-\frac{r_*}{r}\right)=\ln S-\frac{K_e}{r}=\ln\frac{p_v}{p_d(r)},\qquad
p_v-p_d(r)\ (\text{Hertz–Knudsen})
$$

と、Goodheart / Gyarmathy の $(1-r_*/r)\ln S$ も Hertz–Knudsen の $p_v-p_d$ も**同じ量**である。
過飽和 $S\le1$ では $r_*$ が定義されない ($\ln S\le0$) が右辺は定義され常に負、すなわち**蒸発**
(液滴は縮み、潜熱を吸収して蒸気に戻る)。forge の初期実装は $dr/dt<0$ を 0 にクランプし
(凝縮のみ)、過熱域 (圧縮波・衝撃波・高温境界層) で液相が受動スカラーとして凍結していた。
`condEvaporation=1` (既定) で $S\le1$ のセルに次の蒸発モデルを適用する
([計画](../plans/accepted/condensation-evaporation.md)):

- **駆動力**: 既定 (`condEvapKelvin=0`) は Kelvin 項を落とし平面飽和圧で評価する — HK は
  $\frac{\alpha}{\rho_l}\frac{p_v-p_{sat}}{\sqrt{2\pi RT}}$、Goodheart/Gyarmathy は前因子 $\times\ln S\times$ Kn 補正。
  Kelvin 項 ($e^{K_e/r}$, 水 290 K で $K_e\approx1.1$ nm: 100 nm で 1 %, 10 nm で 12 %) は
  「小さいほど速く蒸発する」正帰還で離散的に $Q_1\to0$ 崩壊 ($1/r$, $e^{K_e/r}$ 発散) を招くうえ、
  **質量収支には効かない** ($r_0$ から $r$ まで縮んだ残量は $(r/r_0)^3$; $r_0$≈300 nm なら
  $r<10$ nm 段階の残量は $4\times10^{-5}$)。`condEvapKelvin=1` で Kelvin 込みにできる。
- **評価半径**: 体積平均半径 $r_{30}=\big(g/(\tfrac43\pi\rho_lQ_0)\big)^{1/3}$ を使う
  ($Q_1$ に敏感な $r_{10}=Q_1/Q_0$ でなく、θ 律速される $g$ を通す)。
- **λ スケール更新** (monodisperse 整合): 1 step の縮小比 $\lambda=r_{new}/r_{30}=1+\dot r\,\Delta t/r_{30}$
  で $Q_1\to\lambda Q_1,\ Q_2\to\lambda^2Q_2,\ g\to\lambda^3 g$ ($Q_0$ 不変)、ソースは
  $S_\phi=(\phi_{new}-\phi)/\Delta t$。$Q_1^2\le Q_0Q_2$ 等の実現可能性が自動で保たれる。
  $\lambda\ge\lambda_{\min}=\tfrac12$ (半径半減/step)、$g(1-\lambda^3)\le$ `dg_max`、
  潜熱冷却 $g(1-\lambda^3)L/c_v\le$ `dT_max` (1 K/step) で下から律速。
- **液滴消滅**: $\lambda_{\min}r_{30}<r_{\min}$ すなわち $r_{30}<2r_{\min}$ (`condEvapRmin`=1 nm) なら
  中間状態を作らず $\lambda=0$ (全モーメント消滅、$Q_0$ も 0)。実現可能性クランプ kernel でも
  「$S\le1$ かつ ($Q_0=0$ または $r_{30}<2r_{\min}$) かつ $g\le5\times10^{-7}$」で 4 モーメントを硬 0 化する
  (陰的緩和では厳密 0 に届かないため)。$g=0$ で $Q_n$ だけ残る float アンダーフローの「モーメント塵」も
  $S\le1$ では同時に掃除する。$Y_w$ (総水) は不変なので質量は自動保存、エネルギーは
  二相 EOS $e=e_v-gL$ が $g\to0$ で潜熱吸収を自動処理する。
- **点陰的ヤコビアン**: $g$ 行に蒸気復帰 $\partial S_g/\partial(\rho g)$ (carrier では
  $\partial p_v/\partial(\rho g)=-R_wT$) と潜熱冷却 $\partial S_g/\partial T\cdot\partial T/\partial(\rho g)$ を
  数値微分で入れる (どちらも負帰還 → `sj_g`$\ge0$)。$Q_1$ 行の $\partial S/\partial Q_1$ は
  Kelvin 正帰還側なので入れない。
- **適用範囲**: $S\le1$ のみ。$S>1$ の亜臨界 ($\bar r<r_*$) 蒸発は核生成域で $\bar r$ が $r_*$ に
  張り付く既知の罠 (`COND_RNUC_FAC`) を再燃させるため従来どおり成長 0 (クランプ) のまま。
  核生成 $J$ は $S\le1$ で元々 0。二温度 (`condTwoTemp`) の $T_d\le T$ 側は未対応 (一温度のみ)。

蒸発は負帰還 (蒸発→冷却→$p_{sat}$↓・蒸気復帰→$S\to1$) なので、圧縮帯では**平衡蒸発**
($S\to1$ に漸近、$g$ は僅かに減) になり、高温境界層 (回復温度 $\gg T_{sat}$) では液相が消える。

---

### 4. 液滴温度 $T_d$ — Hill のエネルギーバランス (純蒸気簡約)

$T_d$ は液滴表面のエネルギーバランス (Hill モデル, Johnson 2001 の比較で採用) から反復で解く。
純窒素では不活性キャリアガス寄与 $\beta_c$ 項を落とし ($\zeta=1.0$):

$$
K_v\!\left(1-\frac{\beta_d T_d}{\beta T_v}\right) + \left(1-\frac{\beta_d}{\beta}\right)(L_v-1)\frac{\gamma}{\gamma-1} = 0
$$

$$
K_v=\frac{1}{2}\frac{\gamma+1}{\gamma-1},\quad
L_v=\frac{(\gamma-1)L}{\gamma R T_v},\quad
\beta=\frac{p_v}{\sqrt{2\pi R T_v}},\quad
\beta_d=\frac{p_d}{\sqrt{2\pi R T_d}}
$$

液滴表面の蒸気圧は **Kelvin 効果込み**:

$$
p_d = p_\infty(T_d)\,\exp\!\left\{\frac{2\sigma}{\rho_l R T_d\, r}\right\}
$$

$\beta_d$ が $p_d/\sqrt{T_d}$ を通じて $T_d$ に依存し、$p_d$ も $T_d,r$ に依存するため、$T_v,p_v,r$ を所与に
$T_d$ を反復で求める。$T_d$ は成長則 (式 3) と Kelvin 式にのみ入り、**核生成には直接入らない** (核生成は
バルク $T$ 上昇を介した間接ループでのみ抑制される)。

#### 実装した準定常 $T_d$ 閉包 (carrier + Hertz–Knudsen, `condTwoTemp=1`)

上の Hill 式は純蒸気 (N2) 向けの整理である。希薄水/N2 (carrier) では、**準定常の液滴エネルギーバランス**
「凝縮による潜熱解放 = 液滴から気相 (キャリア) への熱伝導」を直接解く形で実装した。Hertz–Knudsen
(質量律速) 成長に対し、質量フラックス $j$ と熱フラックスを連立する:

$$
\underbrace{L\,j(T_d)}_{\text{潜熱解放}} = \underbrace{h\,(T_d-T_g)}_{\text{気相への熱伝導}},\qquad
j(T_d)=\frac{\alpha\,[\,p_v-p_d(T_d)\,]}{\sqrt{2\pi R T_g}},\qquad
h=\frac{\lambda_g}{r\,(1+3.18\,Kn)}
$$

液滴表面の平衡蒸気圧は **液滴温度 $T_d$ で** Kelvin 補正: $p_d(T_d)=p_{sat}(T_d)\exp\!\{2\sigma/(\rho_l R T_d r)\}$。
$\lambda_g$ はキャリア (N2) の熱伝導率、$Kn=\lambda/2r$ (平均自由行程は全圧基準)。この 1 変数非線形方程式
$F(T_d)=L\,j(T_d)-h(T_d-T_g)=0$ を **Newton 法** ($\partial p_{sat}/\partial T_d$ と Kelvin 項の微分を含む解析勾配,
$T_d\ge T_g$ にクランプ) で解き、得た $T_d$ で成長率 $\dot r=j(T_d)/\rho_l$ を評価する。

- **効果**: 自己加熱で $T_d>T_g$ となり $p_d(T_d)>p_d(T_g)$、実効過飽和 $p_v-p_d$ が下がるため**成長が遅くなる**
  (一温度 $T_d=T_g$ は過冷却を過大評価し成長を過大評価していた)。$Kn$ が大きい (小液滴) ほど $h$ が小さく
  $T_d$ の上振れが大きい。
- **適用範囲**: 本フラグは **Hertz–Knudsen 経路 (`condGrowthModel=0`) のみ**に効かせる。Gyarmathy
  (`condGrowthModel=1`, 熱伝導律速) は元来「液滴が $T_s$ 付近まで自己加熱した極限」を仮定し過冷却
  $(T_s-T_g)$ を駆動力に使う式なので、$T_d$ を二重計上しないため非適用 (一温度のまま)。
- **簡約**: $T_d$ は**成長キネティクスのみ**に使う。二相 EOS のバルクエネルギーは液滴を $T_g$ のままとする
  (液相顕熱差 $\sim g\,c_l\,\Delta T_d \sim 1\%$ オーダーで潜熱 $\sim$ 数% に対し小、既知の近似)。
- **既定 off** (`condTwoTemp=0`, 一温度 $T_d=T_g$) でビット不変。

---

### 5. 二相 EOS と圧力 — 気相分圧近似

窒素を熱量的完全気体として扱う (実在気体効果 $<4\%$)。**液滴体積を無視**し、蒸気を理想気体と
みなすと、圧力に効くのは気相分圧だけ:

$$
\boxed{\ p = (\rho-\rho g)\,R\,T\ } \qquad (R: \text{蒸気の気体定数})
$$

混合の内部エネルギー (潜熱込み) から $T$ を逆算:

$$
\rho e = (\rho-\rho g)\,e_v(T) + \rho g\,e_l(T),\qquad \rho e=\rho E-\tfrac12\rho|\mathbf u|^2
$$

$e_v-e_l = L - RT$ (潜熱)。混合体積比熱 $C_v=(\rho-\rho g)c_{v,v}+\rho g\,c_l$。

#### 圧力微分 (一般 EOS Roe の $(\chi,\kappa)$ + 凝縮新項 $\xi_g$)

$$
\kappa\equiv\left.\frac{\partial p}{\partial(\rho e)}\right|_{\rho,\rho g}=\frac{(\rho-\rho g)R}{C_v}
\quad(\gamma-1\ \text{に相当}),\qquad
\chi\equiv\left.\frac{\partial p}{\partial\rho}\right|_{\rho e,\rho g}=RT-\kappa\,e_v
$$

$$
\boxed{\ \xi_g\equiv\left.\frac{\partial p}{\partial(\rho g)}\right|_{\rho,\rho e}=-RT+\kappa(e_v-e_l)=-RT+\kappa(L-RT)\ }
$$

$\xi_g$ は二相化で**新規に出る圧力微分**で、「$\rho g$↑ = 蒸気→液化 = 潜熱放出 → $T$↑ → $p$↑」を表す
(N2 70 K 概算で $\xi_g\approx+5\times10^4>0$、凝縮で昇圧)。保存量での $\partial p/\partial U$ と flux
Jacobian への入り方、Roe/SLAU への影響は 本ドキュメントの「実装」節 の対流フラックス節を参照。

#### 保存量での $\partial p/\partial U$ と semi-perfect (TP) 一般化

内部エネルギー密度 $\tilde\epsilon=\rho e=\rho E-\tfrac12|\mathbf m|^2/\rho$ ($\mathbf m=\rho\mathbf u$) の
全微分 $dp=\chi\,d\rho+\kappa\,d\tilde\epsilon$ と連鎖則
($\partial\tilde\epsilon/\partial\rho|_{\mathbf m,\rho E}=e_k\equiv\tfrac12|\mathbf u|^2$,
$\partial\tilde\epsilon/\partial m_k=-u_k$, $\partial\tilde\epsilon/\partial(\rho E)=1$) より、保存量
$q=(\rho,\rho u,\rho v,\rho w,\rho E)$ に対して

$$
\frac{\partial p}{\partial q_1}=\chi+\kappa\,e_k,\quad
\frac{\partial p}{\partial q_{2,3,4}}=-\kappa(u,v,w),\quad
\frac{\partial p}{\partial q_5}=\kappa
$$

凝縮種は $\partial p/\partial(\rho g)=\xi_g$ の列が加わる。

**semi-perfect (thermally-perfect) gas**: 比熱が $T$ 依存で $\epsilon(T)=\int_0^T C_v\,dT'\ne C_v(T)\,T$
のとき、$T,\rho,\tilde\epsilon$ の全微分
($(\partial T/\partial\rho)|_{\tilde\epsilon}=-\epsilon/(\rho C_v)$,
$(\partial T/\partial\tilde\epsilon)|_\rho=1/(\rho C_v)$) を用いて

$$
\chi=R\Big(T-\frac{\epsilon(T)}{C_v(T)}\Big)=\frac{p}{\rho}-(\gamma(T)-1)\epsilon(T),\qquad
\kappa=\gamma(T)-1
$$

となり、**$\chi\ne0$** となる (calorically-perfect は $\epsilon=C_v T$ の特殊形で $\chi=0$ に縮約)。
本書の $\chi=RT-\kappa e_v$ は $e_v$ を**積分内部エネルギー** $\epsilon_v(T)$ と取れば
($\kappa e_v=(R/C_v)\epsilon_v=(\gamma-1)\epsilon_v$) この一般形と一致する。forge は既に多成分 TP の陰解法
Jacobian で per-cell $\kappa=\gamma[ic]-1$ を実装済み ([thermophysics plan](../plans/accepted/thermophysics-multicomponent-tpgas.md))
であり、凝縮 Phase 2 は同枠で $\kappa\to(\rho-\rho g)R/C_v$ (混合) と $(\rho g)$ 列の $\xi_g$ を足す最小拡張となる。

---

#### 対流流束の面温度と二相圧力の整合 (2026-08-18 修正)

TP (`thermalMethod 2`) の SLAU 流束は面エンタルピー $h=h_{mix}(T_f)+\tfrac12|u|^2-gL$ を、MUSCL 面値 $(P_f,\rho_f)$ から
$T_f=P_f/(\rho_f R)$ で再構成する。ここで $R$ に**全水分を気相と数えた** $R_{mix}$ を使うと、状態側の圧力が蒸気のみ
$P=\rho T(R_{mix}-gR_w)$ (§5) なので $T_f$ が $gR_w/R_{mix}$ (~2 %) 低く出て、面エンタルピーが状態より $c_p\Delta T\approx6$ kJ/kg
小さくなる。定常解はこの流束で $h_0$ が一定になるよう調整されるため、**凝縮帯で全エンタルピー $(\rho e+P)/\rho$ が +0.6 %
(潜熱の ~16 %) 増え、$T$ が 4–5 K 高い**非保存解になっていた (case/44 va2 M4.75: 軸 $h_0$ 836→841 kJ/kg。出口ノードだけ
境界流束が `Ht` を直読して整合し、$T$ が 5 K 低く $g$ が跳ねる「出口の跳ね」として顕在化)。修正: 面温度の再構成に
$R_{gas}=R_{mix}-gR_w$ (pure-condensible は $(1-g)R$) を使う (`convectiveFlux_slau_d.inc.cuh`)。修正後は軸 $h_0$ が
836.2±0.1 kJ/kg で保存、出口ノードの跳ねは消える。同型の不整合は `outlet_statPress` の TP ghost ($T=P/(\rho R_{mix})$) にも
あり、超音速流出では `roe/T/sonic` を内部値で全量外挿するよう変更 (亜音速流出 × 凝縮は近似のまま、要フォロー)。
ROE 流束の TP 分岐は単成分のみで二相補正未実装 (凝縮 TP は SLAU を使うこと)。

### 6. 数値解法 — fractional-step

支配方程式を「ソース無しの均質部 (対流)」と「凝縮ソース部」に**演算子分割** (fractional-step) する。
均質部は forge 既存の有限体積対流 (SLAU/Roe) + 時間積分で解き、凝縮ソース ($J,\ dr/dt$) は別途
積分する。核生成率 $J\propto\exp(-\Delta G^*/k_BT)$、$dr/dt\propto\ln(p/p_{sat})$ は過飽和に指数的・
対数的で極めて stiff なため、ソース積分は **point-implicit** (源項ヤコビアン $\partial S/\partial U$ を
対角に組む) とする。

#### 精度・無次元化

数密度 $\rho Q_0\sim10^{18}$–$10^{22}$ /m$^3$、$r\sim$ nm のため $Q_0$ と $Q_3$ は $\sim$27 桁跨ぐ。
モーメントを基準量で無次元化して O(1) に抑える:

$$
\mu_n = \frac{Q_n}{N_{ref}\,r_{ref}^{\,n}}\quad(\text{例 } N_{ref}=10^{18}\,/\text{kg},\ r_{ref}=1\,\text{nm})
$$

平均半径 $r=(\mu_3/\mu_0)^{1/3} r_{ref}$ 等のモーメント間比も無次元量どうしで安定する。**まず float
で進め**、無次元化で桁を抑える。float で桁落ちが顕在化したらモーメントの double 化 (混合精度) を
フォールバックとして検討する (詳細は implementation.md)。

---

### 7. 検証

1. **dry 一致**: 凝縮ソースを切った状態 ($g\equiv0$) で case/34 Arthur ノズル run_0003 と一致 (回帰防止)。
2. **N2 凝縮 ON**: case/34 で貯気を振り、中心線静圧が dry 等エントロピー線より**上振れ**することを示す
   (論文 Fig.11)。
3. **H2O 凝縮**: case/16.nozzle_wys + `papers/condensation/wyslouzil2000.pdf` で水蒸気凝縮を検証し、
   N2 と H2O を同じ枠組みで切替できることを実証する。

---

### 8. N2 物性フィット (Lin 2014 Appendix 1)

実装で使う窒素の物性相関。臨界 $T_c=126.192$ K (表面張力は $T_c=126$ K)、臨界密度
$\rho_c=313.3$ kg/m³、三重点 $T_t=63.15$ K、$R_{N_2}=296.8$ J/(kg·K)。

#### 相の方針 — 過冷却液 (supercooled liquid)、固体ではない

**Lin 2014 のモデルは液相 (過冷却液)**。非平衡凝縮は三重点以下でも準安定な**過冷却液滴**を作る
(均質核生成は液滴を生成、結晶化は膨張時間に対して遅い) ため、本実装も**全域で液相フィットを使う**
(固相昇華へは切替えない)。論文本文も "conservation of the **liquid** phase", "**liquid** mass
fraction $g$", "condensation of **liquid droplets**" とし、$r_*=2\sigma/(\rho_l RT\ln(p_v/p_{sat}))$
は**液密度 $\rho_l$** を使う。Appendix 1 が固相フィット ($p_{sat}^s,\rho_s,L_s,\sigma_s$) も併載するのは
図 (Fig.9,10) の比較用であり、**シミュレーションでは使わない**。

> **低温外挿の注意 (実装ガード)**: 液相フィット (Jacobsen $p_{sat}$、潜熱多項式) は **~40K 以下に外挿
> すると破綻**する ($p_{sat}$ が ~33K で非単調、潜熱が崩落)。有効域は概ね **45–126K**。凝縮 ON では
> 潜熱放出で活性域が ~40K 以上に留まる (論文 Ma6 ノズルで出口 40→56K) が、過渡/dry セルが低温に落ちても
> 破綻しないよう、**物性評価温度を $[45\text{K}, T_c)$ にクランプ**する (より低温=より低 $p_{sat}$ の
> 安全側=過剰核生成を起こさない)。既知の近似。実装 `COND_T_PROP_FLOOR`
> ([condensationProperties_d.cuh](../solver_density_cuda/cuda_forge/condensationProperties_d.cuh))。

> **気相 thermo の制約**: 凝縮ケースの気相は **`thermalMethod 0` (熱量的完全気体, $c_p/\gamma$ 一定)**
> を使うこと。NASA-9/CEA は ~200K 未満で無効 (出口 ~27K) で、低温気相物性は CEA に無い。論文も
> calorically perfect gas。

以下のフィット式は参照用に液・固双方を併記するが、**既定実装は液 (過冷却) のみ**を使う。

#### 飽和蒸気圧 $p_{sat}(T)$

- **液 (Jacobsen, 式22)** [atm], $T_c=126.192$:
$$
\ln p_{sat}^l = \frac{n_1}{T}+n_2+n_3T+n_4(T_c-T)^{1.95}+n_5T^3+n_6T^4+n_7T^5+n_8T^6+n_9\ln T
$$
$n_1=8394.409444,\ n_2=-1890.045259,\ n_3=-7.282229165,\ n_4=0.01022850966,$
$n_5=5.556063825\times10^{-4},\ n_6=-5.944544662\times10^{-6},\ n_7=2.715433932\times10^{-8},$
$n_8=-4.879535904\times10^{-11},\ n_9=509.5360824$ (→ ×101325 で Pa)。
- **固 (Frels, 式23)** [mmHg]: $\log_{10} p_{sat}^s = 7.614676 - 356.281/T$ (→ ×133.322 で Pa)。

#### 凝縮相密度 $\rho_l/\rho_s(T)$

- **液 (Nowak, 式24)**, $\tau=1-T/T_c$, $\rho_c=313.3$:
$$
\ln(\rho_l/\rho_c)=n_1\tau^{0.3294}+n_2\tau^{4/6}+n_3\tau^{16/6}+n_4\tau^{35/6}
$$
$n_1=1.48654237,\ n_2=-0.280476066,\ n_3=0.0894143085,\ n_4=-0.119879866$。
- **固 (Scott, 式25)** [kg/m³]: $\rho_s = 1068.49 - 1.97830\,T$。

#### 潜熱 $L(T)$

**H₂O (2026-08-18 更新)**: $L(T)=h_v(T)-h_l(T)$ を **CEA (NASA-9) の気相 H₂O と液相 H₂O(L) の全エンタルピー差**として作る
(`h2o_latent`, `condensationProperties_d.cuh`)。$h_v$ は forge 種 DB と同じ CEA 2002 の H₂O 200–1000 K 係数、$h_l$ は CEA
`thermo.inp` の H₂O(L) 273.15–373.15 K 係数 (Cox 1989 / Haar 1984)。両方とも**生成エンタルピー込みの絶対値**で評価するので
差は datum に依らず、forge の sensible-enthalpy シフト (`thermoHrefTemp`) の影響を受けない (CEA 内での気液差をそのまま保つ)。
**273.15 K 未満は CEA に液相フィットが無い** (CEA は氷 H₂O(cr) に切替) — H₂O(L) 多項式の外挿は 250 K 以下で発散する
($c_{p,l}$ 230 K で 10 kJ/kgK) ため、$h_l$ を 273.15 K の値と勾配 $c_{p,l}(273.15)=4228$ J/kgK で線形外挿する (過冷却水の標準的扱い)。
値: $L$(273.15)=2.501, $L$(250)=2.556, $L$(230)=2.603, $L$(200)=2.674 MJ/kg。旧線形フィット $3.1485\times10^6-2370T$ は
この構成と <0.05 % で一致していた (等価; case/44 M4.75 で T 差 0.002 K, g 差 1e-4)。氷 (昇華熱 2.835 MJ/kg) は使わず、
飽和線も過冷却液 (Murphy–Koop) で統一している。

N2 は従来どおり Lin 2014 のフィット (`n2_latent`)。

#### 表面張力 $\sigma(T)$

- **液 (Stansfield, 式28)** [dyn/cm], $T_c=126$, $\sigma_0=29.06$: $\sigma_l=\sigma_0(1-T/T_c)^{1.247}$ (→ ×$10^{-3}$ で N/m)。
- **固 (Dotson, 式29)** [N/m]: $\sigma_s = 2.0\times10^{-4}(\rho_s)^{2/3}$ ($\rho_s$ は式25)。

#### Kelvin 効果と $T_d$ (Appendix 2)

液滴表面蒸気圧 (式32): $p_d=p_\infty(T_d)\exp\{2\sigma/(\rho_l RT_d r)\}$。$T_d$ の Hill 陰的式は
[4 節](#4-液滴温度-t_d--hill-のエネルギーバランス純蒸気簡約) 参照。$R_{N_2}=296.8$ J/(kg·K)。

## 実装

### 1. 段階実装 (Phase 分け)

| Phase | 内容 | 状態 |
| --- | --- | --- |
| **Phase 1** | 4 モーメントを**受動スカラー** (ソース=0) として輸送する骨格。case/34 dry 回帰一致 | 本セッション |
| **Phase 2** | N2 凝縮物理 (CNT+Iland 核生成 / Goodheart 成長 / Hill $T_d$ / 二相 EOS / point-implicit ソース / $\mu_n$ 無次元化) | 次セッション |
| **Phase 3** | H2O (CNT+Kantrowitz / Hertz–Knudsen) + case/16 Wyslouzil 検証 | 後続 |

Phase 1 では `solverConfig` の `condensation`/`nCondSpecies` で機構を ON/OFF し、モーメントは
気相速度で運ばれる受動スカラー (拡散なし、ソース=0) として移流するだけ。ソース=0 ゆえ
$g\equiv0$ で気相に影響せず、dry 計算と場が一致することを回帰確認する。

---

### 2. 追加保存変数 (多成分一般化)

RANS の固定 2 本 (`roK`,`roOmega`) と違い、モーメントは**凝縮種ごと 4 本 × 可変種数**。化学種輸送の
動的登録 (`variables::registerSpecies`, [variables.cpp](../solver_density_cuda/variables.cpp)) を手本に
`registerCondensation(nCondSpecies)` を新設し、種 $s$ ごとに次を `cellValNames`/`c`/`c_d` へ追加する
([variables.hpp](../solver_density_cuda/variables.hpp) の roK/roOmega 構成を 1:1 踏襲):

- 保存量: `rog_<s>`, `roQ2_<s>`, `roQ1_<s>`, `roQ0_<s>` (Phase 2 で無次元 `roμn_<s>` に読み替え)
- RK ステージ: `*N_<s>`, `*M_<s>`
- 残差: `res_rog_<s>`, `res_rog_<s>_m` ほか
- point-implicit 対角: `src_jac_g_<s>`, `transport_diag_g_<s>` ほか
- primitive: `g_<s>`, `Q0_<s>` ほか (`dependentVariables` で $\rho\phi/\rho$ から復元)

`allocVariables` の前に 1 度だけ呼ぶ (`registerSpecies` の隣)。`output_cellValNames` に
`rog_<s>`/`g_<s>` 等を追加して可視化対象にする。`nCondSpecies<=0` のときは何も登録せず既存経路を保つ。

---

### 3. 移流 — 既存 ScalarTransportDesc の流用

汎用スカラ輸送コア
[`scalarTransport_d.{cu,cuh}`](../solver_density_cuda/cuda_forge/scalarTransport_d.cuh) の
`ScalarTransportDesc` (1 次風上移流 + point-implicit `transport_diag` 対角 + 任意拡散) は物理非依存。
RANS の `ransTransport_d.cu` (`buildScalarDescs`/`ransTransport_d_wrapper`,
[ransTransport_d.cu:83-115](../solver_density_cuda/cuda_forge/ransTransport_d.cu)) と同型の
**`condensationTransport_d.{cu,cuh}`** を新設する:

- 種 × 4 モーメントぶんの `ScalarTransportDesc` を構築し `scalarTransportResidual_d()` を回す。
- モーメントは分子拡散しない受動スカラーなので **`diffusion=0`** (SST の $\sigma_k/\sigma_\omega$ 拡散は不要)、
  `floor=0`。Phase 1 は `src_jac=0`。
- ゲート: `condensationEnabled(cfg) = (cfg.condensation==1 && cfg.nCondSpecies>=1)`。
- `assembleResidual` ([main.cpp](../solver_density_cuda/main.cpp)) の RANS transport 呼出の直後に
  `condensationTransport_d_wrapper` を、explicit RK / segregated point-implicit に
  `condensationTimeIntegration_d_wrapper` を、roK/roOmega と同じ位置で差し込む。

#### 対流フラックスと二相 EOS (Phase 2)

Phase 2 で圧力が $p=(\rho-\rho g)RT$ となり気相へ弱く逆結合する。**SLAU を初手**とし、圧力流束は
$p$ を二相 EOS で評価するだけ、$\rho g$ は質量流束に乗る受動スカラーとして扱う ($\xi_g$ を陽に
Jacobian へ入れない)。保存量での圧力微分:

$$
\frac{\partial p}{\partial\rho}=\chi+\kappa K,\quad
\frac{\partial p}{\partial(\rho u_k)}=-\kappa u_k,\quad
\frac{\partial p}{\partial(\rho E)}=\kappa,\quad
\frac{\partial p}{\partial(\rho g)}=\xi_g
$$

($K=\tfrac12|\mathbf u|^2$、$\kappa,\chi,\xi_g$ は 本ドキュメントの「理論」節 5 節)。flux Jacobian では
運動量行に $n_k\,\partial p/\partial U_j$、エネルギー行に $u_n\,\partial p/\partial U_j$ が入り、
**既存の $(\gamma-1)$ が全て $\kappa$ に置換、$(\rho g)$ 列に $\xi_g$ が新規追加**される。Roe は固有構造が
$\kappa,\chi,\xi_g$ で変わるため一般 EOS Roe (Vinokur–Montagné 流) が必要で後段。

---

### 4. ソース項 (Phase 2 — 初期実装済, 安定性優先)

[condensationSource_d.{cu,cuh}](../solver_density_cuda/cuda_forge/condensationSource_d.cu) に実装。
device 側で核生成 $J$ (CNT × Iland)・臨界半径 $r_*$・成長 $dr/dt$ (Goodheart) を**現在セル状態から一度だけ
評価して freeze** (一温度 $T_v=T_d=T$)、相変化ソースを各モーメント残差へ加算する:

$$
S_{Q_0}=J,\ S_{Q_1}=Jr_{\rm nuc}+\rho Q_0\tfrac{dr}{dt},\ S_{Q_2}=Jr_{\rm nuc}^2+2\rho Q_1\tfrac{dr}{dt},\
S_g=\tfrac{4}{3}\pi\rho_l\big(Jr_{\rm nuc}^3+3\rho Q_2\tfrac{dr}{dt}\big)
$$

**核生成半径 $r_{\rm nuc}=C_{\rm nuc}\,r_*$ ($C_{\rm nuc}$=`COND_RNUC_FAC`=1.01)**: 臨界半径 $r_*$ ちょうどで生むと
全成長則が含む Kelvin 因子 $(1-r_*/\bar r)$ (Hertz–Knudsen では $p_v-p_d$、$p_d$ は Kelvin 平衡蒸気圧) が
$\bar r=r_*$ で **0 (不安定平衡)** となり、平均半径が $r_*$ に張り付いて成長が起動しない。さらに $\bar r$ が
$r_*$ を僅かに割ると $dr/dt<0$ → roQ1 が負帰還で 0 へ暴走崩壊し、$\bar r\to0$ で成長則が発散する。
わずかに超臨界 ($r_{\rm nuc}=1.01r_*$) で生むことで $\bar r>r_*$ から成長が立ち上がる。
**本体 kernel とヤコビアン (`cond_source_vector`) の両方で同一の成長ガード ($\bar r>r_*$ のみ成長・$dr/dt\ge0$) と
$r_{\rm nuc}$ を使う** (ヤコビアンだけガード無しだと亜臨界で過大な $\partial S/\partial Q_1$ が出て roQ1 を 0 に潰す)。

`res_ro<φ> += S_φ·V` で advection 残差へ加算。物性は [condensationProperties_d.cuh](../solver_density_cuda/cuda_forge/condensationProperties_d.cuh)。

#### 安定化 (初期実装の主眼)

- **$p_{sat}$ の Clausius–Clapeyron 外挿**: 液フィットは ~40K 以下で破綻するため $T<50$K は
  $p_{sat}(T)=p_{sat}(50)\exp(-(L/R)(1/T-1/50))$ で物理外挿 (低温で単調減少)。**45K クランプのまま psat を
  凍結すると過飽和 $S=p/p_{sat}$ が潰れ核生成が起きない**ため psat だけ C-C 外挿 (σ,ρ_l,L はクランプ)。
- **clamp/limiter**: $J$ 上限、$dr/dt<0$ は 0 (蒸発なし — `condEvaporation=1` では $S\le1$ で蒸発分岐に置換、理論節参照)、$\bar r\le r_*$ で成長停止、モーメント非負 (floor)、
  $g\le1$。**1 step の $\Delta g$・潜熱 $\Delta T$・蒸気枯渇を $\theta$ で律速**し全モーメントを同 $\theta$ で
  縮小 ($g{=}Q_3$ 整合保持)。
- **src_jac**: 初期実装はソース由来 0。時間項 $V/\Delta\tau$ + 移流 `transport_diag` の対角のみで安定化。
  ソースの T 依存 ($dJ/dT, d(dr/dt)/dT$) を含む完全自己抑制線形化は**後続**。

#### 検証 (case/34 run_0006, CPG)

収束 dry 場から restart、**NaN なし・bounded** ($g\in[0,1]$ max 0.24)・最小 T が dry 27.8→34.6K に上昇
(潜熱)・**g=0 域は dry と完全一致** (上流 P_cond/P_dry=1.00) を確認。g が過大 (paper ~0.75–2.3% に対し
22%) なのは収束 dry からの一斉 onset アーティファクト ($S\sim3\times10^4$ で全冷却セルが同時核生成)。
**定量一致は後続** (膨張流からの physical onset / レート較正 / 完全 src_jac 線形化)。

#### 検証 (case/16 Wyslouzil H2O, Fig.3 1kPa)

正しい Fig.3 ノズル形状 (`mesh/make_nozzle_fig3.py`, throat 全高 5mm・発散 95mm・壁間 1.8°・A/A\*≈1.6;
旧 `nozzle_H` は別形状で不一致だった) + SST。dry SST は実験 isentrope と中心線 p/p0 ±1.5% 一致。
**$r_{\rm nuc}=1.01r_*$ 修正前は液滴が $r_*$ (~0.3nm) に張り付き $g$ が蒸気の <3% で潰れていた** (核生成 $n\sim10^{19}$ は
出るが成長ゼロ; HK/Gyar が同一結果になるのが兆候)。修正後は **$g\to$ 蒸気のほぼ全量 (~90%)**・$\bar r$ が 5→25nm に
成長し潜熱バンプが出る。中心線 p/p0 は **下流 (x≳4cm) で全モデル実験 ±数%**。**onset 域 (x≈2cm) は Kantrowitz 無し
(HK/Gyar) が過早・オーバーシュート (0.44 vs exp 0.36)、Kantrowitz 有り (Kw+HK/Kw+Gyar) が onset を遅らせ実験に最良**
(Fluent UDF で Kw+HK が最適だった知見と一致)。残課題: x≈3cm のピークが実験よりやや高い (onset レート微調整)。

#### モデル切替 (config フラグ)

核生成/成長/補正は **enum + switch (device)** で種ごとに切替 (`CondSpeciesProps.model` で N2/H2O)。
加えて、以下を `solverConfig.yaml > condensation` のフラグで個別 on/off でき、感度評価に使う:

| フラグ | 既定 | 効果 |
| --- | --- | --- |
| `condKantrowitz` | 0 | 1 で核生成に **Kantrowitz 非等温補正** $J\to J/(1+\theta)$, $\theta=\frac{2(\gamma_v-1)}{\gamma_v+1}b(b-\tfrac12)$, $b=L/(R_vT)$ ($\gamma_v$=**凝縮種 (蒸気) の比熱比**: H2O 1.331 / N2 1.4, `CondSpeciesProps.cp/cv`)。0 は等温 CNT。|
| `condKantrowitzGammaMode` | 0 | Kantrowitz の $\gamma$ の取り方。0=蒸気 $\gamma_v$ (既定, 2026-09-10)、1=**旧挙動** (セル気相混合 $c_p/c_v$; H2O–N2 では ≈1.40 で $\theta$ が 17 % 過大)。A/B 用。|
| `condSonicModel` | 自動 | 凝縮セルの音速と $\gamma$。1=**二相 frozen** $c^2=\gamma_{2\phi}\,p/\rho$ (§5「二相 frozen 音速」)、0=旧挙動 (全蒸気気相 $\sqrt{\gamma_{mix}R_{mix}T}$、$g$ を無視)。未指定時は `input/condSonicResolve.hpp` が bcond 読込後に解決: **TP carrier H2O (`thermalMethod 2`, `condGasSpecies>=0`, `condModel 1`) かつ `condEquilibrium 0` かつ全 bcond が `inlet_Pressure`/`outflow`/`wall`/`slip`/`periodic` のとき 1、それ以外 0** (pure / CPG / `condEquilibrium 1,2` / `outlet_statPress`・`wall_isothermal` 等の ghost 再構築が全蒸気 EOS の境界は未検証のため旧式; 2026-09-10)。明示 1 は未検証構成でも従うが起動ログに警告。CPG 分岐 (`thermalMethod 0`) にはキーを 1 にしても効かない。|
| `condGrowthModel` | 0 | 0=既定 (H2O: Hertz–Knudsen 質量律速 / N2: Goodheart)、1=**Gyarmathy** 熱伝導律速 $\frac{dr}{dt}=\frac{kRT^2}{\rho_lL^2}\ln S\frac{1-r_*/r}{r(1+3.18Kn)}$ (極超音速ノズル凝縮で標準)。|
| `condGyarmathyC` | 3.18 | Gyarmathy の Knudsen 補正係数 $1/(1+C\,Kn)$。小さいほど成長速く onset 早、大きいほど遅。標準 3.18。|
| `condTwoTemp` | 0 | 1 で **液滴温度 $T_d$** を準定常 Hill バランスで解き Hertz–Knudsen 成長の駆動力 $p_v-p_d(T_d)$ に反映 (自己加熱で成長↓)。theory.md §4 参照。Gyarmathy 経路には非適用。|
| `condEvaporation` | **1** | 1 (既定, 2026-08-18 検証後に既定化) で **蒸発** ($S\le1$ のセルで負成長 λ スケール更新 + 液滴消滅)。理論節「統一駆動力と蒸発」参照。0 は旧挙動 $dr/dt<0\to0$ (液相が過熱域で凍結)。|
| `condEvapRmin` | 1e-9 | 完全蒸発とみなす $r_{30}$ [m]。$r_{30}<2r_{\min}$ で 4 モーメントを一括 0。|
| `condEvapKelvin` | 0 | 1 で蒸発駆動力に Kelvin 項 $p_d=p_{sat}e^{K_e/r_{30}}$ を含める (既定は平面 $p_{sat}$; 正帰還回避)。|

| `condEquilibrium` | 0 | 1 で **平衡凝縮** (局所瞬時平衡: 各セルで $p_v=p_{sat}(T)$ となる液相分率 $g_{eq}$ へ緩和)。核生成・成長・モーメントは使わない (Q0–Q2 ソース 0)。理論節「平衡凝縮」参照。|
| `condEqRelax` | 1.0 | 平衡緩和の 1 step あたり係数 $\alpha$: $S_g=\alpha\rho(g_{eq}-g)/\Delta t_{loc}$ (θ 律速 $\Delta g\le$5e-3, $\Delta T\le$1 K は共通)。|

蒸発の実装は `cond_evap_rate` / `cond_evap_source` ([condensationSource_d.cuh](../solver_density_cuda/cuda_forge/condensationSource_d.cuh))
と本体 kernel の蒸発分岐 (`condensation_source_d`)、消滅硬クランプは `cond_realizability_clamp_d`
([condensationTransport_d.cu](../solver_density_cuda/cuda_forge/condensationTransport_d.cu))。診断出力
`condS_<s>` (過飽和 $S$)・`condDrdt_<s>` (成長率, 負=蒸発)・`condR30_<s>` ($r_{30}$, 蒸発分岐のみ) を
`res_*.h5` に書く。0-D 単体テスト `tests/unit/test_cond_evaporation.cpp` (HK 解析解・実現可能性・消滅・
律速・Kelvin 質量収支・S>1 不発火・N2 経路) を持つ。

実装は [`condensationSource_d.cuh`](../solver_density_cuda/cuda_forge/condensationSource_d.cuh) の
`cond_nucleation` (Kantrowitz) と `cond_growth` (Gyarmathy 分岐)。Gyarmathy は N2 Goodheart と前因子
$kRT^2\ln S/(\rho_lL^2)$ を共有し Knudsen 内挿のみ $1/(1+3.18Kn)$ に差し替えた形。carrier では Kn 用の
平均自由行程に全圧 $p$ を用いる (希薄蒸気で $p_v$ を使うと Kn 過大になるのを回避)。

---

### 5. 一温度 二相 EOS の温度逆算 (Phase 2)

#### 実装ロードマップ

論文 (Lin 2014) は気相を**完全理想気体 (calorically perfect, $c_p$ 一定)** で扱う。これに合わせ:

1. **初期実装 = CPG 一温度 二相 EOS** (既定)。`thermalMethod 0` 分岐に追加。case/34 の既存初期場
   (CPG 基準 roe) がそのまま使え、TP 初期化問題を回避。**$L(T)$ の温度依存は入れる**。
2. **後続 A = 液滴温度独立** ($T_d$): Hill のエネルギーバランスで $T_v\ne T_d$ の 2 温度モデル。
3. **後続 B = thermally perfect 凝縮**: 気相 $e_v(T)$ を NASA-9 thermo で評価する版
   (`thermalMethod 2` 分岐、`cond_T_from_e_onetemp`、実装・unit 検証済)。$T_d$ 独立と組み合わせ可。

両分岐 (CPG/TP) を [condensationEOS_d.cuh](../solver_density_cuda/cuda_forge/condensationEOS_d.cuh)
に持ち、case/34 は CPG (`thermalMethod 0`) で進める。以下は CPG (初期実装) を基準に記述する。

#### CPG 一温度 二相 EOS (`thermalMethod 0`, 既定)

[dependentVariables_d.cu](../solver_density_cuda/cuda_forge/dependentVariables_d.cu) の
`thermalMethod==0` 分岐に追加 (`cond_T_from_e_cpg`)。気相 $e_v(T)=c_v T$ ($c_v=c_p/\gamma$ 一定)、
液相 $e_l=e_v+R_vT-L(T)=c_pT-L$ (∵ $h_v-h_l=L$, $h_v=e_v+R_vT$, $h_l\approx e_l$) より
$e=(1-g)e_v+ge_l=(c_v+gR_v)T-gL(T)$ (論文 Eq.18, $g_0{=}0$)。$T$ を Newton で反転、
$p=(1-g)\rho R_v T$ ($R_v=(\gamma-1)c_v$)。

#### 一温度近似 (T_v=T_d=T)

初期実装は気相温度 $T_v$ と液滴温度 $T_d$ を分けず $T_v=T_d=T$ とする。**$T_d$ は輸送変数に
追加しない** (Hill $T_d$ 式・$(T_v,T_d)$ 2 変数局所 Newton は後続拡張)。`cond_T_from_e_onetemp` は
T 1 変数 Newton で、拡張時に 2×2 Newton へ置換できる設計。

#### 混合内部エネルギーと温度反転

保存量から $e_{in}=\rho e/\rho-\tfrac12|\mathbf u|^2$ を作り、

$$
G(T)=(c_v+gR_v)T - g\,L(T)-e_{in}=0,\quad G'(T)=(c_v+gR_v) - g\,L'(T)
$$

をセルごとに Newton で解く ($L'$ は数値微分)。$g$ は総液相質量分率 $\sum_s \rho g_s/\rho$ ($g_s$ は
device `rog` 配列、`condensationInit_d` が構築)。$g<10^{-12}$ で従来 CPG 経路 ($T=e_{in}/c_v$) を呼び、
**単相 CPG に bit 同一**を保証 (圧力も $(1-g)=1$ で従来 $\rho R T$ と一致)。$\rho e$・$Ht$ は
$e_{mix}=(c_v+gR_v)T-gL$ で再構成、音速は下記「二相 frozen 音速」。

TP 版 (`thermalMethod 2`, 後続 B) は $e_v(T)$ を NASA-9 thermo で評価する以外は同形
($e=e_v(T)+gR_vT-gL$, `cond_T_from_e_onetemp`)、$g=0$ で `thermo_T_from_e` に厳密縮約。

#### 二相 frozen 音速 (2026-09-10, `condSonicModel: 1`)

流束 (SLAU の $\hat c$)・局所 $\Delta t$・block-DPLUR の音響固有値が使う音速は、上の一温度二相 EOS
$p=p(\rho,e,g)$ と整合した **frozen (相変化なし・液滴は気相と同速度・同温度) の等エントロピー微分**
$c^2=(\partial p/\partial\rho)_{e,g}+(p/\rho^2)(\partial p/\partial e)_{\rho,g}$ で定義する。
$p=\rho R_{eff}T$、$e=e_{gas}^{全蒸気}(T)+g(R_wT-L(T))$ から

$$
c^2=\gamma_{2\phi}\,R_{eff}\,T=\gamma_{2\phi}\frac{p}{\rho},\qquad
\gamma_{2\phi}=\frac{c_{p,2\phi}}{c_{v,2\phi}},\quad
c_{p,2\phi}=c_{p}^{全蒸気}-g\,L'(T),\quad c_{v,2\phi}=c_{p,2\phi}-R_{eff}
$$

- carrier (H2O in N2): $R_{eff}=R_{mix}-gR_w$ ($R_{mix}=\sum_iY_iR_i$ は総水を蒸気と数えた混合)。
- pure TP: $R_{eff}=(1-g)R_v$、$c_p^{全蒸気}=c_p(T)$。**CPG 分岐 (`thermalMethod 0`, pure N2 case/34) は旧式のまま**:
  `n2_latent` の多項式は低温 (例 45 K) で $L'>0$ ($c_l=c_{p,v}-L'<0$, 熱力学的に不整合) となり $g=0.2$ で $c_{v,2\phi}<0$ になる
  (codex 指摘 2026-09-10)。潜熱フィットの整合修正が先 (plan §5.1)。
- $L'(T)=dL/dT=c_{p,v}-c_l$ (Kirchhoff) なので $c_{p,2\phi}$ は「$g$ ぶんの蒸気 $c_{p,v}$ を液 $c_l$ に置き換えた」混合比熱
  ($c_{p,2\phi}=(1-g)c_{p,v}+gc_l$ に pure で一致)。$c_{v,2\phi}$ は温度反転 Newton の $de_{mix}/dT$ と同一。
- $g=0$ で $\gamma_{mix}R_{mix}T$ (従来) に厳密に戻る。`condensation: 0` はコード経路不変。
- 旧実装 (`condSonicModel: 0`) は $g$ を無視した全蒸気気相 $\sqrt{\gamma_{mix}R_{mix}T}$ で、$g$ が大きいほど
  $c$ を過大評価した (H2O–N2 Wyslouzil の出口 $g\approx0.011$ で +1.6 %、pure N2 case/34 の $g\sim0.2$ では +10 % 超)。
- per-cell `gamma` (block-DPLUR の $\kappa=\gamma-1$, TP 出口 BC) も $\gamma_{2\phi}$ にする。固定 $g,Y$ では
  $(\partial p/\partial e)_{\rho,g,Y}/\rho=R_{eff}/c_{v,2\phi}=\gamma_{2\phi}-1$ なので、一般 EOS Jacobian の $\kappa$, $\chi=c^2-\kappa h$ は
  **固定 $g,Y$ の frozen 近似**として整合する。ただし NS 保存量を更新してから `rog` を別更新する分離解法なので、$\rho g$ 列
  ($\partial p/\partial(\rho g)|_{\rho,\rho e}=-R_wT+\kappa(L-R_wT)$) を含む厳密 Jacobian ではない (LHS 近似; 収束経路にのみ効く)。
  `cp` 配列 (粘性・熱伝導の Pr 換算用) は全蒸気 $c_p$ のまま。
- 防御: $c_{v,2\phi}\le0.05\,c_{p,2\phi}$ か $c^2\le0$ のセルは旧式にフォールバック (H2O は $L'<0$ で起きない; pure N2 TP を明示 ON した場合の保険)。
- 境界: `outflow` (超音速) は内部 `sonic` をコピーするので二相音速がそのまま出る。TP 亜音速 `outlet_statPress` の ghost 再構築と
  node 等温壁ピン (`pin_wall_node_temperature_d`) は全蒸気 EOS のままなので内部と ghost の音速が不連続 (既存近似、未対応)。
- 平衡音速 (相変化を許す $S=1$ 拘束下の微分) は別物で未実装。`condEquilibrium: 2` は $g$ を状態量として再決定する経路で
  frozen 音速の妥当性は未検証 → 既定 OFF。

実装 `cond_twophase_sonic` ([condensationEOS_d.cuh](../solver_density_cuda/cuda_forge/condensationEOS_d.cuh))、
単体テスト `tests/unit/test_cond_sonic.cpp` (有限差分 $(\partial p/\partial\rho)_e+(p/\rho^2)(\partial p/\partial e)_\rho$ との一致 ≤1e-7、T×Y_w×g 掃引 84 状態で $c_{v,2\phi}>0$、`resolveCondSonicModel` の選択条件) と `tools/test_eos_jacobian.cpp` mode 2 (固定 $g,Y$ 二相状態で実 `accumulate_split_jacobian_cf` が流束 FD と 2e-8、float32 で 3e-7)。
計画・検証は [plans/active/condensation-kantrowitz-gamma-twophase-sonic.md](../plans/active/condensation-kantrowitz-gamma-twophase-sonic.md)。

**検証 (2026-09-10, case/16 Wyslouzil 2D node/cell SST 凝縮, `run_0331`–`run_0342`)**: 旧キー回帰は反復ノイズ床以内
(node 1.3e-5 vs 床 1.6e-5, cell ≤1.5 倍)、凝縮 OFF はコード経路不変。**γ_v 修正で onset 22.96 → 22.51 mm (−0.45 mm 上流)、
onset 帯の壁 p/p0 偏差 −4.8 → −4.5 %** (下流 +4.8 → +5.1 %)。**二相音速は凝縮帯の $c$ を −1.3 % (最大 −1.63 %) 下げるが、今回の条件では
収束場の差は小さい** (参照との差 $g$ 1.4e-5, $U_y$ 1.1e-5, $P$/$T$ ≤2e-6 = 同一バイナリ反復ノイズと同水準)。凝縮帯は流れ方向に M≥1.3 だが
SLAU は面法線速度で Mach を作るので横向き面は亜音速のまま (codex 計測: 凝縮帯内部面 16,726 のうち 8,304 面が両側とも法線亜音速) で
圧力流束は $c$ に依存する。差が小さい理由は特定していない (「全風上で影響 0」ではない)。
出力の Mach は二相 $c$ で定義されるので出口 M は 1.605 → 1.631 (+1.7 %)。実カーネルの `sonic`/`gamma` は host 式と
3e-7 / 5e-7 で一致。全 node run で `check_convergence` PASS (48000 step) と報告量の時系列 STEADY。

#### 検証 (host unit test 済)

- **(b) g=0 厳密縮約**: CPG `cond_T_from_e_cpg(g=0)` が $e_{in}/c_v$ と完全一致 (diff 0、40–290K)。
  TP `cond_T_from_e_onetemp(g=0)` も `thermo_T_from_e` と一致 (diff 0)。
- **(c) g>0 安定・物理**: $e_{in}=c_v\cdot50\text{K}$ 固定で $g$ を 0→0.08 に上げると潜熱で $T$ 単調上昇
  (50→72.3K、論文 40→56K と同オーダー)、Newton 残差 ~$10^{-9}$、$p=(1-g)\rho RT$ 正。CPG と TP はほぼ同値
  (N2 は $c_p$ 平坦)。
- 既存 `thermalMethod 0/2` 単相経路は未改変 (Phase 1 の run_0004 は不変)。

#### 後続拡張

- **A 液滴温度独立**: $e=(1-g)e_v(T_v)+g\,e_l(T_d)$ と Hill $T_d$ 式を組み、$(T_v,T_d)$ の 2 変数局所
  Newton へ ($T_d$ は局所量で輸送変数にはしない)。
- **B thermally perfect**: 上記 TP 版を既定化。多成分凝縮では $e=e_v-\sum_s g_s L_s(T)$、$R_v$ も気相混合で一般化。

---

### 6. 連成・陰解法 — 分離 (segregated)

モーメント輸送は **NS (気相 5×5 ブロック) と分離して解く**。roK/roOmega・化学種と同じく
block-DPLUR 後段の point-implicit 更新に乗せる (`applySSTPointImplicit` 同型の
`applyCondensationPointImplicit`, [update_d.cu](../solver_density_cuda/cuda_forge/update_d.cu))。
これは論文の fractional-step (均質対流部 → 凝縮ソース部) とも整合する。

point-implicit 対角 (時間 + 移流、Phase 1 は source=0):

$$
D = \underbrace{\frac{V}{\Delta\tau}}_{\text{時間項}} + \underbrace{\sum_f\frac{\max(\dot m_f,0)}{\rho}}_{\text{移流項}\,=\,\texttt{transport\_diag}},\qquad
\delta(\rho\phi)=\frac{\text{Res}}{D}
$$

Phase 2 の二相 EOS による気相逆結合 ($p$ が $g$ 依存) は密結合せず**毎反復の従属変数 $T/P$ 再計算で
渡す loose coupling**。圧力微分 $\xi_g$ の陰解法対角への取込みと気相+モーメント密結合 (block) は
収束不良時のフォールバックとして検討する。

---

### 7. 境界・初期・残差

- **境界**: 入口は dry (液相モーメント=0)。`ransBoundary_d.cu` の Dirichlet(=0)/Neumann(zero-grad)/
  wall(zero-grad) パターン ([ransBoundary_d.cu](../solver_density_cuda/cuda_forge/ransBoundary_d.cu)) を
  手本に `condensationBoundary_d.cu` を新設。inlet=0 固定、outlet/slip/axis=zero-gradient、壁は当面
  droplet 付着無視で zero-gradient。
- **初期**: `setInitial.hpp` の H2D コピー対象にモーメントを追加 (既定 0)。restart 読込
  (`readValueHDF5`) は roK/roOmega 同様、無ければ 0 フォールバック。
- **残差**: `main.cpp` の `residualEquationNames()`/`gatherResidualSnapshot()` に `rms_rog_<s>` 等を
  条件付き追加 → CSV 列に出力。

---

### 7b. 多成分燃焼ガス中の H₂O 凝縮 — carrier を擬似種に畳む運用 (2026-08-17)

燃焼ガス (N₂/CO₂/O₂/H₂O) の H₂O 凝縮では、carrier 全種を独立種にする必要はない (多成分 TP × 陰解法の
結合不安定・種数分の輸送コスト)。**H₂O 以外を NASA-9 の質量分率線形混合で 1 つの擬似種 `MIXDRY` に
畳み、H₂O だけ独立種**とする 2 種 TP (`thermalMethod 2`, `species: [MIXDRY, H2O]`, `condModel 1`,
`condGasSpecies 1`) で解く。混合則が線形なので dry の熱力学は全種独立と厳密に同じ。設計ツール側の
実装は [methods/design/overview.md](design/overview.md#ガスモデル-semi-perfect-nasa-9cea-frozen-組成--forge_designgas)
(`evaluate.tp_species: split_h2o`)、検証は
[plans/accepted/tooling-nozzle-tp-split-h2o-condensation.md](../plans/accepted/tooling-nozzle-tp-split-h2o-condensation.md)
(Wyslouzil fig3 で N₂ 擬似種 + H₂O が既存 `[N2, H2O]` CPG 結果を再現、イソブタン M4.2 H₂O 5 %)。
気相 thermo の低温側は forge の `Tlo` クランプ (cp 凍結) が効く。

### 8. 精度・無次元化 (Phase 2)

Phase 1 は source=0 ゆえモーメントは常時 0 で精度無関係。Phase 2 で
$\mu_n=Q_n/(N_{ref}r_{ref}^n)$ ($N_{ref}=10^{18}$/kg, $r_{ref}=1$ nm) を導入し、保存量を $\rho\mu_n$ と
する。**まず全 float**で進め、無次元化で桁を O(1) に抑える。$r=(\mu_3/\mu_0)^{1/3}r_{ref}$ 等の
モーメント間比も無次元量どうしで安定。float で桁落ちが顕在化したら、その時点でモーメントの
double 化 (混合精度 cond storage、`flow_float` とは別の `cond_float=double`) をフォールバックとして
導入する。
