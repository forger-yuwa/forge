# SST 乱流モデルの理論メモ

## 1. 目的

本章では forge に導入する RANS ベースの **Menter SST**
($k$-$\omega$ Shear Stress Transport) モデルの支配方程式と、
本コードで採用する前提を整理する。

最終的な導入対象は SST 全系だが、初回マイルストーンは
scalar transport 配線確認のため `k`, `\omega` の advection-only とする。

初期スコープは次のとおりとする。

- 低 Re の壁解像型 SST
- 圧縮性 NS ソルバ上の保存形輸送
- まず generic scalar transport 基盤上で `k`, `\omega` を扱う
- 初回は advection-only を対象とし、その後 diffusion / source を追加する
- 軸対称は既存幾何の上で advection を確認し、SST 固有幾何項は別 plan で拡張
- 陽解法 (`timeIntegration = 1, 3, 4`) を対象とし、陰解法は別 plan

## 2. 保存変数

既存の 5 保存変数

$$
Q_{NS} = [\rho, \rho u_x, \rho u_y, \rho u_z, \rho E]^T
$$

に加えて、SST では次の 2 変数を導入する。

$$
Q_{turb} = [\rho k, \rho \omega]^T
$$

ここで $k$ は乱流運動エネルギー、$\omega$ は比散逸率である。

## 3. 輸送方程式

圧縮性保存形では、2 つの追加方程式を generic scalar transport の特殊例として扱う。
最終形は次である。

$$
\frac{\partial (\rho k)}{\partial t}
+ \nabla \cdot (\rho \mathbf{u} k)
= \nabla \cdot \left[
\left(\mu + \sigma_k \mu_t\right) \nabla k
\right]
+ P_k - \beta^* \rho \omega k
$$

$$
\frac{\partial (\rho \omega)}{\partial t}
+ \nabla \cdot (\rho \mathbf{u} \omega)
= \nabla \cdot \left[
\left(\mu + \sigma_\omega \mu_t\right) \nabla \omega
\right]
+ \alpha \frac{\omega}{k} P_k
- \beta \rho \omega^2
+ D_{\mathrm{cross}}
$$

ここで $P_k$ は乱流生成項、$D_{\mathrm{cross}}$ は SST 特有の
cross-diffusion 項である。

初回マイルストーンでは右辺を落とし、

$$
\frac{\partial (\rho \phi)}{\partial t}
+ \nabla \cdot (\rho \mathbf{u} \phi) = 0
$$

で $\phi \in \{k, \omega\}$ を輸送し、保存量配線と境界条件・時間積分を先に確立する。

### 3.1 source term の全体像

本来の SST では、対流・拡散だけでなく生成項と散逸項が入る。実装上は
`k` と `\omega` を独立な scalar と見なしつつ、source は同じ 2 本の式へ加える。

標準形は次である。

$$
S_k = P_k - \beta^* \rho \omega k
$$

$$
S_\omega = \alpha \frac{\omega}{\max(k, k_{\min})} P_k
- \beta \rho \omega^2
+ D_{\mathrm{cross}}
$$

ここで

$$
P_k = \min\left(\tau_{ij} \frac{\partial u_i}{\partial x_j},\; 10\,\beta^*\rho\omega k\right)
$$

であり、$\tau_{ij}$ は Reynolds stress を近似する eddy-viscosity closure から作る。
SST では stagnation point などで $P_k$ が過大になりやすいため、上式のように
production limiter を入れるのが基本である。

cross-diffusion 項は、SST が $k$-$\varepsilon$ 側と $k$-$\omega$ 側を混ぜるために入る補正で、
一般には

$$
D_{\mathrm{cross}} = 2 (1 - F_1) \rho \sigma_{\omega 2}
\frac{1}{\max(\omega, \omega_{\min})}
\nabla k \cdot \nabla \omega
$$

と書ける。`F_1` は壁近傍で 1 に寄る blending function で、壁から離れるほど 0 に近づく。
したがって cross-diffusion は壁近傍では消え、自由流側で効いてくる。

#### 3.1.1 $\omega$ 生産項の 3 つの流儀と、forge 実装の現状 (2026-08-12 追記)

$\omega$ 生産項には文献上 3 つの書き方があり、**forge の実装は上式 (本文の標準形) と一致していない**。

| # | 形 | 出典・位置づけ |
| --- | --- | --- |
| (a) | $P_\omega=\alpha\dfrac{\omega}{\max(k,k_{\min})}P_k$ | 本文の標準形 (上式)。$\nu_t=k/\omega$ の限りで (c) と等価 |
| (b) | $P_\omega=\alpha\rho S^2$ | Menter 2003 論文の記載。**NASA TMR が誤記 (typo) と明記** |
| (c) | $P_\omega=\alpha\dfrac{\dot P_k}{\nu_t}$ | **正しい形** (TMR)。$\dot P_k$ は production limiter 適用後の値、$\nu_t$ は SST の (shear-limited) 実値 |

**一次情報** ([NASA Turbulence Modeling Resource, SST](https://tmbwg.github.io/turbmodels/sst.html)):

> "a typographical error existed in this paper that was subsequently corrected by the authors.
> In the omega equation (2nd part of eqn (1) in the paper), the production term was incorrectly
> given as α ρ S² (using the paper's notation). Instead, it should have read **α Ṗₖ / νₜ**."
>
> "the production limiter is used for **both** k and omega equations, and the constant is changed
> from 20 to 10."

**forge の現状**: [`ransSource_d.cu`](../../solver_density_cuda/cuda_forge/ransSource_d.cu) は
`Pw = alpha * rho * S_prod` = **(b) の誤記形**を実装している。また production limiter は
$k$ 側にのみ適用され $\omega$ 側には適用されていない。

(b) と (c) は $\nu_t$ が shear limiter で抑えられている領域 (近壁・強せん断) で乖離する。
本件の実害は node 壁関数の第一内層で観測されている
([plan: turbulence-node-wf-omega-source](../../plans/accepted/turbulence-node-wf-omega-source.md))。
**一般内部場を (c) へ直す変更は cell/node 全ケースに影響するため、別 plan として起票し
回帰を取ってから行う** (本文のこの節は現状の記述であり、修正されたら更新すること)。

### 3.2 blending による係数切替

SST では係数を 2 組持ち、$F_1$ で混合する。

$$
\sigma_k = F_1 \sigma_{k1} + (1 - F_1) \sigma_{k2},
\qquad
\sigma_\omega = F_1 \sigma_{\omega1} + (1 - F_1) \sigma_{\omega2}
$$

$$
\alpha = F_1 \alpha_1 + (1 - F_1) \alpha_2,
\qquad
\beta = F_1 \beta_1 + (1 - F_1) \beta_2
$$

`F_1` と `F_2` は壁距離 $y$、$k$、$\omega$、勾配を使って定義される。
forge では既存の `wall_dist` をその入力として再利用する。

source と blending の実装上の注意は次のとおりである。

- $k$ と $\omega$ は正値制約を壊しやすいので、source でも下限値を入れる
- $\omega$ の発散は壁近傍で特に起こりやすいため、$\omega_{\min}$ を明示する
- production limiter と blending を分離して実装し、検証しやすくする

## 4. 渦粘性

SST の渦粘性は、$k$-$\omega$ 系をベースに blending function を通して
壁近傍と自由流側を切り替える。

代表形は次である。

$$
\mu_t = \rho \nu_t,
\qquad
\nu_t = \frac{a_1 k}{\max(a_1 \omega, S F_2)}
$$

ここで $S$ はひずみ速度の代表量、$F_2$ は壁距離を含む blending function
である。forge では既存の `wall_dist` を再利用し、追加の壁距離ソルバは
導入しない。

## 5. 係数と blending

SST では $k$-$\omega$ 壁近傍モデルと、自由流側の $k$-$\varepsilon$ 的な
振る舞いを blending function $F_1$, $F_2$ で接続する。

実装では次を守る。

- 係数はコード内に固定値として持つ初期実装とする
- 将来、YAML から係数上書きを許す場合も、既定値は Menter SST を維持する
- $k$, $\omega$, $\mu_t$ の positivity を壊さないために下限値を設ける

## 6. 壁条件

初期スコープは壁関数なしの壁解像型とする。SST の壁条件は、
`k` と `\omega` を壁で別々に閉じる必要がある。

### 6.1 no-slip wall / wall_isothermal

最小構成では、壁面で乱流エネルギーを消し、比散逸率は壁距離から与える。

$$
k_w = 0
$$

$$
\omega_w = \max\left(\frac{60\,\nu_w}{\beta_1 y_w^2},\; \omega_{\min}\right)
$$

ここで $y_w$ は最初のセル中心から壁までの距離、$\nu_w$ は壁面近傍の動粘性係数である。
`wall_isothermal` でも turbulence の考え方は同じで、温度境界は熱方程式側の責務とし、
SST 側は同じ $k_w$, $\omega_w$ を使う。

SST の理論上の ghost-cell closure は、壁値 $\phi_w$ を first-order の線形補間で与える形で書ける。
壁面がセル中心と ghost セルのちょうど中点にあるとき、

$$
\phi_w = \frac{\phi_c + \phi_g}{2}
$$

なので、Dirichlet 値を ghost へ戻す式は

$$
\phi_g = 2\phi_w - \phi_c
$$

である。したがって SST では

$$
k_g = -k_c
$$

$$
\omega_g = 2\omega_w - \omega_c
$$

となる。つまり、`k` は wall value が 0 なら ghost 側は負値になるが、これは「負の乱流エネルギーを物理的に意味づける」ものではなく、壁値 $k_w=0$ を中点補間で満たすための数値的な反射である。

壁からの距離が中点対称でない場合は、一般には距離比で係数が変わる。
現行コードは `kb` / `omegab` を壁値として持ち、ghost にはこの反射式を適用する形に更新した。

### 6.2 slip / axis / periodic

これらの境界は、SST scalar に対しては基本的に zero normal flux とする。

- `slip`: 法線方向勾配を持たせず、ghost へ copy する
- `axis`: 軸対称の対称条件として mirror / copy を使う
- `periodic`: 対向面の値をそのまま使う

このグループでは、壁のような強い散逸条件は入れない。

SST の立場では、`slip` / `axis` / `periodic` は壁面 Dirichlet ではなく、
\(\partial n\) に対する零流束条件として扱う。したがって ghost 値は通常

$$
\phi_g = \phi_c
$$

の copy で十分である。

### 6.3 inlet / outlet

流入は自由流乱流量を直接与える。実装上は boundary config から `k` と `omega` を受け取り、
必要に応じて turbulence intensity と length scale から換算する。

$$
k_\infty = \frac{3}{2} (U_\infty I)^2,
\qquad
\omega_\infty = \frac{\sqrt{k_\infty}}{C_\mu^{1/4} L_t}
$$

流出は 1 次外挿または zero-gradient を基本とし、逆流時のみ安定化のための clamp を入れる。

### 6.4 実装メモ

forge では `wall_dist` を既に持っているので、壁条件と source の両方で再利用する。
SST の壁条件は、単なる ghost copy ではなく、壁距離に基づく `\omega_w` の決定が必要になる。
実装上はまず boundary-plane 側で `kb`, `omegab` を持ち、そこから ghost セルへ

$$
k_g = 2k_w - k_c,
\qquad
\omega_g = 2\omega_w - \omega_c
$$

を与えるか、あるいは wall カーネルの中で ghost 値を直接埋める。
現行のコードは `wall` / `wall_isothermal` の boundary value を `kb` / `omegab` として保持し、
wall kernel 側で ghost へ反射して入れる。
等距離でない場合は、壁からセル中心までの距離を $d_c$、壁から ghost 中心までの距離を $d_g$
として

$$
\phi_g = \frac{(d_c + d_g)\phi_w - d_g\phi_c}{d_c}
$$

を使う。
したがって、壁条件は source よりも先に正しく閉じる必要がある。

### 6.5 automatic (enhanced / y⁺ 非依存) wall treatment

§6.1 の `\omega_w = 60\nu/(\beta_1 y^2)`, `k_w = 0` は **第一セルが粘性低層 (y⁺≲1) に
ある**ことを前提とした wall-resolved 型である。第一セルがバッファ層 (y⁺≈5〜30) や
対数層 (y⁺≳30) に落ちると、この `\omega_w` は過大評価になり、壁せん断応力 $\tau_w$ を
分子勾配 $\mu\,\partial U/\partial n$ から評価する運動量境界も $\tau_w$ を過小評価する。

automatic wall treatment は、粘性低層・バッファ層・対数層を **1 つの定式で滑らかに繋ぎ**、
第一セルの y⁺ 位置に依存せず妥当な $\tau_w$, $\omega_w$ を与える (Menter, *automatic
near-wall treatment*)。**`wallTreatmentSST` の既定は 1 (automatic)** (2026-06-28 に 0→1 へ変更, user 指示)。
`0` で §6.1 の wall-resolved 型 ($60\nu/\beta_1 y^2$, $k_w=0$) に戻せる。**注**: 既定変更により y⁺~1 の
wall-resolved 前提で検証していた cell ケース (`case/26` $C_f$, flat-plate 回帰 等) は automatic に切り替わるため
要再検証。node-centered の Dirichlet 系既定もまとめて ON: `nodeWallDirichlet=1` (壁 no-slip Dirichlet, §6.5(e)・
plan discretization-node-wall-implicit-dirichlet)、`nodeKwfDirichlet=1` (第一内層ノード $k$ Dirichlet, §6.5(e))。

#### (a) 摩擦速度 $u_\tau$ — Reichardt 普遍速度則の逆解き

粘性低層・バッファ・対数を 1 式で繋ぐ Reichardt 則を用いる。

$$
u^+ = \frac{1}{\kappa}\ln(1+\kappa y^+)
      + 7.8\left[1 - e^{-y^+/11} - \tfrac{y^+}{11}e^{-y^+/3}\right]
$$

ここで $u^+ = U_t/u_\tau$, $y^+ = u_\tau y/\nu$, $U_t$ は壁接線速度の大きさ
(no-slip 壁では $U_t = |\mathbf{U}_c - (\mathbf{U}_c\cdot\hat{\mathbf n})\hat{\mathbf n}|$)、
$y$ は第一セル中心-壁距離 (`wall_dist`)、$\kappa = 0.41$。
$U_t$, $y$, $\nu$ は既知なので、

$$
f(u_\tau) = \frac{U_t}{u_\tau} - u^+\!\big(y^+(u_\tau)\big) = 0
$$

を Newton 法数回 (初期値 $u_\tau^{(0)} = \sqrt{\nu U_t / y}$、粘性則) で解く。
$y^+\to 0$ では $u^+ = y^+$ (粘性解) に、大 $y^+$ では対数則に自動収束する。

#### (b) $\omega$ 壁面境界 — Menter automatic ブレンド

粘性低層の解析漸近と対数層の値を二乗和の平方根で滑らかに繋ぐ。

$$
\omega_{vis} = \frac{6\nu}{\beta_1 y^2},
\qquad
\omega_{log} = \frac{u_\tau}{\sqrt{\beta^*}\,\kappa\,y},
\qquad
\omega_w = \sqrt{\omega_{vis}^2 + \omega_{log}^2}
$$

粘性低層では $\omega_{vis}$、対数層では $\omega_{log}$ が支配する。§6.1 の係数 60 ではなく
解析漸近の係数 6 を用いる ($y^+\to0$ の厳密漸近に一致させるため)。$\omega_w$ は全 $y^+$ で
妥当 ($y^+\to0$ で $\omega_{vis}$ 支配=壁解像値) なので、ghost 反射に加えて
**wall-adjacent セル値そのものを $\omega_c = \omega_w$ にピン留め**する (conserved $\rho\omega$ も整合)。
これにより渦粘性 $\mu_t = \rho a_1 k/\max(a_1\omega, S F_2)$ が近壁で適切に上限され、$k$ の暴走
(後述 (b')) を断つ。

#### (b') $k$ 壁面境界 — zero-gradient ($P_k$ の wall-function 化とセット)

wall-function では $k$ を **zero-gradient (Neumann, $\partial k/\partial n = 0$)** とするのが標準
(OpenFOAM `kqRWallFunction`、ANSYS automatic の $k$ zero-flux)。第一セルがバッファ・対数層に
あるとき $k$ は有限の対数層平衡値 $k \approx u_\tau^2/\sqrt{\beta^*}$ を取るべきだからである。

**ただし $k$ zero-gradient は近壁の生産項 $P_k$ の wall-function 化 (下記 (d)) とセットで初めて
成立する。** $P_k$ を標準 SST のまま**解像速度勾配**から計算したまま $k$ を zero-gradient にすると、
粗メッシュでは誤った解像生産で近壁 $k$ が暴走し $\mu_t$ 過大 → $u_\tau$/$C_f$ 過大になる
(実測 `case/26`: $k$ zero-gradient 単独で y⁺≈30 が $C_f/C_{f,\text{corr}}$ 0.91→1.80, y⁺≈10 が
1.10→1.68 と悪化した)。(d) の $P_k$ 置換と (b) の $\omega$ ピン留めを同時に入れると、$k$ は
自律的に平衡値へ収束し暴走しない。したがって automatic モードでは
$$
k_w = k_c,\qquad k_g = k_c \quad(\partial k/\partial n = 0)
$$
とする (wall-resolved モード=既定 0 は $k_w = 0$, $k_g = -k_c$ を保つ)。

#### (c) 運動量壁せん断 — 有効壁粘性

第一セルが粗いと分子勾配は $\tau_w$ を過小評価するため、モデル化した壁せん断
$\boldsymbol{\tau}_w = \rho u_\tau^2\,\hat{\mathbf e}_t$ ($\hat{\mathbf e}_t$ は接線速度方向)
を運動量残差に直接課す。等価な有効壁粘性は

$$
\mu_{eff} = \frac{\rho u_\tau^2\, y}{U_t}
$$

であり、$y^+\to0$ では $u_\tau^2 = \nu U_t/y$ より $\mu_{eff}\to\mu$ (wall-resolved に一致)、
対数層では $\mu_{eff} > \mu$ (wall function 化) となり、全層で連続に繋がる。
no-slip 壁では壁せん断は仕事をしない ($\boldsymbol{\tau}_w\cdot\mathbf{U}_w = 0$) ため、
エネルギー方程式の壁せん断仕事項と熱流束は不変である。

#### (d) wall-adjacent セルの乱流生産 $P_k$ — wall-function 値への置換

(b') の $k$ zero-gradient を成立させる要。wall-adjacent セルでは標準 SST の解像勾配生産
$P_k = \mu_t S^2$ を、普遍速度則から評価した wall-function 生産に置換する。乱流応力 (= 全応力
− 粘性応力) と壁面速度勾配の積として
$$
P_k = (\tau_w - \tau_{vis})\,\frac{\partial U}{\partial y}
    = \rho\frac{u_\tau^4}{\nu}\,g\,(1-g),
\qquad
g \equiv \left.\frac{du^+}{dy^+}\right|_{y^+_1}
$$
ここで $\partial U/\partial y = (u_\tau^2/\nu)\,g$, $\tau_w = \rho u_\tau^2$, $\tau_{vis} = \mu\,\partial U/\partial y
= \tau_w\,g$ を用いた。$g$ は (a) の Reichardt 則の微分 $du^+/dy^+$ を第一セル $y^+_1$ で評価。

- **粘性低層** ($y^+\to0$): $g\to1$ より $P_k\to0$ — せん断は全て分子粘性が担い乱流生産ゼロ。
  標準 SST の wall-resolved 極限に一致する (細メッシュで既定 0 と同結果)。
- **対数層** ($y^+\gg1$): $g\to 1/(\kappa y^+)$ より $P_k\to \rho u_\tau^3/(\kappa y)$ — 標準的な log 則生産。

この $P_k$ と (b) の $\omega$ ピン留めにより、$k$ 方程式の生産・消滅平衡
$P_k = \beta^*\rho k\,\omega_w$ から $k\to u_\tau^2/\sqrt{\beta^*}$ (対数層平衡値) へ自律収束し、
$k$ を直接固定せずとも暴走しない。検証は plan / `case/26` を参照 (3 帯 y⁺≈0.35/10/30 で
$C_f/C_{f,\text{corr}}$ がいずれも 0.9–1.0 に収まる)。

#### (e) 適用先 — cell は壁隣接第一セル / node は壁ノード+第一内層ノード

(b)〜(d) の壁関数処理は離散化で当てる位置が異なる。**生産 $P_k$ 置換 (d) と $\omega$ ピン (b) で扱いが分かれる**
のがポイント。

- **cell-centered**: 壁隣接「第一セル」(セル中心が壁から $\sim\Delta y/2$) が唯一の近壁 DOF。$P_k$ 置換・$\omega$
  ピン・$\omega$ 残差ゼロ化・$\omega$ 行 decouple をすべて第一セルに当てる。SU2 の Normal_Neighbor と等価で
  x_R が SU2 と一致 (検証済み)。
- **node-centered (median-dual)**: 壁面上に壁ノード ($y=0$) と第一内層ノードの **2 つの DOF** がある。
  - **$P_k$ 置換 (d) は壁ノード+第一内層ノードの両方**に当てる。壁ノードだけだと第一内層ノードが標準 SST
    解像生産 $P_k=\mu_t S^2$ のまま残り、淀み凹コーナーで $k$ が暴走する (実測 `case/18.backstep` 段差コーナー
    $\mu_t/\mu\!\approx\!6789$, $k\!\approx\!52$ → 修正後 $77$, $k\!\approx\!3$ = cell 同等)。第一内層ノードは
    $u_\tau$ 逆解きで既に選ぶ代表内部点 (壁内向き法線 $-\hat{\mathbf n}$ との $\cos$ 最大の内部隣接ノード,
    SU2 Normal_Neighbor 流) を使う。
  - **$\omega$ ピン (b)・残差ゼロ化・行 decouple は壁ノードのみ**に当てる (第一内層ノードの $\omega$ は通常 solved、
    壁ノードのピンが近傍を抑える)。$\omega$ を第一内層ノードに移すと、段差面と底壁で同じ第一内層ノードを共有する
    凹コーナーで $\omega$ ピン値が競合し、かつ「壁ノードからの $\omega$ 拡散アンカ」が外れて $\omega$ が崩壊し
    $k$ がさらに暴走する (実測で確認、不採用)。
  - **再付着 $\mu_t$ ピーク (残課題)**: node の壁ノードは cell 第一セルより壁から遠い ($\Delta y$ vs $\Delta y/2$) ため
    $\omega$ ピン値が $\sim1/4$ 低く、再付着せん断層の近壁で $k$ destruction が cell より弱く $\mu_t$ ピークが残る
    (`case/18.backstep` 再付着 $\mu_t/\mu\!\approx\!5800$ vs cell 990、局所 ~27 node)。場平均・$x_R$ は cell/SU2 と整合。
  - **オプション: node $k$ Dirichlet (`nodeKwfDirichlet=1`, SU2 `SetTurbVars_WF` 流)**: 第一内層ノードの $k$ を
    $k_{wf}=\omega_w\,\mu_{t,wall}/\rho$ ($\mu_{t,wall}=\nu(1/g-1)$, $g=du^+/dy^+$) にハード Dirichlet し、上記
    再付着ピークを除去できる (`case/18.backstep` $\mu_t/\mu$ max 5800→890 = cell 990 並)。**ただし副作用**: 非平衡な
    再付着点では $u_\tau\to0$ で $k_{wf}\to0$ となり再付着乱流を抑えすぎ、$x_R$ が伸びる (7.63→8.67、cell 7.95/SU2 7.89
    より長くなる)。場の清浄さ ($\mu_t$ ピーク除去) と $x_R$ 精度のトレードオフ。**既定 ON** (2026-06-28, user 指示;
    `nodeKwfDirichlet=0` で $k$ を解く prod-fix に戻すと $x_R$ が cell/SU2 寄り)。$\omega$ ピンは Dirichlet 化しても
    壁ノードのまま (第一内層へ移すと凹コーナーで $\omega$ ピン値が race し崩壊するため不可)。

  - **★ 既知の不整合 (2026-08-12, 未修正)**: 第一内層ノードは **$P_k$ だけ壁関数値に置換され、
    $\omega$ は解像ひずみ $S$ で解かれる**という非対称な状態にある。$\omega$ 生産は
    $P_\omega=\alpha\rho S^2$ で $P_k$ と結合していない (§3.1.1 の誤記形) ため、y+≫1 の
    未解像な離散勾配がそのまま $\omega$ を駆動する。実測 (case/26 平板 y+30, x=0.6):
    解像 $S$=1.075e5 に対し壁法則整合 $S_{\rm wf}=u_\tau^2g/\nu$=5.146e4 = **2.09 倍**
    ($S^2$ で 4.36 倍)。結果、第一内層 $\omega$ が**壁ノードのピン値より高く跳ねる**
    (1.146e5 → 2.645e5 → 9.635e4)。これは $\alpha\rho S^2$ と $\beta\rho\omega^2$ の局所平衡
    ($\sqrt{\alpha/\beta}S$=2.754e5, 実測比 0.960) そのもの。
    高い $\omega$ は $D_k=\beta^*\rho k\omega$ を増やして $k$ を下げ、第一内層では
    $\mu_t$ が strain-limited 枝 ($SF_2/(a_1\omega)$=1.288) にあるため分子経由で $\mu_t$ も下がる。
    node の $C_f$ が壁解像基準の 0.843 に留まる件の有力仮説。
    **cell では同じ DOF の $\omega$ がピンされていて $S$ に応答しない**のが差の出所。
    分離実験と修正候補は
    [plan: turbulence-node-wf-omega-source](../../plans/accepted/turbulence-node-wf-omega-source.md)。
  - **★ 分離実験の結果 (2026-08-12)**: この構成の寄与を 2×2 factorial で分離した
    ($\omega$ Dirichlet ピン × SST shear limiter 迂回、case/26 平板 y+27、
    外部相関 $C_f/C_{f,\mathrm{KS}}$、W–I 双対面で実際に加わった接線力ベース、
    判定帯 $x$=[0.640, 0.784] で全 run $Re_\theta\ge4000$):

    | | bypass=0 | bypass=1 |
    | --- | --- | --- |
    | pin=0 | 0.7615 | 0.7670 |
    | pin=1 | 0.8594 | **0.8864** |

    平均主効果 pin **+0.1087** / bypass **+0.0163** / 相互作用 **+0.0215**。
    → **$\omega$ 状態が最大因子**。limiter 迂回は単独では小さいが ($+0.0055$)、
    $\omega$ ピン時には $+0.0270$ と有意な正の相互作用を持つ。
    ただし**両方を入れても 0.8864** で壁解像基準 $0.94\pm0.02$ に届かない
    (中心まで 70% / 下限まで 79% の回復) ので、
    **この特定の $\omega$ Dirichlet 値と迂回の組合せでは説明しきれない**。
    残差分の候補は ①ピン値そのもの ②ピンする DOF/範囲 ③$P_\omega$・拡散・交差拡散など
    別の $\omega$ 機構 ④$k$–$\omega$ 結合。
  - **★ 壁応力の伝達は健全 (2026-08-12)**: 「モデル $\tau_w$ が内部ノード残差に届いていない」
    可能性は診断で否定された。AddTauWall 後の接線力は W–I 双対面上で意図した
    $\tau_w A_{\rm WI}$ になり、その同じ力が内部ノードの運動量残差へ加算される
    (壁ノード側の運動量残差は `nodeWallDirichlet_d.cu` でゼロ化されるので反力は内部側に残る)。
    再スケールが traction ベクトル全体に掛かる件も、**発達域では法線成分が接線の
    $2.7\times10^{-5}$** で実害なし (ただし前縁の壁∩対称面コーナーでは局所 0.34 に達するので
    「全域で無害」と一般化しない)。
  - **★ E3: $\omega$ 生産の壁法則整合化 (2026-08-13, 実験実施・不採用)**: 上記の非対称
    ($k$ は壁関数値・$\omega$ は解像ひずみ) を解消する試み。第一内層ノード**のみ**、
    $P_\omega$ の入力ひずみを $S_{\rm wf}=u_\tau^2 g/\nu$ に差し替える
    ($P_{\omega,\rm wf}=\alpha\rho S_{\rm wf}^2$)。$\omega$ の状態は free、limiter も通常式のまま
    (Dirichlet を使わない)。実装は `ransWallFunction_d.cu` が $S_{\rm wf}^2$ を `wf_sprod` へ格納し
    `ransSource_d.cu` が $\omega$ 生産の入力にのみ使う。$\alpha P_{k,\rm wf}/\nu_{t,\rm wf}$ の形は
    低 $y^+$ で 0/0 になるため $S_{\rm wf}^2$ を直接使う。
    **結果**: case/26 平板で $C_f/C_{f,\mathrm{KS}}$ 0.7615→**0.8445** ($\omega$ 2.62e5→1.35e5)。
    **case/40 ベルノズルでは統一評価 (`wall_tau_eval.py`) で $\tau_w$ が生産比 1.1954 倍、
    y+1 内部基準比 1.1414 に過大化**し (旧評価の「応答比 1.1226」は別評価系の値で、
    内部基準比 1.1414 とは**等号で結べない**)、
    事前に定めた採否ゲート (「case/26 の改善が case/40 の過大化を伴わない」) を満たさないため
    **既定化しない**。実験経路として `FORGE_WF_OMEGA_SOURCE=1` で残す
    (**node SST 壁関数のときだけ有効**。そうでないと `wf_sprod` が未初期化のまま参照され
    $\omega$ 生産を全域で消すため、ゲートは必須)。
    振れ幅は**各ケース内で同一定義の before/after** を取ると
    case/26 **+10.9%** ($C_f$/KS 系) / case/40 **+19.5%** (面積重み $\tau_w$ 系) で
    **同程度ではない** (旧記載「両ケース 11–12% 同程度」は撤回)。
    **両者は同一の汎関数ではない**ので大小は参考比較にとどめる。
    基準に対する絶対誤差は統一評価で **case/26 が −23.9%→−15.6% (不足側のまま)**、
    **case/40 が −4.51%→+14.14% (不足→過剰へ)** と**符号が違う**。
    → **今回試した 3 つの $\omega$ 介入 (Dirichlet ピン / limiter 迂回 / E3) では
    両ケースのゲートを同時に満たせなかった**。$\omega$ の輸送・ブレンド・適用範囲・$k$ との結合は
    未検証であり、「$\omega$ 側を出し切った」わけではない。次は独立性の高い代表点経路を診断する
    ([plan: turbulence-node-wf-representative-point](../../plans/accepted/turbulence-node-wf-representative-point.md))。
  - **診断の使い方**: `FORGE_WI_FORCE_DIAG=1` で `wi_ftan`/`wi_fnrm`/`wi_fnrm_abs`/`wi_ftan_res`
    を出力する (OFF なら確保も出力もされない)。`FORGE_WF_LIMITER_BYPASS` (−1/0/1) で
    limiter 迂回を $\omega$ ピンから独立に切れる (node SST 壁関数時のみ有効)。
    第一内層ノードのマスクは `wf_irep_flag` (`wf_pk>=0` は壁ノードにも立つので代用不可、
    cell では常に 0)。`FORGE_OMEGA_BUDGET=1` で $\omega$ 項別収支
    (`omg_prod`/`omg_dest`/`omg_cross`/`omg_trans`/`omg_axisym`)。
    `omg_trans` は**軸対称 $1/y$ ソース加算前**の残差なので、軸対称ケースでは `omg_axisym` と
    併せて見ること。**これらの env はすべて OFF で確保も出力もされない** (既定経路不変)。

#### (f) 断熱壁の熱的閉包 — Crocco 型回復温度 $T_{aw}$ (`sstThermalWallFunction`)

(a)–(e) は**運動量・乱流スカラー側のみ**の壁関数であり、エネルギー側の壁法則を持たない。
第一 DOF が対数層にある粗い壁メッシュでは、粘性散逸熱が乱流伝導で壁へ運ばれて成立する
**断熱回復温度が解像されず、壁温出力が $-200\,\mathrm{K}$ 級に冷える** (case/40 ベルノズル:
wf=1 壁温 1193 K vs SU2 壁関数 1422 K / y+≈1 low-Re 基準 1387–1414 K。
`notes/sessions/wall-temperature-three-way-analysis.md` §7)。SU2 は
`CNSSolver::SetTau_Wall_WF` で断熱壁温を Crocco–Busemann 回復温度へ更新する。

opt-in `turbulence.sstThermalWallFunction` (既定 0) で、代表点 (§6.5(a) の Normal_Neighbor
$\mathrm{rep}$) の状態から

$$T_{aw} = T_{\mathrm{rep}} + r\,\frac{U_{t,\mathrm{rep}}^2}{2 c_p}, \qquad r = \mathrm{Pr}^{1/3}$$

を毎ステップ評価する (`ransWallFunction_d.cu compute_wall_friction_sst_d`、出力配列 `Taw_diag`)。
**断熱壁 (`kind: wall`) かつ `wallTreatmentSST==1`** が対象、`wall_isothermal` には適用しない
(物理壁温が既にある)。

**モード体系 (2026-08-11 確定)**: `sstThermalWallFunction` = 0 (OFF) / 1 (output-only) /
2 (SU2 overlay, 実験・未採用) / **3 (defect-flux 閉包 — 正式採用, 断熱壁 node 生産閉包)**。
以下はまず mode 1 (output-only) の設計を述べる — mode 3 採用後は**比較用 fallback** の位置づけで、
$T_{aw}$ を**境界出力専用** — 「壁関数が再構成した物理壁面温度のモデル値」として `Taw_diag[W]` と
bvar `Tsb` (壁面出力 `res_wall_*` の `Ts`) に書くだけで、**場の解には一切作用させない**:

- $T_{aw}$ を使わないもの: $T[W]$ の状態ピン / GG・LSQ 勾配 (境界勾配は §6.2/§7.2.2 の
  owner-state-only により bvar を読まないので、`Tsb=T_{aw}` は構造的に勾配へ混入しない) /
  W–I 内部粘性流束 / `res_roe[W]` のゼロ化。
- 出力上の区別: **volume field の `VALUE/T` = 壁半 CV の計算状態温度 $T[W]$**、
  **`BCONDS/<wall>/VALUE/Ts` = モデル壁面温度 $T_{aw}$**。両者は別物として報告する
  (SU2 の壁関数壁温出力もこの意味のモデル値であり、モデル値同士の比較は正当だが、
  $T[W]$ と混同して「解いた場の壁温」と報告してはならない)。

**壁の熱経路 (正しい構造)**: 物理壁面 → (境界半割面流束 $q_w$) → 壁ノード $W$ →
(通常の内部粘性拡散: compact 温度差 $T[I]-T[W]$ + DOF 状態勾配) → 内部ノード $I$。
境界種別ごとの $q_w$: **断熱壁 = 厳密 $0$** (`viscousFlux_wall_d` の `adiabaticWall` 分岐) /
指定熱流束壁 = 指定値 (将来) / 熱壁関数 = 壁関数の $q_w$ (将来。`res_roe[W]` へ外部境界流束
として直接加え、内部辺のような $\pm F$ ペアにしない) / 等温 Dirichlet 壁 = $T[W]=T_w$ を強制し
$q_w$ は結果診断 (Dirichlet と Neumann を同時強制しない)。SST 断熱壁では $q_w=0$ なので
境界エネルギー流束はゼロで、壁半 CV の熱収支は W–I 内部流束と粘性仕事だけで通常どおり解かれる
(`res_roe[W]` はゼロ化しない)。

**却下した 2 方式 (実測発散, 復活させない)**:

1. **状態ピン版** ($T_{aw}$ を $T[W]$ にピン + `res_roe[W]` ゼロ化): 壁 CV 自身のエネルギー
   方程式を捨てながら W–I 流束だけ高い $T_{aw}$ で計算する非保存構成 (外部**熱源**) となり、
   $T_{rep}$→$T_{aw}$→強制コピーの正帰還で**過熱暴走** (node 1832 K / cell $T_t$ 飽和,
   run_0038/0039 旧版)。
2. **W–I compact 端点置換版** ($T[W]\to T_{aw}$ を内部辺の compact 温度差に注入。全辺・
   代表辺限定の両方): 置換後の流束が $T[W]$ に依存しなくなり**壁半 CV の復元項を喪失**、
   補償のない外部**吸熱**として壁ノードを EOS 床まで異常冷却して発散 (step1000 で
   $T[W]$ min 1100.7→27.1 K [全辺] / 15.5 K [代表辺限定], y+30 wall-function ケース,
   run_0053/0057)。「注入流束が有界」でも壁半 CV 単体の収支を閉じなければ発散する。

教訓: **W–I 内部辺の流束にモデル温度を混ぜる方式は、内部辺全体の $\pm F$ 保存だけでなく
壁半 CV 単体のエネルギー収支を必ず監査せよ**。壁への熱の出入りは物理壁面の境界半割面流束
$q_w$ だけで表現するのが正しい構造。

**旧「弱閉包」版の訂正 (2026-08-11)**: 旧版は「`Tsb=T_{aw}` は出力専用で場に触れない」と
していたが、当時は境界 GG/LSQ 勾配が bvar を読んでいたため勾配閉包経由で場に触れる経路が実在した
(効果は弱く、実データ [case/40 `run_0038`] では壁ノードの実状態 $T[W]$ [bell 平均 1195.3 K] は
OFF 基準 [1196 K] とほぼ不変 — 「1417.9 K = SU2 と 4 K 一致」の報告は**モデル値同士の一致**であり
「解いた場の壁温」ではなかった)。§6.2/§7.2.2 の owner-state-only 一般化後は bvar が勾配へ混入する
経路自体が消え、本設計 (出力専用) が構造的に保証されるようになった。経緯は
`plans/archived/turbulence-sst-thermal-wall-function.md` (初代・弱閉包) と
[`plans/accepted/turbulence-sst-adiabatic-taw-fluxmodel.md`](../../plans/accepted/turbulence-sst-adiabatic-taw-fluxmodel.md)
(流束注入の試行と棄却、最終方針 §0) を参照。

**mode 3 — 保存的壁層エネルギー流束閉包 (defect-flux, `sstThermalWallFunction: 3`, 2026-08-11 追加・
同日**正式採用** = 断熱壁 × node × `wallTreatmentSST==1` の生産閉包)**:
断熱壁の定応力 Couette 層ではエネルギー式 $\frac{d}{dy}(q+\tau u)=0$ を壁 ($q_w=0,\,u=0$) から
積分すると **任意の高さで $q(y)+\tau u(y)=0$** — 層が Crocco 平衡 ($T_W=T_{aw}$) に乗るとき
W–I 双対面を通る全 (伝導+粘性仕事) エネルギー流束は厳密ゼロになる。そこで W 入射 W–I 辺の
「線形伝導 + $\tau\!\cdot\!u$ 仕事」の合計を、平衡からのずれに比例する **defect 流束**

$$F^{wl}_e = H_T\,\bigl(T_{aw,rep}-T_W\bigr)\,(\vec S_e\!\cdot\!\hat m),\qquad
H_T=\frac{\rho_{rep}\,c_p\,u_\tau}{T^+(y^+_{rep};\mathrm{Pr})}$$

($\hat m$=内向き壁法線, $T^+$=Kader 原式 §10.3, $T_{aw,rep}$=Taw_diag と同式) に置換し、W/I へ
$\pm$ で加える。性質: (i) $T_W=T_{aw}$ で $F=0$、$\partial F/\partial T_W=-H_T A<0$ の**負帰還**で
$T_W\to T_{aw}$ が残差を解いた結果として実現 (ピン・ゼロ化・overlay 一切なし)、(ii) $\pm$ ペアで
厳密保存 (mode 2 の吸熱・ピン方式の湧き出しの両方を構造的に排除)、(iii) 双対 CV の面閉性から
$\sum_e \vec S_e\!\cdot\!\hat m = A_{wall}$ が厳密に成立し多重辺・非構造へ自然に一般化
(case/40 実測 $A_{eff}/A_{wall}=0.994$)、(iv) 平衡点は $T^+$ の値に依存しない ($T^+$ は緩和速度のみ)
— Kader 非圧縮較正の誤差が答えを汚さない、(v) $y^+\!\to\!0$ で $T^+\!\to\!\mathrm{Pr}\,y^+$ →
$H_T\to\lambda/y$ (解像 Fourier へ退化)。淀み ($u_\tau\to0$) は $H_T=\lambda_{rep}/y$ へ退避。
壁側 $\mu_t$ は熱経路から完全に外れる (mode 2 を殺した壁 $\mu_t$ 過大問題を回避)。運動量行
(AddTauWall)・物理壁面 $q_w=0$・GG/LSQ owner-state-only は不変。離散平衡は
$T_W^* = T_{aw} + R_{rest}/(H_T A)$ ($R_{rest}$=W–W' 辺+対流の残り) で、凍結場診断では
221 壁ノード中 212 で $|T_W^*-T_{aw}|<5\,$K。詳細・検証は
[`plans/accepted/turbulence-sst-su2-taw-coupling.md`](../../plans/accepted/turbulence-sst-su2-taw-coupling.md)。

**experimental mode 2 — SU2 式熱結合 (`sstThermalWallFunction: 2`, 2026-08-11 追加, 未採用)**:
モード体系は 0=無効 / 1=**output-only (生産 baseline, 上記の設計)** / 2=experimental SU2 coupled。
mode 2 は SU2 v8.5.0 の処理順 (勾配計算 → 壁 primitive 温度の $T_{aw}$ 上書き → 内部粘性流束が
corrected-gradient で $T_{aw}$ を端点参照) を、状態配列を汚さない per-node overlay 配列
`Taw_Prim_Overlay` で再現する: overlay 端点を持つ W 入射内部辺のみ、熱流束を SU2 式
$\vec g_{corr} = \bar{\vec g} + (\Delta T_{flux} - \bar{\vec g}\cdot\vec d)\vec d/|\vec d|^2$
(算術平均 $k_{eff}$)、粘性仕事の面速度を算術平均に切替え、壁ノードにも `roK_wf` Dirichlet を
拡張する (SU2 `SetTurbVars_WF` の両点書込)。GG/LSQ owner-state-only・$q_w=0$・
`res_roe[W]` 生存は不変。**検証状況: 凍結場収支診断・短尺動的試験の結果、現形は未採用**
(壁 μt を SU2 実効値へ寄せても壁ノードの単調冷却が止まらない。SU2 の安定は「primitive 上書きで
zombie 化した壁保存量の無害化」によるもので、forge は $T[W]$ が実効なため同機構が働かない)。
詳細と設計論点は
[`plans/accepted/turbulence-sst-su2-taw-coupling.md`](../../plans/accepted/turbulence-sst-su2-taw-coupling.md)
§10/§10b。

**cell の既知バイアス**: cell の代表点 = 第一セル (y+~30) の $T_1$ が node/SU2 の同高さより
~100–160 K 高く (cell wf=1 の BL 熱監査は follow-up)、$T_{aw}$ が +70–90 K 過大になる
(`run_0039`: 1490 K)。node を一次対象とする (ユーザ方針)。

#### (g) 等温壁のエネルギー壁関数 — Kader $q_w$ 流束置換 (`sstEnergyWallFunction`)

(f) の弱閉包は**断熱壁の壁温出力**を直すだけで、壁熱流束はモデル化しない。等温壁
(`kind: wall_isothermal`) × 壁関数メッシュ (第一 DOF が対数層) では、壁熱流束を解像伝導
$\lambda\,\partial T/\partial n$ で評価すると sublayer の温度勾配が解像されず $q_w$ を大きく
誤る ($\tau_w$ の解像勾配過小と同型の欠陥のエネルギー版)。冷却壁の熱負荷 $q_w$ が生産量の
ケース (冷却ノズル壁) では、壁隣接の伝導流束自体を Kader 型温度壁法則で**モデル置換**する。

opt-in `turbulence.sstEnergyWallFunction` (既定 0)。適用条件は `wallTreatmentSST==1` かつ
`wall_isothermal`。WMLES の熱壁モデル (§10.3–10.4) と同じ式・同じ流束置換機構を SST 側から
使う:

$$
q_w = \frac{\rho_{\mathrm{rep}}\,c_p\,u_\tau\,\big(T_w - T_{aw,\mathrm{rep}}\big)}{T^+(y^+;\ \mathrm{Pr})},
\qquad
T_{aw,\mathrm{rep}} = T_{\mathrm{rep}} + r\,\frac{U_{t,\mathrm{rep}}^2}{2 c_p},\quad r=\mathrm{Pr}^{1/3}
$$

- $T^+$ は Kader ブレンド (§10.3 の原式, `wallLaw_kader_tplus`)。$u_\tau,\ y^+$ は
  SST automatic wall treatment (§6.5(a)) の Reichardt Newton の収束値を再利用する
  (運動量経路はビット不変)。$\mathrm{Pr}=\mu\,c_p/\lambda$ は代表点で評価。
- 駆動温度差は回復温度 $T_{aw,\mathrm{rep}}$ ((f) の `Taw_diag` と同一式)。$T_m$ 直接では
  高速流の摩擦加熱分を誤る。$q_w>0$ = 壁→流体。
- **淀み/低 $y^+$ 退避**: $U_t\to0$ または $y^+<0.1$ では $T^+$ 評価が退化するため、純伝導式
  $q_w=\lambda\,(T_w-T_{aw,\mathrm{rep}})/y$ へ連続に退避する (WMLES §10.3 のガードと同処方。
  low-Re 極限では Kader の層流極限 $T^+\to\mathrm{Pr}\,y^+$ と一致し連続)。
- **流束の適用先**は WMLES §10.4 と同一機構: node は壁ノードに `Qw_Wall` を格納し
  W↔I 内部双対面の解像伝導を $q_w S$ に xor 置換 (AddQWall)。cell は壁境界面の熱流束を
  $q_w S$ に置換 (運動量側は §6.5 の従来処理のまま)。置換は**加算でなく置換**であること —
  1 面でも解像伝導が残ると (f) の正帰還の再発点になる。
- **状態層**は既存のまま: node 等温壁の壁ノード $T$ ピン + `res_roe` 0 化 + DPLUR roe 行
  decouple (`iso_wall_flag`) は乱流モデルによらず共有 (`nodeWallDirichlet` 系)。
- 断熱壁の $T_{aw}$ **強閉包** ($q_w=0$ 置換 + 壁状態 $T=T_{aw}$ 保持) は本機構の拡張として
  可能だが optional (弱閉包 (f) で壁温出力は足りている)。等温壁の検証が通ってから扱う。

**検証状況 (2026-08-11)**: 平板等温壁 (case/26, 真値 = y+0.35 low-Re 解像 [Colburn ±1.5%])
で現行 EWT の $q_w$ **+42% 過大** (y+30) が、本機構で **y+0.35: −0.1% / y+10: +6.2% /
y+30: +0.1%** に解消 (y+ 依存が実用上消滅)。**適用限界**: 非圧縮較正の Kader $T^+$ は強い
密度変化を持つ極超音速冷却壁 BL で破綻する — case/40 ベル部 (M≈4, $T_w/T_{aw}\approx0.4$)
で +37〜+265% (積分 +87%) の系統過大を実測 (チャンバ/スロートは −9% で実用域)。
$T^+$ の圧縮性補正 (semi-local / Crocco–Busemann / 壁物性評価) が済むまで、超音速ベルの
$q_w$ は y+≈1 low-Re を正とする。

計画: [`turbulence-sst-thermal-flux-model.md`](../../plans/active/turbulence-sst-thermal-flux-model.md)。

## 7. 2D / 3D / 軸対称

2D と 3D は、既存ソルバと同様に同一の離散化カーネルを共有する。
2D は単に $u_z = 0$ と z 方向幾何が退化した特別な 3D として扱う。

軸対称では事情が異なる。既存の軸対称実装は B 流儀

$$
V = \bar{r} A_{planar},
\qquad
\mathbf{S}_f = \bar{r}_f \, dl_f \, \hat{\mathbf{n}}
$$

を用いている。advection-only 段階では既存の体積・面積重みで scalar の対流を
通せるが、SST 特有の diffusion / source には幾何項の検討が必要である。
詳細は子 plan
[`architecture-axisym-sst.md`](../../plans/accepted/architecture-axisym-sst.md)
で扱い、以下に理論的な要点を整理する。

### 7.1 対流・拡散には幾何 source が残らない

$k$, $\omega$ の対流・拡散は円筒座標でも発散形を保つ。
$\phi \in \{k, \omega\}$ について

$$
\frac{1}{r}\partial_r\!\big(r\,\rho u_r \phi\big) + \partial_z(\rho u_z \phi)
=
\frac{1}{r}\partial_r\!\big(r\,\Gamma_\phi\,\partial_r \phi\big)
+ \partial_z(\Gamma_\phi\,\partial_z \phi)
+ S_\phi
$$

であり、両辺に $r$ を掛けて planar セル $A$ で積分すると、対流項・拡散項は
ともに $\partial_r(r\,\cdot) + \partial_z(r\,\cdot)$ の純粋な発散になる。
運動量方程式の $-\partial_r p$ のように「発散形に直せない残差」が無いため、
$r$ 重み付き面積 $\mathbf{S}_f = \bar r_f\,dl_f\,\hat{\mathbf n}$ を使う限り、
**対流・拡散からは軸対称特有の幾何 source は発生しない**。
拡散の $\theta\theta$ 寄与
$\frac{1}{r}\partial_r(r\,\Gamma_\phi\,\partial_r\phi)$ も $r$ 重み面積に
畳み込まれており、別途のソース項は不要である（子 plan で単体検証する）。

### 7.2 生産項のフープひずみ（幾何 source の本体）

幾何 source の本質は **ひずみ速度テンソルの周方向成分** にある。
軸対称流れ ($u_\theta = 0$, $\partial_\theta = 0$) でも、円筒座標のひずみ速度は
平面に現れない対角成分

$$
S_{\theta\theta} = \frac{u_r}{r} = \frac{U_y}{r}
$$

を持つ（フープひずみ, hoop strain）。SST の生産項は
$P_k = \mu_t S^2$, $P_\omega = \alpha\rho S^2$ で
$S = \sqrt{2\,S_{ij}S_{ij}}$ を用いるため、正しい大きさは

$$
S^2 = 2\!\left(S_{rr}^2 + S_{zz}^2 + S_{\theta\theta}^2\right)
+ \left(\partial_z u_r + \partial_r u_z\right)^2,
\qquad
S_{\theta\theta} = \frac{u_r}{r}
$$

となる。planar 速度勾配 ($\partial_x U_x \dots \partial_z U_z$) だけから組んだ
$S^2$ には $2\,S_{\theta\theta}^2 = 2\,(u_r/r)^2$ が欠落しており、軸近傍や急拡大部で
生産が過小評価される。同じ $S$ は渦粘性の制限子
$\mu_t = \rho a_1 k / \max(a_1\omega,\,S F_2)$ にも入るため、両方を一貫して
補正する必要がある。

### 7.3 圧縮性（dilatation）補正

超音速ノズルでは膨張・圧縮による発散 $\nabla\!\cdot\!\mathbf u$ が大きい。
軸対称の完全な発散は

$$
\nabla\!\cdot\!\mathbf u
= \partial_z u_z + \frac{1}{r}\partial_r(r u_r)
= \partial_x U_x + \partial_y U_y + \frac{u_r}{r}
\quad(=\texttt{axisym\_divU})
$$

であり、planar の発散に $u_r/r$ を加えた量となる。

#### 圧縮性 Boussinesq と生産項の正確形

圧縮性では乱流応力の Boussinesq 近似は等方成分を分離した形が正確である。

$$
\tau_{ij} = 2\mu_t\left(S_{ij} - \tfrac13(\nabla\!\cdot\!\mathbf u)\,\delta_{ij}\right)
- \tfrac23\rho k\,\delta_{ij}
$$

生産項 $P_k = \tau_{ij}\,\partial_j u_i = \tau_{ij} S_{ij}$ を展開すると

$$
P_k = \underbrace{2\mu_t\!\left(S_{ij}S_{ij} - \tfrac13(\nabla\!\cdot\!\mathbf u)^2\right)}_{\text{(A) deviatoric}}
\;\underbrace{-\;\tfrac23\rho k\,(\nabla\!\cdot\!\mathbf u)}_{\text{(B) isotropic}}
$$

となる。非圧縮形 $P_k = \mu_t S^2 = 2\mu_t S_{ij}S_{ij}$ (config `dilatationCorrection: 0`)
は (A) の $-\tfrac13(\nabla\!\cdot\!\mathbf u)^2$ と (B) を欠き、正確形との差は

$$
P_k^{\rm 現行} - P_k^{\rm 正確}
= \tfrac23\mu_t(\nabla\!\cdot\!\mathbf u)^2 + \tfrac23\rho k\,(\nabla\!\cdot\!\mathbf u).
$$

- **(A) deviatoric トレース除去**: 等方膨張・圧縮そのものはせん断生産に数えない。
  $S^2 \to S^2 - \tfrac23(\nabla\!\cdot\!\mathbf u)^2$（$S^2=2S_{ij}S_{ij}$ に対して
  $2\cdot\tfrac13(\nabla\!\cdot\!\mathbf u)^2 = \tfrac23(\nabla\!\cdot\!\mathbf u)^2$ を引く）。
  膨張コアでの spurious な $k$ 生成を抑える。低リスク。
- **(B) 等方項**: 膨張 ($\nabla\!\cdot\!\mathbf u>0$) で $k$ のシンク、
  圧縮 ($\nabla\!\cdot\!\mathbf u<0$, 衝撃波) でソース。$k$ 方程式固有で、
  $\omega$ 生産には入れない。膨張で $P_k<0$ になり得るため $P_k \ge 0$ クリップを併用する。

これら (A)(B) はいずれも $\nabla\!\cdot\!\mathbf u$ (= `axisym_divU`) を用いる。
軸対称では $u_r/r$ を含む完全発散を使わなければ強圧縮場で過大/過小評価が残る。

#### 別物: dilatation-dissipation 系（採用しない）

Sarkar / Zeman / Wilcox の乱流マッハ数補正
$\varepsilon = \varepsilon_s(1 + \alpha_1 M_t^2)$, $M_t=\sqrt{2k}/a$ は
**散逸**を増やす別系統の補正で、自由せん断層向けにキャリブレートされている。
壁境界層・ノズル内部では逆効果になりやすく、標準 SST には組み込まない。
本実装では (A)(B) の生産項補正のみを対象とし、$M_t$ 散逸補正は採用しない。

### 7.4 既存の軸対称診断量の再利用

7.2・7.3 で必要な量は、NS 用に
[`axisymmetricSource_d.cu`](../../solver_density_cuda/cuda_forge/axisymmetricSource_d.cu)
が既に計算・保持している。

| 量 | 変数 | 内容 |
| --- | --- | --- |
| フープひずみ $u_r/r$ | `axisym_uy_over_r` | $S_{\theta\theta}$ そのもの |
| 完全発散 $\nabla\!\cdot\!\mathbf u$ | `axisym_divU` | $\partial_x U_x+\partial_y U_y+u_r/r$ |

したがって SST 側は新規物理量を計算せず、これらを kernel 引数で受け取り、
`isAxisymmetric` のとき $S^2$ と発散補正に畳み込むだけでよい。

### 7.5 応力生産アノマリーと Kato–Launder 補正

標準 SST の生産項 $P_k = \mu_t S^2$ ($S=\sqrt{2S_{ij}S_{ij}}$) は **ひずみの大きさ**のみで
生産を見積もる。ところがひずみテンソルは

$$
S_{ij} = \underbrace{\tfrac13(\nabla\!\cdot\!\mathbf u)\delta_{ij}}_{\text{等方(膨張)}}
+ \underbrace{S'_{ij}}_{\text{偏差}},\qquad
S'_{ij} = \underbrace{\text{せん断成分}}_{\Omega\ \text{を伴う}}
+ \underbrace{\text{法線(伸長)異方性}}_{\text{非回転}}
$$

と分解され、偏差ひずみ $S'$ は**せん断**だけでなく**非回転の伸長(加速)ひずみ**も含む。
速度勾配の反対称部=渦度 $\Omega_{ij}=\tfrac12(\partial_j u_i-\partial_i u_j)$,
$\Omega=\sqrt{2\Omega_{ij}\Omega_{ij}}=|\boldsymbol\omega|$ は回転の不変尺度であり、
**単純せん断では $S\approx\Omega$、純伸長(よどみ点・加速)では $\Omega\to0$ かつ $S$ 大**となる。

§7.3 の `dilatationCorrection` は等方(トレース=$\nabla\!\cdot\!\mathbf u$)分のみを除くため、
**非回転の偏差伸長ひずみは残り**、強加速場（ノズル喉中心線・よどみ点）で
$P_k=\mu_t S'^2$ が偽の乱流を生む。これは strain-based SST の古典的欠陥
(**stagnation / round-jet anomaly**) であり、`dilatationCorrection` では消えない別系統の問題。

**Kato–Launder 修正** (Kato & Launder 1993) は生産の片方の $S$ を渦度 $\Omega$ に置換する:

$$
P_k = \mu_t\,S\,\Omega \qquad (\text{標準は } \mu_t S^2)
$$

- せん断層 ($S\approx\Omega$): $P_k\approx\mu_t S^2$ で**従来と一致**(壁 BL の本物の乱流は不変)。
- 非回転加速 ($\Omega\to0$): $P_k\to0$ で**偽生産を除去**。

不変量 $S,\Omega$ を使うため**座標系に依存しない**(「対角項除外」のような非不変操作は不可:
45° 回転でせん断と伸長が入れ替わるため)。

#### case 29 の「軸中心 k 過大」— 機構特定済み: **二層構造** (経緯の誤診は訂正済み)

詳細は [`architecture-axisym-axis-singularity.md`](../../plans/accepted/architecture-axisym-axis-singularity.md) /
[case 29 検証](../../case/29.bell_vs_conical/README.md)。結論:

**層 (a) — 実用格子 (軸セル ≥ 0.25mm) の軽度スパイク = 生産駆動。Kato–Launder で解決。**
平均流は滑らか。軸第1セル $k$=17 vs 核 4.5 は、大きな軸方向ひずみ + フープ + 生産リミタ正帰還 +
軸無フラックス滞留による。`katoLaunder: 1` で $k$ 17→1.9、壁 BL・推力不変 — 3D の平坦
プロファイル (第1セル/核 ~0.85) と整合する妥当な対処。

**層 (b) — 軸を ≤10µm に細分すると平均流の数値不安定 (未根治)。**
2 面クラスタ格子 (壁解像維持 + 軸細分) で **$U_y$ が近軸チェッカーボード (偶奇) 振動**
(物理 ~0.5 m/s のところ ±7 m/s 交番)。勾配計算は場の FD と一致 (正しい) — 場そのものが振動し、
$(\partial_r U_y)^2 \sim 10^{11}$ が $S^2$ 経由で $k$ を爆発させる。軸細分で悪化 (k 収束せず /
rms_roK 1e11 爆発)。機構: 半径方向は近軸で深い低 Mach ($u_r\to0$) → 風上流束の圧力-速度結合が
消え偶奇モードが減衰しない (古典的チェッカーボード)。3D は軸が内部点でこの構造を持たない。
`lowMachPrecond` 1/2 とも失敗 (超音速コア/細軸で発散)。

**実用指針**: 軸対称 RANS ノズルは `katoLaunder: 1` を使い、**軸方向の過細分を避ける**
(軸セル ≥ ~0.1mm; 壁解像は壁側クラスタで)。

> 経緯メモ (誤診の訂正): ①「planar 面積 vs $r$ 重み体積の不整合」は誤り (flux 面積・体積は両方
> $r$ 重み済, [`variables.cpp`](../../solver_density_cuda/variables.cpp))。②「勾配の planar 面積が
> bug」も誤り (勾配は planar 作用素で正しい)。③ uniform 格子のグリッド試験は $k$ 非収束で無効
> だった。最終特定は 2 面クラスタ格子での $U_y$ チェッカーボード観測による。

#### 適用範囲・既定

- config `katoLaunder` (turbulence): `0`=標準 $\mu_t S^2$ (既定, 既存挙動とビット一致),
  `1`=Kato–Launder $\mu_t S\Omega$。
- $\omega$ 生産 $P_\omega=\alpha\rho S^2$ にも同じ置換 ($\alpha\rho S\Omega$) を適用し、$\mu_t$ の
  時間スケール整合を保つ。
- `dilatationCorrection` と**直交・併用可** (等方分は除いた上で $S\to\sqrt{S^2_{\rm corr}}$ と
  $\Omega$ の積を取る)。

## 8. DES / DDES (Detached-Eddy Simulation)

剥離支配の非定常流 (SBLI・ピントルバルブ・後退段) を RANS より高精度・wall-resolved LES より
低コストで解くため、SST に **DDES (Delayed DES, Spalart et al. 2006)** の length scale 修正を加える。
付着境界層は RANS のまま保護し、剥離せん断層・自由せん断層でのみ LES へ切り替える。
実装は [`implementation.md` §3.8](implementation.md) と plan
[`turbulence-iddes-sst.md`](../../plans/active/turbulence-iddes-sst.md) を参照。

### 8.1 RANS length scale と k 消滅項の整合

SST の k 方程式の消滅項は
$$
D_k = \beta^{*}\rho k\omega = \frac{\rho k^{3/2}}{l_\mathrm{RANS}}
$$
である。この恒等式から RANS length scale は
$$
l_\mathrm{RANS} = \frac{\sqrt{k}}{\beta^{*}\omega}\qquad(\beta^{*}=0.09)
$$
**でなければならない**。$\rho k^{3/2}/l_\mathrm{RANS} = \rho k^{3/2}\cdot\beta^{*}\omega/\sqrt{k} = \beta^{*}\rho k\omega$
と一致する。これは Strelets (2001) / Gritskevich (2012) の SST-DES で用いられる定義であり、
DES は **この $l_\mathrm{RANS}$ を $l_\mathrm{DDES}$ に置換するだけ**で、$f_d=0$ の完全シールド域では
$l_\mathrm{DDES}=l_\mathrm{RANS}$ となり $D_k$ が標準 SST に厳密に戻る。

> **注意 (実装で判明した落とし穴)**: $l_\mathrm{RANS}=\sqrt{k}/(\beta^{*1/4}\omega)$ と書くと上の恒等式を
> 満たさず、$D_k=\beta^{*1/4}\rho k\omega$ となって標準 SST より係数が $\beta^{*-3/4}\approx 6.1$ 倍
> 大きい消滅項になる。この場合 $f_d=0$ の付着 BL でも mode1 が標準 SST に縮退せず、
> **モデル応力枯渇 (Modeled Stress Depletion)** で BL の $k$ が崩壊する (検証で `roK` が ~100% 変化)。
> 必ず $\beta^{*}$ (=0.09) を使うこと。

### 8.2 LES グリッドスケール $\Delta_\mathrm{max}$

$$
\Delta = \Delta_\mathrm{max} = \max_{\text{隣接面}} |\mathbf{cc}_{ic} - \mathbf{cc}_{jc}|
$$

隣接 CV 重心間距離の最大値。等方メッシュなら $V^{1/3}$ で代用できるが、**境界層の高アスペクト比
セル**では $V^{1/3}$ が接線スケールを大幅に過小評価し、シールドがあっても DES リミッタが
BL 内で誤発火しうる。$\Delta_\mathrm{max}$ は Spalart 系で標準。幾何セットアップ時に実行時の
`volume`/`ccx` (= CV 重心) から 1 回だけ計算するので、cell / node (median-dual) 双方で自動的に
双対 CV の $\Delta$ になる。

### 8.3 DDES シールド関数 (Spalart et al. 2006)

$$
r_d = \frac{\nu_t + \nu}{\kappa^2 d^2\,\sqrt{\partial U_i/\partial x_j\,\partial U_i/\partial x_j}}\qquad(\kappa=0.41)
$$
$$
f_d = 1 - \tanh\!\bigl([8\,r_d]^3\bigr)
$$

| $r_d$ | 物理的意味 | $f_d$ | 動作 |
|-------|-----------|-------|------|
| $\gtrsim 1$ | BL 内部 (渦粘性大 / 近壁 $d^2$ 小) | $\approx 0$ | 純 RANS (保護) |
| $\to 0$ | 剥離域・自由せん断層 | $\to 1$ | DES limiter 有効 |

分子の $(\nu_t+\nu)$ により、対数層では $\nu_t\!\sim\!\kappa u_\tau y$, $S\!\sim\!u_\tau/(\kappa y)$ から
$r_d\!\to\!1$ となり付着 BL が保護される。粘性低層では $\nu$ 項により $r_d\!\sim\!1/(\kappa^2 y^{+2})$ で
やはり保護される。

### 8.4 ハイブリッド length scale

$$
l_\mathrm{DDES} = l_\mathrm{RANS} - f_d\cdot\max\!\bigl(0,\; l_\mathrm{RANS} - C_\mathrm{DES}\Delta\bigr)
$$
$$
D_k^\mathrm{DDES} = \frac{\rho k^{3/2}}{l_\mathrm{DDES}},\qquad C_\mathrm{DES}=0.78\ (\text{初版 }k\text{-}\omega\text{ 固定})
$$

リミッタは **2 重に保護**される: ① シールド $f_d\approx 0$、または ② $l_\mathrm{RANS}<C_\mathrm{DES}\Delta$
(グリッドが BL 長スケールに対して粗い near-wall) で $\max(0,\cdot)=0$ となり、どちらでも
$l_\mathrm{DDES}=l_\mathrm{RANS}$ に戻る。これにより付着 BL は $f_d$ が多少立っても撹乱されない。
$\omega$ 方程式・渦粘性式は **変更しない** (標準 SST-DES は $k$ 消滅項のみ修正)。

$C_\mathrm{DES}$ は本来 $F_1$ ブレンド $C_\mathrm{DES}=F_1 C_{\mathrm{DES},k\omega}+(1-F_1)C_{\mathrm{DES},k\varepsilon}$
($0.78/0.61$) だが、渦粘性カーネルが $F_1$ を持たないため初版は $C_{\mathrm{DES},k\omega}=0.78$ 固定とする。

### 8.5 float32 ガード

forge は float32。$r_d$ は $d^2$ で効くため近壁・よどみ・一様流で飽和しやすい。
勾配大きさ・分母に floor を入れて 0 割を防ぎ ($\nabla u\to 0$ の一様流で $r_d\to$ 大 $\to f_d\to 0$ =
RANS が正しい極限)、$r_d$ を $[0,10]$ に clamp して $f_d$ の 0/1 張り付きを抑える。
診断のため `r_d`・`f_d`・`l_des`・`delta_les` を HDF5 出力する (一様な $f_d$ 張り付きは NaN より
発見しにくい DES 不全の兆候)。

### 8.6 IDDES (`model: "sst-iddes"`, Shur et al. 2008 / SST 版 Gritskevich et al. 2012)

DDES に **WMLES 分岐** を加えた拡張。解像乱流が存在する近壁領域では、RANS 層を壁モデルとして
薄く残し ($f_B$ で壁直近のみ RANS 強制)、log 層から LES に切り替える。解像乱流が無い場合は
DDES 同等に縮退する。実装は Gritskevich et al. (2012) の SST-IDDES に従う
(関数形状の数値検証: `solver_density_cuda/tools/verify_iddes_functions.py`)。

**壁距離依存グリッドスケール** (Shur 原型の壁法線刻み $h_{wn}$ 項を落とした Gritskevich 簡略形。
非構造格子で壁法線刻みの定義が不要で、$d_w$ = `wall_dist` と $h_\mathrm{max}$ = `delta_les` のみで組める):

$$
\Delta_\mathrm{IDDES} = \min\bigl(C_w \max(d_w, h_\mathrm{max}),\; h_\mathrm{max}\bigr),\qquad C_w=0.15
$$

壁で $0.15h_\mathrm{max}$ → 遠方で $h_\mathrm{max}$ ($[0.15h,h]$ に有界・単調)。近壁で $\Delta$ が
縮むことで $l_\mathrm{LES}$ が下がり、WMLES 分岐が log 層まで LES を入れられる。

**関数群** ($\alpha = 0.25 - d_w/h_\mathrm{max}$、$r_{dt}$: $\nu_t$ 版・$r_{dl}$: $\nu$ 版の $r_d$):

$$
f_B = \min(2e^{-9\alpha^2}, 1),\qquad
f_{e1} = \begin{cases}2e^{-11.09\alpha^2} & \alpha\ge0\\ 2e^{-9\alpha^2} & \alpha<0\end{cases}
$$
$$
f_t = \tanh[(C_t^2 r_{dt})^3],\quad f_l = \tanh[(C_l^2 r_{dl})^{10}],\quad
f_{e2} = 1-\max(f_t,f_l),\quad f_e = f_{e2}\max(f_{e1}-1,0)
$$
$$
f_{dt} = 1-\tanh[(C_{dt1}r_{dt})^3],\qquad \tilde f_d = \max(1-f_{dt},\,f_B)
$$
$$
l_\mathrm{IDDES} = \tilde f_d(1+f_e)\,l_\mathrm{RANS} + (1-\tilde f_d)\,C_\mathrm{DES}\Delta_\mathrm{IDDES}
$$

定数 (SST 較正値): $C_t=1.87$, $C_l=5.0$, $C_{dt1}=20$, $C_w=0.15$。$f_e$ は log-layer mismatch
(RANS→LES 遷移部の速度対数則ずれ) を $l_\mathrm{RANS}$ の増強で補正する項で、付着 log 層
($r_{dt}\approx1$) では $f_t\to1\Rightarrow f_e=0$ と自動停止する。

挙動 (verify スクリプトで数値確認済):

| 状態 | $f_{dt}$ | $\tilde f_d$ | $l_\mathrm{IDDES}$ |
|------|---------|--------------|--------------------|
| 付着 BL (モデル乱流 $r_{dt}\approx1$) | $\to0$ | $\to1$ ($f_e=0$) | $= l_\mathrm{RANS}$ 厳密 (シールド) |
| 解像乱流あり (瞬時 $\nu_t$ SGS 級 $r_{dt}\ll0.02$) | $\to1$ | $\to f_B$ | WMLES 分岐 ($f_B(1+f_e)l_\mathrm{RANS}+(1-f_B)l_\mathrm{LES}$) |
| 乱流流入なし | $\to0$ 側 | DDES 同等 | DDES 縮退 |

$D_k$ の切替 (§8.1)・$\omega$ 方程式不変・float32 ガード (§8.5) は DDES と共通。DDES の
$f_d$ ($C_{d1}=8$, Spalart 2006 較正) は変更せず、IDDES の $f_{dt}$ は $C_{dt1}=20$ の別関数として持つ。
WMLES モードの**定量**評価は解像乱流の存在が前提 (周期チャネルは体積力 config、非周期ケースは
流入乱流生成にゲート。機能確認 = $\tilde f_d/f_e/l_\mathrm{des}$ 分布の妥当性までが現フェーズ)。

---

## 9. 初期実装の範囲外

次は本章の対象外とする。

- 壁関数
- 遷移モデル
- 圧縮性補正の詳細オプション化
- 陰解法 Jacobian への SST 連成

軸対称 SST の幾何 source (§7.1–7.4) は子 plan
[`architecture-axisym-sst.md`](../../plans/accepted/architecture-axisym-sst.md)
で実装・検証する。その他は別 plan または後続フェーズで扱う。
---

## 10. WMLES 代数壁応力モデル (Reichardt + Kader)

計画: [`turbulence-wmles-wall-stress.md`](../../plans/active/turbulence-wmles-wall-stress.md)。
壁モデル LES (WMLES) では第一解点を対数層 ($y^+\approx 50\text{–}100$) に置き、壁面せん断応力
$\tau_w$ と壁面熱流束 $q_w$ を代数壁法則で与えて壁の粘性流束を置き換える。SST automatic wall
treatment (§6.5) と同じ Reichardt 則を速度に、Kader 則を温度に使う平衡壁モデルである。
LES (`model: wale/sigma`) 系の機能であり、SST の $\omega$ ピン・$P_k$ 置換 (§6.5(b)(d)) は関与しない。

### 10.1 マッチング点と壁面状態

- マッチング点 = 壁に最も近い内側の解点。cell: 壁第一セル中心 ($y=$ `wall_dist`)。
  node: SU2 Normal_Neighbor 流に壁ノードの内部双対面から内向き法線 cos 最大の内部ノード
  ($y=(\mathbf x_I-\mathbf x_W)\cdot(-\hat{\mathbf n})$)。§6.5(a) の代表点選択と同一機構。
- 壁面状態 $T_w, p_w$: node は壁ノードの解、cell は ghost との面平均 (断熱ミラーでは内点値に一致)。
  等温壁は指定 $T_\mathrm{wall}$ を優先する。
- 壁面物性: $\rho_w = p_w/(R_w T_w)$、$\mu_w=\mu(T_w)$、$\lambda_w=\lambda(T_w)$、
  $c_{p,w}=c_p(T_w)$。`viscMethod` (定数 / Sutherland / 分子論) と `thermalMethod` (CPG/TP) に
  整合した式で評価し、$\gamma$ 一定・$c_p$ 一定を焼き込まない。
  $\mathrm{Pr}=\mu_w c_{p,w}/\lambda_w$ もその場評価する (層流 Pr の設定キーは存在しない)。

### 10.2 運動量 — Reichardt 則の逆解き (warm start Newton)

壁平行速度 $u_{\parallel,i}=u_{i,m}-(u_{j,m}n_j)n_i$ に対し、Reichardt 則 (§6.5(a) と同一式)
$u^+=f(y^+)$ を

$$
F(u_\tau) = u_\tau f\!\big(y^+(u_\tau)\big) - u_\parallel = 0 ,\qquad
F'(u_\tau) = f(y^+) + y^+ f'(y^+)
$$

の形で Newton 反復して $u_\tau$ を求める。SST 壁関数の残差形 $U_t/u_\tau - u^+$ (§6.5(a)) と根は
同一だが、この形は $u_\tau\to 0$ で正則であり warm start (前 step の面ごと $u_\tau$ を初期値に採用)
と相性が良い。初回・無効値は粘性則 $u_\tau^{(0)}=\sqrt{\mu_w u_\parallel/(\rho_w d)}$。
収束判定は相対 $10^{-6}$・最大 20 回。不収束時は層流応力 $\tau_w=\mu_w u_\parallel/d$ に
フォールバックしカウンタを残す。$u_\parallel$ が機械精度近傍なら Newton をスキップして層流極限。
出力は $\tau_w=\rho_w u_\tau^2$、向きはマッチング点瞬時速度 $\hat e_{\parallel}$ (cell)。

### 10.3 エネルギー — Kader 温度壁法則と回復温度

等温壁のみ計算する (断熱壁は $q_w=0$ で閉じ、壁温は解に含まれる摩擦加熱込みの値)。

$$
T^+ = \mathrm{Pr}\,y^+ e^{-\Gamma} + \left[2.12\,\ln(1+y^+) + \beta(\mathrm{Pr})\right] e^{-1/\Gamma},
\qquad
\Gamma = \frac{0.01\,(\mathrm{Pr}\,y^+)^4}{1 + 5\,\mathrm{Pr}^3 y^+}
$$

$$
\beta(\mathrm{Pr}) = \left(3.85\,\mathrm{Pr}^{1/3} - 1.3\right)^2 + 2.12\ln \mathrm{Pr}
$$

(Kader 1981 の原式。**旧版 (〜2026-08-11) は対数部を $\mathrm{Pr}_t(u^+ + \beta)$ としていた** —
Jayatilleke 型の運動量アナロジー形と Kader の $\beta$ ($2.12\ln(1+y^+)$ と対で較正された定数)
の混用で、log 層の $T^+$ を +30% 級過大評価し等温壁 $q_w$ が系統的に過小になっていた
(case/26 平板 y+30 で −13%)。原式に修正。$u^+$・$\mathrm{Pr}_t$ は $T^+$ に不要になった。)

駆動温度差は回復温度 $T_r = T_m + r\,u_\parallel^2/(2c_{p,w})$ ($r=\mathrm{Pr}^{1/3}$、
$c_p$ ベースで $\gamma$/Mach を使わない) を用い

$$
q_w = \frac{\rho_w\,c_{p,w}\,u_\tau\,(T_r - T_w)}{T^+}
$$

($q_w>0$ が壁→流体)。$T_m$ を直接使うと高速流で摩擦加熱分を誤るため $T_r$ を使う。
$y^+, u^+$ は §10.2 で収束した値を再利用する。

### 10.4 流束の適用先 (cell / node)

- **cell**: 壁境界面の粘性流束を置換する。運動量は接線せん断を $\tau_w$ (向き $-\hat e_\parallel$)
  に、エネルギー熱流束を $q_w S$ に置き換える。壁面 no-slip のため粘性仕事は 0。法線粘性・
  体積項は落とす (§6.5(c) の SST 分岐と同じ処方)。
- **node**: 壁ノードは速度 Dirichlet で運動量残差がゼロ化されるため、壁境界半割面の置換では
  内点に応力が届かない。SST node と同じく壁ノードに $\tau_w=\rho_w u_\tau^2$ を格納し、
  W↔I 内部双対面の解像 traction を大きさ $\tau_w$ に再スケールする (AddTauWall, §6.5(e))。
  向きは面ごとの解像接線方向 (瞬時場に追随)。エネルギーも対称に、W↔I 面の解像伝導熱流束を
  $q_w$ ベースに再スケールする (AddQWall)。等温壁では壁ノード温度を $T_\mathrm{wall}$ に
  Dirichlet ピンする。
- SGS 渦粘性は壁面流束評価には使わない (壁モデルが全応力を与える)。内部領域の SGS・対流
  (KEEP+ES 散逸)・勾配計算は不変。

### 10.5 適用限界

平衡壁法則 (§10.2–10.3) は付着・準平衡境界層が前提。強い非平衡 (衝撃波直下・剥離点近傍) では
$\tau_w$ を誤り、高 Re の緩剥離では標準 SGS 係数との組合せで剥離を見逃す実例がある
([`turbulence-des-wmles-survey.md`](../../notes/investigations/turbulence-des-wmles-survey.md))。
剥離を含むケースでは SGS 係数・ES 散逸 $\sigma$ の感度確認を必須とする。

## 11. 遷移モデル $\gamma$–$Re_{\theta t}$ (Langtry–Menter 2009, `turbulence.transition: "lm2009"`)

低 Re SST は遷移の閉包を持たず、前縁から乱流になる。冷却翼 (case/53, 54) の負圧面層流域で $h$ が +44〜+80 % 過大に出ることの**最優先の仮説**が
これである (層流対照はその領域で実測を下回り、同一メッシュの SU2 も同じだけ過大)。そこで SST に Langtry–Menter の 2 方程式遷移モデル
(AIAA J. 47(12), 2009。相関は同論文の $Re_{\theta c}$・$F_{length}$ = SU2 の `LM` 既定 `MENTER_LANGTRY`) を足す。
**式と定数は SU2 8.x の実装 (`trans_sources.hpp`, `trans_correlations.hpp`, `CTransLMSolver.cpp`) と 1 対 1 に合わせる**
(同一メッシュでの突き合わせを検証の柱にするため)。

### 11.1 輸送方程式

間欠度 $\gamma$ と遷移開始運動量厚さレイノルズ数 $\tilde{Re}_{\theta t}$ を保存形 $\rho\gamma$, $\rho\tilde{Re}_{\theta t}$ で運ぶ。

$$\frac{\partial\rho\gamma}{\partial t}+\nabla\cdot(\rho u\gamma)=P_\gamma-E_\gamma+\nabla\cdot\Big[\Big(\mu+\frac{\mu_t}{\sigma_f}\Big)\nabla\gamma\Big],\qquad
\frac{\partial\rho\tilde{Re}_{\theta t}}{\partial t}+\nabla\cdot(\rho u\tilde{Re}_{\theta t})=P_{\theta t}+\nabla\cdot\big[\sigma_{\theta t}(\mu+\mu_t)\nabla\tilde{Re}_{\theta t}\big]$$

$$P_\gamma=F_{length}\,c_{a1}\,\rho S\,\sqrt{\gamma F_{onset}}\,(1-c_{e1}\gamma),\qquad
E_\gamma=c_{a2}\,\rho\Omega\,\gamma F_{turb}\,(c_{e2}\gamma-1),\qquad
P_{\theta t}=c_{\theta t}\frac{\rho}{t}\,(Re_{\theta t}-\tilde{Re}_{\theta t})(1-F_{\theta t})$$

定数: $c_{e1}=1$, $c_{a1}=2$, $c_{e2}=50$, $c_{a2}=0.06$, $\sigma_f=1$, $c_{\theta t}=0.03$, $\sigma_{\theta t}=2$。$t=500\mu/(\rho U^2)$。

**開始関数**: $Re_v=\rho y^2S/\mu$, $R_T=\rho k/(\mu\omega)$,
$F_{onset1}=Re_v/(2.193\,Re_{\theta c})$, $F_{onset2}=\min(\max(F_{onset1},F_{onset1}^4),2)$, $F_{onset3}=\max(1-(R_T/2.5)^3,0)$,
$F_{onset}=\max(F_{onset2}-F_{onset3},0)$, $F_{turb}=e^{-(R_T/4)^4}$。
$F_{length}=F_{length,1}(1-F_{sub})+40F_{sub}$, $F_{sub}=e^{-(R_\omega/200)^2}$, $R_\omega=\rho y^2\omega/\mu$。

**境界層の外で $\tilde{Re}_{\theta t}$ を自由流の相関値へ引く**ための遮蔽:
$\theta_{BL}=\tilde{Re}_{\theta t}\mu/(\rho U)$, $\delta_{BL}=7.5\theta_{BL}$, $\delta=50\,\Omega y\,\delta_{BL}/U$,
$F_{wake}=e^{-(Re_\omega/10^5)^2}$ ($Re_\omega=\rho\omega y^2/\mu$),
$F_{\theta t}=\min\{\max[F_{wake}e^{-(y/\delta)^4},\,1-((\gamma-1/c_{e2})/(1-1/c_{e2}))^2],\,1\}$。

**経験相関** ($Tu$ [%] $=100\sqrt{2k/3}/U$、下限 0.027):
$Re_{\theta t}=(1173.51-589.428Tu+0.2196/Tu^2)F(\lambda_\theta)$ ($Tu\le1.3$)、$331.5\,(Tu-0.5658)^{-0.671}F(\lambda_\theta)$ ($Tu>1.3$)、下限 20。
$\lambda_\theta=(\rho\theta^2/\mu)\,dU/ds$ を $[-0.1,0.1]$ に制限、$\theta=Re_{\theta t}\mu/(\rho U)$ なので**不動点反復**で解く。
$F(\lambda_\theta)=1-(-12.986\lambda-123.66\lambda^2-405.689\lambda^3)e^{-(Tu/1.5)^{1.5}}$ ($\lambda\le0$)、
$1+0.275(1-e^{-35\lambda})e^{-Tu/0.5}$ ($\lambda>0$)。$dU/ds$ は速度勾配テンソルから $\hat u\cdot\nabla|u|$。
$Re_{\theta c}(\tilde{Re}_{\theta t})$ と $F_{length,1}(\tilde{Re}_{\theta t})$ は Langtry–Menter の区分多項式 (SU2 `MENTER_LANGTRY`)。

### 11.2 SST との結合

剥離誘起遷移の補正 $\gamma_{sep}=\min\{2\max[Re_v/(3.235Re_{\theta c})-1,0]F_{reattach},2\}F_{\theta t}$, $F_{reattach}=e^{-(R_T/20)^4}$、
$\gamma_{eff}=\max(\gamma,\gamma_{sep})$ を使って、$k$ 式だけを変える ($\omega$ 式は不変):

$$\tilde P_k=\gamma_{eff}P_k,\qquad \tilde D_k=\min(\max(\gamma_{eff},0.1),1)\,D_k,\qquad F_1=\max(F_{1,orig},F_3),\ F_3=e^{-(R_y/120)^8},\ R_y=\rho y\sqrt k/\mu$$

### 11.3 境界条件と自由流

壁: $\gamma$, $\tilde{Re}_{\theta t}$ とも**法線勾配 0** (node では壁節点を固定せず、壁半割面の拡散流束を足さない)。
入口: $\gamma=1$、$\tilde{Re}_{\theta t}$ = 入口の $Tu$ から上の相関 ($\lambda_\theta=0$) で決めた値。出口・slip・周期: $k$/$\omega$ と同じ扱い。
**自由流の $Tu$ の減衰が遷移位置を決める**ので、入口の $k$ と $\omega$ (= 粘性比) の与え方は結果に直結する
($k$ の減衰率は $\omega$ で決まる)。翼列では SST だけのときに既に「粘性比が圧力面の $h$ を 23 pt 動かす」ことが分かっており
(報告 §02)、遷移モデルではこの感度がさらに大きくなる。
