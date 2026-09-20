# case/51 — TN D-8233 深いすきま (乱流, M 10.3) の厳密 validation

**T2**。[Throckmorton 1976, NASA TN D-8233] の単一横断すきま (W = 2.29 mm, D = 45.72 mm,
**D/W = 20.0**) を M 10.3 の乱流境界層下で再現する。T3 (case/52) が文献帯の内側に収まらなかった
ため優先度を上げた。T3 との違いは、**この報告だけが h = q/(T_aw − T_w) の形で分母を定義し、
壁温を振り、BL 厚を測っている**こと。

計画: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md)
調査: [`notes/investigations/hypersonic-gap-cavity-heating-survey.md`](../../notes/investigations/hypersonic-gap-cavity-heating-survey.md)

## この case が持つ非自明な入力

`conditions.json` が入力台帳 (本文・Appendix・Fig 4 / Fig 8(a) を 200–260 dpi 描画から目視で読み取り)。

| 量 | 値 | 出所 |
| --- | --- | --- |
| M | 10.3 | 本文 |
| すきま | W 2.29 mm × D 45.72 mm (D/W 19.97) | 本文 |
| $P_t$ | 2.41 / 5.34 / 12.07 MN/m² | 本文 |
| $\delta^*/W$ | **46–56** (δ* 10.5–12.8 cm) | Fig 4 |
| $h_{fp}$ | $7.95\times10^{-5}\,Re'^{0.69}$ W/(m²K) | Appendix (A2) + 図 4 点 |
| 壁温比 | **Re 系列ごとに実測値が違う** (1.47e6: 1.00/1.19/1.37/1.61 …) | Fig 8(a) |
| $T_t$ | **本文に記載なし** → Gate A で解いた | ↓ |

## Gate A — 本文に無い $T_t$ を解いて決めた (**PASS**)

$P_t$・$M$・$Re'$ が 3 系列とも与えられているので $T_t$ は未知 1・条件 1 で解ける。
半完全気体で厳密に解く ($M$ 10.3 では $T_\infty \approx 51$ K まで落ちるため $\gamma$ 一定は使えない)。

```
  Re_inf [1/m]  Pt [MPa]   Tt [K]  T_inf [K]  p_inf [Pa]  rho [kg/m3]   U [m/s]  T_aw [K]
     1.470e+06      2.41   1112.0      52.63        57.0      0.00377    1498.7     995.5
     3.320e+06      5.34   1096.1      51.81       125.6      0.00845    1486.9     981.3
     7.820e+06     12.07   1062.6      50.08       281.0      0.01955    1461.8     951.3

解けた Tt : 1112.0 / 1096.1 / 1062.6 K   平均 1090.3 K, ばらつき 4.5 %
GATE A VERDICT: PASS  (3 系列が同じ Tt を返すか; 閾値 10 %)
```

独立に digitize した 3 組の $(P_t, Re')$ が同じ $T_t$ を返した。これは
**(a) 欠けていた $T_t$ が 1090 K と決まった** (Langley CFHT の公称運転温度帯と一致) ことと、
**(b) digitize した $P_t$ と $Re'$ が互いに整合している** ことの両方を意味する。
ツール: `tools/conditions.py` → `derived.json`。

> **罠 (記録)**: 最初、NASA-9 多項式の下限を避けるため cp と**あわせて粘性も** 200 K で床止め
> した。$T_\infty \approx 50$ K では $\mu$ が真値の **4.7 倍**になり、$Re'$ が 3.5 倍外れて
> Gate A が解けなかった。床止めは cp だけでよい (Chapman–Enskog は $T^*\sim0.5$ まで有効で、
> 45 K でも Sutherland と 16 % 差)。$\mu_\infty$ は $Re'$ に直接乗るので、低温の粘性則は
> この case の一次的な不確かさ要因。

## Gate B — 分母 $h_{fp}$ の独立検算 (**PASS、ただし系統差あり**)

報告の分母は**較正パネルの実測から作った相関**であって理論ではない。Van Driest II
(Karman–Schoenherr + 圧縮性変換) と Reynolds 相似で独立に組み、桁と傾きを比べた。
BL は前縁から発達させず、Fig 4 の $\theta$ を $Re_\theta$ として直接与える
(すきまに近づくのはトンネル壁の厚い乱流 BL)。

| $Re'$ | $Re_\theta$ | $C_f$ (VD-II) | $h_{VDII}$ | $h_{corr}$ | 比 |
| --- | --- | --- | --- | --- | --- |
| 1.47e6 | 2.09e4 | 5.39e-4 | 1.975 | 1.432 | **1.38** |
| 3.32e6 | 3.95e4 | 4.75e-4 | 3.864 | 2.512 | **1.54** |
| 7.82e6 | 8.06e4 | 4.14e-4 | 7.661 | 4.536 | **1.69** |

($T_w$ = 300 K。400 K / 480 K でも比は 1.29–1.57 / 1.23–1.49 とほぼ変わらない。)

```
GATE B VERDICT: PASS  (独立理論が相関の 0.5-2.0 倍に入るか; 中央 1.434)
```

**桁は合うが系統差が残る**。しかもただのオフセットではなく**傾きが違う**:
Van Driest II の実効冪は $d\ln h/d\ln Re' = 0.81$、報告の相関は **0.69**。
digitize した 4 点自体から出る冪も 0.70 なので、これは読み取り誤差ではなく報告側の相関の性質。

**validation への含意 (重要)**: 報告の $h_{fp}$ は圧縮性乱流平板理論より**系統的に小さく、
高 $Re$ ほど小さい**。分母が小さければ報告の $h/h_{fp}$ は**その分だけ大きく出る**。
$Re'$ = 7.82e6 では 1.69 倍、すなわち報告値を理論平板で割り直すと **41 % 下がる**。
これは case/52 (T3) で見えた 1–2 桁の差を埋めるものではないが、
**「文献の比の絶対値」を CFD と突き合わせるときの 40–70 % の系統項**として
不確かさ予算に入れる必要がある。CFD 側は**報告と同じ演算** (実測相関で割る) で比を作り、
理論平板で割った値も併記する。

ツール: `tools/flatplate.py` → `gateB.json`。

## T2-G0 — 流入 BL の再構成と、$\theta$ の定義の特定

すきまに近づくのはトンネル壁の厚い乱流 BL なので、前縁から発達させず入口にプロファイルを
与える。**$C_f$ は仮定しない** (Gate B で相関と VD-II が 1.8 倍違うので、どちらかを入力にすると
流入場がその選択に汚染される)。代わりに実測の $\delta^*,\theta$ から $(\delta, u_\tau)$ を解く。

構成: van Driest 変換した $u_{vd}^+$ に Spalding 内層 + Coles wake、温度は Walz/Crocco-Busemann
($r$=0.89)、密度は等圧。

### 実測の $H$ は圧縮性定義では出ない

報告の形状係数は $H=\delta^*/\theta$ = 12.80/1.44 = **8.89**。しかし $M$ 10.3 では

| | $H$ |
| --- | --- |
| 圧縮性 $\delta^*,\theta$ (質量流束) | **24.0 – 40.4** ($T_w$ 300→995 K、速度冪 $1/7$→$1/20$ の全域) |
| 運動学 $\delta^*,\theta$ (速度のみ) | **1.2 – 1.9** |

$T_w$ も速度形状も振れない。$M$ 10.3 では「冷却壁」 ($T_w/T_{aw}$=0.30) でも $T_w/T_\infty$ = **5.7** と
自由流静温より遥かに高温で、BL 内密度が $\rho_e$ の 0.12–0.18 倍まで落ちるため $H$ が大きくなる。
圧縮性定義で $H$=8.9 を出すには $M_{eq}\approx6.7$ ($T_\infty\approx115$ K) が要り、$M$=10.3 と両立しない。

### 原典 Fig 4 を再読して digitize を確認した

**読み取りは正しかった**。Fig 4 右パネルは $w$ = 0.229 cm を図中に明記したうえで
$\delta^*/w$ と $\theta/w$ を別軸で直接プロットしており、$\delta^*/w$ = 56、$\theta/w$ = 6.3
($\to \theta$ = 1.44 cm) と読める。左パネルの $\delta^*$ [cm]・$\theta$ [cm] とも整合する。
**報告は本当に $H\approx9$ を載せている。**

### 原典 ref.6 を入手した — 運動学説は否定され、別の答えが出た

TN D-8233 の参考文献 6 は **TN D-7939** (`Pressure Gradient Effects on Heat Transfer to RSI
Tile-Array Gaps`, 1975; NTRS 19750020309) だった (ref.3 の TM X-71945 とは別物。両方取得済み)。
ここに BL プロファイルと**分母の平板データの両方**がある。

**(a) 定義は圧縮性だった**。式 (5)(6) は
$\delta^*=\int(1-\frac{\rho u}{\rho_e U_e})dz$、$\theta=\int\frac{\rho u}{\rho_e U_e}(1-\frac{u}{U_e})dz$。
→ **運動学定義という推定は否定された**。

**(b) $\rho$ と $u$ は独立に測られていない**。プロファイルは **11 本ピトー・レーク**の
ピトー圧から還元したもので、密度分布を得るには仮定が要る。その仮定は報告に書かれていない。

**(c) この BL は 2 次元平衡平板ではない**。TN D-7939 本文が明記している:

> "A significant change in the boundary-layer profile shape as a function of transverse position
> is observed. This transverse variation ... is an indication of the **three-dimensional character**
> of the boundary-layer flow in a nozzle of square cross section. The three dimensionality
> results from the **corner interaction** of the tunnel sidewall, floor, and ceiling boundary layers."

> "The tunnel nozzle expands the flow to a point approximately 1 meter upstream of the
> test-section center ... The intersection of the expansion-section and test-section walls then
> constitutes a **compression corner** for the wall boundary-layer flow. The data presented herein
> were obtained **in the compression region** downstream of this corner."

$\delta^*$ は横方向に大きく変わり (中心線で最大)、$Re'$ にはほとんど依らない、とも書いてある。

**(d) 分母の原式は $Re'$ の式ではなく $\delta^*$ の式**だった。TN D-7939 式 (7):

$$N_{St}\,R_{w,\theta}^{0.07} = A\,\delta^* + B,\qquad A=-2.77\times10^{-4}\ \mathrm{cm^{-1}},\ B=4.66\times10^{-3}$$

TN D-8233 の Appendix はこれを $h_{fp}=C\,Re'^{0.69}$ に**還元**したもの。$A\delta^*+B$ は
$\delta^*$=16.8 cm でゼロを切る鋭敏な式で、$\delta^*$ の横方向変化がそのまま効く。

**(e) 著者自身が適用範囲を限定している**:

> "**The correlation has no application to any flow other than this tunnel-wall boundary layer**;
> however, it indicates a good understanding of this flow."

**(f) digitize は原図 6 倍描画で確認した**。Appendix 図の縦軸は $h_{fp}$ [**W/(m²K)**] の対数軸で、
4 点は 1.45 / 2.5 / 3.7 / 4.65。$C$ = 7.95e-5・冪 0.69 は正しい。

### T2 の成立可否への含意 (重要)

$H$=8.9 が 2 次元平衡プロファイルで出ないことと、$h_{fp}$ が VD-II の 1/1.8 しかないことは、
**同じ 1 つの原因で説明がつく**: この基準流れは角干渉で 3 次元化し圧縮領域にある
トンネル壁 BL であって、平板ではない。

- **2D 平板 CFD で $h_{fp}$ の絶対値は再現できない**。T2-G0 の「smooth の熱伝達が合うこと」は
  原理的に通らない。
- 使えるとすれば**比 $h/h_{fp}$ のみ**で、それは「比が局所のエッジ条件と $\delta^*/W$ で決まり、
  3 次元性と圧力勾配が分子分母で相殺する」という**未検証の仮定**に乗る。
  これはまさに T2-G0 が試すはずだった仮定である。
- → 計画 §4.2b の「満たせない場合 T2 は探索に降格し、case/49 への受け渡しは行わない」に
  **該当する可能性が高い**。降格するなら、乱流一次検証は分母の定義が異なる他文献
  ([A78] TP-1187 / [WAC75] TM X-3225 / [A85] TP-2307 / [HN90] TP-2988) から選び直す。

$\theta$ の 2.5–2.8 倍差は、(b) のピトー還元の密度仮定と (c) の 3 次元性のどちらか
(あるいは両方) に帰着すると考えられるが、報告からは分離できない。**$\theta$ は T2 の拘束条件
として使わない**。$\delta^*$ は TM X-71945 本文の独立記述
(「$Re$ = 1×10⁶/ft で $\delta^*\approx$ 4.75 in = 12.07 cm」、Fig 4 の同条件 11.5 cm と 5 % 一致)
で裏が取れているので、こちらを錨にする。

### 訂正後の台帳 (`bl_derived.json`)

| $Re'$ [1/m] | $\delta^*$ [cm] (実測) | $\delta$ [cm] | $\theta_{圧縮性}$ [cm] | $Re_\theta$ | $\delta/W$ |
| --- | --- | --- | --- | --- | --- |
| 1.47e6 | 12.80 | 22.95 | 0.511 | 7.5e3 | 100 |
| 3.32e6 | 11.50 | 20.99 | 0.461 | 1.5e4 | 92 |
| 6.10e6 | 10.80 | 19.92 | 0.431 | 2.6e4 | 87 |
| 7.82e6 | 10.50 | 19.44 | 0.418 | 3.3e4 | 85 |

$\delta$ = 19–23 cm は 78.74 cm 角試験部にコア 33–39 cm を残して成立する
(報告の $\theta$ を圧縮性と読んで固定すると $\delta$ = 49–65 cm となり、対向壁の BL が重なって
$M$ 10.3 の一様コアが消えるため不可能 — これが定義の取り違えを疑った最初の手掛かりだった)。
$H_{圧縮性}$ が 4 系列とも 25.0–25.1 で揃うのは、プロファイル族が素直であることの傍証。

### Gate B への波及

$Re_\theta$ は圧縮性 $\theta$ (0.42–0.51 cm) で作るので、当初の digitize 値より 2.5–2.8 倍小さい。
Karman–Schoenherr の $C_f$ が上がり、$h_{VDII}/h_{corr}$ は **1.57 – 2.04 (中央 1.77)** となる。
分母の系統差は**約 1.8 倍**として扱う。

## 未着手 (T2-G0 以降)

- **T2 を一次検証として残すか、探索に降格するかの決定** ← 最優先 (上節 (a)–(f) を踏まえて)
- 降格する場合: 乱流一次検証の差し替え先を [A78] TP-1187 / [WAC75] TM X-3225 / [A85] TP-2307 /
  [HN90] TP-2988 から選定 (分母の定義が異なるので正規化から組み直す)
- Fig 5–8 の $h/h_{fp}$ 分布の digitize (現状 Fig 5(a) Re 1.47e6 の 34 点のみ)
- ~~2D 中央断面近似の妥当性~~ → **上節 (c) で否定的な答えが出た** (角干渉による 3 次元性 +
  圧縮領域)。残るのは「比なら相殺するか」の検証
- 流入 BL の **CSV 化** (プロファイルは組めた。降格しないと決めた場合に `inletProfile` 形式へ)
- 加熱板の長さと温度分布 (本文未記載)

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| — | **計算はまだ実施していない**。現状は入力台帳・Gate A・Gate B・T2-G0a (流入 BL 再構成) | `derived.json`, `gateB.json`, `bl_derived.json`, `inlet_bl.json` | — |
