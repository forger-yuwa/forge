# 極超音速の「すきま / 深いキャビティ」熱伝達 — 平板基準に対する既往知見のサーベイ

- **作成**: 2026-09-19 (初版)、**全面改訂**: 2026-09-19 (乱流・深いすきまの一次データを大幅追加)
- **目的**: 平板 (smooth wall) の熱伝達に対して、**すきま (gap) / 深いキャビティ内部**の熱伝達がどうなるかを、
  試験・相関・CFD の既往文献で押さえる。あわせて **forge の validation に使えるデータ**を選ぶ。
- **動機**: [`plans/active/case-plate-annular-cavity-m5.md`](../../plans/active/case-plate-annular-cavity-m5.md)
  (case/49: M5 平板 + 環状深キャビティ 幅 2.5 mm × 深さ 50 mm、**深さ/幅 $D/W$ = 20**、キャビティ壁 500 K 等温 ×
  平板断熱) の結果を「妥当」と言うための外部基準。検証計画は
  [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md)。
- **PDF 実体**: `papers/gap_heating/` (git 追跡外)。NTRS は ID で引ける:
  `https://ntrs.nasa.gov/api/citations/<ID>/downloads/<ID>.pdf`。**1970 年代のスキャンでも OCR テキスト層がある**
  (PyMuPDF で本文抽出可、図は 200 dpi レンダで digitize 可)。


## 2026-09-20 追記 (2) — 深部の実測値は対流熱流束ではない

薄肉過渡法 (thin-skin calorimetry) を使う文献 ([W70] TN D-5908, [TH76] TN D-8233,
[A78] TP-1187 …) の**深部の値は、薄板内の横方向伝導を含む**。本文が「表面伝導の補正なし」と
明記しているものもある。304SS・0.305 mm で $\sqrt{\alpha_s t}$ = **1.74 mm** (0.8 s) なので、
すきま幅が mm オーダーだと**拡散長がすきま幅を超える**。

CFD の $q_{conv}$ を熱源に薄板過渡を順方向に解き、実測と同じ $\rho c\tau\,dT/dt$ をかけて
比べると ([`case/50/tools/skin_smear.py`](../../case/50.deep_cavity_wieting_m7/tools/skin_smear.py))、
W70 $w/d$=0.063 の $z/W$=2 で **0.0016 → 0.112** (実測 0.135)、$z/W$=3 で 0.074 (実測 0.078)。
**「CFD が深部で 1–2 桁低い」は、比較していた量が違っただけだった。**
幅が 3.4 倍広い $w/d$=0.211 では拡散長が 0.4 幅しかなく補正も小さい — 補正量が
拡散長/すきま幅で決まるという予測どおりに振る舞う。

**文献の使い方**: 深部 ($z/W\gtrsim1$) の値を対流熱流束として CFD と直接比べない。
比較は必ず「CFD を薄板で滲ませてから」行う。薄板厚・材質・計測時間が要るので、
文献選定時にこれらが書かれているかを確認する。

## 2026-09-20 追記 — TN D-8233 / TN D-7939 の基準流れについて分かったこと

[TH76] TN D-8233 を T2 (乱流一次検証) に立てる前提で原典を読み込んだ結果、**分母 $h_{fp}$ の
素性**が分かった。ここは文献選定に直接効くので survey 側にも残す。

- TN D-8233 の**参考文献 6 は TN D-7939** (`Pressure Gradient Effects...`, NTRS 19750020309)。
  ref.3 の TM X-71945 とは別物。BL プロファイルと平板基準データの**両方**がこちらにある。
- 平板基準の**原式は $Re'$ の式ではない**: $N_{St}R_{w,\theta}^{0.07}=A\delta^*+B$
  ($A$=−2.77e−4 cm⁻¹, $B$=4.66e−3)。TN D-8233 Appendix の $h_{fp}=C\,Re'^{0.69}$ はその還元形。
- **基準流れは平板ではなくトンネル壁 BL** で、TN D-7939 自身が
  (i) 角干渉 (側壁・床・天井) による**3 次元性**でプロファイル形状が横方向に大きく変わること、
  (ii) データがノズル/試験部交差部の作る**圧縮領域**で取られたこと、
  (iii) **"The correlation has no application to any flow other than this tunnel-wall boundary layer"**
  を明記している。
- 独立検算 (Van Driest II + Reynolds 相似) では $h_{fp}$ は圧縮性乱流平板理論の **1/1.8**、
  $Re'$ 依存の冪も 0.69 vs 0.81 と違う。上の 3 点と整合する。
- 報告の $\delta^*/\theta$ = 8.9 は $M$ 10.3 の圧縮性定義では再現できない (どの $T_w$・速度冪でも
  $H\ge24$)。$\delta^*$ は TM X-71945 本文の独立記述で裏が取れるが、$\theta$ は取れない。

さらに、分母は **3 つの独立経路 (Appendix 相関 / Van Driest II / TN D-7939 原式 (7)) が最大 2 倍、
冪も 0.69 / 0.78 / 1.16 とばらつき、第一原理から再構成できない**。原式 (7) を閉じるには
実測の局所壁静圧と局所エッジ速度が要るが、どちらも digitize されていない。

**文献選定への含意**: TN D-8233 は「$h=q/(T_{aw}-T_w)$ で分母を定義し、壁温を振り、BL 厚を測った
唯一の乱流 deep gap データ」という点で依然として価値が高いが、**その分母は 2D 平板 CFD で
再現できる量ではない**。比 $h/h_{fp}$ の再現に賭けるか、一次検証を移すかの判断が要る。

**差し替え候補の一次評価 (2026-09-20)**: 原典を読んだ範囲では **[A78] TP-1187 が有力**。
分母が**同一パネルホルダ上の 2D 平板較正パネルの冷壁実測** (空力フェンスで一様 2 次元、
鋭利前縁 + BL トリップ) なので、**forge が分子も分母も同じ条件で計算でき**、実験側の正規化を
そのまま再現できる。$\delta^*/W$ も 4.7–9.0 と TN D-8233 の 46–56 より 1 桁小さく、
2D 中央断面近似に有利。理論とのずれ (層流 10 % / 乱流 30 % 高) を報告自身が明記している点も
TN D-8233 と対照的。気体 (メタン-空気燃焼ガス) は `case/50` の `CombustionProducts` が流用できる。
難点は、すきまが流れ方向でタイル配列・傾斜壁 (60°/75°/90°) をもち、配列交点の 3 次元性が残ること。
比較表は [`case/51 README`](../../case/51.gap_turbulent_thd8233/README.md)。
詳細は [`case/51 README`](../../case/51.gap_turbulent_thd8233/README.md) と
[計画 §4.2b-1](../../plans/active/case-hypersonic-gap-heating-validation.md)。


## 0. 要約 — 何が分かっているか

1. **すきま内の熱伝達は、開口部近傍を除いて平板より 1〜2 桁小さい**。深さ方向は単調減衰。
   - **層流**: 深さ 30 % で $q/q_{fp}\le0.1$ ([A78] TP-1187)。深いキャビティの面積平均は $0.5\!-\!1.1\,q_{fp}$ ([W70])。
   - **乱流**: 深さ 60 % で $q/q_{fp}\le0.1$ ([A78])。$z/W>4$ で **2 % 未満** ([TH75])。中深部で 10 % 未満 ([HN90])。
   - **乱流の方が深部まで熱が入る** (層流 30 % vs 乱流 60 %)。これは case/49 (乱流) が
     Wieting の層流データより深部が熱くなる方向を意味する。
2. **平板を超えるのは開口部・角・衝突域だけ**。すきま上端 $0.6\!-\!1.5\,q_{fp}$ ([W70])、
   縦すきま流が下流タイル前面に衝突する所で **2.4〜4.5 倍** ([A78])、段差付きタイルで 3.2 倍 ([WAC75])。
   **これらは「配列・交差・段差」があって初めて起きる** — 単独の横すきまでは起きない。
3. **深さ方向の減衰長は「すきま幅 $W$」でスケールする** ($z/W$)。相関はべき乗形
   ($h/h_{ref}\propto (z/W)^{-1}$ 程度, [CK75])。深さそのものでなく $z/W$ で整理するのが慣例。
4. **壁温不連続 (外面が熱く、すきま壁が冷たい) は効く** ([TH76] TN D-8233): $T_{surf}/T_{gap}$ を上げると
   **上端の熱伝達は増え、深部は減る**。反転は $2<z/W<3$。**case/49 は配位としては同型**だが
   (断熱平板 ~1190 K × 500 K 等温すきま、$T_{surf}/T_{gap}\approx2.4$ は試験の範囲外、かつ断熱板 vs 加熱等温板で
   熱的履歴が違う)、**この傾向を case/49 の深部にそのまま移すのは未検証の仮説**である
   (検証計画 §4.11 では「適用不確かさ = 未定量」として扱う)。少なくとも等温模型のデータをそのまま当てるのは誤り。
   **[TH76] の実測温度比は $Re$ 系列ごとに違う** (1.47e6: 1.00/1.19/1.37/1.61、3.32e6: 1.00/1.18/1.34/1.52、
   7.82e6: 1.00/1.17/1.38/1.59)。
5. **Reynolds 数は上端と深部で逆に効く** ([TH76]): $Re$ 増 → 上端の無次元熱伝達は減り、深部は増える
   (高エネルギー流の拡散が深くなる)。
6. **後壁 (下流側) が前壁より熱い**のは $z/W<4$ まで。$z/W>5$ では両壁ほぼ等しい ([TH76])。
   上端角から出るせん断層が下流壁に当たるため。
7. **すきまは通気 (venting) で決まる**: すきま床圧は外面圧に追従し ([TH75])、圧力勾配・通気があると
   流入が増えて加熱が上がる ([HN90] チャイン部、[D83] エレボン cove の leak area 依存)。
   **袋小路のすきま (case/49 のような閉じた環状) と、貫通するすきまは別物**。
8. **深部の絶対値は試験でも信用できない**: $q/q_{fp}\lesssim0.02$ は計測限界 ([A78], [TH75], [CK75] は負値も記録)。
   **CFD をこの帯域で「検証」できない**。

## 1. レジームの整理 (どの文献がどれを見ているか)

| レジーム | 代表文献 | case/49 との関係 |
| --- | --- | --- |
| **2D 矩形深キャビティ** (単独, 層流) | [W70] | ★ 幾何が最も近い ($D/W$ 16)。層流 |
| **単独横すきま** (平板, 乱流, 壁温不連続) | [TH76] | ★★★ **$D/W$=20.0 が完全一致**、乱流、壁温比可変 |
| **タイル配列のすきま** (縦/横/交差, 層流+乱流) | [WAC75] [A78] [A85] [J73] [Q75] | 交差・衝突加熱は case/49 には無い (環状で閉じる) が、深さ減衰は使える |
| **圧力勾配・通気つきすきま** | [TH75] [HN90] [SM79] | 偏心で周方向に圧力差がつく case/49 の派生に効く |
| **wing–elevon cove (貫通・漏れ流れつき深い溝)** | [DB78] [D83] [H80] [K78] [S77] | 貫通流の有無という点で対照。case/49 は閉じた溝 |
| **段差つきタイル / filler bar** | [P83] [WAC75] | 段差があると 3 倍級。case/49 は段差なし |
| **浅い/損傷キャビティ ($L/H$ 大)** | [EV08] [EV10] [EV11] [HO09] | **適用範囲外**の証拠として引く |
| **飛行データ** | [PM82] | 地上試験の外挿検証 |

## 2. 一次データ (層流)

### [W70] Wieting, *Heat-Transfer Distributions in Deep Cavities in Hypersonic Separated Flow*, NASA TN D-5908 (1970) — NTRS 19700027852

- **幾何**: 2D 矩形キャビティ。深さ $d$=20.3 mm 固定、幅 $w$=1.27/4.29/7.77/10.64 mm →
  $w/d$=0.063/0.211/0.383/0.524 (**$d/w$=16.0/4.7/2.6/1.9**)。スパン 76 mm (25 mm 版で 2D 性を確認)。
- **条件**: Langley 7-inch Mach 7 pilot tunnel、**メタン–空気燃焼生成物** ($Pr\approx0.75$)。$M$ 6.78–6.94、
  $T_t$ 1650/1880 K、$p_t$ 10.9–13.9 MPa、$Re_\infty$ 4.5–7.5e6 /m、**壁 294 K 冷壁**、層流。
- **基準**: $q_{fp}$ は**解析値** (Pohlhausen 層流平板 + Eckert 基準温度、キャビティ中点)。Table IV に記載
  (例: $M$ 6.90 / $Re$ 5.71e6 /m / $T_t$ 1650 K で **30.4 kW/m²**)。
- **結果**: $p_c/p_m\approx1$ (Fig 4, $Re$・スパン長に依らず)。後壁上端 $q_s/q_{fp}$=**0.6** ($w/d$ 0.063) 〜
  **1.5** (0.524)。**$x/d>0.5$ で 10 % 未満**。平均 $\bar q_c/q_{fp}$≈**1.1→0.5** (Fig 12; 200 dpi 試読では
  1.05–1.10 / 0.68–0.78 / 0.55–0.62 / 0.47–0.50)。総入熱は $w/d\gtrsim0.15$ で $Q_c<Q_{fp}$ (Fig 11)。
- **理論**: Burggraf (非粘性コア) は $w/d$=0.383/0.524 で分布まで一致。**深い 2 つは粘性コア**で実験が理論より高い。
  Chapman $\bar q/q_{fp}=0.60\,Pr^{0.33}$、Denison–Baum 0.61→0.56。
- **計測**: 薄肉 (0.30 mm SS) 過渡法、$q$ 精度 ±0.68 kW/m² (代表値の ±2.2 %)。
  すきま内定常化時間 $d^2/\nu\approx0.2$ s。

### その他の層流一次データ

| 文献 | 内容 |
| --- | --- |
| [F73] Foster, Lockman, Grifall, *TPS Gap Heating Rates of the Rockwell Flat Plate Model (OH2A/OH2B)*, NASA CR-134077 (1973) | **単一すきま**と交差部の層流データ (配列でない = 2D 再現向き)。未入手 |
| [N66/N70] Nestler, *Laminar Heat Transfer to Cavities in Hypersonic Low Density Flow* (3rd IHTC 1966) / *Hypersonic Laminar Cavity Heat Transfer* (4th IHTC 1970) | 層流キャビティの古典 |
| [NI64] Nicoll, *A Study of Laminar Hypersonic Cavity Flows*, AIAA J. 2(9) 1964 / ARL 63-73 | 熱伝達と回復係数。薄肉過渡法の限界も論じている |
| [H69] Hahn, *Experimental Investigation of Separated Flow over a Cavity at Hypersonic Speed*, AIAA J. 7(6) 1969 | キャビティ分離流 |
| [CH61] Charwat et al., *An Investigation of Separated Flows* I/II, J. Aerospace Sci. 28 (1961) | 空洞流れの原典 (圧力場 / 空洞内流れと熱伝達) |

## 3. 一次データ (**乱流** — 本改訂の主眼)

### [TH76] Throckmorton, *Effect of a Surface-to-Gap Temperature Discontinuity on the Heat Transfer to RSI Tile Gaps*, NASA TN D-8233 (1976) — NTRS 19760019344 ★**最重要**

- **幾何**: **平板 + 単独すきま**模型 (配列でない)。すきま **幅 2.29 mm × 深さ 45.7 mm → $D/W$ = 20.0**
  (case/49 の 2.5 mm × 50 mm = 20.0 と**一致**)。薄肉 0.406 mm SS、midspan に 21 熱電対。
  掃引角 $\Lambda$=0°(横) 〜 60°。
- **条件**: Langley continuous-flow hypersonic tunnel、$M$=10.3、$Re_\infty$=1.5/3.3/7.8e6 /m、
  **厚い乱流トンネル側壁 BL** ($\delta^*,\theta$ は TN D-7939 Fig 4)。
  **$T_{surf}/T_{gap}$ = 1.0 / 1.2 / 1.4 / 1.6** (すきま上流の平板を電熱で加熱して作る)。
- **模型の熱構造 (重要)**: 加熱真鍮板 (3.175 mm, 4 台の電熱ヒータ) → **断熱材 + 水冷通路で熱的に絶縁**された
  薄肉 (0.406 mm, 304SS) すきま → **非加熱** 17-4PH 後板。**異温度の等温面を直接突き合わせた形ではない**。
- **基準と定義**: $h=q/(T_{aw}-T_w)$、$T_{aw}=T_\infty+0.89\,(T_t-T_\infty)$ (乱流 $r$=0.89)。
  無次元化の分母 $h_{fp}$ は**平滑平板の実験相関** ($N_{St}Re^{x}_{w,\theta}$=const を $h_{fp}=C\,(Re'_\infty)^{0.69}$ に還元、
  Appendix)。**test ごとの smooth 実測値ではない**。$T_{surf}/T_{gap}$=1.0 が等温 (従来型) 条件。
- **実際の温度比は 1.00 / 1.18 / 1.34 / 1.52** (Fig 5(b))。公称 1.0/1.2/1.4/1.6 と異なる。
- **測定範囲**: 基本分布 (Fig 5–7) は主に **$z/W\lesssim6$**。黒塗り記号は「極端に低く精度に疑義のある」データと本文が明記。
  → **$z/W\simeq20$ の底面を拘束する根拠にはならない**。
- **結果**:
  - 前壁・後壁とも $h$ は深さとともに単調減少。**後壁 > 前壁 ($z/W<4$)**、$z/W>5$ でほぼ等しい
    (上端角からのせん断層が後壁に衝突するため)。
  - **$T_{surf}/T_{gap}$ ↑ → 上端の $h$ ↑・深部の $h$ ↓、反転は $2<z/W<3$**。
    機構は「壁温のステップ減少 → 壁近傍温度勾配の不連続 → 局所熱伝達増」+ 「上端で熱を取られた分、
    深部に運べるエネルギーが減る」。
  - **$Re$ ↑ → 上端の無次元 $h$ ↓・深部 ↑** (高エネルギー流の拡散が深くなる)。$Re$ が上がるほど
    非等温効果は薄まる。
  - 掃引角 0–60° の影響は小さい。
- **なぜ効くか (case/49 への含意)**: case/49 は**断熱平板 (回復温度 ~1190 K) × 500 K 等温すきま**で
  $T_{surf}/T_{gap}\approx2.4$。[TH76] の傾向をそのまま外挿すると、**上端の熱流束は等温条件より高く、
  深部は低く**出るはず。等温模型ベースの相関 ([CK75] など) をそのまま当てると深部を過大評価する。

### [WAC75] Weinstein, Avery, Chapman, *Aerodynamic Heating to the Gaps and Surfaces of Simulated RSI Tile Arrays in Turbulent Flow at Mach 6.6*, NASA TM X-3225 (1975) — NTRS 19760004291

- **幾何**: 46 cm 角のタイル配列 (中央 1 枚が計測用 304SS、壁厚 0.51 mm、15.2 cm 角)。
  **すきま深さ 6.50 cm**、幅 0.10/0.18/0.30/0.41 cm (**$D/W$ = 16 〜 65**)、エッジ半径 0.25 cm。
  in-line / staggered 配列、タイル突出 (段差) 有無。
- **条件**: Langley 8-ft High-Temperature Structures Tunnel、**メタン–空気燃焼生成物**、$M$=6.6、
  $T_t$≈1750 K、$Re_\infty$=2.0–4.9e6 /m、$\delta^*$=0.81–1.62 cm (**$\delta^*\gg W$**)、乱流。
- **基準**: 計測タイル中央の熱流束 (較正パネル中央と 5 % 以内) を**等価平板値**として無次元化。
- **結果**: in-line 最大 1.8 $q_{fp}$ (前縁 R 直後の上面)。staggered は縦横すきま交差部の前面で **2.9 倍**、
  段差 (δ* の 20 % 突出) で **3.2 倍**。**側壁の熱伝達は上端で平板値程度、深さとともに急減**。
  縦すきま流が横すきまへ溢れる角部では平板の 1.2 倍。0.5 cm 深さで前面中央は平板の 2 倍超 (staggered)。
  垂直壁の積分熱流束は in-line が staggered より 40 % 大、タイル全体の総熱負荷は 13 % 大。
- **含意**: **配列・交差があると上端付近は平板超**。case/49 の環状すきま (交差なし・閉じている) は
  むしろ [TH76] の単独すきまに近い。

### [A78] Avery, *Aerodynamic Heating in Gaps of TPS Tile Arrays in Laminar and Turbulent Boundary Layers*, NASA TP-1187 (1978) — NTRS 19780020430

- **幾何**: staggered 金属タイル配列、前面傾斜 $\theta$=60/75/90°、縦すきま長 15.24/30.48 cm、
  幅 0.10/0.18/0.30/0.41 cm。$\delta^*$=0.36–1.62 cm。
- **条件**: 8-ft HTST、$M$≈7、$T_t$≈1800 K、$Re_\infty$=1.0–4.8e6 /m。**層流と乱流の両方**。
- **結果 (深さ依存の決定的な数値)**:
  - **層流**: 前面の最大でも深さ **30 %** で $0.1\,q_{FP}$ に落ちる。
  - **乱流**: 深さ **60 %** で $0.1\,q_{FP}$ 以下。**上端に「一定値の領域」**があり、
    $\theta$=90° で $W$=0.30/0.41 cm なら深さの上 **15 %**、他は上 **5 %**。
  - 最大衝突加熱は乱流で **4.5 $q_{FP}$** ($\theta$=90°, $W$=0.41 cm)。$\theta$ を 60° に寝かせると大幅低下。
    層流では $W$=0.41 cm で 2.4 倍、0.10/0.18 cm では**平板以下**。
  - 乱流では横すきま内の流れは 3 次元で、影響範囲はすきま幅の 3 倍程度に及ぶ。
  - $q/q_{FP}<0.02$ のデータは**計測系の限界**で信用しない、と明記。
  - **衝突加熱の経験相関** (回帰 + 図式) を層流・乱流それぞれに提示。
- **含意**: 「乱流は層流より深くまで熱が入る」を定量化した唯一級の資料。case/49 の乱流深部評価の目安。

### その他の乱流一次データ

| 文献 | NTRS | 内容 |
| --- | --- | --- |
| [TH75] Throckmorton, *Pressure Gradient Effects on Heat Transfer to RSI Tile-Array Gaps*, **NASA TN D-7939** (1975) | 19750020309 | 実寸模擬タイル配列 (15.25 cm 角)、縦 0.30 / 横 0.20 cm、**深さ 2.86 cm**、$M$=10.3、$Re$=1.6/3.3/6.1e6 /m、厚い乱流壁 BL。**中心線外・$z/W>4$ で非攪乱面の 2 % 未満**。横圧力勾配の系統的影響なし。基準はすきまを埋めた smooth model 実測 |
| [A85] Avery, *Experimental Aerodynamic Heating to Simulated Shuttle Tiles in Laminar and Turbulent BL with Variable Flow Angles at Mach 7*, **NASA TP-2307** (1985) | 19850026038 | 流れ角 0/15/30/45/60°、T-gap。局所ピークと総熱負荷。**流れ角 30–50° で局所加熱最小** |
| [J73] Johnson, *Heat Transfer Data to Cavities Between Simulated RSI Tiles at Mach 8*, NASA CR-128770 / *Interference Heating to Cavities…* (1973) | 19730024766 | すきま幅と BL 厚の効果 (乱流, $M$ 8)。交差部の問題を最初に指摘 |
| [Q75] Quan & Lockman, *Convective Heating Tests of a Longitudinal Gap on the Rockwell Flat Plate Model*, NASA CR-141539 (1975) | 19750022163 | **単独の縦すきま**の対流加熱試験 |
| [HN90] Hunt & Notestine, *Aerodynamic Pressure and Heating-Rate Distributions in Tile Gaps Around Chine Regions with Pressure Gradients*, **NASA TP-2988** (1990) | 19900014354 | 8-ft HTT、$M$=6.6、$T_t$≈1890 K、$\alpha$=7/10/13°。**すきま圧は幅が広いほど低い (通気) → 流入増 → 加熱増**。**周方向すきまの中深部は外面値の約 10 % 未満**。T 字接合と段差が最悪 |
| [H88] Hunt, *Aerodynamic Pressures and Heating Rates on Surfaces Between Split Elevons at Mach 6.6*, NASA TP-2855 (1988) | 19890003451 | 分割エレボン間の深い隙間 |
| [SM79] Scott & Maraia, *Gap Heating with Pressure Gradients*, AIAA 79-1043 (1979) | 19790054016 (PDF なし) | JSC 側の圧力勾配つきすきま加熱 |
| [PM82] Pitts & Murbach, *Flight Measurements of Tile Gap Heating on the Space Shuttle*, AIAA 82-0840 (1982) | 19820050466 (PDF なし) | **飛行実測**。地上データの外挿妥当性 |
| [P83] Petley & Smith, *Analysis of Gap Heating due to Stepped Tiles in the Shuttle TPS*, NASA TP-2209 (1983) | 19830023442 | 段差タイルによるすきま内流れと filler bar 焼損の解析 (流れ + 熱の連成モデル) |
| [HO23] Hollis & Hollingsworth, *Experimental Investigation of Block-TPS Fence/Gap Roughness Effects on Transition Onset and Turbulent Heating*, NASA/TM-20230001707 (2023) | 20230001707 (PDF 直リンク不可) | **現代の計測** (phosphor thermography)。ブロック TPS のすきま・フェンスによる遷移と乱流加熱 |

### wing–elevon cove (貫通流つきの深い溝) — 対照群

| 文献 | NTRS | 内容 |
| --- | --- | --- |
| [D83] Deveikis, *Effects of Flow Separation and Cove Leakage on Pressure and Heat-Transfer Distributions along a Wing–Cove–Elevon Configuration*, **NASA TP-2127** (1983) | 19830020107 | $M$=6.9、実寸ヒートシンク模型、**cove 幅 12.7 mm**、スパン 1.05 m、シール漏れ面積 0/13/50/100 %、ランプ角 15–35°、$Re$=1.15–4.5e6 /m、$T_w/T_t\approx0.17$。準層流/遷移/乱流分離 |
| [DB78] Deveikis & Bartlett, NASA TM-74095 / AIAA 78-39 (1978) | 19790010113 | 同系列の初期試験 (漏れ可変 cove の圧力・熱伝達) |
| [H80] Hunt R. L., *Aerothermal Analysis of a Wing–Elevon Cove with Variable Leakage*, NASA TP-1703 (1980) | 19800024184 | 上の解析版 |
| [K78] Keshock, AIAA 78-40 (1978) / [S77] Scott & Murray, AIAA 77-757 (1977) / [B83] Bey, NASA TM-85711 (1983) | — | cove 内部流れによる加熱の解析・予測法 |

## 4. 相関・予測法

| 文献 | 形 |
| --- | --- |
| **[CK75]** Christensen & Kipp, *Data Correlation and Analysis of Arc Tunnel and Wind Tunnel Tests of RSI Joints and Gaps*: **Phase I = NASA CR-134345/134346 (1974)**、**Phase II = NASA CR-141925 / JSC-09651 / MDC E1248 (1975, NTRS 19750020031)** | 多施設データ (Ames 3.5-ft HWT / JSC 10 MW アークジェット / LaRC 8-ft HTST) の多変量回帰。例: in-line gap 層流 風上側 $\ln(h/h_{ref})=-1.00604-0.96104\ln(z/W)$ → $h/h_{ref}\approx0.366\,(z/W)^{-0.96}$。エッジ半径 $E$・表面距離 $S$ を含む形も。668 点中 316 点が層流。**すきま内 $q$ は逆解析 (HEATRAN) 由来で、負値・tile↔plug 伝導の異常を本文が明記** |
| [CR79] *Prediction of In-Depth Gap Heating Ratios from Wing Glove Model Test Data*, NASA CR-160146 (1977) / *…from NASA/Ames Double Wedge Model Test*, CR-160147 (1978) | JSC 10 MW アークジェット試験から **$q(z)/q_{ref}$ を直接**予測式化 (Rockwell) |
| [A78] の衝突加熱相関 | 縦すきま流の衝突加熱を BL・幾何パラメータで整理 (層流/乱流別) |
| [EV08] Everhart, JSR 2009 (NTRS 20080008346) / [EV10] *Turbulent Supersonic/Hypersonic Heating Correlations for Open and Closed Cavities*, JSR 2010 (NTRS 20090007685) / [EV11] *In-Cavity Transition and Heating Augmentation*, AIAA 2011 (20110013257) | **浅い/損傷キャビティ** ($L/H<10$ open, $>13$ closed) の bump factor 相関と予測限界 (95/99/99.9 %)。**$D/W$=20 のすきまには適用外** |
| Burggraf (1965) / Chapman (1956) / Denison–Baum (1963) | 非粘性コア理論・自由せん断層理論。$\bar q/q_{fp}=0.60\,Pr^{0.33}$ など |

## 5. 数値解析の既往 (比較・注意)

- **PLOS ONE 10(1):e0117012 (2015)** ([open](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0117012)):
  2D、幅 2 mm、$L/D$=1/2〜1/6、$M$=2/3/4、壁 300 K、S-A。すきま壁の $q/q_0$ は **U 字**分布、風上角で最大。
  gap effect coefficient は $M$・深さ増で減、迎角増で増。面取りで 7.03→1.36。$M$=5 の中国側試験と比較。
- **Phys. Fluids 37 (2025) 125166** (横すきま)、**Acta Astronautica (2024)** (縦すきま + 漏れ) — 近年の RANS 研究。
- **Wood/Pulsonetti (LAURA)** ほか: キャビティ流れの CFD 予測は難しいと Everhart が明記。
- **乱流モデル**: SWTBLI で標準 SST はピーク熱流束を ~40 % 過小評価するという最近の報告あり
  ([ScienceDirect, 2024](https://www.sciencedirect.com/science/article/abs/pii/S0094576524003345))。
  すきま深部は**ほぼ淀み**で、渦粘性が残ると熱を運びすぎる恐れがある → **層流帯での A/B が必須**。
- 参考: $L/H$=14–20・$H/\delta$=1.1–5.2 のキャビティで、$Re_\theta$=300 (層流) は実験と一致するが
  $Re_\theta$=503 (遷移/乱流) では加熱増幅が桁で変わる、という報告もある (キャビティ内遷移の敏感さ)。

- **TP-1187 を数値計算で再現した既往は見つからなかった (2026-09-26 検索)**。Web 検索 (NTRS・Google Scholar 相当の一般検索) で、
  TP-1187 の熱電対データと CFD を突き合わせた文献は出てこない。同時代の後続は実験 (TP-2307 [A85]、1983 年の会議論文) か
  工学的解析 (Petley ほか NASA TP-2209 (1983) は段差タイルのすきま流れの解析で、本文に TP-1187 の引用は見当たらない)。
  近年の RANS 研究 (上の PLOS ONE 2015、Acta Astronautica 2024、Phys. Fluids 2025) は中国の風洞試験など自前のデータで検証しており、
  TP-1187 は使っていない。原報自身の解析は経験式 (eq. (3)/(4)) の当てはめまで。**forge の照合が公開の CFD 再現としては初出の可能性が高い**
  (網羅的ではない。引用索引で TP-1187 の被引用を洗い直す余地はある)。

- **TP-1187 の乱流分母の出典は Deveikis & Hunt, NASA TN D-7275 (1973)** (NTRS 19730021511)。同じ 8-ft HTST の大型校正パネルで、中心線スタントン数分布・局所条件表・平均加熱率の carpet plot がある。精読メモ: [`tp1187-test-system-reading.md`](tp1187-test-system-reading.md)。

## 6. 検証に何を使うか (優先順)

| 優先 | 文献 | 役割 | 決め手 |
| --- | --- | --- | --- |
| 1 | **[TH76] TN D-8233** | **乱流・単独すきま・$D/W$=20.0・壁温比可変 (実測値は $Re$ 系列ごとに違う)** | case/49 と幾何比・壁温配位が同型。**分母は平板の実験相関** $h_{fp}=C(Re'_\infty)^{0.69}$ (test ごとの smooth 実測ではない)。基本分布は $z/W\lesssim6$ |
| 2 | **[W70] TN D-5908** | **層流・2D 深キャビティ・4 幅** | 乱流モデル無しで離散化と再循環を素で測れる。平均比・圧力比まで問える |
| 3 | **[A78] TP-1187** | 層流 vs 乱流の**深さ減衰 (30 % vs 60 % で 0.1 $q_{FP}$)** | モデル形式差の外部基準 |
| 4 | [CK75] / [CR79] | 相関帯 | 点一致でなく帯で見る |
| 5 | [WAC75] TM X-3225 / [HN90] TP-2988 | 上端・配列・通気の上限値 | case/49 に交差は無いが、上端の桁を外していないかの確認 |
| — | [EV08/10] | **適用範囲外**の明示 | 浅いキャビティ相関を外挿しない根拠 |

## 7. 未入手 (必要になったら取りに行く)

- Nestler, *The Effects of Surface Discontinuities on Convective Heat Transfer in Hypersonic Flow*, AIAA 85-0971 (1985) — 総説の定番
- Fletcher, Briggs, Page, AIAA 70-767 (1970) — 分離・再付着熱伝達の総説
- Nestler, Saydah, Auxer, *Heat transfer to steps and cavities in hypersonic turbulent flow*, AIAA J. 7(7) 1969 ([10.2514/3.5351](https://doi.org/10.2514/3.5351))
- Dunavant & Throckmorton, JSR 11(6) 1974 (交差部を $\delta^*$・走り長・幅・深さで整理)
- Brewer, Saydah, Nestler, Florence, JSR 10(1) 1973
- [SM79] AIAA 79-1043、[PM82] AIAA 82-0840 (NTRS に PDF 無し)
- Tang G. M., 実験流体力学 14(4) 2000 (**深いすきま**, $M$ 9.85/12/15.5)、Jiang & Zuo, 推进技术 20(1) 1999 ($M$ 5)
- Mori K., AIAA 2012 (TSP による深いキャビティ内熱流束の光学計測)
- DLR: IXV TPS パネル (L3K アーク加熱) のすきま/界面熱流束、SHEFEX の facet 間すきま

## 8. ありがちな誤り (検証計画に落とすべき注意)

1. **分母 (平板基準) の取り違え**: (i) 同一位置の平滑模型実測 ([TH75] [TH76] [WAC75])、(ii) 層流平板解析式 +
   Eckert ([W70])、(iii) タイル中央実測 ([WAC75]) と定義が違う。**比だけ引用すると 2 倍ずれる**。
2. **等温データを非等温配位に当てる**: [TH76] が示す通り、$T_{surf}/T_{gap}$ で上端と深部が逆向きに動く。
   case/49 は $T_{surf}/T_{gap}\approx2.4$ で試験範囲の外。**逆に、[TH76] の傾向を case/49 の深部へ外挿するのも未検証**。
3. **層流データで乱流を語る**: 深さ 30 % vs 60 % ([A78])。**桁が違う**。
4. **相関の外挿**: [EV08/10] は浅いキャビティ、[CK75] は $z/W\lesssim10$。$D/W$=20 は外側。
5. **計測限界以下での「一致」**: $q/q_{fp}<0.02$ は試験側が信用していない帯域。
6. **定常化時間の過小評価**: すきま内は $d^2/\nu$ 程度 ([W70])。外部通過時間の 3 桁上。
7. **通気・貫通の有無を混ぜる**: 閉じた環状すきま (case/49) と、貫通する縦すきま・cove は別物
   ([HN90] [D83])。
