# carrier 中の非等温核生成補正・過冷却水の表面張力・空気 (N2/O2) 凝縮 onset — 文献調査

## メタ

- **種別**: 技術調査 (コード変更は plan 側: [condensation-kantrowitz-carrier](../../plans/active/condensation-kantrowitz-carrier.md), [condensation-air](../../plans/active/condensation-air.md))
- **日付**: 2026-09-12
- **動機**: ユーザ指摘「carrier (N2) 中の H2O 凝縮では Kantrowitz 補正はキャリアによる冷却を考慮すべきでは」「表面張力の式は 50 nm 級でも成り立つか」「空気そのものの凝縮解析機能を検証したい」。

## 1. 非等温核生成補正 (Kantrowitz / Feder) と carrier gas

### 1.1 Feder et al. (1966) の一般形

Feder, Russell, Lothe, Pound, *Adv. Phys.* 15 (1966) 111: 定常核生成率の非等温補正は

$$
\frac{J_{noniso}}{J_{iso}}=\frac{b^2}{b^2+q^2},\qquad
q = h - \frac{k_BT}{2} - \gamma\frac{\partial A(n)}{\partial n},\qquad
b^2 = k_B T^2\Big[\big(c_{v,v}+\tfrac{k_B}{2}\big) + \frac{N_c}{N_v}\sqrt{\frac{m_v}{m_c}}\big(c_{v,c}+\tfrac{k_B}{2}\big)\Big]
$$

($h$: 1 分子あたり潜熱、$c_v$: 1 分子あたり定積熱容量、$N_c/N_v$: キャリア/蒸気の分子数比、$m$: 分子質量)。
単原子蒸気+単原子キャリアで $b^2=2k_B^2T^2\,[1+(N_c/N)\sqrt{m/m_c}]$ (Wedekind et al. 2008 式 9; [arXiv:0804.1516](https://arxiv.org/abs/0804.1516))。
$b^2$ は「衝突 1 回で持ち去れるエネルギー揺らぎ」で、衝突頻度 $\propto p_i/\sqrt{m_i}$ の重みでキャリアも寄与する。
表面項 $\gamma\,\partial A/\partial n$ は臨界核で $2\sigma v_l/r_*=k_BT\ln S$ (Kelvin–Thomson) に等しい。

### 1.2 Kantrowitz (1951) 純蒸気形との関係

$N_c=0$、表面項を落とし $q\to h-k_BT/2$、$c_{v,v}+k_B/2=k_B\frac{\gamma_v+1}{2(\gamma_v-1)}$ とすると
$\theta\equiv q^2/b^2=\frac{2(\gamma_v-1)}{\gamma_v+1}\big(b_L-\tfrac12\big)^2$, $b_L=L/(R_vT)$ — forge の `condKantrowitz: 1`
($\theta=\frac{2(\gamma_v-1)}{\gamma_v+1}b_L(b_L-\tfrac12)$, Kantrowitz の原形) と $b_L\approx25$ で 2 % 差。

### 1.3 Wysłouzil 条件 (H2O 1 % in N2) での見積り

T=230 K, $Y_w$=0.011, $g\approx0$: $N_c/N_v=(0.989/28.01)/(0.011/18.02)=58$, $\sqrt{m_v/m_c}=0.80$,
$\tilde c_{v,v}+\tfrac12=3.5$ (H2O 蒸気 $c_v$=1393 J/kg/K → 3.0 $k_B$), $\tilde c_{v,c}+\tfrac12=3.0$ (N2)。
$b^2/(k_B^2T^2)$: 純蒸気 3.5 → carrier 3.5+58×0.80×3.0=143 (**41 倍**)。
$\theta$: 純蒸気 167 → carrier **4.1** ($J/J_{iso}$ 0.006 → **0.20**)。表面項込み ($\ln S\approx3.4$, $q$ ×0.86) で $\theta$ 3.0 → 0.25。
つまり carrier 形は「等温 (1.0)」と「純蒸気形 (0.006)」の間、等温寄りに落ちる。forge 従来コメント「carrier では θ→~1」は定量的には
θ≈3–4 (J は 1/4–1/5)。

### 1.4 pV 仕事 (Wedekind の "pressure effect")

$\Delta\mu_{eff}=k_BT\ln S - v_l(p+p_c-p_{eq})$: Wysłouzil ($p\approx10$ kPa, $v_l=3\times10^{-29}$ m³) で $v_lp_c=3\times10^{-25}$ J vs
$k_BT\ln S=1.1\times10^{-20}$ J → 3×10⁻⁵ で無視。実装しない。

## 2. 過冷却水の表面張力と小半径 (Kelvin / 臨界核) での妥当性

- forge `h2o_sigma`: IAPWS R1-76 形 $\sigma=0.2358\,\tau^{1.256}(1-0.625\tau)$, $\tau=1-T/647.096$ を 273 K 未満へ外挿。
- 過冷却域の実測: Hrubý et al. 2014 (*J. Phys. Chem. Lett.*) と Vinš et al. 2015 (*JPC B*, 対圧毛管上昇法, −25 °C まで) は **IAPWS の滑らかな
  外挿と一致し、第二変曲点 (歴史的データの −9 °C の折れ) の証拠なし** ([PMC5018856](https://pmc.ncbi.nlm.nih.gov/articles/PMC5018856/),
  [pubmed 26276586](https://pubmed.ncbi.nlm.nih.gov/26276586/))。一方 Vinš et al. 2020 (*JPC Lett.*, −31.4 °C=241.8 K まで,
  [JPC Lett. 2020](https://pubs.acs.org/doi/10.1021/acs.jpclett.0c01163)) は **−20 °C 未満で外挿 IAPWS 式から有意な偏差を検出し、深い過冷却域で
  異常の余地がある**と報告 (2026-09-12 codex 指摘で訂正: 初稿は「一致」と誤読)。MD (TIP4P/2005 等) は 227–235 K 以下で IAPWS からの逸脱の尾を示唆。
  **結論: 平面 σ の IAPWS 外挿は T ≳ 250 K では実測で支持されるが、Wysłouzil の核生成温度 (~215–235 K) は実測域外で精度未保証。**
- 曲率依存 (Tolman): $\sigma(r)=\sigma_\infty(1-2\delta/r)$。水の $\delta$ は Wilhelmsen, Bedeaux, Reguera 2015 (**square-gradient theory + CPA EOS**、
  MD ではない; [JCP 142, 171103](https://pubs.aip.org/aip/jcp/article/142/17/171103/939483)) で負の小さい値 (~−0.05 nm)。
  - $r=50$ nm: 補正 0.2 % → **成長した液滴 (Kelvin 項・成長則) では平面 σ で十分**。
  - $r_*\approx1$ nm (臨界核; $v_l=3.0\times10^{-29}$ m³ → **~139 分子**): 補正 ~+10 %。$J\propto\exp(-\Delta G^*/k_BT)$, $\Delta G^*\propto\sigma^3$ なので
    σ +10 % で障壁は $[(1.1)^3-1]\,\Delta G^*/k_BT\approx17$–23 $k_BT$ 増 ($\Delta G^*/k_BT\sim50$–70; 一定倍率の効果で、引用論文の曲率補正が
    障壁を +5–6 $k_BT$ 動かす報告とは別の計算)。**臨界核サイズでは capillarity 近似そのものが精度限界**で、CNT+平面 σ は「経験的に較正された妥協」
    (Wölk & Strey 2001 の水の経験補正がその例)。
- 実務的含意: σ の ±3 % で $J$ が 10² 桁級に動く → onset 位置の一定倍率局所感度 (`condSigmaScale`) で不確かさの効き方を示す。Tolman 補正は
  対象温度・臨界核サイズでの曲率モデルと係数の検証が不足しているので入れない。

## 3. 空気 (N2/O2) 自体の凝縮 — onset データ

### 3.1 一次データ

- **Daum & Gyarmathy (1968), AIAA J. 6(3) 458–465** "Condensation of air and nitrogen in hypersonic wind tunnels": 各種風洞 (conical / wedge /
  contoured) の onset (p, T) を集約。**低圧域では空気は実質純 N2 として振る舞い、N2 の自発核生成が onset を決める**。膨張率パラメータ
  $\dot P=-(1/p)\,dp/dt$ [1/s] と静圧の比が onset を相関。
- **Grossir & Rambaud (2014), AIAA 2014-1153** ([手元 PDF](../../papers/condensation/Grossir_Rambaud_2014_AIAA-2014-1153_nitrogen_condensation_detection_VKI_Longshot.pdf))
  Fig. 4 に上記を再掲 (p–T 平面、air/N2 の飽和線、N2 の「最小実験 onset 曲線」、$\dot P$ 別の理論 onset 線)。**最低圧 (~1 Pa) で飽和線より
  20 K 低温まで過冷却**、1 kPa で ~8 K。輪郭ノズル ($\dot P<1000$/s) はほぼ過冷却なし、conical/wedge ($\dot P$ 10⁴/s 級) で大きい。
- **Hansen & Nothwang (1952), NACA TN-2690** ([手元 PDF](../../papers/condensation/Hansen_Nothwang_1952_NACA-TN-2690_air_condensation_supersonic_tunnels.pdf)):
  空気の飽和蒸気圧を **O2/N2 の理想溶液** ($x_{O_2}$ を平衡から決め $p_w=x_{O_2}p_{O_2}+(1-x_{O_2})p_{N_2}$, 式 A1–A2) で与え、
  Becker–Döring + **Tolman 補正 σ** で上限 (自発核生成) を、飽和線で下限を与える。閾値は $J=10^{12}$ /cm³/s。
- Arthur (1952, case/34) は $P_0$=844 kPa, $T_0$=290 K の N2 source-flow ノズルで、onset 付近 (M≈5.4) の $p\approx1$ kPa・$T\approx43$ K は
  上の最小 onset 曲線の 1 kPa 点に相当する。

### 3.2 最小実験 onset 曲線 (N2) の読み取り値 (Grossir Fig. 4a の破線; 読み取り誤差 ±1.5 K)

| p [Pa] | T_onset [K] | T_sat,N2 [K] (forge `n2_psat` 逆算) | 過冷却 [K] |
| --- | --- | --- | --- |
| 0.3 | 11 | ~31 | ~20 |
| 1 | 14 | ~34 | ~20 |
| 10 | 21 | ~38.5 | ~17 |
| 100 | 32 | ~44 | ~12 |
| 1000 | 43 | ~51 | ~8 |
| 10000 | 53 | ~62 | ~9 |

CSV: [`case/34.arthur_n2_nozzle/daum_gyarmathy_min_onset_n2.csv`](../../case/34.arthur_n2_nozzle/daum_gyarmathy_min_onset_n2.csv)。
空気の飽和線は N2 より ~1 K 高温側 (Fig. 4a)。

### 3.3 forge の現状と差分

- forge の pure-condensible は **N2 のみ** (`condModel 0`: Lin 2014 の物性、CNT×Iland、Goodheart 成長、CPG)。空気は無い。
- `n2_latent` (Lin 式 26 の 4 次多項式) は **60 K 未満で $dL/dT>0$** ($c_l=c_{p,v}-L'<0$: 45 K で −3100 J/kg/K) と熱力学的に不整合
  (2026-09-10 codex 指摘)。Arthur の最低温 (~28–35 K) はこの域に入る。
- 空気凝縮の検証は「空気 ≈ N2」(Daum & Gyarmathy) を前提に **N2 核生成・成長に空気の飽和線 (O2/N2 理想溶液)・混合物性を与えた擬似種**で行うのが
  一次近似。O2 の物性: NIST Antoine ($\log_{10}p[\mathrm{bar}]=3.9523-340.024/(T-4.144)$, 54–154 K)、$\rho_l$ 1141 kg/m³ @90.19 K、
  σ 13.2 mN/m @90.19 K、$L$ 6820 J/mol @NBP ([NIST WebBook](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=4))。

## 4. 参考文献 (本ノートで参照)

- Feder, Russell, Lothe, Pound, *Adv. Phys.* 15, 111 (1966).
- Kantrowitz, *J. Chem. Phys.* 19, 1097 (1951).
- Wedekind, Hyvärinen, Brus, Reguera, *PRL* 101, 125703 (2008); arXiv:0804.1516.
- Hrubý, Vinš, Mareš, Hykl, Kalová, *J. Phys. Chem. Lett.* 5, 425 (2014); Vinš et al., *J. Phys. Chem. B* 119, 5567 (2015); Vinš et al., *J. Phys. Chem. Lett.* 11, 4443 (2020).
- Wilhelmsen, Bedeaux, Reguera, *J. Chem. Phys.* 142, 171103 (2015) (Tolman length of water).
- Daum & Gyarmathy, *AIAA J.* 6, 458 (1968); Daum, *AIAA J.* 1, 1043 (1963).
- Hansen & Nothwang, NACA TN-2690 (1952); Grossir & Rambaud, AIAA 2014-1153; Lin et al. 2014 (papers/).
