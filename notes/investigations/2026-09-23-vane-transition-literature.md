# 冷却翼 (C3X / Mark II) に遷移モデルを当てた公表例 — やり方の調査 (2026-09-23)

plan [`turbulence-transition-lm2009`](../../plans/active/turbulence-transition-lm2009.md) §6.3/§6.4 で見えた 2 つの症状
(**負圧面の遷移が実測より遅く急** / **正圧面が後縁まで層流のまま**) が、文献でも同じように出ているのか、出ているなら
どう処理しているのかを調べた。**先行調査 [`cht-validation-case-survey.md`](cht-validation-case-survey.md) (2026-09-19) は
「どの試験データを使うか」の棚卸しで、遷移モデルの扱いは調べていない** (plan §3 の「文献に適用例が多い」は根拠を示さずに書いた一行だった)。

## 一次資料 (本文を読んだもの)

**Lin, Kusterer, Haj Ayed (B&B-AGEMA), Bohn (RWTH Aachen), Sugimoto (兵庫県立大),
"Conjugate Heat Transfer Analysis of Convection-cooled Turbine Vanes Using γ-Reθ Transition Model",
Int. J. Gas Turbine, Propulsion and Power Systems, 6(3):9–15, 2014** (J-STAGE 公開・無料,
<https://www.jstage.jst.go.jp/article/jgpp/6/3/6_9/_article>)。Mark II 試験番号 **5411** (超音速出口)、STAR-CCM+ と自社コード CHTflow、
2D/3D、$y^+<1$、3D は 286 万セル。**我々と同じ試験・同じ規約** ($h$ の正規化 $H_0$ = 1135 W/m²K、$T_{ref}$ = 811 K) を使っている。

### 1. 入口乱流条件の与え方 — **我々と同じ「粘性比 10」**

> "After Hylton et al. the inlet turbulence intensity was 6.5 % in the experiment, which has also been used in the numerical
> investigations. **Instead of turbulence length scale the turbulent viscosity ratio has been applied** as another inlet boundary
> condition for turbulence conditions. **The inlet turbulent viscosity ratio has been set at 10.** In the study by Mansour et al.
> the inlet turbulent viscosity ratio has **insignificant impact on the heat transfer coefficient distribution prediction at
> suction side** in the case Mark II."

- 報告が与えない長さスケールの代わりに**粘性比を使うのは標準的な逃げ方**で、**値も 10** (我々の生産設定と同じ)。
- 「粘性比の影響は小さい」という引用は**負圧面に限った話**。我々の実測 (plan §6.3) でも負圧面の遷移位置は粘性比 10/30/100 で
  $s/S$ 0.32〜0.34 とほとんど動かず、**動くのは正圧面** (−33 → +25 %) だった。**文献の主張と矛盾しない**。
- ただし彼らは粘性比を振っていない (Mansour らの結論を引いている)。**粘性比を振った感度は我々の方が広い**。

### 2. 負圧面の遷移位置 — Mark II では「実験と同じ位置」

> "By the calculation with the SST-γ-Reθ model, the turbulence is firstly produced at position by **X/L = 0.45, which is the same
> as in the experiment**." / 結論: "only the SST-γ-Reθ model tends to predict a right onset location of the transition"。

- Mark II は負圧面の**衝撃波で遷移が起きる**ので、モデルの急な立ち上がりがそのまま合う。我々の Mark II も同じ
  (負圧面層流域 +82 % → −5.5 %、rms 95 → 9.3 %、`case/54…/run_0037_lm_1um`)。
- **遷移直後の再付着点で熱流束を過大に出すのも同じ**: 彼らは温度で **+2.5 %** が最大誤差 (熱流束で見れば数十 %)。
  我々は負圧面遷移後の $h$ が +20 % (Mark II) / +2〜6 % (C3X)、rms は悪化 (12.8 → 30〜36 %)。
- C3X (亜音速出口・衝撃なし) で遷移が遅い件は、この論文は扱っていない (Mark II のみ)。

### 3. 正圧面が層流のまま — **同じ症状が出ていて、モデル定数を振って対処している**

> "Although the SST-γ-Reθ model can predict a quite well transition at the suction side, **the temperature distribution at the
> pressure side is even worse predicted than the Realizable k-ε model**. … **the boundary layer at the pressure side, near the
> trailing edge, is still laminar in the calculation with SST-γ-Reθ model**, while the other calculation results are turbulent.
> **The main reason behind this effect is that the value of $Re_{\theta t,min}$ in the transition model, which is suitable for the
> suction side, is too high for the pressure side.**"

- 対処は**輸送変数 $\tilde{Re}_{\theta t}$ の下限 $Re_{\theta t,min}$ を case ごとに振る**こと。標準の下限は 20 だが、彼らは
  **100 / 130 / 200 を与えて比較**した:
  - $Re_{\theta t,min}$ = 130 → 負圧面 (第 1 衝撃より上流) が実測とよく合う。
  - **より小さい値 → 正圧面がよく合う**。
  - 下限を上げると両面とも遷移が遅れる。
- つまり **「負圧面に合う下限」と「正圧面に合う下限」が両立しない**ことを論文自身が認めていて、片方を選んでいる。
  これは我々の「入口粘性比を上げると正圧面が合うが負圧面前縁が悪化する」というトレードオフと**同じ構造**
  (効かせ方が違うだけ: 彼らはモデル定数、我々は入口条件)。
- 結論部: "**none of them (他の著者) has investigated the influence of parameter $Re_{\theta t,min}$** on the transition at
  airfoil suction side" — つまり**他の論文は下限を既定のまま使っている**。

### 4. 実測 $h$ の素性 (我々の比較対象の性格)

> "The latter values (Hylton の $h$) have **not been obtained by measurement, but by application of an FEM tool solving the energy
> equation in the solid, while the measured surface temperatures have been used as boundary conditions**. Thus, the uncertainties
> regarding the accuracy of the experimental heat transfer coefficients are **much higher than for the experimental surface
> temperatures**."

- **我々は $h$ で合否を語っているが、文献は温度で語っている**。報告の $h$ は固体の熱方程式を実測壁温を境界値として解いた**導出量**で、
  壁温より不確かさが大きい。この論文の誤差 2.5 % / 2 % は**温度**の値で、$h$ の % とは直接比べられない。

## 二次資料 (要旨のみ)

- **"Conjugate heat transfer investigations of turbine vane based on transition models", Chinese J. Aeronautics, 2013**
  (<https://www.sciencedirect.com/science/article/pii/S1000936113000757>): Mark II に **SST / SST+AGS / SST+γ-Reθ** の 3 本。
  「遷移モデルを入れると層流域の温度・$h$ の精度が上がり、完全乱流域では 3 モデルとも同じで実験とよく合う」
  「AGS に比べて γ-Reθ は**負圧面の衝撃／境界層干渉による遷移**の予測精度を上げる」。= 我々の Mark II の結果と整合。
- Yan et al. (上の JGPP 論文が引用): 同じ Mark II を CFX/Fluent で。CFX の SST Gamma Theta が温度分布でよく一致、
  ただし**負圧面の再付着点で約 6 % 過大**。

## この調査から読み取れること (plan への反映は §5.1 #15)

1. **我々の 2 症状は既知**。正圧面が層流のまま残るのは、少なくとも Mark II では文献でも再現している (=実装の誤りを疑う材料にはならない)。
   forge/SU2 の一致 (plan §6.4) と合わせて、**「この条件では実装に固有でない」の傍証が 1 つ増えた**。
2. **文献は「合わせる」ために $Re_{\theta t,min}$ を case ごとに動かしている** (標準 20 → 130)。我々はモデル定数を触っていない。
   同じ手を採るかは**方針の決定**なので、ここでは採否を決めない (§5.1 #15 に上げる)。採るなら「どちらの面に合わせたか」を明記する種類の調整であり、
   **実験に合う値を事後に選ばない**という本 plan の立場と衝突しうる。
3. **入口粘性比 10 は文献と同じ**。したがって生産設定を変える理由は今のところ無い。
4. **比較量を $h$ から壁温へ寄せる余地がある** (報告の $h$ は導出量)。CHT を回せば壁温で比較でき、文献の 2〜2.5 % と直接並べられる。

## 出典

- Lin, Kusterer, Haj Ayed, Bohn, Sugimoto (2014), JGPP 6(3):9–15. <https://www.jstage.jst.go.jp/article/jgpp/6/3/6_9/_article>
  (PDF は J-STAGE から無料。**`papers/` には置かない** — git 追跡外でもリポジトリを太らせない方針)
- Chinese J. Aeronautics (2013). <https://www.sciencedirect.com/science/article/pii/S1000936113000757>
- NASA CR-168015 (Hylton et al., 1983). <https://ntrs.nasa.gov/citations/19830020105>
