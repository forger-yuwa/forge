# codex 諮問 (diagnose): conjugate-benchmarks-design

- **brief**: [`notes/reviews/briefs/2026-09-29-conjugate-benchmarks-design.md`](../../notes/reviews/briefs/2026-09-29-conjugate-benchmarks-design.md)
- **plan**: [`plans/active/boundary-cht-conjugate-benchmarks.md`](../../plans/active/boundary-cht-conjugate-benchmarks.md)
- **date**: 2026-09-29
- **commit**: `a27fc754` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 5.1 min, rc=0
- **結論**: **A を有限温度 Robin 加熱へ設計変更し、参照解にも同一境界条件を課す方針を、`plan` §3・§4・§6 に反映することを推奨する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **一様熱流束の Robin 代用：却下** | [`conjugateWall.cpp:985`](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:985) は固体温度を **`min(T_c)−20` 以上、流体全温の最大値＋20 以下**に制限する。`T_c≈10⁷ K` では許容区間が空になり、正常な温度でも停止する。**ソルバ無変更の範囲を維持するなら、A を有限温度・有限 `h` の Robin 加熱ベンチマークに変更する**。参照解にも同じ Robin 条件を課し、一様熱流束とは呼ばない。 |
| **Major** | **C の主判定：案 (b) を採用。ただし速度だけの移植は不十分** | 現案は定密度・散逸なしの式と Blasius を使う（[`plan:43`](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:43)、[`plan:61`](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:61)）。一方、既存の温度再現器は **`ρ,u,v,p` と粘性散逸・圧力仕事**を用いる（[`temp_reproduce.py:52`](/home/sano/work/forge-cht/case/63.graetz_cht/temp_reproduce.py:52)）。主判定はこれらを固定した独立な共役温度計算とし、**界面温度・熱流束は参照解の未知数として解く**。forge の界面温度を境界条件として与えてはいけない。主張は「与えられた流れに対する伝熱・連成の検証」に限定する。Blasius と窓の設定だけでは、上流の差が伝熱場へ伝わる影響を除けない。 |
| **Major** | **対照差し引きで有限 M を処理：主判定として却下** | 前身の正式な温度診断は **`VERDICT: 判定不能`**（[`TEMP_REPRODUCE.txt:9`](/home/sano/work/forge-cht/case/63.graetz_cht/run_0014_g2_dT0_r32/TEMP_REPRODUCE.txt:9)）。今回読み直した最終場でも対照の `x=0` 内部32節点の温度幅は **0.138389 K**だが、これは原因の確定ではない。加熱で `ρu` と源項が変われば、差し引いても同じ線形温度方程式にはならない。**絶対温度・絶対熱流束を、源項を含む参照解と比較する**。M の引下げを最初の対策にはしない。低 M 化だけでは `ΔT/T` による密度変化を除けない。 |
| **Major** | **A の条件：候補として採用、妥当性確定は要再検証** | `r_o/R=2`、`k_s/k_f=10,100`、`Pe_D=720` は条件として定義できるが、**`L_h/R` が未指定**で予熱割合は決まらない（[`plan:50`](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:50)）。有限 Robin に変更後は **`h r_o/k_s`、`L_h/R`、上下流長/R、`(T_c−T_in)/T_in`**も登録する。同じ加熱条件で伝導比を比較し、各条件で温度上昇を10 Kに合わせ直して比較を混ぜない。 |
| **Major** | **「固体が効く」の判定：現案を却下** | 等温固体との差だけでは、熱抵抗と軸方向伝導を分離できない（[`plan:90`](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:90)）。また C の `Br` は **`b/L` の指定なしでは軸方向伝導の強さを決めない**（[`plan:58`](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:58)）。固体方程式を無次元化すると軸方向項と厚さ方向項の係数比は `(b/L)²`。**厚さ方向温度差と、固体断面を通る軸方向熱量を別々に登録・評価する**。軸方向熱量の空間変化が界面熱収支へ与える寄与まで確認し、各効果が測定不確かさの5倍以上かつ対応する主判定許容幅を超える条件を選ぶ。 |
| **Major** | **§6 の許容差：目標値として採用、合否基準としては要再検証** | 現案は参照格子2水準以上と「差が減る」だけで、不確かさの上限がない（[`plan:45`](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:45)、[`plan:87`](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-benchmarks.md:87)）。**参照格子は最低3水準**とし、参照離散化・領域切断・写像・反復の不確かさの合計を許容幅の1/3以下にする。合否は `観測差＋不確かさ ≤ 許容幅`。満たせなければ判定不能とする。C の「q の3%」は分母を明記する。 |
| **Minor** | **Robin の数値的不安定性：`T_c` が大きいだけという説明は却下** | 固体行列・荷重は `double` で、`T_c` は荷重の `h*T_c` に入り、行列には `h` が入る（[`solidFem2d.cpp:187`](/home/sano/work/forge-cht/solver_density_cuda/conjugate/solidFem2d.cpp:187)）。連成行列には界面 `D_f` も加わる（同 [`:217`](/home/sano/work/forge-cht/solver_density_cuda/conjugate/solidFem2d.cpp:217)）。小さい `h` による固体単独系の悪条件化は懸念になるが、**実際の連成行列の悪条件化・収束停滞は未確認**。まず確実に存在する温度ガードの不整合を扱う。`Df_scale` を上げてもこの拒否条件は直らない。 |

結論: **A を有限温度 Robin 加熱へ設計変更し、参照解にも同一境界条件を課す方針を、`plan` §3・§4・§6 に反映することを推奨する。**

第 1 仮説: 提案された巨大 `T_c` の Robin 代用は、数値発散する前に、冷却を前提とした固体温度ガードによって拒否される。 **確度: 高**
  
  根拠: [`conjugateWall.cpp:999`](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:999) の判定式を、`case/63.graetz_cht/run_0016_g2_dT10_r32/` の保存場で読み取り専用再現した。上限は **329.999390 K**。巨大 `T_c` 条件では下限が **10,000,285 K**となり、**固体2115節点すべてが拒否**された。これは温度ガードの再現であり、共役計算の結果ではない。
  
  反証条件: 実際に使用するバイナリでこのガードが存在しない、または実際の `ROBIN/TC` がブリーフと異なり、下限と上限が逆転しないこと。ただし後者は提案条件の再現ではない。

第 2 仮説: 解析速度・定密度・源項なし参照との差を、そのまま共役実装の誤差と誤認する可能性がある。**確度: 中**。前身の原因診断は判定不能であり、今回の誤差量は未確認。

第 3 仮説: C は厚さ方向の熱抵抗が効いても、登録窓内の軸方向伝導は検出限界以下になり得る。**確度: 中**。`b/L` 未指定のため未確認。

判別 A/B: **温度ガードの単体再現を実施済み。追加 forge 計算は0 step。** 保存した流体・固体場を固定し、独立に変えるパラメータを `h` とした。`T_c=305+1000/h` と連動させ、305 Kでの Robin 熱流束だけを1000 W/m²に揃えた。

- **A:** `h=100`、`T_c=315 K` → 下限295 K、拒否 **0/2115節点**。
- **B:** `h=10⁻⁴`、`T_c=10,000,305 K` → 許容区間が空、拒否 **2115/2115節点**。

判別条件は「A が通り B だけ落ちれば境界表現と温度ガードの不整合を支持、B も通れば第1仮説を棄却」。結果は前者。**両境界条件が伝熱問題として等価という試験ではなく、A の連成収束も保証しない。**

やらない方がよいこと: 巨大 `T_c` のまま投入し、停止メッセージに従って CFL・`Df_scale`・warmup を調整すること。ダミーの Robin 辺で `tcMin` を下げてガードをすり抜けること。等温固体との差だけで「軸方向伝導も検証済み」とすること。

呼び出し側の前提への異議: **「外側境界に Robin が実装されている」ことと「任意の加熱 Robin が現在の連成処理で受理される」ことは別**。また、上流熱割合は軸対称なら
`Q_up = R∫up q_i dx`、`Q_total = r_o∫heat h(T_c−T_o) dx`
として同じ単位で定義する必要がある。`T_i−T_b` が不確かさに近い予熱末端では局所 Nu を主判定に使わず、温度・熱流束を使う。

不足情報: A の加熱長と下流長、C の `b/L`・加熱温度差・背面 `h`・遠方境界位置、参照解の境界条件一式、固体効果を評価する尺度。**ファイル変更なし。以上の判断は plan 未反映。**
