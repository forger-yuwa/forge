# codex 諮問 (diagnose): cht-axisym-vax1c-replacement

- **brief**: [`notes/reviews/briefs/2026-09-27-cht-axisym-vax1c-replacement.md`](../../notes/reviews/briefs/2026-09-27-cht-axisym-vax1c-replacement.md)
- **plan**: [`plans/active/boundary-cht-axisymmetric-fem2d.md`](../../plans/active/boundary-cht-axisymmetric-fem2d.md)
- **date**: 2026-09-27
- **commit**: `2e2a90f9` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 2.4 min, rc=0
- **結論**: **(c1) を同方向対角限定に訂正し、(c2) の次数を二段階で判定する置き換え仕様を、実行前に plan §6・§5.1 #11 と試験コードへ登録する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 3 / Minor 1）

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **(c1)「両対角で周期解＝1D FE 解」は却下** | `rect_mesh` の交互対角は市松配置であり、節点ごとに接続する三角形が異なる（[test_solid_fem2d_axisym.py:133](/home/sano/work/forge-cht/solver_density_cuda/tools/test_solid_fem2d_axisym.py:133)）。下記の独立な要素式計算では、**交互対角は内部節点でも 1D FE 解の残差がゼロにならない**。周期化では消えない。(c1) の 1D FE 一致判定は**同方向対角だけ**にする。交互対角は (c2) の解析解への誤差で判定する。 |
| **Major** | **(c2) は修正採用、(c3) は保証範囲を限定して採用** | 旧登録は連続する二段階の次数を要求していた（[plan:158](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:158)）。(c2) は **16→32、32→64 の両方で次数 ≥1.8** にする。誤差上限は提案どおり **0.05 % / 0.1 %** を維持するが、理論から保証された値ではなく、事前に選ぶ精度要求と明記する。(c3) の次数撤去は**旧要件からの緩和**である。旧 FAIL を保持し、新試験を「指定格子での精度・単調減少の検証」と位置付けるなら採用できる。 |
| **Major** | **試験コードの判定への同期を採用・実行前必須** | 現コードは固定軸系列を「判定外」とし（[test:355](/home/sano/work/forge-cht/solver_density_cuda/tools/test_solid_fem2d_axisym.py:355)）、円筒殻の C++–Python 差も表示するだけで判定していない（[test:334](/home/sano/work/forge-cht/solver_density_cuda/tools/test_solid_fem2d_axisym.py:334)）。旧 FAIL と新試験の VERDICT を分離し、新試験では全対象の差 ≤1e-9 K・残差・誤差・次数を実際の合否条件にする。 |
| **Minor** | **V-ax2 の格子登録の具体化を採用** | 生成器は `y` を降順に並べ、各セルで同じ対角を使う（[gen_solid_strip.py:66](/home/sano/work/forge-cht/case/58.conjugate_slot/gen_solid_strip.py:66)、[同:79](/home/sano/work/forge-cht/case/58.conjugate_slot/gen_solid_strip.py:79)）。通常の座標表示では**左上→右下**である。「生成器と同じ」だけでなく、座標順・接続・各段の接線／法線分割数・界面節点の対応を登録する。提示された V-ax2 の感度試験方針は維持する。 |

結論: **(c1) を同方向対角限定に訂正し、(c2) の次数を二段階で判定する置き換え仕様を、実行前に plan §6・§5.1 #11 と試験コードへ登録する。**

第 1 仮説: **交互対角にも周期端面での 1D FE 完全一致を要求すると、正しい組立てを不合格にする。**　確度: **高**

  根拠: [test_solid_fem2d_axisym.py:150](/home/sano/work/forge-cht/solver_density_cuda/tools/test_solid_fem2d_axisym.py:150) の接続に、plan の独立 1D FE 解を代入する要素式計算を行った。新パラメータでの求解はしていない。内部節点半径を \(R\)、半径刻みを \(d\)、軸方向刻みを \(a\)、\(C=q_1r_1\) とすると、剛性行の残差は

\[
(KT_{\rm 1D})_i=
\begin{cases}
0 & \text{同方向対角},\\[2pt]
\displaystyle\pm\frac{CaRd}{3(R^2-d^2/4)} & \text{交互対角}.
\end{cases}
\]

  交互対角では、対角が集まる節点と集まらない節点で符号が反転する。**自然端面から離れた内部節点の残差なので、周期化でも解消しない。**

  反証条件: 同じ市松接続・同じ厳密な \(r\) 重み組立てで、1D FE 場を代入した内部節点残差が上式と一致せず、丸め誤差内でゼロになること。

第 2・第 3 仮説: 追加しない。旧周期試験の成功は同方向対角についての記録であり、交互対角への拡張を裏付けていない（[plan:160](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:160)）。

判別 A/B: **旧パラメータの N=8・軸6分割で、変えるのは対角規則だけ**。A＝同方向、B＝交互。周期自由度へ集約した C++ 行列に独立 1D FE 場を代入し、各1回の行列–ベクトル積で内部節点残差を確認する。求解・forge 起動は不要。予測との差を \(q_1r_1\Delta x\) で正規化し、≤1e-11 を照合条件とする。  
→ **A がゼロ・B が上式なら**「両対角で 1D FE 一致」を棄却する。**この組合せにならなければ**今回の要素式と実際の接続・組立ての対応を再調査し、置き換え試験を走らせない。

やらない方がよいこと: **(c1) の交互対角が FAIL したときに組立てを修正すること。** また、(c3) の合格を「同方向対角でも二次精度を検証した」と表現しない。旧 FAIL を新試験の PASS で上書きしない。

呼び出し側の前提への異議: **新パラメータが未求解であることは、判定条件の数学的妥当性を保証しない。** 0.05 %／0.1 % は事前の精度要求として採用できるが、合格予測としては未確認。「組立ては正しい」という旧記録も、確認した条件・作用に限定すべきである。

不足情報: 周期自由度の集約を C++ 側の直接求解へ適用する具体的経路と、線形相対残差の正規化定義が未提示。旧実測値は plan の記録として確認した範囲であり、今回再測定していない。**ファイル変更なし・新パラメータの求解なし・plan 未反映**。
