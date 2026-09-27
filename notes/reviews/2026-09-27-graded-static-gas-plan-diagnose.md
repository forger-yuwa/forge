# codex 諮問 (diagnose): graded-static-gas-plan

- **brief**: [`notes/reviews/briefs/2026-09-27-graded-static-gas-plan.md`](../../notes/reviews/briefs/2026-09-27-graded-static-gas-plan.md)
- **plan**: [`plans/active/axisymmetric-graded-grid-static-gas.md`](../../plans/active/axisymmetric-graded-grid-static-gas.md)
- **date**: 2026-09-27
- **commit**: `9ad4e5ab` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 4.7 min, rc=0
- **結論**: **同じ非一様格子・BC・バイナリで、一様初期場と、両壁温度に整合した純伝導静止初期場の A/B を最初に行う。**
- **extra**: `case/62.conjugate_disk/README.md`, `notes/investigations/node-slip-tangential-density-spurious-flow.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| Major | §4 の候補順位：**要再検証** | 初期場には既に圧力非平衡がある。下記の実測から、最初は **H2 を優先**し、温度勾配と slip を残して初期場だけを替える。 |
| Major | §4(iii)「slip→等温壁で消えれば H1」：**却下** | 等温壁化は熱境界条件も変更し、`nodeWallDirichlet` は運動量を拘束するため、原因によらず擾乱を抑え得る。[境界仕様:28](/home/sano/work/forge-cht/methods/boundary.md:28)、[残差拘束:56](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:56)。対案は下記の初期場 A/B。 |
| Major | §6「P が上昇すれば整合」：**却下** | `unsteady: 0` は局所擬似時間刻みを使うため、反復中の総質量は物理時間の保存試験にならない。[setDT_d.cu:362](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/setDT_d.cu:362)。対案は圧力の**空間一様性**を静止保持の指標にし、絶対圧力は質量の扱いと分けて評価する。 |
| Minor | §3 のピーク・圧力差の記述：**修正採用** | `run_0015_hoop_A_hoop0` の全保存区間では `max|Uy|=31.07 m/s`、step 140。19.67 Pa は **`x=2.5 mm` 列の半径方向圧力差**であり、領域全体の圧力差ではない。区間と抽出位置を明記する。 |

結論: **同じ非一様格子・BC・バイナリで、一様初期場と、両壁温度に整合した純伝導静止初期場の A/B を最初に行う。**

第 1 仮説: **H2――加熱壁と一様初期場の非平衡が擾乱を与え、非一様格子上の更新で増幅される。** 確度: **中**

  根拠:
  
  - 温度ピンは密度を保持して `P=ρRTw` を設定する。[nodeWallDirichlet_d.cu:80](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:80)。
  - `case/62.conjugate_disk/run_0015_hoop_A_hoop0/` の保存場を直接集計した結果は以下。座標は `(x,r)`、単位 mm。

    | step | 観測 |
    |---:|---|
    | 0 | 全速度ゼロだが、P は **1013.25–1091.19 Pa**。静止平衡ではない |
    | 10 | `max|U|=5.82 m/s`、位置 **(0.3125, 5.1566)**。加熱壁に隣接する内部節点 |
    | 30 | 最大位置が **(0.3125, 6.1885)** に移る |
    | 50 | 同位置で **46.01 m/s**。slip 列の最大 `|Ux|` は **1.42 m/s** |

  - 増幅後の最大値は内部の細格子側にある。ただし、最初の非ゼロ保存場が step 10 なので、**slip 起点を否定する証拠にはならない**。
  - 再実行した `check_convergence.py` は **`NOT CONVERGED`**。`check_quasisteady.py` は step 0–500 の `machmax,pmax` に **`DRIFTING`**。同ツールの `classify_series` に渡した `Umax`・`Uymax`・領域圧力差・壁市松振幅もすべて **`DRIFTING`**。上記は過渡の観測である。

  反証条件: **圧力一様・壁温整合の純伝導初期場でも、500 step 内に同程度の半径速度と壁市松が成長するなら、「大きな加熱起動擾乱が必要」という仮説を棄却する。** H2 全体、特に陰的更新そのものの問題まで棄却するわけではない。

第 2・第 3 仮説:
- **H1、確度低〜中:** slip と接線密度勾配の欠陥。既知ケースは純伝導初期場でも発生したとの記録があり、今回も候補に残る。ただし同一原因という証拠はない。[調査記録:9](/home/sano/work/forge-cht/notes/investigations/node-slip-tangential-density-spurious-flow.md:9)。
- **H3、確度低・未確認:** 閉包補正以外の軸対称幾何・離散化誤差。既存 A/B で下げられるのは閉包欠損の支配性であり、幾何処理全体ではない。

判別 A/B:

**変更因子は初期保存量場だけ。**

- **A:** 現行の一様 325 K・静止初期場。
- **B:** 同じ節点上に  
  \[
  T_i=350-25x_i/H,\quad U_i=0,\quad
  \rho_i=P_*/(RT_i),\quad (\rho E)_i=P_*/(\gamma-1)
  \]
  を与える。`R=287 J/(kg K)`。初期総質量を揃えるため、
  \[
  P_*=\frac{M_A}{\sum_i V_i/(RT_i)},\qquad V_i=A_i\bar r_i
  \]
  とする。入力の双対重心・体積から計算すると **約1051.736 Pa**。ノード半径で体積を作り直さない。
- 別々の新規 `run_NNNN_*` に作り、**各500 step**。`warmup: 20000`、slip、次数、CFL、`hoopAreaFromClosure: 0` は共通。出力は初期10 stepを毎step、その後10 stepごと。
- 見る量は全保存量残差、slip 列／内部別の速度、領域全体と `x=H/2` 列の圧力差、既登録の壁市松振幅。比較窓は **step 250–500** と固定する。

**事前判定:**

- A が再現し、B の `max|Uy|` と壁市松振幅がともに A の **1/10以下**なら、起動非平衡が大振幅化を支配する H2 を支持する。静止保持合格とはしない。
- B でも両量がともに A の **1/2以上**なら、大きな起動擾乱を必要とする仮説を棄却し、平衡場自体を壊す作用素を調べる。H1/H3 の選別はこの結果だけではできない。
- 中間・指標不一致・A 未再現は判別不能。短時間で発現しないことを、欠陥不存在に読み替えない。

**§6 の閾値は、次のように登録することを推奨する。**

- `max|U|≤1e-3 m/s` は暫定採用。この条件では初期密度を使った `Pe_H=ρcpUH/k≈0.0023` となり、伝導への影響を小さくする目安になる。ただし誤差保証ではないため、非連成区間で **全壁節点の熱流束誤差≤0.5%（基準120.5 W/m²）**も直接確認する。
- 圧力は **`(Pmax−Pmin)/〈P〉V≤1e-4`**を提案値として事前登録する。「上昇したか」は使わない。
- 合否には、同一設定区間の残差 **`PASS`** と対象量の **`STEADY`**を併用する。現行の一様格子 step 20000 の値だけでは閾値を正当化しない。

やらない方がよいこと: **slip の修正に直行すること、両壁325 Kで静止しただけでH3を除外すること、warmup延長だけで済ませること。** 無加熱試験は、温度勾配と組み合わさって現れる幾何欠陥を除外できない。[plan:54](/home/sano/work/forge-cht/plans/active/axisymmetric-graded-grid-static-gas.md:54) の二分法は強すぎる。

呼び出し側の前提への異議:
- 閉包 A/B は入力メッシュのハッシュ・バイナリのハッシュが同じで、YAML差分も該当キーだけだった。**閉包補正で大振幅が解消しないという判断は採用**し、この A/B は繰り返さない。
- 「加熱だから圧力が上がる」は、質量を保存する物理過程についての期待。固定質量・線形温度場なら期待圧力は約 **1051.74 Pa**であり、一様側の1027 Paも、その期待を満たす基準ではない。
- run の所在は [case README の一覧](/home/sano/work/forge-cht/case/62.conjugate_disk/README.md:101)。**本回答は plan 未反映**。§3・§4・§6への反映は呼び出し側で行う。

不足情報: **step 1–9 の場、初期の局所残差・流束内訳、実行ビルドの未コミット変更1件の内容。** 現資料では最初の発生位置と増幅するコード経路までは確定できない。
