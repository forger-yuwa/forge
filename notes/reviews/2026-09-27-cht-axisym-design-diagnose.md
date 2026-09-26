# codex 諮問 (diagnose): cht-axisym-design

- **brief**: [`notes/reviews/briefs/2026-09-27-cht-axisym-design.md`](../../notes/reviews/briefs/2026-09-27-cht-axisym-design.md)
- **plan**: [`plans/active/boundary-cht-axisymmetric-fem2d.md`](../../plans/active/boundary-cht-axisymmetric-fem2d.md)
- **date**: 2026-09-27
- **commit**: `174180e4` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 3.4 min, rc=0
- **結論**: **現行の軸対称 CHT 拒否を維持し、上記の欠落を §4・§6 に反映して検証契約を確定することを次の一手とする。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（**Critical なし。現稿のままの実装着手は推奨しない**）

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| §4.1–4.2 の弱形式、積分済み荷重、`D_f` | **採用** | 要素内一定の `k_e` という前提では、重心半径による剛性・提示された Robin 行列は正しい。`Q_f` の直接受け渡しも維持する。現実装も荷重を直接渡し、`D_f` を反復の前処理に使っている（[conjugateWall.cpp:715](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:715)、[同:765](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:765)）。`D_f` と G-if の面積は固体側 `A_i^r` でよい。 |
| §4.4 の軸上判定 | **却下・Major** | `A_i^r≤0` は「軸に接する界面」の判定にならない。[plan:64](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:64) の式では、`L=1, r_i=0, r_j=1` の軸上端点でも `A_i^r=1/6>0`。**界面節点の半径を直接検査**し、軸上端点だけの負例と軸上辺の負例を分ける。剛性単独は正定値ではなく半正定値であり、Robin による拘束は連結成分ごとに必要。 |
| §4.3 の診断面積変更と §8-1 | **修正して採用・Major** | `r` 重みは `axisymMethod: 0` に限られ、`1` の幾何は平面のまま（[variables.cpp:530](/home/sano/work/forge-cht/solver_density_cuda/variables.cpp:530)、[同:593](/home/sano/work/forge-cht/solver_density_cuda/variables.cpp:593)）。CHT の拒否ガードは診断単独には届かない。**分母変更を method 0 に限定**し、未検証の実効熱量は有効な値として出さないこと。§8-1 の影響調査は、この変更と同じ plan で扱う。`interfaceDiag` 全項目が誤りとは限らず、`iface_q_eff`・`iface_q_eff_raw` とその利用先を調べる。 |
| 固体出力と半径が変わる辺の検証 | **追加・Major** | 現在の `q_hole` は `Lh(T−Tc)/2` の平面積分（[solidFem2d.cpp:301](/home/sano/work/forge-cht/solver_density_cuda/conjugate/solidFem2d.cpp:301)）。組立てだけ直しても出力収支は直らない。**Robin の出力熱量にも consistent 行列と同じ積分を使う**。また、円筒側面は半径一定なので、半径が変わる辺の誤実装を検出できない。独立した辺積分の厳密値との照合を V-ax1 に追加する。C++ と Python の相互一致だけでは共通の誤りを排除できない。 |
| V-ax2 の準定常・収束条件 | **却下・Major** | `check_quasisteady.py` は入力系列の平均で規格化する（[check_quasisteady.py:280](/home/sano/work/forge-cht/solver_density_cuda/tools/check_quasisteady.py:280)）。**絶対温度に `0.001` を指定しても「温度降下の 0.1%」にはならない**。下記 A/B で確認した。温度系列の抽出法・基準温度・判定窓を登録し、温度降下基準の許容に換算する。`NOT CONVERGED` の一律免除も不可。機械精度床に限る例外条件と、内部温度場・全保存量の停滞確認を明記する。 |
| §6 の事前登録の不足 | **追加・Major** | `fem2d` の G-if は `tol_solid` 未登録だと `REFUSED`（[check_cht_interface.py:244](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_interface.py:244)）。軸対称の **[W/rad] に対応する値**を登録する。併せて、半径・軸方向長さ・物性・温度・`q_*`・G-cons の絶対許容と床・流体格子感度を固定する。[plan:124](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:124) の `Q'` は単位軸長あたりなので、節点荷重総和との比較には軸方向長さを掛ける。 |
| V-ax3 の平面回帰 | **修正して採用・Minor** | 「同じ式」と「ビット同一」は別。既存の平面演算順序を残し、固体単体の組立て・求解は厳密一致で検査する。GPU 連成のノイズ許容は比較量・ノルム・反復数・基準バイナリを事前登録し、差が出てから許容を決めない（[plan:132](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:132)）。 |

補足すると、流体面積と固体面積の違いを「`D_f` と報告値だけの問題」とするのは不十分です。半径が線形に変わる辺では、半辺積分と FEM の形状関数積分で荷重分布も異なります。例えば `L=1, r_a=1, r_b=2` なら、

- 流体の半辺積分：`(0.625, 0.875)`
- 固体の集中量：`(2/3, 5/6)`

総和は一致しますが、一様な物理熱流束でも `Q_f/A_i^r` は両端で `0.9375q, 1.05q` になります。**保存性は保てても局所整合性は別途検証が必要**です。これを理由に荷重を面積で掛け直す変更は推奨しません。

結論: **現行の軸対称 CHT 拒否を維持し、上記の欠落を §4・§6 に反映して検証契約を確定することを次の一手とする。**

第 1 仮説: 現稿の V-ax2 は、要求する温度精度に対して未定常な系列を `STEADY` として受け入れ得る。 確度: **高**  
  根拠: `check_quasisteady.py:280` の規格化と、同関数を読み取り専用で実行した下記 A/B。  
  反証条件: 実際の評価器が既に温度降下基準の許容を実装し、その許容を超えるドリフトを必ず不合格にすることをコードと負例で示す。

第 2 仮説: 同心円環のみでは、半径が変わる界面・Robin 辺の局所誤差が残っても合格する。 確度: **高**。円環の境界半径は一定であり、その方向の検出力がない。実際のノズル壁温への影響量は未確認。

判別 A/B: **入力温度の基準だけを変更する、流体計算不要の判定器試験**。100 点の線形系列 `Tw=1000→1000.5 K`、`Tc=990 K`、`tail=0.4`、`drift=osc=0.001` を固定し、A は `Tw`、B は `Tw−Tc` を入力する。  
  → **A が `STEADY`、B が不合格なら第 1 仮説を支持**。A も温度降下基準で不合格にするなら、第 1 仮説の前提を再検証する。  
  **今回の実行結果：A は `STEADY`、B は `DRIFTING`。** 末尾窓の変化は `0.19697 K`、10 K の温度降下に対して **1.97%**。これは run の収束結果ではなく、判定器の合成入力試験である。

やらない方がよいこと: 円筒殻の解析解に合った段階で拒否を解除すること、`D_f` を調整して面積定義の不整合を解消したと判断すること、絶対温度の `STEADY` を温度降下基準の精度保証に使うこと。

呼び出し側の前提への異議:

- 「面重心の半径を掛けるのは一次近似」は一般には不正確。直線面上の線形な半径なら重心積分は厳密で、FEM との差は積分領域・重み関数の違いである。
- `case/58.conjugate_slot/run_0012_v6p_df5_i50_nlog/` の保存済み判定は、流体 **`NOT CONVERGED`**、全界面 G-if も **`NOT CONVERGED`**。後者は step `96000–99950` で `res_rel=2.3211e-3`。帯限定の `PASS` と区別する。動作回帰の入力には使えるが、全域収束済みの基準とは扱えない。
- 全モードの CHT 拒否は妥当。ただし、診断単独の正しさを保証するガードではない。

不足情報: 軸対称の数値 run、拒否確認の保存先・実行バイナリ識別情報、V-ax2 の具体的入力・評価器が未提示。ローカルには `case/59` がなく、`case/52` の指定 `run_0003` も確認できなかった。軸対称の実測精度・収束性は判定していない。**ファイル変更・`forge` 起動は実施せず、plan 未反映。反映先は本 plan §4・§6・§8-1。**
