# codex 諮問 (diagnose): sern-b1b-regression

- **brief**: [`notes/reviews/briefs/2026-10-05-sern-b1b-regression.md`](../../notes/reviews/briefs/2026-10-05-sern-b1b-regression.md)
- **plan**: [`plans/active/tooling-sern-mesh-blocking.md`](../../plans/active/tooling-sern-mesh-blocking.md)
- **date**: 2026-10-05
- **commit**: `e734ad8a` (feature/sern-design)
- **codex**: effort `high`, 4.2 min, rc=0
- **結論**: **B1b 修正を採用し、C_M の登録判定は保持したまま、追加 CFD より先に既存2系列のモーメントを基準点の寄与に分解する。**
- **extra**: `notes/investigations/2026-10-05-sern-b1b/B1B_VERDICT.txt`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 判断・採否 | 根拠と対案 |
|---|---|---|
| **Major** | **B1b の修正は採用。ただし「全4量への影響は無視できる」は未成立** | 内部面は節点から double で面重心を計算し、修正後の境界半割面も同じ構成になっている（[gmshReader.hpp:1863](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/gmshReader.hpp:1863)、[同:1999](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/gmshReader.hpp:1999)）。幾何の不整合を直す変更として妥当。**修正の採用と旧結果の同等性認定を分ける。** |
| **Major** | **C_M は「影響無視の判定不能」。回帰悪化の確定でも、合格でもない** | [判定原本:6](/home/sano/work/forge-sern-design/notes/investigations/2026-10-05-sern-b1b/B1B_VERDICT.txt:6) の Δ＝−4.20×10⁻⁵、U＝4.81×10⁻⁴から、Δ±U＝**[−5.23×10⁻⁴, +4.39×10⁻⁴]**。許容帯をまたいでいる。原本の `VERDICT: 閾値超過あり → 記録して諮る` を保持する。U は経験的な比較幅であり、統計的信頼区間ではない。 |
| **Major** | **B4 の受入経路へ進む前に、閉性 FAIL の続行を廃止する** | [run_junction_model.py:76](/home/sano/work/forge-sern-design/case/46.sern_design/cad/run_junction_model.py:76) は現在も警告だけで続行する。**非ゼロ終了なら記録して停止**させる。修正版の接続模型・生産格子について、閉性・正体積・品質の個別判定と変換器の来歴を保存し、今回の再現例を回帰検査に固定する。 |
| **Minor** | **2D に同じ修正を広げる案は不要** | 内部面の中点は double、境界半割面も節点差から double で再構成し、重心は `(3N+O)/4` で求める（[gmshReader.hpp:1440](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/gmshReader.hpp:1440)、[同:1596](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/gmshReader.hpp:1596)）。今回、修正前 commit と現行の **2D 関数部分が完全一致**し、保存済み primal 面重心の参照がないことを確認した。「同じ欠陥」の最小確認はこれで済む。2D 幾何全般の無欠陥保証ではない。 |
| **Minor** | **仕様・残作業表を実装範囲に合わせて訂正する** | [methods/discretization.md:304](/home/sano/work/forge-sern-design/methods/discretization.md:304) は保存済み面重心の流用を記載し、[B1b 行:774](/home/sano/work/forge-sern-design/plans/active/tooling-sern-mesh-blocking.md:774) は座標読込みまで double 化する旧案のまま。**「3D 境界面重心の再計算で内部側と整合」へ更新**し、幾何修正の完了と C_M の未認定を別記する。 |

結論: **B1b 修正を採用し、C_M の登録判定は保持したまま、追加 CFD より先に既存2系列のモーメントを基準点の寄与に分解する。**

第 1 仮説: **C_M の差と比較幅は、長いモーメント腕で増幅された C_L の差・変動が主成分である。** 確度: **中（平均差は支持、変動の共変動は未確認）**

根拠: 対象は `case/46.sern_design/run_1067_b1b_old/` と `case/46.sern_design/run_1068_b1b_new/`。集約原本では U が D の約92%を占め、U_M/U_L≈22。基準点は x_ref＝−20H（[問題定義:42](/home/sano/work/forge-sern-design/case/46.sern_design/problem_3d_prod_m6on_wallres_lswx08.yaml:42)）。[積分式:259](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:259) の符号規約では、同じ場について厳密に

`C_M(x_ref=0) = C_M(x_ref=−20H) + 20 C_L`

となる。報告された平均差から計算すると、原点まわりの差は **−1.8×10⁻⁶**で、元の差の約4.3%。ただし、U の分解には時系列の共変動が必要で、集約値だけでは計算できない。

反証条件: 同じ保存時刻の全系列で分解しても、下記の「平均差・比較幅をともに90%以上縮小する」という診断条件を満たさないこと。

第 2 仮説: **長腕の寄与を除いても、圧力分布の変化によるモーメント差・変動が残る。** 確度: **低・未確認**。第 3 仮説は置かない。

判別 A/B: **変更は後処理の x_ref だけ。CFD は0 step。**

- **A**＝登録済みの x_ref＝−20H、**B**＝診断専用の x_ref＝0。両 run の既存20000 step系列を使い、前後10000 step、同じ標本集合・同じ a、d、U の定義で集計する。丸め済み表示値ではなく、元の `force_history.csv` を用いる。
- **結果A**：B の |Δ|≤4.2×10⁻⁶、U≤4.81×10⁻⁵をともに満たす → 長腕による増幅を主因として支持する。
- **結果B**：いずれかを満たさない → 「差と変動の双方が90%以上長腕で説明できる」という第1仮説を棄却する。
- **これは原因を分ける診断であり、登録済み C_M の合格判定には使わない。** どちらの結果でも元の `D=5.23×10⁻⁴` は保持する。

やらない方がよいこと: **4.6%の超過だから合格へ丸めること、U を削ること、基準点を変更して受入条件を通すこと、振幅が減衰する証拠なしに通るまで延長すること。** また、今回の同等性未認定を理由に閉性不良の旧変換器へ戻さない。B4 の開発は進めてよいが、新形状の生産受入は B5/B6 で取り直す。

呼び出し側の前提への異議: **`GATES PASS` は残差収束の証明ではない。** 評価器は `require_residual_pass=False` が既定である（[runner_sern3d.py:340](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:340)）。今回確認できたのは集約原本の `GATES PASS`・4量 `STEADY` という記録までで、「両 run が収束した」とは認定しない。また、R7a の同型例は判定不能という扱いの先例であり、今回の変動原因を確定する証拠ではない。

不足情報: 対象2 run の生履歴・個別 VERDICT は本 checkout にない。必要なのは `force_history.csv`、全残差の判定と区間、準定常判定の設定、旧新格子の検査原本、バイナリ・実効 config・全保存量の初期コピー照合記録。`roY*` を作成して照合する修正はコードで確認したが、再実行時の成功記録は未確認。恒久索引は [case README:383](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:383)。

**ファイル変更・forge 起動なし。plan 未反映。** 呼び出し側の反映先は `tooling-sern-mesh-blocking.md` §5.1 B1b・R1／§6。
