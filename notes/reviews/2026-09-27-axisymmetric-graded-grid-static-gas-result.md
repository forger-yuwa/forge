# codex レビュー: axisymmetric-graded-grid-static-gas (result)

- **plan**: [`plans/active/axisymmetric-graded-grid-static-gas.md`](../../plans/active/axisymmetric-graded-grid-static-gas.md)
- **stage**: `result` (diff base `94085f33`)
- **date**: 2026-09-27
- **commit**: `77e270c8` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 5.2 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M0/m3
- **focus**: 限定した完了範囲 (§1) で閉じてよいか。コード変更は評価器 (case/62/eval_static_hold.py) のみでソルバは不変。判定の根拠 run は case/62 run_0016〜0020
- **extra**: `case/62.conjugate_disk/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

§1で限定した範囲なら閉じてよい。主要な合格数値・連結残差判定・準定常判定は再現できた。  
ソルバ修正や追加計算は不要だが、移動前に以下の文書・証跡を整える必要がある。

1. **Minor — 撤回済みの判断と古い未完了表示が残っている。**  
   根拠: [plan:39](/home/sano/work/forge-cht/plans/active/axisymmetric-graded-grid-static-gas.md:39) は依然「棄却した仮説: 入力メッシュの float32 閉包欠損」と断定している。A/Bが示すのは今回の大振幅過渡への支配性が支持されなかったことであり、欠損自体の不存在ではない。また、[case README:113](/home/sano/work/forge-cht/case/62.conjugate_disk/README.md:113) の「V-ax2bは未完了」は同ファイル119行の合格結果と衝突する。  
   **対案:** 閉包欠損の表現を修正し、旧未完了記述には「当時の判定・後続結果で解消」と明記する。発注元planの[217行](/home/sano/work/forge-cht/plans/accepted/boundary-cht-axisymmetric-fem2d.md:217)に残る「最終到達状態は未確認」も、非連成の到達確認と連成開始後の未検証を区別する。

2. **Minor — 変更範囲と回帰免除理由の記述が実装と一致しない。**  
   根拠: [plan:23](/home/sano/work/forge-cht/plans/active/axisymmetric-graded-grid-static-gas.md:23) は「コード変更は無い」、[107行](/home/sano/work/forge-cht/plans/active/axisymmetric-graded-grid-static-gas.md:107)は「コード修正はしない」だが、実際には[評価器:82](/home/sano/work/forge-cht/case/62.conjugate_disk/eval_static_hold.py:82)と負例テストが変更されている。§7もソルバ修正箇所が「未定」のまま。  
   **対案:** 「ソルバ・既定値は不変。評価器のみ符号つき比較へ修正」と統一する。ソルバ回帰の再実行免除と、評価器の検証実施を分けて記録する。[methods/boundary.md:500](/home/sano/work/forge-cht/methods/boundary.md:500)の「実証済み」にも、登録されたFP64円板条件という限定を添える。

3. **Minor — 初期500 stepの改善率を独立に再集計できる証跡がローカルにない。**  
   根拠: [plan:78](/home/sano/work/forge-cht/plans/active/axisymmetric-graded-grid-static-gas.md:78) の `B/A = 1.5e-5 / 8.6e-9` は250–500 stepの最大値比較だが、`run_0016_ic_{A_uniform,B_conduct}` にある流体場は各 `res_0.h5` と `res_500.h5` のみで、比較量の時系列CSVもない。終端値は確認できたが、この窓のピーク比は再検証できなかった。  
   **対案:** AWSで集計した比較量の時系列CSVと集計コマンドを回収し、run一覧から参照する。再計算は不要。後続の静止保持合格を覆す問題ではない。

確認できた主要結果は次のとおり。パスはすべて `case/62.conjugate_disk/` 配下で、恒久索引は[READMEの計算run一覧](/home/sano/work/forge-cht/case/62.conjugate_disk/README.md:102)。

| 根拠run・判定窓 | 再確認した実測 | 判定 |
|---|---|---|
| `run_0018_hold_B_ext40k`、累積30000–40000 | max\|U\| **5.17520e-4 m/s**、壁熱流束最大誤差 **0.00565277%**、相対圧力差 **2.87620e-7** | 連結残差 **PASS**、72系列 **ALL STEADY** |
| `run_0020_hold_A_ext100k`、累積90000–100000 | max\|U\| **3.63755e-4 m/s**、壁熱流束最大誤差 **0.00278025%**、相対圧力差 **2.77942e-7** | 連結残差 **PASS**、72系列 **ALL STEADY** |
| `run_0019_disk_r32g_condic` | 壁温誤差 **0.0857738%**、G-cons **8.42781e-7%**、連成保存 **1.65885e-10%** | 評価 **PASS**、残差 **PASS**、66系列 **ALL STEADY** |

連結残差は元CSV全行からメモリ上で再構成し、`check_convergence.py` の判定処理で再確認した。Bは8.0–8.3桁、Aは13.3–14.3桁低下。Bの継続区間単独の `NOT CONVERGED` を隠して合格扱いした形ではない。`run_0017` のBが `qerr_cj: DRIFTING` で判別未了だったことも再現した。

指定diffにソルバ・既定値の変更はない。評価器の熱流束符号・単位・体積重みは登録式と整合し、書込みをメモリへ置き換えた評価器本体の実行で、正常符号は `PASS`、両壁反転・片側反転は `FAIL` を確認した。A/B終端質量の約30.1%差も記載どおりで、同一質量の平衡解を保証対象から外した限定は必要である。

**推奨は、上記1→2→3を処理し、§5.1と§6.1に対応を記録してから、現在の限定範囲のまま `accepted` へ移すこと。** 減衰機構・既知slip欠陥・一様ICからの連成開始は未解決／対象外として残す。ファイルは変更しておらず、指摘・推奨はplan未反映。

指摘数: Critical 0 / Major 0 / Minor 3
