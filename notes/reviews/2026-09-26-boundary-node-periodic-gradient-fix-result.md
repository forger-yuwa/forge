# codex レビュー: boundary-node-periodic-gradient-fix (result)

- **plan**: [`plans/active/boundary-node-periodic-gradient-fix.md`](../../plans/active/boundary-node-periodic-gradient-fix.md)
- **stage**: `result` (diff base `1266aba1`)
- **date**: 2026-09-26
- **commit**: `0977ca72` (feature/sern-design)
- **codex**: effort `high`, 5.6 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m1
- **focus**: (1) dF1_inf の限定合格 (DRIFTING 0.3%/tail、漸近 0.0335、末尾全点 ≤0.05) の妥当性 (2) R2 の KE 収支基準を測定後に符号判定へ訂正したことと、それによる『旧は継ぎ目由来の注入』という帰属の妥当性 (3) R1 抽出の [解釈] 4 点と x 継ぎ目の参考列 (4) §4.2a 面フラグ planePeriodic_d の実装と G1-a/G1-b・R3 の裏付け (5) 旧延長 run の step0 rms_roK スパイクを F1 上書きの仮説とする記述 (6) accepted にしてよいか
- **extra**: `case/09.Taylor-Green/_g0_lsq_seam/README.md`, `case/09.Taylor-Green/_g0_lsq_seam/R3.txt`, `case/09.Taylor-Green/_g0_lsq_seam/G_tgv.txt`, `case/09.Taylor-Green/r2_conservation.txt`, `case/09.Taylor-Green/r2_ke_budget.txt`, `case/39.periodic_hills/run_0039_r1_gradfix_new_ext/QUASISTEADY_1p6M.txt`, `case/39.periodic_hills/run_0038_r1_gradfix_old_ext/CONVERGENCE_VERDICT_fromfloor.txt`, `case/39.periodic_hills/r1_extract.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

合併 LSQ と SST 勾配の gather は支持します。ただし、適用範囲外の軸対称に影響する実装漏れと、訂正後も未達の検証があります。  
R2 の保存性は裏付けられましたが、「注入が消え、数値散逸のみになった」という断定は現在の収支評価では証明できません。

1. **Major — 化学種・受動種の半割面除外が、設計した適用範囲を越えています。**

   **根拠:** [speciesTransport_d.cu:690](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/speciesTransport_d.cu:690) と同ファイル1232行は、`excludePeriodic` に `periodicNodeActive` を渡しています。plan §4.2a が要求する `periodicSeamMergeActive` ではありません。一方、[periodicNode_d.cu:174](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/periodicNode_d.cu:174) は軸対称の勾配 gather を実行しません。

   したがって、**node・軸対称・周期・スカラー再構成有効**では、今回から周期半割面だけが除外され、反対側の寄与は補充されません。定数場でも、除去した面の寄与に対応する `−φ S_periodic/A_planar` が従来勾配との差として生じます。「軸対称の既存経路は変えない」に反します。

   R3 の軸対称ケースは [r3_prepare.py:107](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/r3_prepare.py:107) のとおり入口・出口・壁・軸で、周期境界がありません。既存3ケースの `res_1.h5` を再比較すると **R3: PASS** でしたが、この条件漏れは検出できません。

   **対案:** 両呼出しを `periodicSeamMergeActive` に統一し、軸対称×周期の定数・線形スカラーで旧経路不変を確認してください。回転周期も今回の除外変更から外します。

2. **Major — G2′ のゼロ成分は、訂正後の基準でも一部不合格です。**

   **根拠:** [plan:134](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:134) は、ジッタ格子のゼロ成分にも次数判定を適用しています。保存済み `g_harness_g2p_k0.2_N{16,32,64}_j0.2/res_1.h5` を、ハーネスと同じ抽出条件で再集計しました。

   | 水準 | ゼロ成分の継ぎ目最大誤差 | 収束次数 |
   |---|---:|---:|
   | 16³ | 0.0151796062 | — |
   | 32³ | 0.0084277764 | **0.84891** |
   | 64³ | 0.0042566396 | 0.98544 |

   16³→32³ は要求の **0.9未満**です。誤差は十分に丸め床を超えています。[g2p_jitter.py:99](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/g2p_jitter.py:99) の次数計算は主場だけを対象にし、ゼロ成分には現在も旧基準の「床以下」を適用しています。

   **対案:** ハーネスを訂正後の基準と一致させ、この結果を未達として §5.1 に戻してください。漸近領域への到達を追加水準で確認する検証計画を先に定め、現行結果を保存したうえで再判定すべきです。これは直ちに LSQ 実装の誤りを意味しませんが、「検証完了」は成立しません。

3. **Major — G1 の PASS は、現行仕様の全条件を検査した PASS ではありません。**

   **根拠:** [g_suite.py:265](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/g_suite.py:265) は、G1-b の閾値を `max(1e-5*S, 2*Nmax*εφ/h)` とし、超過時にはさらに `n_member` 倍へ広げます。plan の現行基準には、この自動緩和はありません。また、継ぎ目／内部の誤差比は表示だけで、PASS 判定に入りません。

   保存 HDF5 から7ケースの誤差比を再計算すると、いずれも2以内でした。したがって半割面修正の効果自体は支持できます。しかし、[plan:163](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:163) が明記する **GPU 定数スカラー場の未実施**は残っています。ハーネスの定数場計算は CPU double・半割面込みの診断であり、修正後 GPU の定数保存試験ではありません。

   **対案:** 自動緩和を除き、誤差比と定数場を機械判定に組み込んでください。§5.1 #5a/#5c にある直接の化学種勾配確認も、実施または代理試験で代替する根拠を明示して閉じる必要があります。

4. **Major — R2 の符号判定は再現できますが、「収支の残りは数値散逸のみ」という帰属は過剰です。**

   **根拠:** [r2_ke_budget.py:50](/home/sano/work/forge-sern-design/case/09.Taylor-Green/r2_ke_budget.py:50) は空間微分を独自の周期中心差分で評価し、76行では snapshot 間差分から `dK/dt` を作っています。これはソルバが実際に用いた面応力・圧力仕事・対流流束の離散収支ではありません。

   CSV から次は確認できました。

   - 新の `max r = −0.0074228083`
   - `t≤7` の `min(r_old−r_new) = +0.0175604664`
   - 旧の `K/K0>1` かつ `r>0` は15 snapshot連続
   - 旧の `max K/K0 = 1.0187202`

   しかし `r` には、数値散逸に加えて**空間作用素の不一致、snapshot 時間微分誤差**が含まれます。負の総和から、正の寄与が消失したとは断定できません。旧新差を今回の修正に帰属することと、その差を「継ぎ目でのエネルギー注入」と機構まで確定することは別です。

   さらに [r2_ke_budget.py:105](/home/sano/work/forge-sern-design/case/09.Taylor-Green/r2_ke_budget.py:105) と保存された結果は、今も **`VERDICT (新 |r| ≤ 0.05): FAIL`** です。訂正後の符号条件は判定ツールに反映されていません。

   **対案:** 今回は「修正後に KE 超過が消え、この後処理での収支残差が負になった」という観測までに結論を限定してください。符号条件を正式な保存性ゲートの代替にせず、測定後に追加した診断として実装・記録を同期します。機構の確定は、同一状態に対する旧新の離散残差から KE 仕事を分解する後続検証へ残すのが妥当です。

5. **Major — R1 の限定合格には根拠がありますが、「他の基準はすべて成立」は誤りです。**

   **根拠:** 再実行した `check_quasisteady.py` は、保存記録と同じ結果でした。

   - `case/39.periodic_hills/run_0038_r1_gradfix_old_ext/`：**ALL STEADY**
   - `case/39.periodic_hills/run_0039_r1_gradfix_new_ext/`：**dF1_inf = DRIFTING、NOT ALL STEADY**

   新の末尾64点の `dF1_inf` は **0.03334079～0.03344835**。全点が0.05未満なので、**差の大きさに対する限定合格は支持します**。ただし増加33回・減少30回で、単調ではありません。漸近値0.0335を再現する推定手順も提示されていません。

   また旧延長の [CONVERGENCE_VERDICT_fromfloor.txt:3](/home/sano/work/forge-sern-design/case/39.periodic_hills/run_0038_r1_gradfix_old_ext/CONVERGENCE_VERDICT_fromfloor.txt:3) は **NOT CONVERGED**。`rms_roK` のピーク408.52倍を、末尾0.94倍だけで合格扱いにはできません。[plan:215](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:215) の完了記述と矛盾します。

   F1 初回上書きは旧コードに実在し、スパイク原因の**仮説としては妥当**です。ただし、この restart での原因確認や PASS の代わりにはなりません。

   **対案:** 「単調収束中・漸近値」を「微小な増加傾向と変動が残る、観測窓内で上限成立」に訂正し、未達を §5.1 に残してください。旧側は元 run と延長の同一実効設定区間を正式に連結して判定し、スパイクを含む全履歴を保存します。F1 仮説は同一 restart 入力の短い切り分けで確認してください。新の定常性確認も継続課題として明示します。

6. **Minor — 現在仕様・スコープ・残作業表の同期が不足しています。**

   **根拠:** [methods/gradient.md:156](/home/sano/work/forge-sern-design/methods/gradient.md:156) は現在も「LSQ が2倍」「SST未合算」を現状として記述しています。[plan:47](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:47) は化学種・受動種 GG を対象外としていますが、実装では変更済みです。[plans/README.md:27](/home/sano/work/forge-sern-design/plans/README.md:27) も `draft` のままです。

   §5.1 #5d は古い R2 判定式・閾値を保持しつつ「完了」。G2′ の新基準での未達、旧延長の収束未達、F1 仮説の検証が独立した残作業になっていません。

   **対案:** 現在仕様を実装後の記述へ更新し、過去の欠陥・撤回した判断は履歴として分離してください。上記未解決事項を §5.1 に登録し、case の run 一覧も最終状態へ同期します。

補足すると、R1 抽出の第一内部層との比較、x継ぎ目除外、勾配ノルムの L2 比は、z継ぎ目の診断として妥当です。F1 再計算も今回のモデルの式・定数と対応しています。ただし実際に輸送が読んだ F1 の直接測定ではありません。x継ぎ目列も、今回の格子での参考指標という扱いを支持します。これらを局所勾配全成分の一致証明へ広げるべきではありません。

R2 の保存 CSV は全4本で既定の保存上限を満たしました。`check_convergence.py` の **NOT CONVERGED** は非定常 TGV では想定どおりです。一方、R1 最終 run の残差 CSV と本番 HDF5 はローカルに無いため、収束については保存 VERDICT の照合までで、元履歴からの再判定はできていません。

**推奨は、`active` に留め、1→2・3→5→4→6の順に修正・再判定してから再レビューすることです。** 中核実装を撤回する必要はありませんが、現状での `accepted` 移行は認められません。ファイルは変更しておらず、指摘は plan 未反映です。

指摘数: Critical 0 / Major 5 / Minor 1
