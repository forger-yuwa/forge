# codex レビュー: boundary-node-periodic-gradient-fix (result)

- **plan**: [`plans/active/boundary-node-periodic-gradient-fix.md`](../../plans/active/boundary-node-periodic-gradient-fix.md)
- **stage**: `result` (diff base `1266aba1`)
- **date**: 2026-09-26
- **commit**: `80a659b7` (feature/sern-design)
- **codex**: effort `high`, 4.6 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M4/m1
- **focus**: 2 回目 (1 回目 NO-GO C0/M5/m1 への対応確認): (1) M1 の条件統一と #6a を決定論的試験 (面流束・状態のビット一致) で閉じたことの妥当性 (2) #6c の比 (場ごと) と定数場の壁節点除外が事後の緩和でないか (3) #6b 判定不能を accepted の障害にしない扱い (4) R1 の限定合格 + dF1_inf 定常性を継続課題とする扱い、旧延長の連結区間 PASS と F1 切り分け (5) R2 の結論が観測と帰属に限定されているか (6) accepted にしてよいか
- **extra**: `case/09.Taylor-Green/_g0_lsq_seam/README.md`, `case/09.Taylor-Green/_g0_lsq_seam/R6a_axi_periodic.txt`, `case/09.Taylor-Green/_g0_lsq_seam/R6a_massflux_bitcheck.txt`, `case/09.Taylor-Green/_g0_lsq_seam/G_channel.txt`, `case/09.Taylor-Green/_g0_lsq_seam/G2p_jitter.txt`, `case/09.Taylor-Green/r2_ke_budget.txt`, `case/39.periodic_hills/run_0038_r1_gradfix_old_ext/CONVERGENCE_VERDICT_concat.txt`, `case/39.periodic_hills/run_0039_r1_gradfix_new_ext/CONVERGENCE_VERDICT_concat.txt`, `case/39.periodic_hills/_f1split/F1_SPLIT.txt`, `notes/reviews/2026-09-26-boundary-node-periodic-gradient-fix-result.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

合併 LSQ、SST 勾配の合算位置、`periodicSeamMergeActive` への条件統一は支持します。  
ただし、追加試験2件に検出力の問題があり、G1・G2′の未達を完了扱いにする根拠も不足しています。現状での `accepted` 移行は認められません。

1. **Major — #6a の質量流束ビット一致では、修正対象のスカラー作用素を検証できません。**

   **根拠:** [R6a_massflux_bitcheck.txt:7](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/R6a_massflux_bitcheck.txt:7) が比較したのは `massflux` 5104面と `ro,Ux,Uy,Uz,P,sonic` です。`Y`・`Xi`・その勾配・スカラー面流束は含まれません。

   この試験の `speciesFaceReconstruction: 1` では、組成勾配は面エンタルピーに使われますが、通常の質量流束式には入りません。[convectiveFlux_slau_d.inc.cuh:579](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:579) の `mdot` が一致しても、組成再構成の正しさは保証されません。受動種についても、[speciesTransport_d.cu:1275](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1275) は再構成面値を使う条件を `>=2` としています。

   条件統一自体は正しく、コードから軸対称で旧新とも半割面を除外しないことは確認できます。しかし、このダンプから「スカラー面作用素も同一」「場の差はすべて集積順序」と結論するのは過剰です。

   **対案:** 同一入力で `species_gradient_d` が生成する**集積前の面寄与**を比較する決定論的試験に置き換えてください。化学種・受動種の両呼出しを対象にし、誤った `periodicNodeActive` 条件へ戻すと失敗することも確認します。既存の `R3` 規則上の `FAIL` は履歴として保持します。

2. **Major — G2′の次数未達は、内部との誤差比だけでは閉じられません。**

   **根拠:** [G2p_jitter.txt:20](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/G2p_jitter.txt:20) の本判定は次のとおりです。

   | 成分 | 継ぎ目の次数 | 内部の次数 |
   |---|---:|---:|
   | `φ_y` | **0.855** | 0.976 |
   | `φ_z` | **0.869** | 0.939 |

   どちらも継ぎ目で要求の0.9未満です。128³で誤差が内部より小さいことは支持できますが、**誤差の大小と収束次数は別の条件**です。内部の別成分に0.9未満があることから、この2成分の未達を「継ぎ目固有でない」と確定することもできません。

   [plan:218](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:218) は「判定不能」を記録しながら完了扱いにしていますが、§5.1 #7には次数確認の残作業がありません。

   **対案:** G2′を未達として戻してください。追加試験は、非対称な局所 stencil の無次元形状を固定し、`h,h/2,h/4` に縮小して継ぎ目全成分の次数を測る方式を推奨します。現在のジッタ系列は保存し、格子ごとに評価点集合が変わる最大誤差の診断として併記します。直ちに256³の本計算を追加する必要はありません。

3. **Major — 定数場から壁節点を除く変更は、検証範囲の事後縮小です。**

   **根拠:** [G_channel.txt:98](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/G_channel.txt:98) では、`k/ω` は境界条件で定数場が壊れた stencil を除外しています。一方、`Xi` は**除外0点**で定数場を保ち、最大誤差 `5.31 ε|φ|/h`、正式な判定は **FAIL** です。

   したがって、`Xi` の壁節点除外を「k/ωと同じ扱い」とする [plan:215](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:215) の説明は成立しません。壁∩継ぎ目と非継ぎ目壁で最大値が同じことは既存の壁閉包誤差を示唆しますが、当初対象だった壁∩継ぎ目の検証を不要にはしません。ハーネスも変更後の区分を判定していません。

   **対案:** 壁なし継ぎ目・壁∩継ぎ目・非継ぎ目壁を分けて機械判定してください。壁なしには現行の `4ε` を維持し、壁∩継ぎ目には、格納面ベクトルによる閉包誤差を含む参照との比較を追加します。絶対定数保存の未達は、根拠付きの既知制約として別途登録します。

   なお、**G1-bを場ごとの最大誤差ノルムで評価すること自体は妥当**です。全成分に絶対上限も課しており、成分別の不安定な比を参考値として残す整理は支持します。

4. **Major — #6d の `sstSigmaBlend=0/1` 比較には、F1初回上書き仮説を識別する能力がありません。**

   **根拠:** 旧版 `1266aba1` の `ransTransport_d.cu:105–118` は、`sstSigmaBlend=1` の初回に `F1=1` を書きます。一方、`=0` は固定係数を使います。[scalarTransport_d.cu:213](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/scalarTransport_d.cu:213) の式では、両者とも初回の係数は

   `σ_k=0.85、σ_ω=0.5`

   になります。したがって、[F1_SPLIT.txt:9](/home/sano/work/forge-sern-design/case/39.periodic_hills/_f1split/F1_SPLIT.txt:9) の両条件でスパイクが408.7倍になることは、上書き仮説と矛盾しません。「F1仮説は不支持」とする判断は撤回すべきです。新バイナリとの比較は勾配修正も含むため、この原因だけを分離できません。

   **対案:** 旧版を基準に、**F1初回上書きだけを除いた版**との同一restart入力・初回残差比較を行ってください。それまでは「このA/Bでは識別不能」と記録し、#6dを完了扱いにしないでください。

5. **Minor — 残作業表とrun索引に、撤回済み・更新前の記述が残っています。**

   **根拠:** [plan:124](/home/sano/work/forge-sern-design/plans/active/boundary-node-periodic-gradient-fix.md:124) の #5d は、現在も「単調収束中・漸近値」「新 `|r|≤0.05`」「継ぎ目由来の注入」を現行手順として保持しています。いずれも結果欄では撤回済みです。

   [case/39 README:78](/home/sano/work/forge-sern-design/case/39.periodic_hills/README.md:78) は800kへの回し直し中のままで、延長行にも連結区間の `PASS` が未反映です。#7には壁閉包・G2′の継続課題もありません。

   **対案:** 現行判定を§5.1・§6・run索引で統一し、旧判断は履歴と明示してください。後続へ移す課題には、追跡先planと完了条件を付けます。

指定の `git diff 1266aba1...HEAD -- solver_density_cuda design methods procedures` を取得し、HEAD `80a659b7` まで確認しました。中核実装に新たな符号・単位・適用条件の誤りは見つけていません。`methods/gradient.md` の現在仕様化と、化学種の代理として同じGGカーネルを通るトレーサを用いる説明も支持します。

R1のCSVを `check_quasisteady.py` で再判定した結果は、旧 `run_0038_r1_gradfix_old_ext/` が **ALL STEADY**、新 `run_0039_r1_gradfix_new_ext/` が **dF1_inf: DRIFTING / NOT ALL STEADY** でした。新の末尾64点は `0.03334079–0.03344835`。**観測窓内の上限成立に限った限定合格**と、定常性を継続課題にする扱いは支持します。連結区間の収束は保存VERDICTで旧新とも **PASS** ですが、残差元CSVがローカルにないため再判定はできていません。

R2は `run_0173–0176` の保存CSVから全4本の保存上限成立を再確認しました。符号診断も `max r_new=-0.0074228083`、早期の旧新差最小 `+0.0175604664`、旧の正値15点連続を再現しました。現在の「観測とcommit範囲への帰属」に限定した結論は妥当です。

**推奨は、`active` に留め、上記4件の検証を閉じて文書を同期した後に再レビューすることです。** ファイルは変更しておらず、指摘はplan未反映です。

指摘数: Critical 0 / Major 4 / Minor 1
