# codex レビュー: discretization-node-face-weight-midpoint (result)

- **plan**: [`plans/active/discretization-node-face-weight-midpoint.md`](../../plans/active/discretization-node-face-weight-midpoint.md)
- **stage**: `result` (diff base `18916ed9`)
- **date**: 2026-09-22
- **commit**: `78930149` (feature/sern-design)
- **codex**: effort `high`, 2.0 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M2/m3
- **focus**: 2 回目の result レビュー。1 回目 (notes/reviews/2026-09-22-discretization-node-face-weight-midpoint-result.md、NO-GO C0/M4/m1) の指摘への対応を見てほしい: M1 → §5.1 #2c (平板を手順書どおり冷間起動し基準 run と 0.2 % 以内。§6 の合格条件を check_convergence の PASS から手順書の物理量基準 + SERIES STEADY に変更。理由は前縁特異点の既知の残差床で基準 run 自身が NOT CONVERGED のため)。この合格条件の変更が妥当か、それとも逃げか、を判断してほしい。M2 → #7 (製造解。熱伝導の内部面の式を Python に写したもので、運動量・粘性仕事・SST 拡散は未測定、ソルバ本体でない、という限界を明記)。この範囲の測定で accepted に移してよいか、運動量側の測定を前提にすべきか、推奨を 1 つに絞ってほしい。M3/M4 → #2b と §4.3 の表現を限定。読んでよいのは本 plan、同梱の調査ノート、solver_density_cuda/tools/mms_face_weight.py (sed -n '1,60p' と '95,150p')、procedures/verification/48-flat-plate-cooled.md、case/48.flat_plate_cooled_m4/run_0025_B_tw300_y3_fx05/CONVERGENCE_VERDICT.txt と run_0005_B_tw300_y3/CONVERGENCE_VERDICT.txt。禁止: *.log, residual_history.csv, res_*.h5, *.vtu, *.json (eval), plans/README.md, case/*/README.md の全文 (grep は可)。
- **extra**: `notes/investigations/2026-09-22-mms-node-face-weight.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

平板を既存手順の物理量基準で判定する変更は妥当で、逃げとは判断しません。  
ただし製造解の「壁面熱流束が二次精度」という主張は、コードの誤差定義では証明できていません。  
`accepted` への移行には、その訂正と、今回の主因である粘性応力・粘性仕事の測定を要求します。

閲覧制限を優先し、指定ファイル・指定行範囲だけを確認しました。実装 diff、境界・周期 seam・軸の実装、`methods/` と索引の整合性は今回再検証していません。ファイル変更も行っていません。

1. **Major — 製造解の熱流束誤差は、解析的な壁面熱流束に対する誤差ではない。**

   **根拠:** [mms_face_weight.py:141](/home/sano/work/forge/solver_density_cuda/tools/mms_face_weight.py:141) のコメントは「厳密な `k dT/dr` との相対差」ですが、143–146 行の実装は
   \[
   q_h=k_w\frac{T_h(r_w+h)-T_w}{h},\qquad
   q_{e0}=k_w\frac{T_{\rm exact}(r_w+h)-T_w}{h}
   \]
   の差です。厳密温度にも同じ二点差分を適用し、壁勾配近似そのものの誤差を除去しています。

   指定された製造解では、解析的な壁面流束を \(q_w\) とすると
   \[
   \frac{q_{e0}}{q_w}=\frac{1-e^{-h/\delta_T}}{h/\delta_T}
   =1-\frac{h}{2\delta_T}+O(h^2).
   \]
   \(\delta_T=50\,\mu{\rm m}\)、\(h=2,1,0.5,0.25\,\mu{\rm m}\) でこの差分自体の相対誤差を計算すると、**1.974%、0.993%、0.498%、0.250%**。次数は約 **1** です。これは run の実測ではなく、掲載された製造解からの直接計算です。

   したがって、[調査ノート:39](/home/sano/work/forge/notes/investigations/2026-09-22-mms-node-face-weight.md:39) の「壁面熱流束も `p=2.00`」は過大な主張です。温度誤差の測定まで否定するものではありません。

   **対案:** 現在の指標を「厳密温度の離散壁流束に対する差」と改名し、解析的な \(k_w\partial_nT_{\rm exact}\) に対する誤差を別途測定する。二次精度を要求する対象を内部面流束・温度解・壁流束に分け、結果に合わせて plan と調査ノートを訂正してください。

2. **Major — 熱伝導だけの測定では、今回の主要な変更経路の検証を完了できない。**

   **根拠:** [plan:47](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:47) では、問題のうねり rms は熱伝導 **109 W/m²** に対して粘性仕事 **571 W/m²**。主要因は後者です。一方、[調査ノート:46](/home/sano/work/forge/notes/investigations/2026-09-22-mms-node-face-weight.md:46) は粘性応力・粘性仕事を未測定と明記しています。

   限界の明記は適切ですが、未測定部分を検証済みに変えるものではありません。同じ `fx` を使っていても、テンソル応力と面速度の積について、スカラー熱伝導の次数をそのまま引き継げません。また、別スクリプトによる一つの場の仕事項再現は、この MMS の実装や格子収束の検証にはなりません。

   **対案:** **運動量側の測定を移行の前提にする**ことを推奨します。曲面・高 AR の格子系列で、複数速度成分を持つ製造場について応力流束と \(\tau\cdot U_f\) を測定してください。実カーネルを呼ぶ試験、または実カーネルとの直接照合を伴う試験にし、結果が出るまで §5.1 #7 を未完了に戻すべきです。全 SST 系の二次精度証明まで今回の前提に広げる必要はありません。

3. **Minor — 「`p≥1.8` を満たすのは `fx=0.5` のみ」は数値と矛盾する。**

   **根拠:** [plan:121](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:121) の断定に対して、[調査ノート:34](/home/sano/work/forge/notes/investigations/2026-09-22-mms-node-face-weight.md:34) と 37 行では旧式 `code` の流束指標の次数は **1.86 / 1.85**。いずれも 1.8 以上です。

   **対案:** 「この測定では `half` と `code` が閾値を満たし、`half` の誤差が小さい」と訂正する。Major 1 の指標修正後に合否も再判定してください。

4. **Minor — §4.3 の限定は改善したが、定常 A/B の合格証拠との対応がまだ曖昧。**

   **根拠:** [plan:96](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:96) は旧定義の試行 `run_0122` を `DRIFTING` と明記しています。この訂正は妥当です。しかし直後に補正後定義の別 run `run_0128` の `STEADY` を置き、「1.2 → 0.5」と結論づけています。これは [plan:131](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:131) の「同じ定義」「`STEADY` と確認した量」の合格証拠としては使えません。

   **対案:** §4.3 は「旧定義の有限時間区間における低減傾向」に限定する。合格の根拠は §5.1 #3 にある同一定義の生産 A/B に統一し、各 run の平均・変動幅・比・判定区間・VERDICT を併記してください。

5. **Minor — 合格条件の変更が残作業表に反映しきれていない。**

   **根拠:** [plan:114](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:114) の #2 は、依然として「`PASS` を取るのを #2c として残す」としています。一方、#2c は完了扱いで、§6 は物理量基準へ変更済みです。また、[plan:133](/home/sano/work/forge/plans/active/discretization-node-face-weight-midpoint.md:133) の「原理的に取れない」は、確認した二つの失敗判定からは導けません。

   **対案:** #2 を現行条件と実施結果に同期し、旧条件は変更履歴として残す。「原理的に取れない」は「今回の基準 run と試行 run は `PASS` を満たさない」に訂正し、Major 1・2 を優先順付きの残作業として登録してください。

平板の条件変更については、[手順書:15](/home/sano/work/forge/procedures/verification/48-flat-plate-cooled.md:15) に物理量・系列定常性・3% の回帰基準と既知の残差床が明記されているため、**物理量の回帰合格として採用することは妥当**です。ただし、残差収束とは区別してください。確認した判定は双方とも次のとおりです。

- [run_0025 の判定:2](/home/sano/work/forge/case/48.flat_plate_cooled_m4/run_0025_B_tw300_y3_fx05/CONVERGENCE_VERDICT.txt:2): `NOT CONVERGED (stalled/plateau)`
- [run_0005 の判定:2](/home/sano/work/forge/case/48.flat_plate_cooled_m4/run_0005_B_tw300_y3/CONVERGENCE_VERDICT.txt:2): `NOT CONVERGED (stalled/plateau)`

なお、新 run の `rms_roUy` は **1.38e−4**、基準は **1.60e−5** で **8.63 倍**です。「全残差が従来と同じ床」とも言えません。物理量差と `SERIES VERDICT: STEADY` は今回は plan の記載までの確認で、判定元データの独立検証はしていません。

**推奨は、恒久実装を維持したまま plan を `active` に留め、熱流束誤差の定義を訂正し、粘性応力・粘性仕事の格子系列測定を完了してから再レビューすることです。** 推奨は read-only のため **plan 未反映**です。

指摘数: Critical 0 / Major 2 / Minor 3
