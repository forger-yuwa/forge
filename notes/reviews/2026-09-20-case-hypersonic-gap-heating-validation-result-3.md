# codex レビュー: case-hypersonic-gap-heating-validation (result)

- **plan**: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-20
- **commit**: `44c0c5e4` (feature/sern-design)
- **codex**: effort `high`, 6.0 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M7/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

回帰係数の数値は再現できましたが、「幅依存の傾きが実測と12 %以内で一致する」という検証結論は支持できません。  
合否台帳の出典誤り、SU2 比較条件の不一致、未達ゲートが残っています。accepted への移行は時期尚早です。

指定の diff を取得し、ケース実装・保存場・残差・判定記録を確認しました。diff には別計画の CHT／SERN／limiter 変更も含まれるため、本レビューをブランチ全体の回帰合格とは扱えません。ファイルは変更していません。

1. **Major — 幅スケーリングの回帰を、検証精度と観測誤差の分離に読み替えている。**

   根拠: [width_scaling.py:74](/home/sano/work/forge/case/50.deep_cavity_wieting_m7/tools/width_scaling.py:74) は4点の無重み回帰です。[plan:239](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:239) はそこから「幅に依存しない加算項」と断定し、[plan:270](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:270) で12 %の一致を受け渡しています。

   保存 JSON から再計算すると、全点回帰の傾き **0.88324**、加算モデル RMS **8.521 W/m** は再現します。しかし、隣接幅の増分比 `ΔQmeas/ΔQforge` は **1.180／0.812／0.622**。1点を除いた回帰の傾きも **0.737–0.985** に動きます。差そのものは **5.85–29.59 W/m** です。これでは幅に依存しない誤差を識別できず、12 %を精度上限にもできません。

   **対案:** 結論を「4点への事後的な回帰結果」に限定する。測定・digitize・前壁理論補完・格子誤差を傾きに伝播させ、独立条件で検証するまでは case/49 への定量的な受け渡しを止める。

2. **Major — T4 の事前登録に、別実験の精度限界と異なる正規化が混入している。**

   根拠: [tolerances.json:9](/home/sano/work/forge/case/56.gap_tp1187/tolerances.json:9) は **680 W/m²を TP-1187 の測定精度**としています。しかし、この値の出典は手元の **W70・TN D-5908、本文p.7／PDF第9ページ**です。TP-1187 は本文p.8／PDF第10ページで **`q/qFP < 0.02` の精度に疑義**を記しています。

   さらに [acceptance.json:81](/home/sano/work/forge/case/56.gap_tp1187/acceptance.json:81) は比較量を `h/h_FP` と記載し、同じ680 W/m²と「6点中4点」を合否に使います。plan §4.4a が規定する TP-1187 の比較量は冷壁 **`q/qFP`** です。「4/6」の根拠も用途要求や当該実験の不確かさではありません。

   **対案:** TP-1187 の量・精度・検出限界から台帳を作り直す。精度を確定できない点は合否対象から外し、任意の点数基準で一次検証を成立させない。

3. **Major — 段階起動の失敗検出と結果保存が、修正後も完結していない。**

   根拠: 実際に比較計算を起動した [_xcheck_forge2.sh:40](/home/sano/work/forge/case/56.gap_tp1187/_xcheck_forge2.sh:40) は、依然として `run_case.sh ... || true` で失敗を無視します。44行目は出力が存在した場合だけ引き継ぎ、その後も `stage ... done` と表示します。前段と同名の出力も残ります。

   また、[case/50 make_case.py:193](/home/sano/work/forge/case/50.deep_cavity_wieting_m7/tools/make_case.py:193) と [case/55 make_case.py:217](/home/sano/work/forge/case/55.gap_turbulent_m5/tools/make_case.py:217) は**本段終了後も全 `res_*` を削除**し、標準名の残差・VERDICTを改名します。最終場だけを入力メッシュへ戻しても、時系列による検証は再現できません。

   **対案:** 各段を別ディレクトリに保存し、その段で新規生成された出力・終了コード・有限性を確認する。本段の保存場、標準名の残差、VERDICTは保持する。実際の起動スクリプトにも同じ検査を適用する。

4. **Major — 組み直した SU2 比較も入口乱流条件が一致せず、本段の時系列も不足している。**

   根拠: [forge の入口設定:1](/home/sano/work/forge/case/56.gap_tp1187/run_0008_xcheck_forge2/bcondConfig.yaml:1) は `k=156.0306`, `omega=35363.15`。修正済み Sutherland 則を使うと、入口で

   `μ=1.36754e−5 Pa·s`、`ρk/(ωμ)=9.15157`

   です。一方、[SU2 設定:25](/home/sano/work/forge/case/56.gap_tp1187/run_0007_xcheck_su2/case.cfg:25) は粘性比 **10.0**。輸送則の変更を `omega` に反映していません。

   `run_0008_xcheck_forge2` の本段保存場は **step 0と120000の2点**だけです。残る2500／3000は前段出力です。本段の保存済み判定も **`NOT CONVERGED (stalled/plateau)`** でした。この出力では熱流束・温度・区間熱量の静定を確認できません。

   **対案:** 同じ実効物性から両者の `k/omega` を生成する。今回の条件なら粘性比10に対応する `omega` は約 **32362.8 s⁻¹**。同一設定区間で十分な間隔の保存場を残し、全残差と比較量の定常性を確認してから比較する。

5. **Major — 完了条件、とくに G9 を満たしていない。**

   根拠: [plan:661](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:661) 自体が残余変化評価・IC独立性の未実施を認めています。今回のツール再判定も次のとおりです。run の所在は [case/50 README](/home/sano/work/forge/case/50.deep_cavity_wieting_m7/README.md) を参照できます。

   | run（`case/50.deep_cavity_wieting_m7/` 配下） | `check_convergence.py` | `check_quasisteady.py` の深部量 |
   |---|---|---|
   | `run_0016_T1_wd0063_fine_settle` | `NOT CONVERGED`、末尾上昇 | `DRIFTING` |
   | `run_0008_T1_wd0211_long` | `NOT CONVERGED (stalled/plateau)` | `OSCILLATING` |
   | `run_0009_T1_wd0383_long` | `NOT CONVERGED`、末尾上昇 | `DRIFTING` |
   | `run_0011_T1_wd0524_long` | `NOT CONVERGED (stalled/plateau)` | `DRIFTING` |

   4 run の最終保存場には `VALUE/*` の非有限値はありませんでしたが、これは G9 合格を意味しません。T4 の3D検証も [plan:614](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:614) で未完です。

   **対案:** `active` に残す。残差、報告量、IC独立性、量別格子誤差をそれぞれ閉じる。達成不能なゲートは理由付きの未達結果として整理し、合格主張を維持するために基準を緩めない。

6. **Major — 現行 limiter で生成する計算と、引用している既存結果の対応が確定していない。**

   根拠: [solverConfig.cpp:547](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:547) は既定を `limiterScaled=1`、561行目は `venkatK=0.05` に変更しています。現在の [case/50 make_case.py:53](/home/sano/work/forge/case/50.deep_cavity_wieting_m7/tools/make_case.py:53) は参照値も固定します。

   一方、上記4幅の保存 config はいずれも `space: {convMethod: 1, limiter: 2}` のみです。`case/56/run_0005_gap_t8_deep/forge_run.log` には参照値の **`auto`** 決定が実際に記録されています。生成器を直した事実だけでは、既存値が修正後の作用素の検証結果にはなりません。

   **対案:** 各 run に実行バイナリの識別情報と実効 limiter 設定を紐付ける。旧経路は明示して保存し、採用する現行設定で代表ケースを再評価する。本 plan に依存するソルバ変更と検証範囲を追記する。

7. **Major — 「全件採用」の記録と、現在の本文・残作業表が矛盾している。**

   根拠: [plan:684](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:684) は前回の撤回を全件採用と記録していますが、本文には以下が残っています。

   - [183行目](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:183): 「数値側の検証は全項目クリア」
   - [598行目](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:598): 「分母 `qFP` は forge で作れる」
   - [654行目](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:654): 空気物性でのA1感度計算を使い、G1を「Table IVを0.2 %」として処理

   最後の主張は、生産物性では Gate B 最大 **20.26 %でFAIL**とする同じレビュー記録と両立しません。また §4.8 が要求する**物理収支診断**まで、[660行目](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:660) では離散収支の別 plan 待ちに置き換わっています。

   **対案:** 現在有効な結論を一本化し、旧主張には明示的な撤回表示を付ける。物理収支診断は本 plan の残作業に戻す。レビュー記録の「採用」は、対応箇所と成果物を確認してから完了扱いにする。

8. **Minor — 幅スケーリングが引用する VERDICT は、別の熱量を判定している。**

   根拠: [width_scaling.py:57](/home/sano/work/forge/case/50.deep_cavity_wieting_m7/tools/width_scaling.py:57) は前壁置換後の熱量を使いますが、[cavity_eval.py:200](/home/sano/work/forge/case/50.deep_cavity_wieting_m7/tools/cavity_eval.py:200) の時系列は全面 CFD 熱量です。

   今回、置換後の量を全保存場から再計算し、`check_quasisteady.py` の `classify` に通すと、4幅とも **`STEADY`**、drift は順に **0.4／0.8／0.6／1.6 %/tail** でした。最広幅は単調減少で、推定漸近値 **154.071 W/m** です。今回の結果は支持できますが、保存済みの検証経路は不正確です。

   **対案:** 置換後熱量を専用CSV列として出力し、その列のVERDICTと漸近値を引用する。

9. **Minor — run 索引と残作業表に、古い状態・重複・表構造の破損が残る。**

   根拠: [case/56 README:366](/home/sano/work/forge/case/56.gap_tp1187/README.md:366) は状態列欠落、369行目は余分な状態列があります。本文は `run_0005` を実行中としますが、実測残差は **799999 step** まであります。plan §5.1 には番号18／19／20の重複があり、「未作成」の記載と作成済み `case/56/acceptance.json` も併存しています。

   **対案:** run ごとの実施済み区間・成果物・判定日を同期し、残作業IDを一意にする。

**推奨は、accepted 移行を止め、この plan を `active` に残して検証を継続することです。** 優先順は、①12 %精度の受け渡しを撤回、②T4合否台帳の出典と量を修正、③起動・保存経路とSU2条件を修正、④固定した実効設定で未達ゲートを評価、⑤本文・索引を同期、です。これはレビュー上の提案であり、plan 未反映です。

指摘数: Critical 0 / Major 7 / Minor 2
