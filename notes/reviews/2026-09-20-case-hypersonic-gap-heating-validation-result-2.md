# codex レビュー: case-hypersonic-gap-heating-validation (result)

- **plan**: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-20
- **commit**: `fad2d2d2` (feature/sern-design)
- **codex**: effort `high`, 9.5 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M9/m1
- **focus**: **2 巡目**。1 巡目 (同日, C0/M10/m1) の指摘を全件採用して撤回・修正した後の状態を見てほしい。焦点は 3 つ。(A) **撤回と修正が十分か** — 過大主張がまだ残っていないか、撤回したつもりで別の場所に残っていないか、コード修正 (段間引き継ぎ M9・リミッタ基準値 M10) が正しいか。(B) **残っている肯定的主張が持つか** — Gate A/B の条件再構成 (Table IV と 0.2 %、TP-1187 の理論偏り 1.10/1.30 を 1.072/1.335 で再現)、forge が Burggraf 理論と 4 幅で 1 %、参照温度法の平板と 1-2 %、Qc の幅依存。これらも同じ厳しさで検算してほしい。(C) **新しく追加した 2 つの結論** — (i) 観測モデルを一貫適用したら薄板伝導では説明しきれない (#30)、積分量も区間を切ると壊れる (#36)、という否定的結論の導出が正しいか。(ii) 深部の静穏さが入口条件依存だった (#39) という結論。あわせて SU2 クロスチェック (#38) の設計 — 同一 gmsh モデル・理想気体+Sutherland・一様入口・SST で揃え、判定を残差でなく深部 |U| の静定で行う — に穴がないかも見てほしい。SU2 は現在 200k 反復で継続中で未完了
- **extra**: `case/50.deep_cavity_wieting_m7/README.md`, `case/51.gap_turbulent_thd8233/README.md`, `case/55.gap_turbulent_m5/README.md`, `case/56.gap_tp1187/README.md`, `notes/reviews/2026-09-20-case-hypersonic-gap-heating-validation-result.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

撤回は不完全です。新たに、入口条件比較の交絡、SU2との物性・SST設定不一致、段間引き継ぎの失敗検出漏れ、理論積分・薄板モデルの実装誤りを確認しました。
数値の一部は再現できますが、検証完了として `accepted` へ移す根拠は揃っていません。

指定ファイルを読み、`git diff main...HEAD -- solver_density_cuda design methods procedures`、未追跡の対象case、保存場・残差・原論文を確認しました。ファイルは変更していません。以下の対案は**plan未反映**です。

1. **Major — 撤回した主張が、現在の設計判断と評価プログラムに残っている。**

   [plan:94](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:94)は、依然として「独立検算」「分母はCFDで作れる」としています。[同:546](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:546)では、撤回対象だった「照合先はclean」「未計測断面をcase/49へ渡せる」が3D実施の判断根拠です。[同:345](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:345)にも「3つの独立経路」「第一原理から再構成できない」が残っています。

   文書だけではありません。[gap_eval.py:89](/home/sano/work/forge/case/56.gap_tp1187/tools/gap_eval.py:89)は現在も、実測より上ならforgeが誤りという、撤回済みの大小関係を出力します。

   **対案:** 訂正注記を追加するだけでなく、現役の結論・判断根拠・プログラム出力を置き換える。履歴を残す部分は撤回済みと明示する。T2降格とcase/49への誤差「未定量」は維持する。

2. **Major — #39の「入口を一様に変えただけで60倍」は、比較条件と保存時系列が支えていない。**

   比較するメッシュの座標・接続は一致しました。しかし設定は次のように違います。

   | 項目 | `run_0005_gap_t8_deep/` | `run_0006_xcheck_forge/` |
   | --- | --- | --- |
   | 熱力学 | `thermalMethod: 2`、5成分TP | `thermalMethod: 0`、定比熱 |
   | 輸送 | `viscMethod: 2`、混合気体輸送 | `viscMethod: 1`、Sutherland・固定Pr |
   | リミッタ基準 | 開始場から自動決定 | 自由流値に固定 |
   | 入口 | BLプロファイル | 一様 |

   根拠は両runの [生産設定](/home/sano/work/forge/case/56.gap_tp1187/run_0005_gap_t8_deep/solverConfig.yaml:4)と[比較設定](/home/sano/work/forge/case/56.gap_tp1187/run_0006_xcheck_forge/solverConfig.yaml:4)です。

   後者の120k保存場から、深さ5–10幅の温度幅 **0.033325 K**、最大速度 **0.111454 m/s**は再現できました。ただし本段の保存場は初期場と120k終端だけで、深部量の静定を判定できません。残差の再判定も **`NOT CONVERGED (stalled/plateau)`** です。

   **対案:** 現状は「複数条件を変更した計算で深部量が異なった」に限定する。同じ物性・実効離散化・初期条件で入口だけを変更し、深部速度・温度・熱流束の時系列を取得してから因果を判定する。

3. **Major — #38のSU2比較は、同じSutherland則・同じSSTになっていない。速度だけの静定判定も不足する。**

   [forgeの物性実装:58](/home/sano/work/forge/solver_density_cuda/cuda_forge/gasProperties_d.cu:58)では、Sutherland定数が `μ0=1.716e-5`、`T0=273.0`、`S=111.0` に固定されています。configの `visc` を変えても、この経路の係数は変わりません。一方、[SU2設定:13](/home/sano/work/forge/case/56.gap_tp1187/run_0007_xcheck_su2/case.cfg:13)は `1.570579e-5 / 273.15 / 110.4` です。

   `case/56.gap_tp1187/run_0006_xcheck_forge/res_120000.h5` の自由流粘性は **1.367541×10⁻⁵ Pa·s**。SU2設定が意図する **1.251515×10⁻⁵ Pa·s**より約**9.3%高く**、熱伝導率もずれます。

   またforgeは `dilatationCorrection: 2, katoLaunder: 1`、[SU2実行ログ:34](/home/sano/work/forge/case/56.gap_tp1187/run_0007_xcheck_su2/su2.log:34)は **SST-2003m、production modificationなし**です。これでは差を数値散逸に帰属できません。

   **対案:** 実際に評価された `μ(T)`・`λ(T)` とSST生産項を揃える。全保存量の残差確認を残し、深部 `|U|` に加えて `T`・壁熱流束・区間熱量の準定常判定を行う。速度が静定しても熱場の緩和完了は保証されません。

4. **Major — M9修正は正常終了時には有効だが、case/56では失敗段から古い場を引き継げる。**

   [make_case.py:177](/home/sano/work/forge/case/56.gap_tp1187/tools/make_case.py:177)は `check=False` の戻り値を検査せず、185行以降で `res_{nsteps}.h5` の存在だけを確認します。前段の出力も残しています。

   したがって、同じ2500 stepの `ramp1` が途中で失敗すると、`ramp0`の `res_2500.h5` が存在し、正常な当該段の最終場として再利用されます。出力が無い場合もエラーにせず次段へ進みます。[クロスチェック起動スクリプト:43](/home/sano/work/forge/case/56.gap_tp1187/_xcheck_forge.sh:43)も `|| true` で失敗を握り潰しています。

   **対案:** 段ごとに出力を隔離し、終了状態・当該段の到達step・NaN/Infを確認してから引き継ぐ。正常終了、同番号の旧出力を残した失敗、最終出力欠落の3経路を検証する。

5. **Major — M10の修正は部分的で、既存結果の「同一設定で延長」はまだ保証されていない。**

   case/56生成器の基準値固定自体は妥当です。しかし [case/50生成器:53](/home/sano/work/forge/case/50.deep_cavity_wieting_m7/tools/make_case.py:53)と[case/55生成器:53](/home/sano/work/forge/case/55.gap_turbulent_m5/tools/make_case.py:53)には固定値がありません。case/56も `limiterScaled`・`venkatK` を明示していません。

   既存生産runの[ログ:145](/home/sano/work/forge/case/56.gap_tp1187/run_0005_gap_t8_deep/forge_run.log:145)は今も基準値を **auto** と記録しています。さらに[ソルバ設定読込:547](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:547)では `limiterScaled` の既定が変更されています。生成器を直したことは、既存run間の作用素差やその影響を解消しません。

   **対案:** 対象runの実効設定とバイナリ版を台帳化し、新規比較では方式・基準値を固定する。旧auto系列と固定系列の差を測り、影響未評価の延長を同一設定の収束根拠として扱わない。

6. **Major — 「Burggraf理論と4幅すべて1%」は誤りで、理論積分にも端点欠落がある。**

   [plot_summary.py:18](/home/sano/work/forge/case/50.deep_cavity_wieting_m7/tools/plot_summary.py:18)は積分を `s=1e-6 m` から始めています。しかし実装式はリップで `s⁻¹ᐟ²` の可積分特異性を持ち、切り捨てた区間の寄与は無視できません。同じ式を変数変換して `s=0` から積分すると、次になります。

   | `w/d` | forge平均比 | 現行の理論積分 | 端点を含む理論積分 | forgeの差 |
   | --- | ---: | ---: | ---: | ---: |
   | 0.063 | 0.580428 | 0.575939 | 0.592013 | −1.96% |
   | 0.211 | 0.574892 | 0.575196 | 0.583970 | −1.55% |
   | 0.383 | 0.565330 | 0.567451 | 0.573960 | −1.50% |
   | 0.524 | 0.540272 | 0.560671 | 0.566220 | −4.58% |

   根拠runは `case/50.deep_cavity_wieting_m7/` 配下の `run_0006_T1_wd0063_long/`、`run_0008_T1_wd0211_long/`、`run_0009_T1_wd0383_long/`、`run_0011_T1_wd0524_long/`。残差は全て **`NOT CONVERGED`**、`Qc`系列は全て **`STEADY`**でした。最後の系列は単調減少で、判定器の漸近値は **176.188 W/m**です。

   また[README:371](/home/sano/work/forge/case/50.deep_cavity_wieting_m7/README.md:371)の「`Qc`が幅とともに減少」も量名が違います。総入熱は **21.60→75.03→133.57→176.53 W/m**と増加し、減少するのは `Qc/(w qfp)` です。

   **対案:** 特異端点を解析処理または変数変換して積分し、図表を再生成する。「1%検証済み」を撤回し、適用条件・残差未収束を添えた理論との比較値として報告する。

7. **Major — case/50の「Gate A/BがTable IVと0.2%」は、異なる物性モデルの単一系列を一般化している。**

   [conditions.py:55](/home/sano/work/forge/case/50.deep_cavity_wieting_m7/tools/conditions.py:55)を再実行した結果は、**Gate A: PASS、Gate B: FAIL、最大差20.26%**でした。

   空気物性・`Pr=0.75`・`T*=(Taw+Tw)/2`へ変更した感度計算では、A1の **29.2289対29.3 kW/m²**は再現できます。しかし全12系列の比は **0.96922–0.99757**です。「0.2%」はA1についての値で、燃焼生成物を使う生産モデルのGate Bを閉じた結果ではありません。[plan:601](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:601)はこの区別を失っています。

   一方、TP-1187の再計算では回復率中央値 **0.287481**、Re比中央値 **1.060690**、Gate B中央値 **1.072214／1.334793**を再現しました。`case/56.gap_tp1187/run_0002_fp_t8_long/`の平板も **83.486／77.339 kW/m²**、両量 **`STEADY`**です。ただし残差は **`NOT CONVERGED`**であり、これらは採用した近似式との整合診断です。

   **対案:** 生産物性による判定と代替物性による感度を別表にする。W70の原報の計算規約を確認できない場合は「分母差の説明候補」に留め、TP-1187も実験検証の合格へ読み替えない。

8. **Major — 新しい薄板モデルには、床座標と端点熱流束の実装誤りがある。**

   [cavity_eval.py:80](/home/sano/work/forge/case/50.deep_cavity_wieting_m7/tools/cavity_eval.py:80)が出す床座標は **`y/d`**です。しかし[skin_model.py:53](/home/sano/work/forge/case/50.deep_cavity_wieting_m7/tools/skin_model.py:53)は床への距離を **幅 `W`**で割って補間しています。例えば最狭幅の床データ範囲は約0–0.0625なのに、補間側は0–1へ進み、床の大半が端値の外挿になります。

   また[同:75](/home/sano/work/forge/case/50.deep_cavity_wieting_m7/tools/skin_model.py:75)は、戻り値の熱流束を作る際に断熱端でも `lap[0]=lap[1]` とします。最狭幅の全経路積分を再計算すると、

   - 入力熱量：**17.1733846 W/m**
   - 現行の見かけ熱量：**17.4647000 W/m**
   - 端点に同じ断熱差分を適用：**17.1733846 W/m**

   となります。報告した「1.017倍で保存」は、実際には端点処理の誤差です。

   **対案:** 床座標を `d` 基準に揃え、時間更新と観測量算出で同じ境界演算子を使う。全経路保存を検証してから #30/#36 を再評価する。区間積分の大きな感度まで、この修正で消えるとは主張しない。

9. **Major — #30/#36の数値は再現できるが、「伝導では届かない」「W70全体で定量検証不能」までは導けない。**

   `skin_model.py`の表と、`integral_compare.py`のモデル幅／CFD–実測差 **1.28／1.15**は再現しました。したがって「今回の二つの仮定で比較結果が大きく変わる」は支持できます。

   ただし[skin_model.py:88](/home/sano/work/forge/case/50.deep_cavity_wieting_m7/tools/skin_model.py:88)が調べるのは**2幅だけ**です。一定熱源を0.8秒与え、リップ温度を平板熱流束から直線的に上げる設定は、挿入中の熱履歴を確定したものではありません。断熱端とこのDirichlet端が、未知の接続・熱履歴を包含する上下界である証明もありません。

   実験の精度限界も判定に入っていません。最狭幅の `z/W=4`では、実測 **0.0600**に対する原報の精度限界だけで **±0.0232**あります。モデルの **0.0451**はその範囲内なので、[README:318](/home/sano/work/forge/case/50.deep_cavity_wieting_m7/README.md:318)の「伝導では届かない」は成立しません。

   **対案:** 結論を「評価した2幅・熱履歴・端条件では、観測モデルの感度が大きく定量誤差を確定できない」に限定する。原典追跡を中止する判断と、伝導機構の否定・全幅の比較不能という科学的結論を分ける。case/49へのモデル誤差を現時点で渡さない判断は支持します。

10. **Minor — run索引と残作業表が最新状態に同期していない。**

    [case/56のrun一覧:358](/home/sano/work/forge/case/56.gap_tp1187/README.md:358)に `run_0006_xcheck_forge/` と `run_0007_xcheck_su2/` がありません。残作業番号も[plan:540](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:540)以降で18・19・20が重複しています。

    今回確認時点では `run_0005_gap_t8_deep/` は800k保存場まで存在しました。50k–800k系列の再判定は、浅い3点が **`STEADY`**、8.44／14.11幅が **`DRIFTING`**、残差が **`NOT CONVERGED`**です。なお、`plans/README.md`のT2降格は既に修正されており、その旧指摘は解消しています。

    **対案:** 実行状態・評価区間・成果物・VERDICTをrun一覧へ同期し、残作業番号と完了／未完了状態を一意にする。

**推奨は、計画を `active` に留め、検証途中として継続することです。** 優先順は、撤回漏れの除去 → 引き継ぎ・物性・実効設定の修正 → 理論積分と観測モデルの再評価 → 条件を揃えた入口A/B・SU2比較です。T4本体と受け渡し方針が未完了の現状を、結果検証済みとして閉じることは勧めません。

指摘数: Critical 0 / Major 9 / Minor 1
