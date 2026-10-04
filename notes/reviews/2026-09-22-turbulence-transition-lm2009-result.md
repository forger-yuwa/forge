# codex レビュー: turbulence-transition-lm2009 (result)

- **plan**: [`plans/active/turbulence-transition-lm2009.md`](../../plans/active/turbulence-transition-lm2009.md)
- **stage**: `result` (diff base `ded72e17`)
- **date**: 2026-09-22
- **commit**: `20030477` (feature/sern-design)
- **codex**: effort `high`, 6.7 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m4
- **focus**: 実装の中間レビュー (T3A 平板まで完了、翼は進行中)。特に見てほしい点: (1) solver_density_cuda/cuda_forge/transition_d.cu のソース・陰的対角・上下限が SU2 (.external/su2-src/SU2_CFD/include/numerics/turbulent/transition/trans_sources.hpp, SU2_CFD/src/solvers/CTransLMSolver.cpp) と食い違っていないか。(2) ransSource_d.cu の SST 結合 (Pk_base, 破壊係数, F3) の順序。(3) main.cpp の呼び出し順 (transitionSource が ransSource の前、primitive が dependentVariables の後) に 1 反復ラグや未初期化の読みが無いか。(4) §6.2 で基準 (3) を『層流 run が Blasius ±5 %、遷移つき run の層流域は SU2 LM と ±5 %』に改めたのは逃げでないか。(5) 翼で入口粘性比 10/100 により正圧面が層流/乱流に割れる結果を、モデルの外挿としてどう報告すべきか。巨大ファイル (*.log, residual_history.csv, res_*.h5, plans/README.md) は読まないこと
- **extra**: `case/57.transition_flat_plate/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**  
SST 結合の順序と主要ソース式は妥当ですが、初期出力・restart・収束ゲートに修正すべき問題があります。  
T3A は準定常な比較結果として有用でも、計画の完了条件を満たしていません。`in_progress` を維持してください。

指定の diff を取得して確認しました。`Pk_base`、破壊係数、両方の `F1` 計算への `F3` 適用、残差組立て内の `transitionSource → ransSource` の順序に誤りは見つかりませんでした。周期の部分体積・残差／輸送対角 gather・保存量ミラーも実装されています。軸対称は設定検査で拒否されます。ただし、これらの実装確認と検証合格は別です。

1. **Major — 初期出力が未初期化の遷移変数を読み、restart も保存状態を忠実に引き継がない。**

   **根拠:** [variables.cpp:311](/home/sano/work/forge/solver_density_cuda/variables.cpp:311) の確保後、遷移配列はゼロ初期化対象に含まれません。[main.cpp:1225](/home/sano/work/forge/solver_density_cuda/main.cpp:1225) の初期化経路には `transitionPrimitive` がなく、そのまま `writeInitialOutputs` が呼ばれます。入力に遷移保存量がなければ、初期出力の `roGamma`／`roReth` も未初期化です。入力にあっても `gammaTr`／`reTheta`／`gammaEff` は未生成です。

   通常反復では保存量更新後に原始量を再生成せず出力します。一方、[interp_field.py:133](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:133) は保存済み `roGamma`／`roReth` を使わず、更新前の原始量と更新後の密度から再構成します。

   **対案:** 初期出力前に初期化・境界適用・周期同期を完了させる。restart は保存量を優先し、原始量からの復元は旧形式だけに限定する。初期出力と分割 restart の再現試験を追加する。

2. **Major — `--from-floor` が遷移残差の必須列検査を迂回する。**

   **根拠:** 必須列検査は [check_convergence.py:130](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:130) の `analyze` にしかなく、[同:281](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:281) の `analyze_from_floor` にはありません。

   ファイルを作らない模擬入力で、遷移有効・遷移残差２列欠落に対して **通常判定 `False`、床比較判定 `True`** を再現しました。例えば SST 参照から LM に切り替え、対象 CSV が遷移列を欠いていると、この穴を通ります。

   **対案:** 必須列・有限性の検査を両経路で共通化し、方程式系が違う参照を拒否する。§5.1 #2b は修正まで再開扱いにする。

3. **Major — 「単体検査完了」の証拠が、約束した検査を満たしていない。**

   **根拠:** [check_lm_kernel.py:22](/home/sano/work/forge/solver_density_cuda/tools/check_lm_kernel.py:22) の既定許容値は `1e-3`。判定は点ごとの相対誤差ではなく、**全域最大値で正規化した誤差の99.9百分位**です（同:49）。局所的な大誤差を落とせます。さらに `gammaEff` の誤差と反復上限到達は表示するだけで、PASS/FAIL に入りません（同:53）。陰的対角と更新後の上下限も検査対象外です。

   SU2 との数値上の相違もあります。[transition_d.cu:237](/home/sano/work/forge/solver_density_cuda/cuda_forge/transition_d.cu:237) の対角は SU2 のソース微分そのものではありません。例えば破壊項を無視し `γ=1` とすると、生産項由来の対角は SU2 微分の負部の **1.5倍**です。また SU2 の輸送変数 `Reθt` 下限は **`1e-4`** で、相関値の下限 **20** とは別です（[CTransLMSolver.cpp:102](/home/sano/work/forge/.external/su2-src/SU2_CFD/src/solvers/CTransLMSolver.cpp:102)）。

   **対案:** 同一状態を与える独立した点検査で、ソース・対角・`gammaEff`・ガード・更新制限を検査する。対角の追加減衰は設計として明示できるが、クリップについて「定常解では効かない」とするには、作動量と残差の確認が必要です。

4. **Major — T3A の「合格」と計画完了の判断が、保存済み VERDICT と計画自身の条件に反する。**

   **根拠:** 以下は `case/57.transition_flat_plate/` の保存済み判定です。

   | run | 収束判定 | 準定常判定 | `rms_roe` 最終値 |
   |---|---|---|---:|
   | `run_0005_t3a_lm_cont` | `NOT CONVERGED` | 全7量 `STEADY` | `5.08e-2` |
   | `run_0007_t3a_lm_coarse` | `NOT CONVERGED` | 全7量 `STEADY` | `6.95e-2` |
   | `run_0009_t3a_lm_fine` | `NOT CONVERGED` | 全7量 `STEADY` | `3.51e-2` |

   `cf_series.csv` の最終値は計画の表を裏付け、基準格子の準定常判定も再実行して `STEADY` でした。しかし [plan:153](/home/sano/work/forge/plans/active/turbulence-transition-lm2009.md:153) は未収束なら機構診断と明記しています。SST・層流でもエネルギー残差が同程度という事実だけでは、床の位置が前縁であることや、報告量への影響が無視できることまでは証明できません。

   SU2 T3A には判定ファイル・量の時系列が見当たらず、[README:51](/home/sano/work/forge/case/57.transition_flat_plate/README.md:51) の抜粋は運動量・エネルギー残差を含みません。さらに §5.1 の T3B、圧力勾配、翼、報告更新が未完了です。

   **対案:** T3A は「未収束・報告量は準定常・比較差は基準内」と記載する。残差床の局在と精度影響、SU2 側の判定を揃え、未完項目を閉じてから完了レビューに進む。

5. **Major — 「全格子で局所 `y₁⁺≤1`」を示した集計になっていない。**

   **根拠:** [cf_plate.py:91](/home/sano/work/forge/case/57.transition_flat_plate/tools/cf_plate.py:91) は最大値を **`x>0.02 m` に限定**しています。計画の最大値 0.63／0.45／0.32 は、前縁20 mmを含む全壁の保証には使えません。

   **対案:** 全壁の局所最大・超過面積・位置を報告する。幾何学的特異点を除くなら、その範囲と理由を明記し、除外範囲の格子感度を別に評価する。「全域」という表現は現状では撤回する。

6. **Minor — Blasius 基準の訂正は妥当だが、SU2 比較への置換を「実験検証」と呼ぶのは不適切。**

   **根拠:** `ref/t3a_exp.dat` から再計算すると、実験自身が Blasius に対し `x=0.195 m` で **+5.54%**、`0.295 m` で **+11.51%**です。したがって「遷移つき解にも Blasius ±5%」という旧基準の撤回は、逃げとは判断しません。

   ただし [plan:174](/home/sano/work/forge/plans/active/turbulence-transition-lm2009.md:174) の置換後の条件は、既存の SU2 実装照合と重複します。`x=0.3 m` の一点だけで層流域全体の ±5% も示せません。

   **対案:** 層流 run 対 Blasius、LM 対 SU2、LM 対実験を別の検査として記載する。実験との `Cf` 分布比較を残し、旧基準の不合格と訂正理由も保存する。

7. **Minor — 現在仕様と計画本文に、撤回済みの断定・実装との差が残る。**

   **根拠:** [theory.md:1116](/home/sano/work/forge/methods/turbulence/theory.md:1116) は翼の過大加熱の原因を「確認済み」と断定し、plan §1 の「最優先の仮説」と矛盾します。[solver-settings.md:269](/home/sano/work/forge/procedures/solver-settings.md:269) の「SU2 と同じ下限」も前述のとおり不正確です。

   また、実装は局所 `Tu` による初期化ですが §4.4 は全域入口値。変更ログの `Tu≤100%` は初期化・入口にだけ適用され、ソースの [transition_d.cu:166](/home/sano/work/forge/solver_density_cuda/cuda_forge/transition_d.cu:166) には上限がありません。予定された `recommended-settings.md` の遷移レシピもありません。

   **対案:** §4・§6・現在仕様を実装と検証範囲に合わせて同期する。未完の文書整備を §5.1 に追加する。`methods/index.md` の既存リンクは有効です。`plans/README.md` は指定どおり未確認です。

8. **Minor — 「50000 step 以降6桁不変」は実測と食い違う。**

   **根拠:** `case/57.transition_flat_plate/run_0005_t3a_lm_cont/cf_series.csv` では、`x_end` が50000 step の **0.84279609 m** から100000 step の **0.84964669 m** へ **0.813%**変化しています。`STEADY` 判定とは両立しますが、6桁不変ではありません。

   **対案:** 不変桁数の断定を削除し、評価窓・相対変動・VERDICT を記載する。

9. **Minor — 非退行検査の主張を独立に追跡できる記録が不足する。**

   **根拠:** [case README:53](/home/sano/work/forge/case/57.transition_flat_plate/README.md:53) の `run_0001_regr_none_head`／`run_0002_regr_none_new` は破棄済みで、実際に存在しません。ノイズ床倍率1.83という要約はありますが、入力・バイナリ識別・量別比較結果への恒久的な参照がありません。plan §6／§7 には「ビット不変」も残っています。

   **対案:** 非決定性を考慮した比較方針は維持し、コミット／バイナリ識別、実行条件、量別ノイズ床と差を小さな検証記録として残す。今回、既存機能の破壊を実証したわけではありませんが、非退行 PASS の追認もできません。

**推奨は、計画を `active/in_progress` に留め、指摘1～3の修正、指摘4～5の検証、文書同期の順で完了させることです。**

翼については、入口長さスケール未同定の感度として報告してください。`case/53.c3x_vane_cht/run_0146_lm_1um_cont/` の正圧面 `h` 偏差は準定常平均 **−33.25%**、単調減少の推定漸近値 **−33.30%**。`run_0147_lm_1um_mur100/` は今回の時系列再判定で `STEADY`、平均 **+25.48%**、推定漸近値 **+25.46%**でした。**両者とも収束判定は `NOT CONVERGED`** です。

したがって「遷移モデルで改善した」という一つの予測値にはまとめられません。粘性比10／100で約59ポイント動く条件依存性を示し、実験に近い粘性比を事後選択しないこと。「層流／乱流に分かれた」という機構の断定には、`h` だけでなく `gammaEff`・渦粘性・前縁付近の `Tu` の時系列が必要です。粘性比1、両翼、同一メッシュの ON/OFF、CHT は §5.1 の残作業として維持してください。

ファイル変更は行っていません。以上の対案はレビュー提案であり、plan 未反映です。

指摘数: Critical 0 / Major 5 / Minor 4
