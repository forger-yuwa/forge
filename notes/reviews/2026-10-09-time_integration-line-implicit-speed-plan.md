# codex レビュー: time_integration-line-implicit-speed (plan)

- **plan**: [`plans/active/time_integration-line-implicit-speed.md`](../../plans/active/time_integration-line-implicit-speed.md)
- **stage**: `plan`
- **date**: 2026-10-09
- **commit**: `dc3bf1cc` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 5.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m1

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
課題の同定と A・D の方向は妥当です。既存 v2 とも重複しません。  
ただし、B の線形解検証、C の既定動作、性能・回帰試験の合格条件を直してから進めるべきです。

確認範囲：コード、accepted plan、保存済み JSON を照合しました。`run_0249`〜`run_0251` の計測原本と対象 run の残差 CSV はこの作業コピーに無く、18.97／36.37／37.28 ms/step は独立検証できていません。また、調査中の HEAD `4d394a71` には A が既に実装されています。以下は貼付された計画へのレビューです。

1. **Major — B は並列化と逆行列化を分離すべき。現在の検証では誤ったライン解を検出できない。**

   **根拠:** [timeIntegration_d.cu:1855](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1855) は Schur 補行列 `M = D − Kprev·W` を順次構築し、部分ピボット付き LU で解きます。[同:1691](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1691) には、ピボット適用順序で過去に誤解を返した記録もあります。

   逆行列を保存して掛ける方法は厳密算術では同じ系を解きますが、「丸め順序が変わる」だけの説明では、悪条件の Schur 補行列や長いラインでの誤差を評価できません。さらに [同:1895](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1895) は分解失敗時に補正を据え置くため、失敗しても NaN が出ず、全場 RMS では見逃す可能性があります。

   **対案:** B はまず **double・部分ピボット付き LU を維持した並列化**に限定する。凍結した実際の `D/K/rhs` から、緩和前のライン解について成分別にスケーリングした後退誤差  
   `η = ‖b − Ax‖∞ / (‖A‖∞‖x‖∞ + ‖b‖∞)`  
   を測り、独立した倍精度解との比較、拘束行、複数回の行交換、分解失敗件数を合格条件にする。逆行列化は別試験に分けるべきです。

2. **Major — C は「精度の選択肢追加」ではなく、既存設定の暗黙の精度低下になる。**

   **根拠:** [計画:66](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:66) は既定 `implicitSolvePrecision=0` で Thomas を float に変更します。しかし現行の Thomas 因子・前進代入は、このキーによらず double です（[mesh.hpp:199](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/mesh.hpp:199)、[timeIntegration_d.cu:1883](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1883)）。

   point の局所 5×5 解が float で成立することは、121 節点にわたる Schur 補行列の再帰も float で成立する根拠にはなりません。既存 YAML をそのまま使った再開計算まで動作が変わります。

   **対案:** 今回は **Thomas の既定 double を維持**する。float は明示的な opt-in とし、組立精度・係数保存精度・因子精度・前進代入精度の組合せを表で定義する。C の採否は、線形残差と長時間の収束性能を測ってから別途決めるべきです。

3. **Major — §6 は短期の回帰確認しか定義しておらず、目的とする「同じ品質までの高速化」を判定できない。**

   **根拠:** [計画:86](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:86) は 1／20 step の場の差と、C の 5000 step 軌道を挙げるだけで、C の許容値も未定です。

   保存済み `V0_repeat.json` の directional・密度の同一バイナリ再実行差は、1 step で **1.36〜1.61×10⁻¹⁴**、20 step で **8.48×10⁻⁸〜1.25×10⁻⁷**です。20 step の幅を一般的な許容差にはできません。

   また、`case/45.isobutane_m6_d155/run_0224_ns_coldmesh_tw300_linedir_tj5_cap50_ext/` は、[関連 plan:226](/home/sano/work/forge-integ-1005/plans/active/time_integration-implicit-thermal-jacobian.md:226) に記録された `check_convergence` 判定が **NOT CONVERGED**。保存済み系列 JSON を今回 `check_quasisteady.py --series-csv /dev/stdin` に渡した結果も、末尾 40000〜60000 step、`--drift 0.0005` で次のとおりです。

   ```text
   theta_r_40 / theta_r_70 / theta_r_94 / Q_w: DRIFTING
   OVERALL: NOT ALL STEADY
   ```

   これは保存済み派生量系列の再判定であり、HDF5 からの再抽出ではありません。それでも、通算 75000 step の結果が未達の水準を、5000 step で保証できないことは明確です。

   **対案:** 短期の数値回帰と、長期の到達品質を別ゲートにする。全保存量・`P/T`・壁／近軸の局所最大差を段数・設定別に登録し、長期評価には全残差の `check_convergence` と対象量の `check_quasisteady` を必須化する。短期試験を「収束解の一致」の証拠にしないこと。

4. **Major — 全節点がラインに載る case/45 と、条件未指定の case/39 一本では、共有経路を検証できない。**

   **根拠:** [main.cpp:2070](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2070) では `lineKFreeze=1` の後続 subiteration は **sweep 0 でも `storeLU=0`** です。[同:2094](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2094) は Thomas・swap 後の周期補正ミラーを要求します。一方、[mesh.cpp:1129](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/mesh.cpp:1129) は短い鎖を point に戻すため、ライン／point 混在が正規の動作です。

   **対案:** 最小限、`lineImplicit=0`、被覆ゼロ／部分被覆、可変ライン長、`lineKFreeze=0/1`、周期ミラー、軸・壁の拘束行を検証表に追加する。A は `loop>0` ではなく **`onLine && !storeLU`** を契約とし、RHS の拘束処理は毎 sweep 残す。D は面積係数、TP の固有ベクトル、`rowDec`、既存粘性加算を照合する。

   実行試験は [verification/README.md:63](/home/sano/work/forge-integ-1005/procedures/verification/README.md:63) に従って node のみでよく、cell は未検証と明記すればよいです。

5. **Major — 現在の計測手順だけでは、生産時の壁時計短縮を合格判定できない。**

   **根拠:** [prof_runs.sh:7](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/prof_runs.sh:7) は `FORGE_PROFILE=1` を設定し、[main.cpp:929](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:929) は計測区間ごとに `cudaEventSynchronize` します。line 経路には Thomas 用の追加計測区間もあります。したがって、計測下の壁時計には通常実行と異なる同期費用が入ります。

   また、`ic*25+i*5+j` の係数配置に対してスレッドを増やすだけでは、帯域下限への接近は保証されません。実転送量は warp 内のアクセス配置に依存します（[NVIDIA CUDA Best Practices Guide](https://docs.nvidia.com/cuda/archive/13.0.1/cuda-c-best-practices-guide/index.html)）。

   **対案:** カーネル内訳は `nsys`、採否は **native・専有 GPU・`FORGE_PROFILE=0`・profiler なし**の反復計測で決める。初期化／出力を除いた同一区間の中央値・ばらつきと、必要な品質までの総時間を登録する。B の lane 配置、係数ロード、レジスタ／spill を具体化し、22〜24 ms は未検証の目標として扱うべきです。

6. **Minor — C の実装範囲と、並行する粘性 Jacobian 計画との接続が不足している。**

   **根拠:** §7 は CUDA の二ファイルだけを挙げますが、因子の型は [mesh.hpp:201](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/mesh.hpp:201)、確保サイズは [mesh.cpp:1155](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/mesh.cpp:1155) にあります。カーネルのテンプレート化だけでは、計画した保存容量・帯域削減になりません。また、`time_integration-line-viscous-jacobian.md` は同じ `D/K` の組立を変更します。

   **対案:** 型・確保・wrapper・設定・退避経路を影響範囲に追加し、粘性 Jacobian 計画を相互参照する。高速化の A/B では `D/K` の仕様を固定し、双方を統合した後の再検証を残作業表に置く。

**推奨は、A → D → B（double・部分ピボット付き LU 維持）の段階導入です。C と逆行列化は別の opt-in 実験に分離してください。** 既存 v2 が解決したのは毎 sweep の再分解費用であり、今回の不要対角組立・Thomas 実行効率の改善には独立した価値があります。

実装前の修正順は、①既定 double の維持、②線形解・回帰・長期品質ゲートの定量化、③通常実行の性能判定追加、④影響範囲と plan 間の統合順の明記です。本提案は **plan 未反映**です。依頼どおりファイルは変更していません。

指摘数: Critical 0 / Major 5 / Minor 1
