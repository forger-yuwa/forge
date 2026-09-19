# codex レビュー: tooling-nozzle-sern-startup (plan)

- **plan**: [`plans/active/tooling-nozzle-sern-startup.md`](../../plans/active/tooling-nozzle-sern-startup.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `d1c42b8e` (feature/sern-design)
- **codex**: effort `high`, 5.7 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M6/m1
- **focus**: 起動レシピの妥当性。(1) 収束判定 5 条件、特に『全残差列がプラトー』を要求する形は妥当か。見落としはないか。(2) 本段 cfl 5 / 500 step は力係数の頭打ちで決めたが、遅い現象の発達を見逃していないか。(3) 設計 B だけ rms_roY1 が上昇する件の切り分け方。(4) MOC 初期値の作り方 (等エントロピー展開 + NASA-9、被覆外は端値クランプ外挿) に誤りはないか。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
起動時間を短縮する目的は妥当ですが、「全残差プラトー＝収束」と「1500×3＋500 step の生産レシピ確定」は支持できません。  
実測から500 step時点の未発達を確認しました。MOC 初期値にも、一定 γ の展開と NASA-9 熱力学の不整合があります。

1. **Major — プラトーは収束の十分条件ではなく、実装も「全列プラトー」を検査していません。**

   **根拠:** [sern_gates.py:217](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:217) は、標準収束判定が不合格でも、上昇・`still converging` がなければ受理します。実際、`case/46.sern_design/run_0179_ramp_B/` を再判定すると、

   ```text
   check_convergence: NOT CONVERGED (stalled/plateau)
   rms_roOmega: drop=0.8dec
   check_quasisteady（力係数）: ALL STEADY
   evaluate_gates: PASS
   ```

   でした。残差の大きさを問わず停滞を受理するため、離散方程式を十分満たしていない状態も通ります。float32 の更新消失やクランプで動かなくなった状態も区別できません。

   さらに [check_convergence.py:107](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:107) は、3桁以上低下した列には `still converging` を付けません。したがって **3桁以上低下してなお下降中の列は、現在のプラトーゲートを通ります**。仕様と実装が不一致です。

   **対案:** 標準収束 VERDICT と評価量の準定常 VERDICT を分離してください。プラトーは停滞の診断情報とし、受理には方程式別の無次元残差上限、保存収支、場・目的量の定常性を要求する。十分小さい残差が下降中であること自体を失敗理由にする必要はありません。基準を実証するまでは、現在の結果を「未収束・力係数は準定常」と扱うべきです。

2. **Major — 「1500×3＋500 step」は直接検証されておらず、500 stepで遅い変化を実際に見逃します。**

   **根拠:** ramp 成功例 `run_0178_ramp_A/`～`run_0180_ramp_C/` は、本段2000 stepです。[problem_ramp_B.yaml:111](/home/sano/work/forge/case/46.sern_design/problem_ramp_B.yaml:111) と同ファイル121行以降では、暖機・soft・mid も各2000 stepを指定しています。暖機は実装の整数除算により666×3です。これらを組み合わせても、計画のレシピを検証したことにはなりません。

   `case/46.sern_design/run_0179_ramp_B/` の先頭500 stepをメモリ上で切り出して再評価した結果は次のとおりです。ファイルは変更していません。

   | 対象 | 判定・実測 |
   |---|---|
   | `C_T_with_shear` | 正式分類関数で `STEADY` |
   | `rms_roK` | `falling / still converging`、低下1.0桁 |
   | `rms_roY1` | `falling / still converging`、低下2.6桁 |
   | 局所 `vis_turb` | `check_quasisteady.py`: **`DRIFTING`、末尾窓18.3%** |
   | ノード平均 `vis_turb` | 同ツール: **`TRANSIENT-UNSETTLED`** |

   局所点は500→2000 stepの `vis_turb` 絶対変化が最大のノード47159、座標約 `(1.526918, 0.340371, 0)` mです。`vis_turb` は500 stepの `1.58152e-3` から2000 stepの約 `1.63570e-3` へ変化しています。根拠は同runの `res_100.h5`～`res_2000.h5`。メッシュ品質は `MESH_QUALITY.txt` の **`VERDICT: PASS`**（最大AR502、skew0.313）です。

   **対案:** 本段500 stepの固定打ち切りを撤回し、まず提案レシピそのものをA/B/Cで検証する。長時間参照との比較には、力係数に加え、剥離位置・せん断層厚さ・組成断面・乱流粘性・保存収支を含めてください。合否許容差を先に決め、複数の連続窓で条件を満たすまで延長する方式を推奨します。`24000/5000=4.8` はstep数比であり、実時間の高速化率とも区別が必要です。

3. **Major — 異なる次元の残差を中央値で比較する条件③は、数値的な根拠がありません。**

   **根拠:** [sern_gates.py:168](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:168) は、生の `rms_ro`・運動量・エネルギー・`roOmega` 等を同じ集合に入れ、中央値の `1e6` 倍を閾値にしています。[residualMonitor_d.cu:124](/home/sano/work/forge/solver_density_cuda/cuda_forge/residualMonitor_d.cu:124) のRMS計算には、方程式間の物理スケールを揃える正規化がありません。

   この比較は単位・入口状態・化学種数に依存します。極端な異常を検知できても、「残差列の桁が揃う」ことは収束や方程式間の釣り合いを意味しません。

   **対案:** 各方程式を代表保存量・代表時間・CV体積に対応する尺度で無次元化し、それぞれの上限を設定する。局所異常は最大値・分位値・位置で検出してください。既存の中央値ゲートは補助警報に留め、収束条件から外すべきです。

4. **Major — MOC 初期値は NASA-9 に対する等エントロピー展開になっていません。**

   **根拠:** [runner_sern.py:359](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:359) と380行以降は、一定 γ の式で `T/P/q` を計算し、391行で内部エネルギーだけ NASA-9 に置き換えています。EOSに整合する `roe` は作れますが、入口と同じ全エンタルピー・エントロピーは維持しません。

   `problem_ramp_B.yaml` の `m6_on` 組成と、実装と同じ `FrozenGas` で再計算すると、指定 `M=3` で、

   ```text
   EOSからの実Mach = 2.96704
   入口に対する全エンタルピー差 = +1.20093%
   入口に対するエントロピー差 = +30.7478 J/(kg K)
   ```

   となります。これは初期値生成式の検算であり、CFDの収束解の誤差を示すものではありません。

   **対案:** MOC の `M,θ` を近似初期分布として利用しつつ、熱力学変換は同じ NASA-9 物性で、

   \[
   h_{\rm sens}(T)+\tfrac12 M^2\gamma(T)RT=h_{0,\rm in},
   \qquad
   p=p_{\rm in}\exp\!\left[\frac{s^\circ(T)-s^\circ(T_{\rm in})}{R}\right]
   \]

   を解き、`q=M a(T)` とする。[frozen.py:71](/home/sano/work/forge/design/forge_design/gas/frozen.py:71) に必要な関数はあります。

   被覆外の端値クランプ自体は初期値の近似として許容できます。ただし、それを被覆外のMOC解とは呼べません。クランプ領域の割合・位置、入口との整合、壁法線速度を記録し、別初期値から同じ最終状態に到達することを検証してください。

5. **Major — 設計Bを「形状固有・せん断層の乱流輸送」と絞り込む根拠が不足しています。**

   **根拠:** `rms_roY1` は領域全体のRMSであり、上昇だけでは発生場所を特定できません。`run_0182_sm1_B/` と `run_0185_sr_B/` に残る上昇は、いずれも**本段CFL5の履歴**です。soft/mid 内で発生したことを直接示していません。

   また、[runner_sern.py:744](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:744) 以降は段終了を主に終了コードと出力存在で判断し、restart後に段の `res_*` を削除します。段別診断の証拠が残りません。

   見落とせない候補は化学種更新です。該当configでは `speciesImplicitCoupling` 未指定＝0で、流れのblock-DPLURと異なる点陰的更新を使います。化学種の緩和は別の `speciesImplicitRelax`＝1です（[solverConfig.cpp:343](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:343)、[speciesTransport_d.cu:835](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:835)）。さらに [speciesTransport_d.cu:152](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:152) で非負化・再正規化します。これは `implicitRelax` 試験の解差についても調べるべき候補です。

   **対案:** 同一Bメッシュ・同一保存場から、softとmidを独立に変更する試験を先に行う。各段の履歴を保存し、局所種残差、対流・拡散寄与、再正規化の修正量、更新量とfloat32のULP比を測る。その後に化学種結合方式を一因子で比較してください。乱流拡散はコード上存在しますが、今回の原因と断定できる段階ではありません。

6. **Major — 健全性・床ゲートには検査漏れと正当な境界条件の誤検出があります。**

   **根拠:** [sern_gates.py:30](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:30) の対象には `roY*` がなく、必要データが欠けても51行以降でスキップします。最終場の種非負性・`ΣρY=ρ` はゲートに含まれません。

   [同ファイル:134](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:134) の `floor_gate` は、実configの `tMin/roMin` を読まず、読取り例外でも `ok=True` を返します。また、全ノードで `k<=0` を禁止すると、低Reのnode壁で正当に `k=0` をピンする仕様と衝突します（[solver-settings.md:250](/home/sano/work/forge/procedures/solver-settings.md:250)）。

   なお、今回調べたBの最終場では種分率範囲は `[0,1]`、和の最大誤差は約 `1.19e-7` でした。指摘は現runの種破壊ではなく、受理ゲートの保証不足です。

   **対案:** configから必要変数と実効下限を決定し、欠損・読取り不能は判定不能として不合格にする。種分率と保存量の整合を検査し、正当な壁ピンを除いた領域で床を評価する。途中のクランプ発動は、保存時点の床への張り付きとは別に記録してください。

7. **Minor — 計画の対象・検証段階・恒久記録を整理する必要があります。**

   **根拠:** [plans/README.md:51](/home/sano/work/forge/plans/README.md:51) は本体と切り出しの関係を明示しており、重複そのものは問題ありません。ただし「実装前レビュー」に対し、対象機能は既に実装・試験済みです。§1/§2相当の適用範囲と、今後変更する項目が不明確です。

   実測範囲は `m6_on` の3形状ですが、記述は「SERNでは」「生産レシピ」と広がっています。NASA-9の物性・等エントロピー処理は [accepted/tooling-nozzle-semiperfect-gas.md](/home/sano/work/forge/plans/accepted/tooling-nozzle-semiperfect-gas.md) に既存資産があります。また、[case READMEのrun一覧](/home/sano/work/forge/case/46.sern_design/README.md:23) には今回の `run_0123`～`0186` が登録されていません。

   **対案:** 本体planを正本として、対象を当面「平面node・frozen_tp・m6_on」に限定する。残作業の先頭を受理基準と500 step検証に変更し、各試験の設定・参照run・定量閾値を明記する。検証は現行手順どおりnodeを基本とし、cell・軸対称・周期への保証は付けない。run索引も同期してください。

**推奨は、CFL上限探索を止め、受理基準と初期値を修正したうえで起動レシピを再検証することです。** 優先順は、①プラトーと収束の分離・ゲート修正、②NASA-9整合IC、③同一保存場からのBの切り分け、④提案レシピと長時間参照の比較です。現時点では生産採用へ進めません。

ファイル変更・新規CFD計算は行っていません。上記はレビュー提案であり、**plan未反映**です。

指摘数: Critical 0 / Major 6 / Minor 1
