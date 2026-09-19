# codex レビュー: tooling-nozzle-sern-chain (plan)

- **plan**: [`plans/active/tooling-nozzle-sern-chain.md`](../../plans/active/tooling-nozzle-sern-chain.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `8b2668a6` (feature/sern-design)
- **codex**: effort `high`, 7.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m1
- **focus**: §4.15.2 の R4e 案 (d) (機体後縁に有限の厚みを持たせてベース面で終わらせる) と案 (e) (バンド上端を直線にして出口まで通し上バンド下端を z 依存にする) のどちらを採るべきか。論点 1-4 と、2D run_0034 の base+wake 発散が 3D で再発しないかを見てほしい。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
**案 (d) を推奨します。ただし、有限ベースと後流ブロックを一体で設計することが条件です。**  
`run_0034` は有限ベースの不成立を証明していません。実行ログは「soft 段」ではなく、2 次・`cfl_pseudo=4` でした。

ファイルは変更していません。目的である「幅外の人工固体・閉じ壁の除去」は未解決で、既存の accepted plan と重複しません。既存の node 壁境界・陰解法を再利用する方針は妥当です。ただし、今回の修正で 3D SST 評価器全体が成立するとはまだ言えません。

再判定した根拠は次のとおりです。

| 対象 run（`case/46.sern_design/` 配下） | メッシュ品質 | `check_convergence.py` | `check_quasisteady.py` |
|---|---|---|---|
| `run_0034_cycle_euler_m4off_exttop/` | `PASS`、AR 443.4、skew 0.464 | **DIVERGED (NaN/Inf)** | 有効な定常判定なし |
| `run_0122_r5_outflow/` | `PASS`、AR 789.2、skew 0.401 | **NOT CONVERGED (stalled/plateau)** | 力係数3成分は **STEADY**、`pmax` は **DRIFTING** |

保存先・履歴の索引は [case README](/home/sano/work/forge/case/46.sern_design/README.md:23) です。後者の最終場では `Pmax=5,721,807 Pa`、`rms_roOmega≈8.95e15`。性能評価には使えません。

**指摘一覧**

1. **Major — 案 (d) の論点1は、後縁でブロック構成を切り替えるという不要な前提から生じています。**

   **根拠:** [plan:531](/home/sano/work/forge/plans/active/tooling-nozzle-sern-chain.md:531) は、後縁下流で上バンド下端を旧ランプ線へ落とすとしています。現行コードも [mesh_sern3d.py:315](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:315) で上バンド下端をランプ線へ切り替え、側面バンドは [同:413](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:413) で閉じています。実データでも、幅外を塞ぐ **224 面**を確認しました。

   **対案:** 案 (d) では、幅外バンドを出口まで維持してください。下端を `a(x)=yt(x)`、上端を `b(x)` とし、後縁で
   \[
   b(L_{\rm ramp})=y_e+t_{\rm base}
   \]
   とします。下流の `b(x)` は後流ブロック上端として連続に延長します。

   - 幅外：`a..b` の流体バンドを全長に置く。
   - 幅内：同じ `a..b` のブロックを後縁下流だけに置く。
   - 上バンド：全幅で `b..天井` を担当する。
   - ベース壁：`x=L_ramp`、幅内の `a..b` だけに置く。

   共通の節点列 `y_j=a+η_j(b−a)` と**同じ節点ID**を使えば、側面バンドと後流の接続は適合します。上バンドへの担当変更がなく、案 (e) の z 依存の分割線も不要です。メモリ内の接続骨格検算では、92 hex、非多様体面0、流体域連結、幅外の閉じ面0を確認しました。これは接続の成立確認であり、CFD安定性の実証ではありません。

2. **Major — `run_0034` の実行条件と原因診断が誤記され、案の選択根拠を歪めています。**

   **根拠:** [forge_run.log:13](/home/sano/work/forge/case/46.sern_design/run_0034_cycle_euler_m4off_exttop/forge_run.log:13) は `cfl_pseudo=4`、[同:24](/home/sano/work/forge/case/46.sern_design/run_0034_cycle_euler_m4off_exttop/forge_run.log:24) は `convMethod=1`、すなわち2次です。`res_nan_5.h5` の密度非有限値79点がベース近傍に集中することは確認できましたが、これだけでは角の処理を真因と特定できません。

   現行の [slip_d:68](/home/sano/work/forge/solver_density_cuda/cuda_forge/boundaryCond_d.cu:68) は面ごとの ghost・境界値を構成します。「共有ノードを二方向へ順次射影するので必ず破綻する」という説明を、そのまま裏付ける実装ではありません。

   **対案:** 「soft 段」「90°二重 slip が原因」を撤回し、**ベース近傍で発散、根因未確定**に直してください。同じ有限ベース形状について、現行バイナリ・1次・低 `cfl_pseudo` の診断を先に置きます。

   3Dでも再発リスクはあります。ただし、生産の等温 no-slip は `u=0` を共通に課し、[timeIntegration_d.cu:1013](/home/sano/work/forge/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1013) で運動量行を拘束するため、旧 Euler/slip と同一条件ではありません。**肩の丸めは検証候補であって、発散回避の保証ではありません。** 丸めを採るなら半径・接点・変更するランプ区間まで定義してください。

3. **Major — 物理形状が格子間隔に依存しています。これを残すと格子独立性試験が成立しません。**

   **根拠:** [mesh_sern3d.py:274](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:274) は
   `clearance=max(vehicle_clearance,3*first_top_frac)`、
   [同:288](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:288) は `first_wall_frac` で機体終端を決めます。

   `run_0122` の輪郭を使った関数検算では、`first_top_frac` を `0.02→0.01` にするだけで、機体上面基準高さが **3.840894H→3.810894H**、つまり **0.03H** 動きました。

   **対案:** `t_base/H`、機体上面輪郭、丸め半径を物理入力として固定し、格子不足ならメッシュ生成を失敗させてください。格子に合わせて形を動かしてはいけません。

   `vehicle_clearance=0.02H` は `t_base=0.02H` の根拠になりません。構造寸法が未定なら、例えば `t_base/H=0.02` を**暫定モデル値**として明記し、`0.01/0.02/0.04` の形状感度を別試験にします。ベースの壁法線は x 方向なので、厚み方向の分割数だけでなく、直後の Δx・側端の Δz・実測 y⁺も必要です。

4. **Major — ベースの境界条件と帳簿を、既存タグへの追加だけで実装すると誤集計します。**

   **根拠:** [runner_sern3d.py:75](/home/sano/work/forge/design/forge_design/evaluate/runner_sern3d.py:75) は機体上面を slip 固定にしています。また [forces3d:243](/home/sano/work/forge/design/forge_design/evaluate/runner_sern3d.py:243) の集計対象に `vehicle_side` はなく、[同:262](/home/sano/work/forge/design/forge_design/evaluate/runner_sern3d.py:262) は特別扱いしない追加面をノズル力へ加算します。新しい `vehicle_base` を単純追加すると同じ罠になります。

   **対案:** 生産では `vehicle_top/side/base` を同一の `wall_isothermal`・`Tw` に統一し、変換時から壁距離を生成してください。R4fを後回しにせず、案 (d) の生産仕様に含めます。slip は診断用に限定します。

   集計は明示的な面集合で分けてください。

   - ノズル力：ランプ・カウル・ダクト側壁。
   - 機体力：上面・機体側面・ベース。
   - 全制御体積の収支：両者を含む全境界。

   ベースをノズル力から除外する帳簿は妥当ですが、**機体力として保存する必要があります**。また、帳簿から除外してもベース流れによるノズル圧力への影響は残ります。

5. **Major — §6には、今回の失敗を検出する定量的な検証条件が不足しています。**

   **根拠:** 現行の `run_sern_mesh3d_tests.py` は今回も **ALL PASS** でした。しかし [同:87](/home/sano/work/forge/design/tests/run_sern_mesh3d_tests.py:87) は機体上面の z 座標を中央値で検査し、[同:114](/home/sano/work/forge/design/tests/run_sern_mesh3d_tests.py:114) は角の行列式の**絶対値**しか検査しません。人工壁や負向き要素を十分に排除できません。

   また `run_0122` は力係数が `STEADY` でも、`pmax` は末尾区間で **3.1% drift、DRIFTING**。有限値・非上昇だけを認める [sern_gates.py:104](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:104) の限界も残っています。

   **対案:** 実装前に、検証を次の順で固定してください。

   1. **接続:** 内部面は同じIDで2セル共有、境界面は1セル所有、重複タグ0、幅外の横断壁0。ベース面積は設計値と照合する。
   2. **幾何:** 符号付き Jacobian、双対体積・閉性、float32変換後の節点衝突、メッシュ品質 `PASS`。
   3. **起動:** 2D有限ベース診断→3D層流→SST。追加流体には領域別の組成・状態を与え、変更前メッシュからの無条件indexコピーを禁止する。
   4. **受入:** R5bを先に有効化し、全残差 `PASS`、目的量 `STEADY`、持続する内部の床・更新クリップなし。壁の正当な `k=0` ピンは異常クリップと区別する。
   5. **独立性:** 生産3作動点で、固定形状の格子・領域感度を測る。暫定基準として `|ΔC_T|≤0.002`、`|ΔC_L|≤0.002`、`|ΔC_M|≤0.02` を明記し、設計上の識別幅が小さければさらに厳しくする。

   ベース圧・再循環長も時系列監視対象にしてください。定常擬似時間で振動する場合は、その振動を物理的な後流変動と解釈できません。

6. **Minor — 現行仕様と履歴が混在し、今回の実装範囲が曖昧です。**

   **根拠:** §4.2の5変数・2目的に対し、[§4.8:193](/home/sano/work/forge/plans/active/tooling-nozzle-sern-chain.md:193) は4目的・`10d=60` のままです。Raoの [テスト6b:135](/home/sano/work/forge/design/tests/run_sern_moc_tests.py:135) は等長拘束下の局所比較を追加していますが、Shyneの全条件との独立照合や3D粘性最適性の証明にはなりません。

   **対案:** 現行仕様を冒頭で一本化し、旧記述は履歴と明示してください。今回の対象は「R4eの領域修復と評価器成立の前提整備」。MOCは形状パラメータ化として扱い、SST交差拡散などソルバ変更が必要なら既存の関連planまたは別planへ分離します。検証は現行の [verificationルール](/home/sano/work/forge/procedures/verification/README.md:40) に従い node 主体で組み、cellへ切り替えて成立扱いにしないでください。

**推奨は、案 (d)「全長機体＋有限ベース＋適合した側面・後流バンド」です。**

| 比較軸 | 案 (d) | 案 (e) |
|---|---|---|
| 実装量 | 共通断面分布で側面・後流・上バンドを接続できる | z依存の分割線と機体上面との接続設計が追加で必要 |
| 数値面 | ベース近傍の剥離・再循環を検証する必要あり | ベース追加は避けられるが、現案では後縁の縮退処理が未解決 |
| 物理面 | 有限厚の機体という今回の要件に直接対応 | 薄い後縁モデルのままで、その要件を満たさない |

実装前の優先順は、**①接続・物理寸法の確定 → ②ベースを含むBCと帳簿の確定 → ③`run_0034` の診断訂正と起動試験の定義 → ④R5b・定量検証条件の確定**です。この修正後なら実装に進めます。案 (d) が発散しないという保証は、現時点ではありません。

本レビューの提案は **plan 未反映**です。読み取り専用の依頼に従い、文書・コード・runは変更していません。

指摘数: Critical 0 / Major 5 / Minor 1
