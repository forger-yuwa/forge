# codex レビュー: tooling-nozzle-sern-chain (plan)

- **plan**: [`plans/active/tooling-nozzle-sern-chain.md`](../../plans/active/tooling-nozzle-sern-chain.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `41068a06` (feature/sern-design)
- **codex**: effort `high`, 8.3 min, rc=0
- **判定**: **NO-GO**, 指摘 C2/M5/m1
- **focus**: §4.15 の 3 件に集中。(1) rms_roOmega が 8.95e15 で張り付く件 — 入口面 x/H -0.5 と側壁 z/H 1.0 の稜線が疑わしい。BC 優先順位と残差集計のどこを読むべきか、切り分けの手順。(2) vehicle_top の probe 出力と、slip のままか粘性壁にするかの判断。(3) 対称面の鏡像フィルタ — 反転すべき成分の網羅に漏れがないか。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
巨大な `rms_roOmega` は、入口稜線よりも「対称面上の1ノードでの下限張り付きと交差拡散」で説明できます。  
さらに幅外に人工的な壁が残っています。§4.15 の診断方針と評価ゲートを直してから進めるべきです。

ファイルは変更していません。コード、既存の accepted plan、`case/46.sern_design/run_0122_r5_outflow/` の全26保存場・残差履歴を確認しました。再判定結果は次のとおりです。

| 検査 | 結果 |
|---|---|
| `check_convergence.py` | **NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)** |
| `check_quasisteady.py`：`C_T_with_shear,C_T,C_L,C_M` | **OVERALL: ALL STEADY** |
| 全26保存場の `VALUE/*` | NaN/Inf なし |
| 保存済み `MESH_QUALITY.txt` | **PASS**：AR 789.2、skew 0.401 |
| 現行 `metrics.json` の評価ゲート | **PASS**――以下の問題を見逃している |

1. **Critical — 巨大残差の診断が誤った方向を向き、ゲートが異常解を受理している**

   `roOmega` の最大値が一定でも、残差の発生場所は分かりません。実際、残差は開始時から張り付いていません。

   ```
   step 1450: rms_roOmega = 0.223809
   step 1451: rms_roOmega = 1.235441
   step 1452: rms_roOmega = 8.28907e15
   ```

   `res_1500.h5` 以降、**ノードID 646932** が `roOmega=1e-20`、`k=0` に張り付いています。位置は **(x/H,y/H,z/H)=(10.5397,3.76611,0)**、所属は **`sym` のみ**です。入口でも側壁でもありません。

   `res_12000.h5` とメッシュから、実装と同じ面補間による勾配・交差拡散項を独立再計算すると、

   - `∇k·∇omega = −1.07048e16`
   - `k=0` により `F1=0`
   - 実装の `omega` 分母下限 `1e-12` を使った **交差拡散項×体積 = −9.58791e18**
   - 全 **1,146,866** ノードに対する、この1点だけの RMS 寄与 = **8.95298e15**

   となり、報告値をほぼ説明します。**保存場が有限なのは、負の更新を下限で切っているため**です。最初に下限へ落ちる原因は未確定ですが、その後の巨大残差の維持機構は入口ピン説よりはるかに強く裏付けられています。

   根拠は [ransSource_d.cu:136](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransSource_d.cu:136)、[同:212](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransSource_d.cu:212)、更新クリップの [update_d.cu:343](/home/sano/work/forge/solver_density_cuda/cuda_forge/update_d.cu:343)。一方、[sern_gates.py:104](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:104) は有限かつ rising でなければ受理し、今回も `require_residual_pass=false` で **PASS** です。

   **対案・切り分け手順：**

   - R5 の受入では `require_residual_pass` を必須にし、R5b に **`roOmega` 下限到達と更新クリップ**も追加する。`rms_roOmega` を除外して通してはいけません。他の保存量も低下は約1～2桁で、未収束です。
   - 同一メッシュの `res_1000.h5` から新規診断 run を作り、step 1452 相当までを細かく保存する。`res_roOmega`、`omg_trans/prod/dest/cross`、`dK*`、`dOmega*`、`sstF1`、`src_jac_omega`、`transport_diag_omega`、`dt_local`、更新前後の `roK/roOmega` を取得する。
   - 閾値超過**点数だけでなく二乗和寄与率**で局在化し、下限到達前の輸送・負の交差拡散・更新量を追う。ソルバ変更が必要なら、正値性を保つ SST 更新を別 plan にする。

   入口側も確認する場合、読む順序は次です。

   | 責務 | 確認先 |
   |---|---|
   | 壁所属フラグ | [mesh.cpp:922](/home/sano/work/forge/solver_density_cuda/mesh/mesh.cpp:922) |
   | BC 適用順序 | [boundaryCond.cpp:450](/home/sano/work/forge/solver_density_cuda/boundaryCond.cpp:450) |
   | 壁・入口の保存量ピン | [ransBoundary_d.cu:53](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransBoundary_d.cu:53)、[同:116](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransBoundary_d.cu:116) |
   | SST 残差除外 | [ransSource_d.cu:268](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransSource_d.cu:268) |
   | 最終壁残差射影→集計 | [main.cpp:1383](/home/sano/work/forge/solver_density_cuda/main.cpp:1383)、[residualMonitor_d.cu:124](/home/sano/work/forge/solver_density_cuda/cuda_forge/residualMonitor_d.cu:124) |

   入口ピンと残差除外は [既存 accepted plan](/home/sano/work/forge/plans/accepted/turbulence-node-inlet-dirichlet-conserved.md) で実装済みです。BC は `msh.bconds` 順の上書きですが、それを今回の巨大残差の原因とする根拠はありません。

2. **Critical — R4c は未完了。幅外流体の下流端を壁で塞いでいる**

   [mesh_sern3d.py:407](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:407) は、追加した流体バンドの末端を `vehicle_side` として閉じます。この面は **`z=W/2` に限定されず遠方境界まで延びます**。

   `run_0122` の実メッシュには、その壁が **224面**あり、範囲は

   ```
   x/H = 10.5397
   y/H = 3.77011～3.78617
   z/H = 1.0～2.5
   ```

   です。最終保存場の最大圧力 **5.721807 MPa** は、この壁と `side_far` の交線にあります。生産条件の外気圧は **2851 Pa**。幾何品質 PASS は、意図した流体領域ができた証拠になっていません。

   **対案：** 幅外の追加流体を下流流体へ接続する末端ブロックを設計し、内部面として共有する。本物の機体ベースを設けるなら機体幅内に限定する。断面・タグ試験に「幅外に横断する固体面がない」「追加流体から下流への接続がある」を追加し、R4c を再検証してください。

3. **Major — `outflow` の効果を「出口だけの修正」と断定できない**

   [runner_sern3d.py:47](/home/sano/work/forge/design/forge_design/evaluate/runner_sern3d.py:47) の共通関数は `outlet` と `bottom` の両方に使われています。実際、`run_0121_r5_vehicleside_lsw08` → `run_0122_r5_outflow` では、**両境界が同時に** `outlet_statPress` から `outflow` へ変わっています。

   最終保存場の `bottom` の面法線 Mach 数は **−5.20e-6～4.31e-6**。主流 Mach が高くても、この面は超音速流出面ではありません。全量外挿の適用根拠が不足しています。また、2851 Pa と6.5 kPaの差は約2.3倍で「桁違い」でもありません。

   **対案：** 下流出口と遠方境界の設定を分離し、修正後の領域で出口だけを変更した比較を行う。生産3作動点について面法線 Mach、境界距離感度、運動量収支を確認し、`C_T/C_L/C_M` ごとの許容差を事前に固定してください。

4. **Major — `probe.yaml` の `surfaces:` を埋めても面プローブは動かない**

   ビルド対象は [probe/CMakeLists.txt:3](/home/sano/work/forge/solver_density_cuda/probe/CMakeLists.txt:3) の `point_probes.cu` です。その [readYAML():64](/home/sano/work/forge/solver_density_cuda/probe/point_probes.cu:64) は **`points` しか読みません**。`surfaces` の読取・出力実装はありません。

   一方、`run_0122` には既に `res_vehicle_top_15_<step>.h5` があり、圧力・温度・速度を保持しています。面出力機能を新設する必要はありません。

   **対案：** 既存の面 HDF5 系列から、タグと面接続に基づく固定位置・面積重み統計を CSV 化する後処理に絞る。Mach 数は対応する体積場の `sonic` を使い、凍結 TP を一定 γ で近似しない。定常計算なので横軸は擬似時間の **step** と明示し、CSV を正式な準定常判定器へ接続してください。

5. **Major — 「帳簿外だから力には効かない」は誤り。上面の壁モデルを決める根拠にならない**

   [plan:483](/home/sano/work/forge/plans/active/tooling-nozzle-sern-chain.md:483) の主張は、**直接の面積分対象**と**流れ場への影響**を混同しています。上面境界層・熱伝達・後縁流れが変われば、帳簿に含まれるランプやカウルの圧力も変わり得ます。

   現在は [runner_sern3d.py:75](/home/sano/work/forge/design/forge_design/evaluate/runner_sern3d.py:75) で上面だけ `slip`、機体側面は等温粘性壁です。この組合せのモデル誤差は測定されていません。

   **対案：** **生産評価の `vehicle_top` は、機体側面と同じ等温粘性壁に統一する**ことを推奨します。暫定壁温1000 Kも明記する。現在の slip は比較用として残し、領域・SST 更新の問題を解消してから同一幾何で差を測る。上面の壁距離を再生成し、近壁解像度と壁関数の適用状態も検証対象にしてください。

6. **Major — 鏡像の符号規則は型で定義すべき。名前の末尾判定では誤変換する**

   レビュー中に確認した未追跡の [mirror_symmetry.py:32](/home/sano/work/forge/solver_density_cuda/tools/mirror_symmetry.py:32) も点検しました。実行による関数確認では、

   - `limiter_Uz` を反転する――**誤り**。これは速度ではなく制限係数。
   - `dq_block_old_3`、`rhs_block_3` を反転しない――z 運動量成分なので漏れ。
   - `diag_block_03` を反転しない――流れ変数の変換則を適用するなら漏れ。
   - `dUzdz` を反転しない――こちらは正しい。

   **対案：** 変数の意味を明示した登録表を作る。`R=diag(1,1,−1)` として、極性ベクトルは `Rv`、速度勾配は `RGRᵀ` とする。

   | 対象 | z 鏡像での扱い |
   |---|---|
   | `Uz, roUz, twall_z`、スカラーの z 勾配 | 反転 |
   | `dUxdz,dUydz,dUzdx,dUzdy` | 反転 |
   | `dUzdz`、組成、SST の `omega`、制限係数 | 維持 |
   | 渦度ベクトル | 軸性ベクトルなので **x,y を反転、z を維持** |
   | block 補正・残差・行列 | 保存変数の符号行列から変換 |

   未知フィールドは黙って複製せず検出する。VTK 自体もベクトル・テンソルの反転を提供しており、「Reflect は符号を扱えない」という説明も不正確です。ただし forge の分離されたスカラー配列には型付けが必要です。[VTK公式仕様](https://vtk.org/doc/nightly/html/classvtkReflectionFilter.html)

7. **Major — 鏡像の接続処理が hex の向きを壊す。対称面での値の整合確認もない**

   [mirror_symmetry.py:84](/home/sano/work/forge/solver_density_cuda/tools/mirror_symmetry.py:84) の `ids[::-1]` は、3D要素一般の向き修正になりません。単位 hex で同じ操作を検算すると、**符号付き Jacobian が +1 → −1** になります。

   また [同:91](/home/sano/work/forge/solver_density_cuda/tools/mirror_symmetry.py:91) は節点数と違う配列を無条件に複製するため、cell/面中心ベクトルの符号反転が抜けます。保存時には `h0_includes_k` 属性も失われます。

   対称面の共有も、座標だけでは不十分です。`run_0122/res_12000.h5` の `sym` ノードでは **`Uz=−239.85～329.22 m/s`**。これをそのまま共有して「対称な全幅場」と扱うと、奇関数成分の継ぎ目が不整合になります。

   **対案：** 要素種別ごとの向き反転置換、XDMF の Node/Cell 所属、HDF5属性を保持する。対称面は元ID単位で共有し、既存の座標一致双子ノードを一括統合しない。検証には正 Jacobian、二重反転、属性保存、反対称成分の面上誤差を含める。面上の非零 `Uz` は可視化側で黙ってゼロ化せず、元解の対称条件誤差として報告してください。

8. **Minor — 現行仕様・残作業表・検証条件が同期していない**

   目的である「成立した3D評価器を先に作る」は妥当です。ただし次が混在しています。

   - §4.15 の3件に対応する作業・受入条件が §5.1 に整理されていない。
   - 「2D の入口角問題は解決済み」は、同 plan の R4b「測定して却下」と矛盾する。
   - `nodeInletCornerWall` 自体は **3Dにも実装済み**。[solver-settings.md:293](/home/sano/work/forge/procedures/solver-settings.md:293)
   - §6 の「全作動点 PASS」と §4.7 の「PASS 任意」が不一致。
   - [case README の run 索引](/home/sano/work/forge/case/46.sern_design/README.md:23) に、確認時点で `run_0122` の行がない。

   **対案：** §5.1 を現行作業の正本として更新し、ソルバ変更は別 plan と関連付ける。回帰計画は現行の [verification/README.md](/home/sano/work/forge/procedures/verification/README.md) に従い node で組み、SST・等温壁変更には `case/48` を含める。鏡像は可視化用に限定し、周期境界・軸対称の再構成や restart への対応を暗黙に含めないでください。

**推奨は、3D評価器の成立確認に作業を戻すことです。** 優先順は、**異常解を落とすゲート → 幅外末端の流体接続 → step 1452 前後の SST 下限到達診断 → 境界条件・上面粘性壁の検証 → 面時系列・鏡像ツール**です。生産3作動点で全残差と対象量の判定が揃うまで、MOO・3D性能比較には進めません。

これはレビュー提案であり、ユーザー指定の read-only に従って **plan 未反映**です。

指摘数: Critical 2 / Major 5 / Minor 1
