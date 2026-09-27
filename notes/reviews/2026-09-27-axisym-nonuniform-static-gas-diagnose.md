# codex 諮問 (diagnose): axisym-nonuniform-static-gas

- **brief**: [`notes/reviews/briefs/2026-09-27-axisym-nonuniform-static-gas.md`](../../notes/reviews/briefs/2026-09-27-axisym-nonuniform-static-gas.md)
- **plan**: [`plans/active/boundary-cht-axisymmetric-fem2d.md`](../../plans/active/boundary-cht-axisymmetric-fem2d.md)
- **date**: 2026-09-27
- **commit**: `8e4e7ef3` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 3.6 min, rc=0
- **結論**: **同じ非一様格子・同じ FP64 バイナリで、`mesh.hoopAreaFromClosure: 0 / 1` だけを変える、連成前 500 step の A/B を行う。**
- **extra**: `case/62.conjugate_disk/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **同じ非一様格子・同じ FP64 バイナリで、`mesh.hoopAreaFromClosure: 0 / 1` だけを変える、連成前 500 step の A/B を行う。**

第 1 仮説: **入力メッシュに残る FP32 の幾何閉包欠損が偽の圧力力を供給し、非一様格子で擾乱が増幅される。** 確度: **中**  
  根拠:
  - ローカルの `case/62.conjugate_disk/mesh/disk_r32_{u,g1p1}.h5` は、`volume`・`centCoords`・`surfVect` がすべて **float32 保存**。FP64 で読み込んでも入力時の丸めは戻らない。
  - 実ファイルを倍精度で読み、各 CV の符号付き面和から
    \[
    E_i=\sum_f r_fS_{if}-(0,A_i,0)
    \]
    を再計算した。**境界を除く内部 CV** の最大値は以下。ソルバ実行結果ではなく、入力幾何の直接検査である。

    | 入力格子 | `max|E_x|/A` | `max|E_r|/A` |
    |---|---:|---:|
    | 一様 r32 | 0 | `1.404e-6` |
    | 非一様 r32 | `5.960e-6` | `3.115e-6` |

  - 非一様では**内部にも軸方向の閉包欠損がある**。面を `r_face` 倍する処理と、元の面積を hoop に使う処理は別である（[`variables.cpp:540`](/home/sano/work/forge-cht/solver_density_cuda/variables.cpp:540)、[同:568](/home/sano/work/forge-cht/solver_density_cuda/variables.cpp:568)）。
  - `hoopAreaFromClosure: 1` は同じ面ベクトルから閉包を作り、半径方向ソースの置換と軸方向圧力補正を有効にする（[同:575](/home/sano/work/forge-cht/solver_density_cuda/variables.cpp:575)、[`axisymmetricSource_d.cu:111`](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/axisymmetricSource_d.cu:111)）。読み取り階層は **`mesh`**（[`solverConfig.cpp:222`](/home/sano/work/forge-cht/solver_density_cuda/input/solverConfig.cpp:222)）。
  - **閉包欠損の存在は確認できたが、それが速度差約 500 倍・壁熱流束の市松を生む主因かは未確認**。一様格子にも欠損があり、欠損量だけで症状の差を説明してはいけない。

  反証条件: AWS 入力でも欠損とスイッチの有効化を確認したうえで、補正側でも基準側と同程度の市松・半径速度が発生すれば、「閉包欠損が支配的な駆動源」を棄却する。

第 2・第 3 仮説:
- **第 2〔未確認〕:** 加熱起動と低 Mach の再構成・陰的更新の組合せが、閉包欠損とは独立に市松を増幅する。初期温度は 325 K だが、熱壁の状態ピンは密度を保持して `T`・`P`・`roe` を変更するため、起動直後は静止平衡ではない（[`nodeWallDirichlet_d.cu:80`](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:80)）。
- **第 3〔未確認〕:** slip 端・等温壁との角の境界処理。境界質量流束や初発位置の証拠がなく、除外も主因認定もできない。

判別 A/B:  
  **A = `mesh.hoopAreaFromClosure: 0`、B = `1`。** `run_0014` の入力を別々の新規 run に複製し、同じ IC・バイナリ・BC・CFL・再構成で各 **500 step**。両方とも `warmup: 20000` を維持し、連成を開始しない。出力間隔は両方共通で 10 step とする。

  見る量は全保存量残差、`max|Uy|`、同じ軸位置での半径方向圧力差、壁の `iface_q_eff` の市松振幅。市松振幅は「各内部壁節点の値と、左右隣接節点から半径座標で線形補間した値との差」の最大値で定義する。warmup 中の伝導基準は **`0.0241×(350−325)/0.005 = 120.5 W/m²`**。

  → **A で異常が再現し、B の半径速度・市松振幅の両方が step 250–500 の最大値で A の 1/10 以下**なら、第 1 仮説を支持する。  
  → **B でも両方が A の 1/2 以上**なら、閉包補正では症状が消えないため第 1 仮説を主因候補から下げ、第 2 仮説を優先する。  
  → 中間、または A が再現しなければ判別不能。**この短時間比較は発生機構の診断であり、収束・定常化・V-ax2b 合格の判定には使わない。**

やらない方がよいこと:
- 平面化・`pRef`・CFL・次数を同時に変更すること。今回の閉包仮説には上記スイッチが直接対応する。
- 固体への荷重再配分、非一様格子の撤回、warmup の再延長だけで処置すること。
- 現コードで修正済みの `pRef`–hoop ゲージ不整合を再び真因とすること。ソースは既に `P−pRef` を使う（[`axisymmetricSource_d.cu:51`](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/axisymmetricSource_d.cu:51)）。

呼び出し側の前提への異議:
- **Major — 「FP64 だから幾何の丸めは除外できる」は成立しない。** 上記入力ファイルが反証。対案は AWS 側入力の同一性・保存型も確認し、閉包補正の A/B で因果を調べること。
- **Major — 総質量指標の低下から境界漏れを推論できない。** `unsteady: 0` は局所擬似時間刻みを使い、block-DPLUR の時間対角も `V_i/Δτ_i` である（[`setDT_d.cu:362`](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/setDT_d.cu:362)、[`timeIntegration_d.cu:827`](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/timeIntegration_d.cu:827)）。反復途中の `ΣρV` は物理時間の保存試験にならない。対案は、離散連続式の残差総和と境界質量流束を照合すること。質量そのものは**入力メッシュの双対重心半径**を用いる `Σρ A_planar r̄`。ノード半径への置換は不可（[`mesh.cpp:363`](/home/sano/work/forge-cht/solver_density_cuda/mesh/mesh.cpp:363)）。
- **Major — 「流体側の異常」は妥当な切り分け候補だが、「軸対称の流体実装が原因」は未確定。** 推奨は、**流体診断・修正を別 plan に切り出し、CHT plan では V-ax2b の依存障害として追跡する**こと。§5.1 #6・§6 V-ax2b に「`run_0014` 安全停止、非一様検証未達、原因未確定」と記録し、失敗履歴・登録閾値・マージ禁止を維持する（[plan:142](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:142)）。**本回答は plan 未反映**。反映は呼び出し側で行う。

不足情報:  
AWS の `run_0013`・`run_0014` の実効 YAML、入力メッシュのハッシュ、FP64 変更差分・バイナリ識別情報、全残差 CSV、連成開始前までの複数スナップショット、判定区間付き `CONVERGENCE_VERDICT.txt`・準定常 VERDICT。これらは未取得のため、報告された run 数値は独立検証できておらず、一様側も「正常に収束」とは認定しない。
