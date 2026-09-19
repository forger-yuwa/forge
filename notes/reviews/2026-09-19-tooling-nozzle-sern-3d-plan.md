# codex レビュー: tooling-nozzle-sern-3d (plan)

- **plan**: [`plans/active/tooling-nozzle-sern-3d.md`](../../plans/active/tooling-nozzle-sern-3d.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `0280cd26` (feature/sern-design)
- **codex**: effort `high`, 6.7 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M6/m1
- **focus**: §4.15.4 の 2D 有限ベース診断を踏まえ、**2 次再構成がベース角で P を 0 に潰す問題をどう解くべきか**。肩の丸め (半径・接点・変更するランプ区間の決め方) / ベース近傍の格子分布・nj_wake / ベース高さ t_base の下限 / その他の手段を比較し、次に回すべき試験を順序付きで示してほしい。bndFirstOrder は使用禁止

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

有限ベース＋適合バンドの案 (d) は維持してよいです。ただし、§4.15.4 の「ベース角の2次再構成が `P=0` を作る」という原因認定は成立していません。  
**丸めを先行させず、再構成点とリミッタ評価点の不一致、および陰的更新の破綻順序を先に検証すること**を実装着手の条件とします。

ファイルは変更していません。以下の提案は **plan 未反映**です。

**指摘一覧**

1. **Major — `P=0` は再構成破綻の証拠になっていない。しかも、その点はベース角ではない。**

   根拠: `case/46.sern_design/run_0205_r4e_2d_base_wall/res_nan_28.h5` を再集計すると、`P=0` はノード **62103 の1点**です。位置は `(x/H,y/H)=(10.7732117,2.7181569)`、ベース下端 `2.7081569`・上端 `2.7281569` の**中央**です。その点の `ro`・`roe`・`roY0`・`roY1` は既に NaN。`ro`・`roe` が有限かつ `ro>0` の点に限定した最小圧力は **584.09 Pa** でした。

   また、保存される節点 `P` と再構成した面圧力 `P_L/P_R` は別物です。等温壁では [nodeWallDirichlet_d.cu:78](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:78) が `P=ρRTw` を再代入します。ダンプには面状態も更新途中の状態もありません。

   **対案:** 記述を「2次切替後、ベース近傍で発散。最初に破綻する演算は未特定」に訂正する。正常な `res_0.h5` から、面再構成→流束→block-DPLUR→化学種更新→EOS→壁ピンの各境界で非有限・非許容状態を捕捉してください。NaN ダンプから原因の時間順序を逆算してはいけません。

2. **Major — リミッタの評価点と、node の実際の再構成点が一致していない。最優先のコード検証対象である。**

   根拠: SLAU は [convectiveFlux_slau_d.inc.cuh:134](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:134) で再構成先を**ノード間のエッジ中点**にします。一方、実行される融合リミッタは [limiter_d.cu:277](/home/sano/work/forge/solver_density_cuda/cuda_forge/limiter_d.cu:277) で**双対面重心 `pcx/pcy/pcz`**への増分を制限しています。Barth と Venkat の双方がこの経路です。

   `run_0205/sern.h5` のベース近傍では、確認した内部面114面中103面で両点が異なりました。例えば面123581では、

   - エッジ中点: `(10.7732117, 2.7169070)H`
   - 双対面重心: `(10.7936950, 2.7163737)H`

   です。重心で制限した値から、異なる位置の再構成値の有界性は保証できません。**「Barth でも落ちるのでリミッタでは直らない」は、この共通不整合を除外できていません。**

   **対案:** node のリミッタ評価点を実際の再構成点と揃える修正を、別のソルバ plan として検証する。まず現行実装で両点の制限後 `ρ/P` と近傍極値逸脱を記録し、原因への寄与を確認してください。これは今回の発散原因と確定したわけではありませんが、丸めより直接的な検証対象です。共有勾配をゼロ化する必要はありません。

3. **Major — CFL・厚み・BC の切り分けから、等温壁・2次段について結論を出している。**

   根拠: `run_0202`–`0204` は **slip ベース・SST 1次**です。CFL低減と厚み低減を試したのはこの条件だけです。`run_0205`・`0206` の **等温壁・2次**はともに `cfl_pseudo=1`、ログ上 `implicitRelax=1`・`nStepInner=5` です。こちらの低CFL試験はありません。

   また、[runner_sern.py:278](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:278) の `vehicle_kind` は**上面とベースをまとめて**変更します。「ベースだけを壁にした効果」でもありません。slip 版のベース9節点の `wall_dist` は、下端の0から上端の **0.002 m** まで増加しており、ベース全体で0になるわけでもありません。

   **対案:** `run_0205/res_0.h5` を共通初期場として、同一形状・同一BCで2次段の CFL と `implicitRelax` を一つずつ変える。step数×CFLを「同じ擬似時間」の証明にせず、局所 `dt_local`・更新量・内部反復の増幅を比較してください。既存の [update-positivity-guard plan](/home/sano/work/forge/plans/accepted/time_integration-update-positivity-guard.md) にも、床到達が陰的反復不安定の終端症状だった負の結果があります。

4. **Major — 2D 診断は、3D 生産メッシュと同じ形状ではない。形状生成式も異なる。**

   根拠: 両 run の `prepare_info.json` は次の値です。

   | 項目 | `run_0200_r4e_mesh` | `run_0205_r4e_2d_base_wall` |
   |---|---:|---:|
   | `L_ramp/H` | 10.64344 | 10.77321 |
   | `theta_e_deg` | +5.85615° | −1.49028° |
   | `L_cowl/H` | 1.60178 | 1.20000 |
   | `y_veh/H` | 3.84089 | 2.78095 |

   さらに、2D は [mesh_sern.py:288](/home/sano/work/forge/design/forge_design/meshing/mesh_sern.py:288) で既存上面に `t_base*s³` を加えますが、3D は [mesh_sern3d.py:305](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:305) でエルミート終点を `y_e+t_base` にします。下限処理前でも同じ曲線ではなく、2D側ではこの加算が後縁勾配も変えます。

   **対案:** 既存2D結果は「同種トポロジの診断」と位置づける。3Dへ進む判定には、同じ設計入力と共通の上下面輪郭関数を使った2D断面を用意する。丸めを追加する場合も、この共通化を先に行ってください。

5. **Major — `nj_wake` だけを増やす試験では、ベースの壁法線解像を改善できない。**

   根拠: [mesh_sern.py:321](/home/sano/work/forge/design/forge_design/meshing/mesh_sern.py:321) の後流分布は、ベース高さ方向、すなわち鉛直ベースに対して**接線方向**です。`run_0205` の実測は、

   - ベース接線間隔: `Δy/H ≈ 0.0025`
   - ベース直後の壁法線間隔: `Δx/H ≈ 0.08194`
   - 上側外部バンドの第一間隔: `0.02H`

   でした。後流は高さ方向に8分割される一方、最初の流れ方向間隔はベース厚の約4.1倍です。`nj_wake` 増加だけでは、この壁法線間隔は変わりません。

   また、[_geom_start:89](/home/sano/work/forge/design/forge_design/meshing/mesh_sern.py:89) により今回の後流は一様分布になります。厚みを `0.02→0.005` にすると `Δy/H` も `0.0025→0.000625` に変わるため、純粋な形状感度ではありません。

   **対案:** 固定 `t_base` で、①後縁前後の `x` 集中、②上下せん断層と後流バンドの間隔接続、③`nj_wake` の順に検証する。`t_base ≥ 3*first_wall_frac` のような数値都合の物理下限は設けず、形状寸法と格子の成立条件を分離してください。

6. **Major — 受入条件は概ね正しいが、壁処理・履歴保存・ゲート実装が条件を満たす仕様になっていない。**

   根拠: [run_0205/solverConfig.yaml:16](/home/sano/work/forge/case/46.sern_design/run_0205_r4e_2d_base_wall/solverConfig.yaml:16) は `wallTreatmentSST: 1`、すなわち壁関数です。等温no-slipへの変更だけで低Re壁解像が保証されるわけではありません。ベースの局所壁解像と壁関数の適用妥当性を判定する条件が必要です。

   現在の [floor_gate:167](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:167) は最終保存場を調べるため、§4.15.3の「持続する床／更新クリップなし」を確認できません。3Dの `require_residual_pass` も [runner_sern3d.py:348](/home/sano/work/forge/design/forge_design/evaluate/runner_sern3d.py:348) では既定 `False` です。

   **対案:** 生産では残差PASSを必須化し、床・更新クリップの時系列、ベース圧・再循環長、局所壁解像を受入表に追加する。壁解像は `ypls` ではなく第一内部点距離と接線壁応力で確認する。2次切替直前と切替後の履歴を保存し、段階の異なる残差を混ぜないこと。ソルバ修正を伴う場合はSERNだけでなく、`procedures/verification/README.md` に従う標準ケース回帰も追加してください。

7. **Minor — 残作業表と現行コードの説明が同期していない。**

   根拠: planのR5bは `k≤0` 判定と残差列中央値による拒否を実装済みとしていますが、現行 [sern_gates.py:145](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:145) では前者を撤去、[同:223](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:223) では後者を補助警報へ降格しています。また、`run_0205` のsoft終端 `rms_ro` は実CSVでは **7.00e-7** で、表の `1.0e-5` と異なります。

   **対案:** §5.1を実行順に並べ直し、完了した実装と未完了の受入条件を区別する。起動planが持つ判定仕様への参照を一本化してください。

**再確認した検証結果**

`check_convergence.py` を再実行し、以下はすべて **`DIVERGED (NaN/Inf)`** でした。パスはすべて `case/46.sern_design/` 配下です。

| run | NaNダンプ |
|---|---|
| `run_0202_r4e_2d_base` | `res_nan_7.h5` |
| `run_0203_r4e_2d_base_cfl01` | `res_nan_57.h5` |
| `run_0204_r4e_2d_base_t005` | `res_nan_20.h5` |
| `run_0205_r4e_2d_base_wall` | `res_nan_28.h5` |
| `run_0206_r4e_2d_base_wall_barth` | `res_nan_27.h5` |

`run_0205` はメッシュ品質 **`VERDICT: PASS`**（AR最大677.1、skew最大0.439）ですが、`check_quasisteady.py` は **`TRANSIENT-UNSETTLED (only 1 res_*.h5)`**。3Dメッシュ試験は再実行して **`ALL PASS`** でした。接続設計の支持にはなりますが、CFDの成立証明にはなりません。run索引は [case README](/home/sano/work/forge/case/46.sern_design/README.md:253) です。

**推奨は一つ: 案 (d) を維持し、固定形状で原因を特定してから形状感度へ進む。**

次の試験順を推奨します。

1. **現行失敗の再現と演算段階の特定。** `run_0205/res_0.h5` を複製runの共通初期場にし、1次継続と2次切替を比較。最初の異常まで毎step保存し、問題点周辺では内部反復も捕捉する。
2. **再構成点の整合試験。** 現行と、nodeのリミッタ評価点をエッジ中点に揃えた実装を比較する。有限な入力から負の面状態が実際に出る場合に、面状態の許容性を保つ再構成制限を検討する。有界再構成と正値性の関係には一次資料がありますが、forgeのTP・陰的更新全体の正値性まで保証するものではありません。[NASA技術報告](https://ntrs.nasa.gov/citations/19990047773)
3. **同じ等温壁・2次条件で時間積分を分離。** `cfl_pseudo=1→0.1`、次に `implicitRelax=1→0.7` を個別比較。SSTが先に異常化する場合だけ、乱流更新凍結を診断対照に加える。単なる延命を成功扱いしない。
4. **固定形状の局所格子試験。** `t_base/H=0.02` を維持し、まず後縁前後の `Δx` を縮める。その後、上下接続の間隔比を滑らかにし、`nj_wake=9→17` を比較する。格子変更時は場を移植して適応段を置く。
5. **最後に肩の丸めと厚みの形状感度。** 丸めは次の仕様を先に固定する。

   | 候補 | 採用条件・試験仕様 |
   |---|---|
   | 肩の丸め | 原輪郭と鉛直ベースへ接する円弧を固体側に置く。上下の半径、接点、下側ランプの変更区間を明記する |
   | 半径 | 暫定比較は `r/H=0、0.0025、0.005`。物理入力として固定し、円弧は少なくとも8区間で解像して再細分化する |
   | ベース残存高さ | 各円弧がベース上で消費する長さを `d_lower,d_upper` として、`t_flat=t_base−d_lower−d_upper>0` を必須にする。非直角では単純に `t_base−2r` としない |
   | `t_base` | `0.01/0.02/0.04H` は別の形状感度。現状から普遍的な下限は決められない。薄くして安定化する方針は採らない |

丸めた下側肩はノズル輪郭そのものを変更します。変更区間の力もノズル帳簿へ含め、従来の「鉛直ベース面積=`t_base×幅`」という試験は、円弧を含む新しい面集合に合わせて更新する必要があります。

指摘数: Critical 0 / Major 6 / Minor 1
