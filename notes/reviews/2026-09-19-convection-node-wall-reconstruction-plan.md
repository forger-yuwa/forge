# codex レビュー: convection-node-wall-reconstruction (plan)

- **plan**: [`plans/active/convection-node-wall-reconstruction.md`](../../plans/active/convection-node-wall-reconstruction.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `cfc221cf` (feature/sern-design)
- **codex**: effort `high`, 6.7 min, rc=0
- **判定**: **NO-GO**, 指摘 C1/M6/m1
- **focus**: SU2 v8.5.0 のソースと突き合わせた結果 (§4.4/§4.5/§4.6) の妥当性。特に (1) 当初仮説 (壁ノードからの再構成が異常) の棄却は正しいか、(2) 残る実質差 = リミッタ評価点の不一致 を W1 として最優先にする判断、(3) SU2 の汎用セーフガード 3 件 (非物理再構成のエッジ単位1次化 / 更新率上限 / 点ロールバック) を forge に入れるべきか・既定 on にすべきか、(4) ベース発散を SERN 側の未解像問題に戻す判断

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
W1 の評価点整合は必要ですが、§4.5–4.6 の判断根拠に使った再構成診断が、ソルバと異なる座標を使っています。  
「セーフガードは発火しない」「ベース発散はソルバ側の課題ではない」という結論は撤回が必要です。

ローカル SU2 が `v8.5.0` であることを確認し、ソース照合、既存 HDF5 の再計算、`check_convergence.py` の再実行を行いました。新しい CFD 計算は実行していません。

1. **Critical — 測定器の座標が間違っており、仮説棄却と原因除外の根拠が崩れている。**

   根拠: [`check_wall_cv_drain.py:54`](/home/sano/work/forge/solver_density_cuda/tools/check_wall_cv_drain.py:54) は入力 HDF5 の `CELLS/centCoords` を読み、70 行で再構成距離に使います。しかし実行時は [`main.cpp:1095`](/home/sano/work/forge/solver_density_cuda/main.cpp:1095) と [`mesh.cpp:353`](/home/sano/work/forge/solver_density_cuda/mesh/mesh.cpp:353) により、これを**ノード座標へ置換**します。`run_0208` のログにも置換と最大移動量 `0.03593721241 m` が記録されています。

   `case/46.sern_design/run_0208_r4e_diag_2nd/` の保存された原始量・勾配・リミッタを、実行時のノード座標で再評価すると次の結果です。

   | 診断量 | plan の値 | 再評価 |
   |---|---:|---:|
   | `res_1.h5`、62100→62107 の `Ux_L` | 457.42 m/s | **651.765 m/s** |
   | 同、62104→62111 の `Ux_L` | 365.77 m/s | **521.171 m/s** |
   | `Ux_L / (0.5 Ux_neighbor)` | 0.70 | **約 0.999996** |
   | 中点とのずれが辺長の 1% 超の内部面 | 58899、47.7% | **48328、39.1%** |

   負圧の診断も変わります。正しい座標では、`run_0208/res_1.h5` の再構成圧力に非正値はなく、最初に検出した保存場は `res_6.h5`、47702→47703 の **−66.78 Pa** でした。`res_27.h5` では**ベース自身**の62104→62103 が **−15.85 Pa**、62102→62103 が **−38.79 Pa**。一方、対照の `run_0209_r4e_diag_nobase/res_1.h5`〜`res_40.h5` では非正値を検出しませんでした。

   したがって「ベース無しでも同じ負圧が毎 step 出るから原因ではない」「ベースでは物理性検査が発火しない」は支持できません。SU2 の検査は左右いずれかの負圧で発火します（[`CEulerSolver.cpp:1947`](/home/sano/work/forge/.external/su2-src/SU2_CFD/src/solvers/CEulerSolver.cpp:1947)）。

   **対案:** 最初に測定器を修正し、同一演算時点の状態で GPU の再構成値・実質量流束と照合する。その後に W1 の A/B を行う。「壁面外の中点に有限速度があるだけでは異常といえない」は妥当ですが、今回の数値から再構成全体を無罪にはできません。

2. **Major — 「1 次では壁 CV が質量を失えない」という説明は SLAU の式に反する。**

   根拠: [`convectiveFlux_slau_d.inc.cuh:493`](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:493) の質量流束には、両側速度に加えて `−chi*(P_R−P_L)/c_diss` が入ります。左右速度がともにゼロでも、

   \[
   \dot m=-\frac{S}{2c_{\mathrm{diss}}}(P_R-P_L)
   \]

   となり、壁ノード側の圧力が高ければ内部面へ質量が流出します。壁半割面のゼロ流束と、壁 CV の内部面流束は別です。

   **対案:** §4.1・§4.6 と測定器の「一方弁」説明を削除する。1 次の短時間完走は観測事実として残し、原因説明は実際の SLAU 面流束の符号付き総和、密度残差、陰的補正量で行う。

3. **Major — SU2 との「唯一の実質差」という整理は、対象の等温壁に対して成立しない。**

   根拠: plan が引用する `CNSSolver.cpp:508` 付近は熱流束壁の処理です。対象 run は等温壁であり、SU2 の [`BC_Isothermal_Wall_Generic:724`](/home/sano/work/forge/.external/su2-src/SU2_CFD/src/solvers/CNSSolver.cpp:724) は温度差から熱流束を作り、エネルギー残差と Jacobian に加えます。行削除は運動量だけです。

   forge は [`nodeWallDirichlet_d.cu:78`](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:78) で `T/roe/P` をピンし、[`同:163`](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:163) で等温壁のエネルギー残差をゼロ化します。**熱的閉包まで同型ではありません。**

   また、SU2 の Venkatakrishnan の平滑化係数は参照長さと設定係数に依存します（[`CLimiterDetails.hpp:185`](/home/sano/work/forge/.external/su2-src/SU2_CFD/include/limiters/CLimiterDetails.hpp:185)）。forge の局所 `volume` 依存とは異なり、勾配も未比較です。「同条件なら SU2 のリミッタは1」「SU2 の値より小さい」はソースの形式だけから導けません。

   **対案:** 「速度 Dirichlet・連続式を解く・内部辺を再構成する」という構造上の一致に限定する。等温壁、勾配、リミッタ、熱物性、陰的更新の差を比較表に追加し、同一メッシュ・対応する BC で実行するまでは SU2 の成否を予測しない。

4. **Major — W2 は既存ガードと重複し、連成状態を安全に更新する設計が欠けている。**

   根拠: [`update_d.cu:192`](/home/sano/work/forge/solver_density_cuda/cuda_forge/update_d.cu:192) に `updateGuardScale` が存在し、定常・dual-time の両方で使われています。[既存 accepted plan §5](/home/sano/work/forge/plans/accepted/time_integration-update-positivity-guard.md:58) は、このガードが試験した CFL 上限を改善しなかった結果を記録しています。SU2 の20%制限と同一ではありませんが、「更新の汎用セーフガードが無い」は誤りです。

   さらに forge は流れ更新後に SST、化学種、凝縮を別々に更新します（[`main.cpp:1581`](/home/sano/work/forge/solver_density_cuda/main.cpp:1581)）。流れ5変数だけを制限・巻き戻すと、`ΣρY=ρ`、組成依存 EOS、乱流エネルギー、周期共有 DOF との整合を失う可能性があります。TP のエネルギー基準値にも依存するため、SU2 の `|ΔρE|/|ρE|` をそのまま移植する設計では不足です。

   **対案:** W2を分割する。再構成の非物理状態に対する**面単位の左右同時フォールバック**は導入候補とし、組成・エンタルピーも整合して再計算する。更新率制限と点ロールバックは既存ガードとの差分、全連成状態の受理・復旧範囲、EOS 床より前の検査位置を別 plan で確定する。**現段階で3件を一括して既定 on にすることは推奨しません。**

5. **Major — W1 の実装範囲では周期経路と3次再構成の不整合が残る。**

   根拠: 通常経路の実体は [`limiter_r1_fused5_d:273`](/home/sano/work/forge/solver_density_cuda/cuda_forge/limiter_d.cu:273) ですが、周期 node は [`limiter_d.cu:329`](/home/sano/work/forge/solver_density_cuda/cuda_forge/limiter_d.cu:329) で別経路へ分岐し、[`limiterPeriodic_d.cuh:87`](/home/sano/work/forge/solver_density_cuda/cuda_forge/limiterPeriodic_d.cuh:87) も面重心を使います。

   また3次 MUSCL の増分は、勾配射影だけでなく隣接値差を含みます（[`convectiveFlux_common_d.cuh:99`](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_common_d.cuh:99)）。座標だけ揃えても、リミッタが評価する増分と実際の増分は一致しません。SU2 はリミッタ側でも `umusclProjection` を適用しています（[`computeLimiters_impl.hpp:192`](/home/sano/work/forge/.external/su2-src/SU2_CFD/include/limiters/computeLimiters_impl.hpp:192)）。

   **対案:** 通常・周期・次数ごとの再構成増分を列挙し、同じ増分を評価する共通処理を設計する。float32 では `Qt−Qc` の桁落ちも確認する。cell は対象外の分岐を維持し、現行ルールどおり実行検証は node に限定する。

6. **Major — W1 の「近傍 min/max 違反ゼロ」は既定 Venkatakrishnan の保証ではない。**

   根拠: [`limiterFunctions_d.cuh:13`](/home/sano/work/forge/solver_density_cuda/cuda_forge/limiterFunctions_d.cuh:13) は正の `eps2=volume` を使います。例えば `δmin=0, δm=−1, volume=1` なら `ψ=1/3`。評価点が完全一致していても、近傍最小値から負方向への再構成が残ります。

   **対案:** Barth の厳密有界性試験と、Venkatakrishnan の平滑化を含む許容逸脱試験を分ける。密度・圧力の正値性は別ゲートにする。W1 単独に厳密有界性まで要求すると、幾何整合の修正に別のリミッタ変更が混入します。

7. **Major — 検証条件が実装前に固定されておらず、ベース問題の切り離しも早すぎる。**

   根拠: [対象 plan:182](/home/sano/work/forge/plans/active/convection-node-wall-reconstruction.md:182) は許容差を後決めにしています。また「ベース発散は対象外」としながら、198行の完了条件は「壁 CV の `ro` が保たれる」のままです。

   今回、各診断 run の既存 CSV 全区間で再実行した判定は次のとおりです。

   | run（`case/46.sern_design/` 配下） | `check_convergence.py` |
   |---|---|
   | `run_0207_r4e_diag_1st/` | `NOT CONVERGED (stalled/plateau …)` |
   | `run_0208_r4e_diag_2nd/` | `DIVERGED (NaN/Inf)` |
   | `run_0209_r4e_diag_nobase/` | `NOT CONVERGED (stalled/plateau …)` |

   所在は [case README の run 一覧](/home/sano/work/forge/case/46.sern_design/README.md:258)。これらは短時間の発散診断には使えますが、離散定常解の不存在や未解像主因説の証明にはなりません。

   **対案:** 実装前に、共通 IC、設定差分、メッシュ品質、判定区間、全残差の `PASS`、比較対象量の `check_quasisteady.py` 判定、誤差許容値を固定する。標準回帰は [`verification/README.md`](/home/sano/work/forge/procedures/verification/README.md:5) のケースを具体的に指定し、周期・軸対称・TP を変更範囲に応じて追加する。SERN 側への担当移管はよいものの、原因は**未確定**として残し、W1 の実行 A/B を局所格子試験より先に行う。

8. **Minor — タイトル・索引・影響範囲が旧方針のまま。**

   根拠: [`plans/README.md:55`](/home/sano/work/forge/plans/README.md:55) は壁速度再構成の制限問題として掲載し、対象 plan の193行には存在しない「W5」が残っています。更新ガードを含むのに `update_d.cu` も影響範囲にありません。

   **対案:** スコープ確定後、タイトル、索引、§5、§7、§8、関連する `methods/` を同期する。

**推奨は、診断を修正してから W1 に絞って進めることです。** 優先順は、①座標・演算時点を揃えた診断と§4.5–4.6の訂正、②通常・周期・次数別の W1 設計と定量ゲートの確定、③共通 IC での W1 A/B、④面フォールバックの独立検証です。更新率制限・点ロールバックの既定化と「ベースは未解像が主因」という結論は、その結果が出るまで保留してください。

ファイルは変更していません。指摘・推奨は **plan 未反映**です。

指摘数: Critical 1 / Major 6 / Minor 1
