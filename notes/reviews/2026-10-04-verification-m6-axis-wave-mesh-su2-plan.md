# codex レビュー: verification-m6-axis-wave-mesh-su2 (plan)

- **plan**: [`plans/active/verification-m6-axis-wave-mesh-su2.md`](../../plans/active/verification-m6-axis-wave-mesh-su2.md)
- **stage**: `plan`
- **date**: 2026-10-04
- **commit**: `bc1c84d1` (feature/gap-heating-precision)
- **codex**: effort `high`, 5.2 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M7/m1

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
軸細分と SU2 比較は有用ですが、現状の判定基準では「物理」を誤認できます。  
実測で、基準 run の未収束、TP Mach 算出の誤り、CPG 化で形状まで変わる実装経路を確認しました。実装前に修正が必要です。

目的は妥当です。親 plan は既に「軸ノードの数値スパイク」と「壁の δ 差への集束応答」を別々に記載していますが、後者の因果検証は未完了です。本計画は重複ではなく、その未検証部分を扱います。[親 plan:227](/home/sano/work/forge/plans/accepted/tooling-nozzle-deltastar-core-matched-euler.md:227)  
なお、過去の cell 近軸固着は「float32 一般」ではなく粘性対角の幾何不整合へ原因が訂正されています。今回それを既知の真因として再利用してはいけません。[既存の原因訂正:20](/home/sano/work/forge/plans/accepted/architecture-axisym-axis-singularity.md:20)

1. **Major — §6 の観測から「物理／数値」を決める論理が成立していません。**

   **根拠:** 「半径方向細分で変わらない」「細分で山が増える」「同一粗格子で二つのソルバが似る」は、物理の証明になりません。固定した x 解像度の誤差、両者に共通する格子誤差、反復誤差が残ります。特に「増えれば物理の集束」は根拠のない一方向の解釈です。[対象 plan:105](/home/sano/work/forge/plans/active/verification-m6-axis-wave-mesh-su2.md:105)

   閾値も排他的ではありません。例えば b₀=0.236、b₁=0.180、b₂=0.110、ε=0.05 %pt なら、緩和後の許容差 0.15 %pt により **H_num と H_phys が同時成立**します。また、B は TP から CPG への変更に加えて SST 補正も変えるため、A と同じ仮説の独立検証ではありません。

   **対案:** 今回の結論を「半径方向格子感度あり／検出できず」「CPG 条件でソルバ間差あり／検出できず」に限定してください。「物理」の確定は、x 方向細分・反復誤差の切り分けを経た後、親 plan が原因とする壁差を制御した比較へ持ち越すべきです。ε は絶対差で定義し、検出限界が大きければ許容差を広げず判定不能にします。

2. **Major — 基準 run が未収束で、指定した収束ゲートは実行すると拒否されます。**

   **根拠:** 再実行した判定は次のとおりです。

   | run（`case/45.isobutane_m6_d155/` 配下） | `check_convergence.py` の結果 |
   |---|---|
   | `run_0037_euler_rt77p02/` | `PASS (converged)` |
   | `run_0038_ns_final_rt77p02/` | `NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)` |
   | `run_0013_cpg_ns_plainsst/` | `NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)` |

   `run_0038` の最終残差は `rms_roUy=3.50e-4`、`rms_roe=6.28e-1` で、ともに plateau 判定です。計画どおりの `--from-floor run_0038` は **`REFUSED`、exit 2** になります。参照自身の通常判定 PASS が必須だからです。[判定コード:353](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:353)

   **対案:** §5.1 の「非 DIVERGED」を合格条件から除外してください。まず比較基準の収束を確立し、その後に細分比較へ進めます。未収束のまま実施する場合は探索的診断に限定し、「物理」「一致」の根拠に使わないことを明記してください。単に 12000 step を追加すれば解決するとは限りません。

3. **Major — 「同じ設定」だけでは旧 run と同じ離散化を再現できません。**

   **根拠:** `run_0038` は `convMethod: 1, limiter: 2` だけを指定しています。[保存設定:17](/home/sano/work/forge/case/45.isobutane_m6_d155/run_0038_ns_final_rt77p02/solverConfig.yaml:17)  
   現行コードでは省略された `limiterScaled` が **1** になり、評価点と無次元化を含む修正版へ切り替わります。既定変更日は元の最終形作成より後の 2026-09-20 です。[設定コード:619](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:619)

   したがって現在のバイナリで旧 YAML を使った A0′ は、再開だけでなく実装変更の影響も拾い得ます。それを ε として吸収すると、格子感度を検出できなくなります。

   **対案:** 全腕のソルバ・変換器の commit／バイナリ hash、精度、実効設定を固定してください。現行版を使うなら、同じ現行版で A0 を再確立してから A1/A2 と比較し、旧 `run_0038` との差は別の回帰比較として扱ってください。

4. **Major — 予備観測の TP Mach 算出が誤っています。**

   **根拠:** `axis_wave_compare.py` は `Mach` と `gamma` が出力されていない場合、TP でも定数 `G=1.27354` を使います。引数 `cpg` は判定に使われていません。[解析コード:12](/home/sano/work/forge/case/45.isobutane_m6_d155/axis_wave_compare.py:12)

   実際の `run_0037`／`run_0038` には両フィールドがなく、`sonic` は存在します。最終スナップショットを読み直すと、計画の `b(0)=0.23601 %pt` は再現できますが、`hypot(Ux,Uy)/sonic` では **0.23247 %pt** です。`run_0038` の全節点で両 Mach 定義の最大絶対差は **0.24555** でした。

   正しい音速による b の step 0/4000/8000/12000 系列を、`check_quasisteady.py` の `classify` に渡した結果は `STEADY`、末尾 drift 約 0.6% です。ただし、指摘2の残差判定は未収束のままであり、これは算出誤りの検証値です。

   **対案:** forge の Mach は `VALUE/sonic` から統一して算出し、欠落時は明示的に停止してください。§3 の数値を再計算し、事前登録基準も修正後の定義で固定してください。

5. **Major — 共通 Euler 参照は、現在の山指標では相殺しません。**

   **根拠:** §4.1 は参照の軸解像度差が相殺するとしていますが、`b=max(ΔM)−median(ΔM)` は非線形です。A1/A2 で最大位置や中央値を与える位置が変われば、**b(A2)−b(A0) から Euler の x 依存は消えません**。[対象 plan:61](/home/sano/work/forge/plans/active/verification-m6-axis-wave-mesh-su2.md:61)、[指標実装:55](/home/sano/work/forge/case/45.isobutane_m6_d155/axis_wave_compare.py:55)

   **対案:** 格子感度の主指標を、共通座標での NS 同士の直接差と、その区間最大ノルムにしてください。山の変化は、固定した評価位置・重みを持つ線形な指標でも確認します。Euler 比の b は副指標として残せますが、「参照誤差は相殺」とは書けません。

6. **Major — 腕 B の準備手順は、同一形状・同一条件を保証しません。**

   **根拠:** §5 の「CPG problem で prepare」は危険です。`prepare_ns` は `design_chain(p)` を呼び、`cfd_gas: cpg` なら **設計用ガスまで CPG に変更して MOC 壁を再生成**します。[設計コード:215](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:215)、[準備コード:699](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:699)

   さらに、流用予定の SU2 設定は出口圧 **1389 Pa**、元の `run_0038` は **2237 Pa** です。入口乱流も SU2 の強度・粘性比指定と forge の `k=1, omega=18000` は同一性が未確認です。[SU2 設定:26](/home/sano/work/forge/case/45.isobutane_m6_d155/run_0012_su2_sst/sst.cfg:26)、[forge 境界設定:1](/home/sano/work/forge/case/45.isobutane_m6_d155/run_0038_ns_final_rt77p02/bcondConfig.yaml:1)

   **対案:** B は `run_0038` の既存メッシュを直接複製・変換し、設計チェーンを再実行しない手順にしてください。座標・接続・境界タグの対応を検査し、両ソルバの EOS、輸送係数、入口乱流、出口圧を対照表で固定します。TP→CPG の初期場は CPG EOS でエネルギーを再構成する手順も必要です。

7. **Major — 軸の山を準定常判定する実装・精度条件が欠けています。**

   **根拠:** `check_quasisteady.py` の標準量は `shock, asym, machmax, pmax` で、軸 M、b、山の位置、B のトレンド除去残差はありません。未知の量はエラーになります。[抽出器:377](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:377)  
   既定許容差も drift **5%**、fluctuation **10%** です。Mach 自体へ適用すれば、本件の 0.1～0.2% 級構造を判別するには粗すぎます。[既定値:542](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:542)

   **対案:** 指標スクリプトの時系列対応を、run 投入より先に実装してください。b、固定位置の M、ピーク位置、B の比較量を保存し、`--series-csv` で判定します。時間変動の許容幅は空間比較の閾値より十分小さく事前指定し、SU2 側も出口 massflux だけでなく同じ局所量を検査してください。出力間隔と最低スナップショット数も §6 に必要です。

8. **Minor — 「軸側だけ細分」と品質 PASS の条件を具体化してください。**

   **根拠:** 現在見えている `_radial_fracs_capped` を実行すると、A1 は **nj=110、η≈0.691 まで一様化**、A2 は **nj=141、η≈0.835 まで一様化**します。変更領域は山がある η≲0.15 より大幅に広く、コア内の波の伝播誤差も変えます。[メッシュコード:87](/home/sano/work/forge/design/forge_design/meshing/mesh2d.py:87)  
   また、A0 メッシュの再検査結果は `VERDICT: SOFT-PASS`、最大 AR **1207.8**、最大 skew **0.439** で、計画記載の厳密な PASS ではありません。

   **対案:** 「壁近傍を保存したコア細分」と範囲を正確に記載してください。残す壁側節点の座標、変更範囲、導出 nj を検証項目に追加し、AR 緩和を使う場合は対象セルが適用条件を満たすことを確認して、全腕共通の品質基準として事前登録してください。

**推奨は、計画を「数値感度の診断」に修正して進めることです。** 実装前の優先順は、①結論と排他的な判定基準の修正、②バージョン固定と基準 run の収束確立、③正しい Mach・直接差・時系列ゲートの整備、④B の形状固定と物理条件照合、です。「壁差への物理的集束」の確定は、この診断を通過した後の因果検証に分けてください。

ファイルは変更していません。上記提案は **plan 未反映**です。

指摘数: Critical 0 / Major 7 / Minor 1
