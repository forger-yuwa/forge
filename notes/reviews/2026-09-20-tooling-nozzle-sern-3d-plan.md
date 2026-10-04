# codex レビュー: tooling-nozzle-sern-3d (plan)

- **plan**: [`plans/active/tooling-nozzle-sern-3d.md`](../../plans/active/tooling-nozzle-sern-3d.md)
- **stage**: `plan`
- **date**: 2026-09-20
- **commit**: `a1d19222` (feature/sern-design)
- **codex**: effort `high`, 5.3 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M6/m1
- **focus**: §4.21 のやり直し設計だけを見てほしい。(a) O1 の『出力が解を変えるか』の切り分け設計は十分か (出力経路が状態を触りうる箇所をソースで指摘してほしい: outputH5_XDMF / dependentVariables / EOS 床の再適用など)、(b) O2 の指標 (排出率・成長率・残差・定常性) で機序を識別できるか、他に取るべき量はあるか、(c) O3 で粗/細の初期場を interp_field で揃える方針の落とし穴、(d) そもそも『同一入力で NaN が step 34/38/53 とばらつく』ことの説明として、指数成長による増幅と実際の非決定性のどちらを先に疑うべきか、判別方法は

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
§4.21 の再調査は必要ですが、O1 の二択と O3 の「移植すれば同じ初期場」という前提は修正が必要です。  
優先する説明は **実行ごとの微小な非決定性と、その増幅**です。出力による状態変更の併存は、別試験で排除してください。

1. **Major — O1 の判別論理は成立しない。非決定性は既存 run ですでに観測できる。**

   **根拠:** `case/46.sern_design/` 配下の `run_0292_wake_cfl5/` と `run_0323_repro_rep/` を独立に比較しました。設定は同一で、`sern.h5`・境界条件・種 DB・probe の SHA-256 も同一です。

   | 比較対象 | 実測 |
   |---|---|
   | `res_0.h5` | 全 `VALUE` 配列がビット同一 |
   | `res_1.h5` の `ro` | 18,811 点で相違、最大絶対差 `3.72529e-8 kg/m³` |
   | `res_1.h5` の `roe` | 13,613 点で相違、最大絶対差 `0.0625 J/m³` |

   差の発生源として、SLAU の残差集積には浮動小数点 `atomicAdd` が存在します。[`convectiveFlux_slau_d.inc.cuh:563`](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:563)。SST 勾配にも同じ構造があります。[`ransTransport_d.cu:55`](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransTransport_d.cu:55)。加算順序による run 間差は NVIDIA も説明しています。[NVIDIA の解説](https://developer.nvidia.com/blog/?p=113316)

   **指数成長は差を増幅する説明であり、同一初期状態から差が生じる説明ではありません。** また、出力の副作用と非決定性は併存できるため、「反復にも差があれば出力バグではない」とは判定できません。今回の実測だけで `atomicAdd` を唯一の原因と断定することもできません。

   **対案:** 最終場の二回比較を主試験にせず、初期化→残差組立→各 sweep→更新→EOS の順に**最初の差**を探してください。差が最初に出た集積箇所を固定順序で評価する診断経路で再試験し、その後、既知の微小摂動を与えて局所モードの増幅率を測る順序を推奨します。

2. **Major — O1 に、出力呼出しそのものの無作用性を検査する試験がない。**

   **根拠:** 現行ソースで確認できた経路は次のとおりです。

   | 箇所 | 実際の作用 |
   |---|---|
   | [`output.cpp:65`](/home/sano/work/forge/solver_density_cuda/output/output.cpp:65) | `outputH5_XDMF` の書出し処理は D2H コピー。ホストの `var.c` を更新するが、EOS は呼ばない |
   | [`output.cpp:278`](/home/sano/work/forge/solver_density_cuda/output/output.cpp:278) | 境界出力はホストの `bc.bvar` と出力用構造を更新する |
   | [`main.cpp:1374`](/home/sano/work/forge/solver_density_cuda/main.cpp:1374) | 通常の `dependentVariables` は残差組立側で呼ぶ |
   | [`dependentVariables_d.cu:76`](/home/sano/work/forge/solver_density_cuda/cuda_forge/dependentVariables_d.cu:76)、[`:209`](/home/sano/work/forge/solver_density_cuda/cuda_forge/dependentVariables_d.cu:209) | 密度床を適用し、TP 経路で保存エネルギー `roe` を再構成する。単なる読取りではない |

   `roe` の再構成は床に到達した点だけの処理でもありません。**診断のために EOS を追加呼出しすると、診断自体が介入になります。** 一方、現行の通常 HDF5 出力から EOS 再適用へ至る直接経路は確認できませんでした。

   また、名前が出力専用の `wallStressForOutput_node_d` は、出力間隔の分岐ではなく粘性流束 wrapper 内で実行されます。[`viscousFlux_d.cu:1031`](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:1031)。名前だけで O1 の原因候補にしてはいけません。

   **対案:** 時間前進を止め、同じメモリ状態に対して出力呼出し前後を比較してください。GPU の保存量・次反復で読む補助状態と、更新されたホスト配列の後続読者を対象にします。そのうえで「出力なし／同期と退避コピーのみ／通常出力」を分けます。原始量を揃えるための EOS 評価は**複製した状態上だけ**で行い、計算状態へ戻さない設計にしてください。

3. **Major — O2 の指標は症状を分類できるが、排出・線形反復・EOS 補正の機序を識別できない。**

   **根拠:** `dt_local` は CFL と局所状態に依存します。[`setDT_d.cu:216`](/home/sano/work/forge/solver_density_cuda/cuda_forge/setDT_d.cu:216)。したがって `%/step` は粗細・CFL 間で同じ尺度ではありません。陰解法では、さらに残差を直接 `dt_local` 倍した量が密度更新になるわけでもありません。[`main.cpp:1614`](/home/sano/work/forge/solver_density_cuda/main.cpp:1614)

   SLAU の実際の `massflux` には圧力差項が含まれます。[`convectiveFlux_slau_d.inc.cuh:526`](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:526)。密度低下だけから「再構成速度による排出」とは特定できません。符号反転と包絡の成長だけでも、DPLUR の線形反復不安定とは特定できません。

   **対案:** 現行指標に以下を追加してください。

   - 固定したベース節点・物理領域での質量 `ΣρV`、実際の面 `massflux` の流入・流出別収支、圧力差項の寄与。
   - `dt_local`、DPLUR の密度補正、EOS・床・壁ピンによる保存量変更を**別々の時相**で記録。
   - 固定した線形系に対する各 sweep の線形残差と補正量。これで線形解法と外側の非線形更新を分ける。
   - 全残差の正式判定に加え、ベース領域の残差二乗和寄与率、体積・基準量で規格化した残差。

   成長率は固定節点または固定空間モードで、床・クリップ発動前の窓を使って評価します。移動する `ρ_min` や符号反転回数だけを成長率の根拠にしないでください。

4. **Major — O2 の「すべて同じ短い step 数で完走」は、O1 の観測と長期受入の両方を満たさない。**

   **根拠:** `run_0324_repro_out50/` は 40 step 実行ですが、通常の場出力は `res_0.h5` だけです。出力条件は単純な剰余判定で、最終 step の強制保存はありません。[`output.cpp:264`](/home/sano/work/forge/solver_density_cuda/output/output.cpp:264)

   今回の再判定は次のとおりです。

   | run（`case/46.sern_design/` 配下） | ツールの VERDICT |
   |---|---|
   | `run_0292_wake_cfl5/` | `DIVERGED (NaN/Inf)` |
   | `run_0323_repro_rep/` | `DIVERGED (NaN/Inf)` |
   | `run_0324_repro_out50/` | `NOT CONVERGED (stalled/plateau)` |
   | `run_0325_repro_rep2/` | `DIVERGED (NaN/Inf)` |

   `run_0324` の準定常判定も **`TRANSIENT-UNSETTLED (only 1 res_*.h5)`** です。run 索引は [`case/46.sern_design/README.md:300`](/home/sano/work/forge/case/46.sern_design/README.md:300)。

   また、`check_quasisteady.py` の既定量は `shock/asym/machmax/pmax` で、ベース圧・再循環長・ノズル力ではありません。[`check_quasisteady.py:377`](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:377)

   **対案:** 試験を二段にします。短期診断は共通停止点と独立した終端保存を用意し、反復間ばらつき込みで比較。長期受入は別に実施し、全残差と**実際の対象量**の時系列を判定します。ベース圧・領域質量・力係数などは `--series-csv` 経由で判定可能です。途中で異常化したケースは失敗として残し、完走例だけを比較しないでください。

5. **Major — O3 の `interp_field.py` は、そのままでは「共通初期場」を保証しない。**

   **根拠:** [`interp_field.py:119`](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:119) は、`P/Ux` がある出力ファイルでは保存済みの `roUx/roK/roOmega` を使わず、`ro×Ux`、`ro×k`、`ro×omega` から作り直します。

   定常陰解法の出力前には、更新後の速度・乱流原始量を一括再評価していません。[`main.cpp:1769`](/home/sano/work/forge/solver_density_cuda/main.cpp:1769)。実際に `run_0292_wake_cfl5/res_10.h5` では、保存された `roUx` と `ro×Ux` の最大差が **`0.734817 kg/(m²·s)`** あり、該当点では **`9.85081 → 10.58563`** に変わります。これは補間前に生じる状態変更です。

   さらに、最近傍探索は **x,y のみ**です。[`interp_field.py:34`](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:34)。現在の平面 2D 診断では z の省略は問題になりませんが、3D へ流用するとスパン方向の異なる点を混同します。最近傍移植には保存積分・境界所属を保つ制約もありません。[同`:159`](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:159)

   **対案:** 保存量が存在する場合は全保存量を直接移すことを前提にしてください。共通 donor から両格子を作り、移植距離、壁／内部の誤対応、ベース近傍分布、質量・エネルギー積分差を記録します。さらに**初回 EOS・BC 適用後**の状態を検査し、各格子内の CFL 1/5 は同じ初期ファイルを使います。「同一物理場」ではなく「同一 donor と検証済み写像」と表現するのが正確です。

6. **Major — O4 はキー確認だけでは不十分。実効作用素とバイナリを固定する必要がある。**

   **根拠:** `implicitRelax` は `time.deltaT` から読みます。[`solverConfig.cpp:327`](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:327)。また SST 更新にも同じ値が掛かるため、変更は流れの block-DPLUR だけへの介入ではありません。[`update_d.cu:336`](/home/sano/work/forge/solver_density_cuda/cuda_forge/update_d.cu:336)

   今回の checkout には `main.cpp` の等温壁ピン順序・step 末ピンに未コミット変更もあります。現行ソースと過去 run の実行バイナリを同一視できません。加えて、scaled limiter を使う場合の基準値は初期場から自動決定されます。[`main.cpp:1269`](/home/sano/work/forge/solver_density_cuda/main.cpp:1269)

   **対案:** YAML を構造として編集し、読込後の実効設定を期待値と照合して、不一致なら投入を失敗させます。各 run にバイナリ SHA-256、対応ソースと差分、入力ファイル署名、GPU・実行に影響する環境変数を残してください。O4 はまず**連成系全体の緩和 A/B**として扱い、結果だけから流れ Jacobian の原因認定をしないこと。自動基準値が有効なら、それも全ケースで明示固定します。

7. **Minor — §4.21 が残作業の正本に反映されていない。**

   **根拠:** [`plan:715`](/home/sano/work/forge/plans/active/tooling-nozzle-sern-3d.md:715) 以降の §5.1 は依然として W2 を先頭に置き、R4e 行にも撤回済みの原因説明が残っています。run 索引にも「粗は CFL 非依存」「細分で排出が消える」が残っています。[`README.md:292`](/home/sano/work/forge/case/46.sern_design/README.md:292)

   **対案:** O1→測定器・初期場の検証→O3/O4→長期受入を §5.1 と §6 に反映し、旧結果は撤回表示を付けて履歴として残してください。

**推奨は、まず「状態を変えない測定と再現性の監査」を完了してから 2×2 を再実施することです。** 実装前の修正優先順は、①O1 の二択撤回と出力前後検査、②時相・保存量移植・実効設定の固定、③収支と線形残差を含む O2、④短期診断と長期受入の分離です。既存の出力 plan はこの監査を完了した根拠にはなりません。

ファイル変更なし。以上の修正案は **plan 未反映**です。

指摘数: Critical 0 / Major 6 / Minor 1
