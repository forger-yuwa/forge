# codex レビュー: tooling-energy-balance-diagnostics (plan)

- **plan**: [`plans/active/tooling-energy-balance-diagnostics.md`](../../plans/active/tooling-energy-balance-diagnostics.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `8e0523eb` (feature/sern-design)
- **codex**: effort `high`, 4.5 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M7/m1
- **focus**: 発注元 (case-hypersonic-gap-heating-validation §4.8) が求める『壁熱量を離散スキームの言葉で検算する』ために、この診断出力の設計 (§4) とスコープが妥当か。対象 CV・拘束前 res_roe・境界数値流束・符号/単位・解除試験 (1D 熱伝導) の定義で足りるか。node の Dirichlet 拘束 (nodeWallDirichlet_d.cu) と出力登録 (variables.hpp / output.cpp) の実装を読んで、実現可能性と既定挙動・性能への影響を評価してほしい
- **extra**: `plans/active/case-hypersonic-gap-heating-validation.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
目的は妥当で、離散収支診断そのものが解決済みの計画は確認できませんでした。  
しかし現稿は、保存する配列の候補までで、拘束による熱授受・評価時点・面の対応関係が未定義です。このまま実装すると、出力は増えても壁熱量の検算が成立しません。

1. **Major — 拘束前 `res_roe` は、そのまま壁熱量ではない**

   **根拠:** [`nodeWallDirichlet_d.cu:80`](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:80) は `roe=ρ(e(Tw,Y)+ek)` を上書きし、同ファイル:90 は残差をゼロ化します。一方、[`timeIntegration_d.cu:1023`](/home/sano/work/forge/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1023) はエネルギー補正をゼロにしても**密度の更新は残します**。したがって壁温一定でも、密度・組成の変化に伴う壁 CV の蓄積エネルギー変化があります。また、拘束前残差には既に壁面の伝導流束が入っています。

   **対案:** §4 に次の収支を定義してください。`F` は外向き、`S` と拘束反力 `C` は流体への供給を正とします。
   $$
   R_i^{raw}=-\sum_f F_{if}+S_i,\qquad
   D_t(V_iE_i)=-\sum_f F_{if}+S_i+C_i.
   $$
   定常の拘束行なら `C_i=-R_i^{raw}` ですが、**物理境界流束と拘束反力を組み合わせて**実効壁熱量を求める必要があります。過渡まで扱うなら、温度ピン・no-slip の状態上書きによる `Δ(VE)` も記録してください。壁 CV を含む領域と、第一内部列から始まる領域では、比較する熱量の定義を分ける必要があります。

2. **Major — 残差・流束・保存場の評価時点が一致しない**

   **根拠:** [`main.cpp:1583`](/home/sano/work/forge/solver_density_cuda/main.cpp:1583) で残差を組み、その後:1610–1615で DPLUR 補正を適用し、:1744で更新後の場を出力します。単に退避配列を追加すると、**更新前の流束・残差と更新後の `roe`** が同じ出力ファイルに入ります。RK も各ステージの評価後に状態を更新します。dual-time では:1799–1800で空間残差の後に BDF 項が追加されます。

   **対案:** 最初の対象を定常診断に限定し、既存の残差組立て中に、同じ状態のエネルギー・流束・拘束前残差を一組として保存してください。`step`、評価位相、RK stage／dual-time subiteration、空間残差か BDF 込みかを明記します。**定常 DPLUR の反復差分を物理的なエネルギー変化率に換算してはいけません。** 出力時の再評価を選ぶ場合も、状態ピンを含む `assembleResidual` を無条件に再実行すると診断が計算を変えます。

3. **Major — エネルギー流束・ソースの収集範囲が不足している**

   **根拠:** エネルギー残差の供給元は対流・粘性だけではありません。

   - [`speciesTransport_d.cu:271`](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:271): 化学種拡散のエンタルピー輸送 `Σh_sJ_s`。
   - [`viscousFlux_d.cu:285`](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:285): `sstEnergyIncludesK` の `k` 拡散。
   - [`ransSource_d.cu:230`](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransSource_d.cu:230): SST のエネルギーソース。
   - [`axisymmetricSource_d.cu:232`](/home/sano/work/forge/solver_density_cuda/cuda_forge/axisymmetricSource_d.cu:232): `axisymMethod: 1` のエネルギー幾何ソース。
   - [`bodyForce_d.cu:43`](/home/sano/work/forge/solver_density_cuda/cuda_forge/bodyForce_d.cu:43): 体積力の仕事。

   発注元には多成分 TP・SST が含まれるため、純伝導だけ閉じても要求を満たしません。

   **対案:** `res_roe` への全加算経路を一覧化し、面流束と体積ソースを別々に記録してください。流束は各カーネルが**実際に残差へ加算した値**を採取し、後処理で再計算しないこと。`roe` と `roe+roK` のどちらを収支対象にするかも属性に残します。初版で未対応の物理・スキームは、診断要求時に明示的に拒否してください。他方程式の診断を対象外にしても、そこからエネルギー式へ入る寄与は省略できません。

4. **Major — 面配列は既存の `extraFields` 登録だけでは出力できない**

   **根拠:** [`output.cpp:69`](/home/sano/work/forge/solver_density_cuda/output/output.cpp:69) は選択した量を `copyVariables_cell_D2H` に渡し、:144–158では `var.c` を `nCells` 個として出力します。面配列は別の `var.p`／`nPlanes` 管理です（[`variables.cpp:303`](/home/sano/work/forge/solver_density_cuda/variables.cpp:303)）。`eflux_face` を `output_cellValNames` に足すだけでは、面配列の参照・長さ・接続情報を扱えません。

   **対案:** CV 診断と面診断の HDF5 スキーマを分けてください。面側には最低限、`face_id`、owner／neighbor、境界種別・`physID`、向き、使用した面積ベクトル、流束を保存します。node の可視化用 primal 面と、残差を組む dual 面を混同しないこと。§5 は「登録」より先に**データモデルと採取位置の確定**を置き、影響範囲に流束カーネル、`main.cpp`、変数確保・転送を追加すべきです。

5. **Major — `extraFields` は計算・メモリ確保の opt-in を保証しない**

   **根拠:** [`output.cpp:39`](/home/sano/work/forge/solver_density_cuda/output/output.cpp:39) は `output.level>=2` なら登録量を全出力します。また [`variables.cpp:282`](/home/sano/work/forge/solver_density_cuda/variables.cpp:282) は登録された CV 配列を一律確保します。静的登録と毎回の退避では、診断を要求していない実行にもメモリ・帯域コストが発生します。

   **対案:** 明示要求を初期化時に解決し、有効時だけ配列登録・確保・採取・転送を行う契約にしてください。`level: 2` 単独で有効化するかも明記し、今回の「opt-in」なら有効化しない設計を推奨します。OFF の回帰だけでなく、**ON/OFF で解が変わらないこと**と、OFF の追加確保・カーネル・同期がないことを確認してください。

6. **Major — CV 選択・周期同一視・幾何単位の定義がない**

   **根拠:** [`mesh.cpp:699`](/home/sano/work/forge/solver_density_cuda/mesh/mesh.cpp:699) は周期合併前の部分体積を保存し、:710–717では合併体積を全 member に複写します。さらに [`main.cpp:1423`](/home/sano/work/forge/solver_density_cuda/main.cpp:1423) は壁拘束の**後**に周期残差を合算・broadcast します。出力された `volume` と各ノードの値を単純積分すると周期部分を重複計上します。軸対称 method 0 は面積に半径を掛ける方式です（[`variables.cpp:505`](/home/sano/work/forge/solver_density_cuda/variables.cpp:505)）。

   **対案:** 対象領域はまず「指定した solver CV の和集合」と定義し、所属が片側だけの dual 面を領域境界にしてください。任意の幾何断面で CV を切る機能は初版から外すのが妥当です。周期は root 単位の一回集計、または部分 CV 単位の集計に統一し、対応表を出力します。単位も、平面 2D の単位スパン当たり `W/m`、3D の `W`、軸対称 method 0 の単位角度当たり量と全周換算を区別してください。未対応なら周期・軸対称を診断対象として受理しないことです。

7. **Major — 解除試験が符号・精度・対象現象の点で成立していない**

   **根拠:** plan §4 は外向き正ですが、§6 は「流束和＝内部エネルギー変化率」としています。ソース・拘束なしなら正しくは**負の流束和**です。実装も対流では owner に負、粘性では owner に正を加えています（[`convectiveFlux_slau_d.inc.cuh:528`](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:528)、[`viscousFlux_d.cu:320`](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:320)）。float の `atomicAdd` で組んだ残差に「double 相当」を要求する根拠もありません。

   既存の純伝導例 `case/24.laminar_channel_bl/run_isoT_condN_node/` に判定器を再実行した結果は、**`NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)`**。最終代表値は `rms_roe=4.42e-4`、`rms_roUy=8.27e-7` で、ともに rising 判定でした。既存の線形場検証を、そのまま収束済みの解除試験として継承できません。

   **対案:** §6 を以下の独立した確認に分けてください。

   - **組立て恒等式:** 同一状態の、独立に保存した残差と面流束・ソースを比較。内部面相殺、壁列、第一内部列、領域和を検査する。
   - **壁の拘束収支:** 既知の純伝導問題で、物理境界流束・拘束反力・実効壁熱量を別々に検査する。全領域の正味ゼロだけでは、上下壁を両方誤ってゼロにしても通ります。
   - **実用途:** SLAU の対流・再構成・粘性仕事を含む node ケース、および発注元の低 Re SST／TP 経路を追加する。標準の `case/48` と整合させる。
   - **定量判定:** float の加算数と `Σ|F|` に基づく丸め誤差限界を事前定義し、解析解との差は別の離散化誤差として判定する。定常解の主張には `check_convergence.py`、熱量系列には `check_quasisteady.py` の VERDICT を要求する。

   新規 run のメッシュ品質、IC、段階起動、run 索引も検証手順に含めてください。一般の収支自動判定ツールを対象外にすることは妥当ですが、解除試験の自動検算は必要です。

8. **Minor — 「数値流束の取得経路がない」という現状説明は広すぎる**

   **根拠:** [`viscousFlux_d.cu:560`](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:560) は、実際に残差へ入れた壁伝導寄与を `qwall=heatflux/sss` として保存し、[`output.cpp:355`](/home/sano/work/forge/solver_density_cuda/output/output.cpp:355) は境界変数を出力します。ただし、これは拘束後の実効熱量ではなく、内部 dual 面の流束も含みません。

   **対案:** §1 を「壁面の局所伝導寄与は取得可能だが、拘束反力と領域境界の全数値流束を対応づけられない」に修正してください。`output-level-and-h0`、`output-node-wall-surface-viz` は再利用する基盤であり、本計画を不要にする重複ではありません。なお調べた旧 `run_isoT_condN_node` の最終壁出力は上下とも `qwall` 全17点がゼロでした。フィールドの存在だけで既存 run を診断可能と扱わないことです。

**推奨は、定常 node の「同一状態における離散収支」に初版を絞り、§4 を書き直してから再レビューすることです。** 優先順は、①収支式と拘束反力、②採取位相、③全加算経路と対応機能、④CV／面スキーマと opt-in、⑤解除試験です。`cell` 機能を対象外にする判断は現行検証方針と整合しています。

ファイル変更・新規計算は行っていません。上記提案は **plan 未反映**です。

指摘数: Critical 0 / Major 7 / Minor 1
