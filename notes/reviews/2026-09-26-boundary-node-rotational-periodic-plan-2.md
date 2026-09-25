# codex レビュー: boundary-node-rotational-periodic (plan)

- **plan**: [`plans/active/boundary-node-rotational-periodic.md`](../../plans/active/boundary-node-rotational-periodic.md)
- **stage**: `plan`
- **date**: 2026-09-26
- **commit**: `3de26cba` (feature/sern-design)
- **codex**: effort `high`, 6.5 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M6/m2
- **focus**: 2 回目の plan 段。**実装前に §4.0 の修正方針 (LSQ seam 勾配の 2 倍の修正) を見てほしい**。G0 で並進周期の node 継ぎ目の LSQ 勾配が線形場で正確に 2 倍と確認した (§6.2、case/09.Taylor-Green/_g0_lsq_seam/)。現行は gradLSQ=2 の係数事前計算 (solver_density_cuda/cuda_forge/calcGradient_d.cu:577-920、擬似逆 M⁺ と spectral truncation) で各部分 CV が自分側の隣接だけで勾配を出し、periodicGradientGather (cuda_forge/periodicNode_d.cu:166-212) が和を取る。提案: 事前計算で正規方程式 M=ΣΔxΔx^T を seam 越しに合算してから擬似逆を取り、各側の隣接に合併 M⁺ から重みを配って毎 step の gather は現行の和のまま、回転では M_m → R_m^T M_m R_m、Δx を root 系に回す。(i) この修正で合併 CV の LSQ 勾配 (両側の隣接を 1 つの stencil とした LSQ) と厳密に一致するか、重み (距離重み等) や spectral truncation との整合、境界ノード (壁∩seam) の扱い。(ii) 代替 (gather を和でなく root の値で上書き、あるいは 0.5 倍) との比較と、どちらが正しいか。(iii) 勾配を使う全経路 (粘性・再構成・リミタ・SST 生成・化学種・受動種) で他に seam 2 重計上が無いか。(iv) 修正の検証 (G0 の線形場で 1±1e-5、2 次場での精度、case/39 と case/09 の回帰で何を見るか)。§4.1 以降の回転設計の改訂も見てよいが主眼は §4.0。
- **extra**: `case/09.Taylor-Green/_g0_lsq_seam/G0_translational.txt`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
§4.0 の「正規行列を合併してから擬似逆を取る」方針を支持します。ただし、距離重み・隣接点の重複・打ち切りの定義が不足しています。  
さらに、SST 勾配の合併漏れと、全周比較を丸め誤差内に収められない構造的問題があります。以下を実装前に修正してください。

目的の同定は妥当です。`plans/README.md` と `accepted/` を確認した限り、既存の完了済み周期修正は cell の partner 転送が対象で、本件は未解決です。以下の数値実験は読み取り専用の Python 計算であり、forge の検証 run とは区別しています。ファイルは変更していません。

1. **Major — 単純な \(M\) の合算では、全周メッシュと同じ LSQ にはならない**

   **根拠:** 現行は無重み LSQ ではなく、\(w_{ij}=1/|\Delta x_{ij}|^2\) を使います。[`calcGradient_d.cu:668`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:668)  
   また、全周の双対面は primal edge を一意化して生成しますが、分割された seam の接線方向エッジは両側に存在します。そのまま部分行列を足すと、その隣接点だけ重複して数えます。[`gmshReader.hpp:1841`](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/gmshReader.hpp:1841)

   したがって、提案は**部分 stencil を重複込みで連結した LSQ**には代数的に一致します。しかし、**全周の一意な隣接 stencil の LSQ**とは一般に異なります。線形場ではどちらも正解になるため、G0 だけでは検出できません。

   反例として、中心を原点、両側の隣接点を
   \(A=\{(1,0.3),(0,1),(0,-1)\}\)、
   \(B=\{(-1,0.7),(0,1),(0,-1)\}\)、
   \(\phi=x^2+xy+y^2\) として現行の距離重みで計算すると：

   | 演算 | 勾配 |
   |---|---|
   | root 側だけ採用 | \((1.390000,\,0)\) |
   | 両側勾配を平均 | \((0.300000,\,0)\) |
   | 部分 \(M,b\) を単純合算 | \((0.492580,\,0.192580)\) |
   | 重複を除いた全 stencil | \((0.512338,\,0.353896)\) |

   **対案:** §4.0 に、重複 incidence の配分係数 \(\alpha_{mj}\) を明記してください。同じ物理隣接点への係数の総和を 1 にし、

   \[
   d_{mj}^{r}=R_m^{\mathsf T}d_{mj},\qquad
   M_r=\sum_{m,j}\alpha_{mj}w_{mj}d_{mj}^{r}(d_{mj}^{r})^{\mathsf T},
   \]
   \[
   c_{mj}^{m}=R_mM_{r,\tau}^{+}R_m^{\mathsf T}
                 \alpha_{mj}w_{mj}d_{mj}
   \]

   を焼き込む方式を推奨します。これなら毎 step の部分勾配の和を維持できます。隣接点の同定には周期像も含め、`periodicRoot` が同じでも異なる方向の隣接点を潰さないことが必要です。

   spectral truncation は**合併した \(M_r\) に一度だけ**適用します。現行の閾値処理は [`calcGradient_d.cu:671`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:671)。部分行列ごとの打ち切りや擬似逆の合算は不可です。打ち切りがある場合の線形場の期待値は真勾配そのものではなく、保持部分空間への射影です。

   壁∩seam でも実在する内部隣接ノードだけを合併し、壁の疑似点は追加しません。これは現行仕様と整合します。[`methods/discretization.md:640`](/home/sano/work/forge-sern-design/methods/discretization.md:640)

   **root 上書き・一律 0.5 倍は不採用**です。上の反例に加え、4・8 member の角、非対称 stencil、部分的な rank 欠損を処理できません。

2. **Major — SST の `k/ω` 勾配は、合併した直後に上書きされている**

   **根拠:** `periodicGradientGather` は [`main.cpp:1459`](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1459) で実行されますが、後段の [`main.cpp:1481`](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1481) が `ransGradient` を呼び、勾配をゼロクリアして再生成します。[`ransTransport_d.cu:211`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransTransport_d.cu:211)  
   その後の gather はなく、`F1` はこの未合併勾配を読みます。[`ransSource_d.cu:344`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransSource_d.cu:344)

   これは LSQ の「2 倍」と別の既存欠陥です。また `ransGradient` は全 `nPlanes` を走査し、化学種経路にある周期半割面の除外もありません。

   **対案:** `k/ω` 専用の勾配合併を、`ransGradient` の直後、`ransBlendF1` の前へ移してください。周期半割面を除外し、GG の部分寄与を一度だけ回転合算します。既存の早い gather から `k/ω` を外し、NS・化学種まで再合算しない構成にします。

   経路の区別も plan に必要です。NS の速度・密度・圧力・温度は LSQ、化学種・受動種・凝縮モーメントの勾配は GG です。後者は合併体積で正規化してから和を取る現在の方式が整合しており、**0.5 倍してはいけません**。[`speciesTransport_d.cu:639`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/speciesTransport_d.cu:639)、[`speciesTransport_d.cu:1220`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/speciesTransport_d.cu:1220)

3. **Major — §4.4 は回転との可換性を満たしても、実際の面値の bound を保証しない**

   **根拠:** 円筒成分の勾配で求めた制限係数を、Cartesian 速度の再構成へ掛けています。[`plan:75`](/home/sano/work/forge-sern-design/plans/active/boundary-node-rotational-periodic.md:75)  
   円筒成分の勾配には基底の微分項があるため、この二つの再構成増分は異なります。既存コードも「リミタと流束が同じ増分を使う」ことを必須にしています。[`reconIncrement_d.cuh:2`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/reconIncrement_d.cuh:2)

   解析的な反例は \(u_r=1,u_\theta=0\)。円筒成分の勾配はゼロですが、半径 1、角度差 0.2 rad の辺の中点へ厳密な Cartesian 勾配で再構成すると、
   \(u_{r,f}=1.00492108,\ u_{\theta,f}=-0.00099501\)。
   円筒成分の一定値 bound を外れます。

   **対案:** 実際に流束へ渡す `recon_increment` の横速度ベクトルを、評価点ごとの共通基底へ射影し、同じ基底へ射影した近傍値から bound を作って共通 \(\psi_\perp\) を求めてください。診断も最終面値に対して行います。V6 の「Cartesian モードとの差がノイズ床以内」は削除し、回転共変性・最終面値の有界性・面値誤差の格子収束で判定すべきです。

4. **Major — seam の部分双対面を別々に評価する限り、V2-i の全周一致は一般に成立しない**

   **根拠:** 全周側は同じ edge の双対面ベクトルを合算し、そのノルムを面積にします。[`gmshReader.hpp:1903`](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/gmshReader.hpp:1903)、[`gmshReader.hpp:1938`](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/gmshReader.hpp:1938)  
   一方、seam 側は各部分面で流束を評価します。SLAU の質量流束には面積 `sss` と単位法線に依存する非線形項があります。[`convectiveFlux_slau_d.inc.cuh:579`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:579)

   したがって、部分法線が平行でなければ
   \[
   \sum_m F(Q_L,Q_R,S_m)\ne F(Q_L,Q_R,\sum_m S_m).
   \]
   これは勾配を直しても残ります。環状 quad を単位厚さに押し出した例で、径方向 \(r=1\ldots2\)、隣接セル角幅 0.2 rad とすると、部分面積和／合併面積は **1.005020918**。静止・非零圧力差の SLAU 圧力散逸項だけで約 **0.502%** 違い、丸め誤差ではありません。

   **対案:** 全周の作用素との一致を受入条件にするなら、seam 接線エッジの部分双対面も同定・合併し、全周と同じ法線・面積で流束評価する設計を追加してください。V2-i の前提検査には、対応 ID だけでなく、この面集約まで含める必要があります。

5. **Major — 全周／セクタでリミタの基準値が自動的に変わる**

   **根拠:** `limiterRefLength` の既定値は領域の境界箱対角です。全周とセクタでは異なります。[`main.cpp:1328`](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1328)  
   `limiterRoRef/PRef/ARef` の自動平均も、合併体積を持つ seam member を全員集計しています。[`main.cpp:1312`](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1312)  
   基準長は実際のリミタ式に入るため、同じ円筒モードにするだけでは同じ作用素になりません。[`limiter_d.cu:86`](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/limiter_d.cu:86)

   **対案:** V2・V4・V6 では四つの基準値を共通値として明示固定してください。自動算出を回転周期でも保証するなら、代表領域に依存しない基準長と、一意 DOF による平均へ修正する必要があります。

6. **Major — G0 と V5 だけでは、修正の正しさも回帰も判定できない**

   **根拠:** G0 は線形場の比だけ、V5 は「面流束ビット同一」が中心です。[`plan:119`](/home/sano/work/forge-sern-design/plans/active/boundary-node-rotational-periodic.md:119)、[`plan:128`](/home/sano/work/forge-sern-design/plans/active/boundary-node-rotational-periodic.md:128)  
   しかし LSQ 修正は粘性・再構成・SST 生産を変え、時間発展後には勾配を直接読まない流束まで変わります。また `case/39/run_0007` の設定には `wallTreatmentSST: 1` と旧 `kInf/omegaInf` が残っています。[`solverConfig.yaml:49`](/home/sano/work/forge-sern-design/case/39.periodic_hills/run_0007_coarse_rans/solverConfig.yaml:49)

   **対案:** 実装前に次を §6 へ追加してください。

   - **G0:** 2・4・8 member、非対称 stencil、壁∩seam、root 交換、60°回転を含める。非退化方向は最大誤差 \(10^{-5}\)、ゼロ成分は絶対誤差で判定。退化方向は同じ打ち切りを施した参照解と比較する。
   - **二次場:** CPU double の一意 stencil を参照とする作用素試験と、格子細分試験を分ける。一般の非対称 stencil では二次場の**勾配誤差は通常 \(O(h)\)** であり、面再構成の \(O(h^2)\) と混同しない。
   - **周期整合:** 任意の線形スカラー場は回転周期条件を満たさない。局所作用素試験として定義するか、周期条件に整合する場を選ぶ。速度テンソルは剛体回転などで検査する。
   - **case/39:** 現行設定へ整備し、壁∩seam の速度勾配、壁応力、`k/ω` 勾配、`F1` を確認。6000 step 固定を合格条件にせず、定常結果には全残差と対象量の VERDICT を要求する。
   - **case/09:** 一意 DOF で質量・運動量・全エネルギー、KE・エントロピー履歴を確認する。TGV に定常収束 PASS は要求しない。[検証手順:69](/home/sano/work/forge-sern-design/procedures/verification/09-taylor-green.md:69)
   - **閾値:** V2-i の勾配差を流束の次元で判定しない。量ごとの尺度・ノルムを定義する。V3 のソース込み収支にはソース積分も尺度に含め、流束ゼロでも判定可能にする。

   V4 の起動条件には、現行レシピに沿った段階起動・初期 `k/ω`・実効設定の区間記録も必要です。

7. **Minor — 改訂前後のスコープが併存している**

   **根拠:** §2 は円筒リミタを対象外、§4.4 は今回実装。§4.3 は陰解法 fold を別計画へ分離、§5 は実装対象。§7 は並進ビット不変を残しています。[`plan:28`](/home/sano/work/forge-sern-design/plans/active/boundary-node-rotational-periodic.md:28)、[`plan:95`](/home/sano/work/forge-sern-design/plans/active/boundary-node-rotational-periodic.md:95)、[`plan:150`](/home/sano/work/forge-sern-design/plans/active/boundary-node-rotational-periodic.md:150)  
   `methods/boundary.md` も旧リミタ方針です。

   **対案:** §2・§5・§7・現在仕様を同期してください。陰解法の行縮約を今回は分離する判断は支持しますが、V7 は同じ反復予算と判定区間で比較することを明記してください。

8. **Minor — G0 の保存成果物だけでは実測を再検証できない**

   **根拠:** この checkout の `_g0_lsq_seam/` には結果テキストと YAML 二つだけがあり、入力 HDF5、勾配出力、生成・集計スクリプトがありません。YAML の `initial` も線形場を定義していません。[`solverConfig.yaml:19`](/home/sano/work/forge-sern-design/case/09.Taylor-Green/_g0_lsq_seam/solverConfig.yaml:19)

   テキストの **seam 平均 2.000000／内部平均 1.000000** はコード上の欠陥と整合しますが、今回確認できたのはその記録までです。

   **対案:** 修正前後の再現用 `run_NNNN_*`、入力生成・集計スクリプト、更新前の勾配出力、使用 revision／バイナリを保存し、case README に登録してください。G0 は作用素試験として扱い、1 step の結果を収束結果とは呼ばないでください。

**推奨は、全周の一意 stencil を基準にした合併 LSQ を採用し、まず並進周期の修正を独立に通すことです。** 実装前の優先順は、①指摘1の式・重複処理と指摘6の試験を確定、②指摘2の SST 勾配順序を修正対象へ追加、③回転実装前に指摘4・3・5を設計へ反映、④文書と証拠を整備、です。レビューのため **plan 未反映**です。

指摘数: Critical 0 / Major 6 / Minor 2
