# codex レビュー: tooling-sern-mesh-blocking (plan)

- **plan**: [`plans/active/tooling-sern-mesh-blocking.md`](../../plans/active/tooling-sern-mesh-blocking.md)
- **stage**: `plan`
- **date**: 2026-09-22
- **commit**: `2153b6d3` (feature/sern-design)
- **codex**: effort `high`, 5.7 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m1
- **focus**: §4.13 (L446–502) の診断と判断だけを見てほしい。問い: (1) 『鈍頭後端面の壁 CV の排出 (nodeWallDirichlet=1 の壁節点は u=0・T_w 固定で密度だけ自由、内点の外向き流れで質量が抜け、圧力勾配では戻れない)』という真因の推論は、示した証拠 (200 step 刻みの順序表、排出率 1.5 %/step の見積もりと実測の一致、nodeWallDirichlet=0 の A/B で完走、0→1 の段階起動が 55 step で NaN) で成立しているか。別の説明 (局所時間刻みの歪み、等温壁 T_w 1000 K と静止 IC の組合せ、後端面の x 第一間隔 0.25 mm と隅の 16 µm 層の AR、EOS 床の洗浄) を排除できているか。(2) 選択肢 (a)〜(e) の採否と『後端テーパを推奨 (L_t = 5t, t_base = t/5)』は妥当か。テーパでも残るベースに同じ機序があるなら、t_base を決める A/B の設計 (何を測るか) に不足はないか。他に安価な案 (壁 CV の連続式だけ弱形式にする等のソルバ側の小変更、後端面だけ slip、丸い後縁) を見落としていないか。(3) IC (run_junction_model.py、添付) の 2 つの修正 (端面からの距離での速度立ち上げ、境目の 5 mm 混合) は妥当か。診断の run はすべて AWS 上で、ここからは h5 を読めない。ファイル:行 で根拠を示し、推奨は 1 つに絞ってほしい。
- **extra**: `case/46.sern_design/cad/run_junction_model.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
§4.13 の「壁 CV の排出」は有力な仮説ですが、「Dirichlet 固定が真因と確定」は成立していません。IC には実際に不連続が残っています。  
推奨は **鈍頭形状を維持し、IC を修正したうえで局所収支を切り分けること**です。現時点でテーパを実装する根拠は不足しています。

対象は §4.13 と関連する B1d に限定しました。AWS の HDF5・実効設定・VERDICT は未確認であり、run の数値は提示記録として評価しています。ファイルは変更していません。

1. **Major — 「圧力勾配では戻れない」という説明は、実装されている質量流束と一致しません。**

   壁の運動量・温度を固定し、連続式を残す構造は確認できました。ただし、**壁速度を固定しても、内点側の応答と数値質量流束の圧力差項は残ります**。

   SLAU は `mdot` に `−χ(P_R−P_L)/ĉ` を含みます。また、`χ` を決める Mach 数は面法線速度だけでなく、**接線速度を含む速度の大きさ**から計算されます。[convectiveFlux_slau_d.inc.cuh:533](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:533)

   提示された内点速度 `(850, 0, 1290) m/s` と壁速度ゼロでは、`ĉ ≤ 約1092 m/s` なら `χ=0` になります。したがって、まず検証すべき具体的な候補は **「高速の接線流によって SLAU の圧力差項が消失・減衰している」**です。これは現時点では仮説ですが、「node の固定壁だから戻れない」より実装に即しています。

   排出率の見積もりも確証にはなりません。`χ=0`、壁側法線速度ゼロ、内点側法線速度 `u_n>0` とすると、実装式は

   \[
   \frac{\dot m}{A}
   =\frac{\rho_w\rho_i}{\rho_w+\rho_i}u_n
   \]

   となります。壁 CV の質量で割った排出速度は、単純な `u_n/2=425 m/s` ではありません。さらに block-DPLUR の更新には対角行列・近傍補正・緩和が入ります。[timeIntegration_d.cu:1101](/home/sano/work/forge/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1101)

   提示式を再計算すると排出率は **1.564 %/step**。一定率の陽的減少なら `2000 → 約85 Pa / 200 step` で、提示された `164 Pa` に対応する率は **約1.243 %/step**です。桁は近くても「一致して機序を確定」とはいえません。

   **対案:** 問題 CV と全隣接面について、実際の `massflux`、移流項と圧力差項、`χ`、両側音速、面積ベクトル、体積、`Δρ` を毎 step 記録してください。200 step 刻みの最小値表は発生域と進行を示しますが、面別収支と更新順序を証明しません。

2. **Major — `nodeWallDirichlet=0` の A/B は、運動量固定だけを切り分けていません。**

   このスイッチは、速度固定に加えて **温度ピン、エネルギー残差ゼロ化、陰解法のエネルギー行固定**も外します。[nodeWallDirichlet_d.cu:147](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:147)、[timeIntegration_d.cu:1511](/home/sano/work/forge/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1511)

   したがって `run_0428` の完走が支持するのは、「壁の強制処理一式を変えると発散経路が変わる」までです。**運動量固定、熱的固定、音速・局所刻みの変化のどれが支配的かは未分離**です。別の壁で床に達しており、成立した基準解にもなっていません。

   `run_0429` についても、plan 自身が「再循環が立った場ではなかった」と記しています。その場から戻して55 stepで失敗しても、「再循環を成立させてから戻す」という案 (e) 全体の反証にはなりません。[plan:493](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:493)

   **対案:** 「機序確定」「段階起動は不成立」を撤回し、運動量拘束と熱的拘束を分離した診断を先に置いてください。既存の `boundary-weak-isothermal-wall` が関連しますが、実装済みの `nodeIsothermalEnergyBC=1` は **CPG 限定で TP を拒否**するため、この run にそのまま指定できません。[weakIsothermalWall_d.cu:31](/home/sano/work/forge/solver_density_cuda/cuda_forge/weakIsothermalWall_d.cu:31)

3. **Major — IC の二つの修正は方向として妥当ですが、現在のコードは滑らかな IC を生成していません。**

   **速度:** `dist_face()` は `x<xf` で距離を突然 `∞` にします。しかし、端面の矩形の横には上流側にも流体があります。その流体で `x=xf` を横切ると、速度係数が不連続になります。[run_junction_model.py:90](/home/sano/work/forge/case/46.sern_design/cad/run_junction_model.py:90)

   ファイルから式を抽出して評価したところ、側壁内側の流体点 `y=0.05 m, z=ZW−16 µm` で、`x=LSW±1 nm` の速度係数は **`1 → 3.071×10⁻⁷`** でした。これは実メッシュの隣接点差の測定ではありませんが、**関数そのものの不連続を示す反例**です。「不連続は作らない」「速度不連続は原因でなかった」という結論は維持できません。

   **排気／外気混合:** `wz` と `wy` も、`x>LSW−dl`／`x>LCOWL−dl` で二値マスクから混合値へ突然切り替わります。[run_junction_model.py:80](/home/sano/work/forge/case/46.sern_design/cad/run_junction_model.py:80) 同様の評価で、流体内に `wz: 1→0.84555`、`wy: 1→0.65030` の段差が出ました。また、混合域の全幅はコード上 **`2dl=10 mm`** です。

   一方、保存量と `ρY_s` の共通重みによる混合、および速度低減分の運動エネルギーを `roe` から引く処理は、現在の `Uy=Uz=0` の IC では整合しています。[runner_sern.py:314](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:314)、[run_junction_model.py:97](/home/sano/work/forge/case/46.sern_design/cad/run_junction_model.py:97)

   **対案:** 端面距離の上流側切断をなくすか、流体領域で連続な軸方向の重みを掛けてください。混合重みも x 方向に連続化し、生成後は**実際の接続辺ごとの速度・圧力・温度・組成ジャンプと、壁ピン適用前後の差**を検査する必要があります。テーパ導入時には、固定値の `TSW/TC` による領域判定もメッシュ形状と共通化してください。

4. **Major — EOS 床、局所時間刻み、格子の影響を排除できていません。床については説明に誤りがあります。**

   TP の EOS は `ro_temp=max(ro,roMin)` を計算し、**`ro[ic]=ro_temp` と保存量へ書き戻します**。[dependentVariables_d.cu:76](/home/sano/work/forge/solver_density_cuda/cuda_forge/dependentVariables_d.cu:76)、[dependentVariables_d.cu:203](/home/sano/work/forge/solver_density_cuda/cuda_forge/dependentVariables_d.cu:203)

   一方、定常陰解法は最後の更新後に温度ピンを掛け、その後に出力します。温度ピンの圧力は床なしの `P=ρRTw` です。[main.cpp:1793](/home/sano/work/forge/solver_density_cuda/main.cpp:1793)、[nodeWallDirichlet_d.cu:80](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:80)

   つまり、床未満の出力は **「床が作用していない」証拠ではありません**。前段で床が密度を補い、更新で再び減り、壁温ピンが圧力床を上書きする経路と両立します。床による補正と排出が繰り返されている可能性が残ります。

   局所刻みも「受け側の音速だけ」では決まりません。面の速度、`sonic[ic0]`、両 CV の `V/A`、粘性項を用い、CV ごとに面 CFL の最大値を取っています。[setDT_d.cu:64](/home/sano/work/forge/solver_density_cuda/cuda_forge/setDT_d.cu:64)、[setDT_d.cu:181](/home/sano/work/forge/solver_density_cuda/cuda_forge/setDT_d.cu:181) 「過大な流入で負密度」という説明も、正味の流入だけなら密度は増えるため、流束の符号と陰解法補正の確認が必要です。

   格子についても AR 合格は局所更新の健全性を保証しません。双対閉性は **max `1.7e-4`、6,410 CV が規定値超過**で、ドライバは閉性検査の失敗を無視して続行します。[plan:440](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:440)、[run_junction_model.py:68](/home/sano/work/forge/case/46.sern_design/cad/run_junction_model.py:68)

   **対案:** 同一の鈍頭形状で、実効 `cfl_pseudo`・緩和・内反復数の感度を分離し、step 数だけで比較せず対象 CV の累積擬似時間も記録してください。併せて、EOS 前後・更新前後・温度ピン後の `ρ/P/T/ρe`、床補正量、局所閉性、面別 `V/A` を追跡します。`T_w=1000 K` と低速 IC の影響、端面間隔と隅の細層の影響は、この確認まで未排除です。

5. **Major — テーパの寸法選定と B1d の A/B は、形状効果と格子効果を混同します。**

   テーパ側面が流れに沿っても、**残るベース面は依然として x に垂直**です。`L_t=5t, t_base=t/5` は試験候補にはできますが、成立を保証する根拠はありません。機体ベースの実績は、この三重隅・加熱条件・排気／外気接続の検証を代替しません。[plan:497](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:497)

   さらに `Δx₁=t_base/5` とすると、`t_base` の A/B で CV 体積・局所刻み・解像度まで同時に変わります。「薄いベースだから改善した」のか「格子と更新が変わったから改善した」のか分離できません。[plan:500](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:500)

   **対案:** 将来テーパを試す場合も、鈍頭対照を残し、`t_base` と `Δx₁` を独立に振ってください。合否には、全域・全壁タグの床到達数、固定 CV 群の密度履歴、ベース周囲の符号付き質量収支、逆流域、壁圧・熱流束、床補正量を含めます。`check_convergence.py` の同一設定区間の VERDICT と、派生量時系列に対する `check_quasisteady.py` の VERDICT が必要です。「1800 step 完走」では足りません。

   安価な代案の評価も修正が必要です。**連続式は既に残されている**ため、「連続式だけ弱形式にする」では変更内容が定まりません。[timeIntegration_d.cu:1053](/home/sano/work/forge/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1053) 後端面だけ slip にする案は全壁の固定解除とは別ですが、隣接 no-slip 壁との共有節点には壁固定が残ります。[mesh.cpp:962](/home/sano/work/forge/solver_density_cuda/mesh/mesh.cpp:962) 丸い後縁も曲率・剥離・格子を変える形状案であり、安価な根治策とはまだ判断できません。

6. **Minor — 根拠となる場の範囲と時刻が曖昧です。**

   `run_0428` は「後端面は健全、P 2.9〜25 kPa」とされる一方、同じ restart 元について「その隅は P 54 Pa」と記されています。[plan:490](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:490)、[plan:493](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:493) 中央部と隅、または時刻の違いなら、範囲を明記すべきです。

   **対案:** 比較対象を固定の node ID・壁タグ・時刻で指定し、AWS 実行時の commit、バイナリ、実効 config、IC 作成条件を診断記録に添えてください。特に壁温ピンの順序が結論に関わるため、手元コードと AWS バイナリの対応が必要です。

**推奨は一つです。鈍頭形状を維持した診断を先行してください。**

優先順は **① IC の不連続を除去 → ② 問題 CV の毎 step の流束・更新・床補正を記録 → ③ 時間更新と熱的拘束を分離した A/B** です。これで原因を絞る前にテーパを入れると、形状変更で症状を隠しただけかどうか判定できません。

有限厚板で計算を成立させる目的は妥当ですが、§4.13・B1d の「真因確定」と「テーパ推奨」は現状では撤回すべきです。本レビューの提案は、read-only 指定に従い **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
