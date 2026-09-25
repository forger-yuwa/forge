# codex レビュー: convection-slau-wall-normal-chi (plan)

- **plan**: [`plans/accepted/convection-slau-wall-normal-chi.md`](../../plans/accepted/convection-slau-wall-normal-chi.md)
- **stage**: `plan`
- **date**: 2026-09-23
- **commit**: `6156606b` (feature/sern-design)
- **codex**: effort `high`, 5.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
壁隣接面の圧力–質量結合を回復する目的は妥当です。ただし、圧力束まで同時に変更する根拠と、V1・V3 の合否基準に重大な問題があります。  
推奨は、**質量流束の係数だけを変更する opt-in 実験へ絞り、以下を実装前に修正すること**です。

確認範囲は現行コード・関連計画・台帳・一次論文です。`run_0430–0434` は台帳上 AWS にあり、この環境には実データがありません。以下の run 数値は**記録値**、式への代入値は**今回の再計算**です。`check_convergence.py` / `check_quasisteady.py` の再判定は行えておらず、収束・定常性を独立確認したとは扱いません。

既存の `convection-slau2-lowmach` は質量流束を変更せず、`discretization-node-wall-implicit-dirichlet` は運動量行の拘束を扱うため、本件は解決済み機能との重複ではありません。EOS 床の修正を別判断にする境界も妥当です。

1. **Major — 質量補充の修正に、既存の圧力束散逸を弱める変更が混入している。**

   根拠: [plan:82](/home/sano/work/forge/plans/accepted/convection-slau-wall-normal-chi.md:82) は `mdot` と `p_third` の双方へ `χ_n` を適用します。しかし、`χ_n ≥ χ` なので、質量流束の圧力拡散を増やす一方、SLAU の圧力束第3項の係数 `1−χ` は小さくなります。[実装:538](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:538)

   面 575140 の記録状態を現行式へ代入すると、`χ_n=0.5201`、質量流束は `+4.35e−9 → −1.17e−7 kg/s` と期待どおり反転します。一方、面圧力 `p_tilde` も **807 → 1070 Pa、約32.5%増**となります。これは凍結した状態での代数評価ですが、補充経路の回復と運動量側の変更が別作用であることを示します。

   また、原典は面法線 Mach への置換を実際に比較しています。Shima–Kitamura の §III.K・Fig.20 に加え、Furusawa–Kitamura 2023 は `mSLAU` の格子感度と収束悪化を報告しています。局所適用の失敗を直接証明するものではありませんが、「AUSM⁺-up と同型だから妥当」「出典確認は実装をブロックしない」という根拠にはできません。[原典](https://www.researchgate.net/publication/258474939_Parameter-Free_simple_Low-Dissipation_AUSM-Family_scheme_for_all_speeds)、[直接比較論文](https://doi.org/10.1002/fld.5183)

   **対案:** 初回は `chi_mass=χ_n`、`chi_pressure=χ` と分離し、圧力束は現行のままにする。変更した `mdot` は既存の運動量・エネルギー・species 移流へ一貫して渡す。これは検証すべき新しい組合せであり、既存スキームの保証を継承したとは主張しない。出典確認は実装前へ移す。

2. **Major — 単一面の平衡式を、壁 CV 全体の合否判定に使えない。**

   根拠: [plan:88](/home/sano/work/forge/plans/accepted/convection-slau-wall-normal-chi.md:88) は W↔W′ を無視します。しかし、`case/46.sern_design/run_0430_3d_junction_diag_icfix/` の記録では、内点面の排出 `4.345e−9 kg/s` に対して、主要な壁同士の補充は合計 `2.841e−9 kg/s`、**排出の約65%**です。無視できる項ではありません。CV 189814 には、壁同士の補充が内点面の排出を上回る記録もあります。[面別収支:558](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:558)

   一次・壁速度ゼロ・外向き `V_n,i>0` に限定しても、内点面の正確な移流項は
   \[
   A\frac{\rho_w\rho_i}{\rho_w+\rho_i}V_{n,i}
   \]
   です。記録状態の隣接量を固定した単一面モデルでも、平衡圧力比は近似式の **0.313** に対して、調和因子込みでは **0.406** になります。さらに壁同士の流束と内点状態の応答が加わります。

   **対案:** §4.1 は機序を説明する近似モデルに降格する。V1 の主判定は、**全接続面の実際の `mdot` を合計した壁 CV の連続式収支**、密度の回復、床非到達、時系列の定常性にする。単一面モデルから外れただけで「機序の理解が誤り」と判定しない。陰解法では `Δρ=−ΔtΣṁ/V` と直接同一視せず、RHS 収支と `dq` を区別する。

3. **Major — V1 は指定した初期状態で既に不合格になる。**

   根拠: [plan:151](/home/sano/work/forge/plans/accepted/convection-slau-wall-normal-chi.md:151) は `run_0430/res_1800` から開始し、CV 153797 の密度を「常に `10ρMin` 以上」と要求します。一方、記録は初期密度 **`1.711e−4`**、`roMin=1e−4` です。要求値 `1e−3` の **17.1%** しかなく、初期時点で条件を満たしません。[初期状態と床:545](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:545)

   また、後続診断で最初の破綻候補とされた **CV 153880** が監視対象から抜けています。[診断:730](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:730)

   **対案:** 同じ初期状態を維持し、「回復区間では床非到達」「事前に定めた期限までに `10ρMin` を超え、その後維持」を分ける。153880 と全壁ノードの最小密度・床到達数も監視する。入力の密度や床を変更して条件に合わせることは、この比較の交絡になるため避ける。

4. **Major — V3 は機能の受入ゲートになっておらず、比較精度も保証されていない。**

   根拠: [plan:153](/home/sano/work/forge/plans/accepted/convection-slau-wall-normal-chi.md:153) の不合格時処置は「既定化しない」だけです。既定化はもともとスコープ外なので、回帰がどれほど悪化しても opt-in を完成扱いにできてしまいます。また、[plan:157](/home/sano/work/forge/plans/accepted/convection-slau-wall-normal-chi.md:157) は収束 VERDICT の「併記」しか要求しておらず、V3 の各量に準定常判定を要求していません。

   `check_quasisteady.py` の既定許容値は drift **5%**、oscillation **10%**です。これをそのまま使うだけでは、`C_T ±0.1%` の比較精度を担保できません。[判定引数:542](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:542)

   **対案:** V3 不合格は「受入保留・設計へ戻る」とする。基準 run・commit・比較区間・壁圧分布のノルムを固定し、両側の収束判定と、**比較許容差より十分小さい変動幅・ドリフト**を各報告量で要求する。V1 の診断成功と、生産利用可能な検証完了を区別する。`OSCILLATING` を許す場合は平均±振幅で比較し、終端 dump を代表値にしない。

5. **Major — 一次の壁ノード値に関する説明を、高次の面状態へ一般化している。**

   根拠: [plan:67](/home/sano/work/forge/plans/accepted/convection-slau-wall-normal-chi.md:67) の「W↔W′ は変化なし」、同84行の「付着境界層では第3項が0」は一般には成立しません。流束が使う速度は、ピン留めした節点値そのものではなく **`interp_dispatch` 後の面状態**です。[速度再構成:185](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:185)

   亜音速域では
   \[
   \beta_++\beta_--1
   =\tfrac34(M_L-M_R)-\tfrac14(M_L^3-M_R^3)
   \]
   なので、`M_n≈0` は厳密なゼロを保証しません。さらに内部双対面の法線は、物理壁の法線と一致するとは限りません。

   V1 は一次暖機の再現だけで、SLAU2・対象形状の二次化・周期壁接合・軸対称について有効経路を確認する試験が明示されていません。

   **対案:** `χ_n` の入力を「流束が実際に消費する再構成・補正後の面速度」と明記し、上記の無変化の主張を条件付きに直す。SLAU/SLAU2 の双方で面向き反転・等状態保存・非対象面不変を確認し、V1 成功後に同形状の二次本段まで検証する。周期・軸対称は適用マスクと保存性を確認する小規模 node 試験を加える。cell は明示的に無効化すればよく、現行運用に反する cell 生産回帰までは不要です。

6. **Minor — V2 の「場がビット不一致なら実装誤り」は強すぎる。**

   根拠: 残差は面スレッドから浮動小数点 `atomicAdd` で集積します。加算順序を固定していないため、演算式を変えなくても実行順序による差が生じ得ます。[残差集積:623](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:623)

   **対案:** flag 0 の**面流束演算のビット同一性**と、並列集積後の場の回帰を分ける。場については同一バイナリ反復の再現性を先に測り、比較条件を固定する。`check_field_regression` は今回の探索では実体を確認できなかったため、実行可能なツール名・引数・比較対象も plan に書く。

7. **Minor — 実装の変更ファイル一覧と docs 更新順が不足している。**

   根拠: `SLAU_d` の現在の引数には壁フラグがなく、呼出し元は [convectiveFlux_d.cu:260](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:260) です。したがって、§5.1・§7 に列挙した3ファイルだけでは配線できません。docs 更新も実装後になっています。

   **対案:** wrapper を変更対象へ追加し、`node && nodeWallDirichlet==1` の実効条件、無効構成の扱い、`wall_flag_d` の受渡しを明記する。`methods/convection/{theory,implementation}.md` と設定リファレンスを実装前に同期する。block-DPLUR は既に近似 FVS Jacobian と壁行拘束を使うため、今回ただちに Jacobian の厳密微分へ変更する必要はありませんが、CFL・内反復数による固定点と収束性の確認は検証計画に入れる。[Jacobian 定義:22](/home/sano/work/forge/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:22)

**推奨は、圧力束を保持した質量流束限定の opt-in 修正です。** 実装前の優先順は、①文献を踏まえた作用項の限定、②全 CV 収支による V1 と初期密度条件の修正、③V3 の受入条件・比較精度・試験範囲の固定、④配線と docs の明記です。格子品質・局所壁解像の既存 VERDICT、基準 run の実データも検証投入前に揃えてください。

ファイルは変更していません。上記提案は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 2
