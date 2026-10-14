# codex 諮問 (diagnose): core-grid-twin

- **brief**: [`notes/reviews/briefs/2026-10-12-core-grid-twin.md`](../../notes/reviews/briefs/2026-10-12-core-grid-twin.md)
- **plan**: [`plans/active/tooling-nozzle-core-grid.md`](../../plans/active/tooling-nozzle-core-grid.md)
- **date**: 2026-10-12
- **commit**: `71eb2fef` (feature/nozzle-core-grid)
- **codex**: effort `high`, 3.8 min, rc=0
- **結論**: **まず `Gref`・同一起点・同じ微小摂動の双子を、`lineDtDirectionalCap` 50 と 1 の２腕、計４本で比較し、上限に依存する持続的な差の増幅があるかを判別する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 提案・解釈 | 採否と理由・対案 |
|---|---|---|
| Major | 双子の差から `|g| ≈ 1.004` の不安定固有モードを認定する | **要再検証**。§4.60 の正式判定は「判定不能」。動いている背景場での差の増幅は、固定点の固有値の測定ではない。まず「有限時間の摂動増幅」として検証する。[plan:974](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:974) |
| Major | 丸めを種にして、差が育たなければ `λ ≤ 0` とする | **却下**。丸めの差は初回だけでなく反復中にも入る。表示値 `0.000` も減衰の証拠にならない。両上限に同じ有限の初期摂動を与え、丸め・保存精度より十分大きい区間で比較する。[plan:841](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:841)、[plan:975](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:975) |
| Major | 別々に発達した `Gref`・`Grefc` の成長率差を格子幅の効果とする | **今回の主比較から外す**。既存比較は両者 `NOT CONVERGED` で、状態も違う。移送した第４腕を足しても補間と再調整の影響が入る。まず同一格子・同一起点の上限比較に絞る。[plan:758](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:758) |
| Major | 増幅域 `x/r_t = −8.4〜−7.3` は粗化区間の外にある | **却下**。区間 `−8〜−2` と重なる。`Grefc` には端の遷移区間もある。実際の局所 `Δx` と対応させるまで「細かい所だけで増幅」としない。[plan:748](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:748)、[plan:975](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:975) |

結論: **まず `Gref`・同一起点・同じ微小摂動の双子を、`lineDtDirectionalCap` 50 と 1 の２腕、計４本で比較し、上限に依存する持続的な差の増幅があるかを判別する。**

第 1 仮説: **`Gref` の停滞軌道では、上限50の反復が境界層の縁の微小差を持続的に増幅し、上限1ではその増幅が抑えられる。** 確度: **中**。ただし、流れ方向の格子幅が関与する経路は未確定。

根拠:

- 方向別の刻みはライン面を除いた面の CFL から決まり、上限は全面の CFL を使って刻みを制限する。上限変更は局所の時間項を変える介入である。[setDT_d.cu:167](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/setDT_d.cu:167)、[setDT_d.cu:194](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/setDT_d.cu:194)
- 残差の集中する壁から70〜79節点目では、上限50の刻みは point の3.3〜5.3倍。全域中央値 `0.977` を根拠に局所の刻みの影響を除外できない。[plan:809](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:809)
- 別セッションで報告された250 stepごとの2.4〜2.7倍は、有限時間の傾きに換算すると `λ ≈ 0.0035〜0.0040 /step`。これは試験を行う根拠にはなるが、`Gref` での再現証拠ではない。[plan:975](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:975)

反証条件: **十分な信号がある小振幅区間で、上限50に登録した持続増幅が見つからない、または上限1でも同程度に増幅する場合、この起点・摂動・6,000 stepについての第1仮説を退ける。** 信号が丸めの水準に埋もれた場合は反証ではなく判定不能。

第 2 仮説: 再開後の背景場の変化や一時的な増幅を、持続する増幅と取り違えている。確度: 低、未確認。

第 3 仮説: 今回は置かない。大きな内部反復の解き残しと、熱伝導の近傍結合を足すだけで停滞が緩和するという説明は、既存試験の範囲では再提案しない。[plan:841](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:841)、[plan:983](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:983)

判別 A/B: **変更する設定は `lineDtDirectionalCap: 50 ↔ 1` の１点。以下を事前登録する。**

1. **構成**  
   起点は両腕とも `case/45.isobutane_m6_d155/run_0619_cg_gref_b0_ctl/res_60000.h5`。各腕に基準場 `Q₀` と摂動場 `Q₀ + δQ₀` の２本を置く。各6,000 step、250 stepごとに保存し、初期値と最初の1 stepも記録する。同一バイナリ・同一格子・同一の実効設定を使い、SST更新・内部反復数・緩和は元のB0のままとする。

2. **摂動**  
   拘束節点を除いた既知の島に、複数の空間波長を含む固定パターンの微小な `ρE` 摂動を与える。両上限で**同じ保存量配列の組**を使う。初期振幅の案は `max|δρE/ρE| = 10⁻⁷`。これは提案値であり、起動時の同一入力・1 step再評価で測る差と保存精度を十分上回ることを確認する。上回らなければ、無信号を減衰と判定しない。

3. **測定量**  
   主指標は、固定した領域での体積重み付き温度差  
   `D_T(n) = √[Σ Vᵢ(Tᵢ⁺ − Tᵢ)² / Σ Vᵢ]`。  
   `ρ`・`ρE` も成分別に同じ定義で測る。状態差なので、残差用の `R/V` をそのまま転用しない。[cg_region_res.py:16](/home/sano/work/forge-coregrid/case/45.isobutane_m6_d155/cg_region_res.py:16)

   領域は初期状態で固定した島、縮流部、全域を保存する。最大差とその位置、壁からの物理距離、局所 `Δx`・`Δτ` は補助指標とする。最大点の移動や、外から差が流入したことを「その場所固有の成長率」と呼ばない。

4. **傾きと閾値**  
   `ln D_T` の step に対する傾きを使う。最初の1,000 stepは再開過渡として別表示する。その後の**連続1,000 step以上・5点以上**で、差が測定限界の100倍以上、かつ `max|δT/T| < 10⁻³` の区間だけを候補にする。適格区間をすべて報告し、都合のよい窓だけを選ばない。

   持続増幅の運用上の基準は、`λ ≥ 10⁻³ /step`、対数直線の `R² ≥ 0.95`、区間前半・後半の傾きがともに正で２倍以内。`ρ`・`ρE` でも同時期に増幅することを確認する。これは今回提案する検出基準であり、不安定固有値の証明条件ではない。境界付近、成分不一致、適格区間不足は判定不能とする。

→ **結果A**: 上限50だけが持続増幅の条件を満たし、上限1は同じ比較区間で信号を保ったまま減衰する。第1仮説を支持し、「単に収束が遅いだけ」という説明を退ける。  
→ **結果B**: 十分な信号の下で上限50にも持続増幅がない、または両上限で同程度に増幅する。第1仮説を退ける。  
→ 差が保存精度以下、飽和済み、短い過渡しかない場合は**判定不能**。`λ = 0` や `λ < 0` を割り当てない。

やらない方がよいこと: `Grefc` と移送した追加腕を同時投入すること、再構成・熱伝導・内部反復数も変えること、6,000 step末尾の頭打ちを飽和と断定すること。解析を再現できるよう、初期場・適合区間の場・最終場は残す。

呼び出し側の前提への異議: **双子の差の増幅、残差の停滞、格子幅依存は別々の主張である。** 今回判別するのは最初の主張の上限依存まで。成功しても「細かい格子でなぜ増幅するか」は未解明として残す。また、L-θの結果は熱伝導追加による停滞緩和を否定するが、双子の差を測っていないので「増幅が止まらなかった」とはまだ言えない。記録された判定は0〜20kで３本とも `NOT CONVERGED`。[plan:984](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:984)

不足情報: 別セッションの正式判定を止めた資格条件の詳細、`0.000` の未丸め値、双子の状態差の保存精度・検出限界。本レビューでは禁止されたrunデータを読まず、数値はplanの記録に基づく。ファイル変更・計算実行はしていない。**提案はplan未反映で、呼び出し側による事前登録が必要。**
