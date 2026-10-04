# codex 諮問 (diagnose): twophase-diffusion-4c

- **brief**: [`notes/reviews/briefs/2026-10-02-twophase-diffusion-4c.md`](../../notes/reviews/briefs/2026-10-02-twophase-diffusion-4c.md)
- **plan**: [`plans/active/condensation-two-phase-transport.md`](../../plans/active/condensation-two-phase-transport.md)
- **date**: 2026-10-02
- **commit**: `09299642` (feature/species-transport)
- **codex**: effort `high`, 4.1 min, rc=0
- **結論**: **独立残差監査と変更可能な緩和係数を備えた定常専用初版カーネルに着手し、上記の判定修正と制限作動A/Bを受入条件にする。**
- **extra**: `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 論点 | 重大度・採否 | 根拠と対案 |
|---|---|---|
| Q1：定常専用初版への着手 | **採用。ω=1、2×2は必須にしない** | [設計メモ:453](/home/sano/work/forge-species/notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md:453)の記録は事前比較条件を満たす。ただしω=1は初期選択とし、緩和係数を変更できる実装にする。高CFLでの成功やCFDへの適用まで保証したとは扱わない。dual-time併用拒否は維持する。 |
| Q2：制限作動試験を着手前の追加ゲートにするか | **却下。実装の受入試験には採用** | [設計メモ:460](/home/sano/work/forge-species/notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md:460)ではθ・θ_srcが全期間1。これは未検証範囲であり、既定条件の失敗ではない。下記A/Bを受入条件に加える。CFL増大より、更新制限の閾値だけを変える方が切り分けやすい。 |
| Q3：核生成の自己Jacobian追加 | **Major／却下。初版は緩和で扱う** | [ソース:130](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/condensationSourceF_d.cuh:130)では核生成率をモーメントによらず評価し、139行で`SQ0 = J`としている。核生成分の∂S_Q0/∂ρQ0は0であり、`sj_Q0=0`は欠落ではない。T・蒸気・他モーメントを介する結合は非対角項である。必要なら後続で連成を検討し、架空の自己微分を追加しない。 |
| Q4：独立残差による停止契約 | **Major／採用** | [設計メモ:494](/home/sano/work/forge-species/notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md:494)ではf32停止が蒸気残差の最大12%超過を見逃している。格納状態からEOS・流束・ソースを再評価し、蒸気を含む全成分で合否を決める。既に丸めたf32残差をdoubleへ変換するだけでは足りない。 |
| 「最後の10%」の補正監視 | **Major／修正を採用** | [ハーネス:323](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_real_source.cpp:323)は`0.9*o.cap`を使用するため、上限20000では18000回以降しか集計しない。69〜161回で終了すると監視窓は空になる。反復ごとの補正量を記録し、**実際の終了反復数**から末尾10%を集計する。今回の「全期間補正0」という記録が正しければ、今回の合否は変わらない。 |
| 異常値の合否判定 | **Major／修正を採用** | [ハーネス:304](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_real_source.cpp:304)の`std::max`による集計はNaNを取り落とし得る。また[100行](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_real_source.cpp:100)では表範囲外でゼロソースを返し、失敗を伝播しない。状態・残差・尺度・増分の非有限値と、未対応の表範囲外を明示的にFAILにする。 |

結論: **独立残差監査と変更可能な緩和係数を備えた定常専用初版カーネルに着手し、上記の判定修正と制限作動A/Bを受入条件にする。**

第 1 仮説: **実ソースと非分割更新の組合せは、更新制限が途中で作動しても、最後に制限を解除して原方程式残差の基準へ到達できる。** 確度: **中。制限作動時は未確認。**

  根拠: ブリーフの固定条件では全4 variantが残差基準を満たし、補正0。[ハーネス:348](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_real_source.cpp:348)以降では制限を増分に掛け、停止は[298行](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_real_source.cpp:298)以降の独立残差で判断する。ただし、この構造だけでは反復の到達性は保証されない。

  反証条件: 制限を作動させたBだけが、θ=0の停止、残差未達、末尾の補正継続、または制限未解除になること。

第 2・第 3 仮説: **高CFLの未収束は、連成反復またはfloatソース評価の量子化による振動である可能性がある。確度: 中、機構は未確認。** Q0の110〜160 ULPという増分は「全増分が格納時に消える停滞」を支持しないが、浮動小数点由来の周期運動までは除外しない。

判別 A/B: **変更は`DT_MAX`だけ**。同じ初期値、CFL=5、ω=1、対角前処理、上限20000回で、A=1 K／B=0.01 K。先に補正監視と異常値判定を直し、両者で同じ判定を使う。

- Bでθ<1の発生を確認し、蒸気を含む全成分の独立残差比≤1、非負、実反復数の末尾10%で補正0、最後のθ・θ_src=1を確認する。
- **A・Bとも達成なら、当該条件で更新制限が反復を妨げる懸念を除外する。Aだけ達成なら、第1仮説を棄却して制限・commitを調べる。**
- Bでもθが常に1なら判別不成立。この試験はθ_srcの作動検証にはならない。

やらない方がよいこと: **`sj_Q0`へ根拠のない対角を足すこと、停止許容を広げて高CFLを合格させること、ω=1を全条件へ固定すること。**

呼び出し側の前提への異議:

- **Major：「2×2で差がないから蒸気側結合が弱い」は要再検証。** [ハーネス:158](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_real_source.cpp:158)の蒸気摂動はT固定であり、液側は[149行](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_real_source.cpp:149)の温度応答近似を使う。EOSを含む蒸気・液座標での完全なJacobianではない。「試した近似ブロックでは改善しなかった」と記録する。対角版が合格しているため、初版で必須にしない判断は維持できる。
- **Minor：1DのCFLとCFDの`cfl_pseudo`は数値を直接比較できない。** [ハーネス:285](/home/sano/work/forge-species/solver_density_cuda/tests/unit/test_twophase_real_source.cpp:285)はM=U/CFL。CFDへ対応づける際は、名称でなく実際のΔτと輸送・ソース対角の比を比較する。

不足情報: 制限作動時の実ソース接続結果、CUDAへの移植一致、高CFL周期運動のセル別状態履歴。今回はコードと提示記録を監査し、数値試験は再実行していない。禁止されたrun出力は未読、ファイル変更なし。**plan未反映**。呼び出し側で[plan](/home/sano/work/forge-species/plans/active/condensation-two-phase-transport.md) §4.2・§5.1 #4/#4c・§6へ採否と受入条件を反映すること。
