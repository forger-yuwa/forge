# codex レビュー: boundary-node-farfield-characteristic (plan)

- **plan**: [`../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md`](../../../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md)
- **stage**: `plan`
- **date**: 2026-09-28
- **commit**: `2677dc85` (feature/gap-heating-precision)
- **codex**: effort `high`, 6.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M2/m1
- **extra**: `../forge-sern-design/methods/boundary.md`, `../forge-sern-design/solver_density_cuda/tools/farfield_proto1d.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

境界専用 HLLC・同時刻のスカラー風上化・面流束診断という構造は妥当です。
ただし、温度差のある内外状態で音響反射が受入れ基準を超え、TP 用の波速推定にも有限値のまま流束が不連続になる反例があります。
この2点をホスト段階で解消するまで、本体実装への着手は推奨しません。

対象は `/home/sano/work/forge-sern-design`、HEAD `a0b59649` です。既存の静圧出口・ghostless 境界計画との重複は認めません。`side_far` が slip に固定されていること、`top_out_kind: outflow` が slip に落ちることも[runner:60](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern3d.py:60)で確認しました。node 限定・非対応構成の拒否、ghostless A⁺ を近似前処理として維持する方針は既存構造と整合します。

以下の数値は、添付試作と既存物性コードを使った**読み取り専用のホスト評価**です。forge 本体の検証結果ではありません。SERN の `run_0986/0988/0989` は台帳までの確認で、残差・場・VERDICT 原本による独立再認定はしていません。

1. **Major — 固定自由流＋HLLC は、V2d-2 の温度差を伴う音響5%基準を満たさない**

   根拠は、外側状態を固定する[plan:69](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:69)、異なる内外温度で反射率5%以下を要求する[plan:154](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:154)、同じ背景状態だけを音響試験に使う[`farfield_proto1d.py`:157](/home/sano/work/forge-sern-design/solver_density_cuda/tools/farfield_proto1d.py:157)です。

   試作の `run()` と `hllc_ghost` をそのまま使い、初期状態だけを次に変更しました。

   - CPG、γ = 1.4、外側 ρ∞ = 1、P∞ = 1/γ、c∞ = 1。
   - 内部 ρᵢ = 220/600、Pᵢ = P∞。同一気体で Tᵢ/T∞ = 600/220 に相当。
   - 内外とも u = 0.5。内部音速に整合する右向き純音波、振幅 10⁻³P∞、半値幅40 mm。
   - 試作と同じ短領域 `[-1,1]`・長領域 `[-1,3]`、観測点 x = 0.8。

   | Δx | 短領域−長領域の反射振幅／入射振幅 |
   |---|---:|
   | 5 mm | 24.40% |
   | 2.5 mm | 29.79% |
   | 1.25 mm | 33.48% |

   2.5 mm で時間刻みを半減しても **29.793%** です。パルスなし対照の圧力誤差は、この評価では **0** でした。

   したがって、**接触状態が無擾乱で保たれることと、その状態を通る音波が無反射で出ることは別です。** この反例は CPG なので、既知の TP 保存形混合誤差では説明できません。長領域との差を取る現在の判定でも残ります。均一背景で得た0.04〜0.5%を、温度・組成が違う背景へ一般化できません。

   **対案:** V2d-2 を本体実装前のホストゲートへ移してください。5%を受入れ要件として維持するなら、局所の出射音響特性に整合する境界閉包が必要です。固定 U∞ をそのまま渡す方式を、その要件を満たすものとして確定してはいけません。CPG の上記反例から修正を検証し、その後に単成分・多成分 TP へ進める順序を推奨します。

2. **Major — TP 用の音速算術平均は HLLC の波速順序を壊し、流束の連続性も保証しない**

   [plan:73](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:73)では c̃ = (cL+cR)/2 としています。しかし、試作の合格結果は[Roe エンタルピー平均から音速を求める別の式](/home/sano/work/forge-sern-design/solver_density_cuda/tools/farfield_proto1d.py:61)によるものです。

   既存の[`FrozenGas.air()`](/home/sano/work/forge-sern-design/design/forge_design/gas/frozen.py:57)で次の TP 空気を評価しました。温度・圧力とも正で有限です。

   | 状態 | T [K] | P [Pa] | ρ [kg/m³] | c [m/s] |
   |---|---:|---:|---:|---:|
   | L | 2200 | 100000 | 0.15832211 | 904.26248 |
   | R | 220 | 10 | 0.00015832211 | 297.49478 |

   uL = uR = 0 として計画の式を適用すると、

   `SL = −904.26248, S* = 697.96235, SR = 600.87863 m/s`

   となり、**S* > SR** です。右星状態の密度も **−0.00097990035 kg/m³** になります。「波速の見積もりだけ」という変更でも、HLLC の成立条件に直接影響します。HLLC の正値性が波速選択に依存することは、[Batten らの原論文](https://doi.org/10.1137/S1064827593260140)でも明示されています。

   さらに、両側へ同じ速度を加え、添付試作の[分岐:82](/home/sano/work/forge-sern-design/solver_density_cuda/tools/farfield_proto1d.py:82)にこの TP 波速を適用すると、次の不連続が出ます。

   | uL = uR [m/s] | 質量流束 [kg/(m²·s)] |
   |---|---:|
   | −697.962352651736 | −0.110502872 |
   | −697.962350651736 | +8.94×10⁻⁸ |

   速度差 **2×10⁻⁶ m/s** に対し、質量流束が約0.1105跳びます。**値は有限なので、非有限流束の置換カウンタでは検出できません。** 分岐順だけを変えても、波速順序の破綻は解消しません。

   **対案:** 算術平均の採用を撤回し、TP に対応する波速評価を選び直してください。ホストゲートには `SL ≤ S* ≤ SR`、星状態の密度・EOS 上の admissibility、共通速度を加えた掃引を追加する必要があります。不適格時は検証済みの保守的な退避流束へ切り替えるか明示的に失敗させ、有限値も検査対象にしてください。CPG で試した波速と TP 本番の波速を分けたまま「全項目合格」と扱うべきではありません。

3. **Minor — 仕様一覧と影響範囲に旧方式・記載漏れが残る**

   [methods/boundary.md:37](/home/sano/work/forge-sern-design/methods/boundary.md:37)の一覧は、現在も「Riemann 不変量を混合して境界状態を構成する」方式です。一方、本文は固定自由流＋HLLC に変わっています。また、[plan §7:176](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:176)には、§5で変更する対流カーネル、`scalarTransport_d.cu`、`passiveKernels_d.cuh`、追加ツールが載っていません。

   **対案:** 閉包修正と同時に、仕様一覧・§5・§7を同期してください。[methods:142](/home/sano/work/forge-sern-design/methods/boundary.md:142)の「接触波は擾乱なしに通る」も、面流束の整合性と TP の時間発展後の圧力誤差を区別する記述に直すべきです。

推奨は一つです。**専用 HLLC カーネルと同時刻スカラー輸送の構造は維持し、境界閉包の成立性確認を本体実装より先に完了する。** 実施順は、①温度差を伴う音響5%基準、②TP の波速順序・連続性・状態妥当性、③仕様同期、④CUDA・収支診断・runner の実装です。

保存収支の絶対和による規格化、SST の `res_roe` に `res_roK` を重ねて足さない扱い、dual-time の反復数・時間刻み感度、V1–V2の後にV3を行う順序、広幅対照を必須とする幅採否規則は維持して構いません。

ファイルは変更していません。**plan 未反映**です。

指摘数: Critical 0 / Major 2 / Minor 1
