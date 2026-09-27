# codex レビュー: boundary-node-farfield-characteristic (plan)

- **plan**: [`../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md`](../../../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md)
- **stage**: `plan`
- **date**: 2026-09-27
- **commit**: `e4ea8df8` (feature/gap-heating-precision)
- **codex**: effort `high`, 6.9 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m2
- **extra**: `../forge-sern-design/methods/boundary.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

専用流束・同時刻のスカラー面値という構造は妥当です。しかし、境界構成が接触波の受入れ条件と両立せず、逆流時の外気スカラー供給にも欠陥があります。GPU 実装前に境界閉包を修正すべきです。

対象は `/home/sano/work/forge-sern-design`。`plans/README.md`・`accepted/` を確認し、既存の静圧出口計画との重複は認めません。node 限定・非対応構成の拒否、ghostless A⁺ を近似前処理として残す方針も妥当です。

対応する `run_0986/0988/0989/0990` の残差・VERDICT 原本は手元に見つからず、既報の収束・準定常性・力係数差は独立に再認定していません。以下の数値反例は、計画式と既存 SLAU 式を読み取り専用の Python で評価したものです。

1. **Major — 異なるエントロピーの状態から Riemann 不変量を混ぜるため、等圧の接触波を音響擾乱に変える**

   根拠: [plan:57](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:57) の構成では、亜音速流出の `R⁻` を外気から固定します。内部と外気の温度が異なると、圧力・速度が等しくても音速差が境界速度・圧力へ入ります。これは frozen-γ の近似誤差以前の問題で、**CPG でも発生します**。

   V2d の条件に合わせ、γ=1.4、R=287 J/(kg·K)、Pᵢ=P∞=2851 Pa、Tᵢ=600 K、T∞=220 K、Uₙ,ᵢ=Qₙ=0.5a∞ としました。float32 での評価結果は次のとおりです。

   | 量 | 内部の一様流・物理流束 | 計画式＋SLAU |
   |---|---:|---:|
   | 法線速度 | 148.657 m/s | 構成状態 632.868 m/s |
   | 圧力 | 2851 Pa | 構成状態 612.514 Pa |
   | 単位面積の質量流束 | 2.46122 kg/(m²·s) | 4.46563 kg/(m²·s) |
   | 圧力流束成分 p̃ | 2851 Pa | 2116.706 Pa |

   数値流束まで変わるため、「構成状態だけの違い」で済みません。流束式の根拠は [SLAU:572](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:572)、[同:584](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:584) です。

   このままでは、[V2d-1/2:139](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:139) が要求する「接触波の無圧力擾乱」「微小音響反射の測定」と整合しません。特に V2d-2 はパルスを与える前から境界が大きな擾乱を発生させます。SU2 の同型実装を参照したこと自体は正しいものの、それは異なるエントロピーを持つ流れの無反射通過を保証しません。[SU2 公式ソース](https://raw.githubusercontent.com/su2code/SU2/master/SU2_CFD/src/solvers/CEulerSolver.cpp)

   **対案:** 流出するエントロピー・組成と入射音響を分離する、局所 EOS の特性振幅に基づく閉包へ §4.2 を改訂する。まずホスト試験で「P・u 一定、T・Y のみ変化」の接触波適合性を確認する。V2d-2 は無パルス対照を必須とし、基底状態の適合性を確認してから反射率を測る。許容値を緩めて通すべきではありません。

2. **Major — 自由流では流出でも、実際に逆流した面から外気の組成・k・ω が入らない**

   根拠: [plan:64](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:64) は `Qₙ>0` なら構成状態の組成・k・ω を内部値にします。[plan:75](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:75) が負の質量流束で選ぶのも、この構成状態です。したがって、**質量流束の符号で風上化しても、外気値には切り替わりません**。

   CPG、γ=1.4、ρᵢ=ρ∞=1、Pᵢ=P∞=1/γ、aᵢ=a∞=1、Qₙ=+0.2、Uₙ,ᵢ=−0.6 では、

   - Uₙ,b=−0.2、c_b=0.92
   - SLAU 質量流束 = **−0.243786**
   - `sstEnergyIncludesK: 0`、kᵢ=10k∞ なら、流入する k の面値は **10k∞**

   となります。置換は発動せず、外気供給の失敗を検出できません。V2b の「流入配置／流出配置」だけでは、この自由流分類と実流束の不一致を網羅していません。

   **対案:** 実際の入側輸送値を外気、出側を内部とする契約を、エネルギー・組成の整合を含めて閉じる。単純に判定を `Uₙ,ᵢ` の硬い分岐へ戻すのではなく、流束の連続性も検証する。`Qₙ` と質量流束の符号の全組合せ、およびゼロ通過を V0u/V2b に追加してください。

3. **Major — V3b の隣接差だけでは、選んだ最小幅が試験した幅系列に対して許容内とは判定できない**

   根拠: [plan:146](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:146) は隣接する幅の `D≤ε` から最小幅を選びます。しかし、差は累積し得ます。

   例えば、変動幅がゼロで、幅 2.50／3.42／4.35 H の C_M の差分が順に **0／0.004／0.008** なら、隣接差は両方とも許容 0.005 内です。それでも最小幅と最大幅の差は **0.008** で不合格です。これは判定規則の反例であり、既存 run の実測値ではありません。

   **対案:** 候補幅を、それより大きい**すべての試験済み幅**と比較し、全対象量で `D≤ε` を満たす最小幅を採用する。R4d の窓条件・準定常ゲートはそのまま維持してください。

4. **Minor — 指定された熱物性 helper は、構成状態の組成を直接評価できない**

   根拠: [plan:58](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:58) が挙げる [`thermo_state_at_T`:913](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/thermo_d.cuh:913) は、`roY` とセル添字 `ic` から組成を読みます。明示的な `Y_b` を受け取らず、`__device__` 専用です。

   内部セル添字でそのまま呼ぶと、流入で必要な外気組成の内部エネルギーを評価できません。また、ホスト単体試験からの共用にもそのままでは使えません。

   **対案:** 明示的な `Y,T` を受けるホスト／デバイス共用 helper を実装前提に追加する。内部と外気の組成が異なる流入状態で、R_mix・e・h・音速の整合を検査してください。

5. **Minor — 残差試験の規格化が一部未定義で、ゼロ収支・異なる保存量を正しく判定できない**

   根拠: [V2b:132](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:132) の離散恒等式は「相対 1e-6」ですが、比較する正味流束・残差和はゼロになり得ます。続く定常収支の代表量は定義されていますが、この恒等式への適用は書かれていません。

   また、[V1(d):124](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:124) は異なる残差成分を一律に `ṁk∞` 規模で判定しており、ω や流体保存量とは次元が合いません。

   **対案:** 保存量ごとに非零の規格化量と絶対許容を固定する。離散恒等式は正味和で割らず、流束の絶対和と加算丸めの上限を使う。ホスト集計を double にしても、GPU の float32 残差積算誤差は残ることを織り込んでください。

推奨は一つです。**専用流束・同時刻スカラー面値の構造を維持し、境界閉包のホスト検証を GPU 実装より先に置く。** 優先順は **①接触波に適合する閉包、②実逆流時の外気供給、③幅の採否規則、④熱物性 API と試験規格化**です。V1–V2 を機能受入れ、V3 を配置評価と分ける順序は維持してよく、cell 回帰を追加する必要もありません。

ファイルは変更していません。**plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 2
