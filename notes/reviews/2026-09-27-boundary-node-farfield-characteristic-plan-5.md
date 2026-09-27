# codex レビュー: boundary-node-farfield-characteristic (plan)

- **plan**: [`../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md`](../../../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md)
- **stage**: `plan`
- **date**: 2026-09-27
- **commit**: `2677dc85` (feature/gap-heating-precision)
- **codex**: effort `high`, 5.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M2/m3
- **extra**: `../forge-sern-design/methods/boundary.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

目的と専用境界流束の構造は妥当です。しかし、流向を一度だけ選び直す閉包には、流入組成の誤選択と有限の流束ジャンプが残ります。GPU 実装前に §4.2–4.3 の閉包を再設計すべきです。

対象は `/home/sano/work/forge-sern-design`、HEAD `84e13d9b`。`plans/README.md`・`accepted/` を確認し、既存の静圧出口計画との重複は認めません。node 限定と非対応構成の拒否は妥当で、[検証方針](/home/sano/work/forge-sern-design/procedures/verification/README.md:7)とも整合します。ghostless A⁺ を近似前処理として維持する方針も、[既存実装](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/timeIntegration_d.cu:857)と整合します。

`run_0986/0988/0989` の残差・VERDICT 原本はこの環境に見つからず、既報の係数差や準定常性は独立に再認定していません。以下の数値は、**計画式と実装の SLAU 式を Python で直接評価した反例**であり、solver run の結果ではありません。

1. **Major — 「流向を一度だけ選び直す」処理が自己整合せず、流束も不連続になる**

   根拠は [plan:74](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:74)と、[SLAU の質量流束式:584](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:584)です。SLAU の質量流束は Uₙ,b だけでは決まらず、選んだ密度が音速・速度平均・圧力差項を変えます。

   CPG、γ=1.4、面積1、接線速度0、内部状態 `(ρᵢ,Pᵢ,Uₙ,ᵢ)=(1,1,0.6)`、自由流 `(ρ∞,P∞,Qₙ)=(0.7,2.7,0.5)` を float32 で評価すると、a∞=2.32379 なので亜音速分岐に入り、P_b=1.9091609、Uₙ,b=−0.1683811 となります。

   | 密度・組成を取る側 | ρ_b | SLAU 質量流束 |
   |---|---:|---:|
   | 最初の判定：外気側 | 0.13511491 | +0.4109993 |
   | 一度選び直した後：内部側 | 1.6494006 | −0.09861247 |

   **最終流束は流入なのに、構成状態の組成・k・ω は内部値です。** §4.3 の風上化でも、この誤った境界状態を選ぶため修復されません。密度・圧力は正で、置換カウンタも検出しません。「符号の食い違いは微小」という前提も成立していません。

   同じ条件で ρ∞ だけを変えると、最終流束は次のように跳びます。

   | ρ∞ | 最初の質量流束 | 一度選び直した最終流束 |
   |---:|---:|---:|
   | 1.8027097 | +2.3842e−7 | −0.09861247 |
   | 1.8027117 | −2.3842e−7 | −2.3842e−7 |

   入力差約2e−6に対して約0.099のジャンプです。SLAU2 も質量流束は同式なので該当します。

   **対案:** 候補流束の符号で密度側を循環選択する閉包を撤去し、入出射特性の分離から流束と輸送組成を一貫して決める構成へ改訂してください。反復回数を増やすだけでは、上の二候補に自己整合する選択がない問題を解消できません。この反例と密度・圧力・接線速度の掃引を V0u に追加し、CUDA 実装前の必須ゲートにすべきです。

2. **Major — 自由流による永久固定分類は、局所的な亜音速化・逆流時に必要な境界情報を失う**

   [plan:62](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:62)では、自由流が超音速流出なら常に `U_b=U_i` です。したがって、亜音速分岐に追加した流向修正は、この面には到達しません。

   例えば γ=1.4、ρᵢ=ρ∞=1、Pᵢ=P∞=1/γ、aᵢ=a∞=1、Qₙ=+2、実際の Uₙ,ᵢ=−0.2 とすると、

   - 分類は超音速流出のまま。
   - `U_b=U_i` なので、一様状態に対する整合性から質量流束は −0.2。
   - 流入組成は外気値ではなく内部値。

   これは [V0u(vii):137](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:137)の「質量流束が負なら外気値」という契約にも反します。逆流まで至らなくても、局所亜音速化で生じる入射音響特性を全量外挿では指定できません。境界情報の数は特性速度の向きで決まります。[NASA Wind-US の説明](https://www.grc.nasa.gov/www/winddocs/user/bc.html)

   **対案:** 局所 EOS と状態に対応する特性速度を使い、音速通過では λ⁺・λ⁻ による流束寄与が連続になる構成にしてください。「状態を硬く切り替えると不連続になる」ことへの対処を、分類の永久固定で代用すべきではありません。V0u には自由流分類と局所流況が異なる場合、V2 には局所音速通過・逆流を追加してください。

3. **Minor — ホスト単体試験から、計画中の流束関数を直接検証できない**

   [plan:107](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:107)では流束関数を `__device__` とし、[V0u:135](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:135)はホスト試験で面流束の連続性と内部面実装との一致を検証するとしています。この呼出し契約は成立していません。

   **対案:** 状態から流束までの純粋計算部分を `__host__ __device__` にし、GPU 配列・診断・残差加算を外側へ分離してください。そのうえで、ホストの閉包試験と、実際の CUDA カーネルを使う一致試験を区別する必要があります。

4. **Minor — 保存収支試験が依存する既存帳簿は、汎用の全成分診断になっていない**

   [V2b:147](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:147)は全節点・毎評価の帳簿を要求します。しかし、既存の [`ledgerCapture`:550](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:550)は化学種を `roY0/roY1` に固定し、節点・成分ごとに `cudaMemcpy` します。節点指定も[列挙式](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:497)です。

   小規模な二成分試験には使えても、計画が受け入れる多成分全般の収支を保証する診断には不足しています。

   **対案:** §5 に、実際の種数からの成分列挙、全節点指定、一括転送、評価回・段階の対応付けを追加してください。既存のソース前後の採取位置は利用できるため、診断機構を全面的に作り直す必要はありません。

5. **Minor — V3b は最大試験幅を無条件に採用できる規則になっている**

   [plan:161](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:161)の「それより広いすべての試験済み幅と比較」には、比較対象が空の場合の扱いがありません。

   例えば幅2.50／3.42／4.35 Hの C_M 差分が0／0.01／0.02なら、最初の二候補は不合格ですが、最大幅は比較相手がなく条件を満たしてしまいます。これでは「必要幅4.35 H」を根拠なく採用できます。

   **対案:** 採用候補には少なくとも一つの広幅対照を必須としてください。最大幅しか残らなければ「必要幅未確定」とし、系列を延長する規則を加えるべきです。

推奨は一つです。**専用カーネル・同時刻のスカラー面値・面流束診断という構造を維持し、境界閉包の再設計とホスト検証を GPU 実装より先に完了する。** 優先順は **①流向選択の自己整合性と連続性、②局所特性に対応する分類、③試験 API・帳簿・幅採否規則**です。

V1–V2 を機能受入れ、V3 を配置評価とする順序、dual-time の `nSub` 倍増・Δt 半減比較は維持して構いません。今回の阻害要因は試験数ではなく、受入れ条件と両立しない閉包です。

ファイルは変更していません。**plan 未反映**です。

指摘数: Critical 0 / Major 2 / Minor 3
