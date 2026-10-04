# codex レビュー: boundary-node-farfield-characteristic (plan)

- **plan**: [`../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md`](../../../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md)
- **stage**: `plan`
- **date**: 2026-09-27
- **commit**: `2677dc85` (feature/gap-heating-precision)
- **codex**: effort `high`, 7.5 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m1
- **extra**: `../forge-sern-design/methods/boundary.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

専用境界流束・同時刻のスカラー面値・面流束診断という構造は妥当です。
しかし、自由流をそのまま外側に置く SLAU は、計画が主張する特性の分離を満たしません。音響反射と TP 接触波の受入れにも、実装前に解消すべき問題があります。

対象は `/home/sano/work/forge-sern-design`、HEAD `b953c828` です。`plans/README.md` と `accepted/` に同機能の実装済み計画は見つかりません。node 限定・非対応構成の拒否は[検証方針](/home/sano/work/forge-sern-design/procedures/verification/README.md:7)に沿い、ghostless A⁺ を近似前処理として維持する方針も[既存実装](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/timeIntegration_d.cu:858)と整合します。

ただし、`run_0986/0988/0989` の場・残差・VERDICT 原本はローカルにありません。SERN の係数差と準定常性は台帳の記録までの確認で、独立に再認定していません。以下の数値は、**コードの流束式・既存物性コードを読み取り専用の Python で評価したもの**です。forge 本体の run 結果とは区別します。

1. **Major — 固定自由流＋SLAU は、出射特性を内部状態だけで決める境界になっていない**

   根拠は [plan:61](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:61) の「局所の波速で振り分ける」という主張と、[SLAU の速度平均:541](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:541)、[質量流束:584](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:584)です。

   γ = 1.4、面積1、接線速度0、左右とも ρ = 1、P = 1/γ、音速1 とし、内部 Uₙ,L = 2、外側 Uₙ,R = 3 とします。両側とも超音速流出で、この Riemann 問題の波はすべて外向きです。それでもコードの式は次を返します。

   | 流束成分 | 内部状態の物理流束 | 計画の SLAU 流束 |
   |---|---:|---:|
   | 質量 | 2.000000 | 2.500000 |
   | 法線運動量 | 4.714286 | 5.714286 |
   | エネルギー | 9.000000 | 11.250000 |

   **外から情報が入らない条件でも、外側速度が質量・エネルギー流束を25%変えます。** SLAU2 も質量流束は同じです。流束の連続性と、入出射特性の正しい分離は別の性質です。

   また、[V2a:128](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:128) の5%反射基準にも強い不適合の兆候があります。コードの流束を M = 0.3 の一様状態で線形化し、境界半 CV を持つ1次元 node モデルで確認しました。内部は2次再構成、境界は計画どおり1次、RK3、音響 CFL 0.25。長さ1 mと3 mの差を、同じ観測点の入射振幅で規格化した結果です。

   | Δx | SLAU | SLAU2 |
   |---|---:|---:|
   | 5 mm | 16.96% | 14.01% |
   | 2.5 mm | 17.78% | 14.64% |
   | 1.25 mm | 17.94% | 14.77% |

   最細格子で時間刻みを半減しても、反射率の変化は約1e−6でした。これは forge の V2a 不合格を直接認定する結果ではありませんが、**空間・時間精度を上げれば5%に入るという見込みを支持しません。**

   **対案:** 内部面の SLAU/SLAU2 は維持し、遠方境界だけを EOS 整合な HLLC 流束に変更することを推奨します。超音速極限・接触波流束・小振幅音響をホストで先に検証し、その合格を CUDA・runner・大規模診断実装の前提にしてください。現在の「全体を実装してから V2a で判断する」順序は費用対効果が悪いです。

2. **Major — TP の「接触波流束が正しい」ことから「圧力擾乱を作らない」は導けない**

   [plan:63](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:63) の質量流束の説明は、基本 SLAU の面評価として成立します。しかし、TP で保存量を更新し、非線形 EOS から圧力を復元するところまでの圧力平衡は保証しません。

   根拠は [NASA エンタルピーを使う流束](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:400)と、[e = h − RT の物性定義](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/thermo_d.cuh:339)です。既存の [`FrozenGas`](/home/sano/work/forge-sern-design/design/forge_design/gas/frozen.py:75) と SERN の実組成で、等圧・等速度の二状態を保存量として混合して確認しました。

   - 高温側：P = 2851 Pa、T = 600 K、Y_EXH = 0.13
   - 外気側：P = 2851 Pa、T = 220 K、Y_EXH = 0
   - 更新：U_new = 0.9 U_hot + 0.1 U_air

   復元結果は **T = 514.973 K、Y_EXH = 0.0992112、P = 2869.641 Pa**。圧力誤差は **+0.65385%** です。これは等速度の接触面を1次風上で移流した際に生じる保存量更新の形であり、面流束が風上側の物理流束でも発生します。単一の TP 空気でも同じ計算で **+0.51328%** でした。多成分保存形での同種の問題は、[Abgrall の原論文](https://www.sciencedirect.com/science/article/abs/pii/S0021999196900856)でも扱われています。

   [V2d:142](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:142)には長領域との比較が既にありますが、絶対圧力誤差の合格条件は別に残っています。この単一更新の反例だけでガウス塊の試験結果までは断定できません。しかし、**境界閉包だけで TP の圧力平衡まで保証できるという前提は誤り**です。HLLC への変更だけでも、この問題は解消しません。

   **対案:** V0u の「面流束の接触波整合」と、時間発展後の「圧力平衡」を明確に分離してください。単一 CV の保存量更新＋EOS 復元を前提試験に加え、V2d は長領域自身の圧力誤差と短領域の追加誤差を別々に判定するべきです。長領域が絶対基準を満たさなければ、既存 TP 輸送の精度問題を依存課題として扱い、farfield の変更だけで達成済みにしない契約が必要です。

3. **Major — V2f の逆流判定が、採用した質量流束の契約と食い違う**

   [V2f:137](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:137) は、内部法線速度を負にした帯で外気組成が使われることを要求しています。一方、[§4.3:74](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:74) は、**数値質量流束**が負の場合だけ外気組成を使う契約です。

   V2f に対応する γ = 1.4、ρᵢ = ρ∞ = 1、Pᵢ = P∞ = 1/γ、a = 1、Qₙ = 2、Uₙ,ᵢ = −0.2 を評価すると、

   **ṁ = +0.71999997**

   です。内部速度は逆向きでも、数値流束は流出です。Yᵢ = 0.13、Y∞ = 0 なら、§4.3 に従う種流束は **+0.09359999** で、内部組成を使います。外気組成を強制すれば、むしろ採用した風上化契約を破ります。

   **対案:** V2f の判定対象を「内部速度が負の面」から「実際の ṁ が負の面」に修正してください。上の入力は「内部速度と数値流束の符号が異なる試験」として残し、別に ṁ < 0 が確認できる入力を用意して外気組成の注入を検証するべきです。

4. **Minor — 最新方針が実装手順と仕様一覧へ反映されていない**

   次の不整合が残っています。

   - [§5:93](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:93)：撤去した「境界状態の構成」が実装手順に残り、流束関数も `__device__` 表記のまま。
   - [§6:119](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:119)：存在しなくなった c_b ≤ 0 の置換を参照。
   - [methods/boundary.md:37](/home/sano/work/forge-sern-design/methods/boundary.md:37)、[plans/README.md:28](/home/sano/work/forge-sern-design/plans/README.md:28)：旧 Riemann 不変量方式を記載。

   **対案:** 実装前に、最新の流束定義、`__host__ __device__` API、非有限流束カウンタ、影響ファイル一覧へ統一してください。帳簿拡張も §5 の独立した実装項目として明記すべきです。

推奨は一つです。**境界専用 HLLC を候補として、ホストでの成立性確認を全体実装より先に完了する。** 優先順は、①超音速極限と音響5%基準、②TP の面評価と時間発展の誤差分離、③V2f の判定修正、④仕様・実装手順の同期です。

同時刻のスカラー面値、SST ソース込み収支、dual-time の反復数・時間刻み感度、V1–V2 の後に V3 を行う順序、広い全候補との幅比較は維持して構いません。ファイルは変更していません。**plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 1
