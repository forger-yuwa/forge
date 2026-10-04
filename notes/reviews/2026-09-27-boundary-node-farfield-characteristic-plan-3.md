# codex レビュー: boundary-node-farfield-characteristic (plan)

- **plan**: [`../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md`](../../../forge-sern-design/plans/active/boundary-node-farfield-characteristic.md)
- **stage**: `plan`
- **date**: 2026-09-27
- **commit**: `e4ea8df8` (feature/gap-heating-precision)
- **codex**: effort `high`, 6.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M6/m1
- **extra**: `../forge-sern-design/methods/boundary.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

専用境界流束と同時刻のスカラー面値という構造は妥当です。しかし、TP の分岐に不連続が残り、SLAU・SST と両立しない合格条件があります。現版のまま実装すると、正しい実装を不合格にする検証になります。

対象 worktree は `/home/sano/work/forge-sern-design`。ファイルは変更していません。`plans/README.md` と `accepted/` を確認した限り、既存の静圧出口計画とは目的が異なり、重複ではありません。node 限定・非対応構成の起動拒否も現行方針と整合します。

手元の `case/46.sern_design/r4d_view/z2p50H/res_20000.h5` では、最外面 51,143 節点の **T=176.346–596.966 K、Y0=0–0.128410、a²ρ/P=1.372255–1.404682** を確認しました。ただし、対応する `run_0986/0988/0989/0990` の残差履歴・VERDICT 原本は見つからず、既報の収束・準定常性・力係数差は独立に再認定していません。以下の数値反例は、計画式と既存流束式を Python/float32 で評価したものです。

1. **Major — TP では特性の分岐が時間変化し、有限の流束ジャンプが残る**

   根拠: [plan:57](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:57) は γ*=γᵢ、c∞=√(γ*P∞/ρ∞) と定義する一方、[plan:60](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:60) は判定を「面ごとに固定」としています。**Qₙ が固定でも、比較対象の c∞ は内部温度・組成で変わります。**

   ρᵢ=ρ∞=1、Pᵢ=P∞=1、Uₙ,ᵢ=0.5、Qₙ=√1.39、面積=1 とし、計画式と既存 SLAU の質量流束式を評価しました。

   | γ* | 分岐 | Uₙ,b | P_b | 質量流束 |
   |---:|---|---:|---:|---:|
   | 1.389999 | 超音速流出 | 0.500000 | 1.000000 | 0.500000 |
   | 1.390001 | 亜音速 | 0.839491 | 0.662370 | 0.668000 |

   γ の差約 2×10⁻⁶ に対し、質量流束は約 **34%** 跳びます。近似 Riemann 流束を挟んでも解消しません。式の根拠は [convectiveFlux_slau_d.inc.cuh:541](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:541)。

   SU2 の引用箇所は固定の `Gamma` を使っています。その分岐を内部 γ が変動する TP に拡張しても、同じ性質は保証されません。[SU2 公式ソース](https://raw.githubusercontent.com/su2code/SU2/master/SU2_CFD/src/solvers/CEulerSolver.cpp)

   **対案:** 分類用音速を自由流の実物性から計算した固定値として保持し、再構成用の frozen-γ 音速と明確に分ける。V0u には Uₙ,ᵢ だけでなく、**Tᵢ・Yᵢ の変化で分類閾値を横切る試験**を追加する。

2. **Major — V0u の超音速流入極限は、採用する SLAU の式と両立しない**

   根拠: [plan:120](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:120) は超音速流入で F=F(U∞) を要求します。しかし SLAU は、両側が超音速流入でも一般に自由流の物理流束そのものにはなりません。[質量流束:584](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:584)

   CPG、γ=1.4、ρᵢ=ρ∞=1、Pᵢ=P∞=1/γ、cᵢ=c∞=1、Uₙ,ᵢ=−2、Uₙ,∞=−3 では、構成状態は正しく U_b=U∞ になります。それでも、

   - SLAU の質量流束: **−2.5**
   - 自由流の物理質量流束: **−3.0**
   - SLAU のエネルギー流束: **−17.5**
   - 自由流の物理エネルギー流束: **−21.0**

   となります。相対 10⁻⁶ の条件には原理的に届きません。

   **対案:** 弱形式 SLAU 境界という設計を維持し、期待値を **U_b=U∞ および F_num(Uᵢ,U∞)** に修正する。物理流束との一致は Uᵢ=U∞ の整合性試験で要求する。「構成状態の極限」と「数値流束の極限」を分けてください。

3. **Major — V1(d) は、消滅項を持つ SST に一様定常解を要求している**

   根拠: [plan:123](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:123) は 2000 step 後も k・ω が自由流値から相対 10⁻⁵ 以内とします。しかし一様流では生成・交差拡散が消えても、

   `Sₖ = −0.09ρkω`、`Sω = −βρω²`

   が残ります。[ransSource_d.cu:210](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransSource_d.cu:210)、[同:232](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransSource_d.cu:232)、[同:248](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransSource_d.cu:248)。

   例えば ρ=k=V=1、ω=10、F1=0 なら、初回から残差に **−0.9、−8.28** が加算されます。境界が正しくても k・ω は変化します。`sstEnergyIncludesK: 1` なら、その変化は平均流エネルギーにも戻ります。

   **対案:** V1(d) は一様状態での**輸送残差の打消し**をソース前に検査する試験に変更する。完全な SST の時間発展・定常解は、ソース込みの V2b で検証し、k・ω の一様値保持を要求しない。

4. **Major — V2b のソース符号が逆で、規格化にもゼロ除算がある**

   根拠: [plan:132](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:132) は `Σ面流束 + Σ体積ソース = 0` とします。コードは外向き流束を残差へ **−F**、生成を正とするソースを **+S** で加算しています。[scalarTransport_d.cu:174](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/scalarTransport_d.cu:174)、[ransSource_d.cu:247](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransSource_d.cu:247)。

   したがって定常収支は **ΣF − ΣS = 0** です。帳簿の段間差をそのまま S として使うなら、現計画は符号が逆です。

   また、自由流で Y_EXH=0 の試験では、その化学種の流入流束もゼロです。[plan:133](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:133) の分母では判定できません。横方向運動量にも同様の問題があります。

   **対案:** 保存量ごとに `R = −ΣF + ΣS` を正本として固定する。規格化は事前登録した非零の代表量と絶対許容差を併用し、ゼロ成分を判定から落とさない。`sstEnergyIncludesK: 1` の `res_roe` は既に全エネルギーの残差なので、**`res_roK` をさらに加えない**ことも帳簿仕様に明記する。[分割保持仕様:75](/home/sano/work/forge-sern-design/plans/accepted/turbulence-sst-energy-includes-k.md:75)

5. **Major — V2d は初期値と誤差予算が未定義で、TP 近似の受入れ試験として再現できない**

   根拠: [plan:137](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:137) は塊の T・Y を指定していますが、圧力・速度・密度の整合条件、形状、格子幅、参照解の収束性を定義していません。

   圧力・速度を連続にした温度／組成の塊なら、連続 Euler 系では移流する接触波であり、**物理的な圧力擾乱はゼロ**です。その場合 `ΔP_scale` はゼロ、または参照コードの数値誤差になります。逆に圧力差を持たせれば別の問題です。さらに T を比較量に挙げながら、T の合否条件がありません。

   **対案:** 完全な初期状態と評価窓を先に固定し、接触波と音響擾乱を別試験にする。接触波の圧力誤差は P∞ で規格化し、T・Y にも許容値を置く。独立 HLLC 参照は格子・時間刻みを細分化して誤差を評価し、参照誤差を受入れ許容の 1/5 以下にする。共通領域の forge 長領域対照も使い、内部離散化誤差と境界誤差を分離してください。

6. **Major — 非定常試験の内部反復条件が現行手順に反し、宣言した時間積分経路も検証されない**

   根拠: [plan:128](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:128) の「3 桁低下または 30 回」は、残差が十分下がらなくても 30 回で評価を進めます。現行手順は、**残差低下だけでは不十分で、`nSub` 倍増時に同じ物理時刻の解が変わらないことを確認する**よう明記しています。[recommended-settings.md:256](/home/sano/work/forge-sern-design/procedures/recommended-settings.md:256)

   `check_convergence.py` も通常は `outer_end` を選ぶため、その VERDICT だけで各物理ステップの内部反復を保証できません。[check_convergence.py:40](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_convergence.py:40)

   また、定常陽解法と SST を含む dual-time の試験が明示されていません。後者は `sstEnergyIncludesK` の処理が定常系と異なる実経路です。[main.cpp:1979](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1979)

   **対案:** V2a/V2d に `nSub` 倍増・Δt 半減の比較を追加し、観測量の変化を許容の 1/5 以下にする。非定常の精度判定を定常用 VERDICT と区別して保存する。対応範囲を維持するなら、定常陽解法と SST dual-time の `sstEnergyIncludesK` 0/1 を検証表に追加する。

7. **Minor — 診断流束と SERN 帳簿の圧力基準が未定義**

   根拠: 既存 SLAU が残差に投入する運動量流束は **p̃−`d_pRef`** を使います。[convectiveFlux_slau_d.inc.cuh:615](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:615)。一方、SERN 帳簿は **P−p_a** です。[sern_momentum.py:62](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_momentum.py:62)

   §4.3 のダンプを両用途へ直接流用すると、`pRef ≠ p_a` の場合に境界群ごとの寄与がずれます。閉曲面全体で定数圧力が相殺することは、部分境界の混在を正当化しません。

   **対案:** ダンプに面積ベクトル、圧力基準、残差投入流束を記録する。SERN 側へ渡す際は運動量流束に `(pRef−p_a)S` を加え、基準を統一する。

推奨は一つです。**専用流束・同時刻スカラー面値の構造を維持し、上記を計画へ反映してから再レビューする。** 優先順は **①TP 分岐、②SLAU 極限、③SST 自由流試験、④保存収支、⑤独立参照と時間精度、⑥診断仕様**です。ghostless A⁺ を近似前処理として残す方針自体は許容できますが、今回の不連続や合否条件の矛盾を解消するものではありません。

ユーザー指定に従い、**plan 未反映**です。

指摘数: Critical 0 / Major 6 / Minor 1
