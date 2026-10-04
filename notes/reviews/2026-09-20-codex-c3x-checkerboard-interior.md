# codex 切り分け — 市松は内部離散化か (2026-09-20)

対象: `case/53.c3x_vane_cht` の壁熱流束 2 節点交番。弱形式 A/B (`run_0039_weakbc`) が悪化した後の続き。
実行: read-only、model `gpt-6-astra`。生ログは同名 `.log`。

## 本レビューで確定したこと

### 1. 「強制側では壁半割面の流束が解に入らない」— **正しい** (実験でも確認)

`zero_res_roe_bplane_d` ([`nodeWallDirichlet_d.cu:89`](../../solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu)) と
エネルギー行の単位行化 ([`timeIntegration_d.cu:1042`](../../solver_density_cuda/cuda_forge/timeIntegration_d.cu)) で
更新から除かれる。**直接寄与は交番の原因にならない。**

**当方の独立実験** (`run_0040_wallflux0`, `mesh.nodeIsothermalEnergyBC: 2` = 強制のまま壁半割面の熱流束だけ 0):
場の変化は run-to-run ノイズ床の **0.70–0.85 倍** = **変わらない**。

### 2. ただし「壁 BC 全体が無関係」「W–I 熱拡散が原因」は**言い過ぎ** (codex 指摘、採用)

壁温ピンは内部勾配に効き、第一内部点には他の内部面・粘性仕事・乱流熱伝導も作用する。
一般化への**反例**: **cell** (ゼロ化対象外)、**`Qw_Wall` による W–I 流束置換** (内部点側には残る、
[`viscousFlux_d.cu:261`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu))。
周期集約はゼロ化後だが、消した壁熱流束自体を復活させない。

### 3. forge と SU2 の非直交補正の差 — **有力だが 20 倍は説明できない**

同じ面平均勾配 $\bar g$、辺ベクトル $d$、面積ベクトル $S$ で

$$F/\kappa=\bar g\cdot S+a(\Delta T-\bar g\cdot d),\qquad
a_{\rm forge}=\frac{|S|^2}{|d\cdot S|},\quad a_{\rm SU2}=\frac{d\cdot S}{|d|^2}$$

([`viscousFlux_d.cu:153`](../../solver_density_cuda/cuda_forge/viscousFlux_d.cu) /
[`flow_diffusion.cpp:101`](../../.external/su2-src/SU2_CFD/src/numerics/flow/flow_diffusion.cpp))。
**直交面では一致**し、正向きの非直交面では係数比が $1/\cos^2\theta$。
対称配置で LSQ が交番を拾わず $\bar g\approx0$ になる場合も、両者の直接差分は交番を減衰させ、
**forge 側の係数はむしろ大きい**。

> **したがって「forge の減衰が弱いから 20 倍」という説明は成立しない** (Claude の筋書きを否定)。
> 実メッシュでの補正項の相殺と交番の励起を調べる必要がある。

## 推奨 A/B (1 つ)

A = 現行、B = **内部面の熱伝導補正だけ**を SU2 型へ (LSQ・物性補間・強制壁は固定)。

- 十分定常化した第一内部点 $T$ と $q$ の交番が**低下すれば**、この補正式が振幅に因果的に寄与すると確定。
- **不変・悪化なら**、この置換だけでは解消しないと確定 (LSQ 自体や他の内部結合までは棄却しない)。

→ **採用**。`space.heatCorrSU2` (既定 0) として実装。

---

## codex の回答 (原文)

**A. X は前半が正しく、原因を限定する後半は言い過ぎです。**  
対象の node・強制等温壁では、壁半割面の熱流束は壁ノードだけに加算され、[残差ゼロ化](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:89)と[エネルギー行の単位行化](/home/sano/work/forge/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1042)で更新から除かれます。この直接寄与は交番の原因になりません。

ただし「壁 BC 全体が無関係」「W–I 熱拡散が原因」とは確定しません。壁温ピンは内部勾配に効き、第一内部点には他の内部面・粘性仕事・乱流熱伝導も作用します。一般化への反例は **cell**（ゼロ化対象外）と、**`Qw_Wall` による W–I 流束置換**（内部点側には残る、[コード](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:261)）。周期集約はゼロ化後ですが、消した壁熱流束自体を復活させません。

**B. 有力候補ですが、20 倍を説明できるとは断定できません。**  
同じ面平均勾配 \(\bar g\)、辺ベクトル \(d\)、面積ベクトル \(S\) なら、熱伝導の幾何部分は

\[
F/\kappa=\bar g\cdot S+a(\Delta T-\bar g\cdot d),\qquad
a_{\rm forge}=\frac{|S|^2}{|d\cdot S|},\quad
a_{\rm SU2}=\frac{d\cdot S}{|d|^2}.
\]

[forge](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:153)と[SU2](/home/sano/work/forge/.external/su2-src/SU2_CFD/src/numerics/flow/flow_diffusion.cpp:101)は直交面で一致します。正向きの非直交面では係数比が \(1/\cos^2\theta\)。対称配置で LSQ が交番を拾わず \(\bar g\approx0\) になる場合も、両者の直接差分は交番を減衰させ、forge 側の係数はむしろ大きいです。したがって「forge の減衰が弱いから20倍」という説明は成立せず、実メッシュでの補正項の相殺と交番の励起を調べる必要があります。

**C. 次の A/B は一つだけ。**（read-only のため **plan 未反映**）  
A＝現行、B＝内部面の**熱伝導補正だけ**を SU2 型へ変更し、LSQ・物性補間・強制壁を固定：十分定常化した第一内部点 T と q の交番が低下すれば、この補正式が振幅に因果的に寄与すると確定。  
不変・悪化なら、この置換だけでは解消しないと確定し、LSQ 自体や他の内部結合までは棄却しません。
