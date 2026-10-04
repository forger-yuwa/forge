# codex 診断: C3X 翼列 (case/53) の前縁発散 — node 等温壁の保存量が凍結している

- 実行: 2026-09-20, `codex exec --sandbox read-only` (自由形式の質問。plan レビューではない)
- 質問文: scratchpad/codex_q.txt (要旨は本文冒頭)
- 生ログ: この応答の全文 (以下)

---

**1. 実際に何が壊れているか**

まず、`convMethod: 1, limiter: 0` は**リミタなしの2次再構成**です。1次ではありません。

[再構成の分岐](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_common_d.cuh:171)は `convMethod=1` を2次へ送り、[limiter の処理](/home/sano/work/forge/solver_density_cuda/cuda_forge/limiter_d.cu:450)は `limiter=0` で係数を **1** にします。[setup_run.py](/home/sano/work/forge/case/53.c3x_vane_cht/tools/setup_run.py:47)は warm・sst_soft・sst_mid をすべてこの設定にしています。

したがって、これまでの結果は「1次でも前縁が発散した」という証拠ではありません。

保存結果も次のようになります。

| 保存場・記録 | 再確認した事実 |
|---|---|
| `CONVERGENCE_VERDICT_warm.txt` | **NOT CONVERGED**。`rms_roUx` は低下0.1桁、`rms_roUy` は上昇 |
| `_warm_res_3000.h5` | 非有限なし。ただし最大静圧411 kPa、最大温度915 K。良好な収束場ではない |
| `res_0.h5` | SST開始時は正密度。壁 ω は約 `1.17e9–6.09e9 s⁻¹` |
| `residual_history.csv` | ω残差が初めて `1e10` を超えるのは **step 518**。600 step完走を安定判定に使えない |
| `res_1000.h5` | 前縁壁ノード0で **ρ = −6.061 kg/m³**。全域最大 k は140程度 |
| `res_nan_1762.h5` | 保存量5成分で各70節点が非有限。位置は**出口域**、壁ノードではない |

再実行した収束判定は、SST区間について：

```text
[last step 1761] -> DIVERGED (NaN/Inf)
rms_roe : NaN/Inf present
```

CSVの最初の非有限残差は `rms_roe` の **1761**、ログの検出・ダンプは **1762** です。

最終ダンプの非有限節点は `x=0.18004–0.19533 m, y=−0.42783–−0.29766 m`。有限値の最大値も `roK≈7.91e12`、`roOmega≈1.74e21` で、ご提示の指紋とは異なります。**前縁異常は確認できますが、最終 NaN を「鼻の数十節点だけ」とする説明は、このファイルには当てはまりません。**

**実装上の重要な欠陥：等温壁で `roe` が初期値に固定される**

コードの組み合わせはこうなっています。

- [等温壁の陰解法行](/home/sano/work/forge/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1023)は `dq_roe=0`。
- [更新](/home/sano/work/forge/solver_density_cuda/cuda_forge/update_d.cu:241)は `roe=roeN+dq_roe`。
- [壁温ピン](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:74)は残差評価中に `roe=ρe(Tw)` とするものの、その後の更新が古い `roeN` を戻す。
- [残差評価の順序](/home/sano/work/forge/solver_density_cuda/main.cpp:1350)は **EOS → 物性 → SST壁境界 → 壁温ピン**。物性と壁 ω は、壁温に戻す前の温度を使う。

実測で、**壁320節点すべての `roe` が `res_0` と `res_1000` でビット一致**しました。密度は変わっているため、正密度の壁節点でも `roe/[ρcvTw]` は **0.579–1.399**。保存量は等温条件を満たしていません。

前縁ノード0、`(x,y)=(0.001097,0.116540) m` は特に明瞭です。

| 量 | `res_1000.h5` |
|---|---:|
| 保存密度 | −6.061 |
| 保存 `roe` | 659599.5 |
| 表示温度 | 574.013 K |
| 保存 `vis_lam` | **0.00442923 Pa·s** |
| 574.013 Kに対応する粘性 | **0.0000293284 Pa·s** |
| 壁 ω | **8.11823e15 s⁻¹** |

EOSが密度を `roMin=1e-4` に戻した時の温度は、

\[
T_{\rm EOS}=\frac{659599.5}{10^{-4}(1004.5/1.4)}
          \simeq9.193\times10^6\ {\rm K}.
\]

この温度を[実装のSutherland式](/home/sano/work/forge/solver_density_cuda/cuda_forge/gasProperties_d.cu:58)に入れると、保存されている **0.00442923** を再現します。壁温から期待される値の **151倍**です。

さらに[壁 ω の式](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransBoundary_d.cu:42)に、実測距離 `2.08919 µm` を入れると、

\[
\omega_w=\frac{60\mu}{0.075\,\rho\,y_1^2}
         \simeq8.11823\times10^{15},
\]

も再現します。**異常な ω は、少なくともこの時点では「密度床＋誤った壁物性」の直接的な結果です。**

ただし、step 518近辺の場・勾配・流束が保存されていないため、**最初に負密度を作った寄与が再構成、陰解法更新、乱流輸送のどれかを一意には確定できません**。前縁異常から最終出口破綻までの因果も、現在の疎な保存場だけでは証明できません。

**SSTそのものについて**

[ransSource_d.cu](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransSource_d.cu:124)には KL の `SΩ`、`Pk≤10β*ρkω` の制限、制限後 `Pk` に基づく ω生成が実装されています。[渦粘性](/home/sano/work/forge/solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu:168)にも SST のせん断制限があります。壁 `k=0` も保存場で確認しました。

よどみ点の偽生成への対策として KL は合理的ですが、**今回「kの暴走が最初」とする証拠はありません**。壁 ω の `60ν/(β₁y₁²)` 自体は標準的な境界条件で、2 µmなら `10⁹` 級でも異常とは限りません。[NASA TMRのSST定義](https://tmbwg.github.io/turbmodels/sst.html)

`inlet_Pressure_dir` にも別の弱点があります。[CPG分岐](/home/sano/work/forge/solver_density_cuda/cuda_forge/boundaryCond_d.cu:1381)は `Pcell>Pt` で平方根が負になり、方向ベクトルのゼロ割ガードもありません。ただし今回の入口方向は非ゼロ、保存時点の入口圧力は319.5 kPa未満、非有限位置も出口なので、**今回の主因とする根拠はありません**。

**2. 推奨する起動と対策**

**優先順位は「等温壁状態の整合性修正 → 本当の1次で再起動 → SST → 制限付き2次」です。** 現行コードのまま、CFLだけで解決済みにしてはいけません。

実装修正の要件は、壁で常に `roe=ρe(Tw)` を満たし、その状態から EOS・物性・壁 ω を評価することです。陰解法で密度を更新するなら、静止等温壁の拘束は `δroe=e(Tw)δρ`。単にエネルギー増分をゼロにする現在の処理との整合を取る必要があります。修正後はまず「壁の保存量から算出した温度＝Tw」「μ＝μ(Tw)」を検証します。

以下は**未検証の保守的な起動案**です。

```yaml
space:
  convMethod: 0
  limiter: 2
  pRef: 250000.0

time:
  unsteady: 0
  timeIntegration: 11
  nStepInner: 10
  deltaT:
    control: 1
    dt: 1.0e-9
    cfl: 0.01
    cfl_pseudo: 0.01
    implicitRelax: 0.5
    blockDPLUR: 1
    lowMachPrecond: 0
    dt_min: 1.0e-12
    dt_max: 1.0
    detectNaN: 1
```

- まず `model: none`。現在のwarm場は未収束なので、引き継ぐ場合も「初期推定」として扱う。
- 全保存量の収束と壁状態の整合を確認してから SST。初期は同じCFL・次数を保持する。
- 正値性・残差を確認して `cfl_pseudo: 0.02 → 0.05`。固定600 stepで次段へ進めない。
- SSTの1次場を得てから `convMethod: 1, limiter: 2`。切替時はCFLを戻す。`implicitRelax: 0.7`、`nStepInner: 4` はその後の効率化候補。
- `lowMachPrecond: 0` を維持。これは壁の状態不整合を直す機能ではなく、既存検証にもSSTとの発散例がある。

根拠は [起動手順 §0・§2](/home/sano/work/forge/procedures/divergence-and-startup.md)、[推奨設定 §1.2・§2](/home/sano/work/forge/procedures/recommended-settings.md)、[更新ガード検証 §5](/home/sano/work/forge/plans/accepted/time_integration-update-positivity-guard.md)。一般レシピのCFL 0.5以上を、この未検証翼列へそのまま移すことは勧めません。

乱流設定はまず次を維持します。

```yaml
turbulence:
  model: sst
  scalarDiffusion: 1
  wallTreatmentSST: 0
  katoLaunder: 1
  dilatationCorrection: 0
  turbulentPrandtl: 0.9
  kInit: 57.4
  omegaInit: 1384.0
```

**`kInf` / `omegaInf` は存在しない設定キー**です。[読込実装](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:693)は `kInit` / `omegaInit`。ただし今回のHDF5には非ゼロの保存量があるので、これを今回の直接原因にはしません。既存 `roK/roOmega` がある場合、初期化キーだけでは上書きされない点にも注意が必要です。

`57.4` は `M₁=0.17, Tu=6.5%` からの `k=1.5(Tu U)²≈57.5` と整合します。`1384` は

\[
L=\frac{\sqrt{k}}{0.09^{1/4}\omega}\simeq10\ {\rm mm}
\]

という長さ尺度の仮定です。**Tuだけでは ω は決まりません**。安定化目的で ω を恣意的に増やす前に、この尺度の根拠を確認すべきです。今回の未収束入口速度は約188 m/sなので、同じkはTu約3.3%に相当します。

次の診断では `output.level: 2` と短い保存間隔を使い、特にstep 500付近の**更新前後のρ・roe、面再構成、密度残差、補正量、ω項別収支**を追う必要があります。NaN検出だけでは負密度を見逃します。

**3. メッシュについて**

再実行した品質判定は：

```text
VERDICT: SOFT-PASS (<0.1% outliers)
AR max = 822.3
skew max = 0.964
```

| 調査対象 | 結果 |
|---|---|
| 双対CV体積 | 全22442 CVで正。最小 `2.275e-10 m²` |
| 壁ω用の第一内部距離 | 全320壁節点で取得可能。`1.987–2.259 µm` |
| 周期節点 | 294対、一対一。平行移動後の最大誤差 `2.15e-8 m` |
| 周期対の保存値 | 調べたρ・P・ρk・ρωで差ゼロ |
| 壁と周期境界の共有節点 | なし |
| skew超過4セル | `x≈0.00458, y≈0.10221–0.10239 m`。問題の前縁ノード0から約14.6 mm |

したがって、**ωの距離ゼロ、双対体積ゼロ、周期∩壁角、周期節点不一致を主因とする証拠はありません**。

高skewの4セルは primal ID **0、1、23420、23440**。修正するなら、この位置の境界層外縁と三角形領域の接続を局所再分割・平滑化します。前縁だけを一律細分化する対策とは別です。

第一層厚は、最初の切り分けでは **2 µmを維持**します。4 µm化は壁ωを概ね1/4にするため診断には使えますが、今回の熱力学的不整合を隠す可能性があります。また、**2 µmだからy⁺≈1と確認されたわけではありません**。正常な場を得てから第一内部距離と接線壁応力で局所 `y₁⁺` を評価すべきです。

**4. forgeでは無理なのか、代替は何か**

**node median-dual＋低Re SSTで翼列を解くこと自体に、今回示された原理的な無理はありません。ただし現行の等温壁処理と今回のrunを、熱伝達の定量検証に使うのは無理です。** 壁温表示が正しくても、物性が151倍違う状態を通っています。

代替を一つに絞るなら、**同一C3Xメッシュ・同じ実測壁温分布によるSU2の低Re SST対照計算**を推奨します。[比較手順](/home/sano/work/forge/procedures/su2-cross-check.md)に沿い、乱流モデルの変種・入口乱流量も合わせる。壁関数やhwall緩和へ先に変えると、解きたい熱伝達問題そのものが変わり、forge固有の欠陥を切り分けにくくなります。

今回の提案は **plan未反映**です。read-only指定に従っています。また、[case READMEのrun一覧](/home/sano/work/forge/case/53.c3x_vane_cht/README.md:157)はまだ「無し」で、実在する `run_0001_measTw/` を反映していません。
tokens used
150,742
**「前縁よどみ点の SST 過大生成が真因」とは、今回の保存データからは言えません。確認できたのは、①起動設定の次数の誤認、②等温壁の保存量・物性更新の不整合、③前縁で既に負密度、その後に出口域で NaN、です。** 特に②は設定調整だけでは解消しない実装上の問題です。

ファイル変更・新規計算は行っていません。以下の数値はすべて `case/53.c3x_vane_cht/run_0001_measTw/` の保存場・ログを読み直した結果です。

**1. 実際に何が壊れているか**

まず、`convMethod: 1, limiter: 0` は**リミタなしの2次再構成**です。1次ではありません。

[再構成の分岐](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_common_d.cuh:171)は `convMethod=1` を2次へ送り、[limiter の処理](/home/sano/work/forge/solver_density_cuda/cuda_forge/limiter_d.cu:450)は `limiter=0` で係数を **1** にします。[setup_run.py](/home/sano/work/forge/case/53.c3x_vane_cht/tools/setup_run.py:47)は warm・sst_soft・sst_mid をすべてこの設定にしています。

したがって、これまでの結果は「1次でも前縁が発散した」という証拠ではありません。

保存結果も次のようになります。

| 保存場・記録 | 再確認した事実 |
|---|---|
| `CONVERGENCE_VERDICT_warm.txt` | **NOT CONVERGED**。`rms_roUx` は低下0.1桁、`rms_roUy` は上昇 |
| `_warm_res_3000.h5` | 非有限なし。ただし最大静圧411 kPa、最大温度915 K。良好な収束場ではない |
| `res_0.h5` | SST開始時は正密度。壁 ω は約 `1.17e9–6.09e9 s⁻¹` |
| `residual_history.csv` | ω残差が初めて `1e10` を超えるのは **step 518**。600 step完走を安定判定に使えない |
| `res_1000.h5` | 前縁壁ノード0で **ρ = −6.061 kg/m³**。全域最大 k は140程度 |
| `res_nan_1762.h5` | 保存量5成分で各70節点が非有限。位置は**出口域**、壁ノードではない |

再実行した収束判定は、SST区間について：

```text
[last step 1761] -> DIVERGED (NaN/Inf)
rms_roe : NaN/Inf present
```

CSVの最初の非有限残差は `rms_roe` の **1761**、ログの検出・ダンプは **1762** です。

最終ダンプの非有限節点は `x=0.18004–0.19533 m, y=−0.42783–−0.29766 m`。有限値の最大値も `roK≈7.91e12`、`roOmega≈1.74e21` で、ご提示の指紋とは異なります。**前縁異常は確認できますが、最終 NaN を「鼻の数十節点だけ」とする説明は、このファイルには当てはまりません。**

**実装上の重要な欠陥：等温壁で `roe` が初期値に固定される**

コードの組み合わせはこうなっています。

- [等温壁の陰解法行](/home/sano/work/forge/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1023)は `dq_roe=0`。
- [更新](/home/sano/work/forge/solver_density_cuda/cuda_forge/update_d.cu:241)は `roe=roeN+dq_roe`。
- [壁温ピン](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:74)は残差評価中に `roe=ρe(Tw)` とするものの、その後の更新が古い `roeN` を戻す。
- [残差評価の順序](/home/sano/work/forge/solver_density_cuda/main.cpp:1350)は **EOS → 物性 → SST壁境界 → 壁温ピン**。物性と壁 ω は、壁温に戻す前の温度を使う。

実測で、**壁320節点すべての `roe` が `res_0` と `res_1000` でビット一致**しました。密度は変わっているため、正密度の壁節点でも `roe/[ρcvTw]` は **0.579–1.399**。保存量は等温条件を満たしていません。

前縁ノード0、`(x,y)=(0.001097,0.116540) m` は特に明瞭です。

| 量 | `res_1000.h5` |
|---|---:|
| 保存密度 | −6.061 |
| 保存 `roe` | 659599.5 |
| 表示温度 | 574.013 K |
| 保存 `vis_lam` | **0.00442923 Pa·s** |
| 574.013 Kに対応する粘性 | **0.0000293284 Pa·s** |
| 壁 ω | **8.11823e15 s⁻¹** |

EOSが密度を `roMin=1e-4` に戻した時の温度は、

\[
T_{\rm EOS}=\frac{659599.5}{10^{-4}(1004.5/1.4)}
          \simeq9.193\times10^6\ {\rm K}.
\]

この温度を[実装のSutherland式](/home/sano/work/forge/solver_density_cuda/cuda_forge/gasProperties_d.cu:58)に入れると、保存されている **0.00442923** を再現します。壁温から期待される値の **151倍**です。

さらに[壁 ω の式](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransBoundary_d.cu:42)に、実測距離 `2.08919 µm` を入れると、

\[
\omega_w=\frac{60\mu}{0.075\,\rho\,y_1^2}
         \simeq8.11823\times10^{15},
\]

も再現します。**異常な ω は、少なくともこの時点では「密度床＋誤った壁物性」の直接的な結果です。**

ただし、step 518近辺の場・勾配・流束が保存されていないため、**最初に負密度を作った寄与が再構成、陰解法更新、乱流輸送のどれかを一意には確定できません**。前縁異常から最終出口破綻までの因果も、現在の疎な保存場だけでは証明できません。

**SSTそのものについて**

[ransSource_d.cu](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransSource_d.cu:124)には KL の `SΩ`、`Pk≤10β*ρkω` の制限、制限後 `Pk` に基づく ω生成が実装されています。[渦粘性](/home/sano/work/forge/solver_density_cuda/cuda_forge/turbulent_viscosity_d.cu:168)にも SST のせん断制限があります。壁 `k=0` も保存場で確認しました。

よどみ点の偽生成への対策として KL は合理的ですが、**今回「kの暴走が最初」とする証拠はありません**。壁 ω の `60ν/(β₁y₁²)` 自体は標準的な境界条件で、2 µmなら `10⁹` 級でも異常とは限りません。[NASA TMRのSST定義](https://tmbwg.github.io/turbmodels/sst.html)

`inlet_Pressure_dir` にも別の弱点があります。[CPG分岐](/home/sano/work/forge/solver_density_cuda/cuda_forge/boundaryCond_d.cu:1381)は `Pcell>Pt` で平方根が負になり、方向ベクトルのゼロ割ガードもありません。ただし今回の入口方向は非ゼロ、保存時点の入口圧力は319.5 kPa未満、非有限位置も出口なので、**今回の主因とする根拠はありません**。

**2. 推奨する起動と対策**

**優先順位は「等温壁状態の整合性修正 → 本当の1次で再起動 → SST → 制限付き2次」です。** 現行コードのまま、CFLだけで解決済みにしてはいけません。

実装修正の要件は、壁で常に `roe=ρe(Tw)` を満たし、その状態から EOS・物性・壁 ω を評価することです。陰解法で密度を更新するなら、静止等温壁の拘束は `δroe=e(Tw)δρ`。単にエネルギー増分をゼロにする現在の処理との整合を取る必要があります。修正後はまず「壁の保存量から算出した温度＝Tw」「μ＝μ(Tw)」を検証します。

以下は**未検証の保守的な起動案**です。

```yaml
space:
  convMethod: 0
  limiter: 2
  pRef: 250000.0

time:
  unsteady: 0
  timeIntegration: 11
  nStepInner: 10
  deltaT:
    control: 1
    dt: 1.0e-9
    cfl: 0.01
    cfl_pseudo: 0.01
    implicitRelax: 0.5
    blockDPLUR: 1
    lowMachPrecond: 0
    dt_min: 1.0e-12
    dt_max: 1.0
    detectNaN: 1
```

- まず `model: none`。現在のwarm場は未収束なので、引き継ぐ場合も「初期推定」として扱う。
- 全保存量の収束と壁状態の整合を確認してから SST。初期は同じCFL・次数を保持する。
- 正値性・残差を確認して `cfl_pseudo: 0.02 → 0.05`。固定600 stepで次段へ進めない。
- SSTの1次場を得てから `convMethod: 1, limiter: 2`。切替時はCFLを戻す。`implicitRelax: 0.7`、`nStepInner: 4` はその後の効率化候補。
- `lowMachPrecond: 0` を維持。これは壁の状態不整合を直す機能ではなく、既存検証にもSSTとの発散例がある。

根拠は [起動手順 §0・§2](/home/sano/work/forge/procedures/divergence-and-startup.md)、[推奨設定 §1.2・§2](/home/sano/work/forge/procedures/recommended-settings.md)、[更新ガード検証 §5](/home/sano/work/forge/plans/accepted/time_integration-update-positivity-guard.md)。一般レシピのCFL 0.5以上を、この未検証翼列へそのまま移すことは勧めません。

乱流設定はまず次を維持します。

```yaml
turbulence:
  model: sst
  scalarDiffusion: 1
  wallTreatmentSST: 0
  katoLaunder: 1
  dilatationCorrection: 0
  turbulentPrandtl: 0.9
  kInit: 57.4
  omegaInit: 1384.0
```

**`kInf` / `omegaInf` は存在しない設定キー**です。[読込実装](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:693)は `kInit` / `omegaInit`。ただし今回のHDF5には非ゼロの保存量があるので、これを今回の直接原因にはしません。既存 `roK/roOmega` がある場合、初期化キーだけでは上書きされない点にも注意が必要です。

`57.4` は `M₁=0.17, Tu=6.5%` からの `k=1.5(Tu U)²≈57.5` と整合します。`1384` は

\[
L=\frac{\sqrt{k}}{0.09^{1/4}\omega}\simeq10\ {\rm mm}
\]

という長さ尺度の仮定です。**Tuだけでは ω は決まりません**。安定化目的で ω を恣意的に増やす前に、この尺度の根拠を確認すべきです。今回の未収束入口速度は約188 m/sなので、同じkはTu約3.3%に相当します。

次の診断では `output.level: 2` と短い保存間隔を使い、特にstep 500付近の**更新前後のρ・roe、面再構成、密度残差、補正量、ω項別収支**を追う必要があります。NaN検出だけでは負密度を見逃します。

**3. メッシュについて**

再実行した品質判定は：

```text
VERDICT: SOFT-PASS (<0.1% outliers)
AR max = 822.3
skew max = 0.964
```

| 調査対象 | 結果 |
|---|---|
| 双対CV体積 | 全22442 CVで正。最小 `2.275e-10 m²` |
| 壁ω用の第一内部距離 | 全320壁節点で取得可能。`1.987–2.259 µm` |
| 周期節点 | 294対、一対一。平行移動後の最大誤差 `2.15e-8 m` |
| 周期対の保存値 | 調べたρ・P・ρk・ρωで差ゼロ |
| 壁と周期境界の共有節点 | なし |
| skew超過4セル | `x≈0.00458, y≈0.10221–0.10239 m`。問題の前縁ノード0から約14.6 mm |

したがって、**ωの距離ゼロ、双対体積ゼロ、周期∩壁角、周期節点不一致を主因とする証拠はありません**。

高skewの4セルは primal ID **0、1、23420、23440**。修正するなら、この位置の境界層外縁と三角形領域の接続を局所再分割・平滑化します。前縁だけを一律細分化する対策とは別です。

第一層厚は、最初の切り分けでは **2 µmを維持**します。4 µm化は壁ωを概ね1/4にするため診断には使えますが、今回の熱力学的不整合を隠す可能性があります。また、**2 µmだからy⁺≈1と確認されたわけではありません**。正常な場を得てから第一内部距離と接線壁応力で局所 `y₁⁺` を評価すべきです。

**4. forgeでは無理なのか、代替は何か**

**node median-dual＋低Re SSTで翼列を解くこと自体に、今回示された原理的な無理はありません。ただし現行の等温壁処理と今回のrunを、熱伝達の定量検証に使うのは無理です。** 壁温表示が正しくても、物性が151倍違う状態を通っています。

代替を一つに絞るなら、**同一C3Xメッシュ・同じ実測壁温分布によるSU2の低Re SST対照計算**を推奨します。[比較手順](/home/sano/work/forge/procedures/su2-cross-check.md)に沿い、乱流モデルの変種・入口乱流量も合わせる。壁関数やhwall緩和へ先に変えると、解きたい熱伝達問題そのものが変わり、forge固有の欠陥を切り分けにくくなります。

今回の提案は **plan未反映**です。read-only指定に従っています。また、[case READMEのrun一覧](/home/sano/work/forge/case/53.c3x_vane_cht/README.md:157)はまだ「無し」で、実在する `run_0001_measTw/` を反映していません。
