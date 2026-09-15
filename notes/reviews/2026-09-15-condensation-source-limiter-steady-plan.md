# codex レビュー: condensation-source-limiter-steady (plan)

- **plan**: [`plans/active/condensation-source-limiter-steady.md`](../../plans/active/condensation-source-limiter-steady.md)
- **stage**: `plan`
- **date**: 2026-09-15
- **commit**: `3a065ca8` (feature/sern-design)
- **codex**: effort `high`, 5.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M6/m1
- **focus**: §4.2 の更新クランプ化 (θ を残差でなく更新量に掛ける) が定常固定点を Δτ 非依存にしつつ起動安定性を保てるか、§4.3 平衡緩和形を触らない判断、§6 の判定ゲート
- **extra**: `case/44.vitiated_air_wt/README.md`, `solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

## 判定: **GO-with-changes**

残差から θ を外す方向は正しい。ただし、現在の計画では固定点の Δτ 非依存と起動安定性を保証できません。  
実装前に、蒸発経路、更新クランプの境界条件、平衡緩和形の説明、検証ゲートを修正してください。ファイルは変更していません。

## 指摘一覧

### 1. **Major — §4.2 の固定点保存は条件付き。現在のクランプ設計には更新停止の穴がある**

**根拠:** [対象 plan:57](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:57)、[scalarTransport_d.cu:285](/home/sano/work/forge/solver_density_cuda/cuda_forge/scalarTransport_d.cu:285)、[condensationTransport_d.cu:82](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationTransport_d.cu:82)。

更新を
\[
U^{k+1}=U^k+\theta_u B R(U^k)
\]
とすれば、**正の \(\theta_u\)、非特異な \(B\)、追加の射影なし**という条件では固定点は \(R=0\) のままです。しかし計画の式は `avail=0, Δg>0` で \(\theta_u=0\) となり、残差が残っていても４モーメントを全部停止できます。

また、現在の更新には既に次が入っています。

- スカラー更新内部の `max(updated, floor)`。
- 後続の `rog≤roY_w`／`rog≤0.99ρ` と各モーメント非負の硬クランプ。

したがって、**更新後の差分がゼロでも未加工の残差がゼロとは限りません**。４増分を同率に縮めても、先行する個別クリップで失われた増分方向は戻りません。

**対案:** 未クリップの候補増分を取得し、４モーメントをまとめて制限してから確定する。`avail=0`、負の候補、純物質の `0.99ρ` 上限、float32 の丸め停止を明示的に扱い、未加工残差と硬クランプの補正量を監視する。「\(\theta_u\to1\)」だけを固定点の証拠にしないこと。

### 2. **Major — 蒸発側を残すと、非平衡ソースの Δτ 依存は解消しない**

**根拠:** [condensationSource_d.cuh:266](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationSource_d.cuh:266)、[condensationSourceKernels_d.cuh:133](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:133)。

蒸発は `cond_evap_source` 内で
\[
\lambda=1+\dot r\,\Delta\tau/r_{30},\qquad
S_g=\rho g(\lambda^3-1)/\Delta\tau
\]
を作り、さらに半径半減・Δg・ΔT制限を掛けています。**制限が非作動でも**
\[
S_g=\rho g(3a+3a^2\Delta\tau+a^3\Delta\tau^2),\quad a=\dot r/r_{30}
\]
なので Δτ が残ります。`SQ2` も同様です。

蒸発分岐では `diagLim` が初期値１のまま返るため、`condLim≈1` はこの依存を検出できません。`condEvaporation` は既定で有効です。

**対案:** 定常の非平衡モデル全体を Δτ 非依存にするなら、蒸発も瞬間速度のソースと更新制限に分離する設計を追加する。液滴消滅処理も含め、**状態を固定して Δτ だけ変えても残差が変わらない単体試験**を必須にする。凝縮側だけの修正で「非平衡定常解が Δτ 非依存」と完了宣言してはいけません。

### 3. **Major — §4.3 は輸送を落としており、平衡緩和形を据え置く理由が誤っている**

**根拠:** [main.cpp:1205](/home/sano/work/forge/solver_density_cuda/main.cpp:1205)、[condensationSourceKernels_d.cuh:115](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:115)。

`eq=1` は輸送残差を保持したまま緩和ソースを加算します。定常条件は
\[
R_{\mathrm{transport}}
+V\alpha\theta\rho(g_{eq}-g)/\Delta\tau=0
\]
であり、一般には \(g=g_{eq}\) ではありません。

例えば、固定した上流値を持つ１セルの移流緩和モデルでは、
\[
g=\frac{\beta g_{\mathrm{in}}+kg_{eq}}{\beta+k},
\qquad k=\alpha\theta/\Delta\tau
\]
となり、固定点自体が Δτ に依存します。「θ は接近速度だけを変える」は輸送のない局所緩和に限った説明です。

**対案:** `eq=1` を今回変更しない判断は維持してよい。ただし理由を「旧モデルの互換性維持とスコープ分離」に変更し、**輸送との釣り合いに擬似時間依存がある既知の制約**として記録する。平衡計算の推奨を `eq=2` とする説明と、非平衡修正の達成範囲を分けてください。

### 4. **Major — dual-time の前提が成立していない。凝縮モーメントに物理時間項がない**

**根拠:** [update_d.cu:423](/home/sano/work/forge/solver_density_cuda/cuda_forge/update_d.cu:423)、[update_d.cu:485](/home/sano/work/forge/solver_density_cuda/cuda_forge/update_d.cu:485)、[main.cpp:1555](/home/sano/work/forge/solver_density_cuda/main.cpp:1555)。

確認した現行実装では、BDF物理時間項と時間レベルのシフトは平均流５変数と `roK/roOmega` が対象です。`rog`・`roQ0..2` は含まれません。

一方、凝縮モーメントはサブ反復ごとに `N/M` を現在値へコピーし、定常と同じ point-implicit 更新を行っています。その対角にも物理時間項がありません。つまり、**小さい物理 Δt を指定しただけでは、モーメントがその物理時間で積分される保証がありません**。

これは既存の構造的問題ですが、§2・§6の「物理 Δt が小さいので実害なし」「１本変わらなければ確認完了」を無効にします。

**対案:** 新経路の初回適用を定常に限定する。凝縮 dual-time は独立した前提課題として、物理履歴・BDF残差・対角を整備してから、物理 Δt 半減とサブ反復数変更で検証する。今回の新旧一致を時間精度の保証に使わないこと。

### 5. **Major — §6 は未収束の参照場を固定点と扱い、VERDICT の添付を合否条件にしていない**

**実測:** 以下を再実行しました。

- `check_convergence.py`：３ run とも **`NOT CONVERGED (stalled/plateau)`**
- `check_quasisteady.py --quantity machmax,pmax`：３ run とも **`STEADY`**
- 保存済み全スナップショットの `VALUE/*`：NaN/Inf **０**

根拠 run は次の３つです。

- `case/44.vitiated_air_wt/run_0127_va3_M4.19_Lc8_noneq_inletTt_merged/`
- `case/44.vitiated_air_wt/run_0130_va3_M4.19_Lc8_noneq_inletTt_cfl1/`
- `case/44.vitiated_air_wt/run_0131_va3_M4.19_Lc8_noneq_inletTt_cfl05/`

| run | 最終 `condLim_0` 最小値 | 出口 g〔質量流束平均〕 | `rms_roQ0_0` 低下 |
|---|---:|---:|---:|
| `0127` | 0.241140 | 0.436890 % | 2.5桁 |
| `0130` | 0.496454 | 0.573848 % | 2.0桁 |
| `0131` | **0.997275** | 0.584432 % | 2.0桁 |

出口断面は既存後処理と同じ \(x/r_t=22.2561\)。出口 g/M/T と g最大値の時系列を、`check_quasisteady.py` の `classify_series` に直接渡した追加確認では、変動・ドリフト閾値各0.2%で **`STEADY`** でした。これはCFL依存の数値差を支持しますが、残差収束の証明にはなりません。

さらに `run_0131` は最終場の１点で θ<1 です。「θ≡1」「既に新固定点」は撤回すべきです。

**対案:** [plan:107](/home/sano/work/forge/plans/active/condensation-source-limiter-steady.md:107) のゲートを以下に変更する。

- 比較する両 run の `check_convergence` **PASSを必須**とし、未達は未完了とする。
- 凝縮固有量をCSV時系列化し、`--series-csv` で判定する。既定のドリフト5%・変動10%は、場比較1%のゲートより緩すぎる。
- gの相対差について、ノルム・分母・凝縮領域の和集合を定義する。出口平均の重み、onsetの流線・閾値、T/Mの許容差も明記する。
- 「24000 stepで終了」ではなく、収束・定常性ゲートで終了する。
- `run_0131` は比較参考とし、無制限ソースの固定点の正解データとは扱わない。

### 6. **Major — 更新位置・状態の基準・試験対象が不足し、起動安定性を判定できない**

**根拠:** [main.cpp:1348](/home/sano/work/forge/solver_density_cuda/main.cpp:1348)、[main.cpp:1400](/home/sano/work/forge/solver_density_cuda/main.cpp:1400)、[scalarTransport_d.cu:250](/home/sano/work/forge/solver_density_cuda/cuda_forge/scalarTransport_d.cu:250)、[verification/README.md:40](/home/sano/work/forge/procedures/verification/README.md:40)。

- 凝縮はNSのblock-DPLUR sweep内部ではなく、**流れ・化学種更新後の独立したpoint-implicit更新**です。「陰的 sweep 後」では呼び出し粒度が曖昧です。
- 密度が変わると、実際の質量分率変化は  
  \[
  g^{new}-g^{old}
  =\frac{\Delta(\rho g)-g^{old}\Delta\rho}{\rho^{new}}
  \]
  です。計画の `Δ(ρg)/ρ` を使うなら、「更新済み流れを固定したモーメント修正による変化」と定義する必要があります。
- `L/cv` はEOS温度変化の近似です。現行EOSは `cv + g(Rw−dL/dT)` を使います。[condensationEOS_d.cuh:280](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:280)
- RK4は未制限残差の累積バッファを持ちます。段の状態だけ制限しても、最終段で累積増分が再投入されます。RK3の段基準とも区別が必要です。
- §6-3は陰解法の軸対称case/44であり、RKの呼び出し粒度を決める試験になっていません。
- node/cell双方の試験指定がありません。周期モーメント輸送には既存課題F-cf1もあります。[condensation-followups.md:53](/home/sano/work/forge/plans/active/condensation-followups.md:53)

**対案:** §5の実装前に、定常の候補作成→制限→確定→primitive同期の順序と、ρ・組成・温度の評価時点を固定する。単体試験は提案の３項目に加え、**輸送とソースが非ゼロで釣り合う状態、蒸気枯渇、負増分、密度・組成変化、float32極小増分**を含める。

回帰は少なくともnode/cell、`condFloat=0/1`とdouble退避を明示する。RK・周期を保証範囲に含めるなら専用試験を追加し、未検証のまま既定適用しないこと。

### 7. **Minor — `condEquilibrium=2` は設定の既定値ではない**

**根拠:** [solverConfig.hpp:436](/home/sano/work/forge/solver_density_cuda/input/solverConfig.hpp:436)、[solverConfig.cpp:721](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:721)。

宣言と設定読込の既定値はともに **0＝非平衡**です。

**対案:** §4.3の「`condEquilibrium 2`, 既定」を「平衡凝縮を選ぶ場合の推奨方式」に直す。

## 推奨

**非平衡・定常に絞り、凝縮／蒸発の残差をΔτ非依存にして、未クリップの４モーメント候補増分をまとめて制限する方針を推奨します。** 既存の分離解法を利用でき、密結合ブロック化より変更を限定できます。

実装前の優先順位は次のとおりです。

1. 蒸発を含む残差のΔτ非依存条件と、クランプの境界・停止条件を定義する。
2. `eq=1`据え置きの説明を訂正し、dual-time・RKの保証を初回スコープから外す。
3. 更新位置・評価状態を固定し、収束PASS／凝縮量STEADY／node・cell回帰を完了ゲートにする。

`plans/README.md` と `accepted/` の確認では、この修正が完了済みの計画はありません。F-cf7からの独立化は妥当です。上記はレビュー提案であり、指定どおり**plan未反映**です。

指摘数: Critical 0 / Major 6 / Minor 1
