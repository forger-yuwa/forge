# codex レビュー: condensation-h2o-latent-supercritical (plan)

- **plan**: [`plans/active/condensation-h2o-latent-supercritical.md`](../../plans/active/condensation-h2o-latent-supercritical.md)
- **stage**: `plan`
- **date**: 2026-09-14
- **commit**: `580dd8be` (feature/perf-3d-speedup)
- **codex**: effort `high`, 6.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

課題の同定と Watson 外挿の採用方針は妥当です。ただし、`h2o_latent` だけの変更では、ゼロ除算・float 表の大誤差・温度反転の失敗を持ち込みます。  
以下の Major を計画に反映するまで、実装着手は勧めません。

現行の液相エンタルピーのクランプと潜熱の上下限制限は、[condensationProperties_d.cuh:283](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:283) で確認しました。アンカー値も **2,269,078.998 J/kg** と再計算でき、目的の説明は正確です。[F-cf5](/home/sano/work/forge-perf/plans/active/condensation-followups.md:57) は未解決で、accepted の float 化計画も旧物性を表にしたものです。解決済み課題との重複ではありません。

以下の変更後の数値は、**提案式を既存コードの計算手順に代入した Python／NumPy による再現計算**です。変更後 CUDA バイナリの実測ではありません。既存 run は実ファイルを読み、指定ツールで検査しました。

1. **Major — `L=0` は、既存の蒸発式の定義域を破る**

   根拠：[condensationSource_d.cuh:247](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSource_d.cuh:247) の Gyarmathy 蒸発式は `L*L` で除算します。`T>=Tc` の停止条件はありません。

   H2O、`growthModel=1`、`T=650 K`、`p_v=p_gas=100 kPa`、`r=100 nm`、Kelvin 無効で再現すると、次になります。

   | 潜熱 | `dr/dt` |
   |---|---:|
   | 現行 | −0.0136006 m/s |
   | 提案式 | **−Inf** |

   蒸発量は後段の制限で有限になり得ますが、`drdt` はそのまま診断配列に書かれます（[condensationSourceKernels_d.cuh:140](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSourceKernels_d.cuh:140)）。したがって「物性関数だけ直せばよい」は成立しません。

   **対案：** `L=0` を許容する変更と同時に、成長・蒸発側の定義域を明示し、臨界点以上では `1/L²` の式を評価しない分岐を設計する。域外の湿潤状態をどう扱うか、保存量を含めて規約化し、`condFloat=0/1` 双方で有限性を試験してください。

2. **Major — §4.6 の表精度評価は、最も危険な使用区間を除外している**

   根拠：[condensationTables_d.cuh:80](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationTables_d.cuh:80) の格子では、臨界点 `647.096 K` が **646.90–647.15 K の区間内部**に入ります。同ファイルの片側差分・Hermite 係数・float 評価を再現した結果です。

   | T [K] | Watson [kJ/kg] | float 表 [kJ/kg] | 相対誤差 |
   |---|---:|---:|---:|
   | 646.95 | 129.417 | 120.715 | −6.72% |
   | 647.00 | 110.357 | 83.708 | **−24.15%** |

   `647.0 K` は `cond_tab_wet_ok` の許容範囲内で、[面潜熱](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSourceF_d.cuh:253)にも使われます。既存の表試験の高温側許容値 `1e-4` にも明確に違反します（[test_cond_float.cpp:64](/home/sano/work/forge-perf/solver_density_cuda/tests/unit/test_cond_float.cpp:64)）。

   また、§4.6 の「表面張力も同じ性質」は数学的に誤りです。`σ∝τ^1.256` の一階微分はゼロへ向かいますが、Watson の一階微分は発散します。

   **対案：** 今回は既存の double 退避を利用し、臨界点付近を表の使用域から外す方針を推奨します。退避温度は誤差試験で決定し、面・Newton・ソース摂動点を含めて統一する。試験範囲を `646 K` で打ち切らず、切替点両側と `Tc` まで含めてください。

3. **Major — 未記載の利用先 `cond_Tsat` が壊れる**

   根拠：[condensationProperties_d.cuh:353](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:353) は、Murphy–Koop の飽和圧を反転するとき、勾配を `L/(R*T*T)` で近似しています。潜熱と飽和圧は、実装上も完全に独立ではありません。

   `p_v=h2o_psat(640 K)`、初期推定 `500 K` で同じ25反復を再現すると、

   - 現行：`639.999995 K`
   - 提案式：**`471.355064 K`**、`ln(psat/p_v)=−2.775736`

   となります。失敗判定を返さず、この温度が診断値として採用されます。一方、[float 版](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSourceF_d.cuh:230) は飽和圧表自身の微分を使うため、両経路の乖離も拡大します。

   **対案：** double 版も、実際の `ln(psat)` の微分と括弧付き反転へ変更する。高温を含む往復試験 `Tsat(psat(T))` を、根から離れた初期推定値で追加してください。

4. **Major — `G` の単調性だけでは、現行温度反転の成功を保証できない**

   根拠：[condensationEOS_d.cuh:124](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:124) は固定幅 `±0.1 K` で微分し、129行では温度の括弧幅だけで反復を終了します。Watson の急勾配下では、温度幅が小さくてもエネルギー残差が残ります。

   `cv=1393.5`、`R=461.5`、`g=0.9`、初期推定 `300 K`、真の温度 `647.095 K` から生成したエネルギーを反転すると、**残差約1.095 J/kg、`ok=false`** でした。真の温度 `647.0959 K` では約3.201 J/kgです。

   保存量保護はありますが、「全域で堅牢性が改善する」という主張の根拠にはなりません。また、`T>Tc` では `L'=0` なので、§4.5 の厳密不等式 `G'>a` も成立しません。

   **対案：** 停止条件をエネルギー残差と整合させ、急勾配下で温度幅だけを理由に終了しないよう修正する。CPG と TP hybrid の実際の反転関数について、`373.15 K`、表退避点、`Tc` の両側で往復誤差・成功フラグ・失敗時の `roe` 保護を検証してください。

5. **Major — §6 は低温回帰にはなるが、変更部分の統合検証になっていない**

   根拠：[計画:108](/home/sano/work/forge-perf/plans/active/condensation-h2o-latent-supercritical.md:108) の唯一の場試験は、新しい枝を使わない低温 run です。また、共有コードには node／cell 双方の検証が必要です（[verification/README.md:39](/home/sano/work/forge-perf/procedures/verification/README.md:39)）。

   `case/16.nozzle_wys/run_0462_cond3d_local_spinup/` の実ファイルでは、`res_3000.h5` の温度は **201.9272–295.386 K、373.15 K 超は0/2,370,000**。保存された2枚の `VALUE/*` に NaN/Inf はありません。ただし、実行した判定は次のとおりです。

   ```text
   check_convergence.py:
   NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)

   check_quasisteady.py --quantity pmax:
   TRANSIENT-UNSETTLED
   only 2 snapshot(s) (<4)
   ```

   run の索引は [case/16 の README](/home/sano/work/forge-perf/case/16.nozzle_wys/README.md) にあります。この run は未収束状態の回帰対照には使えますが、高温モデルの妥当性を立証しません。

   **対案：** 低温のビット回帰を残し、高温湿潤状態の device 試験を追加する。`condFloat=0/1`、成長・蒸発、ソース Jacobian、面潜熱、音速を対象にする。場の回帰は node／cell 双方について、基準バイナリ・同一入力・比較する配列とステップ数を固定し、収束・準定常の判定と回帰差分の判定を明記してください。

6. **Minor — 蒸気表との精度主張が追試できない**

   根拠：[計画:50](/home/sano/work/forge-perf/plans/active/condensation-h2o-latent-supercritical.md:50) は相対誤差と指数探索結果だけを示し、参照潜熱値、出典、探索刻み、計算手順を示していません。

   **対案：** [IAPWS の飽和物性リリース](https://iapws.org/technical-guidance/release/Supp-sat)など、採用した一次資料・版・5点の参照値を固定し、誤差と RMS を再計算できる試験にする。飽和潜熱の近似精度と、既存の理想気体ベース二相 EOS の高温精度は区別してください。

7. **Minor — 実装ステップがリポジトリの文書更新順序と逆**

   根拠：[計画:89](/home/sano/work/forge-perf/plans/active/condensation-h2o-latent-supercritical.md:89) は実装・単体試験の後に `methods` 更新を置いていますが、[AGENTS.md:117](/home/sano/work/forge-perf/AGENTS.md:117) は現在仕様文書と計画を揃えてから実装する順序です。

   **対案：** 今回の採否、適用域、退避・反転・ゼロ潜熱の規約を先に `methods` と plan に反映し、影響ファイル一覧と残作業表を更新してください。

**推奨は、Watson 外挿を採用し、利用先の修正まで含む計画へ改訂することです。** 実装前の優先順は、①ゼロ潜熱の定義域処理、②臨界近傍の表退避、③`Tsat` 反転、④エネルギー反転、⑤検証計画の具体化です。373.15 K での傾きの不連続を許容する判断自体は支持します。

ファイルは変更していません。以上の指摘・推奨は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 2
