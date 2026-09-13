# codex レビュー: condensation-float-speedup (plan)

- **plan**: [`plans/active/condensation-float-speedup.md`](../../plans/active/condensation-float-speedup.md)
- **stage**: `plan`
- **date**: 2026-09-13
- **commit**: `cfcda6d1` (feature/perf-3d-speedup)
- **codex**: effort `high`, 8.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M7/m2
- **focus**: 凝縮経路 (相変化ソース・二相 EOS 反転・SLAU 面潜熱) の float 化計画。特に (1) 対数空間 CNT の float 範囲・丸めの安全性、(2) dry セル早期退出の等価性 (S<=1, g=0, Q0=0 で本当に全ソース・src_jac が 0 か)、(3) 二相ハイブリッド反転の設計、(4) 判定基準 §4.3 と回帰ケース §6 の妥当性 (凝縮 onset・h0 保存)、(5) スコープ外に残す double 経路の妥当性、を見てください。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

高速化の対象同定は妥当で、既存の dry 高速化計画とは重複していません。ただし、物性の単純な float 化、CNT の範囲保証、回帰基準には修正が必要です。**現状の §4・§6 のまま実装へ進むことは推奨しません。**

1. **Major — 物性係数の float コピーでは、要求精度を満たせない**

   根拠は [condensationProperties_d.cuh:114](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:114) の N2 飽和圧と、[同ファイル:231](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:231) の H2O 潜熱評価です。いずれも大きな項の相殺を含みます。特に液相 NASA-9 の値を差分して `cpl` を作る処理は誤差を増幅します。

   現行式を NumPy の float32 演算で再評価した結果は次のとおりです。これは **CUDA 実測ではなく、演算精度を落とした場合の数値実験**です。

   | 評価対象 | double | float32 | 相対差 |
   |---|---:|---:|---:|
   | H2O `L(200 K)` | 2,674,398.6 J/kg | 2,699,105.0 J/kg | `9.24e-3` |
   | N2 `psat(100 K)` | 776,276.8 Pa | 776,053.8 Pa | `2.87e-4` |

   さらに、**係数だけを float に丸め、演算を double に保っても**、H2O `L(300 K)` の相対差は `1.71e-4`。FMA の有無だけでは解決しません。§4.3 の `2e-6` と、§4.2 の「潜熱誤差は約 `1e-6`」という前提は成立していません。

   **対案:** 単体試験を実装順の先頭へ移し、定数 `h_l(273.15)`・`cpl` の double 前計算、変数シフト等による安定な式評価を先に設計してください。許容値を満たせない物性評価は、当面 double に残すべきです。

2. **Major — 対数 CNT は必要だが、`J≤1e35` を全経路の範囲保証に使えない**

   [condensationSource_d.cuh:92](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSource_d.cuh:92) は無制限の `J` を返します。`Jmax` を適用するのは [condensationSource_d.cu:187](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSource_d.cu:187) の本体だけで、[同ファイル:224](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSource_d.cu:224) のヤコビアン用 `cond_source_vector` には適用されません。

   現行式の double 再評価では、計画の試験範囲内である **N2、`T=60.2 K`、`S=100`、`condKantrowitz=0`** で、
   `lnJ≈93.2733`、`J≈3.22e40`。float 最大値の対数 `88.7228` を超えます。下限 `-80` だけでは安全になりません。

   また、最終的な `lnJ` の ULP から式全体の誤差を見積もる説明は不適切です。飽和圧・表面張力・障壁項の誤差と相殺を伝播させる必要があります。`∂lnSg/∂T≥0.05` も、成長停止・律速・蒸発消滅の分岐を含む全域の保証ではありません。

   **対案:** 本体の上限付き評価と、上限のない摂動評価を明確に分け、指数化前の上限処理、対数での積評価、必要箇所の double 退避を設計してください。摂動側にも単純に `Jmax` を掛けると既存ヤコビアンを変更するため、別の設計判断になります。`src_jac` 自体の有限性・符号・double 差分との比較を試験対象に追加すべきです。

3. **Major — 二相ハイブリッド反転に、失敗判定と保存量保護がない**

   [condensationEOS_d.cuh:22](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:22) の TP carrier 反転は、反復上限後も成功フラグなしで温度を返します。[dependentVariables_d.cu:197](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/dependentVariables_d.cu:197) は、その結果で `roe` を再構成します。

   したがって「float 12 回＋double 1〜3 回」「`|ΔT|<1e-3+1e-6T`」だけでは、`1e-8·T` の精度やエネルギー保存を保証できません。物性接続点や温度上限・下限では、滑らかな Newton の二次収束を前提にできません。CPG で既に同種の問題を修正した経緯が [condensationEOS_d.cuh:94](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationEOS_d.cuh:94) にあります。

   **対案:** 二相エネルギー残差による成功条件、非有限値検出、未達時の double 継続／括弧付き反転、失敗時に `roe` を改変しない経路を明記してください。`thermoFloat=0`・datum 無効時の扱いも固定し、温度精度に加えて保存量の反復往復ドリフトを検証すべきです。

4. **Major — `condL` キャッシュは、現状の配置では評価時点が変わる**

   [main.cpp:1145](/home/sano/work/forge-perf/solver_density_cuda/main.cpp:1145) の `dependentVariables` 後に境界処理が入り、[nodeWallDirichlet_d.cu:79](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:79) は実セルである壁ノードの `T` を変更します。一方、現行 SLAU は [convectiveFlux_slau_d.inc.cuh:340](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:340) で、その時点の `T_cell` を読みます。

   `dependentVariables` だけで `condL` を生成すると、こうしたセルでは旧温度の潜熱を使います。「セル温度評価だから等価」は、更新順まで一致して初めて成立します。

   **対案:** 最後の温度・凝縮状態更新後にキャッシュを確定するか、状態を変更する境界処理で更新してください。配列登録・初期化を含め、node/cell、ghost、周期相手セルの参照範囲を明記し、旧面評価との凍結状態比較を追加すべきです。

5. **Major — `condFloat: 0` の互換性設計と、ビット一致条件が整合していない**

   計画はソース kernel の切替を明記していますが、SLAU、clamp、移流融合、温度反転まで旧経路に戻す条件が不明です。特に [dependentVariables_d.cu:89](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/dependentVariables_d.cu:89) の `useHybrid` は、反転だけでなく組成 `Y` の計算精度も切り替えます。単に `condensation==0` を外すと、double に残す平衡反転への入力も変わります。

   また、指定された `case/16.nozzle_wys/run_0456_perf_regress_node2d_cond/` の `base_r1` と `base_r2` を実際に比較すると、`res_300.h5` の `ro` は **29,163 要素が非一致**、正規化最大差は `9.55e-7`。同一バイナリでも全場ビット一致は成立していません。

   **対案:** `condFloat` と各物理モードの分岐表を作り、旧経路を残す範囲を全変更に適用してください。ビット一致は決定的な単体・凍結状態評価で確認し、反復 run はノイズ床で判定すべきです。

   重点の dry 退出については、**非平衡の通常物性域ではソース・ヤコビアンがゼロになる判断は妥当**です。ただし [condensationSource_d.cu:44](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSource_d.cu:44) の初期化を維持し、輸送残差を消さないこと。[同ファイル:113](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/condensationSource_d.cu:113) の `condEquilibrium=2` に必要な残差凍結を、早期退出で飛ばしてはいけません。

6. **Major — 短期回帰から onset・`h0` の物理的非劣化を判定する手順が不足している**

   指定 run に判定ツールを実行した結果は以下です。

   | run（リポジトリ相対パス） | `check_convergence.py` | `check_quasisteady.py` |
   |---|---|---|
   | `case/16.nozzle_wys/run_0456_perf_regress_node2d_cond/base_r1/` | `NOT CONVERGED (stalled/plateau)` | `TRANSIENT-UNSETTLED`：2枚 |
   | `case/34.arthur_n2_nozzle/run_0100_merge_regress_cell_air/merged/` | `NOT CONVERGED (stalled/plateau)` | 今回未実行 |
   | `case/44.vitiated_air_wt/run_0201_perf_regress_node_axisym_euler_cond/merged/` | `NOT CONVERGED (stalled/plateau)` | `TRANSIENT-UNSETTLED`：2枚 |

   短い restart 区間の判定だけで元の場を否定するものではありませんが、これらを収束解の保証には使えません。

   さらに、[check_quasisteady.py:280](/home/sano/work/forge-perf/solver_density_cuda/tools/check_quasisteady.py:280) の標準量には onset・`h0` がありません。統合は [condensation-followups.md:52](/home/sano/work/forge-perf/plans/active/condensation-followups.md:52) の未完了項目です。`h0±0.1%` も分母が未定義で、[compare_condfix.py:69](/home/sano/work/forge-perf/case/16.nozzle_wys/compare_condfix.py:69) は sensible datum のため相対値を避けています。

   **対案:** 短期回帰と物理検証を分け、後者には元 run を含む残差履歴、報告量の十分な時系列、対象量そのものの VERDICT を要求してください。`h0` は `VALUE/h0` に基づく J/kg の絶対偏差と許容増分を定義し、onset とともに場比較の「ノイズまたは truth」救済から独立した必須条件にすべきです。

7. **Major — 単体範囲と回帰ケースの実体が、変更対象を覆っていない**

   `case/34.arthur_n2_nozzle/run_0100_merge_regress_cell_air/merged/res_12000.h5` の保存場を読んだところ、**`T_min=39.1903 K`、`g_max=0.0603974`**。計画の N2 `45–120 K`、`g≤0.05` の外です。`run_0456/.../base_r1/res_300.h5` の `condS_0` 最大値も **170.635** で、`S≤100` を超えます。いずれも保存場の範囲監査値であり、定常値の主張ではありません。

   また [case/34 の README:78](/home/sano/work/forge-perf/case/34.arthur_n2_nozzle/README.md:78) によれば、`run_0100–0102` の実体は **cell 空気／node 空気／cell dry**。計画が記す pure N2 の cell/node は含まれていません。host 単体だけでは、device の物性・`src_jac`・分岐の実装を直接検証できません。

   **対案:** 実測範囲、物性接続点、`S=1`、`g≈Yw`、消滅閾値、`Jmax` 前後を含めて試験を拡張してください。pure N2 の両離散化、Feder 2/3、Gyarmathy、蒸発の device 試験を具体化し、移流融合には周期境界を含む4モーメントの保存性試験を加えるべきです。各 run のメッシュ品質・IC 出所・段階起動条件も §6 に固定してください。

8. **Minor — 速度目標は dry 起動区間に限定され、実測根拠を今回追試できない**

   [plan:47](/home/sano/work/forge-perf/plans/active/condensation-float-speedup.md:47) の計測対象はモーメントゼロからの短期 run です。指定された `run_0420`・`run_0421` はこの workspace に存在せず、`89.5 ms/step` と kernel 内訳は今回独立検証できていません。

   **対案:** バイナリ・入力ハッシュ付き生ログと nsys 結果を参照可能にし、発達した凝縮場からの速度測定も追加してください。そこで湿潤セル比率・反転回数・fallback 頻度を記録すれば、CPG・平衡・二温度・音速を double に残す費用対効果も判断できます。これらを double に残す方向自体は妥当です。

9. **Minor — 検証ツールの仕様が計画の記述と違う**

   [perf_regress.py:72](/home/sano/work/forge-perf/solver_density_cuda/tools/perf_regress.py:72) の除外リストに `condTheta`・`condLim` はなく、現状では B 判定になります。また、`DFMA=0` だけでは `DADD`・`DMUL` 等の FP64 混入を排除できません。

   **対案:** 判定対象をツールと同期させてください。`condLim` は律速規約を検証する別ゲートとして残し、命令確認は kernel ごとの FP64 命令群と呼出先まで対象にすべきです。

**推奨は、物性の安定評価を先に確立する選択的混合精度化です。** 実装前の修正順は、①物性・CNT の範囲／誤差設計、②反転の成功判定・互換分岐・キャッシュ更新順、③実測範囲を覆う device 試験、④定常性と保存性の検証条件、⑤発達した凝縮場での性能条件、としてください。モーメント融合には既存の [scalarTransportResidualMulti_d](/home/sano/work/forge-perf/solver_density_cuda/cuda_forge/scalarTransport_d.cu:344) を再利用する方針が妥当です。

ファイル変更禁止の依頼に従い、**plan 未反映**です。

指摘数: Critical 0 / Major 7 / Minor 2
