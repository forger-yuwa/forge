# codex レビュー: tooling-nozzle-sern-chain (plan)

- **plan**: [`plans/active/tooling-nozzle-sern-chain.md`](../../plans/active/tooling-nozzle-sern-chain.md)
- **stage**: `plan`
- **date**: 2026-09-09
- **commit**: `1a9d75c2` (feature/sern-design)
- **codex**: effort `high`, 8.8 min, rc=0
- **判定**: **NO-GO**, 指摘 C2/M7/m0
- **focus**: §5.1 残作業表の最優先 3 件 (作動点のサイクル値化・外部流ブロック・3D SST 後縁 3 重点) の方針と順序が妥当かに重点を置く
- **extra**: `case/46.sern_design/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
設計チェーンの目的は妥当ですが、評価器が発散・物性誤差・力の集計範囲の違いを設計性能として取り込む状態です。  
最優先は MOO 再取得ではなく、評価条件と判定の修復です。「3D SST 解決済み」という扱いも撤回すべきです。

1. **Critical — §4.13 の発散結果採用は撤回すべき。判定コードには NaN を `STEADY` にする欠陥もある。**

   **根拠:** [driver_sern.py:118](/home/sano/work/forge/design/forge_design/opt/driver_sern.py:118) は、力係数が `STEADY` なら `rc != 0` でも採用します。標準レシピでこの条件を満たすと再試行さえ省略します。さらに [sern_forces.py:108](/home/sano/work/forge/design/forge_design/metrics/sern_forces.py:108) を実行したところ、`steadiness([1,1,1,NaN])` は **`STEADY`、平均・傾き・変動幅すべて NaN** を返しました。

   判定対象も不一致です。最適化には `C_T_with_shear` を使う一方、定常判定は圧力だけの `C_T,C_L,C_M` です。`pareto.json` への書き出しでは `degraded` が脱落します。[driver_sern.py:139](/home/sano/work/forge/design/forge_design/opt/driver_sern.py:139)、[同:204](/home/sano/work/forge/design/forge_design/opt/driver_sern.py:204)

   **対案:** 発散・非有限値を含む評価はサロゲート学習と Pareto から除外する。保存場の健全性、全残差の収束、**実際に最適化する量**の定常性を独立した必須ゲートにする。壁ピン残差が判定を妨げるなら、その診断を修正する。失敗を物理的な `INFEASIBLE` と混同しない。既存の [accepted MOO 計画:33](/home/sano/work/forge/plans/accepted/tooling-nozzle-moo-loop.md:33) も、ゲート不合格を学習対象外としています。

2. **Critical — 作動点のサイクル値化は未完了。現在の「生産構成」は定動圧経路を再現していない。**

   **根拠:** [runner_sern.py:80](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:80) は排気と外気に同じ `gamma,R` を使います。現行 YAML をこの関数に通して再計算すると、次の値になります。

   | 作動点 | 外部流動圧 \(½ρu²\) | 計画の 71.85 kPa に対する差 |
   |---|---:|---:|
   | `m6_on` | 60.699 kPa | −15.5% |
   | `m10_on` | 62.909 kPa | −12.4% |
   | `m4_off` | 71.776 kPa | −0.1% |

   §8-7 の「M6 は誤差 0.1%」は音速・速度についての話で、密度は空気比 **−15.7%**。外力・せん断層の評価を正当化しません。

   熱力学も凍結 semi-perfect ではありません。[tmx_operating_points.py:60](/home/sano/work/forge/case/46.sern_design/cea/tmx_operating_points.py:60) は平衡 `GAMMAs` から CPG の `cp` を作っています。M6 の生成値は約 2202 J/kg/K、CEA 出力の平衡 `Cp` は 2551 J/kg/K です。これは音速に合わせた代用 CPG であり、凍結組成の熱力学を検証したことにはなりません。

   **対案:** **外気＝空気、排気＝CEA 組成を凍結した TP** に統一し、入口状態・エネルギー・理想推力の正規化まで同じ物性で計算する。既存の [多成分 TP 実装](/home/sano/work/forge/plans/accepted/thermophysics-multicomponent-tpgas.md:31) を再利用でき、有限速度化学の実装を待つ必要はありません。この修正より MOO 再取得を先に置くべきではありません。

3. **Major — §5.1-3 の「3D SST 解決済み」「三重点が真因」は証拠を超えている。**

   **根拠:** 今回 `check_convergence.py` を再実行しました。パスはすべて `case/46.sern_design/` 配下です。

   | run | 再判定 |
   |---|---|
   | `run_0082_3d_sst_lsw08/` | **NOT CONVERGED — stalled/plateau** |
   | `run_0083_3d_sst_cycle_m6on/` | **DIVERGED — NaN/Inf** |
   | `run_0084_3d_sst_cycle_pmin60/` | **DIVERGED — NaN/Inf** |
   | `run_0087_3d_sst_cycle_thick_ztaper/` | **DIVERGED — NaN/Inf** |
   | `run_0088_3d_sst_cycle_cfl01/` | **DIVERGED — NaN/Inf** |

   `run_0082` の全9体積スナップショットの `VALUE/*` は有限でした。しかし終端の `rms_roUy=3.97e-5`、`rms_roOmega=7.58e16`。完走の証拠であって、収束の証拠ではありません。

   また、9月8日に [SST 既定値](/home/sano/work/forge/solver_density_cuda/input/solverConfig.hpp:202) と [スカラー拡散の相対ガード](/home/sano/work/forge/solver_density_cuda/cuda_forge/scalarTransport_d.cu:104) が変更されています。9月6日の結果から現行ソルバの限界を断定できません。`L_sw` を変えれば流れ自体も変わるため、成功だけで原因は確定しません。

   **対案:** 状態を「加速点で完走、生産点は未成立」に戻す。バイナリ・実効設定を固定し、現行実装で同一幾何の再現試験を先行する。`L_sw` 分離や鈍頭化は、その後の比較対象とする。「同位置終了を除外しても設計損失なし」は未検証として撤回する。

4. **Major — 擬似時間の履歴を物理時間として解釈している。**

   **根拠:** [runner_sern.py:152](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:152) は `unsteady:0, dualTime:0`。`run_0088` の実 config も同じです。したがって、[plan:362](/home/sano/work/forge/plans/active/tooling-nozzle-sern-chain.md:362) の「CFL 半減で破綻 step が倍＝同じ物理時刻」は成立しません。§4.7 の RSS/FSS の平均・振幅も、この反復履歴から物理的な統計として取得できません。

   圧力床への到達が先に観測されたことも、床が根因である証明ではありません。既存の [正値性ガード計画:55](/home/sano/work/forge/plans/accepted/time_integration-update-positivity-guard.md:55) は、別ケースで床 NaN を陰的反復不安定の終端症状と判定しています。

   **対案:** 定常反復の安定性と物理的非定常を分離する。定常解が得られず振動を採用する場合は、物理時間を持つ dual-time に移行し、時間刻み・内部反復・統計窓の独立性を確認する。床到達は診断指標として扱い、更新前後の保存量・流束・乱流源項で原因を切り分ける。

5. **Major — 「3D で推力 −4.5%」には集計対象の変更が混入しており、2% 基準の判断に使えない。**

   **根拠:** [runner_sern3d.py:197](/home/sano/work/forge/design/forge_design/evaluate/runner_sern3d.py:197) は幅外の機体下面も総推力へ加算します。2D にはその面がありません。保存された全力履歴を `check_quasisteady.py` の `classify` 関数へ渡して再判定すると、次の各系列は **`STEADY`** でした。

   | 帳簿の検算 | 末尾平均 \(C_T\) |
   |---|---:|
   | `run_0027_diag3d_euler_accel_slip/` 総計 | 0.932472 |
   | 同 run の幅外ランプ寄与 | −0.025349 |
   | 同寄与を除いた値 | 0.957821 |
   | `run_0029_ref2d_euler_node_accel_slip/` | 0.976214 |

   差は **−4.48% → −1.88%** に変わります。ただし両 run の収束 VERDICT は **NOT CONVERGED**。これは帳簿の不一致を示す検算であり、2D 設計の十分性を証明する値ではありません。

   さらに [mesh_sern3d.py:186](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:186) は外側計算領域の端まで `ramp` を生成します。`Z_ext` を変えると集計する機体面積まで変わります。

   **対案:** ノズル力と機体力をタグ・積分とも分離し、物理的な機体幅を遠方境界位置から独立させる。同じ対象面・基準点・正規化で比較し、閉じた制御体積の運動量収支でも検算する。

6. **Major — 外部流ブロックの追加は必要だが、「実装・検証済み」とする条件が不足している。**

   **根拠:** 2D は `ext_top` を持つ一方、[3D メッシャ:192](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:192) はランプ後縁以降も `top_out` で閉じています。2D で問題視した構成が3Dに残っています。

   また、[plan:452](/home/sano/work/forge/plans/active/tooling-nozzle-sern-chain.md:452) は2形状の剥離ゼロから「作動点の性質なので形状によらず剥離しない」と一般化しています。カウル後縁の圧力だけでは、ランプの逆圧力勾配、入口境界層、側壁干渉を拘束できません。

   **対案:** 3Dにも対応する上側外部領域を設け、固定した機体形状について遠方境界・出口距離・格子の感度を確認する。剥離制約を要求外として外す判断は維持してよいものの、「設計箱全域で剥離しない」とする根拠には使わない。メッシュ品質 PASS と領域独立性を別ゲートにする。

7. **Major — §6 は MOC の整合性検証と、最大推力・CFD 妥当性の検証を混同している。**

   **根拠:** `run_sern_moc_tests.py` は、CSV 保存だけをメモリ上で除外して再実行し **ALL PASS** でした。MLN の面積比は `3.60064 / 3.60031`、推力検算は `2.26635 / 2.26702`。MOC の基礎部分は支持できます。

   しかし [同テスト:117](/home/sano/work/forge/design/tests/run_sern_moc_tests.py:117) の最適性掃引は **assert なし**です。[plan:387](/home/sano/work/forge/plans/active/tooling-nozzle-sern-chain.md:387) の「Rao 点が多作動点・粘性 MOO の Pareto 端点に出なければ実装誤り」も不適切です。固定長・単一作動点の非粘性問題とは目的と制約が違います。Shyne の対象も外部流を含む特定の定式化です。[NASA TM-103175](https://ntrs.nasa.gov/citations/19900015790)

   加えて、現行 `check_quasisteady.py --quantity C_T,C_L,C_M` は実行すると **`ERROR: unknown quantity`**。独自の `steadiness` と正式ツールは同じ判定器ではありません。

   **対案:** Rao 検証は同一ガス・作動点・長さ拘束で独立に行い、式11–15との対応と制約付き微小摂動を検算する。CFD は力係数抽出を正式ツールへ接続し、摩擦込みの力・モーメント・剥離量を判定する。格子・領域誤差には、例えば **\(|\Delta C_T|<0.002\)** のように最適化の差より小さい許容値を明記する。

8. **Major — 最適化する問題が文書とコードで違い、制約にも抜けがある。**

   **根拠:** [plan:77](/home/sano/work/forge/plans/active/tooling-nozzle-sern-chain.md:77) は6変数・4目的ですが、[driver_sern.py:33](/home/sano/work/forge/design/forge_design/opt/driver_sern.py:33) は5変数、目的は重み付き推力と長さの2個です。`theta_c` は [moc_sern.py:387](/home/sano/work/forge/design/forge_design/geometry/moc_sern.py:387) で場から決まります。

   モーメント制約も作動点ごとの値ではなく加重平均だけです。ある作動点のトリム不能を別作動点で相殺できます。さらに `L_ramp_max` は粗い probe の検査に留まり、本評価で得た長さには再適用していません。[driver_sern.py:80](/home/sano/work/forge/design/forge_design/opt/driver_sern.py:80)、[同:138](/home/sano/work/forge/design/forge_design/opt/driver_sern.py:138)

   **対案:** 初版を**5変数・推力効率と長さの2目的**として明文化する。モーメントは作動点別の許容窓にし、許容値を「候補の何割を残すか」から決めない。最終輪郭生成後に包絡・長さを厳密に再検査する。

9. **Major — 「3D 最適化には MOC への帰還が必須」という追加方針には根拠がない。**

   **根拠:** [plan:497](/home/sano/work/forge/plans/active/tooling-nozzle-sern-chain.md:497) は、3D 排除効果の帰還まで実装して初めて3D最適化としています。しかし形状を \(G(d)\)、評価を \(J_{\mathrm{3D}}(G(d))\) とすれば、3D評価で \(d\) を更新する時点で形状は3D環境に応答します。MOCへの帰還は必要条件ではありません。

   横方向膨張を「等価2D面積」へ縮約する案には、どの保存則・圧力分布を保つかの定義もありません。現在仕様の [methods/design/overview.md:799](/home/sano/work/forge/methods/design/overview.md:799) は無帰還・3D確認のままで、残作業表と完了条件が別の計画になっています。

   **対案:** MOC は形状パラメータ化として維持する。検証済み3D評価器を得た後、まず `L_sw` を追加した小規模探索で設計改善を測る。帰還は、既存パラメータ族の不足を実測してから別計画で扱う。多成分・MOO・角部処理は既存 accepted 計画を再利用し、再開発しない。

推奨は一つです。**本計画を「評価器の成立確認を先に完了し、その後に最適化を再開する計画」へ組み直してください。** 優先順は次のとおりです。

| 優先 | 実装再開前に確定すべき内容 |
|---|---|
| 1 | 発散採用の撤回、有限値・全残差・実目的量のゲート、力の集計範囲 |
| 2 | サイクル作動点の凍結 TP 化と外気の空気物性 |
| 3 | 固定した機体形状に対する2D/3D外部領域と格子・領域独立性 |
| 4 | 現行バイナリで生産3作動点の3D SSTを検証。必要ならソルバ修正を別計画化 |
| 5 | 小規模探索で判別能力を確認してから MOO 再取得。早期停止・帰還は後段 |

ファイル変更・新規 CFD 投入は行っていません。上記方針は **plan 未反映**です。既存 run の索引は [case/46.sern_design/README.md](/home/sano/work/forge/case/46.sern_design/README.md)。`run_0082` のメッシュ再判定は **PASS（AR 535.8、skew 0.365）**ですが、収束判定を代替しません。

指摘数: Critical 2 / Major 7 / Minor 0
