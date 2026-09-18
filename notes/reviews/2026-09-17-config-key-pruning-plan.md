# codex レビュー: config-key-pruning (plan)

- **plan**: [`plans/active/config-key-pruning.md`](../../plans/active/config-key-pruning.md)
- **stage**: `plan`
- **date**: 2026-09-17
- **commit**: `4bd1c876` (feature/sern-design)
- **codex**: effort `high`, 5.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M9/m0
- **focus**: 初回 (plan 段)。ユーザ指示「パラメータが多すぎるので採用があり得ないものは消す」で起票。§4 に削除基準 (A 実測で効果なし/有害, B 後継に置換された旧経路で A/B 終了, C 決定で不採用のモード値, D 未使用の内部定数) と残す基準 (P 生産で選ぶ, V 診断・A/B で要る) を置き、§5.2 に全 155 キーの分類を根拠つき (ファイル:行 か run 実測) で作った。削除候補 31・保留 7・残す 117。問い: (A) 分類基準そのものは妥当か、特に「V (診断で要る)」の線引き (§10: 切り分け手順が plan/notes に書かれているものだけ残し、他は git にあるので消す) を採ってよいか。(B) 削除候補 31 の中に「消してはいけないもの」はあるか。根拠が弱い/読み違えているものを指摘してほしい。(C) 保留 7 件の扱い (ducrosLimiter は enable=0 で ducros 場を恒久 0 にするため単純削除できない、condLimiterMode は RK/dual-time が自動降格で旧経路に依存、nodeOmegaWfDirichlet は accepted plan が「残置」と明示決定済み) は妥当か。(D) §6 の検証 (削除キーを既定値で書いていた run のビット不変、旧キーの起動時エラー、config lint、単体試験) で十分か。実装はまだしていない。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

目的は妥当ですが、現在の削除表をそのまま実装してはいけません。時間積分方式の同一視、SST の最終決定の読み違い、棚卸しの取りこぼしがあります。推奨は、**根拠が確定した削除だけに対象を絞り、共有パーサと設定生成器を含めて検証すること**です。

1. **Major — `timeIntegration: 1` は `3` の別名ではない。削除根拠が誤り。**

   **根拠:** [solverConfig.cpp:907](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:907) は `1` を1段 Euler、[同:919](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:919) は `3` を3段 TVD Runge–Kutta として、段数と係数を別々に設定しています。共通カーネルも [timeIntegration_d.cu:464](/home/sano/work/forge/solver_density_cuda/cuda_forge/timeIntegration_d.cu:464) でその係数を使います。同じ分岐に入ることは、時間積分の同一性を意味しません。

   **対案:** 今回の削除対象から外す。Euler を廃止するなら、次数・安定性・計算量の変更を伴う機能廃止として別途判断する。

2. **Major — SST 分離型2キーの削除理由は、引用先の最終決定と逆。**

   **根拠:** [accepted plan §4.0:41](/home/sano/work/forge/plans/accepted/turbulence-sst-energy-includes-k.md:41) は `sstEnergyIncludesK` の既定を **0** と決定し、分離型を個別利用可能としています。[残作業表:111](/home/sano/work/forge/plans/accepted/turbulence-sst-energy-includes-k.md:111) にも「撤去せず現状維持」とあります。現行 [solverConfig.cpp:527](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:527) もこの決定どおりです。本計画は同文書に残った旧実装ステップを採用しています。

   **対案:** `sstIsotropicStress` / `sstEnergyKSource` は今回保留。「非推奨なので廃止する」という新しい判断は可能ですが、「既に後継へ置換済み」として処理しないこと。

3. **Major — 「全155キーを分類済み」は成立していない。**

   **根拠:** [config_key_inventory.py:8](/home/sano/work/forge/solver_density_cuda/tools/config_key_inventory.py:8) の抽出は単純な変数参照しか扱わず、`config["mesh"]["…"]` などを取りこぼします。抽出処理を読み取り専用で再現すると155名称になりましたが、実際に読む `timeIntegration`、`discretization`、`nodeWallDirichlet`、`gradLSQDegenThresh`、`primPack` などが含まれません。同名末端キーも `setdefault(key, …)` で統合されます。

   また、提示表を数えると **A=14キー、B=14キー、D=10キーだけで38キー**。これに C のキー削除・モード廃止が加わり、「削除31」と一致しません。「保留7」は7行であり、キー数ではありません。残す117キーにも個別分類がありません。

   **対案:** 完全修飾パスを識別子にし、セクション・キー・別名・モード値を分離して再集計する。各項目に実効既定値、利用箇所、分類、根拠、移行先、検証を持たせる。未知キー検出の許可リストに現在の抽出結果を使わないこと。

4. **Major — 共有パーサと設定生成器が影響範囲から漏れている。**

   **根拠:** [convertGmshToForge.cpp:26](/home/sano/work/forge/solver_density_cuda/mesh/convertGmshToForge.cpp:26) は同じ `solverConfig::read()` を使い、[同:58](/home/sano/work/forge/solver_density_cuda/mesh/convertGmshToForge.cpp:58) は `cfg.gpu` を確保処理へ渡します。[variables.cpp:288](/home/sano/work/forge/solver_density_cuda/variables.cpp:288) では値が1なら `cudaMalloc` します。solver の CPU 境界処理が未対応でも、変換器の CPU のみの経路まで不要とは言えません。

   さらに [runner.py:80](/home/sano/work/forge/design/forge_design/evaluate/runner.py:80)、[runner_wt.py:256](/home/sano/work/forge/design/forge_design/evaluate/runner_wt.py:256)、[runner_sern.py:214](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:214) は削除予定の `gpu` / `nodeWallDirichlet` を生成します。起動時拒否だけ先に入れると、新規設計ケースも起動不能になります。

   **対案:** 設計用 YAML の仕様変更は対象外のままでよいが、`solverConfig.yaml` の生成器・テンプレート・変換器は対象に含める。solver の GPU 固定と変換器のホスト実行を分離し、`initial` の必須性も呼出側別に定義する。

5. **Major — `condEquilibrium: 1 → 2` の「固定点同一」は一般には誤り。**

   **根拠:** より新しい [accepted plan:83](/home/sano/work/forge/plans/accepted/condensation-source-limiter-steady.md:83) は、緩和形の定常条件が
   `R_transport + V α θ ρ(g_eq−g)/Δτ = 0`
   となり、**固定点も `Δτ` に依存する**と明記しています。輸送のある問題では代数拘束 `g=g_eq` と同一ではありません。また [solverConfig.cpp:876](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:876) は mode 2 を単一凝縮種に制限しています。

   **対案:** `condEquilibrium: 1` と `condEq*` は今回保留。廃止する場合は「同値な整理」ではなく、旧モデルとその受付範囲の廃止として移行可能範囲を明記する。下流の一例の近接だけを互換性の証明にしない。

6. **Major — `mesh.axisymMethod: 1` は「不採用決定済み」ではない。**

   **根拠:** [対応する active plan:77](/home/sano/work/forge/plans/active/axisymmetric-su2-source-formulation.md:77) は陰解法の収束深さを残課題とし、[同:84](/home/sano/work/forge/plans/active/axisymmetric-su2-source-formulation.md:84) では method 1 を opt-in として保持し、線形化改善後に再評価するとしています。特定の陰解法・条件での不具合は、定式化全体の不採用決定ではありません。

   **対案:** 今回は保留し、既存 plan の責務と残作業を維持する。`bndFirstOrder` を既存の削除 plan に委ねる扱いは妥当です。

7. **Major — A・D・V の基準が、限定的な実測や利用件数から機能廃止へ飛躍している。**

   **根拠:** [line-implicit v2 plan:49](/home/sano/work/forge/plans/accepted/time_integration-line-implicit-viscous-v2.md:49) は粘性結合の効果を「**このケースでは**僅差」と限定し、粘性律速になる条件も記載しています。[凝縮 plan:192](/home/sano/work/forge/plans/accepted/condensation-nonequilibrium.md:192) の二温度評価も希薄水/N2に限定されています。これらは全用途で不要という証拠ではありません。

   D の `tMaxReaction` / `freezeBelowT` も単なる内部反復定数ではなく、[chemistry_d.cu:86](/home/sano/work/forge/solver_density_cuda/cuda_forge/chemistry_d.cu:86) で反応源の温度評価・停止条件を変えます。設定実績ゼロだけでは、6000 K固定・低温凍結機能廃止の妥当性を証明できません。

   **対案:** A は「検証した適用範囲」と「廃止する対応範囲」を対応づける。`lineViscCoupling` / `lineViscousDtRelief` / `condTwoTemp` と上記化学2キーは今回保留する。`condGyarmathyC` / `condN2LiquidCp` の固定化も、標準モデルとして固定する判断と感度試験機能の廃止を明記する。

   V は「文書の有無」ではなく、**現在も必要な切り分け目的・実行手順・解除条件**で判定すべきです。文書がないものは用途を確認して分類し、過去の手順が残っているだけのものは残置理由にしない。旧 commit の復元では、その後の修正を含む現行コード上の A/B を代替できません。

8. **Major — §6 のケースでは、削除する主要経路を検証できない。**

   **根拠:** 指定された3 run の実 config は、すべて **node・SLAU・`turbulence.model: none`・`viscMethod: 0`** です。SST、粘性壁、line-implicit、cell、陽解法、KEEP/Roe の削除回帰を通せません。[検証手順:運用ルール](/home/sano/work/forge/procedures/verification/README.md:32) の node/cell 両系統確認も不足します。

   `check_convergence.py` を実行した結果は次のとおりで、保存済み判定とも整合しました。

   | run | VERDICT |
   |---|---|
   | `case/16.nozzle_wys/run_0476_tp5_node_euler_profile_base/` | `PASS (converged)` |
   | `case/44.vitiated_air_wt/run_0312_passiveG_order1_bdf2_dt8e-6_nsub40/` | `NOT CONVERGED` |
   | `case/09.Taylor-Green/run_0156_passiveG_fct_step_seam/` | `NOT CONVERGED` |

   後二者は非定常試験なので、この判定だけで試験失格とはしません。ただし定常収束の証拠には使えません。

   **対案:** 削除項目ごとの実行経路に対応する最小回帰表を作る。SST・壁出力・cell/node・対象時間積分を実際に作動させる。定常比較は全残差の `check_convergence` と対象量の `check_quasisteady`、非定常は同一物理時刻・保存収支・内部反復依存を判定基準にする。新規複製 run、品質確認、IC、出力先・台帳・残差図も検証手順に含める。

9. **Major — 「長時間 run の保存量ビット一致」だけでは、実行可能で十分な非退行判定にならない。**

   **根拠:** [viscousFlux_d.cu:320](/home/sano/work/forge/solver_density_cuda/cuda_forge/viscousFlux_d.cu:320) などは浮動小数点 `atomicAdd` で残差を集積します。既存 [回帰記録:192](/home/sano/work/forge/plans/accepted/condensation-float-speedup.md:192) でも同一バイナリ反復ノイズを評価しています。長時間ビット一致を唯一の条件にすると、無変更でも通らない可能性があります。

   また、`nodeWallStressEdgeKernel` のような出力経路は保存量だけでは検査できません。`implicitRelaxSST=-1` も [solverConfig.cpp:344](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:344) で `implicitRelax` を継承する指定であり、固定数値への置換とは異なります。

   **対案:** 旧キーありの旧 config と、当該キーを除いた新 config を明示的に対応づける。決定的な局所試験ではビット一致、CFDでは同一バイナリ反復幅と事前に定めた物理量許容値の両方を使う。全保存量に加え、壁応力・熱流束・診断量も対象にする。旧キー拒否試験は完全修飾パスごとに行い、残す同名キーや残すモード値の受理も確認する。

保留方針については、特に `ducrosLimiter`、`condLimiterMode`、`nodeOmegaWfDirichlet` の保持に同意します。それぞれ [センサのゼロ化](/home/sano/work/forge/solver_density_cuda/cuda_forge/ducrosSensor_d.cu:32)、[旧経路への自動降格](/home/sano/work/forge/solver_density_cuda/main.cpp:1122)、[明示的な残置決定](/home/sano/work/forge/plans/accepted/turbulence-sst-su2-taw-coupling.md:380) を確認しました。残りも今回は保持し、`draft` かどうかではなく保留解除条件を記載してください。

**推奨は、根拠が確定した削除に限定した段階実施です。** 実装前の修正順は、①誤分類・廃止根拠の修正、②完全修飾キーによる全分類の再作成、③共有パーサ・生成器の移行設計、④経路別回帰と定量ゲートの確定、です。`lineKFreeze` などの既定値変更は、挙動不変の削除と分けてください。

ファイル変更・新規CFD実行はしていません。レビュー提案は **plan 未反映**です。

指摘数: Critical 0 / Major 9 / Minor 0
