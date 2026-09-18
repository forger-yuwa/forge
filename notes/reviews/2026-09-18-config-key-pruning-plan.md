# codex レビュー: config-key-pruning (plan)

- **plan**: [`plans/active/config-key-pruning.md`](../../plans/active/config-key-pruning.md)
- **stage**: `plan`
- **date**: 2026-09-18
- **commit**: `f8a6d243` (feature/sern-design)
- **codex**: effort `high`, 5.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M7/m1
- **focus**: 2 回目 (1 回目 GO-with-changes C0/M9/m0 を全採用)。1 回目の指摘はすべて §5.2/§6/§2 に反映した。加えて**実装後に自分で 2 度誤りを見つけて訂正**している: (a) 最初に 17 件削除したが、updateGuardAlpha / lowMachThornber / multispeciesRhoYCommonLimiter は元 plan の「処置」が明示的に「opt-in のまま残置」で、M2 と同じ読み違いだった; passiveFctTolAbs は同一 plan・同一機構の 2 キーを保留にしながら 1 つだけ削っていた。(b) turbulentSchmidt は「使用 0」が誤りで実際は 162 run が使用中 (37 run は 0.5 で移行先 Sc_t の既定 0.7 と別値)、physProp.isAxisymmetric も 625 run が使用中で、拒否にすると既存 run の再実行が落ちる。原因は棚卸しツールの使用数がキー名ベースで節をまたいで合算していたこと。ツールを完全修飾パスの PyYAML 集計に作り替え、run は case/run で重複排除し、既定値と数値比較するようにした (live 181 パス / 拒否済み 11 / run config 2796 本)。**最終的に削除したのは 5 パス (blockDPLURDiagCache, blockDPLURDqPack, primPack, lineDtWallRelief, implicitRelaxSST) と定数化 6 件**。検証は (i) 削除キーを書いた config が起動時エラーになること、(ii) 削除前後の場の差が同一設定の反復ノイズ以下、(iii) 単体試験 7 本 PASS、(iv) 残置決定のあるキーが受理されること。問い: (A) 残った削除 5 パス + 定数化 6 件に、まだ「消してはいけないもの」はあるか。(B) 棚卸しツールの数え方に穴は残っていないか。(C) 第 2 陣の候補 (physProp.isCompressible と physProp.ro は必須キーなのに読むコードが無い、time.last.control は消費者ゼロ、time.last.time は読む場所が無いのに 1988 run が書いている、mesh.meshFormat は合法値が 1 つ) の進め方は妥当か。設計 §4 の基準そのものの再指摘は不要。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
削除 5 パスのうち 3 パスに、まだ「opt-in 残置」という既存決定との矛盾があります。  
棚卸しには探索漏れ・誤抽出・既定値比較の欠陥が残り、検証結果も現行の合格基準を満たしていません。

1. **Major — 性能スイッチ 3 件で、残置決定の読み違いが再発している。**

   **根拠:** [performance-3d-node-sst-speedup.md:160](/home/sano/work/forge/plans/accepted/performance-3d-node-sst-speedup.md:160) は、`blockDPLURDqPack` と `primPack` を「不採用で確定、**opt-in 残置**」としています。`blockDPLURDiagCache` も [solverConfig.hpp:118](/home/sano/work/forge/solver_density_cuda/input/solverConfig.hpp:118) に「opt-in 記録用に残す」と明記されています。「既定に採用しない」と「設定機能を廃止する」は別の決定です。

   **対案:** この 3 件は今回の削除から戻してください。廃止するなら、既存の残置決定を置き換える判断として、再現用途を終了する理由と影響する run を明記する必要があります。自己訂正済みの `updateGuardAlpha` 等と同じ扱いが必要です。

2. **Major — 棚卸しが入れ子の検証 run を探索せず、「使用 0」を再び誤判定している。**

   **根拠:** [config_key_inventory.py:129](/home/sano/work/forge/solver_density_cuda/tools/config_key_inventory.py:129) の探索は `case/*/run_*/solverConfig.yaml` に限定されています。現在のツリーだけでも、`run_*` を祖先に持つ **826 config** が探索対象外でした。実際に漏れている使用例は以下です。

   | パス | 指定値 |
   |---|---:|
   | `forge-perf/case/44.vitiated_air_wt/run_0200_perf_regress_node_axisym_sst_tp/final_dcache/solverConfig.yaml:12` | `blockDPLURDiagCache: 1` |
   | `forge-perf/case/16.nozzle_wys/run_0452_perf_bench3d_coarse_bench/cb_dqpack_r1_n200/solverConfig.yaml:29` | `blockDPLURDqPack: 1` |
   | `forge-perf/case/16.nozzle_wys/run_0452_perf_bench3d_coarse_bench/cb_prim_r2_n200/solverConfig.yaml:7` | `primPack: 1` |

   **対案:** 入れ子を含めて探索し、`case` 以下の config 相対パスを識別子にしてください。同一パスでも内容が違う worktree は別設定として残すべきです。今回、現行探索範囲では内容の異なる重複は見つかりませんでしたが、現在の「先勝ち」方式には検出処理がありません。読み取り失敗も黙って除外せず集計不完全として報告してください。

3. **Major — 「live 181 パス」は、実際に受理される設定キーの一覧になっていない。**

   **根拠:** [config_key_inventory.py:46](/home/sano/work/forge/solver_density_cuda/tools/config_key_inventory.py:46) と [同:59](/home/sano/work/forge/solver_density_cuda/tools/config_key_inventory.py:59) はコメントを含むソース全体を走査します。実行確認では、コメントを除くだけで抽出数が **181→177** になり、`time.last.time`、`physProp.isCompressible.ro`、`turbulence.LESmodel`、`turbulence.LESorRANS` が消えました。

   さらに `mesh`・`time.deltaT` 等の節自体や、既に拒否される `time.last.nStep`・`time.nInnerLoop` も live に入ります。必須ヘルパの解析では、`mesh.meshFormat` の既定値が `"mesh"`、`time.last.control` が `"time.last"` になります。これは既定値ではなくエラー表示用の節名です。

   **対案:** コメント、節、値の読み取り、拒否専用参照を区別し、必須・任意ヘルパを別々に解析してください。必須キーの既定値は「なし」とすること。分類表の網羅性確認は、この修正後の一覧でやり直す必要があります。

4. **Major — 非既定値の比較が、小さい値の重大な変更を「既定と同じ」にする。**

   **根拠:** [config_key_inventory.py:105](/home/sano/work/forge/solver_density_cuda/tools/config_key_inventory.py:105) は絶対許容差に最低 `1e-12` を与えます。実際に `same_as_default(0, '1.0e-30')` と `same_as_default(1e-20, '1.0e-30')` がともに `True` でした。`passiveFctTolAbs` のゼロ化や **10 桁の変更**が非既定使用から消えます。

   また [同:121](/home/sano/work/forge/solver_density_cuda/tools/config_key_inventory.py:121) の末端名による初期値補完では、`mesh.renumber → meshRenumber` や `turbulence.turbulentSchmidt → Sc_t` を解決できません。

   **対案:** 数値表記を正規化して比較し、単位に依存する一律の絶対許容差を使わないこと。別名・継承・条件付き既定値は明示的に対応づけ、解決できないものは「不明」としてください。「記載 run 数」と「非既定 run 数」も別々に出すべきです。

5. **Major — 未知キー検出は、plan が「拾える」としている誤配置を検出できない。**

   **根拠:** [check_solver_config.py:44](/home/sano/work/forge/solver_density_cuda/tools/check_solver_config.py:44) は完全修飾パスではなく、末端名の文字列がソースにあるかだけを調べます。以下はいずれも `unknown_keys()` が空リストを返しました。

   - トップレベルの `lowMachPrecond`
   - `physProp.speciesFaceReconstruction`
   - `space.keepDissType`
   - 拒否済みの `mesh.primPack`
   - 読まれない `time.last.time`

   したがって [plan:162](/home/sano/work/forge/plans/active/config-key-pruning.md:162) の「未知キー検出で拾える」は成立していません。既存の Python 単体試験は実行して `ALL PASS` でしたが、この欠陥を検出しません。

   **対案:** 修正した完全修飾パス一覧で検査し、誤配置は WARN、solver が拒否するキーは FAIL にしてください。上記の実例を回帰試験に追加する必要があります。

6. **Major — 自己訂正が設定の正本文書と実施表に反映されていない。**

   **根拠:** [recommended-settings.md:214](/home/sano/work/forge/procedures/recommended-settings.md:214) 以降は、復元した `updateGuardAlpha`・`lowMachThornber`・`multispeciesRhoYCommonLimiter`・`turbulentSchmidt`・旧軸対称キー・`passiveFctTolAbs` を現在も「書くと起動時エラー」と案内しています。`lowMachThornber` の移行先を `lowMachPrecond` とする記述も残っています。

   [plan:94](/home/sano/work/forge/plans/active/config-key-pruning.md:94) の確定削除表と [plan:174](/home/sano/work/forge/plans/active/config-key-pruning.md:174) の検証表にも復元済みキーが残り、「定数化 7 件」のままです。

   **対案:** 変更ログだけでなく、確定表・残作業表・検証表・設定文書を最終対象に同期してください。現状は、次の担当者が正本文書に従うと、復元した設定を再び削除・置換してしまいます。

7. **Major — 非退行の主張が全量では成立せず、最終対象の受入試験も不足している。**

   **根拠:** 以下の `res_200.h5` を再比較しました。

   - 旧: `case/44.vitiated_air_wt/run_0467_sweep_cflp12_nsub10_float/`
   - 新: `case/44.vitiated_air_wt/run_0482_prune_regress_float/`
   - 新の反復: `case/44.vitiated_air_wt/run_0483_prune_regress_float_rep/`

   | 量の最大絶対差 | 旧→新 | 新→反復 |
   |---|---:|---:|
   | `ro` | `1.19209e-5` | `1.28746e-5` |
   | `roe` | `6.11719` | `9.72656` |
   | **`roUy`** | **`3.83759e-3`** | **`3.37625e-3`** |
   | **`P`** | **`2.125`** | **`2.0`** |

   引用された `ro`・`roe` の数値は再現しましたが、**全量が反復幅以内ではありません**。これだけで劣化とは断定できないものの、§6.2 の現行基準では合格にできません。

   3 run の再判定はすべて `check_convergence: NOT CONVERGED`、`check_quasisteady --quantity pmax: TRANSIENT-UNSETTLED` でした。非定常の同時刻回帰には利用できますが、定常解の一致を示す結果ではありません。また全て `model: none`・`viscMethod: 0` で、SST・DES・WMLES の定数化を検証できません。拒否試験の `run_0481_prune_reject_test/forge_run.log` は、今は残置する `lowMachThornber` を拒否した旧実装の記録です。

   **対案:** 比較量・ノルム・物理的許容差・反復数を事前に固定し、全対象量を判定してください。最終削除パスごとの拒否／残置キーの受理試験と、§6 の SST・DES・WMLES・cell/node 経路の具体的な試験を受入条件にすること。定常量を報告する試験には、収束と準定常の両 VERDICT を要求してください。

8. **Minor — `wmlesPrt` の削除理由と移行案内が実装に合っていない。**

   **根拠:** [solverConfig.cpp:450](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:450) は `use turbulentPrandtl` と案内します。しかし [wallLaw_d.cuh:130](/home/sano/work/forge/solver_density_cuda/cuda_forge/wallLaw_d.cuh:130) は、Kader 原式への修正で `Pr_t` が不要になったと明記し、実際の [wmlesWallModel_d.cu:223](/home/sano/work/forge/solver_density_cuda/cuda_forge/wmlesWallModel_d.cu:223) も `wallLaw_kader_tplus(Pr, yp)` を呼んでいます。

   **対案:** これは「定数化」よりも「消費者のない引数の撤去」と分類するのが正確です。`turbulentPrandtl` への移送は案内せず、不要なメンバ・引数まで整理してください。

**推奨は、第 2 陣を止め、第 1 陣と棚卸しの整合を先に取り直すことです。** 実装前の優先順は、①指摘 1・6 の残置判断と正本文書、②指摘 2〜5 の棚卸し・検査、③指摘 7 の受入試験です。

A について、残る `lineDtWallRelief` の既定 0 への固定と `implicitRelaxSST` の継承統合には、確認した差分上、新たな既定挙動の破壊は見つかりませんでした。定数化 6 件にも即時撤回を要する数値変更は見つかりませんが、検証完了とは評価できません。

C の段階移行方針は妥当です。ただし `isCompressible`・`ro`・`time.last.control` は「読まれない」ではなく、**パーサで読まれるが下流の消費者がない**キーです（[solverConfig.cpp:312](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:312)、[同:615](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:615)、[同:646](/home/sano/work/forge/solver_density_cuda/input/solverConfig.cpp:646)）。まず任意化と無効である旨の警告を行い、`mesh.meshFormat` は省略時 `hdf5` としつつ、現在の不正値拒否を維持してください（[main.cpp:1094](/home/sano/work/forge/solver_density_cuda/main.cpp:1094)）。

ファイルは変更していません。本レビューの推奨は **plan 未反映**です。

指摘数: Critical 0 / Major 7 / Minor 1
