# codex レビュー: tooling-rerun-conditions (result)

- **plan**: [`plans/active/tooling-rerun-conditions.md`](../../plans/active/tooling-rerun-conditions.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-06
- **commit**: `8ef55aea` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 8.3 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m3

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**

主要 run の準定常判定と数値は概ね再現できました。ただし、無変更時の回帰、外部依存の検査漏れ、Euler 比較元の準定常性の裏付け不足が残っています。現状のまま `accepted` へ移すことには反対します。

1. **Major — 無変更の再実行が、正当な `cfl ≠ cfl_pseudo` 設定を拒否する。**

   [rerun_conditions.py:820](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:820) は参照の `cfl` を推奨値にし、[同:846](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:846) で両 CFL に一致を要求します。

   `run_0094` の config 読み込みだけをメモリ上で `cfl: 4.0, cfl_pseudo: 5.0` に変更すると、**条件変更なしでも停止し、`--cfl 4.0` を提案**しました。この提案は実効 CFL を 5→4 に変えます。定常計算で効くのは `cfl_pseudo` です（[solver-settings.md:38](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:38)）。

   **対案:** 無変更時は両値をそれぞれ保持する。推奨との照合は、条件変更に対して指定した実効 CFL に適用し、異なる２値を持つ参照の無変更試験を追加する。

2. **Major — 外部依存の拒否がファイル拡張子に依存し、入力契約を迂回できる。**

   [_file_refs:295](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:295) は列挙した拡張子で終わる文字列しか検査しません。そのため、以下を参照 config にメモリ上で追加すると、いずれも `build_plan` が受理しました。

   - `physProp.speciesDBFile: /tmp/review_species_db`
   - `physProp.chemistry.mechanismFile: /tmp/review_mechanism`

   solver はこれらを拡張子に関係なくファイル参照として読みます（[solverConfig.cpp:1061](/home/sano/work/forge-integ-1005/solver_density_cuda/input/solverConfig.cpp:1061)、[同:1090](/home/sano/work/forge-integ-1005/solver_density_cuda/input/solverConfig.cpp:1090)）。相対参照なら複製漏れ、絶対参照なら外部依存が残り、§4.2 の「作成前に拒否」を満たしません。

   **対案:** ファイル参照は設定キーの完全修飾パスで検査する。未対応の依存キーは、拡張子・存在の有無によらず作成前に拒否する。

3. **Major — 腕 E の比較元について、準定常性を確認できない。**

   `case/45.isobutane_m6_d155/run_0086_euler_wallfit_pincal_r1_ext6k/` の [quantities_series.csv](/home/sano/work/forge/case/45.isobutane_m6_d155/run_0086_euler_wallfit_pincal_r1_ext6k/quantities_series.csv:1) は **step 6000 の１行だけ**です。§6 の規約で再実行した結果は、

   ```text
   TRANSIENT-UNSETTLED
   only 1 snapshot(s) (<10)
   ```

   残差判定も `NOT CONVERGED (stalled/plateau)` です。これは実流れが過渡だという証明ではなく、**比較元の準定常性が監査資料から判定できない**という問題です。腕 E の２量が基準スナップショットに近いことだけでは、[plan:183](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:183) の「状態変換は Euler 解を保つ」を裏付けきれません。

   **対案:** 比較元の時系列を回収し、出口 M・流量を同じ規約で判定する。資料がなければ元条件の対照計算で確認し、両者の末尾平均を比較する。それまでは「基準スナップショットに対する２量の差が閾値内」に主張を限定する。

4. **Minor — `--Ps` が参照 Ps と同値だと、出口 Pt の同期を省く。**

   [rerun_conditions.py:616](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:616) は Ps が変わった場合だけ Pt も更新します。参照出口を `Ps: 2237, Pt: 3000` として `--Ps 2237` を指定すると、**Pt は 3000 のまま、変更記録は空、推奨は `none`**になりました。

   **対案:** `--Ps` 指定時は Ps・Pt を個別に照合し、いずれかが変われば変更として記録する。なお現行出口カーネルは逆流も静圧アンカーで扱っているため、これを直ちに流れの破壊とは評価しません（[boundaryCond_d.cu:593](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/boundaryCond_d.cu:593)）。「逆流用 Pt」という文書説明も現行実装に合わせるべきです。

5. **Minor — 現在仕様・手順に旧説明が残る。**

   [methods/design/overview.md:977](/home/sano/work/forge-integ-1005/methods/design/overview.md:977) は Euler を「乱流なし、`prepare_info` に `viscous` なし」と説明していますが、実装では壁・輸送設定で分類し、`prepare_info` は照合専用です。また [procedures/nozzle-design-workflow.md:143](/home/sano/work/forge-integ-1005/procedures/nozzle-design-workflow.md:143) の Pt 変更推奨には、Tt・組成・凝縮との併用禁止条件がありません。

   **対案:** 分類説明を一本化し、手順の推奨行にも適用条件と複合変更の未検証範囲を明記する。

6. **Minor — 発散の step 表記と NaN の対象が原ログと対応していない。**

   全残差 CSV を確認すると、最初の非有限値はそれぞれ次のとおりでした。共通の接頭辞は `case/45.isobutane_m6_d155/` です。

   | run | CSV の最初の非有限 step | 列 |
   |---|---:|---|
   | `run_0121_rerun_pt08_scale` | 333 | `rms_roe` |
   | `run_0122_rerun_pt08_noscale` | 119 | `rms_roe` |
   | `run_0129_rerun_pt08_scale_full` | 467 | `rms_roe` |
   | `run_0130_rerun_pt08_noscale_full` | 1 | `rms_roe` |

   [plan:183](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:183)・[同:185](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:185) は 334／120／468／2 と記載しています。更新回数なら説明が必要です。また `rms_ro` は全行有限で、場の「ro NaN・54節点」を確認できるスナップショットは取得済み資料にありません。

   **対案:** CSV の step と更新回数を区別し、残差 NaN と場の NaN を分けて記録する。位置・節点数の主張には元スナップショットまたは診断出力を添付する。

確認できた検証結果は以下です。数値は CSV の末尾５枚平均、判定は現行ツールで再実行しました。各 run は `case/45.isobutane_m6_d155/` 配下です。

| run | 残差 VERDICT | 対象量の VERDICT・数値 |
|---|---|---|
| `run_0119_rerun_ctrl` | NOT CONVERGED | δ_E・M・流量 STEADY。δ_E = 0.7252822 |
| `run_0120_rerun_euler_pt08` | NOT CONVERGED | M・流量 STEADY。M = 5.9999977 |
| `run_0132_rerun_pt08_scale_cfl1_ext` | NOT CONVERGED | ３量 STEADY。δ_E = 0.7494484、M = 5.9921242 |
| `run_0134_rerun_euler_tt1500` | NOT CONVERGED | M・流量 DRIFTING |
| `run_0138_rerun_fullpath_blk4` | NOT CONVERGED | M は STEADY、δ_E・流量は DRIFTING |
| `run_0139_rerun_pt08_noscale_cfl1_ext` | NOT CONVERGED | M・流量 DRIFTING。M = 5.9939767 |

判定資料の所在は [case README の run 台帳](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:96) に記載されています。取得済み主要最終場５本では `VALUE/*` の NaN/Inf はありませんでした。A3 の結果を**対象量の準定常性**として扱い、Tt/Y の未解決事項を [§5.1 #11](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:94) に残した判断は妥当です。

指定の `main...HEAD` 差分に CUDA 数値カーネルの変更はありません。実行した `run_mesh_params_tests.py` と `run_physical_wall_analytic_tests.py` はともに FAIL 0。ファイル作成を伴う `test_rerun_conditions.py` 全体は、今回の変更禁止条件に従い再実行していません。

**推奨は、現行の `scale-ic pt + full + 本段 cfl_pseudo 1` 方針を維持し、上記の入力契約修正と比較元の監査を済ませてから再レビューすることです。** 新しい指摘を §5.1 に追加し、移動はその後にしてください。今回は read-only のため **plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 3
