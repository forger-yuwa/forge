# codex レビュー: boundary-conjugate-heat-transfer (result)

- **plan**: [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-23
- **commit**: `d48a3d48` (feature/cht-phase2-fem2d)
- **codex**: effort `high`, 6.6 min, rc=0
- **判定**: **GO**, 指摘 C0/M0/m2
- **focus**: 5 巡目 (notes/reviews/2026-09-23-boundary-conjugate-heat-transfer-result-5.md) の Major 2 + Minor 1 への対処を検証する。(1) run の solverConfig.yaml から conjugate.mode: fem2d を判別して固体チェックポイントを必須にした、(2) step と content_sha1 を両側で必須にし欠落・不一致を REFUSED、(3) 固体ファイルを conjugate.solid から解決。codex が示した 4 通り (チェックポイント削除 / state の step 削除 / state の content_sha1 削除 / solid.h5 の content_sha1 削除) はすべて REFUSED になることを当方で確認済み。**C3X の accepted 移行可否を明確に述べてほしい**。Mark II は報告項目 (局所は gate にしない)。残る指摘が Minor だけなら GO としてよいか

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO**

**C3X の Phase 2／V4b は `accepted` に移してよい**と判断します。前回の Major 2件は解消し、登録済みゲートを再現しました。  
残件は Minor のみです。承認範囲は C3X V4b とし、V2・V3・V6 等の未完了作業は残してください。

1. **Minor — G-cons に不要な JSON／NPZ 依存が残っています。**

   [check_cht_balance.py:130](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:130) は、ソルバ内 `fem2d` でも先に `--solid` の JSON と NPZ を要求します。必要な HDF5 とチェックポイントが揃った `run_0161_r4_avg42` でも、JSON 指定を省くと終了しました。

   **対案:** config による経路判別を先に行い、ソルバ内経路では固体 HDF5 の界面情報・Robin 条件と保存温度だけで検査してください。これは可搬性の問題で、今回の数値判定を覆しません。

2. **Minor — plan 本文に旧求解方式が残っています。**

   [plan §4.6a:293](/home/sano/work/forge-cht/plans/active/boundary-conjugate-heat-transfer.md:293) は直接求解式、[同:301](/home/sano/work/forge-cht/plans/active/boundary-conjugate-heat-transfer.md:301) は「右辺だけ更新」「最終更新では必ず再分解」のままです。実装は [conjugateWall.cpp:751](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:751) の残差補正形で、再分解条件にも最終更新の分岐はありません。

   **対案:** 更新済みの [methods/boundary.md:465](/home/sano/work/forge-cht/methods/boundary.md:465) に合わせて本文を同期してください。C3X V4b の承認と、計画全体の完了も明記して区別します。

前回指摘の再検証結果は以下です。ファイル欠落・属性変更は**メモリ上で模擬**しました。

- 指定の4通り――チェックポイント欠落、状態の `step` 欠落、状態・固体それぞれの `content_sha1` 欠落――はすべて **`REFUSED`**。
- 負の step、時刻不一致、ハッシュ不一致、壁の `step_abs` 欠落、積分荷重1点の NaN も **`REFUSED`**。
- `conjugate.solid` を別名にした試験では、その指定先を読み **`PASS`**。

保存済み結果も再判定しました。準定常の `drift`・`osc` は壁温 **0.0006**、熱量 **0.002** です。

| run | G-if | G-cons | 壁温平均／最大 | 瞬時 `Qf_eff_total` |
|---|---|---|---|---|
| `case/53.c3x_vane_cht/run_0161_r4_avg42/` | **PASS** | **PASS：0.0016%** | **STEADY／STEADY** | **STEADY** |
| `case/54.markii_vane_cht/run_0038_r4_avg42/` | **NOT CONVERGED** | **PASS：0.0232%** | **STEADY／TRANSIENT-UNSETTLED** | **DRIFTING** |

C3X の G-if は末尾80更新で、絶対残差 **52.13 W/m²**、相対残差 **2.273×10⁻⁴**、温度更新 **0.002422 K**、固体内部残差 **2.410×10⁻⁴ W/m**。すべて登録許容内です。Phase 1 比も壁温平均 **+0.20914 K**、局所最大差 **0.51426 K** で許容内でした。

再開回帰の `case/53.c3x_vane_cht/run_0162_rs4_cont/` と `run_0164_rs4_legB/` は累積8000 stepで揃い、固体界面温度差は RMS **0.000677 K**、最大 **0.003247 K**。再開側 G-cons は **`PASS：0.0004%`** です。

両翼の流体残差は **`NOT CONVERGED (stalled/plateau)`**。今回の GO は登録どおり界面ゲートと準定常性によるものです。両翼各9枚の保存流体場の `VALUE/*` と最終固体温度に非有限値はなく、`test_solid_shell.py`・`test_solid_fem2d.py` はともに **`PASS (all)`** でした。根拠の恒久索引は [C3X 台帳](/home/sano/work/forge-cht/case/53.c3x_vane_cht/README.md:781) と [Mark II 台帳](/home/sano/work/forge-cht/case/54.markii_vane_cht/README.md:84) です。

推奨は、**C3X V4b を承認し、Minor は残作業として追跡すること**です。Mark II の局所未収束は指定どおり承認阻害要因にしません。計画全体を `done` にする際は、未完了項目の移管先を残してください。ファイルは変更しておらず、本レビューは **plan 未反映**です。

指摘数: Critical 0 / Major 0 / Minor 2
