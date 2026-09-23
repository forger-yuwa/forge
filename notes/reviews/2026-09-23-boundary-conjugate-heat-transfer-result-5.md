# codex レビュー: boundary-conjugate-heat-transfer (result)

- **plan**: [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-23
- **commit**: `6725c5dc` (feature/cht-phase2-fem2d)
- **codex**: effort `high`, 5.9 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M2/m1
- **focus**: 4 巡目 (notes/reviews/2026-09-23-boundary-conjugate-heat-transfer-result-4.md) の Major 3 + Minor 1 への対処を検証する。(1) 壁ダンプの最終強制保存と step_abs、G-cons の時刻一致 REFUSED、(2) 固体側を solid.h5 の ROBIN/* と保存 SOLID/T から直接集計する経路、(3) 集計する荷重配列の NaN 検査。再判定は run_0161_r4_avg42 (C3X)・run_0162/0163/0164 (再開回帰)・run_0038_r4_avg42 (Mark II)。Mark II は報告項目。**C3X の accepted 移行可否**を明確に述べてほしい。残る指摘が Minor だけなら GO としてよいか、Major が残るなら何か

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**

C3X の V4b 数値条件は再現できましたが、G-cons が検証に必要な情報の欠落を `PASS` にする経路が残っています。  
**C3X の `accepted` 移行は現時点では不可**です。Mark II の未収束は、指定どおり承認阻害要因にしていません。

1. **Major — 固体チェックポイントが無いと、初期壁温から復元した別状態を合格させます。**

   [check_cht_balance.py:165](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:165) はチェックポイントの存在だけで検査経路を選び、無い場合は [同:211](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:211) で `wall_profile` から固体を復元します。ソルバ内連成でも、この外部ループ用経路に入ります。

   **再現:** `run_0161_r4_avg42` のチェックポイントを読み込み時だけ「存在しない」と模擬すると、固体放熱量が **43474.1232 → 43430.6385 W/m** に変わり、**G-cons 0.0984%／`PASS`** になりました。保存固体も時刻一致も確認していません。

   **対案:** `solverConfig.yaml` からソルバ内連成かを判別し、`fem2d` ではチェックポイントを必須にしてください。欠落時は `REFUSED` とし、初期プロファイルへの代替を禁止します。

2. **Major — 時刻・固体同一性の必須属性が欠けても照合を省略します。**

   [check_cht_balance.py:168](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:168) は欠落した `step` を `-1` とし、[同:206](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:206) は非負の場合だけ時刻を比較します。ハッシュも [同:180](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:180) で**双方が存在する場合だけ**比較しています。

   **再現:** 同じ C3X run について、次の各変更を独立にメモリ上で模擬すると、すべて **0.0016%／`PASS`** でした。

   - チェックポイントの `step` を削除、または `-1` に変更。
   - チェックポイントの `content_sha1` を削除。
   - 固体 HDF5 の `content_sha1` を削除。

   また、[同:174](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:174) は依然として固定名 `solid.h5` を読み、設定の `conjugate.solid` を参照しません。前回要求した入力照合は未完了です。

   **対案:** 両側の有効な時刻と非空のハッシュを必須にし、欠落・不正値・不一致をすべて `REFUSED` にしてください。固体ファイルは `conjugate.solid` から解決し、ソルバ内経路はその HDF5 の界面・Robin データと保存温度だけで検査できる形にします。

3. **Minor — 文書同期と「全件修正済み」の記録がまだ不正確です。**

   `tol_solid` の設定例と今回の run 台帳は更新されています。一方、[methods/boundary.md:465](/home/sano/work/forge-cht/methods/boundary.md:465) は直接求解式のままなのに、後段で「上の残差補正形」と説明しています。[同:491](/home/sano/work/forge-cht/methods/boundary.md:491) の再開手順も CSV のコピーだけで、`fem2d` のチェックポイント引き継ぎが抜けています。[plans/README.md:30](/home/sano/work/forge-cht/plans/README.md:30) には `draft` が残っています。

   **対案:** 現在の残差補正式と完全な再開手順に更新し、[plan §5.1 #77](/home/sano/work/forge-cht/plans/active/boundary-conjugate-heat-transfer.md:515) の「全件採用・修正済み」を今回の残件に合わせて訂正してください。C3X V4b の承認範囲と V2・V3・V6 等の残作業も分けて明記します。

再検証結果は以下です。準定常判定は `check_quasisteady.py` を使い、`drift`・`osc` とも壁温 **0.0006**、熱量 **0.002** としました。

| run | G-if | G-cons | 壁ダンプの `Tw_max` | 瞬時 `Qf_eff_total` |
|---|---|---|---|---|
| `case/53.c3x_vane_cht/run_0161_r4_avg42/` | **PASS** | **PASS：0.0016%** | **STEADY** | **STEADY** |
| `case/54.markii_vane_cht/run_0038_r4_avg42/` | **NOT CONVERGED** | **PASS：0.0232%** | **TRANSIENT-UNSETTLED** | **DRIFTING** |

42更新平均の連成熱量は両翼とも **`STEADY`**。両翼の流体残差は **`NOT CONVERGED (stalled/plateau)`** です。全保存流体場の `VALUE/*` と最終固体温度に NaN/Inf はありませんでした。

確認できた改善と数値根拠は次のとおりです。

- C3X の Phase 1 比は壁温平均 **+0.2091 K**、局所最大差 **0.5143 K**。登録許容内です。
- `case/53.c3x_vane_cht/run_0162_rs4_cont/` と、`run_0163_rs4_legA/` → `run_0164_rs4_legB/` の最終出力は累積 **8000 step** で一致。保存固体界面温度差は RMS **0.000677 K**、最大 **0.003247 K**。再開側 G-cons は **`PASS：0.0004%`**。
- `iface_Qf_eff[269]` の NaN、および明示的な時刻不一致は **`REFUSED`**。JSON の `h` を10%変更しても、保存 HDF5 から求める固体放熱量は変わりません。
- `test_solid_shell.py`・`test_solid_fem2d.py` はともに **`PASS (all)`**。

根拠成果物は各 run の `conjugate_history.csv`、`conjugate_state_5.h5`、最終 `res_wall_5_*.h5`・`res_solid_5_*.h5`・`res_*.h5` です。恒久索引は [C3X 台帳](/home/sano/work/forge-cht/case/53.c3x_vane_cht/README.md:781) と [Mark II 台帳](/home/sano/work/forge-cht/case/54.markii_vane_cht/README.md:84) にあります。

推奨は、**`active` に留め、1・2の欠落時拒否を回帰試験付きで修正し、3の文書同期後に C3X V4b を再判定すること**です。ファイルは変更しておらず、指摘は **plan 未反映**です。

指摘数: Critical 0 / Major 2 / Minor 1
