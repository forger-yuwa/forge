# codex レビュー: boundary-conjugate-heat-transfer (result)

- **plan**: [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-23
- **commit**: `4ac92a4c` (feature/cht-phase2-fem2d)
- **codex**: effort `high`, 5.9 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m1
- **focus**: 3 巡目 (notes/reviews/2026-09-23-boundary-conjugate-heat-transfer-result-3.md) の Major 4 + Minor 1 への対処を最優先で検証する。(1) 最終保存の時刻一致・累積 step と暖機の扱い・flux_avg=1 での位相復元、(2) check_cht_balance が保存 SOLID/T を RCM 並べ替えを戻して使っているか、(3) cht_loop の --flux 尊重、(4) cht_wall_series の q_eff 既定と Qf_eff_total 列。再判定は run_0157_r3_avg42 (C3X)・run_0158_rs_cont / run_0159_rs_legA / run_0160_rs_legB (再開回帰 flux_avg=1 warmup=100 非倍数分割)・run_0037_r3_avg42 (Mark II)。Mark II は報告項目 (局所は gate にしない) であることを前提に、C3X の accepted 移行可否を判断してほしい

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**

C3X の G-if・G-cons・準定常判定は再現できました。  
ただし、最終壁出力と収支ゲートに欠陥が残っています。Mark II の局所未収束は、指定どおり承認阻害要因にはしていません。

1. **Major — 最終壁ダンプが保存されず、異なる時刻の収支を PASS にできます。**

   [output.cpp:286](/home/sano/work/forge-cht/solver_density_cuda/output/output.cpp:286) は、壁出力を依然として出力間隔だけで制御しています。流体・固体に追加された最終ステップの強制保存がありません。また、[check_cht_balance.py:156](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:156) は固体の `step` を表示するだけで、壁ダンプとの一致を検査しません。

   **実測:** `case/53.c3x_vane_cht/run_0160_rs_legB/` の流体・固体出力はローカル3975、累積8000ステップまでありますが、壁は `res_wall_5_3015.h5`、すなわち累積7040までです。この960ステップずれた組み合わせを G-cons は **0.0040%／PASS** と判定しました。

   **対案:** 壁にも最終強制保存を適用し、出力に累積ステップを記録してください。G-cons は時刻不一致を `REFUSED` にし、今回の非倍数分割回帰で確認してください。

2. **Major — 固体ハッシュを照合しても、収支計算には別入力の冷却条件・メッシュを使っています。**

   [check_cht_balance.py:132](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:132) は `--solid` の JSON・外部 NPZ から作用素を作ります。一方、ハッシュ照合対象はチェックポイントと固定名 `solid.h5` だけです。[同:187](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:187) の Robin 放熱計算には、照合していない JSON・NPZ が使われます。

   **再現:** `run_0157_r3_avg42` の保存ファイルを変えず、JSON の `h` だけをメモリ上で10%増やすと、ハッシュ不一致による拒否なしに、固体放熱量が **43474.02 → 47821.43 W/m** に変わりました。

   **対案:** ソルバ内連成の検査では `solverConfig.yaml` の `conjugate.solid` を読み、その HDF5 の座標・Robin 辺・`H`・`TC` と保存温度から直接集計してください。必須ハッシュの欠落も拒否対象にします。

3. **Major — 積分済み荷重の NaN を除外して PASS にします。**

   [check_cht_balance.py:92](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:92) は `iface_q_eff` の有限性だけを検査します。しかし実際の集計対象は `iface_Qf_eff` で、[同:143](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:143) が非有限要素を除外します。

   **再現:** C3X 最終ダンプの `iface_Qf_eff[269]`、元の値 **3.18544 W/m** だけをメモリ上で NaN にすると、**「NaN 節点 0/480」「0.0132%」「VERDICT: PASS」**になりました。

   **対案:** 実際に採用する荷重配列の全要素を検査し、NaN/Inf が1点でもあれば集計前に不合格にしてください。この入力を回帰試験に追加してください。

4. **Minor — 前回指摘した仕様・台帳の同期が残っています。**

   [methods/boundary.md:462](/home/sano/work/forge-cht/methods/boundary.md:462) の設定例には `tol_solid` がなく、[同:476](/home/sano/work/forge-cht/methods/boundary.md:476) は現在の残差補正方式と異なる `D_f` 凍結説明です。[plans/README.md:30](/home/sano/work/forge-cht/plans/README.md:30) も `draft`・ソルバ内薄肉シェルのままです。

   [C3X の run 台帳](/home/sano/work/forge-cht/case/53.c3x_vane_cht/README.md:779)には今回の `run_0157`～`run_0160`、[Mark II の台帳](/home/sano/work/forge-cht/case/54.markii_vane_cht/README.md:81)には `run_0037_r3_avg42` がありません。

   **対案:** 設定例・更新式・チェックポイントを含む再開手順・今回の成果物を同期し、C3X V4b の承認範囲と V2・V3・V6 等の残作業を分けて明記してください。

再実行した判定は以下です。壁温の準定常許容は `drift=osc=0.0006`、熱量は `0.002` です。

| 対象 run | G-if | G-cons¹ | 壁温最大 | 瞬時 `Qf_eff_total` | 42更新平均の連成熱量 |
|---|---|---|---|---|---|
| `case/53.c3x_vane_cht/run_0157_r3_avg42/` | **PASS** | **PASS：0.0059%** | **STEADY** | **STEADY** | **STEADY** |
| `case/54.markii_vane_cht/run_0037_r3_avg42/` | **NOT CONVERGED** | **PASS：0.1194%** | **OSCILLATING：663.3 ± 0.28 K** | **DRIFTING** | **STEADY** |

¹ この2 run は壁・固体とも40000ステップです。保存固体 HDF5 から直接計算した Robin 放熱量でも照合しました。

前回の `--flux` 問題は、平均なし・4枚平均とも修正を確認しました。RCM の逆変換、`flux_avg=1` の位相復元、累積暖機も改善しています。再開回帰の保存固体界面温度差は RMS **0.000507 K**、最大 **0.001759 K**でした。

両翼の流体残差は **`NOT CONVERGED (stalled/plateau)`**。確認した保存場に NaN/Inf はありません。`test_solid_shell.py`・`test_solid_fem2d.py` はともに **`PASS (all)`**です。

推奨は、**`active` に留め、1～3を修正・回帰検証し、4の文書同期後に C3X V4b を再判定すること**です。ファイルは変更しておらず、提案は **plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 1
