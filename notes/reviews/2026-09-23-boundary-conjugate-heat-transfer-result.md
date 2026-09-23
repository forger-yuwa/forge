# codex レビュー: boundary-conjugate-heat-transfer (result)

- **plan**: [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-23
- **commit**: `9525377e` (feature/cht-phase2-fem2d)
- **codex**: effort `high`, 8.8 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M7/m2
- **focus**: Phase 2 (ソルバ内 fem2d 連成) の §4.6a 設計・§5.1 #66-#72 の実測・§6 V4b の合否に集中する。特に (1) 流束の 42 更新移動平均が固定点や遅れの点で妥当か、(2) D_f を分解と一緒に凍結する実装が正しいか、(3) G-if が壁温の単調ドリフト中でも通り得ることの扱い、(4) Mark II を報告項目にした判断

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**

`D_f` の左右辺凍結と C++ 固体単独求解は確認できましたが、連成時の固定点・界面荷重・合否判定に未解決の欠陥があります。  
C3X の「全登録項目 PASS」は現状では成立せず、`accepted` への移動は見送るべきです。

1. **Major — 分解再利用で `K_s(T)` が凍結され、現在の物性に対する固定点を保存していません。**

   [conjugateWall.cpp:633](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:633) は温度差だけで再分解を判断します。毎更新の残差は現在の `K_s(u)` で計算しますが、求解は古い分解を使います。`D_f` の不一致は直っていますが、実際の固定点は `K_s(u_fact)u=b+Q` のままです。[plan §4.6a:303](/home/sano/work/forge-cht/plans/active/boundary-conjugate-heat-transfer.md:303) にある「最終出力前は必ず再分解」も実装されていません。

   **実測:** `case/53.c3x_vane_cht/run_0150_dffrozen_avg42/res_solid_5_40000.h5` を Python の現在物性で再組立てすると、内部残差最大は **0.0317565 W/m**。同じ保存荷重を C++ の Picard 求解で解き直すと **4.55e−11 W/m** になり、温度は最大 **0.03679 K** 動きます。固体単独の組立て自体は Python と相対差 **3.59e−14** で整合しており、問題は連成側の再利用方法です。

   **対案:** 古い分解は現在の方程式の残差補正に使い、絶対温度を古い行列で解き続けない構成にする。最終判定前には現在物性で再分解・再求解し、内部残差を確認してください。

2. **Major — G-if が固体内部残差と平均バッファの待機条件を検査せず、未完成の状態を PASS にします。**

   [check_cht_interface.py:117](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_interface.py:117) は `tol_solid` と `res_solid` の両方がある場合だけ内部残差を判定します。実際の `run_0150` の登録 JSON に `tol_solid` はなく、末尾80更新の内部残差最大 **0.0346393 W/m** は判定から外れています。許容を指定しても列が欠けていれば検査を省略します。

   また、[同:95](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_interface.py:95) は行数しか見ず、V4b(g) の「バッファ充填後、さらに `2N` 更新待つ」を実装していません。メモリ上の合成入力で、**開始直後の80更新・内部残差 `1e9 W/m`** が `VERDICT: PASS` になることを再現しました。

   **対案:** `fem2d` では内部残差の許容と列を必須にし、欠落は `REFUSED`。平均窓長・充填数・更新番号も記録して、待機期間を除いた連続80更新だけを判定してください。

3. **Major — `q_eff` を作る面積と固体荷重へ戻す辺長が異なり、界面熱量が保存されません。**

   [conjugateWall.cpp:225](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:225) は節点積分熱量を流体の合成面積 `surfArea` で割ります。一方、[同:602](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:602) は固体の集中辺長を掛けます。曲面・角では両者が異なるため、転送荷重に `L_i/A_i` が混入します。これは §4.4b が警告した幾何の問題そのものです。

   **実測:** `case/54.markii_vane_cht/run_0034_insolver_fem2d_avg42/` の角 `(0.064912, −0.000686) m` では **`L/A=1.29891`**。最終壁ダンプの節点荷重は、本来 **2.02762 W/m** のところ **2.63370 W/m** になります。全周積分への影響は約 **0.0052%** と小さく、積分ゲートだけでは局所の約30%増幅を見逃します。

   **対案:** 流体側の積分済み荷重 `Rraw−Fw−e_w Rrho` を直接、一度だけ固体へ渡す。局所残差を熱流束に換算する段階で固体の集中辺長を使い、角を含む保存試験を追加してください。

4. **Major — 再開用 HDF5 は書くだけで、固体状態を復元する経路がありません。**

   [conjugateWall.cpp:451](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:451) は起動時に内部温度を壁温平均で初期化します。`conjugate_state` の読込み処理はありません。[同:879](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:879) の保存内容も `SOLID/T`・step・界面ハッシュだけで、`D_f`、分解基準温度、42更新の荷重履歴は保存していません。

   したがって、CSV で壁温だけを引き継いでも、連成状態と平均の履歴は失われます。仕様の「再開」としては未実装です。

   **対案:** 固体温度・凍結係数・平均バッファ・更新位相を保存／復元し、固体入力ハッシュを検査する。連続実行と途中再開の同値試験を通してください。

5. **Major — 準定常判定の百分率指定が100倍緩く、判定対象の熱量も異なっています。**

   [check_quasisteady.py:543](/home/sano/work/forge-cht/solver_density_cuda/tools/check_quasisteady.py:543) の引数は割合です。plan #71 の **`--drift 0.06 --osc 0.06` は6%** であり、登録した壁温 **0.06%** には `0.0006` が必要です。

   登録値で `conjugate_history.csv` を再判定した結果は次のとおりです。

   | run | 壁温平均 | 最大壁温 | 総界面熱量 |
   |---|---|---|---|
   | `case/53.c3x_vane_cht/run_0150_dffrozen_avg42/` | `STEADY` | `STEADY` | `STEADY` |
   | `case/54.markii_vane_cht/run_0034_insolver_fem2d_avg42/` | `STEADY` | **`OSCILLATING`** | `STEADY` |

   Mark II の最大壁温は末尾で約 **663.6 ± 0.4 K**。README の **ALL STEADY は訂正が必要**です。また、[cht_wall_series.py:40](/home/sano/work/forge-cht/solver_density_cuda/tools/cht_wall_series.py:40) は既定が `q_compact` で、`q_eff` を選択できません。`run_0150/wall_series.csv` の最終熱量 **42245.8 W/m** は、連成履歴の **43471.9 W/m** と別の量です。

   **対案:** 壁温は `0.0006`、熱量は `0.002` で再判定し、連成に使った平均 `q_eff` を対象として記録を更新してください。

6. **Major — 方針上の安全停止条件が実装されていません。**

   [conjugateWall.cpp:678](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:678) は非有限と温度下限だけを検査し、§4.6a のガス全温に基づく上限を検査しません。

   また、[同:704](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:704) は「10更新連続で増加」に加え、**直前の1更新から2倍超**を要求します。更新量が毎回1.1倍なら10更新で2.59倍になりますが、この停止条件には入りません。

   **対案:** 増加区間の始点を保存して累積増幅を判定し、登録した温度上限も実装してください。

7. **Major — 区間識別が `local1d` の冷却条件変更を見逃します。**

   [stage_manifest.py:97](/home/sano/work/forge-cht/solver_density_cuda/tools/stage_manifest.py:97) に `conjugate.back` と `conjugate.h_c` がありません。`stage_key()` で **`h_c: 100` と `10000` が同一キーになる**ことを確認しました。固体抵抗が変わる別問題を、同じ収束区間として連結できます。

   **対案:** 固体方程式を変える全入力を実効値で区間キーに含め、冷却条件の変更で区間が分かれる試験を追加してください。

8. **Minor — 現在仕様・索引・残作業表が同期していません。**

   [methods/boundary.md:330](/home/sano/work/forge-cht/methods/boundary.md:330) は `local1d` のみ実装済み、`fem2d` 未実装と記載する一方、後段に実装済みの説明があります。[methods/index.md:22](/home/sano/work/forge-cht/methods/index.md:22) は CHT 全体を未実装としています。plan も `draft` のまま、過去の未実装記述と完了記述が併存しています。

   **対案:** §5.1 を現在の完了／未完了に整理し、今回の指摘と、V2・V3・V6 等の残作業を明示する。Phase 2 の結果承認と計画全体の完了を区別してください。

9. **Minor — `conjugate` の未知キー拒否が未実装です。**

   [solverConfig.cpp:500](/home/sano/work/forge-cht/solver_density_cuda/input/solverConfig.cpp:500) は既知キーを個別に読むだけです。例えば `flux_avgg: 42` は拒否されず、平均なしの既定値へ落ちます。§4.6a の「未知キーは拒否する」と一致しません。

   **対案:** `conjugate` と `conjugate.gate` に許可キー検査を実装してください。

推奨は、**`active` に留め、固定点・保存荷重・ゲートを修正してから V4b を再判定すること**です。

42更新の後方平均は、流体が定常化した極限では一定荷重を変えません。ただし平均の遅れは残ります。`run_0148_coldstart_avg42` の G-if 初回到達 **34150 step** では、登録許容による壁温平均・最大壁温・熱量はすべて **`DRIFTING`** でした。plan の「G-if 到達＝連成完了ではない」という注意は正しく、速度比較も全ゲートを満たした時点で行う必要があります。

Mark II を局所同値試験の合否対象から外す判断は、Phase 1 自身が局所残差を満たしていない以上、妥当です。ただし、これは超音速 CHT の検証完了を意味しません。再実行した判定は **G-if: `NOT CONVERGED`**。流体残差は C3X・Mark II とも **`NOT CONVERGED (stalled/plateau)`** でした。両 run の全残差系列と最終場には NaN/Inf はありません。

壁温差の記録値は再現できました。C3X は Phase 1 比で平均 **+0.19912 K**・局所最大 **0.49649 K**、Mark II は平均 **+0.36956 K**・局所最大 **5.38214 K** です。これらの差が小さいことだけでは、上記の保存・内部残差の問題は解消しません。

ファイルは変更していません。判定器の結果ファイル書込みも抑止しました。提案は **plan 未反映**です。

指摘数: Critical 0 / Major 7 / Minor 2
