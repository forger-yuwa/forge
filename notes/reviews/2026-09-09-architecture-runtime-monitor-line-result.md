# codex レビュー: architecture-runtime-monitor-line (result)

- **plan**: [`plans/accepted/architecture-runtime-monitor-line.md`](../../plans/accepted/architecture-runtime-monitor-line.md)
- **stage**: `result` (diff base `HEAD`)
- **date**: 2026-09-09
- **commit**: `1a9d75c2` (feature/sern-design)
- **codex**: effort `high`, 5.7 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m3
- **focus**: 実装 diff は未コミットで作業ツリーにある。git diff HEAD -- solver_density_cuda methods procedures plans case/*/README.md で取ること。特に (1) ResidualCsvLogger::latestSummary が flush 直後に正しい outer_end 行を返すか (GPU 経路の rows_/slot と CPU 経路の両方)、(2) setDT_d_wrapper の printCfl=false 化が dt 適応 (adaptDt, dtControl==1) の挙動を変えていないか (定常・非定常・dual-time の 3 経路)、(3) unsteady==1 && dtControl==1 の explicit で monitorCflMax が monitor step 以外に古い値のまま表示されないか、(4) 3 run の CSV 差が run-to-run ノイズと同じ桁という主張の妥当性、を見てほしい。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**

通常の残差表示は正しく、`printCfl` の変更による `dt` 適応の回帰も見つかりませんでした。
ただし、適応時の `dt`／CFL 表示が対応しておらず、「数値変化なし」「速度非劣化」の検証も成立していません。

指定どおり作業ツリーの `git diff HEAD -- ...` を確認しました。ファイル変更・ビルド・新規計算は行っていません。

1. **Major — 適応時に表示する `dt` と `maxCFL` が別の刻みに対応する。**

   [setDT_d.cu:414](/home/sano/work/forge/solver_density_cuda/cuda_forge/setDT_d.cu:414) は適応前の `dt` で求めた `cfl_max` を保持したまま、`cfg.dt` を更新します。[main.cpp:513](/home/sano/work/forge/solver_density_cuda/main.cpp:513) は、その更新後の `dt` と適応前の CFL を並べて表示します。

   例えばコード上、適応前 `dt=1e-6`、`cfl_max=2`、目標 CFL=1、上下限非拘束なら、表示は `dt=5e-7 / maxCFL=2`。表示した刻みに換算した CFL は 1 です。「`cfg.dt` 基準の物理 CFL」という仕様を満たしません。

   **重点③の古い値の再表示はありません。** 非 monitor step は [main.cpp:498](/home/sano/work/forge/solver_density_cuda/main.cpp:498) で戻り、表示する step では explicit 側も同じ条件で更新します。問題は更新頻度ではなく、刻みの対応です。適応・時間加算の順序自体は既存由来で、今回の数値回帰とは区別します。

   **対案:** CFL を評価した刻みも保持し、`dt_CFL`／`maxCFL` を対応づけ、適応後の値は `dt_next` と明示する。`unsteady=1, dtControl=1, monitorInterval>1` と上下限拘束を検証に追加する。現行 Sod は固定 `dt` なので、この問題を検証できていません。

2. **Major — 「run-to-run ノイズと同じ桁＝数値変化なし」は証明になっておらず、一部の掲載値も再現できない。**

   根拠は [plan §9:110](/home/sano/work/forge/plans/accepted/architecture-runtime-monitor-line.md:110)。保存済み新旧 CSV の全行・全 `rms_*` を、`|新−旧| / max(|旧|, 1e-30)` で再計算しました。

   | 対象 run（比較先は各 `_baseline/`） | 最大相対差 | 確認結果 |
   |---|---:|---|
   | `case/16.nozzle_wys/run_0330_monitor_steady/` | `4.4933e-3` | 掲載値を再現。ただし掲載された旧同士の差 `2.2e-3` の約2倍 |
   | `case/05.sod_shock_tube/run_0013_monitor_explicit/` | `1.0` | 掲載値 `0.77` を再現できない |
   | `case/09.Taylor-Green/run_0051_monitor_dualtime/` | `0.19210` | 掲載値を再現。最大絶対差も `5.8280e-9` |

   Sod の最大差は step 2 `outer_begin` の `rms_roUy`：新 `8.36e-28`、旧 `9.18e-11`。ゼロ近傍の相対差です。一方、`rms_ro` の最大相対差は約 `6.58e-7`、`rms_roUx` は `3.50e-6` であり、列をまとめた最大値だけでは比較の意味が失われます。

   指定 run 内には旧版の再試行 CSV が保存されておらず、ノイズ床を追試できません。また、[§6:81](/home/sano/work/forge/plans/accepted/architecture-runtime-monitor-line.md:81) の node 側ビット一致条件を満たさない理由・判定基準変更が整理されていません。

   **対案:** 旧旧・新新・旧新の反復結果を別ディレクトリに保存し、比較式と列別の絶対・相対許容差を明記する。それまでは「数値変化なし」を「保存済み新旧比較で上記差を観測」に修正する。

   なお、3 run と各 baseline の `check_convergence.py` はすべて **`NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)`**。保存済み `CONVERGENCE_VERDICT.txt` とも整合します。`check_quasisteady.py --quantity pmax` は3 run とも **`TRANSIENT-UNSETTLED`**。これは非定常ログ回帰試験そのものの失敗ではありませんが、収束解・定常量の一致を裏付ける結果ではありません。

3. **Major — 壁時計の非劣化という受入条件が未検証。**

   [plan:117](/home/sano/work/forge/plans/accepted/architecture-runtime-monitor-line.md:117) は新版 `3.99 ms/step` と旧版 `1.084 s / 200` を比較しています。しかし旧版の `clock()` は `main()` 冒頭から、新版の時計は初期化・初期出力後の [main.cpp:1688](/home/sano/work/forge/solver_density_cuda/main.cpp:1688) からです。時計の種類と測定区間が両方異なります。

   **対案:** native で新旧を同じ開始・終了位置の壁時計により複数回測定し、中央値とばらつきを保存する。現状の「同等以上」を撤回し、§6(c) を未完了に戻す。

4. **Minor — 初期残差ゼロから増えた列が `worst` 候補から永久に除外される。**

   [main.cpp:531](/home/sano/work/forge/solver_density_cuda/main.cpp:531) は `rms0 > 0` を必須にしています。Sod の `rms_roUy`／`rms_roUz` は初期ゼロから最終約 `3.74e-6`／`3.50e-6` に増えていますが、`worst` には選ばれません。今回の値は小さくても、同じ条件なら大きな増大も有限である限り隠れます。

   **対案:** `初期ゼロ → 現在非ゼロ` を通常の低下桁数とは別に `from-zero` として表示する。無次元化されていない絶対 epsilon を分母に足して順位を作らないこと。

5. **Minor — 定常時の同期削減について文書が実装以上の効果を述べている。**

   [methods/architecture/overview.md:277](/home/sano/work/forge/methods/architecture/overview.md:277) は `printCfl=false` により host 同期がなくなるとしています。しかし実際の条件は `needHostRead = (adaptDt && dtControl==1) || printCfl`。`run_0330` は `dtControl=1, monitorInterval=1` なので、定常でも毎 step の CFL 読み出しが残ります。

   **対案:** 「表示目的の読み出しを除去する。適応目的の読み出しは残る」に修正する。dual-time は [main.cpp:1443](/home/sano/work/forge/solver_density_cuda/main.cpp:1443) で固定 `dt` を強制するため、同経路の「物理時間を適応する」というコメントも訂正する。

6. **Minor — ω 過渡の数値と「1 行目から」という記述が実測と異なる。**

   [plan:118](/home/sano/work/forge/plans/accepted/architecture-runtime-monitor-line.md:118) の `run_0317` の推移は `24→1.8e4→330` ですが、`case/16.nozzle_wys/run_0317_inletprof_tt_h2o/residual_history.csv` の `outer_end` は初期 `24.5926`、最大 `41225.1367`、最終 `330.6501` です。また `run_0330` の最初のモニタ行には `worst rms_roOmega` はありません。

   **対案:** 実測値と表示開始時点に修正する。restart 直後の ω 過渡が既存 run にもある、という定性的説明は支持できます。

確認できた正常点も明記します。GPU の `rows_`／slot は flush 中に値をコピーし、CPU は即時書き出し後に同じ値を保持するため、通常経路の `latestSummary()` は正しい構造です。保存された全500モニタ行で、残差値・`worst` 列・変化桁数の符号を CSV と照合し、不整合はありませんでした。全 CSV と保存済み HDF5 の数値 dataset に NaN/Inf はありません。境界半割面・周期 seam・軸・数値カーネルの閾値を変更する差分もありませんでした。

**推奨は、受入れを保留して修正・再検証することです。** 優先順は①刻みと CFL の対応、②再現可能な数値比較、③同条件の速度測定、④Minor 修正です。§5.1 では検証済み扱いの #1 を再開し、これらを残作業として記録してください。EOS 床・境界質量流量の将来課題は既に残っています。`methods/index.md` と plan 索引のリンクは整合していますが、`check_plans.py` は result レビュー行不足で **`VERDICT: FAIL`**。採否を §6.1 に記録してから完了扱いにしてください。本レビューの提案は **plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 3
