# console モニタ行 — モード別に意味のある実行時情報を 1 行/monitor で出す

## メタ

- **area**: `architecture`
- **status**: `done`
- **related_docs**:
  - [`methods/architecture/overview.md`](../../methods/architecture/overview.md) §8.4 (残差モニタ) / §8.5 (console モニタ行)
- **related_plans**:
  - [`architecture-residual-monitor-async.md`](../accepted/architecture-residual-monitor-async.md) (残差 CSV の device 常駐化・`monitorInterval` の導入)
  - [`architecture-detect-nan.md`](../accepted/architecture-detect-nan.md)
- **created**: `2026-09-09`
- **owner**: `sano`

## 1. 目的

`forge_run.log` に毎ステップ出ている `max cfl` / `dt` は、定常計算 (`unsteady: 0`, 陰解法・陽解法とも) では
計算に使われていない値である (局所擬似時間 `dt_local = cfl_pseudo·V/λ` で `cfg.dt` が打ち消される。
`dtControl: 1` だと `cfg.dt` は `dt_min` に貼り付き、`max cfl` は定数になる。例: case/16 run_0317 は
1000 step 全て `max cfl 22.24 / dt 1e-08`)。一方 log には残差も速度も出ておらず、`tail -f` では何も分からない。
本計画は console 出力を **モード (定常/非定常, 陽/陰) ごとに意味のある量だけ**を **1 行/monitor** で出す形に置き換える。

## 2. スコープ

- **やる**:
  - `StepMonitor` (main.cpp): `monitorInterval` ごとに 1 行。壁時計 `ms/step`・経過・ETA、残差要約
    (`rms_ro` + 低下桁数最小の列)、NaN フラグ、unsteady のみ `t`/`dt`/`maxCFL`。
  - 旧 `Step :` / `max cfl :` / `dt :` / `Stage : n` / `Inner : n` 行の廃止。`setDT_d_wrapper` は印字せず
    `cfg.monitorCflMax` に格納。定常では `printCfl=false` で host 同期そのものを省く。
  - 起動時にモード要約を 1 回印字 (`mode`, `cfl_pseudo`, `implicitRelax`, `nStepInner`, `monitorInterval`, `nStepOuter`)。
  - 終了時 `Time = ... s` を `clock()` (CPU 時間) から壁時計に変更。
- **やらない** (残作業表へ):
  - EOS 床 (`pMin`/`roMin`/`tMin`) ヒット数・`min P`/`max T` の印字 (reduction 追加が要る)。
  - 境界質量流量 (入口/出口/収支) の印字。
  - `residual_history.csv` の列・行構成の変更 (不変を保つ)。

## 3. 関連 docs と前提

- モード判定の軸は `time.unsteady`: 0 なら陽/陰を問わず局所擬似時間 (`setDT_d.cu` の `if (dualTime==1 or unsteady==0)` 分岐)、
  1 なら `cfg.dt` が物理 dt。詳細は methods §8.5 の表。
- 残差値は `ResidualCsvLogger` が flush 時に保持する最新 `outer_end` 行を使う (CSV と同一値・追加 reduction 無し)。
- 既存の `monitorInterval` の意味 (CSV flush 間隔 = console 出力間隔) は変えない。

## 4. 設計方針

- `ResidualCsvLogger` に `latestSummary()` を追加: `flushDevice()` / `writeRowImmediate()` で step 0 `outer_begin` の
  rms を基準 `rms0` として保持し、最新 `outer_end` の rms と step を保持する。低下桁数 = `log10(rms0/rms)`。
- `StepMonitor::report(iStep)` は `iStep % monitorInterval == 0` の step 末尾で `residual_logger.flush()` → `latestSummary()`
  → 1 行印字。flush は monitor step でいずれ起きる D2H を前倒しするだけ (行構成不変)。
- `setDT_d_wrapper(…, printCfl)` の `printCfl` は「host 読みして `cfg.monitorCflMax` に格納」の意味に変える (名前は互換のため維持)。
  呼び出し側 3 箇所 (steady implicit / dual-time / explicit) で `unsteady==1 && onMonitor` のときだけ true。
- 壁時計は `std::chrono::steady_clock`。monitor step では flush の D2H 同期が入るので区間平均 `ms/step` は GPU 実働を含む。

## 5. 実装ステップ

1. `solver_density_cuda/input/solverConfig.hpp`: `flow_float monitorCflMax = -1` を追加。
2. `solver_density_cuda/cuda_forge/setDT_d.cu`: printf 2 行を削除し `cfg.monitorCflMax = cfl_max` に置換。
3. `solver_density_cuda/main.cpp`: `ResidualCsvLogger::latestSummary()`、`StepMonitor` クラス、`advanceOneStep` 冒頭の
   `Step :` 印字と `advanceExplicitRK` の `Stage :` 印字を削除、3 経路の `printCfl` 条件を `unsteady==1 && onMonitor` に、
   `main()` の起動時モード要約と終了時壁時計。
4. `methods/architecture/overview.md` §8.5 (済)、`procedures/solver-settings.md` / `calculation-workflow.md` の
   `max cfl` 記述を新形式に合わせる。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | ~~実装 + 3 モード検証~~ | 決着 (2026-09-09, §9): 4 run (定常陰/非定常陽 固定 dt/非定常陽 CFL 適応/dual-time) で行形式確認・CSV は旧×2/新×2 の列別比較で run-to-run ノイズ内 |
| 1b | ~~codex result Major 1〜3 / Minor 4〜6~~ | 決着 (2026-09-09, §6.1): dt/maxCFL の刻み対応 (`monitorCflDt`, `dt_next`)、再現可能な CSV 比較 (`csv_compare.txt`)、同条件壁時計 (`wall_compare.txt`)、`from0` 表示、文書訂正を全て採用 |
| 1c | 定常 `dtControl: 1` の adapt 読み | 定常では `cfg.dt` の適応は解に効かないのに monitor step ごとに `thrust::max_element` の host 読みが残る (§6.1 Minor 5)。`unsteady==0` では adaptDt 自体を切ってよいか (dt_local に cfg.dt が入らないことの再確認) を別途判断 |
| 2 | EOS 床ヒット数 / min P / max T | monitor step だけ fused reduction で `floorHits N  minP  maxT` を行に追加。発散の早期指標 ([detectnan-reports-symptom-not-cause] の指紋) |
| 3 | 境界質量流量 | 入口/出口 bcond ごとの ṁ と収支をモニタ行へ。`check_quasisteady.py` の前段の目安 |
| 4 | ~~docs 追従~~ | 決着 (2026-09-09): `procedures/solver-settings.md` / `calculation-workflow.md` の `max cfl` 記述を更新済み |

## 6. 検証

- **ビルド**: `solver_density_cuda/tools/build_native_wsl.sh` (native)。
- **検証ケース** (いずれも新規 run ディレクトリに複製して短く回す):
  - 定常陰解法: `case/16.nozzle_wys/run_0330_monitor_steady/` (run_0317 複製, 200 step)
  - 非定常陽解法 (固定 dt): `case/05.sod_shock_tube/run_0013_monitor_explicit/` (run_0012 複製)
  - 非定常陽解法 (CFL 適応 `deltaT.control: 1`, `monitorInterval: 5`): `case/05.sod_shock_tube/run_0014_monitor_explicit_adaptdt/`
  - dual-time: `case/09.Taylor-Green/run_0051_monitor_dualtime/` (run_0050 複製, 短縮)
- **判定基準**: (a) 各モードで表の通りの列が出る (定常に `t/dt/maxCFL` が出ない、unsteady には出る。CFL 適応の
  陽解法では `dt` (使用) と `maxCFL` が同じ刻みで対になり `dt_next` が別表示)。
  (b) `residual_history.csv` の旧 vs 新の差が run-to-run ノイズ内: 旧×2・新×2 を `_baseline/ _baseline2/ (run 本体) _new2/`
  に保存し、列ごとに `max|x−y| / max(|x|,|y|)_列最大` を旧-旧・新-新・旧-新で並べ、旧-新 ≤ 3×max(旧-旧, 新-新) を合格とする
  (`csv_compare.txt`)。node でもビット一致は要求しない (残差 reduction の atomicAdd と場の run-to-run 揺らぎが元からある。
  改定 2026-09-09, codex Major 2)。
  (c) 同一区間の壁時計: `/usr/bin/time %e` でプロセス全体を旧/新交互に 3 回ずつ計測し中央値で比較 (`wall_compare.txt`)。
  - 注: これらの run はログ形式の回帰試験であり、収束・定常量の主張はしない (`check_convergence.py` は全て NOT CONVERGED
    で当然。200/100 step の短縮 run)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan 免除 | `2026-09-09` | — | — | 起票・実装・検証を同日同セッションで完了した後に AGENTS.md の codex レビュー節が追加された (2026-09-09 以前起票・実装進行中の免除条項に該当)。result 段は実施 |
| result | `2026-09-09` | [2026-09-09-architecture-runtime-monitor-line-result.md](../../notes/reviews/2026-09-09-architecture-runtime-monitor-line-result.md) | NO-GO, C0/M3/m3 | 全件採用。M1 dt/maxCFL 刻み不一致→`cfg.monitorCflDt` を適応前に対で保存し `dt_next` 別表示 + run_0014 (CFL 適応) で検証。M2 CSV 比較→旧×2/新×2 を保存し列別比較 (§6(b) 改定)。M3 壁時計→`/usr/bin/time` 同一区間 3 回 (§6(c))。m4 初期ゼロ列→`from0` 表示。m5 methods §8.5 の同期記述と dual-time コメントを訂正。m6 ω 最大値 4.1e4 に訂正。副次: レビュー中に露見した陽解法 `totalTime` の加算順バグ (適応後 dt を加算) を修正。再レビューは回していない (指摘は機械的に確認可能な範囲) |

## 7. 影響範囲

- `main.cpp`, `cuda_forge/setDT_d.cu`, `input/solverConfig.hpp`
- `forge_run.log` を grep するスクリプト: `case/39.periodic_hills/run_lineimp2_arms.sh` (`Time = ` は維持するので影響なし)
- docs: methods §8.4/§8.5 (済)、procedures 2 箇所 (残作業 #4)

## 8. 完了条件

- [x] 関連 `methods/<area>/` の現在仕様を更新済み
- [x] 実装・検証完了 (本計画の §6 を満たす)
- [x] 本計画の `status` を `done` に変更し、§9 に変更ログを記載
- [x] ファイルを `plans/active/` → `plans/accepted/` へ移動
- [x] [`plans/README.md`](../README.md) の一覧を同期

## 9. 変更ログ

- `2026-09-09` — 初稿。対話で決めた方針 (モード別の時間量・1 行/monitor・残差要約・壁時計) を反映。
- `2026-09-09` — 実装・検証完了 (`main.cpp` `StepMonitor` / `ResidualCsvLogger::latestSummary`, `setDT_d.cu` の printf → `cfg.monitorCflMax`/`monitorCflDt`)。
  codex result レビュー (§6.1, NO-GO → 全指摘採用) 後の最終検証 (native RTX 3060, 旧バイナリ = HEAD 1a9d75c2 を
  `FORGE_CUDA_ARCHITECTURES=86` で別ビルド。各 run に旧×2 `_baseline/ _baseline2/`・新×2 `run 本体/ _new2/` を保存):

  | run | モード | 行形式 | CSV 列別 旧-新 最大 / ノイズ床 max(旧-旧, 新-新) | log 行数 旧→新 |
  | --- | --- | --- | --- | --- |
  | `case/16.nozzle_wys/run_0330_monitor_steady` (200 step) | 定常陰解法 node TP+SST | `t/dt/maxCFL` 無し | 1.6e-4 / 1.5e-4 (rms_roOmega 列) | 907→307 |
  | `case/05.sod_shock_tube/run_0013_monitor_explicit` (200 step) | 非定常陽解法 RK3 cell 固定 dt | `t dt maxCFL 0.07`, `from0 rms_roUy` | 1.5e-1 / 1.5e-1 (rms_roUy ≈1e-6 のゼロ近傍列; rms_ro は 1e-6 / 9e-7) | 1512→312 |
  | `case/05.sod_shock_tube/run_0014_monitor_explicit_adaptdt` (200 step, monitorInterval 5) | 非定常陽解法 CFL 適応 | `t 1.0969e-06 dt 1.10e-06 maxCFL 0.20 dt_next 1.09e-06` (t = 使用 dt の累積) | 2.1e-1 / 3.3e-1 (同上ゼロ近傍列) | — |
  | `case/09.Taylor-Green/run_0051_monitor_dualtime` (100 step) | dual-time node SST | `t dt maxCFL 0.13` (旧は 20 sub-iter × 2 行/step) | 1.2e-4 / 1.3e-4 | 4330→230 |

  4 run とも `csv_compare.txt` の VERDICT は「旧-新 ≤ 3×ノイズ床」で合格 = 数値変化なし (コード変更は印字・host 読みの
  間引き・`totalTime` 加算順のみ)。壁時計 (`wall_compare.txt`, run_0330 プロセス全体, 旧/新交互 3 回): 旧 1.13/1.19/1.19 s,
  新 1.13/1.21/1.24 s (中央値 1.19 vs 1.21 s, +1.7 %, 反復ばらつき 0.06〜0.11 s の内)。`monitorInterval: 1` での
  per-step ostringstream/cout が上限で、`monitorInterval>1` では消える。
  副産物: run_0330 の `worst rms_roOmega (+2.4)` (step 2 の行から表示) は run_0317 にも元からある restart 直後の SST ω
  残差過渡 (outer_end 24.6→最大 4.1e4→1000 step で 331) で、モニタ行がこれを冒頭から見せる。
  副次修正: 陽解法 `advanceExplicitRK` で `totalTime += dt` を setDT (適応) の前に移動 (旧は次 step 用 dt を加算していた。
  固定 dt では無影響、CFL 適応では t が step ごとに (dt_next − dt) ずれていた)。
  ビルドの注意: 新規 build dir は `FORGE_CUDA_ARCHITECTURES=86` を付けないと `compute_52` 既定で
  `residualMonitor_d.cu` の `atomicAdd(double)` がコンパイルエラーになる。
