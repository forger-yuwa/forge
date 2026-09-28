# case/60 — Deveikis & Hunt TN D-7275 大型校正パネル (M7, 8-ft HTST) = **T4-0b**

計画: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md) §4.12 / §6 G15 / §5.1 #61
精読メモ: [`notes/investigations/tp1187-test-system-reading.md`](../../notes/investigations/tp1187-test-system-reading.md)

TP-1187 (case/56) の乱流分母 q_FP の出所。フェンス付きパネルホルダー・鋭い前縁・0.24 cm 球トリップ (前縁から 12.7 cm)。
T4-0a (case/59、Cary) の比較方法が固まってから計算する。

## 実験データ

- `digitize_d7275_fig20_test26.json`: 試験 26 (Tt 1867 K、Pt 18.06 MPa、α 0°、M∞ 6.64、Re∞ 4.757×10⁶/m) の中心線 St*_l (Fig 20 の □、10 点)、
  Table II/III の試験 26 行、換算式と読み取り幅 (St ±1.3 %、x 対応 ±3 %)。

## 準備 (2026-09-27)

- 自由流の再構成: 燃焼ガスモデル (case/50 の CombustionProducts、φ 0.711) で M∞ 6.64・p∞ 2117 Pa・Tt 1867 K から
  T∞ 227.7 K (表 228.3)、Re∞ 4.71×10⁶/m (表 4.757、−1 %)、(ρVcp)∞ 67.95 kW/m²K (表 68.89、−1.4 %)。
- メッシュ `mesh/fp_d7275_y3` (case/56 の `gen_mesh.py`、前縁から 2.6 m、y1 3 µm、H 0.5 m、37.6 万節点): `VERDICT: PASS` (AR 最大 854)。
- 生成器 `gen_runs.py` (段の引き継ぎは restart_field.py、FP64 アキュムレータ、stage manifest)。
- トリップ (12.7 cm) は 2D で入れられないので、乱流の起点は前縁と、同じ解を後処理で 12.7 cm ずらした場合の 2 通りで示す (有効長の感度、上下界ではない)。

## 計算 run 一覧

| run | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
|---|---|---|---|
| `run_0001_t26` | 試験 26 (α 0°)。燃焼ガス、壁 300 K、前縁から乱流 (SST、Tu 0.5 %、μt/μ 10)、`fp_d7275_y3`、段階起動 lam→soft→mid→2 次ランプ→本段 60k (cfl 1.5)、**全域 FP64** (`~/forge56-double`)、段の引き継ぎは restart_field | 本段 60k 完走、NaN 0。判定区間 (`--segment`) で `NOT CONVERGED (still converging)` (全列 1.6–3.3 dec、下降中)。壁解像: y1+>1 は前縁の 2 列のみ (最大 3.54、長さ重み 0.03 %)、パネル域 1.07–2.46 m は 0.26–0.28。暫定 R (St*_l、前縁起点) = 1.27–1.43 (平均 1.31)、トリップ起点で平均 1.32、x 対応 ±3 % で 1.307–1.314 (**未収束につき解釈しない**) | 延長元 |
| `run_0002_t26_ext` | run_0001 から同一設定 +60k (restart_field、倍精度のまま 12 量ビット一致)。1 回目は `species_db.yaml` を複製し忘れて forge が起動せず → スクリプトを直して空の run を消して再投入 | 延長区間 `NOT CONVERGED (stalled/plateau)`、run_0001 の判定区間と接続しても 2.2–3.9 dec で `still converging`/プラトー (**全域 FP64 なので丸め床ではない**)。St **ALL STEADY**。**R (St*_l、前縁起点) = 1.266–1.442、平均 1.313** (トリップ起点 1.328、x ±3 % で 1.313–1.319)。forge の q_w: パネル平均 89.0 kW/m²、位置 II 87.8 kW/m² (Tw 300 K) | active (**T4-0b 試験 26 の結果、収束ゲート未達**) |
| `run_0003_t26_probe` | run_0002 の最終場から同一設定 4000 step、1000 毎に出力 (残差の出どころを探す短い継続) | 2000→4000 step の場の相対変化は**下流の slip バッファ (x > 2.55 m、平板の後ろの後流) だけ** (ρ 7e-3、P 8e-3、roK 0.15)。比較域の壁近傍 (wd < 1 mm) は 2e-6 以下、前縁 1e-6 → **残差プラトーは比較域の外に局在** | ref (診断) |
| `run_0004_t26_wdA` / `run_0005_t26_wdB` | **T4-0b-WD** (事前登録 commit fa55f89a): 下流 slip を physID 7 に分けたメッシュで、A = 壁距離そのまま / B = 変換時 `wallDistExtraPhysIDs: [7]` (x ≥ 2.6 m の壁距離だけ変わる)。run_0002 の最終場から restart_field --keep-src-dtype (化学種を含む 12 量ビット一致)、同一設定 20k、1000 毎。**1 回目は新メッシュに roY が無く 7 量しか移らなかったので起動 1 分で止めて作り直した** | **1 回目は AWS ディスク満杯で step 1900 に Write failed → 整理して 2 回目を投入** (出力 2000 毎)。2 回目の結果: 両腕 `NOT CONVERGED (stalled/plateau)`。B の末尾残差は A の 0.53–0.73 倍で **1/10 に届かない → 事前登録の読み = 壁距離単独主因説を退ける** (codex 再読: 『壁距離変更だけでは解消しないが寄与はある』)。熱流束の差は比較点で 1.9e-7・パネル域で 2.1e-6 → **比較量は下流の変更に鈍い**。St 両腕 ALL STEADY、R 平均 1.313 | ref (A/B) |

> **AWS ディスク整理 (2026-09-27)**: run_0001 の中間スナップショット、run_0002 の res_0 (run_0001 最終場の複製)、run_0003_t26_probe の場のファイル (結論は上の行に記録済み) を削除。

> **AWS ディスク圧縮 (2026-09-27、ユーザ指示「貴方の関連のディスク容量圧縮の為、ファイルをガンガン削除していい」)**: 判定の根拠 (壁・表面出力 `res_wall_*`/`res_plate_*`/`res_gap_*`、残差履歴、VERDICT、時系列 CSV、設定) と最新の継続用最終場 (Cary: run_0005・0014–0024、D-7275: run_0002・run_0005) だけを残し、それ以外の全場スナップショット・run 内の `mesh.h5`・旧 3D メッシュ `gap3d_prod.h5` を削除した (計約 13 GB)。削除した場が要る再解析は再計算になる。
| `run_0006_t26_resmap` | run_0002 から 2000 step、`FORGE_OUT_RESIDUALS=1` で残差場を狙ったが、`output.level: 1` のフィルタで残差場が落ちて取れなかった | — | 破棄 |
| `run_0007_t26_resmap2` | 同上、`extraFields` に `res_ro/res_roUx/res_roUy/res_roe` を明示 | **残差場の二乗和の 99.99–100 % が平板の後ろ (x ≥ 2.6 m、節点の 6.3 %) にある**。最大は出口∩上端 slip の角 (x≈2.79 m、y≈0.47–0.50 m) と出口∩下流 slip の角 (y=0)。比較域 1.07–2.46 m と境界層 (wd < 3 cm) は 0.00 % → **残差プラトーは出口の角に局在** (状態差の測定 run_0003 と整合、今度は残差そのもの) | ref (診断) |
| `run_0008_t26_outflow` | **T4-0b-OUT** (事前登録 commit 57b94f2e): run_0002 の最終場から出口だけ `outlet_statPress` → `outflow` (現行レシピの超音速 node 出口)。20k、2000 毎、残差場も出力 | `NOT CONVERGED (stalled/plateau)`。末尾残差は同じ長さの対照 run_0004 (statPress) と 1.00 倍 → **出口の種類を替えても改善しなかった** (超音速流出では両者とも内部状態の外挿で実質同じ作用なので、出口閉包そのものは除外していない — codex)。残差は依然 100 % が x ≥ 2.6 m。熱流束の差は比較点 3.9e-6・パネル域 4.9e-6 → **比較量は鈍い**。St ALL STEADY、R 平均 1.313 | ref (A/B) |
| `run_0009_t26_relaxA` / `run_0010_t26_relaxB` | **T4-0b-RELAX** (事前登録 commit cd395f01): run_0002 の最終場から implicitRelax 1.0 / 0.5 だけを変える。20k、残差場 (k/ω 含む) を出力 | 両方 `NOT CONVERGED (stalled/plateau)`。末尾残差 B/A = 0.86–1.19 (1/10 に遠い)、出口近傍 (x ≥ 2.7 m) の 2000 step 状態差は B が A の約 0.5 倍 → **事前登録の読み = 判別保留** (中間的。振れ幅の半減は更新量半減の見かけの可能性)。A の roK は落ち着き条件も未達 | ref (A/B) |
| `run_0011_t26_ext` | **T4-0b-EXT** (ユーザ指示「出口側を slip にして延長したらいいよ」、事前登録 commit 820d20a5): 下流 slip を 0.2 → 1.0 m (出口 x 3.6 m) に延長したメッシュ `fp_d7275_y3_ext` (上流・平板の分割は元と一致、42.5 万節点、`SOFT-PASS` AR 最大 ≈ 1043)。run_0002 の最終場から interp_field (化学種を含む 12 量、ΣY = 1)、全域 FP64、本段設定 40k、5000 毎、残差場も出力 | `NOT CONVERGED (stalled/plateau)`、未達残差は同じ長さの対照の **1.9–14.6 倍に悪化** → 事前登録の読み = 延長だけでは解消しない。**残差の 96–98 % は延長した slip 区間の内部 (2.6–3.5 m)、出口付近は 2–4 %** → 源は出口の角でなく平板の後ろの slip 後流。St ALL STEADY、R 平均 1.313 (延長前と同じ) | ref |
| `run_0012_t26_wall` | **T4-0b-WALL** (ユーザ承認「そのやり方でOK」、事前登録 commit 679eb81e): 平板 (等温壁 300 K) を出口 x = 2.8 m まで延ばし slip 後流を無くす。メッシュ `fp_d7275_y3_wall` (元と形状・分割が同一、下流 2.6–2.8 m の底辺だけ slip → plate、`VERDICT: PASS`)。run_0002 の最終場から interp_field (化学種を含む 12 量、ΣY = 1、NaN 0)、全域 FP64、本段設定 40k、5000 毎、残差場 (k/ω 含む) も出力 | `NOT CONVERGED (stalled/plateau)` (流れ 2.0–2.7 dec 横ばい、k/ω 6.5–7.4 dec 下降中)。step 16–20k の残差は対照 run_0009 の 0.26–0.29 倍 (roUy 0.94、roOmega 2.8) → **読み = not_resolved** (1/10 に届かず)。**slip 後流の残差は消え、残りの 99.7–99.9 % は上端 slip の出口直前 (y 0.45–0.50、x 2.70–2.81 m)** — 前縁衝撃波 (約 9.5°) が x = 2.8 m で y ≈ 0.47 に達し上端∩出口の角に掛かる。St 比較 10 点は run_0002 と表示桁で一致 (R_A 平均 1.313)、St ALL STEADY | active |
| `run_0013_t26_h05` / `run_0014_t26_h079` | **T4-0b-H** (codex diagnose 2026-09-27、事前登録 commit 2904655a): 上端 slip の位置だけを変える A/B。A = `fp_d7275_y3_wall` (H 0.5 m、run_0012 最終場から restart_field --keep-src-dtype、12 量ビット一致) / B = `fp_d7275_y3_wall_top8` (A の全節点を座標完全一致で含み上に同じ等比で 8 層、H 0.794506 m、39.4 万節点、`PASS`。`tools/stack_init.py` で共通節点ビット一致コピー、追加節点は入口自由流)。各 20k、1000 毎、残差場 (k/ω 含む)、旧上端付近の点プローブ (T・P・U、10 step 毎) | 両腕 `NOT CONVERGED (stalled/plateau)`。**旧上端域の √Σres² B/A = 7e-4–5e-3 (全 6 成分 1/10 以下)**、旧上端直下の P 振幅 A 3.8e-5–1.5e-4 → B は出力分解能以下、B の残りの流れ残差は自由流域に薄く分布 (同水準の移動なし、全域 √Σres² B/A ro 0.032・roe 0.10) → 機械的な読み = 上端近接説を支持 (解釈は諮問中)。rms_roOmega は両腕 0.6 で横ばい (壁第一内部節点列に散在、上端と無関係)。比較域 q の B/A−1 ≤ 1.5e-6、St ALL STEADY、R_A 平均 1.313 | active |
| `run_0015_t26_tf1` / `run_0016_t26_tf0` → `run_0017_t26_tf1_ext` / `run_0018_t26_tf0_ext` | **T4-0b-TF** (codex diagnose 2026-09-27、事前登録 commit 4f6cf102): `physProp.thermoFloat` 1 / 0 だけを変える。run_0014 最終場から restart_field --keep-src-dtype (12 量ビット一致)、`fp_d7275_y3_wall_top8`、各 5k + 化学種残差場・ω 項別収支 (`FORGE_OMEGA_BUDGET=1`)・step 4000–5000 を 100 毎。5k で k/ω の窓間中央値が 5 % 超動いたので事前登録どおり各 +5k 延長 (0017/0018、restart_field) | 5k 時点 (暫定): 自由流域の流れ・化学種残差 B/A 0.92–1.06 (不変)、第一内部列は roUy 0.028・roe 0.035・ro 0.10・roOmega 0.38。第一内部列 ω 収支は trans ≈ dest ≈ 6.5e8 の相殺残りが A ±1e2 → B ±3e-2。**10k の結果**: 両腕 `NOT CONVERGED (stalled/plateau)`。自由流域の √Σres² B/A 0.91–1.05 → **float 熱力学経路は自由流域の床の主因でない (棄却)**。第一内部列は roUy 0.028・roe 0.042・ro 0.10 (float 経路の寄与大)、ω 0.46・k 0.84 は下降中で判別保留。第一内部列 ω 収支の相殺残り A ±1e2 → B ±1e-3–3。新たな観測: 自由流域の |res_ro|/(ρU√vol) 中央値 1.2e-8 (float32 水準)、変換メッシュの幾何は float32 | active |
| `run_0019_t26_hf_f` / `run_0020_t26_hf_d` | **T4-0b-HF** (codex diagnose 2026-09-27、事前登録 commit 19d1553a): 両腕 thermoFloat 0、SLAU の面エンタルピー評価だけを float (A、`~/forge56-double`) / double (B、`~/forge56-hface` = 同じ 4687c3c + `tools/hface_double.patch` の 2 行) に切り替える。run_0018 最終場から restart_field --keep-src-dtype (12 量ビット一致)、各 5k、場は 3000–5000 を 100 毎、点プローブ 1 step 毎。B の窓間変化が 5 % 超 (−4〜−40 %) なので事前登録どおり各 +5k 延長 (`run_0021_t26_hf_f_ext` / `run_0022_t26_hf_d_ext`、restart_field) | 5k 時点 (暫定): **自由流域の √Σres² B/A = 3.4e-4〜1.2e-3 (流れ・化学種)**、ω 4.9e-3、k 0.12 (時系列 21 枚の末尾窓中央値)。B の全域履歴は 3.3 dec 低下して下降中 (A は横ばい)。第一内部列は ω 0.53・k 0.72・roUy 0.11・ro 0.30 (両腕とも下降中)。**10k**: 自由流域 B/A は同じ ~1e-3 (roe 3.5e-4) だが B の窓間変化 −13〜−30 % (再開直後の残差の跳ねからの戻りを含む) → 字面では判別保留。両腕 `NOT CONVERGED` (B は流れ falling・ω plateau、A は全列横ばい)。第一内部列 ω は B/A 0.58 (42 % 低下、停滞は解消せず)。比較域 q の差 ≤ 4.4e-7、R_A 平均 1.313。**読み = 判別保留** (codex 2026-09-27: 残差比は満たすが窓間変化 5 % 超・振幅は測定不能)。観測としては面エンタルピー精度への大きな感度。stop_rule により追跡はここで区切る | ref (A/B) |
| `run_0023_t28` | **試験 28** (α 9.8°、Tt 1817 K、M∞ 6.60、p∞ 2.144 kPa。事前登録 T4-0b-T28、commit dffa3516): 入口は forge の燃焼ガスで解いた斜め衝撃波 (β 16.53°) の背後の局所状態を平板に平行に (T 348.65 K、U 1923.95 m/s、p 8451 Pa、ρ 0.0819 kg/m³、M 5.11、Tu 0.5 %・μt/μ 10)。`fp_d7275_y3_wall_top8` (壁を出口まで・上端 0.7945 m)、`gen_runs.py --test 28` の段階起動 + 本段 60k、全域 FP64。実験点は Fig 19(c) △ 11 点を Table I の熱電対位置に対応 | 本段 `PASS (converged)` (4.0–4.6 dec) だが **St は下流 8 点が DRIFTING/TRANSIENT-UNSETTLED** (層流水準から上昇中、最下流は 50k→60k でも +2.4 %)、**壁解像 FAIL** (y1+ 平均 0.89、>1 が 14.8 %)。暫定 R28 = 1.377 (読みに使わない) | 延長元 |
| `run_0024_t28_ext` | run_0023 から同一設定 +60k (restart_field --keep-src-dtype、12 量ビット一致)。St の準定常の確認 | 実行中 | active |
| `run_0025_t28_y2` | 第一層 2 µm のメッシュ `fp_d7275_y2_wall_top8` (41.1 万節点、AR 最大 1281 = 構造格子の境界層セルなので **AR 緩和 (≤5000)** で PASS) へ run_0023 最終場を interp_field (12 量)、60k。壁解像の確認 | 実行中 | active |

- 2026-09-27 AWS 整理: run_0006〜0011・0015〜0020 の中間の全場スナップショットを削除 (最終場・壁出力・残差履歴・VERDICT・点プローブは保持)。T4-0b-HF の領域別時系列は削除前に `tools/hf_ab.py` で集計済み。

## Eckert 参照温度法との比較 (2026-09-27)

`tools/eckert_compare.py` (論文 p.14 の乱流の式 St*ₗ = 0.0296·Pr*^(−2/3)·R*ₗ^(−1/5)、Table III の換算係数、Pr*(600 K) = 0.6905) で `run_0002_t26_ext` の壁熱流束 (step 60000) を実験点と同じ座標に置いた。比較 10 点で **forge/Eckert 0.969–1.009 (平均 0.991)**、**実験/Eckert 0.700–0.777 (平均 0.755)**、forge/実験 1.313。感度: 乱流起点をトリップにすると forge/Eckert 0.976・実験/Eckert 0.744、T* を Tw 300 K から組み直すと forge/Eckert 1.001。forge/Eckert は前縁近くで 1 を下回る (x 0.2 m で 0.89)。まとめのページ: https://claude.ai/artifact/W6RzaLpfYWGopBEdLd8aDx (非公開)。

**訂正 (同日)**: 上の forge/Eckert 0.991 は Eckert に比較用の Python 気体モデル (case/56 `gas_htst`) の Pr*(600 K) = 0.6905 を入れた値。forge 自身の輸送物性 (kinetic theory、`run_0002_t26_ext/res_60000.h5` の `vis_lam`・`thermCond` から T* ±2 K の中央値: μ* 2.95e-5 Pa·s、Pr* 0.7446、`forge_props_run0002_Tstar.json`) で Eckert の q を組むと **forge/Eckert 1.060 (1.036–1.080)**、T* 594 K で 1.054。同じ Python 物性で q を直接組むと 1.020。→ forge/Eckert は Eckert に入れる物性で 0.99〜1.06。実験/Eckert (0.755) と forge/実験 (1.313) は同じ換算係数で割るので変わらない。初版ページは物性の表を持たず、この違いに気づかなかった (ユーザ指摘で `forge-report` skill に物性の表を追加)。ページは v2 で解析条件・物性・判定ゲートを追加して差し替え。
- 2026-09-28 AWS 整理 (ユーザ指示「終わったやつはガンガン消していい」): run_0001〜0022 (0002 の res_60000 は保持) の全場スナップショット `res_<n>.h5` と run 内の `mesh.h5` の複製を削除、使い終わったメッシュも削除。残差履歴・ログは gzip (`*.csv.gz`)、壁出力・VERDICT・時系列 CSV は保持。
