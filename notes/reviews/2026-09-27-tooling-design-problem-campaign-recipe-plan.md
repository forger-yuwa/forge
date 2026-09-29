# codex レビュー: tooling-design-problem-campaign-recipe (plan)

- **plan**: [`plans/active/tooling-design-problem-campaign-recipe.md`](../../plans/active/tooling-design-problem-campaign-recipe.md)
- **stage**: `plan`
- **date**: 2026-09-27
- **commit**: `e79bad40` (feature/gap-heating-precision)
- **codex**: effort `high`, 5.8 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2
- **extra**: `notes/reviews/2026-09-27-design-problem-schema-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
課題の同定と「機種別 runner＋小さな共通契約」という方向は妥当です。  
ただし、熱力学の修正方針、δ* 反復の実行契約、評価量・合否の定義を以下のように直してから実装すべきです。

1. **Major — §4.6 の γ*/cp(Tt) は、整合した CPG 定数対になっていない。NS 初期化には別の熱力学的不整合もある。**

   `Problem.R_gas` は R = cp(γ−1)/γ を計算します（[probdef.py:33](/home/sano/work/forge/design/forge_design/probdef.py:33)）。va3 のガスモデルを実行して確認すると、モデルの R は **285.2704 J/(kg·K)**、提案された γ* と cp(Tt) の組から得る R は **292.5894 J/(kg·K)**、**2.57% 過大**です。異なる温度で評価した定数を組み合わせています。

   さらに NS の ω 床は、`roe` から **CPG の式で温度を逆算**しています（[runner_axismach.py:820](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:820)）。TP の IC は `thermoHrefTemp` を反映した内部エネルギーなので、この逆算は不整合です（[ic.py:69](/home/sano/work/forge/design/forge_design/evaluate/ic.py:69)）。`run_0509` の Euler 参照場にこの初期化式を適用した読取専用検算では、実際の T が **282.03–1160.84 K** なのに、逆算値は **50–685.33 K**、**57.18% のノードが 50 K 下限に張り付きました**。これは初期化式の検算であり、NS 計算を実行した結果ではありません。

   **対案:** TP の R は組成から直接取得し、温度は組成・エネルギー基準と整合した EOS で復元する。CPG 近似には同一参照温度の γ/cp、または R と一方の定数を使う。§5.1 #6 は「旧結果不変」だけでなく、**熱力学整合性の試験と、不具合修正に伴う結果変化の評価**を分離してください。

2. **Major — `ExecutionPlan` が、事前に決まる契約と実行後に確定する成果物を区別していない。**

   [plan:89](/home/sano/work/forge/plans/active/tooling-design-problem-campaign-recipe.md:89) は事前解決・必須入力の充足確認・計画だけによる実行を要求します。一方、δ* の次の壁は前 pass の NS 場から生成されます（[deltastar_loop.py:135](/home/sano/work/forge/design/forge_design/feedback/deltastar_loop.py:135)）。campaign 開始時には、その場も補正壁も存在しません。ここを未定義にすると、事前解決が後編集依存に戻ります。

   また、現在の初期場指定は run ディレクトリから最大番号の `res_*.h5` を選びます（[runner_axismach.py:805](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:805)）。同じ入力パスでも継続計算後には別の場を読むため、計画の決定性と再現性を保証できません。

   **対案:** campaign の事前検証では「将来生成される参照」と「欠落入力」を区別し、recipe が各 pass の直前に具体的な実行計画を確定する契約にする。参照には candidate・形状・メッシュ・輸送種配置・エネルギー基準・**特定 snapshot と内容ハッシュ**を含め、不適合を拒否する。各段の設定と履歴は既存 `StageManifest` に保存し、最終 YAML だけでなく**全段**を計画と照合してください。汎用 DAG は不要です。

3. **Major — V2 は qualification が失敗しても合格でき、補正後の順位付けを検証していない。**

   [plan:151](/home/sano/work/forge/plans/active/tooling-design-problem-campaign-recipe.md:151) の「ゲートを満たすか、満たさなければ pass 1 まで回して結果を記録」は、**ゲート不合格の記録だけで V2 を通せる**条件です。最終順位の使用量、失敗候補の除外、全候補失敗時の扱いもありません。

   既存 `run_pass` は solver の戻り値を表示した後に `collect` と次の δ* 抽出へ進み、収束検査も `check=False` です（[deltastar_loop.py:156](/home/sano/work/forge/design/forge_design/feedback/deltastar_loop.py:156)）。したがって既存関数を呼ぶだけでは qualification の合否は成立しません。

   **対案:** 生の判定とは別に、`accepted / rejected / incomplete` を定義する。NS/Euler 比と出口コア M は**両系列の定常性**を確認し、反復上限でも不合格なら最終順位から除外する。全候補不合格なら「勝者なし」とする。V2 には補正後の目的量・制約の再評価を追加し、CFD 不要の試験で「Euler 順位が NS で逆転」「判定欠落」「全候補失敗」を検証してください。

4. **Major — V0/V2 の「軸 M 目標差」「軸 M 出口」は、比較元と runner で定義が異なる。**

   runner は軸方向 **200 点へ補間**して最大誤差を計算し、`M_axis_exit` を **x_E** で評価します（[runner_axismach.py:572](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:572)）。case 側は実際の軸ノード上で最大誤差を取り、出口側の **x_max−2r_t に近い断面**を使います（[lumpX_series_csv.py:13](/home/sano/work/forge/case/44.vitiated_air_wt/lumpX_series_csv.py:13)、[同:30](/home/sano/work/forge/case/44.vitiated_air_wt/lumpX_series_csv.py:30)）。

   `case/44.vitiated_air_wt/run_0509_va3_M4.19_Lc8_dry_lumpX/res_24000.h5` を読み直した結果です。

   | 量 | case 側の保存系列 | runner の計算式 |
   |---|---:|---:|
   | 軸 M 目標差 max | 0.002129031659 | 0.002120321373 |
   | 軸 M 出口 | 4.190995693 | 4.190284303 |

   前者は同じ場でも **約0.41%** 異なります。これは float32 の再実行ノイズではなく測定定義の差です。

   判定も再確認しました。対象区間は **S2_main、本段 step 0–24000**。`check_convergence.py` は **`NOT CONVERGED (stalled/plateau)`**、保存6量の `check_quasisteady.py --series-csv` は **`OVERALL: ALL STEADY`**。runner 定義で再抽出した2量も同ツールの `classify_series` で **`STEADY`** でした。全7 snapshot の `VALUE/*` に NaN/Inf はありません。run の索引は [case README:783](/home/sano/work/forge/case/44.vitiated_air_wt/README.md:783) です。

   **対案:** 評価量ごとに抽出位置・補間格子・積分重み・正規化・版を固定し、探索、時系列判定、旧 run 比較を同じ抽出関数で行う。その定義で V0 を取り直してください。`run_0509` は**未収束の準定常回帰参照**として明記すべきです。

5. **Major — Euler→NS のメッシュ切替と壁解像が、V1/V2 の必須検証になっていない。**

   現在の `prepare_ns` は `problem.mesh` を優先し、無指定の場合だけ `wall_first_frac=4.5e-5` を採ります（[runner_axismach.py:746](/home/sano/work/forge/design/forge_design/evaluate/runner_axismach.py:746)）。対象 Euler YAML は `wall_first_frac=0.005` を指定しています。単に同じ problem を NS 関数に渡すと、**NS 既定の約111倍の値**を引き継ぎます。slip→no-slip と solverConfig の一致だけでは、この誤配線を検出できません。

   **対案:** recipe ごとのメッシュ条件を V1 の比較対象として数値で固定する。V2 では各補正壁のメッシュ品質、cross-mesh IC の適合性、低 Re SST の局所 y₁⁺ を必須検証にする。品質不合格なら投入停止、壁解像不足なら qualification 不合格とする。`ypls` を根拠にしない規則も既にあります（[AGENTS.md:104](/home/sano/work/forge/AGENTS.md:104)）。

6. **Minor — 1 回の再実行との差を、そのまま許容差にする V0 は不安定。**

   [plan:142](/home/sano/work/forge/plans/active/tooling-design-problem-campaign-recipe.md:142) では、偶然ゼロだった差は許容差ゼロになり、偶然大きかった差は過大な許容になります。また `run_0509` は既に `restart_field.py` を使っているため、「誤った段間移植の修正前結果」ではありません（[run_lumpX_staged.py:96](/home/sano/work/forge/case/44.vitiated_air_wt/run_lumpX_staged.py:96)）。

   **対案:** バイナリ・初期場・実効設定を固定した反復測定から、量別の絶対/相対許容差と安全係数を事前定義する。V0 の取得を数値変更前の作業として残作業表に置き、`run_0509` は「正しい段間移植を行う参照経路」と記述してください。

7. **Minor — §3 の「出口診断は定数 γ の特性線追跡」はコードと一致しない。**

   `gamma` 引数は渡されていますが、`core_radius_traced` は **出力の `sonic`** から M を作り、μ = asin(1/M) で追跡しています。関数内で `gamma` を使用していません（[extract.py:87](/home/sano/work/forge/design/forge_design/metrics/extract.py:87)、[同:102](/home/sano/work/forge/design/forge_design/metrics/extract.py:102)）。出口 M 自体も `sonic` を使います。

   **対案:** §3・§4.6 の用途一覧を訂正する。この診断を新たにガスモデル経由へ置き換える必要はなく、未使用引数の整理に留め、実際に誤っている NS 温度逆算を優先してください。

推奨は、**現方針を維持し、上記契約を確定した axismach の縦断実装を1本作ること**です。実装前の修正順は **①熱力学、②各 pass の実行計画と成果物参照、③合否・最終順位、④評価量と V0、⑤メッシュ・壁解像ゲート**とします。

目的は未解決です。既存ベル driver は固定の `DV_ORDER` を持ち、汎用 campaign は実現していません（[driver.py:48](/home/sano/work/forge/design/forge_design/opt/driver.py:48)）。case/44 の選択は [検証手順](/home/sano/work/forge/procedures/verification/README.md:11) と整合し、今回 cell・周期境界まで CFD 検証を広げる必要はありません。共有 `probdef` の互換性は CPG・semiperfect・`frozen_tp` の読込試験で確認してください。

ファイル変更・新規 CFD 実行は行っていません。以上の修正提案は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 2
