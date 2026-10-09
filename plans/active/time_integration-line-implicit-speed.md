# ライン陰解法 (`lineImplicit`) の 1 step あたりの費用を下げる

## メタ

- **area**: `time_integration`
- **status**: `in_progress`
- **related_docs**:
  - `methods/time_integration/implementation.md` の「line-implicit (`lineImplicit`, 2026-09-02)」と「v2」
- **related_plans**:
  - [`../accepted/time_integration-line-implicit.md`](../accepted/time_integration-line-implicit.md) (block-Thomas の実装、1 ライン 1 スレッド・内部 double)
  - [`../accepted/time_integration-line-implicit-viscous-v2.md`](../accepted/time_integration-line-implicit-viscous-v2.md) (factor/solve の分離)
  - [`time_integration-implicit-thermal-jacobian.md`](time_integration-implicit-thermal-jacobian.md) §6.0 (case/45 で方向別の擬似 dt + キー 5 + 上限 50 を本線の候補にしている)
- **created**: `2026-10-09`
- **owner**: Claude (ユーザ指示 2026-10-09「現状 LINE にすると計算速度がほぼ半減するのは課題だと思うので、改善策がないか調べてほしい」)

## 1. 目的

case/45 の冷却ノズル (FP64 ビルド、AWS g5 の A10G) で、ライン陰解法を入れると 1 step の壁時計が point の約 1.9 倍になる (§4.1)。
解く線形系を変えずに (または精度の選択を明示した上で) この増分を減らし、方向別の擬似 dt で得た step 数の短縮 (§6.0 の run_0224) を壁時計の短縮として取り出す。

## 2. スコープ

- **やる**: 計測 (§4.1)、ライン上の節点の対角の組み直しの省略、block-Thomas の factor/solve のライン内の並列化、
  Thomas の内部精度を `implicitSolvePrecision` に従わせる選択肢、ライン面の近傍行列 K の閉形式での直接計算。
- **やらない**: K に粘性・熱伝導の結合を入れること (線形系が変わる数値の変更。別の議題として扱う)、sweep 回数 (`nStepInner`) の変更、
  streamwise のライン族。

## 3. 関連 docs と前提

- 実装: `solver_density_cuda/cuda_forge/timeIntegration_d.cu` の `implicit_defect_correction_block_d` (sweep)、`lineThomasFactor_d`・`lineThomasSolve_d`、
  `block_dplur_jacobian_d.cuh` の `accumulate_split_jacobian_cf`。
- case/45 の格子: ライン 4719 本 × 最大 121 節点、被覆 570999/570999 CV (100 %)。sweep 5 回/step (`nStepInner 5`)。

## 4. 設計方針 (案。codex plan 段の前)

### 4.1 計測 (2026-10-09、nsys、GPU 専有、同じ新バイナリ 35e498b1…、run_0183 の res_100000 から 500 step)

| run | 設定 | 壁時計 ms/step (nsys の下) | GPU カーネル合計 ms/step (初期化込み) |
| --- | --- | --- | --- |
| `case/45.isobutane_m6_d155/run_0249_prof_point` | point cfl 4 | 18.97 | 17.09 |
| `case/45.isobutane_m6_d155/run_0251_prof_lineonly` | `lineImplicit 1` だけ | 36.37 | 34.32 |
| `case/45.isobutane_m6_d155/run_0250_prof_linedir_tj5_cap50` | ライン + 方向別 + キー 5 + 上限 50 | 37.28 | 35.26 |

ライン解法で増えるカーネル (ms/step): `lineThomasSolve_d` 9.85 (1 回 1.97 ms × 5 sweep)、`lineThomasFactor_d` 5.41 (1 回)、
`implicit_defect_correction_block_d` 4.14 → 6.12 (ラインだけ) → 7.00 (+キー 5)。その他のカーネルは ±0.05 以内。
`ncu` は性能カウンタが管理者権限の設定 (`RmProfilingAdminOnly: 1`) で取れなかった (sudo は使っていない)。

### 4.2 遅さの理由 (コードから。測定で確かめていない部分は見積もりと書く)

1. **ライン上の節点は sweep ごとに対角を組み直している**: `cached = useDiagCache && loop > 0 && !onLine` なので、ライン上の節点は loop > 0 でも
   全面の FVS の対角・粘性の対角・キー 5 の熱伝導を計算する。ライン上の節点の対角を使うのは `storeLU` の sweep (loop 0) だけで、loop > 0 では捨てている。
   case/45 は全節点がライン上なので、point にある対角キャッシュの効果がまるごと消えている。block カーネルの増分 +1.98 ms (キー 5 で +2.86 ms) の主因と見ている。
2. **Thomas は 1 ライン 1 スレッドで 4719 スレッド**: 64 スレッド/ブロック × 74 ブロックで、80 SM の各々に 2 warp 程度しか載らない。
   各スレッドが 121 節点を直列に処理し、節点ごとに依存した大域メモリの読みと 5×5 の LU・代入を double で行う。
3. **double の演算**: A10G (GA102) の FP64 は FP32 の 1/64 の速さで、除算はさらに重い。solve は節点ごとに約 75 回の積和と 5 回の除算、
   factor は約 265 回の積和と 30 回の除算 (5 列の代入を含む)。
4. **メモリ量の下限 (見積もり)**: solve は節点あたり約 0.77 kB (Kprev・LU・W は double で各 200 B) を読み書きするので、帯域 600 GB/s でも 1 回 0.73 ms、5 sweep で 3.6 ms/step。
   factor は約 1.2 kB/節点で 1.15 ms。現状 (1.97 ms・5.41 ms) はこの下限の 2.7 倍・4.7 倍で、並列度の不足 (2.) が主と見ている。

### 4.3 改善の候補 (見込みはすべて見積もり、未測定)

| 案 | 内容 | 線形系・丸め | 見込み (ms/step) |
| --- | --- | --- | --- |
| A | ライン上の節点は `storeLU` の sweep 以外で対角を組まない (`DoDiag=false` の経路と粘性の対角の省略) | 不変 (使っていない量を計算しないだけ) | block カーネル 7.0 → 約 4.5 |
| B | Thomas をライン内で並列化: 1 ラインを 5〜8 スレッド (行ごと) か 1 warp に割り、5×5 の積・代入を行で並べる。solve は LU と行交換の代わりに逆行列 M̃⁻¹ を保存して行列ベクトル積だけにし、除算をなくす | 線形系は不変、丸めの順序が変わる | Thomas 15.3 → 約 5 (double のまま、メモリの下限に近づく) |
| C | Thomas の内部精度を `implicitSolvePrecision` に従わせる (0 = float で演算・LU/W/K を float で保存、1 = 従来の double) | 既定の ISP 0 で丸めが float になる (sweep の点解は ISP 0 で既に float) | B と合わせて約 2〜2.5 |
| D | K を単位ベクトルで 5 回呼んで列を抜く代わりに、閉形式 −A⁻ = n₂I + n_a1 r₁⊗l₁ + n_a5 r₅⊗l₅ を 1 回で 25 成分組む | 不変 (同じ式、丸めの順序が変わる) | loop 0 だけなので小 |

A〜D をすべて入れると、ラインの 1 step は 37 ms → 22〜24 ms (point 19 ms の約 1.2 倍) の見込み (未検証の目標)。

**codex plan 段 (2026-10-09、§6.1) の反映**: 順序は A → D → B (double・部分ピボット付き LU のまま並列化)。**Thomas の既定は double を維持**し、
float (C) と逆行列の保存は別の opt-in の実験に分ける (組立・係数保存・因子・前進代入の精度の組み合わせを表で定義してから)。
A の契約は「`onLine && storeLU == 0` の節点は対角を組まない、RHS の拘束 (壁・軸の rhs 0) は毎 sweep 残す」(`lineKFreeze 1` の後続 subiteration は sweep 0 でも `storeLU 0`)。

## 5. 実装ステップ

1. A (sweep カーネル、4d394a71 で実装済み)。2. D。3. B (`lineThomasFactor_d`・`lineThomasSolve_d` の並列化、double・LU のまま、旧版は退避スイッチで残す)。
4. (別の実験) C と逆行列の保存。各段で §6 の確認を行う。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | codex plan 段の採否 | 判断: 2026-10-09 GO-with-changes (C0/M5/m1) を全件採用 (§6.1)。順序 A → D → B、既定 double 維持 | F |
| 2 | A の確認 | 完了: 短期の一致は合格、性能 −1.24 ms/step (§6.0) | O |
| 3 | D | 閉形式の K。合格: 面積係数・TP の固有ベクトル・`rowDec`・既存の粘性の加算を含めて、旧の列抽出の K と同じ状態で成分の相対差 ≤ 1e-12 (FP64) | O |
| 4 | B (double・LU) | 不採用 (2026-10-09、§6.0): 3 つの書き方とも案 A より遅い。並列版は opt-in で残す | O |
| 5 | 検証の表 (§6) の他経路 | `lineImplicit 0`・部分被覆・可変長・`lineKFreeze 0/1`・周期ミラー・軸と壁の拘束行 (node のみ、cell は未検証と明記) | O |
| 6 | 粘性 Jacobian plan との統合の後の再検証 | [time_integration-line-viscous-jacobian](time_integration-line-viscous-jacobian.md) は同じ D/K の組立を変える。速度の A/B では D/K の仕様を固定し、統合後に §6 を回し直す。2026-10-09: 粘性 Jacobian は本線不採用 (値 2・3 は診断として保留、既定 0 の D/K は不変) なので、既定の経路の再検証は不要 | O |
| 7 | 逆行列の保存 (double、別の opt-in 実験) | 実装: `FORGE_LINE_INV=1` (§6.2 に事前登録)。判断: 2026-10-09 codex 諮問 — 次の候補は double の逆行列の保存、C (float) は後回し。部分ピボット付き double LU から逆行列を作り、代入だけを行列ベクトル積に替える (D/K の組立と保存の精度は固定)。候補の基準 (実装前に確定): 緩和前のスケーリングした後退誤差 η ≤ 1e-11、物理尺度で正規化した補正の差 ≤ 1e-8、新しい分解の失敗 0。性能は逆行列の生成の費用を含め、native・専有 GPU・profiler なしの反復計測と、同じ品質までの総時間。B の失敗は「その実装の不採用」で、律速の原因の確定ではない | F |
| 8 | C (Thomas を float、別の opt-in 実験) | #7 の後。精度の組み合わせの表と、線形残差・長期の収束性能で採否 | F |
| 9 | 本線の評価 (方向別 dt + 上限 + point 仕上げの総壁時計) | 粘性 Jacobian の探索を終えた後の本線 (codex 2026-10-09)。ライン陰解法の速度の改善はこの総壁時計で採否を決める。**粘性入りの腕を加える (2026-10-09、ユーザの指摘「粘性をライン向けのヤコビアンに入れたほうが良くね？」を受けた方針)**: ノズルで 2000 step 壊れなかった唯一の粘性入りの形 (値 3・マスク 5 = 運動量の薄層の結合と仕事を入れ、熱伝導の近傍 K を外す) を、値 0 と point だけの腕と同じ終了条件・総壁時計で比べる (viscous-jacobian plan §6.12 の条件: 仕上げ込みで速くなければ終える)。§6 に事前登録し codex plan 段に諮ってから回す。残差の停滞に面エンタルピーの精度が効くかの監査は [time_integration-implicit-thermal-jacobian](time_integration-implicit-thermal-jacobian.md) §5.1 #5 (未着手) | F |

## 6. 検証 (codex plan 段の反映後、2026-10-09)

- **短期の一致 (数値を変えない A・D)**: 旧 (現行、sha256 35e498b1…) と新で、run_0183 の res_100000 から run_0223 の設定 (ライン + 方向別 + キー 5 + 上限 50) を
  1 step と 20 step。旧・新それぞれ 2 本 (20 step) と 1 本 (1 step)。量: ρ・ρu・ρv・ρE・ρk・ρω・P・T の場の差の相対 RMS と最大絶対差。
  **合格: 旧 × 新の各量の相対 RMS・最大絶対差が、同じ設定の同じバイナリの再実行 (旧 × 旧・新 × 新) の最大の 3 倍以内** (1 step は再実行の差が丸めの桁なので、
  実装の誤りは桁違いの差として出る)。短期の一致は「収束解の一致」の証拠にしない。
- **線形解 (B)**: 凍結した実際の D・K・rhs から、緩和前のライン解の成分別にスケーリングした後退誤差 η = ‖b − Ax‖∞/(‖A‖∞‖x‖∞ + ‖b‖∞) と、
  独立した倍精度の解との差、行交換の回数、分解の失敗の件数を旧・新で比べる (合格値は B の実装前に登録する)。
- **長期の品質 (B 以降)**: run_0223 の設定で同じ step 数を回し、全残差の `check_convergence` と θ_r・Q_w・欠損の `check_quasisteady` を旧・新で比べる。
- **性能**: カーネルの内訳は nsys (§4.1)。**採否は native・専有 GPU・`FORGE_PROFILE` なし・profiler なしの壁時計**で、初期化と出力を除いた区間の
  ms/step を 3 回以上測って中央値とばらつきを出す。
- **他の経路**: §5.1 #5。

### 6.0 結果

- **A の短期の一致** (`case/45.isobutane_m6_d155/run_0254〜0259_abA_*`、`_band_ab/cold_pair/AB_abA.json`): **合格**。20 step の旧 × 新の相対 RMS・最大絶対差が 8 量とも再実行の最大の 3 倍以内
  (ρ の RMS 4.2e-9 vs 再実行 3.8e-9、ρω の最大 3.8e4 vs 2.1e4)。1 step の旧 × 新は ρ の RMS 1.2e-14 (V0 の directional・1 step の再実行 1.4〜1.6e-14 と同じ桁)。
- **A の性能** (`run_0268〜0274_timeA_*`、専有 GPU・`FORGE_PROFILE` なし・profiler なし・1000 step、forge_run.log の Time): 旧 36.29 / 36.22 / 36.22、新 34.99 / 34.96 / 34.98 ms/step
  → 中央値 36.22 → 34.98 (**−1.24 ms/step、−3.4 %**、ばらつき ±0.04)。point (新バイナリ) は 18.57 ms/step、ラインは point の 1.88 倍。見積もり (−2.5 ms) の半分。
- **順序の変更 (2026-10-09)**: D の見込みは loop 0 だけで 0.3 ms 程度なので、B (Thomas 15.3 ms) を先にする。B の合格の条件 (事前登録): 新旧の Thomas を同じ入力で両方解くデバッグの経路
  (`FORGE_LINE_COMPARE=1`) で、run_0223 の設定 20 step の全 sweep について factor の LU・W・ピボットと solve の補正 dq の新旧の最大絶対差 / 旧の最大 ≤ 1e-12 (double)、分解の失敗のライン数が一致。
  加えて A と同じ短期の一致 (旧 × 新が再実行の 3 倍以内) と、性能 (3 回の中央値)。
- **B の結果 (2026-10-09)**: 3 つの書き方を試し、**どれも案 A より遅かった**ので既定は 1 ライン 1 スレッドのまま (並列版は `FORGE_LINE_PAR=1` の opt-in で残す)。
  | 版 | commit | 内容 | 新旧の比較 (`FORGE_LINE_COMPARE=1`、20 step・120 回) | 壁時計 ms/step (中央値、3 回) | nsys の代入 / 分解 [ms/step] |
  | --- | --- | --- | --- | --- | --- |
  | 案 A のみ | 4d394a71 | 1 ライン 1 スレッド | — | 34.98 (取り直し 34.98) | 9.85 / 5.41 (run_0250、A 前) |
  | B1 | 64ed2cd6 | 8 レーン/ライン、各レーンが 5×5 の代入を重複 | 全回ビット一致 | 40.53 | 15.46 / 5.52 |
  | B2 | 874ae90a | 代入を行の分担、行のポインタを 1 回だけ選ぶ | 全回ビット一致 | 37.84・38.05 (3 本目はディスク満杯で失敗) | — |
  | B3 | 3764852d | 1 ライン 1 スレッドで lu5 の行の入れ替えを静的な添字に (スタック 240/40 B → 0 B) | 全回ビット一致 | 39.59 | 13.90 / 5.86 |
  | 戻し | d5538001 | lu5 を元に戻す (既定 = 案 A + 値 2 の経路) | — | 35.06・34.99 | — |
  読み (仮説、codex 諮問で「実行時間だけでは確定しない」): 旧版でも 4719 本のラインは全部同時に走っており、カーネルの時間は「1 本のラインの 121 節点の直列の連鎖」で決まると見ている。レーンを増やしても連鎖は短くならず、shfl の分だけ伸びた。
  レジスタ化は命令数が増えて代入を 9.85 → 13.90 ms に遅くした (原因の内訳は ncu が権限で使えず未確認)。連鎖そのものを短くするには、精度を変える C (float) か
  逆行列の保存 (1 節点の代入が 1 回の行列ベクトル積になる) が要る。どちらも数値が変わるので別の opt-in の実験 (§5.1 #7)。
  run: `case/45.isobutane_m6_d155/run_0275_cmpB`・`run_0276〜0281_abB_*` (短期の一致、1 step の ρ の RMS 1.2e-14)・`run_0283〜0285_timeB_*`・`run_0286_prof_lineB`・`run_0287_cmpB2`・`run_0289〜0291_timeB2_*`・
  `run_0292_cmpB3`・`run_0293〜0295_timeB3_*`・`run_0296_timeA_recheck`・`run_0297_prof_lineB3`・`run_0298・0299_timeD_*`。

### 6.2 事前登録: §5.1 #7 逆行列の保存 (`FORGE_LINE_INV=1`、2026-10-09、codex plan 段 [記録](../../notes/reviews/2026-10-09-time_integration-line-implicit-speed-plan-2.md) の反映後)

- **実装**: `lineThomasFactor_d<true>` が部分ピボット付き double LU (従来と同じ `lu5_factor`) から単位ベクトルの代入 5 回で M̃⁻¹ を作り、LU の代わりに保存する。
  W_k = M̃⁻¹Knext_k は従来どおり LU の代入 (不変)。`lineThomasSolve_d<true>` の前進は `lu5_solve` の代わりに 5×5 の行列ベクトル積。後退代入・D/K の組立・保存の精度は不変。
  既定 (変数なし) の経路は変えない。`FORGE_LINE_PAR`・`FORGE_LINE_MONO`・`FORGE_LINE_DEBUG_POINT`・`FORGE_LINE_NOOP` とは組み合わせない (最初の呼び出しで止める)。
  実効のモード (LU / INV / PAR / …) を factor・solve の 1 回目と 1000 回ごとにログへ書く。
- **バイナリ**: lineI = この節を書いた commit + typedef double (AWS、`~/forge-linevisc-fp64` を上書き)。設定は run_0223 と同じ (値 0・方向別・キー 5・上限 50・cfl 4・緩和 0.7・sweep 5・ISP 0)、
  初期場は run_0183 の res_100000。尺度は run_0183 のリミッタの基準値 (ρ_ref 0.8740、a_ref 359.77 m/s): ρ → ρ_ref、ρu → ρ_ref·a_ref = 314.4、ρE → ρ_ref·a_ref² = 1.131e5。
- **(1) 全ラインの照合と後退誤差** `run_0324_invcmp`: `FORGE_LINE_COMPARE=1` + `FORGE_LINE_INV=1` で 20 step (factor 20 回・solve 100 回)。従来の LU を別のバッファで解き、逆行列の解で進める。
  各 solve の直後に、保存した D・K・rhs と緩和前の解から、**全ライン**の η = ‖b̂ − Âx̂‖∞ / (‖Â‖∞‖x̂‖∞ + ‖b̂‖∞) (尺度で無次元化) を逆行列・従来の両方で出す (`[lineEta]` の行: 最大・ライン・評価本数・1e-11 超の本数・分解の失敗の本数)。
  合格 (すべて): 両方の腕の係数・因子・補正に非有限 0 件、全 solve で成分ごとの max|Δdq_c| / 尺度_c ≤ 1e-8 (dq は緩和後)、W の差 0、分解の失敗の不一致 0 (各腕の失敗の本数も記録)、
  逆行列の η ≤ 1e-11 (評価した全ライン・全 solve)。**従来の η も 1e-11 を超える solve があれば**、逆行列の精度不足とは断定せず、その solve は共通の系の問題として判別不能とする。
  成分ごとの max|Δdq_c| / max|dq_c| も記録する。
- **(2) 独立の照合 (CPU)** `run_0325_invdump_c1`・`run_0326_invdump_c20` (逆行列) と `run_0327_ludump_c1`・`run_0328_ludump_c20` (従来): factor の 1・20 回目の書き出し。
  `FORGE_LINE_DUMP_NODES` は節点番号で、節点 1572・4113・5321・7985・198560・264263 を含むライン (S0 の書き出しでライン 12・33・43・65・1640・2183) を選ぶ。
  `line_eta.py` (密行列の CPU 計算) の η が (1) の GPU の η と桁で一致することを確かめる (独立した実装の照合。合否は (1) で決める)。
- **(3) 性能** `run_0329〜0331_timeLU_*`・`run_0332〜0334_timeINV_*`: 1000 step、交互に 3 組。比較・書き出し・profiler の環境変数は無効。投入前に GPU 上の計算プロセス
  (`nvidia-smi --query-compute-apps`) が自分のもの以外に無いことを記録する。量は各 run の **101〜999 step の 1 step ごとの ms の平均** (暖機と最後の出力の step を除く、ログの各 step の行から)。
  判定: 組ごとの差 (INV − LU) が 3 組とも負で、その中央値の大きさが従来の 3 本の幅 (最大 − 最小) を超える → 速い。3 組とも正で同じ条件 → 遅い。それ以外 → 判別不能。
- **(4) 途中の到達時間 (予備の選別)** `run_0335_qualLU`・`run_0336_qualINV`: 10000 step・500 ごと (extraFields `res_ro`)。到達 = |2π Σ res_ro| ≤ 1.0 kg/s が step 0 を除く出力で
  **2 回続けて**成り立つ最初の区間で、その手前の出力との線形補間の step。到達までの累積壁時計はログの各 step の ms の和で直接出す (単価 × step の推定にしない)。
  片方・両方が届かなければ打ち切りとして記録し、採用の候補の判定は保留。非有限・非物理値なし、`check_convergence` の VERDICT を記録 (短い試験の NOT CONVERGED は許す)。
  **これは品質の判定ではない**: 全残差・θ_r・Q_w を含む品質と、point 仕上げまでの総壁時計は §5.1 #9 で判定する。
- **分岐**: (1) の合格条件を外す → 不採用 (従来の η も超えた solve は判別不能として別に書く)。(3) で「遅い」→ 不採用、「判別不能」→ 速さの根拠なしとして保留。
  (1) 合格・(3) で速い・(4) で両方到達して到達の差が ±10 % 以内かつ累積壁時計が短い → **#9 の比較に逆行列の腕を使う候補** (既定にするかは #9 の後)。
  (4) の到達の差が 10 % を超える → 判別不能 (1 本ずつでは軌道の違いが性能差を上回る)。
- 監査: 中間の res は判定と監査が終わるまで残す。消すときは README に書く。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan (#7 の事前登録) | 2026-10-09 | [2026-10-09-time_integration-line-implicit-speed-plan-2.md](../../notes/reviews/2026-10-09-time_integration-line-implicit-speed-plan-2.md) | GO-with-changes, C0/M5/m1 | 全件採用: M1 非有限の数え上げ、M2 全ラインの η を比較の経路で、M3 到達を |Σ|・2 回連続・step 0 除外・500 ごと、M4 (4) は予備の選別で品質は #9、M5 時計を 1 step ごとの ms (101〜999) と組ごとの差で、m6 スイッチの併用を止める (§6.2) |
| plan | 2026-10-09 | [2026-10-09-time_integration-line-implicit-speed-plan.md](../../notes/reviews/2026-10-09-time_integration-line-implicit-speed-plan.md) | GO-with-changes, C0/M5/m1 | 全件採用: M1 B は double・LU のまま並列化し逆行列化を分離、後退誤差 η で判定 (§6)。M2 既定 double 維持、float は opt-in (§4.3 末尾、§5.1 #7)。M3 短期の一致と長期の品質を別ゲート (§6)。M4 他経路の表 (§5.1 #5)、A の契約を onLine && !storeLU と明記。M5 採否は profiler なしの反復計測 (§6)。m6 影響範囲と粘性 Jacobian plan との統合順 (§5.1 #6、§7) |

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/timeIntegration_d.cu` (sweep・Thomas・wrapper)、`block_dplur_jacobian_d.cuh` (D)、`solver_density_cuda/mesh/mesh.hpp`・`mesh.cpp` (因子の型と確保、C のとき)、
  `solver_density_cuda/input/solverConfig.*` (C の opt-in のとき)。`lineImplicit 0` の経路は変えない。

## 変更ログ

- 2026-10-09: 起票 (計測と候補の整理、draft)。
- 2026-10-09: codex plan 段 (GO-with-changes) を全件採用し in_progress。A を実装 (4d394a71)。
- 2026-10-09: A の性能 −1.24 ms/step。B は 3 つの書き方とも遅く不採用 (既定は 1 ライン 1 スレッド、並列版は opt-in)。
