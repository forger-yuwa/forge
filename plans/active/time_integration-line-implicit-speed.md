# ライン陰解法 (`lineImplicit`) の 1 step あたりの費用を下げる

## メタ

- **area**: `time_integration`
- **status**: `draft`
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

A〜D をすべて入れると、ラインの 1 step は 37 ms → 22〜24 ms (point 19 ms の約 1.2 倍) の見込み。A は数値を変えないので先に入れられる。

## 5. 実装ステップ

1. A (sweep カーネル)。2. B (`lineThomasFactor_d`・`lineThomasSolve_d` の書き直し、旧版は退避スイッチで残す)。3. C (精度のテンプレート化)。4. D。
各段で §6 の計測と一致の確認を行う。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | 方針の決定 (ユーザ) と codex plan 段 | A〜D のどこまでやるか。cuda_forge の数値カーネルの変更なので plan 段のレビューと諮問を経てから実装する | F |
| 2 | A の実装と確認 | `timeIntegration_d.cu` の sweep。合格: 1 step・20 step の場の差が同じバイナリの再実行の差 (V0_repeat.json の幅) に入る、block カーネルの ms/step の減少を nsys で確認 | O |
| 3 | B・C・D | 合格は §6 | O |

## 6. 検証 (案)

- **一致**: 旧 (現行) と新で、run_0183 の res_100000 から 1 step と 20 step の場の差の相対 RMS を、同じバイナリの再実行の差
  (`_band_ab/cold_pair/V0_repeat.json`、directional・20 step の ρ で 0.85〜1.2e-7) と比べる。C (float) は丸めが変わるので別の基準
  (同じ状態でのライン解の相対差と、run_0223 の設定で 5000 step の欠損・θ_r の軌道) を事前に決める。
- **速さ**: §4.1 と同じ手順 (nsys、GPU 専有、500 step) でカーネル別の ms/step。
- **他のケース**: ライン陰解法を使う case/39 (DDES、dual-time) の FP32 ビルドでも 1 本確かめる。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |

## 7. 影響範囲

- `solver_density_cuda/cuda_forge/timeIntegration_d.cu`、`block_dplur_jacobian_d.cuh`。`lineImplicit 0` の経路は変えない。

## 変更ログ

- 2026-10-09: 起票 (計測と候補の整理、draft)。
