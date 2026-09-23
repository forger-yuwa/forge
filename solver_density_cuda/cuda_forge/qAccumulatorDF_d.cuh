// qAccumulatorDF_d.cuh — 保存量アキュムレータの **double-float (float 2 本)** 版
//   (plans/active/time_integration-fp64-accumulator.md §4.3 の代替案。FP64 影と同じ役目)
//
// **狙い**: FP64 影は `Qacc` に 8 B/変数を要る (5 保存量で 40 B/CV)。double-float なら
// `Q` (既にある float ミラー) を上位、残余 `lo` を下位に置けるので **追加は 4 B/変数 = 20 B/CV**。
// 実効精度は ~48 bit で、必要なのは「1 ULP 未満の残余を溜める」ことだけなので十分。
// 加算は two-sum (Knuth) 2 回 + 加算 1 回 = **加減算 13 回 + scale 乗算**
// (`twoSum` は 6 加減算。2026-09-24 に 7 → 13 へ再訂正、codex result-2 m2。
// 「3 flop」は `fastTwoSum` 1 回分の値で、commit 1 回の費用ではなかった)。
//
// **罠とその実測** (2026-09-24): two-sum は「加算の丸め誤差をちょうど拾う」ことに依存する。
//     s = a + b;  e = b - (s - a);
// コンパイラが融合・再結合すると `e` が 0 に潰れ、黙って素の float に劣化しうる。
// **ただし nvcc では実際には起きなかった** (`tests/unit/test_qacc_df.cu` で実測):
//   - 素の `+`/`-` 版 (`-DFORGE_DF_NAIVE`)、`--use_fast_math` 併用、`scale` を 1/0.5/0.375/0.123 と
//     変えた版、`-fmad=false` — **どれも intrinsic 版と同一の結果** (FP64 影との差 9.437e-16)。
//   - 理由は 2 つで、いずれも構造的:
//     (1) nvcc の `--use_fast_math` は除算・sqrt・超越関数・非正規化数に効くもので、
//         **FP 加算の再結合はしない** (CPU の `-ffast-math` とは違う)。
//     (2) `dq * scale` は `s` と `e` の**両方に現れる共通部分式**なので積が値として確定する。
//         `s` 側だけ fma に融合すると積を二度違う丸めで計算することになるため、コンパイラはやらない。
// したがって intrinsic は**保険**であって必須ではない。ただし将来の書き換えで (2) の性質が
// 崩れると**黙って劣化する**ので、`test_qacc_df.cu` の (e) (FP64 影との一致) を番人にする。
#ifndef Q_ACCUMULATOR_DF_D_CUH
#define Q_ACCUMULATOR_DF_D_CUH

#include "flowFormat.hpp"

#ifdef FORGE_DF_NAIVE
__device__ inline float dfAdd(float a, float b) { return a + b; }
__device__ inline float dfSub(float a, float b) { return a - b; }
#else
__device__ inline float dfAdd(float a, float b) { return __fadd_rn(a, b); }
__device__ inline float dfSub(float a, float b) { return __fsub_rn(a, b); }
#endif

// fast-two-sum: **|a| >= |b| のときだけ**丸め誤差 e を厳密に返す (3 flop)。
// ⚠ **前提が破れると誤差を丸ごと落とす** (2026-09-24, codex plan M1 の反例を自分で再現):
//     a = 2^-25, b = 1  ->  s = 1.0, e = 0.0   (2^-25 が消える。正しくは e = 2.98e-08)
// 「commit では |Q| >> |dq| なので常に成立」と書いていたが**誤り**: 運動量はゼロ近傍・符号反転・
// 初期過渡で |Q| < |dq| になりうるし、正値性ガードも運動量成分にこの条件を課していない。
// **したがって本体では使わない。** 残すのは「前提が成り立つ場所での比較用」としてのみ。
__device__ inline void fastTwoSum(float a, float b, float& s, float& e)
{
    s = dfAdd(a, b);
    e = dfSub(b, dfSub(s, a));
}

// two-sum (Knuth): **大小関係を問わず**丸め誤差 e を厳密に返す (6 flop)。
// fast-two-sum の 2 倍の演算だが前提が要らない。保存量 5 本すべてに使えるのはこちらだけ。
__device__ inline void twoSum(float a, float b, float& s, float& e)
{
    s = dfAdd(a, b);
    const float bb = dfSub(s, a);
    e = dfAdd(dfSub(a, dfSub(s, bb)), dfSub(b, bb));
}

// commit: 上位 q と下位 lo の対に dq を積む。q は FP64 影版と同じく「そのまま使える値」。
__device__ inline void qaccDFCommitCell(float& q, float& lo, float dq, float scale)
{
    float s, e;
    // FORGE_DF_FASTTWOSUM: **旧実装を忠実に再現**する (試験が本当に捕まえるかの確認用)。
    // ⚠ 以前は `lo = dfAdd(lo, e)` まで #else に入れてしまい、**旧実装とは別の粗い壊れ方**を
    // 捕まえていた (2026-09-24, codex result-2 m2)。差し替えるのは two-sum の種類だけにする。
#ifdef FORGE_DF_FASTTWOSUM
    fastTwoSum(q, dq * scale, s, e);
#else
    twoSum(q, dq * scale, s, e);     // **大小関係を仮定しない** (運動量のゼロ近傍・符号反転に備える)
#endif
    lo = dfAdd(lo, e);               // 取りこぼした分を下位へ
    // 正規化: 下位が 1 ULP を超えたら上位へ畳む
    float s2, e2;
#ifdef FORGE_DF_FASTTWOSUM
    fastTwoSum(s, lo, s2, e2);
#else
    twoSum(s, lo, s2, e2);
#endif
    q  = s2;
    lo = e2;
}

// reconcile: FP32 側が別の writer に書き換えられていたらそれを採用し、下位を捨てる。
// **FP64 影より単純**: 正本が q 自身なので「q が変わったか」を別に持つ必要がある。
// ここでは呼び出し側が変更前の q を渡す。
__device__ inline bool qaccDFReconcileCell(float qNow, float qExpected, float& lo)
{
    if (qNow != qExpected) { lo = 0.0f; return true; }
    return false;
}

#endif  // Q_ACCUMULATOR_DF_D_CUH
