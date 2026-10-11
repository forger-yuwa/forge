// test_quantize_state.cu — 診断 FORGE_DIAG_QUANTIZE_STATE の 1 値の処理と自己点検の集計・書式の単体試験
//   (plans/active/architecture-float-state-double-geometry.md §5.1 #19・§6.36)
//   build (FP64 のビルド専用。作業ツリーの flowFormat.hpp は float の既定のままにし、typedef を double にした写しを
//   別のディレクトリ <fp64> に置いて先に探させる):
//     cd solver_density_cuda
//     nvcc --expt-relaxed-constexpr -I<fp64> -I. -o test_quantize_state tests/unit/test_quantize_state.cu
//   (<fp64>/flowFormat.hpp は flowFormat.hpp の float の 4 行をコメントにし、double の 4 行を生かしたもの。
//    FP64 のビルドの写しの solver_density_cuda で -I. だけでもよい)
//
// **製品の関数を直接使う** (cuda_forge/quantizeStateDiag_d.cuh):
//   - quantizeStateValue: カーネル quantizeStateToFloat32_d (quantizeStateDiag_d.cu) が 1 CV の 1 量に当てる処理そのもの。
//     下の試験のカーネルは製品のカーネルと同じ 1 行 (ic < n なら q[ic] = quantizeStateValue(q[ic]))。
//   - quantizeStateRepresentable / quantizeStateCountNotRepresentable / quantizeStateFormat / quantizeStateCheckLine /
//     quantizeStateInteriorLine / quantizeStateWhenOuter / quantizeStateWhenSpecies: 自己点検の判定・集計・行。
// 項目:
//   (1) 乱数の double 20 万個 (float の範囲内。非正規化数の範囲・0・負・丁度中間 (偶数丸め)・float で表せる値を含む) を
//       デバイスで量子化し、結果が host の (double)(float)x と生のビットで一致する。host 側の quantizeStateValue も一致。
//   (2) 2 回当てても変わらない (冪等)。量子化の後は全部 float32 で表せる。
//   (3) 量子化の前に表せない値の個数 (quantizeStateCountNotRepresentable) = 量子化で値が変わった個数。
//   (4) 自己点検の集計 (区間・NaN・入口のピンの除外) と行の書式が期待の文字列と一致する。
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <string>
#include <vector>
#include <cuda_runtime.h>
#include "flowFormat.hpp"
#include "cuda_forge/quantizeStateDiag_d.cuh"

static_assert(sizeof(flow_float) == 8, "test_quantize_state は FP64 のビルド専用 (flowFormat.hpp の typedef double の写しで build する)");

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++g_fail; printf("  FAIL: " __VA_ARGS__); printf("\n"); } } while (0)

#define CUDA_OK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); return 2; } } while (0)

// 製品のカーネル (quantizeStateDiag_d.cu の quantizeStateToFloat32_d) と同じ本文。
__global__ void quantizeKernel(geom_int n, flow_float* q)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < n) q[ic] = quantizeStateValue(q[ic]);
}

static uint64_t bits(double x) { uint64_t u; std::memcpy(&u, &x, sizeof(u)); return u; }

// host の参照: (double)(float)x (実行時に変換する。最近接偶数、非正規化数を保つ)。
static double hostRef(double x)
{
    volatile double vx = x;
    volatile float f = static_cast<float>(vx);
    return static_cast<double>(f);
}

// デバイスで q を量子化する (製品と同じ 1 本のカーネル)。
static int runDevice(std::vector<double>& q, int block)
{
    const geom_int n = (geom_int)q.size();
    double* d = nullptr;
    CUDA_OK(cudaMalloc(&d, q.size()*sizeof(double)));
    CUDA_OK(cudaMemcpy(d, q.data(), q.size()*sizeof(double), cudaMemcpyHostToDevice));
    const int grid = (n + block - 1) / block;
    quantizeKernel<<<grid, block>>>(n, d);
    CUDA_OK(cudaGetLastError());
    CUDA_OK(cudaMemcpy(q.data(), d, q.size()*sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_OK(cudaFree(d));
    return 0;
}

int main()
{
    // ---- 入力 (20 万個) ----
    const size_t N = 200000;
    std::vector<double> x;
    x.reserve(N);
    const float fden = std::numeric_limits<float>::denorm_min();   // 2^-149
    // 特別な値: ±0、float の最小の非正規化数・最小の正規化数・最大、その半分・中間 (偶数丸めの両向き)
    const double specials[] = {
        0.0, -0.0,
        (double)fden, -(double)fden, (double)FLT_MIN, -(double)FLT_MIN, (double)FLT_MAX, -(double)FLT_MAX,
        std::ldexp(1.0, -150),            // 最小の非正規化数の丁度半分 → 偶数丸めで 0
        -std::ldexp(1.0, -150),           // → −0
        std::ldexp(3.0, -151),            // 1.5·2^-150 → 2^-149 (切り上げ)
        std::ldexp(3.0, -150),            // 1.5·2^-149 (丁度中間) → 2^-148 (偶数)
        std::ldexp(5.0, -150),            // 2.5·2^-149 (丁度中間) → 2^-148 (偶数)
        std::ldexp(1.0, -160),            // 非正規化数より小さい → 0
        1.0 + std::ldexp(1.0, -24),       // 1 と次の float の丁度中間 → 1 (偶数)
        1.0 + std::ldexp(3.0, -24),       // (1+2^-23) と (1+2^-22) の丁度中間 → 1+2^-22 (偶数)
        -(1.0 + std::ldexp(1.0, -24)),
        1.0 + std::ldexp(1.0, -30),       // 1 に近い表せない値
        (double)FLT_MAX - std::ldexp(1.0, 102),   // FLT_MAX の少し下 (表せない) → FLT_MAX
        101325.0, 1.2345678901234567e-3, -2.5e-7, 300.0, 1.0e6,
    };
    for (double v : specials) x.push_back(v);
    std::mt19937_64 rng(20261011ULL);
    std::uniform_int_distribution<int> expDist(-160, 127);     // 2^-160 .. 2^127 (非正規化数の範囲と、それより小さい値を含む)
    std::uniform_int_distribution<uint64_t> manDist(0, (1ULL << 52) - 1);
    std::uniform_int_distribution<int> kind(0, 9);
    std::uniform_real_distribution<float> fUni(-1.0e6f, 1.0e6f);
    while (x.size() < N) {
        const int k = kind(rng);
        double v;
        if (k == 0) {
            v = (double)fUni(rng);                                   // float で表せる値 (変わらないこと)
        } else if (k == 1) {
            // float の丁度中間: f と次の float の中点 (double で丁度表せる)
            const float f = fUni(rng);
            const float g = std::nextafter(f, 2.0f*std::fabs(f) + 1.0f);
            v = 0.5*((double)f + (double)g);
        } else {
            const double m = 1.0 + (double)manDist(rng) * std::ldexp(1.0, -52);   // [1, 2)
            v = std::ldexp(m, expDist(rng));
            if (v > (double)FLT_MAX) v = (double)FLT_MAX;            // float の範囲内に留める
            if (k >= 6) v = -v;                                       // 負
        }
        x.push_back(v);
    }
    // 網羅の記録
    long long nZero = 0, nNeg = 0, nSub = 0, nBelow = 0, nRep = 0;
    for (double v : x) {
        const double a = std::fabs(v);
        nZero += (a == 0.0); nNeg += std::signbit(v);
        nSub += (a > 0.0 && a < (double)FLT_MIN && a >= (double)fden);
        nBelow += (a > 0.0 && a < (double)fden);
        nRep += quantizeStateRepresentable(v);
    }
    printf("入力 %zu 個: 0 が %lld、負 %lld、float の非正規化数の範囲 %lld、それより小さい %lld、float で表せる %lld\n",
           x.size(), nZero, nNeg, nSub, nBelow, nRep);
    CHECK(nZero >= 2 && nNeg > 50000 && nSub > 1000 && nBelow > 1000 && nRep > 10000, "入力の網羅が足りない");

    // ---- (1) デバイスの量子化 = host の (double)(float)x (ビット一致) ----
    std::vector<double> y = x;
    if (runDevice(y, 128) != 0) return 2;   // ブロックの大きさで割り切れない長さ (200000 = 1562·128 + 64)
    long long nMis = 0, nMisHost = 0;
    for (size_t i = 0; i < N; ++i) {
        const double r = hostRef(x[i]);
        if (bits(y[i]) != bits(r)) { if (nMis < 5) printf("  不一致 x=%.17g device=%.17g host=%.17g\n", x[i], y[i], r); ++nMis; }
        if (bits(quantizeStateValue(x[i])) != bits(r)) ++nMisHost;
    }
    printf("(1) デバイスと host の (double)(float)x のビット不一致 %lld、host の quantizeStateValue との不一致 %lld\n", nMis, nMisHost);
    CHECK(nMis == 0, "(1) デバイスの量子化が (double)(float)x とビット一致しない (%lld 個)", nMis);
    CHECK(nMisHost == 0, "(1) host の quantizeStateValue が (double)(float)x とビット一致しない (%lld 個)", nMisHost);
    // 特別な値の期待 (偶数丸め・非正規化数・符号)
    {
        const double want[] = {
            0.0, -0.0,
            (double)fden, -(double)fden, (double)FLT_MIN, -(double)FLT_MIN, (double)FLT_MAX, -(double)FLT_MAX,
            0.0, -0.0, (double)fden, std::ldexp(1.0, -148), std::ldexp(1.0, -148), 0.0,
            1.0, 1.0 + std::ldexp(1.0, -22), -1.0, 1.0, (double)FLT_MAX,
        };
        const size_t nw = sizeof(want)/sizeof(want[0]);
        for (size_t i = 0; i < nw; ++i)
            CHECK(bits(y[i]) == bits(want[i]), "(1) 特別な値 #%zu: x=%.17g → %.17g (期待 %.17g)", i, x[i], y[i], want[i]);
    }

    // ---- (2) 冪等・量子化の後は全部表せる ----
    std::vector<double> z = y;
    if (runDevice(z, 512) != 0) return 2;
    long long nIdem = 0, nNotRepAfter = 0;
    for (size_t i = 0; i < N; ++i) {
        if (bits(z[i]) != bits(y[i])) ++nIdem;
        if (!quantizeStateRepresentable(y[i])) ++nNotRepAfter;
    }
    printf("(2) 2 回目で変わった値 %lld、1 回目の後に float32 で表せない値 %lld\n", nIdem, nNotRepAfter);
    CHECK(nIdem == 0, "(2) 2 回当てると変わる (%lld 個)", nIdem);
    CHECK(nNotRepAfter == 0, "(2) 量子化の後に表せない値が残る (%lld 個)", nNotRepAfter);

    // ---- (3) 表せない値の個数 = 量子化で変わった個数 ----
    long long nChanged = 0;
    for (size_t i = 0; i < N; ++i) nChanged += (bits(y[i]) != bits(x[i]));
    const long long nNotRep = quantizeStateCountNotRepresentable(x, 0, N);
    printf("(3) 量子化の前に表せない値 %lld、量子化で変わった値 %lld\n", nNotRep, nChanged);
    CHECK(nNotRep == nChanged, "(3) 表せない値の個数 %lld と量子化で変わった個数 %lld が違う", nNotRep, nChanged);
    CHECK(quantizeStateCountNotRepresentable(y, 0, N) == 0, "(3) 量子化の後の集計が 0 でない");

    // ---- (4) 自己点検の集計 (区間・NaN・入口のピンの除外) と行の書式 ----
    {
        // 10 CV (内点 6・ghost 4)。roUx に内点 2・ghost 1 の表せない値と、内点に NaN 1 (表せないに数える)。
        std::vector<double> ro(10, 1.25), roUx(10, -3.0), roN(10, 0.5);
        roUx[1] = 1.0 + std::ldexp(1.0, -30);
        roUx[4] = std::nan("");
        roUx[5] = -7.0e-3;                        // 0.007 は float で表せない
        roUx[8] = std::ldexp(1.0, -160);          // 非正規化数より小さい (表せない)
        roN[9]  = std::ldexp(1.0, -149);          // 最小の非正規化数 (表せる)
        CHECK(!quantizeStateRepresentable(std::nan("")), "(4) NaN を表せないに数えない");
        CHECK(quantizeStateRepresentable(-0.0) && quantizeStateRepresentable(0.0), "(4) ±0 を表せないに数える");
        const long long cAll = quantizeStateCountNotRepresentable(roUx, 0, 10), cInt = quantizeStateCountNotRepresentable(roUx, 0, 6);
        CHECK(cAll == 4 && cInt == 3, "(4) 区間の集計: 全 %lld (期待 4)、内点 %lld (期待 3)", cAll, cInt);
        CHECK(quantizeStateCountNotRepresentable(roUx, 0, 99) == 4, "(4) 区間の上端が配列の長さを超えたときの集計");
        CHECK(quantizeStateCountNotRepresentable(roN, 0, 10) == 0, "(4) 非正規化数を表せないに数える");

        std::vector<QuantizeStateCountGroup> all = {
            {"Q",   {"ro", "roUx"}, {quantizeStateCountNotRepresentable(ro, 0, 10), cAll}},
            {"Q_N", {"roN"},        {quantizeStateCountNotRepresentable(roN, 0, 10)}},
        };
        std::vector<QuantizeStateCountGroup> in = {
            {"Q",   {"ro", "roUx"}, {quantizeStateCountNotRepresentable(ro, 0, 6), cInt}},
            {"Q_N", {"roN"},        {quantizeStateCountNotRepresentable(roN, 0, 6)}},
        };
        // 入口のピン (scalarDirichletPin、内点 [nCells] だけ): 内点 5 をピン → ピンを除く内点の roUx は 2 (1・4)
        std::vector<double> pin(6, 0.0);
        pin[5] = 1.0;
        long long nPin = 0;
        for (double p : pin) nPin += (p == 1.0);
        CHECK(quantizeStateCountNotRepresentable(roUx, 0, 6, &pin) == 2, "(4) ピンを除く集計 (期待 2)");
        CHECK(quantizeStateCountNotRepresentable(roUx, 0, 10, &pin) == 2, "(4) ピンの配列 (内点の長さ) より長い区間はピンの長さで切る (期待 2)");
        std::vector<QuantizeStateCountGroup> np = {
            {"Q",   {"ro", "roUx"}, {quantizeStateCountNotRepresentable(ro, 0, 6, &pin), quantizeStateCountNotRepresentable(roUx, 0, 6, &pin)}},
            {"Q_N", {"roN"},        {quantizeStateCountNotRepresentable(roN, 0, 6, &pin)}},
        };
        long long tAll = -1, tIn = -1, tNp = -1;
        const std::string body = quantizeStateFormat(all, &tAll), bodyIn = quantizeStateFormat(in, &tIn);
        const std::string bodyNp = quantizeStateFormat(np, &tNp);
        const std::string line = quantizeStateCheckLine(7, body, 10, quantizeStateWhenOuter());
        const std::string lineIn = quantizeStateInteriorLine(7, 6, bodyIn, nPin, bodyNp);
        const std::string lineInNoPinArray = quantizeStateInteriorLine(7, 6, bodyIn);   // ピンの配列が無いとき
        const std::string wantLine = "[quantizeState] step 7: 表せない値 ro 0 roUx 4 (Q) / roN 0 (Q_N) / 合計 4 "
                                     "(全 CV 10、updateVariablesOuter の直後)";
        const std::string wantInHead = "[quantizeState] step 7: 注意: 表せない値がある (0 を期待)。内点 (ic < nCells 6) に限ると "
                                       "ro 0 roUx 3 (Q) / roN 0 (Q_N) / 合計 3";
        const std::string wantIn = wantInHead + "。入口のピン (scalarDirichletPin 1、1 節点) を除く内点では "
                                   "ro 0 roUx 2 (Q) / roN 0 (Q_N) / 合計 2";
        printf("(4) %s\n    %s\n", line.c_str(), lineIn.c_str());
        CHECK(tAll == 4 && tIn == 3 && tNp == 2, "(4) 合計 %lld/%lld/%lld (期待 4/3/2)", tAll, tIn, tNp);
        CHECK(line == wantLine, "(4) 行の書式が違う:\n    得た '%s'\n    期待 '%s'", line.c_str(), wantLine.c_str());
        CHECK(lineIn == wantIn, "(4) 内点の行の書式が違う:\n    得た '%s'\n    期待 '%s'", lineIn.c_str(), wantIn.c_str());
        CHECK(lineInNoPinArray == wantInHead, "(4) ピンの配列が無いときの内点の行が違う:\n    得た '%s'\n    期待 '%s'",
              lineInNoPinArray.c_str(), wantInHead.c_str());

        // 期待どおり (全部 0) の行と、(b) の roY_N・roY_M の行
        std::vector<QuantizeStateCountGroup> zero = {
            {"roY_N", {"roY0N", "roY1N"}, {0, 0}},
            {"roY_M", {"roY0M", "roY1M"}, {0, 0}},
        };
        long long tz = -1;
        const std::string lineSp = quantizeStateCheckLine(1, quantizeStateFormat(zero, &tz), 576000, quantizeStateWhenSpecies(2));
        const std::string wantSp = "[quantizeState] step 1: 表せない値 roY0N 0 roY1N 0 (roY_N) / roY0M 0 roY1M 0 (roY_M) / 合計 0 "
                                   "(全 CV 576000、step 2 の化学種の更新の頭 (speciesUpdateOuter の直後))";
        printf("(4) %s\n", lineSp.c_str());
        CHECK(tz == 0 && lineSp == wantSp, "(4) (b) の行の書式が違う:\n    得た '%s'\n    期待 '%s'", lineSp.c_str(), wantSp.c_str());
    }

    printf("%s (失敗 %d)\n", g_fail == 0 ? "PASS" : "FAIL", g_fail);
    return g_fail == 0 ? 0 : 1;
}
