// test_qacc_df.cu — double-float (float 2 本) アキュムレータを FP64 影と突き合わせる
//   (plans/active/time_integration-fp64-accumulator.md §4.3 の代替案検証)
//   build: nvcc --expt-relaxed-constexpr -I. -o test_qacc_df tests/unit/test_qacc_df.cu
//   劣化確認: 上に -DFORGE_DF_NAIVE を付けると intrinsic を外した版になる
//   縮約確認: さらに --use_fast_math を付ける
//
// S0 (FP64 影) と同じ 3 セルを見る:
//   (a) 通常セル ro=0.017755, dq=2.97e-10 (0.159 ULP), N=100
//   (b) 壁ピンセル  (c) OFF 経路  (d) 残余ゼロのビット同一性
// さらに (e) FP64 影との一致、(f) 追加メモリの比較。
#include <cstdio>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "flowFormat.hpp"
#include "cuda_forge/qAccumulator_d.cuh"
#include "cuda_forge/qAccumulatorDF_d.cuh"

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++g_fail; printf("  FAIL: " __VA_ARGS__); printf("\n"); } } while (0)

// (a)(c)(e): N step 回して double-float / FP64 影 / OFF の到達点を返す
__global__ void k_acc(float* q_df, float* lo_df, double* qacc, float* q_fp64mirror,
                      float* q_off, float dq, int nStep)
{
    if (blockDim.x * blockIdx.x + threadIdx.x != 0) return;
    for (int n = 0; n < nStep; ++n) {
        qaccDFCommitCell(q_df[0], lo_df[0], dq, 1.0f);            // double-float
        qaccCommitCell(qacc[0], q_fp64mirror[0], dq, 1.0f);        // FP64 影
        q_off[0] = q_off[0] + dq;                                  // OFF (現行構造)
    }
}

// (b): 毎 step 別の writer が上位を書き換える (壁ピン相当)
__global__ void k_pinned(float* q, float* lo, const float* pinVal, float dq, int nStep, int* nAdopt)
{
    if (blockDim.x * blockIdx.x + threadIdx.x != 0) return;
    int c = 0;
    for (int n = 0; n < nStep; ++n) {
        qaccDFCommitCell(q[0], lo[0], dq, 1.0f);
        const float before = q[0];
        q[0] = pinVal[0];                                          // 別の writer
        if (qaccDFReconcileCell(q[0], before, lo[0])) ++c;
    }
    nAdopt[0] = c;
}

// (d): 残余ゼロの 1 step は OFF とビット一致するか
__global__ void k_bitwise(const float* a, const float* b, int* nDiff, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= n) return;
    float q = a[i], lo = 0.0f;
    qaccDFCommitCell(q, lo, b[i], 1.0f);
    const float off = a[i] + b[i];
    if (q != off) atomicAdd(nDiff, 1);
}

int main()
{
    const double RO0 = 0.0177550;
    const float dq = (float)2.97e-10;
    const int N = 100;
    const double RO0f = (double)(float)RO0;
    const double ulp = (double)std::nextafter((float)RO0, 1.0f) - RO0f;
#ifdef FORGE_DF_NAIVE
    printf("*** FORGE_DF_NAIVE: intrinsic を外した版 (縮約・再結合の影響を見る) ***\n");
#else
    printf("intrinsic 版 (__fadd_rn / __fsub_rn)\n");
#endif
    printf("ro=%.7f  1 ULP=%.4e  dq=%.4e (= %.3f ULP)  N=%d\n\n",
           RO0, ulp, (double)dq, (double)dq / ulp, N);

    float *d_qdf, *d_lodf, *d_mir, *d_off, *d_pin; double* d_acc; int* d_n;
    cudaMalloc(&d_qdf, 4); cudaMalloc(&d_lodf, 4); cudaMalloc(&d_mir, 4);
    cudaMalloc(&d_off, 4); cudaMalloc(&d_pin, 4); cudaMalloc(&d_acc, 8); cudaMalloc(&d_n, 4);

    {   // --- (a)(c)(e) ---
        float q = (float)RO0, lo = 0.0f; double acc = RO0;
        cudaMemcpy(d_qdf, &q, 4, cudaMemcpyHostToDevice);
        cudaMemcpy(d_lodf, &lo, 4, cudaMemcpyHostToDevice);
        cudaMemcpy(d_mir, &q, 4, cudaMemcpyHostToDevice);
        cudaMemcpy(d_off, &q, 4, cudaMemcpyHostToDevice);
        cudaMemcpy(d_acc, &acc, 8, cudaMemcpyHostToDevice);
        k_acc<<<1, 32>>>(d_qdf, d_lodf, d_acc, d_mir, d_off, dq, N);
        cudaDeviceSynchronize();
        float qOut, loOut, mirOut, offOut; double accOut;
        cudaMemcpy(&qOut, d_qdf, 4, cudaMemcpyDeviceToHost);
        cudaMemcpy(&loOut, d_lodf, 4, cudaMemcpyDeviceToHost);
        cudaMemcpy(&mirOut, d_mir, 4, cudaMemcpyDeviceToHost);
        cudaMemcpy(&offOut, d_off, 4, cudaMemcpyDeviceToHost);
        cudaMemcpy(&accOut, d_acc, 8, cudaMemcpyDeviceToHost);

        const double want = N * (double)dq;
        const double gotDF = ((double)qOut + (double)loOut) - RO0f;   // 対の値 = hi + lo
        const double gotFP64 = accOut - RO0;
        printf("(a) 通常セル N=%d step:\n", N);
        printf("    double-float  (hi+lo)-Q0 = %.6e   期待 %.6e   相対差 %.2e\n", gotDF, want, fabs(gotDF/want - 1.0));
        printf("    FP64 影        Qacc-Q0   = %.6e   相対差 %.2e\n", gotFP64, fabs(gotFP64/want - 1.0));
        printf("    ミラーの移動: double-float %.0f ULP / FP64 影 %.0f ULP\n",
               ((double)qOut - RO0f)/ulp, ((double)mirOut - RO0f)/ulp);
        printf("(c) OFF 経路 Q-Q0 = %.4e  <- 現行構造の症状\n", (double)offOut - RO0f);
        printf("(e) double-float と FP64 影の差 = %.3e (相対 %.2e)\n\n",
               fabs(gotDF - gotFP64), fabs(gotDF - gotFP64)/want);
        CHECK(fabs(gotDF/want - 1.0) < 1e-6, "(a) double-float が N*dq を積めていない (相対差 %.2e)", fabs(gotDF/want-1.0));
        CHECK(fabs((double)qOut - (double)mirOut) == 0.0, "(a) 上位が FP64 影のミラーと一致しない");
        CHECK((double)offOut - RO0f == 0.0, "(c) OFF 経路が 0 でない");
        CHECK(fabs(gotDF - gotFP64)/want < 1e-6, "(e) double-float が FP64 影とずれた (相対 %.2e)", fabs(gotDF-gotFP64)/want);
    }

    {   // --- (b) 壁ピン ---
        const float pin = (float)(RO0 * 1.001);
        float q = (float)RO0, lo = 0.0f;
        cudaMemcpy(d_qdf, &q, 4, cudaMemcpyHostToDevice);
        cudaMemcpy(d_lodf, &lo, 4, cudaMemcpyHostToDevice);
        cudaMemcpy(d_pin, &pin, 4, cudaMemcpyHostToDevice);
        cudaMemset(d_n, 0, 4);
        k_pinned<<<1, 32>>>(d_qdf, d_lodf, d_pin, dq, N, d_n);
        cudaDeviceSynchronize();
        float qOut, loOut; int nAdopt;
        cudaMemcpy(&qOut, d_qdf, 4, cudaMemcpyDeviceToHost);
        cudaMemcpy(&loOut, d_lodf, 4, cudaMemcpyDeviceToHost);
        cudaMemcpy(&nAdopt, d_n, 4, cudaMemcpyDeviceToHost);
        const double ulpPin = (double)std::nextafter(pin, 1.0f) - (double)pin;
        const double drift = fabs(((double)qOut + (double)loOut) - (double)pin) / ulpPin;
        printf("(b) 壁ピン: hi+lo はピンから %.3f ULP、採用 %d/%d step\n\n", drift, nAdopt, N);
        CHECK(drift < 1.0, "(b) ピンから 1 ULP 以上離れた (%.3f)", drift);
    }

    {   // --- (d) 残余ゼロのビット同一性 ---
        const int M = 200000;
        std::vector<float> a(M), b(M);
        unsigned s = 12345u;
        auto rnd = [&s]() { s = s*1664525u + 1013904223u; return (double)(s >> 8) / 16777216.0; };
        for (int i = 0; i < M; ++i) { a[i] = (float)(1e-3 + rnd()*1e3); b[i] = (float)((rnd()-0.5)*2e-6); }
        float *da, *db; cudaMalloc(&da, M*4); cudaMalloc(&db, M*4);
        cudaMemcpy(da, a.data(), M*4, cudaMemcpyHostToDevice);
        cudaMemcpy(db, b.data(), M*4, cudaMemcpyHostToDevice);
        cudaMemset(d_n, 0, 4);
        k_bitwise<<<(M+255)/256, 256>>>(da, db, d_n, M);
        cudaDeviceSynchronize();
        int nDiff; cudaMemcpy(&nDiff, d_n, 4, cudaMemcpyDeviceToHost);
        printf("(d) 残余ゼロの 1 step: 上位が OFF と違うサンプル %d / %d\n\n", nDiff, M);
        CHECK(nDiff == 0, "(d) 上位が OFF とビット一致しない (%d 件)", nDiff);
    }

    printf("(f) 追加メモリ (5 保存量/CV): double-float %zu B  vs  FP64 影 %zu B\n",
           5*sizeof(float), 5*sizeof(double));
    printf("\n%s (失敗 %d)\n", g_fail ? "VERDICT: FAIL" : "VERDICT: PASS", g_fail);
    return g_fail ? 1 : 0;
}
