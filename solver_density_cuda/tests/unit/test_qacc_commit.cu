// test_qacc_commit.cu — 保存量 FP64 影アキュムレータの単体試験
//   (plans/active/time_integration-fp64-accumulator.md §5.1 S0 / §6.0 G0)
//   build: nvcc --expt-relaxed-constexpr -I. -o test_qacc_commit tests/unit/test_qacc_commit.cu
//
// **製品コードにフックを入れない**ため、commit / reconcile の device 関数を直接駆動する。
// plan に**事前登録した 3 セル**:
//   (a) 通常セル   ro=0.017755, dq=2.97e-10 (= 0.159 ULP), N=100 → Qacc-Q0 = N*dq (相対 1e-6)、
//                  ミラーは 16 ULP 動く。**現行構造 (FP32 累積) は 0 のまま**。
//   (b) 壁ピンセル 毎 step FP32 writer が Q を書き換える → reconcile が採用し Qacc == (double)Q になる。
//   (c) OFF 経路   影アキュムレータを通さない従来の加算 → 0.0 のまま (現行の症状の再現)。
// さらに (d) ビット同一性: 残余ゼロの 1 step は ON と OFF がビット一致する。
#include <cstdio>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "flowFormat.hpp"
#include "cuda_forge/qAccumulator_d.cuh"

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++g_fail; printf("  FAIL: " __VA_ARGS__); printf("\n"); } } while (0)

// (a)(c): N step 回して、ON (影アキュムレータ) と OFF (FP32 累積) の到達点を返す。
__global__ void k_accumulate(double* qacc, flow_float* q_on, flow_float* q_off,
                             flow_float dq, int nStep)
{
    if (blockDim.x * blockIdx.x + threadIdx.x != 0) return;
    for (int n = 0; n < nStep; ++n) {
        qaccCommitCell(qacc[0], q_on[0], dq, (flow_float)1.0);   // ON
        q_off[0] = q_off[0] + dq;                                 // OFF (現行構造と同じ)
    }
}

// (b): 毎 step、commit のあとに別の writer が Q を上書きし、reconcile で採用させる。
__global__ void k_pinned(double* qacc, flow_float* q, const flow_float* pinVal,
                         flow_float dq, int nStep, int* nAdopt)
{
    if (blockDim.x * blockIdx.x + threadIdx.x != 0) return;
    int c = 0;
    for (int n = 0; n < nStep; ++n) {
        qaccCommitCell(qacc[0], q[0], dq, (flow_float)1.0);
        q[0] = pinVal[0];                                         // 壁ピン相当の FP32 writer
        if (qaccReconcileCell(qacc[0], q[0])) ++c;
    }
    nAdopt[0] = c;
}

// (d): 残余ゼロの 1 step が ON/OFF でビット一致するか (多数サンプル)
__global__ void k_bitwise(const flow_float* a, const flow_float* b, int* nDiff, int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= n) return;
    double acc = (double)a[i];
    flow_float on = a[i];
    qaccCommitCell(acc, on, b[i], (flow_float)1.0);
    const flow_float off = a[i] + b[i];
    if (on != off) atomicAdd(nDiff, 1);
}

int main()
{
    const double RO0 = 0.0177550;
    const flow_float dq = (flow_float)2.97e-10;
    const int N = 100;
    // **float32 の値を基準に取る**。RO0 (double) との差 (-7.9e-11) を累積と取り違えないため。
    const double RO0f = (double)(flow_float)RO0;
    const double ulp = (double)std::nextafter((flow_float)RO0, (flow_float)1.0) - RO0f;
    printf("ro=%.7f  1 ULP=%.4e  dq=%.4e (= %.3f ULP)  N=%d\n",
           RO0, ulp, (double)dq, (double)dq / ulp, N);

    double *d_acc; flow_float *d_on, *d_off, *d_pin; int* d_n;
    cudaMalloc(&d_acc, sizeof(double)); cudaMalloc(&d_on, sizeof(flow_float));
    cudaMalloc(&d_off, sizeof(flow_float)); cudaMalloc(&d_pin, sizeof(flow_float));
    cudaMalloc(&d_n, sizeof(int));

    // --- (a) 通常セル / (c) OFF 経路 ---
    {
        double acc = RO0; flow_float q = (flow_float)RO0;
        cudaMemcpy(d_acc, &acc, sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_on, &q, sizeof(flow_float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_off, &q, sizeof(flow_float), cudaMemcpyHostToDevice);
        k_accumulate<<<1, 32>>>(d_acc, d_on, d_off, dq, N);
        cudaDeviceSynchronize();
        double accOut; flow_float onOut, offOut;
        cudaMemcpy(&accOut, d_acc, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&onOut, d_on, sizeof(flow_float), cudaMemcpyDeviceToHost);
        cudaMemcpy(&offOut, d_off, sizeof(flow_float), cudaMemcpyDeviceToHost);
        const double want = N * (double)dq;
        const double got = accOut - RO0;          // Qacc は double の RO0 から始めたので RO0 基準
        const double mirrorUlp = ((double)onOut - RO0f) / ulp;
        printf("\n(a) 通常セル: Qacc-Q0 = %.4e (期待 %.4e, 相対差 %.2e)\n", got, want, fabs(got / want - 1.0));
        printf("    ミラー (float32) は %.0f ULP 動いた\n", mirrorUlp);
        printf("(c) OFF 経路: Q-Q0 = %.4e  <- 現行構造の症状 (float32 の初期値基準)\n",
               (double)offOut - RO0f);
        CHECK(fabs(got / want - 1.0) < 1e-6, "(a) 影アキュムレータが N*dq を積めていない");
        CHECK(fabs(mirrorUlp - 16.0) <= 1.0, "(a) ミラーの移動が 16 ULP でない (%.1f)", mirrorUlp);
        CHECK((double)offOut - RO0f == 0.0, "(c) OFF 経路が 0 でない (症状が再現しない)");
    }

    // --- (b) 壁ピンセル ---
    {
        const flow_float pin = (flow_float)(RO0 * 1.001);   // 別の writer が毎 step 書く値
        double acc = RO0; flow_float q = (flow_float)RO0;
        cudaMemcpy(d_acc, &acc, sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_on, &q, sizeof(flow_float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_pin, &pin, sizeof(flow_float), cudaMemcpyHostToDevice);
        cudaMemset(d_n, 0, sizeof(int));
        k_pinned<<<1, 32>>>(d_acc, d_on, d_pin, dq, N, d_n);
        cudaDeviceSynchronize();
        double accOut; flow_float qOut; int nAdopt;
        cudaMemcpy(&accOut, d_acc, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&qOut, d_on, sizeof(flow_float), cudaMemcpyDeviceToHost);
        cudaMemcpy(&nAdopt, d_n, sizeof(int), cudaMemcpyDeviceToHost);
        // **期待するのは「ピンから離れられないこと」**。毎 step 採用されるわけではない:
        // 残余が 1/2 ULP を超えるまでは (float)Qacc == pin なので reconcile は発火しない。
        // 理論上の発火回数 ≈ N * (dq/ULP) / 0.5 = 100 * 0.159 / 0.5 ≈ 32 回。
        const double ulpPin = (double)std::nextafter(pin, (flow_float)1.0) - (double)pin;
        const double drift = fabs(accOut - (double)pin) / ulpPin;
        printf("\n(b) 壁ピンセル: Qacc = %.9f, Q(pin) = %.9f, ピンからの隔たり %.3f ULP, 採用 %d/%d step\n",
               accOut, (double)pin, drift, nAdopt, N);
        CHECK(drift < 1.0, "(b) Qacc がピンから 1 ULP 以上離れた (%.3f ULP) = 累積が逃げている", drift);
        CHECK(nAdopt > N / 8 && nAdopt < N, "(b) 採用回数 %d が想定域 (N/8, N) の外", nAdopt);
    }

    // --- (d) 残余ゼロの 1 step は ON/OFF でビット一致 ---
    {
        const int M = 200000;
        std::vector<flow_float> a(M), b(M);
        unsigned s = 12345u;
        auto rnd = [&s]() { s = s * 1664525u + 1013904223u; return (double)(s >> 8) / 16777216.0; };
        for (int i = 0; i < M; ++i) {
            a[i] = (flow_float)(1e-3 + rnd() * 1e3);
            b[i] = (flow_float)((rnd() - 0.5) * 2e-6);
        }
        flow_float *da, *db; cudaMalloc(&da, M * sizeof(flow_float)); cudaMalloc(&db, M * sizeof(flow_float));
        cudaMemcpy(da, a.data(), M * sizeof(flow_float), cudaMemcpyHostToDevice);
        cudaMemcpy(db, b.data(), M * sizeof(flow_float), cudaMemcpyHostToDevice);
        cudaMemset(d_n, 0, sizeof(int));
        k_bitwise<<<(M + 255) / 256, 256>>>(da, db, d_n, M);
        cudaDeviceSynchronize();
        int nDiff; cudaMemcpy(&nDiff, d_n, sizeof(int), cudaMemcpyDeviceToHost);
        printf("\n(d) 残余ゼロの 1 step: ON と OFF が違うサンプル %d / %d\n", nDiff, M);
        CHECK(nDiff == 0, "(d) 残余ゼロでも ON/OFF がビット一致しない (%d 件)", nDiff);
        cudaFree(da); cudaFree(db);
    }

    // --- (e) restart の代償: 残余を捨てる頻度で効果がどこまで消えるか ---
    //     (plan §5.1 S1b-⑤。codex result M5 の反例を試験に残す)
    //     残余は常に ½ ULP 以下だが、**捨てる頻度が上がると累積が丸ごと消える**。
    {
        const flow_float q0 = (flow_float)1.0;
        const flow_float dqs = (flow_float)ldexp(1.0, -26);   // 0.125 ULP @ Q=1
        const int N = 100;
        const double ulp1 = (double)std::nextafter(q0, (flow_float)2.0) - (double)q0;
        printf("\n(e) restart の代償 (Q=1, dq=%.3f ULP, N=%d step):\n", (double)dqs / ulp1, N);
        printf("    %-24s %14s\n", "残余を捨てる間隔", "ミラーの移動 [ULP]");
        struct R { const char* name; int every; };
        const R cases[] = {{"捨てない (連続)", 0}, {"100 step ごと", 100}, {"10 step ごと", 10},
                           {"1 step ごと (毎回)", 1}};
        double moved[4];
        for (int k = 0; k < 4; ++k) {
            double acc = (double)q0; flow_float q = q0;
            for (int n = 0; n < N; ++n) {
                acc += (double)dqs; q = (flow_float)acc;
                // restart 相当: FP32 の Q から正本を作り直す = 残余 (≤½ ULP) を捨てる
                if (cases[k].every > 0 && (n + 1) % cases[k].every == 0) acc = (double)q;
            }
            moved[k] = ((double)q - (double)q0) / ulp1;
            printf("    %-24s %14.1f\n", cases[k].name, moved[k]);
        }
        CHECK(moved[0] >= 11.0, "(e) 連続で 12 ULP 動かない (%.1f)", moved[0]);
        CHECK(moved[3] == 0.0, "(e) 毎 step restart なのに動いた (%.1f) — 反例が再現していない", moved[3]);
        CHECK(moved[1] > moved[3], "(e) 100 step ごとの方が毎回より動かない");
        printf("    -> **残余は常に ½ ULP 以下でも、捨てる頻度が上がると効果が丸ごと消える**。\n");
        printf("       「½ ULP しか失わないから実害なし」と一般化しないこと。\n");
    }

    printf("\n%s (失敗 %d)\n", g_fail ? "VERDICT: FAIL" : "VERDICT: PASS", g_fail);
    return g_fail ? 1 : 0;
}
