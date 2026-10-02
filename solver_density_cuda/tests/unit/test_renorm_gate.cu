// test_renorm_gate.cu — 再正規化の受入ゲートの計測 (cuda_forge/renormGate_d.cuh) の集計間隔不変性 (plan §5.1 #1b-pre;
//   codex diagnose notes/reviews/2026-10-02-twophase-1b-final-diagnose.md)。
//   同じ 100 更新の補正列 (丸め規模の合成入力) を、A = 毎更新でログ / B = 10 更新ごとにログ の 2 通りで集計し、
//   末尾窓 (ceil(0.1N) = 10 更新) の判定 (κ = 2 n_s ε₃₂, n_s = 2) が一致することを確かめる。比較として却下された計り方
//   (区間累積の Σ|Δq|V を現在総量で割る) が集計間隔で変わることも出す。境界例: 窓内 1 更新だけ κ 超過・総量 0 の成分・NaN。
// build: nvcc -O2 -arch=sm_86 -std=c++17 -I. -o test_renorm_gate tests/unit/test_renorm_gate.cu
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>
#include "cuda_forge/renormGate_d.cuh"

static int g_fail = 0;
static void verdict(bool ok, const char* m) { printf("[%s] %s\n", ok ? "PASS" : "FAIL", m); if (!ok) ++g_fail; }
#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { printf("CUDA %s\n", cudaGetErrorString(e_)); return 2; } } while (0)

// 再正規化の写し: ΣρY = ρ への係数 f を化学種 (2 本) と液・Q に掛け、float に戻す。直前・直後の格納値で rng_accumulate (本番と同じ関数)。
__global__ void k_renorm(int nC, float* rY0, float* rYw, float* rg, float* q2, float* q1, float* q0, const float* ro, double* acc, int liquidToo)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= nC) return;
    const double sum = (double)rY0[i] + (double)rYw[i];
    const double f = (double)ro[i]/sum;
    double qm[RNG_NC] = {rYw[i], rg[i], q2[i], q1[i], q0[i]};
    rY0[i] = (float)((double)rY0[i]*f); rYw[i] = (float)((double)rYw[i]*f);
    if (liquidToo) { rg[i] = (float)((double)rg[i]*f); q2[i] = (float)((double)q2[i]*f); q1[i] = (float)((double)q1[i]*f); q0[i] = (float)((double)q0[i]*f); }
    double qp[RNG_NC] = {rYw[i], rg[i], q2[i], q1[i], q0[i]};
    rng_accumulate(acc, qm, qp, 1.0, fabs(f - 1.0));
}

struct Seq { std::vector<double> hist; std::vector<double> cumNum, totNow; };   // 履歴 (100×6) と却下された計り方用の Σ|Δq|V・現在総量

// 100 更新の合成入力: 各更新の前に化学種を丸め規模 (±scale ULP) 乱し、再正規化する。spike: その更新だけ ρY_w を相対 dev ずらす。
static int run_sequence(Seq& S, int N, double ulps, int spikeAt, double spikeDev, bool liquidZero, bool nanAt)
{
    const int nC = 4096; const float rho = 0.36f;
    std::vector<float> rY0(nC), rYw(nC), rg(nC), q2(nC), q1(nC), q0(nC), ro(nC, rho);
    for (int i = 0; i < nC; ++i) { rYw[i] = rho*0.01095f; rY0[i] = rho - rYw[i]; rg[i] = liquidZero ? 0.0f : rho*1e-3f; q2[i] = liquidZero ? 0.0f : 2e-2f; q1[i] = liquidZero ? 0.0f : 7e7f; q0[i] = liquidZero ? 0.0f : 3e15f; }
    float *d[7]; double *acc, *hist;
    for (int k = 0; k < 7; ++k) CK(cudaMalloc(&d[k], nC*sizeof(float)));
    CK(cudaMalloc(&acc, RNG_ACC*sizeof(double))); CK(cudaMalloc(&hist, (size_t)N*RNG_ENTRY*sizeof(double)));
    std::mt19937 rng(7); std::uniform_int_distribution<int> U(-(int)ulps, (int)ulps);
    S.cumNum.assign((size_t)N*RNG_NC, 0.0); S.totNow.assign((size_t)N*RNG_NC, 0.0);
    for (int n = 0; n < N; ++n) {
        for (int i = 0; i < nC; ++i) {   // 丸め規模の乱れ (ulps 単位)
            for (int t = U(rng); t != 0; t += (t > 0 ? -1 : 1)) rY0[i] = std::nextafter(rY0[i], t > 0 ? 1.0f : 0.0f);
            for (int t = U(rng); t != 0; t += (t > 0 ? -1 : 1)) rYw[i] = std::nextafter(rYw[i], t > 0 ? 1.0f : 0.0f);
        }
        if (n == spikeAt) for (int i = 0; i < nC; ++i) rYw[i] = (float)(rYw[i]*(1.0 + spikeDev));
        if (nanAt && n == N - 1) rYw[0] = NAN;
        std::vector<float> wpre = rYw, gpre = rg;
        CK(cudaMemcpy(d[0], rY0.data(), nC*4, cudaMemcpyHostToDevice)); CK(cudaMemcpy(d[1], rYw.data(), nC*4, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d[2], rg.data(), nC*4, cudaMemcpyHostToDevice)); CK(cudaMemcpy(d[3], q2.data(), nC*4, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d[4], q1.data(), nC*4, cudaMemcpyHostToDevice)); CK(cudaMemcpy(d[5], q0.data(), nC*4, cudaMemcpyHostToDevice));
        CK(cudaMemcpy(d[6], ro.data(), nC*4, cudaMemcpyHostToDevice));
        CK(cudaMemset(acc, 0, RNG_ACC*sizeof(double)));
        k_renorm<<<(nC + 255)/256, 256>>>(nC, d[0], d[1], d[2], d[3], d[4], d[5], d[6], acc, liquidZero ? 0 : 1);
        rng_finalize_d<<<1, 1>>>(acc, hist + (size_t)n*RNG_ENTRY);
        CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
        CK(cudaMemcpy(rY0.data(), d[0], nC*4, cudaMemcpyDeviceToHost)); CK(cudaMemcpy(rYw.data(), d[1], nC*4, cudaMemcpyDeviceToHost));
        CK(cudaMemcpy(rg.data(), d[2], nC*4, cudaMemcpyDeviceToHost)); CK(cudaMemcpy(q2.data(), d[3], nC*4, cudaMemcpyDeviceToHost));
        CK(cudaMemcpy(q1.data(), d[4], nC*4, cudaMemcpyDeviceToHost)); CK(cudaMemcpy(q0.data(), d[5], nC*4, cudaMemcpyDeviceToHost));
        // 却下された計り方の材料 (ρY_w と ρg だけ): 更新の Σ|Δq|V と更新後の総量
        double dw = 0, dg = 0, tw = 0, tg = 0;
        for (int i = 0; i < nC; ++i) { dw += fabs((double)rYw[i] - wpre[i]); dg += fabs((double)rg[i] - gpre[i]); tw += rYw[i]; tg += rg[i]; }
        S.cumNum[(size_t)n*RNG_NC + 0] = dw; S.cumNum[(size_t)n*RNG_NC + 1] = dg; S.totNow[(size_t)n*RNG_NC + 0] = tw; S.totNow[(size_t)n*RNG_NC + 1] = tg;
    }
    S.hist.resize((size_t)N*RNG_ENTRY);
    CK(cudaMemcpy(S.hist.data(), hist, S.hist.size()*sizeof(double), cudaMemcpyDeviceToHost));
    for (int k = 0; k < 7; ++k) cudaFree(d[k]);
    cudaFree(acc); cudaFree(hist);
    return 0;
}

// ログ間隔 L で集計: 区間ログ (各区間の max) と終了時の末尾窓。返り値: 末尾窓の判定と、区間ログから末尾窓を覆う区間の max。
struct Agg { RngWindow fin; bool pass; double intervalCoverMax; double rejected; };
static Agg aggregate(const Seq& S, int N, int L, double kappa)
{
    Agg a{};
    const size_t b = rng_final_window_begin((size_t)N);
    double cover = 0.0, rej = 0.0;
    for (int s = 0; s < N; s += L) {
        const int e = std::min(N, s + L);
        const RngWindow w = rng_window_max(S.hist, (size_t)s, (size_t)e);
        if ((size_t)e > b) for (int k = 0; k < RNG_NC; ++k) cover = std::max(cover, w.C[k]);
        // 却下された計り方: 区間の Σ|Δρg|V の和を区間末の現在総量で割る
        double num = 0; for (int n = s; n < e; ++n) num += S.cumNum[(size_t)n*RNG_NC + 1];
        const double tot = S.totNow[(size_t)(e - 1)*RNG_NC + 1];
        if ((size_t)e == (size_t)N) rej = (tot > 0) ? num/tot : 0.0;
    }
    a.fin = rng_window_max(S.hist, b, (size_t)N);
    a.pass = rng_judge(a.fin, kappa);
    a.intervalCoverMax = cover; a.rejected = rej;
    return a;
}

static void report(const char* tag, const Seq& S, int N, double kappa, bool expectPass)
{
    const Agg A = aggregate(S, N, 1, kappa), B = aggregate(S, N, 10, kappa);
    auto eqv = [](double x, double y) { return x == y || (x != x && y != y); };   // NaN 同士も一致とみなす
    bool same = (A.pass == B.pass) && eqv(A.fin.F, B.fin.F) && (A.fin.nonfinite == B.fin.nonfinite);
    for (int k = 0; k < RNG_NC; ++k) same = same && eqv(A.fin.C[k], B.fin.C[k]);
    char m[640];
    snprintf(m, sizeof m, "%s: 末尾窓 [%zu,%d) A(毎更新)/B(10 更新ごと) とも %s/%s、max C: w %.3e g %.3e Q2 %.3e Q1 %.3e Q0 %.3e、max|f-1| %.3e、非有限 %ld、"
             "κ %.7e、区間ログで窓を覆う max C A %.3e / B %.3e (参考: 却下された区間累積/現在総量 [g] A %.3e / B %.3e)",
             tag, A.fin.begin, N, A.pass ? "PASS" : "FAIL", B.pass ? "PASS" : "FAIL", A.fin.C[0], A.fin.C[1], A.fin.C[2], A.fin.C[3], A.fin.C[4],
             A.fin.F, A.fin.nonfinite, kappa, A.intervalCoverMax, B.intervalCoverMax, A.rejected, B.rejected);
    verdict(same && A.pass == expectPass, m);
}

int main()
{
    const int N = 100;
    const double kappa = rng_kappa(2);
    printf("=== 再正規化ゲートの集計間隔不変性 (N = %d, 末尾窓 %zu 更新, κ = 2·2·ε₃₂ = %.7e) ===\n", N, (size_t)N - rng_final_window_begin(N), kappa);
    Seq s1; if (run_sequence(s1, N, 1, -1, 0.0, false, false)) return 2;
    report("(a) 丸め規模 ±1 ULP (期待 PASS)", s1, N, kappa, true);
    Seq s2; if (run_sequence(s2, N, 1, 95, 1e-3, false, false)) return 2;
    report("(b) 窓内の更新 95 だけ ρY_w を 1e-3 ずらす (f−1 ≈ 1.1e-5; 期待 FAIL)", s2, N, kappa, false);
    Seq s3; if (run_sequence(s3, N, 1, 50, 1e-3, false, false)) return 2;
    report("(c) 窓の外 (更新 50) だけ 1e-3 (期待 PASS)", s3, N, kappa, true);
    Seq s4; if (run_sequence(s4, N, 1, -1, 0.0, true, false)) return 2;
    report("(d) 液・Q の総量 0 で補正 0 (期待 PASS: C = 0)", s4, N, kappa, true);
    Seq s5; if (run_sequence(s5, N, 1, -1, 0.0, false, true)) return 2;
    report("(e) 最終更新に NaN (期待 FAIL)", s5, N, kappa, false);
    // 総量 0 で補正 > 0 (分母 0) は inf → FAIL (finalize の規則を直接)
    {
        double hacc[RNG_ACC] = {0}; hacc[1] = 1e-30; hacc[RNG_NC + 1] = 0.0;   // ρg: 分子 > 0, 分母 0
        double *dacc, *dent; cudaMalloc(&dacc, sizeof hacc); cudaMalloc(&dent, RNG_ENTRY*sizeof(double));
        cudaMemcpy(dacc, hacc, sizeof hacc, cudaMemcpyHostToDevice);
        rng_finalize_d<<<1, 1>>>(dacc, dent); cudaDeviceSynchronize();
        std::vector<double> h(RNG_ENTRY); cudaMemcpy(h.data(), dent, RNG_ENTRY*sizeof(double), cudaMemcpyDeviceToHost);
        const RngWindow w = rng_window_max(h, 0, 1);
        char m[200]; snprintf(m, sizeof m, "(f) 分母 0・分子 1e-30 の成分: C_g = %g、判定 %s (期待 FAIL)", h[1], rng_judge(w, kappa) ? "PASS" : "FAIL");
        verdict(!rng_judge(w, kappa) && std::isinf(h[1]), m);
        cudaFree(dacc); cudaFree(dent);
    }
    printf("\n%s\n", g_fail ? "FAIL あり" : "ALL PASS");
    return g_fail ? 1 : 0;
}
