// test_farfield_flux.cu — 遠方境界 `farfield` の外側状態と HLLC の単体試験 (plans/active/boundary-node-farfield-characteristic.md §6 V0u・V0k)。
//   build: nvcc -std=c++17 --expt-relaxed-constexpr -I. -Icuda_forge -o test_farfield_flux tests/unit/test_farfield_flux.cu
//   V0u (ホスト、__host__ __device__ の ff_outer_state / ff_hllc を直接呼ぶ):
//     (iv)   一様 (U_i = U_∞) で数値流束 = 物理流束 (相対 1e-6)、CPG と TP、亜音速/超音速の流入・流出・接線流
//     (ii)   接触波 (P・U_n 同じ、T 600/220 K、Y 0.13/0): 流出 = 内部の物理流束、流入 = 自由流の物理流束 (相対 1e-6)
//     (i)    連続性: U_n,i・ρ_i・P_i・U_n,∞ を法線 Mach 0・±0.9・±1 付近で 1e-6 刻みに掃引、隣接点の流束差 ≤ 1e-4 × 流束規模
//     (iii)  流向の組合せ: 混合帯の外では外側状態の組成・k・ω = 自由流、常に内部値と自由流値の間 (有界)
//     (iii'') 真空と float32: 独立掃引で非有限 0・外側状態の負値 0、真空置換と HLL 退避を別々に数える
//     (iii''') 混合帯の逆流: Y_R = 0.065・k_R = 0.075 (CPG/TP、SST k 込み)
//   V0k (GPU): 実カーネル farfield_flux_d に同じ入力を与え、残差・massflux・面の値をホスト計算と相対 1e-6 で照合
//   (v) 非対応構成の起動拒否は forge 本体の起動で確かめる (tests/unit/test_farfield_reject.sh)。
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cuda_runtime.h>
#include "flowFormat.hpp"
#include "cuda_forge/thermo_d.cuh"
#include "cuda_forge/convection/farfieldFlux_d.inc.cuh"

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++g_fail; printf("  FAIL: " __VA_ARGS__); printf("\n"); } } while (0)
template <class T> static T* up(const std::vector<T>& v) { T* d = nullptr; cudaMalloc((void**)&d, v.size() * sizeof(T)); cudaMemcpy(d, v.data(), v.size() * sizeof(T), cudaMemcpyHostToDevice); return d; }
template <class T> static std::vector<T> down(const T* d, size_t n) { std::vector<T> v(n); cudaMemcpy(v.data(), d, n * sizeof(T), cudaMemcpyDeviceToHost); return v; }

// 2 成分 (0 = EXH 相当、1 = AIR 相当)。cp が T に依存する NASA-9 (a2 定数 + a3 T) で TP の経路を通す
static SpeciesThermo mkSp(double MW, double a2, double a3)
{
    // 区間可変の NASA-9 (2026-10-01 #13-1) に合わせ、旧 low/high (同じ係数・Tmid 1000 K) と同値の 2 区間で組む。
    SpeciesThermo s{}; s.MW = MW; s.invMW = 1.0 / MW;
    const double Tb[3] = {200.0, 1000.0, 6000.0};
    double a[2][9] = {};
    a[0][2] = a[1][2] = a2; a[0][3] = a[1][3] = a3;
    thermo_set_intervals(s, 2, Tb, a);
    return s;
}
static SpeciesThermo g_sp[2] = { mkSp(0.0250, 4.0, 2.0e-4), mkSp(0.02896, 3.5, 1.0e-4) };

static FfGas gasCPG() { FfGas g{}; g.thermalMethod = 0; g.ga = 1.4; g.cp = 1004.5; g.sp = nullptr; g.nSpecies = 0; return g; }
static FfGas gasTP()  { FfGas g{}; g.thermalMethod = 2; g.ga = 1.4; g.cp = 1004.5; g.sp = g_sp; g.nSpecies = 2; return g; }

// 状態を (T, P, 速度, Y_EXH, k, ω) から作る
static FfPrim mk(const FfGas& g, double T, double P, double ux, double uy, double uz, double yExh, double k = 0.0, double om = 0.0)
{
    FfPrim s{};
    s.p = P; s.u[0] = ux; s.u[1] = uy; s.u[2] = uz; s.k = k; s.om = om;
    s.Y[0] = yExh; s.Y[1] = 1.0 - yExh;
    const double R = (g.thermalMethod == 2) ? thermo_R_mix(g.sp, 2, s.Y) : g.cp - g.cp / g.ga;
    s.r = P / (R * T);
    ff_thermo(g, s);
    return s;
}
static double un(const FfPrim& s, const double n[3]) { return s.u[0] * n[0] + s.u[1] * n[1] + s.u[2] * n[2]; }
// 物理流束 (単位面積)
static void physFlux(const FfPrim& s, const double n[3], int energyK, double F[5])
{
    const double q = un(s, n);
    const double p = s.p + (energyK ? 2.0 / 3.0 * s.r * s.k : 0.0);
    const double E = s.e + 0.5 * (s.u[0] * s.u[0] + s.u[1] * s.u[1] + s.u[2] * s.u[2]) + (energyK ? s.k : 0.0);
    F[0] = s.r * q;
    for (int d = 0; d < 3; ++d) F[1 + d] = s.r * s.u[d] * q + p * n[d];
    F[4] = (s.r * E + p) * q;
}
// 流束の規模 (保存量別): 質量 ρ(|u|+c)、運動量 ρ(|u|+c)² + P、エネルギー ρ(|u|+c)(|e| + (|u|+c)²)
static void scale(const FfPrim& a, const FfPrim& b, double sc[5])
{
    for (int q = 0; q < 5; ++q) sc[q] = 0.0;
    for (const FfPrim* s : { &a, &b }) {
        const double v = sqrt(s->u[0] * s->u[0] + s->u[1] * s->u[1] + s->u[2] * s->u[2]) + s->c;
        sc[0] = fmax(sc[0], s->r * v);
        for (int d = 1; d <= 3; ++d) sc[d] = fmax(sc[d], s->r * v * v + s->p);
        sc[4] = fmax(sc[4], s->r * v * (fabs(s->e) + v * v));
    }
}
static double relErr(const double A[5], const double B[5], const double sc[5])
{
    double m = 0.0;
    for (int q = 0; q < 5; ++q) m = fmax(m, fabs(A[q] - B[q]) / sc[q]);
    return m;
}
static void unitN(double n[3], double a, double b, double c) { const double s = sqrt(a * a + b * b + c * c); n[0] = a / s; n[1] = b / s; n[2] = c / s; }

// 1 面の評価 (外側状態 + HLLC)
struct Eval { FfPrim R; int vac, hll; double F[5]; };
static Eval evalFace(const FfGas& g, const FfPrim& i, const FfPrim& f, const double n[3], int energyK = 0)
{
    Eval e{};
    e.vac = ff_outer_state(g, i, f, n, e.R);
    e.hll = ff_hllc(i, e.R, n, 0.0, energyK, e.F);
    return e;
}

// ---------------------------------------------------------------------------------------------------------
static void test_uniform()
{
    printf("[V0u (iv)] 一様状態で数値流束 = 物理流束\n");
    double n[3]; unitN(n, 0.6, -0.48, 0.64);
    const double Ms[] = { 0.5, -0.5, 2.0, -2.0, 0.0, 0.95, -0.95, 1.0, -1.0 };
    for (int gi = 0; gi < 2; ++gi) {
        const FfGas g = gi ? gasTP() : gasCPG();
        double worst = 0.0; int bad = 0;
        for (double M : Ms) {
            FfPrim s = mk(g, 400.0, 3.0e4, 0.0, 0.0, 0.0, 0.3, 50.0, 1.0e4);
            // 法線成分 M c と接線成分 0.7 c
            double t[3] = { n[1], -n[0], 0.0 }; const double tl = sqrt(t[0] * t[0] + t[1] * t[1]); for (double& x : t) x /= tl;
            for (int d = 0; d < 3; ++d) s.u[d] = M * s.c * n[d] + 0.7 * s.c * t[d];
            ff_thermo(g, s);
            for (int eK = 0; eK < 2; ++eK) {
                Eval e = evalFace(g, s, s, n, eK);
                double Fp[5]; physFlux(s, n, eK, Fp);
                double sc[5]; scale(s, s, sc);
                const double r = relErr(e.F, Fp, sc);
                worst = fmax(worst, r);
                if (!(r <= 1e-6) || e.vac || e.hll) { ++bad; printf("    M %g eK %d: rel %.3e vac %d hll %d\n", M, eK, r, e.vac, e.hll); }
            }
        }
        printf("  %s: 最大相対差 %.3e\n", gi ? "TP" : "CPG", worst);
        CHECK(bad == 0, "(iv) %s 一様で %d 件が不一致", gi ? "TP" : "CPG", bad);
    }
}

static void test_contact()
{
    printf("[V0u (ii)] 接触波 (P・U_n 同じ、T 600/220 K、Y 0.13/0)\n");
    double n[3]; unitN(n, 1.0, 0.2, -0.1);
    for (int gi = 0; gi < 2; ++gi) {
        const FfGas g = gi ? gasTP() : gasCPG();
        for (int dir = -1; dir <= 1; dir += 2) {
            FfPrim i = mk(g, 600.0, 2851.0, 0, 0, 0, 0.13);
            FfPrim f = mk(g, 220.0, 2851.0, 0, 0, 0, 0.0);
            const double q = dir * 0.3 * f.c;   // 両側とも亜音速
            for (int d = 0; d < 3; ++d) { i.u[d] = q * n[d]; f.u[d] = q * n[d]; }
            Eval e = evalFace(g, i, f, n);
            double Fp[5]; physFlux(dir > 0 ? i : f, n, 0, Fp);
            double sc[5]; scale(i, f, sc);
            const double r = relErr(e.F, Fp, sc);
            printf("  %s %s: rel %.3e (Y_R %.4f)\n", gi ? "TP" : "CPG", dir > 0 ? "流出" : "流入", r, e.R.Y[0]);
            CHECK(r <= 1e-6 && !e.vac && !e.hll, "(ii) %s dir %d rel %.3e", gi ? "TP" : "CPG", dir, r);
            if (gi) CHECK(fabs(e.R.Y[0] - 0.0) < 1e-12, "(ii) TP 外側組成が自由流でない (%.3e)", e.R.Y[0]);
        }
    }
}

// 1 つの量を掃引して隣接差を見る
template <class Build>
static double sweepMaxJump(const FfGas& g, Build build, double x0, double dx, int N, const double n[3], int energyK, int& nVac, int& nHll)
{
    double worst = 0.0; double Fprev[5]; bool have = false; double sc[5];
    for (int k = -N; k <= N; ++k) {
        FfPrim i, f; build(x0 + k * dx, i, f);
        Eval e = evalFace(g, i, f, n, energyK);
        nVac += e.vac ? 1 : 0; nHll += e.hll ? 1 : 0;
        if (have) worst = fmax(worst, relErr(e.F, Fprev, sc));
        scale(i, f, sc);
        for (int q = 0; q < 5; ++q) Fprev[q] = e.F[q];
        have = true;
    }
    return worst;
}

static void test_continuity()
{
    printf("[V0u (i)] 連続性 (1e-6 刻み、隣接差 ≤ 1e-4 × 規模)\n");
    double n[3]; unitN(n, 0.8, 0.0, 0.6);
    const double centers[] = { 0.0, 0.9, 1.0, -0.9, -1.0, 1.1, -1.1 };
    for (int gi = 0; gi < 2; ++gi) {
        const FfGas g = gi ? gasTP() : gasCPG();
        const FfPrim i0 = mk(g, 700.0, 4.0e3, 0, 0, 0, 0.13, 20.0, 5e3);
        const FfPrim f0 = mk(g, 220.0, 2851.0, 0, 0, 0, 0.0, 1.0, 1e3);
        double worst = 0.0; int nVac = 0, nHll = 0;
        for (double Mc : centers) {
            // (a) 内部の法線速度 (自由流は亜音速で静止に近い)
            auto bUi = [&](double M, FfPrim& i, FfPrim& f) {
                i = i0; f = f0; for (int d = 0; d < 3; ++d) i.u[d] = M * i0.c * n[d] + (d == 1 ? 0.3 * i0.c : 0.0);
                for (int d = 0; d < 3; ++d) f.u[d] = 0.2 * f0.c * n[d]; ff_thermo(g, i); ff_thermo(g, f); };
            // (b) 自由流の法線速度 (内部は亜音速流出)
            auto bUf = [&](double M, FfPrim& i, FfPrim& f) {
                i = i0; f = f0; for (int d = 0; d < 3; ++d) i.u[d] = 0.3 * i0.c * n[d];
                for (int d = 0; d < 3; ++d) f.u[d] = M * f0.c * n[d] + (d == 1 ? 0.5 * f0.c : 0.0); ff_thermo(g, i); ff_thermo(g, f); };
            for (int eK = 0; eK < 2; ++eK) {
                worst = fmax(worst, sweepMaxJump(g, bUi, Mc, 1e-6, 200, n, eK, nVac, nHll));
                worst = fmax(worst, sweepMaxJump(g, bUf, Mc, 1e-6, 200, n, eK, nVac, nHll));
            }
            // (c) 内部の密度・圧力 (法線 Mach を Mc に置いたまま)
            auto bRho = [&](double s, FfPrim& i, FfPrim& f) {
                i = i0; f = f0; i.r = i0.r * s; for (int d = 0; d < 3; ++d) i.u[d] = Mc * i0.c * n[d]; ff_thermo(g, i); ff_thermo(g, f); };
            auto bP = [&](double s, FfPrim& i, FfPrim& f) {
                i = i0; f = f0; i.p = i0.p * s; for (int d = 0; d < 3; ++d) i.u[d] = Mc * i0.c * n[d]; ff_thermo(g, i); ff_thermo(g, f); };
            worst = fmax(worst, sweepMaxJump(g, bRho, 1.0, 1e-6, 200, n, 0, nVac, nHll));
            worst = fmax(worst, sweepMaxJump(g, bP, 1.0, 1e-6, 200, n, 0, nVac, nHll));
        }
        // (d) 混合帯の逆流の反例状態 (iii''') の近傍
        auto bMix = [&](double M, FfPrim& i, FfPrim& f) {
            i = mk(g, 600.0, 1.0e5, 0, 0, 0, 0.13, 0.14, 1.0); f = mk(g, 220.0, 1.0e6, 0, 0, 0, 0.0, 0.01, 1.0);
            for (int d = 0; d < 3; ++d) { i.u[d] = M * i.c * n[d]; f.u[d] = -3.0 * f.c * n[d]; } ff_thermo(g, i); ff_thermo(g, f); };
        worst = fmax(worst, sweepMaxJump(g, bMix, 0.95, 1e-6, 200, n, 1, nVac, nHll));
        printf("  %s: 最大隣接差 %.3e (置換 %d、HLL 退避 %d)\n", gi ? "TP" : "CPG", worst, nVac, nHll);
        CHECK(worst <= 1e-4, "(i) %s 連続性 %.3e", gi ? "TP" : "CPG", worst);
        CHECK(nVac == 0 && nHll == 0, "(i) %s 掃引中に置換 %d / 退避 %d", gi ? "TP" : "CPG", nVac, nHll);
    }
}

static void test_direction()
{
    printf("[V0u (iii)] 流向の組合せ: 帯の外で外側組成・k・ω = 自由流、常に有界\n");
    double n[3]; unitN(n, 0.0, 1.0, 0.0);
    const FfGas g = gasTP();
    int bad = 0, nIn = 0, nOut = 0;
    const double Mi[] = { -1.5, -0.95, -0.5, -1e-4, 0.0, 1e-4, 0.5, 0.95, 1.5 };
    const double Mf[] = { -1.5, -0.95, -0.5, -1e-4, 0.0, 1e-4, 0.5, 0.95, 1.5 };
    const double Pr[] = { 0.5, 1.0, 2.0 };
    for (double a : Mi) for (double b : Mf) for (double pr : Pr) {
        FfPrim i = mk(g, 600.0, 2851.0 * pr, 0, 0, 0, 0.13, 0.14, 2.0);
        FfPrim f = mk(g, 220.0, 2851.0, 0, 0, 0, 0.0, 0.01, 1.0);
        for (int d = 0; d < 3; ++d) { i.u[d] = a * i.c * n[d]; f.u[d] = b * f.c * n[d]; }
        i.u[0] = 0.2 * i.c; f.u[0] = 0.4 * f.c;
        ff_thermo(g, i); ff_thermo(g, f);
        Eval e = evalFace(g, i, f, n, 1);
        (e.F[0] < 0.0 ? nIn : nOut)++;
        const bool outBand = (a < 0.9) && (b > -0.9);
        auto within = [](double v, double x, double y) { return v >= fmin(x, y) - 1e-14 && v <= fmax(x, y) + 1e-14; };
        bool ok = within(e.R.Y[0], i.Y[0], f.Y[0]) && within(e.R.k, i.k, f.k) && within(e.R.om, i.om, f.om);
        if (outBand) ok = ok && e.R.Y[0] == f.Y[0] && e.R.k == f.k && e.R.om == f.om && e.R.u[0] == f.u[0];
        if (a >= 1.0) ok = ok && e.R.Y[0] == i.Y[0] && e.R.k == i.k;          // 内部の超音速流出 = 内部 (優先)
        if (b <= -1.0 && a < 0.9) ok = ok && e.R.p == f.p && e.R.r == f.r;      // 自由流の超音速流入 = 自由流
        if (!ok || e.vac || e.hll) { ++bad; printf("    Mi %g Mf %g pr %g: Y_R %.4f k_R %.4f vac %d hll %d\n", a, b, pr, e.R.Y[0], e.R.k, e.vac, e.hll); }
    }
    printf("  組合せ %d 件 (ṁ<0 %d、ṁ≥0 %d)\n", nIn + nOut, nIn, nOut);
    CHECK(bad == 0, "(iii) %d 件が規則に反する", bad);
    CHECK(nIn > 0 && nOut > 0, "(iii) 流入・流出の両方が現れていない");
}

static void test_vacuum_float()
{
    printf("[V0u (iii'')] 真空と float32 の独立掃引\n");
    double n[3]; unitN(n, 1.0, 0.0, 0.0);
    for (int gi = 0; gi < 2; ++gi) {
        const FfGas g = gi ? gasTP() : gasCPG();
        for (int fl = 0; fl < 2; ++fl) {
            long nAll = 0, nNonfin = 0, nNeg = 0, nVac = 0, nHll = 0;
            for (int a = 0; a <= 24; ++a) for (int b = 0; b <= 26; ++b)
            for (int rr = -2; rr <= 2; ++rr) for (int pr = -2; pr <= 2; ++pr) {
                const double Mi = -3.0 + 0.25 * a, Mf = -3.0 + 0.5 * b;
                FfPrim i = mk(g, 400.0, 3.0e3 * pow(10.0, pr), 0, 0, 0, 0.13);
                i.r *= pow(10.0, 0.5 * rr);
                ff_thermo(g, i);
                if (gi && !(i.T > 200.0 && i.T < 6000.0)) continue;   // NASA の定義域の外は掃引しない
                FfPrim f = mk(g, 220.0, 3.0e3, 0, 0, 0, 0.0);
                i.u[0] = Mi * i.c; f.u[0] = Mf * f.c; i.u[1] = 0.3 * i.c;
                if (fl) {   // float32 の入力 (ソルバの保存場)
                    i.r = (float)i.r; i.p = (float)i.p; for (double& x : i.u) x = (float)x;
                    f.r = (float)f.r; f.p = (float)f.p; for (double& x : f.u) x = (float)x;
                }
                ff_thermo(g, i); ff_thermo(g, f);
                Eval e = evalFace(g, i, f, n);
                ++nAll;
                bool fin = true; for (double x : e.F) fin = fin && std::isfinite(x);
                nNonfin += fin ? 0 : 1;
                nNeg += (e.R.r > 0.0 && e.R.p > 0.0 && e.R.c > 0.0) ? 0 : 1;
                nVac += e.vac ? 1 : 0; nHll += e.hll ? 1 : 0;
            }
            printf("  %s %s: %ld 件、非有限 %ld、負値 %ld、真空置換 %ld、HLL 退避 %ld\n",
                   gi ? "TP" : "CPG", fl ? "float32" : "float64", nAll, nNonfin, nNeg, nVac, nHll);
            CHECK(nNonfin == 0 && nNeg == 0, "(iii'') %s %s 非有限 %ld / 負値 %ld", gi ? "TP" : "CPG", fl ? "f32" : "f64", nNonfin, nNeg);
        }
    }
    // 真空の反例: 強い 2 膨張 (内部が外向きに、自由流が離れる向きに速い)
    const FfGas g = gasCPG();
    FfPrim i = mk(g, 300.0, 1.0e3, 0, 0, 0, 0.0), f = mk(g, 300.0, 1.0e3, 0, 0, 0, 0.0);
    i.u[0] = -8.0 * i.c; f.u[0] = 8.0 * f.c; ff_thermo(g, i); ff_thermo(g, f);
    Eval e = evalFace(g, i, f, n);
    bool fin = true; for (double x : e.F) fin = fin && std::isfinite(x);
    printf("  真空の反例: vac %d、HLL %d、有限 %d\n", e.vac, e.hll, (int)fin);
    CHECK(e.vac == 1 && fin, "(iii'') 真空の反例で置換されない (vac %d)", e.vac);
}

static void test_mixband_backflow()
{
    printf("[V0u (iii''')] 混合帯の逆流 (Y_R 0.065、k_R 0.075)\n");
    double n[3]; unitN(n, 1.0, 0.0, 0.0);
    for (int gi = 0; gi < 2; ++gi) {
        const FfGas g = gi ? gasTP() : gasCPG();
        FfPrim i, f;
        if (gi) { i = mk(g, 600.0, 1.0e5, 0, 0, 0, 0.13, 0.14, 1.0); f = mk(g, 220.0, 1.0e6, 0, 0, 0, 0.0, 0.01, 1.0); }
        else {  // plan の無次元状態: 内部 (1, 0.95, 1/γ)・自由流 (10, −3, 10/γ)
            i = FfPrim{}; i.r = 1.0; i.p = 1.0 / 1.4; i.Y[0] = 0.13; i.Y[1] = 0.87; i.k = 0.14; i.om = 1.0;
            f = FfPrim{}; f.r = 10.0; f.p = 10.0 / 1.4; f.Y[0] = 0.0; f.Y[1] = 1.0; f.k = 0.01; f.om = 1.0;
            ff_thermo(g, i); ff_thermo(g, f);
        }
        i.u[0] = 0.95 * i.c; f.u[0] = -3.0 * f.c; ff_thermo(g, i); ff_thermo(g, f);
        Eval e = evalFace(g, i, f, n, 1);
        const double yR = gi ? e.R.Y[0] : 0.065;   // CPG は組成を持たない (nY = 0)
        printf("  %s: Y_R %.6f k_R %.6f ṁ %.4e vac %d hll %d\n", gi ? "TP" : "CPG", e.R.Y[0], e.R.k, e.F[0], e.vac, e.hll);
        CHECK(fabs(yR - 0.065) < 1e-12 && fabs(e.R.k - 0.075) < 1e-12, "(iii''') %s Y_R %.6f k_R %.6f", gi ? "TP" : "CPG", e.R.Y[0], e.R.k);
        // SST k 込みの流束が k_R から作られている: R 側の k を別の値にすると流束が変わる (k_R を実際に読んでいる)
        FfPrim R2 = e.R; R2.k = f.k; double F2[5]; ff_hllc(i, R2, n, 0.0, 1, F2);
        double sc[5]; scale(i, f, sc);
        CHECK(relErr(e.F, F2, sc) > 1e-8, "(iii''') %s HLLC が k_R を読んでいない", gi ? "TP" : "CPG");
    }
}

// ---------------------------------------------------------------------------------------------------------
// V0k: 実カーネルとホストの照合
struct KCase { FfPrim i; double n[3]; double S; };
static void test_kernel()
{
    printf("[V0k] GPU カーネル farfield_flux_d とホスト計算の照合\n");
    for (int gi = 0; gi < 2; ++gi) {
        FfGas g = gi ? gasTP() : gasCPG();
        const int nY = gi ? 2 : 0;
        // 自由流 (M 2 で +x)
        FfPrim f = mk(g, 220.0, 2851.0, 0, 0, 0, 0.0, 0.01, 1.0e3);
        f.u[0] = 2.0 * f.c; f.u[1] = 0.1 * f.c; ff_thermo(g, f);
        std::vector<KCase> cs;
        const double Ms[] = { -1.5, -0.95, -0.3, 0.0, 0.3, 0.95, 1.0, 1.5 };
        const double dirs[][3] = { {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0.6,0.8,0}, {-0.6,0,0.8} };
        for (double M : Ms) for (auto& dv : dirs) {
            KCase k{}; unitN(k.n, dv[0], dv[1], dv[2]); k.S = 1.0e-4 * (1.0 + 0.1 * cs.size());
            k.i = mk(g, 400.0 + 10.0 * cs.size(), 2851.0 * (0.5 + 0.02 * cs.size()), 0, 0, 0, 0.13, 0.05, 2.0e3);
            for (int d = 0; d < 3; ++d) k.i.u[d] = M * k.i.c * k.n[d] + (d == 2 ? 0.2 * k.i.c : 0.0);
            ff_thermo(g, k.i);
            cs.push_back(k);
        }
        const int nb = (int)cs.size();
        // float の入力場 (ソルバと同じ型)。面 ib は節点 ib に 1 対 1
        std::vector<int> pl(nb), ce(nb);
        std::vector<float> sx(nb), sy(nb), sz(nb), ss(nb), ro(nb), ux(nb), uy(nb), uz(nb), P(nb), kk(nb), om(nb), Y0(nb), Y1(nb);
        for (int b = 0; b < nb; ++b) {
            pl[b] = b; ce[b] = b;
            sx[b] = (float)(cs[b].n[0] * cs[b].S); sy[b] = (float)(cs[b].n[1] * cs[b].S); sz[b] = (float)(cs[b].n[2] * cs[b].S); ss[b] = (float)cs[b].S;
            ro[b] = (float)cs[b].i.r; ux[b] = (float)cs[b].i.u[0]; uy[b] = (float)cs[b].i.u[1]; uz[b] = (float)cs[b].i.u[2]; P[b] = (float)cs[b].i.p;
            kk[b] = (float)cs[b].i.k; om[b] = (float)cs[b].i.om; Y0[b] = (float)cs[b].i.Y[0]; Y1[b] = (float)cs[b].i.Y[1];
        }
        FfInf inf{}; inf.r = f.r; for (int d = 0; d < 3; ++d) inf.u[d] = f.u[d]; inf.p = f.p; inf.k = f.k; inf.om = f.om; inf.Y[0] = f.Y[0]; inf.Y[1] = f.Y[1];
        FfGas gd = g;
        SpeciesThermo* dsp = nullptr;
        if (gi) { dsp = up(std::vector<SpeciesThermo>(g_sp, g_sp + 2)); gd.sp = dsp; }
        int *dpl = up(pl), *dce = up(ce);
        float *dsx = up(sx), *dsy = up(sy), *dsz = up(sz), *dss = up(ss), *dro = up(ro), *dux = up(ux), *duy = up(uy), *duz = up(uz), *dP = up(P), *dk = up(kk), *dom = up(om);
        float *dY0 = up(Y0), *dY1 = up(Y1);
        std::vector<float> zero(nb, 0.0f);
        float *r0 = up(zero), *r1 = up(zero), *r2 = up(zero), *r3 = up(zero), *r4 = up(zero), *mf = up(zero), *fY0 = up(zero), *fY1 = up(zero), *fK = up(zero), *fO = up(zero);
        float** dRoY = up(std::vector<float*>{ dY0, dY1 });
        float** dFY = up(std::vector<float*>{ fY0, fY1 });
        farfield_flux_d<<<(nb + 127) / 128, 128>>>(nb, dpl, dce, dsx, dsy, dsz, dss, dro, dux, duy, duz, dP,
            nY ? dRoY : nullptr, dk, dom, gd, inf, (flow_float)0.0f, 1, r0, r1, r2, r3, r4, mf, nY ? dFY : nullptr, fK, fO, nullptr);
        const cudaError_t err = cudaDeviceSynchronize();
        CHECK(err == cudaSuccess, "V0k カーネル実行エラー %s", cudaGetErrorString(err));
        auto R0 = down(r0, nb), R1 = down(r1, nb), R2 = down(r2, nb), R3 = down(r3, nb), R4 = down(r4, nb), MF = down(mf, nb);
        auto FY0 = down(fY0, nb), FK = down(fK, nb), FO = down(fO, nb);
        unsigned long long hv = 0, hh = 0;
        cudaMemcpyFromSymbol(&hv, g_ffVac, sizeof(hv)); cudaMemcpyFromSymbol(&hh, g_ffHll, sizeof(hh));
        // ホストで同じ手順 (カーネルと同じ float 入力から)
        double worst = 0.0, worstS = 0.0; int bad = 0;
        for (int b = 0; b < nb; ++b) {
            const double S = ss[b]; const double n[3] = { sx[b] / S, sy[b] / S, sz[b] / S };
            FfPrim I{}; I.r = ro[b]; I.u[0] = ux[b]; I.u[1] = uy[b]; I.u[2] = uz[b]; I.p = P[b]; I.k = kk[b]; I.om = om[b];
            if (nY) { const double s0 = fmax((double)Y0[b], 0.0), s1 = fmax((double)Y1[b], 0.0), sm = s0 + s1; I.Y[0] = s0 / sm; I.Y[1] = s1 / sm; }
            ff_thermo(g, I);
            FfPrim Fs{}; Fs.r = inf.r; for (int d = 0; d < 3; ++d) Fs.u[d] = inf.u[d]; Fs.p = inf.p; Fs.k = inf.k; Fs.om = inf.om; Fs.Y[0] = inf.Y[0]; Fs.Y[1] = inf.Y[1];
            ff_thermo(g, Fs);
            Eval e = evalFace(g, I, Fs, n, 1);
            double sc[5]; scale(I, Fs, sc);
            const double dev[5] = { -R0[b] / S, -R1[b] / S, -R2[b] / S, -R3[b] / S, -R4[b] / S };
            const double r = relErr(dev, e.F, sc);
            const double rm = fabs(MF[b] / S - e.F[0]) / sc[0];
            double rs = fabs(FK[b] - e.R.k) / fmax(fabs(e.R.k), 1e-30) + fabs(FO[b] - e.R.om) / fmax(fabs(e.R.om), 1e-30);
            if (nY) rs += fabs(FY0[b] - e.R.Y[0]);
            worst = fmax(worst, fmax(r, rm)); worstS = fmax(worstS, rs);
            if (!(r <= 1e-6 && rm <= 1e-6 && rs <= 1e-6)) { ++bad; printf("    面 %d: 流束 %.3e ṁ %.3e 面値 %.3e (vac %d)\n", b, r, rm, rs, e.vac); }
        }
        printf("  %s: %d 面、流束の最大相対差 %.3e、面の値 %.3e、GPU の置換 %llu・退避 %llu\n", gi ? "TP" : "CPG", nb, worst, worstS, hv, hh);
        CHECK(bad == 0, "V0k %s %d 面が不一致", gi ? "TP" : "CPG", bad);
        CHECK(hv == 0 && hh == 0, "V0k %s GPU で置換 %llu / 退避 %llu", gi ? "TP" : "CPG", hv, hh);
        const unsigned long long z = 0; cudaMemcpyToSymbol(g_ffVac, &z, sizeof(z)); cudaMemcpyToSymbol(g_ffHll, &z, sizeof(z));
        for (void* p : std::vector<void*>{ dpl, dce, dsx, dsy, dsz, dss, dro, dux, duy, duz, dP, dk, dom, dY0, dY1, r0, r1, r2, r3, r4, mf, fY0, fY1, fK, fO, dRoY, dFY, dsp }) cudaFree(p);
    }
}

int main()
{
    test_uniform();
    test_contact();
    test_continuity();
    test_direction();
    test_vacuum_float();
    test_mixband_backflow();
    test_kernel();
    printf(g_fail ? "RESULT: FAIL (%d)\n" : "RESULT: PASS\n", g_fail);
    return g_fail ? 1 : 0;
}
