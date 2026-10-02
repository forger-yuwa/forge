// test_twophase_kernel.cu — 二相拡散の CUDA 実装 (cuda_forge/twoPhaseDiffusion_d.cuh) の GPU 単体試験
//   plans/active/condensation-two-phase-transport.md §5.1 #4e。受入は notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md §14.2
//   (実行前に固定): T1 面流束 (GPU vs 同じ float 入力の double 参照)、T2 §6 の単体条件、T3 #4c/#4d の 1D 問題で GPU とホストの一致。
//   面の代数 tp_face_flux と 1 セルの更新 tp_vl_update は本番カーネルと同じ関数 (atomicAdd の前の値を probe で書き出す)。
//
// build (GPU 1 本; test_twophase_real_source.cpp と同じ埋め込みヘッダと種 DB の .o を使う):
//   cp tests/unit/test_twophase_kernel.cu <OUT>/ ; cd solver_density_cuda
//   nice -n 19 nvcc -O2 -arch=sm_86 --expt-relaxed-constexpr -std=c++17 -I . -I <GEN> -o <OUT>/test_twophase_kernel \
//       tests/unit/test_twophase_kernel.cu <OUT>/speciesDB.o <OUT>/speciesTransportDB.o -lyaml-cpp
//   <OUT>/test_twophase_kernel [--quick]
#define main twophase_real_source_main_unused
#include "tests/unit/test_twophase_real_source.cpp"
#undef main
#include "cuda_forge/twoPhaseDiffusion_d.cuh"
#include "cuda_forge/passiveKernels_d.cuh"          // passive_bounds_d (#4h の下流の床)
#include "cuda_forge/condensationRealizability_d.cuh"   // cond_realizability_clamp_f_d (#4h の液の射影)
#include <random>

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); exit(2); } } while (0)

// ------------------------------------------------------------------ device probes (atomicAdd の前の値)
// 本番カーネル twophase_diffusion_d と同じ式で面の h_k(T_f)・L(T_f) を作る (T_f = f T0 + (1−f) T1 を float で)。
__global__ void k_face(int nF, TpFaceIn* in, const float* T0, const float* T1, const SpeciesThermoF* sp, CondSpeciesProps cp, int setProps, TpFaceOut* out)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= nF) return;
    TpFaceIn a = in[i];
    if (setProps) {
        const float g = 1.0f - a.f;
        const float Tf = a.f*T0[i] + g*T1[i];
        for (int s = 0; s < a.n; ++s) a.h[s] = thermo_h_mass_f(sp[s], Tf);
        a.L = (float)cond_latent(cp, (double)Tf);
        in[i] = a;   // 参照側が同じ h・L を読めるように書き戻す (h・L の評価誤差は別に見る)
    }
    tp_face_flux(a, out[i]);
}
__global__ void k_update(int n, const TpCellIn* in, TpCellOut* out, int thetaRound)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n) tp_vl_update(in[i], out[i], thetaRound);
}
// 液の更新値を bounds の前で見るための 1 セル版 (#4f (1)): 更新本体と同じ関数を通し、θ の float 値も返す。
__global__ void k_update_probe(TpCellIn in, int thetaRound, TpCellOut* out, float* Th)
{
    tp_vl_update(in, *out, thetaRound);
    // 実効 θ (float) を同じ規則で再現して返す
    float t = (float)out->theta;
    if (thetaRound != 0 && (double)t > out->theta) t = nextafterf(t, 0.0f);
    *Th = t;
}

struct Gpu {
    SpeciesThermoF* sp = nullptr; int nsp = 0; CondSpeciesProps cp{};
    TpFaceIn* fin = nullptr; TpFaceOut* fout = nullptr; float *T0 = nullptr, *T1 = nullptr; int capF = 0;
    TpCellIn* cin = nullptr; TpCellOut* cout_ = nullptr; int capC = 0;
    void faces(std::vector<TpFaceIn>& in, const std::vector<float>& t0, const std::vector<float>& t1, bool props, std::vector<TpFaceOut>& out) {
        const int n = (int)in.size();
        if (n > capF) { cudaFree(fin); cudaFree(fout); cudaFree(T0); cudaFree(T1); CK(cudaMalloc(&fin, n*sizeof(TpFaceIn))); CK(cudaMalloc(&fout, n*sizeof(TpFaceOut)));
                        CK(cudaMalloc(&T0, n*sizeof(float))); CK(cudaMalloc(&T1, n*sizeof(float))); capF = n; }
        CK(cudaMemcpy(fin, in.data(), n*sizeof(TpFaceIn), cudaMemcpyHostToDevice));
        if (props) { CK(cudaMemcpy(T0, t0.data(), n*sizeof(float), cudaMemcpyHostToDevice)); CK(cudaMemcpy(T1, t1.data(), n*sizeof(float), cudaMemcpyHostToDevice)); }
        k_face<<<(n + 127)/128, 128>>>(n, fin, T0, T1, sp, cp, props ? 1 : 0, fout);
        CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
        out.resize(n);
        CK(cudaMemcpy(out.data(), fout, n*sizeof(TpFaceOut), cudaMemcpyDeviceToHost));
        if (props) CK(cudaMemcpy(in.data(), fin, n*sizeof(TpFaceIn), cudaMemcpyDeviceToHost));
    }
    void cells(const std::vector<TpCellIn>& in, std::vector<TpCellOut>& out) {
        const int n = (int)in.size();
        if (n > capC) { cudaFree(cin); cudaFree(cout_); CK(cudaMalloc(&cin, n*sizeof(TpCellIn))); CK(cudaMalloc(&cout_, n*sizeof(TpCellOut))); capC = n; }
        CK(cudaMemcpy(cin, in.data(), n*sizeof(TpCellIn), cudaMemcpyHostToDevice));
        k_update<<<(n + 127)/128, 128>>>(n, cin, cout_, TP_THETA_ROUND_DEFAULT);
        CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
        out.resize(n);
        CK(cudaMemcpy(out.data(), cout_, n*sizeof(TpCellOut), cudaMemcpyDeviceToHost));
    }
};
static Gpu G;
static SpeciesThermo g_sp3[3];   // T1 の 3 種 (N2, O2, H2O) の double 係数 (h の評価差の情報用)

// ------------------------------------------------------------------ 参照 (C++ 1D ハーネスの拡散式を一般化: 向き L→R が正、Σz 正規化)
struct RefFace {
    double J[THERMO_MAX_SPECIES], Jv, Jl, JQ[TP_NQ], q, diag0[THERMO_MAX_SPECIES], diag1[THERMO_MAX_SPECIES], diagt0, diagt1;
    double j0[THERMO_MAX_SPECIES], jc[THERMO_MAX_SPECIES];   // 分子 (補正前・補正後; forge 向き)
    double tol[THERMO_MAX_SPECIES], tolv, toll, tolQ[TP_NQ], tolq;   // §14.2 T1 の許容
    double zL[THERMO_MAX_SPECIES], zR[THERMO_MAX_SPECIES], rgf;
};
// 入力は TpFaceIn (float) を double に上げたもの。h・L も in のもの (device が作った値) を使う → 演算誤差だけを見る。
static void ref_face(const TpFaceIn& in, RefFace& r)
{
    const int n = in.n, iw = in.iw;
    const double rhoL = in.rho0, rhoR = in.rho1, gL = in.rg0, gR = in.rg1;
    const double rgL = std::max(rhoL - gL, 1e-30), rgR = std::max(rhoR - gR, 1e-30);
    double sL = 0, sR = 0;
    for (int k = 0; k < n; ++k) {
        r.zL[k] = ((k == iw) ? ((double)in.rY0[k] - gL) : (double)in.rY0[k])/rgL;
        r.zR[k] = ((k == iw) ? ((double)in.rY1[k] - gR) : (double)in.rY1[k])/rgR;
        sL += r.zL[k]; sR += r.zR[k];
    }
    for (int k = 0; k < n; ++k) { r.zL[k] /= sL; r.zR[k] /= sR; }
    const double f = in.f, a = in.geo, aa = in.geo_abs, Gam = (double)in.ct;
    const double rgf = f*rgL + (1 - f)*rgR; r.rgf = rgf;
    // ハーネスの形 (物理の向き L→R): j0 = −ρ_g,f D a (z_R − z_L)、S = Σ j0、upL = (S ≤ 0)、j = j0 − z_up S
    double j0p[THERMO_MAX_SPECIES], S = 0;
    for (int k = 0; k < n; ++k) { j0p[k] = -rgf*(double)in.D[k]*a*(r.zR[k] - r.zL[k]); S += j0p[k]; }
    const bool upL = (S <= 0);
    auto Jt = [&](double phiL, double phiR) { return -Gam*a*(phiR - phiL); };
    double dif[THERMO_MAX_SPECIES];
    double AS = 0, Amol[THERMO_MAX_SPECIES];
    for (int k = 0; k < n; ++k) { Amol[k] = rgf*(double)in.D[k]*a*(fabs(r.zL[k]) + fabs(r.zR[k])); AS += Amol[k]; }
    const double E = EPS32;
    for (int k = 0; k < n; ++k) {
        const double zu = upL ? r.zL[k] : r.zR[k];
        const double jp = j0p[k] - zu*S;
        r.j0[k] = -j0p[k]; r.jc[k] = -jp;
        if (k == iw) {
            const double vL = ((double)in.rY0[k] - gL)/rhoL, vR = ((double)in.rY1[k] - gR)/rhoR;
            dif[k] = jp + Jt(vL, vR);
            r.Jv = -dif[k];
            r.tolv = 16*E*(Amol[k] + fabs(zu)*AS + Gam*a*(fabs(vL) + fabs(vR)));
        } else {
            const double yL = in.rY0[k]/rhoL, yR = in.rY1[k]/rhoR;
            dif[k] = jp + Jt(yL, yR);
            r.J[k] = -dif[k];
            r.tol[k] = 16*E*(Amol[k] + fabs(zu)*AS + Gam*a*(fabs(yL) + fabs(yR)));
        }
    }
    const double difl = Jt(gL/rhoL, gR/rhoR);
    r.Jl = -difl;
    r.toll = 16*E*Gam*a*(fabs(gL/rhoL) + fabs(gR/rhoR));
    if (iw >= 0) { r.J[iw] = r.Jv + r.Jl; r.tol[iw] = r.tolv + r.toll; }
    for (int m = 0; m < TP_NQ; ++m) {
        const double qL = in.rQ0[m]/rhoL, qR = in.rQ1[m]/rhoR;
        r.JQ[m] = -Jt(qL, qR);
        r.tolQ[m] = 16*E*Gam*a*(fabs(qL) + fabs(qR));
    }
    // エネルギー Σ_{k≠水} h_k J_k + h_v J_w − L J_l
    double q = 0, qa = 0, qt = 0;
    for (int k = 0; k < n; ++k) { q += (double)in.h[k]*r.J[k]; qa += fabs((double)in.h[k]*r.J[k]); qt += fabs((double)in.h[k])*r.tol[k]; }
    q -= (double)in.L*r.Jl; qa += fabs((double)in.L*r.Jl); qt += fabs((double)in.L)*r.toll;
    r.q = q;
    r.tolq = 4*qt + 64*E*qa;
    // 対角 (ハーネス: 分子 + 乱流、補正の流出側に |S|/ρ_g)
    for (int k = 0; k < n; ++k) {
        const double md = rgf*(double)in.D[k]*aa;
        r.diag0[k] = md/rgL + Gam*aa/rhoL; r.diag1[k] = md/rgR + Gam*aa/rhoR;
        if (upL) r.diag0[k] += fabs(S)/rgL; else r.diag1[k] += fabs(S)/rgR;
    }
    r.diagt0 = Gam*aa/rhoL; r.diagt1 = Gam*aa/rhoR;
}
static bool finite_out(const TpFaceOut& o, int n)
{
    bool ok = std::isfinite(o.Jv) && std::isfinite(o.Jl) && std::isfinite(o.q) && std::isfinite(o.diagt0) && std::isfinite(o.diagt1);
    for (int k = 0; k < n; ++k) ok = ok && std::isfinite(o.J[k]) && std::isfinite(o.diag0[k]) && std::isfinite(o.diag1[k]);
    for (int m = 0; m < TP_NQ; ++m) ok = ok && std::isfinite(o.JQ[m]);
    return ok;
}

// ------------------------------------------------------------------ T1 面流束 (受入 2)
static void test_T1(const Ctx& c, int nrep)
{
    printf("\n=== T1 面流束: GPU tp_face_flux と同じ float 入力の double 参照 (3 種 N2/O2/H2O、D_k 固定、乱流あり) ===\n");
    std::mt19937_64 rng(20261002);
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    for (double dvr : {1.0, 3.0}) {
        double worst[6] = {0, 0, 0, 0, 0, 0};   // 比 (species / v / l / Q / E / diag)
        long nF = 0; bool fin = true; long nup0 = 0;
        for (int rep = 0; rep < nrep; ++rep) {
            // 40 セルの状態 (入口状態のまわり)
            const int N = NC;
            std::vector<double> rho(N), g(N), Yw(N), YO2(N), T(N), q0(N), q1(N), q2(N);
            for (int i = 0; i < N; ++i) {
                rho[i] = RHO*(0.95 + 0.1*U01(rng));
                g[i] = (U01(rng) < 0.2) ? 0.0 : 5e-3*U01(rng);
                const double r = U01(rng);
                Yw[i] = (r < 0.1) ? g[i] : g[i] + (YW_IN - g[i] > 0 ? (YW_IN*1.2 - g[i])*U01(rng) : 1e-4*U01(rng));   // 1 割は蒸気 0
                YO2[i] = 0.2*(0.9 + 0.2*U01(rng));
                T[i] = 215.0 + 20.0*U01(rng);
                q0[i] = Q0_IN*pow(10.0, -1 + 2*U01(rng)); q1[i] = Q1_IN*pow(10.0, -1 + 2*U01(rng)); q2[i] = Q2_IN*pow(10.0, -1 + 2*U01(rng));
            }
            std::vector<TpFaceIn> in(N - 1); std::vector<float> t0(N - 1), t1(N - 1);
            for (int i = 0; i + 1 < N; ++i) {
                TpFaceIn& a = in[i]; memset(&a, 0, sizeof a);
                a.n = 3; a.iw = 2; a.f = 0.5f; a.geo = (float)(1.0/DX); a.geo_abs = a.geo;
                const int L = i, R = i + 1;
                a.rho0 = (float)rho[L]; a.rho1 = (float)rho[R];
                a.rY0[0] = (float)(rho[L]*(1 - Yw[L] - YO2[L])); a.rY0[1] = (float)(rho[L]*YO2[L]); a.rY0[2] = (float)(rho[L]*Yw[L]);
                a.rY1[0] = (float)(rho[R]*(1 - Yw[R] - YO2[R])); a.rY1[1] = (float)(rho[R]*YO2[R]); a.rY1[2] = (float)(rho[R]*Yw[R]);
                a.rg0 = (float)(rho[L]*g[L]); a.rg1 = (float)(rho[R]*g[R]);
                a.rQ0[0] = (float)(rho[L]*q2[L]); a.rQ0[1] = (float)(rho[L]*q1[L]); a.rQ0[2] = (float)(rho[L]*q0[L]);
                a.rQ1[0] = (float)(rho[R]*q2[R]); a.rQ1[1] = (float)(rho[R]*q1[R]); a.rQ1[2] = (float)(rho[R]*q0[R]);
                a.D[0] = (float)DMOL; a.D[1] = (float)(1.3*DMOL); a.D[2] = (float)(dvr*DMOL);
                a.ct = (float)(MUT/SCT);
                t0[i] = (float)T[L]; t1[i] = (float)T[R];
            }
            std::vector<TpFaceOut> out;
            G.faces(in, t0, t1, true, out);
            for (size_t i = 0; i < in.size(); ++i) {
                RefFace r; ref_face(in[i], r);
                const TpFaceOut& o = out[i];
                if (!finite_out(o, 3)) fin = false;
                nup0 += o.up0;
                for (int k = 0; k < 3; ++k) worst[0] = std::max(worst[0], fabs((double)o.J[k] - r.J[k])/std::max(r.tol[k], 1e-300));
                worst[1] = std::max(worst[1], fabs((double)o.Jv - r.Jv)/std::max(r.tolv, 1e-300));
                worst[2] = std::max(worst[2], fabs((double)o.Jl - r.Jl)/std::max(r.toll, 1e-300));
                for (int m = 0; m < TP_NQ; ++m) worst[3] = std::max(worst[3], fabs((double)o.JQ[m] - r.JQ[m])/std::max(r.tolQ[m], 1e-300));
                worst[4] = std::max(worst[4], fabs((double)o.q - r.q)/std::max(r.tolq, 1e-300));
                for (int k = 0; k < 3; ++k) {
                    worst[5] = std::max(worst[5], fabs((double)o.diag0[k] - r.diag0[k])/(16*EPS32*r.diag0[k]));
                    worst[5] = std::max(worst[5], fabs((double)o.diag1[k] - r.diag1[k])/(16*EPS32*r.diag1[k]));
                }
                worst[5] = std::max(worst[5], fabs((double)o.diagt0 - r.diagt0)/(16*EPS32*r.diagt0));
                worst[5] = std::max(worst[5], fabs((double)o.diagt1 - r.diagt1)/(16*EPS32*r.diagt1));
                ++nF;
            }
        }
        double wm = 0; for (double w : worst) wm = std::max(wm, w);
        char b[512];
        snprintf(b, sizeof b, "T1 D_v/D_N2 = %.0f: %ld 面 (補正 z をセル 0 から %ld 面)、|GPU − 参照|/許容 の最大: 種 %.3f・蒸気 %.3f・液 %.3f・Q %.3f・エネルギー %.3f・対角 %.3f、非有限 %s",
                 dvr, nF, nup0, worst[0], worst[1], worst[2], worst[3], worst[4], worst[5], fin ? "0" : "あり");
        verdict(fin && wm <= 1.0, b);
    }
    // h・L の評価 (float の thermo_h_mass_f・double の cond_latent) と double 参照の差 (情報)
    {
        std::vector<TpFaceIn> in(1); memset(&in[0], 0, sizeof(TpFaceIn)); in[0].n = 3; in[0].iw = 2; in[0].f = 0.5f; in[0].rho0 = in[0].rho1 = 1; in[0].rY0[0] = in[0].rY1[0] = 1;
        double eh = 0, eL = 0;
        for (double Tq : {215.0, 221.35, 230.0, 250.0, 290.0}) {
            std::vector<float> t0(1, (float)Tq), t1(1, (float)Tq); std::vector<TpFaceOut> out;
            G.faces(in, t0, t1, true, out);
            const double Tf = (double)(float)Tq;
            for (int k = 0; k < 3; ++k) eh = std::max(eh, fabs((double)in[0].h[k] - thermo_h_mass(g_sp3[k], Tf))/std::max(fabs(thermo_h_mass(g_sp3[k], Tf)), 1.0));
            eL = std::max(eL, fabs((double)in[0].L - cond_latent(c.cp, Tf))/cond_latent(c.cp, Tf));
        }
        printf("[INFO] T1 面の h_k (float) と double の相対差 最大 %.2e、L (double 評価の float 化) %.2e\n", eh, eL);
    }
}

// ------------------------------------------------------------------ T2 §6 の単体条件 (受入 1)
static void face_set_gas(TpFaceIn& a, double rho, double gL, double gR, const double* zg, int nG, double dx)
{
    memset(&a, 0, sizeof a);
    a.n = nG; a.iw = nG - 1; a.f = 0.5f; a.geo = (float)(1.0/dx); a.geo_abs = a.geo;
    a.rho0 = a.rho1 = (float)rho;
    for (int k = 0; k < nG; ++k) {
        a.rY0[k] = (float)((1 - gL)*zg[k]*rho + ((k == nG - 1) ? gL*rho : 0.0));
        a.rY1[k] = (float)((1 - gR)*zg[k]*rho + ((k == nG - 1) ? gR*rho : 0.0));
    }
    a.rg0 = (float)(gL*rho); a.rg1 = (float)(gR*rho);
}

static void test_T2_structure()
{
    printf("\n=== T2 (i) 分子流束の構造: 三成分 (等分子量)・気相一様・g のみ勾配 (∇g = 0.1)、D_k 固定 [1,2,3]e-5 ===\n");
    const double zg[3] = {0.2, 0.3, 0.5}, Dk[3] = {1e-5, 2e-5, 3e-5};
    for (double dx : {1.0, 1e-3}) {
        const double gL = 0.1, gR = gL + 0.1*dx, rho = 1.0;
        std::vector<TpFaceIn> in(1); face_set_gas(in[0], rho, gL, gR, zg, 3, dx);
        for (int k = 0; k < 3; ++k) in[0].D[k] = (float)Dk[k];
        in[0].ct = 0.0f;
        std::vector<float> t; std::vector<TpFaceOut> out;
        G.faces(in, t, t, false, out);
        RefFace r; ref_face(in[0], r);   // 同じ float 入力の double
        // 厳密入力 (double) の物理許容
        const double rgf = rho*(1 - 0.5*(gL + gR)), lim = 1e-6*rgf*3e-5*0.1;
        double zl[3], zr[3], s0 = 0, S = 0, jex[3];
        for (int k = 0; k < 3; ++k) { zl[k] = zg[k]; zr[k] = zg[k]; s0 += zl[k]; }
        for (int k = 0; k < 3; ++k) { jex[k] = -rgf*Dk[k]*(zr[k] - zl[k])/dx; S += jex[k]; }
        double jexm = 0; for (int k = 0; k < 3; ++k) jexm = std::max(jexm, fabs(jex[k] - zl[k]*S));
        // 演算誤差の許容 (Python ハーネス s1_tolerances の "a")
        const double a = 1.0/dx, ca = 3.0 + 3;
        double A[3], As = 0; for (int k = 0; k < 3; ++k) { A[k] = rgf*Dk[k]*a*EPS32*(r.zL[k]*ca + r.zR[k]*ca); As += A[k]; }
        double worst = 0, jmax = 0;
        for (int k = 0; k < 3; ++k) {
            const double tol = A[k] + 0.5*(r.zL[k] + r.zR[k])*As;
            const double gpu = (k == 2) ? (double)out[0].Jv : (double)out[0].J[k];
            const double ref = (k == 2) ? r.Jv : r.J[k];
            worst = std::max(worst, fabs(gpu - ref)/tol); jmax = std::max(jmax, fabs(gpu));
        }
        char b[384];
        snprintf(b, sizeof b, "T2 (i) dx=%g m: 厳密入力 double max|j| %.3e ≤ 物理許容 %.3e、GPU 演算誤差/許容 (ハーネス s1 \"a\") 最大 %.3f (GPU max|j| %.3e)",
                 dx, jexm, lim, worst, jmax);
        verdict(jexm <= lim && worst <= 1.0, b);
    }
}

static void test_T2_identities(int nF)
{
    printf("\n=== T2 (ii) 面恒等式 (GPU の float 出力から double で評価) ===\n");
    std::mt19937_64 rng(20261003);
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    for (int nG : {2, 3}) {
        std::vector<TpFaceIn> in(nF);
        for (int i = 0; i < nF; ++i) {
            TpFaceIn& a = in[i]; memset(&a, 0, sizeof a);
            a.n = nG; a.iw = nG - 1; a.f = (float)(0.3 + 0.4*U01(rng)); a.geo = (float)pow(10.0, -1 + 4*U01(rng)); a.geo_abs = a.geo;
            const double vch[4] = {0.0, 1e-8, 1e-4, 1e-2};
            for (int side = 0; side < 2; ++side) {
                const double rho = 0.05 + 1.95*U01(rng);
                const double g = (U01(rng) < 0.8) ? 0.3*U01(rng) : 0.0;
                const double vfrac = vch[(int)(4*U01(rng)) & 3]*U01(rng);
                double z[3]; double zs = 0; for (int k = 0; k < nG - 1; ++k) { z[k] = -log(std::max(U01(rng), 1e-300)); zs += z[k]; }
                for (int k = 0; k < nG - 1; ++k) z[k] = z[k]/zs*(1 - vfrac);
                float* rY = side ? a.rY1 : a.rY0;
                for (int k = 0; k < nG - 1; ++k) rY[k] = (float)((1 - g)*z[k]*rho);
                rY[nG - 1] = (float)(((1 - g)*vfrac + g)*rho);
                (side ? a.rho1 : a.rho0) = (float)rho; (side ? a.rg1 : a.rg0) = (float)(g*rho);
            }
            for (int k = 0; k < nG; ++k) a.D[k] = (float)(1e-5*(1 + 4*U01(rng)));
            a.ct = (float)(1e-5*U01(rng));
        }
        std::vector<float> t; std::vector<TpFaceOut> out;
        G.faces(in, t, t, false, out);
        double w1 = 0, w2 = 0;
        for (int i = 0; i < nF; ++i) {
            // Σ_気相 j_k (分子のみ) は乱流を除いて見る: j_k = J_k − J_t(Y_k)、j_v = J_v − J_t(v) を double で戻す
            const TpFaceIn& a = in[i]; const TpFaceOut& o = out[i];
            RefFace r; ref_face(a, r);
            double sj = 0, sa = 0;
            for (int k = 0; k < nG; ++k) {
                const double gpuJ = (k == a.iw) ? (double)o.Jv : (double)o.J[k];
                const double phiL = (k == a.iw) ? ((double)a.rY0[k] - a.rg0)/a.rho0 : (double)a.rY0[k]/a.rho0;
                const double phiR = (k == a.iw) ? ((double)a.rY1[k] - a.rg1)/a.rho1 : (double)a.rY1[k]/a.rho1;
                sj += gpuJ - (double)a.ct*a.geo*(phiR - phiL);
                sa += fabs(r.j0[k]);
            }
            // 乱流分の打ち消し誤差を許容に含める (被演算子の大きさ)
            double ta = 0; for (int k = 0; k < nG; ++k) ta += (double)a.ct*a.geo*(fabs(a.rY0[k]/a.rho0) + fabs(a.rY1[k]/a.rho1) + fabs(a.rg0/a.rho0) + fabs(a.rg1/a.rho1));
            w1 = std::max(w1, fabs(sj)/(8*EPS32*(sa + ta) + 1e-300));
            const double Jw = o.J[a.iw], Jv = o.Jv, Jl = o.Jl;
            w2 = std::max(w2, fabs(Jw - Jv - Jl)/(8*EPS32*(fabs(Jw) + fabs(Jv) + fabs(Jl)) + 1e-300));
        }
        // #4f (3): 分子の恒等式は同じ状態で ct = 0 にして直接検査する (許容は §6 の 8ε₃₂ Σ|j_k⁰|; 乱流の差し引きを混ぜない)。
        //   上の w1 は乱流を double で差し引き、分母を 8ε(Σ|j⁰| + 乱流の被演算子) に広げていた (#4e の記録; 判定に使わない)。
        std::vector<TpFaceIn> in0 = in; for (auto& a : in0) a.ct = 0.0f;
        std::vector<TpFaceOut> out0; G.faces(in0, t, t, false, out0);
        double w3 = 0; long nz = 0;
        for (int i = 0; i < nF; ++i) {
            RefFace r; ref_face(in0[i], r);
            double sj = 0, sa = 0;
            for (int k = 0; k < nG; ++k) { sj += (k == in0[i].iw) ? (double)out0[i].Jv : (double)out0[i].J[k]; sa += fabs(r.j0[k]); }
            if (sa == 0.0) { if (sj != 0.0) w3 = INFINITY; ++nz; continue; }
            w3 = std::max(w3, fabs(sj)/(8*EPS32*sa));
        }
        printf("[INFO] T2 (ii) #4e の記録 (乱流を差し引き・分母を広げた形): 気相 %d 種 |Σj_k|/(8ε(Σ|j⁰|+乱流)) 最大 %.3f\n", nG, w1);
        char b[320];
        snprintf(b, sizeof b, "T2 (ii) 気相 %d 種 × %d 面: ct = 0 で |Σ_気相 j_k|/(8ε₃₂ Σ|j_k⁰|) 最大 %.3f (Σ|j⁰| = 0 の面 %ld)、|J_w − J_v − J_l|/(8ε Σ|J|) 最大 %.3f", nG, nF, w3, nz, w2);
        verdict(w3 <= 1.0 && w2 <= 1.0, b);
    }
}

// 閉じた直列セルの BE 1 物理更新を GPU の面流束と GPU の更新で解き切る (乱流のみ、分子 D 指定、ρ = V = Δt = 1、V/Δτ = 1)。
//   非水種 (1 本) と BE 項はホストで同じ点対角の形 (BE の V/Δt を対角に)。停止: 全成分 (N, w, v, l) で max|r| ≤ max(1e-7 r0, 6ε max|Mρφ|)。
struct Box { std::vector<float> rN, rW, rG; int N; double ct, D; };
static double g_beMp = 1.0, g_beFloor = 6.0;   // 擬似時間の V/Δτ と停止の床 (ε の倍数)。判定は 1.0 / 6.0、他は診断 (判定外)
static int be_solve(Box& b, const Box& bn, int cap, double omega, double* minrv, double* corr, bool* capHit)
{
    const int N = b.N; double r0[4] = {0, 0, 0, 0};
    *capHit = false;
    for (int it = 0; it <= cap; ++it) {
        std::vector<TpFaceIn> in(N - 1);
        for (int i = 0; i + 1 < N; ++i) {
            TpFaceIn& a = in[i]; memset(&a, 0, sizeof a);
            a.n = 2; a.iw = 1; a.f = 0.5f; a.geo = 1.0f; a.geo_abs = 1.0f; a.rho0 = a.rho1 = 1.0f;
            a.rY0[0] = b.rN[i]; a.rY0[1] = b.rW[i]; a.rY1[0] = b.rN[i + 1]; a.rY1[1] = b.rW[i + 1];
            a.rg0 = b.rG[i]; a.rg1 = b.rG[i + 1]; a.D[0] = a.D[1] = (float)b.D; a.ct = (float)b.ct;
        }
        std::vector<float> t; std::vector<TpFaceOut> out;
        G.faces(in, t, t, false, out);
        std::vector<float> RN(N, 0.f), RW(N, 0.f), RG(N, 0.f), dN(N, 0.f), dW(N, 0.f), dG(N, 0.f);
        for (int i = 0; i + 1 < N; ++i) {
            RN[i] += out[i].J[0]; RN[i + 1] -= out[i].J[0]; RW[i] += out[i].J[1]; RW[i + 1] -= out[i].J[1]; RG[i] += out[i].Jl; RG[i + 1] -= out[i].Jl;
            dN[i] += out[i].diag0[0]; dN[i + 1] += out[i].diag1[0]; dW[i] += out[i].diag0[1]; dW[i + 1] += out[i].diag1[1];
            dG[i] += out[i].diagt0; dG[i + 1] += out[i].diagt1;
        }
        // BE の物理時間項 −V(ρφ − ρφⁿ)/Δt
        double rmax[4] = {0, 0, 0, 0}, Mmax[4] = {0, 0, 0, 0};
        for (int i = 0; i < N; ++i) {
            RN[i] -= (b.rN[i] - bn.rN[i]); RW[i] -= (b.rW[i] - bn.rW[i]); RG[i] -= (b.rG[i] - bn.rG[i]);
            const double rv = (double)RW[i] - (double)RG[i];
            rmax[0] = std::max(rmax[0], fabs((double)RN[i])); rmax[1] = std::max(rmax[1], fabs((double)RW[i]));
            rmax[2] = std::max(rmax[2], fabs(rv)); rmax[3] = std::max(rmax[3], fabs((double)RG[i]));
            Mmax[0] = std::max(Mmax[0], fabs((double)b.rN[i])); Mmax[1] = std::max(Mmax[1], fabs((double)b.rW[i]));
            Mmax[2] = std::max(Mmax[2], fabs((double)b.rW[i] - b.rG[i])); Mmax[3] = std::max(Mmax[3], fabs((double)b.rG[i]));
        }
        bool pass = true;
        for (int q = 0; q < 4; ++q) {
            if (it == 0) r0[q] = rmax[q];
            // 停止の床は #4a で固定した 6ε·max|M ρφ| (M = V/Δt = 1; 擬似時間の V/Δτ は前処理にだけ入る)。
            // 2026-10-02 初回は M = V/Δτ + V/Δt = 2 (12ε 相当) で回して 100 更新 6.72e-6 の FAIL (設計メモ §14.3)。
            const double tol = std::max(1e-7*r0[q], g_beFloor*EPS32*1.0*Mmax[q]);
            if (!(rmax[q] <= tol)) pass = false;
        }
        if (pass && it > 0) return it;
        if (it == cap) { *capHit = true; return it; }
        std::vector<TpCellIn> ci(N); std::vector<TpCellOut> co;
        for (int i = 0; i < N; ++i) {
            TpCellIn& c = ci[i]; memset(&c, 0, sizeof c);
            c.M = (float)(g_beMp + 1.0); c.V = 1.0f; c.Rw = RW[i]; c.Rg = RG[i]; c.Dv = dW[i]; c.Dg = dG[i];
            c.DQ[0] = c.DQ[1] = c.DQ[2] = 1.0f;
            c.rYw = b.rW[i]; c.rg = b.rG[i]; c.rho = 1.0f; c.omega = (float)omega; c.dg_max = 1e30; c.dT_max = 1e30; c.L = 0; c.cveff = 1;
        }
        G.cells(ci, co);
        for (int i = 0; i < N; ++i) {
            b.rN[i] += (float)omega*(RN[i]/((float)(g_beMp + 1.0) + dN[i]));
            b.rG[i] = co[i].rg; b.rW[i] = co[i].rYw;
            *corr += co[i].qcut + co[i].vround;
            *minrv = std::min(*minrv, (double)b.rW[i] - (double)b.rG[i]);
        }
    }
    return cap;
}

static void test_T2_nonneg_and_conservation(bool quick)
{
    printf("\n=== T2 非負 (2 セル) と 3 セル判別 B (閉じた箱・乱流のみ・BE を GPU の面流束と更新で解き切る) ===\n");
    {
        Box b; b.N = 2; b.ct = 1.0; b.D = 0.0;
        b.rG = {0.1f, 0.2f}; b.rW = {0.1f, 0.2f}; b.rN = {0.9f, 0.8f};
        Box bn = b; double minrv = INFINITY, corr = 0; bool capHit;
        const int it = be_solve(b, bn, 5000, 0.5, &minrv, &corr, &capHit);
        char m[256]; snprintf(m, sizeof m, "T2 非負 2 セル (蒸気 0, Y_w = g = 0.1/0.2): 反復 %d (上限到達 %s)、全反復の min ρv %.3e、状態補正 %.2e、解 ρg = [%.6f, %.6f]",
                              it, capHit ? "あり" : "なし", minrv, corr, b.rG[0], b.rG[1]);
        verdict(!capHit && minrv >= 0.0 && corr == 0.0, m);
    }
    for (int nupd : {1, quick ? 100 : 1000}) {
        Box b; b.N = 3; b.ct = 1.0; b.D = 0.0;
        b.rG = {0.1f, 0.2f, 0.2f}; b.rW = {0.2f, 0.3f, 0.3f}; b.rN = {0.8f, 0.7f, 0.7f};
        double t0[3] = {0, 0, 0}; for (int i = 0; i < 3; ++i) { t0[0] += b.rG[i]; t0[1] += b.rW[i]; t0[2] += b.rN[i]; }
        double minrv = INFINITY, corr = 0; int caps = 0, maxit = 0;
        for (int u = 0; u < nupd; ++u) {
            Box bn = b; bool capHit;
            const int it = be_solve(b, bn, 5000, 1.0, &minrv, &corr, &capHit);
            caps += capHit ? 1 : 0; maxit = std::max(maxit, it);
        }
        double t1[3] = {0, 0, 0}; for (int i = 0; i < 3; ++i) { t1[0] += b.rG[i]; t1[1] += b.rW[i]; t1[2] += b.rN[i]; }
        double e = 0; for (int k = 0; k < 3; ++k) e = std::max(e, fabs(t1[k] - t0[k])/t0[k]);
        double sy = 0; for (int i = 0; i < 3; ++i) sy = std::max(sy, fabs((double)b.rW[i] + b.rN[i] - 1.0));
        char m[320]; snprintf(m, sizeof m, "T2 3 セル判別 B %d 更新: 総液量・各種総量の相対変化 最大 %.2e (液 %.2e)、min ρv %.3e、上限到達 %d、最大反復 %d、状態補正 %.2e、|ΣρY−ρ| %.1e",
                              nupd, e, fabs(t1[0] - t0[0])/t0[0], minrv, caps, maxit, corr, sy);
        // 1000 更新は #4e の記録 (試験側 BE が擬似 V/Δτ=1・float 停止のとき FAIL 3.14e-6)。判定は #4f (2) の独立残差版 (test_4f_box) に移した。
        if (nupd == 1) verdict(e <= 1e-6 && minrv >= 0.0 && caps == 0, m);
        else printf("[INFO] #4e の記録 (旧判定 FAIL, 判定は #4f (2) へ): %s\n", m);
    }
    // 診断 (判定外): 1000 更新の累積が停止の床と前処理 (擬似時間の V/Δτ) のどちらで決まるか
    for (const double* v : (const double[][2]){{1.0, 2.0}, {0.0, 6.0}}) {
        g_beMp = v[0]; g_beFloor = v[1];
        Box b; b.N = 3; b.ct = 1.0; b.D = 0.0;
        b.rG = {0.1f, 0.2f, 0.2f}; b.rW = {0.2f, 0.3f, 0.3f}; b.rN = {0.8f, 0.7f, 0.7f};
        double t0 = 0; for (int i = 0; i < 3; ++i) t0 += b.rG[i];
        double minrv = INFINITY, corr = 0; int caps = 0, maxit = 0;
        for (int u = 0; u < (quick ? 100 : 1000); ++u) { Box bn = b; bool ch; const int it = be_solve(b, bn, 5000, 1.0, &minrv, &corr, &ch); caps += ch; maxit = std::max(maxit, it); }
        double t1 = 0; for (int i = 0; i < 3; ++i) t1 += b.rG[i];
        printf("[INFO] 診断 3 セル %d 更新: 擬似 V/Δτ %.0f・床 %.1fε → 総液量 %.2e、上限到達 %d、最大反復 %d\n", quick ? 100 : 1000, v[0], v[1], (t1 - t0)/t0, caps, maxit);
    }
    g_beMp = 1.0; g_beFloor = 6.0;
}

static void test_T2_energy(const Ctx& c)
{
    printf("\n=== T2 エネルギー: 分子 0・総組成一様・乱流のみ・g のみ勾配 — 1 更新の温度変化の EOS 接線 ===\n");
    const int N = 6; const double rho = RHO, Yw = YW_IN, T = 230.0;
    std::vector<TpFaceIn> in(N - 1); std::vector<float> t0(N - 1, (float)T), t1(N - 1, (float)T);
    // g は 2 次の分布 (1 次だと内部セルの r_g が丸めだけになり相対差が定義できない; 2026-10-02 初回は 1 次で FAIL 2.8e3、設計メモ §14.3)
    std::vector<double> g(N); for (int i = 0; i < N; ++i) g[i] = 2e-3*(double)(i*i)/((N - 1)*(N - 1));
    for (int i = 0; i + 1 < N; ++i) {
        TpFaceIn& a = in[i]; memset(&a, 0, sizeof a);
        a.n = 2; a.iw = 1; a.f = 0.5f; a.geo = (float)(1.0/DX); a.geo_abs = a.geo; a.rho0 = a.rho1 = (float)rho;
        a.rY0[0] = a.rY1[0] = (float)(rho*(1 - Yw)); a.rY0[1] = a.rY1[1] = (float)(rho*Yw);
        a.rg0 = (float)(rho*g[i]); a.rg1 = (float)(rho*g[i + 1]); a.ct = (float)(MUT/SCT);
    }
    // 2 種 (N2, H2O) の h は Ctx の sp[0], sp[1] (k_face の setProps)
    std::vector<TpFaceOut> out; G.faces(in, t0, t1, true, out);
    std::vector<double> rE(N, 0.0), rg(N, 0.0);
    for (int i = 0; i + 1 < N; ++i) { rE[i] += out[i].q; rE[i + 1] -= out[i].q; rg[i] += out[i].Jl; rg[i + 1] -= out[i].Jl; }
    double worst = 0;
    for (int i = 0; i < N; ++i) {
        if (rg[i] == 0.0) continue;
        const double Yd[2] = {1 - Yw, Yw};
        double cpm, hm; thermo_cph_mix(c.sp, 2, Yd, T, &cpm, &hm);
        const double cv = cpm - thermo_R_mix(c.sp, 2, Yd);
        const double L = cond_latent(c.cp, T), dL = (cond_latent(c.cp, T + 0.1) - cond_latent(c.cp, T - 0.1))/0.2;
        const double cveff = cv + g[i]*(c.Rw - dL);
        const double dTt = (rE[i] - (c.Rw*T - L)*rg[i])/(rho*cveff);
        const double dTa = -c.Rw*T*rg[i]/(rho*cveff);
        worst = std::max(worst, fabs(dTt - dTa)/fabs(dTa));
    }
    char m[160]; snprintf(m, sizeof m, "T2 エネルギー: 接線 ΔT と解析値 −R_w T r_g Δt/(ρ c_v,eff) の相対差 最大 %.2e (≤1 %%)", worst);
    verdict(worst <= 1e-2, m);
}

// ------------------------------------------------------------------ T3 1D 問題 (#4c/#4d) を GPU の面流束・更新で回す (受入 2)
// 移流・ソース・独立残差 (停止)・N2/E の更新はハーネスのまま。拡散 (面流束と対角) と蒸気・液・Q の更新だけ GPU。
// #4h: 非負制限を外した更新の後段 (本番の受動種の床 → 実現可能性クランプの下限/上限) を 1D で写す (float; 射影・消滅はこの 1D に無い)。
//   段ごとの補正 |Δ|V (成分 w, v, g, Q2, Q1, Q0) を更新開始時の格納値の総量で割って、更新ごとの C を返す (本番 [twophase-corr-gate] と同じ定義)。
struct TcUpd { double C[6]; };
static std::vector<TcUpd> g_t3_tc;
static RunOut run_gpu(const Ctx& c, double omega, int cap, int sweeps = 0, int noNonneg = 0)
{
    RunOut ro;
    State<float> st, inlet; for (auto& a : inlet.v) a.assign(1, 0.0f);
    const double rYw = RHO*YW_IN, rg = RHO*G_IN, rN2 = RHO - rYw;
    const double rE = e_of_T(c, rN2, rYw, rg, T_IN);
    const double init[NVAR] = {rN2, rYw, rg, RHO*Q0_IN, RHO*Q1_IN, RHO*Q2_IN, rE};
    for (int q = 0; q < NVAR; ++q) { inlet.v[q][0] = (float)init[q]; for (int i = 0; i < NC; ++i) st.v[q][i] = (float)init[q]; }
    std::vector<double> T(NC, T_IN);
    const double M = DX/(CFL*DX/U);
    double r0[NVAR + 1] = {};
    const int qlist[] = {K_N2, K_W, K_V, K_L, K_Q0, K_Q1, K_Q2, K_E};
    auto cell = [&](int i, int q) -> float { return (i < 0) ? inlet.v[q][0] : st.v[q][i]; };
    for (int it = 0; it <= cap; ++it) {
        for (int i = 0; i < NC; ++i) {
            bool ok; T[i] = eos_T(c, (double)st.v[K_N2][i], (double)st.v[K_W][i], (double)st.v[K_L][i], (double)st.v[K_E][i], T[i], &ok);
            if (!ok || !std::isfinite(T[i])) { ro.bad = "EOS"; ro.iters = it; return ro; }
        }
        State<double> s64; for (int q = 0; q < NVAR; ++q) for (int i = 0; i < NC; ++i) s64.v[q][i] = (double)st.v[q][i];
        State<double> in64; for (auto& v_ : in64.v) v_.assign(1, 0.0); for (int q = 0; q < NVAR; ++q) in64.v[q][0] = (double)inlet.v[q][0];
        Assembled a64; assemble<double>(c, s64, T, in64, a64, false);
        if (!scan_assembled(a64, "独立残差", it, ro)) { ro.iters = it; return ro; }
        bool pass = true; double rmax = 0.0;
        for (int q : qlist) {
            double m = 0.0, Am = 0.0; for (int i = 0; i < NC; ++i) { m = std::max(m, fabs(a64.r[q][i])); Am = std::max(Am, a64.A[q][i]); }
            if (it == 0) r0[q] = m;
            const double tol = std::max(1e-7*r0[q], 6.0*EPS32*Am);
            const double ratio = (tol > 0) ? m/tol : (m == 0 ? 0.0 : INFINITY);
            ro.comp_ratio[q] = ratio; rmax = std::max(rmax, ratio);
            if (!(ratio <= 1.0)) pass = false;
        }
        ro.ratio_max = rmax;
        if (pass && it > 0) { ro.converged = true; ro.iters = it; break; }
        if (it == cap) { ro.iters = it; break; }
        // 反復の残差 (float): 移流 (ハーネスと同じ) + 拡散 (GPU) + ソース (ハーネスの写し)
        std::vector<float> R[NVAR + 1]; std::vector<double> diag[NVAR + 1];
        for (int q = 0; q <= NVAR; ++q) { R[q].assign(NC, 0.0f); diag[q].assign(NC, 0.0); }
        const float rho = (float)RHO, mdot = (float)(RHO*U);
        std::vector<TpFaceIn> fin(NC); std::vector<float> t0(NC), t1(NC);
        for (int i = 0; i < NC; ++i) {   // 拡散のある面 i = 0..NC-1 (左 = i−1、i = 0 は入口ゴースト)
            TpFaceIn& a = fin[i]; memset(&a, 0, sizeof a);
            const int L = i - 1, Rr = i;
            a.n = 2; a.iw = 1; a.f = 0.5f; a.geo = (float)(1.0/DX); a.geo_abs = a.geo; a.rho0 = a.rho1 = rho;
            a.rY0[0] = cell(L, K_N2); a.rY0[1] = cell(L, K_W); a.rY1[0] = cell(Rr, K_N2); a.rY1[1] = cell(Rr, K_W);
            a.rg0 = cell(L, K_L); a.rg1 = cell(Rr, K_L);
            a.rQ0[0] = cell(L, K_Q2); a.rQ0[1] = cell(L, K_Q1); a.rQ0[2] = cell(L, K_Q0);
            a.rQ1[0] = cell(Rr, K_Q2); a.rQ1[1] = cell(Rr, K_Q1); a.rQ1[2] = cell(Rr, K_Q0);
            a.D[0] = a.D[1] = (float)DMOL; a.ct = (float)(MUT/SCT);
            t0[i] = (float)((L < 0) ? T_IN : T[L]); t1[i] = (float)T[Rr];
        }
        std::vector<TpFaceOut> fo; G.faces(fin, t0, t1, true, fo);
        for (int i = 0; i <= NC; ++i) {
            const int L = i - 1, Rr = i;
            float flux[NVAR + 1] = {};
            flux[K_N2] = mdot*(cell(L, K_N2)/rho); flux[K_L] = mdot*(cell(L, K_L)/rho);
            const float yv = (cell(L, K_W) - cell(L, K_L))/rho;
            flux[K_V] = mdot*yv; flux[K_W] = flux[K_V] + flux[K_L];
            flux[K_Q0] = mdot*(cell(L, K_Q0)/rho); flux[K_Q1] = mdot*(cell(L, K_Q1)/rho); flux[K_Q2] = mdot*(cell(L, K_Q2)/rho);
            flux[K_E] = mdot*(cell(L, K_E)/rho);
            if (Rr < NC) {   // GPU の面流束 (セル 0 = 左へ入る向き) を物理の向き (左→右) に直して足す
                const TpFaceOut& o = fo[i];
                flux[K_N2] += -o.J[0]; flux[K_W] += -o.J[1]; flux[K_L] += -o.Jl;
                flux[K_Q2] += -o.JQ[0]; flux[K_Q1] += -o.JQ[1]; flux[K_Q0] += -o.JQ[2]; flux[K_E] += -o.q;
                if (L >= 0) { diag[K_N2][L] += o.diag0[0]; diag[K_V][L] += o.diag0[1]; diag[K_L][L] += o.diagt0; for (int q : {K_Q0, K_Q1, K_Q2}) diag[q][L] += o.diagt0; }
                diag[K_N2][Rr] += o.diag1[0]; diag[K_V][Rr] += o.diag1[1]; diag[K_L][Rr] += o.diagt1; for (int q : {K_Q0, K_Q1, K_Q2}) diag[q][Rr] += o.diagt1;
            }
            if (L >= 0) for (int q = 0; q <= NVAR; ++q) diag[q][L] += RHO*U/RHO;
            for (int q = 0; q < NVAR; ++q) {
                if (L >= 0) R[q][L] -= flux[q];
                if (Rr < NC) R[q][Rr] += flux[q];
            }
        }
        std::vector<Src> src(NC);
        for (int i = 0; i < NC; ++i) {
            src[i] = cond_source_cell_f(c, (float)RHO, (float)T[i], st.v[K_W][i], st.v[K_L][i], st.v[K_Q0][i], st.v[K_Q1][i], st.v[K_Q2][i], false);
            if (src[i].oot) { ro.oot++; ro.bad = "表範囲外"; ro.iters = it; return ro; }
            const float V = (float)DX;
            R[K_L][i] += src[i].Sg*V; R[K_Q0][i] += src[i].SQ0*V; R[K_Q1][i] += src[i].SQ1*V; R[K_Q2][i] += src[i].SQ2*V;
        }
        // 蒸気・液・Q の更新 (GPU tp_vl_update)
        std::vector<TpCellIn> ci(NC); std::vector<TpCellOut> co;
        for (int i = 0; i < NC; ++i) {
            TpCellIn& cc = ci[i]; memset(&cc, 0, sizeof cc);
            cc.M = (float)M; cc.V = (float)DX;
            cc.Rw = R[K_W][i]; cc.Rg = R[K_L][i]; cc.RQ[0] = R[K_Q2][i]; cc.RQ[1] = R[K_Q1][i]; cc.RQ[2] = R[K_Q0][i];
            cc.Dv = (float)diag[K_V][i]; cc.Dg = (float)diag[K_L][i]; cc.DQ[0] = (float)diag[K_Q2][i]; cc.DQ[1] = (float)diag[K_Q1][i]; cc.DQ[2] = (float)diag[K_Q0][i];
            cc.sjg = src[i].sjg; cc.sjQ[0] = 0.0f; cc.sjQ[1] = src[i].sjq1; cc.sjQ[2] = 0.0f;
            cc.rYw = st.v[K_W][i]; cc.rg = st.v[K_L][i]; cc.rQ[0] = st.v[K_Q2][i]; cc.rQ[1] = st.v[K_Q1][i]; cc.rQ[2] = st.v[K_Q0][i];
            cc.rho = (float)RHO; cc.omega = (float)omega; cc.dg_max = DG_MAX; cc.dT_max = DT_MAX;
            double cpm, hm; const double Yd[2] = {(double)st.v[K_N2][i]/RHO, (double)st.v[K_W][i]/RHO};
            thermo_cph_mix(c.sp, 2, Yd, T[i], &cpm, &hm);
            const double gg = (double)st.v[K_L][i]/RHO;
            cc.L = cond_latent(c.cp, T[i]);
            cc.cveff = cpm - thermo_R_mix(c.sp, 2, Yd) + gg*(c.Rw - (cond_latent(c.cp, T[i] + 0.1) - cond_latent(c.cp, T[i] - 0.1))/0.2);
        }
        if (sweeps > 0) {   // #4g: 緩和整合 scalar-DPLUR の増分 (本番と同じ tp_denoms / tp_dplur_solve; 非対角は 1 次風上の流入 ṁ/ρ_左、ω = implicitRelax = 1)
            std::vector<float> R5(5*NC), D5(5*NC), a(5*NC, 0.0f), b(5*NC, 0.0f);
            for (int i = 0; i < NC; ++i) {
                float Dd[5]; tp_denoms(ci[i], Dd);
                const float Rr[5] = {ci[i].Rw - ci[i].Rg, ci[i].Rg, ci[i].RQ[0], ci[i].RQ[1], ci[i].RQ[2]};
                for (int q = 0; q < 5; ++q) { R5[q*NC + i] = Rr[q]; D5[q*NC + i] = Dd[q]; }
            }
            const float off = mdot/rho;
            for (int k = 0; k < sweeps; ++k) {
                for (int i = 0; i < NC; ++i) for (int q = 0; q < 5; ++q) {
                    const float nb = (i > 0) ? off*a[q*NC + i - 1] : 0.0f;
                    b[q*NC + i] = tp_dplur_solve(1.0f, R5[q*NC + i], nb, D5[q*NC + i]);
                }
                std::swap(a, b);
            }
            for (int i = 0; i < NC; ++i) { ci[i].useInc = 1; for (int q = 0; q < 5; ++q) ci[i].inc[q] = a[q*NC + i]; }
        }
        for (int i = 0; i < NC; ++i) ci[i].noNonneg = noNonneg;
        G.cells(ci, co);
        if (noNonneg) {   // 後段: 液と Q の床 (受動種の床) → 液の上限 ρg ≤ ρY_w (実現可能性クランプ; 総水分は固定)。補正を段ごとに計上
            double num[6] = {0, 0, 0, 0, 0, 0}, den[6] = {0, 0, 0, 0, 0, 0};
            for (int i = 0; i < NC; ++i) {
                TpCellOut& o = co[i];
                const double d0[6] = {ci[i].rYw, (double)ci[i].rYw - ci[i].rg, ci[i].rg, ci[i].rQ[0], ci[i].rQ[1], ci[i].rQ[2]};
                for (int k = 0; k < 6; ++k) den[k] += d0[k]*DX;
                num[0] += (o.vround + o.wfloor)*DX; num[1] += (o.vround + o.wfloor)*DX;
                for (int m = 0; m < 3; ++m) num[3 + m] += o.qc[m]*DX;
                float g = o.rg; if (g < 0.0f) { num[2] += -(double)g*DX; num[1] += -(double)g*DX; g = 0.0f; }
                for (int m = 0; m < 3; ++m) if (o.rQ[m] < 0.0f) { num[3 + m] += -(double)o.rQ[m]*DX; o.rQ[m] = 0.0f; }
                if (g > o.rYw) { num[2] += (double)(g - o.rYw)*DX; num[1] += (double)(g - o.rYw)*DX; g = o.rYw; }
                o.rg = g;
            }
            TcUpd t; for (int k = 0; k < 6; ++k) t.C[k] = (den[k] == 0.0) ? (num[k] == 0.0 ? 0.0 : INFINITY) : num[k]/den[k];
            g_t3_tc.push_back(t);
        }
        double qcut_this = 0, vround_this = 0; int thetaLastLt1 = 0;
        for (int i = 0; i < NC; ++i) {
            const TpCellOut& o = co[i];
            if (!std::isfinite(o.theta)) { ro.bad = "非有限 θ"; ro.iters = it; return ro; }
            ro.theta_min = std::min(ro.theta_min, o.theta); if (o.theta < 1.0) { ++ro.theta_lt1; ++thetaLastLt1; }
            ro.theta_src_min = std::min(ro.theta_src_min, (double)src[i].theta);
            ro.withheld_v += o.withheld_v; ro.withheld_g += o.withheld_g;
            ro.qcut += o.qcut; ro.vround += o.vround; qcut_this += o.qcut; vround_this += o.vround;
            float dN2 = R[K_N2][i]/((float)M + (float)diag[K_N2][i]);
            float dE  = R[K_E][i]/((float)M + (float)U);
            if (omega != 1.0) { dN2 *= (float)omega; dE *= (float)omega; }
            st.v[K_L][i] = o.rg; st.v[K_W][i] = o.rYw; st.v[K_Q2][i] = o.rQ[0]; st.v[K_Q1][i] = o.rQ[1]; st.v[K_Q0][i] = o.rQ[2];
            st.v[K_N2][i] += dN2; st.v[K_E][i] += dE;
            ro.min_rv = std::min(ro.min_rv, (double)(st.v[K_W][i] - st.v[K_L][i]));
            ro.min_rg = std::min(ro.min_rg, (double)st.v[K_L][i]);
            for (int q : {K_Q0, K_Q1, K_Q2}) ro.min_Q = std::min(ro.min_Q, (double)st.v[q][i]);
        }
        ro.theta_final_lt1 = thetaLastLt1;
        if (thetaLastLt1 > 0) { ++ro.theta_lt1_updates; ro.theta_last_it = it; }
        ro.qcut_it.push_back(qcut_this); ro.vround_it.push_back(vround_this);
    }
    {
        const int n = (int)ro.qcut_it.size();
        ro.tail_n = std::max(1, (int)std::ceil(0.1*n));
        for (int k = std::max(0, n - ro.tail_n); k < n; ++k) { ro.qcut_last10 += ro.qcut_it[k]; ro.vround_last10 += ro.vround_it[k]; }
    }
    {
        Assembled af; State<double> s64; for (int q = 0; q < NVAR; ++q) for (int i = 0; i < NC; ++i) s64.v[q][i] = (double)st.v[q][i];
        State<double> in64; for (auto& v_ : in64.v) v_.assign(1, 0.0); for (int q = 0; q < NVAR; ++q) in64.v[q][0] = (double)inlet.v[q][0];
        assemble<double>(c, s64, T, in64, af, false);
        for (int i = 0; i < NC; ++i) if (af.src[i].theta < 1.0f) ++ro.theta_src_final_lt1;
        ro.g_out = s64.v[K_L][NC - 1]/RHO; ro.T_out = T[NC - 1];
    }
    return ro;
}

static void test_T3(const Ctx& c, int cap)
{
    printf("\n=== T3 #4c/#4d の 1D 問題: 拡散と蒸気・液・Q の更新を GPU に置き換えてホスト版と比べる (float32, CFL 5, ω 1, 対角, 上限 %d) ===\n", cap);
    for (double dtm : {1.0, 0.01}) {
        DT_MAX = dtm;
        const RunOut h = run<float>(c, RunOpt{1.0, false, cap, "host"});
        const RunOut d = run_gpu(c, 1.0, cap);
        const bool okH = accept(h), okD = accept(d);
        const double dg = fabs(d.g_out - h.g_out)/h.g_out, dT = fabs(d.T_out - h.T_out);
        const double dit = fabs((double)d.iters - h.iters)/h.iters;
        const bool act = (dtm < 1.0) ? (d.theta_lt1 > 0) : true;
        char m[640];
        snprintf(m, sizeof m, "T3 DT_MAX %.2g K: GPU %s 反復 %d (最大比 %.3f [v %.2f], θ 最小 %.3f, θ<1 %ld, 最後の θ<1 %d, θ_src 最終<1 %d, 末尾補正 %.1e/%.1e, min ρv %.2e) / "
                 "ホスト %s 反復 %d、出口 g %.6e / %.6e (相対 %.1e)、出口 T %.4f / %.4f K (差 %.1e)、反復差 %.1f %%",
                 dtm, okD ? "合格" : "不合格", d.iters, d.ratio_max, d.comp_ratio[K_V], d.theta_min, d.theta_lt1, d.theta_final_lt1, d.theta_src_final_lt1,
                 d.qcut_last10, d.vround_last10, d.min_rv, okH ? "合格" : "不合格", h.iters, d.g_out, h.g_out, dg, d.T_out, h.T_out, dT, 100*dit);
        verdict(okD && okH && dg <= 1e-4 && dT <= 0.01 && dit <= 0.10 && act, m);
        if (!d.bad.empty()) printf("    GPU 打ち切り: %s\n", d.bad.c_str());
    }
    DT_MAX = 1.0;
}

// ------------------------------------------------------------------ #4f (1) θ 丸めの判別 A/B (plan §5.1 #4f の入力そのまま)
static void test_4f_theta()
{
    printf("\n=== #4f (1) θ の float 丸め A/B: 1 セル (M=V=ρ=ω=1、輸送対角・ソース Jacobian 0、rYw 0.01f・rg 3e-6f・Rw 0・Rg −1.45e-4f) ===\n");
    TpCellIn c; memset(&c, 0, sizeof c);
    c.M = 1.0f; c.V = 1.0f; c.rho = 1.0f; c.omega = 1.0f;
    c.rYw = 0.01f; c.rg = 3e-6f; c.Rw = 0.0f; c.Rg = -1.45e-4f;
    c.dg_max = 0.005; c.dT_max = 1.0; c.L = 2.5e6; c.cveff = 1000.0;
    TpCellOut* d_o; float* d_t; CK(cudaMalloc(&d_o, sizeof(TpCellOut))); CK(cudaMalloc(&d_t, sizeof(float)));
    TpCellOut o[2]; float Th[2];
    for (int m = 0; m < 2; ++m) {
        k_update_probe<<<1, 1>>>(c, m, d_o, d_t); CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
        CK(cudaMemcpy(&o[m], d_o, sizeof(TpCellOut), cudaMemcpyDeviceToHost)); CK(cudaMemcpy(&Th[m], d_t, sizeof(float), cudaMemcpyDeviceToHost));
        const double rv = (double)o[m].rYw - (double)o[m].rg;
        printf("[INFO] %s: θ (double) %.17g → float %.9g (%s)、液 ρg %.9g、蒸気 ρv %.9g、Q [%g %g %g]、qcut %.3g、vround %.3g\n",
               m ? "B (切り上がったら 0 側の隣)" : "A (最近接)", o[m].theta, (double)Th[m], ((double)Th[m] > o[m].theta) ? "切り上がり" : "切り下がりか一致",
               (double)o[m].rg, rv, o[m].rQ[0], o[m].rQ[1], o[m].rQ[2], o[m].qcut, o[m].vround);
    }
    auto nonneg = [](const TpCellOut& x) { return x.rg >= 0.0f && (double)x.rYw - (double)x.rg >= 0.0 && x.rQ[0] >= 0.0f && x.rQ[1] >= 0.0f && x.rQ[2] >= 0.0f; };
    const bool aNeg = !(o[0].rg >= 0.0f), bOk = nonneg(o[1]) && o[1].qcut == 0.0 && o[1].vround == 0.0;
    const char* rule;
    if (aNeg && bOk) rule = "A だけ液が負・B は全量非負・補正 0 → B を採用";
    else if (aNeg) rule = "A・B とも負 → commit を調べる";
    else rule = "A も非負 → GPU での反例不成立 (演算順との差を確認)";
    printf("[判定] #4f (1): %s\n", rule);
    verdict(nonneg(o[1]) && o[1].qcut == 0.0 && o[1].vround == 0.0, "#4f (1) B (安全側の丸め) で液・蒸気・Q が非負、qcut・vround 0");
    CK(cudaFree(d_o)); CK(cudaFree(d_t));
}

// ------------------------------------------------------------------ #4f (2) 3 セル 1000 更新: 試験側 BE を V/Δτ = 0 で、全成分を独立残差で
// 停止・判定とも、格納 float 状態を double に上げて ref_face (double) で組んだ BE 残差 (独立残差) で行う。床 6ε·max|(V/Δt)ρφ|、相対 1e-7·r0。
static void box_indep_resid(const Box& b, const Box& bn, double r[4][8], double sc[4])
{
    const int N = b.N;
    for (int q = 0; q < 4; ++q) { for (int i = 0; i < N; ++i) r[q][i] = 0.0; sc[q] = 0.0; }
    for (int i = 0; i + 1 < N; ++i) {
        TpFaceIn a; memset(&a, 0, sizeof a);
        a.n = 2; a.iw = 1; a.f = 0.5f; a.geo = 1.0f; a.geo_abs = 1.0f; a.rho0 = a.rho1 = 1.0f;
        a.rY0[0] = b.rN[i]; a.rY0[1] = b.rW[i]; a.rY1[0] = b.rN[i + 1]; a.rY1[1] = b.rW[i + 1];
        a.rg0 = b.rG[i]; a.rg1 = b.rG[i + 1]; a.D[0] = a.D[1] = (float)b.D; a.ct = (float)b.ct;
        RefFace f; ref_face(a, f);
        const double F[4] = {f.J[0], f.J[1], f.Jv, f.Jl};   // N, w, v, l (セル i へ入る向き)
        for (int q = 0; q < 4; ++q) { r[q][i] += F[q]; r[q][i + 1] -= F[q]; }
    }
    for (int i = 0; i < N; ++i) {
        const double x[4] = {b.rN[i], b.rW[i], (double)b.rW[i] - b.rG[i], b.rG[i]}, xn[4] = {bn.rN[i], bn.rW[i], (double)bn.rW[i] - bn.rG[i], bn.rG[i]};
        for (int q = 0; q < 4; ++q) { r[q][i] -= (x[q] - xn[q]); sc[q] = std::max(sc[q], fabs(x[q])); }
    }
}
static int be_solve_indep(Box& b, const Box& bn, int cap, double* minrv, double* corr, bool* capHit, double ratio_out[4])
{
    const int N = b.N; double r0[4] = {0, 0, 0, 0};
    *capHit = false;
    for (int it = 0; it <= cap; ++it) {
        double r[4][8], sc[4]; box_indep_resid(b, bn, r, sc);
        bool pass = true;
        for (int q = 0; q < 4; ++q) {
            double m = 0; for (int i = 0; i < N; ++i) m = std::max(m, fabs(r[q][i]));
            if (it == 0) r0[q] = m;
            const double tol = std::max(1e-7*r0[q], 6.0*EPS32*sc[q]);
            ratio_out[q] = (tol > 0) ? m/tol : (m == 0 ? 0.0 : INFINITY);
            if (!(ratio_out[q] <= 1.0)) pass = false;
        }
        if (pass && it > 0) return it;
        if (it == cap) { *capHit = true; return it; }
        // 反復 (前処理): GPU の面流束と GPU 更新、擬似時間 V/Δτ = 0 (対角は BE の V/Δt = 1 だけ)
        std::vector<TpFaceIn> in(N - 1);
        for (int i = 0; i + 1 < N; ++i) {
            TpFaceIn& a = in[i]; memset(&a, 0, sizeof a);
            a.n = 2; a.iw = 1; a.f = 0.5f; a.geo = 1.0f; a.geo_abs = 1.0f; a.rho0 = a.rho1 = 1.0f;
            a.rY0[0] = b.rN[i]; a.rY0[1] = b.rW[i]; a.rY1[0] = b.rN[i + 1]; a.rY1[1] = b.rW[i + 1];
            a.rg0 = b.rG[i]; a.rg1 = b.rG[i + 1]; a.D[0] = a.D[1] = (float)b.D; a.ct = (float)b.ct;
        }
        std::vector<float> t; std::vector<TpFaceOut> out; G.faces(in, t, t, false, out);
        std::vector<float> RN(N, 0.f), RW(N, 0.f), RG(N, 0.f), dN(N, 0.f), dW(N, 0.f), dG(N, 0.f);
        for (int i = 0; i + 1 < N; ++i) {
            RN[i] += out[i].J[0]; RN[i + 1] -= out[i].J[0]; RW[i] += out[i].J[1]; RW[i + 1] -= out[i].J[1]; RG[i] += out[i].Jl; RG[i + 1] -= out[i].Jl;
            dN[i] += out[i].diag0[0]; dN[i + 1] += out[i].diag1[0]; dW[i] += out[i].diag0[1]; dW[i + 1] += out[i].diag1[1];
            dG[i] += out[i].diagt0; dG[i + 1] += out[i].diagt1;
        }
        for (int i = 0; i < N; ++i) { RN[i] -= (b.rN[i] - bn.rN[i]); RW[i] -= (b.rW[i] - bn.rW[i]); RG[i] -= (b.rG[i] - bn.rG[i]); }
        std::vector<TpCellIn> ci(N); std::vector<TpCellOut> co;
        for (int i = 0; i < N; ++i) {
            TpCellIn& c = ci[i]; memset(&c, 0, sizeof c);
            c.M = 1.0f; c.V = 1.0f; c.Rw = RW[i]; c.Rg = RG[i]; c.Dv = dW[i]; c.Dg = dG[i];
            c.DQ[0] = c.DQ[1] = c.DQ[2] = 1.0f;
            c.rYw = b.rW[i]; c.rg = b.rG[i]; c.rho = 1.0f; c.omega = 1.0f; c.dg_max = 1e30; c.dT_max = 1e30; c.L = 0; c.cveff = 1;
        }
        G.cells(ci, co);
        for (int i = 0; i < N; ++i) {
            b.rN[i] += RN[i]/(1.0f + dN[i]);
            b.rG[i] = co[i].rg; b.rW[i] = co[i].rYw;
            *corr += co[i].qcut + co[i].vround;
            *minrv = std::min(*minrv, (double)b.rW[i] - (double)b.rG[i]);
        }
    }
    return cap;
}
static void test_4f_box(bool quick)
{
    const int nupd = quick ? 100 : 1000;
    printf("\n=== #4f (2) 3 セル判別 B %d 更新: 試験側 BE を V/Δτ = 0、停止・判定とも独立残差 (double, 全成分 N・w・v・l)、床 6ε ===\n", nupd);
    Box b; b.N = 3; b.ct = 1.0; b.D = 0.0;
    b.rG = {0.1f, 0.2f, 0.2f}; b.rW = {0.2f, 0.3f, 0.3f}; b.rN = {0.8f, 0.7f, 0.7f};
    double t0[4] = {0, 0, 0, 0}; for (int i = 0; i < 3; ++i) { t0[0] += b.rN[i]; t0[1] += b.rW[i]; t0[2] += (double)b.rW[i] - b.rG[i]; t0[3] += b.rG[i]; }
    double minrv = INFINITY, corr = 0, rmax[4] = {0, 0, 0, 0}; int caps = 0, maxit = 0;
    for (int u = 0; u < nupd; ++u) {
        Box bn = b; bool ch; double rr[4];
        const int it = be_solve_indep(b, bn, 5000, &minrv, &corr, &ch, rr);
        caps += ch; maxit = std::max(maxit, it);
        for (int q = 0; q < 4; ++q) rmax[q] = std::max(rmax[q], rr[q]);
    }
    double t1[4] = {0, 0, 0, 0}; for (int i = 0; i < 3; ++i) { t1[0] += b.rN[i]; t1[1] += b.rW[i]; t1[2] += (double)b.rW[i] - b.rG[i]; t1[3] += b.rG[i]; }
    double e[4], em = 0; for (int q = 0; q < 4; ++q) { e[q] = fabs(t1[q] - t0[q])/t0[q]; em = std::max(em, e[q]); }
    double rm = 0; for (int q = 0; q < 4; ++q) rm = std::max(rm, rmax[q]);
    char m[480];
    snprintf(m, sizeof m, "#4f (2) %d 更新: 総量の相対変化 N %.2e・w %.2e・v %.2e・l %.2e (≤1e-6)、各更新の終了時の独立残差比 最大 N %.2f・w %.2f・v %.2f・l %.2f、上限到達 %d、最大反復 %d、min ρv %.3e、状態補正 %.2e",
             nupd, e[0], e[1], e[2], e[3], rmax[0], rmax[1], rmax[2], rmax[3], caps, maxit, minrv, corr);
    verdict(em <= 1e-6 && rm <= 1.0 && caps == 0 && minrv >= 0.0, m);
}

// ------------------------------------------------------------------ #4g (1): R = 0 で全増分・補正 0、1 sweep・ω 1 で点対角とビット一致
__global__ void k_dplur1(int n, TpCellIn* in)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= n) return;
    TpCellIn c = in[i];
    float D[5]; tp_denoms(c, D);
    const float R[5] = {c.Rw - c.Rg, c.Rg, c.RQ[0], c.RQ[1], c.RQ[2]};
    for (int q = 0; q < 5; ++q) c.inc[q] = tp_dplur_solve(1.0f, R[q], 0.0f, D[q]);   // ゼロ開始の 1 sweep (近傍の δ = 0)
    c.useInc = 1;
    in[i] = c;
}
static bool same_out(const TpCellOut& a, const TpCellOut& b)
{
    return !memcmp(&a.rYw, &b.rYw, sizeof(float)) && !memcmp(&a.rg, &b.rg, sizeof(float)) && !memcmp(a.rQ, b.rQ, sizeof(a.rQ))
        && a.theta == b.theta && a.qcut == b.qcut && a.vround == b.vround && !memcmp(&a.dv, &b.dv, sizeof(float)) && !memcmp(&a.dg, &b.dg, sizeof(float))
        && !memcmp(a.dq, b.dq, sizeof(a.dq)) && a.reason == b.reason;
}
static void test_4g_unit()
{
    printf("\n=== #4g (1) DPLUR の単体: R = 0 で増分・補正 0、1 sweep・ω 1 で点対角とビット一致 (GPU, 乱数 2 万セル) ===\n");
    std::mt19937_64 rng(20261004); std::uniform_real_distribution<double> U(0.0, 1.0);
    const int n = 20000;
    std::vector<TpCellIn> base(n);
    for (int i = 0; i < n; ++i) {
        TpCellIn& c = base[i]; memset(&c, 0, sizeof c);
        c.M = (float)(1e-3 + U(rng)); c.V = (float)(1e-9 + 1e-6*U(rng)); c.rho = (float)(0.1 + U(rng));
        c.Rw = (float)((U(rng) - 0.5)*1e-6); c.Rg = (float)((U(rng) - 0.5)*1e-7); for (int m = 0; m < 3; ++m) c.RQ[m] = (float)((U(rng) - 0.5)*pow(10.0, 4*m));
        c.Dv = (float)U(rng); c.Dg = (float)U(rng); c.sjg = (float)(U(rng)*1e3); for (int m = 0; m < 3; ++m) { c.DQ[m] = (float)U(rng); c.sjQ[m] = (float)(U(rng)*1e3); }
        c.rg = (U(rng) < 0.3) ? 0.0f : (float)(1e-3*U(rng)); c.rYw = c.rg + ((U(rng) < 0.2) ? 0.0f : (float)(1e-2*U(rng)));
        for (int m = 0; m < 3; ++m) c.rQ[m] = (U(rng) < 0.3) ? 0.0f : (float)pow(10.0, 4*m)*(float)U(rng);
        c.omega = (U(rng) < 0.5) ? 1.0f : 0.5f; c.dg_max = 5e-3; c.dT_max = (U(rng) < 0.5) ? 1.0 : 0.01; c.L = 2.5e6; c.cveff = 1000.0;
    }
    // 非零残差: 点対角 vs DPLUR 1 sweep
    std::vector<TpCellOut> oP, oD;
    G.cells(base, oP);
    TpCellIn* d; CK(cudaMalloc(&d, n*sizeof(TpCellIn))); CK(cudaMemcpy(d, base.data(), n*sizeof(TpCellIn), cudaMemcpyHostToDevice));
    k_dplur1<<<(n + 127)/128, 128>>>(n, d); CK(cudaDeviceSynchronize());
    std::vector<TpCellIn> dp(n); CK(cudaMemcpy(dp.data(), d, n*sizeof(TpCellIn), cudaMemcpyDeviceToHost));
    G.cells(dp, oD);
    long diff = 0, thlt1 = 0; for (int i = 0; i < n; ++i) { if (!same_out(oP[i], oD[i])) ++diff; if (oP[i].theta < 1.0) ++thlt1; }
    char m[256]; snprintf(m, sizeof m, "#4g (1) 非零残差 %d セル (θ<1 %ld): 1 sweep・ω 1 の DPLUR と点対角の出力 (状態・θ・補正・増分・制限) が異なるセル %ld", n, thlt1, diff);
    verdict(diff == 0, m);
    // R = 0
    for (auto& c : base) { c.Rw = c.Rg = 0.0f; for (int m2 = 0; m2 < 3; ++m2) c.RQ[m2] = 0.0f; }
    CK(cudaMemcpy(d, base.data(), n*sizeof(TpCellIn), cudaMemcpyHostToDevice));
    k_dplur1<<<(n + 127)/128, 128>>>(n, d); CK(cudaDeviceSynchronize());
    CK(cudaMemcpy(dp.data(), d, n*sizeof(TpCellIn), cudaMemcpyDeviceToHost));
    G.cells(dp, oD);
    long bad = 0;
    for (int i = 0; i < n; ++i) {
        const TpCellOut& o = oD[i]; const TpCellIn& c = dp[i];
        bool z = (o.dv == 0.0f && o.dg == 0.0f && o.dq[0] == 0.0f && o.dq[1] == 0.0f && o.dq[2] == 0.0f && o.qcut == 0.0 && o.vround == 0.0 && o.theta == 1.0
                  && !memcmp(&o.rg, &c.rg, 4) && !memcmp(&o.rYw, &c.rYw, 4) && !memcmp(o.rQ, c.rQ, sizeof(o.rQ)));
        for (int q = 0; q < 5; ++q) z = z && (c.inc[q] == 0.0f);
        if (!z) ++bad;
    }
    snprintf(m, sizeof m, "#4g (1) R = 0 の %d セル: 増分・補正が 0 でない、または状態が変わったセル %ld", n, bad);
    verdict(bad == 0, m);
    cudaFree(d);
}

static void test_4g_T3(const Ctx& c, int cap)
{
    printf("\n=== #4g (3)(4) T3 の 1D 問題を DPLUR (5 sweep, ω = implicitRelax = 1, condTwoPhaseRelax 1) で (float32, CFL 5, 上限 %d) ===\n", cap);
    for (double dtm : {1.0, 0.01}) {
        DT_MAX = dtm;
        const RunOut h = run<float>(c, RunOpt{1.0, false, cap, "host 点対角"});
        const RunOut d = run_gpu(c, 1.0, cap, 5);
        long negs = (d.min_rv < 0) + (d.min_rg < 0) + (d.min_Q < 0);
        bool fin = d.bad.empty(); for (int q = 0; q <= NVAR; ++q) fin = fin && std::isfinite(d.comp_ratio[q]);
        const bool ok3 = d.converged && fin && d.ratio_max <= 1.0 && negs == 0 && d.qcut_last10 == 0 && d.vround_last10 == 0
                         && d.theta_final_lt1 == 0 && d.theta_src_final_lt1 == 0 && (dtm >= 1.0 || d.theta_lt1 > 0);
        const double dg = fabs(d.g_out - h.g_out)/h.g_out, dT = fabs(d.T_out - h.T_out);
        char m[640];
        snprintf(m, sizeof m, "#4g (3) DT_MAX %.2g K DPLUR: %s 反復 %d、独立残差の最大比 %.3f [w %.2f v %.2f l %.2f Q0 %.2f Q1 %.2f Q2 %.2f E %.2f]、min ρv %.2e ρg %.2e ρQ %.2e、"
                 "末尾 10%% の Qcut/vround %.1e/%.1e、最後の θ<1 %d・θ_src<1 %d、θ<1 のセル·反復 %ld (θ 最小 %.3f)%s",
                 dtm, d.converged ? "収束" : (d.bad.empty() ? "上限到達" : d.bad.c_str()), d.iters, d.ratio_max, d.comp_ratio[K_W], d.comp_ratio[K_V], d.comp_ratio[K_L],
                 d.comp_ratio[K_Q0], d.comp_ratio[K_Q1], d.comp_ratio[K_Q2], d.comp_ratio[K_E], d.min_rv, d.min_rg, d.min_Q, d.qcut_last10, d.vround_last10,
                 d.theta_final_lt1, d.theta_src_final_lt1, d.theta_lt1, d.theta_min, (dtm < 1.0 && d.theta_lt1 == 0) ? " — θ<1 が起きない (判別不成立)" : "");
        verdict(ok3, m);
        snprintf(m, sizeof m, "#4g (4) DT_MAX %.2g K: 出口 g %.6e / 点対角 %.6e (相対 %.1e ≤ 1e-4)、出口 T %.4f / %.4f K (差 %.1e ≤ 0.01)、反復 DPLUR %d / 点対角 %d (記録のみ)",
                 dtm, d.g_out, h.g_out, dg, d.T_out, h.T_out, dT, d.iters, h.iters);
        verdict(dg <= 1e-4 && dT <= 0.01, m);
    }
    DT_MAX = 1.0;
}

// ------------------------------------------------------------------ #4h 非負制限 θ_vg を外す (condTwoPhaseNonnegLimit 0)
static void test_4h(const Ctx& c, int cap)
{
    printf("\n=== #4h (a) codex の反例 (θ_thr = 1, ρY_w 0.01, ρg 0.0099, δρv −0.0002, δρg 0): 本番の更新 → 受動種の床 → 実現可能性クランプ (GPU) ===\n");
    TpCellIn in; memset(&in, 0, sizeof in);
    in.M = 1.0f; in.V = 1.0f; in.rho = 1.0f; in.omega = 1.0f; in.dg_max = 1e30; in.dT_max = 1e30; in.L = 2.5e6; in.cveff = 1000.0;
    in.rYw = 0.01f; in.rg = 0.0099f; in.useInc = 1; in.inc[0] = -0.0002f; in.inc[1] = 0.0f;
    for (int mode = 0; mode < 2; ++mode) {
        TpCellIn ci = in; ci.noNonneg = mode;   // 0: 従来 (キー 1)、1: 非負制限なし (キー 0)
        std::vector<TpCellIn> v1(1, ci); std::vector<TpCellOut> o1; G.cells(v1, o1);
        const TpCellOut o = o1[0];
        // 後段 (本番のカーネルをそのまま): 受動種の床 (液) → 実現可能性クランプ (evap 0, 射影 1)
        float h[7] = {1.0f /*ro*/, o.rYw, o.rg, o.rQ[2] /*Q0*/, o.rQ[1] /*Q1*/, o.rQ[0] /*Q2*/, 230.0f /*T*/};
        float* d; CK(cudaMalloc(&d, sizeof h)); CK(cudaMemcpy(d, h, sizeof h, cudaMemcpyHostToDevice));
        double* acc; CK(cudaMalloc(&acc, 36*sizeof(double))); CK(cudaMemset(acc, 0, 36*sizeof(double)));
        geom_float vol = 1.0; geom_float* dv; CK(cudaMalloc(&dv, sizeof vol)); CK(cudaMemcpy(dv, &vol, sizeof vol, cudaMemcpyHostToDevice));
        double* st; CK(cudaMalloc(&st, 8*sizeof(double))); CK(cudaMemset(st, 0, 8*sizeof(double)));
        passive_bounds_d<<<1, 1>>>(1, d + 2, 0, d, dv, nullptr, st, nullptr, acc, 2);
        CondTablesF tbz{};
        cond_realizability_clamp_f_d<<<1, 1>>>(1, d, d + 1, d + 2, d + 3, d + 4, d + 5, 0, (float)c.Rw, 1e-9f, 5e-7f, 0.0f, d + 6, d + 6, tbz, c.cp,
                                               nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 1, nullptr, acc);
        CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
        float r[7]; CK(cudaMemcpy(r, d, sizeof r, cudaMemcpyDeviceToHost));
        double ha[36]; CK(cudaMemcpy(ha, acc, sizeof ha, cudaMemcpyDeviceToHost));
        const double w0 = in.rYw, g0 = in.rg, wtrial = (double)in.rYw + (double)in.inc[0];
        const double cw = (o.vround + o.wfloor)/w0, cgFloor = ha[1*6 + 2]/g0, cgClamp = ha[2*6 + 2]/g0;
        char m[480];
        snprintf(m, sizeof m, "%s: θ %.6g、commit 後 ρY_w %.7g ρg %.7g (vround %.2e、総水分の下限 %.2e)、後段の後 ρY_w %.7g ρg %.7g ρv %.3e、"
                 "総水分 初期比 %+.4e (試行 %.7g)、補正 C: 総水分 commit %.3e、液 床 %.3e・上限 %.3e",
                 mode ? "キー 0 (非負制限なし)" : "キー 1 (従来)", o.theta, (double)o.rYw, (double)o.rg, o.vround, o.wfloor, (double)r[1], (double)r[2],
                 (double)r[1] - (double)r[2], ((double)r[1] - w0)/w0, wtrial, cw, cgFloor, cgClamp);
        // 判定: 総水分は commit の試行値 (float の格納値) から後段で変わらない (増やさない)・vround を使わない・0 ≤ ρg ≤ ρY_w
        if (mode == 1) verdict(r[1] == o.rYw && o.rYw <= in.rYw && o.vround == 0.0 && r[2] >= 0.0f && r[2] <= r[1],
                               (std::string("#4h (a) ") + m + " — 総水分を増やさない・0 ≤ ρg ≤ ρY_w").c_str());
        else printf("[INFO] #4h (a) %s\n", m);
        cudaFree(d); cudaFree(acc); cudaFree(dv); cudaFree(st);
    }
    printf("\n=== #4h (b) T3 の 1D 問題を DPLUR (5 sweep) + キー 0 で (float32, CFL 5, 上限 %d; 後段の床・上限を 1D に写して補正を計上) ===\n", cap);
    const double kappa = 4.0*EPS32;
    for (double dtm : {1.0, 0.01}) {
        DT_MAX = dtm;
        g_t3_tc.clear();
        const RunOut k1 = run_gpu(c, 1.0, cap, 5, 0);
        g_t3_tc.clear();
        const RunOut d = run_gpu(c, 1.0, cap, 5, 1);
        long negs = (d.min_rv < 0) + (d.min_rg < 0) + (d.min_Q < 0);
        bool fin = d.bad.empty(); for (int q = 0; q <= NVAR; ++q) fin = fin && std::isfinite(d.comp_ratio[q]);
        const size_t N = g_t3_tc.size(), W = (N == 0) ? 0 : (size_t)std::ceil(0.1*(double)N);
        double mx[6] = {0, 0, 0, 0, 0, 0}; bool cfin = true;
        for (size_t u = N - W; u < N; ++u) for (int k = 0; k < 6; ++k) { const double v = g_t3_tc[u].C[k]; if (!std::isfinite(v)) cfin = false; else mx[k] = std::max(mx[k], v); }
        bool cok = cfin && W > 0; for (int k = 0; k < 6; ++k) cok = cok && (mx[k] <= kappa);
        const bool ok = d.converged && fin && d.ratio_max <= 1.0 && negs == 0 && d.theta_final_lt1 == 0 && d.theta_src_final_lt1 == 0 && cok;
        char m[640];
        snprintf(m, sizeof m, "#4h (b) DT_MAX %.2g K キー 0: %s 反復 %d (キー 1 は %d)、独立残差の最大比 %.3f、min ρv %.2e ρg %.2e ρQ %.2e、最後の θ<1 %d・θ_src<1 %d、θ<1 のセル·反復 %ld、"
                 "末尾 %zu/%zu 更新の max C [w %.2e v %.2e g %.2e Q2 %.2e Q1 %.2e Q0 %.2e] (κ %.3e)、出口 g %.6e (キー 1 %.6e)",
                 dtm, d.converged ? "収束" : (d.bad.empty() ? "上限到達" : d.bad.c_str()), d.iters, k1.iters, d.ratio_max, d.min_rv, d.min_rg, d.min_Q,
                 d.theta_final_lt1, d.theta_src_final_lt1, d.theta_lt1, W, N, mx[0], mx[1], mx[2], mx[3], mx[4], mx[5], kappa, d.g_out, k1.g_out);
        verdict(ok, m);
    }
    DT_MAX = 1.0;
}

int main(int argc, char** argv)
{
    const bool quick = (argc > 1 && !strcmp(argv[1], "--quick"));
    Ctx c;
    {
        ResolvedSpeciesDB db = speciesDB_resolve(std::vector<std::string>{"N2", "H2O"}, "");
        c.sp[0] = db.species[0]; c.sp[1] = db.species[1];
        for (int k = 0; k < 2; ++k) if (!(c.sp[k].invMW > 0.0)) c.sp[k].invMW = 1.0/c.sp[k].MW;
        c.cp = cond_test_props_H2O(true);   // device 用の潜熱ペアも作る
        cond_tables_build_host(c.cp, c.ht); c.tb = cond_tables_view_host(c.ht);
        c.cpf = condProps_to_f(c.cp);
        c.Rw = c.cp.R; c.RN2 = THERMO_RU/c.sp[0].MW;
    }
    // device の種 (T1 は N2, O2, H2O; T2 エネルギー・T3 は N2, H2O を index 0, 1 で使う → T1 だけ O2 を挟んだ別配列)
    SpeciesThermoF* sp2 = nullptr; SpeciesThermoF* sp3 = nullptr;
    {
        ResolvedSpeciesDB db3 = speciesDB_resolve(std::vector<std::string>{"N2", "O2", "H2O"}, "");
        SpeciesThermoF h3[3], h2[2];
        for (int k = 0; k < 3; ++k) { SpeciesThermo s = db3.species[k]; if (!(s.invMW > 0.0)) s.invMW = 1.0/s.MW; g_sp3[k] = s; h3[k] = thermo_to_float(s); }
        for (int k = 0; k < 2; ++k) h2[k] = thermo_to_float(c.sp[k]);
        CK(cudaMalloc(&sp3, 3*sizeof(SpeciesThermoF))); CK(cudaMemcpy(sp3, h3, 3*sizeof(SpeciesThermoF), cudaMemcpyHostToDevice));
        CK(cudaMalloc(&sp2, 2*sizeof(SpeciesThermoF))); CK(cudaMemcpy(sp2, h2, 2*sizeof(SpeciesThermoF), cudaMemcpyHostToDevice));
    }
    G.cp = c.cp;
    printf("=== #4e 二相拡散 CUDA 実装の単体試験 (tp_face_flux / tp_vl_update) ===\n");
    G.sp = sp3; test_T1(c, quick ? 20 : 200);
    test_T2_structure();
    test_T2_identities(quick ? 2000 : 20000);
    G.sp = sp2;
    test_T2_nonneg_and_conservation(quick);
    test_4f_box(quick);
    test_T2_energy(c);
    test_4f_theta();
    test_T3(c, quick ? 3000 : 20000);
    test_4g_unit();
    test_4g_T3(c, quick ? 3000 : 20000);
    test_4h(c, quick ? 3000 : 20000);
    printf("\n%s\n", g_fail ? "FAIL あり" : "ALL PASS");
    return g_fail ? 1 : 0;
}
