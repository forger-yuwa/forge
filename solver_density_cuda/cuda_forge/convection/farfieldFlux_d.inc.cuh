// =============================================================================
// 遠方境界 `farfield` の境界半割面流束 (node 専用)。plans/active/boundary-node-farfield-characteristic.md §4.2–4.3。
//   - convectiveFlux_d.cu に textual include される断片 (単一 TU を保つため、CMake の source に足さない)。
//   - node にゴーストセルは無い。境界節点 i の状態 U_i と、その場で作る外側状態 U_R から、半割面の流束を HLLC で作り
//     境界節点の残差に加える (弱形式)。外側状態は保存しない。
//   - 外側状態 (§4.2): 圧力・法線速度 = 内部エントロピーの 2 膨張波近似 (TRRS、frozen γ_i)、密度・接線速度・組成・k・ω = 自由流、
//     超音速の境目は原始変数の滑らかな混合 (自由流の超音速流入 → 自由流、内部の超音速流出 → 内部 [優先])、
//     真空・非物理は U_R = U_i に置換して数える。
//   - 面流束 (§4.3): HLLC (Toro)、波速 Davis。S_L ≤ S* ≤ S_R と星密度 > 0 を検査し、外れたら HLL に退避して数える。
//   - 化学種・k・ω は本カーネルでは残差に足さない (各輸送カーネルが massflux[ip] と下の面の値配列で運ぶ)。
//   - 外側状態と HLLC は __host__ __device__ の純粋関数 (ホスト単体試験 tests/unit/test_farfield_flux.cpp から呼ぶ)。
// =============================================================================
#pragma once

#ifndef FF_HD
#  if defined(__CUDACC__)
#    define FF_HD __host__ __device__ inline
#  else
#    define FF_HD inline
#  endif
#endif

#include <cmath>

// 1 側の原始状態 (倍精度で扱う。境界面は数が少ないので精度を優先)
struct FfPrim {
    double r, u[3], p;           // 密度・速度・静圧
    double Y[THERMO_MAX_SPECIES]; // 質量分率 (nSpecies >= 2 のときだけ使う)
    double k, om;                // 乱流 (無いときは 0)
    double T, e, c;              // 熱力学 (ff_thermo が埋める): 温度・内部エネルギー (datum はソルバと同じ)・音速
};

// 熱力学モデル (CPG / 単成分 TP / 多成分 TP) を ff_thermo に渡す
struct FfGas {
    int thermalMethod;           // 2 = TP (NASA)、それ以外 = CPG
    double ga, cp;               // CPG
    const SpeciesThermo* sp;     // TP
    int nSpecies;
};

FF_HD void ff_thermo(const FfGas& g, FfPrim& s)
{
    if (g.thermalMethod == 2) {
        const bool mix = (g.nSpecies >= 2);
        const int n = mix ? g.nSpecies : 1;
        const double R = mix ? thermo_R_mix(g.sp, n, s.Y) : thermo_R_species(g.sp[0]);
        s.T = s.p / (s.r * R);
        s.e = mix ? thermo_e_mix(g.sp, n, s.Y, s.T) : (thermo_h_mass(g.sp[0], s.T) - R * s.T);
        const double cpv = mix ? thermo_cp_mix(g.sp, n, s.Y, s.T) : thermo_cp_mass(g.sp[0], s.T);
        const double gam = cpv / fmax(cpv - R, 1.0e-6);
        s.c = sqrt(gam * s.p / s.r);
    } else {
        const double R = g.cp - g.cp / g.ga;
        s.T = s.p / (s.r * R);
        s.e = s.p / ((g.ga - 1.0) * s.r);
        s.c = sqrt(g.ga * s.p / s.r);
    }
}

FF_HD double ff_smooth01(double x)
{
    const double t = fmin(fmax(x, 0.0), 1.0);
    return t * t * (3.0 - 2.0 * t);
}

FF_HD void ff_blend(FfPrim& a, const FfPrim& b, double w, int nY)
{
    if (w <= 0.0) return;
    a.r = (1.0 - w) * a.r + w * b.r;
    for (int d = 0; d < 3; ++d) a.u[d] = (1.0 - w) * a.u[d] + w * b.u[d];
    a.p = (1.0 - w) * a.p + w * b.p;
    for (int s = 0; s < nY; ++s) a.Y[s] = (1.0 - w) * a.Y[s] + w * b.Y[s];
    a.k = (1.0 - w) * a.k + w * b.k;
    a.om = (1.0 - w) * a.om + w * b.om;
}

// 原始状態の複写。構造体代入 (R = f) は device で組成 Y が 0 になる事例があった (2026-09-29, V1b 側面) ので、使う成分を明示的に写す
FF_HD void ff_copy(FfPrim& a, const FfPrim& b, int nY)
{
    a.r = b.r; a.p = b.p; a.k = b.k; a.om = b.om; a.T = b.T; a.e = b.e; a.c = b.c;
    for (int d = 0; d < 3; ++d) a.u[d] = b.u[d];
    for (int s = 0; s < nY; ++s) a.Y[s] = b.Y[s];
}

// 外側状態 (§4.2)。i = 境界節点 (ff_thermo 済み)、f = 自由流 (ff_thermo 済み)、n = 外向き単位法線。
// 戻り値: 0 = 通常、非 0 = U_R = U_i に置換した (1 = 真空判定、2 = 密度・圧力が非正/非有限、3 = 熱物性が非有限)。R は ff_thermo 済みで返す。
FF_HD int ff_outer_state(const FfGas& g, const FfPrim& i, const FfPrim& f, const double n[3], FfPrim& R, double band = 0.1)
{
    const int nY = (g.thermalMethod == 2 && g.nSpecies >= 2) ? g.nSpecies : 0;
    const double uni = i.u[0] * n[0] + i.u[1] * n[1] + i.u[2] * n[2];
    const double unf = f.u[0] * n[0] + f.u[1] * n[1] + f.u[2] * n[2];
    const double gi = i.c * i.c * i.r / i.p;          // frozen γ_i
    const double gf = f.c * f.c * f.r / f.p;          // frozen γ_∞
    const double z = (gi - 1.0) / (2.0 * gi);
    const double cpo = i.c * pow(f.p / i.p, z);        // 擬似外側 (P∞、内部エントロピー) の音速
    const double B = i.c + cpo - 0.5 * (gi - 1.0) * (unf - uni);
    bool bad = !(B > 1.0e-3 * (i.c + cpo));
    // 圧力・法線速度 (TRRS)。密度・接線速度・組成・k・ω は自由流側
    ff_copy(R, f, nY);
    if (!bad) {
        const double ps = pow(B / (i.c * pow(i.p, -z) + cpo * pow(f.p, -z)), 1.0 / z);
        const double us = uni + 2.0 * i.c / (gi - 1.0) * (1.0 - pow(ps / i.p, z));
        R.p = ps;
        R.r = f.r * pow(ps / f.p, 1.0 / gf);
        for (int d = 0; d < 3; ++d) R.u[d] = f.u[d] + (us - unf) * n[d];
    }
    // 超音速の境目: 自由流の超音速流入へ、続いて内部の超音速流出へ (後者が優先)
    ff_blend(R, f, ff_smooth01(((-1.0 + band) - unf / f.c) / band), nY);
    ff_blend(R, i, ff_smooth01((uni / i.c - (1.0 - band)) / band), nY);
    if (bad) { ff_copy(R, i, nY); return 1; }
    if (!(R.r > 0.0 && R.p > 0.0 && std::isfinite(R.r) && std::isfinite(R.p))) { ff_copy(R, i, nY); return 2; }
    ff_thermo(g, R);
    if (!(R.c > 0.0 && std::isfinite(R.c) && std::isfinite(R.e))) { ff_copy(R, i, nY); return 3; }
    return 0;
}

// HLLC 流束 (Toro、Davis 波速)。F[0] = 質量、F[1..3] = 運動量 (圧力は p − pRef)、F[4] = エネルギー、いずれも単位面積あたり。
// energyK: SST 全エネルギー (E* = E + k、p* = p + 2/3 ρk)。戻り値: 0 = HLLC、1 = HLL に退避。
FF_HD int ff_hllc(const FfPrim& L, const FfPrim& R, const double n[3], double pRef, int energyK, double F[5])
{
    const double unL = L.u[0] * n[0] + L.u[1] * n[1] + L.u[2] * n[2];
    const double unR = R.u[0] * n[0] + R.u[1] * n[1] + R.u[2] * n[2];
    const double pL = L.p + (energyK ? 2.0 / 3.0 * L.r * L.k : 0.0);
    const double pR = R.p + (energyK ? 2.0 / 3.0 * R.r * R.k : 0.0);
    const double EL = L.e + 0.5 * (L.u[0] * L.u[0] + L.u[1] * L.u[1] + L.u[2] * L.u[2]) + (energyK ? L.k : 0.0);
    const double ER = R.e + 0.5 * (R.u[0] * R.u[0] + R.u[1] * R.u[1] + R.u[2] * R.u[2]) + (energyK ? R.k : 0.0);
    const double SL = fmin(unL - L.c, unR - R.c);
    const double SR = fmax(unL + L.c, unR + R.c);
    const double den = L.r * (SL - unL) - R.r * (SR - unR);
    const double Sm = (pR - pL + L.r * unL * (SL - unL) - R.r * unR * (SR - unR)) / den;
    const double rsL = L.r * (SL - unL) / (SL - Sm);
    const double rsR = R.r * (SR - unR) / (SR - Sm);
    const bool ok = (SL <= Sm) && (Sm <= SR) && (rsL > 0.0) && (rsR > 0.0) && std::isfinite(Sm);
    double FL[5], FR[5], UL[5], UR[5];
    FL[0] = L.r * unL; FR[0] = R.r * unR;
    for (int d = 0; d < 3; ++d) { FL[1 + d] = L.r * L.u[d] * unL + (pL - pRef) * n[d]; FR[1 + d] = R.r * R.u[d] * unR + (pR - pRef) * n[d]; }
    FL[4] = (L.r * EL + pL) * unL; FR[4] = (R.r * ER + pR) * unR;
    UL[0] = L.r; UR[0] = R.r;
    for (int d = 0; d < 3; ++d) { UL[1 + d] = L.r * L.u[d]; UR[1 + d] = R.r * R.u[d]; }
    UL[4] = L.r * EL; UR[4] = R.r * ER;
    if (SL >= 0.0) { for (int q = 0; q < 5; ++q) F[q] = FL[q]; return ok ? 0 : 1; }
    if (SR <= 0.0) { for (int q = 0; q < 5; ++q) F[q] = FR[q]; return ok ? 0 : 1; }
    if (!ok) {
        for (int q = 0; q < 5; ++q) F[q] = (SR * FL[q] - SL * FR[q] + SL * SR * (UR[q] - UL[q])) / (SR - SL);
        return 1;
    }
    const bool left = (Sm >= 0.0);
    const FfPrim& K = left ? L : R;
    const double S = left ? SL : SR, unK = left ? unL : unR, pK = left ? pL : pR, EK = left ? EL : ER;
    const double* FK = left ? FL : FR; const double* UK = left ? UL : UR;
    const double f = K.r * (S - unK) / (S - Sm);
    double Us[5];
    Us[0] = f;
    for (int d = 0; d < 3; ++d) Us[1 + d] = f * (K.u[d] + (Sm - unK) * n[d]);
    Us[4] = f * (EK + (Sm - unK) * (Sm + pK / (K.r * (S - unK))));
    for (int q = 0; q < 5; ++q) F[q] = FK[q] + S * (Us[q] - UK[q]);
    return 0;
}

#if defined(__CUDACC__)
__device__ unsigned long long g_ffHll = 0;      // HLL 退避の面数 (累積)
__device__ unsigned long long g_ffVac = 0;      // 真空・非物理で U_R = U_i にした面数 (累積)
#define FF_DUMP_NF 24

// 自由流 (bcond floats、面によらず一定)
struct FfInf {
    double r, u[3], p, k, om;
    double Y[THERMO_MAX_SPECIES];
};

__global__ void farfield_flux_d(
    geom_int nb, const geom_int* bplane_plane, const geom_int* bplane_cell,
    const geom_float* sx, const geom_float* sy, const geom_float* sz, const geom_float* ss,
    const flow_float* ro, const flow_float* Ux, const flow_float* Uy, const flow_float* Uz, const flow_float* P,
    flow_float* const* roY, const flow_float* kturb, const flow_float* omega,
    FfGas gas, FfInf inf, flow_float pRef, int energyK,
    flow_float* res_ro, flow_float* res_roUx, flow_float* res_roUy, flow_float* res_roUz, flow_float* res_roe,
    flow_float* massflux,
    flow_float* const* faceY, flow_float* faceK, flow_float* faceOm,
    float* dumpBuf)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib >= nb) return;
    const geom_int ip = bplane_plane[ib];
    const geom_int ic = bplane_cell[ib];
    const double S = ss[ip];
    const double nrm[3] = {sx[ip] / S, sy[ip] / S, sz[ip] / S};
    const int nY = (gas.thermalMethod == 2 && gas.nSpecies >= 2) ? gas.nSpecies : 0;

    FfPrim I = {}, Fs = {}, Rs = {};
    I.r = ro[ic]; I.u[0] = Ux[ic]; I.u[1] = Uy[ic]; I.u[2] = Uz[ic]; I.p = P[ic];
    if (nY > 0) {
        double sum = 0.0;
        for (int s = 0; s < nY; ++s) { I.Y[s] = fmax((double)roY[s][ic], 0.0); sum += I.Y[s]; }
        for (int s = 0; s < nY; ++s) I.Y[s] /= fmax(sum, 1.0e-300);
    }
    I.k = (kturb != nullptr) ? (double)kturb[ic] : 0.0;
    I.om = (omega != nullptr) ? (double)omega[ic] : 0.0;
    ff_thermo(gas, I);
    Fs.r = inf.r; Fs.u[0] = inf.u[0]; Fs.u[1] = inf.u[1]; Fs.u[2] = inf.u[2]; Fs.p = inf.p; Fs.k = inf.k; Fs.om = inf.om;
    for (int s = 0; s < nY; ++s) Fs.Y[s] = inf.Y[s];
    ff_thermo(gas, Fs);

    const int vac = ff_outer_state(gas, I, Fs, nrm, Rs);
    double F[5];
    const int hll = ff_hllc(I, Rs, nrm, (double)pRef, energyK, F);
    bool finite = true;
    for (int q = 0; q < 5; ++q) finite = finite && std::isfinite(F[q]);
    if (!finite) { for (int q = 0; q < 5; ++q) F[q] = 0.0; }
    if (vac || hll || !finite) {
        if (vac) atomicAdd(&g_ffVac, 1ULL);
        if (hll || !finite) atomicAdd(&g_ffHll, 1ULL);
    }
    atomicAdd(&res_ro[ic],   (flow_float)(-F[0] * S));
    atomicAdd(&res_roUx[ic], (flow_float)(-F[1] * S));
    atomicAdd(&res_roUy[ic], (flow_float)(-F[2] * S));
    atomicAdd(&res_roUz[ic], (flow_float)(-F[3] * S));
    atomicAdd(&res_roe[ic],  (flow_float)(-F[4] * S));
    massflux[ip] = (flow_float)(F[0] * S);
    // 流入 (質量流束 < 0) で運ぶスカラーの面の値 = 外側状態の値 (流出では各輸送カーネルが内部値を使う)
    for (int s = 0; s < nY; ++s) faceY[s][ip] = (flow_float)Rs.Y[s];
    if (faceK != nullptr) faceK[ip] = (flow_float)Rs.k;
    if (faceOm != nullptr) faceOm[ip] = (flow_float)Rs.om;
    if (dumpBuf != nullptr) {
        float* o = dumpBuf + (size_t)ib * FF_DUMP_NF;
        o[0] = (float)ip; o[1] = (float)ic; o[2] = (float)nrm[0]; o[3] = (float)nrm[1]; o[4] = (float)nrm[2]; o[5] = (float)S;
        for (int q = 0; q < 5; ++q) o[6 + q] = (float)(F[q] * S);
        o[11] = (float)Rs.r; o[12] = (float)Rs.u[0]; o[13] = (float)Rs.u[1]; o[14] = (float)Rs.u[2]; o[15] = (float)Rs.p;
        o[16] = (float)Rs.k; o[17] = (float)Rs.om; o[18] = (float)(nY > 0 ? Rs.Y[0] : 0.0);
        o[19] = (float)vac;   // 置換の理由 (0 = なし、1 = 真空判定、2 = 密度・圧力、3 = 熱物性)
        o[20] = (float)hll; o[21] = (float)pRef; o[22] = (float)I.c; o[23] = (float)Rs.c;
    }
}
#endif
