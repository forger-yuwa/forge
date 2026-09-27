#pragma once

// =============================================================================
// transportMix_d.cuh
//   種ごとに出所を選ぶ輸送物性 (粘性 μ・熱伝導率 λ) と CEA 形の frozen 混合則 (double)。
//   plans/active/thermophysics-solver-owned-species-db.md §4.3b・§4.3c・§5.1 #5t2、仕様 methods/thermophysics.md。
//
//   **段 1 (#5t2): host だけが使う**。GPU カーネルからはまだ呼ばない (現行カーネルは thermo_d.cuh の
//   thermo_mu_mix / thermo_lambda_mix のまま)。段 2 で同じ式を GPU から呼べるよう THERMO_HD にしてある。
//   係数の解決・検査・記録は input/speciesTransportDB.{hpp,cpp}。
//
//   種ごとのモデル (SpeciesTransportD::model):
//     CEA      CEA trans.inp の種別フィット ln f = A ln T + B/T + C/T² + D (f: μ [μP]・λ [μW/(cm K)])。
//              区間は cea2.f TRANIN の kt と同じ選び方 (T ≤ Thi となる最初の区間、無ければ最後 = 範囲外は最寄り区間で外挿)。
//     FIT      ユーザ指定フィット (speciesDBFile の transport_fit)。式・単位・区間の選び方は CEA と同じ。
//     KINETIC  μ: Chapman–Enskog (Neufeld 1972 の Ω(2,2)*、T* を [0.3, 100] にクランプ; thermo_d.cuh と同じ範囲)。
//              双極子モーメントがあれば Brokaw の極性補正 Ω(2,2)* += 0.2 δ*²/T* (δ* = μ_d²/(2 ε σ³), 同じクランプ済み T*)。
//              λ: 修正 Eucken λ = μ (c_p + 1.25 R_s) (現行と同じ式; Warnatz 式は未決 §10)。
//     H2O_IAPWS_CEA_V1  (config 名 custom:h2o_iapws_cea_v1) H2O 専用:
//              T ≤ 500 K は IAPWS 希薄気体 (粘性 IAPWS 2008 μ₀、熱伝導 IAPWS 2011 λ₀)、T ≥ 700 K は CEA (上の CEA 式)、
//              500 < T < 700 K は ln μ・ln λ を smoothstep w = 3s² − 2s³ (s = (T−500)/200) で IAPWS → CEA。
//              150 K 未満は 150 K での対数勾配 n = d ln f/d ln T に合わせた冪 f(150)(T/150)^n (C¹)。
//              253.15 K 未満は IAPWS の公式適用域外 (起動ログに出す)。
//   相互作用粘性 η_ij (TransportPairD::kind; 規約は対称):
//     両種が KINETIC → 二元 Chapman–Enskog (σ_ij = (σ_i+σ_j)/2, ε_ij = √(ε_i ε_j), 換算質量 2M_iM_j/(M_i+M_j); 極性補正なし)
//     それ以外 → CEA trans.inp の相互作用フィット、無ければ CEA の剛体球近似 (cea2.f 行 5565–5570)。
//   混合則 (cea2.f TRANP, 行 5599–5634):
//     φ_ij = 2 M_j η_i / ((M_i+M_j) η_ij)、ψ_ij = φ_ij [1 + 2.41 (M_i−M_j)(M_i−0.142 M_j)/(M_i+M_j)²]、φ_ii = ψ_ii = 1
//     μ = Σ_i X_i η_i / Σ_j φ_ij X_j,   λ = Σ_i X_i λ_i / Σ_j ψ_ij X_j   (λ は frozen)
//   単位: 返り値は SI (μ [Pa s], λ [W/(m K)])。CEA 単位からは μP → ×1e-7、μW/(cm K) → ×1e-4。
// =============================================================================

#include "cuda_forge/thermo_d.cuh"   // THERMO_HD, SpeciesThermo, thermo_cp_mass, THERMO_RU

// 上限 (超過は起動時に拒否; input/speciesTransportDB.cpp)
#define TRANSPORT_MAX_FIT_INTERVALS 3    // CEA trans.inp と同じ (cea2.f UTRAN: nV, nC <= 3)
#define TRANSPORT_MAX_REAL_SPECIES  32   // lump 展開後の実種数

#define TRANSPORT_MODEL_CEA              1
#define TRANSPORT_MODEL_KINETIC          2
#define TRANSPORT_MODEL_FIT              3
#define TRANSPORT_MODEL_H2O_IAPWS_CEA_V1 4

#define TRANSPORT_PAIR_CE           1   // 二元 Chapman–Enskog (両種 kinetic)
#define TRANSPORT_PAIR_CEA          2   // CEA trans.inp の相互作用フィット
#define TRANSPORT_PAIR_RIGID_SPHERE 3   // CEA の剛体球近似

// 単位換算
#define TRANSPORT_MICROPOISE_TO_PAS 1.0e-7   // μP → Pa s
#define TRANSPORT_CEA_COND_TO_SI    1.0e-4   // μW/(cm K) → W/(m K)
// Chapman–Enskog 粘性の定数 (thermo_mu_species と同じ): μ[Pa s] = 2.6693e-6 √(M[g/mol] T) / (σ[Å]² Ω)
#define TRANSPORT_CE_MU_CONST 2.6693e-6

// custom:h2o_iapws_cea_v1 の定数 (版 v1。変えたら版を上げる; 記録と互換性ハッシュに入る)
#define TRANSPORT_H2O_TC          647.096   // IAPWS 臨界温度 [K]
#define TRANSPORT_H2O_BLEND_LO    500.0     // これ以下は IAPWS
#define TRANSPORT_H2O_BLEND_HI    700.0     // これ以上は CEA
#define TRANSPORT_H2O_POWER_BELOW 150.0     // これ未満は冪外挿
#define TRANSPORT_H2O_IAPWS_TMIN  253.15    // IAPWS の公式適用域の下端 (ログのみ)

// CEA 形フィット (1 物性ぶん)。行 k: Tlo[k] ≤ T ≤ Thi[k] で ln f = A ln T + B/T + C/T² + D。
struct TransportFitD {
    int    n;                                   // 区間数 (0 = データなし)
    double Tlo[TRANSPORT_MAX_FIT_INTERVALS];
    double Thi[TRANSPORT_MAX_FIT_INTERVALS];
    double A[TRANSPORT_MAX_FIT_INTERVALS];
    double B[TRANSPORT_MAX_FIT_INTERVALS];
    double C[TRANSPORT_MAX_FIT_INTERVALS];
    double D[TRANSPORT_MAX_FIT_INTERVALS];
};

// 実種 1 つの輸送モデル
struct SpeciesTransportD {
    int           model;       // TRANSPORT_MODEL_*
    double        MW;          // [kg/mol] (熱物性と同じ値)
    TransportFitD V;           // CEA / FIT / H2O_IAPWS_CEA_V1: 粘性フィット (μP)
    TransportFitD C;           // 同: 熱伝導フィット (μW/(cm K))
    double        sigma_LJ;    // KINETIC: [Å]
    double        eps_kB;      // KINETIC: [K]
    double        deltaStar;   // KINETIC: 換算双極子 δ* (0 = 非極性)
    SpeciesThermo thermo;      // KINETIC の修正 Eucken の c_p 用 (絶対基準係数; c_p は datum に依らない)
};

// 異種の組 (a < b は resolver が決めた向き; 剛体球近似の式は向きで丸めが変わるので固定する)
struct TransportPairD {
    int           kind;        // TRANSPORT_PAIR_*
    int           a, b;        // 実種 index
    TransportFitD V;           // CEA: 相互作用粘性フィット (μP)
    double        sigma_ab;    // CE: (σ_a+σ_b)/2 [Å]
    double        eps_ab;      // CE: √(ε_a ε_b) [K]
};

// 上三角 (a < b) の組の通し番号
THERMO_HD int transport_pair_index(int a, int b, int n)
{
    return a*n - a*(a + 1)/2 + (b - a - 1);
}

// CEA 形フィットの評価 (単位はフィットの単位のまま)。区間は cea2.f TRANIN と同じ選び方。
THERMO_HD double transport_fit_eval(const TransportFitD& f, double T)
{
    int k = f.n - 1;
    for (int i = 0; i < f.n; ++i) { if (T <= f.Thi[i]) { k = i; break; } }
    const double lnT = log(T);
    return exp(f.A[k]*lnT + f.B[k]/T + f.C[k]/(T*T) + f.D[k]);
}

// Neufeld (1972) の Ω(2,2)* (double; T* は [0.3, 100] にクランプ済みを渡す)
THERMO_HD double transport_omega22(double Ts)
{
    return 1.16145*pow(Ts, -0.14874) + 0.52487*exp(-0.77320*Ts) + 2.16178*exp(-2.43787*Ts);
}

THERMO_HD double transport_clamp_Tstar(double Ts)
{
    if (Ts < 0.3)   Ts = 0.3;
    if (Ts > 100.0) Ts = 100.0;
    return Ts;
}

// KINETIC 単成分粘性 [Pa s]
THERMO_HD double transport_mu_kinetic(const SpeciesTransportD& s, double T)
{
    const double Ts = transport_clamp_Tstar(T/s.eps_kB);
    const double om = transport_omega22(Ts) + 0.2*s.deltaStar*s.deltaStar/Ts;
    return TRANSPORT_CE_MU_CONST*sqrt(s.MW*1000.0*T)/(s.sigma_LJ*s.sigma_LJ*om);
}

// IAPWS 2008 希薄気体粘性 μ₀ [Pa s] と IAPWS 2011 希薄気体熱伝導率 λ₀ [W/(m K)]、およびその対数勾配 d ln f/d ln T
THERMO_HD void transport_iapws_h2o(double T, double* mu, double* lam, double* nmu, double* nlam)
{
    const double H[4] = {1.67752, 2.20462, 0.6366564, -0.241605};
    const double L[5] = {2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4};
    const double Tb = T/TRANSPORT_H2O_TC;
    double sH = 0.0, kH = 0.0, sL = 0.0, kL = 0.0, p = 1.0;   // p = Tb^-k
    for (int k = 0; k < 5; ++k) {
        if (k < 4) { sH += H[k]*p; kH += k*H[k]*p; }
        sL += L[k]*p; kL += k*L[k]*p;
        p /= Tb;
    }
    *mu  = 100.0*sqrt(Tb)/sH*1.0e-6;     // μPa s → Pa s
    *lam = sqrt(Tb)/sL*1.0e-3;           // mW/(m K) → W/(m K)
    *nmu  = 0.5 + kH/sH;                 // d ln μ₀/d ln T
    *nlam = 0.5 + kL/sL;
}

// custom:h2o_iapws_cea_v1 [SI]
THERMO_HD void transport_h2o_iapws_cea_v1(const SpeciesTransportD& s, double T, double* mu, double* lam)
{
    double mI, lI, nm, nl;
    if (T >= TRANSPORT_H2O_BLEND_HI) {
        *mu  = transport_fit_eval(s.V, T)*TRANSPORT_MICROPOISE_TO_PAS;
        *lam = transport_fit_eval(s.C, T)*TRANSPORT_CEA_COND_TO_SI;
        return;
    }
    if (T < TRANSPORT_H2O_POWER_BELOW) {
        const double T0 = TRANSPORT_H2O_POWER_BELOW;
        transport_iapws_h2o(T0, &mI, &lI, &nm, &nl);
        const double r = log(T/T0);
        *mu  = mI*exp(nm*r);
        *lam = lI*exp(nl*r);
        return;
    }
    transport_iapws_h2o(T, &mI, &lI, &nm, &nl);
    if (T <= TRANSPORT_H2O_BLEND_LO) { *mu = mI; *lam = lI; return; }
    const double s01 = (T - TRANSPORT_H2O_BLEND_LO)/(TRANSPORT_H2O_BLEND_HI - TRANSPORT_H2O_BLEND_LO);
    const double w   = s01*s01*(3.0 - 2.0*s01);
    const double mC  = transport_fit_eval(s.V, T)*TRANSPORT_MICROPOISE_TO_PAS;
    const double lC  = transport_fit_eval(s.C, T)*TRANSPORT_CEA_COND_TO_SI;
    *mu  = exp((1.0 - w)*log(mI) + w*log(mC));
    *lam = exp((1.0 - w)*log(lI) + w*log(lC));
}

// 実種 1 つの μ [Pa s]・λ [W/(m K)]
THERMO_HD void transport_species(const SpeciesTransportD& s, double T, double* mu, double* lam)
{
    switch (s.model) {
    case TRANSPORT_MODEL_KINETIC: {
        const double m = transport_mu_kinetic(s, T);
        *mu  = m;
        *lam = m*(thermo_cp_mass(s.thermo, T) + 1.25*THERMO_RU/s.MW);
        return;
    }
    case TRANSPORT_MODEL_H2O_IAPWS_CEA_V1:
        transport_h2o_iapws_cea_v1(s, T, mu, lam);
        return;
    default:   // CEA / FIT
        *mu  = transport_fit_eval(s.V, T)*TRANSPORT_MICROPOISE_TO_PAS;
        *lam = transport_fit_eval(s.C, T)*TRANSPORT_CEA_COND_TO_SI;
        return;
    }
}

// 相互作用粘性 η_ab [Pa s]。eta_a, eta_b は各種の単成分粘性 [Pa s] (剛体球近似で使う)。
THERMO_HD double transport_eta_pair(const TransportPairD& p, const SpeciesTransportD* sp, const double* eta, double T)
{
    const SpeciesTransportD& A = sp[p.a];
    const SpeciesTransportD& B = sp[p.b];
    if (p.kind == TRANSPORT_PAIR_CE) {
        const double Ts = transport_clamp_Tstar(T/p.eps_ab);
        const double Mr = 2.0*A.MW*B.MW/(A.MW + B.MW)*1000.0;   // 換算質量 ×2 [g/mol]
        return TRANSPORT_CE_MU_CONST*sqrt(Mr*T)/(p.sigma_ab*p.sigma_ab*transport_omega22(Ts));
    }
    if (p.kind == TRANSPORT_PAIR_CEA) {
        return transport_fit_eval(p.V, T)*TRANSPORT_MICROPOISE_TO_PAS;
    }
    // 剛体球近似 (cea2.f 5565–5570; i = a, j = b)
    const double ea = eta[p.a], eb = eta[p.b];
    const double ratio = sqrt(B.MW/A.MW);
    const double d = 1.0 + sqrt(ratio*ea/eb);
    return 5.656854*ea*sqrt(B.MW/(A.MW + B.MW))/(d*d);
}

// CEA frozen 混合則。X: 実種のモル分率 (n 個)。pairs: 上三角 n(n-1)/2 個 (transport_pair_index の順)。
// mu_i, lam_i (任意, nullptr 可): 実種ごとの単成分値の出力。
THERMO_HD void transport_mix(const SpeciesTransportD* sp, const TransportPairD* pairs, int n, const double* X, double T,
                             double* mu, double* lam, double* mu_i = nullptr, double* lam_i = nullptr)
{
    double eta[TRANSPORT_MAX_REAL_SPECIES], con[TRANSPORT_MAX_REAL_SPECIES];
    for (int i = 0; i < n; ++i) transport_species(sp[i], T, &eta[i], &con[i]);
    double m = 0.0, l = 0.0;
    for (int i = 0; i < n; ++i) {
        if (mu_i)  mu_i[i]  = eta[i];
        if (lam_i) lam_i[i] = con[i];
        if (!(X[i] > 0.0)) continue;
        double sv = 0.0, sc = 0.0;
        for (int j = 0; j < n; ++j) {
            double phi = 1.0, psi = 1.0;
            if (j != i) {
                const int a = (i < j) ? i : j, b = (i < j) ? j : i;
                const double eij = transport_eta_pair(pairs[transport_pair_index(a, b, n)], sp, eta, T);
                const double Mi = sp[i].MW, Mj = sp[j].MW;
                phi = 2.0*Mj*eta[i]/((Mi + Mj)*eij);
                psi = phi*(1.0 + 2.41*(Mi - Mj)*(Mi - 0.142*Mj)/((Mi + Mj)*(Mi + Mj)));
            }
            sv += phi*X[j];
            sc += psi*X[j];
        }
        m += eta[i]*X[i]/sv;
        l += con[i]*X[i]/sc;
    }
    *mu = m;
    *lam = l;
}
