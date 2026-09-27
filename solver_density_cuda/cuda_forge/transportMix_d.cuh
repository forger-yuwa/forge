#pragma once

// =============================================================================
// transportMix_d.cuh
//   種ごとに出所を選ぶ輸送物性 (粘性 μ・熱伝導率 λ) と CEA 形の frozen 混合則 (double)。
//   plans/active/thermophysics-solver-owned-species-db.md §4.3b・§4.3c・§5.1 #5t2、仕様 methods/thermophysics.md。
//
//   host (解決・記録・試験) と GPU (段 2, #5t2-2) が同じ式を使う (THERMO_HD)。GPU では physProp.transport があり
//   viscMethod 2 のときだけ、セル (gasProperties_d) と壁 (wmlesWallModel_d) が transport_mix_Y を呼ぶ。
//   書かれていなければ現行 (thermo_d.cuh の thermo_mu_mix / thermo_lambda_mix) のまま。
//   係数の解決・検査・記録は input/speciesTransportDB.{hpp,cpp}、device への転送は thermo_d.cu (thermo_init_db)。
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
//              253.15 K (IAPWS の公式適用域の下端) 未満は、253.15 K での IAPWS の対数勾配 n = d ln f/d ln T に合わせた冪
//              f(253.15)(T/253.15)^n (C¹; n は式から計算、μ ≈ 0.749・λ ≈ 0.975)。IAPWS μ₀ は 202 K で最小・134 K に極があり、式のままの外挿は非物理 (2026-09-27 ユーザ決定)。
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
#define TRANSPORT_H2O_POWER_BELOW 253.15    // これ未満は冪外挿 (IAPWS の公式適用域の下端)
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

// 温度とその関数 (混合則の 1 回の評価で使い回す; GPU の double の log・除算を減らす #5t2-2)
struct TransportTemp {
    double T, lnT, iT;   // T, ln T, 1/T
};

THERMO_HD TransportTemp transport_temp(double T)
{
    TransportTemp t;
    t.T = T;
    t.lnT = log(T);
    t.iT = 1.0/T;
    return t;
}

// CEA 形フィットの評価 (単位はフィットの単位のまま)。区間は cea2.f TRANIN と同じ選び方。
THERMO_HD double transport_fit_eval(const TransportFitD& f, const TransportTemp& t)
{
    int k = f.n - 1;
    for (int i = 0; i < f.n; ++i) { if (t.T <= f.Thi[i]) { k = i; break; } }
    return exp(f.A[k]*t.lnT + (f.B[k] + f.C[k]*t.iT)*t.iT + f.D[k]);
}

THERMO_HD double transport_fit_eval(const TransportFitD& f, double T)
{
    return transport_fit_eval(f, transport_temp(T));
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
THERMO_HD void transport_h2o_iapws_cea_v1(const SpeciesTransportD& s, const TransportTemp& tt, double* mu, double* lam)
{
    const double T = tt.T;
    double mI, lI, nm, nl;
    if (T >= TRANSPORT_H2O_BLEND_HI) {
        *mu  = transport_fit_eval(s.V, tt)*TRANSPORT_MICROPOISE_TO_PAS;
        *lam = transport_fit_eval(s.C, tt)*TRANSPORT_CEA_COND_TO_SI;
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
    const double mC  = transport_fit_eval(s.V, tt)*TRANSPORT_MICROPOISE_TO_PAS;
    const double lC  = transport_fit_eval(s.C, tt)*TRANSPORT_CEA_COND_TO_SI;
    *mu  = exp((1.0 - w)*log(mI) + w*log(mC));
    *lam = exp((1.0 - w)*log(lI) + w*log(lC));
}

THERMO_HD void transport_h2o_iapws_cea_v1(const SpeciesTransportD& s, double T, double* mu, double* lam)
{
    transport_h2o_iapws_cea_v1(s, transport_temp(T), mu, lam);
}

// 実種 1 つの μ [Pa s]・λ [W/(m K)]
THERMO_HD void transport_species(const SpeciesTransportD& s, const TransportTemp& tt, double* mu, double* lam)
{
    const double T = tt.T;
    switch (s.model) {
    case TRANSPORT_MODEL_KINETIC: {
        const double m = transport_mu_kinetic(s, T);
        *mu  = m;
        *lam = m*(thermo_cp_mass(s.thermo, T) + 1.25*THERMO_RU/s.MW);
        return;
    }
    case TRANSPORT_MODEL_H2O_IAPWS_CEA_V1:
        transport_h2o_iapws_cea_v1(s, tt, mu, lam);
        return;
    default:   // CEA / FIT
        *mu  = transport_fit_eval(s.V, tt)*TRANSPORT_MICROPOISE_TO_PAS;
        *lam = transport_fit_eval(s.C, tt)*TRANSPORT_CEA_COND_TO_SI;
        return;
    }
}

THERMO_HD void transport_species(const SpeciesTransportD& s, double T, double* mu, double* lam)
{
    transport_species(s, transport_temp(T), mu, lam);
}

// 相互作用粘性 η_ab [Pa s]。eta_a, eta_b は各種の単成分粘性 [Pa s] (剛体球近似で使う)。
THERMO_HD double transport_eta_pair(const TransportPairD& p, const SpeciesTransportD* sp, const double* eta, const TransportTemp& tt)
{
    const double T = tt.T;
    const SpeciesTransportD& A = sp[p.a];
    const SpeciesTransportD& B = sp[p.b];
    if (p.kind == TRANSPORT_PAIR_CE) {
        const double Ts = transport_clamp_Tstar(T/p.eps_ab);
        const double Mr = 2.0*A.MW*B.MW/(A.MW + B.MW)*1000.0;   // 換算質量 ×2 [g/mol]
        return TRANSPORT_CE_MU_CONST*sqrt(Mr*T)/(p.sigma_ab*p.sigma_ab*transport_omega22(Ts));
    }
    if (p.kind == TRANSPORT_PAIR_CEA) {
        return transport_fit_eval(p.V, tt)*TRANSPORT_MICROPOISE_TO_PAS;
    }
    // 剛体球近似 (cea2.f 5565–5570; i = a, j = b)
    const double ea = eta[p.a], eb = eta[p.b];
    const double ratio = sqrt(B.MW/A.MW);
    const double d = 1.0 + sqrt(ratio*ea/eb);
    return 5.656854*ea*sqrt(B.MW/(A.MW + B.MW))/(d*d);
}

THERMO_HD double transport_eta_pair(const TransportPairD& p, const SpeciesTransportD* sp, const double* eta, double T)
{
    return transport_eta_pair(p, sp, eta, transport_temp(T));
}

// CEA frozen 混合則。X: 実種のモル分率 (n 個)。pairs: 上三角 n(n-1)/2 個 (transport_pair_index の順)。
// mu_i, lam_i (任意, nullptr 可): 実種ごとの単成分値の出力。
//   η_ij は各組 (a<b) を 1 回だけ評価し、φ_ab・φ_ba (ψ も) を両方向の分母へ加える (#5t2-2)。
//   分母 Σ_j φ_ij X_j の加算順は j = 0..n−1 (自己項 φ_ii X_i は j = i の位置) で、組ごとに再評価していた段 1 の形と同じ。
//   X_i ≤ 0 の種は分子に入らないので分母も作らない。分母へ入る項が φ·0 だけの組 (相手の X が 0) は η_ij を評価しない
//   (0 を足すだけなので値は変わらない)。ln T・1/T は 1 回だけ計算して全フィットで使い回す。
THERMO_HD void transport_mix(const SpeciesTransportD* sp, const TransportPairD* pairs, int n, const double* X, double T,
                             double* mu, double* lam, double* mu_i = nullptr, double* lam_i = nullptr)
{
    const TransportTemp tt = transport_temp(T);
    double eta[TRANSPORT_MAX_REAL_SPECIES], con[TRANSPORT_MAX_REAL_SPECIES];
    double sv[TRANSPORT_MAX_REAL_SPECIES], sc[TRANSPORT_MAX_REAL_SPECIES];
    for (int i = 0; i < n; ++i) {
        // X_i = 0 の種の単成分値はどこにも入らない (分子にも、0 でない分母項を持つ組にも) ので評価しない
        if (X[i] != 0.0 || mu_i || lam_i) transport_species(sp[i], tt, &eta[i], &con[i]);
        else { eta[i] = 0.0; con[i] = 0.0; }
        if (mu_i)  mu_i[i]  = eta[i];
        if (lam_i) lam_i[i] = con[i];
        sv[i] = 0.0;
        sc[i] = 0.0;
    }
    for (int a = 0; a < n; ++a) {
        const bool pa = (X[a] > 0.0);
        if (pa) { sv[a] += 1.0*X[a]; sc[a] += 1.0*X[a]; }   // φ_aa = ψ_aa = 1 (j = a の位置)
        for (int b = a + 1; b < n; ++b) {
            const bool pb = (X[b] > 0.0);
            const bool ua = pa && X[b] != 0.0;   // 分母 a に 0 でない項が入る
            const bool ub = pb && X[a] != 0.0;   // 分母 b に 0 でない項が入る
            if (!ua && !ub) continue;
            const double eab = transport_eta_pair(pairs[transport_pair_index(a, b, n)], sp, eta, tt);
            const double Ma = sp[a].MW, Mb = sp[b].MW;
            const double is = 1.0/(Ma + Mb);
            const double ie = 1.0/eab;
            if (ua) {   // i = a, j = b
                const double phi = 2.0*Mb*eta[a]*is*ie;
                const double psi = phi*(1.0 + 2.41*(Ma - Mb)*(Ma - 0.142*Mb)*is*is);
                sv[a] += phi*X[b];
                sc[a] += psi*X[b];
            }
            if (ub) {   // i = b, j = a
                const double phi = 2.0*Ma*eta[b]*is*ie;
                const double psi = phi*(1.0 + 2.41*(Mb - Ma)*(Mb - 0.142*Ma)*is*is);
                sv[b] += phi*X[a];
                sc[b] += psi*X[a];
            }
        }
    }
    double m = 0.0, l = 0.0;
    for (int i = 0; i < n; ++i) {
        if (!(X[i] > 0.0)) continue;
        m += eta[i]*X[i]/sv[i];
        l += con[i]*X[i]/sc[i];
    }
    *mu = m;
    *lam = l;
}

// -----------------------------------------------------------------------------
// GPU で使う表 (thermo_d.cu の thermo_init_db が device へ上げる; 値渡しでカーネルへ)。
//   expand: 輸送する種 s (physProp.species の順, nTransported 個) × 実種 r の密行列 (行優先)。
//           値は lump 内モル分率 (lump でない種は 1、含まない実種は 0)。重複実種は同じ列に入る。
// -----------------------------------------------------------------------------
struct TransportTableD {
    int                      nReal;
    int                      nTransported;
    const SpeciesTransportD* sp;      // [nReal]
    const TransportPairD*    pairs;   // [nReal(nReal-1)/2]
    const double*            expand;  // [nTransported*nReal]
};

// 輸送種のモル分率 Xs → 実種のモル分率 Xreal (Xreal_r = Σ_s Xs_s·expand[s,r]; 加算順は s 昇順 = host の疎な展開と同じ値)
THERMO_HD void transport_expand_X(const TransportTableD& t, const double* Xs, double* Xreal)
{
    for (int r = 0; r < t.nReal; ++r) Xreal[r] = 0.0;
    for (int s = 0; s < t.nTransported; ++s) {
        const double* e = t.expand + s*t.nReal;
        for (int r = 0; r < t.nReal; ++r) Xreal[r] += Xs[s]*e[r];
    }
}

// セル・壁の共通入口 (#5t2-2)。組成は**モル基底で展開**する:
//   正規化した輸送種 Y (負値は 0 に切る) → 輸送種 X = (Y/M)/Σ(Y/M) (M: 輸送種の分子量 = 熱物性の MW; lump は合成 MW)
//   → Xreal_r = Σ_s X_s·expand[s,r] → transport_mix (CEA frozen)。
//   Y を直接展開行列に掛けてはいけない (行列は lump 内モル分率; N2/He lump で X_He 0.6 が 0.887 になる — codex 検算)。
//   th: 輸送種の熱物性 (MW だけ使う)。nY < 2 (単成分) のときは Y を読まず X = {1}。
//   Xreal_out (任意): 展開後の実種モル分率 (試験用)。
THERMO_HD void transport_mix_Y(const SpeciesThermo* th, const TransportTableD& t, int nY, const double* Y, double T,
                               double* mu, double* lam, double* Xreal_out = nullptr)
{
    double Xs[THERMO_MAX_SPECIES];
    double Xr[TRANSPORT_MAX_REAL_SPECIES];
    if (nY < 2) {
        Xs[0] = 1.0;
        for (int s = 1; s < t.nTransported; ++s) Xs[s] = 0.0;
    } else {
        double Yc[THERMO_MAX_SPECIES];
        for (int s = 0; s < nY; ++s) Yc[s] = (Y[s] > 0.0) ? Y[s] : 0.0;
        thermo_X_from_Y(th, nY, Yc, Xs);
    }
    transport_expand_X(t, Xs, Xr);
    transport_mix(t.sp, t.pairs, t.nReal, Xr, T, mu, lam);
    if (Xreal_out) for (int r = 0; r < t.nReal; ++r) Xreal_out[r] = Xr[r];
}

// device 側の輸送表 (thermo_d.cu が所有; thermo_init_db が physProp.transport を解決した DB から上げる)。
//   physProp.transport が無い run では nullptr (カーネルは現行の thermo_mu_mix / thermo_lambda_mix を使う)。
const TransportTableD* thermo_transport_table();
