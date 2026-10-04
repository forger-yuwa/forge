#pragma once

#include <cuda_runtime.h>
#include <math.h>   // host コンパイル時の pow/log/exp/sqrt (device は組込み)
#include "thermo_d.cuh"   // SpeciesThermo / thermo_h_mass (H2O 潜熱の気相 = 種 DB と同じ評価; plan thermophysics-solver-owned-species-db #10)

// 非平衡凝縮: 凝縮種ごとの物性相関 (飽和蒸気圧・凝縮相密度・潜熱・表面張力)。
// 物性評価は exp/log で桁が飛ぶため内部は **double** で計算する (flow_float が float でも安全)。
//
// 【相の方針 (重要)】Lin 2014 のモデルは **液相 (過冷却=supercooled liquid)**。非平衡凝縮は
// 三重点以下でも準安定な過冷却液滴を作るため、本実装も**全域で液相フィットを使う** (固相昇華へは
// 切替えない)。r* = 2σ/(ρ_l RT ln(p_v/p_sat)) も液密度 ρ_l を使う。詳細は methods/condensation.md 8 節。
//
// 【低温外挿ガード (重要)】液相フィット (Jacobsen p_sat, 潜熱多項式) は ~40K 以下に外挿すると破綻する
// (p_sat が ~33K で非単調、潜熱が崩落)。有効域は概ね 45–126K。凝縮 ON では潜熱放出で活性域が ~40K 以上に
// 留まるが、過渡/dry セルが低温に落ちても破綻しないよう、物性評価温度を [T_PROP_FLOOR, Tc) にクランプする。
// クランプは「より低温=より低 p_sat」側に倒れる安全側 (過剰核生成を起こさない)。既知の近似として docs に明記。
//
// 【気相 thermo の制約】凝縮ケースの気相は **thermalMethod 0 (熱量的完全気体, cp/γ 一定)** を使うこと。
// NASA-9/CEA は ~200K 未満で無効 (出口 ~27K)。論文も calorically perfect gas。
//
// 種の切替は enum + 係数構造体 CondSpeciesProps で行う。N2 を最初に実装、H2O は係数差し替えで足す。

// 液相フィットの有効下限 (これ未満は外挿破綻するためクランプ)。
#define COND_T_PROP_FLOOR 45.0

enum CondPropModel {
    COND_MODEL_N2  = 0,
    COND_MODEL_H2O = 1,   // Phase 3
};

// H2O の潜熱モデル = 気液ペア (plans/active/thermophysics-solver-owned-species-db.md §4.8, #10; methods/condensation.md §4.3)。
//   L(T) = h_v(T) − h_l(T)。h_v はペアの気相種そのもの (種 DB の係数・区間・外挿規約、datum オフセット焼き込み後 =
//   thermo_init_db が device に上げる係数と同一)、h_l は共通データの液相 H2O(L) を気相と同じ MW で質量換算し**気相と同じ datum 定数** (R_u Δa7/MW) を足す。
//   起動時に host (condensationTransport_d.cu cond_latent_pair_for) が種 DB の解決結果から cond_latent_pair_make で作り、
//   host 側と device 側に 1 つずつ置いて **ポインタ** (CondLatentRef) で渡す。気相係数のハードコードは持たない。
//   値渡しにしないのは、係数 ~40 double を kernel のローカル構造体に持つとレジスタが溢れるため
//   (実測: dependentVariables_d REG 126→166、512 スレッド/ブロックの上限 128 を超えて起動不能になる)。
struct CondLatentPair {
    int           valid = 0;   // 1: 構築済み。0: 未構築 → H2O の cond_latent は NaN (種 DB を経ない L は作らない)
    SpeciesThermo gas;         // ペアの気相 (datum 焼き込み後)
    double liq[8];             // 液相 NASA-9 a0..a7 (絶対基準のまま; h だけなので a8 は持たない)
    double liqTlo, liqThi;     // 液相フィットの区間 [K] (H2O(L): 273.15–373.15)
    double liqMW;              // 質量換算の MW (= gas.MW; 共通データで一致を検査済み)
    double liqR;               // R_u/liqMW [J/(kg K)]
    double hShift;             // 気相と同じ datum オフセット [J/kg] = R_u Δa7/MW (Tref=0 なら 0)。h_l = h_l,abs + hShift
    double hlLo, cplLo, hlHi;  // 延長規約の前計算 (絶対基準) [J/kg, J/(kg K)]: h_l(Tlo)、T<Tlo の一定 c_p,l (1 K 差分)、h_l(Thi) (T>Thi は頭打ち)
};
// 気液ペアへの参照。host コード (表の生成・起動ログ・単体試験) は h を、kernel は d (device メモリ) を読む。
// どちらかが nullptr のまま H2O の潜熱を評価すると NaN (設定漏れを露呈させる)。
struct CondLatentRef {
    const CondLatentPair* h = nullptr;
    const CondLatentPair* d = nullptr;
};

// 凝縮種の物性パラメータ (device へ value 渡しできる POD)。
struct CondSpeciesProps {
    int    model;   // CondPropModel
    double R;       // 蒸気の気体定数 [J/(kg·K)]
    double cv;      // 蒸気の定積比熱 [J/(kg·K)] (calorically perfect)
    double cp;      // 蒸気の定圧比熱 [J/(kg·K)]
    double Tt;      // 三重点温度 [K] (液/固 切替)
    double Tc;      // 臨界温度 [K]
    double M;       // 分子量 [kg/mol]
    double sigmaScale; // 表面張力の倍率 (感度試験専用 condSigmaScale, 既定 1.0; 核生成・Kelvin・蒸発に一貫)
    // N2 低温物性の切替 (plans/accepted/condensation-air.md §4.2)。H2O では未使用。
    int    latentLowT;   // 1: 70 K 未満の潜熱を c_l 一定の線形外挿 (既定), 0: 旧 4 次多項式 (60 K 未満で L'>0)
    int    psatLowT;     // 1: 50 K 未満の飽和圧 C–C 外挿を新 L(T) の積分で再構成 (既定), 0: 旧 (L_poly(50) 一定の C–C; 診断用)
    double liquidCp;     // 液 N2 の比熱 c_l [J/(kg K)] (線形外挿の傾き c_p,v − c_l に使う; 既定 2000)
    int    gasKgasModel; // 成長則の気相熱伝導率: 0=N2 (n2_kgas), 1=空気 Sutherland (CPG carrier 空気)
    CondLatentRef lat;   // H2O の潜熱 (気液ペアへの参照; #10)。N2 では未使用
};

// kernel に値渡しする物性オプション (config 由来)。condProps_make() で CondSpeciesProps に反映する。
struct CondPropOpts {
    int    latentLowT;
    int    psatLowT;
    double liquidCp;
    int    gasKgasModel;
    double sigmaScale;
    double Yw;           // CPG carrier 形の凝縮種質量分率 (condVaporMassFraction; <=0 で pure)
    CondLatentRef h2oLatent;    // H2O の潜熱モデル (cond_prop_opts が種 DB から設定; 未設定なら H2O の L は NaN)
};

// N2 既定パラメータ。R=296.8, γ=1.4 → cv=R/(γ-1)=742, cp=γcv=1038.8。
__host__ __device__ inline CondSpeciesProps condProps_N2()
{
    CondSpeciesProps s;
    s.model = COND_MODEL_N2;
    s.R  = 296.8;
    s.cv = 742.0;
    s.cp = 1038.8;
    s.Tt = 63.15;
    s.Tc = 126.192;
    s.M  = 0.0280134;
    s.sigmaScale = 1.0;
    s.latentLowT = 1; s.psatLowT = 1; s.liquidCp = 2000.0; s.gasKgasModel = 0;
    return s;
}

// 物性評価温度を液相フィット有効域 [T_PROP_FLOOR, Tc) にクランプ。
__host__ __device__ inline double cond_clamp_Tprop(double T, double Tc)
{
    double Tcl = T;
    if (Tcl < COND_T_PROP_FLOOR) Tcl = COND_T_PROP_FLOOR;
    const double Tmax = Tc - 0.5; // 臨界点直下で発散させない
    if (Tcl > Tmax) Tcl = Tmax;
    return Tcl;
}

// --- N2 潜熱 (蒸発) L(T) [J/kg] --- 式26 (4次多項式 MJ/kg)。psat の C-C 外挿でも使うため前方に置く。
// 【注意】この多項式は 60 K 未満で dL/dT>0 (液比熱 c_l=c_p,v−L' が負) となり熱力学的に不整合 (2026-09-10/12 codex 指摘)。
//   低温整合版は n2_latent_ex(T, lowT=1, c_l) (70 K 未満を c_l 一定の線形外挿, C0 接続)。本関数 (旧) は A/B 用に残す。
__host__ __device__ inline double n2_latent_poly(double T)
{
    const double Tcl = cond_clamp_Tprop(T, 126.192);
    const double p1 = -2.137e-8, p2 = 7.18e-6, p3 = -9.142e-4, p4 = 0.05069, p5 = -0.809;
    const double L = p1*Tcl*Tcl*Tcl*Tcl + p2*Tcl*Tcl*Tcl + p3*Tcl*Tcl + p4*Tcl + p5; // MJ/kg
    const double Lj = L * 1.0e6;
    return (Lj > 0.0) ? Lj : 0.0;
}
__host__ __device__ inline double n2_latent(double T) { return n2_latent_poly(T); }   // 旧 (A/B 用)

#define COND_N2_LATENT_TA 70.0   // 線形外挿の接続温度 [K] (C0 接続; L'(70) は多項式 −1072 vs 線形 c_p,v−c_l で不連続)
#define COND_N2_CPV 1038.8       // N2 蒸気 c_p [J/(kg K)] (CPG)
// 低温整合版: T>=Ta は多項式、T<Ta は L(T)=L(Ta)+(c_p,v−c_l)(T−Ta) (c_l>0 なら L'<0 が保証される)。lowT=0 で旧多項式。
__host__ __device__ inline double n2_latent_ex(double T, int lowT, double cl)
{
    if (!lowT || T >= COND_N2_LATENT_TA) return n2_latent_poly(T);
    return n2_latent_poly(COND_N2_LATENT_TA) + (COND_N2_CPV - cl)*(T - COND_N2_LATENT_TA);
}

// Jacobsen 液飽和圧 (式22, atm→Pa)。有効域内の生評価。
__host__ __device__ inline double n2_psat_jacobsen(double Tcl)
{
    const double Tc = 126.192;
    const double n1 = 8394.409444, n2 = -1890.045259, n3 = -7.282229165;
    const double n4 = 0.01022850966, n5 = 5.556063825e-4, n6 = -5.944544662e-6;
    const double n7 = 2.715433932e-8, n8 = -4.879535904e-11, n9 = 509.5360824;
    const double dTc = Tc - Tcl;
    const double lnP = n1/Tcl + n2 + n3*Tcl + n4*pow(dTc, 1.95)
                     + n5*Tcl*Tcl*Tcl + n6*Tcl*Tcl*Tcl*Tcl + n7*Tcl*Tcl*Tcl*Tcl*Tcl
                     + n8*Tcl*Tcl*Tcl*Tcl*Tcl*Tcl + n9*log(Tcl);
    return exp(lnP) * 101325.0;
}

// --- N2 飽和蒸気圧 p_sat(T) [Pa] (過冷却液) ---
// 有効域 [T_switch, Tc) は Jacobsen (式22)。**T_switch 未満は Clausius–Clapeyron で物理外挿**:
//   p_sat(T) = p_sat(T_sw)·exp(-(L/R)(1/T - 1/T_sw))   (低温で単調に減少、過飽和を正しく与える)。
// 45K クランプ (cond_clamp_Tprop) では psat を凍結し過飽和 S=p/psat を潰してしまい核生成が起きないため、
// psat だけは C-C 外挿を使う (σ,ρ_l,L は緩変化なのでクランプのまま)。L_ref=L(T_sw)。
#define COND_PSAT_TSWITCH 50.0
// psatLowT=1 (既定): 50 K 未満の C–C 外挿を低温整合潜熱 L(T)=L_a+c'(T−T_a) の積分 (閉形式) で再構成:
//   ln(p/p_s) = (1/R)[(L_a − c' T_a)(1/T_s − 1/T) + c' ln(T/T_s)],  c' = c_p,v − c_l。接続 (50 K) は C0 (微分は不連続)。
//   38 K で旧 (L_poly(50)=204 kJ/kg 一定の C–C) の 0.518 倍 (plans/accepted/condensation-air.md §4.2)。
// psatLowT=0: 旧 C–C (診断「潜熱だけ新」の A/B 用)。latentLowT=0 のときは閉形式の L も多項式側では定義できないので旧 C–C に落とす。
__host__ __device__ inline double n2_psat_ex(double T, int psatLowT, int latentLowT, double cl)
{
    const double Tc = 126.192;
    if (T >= COND_PSAT_TSWITCH) {
        double Tcl = (T < Tc - 0.5) ? T : (Tc - 0.5);  // 臨界直下のみ上クランプ
        return n2_psat_jacobsen(Tcl);
    }
    const double Tsw  = COND_PSAT_TSWITCH;
    const double psw  = n2_psat_jacobsen(Tsw);
    const double Rv   = 296.8;
    double Tlo = (T > 5.0) ? T : 5.0;     // 0 割回避
    if (psatLowT && latentLowT) {
        const double Ta = COND_N2_LATENT_TA, La = n2_latent_poly(Ta), cpr = COND_N2_CPV - cl;
        const double lnr = ((La - cpr*Ta)*(1.0/Tsw - 1.0/Tlo) + cpr*log(Tlo/Tsw))/Rv;
        return psw * exp(lnr);
    }
    const double Lref = n2_latent_poly(Tsw);   // 旧: ~2.04e5 J/kg 一定
    return psw * exp(-(Lref/Rv)*(1.0/Tlo - 1.0/Tsw));
}
__host__ __device__ inline double n2_psat(double T) { return n2_psat_ex(T, 0, 0, 2000.0); }   // 旧 (A/B 用)

// --- N2 凝縮相 (液) 密度 ρ_l(T) [kg/m³] --- Nowak (式24)。
__host__ __device__ inline double n2_rho_cond(double T)
{
    const double Tc = 126.192, rhoc = 313.3;
    const double Tcl = cond_clamp_Tprop(T, Tc);
    double tau = 1.0 - Tcl/Tc;
    if (tau < 0.0) tau = 0.0;
    const double n1 = 1.48654237, n2 = -0.280476066, n3 = 0.0894143085, n4 = -0.119879866;
    const double lnr = n1*pow(tau, 0.3294) + n2*pow(tau, 4.0/6.0)
                     + n3*pow(tau, 16.0/6.0) + n4*pow(tau, 35.0/6.0);
    return rhoc * exp(lnr);
}

// --- N2 表面張力 (液) σ(T) [N/m] --- Stansfield (式28, dyn/cm)。
__host__ __device__ inline double n2_sigma(double T)
{
    const double Tc = 126.0, sig0 = 29.06; // dyn/cm
    const double Tcl = cond_clamp_Tprop(T, 126.192);
    double x = 1.0 - Tcl/Tc;
    if (x < 0.0) x = 0.0;
    return sig0 * pow(x, 1.247) * 1.0e-3; // dyn/cm -> N/m
}

// --- N2 気相輸送物性 (成長則の熱伝導 k と 平均自由行程 λ 用) ---
// 成長則 (Goodheart) は気相熱伝導 k(T) と Knudsen 数 Kn=λ/2r を要する。case/34 は非粘性 (viscMethod 0)
// で流れ場の μ/k を持たないため、凝縮ソース内で N2 の簡易モデルから別途評価する。低温外挿のため近似。

// N2 気相粘性 μ(T) [Pa·s] — Sutherland (μ0=1.663e-5 @273K, C=107K)。低温は外挿 (近似)。
__host__ __device__ inline double n2_mu_gas(double T)
{
    const double mu0 = 1.663e-5, T0 = 273.0, C = 107.0;
    double Tc = (T > 5.0) ? T : 5.0;
    return mu0 * (T0 + C)/(Tc + C) * pow(Tc/T0, 1.5);
}

// N2 気相熱伝導 k(T) [W/(m·K)] — k = μ cp/Pr (Pr=0.72, cp=1038.8)。
__host__ __device__ inline double n2_kgas(double T)
{
    const double cp = 1038.8, Pr = 0.72;
    return n2_mu_gas(T) * cp / Pr;
}

// 気相 平均自由行程 λ(T,p) [m] — kinetic theory λ = (μ/p)√(π R T/2)。R_N2=296.8。
__host__ __device__ inline double cond_mean_free_path(double T, double p, double R)
{
    const double pf = (p > 1.0) ? p : 1.0;
    return (n2_mu_gas(T)/pf) * sqrt(3.14159265358979*R*T/2.0);
}

// =====================================================================================
// H2O (水) 物性 — 過冷却液 (Wyslouzil ノズルは ~207K まで膨張、過冷却水滴)。
// psat は Murphy-Koop (2005) の過冷却液式 (123–332K で妥当)、他は標準フィット。
// =====================================================================================

// 水 飽和蒸気圧 p_sat(T) [Pa] — Murphy & Koop (2005) liquid (過冷却水, 123<T<332K)。
__host__ __device__ inline double h2o_psat(double T)
{
    double Tc = (T > 120.0) ? T : 120.0;
    const double lnp = 54.842763 - 6763.22/Tc - 4.210*log(Tc) + 0.000367*Tc
        + tanh(0.0415*(Tc - 218.8)) * (53.878 - 1331.22/Tc - 9.44523*log(Tc) + 0.014025*Tc);
    return exp(lnp); // Pa
}

// 水 液密度 ρ_l(T) [kg/m³] — 過冷却水の簡易フィット (200–300K で ~%)。
__host__ __device__ inline double h2o_rho_cond(double T)
{
    // 過冷却水: 277K で最大 ~1000、低温で僅かに低下。簡易: 線形近似 + フロア。
    double r = 1000.0 - 0.12*(277.0 - T); // ゆるい近似
    if (r < 920.0) r = 920.0;
    return r;
}

// 水 蒸発潜熱 L(T) [J/kg] = h_v(T) − h_l(T) = **気相 H2O と液相 H2O(L) の全エンタルピー差** (CEA の絶対基準のペア)。
//   2026-08-18 のユーザ指示「CEA 式で L を逆算」を、2026-09-27 に種 DB の気液ペアへ移した
//   (plans/active/thermophysics-solver-owned-species-db.md §4.8, #10; 旧実装は気相 H2O の 200–1000 K 係数を再ハードコードし、
//    200 K 未満へ多項式のまま外挿・1000 K で頭打ちしていた。種 DB は 200 K 未満を c_p(200 K) 一定の線形外挿にするので
//    120 K で −2.39 kJ/kg、150 K で −0.47 kJ/kg、200 K 以上は丸め程度で変わる)。
//   h_v: ペアの気相種 = 種 DB と同じ係数・区間・外挿規約・datum (thermo_h_mass そのもの)。
//   h_l: 共通データの H2O(L) 273.15–373.15 K 係数 (CEA; Cox 1989 / Haar 1984) を気相と同じ MW で質量換算し、気相と同じ datum 定数を足す。
//        **273.15 K 未満は CEA に液相フィットが無く** (CEA は氷 H2O(cr))、多項式外挿は 250 K 以下で発散 (cp_l 230 K で 10 kJ/kgK,
//        200 K で 50 kJ/kgK) するので、h_l(273.15) から cp_l ≈ 4228 J/kgK (1 K 差分) 一定で線形外挿する (過冷却水の標準的な扱い)。
//        373.15 K 超は h_l(373.15) で頭打ち (従来と同じ延長規約)。
//   評価温度の下限は COND_T_PROP_FLOOR (45 K; 他の凝縮物性と同じ入力クランプ、NaN もここに落ちる)。L は [1.5, 3.5] MJ/kg にクランプ。
//   L(250 K)=2.563 MJ/kg, L(273.15)=2.501 MJ/kg。氷 (昇華熱 2.835 MJ/kg) は使わない (飽和線も過冷却液 Murphy–Koop で統一)。
//   datum 不変性 (thermoHrefTemp を変えても L が不変) と既知の気液差は tests/unit/test_cond_latent_pair.cu。

// 液相の NASA-9 h [J/kg] (絶対基準)。H2O(L) の係数は a0=1.3e9, a1=−2.4e7 と大きく項どうしが強く打ち消すので、演算順は
// #10 以前の h2o_latent (h2o_nasa9_h_mass) と同じにする (液相部分は旧実装とビット一致; 丸めの再配置で ~5e-4 J/kg 動くのを避ける)。
__host__ __device__ inline double cond_liquid_h_abs_poly(const double* a, double R, double T)
{
    const double hRT = -a[0]/(T*T) + a[1]*log(T)/T + a[2] + a[3]*T/2.0 + a[4]*T*T/3.0
                     + a[5]*T*T*T/4.0 + a[6]*T*T*T*T/5.0 + a[7]/T;
    return hRT*R*T;
}

// 気液ペアを作る (host/device 共用; 起動時に 1 回)。gasAbs・liqAbs は datum 前の絶対基準。Tref>0 なら thermo_init_db と
// 同じ演算で気相の全区間の a7 に Δa7 = −h_abs,gas(Tref)/Ru を足し (device の種 DB 係数とビット一致)、液相には**同じ定数**
// R_u Δa7/MW を h に足す (a7 に Δa7 を足すのと数学的に同じ。液相多項式の打ち消しの丸めを datum から切り離すため、係数は絶対基準のまま)。
__host__ __device__ inline CondLatentPair cond_latent_pair_make(const SpeciesThermo& gasAbs, const double* liqAbs,
                                                                double liqTlo, double liqThi, double Tref)
{
    CondLatentPair p;
    p.valid = 1;
    p.gas = gasAbs;
    p.gas.invMW = 1.0/p.gas.MW;
    double da7 = 0.0;
    if (Tref > 0.0) {
        const double h_ref = thermo_h_molar(gasAbs, Tref);   // 移動前の絶対 h [J/mol] (thermo_init_db と同じ)
        da7 = -h_ref / THERMO_RU;
        thermo_add_a7(p.gas, da7);   // 全区間 (thermo_init_db と同じ)
        p.gas.h_datum  = h_ref;
    }
    for (int k = 0; k < 8; ++k) p.liq[k] = liqAbs[k];
    p.liqTlo = liqTlo; p.liqThi = liqThi; p.liqMW = gasAbs.MW;
    p.liqR   = THERMO_RU/p.liqMW;
    p.hShift = THERMO_RU*da7/p.liqMW;
    p.hlLo  = cond_liquid_h_abs_poly(p.liq, p.liqR, liqTlo);
    p.cplLo = (cond_liquid_h_abs_poly(p.liq, p.liqR, liqTlo + 0.5) - cond_liquid_h_abs_poly(p.liq, p.liqR, liqTlo - 0.5 + 1.0e-9))/1.0;   // ≈4228 J/kgK (従来と同じ 1 K 差分)
    p.hlHi  = cond_liquid_h_abs_poly(p.liq, p.liqR, liqThi);
    return p;
}

// 気相 h_v(T) [J/kg]: thermo_h_mass (thermo_d.cuh) と同じ式・同じ分岐 (区間は thermo_interval と同じ規約
//   (区切りちょうどは上の区間)、[Tlo,Thi] の外は端の c_p で線形外挿)。
//   thermo_pick_coeffs のポインタ選択を kernel 内のローカル構造体に使うとローカルメモリへ落ちる (dependentVariables_d で
//   REG 126→178・STACK 360→1608 を実測) ので、係数を値で選ぶ (区間番号 k の 3 択。THERMO_MAX_INTERVALS=3 に合わせる)。
//   thermo_h_mass とのビット一致は tests/unit/test_cond_latent_pair.cu。
#if THERMO_MAX_INTERVALS != 3
#error "cond_gas_coef assumes THERMO_MAX_INTERVALS == 3 (select by value)"
#endif
__host__ __device__ inline double cond_gas_coef(const SpeciesThermo& sp, int k, int i)
{
    return (k == 0) ? sp.coef[0][i] : ((k == 1) ? sp.coef[1][i] : sp.coef[2][i]);
}
__host__ __device__ inline double cond_gas_h_molar_clamped(const SpeciesThermo& sp, double Tc)
{
    const int k = thermo_interval(sp, Tc);
    const double a0 = cond_gas_coef(sp, k, 0), a1 = cond_gas_coef(sp, k, 1);
    const double a2 = cond_gas_coef(sp, k, 2), a3 = cond_gas_coef(sp, k, 3);
    const double a4 = cond_gas_coef(sp, k, 4), a5 = cond_gas_coef(sp, k, 5);
    const double a6 = cond_gas_coef(sp, k, 6), a7 = cond_gas_coef(sp, k, 7);
    const double Ti  = 1.0/Tc;
    const double Ti2 = Ti*Ti;
    const double lnT = log(Tc);
    const double hRT = -a0*Ti2 + a1*lnT*Ti + a2
                     + a3*Tc/2.0 + a4*Tc*Tc/3.0 + a5*Tc*Tc*Tc/4.0
                     + a6*Tc*Tc*Tc*Tc/5.0 + a7*Ti;
    return THERMO_RU * Tc * hRT;
}
__host__ __device__ inline double cond_gas_cp_molar_clamped(const SpeciesThermo& sp, double Tc)
{
    const int k = thermo_interval(sp, Tc);
    const double a0 = cond_gas_coef(sp, k, 0), a1 = cond_gas_coef(sp, k, 1);
    const double a2 = cond_gas_coef(sp, k, 2), a3 = cond_gas_coef(sp, k, 3);
    const double a4 = cond_gas_coef(sp, k, 4), a5 = cond_gas_coef(sp, k, 5);
    const double a6 = cond_gas_coef(sp, k, 6);
    const double Ti  = 1.0/Tc;
    const double Ti2 = Ti*Ti;
    return THERMO_RU * ( a0*Ti2 + a1*Ti + a2
                       + a3*Tc + a4*Tc*Tc + a5*Tc*Tc*Tc + a6*Tc*Tc*Tc*Tc );
}
__host__ __device__ inline double cond_gas_h_mass(const SpeciesThermo& sp, double T)
{
    double h;
    if (T < sp.Tlo)      h = cond_gas_h_molar_clamped(sp, sp.Tlo) + cond_gas_cp_molar_clamped(sp, sp.Tlo)*(T - sp.Tlo);
    else if (T > sp.Thi) h = cond_gas_h_molar_clamped(sp, sp.Thi) + cond_gas_cp_molar_clamped(sp, sp.Thi)*(T - sp.Thi);
    else                 h = cond_gas_h_molar_clamped(sp, T);
    return h / sp.MW;
}

// 液相 h_l(T) [J/kg] (延長規約込み、気相と同じ datum オフセット hShift を足す)
__host__ __device__ inline double cond_liquid_h_mass(const CondLatentPair& p, double T)
{
    double h;
    if (T < p.liqTlo)      h = p.hlLo - p.cplLo*(p.liqTlo - T);
    else if (T > p.liqThi) h = p.hlHi;
    else                   h = cond_liquid_h_abs_poly(p.liq, p.liqR, T);
    return h + p.hShift;
}

__host__ __device__ inline double h2o_latent_pair(const CondLatentPair& p, double T)
{
    if (!p.valid) return NAN;   // 種 DB を経ない潜熱は作らない (起動時の設定漏れを NaN で露呈させる)
    const double Tg = (T > COND_T_PROP_FLOOR) ? T : COND_T_PROP_FLOOR;
    const double hv = cond_gas_h_mass(p.gas, Tg);   // 種 DB と同じ評価 (= thermo_h_mass; 200 K 未満は c_p(200 K) 一定の線形外挿)
    double L = hv - cond_liquid_h_mass(p, Tg);
    if (L < 1.5e6) L = 1.5e6;
    if (L > 3.5e6) L = 3.5e6;
    return L;
}
// 参照から評価 (host は r.h、device は r.d)。
__host__ __device__ inline double h2o_latent(const CondLatentRef& r, double T)
{
#ifdef __CUDA_ARCH__
    const CondLatentPair* p = r.d;
#else
    const CondLatentPair* p = r.h;
#endif
    return (p != nullptr) ? h2o_latent_pair(*p, T) : NAN;
}

// 水 表面張力 σ(T) [N/m] — IAPWS 形を過冷却へ外挿 (Tc=647.096K)。
__host__ __device__ inline double h2o_sigma(double T)
{
    const double Tc = 647.096;
    double tau = (Tc - T)/Tc;
    if (tau < 0.0) tau = 0.0;
    double s = 0.2358 * pow(tau, 1.256) * (1.0 - 0.625*tau); // N/m
    if (s < 0.0) s = 0.0;
    return s;
}

// --- 種ディスパッチ (model で N2 / H2O を切替) ---
__host__ __device__ inline double cond_psat(const CondSpeciesProps& s, double T)
{
    return (s.model == COND_MODEL_H2O) ? h2o_psat(T) : n2_psat_ex(T, s.psatLowT, s.latentLowT, s.liquidCp);
}
// 空気 (CPG carrier) の気相熱伝導率 [W/(m K)]: Sutherland μ (μ0=1.716e-5 @273 K, C=111) × c_p/Pr (c_p 1008.7, Pr 0.72)。低温外挿 (近似)。
__host__ __device__ inline double air_kgas(double T)
{
    const double mu0 = 1.716e-5, T0 = 273.0, C = 111.0, cp = 1008.7, Pr = 0.72;
    double Tc = (T > 5.0) ? T : 5.0;
    return mu0 * (T0 + C)/(Tc + C) * pow(Tc/T0, 1.5) * cp/Pr;
}
// 成長則・蒸発・二温度が使う気相 (キャリア) 熱伝導率のディスパッチ (gasKgasModel: 0=N2, 1=空気)。
__host__ __device__ inline double cond_kgas(const CondSpeciesProps& s, double T)
{
    return (s.gasKgasModel == 1) ? air_kgas(T) : n2_kgas(T);
}
__host__ __device__ inline double cond_rho_cond(const CondSpeciesProps& s, double T)
{
    return (s.model == COND_MODEL_H2O) ? h2o_rho_cond(T) : n2_rho_cond(T);
}
__host__ __device__ inline double cond_latent(const CondSpeciesProps& s, double T)
{
    return (s.model == COND_MODEL_H2O) ? h2o_latent(s.lat, T) : n2_latent_ex(T, s.latentLowT, s.liquidCp);
}
__host__ __device__ inline double cond_sigma(const CondSpeciesProps& s, double T)
{
    // sigmaScale は感度試験専用の一定倍率 (既定 1.0: ×1.0 は IEEE 恒等でビット不変)。
    return s.sigmaScale * ((s.model == COND_MODEL_H2O) ? h2o_sigma(T) : n2_sigma(T));
}

// 飽和温度 T_sat(p_v): ln p_sat(T) = ln p_v を Clausius-Clapeyron 勾配 d ln p/dT = L/(R T^2) の Newton で反転。
// 過冷却度 T_sat - T の診断用 (plans/accepted/condensation-equilibrium.md)。p_v <= 1e-6 Pa は 0 を返す。
__host__ __device__ inline double cond_Tsat(const CondSpeciesProps& s, double pv, double T_guess)
{
    if (!(pv > 1.0e-6)) return 0.0;
    const double lnpv = log(pv);
    double T = (T_guess > 50.0 && T_guess < 1000.0) ? T_guess : 250.0;
    #pragma unroll 1
    for (int it = 0; it < 25; ++it) {
        const double ps = cond_psat(s, T);
        if (!(ps > 1.0e-300)) { T *= 1.2; continue; }
        const double f  = log(ps) - lnpv;
        double dfdT = cond_latent(s, T)/(s.R*T*T);
        if (dfdT < 1.0e-6) dfdT = 1.0e-6;
        double dT = f/dfdT;
        if (dT >  0.3*T) dT =  0.3*T;
        if (dT < -0.3*T) dT = -0.3*T;
        T -= dT;
        if (T < 50.0) T = 50.0;
        if (T > 1000.0) T = 1000.0;
        if (fabs(dT) < 1.0e-4) break;
    }
    return T;
}

// H2O 既定パラメータ (キャリア+凝縮種)。R_H2O=461.5, M=0.0180153。
__host__ __device__ inline CondSpeciesProps condProps_H2O()
{
    CondSpeciesProps s;
    s.model = COND_MODEL_H2O;
    s.R  = 461.5;
    // 蒸気 cp/cv: NASA-9 (species_db.yaml H2O) の 200–300 K 平均 (cp 1851–1865, γ_v 1.329–1.332)。
    // Kantrowitz 非等温補正の γ_v=cp/cv=1.331 に使う (純蒸気形の係数 2(γ_v−1)/(γ_v+1)=R_v/(c_v,v+R_v/2))。
    // 気相 EOS は TP では NASA-9 を使うのでここは核生成補正専用 (plans/active/condensation-kantrowitz-gamma-twophase-sonic.md)。
    s.cv = 1393.5;
    s.cp = 1855.0;
    s.Tt = 273.16;
    s.Tc = 647.096;
    s.M  = 0.0180153;
    s.sigmaScale = 1.0;
    s.latentLowT = 1; s.psatLowT = 1; s.liquidCp = 2000.0; s.gasKgasModel = 0;   // H2O では未使用
    // 潜熱の気液ペアへの参照 s.lat は未設定 (nullptr → L は NaN)。種 DB から作ったもの (CondPropOpts::h2oLatent) を condProps_make で入れる。
    return s;
}

// config 由来のオプションを反映した凝縮種物性 (kernel 内で構築する場合はこれを使う)。
__host__ __device__ inline CondSpeciesProps condProps_make(int model, const CondPropOpts& o)
{
    CondSpeciesProps s = (model == COND_MODEL_H2O) ? condProps_H2O() : condProps_N2();
    s.sigmaScale = o.sigmaScale; s.latentLowT = o.latentLowT; s.psatLowT = o.psatLowT; s.liquidCp = o.liquidCp; s.gasKgasModel = o.gasKgasModel;
    if (model == COND_MODEL_H2O) s.lat = o.h2oLatent;   // 潜熱の気液ペア (#10)
    return s;
}
