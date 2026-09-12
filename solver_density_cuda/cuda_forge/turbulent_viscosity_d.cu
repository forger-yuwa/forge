#include "calcGradient_d.cuh"
#include <cstdlib>
#include <cstdio>
#include "cuda_forge/cudaWrapper.cuh"


__global__ void WALE_d
( 
 // mesh structure
 geom_int nCells,
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro  ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* T   ,
 flow_float* Ht  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* vis_turb , flow_float* wall_dist
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells ) {

        flow_float karman = 0.41;
        //flow_float Cw = 0.5;
        flow_float Cw = 0.325;

        geom_float volume = vol[ic];


        flow_float g11 = dUxdx[ic];
        flow_float g12 = dUxdy[ic];
        flow_float g13 = dUxdz[ic];
        flow_float g21 = dUydx[ic];
        flow_float g22 = dUydy[ic];
        flow_float g23 = dUydz[ic];
        flow_float g31 = dUzdx[ic];
        flow_float g32 = dUzdy[ic];
        flow_float g33 = dUzdz[ic];

        // Sd_ij = ½(ḡ²_ij + ḡ²_ji) − ⅓δ_ij ḡ²_kk, ḡ²_ij = g_ik g_kj (**行列 2 乗**; Nicoud-Ducros 1999)。
        // 旧実装は成分ごとの 2 乗 (½(g12²+g21²) 等) で本来の WALE 演算子ではなかった
        // (純せん断で Sd=0 になるべきところ偽値を作る)。式は tools/ 相当の python 検証済
        // (plans/active/turbulence-wale-fix.md)。
        flow_float q11 = g11*g11 + g12*g21 + g13*g31;
        flow_float q12 = g11*g12 + g12*g22 + g13*g32;
        flow_float q13 = g11*g13 + g12*g23 + g13*g33;
        flow_float q21 = g21*g11 + g22*g21 + g23*g31;
        flow_float q22 = g21*g12 + g22*g22 + g23*g32;
        flow_float q23 = g21*g13 + g22*g23 + g23*g33;
        flow_float q31 = g31*g11 + g32*g21 + g33*g31;
        flow_float q32 = g31*g12 + g32*g22 + g33*g32;
        flow_float q33 = g31*g13 + g32*g23 + g33*g33;
        flow_float qtr3 = (q11 + q22 + q33)/3.0;
        flow_float Sd11 = q11 - qtr3;
        flow_float Sd12 = 0.5*(q12+q21);
        flow_float Sd13 = 0.5*(q13+q31);
        flow_float Sd21 = Sd12;
        flow_float Sd22 = q22 - qtr3;
        flow_float Sd23 = 0.5*(q23+q32);
        flow_float Sd31 = Sd13;
        flow_float Sd32 = Sd23;
        flow_float Sd33 = q33 - qtr3;


        flow_float S11 = 0.5*(g11+g11);
        flow_float S12 = 0.5*(g12+g21);
        flow_float S13 = 0.5*(g13+g31);
        flow_float S21 = 0.5*(g21+g12);
        flow_float S22 = 0.5*(g22+g22);
        flow_float S23 = 0.5*(g23+g32);
        flow_float S31 = 0.5*(g31+g13);
        flow_float S32 = 0.5*(g32+g23);
        flow_float S33 = 0.5*(g33+g33);

        flow_float SdijSdij = Sd11*Sd11 + Sd12*Sd12 + Sd13*Sd13
                            + Sd21*Sd21 + Sd22*Sd22 + Sd23*Sd23
                            + Sd31*Sd31 + Sd32*Sd32 + Sd33*Sd33;

        flow_float SijSij = S11*S11 + S12*S12 + S13*S13
                          + S21*S21 + S22*S22 + S23*S23
                          + S31*S31 + S32*S32 + S33*S33;

        geom_float d = wall_dist[ic];
        flow_float Ls = min(karman*d , Cw*pow(vol[ic],1.0/3.0));

        // WALE は一様流（勾配ゼロ）で分子分母とも 0 になり 0/0=NaN を生む。標準的に分母へ
        // 小さな正則化項を加える（無勾配では vis_turb=0 となり物理的に正しい）。
        const flow_float wale_denom = pow(SijSij, 5.0/2.0) + pow(SdijSdij, 5.0/4.0)
                                    + static_cast<flow_float>(1.0e-30);
        vis_turb[ic] =  ro[ic]*Ls*Ls*pow(SdijSdij, 3.0/2.0)/wale_denom;
    }
}


// σ-model (Nicoud, Toda, Cabrit, Bose, Lee, PoF 2011): 静的 SGS の改良版。
//   ν_t = (C_σ Δ)² D_σ,  D_σ = σ3(σ1−σ2)(σ2−σ3)/σ1²,  Δ = V^{1/3},  C_σ = 1.35。
//   σ1≥σ2≥σ3 は速度勾配 g の特異値 (= G=gᵀg の固有値の平方根、SPD 3×3 の三角関数閉形式)。
//   設計性質: 2成分流 (TG 初期場含む)・純せん断・剛体回転・等方/軸対称伸長で **厳密に 0**
//   → WALE の遷移前倒し (層流勾配で偽 ν_t) への構造的耐性。壁減衰も D_σ 自身が担うため
//   wall_dist 不要。式は tools/verify_sigma_model.py で SVD と突き合わせ検証済
//   (plans/active/turbulence-sigma-model.md)。固有値は double で評価 (縮退近傍の丸め対策)。
__global__ void SIGMA_d
(
 geom_int nCells,
 geom_float* vol,
 flow_float* ro,
 flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
 flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
 flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
 flow_float* vis_turb
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const double Csig = 1.35;
        const double g11=dUxdx[ic], g12=dUxdy[ic], g13=dUxdz[ic];
        const double g21=dUydx[ic], g22=dUydy[ic], g23=dUydz[ic];
        const double g31=dUzdx[ic], g32=dUzdy[ic], g33=dUzdz[ic];
        // G = gᵀ g (SPD)
        const double G11 = g11*g11 + g21*g21 + g31*g31;
        const double G12 = g11*g12 + g21*g22 + g31*g32;
        const double G13 = g11*g13 + g21*g23 + g31*g33;
        const double G22 = g12*g12 + g22*g22 + g32*g32;
        const double G23 = g12*g13 + g22*g23 + g32*g33;
        const double G33 = g13*g13 + g23*g23 + g33*g33;
        const double I1 = G11 + G22 + G33;
        const double I2 = G11*G22 + G11*G33 + G22*G33 - (G12*G12 + G13*G13 + G23*G23);
        const double I3 = G11*(G22*G33 - G23*G23) - G12*(G12*G33 - G23*G13) + G13*(G12*G23 - G22*G13);
        double a1 = I1*I1/9.0 - I2/3.0;      a1 = fmax(a1, 0.0);
        double a2 = I1*I1*I1/27.0 - I1*I2/6.0 + I3/2.0;
        const double s  = sqrt(a1);
        const double dn = a1*s;              // a1^{3/2}
        double arg = (dn > 0.0) ? a2/dn : 0.0;
        arg = fmin(1.0, fmax(-1.0, arg));
        const double a3 = acos(arg)/3.0;
        const double PI3 = 1.0471975511965976; // π/3
        double l1 = I1/3.0 + 2.0*s*cos(a3);
        double l2 = I1/3.0 - 2.0*s*cos(PI3 + a3);
        double l3 = I1/3.0 - 2.0*s*cos(PI3 - a3);
        // 閉形式は l1 >= l3 >= l2 の順で返る (検証済) が、丸め対策で明示ソート
        double t;
        if (l1 < l2) { t=l1; l1=l2; l2=t; }
        if (l1 < l3) { t=l1; l1=l3; l3=t; }
        if (l2 < l3) { t=l2; l2=l3; l3=t; }
        const double s1 = sqrt(fmax(l1, 0.0));
        const double s2 = sqrt(fmax(l2, 0.0));
        const double s3 = sqrt(fmax(l3, 0.0));
        const double Dsig = (s1 > 0.0) ? s3*(s1-s2)*(s2-s3)/(s1*s1) : 0.0;
        const double delta = cbrt((double)vol[ic]);
        vis_turb[ic] = (flow_float)((double)ro[ic] * (Csig*delta)*(Csig*delta) * fmax(Dsig, 0.0));
    }
}


// SST 渦粘性: mu_t = rho * a1 * k / max(a1 * omega, S * F2)
// 軸対称 (isAxisymmetric=1) では planar 速度勾配に現れないフープひずみ
// S_thetatheta = u_r/r (= axisym_uy_over_r) を S^2 に加える (methods/turbulence/theory.md §7.2)。
__global__ void sst_eddy_viscosity_d(
    geom_int nCells,
    flow_float* ro,
    flow_float* k,
    flow_float* omega,
    flow_float* vis_lam,
    flow_float* wall_dist,
    flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
    flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
    flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
    int isAxisymmetric,
    flow_float* axisym_uy_over_r,
    flow_float* vis_turb,
    // node ω Dirichlet (nodeOmegaWfDirichlet, 第一内層ノード roOmega_wf>=0): Bradshaw/strain
    // リミッタを迂回し νt=ρk/ω (=構成上 ρν_t,wall の mixing-length 値) を使う。ジャンプの巨大 S が
    // リミッタ経由で νt を頭打ちにし「S 大→νt 小→混合無→ジャンプ維持」の正帰還ロックを作るため
    // (case/38 run_0011/0012 で実測・log 則張替でも再形成 = 構造不安定)。nullptr で無効。
    flow_float* roOmega_wf,
    // 2×2 分離実験 (plan §3): ω ピンと strain limiter 迂回を独立に切る。
    //   wf_irep_flag : 第一内層ノードのマスク (bypass を pin と独立に当てるため)
    //   bypassMode   : -1=従来 (roOmega_wf>=0 に追随) / 0=強制 OFF / 1=強制 ON (irep のみ)
    flow_float* wf_irep_flag, int bypassMode)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    constexpr flow_float a1       = static_cast<flow_float>(0.31);
    constexpr flow_float betaStar = static_cast<flow_float>(0.09);
    constexpr flow_float kSmall   = static_cast<flow_float>(1.0e-12);

    const flow_float rho    = max(ro[ic], kSmall);
    const flow_float k_c    = max(k[ic],  static_cast<flow_float>(0.0));
    const flow_float w_c    = max(omega[ic], kSmall);
    const flow_float mu_lam = vis_lam[ic];
    const flow_float y      = max(wall_dist[ic], kSmall);
    const flow_float nu     = mu_lam / rho;

    // ひずみ速度の大きさ S = sqrt(2 Sij Sij)
    const flow_float S11 = dUxdx[ic];
    const flow_float S22 = dUydy[ic];
    const flow_float S33 = dUzdz[ic];
    const flow_float S12 = static_cast<flow_float>(0.5) * (dUxdy[ic] + dUydx[ic]);
    const flow_float S13 = static_cast<flow_float>(0.5) * (dUxdz[ic] + dUzdx[ic]);
    const flow_float S23 = static_cast<flow_float>(0.5) * (dUydz[ic] + dUzdy[ic]);
    flow_float S_sq = static_cast<flow_float>(2.0) *
        (S11*S11 + S22*S22 + S33*S33 +
         static_cast<flow_float>(2.0) * (S12*S12 + S13*S13 + S23*S23));
    // 軸対称: フープひずみ S_thetatheta = u_r/r を加算 (2 * S_tt^2)
    if (isAxisymmetric) {
        const flow_float S_tt = axisym_uy_over_r[ic];
        S_sq += static_cast<flow_float>(2.0) * S_tt * S_tt;
    }
    const flow_float S_mag = sqrt(max(S_sq, static_cast<flow_float>(0.0)));

    // F2 ブレンド関数
    const flow_float arg2_a = static_cast<flow_float>(2.0) * sqrt(k_c) / (betaStar * w_c * y);
    const flow_float arg2_b = static_cast<flow_float>(500.0) * nu / (w_c * y * y);
    const flow_float arg2   = max(arg2_a, arg2_b);
    const flow_float F2     = tanh(arg2 * arg2);

    const bool wfPin = (roOmega_wf != nullptr && roOmega_wf[ic] >= static_cast<flow_float>(0.0));
    // limiter 迂回の可否を ω ピンから切り離す (既定 bypassMode=-1 は従来と完全同一)。
    bool bypass;
    if (bypassMode < 0)       bypass = wfPin;                       // 従来
    else if (bypassMode == 0) bypass = false;                       // 強制 OFF (E1: pin のみ)
    else                      bypass = (wf_irep_flag != nullptr
                                        && wf_irep_flag[ic] > static_cast<flow_float>(0.5)); // 強制 ON (E2)
    vis_turb[ic] = bypass ? rho * k_c / w_c
                          : rho * a1 * k_c / max(a1 * w_c, S_mag * F2);
}


// ===========================================================================
// SST-DDES length scale (methods/turbulence §8, plan turbulence-iddes-sst.md §4.4)
// ---------------------------------------------------------------------------
// Spalart et al. (2006) シールド関数 f_d で剥離域のみ LES、付着 BL は RANS に保つ。
//   r_d = (ν_t + ν) / (κ² d² √(∂Ui/∂xj·∂Ui/∂xj))      κ=0.41
//   f_d = 1 - tanh((8 r_d)³)
//   l_RANS = √k / (β*^{1/4} ω),  l_LES = C_DES·Δ
//   l_DDES = l_RANS - f_d·max(0, l_RANS - C_DES·Δ)
//
// ★ float32 ガード (plan §4.4): r_d は d² で効くため近壁・よどみ・一様流で飽和しやすい。
//   - 勾配大きさ / 分母に floor を入れ 0 割を防ぐ (一様流 ∇u→0 では r_d→大→f_d→0=RANS が正しい極限)。
//   - r_d を [0, 10] に clamp し f_d が 0/1 へ過度に張り付くのを防ぐ。
//   - r_d / f_d を必ず出力 (NaN が無くても全域張り付きは DES 不全の兆候・l_des だけでは見抜けない)。
//
// 依存順 (plan §4.7・最重要): r_d は ν_t=vis_turb/ρ を読むため、必ず sst_eddy_viscosity_d
// (vis_turb 計算) の後に呼ぶ。本カーネルは CV 抽象 (volume→Δ, wall_dist→d, ∇u, vis_turb) のみに
// 依存し cfg.discretization を参照しないので、node-centered (median-dual) でも無変更で成立する (§5.6)。
__device__ static flow_float compute_rd_ddes(
    flow_float nu_t, flow_float nu, flow_float y, flow_float grad2)
{
    constexpr flow_float kKarman    = static_cast<flow_float>(0.41);
    constexpr flow_float kGradFloor = static_cast<flow_float>(1.0e-30);
    constexpr flow_float kDenomFloor= static_cast<flow_float>(1.0e-30);
    const flow_float gradMag = sqrt(max(grad2, kGradFloor));
    flow_float denom = kKarman * kKarman * y * y * gradMag;
    denom = max(denom, kDenomFloor);
    const flow_float rd = (nu_t + nu) / denom;
    return min(rd, static_cast<flow_float>(10.0));
}

__device__ static flow_float compute_fd_ddes(flow_float rd)
{
    const flow_float arg = static_cast<flow_float>(8.0) * rd;
    return static_cast<flow_float>(1.0) - tanh(arg * arg * arg);
}

// ===========================================================================
// SST-IDDES (DESmode==2, methods/turbulence §8.6, plan §3.4/§4.6)
// ---------------------------------------------------------------------------
// Shur et al. (2008) の IDDES を SST 用に定数調整した Gritskevich et al. (2012)。
// DDES に WMLES 分岐を追加: 解像乱流があれば (瞬時 νt が SGS 級 → r_dt≪1 → f_dt→1)
// f̃_d→f_B となり壁直近のみ RANS・log 層から LES。解像乱流が無ければ DDES 同等に縮退。
//   Δ_IDDES = min(C_w·max(d_w, h_max), h_max)   (h_wn 不要の Gritskevich 簡略形;
//             wall_dist と delta_les だけで組める。壁で 0.15h_max → 遠方で h_max)
//   l_IDDES = f̃_d(1+f_e)·l_RANS + (1−f̃_d)·C_DES·Δ_IDDES
// 定数 (SST 較正): C_t=1.87, C_l=5.0, C_dt1=20 (DDES f_d の 8 と異なる), C_w=0.15。
// 関数形状・RANS 縮退・WMLES 分岐・float32 ガードは tools/verify_iddes_functions.py で
// 検証済 (VERDICT: PASS)。DDES (DESmode==1) 経路はビット不変で共存する。
__device__ static flow_float compute_l_iddes(
    flow_float l_rans, flow_float C_DES, flow_float hmax,
    flow_float y, flow_float rdt, flow_float rdl,
    flow_float* fdtil_out, flow_float* fe_out)
{
    constexpr flow_float kCt   = static_cast<flow_float>(1.87);
    constexpr flow_float kCl   = static_cast<flow_float>(5.0);
    constexpr flow_float kCdt1 = static_cast<flow_float>(20.0);
    constexpr flow_float kCw   = static_cast<flow_float>(0.15);
    constexpr flow_float one   = static_cast<flow_float>(1.0);

    const flow_float h  = max(hmax, static_cast<flow_float>(1.0e-20));
    const flow_float al = static_cast<flow_float>(0.25) - y / h;
    const flow_float al2 = al * al;

    // 近壁ブレンド f_B と log-layer mismatch 補正 f_e (付着 log 層 r_dt≈1 では f_t→1 で自動停止)
    const flow_float fb  = min(static_cast<flow_float>(2.0)
                               * exp(static_cast<flow_float>(-9.0) * al2), one);
    const flow_float fe1 = (al >= static_cast<flow_float>(0.0))
        ? static_cast<flow_float>(2.0) * exp(static_cast<flow_float>(-11.09) * al2)
        : static_cast<flow_float>(2.0) * exp(static_cast<flow_float>(-9.0)  * al2);
    const flow_float argt = kCt * kCt * rdt;
    const flow_float ft   = tanh(argt * argt * argt);
    const flow_float argl  = kCl * kCl * rdl;                 // (C_l² r_dl)^10 を乗算で構成
    const flow_float argl2 = argl * argl;
    const flow_float argl4 = argl2 * argl2;
    const flow_float fl   = tanh(argl4 * argl4 * argl2);
    const flow_float fe2  = one - max(ft, fl);
    const flow_float fe   = fe2 * max(fe1 - one, static_cast<flow_float>(0.0));

    // WMLES/DDES 自動切替 f̃_d
    const flow_float argdt = kCdt1 * rdt;
    const flow_float fdt   = one - tanh(argdt * argdt * argdt);
    const flow_float fdtil = max(one - fdt, fb);

    // 壁距離依存グリッドスケール Δ_IDDES と凸ブレンド
    const flow_float delta_i = min(kCw * max(y, hmax), hmax);
    const flow_float l_les   = C_DES * delta_i;
    *fdtil_out = fdtil;
    *fe_out    = fe;
    return fdtil * (one + fe) * l_rans + (one - fdtil) * l_les;
}

__global__ void compute_des_length_d(
    geom_int nCells,
    int DESmode,
    flow_float* ro,
    flow_float* k,
    flow_float* omega,
    flow_float* vis_lam,
    flow_float* vis_turb,
    flow_float* wall_dist,
    flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
    flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
    flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
    flow_float* delta_les,
    flow_float C_DES,
    flow_float* l_des,
    flow_float* fd_shield,
    flow_float* rd_des,
    flow_float* fe_iddes)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    constexpr flow_float kSmall    = static_cast<flow_float>(1.0e-12);
    constexpr flow_float kBetaStar = static_cast<flow_float>(0.09);

    const flow_float rho  = max(ro[ic], kSmall);
    const flow_float k_c  = max(k[ic],  static_cast<flow_float>(0.0));
    const flow_float w_c  = max(omega[ic], kSmall);
    const flow_float nu   = vis_lam[ic] / rho;
    const flow_float nu_t = max(vis_turb[ic], static_cast<flow_float>(0.0)) / rho;
    const flow_float y    = max(wall_dist[ic], kSmall);

    // |∂Ui/∂xj|² = 全 9 成分の二乗和 (= S²+Ω²)。軸対称項は含めない (検証は 3D 直交・plan §4.7)。
    const flow_float g11 = dUxdx[ic], g12 = dUxdy[ic], g13 = dUxdz[ic];
    const flow_float g21 = dUydx[ic], g22 = dUydy[ic], g23 = dUydz[ic];
    const flow_float g31 = dUzdx[ic], g32 = dUzdy[ic], g33 = dUzdz[ic];
    const flow_float grad2 = g11*g11 + g12*g12 + g13*g13
                           + g21*g21 + g22*g22 + g23*g23
                           + g31*g31 + g32*g32 + g33*g33;

    // l_RANS = √k/(β* ω)。これは SST k 消滅項 D_k = β*ρkω = ρk^{3/2}/l_RANS と厳密に整合する
    // 唯一の定義 (Strelets 2001 / Gritskevich 2012)。f_d=0 の完全シールド域で l_DDES=l_RANS となり
    // D_k が標準 SST に一致する (これを満たさないと mode1 はシールド域でも標準 SST に縮退せず MSD を起こす)。
    const flow_float l_rans = sqrt(k_c) / (kBetaStar * w_c);

    if (DESmode == 2) {
        // IDDES: r_dt (νt 版)・r_dl (ν 版) は compute_rd_ddes を第 2 引数 0 で流用
        // (floor + clamp[0,10] 込み。f_t/f_l/f_dt の tanh は clamp 上限で完全飽和するので影響なし)。
        const flow_float rdt = compute_rd_ddes(nu_t, static_cast<flow_float>(0.0), y, grad2);
        const flow_float rdl = compute_rd_ddes(nu,   static_cast<flow_float>(0.0), y, grad2);
        flow_float fdtil, fe;
        const flow_float l_i = compute_l_iddes(l_rans, C_DES, delta_les[ic], y, rdt, rdl,
                                               &fdtil, &fe);
        l_des[ic]     = max(l_i, static_cast<flow_float>(1.0e-20));
        fd_shield[ic] = fdtil;   // DESmode==2 では f̃_d (1=RANS...0=LES の向きは f_d と逆なので注意)
        rd_des[ic]    = rdt;
        fe_iddes[ic]  = fe;
        return;
    }

    // DDES (DESmode==1・従来経路ビット不変)
    const flow_float rd = compute_rd_ddes(nu_t, nu, y, grad2);
    const flow_float fd = compute_fd_ddes(rd);
    const flow_float l_les  = C_DES * delta_les[ic];
    const flow_float l_ddes = l_rans - fd * max(static_cast<flow_float>(0.0), l_rans - l_les);

    l_des[ic]     = max(l_ddes, static_cast<flow_float>(1.0e-20));
    fd_shield[ic] = fd;
    rd_des[ic]    = rd;
    fe_iddes[ic]  = static_cast<flow_float>(0.0);
}


void turbulent_viscosity_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    // node SST 壁関数が有効か (2×2 分離ゲートと wf_irep_flag の受け渡しに使う)
    const bool sstWfActive = (cfg.discretization == "node" && cfg.LESorRANS == 2
                              && cfg.RANSmodel == 1 && cfg.wallTreatmentSST == 1);

    if (cfg.LESorRANS == 1) { // LES

        if (cfg.LESmodel == 1) {
            // node-centered (median-dual) でも成立: WALE_d は DOF (cell/node) ごとに勾配・体積・wall_dist を
            // 読むだけで cfg.discretization に依存しないため、SST と同じ normalcell グリッドで node も網羅する。
            WALE_d<<<cuda_cfg.dimGrid_normalcell , cuda_cfg.dimBlock>>> (
                // mesh structure
                msh.nCells,
                var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
                var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

                // basic variables
                var.c_d["ro"] ,
                var.c_d["Ux"] ,
                var.c_d["Uy"] ,
                var.c_d["Uz"] ,
                var.c_d["P"]  , 
                var.c_d["T"]  ,
                var.c_d["Ht"] ,

                // gradient
                var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
                var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
                var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
                var.c_d["vis_turb"] , var.c_d["wall_dist"]
            ) ;
        } else if (cfg.LESmodel == 2) {
            // σ-model (Nicoud 2011)。WALE と同じ normalcell グリッドで cell/node 両対応。
            SIGMA_d<<<cuda_cfg.dimGrid_normalcell , cuda_cfg.dimBlock>>> (
                msh.nCells,
                var.c_d["volume"],
                var.c_d["ro"],
                var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
                var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
                var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
                var.c_d["vis_turb"]
            );
        }

    } else if (cfg.LESorRANS == 2 && cfg.RANSmodel == 1) { // RANS SST
        sst_eddy_viscosity_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
            msh.nCells,
            var.c_d["ro"],
            var.c_d["k"],
            var.c_d["omega"],
            var.c_d["vis_lam"],
            var.c_d["wall_dist"],
            var.c_d["dUxdx"], var.c_d["dUxdy"], var.c_d["dUxdz"],
            var.c_d["dUydx"], var.c_d["dUydy"], var.c_d["dUydz"],
            var.c_d["dUzdx"], var.c_d["dUzdy"], var.c_d["dUzdz"],
            cfg.isAxisymmetric,
            var.c_d["axisym_uy_over_r"],
            var.c_d["vis_turb"],
            (cfg.discretization == "node" && cfg.wallTreatmentSST == 1 && cfg.nodeOmegaWfDirichlet == 1)
                ? var.c_d["roOmega_wf"] : nullptr,
                // wf_irep_flag は ransWallFunction の init でのみ 0 クリアされるので、
                // node SST 壁関数が有効なときだけ渡す (未初期化値を読ませない)。
                sstWfActive ? var.c_d["wf_irep_flag"] : nullptr,
                // 2×2 分離: FORGE_WF_LIMITER_BYPASS (未設定=-1 従来 / 0=OFF / 1=ON)。
                // 壁関数が非有効な構成では常に -1 (従来経路) に固定する。
                !sstWfActive ? -1 : [] {
                    static const int v = [] {
                        const char* s = std::getenv("FORGE_WF_LIMITER_BYPASS");
                        const int x = s ? std::atoi(s) : -1;
                        if (s) printf("[EXPERIMENT] SST shear limiter bypass mode = %d "
                                      "(FORGE_WF_LIMITER_BYPASS)\n", x);
                        return x;
                    }();
                    return v;
                }());

        // SST-DES: vis_turb 計算直後に l_DDES とシールド診断を計算する (依存順・plan §4.7)。
        // DESmode==0 (既定) では呼ばない → 既存 SST と完全ビット不変。
        // C_DES は初版で C_DES_kw 固定 (F1 ブレンドは sst_eddy が F1 を持たないため Phase 2・plan §3.2)。
        if (cfg.DESmode > 0) {
            compute_des_length_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
                msh.nCells,
                cfg.DESmode,
                var.c_d["ro"],
                var.c_d["k"],
                var.c_d["omega"],
                var.c_d["vis_lam"],
                var.c_d["vis_turb"],
                var.c_d["wall_dist"],
                var.c_d["dUxdx"], var.c_d["dUxdy"], var.c_d["dUxdz"],
                var.c_d["dUydx"], var.c_d["dUydy"], var.c_d["dUydz"],
                var.c_d["dUzdx"], var.c_d["dUzdy"], var.c_d["dUzdz"],
                var.c_d["delta_les"],
                cfg.C_DES_kw,
                var.c_d["l_des"],
                var.c_d["fd_shield"],
                var.c_d["rd_des"],
                var.c_d["fe_iddes"]);
        }

    } else {
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["vis_turb"], 0.0, msh.nCells_all*sizeof(flow_float)));
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

}