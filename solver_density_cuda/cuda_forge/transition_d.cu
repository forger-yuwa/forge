#include "transition_d.cuh"

#include "scalarTransport_d.cuh"
#include "periodicNode_d.cuh"

#include <array>
#include <cstdlib>
#include <stdexcept>

// Langtry–Menter 2009 γ–Re_θt (methods/turbulence/theory.md §11)。参照実装は tools/lm2009_reference.py (numpy 倍精度)。
// ここを直すときは SU2 の trans_sources.hpp / trans_correlations.hpp / CTransLMSolver.cpp の該当行と参照実装の両方を見ること。

namespace {

constexpr flow_float kCe1    = static_cast<flow_float>(1.0);
constexpr flow_float kCa1    = static_cast<flow_float>(2.0);
constexpr flow_float kCe2    = static_cast<flow_float>(50.0);
constexpr flow_float kCa2    = static_cast<flow_float>(0.06);
constexpr flow_float kCtheta = static_cast<flow_float>(0.03);
constexpr flow_float kGammaMin = static_cast<flow_float>(1.0e-4);   // SU2 CTransLMSolver.cpp lowerlimit[0]
constexpr flow_float kGammaMax = static_cast<flow_float>(1.0);
// $Re_{\theta t,min}$: 相関値と輸送変数 $\tilde{Re}_{\theta t}$ の下限。**既定 20 は LM2009/SU2 の推奨値** (`Corr_Ret_lim`)。
// `turbulence.transitionRethMin` で case ごとに変えられる (文献では 100〜200 に上げて負圧面の遷移を遅らせる調整例がある:
// Lin ら JGPP 6(3):9-15, 2014。notes/investigations/2026-09-23-vane-transition-literature.md)。下限を上げると遷移が遅れる。
constexpr flow_float kRethMinDefault = static_cast<flow_float>(20.0);
constexpr flow_float kTuMin    = static_cast<flow_float>(0.027);
constexpr flow_float kTuMax    = static_cast<flow_float>(100.0);    // 停滞点の保護 (U→0 で Tu→∞)。相関の適用範囲外なので値に意味は無い

inline bool transitionOn(const variables& var) { return var.transitionRegistered != 0; }

// Re_θc(Re_θt) (Langtry–Menter 2009 式 (41))。4 次多項式は桁落ちするので倍精度で評価する (節点あたり 5 回の乗算)。
__device__ __forceinline__ flow_float lm_re_theta_c(flow_float ret_f)
{
    const double r = static_cast<double>(ret_f);
    const double v = (r <= 1870.0)
        ? (-396.035e-2 + 10120.656e-4 * r - 868.230e-6 * r * r + 696.506e-9 * r * r * r - 174.105e-12 * r * r * r * r)
        : (r - (593.11 + 0.482 * (r - 1870.0)));
    return static_cast<flow_float>(v);
}

// F_length1(Re_θt) (同 式 (40))。
__device__ __forceinline__ flow_float lm_f_length1(flow_float ret_f)
{
    const double r = static_cast<double>(ret_f);
    double v;
    if      (r < 400.0)  v = 39.8189 - 119.270e-4 * r - 132.567e-6 * r * r;
    else if (r < 596.0)  v = 263.404 - 123.939e-2 * r + 194.548e-5 * r * r - 101.695e-8 * r * r * r;
    else if (r < 1200.0) v = 0.5 - 3.0e-4 * (r - 596.0);
    else                 v = 0.3188;
    return static_cast<flow_float>(v);
}

// 自由流 (λ_θ=0) の Re_θt(Tu)。入口値と初期値に使う (SU2 CTransLMSolver.cpp:112-124 と同じ形)。
__device__ __forceinline__ flow_float lm_re_theta_t_freestream(flow_float tu, flow_float rethMin)
{
    flow_float v;
    if (tu <= static_cast<flow_float>(1.3)) {
        v = static_cast<flow_float>(1173.51) - static_cast<flow_float>(589.428) * tu + static_cast<flow_float>(0.2196) / (tu * tu);
    } else {
        v = static_cast<flow_float>(331.5) * pow(tu - static_cast<flow_float>(0.5658), static_cast<flow_float>(-0.671));
    }
    return max(v, rethMin);
}

__device__ __forceinline__ flow_float lm_local_tu(flow_float k, flow_float umag)
{
    const flow_float tu = static_cast<flow_float>(100.0) * sqrt(static_cast<flow_float>(2.0 / 3.0) * max(k, static_cast<flow_float>(0.0)))
                        / max(umag, static_cast<flow_float>(1.0e-30));
    return min(max(tu, kTuMin), kTuMax);
}

// 原始量の復元 + 初回初期化。
__global__ void transition_primitive_d(
    geom_int nCells, int doInit, flow_float rethMin,
    flow_float* ro, flow_float* Ux, flow_float* Uy, flow_float* Uz, flow_float* k,
    flow_float* roGamma, flow_float* roReth, flow_float* gammaTr, flow_float* reTheta, flow_float* gammaEff)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const flow_float r = max(ro[ic], static_cast<flow_float>(1.0e-30));
    if (doInit != 0) {
        const flow_float umag = sqrt(Ux[ic] * Ux[ic] + Uy[ic] * Uy[ic] + Uz[ic] * Uz[ic]);
        roGamma[ic] = r;
        roReth[ic]  = r * lm_re_theta_t_freestream(lm_local_tu(k[ic], umag), rethMin);
        gammaEff[ic] = static_cast<flow_float>(1.0);
    }
    gammaTr[ic] = min(max(roGamma[ic] / r, kGammaMin), kGammaMax);
    reTheta[ic] = max(roReth[ic] / r, rethMin);
    // ソース評価の前 (初期出力) は γ_sep が未定なので γ を入れておく。評価後は transition_lm_source_d が毎反復上書きする。
    if (gammaEff[ic] <= static_cast<flow_float>(0.0)) gammaEff[ic] = gammaTr[ic];
}

// node 入口ピン: γ=1, Re_θt = 自由流相関 (局所 Tu)。
__global__ void transition_inlet_pin_d(
    geom_int nCells, flow_float rethMin, flow_float* pin,
    flow_float* ro, flow_float* Ux, flow_float* Uy, flow_float* Uz, flow_float* k,
    flow_float* roGamma, flow_float* roReth, flow_float* gammaTr, flow_float* reTheta)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    if (pin[ic] != static_cast<flow_float>(1.0)) return;
    const flow_float umag = sqrt(Ux[ic] * Ux[ic] + Uy[ic] * Uy[ic] + Uz[ic] * Uz[ic]);
    const flow_float ret  = lm_re_theta_t_freestream(lm_local_tu(k[ic], umag), rethMin);
    gammaTr[ic] = static_cast<flow_float>(1.0);
    reTheta[ic] = ret;
    roGamma[ic] = ro[ic];
    roReth[ic]  = ro[ic] * ret;
}

// ソース・陰的対角・γ_eff。体積 vol は周期 node では部分体積 (ransSource と同じ理由)。
__global__ void transition_lm_source_d(
    geom_int nCells, geom_float* vol, flow_float rethMin,
    flow_float* ro, flow_float* Ux, flow_float* Uy, flow_float* Uz, flow_float* sonic,
    flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
    flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
    flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
    flow_float* k, flow_float* omega, flow_float* vis_lam, flow_float* wall_dist,
    flow_float* gammaTr, flow_float* reTheta,
    geom_int* wall_flag, flow_float* pin,
    flow_float* res_roGamma, flow_float* res_roReth,
    flow_float* src_jac_gamma, flow_float* src_jac_reth,
    flow_float* gammaEff,
    // 診断 (output.level 2 / extraFields で確保されたときだけ非 null)
    flow_float* d_fonset, flow_float* d_flength, flow_float* d_ftheta, flow_float* d_rethcorr, flow_float* d_gammasep,
    flow_float* d_pgamma, flow_float* d_egamma, flow_float* d_ptheta, flow_float* d_iter)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    const flow_float gam = gammaTr[ic];
    const flow_float ret = reTheta[ic];
    // 既定 (ソースを評価しない節点): γ_eff = γ、ソース・対角 0。未初期化値を残さない。
    gammaEff[ic]      = gam;
    src_jac_gamma[ic] = static_cast<flow_float>(0.0);
    src_jac_reth[ic]  = static_cast<flow_float>(0.0);
    if (d_fonset != nullptr) {
        d_fonset[ic] = 0; d_flength[ic] = 0; d_ftheta[ic] = 0; d_rethcorr[ic] = 0; d_gammasep[ic] = 0;
        d_pgamma[ic] = 0; d_egamma[ic] = 0; d_ptheta[ic] = 0; d_iter[ic] = 0;
    }

    const bool pinned = (pin != nullptr && pin[ic] == static_cast<flow_float>(1.0));
    if (pinned) {
        res_roGamma[ic] = static_cast<flow_float>(0.0);
        res_roReth[ic]  = static_cast<flow_float>(0.0);
        return;
    }

    const flow_float rho  = max(ro[ic], static_cast<flow_float>(1.0e-30));
    const flow_float mu   = vis_lam[ic];
    const flow_float dist = wall_dist[ic];
    const flow_float ux = Ux[ic], uy = Uy[ic], uz = Uz[ic];
    const flow_float umag = sqrt(ux * ux + uy * uy + uz * uz);
    const bool onWall = (wall_flag != nullptr && wall_flag[ic] == 1);
    // 壁ノード (u=0)・壁距離 0・停滞点 (U<1e-6 a) はソースを評価しない (1/U, 1/U² のゼロ割)。輸送だけが効く = 法線勾配ゼロ。
    if (onWall || dist <= static_cast<flow_float>(1.0e-10) || umag < static_cast<flow_float>(1.0e-6) * sonic[ic]) return;

    // ひずみ・渦度 (SU2: StrainMag = sqrt(2 SijSij), VorticityMag)
    const flow_float S11 = dUxdx[ic], S22 = dUydy[ic], S33 = dUzdz[ic];
    const flow_float S12 = static_cast<flow_float>(0.5) * (dUxdy[ic] + dUydx[ic]);
    const flow_float S13 = static_cast<flow_float>(0.5) * (dUxdz[ic] + dUzdx[ic]);
    const flow_float S23 = static_cast<flow_float>(0.5) * (dUydz[ic] + dUzdy[ic]);
    const flow_float S = sqrt(static_cast<flow_float>(2.0) * (S11 * S11 + S22 * S22 + S33 * S33
                       + static_cast<flow_float>(2.0) * (S12 * S12 + S13 * S13 + S23 * S23)));
    const flow_float wx = dUzdy[ic] - dUydz[ic];
    const flow_float wy = dUxdz[ic] - dUzdx[ic];
    const flow_float wz = dUydx[ic] - dUxdy[ic];
    const flow_float Om = sqrt(wx * wx + wy * wy + wz * wz);

    const flow_float k_c = max(k[ic], static_cast<flow_float>(0.0));
    const flow_float w_c = max(omega[ic], static_cast<flow_float>(1.0e-12));
    const flow_float tu  = max(static_cast<flow_float>(100.0) * sqrt(static_cast<flow_float>(2.0 / 3.0) * k_c) / umag, kTuMin);

    const flow_float rec     = max(lm_re_theta_c(ret), static_cast<flow_float>(1.0e-6));
    const flow_float r_omega = rho * dist * dist * w_c / mu;
    const flow_float f_sub   = exp(-(r_omega / static_cast<flow_float>(200.0)) * (r_omega / static_cast<flow_float>(200.0)));
    const flow_float f_len   = lm_f_length1(ret) * (static_cast<flow_float>(1.0) - f_sub) + static_cast<flow_float>(40.0) * f_sub;

    const flow_float r_t  = rho * k_c / mu / w_c;
    const flow_float re_v = rho * dist * dist * S / mu;
    // F_onset2 = min(max(F1, F1^4), 2)。F1^4 の overflow を避けるため先に 2 で頭打ち (結果は同じ)。
    const flow_float fo1  = min(re_v / (static_cast<flow_float>(2.193) * rec), static_cast<flow_float>(2.0));
    const flow_float fo2  = min(max(fo1, fo1 * fo1 * fo1 * fo1), static_cast<flow_float>(2.0));
    const flow_float rt25 = min(r_t / static_cast<flow_float>(2.5), static_cast<flow_float>(1.0e3));
    const flow_float fo3  = max(static_cast<flow_float>(1.0) - rt25 * rt25 * rt25, static_cast<flow_float>(0.0));
    const flow_float f_onset = max(fo2 - fo3, static_cast<flow_float>(0.0));

    // 流線方向の加速 dU/ds = (U_i U_j / |U|^2) dU_i/dx_j
    const flow_float inv_u = static_cast<flow_float>(1.0) / umag;
    const flow_float dUdx = (ux * dUxdx[ic] + uy * dUydx[ic] + uz * dUzdx[ic]) * inv_u;
    const flow_float dUdy = (ux * dUxdy[ic] + uy * dUydy[ic] + uz * dUzdy[ic]) * inv_u;
    const flow_float dUdz = (ux * dUxdz[ic] + uy * dUydz[ic] + uz * dUzdz[ic]) * inv_u;
    const flow_float duds = (ux * dUdx + uy * dUdy + uz * dUdz) * inv_u;

    const flow_float nu      = mu / rho;
    const flow_float tscale  = static_cast<flow_float>(500.0) * nu * inv_u * inv_u;
    const flow_float theta_bl = ret * nu * inv_u;
    const flow_float delta   = static_cast<flow_float>(50.0) * Om * dist * inv_u * (static_cast<flow_float>(7.5) * theta_bl)
                             + static_cast<flow_float>(1.0e-20);
    const flow_float rw5     = min(r_omega / static_cast<flow_float>(1.0e5), static_cast<flow_float>(1.0e3));
    const flow_float f_wake  = exp(-rw5 * rw5);
    const flow_float dd      = min(dist / delta, static_cast<flow_float>(1.0e3));   // (d/δ)^4 の overflow 保護 (exp(-1e12)=0 で結果不変)
    const flow_float var1    = (gam - static_cast<flow_float>(1.0) / kCe2) / (static_cast<flow_float>(1.0) - static_cast<flow_float>(1.0) / kCe2);
    const flow_float f_theta = min(max(f_wake * exp(-dd * dd * dd * dd), static_cast<flow_float>(1.0) - var1 * var1), static_cast<flow_float>(1.0));
    const flow_float rt4     = min(r_t / static_cast<flow_float>(4.0), static_cast<flow_float>(1.0e3));
    const flow_float f_turb  = exp(-rt4 * rt4 * rt4 * rt4);

    // Re_θt(Tu, λ_θ) の不動点反復 (SU2 と同じ初期値 20・上限 100 回)。
    flow_float corr = rethMin, corr_old = static_cast<flow_float>(0.0);
    const flow_float tu15 = exp(-pow(tu / static_cast<flow_float>(1.5), static_cast<flow_float>(1.5)));
    const flow_float tu05 = exp(-tu / static_cast<flow_float>(0.5));
    const flow_float base = (tu <= static_cast<flow_float>(1.3))
        ? (static_cast<flow_float>(1173.51) - static_cast<flow_float>(589.428) * tu + static_cast<flow_float>(0.2196) / (tu * tu))
        : (static_cast<flow_float>(331.5) * pow(tu - static_cast<flow_float>(0.5658), static_cast<flow_float>(-0.671)));
    int it = 0;
    for (; it < 100; ++it) {
        const flow_float th  = corr * nu * inv_u;
        flow_float lam = th * th / nu * duds;
        lam = min(max(lam, static_cast<flow_float>(-0.1)), static_cast<flow_float>(0.1));
        flow_float fl;
        if (lam <= static_cast<flow_float>(0.0)) {
            fl = static_cast<flow_float>(1.0) - (static_cast<flow_float>(-12.986) * lam - static_cast<flow_float>(123.66) * lam * lam
                 - static_cast<flow_float>(405.689) * lam * lam * lam) * tu15;
        } else {
            fl = static_cast<flow_float>(1.0) + static_cast<flow_float>(0.275) * (static_cast<flow_float>(1.0) - exp(static_cast<flow_float>(-35.0) * lam)) * tu05;
        }
        corr = max(base * fl, rethMin);
        if (fabs(corr - corr_old) <= static_cast<flow_float>(1.0e-6) * corr) { ++it; break; }
        corr_old = corr;
    }

    const flow_float Pg = f_len * kCa1 * rho * S * sqrt(max(f_onset * gam, static_cast<flow_float>(0.0))) * (static_cast<flow_float>(1.0) - kCe1 * gam);
    const flow_float Eg = kCa2 * rho * Om * gam * f_turb * (kCe2 * gam - static_cast<flow_float>(1.0));
    const flow_float Pt = kCtheta * rho / tscale * (corr - ret) * (static_cast<flow_float>(1.0) - f_theta);

    const geom_float v = vol[ic];
    atomicAdd(&res_roGamma[ic], (Pg - Eg) * v);
    atomicAdd(&res_roReth[ic],  Pt * v);

    // 陰的対角 = 体積あたりのソース微分の負の部分 (項ごと)。∂P_γ/∂γ の +0.5/√γ は陽的に残す (γ→0 で正側に発散するため対角に入らない)。
    //   ∂P_γ/∂(ργ) ∋ -1.5 c_e1 F_len c_a1 S √(F_onset γ)      ∂E_γ/∂(ργ) = c_a2 Ω f_turb (2 c_e2 γ - 1)
    //   ∂P_θt/∂(ρRe_θt) = -(c_θt/t)(1 - F_θt)
    const flow_float jg_p = static_cast<flow_float>(1.5) * kCe1 * f_len * kCa1 * S * sqrt(max(f_onset * gam, static_cast<flow_float>(0.0)));
    const flow_float jg_e = kCa2 * Om * f_turb * (static_cast<flow_float>(2.0) * kCe2 * gam - static_cast<flow_float>(1.0));
    src_jac_gamma[ic] = jg_p + max(jg_e, static_cast<flow_float>(0.0));
    src_jac_reth[ic]  = kCtheta / tscale * (static_cast<flow_float>(1.0) - f_theta);

    // 剥離誘起遷移 γ_sep (SU2 CTransLMSolver::Postprocessing)
    const flow_float rt20 = min(r_t / static_cast<flow_float>(20.0), static_cast<flow_float>(1.0e3));
    const flow_float f_reattach = exp(-rt20 * rt20 * rt20 * rt20);
    flow_float gsep = static_cast<flow_float>(2.0) * max(re_v / (static_cast<flow_float>(3.235) * rec) - static_cast<flow_float>(1.0), static_cast<flow_float>(0.0)) * f_reattach;
    gsep = min(gsep, static_cast<flow_float>(2.0)) * f_theta;
    gsep = min(max(gsep, static_cast<flow_float>(0.0)), static_cast<flow_float>(2.0));
    gammaEff[ic] = max(gam, gsep);

    if (d_fonset != nullptr) {
        d_fonset[ic] = f_onset; d_flength[ic] = f_len; d_ftheta[ic] = f_theta; d_rethcorr[ic] = corr; d_gammasep[ic] = gsep;
        d_pgamma[ic] = Pg; d_egamma[ic] = Eg; d_ptheta[ic] = Pt; d_iter[ic] = static_cast<flow_float>(it);
    }
}

// point-implicit 更新: D = V/Δτ + V·src_jac + transport_diag、δ = relax·res/D。上下限は更新後の保存量に直接掛ける。
__global__ void transition_point_implicit_d(
    geom_int nCells, geom_float* vol, flow_float* dt_local, flow_float relax, flow_float rethMin,
    flow_float* ro,
    flow_float* roGamma, flow_float* roReth,
    flow_float* res_roGamma, flow_float* res_roReth,
    flow_float* src_jac_gamma, flow_float* src_jac_reth,
    flow_float* transport_diag_gamma, flow_float* transport_diag_reth)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const flow_float v    = vol[ic];
    const flow_float dt_l = max(dt_local[ic], static_cast<flow_float>(1.0e-30));
    const flow_float r    = max(ro[ic], static_cast<flow_float>(1.0e-30));
    const flow_float Dg = v / dt_l + v * src_jac_gamma[ic] + transport_diag_gamma[ic];
    const flow_float Dt = v / dt_l + v * src_jac_reth[ic]  + transport_diag_reth[ic];
    const flow_float dg = relax * res_roGamma[ic] / max(Dg, static_cast<flow_float>(1.0e-30));
    const flow_float dt = relax * res_roReth[ic]  / max(Dt, static_cast<flow_float>(1.0e-30));
    roGamma[ic] = min(max(roGamma[ic] + dg, r * kGammaMin), r * kGammaMax);
    roReth[ic]  = max(roReth[ic] + dt, r * rethMin);
}

std::array<ScalarTransportDesc, 2> buildTransitionDescs(variables& var)
{
    // γ: σ_f=1 → μ + μ_t。Re_θt: σ_θt=2 → 2(μ + μ_t)。勾配・RK 用バッファは使わない (node 境界半割面の拡散は skip、更新は point-implicit のみ)。
    std::array<ScalarTransportDesc, 2> d = {{
        {var.c_d["gammaTr"], nullptr, nullptr, nullptr, var.c_d["roGamma"], nullptr, nullptr, var.c_d["res_roGamma"], nullptr,
         var.c_d["src_jac_gamma"], var.c_d["transport_diag_gamma"], static_cast<flow_float>(0.0), static_cast<flow_float>(1.0), 1},
        {var.c_d["reTheta"], nullptr, nullptr, nullptr, var.c_d["roReth"], nullptr, nullptr, var.c_d["res_roReth"], nullptr,
         var.c_d["src_jac_reth"], var.c_d["transport_diag_reth"], static_cast<flow_float>(0.0), static_cast<flow_float>(2.0), 1}
    }};
    d[1].sigma_lam = static_cast<flow_float>(2.0);
    return d;
}

}  // namespace

void transitionValidateConfig(const solverConfig& cfg)
{
    if (!cfg.transitionEnabled()) return;
    auto fail = [](const std::string& why) {
        throw std::runtime_error("turbulence.transition: lm2009 is not accepted with this configuration: " + why
                                 + " (plan turbulence-transition-lm2009 §4.1)");
    };
    if (cfg.discretization != "node")                           fail("mesh.discretization must be node");
    if (!(cfg.LESorRANS == 2 && cfg.RANSmodel == 1))            fail("turbulence.model must be sst");
    if (cfg.DESmode != 0)                                       fail("DES is not supported");
    if (cfg.wallTreatmentSST != 0)                              fail("wallTreatmentSST must be 0 (low-Re wall)");
    if (cfg.isAxisymmetric != 0)                                fail("axisymmetric runs are not verified");
    if (cfg.sstEnergyIncludesK != 0)                            fail("sstEnergyIncludesK must be 0");
    if (cfg.scalarDiffusion != 1)                               fail("scalarDiffusion must be 1");
    if (!(cfg.isImplicit == 1 && cfg.unsteady == 0))            fail("only the steady implicit (point-implicit scalar) path is implemented");
}

void transitionPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!transitionOn(var)) return;
    const int doInit = (var.transitionNeedsInit != 0) ? 1 : 0;
    transition_primitive_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, doInit, static_cast<flow_float>(cfg.transitionRethMin),
        var.c_d["ro"], var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"], var.c_d["k"],
        var.c_d["roGamma"], var.c_d["roReth"], var.c_d["gammaTr"], var.c_d["reTheta"], var.c_d["gammaEff"]);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
    if (doInit != 0) {
        printf("[transition] no roGamma/roReth in the input: initialised gamma=1, Re_theta_t from the local-Tu free-stream correlation\n");
        var.transitionNeedsInit = 0;
    }
}

void applyTransitionBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!transitionOn(var)) return;
    transition_inlet_pin_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, static_cast<flow_float>(cfg.transitionRethMin), var.c_d["scalarDirichletPin"],
        var.c_d["ro"], var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"], var.c_d["k"],
        var.c_d["roGamma"], var.c_d["roReth"], var.c_d["gammaTr"], var.c_d["reTheta"]);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void transitionTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!transitionOn(var)) return;
    for (const char* nm : {"res_roGamma", "res_roReth", "transport_diag_gamma", "transport_diag_reth"}) {
        CHECK_CUDA_ERROR(cudaMemset(var.c_d[nm], 0, msh.nCells * sizeof(flow_float)));
    }
    const auto descs = buildTransitionDescs(var);
    scalarTransportResidualMulti_d(cfg, cuda_cfg, msh, var, descs.data(), 2);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void transitionSource_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!transitionOn(var)) return;
    auto opt = [&](const char* nm) -> flow_float* { auto it = var.c_d.find(nm); return (it == var.c_d.end()) ? nullptr : it->second; };
    transition_lm_source_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        msh.nCells,
        (msh.volumePartial_d != nullptr) ? msh.volumePartial_d : var.c_d["volume"],
        static_cast<flow_float>(cfg.transitionRethMin),
        var.c_d["ro"], var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"], var.c_d["sonic"],
        var.c_d["dUxdx"], var.c_d["dUxdy"], var.c_d["dUxdz"],
        var.c_d["dUydx"], var.c_d["dUydy"], var.c_d["dUydz"],
        var.c_d["dUzdx"], var.c_d["dUzdy"], var.c_d["dUzdz"],
        var.c_d["k"], var.c_d["omega"], var.c_d["vis_lam"], var.c_d["wall_dist"],
        var.c_d["gammaTr"], var.c_d["reTheta"],
        (cfg.discretization == "node") ? msh.wall_flag_d : nullptr,
        var.c_d["scalarDirichletPin"],
        var.c_d["res_roGamma"], var.c_d["res_roReth"],
        var.c_d["src_jac_gamma"], var.c_d["src_jac_reth"],
        var.c_d["gammaEff"],
        opt("lmFonset"), opt("lmFlength"), opt("lmFtheta"), opt("lmRethCorr"), opt("lmGammaSep"),
        opt("lmPgamma"), opt("lmEgamma"), opt("lmPtheta"), opt("lmCorrIter"));
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void applyTransitionPointImplicit_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!transitionOn(var)) return;
    transition_point_implicit_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, var.c_d["volume"], var.c_d["dt_local"], static_cast<flow_float>(cfg.implicitRelax), static_cast<flow_float>(cfg.transitionRethMin),
        var.c_d["ro"],
        var.c_d["roGamma"], var.c_d["roReth"],
        var.c_d["res_roGamma"], var.c_d["res_roReth"],
        var.c_d["src_jac_gamma"], var.c_d["src_jac_reth"],
        var.c_d["transport_diag_gamma"], var.c_d["transport_diag_reth"]);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
    periodicMirrorTransitionState_d_wrapper(cfg, cuda_cfg, msh, var);
}
