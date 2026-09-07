#include "speciesTransport_d.cuh"

#include "scalarTransport_d.cuh"
#include "thermo_d.cuh"
#include "species_eos_coupling_d.cuh"   // 案C EOS クロス応答 (解析 JVP)

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace {

// device の roY / res_roY / transport_diag ポインタ配列 (flow_float*[nSpecies])。
// speciesInit_d で 1 度だけ構築 (kinetic 拡散カーネルが化学種ループするため)。
flow_float** g_roY_dev = nullptr;
flow_float** g_roYN_dev = nullptr;   // roYN: species ステップ始点ベースライン (speciesUpdateOuter で roY に一致)
flow_float** g_resroY_dev = nullptr;
flow_float** g_transdiag_dev = nullptr;
int          g_nSpecies = 0;

// face 整合再構成 (speciesFaceReconstruction==1) 用: Y{s} と セル勾配 ∇Y{s} のポインタ配列。
flow_float** g_Y_dev    = nullptr;
flow_float** g_dYdx_dev = nullptr;
flow_float** g_dYdy_dev = nullptr;
flow_float** g_dYdz_dev = nullptr;
flow_float** g_limiterY_dev = nullptr;   // ψ_Y[s] (Venkat on Y) ポインタ配列
// S3: convectiveFlux が書き出す upwind 再構成 face 組成 (layout [ip*nSpecies+s])。species 移流が同一面組成で読む。
flow_float*  g_Yface_dev = nullptr;
int          g_Yface_nPlanes = 0;

// 案C 用 scratch: dq_roY{s}_old のポインタ配列 (sweep の swap でポインタが変わるため毎回再構築) と
// セルごとの δp_Y [Pa]。speciesInit_d でなく案C 経路の初回に遅延確保。
flow_float** g_dqYold_dev = nullptr;   // device 上の flow_float*[nSpecies]
flow_float*  g_dpY_eos    = nullptr;   // device 上の flow_float[nCells_all]
int          g_dqYold_cap = 0;         // 確保済みポインタ数
geom_int     g_dpY_cap    = 0;         // 確保済みセル数

constexpr flow_float kSmall = static_cast<flow_float>(1.0e-30);

inline bool speciesEnabled(const variables& var)
{
    return var.nSpeciesRegistered >= 2;
}

// 化学種 s 用のスカラ輸送記述子を構築する (RANS buildScalarDescs と同形)。
// floor=0 (Y_s>=0)、sigma=0 (M2 は移流のみ。拡散は M4)。
ScalarTransportDesc buildSpeciesDesc(variables& var, int s)
{
    const std::string i = std::to_string(s);
    return ScalarTransportDesc{
        var.c_d["Y"+i], nullptr, nullptr, nullptr,  // dphidx/y/z: 化学種は汎用拡散(diffusion=0)未使用
        var.c_d["roY"+i], var.c_d["roY"+i+"N"], var.c_d["roY"+i+"M"],
        var.c_d["res_roY"+i], var.c_d["res_roY"+i+"_m"],
        var.c_d["src_jac_Y"+i], var.c_d["transport_diag_Y"+i],
        static_cast<flow_float>(0.0), static_cast<flow_float>(0.0),
        0  // 汎用拡散は使わない (化学種 Fick 拡散は species_diffusion_d で別途)
    };
}

__global__ void species_primitive_d(
    geom_int nCells_all,
    flow_float* ro,
    flow_float* roY,
    flow_float* Y)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells_all) {
        Y[ic] = roY[ic] / max(ro[ic], kSmall);
    }
}

// Neumann (zero-gradient) ghost 充填: roY[ig]=roY[ic], Y[ig]=Y[ic]。
__global__ void species_neumann_boundary_d(
    geom_int nb,
    geom_int* bplane_cell,
    geom_int* bplane_cell_ghst,
    flow_float* roY,
    flow_float* Y)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib < nb) {
        const geom_int ic = bplane_cell[ib];
        const geom_int ig = bplane_cell_ghst[ib];
        roY[ig] = roY[ic];
        Y[ig]   = Y[ic];
    }
}

// Dirichlet ghost 充填 (M5: 組成依存超音速入口): Y[ig]=Y_s^in, roY[ig]=ρ[ig]·Y_s^in。
// ρ[ig] は applyBconds (inlet カーネル) が設定済みの ghost 密度。Yb は per-face 入口組成 bvar。
__global__ void species_dirichlet_boundary_d(
    geom_int nb,
    geom_int* bplane_cell_ghst,
    flow_float* ro,
    flow_float* Yb,
    flow_float* roY,
    flow_float* Y)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib < nb) {
        const geom_int ig = bplane_cell_ghst[ib];
        const flow_float Yin = Yb[ib];
        Y[ig]   = Yin;
        roY[ig] = ro[ig] * Yin;
    }
}

// 実現可能性 + 再正規化: 各 ρY_s>=0 にクランプ後、Σ_s ρY_s = ρ となるよう再スケール (ΣY_s=1)。
__global__ void species_renormalize_d(
    geom_int nCells,
    int nSpecies,
    flow_float** roY,
    flow_float* ro)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        double sum = 0.0;
        for (int s = 0; s < nSpecies; s++) {
            flow_float v = roY[s][ic];
            if (v < 0.0) v = 0.0;
            roY[s][ic] = v;
            sum += (double)v;
        }
        const double factor = (double)ro[ic] / (sum > (double)kSmall ? sum : (double)kSmall);
        for (int s = 0; s < nSpecies; s++) {
            roY[s][ic] = (flow_float)((double)roY[s][ic] * factor);
        }
    }
}

// M4: 化学種 Fick 拡散 + 質量保存補正 (ΣJ=0) + エンタルピー拡散のエネルギー結合。
//   J_s = ρ D_s ∇Y_s (over-relaxed 法線), 補正 J_s* = J_s - Y_s Σ_k J_k (Σ J_s*=0),
//   res_roY_s += J_s*, res_roe += Σ_s h_s(T) J_s*。
//   D_s: diffMethod==1 で混合平均 (Chapman-Enskog), ==0 で定数 Schmidt D=μ/(ρ Sc)。
__global__ void species_diffusion_d(
    geom_int nCells,
    geom_int nNormalHaloPlanes,
    geom_int* normal_halo_planes,
    geom_int* plane_cells,
    geom_float* ccx, geom_float* ccy, geom_float* ccz,
    geom_float* fx, geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,
    const SpeciesThermo* sp, int nSpecies,
    flow_float** roY, flow_float** res_roY, flow_float** transport_diag,
    flow_float* ro, flow_float* T, flow_float* P, flow_float* vis_lam, flow_float* vis_turb,
    flow_float* res_roe,
    int diffMethod, flow_float Sc, flow_float Sc_t,
    int isNode, flow_float** dYdx, flow_float** dYdy, flow_float** dYdz)
{
    geom_int ih = blockDim.x * blockIdx.x + threadIdx.x;
    if (ih >= nNormalHaloPlanes) return;

    const geom_int ip  = normal_halo_planes[ih];
    const geom_int ic0 = plane_cells[2 * ip + 0];
    const geom_int ic1 = plane_cells[2 * ip + 1];

    const geom_float f   = fx[ip];
    const geom_float sxx = sx[ip], syy = sy[ip], szz = sz[ip], sss = ss[ip];

    // node モードの境界半割面 (ghost を含む面): 拡散流束を加えない (skip)。
    // 根拠 (plan diffusion-node-boundary-real-distance.md §3 (c), scalar_diffusion_first_order_d と同方針):
    //  - Dirichlet (固定組成入口) では境界ノード Y はピンで上書きされ無意味、内部ノードへの拡散は
    //    内部双対面 W↔I が実距離で運ぶ (本ループの非境界 plane で計算済)。
    //  - Neumann (壁 zero-grad / slip) は ∂Y/∂n=0 ＝半割面フラックス 0。
    //  → ghost mirror の dcc≈0 退化も ∇Y·S 弱形式の境界閉包依存も不要。エネルギー結合 (Σh_sJ_s) も
    //    半割面では 0。cell は ghost で正しく閉じるので従来どおり。
    if (isNode != 0 && (ic0 >= nCells || ic1 >= nCells)) {
        return;
    }

    const flow_float dccx = ccx[ic1] - ccx[ic0];
    const flow_float dccy = ccy[ic1] - ccy[ic0];
    const flow_float dccz = ccz[ic1] - ccz[ic0];
    const flow_float dcc  = sqrt(dccx*dccx + dccy*dccy + dccz*dccz);
    const flow_float denom = dccx*sxx + dccy*syy + dccz*szz;
    const flow_float Dsafe = (fabs(denom) < 1.0e-30) ? ((denom>=0)?1.0e-30:-1.0e-30) : denom;
    const flow_float delta = dcc * sss * sss / Dsafe;       // over-relaxed 法線

    const double ro0 = (double)max(ro[ic0], (flow_float)1.0e-30);
    const double ro1 = (double)max(ro[ic1], (flow_float)1.0e-30);
    const double ro_face = (double)f*ro0 + (1.0-(double)f)*ro1;
    const double T_face  = (double)f*(double)T[ic0] + (1.0-(double)f)*(double)T[ic1];
    const double P_face  = (double)f*(double)P[ic0] + (1.0-(double)f)*(double)P[ic1];

    // 面の質量分率 Yf (正規化) -> モル分率 X
    double Yf[THERMO_MAX_SPECIES], X[THERMO_MAX_SPECIES];
    double ysum = 0.0;
    for (int s=0;s<nSpecies;s++){
        double y = (double)f*((double)roY[s][ic0]/ro0) + (1.0-(double)f)*((double)roY[s][ic1]/ro1);
        if (y<0.0) y=0.0; Yf[s]=y; ysum+=y;
    }
    const double yinv = 1.0/(ysum>1.0e-30?ysum:1.0e-30);
    for (int s=0;s<nSpecies;s++) Yf[s]*=yinv;
    thermo_X_from_Y(sp, nSpecies, Yf, X);

    const double mu_face = (double)f*(double)vis_lam[ic0] + (1.0-(double)f)*(double)vis_lam[ic1];
    const double mut_face = (double)f*(double)vis_turb[ic0] + (1.0-(double)f)*(double)vis_turb[ic1];
    // 乱流化学種拡散 D_t = μ_t/(ρ Sc_t) は全種共通で加える。
    const double Dt = (mut_face > 0.0) ? mut_face/(ro_face*(double)Sc_t) : 0.0;

    // 各化学種の非補正 Fick flux J_s と Σ
    double Js[THERMO_MAX_SPECIES];
    double sumJ = 0.0;
    for (int s=0;s<nSpecies;s++){
        const double Ys0 = (double)roY[s][ic0]/ro0;
        const double Ys1 = (double)roY[s][ic1]/ro1;
        double D;
        if (diffMethod == 1) D = thermo_Dmix_species(sp, nSpecies, X, s, T_face, P_face);
        else                 D = mu_face/(ro_face*(double)Sc);
        D += Dt;  // 層流 (Fick/Sc) + 乱流 (μ_t/Sc_t)
        Js[s] = ro_face * D * ((Ys1 - Ys0)/(double)dcc) * (double)delta;
        sumJ += Js[s];
        // point-implicit 拡散対角 (各セル ρ で正規化)
        const double diag = ro_face * D * fabs((double)delta) / (double)max(dcc,(flow_float)1.0e-30);
        if (ic0 < nCells) atomicAdd(&transport_diag[s][ic0], (flow_float)(diag/ro0));
        if (ic1 < nCells) atomicAdd(&transport_diag[s][ic1], (flow_float)(diag/ro1));
    }

    // 補正 J_s* = J_s - Y_s Σ_k J_k (Σ J_s* = 0) と エネルギー結合 Σ h_s J_s*
    double q = 0.0;
    for (int s=0;s<nSpecies;s++){
        const double Jc = Js[s] - Yf[s]*sumJ;
        if (ic0 < nCells) atomicAdd(&res_roY[s][ic0],  (flow_float)Jc);
        if (ic1 < nCells) atomicAdd(&res_roY[s][ic1], -(flow_float)Jc);
        const double hs = thermo_h_mass(sp[s], T_face);   // NASA 絶対エンタルピー [J/kg]
        q += hs * Jc;
    }
    if (ic0 < nCells) atomicAdd(&res_roe[ic0],  (flow_float)q);
    if (ic1 < nCells) atomicAdd(&res_roe[ic1], -(flow_float)q);
}

// TP 多成分気体の組成-エネルギー整合補正 (M-stagger 修正)。
// speciesTimeIntegration で roY が更新されたとき、roe を同じ温度が維持されるよう補正する:
//   roe[ic] += Σ_s (roY_s[ic] - roYN_s[ic]) * h_s(T[ic])
// 背景: roY が変化しても T 一定を保つには roe を Σ_s Δ(roY_s)*h_s(T) だけ変化させる必要がある
// (等温条件下でのエンタルピー差)。staggered update では roe が未補正のまま speciesPrimitive
// の Newton 反転に入るため、TP h_mix が大きく負の成分 (H2O: ~-13.4 MJ/kg) では T が不物理的
// にジャンプし発散を引き起こす。本カーネルはその inconsistency を 1 次で除去する。
// 適用条件: thermalMethod==2 (TP) かつ nSpecies>=2。単成分/CPG では no-op。
__global__ void species_energy_correction_kernel(
    geom_int nCells,
    int nSpecies,
    flow_float* const* roY_dev,
    flow_float* const* roYN_dev,
    const flow_float* T,
    flow_float* roe,
    const SpeciesThermo* sp)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    const double Tic = (double)T[ic];
    double droe = 0.0;
    for (int s = 0; s < nSpecies; s++) {
        const double droY_s = (double)roY_dev[s][ic] - (double)roYN_dev[s][ic];
        if (droY_s != 0.0) {
            droe += droY_s * thermo_h_mass(sp[s], Tic);
        }
    }
    roe[ic] += (flow_float)droe;
}

// 緩和整合 scalar-DPLUR の 1 Jacobi sweep (speciesImplicitCoupling==1)。
// 流れ scalar-DPLUR (timeIntegration_d.cu implicit_defect_correction_d) を化学種 1 本にミラーする。
//   D·δ(ρY_s) = res_roY_s + Σ_f offdiag_f·δ(ρY_s)_old[nbr]
//   D    = V/Δτ + transport_diag[ic]              (対角=流出質量流束/ρ + 粘性時 Fick 拡散対角)
//   offdiag_f = max(∓ṁ_f,0)/ρ_nbr                 (非対角=流入質量流束/ρ_nbr, 1次風上の凍結 Jacobian)
// massflux[ip]>0 は ic0→ic1。ic=ic0 の流入は ṁ<0、ic=ic1 の流入は ṁ>0。ρ_nbr は基準 roN(=ρⁿ) で
// transport_diag の ρ と整合させる (flow commit 後の ρ ではなく)。dq_old は呼び出し側で memset 済み。
__global__ void species_dplur_sweep_d(
    flow_float implicit_relax,
    flow_float* dt_local,
    geom_int nCells,
    geom_float* vol,
    geom_int* plane_cells,
    geom_int* cell_planes_index,
    geom_int* cell_planes,
    flow_float* massflux,
    flow_float* roN,
    flow_float* res_roY,
    flow_float* transport_diag,
    flow_float* src_jac,
    flow_float* dq_old,
    flow_float* dq_new)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float dt_l = dt_local[ic];
        const geom_float v = vol[ic];

        flow_float neighbor = 0.0;
        const geom_int plane_begin = cell_planes_index[ic];
        const geom_int plane_end   = cell_planes_index[ic + 1];
        for (geom_int po = plane_begin; po < plane_end; ++po) {
            const geom_int ip  = cell_planes[po];
            const geom_int ic0 = plane_cells[2 * ip + 0];
            const geom_int ic1 = plane_cells[2 * ip + 1];
            const geom_int other_ic = (ic0 == ic) ? ic1 : ic0;
            const flow_float mdot = massflux[ip];
            // 流入 (other→ic) の質量流束。ic=ic0: ṁ<0 で流入、ic=ic1: ṁ>0 で流入。
            const flow_float inflow = (ic0 == ic) ? max(-mdot, static_cast<flow_float>(0.0))
                                                  : max( mdot, static_cast<flow_float>(0.0));
            if (other_ic < nCells) {
                const flow_float offdiag = inflow / max(roN[other_ic], static_cast<flow_float>(1.0e-30));
                neighbor += offdiag * dq_old[other_ic];
            }
        }

        // D = V/Δτ + V·src_jac + transport_diag (点陰的対角と同形。src_jac は化学種では 0)。
        const flow_float diag = max(
            static_cast<flow_float>(v / max(dt_l, static_cast<flow_float>(1.0e-30)))
                + static_cast<flow_float>(v) * src_jac[ic]
                + transport_diag[ic],
            static_cast<flow_float>(1.0e-30));

        dq_new[ic] = implicit_relax * (res_roY[ic] + neighbor) / diag;
    }
}

// 緩和整合 scalar-DPLUR の commit: ρY_s = ρY_s^N + δ(ρY_s)。実現可能性フロア (ρY_s>=0)。
// Σ_s ρY_s = ρ の再正規化は呼び出し側 (speciesRenormalize_d_wrapper) が行う。
__global__ void species_commit_correction_d(
    geom_int nCells,
    flow_float* roY,
    flow_float* roYN,
    flow_float* dq)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        roY[ic] = max(roYN[ic] + dq[ic], static_cast<flow_float>(0.0));
    }
}

// ===== 案C: block-triangular roe↔roY coupling (EOS クロス応答) =====

// 組成接空間への射影: z_s = δz*_s - Y_s Σ_r δz*_r (Σ_s z_s = 0)。dq[s] を in-place で z_s に書換える。
// Y_s = ρY_s^N/ρ^N (現組成で補正を分配; 単純等配分でない)。
__global__ void species_eos_project_tangent_d(
    geom_int nCells, int nSpecies,
    flow_float** dq, flow_float** roYN, flow_float* roN)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const double roc = max((double)roN[ic], (double)kSmall);
        double sumdq = 0.0;
        for (int s = 0; s < nSpecies; ++s) sumdq += (double)dq[s][ic];
        for (int s = 0; s < nSpecies; ++s) {
            const double Ys = (double)roYN[s][ic] / roc;
            dq[s][ic] = (flow_float)((double)dq[s][ic] - Ys * sumdq);
        }
    }
}

// セルごとの解析 EOS-JVP δp_Y = Σ_s (∂p/∂(ρY_s)) z_s を評価して dpY[ic] に格納。
// Y は ρY_s^N/ρ^N の正規化, T は現状態, ρ は ρ^N, z は射影済み dq。
__global__ void species_eos_dp_cell_d(
    geom_int nCells, int nSpecies, const SpeciesThermo* sp,
    flow_float** dq, flow_float** roYN, flow_float* roN, flow_float* T,
    flow_float* dpY)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const double roc = max((double)roN[ic], (double)kSmall);
        double Y[THERMO_MAX_SPECIES];
        double z[THERMO_MAX_SPECIES];
        double ysum = 0.0;
        for (int s = 0; s < nSpecies; ++s) {
            double y = (double)roYN[s][ic] / roc;
            if (y < 0.0) y = 0.0;
            Y[s] = y; ysum += y;
            z[s] = (double)dq[s][ic];
        }
        const double inv = 1.0 / (ysum > (double)kSmall ? ysum : (double)kSmall);
        for (int s = 0; s < nSpecies; ++s) Y[s] *= inv;
        double dp = 0.0, dT = 0.0;
        species_eos_cross_response(sp, nSpecies, Y, (double)T[ic], roc, z, &dp, &dT);
        dpY[ic] = (flow_float)dp;
    }
}

// クロスエネルギー流束 A_QY δY を res_roe へ移項する。エネルギー流束 mdot·H に対し
// δH = δp_Y/ρ (保存 ρe・ρ 固定) を上流差分で載せる (対流流束 res_roe_temp と同符号規約)。
__global__ void species_eos_cross_flux_d(
    geom_int nNormalPlanes, geom_int* plane_cells,
    flow_float* massflux, flow_float* dpY, flow_float* roN,
    flow_float* res_roe)
{
    const geom_int ip = blockDim.x * blockIdx.x + threadIdx.x;
    if (ip < nNormalPlanes) {
        const geom_int ic0 = plane_cells[2*ip+0];
        const geom_int ic1 = plane_cells[2*ip+1];
        const flow_float mdot = massflux[ip];
        const flow_float dH0  = dpY[ic0] / max(roN[ic0], kSmall);
        const flow_float dH1  = dpY[ic1] / max(roN[ic1], kSmall);
        const flow_float cross = 0.5f*(mdot + fabsf(mdot))*dH0
                               + 0.5f*(mdot - fabsf(mdot))*dH1;
        atomicAdd(&res_roe[ic0], -cross);
        atomicAdd(&res_roe[ic1],  cross);
    }
}

// 案C 最終 commit: ρY_s = ρY_s^N + z_s + Y_s^N δρ (δρ=ρ-ρ^N)。Σ_s δ(ρY_s)=δρ を満たす。
__global__ void species_eos_final_commit_d(
    geom_int nCells, int nSpecies,
    flow_float** roY, flow_float** roYN, flow_float** dq,
    flow_float* ro, flow_float* roN)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float dro = ro[ic] - roN[ic];
        const flow_float roc = max(roN[ic], kSmall);
        for (int s = 0; s < nSpecies; ++s) {
            const flow_float Ys = roYN[s][ic] / roc;
            const flow_float v  = roYN[s][ic] + dq[s][ic] + Ys * dro;
            roY[s][ic] = max(v, static_cast<flow_float>(0.0));
        }
    }
}

}  // namespace

static flow_float** g_srcjac_dev = nullptr;   // src_jac_Y{s} の device ポインタ配列 (化学反応ソース用)
flow_float** species_resroY_device_ptr() { return g_resroY_dev; }
flow_float** species_srcjac_device_ptr() { return g_srcjac_dev; }

void speciesInit_d(solverConfig& cfg, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var)) {
        g_roY_dev = nullptr;
        g_nSpecies = 0;
        return;
    }
    g_nSpecies = var.nSpeciesRegistered;

    std::vector<flow_float*> hroY(g_nSpecies), hroYN(g_nSpecies), hres(g_nSpecies), htd(g_nSpecies), hsj(g_nSpecies);
    for (int s = 0; s < g_nSpecies; s++) {
        const std::string i = std::to_string(s);
        hroY[s]  = var.c_d["roY"+i];
        hroYN[s] = var.c_d["roY"+i+"N"];
        hres[s]  = var.c_d["res_roY"+i];
        htd[s]   = var.c_d["transport_diag_Y"+i];
        hsj[s]   = var.c_d["src_jac_Y"+i];
    }
    const size_t pbytes = g_nSpecies*sizeof(flow_float*);
    gpuErrchk( cudaMalloc((void**)&g_roY_dev,      pbytes) );
    gpuErrchk( cudaMalloc((void**)&g_roYN_dev,     pbytes) );
    gpuErrchk( cudaMalloc((void**)&g_resroY_dev,   pbytes) );
    gpuErrchk( cudaMalloc((void**)&g_transdiag_dev, pbytes) );
    gpuErrchk( cudaMemcpy(g_roY_dev,       hroY.data(),  pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(g_roYN_dev,      hroYN.data(), pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(g_resroY_dev,    hres.data(),  pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(g_transdiag_dev, htd.data(),   pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMalloc((void**)&g_srcjac_dev, pbytes) );
    gpuErrchk( cudaMemcpy(g_srcjac_dev, hsj.data(), pbytes, cudaMemcpyHostToDevice) );

    // face 整合再構成用 Y / ∇Y ポインタ配列。
    std::vector<flow_float*> hY(g_nSpecies), hdx(g_nSpecies), hdy(g_nSpecies), hdz(g_nSpecies);
    for (int s = 0; s < g_nSpecies; s++) {
        const std::string i = std::to_string(s);
        hY[s]  = var.c_d["Y"+i];
        hdx[s] = var.c_d["dY"+i+"dx"]; hdy[s] = var.c_d["dY"+i+"dy"]; hdz[s] = var.c_d["dY"+i+"dz"];
    }
    gpuErrchk( cudaMalloc((void**)&g_Y_dev,    pbytes) ); gpuErrchk( cudaMemcpy(g_Y_dev,    hY.data(),  pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMalloc((void**)&g_dYdx_dev, pbytes) ); gpuErrchk( cudaMemcpy(g_dYdx_dev, hdx.data(), pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMalloc((void**)&g_dYdy_dev, pbytes) ); gpuErrchk( cudaMemcpy(g_dYdy_dev, hdy.data(), pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMalloc((void**)&g_dYdz_dev, pbytes) ); gpuErrchk( cudaMemcpy(g_dYdz_dev, hdz.data(), pbytes, cudaMemcpyHostToDevice) );
    std::vector<flow_float*> hlim(g_nSpecies);
    for (int s = 0; s < g_nSpecies; s++) hlim[s] = var.c_d["limiter_Y"+std::to_string(s)];
    gpuErrchk( cudaMalloc((void**)&g_limiterY_dev, pbytes) ); gpuErrchk( cudaMemcpy(g_limiterY_dev, hlim.data(), pbytes, cudaMemcpyHostToDevice) );

    std::cout << "speciesInit_d: built device roY/roYN/res/diag/Y/dY/limiterY[] for nSpecies=" << g_nSpecies << "\n";
}

// face 整合再構成用アクセサ
flow_float** species_Y_device_ptr()    { return g_Y_dev; }
flow_float** species_dYdx_device_ptr() { return g_dYdx_dev; }
flow_float** species_dYdy_device_ptr() { return g_dYdy_dev; }
flow_float** species_dYdz_device_ptr() { return g_dYdz_dev; }
flow_float** species_limiterY_device_ptr() { return g_limiterY_dev; }

// S3: convectiveFlux 用 face 組成バッファを確保し device ポインタを返す。
flow_float* species_Yface_alloc(int nPlanes)
{
    if (g_nSpecies < 2) return nullptr;
    if (g_Yface_dev == nullptr || g_Yface_nPlanes < nPlanes) {
        if (g_Yface_dev) cudaFree(g_Yface_dev);
        gpuErrchk( cudaMalloc((void**)&g_Yface_dev, (size_t)nPlanes*g_nSpecies*sizeof(flow_float)) );
        g_Yface_nPlanes = nPlanes;
    }
    return g_Yface_dev;
}

// S3: species 移流流束を **convectiveFlux が書いた face 組成** で組む (energy 流束と同一面組成)。
// 対角 transport_diag は 1 次風上のまま (defect-correction)。ΣY_face=1 なので Σ res_roY = res_ro。
__global__ void species_advection_faceY_d(
    geom_int nCells, geom_int nNormalHaloPlanes, geom_int* normal_halo_planes, geom_int* plane_cells,
    flow_float* ro, flow_float* massflux, int nSpecies, flow_float* Yface,
    flow_float** res_roY, flow_float** transport_diag,
    int isNode, flow_float** roY)
{
    geom_int ih = blockDim.x*blockIdx.x + threadIdx.x;
    if (ih < nNormalHaloPlanes) {
        const geom_int ip  = normal_halo_planes[ih];
        const geom_int ic0 = plane_cells[2*ip+0];
        const geom_int ic1 = plane_cells[2*ip+1];
        const flow_float mdot = massflux[ip];
        const flow_float d0 = max(mdot, (flow_float)0.0) / max(ro[ic0], (flow_float)1.0e-30);
        const flow_float d1 = max(-mdot,(flow_float)0.0) / max(ro[ic1], (flow_float)1.0e-30);
        // node 境界半割面 (ic1=ghost): 主ループ (SLAU/ROE) は境界半割面を除外するため Yface[ip] が
        // 未書込 (stale)。node は ghost を読まない設計なので、境界ノード ic0 自身の組成を面組成に使う。
        const bool nodeBnd = (isNode != 0 && ic1 >= nCells);
        for (int s = 0; s < nSpecies; ++s) {
            const flow_float Yf = nodeBnd
                ? (roY[s][ic0] / max(ro[ic0], (flow_float)1.0e-30))
                : Yface[(size_t)ip*nSpecies + s];   // 内部面 upwind は convectiveFlux 側で確定済み
            const flow_float flux = mdot * Yf;
            if (ic0 < nCells) { atomicAdd(&res_roY[s][ic0], -flux); atomicAdd(&transport_diag[s][ic0], d0); }
            if (ic1 < nCells) { atomicAdd(&res_roY[s][ic1],  flux); atomicAdd(&transport_diag[s][ic1], d1); }
        }
    }
}

void speciesAdvectionFaceY_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var) || g_Yface_dev == nullptr) return;
    dim3 dimGrid_nh = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));
    species_advection_faceY_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(
        msh.nCells, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
        var.c_d["ro"], var.p_d["massflux"], g_nSpecies, g_Yface_dev, g_resroY_dev, g_transdiag_dev,
        (cfg.discretization == "node") ? 1 : 0, g_roY_dev);
    gpuErrchk( cudaPeekAtLastError() );
}

// 化学種セル勾配 ∇Y{s} を Green-Gauss で計算する (calcGradient と同形)。speciesFaceReconstruction==1 のみ。
// 境界は Neumann ghost (applySpeciesBoundaries 済) を用い、内部面と同様に集計する。
__global__ void species_gradient_d(
    geom_int nCells, geom_int nPlanes, geom_int* plane_cells,
    geom_float* vol, geom_float* fx, geom_float* sx, geom_float* sy, geom_float* sz,
    int nSpecies, flow_float** Y, flow_float** dYdx, flow_float** dYdy, flow_float** dYdz)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;
    if (ip < nPlanes) {
        geom_int ic0 = plane_cells[2*ip+0];
        geom_int ic1 = plane_cells[2*ip+1];
        geom_float f = fx[ip];
        const geom_float sxx = sx[ip], syy = sy[ip], szz = sz[ip];
        for (int s = 0; s < nSpecies; ++s) {
            const flow_float Yf = f*Y[s][ic0] + (1.0-f)*Y[s][ic1];
            atomicAdd(&dYdx[s][ic0],  sxx*Yf); atomicAdd(&dYdy[s][ic0],  syy*Yf); atomicAdd(&dYdz[s][ic0],  szz*Yf);
            if (ic1 < nCells) { atomicAdd(&dYdx[s][ic1], -sxx*Yf); atomicAdd(&dYdy[s][ic1], -syy*Yf); atomicAdd(&dYdz[s][ic1], -szz*Yf); }
        }
    }
}

__global__ void species_gradient_normalize_d(
    geom_int nCells, geom_float* vol, int nSpecies,
    flow_float** dYdx, flow_float** dYdy, flow_float** dYdz)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float invv = 1.0/max(vol[ic], (flow_float)1.0e-30);
        for (int s = 0; s < nSpecies; ++s) {
            dYdx[s][ic] *= invv; dYdy[s][ic] *= invv; dYdz[s][ic] *= invv;
        }
    }
}

void speciesGradient_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;
    const int n = var.nSpeciesRegistered;
    for (int s = 0; s < n; ++s) {
        const std::string i = std::to_string(s);
        cudaMemset(var.c_d["dY"+i+"dx"], 0, msh.nCells*sizeof(flow_float));
        cudaMemset(var.c_d["dY"+i+"dy"], 0, msh.nCells*sizeof(flow_float));
        cudaMemset(var.c_d["dY"+i+"dz"], 0, msh.nCells*sizeof(flow_float));
    }
    flow_float* gvol = (cfg.isAxisymmetric == 1) ? var.c_d["A_planar"] : var.c_d["volume"];
    flow_float* gsx = (cfg.isAxisymmetric == 1) ? var.p_d["sx_planar"] : var.p_d["sx"];
    flow_float* gsy = (cfg.isAxisymmetric == 1) ? var.p_d["sy_planar"] : var.p_d["sy"];
    flow_float* gsz = (cfg.isAxisymmetric == 1) ? var.p_d["sz_planar"] : var.p_d["sz"];
    species_gradient_d<<<cuda_cfg.dimGrid_plane, cuda_cfg.dimBlock>>>(
        msh.nCells, msh.nPlanes, msh.map_plane_cells_d, gvol, var.p_d["fx"], gsx, gsy, gsz,
        n, g_Y_dev, g_dYdx_dev, g_dYdy_dev, g_dYdz_dev);
    species_gradient_normalize_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, gvol, n, g_dYdx_dev, g_dYdy_dev, g_dYdz_dev);
    gpuErrchk( cudaPeekAtLastError() );
}

flow_float** species_roY_device_ptr()
{
    return g_roY_dev;
}

void speciesEnergyCorrection_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (cfg.thermalMethod != 2) return;
    if (!speciesEnabled(var) || g_roY_dev == nullptr || g_roYN_dev == nullptr) return;

    dim3 dimGrid_c = dim3((geom_int)ceil(msh.nCells / (flow_float)cuda_cfg.blocksize));
    species_energy_correction_kernel<<<dimGrid_c, cuda_cfg.dimBlock>>>(
        msh.nCells,
        g_nSpecies,
        g_roY_dev,
        g_roYN_dev,
        var.c_d["T"],
        var.c_d["roe"],
        thermo_species_device_ptr());

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void speciesPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var)) return;

    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        species_primitive_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
            msh.nCells_all,
            var.c_d["ro"],
            var.c_d["roY"+i],
            var.c_d["Y"+i]);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void speciesBoundary_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, bcond& bc, mesh& msh, variables& var)
{
    (void)msh;
    if (!speciesEnabled(var)) return;
    if (bc.iPlanes.empty()) return;

    const geom_int nb = static_cast<geom_int>(bc.iPlanes.size());
    // M5: 入口種別は組成依存 Dirichlet (roY[ig]=ρ[ig]·Y_s^in)。入口組成 bvar_d["Y{s}"] は
    // readBcondConfig が inlet_* に対して登録済み。他種別は従来通り Neumann (zero-gradient)。
    const bool isInlet = bc.bcondKind.rfind("inlet_", 0) == 0;
    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        const auto ybIt = bc.bvar_d.find("Y"+i);
        if (isInlet && ybIt != bc.bvar_d.end()) {
            species_dirichlet_boundary_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(
                nb,
                bc.map_bplane_cell_ghst_d,
                var.c_d["ro"],
                ybIt->second,
                var.c_d["roY"+i],
                var.c_d["Y"+i]);
        } else {
            species_neumann_boundary_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(
                nb,
                bc.map_bplane_cell_d,
                bc.map_bplane_cell_ghst_d,
                var.c_d["roY"+i],
                var.c_d["Y"+i]);
        }
    }
}

void applySpeciesBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;

    for (auto& bc : msh.bconds) {
        speciesBoundary_d_wrapper(cfg, cuda_cfg, bc, msh, var);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void speciesTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;

    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roY"+i], 0, msh.nCells * sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["transport_diag_Y"+i], 0, msh.nCells * sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["src_jac_Y"+i], 0, msh.nCells * sizeof(flow_float)));
    }

    if (cfg.speciesFaceReconstruction >= 2 && g_Yface_dev != nullptr) {
        // S3: convectiveFlux が書いた同一 face 組成で移流 (energy 流束と整合)。diag は 1 次のまま。
        speciesAdvectionFaceY_d_wrapper(cfg, cuda_cfg, msh, var);
    } else {
        for (int s = 0; s < var.nSpeciesRegistered; s++) {
            const ScalarTransportDesc desc = buildSpeciesDesc(var, s);
            scalarTransportResidual_d(cfg, cuda_cfg, msh, var, desc);
        }
    }

    // M4: 粘性ケースのみ Fick 拡散 + ΣJ=0 補正 + エンタルピー拡散 (res_roe へ加算)。
    if (cfg.viscMethod != 0 && g_roY_dev != nullptr) {
        dim3 dimGrid_nh = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));
        species_diffusion_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(
            msh.nCells, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
            var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
            var.p_d["fx"], var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
            thermo_species_device_ptr(), g_nSpecies,
            g_roY_dev, g_resroY_dev, g_transdiag_dev,
            var.c_d["ro"], var.c_d["T"], var.c_d["P"], var.c_d["vis_lam"], var.c_d["vis_turb"],
            var.c_d["res_roe"],
            cfg.speciesDiffusionMethod, cfg.Sc, cfg.Sc_t,
            (cfg.discretization == "node") ? 1 : 0, g_dYdx_dev, g_dYdy_dev, g_dYdz_dev);
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void speciesTimeIntegration_d_wrapper(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;

    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const ScalarTransportDesc desc = buildSpeciesDesc(var, s);
        scalarTimeIntegration_d(loop, cfg, cuda_cfg, msh, var, desc);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void speciesRenormalize_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var) || g_roY_dev == nullptr) return;

    species_renormalize_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells,
        g_nSpecies,
        g_roY_dev,
        var.c_d["ro"]);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

bool speciesImplicitCoupled(solverConfig& cfg, variables& var)
{
    return cfg.speciesImplicitCoupling == 1 && speciesEnabled(var);
}

// 緩和整合 scalar-DPLUR ソルバ (speciesImplicitCoupling==1)。
// 凍結 res_roY/transport_diag/massflux/dt_local (assembleResidual で確定) に対し、各化学種の補正
// δ(ρY_s) を dq=0 から nStepInner 回 Jacobi sweep で緩和し (流れ block と同一 implicitRelax/sweep 回数)、
// ρY_s = ρY_s^N + δ(ρY_s) を commit する。呼び出し前に speciesUpdateOuter で ρY_s^N=ρY_s を取ること。
// commit 後の Σ_s ρY_s=ρ 再正規化・Y=ρY/ρ 同期は呼び出し側 (main.cpp) で行う。
void speciesImplicitDPLURSolve_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;

    const int nSpecies = var.nSpeciesRegistered;
    const size_t bytes = static_cast<size_t>(msh.nCells_all) * sizeof(flow_float);

    // dq (new/old) を 0 初期化 (古典 DPLUR は dq=0 から開始)。
    for (int s = 0; s < nSpecies; s++) {
        const std::string i = std::to_string(s);
        cudaMemset(var.c_d["dq_roY"+i],        0, bytes);
        cudaMemset(var.c_d["dq_roY"+i+"_old"], 0, bytes);
    }

    const int nSweep = std::max(1, cfg.nStepInner);
    for (int iSweep = 0; iSweep < nSweep; ++iSweep) {
        for (int s = 0; s < nSpecies; s++) {
            const std::string i = std::to_string(s);
            species_dplur_sweep_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
                cfg.implicitRelax,
                var.c_d["dt_local"],
                msh.nCells,
                var.c_d["volume"],
                msh.map_plane_cells_d,
                msh.map_cell_planes_index_d,
                msh.map_cell_planes_d,
                var.p_d["massflux"],
                var.c_d["roN"],
                var.c_d["res_roY"+i],
                var.c_d["transport_diag_Y"+i],
                var.c_d["src_jac_Y"+i],
                var.c_d["dq_roY"+i+"_old"],
                var.c_d["dq_roY"+i]);
        }
        // sweep 後に old↔new を swap (最終補正は dq_old 側に残る)。
        for (int s = 0; s < nSpecies; s++) {
            const std::string i = std::to_string(s);
            std::swap(var.c_d["dq_roY"+i], var.c_d["dq_roY"+i+"_old"]);
        }
    }

    // commit: ρY_s = ρY_s^N + δ(ρY_s) (最終補正は dq_roY{s}_old)。
    for (int s = 0; s < nSpecies; s++) {
        const std::string i = std::to_string(s);
        species_commit_correction_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
            msh.nCells,
            var.c_d["roY"+i],
            var.c_d["roY"+i+"N"],
            var.c_d["dq_roY"+i+"_old"]);
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

// ===== 案C: block-triangular roe↔roY coupling のオーケストレーション =====

bool speciesEOSCoupled(solverConfig& cfg, variables& var)
{
    return cfg.speciesImplicitCoupling == 2 && speciesEnabled(var);
}

// dq_roY{s}_old のポインタ配列を (再)構築する。sweep の swap でポインタが変わるため毎回呼ぶ。
static void rebuildDqYoldPtrs(variables& var, int nSpecies)
{
    if (g_dqYold_cap < nSpecies) {
        if (g_dqYold_dev) cudaFree(g_dqYold_dev);
        gpuErrchk( cudaMalloc((void**)&g_dqYold_dev, nSpecies*sizeof(flow_float*)) );
        g_dqYold_cap = nSpecies;
    }
    std::vector<flow_float*> h(nSpecies);
    for (int s = 0; s < nSpecies; ++s) h[s] = var.c_d["dq_roY"+std::to_string(s)+"_old"];
    gpuErrchk( cudaMemcpy(g_dqYold_dev, h.data(), nSpecies*sizeof(flow_float*), cudaMemcpyHostToDevice) );
}

void speciesEOSCrossPredictInject_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;
    const int nSpecies = var.nSpeciesRegistered;
    const size_t bytes = static_cast<size_t>(msh.nCells_all) * sizeof(flow_float);

    // --- species scalar-DPLUR sweep で δ(ρY_s)* を予測 (dq=0 から, commit しない) ---
    for (int s = 0; s < nSpecies; ++s) {
        const std::string i = std::to_string(s);
        cudaMemset(var.c_d["dq_roY"+i],        0, bytes);
        cudaMemset(var.c_d["dq_roY"+i+"_old"], 0, bytes);
    }
    const int nSweep = std::max(1, cfg.nStepInner);
    for (int iSweep = 0; iSweep < nSweep; ++iSweep) {
        for (int s = 0; s < nSpecies; ++s) {
            const std::string i = std::to_string(s);
            species_dplur_sweep_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
                cfg.implicitRelax, var.c_d["dt_local"], msh.nCells, var.c_d["volume"],
                msh.map_plane_cells_d, msh.map_cell_planes_index_d, msh.map_cell_planes_d,
                var.p_d["massflux"], var.c_d["roN"],
                var.c_d["res_roY"+i], var.c_d["transport_diag_Y"+i], var.c_d["src_jac_Y"+i],
                var.c_d["dq_roY"+i+"_old"], var.c_d["dq_roY"+i]);
        }
        for (int s = 0; s < nSpecies; ++s) {
            const std::string i = std::to_string(s);
            std::swap(var.c_d["dq_roY"+i], var.c_d["dq_roY"+i+"_old"]);
        }
    }
    // 最終補正 δ(ρY_s)* は dq_roY{s}_old に残る。ポインタ配列を再構築。
    rebuildDqYoldPtrs(var, nSpecies);

    // --- 接空間射影 z_s = δz*_s - Y_s Σδz* (dq_old を in-place で z に) ---
    species_eos_project_tangent_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, nSpecies, g_dqYold_dev, g_roYN_dev, var.c_d["roN"]);

    // --- 解析 δp_Y をセルごとに評価 ---
    if (g_dpY_cap < msh.nCells_all) {
        if (g_dpY_eos) cudaFree(g_dpY_eos);
        gpuErrchk( cudaMalloc((void**)&g_dpY_eos, bytes) );
        g_dpY_cap = msh.nCells_all;
    }
    cudaMemset(g_dpY_eos, 0, bytes);
    species_eos_dp_cell_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, nSpecies, thermo_species_device_ptr(),
        g_dqYold_dev, g_roYN_dev, var.c_d["roN"], var.c_d["T"], g_dpY_eos);

    // --- クロスエネルギー流束を res_roe へ移項 (内部 normal plane のみ) ---
    species_eos_cross_flux_d<<<cuda_cfg.dimGrid_plane, cuda_cfg.dimBlock>>>(
        msh.nNormalPlanes, msh.map_plane_cells_d,
        var.p_d["massflux"], g_dpY_eos, var.c_d["roN"], var.c_d["res_roe"]);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void speciesEOSFinalCommit_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var)) return;
    const int nSpecies = var.nSpeciesRegistered;
    rebuildDqYoldPtrs(var, nSpecies);   // PredictInject 後と同じ z (dq_old) を指す
    species_eos_final_commit_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, nSpecies, g_roY_dev, g_roYN_dev, g_dqYold_dev,
        var.c_d["ro"], var.c_d["roN"]);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void speciesUpdateOuter_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!speciesEnabled(var)) return;

    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        gpuErrchk( cudaMemcpy(var.c_d["roY"+i+"N"], var.c_d["roY"+i], bytes, cudaMemcpyDeviceToDevice) );
        gpuErrchk( cudaMemcpy(var.c_d["roY"+i+"M"], var.c_d["roY"+i], bytes, cudaMemcpyDeviceToDevice) );
    }
}

void speciesUpdateInner_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!speciesEnabled(var)) return;

    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        gpuErrchk( cudaMemcpy(var.c_d["roY"+i+"M"], var.c_d["roY"+i], bytes, cudaMemcpyDeviceToDevice) );
    }
}
