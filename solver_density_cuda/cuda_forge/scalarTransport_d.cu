#include "scalarTransport_d.cuh"

namespace {

__global__ void scalar_advection_first_order_d(
    geom_int nCells,
    geom_int nNormalHaloPlanes,
    int isNode,
    geom_int* normal_halo_planes,
    geom_int* plane_cells,
    flow_float* phi,
    flow_float* ro,
    flow_float* massflux,
    flow_float* res_rho_phi,
    flow_float* transport_diag)
{
    geom_int ih = blockDim.x * blockIdx.x + threadIdx.x;

    if (ih < nNormalHaloPlanes) {
        const geom_int ip = normal_halo_planes[ih];
        const geom_int ic0 = plane_cells[2 * ip + 0];
        const geom_int ic1 = plane_cells[2 * ip + 1];

        const flow_float mdot = massflux[ip];

        // node モード境界半割面 (ic1=ghost): node は ghost を読まない設計。流入側の上流値は
        // 境界ノード ic0 自身の値 (BC が Dirichlet/壁ではそこへピン、Neumann はゼロ勾配で
        // 境界値=ic0 値) を使い、ghost phi[ic1] への依存を断つ (cell は従来どおり ghost)。
        flow_float phi_ext = phi[ic1];
        if (isNode != 0 && ic1 >= nCells) phi_ext = phi[ic0];
        const flow_float phi_upwind = (mdot >= 0.0) ? phi[ic0] : phi_ext;
        const flow_float flux = mdot * phi_upwind;

        // point-implicit 移流対角: 1次風上の -∂res/∂(ρφ) は流出側セルに max(±ṁ,0)/ρ [m³/s]。
        // 流出 (outflow) のみ正の対角を生み、φ の陽的移流 stiff 性を陰化する（拡散と同じ対角に集約）。
        if (ic0 < nCells) {
            atomicAdd(&res_rho_phi[ic0], -flux);
            atomicAdd(&transport_diag[ic0], max(mdot, static_cast<flow_float>(0.0)) / max(ro[ic0], static_cast<flow_float>(1.0e-30)));
        }
        if (ic1 < nCells) {
            atomicAdd(&res_rho_phi[ic1], flux);
            atomicAdd(&transport_diag[ic1], max(-mdot, static_cast<flow_float>(0.0)) / max(ro[ic1], static_cast<flow_float>(1.0e-30)));
        }
    }
}

__global__ void scalar_diffusion_first_order_d(
    geom_int nCells,
    geom_int nNormalHaloPlanes,
    geom_int nWallHaloPlanes,
    int isNode,
    geom_int* normal_halo_planes,
    geom_int* plane_cells,
    geom_float* ccx,
    geom_float* ccy,
    geom_float* ccz,
    geom_float* fx,
    geom_float* sx,
    geom_float* sy,
    geom_float* sz,
    geom_float* ss,
    flow_float* phi,
    flow_float* dphidx,
    flow_float* dphidy,
    flow_float* dphidz,
    flow_float* ro,
    flow_float* vis_lam,
    flow_float* vis_turb,
    flow_float sigma,
    flow_float* res_rho_phi,
    flow_float* transport_diag,
    flow_float sigma2,
    flow_float* F1blend)
{
    geom_int ih = blockDim.x * blockIdx.x + threadIdx.x;

    if (ih < nNormalHaloPlanes) {
        const geom_int ip = normal_halo_planes[ih];
        const geom_int ic0 = plane_cells[2 * ip + 0];
        const geom_int ic1 = plane_cells[2 * ip + 1];

        const geom_float f = fx[ip];
        const geom_float sxx = sx[ip];
        const geom_float syy = sy[ip];
        const geom_float szz = sz[ip];
        const geom_float sss = ss[ip];

        // node モードの境界半割面 (ghost を含む面: ic0 か ic1 が nCells 以上) は拡散流束を加えない (skip)。
        // 根拠 (plan diffusion-node-boundary-real-distance.md §3 (c)):
        //  - Dirichlet (入口固定値 / k=0・ω ピン) では境界ノードはピンで上書きされ、半割面に何を足しても無意味。
        //    内部ノード I への拡散は内部双対面 W↔I が実距離で運ぶ (本ループの非境界 plane で計算済)。
        //  - Neumann (slip / zero-grad 出口) は物理的に ∂φ/∂n=0 ＝半割面フラックス 0 そのもの。
        //  → ghost mirror の dcc≈0 退化も、∇φ·S 弱形式の境界閉包依存も不要。skip が正しく最シンプル。
        // periodic 面は partner 実セルで両側 < nCells なので該当せず通常処理。(void)nWallHaloPlanes,dphi*。
        if (isNode != 0 && (ic0 >= nCells || ic1 >= nCells)) {
            return;
        }

        const flow_float dcc_x = ccx[ic1] - ccx[ic0];
        const flow_float dcc_y = ccy[ic1] - ccy[ic0];
        const flow_float dcc_z = ccz[ic1] - ccz[ic0];
        const flow_float dcc = sqrt(dcc_x * dcc_x + dcc_y * dcc_y + dcc_z * dcc_z);

        const flow_float denom = dcc_x * sxx + dcc_y * syy + dcc_z * szz;
        // ゼロ割ガードは相対値 (|d||S| の 1e-6) にする。旧 1.0e-12 は [m³] の絶対値で、3D の µm 級双対面
        // (|d·S| ~ 1e-15) では拡散コンダクタンスそのものを 1/10〜1/1000 に削っていた (case/16 3D 角線で ω が
        // 壁漸近解の 1/30 に落ち k 未減衰→角部加熱; codex レビュー 2026-09-08, plan turbulence-sst-node-corner-heating)。
        const flow_float denom_floor = static_cast<flow_float>(1.0e-6) * dcc * sss;
        const flow_float safe_denom = (abs(denom) < denom_floor) ? ((denom >= 0.0) ? denom_floor : -denom_floor) : denom;
        const flow_float delta = dcc * sss * sss / safe_denom;

        // sstSigmaBlend: σ = F1·σ1 + (1−F1)·σ2 をノードごとに (F1 は前 step の ransSource 値, 初期 1)
        // cell モードの境界 ghost (ic >= nCells) は F1 を持たないので内部側の F1 を使う (codex 2026-09-08)
        const flow_float F1a = (F1blend != nullptr) ? F1blend[(ic0 < nCells) ? ic0 : ic1] : static_cast<flow_float>(1.0);
        const flow_float F1b = (F1blend != nullptr) ? F1blend[(ic1 < nCells) ? ic1 : ic0] : static_cast<flow_float>(1.0);
        const flow_float sig0 = (F1blend != nullptr) ? (F1a * sigma + (static_cast<flow_float>(1.0) - F1a) * sigma2) : sigma;
        const flow_float sig1 = (F1blend != nullptr) ? (F1b * sigma + (static_cast<flow_float>(1.0) - F1b) * sigma2) : sigma;
        const flow_float mu0 = vis_lam[ic0] + sig0 * max(vis_turb[ic0], static_cast<flow_float>(0.0));
        const flow_float mu1 = vis_lam[ic1] + sig1 * max(vis_turb[ic1], static_cast<flow_float>(0.0));
        const flow_float mu_face = f * mu0 + (1.0 - f) * mu1;

        const flow_float dphi = phi[ic1] - phi[ic0];
        const flow_float flux = mu_face * (dphi / dcc) * delta;

        // point-implicit 拡散対角: -∂res/∂(ρφ) = (μ_face/ρ)·|δ|/dcc ≥ 0 [m³/s]。
        // 各セルは自身の ρ で正規化（φ=ρφ/ρ）。壁近傍の細セルで大きくなり陽的拡散 stiff 性を陰化する。
        const flow_float diag_face = mu_face * fabs(delta) / max(dcc, static_cast<flow_float>(1.0e-30));

        if (ic0 < nCells) {
            atomicAdd(&res_rho_phi[ic0], flux);
            atomicAdd(&transport_diag[ic0], diag_face / max(ro[ic0], static_cast<flow_float>(1.0e-30)));
        }
        if (ic1 < nCells) {
            atomicAdd(&res_rho_phi[ic1], -flux);
            atomicAdd(&transport_diag[ic1], diag_face / max(ro[ic1], static_cast<flow_float>(1.0e-30)));
        }
    }
}

__global__ void runge_kutta_exp_scalar_4th_d(
    int loop,
    flow_float coef_DT,
    flow_float coef_Res,
    flow_float* dt_local,
    geom_int nCells,
    geom_float* vol,
    flow_float* rho_phi,
    flow_float* rho_phi_N,
    flow_float* rho_phi_M,
    flow_float* res_rho_phi,
    flow_float* res_rho_phi_m)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        const flow_float dt_l = dt_local[ic];
        const geom_float v = vol[ic];

        if (loop == 0) {
            res_rho_phi_m[ic] = 0.0;
        }

        res_rho_phi_m[ic] += coef_Res * res_rho_phi[ic] * dt_l / v;

        if (loop < 3) {
            rho_phi[ic] = rho_phi_N[ic] + coef_DT * res_rho_phi[ic] * dt_l / v;
        } else {
            res_rho_phi[ic] = res_rho_phi_m[ic] * v / dt_l;
            rho_phi[ic] = rho_phi_N[ic] + res_rho_phi_m[ic];
        }
    }
}

__global__ void runge_kutta_exp_scalar_d(
    flow_float coef_N,
    flow_float coef_M,
    flow_float coef_Res,
    flow_float* dt_local,
    geom_int nCells,
    geom_float* vol,
    flow_float* rho_phi,
    flow_float* rho_phi_N,
    flow_float* rho_phi_M,
    flow_float* res_rho_phi,
    flow_float* src_jac,
    flow_float* transport_diag,
    flow_float floor)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        const flow_float dt_l = dt_local[ic];
        const geom_float v = vol[ic];
        // 源項（消散）+ 輸送項（移流+拡散）の stiff 性を point-implicit で減衰。
        //   src_jac=∂D/∂(ρφ)≥0、transport_diag=Σ_f[max(±ṁ,0)+μ_face|δ|/dcc]/ρ≥0 [m³/s]（/v で 1/s 化）。
        // 源項・輸送対角が双方 0 なら従来の純陽的更新と一致。減衰係数は陰解法対角
        // D_φ=V/Δτ+V·src_jac+transport_diag に Δτ/V を掛けた形と整合。理論は methods/time_integration。
        const flow_float fac = static_cast<flow_float>(1.0)
            + coef_Res * dt_l * (src_jac[ic] + transport_diag[ic] / v);
        const flow_float updated = coef_N * rho_phi_N[ic] + coef_M * rho_phi_M[ic]
                    + (coef_Res * res_rho_phi[ic] * dt_l / v) / fac;
        // realizability 下限。source 側 point-implicit と整合。下限に達しない範囲では無影響。
        rho_phi[ic] = max(updated, floor);
    }
}

}

void scalarTransportResidual_d(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                               const ScalarTransportDesc& desc)
{
    dim3 dimGrid_normal_halo = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));

    scalar_advection_first_order_d<<<dimGrid_normal_halo , cuda_cfg.dimBlock>>>(
        msh.nCells,
        msh.nNormal_halo_Planes,
        (cfg.discretization == "node") ? 1 : 0,
        msh.normal_halo_planes_d,
        msh.map_plane_cells_d,
        desc.phi,
        var.c_d["ro"],
        var.p_d["massflux"],
        desc.res_rho_phi,
        desc.transport_diag);

    if (cfg.scalarDiffusion == 1 && desc.diffusion == 1) {
        scalar_diffusion_first_order_d<<<dimGrid_normal_halo , cuda_cfg.dimBlock>>>(
            msh.nCells,
            msh.nNormal_halo_Planes,
            (cfg.discretization == "node") ? msh.nWallHaloPlanes : 0,
            (cfg.discretization == "node") ? 1 : 0,
            msh.normal_halo_planes_d,
            msh.map_plane_cells_d,
            var.c_d["ccx"],
            var.c_d["ccy"],
            var.c_d["ccz"],
            var.p_d["fx"],
            var.p_d["sx"],
            var.p_d["sy"],
            var.p_d["sz"],
            var.p_d["ss"],
            desc.phi,
            desc.dphidx,
            desc.dphidy,
            desc.dphidz,
            var.c_d["ro"],
            var.c_d["vis_lam"],
            var.c_d["vis_turb"],
            desc.sigma,
            desc.res_rho_phi,
            desc.transport_diag,
            desc.sigma2,
            desc.F1);
    }
}

void scalarTimeIntegration_d(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                             const ScalarTransportDesc& desc)
{
    if (cfg.timeIntegration == 4) {
        runge_kutta_exp_scalar_4th_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
            loop,
            cfg.coef_DT_4thRunge[loop],
            cfg.coef_Res_4thRunge[loop],
            var.c_d["dt_local"],
            msh.nCells,
            var.c_d["volume"],
            desc.rho_phi,
            desc.rho_phi_N,
            desc.rho_phi_M,
            desc.res_rho_phi,
            desc.res_rho_phi_m);
    } else if (cfg.timeIntegration == 1 || cfg.timeIntegration == 3) {
        runge_kutta_exp_scalar_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
            cfg.coef_N[loop],
            cfg.coef_M[loop],
            cfg.coef_Res[loop],
            var.c_d["dt_local"],
            msh.nCells,
            var.c_d["volume"],
            desc.rho_phi,
            desc.rho_phi_N,
            desc.rho_phi_M,
            desc.res_rho_phi,
            desc.src_jac,
            desc.transport_diag,
            desc.floor);
    } else if (cfg.timeIntegration == 11) {
        // 陰解法 (block-DPLUR) ステップでの化学種更新。1 回の point-implicit forward-Euler:
        //   ρφ = ρφ_N + (res·Δτ/V) / (1 + Δτ(src_jac + transport_diag/V))。
        // coef_N=1, coef_M=0, coef_Res=1 (loop 無関係)。点陰的対角が高 CFL_pseudo を安定化し、
        // Δτ→∞ で局所定常増分 res/diag に収束 (SST point-implicit と同形)。
        runge_kutta_exp_scalar_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
            static_cast<flow_float>(1.0),
            static_cast<flow_float>(0.0),
            static_cast<flow_float>(1.0),
            var.c_d["dt_local"],
            msh.nCells,
            var.c_d["volume"],
            desc.rho_phi,
            desc.rho_phi_N,
            desc.rho_phi_M,
            desc.res_rho_phi,
            desc.src_jac,
            desc.transport_diag,
            desc.floor);
    }
}
