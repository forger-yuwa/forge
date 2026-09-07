#include "ransBoundary_d.cuh"
#include "ransWallFunction_d.cuh"

#include "cuda_forge/cudaWrapper.cuh"

namespace {

constexpr flow_float kSstBeta1  = static_cast<flow_float>(0.075);
constexpr flow_float kSstBetaSt = static_cast<flow_float>(0.09);   // β*
constexpr flow_float kKappa     = static_cast<flow_float>(0.41);   // von Karman 定数
constexpr flow_float kOmegaMin  = static_cast<flow_float>(1.0e-10);
constexpr flow_float kSmall = static_cast<flow_float>(1.0e-12);

// SST 壁面 ω 境界。
// wallTreatment==0: low-Re 壁解像型 ω_w = 60ν/(β₁ y²) (methods/turbulence §6.1)。
// wallTreatment==1: Menter automatic ブレンド ω_w = √(ω_vis² + ω_log²) (§6.5 (b))。
//   ω_vis = 6ν/(β₁ y²), ω_log = u_τ/(√β* κ y)。u_τ は事前に utau に格納済み。
__global__ void rans_wall_scalar_boundary_d(
    geom_int nb,
    geom_int* bplane_cell,
    geom_int* bplane_cell_ghst,
    flow_float* ro,
    flow_float* vis_lam,
    flow_float* wall_dist,
    flow_float* k,
    flow_float* omega,
    flow_float* kb,
    flow_float* omegab,
    int wallTreatment,
    flow_float* utau,
    flow_float* roOmega,
    int isNode)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;

    if (ib < nb) {
        const geom_int ic = bplane_cell[ib];
        const geom_int ig = bplane_cell_ghst[ib];

        const flow_float rho_c = max(ro[ic], kSmall);
        const flow_float nu_c = vis_lam[ic] / rho_c;
        const flow_float y_w = max(wall_dist[ic], kSmall);

        if (wallTreatment == 1) {
            // automatic wall treatment (methods/turbulence §6.5):
            //   ω: Menter ブレンド ω_w=√(ω_vis²+ω_log²) を **wall-adjacent セル値にピン留め**する。
            //      ω_w は全 y⁺ で妥当 (y⁺→0 で ω_vis 支配=壁解像値) なので μ_t=ρa₁k/max(a₁ω,SF₂) を
            //      適切に上限し、k runaway を断つ。conserved roOmega も合わせて更新する。
            //   k: zero-gradient (Neumann)。近壁 P_k は wall-function 値 (ransSource の wf_pk) に
            //      置換済みなので、k は平衡値 u_τ²/√β* に自律収束する (b' の runaway は起きない)。
            const flow_float omega_vis = static_cast<flow_float>(6.0) * nu_c / (kSstBeta1 * y_w * y_w);
            const flow_float omega_log = utau[ib] / (sqrt(kSstBetaSt) * kKappa * y_w);
            const flow_float omega_w   = max(sqrt(omega_vis * omega_vis + omega_log * omega_log), kOmegaMin);
            omegab[ib]  = omega_w;
            omega[ic]   = omega_w;                       // wall-adjacent セルにピン留め
            roOmega[ic] = rho_c * omega_w;               // conserved も整合
            kb[ib] = k[ic];
            k[ig]     = k[ic];                           // ∂k/∂n = 0 (zero-gradient)
            omega[ig] = static_cast<flow_float>(2.0) * omega_w - omega[ic];
        } else {
            // wall-resolved (low-Re): k_w = 0 Dirichlet, ω_w = 60ν/(β₁ y²)。
            kb[ib] = static_cast<flow_float>(0.0);
            omegab[ib] = max(static_cast<flow_float>(60.0) * nu_c / (kSstBeta1 * y_w * y_w), kOmegaMin);
            k[ig]     = static_cast<flow_float>(2.0) * kb[ib] - k[ic];
            omega[ig] = static_cast<flow_float>(2.0) * omegab[ib] - omega[ic];
            // node-centered: 壁ノードが CV 中心 (ic) なので、ゴーストミラー (dcc≈0 退化) では ω[ic] が
            // ω_w にピン留めされない → ω が過小 → k 過小消散 → 過剰乱流。ω[ic] を直接 ω_w に固定し、
            // 保存量 roOmega も ρ·ω_w に整合させる (point-implicit で dω=0 decouple して保持)。cell は
            // ゴーストで正しくピンされるので変更しない。
            if (isNode != 0) {
                omega[ic]   = omegab[ib];
                roOmega[ic] = rho_c * omegab[ib];
            }
        }
    }
}

__global__ void rans_dirichlet_scalar_boundary_d(
    geom_int nb,
    geom_int* bplane_cell,
    geom_int* bplane_cell_ghst,
    flow_float* ro,
    flow_float* k,
    flow_float* omega,
    flow_float* roK,
    flow_float* roOmega,
    flow_float* scalarDirichletPin,
    flow_float* kb,
    flow_float* omegab,
    int isNode)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;

    if (ib < nb) {
        const geom_int ic = bplane_cell[ib];
        const geom_int ig = bplane_cell_ghst[ib];

        k[ig] = static_cast<flow_float>(2.0) * kb[ib] - k[ic];
        omega[ig] = static_cast<flow_float>(2.0) * omegab[ib] - omega[ic];

        // node-centered: 境界面の流束/勾配 (ransTransport) は ghost を使わず境界ノード値 k[ic0] を
        // そのまま境界値に使う設計。壁は omega[ic]=omega_w をピンするが、入口 Dirichlet が ghost しか
        // 書かないと境界ノードがピンされず、入口 k/ω が freestream から drift する (実測 1000→918)。
        // node では境界ノードを Dirichlet 値に直接ピンする (cell は ghost 経由で正しく課されるので不変)。
        // primitive だけでなく保存量 roK/roOmega も ρ·Dirichlet 値に整合させ (放置すると roOmega/ro vs omega が
        // 最大 888727 不整合し残差が rms を汚染、入口×壁コーナーで最大化)、scalarDirichletPin を立てて
        // ransSource で残差除外する。詳細: plans/active/turbulence-node-inlet-dirichlet-conserved.md。
        if (isNode != 0) {
            k[ic]       = kb[ib];
            omega[ic]   = omegab[ib];
            roK[ic]     = ro[ic] * kb[ib];
            roOmega[ic] = ro[ic] * omegab[ib];
            scalarDirichletPin[ic] = static_cast<flow_float>(1.0);
        }
    }
}

__global__ void rans_neumann_scalar_boundary_d(
    geom_int nb,
    geom_int* bplane_cell,
    geom_int* bplane_cell_ghst,
    flow_float* k,
    flow_float* omega,
    flow_float* kb,
    flow_float* omegab)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;

    if (ib < nb) {
        const geom_int ic = bplane_cell[ib];
        const geom_int ig = bplane_cell_ghst[ib];

        k[ig] = k[ic];
        omega[ig] = omega[ic];
        kb[ib] = k[ic];
        omegab[ib] = omega[ic];
    }
}

}

// node-centered: 壁ノードは壁面上 (wall_dist=0) なので SST 壁 omega BC が overflow する。
// 各壁ノードの実効 y を「内部隣接ノード (壁面ノードを除く) の wall_dist 平均」に置換する。
// cell モード・非壁ノードは wall_dist をそのままコピーする (挙動不変)。
__global__ void compute_wall_y_eff_d(
    geom_int nCells, geom_int* wall_flag, flow_float* wall_dist,
    geom_int* cell_planes_index, geom_int* cell_planes, geom_int* plane_cells,
    flow_float* wall_y_eff)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        if (wall_flag != nullptr && wall_flag[ic] == 1) {
            // 最近接の内部隣接ノードの wall_dist = 壁→第1オフ壁点の法線距離 Δy。ω_w=60ν/(β₁Δy²) の正しい値。
            // (平均は Δy を 1 桁過大評価し ω を ~100x 過小→k 過小消散→過剰乱流。最近接が物理的に正しい。)
            flow_float ymin = static_cast<flow_float>(1.0e30);
            int cnt = 0;
            for (geom_int j = cell_planes_index[ic]; j < cell_planes_index[ic + 1]; ++j) {
                const geom_int ip = cell_planes[j];
                const geom_int a = plane_cells[2 * ip + 0];
                const geom_int b = plane_cells[2 * ip + 1];
                const geom_int other = (a == ic) ? b : a;
                // wall_dist==0 の非壁ノード (slip 壁が wallDistExtraPhysIDs で壁点集合に入っている場合の
                // slip 壁ノード) は第一オフ壁点ではないので除外する (2026-09-08: no-slip→slip 接合の壁ノードで
                // y_eff=kSmall→ω=1e23→NaN になった, case/16 run_0225)
                if (other < nCells && wall_flag[other] == 0 && wall_dist[other] > kSmall) {
                    ymin = min(ymin, max(wall_dist[other], kSmall));
                    ++cnt;
                }
            }
            // 凹コーナー (段差リップ・2壁交線) では全エッジ隣接が壁ノードで cnt==0 になる。
            // 構造化要素では対角の内部ノードへ直接エッジが無いため、2-ring 探索 (隣接壁ノードの隣接 =
            // 対角内部ノード) で第一オフ壁距離を拾う。これをしないと wall_dist=0 のコーナー (3D 段差) で
            // y=kSmall→ω=60ν/(β₁kSmall²) が float overflow し ω=inf になる (2D は偶々 wall_dist>0 で救われていた)。
            if (cnt == 0) {
                for (geom_int j = cell_planes_index[ic]; j < cell_planes_index[ic + 1]; ++j) {
                    const geom_int ip = cell_planes[j];
                    const geom_int a = plane_cells[2 * ip + 0];
                    const geom_int b = plane_cells[2 * ip + 1];
                    const geom_int nb = (a == ic) ? b : a;          // 隣接壁ノード
                    if (nb >= nCells || wall_flag[nb] != 1) continue;
                    for (geom_int j2 = cell_planes_index[nb]; j2 < cell_planes_index[nb + 1]; ++j2) {
                        const geom_int ip2 = cell_planes[j2];
                        const geom_int a2 = plane_cells[2 * ip2 + 0];
                        const geom_int b2 = plane_cells[2 * ip2 + 1];
                        const geom_int other2 = (a2 == nb) ? b2 : a2;
                        if (other2 < nCells && other2 != ic && wall_flag[other2] == 0 && wall_dist[other2] > kSmall) {
                            ymin = min(ymin, max(wall_dist[other2], kSmall));
                            ++cnt;
                        }
                    }
                }
            }
            wall_y_eff[ic] = (cnt > 0) ? ymin : max(wall_dist[ic], kSmall);
        } else {
            wall_y_eff[ic] = wall_dist[ic];
        }
    }
}

void ransBoundary_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var)
{
    if (!(cfg.LESorRANS == 2 && cfg.RANSmodel == 1)) {
        return;
    }

    if (bc.iPlanes.empty()) {
        return;
    }

    if (bc.bcondKind == "wall" || bc.bcondKind == "wall_isothermal") {
        // automatic wall treatment (wallTreatmentSST==1) では先に摩擦速度 u_τ を Reichardt 逆解きで
        // 求め bc.bvar_d["utau"] に格納する。ω 壁 BC と viscousFlux 壁せん断が同じ u_τ を共有する。
        computeWallFrictionSST_d_wrapper(cfg , cuda_cfg , bc , msh , var);

        // node-centered: 壁ノードの y=wall_dist=0 で omega 壁 BC が overflow するため、内部隣接ノードの
        // wall_dist 平均 (wall_y_eff) を omega BC の y に使う。cell モードは wall_dist をそのまま使う。
        flow_float* y_for_omega = var.c_d["wall_dist"];
        if (cfg.discretization == "node" && msh.wall_flag_d != nullptr) {
            compute_wall_y_eff_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
                msh.nCells, msh.wall_flag_d, var.c_d["wall_dist"],
                msh.map_cell_planes_index_d, msh.map_cell_planes_d, msh.map_plane_cells_d,
                var.c_d["wall_y_eff"]);
            y_for_omega = var.c_d["wall_y_eff"];
        }

        rans_wall_scalar_boundary_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>(
            static_cast<geom_int>(bc.iPlanes.size()),
            bc.map_bplane_cell_d,
            bc.map_bplane_cell_ghst_d,
            var.c_d["ro"],
            var.c_d["vis_lam"],
            y_for_omega,
            var.c_d["k"],
            var.c_d["omega"],
            bc.bvar_d["kb"],
            bc.bvar_d["omegab"],
            cfg.wallTreatmentSST,
            bc.bvar_d["utau"],
            var.c_d["roOmega"],
            (cfg.discretization == "node") ? 1 : 0);
        return;
    }

    if (bc.bcondKind == "inlet_uniformVelocity" ||
        bc.bcondKind == "inlet_fluctVelocity" ||
        bc.bcondKind == "inlet_Pressure" ||
        bc.bcondKind == "inlet_Pressure_dir") {
        rans_dirichlet_scalar_boundary_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>(
            static_cast<geom_int>(bc.iPlanes.size()),
            bc.map_bplane_cell_d,
            bc.map_bplane_cell_ghst_d,
            var.c_d["ro"],
            var.c_d["k"],
            var.c_d["omega"],
            var.c_d["roK"],
            var.c_d["roOmega"],
            var.c_d["scalarDirichletPin"],
            bc.bvar_d["k"],
            bc.bvar_d["omega"],
            (cfg.discretization == "node") ? 1 : 0);
        return;
    }

    if (bc.bcondKind == "slip" ||
        bc.bcondKind == "axis" ||
        bc.bcondKind == "outlet_statPress" ||
        bc.bcondKind == "outflow" ||
        bc.bcondKind == "periodic") {
        rans_neumann_scalar_boundary_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>(
            static_cast<geom_int>(bc.iPlanes.size()),
            bc.map_bplane_cell_d,
            bc.map_bplane_cell_ghst_d,
            var.c_d["k"],
            var.c_d["omega"],
            bc.bvar_d["kb"],
            bc.bvar_d["omegab"]);
        return;
    }
}