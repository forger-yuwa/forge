#include "wmlesWallModel_d.cuh"

#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/wallLaw_d.cuh"
#include "cuda_forge/thermo_d.cuh"
#include "cuda_forge/speciesTransport_d.cuh"  // species_roY_device_ptr()
#include "cuda_forge/nodeWallDirichlet_d.cuh"  // pinWallNodeTemperature_bcond (状態層・共有)

#include <unordered_set>

// WMLES 代数壁応力モデル (theory §10)。壁境界面並列で
//   取得層: マッチング点 (cell=第一セル / node=Normal_Neighbor 内点) と壁面状態 (T_w, p_w)
//   壁モデル: Reichardt warm start Newton → u_τ, τ_w / Kader → q_w (等温壁のみ)
//   書き込み: bvar (utau/ypls/twall_*/qwall), node は Tau_Wall/Qw_Wall (per-CV)
// を行う。流束への適用は viscousFlux_d.cu 側 (cell: wallTreatment==2 置換 / node: AddTauWall+AddQWall)。

namespace {

constexpr flow_float kSmall = static_cast<flow_float>(1.0e-12);

// 壁温 T_w での層流物性 (μ, λ, cp, R)。gasProperties_d (セル評価) と同式で壁温評価する。
//   viscMethod 0: 定数 / 1: Sutherland (熱伝導は定数) / 2: kinetic theory (Chapman-Enskog+Wilke)
//   thermalMethod 2 (TP): NASA 多項式 (組成 Y は呼び出し側で正規化済み)。それ以外 (CPG): 定数 γ, cp。
__device__ inline void wmles_wall_props(
    int thermalMethod, int viscMethod,
    flow_float ga, flow_float cp_const, flow_float visc_const, flow_float thermCond_const,
    const SpeciesThermo* sp, int nSpecies, const double* Y,
    flow_float Tw,
    flow_float& mu_w, flow_float& lam_w, flow_float& cp_w, flow_float& R_w)
{
    const int n = (nSpecies >= 1) ? nSpecies : 1;
    if (thermalMethod == 2) {
        const double Twd = (double)Tw;
        R_w  = (flow_float)((n > 1) ? thermo_R_mix(sp, n, Y) : thermo_R_species(sp[0]));
        cp_w = (flow_float)((n > 1) ? thermo_cp_mix(sp, n, Y, Twd) : thermo_cp_mass(sp[0], Twd));
    } else {
        R_w  = cp_const - cp_const/ga;   // wall_isothermal_d と同式
        cp_w = cp_const;
    }
    if (viscMethod == 0) {
        mu_w  = visc_const;
        lam_w = thermCond_const;
    } else if (viscMethod == 1) {   // Sutherland (gasProperties_d と同定数)。熱伝導は定数
        const flow_float T0  = 273.0;
        const flow_float mu0 = 1.716e-5;
        const flow_float Smu = 111.0;
        mu_w  = mu0*pow(Tw/T0, static_cast<flow_float>(3.0/2.0))*(T0+Smu)/(Tw+Smu);
        lam_w = thermCond_const;
    } else {                        // kinetic theory
        double X[THERMO_MAX_SPECIES];
        thermo_X_from_Y(sp, n, Y, X);
        const double Twd = (double)Tw;
        mu_w  = (flow_float)thermo_mu_mix(sp, n, X, Twd);
        lam_w = (flow_float)thermo_lambda_mix(sp, n, X, Twd);
    }
}

// node 用: Tau_Wall / Qw_Wall を全 CV -1 (inactive) に初期化 (毎 step、bcond ループ前)。
__global__ void wmles_init_wall_arrays_d(geom_int nCells, flow_float* Tau_Wall, flow_float* Qw_Wall)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        Tau_Wall[ic] = static_cast<flow_float>(-1.0);
        Qw_Wall[ic]  = static_cast<flow_float>(-1.0);
    }
}

__global__ void wmles_wall_model_d(
    geom_int nb,
    geom_int* bplane_plane,
    geom_int* bplane_cell,
    geom_int* bplane_cell_ghst,
    geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,
    // 取得層 (node): SU2 Normal_Neighbor 内点選択用 CSR (ransWallFunction_d と同機構)
    int isNode,
    geom_int nNormalPlanes,
    geom_int* cell_planes_index, geom_int* cell_planes, geom_int* plane_cells,
    geom_int* wall_flag,
    geom_float* ccx, geom_float* ccy, geom_float* ccz,
    flow_float* wall_dist,
    // 場
    flow_float* ro, flow_float* Ux, flow_float* Uy, flow_float* Uz,
    flow_float* T, flow_float* P,
    // 物性 config
    int thermalMethod, flow_float ga, flow_float cp_const,
    const SpeciesThermo* sp, flow_float* const* roY, int nSpecies,
    int viscMethod, flow_float visc_const, flow_float thermCond_const,
    // モデルパラメータ
    flow_float tol, int maxIt,
    // 熱条件: 0=断熱 (q_w=0) / 1=等温 (Tsb=指定壁温)
    int isothermal, flow_float* Tsb,
    // 出力 (bvar)。utau_b は warm start を兼ねる (前 step 値を初期値に使う)
    flow_float* utau_b, flow_float* ypls_b,
    flow_float* twall_x_b, flow_float* twall_y_b, flow_float* twall_z_b,
    flow_float* qwall_b,
    // 出力 (node のみ): 壁ノードの τ_w / q_w。viscousFlux_d の AddTauWall / AddQWall が消費
    flow_float* Tau_Wall, flow_float* Qw_Wall,
    // Newton 不収束カウンタ (デバッグ。恒常的に増えるならメッシュ/場の異常)
    unsigned int* nonconvCount)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib >= nb) return;

    const geom_int ip = bplane_plane[ib];
    const geom_int ic = bplane_cell[ib];
    const geom_int ig = bplane_cell_ghst[ib];

    // 壁単位法線 (流体→壁の外向き)
    const flow_float sss = max(ss[ip], kSmall);
    const flow_float nx = sx[ip] / sss;
    const flow_float ny = sy[ip] / sss;
    const flow_float nz = sz[ip] / sss;

    // ---- 取得層: マッチング点 irep と壁面垂直距離 d、壁面状態 (T_w_in, p_w) ----
    //   cell: irep=第一セル, d=wall_dist, 壁面値=ghost との面平均 (断熱ミラーでは内点値)。
    //   node: irep=Normal_Neighbor 内点 (壁ノードは u=0 Dirichlet で退化するため),
    //         d=(x_I-x_W)·(-n), 壁面値=壁ノードの解そのもの。
    geom_int   irep = ic;
    flow_float d    = max(wall_dist[ic], kSmall);
    flow_float Tw_in, p_w;
    if (isNode != 0) {
        flow_float best_cos = static_cast<flow_float>(-2.0);
        geom_int   bestI = -1;
        flow_float bestDn = kSmall;
        const flow_float xw = ccx[ic], yw = ccy[ic], zw = ccz[ic];
        for (geom_int j = cell_planes_index[ic]; j < cell_planes_index[ic + 1]; ++j) {
            const geom_int ipn = cell_planes[j];
            if (ipn >= nNormalPlanes) continue;                 // 内部双対面のみ
            const geom_int a = plane_cells[2 * ipn + 0];
            const geom_int b = plane_cells[2 * ipn + 1];
            const geom_int cand = (a == ic) ? b : a;
            if (wall_flag[cand] != 0) continue;                 // 内部ノードのみ
            const flow_float dx = ccx[cand] - xw;
            const flow_float dy = ccy[cand] - yw;
            const flow_float dz = ccz[cand] - zw;
            const flow_float dist = sqrt(dx * dx + dy * dy + dz * dz);
            if (dist <= kSmall) continue;
            const flow_float dn_in = -(dx * nx + dy * ny + dz * nz);
            if (dn_in <= static_cast<flow_float>(0.0)) continue;
            const flow_float cosv = dn_in / dist;
            if (cosv > best_cos) { best_cos = cosv; bestI = cand; bestDn = dn_in; }
        }
        if (bestI < 0) {   // 代表点なし: 層流退避 (τ_w=0 扱い)
            utau_b[ib] = 0.0; ypls_b[ib] = 0.0;
            twall_x_b[ib] = 0.0; twall_y_b[ib] = 0.0; twall_z_b[ib] = 0.0;
            qwall_b[ib] = 0.0;
            Tau_Wall[ic] = 0.0;
            if (isothermal != 0) Qw_Wall[ic] = 0.0;
            return;
        }
        irep = bestI;
        d    = max(bestDn, kSmall);
        Tw_in = T[ic];
        p_w   = P[ic];
    } else {
        Tw_in = static_cast<flow_float>(0.5)*(T[ig] + T[ic]);
        p_w   = static_cast<flow_float>(0.5)*(P[ig] + P[ic]);
    }
    const flow_float Tw = (isothermal != 0) ? Tsb[ib] : Tw_in;   // 等温壁は指定壁温を優先

    // ---- 壁面物性 (T_w 評価; γ/cp 一定を焼き込まない) ----
    //   組成: node は壁ノード, cell は第一セル (ghost 組成=内部と同一)。
    double Y[THERMO_MAX_SPECIES];
    if (!thermo_cell_Y(roY, nSpecies, ic, Y)) Y[0] = 1.0;   // 共通 helper (thermo_d.cuh)
    flow_float mu_w, lam_w, cp_w, R_w;
    wmles_wall_props(thermalMethod, viscMethod, ga, cp_const, visc_const, thermCond_const,
                     sp, nSpecies, Y, Tw, mu_w, lam_w, cp_w, R_w);
    const flow_float rho_w = p_w / (R_w * Tw);
    const flow_float Pr    = mu_w * cp_w / max(lam_w, kSmall);

    // ---- マッチング点の壁平行速度 ----
    const flow_float uc = Ux[irep], vc = Uy[irep], wc = Uz[irep];
    const flow_float un = uc*nx + vc*ny + wc*nz;
    const flow_float upx = uc - un*nx;
    const flow_float upy = vc - un*ny;
    const flow_float upz = wc - un*nz;
    const flow_float upar = sqrt(upx*upx + upy*upy + upz*upz);
    const flow_float Tm   = T[irep];

    if (upar <= kSmall) {
        // よどみ: Newton スキップ、層流極限 (τ_w→0)。熱は純伝導。
        utau_b[ib] = 0.0; ypls_b[ib] = 0.0;
        twall_x_b[ib] = 0.0; twall_y_b[ib] = 0.0; twall_z_b[ib] = 0.0;
        const flow_float qw = (isothermal != 0) ? lam_w*(Tw - Tm)/d : static_cast<flow_float>(0.0);
        qwall_b[ib] = qw;
        if (isNode != 0) {
            Tau_Wall[ic] = 0.0;
            if (isothermal != 0) Qw_Wall[ic] = qw;
        }
        return;
    }

    // ---- u_τ: warm start Newton (wallLaw_solve_utau, theory §10.2) ----
    flow_float utau = utau_b[ib];   // 前 step 値 (初回/無効値はソルバ内で層流推定に置換)
    const bool converged = wallLaw_solve_utau(upar, d, rho_w, mu_w, tol, maxIt, utau);
    flow_float tauw;
    if (converged) {
        tauw = rho_w * utau * utau;
    } else {
        // 不収束: 層流応力にフォールバック (黙って握りつぶさない: カウンタ++)
        tauw = mu_w * upar / d;
        utau = sqrt(tauw / rho_w);
        if (nonconvCount != nullptr) atomicAdd(nonconvCount, 1u);
    }
    const flow_float yp = rho_w * utau * d / mu_w;
    const flow_float up = upar / max(utau, kSmall);

    // ---- τ_w ベクトル (流束/面積, 流体を減速させる向き -ê_∥) ----
    const flow_float inv_up = static_cast<flow_float>(1.0) / upar;
    twall_x_b[ib] = -tauw * upx * inv_up;
    twall_y_b[ib] = -tauw * upy * inv_up;
    twall_z_b[ib] = -tauw * upz * inv_up;

    // ---- q_w (等温壁のみ; 壁→流体を正, theory §10.3) ----
    //   駆動温度差は回復温度 T_r = T_m + Pr^{1/3} u_∥²/(2 cp_w)。
    //   y⁺ が十分小さいときは Kader の層流極限 (T⁺→Pr y⁺) と等価な純伝導式で評価 (0/0 回避)。
    flow_float qw = 0.0;
    if (isothermal != 0) {
        const flow_float Tr = Tm + cbrt(Pr) * upar * upar / (static_cast<flow_float>(2.0) * cp_w);
        if (yp < static_cast<flow_float>(1.0e-1)) {
            qw = lam_w * (Tw - Tr) / d;
        } else {
            const flow_float Tp = wallLaw_kader_tplus(Pr, yp);
            qw = rho_w * cp_w * utau * (Tw - Tr) / max(Tp, kSmall);
        }
    }
    qwall_b[ib] = qw;

    utau_b[ib] = utau;
    ypls_b[ib] = yp;

    if (isNode != 0) {
        Tau_Wall[ic] = tauw;                       // AddTauWall (W↔I 接線 traction 再スケール)
        if (isothermal != 0) Qw_Wall[ic] = qw;     // AddQWall (W↔I 熱流束置換) + res_roe ピンマーカ
    }
}


} // namespace

bool wmlesActiveForBcond(const solverConfig& cfg, const bcond& bc)
{
    if (cfg.LESorRANS == 2) return false;   // RANS は SST 壁関数側 (共存しない)
    if (!(bc.bcondKind == "wall" || bc.bcondKind == "wall_isothermal")) return false;
    const auto it = bc.inputInts.find("wallModelLES");
    return (it != bc.inputInts.end() && it->second == 1);
}

bool wmlesAnyActive(const solverConfig& cfg, const mesh& msh)
{
    for (const auto& bc : msh.bconds) if (wmlesActiveForBcond(cfg, bc)) return true;
    return false;
}

bool wmlesNodeActive(const solverConfig& cfg, const mesh& msh)
{
    return cfg.discretization == "node" && msh.wall_flag_d != nullptr && wmlesAnyActive(cfg, msh);
}

bool wmlesNodeIsothermalActive(const solverConfig& cfg, const mesh& msh)
{
    if (cfg.discretization != "node" || msh.wall_flag_d == nullptr) return false;
    for (const auto& bc : msh.bconds)
        if (wmlesActiveForBcond(cfg, bc) && bc.bcondKind == "wall_isothermal") return true;
    return false;
}

void applyWmlesWallModel(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (!wmlesAnyActive(cfg, msh)) return;

    const int isNode = (cfg.discretization == "node" && msh.wall_flag_d != nullptr) ? 1 : 0;

    // Newton 不収束カウンタ (一度だけ確保)
    static unsigned int* nonconv_d = nullptr;
    if (nonconv_d == nullptr) {
        gpuErrchk(cudaMalloc(&nonconv_d, sizeof(unsigned int)));
        gpuErrchk(cudaMemset(nonconv_d, 0, sizeof(unsigned int)));
    }

    // warm start 用 utau bvar は type-0 (未初期化 device メモリ) のため、初回のみ 0 埋めする
    // (0 → ソルバ内で層流推定に置換。ゴミ値からの Newton 発散を防ぐ)。
    static std::unordered_set<const void*> utauInitialized;

    // node: Tau_Wall / Qw_Wall を毎 step -1 (inactive) 初期化 (SST の init_wf_pk_d と同位相)
    if (isNode != 0) {
        wmles_init_wall_arrays_d<<<cuda_cfg.dimGrid_normalcell , cuda_cfg.dimBlock>>>(
            msh.nCells, var.c_d["Tau_Wall"], var.c_d["Qw_Wall"]);
    }

    for (auto& bc : msh.bconds) {
        if (!wmlesActiveForBcond(cfg, bc)) continue;
        if (bc.iPlanes.empty()) continue;

        if (utauInitialized.insert(bc.bvar_d["utau"]).second) {
            gpuErrchk(cudaMemset(bc.bvar_d["utau"], 0, bc.iPlanes.size()*sizeof(flow_float)));
        }

        const int isothermal = (bc.bcondKind == "wall_isothermal") ? 1 : 0;

        // node 等温壁: 壁ノード温度状態を指定壁温にピン (流束評価の前・状態層の共有関数)
        if (isNode != 0 && isothermal != 0) {
            pinWallNodeTemperature_bcond(cfg, cuda_cfg, bc, var);
        }

        wmles_wall_model_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>(
            static_cast<geom_int>(bc.iPlanes.size()),
            bc.map_bplane_plane_d,
            bc.map_bplane_cell_d,
            bc.map_bplane_cell_ghst_d,
            var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
            isNode,
            msh.nNormalPlanes,
            msh.map_cell_planes_index_d, msh.map_cell_planes_d, msh.map_plane_cells_d,
            msh.wall_flag_d,
            var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
            var.c_d["wall_dist"],
            var.c_d["ro"], var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"],
            var.c_d["T"], var.c_d["P"],
            cfg.thermalMethod, cfg.gamma, cfg.cp,
            thermo_species_device_ptr(), species_roY_device_ptr(), cfg.nSpecies,
            cfg.viscMethod, cfg.visc, cfg.thermCond,
            cfg.wmlesNewtonTol, cfg.wmlesNewtonMaxIt,
            isothermal, bc.bvar_d["Ts"],
            bc.bvar_d["utau"], bc.bvar_d["ypls"],
            bc.bvar_d["twall_x"], bc.bvar_d["twall_y"], bc.bvar_d["twall_z"],
            bc.bvar_d["qwall"],
            var.c_d["Tau_Wall"], var.c_d["Qw_Wall"],
            nonconv_d);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

    // 不収束カウンタの定期報告 (恒常的に出る場合はメッシュ/スキームの問題: theory §10.2)
    static long callCount = 0;
    static unsigned int lastReported = 0;
    callCount++;
    if (callCount % 1000 == 0) {
        unsigned int h = 0;
        gpuErrchk(cudaMemcpy(&h, nonconv_d, sizeof(unsigned int), cudaMemcpyDeviceToHost));
        if (h != lastReported) {
            printf("[WMLES] u_tau Newton non-converged (laminar fallback) count: %u (cumulative)\n", h);
            lastReported = h;
        }
    }
}
