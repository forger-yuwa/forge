#include <cstdlib>
#include "convectiveFlux_d.cuh"
#include "lowMachPrecond_d.cuh"
#include "speciesTransport_d.cuh"  // species_roY_device_ptr()
#include "condensationProperties_d.cuh"  // n2_latent (二相エネルギー流束の潜熱補正)
#include "convectiveFlux_common_d.cuh"

#include <stdexcept>
#include "convectiveFlux_slau_d.inc.cuh"

#include "legacy/convectiveFlux_ausm_keep_d.inc.cuh"



#include "convectiveFlux_hlle_d.inc.cuh"



#include "convectiveFlux_roe_d.inc.cuh"

#include "convectiveFlux_keep_d.inc.cuh"

#include "legacy/convectiveFlux_keepslau_d.inc.cuh"

#include "convectiveFlux_boundary_d.inc.cuh"




void convectiveFlux_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_ns)
{
    // free-stream 保存: 基準静圧 pRef を device 定数へ転送 (既定 0.0 でビット不変)
    {
        flow_float pRef_h = static_cast<flow_float>(cfg.pRef);
        CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_pRef, &pRef_h, sizeof(flow_float)));
    }
    // U∞≠0 の移流基準差分 (KEEP CPG): 参照一様流を device 定数へ (既定 0.0 = off でビット不変)
    {
        flow_float roRef_h = static_cast<flow_float>(cfg.roRef);
        flow_float uRefX_h = static_cast<flow_float>(cfg.uRefX);
        flow_float uRefY_h = static_cast<flow_float>(cfg.uRefY);
        flow_float uRefZ_h = static_cast<flow_float>(cfg.uRefZ);
        CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_roRef, &roRef_h, sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_uRefX, &uRefX_h, sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_uRefY, &uRefY_h, sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemcpyToSymbol(d_uRefZ, &uRefZ_h, sizeof(flow_float)));
    }

    // D2a/D4 診断フラグを env から 1 度だけ device へ設定 (既定 off = ビット不変)。
    {
        static bool s_init = false;
        if (!s_init) {
            int c1 = 0, clog = 0; flow_float th = 0.05f, lth = 0.3f, blend = 0.0f;
            if (const char* e = getenv("FORGE_CONTACT_1ST"))       c1   = atoi(e);
            if (const char* e = getenv("FORGE_CONTACT_THRESH"))    th   = (flow_float)atof(e);
            if (const char* e = getenv("FORGE_CONTACT_LOG"))       clog = atoi(e);
            if (const char* e = getenv("FORGE_CONTACT_LOG_THRESH"))lth  = (flow_float)atof(e);
            if (const char* e = getenv("FORGE_CONTACT_BLEND"))     blend= (flow_float)atof(e);
            int fthy = 0;
            if (const char* e = getenv("FORGE_FACE_THERMOY"))      fthy = atoi(e);
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_faceThermoY,      &fthy, sizeof(int)));
            const int rem = (cfg.discretization == "node") ? 1 : 0;   // node は常にエッジ中点再構成
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_reconEdgeMid,     &rem,  sizeof(int)));
            const int sfr = cfg.speciesFaceReconstruction;   // config 由来 (env でない)
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_speciesFaceRecon, &sfr,  sizeof(int)));
            const int rycl = cfg.multispeciesRhoYCommonLimiter;   // config 由来 (opt-in 診断)
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_rhoYCommonLim,    &rycl, sizeof(int)));
            const unsigned long long zero = 0ULL;
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_speciesOvershoot, &zero, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contact1st,       &c1,   sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactThresh,    &th,   sizeof(flow_float)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactLog,       &clog, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactLogThresh, &lth,  sizeof(flow_float)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_contactBlend,     &blend,sizeof(flow_float)));
            s_init = true;
        }
    }

    // rho-Y 共通リミタ診断: 一定間隔で device カウンタを読み出して 1 行印字しリセット (opt-in 時のみ)。
    if (cfg.multispeciesRhoYCommonLimiter == 1) {
        static int s_rycl_call = 0;
        const int interval = 200;
        if ((s_rycl_call % interval) == 0) {
            int psimin=0, ymin=0, ymax=0;
            unsigned long long lt001=0, lt01=0, byrho=0, bysp=0, fb=0, ovs=0;
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&psimin, g_psiRhoY_min_scaled, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&lt001,  g_psiRhoY_lt001, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&lt01,   g_psiRhoY_lt01,  sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&byrho,  g_rhoYMinByRho,  sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&bysp,   g_rhoYMinBySpecies, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&fb,     g_rhoYFallback,  sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&ovs,    g_speciesOvershoot, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&ymin,   g_Yface_min_scaled, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&ymax,   g_Yface_max_scaled, sizeof(int)));
            printf("RHOYLIM call=%d minPsiRhoY=%.4f lt0.01=%llu lt0.1=%llu minBy[rho=%llu sp=%llu] overshoot=%llu fallback=%llu Yface=[%.5f,%.5f]\n",
                   s_rycl_call, (double)psimin*1e-6, lt001, lt01, byrho, bysp, ovs, fb,
                   (double)ymin*1e-6, (double)ymax*1e-6);
            // 次区間用にリセット (min/max は両端へ)。
            const int big=2000000, nbig=-2000000; const unsigned long long z=0ULL;
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_psiRhoY_min_scaled, &big, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_Yface_min_scaled,   &big, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_Yface_max_scaled,   &nbig, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_psiRhoY_lt001, &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_psiRhoY_lt01,  &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_rhoYMinByRho,  &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_rhoYMinBySpecies, &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_rhoYFallback,  &z, sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_speciesOvershoot, &z, sizeof(unsigned long long)));
        }
        s_rycl_call++;
    }

    // initialize
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_ro"]  , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roUx"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roUy"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roUz"], 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roe"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.p_d["massflux"], 0.0, msh.nPlanes*sizeof(flow_float)));

    dim3 dimGrid_normal_halo = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));

    // node-centered 弱形式 (Phase 2): 主対流ループは「内部双対面のみ」を処理し、末尾に並ぶ全非 periodic 境界 plane
    // (nBoundaryHaloPlanes = wall+inlet+outlet+slip...) を除外する。全境界は別途 convectiveFlux_boundary_d が bvar を
    // R 状態とする弱形式で担う。これは境界ノードが物理境界上に乗る node-centered で ghost を主ループの右状態に食わせる
    // と退化幾何 (d_along_n=0) により出口列/コーナーで近壁 BL が崩壊するため (case/26 で実証)。
    // **periodic も除外する** (median-dual M4, §4.5): node の周期境界は DOF 同一視 (periodicNodeGather + 合併体積)
    // で扱い、継ぎ目に双対面を作らない。周期半割面を主ループに入れると外向き境界パッチを CV 間面と誤用し発散する。
    // よって node の主ループ境界 = 純内部双対面 (nNormalPlanes)。cell モードは全 plane を主ループでゴースト処理 (従来どおり)。
    const geom_int convPlaneBound = (cfg.discretization == "node")
                                  ? msh.nNormalPlanes : msh.nNormal_halo_Planes;

    // (検証A) node 弱形式の plane 振り分けを 1 度だけログ: 主ループ面 = 内部(+periodic)、弱形式面 = 全境界半割面。
    if (cfg.discretization == "node") {
        static bool s_planeLogDone = false;
        if (!s_planeLogDone) {
            s_planeLogDone = true;
            const geom_int mainPlanes = convPlaneBound;
            const geom_int periodicInMain = mainPlanes - msh.nNormalPlanes;  // 主ループのうち内部以外 = periodic
            printf("[node weak-form] conv main planes = %ld (internal %ld + periodic %ld), "
                   "boundary weak planes = %ld (wall %ld + non-wall %ld)\n",
                   (long)mainPlanes, (long)msh.nNormalPlanes, (long)periodicInMain,
                   (long)msh.nBoundaryHaloPlanes, (long)msh.nWallHaloPlanes,
                   (long)(msh.nBoundaryHaloPlanes - msh.nWallHaloPlanes));
        }
    }

    // -----------------------
    // *** sum over planes ***
    // -----------------------
    // 非平衡凝縮: 二相エネルギー流束補正用の総液相分率 g。off / 未登録なら nullptr (補正なし)。
    // 現状 nCondSpecies=1 (g_0 が総 g)。多成分は総和 g 配列を別途用意する (TODO)。
    flow_float* cond_g = nullptr;
    if (cfg.condensation == 1 && var.nCondSpeciesRegistered >= 1) cond_g = var.c_d["g_0"];

    // 全スキーム共通の引数バンドルを 1 度だけ構築 (SLAU/HLLE/ROE で使い回す)。
    FaceGeom geom {
        msh.nCells, msh.nPlanes, msh.nNormalPlanes, msh.map_plane_cells_d,
        convPlaneBound, msh.normal_halo_planes_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],
        var.p_d["massflux"] };
    PrimState st {
        var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"],
        var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"], var.c_d["P"], var.c_d["Ht"], var.c_d["sonic"] };
    ResidualOut reso {
        var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"] };
    LimiterFields lim {
        var.c_d["limiter_ro"], var.c_d["limiter_Ux"], var.c_d["limiter_Uy"], var.c_d["limiter_Uz"], var.c_d["limiter_P"],
        var.c_d["ducros"] };
    GradFields grd {
        var.c_d["drodx"], var.c_d["drody"], var.c_d["drodz"],
        var.c_d["dUxdx"], var.c_d["dUxdy"], var.c_d["dUxdz"],
        var.c_d["dUydx"], var.c_d["dUydy"], var.c_d["dUydz"],
        var.c_d["dUzdx"], var.c_d["dUzdy"], var.c_d["dUzdz"],
        var.c_d["dPdx"] , var.c_d["dPdy"] , var.c_d["dPdz"] };
    // SST 全エネルギー E_t = E_m + ρk (sstEnergyIncludesK): 面エンタルピー +(5/3)k, 圧力流束 p* = p + (2/3)ρk。
    const bool sstEnergyK = (cfg.sstEnergyIncludesK != 0 && cfg.LESorRANS == 2 && cfg.RANSmodel == 1);
    if (sstEnergyK && !(cfg.solver == "SLAU" || cfg.solver == "SLAU2")) {
        throw std::runtime_error("turbulence.sstEnergyIncludesK=1 is implemented for solver SLAU/SLAU2 only");
    }
    CondArgs cnd { cfg.cp, cond_g, var.c_d["T"], cfg.condModel, sstEnergyK ? var.c_d["k"] : nullptr, sstEnergyK ? 1 : 0 };

    if (cfg.solver == "SLAU" || cfg.solver == "SLAU2") {
        int slauVariant = (cfg.solver == "SLAU2") ? 2 : 1;
        SpeciesArgs spA {
            cfg.thermalMethod, thermo_species_device_ptr(), cfg.nSpecies,
            species_roY_device_ptr(),
            species_Y_device_ptr(), species_dYdx_device_ptr(), species_dYdy_device_ptr(), species_dYdz_device_ptr(),
            species_limiterY_device_ptr(),
            (cfg.speciesFaceReconstruction >= 2) ? species_Yface_alloc(msh.nPlanes) : nullptr,
            var.c_d["Rmix"] };

        SLAU_d<<<dimGrid_normal_halo , cuda_cfg.dimBlock>>> (
            cfg.convMethod, cfg.limiter, slauVariant,
            cfg.lowMachPrecond, cfg.precondEps,
            cfg.lowMachThornber,
            cfg.gamma,
            spA, cnd, geom, st, reso, lim, grd
        ) ;

    } else if (cfg.solver == "HLLE") {
        HLLE_d<<<dimGrid_normal_halo , cuda_cfg.dimBlock>>> (
            cfg.convMethod, cfg.limiter,
            cfg.gamma,
            cnd, geom, st, reso, lim, grd
        );

    } else if (cfg.solver == "ROE") {
        //ROE_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>> ( 
        ROE_d<<<dimGrid_normal_halo , cuda_cfg.dimBlock>>> (
            cfg.convMethod, cfg.limiter,
            cfg.gamma,
            cfg.thermalMethod, thermo_species_device_ptr(), cfg.nSpecies,
            cnd, geom, st, reso, lim, grd
        );

    } else if (cfg.solver == "KEEP") {
        // 純粋 KEEP 中心流束 + opt-in ES 散逸レイヤ (keepDissType, 既定 0=散逸なし・ビット不変)。
        // LES/ILES 向け低散逸対流。SGS 散逸は WALE が担う。
        // TP (thermalMethod==2): 単成分 (Step 3)・多成分 (Step 4) とも scalar/matrix 対応。
        // 多成分 matrix の市松プラトーバグは根治済 (真因: ΣρY≠ρ 共通モードノイズの w 増幅。
        // カーネル内 Y 正規化で除去。plan 変更ログ 2026-07-19 深掘り参照)。
        KEEP_d<<<dimGrid_normal_halo , cuda_cfg.dimBlock>>> (
            cfg.gamma,
            cfg.thermalMethod, thermo_species_device_ptr(),
            cfg.nSpecies, species_roY_device_ptr(),
            cfg.keepDissType, cfg.keepDissCoeff,
            cfg.keepDissCprime, cfg.precondEps,
            cfg.keepDissJump, cfg.keepDissPrecond,
            // f_d 駆動 σ ブレンド (keepDissFdBlend=1 かつ DES のときのみ fd_shield を渡す。
            // config 検証で DESmode>0 は保証済み → fd_shield は毎サブ反復 turbulent_viscosity が更新)
            (cfg.keepDissFdBlend == 1) ? var.c_d["fd_shield"] : nullptr,
            (cfg.DESmode == 1) ? 1 : 0,   // fd_shield の向き (DDES: 1=LES / IDDES: 1=RANS)
            cfg.keepDissCoeffMax,
            cfg.keepDissCbCoeff, cfg.keepDissCbEps,
            cfg.keepDissOpBlendRaw,
            grd,
            geom, st, reso
        );

    } else {
        std::cerr << "Error: unsupported solver name " << cfg.solver
                  << " (enabled: SLAU, HLLE, ROE, KEEP)" << std::endl;
        exit(EXIT_FAILURE);
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();


    // cell モード: SLAU/ROE/HLLE は主ループが normal_halo_planes_d 経由で全境界 plane をゴースト処理するため、
    // この dedicated 境界カーネルは二重計上になるのでスキップする。
    // node モード (弱形式・Phase 2): 主ループは全非 periodic 境界 plane を除外した (convPlaneBound=内部+periodic) ので、
    // 全境界 (wall/inlet/outlet/slip...) の境界寄与を本 convectiveFlux_boundary_d が bvar を R 状態とする弱形式で担う。
    // (壁: bvar u=0 → mdot=0, pressure-only。出口: bvar Uxb/Psb。slip: Uxb=0+Psb=P[ic]。入口: 固定状態。)
    const bool nodeMode = (cfg.discretization == "node");
    bool skipBoundaryFluxKernel = (!nodeMode)
                               && (cfg.solver == "SLAU"
                                || cfg.solver == "SLAU2"
                                || cfg.solver == "ROE"
                                || cfg.solver == "HLLE"
                                || cfg.solver == "KEEP");

    for (auto& bc : msh.bconds)
    {
        if (bc.bcondKind == "periodic") {
            continue;
        }
        if (skipBoundaryFluxKernel) {
            continue;
        }
        convectiveFlux_boundary_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> (
            cfg.gamma,
            // 移流基準差分 (plan §8.5): 主ループ KEEP_d の advGauge (d_roRef>0 && CPG) と厳密同条件。
            // CV の全面に載らないと telescoping が破れるため、境界だけの on/off は不可。
            ((cfg.solver == "KEEP" && cfg.roRef > 0.0 && cfg.thermalMethod != 2) ? 1 : 0),
            // node × outlet_statPress: p_tilde を内部値 Ps[ic] で評価 (§2.11)。cell 経路は不変。
            ((nodeMode && bc.bcondKind == "outlet_statPress") ? 1 : 0),
            // mesh structure
            bc.iPlanes.size(),
            bc.map_bplane_plane_d,  
            bc.map_bplane_cell_d,  
            bc.map_bplane_cell_ghst_d,

            var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
            var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
            var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],
            (nodeMode ? var.p_d["massflux"] : nullptr),   // node 弱形式境界のみ massflux 書き戻し (cell は主ループが書く)

            // basic variables
            var.c_d["ro"] ,
            var.c_d["roUx"] ,
            var.c_d["roUy"] ,
            var.c_d["roUz"] ,
            var.c_d["roe"] ,
            var.c_d["Ux"]  ,
            var.c_d["Uy"]  ,
            var.c_d["Uz"]  ,
            var.c_d["P"]  ,
            var.c_d["Ht"]  ,
            var.c_d["sonic"]  ,
            var.c_d["T"]  ,

            bc.bvar_d["ro"],
            bc.bvar_d["roUx"],
            bc.bvar_d["roUy"],
            bc.bvar_d["roUz"],
            bc.bvar_d["roe"],
            bc.bvar_d["Ux"],
            bc.bvar_d["Uy"],
            bc.bvar_d["Uz"],
            bc.bvar_d["Tt"],
            bc.bvar_d["Pt"],
            bc.bvar_d["Ts"],
            bc.bvar_d["Ps"],
 
            var.c_d["res_ro"] ,
            var.c_d["res_roUx"] ,
            var.c_d["res_roUy"] ,
            var.c_d["res_roUz"] ,
            var.c_d["res_roe"]  ,
            // sstEnergyIncludesK: 内部側 k と境界側 k (bvar kb / 入口 k / 無ければ内部値)
            sstEnergyK ? var.c_d["k"] : nullptr,
            // 境界側 k: 入口 (Dirichlet, bvar "k" = 指定値) のみ bvar を使う。壁/出口/slip の bvar kb は node では未充填になり得るため内部値 (Neumann)
            (sstEnergyK && bc.bcondKind.rfind("inlet", 0) == 0 && bc.bvar_d.count("k")) ? bc.bvar_d["k"] : nullptr,
            sstEnergyK ? 1 : 0
        ) ;
    }


    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();


}