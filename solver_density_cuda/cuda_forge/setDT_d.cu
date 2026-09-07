//#if defined(__CUDA_ARCH__)
//#pragma message "content of __CUDA_ARCH__: " __CUDA_ARCH__
//#endif

#include "setDT_d.cuh"
#include "lowMachPrecond_d.cuh"   // Phase 4: β・c' で Δτ' を前処理スペクトル半径に合わせて拡大


__global__ void setCFL_pln_d
( 
 flow_float dt,
 //flow_float dt_local,
 flow_float visc,
 flow_float* vis_turb  ,

 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* cfl ,
 //flow_float* cfl_pseudo ,
 flow_float* sonic,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,

 //plane variables
 flow_float* cfl_pln
 //flow_float* cfl_pseudo_pln

)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;


    if (ip < nPlanes) {

        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        //__syncthreads();

        geom_float f = fx[ip];

        geom_float vol0 = vol[ic0];
        geom_float vol1 = vol[ic1];
        
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];


        geom_float dx0 = vol0/sss;
        geom_float dx1 = vol1/sss;
        geom_float dx_min = min(dx0,dx1);

        flow_float Ux0 = Ux[ic0];
        flow_float Uy0 = Uy[ic0];
        flow_float Uz0 = Uz[ic0];

        flow_float Ux1 = Ux[ic1];
        flow_float Uy1 = Uy[ic1];
        flow_float Uz1 = Uz[ic1];

        flow_float US  = (f*Ux0 + (1.0-f)*Ux1)*sxx
                        +(f*Uy0 + (1.0-f)*Uy1)*syy
                        +(f*Uz0 + (1.0-f)*Uz1)*szz;

        flow_float rof = f*ro[ic0] + (1.0-f)*ro[ic1];
        flow_float v_turb = f*vis_turb[ic0] + (1.0-f)*vis_turb[ic1];
        flow_float lambda = abs(US)/sss + sonic[ic0] + 2.0*(visc+v_turb)/(rof*dx_min);

        cfl_pln[ip] = dt*lambda/dx_min;

        //cfl_pseudo_pln[ip] = dt_local*lambda/dx_min;
    }
}

__global__ void setCFL_cell_d
( 
 int dtControl, flow_float cfl_target, flow_float cfl_pseudo_target,
 int dualTime, int unsteady, 
 flow_float dt,

 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* cell_planes_index, geom_int* cell_planes,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* cfl ,
 flow_float* dt_local  ,
 flow_float* sonic,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,

 //plane variables
 flow_float* cfl_pln,
 //flow_float* cfl_pseudo_pln

 // 軸対称 near-axis 安定化 (isAxisymmetric==1 かつ axisBeta>0 のときのみ作用)
 int isAxisymmetric,
 flow_float* A_planar,
 flow_float axisBeta,

 // v2 line-implicit 粘性 CFL 割引 (plans/active/time_integration-line-implicit-viscous-v2.md):
 // on-line セルに限り、面 λ の粘性項 2ν_eff/(ρ·dx_min) 由来の CFL を (1−θ) 倍に割引く。
 // 壁近傍セルの dx_min は壁法線 (=line 方向) なので方向整合は近似的に成立。θ=0 で完全不変。
 geom_int* plane_cells,
 flow_float visc_lam,
 flow_float* vis_turb,
 flow_float lineReliefTheta,
 const geom_int* line_prev,
 const geom_int* line_next,
 // 方向別 dt (lineDtDirectional==1): line 面 (Thomas が厳密に解く結合) の λ を CFL の max から
 // 完全除外し、Δτ を off-line 面 (lag 側) の λ だけで決める。壁セルの Δτ は streamwise 基準 (×AR)。
 int lineDtDirectional,
 // 診断 (lineDtWallRelief==1): wall 種境界半割面の λ も on-line セルの max から除外 (壁端点律速の切り分け)。
 const unsigned char* plane_wall,
 int lineDtWallRelief
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {

        geom_int index_st = cell_planes_index[ic];
        geom_int index_en = cell_planes_index[ic+1];
        geom_int np = index_en - index_st;

        cfl[ic] = 0.0;
        //cfl_pseudo[ic] = 0.0;

        const bool onLine = (line_prev != nullptr && (line_prev[ic] >= 0 || line_next[ic] >= 0));
        const bool lineRelief = (lineReliefTheta > (flow_float)0.0 && onLine);

        for (geom_int ilp=index_st; ilp<index_en; ilp++) {
            geom_int ip = cell_planes[ilp];

            if (lineDtWallRelief != 0 && onLine && plane_wall != nullptr && plane_wall[ip] != 0) continue;
            if (lineDtDirectional != 0 && onLine) {
                const geom_int ic0 = plane_cells[2*ip+0];
                const geom_int ic1 = plane_cells[2*ip+1];
                const geom_int other = (ic0 == ic) ? ic1 : ic0;
                if (other == line_prev[ic] || other == line_next[ic]) continue;  // line 面は除外
            }
            flow_float cfl_p = cfl_pln[ip];
            if (lineRelief) {
                // setCFL_pln_d と同一式で面の粘性 CFL を再計算し θ 分を引く (対流+音響分は必ず残る)。
                const geom_int ic0 = plane_cells[2*ip+0];
                const geom_int ic1 = plane_cells[2*ip+1];
                const flow_float f = fx[ip];
                const flow_float sss = ss[ip];
                const flow_float dx_min = min(vol[ic0]/sss, vol[ic1]/sss);
                const flow_float rof = f*ro[ic0] + ((flow_float)1.0-f)*ro[ic1];
                const flow_float v_turb = f*vis_turb[ic0] + ((flow_float)1.0-f)*vis_turb[ic1];
                const flow_float cfl_visc = dt * (flow_float)2.0 * (visc_lam + v_turb)
                                            / (rof * dx_min * dx_min);
                cfl_p = max(cfl_p - lineReliefTheta * cfl_visc, (flow_float)0.0);
            }
            cfl[ic] = max(cfl[ic], cfl_p);
        }

        // 軸対称 near-axis 半径音響スペクトル半径を加える: λ_axis = β·(|u_r|+c)·A_planar。
        // revolved 軸面積 (r_f·S→0) が落とす半径音響モードを planar 面積で補う。face 項と同じ
        // 無次元化 (dt·λ/V) で cfl へ加算 → 近軸で Δτ=cfl_pseudo·dt/cfl ∝ cfl_pseudo·r/(|u_r|+c)。
        // V=r·A_planar (per-radian) なので A_planar/V=1/r。LHS 時間項 v/Δτ が全保存式で自動的に強まる。
        if (isAxisymmetric == 1 && axisBeta > (flow_float)0.0) {
            const flow_float u_r      = fabsf(Uy[ic]);                                 // 半径方向速度 (axisym: y=radial)
            const flow_float lam_axis = axisBeta * (u_r + sonic[ic]) * A_planar[ic];   // [V/time]
            cfl[ic] += dt * lam_axis / max(vol[ic], (flow_float)1.0e-30);
        }

        //if (dualTime == 1 or unsteady == 0) {
        //    dt_local[ic] = cfl_pseudo_target*dt/cfl[ic];
        //} else {
        //    dt_local[ic] = dt;
        //}
    }
}

__global__ void setDTlocal_pseudo_cell_d
( 
 flow_float cfl_pseudo_target,
 flow_float dt,

 // mesh structure
 geom_int nCells,
 flow_float* cfl ,
 flow_float* dt_local
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        const flow_float local_cfl_rate = cfl[ic] / max(dt, static_cast<flow_float>(1.0e-30));
        dt_local[ic] = cfl_pseudo_target / max(local_cfl_rate, static_cast<flow_float>(1.0e-30));
    }
}

__global__ void setDTlocal_uniform_cell_d
( 
 flow_float dt,
 geom_int nCells,

 flow_float* dt_local

)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        dt_local[ic] = dt;
    }
}

// Phase 4 (lowMachPrecond>=2): 擬似時間刻みを前処理スペクトル半径に合わせて拡大する。
// 既存 dt_local は物理スペクトル半径 λ_phys=|u|+c 基準。前処理後の律速は
// λ'=½(1+β)|u|+c' (低マッハで小) なので、dt_local ×= λ_phys/λ' で擬似 CFL を λ' 基準に合わせる。
// ε フロアにより λ'≥ε·c 程度に下限が付き、倍率は ~1/ε で有界 (発散しない)。
// LHS を前処理した本モードでのみ整合的 (Phase 2 で LHS 非前処理のまま setDT 前処理して破綻した轍を踏まない)。
__global__ void setDTlocal_precond_scale_d
(
 flow_float precondEps,
 geom_int nCells,
 flow_float* ro,
 flow_float* Ux, flow_float* Uy, flow_float* Uz,
 flow_float* sonic,
 flow_float* dt_local
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float c = max(sonic[ic], static_cast<flow_float>(1.0e-8));
        const flow_float vmag = sqrt(Ux[ic]*Ux[ic] + Uy[ic]*Uy[ic] + Uz[ic]*Uz[ic]);
        const flow_float beta = lowMachBeta(c, vmag, precondEps);
        const flow_float cprime = lowMachCprime(c, vmag, vmag, precondEps);
        const flow_float lam_phys = vmag + c;
        const flow_float lam_prec = static_cast<flow_float>(0.5)*(static_cast<flow_float>(1.0)+beta)*vmag + cprime;
        const flow_float ratio = lam_phys / max(lam_prec, static_cast<flow_float>(1.0e-30));
        dt_local[ic] *= max(ratio, static_cast<flow_float>(1.0));
    }
}


// 診断: 出口近傍の局所 dt キャップ (壁∩出口コーナー不安定の切り分け用, env ゲート)。
//   FORGE_DT_OUTLET_SCALE (0<s<1) と FORGE_DT_OUTLET_XMIN [m] の両方を与えると、
//   x > XMIN のセルの dt_local を s 倍に縮める。既定 (未設定) は完全に不変。
__global__ void scaleOutletDt_d(geom_int nCells, const flow_float* ccx,
                                flow_float xmin, flow_float scale, flow_float* dt_local)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells && ccx[ic] > xmin) dt_local[ic] *= scale;
}

void setDT_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , bool adaptDt , bool printCfl)
{
    setCFL_pln_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>> ( 
        cfg.dt,
        //cfg.dt_local,
        cfg.visc,
        var.c_d["vis_turb"] ,

        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

        // basic variables
        //var.c_d["convx"] , var.c_d["convy"] , var.c_d["convz"] ,
        //var.c_d["diffx"] , var.c_d["diffy"] , var.c_d["diffz"] ,
        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["cfl"]  , 
        //var.c_d["cfl_pseudo"] , 
        var.c_d["sonic"]  , 
        var.c_d["Ux"]  , 
        var.c_d["Uy"]  , 
        var.c_d["Uz"]  ,

        var.p_d["cfl_pln"]
        //var.p_d["cfl_pseudo_pln"]
    ) ;
    gpuErrchk( cudaPeekAtLastError() );
    // 注: 同一 default stream なので後続カーネルは順序保証される。per-step 同期削減のため
    // 中間の cudaDeviceSynchronize は撤去し peek のみとする (誤差チェックは末尾/次同期点で担保)。

    setCFL_cell_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
        cfg.dtControl, cfg.cfl, cfg.cfl_pseudo,
        cfg.dualTime, cfg.unsteady,
        cfg.dt,
        //cfg.dt_local,

        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes, msh.map_cell_planes_index_d , msh.map_cell_planes_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

        // basic variables
        //var.c_d["convx"] , var.c_d["convy"] , var.c_d["convz"] ,
        //var.c_d["diffx"] , var.c_d["diffy"] , var.c_d["diffz"] ,
        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["cfl"]  , 
        var.c_d["dt_local"] ,
        var.c_d["sonic"]  ,
        var.c_d["Ux"]  ,
        var.c_d["Uy"]  ,
        var.c_d["Uz"]  ,

        var.p_d["cfl_pln"] ,
        //var.p_d["cfl_pseudo_pln"]

        // 軸対称 near-axis 安定化
        cfg.isAxisymmetric,
        (cfg.isAxisymmetric == 1) ? var.c_d["A_planar"] : var.c_d["volume"],
        cfg.axisTimestepBeta,

        // v2 line-implicit 粘性 CFL 割引 (lineViscousDtRelief=0 で完全不変)
        msh.map_plane_cells_d,
        cfg.visc,
        var.c_d["vis_turb"],
        (cfg.lineImplicit == 1) ? cfg.lineViscousDtRelief : (flow_float)0.0,
        (cfg.lineImplicit == 1) ? msh.line_prev_d : nullptr,
        (cfg.lineImplicit == 1) ? msh.line_next_d : nullptr,
        (cfg.lineImplicit == 1) ? cfg.lineDtDirectional : 0,
        (cfg.lineImplicit == 1) ? msh.plane_wall_flag_d : nullptr,
        (cfg.lineImplicit == 1) ? cfg.lineDtWallRelief : 0
    ) ;
    gpuErrchk( cudaPeekAtLastError() );

    if (cfg.dualTime == 1 or cfg.unsteady == 0) {
        setDTlocal_pseudo_cell_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
            cfg.cfl_pseudo,
            cfg.dt,
            msh.nCells,
            var.c_d["cfl"],
            var.c_d["dt_local"]
        );

    } else {
        setDTlocal_uniform_cell_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
            cfg.dt,
            msh.nCells,
            var.c_d["dt_local"]
        );
    }

    // Phase 4: 完全前処理モード (LHS) では擬似時間刻みを前処理スペクトル半径基準に拡大する。
    // lowMachPrecond==2 (RHS+LHS) / ==3 (LHS のみ) の双方で擬似時間項は前処理されるため両方で適用する。
    if (cfg.lowMachPrecond >= 2) {
        setDTlocal_precond_scale_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
            cfg.precondEps,
            msh.nCells,
            var.c_d["ro"],
            var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"],
            var.c_d["sonic"],
            var.c_d["dt_local"]
        );
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchkKernelSync();
    }

    // max cfl の host 読み出し (thrust::max_element は D2H + 同期)。
    //  - printCfl: printf が値を host で要するため、本質的に host read が必要。
    //  - dt 適応 (adaptDt && dtControl==1): 原理的には device 上で完結できる (cfl_max も cfg.dt も device 化すれば
    //    host 同期不要) が、現状 cfg.dt が host スカラ (多数カーネルに値渡し・dual-time BDF 係数等で host 使用) の
    //    ため host で計算している。よってこの条件は「現実装の都合」であり本質的要請ではない。
    //    device-resident dt 化は explicit/dual-time の毎ステップ適応 (dtControl==1) で per-step 同期を消せるが、
    //    定常 implicit は dt_local=cfl_pseudo·dx/λ で cfg.dt が打ち消され不影響なので利得なし (別 plan 候補)。
    // 両条件とも false なら host 同期は一切発生しない (後続カーネルは同一 default stream で順序保証)。
    // 出口近傍 dt キャップ (診断, env ゲート): dt_local 確定後に縮める。既定 (未設定) は不変。
    {
        static const double dtOutletScale = [](){ const char* e = getenv("FORGE_DT_OUTLET_SCALE"); return e ? atof(e) : 0.0; }();
        static const double dtOutletXmin  = [](){ const char* e = getenv("FORGE_DT_OUTLET_XMIN");  return e ? atof(e) : 0.0; }();
        if (dtOutletScale > 0.0 && dtOutletScale < 1.0) {
            scaleOutletDt_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
                msh.nCells, var.c_d["ccx"], (flow_float)dtOutletXmin, (flow_float)dtOutletScale, var.c_d["dt_local"]);
            gpuErrchk( cudaPeekAtLastError() );
        }
    }

    const bool needHostRead = (adaptDt && cfg.dtControl == 1) || printCfl;
    if (needHostRead) {
        thrust::device_ptr<flow_float> d_ptr = thrust::device_pointer_cast(var.c_d["cfl"]);
        const flow_float cfl_max = *(thrust::max_element(d_ptr, d_ptr + msh.nCells));

        if (adaptDt && cfg.dtControl == 1) { // cfl based time control
            flow_float cfl_target = cfg.cfl;
            cfg.dt = cfg.dt*cfl_target/cfl_max;

            cfg.dt = max(cfg.dt, cfg.dt_min);
            cfg.dt = min(cfg.dt, cfg.dt_max);
        }

        if (printCfl) {
            printf("  max cfl : %f    \n", cfl_max);
            printf("  dt      : %e [s]\n", cfg.dt);
        }
    }

    gpuErrchk( cudaPeekAtLastError() );
}