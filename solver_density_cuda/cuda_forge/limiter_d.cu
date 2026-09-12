#include "limiter_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"

__global__ void fill_limiter_d(flow_float* values, geom_int nValues, flow_float value)
{
    geom_int index = blockDim.x * blockIdx.x + threadIdx.x;

    if (index < nValues) {
        values[index] = value;
    }
}
__global__ void fill_limiter5_d(flow_float* a, flow_float* b, flow_float* c, flow_float* d, flow_float* e,
                                geom_int nValues, flow_float value)
{
    geom_int index = blockDim.x * blockIdx.x + threadIdx.x;
    if (index < nValues) { a[index] = value; b[index] = value; c[index] = value; d[index] = value; e[index] = value; }
}

// Limiters for Unstructured Higher-Order Accurate Solutions of the Euler Equations
// Krzysztof Michalak

__device__ flow_float venkata_limiter(flow_float delta_p_max, flow_float delta_p_min, 
                                      flow_float delta_m, flow_float volume) {

    flow_float K = 1.f;
    flow_float eps2 = K*K*K*volume;
    //return (x*x + 2.0*x + eps*eps)/(x*x + x + 2.0 + eps*eps);
    flow_float res;

    // K11: 元は /(...)/delta_m と除算2回。/(denom*delta_m) に統合して除算1回に（compute律速の limiter 向け）。
    if (delta_m > 1e-20f) {
        flow_float delta_p = delta_p_max;
        res = ((delta_p*delta_p+eps2)*delta_m +2*delta_m*delta_m*delta_p)
              /((delta_p*delta_p +2.0f*delta_m*delta_m +delta_p*delta_m +eps2)*delta_m);
    } else if (delta_m < -1e-20f) {
        flow_float delta_p = delta_p_min;
        res = ((delta_p*delta_p+eps2)*delta_m +2*delta_m*delta_m*delta_p)
              /((delta_p*delta_p +2.0f*delta_m*delta_m +delta_p*delta_m +eps2)*delta_m);
    } else {
        res = 1.0f;
    }

    return res;
}

__device__ flow_float barth_Jespersen_limiter(flow_float delta_p_max, flow_float delta_p_min, 
                                              flow_float delta_m, flow_float volume) {

    flow_float res;

    if (delta_m > 1e-20f) {
        res = min(1.0f, delta_p_max/delta_m);
    } else if (delta_m < -1e-20f) {
        res = min(1.0f, delta_p_min/delta_m);
    } else {
        res = 1.0f;
    }

    return min(res, 1.0f);
}


//__device__ flow_float nishikawa_r1_limiter(deltas delta_dash) {
//    flow_float res;
//    flow_float del = delta_dash.del;
//    flow_float rik = delta_dash.rik;
//
//    if (del < 0.0) {
//        res = 0.0;
//    } else if (del >= 0.0 && del < 1.0/rik) {
//        res = rik*del*(1.0+0.5*(1.0-rik)*del);
//    } else if (del < 1.0 && del >= 1.0/rik) {
//        res = 1.0 + rik*0.5/(1.0-rik)*(1.0-del)*(1.0-del);
//    } else if (del >= 1.0) {
//        res = 1.0;
//    }
//
//    return res;
//}

// Modified multi-dimensional limiting process with enhanced shock stability on unstructured grids
__global__ void limiter_r1_d
( 
 int limiter_scheme,
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells, 
 geom_int* cell_planes_index, geom_int* cell_planes,  

 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,

 // variables
 flow_float* Q  ,

 flow_float* limiter_Q  ,

 flow_float* dQdx  , flow_float* dQdy , flow_float* dQdz
) 
{
    geom_int ic0 = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic0 < nCells) {
        if (limiter_scheme == 0) {
            limiter_Q[ic0] = 1.0f;
            return;
        }

        geom_int ip; 
        geom_int ic1; 
        geom_float dcp_x; 
        geom_float dcp_y; 
        geom_float dcp_z; 

        geom_float dcp2_x; 
        geom_float dcp2_y; 
        geom_float dcp2_z; 

        geom_int index_st = cell_planes_index[ic0];
        geom_int index_en = cell_planes_index[ic0+1];
        geom_int np = index_en - index_st;

        flow_float Q_max =Q[ic0];
        flow_float Q_min =Q[ic0];

        flow_float denomi;

        flow_float limiter_Q_temp  = 1.0f;
        flow_float limiter_Q_temp2 = 1.0f;

        deltas delta;
        deltas delta_dash;
        flow_float volume = vol[ic0];

        int int_one=1;
        int int_two=2;
        int int_three=3;

        flow_float (*limiter_function)(flow_float, flow_float, flow_float, flow_float);

        if (limiter_scheme == int_one) { //barth
            limiter_function = barth_Jespersen_limiter;
        } else if (limiter_scheme == int_two or limiter_scheme == -1) { //venkata
            limiter_function = venkata_limiter;
        //} else if (limiter_scheme == int_three or limiter_scheme == -1) { //
        //    limiter_function = nishikawa_r1_limiter;
        } else {
            printf("Error: something wrong");
        }

        for (geom_int ilp=index_st; ilp<index_en; ilp++) {
            ip  = cell_planes[ilp];

            if (ip >= nNormalPlanes) continue;

            ic1 = plane_cells[2*ip+0] + plane_cells[2*ip+1] -ic0;

            Q_max = max(Q_max, Q[ic1]);
            Q_min = min(Q_min, Q[ic1]);
        }

        for (geom_int ilp=index_st; ilp<index_en; ilp++) {
            ip  = cell_planes[ilp];

            if (ip >= nNormalPlanes) continue;

            ic1 = plane_cells[2*ip+0] + plane_cells[2*ip+1] -ic0;

            dcp_x = pcx[ip] - ccx[ic0];
            dcp_y = pcy[ip] - ccy[ic0];
            dcp_z = pcz[ip] - ccz[ic0];

            // K9: dcp2_*, ri, rk, rik は未使用（dead）。compute律速の limiter から sqrt 2回/面/変数を除去。

            flow_float delta_p_max;
            flow_float delta_p_min;
            flow_float delta_m;

            flow_float Qt = Q[ic0] + dQdx[ic0]*dcp_x + dQdy[ic0]*dcp_y + dQdz[ic0]*dcp_z;

            delta_p_max = Q_max - Q[ic0];
            delta_p_min = Q_min - Q[ic0];
            delta_m     = Qt    - Q[ic0];

            limiter_Q_temp = limiter_function(delta_p_max, delta_p_min, delta_m, volume);

            limiter_Q_temp2 = min(limiter_Q_temp2, limiter_Q_temp);
        }

        limiter_Q[ic0] = min(max(limiter_Q_temp2, 0.0f),1.0f);

            //ic1 = plane_cells[2*ip+0] + plane_cells[2*ip+1] -ic0;
        //limiter_Q[ic1] = dQdx[ic0];
    }
}


// limiter_r1_d の 5 変数 (ro,Ux,Uy,Uz,P) 融合版。connectivity/geometry を 1 回読み・plane ループを共有して
// 5 変数の min/max と limiter を同時計算する (per-variable 5 回 launch の冗長 geometry 読みを除去)。
// 数式は limiter_r1_d と同一 (#pragma unroll で k は定数化 → Qk 等は constant param のまま, ポインタ配列は消える)。
// limiter_scheme は template 引数 (SCHEME) にして、リミッタ関数をコンパイル時に確定させる。
// 旧実装の関数ポインタ経由の呼び出しは device 上で間接 CALL になり (インライン化不可・レジスタ退避)、
// 2 パス × 5 変数 × 近傍面のループで 1 セルあたり数十回の CALL が入って 3D 2.37 M 節点で 9.5 ms/step を
// 食っていた (plans/active/performance-3d-node-sst-speedup.md §4.1)。数式・演算順序は同一 (ビット同一)。
template<int SCHEME>
__global__ void limiter_r1_fused5_d
(
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,
 geom_int* cell_planes_index, geom_int* cell_planes,
 geom_float* vol, geom_float* ccx, geom_float* ccy, geom_float* ccz,
 geom_float* pcx, geom_float* pcy, geom_float* pcz, geom_float* fx,
 flow_float* Q0, flow_float* Q1, flow_float* Q2, flow_float* Q3, flow_float* Q4,
 flow_float* L0, flow_float* L1, flow_float* L2, flow_float* L3, flow_float* L4,
 flow_float* d0x, flow_float* d0y, flow_float* d0z,
 flow_float* d1x, flow_float* d1y, flow_float* d1z,
 flow_float* d2x, flow_float* d2y, flow_float* d2z,
 flow_float* d3x, flow_float* d3y, flow_float* d3z,
 flow_float* d4x, flow_float* d4y, flow_float* d4z
)
{
    geom_int ic0 = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic0 >= nCells) return;

    flow_float* Q[5]   = {Q0,Q1,Q2,Q3,Q4};
    flow_float* Lim[5] = {L0,L1,L2,L3,L4};
    flow_float* dQx[5] = {d0x,d1x,d2x,d3x,d4x};
    flow_float* dQy[5] = {d0y,d1y,d2y,d3y,d4y};
    flow_float* dQz[5] = {d0z,d1z,d2z,d3z,d4z};

    if (SCHEME == 0) {
        #pragma unroll
        for (int k=0;k<5;k++) Lim[k][ic0] = 1.0f;
        return;
    }

    const geom_int index_st = cell_planes_index[ic0];
    const geom_int index_en = cell_planes_index[ic0+1];
    const flow_float volume  = vol[ic0];
    const geom_float cx0 = ccx[ic0], cy0 = ccy[ic0], cz0 = ccz[ic0];

    flow_float qc[5], qmax[5], qmin[5], gx[5], gy[5], gz[5], ltmp[5];
    #pragma unroll
    for (int k=0;k<5;k++){
        qc[k]=Q[k][ic0]; qmax[k]=qc[k]; qmin[k]=qc[k];
        gx[k]=dQx[k][ic0]; gy[k]=dQy[k][ic0]; gz[k]=dQz[k][ic0]; ltmp[k]=1.0f;
    }

    // pass1: neighbor min/max (geometry/connectivity を 1 回だけ読む)
    for (geom_int ilp=index_st; ilp<index_en; ilp++) {
        geom_int ip = cell_planes[ilp];
        if (ip >= nNormalPlanes) continue;
        geom_int ic1 = plane_cells[2*ip+0] + plane_cells[2*ip+1] - ic0;
        #pragma unroll
        for (int k=0;k<5;k++){ flow_float qn=Q[k][ic1]; qmax[k]=max(qmax[k],qn); qmin[k]=min(qmin[k],qn); }
    }

    // pass2: limiter
    for (geom_int ilp=index_st; ilp<index_en; ilp++) {
        geom_int ip = cell_planes[ilp];
        if (ip >= nNormalPlanes) continue;
        flow_float dcp_x = pcx[ip]-cx0;
        flow_float dcp_y = pcy[ip]-cy0;
        flow_float dcp_z = pcz[ip]-cz0;
        #pragma unroll
        for (int k=0;k<5;k++){
            flow_float Qt = qc[k] + gx[k]*dcp_x + gy[k]*dcp_y + gz[k]*dcp_z;
            const flow_float lk = (SCHEME == 1)
                ? barth_Jespersen_limiter(qmax[k]-qc[k], qmin[k]-qc[k], Qt-qc[k], volume)
                : venkata_limiter        (qmax[k]-qc[k], qmin[k]-qc[k], Qt-qc[k], volume);
            ltmp[k] = min(ltmp[k], lk);
        }
    }

    #pragma unroll
    for (int k=0;k<5;k++) Lim[k][ic0] = min(max(ltmp[k], (flow_float)0.0), (flow_float)1.0);
}


void limiter_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    // 5 配列の 1.0 充填を 1 カーネルに (起動 5→1)。
    fill_limiter5_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        var.c_d["limiter_ro"], var.c_d["limiter_Ux"], var.c_d["limiter_Uy"], var.c_d["limiter_Uz"], var.c_d["limiter_P"],
        msh.nCells_all, 1.0f);

    // limiter<=0: 0=明示 off、-1=「リミタ off」(solverConfig.cpp が受理する正式値。KEEP は lim を
    // 一切参照しないため計算自体が無駄 — 修正前は == 0 のみ早期 return しており、-1 が Venkatakrishnan
    // 分岐 (venkata_limiter, limiter_scheme==2 と同一コスト) に落ちて全 cell の全変数を計算していた
    // = KEEP 系 run の実測 ~22% の GPU 時間が丸ごと無駄だった)。fill 済みの 1.0 がそのまま「無制限」を表す。
    if (cfg.limiter <= 0) {
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchkKernelSync();
        return;
    }

    // ro,Ux,Uy,Uz,P を 1 カーネルに融合 (connectivity/geometry の 5 重読みを除去)。数式は per-variable と同一。
    // SCHEME: 1=Barth-Jespersen, それ以外 (2 / -1 は上で return 済) = Venkatakrishnan。
    #define FORGE_LIMITER_FUSED5_ARGS \
        msh.nCells, \
        msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d, \
        msh.map_cell_planes_index_d , msh.map_cell_planes_d , \
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"], \
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"], \
        var.c_d["ro"], var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"], var.c_d["P"], \
        var.c_d["limiter_ro"], var.c_d["limiter_Ux"], var.c_d["limiter_Uy"], var.c_d["limiter_Uz"], var.c_d["limiter_P"], \
        var.c_d["drodx"], var.c_d["drody"], var.c_d["drodz"], \
        var.c_d["dUxdx"], var.c_d["dUxdy"], var.c_d["dUxdz"], \
        var.c_d["dUydx"], var.c_d["dUydy"], var.c_d["dUydz"], \
        var.c_d["dUzdx"], var.c_d["dUzdy"], var.c_d["dUzdz"], \
        var.c_d["dPdx"] , var.c_d["dPdy"] , var.c_d["dPdz"]
    if (cfg.limiter == 1)
        limiter_r1_fused5_d<1><<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>> (FORGE_LIMITER_FUSED5_ARGS);
    else
        limiter_r1_fused5_d<2><<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>> (FORGE_LIMITER_FUSED5_ARGS);
    #undef FORGE_LIMITER_FUSED5_ARGS

    // 多成分 face 整合再構成: 各化学種 Y_s に Venkat リミタ ψ_Y を計算 (∇Y は speciesGradient 済)。
    // speciesFaceReconstruction==1 のみ。flux では min(ψ_ρ, ψ_Y) を Y 再構成に使う (boundedness)。
    if (cfg.speciesFaceReconstruction >= 1 && var.nSpeciesRegistered >= 2) {
        for (int s = 0; s < var.nSpeciesRegistered; ++s) {
            const std::string i = std::to_string(s);
            fill_limiter_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(var.c_d["limiter_Y"+i], msh.nCells_all, 1.0f);
            limiter_r1_d<<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>> (
                cfg.limiter, msh.nCells, msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
                msh.map_cell_planes_index_d , msh.map_cell_planes_d ,
                var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
                var.c_d["Y"+i] , var.c_d["limiter_Y"+i] ,
                var.c_d["dY"+i+"dx"] , var.c_d["dY"+i+"dy"] , var.c_d["dY"+i+"dz"]
            ) ;
        }
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

}