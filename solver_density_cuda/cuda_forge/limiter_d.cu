#include "limiter_d.cuh"
#include "calcGradient_d.cuh"
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

#include "limiterFunctions_d.cuh"
#include "passiveLimiter_d.cuh"
#include "limiterPeriodic_d.cuh"
#include "passiveTransport_d.cuh"
#include "periodicNode_d.cuh"

// 周期 node (合併 CV) の 2 段リミッタ用スクラッチ (極値 Q_max / Q_min)。nCells_all で 1 回だけ確保し使い回す
// (plans/active/species-passive-scalar-unification.md §4.8)。
static flow_float* s_lim_qmax = nullptr;
static flow_float* s_lim_qmin = nullptr;
static geom_int    s_lim_scratch_n = 0;
static void periodicLimiterScratch(const mesh& msh)
{
    if (s_lim_qmax != nullptr && s_lim_scratch_n == msh.nCells_all) return;
    if (s_lim_qmax != nullptr) { cudaFree(s_lim_qmax); cudaFree(s_lim_qmin); }
    gpuErrchk( cudaMalloc((void**)&s_lim_qmax, sizeof(flow_float)*msh.nCells_all) );
    gpuErrchk( cudaMalloc((void**)&s_lim_qmin, sizeof(flow_float)*msh.nCells_all) );
    s_lim_scratch_n = msh.nCells_all;
}

// 周期 node の 1 変数 2 段リミッタ: 極値 → group max/min gather → 合併極値で自分の面の ψ → group min gather。
// SCALED=true は受動種の無次元化版 (limiter_r1_scaled_d 相当)、false は limiter_r1_d 相当。
// 非周期・cell では呼ばない (従来の 1 段 kernel がビット不変で残る)。
template<bool SCALED>
static void limiter_periodic_merged
(
 solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
 flow_float phi_floor, flow_float* Q, flow_float* limiter_Q,
 flow_float* dQdx, flow_float* dQdy, flow_float* dQdz
)
{
    periodicLimiterScratch(msh);
    limiter_extrema_d<<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>>(
        msh.nCells, msh.nNormalPlanes, msh.map_plane_cells_d,
        msh.map_cell_planes_index_d, msh.map_cell_planes_d,
        Q, s_lim_qmax, s_lim_qmin);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    periodicGatherMaxArray_d_wrapper(cfg, cuda_cfg, msh, s_lim_qmax);
    periodicGatherMinArray_d_wrapper(cfg, cuda_cfg, msh, s_lim_qmin);
    limiter_psi_merged_d<SCALED><<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>>(
        cfg.limiter, msh.nCells, msh.nNormalPlanes,
        msh.map_cell_planes_index_d, msh.map_cell_planes_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"], var.p_d["pcy"], var.p_d["pcz"],
        phi_floor, Q, s_lim_qmax, s_lim_qmin, limiter_Q, dQdx, dQdy, dQdz);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    periodicGatherMinArray_d_wrapper(cfg, cuda_cfg, msh, limiter_Q);
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
 flow_float* d4x, flow_float* d4y, flow_float* d4z,
 const flow_float* __restrict__ prim   // 原始量 AoS パック [ro,Ux,Uy,Uz,P,...] (nullptr なら Q0..Q4 から gather)
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
        if (prim != nullptr) {
            const float4 a = *reinterpret_cast<const float4*>(prim + (size_t)ic1*8);
            const flow_float qn5 = prim[(size_t)ic1*8 + 4];
            const flow_float qn[5] = {a.x, a.y, a.z, a.w, qn5};
            #pragma unroll
            for (int k=0;k<5;k++){ qmax[k]=max(qmax[k],qn[k]); qmin[k]=min(qmin[k],qn[k]); }
        } else {
        #pragma unroll
        for (int k=0;k<5;k++){ flow_float qn=Q[k][ic1]; qmax[k]=max(qmax[k],qn); qmin[k]=min(qmin[k],qn); }
        }
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
        passiveLimiter_d_wrapper(cfg, cuda_cfg, msh, var);   // 受動種 ψ_P (limiter<=0 では 1.0 充填)
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
        var.c_d["dPdx"] , var.c_d["dPdy"] , var.c_d["dPdz"], \
        (const flow_float*)nullptr   /* 原始量パックは不採用 (実測で利得なし) */
    // 周期 node (合併 CV) は 2 段 (極値の group max/min → ψ の group min) で周期対の ψ を一致させる (§4.8)。
    const bool perNode = periodicNodeActive(cfg, msh);
    if (perNode) {
        const char* qn[5]  = {"ro","Ux","Uy","Uz","P"};
        const char* ln[5]  = {"limiter_ro","limiter_Ux","limiter_Uy","limiter_Uz","limiter_P"};
        const char* gxn[5] = {"drodx","dUxdx","dUydx","dUzdx","dPdx"};
        const char* gyn[5] = {"drody","dUxdy","dUydy","dUzdy","dPdy"};
        const char* gzn[5] = {"drodz","dUxdz","dUydz","dUzdz","dPdz"};
        for (int k = 0; k < 5; ++k)
            limiter_periodic_merged<false>(cfg, cuda_cfg, msh, var, 0.0f,
                var.c_d[qn[k]], var.c_d[ln[k]], var.c_d[gxn[k]], var.c_d[gyn[k]], var.c_d[gzn[k]]);
    } else if (cfg.limiter == 1)
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
            if (perNode) {
                limiter_periodic_merged<false>(cfg, cuda_cfg, msh, var, 0.0f,
                    var.c_d["Y"+i], var.c_d["limiter_Y"+i], var.c_d["dY"+i+"dx"], var.c_d["dY"+i+"dy"], var.c_d["dY"+i+"dz"]);
                continue;
            }
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

    // 受動種 (passiveScalarScheme 1, speciesFaceReconstruction>=1): 受動種ごとの無次元化 Venkat ψ_P。
    passiveLimiter_d_wrapper(cfg, cuda_cfg, msh, var);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

// 受動種の無次元化 Venkat リミッタ ψ_P (passiveLimiter_d.cuh)。limiter<=0 または speciesFaceReconstruction<1 では 1.0 充填のみ。
void passiveLimiter_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (!passiveSchemeEnabled(cfg)) return;
    const int n = passive_count();
    const auto& prims = passive_prim_names();
    for (int q = 0; q < n; ++q) {
        fill_limiter_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(var.c_d["limiter_"+prims[q]], msh.nCells_all, 1.0f);
    }
    if (cfg.limiter <= 0 || cfg.speciesFaceReconstruction < 1) return;
    const bool perNode = periodicNodeActive(cfg, msh);
    for (int q = 0; q < n; ++q) {
        const std::string& pn = prims[q];
        if (perNode) {   // 周期 node: 合併極値 (φ_ref も合併極値) の 2 段 ψ_P (§4.8)
            limiter_periodic_merged<true>(cfg, cuda_cfg, msh, var, static_cast<flow_float>(1.0e-30),
                var.c_d[pn], var.c_d["limiter_"+pn], var.c_d["d"+pn+"dx"], var.c_d["d"+pn+"dy"], var.c_d["d"+pn+"dz"]);
            continue;
        }
        limiter_r1_scaled_d<<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>> (
            cfg.limiter, msh.nCells, msh.nNormalPlanes, msh.map_plane_cells_d,
            msh.map_cell_planes_index_d, msh.map_cell_planes_d,
            var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
            var.p_d["pcx"], var.p_d["pcy"], var.p_d["pcz"],
            static_cast<flow_float>(1.0e-30),
            var.c_d[pn], var.c_d["limiter_"+pn],
            var.c_d["d"+pn+"dx"], var.c_d["d"+pn+"dy"], var.c_d["d"+pn+"dz"]);
    }
}
