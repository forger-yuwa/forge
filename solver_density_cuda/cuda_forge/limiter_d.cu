#include "limiter_d.cuh"
#include "cuda_forge/reconIncrement_d.cuh"   // 再構成増分の唯一の定義 (流束と共有)
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
// 有界性診断 (§4.20 / W1h) の device カウンタ。周期経路の検査カーネルからも使うので前方に置く。
__device__ int g_limDiag = 0;
__device__ unsigned long long g_limG1[5] = {0,0,0,0,0};   // ro, Ux, Uy, Uz, P
__device__ unsigned long long g_limSides = 0;
__device__ unsigned long long g_limNonFinite[5] = {0,0,0,0,0};   // 非有限を「逸脱なし」にしない (codex C2)
__device__ unsigned long long g_limMaxExcess[5] = {0,0,0,0,0};   // 許容幅に対する最大逸脱倍率 x1e3

template<bool SCALED>
static void limiter_periodic_merged
(
 solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
 // matchRecon: 対象は**流れ 5 変数・node のみ** (plan §4.14)。化学種・受動スカラーは 0 を渡すこと
 flow_float phi_floor, int matchRecon, flow_float* Q, flow_float* limiter_Q,
 flow_float* dQdx, flow_float* dQdy, flow_float* dQdz,
 // 無次元化 Venkatakrishnan (codex plan-3 Critical 1)。既定 limScaled=0 で式は変更前と同一。
 // 化学種・受動スカラーは対象外なので 0 を渡すこと。
 int limScaled = 0, flow_float qRef = (flow_float)1.0,
 // W1h: G1 検査の変数番号 (0..4 = ro/Ux/Uy/Uz/P)。-1 = 検査しない (化学種・受動スカラー)。
 int kVar = -1
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
        cfg.limiter, msh.nCells, msh.nNormalPlanes, msh.map_plane_cells_d,
        msh.map_cell_planes_index_d, msh.map_cell_planes_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"], var.p_d["pcy"], var.p_d["pcz"],
        phi_floor, Q, s_lim_qmax, s_lim_qmin, limiter_Q, dQdx, dQdy, dQdz,
        matchRecon, (cfg.discretization == "node" ? 1 : 0), cfg.convMethod,
        limScaled, qRef,
        (limScaled == 2 ? (flow_float)cfg.venkatK
            : (flow_float)(cfg.venkatK*cfg.venkatK*cfg.venkatK/(cfg.limiterRefLength*cfg.limiterRefLength*cfg.limiterRefLength))),
        cfg.limiterLengthFromArea,
        (var.c_d.count("A_planar") ? var.c_d["A_planar"] : var.c_d["volume"]));
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    periodicGatherMinArray_d_wrapper(cfg, cuda_cfg, msh, limiter_Q);
    // W1h: ψ は gather の**後**に確定するので、G1 検査はここで別カーネルにする (plan §5.1 W1h)。
    // 流れ 5 変数のみ (kVar 0..4)。化学種・受動スカラーは kVar<0 で呼ばれないこと。
    if (cfg.limiterDiag > 0 && kVar >= 0) {
        // 許容幅の絶対床に使う基準 (limiterScaled に依らず main.cpp で計算済み)
        const flow_float qRefDiag = (kVar == 0) ? (flow_float)cfg.limiterRoRef
                                  : ((kVar == 4) ? (flow_float)cfg.limiterPRef : (flow_float)cfg.limiterARef);
        unsigned long long *g1_d=nullptr, *nf_d=nullptr, *mx_d=nullptr, *sd_d=nullptr;
        CHECK_CUDA_ERROR(cudaGetSymbolAddress((void**)&g1_d, g_limG1));
        CHECK_CUDA_ERROR(cudaGetSymbolAddress((void**)&nf_d, g_limNonFinite));
        CHECK_CUDA_ERROR(cudaGetSymbolAddress((void**)&mx_d, g_limMaxExcess));
        CHECK_CUDA_ERROR(cudaGetSymbolAddress((void**)&sd_d, g_limSides));
        limiter_g1_check_periodic_d<<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>>(
            msh.nCells, msh.nNormalPlanes, msh.map_plane_cells_d,
            msh.map_cell_planes_index_d, msh.map_cell_planes_d,
            var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
            var.p_d["pcx"], var.p_d["pcy"], var.p_d["pcz"],
            Q, s_lim_qmax, s_lim_qmin, limiter_Q, dQdx, dQdy, dQdz,
            (cfg.discretization == "node" ? 1 : 0), cfg.convMethod, kVar, qRefDiag,
            g1_d, nf_d, mx_d, sd_d);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    }
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
// --- 有界性診断 (plan convection-node-wall-reconstruction §4.20) ---
// **後処理で測れない**ので流束が使う値をカーネル内で直接数える。
// 後処理は res_*.h5 から再構成を組み直すが、保存量は更新後・勾配は更新前で 1 step ずれており、
// さらに境界条件が壁速度・壁温由来の P を書き換えるため、ro 以外は流束時の状態を再現できない
// (実測: Uy の逸脱 1.3 万〜1.9 万が全構成で出て、厳密有界のはずの Barth でも消えない = 判定不能だった)。
// ここでは psi 確定後に**流束と同じ増分**で再構成値を作り、そのノードの近傍 min/max を外れた
// face-side を変数ごとに数える。既定 off (g_limDiag=0) で atomicAdd は一切走らない。

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
 const flow_float* __restrict__ prim,  // 原始量 AoS パック [ro,Ux,Uy,Uz,P,...] (nullptr なら Q0..Q4 から gather)
 // リミッタの試行増分を流束と一致させる (plan convection-node-wall-reconstruction §4.8)。
 //   matchRecon=0: 従来 (双対面重心で g·d のみ。式は変更前と同一)
 //   matchRecon=1: 流束と同じ点・同じ形 — edgeMid なら目標点はエッジ中点、convM==2 は隣接値差の項も含む
 int matchRecon, int edgeMid, int convM,
 // 無次元化 Venkatakrishnan (plan §4.13)。scaled=0 で従来の式 (ビット変化なし)。
 //   qref[5]  : ro_ref, a_ref, a_ref, a_ref, p_ref (速度 3 成分は共通の a_ref)
 //   eps2Coef : (K / L_ref)^3 。eps2 = eps2Coef * h_i^3、h_i は下の lenArea で決まる
 //   lenArea  : 1 = h_i = sqrt(A_planar) (2D/軸対称) / 0 = cbrt(volume) (3D)
 int scaled, flow_float qr0, flow_float qr1, flow_float qr4,
 flow_float eps2Coef, int lenArea, geom_float* A_planar
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
        if (matchRecon == 0) {
            // 従来経路 (式は変更前と同一): 双対面重心で g·d を評価し、Qt を作ってから差を取る
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
        } else {
            // 流束一致経路: 目標点と増分の形を convectiveFlux 側と揃える。
            // 増分 delta を直接作る (Qt を経由しないので float32 の足して引く桁落ちが無い)。
            const geom_int ic1 = plane_cells[2*ip+0] + plane_cells[2*ip+1] - ic0;
            flow_float dcp_x, dcp_y, dcp_z;
            if (edgeMid != 0) {                      // node: 目標点 = エッジ中点 (g_reconEdgeMid と同じ)
                dcp_x = (flow_float)0.5*(ccx[ic1]-cx0);
                dcp_y = (flow_float)0.5*(ccy[ic1]-cy0);
                dcp_z = (flow_float)0.5*(ccz[ic1]-cz0);
            } else {                                  // cell: 双対面重心のままで流束と整合している
                dcp_x = pcx[ip]-cx0; dcp_y = pcy[ip]-cy0; dcp_z = pcz[ip]-cz0;
            }
            #pragma unroll
            for (int k=0;k<5;k++){
                // 流束と**同じ関数**で増分を作る (reconIncrement_d.cuh)。Qt を経由しないので桁落ちも無い
                const flow_float delta = recon_increment(convM, qc[k], Q[k][ic1],
                                                         gx[k], gy[k], gz[k], dcp_x, dcp_y, dcp_z);
                flow_float lk;
                if (scaled == 2 && SCHEME != 1) {
                    // 比の形: 無次元なので基準値も長さも要らない (venkatK を eps として使う)
                    lk = (fabsf(delta) > (flow_float)1.0e-20)
                       ? venkata_limiter_ratio(qmax[k]-qc[k], qmin[k]-qc[k], delta, eps2Coef)
                       : (flow_float)1.0;
                } else if (scaled == 1 && SCHEME != 1) {
                    // 変数ごとの固定参照で無次元化してから Venkatakrishnan
                    const flow_float qr = (k == 0) ? qr0 : ((k == 4) ? qr4 : qr1);
                    const flow_float inv = (flow_float)1.0/qr;
                    const flow_float hi = (lenArea != 0) ? sqrtf(A_planar[ic0]) : cbrtf(volume);
                    const flow_float e2 = eps2Coef * hi*hi*hi;
                    lk = venkata_limiter_scaled((qmax[k]-qc[k])*inv, (qmin[k]-qc[k])*inv, delta*inv, e2);
                } else {
                    lk = (SCHEME == 1)
                        ? barth_Jespersen_limiter(qmax[k]-qc[k], qmin[k]-qc[k], delta, volume)
                        : venkata_limiter        (qmax[k]-qc[k], qmin[k]-qc[k], delta, volume);
                }
                ltmp[k] = min(ltmp[k], lk);
            }
        }
    }

    #pragma unroll
    for (int k=0;k<5;k++) Lim[k][ic0] = min(max(ltmp[k], (flow_float)0.0), (flow_float)1.0);

    // pass3 (診断・既定 off): psi 確定後に流束と同じ再構成を作り、近傍 min/max を外れた面側を数える。
    if (g_limDiag != 0) {
        flow_float psi[5];
        #pragma unroll
        for (int k=0;k<5;k++) psi[k] = min(max(ltmp[k], (flow_float)0.0), (flow_float)1.0);
        for (geom_int ilp=index_st; ilp<index_en; ilp++) {
            geom_int ip = cell_planes[ilp];
            if (ip >= nNormalPlanes) continue;
            const geom_int ic1 = plane_cells[2*ip+0] + plane_cells[2*ip+1] - ic0;
            // 評価点は**流束側の規則だけ**で決める (codex plan-3 Critical 2)。node 流束は
            // `g_reconEdgeMid`=1 で `matchRecon` に依らず常にエッジ中点 (`convectiveFlux_d.cu:69`,
            // `convectiveFlux_slau_d.inc.cuh:134`)。診断を `matchRecon` で分岐させると
            // `mr0` の構成だけ双対面重心を測ることになり、A/B が別々の点の比較になる。
            // 増分の形 (convM) も流束と同じにする。
            flow_float dx, dy, dz;
            if (edgeMid != 0) {
                dx = (flow_float)0.5*(ccx[ic1]-cx0); dy = (flow_float)0.5*(ccy[ic1]-cy0); dz = (flow_float)0.5*(ccz[ic1]-cz0);
            } else {
                dx = pcx[ip]-cx0; dy = pcy[ip]-cy0; dz = pcz[ip]-cz0;
            }
            atomicAdd(&g_limSides, 1ULL);
            #pragma unroll
            for (int k=0;k<5;k++) {
                const flow_float d = recon_increment(convM, qc[k], Q[k][ic1], gx[k], gy[k], gz[k], dx, dy, dz);
                const flow_float qf = qc[k] + psi[k]*d;
                // 非有限は「逸脱なし」ではない (codex plan-3 Critical 2)。NaN は大小比較が両方 false に
                // なるので明示的に数える。入力 (qc/ψ/増分) の非有限も同じ枠で拾う。
                if (!isfinite(qf) || !isfinite(d) || !isfinite(psi[k]) || !isfinite(qc[k])) {
                    atomicAdd(&g_limNonFinite[k], 1ULL);
                    continue;
                }
                // 丸め許容は「近傍レンジ」と「値の大きさ」の大きい方 (自由流でレンジ 0・0 を跨ぐ速度の両方に耐える)
                // 許容幅には**構成によらない絶対床** (1e-6·q_ref) を敷く。恒等 0 の変数 (疑似 2D の Uy/Uz 等)
                // では近傍レンジも値の大きさも 0 になり、float ノイズが全部「逸脱」に化けるため。
                const flow_float mag = max(fabsf(qmax[k]), fabsf(qmin[k]));
                const flow_float qrk = (k == 0) ? qr0 : ((k == 4) ? qr4 : qr1);
                const flow_float sc  = max(max(qmax[k]-qmin[k], mag) * (flow_float)1.0e-5,
                                           qrk * (flow_float)1.0e-6);
                if (qf > qmax[k] + sc || qf < qmin[k] - sc) {
                    atomicAdd(&g_limG1[k], 1ULL);
                    // 逸脱量そのものも残す (件数 0 でも「どれだけ外れたか」を言えるように)
                    const flow_float ex = max(qf - qmax[k], qmin[k] - qf) / max(sc, (flow_float)1.0e-30);
                    atomicMax(&g_limMaxExcess[k], (unsigned long long)(ex * 1.0e3));
                }
            }
        }
    }
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

    // 有界性診断 (§4.20)。既定 off。on のとき一定間隔でカウンタを 1 行印字してリセットする。
    {
        const int diag = (cfg.limiterDiag > 0) ? 1 : 0;
        CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_limDiag, &diag, sizeof(int)));
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
        ((cfg.primPack != 0 && cfg.gradLSQ == 2) ? prim_pack_device_ptr() : nullptr), \
        cfg.limiterMatchRecon, (cfg.discretization == "node" ? 1 : 0), cfg.convMethod, \
        cfg.limiterScaled, (flow_float)cfg.limiterRoRef, (flow_float)cfg.limiterARef, (flow_float)cfg.limiterPRef, \
        (cfg.limiterScaled == 2 ? (flow_float)cfg.venkatK \
            : (flow_float)(cfg.venkatK*cfg.venkatK*cfg.venkatK/(cfg.limiterRefLength*cfg.limiterRefLength*cfg.limiterRefLength))), \
        cfg.limiterLengthFromArea, \
        (var.c_d.count("A_planar") ? var.c_d["A_planar"] : var.c_d["volume"])
    // 周期 node (合併 CV) は 2 段 (極値の group max/min → ψ の group min) で周期対の ψ を一致させる (§4.8)。
    const bool perNode = periodicNodeActive(cfg, msh);
    if (perNode) {
        const char* qn[5]  = {"ro","Ux","Uy","Uz","P"};
        const char* ln[5]  = {"limiter_ro","limiter_Ux","limiter_Uy","limiter_Uz","limiter_P"};
        const char* gxn[5] = {"drodx","dUxdx","dUydx","dUzdx","dPdx"};
        const char* gyn[5] = {"drody","dUxdy","dUydy","dUzdy","dPdy"};
        const char* gzn[5] = {"drodz","dUxdz","dUydz","dUzdz","dPdz"};
        // 無次元化の基準は通常経路 (limiter_r1_fused5_d) と同じ割り当て: k=0→ρ, k=4→P, 速度→音速
        const flow_float qref5[5] = {(flow_float)cfg.limiterRoRef, (flow_float)cfg.limiterARef,
                                     (flow_float)cfg.limiterARef,  (flow_float)cfg.limiterARef,
                                     (flow_float)cfg.limiterPRef};
        for (int k = 0; k < 5; ++k)
            limiter_periodic_merged<false>(cfg, cuda_cfg, msh, var, 0.0f, cfg.limiterMatchRecon,
                var.c_d[qn[k]], var.c_d[ln[k]], var.c_d[gxn[k]], var.c_d[gyn[k]], var.c_d[gzn[k]],
                cfg.limiterScaled, qref5[k], k);
    } else if (cfg.limiter == 1)
        limiter_r1_fused5_d<1><<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>> (FORGE_LIMITER_FUSED5_ARGS);
    else
        limiter_r1_fused5_d<2><<<cuda_cfg.dimGrid_normalcell_small , cuda_cfg.dimBlock_small>>> (FORGE_LIMITER_FUSED5_ARGS);
    #undef FORGE_LIMITER_FUSED5_ARGS

    if (cfg.limiterDiag > 0) {
        static int s_lim_call = 0;
        const int interval = cfg.limiterDiag;
        if ((s_lim_call % interval) == 0) {
            gpuErrchkKernelSync();
            unsigned long long g1[5] = {0,0,0,0,0}, sides = 0;
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(g1,    g_limG1,    5*sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(&sides, g_limSides, sizeof(unsigned long long)));
            unsigned long long nf[5] = {0,0,0,0,0}, ex[5] = {0,0,0,0,0};
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(nf, g_limNonFinite, 5*sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyFromSymbol(ex, g_limMaxExcess, 5*sizeof(unsigned long long)));
            // 検査件数 0 は「逸脱なし」ではない (周期 node には pass3 が無い等)。明示的に FAIL と印字する。
            const char* verdict = (sides == 0ULL) ? " VERDICT=NO-CHECKS(FAIL)"
                : ((g1[0]|g1[1]|g1[2]|g1[3]|g1[4]|nf[0]|nf[1]|nf[2]|nf[3]|nf[4]) ? " VERDICT=VIOLATIONS" : " VERDICT=CLEAN");
            printf("LIMG1 call=%d sides=%llu out[ro=%llu Ux=%llu Uy=%llu Uz=%llu P=%llu]"
                   " nonfinite[%llu %llu %llu %llu %llu] maxexcess_x1e3[%llu %llu %llu %llu %llu]%s\n",
                   s_lim_call, sides, g1[0], g1[1], g1[2], g1[3], g1[4],
                   nf[0], nf[1], nf[2], nf[3], nf[4], ex[0], ex[1], ex[2], ex[3], ex[4], verdict);
            const unsigned long long z5b[5] = {0,0,0,0,0};
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_limNonFinite, z5b, 5*sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_limMaxExcess, z5b, 5*sizeof(unsigned long long)));
            const unsigned long long z5[5] = {0,0,0,0,0}, z = 0ULL;
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_limG1,    z5, 5*sizeof(unsigned long long)));
            CHECK_CUDA_ERROR(cudaMemcpyToSymbol(g_limSides, &z, sizeof(unsigned long long)));
        }
        s_lim_call++;
    }

    // 多成分 face 整合再構成: 各化学種 Y_s に Venkat リミタ ψ_Y を計算 (∇Y は speciesGradient 済)。
    // speciesFaceReconstruction==1 のみ。flux では min(ψ_ρ, ψ_Y) を Y 再構成に使う (boundedness)。
    if (cfg.speciesFaceReconstruction >= 1 && var.nSpeciesRegistered >= 2) {
        for (int s = 0; s < var.nSpeciesRegistered; ++s) {
            const std::string i = std::to_string(s);
            fill_limiter_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(var.c_d["limiter_Y"+i], msh.nCells_all, 1.0f);
            if (perNode) {
                // 化学種は対象外 (matchRecon=0)。対象は流れ 5 変数・node のみ (plan §4.14)
                limiter_periodic_merged<false>(cfg, cuda_cfg, msh, var, 0.0f, 0,
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
            // 受動スカラーは対象外 (matchRecon=0)
            limiter_periodic_merged<true>(cfg, cuda_cfg, msh, var, static_cast<flow_float>(1.0e-30), 0,
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
