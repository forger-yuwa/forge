#include "timeIntegration_d.cuh"
#include "weakIsothermalWall_d.cuh"
#include "lowMachPrecond_d.cuh"   // Phase 4: β (lowMachBeta)・c' (lowMachCprime) device ヘルパ
#include "cuda_forge/eos_jacobian_d.cuh"  // 一般EOS固有系 (eos_split_jacobian_general_closed)。precond 経路で使用
#include "cuda_forge/block_dplur_jacobian_d.cuh"  // block_dplur::accumulate_split_jacobian_cf (共有ヘッダ)
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cooperative_groups.h>
#include <cstring>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

// 診断 (env FORGE_AXIS_DIAG_ALPHA, 既定 0=不変): 近軸 (r→0) の半径方向音響モード安定化。
// roUy 対角に α·A_planar·c を加える。revolved 軸面積 (r_f·S→0) が落とす半径音響スペクトル半径を
// planar 面積 A_pl で補うもの。near-axis radial-momentum 不安定 (case/28 TP) の最小修正診断。
__device__ float g_axisDiagAlpha = 0.0f;

namespace block_dplur {

template<typename T>
__device__ __forceinline__ void zero5(T* vec)
{
    #pragma unroll
    for (int i = 0; i < 5; ++i) {
        vec[i] = 0.0;
    }
}

template<typename T>
__device__ __forceinline__ void add_identity_scaled(T mat[5][5], T scale)
{
    #pragma unroll
    for (int row = 0; row < 5; ++row) {
        mat[row][row] += scale;
    }
}

// accumulate_split_jacobian_cf は cuda_forge/block_dplur_jacobian_d.cuh へ移設 (host/device 共有・Level3
// 単体テスト tools/test_eos_jacobian.cpp から呼ぶため)。block_dplur:: で参照する。

__device__ __forceinline__ void multiply_add_5x5_vec(
    const flow_float mat[5][5],
    const flow_float vec[5],
    flow_float out[5]
)
{
    #pragma unroll
    for (int row = 0; row < 5; ++row) {
        flow_float sum = 0.0;
        #pragma unroll
        for (int col = 0; col < 5; ++col) {
            sum += mat[row][col] * vec[col];
        }
        out[row] += sum;
    }
}

__device__ __forceinline__ void copy5(const flow_float* src, flow_float* dst)
{
    #pragma unroll
    for (int i = 0; i < 5; ++i) {
        dst[i] = src[i];
    }
}

template<typename T>
__device__ __forceinline__ void zero5x5(T mat[5][5])
{
    #pragma unroll
    for (int row = 0; row < 5; ++row) {
        #pragma unroll
        for (int col = 0; col < 5; ++col) {
            mat[row][col] = 0.0;
        }
    }
}

__device__ __forceinline__ void add_scaled_5x5(flow_float dst[5][5], const flow_float src[5][5], flow_float scale)
{
    #pragma unroll
    for (int row = 0; row < 5; ++row) {
        #pragma unroll
        for (int col = 0; col < 5; ++col) {
            dst[row][col] += scale * src[row][col];
        }
    }
}

__device__ __forceinline__ void build_abs_jacobian(
    flow_float gamma,
    flow_float nx,
    flow_float ny,
    flow_float nz,
    flow_float u,
    flow_float v,
    flow_float w,
    flow_float H,
    flow_float c,
    flow_float abs_jac[5][5]
)
{
    zero5x5(abs_jac);

    const flow_float sonic = max(c, static_cast<flow_float>(1.0e-8));
    const flow_float ek = 0.5 * (u * u + v * v + w * w);
    const flow_float U = u * nx + v * ny + w * nz;
    const flow_float chi = (gamma - 1.0) / sonic;
    const flow_float inv_sqrt2 = static_cast<flow_float>(0.7071067811865475244);

    flow_float lambda[5] = {
        fabs(U + sonic),
        fabs(U),
        fabs(U),
        fabs(U),
        fabs(U - sonic)
    };

    flow_float R[5][5];
    flow_float L[5][5];

    R[0][0] = inv_sqrt2 / sonic;             R[0][1] = ny / sonic;               R[0][2] = nz / sonic;               R[0][3] = nx / sonic;               R[0][4] = inv_sqrt2 / sonic;
    R[1][0] = (u / sonic + nx) * inv_sqrt2; R[1][1] = u * ny / sonic + nz;      R[1][2] = u * nz / sonic - ny;      R[1][3] = u * nx / sonic;           R[1][4] = (u / sonic - nx) * inv_sqrt2;
    R[2][0] = (v / sonic + ny) * inv_sqrt2; R[2][1] = v * ny / sonic;           R[2][2] = v * nz / sonic + nx;      R[2][3] = v * nx / sonic - nz;      R[2][4] = (v / sonic - ny) * inv_sqrt2;
    R[3][0] = (w / sonic + nz) * inv_sqrt2; R[3][1] = w * ny / sonic - nx;      R[3][2] = w * nz / sonic;           R[3][3] = w * nx / sonic + ny;      R[3][4] = (w / sonic - nz) * inv_sqrt2;
    R[4][0] = (ek / sonic + 1.0 / chi + U) * inv_sqrt2;
    R[4][1] = ek * ny / sonic + nz * u - nx * w;
    R[4][2] = ek * nz / sonic + nx * v - ny * u;
    R[4][3] = ek * nx / sonic + ny * w - nz * v;
    R[4][4] = (ek / sonic + 1.0 / chi - U) * inv_sqrt2;

    L[0][0] = (chi * ek - U) * inv_sqrt2;   L[0][1] = (-chi * u + nx) * inv_sqrt2; L[0][2] = (-chi * v + ny) * inv_sqrt2; L[0][3] = (-chi * w + nz) * inv_sqrt2; L[0][4] = chi * inv_sqrt2;
    L[1][0] = ny * (-chi * ek + sonic) - nz * u + nx * w;
    L[1][1] = ny * chi * u + nz;            L[1][2] = ny * chi * v;              L[1][3] = ny * chi * w - nx;         L[1][4] = -ny * chi;
    L[2][0] = nz * (-chi * ek + sonic) - nx * v + ny * u;
    L[2][1] = nz * chi * u - ny;            L[2][2] = nz * chi * v + nx;         L[2][3] = nz * chi * w;              L[2][4] = -nz * chi;
    L[3][0] = nx * (-chi * ek + sonic) - ny * w + nz * v;
    L[3][1] = nx * chi * u;                 L[3][2] = nx * chi * v - nz;         L[3][3] = nx * chi * w + ny;         L[3][4] = -nx * chi;
    L[4][0] = (chi * ek + U) * inv_sqrt2;   L[4][1] = (-chi * u - nx) * inv_sqrt2; L[4][2] = (-chi * v - ny) * inv_sqrt2; L[4][3] = (-chi * w - nz) * inv_sqrt2; L[4][4] = chi * inv_sqrt2;

    #pragma unroll
    for (int row = 0; row < 5; ++row) {
        #pragma unroll
        for (int col = 0; col < 5; ++col) {
            flow_float sum = 0.0;
            #pragma unroll
            for (int k = 0; k < 5; ++k) {
                sum += R[row][k] * lambda[k] * L[k][col];
            }
            abs_jac[row][col] = sum;
        }
    }
}

// LU-SGS の通量分割 Jacobian を同時構築する。
//   a_plus = A^+ = R Λ^+ L,  k_off = -A^- = R(-Λ^-)L = ½(|A|-A)。
// 検証済みの一般EOS閉形式 eos_split_jacobian_general_closed (eos_jacobian_d.cuh, Level1〜3 検証済) を流用し
// CPG/TP を統一 (旧来の CPG 専用べた打ち R/L=RL≠I の近似を廃止)。H は実全エンタルピー Ht[ic]、
// κ=γ−1、χ_eos=c²−κh。CPG では χ_eos≈0 で従来の CPG 固有系に簡約 (収束先は不変・厳密 Jacobian で僅かに向上)。
// double で組み立て float へ格納 (precond カーネルは元々 double 相当)。標準経路 accumulate_split_jacobian_cf は
// 別実装 (CPG ビット不変) のまま。
__device__ __forceinline__ void build_jacobian_split(
    flow_float gamma,
    flow_float nx, flow_float ny, flow_float nz,
    flow_float u, flow_float v, flow_float w,
    flow_float H, flow_float c,
    flow_float a_plus[5][5], flow_float k_off[5][5]
)
{
    const double sonic = (double)max(c, static_cast<flow_float>(1.0e-8));
    const double ek    = 0.5*((double)u*u + (double)v*v + (double)w*w);
    const double kappa = (double)gamma - 1.0;
    const double chi   = sonic*sonic - kappa*((double)H - ek);   // χ_eos = c²−κh (CPG で ≈0)
    double Ap[5][5], Am[5][5];
    eos_split_jacobian_general_closed((double)u,(double)v,(double)w, (double)nx,(double)ny,(double)nz,
                                      sonic, (double)H, kappa, chi, Ap, Am);   // Ap=A⁺, Am=A⁻
    #pragma unroll
    for (int i=0;i<5;++i)
        #pragma unroll
        for (int j=0;j<5;++j){ a_plus[i][j]=(flow_float)Ap[i][j]; k_off[i][j]=(flow_float)(-Am[i][j]); }
}

template<typename T>
__device__ __forceinline__ bool solve_5x5(T mat[5][5], T rhs[5], T sol[5])
{
    #pragma unroll
    for (int col = 0; col < 5; ++col) {
        int pivot = col;
        T pivot_abs = fabs(mat[col][col]);
        #pragma unroll
        for (int row = col + 1; row < 5; ++row) {
            const T candidate = fabs(mat[row][col]);
            if (candidate > pivot_abs) {
                pivot = row;
                pivot_abs = candidate;
            }
        }

        if (pivot_abs < static_cast<T>(1.0e-20)) {
            zero5(sol);
            return false;
        }

        if (pivot != col) {
            #pragma unroll
            for (int k = 0; k < 5; ++k) {
                const T tmp = mat[col][k];
                mat[col][k] = mat[pivot][k];
                mat[pivot][k] = tmp;
            }
            const T rhs_tmp = rhs[col];
            rhs[col] = rhs[pivot];
            rhs[pivot] = rhs_tmp;
        }

        const T inv_pivot = static_cast<T>(1.0) / mat[col][col];
        #pragma unroll
        for (int row = col + 1; row < 5; ++row) {
            const T factor = mat[row][col] * inv_pivot;
            mat[row][col] = 0.0;
            #pragma unroll
            for (int k = col + 1; k < 5; ++k) {
                mat[row][k] -= factor * mat[col][k];
            }
            rhs[row] -= factor * rhs[col];
        }
    }

    for (int row = 4; row >= 0; --row) {
        T sum = rhs[row];
        #pragma unroll
        for (int col = row + 1; col < 5; ++col) {
            sum -= mat[row][col] * sol[col];
        }
        sol[row] = sum / mat[row][row];
    }

    return true;
}

// 5×5 を 2 つの RHS について同時に解く (部分ピボット Gauss 消去を 1 回・float)。
// Phase 4 (lowMachPrecond=2) の Sherman-Morrison 解法で D0⁻¹b と D0⁻¹g を同じ分解で得るのに使う。
// D0 は物理ブロックで良条件 (既存 0/1 カーネルが float で解いているのと同形) ゆえ float で十分。
__device__ __forceinline__ bool solve_5x5_2rhs(flow_float mat[5][5],
                                               flow_float r1[5], flow_float r2[5],
                                               flow_float s1[5], flow_float s2[5])
{
    #pragma unroll
    for (int col = 0; col < 5; ++col) {
        int pivot = col;
        flow_float pivot_abs = fabs(mat[col][col]);
        #pragma unroll
        for (int row = col + 1; row < 5; ++row) {
            const flow_float candidate = fabs(mat[row][col]);
            if (candidate > pivot_abs) { pivot = row; pivot_abs = candidate; }
        }

        if (pivot_abs < static_cast<flow_float>(1.0e-20)) {
            zero5(s1); zero5(s2);
            return false;
        }

        if (pivot != col) {
            #pragma unroll
            for (int k = 0; k < 5; ++k) {
                const flow_float tmp = mat[col][k]; mat[col][k] = mat[pivot][k]; mat[pivot][k] = tmp;
            }
            flow_float t = r1[col]; r1[col] = r1[pivot]; r1[pivot] = t;
            t = r2[col]; r2[col] = r2[pivot]; r2[pivot] = t;
        }

        const flow_float inv_pivot = static_cast<flow_float>(1.0) / mat[col][col];
        #pragma unroll
        for (int row = col + 1; row < 5; ++row) {
            const flow_float factor = mat[row][col] * inv_pivot;
            mat[row][col] = 0.0;
            #pragma unroll
            for (int k = col + 1; k < 5; ++k) mat[row][k] -= factor * mat[col][k];
            r1[row] -= factor * r1[col];
            r2[row] -= factor * r2[col];
        }
    }

    for (int row = 4; row >= 0; --row) {
        flow_float sum1 = r1[row], sum2 = r2[row];
        #pragma unroll
        for (int col = row + 1; col < 5; ++col) { sum1 -= mat[row][col] * s1[col]; sum2 -= mat[row][col] * s2[col]; }
        const flow_float inv = static_cast<flow_float>(1.0) / mat[row][row];
        s1[row] = sum1 * inv;
        s2[row] = sum2 * inv;
    }

    return true;
}

__device__ __forceinline__ void load_block_vec(
    geom_int ic,
    flow_float* v0,
    flow_float* v1,
    flow_float* v2,
    flow_float* v3,
    flow_float* v4,
    flow_float out[5]
)
{
    out[0] = v0[ic];
    out[1] = v1[ic];
    out[2] = v2[ic];
    out[3] = v3[ic];
    out[4] = v4[ic];
}

__device__ __forceinline__ void store_block_vec(
    geom_int ic,
    const flow_float in[5],
    flow_float* v0,
    flow_float* v1,
    flow_float* v2,
    flow_float* v3,
    flow_float* v4
)
{
    v0[ic] = in[0];
    v1[ic] = in[1];
    v2[ic] = in[2];
    v3[ic] = in[3];
    v4[ic] = in[4];
}

}

__global__ void runge_kutta_exp_4th_d
// see https://sci-hub.se/https://doi.org/10.1016/j.compfluid.2003.10.004
// N: previous outer step , M: previous inner loop
( 
 int loop, 
 flow_float coef_DT,
 flow_float coef_Res,

 flow_float dt ,
 flow_float* dt_local ,

 // mesh structure
 geom_int nCells_all , geom_int nCells,
 geom_float* vol ,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* roN ,
 flow_float* roUxN ,
 flow_float* roUyN ,
 flow_float* roUzN ,
 flow_float* roeN ,

 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,

 flow_float* res_ro,
 flow_float* res_roUx,
 flow_float* res_roUy,
 flow_float* res_roUz,
 flow_float* res_roe,

 flow_float* res_ro_m,
 flow_float* res_roUx_m,
 flow_float* res_roUy_m,
 flow_float* res_roUz_m,
 flow_float* res_roe_m

)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    geom_float v = vol[ic];

    if (ic < nCells) {
        flow_float dt_l = dt_local[ic];

        if (loop == 0) {
            res_ro_m[ic]   = 0.0;
            res_roUx_m[ic] = 0.0;
            res_roUy_m[ic] = 0.0;
            res_roUz_m[ic] = 0.0;
            res_roe_m[ic]  = 0.0;
        }
        // N: previous outer step , M: previous inner loop
        res_ro_m[ic]   += coef_Res*res_ro[ic]*dt_l/v;
        res_roUx_m[ic] += coef_Res*res_roUx[ic]*dt_l/v;
        res_roUy_m[ic] += coef_Res*res_roUy[ic]*dt_l/v;
        res_roUz_m[ic] += coef_Res*res_roUz[ic]*dt_l/v;
        res_roe_m[ic]  += coef_Res*res_roe[ic]*dt_l/v;

        if (loop < 3) {
            ro[ic]   = roN[ic]   +coef_DT*res_ro[ic]*dt_l/v;
            roUx[ic] = roUxN[ic] +coef_DT*res_roUx[ic]*dt_l/v;
            roUy[ic] = roUyN[ic] +coef_DT*res_roUy[ic]*dt_l/v;
            roUz[ic] = roUzN[ic] +coef_DT*res_roUz[ic]*dt_l/v;
            roe[ic]  = roeN[ic]  +coef_DT*res_roe[ic]*dt_l/v;
        } else {
            res_ro[ic]   = res_ro_m[ic]*v/dt_l ;
            res_roUx[ic] = res_roUx_m[ic]*v/dt_l ;
            res_roUy[ic] = res_roUy_m[ic]*v/dt_l ;
            res_roUz[ic] = res_roUz_m[ic]*v/dt_l ;
            res_roe[ic]  = res_roe_m[ic]*v/dt_l ;

            ro[ic]   = roN[ic]  + res_ro_m[ic] ;
            roUx[ic] = roUxN[ic]+ res_roUx_m[ic] ;
            roUy[ic] = roUyN[ic]+ res_roUy_m[ic] ;
            roUz[ic] = roUzN[ic]+ res_roUz_m[ic] ;
            roe[ic]  = roeN[ic] + res_roe_m[ic] ;
        }
    }
}

__global__ void runge_kutta_exp_d
// see https://sci-hub.se/https://doi.org/10.1016/j.compfluid.2003.10.004
// N: previous outer step , M: previous inner loop
( 
 int loop,
 flow_float coef_N,
 flow_float coef_M,
 flow_float coef_Res,

 flow_float dt ,
 flow_float* dt_local ,

 // mesh structure
 geom_int nCells_all , geom_int nCells,
 geom_float* vol ,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* roN ,
 flow_float* roUxN ,
 flow_float* roUyN ,
 flow_float* roUzN ,
 flow_float* roeN ,

 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,

 flow_float* res_ro,
 flow_float* res_roUx,
 flow_float* res_roUy,
 flow_float* res_roUz,
 flow_float* res_roe
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        geom_float v = vol[ic];
        flow_float dt_l = dt_local[ic];
        // N: previous outer step , M: previous inner loop
        ro[ic]   = coef_N*roN[ic]   + coef_M*roM[ic]   + coef_Res*res_ro[ic]*dt_l/v;
        roUx[ic] = coef_N*roUxN[ic] + coef_M*roUxM[ic] + coef_Res*res_roUx[ic]*dt_l/v;
        roUy[ic] = coef_N*roUyN[ic] + coef_M*roUyM[ic] + coef_Res*res_roUy[ic]*dt_l/v;
        roUz[ic] = coef_N*roUzN[ic] + coef_M*roUzM[ic] + coef_Res*res_roUz[ic]*dt_l/v;
        roe[ic]  = coef_N*roeN[ic]  + coef_M*roeM[ic]  + coef_Res*res_roe[ic]*dt_l/v;
    }
}
__global__ void implicit_defect_correction_d
(
 int loop,
 flow_float dt,
 flow_float* dt_local,
 flow_float implicit_relax,
 flow_float* gamma_arr,   // per-cell γ (TP: γ_mix(T), CPG: cfg.gamma)。軸対称ソース Jacobian 用

 // mesh structure
 geom_int nCells_all , geom_int nCells,
 geom_float* vol,
 geom_int* plane_cells,
 geom_int* cell_planes_index,
 geom_int* cell_planes,
 geom_float* ccx,
 geom_float* ccy,
 geom_float* ccz,
 // 面ごとの差 e = cc[ic1] − cc[ic0] (ic0/ic1 = plane_cells[2*ip+0/1]、double の座標から 1 回だけ丸めた値、var.p_d["ge_*"]。
 // plans/active/architecture-float-state-double-geometry.md §4.2a、段 ③)。粘性の対角の cc[other] − cc[ic] に使う。
 const flow_float* ge_x,
 const flow_float* ge_y,
 const flow_float* ge_z,
 geom_float* sx,
 geom_float* sy,
 geom_float* sz,
 geom_float* ss,

 // variables
 flow_float* ro,
 flow_float* roUx,
 flow_float* roUy,
 flow_float* roUz,
 flow_float* roe,

 flow_float* roN,
 flow_float* roUxN,
 flow_float* roUyN,
 flow_float* roUzN,
 flow_float* roeN,

 flow_float laminar_visc,
 flow_float* vis_turb,
 flow_float* sonic,
 flow_float* Ux,
 flow_float* Uy,
 flow_float* Uz,

 flow_float* res_ro,
 flow_float* res_roUx,
 flow_float* res_roUy,
 flow_float* res_roUz,
 flow_float* res_roe,

 flow_float* corr_ro_old,
 flow_float* corr_roUx_old,
 flow_float* corr_roUy_old,
 flow_float* corr_roUz_old,
 flow_float* corr_roe_old,

 flow_float* corr_ro_new,
 flow_float* corr_roUx_new,
 flow_float* corr_roUy_new,
 flow_float* corr_roUz_new,
 flow_float* corr_roe_new,

 // 軸対称ソース項 Jacobian 用 (CPG/TP 共通)
 int isAxisymmetric,
 flow_float* A_planar,
 flow_float axisRFloor,
 // dual-time 物理時間項の対角係数 a/Δt（定常は 0）
 flow_float unsteady_diag
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        geom_float v = vol[ic];
        flow_float dt_l = dt_local[ic];
        const flow_float density = max(ro[ic], static_cast<flow_float>(1.0e-30));
        const flow_float velocity_x = Ux[ic];
        const flow_float velocity_y = Uy[ic];
        const flow_float velocity_z = Uz[ic];
        const flow_float local_sonic = max(sonic[ic], static_cast<flow_float>(0.0));
        const flow_float nu_eff = (laminar_visc + max(vis_turb[ic], static_cast<flow_float>(0.0))) / density;

        if (loop == 0) {
            corr_ro_old[ic] = 0.0;
            corr_roUx_old[ic] = 0.0;
            corr_roUy_old[ic] = 0.0;
            corr_roUz_old[ic] = 0.0;
            corr_roe_old[ic] = 0.0;
        }

        flow_float diag_face_sum = 0.0;
        flow_float neighbor_ro = 0.0;
        flow_float neighbor_roUx = 0.0;
        flow_float neighbor_roUy = 0.0;
        flow_float neighbor_roUz = 0.0;
        flow_float neighbor_roe = 0.0;
        const geom_int plane_begin = cell_planes_index[ic];
        const geom_int plane_end = cell_planes_index[ic + 1];
        for (geom_int plane_offset = plane_begin; plane_offset < plane_end; ++plane_offset) {
            const geom_int ip = cell_planes[plane_offset];
            const flow_float face_area = max(ss[ip], static_cast<flow_float>(1.0e-30));
            const flow_float advective_radius = fabs(
                velocity_x * sx[ip] + velocity_y * sy[ip] + velocity_z * sz[ip]
            ) / face_area + local_sonic;

            const geom_int ic0 = plane_cells[2 * ip + 0];
            const geom_int ic1 = plane_cells[2 * ip + 1];
            const geom_int other_ic = (ic0 == ic) ? ic1 : ic0;
            // cc[other] − cc[ic]: ic が ic0 側なら e、ic1 側なら −e (符号反転は厳密。§4.2a、段 ③)
            const flow_float dcc_x = (ic0 == ic) ? ge_x[ip] : -ge_x[ip];
            const flow_float dcc_y = (ic0 == ic) ? ge_y[ip] : -ge_y[ip];
            const flow_float dcc_z = (ic0 == ic) ? ge_z[ip] : -ge_z[ip];
            const flow_float dcc = max(
                sqrt(dcc_x * dcc_x + dcc_y * dcc_y + dcc_z * dcc_z),
                static_cast<flow_float>(1.0e-30)
            );
            const flow_float dcc_dot_s = max(
                fabs(dcc_x * sx[ip] + dcc_y * sy[ip] + dcc_z * sz[ip]),
                static_cast<flow_float>(1.0e-30)
            );
            const flow_float delta = max(dcc * face_area * face_area / dcc_dot_s, static_cast<flow_float>(1.0e-30));
            // 粘性対角は residual の粘性流束 Jacobian (2ν·ss²/dcc_dot_s = 2ν·delta/dcc) と整合させる。
            // 旧 face_area·(2ν/delta) は面積 delta を長さ²扱いし ≈2ν に潰れ、軸対称近軸で r 重みを失い、
            // ゼロ面積(対称)面にもスプリアス項を載せていた (residual は ip<nNormalPlanes で軸面を除外)。
            const flow_float viscous_diag = static_cast<flow_float>(2.0) * nu_eff * delta / dcc;
            const flow_float face_coeff = face_area * advective_radius + viscous_diag;
            const flow_float offdiag_coeff = static_cast<flow_float>(0.5) * face_coeff;
            diag_face_sum += face_coeff;

            if (other_ic < nCells) {
                neighbor_ro += offdiag_coeff * corr_ro_old[other_ic];
                neighbor_roUx += offdiag_coeff * corr_roUx_old[other_ic];
                neighbor_roUy += offdiag_coeff * corr_roUy_old[other_ic];
                neighbor_roUz += offdiag_coeff * corr_roUz_old[other_ic];
                neighbor_roe += offdiag_coeff * corr_roe_old[other_ic];
            }
        }

        const flow_float diag = max(
            static_cast<flow_float>(v / max(dt_l, static_cast<flow_float>(1.0e-30)))
                + static_cast<flow_float>(v) * unsteady_diag + diag_face_sum,
            static_cast<flow_float>(1.0e-30)
        );
        const flow_float inv_diag = 1.0 / diag;

        // 軸対称フープ源 (res_roUy += (P-τθθ)·A_planar, axisymmetricSource_d.cu) の Jacobian 対角成分を
        // roUy 方程式の対角に陰化する。これが無いと近軸の剛フープ源 (∝1/r) が陽 (lagged) 扱いになり、
        // block-DPLUR が安定な CFL でも scalar が発散する (block は diag_block[2][2] で陰化済。
        // 切り分けは case/29 README / plan time_integration-scalar-dplur-axisym-source.md)。
        // block 版 diag_block[2][2] と同形: A_pl·((γ-1)u_y + 2μ/(ρ r_eff))。
        // γ は per-cell gamma_arr[ic] (TP=γ_mix(T) / CPG=cfg.gamma) を使い thermally perfect でも整合。
        // scalar 対角の正値性 (対角優位) を保つため非負側のみ加える (defect-correction の不動点は不変)。
        flow_float diag_roUy = diag;
        if (isAxisymmetric == 1 &&
            !(axisRFloor > (flow_float)0.0 && ccy[ic] < axisRFloor)) {
            const flow_float A_pl = max(A_planar[ic], static_cast<flow_float>(1.0e-30));
            const flow_float r_eff = max(static_cast<flow_float>(v) / A_pl, static_cast<flow_float>(1.0e-30));
            const flow_float g1 = gamma_arr[ic] - static_cast<flow_float>(1.0);
            const flow_float mu_total = laminar_visc + max(vis_turb[ic], static_cast<flow_float>(0.0));
            const flow_float hoop = static_cast<flow_float>(2.0) * mu_total / (density * r_eff);
            const flow_float src_diag = A_pl * (g1 * velocity_y + hoop);
            diag_roUy = diag + max(src_diag, static_cast<flow_float>(0.0));
        }
        const flow_float inv_diag_roUy = 1.0 / diag_roUy;

        const flow_float jacobi_ro = (res_ro[ic] + neighbor_ro) * inv_diag;
        const flow_float jacobi_roUx = (res_roUx[ic] + neighbor_roUx) * inv_diag;
        const flow_float jacobi_roUy = (res_roUy[ic] + neighbor_roUy) * inv_diag_roUy;
        const flow_float jacobi_roUz = (res_roUz[ic] + neighbor_roUz) * inv_diag;
        const flow_float jacobi_roe = (res_roe[ic] + neighbor_roe) * inv_diag;

        corr_ro_new[ic] = implicit_relax * jacobi_ro;
        corr_roUx_new[ic] = implicit_relax * jacobi_roUx;
        corr_roUy_new[ic] = implicit_relax * jacobi_roUy;
        corr_roUz_new[ic] = implicit_relax * jacobi_roUz;
        corr_roe_new[ic] = implicit_relax * jacobi_roe;
    }
}

// 5×5 行列を複数保持しレジスタ消費が大きいため、block 上限を超えないよう __launch_bounds__ で
// 1 block あたりスレッド数を 128 に制限する（起動時の "too many resources" を回避）。
#define BLOCK_DPLUR_THREADS 128
// 占有率実験用: __launch_bounds__ の最小常駐ブロック数 (既定 1 = 従来どおり制限なし, 128 regs → 占有率 ~28 %)。
// -DBLOCK_DPLUR_MINBLOCKS=n でビルドすると regs ≤ 65536/(128 n) に制限され (spill と引き換えに) 占有率が上がる。
#ifndef BLOCK_DPLUR_MINBLOCKS
#define BLOCK_DPLUR_MINBLOCKS 1
#endif
// block-DPLUR の閉形式 FVS 版。線形 solve の内部精度を ST (float 既定 / double で軸対称近軸を根治) で
// テンプレート化。残差/状態 (flow_float=float) を ST へキャストして取り込み、R/L を作らず閉形式で
// diag/nbr を畳み、ST で in-place 5×5 solve、補正を float dq_new へ書戻す (混合精度 iterative refinement)。
// 詳細: plans/archived/precision-mixed-axisym.md。
template<typename ST>
__global__ void __launch_bounds__(BLOCK_DPLUR_THREADS, BLOCK_DPLUR_MINBLOCKS) implicit_defect_correction_block_d
(
 int loop,
 flow_float dt,
 const flow_float* __restrict__ dt_local,
 flow_float implicit_relax,
 const flow_float* __restrict__ gamma_arr,   // per-cell γ (TP: γ_mix(T), CPG: cfg.gamma)。frozen-coefficient Jacobian 用
 int thermallyPerfect,    // 1: TP 固有系 (実 Ht・χ_eos=c²−κh, κ=γ−1), 0: CPG 閉形式 (従来・ビット不変)

 geom_int nCells_all , geom_int nCells,
 const geom_float* __restrict__ vol,
 const geom_int* __restrict__ plane_cells,
 const geom_int* __restrict__ cell_planes_index,
 const geom_int* __restrict__ cell_planes,
 const geom_float* __restrict__ ccx,
 const geom_float* __restrict__ ccy,
 const geom_float* __restrict__ ccz,
 // 面ごとの差 e = cc[ic1] − cc[ic0] (ic0/ic1 = plane_cells[2*ip+0/1]、double の座標から 1 回だけ丸めた値、var.p_d["ge_*"]。
 // plans/active/architecture-float-state-double-geometry.md §4.2a、段 ③)。粘性の対角の cc[other] − cc[ic] に使う。
 const flow_float* __restrict__ ge_x,
 const flow_float* __restrict__ ge_y,
 const flow_float* __restrict__ ge_z,
 const geom_float* __restrict__ sx,
 const geom_float* __restrict__ sy,
 const geom_float* __restrict__ sz,
 const geom_float* __restrict__ ss,

 const flow_float* __restrict__ ro,
 const flow_float* __restrict__ roUx,
 const flow_float* __restrict__ roUy,
 const flow_float* __restrict__ roUz,
 const flow_float* __restrict__ roe,

 flow_float laminar_visc,
 const flow_float* __restrict__ vis_turb,
 const flow_float* __restrict__ sonic,
 const flow_float* __restrict__ Ux,
 const flow_float* __restrict__ Uy,
 const flow_float* __restrict__ Uz,
 const flow_float* __restrict__ Ht,

 const flow_float* __restrict__ res_ro,
 const flow_float* __restrict__ res_roUx,
 const flow_float* __restrict__ res_roUy,
 const flow_float* __restrict__ res_roUz,
 const flow_float* __restrict__ res_roe,

 const flow_float* __restrict__ dq_old_0,
 const flow_float* __restrict__ dq_old_1,
 const flow_float* __restrict__ dq_old_2,
 const flow_float* __restrict__ dq_old_3,
 const flow_float* __restrict__ dq_old_4,

 flow_float* dq_new_0,
 flow_float* dq_new_1,
 flow_float* dq_new_2,
 flow_float* dq_new_3,
 flow_float* dq_new_4,

 flow_float* rhs_0,
 flow_float* rhs_1,
 flow_float* rhs_2,
 flow_float* rhs_3,
 flow_float* rhs_4,

 flow_float* diag_00, flow_float* diag_01, flow_float* diag_02, flow_float* diag_03, flow_float* diag_04,
 flow_float* diag_10, flow_float* diag_11, flow_float* diag_12, flow_float* diag_13, flow_float* diag_14,
 flow_float* diag_20, flow_float* diag_21, flow_float* diag_22, flow_float* diag_23, flow_float* diag_24,
 flow_float* diag_30, flow_float* diag_31, flow_float* diag_32, flow_float* diag_33, flow_float* diag_34,
 flow_float* diag_40, flow_float* diag_41, flow_float* diag_42, flow_float* diag_43, flow_float* diag_44,

 // 軸対称ソースヤコビアン用（isAxisymmetric==1 のときのみ使用）
 int isAxisymmetric,
 const flow_float* __restrict__ A_planar,

 // 軸対称 r 床 (axisymMethod==0): ccy < axisRFloor の帯は hoop ソース不課につき Jacobian も課さない。
 flow_float axisRFloor,

 // dual-time 物理時間項の対角係数 a/Δt（定常は 0）
 flow_float unsteady_diag,

 // node-centered 軸対称: 軸上 CV で半径方向運動量 (roUy, index2) 行を decouple する (nullptr 可)。
 // SU2 流の対称面を Jacobian 内で課す = solve の外で状態を手術せず一貫して dq_roUy=0 を得る。
 const geom_int* __restrict__ axis_flag,     // (未使用: 旧 nodeAxisDirichlet の全 5 行 decouple。常に nullptr)
 // node × 軸対称: 軸ノードで roUy 行 (index 2) のみ単位行化 (nullptr で無効)。
 const geom_int* __restrict__ axis_ur_flag,

 // axisymMethod==1 (isAxisymmetric enc==2) の軸ソース Jacobian ガード: 軸上ノード (==1) はソース 0 なので
 // Jacobian も加えない。decouple 用 axis_flag (nodeAxisDirichlet ゲート) とは独立に渡す (nullptr 可)。
 const geom_int* __restrict__ axis_flag_src,

 // node-centered 壁 no-slip: 壁ノードで運動量3行 (index1=roUx,2=roUy,3=roUz) を decouple する (nullptr 可)。
 // SU2 `DeleteValsRowi` 相当。残差射影だけでは block-DPLUR が壁運動量を連成したまま dq≠0 を返し速度 drift
 // するのを防ぐ。連続(行0)・エネルギー(行4)は保持。methods/discretization.md §7.2.1。
 const geom_int* __restrict__ wall_flag,

 // node-centered 等温壁: 壁ノードでエネルギー行 (index4=roe) を decouple する (nullptr 可)。
 // 壁ノード T ピン (applyNodeIsothermalWallPin / WMLES 等温 pin) と対。ピンで状態を上書きしながら
 // エネルギー行を連成したまま解くと Jacobian 不整合で発散する (2026-07-20 純伝導検証で実測)。
 const geom_int* __restrict__ iso_wall_flag,
 // 等温壁エネルギー境界の弱形式 (mesh.nodeIsothermalEnergyBC=1) の近似対角項の素材
 // g = Σ_{壁半割面} k_eff A_half / d_1 [W/K] (viscousFlux_wall_d が残差と同じ k_eff・幾何で積む)。
 // 厳密微分は対角 0 で温度微分は内部点の列にある。ここで足すのは SU2 型の**近似対角**である
 // (codex plan レビュー M1)。forge は残差微分の符号を反転して行列を組むので **+** で足す。
 // CPG: ΔA[4][4] = + g/(ρ c_v)。nullptr なら何もしない (既定はビット不変)。
 const flow_float* __restrict__ weakIsoDiag,
 flow_float weakIsoCv,

 // node-centered 弱形式 (Phase 2, 5e): node モードはゴーストセルを使わない。境界半割面 (has_nbr=false=ゴースト
 // 側) をこの node-to-node Jacobian ループから完全に除外する (continue)。境界ノードは物理境界上に乗るため
 // node→ghost が退化 (dcc≈0) し粘性対角 2ν·delta/dcc が爆発→対角巨大→dq≈0 で境界ノードが凍結する (出口 BL
 // 崩壊・残差プラトーの真因)。境界の対流/粘性は弱形式カーネルが残差側で担う。cell (isNode=0) は境界ゴースト
 // が法線方向に正しく置かれ非退化なので従来どおり境界面も処理する。
 int isNode,

 // --- line-implicit (plans/active/time_integration-line-implicit.md) ---
 // line_prev/next != nullptr で有効。ライン CV は点解せず、diag (storeLU 時)・rhs (毎 sweep)・
 // 近傍行列 K (storeLU 時) を保存して Thomas カーネルに委ねる。
 const geom_int* line_prev,
 const geom_int* line_next,
 flow_float* Kprev,
 flow_float* Knext,
 int storeLU,
 int lineViscCoupling,
 // 対角キャッシュ (plans/active/performance-3d-node-sst-speedup.md §4.2-4): 1 のとき loop==0 で組んだ対角 5×5
 // (拘束行の単位行化込み) を diag_** に保存し、loop>0 は対角組立・粘性対角・軸対称 Jacobian・近傍幾何読みを省略して
 // 保存値を読む。状態は sweep 中凍結なので結果はビット同一。線形 solve・rhs 拘束・近傍積は毎 sweep 従来どおり。
 // 呼び出し側で float・point 経路 (line_prev==nullptr) に限定する。
 int useDiagCache,
 // 近傍 dq の AoS 版 (stride 8 floats = 32 B セクタ整列, [0..4] を使用)。非 nullptr のとき近傍 gather は
 // dq_pack_old から float4+float の 2 ロード (SoA 5 配列の 5 ロード = 5 セクタから 1 セクタへ)。dq_pack_new には
 // dq_new と同じ値を書く。SoA の dq_new_* も従来どおり書く (commit・周期ミラー・診断が読む)。
 // line-implicit / node 周期 (SoA だけを書き換える経路) では呼び出し側が nullptr を渡す。
 const flow_float* __restrict__ dq_pack_old,
 flow_float* dq_pack_new,
 // エネルギー行の熱伝導 Jacobian (plans/active/time_integration-implicit-thermal-jacobian.md、ビットマスク、0 で従来どおり):
 //   ビット 1: 内部の node 間面でエネルギー行の粘性対角を k_face·δ/dcc·(γ/cp)·∂e/∂Q に置き換える (k_face は残差と同じ式)。
 //   ビット 2: 等温壁の節点のエネルギー行を拘束の行 [−e_w,0,0,0,1] にする。
 //   ビット 4 (ビット 1 と併用): 行 4 の従来のスカラー 2ν_eff·δ/dcc も残し、温度の項はその上に足す。
 // thermCondArr・cpArr・fxArr は thermalJac のビット 1 が立っているときだけ読む (それ以外は nullptr でよい)。
 int thermalJac,
 const flow_float* __restrict__ thermCondArr,
 const flow_float* __restrict__ cpArr,
 const flow_float* __restrict__ fxArr,
 flow_float Prt,
 const flow_float* __restrict__ visLamArr,   // 節点ごとの層流粘性 (lineViscCoupling 2 のときだけ非 nullptr)
 flow_float* dbgLineVisc,                    // デバッグ (FORGE_LINE_DUMP_DIR): ライン面のスカラー 2ν_eff·δ/dcc の和を storeLU の sweep で書く。通常 nullptr
 int lineViscTerms                           // 診断 (FORGE_LVC_TERMS、既定 7): 値 2/3 の薄層の項のマスク (plan time_integration-line-viscous-jacobian §6.8)
)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        // 状態/残差 (flow_float=float) を solve 精度 ST へキャストして取り込む (混合精度)。
        const ST gamma = static_cast<ST>(gamma_arr[ic]);   // 局所 γ
        const ST v = static_cast<ST>(vol[ic]);
        const ST dt_l = static_cast<ST>(dt_local[ic]);
        const ST density = max(static_cast<ST>(ro[ic]), static_cast<ST>(1.0e-30));
        const ST velocity_x = static_cast<ST>(Ux[ic]);
        const ST velocity_y = static_cast<ST>(Uy[ic]);
        const ST velocity_z = static_cast<ST>(Uz[ic]);
        const ST local_sonic = max(static_cast<ST>(sonic[ic]), static_cast<ST>(1.0e-8));
        const ST local_Ht = static_cast<ST>(Ht[ic]);   // 一般EOS固有系のエネルギー成分 (TP)
        // line-implicit: 自 CV の役割 (ラインに載るか) と decouple 行マスク (保存 K の行ゼロ化用)
        const geom_int lp = (line_prev != nullptr) ? line_prev[ic] : (geom_int)(-1);
        const geom_int ln_ = (line_next != nullptr) ? line_next[ic] : (geom_int)(-1);
        const bool onLine = (lp >= 0) || (ln_ >= 0);
        // 対角キャッシュを読む sweep か (loop>0 かつ line に載らない CV)。
        const bool cached = (useDiagCache != 0) && (loop > 0) && !onLine;
        // ライン上の節点の対角は storeLU の sweep (Thomas の因子を作る sweep) でしか使わない。それ以外の sweep では
        // 組んでも捨てるだけなので組まない (plan time_integration-line-implicit-speed §4.3 案 A、数値は不変)。
        const bool skipDiag = cached || (onLine && storeLU == 0);
        ST nu_eff = static_cast<ST>(0.0);
        if (!skipDiag) nu_eff = (static_cast<ST>(laminar_visc) + max(static_cast<ST>(vis_turb[ic]), static_cast<ST>(0.0))) / density;
        bool rowDec[5] = {false, false, false, false, false};
        if (onLine) {
            if (axis_ur_flag != nullptr && axis_ur_flag[ic] == 1) rowDec[2] = true;
            if (wall_flag != nullptr && wall_flag[ic] == 1) { rowDec[1] = true; rowDec[2] = true; rowDec[3] = true; }
            if (iso_wall_flag != nullptr && iso_wall_flag[ic] == 1) rowDec[4] = true;
        }

        // (loop==0 の dq_old ゼロ化は blockDPLURSolve の cudaMemset が担う。dq_old は本カーネルでは読み取り専用
        //  (const __restrict__) にして read-only キャッシュ経路を許す。dq_new とは別バッファ = 別名無し。)

        ST diag_block[5][5];
        block_dplur::zero5x5(diag_block);
        if (!skipDiag) {
            block_dplur::add_identity_scaled(diag_block, static_cast<ST>(v / max(dt_l, static_cast<ST>(1.0e-30))));
            // dual-time: 物理時間項 a·V/Δt を対角へ（定常は unsteady_diag==0）。
            block_dplur::add_identity_scaled(diag_block, v * static_cast<ST>(unsteady_diag));
        }

        ST rhs[5] = {
            static_cast<ST>(res_ro[ic]),
            static_cast<ST>(res_roUx[ic]),
            static_cast<ST>(res_roUy[ic]),
            static_cast<ST>(res_roUz[ic]),
            static_cast<ST>(res_roe[ic])
        };
        ST neighbor_accum[5];
        block_dplur::zero5(neighbor_accum);
        ST dbgLineViscSum = static_cast<ST>(0.0);   // デバッグ: ライン面のスカラー粘性の対角の和 (値 2 で外した量)

        const geom_int plane_begin = cell_planes_index[ic];
        const geom_int plane_end = cell_planes_index[ic + 1];
        for (geom_int plane_offset = plane_begin; plane_offset < plane_end; ++plane_offset) {
            const geom_int ip = cell_planes[plane_offset];
            const ST face_area = max(static_cast<ST>(ss[ip]), static_cast<ST>(1.0e-30));
            const geom_int ic0 = plane_cells[2 * ip + 0];
            const geom_int ic1 = plane_cells[2 * ip + 1];
            const geom_int other_ic = (ic0 == ic) ? ic1 : ic0;

            // 格納法線 (sx,sy,sz) は ic0→ic1。ic が neighbor 側 (ic1) のとき符号反転。
            const ST nsign = (ic0 == ic) ? static_cast<ST>(1.0) : static_cast<ST>(-1.0);
            const ST nx = nsign * static_cast<ST>(sx[ip]) / face_area;
            const ST ny = nsign * static_cast<ST>(sy[ip]) / face_area;
            const ST nz = nsign * static_cast<ST>(sz[ip]) / face_area;

            // 閉形式 FVS (R/L を作らず rank-2 外積で A⁺S を diag・k_off·sdq を nbr へ)。
            // 対流項: has_nbr=false (境界半割面) のとき自セル状態+面法線から A⁺S を対角に積むだけで、ゴースト
            // セルの状態/中心は一切読まない (ghostless)。よって node モードでもこの対流寄与は残す (境界ノードの
            // 流出 Jacobian=陰的安定化に必要。除くと rms_roUy 等が発散した)。
            const bool has_nbr = (other_ic < nCells);
            // ライン面: dq_old の lag 参照をスキップ (Thomas が厳密連成) — sdq=0 で対角 A⁺ だけ積む
            const bool isLineFace = onLine && has_nbr && (other_ic == lp || other_ic == ln_);
            ST sdq[5] = {static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0)};
            // loop==0 は dq_old≡0 (blockDPLURSolve の memset) なので gather を省く (寄与は厳密に 0 = ビット同一)。
            if (has_nbr && !isLineFace && loop > 0) {
                if (dq_pack_old != nullptr) {
                    const float4 q4 = *reinterpret_cast<const float4*>(dq_pack_old + (size_t)other_ic * 8);
                    const flow_float q5 = dq_pack_old[(size_t)other_ic * 8 + 4];
                    sdq[0] = face_area * static_cast<ST>(q4.x);
                    sdq[1] = face_area * static_cast<ST>(q4.y);
                    sdq[2] = face_area * static_cast<ST>(q4.z);
                    sdq[3] = face_area * static_cast<ST>(q4.w);
                    sdq[4] = face_area * static_cast<ST>(q5);
                } else {
                    sdq[0] = face_area * static_cast<ST>(dq_old_0[other_ic]);
                    sdq[1] = face_area * static_cast<ST>(dq_old_1[other_ic]);
                    sdq[2] = face_area * static_cast<ST>(dq_old_2[other_ic]);
                    sdq[3] = face_area * static_cast<ST>(dq_old_3[other_ic]);
                    sdq[4] = face_area * static_cast<ST>(dq_old_4[other_ic]);
                }
            }
            if (skipDiag) {
                block_dplur::accumulate_split_jacobian_cf<ST, false>(
                    gamma, nx, ny, nz, velocity_x, velocity_y, velocity_z,
                    local_sonic, local_Ht, thermallyPerfect != 0,
                    face_area, has_nbr, sdq, diag_block, neighbor_accum
                );
            } else {
                block_dplur::accumulate_split_jacobian_cf<ST>(
                    gamma, nx, ny, nz, velocity_x, velocity_y, velocity_z,
                    local_sonic, local_Ht, thermallyPerfect != 0,
                    face_area, has_nbr, sdq, diag_block, neighbor_accum
                );
            }
            // K 行列の列抽出 (状態凍結ゆえ storeLU=loop0 のみ): nbr 寄与は sdq に線形なので
            // 単位ベクトル×face_area で列が得られる (対角へは dummy に捨てる)。decouple 行は 0。
            if (isLineFace && storeLU != 0) {
                flow_float* Kdst = (other_ic == lp) ? Kprev : Knext;
                ST ddum[5][5];
                ST kcol[5];
                ST evec[5];
                for (int j = 0; j < 5; ++j) {
                    block_dplur::zero5x5(ddum);
                    block_dplur::zero5(kcol);
                    #pragma unroll
                    for (int q = 0; q < 5; ++q) evec[q] = static_cast<ST>(0.0);
                    evec[j] = face_area;
                    block_dplur::accumulate_split_jacobian_cf<ST>(
                        gamma, nx, ny, nz, velocity_x, velocity_y, velocity_z,
                        local_sonic, local_Ht, thermallyPerfect != 0,
                        face_area, true, evec, ddum, kcol);
                    #pragma unroll
                    for (int i = 0; i < 5; ++i)
                        Kdst[(size_t)ic * 25 + i * 5 + j] =
                            rowDec[i] ? (flow_float)0.0 : static_cast<flow_float>(kcol[i]);
                }
            }

            // 粘性対角: node モードはゴーストセルを使わない。境界半割面 (has_nbr=false=ゴースト側) では
            // dcc 計算 (ccx[ghost] 読み) も viscous_diag も行わない。境界ノードは境界面上に乗るため node→ghost が
            // 退化 (dcc≈0) し 2ν·delta/dcc が爆発→対角巨大→dq≈0 で境界ノードが凍結する (出口 BL 崩壊・残差
            // プラトーの真因)。境界粘性は弱形式カーネルが残差側で担う。内部 node-to-node 面のみ粘性対角を課す。
            // cell モード (isNode=0) は境界ゴーストが法線方向に正しく置かれ非退化なので従来どおり境界面も課す。
            if (!skipDiag && !(isNode != 0 && !has_nbr)) {
                // cc[other] − cc[ic]: ic が ic0 側なら e、ic1 側なら −e を ST にする (§4.2a、段 ③)。旧は ST にしてから
                // 座標の差を取っていたので、ST = float (implicitSolvePrecision 0) では FP64 のビルドでも値が変わる (左辺だけ)。
                const ST dcc_x = static_cast<ST>((ic0 == ic) ? ge_x[ip] : -ge_x[ip]);
                const ST dcc_y = static_cast<ST>((ic0 == ic) ? ge_y[ip] : -ge_y[ip]);
                const ST dcc_z = static_cast<ST>((ic0 == ic) ? ge_z[ip] : -ge_z[ip]);
                const ST dcc = max(sqrt(dcc_x * dcc_x + dcc_y * dcc_y + dcc_z * dcc_z), static_cast<ST>(1.0e-30));
                const ST dcc_dot_s = max(
                    fabs(dcc_x * static_cast<ST>(sx[ip]) + dcc_y * static_cast<ST>(sy[ip]) + dcc_z * static_cast<ST>(sz[ip])),
                    static_cast<ST>(1.0e-30)
                );
                const ST delta = max(dcc * face_area * face_area / dcc_dot_s, static_cast<ST>(1.0e-30));
                // 粘性対角は residual の粘性流束 Jacobian (2ν·ss²/dcc_dot_s = 2ν·delta/dcc) と整合させる
                // (旧 face_area·(2ν/delta) は ≈2ν に潰れ近軸で r 重み喪失・ゼロ面積面にスプリアス。詳細は site1 コメント)。
                const ST viscous_diag = static_cast<ST>(2.0) * nu_eff * delta / dcc;
                if (isLineFace) dbgLineViscSum += viscous_diag;
                if (isLineFace && lineViscCoupling >= 2) {
                    // 薄層の粘性・熱伝導の Jacobian (plan time_integration-line-viscous-jacobian §4.1)。ライン面では従来の
                    // スカラー 2ν·δ/dcc の代わりに、残差と同じ面の μ_f・k_f で D (自節点) と K (ライン上の隣) を組む。
                    // 値 3 (診断、§6.4): 従来のスカラーも全行の対角に残す (拘束の行は後で単位行に上書きされる)。
                    if (lineViscCoupling == 3) block_dplur::add_identity_scaled(diag_block, viscous_diag);
                    const ST f0 = static_cast<ST>(fxArr[ip]);
                    const ST fi = (ic0 == ic) ? f0 : static_cast<ST>(1.0) - f0;    // 自節点の補間の重み
                    const ST omfi = static_cast<ST>(1.0) - fi;
                    const ST mlam_i = (visLamArr != nullptr) ? static_cast<ST>(visLamArr[ic]) : static_cast<ST>(laminar_visc);
                    const ST mlam_j = (visLamArr != nullptr) ? static_cast<ST>(visLamArr[other_ic]) : static_cast<ST>(laminar_visc);
                    const ST mut_i = static_cast<ST>(vis_turb[ic]);
                    const ST mut_j = static_cast<ST>(vis_turb[other_ic]);
                    const ST cp_i = static_cast<ST>(cpArr[ic]);
                    const ST cp_j = static_cast<ST>(cpArr[other_ic]);
                    const ST mu_f = fi * (mlam_i + mut_i) + omfi * (mlam_j + mut_j);
                    const ST k_f = fi * static_cast<ST>(thermCondArr[ic]) + omfi * static_cast<ST>(thermCondArr[other_ic])
                                 + (fi * cp_i + omfi * cp_j) * (fi * mut_i + omfi * mut_j) / static_cast<ST>(Prt);
                    const ST beta = max(mu_f, static_cast<ST>(0.0)) * delta / dcc;
                    const ST kappa = max(k_f, static_cast<ST>(0.0)) * delta / dcc;
                    const bool jVel = (wall_flag != nullptr && wall_flag[other_ic] == 1);
                    const bool jTemp = (iso_wall_flag != nullptr && iso_wall_flag[other_ic] == 1);
                    ST Kv[5][5];
                    block_dplur::zero5x5(Kv);
                    block_dplur::accumulate_thinlayer_visc_jacobian<ST>(
                        beta, kappa, nx, ny, nz, fi,
                        density, velocity_x, velocity_y, velocity_z, static_cast<ST>(roe[ic]), gamma, max(cp_i, static_cast<ST>(1.0e-30)),
                        max(static_cast<ST>(ro[other_ic]), static_cast<ST>(1.0e-30)),
                        static_cast<ST>(Ux[other_ic]), static_cast<ST>(Uy[other_ic]), static_cast<ST>(Uz[other_ic]),
                        static_cast<ST>(roe[other_ic]), static_cast<ST>(gamma_arr[other_ic]), max(cp_j, static_cast<ST>(1.0e-30)),
                        jVel, jTemp, diag_block, (storeLU != 0) ? Kv : nullptr, lineViscTerms);
                    if (storeLU != 0) {
                        flow_float* Kdst = (other_ic == lp) ? Kprev : Knext;
                        for (int i = 0; i < 5; ++i)
                            if (!rowDec[i])
                                for (int j = 0; j < 5; ++j) Kdst[(size_t)ic * 25 + i * 5 + j] += static_cast<flow_float>(Kv[i][j]);
                    }
                } else if (isLineFace && lineViscCoupling == 1) {
                    // v2 (plans/active/time_integration-line-implicit-viscous-v2.md): line 面は
                    // スカラー粘性結合 K += α·I (α=ν_eff·δ/dcc) と対にし、対角は 2α→α に置換して
                    // 真の 1D 拡散行 [−α, 2α, −α] を line 内で完成させる (off-line 面は従来 2α のまま)。
                    const ST alpha = static_cast<ST>(0.5) * viscous_diag;
                    block_dplur::add_identity_scaled(diag_block, alpha);
                    if (storeLU != 0) {
                        flow_float* Kdst = (other_ic == lp) ? Kprev : Knext;
                        #pragma unroll
                        for (int i = 0; i < 5; ++i)
                            if (!rowDec[i]) Kdst[(size_t)ic * 25 + i * 5 + i] += static_cast<flow_float>(alpha);
                    }
                } else if ((thermalJac & 1) != 0 && isNode != 0 && has_nbr) {
                    // 熱伝導の Jacobian (implicitThermalJacobian ビット 1): 連続・運動量の行は従来どおりスカラー、
                    // エネルギー行は熱伝導の残差 k_face·(T_j − T_i)·δ/dcc の Q_i による微分の符号反転 = Λ^T·(γ/cp)·∂e/∂Q。
                    // k_face は viscousFlux_d.cu の tc_face と同じ式 (f 補間の層流 k + 面 cp × 面 μ_t / Pr_t)、物性・γ は凍結。
                    // ∂e/∂ρ = −(e − ½|u|²)/ρ、∂e/∂(ρu_k) = −u_k/ρ、∂e/∂(ρE) = 1/ρ (e = ρE/ρ − ½|u|²; TP でも Y は正規化済みなので厳密)。
                    #pragma unroll
                    // ビット 4: 行 4 にも従来のスカラー (スペクトル半径の近似) を残し、その上に温度の項を足す。
                    // 行 0〜3 と同じ一律の減衰を行 4 にも保つので、行ごとに歩幅の縮め方が違う状態 (温度を変えない
                    // 補正でエネルギー行だけ減衰が抜け、高速域で偽の ΔT を作る) を避ける (plan §4.3)。
                    const int nScalarRows = ((thermalJac & 4) != 0) ? 5 : 4;
                    for (int r = 0; r < nScalarRows; ++r) diag_block[r][r] += viscous_diag;
                    const ST f = static_cast<ST>(fxArr[ip]);
                    const ST omf = static_cast<ST>(1.0) - f;
                    const ST cp_face = f * static_cast<ST>(cpArr[ic0]) + omf * static_cast<ST>(cpArr[ic1]);
                    const ST mut_face = f * static_cast<ST>(vis_turb[ic0]) + omf * static_cast<ST>(vis_turb[ic1]);
                    const ST k_face = f * static_cast<ST>(thermCondArr[ic0]) + omf * static_cast<ST>(thermCondArr[ic1])
                                    + cp_face * mut_face / static_cast<ST>(Prt);
                    const ST q2 = velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z;
                    const ST e_int = static_cast<ST>(roe[ic]) / density - static_cast<ST>(0.5) * q2;
                    const ST cfac = (k_face * delta / dcc) * (gamma / max(static_cast<ST>(cpArr[ic]), static_cast<ST>(1.0e-30))) / density;
                    diag_block[4][0] += -cfac * (e_int - static_cast<ST>(0.5) * q2);
                    diag_block[4][1] += -cfac * velocity_x;
                    diag_block[4][2] += -cfac * velocity_y;
                    diag_block[4][3] += -cfac * velocity_z;
                    diag_block[4][4] += cfac;
                } else {
                    block_dplur::add_identity_scaled(diag_block, viscous_diag);
                    // **診断専用 A/B** (codex 2026-09-21、既定はコンパイルから除外されビット不変)。
                    // 粘性対角 2ν_eff·delta/dcc は `add_identity_scaled` で**5 行すべてに同じ量**が入る。
                    // エネルギー行の真の拡散 Jacobian はこれと違う (完全気体で rho,rho u を固定すると
                    // ∂T/∂(ρE)=1/(ρc_v) なので k∂T/∂(ρE)=γα、さらに交差微分がある) が、
                    // **不足率は α/ν=1.11-1.39 からは証明できない** (codex 指摘)。
                    // ここで調べるのは「非収束状態の 2 節点指標がエネルギー対角に感度を持つか」だけで、
                    // 1.4 倍は物理係数の再現ではなく試験強度である。
                    // 効いた場合も「Pr 補正が正しい」ではなく「指標が陰解法作用素に依存する」と結論する。
                    // plan boundary-conjugate-heat-transfer §5.1 #43。
#if defined(FORGE_TEST_ENERGY_VISCOUS_DIAG)
                    diag_block[4][4] += static_cast<ST>(0.4) * viscous_diag;
#endif
                }
            }
        }

        #pragma unroll
        for (int i = 0; i < 5; ++i) {
            rhs[i] += neighbor_accum[i];
        }
        if (dbgLineVisc != nullptr && storeLU != 0) dbgLineVisc[ic] = static_cast<flow_float>(dbgLineViscSum);

        // 軸対称ソース項のヤコビアンを対角ブロックに加える（roUy 行 = index 2）。詳細は実装ドキュメント参照。
        // axisRFloor 帯 (r 床, ソース不課) は Jacobian も課さない。
        if (!skipDiag && isAxisymmetric == 1 &&
            !(static_cast<ST>(axisRFloor) > static_cast<ST>(0.0) && static_cast<ST>(ccy[ic]) < static_cast<ST>(axisRFloor))) {
            const ST A_pl = static_cast<ST>(A_planar[ic]);
            const ST r_eff = max(v / max(A_pl, static_cast<ST>(1.0e-30)), static_cast<ST>(1.0e-30));
            const ST g1 = gamma - static_cast<ST>(1.0);
            const ST q2 = velocity_x*velocity_x + velocity_y*velocity_y + velocity_z*velocity_z;
            const ST mu_total = static_cast<ST>(laminar_visc) + max(static_cast<ST>(vis_turb[ic]), static_cast<ST>(0.0));
            const ST hoop = static_cast<ST>(2.0) * mu_total / (density * r_eff);
            diag_block[2][0] += -A_pl * (static_cast<ST>(0.5)*g1*q2 + hoop * velocity_y);
            diag_block[2][1] += A_pl * (g1 * velocity_x);
            diag_block[2][2] += A_pl * (g1 * velocity_y + hoop);
            diag_block[2][3] += A_pl * (g1 * velocity_z);
            diag_block[2][4] += -A_pl * g1;
            // 診断: 近軸半径音響スペクトル半径 α·A_pl·c を roUy 対角に補う (FORGE_AXIS_DIAG_ALPHA>0 のみ)。
            diag_block[2][2] += static_cast<ST>(g_axisDiagAlpha) * A_pl * local_sonic;
        } else if (!skipDiag && isAxisymmetric == 2) {
            // SU2 流 (axisymMethod==1) 非粘性軸対称ソースの解析 Jacobian (CSourceAxisymmetric_Flow 移植,
            // 行/列 = [ro, roUx, roUy, roe] → forge [0,1,2,4])。forge 対角は -∂S/∂U = +SU2 jacobian。
            // 軸ノード (axis_flag_src==1) と y≤eps はソース 0 のためスキップ。γ は frozen (gamma_arr)。
            const ST y = static_cast<ST>(ccy[ic]);
            const bool onAxisSrc = (axis_flag_src != nullptr && axis_flag_src[ic] == 1);
            if (!onAxisSrc && y > static_cast<ST>(1.0e-12)) {
                const ST yv = v / y;
                const ST g1 = gamma - static_cast<ST>(1.0);
                const ST uu = velocity_x, ww = velocity_y;
                const ST q2d = uu*uu + ww*ww;
                const ST et = static_cast<ST>(roe[ic]) / density;   // 比全エネルギー
                diag_block[0][2] += yv;
                diag_block[1][0] += yv * (-uu * ww);
                diag_block[1][1] += yv * ww;
                diag_block[1][2] += yv * uu;
                diag_block[2][0] += yv * (-ww * ww);
                diag_block[2][2] += yv * static_cast<ST>(2.0) * ww;
                diag_block[4][0] += yv * (-gamma * ww * et + g1 * ww * q2d);
                diag_block[4][1] += yv * (-g1 * uu * ww);
                diag_block[4][2] += yv * (gamma * et - static_cast<ST>(0.5) * g1 * (q2d + static_cast<ST>(2.0) * ww * ww));
                diag_block[4][4] += yv * (gamma * ww);
                // 粘性軸対称ソースの stiff 主対角: S_roUy ∋ -V·2μ_tot·v/y² → -∂S/∂(ρv) = +V·2μ/(ρy²)。
                // 近軸第一列 (y~1e-4) で極めて stiff で、これを lag すると implicit が喉部近軸で
                // limit cycle 化し rms_ro ~1e-5 で頭打ちになる (explicit は 3e-7 到達 = 空間は健全)。
                const ST mu_tot_ax = static_cast<ST>(laminar_visc) + max(static_cast<ST>(vis_turb[ic]), static_cast<ST>(0.0));
                diag_block[2][2] += yv * static_cast<ST>(2.0) * mu_tot_ax / (density * y);
            }
        }

        // node × 軸対称: 軸ノードの半径運動量行のみ decouple (dq_roUy=0)。状態は enforceAxisSymmetry がピン。
        if (axis_ur_flag != nullptr && axis_ur_flag[ic] == 1) {
            if (!skipDiag) {
                for (int jj = 0; jj < 5; ++jj) diag_block[2][jj] = static_cast<ST>(0.0);
                diag_block[2][2] = static_cast<ST>(1.0);
            }
            rhs[2] = static_cast<ST>(0.0);
        }

        // SU2 `DeleteValsRowi` 相当の壁 no-slip Dirichlet: 壁ノードで運動量3行 (index 1,2,3) を単位行に
        // 置換し rhs=0 → solve が一貫して dq_roUx=dq_roUy=dq_roUz=0 を返す。連続(0)・エネルギー(4)行は
        // 保持され ρ,ρe は保存式で発展、圧力は EOS が復元 (CPG/TP 共通)。残差射影だけでは block-DPLUR が
        // 壁運動量を連成し dq≠0 を返して壁速度が drift する問題を Jacobian 整合で根治する。
        if (wall_flag != nullptr && wall_flag[ic] == 1) {
            for (int row = 1; row <= 3; ++row) {
                if (!skipDiag) {
                    for (int jj = 0; jj < 5; ++jj) diag_block[row][jj] = static_cast<ST>(0.0);
                    diag_block[row][row] = static_cast<ST>(1.0);
                }
                rhs[row] = static_cast<ST>(0.0);
            }
        }

        // 弱形式の等温壁 (nodeIsothermalEnergyBC=1): エネルギー行は残したまま、壁寄与の近似対角を足す。
        // iso_wall_flag が nullptr になっているので下の単位行化とは排他。
        if (weakIsoDiag != nullptr && !skipDiag) {
            const ST g = static_cast<ST>(weakIsoDiag[ic]);
            if (g > static_cast<ST>(0.0)) {
                const ST rho = static_cast<ST>(max(ro[ic], (flow_float)1.0e-30));
                diag_block[4][4] += g / (rho * static_cast<ST>(weakIsoCv));
            }
        }

        // 等温壁ノード: エネルギー行 (4) も単位行に置換し dq_roe=0 → 壁ノード T は pin (applyBconds 位相) が
        // 一意に決める。連続 (0) 行は保持 (ρ は保存式で発展し P=ρRTw が追従)。
        if (iso_wall_flag != nullptr && iso_wall_flag[ic] == 1) {
            if (!skipDiag) {
                for (int jj = 0; jj < 5; ++jj) diag_block[4][jj] = static_cast<ST>(0.0);
                diag_block[4][4] = static_cast<ST>(1.0);
                // implicitThermalJacobian ビット 2: 拘束の行 Δ(ρE)_w − e_w·Δρ_w = 0 (壁温のピン ρE = ρ·e(T_w) と一致、壁は u = 0)。
                // lineViscCoupling 2 も同じ拘束の行にする (隣の熱伝導の K を消す前提の ΔT_w = 0、plan time_integration-line-viscous-jacobian §4.1)。
                if ((thermalJac & 2) != 0 || lineViscCoupling >= 2) diag_block[4][0] = -static_cast<ST>(roe[ic]) / density;
            }
            rhs[4] = static_cast<ST>(0.0);
        }

        // 対角キャッシュ: loop>0 は保存値を読む / loop==0 (useDiagCache かつ line 外) は組んだ対角を保存する。
        // ST=float・point 経路に限定して呼ばれる (呼び出し側ゲート) ので、保存/読込で丸めは発生しない (ビット同一)。
        if (cached) {
            diag_block[0][0]=static_cast<ST>(diag_00[ic]); diag_block[0][1]=static_cast<ST>(diag_01[ic]); diag_block[0][2]=static_cast<ST>(diag_02[ic]); diag_block[0][3]=static_cast<ST>(diag_03[ic]); diag_block[0][4]=static_cast<ST>(diag_04[ic]);
            diag_block[1][0]=static_cast<ST>(diag_10[ic]); diag_block[1][1]=static_cast<ST>(diag_11[ic]); diag_block[1][2]=static_cast<ST>(diag_12[ic]); diag_block[1][3]=static_cast<ST>(diag_13[ic]); diag_block[1][4]=static_cast<ST>(diag_14[ic]);
            diag_block[2][0]=static_cast<ST>(diag_20[ic]); diag_block[2][1]=static_cast<ST>(diag_21[ic]); diag_block[2][2]=static_cast<ST>(diag_22[ic]); diag_block[2][3]=static_cast<ST>(diag_23[ic]); diag_block[2][4]=static_cast<ST>(diag_24[ic]);
            diag_block[3][0]=static_cast<ST>(diag_30[ic]); diag_block[3][1]=static_cast<ST>(diag_31[ic]); diag_block[3][2]=static_cast<ST>(diag_32[ic]); diag_block[3][3]=static_cast<ST>(diag_33[ic]); diag_block[3][4]=static_cast<ST>(diag_34[ic]);
            diag_block[4][0]=static_cast<ST>(diag_40[ic]); diag_block[4][1]=static_cast<ST>(diag_41[ic]); diag_block[4][2]=static_cast<ST>(diag_42[ic]); diag_block[4][3]=static_cast<ST>(diag_43[ic]); diag_block[4][4]=static_cast<ST>(diag_44[ic]);
        } else if (useDiagCache != 0 && !onLine) {
            diag_00[ic]=static_cast<flow_float>(diag_block[0][0]); diag_01[ic]=static_cast<flow_float>(diag_block[0][1]); diag_02[ic]=static_cast<flow_float>(diag_block[0][2]); diag_03[ic]=static_cast<flow_float>(diag_block[0][3]); diag_04[ic]=static_cast<flow_float>(diag_block[0][4]);
            diag_10[ic]=static_cast<flow_float>(diag_block[1][0]); diag_11[ic]=static_cast<flow_float>(diag_block[1][1]); diag_12[ic]=static_cast<flow_float>(diag_block[1][2]); diag_13[ic]=static_cast<flow_float>(diag_block[1][3]); diag_14[ic]=static_cast<flow_float>(diag_block[1][4]);
            diag_20[ic]=static_cast<flow_float>(diag_block[2][0]); diag_21[ic]=static_cast<flow_float>(diag_block[2][1]); diag_22[ic]=static_cast<flow_float>(diag_block[2][2]); diag_23[ic]=static_cast<flow_float>(diag_block[2][3]); diag_24[ic]=static_cast<flow_float>(diag_block[2][4]);
            diag_30[ic]=static_cast<flow_float>(diag_block[3][0]); diag_31[ic]=static_cast<flow_float>(diag_block[3][1]); diag_32[ic]=static_cast<flow_float>(diag_block[3][2]); diag_33[ic]=static_cast<flow_float>(diag_block[3][3]); diag_34[ic]=static_cast<flow_float>(diag_block[3][4]);
            diag_40[ic]=static_cast<flow_float>(diag_block[4][0]); diag_41[ic]=static_cast<flow_float>(diag_block[4][1]); diag_42[ic]=static_cast<flow_float>(diag_block[4][2]); diag_43[ic]=static_cast<flow_float>(diag_block[4][3]); diag_44[ic]=static_cast<flow_float>(diag_block[4][4]);
        }

        if (onLine) {
            // line-implicit: 点解せず Thomas カーネル用に保存する。
            //   diag: 状態凍結ゆえ storeLU (loop==0) のみ / rhs: ライン外 lag 込みなので毎 sweep。
            //   dq_new は Thomas が上書きする (保険で前回反復値を置く)。
            if (storeLU != 0) {
                diag_00[ic]=static_cast<flow_float>(diag_block[0][0]); diag_01[ic]=static_cast<flow_float>(diag_block[0][1]); diag_02[ic]=static_cast<flow_float>(diag_block[0][2]); diag_03[ic]=static_cast<flow_float>(diag_block[0][3]); diag_04[ic]=static_cast<flow_float>(diag_block[0][4]);
                diag_10[ic]=static_cast<flow_float>(diag_block[1][0]); diag_11[ic]=static_cast<flow_float>(diag_block[1][1]); diag_12[ic]=static_cast<flow_float>(diag_block[1][2]); diag_13[ic]=static_cast<flow_float>(diag_block[1][3]); diag_14[ic]=static_cast<flow_float>(diag_block[1][4]);
                diag_20[ic]=static_cast<flow_float>(diag_block[2][0]); diag_21[ic]=static_cast<flow_float>(diag_block[2][1]); diag_22[ic]=static_cast<flow_float>(diag_block[2][2]); diag_23[ic]=static_cast<flow_float>(diag_block[2][3]); diag_24[ic]=static_cast<flow_float>(diag_block[2][4]);
                diag_30[ic]=static_cast<flow_float>(diag_block[3][0]); diag_31[ic]=static_cast<flow_float>(diag_block[3][1]); diag_32[ic]=static_cast<flow_float>(diag_block[3][2]); diag_33[ic]=static_cast<flow_float>(diag_block[3][3]); diag_34[ic]=static_cast<flow_float>(diag_block[3][4]);
                diag_40[ic]=static_cast<flow_float>(diag_block[4][0]); diag_41[ic]=static_cast<flow_float>(diag_block[4][1]); diag_42[ic]=static_cast<flow_float>(diag_block[4][2]); diag_43[ic]=static_cast<flow_float>(diag_block[4][3]); diag_44[ic]=static_cast<flow_float>(diag_block[4][4]);
            }
            rhs_0[ic] = static_cast<flow_float>(rhs[0]);
            rhs_1[ic] = static_cast<flow_float>(rhs[1]);
            rhs_2[ic] = static_cast<flow_float>(rhs[2]);
            rhs_3[ic] = static_cast<flow_float>(rhs[3]);
            rhs_4[ic] = static_cast<flow_float>(rhs[4]);
            dq_new_0[ic] = dq_old_0[ic];
            dq_new_1[ic] = dq_old_1[ic];
            dq_new_2[ic] = dq_old_2[ic];
            dq_new_3[ic] = dq_old_3[ic];
            dq_new_4[ic] = dq_old_4[ic];
        } else {
        // diag_block を破壊して in-place で解く (solve_mat コピー排除)。
        ST correction[5] = {static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0), static_cast<ST>(0.0)};
        const bool ok = block_dplur::solve_5x5(diag_block, rhs, correction);
        if (!ok) {
            block_dplur::zero5(correction);
        }

        const ST relax = static_cast<ST>(implicit_relax);
        #pragma unroll
        for (int i = 0; i < 5; ++i) {
            correction[i] *= relax;
        }

        // 古典 DPLUR: float dq_new へ書戻し。Q への commit は applyBlockImplicitCorrection。
        dq_new_0[ic] = static_cast<flow_float>(correction[0]);
        dq_new_1[ic] = static_cast<flow_float>(correction[1]);
        dq_new_2[ic] = static_cast<flow_float>(correction[2]);
        dq_new_3[ic] = static_cast<flow_float>(correction[3]);
        dq_new_4[ic] = static_cast<flow_float>(correction[4]);
        if (dq_pack_new != nullptr) {
            *reinterpret_cast<float4*>(dq_pack_new + (size_t)ic * 8) =
                make_float4(static_cast<flow_float>(correction[0]), static_cast<flow_float>(correction[1]),
                            static_cast<flow_float>(correction[2]), static_cast<flow_float>(correction[3]));
            dq_pack_new[(size_t)ic * 8 + 4] = static_cast<flow_float>(correction[4]);
        }
        // rhs_** の診断書き出しは撤去 (読者なし。5 配列×sweep の書込 ≈240 MB/step を節約, 2026-09-12)。
        // line 経路 (上の onLine 分岐) は Thomas カーネルが rhs を読むので従来どおり書く。
        }
    }
}

// =============================================================================
// Phase 4 (a): 完全 Γ⁻¹A 低マッハ前処理の block DPLUR (lowMachPrecond>=2; 2=RHS+LHS / 3=LHS-only)。
// 既存 implicit_defect_correction_block_d (lowMachPrecond 0/1) とは別カーネルにし、
// 0/1 経路のレジスタ・ビットを一切変えない。
// 保存形は (Γ_c·V/Δτ' + A_c)ΔQ = -R で、**前処理は擬似時間項 Γ_c のみ**。フラックス A_c は
// 既存と同じ物理厳密 FVS (a_plus=A_c⁺・k_off=-A_c⁻) をそのまま使う (非前処理が正しい)。収束は
// (Δτ'/V)Γ_c⁻¹A_c の固有値 λ' で一様に前処理され、スカラー Δτ'=cell/ρ' で効く。
//   - 擬似時間項: Γ_c·V/Δτ' (Δτ' は setDT 側で前処理スペクトル半径から拡大した dt_local)。
//   - フラックス: 物理 a_plus を対角・k_off を近傍 (既存 block と同一)。
//   - 物理 BDF 項 a·V/Δt·I は非前処理 (dual-time 所有)。
// **Sherman-Morrison 解法**: Γ_c=I+α g rᵀ がランク1なので D=D0+γ g rᵀ (D0=物理ブロック・良条件)。
//   D0 を float で 2 RHS 同時 (solve_5x5_2rhs) に解き、悪条件 ~1/β は分母スカラーのみ double に隔離。
//   FP64 を回避して 0/1 カーネルに近い速度。β=1 (超音速) で Γ_c=I・Δτ'=Δτ・フラックス同一ゆえ現行と解一致。
// 理論: methods/time_integration/theory.md「低マッハ前処理固有系」、計画 §5 Phase 4。
__global__ void __launch_bounds__(BLOCK_DPLUR_THREADS) implicit_defect_correction_block_precond_d
(
 int loop,
 flow_float dt,
 flow_float* dt_local,
 flow_float implicit_relax,
 flow_float* gamma_arr,   // per-cell γ (TP: γ_mix(T), CPG: cfg.gamma)
 flow_float precondEps,

 geom_int nCells_all , geom_int nCells,
 geom_float* vol,
 geom_int* plane_cells,
 geom_int* cell_planes_index,
 geom_int* cell_planes,
 geom_float* ccx, geom_float* ccy, geom_float* ccz,
 // 面ごとの差 e = cc[ic1] − cc[ic0] (var.p_d["ge_*"]、plans/active/architecture-float-state-double-geometry.md §4.2a、段 ③)
 const flow_float* ge_x, const flow_float* ge_y, const flow_float* ge_z,
 geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,

 flow_float* ro, flow_float* roUx, flow_float* roUy, flow_float* roUz, flow_float* roe,

 flow_float laminar_visc,
 flow_float* vis_turb,
 flow_float* sonic,
 flow_float* Ux, flow_float* Uy, flow_float* Uz, flow_float* Ht,

 flow_float* res_ro, flow_float* res_roUx, flow_float* res_roUy, flow_float* res_roUz, flow_float* res_roe,

 flow_float* dq_old_0, flow_float* dq_old_1, flow_float* dq_old_2, flow_float* dq_old_3, flow_float* dq_old_4,
 flow_float* dq_new_0, flow_float* dq_new_1, flow_float* dq_new_2, flow_float* dq_new_3, flow_float* dq_new_4,

 int isAxisymmetric,
 flow_float* A_planar,
 flow_float axisRFloor,
 flow_float unsteady_diag
)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        const flow_float gamma = gamma_arr[ic];   // 局所 γ (TP の γ_mix。CPG は cfg.gamma で不変)
        const geom_float v = vol[ic];
        const flow_float dt_l = dt_local[ic];
        const flow_float density = max(ro[ic], static_cast<flow_float>(1.0e-30));
        const flow_float vx = Ux[ic];
        const flow_float vy = Uy[ic];
        const flow_float vz = Uz[ic];
        const flow_float local_sonic = max(sonic[ic], static_cast<flow_float>(1.0e-8));
        const flow_float local_enthalpy = max(Ht[ic], static_cast<flow_float>(1.0e-8));
        const flow_float nu_eff = (laminar_visc + max(vis_turb[ic], static_cast<flow_float>(0.0))) / density;
        const flow_float velMag = sqrt(vx*vx + vy*vy + vz*vz);
        const flow_float beta = lowMachBeta(local_sonic, velMag, precondEps);

        if (loop == 0) {
            dq_old_0[ic] = 0.0; dq_old_1[ic] = 0.0; dq_old_2[ic] = 0.0; dq_old_3[ic] = 0.0; dq_old_4[ic] = 0.0;
        }

        // 対角ブロックは Γ_c=I+α g rᵀ がランク1ゆえ D = D0 + γ g rᵀ と書ける:
        //   D0 = V/Δτ'·I + a·V/Δt·I + Σ A_c⁺ S + 粘性 + 軸対称  (物理ブロック・良条件・float 可)
        //   γ g rᵀ = (V/Δτ')·α g rᵀ                             (Γ_c 前処理寄与、悪条件 ~1/β の源)
        // Sherman-Morrison: x = y - [γ(rᵀy)/(1+γ(rᵀz))] z,  y=D0⁻¹b, z=D0⁻¹g。
        // D0 を float で 2 RHS 同時に解き、悪条件は分母スカラー(double)に隔離 → FP64 を回避 (RTX 等で高速)。
        flow_float D0[5][5];
        block_dplur::zero5x5(D0);
        const flow_float v_over_dtau = static_cast<flow_float>(v / max(dt_l, static_cast<flow_float>(1.0e-30)));
        block_dplur::add_identity_scaled(D0, v_over_dtau);
        block_dplur::add_identity_scaled(D0, static_cast<flow_float>(v) * unsteady_diag);  // dual-time BDF (非前処理)

        flow_float b[5] = { res_ro[ic], res_roUx[ic], res_roUy[ic], res_roUz[ic], res_roe[ic] };
        flow_float nbr[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

        const geom_int plane_begin = cell_planes_index[ic];
        const geom_int plane_end = cell_planes_index[ic + 1];
        for (geom_int plane_offset = plane_begin; plane_offset < plane_end; ++plane_offset) {
            const geom_int ip = cell_planes[plane_offset];
            const flow_float face_area = max(ss[ip], static_cast<flow_float>(1.0e-30));
            const geom_int ic0 = plane_cells[2 * ip + 0];
            const geom_int ic1 = plane_cells[2 * ip + 1];
            const geom_int other_ic = (ic0 == ic) ? ic1 : ic0;
            const flow_float nsign = (ic0 == ic) ? static_cast<flow_float>(1.0) : static_cast<flow_float>(-1.0);
            const flow_float nx = nsign * sx[ip] / face_area;
            const flow_float ny = nsign * sy[ip] / face_area;
            const flow_float nz = nsign * sz[ip] / face_area;

            // フラックスは物理の厳密 FVS。前処理は時間項のみ (保存形 (Γ_c V/Δτ'+A_c)ΔQ=-R)。
            flow_float a_plus[5][5];
            flow_float k_off[5][5];
            block_dplur::build_jacobian_split(gamma, nx, ny, nz, vx, vy, vz,
                                              local_enthalpy, local_sonic, a_plus, k_off);
            block_dplur::add_scaled_5x5(D0, a_plus, face_area);

            // cc[other] − cc[ic]: ic が ic0 側なら e、ic1 側なら −e (§4.2a、段 ③)
            const flow_float dcc_x = (ic0 == ic) ? ge_x[ip] : -ge_x[ip];
            const flow_float dcc_y = (ic0 == ic) ? ge_y[ip] : -ge_y[ip];
            const flow_float dcc_z = (ic0 == ic) ? ge_z[ip] : -ge_z[ip];
            const flow_float dcc = max(sqrt(dcc_x*dcc_x + dcc_y*dcc_y + dcc_z*dcc_z), static_cast<flow_float>(1.0e-30));
            const flow_float dcc_dot_s = max(fabs(dcc_x*sx[ip] + dcc_y*sy[ip] + dcc_z*sz[ip]), static_cast<flow_float>(1.0e-30));
            const flow_float delta = max(dcc * face_area * face_area / dcc_dot_s, static_cast<flow_float>(1.0e-30));
            // 粘性対角は residual の粘性流束 Jacobian (2ν·ss²/dcc_dot_s = 2ν·delta/dcc) と整合させる
            // (旧 face_area·(2ν/delta) は ≈2ν に潰れ近軸で r 重み喪失・ゼロ面積面にスプリアス。詳細は site1 コメント)。
            const flow_float viscous_diag = static_cast<flow_float>(2.0) * nu_eff * delta / dcc;
            block_dplur::add_identity_scaled(D0, viscous_diag);

            // 近傍 += k_off S ΔQ_nbr (= -A_c⁻ S ΔQ_nbr)。
            if (other_ic < nCells) {
                flow_float dqn[5];
                block_dplur::load_block_vec(other_ic, dq_old_0, dq_old_1, dq_old_2, dq_old_3, dq_old_4, dqn);
                #pragma unroll
                for (int i = 0; i < 5; ++i) dqn[i] *= face_area;
                block_dplur::multiply_add_5x5_vec(k_off, dqn, nbr);
            }
        }

        #pragma unroll
        for (int i = 0; i < 5; ++i) b[i] += nbr[i];

        // 軸対称ソースヤコビアン (物理ブロック D0 へ float で加算、既存 block と同式)。
        // axisRFloor 帯 (r 床, ソース不課) は Jacobian も課さない。
        if (isAxisymmetric == 1 &&
            !(axisRFloor > (flow_float)0.0 && ccy[ic] < axisRFloor)) {
            const flow_float A_pl = A_planar[ic];
            const flow_float r_eff = max(v / max(A_pl, static_cast<flow_float>(1.0e-30)), static_cast<flow_float>(1.0e-30));
            const flow_float g1 = gamma - static_cast<flow_float>(1.0);
            const flow_float q2 = vx*vx + vy*vy + vz*vz;
            const flow_float mu_total = laminar_visc + max(vis_turb[ic], static_cast<flow_float>(0.0));
            const flow_float hoop = static_cast<flow_float>(2.0) * mu_total / (density * r_eff);
            // ∂P/∂Q の第 1 成分は一般 EOS で χ_eos + κ e_k (χ_eos = c² − κ h)。CPG では χ_eos=0 で
            // 従来式 ½κq² にビット一致。TP (thermalMethod 2) では χ_eos≠0 で、これを落とすと
            // 対流 Jacobian (rvec, 下記) と軸ソース Jacobian が不整合になり、ホップ項が支配する
            // 軸近傍で block-DPLUR が発散する (case/42 run_0020: 一定 cp 種・陽解法では完走、
            // 実 NASA-9 + 陰解法のみ発散 → 2026-08-17 に特定)。
            const flow_float chi_hoop = local_sonic*local_sonic
                                      - g1*(local_enthalpy - static_cast<flow_float>(0.5)*q2);
            D0[2][0] += -A_pl * (chi_hoop + static_cast<flow_float>(0.5)*g1*q2 + hoop*vy);
            D0[2][1] +=  A_pl * (g1*vx);
            D0[2][2] +=  A_pl * (g1*vy + hoop);
            D0[2][3] +=  A_pl * (g1*vz);
            D0[2][4] += -A_pl * g1;
        }

        // Γ_c のランク1寄与: g=(1,u,v,w,H_t), r=∂p/∂Q=(χ_eos+κek,-κu,-κv,-κw,κ), γ=(V/Δτ')·(1-β)/(βc²)。
        // CPG では H_t=c²/(γ-1)+ek・χ_eos=0 で従来式に簡約 (ビット不変)。TP は実 H_t(=local_enthalpy) と
        // χ_eos=c²−κh を使う (build_jacobian_split と同じ一般EOS整合)。κ=γ-1。
        const flow_float ek = static_cast<flow_float>(0.5) * velMag * velMag;
        const flow_float gm1 = gamma - static_cast<flow_float>(1.0);
        const flow_float Htot = local_enthalpy;   // 実 Ht[ic] (CPG/TP 統一。CPG も Ht[ic]=ek+c²/(γ-1))
        const flow_float chi_eos = local_sonic*local_sonic - gm1*(local_enthalpy - ek);  // c²−κh (CPG で ≈0)
        flow_float gvec[5] = { static_cast<flow_float>(1.0), vx, vy, vz, Htot };
        const flow_float rvec[5] = { chi_eos + gm1*ek, -gm1*vx, -gm1*vy, -gm1*vz, gm1 };
        const double dbeta = static_cast<double>(beta);
        const double alpha = (1.0 - dbeta) / (dbeta * static_cast<double>(local_sonic) * static_cast<double>(local_sonic));
        const double gam = static_cast<double>(v_over_dtau) * alpha;   // = (V/Δτ')·α

        // D0 を float で 2 RHS 同時に解く: y=D0⁻¹b, z=D0⁻¹g。
        flow_float y[5], z[5];
        const bool ok = block_dplur::solve_5x5_2rhs(D0, b, gvec, y, z);

        // Sherman-Morrison スカラー (悪条件 1/β はここだけ double): x = y - [γ(rᵀy)/(1+γ(rᵀz))] z。
        double ry = 0.0, rz = 0.0;
        #pragma unroll
        for (int i = 0; i < 5; ++i) { ry += static_cast<double>(rvec[i]) * y[i]; rz += static_cast<double>(rvec[i]) * z[i]; }
        const double denom = 1.0 + gam * rz;
        const double sfac = (fabs(denom) > 1.0e-300) ? gam * ry / denom : 0.0;

        flow_float correction[5];
        #pragma unroll
        for (int i = 0; i < 5; ++i)
            correction[i] = ok ? static_cast<flow_float>((static_cast<double>(y[i]) - sfac * static_cast<double>(z[i]))
                                                         * static_cast<double>(implicit_relax))
                               : static_cast<flow_float>(0.0);

        block_dplur::store_block_vec(ic, correction, dq_new_0, dq_new_1, dq_new_2, dq_new_3, dq_new_4);
    }
}

// block DPLUR の sweep 間バッファ入れ替え。ドライバ側から各 sweep 後に明示的に呼ぶ
// （旧実装は wrapper 内部で暗黙に swap していたが、古典 DPLUR では制御フローを明示化する）。
// 近傍 dq の AoS バッファ (stride 8)。wrapper で nCells_all に合わせて確保し、swap で old/new を入れ替える。
static flow_float* g_dqPackOld = nullptr;
static flow_float* g_dqPackNew = nullptr;
static geom_int    g_dqPackN   = 0;

void swapBlockImplicitCorrectionBuffers(variables& var)
{
    std::swap(g_dqPackOld, g_dqPackNew);
    std::swap(var.c_d["dq_block_old_0"], var.c_d["dq_block_new_0"]);
    std::swap(var.c_d["dq_block_old_1"], var.c_d["dq_block_new_1"]);
    std::swap(var.c_d["dq_block_old_2"], var.c_d["dq_block_new_2"]);
    std::swap(var.c_d["dq_block_old_3"], var.c_d["dq_block_new_3"]);
    std::swap(var.c_d["dq_block_old_4"], var.c_d["dq_block_new_4"]);
}

// scalar 対角版 (blockDPLUR==0) の sweep 間バッファ入れ替え。block 版と同様にドライバ側から呼ぶ。
void swapScalarImplicitCorrectionBuffers(variables& var)
{
    std::swap(var.c_d["dq_ro_old"],   var.c_d["dq_ro_new"]);
    std::swap(var.c_d["dq_roUx_old"], var.c_d["dq_roUx_new"]);
    std::swap(var.c_d["dq_roUy_old"], var.c_d["dq_roUy_new"]);
    std::swap(var.c_d["dq_roUz_old"], var.c_d["dq_roUz_new"]);
    std::swap(var.c_d["dq_roe_old"],  var.c_d["dq_roe_new"]);
}

//TODO __global__ void runge_kutta_dual_explicit_d
//TODO // see https://sci-hub.se/https://doi.org/10.1016/j.compfluid.2003.10.004
//TODO // N: previous outer step , M: previous inner loop
//TODO ( 
//TODO  geom_int dt ,
//TODO 
//TODO  // mesh structure
//TODO  geom_int nCells_all , geom_int nCells,
//TODO  geom_float* vol ,
//TODO 
//TODO  // variables
//TODO  flow_float* ro  ,
//TODO  flow_float* roUx  ,
//TODO  flow_float* roUy  ,
//TODO  flow_float* roUz  ,
//TODO  flow_float* roe  ,
//TODO 
//TODO  flow_float* roN ,
//TODO  flow_float* roUxN ,
//TODO  flow_float* roUyN ,
//TODO  flow_float* roUzN ,
//TODO  flow_float* roeN ,
//TODO 
//TODO  flow_float* roM ,
//TODO  flow_float* roUxM ,
//TODO  flow_float* roUyM ,
//TODO  flow_float* roUzM ,
//TODO  flow_float* roeM ,
//TODO 
//TODO  flow_float* res_ro,
//TODO  flow_float* res_roUx,
//TODO  flow_float* res_roUy,
//TODO  flow_float* res_roUz,
//TODO  flow_float* res_roe,
//TODO 
//TODO  flow_float* res_ro_dual ,
//TODO  flow_float* res_roUx_dual ,
//TODO  flow_float* res_roUy_dual ,
//TODO  flow_float* res_roUz_dual ,
//TODO  flow_float* res_roe_dual 
//TODO 
//TODO )
//TODO {
//TODO     geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
//TODO 
//TODO     geom_float v = vol[ic];
//TODO 
//TODO     if (ic < nCells) {
//TODO         // N: previous outer step , M: previous inner loop
//TODO         res_ro_dual[ic]   = -(ro[ic]-roN[ic])*v/dt     + res_ro[ic];
//TODO         res_roUx_dual[ic] = -(roUx[ic]-roUxN[ic])*v/dt + res_roUx[ic];
//TODO         res_roUy_dual[ic] = -(roUy[ic]-roUyN[ic])*v/dt + res_roUy[ic];
//TODO         res_roUz_dual[ic] = -(roUz[ic]-roUzN[ic])*v/dt + res_roUz[ic];
//TODO         res_roe_dual[ic]  = -(roe[ic]-roeN[ic])*v/dt   + res_roe[ic];
//TODO     }
//TODO     __syncthreads();
//TODO }

// ---- ライン行列の書き出し (plan time_integration-line-viscous-jacobian §6.2、デバッグ) ----
// FORGE_LINE_DUMP_DIR=<dir> で有効。FORGE_LINE_DUMP_CALL (既定 1) 回目の factor の直前に、FORGE_LINE_DUMP_NODES (節点番号のカンマ区切り) を含む
// ラインの全節点について D・Kprev・Knext・ライン面のスカラー粘性の和・dt_local・体積・拘束のフラグ・状態を書き、続く各 sweep の solve の後に rhs・dq を書く。
// 形式: <dir>/meta.txt (名前・行数・列数) と <dir>/<name>.f64 (float64、(節点, 列) の行優先)。
namespace line_dump {
struct State {
    bool init = false, on = false;
    std::string dir;
    int targetCall = 1, factorCalls = 0, sweep = -1;
    std::vector<int> nodes;          // 書き出す節点 (ラインの順)
    std::vector<int> lineOf;         // 各節点のライン番号
    flow_float* viscBuf = nullptr;
};
static State g;
static void initOnce(mesh& msh) {
    if (g.init) return;
    g.init = true;
    const char* d = getenv("FORGE_LINE_DUMP_DIR");
    if (!d || !*d || msh.nImplicitLines <= 0) return;
    g.on = true; g.dir = d;
    if (const char* c = getenv("FORGE_LINE_DUMP_CALL")) g.targetCall = atoi(c);
    std::vector<int> want;
    if (const char* n = getenv("FORGE_LINE_DUMP_NODES")) {
        std::stringstream ss(n); std::string t;
        while (std::getline(ss, t, ',')) if (!t.empty()) want.push_back(atoi(t.c_str()));
    }
    std::vector<geom_int> off(msh.nImplicitLines + 1);
    gpuErrchk(cudaMemcpy(off.data(), msh.line_offsets_d, sizeof(geom_int) * (msh.nImplicitLines + 1), cudaMemcpyDeviceToHost));
    std::vector<geom_int> cells(off.back());
    gpuErrchk(cudaMemcpy(cells.data(), msh.line_cells_d, sizeof(geom_int) * off.back(), cudaMemcpyDeviceToHost));
    for (geom_int l = 0; l < msh.nImplicitLines; ++l) {
        bool hit = false;
        for (geom_int p = off[l]; p < off[l + 1] && !hit; ++p)
            for (int w : want) if (cells[p] == w) { hit = true; break; }
        if (!hit) continue;
        for (geom_int p = off[l]; p < off[l + 1]; ++p) { g.nodes.push_back(cells[p]); g.lineOf.push_back(l); }
    }
    gpuErrchk(cudaMalloc((void**)&g.viscBuf, sizeof(flow_float) * msh.nCells_all));
    gpuErrchk(cudaMemset(g.viscBuf, 0, sizeof(flow_float) * msh.nCells_all));
    printf("[lineDump] %s: factor %d 回目、%zu 節点 (要求 %zu 節点を含むライン)\n", g.dir.c_str(), g.targetCall, g.nodes.size(), want.size());
}
template<typename T>
static void put(const std::string& name, const std::vector<const T*>& cols, size_t nAll) {
    std::vector<double> out(g.nodes.size() * cols.size());
    std::vector<T> h(nAll);
    for (size_t c = 0; c < cols.size(); ++c) {
        if (cols[c] == nullptr) { for (size_t i = 0; i < g.nodes.size(); ++i) out[i * cols.size() + c] = 0.0; continue; }
        gpuErrchk(cudaMemcpy(h.data(), cols[c], sizeof(T) * nAll, cudaMemcpyDeviceToHost));
        for (size_t i = 0; i < g.nodes.size(); ++i) out[i * cols.size() + c] = (double)h[g.nodes[i]];
    }
    std::ofstream f(g.dir + "/" + name + ".f64", std::ios::binary);
    f.write((const char*)out.data(), sizeof(double) * out.size());
    std::ofstream m(g.dir + "/meta.txt", std::ios::app);
    m << name << " " << g.nodes.size() << " " << cols.size() << "\n";
}
// AoS の [n*25] 配列 (Kprev/Knext) を 25 列に分けて書く
static void putAoS25(const std::string& name, const flow_float* a, size_t nAll) {
    std::vector<flow_float> h(nAll * 25);
    gpuErrchk(cudaMemcpy(h.data(), a, sizeof(flow_float) * nAll * 25, cudaMemcpyDeviceToHost));
    std::vector<double> out(g.nodes.size() * 25);
    for (size_t i = 0; i < g.nodes.size(); ++i) for (int k = 0; k < 25; ++k) out[i * 25 + k] = (double)h[(size_t)g.nodes[i] * 25 + k];
    std::ofstream f(g.dir + "/" + name + ".f64", std::ios::binary);
    f.write((const char*)out.data(), sizeof(double) * out.size());
    std::ofstream m(g.dir + "/meta.txt", std::ios::app);
    m << name << " " << g.nodes.size() << " 25\n";
}
static void atFactor(solverConfig& cfg, mesh& msh, variables& var) {
    initOnce(msh);
    if (!g.on) return;
    ++g.factorCalls;
    g.sweep = -1;
    if (g.factorCalls != g.targetCall) return;
    g.sweep = 0;
    const size_t n = msh.nCells_all;
    { std::ofstream m(g.dir + "/meta.txt"); m << "# name rows cols (float64、行優先)。implicitRelax " << cfg.implicitRelax
                                              << " lineViscCoupling " << cfg.lineViscCoupling << " implicitThermalJacobian " << cfg.implicitThermalJacobian << "\n"; }
    std::vector<double> ids(g.nodes.size() * 2);
    for (size_t i = 0; i < g.nodes.size(); ++i) { ids[2 * i] = g.nodes[i]; ids[2 * i + 1] = g.lineOf[i]; }
    { std::ofstream f(g.dir + "/node_line.f64", std::ios::binary); f.write((const char*)ids.data(), sizeof(double) * ids.size());
      std::ofstream m(g.dir + "/meta.txt", std::ios::app); m << "node_line " << g.nodes.size() << " 2\n"; }
    std::vector<const flow_float*> D;
    for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) D.push_back(var.c_d["diag_block_" + std::to_string(i) + std::to_string(j)]);
    put<flow_float>("D", D, n);
    putAoS25("Kprev", msh.line_Kprev_d, n);
    putAoS25("Knext", msh.line_Knext_d, n);
    put<flow_float>("scalar_visc_line", {g.viscBuf}, n);
    put<flow_float>("dt_vol", {var.c_d["dt_local"], var.c_d["volume"]}, n);
    put<geom_int>("flags_wall_iso_axis", {msh.wall_flag_d, msh.iso_wall_flag_d, msh.axis_flag_d}, (size_t)msh.nCells);   // フラグは節点の範囲 (nCells) で確保
    put<flow_float>("state_ro_roU_roe_cp_gamma", {var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"], var.c_d["cp"], var.c_d["gamma"]}, n);
    printf("[lineDump] factor の直前を書いた (%zu 節点)\n", g.nodes.size());
}
static void afterSolve(mesh& msh, variables& var) {
    if (!g.on || g.sweep < 0) return;
    const size_t n = msh.nCells_all;
    const std::string k = std::to_string(g.sweep);
    put<flow_float>("rhs_s" + k, {var.c_d["rhs_block_0"], var.c_d["rhs_block_1"], var.c_d["rhs_block_2"], var.c_d["rhs_block_3"], var.c_d["rhs_block_4"]}, n);
    put<flow_float>("dqnew_s" + k, {var.c_d["dq_block_new_0"], var.c_d["dq_block_new_1"], var.c_d["dq_block_new_2"], var.c_d["dq_block_new_3"], var.c_d["dq_block_new_4"]}, n);
    put<flow_float>("dqold_s" + k, {var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"], var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"]}, n);
    ++g.sweep;
}
} // namespace line_dump

// 診断: 値 2/3 の薄層の項のマスク (FORGE_LVC_TERMS、既定 7 = 全部。plan time_integration-line-viscous-jacobian §6.8)
static int lineViscTermsEnv() {
    static const int v = [](){ const char* e = getenv("FORGE_LVC_TERMS"); const int t = (e && *e) ? atoi(e) : 7;
                               if (t != 7) printf("[lineViscCoupling] 診断のマスク FORGE_LVC_TERMS=%d (1 運動量 D/K、2 熱伝導の近傍 K、4 仕事 D/K、8 熱伝導の K の密度の列を外す)\n", t);
                               return t; }();
    return v;
}

static flow_float* lineDumpViscBuf(mesh& msh) {
    line_dump::initOnce(msh);
    return line_dump::g.on ? line_dump::g.viscBuf : nullptr;
}

void timeIntegration_d_wrapper(int loop , solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , int lineStoreK)
{
    // 軸対称エンコード: 0=非軸対称 / 1=r 重み方式 (hoop Jacobian) / 2=SU2 流 planar+ソース (SU2 4x4 Jacobian)。
    // ==1 判定しかしない旧カーネル (scalar/lowmach) は 2 のとき軸対称 Jacobian を持たない (source lag, 定常解不変)。
    // 診断 (env): FORGE_DIAG_SU2JAC_OFF=1 で SU2 ソース Jacobian を落とす (ソースは残す = 完全 lag)。
    static const bool diagSu2JacOff = (getenv("FORGE_DIAG_SU2JAC_OFF") != nullptr);
    int axisymEnc = (cfg.isAxisymmetric == 1) ? ((cfg.axisymMethod == 1) ? 2 : 1) : 0;
    if (diagSu2JacOff && axisymEnc == 2) axisymEnc = 0;
    // 診断 near-axis 安定化係数を env から 1 度だけ device へ設定 (既定 0 = 不変)。
    static bool s_axisAlphaInit = false;
    if (!s_axisAlphaInit) {
        float a = 0.0f;
        if (const char* e = getenv("FORGE_AXIS_DIAG_ALPHA")) a = static_cast<float>(atof(e));
        cudaMemcpyToSymbol(g_axisDiagAlpha, &a, sizeof(float));
        s_axisAlphaInit = true;
    }
    if (cfg.timeIntegration == 4) { // 4th order runge kutta
        runge_kutta_exp_4th_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
            loop, 
            cfg.coef_DT_4thRunge[loop],
            cfg.coef_Res_4thRunge[loop],
            cfg.dt ,
            var.c_d["dt_local"],

            // mesh structure
            msh.nCells_all , msh.nCells ,
            var.c_d["volume"],

            // basic variables
            var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"] , var.c_d["roe"] ,
            var.c_d["roN"] , var.c_d["roUxN"], var.c_d["roUyN"] , var.c_d["roUzN"], var.c_d["roeN"] ,
            var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"], var.c_d["roeM"] ,
            var.c_d["res_ro"]  , var.c_d["res_roUx"]  , var.c_d["res_roUy"]  , var.c_d["res_roUz"] , var.c_d["res_roe"] ,
            var.c_d["res_ro_m"], var.c_d["res_roUx_m"], var.c_d["res_roUy_m"], var.c_d["res_roUz_m"] , var.c_d["res_roe_m"] 
        ) ;

    } else if (cfg.timeIntegration == 1 or cfg.timeIntegration == 3) { // explicit
        runge_kutta_exp_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
            loop,
            cfg.coef_N[loop],
            cfg.coef_M[loop],
            cfg.coef_Res[loop],
            cfg.dt , 
            var.c_d["dt_local"],

            // mesh structure
            msh.nCells_all , msh.nCells ,
            var.c_d["volume"],

            // basic variables
            var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"] , var.c_d["roe"] ,
            var.c_d["roN"] , var.c_d["roUxN"], var.c_d["roUyN"] , var.c_d["roUzN"], var.c_d["roeN"] ,
            var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"], var.c_d["roeM"] ,
            var.c_d["res_ro"]  , var.c_d["res_roUx"]  , var.c_d["res_roUy"]  , var.c_d["res_roUz"] , var.c_d["res_roe"] 
        ) ;
    } else if (cfg.timeIntegration == 11) { // implicit defect-correction with diagonal Jacobian approximation
        if (cfg.blockDPLUR == 1) {
            // レジスタ過多のため専用の小さい block サイズで起動（__launch_bounds__ と整合）。
            const int block_threads = BLOCK_DPLUR_THREADS;
            const int block_grid = (msh.nCells_all + block_threads - 1) / block_threads;
            if (cfg.lowMachPrecond >= 2) {
              // Phase 4: 完全 Γ⁻¹A 前処理の倍精度カーネル (dt_local は前処理 Δτ' に拡大済)。
              // lowMachPrecond==2: RHS 散逸 c' (Phase 1) + 本 LHS 前処理。
              // lowMachPrecond==3: LHS 前処理のみ (RHS 散逸は c_hat で =0 と不変)。本カーネルは
              //   擬似時間項 Γ_c のみ前処理しフラックス A_c は非前処理ゆえ、収束解は前処理なしと一致する。
              implicit_defect_correction_block_precond_d<<<block_grid , block_threads>>>(
                loop, cfg.dt, var.c_d["dt_local"], cfg.implicitRelax, var.c_d["gamma"], cfg.precondEps,
                msh.nCells_all, msh.nCells, var.c_d["volume"],
                msh.map_plane_cells_d, msh.map_cell_planes_index_d, msh.map_cell_planes_d,
                var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["ge_x"], var.p_d["ge_y"], var.p_d["ge_z"],
                var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
                var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"],
                cfg.visc, var.c_d["vis_turb"], var.c_d["sonic"],
                var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"], var.c_d["Ht"],
                var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"],
                var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"], var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"],
                var.c_d["dq_block_new_0"], var.c_d["dq_block_new_1"], var.c_d["dq_block_new_2"], var.c_d["dq_block_new_3"], var.c_d["dq_block_new_4"],
                axisymEnc,
                (cfg.isAxisymmetric == 1) ? ((cfg.axisRFloor > (flow_float)0.0 || cfg.hoopAreaFromClosure == 1) ? var.c_d["A_closure_y"] : var.c_d["A_planar"]) : var.c_d["volume"],
                cfg.axisRFloor,
                cfg.unsteadyDiagCoef
              );
            } else {
            // implicitSolvePrecision: 0=float (既定・高速), 1=double (軸対称近軸の根治, 遅い)。
            // 同じテンプレートカーネルを ST=float/double で起動。引数は共通 (FORGE_BDPLUR_ARGS)。
            #define FORGE_BDPLUR_ARGS \
                loop, cfg.dt, var.c_d["dt_local"], cfg.implicitRelax, var.c_d["gamma"], \
                (cfg.thermalMethod == 2 ? 1 : 0), \
                msh.nCells_all, msh.nCells, var.c_d["volume"], \
                msh.map_plane_cells_d, msh.map_cell_planes_index_d, msh.map_cell_planes_d, \
                var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"], \
                var.p_d["ge_x"], var.p_d["ge_y"], var.p_d["ge_z"],  /* 面ごとの差 e (段 ③) */ \
                var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"], \
                var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"], \
                cfg.visc, var.c_d["vis_turb"], var.c_d["sonic"], \
                var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"], var.c_d["Ht"], \
                var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"], \
                var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"], var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"], \
                var.c_d["dq_block_new_0"], var.c_d["dq_block_new_1"], var.c_d["dq_block_new_2"], var.c_d["dq_block_new_3"], var.c_d["dq_block_new_4"], \
                var.c_d["rhs_block_0"], var.c_d["rhs_block_1"], var.c_d["rhs_block_2"], var.c_d["rhs_block_3"], var.c_d["rhs_block_4"], \
                var.c_d["diag_block_00"], var.c_d["diag_block_01"], var.c_d["diag_block_02"], var.c_d["diag_block_03"], var.c_d["diag_block_04"], \
                var.c_d["diag_block_10"], var.c_d["diag_block_11"], var.c_d["diag_block_12"], var.c_d["diag_block_13"], var.c_d["diag_block_14"], \
                var.c_d["diag_block_20"], var.c_d["diag_block_21"], var.c_d["diag_block_22"], var.c_d["diag_block_23"], var.c_d["diag_block_24"], \
                var.c_d["diag_block_30"], var.c_d["diag_block_31"], var.c_d["diag_block_32"], var.c_d["diag_block_33"], var.c_d["diag_block_34"], \
                var.c_d["diag_block_40"], var.c_d["diag_block_41"], var.c_d["diag_block_42"], var.c_d["diag_block_43"], var.c_d["diag_block_44"], \
                axisymEnc, (cfg.isAxisymmetric == 1) ? ((cfg.axisRFloor > (flow_float)0.0 || cfg.hoopAreaFromClosure == 1) ? var.c_d["A_closure_y"] : var.c_d["A_planar"]) : var.c_d["volume"], cfg.axisRFloor, cfg.unsteadyDiagCoef, \
                nullptr,  /* axis_flag: 旧 nodeAxisDirichlet の全 5 行 decouple (撤去) */ \
                ((cfg.discretization == "node" && cfg.isAxisymmetric == 1) ? msh.axis_flag_d : nullptr),  /* axis_ur_flag: 軸ノードの roUy 行 decouple (常時) */ \
                ((cfg.discretization == "node" && cfg.isAxisymmetric == 1) ? msh.axis_flag_d : nullptr),  /* axis_flag_src: SU2 流 (enc==2) の軸ソース Jacobian ガード (軸ノードはソース 0) */ \
                ((cfg.discretization == "node" && cfg.nodeWallDirichlet == 1) ? msh.wall_flag_d : nullptr),  /* wall_flag: 壁運動量3行 decouple */ \
                ((cfg.discretization == "node" && cfg.nodeWallDirichlet == 1 && cfg.nodeIsothermalEnergyBC != 1) ? msh.iso_wall_flag_d : nullptr),  /* iso_wall_flag: 等温壁 roe 行 decouple (T ピンと対)。弱形式 (nodeIsothermalEnergyBC=1) では単位行化も rowDec も外す (plan boundary-weak-isothermal-wall §4.3) */ \
                (weakIsoWall::active(cfg, msh) ? weakIsoWall::diagBuf(msh) : nullptr),  /* weakIsoDiag: 弱形式の近似対角 */ \
                (flow_float)(cfg.cp / max(cfg.gamma, 1.0e-30)),  /* weakIsoCv = cp/gamma = c_v (CPG) */ \
                ((cfg.discretization == "node") ? 1 : 0),  /* isNode: 5e 境界半割面の粘性対角スキップ */ \
                ((cfg.lineImplicit == 1) ? msh.line_prev_d : nullptr), \
                ((cfg.lineImplicit == 1) ? msh.line_next_d : nullptr), \
                msh.line_Kprev_d, msh.line_Knext_d, (((loop == 0) && (lineStoreK != 0)) ? 1 : 0), cfg.lineViscCoupling,  /* line-implicit */ \
                ((cfg.implicitSolvePrecision == 0 && cfg.lineImplicit == 0 && cfg.blockDPLURDiagCache != 0) ? 1 : 0),  /* useDiagCache: float・point 経路のみ */ \
                (usePack ? (const flow_float*)g_dqPackOld : nullptr), (usePack ? g_dqPackNew : nullptr),  /* 近傍 dq の AoS 版 */ \
                cfg.implicitThermalJacobian,  /* エネルギー行の熱伝導 Jacobian / 等温壁の拘束の行 (ビットマスク) */ \
                (((cfg.implicitThermalJacobian & 1) || cfg.lineViscCoupling >= 2) ? var.c_d["thermCond"] : nullptr), \
                (((cfg.implicitThermalJacobian & 1) || cfg.lineViscCoupling >= 2) ? var.c_d["cp"] : nullptr), \
                (((cfg.implicitThermalJacobian & 1) || cfg.lineViscCoupling >= 2) ? var.p_d["fx"] : nullptr), \
                cfg.turbulentPrandtl, \
                ((cfg.lineViscCoupling >= 2) ? var.c_d["vis_lam"] : nullptr),  /* 節点ごとの層流粘性 (残差の μ_f と揃える) */ \
                lineDumpViscBuf(msh),  /* デバッグの書き出し (FORGE_LINE_DUMP_DIR のときだけ非 nullptr) */ \
                lineViscTermsEnv()     /* 診断のマスク (FORGE_LVC_TERMS、既定 7) */
            // 近傍 dq の AoS 経路: line-implicit と node 周期 (SoA だけを直接書き換える) では使わない。
            const bool usePack = (cfg.lineImplicit == 0) && (cfg.blockDPLURDqPack != 0) &&
                                 !(cfg.discretization == "node" && msh.periodicRoot_d != nullptr && msh.nPeriodicMembers > 0);
            if (usePack && (g_dqPackOld == nullptr || g_dqPackN != msh.nCells_all)) {
                if (g_dqPackOld) { cudaFree(g_dqPackOld); cudaFree(g_dqPackNew); }
                const size_t nb = (size_t)msh.nCells_all * 8 * sizeof(flow_float);
                gpuErrchk(cudaMalloc((void**)&g_dqPackOld, nb)); gpuErrchk(cudaMalloc((void**)&g_dqPackNew, nb));
                gpuErrchk(cudaMemset(g_dqPackOld, 0, nb)); gpuErrchk(cudaMemset(g_dqPackNew, 0, nb));
                g_dqPackN = msh.nCells_all;
            }
            if (cfg.implicitSolvePrecision == 1)
                implicit_defect_correction_block_d<double><<<block_grid , block_threads>>>(FORGE_BDPLUR_ARGS);
            else
                implicit_defect_correction_block_d<float><<<block_grid , block_threads>>>(FORGE_BDPLUR_ARGS);
            #undef FORGE_BDPLUR_ARGS
            }
        } else {
            implicit_defect_correction_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
                loop,
                cfg.dt,
                var.c_d["dt_local"],
                cfg.implicitRelax,
                var.c_d["gamma"],
                msh.nCells_all,
                msh.nCells,
                var.c_d["volume"],
                msh.map_plane_cells_d,
                msh.map_cell_planes_index_d,
                msh.map_cell_planes_d,
                var.c_d["ccx"],
                var.c_d["ccy"],
                var.c_d["ccz"],
                var.p_d["ge_x"],
                var.p_d["ge_y"],
                var.p_d["ge_z"],
                var.p_d["sx"],
                var.p_d["sy"],
                var.p_d["sz"],
                var.p_d["ss"],
                var.c_d["ro"],
                var.c_d["roUx"],
                var.c_d["roUy"],
                var.c_d["roUz"],
                var.c_d["roe"],
                var.c_d["roN"],
                var.c_d["roUxN"],
                var.c_d["roUyN"],
                var.c_d["roUzN"],
                var.c_d["roeN"],
                cfg.visc,
                var.c_d["vis_turb"],
                var.c_d["sonic"],
                var.c_d["Ux"],
                var.c_d["Uy"],
                var.c_d["Uz"],
                var.c_d["res_ro"],
                var.c_d["res_roUx"],
                var.c_d["res_roUy"],
                var.c_d["res_roUz"],
                var.c_d["res_roe"],
                var.c_d["dq_ro_old"],
                var.c_d["dq_roUx_old"],
                var.c_d["dq_roUy_old"],
                var.c_d["dq_roUz_old"],
                var.c_d["dq_roe_old"],
                var.c_d["dq_ro_new"],
                var.c_d["dq_roUx_new"],
                var.c_d["dq_roUy_new"],
                var.c_d["dq_roUz_new"],
                var.c_d["dq_roe_new"],
                axisymEnc,
                (cfg.isAxisymmetric == 1) ? ((cfg.axisRFloor > (flow_float)0.0 || cfg.hoopAreaFromClosure == 1) ? var.c_d["A_closure_y"] : var.c_d["A_planar"]) : var.c_d["volume"],
                cfg.axisRFloor,
                cfg.unsteadyDiagCoef
            );
        }
        // 古典 DPLUR: buffer swap と Q への commit はドライバ側 (main.cpp blockDPLURSolve /
        // applyBlockImplicitCorrection) で明示的に行う。ここでは sweep カーネルの起動のみ。
//TODO    } else if (cfg.timeIntegration == 10) { // implicit (m-time stepping & explicit scheme)
//TODO        runge_kutta_dual_explicit_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
//TODO            cfg.dt , 
//TODO
//TODO            // mesh structure
//TODO            msh.nCells_all , msh.nCells ,
//TODO            var.c_d["volume"],
//TODO
//TODO            // basic variables
//TODO            var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"] , var.c_d["roe"] ,
//TODO            var.c_d["roN"] , var.c_d["roUxN"], var.c_d["roUyN"] , var.c_d["roUzN"], var.c_d["roeN"] ,
//TODO            var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"], var.c_d["roeM"] ,
//TODO            var.c_d["res_ro"]  , var.c_d["res_roUx"]  , var.c_d["res_roUy"]  , var.c_d["res_roUz"] , var.c_d["res_roe"] ,
//TODO            var.c_d["res_ro_m"], var.c_d["res_roUx_m"], var.c_d["res_roUy_m"], var.c_d["res_roUz_m"] , var.c_d["res_roe_m"] 
//TODO        ) ;
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

}
// =============================================================================
// line-implicit: ライン block-Thomas (plans/active/time_integration-line-implicit.md)。
// sweep カーネルが保存した diag (loop0 凍結)・K (loop0)・rhs (毎 sweep, ライン外 lag 込み) から、
// 各ラインの block 三重対角系
//   D_k ΔQ_k − Kprev_k ΔQ_{k-1} − Knext_k ΔQ_{k+1} = rhs_k
// を前進消去+後退代入で厳密に解き dq_new を上書きする。1 ライン = 1 スレッド (v1)、内部 double。
// scratch: W_k = D̃_k⁻¹ Knext_k (25/cell), y_k = D̃_k⁻¹ b̃_k (5/cell)。
namespace line_implicit {

__device__ __forceinline__ bool lu5_factor(double A[5][5], int piv[5])
{
    for (int col = 0; col < 5; ++col) {
        int pv = col; double pa = fabs(A[col][col]);
        for (int r = col + 1; r < 5; ++r) { const double c = fabs(A[r][col]); if (c > pa) { pv = r; pa = c; } }
        if (pa < 1.0e-30) return false;
        piv[col] = pv;
        if (pv != col) for (int k = 0; k < 5; ++k) { const double tmp = A[col][k]; A[col][k] = A[pv][k]; A[pv][k] = tmp; }
        const double inv = 1.0 / A[col][col];
        for (int r = col + 1; r < 5; ++r) {
            const double f = A[r][col] * inv;
            A[r][col] = f;                      // L を下三角に格納
            for (int k = col + 1; k < 5; ++k) A[r][k] -= f * A[col][k];
        }
    }
    return true;
}

__device__ __forceinline__ void lu5_solve(const double A[5][5], const int piv[5], double x[5])
{
    // LAPACK getrs 流: ① 行交換を全て先に適用 (LASWP) ② 単位下三角 L 前進代入 ③ U 後退代入。
    // 交換と代入をインタリーブする書き方は、後段ピボットが L 部分も行交換する getrf 形格納と
    // 非整合で誤解を返す (2026-09-02 に numpy 照合で確認済みの罠)。
    for (int col = 0; col < 5; ++col) {
        if (piv[col] != col) { const double tmp = x[col]; x[col] = x[piv[col]]; x[piv[col]] = tmp; }
    }
    for (int col = 0; col < 5; ++col) {
        for (int r = col + 1; r < 5; ++r) x[r] -= A[r][col] * x[col];
    }
    for (int r = 4; r >= 0; --r) {
        double s = x[r];
        for (int c = r + 1; c < 5; ++c) s -= A[r][c] * x[c];
        x[r] = s / A[r][r];
    }
}

} // namespace line_implicit

__global__ void lineThomas_d
(
 geom_int nLines,
 const geom_int* line_offsets,
 const geom_int* line_cells,
 const flow_float* Kprev, const flow_float* Knext,
 const flow_float* d00, const flow_float* d01, const flow_float* d02, const flow_float* d03, const flow_float* d04,
 const flow_float* d10, const flow_float* d11, const flow_float* d12, const flow_float* d13, const flow_float* d14,
 const flow_float* d20, const flow_float* d21, const flow_float* d22, const flow_float* d23, const flow_float* d24,
 const flow_float* d30, const flow_float* d31, const flow_float* d32, const flow_float* d33, const flow_float* d34,
 const flow_float* d40, const flow_float* d41, const flow_float* d42, const flow_float* d43, const flow_float* d44,
 const flow_float* rhs0, const flow_float* rhs1, const flow_float* rhs2, const flow_float* rhs3, const flow_float* rhs4,
 const flow_float* dq_old_0, const flow_float* dq_old_1, const flow_float* dq_old_2, const flow_float* dq_old_3, const flow_float* dq_old_4,
 flow_float* dq_new_0, flow_float* dq_new_1, flow_float* dq_new_2, flow_float* dq_new_3, flow_float* dq_new_4,
 flow_float implicit_relax,
 double* Wd, double* yd
)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    const geom_int b = line_offsets[l];
    const geom_int e = line_offsets[l + 1];
    bool fail = false;

    // ---- 前進消去 ----
    for (geom_int p = b; p < e && !fail; ++p) {
        const geom_int ic = line_cells[p];
        double M[5][5] = {
            {(double)d00[ic],(double)d01[ic],(double)d02[ic],(double)d03[ic],(double)d04[ic]},
            {(double)d10[ic],(double)d11[ic],(double)d12[ic],(double)d13[ic],(double)d14[ic]},
            {(double)d20[ic],(double)d21[ic],(double)d22[ic],(double)d23[ic],(double)d24[ic]},
            {(double)d30[ic],(double)d31[ic],(double)d32[ic],(double)d33[ic],(double)d34[ic]},
            {(double)d40[ic],(double)d41[ic],(double)d42[ic],(double)d43[ic],(double)d44[ic]}};
        double bk[5] = {(double)rhs0[ic],(double)rhs1[ic],(double)rhs2[ic],(double)rhs3[ic],(double)rhs4[ic]};
        if (p > b) {
            const geom_int icm = line_cells[p - 1];
            // M -= Kprev·W_{k-1},  b̃_k = b_k + Kprev·y_{k-1}
            // (標準形 L=−Kprev, U=−Knext につき符号は加算側に出る)
            double Kp[5][5];
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < 5; ++j)
                    Kp[i][j] = (double)Kprev[(size_t)ic * 25 + i * 5 + j];
            for (int i = 0; i < 5; ++i) {
                double bacc = 0.0;
                for (int m = 0; m < 5; ++m) bacc += Kp[i][m] * yd[(size_t)icm * 5 + m];
                bk[i] += bacc;
                for (int j = 0; j < 5; ++j) {
                    double macc = 0.0;
                    for (int m = 0; m < 5; ++m) macc += Kp[i][m] * Wd[(size_t)icm * 25 + m * 5 + j];
                    M[i][j] -= macc;
                }
            }
        }
        int piv[5];
        if (!line_implicit::lu5_factor(M, piv)) { fail = true; break; }
        line_implicit::lu5_solve(M, piv, bk);                  // y_k
        for (int i = 0; i < 5; ++i) yd[(size_t)ic * 5 + i] = bk[i];
        if (p + 1 < e) {                                        // W_k = M⁻¹·Knext_k
            for (int j = 0; j < 5; ++j) {
                double col[5];
                for (int i = 0; i < 5; ++i) col[i] = (double)Knext[(size_t)ic * 25 + i * 5 + j];
                line_implicit::lu5_solve(M, piv, col);
                for (int i = 0; i < 5; ++i) Wd[(size_t)ic * 25 + i * 5 + j] = col[i];
            }
        }
    }

    // ---- 後退代入 (relax を掛けて dq_new へ) / 失敗時は前回反復値を保持 ----
    if (fail) {
        for (geom_int p = b; p < e; ++p) {
            const geom_int ic = line_cells[p];
            dq_new_0[ic] = dq_old_0[ic]; dq_new_1[ic] = dq_old_1[ic]; dq_new_2[ic] = dq_old_2[ic];
            dq_new_3[ic] = dq_old_3[ic]; dq_new_4[ic] = dq_old_4[ic];
        }
        return;
    }
    double dq[5];
    {
        const geom_int ic = line_cells[e - 1];
        for (int i = 0; i < 5; ++i) dq[i] = yd[(size_t)ic * 5 + i];
        dq_new_0[ic] = (flow_float)(implicit_relax * dq[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * dq[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * dq[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * dq[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * dq[4]);
    }
    for (geom_int p = e - 2; p >= b; --p) {
        const geom_int ic = line_cells[p];
        double nx[5];
        for (int i = 0; i < 5; ++i) {
            double acc = yd[(size_t)ic * 5 + i];
            for (int m = 0; m < 5; ++m) acc += Wd[(size_t)ic * 25 + i * 5 + m] * dq[m];
            nx[i] = acc;
        }
        for (int i = 0; i < 5; ++i) dq[i] = nx[i];
        dq_new_0[ic] = (flow_float)(implicit_relax * dq[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * dq[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * dq[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * dq[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * dq[4]);
        if (p == b) break;   // geom_int が unsigned の場合の負回りガード
    }
}

// v2 (plans/active/time_integration-line-implicit-viscous-v2.md): factor/solve 分離。
// 前進消去の M̃_k = D_k − Kprev·W_{k−1} の構築・LU 分解・W_k = M̃⁻¹Knext は rhs に依存しない
// (D, K は storeLU 時に凍結) ので、storeLU のタイミングで 1 回だけ行い LU/piv/W を保存する。
// sweep 毎の solve は保存済み因子での代入 (Kp·y 25 積 + LASWP 前進/後退) だけになる。
// モノリシック版 (lineThomas_d) は毎 sweep この 5 列 solve + Kp·W (625 積) を再計算しており、
// これが DDES A/B での step 単価 2.44 倍の主犯 — 分離は厳密 (近似ゼロ) の最適化。
// INV (plan time_integration-line-implicit-speed §5.1 #7、FORGE_LINE_INV=1 の opt-in 実験): LU の代わりに逆行列 M̃⁻¹ を LUd に保存する。
// 逆行列は同じ部分ピボット付き LU で単位ベクトルを 5 回代入して作る (double)。W_k = M̃⁻¹Knext_k は従来どおり LU の代入で作る (不変)。
// solve の前進は lu5_solve (行交換・除算を含む直列の代入) の代わりに 5×5 の行列ベクトル積になる。pivd は使わない (恒等を書く)。
template<bool INV>
__global__ void lineThomasFactor_d
(
 geom_int nLines,
 const geom_int* line_offsets,
 const geom_int* line_cells,
 const flow_float* Kprev, const flow_float* Knext,
 const flow_float* d00, const flow_float* d01, const flow_float* d02, const flow_float* d03, const flow_float* d04,
 const flow_float* d10, const flow_float* d11, const flow_float* d12, const flow_float* d13, const flow_float* d14,
 const flow_float* d20, const flow_float* d21, const flow_float* d22, const flow_float* d23, const flow_float* d24,
 const flow_float* d30, const flow_float* d31, const flow_float* d32, const flow_float* d33, const flow_float* d34,
 const flow_float* d40, const flow_float* d41, const flow_float* d42, const flow_float* d43, const flow_float* d44,
 double* Wd, double* LUd, signed char* pivd, unsigned char* faild
)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    const geom_int b = line_offsets[l];
    const geom_int e = line_offsets[l + 1];
    faild[l] = 0;

    for (geom_int p = b; p < e; ++p) {
        const geom_int ic = line_cells[p];
        double M[5][5] = {
            {(double)d00[ic],(double)d01[ic],(double)d02[ic],(double)d03[ic],(double)d04[ic]},
            {(double)d10[ic],(double)d11[ic],(double)d12[ic],(double)d13[ic],(double)d14[ic]},
            {(double)d20[ic],(double)d21[ic],(double)d22[ic],(double)d23[ic],(double)d24[ic]},
            {(double)d30[ic],(double)d31[ic],(double)d32[ic],(double)d33[ic],(double)d34[ic]},
            {(double)d40[ic],(double)d41[ic],(double)d42[ic],(double)d43[ic],(double)d44[ic]}};
        if (p > b) {
            const geom_int icm = line_cells[p - 1];
            double Kp[5][5];
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < 5; ++j)
                    Kp[i][j] = (double)Kprev[(size_t)ic * 25 + i * 5 + j];
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < 5; ++j) {
                    double macc = 0.0;
                    for (int m = 0; m < 5; ++m) macc += Kp[i][m] * Wd[(size_t)icm * 25 + m * 5 + j];
                    M[i][j] -= macc;
                }
        }
        int piv[5];
        if (!line_implicit::lu5_factor(M, piv)) { faild[l] = 1; return; }
        if (INV) {
            for (int j = 0; j < 5; ++j) {                       // M̃⁻¹ の列 j = M̃⁻¹ e_j
                double col[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
                col[j] = 1.0;
                line_implicit::lu5_solve(M, piv, col);
                for (int i = 0; i < 5; ++i) LUd[(size_t)ic * 25 + i * 5 + j] = col[i];
            }
            for (int i = 0; i < 5; ++i) pivd[(size_t)ic * 5 + i] = (signed char)i;
        } else {
            for (int i = 0; i < 5; ++i) {
                pivd[(size_t)ic * 5 + i] = (signed char)piv[i];
                for (int j = 0; j < 5; ++j) LUd[(size_t)ic * 25 + i * 5 + j] = M[i][j];
            }
        }
        if (p + 1 < e) {                                        // W_k = M̃⁻¹·Knext_k
            for (int j = 0; j < 5; ++j) {
                double col[5];
                for (int i = 0; i < 5; ++i) col[i] = (double)Knext[(size_t)ic * 25 + i * 5 + j];
                line_implicit::lu5_solve(M, piv, col);
                for (int i = 0; i < 5; ++i) Wd[(size_t)ic * 25 + i * 5 + j] = col[i];
            }
        }
    }
}

template<bool INV>   // INV: LUd は逆行列 M̃⁻¹ (lineThomasFactor_d<true> が保存)、前進は行列ベクトル積
__global__ void lineThomasSolve_d
(
 geom_int nLines,
 const geom_int* line_offsets,
 const geom_int* line_cells,
 const flow_float* Kprev,
 const double* Wd, const double* LUd, const signed char* pivd, const unsigned char* faild,
 const flow_float* rhs0, const flow_float* rhs1, const flow_float* rhs2, const flow_float* rhs3, const flow_float* rhs4,
 const flow_float* dq_old_0, const flow_float* dq_old_1, const flow_float* dq_old_2, const flow_float* dq_old_3, const flow_float* dq_old_4,
 flow_float* dq_new_0, flow_float* dq_new_1, flow_float* dq_new_2, flow_float* dq_new_3, flow_float* dq_new_4,
 flow_float implicit_relax,
 double* yd
)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    const geom_int b = line_offsets[l];
    const geom_int e = line_offsets[l + 1];
    if (faild[l] != 0) {
        for (geom_int p = b; p < e; ++p) {
            const geom_int ic = line_cells[p];
            dq_new_0[ic] = dq_old_0[ic]; dq_new_1[ic] = dq_old_1[ic]; dq_new_2[ic] = dq_old_2[ic];
            dq_new_3[ic] = dq_old_3[ic]; dq_new_4[ic] = dq_old_4[ic];
        }
        return;
    }

    // ---- 前進 (保存因子で代入のみ) ----
    for (geom_int p = b; p < e; ++p) {
        const geom_int ic = line_cells[p];
        double bk[5] = {(double)rhs0[ic],(double)rhs1[ic],(double)rhs2[ic],(double)rhs3[ic],(double)rhs4[ic]};
        if (p > b) {
            const geom_int icm = line_cells[p - 1];
            for (int i = 0; i < 5; ++i) {
                double bacc = 0.0;
                for (int m = 0; m < 5; ++m)
                    bacc += (double)Kprev[(size_t)ic * 25 + i * 5 + m] * yd[(size_t)icm * 5 + m];
                bk[i] += bacc;
            }
        }
        if (INV) {
            for (int i = 0; i < 5; ++i) {
                double s = 0.0;
                for (int j = 0; j < 5; ++j) s += LUd[(size_t)ic * 25 + i * 5 + j] * bk[j];
                yd[(size_t)ic * 5 + i] = s;
            }
        } else {
            double M[5][5];
            int piv[5];
            for (int i = 0; i < 5; ++i) {
                piv[i] = (int)pivd[(size_t)ic * 5 + i];
                for (int j = 0; j < 5; ++j) M[i][j] = LUd[(size_t)ic * 25 + i * 5 + j];
            }
            line_implicit::lu5_solve(M, piv, bk);
            for (int i = 0; i < 5; ++i) yd[(size_t)ic * 5 + i] = bk[i];
        }
    }

    // ---- 後退代入 (relax を掛けて dq_new へ) ----
    double dq[5];
    {
        const geom_int ic = line_cells[e - 1];
        for (int i = 0; i < 5; ++i) dq[i] = yd[(size_t)ic * 5 + i];
        dq_new_0[ic] = (flow_float)(implicit_relax * dq[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * dq[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * dq[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * dq[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * dq[4]);
    }
    for (geom_int p = e - 2; p >= b; --p) {
        const geom_int ic = line_cells[p];
        double nx[5];
        for (int i = 0; i < 5; ++i) {
            double acc = yd[(size_t)ic * 5 + i];
            for (int m = 0; m < 5; ++m) acc += Wd[(size_t)ic * 25 + i * 5 + m] * dq[m];
            nx[i] = acc;
        }
        for (int i = 0; i < 5; ++i) dq[i] = nx[i];
        dq_new_0[ic] = (flow_float)(implicit_relax * dq[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * dq[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * dq[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * dq[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * dq[4]);
        if (p == b) break;   // geom_int が unsigned の場合の負回りガード
    }
}

// v3 (plan time_integration-line-implicit-speed §4.3 案 B): ライン内の並列化。1 ライン = 8 レーンの tile、レーン 0〜4 が行 0〜4 を担当し、
// レーン 5〜7 は shfl に加わるだけ (書かない)。演算の順序は v2 (lineThomasFactor_d / lineThomasSolve_d) と要素ごとに同じにしてある
// (積和の順序・ピボットの選び方・除算の位置)。並列にしたのは行どうしで独立な部分 (行ごとの積和・消去) と、メモリの読み書きの同時発行。
// LU の代入 (前進・後退) は行どうしが直列なので、各レーンが同じ 5×5 の代入を重複して行い、全レーンが同じ解を持つ。
namespace line_implicit_par {
constexpr int TILE = 8;
}

// レーンの行に対応するポインタを 1 回だけ選ぶ (実行時の添字で引くポインタの配列はローカルメモリに置かれて遅い)
template<typename P>
__device__ __forceinline__ P linePick5(int r, P a0, P a1, P a2, P a3, P a4) {
    return (r == 0) ? a0 : (r == 1) ? a1 : (r == 2) ? a2 : (r == 3) ? a3 : a4;
}

// ---- Thomas の並べ替え版 (plan time_integration-line-implicit-speed §5.1 #10・§6.7、FORGE_LINE_LAYOUT=1 の opt-in) ----
// ncu (§6.6): 既定の Thomas は隣のスレッド (隣のライン) が 121 節点離れた番地を読み、1 回の読み込み命令が約 31 セクタに散る。
// Thomas の中だけで使う W・LU・ピボット・y と、代入で毎 sweep 読む Kprev の写しを (ライン内の位置 k, 成分 e, ライン l) の並び
//   idx = (k·B + e)·nLines + l   (B = 25 / 5)
// に置き、同じ k・e を読む隣のスレッドが隣の番地を読むようにする。D・Knext (分解で 1 回)・rhs・dq (毎 sweep 5 個ずつ) は節点番号の並びのまま。
// 演算 (読む値・積和の順序・ピボットの規則) は既定の lineThomasFactor_d<false> / lineThomasSolve_d<false> と要素ごとに同じなので、解はビットで一致するはず。
__device__ __forceinline__ size_t lineLIdx(geom_int k, int e, int B, geom_int nLines, geom_int l) {
    return ((size_t)k * B + e) * (size_t)nLines + (size_t)l;
}
__global__ void lineThomasFactorL_d
(
 geom_int nLines, const geom_int* line_offsets, const geom_int* line_cells,
 const flow_float* Kprev, const flow_float* Knext,
 const flow_float* d00, const flow_float* d01, const flow_float* d02, const flow_float* d03, const flow_float* d04,
 const flow_float* d10, const flow_float* d11, const flow_float* d12, const flow_float* d13, const flow_float* d14,
 const flow_float* d20, const flow_float* d21, const flow_float* d22, const flow_float* d23, const flow_float* d24,
 const flow_float* d30, const flow_float* d31, const flow_float* d32, const flow_float* d33, const flow_float* d34,
 const flow_float* d40, const flow_float* d41, const flow_float* d42, const flow_float* d43, const flow_float* d44,
 double* Wt, double* LUt, signed char* pivt, double* Kpt, unsigned char* faild
)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    const geom_int b = line_offsets[l];
    const geom_int e = line_offsets[l + 1];
    faild[l] = 0;
    for (geom_int p = b; p < e; ++p) {
        const geom_int k = p - b;
        const geom_int ic = line_cells[p];
        double M[5][5] = {
            {(double)d00[ic],(double)d01[ic],(double)d02[ic],(double)d03[ic],(double)d04[ic]},
            {(double)d10[ic],(double)d11[ic],(double)d12[ic],(double)d13[ic],(double)d14[ic]},
            {(double)d20[ic],(double)d21[ic],(double)d22[ic],(double)d23[ic],(double)d24[ic]},
            {(double)d30[ic],(double)d31[ic],(double)d32[ic],(double)d33[ic],(double)d34[ic]},
            {(double)d40[ic],(double)d41[ic],(double)d42[ic],(double)d43[ic],(double)d44[ic]}};
        if (p > b) {
            double Kp[5][5];
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < 5; ++j) {
                    Kp[i][j] = (double)Kprev[(size_t)ic * 25 + i * 5 + j];
                    Kpt[lineLIdx(k, i * 5 + j, 25, nLines, l)] = Kp[i][j];          // 代入で毎 sweep 読む写し
                }
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < 5; ++j) {
                    double macc = 0.0;
                    for (int m = 0; m < 5; ++m) macc += Kp[i][m] * Wt[lineLIdx(k - 1, m * 5 + j, 25, nLines, l)];
                    M[i][j] -= macc;
                }
        }
        int piv[5];
        if (!line_implicit::lu5_factor(M, piv)) { faild[l] = 1; return; }
        for (int i = 0; i < 5; ++i) {
            pivt[lineLIdx(k, i, 5, nLines, l)] = (signed char)piv[i];
            for (int j = 0; j < 5; ++j) LUt[lineLIdx(k, i * 5 + j, 25, nLines, l)] = M[i][j];
        }
        if (p + 1 < e) {
            for (int j = 0; j < 5; ++j) {
                double col[5];
                for (int i = 0; i < 5; ++i) col[i] = (double)Knext[(size_t)ic * 25 + i * 5 + j];
                line_implicit::lu5_solve(M, piv, col);
                for (int i = 0; i < 5; ++i) Wt[lineLIdx(k, i * 5 + j, 25, nLines, l)] = col[i];
            }
        }
    }
}

__global__ void lineThomasSolveL_d
(
 geom_int nLines, const geom_int* line_offsets, const geom_int* line_cells,
 const double* Kpt, const double* Wt, const double* LUt, const signed char* pivt, const unsigned char* faild,
 const flow_float* rhs0, const flow_float* rhs1, const flow_float* rhs2, const flow_float* rhs3, const flow_float* rhs4,
 const flow_float* dq_old_0, const flow_float* dq_old_1, const flow_float* dq_old_2, const flow_float* dq_old_3, const flow_float* dq_old_4,
 flow_float* dq_new_0, flow_float* dq_new_1, flow_float* dq_new_2, flow_float* dq_new_3, flow_float* dq_new_4,
 flow_float implicit_relax,
 double* yt
)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    const geom_int b = line_offsets[l];
    const geom_int e = line_offsets[l + 1];
    if (faild[l] != 0) {
        for (geom_int p = b; p < e; ++p) {
            const geom_int ic = line_cells[p];
            dq_new_0[ic] = dq_old_0[ic]; dq_new_1[ic] = dq_old_1[ic]; dq_new_2[ic] = dq_old_2[ic];
            dq_new_3[ic] = dq_old_3[ic]; dq_new_4[ic] = dq_old_4[ic];
        }
        return;
    }
    for (geom_int p = b; p < e; ++p) {
        const geom_int k = p - b;
        const geom_int ic = line_cells[p];
        double bk[5] = {(double)rhs0[ic],(double)rhs1[ic],(double)rhs2[ic],(double)rhs3[ic],(double)rhs4[ic]};
        if (p > b) {
            for (int i = 0; i < 5; ++i) {
                double bacc = 0.0;
                for (int m = 0; m < 5; ++m)
                    bacc += Kpt[lineLIdx(k, i * 5 + m, 25, nLines, l)] * yt[lineLIdx(k - 1, m, 5, nLines, l)];
                bk[i] += bacc;
            }
        }
        double M[5][5];
        int piv[5];
        for (int i = 0; i < 5; ++i) {
            piv[i] = (int)pivt[lineLIdx(k, i, 5, nLines, l)];
            for (int j = 0; j < 5; ++j) M[i][j] = LUt[lineLIdx(k, i * 5 + j, 25, nLines, l)];
        }
        line_implicit::lu5_solve(M, piv, bk);
        for (int i = 0; i < 5; ++i) yt[lineLIdx(k, i, 5, nLines, l)] = bk[i];
    }
    double dq[5];
    {
        const geom_int k = e - 1 - b;
        const geom_int ic = line_cells[e - 1];
        for (int i = 0; i < 5; ++i) dq[i] = yt[lineLIdx(k, i, 5, nLines, l)];
        dq_new_0[ic] = (flow_float)(implicit_relax * dq[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * dq[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * dq[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * dq[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * dq[4]);
    }
    for (geom_int p = e - 2; p >= b; --p) {
        const geom_int k = p - b;
        const geom_int ic = line_cells[p];
        double nx[5];
        for (int i = 0; i < 5; ++i) {
            double acc = yt[lineLIdx(k, i, 5, nLines, l)];
            for (int m = 0; m < 5; ++m) acc += Wt[lineLIdx(k, i * 5 + m, 25, nLines, l)] * dq[m];
            nx[i] = acc;
        }
        for (int i = 0; i < 5; ++i) dq[i] = nx[i];
        dq_new_0[ic] = (flow_float)(implicit_relax * dq[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * dq[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * dq[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * dq[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * dq[4]);
        if (p == b) break;
    }
}

// ---- 並べ替え版 + 連鎖の短縮 (plan time_integration-line-implicit-speed §5.1 #12、FORGE_LINE_LAYOUT=2 の opt-in) ----
// lineThomasFactorL_d / lineThomasSolveL_d と同じ並び・同じ演算 (読む値・積和の順序・ピボットの規則) で、次の 3 点だけを変える:
//   (1) 前の節点の W (分解) と y (代入の前進) を書いてから読み直さず、レジスタに持ち回す (値は書いたものと同じ)
//   (2) ポインタに __restrict__ を付け、別名の可能性で次の読みが書きの後ろへ下がらないようにする
//   (3) 節点のループを 2 回ずつ展開し、次の節点の読みを早く出せるようにする
// 補正はビットで一致するはず (§6.7 と同じ照合で確かめる)。
__global__ void lineThomasFactorLP_d
(
 geom_int nLines, const geom_int* __restrict__ line_offsets, const geom_int* __restrict__ line_cells,
 const flow_float* __restrict__ Kprev, const flow_float* __restrict__ Knext,
 const flow_float* __restrict__ d00, const flow_float* __restrict__ d01, const flow_float* __restrict__ d02, const flow_float* __restrict__ d03, const flow_float* __restrict__ d04,
 const flow_float* __restrict__ d10, const flow_float* __restrict__ d11, const flow_float* __restrict__ d12, const flow_float* __restrict__ d13, const flow_float* __restrict__ d14,
 const flow_float* __restrict__ d20, const flow_float* __restrict__ d21, const flow_float* __restrict__ d22, const flow_float* __restrict__ d23, const flow_float* __restrict__ d24,
 const flow_float* __restrict__ d30, const flow_float* __restrict__ d31, const flow_float* __restrict__ d32, const flow_float* __restrict__ d33, const flow_float* __restrict__ d34,
 const flow_float* __restrict__ d40, const flow_float* __restrict__ d41, const flow_float* __restrict__ d42, const flow_float* __restrict__ d43, const flow_float* __restrict__ d44,
 double* __restrict__ Wt, double* __restrict__ LUt, signed char* __restrict__ pivt, double* __restrict__ Kpt, unsigned char* __restrict__ faild
)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    const geom_int b = line_offsets[l];
    const geom_int e = line_offsets[l + 1];
    faild[l] = 0;
    double Wp[5][5];                                       // 前の節点の W (レジスタに持ち回す)
    #pragma unroll 2
    for (geom_int p = b; p < e; ++p) {
        const geom_int k = p - b;
        const geom_int ic = line_cells[p];
        double M[5][5] = {
            {(double)d00[ic],(double)d01[ic],(double)d02[ic],(double)d03[ic],(double)d04[ic]},
            {(double)d10[ic],(double)d11[ic],(double)d12[ic],(double)d13[ic],(double)d14[ic]},
            {(double)d20[ic],(double)d21[ic],(double)d22[ic],(double)d23[ic],(double)d24[ic]},
            {(double)d30[ic],(double)d31[ic],(double)d32[ic],(double)d33[ic],(double)d34[ic]},
            {(double)d40[ic],(double)d41[ic],(double)d42[ic],(double)d43[ic],(double)d44[ic]}};
        if (p > b) {
            double Kp[5][5];
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < 5; ++j) {
                    Kp[i][j] = (double)Kprev[(size_t)ic * 25 + i * 5 + j];
                    Kpt[lineLIdx(k, i * 5 + j, 25, nLines, l)] = Kp[i][j];
                }
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < 5; ++j) {
                    double macc = 0.0;
                    for (int m = 0; m < 5; ++m) macc += Kp[i][m] * Wp[m][j];
                    M[i][j] -= macc;
                }
        }
        int piv[5];
        if (!line_implicit::lu5_factor(M, piv)) { faild[l] = 1; return; }
        for (int i = 0; i < 5; ++i) {
            pivt[lineLIdx(k, i, 5, nLines, l)] = (signed char)piv[i];
            for (int j = 0; j < 5; ++j) LUt[lineLIdx(k, i * 5 + j, 25, nLines, l)] = M[i][j];
        }
        if (p + 1 < e) {
            for (int j = 0; j < 5; ++j) {
                double col[5];
                for (int i = 0; i < 5; ++i) col[i] = (double)Knext[(size_t)ic * 25 + i * 5 + j];
                line_implicit::lu5_solve(M, piv, col);
                for (int i = 0; i < 5; ++i) { Wt[lineLIdx(k, i * 5 + j, 25, nLines, l)] = col[i]; Wp[i][j] = col[i]; }
            }
        }
    }
}

__global__ void lineThomasSolveLP_d
(
 geom_int nLines, const geom_int* __restrict__ line_offsets, const geom_int* __restrict__ line_cells,
 const double* __restrict__ Kpt, const double* __restrict__ Wt, const double* __restrict__ LUt, const signed char* __restrict__ pivt, const unsigned char* __restrict__ faild,
 const flow_float* __restrict__ rhs0, const flow_float* __restrict__ rhs1, const flow_float* __restrict__ rhs2, const flow_float* __restrict__ rhs3, const flow_float* __restrict__ rhs4,
 const flow_float* __restrict__ dq_old_0, const flow_float* __restrict__ dq_old_1, const flow_float* __restrict__ dq_old_2, const flow_float* __restrict__ dq_old_3, const flow_float* __restrict__ dq_old_4,
 flow_float* __restrict__ dq_new_0, flow_float* __restrict__ dq_new_1, flow_float* __restrict__ dq_new_2, flow_float* __restrict__ dq_new_3, flow_float* __restrict__ dq_new_4,
 flow_float implicit_relax,
 double* __restrict__ yt
)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    const geom_int b = line_offsets[l];
    const geom_int e = line_offsets[l + 1];
    if (faild[l] != 0) {
        for (geom_int p = b; p < e; ++p) {
            const geom_int ic = line_cells[p];
            dq_new_0[ic] = dq_old_0[ic]; dq_new_1[ic] = dq_old_1[ic]; dq_new_2[ic] = dq_old_2[ic];
            dq_new_3[ic] = dq_old_3[ic]; dq_new_4[ic] = dq_old_4[ic];
        }
        return;
    }
    double yp[5] = {0.0, 0.0, 0.0, 0.0, 0.0};              // 前の節点の y (レジスタに持ち回す)
    #pragma unroll 2
    for (geom_int p = b; p < e; ++p) {
        const geom_int k = p - b;
        const geom_int ic = line_cells[p];
        double bk[5] = {(double)rhs0[ic],(double)rhs1[ic],(double)rhs2[ic],(double)rhs3[ic],(double)rhs4[ic]};
        if (p > b) {
            for (int i = 0; i < 5; ++i) {
                double bacc = 0.0;
                for (int m = 0; m < 5; ++m) bacc += Kpt[lineLIdx(k, i * 5 + m, 25, nLines, l)] * yp[m];
                bk[i] += bacc;
            }
        }
        double M[5][5];
        int piv[5];
        for (int i = 0; i < 5; ++i) {
            piv[i] = (int)pivt[lineLIdx(k, i, 5, nLines, l)];
            for (int j = 0; j < 5; ++j) M[i][j] = LUt[lineLIdx(k, i * 5 + j, 25, nLines, l)];
        }
        line_implicit::lu5_solve(M, piv, bk);
        for (int i = 0; i < 5; ++i) { yt[lineLIdx(k, i, 5, nLines, l)] = bk[i]; yp[i] = bk[i]; }
    }
    double dq[5];
    {
        const geom_int ic = line_cells[e - 1];
        for (int i = 0; i < 5; ++i) dq[i] = yp[i];
        dq_new_0[ic] = (flow_float)(implicit_relax * dq[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * dq[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * dq[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * dq[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * dq[4]);
    }
    #pragma unroll 2
    for (geom_int p = e - 2; p >= b; --p) {
        const geom_int k = p - b;
        const geom_int ic = line_cells[p];
        double nx[5];
        for (int i = 0; i < 5; ++i) {
            double acc = yt[lineLIdx(k, i, 5, nLines, l)];
            for (int m = 0; m < 5; ++m) acc += Wt[lineLIdx(k, i * 5 + m, 25, nLines, l)] * dq[m];
            nx[i] = acc;
        }
        for (int i = 0; i < 5; ++i) dq[i] = nx[i];
        dq_new_0[ic] = (flow_float)(implicit_relax * dq[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * dq[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * dq[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * dq[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * dq[4]);
        if (p == b) break;
    }
}

// ---- Thomas の float 版 (plan time_integration-line-implicit-speed §5.1 #8・§6.4、FORGE_LINE_F32=1/2 の opt-in 実験) ----
// 既定の経路 (lineThomasFactor_d / lineThomasSolve_d、double) とは別のカーネル。演算の型 T と、Thomas の中だけで使う W・LU・y の保存の型 S を選ぶ:
//   FORGE_LINE_F32=1: T = float、S = double (演算だけ float。律速が演算かメモリかの切り分け用)
//   FORGE_LINE_F32=2: T = float、S = float  (演算と保存の両方 float)
// D・K・rhs・dq の保存 (flow_float) と組立 (ISP 0 なら float) は不変。式・ピボットの規則・失敗の判定 (|pivot| < 1e-30) は double 版と同じ。
namespace line_implicit_t {
template<typename T>
__device__ __forceinline__ bool lu5_factor(T A[5][5], int piv[5])
{
    for (int col = 0; col < 5; ++col) {
        int pv = col; T pa = fabs(A[col][col]);
        for (int r = col + 1; r < 5; ++r) { const T c = fabs(A[r][col]); if (c > pa) { pv = r; pa = c; } }
        if (!(pa >= (T)1.0e-30) || !isfinite(pa)) return false;   // NaN・Inf のピボットも失敗にする (float 版だけ、codex 2026-10-09 M2)
        piv[col] = pv;
        if (pv != col) for (int k = 0; k < 5; ++k) { const T tmp = A[col][k]; A[col][k] = A[pv][k]; A[pv][k] = tmp; }
        const T inv = (T)1.0 / A[col][col];
        for (int r = col + 1; r < 5; ++r) {
            const T f = A[r][col] * inv;
            A[r][col] = f;
            for (int k = col + 1; k < 5; ++k) A[r][k] -= f * A[col][k];
        }
    }
    for (int i = 0; i < 5; ++i) for (int k = 0; k < 5; ++k) if (!isfinite(A[i][k])) return false;   // 消去で溢れた因子も失敗にする
    return true;
}
template<typename T>
__device__ __forceinline__ void lu5_solve(const T A[5][5], const int piv[5], T x[5])
{
    for (int col = 0; col < 5; ++col) if (piv[col] != col) { const T tmp = x[col]; x[col] = x[piv[col]]; x[piv[col]] = tmp; }
    for (int col = 0; col < 5; ++col) for (int r = col + 1; r < 5; ++r) x[r] -= A[r][col] * x[col];
    for (int r = 4; r >= 0; --r) {
        T s = x[r];
        for (int c = r + 1; c < 5; ++c) s -= A[r][c] * x[c];
        x[r] = s / A[r][r];
    }
}
} // namespace line_implicit_t

template<typename T, typename S>
__global__ void lineThomasFactorT_d
(
 geom_int nLines, const geom_int* line_offsets, const geom_int* line_cells,
 const flow_float* Kprev, const flow_float* Knext,
 const flow_float* d00, const flow_float* d01, const flow_float* d02, const flow_float* d03, const flow_float* d04,
 const flow_float* d10, const flow_float* d11, const flow_float* d12, const flow_float* d13, const flow_float* d14,
 const flow_float* d20, const flow_float* d21, const flow_float* d22, const flow_float* d23, const flow_float* d24,
 const flow_float* d30, const flow_float* d31, const flow_float* d32, const flow_float* d33, const flow_float* d34,
 const flow_float* d40, const flow_float* d41, const flow_float* d42, const flow_float* d43, const flow_float* d44,
 S* Wd, S* LUd, signed char* pivd, unsigned char* faild
)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    const geom_int b = line_offsets[l];
    const geom_int e = line_offsets[l + 1];
    faild[l] = 0;
    for (geom_int p = b; p < e; ++p) {
        const geom_int ic = line_cells[p];
        T M[5][5] = {
            {(T)d00[ic],(T)d01[ic],(T)d02[ic],(T)d03[ic],(T)d04[ic]},
            {(T)d10[ic],(T)d11[ic],(T)d12[ic],(T)d13[ic],(T)d14[ic]},
            {(T)d20[ic],(T)d21[ic],(T)d22[ic],(T)d23[ic],(T)d24[ic]},
            {(T)d30[ic],(T)d31[ic],(T)d32[ic],(T)d33[ic],(T)d34[ic]},
            {(T)d40[ic],(T)d41[ic],(T)d42[ic],(T)d43[ic],(T)d44[ic]}};
        if (p > b) {
            const geom_int icm = line_cells[p - 1];
            T Kp[5][5];
            for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) Kp[i][j] = (T)Kprev[(size_t)ic * 25 + i * 5 + j];
            for (int i = 0; i < 5; ++i)
                for (int j = 0; j < 5; ++j) {
                    T macc = (T)0.0;
                    for (int m = 0; m < 5; ++m) macc += Kp[i][m] * (T)Wd[(size_t)icm * 25 + m * 5 + j];
                    M[i][j] -= macc;
                }
        }
        int piv[5];
        if (!line_implicit_t::lu5_factor<T>(M, piv)) { faild[l] = 1; return; }
        for (int i = 0; i < 5; ++i) {
            pivd[(size_t)ic * 5 + i] = (signed char)piv[i];
            for (int j = 0; j < 5; ++j) LUd[(size_t)ic * 25 + i * 5 + j] = (S)M[i][j];
        }
        if (p + 1 < e) {
            for (int j = 0; j < 5; ++j) {
                T col[5];
                for (int i = 0; i < 5; ++i) col[i] = (T)Knext[(size_t)ic * 25 + i * 5 + j];
                line_implicit_t::lu5_solve<T>(M, piv, col);
                for (int i = 0; i < 5; ++i) Wd[(size_t)ic * 25 + i * 5 + j] = (S)col[i];
            }
        }
    }
}

template<typename T, typename S>
__global__ void lineThomasSolveT_d
(
 geom_int nLines, const geom_int* line_offsets, const geom_int* line_cells,
 const flow_float* Kprev,
 const S* Wd, const S* LUd, const signed char* pivd, const unsigned char* faild,
 const flow_float* rhs0, const flow_float* rhs1, const flow_float* rhs2, const flow_float* rhs3, const flow_float* rhs4,
 const flow_float* dq_old_0, const flow_float* dq_old_1, const flow_float* dq_old_2, const flow_float* dq_old_3, const flow_float* dq_old_4,
 flow_float* dq_new_0, flow_float* dq_new_1, flow_float* dq_new_2, flow_float* dq_new_3, flow_float* dq_new_4,
 flow_float implicit_relax,
 S* yd
)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    const geom_int b = line_offsets[l];
    const geom_int e = line_offsets[l + 1];
    if (faild[l] != 0) {
        for (geom_int p = b; p < e; ++p) {
            const geom_int ic = line_cells[p];
            dq_new_0[ic] = dq_old_0[ic]; dq_new_1[ic] = dq_old_1[ic]; dq_new_2[ic] = dq_old_2[ic];
            dq_new_3[ic] = dq_old_3[ic]; dq_new_4[ic] = dq_old_4[ic];
        }
        return;
    }
    for (geom_int p = b; p < e; ++p) {
        const geom_int ic = line_cells[p];
        T bk[5] = {(T)rhs0[ic],(T)rhs1[ic],(T)rhs2[ic],(T)rhs3[ic],(T)rhs4[ic]};
        if (p > b) {
            const geom_int icm = line_cells[p - 1];
            for (int i = 0; i < 5; ++i) {
                T bacc = (T)0.0;
                for (int m = 0; m < 5; ++m) bacc += (T)Kprev[(size_t)ic * 25 + i * 5 + m] * (T)yd[(size_t)icm * 5 + m];
                bk[i] += bacc;
            }
        }
        T M[5][5]; int piv[5];
        for (int i = 0; i < 5; ++i) {
            piv[i] = (int)pivd[(size_t)ic * 5 + i];
            for (int j = 0; j < 5; ++j) M[i][j] = (T)LUd[(size_t)ic * 25 + i * 5 + j];
        }
        line_implicit_t::lu5_solve<T>(M, piv, bk);
        for (int i = 0; i < 5; ++i) yd[(size_t)ic * 5 + i] = (S)bk[i];
    }
    T dq[5];
    {
        const geom_int ic = line_cells[e - 1];
        for (int i = 0; i < 5; ++i) dq[i] = (T)yd[(size_t)ic * 5 + i];
        dq_new_0[ic] = (flow_float)(implicit_relax * (flow_float)dq[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * (flow_float)dq[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * (flow_float)dq[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * (flow_float)dq[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * (flow_float)dq[4]);
    }
    for (geom_int p = e - 2; p >= b; --p) {
        const geom_int ic = line_cells[p];
        T nx[5];
        for (int i = 0; i < 5; ++i) {
            T acc = (T)yd[(size_t)ic * 5 + i];
            for (int m = 0; m < 5; ++m) acc += (T)Wd[(size_t)ic * 25 + i * 5 + m] * dq[m];
            nx[i] = acc;
        }
        for (int i = 0; i < 5; ++i) dq[i] = nx[i];
        dq_new_0[ic] = (flow_float)(implicit_relax * (flow_float)dq[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * (flow_float)dq[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * (flow_float)dq[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * (flow_float)dq[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * (flow_float)dq[4]);
        if (p == b) break;
    }
}

__global__ void lineThomasFactorPar_d
(
 geom_int nLines,
 const geom_int* line_offsets,
 const geom_int* line_cells,
 const flow_float* Kprev, const flow_float* Knext,
 const flow_float* d00, const flow_float* d01, const flow_float* d02, const flow_float* d03, const flow_float* d04,
 const flow_float* d10, const flow_float* d11, const flow_float* d12, const flow_float* d13, const flow_float* d14,
 const flow_float* d20, const flow_float* d21, const flow_float* d22, const flow_float* d23, const flow_float* d24,
 const flow_float* d30, const flow_float* d31, const flow_float* d32, const flow_float* d33, const flow_float* d34,
 const flow_float* d40, const flow_float* d41, const flow_float* d42, const flow_float* d43, const flow_float* d44,
 double* Wd, double* LUd, signed char* pivd, unsigned char* faild
)
{
    namespace cg = cooperative_groups;
    cg::thread_block_tile<line_implicit_par::TILE> tile = cg::tiled_partition<line_implicit_par::TILE>(cg::this_thread_block());
    const geom_int l = (geom_int)((blockDim.x * blockIdx.x + threadIdx.x) / line_implicit_par::TILE);
    if (l >= nLines) return;                       // tile の 8 レーンは同じ l なので tile ごと抜ける
    const int r = (int)tile.thread_rank();
    const bool act = (r < 5);
    const int rr = act ? r : 0;                    // 補助レーンは行 0 の値を持つだけ
    const flow_float* const D0 = linePick5(rr, d00, d10, d20, d30, d40);
    const flow_float* const D1 = linePick5(rr, d01, d11, d21, d31, d41);
    const flow_float* const D2 = linePick5(rr, d02, d12, d22, d32, d42);
    const flow_float* const D3 = linePick5(rr, d03, d13, d23, d33, d43);
    const flow_float* const D4 = linePick5(rr, d04, d14, d24, d34, d44);
    const geom_int b = line_offsets[l];
    const geom_int e = line_offsets[l + 1];
    if (r == 0) faild[l] = 0;
    double Wc[5] = {0.0, 0.0, 0.0, 0.0, 0.0};      // 前の節点の W の列 (レーン j が W[:, j])
    for (geom_int p = b; p < e; ++p) {
        const geom_int ic = line_cells[p];
        double M[5] = {(double)D0[ic], (double)D1[ic], (double)D2[ic], (double)D3[ic], (double)D4[ic]};
        if (p > b) {                               // M_i −= Σ_m Kprev[i][m] W_{k−1}[m][j] (v2 と同じ順序)
            double Kp[5];
            for (int m = 0; m < 5; ++m) Kp[m] = (double)Kprev[(size_t)ic * 25 + rr * 5 + m];
            for (int j = 0; j < 5; ++j) {
                double macc = 0.0;
                for (int m = 0; m < 5; ++m) macc += Kp[m] * tile.shfl(Wc[m], j);
                M[j] -= macc;
            }
        }
        // 部分ピボット付き LU (line_implicit::lu5_factor と同じ演算、行を各レーンが持つ)
        int piv[5];
        bool fail = false;
        for (int col = 0; col < 5; ++col) {
            double a[5];
            for (int q = 0; q < 5; ++q) a[q] = fabs(tile.shfl(M[col], q));
            int pv = col; double pa = a[col];
            for (int q = col + 1; q < 5; ++q) { if (a[q] > pa) { pv = q; pa = a[q]; } }
            if (pa < 1.0e-30) { fail = true; break; }
            piv[col] = pv;
            if (pv != col) {
                for (int k = 0; k < 5; ++k) {
                    const double vc = tile.shfl(M[k], col);
                    const double vp = tile.shfl(M[k], pv);
                    if (r == col) M[k] = vp; else if (r == pv) M[k] = vc;
                }
            }
            double Pc[5];
            for (int k = 0; k < 5; ++k) Pc[k] = tile.shfl(M[k], col);
            const double inv = 1.0 / Pc[col];
            if (act && r > col) {
                const double f = M[col] * inv;
                M[col] = f;
                for (int k = col + 1; k < 5; ++k) M[k] -= f * Pc[k];
            }
        }
        if (fail) { if (r == 0) faild[l] = 1; return; }
        if (act) {
            pivd[(size_t)ic * 5 + r] = (signed char)piv[r];
            for (int k = 0; k < 5; ++k) LUd[(size_t)ic * 25 + r * 5 + k] = M[k];
        }
        if (p + 1 < e) {                           // W_k = M̃⁻¹ Knext_k: レーン j が列 j を代入 (v2 と同じ lu5_solve)
            double Mf[5][5];
            for (int q = 0; q < 5; ++q)
                for (int k = 0; k < 5; ++k) Mf[q][k] = tile.shfl(M[k], q);
            if (act) {
                double col[5];
                for (int q = 0; q < 5; ++q) col[q] = (double)Knext[(size_t)ic * 25 + q * 5 + r];
                line_implicit::lu5_solve(Mf, piv, col);
                for (int q = 0; q < 5; ++q) { Wc[q] = col[q]; Wd[(size_t)ic * 25 + q * 5 + r] = col[q]; }
            }
        }
    }
}

// 行を分担した LU の代入 (line_implicit::lu5_solve と同じ演算・同じ順序): レーン r は U/L の行 r と x_r を持つ。
// ① 行交換を全て先に (LASWP) ② 単位下三角の前進 ③ 上三角の後退 (行 r は x_{r+1..4} を昇順で引いてから U_rr で割る)。
__device__ __forceinline__ double line_lu5_solve_rows(cooperative_groups::thread_block_tile<line_implicit_par::TILE>& tile,
                                                      int r, const double Lr[5], const int piv[5], double x)
{
    for (int col = 0; col < 5; ++col) {
        const int pc = piv[col];
        if (pc != col) {
            const double xc = tile.shfl(x, col);
            const double xp = tile.shfl(x, pc);
            if (r == col) x = xp; else if (r == pc) x = xc;
        }
    }
    for (int col = 0; col < 4; ++col) {
        const double xc = tile.shfl(x, col);
        if (r > col && r < 5) x -= Lr[col] * xc;
    }
    double xs[5];
    for (int row = 4; row >= 0; --row) {
        double s = x;
        if (row < 4) {
            for (int c = row + 1; c < 5; ++c) xs[c] = tile.shfl(x, c);
            if (r == row) for (int c = row + 1; c < 5; ++c) s -= Lr[c] * xs[c];
        }
        if (r == row) x = s / Lr[row];
        // 行 row の x が確定したので、次の行 (row−1) は shfl で読む (shfl は全レーンで実行する)
    }
    return x;
}

__global__ void lineThomasSolvePar_d
(
 geom_int nLines,
 const geom_int* line_offsets,
 const geom_int* line_cells,
 const flow_float* Kprev,
 const double* Wd, const double* LUd, const signed char* pivd, const unsigned char* faild,
 const flow_float* rhs0, const flow_float* rhs1, const flow_float* rhs2, const flow_float* rhs3, const flow_float* rhs4,
 const flow_float* dq_old_0, const flow_float* dq_old_1, const flow_float* dq_old_2, const flow_float* dq_old_3, const flow_float* dq_old_4,
 flow_float* dq_new_0, flow_float* dq_new_1, flow_float* dq_new_2, flow_float* dq_new_3, flow_float* dq_new_4,
 flow_float implicit_relax,
 double* yd
)
{
    namespace cg = cooperative_groups;
    cg::thread_block_tile<line_implicit_par::TILE> tile = cg::tiled_partition<line_implicit_par::TILE>(cg::this_thread_block());
    const geom_int l = (geom_int)((blockDim.x * blockIdx.x + threadIdx.x) / line_implicit_par::TILE);
    if (l >= nLines) return;
    const int r = (int)tile.thread_rank();
    const bool act = (r < 5);
    const int rr = act ? r : 0;
    const flow_float* const RHS = linePick5(rr, rhs0, rhs1, rhs2, rhs3, rhs4);
    const flow_float* const DQO = linePick5(rr, dq_old_0, dq_old_1, dq_old_2, dq_old_3, dq_old_4);
    flow_float* const DQN = linePick5(rr, dq_new_0, dq_new_1, dq_new_2, dq_new_3, dq_new_4);
    const geom_int b = line_offsets[l];
    const geom_int e = line_offsets[l + 1];
    if (faild[l] != 0) {
        if (act) for (geom_int p = b; p < e; ++p) { const geom_int ic = line_cells[p]; DQN[ic] = DQO[ic]; }
        return;
    }
    // ---- 前進 (保存因子で代入のみ) ----
    double y[5] = {0.0, 0.0, 0.0, 0.0, 0.0};      // 前の節点の y (全レーンが全成分を持つ)
    for (geom_int p = b; p < e; ++p) {
        const geom_int ic = line_cells[p];
        double x = (double)RHS[ic];
        if (p > b) {
            double bacc = 0.0;
            for (int m = 0; m < 5; ++m) bacc += (double)Kprev[(size_t)ic * 25 + rr * 5 + m] * y[m];
            x += bacc;
        }
        double Lr[5];
        int piv[5];
        for (int k = 0; k < 5; ++k) Lr[k] = LUd[(size_t)ic * 25 + rr * 5 + k];
        for (int q = 0; q < 5; ++q) piv[q] = (int)pivd[(size_t)ic * 5 + q];
        x = line_lu5_solve_rows(tile, r, Lr, piv, x);
        for (int q = 0; q < 5; ++q) y[q] = tile.shfl(x, q);
        if (act) yd[(size_t)ic * 5 + r] = x;
    }
    // ---- 後退代入 (relax を掛けて dq_new へ) ----
    double dq[5];
    for (int q = 0; q < 5; ++q) dq[q] = y[q];
    {
        const geom_int ic = line_cells[e - 1];
        if (act) DQN[ic] = (flow_float)(implicit_relax * dq[r]);
    }
    for (geom_int p = e - 1; p > b; --p) {
        const geom_int ic = line_cells[p - 1];
        double acc = yd[(size_t)ic * 5 + rr];
        for (int m = 0; m < 5; ++m) acc += Wd[(size_t)ic * 25 + rr * 5 + m] * dq[m];
        for (int q = 0; q < 5; ++q) dq[q] = tile.shfl(acc, q);
        if (act) DQN[ic] = (flow_float)(implicit_relax * dq[r]);
    }
}

// 診断: ライン CV を「保存済み diag/rhs の点解」だけで更新する (K/Thomas 不使用)。
// FORGE_LINE_DEBUG_POINT=1 で有効。格納 (diag/rhs) の正しさと Thomas 本体の切り分け用。
__global__ void lineDebugPoint_d
(
 geom_int nLines, const geom_int* line_offsets, const geom_int* line_cells,
 const flow_float* d00, const flow_float* d01, const flow_float* d02, const flow_float* d03, const flow_float* d04,
 const flow_float* d10, const flow_float* d11, const flow_float* d12, const flow_float* d13, const flow_float* d14,
 const flow_float* d20, const flow_float* d21, const flow_float* d22, const flow_float* d23, const flow_float* d24,
 const flow_float* d30, const flow_float* d31, const flow_float* d32, const flow_float* d33, const flow_float* d34,
 const flow_float* d40, const flow_float* d41, const flow_float* d42, const flow_float* d43, const flow_float* d44,
 const flow_float* rhs0, const flow_float* rhs1, const flow_float* rhs2, const flow_float* rhs3, const flow_float* rhs4,
 flow_float* dq_new_0, flow_float* dq_new_1, flow_float* dq_new_2, flow_float* dq_new_3, flow_float* dq_new_4,
 flow_float implicit_relax
)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    for (geom_int p = line_offsets[l]; p < line_offsets[l + 1]; ++p) {
        const geom_int ic = line_cells[p];

        double M[5][5] = {
            {(double)d00[ic],(double)d01[ic],(double)d02[ic],(double)d03[ic],(double)d04[ic]},
            {(double)d10[ic],(double)d11[ic],(double)d12[ic],(double)d13[ic],(double)d14[ic]},
            {(double)d20[ic],(double)d21[ic],(double)d22[ic],(double)d23[ic],(double)d24[ic]},
            {(double)d30[ic],(double)d31[ic],(double)d32[ic],(double)d33[ic],(double)d34[ic]},
            {(double)d40[ic],(double)d41[ic],(double)d42[ic],(double)d43[ic],(double)d44[ic]}};
        double bk[5] = {(double)rhs0[ic],(double)rhs1[ic],(double)rhs2[ic],(double)rhs3[ic],(double)rhs4[ic]};
        int piv[5];
        if (line_implicit::lu5_factor(M, piv)) {
            line_implicit::lu5_solve(M, piv, bk);
        } else {
            for (int i = 0; i < 5; ++i) bk[i] = 0.0;
        }
        dq_new_0[ic] = (flow_float)(implicit_relax * bk[0]);
        dq_new_1[ic] = (flow_float)(implicit_relax * bk[1]);
        dq_new_2[ic] = (flow_float)(implicit_relax * bk[2]);
        dq_new_3[ic] = (flow_float)(implicit_relax * bk[3]);
        dq_new_4[ic] = (flow_float)(implicit_relax * bk[4]);
    }
}

// v2: モノリシック版へ戻す退避スイッチ (FORGE_LINE_MONO=1)。既定は factor/solve 分離。
static bool lineMonoEnabled() {
    static const bool v = [](){ const char* e = getenv("FORGE_LINE_MONO"); return e && atoi(e) != 0; }();
    return v;
}

// ---- v3 の切り替えと判定の経路 (plan time_integration-line-implicit-speed §4.3 案 B・§6) ----
// 既定は 1 ライン 1 スレッド (lineThomasFactor_d / lineThomasSolve_d)、FORGE_LINE_PAR=1 でライン内の並列化 (lineThomasFactorPar_d / lineThomasSolvePar_d)。
// FORGE_LINE_COMPARE=1: 1 ライン 1 スレッドを別のバッファで、並列版を本来のバッファで同じ入力から解き、因子 (LU・W・ピボット・失敗) と補正 dq の差を出力する
// (判定用で遅い。解は並列版を使う。2 つは独立に書いた実装なので、片方を変えたときの照合に使う)。
// 2026-10-09: 並列版は 1 ライン 1 スレッドより遅かった (case/45、40.5 / 37.8 ms/step vs 35.0、plan §6.0) ので既定は 1 ライン 1 スレッド、
// (lu5 の行の入れ替えを静的な添字にしてレジスタに置く書き直しも代入を 9.85 → 13.90 ms/step に遅くしたので戻した)、
// 並列版は FORGE_LINE_PAR=1 の opt-in。lineSerialEnabled() は「1 ライン 1 スレッドを使う」の意味のまま残す。
static bool lineSerialEnabled() {
    static const bool v = [](){ const char* e = getenv("FORGE_LINE_PAR"); return !(e && atoi(e) != 0); }();
    return v;
}
static bool lineCompareEnabled() {
    static const bool v = [](){ const char* e = getenv("FORGE_LINE_COMPARE"); return e && atoi(e) != 0; }();
    return v;
}
// FORGE_LINE_INV=1 (plan time_integration-line-implicit-speed §5.1 #7、opt-in の実験): 1 ライン 1 スレッドの factor で逆行列を保存し、solve の前進を行列ベクトル積にする。
// 並列版 (FORGE_LINE_PAR=1) とは組み合わせない (起動時に止める)。FORGE_LINE_COMPARE=1 と組むと、従来の LU (別のバッファ) と逆行列 (本来のバッファ) を比べる。
static bool lineInvEnabled() {
    static const bool v = [](){
        const char* e = getenv("FORGE_LINE_INV"); const bool on = e && atoi(e) != 0;
        if (on) {   // 逆行列の経路を通らない診断の経路とは組み合わせない (codex 2026-10-09 m6: 黙って無効になるのを防ぐ)
            for (const char* k : {"FORGE_LINE_PAR", "FORGE_LINE_MONO", "FORGE_LINE_DEBUG_POINT", "FORGE_LINE_NOOP"}) {
                const char* p = getenv(k);
                if (p && atoi(p) != 0) { fprintf(stderr, "[line] FORGE_LINE_INV と %s は組み合わせない — 止める\n", k); exit(EXIT_FAILURE); }
            }
            printf("[line] FORGE_LINE_INV=1: 逆行列を保存して前進の代入を行列ベクトル積にする (opt-in の実験)\n");
        }
        return on; }();
    return v;
}
// FORGE_LINE_F32=1/2 (plan time_integration-line-implicit-speed §5.1 #8・§6.4、opt-in の実験): Thomas を float で (1 = 演算だけ、2 = 演算と W・LU・y の保存)。
// 逆行列・並列版・診断の経路とは組み合わせない (最初の呼び出しで止める)。FORGE_LINE_COMPARE=1 と組むと従来の double の LU (別のバッファ) と比べる。
static int lineF32Mode() {
    static const int v = [](){
        const char* e = getenv("FORGE_LINE_F32"); const int m = e ? atoi(e) : 0;
        if (m < 0 || m > 2) { fprintf(stderr, "[line] FORGE_LINE_F32=%d は 0/1/2 のどれか — 止める\n", m); exit(EXIT_FAILURE); }
        if (m != 0) {
            for (const char* k : {"FORGE_LINE_INV", "FORGE_LINE_PAR", "FORGE_LINE_MONO", "FORGE_LINE_DEBUG_POINT", "FORGE_LINE_NOOP"}) {
                const char* q = getenv(k);
                if (q && atoi(q) != 0) { fprintf(stderr, "[line] FORGE_LINE_F32 と %s は組み合わせない — 止める\n", k); exit(EXIT_FAILURE); }
            }
            printf("[line] FORGE_LINE_F32=%d: Thomas を float で (%s、opt-in の実験)\n", m, m == 1 ? "演算だけ float、W・LU・y は double" : "演算と W・LU・y の保存を float");
        }
        return m; }();
    return v;
}
// Thomas の配列の並び (plan time_integration-line-implicit-speed §5.1 #10・#12・#13、§6.7・§6.10・§6.19)。数値は不変 (ビット一致を確認済み)。
//   0 = 従来 (節点番号の並び)、1 = (位置, 成分, ライン) の並び、2 = 1 + 前の節点の W・y をレジスタに持つ (既定)。
// 実効の並びは、メモリを確保する前に、環境変数・診断のスイッチ・メッシュ (区画の数と被覆) から一度だけ決め (lineLayoutResolve)、
// factor・solve・比較・表示で共有する (codex plan-8 M2)。
//   FORGE_LINE_LAYOUT=0/1/2: 明示。1・2 と診断のスイッチ (逆行列・float・並列版・一体型・点の診断・NOOP) の組み合わせは止める。確保の失敗も止める。
//   未指定: 診断のスイッチがあれば 0。区画 (最長 × 本数) が被覆 CV の 1.5 倍を超えれば 0 (長さのばらつきによる無駄を避ける)。
//           並べ替えた配列の確保に失敗すれば、確保済みの分を解放して 0 (codex plan-8 M3)。それ以外は 2。
//   試験用: FORGE_LINE_LAYOUT_SLOT_RATIO (区画比の上限、既定 1.5)、FORGE_LINE_LAYOUT_FAKE_OOM=1 (確保の失敗を模擬)。
namespace line_layout { static int mode = -1; }   // −1 = 未決定
static bool lineLayoutDiagSwitch(const char** which) {
    for (const char* k : {"FORGE_LINE_INV", "FORGE_LINE_F32", "FORGE_LINE_PAR", "FORGE_LINE_MONO", "FORGE_LINE_DEBUG_POINT", "FORGE_LINE_NOOP"}) {
        const char* q = getenv(k);
        if (q && atoi(q) != 0) { if (which) *which = k; return true; }
    }
    return false;
}
static int lineLayoutMode() { return line_layout::mode < 0 ? 0 : line_layout::mode; }
static bool lineLayoutEnabled() { return lineLayoutMode() != 0; }
namespace line_layout {
static double* W = nullptr; static double* LU = nullptr; static double* Kp = nullptr; static double* y = nullptr; static signed char* piv = nullptr;
static size_t slots = 0;   // maxLen · nLines
static void geometry(mesh& msh, geom_int& maxLen, geom_int& nOn) {
    std::vector<geom_int> off(msh.nImplicitLines + 1);
    gpuErrchk(cudaMemcpy(off.data(), msh.line_offsets_d, sizeof(geom_int) * off.size(), cudaMemcpyDeviceToHost));
    maxLen = 0;
    for (geom_int l = 0; l < msh.nImplicitLines; ++l) maxLen = std::max(maxLen, off[l + 1] - off[l]);
    nOn = off[msh.nImplicitLines];
}
static void release() {
    for (double** q : {&W, &LU, &Kp, &y}) { if (*q) cudaFree(*q); *q = nullptr; }
    if (piv) cudaFree(piv);
    piv = nullptr; slots = 0;
}
// 確保できれば true。失敗したら確保済みの分を解放し、CUDA のエラー状態を消して false (自動選択の退避用)。
static bool tryEnsure(mesh& msh) {
    if (W) return true;
    geom_int maxLen = 0, nOn = 0;
    geometry(msh, maxLen, nOn);
    const size_t s = (size_t)maxLen * (size_t)msh.nImplicitLines;
    const char* fake = getenv("FORGE_LINE_LAYOUT_FAKE_OOM");
    const bool fakeOOM = fake && atoi(fake) != 0;
    bool ok = cudaMalloc((void**)&W, sizeof(double) * 25 * s) == cudaSuccess;
    ok = ok && !fakeOOM;                                               // 試験: 1 つ目を確保した後で失敗したことにする (部分の解放を通す)
    ok = ok && cudaMalloc((void**)&LU, sizeof(double) * 25 * s) == cudaSuccess;
    ok = ok && cudaMalloc((void**)&Kp, sizeof(double) * 25 * s) == cudaSuccess;
    ok = ok && cudaMalloc((void**)&y, sizeof(double) * 5 * s) == cudaSuccess;
    ok = ok && cudaMalloc((void**)&piv, 5 * s) == cudaSuccess;
    if (!ok) { release(); (void)cudaGetLastError(); return false; }
    gpuErrchk(cudaMemset(W, 0, sizeof(double) * 25 * s));  gpuErrchk(cudaMemset(LU, 0, sizeof(double) * 25 * s));
    gpuErrchk(cudaMemset(Kp, 0, sizeof(double) * 25 * s)); gpuErrchk(cudaMemset(y, 0, sizeof(double) * 5 * s));
    gpuErrchk(cudaMemset(piv, 0, 5 * s));
    slots = s;
    printf("[line] Thomas の並べ替えた配列: 最長 %ld 節点 × %ld 本 = %zu 区画 (被覆 %ld CV、%.1f MB)\n", (long)maxLen, (long)msh.nImplicitLines, slots,
           (long)nOn, (double)slots * (25 * 3 + 5) * sizeof(double) / 1.0e6);
    return true;
}
static void ensure(mesh& msh) {
    if (!tryEnsure(msh)) { fprintf(stderr, "[line] Thomas の並べ替えた配列を確保できない (FORGE_LINE_LAYOUT 明示) — 止める\n"); exit(EXIT_FAILURE); }
}
} // namespace line_layout
static void lineLayoutResolve(mesh& msh) {
    if (line_layout::mode >= 0) return;
    const char* which = nullptr;
    const bool diag = lineLayoutDiagSwitch(&which);
    const char* e = getenv("FORGE_LINE_LAYOUT");
    if (e) {
        const int v = atoi(e);
        if (v < 0 || v > 2) { fprintf(stderr, "[line] FORGE_LINE_LAYOUT=%s は 0/1/2 のどれか — 止める\n", e); exit(EXIT_FAILURE); }
        if (v != 0 && diag) { fprintf(stderr, "[line] FORGE_LINE_LAYOUT=%d と %s は組み合わせない — 止める\n", v, which); exit(EXIT_FAILURE); }
        line_layout::mode = v;
        if (v != 0) line_layout::ensure(msh);
        printf("[line] Thomas の配列の並び: %s (FORGE_LINE_LAYOUT=%d、明示)\n", v == 2 ? "LAYOUT2" : v == 1 ? "LAYOUT" : "従来 (節点番号の並び)", v);
        return;
    }
    if (diag) {
        line_layout::mode = 0;
        printf("[line] Thomas の配列の並び: 従来 (診断のスイッチ %s があるので既定の LAYOUT2 を使わない)\n", which);
        return;
    }
    geom_int maxLen = 0, nOn = 0;
    line_layout::geometry(msh, maxLen, nOn);
    const char* rl = getenv("FORGE_LINE_LAYOUT_SLOT_RATIO");
    const double limit = rl ? atof(rl) : 1.5;
    const double ratio = (double)maxLen * (double)msh.nImplicitLines / (double)std::max(nOn, (geom_int)1);
    if (ratio > limit) {
        line_layout::mode = 0;
        printf("[line] Thomas の配列の並び: 従来 (区画 %ld × %ld が被覆 %ld CV の %.2f 倍 > %.2f、ラインの長さのばらつきが大きい)\n",
               (long)maxLen, (long)msh.nImplicitLines, (long)nOn, ratio, limit);
        return;
    }
    if (!line_layout::tryEnsure(msh)) {
        line_layout::mode = 0;
        printf("[line] Thomas の配列の並び: 従来 (並べ替えた配列を確保できない)\n");
        return;
    }
    line_layout::mode = 2;
    printf("[line] Thomas の配列の並び: LAYOUT2 (既定。区画比 %.2f、FORGE_LINE_LAYOUT=0 で従来)\n", ratio);
}
namespace line_layout {
// 比較 (FORGE_LINE_COMPARE=1): 並べ替えた因子 (LU・W・ピボット) を従来の節点番号の並びの因子とビット列で比べる (codex 2026-10-10 m3)
__global__ void cmp_d(geom_int nLines, const geom_int* off, const geom_int* cells, const double* LUt, const double* Wt, const signed char* pivt,
                      const double* LU, const double* W, const signed char* piv, const unsigned char* fail, const unsigned char* failRef,
                      unsigned long long* cnt)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    if (fail[l] != 0 || failRef[l] != 0) { if (fail[l] != failRef[l]) atomicAdd(&cnt[3], 1ULL); return; }
    const geom_int b = off[l], e = off[l + 1];
    for (geom_int p = b; p < e; ++p) {
        const geom_int k = p - b, ic = cells[p];
        for (int q = 0; q < 25; ++q) {
            if (__double_as_longlong(LUt[lineLIdx(k, q, 25, nLines, l)]) != __double_as_longlong(LU[(size_t)ic * 25 + q])) atomicAdd(&cnt[0], 1ULL);
            if (p + 1 < e && __double_as_longlong(Wt[lineLIdx(k, q, 25, nLines, l)]) != __double_as_longlong(W[(size_t)ic * 25 + q])) atomicAdd(&cnt[1], 1ULL);
        }
        for (int q = 0; q < 5; ++q) if (pivt[lineLIdx(k, q, 5, nLines, l)] != piv[(size_t)ic * 5 + q]) atomicAdd(&cnt[2], 1ULL);
    }
}
} // namespace line_layout
namespace line_f32 {
static float* W = nullptr; static float* LU = nullptr; static float* y = nullptr;
static void ensure(mesh& msh) {
    if (W) return;
    const size_t n = msh.nCells_all;
    gpuErrchk(cudaMalloc((void**)&W, sizeof(float) * 25 * n));
    gpuErrchk(cudaMalloc((void**)&LU, sizeof(float) * 25 * n));
    gpuErrchk(cudaMalloc((void**)&y, sizeof(float) * 5 * n));
    gpuErrchk(cudaMemset(W, 0, sizeof(float) * 25 * n));
    gpuErrchk(cudaMemset(LU, 0, sizeof(float) * 25 * n));   // ライン外 (ゴースト等) の要素が未初期化のまま非有限の検査に入らないように
    gpuErrchk(cudaMemset(y, 0, sizeof(float) * 5 * n));
}
} // namespace line_f32
// 実効のモードを記録する (factor・solve の 1 回目と 1000 回ごと)
static void lineModeNote(const char* phase, long& count) {
    ++count;
    if (count == 1 || count % 1000 == 0) {
        const char* m = getenv("FORGE_LINE_MONO"); const char* dp = getenv("FORGE_LINE_DEBUG_POINT"); const char* nn = getenv("FORGE_LINE_NOOP");
        const char* mode = (nn && atoi(nn)) ? "NOOP" : (dp && atoi(dp)) ? "DEBUG_POINT" : (m && atoi(m)) ? "MONO" :
                           lineLayoutEnabled() ? (lineLayoutMode() == 2 ? "LAYOUT2" : "LAYOUT") : lineF32Mode() == 1 ? "F32c" : lineF32Mode() == 2 ? "F32cs" :
                           lineInvEnabled() ? "INV" : (!lineSerialEnabled() ? "PAR" : "LU");
        printf("[line] %s %ld 回目: モード %s%s\n", phase, count, mode, lineCompareEnabled() ? " (比較あり)" : "");
    }
}
namespace line_cmp {
struct Alt { double* W = nullptr; double* LU = nullptr; double* y = nullptr; signed char* piv = nullptr; unsigned char* fail = nullptr;
             flow_float* dq[5] = {nullptr, nullptr, nullptr, nullptr, nullptr}; unsigned long long* red = nullptr; double* eta = nullptr; };
static Alt a;
static int factorCalls = 0, solveCalls = 0;
static void ensure(mesh& msh) {
    if (a.W) return;
    const size_t n = msh.nCells_all;
    gpuErrchk(cudaMalloc((void**)&a.W, sizeof(double) * 25 * n));
    gpuErrchk(cudaMalloc((void**)&a.LU, sizeof(double) * 25 * n));
    gpuErrchk(cudaMalloc((void**)&a.y, sizeof(double) * 5 * n));
    gpuErrchk(cudaMalloc((void**)&a.piv, sizeof(signed char) * 5 * n));
    gpuErrchk(cudaMalloc((void**)&a.fail, sizeof(unsigned char) * std::max(msh.nImplicitLines, (geom_int)1)));
    for (int k = 0; k < 5; ++k) gpuErrchk(cudaMalloc((void**)&a.dq[k], sizeof(flow_float) * n));
    gpuErrchk(cudaMalloc((void**)&a.red, sizeof(unsigned long long) * 4));
    gpuErrchk(cudaMalloc((void**)&a.eta, sizeof(double) * std::max(msh.nImplicitLines, (geom_int)1)));
}
// 最大絶対差・最大絶対値 (非負の double のビット列は大小の順を保つので unsigned long long の atomicMax で取る)・一致しない件数・非有限の件数
// (codex 2026-10-09 M1: NaN は atomicMax・std::max の比較で落ちるので、非有限は別に数えて比較の失格にする)
template<typename T>
__global__ void diff_d(size_t n, const T* x, const T* ref, unsigned long long* red) {
    for (size_t i = blockDim.x * (size_t)blockIdx.x + threadIdx.x; i < n; i += (size_t)blockDim.x * gridDim.x) {
        const double xv = (double)x[i], rv = (double)ref[i];
        if (!isfinite(xv) || !isfinite(rv)) { atomicAdd(&red[3], 1ULL); atomicAdd(&red[2], 1ULL); continue; }
        const double d = fabs(xv - rv);
        const double m = fabs(rv);
        atomicMax(&red[0], (unsigned long long)__double_as_longlong(d));
        atomicMax(&red[1], (unsigned long long)__double_as_longlong(m));
        bool ne = false;                                   // ビット列で比べる (+0 と −0 も区別する、codex 2026-10-10 m3)
        const unsigned char* px = reinterpret_cast<const unsigned char*>(&x[i]);
        const unsigned char* pr = reinterpret_cast<const unsigned char*>(&ref[i]);
        for (size_t b = 0; b < sizeof(T); ++b) ne |= (px[b] != pr[b]);
        if (ne) atomicAdd(&red[2], 1ULL);
    }
}
static unsigned long long nonfiniteTotal = 0;   // 比較の全体で見つけた非有限の件数 (1 件でも比較は失格)
// 実際に使った因子・中間 (float の版は float のバッファ) の非有限を数える (codex 2026-10-09 M2: 比較の差だけでは使っていない double の配列を見てしまう)
template<typename S>
__global__ void nonfinite_d(size_t n, const S* x, unsigned long long* cnt) {
    for (size_t i = blockDim.x * (size_t)blockIdx.x + threadIdx.x; i < n; i += (size_t)blockDim.x * gridDim.x)
        if (!isfinite((double)x[i])) atomicAdd(cnt, 1ULL);
}
template<typename S>
static unsigned long long countNonfinite(size_t n, const S* x) {
    gpuErrchk(cudaMemset(a.red, 0, sizeof(unsigned long long)));
    nonfinite_d<S><<<256, 256>>>(n, x, a.red);
    gpuErrchk(cudaPeekAtLastError());
    unsigned long long h = 0;
    gpuErrchk(cudaMemcpy(&h, a.red, sizeof(h), cudaMemcpyDeviceToHost));
    return h;
}
template<typename T>
static void diff(const char* what, size_t n, const T* x, const T* ref, double& md, double& mx, unsigned long long& nd) {
    gpuErrchk(cudaMemset(a.red, 0, sizeof(unsigned long long) * 4));
    diff_d<T><<<256, 256>>>(n, x, ref, a.red);
    gpuErrchk(cudaPeekAtLastError());
    unsigned long long h[4];
    gpuErrchk(cudaMemcpy(h, a.red, sizeof(h), cudaMemcpyDeviceToHost));
    double d, m; memcpy(&d, &h[0], 8); memcpy(&m, &h[1], 8);
    md = std::max(md, d); mx = std::max(mx, m); nd += h[2];
    if (h[3] != 0) { nonfiniteTotal += h[3]; printf("[lineCompare] 非有限 %llu 件 (%s)\n", h[3], what); }
}
// 全ラインの後退誤差 (plan time_integration-line-implicit-speed §6.2 (2)、codex 2026-10-09 M2): 保存した D (storeLU の sweep の値が残る)・Kprev・Knext・rhs と
// 緩和前の解 x = dq/relax から、尺度 S = diag(ρ_ref, ρ_ref a_ref ×3, ρ_ref a_ref²) で無次元化した η = ‖b̂ − Âx̂‖∞ / (‖Â‖∞‖x̂‖∞ + ‖b̂‖∞) をラインごとに出す。
// 分解に失敗したラインは −1 (評価しない)、非有限は +inf。
struct EtaArgs { const flow_float* d[25]; const flow_float* rhs[5]; const flow_float* dq[5]; };
__global__ void eta_d(geom_int nLines, const geom_int* off, const geom_int* cells, const flow_float* Kprev, const flow_float* Knext,
                      EtaArgs A, const unsigned char* fail, double relax, double sc0, double sc1, double sc4, double* eta)
{
    const geom_int l = blockDim.x * blockIdx.x + threadIdx.x;
    if (l >= nLines) return;
    if (fail[l] != 0) { eta[l] = -1.0; return; }
    const double sc[5] = {sc0, sc1, sc1, sc1, sc4};
    const geom_int b = off[l], e = off[l + 1];
    double mA = 0.0, mX = 0.0, mB = 0.0, mR = 0.0;
    bool bad = false;
    for (geom_int p = b; p < e; ++p) {
        const geom_int ic = cells[p];
        double x[5], xm[5], xp[5];
        for (int j = 0; j < 5; ++j) {
            x[j]  = (double)A.dq[j][ic] / relax;
            xm[j] = (p > b)     ? (double)A.dq[j][cells[p - 1]] / relax : 0.0;
            xp[j] = (p + 1 < e) ? (double)A.dq[j][cells[p + 1]] / relax : 0.0;
        }
        for (int i = 0; i < 5; ++i) {
            double r = (double)A.rhs[i][ic], rowA = 0.0;
            for (int j = 0; j < 5; ++j) {
                const double dij = (double)A.d[i * 5 + j][ic];
                r -= dij * x[j]; rowA += fabs(dij) * sc[j];
                if (p > b)     { const double k = (double)Kprev[(size_t)ic * 25 + i * 5 + j]; r += k * xm[j]; rowA += fabs(k) * sc[j]; }
                if (p + 1 < e) { const double k = (double)Knext[(size_t)ic * 25 + i * 5 + j]; r += k * xp[j]; rowA += fabs(k) * sc[j]; }
            }
            const double rh = fabs(r) / sc[i], bh = fabs((double)A.rhs[i][ic]) / sc[i], xh = fabs(x[i]) / sc[i], ah = rowA / sc[i];
            if (!isfinite(rh) || !isfinite(bh) || !isfinite(xh) || !isfinite(ah)) bad = true;
            mR = fmax(mR, rh); mB = fmax(mB, bh); mX = fmax(mX, xh); mA = fmax(mA, ah);
        }
    }
    const double den = mA * mX + mB;
    eta[l] = bad ? INFINITY : (den > 0.0 ? mR / den : 0.0);
}
static void eta(const char* arm, int call, solverConfig& cfg, mesh& msh, variables& var, flow_float* const dq[5], const unsigned char* fail) {
    EtaArgs A;
    static const char* dn[25] = {"diag_block_00","diag_block_01","diag_block_02","diag_block_03","diag_block_04","diag_block_10","diag_block_11","diag_block_12","diag_block_13","diag_block_14",
                                 "diag_block_20","diag_block_21","diag_block_22","diag_block_23","diag_block_24","diag_block_30","diag_block_31","diag_block_32","diag_block_33","diag_block_34",
                                 "diag_block_40","diag_block_41","diag_block_42","diag_block_43","diag_block_44"};
    for (int k = 0; k < 25; ++k) A.d[k] = var.c_d[dn[k]];
    static const char* rn[5] = {"rhs_block_0","rhs_block_1","rhs_block_2","rhs_block_3","rhs_block_4"};
    for (int k = 0; k < 5; ++k) { A.rhs[k] = var.c_d[rn[k]]; A.dq[k] = dq[k]; }
    const double ro = cfg.limiterRoRef > 0.0 ? cfg.limiterRoRef : 1.0, ar = cfg.limiterARef > 0.0 ? cfg.limiterARef : 1.0;
    const int threads = 64, grid = (int)((msh.nImplicitLines + threads - 1) / threads);
    eta_d<<<grid, threads>>>(msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d, msh.line_Kprev_d, msh.line_Knext_d, A, fail,
                             (double)cfg.implicitRelax, ro, ro * ar, ro * ar * ar, a.eta);
    gpuErrchk(cudaPeekAtLastError());
    std::vector<double> h(msh.nImplicitLines);
    gpuErrchk(cudaMemcpy(h.data(), a.eta, sizeof(double) * h.size(), cudaMemcpyDeviceToHost));
    double mx = 0.0; long arg = -1; size_t nEval = 0, nFail = 0, nOver = 0, nBad = 0;
    for (size_t i = 0; i < h.size(); ++i) {
        if (h[i] < 0.0) { ++nFail; continue; }
        ++nEval;
        if (!std::isfinite(h[i])) { ++nBad; continue; }
        if (h[i] > 1.0e-11) ++nOver;
        if (h[i] > mx) { mx = h[i]; arg = (long)i; }
    }
    printf("[lineEta] solve %d %s: η 最大 %.6e (ライン %ld)、評価 %zu 本、1e-11 超 %zu 本、非有限 %zu 本、分解の失敗 %zu 本 (尺度 ρ_ref %.6g・a_ref %.6g)\n",
           call, arm, mx, arg, nEval, nOver, nBad, nFail, ro, ar);
    if (nBad) nonfiniteTotal += nBad;
}
} // namespace line_cmp

static void launchLineFactor(bool par, bool inv, int f32, mesh& msh, variables& var, double* W, double* LU, signed char* piv, unsigned char* fail, bool layout = false)
{
    const int threads = 64;
    const int grid = par ? (int)(((size_t)msh.nImplicitLines * line_implicit_par::TILE + threads - 1) / threads)
                         : (int)((msh.nImplicitLines + threads - 1) / threads);
    #define FORGE_LINE_FACTOR_ARGS \
        msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d, msh.line_Kprev_d, msh.line_Knext_d, \
        var.c_d["diag_block_00"], var.c_d["diag_block_01"], var.c_d["diag_block_02"], var.c_d["diag_block_03"], var.c_d["diag_block_04"], \
        var.c_d["diag_block_10"], var.c_d["diag_block_11"], var.c_d["diag_block_12"], var.c_d["diag_block_13"], var.c_d["diag_block_14"], \
        var.c_d["diag_block_20"], var.c_d["diag_block_21"], var.c_d["diag_block_22"], var.c_d["diag_block_23"], var.c_d["diag_block_24"], \
        var.c_d["diag_block_30"], var.c_d["diag_block_31"], var.c_d["diag_block_32"], var.c_d["diag_block_33"], var.c_d["diag_block_34"], \
        var.c_d["diag_block_40"], var.c_d["diag_block_41"], var.c_d["diag_block_42"], var.c_d["diag_block_43"], var.c_d["diag_block_44"], \
        W, LU, piv, fail
    if (layout) {
        line_layout::ensure(msh);
        if (lineLayoutMode() == 2) lineThomasFactorLP_d<<<grid, threads>>>(
            msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d, msh.line_Kprev_d, msh.line_Knext_d,
            var.c_d["diag_block_00"], var.c_d["diag_block_01"], var.c_d["diag_block_02"], var.c_d["diag_block_03"], var.c_d["diag_block_04"],
            var.c_d["diag_block_10"], var.c_d["diag_block_11"], var.c_d["diag_block_12"], var.c_d["diag_block_13"], var.c_d["diag_block_14"],
            var.c_d["diag_block_20"], var.c_d["diag_block_21"], var.c_d["diag_block_22"], var.c_d["diag_block_23"], var.c_d["diag_block_24"],
            var.c_d["diag_block_30"], var.c_d["diag_block_31"], var.c_d["diag_block_32"], var.c_d["diag_block_33"], var.c_d["diag_block_34"],
            var.c_d["diag_block_40"], var.c_d["diag_block_41"], var.c_d["diag_block_42"], var.c_d["diag_block_43"], var.c_d["diag_block_44"],
            line_layout::W, line_layout::LU, line_layout::piv, line_layout::Kp, fail);
        else lineThomasFactorL_d<<<grid, threads>>>(
            msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d, msh.line_Kprev_d, msh.line_Knext_d,
            var.c_d["diag_block_00"], var.c_d["diag_block_01"], var.c_d["diag_block_02"], var.c_d["diag_block_03"], var.c_d["diag_block_04"],
            var.c_d["diag_block_10"], var.c_d["diag_block_11"], var.c_d["diag_block_12"], var.c_d["diag_block_13"], var.c_d["diag_block_14"],
            var.c_d["diag_block_20"], var.c_d["diag_block_21"], var.c_d["diag_block_22"], var.c_d["diag_block_23"], var.c_d["diag_block_24"],
            var.c_d["diag_block_30"], var.c_d["diag_block_31"], var.c_d["diag_block_32"], var.c_d["diag_block_33"], var.c_d["diag_block_34"],
            var.c_d["diag_block_40"], var.c_d["diag_block_41"], var.c_d["diag_block_42"], var.c_d["diag_block_43"], var.c_d["diag_block_44"],
            line_layout::W, line_layout::LU, line_layout::piv, line_layout::Kp, fail);
    }
    else if (f32 == 1) lineThomasFactorT_d<float, double><<<grid, threads>>>(FORGE_LINE_FACTOR_ARGS);
    else if (f32 == 2) {
        line_f32::ensure(msh);
        #define FORGE_LINE_FACTOR_ARGS_F \
            msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d, msh.line_Kprev_d, msh.line_Knext_d, \
            var.c_d["diag_block_00"], var.c_d["diag_block_01"], var.c_d["diag_block_02"], var.c_d["diag_block_03"], var.c_d["diag_block_04"], \
            var.c_d["diag_block_10"], var.c_d["diag_block_11"], var.c_d["diag_block_12"], var.c_d["diag_block_13"], var.c_d["diag_block_14"], \
            var.c_d["diag_block_20"], var.c_d["diag_block_21"], var.c_d["diag_block_22"], var.c_d["diag_block_23"], var.c_d["diag_block_24"], \
            var.c_d["diag_block_30"], var.c_d["diag_block_31"], var.c_d["diag_block_32"], var.c_d["diag_block_33"], var.c_d["diag_block_34"], \
            var.c_d["diag_block_40"], var.c_d["diag_block_41"], var.c_d["diag_block_42"], var.c_d["diag_block_43"], var.c_d["diag_block_44"], \
            line_f32::W, line_f32::LU, piv, fail
        lineThomasFactorT_d<float, float><<<grid, threads>>>(FORGE_LINE_FACTOR_ARGS_F);
        #undef FORGE_LINE_FACTOR_ARGS_F
    }
    else if (par) lineThomasFactorPar_d<<<grid, threads>>>(FORGE_LINE_FACTOR_ARGS);
    else if (inv) lineThomasFactor_d<true><<<grid, threads>>>(FORGE_LINE_FACTOR_ARGS);
    else          lineThomasFactor_d<false><<<grid, threads>>>(FORGE_LINE_FACTOR_ARGS);
    #undef FORGE_LINE_FACTOR_ARGS
    gpuErrchk( cudaPeekAtLastError() );
}

static void launchLineSolve(bool par, bool inv, int f32, solverConfig& cfg, mesh& msh, variables& var, double* W, double* LU, signed char* piv,
                            unsigned char* fail, double* y, flow_float* const dqn[5], bool layout = false)
{
    const int threads = 64;
    const int grid = par ? (int)(((size_t)msh.nImplicitLines * line_implicit_par::TILE + threads - 1) / threads)
                         : (int)((msh.nImplicitLines + threads - 1) / threads);
    #define FORGE_LINE_SOLVE_ARGS \
        msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d, msh.line_Kprev_d, W, LU, piv, fail, \
        var.c_d["rhs_block_0"], var.c_d["rhs_block_1"], var.c_d["rhs_block_2"], var.c_d["rhs_block_3"], var.c_d["rhs_block_4"], \
        var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"], var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"], \
        dqn[0], dqn[1], dqn[2], dqn[3], dqn[4], cfg.implicitRelax, y
    if (layout) {
        line_layout::ensure(msh);
        if (lineLayoutMode() == 2) lineThomasSolveLP_d<<<grid, threads>>>(
            msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d, line_layout::Kp, line_layout::W, line_layout::LU, line_layout::piv, fail,
            var.c_d["rhs_block_0"], var.c_d["rhs_block_1"], var.c_d["rhs_block_2"], var.c_d["rhs_block_3"], var.c_d["rhs_block_4"],
            var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"], var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"],
            dqn[0], dqn[1], dqn[2], dqn[3], dqn[4], cfg.implicitRelax, line_layout::y);
        else lineThomasSolveL_d<<<grid, threads>>>(
            msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d, line_layout::Kp, line_layout::W, line_layout::LU, line_layout::piv, fail,
            var.c_d["rhs_block_0"], var.c_d["rhs_block_1"], var.c_d["rhs_block_2"], var.c_d["rhs_block_3"], var.c_d["rhs_block_4"],
            var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"], var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"],
            dqn[0], dqn[1], dqn[2], dqn[3], dqn[4], cfg.implicitRelax, line_layout::y);
    }
    else if (f32 == 1) lineThomasSolveT_d<float, double><<<grid, threads>>>(FORGE_LINE_SOLVE_ARGS);
    else if (f32 == 2) {
        line_f32::ensure(msh);
        lineThomasSolveT_d<float, float><<<grid, threads>>>(
            msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d, msh.line_Kprev_d, line_f32::W, line_f32::LU, piv, fail,
            var.c_d["rhs_block_0"], var.c_d["rhs_block_1"], var.c_d["rhs_block_2"], var.c_d["rhs_block_3"], var.c_d["rhs_block_4"],
            var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"], var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"],
            dqn[0], dqn[1], dqn[2], dqn[3], dqn[4], cfg.implicitRelax, line_f32::y);
    }
    else if (par) lineThomasSolvePar_d<<<grid, threads>>>(FORGE_LINE_SOLVE_ARGS);
    else if (inv) lineThomasSolve_d<true><<<grid, threads>>>(FORGE_LINE_SOLVE_ARGS);
    else          lineThomasSolve_d<false><<<grid, threads>>>(FORGE_LINE_SOLVE_ARGS);
    #undef FORGE_LINE_SOLVE_ARGS
    gpuErrchk( cudaPeekAtLastError() );
}

// factor 位相: storeLU した sweep の直後に 1 回だけ呼ぶ (blockDPLURSolve が管理)。
void lineThomasFactor_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (msh.nImplicitLines <= 0) return;
    static long nFactorCalls = 0;
    (void)lineInvEnabled(); (void)lineF32Mode(); lineLayoutResolve(msh);                   // スイッチの整合と並びを最初の呼び出しで決める (診断の早期 return より前)
    lineModeNote("factor", nFactorCalls);
    if (lineMonoEnabled()) return;   // モノリシック時は毎 sweep の lineThomas_d が全てやる
    static const bool dbgPoint = [](){ const char* e = getenv("FORGE_LINE_DEBUG_POINT"); return e && atoi(e) != 0; }();
    static const bool dbgNoop = [](){ const char* e = getenv("FORGE_LINE_NOOP"); return e && atoi(e) != 0; }();
    if (dbgNoop || dbgPoint) return;
    line_dump::atFactor(cfg, msh, var);
    if (!lineCompareEnabled()) {
        launchLineFactor(!lineSerialEnabled(), lineInvEnabled(), lineF32Mode(), msh, var, msh.line_W_d, msh.line_LU_d, msh.line_piv_d, msh.line_fail_d, lineLayoutEnabled());
        return;
    }
    // 比較の相手: FORGE_LINE_F32=1/2 なら float 版、FORGE_LINE_INV=1 なら逆行列 (1 ライン 1 スレッド)、どちらでもなければ並列版。基準は常に従来の LU (1 ライン 1 スレッド) を別のバッファで
    const bool cmpInv = lineInvEnabled();
    const int cmpF32 = lineF32Mode();
    const bool cmpLay = lineLayoutEnabled();
    // 比較: 書かれない要素 (各ラインの最後の節点の W、失敗したラインの残り) で差が出ないよう、両方のバッファを 0 にしてから解く
    line_cmp::ensure(msh);
    const size_t n = msh.nCells_all;
    gpuErrchk(cudaMemset(msh.line_W_d, 0, sizeof(double) * 25 * n)); gpuErrchk(cudaMemset(line_cmp::a.W, 0, sizeof(double) * 25 * n));
    gpuErrchk(cudaMemset(msh.line_LU_d, 0, sizeof(double) * 25 * n)); gpuErrchk(cudaMemset(line_cmp::a.LU, 0, sizeof(double) * 25 * n));
    gpuErrchk(cudaMemset(msh.line_piv_d, 0, 5 * n)); gpuErrchk(cudaMemset(line_cmp::a.piv, 0, 5 * n));
    gpuErrchk(cudaMemset(msh.line_y_d, 0, sizeof(double) * 5 * n)); gpuErrchk(cudaMemset(line_cmp::a.y, 0, sizeof(double) * 5 * n));
    if (cmpF32 == 2) {
        line_f32::ensure(msh);
        gpuErrchk(cudaMemset(line_f32::W, 0, sizeof(float) * 25 * n)); gpuErrchk(cudaMemset(line_f32::LU, 0, sizeof(float) * 25 * n));
        gpuErrchk(cudaMemset(line_f32::y, 0, sizeof(float) * 5 * n));
    }
    launchLineFactor(false, false, 0, msh, var, line_cmp::a.W, line_cmp::a.LU, line_cmp::a.piv, line_cmp::a.fail);
    launchLineFactor(!cmpInv && cmpF32 == 0 && !cmpLay, cmpInv, cmpF32, msh, var, msh.line_W_d, msh.line_LU_d, msh.line_piv_d, msh.line_fail_d, cmpLay);
    double dLU = 0, mLU = 0, dW = 0, mW = 0, dp = 0, mp = 0, df = 0, mf = 0; unsigned long long nLU = 0, nW = 0, np = 0, nf = 0;
    line_cmp::diff("LU", 25 * n, msh.line_LU_d, (const double*)line_cmp::a.LU, dLU, mLU, nLU);
    line_cmp::diff("W", 25 * n, msh.line_W_d, (const double*)line_cmp::a.W, dW, mW, nW);
    line_cmp::diff("piv", 5 * n, msh.line_piv_d, (const signed char*)line_cmp::a.piv, dp, mp, np);
    line_cmp::diff("fail", (size_t)msh.nImplicitLines, msh.line_fail_d, (const unsigned char*)line_cmp::a.fail, df, mf, nf);
    printf("[lineCompare] factor %d%s: LU 最大差 %.3e / 最大 %.3e (不一致 %llu)、W %.3e / %.3e (不一致 %llu)、ピボットの不一致 %llu、失敗の不一致 %llu\n",
           ++line_cmp::factorCalls, cmpLay ? " (LAYOUT: 因子は並べ替えたバッファにあるので LU・W・ピボットの差は意味がない)" :
           cmpF32 == 2 ? " (F32cs: 因子は float のバッファにあるので LU・W・ピボットの差は意味がない)" : cmpF32 == 1 ? " (F32c: 因子は float で計算して double に保存)" :
           cmpInv ? " (逆行列 vs LU: LU・ピボットは中身が違うので比べない)" : "", dLU, mLU, nLU, dW, mW, nW, np, nf);
    {   // 実際に使った因子の非有限 (腕と従来の両方)
        const unsigned long long aLU = cmpLay ? line_cmp::countNonfinite(25 * line_layout::slots, (const double*)line_layout::LU) :
                                       (cmpF32 == 2) ? line_cmp::countNonfinite(25 * n, line_f32::LU) : line_cmp::countNonfinite(25 * n, (const double*)msh.line_LU_d);
        const unsigned long long aW  = cmpLay ? line_cmp::countNonfinite(25 * line_layout::slots, (const double*)line_layout::W) :
                                       (cmpF32 == 2) ? line_cmp::countNonfinite(25 * n, line_f32::W)  : line_cmp::countNonfinite(25 * n, (const double*)msh.line_W_d);
        const unsigned long long rLU = line_cmp::countNonfinite(25 * n, (const double*)line_cmp::a.LU);
        const unsigned long long rW  = line_cmp::countNonfinite(25 * n, (const double*)line_cmp::a.W);
        printf("[lineNonfinite] factor %d: 腕 LU %llu・W %llu、従来 LU %llu・W %llu\n", line_cmp::factorCalls, aLU, aW, rLU, rW);
    }
    if (cmpLay) {   // 並べ替えた因子と従来の因子のビット列の比較
        unsigned long long* cnt = nullptr; gpuErrchk(cudaMalloc((void**)&cnt, sizeof(unsigned long long) * 4)); gpuErrchk(cudaMemset(cnt, 0, sizeof(unsigned long long) * 4));
        const int th = 64, gr = (int)((msh.nImplicitLines + th - 1) / th);
        line_layout::cmp_d<<<gr, th>>>(msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d, line_layout::LU, line_layout::W, line_layout::piv,
                                       line_cmp::a.LU, line_cmp::a.W, line_cmp::a.piv, msh.line_fail_d, line_cmp::a.fail, cnt);
        gpuErrchk(cudaPeekAtLastError());
        unsigned long long h[4]; gpuErrchk(cudaMemcpy(h, cnt, sizeof(h), cudaMemcpyDeviceToHost)); gpuErrchk(cudaFree(cnt));
        printf("[lineLayoutCmp] factor %d: LU 不一致 %llu・W 不一致 %llu・ピボット 不一致 %llu・失敗 不一致 %llu\n", line_cmp::factorCalls, h[0], h[1], h[2], h[3]);
    }
}

void lineThomas_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (msh.nImplicitLines <= 0) return;
    static long nSolveCalls = 0;
    (void)lineInvEnabled(); (void)lineF32Mode(); lineLayoutResolve(msh);
    lineModeNote("solve", nSolveCalls);
    static const bool dbgPoint = [](){ const char* e = getenv("FORGE_LINE_DEBUG_POINT"); return e && atoi(e) != 0; }();
    static const bool dbgNoop = [](){ const char* e = getenv("FORGE_LINE_NOOP"); return e && atoi(e) != 0; }();
    if (dbgNoop) return;   // 切り分け: ライン CV は dq 据え置き (sweep 内 placeholder のまま)
    if (dbgPoint) {
        const int threads0 = 64;
        const int grid0 = (int)((msh.nImplicitLines + threads0 - 1) / threads0);
        lineDebugPoint_d<<<grid0, threads0>>>(
            msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d,
            var.c_d["diag_block_00"], var.c_d["diag_block_01"], var.c_d["diag_block_02"], var.c_d["diag_block_03"], var.c_d["diag_block_04"],
            var.c_d["diag_block_10"], var.c_d["diag_block_11"], var.c_d["diag_block_12"], var.c_d["diag_block_13"], var.c_d["diag_block_14"],
            var.c_d["diag_block_20"], var.c_d["diag_block_21"], var.c_d["diag_block_22"], var.c_d["diag_block_23"], var.c_d["diag_block_24"],
            var.c_d["diag_block_30"], var.c_d["diag_block_31"], var.c_d["diag_block_32"], var.c_d["diag_block_33"], var.c_d["diag_block_34"],
            var.c_d["diag_block_40"], var.c_d["diag_block_41"], var.c_d["diag_block_42"], var.c_d["diag_block_43"], var.c_d["diag_block_44"],
            var.c_d["rhs_block_0"], var.c_d["rhs_block_1"], var.c_d["rhs_block_2"], var.c_d["rhs_block_3"], var.c_d["rhs_block_4"],
            var.c_d["dq_block_new_0"], var.c_d["dq_block_new_1"], var.c_d["dq_block_new_2"], var.c_d["dq_block_new_3"], var.c_d["dq_block_new_4"],
            cfg.implicitRelax);
        gpuErrchk( cudaPeekAtLastError() );
        return;
    }
    const int threads = 64;
    const int grid = (int)((msh.nImplicitLines + threads - 1) / threads);
    if (!lineMonoEnabled()) {
        // 保存済み LU/piv/W での代入のみ (factor は lineThomasFactor_d_wrapper が実施済み)。既定は 1 ライン 1 スレッド、FORGE_LINE_PAR=1 で並列版。
        flow_float* const dqn[5] = {var.c_d["dq_block_new_0"], var.c_d["dq_block_new_1"], var.c_d["dq_block_new_2"], var.c_d["dq_block_new_3"], var.c_d["dq_block_new_4"]};
        if (!lineCompareEnabled()) {
            launchLineSolve(!lineSerialEnabled(), lineInvEnabled(), lineF32Mode(), cfg, msh, var, msh.line_W_d, msh.line_LU_d, msh.line_piv_d, msh.line_fail_d, msh.line_y_d, dqn, lineLayoutEnabled());
        } else {
            // 比較: ライン外の CV の dq_new (点解) を写してから v2 を別のバッファへ、v3 を本来のバッファへ解き、dq_new の 5 成分を比べる
            line_cmp::ensure(msh);
            const size_t n = msh.nCells_all;
            for (int k = 0; k < 5; ++k) gpuErrchk(cudaMemcpy(line_cmp::a.dq[k], dqn[k], sizeof(flow_float) * n, cudaMemcpyDeviceToDevice));
            const bool cmpInv = lineInvEnabled();
            const int cmpF32 = lineF32Mode();
            const bool cmpLay = lineLayoutEnabled();
            launchLineSolve(false, false, 0, cfg, msh, var, line_cmp::a.W, line_cmp::a.LU, line_cmp::a.piv, line_cmp::a.fail, line_cmp::a.y, line_cmp::a.dq);
            launchLineSolve(!cmpInv && cmpF32 == 0 && !cmpLay, cmpInv, cmpF32, cfg, msh, var, msh.line_W_d, msh.line_LU_d, msh.line_piv_d, msh.line_fail_d, msh.line_y_d, dqn, cmpLay);
            double d = 0, m = 0; unsigned long long nd = 0;
            double dk[5], mk[5];                       // 成分ごと (緩和後の dq_new の最大絶対差と最大絶対値)
            for (int k = 0; k < 5; ++k) {
                double d1 = 0, m1 = 0; unsigned long long n1 = 0;
                line_cmp::diff("dq", n, dqn[k], (const flow_float*)line_cmp::a.dq[k], d1, m1, n1);
                dk[k] = d1; mk[k] = m1; d = std::max(d, d1); m = std::max(m, m1); nd += n1;
            }
            printf("[lineCompare] solve %d: dq 最大差 %.3e / 最大 %.3e (不一致 %llu); 成分ごと 差/最大 = %.3e/%.3e %.3e/%.3e %.3e/%.3e %.3e/%.3e %.3e/%.3e (非有限の累計 %llu)\n",
                   ++line_cmp::solveCalls, d, m, nd, dk[0], mk[0], dk[1], mk[1], dk[2], mk[2], dk[3], mk[3], dk[4], mk[4], line_cmp::nonfiniteTotal);
            {   // 実際に使った前進の中間 y の非有限 (腕と従来の両方)
                const unsigned long long ay = cmpLay ? line_cmp::countNonfinite(5 * line_layout::slots, (const double*)line_layout::y) :
                                              (cmpF32 == 2) ? line_cmp::countNonfinite(5 * n, line_f32::y) : line_cmp::countNonfinite(5 * n, (const double*)msh.line_y_d);
                const unsigned long long ry = line_cmp::countNonfinite(5 * n, (const double*)line_cmp::a.y);
                printf("[lineNonfinite] solve %d: 腕 y %llu、従来 y %llu\n", line_cmp::solveCalls, ay, ry);
            }
            if (cmpInv || cmpF32 != 0 || cmpLay) {     // 後退誤差は比べる版 (本来のバッファ) と従来の LU (別のバッファ) の両方
                line_cmp::eta(cmpLay ? (lineLayoutMode() == 2 ? "LAYOUT2" : "LAYOUT") : cmpF32 == 1 ? "F32c" : cmpF32 == 2 ? "F32cs" : "INV", line_cmp::solveCalls, cfg, msh, var, dqn, msh.line_fail_d);
                line_cmp::eta("LU", line_cmp::solveCalls, cfg, msh, var, line_cmp::a.dq, line_cmp::a.fail);
            }
        }
        line_dump::afterSolve(msh, var);
        return;
    }
    lineThomas_d<<<grid, threads>>>(
        msh.nImplicitLines, msh.line_offsets_d, msh.line_cells_d,
        msh.line_Kprev_d, msh.line_Knext_d,
        var.c_d["diag_block_00"], var.c_d["diag_block_01"], var.c_d["diag_block_02"], var.c_d["diag_block_03"], var.c_d["diag_block_04"],
        var.c_d["diag_block_10"], var.c_d["diag_block_11"], var.c_d["diag_block_12"], var.c_d["diag_block_13"], var.c_d["diag_block_14"],
        var.c_d["diag_block_20"], var.c_d["diag_block_21"], var.c_d["diag_block_22"], var.c_d["diag_block_23"], var.c_d["diag_block_24"],
        var.c_d["diag_block_30"], var.c_d["diag_block_31"], var.c_d["diag_block_32"], var.c_d["diag_block_33"], var.c_d["diag_block_34"],
        var.c_d["diag_block_40"], var.c_d["diag_block_41"], var.c_d["diag_block_42"], var.c_d["diag_block_43"], var.c_d["diag_block_44"],
        var.c_d["rhs_block_0"], var.c_d["rhs_block_1"], var.c_d["rhs_block_2"], var.c_d["rhs_block_3"], var.c_d["rhs_block_4"],
        var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"], var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"],
        var.c_d["dq_block_new_0"], var.c_d["dq_block_new_1"], var.c_d["dq_block_new_2"], var.c_d["dq_block_new_3"], var.c_d["dq_block_new_4"],
        cfg.implicitRelax,
        msh.line_W_d, msh.line_y_d);
    gpuErrchk( cudaPeekAtLastError() );
}
