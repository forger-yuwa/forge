#pragma once

#include "cuda_forge/cudaConfig.cuh"

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

__global__ void calcGradient_1_d
( 
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro  ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* T   ,
 flow_float* roe  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,

 //flow_float* droUxdx , flow_float* droUxdy, flow_float* droUxdz,
 //flow_float* droUydx , flow_float* droUydy, flow_float* droUydz,
 //flow_float* droUzdx , flow_float* droUzdy, flow_float* droUzdz,
 //flow_float* droedx  , flow_float* droedy , flow_float* droedz,

 flow_float* divU
);


// calcGradient_b_d は calcGradient_d.cu 内でのみ定義・呼出される (単一 TU)。
// 宣言はここに置かない (シグネチャは同ファイル内の定義を正とする)。

__global__ void calcGradient_2_d
( 
 // mesh structure
 geom_int nCells,
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro  ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* T   ,
 flow_float* roe  ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* drodx  , flow_float* drody , flow_float* drodz,
 flow_float* dPdx   , flow_float* dPdy  , flow_float* dPdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,

 //flow_float* droUxdx  , flow_float* droUxdy , flow_float* droUxdz,
 //flow_float* droUydx  , flow_float* droUydy , flow_float* droUydz,
 //flow_float* droUzdx  , flow_float* droUzdy , flow_float* droUzdz,
 //flow_float* droedx  , flow_float* droedy , flow_float* droedz,

 flow_float* divU
);

void calcGradient_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

// 原始量の AoS パック (stride 8: ro,Ux,Uy,Uz,P,T,0,0)。calcGradient_d_wrapper (gradLSQ==2) が applyBconds 後に組み、
// LSQ 勾配とリミッタの近傍 gather が 6 配列 (6 セクタ) の代わりに 1 セクタで読む。未構築なら nullptr。
const flow_float* prim_pack_device_ptr();

// gradLSQ=2 の事前計算 LSQ 係数 (cInt) の読み取り専用ビュー (plan gradient-scalar-lsq-unification §4.1)。
// cInt は cell_planes CSR の incidence ごとに 3 個 (境界 incidence は 0、継ぎ目は合併済み)。calcGradient_d_wrapper が
// 初回に構築する。nCells・nInc・msh は構築時の値 (単一メッシュ・単一プロセス前提。利用側が一致を検査する)。
struct LsqCoefView {
    const flow_float* cInt = nullptr;
    geom_int nCells = 0;
    geom_int nInc = 0;
    const mesh* msh = nullptr;
};
LsqCoefView lsq_coef_view();

// 段 ② の診断 (FORGE_DIAG_GEOM_STAGE2_DUMP、plans/active/architecture-float-state-double-geometry.md §6.3) 用の LSQ の係数の作り直し。
// 本番の係数・勾配は触らない。legacy = false は本番と同じ変位 (面ごとの差 ge_x..ge_z = double の値の位置の差を 1 回丸めた値)、
// true は段 ① までの変位 (flow_float の座標 ccx..ccz の差)。node の gradLSQ 1/2 だけ (それ以外は空)。
struct LsqDiagCoef {
    int gradLSQ = 0;
    std::vector<flow_float> cInt;       // gradLSQ=2: [3*nInc] (cell_planes CSR の incidence ごと。境界 incidence は 0)
    std::vector<int> seamGroup;         // gradLSQ=2: [nInc] 継ぎ目の group の root (継ぎ目でない incidence は -1)
    std::vector<int> seamClass;         // gradLSQ=2: [nInc] group 内の同値類の番号 (同上 -1)
    int nDegen = 0;                     // gradLSQ=2: スペクトル打ち切りした節点・group の数
    std::vector<flow_float> M6;         // gradLSQ=1: [6*nCells] Mxx,Mxy,Mxz,Myy,Myz,Mzz の順に nCells ずつ (flow_float に格納した値)
};
void lsqCoefForDiag(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, bool legacy, LsqDiagCoef& out);

// スカラー勾配の LSQ 経路 (mesh.scalarGradient: lsq、node のみ)。
inline bool scalarGradientLsqActive(const solverConfig& cfg)
{
    return cfg.discretization == "node" && cfg.scalarGradient == "lsq";
}

// 多変数 LSQ gather: phi_dev[q] (q < nVar、device 上のポインタ配列) の勾配を gx/gy/gz_dev[q] に書く。
// NS と同じ係数・同じ走査順の差分形。最大 4 変数ずつのチャンクで回す。ic < nCells だけ書く (ghost は触らない)。
// 周期の和→broadcast は呼び出し側 (periodicSeamMergeActive のとき)。
void lsqScalarGradient_d_wrapper(cudaConfig& cuda_cfg, mesh& msh, int nVar,
                                 flow_float** phi_dev, flow_float** gx_dev, flow_float** gy_dev, flow_float** gz_dev);
