#include "cuda_forge/calcStructualVariables_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include "flowFormat.hpp"
#include "iostream"
#include <cstdlib>
#include <cstdio>

__global__ 
void calcStructualVariables_d 
( 
 geom_int nPlanes, geom_int nNormalPlanes,
 geom_int* plane_cells ,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz ,  geom_float* ss,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz ,
 geom_float* ccx ,  geom_float* ccy ,  geom_float* ccz ,
 geom_float* fx  ,  geom_float* dcc ,
 int nodeMode   // 1: node-centered (median-dual) → 内部双対面は fx=0.5 (ノード間中点) に固定
)
{
    geom_int  ip  = blockDim.x*blockIdx.x + threadIdx.x;

    //if (ip < nNormalPlanes) {
    if (ip < nPlanes) { // include ghost cell
        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        geom_float ccx0 = ccx[ic0];
        geom_float ccy0 = ccy[ic0];
        geom_float ccz0 = ccz[ic0];

        geom_float ccx1 = ccx[ic1];
        geom_float ccy1 = ccy[ic1];
        geom_float ccz1 = ccz[ic1];

        geom_float dcx = ccx1 - ccx0;
        geom_float dcy = ccy1 - ccy0;
        geom_float dcz = ccz1 - ccz0;
        geom_float dc  = sqrt( pow(dcx, 2.0) + pow(dcy, 2.0) + pow(dcz, 2.0));

        geom_float dc0px = sx[ip]*(pcx[ip] - ccx0)/ss[ip];
        geom_float dc0py = sy[ip]*(pcy[ip] - ccy0)/ss[ip];
        geom_float dc0pz = sz[ip]*(pcz[ip] - ccz0)/ss[ip];
        geom_float dc0p  = sqrt( pow(dc0px, 2.0) + pow(dc0py, 2.0) + pow(dc0pz, 2.0));

        geom_float dc1px = sx[ip]*(pcx[ip] - ccx1)/ss[ip];
        geom_float dc1py = sy[ip]*(pcy[ip] - ccy1)/ss[ip];
        geom_float dc1pz = sz[ip]*(pcz[ip] - ccz1)/ss[ip];
        geom_float dc1p  = sqrt( pow(dc1px, 2.0) + pow(dc1py, 2.0) + pow(dc1pz, 2.0));

        // 退化ガード: 面上に両セル中心が乗る (node 値位置=ノード座標の境界半割面: ノードとその鏡映ゴースト)
        // と 0/0 になるため中点 fx=0.5 にする。非退化面は従来どおり (ビット不変)。
        fx [ip] = ((dc0p + dc1p) > (geom_float)0.0) ? dc1p/(dc0p + dc1p) : (geom_float)0.5;
        dcc[ip] = dc;

        // node-centered: 内部双対面はノード–ノード中点で値を取る (φ_f=½(φ0+φ1)) のが標準的な
        // median-dual エッジ補間。歪み面で重心射影 fx が 0.5 からずれる (cell では従来どおり幾何 fx)。
        // 境界半割面 (ip>=nNormalPlanes) はゴースト/弱形式扱いのため対象外。
        if (nodeMode == 1 && ip < nNormalPlanes) {
            fx[ip] = (geom_float)0.5;
        }
    }
};

void calcStructualVariables_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh,  variables& v)
{
    // (nodeMidpointFx 撤去 2026-08-16 → 2026-09-21 に固定スキームとして復活。下記)
    // node の内部双対面は fx=0.5 (辺中点) の固定スキーム (2026-09-21, plan discretization-node-face-weight-midpoint)。
    // 幾何 fx は (i) 下の式が法線射影でなく成分ごとの積のノルムなので**回転不変でなく**、(ii) 射影に直しても
    // 高 AR の曲面壁層では弦のたるみ κΔs²/8 ~ d1 のため 0.03–1.00 に散る。「値=ノード座標なら幾何 fx は中点相当」
    // (2026-08-16 の nodeMidpointFx 撤去の前提) は等方セルでしか成り立たない。cell と境界半割面は従来どおり。
    const int nodeMode = (cfg.discretization == "node") ? 1 : 0;
    calcStructualVariables_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>>(
        msh.nPlanes, msh.nNormalPlanes,
        msh.map_plane_cells_d,
        v.p_d["sx"] , v.p_d["sy"] , v.p_d["sz"], v.p_d["ss"],
        v.p_d["pcx"], v.p_d["pcy"], v.p_d["pcz"],
        v.c_d["ccx"], v.c_d["ccy"], v.c_d["ccz"],
        v.p_d["fx"] , v.p_d["dcc"],
        nodeMode
    );

    //for (auto& bc : msh.bconds)
    //{
    //    calcStructualVariables_bp_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>(
    //        bc.iPlanes.size() ,
    //        bc.map_bplane_plane_d,  bc.map_bplane_cell_d,
    //        v.p_d["sx"] , v.p_d["sy"] , v.p_d["sz"], v.p_d["ss"],
    //        v.p_d["pcx"], v.p_d["pcy"], v.p_d["pcz"],
    //        v.c_d["ccx"], v.c_d["ccy"], v.c_d["ccz"],
    //        v.p_d["fx"] , v.p_d["dcc"]
    //    );
    //}

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

};

