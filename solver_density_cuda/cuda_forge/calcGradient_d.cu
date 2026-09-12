#include "calcGradient_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"
#include <vector>
#include <cstdio>

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
 flow_float* Ht  ,

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
 //flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz,

 flow_float* divU
)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;

    __syncthreads();

    //if (ip < nNormalPlanes) {
    if (ip < nPlanes) {
        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        geom_float f   = fx[ip];
        flow_float Uxf = f*Ux[ic0] + (1.0f-f)*Ux[ic1];
        flow_float Uyf = f*Uy[ic0] + (1.0f-f)*Uy[ic1];
        flow_float Uzf = f*Uz[ic0] + (1.0f-f)*Uz[ic1];

        flow_float Pf  = f*P[ic0]  + (1.0f-f)*P[ic1];
        flow_float Tf  = f*T[ic0]  + (1.0f-f)*T[ic1];
        flow_float rof = f*ro[ic0] + (1.0f-f)*ro[ic1];
        flow_float roef= f*roe[ic0]+ (1.0f-f)*roe[ic1];
        flow_float Htf= f*Ht[ic0]+ (1.0f-f)*Ht[ic1];
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];

        flow_float US  = Uxf*sxx +Uyf*syy +Uzf*szz;

//        if (ic0 == 0 or ic1 == 0) {
//            printf("ic0 = %d, ic1 = %d, f = %f, Uxf = %f, Uyf = %f, Uzf = %f\n", ic0, ic1, f, Uxf, Uyf, Uzf);
//            printf("sxx = %f, syy = %f, szz = %f\n", sxx, syy, szz);
//            printf("Pf  = %f\n", Pf);    
//            printf("US  = %f\n", US);    
//        }
//
//        if (ic0 == 1 or ic1 == 1) {
//            printf("ic0 = %d, ic1 = %d, f = %f, Uxf = %f, Uyf = %f, Uzf = %f\n", ic0, ic1, f, Uxf, Uyf, Uzf);
//        }



        atomicAdd(&dUxdx[ic0], sxx*Uxf);
        atomicAdd(&dUxdy[ic0], syy*Uxf);
        atomicAdd(&dUxdz[ic0], szz*Uxf);

        atomicAdd(&dUydx[ic0], sxx*Uyf);
        atomicAdd(&dUydy[ic0], syy*Uyf);
        atomicAdd(&dUydz[ic0], szz*Uyf);

        atomicAdd(&dUzdx[ic0], sxx*Uzf);
        atomicAdd(&dUzdy[ic0], syy*Uzf);
        atomicAdd(&dUzdz[ic0], szz*Uzf);

        atomicAdd(&dTdx[ic0], sxx*Tf);
        atomicAdd(&dTdy[ic0], syy*Tf);
        atomicAdd(&dTdz[ic0], szz*Tf);

        atomicAdd(&dPdx[ic0], sxx*Pf);
        atomicAdd(&dPdy[ic0], syy*Pf);
        atomicAdd(&dPdz[ic0], szz*Pf);

        //atomicAdd(&dHtdx[ic0], sxx*Htf);
        //atomicAdd(&dHtdy[ic0], syy*Htf);
        //atomicAdd(&dHtdz[ic0], szz*Htf);



        atomicAdd(&drodx[ic0], sxx*rof);
        atomicAdd(&drody[ic0], syy*rof);
        atomicAdd(&drodz[ic0], szz*rof);


        //atomicAdd(&droUxdx[ic0], sxx*rof*Uxf);
        //atomicAdd(&droUxdy[ic0], syy*rof*Uxf);
        //atomicAdd(&droUxdz[ic0], szz*rof*Uxf);

        //atomicAdd(&droUydx[ic0], sxx*rof*Uyf);
        //atomicAdd(&droUydy[ic0], syy*rof*Uyf);
        //atomicAdd(&droUydz[ic0], szz*rof*Uyf);

        //atomicAdd(&droUzdx[ic0], sxx*rof*Uzf);
        //atomicAdd(&droUzdy[ic0], syy*rof*Uzf);
        //atomicAdd(&droUzdz[ic0], szz*rof*Uzf);

        //atomicAdd(&droedx[ic0], sxx*roef);
        //atomicAdd(&droedy[ic0], syy*roef);
        //atomicAdd(&droedz[ic0], szz*roef);


        atomicAdd(&divU[ic0], US);


        atomicAdd(&dUxdx[ic1], -sxx*Uxf);
        atomicAdd(&dUxdy[ic1], -syy*Uxf);
        atomicAdd(&dUxdz[ic1], -szz*Uxf);

        atomicAdd(&dUydx[ic1], -sxx*Uyf);
        atomicAdd(&dUydy[ic1], -syy*Uyf);
        atomicAdd(&dUydz[ic1], -szz*Uyf);

        atomicAdd(&dUzdx[ic1], -sxx*Uzf);
        atomicAdd(&dUzdy[ic1], -syy*Uzf);
        atomicAdd(&dUzdz[ic1], -szz*Uzf);

        atomicAdd(&dTdx[ic1], -sxx*Tf);
        atomicAdd(&dTdy[ic1], -syy*Tf);
        atomicAdd(&dTdz[ic1], -szz*Tf);

        atomicAdd(&dPdx[ic1], -sxx*Pf);
        atomicAdd(&dPdy[ic1], -syy*Pf);
        atomicAdd(&dPdz[ic1], -szz*Pf);

        //atomicAdd(&dHtdx[ic1], -sxx*Pf);
        //atomicAdd(&dHtdy[ic1], -syy*Pf);
        //atomicAdd(&dHtdz[ic1], -szz*Pf);

        atomicAdd(&drodx[ic1], -sxx*rof);
        atomicAdd(&drody[ic1], -syy*rof);
        atomicAdd(&drodz[ic1], -szz*rof);

        //atomicAdd(&droUxdx[ic1], -sxx*rof*Uxf);
        //atomicAdd(&droUxdy[ic1], -syy*rof*Uxf);
        //atomicAdd(&droUxdz[ic1], -szz*rof*Uxf);

        //atomicAdd(&droUydx[ic1], -sxx*rof*Uyf);
        //atomicAdd(&droUydy[ic1], -syy*rof*Uyf);
        //atomicAdd(&droUydz[ic1], -szz*rof*Uyf);

        //atomicAdd(&droUzdx[ic1], -sxx*rof*Uzf);
        //atomicAdd(&droUzdy[ic1], -syy*rof*Uzf);
        //atomicAdd(&droUzdz[ic1], -szz*rof*Uzf);

        //atomicAdd(&droedx[ic1], -sxx*roef);
        //atomicAdd(&droedy[ic1], -syy*roef);
        //atomicAdd(&droedz[ic1], -szz*roef);

        atomicAdd(&divU[ic1], -US);
    }
    __syncthreads();
}

// node 境界 Green-Gauss 寄与 (methods/discretization.md §6.2/§7.2.2, 2026-08-11 全面改訂)。
// 全 non-periodic bcond・全 primitive (ro,Ux,Uy,Uz,P,T) について、境界面値を owner ノードの
// **状態値そのもの** (bvar ではない) として加算する。no-slip/等温壁/axis のように状態が既に
// BC 値へピン済みの境界ではビット不変 (状態=bvar 相当)、outlet や SST 断熱壁 T_aw のように
// 状態と bvar が乖離し得る境界でのみ実効的に変わる (旧 outlet 専用の `interiorPT` 特例を
// 一般原則に格上げし、bvar 引数ごと撤去した)。
__global__ void calcGradient_b_d
(
  // mesh structure
 geom_int nb,
 geom_int* bplane_plane,
 geom_int* bplane_cell,
 geom_int* bplane_cell_ghst,

 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

// variables (owner ノードの状態; ib==半割面インデックスだが値は ic で引く)
 flow_float* ro  ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* T   ,
 flow_float* roe  ,
 flow_float* Ht   ,

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
 //flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz,

 flow_float* divU
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

        flow_float rof = ro[ic];
        flow_float Uxf = Ux[ic];
        flow_float Uyf = Uy[ic];
        flow_float Uzf = Uz[ic];
        flow_float ek  = 0.5f*(Uxf*Uxf + Uyf*Uyf + Uzf*Uzf);
        flow_float Pf  = P[ic];
        flow_float Tf  = T[ic];
        flow_float roef= roe[ic];
        flow_float Htf = roef/rof + Pf/rof;

        flow_float US  = Uxf*sxx +Uyf*syy +Uzf*szz;

        atomicAdd(&dUxdx[ic], sxx*Uxf);
        atomicAdd(&dUxdy[ic], syy*Uxf);
        atomicAdd(&dUxdz[ic], szz*Uxf);

        atomicAdd(&dUydx[ic], sxx*Uyf);
        atomicAdd(&dUydy[ic], syy*Uyf);
        atomicAdd(&dUydz[ic], szz*Uyf);

        atomicAdd(&dUzdx[ic], sxx*Uzf);
        atomicAdd(&dUzdy[ic], syy*Uzf);
        atomicAdd(&dUzdz[ic], szz*Uzf);

        atomicAdd(&dTdx[ic], sxx*Tf);
        atomicAdd(&dTdy[ic], syy*Tf);
        atomicAdd(&dTdz[ic], szz*Tf);

        atomicAdd(&dPdx[ic], sxx*Pf);
        atomicAdd(&dPdy[ic], syy*Pf);
        atomicAdd(&dPdz[ic], szz*Pf);

        //atomicAdd(&dHtdx[ic], sxx*Htf);
        //atomicAdd(&dHtdy[ic], syy*Htf);
        //atomicAdd(&dHtdz[ic], szz*Htf);

        atomicAdd(&drodx[ic], sxx*rof);
        atomicAdd(&drody[ic], syy*rof);
        atomicAdd(&drodz[ic], szz*rof);

        //atomicAdd(&droUxdx[ic], sxx*rof*Uxf);
        //atomicAdd(&droUxdy[ic], syy*rof*Uxf);
        //atomicAdd(&droUxdz[ic], szz*rof*Uxf);

        //atomicAdd(&droUydx[ic], sxx*rof*Uyf);
        //atomicAdd(&droUydy[ic], syy*rof*Uyf);
        //atomicAdd(&droUydz[ic], szz*rof*Uyf);

        //atomicAdd(&droUzdx[ic], sxx*rof*Uzf);
        //atomicAdd(&droUzdy[ic], syy*rof*Uzf);
        //atomicAdd(&droUzdz[ic], szz*rof*Uzf);

        //atomicAdd(&droedx[ic], sxx*roef);
        //atomicAdd(&droedy[ic], syy*roef);
        //atomicAdd(&droedz[ic], szz*roef);

        atomicAdd(&divU[ic], US);
    }
    __syncthreads();
}


// K5: atomicAdd を使わない cell-gather 版（face並列+atomic を cell並列に置換）。
// 各セルが自分の面を走査し、面値 φf=f*φ[ic0]+(1-f)*φ[ic1] と外向き符号で寄与を集約。
// 内部面は両側セルが各自で再計算（face値補間は安価）。raw sum を書き、正規化は calcGradient_2_d が行う。
__global__ void calcGradient_cellgather_d
(
 geom_int nCells,
 geom_int* plane_cells,
 geom_int* cell_planes_index, geom_int* cell_planes,
 geom_float* sx, geom_float* sy, geom_float* sz,
 geom_float* fx,
 flow_float* ro, flow_float* Ux, flow_float* Uy, flow_float* Uz, flow_float* P, flow_float* T,
 flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
 flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
 flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
 flow_float* drodx, flow_float* drody, flow_float* drodz,
 flow_float* dPdx,  flow_float* dPdy,  flow_float* dPdz,
 flow_float* dTdx,  flow_float* dTdy,  flow_float* dTdz,
 flow_float* divU,
 geom_int gradPlaneSkipFrom   // node 弱形式: ip>=この値 (=nNormalPlanes) の境界面はスキップし
                              // calcGradient_b_d (owner-state) に委ねる。cell では nPlanes 以上を渡し従来どおり全面処理。
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    flow_float gUxx=0,gUxy=0,gUxz=0, gUyx=0,gUyy=0,gUyz=0, gUzx=0,gUzy=0,gUzz=0;
    flow_float grox=0,groy=0,groz=0, gPx=0,gPy=0,gPz=0, gTx=0,gTy=0,gTz=0, gdiv=0;

    const geom_int st = cell_planes_index[ic];
    const geom_int en = cell_planes_index[ic+1];
    for (geom_int ilp=st; ilp<en; ilp++) {
        geom_int ip  = cell_planes[ilp];
        if (ip >= gradPlaneSkipFrom) continue;   // node: 境界半割面は calcGradient_b_d で処理
        geom_int ic0 = plane_cells[2*ip+0];
        geom_int ic1 = plane_cells[2*ip+1];
        flow_float sgn = (ic0==ic) ? 1.0f : -1.0f;   // 外向き法線符号（owner ic0 が +）
        geom_float f = fx[ip];
        flow_float g = 1.0f - f;

        flow_float Uxf=f*Ux[ic0]+g*Ux[ic1];
        flow_float Uyf=f*Uy[ic0]+g*Uy[ic1];
        flow_float Uzf=f*Uz[ic0]+g*Uz[ic1];
        flow_float Pf =f*P[ic0] +g*P[ic1];
        flow_float Tf =f*T[ic0] +g*T[ic1];
        flow_float rof=f*ro[ic0]+g*ro[ic1];

        flow_float sxx=sgn*sx[ip], syy=sgn*sy[ip], szz=sgn*sz[ip];

        gUxx+=sxx*Uxf; gUxy+=syy*Uxf; gUxz+=szz*Uxf;
        gUyx+=sxx*Uyf; gUyy+=syy*Uyf; gUyz+=szz*Uyf;
        gUzx+=sxx*Uzf; gUzy+=syy*Uzf; gUzz+=szz*Uzf;
        gTx +=sxx*Tf;  gTy +=syy*Tf;  gTz +=szz*Tf;
        gPx +=sxx*Pf;  gPy +=syy*Pf;  gPz +=szz*Pf;
        grox+=sxx*rof; groy+=syy*rof; groz+=szz*rof;
        gdiv+= Uxf*sxx + Uyf*syy + Uzf*szz;
    }

    dUxdx[ic]=gUxx; dUxdy[ic]=gUxy; dUxdz[ic]=gUxz;
    dUydx[ic]=gUyx; dUydy[ic]=gUyy; dUydz[ic]=gUyz;
    dUzdx[ic]=gUzx; dUzdy[ic]=gUzy; dUzdz[ic]=gUzz;
    drodx[ic]=grox; drody[ic]=groy; drodz[ic]=groz;
    dPdx[ic]=gPx;   dPdy[ic]=gPy;   dPdz[ic]=gPz;
    dTdx[ic]=gTx;   dTdy[ic]=gTy;   dTdz[ic]=gTz;
    divU[ic]=gdiv;
}

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
 flow_float* Ht  ,

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
 //flow_float* dHtdx  , flow_float* dHtdy , flow_float* dHtdz,

 flow_float* divU
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells ) {

        geom_float volume = vol[ic];

        dUxdx[ic] = dUxdx[ic]/volume;
        dUxdy[ic] = dUxdy[ic]/volume;
        dUxdz[ic] = dUxdz[ic]/volume;
                                    
        dUydx[ic] = dUydx[ic]/volume;
        dUydy[ic] = dUydy[ic]/volume;
        dUydz[ic] = dUydz[ic]/volume;

        dUzdx[ic] = dUzdx[ic]/volume;
        dUzdy[ic] = dUzdy[ic]/volume;
        dUzdz[ic] = dUzdz[ic]/volume;

        dPdx[ic] = dPdx[ic]/volume;
        dPdy[ic] = dPdy[ic]/volume;
        dPdz[ic] = dPdz[ic]/volume;

        drodx[ic] = drodx[ic]/volume;
        drody[ic] = drody[ic]/volume;
        drodz[ic] = drodz[ic]/volume;

        dTdx[ic] = dTdx[ic]/volume;
        dTdy[ic] = dTdy[ic]/volume;
        dTdz[ic] = dTdz[ic]/volume;

        //droUxdx[ic] = droUxdx[ic]/volume;
        //droUxdy[ic] = droUxdy[ic]/volume;
        //droUxdz[ic] = droUxdz[ic]/volume;

        //droUydx[ic] = droUydx[ic]/volume;
        //droUydy[ic] = droUydy[ic]/volume;
        //droUydz[ic] = droUydz[ic]/volume;

        //droUzdx[ic] = droUzdx[ic]/volume;
        //droUzdy[ic] = droUzdy[ic]/volume;
        //droUzdz[ic] = droUzdz[ic]/volume;

        //droedx[ic] = droedx[ic]/volume;
        //droedy[ic] = droedy[ic]/volume;
        //droedz[ic] = droedz[ic]/volume;

        //dHtdx[ic] = dHtdx[ic]/volume;
        //dHtdy[ic] = dHtdy[ic]/volume;
        //dHtdz[ic] = dHtdz[ic]/volume;

        divU[ic] = divU[ic]/volume;
    }
}


// bndFirstOrder 診断/ロバスト化: 境界隣接 CV の再構成勾配 (ro,Ux,Uy,Uz,P) を 0 にし、
// その CV からの 2 次 MUSCL 外挿を 1 次に落とす (壁近傍の過大再構成切り分け/抑制)。
// divU・dT は粘性で使うため残す (本オプションは対流再構成のみを対象とする)。
__global__ void zeroBndNodeGradient_d
(
 geom_int nCells, geom_int* bnode_flag,
 flow_float* drodx, flow_float* drody, flow_float* drodz,
 flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
 flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
 flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
 flow_float* dPdx,  flow_float* dPdy,  flow_float* dPdz
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells || bnode_flag[ic] == 0) return;
    drodx[ic]=0; drody[ic]=0; drodz[ic]=0;
    dUxdx[ic]=0; dUxdy[ic]=0; dUxdz[ic]=0;
    dUydx[ic]=0; dUydy[ic]=0; dUydz[ic]=0;
    dUzdx[ic]=0; dUzdy[ic]=0; dUzdz[ic]=0;
    dPdx[ic]=0;  dPdy[ic]=0;  dPdz[ic]=0;
}

// ===================== 最小二乗 (LSQ) 勾配 (node-centered 用) =====================
// median-dual の近壁で Green-Gauss 面勾配が checkerboard を持ち粘性 BL を振動させるため、
// LSQ 勾配 (近傍 CV 中心の重み付き最小二乗) で checkerboard-free に求める。
// b 配列は既存の勾配配列を流用 (accum で b を貯め、solve で in-place に勾配へ上書き)。
// M (対称 3x3) は scratch。境界は LSQ 点にしない (実在する内部隣接 node のみ。bvar 疑似点は
// M/RHS とも不使用 — methods/discretization.md §7.3, 2026-08-11 改訂)。

// 対称 3x3 を解く。2D (mzz≈0) は xy の 2x2 を解き gz=0。
__device__ __forceinline__ void lsq_solve_sym3(
    double mxx,double mxy,double mxz,double myy,double myz,double mzz,
    double bx,double by,double bz,
    flow_float& gx, flow_float& gy, flow_float& gz)
{
    const double tol = 1.0e-12 * (mxx + myy + 1.0e-300);
    if (mzz <= tol) { // 2D (z 方向に格子が無い)
        double det = mxx*myy - mxy*mxy;
        double idet = (fabs(det) > 1.0e-300) ? 1.0/det : 0.0;
        gx = (flow_float)((myy*bx - mxy*by)*idet);
        gy = (flow_float)((mxx*by - mxy*bx)*idet);
        gz = (flow_float)0.0;
    } else {
        double c00=myy*mzz-myz*myz, c01=mxz*myz-mxy*mzz, c02=mxy*myz-mxz*myy;
        double det=mxx*c00+mxy*c01+mxz*c02;
        double idet=(fabs(det)>1.0e-300)?1.0/det:0.0;
        double c11=mxx*mzz-mxz*mxz, c12=mxz*mxy-mxx*myz, c22=mxx*myy-mxy*mxy;
        gx=(flow_float)((c00*bx+c01*by+c02*bz)*idet);
        gy=(flow_float)((c01*bx+c11*by+c12*bz)*idet);
        gz=(flow_float)((c02*bx+c12*by+c22*bz)*idet);
    }
}

// 内部面 (ip<nNormalPlanes) について per-cell gather で M と b を貯める。
__global__ void lsqGrad_accumInternal_d(
    geom_int nCells, geom_int* plane_cells,
    geom_int* cell_planes_index, geom_int* cell_planes, geom_int nNormalPlanes,
    geom_float* ccx, geom_float* ccy, geom_float* ccz,
    flow_float* ro, flow_float* Ux, flow_float* Uy, flow_float* Uz, flow_float* P, flow_float* T,
    flow_float* Mxx, flow_float* Mxy, flow_float* Mxz, flow_float* Myy, flow_float* Myz, flow_float* Mzz,
    flow_float* drodx, flow_float* drody, flow_float* drodz,
    flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
    flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
    flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
    flow_float* dPdx,  flow_float* dPdy,  flow_float* dPdz,
    flow_float* dTdx,  flow_float* dTdy,  flow_float* dTdz)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    double mxx=0,mxy=0,mxz=0,myy=0,myz=0,mzz=0;
    double bro[3]={0,0,0}, bux[3]={0,0,0}, buy[3]={0,0,0}, buz[3]={0,0,0}, bp[3]={0,0,0}, bt[3]={0,0,0};
    const double cx=ccx[ic], cy=ccy[ic], cz=ccz[ic];
    const double r0=ro[ic], ux0=Ux[ic], uy0=Uy[ic], uz0=Uz[ic], p0=P[ic], t0=T[ic];
    const geom_int st=cell_planes_index[ic], en=cell_planes_index[ic+1];
    for (geom_int ilp=st; ilp<en; ++ilp) {
        const geom_int ip=cell_planes[ilp];
        if (ip>=nNormalPlanes) continue;          // 境界・periodic は LSQ 点にしない (§7.3)
        const geom_int ic0=plane_cells[2*ip+0], ic1=plane_cells[2*ip+1];
        const geom_int jc=(ic0==ic)?ic1:ic0;
        const double dx=ccx[jc]-cx, dy=ccy[jc]-cy, dz=ccz[jc]-cz;
        const double w=1.0/fmax(dx*dx+dy*dy+dz*dz,1.0e-300);
        mxx+=w*dx*dx; mxy+=w*dx*dy; mxz+=w*dx*dz; myy+=w*dy*dy; myz+=w*dy*dz; mzz+=w*dz*dz;
        const double dr=ro[jc]-r0, dux=Ux[jc]-ux0, duy=Uy[jc]-uy0, duz=Uz[jc]-uz0, dp=P[jc]-p0, dt=T[jc]-t0;
        bro[0]+=w*dx*dr; bro[1]+=w*dy*dr; bro[2]+=w*dz*dr;
        bux[0]+=w*dx*dux; bux[1]+=w*dy*dux; bux[2]+=w*dz*dux;
        buy[0]+=w*dx*duy; buy[1]+=w*dy*duy; buy[2]+=w*dz*duy;
        buz[0]+=w*dx*duz; buz[1]+=w*dy*duz; buz[2]+=w*dz*duz;
        bp[0]+=w*dx*dp;  bp[1]+=w*dy*dp;  bp[2]+=w*dz*dp;
        bt[0]+=w*dx*dt;  bt[1]+=w*dy*dt;  bt[2]+=w*dz*dt;
    }
    Mxx[ic]=mxx; Mxy[ic]=mxy; Mxz[ic]=mxz; Myy[ic]=myy; Myz[ic]=myz; Mzz[ic]=mzz;
    drodx[ic]=bro[0]; drody[ic]=bro[1]; drodz[ic]=bro[2];
    dUxdx[ic]=bux[0]; dUxdy[ic]=bux[1]; dUxdz[ic]=bux[2];
    dUydx[ic]=buy[0]; dUydy[ic]=buy[1]; dUydz[ic]=buy[2];
    dUzdx[ic]=buz[0]; dUzdy[ic]=buz[1]; dUzdz[ic]=buz[2];
    dPdx[ic]=bp[0];   dPdy[ic]=bp[1];   dPdz[ic]=bp[2];
    dTdx[ic]=bt[0];   dTdy[ic]=bt[1];   dTdz[ic]=bt[2];
}

// 境界半割面は LSQ 点にしない (methods/discretization.md §7.3, 2026-08-11 撤去)。
// `lsqGrad_accumInternal_d` が既に境界面 (`ip>=nNormalPlanes`) を skip しているため、
// 境界寄与を追加しない = M/b とも境界点ゼロ。旧 `lsqGrad_accumBoundary_d` (bvar 疑似点) は削除した。

// M g = b を変数ごとに解いて勾配を in-place に上書きし、divU を更新する。
__global__ void lsqGrad_solve_d(
    geom_int nCells,
    flow_float* Mxx, flow_float* Mxy, flow_float* Mxz, flow_float* Myy, flow_float* Myz, flow_float* Mzz,
    flow_float* drodx, flow_float* drody, flow_float* drodz,
    flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
    flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
    flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
    flow_float* dPdx,  flow_float* dPdy,  flow_float* dPdz,
    flow_float* dTdx,  flow_float* dTdy,  flow_float* dTdz,
    flow_float* divU)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const double mxx=Mxx[ic], mxy=Mxy[ic], mxz=Mxz[ic], myy=Myy[ic], myz=Myz[ic], mzz=Mzz[ic];
    lsq_solve_sym3(mxx,mxy,mxz,myy,myz,mzz, drodx[ic],drody[ic],drodz[ic], drodx[ic],drody[ic],drodz[ic]);
    lsq_solve_sym3(mxx,mxy,mxz,myy,myz,mzz, dUxdx[ic],dUxdy[ic],dUxdz[ic], dUxdx[ic],dUxdy[ic],dUxdz[ic]);
    lsq_solve_sym3(mxx,mxy,mxz,myy,myz,mzz, dUydx[ic],dUydy[ic],dUydz[ic], dUydx[ic],dUydy[ic],dUydz[ic]);
    lsq_solve_sym3(mxx,mxy,mxz,myy,myz,mzz, dUzdx[ic],dUzdy[ic],dUzdz[ic], dUzdx[ic],dUzdy[ic],dUzdz[ic]);
    lsq_solve_sym3(mxx,mxy,mxz,myy,myz,mzz, dPdx[ic],dPdy[ic],dPdz[ic], dPdx[ic],dPdy[ic],dPdz[ic]);
    lsq_solve_sym3(mxx,mxy,mxz,myy,myz,mzz, dTdx[ic],dTdy[ic],dTdz[ic], dTdx[ic],dTdy[ic],dTdz[ic]);
    divU[ic] = dUxdx[ic] + dUydy[ic] + dUzdz[ic];
}

// ===================== gradLSQ=2: 係数事前計算 LSQ 勾配 (§7.3.1) =====================
// メッシュ静的性を使い、正規方程式の解を幾何のみの線形演算子 g_i = Σ_j c_ij (φ_j−φ_i),
// c_ij = M⁺ w_ij d_ij に畳む。setup 1 回 (double): M 組立 → 解析固有分解 → スペクトル打ち切り
// 擬似逆 M⁺ (λ < gradLSQDegenThresh·λ_max のモードを落とす = 退化方向 1 次化) → float32 係数
// テーブル焼き込み。runtime: float32 gather のみ (M 組立・solve 消滅、係数は 6 変数共有)。
// gradLSQ=1 の「M を float32 格納して毎ステップ solve」の条件数 2 乗増幅と、退化ノード
// (近傍方向が共線/共面) の発散 (case/29, plan §9) を同時に解決する。periodic 境界の半割面は
// LSQ 点に含めない (GG 経路と同方針: DOF 合併機構に委ねる)。
// plans/active/discretization-lsq-gradient.md 候補④ / methods/discretization.md §7.3.1。

// 対称 3×3 の解析固有値 (Smith のトリゴ法)。l0 <= l1 <= l2 で返す。
__device__ __forceinline__ static void lsqPre_sym3_eigvals(
    double m00,double m01,double m02,double m11,double m12,double m22,
    double& l0,double& l1,double& l2)
{
    const double p1 = m01*m01 + m02*m02 + m12*m12;
    const double q  = (m00 + m11 + m22)/3.0;
    const double p2 = (m00-q)*(m00-q) + (m11-q)*(m11-q) + (m22-q)*(m22-q) + 2.0*p1;
    const double p  = sqrt(fmax(p2/6.0, 0.0));
    if (p < 1.0e-300) { l0 = l1 = l2 = q; return; }
    const double b00=(m00-q)/p, b01=m01/p, b02=m02/p, b11=(m11-q)/p, b12=m12/p, b22=(m22-q)/p;
    double detB = b00*(b11*b22-b12*b12) - b01*(b01*b22-b12*b02) + b02*(b01*b12-b11*b02);
    double r = detB/2.0;
    r = fmin(fmax(r, -1.0), 1.0);
    const double phi = acos(r)/3.0;
    const double twoPiOver3 = 2.0943951023931953;
    l2 = q + 2.0*p*cos(phi);
    l0 = q + 2.0*p*cos(phi + twoPiOver3);
    l1 = 3.0*q - l2 - l0;
}

// 固有値 mu に対応する固有ベクトル: (M−μI) の 2 行のクロス積のうち最大ノルムを採用。
// mu が他の固有値からよく分離しているときのみ呼ぶ (呼び出し側で保証)。
__device__ __forceinline__ static void lsqPre_sym3_eigvec(
    double m00,double m01,double m02,double m11,double m12,double m22, double mu,
    double& vx,double& vy,double& vz)
{
    const double r0x=m00-mu, r0y=m01,   r0z=m02;
    const double r1x=m01,    r1y=m11-mu,r1z=m12;
    const double r2x=m02,    r2y=m12,   r2z=m22-mu;
    double c01x=r0y*r1z-r0z*r1y, c01y=r0z*r1x-r0x*r1z, c01z=r0x*r1y-r0y*r1x;
    double c02x=r0y*r2z-r0z*r2y, c02y=r0z*r2x-r0x*r2z, c02z=r0x*r2y-r0y*r2x;
    double c12x=r1y*r2z-r1z*r2y, c12y=r1z*r2x-r1x*r2z, c12z=r1x*r2y-r1y*r2x;
    const double n01=c01x*c01x+c01y*c01y+c01z*c01z;
    const double n02=c02x*c02x+c02y*c02y+c02z*c02z;
    const double n12=c12x*c12x+c12y*c12y+c12z*c12z;
    double nx,ny,nz,nn;
    if (n01 >= n02 && n01 >= n12)      { nx=c01x; ny=c01y; nz=c01z; nn=n01; }
    else if (n02 >= n12)               { nx=c02x; ny=c02y; nz=c02z; nn=n02; }
    else                               { nx=c12x; ny=c12y; nz=c12z; nn=n12; }
    const double inv = (nn > 1.0e-300) ? 1.0/sqrt(nn) : 0.0;
    vx=nx*inv; vy=ny*inv; vz=nz*inv;
}

// 対称 3×3 の直接逆行列 (adjugate/det)。呼び出し側で det の健全性を保証。
__device__ __forceinline__ static void lsqPre_sym3_inv(
    double m00,double m01,double m02,double m11,double m12,double m22,
    double& i00,double& i01,double& i02,double& i11,double& i12,double& i22)
{
    const double c00=m11*m22-m12*m12, c01=m02*m12-m01*m22, c02=m01*m12-m02*m11;
    const double det=m00*c00+m01*c01+m02*c02;
    const double idet=(fabs(det)>1.0e-300)?1.0/det:0.0;
    i00=c00*idet; i01=c01*idet; i02=c02*idet;
    i11=(m00*m22-m02*m02)*idet; i12=(m02*m01-m00*m12)*idet;
    i22=(m00*m11-m01*m01)*idet;
}

// setup 1/3: ノードごとに M (double) を組み、スペクトル打ち切り擬似逆 M⁺ を Minv6 に格納。
// 近傍 = 内部双対面の相手 CV 中心 + 非 periodic 境界半割面の重心 pc (planePeriodic[ip]==0)。
// M は**内部双対面の隣接ノードのみ**で組む (methods/discretization.md §7.3.1, 2026-08-11 改訂:
// 非 periodic 境界半割面重心を LSQ 点として加える旧処理を撤去 — ジャンプが常に 0 の疑似点でも
// M に入れば「その方向の勾配をゼロへ向かわせる」暗黙の正則化になるため)。境界・periodic とも
// `ip>=nNormalPlanes` は一様に skip する (`planePeriodic`/`pcx,pcy,pcz` 引数は不要になった)。
__global__ void lsqPre_setup_d(
    geom_int nCells, geom_int* plane_cells,
    geom_int* cell_planes_index, geom_int* cell_planes, geom_int nNormalPlanes,
    flow_float* ccx, flow_float* ccy, flow_float* ccz,
    double thresh,
    double* Minv6, int* degenCount)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    double m00=0,m01=0,m02=0,m11=0,m12=0,m22=0;
    const double cx=ccx[ic], cy=ccy[ic], cz=ccz[ic];
    const geom_int st=cell_planes_index[ic], en=cell_planes_index[ic+1];
    for (geom_int ilp=st; ilp<en; ++ilp) {
        const geom_int ip=cell_planes[ilp];
        if (ip >= nNormalPlanes) continue;   // 境界・periodic は LSQ 点にしない
        const geom_int ic0=plane_cells[2*ip+0], ic1=plane_cells[2*ip+1];
        const geom_int jc=(ic0==ic)?ic1:ic0;
        const double dx=ccx[jc]-cx, dy=ccy[jc]-cy, dz=ccz[jc]-cz;
        const double w=1.0/fmax(dx*dx+dy*dy+dz*dz,1.0e-300);
        m00+=w*dx*dx; m01+=w*dx*dy; m02+=w*dx*dz; m11+=w*dy*dy; m12+=w*dy*dz; m22+=w*dz*dz;
    }
    double l0,l1,l2;
    lsqPre_sym3_eigvals(m00,m01,m02,m11,m12,m22, l0,l1,l2);
    const double cut = fmax(thresh*l2, 1.0e-300);
    double i00,i01,i02,i11,i12,i22;
    if (l2 <= 1.0e-300) {
        // 近傍ゼロ (来ないはず): 勾配 0
        i00=i01=i02=i11=i12=i22=0.0;
        atomicAdd(degenCount, 1);
    } else if (l0 >= cut) {
        // 全モード健全: 直接逆行列
        lsqPre_sym3_inv(m00,m01,m02,m11,m12,m22, i00,i01,i02,i11,i12,i22);
    } else if (l1 >= cut) {
        // 1 モード退化 (l0 のみ落とす)。l0 は kept (>=cut) からよく分離 → 固有ベクトル健全。
        // M⁺ = inv(M + l2·v0v0ᵀ)·(I − v0v0ᵀ) (シフトで退化方向を持ち上げてから射影で消す)
        double vx,vy,vz;
        lsqPre_sym3_eigvec(m00,m01,m02,m11,m12,m22, l0, vx,vy,vz);
        const double K=l2;
        double s00,s01,s02,s11,s12,s22;
        lsqPre_sym3_inv(m00+K*vx*vx, m01+K*vx*vy, m02+K*vx*vz,
                        m11+K*vy*vy, m12+K*vy*vz, m22+K*vz*vz,
                        s00,s01,s02,s11,s12,s22);
        // P = I − v vᵀ を右から掛ける (両者は同じ固有基底を持ち可換 → 結果は対称)
        const double p00=1.0-vx*vx, p01=-vx*vy, p02=-vx*vz;
        const double p11=1.0-vy*vy, p12=-vy*vz, p22=1.0-vz*vz;
        i00=s00*p00+s01*p01+s02*p02; i01=s00*p01+s01*p11+s02*p12; i02=s00*p02+s01*p12+s02*p22;
        i11=s01*p01+s11*p11+s12*p12; i12=s01*p02+s11*p12+s12*p22;
        i22=s02*p02+s12*p12+s22*p22;
        atomicAdd(degenCount, 1);
    } else {
        // 2 モード退化: 最大固有値方向のみ残す。l2 は dropped (<cut) からよく分離。
        double vx,vy,vz;
        lsqPre_sym3_eigvec(m00,m01,m02,m11,m12,m22, l2, vx,vy,vz);
        const double il=1.0/l2;
        i00=il*vx*vx; i01=il*vx*vy; i02=il*vx*vz;
        i11=il*vy*vy; i12=il*vy*vz; i22=il*vz*vz;
        atomicAdd(degenCount, 1);
    }
    Minv6[6*ic+0]=i00; Minv6[6*ic+1]=i01; Minv6[6*ic+2]=i02;
    Minv6[6*ic+3]=i11; Minv6[6*ic+4]=i12; Minv6[6*ic+5]=i22;
}

// setup 2/3: 内部双対面 incidence の係数 c = M⁺·(w d) を cell_planes CSR 位置に float32 で格納。
__global__ void lsqPre_coefInternal_d(
    geom_int nCells, geom_int* plane_cells,
    geom_int* cell_planes_index, geom_int* cell_planes, geom_int nNormalPlanes,
    flow_float* ccx, flow_float* ccy, flow_float* ccz,
    const double* Minv6, flow_float* cInt)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const double cx=ccx[ic], cy=ccy[ic], cz=ccz[ic];
    const double i00=Minv6[6*ic+0], i01=Minv6[6*ic+1], i02=Minv6[6*ic+2];
    const double i11=Minv6[6*ic+3], i12=Minv6[6*ic+4], i22=Minv6[6*ic+5];
    const geom_int st=cell_planes_index[ic], en=cell_planes_index[ic+1];
    for (geom_int ilp=st; ilp<en; ++ilp) {
        const geom_int ip=cell_planes[ilp];
        if (ip >= nNormalPlanes) { cInt[3*ilp+0]=0; cInt[3*ilp+1]=0; cInt[3*ilp+2]=0; continue; }
        const geom_int ic0=plane_cells[2*ip+0], ic1=plane_cells[2*ip+1];
        const geom_int jc=(ic0==ic)?ic1:ic0;
        const double dx=ccx[jc]-cx, dy=ccy[jc]-cy, dz=ccz[jc]-cz;
        const double w=1.0/fmax(dx*dx+dy*dy+dz*dz,1.0e-300);
        const double bx=w*dx, by=w*dy, bz=w*dz;
        cInt[3*ilp+0]=(flow_float)(i00*bx+i01*by+i02*bz);
        cInt[3*ilp+1]=(flow_float)(i01*bx+i11*by+i12*bz);
        cInt[3*ilp+2]=(flow_float)(i02*bx+i12*by+i22*bz);
    }
}

// runtime 1/2: 内部 incidence の gather (per-node, atomic なし)。6 変数 × 3 成分を直接書く。
__global__ void lsqPreGrad_internal_d(
    geom_int nCells, geom_int* plane_cells,
    geom_int* cell_planes_index, geom_int* cell_planes, geom_int nNormalPlanes,
    const flow_float* cInt,
    flow_float* ro, flow_float* Ux, flow_float* Uy, flow_float* Uz, flow_float* P, flow_float* T,
    flow_float* drodx, flow_float* drody, flow_float* drodz,
    flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
    flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
    flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
    flow_float* dPdx,  flow_float* dPdy,  flow_float* dPdz,
    flow_float* dTdx,  flow_float* dTdy,  flow_float* dTdz)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    flow_float gro[3]={0,0,0}, gux[3]={0,0,0}, guy[3]={0,0,0}, guz[3]={0,0,0}, gp[3]={0,0,0}, gt[3]={0,0,0};
    const flow_float r0=ro[ic], ux0=Ux[ic], uy0=Uy[ic], uz0=Uz[ic], p0=P[ic], t0=T[ic];
    const geom_int st=cell_planes_index[ic], en=cell_planes_index[ic+1];
    for (geom_int ilp=st; ilp<en; ++ilp) {
        const geom_int ip=cell_planes[ilp];
        if (ip >= nNormalPlanes) continue;      // 境界・periodic は LSQ 点にしない (§7.3)
        const geom_int ic0=plane_cells[2*ip+0], ic1=plane_cells[2*ip+1];
        const geom_int jc=(ic0==ic)?ic1:ic0;
        const flow_float c0=cInt[3*ilp+0], c1=cInt[3*ilp+1], c2=cInt[3*ilp+2];
        flow_float d;
        d=ro[jc]-r0;  gro[0]+=c0*d; gro[1]+=c1*d; gro[2]+=c2*d;
        d=Ux[jc]-ux0; gux[0]+=c0*d; gux[1]+=c1*d; gux[2]+=c2*d;
        d=Uy[jc]-uy0; guy[0]+=c0*d; guy[1]+=c1*d; guy[2]+=c2*d;
        d=Uz[jc]-uz0; guz[0]+=c0*d; guz[1]+=c1*d; guz[2]+=c2*d;
        d=P[jc]-p0;   gp[0]+=c0*d;  gp[1]+=c1*d;  gp[2]+=c2*d;
        d=T[jc]-t0;   gt[0]+=c0*d;  gt[1]+=c1*d;  gt[2]+=c2*d;
    }
    drodx[ic]=gro[0]; drody[ic]=gro[1]; drodz[ic]=gro[2];
    dUxdx[ic]=gux[0]; dUxdy[ic]=gux[1]; dUxdz[ic]=gux[2];
    dUydx[ic]=guy[0]; dUydy[ic]=guy[1]; dUydz[ic]=guy[2];
    dUzdx[ic]=guz[0]; dUzdy[ic]=guz[1]; dUzdz[ic]=guz[2];
    dPdx[ic]=gp[0];   dPdy[ic]=gp[1];   dPdz[ic]=gp[2];
    dTdx[ic]=gt[0];   dTdy[ic]=gt[1];   dTdz[ic]=gt[2];
}

// 境界寄与は無い (M に境界点を入れていないため境界の係数テーブル自体が存在しない, §7.3)。
// runtime 2/2: divU。
__global__ void lsqPreGrad_finish_d(
    geom_int nCells,
    flow_float* dUxdx, flow_float* dUydy, flow_float* dUzdz, flow_float* divU)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    divU[ic] = dUxdx[ic] + dUydy[ic] + dUzdz[ic];
}

void calcGradient_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    flow_float* grad_volume = (cfg.isAxisymmetric == 1) ? var.c_d["A_planar"] : var.c_d["volume"];
    flow_float* grad_sx = (cfg.isAxisymmetric == 1) ? var.p_d["sx_planar"] : var.p_d["sx"];
    flow_float* grad_sy = (cfg.isAxisymmetric == 1) ? var.p_d["sy_planar"] : var.p_d["sy"];
    flow_float* grad_sz = (cfg.isAxisymmetric == 1) ? var.p_d["sz_planar"] : var.p_d["sz"];
    flow_float* grad_ss = (cfg.isAxisymmetric == 1) ? var.p_d["ss_planar"] : var.p_d["ss"];

    // initialize
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUxdx"], 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUxdy"], 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUxdz"], 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUydx"], 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUydy"], 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUydz"], 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUzdx"], 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUzdy"], 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dUzdz"], 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["drodx"], 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["drody"], 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["drodz"], 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dPdx"] , 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dPdy"] , 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dPdz"] , 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dTdx"] , 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dTdy"] , 0.0f, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["dTdz"] , 0.0f, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["dHtdx"] , 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["dHtdy"] , 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["dHtdz"] , 0.0, msh.nCells*sizeof(flow_float)));
    CHECK_CUDA_ERROR(cudaMemset(var.c_d["divU"] , 0.0f, msh.nCells*sizeof(flow_float)));

    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUxdx"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUxdy"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUxdz"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUydx"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUydy"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUydz"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUzdx"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUzdy"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droUzdz"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droedx"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droedy"], 0.0, msh.nCells*sizeof(flow_float)));
    //CHECK_CUDA_ERROR(cudaMemset(var.c_d["droedz"], 0.0, msh.nCells*sizeof(flow_float)));


    // ---- LSQ 勾配 係数事前計算版 (node-centered, gradLSQ=2, methods §7.3.1) ----
    // setup 1 回 (double 固有分解 + スペクトル打ち切り擬似逆) で係数を焼き込み、
    // runtime は float32 gather のみ。退化ノード数は起動時に 1 度だけログする。
    if (cfg.gradLSQ == 2 && cfg.discretization == "node") {
        // 境界 (periodic 含む) は LSQ 点にしない (§7.3.1, 2026-08-11 改訂) — 係数テーブルは
        // 内部双対面 incidence のみ。
        static flow_float *cInt=nullptr;
        static geom_int pre_n=0;
        if (cInt==nullptr || pre_n!=msh.nCells) {
            if (cInt){cudaFree(cInt); cInt=nullptr;}
            geom_int nInc=0;
            CHECK_CUDA_ERROR(cudaMemcpy(&nInc, msh.map_cell_planes_index_d + msh.nCells,
                                        sizeof(geom_int), cudaMemcpyDeviceToHost));
            gpuErrchk(cudaMalloc(&cInt, sizeof(flow_float)*3*nInc));
            double* Minv6=nullptr; int* degD=nullptr;
            gpuErrchk(cudaMalloc(&Minv6, sizeof(double)*6*msh.nCells));
            gpuErrchk(cudaMalloc(&degD, sizeof(int)));
            CHECK_CUDA_ERROR(cudaMemset(degD, 0, sizeof(int)));
            lsqPre_setup_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
                msh.nCells, msh.map_plane_cells_d,
                msh.map_cell_planes_index_d, msh.map_cell_planes_d, msh.nNormalPlanes,
                var.c_d["ccx"],var.c_d["ccy"],var.c_d["ccz"],
                (double)cfg.gradLSQDegenThresh, Minv6, degD);
            lsqPre_coefInternal_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
                msh.nCells, msh.map_plane_cells_d,
                msh.map_cell_planes_index_d, msh.map_cell_planes_d, msh.nNormalPlanes,
                var.c_d["ccx"],var.c_d["ccy"],var.c_d["ccz"], Minv6, cInt);
            gpuErrchk( cudaPeekAtLastError() );
            gpuErrchkKernelSync();
            int degH=0;
            CHECK_CUDA_ERROR(cudaMemcpy(&degH, degD, sizeof(int), cudaMemcpyDeviceToHost));
            printf("gradLSQ=2 precomp: %d/%lld nodes spectral-truncated (thresh=%.2e)\n",
                   degH, (long long)msh.nCells, (double)cfg.gradLSQDegenThresh);
            cudaFree(Minv6); cudaFree(degD);
            pre_n = msh.nCells;
        }
        lsqPreGrad_internal_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
            msh.nCells, msh.map_plane_cells_d,
            msh.map_cell_planes_index_d, msh.map_cell_planes_d, msh.nNormalPlanes, cInt,
            var.c_d["ro"],var.c_d["Ux"],var.c_d["Uy"],var.c_d["Uz"],var.c_d["P"],var.c_d["T"],
            var.c_d["drodx"],var.c_d["drody"],var.c_d["drodz"],
            var.c_d["dUxdx"],var.c_d["dUxdy"],var.c_d["dUxdz"],
            var.c_d["dUydx"],var.c_d["dUydy"],var.c_d["dUydz"],
            var.c_d["dUzdx"],var.c_d["dUzdy"],var.c_d["dUzdz"],
            var.c_d["dPdx"],var.c_d["dPdy"],var.c_d["dPdz"],
            var.c_d["dTdx"],var.c_d["dTdy"],var.c_d["dTdz"]);
        lsqPreGrad_finish_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
            msh.nCells, var.c_d["dUxdx"], var.c_d["dUydy"], var.c_d["dUzdz"], var.c_d["divU"]);
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchkKernelSync();
        return;
    }

    // ---- LSQ 勾配 (node-centered, gradLSQ=1): 毎ステップ solve (回帰対照として残置) ----
    if (cfg.gradLSQ == 1 && cfg.discretization == "node") {
        static flow_float *Mxx=nullptr,*Mxy=nullptr,*Mxz=nullptr,*Myy=nullptr,*Myz=nullptr,*Mzz=nullptr;
        static geom_int M_n=0;
        if (Mxx==nullptr || M_n!=msh.nCells) {
            if (Mxx){cudaFree(Mxx);cudaFree(Mxy);cudaFree(Mxz);cudaFree(Myy);cudaFree(Myz);cudaFree(Mzz);}
            size_t sz=sizeof(flow_float)*msh.nCells;
            gpuErrchk(cudaMalloc(&Mxx,sz)); gpuErrchk(cudaMalloc(&Mxy,sz)); gpuErrchk(cudaMalloc(&Mxz,sz));
            gpuErrchk(cudaMalloc(&Myy,sz)); gpuErrchk(cudaMalloc(&Myz,sz)); gpuErrchk(cudaMalloc(&Mzz,sz));
            M_n=msh.nCells;
        }
        lsqGrad_accumInternal_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
            msh.nCells, msh.map_plane_cells_d, msh.map_cell_planes_index_d, msh.map_cell_planes_d, msh.nNormalPlanes,
            var.c_d["ccx"],var.c_d["ccy"],var.c_d["ccz"],
            var.c_d["ro"],var.c_d["Ux"],var.c_d["Uy"],var.c_d["Uz"],var.c_d["P"],var.c_d["T"],
            Mxx,Mxy,Mxz,Myy,Myz,Mzz,
            var.c_d["drodx"],var.c_d["drody"],var.c_d["drodz"],
            var.c_d["dUxdx"],var.c_d["dUxdy"],var.c_d["dUxdz"],
            var.c_d["dUydx"],var.c_d["dUydy"],var.c_d["dUydz"],
            var.c_d["dUzdx"],var.c_d["dUzdy"],var.c_d["dUzdz"],
            var.c_d["dPdx"],var.c_d["dPdy"],var.c_d["dPdz"],
            var.c_d["dTdx"],var.c_d["dTdy"],var.c_d["dTdz"]);
        // 境界 (periodic 含む) は LSQ 点にしない (§7.3, 2026-08-11 改訂) — accumInternal が
        // 既に境界面を skip しているため追加の除外処理は不要 (旧 lsqGrad_accumBoundary_d は削除)。
        lsqGrad_solve_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
            msh.nCells, Mxx,Mxy,Mxz,Myy,Myz,Mzz,
            var.c_d["drodx"],var.c_d["drody"],var.c_d["drodz"],
            var.c_d["dUxdx"],var.c_d["dUxdy"],var.c_d["dUxdz"],
            var.c_d["dUydx"],var.c_d["dUydy"],var.c_d["dUydz"],
            var.c_d["dUzdx"],var.c_d["dUzdy"],var.c_d["dUzdz"],
            var.c_d["dPdx"],var.c_d["dPdy"],var.c_d["dPdz"],
            var.c_d["dTdx"],var.c_d["dTdy"],var.c_d["dTdz"],
            var.c_d["divU"]);
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchkKernelSync();
        return;   // LSQ は正規化不要 (直接勾配を解く)。GG 経路はスキップ。
    }

    // sum over planes
    // K5: atomicAdd 版 calcGradient_1_d を cell-gather 版へ置換（raw sum を書く。正規化は下の _2_d）。
    calcGradient_cellgather_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
        msh.nCells,
        msh.map_plane_cells_d,
        msh.map_cell_planes_index_d , msh.map_cell_planes_d ,
        grad_sx , grad_sy , grad_sz ,
        var.p_d["fx"],
        var.c_d["ro"] , var.c_d["Ux"] , var.c_d["Uy"] , var.c_d["Uz"] , var.c_d["P"] , var.c_d["T"] ,
        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
        var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
        var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
        var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],
        var.c_d["divU"],
        // node 弱形式: 境界半割面 (ip>=nNormalPlanes) を cellgather でスキップし calcGradient_b_d
        // (owner-state 境界面値, §6.2/§7.2.2) に委ねる。壁ノードは壁上に乗りミラーゴーストが退化するため、
        // ミラー経由の cellgather 勾配が不正確 → 粘性 BL 振動の原因。cell は nPlanes (全面処理, 従来)。
        (cfg.discretization == "node") ? msh.nNormalPlanes : msh.nPlanes
    ) ;

    // node 弱形式 境界勾配: 各境界半割面の φ_b·S_b を owner ノードの状態値 φ[W] で加える
    // (cellgather は内部のみ計算済み。bvar は参照しない — §6.2/§7.2.2, 2026-08-11 改訂)。
    // 軸対称は cellgather と同じ planar 面積 (grad_sx=sx_planar) を使い r 重みを一致させる。
    if (cfg.discretization == "node") {
        for (auto& bc : msh.bconds)
        {
            // periodic は DOF 同一視 (§4.5) で扱うため境界半割面の勾配寄与を加えない。
            // 勾配は内部双対面のみ (片側) で計算し、後段 periodicGradientGather で両側を厳密合併する。
            if (bc.bcondKind == "periodic") continue;
            calcGradient_b_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> (
                bc.iPlanes.size(),
                bc.map_bplane_plane_d, bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
                var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
                grad_sx          , grad_sy       , grad_sz       , grad_ss,
                // owner ノードの状態のみ (bvar は参照しない, §6.2/§7.2.2)
                var.c_d["ro"], var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"],
                var.c_d["P"], var.c_d["T"], var.c_d["roe"], var.c_d["Ht"],
                var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
                var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
                var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
                var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
                var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
                var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],
                var.c_d["divU"]
            ) ;
        }
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchkKernelSync();
    }


    calcGradient_2_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
        // mesh structure
        msh.nCells,
        grad_volume, var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        grad_sx , grad_sy , grad_sz , grad_ss,

        // basic variables
        var.c_d["ro"] ,
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["P"]  , 
        var.c_d["T"]  ,
        var.c_d["roe"] ,
        var.c_d["Ht"] ,

        // gradient
        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
        var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
        var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"],
        var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],

        //var.c_d["droUxdx"] , var.c_d["droUxdy"] , var.c_d["droUxdz"],
        //var.c_d["droUydx"] , var.c_d["droUydy"] , var.c_d["droUydz"],
        //var.c_d["droUzdx"] , var.c_d["droUzdy"] , var.c_d["droUzdz"],
        //var.c_d["droedx"]  , var.c_d["droedy"]  , var.c_d["droedz"],
        //var.c_d["dHtdx"]  , var.c_d["dHtdy"]  , var.c_d["dHtdz"],
 
        var.c_d["divU"]
    ) ;

    // 境界隣接 CV の再構成勾配を 0 に (境界側 1 次化)。bndFirstOrder=1 のときのみ。
    if (cfg.bndFirstOrder == 1 && msh.bnode_flag_d != nullptr) {
        zeroBndNodeGradient_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
            msh.nCells, msh.bnode_flag_d,
            var.c_d["drodx"] , var.c_d["drody"] , var.c_d["drodz"],
            var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
            var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
            var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
            var.c_d["dPdx"]  , var.c_d["dPdy"]  , var.c_d["dPdz"]
        );
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

}