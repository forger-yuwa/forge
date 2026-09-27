#include "gasProperties_d.cuh"
#include "thermo_d.cuh"
#include "transportMix_d.cuh"      // 種ごとの輸送物性 (physProp.transport; plan #5t2-2)
#include "speciesTransport_d.cuh"  // species_roY_device_ptr()

// 種ごとの輸送物性の経路 (physProp.transport があり viscMethod 2) で使うセル組成。
//   現行 viscMethod 2 と同じ作り方 (Y_s = max(ρY_s/ρ, 0) を和で正規化)。単成分・roY 無しは {1, 0, ...}。
//   返り値は transport_mix_Y へ渡す種数 (roY 無しは 1 = X {1})。
__device__ inline int gas_transport_cell_Y(int nSpecies, flow_float* const* roY, const flow_float* ro, geom_int ic, double* Y)
{
    if (nSpecies <= 1 || roY == nullptr) {
        Y[0] = 1.0;
        for (int s = 1; s < nSpecies && s < THERMO_MAX_SPECIES; s++) Y[s] = 0.0;
        return 1;
    }
    const double ro_d = (double)max(ro[ic], (flow_float)1.0e-30);
    double ysum = 0.0;
    for (int s=0;s<nSpecies;s++){ double y=(double)roY[s][ic]/ro_d; if(y<0.0)y=0.0; Y[s]=y; ysum+=y; }
    const double inv = 1.0/(ysum>1.0e-30?ysum:1.0e-30);
    for (int s=0;s<nSpecies;s++) Y[s]*=inv;
    return nSpecies;
}

__global__ void gasProperties_d
(
 // gas properties
 int thermalMethod , int viscMethod ,
 // 層流熱伝導モデル: 0=一定 / 1=constant-Pr (k=μ(T)·cp/Pr_lam)。viscMethod 0/1 で有効 (2 は kinetic theory が λ を持つ)。
 int thermCondMethod , flow_float prandtlLam ,

 // gas properties
 flow_float gamma , flow_float cp , flow_float visc_lam, flow_float thermCond_const,

 // thermally-perfect 化学種データ (kinetic theory 輸送 viscMethod==2 用)
 const SpeciesThermo* sp , int nSpecies , flow_float** roY ,

 // 種ごとの輸送物性 (physProp.transport; transportOn=1 のとき viscMethod 2 は CEA frozen 混合則 transport_mix_Y)
 int transportOn , TransportTableD ttab ,

 // mesh structure
 geom_int nCells_all , geom_int nCells,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,

 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,
 flow_float* T   ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,

 flow_float* vis_lam_array,
 flow_float* thermCond_array,
 flow_float* gam_array,
 flow_float* cp_array

)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells_all) {

        if (thermalMethod == 0) {
            gam_array[ic] = gamma;
            cp_array[ic] = cp;
        }

        if (viscMethod == 0) {
            vis_lam_array[ic]   = visc_lam;
            thermCond_array[ic] = (thermCondMethod == 1)
                ? vis_lam_array[ic]*cp_array[ic]/prandtlLam    // constant-Pr
                : thermCond_const;                             // 既存どおり一定 (face 平均で同値)

        } else if (viscMethod == 1) { // sutherland (熱伝導は thermCondMethod で選択)
            flow_float T0  = 273.0f;
            flow_float mu0 = 1.716e-5f;
            flow_float Smu = 111.0f;
            vis_lam_array[ic]   = mu0*pow(T[ic]/T0,3.0f/2.0f)*(T0+Smu)/(T[ic]+Smu);
            // constant-Pr: k=μ(T)·cp/Pr で分子 Pr を一定に保つ (一定 k は高温で Pr~1.5-1.9 に漂う)。
            thermCond_array[ic] = (thermCondMethod == 1)
                ? vis_lam_array[ic]*cp_array[ic]/prandtlLam
                : thermCond_const;

        } else if (viscMethod == 2 && transportOn != 0) { // 種ごとの出所 + CEA frozen 混合則 (double 評価・float 格納)
            double Y[THERMO_MAX_SPECIES];
            const int nY = gas_transport_cell_Y(nSpecies, roY, ro, ic, Y);
            double mu, lam;
            transport_mix_Y(sp, ttab, nY, Y, (double)T[ic], &mu, &lam);
            vis_lam_array[ic]   = (flow_float)mu;
            thermCond_array[ic] = (flow_float)lam;

        } else if (viscMethod == 2) { // kinetic theory (Chapman-Enskog + Wilke/Mason-Saxena)
            // 組成 Y (単成分 or roY=nullptr のときは Y={1}) を構築し double で評価。
            double Y[THERMO_MAX_SPECIES];
            double X[THERMO_MAX_SPECIES];
            const double ro_d = (double)max(ro[ic], (flow_float)1.0e-30);
            if (nSpecies <= 1 || roY == nullptr) {
                Y[0] = 1.0f;
            } else {
                double ysum = 0.0;
                for (int s=0;s<nSpecies;s++){ double y=(double)roY[s][ic]/ro_d; if(y<0.0)y=0.0; Y[s]=y; ysum+=y; }
                const double inv = 1.0/(ysum>1.0e-30?ysum:1.0e-30);
                for (int s=0;s<nSpecies;s++) Y[s]*=inv;
            }
            const int n = (nSpecies >= 1) ? nSpecies : 1;
            thermo_X_from_Y(sp, n, Y, X);
            const double Td = (double)T[ic];
            vis_lam_array[ic]   = (flow_float)thermo_mu_mix(sp, n, X, Td);
            thermCond_array[ic] = (flow_float)thermo_lambda_mix(sp, n, X, Td);
        }
    }
}


void gasProperties_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    const TransportTableD* ttab = thermo_transport_table();   // physProp.transport が無ければ nullptr (現行経路)
    gasProperties_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
        cfg.thermalMethod,
        cfg.viscMethod ,
        cfg.thermCondMethod , cfg.prandtlLam ,

        // gas properties
        cfg.gamma , cfg.cp , cfg.visc, cfg.thermCond,

        // 化学種データ (kinetic theory)
        thermo_species_device_ptr() , cfg.nSpecies , species_roY_device_ptr() ,

        // 種ごとの輸送物性 (physProp.transport)
        ttab ? 1 : 0 , ttab ? *ttab : TransportTableD{0, 0, nullptr, nullptr, nullptr} ,

        // mesh structure
        msh.nCells_all , msh.nCells ,

        // basic variables
        var.c_d["ro"]  , var.c_d["roUx"], var.c_d["roUy"] , var.c_d["roUz"], var.c_d["roe"] ,
        var.c_d["P"]   , var.c_d["Ht"]  , var.c_d["sonic"], var.c_d["T"],
        var.c_d["Ux"]  , var.c_d["Uy"]  , var.c_d["Uz"] ,

        var.c_d["vis_lam"] , var.c_d["thermCond"] , var.c_d["gamma"] , var.c_d["cp"]
    ) ;
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}


// -----------------------------------------------------------------------------
// 試験用 (FORGE_TRANSPORT_PROBE; tests/unit/test_transport_gpu.py)。計算経路からは呼ばない。
//   セル経路と同じ組成の作り方 (gas_transport_cell_Y) と同じ評価関数 (transport_mix_Y) で、double の μ・λ と
//   展開後の実種モル分率を全セル (ghost 込み nCells_all) について書き出す。
// -----------------------------------------------------------------------------
__global__ void gasPropertiesTransportProbe_d(const SpeciesThermo* sp, int nSpecies, flow_float** roY, TransportTableD ttab,
                                              geom_int nCells_all, flow_float* ro, flow_float* T,
                                              double* mu, double* lam, double* Xreal)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells_all) return;
    double Y[THERMO_MAX_SPECIES];
    const int nY = gas_transport_cell_Y(nSpecies, roY, ro, ic, Y);
    transport_mix_Y(sp, ttab, nY, Y, (double)T[ic], &mu[ic], &lam[ic], Xreal + (size_t)ic*ttab.nReal);
}

// 試験用: 与えた double の輸送種 Y (K×nSpecies, 行優先) と T で同じ評価関数を呼ぶ (float の roY を経由しない比較用)。
__global__ void gasPropertiesTransportProbeStates_d(const SpeciesThermo* sp, int nSpecies, TransportTableD ttab, int K,
                                                    const double* Y, const double* T, double* mu, double* lam, double* Xreal)
{
    const int k = blockDim.x*blockIdx.x + threadIdx.x;
    if (k >= K) return;
    transport_mix_Y(sp, ttab, nSpecies, Y + (size_t)k*nSpecies, T[k], &mu[k], &lam[k], Xreal + (size_t)k*ttab.nReal);
}

bool gasPropertiesTransportProbe_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                                           double* mu_d, double* lam_d, double* Xreal_d)
{
    const TransportTableD* ttab = thermo_transport_table();
    if (ttab == nullptr) return false;
    gasPropertiesTransportProbe_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        thermo_species_device_ptr(), cfg.nSpecies, species_roY_device_ptr(), *ttab,
        msh.nCells_all, var.c_d["ro"], var.c_d["T"], mu_d, lam_d, Xreal_d);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
    return true;
}

bool gasPropertiesTransportProbeStates_d_wrapper(solverConfig& cfg, int K, const double* Y_d, const double* T_d,
                                                 double* mu_d, double* lam_d, double* Xreal_d)
{
    const TransportTableD* ttab = thermo_transport_table();
    if (ttab == nullptr) return false;
    const int bs = 128;
    gasPropertiesTransportProbeStates_d<<<(K + bs - 1)/bs, bs>>>(
        thermo_species_device_ptr(), cfg.nSpecies, *ttab, K, Y_d, T_d, mu_d, lam_d, Xreal_d);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
    return true;
}
