#include "dependentVariables_d.cuh"

__global__ void updateVariablesOuter_d
( 
 // mesh structure
 geom_int nCells_all , geom_int nCells,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,
 flow_float* roK  ,
 flow_float* roOmega  ,

 flow_float* roN ,
 flow_float* roUxN ,
 flow_float* roUyN ,
 flow_float* roUzN ,
 flow_float* roeN ,
 flow_float* roKN ,
 flow_float* roOmegaN ,
 
 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,
 flow_float* roKM ,
 flow_float* roOmegaM 

)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells_all) {

        roN[ic]   = ro[ic];
        roUxN[ic] = roUx[ic];
        roUyN[ic] = roUy[ic];
        roUzN[ic] = roUz[ic];
        roeN[ic]  = roe[ic];
        roKN[ic] = roK[ic];
        roOmegaN[ic] = roOmega[ic];

        roM[ic]   = ro[ic];
        roUxM[ic] = roUx[ic];
        roUyM[ic] = roUy[ic];
        roUzM[ic] = roUz[ic];
        roeM[ic]  = roe[ic];
        roKM[ic] = roK[ic];
        roOmegaM[ic] = roOmega[ic];

    }
}


void updateVariablesOuter_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    updateVariablesOuter_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
        // mesh structure
        msh.nCells_all , msh.nCells ,

        // basic variables
        var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"]  , var.c_d["roe"]  , var.c_d["roK"] , var.c_d["roOmega"] ,
        var.c_d["roN"] , var.c_d["roUxN"], var.c_d["roUyN"] , var.c_d["roUzN"] , var.c_d["roeN"] , var.c_d["roKN"] , var.c_d["roOmegaN"] ,
        var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"] , var.c_d["roeM"] , var.c_d["roKM"] , var.c_d["roOmegaM"] 
    ) ;

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

__global__ void updateVariablesInner_d
( 
 // mesh structure
 geom_int nCells_all , geom_int nCells,

 // variables
 flow_float* ro  ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe  ,
 flow_float* roK  ,
 flow_float* roOmega  ,

 flow_float* roM ,
 flow_float* roUxM ,
 flow_float* roUyM ,
 flow_float* roUzM ,
 flow_float* roeM ,
 flow_float* roKM ,
 flow_float* roOmegaM 

)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells_all) {

        roM[ic]   = ro[ic];
        roUxM[ic] = roUx[ic];
        roUyM[ic] = roUy[ic];
        roUzM[ic] = roUz[ic];
        roeM[ic]  = roe[ic];
        roKM[ic] = roK[ic];
        roOmegaM[ic] = roOmega[ic];
    }
}


void updateVariablesInner_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    updateVariablesInner_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> ( 
        // mesh structure
        msh.nCells_all , msh.nCells ,

        // basic variables
        var.c_d["ro"]  , var.c_d["roUx"] , var.c_d["roUy"]  , var.c_d["roUz"]  , var.c_d["roe"]  , var.c_d["roK"] , var.c_d["roOmega"] ,
        var.c_d["roM"] , var.c_d["roUxM"], var.c_d["roUyM"] , var.c_d["roUzM"] , var.c_d["roeM"] , var.c_d["roKM"] , var.c_d["roOmegaM"] 
    ) ;

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

__global__ void applyScalarImplicitCorrection_d
(
 geom_int nCells_all , geom_int nCells,

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

 flow_float* dq_ro,
 flow_float* dq_roUx,
 flow_float* dq_roUy,
 flow_float* dq_roUz,
 flow_float* dq_roe
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        ro[ic] = roN[ic] + dq_ro[ic];
        roUx[ic] = roUxN[ic] + dq_roUx[ic];
        roUy[ic] = roUyN[ic] + dq_roUy[ic];
        roUz[ic] = roUzN[ic] + dq_roUz[ic];
        roe[ic] = roeN[ic] + dq_roe[ic];
    }
}

void applyScalarImplicitCorrection_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    applyScalarImplicitCorrection_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells_all , msh.nCells,
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
        var.c_d["dq_ro_old"],
        var.c_d["dq_roUx_old"],
        var.c_d["dq_roUy_old"],
        var.c_d["dq_roUz_old"],
        var.c_d["dq_roe_old"]
    );

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

// 陰的更新の正値性ガード (plans/active/time_integration-update-positivity-guard.md):
// commit 予定の Δq が「ro または内部エネルギー e_i を 1 step で alpha 倍未満に落とす」場合、
// 半減列 s∈{1,1/2,...,1/32} で Δq を縮小する (床に当てる前の局所 under-relax)。
// P でなく (ro, e_i) を見るのは TP の EOS 反転を避けるため (CPG では P=(γ-1)e_i と等価)。
__device__ inline flow_float updateGuardScale
(
 flow_float alpha,
 flow_float ro0, flow_float rux0, flow_float ruy0, flow_float ruz0, flow_float roe0,
 flow_float d0, flow_float d1, flow_float d2, flow_float d3, flow_float d4
)
{
    const flow_float tiny = (flow_float)1e-30;
    const flow_float ei0 = roe0 - (flow_float)0.5*(rux0*rux0 + ruy0*ruy0 + ruz0*ruz0)/fmaxf(ro0, tiny);
    if (!(ro0 > (flow_float)0.0) || !(ei0 > (flow_float)0.0)) return (flow_float)1.0; // 病的セルは対象外
    flow_float s = (flow_float)1.0;
    #pragma unroll
    for (int k = 0; k < 6; ++k) {
        const flow_float ro_n = ro0 + s*d0;
        const flow_float ux = rux0 + s*d1, uy = ruy0 + s*d2, uz = ruz0 + s*d3;
        const flow_float ei_n = (roe0 + s*d4) - (flow_float)0.5*(ux*ux + uy*uy + uz*uz)/fmaxf(ro_n, tiny);
        if (ro_n >= alpha*ro0 && ei_n >= alpha*ei0) return s;
        s *= (flow_float)0.5;
    }
    return s; // 6 回半減しても不成立: s=1/64 で commit (EOS 床が最後の防波堤)
}

__global__ void applyBlockImplicitCorrection_d
(
 geom_int nCells_all , geom_int nCells,

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

 flow_float* dq_block_0,
 flow_float* dq_block_1,
 flow_float* dq_block_2,
 flow_float* dq_block_3,
 flow_float* dq_block_4,
 geom_int* axis_flag , // node-centered 軸対称: 軸上 CV で roUy=0 (対称条件)。nullptr 可。
 flow_float updateGuardAlpha   // >0 で正値性ガード有効 (0 = 既存とビット同一)
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    if (ic < nCells) {
        flow_float d0 = dq_block_0[ic], d1 = dq_block_1[ic], d2 = dq_block_2[ic],
                   d3 = dq_block_3[ic], d4 = dq_block_4[ic];
        if (updateGuardAlpha > (flow_float)0.0) {
            const flow_float s = updateGuardScale(updateGuardAlpha,
                roN[ic], roUxN[ic], roUyN[ic], roUzN[ic], roeN[ic], d0, d1, d2, d3, d4);
            d0 *= s; d1 *= s; d2 *= s; d3 *= s; d4 *= s;
        }
        ro[ic]   = roN[ic]   + d0;
        roUx[ic] = roUxN[ic] + d1;
        // 軸上 CV は半径方向運動量を 0 に保つ (SU2 MARKER_SYM 流。block-DPLUR の連成で
        // 補正が漏れるのを commit 時に射影。roUy~0 なので KE 不整合は無視できる)。
        roUy[ic] = (axis_flag != nullptr && axis_flag[ic] == 1)
                       ? (flow_float)0.0 : roUyN[ic] + d2;
        roUz[ic] = roUzN[ic] + d3;
        roe[ic]  = roeN[ic]  + d4;
    }
}

void applyBlockImplicitCorrection_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    // 古典 DPLUR の sweep ループ後、最終補正は dq_block_old に入っている（最後の swap 後）。
    applyBlockImplicitCorrection_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells_all , msh.nCells,
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
        var.c_d["dq_block_old_0"],
        var.c_d["dq_block_old_1"],
        var.c_d["dq_block_old_2"],
        var.c_d["dq_block_old_3"],
        var.c_d["dq_block_old_4"],
        // 注: 軸上 roUy=0 を commit で強制すると block-DPLUR の線形 solve と不整合になり発散級に悪化
        // (Mach~1000)。SU2 流の対称面は Jacobian を整合的に修正する必要がある (block-DPLUR の row 修正)。
        // 暫定で無効 (nullptr=baseline)。near-axis corner は open issue (docs §7.1)。
        nullptr,
        cfg.updateGuardAlpha
    );

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

// SST (k-ω) の segregated point-implicit 更新。
// 平均流 block DPLUR の commit 後に呼び、凍結していたスカラーを陰的に更新する。
// 消散項 + 輸送項(移流+拡散) を point-implicit 化: D_φ = V/Δτ + V·(∂D/∂φ) + transport_diag、
// δ(ρφ) = relax·res_roφ / D_φ。生産・近傍 ΔQ は res_roφ に含む lagged 扱い。ρφ = ρφ_baseline(N) + δ。
__global__ void applySSTPointImplicit_d
(
 geom_int nCells,
 geom_float* vol,
 flow_float* dt_local,
 flow_float implicit_relax,

 flow_float* roK,
 flow_float* roOmega,
 flow_float* roKN,
 flow_float* roOmegaN,

 flow_float* res_roK,
 flow_float* res_roOmega,
 flow_float* src_jac_k,
 flow_float* src_jac_omega,
 flow_float* transport_diag_k,
 flow_float* transport_diag_omega,
 flow_float unsteady_diag,
 // node-centered 壁 ω Dirichlet: 壁ノードで dω=0 にし rans_wall_scalar_boundary_d がピンした ω_w を保つ
 // (point-implicit の行=1/他=0)。k はノイマンなので通常更新。nullptr/0 で無効 (cell 不変)。
 geom_int* wall_flag, int decouple_wall_omega,
 // node-centered k Dirichlet (SU2 SetTurbVars_WF 流): roK_wf[ic]>=0 の第一内層ノードで roK=roK_wf に固定
 // (dk=0)。near-wall k 蓄積 (再付着 μ_t ピーク) を断つ。nullptr/全-1 で無効 (cell 不変)。
 flow_float* roK_wf,
 flow_float* roOmega_wf
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float v    = vol[ic];
        const flow_float dt_l = max(dt_local[ic], static_cast<flow_float>(1.0e-30));
        const bool omegaWall = (decouple_wall_omega != 0 && wall_flag != nullptr && wall_flag[ic] == 1);
        const bool kPinned   = (roK_wf != nullptr && roK_wf[ic] >= static_cast<flow_float>(0.0));
        const bool omgPinnedWf = (roOmega_wf != nullptr && roOmega_wf[ic] >= static_cast<flow_float>(0.0));
        // D_φ = V/Δτ + V·src_jac（消散）+ transport_diag（移流+拡散、既に [m³/s]）+ V·unsteady_diag。
        // transport_diag は壁近傍の陽的 k/ω 輸送 stiff 性を陰化し安定 cfl_pseudo を一桁以上引き上げる。
        const flow_float Dk = v / dt_l + v * src_jac_k[ic]     + transport_diag_k[ic]     + v * unsteady_diag;
        const flow_float Dw = v / dt_l + v * src_jac_omega[ic] + transport_diag_omega[ic] + v * unsteady_diag;

        // k Dirichlet (node 第一内層ノード): dk で更新せず roK_wf に固定。それ以外は通常 point-implicit。
        const flow_float dk = implicit_relax * res_roK[ic]     / max(Dk, static_cast<flow_float>(1.0e-30));
        // 壁 ω decouple: 壁ノードでは dω=0 とし、ピンした roOmega=ρ·ω_w を保つ (k はノイマンで通常更新)。
        const flow_float dw = omegaWall ? static_cast<flow_float>(0.0)
                                        : implicit_relax * res_roOmega[ic] / max(Dw, static_cast<flow_float>(1.0e-30));

        // realizability: ρk ≥ 0, ρω > 0。dual-time でも正しいよう現在反復値への in-place 加算
        // （定常では roKN==roK のため roKN+dk と等価）。dual-time の roKN/roKNN は BDF 残差項側で使用。
        roK[ic]     = kPinned ? roK_wf[ic] : max(roK[ic] + dk, static_cast<flow_float>(0.0));
        // ω Dirichlet (node 第一内層ノード, roOmega_wf): dw で更新せず ρω_w に固定 (k と対)。
        roOmega[ic] = omgPinnedWf ? roOmega_wf[ic]
                                  : max(roOmega[ic] + dw, static_cast<flow_float>(1.0e-20));
    }
}

void applySSTPointImplicit_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    applySSTPointImplicit_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells,
        var.c_d["volume"],
        var.c_d["dt_local"],
        cfg.implicitRelaxSST,
        var.c_d["roK"],
        var.c_d["roOmega"],
        var.c_d["roKN"],
        var.c_d["roOmegaN"],
        var.c_d["res_roK"],
        var.c_d["res_roOmega"],
        var.c_d["src_jac_k"],
        var.c_d["src_jac_omega"],
        var.c_d["transport_diag_k"],
        var.c_d["transport_diag_omega"],
        cfg.unsteadyDiagCoef,
        (cfg.discretization == "node") ? msh.wall_flag_d : nullptr,
        (cfg.discretization == "node" && msh.wall_flag_d != nullptr) ? 1 : 0,
        // k Dirichlet は wallTreatmentSST==1 && nodeKwfDirichlet==1 の node のみ。それ以外は nullptr で無効化
        // (roK_wf が init されない経路で誤発火しないよう堅牢化)。
        (cfg.discretization == "node" && cfg.wallTreatmentSST == 1 && cfg.nodeKwfDirichlet == 1) ? var.c_d["roK_wf"] : nullptr,
        (cfg.discretization == "node" && cfg.wallTreatmentSST == 1 && cfg.nodeOmegaWfDirichlet == 1) ? var.c_d["roOmega_wf"] : nullptr
    );

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

// dual-time: 物理時間 BDF 項を残差に加える（assembleResidual の末尾で呼ぶ）。
// res は -R_spatial（陽解法 Q=Q_N+res·dt/V と整合）。dual-time の残差 R* = R_spatial + (V/Δt)(a Q - b Q^n + c Q^{n-1})
// を 0 にするので res* = res - (V/Δt)(a Q - b Q^n + c Q^{n-1}) とする（roN=Q^n, roNN=Q^{n-1}）。
// (a,b,c): BDF1=(1,1,0), BDF2=(1.5,2,0.5)。include_scalar==1 で k/ω も処理。
__global__ void addUnsteadyTimeTerm_d
(
 geom_int nCells,
 geom_float* vol,
 flow_float dt_phys,
 flow_float a, flow_float b, flow_float c,

 flow_float* ro,   flow_float* roUx,   flow_float* roUy,   flow_float* roUz,   flow_float* roe,
 flow_float* roN,  flow_float* roUxN,  flow_float* roUyN,  flow_float* roUzN,  flow_float* roeN,
 flow_float* roNN, flow_float* roUxNN, flow_float* roUyNN, flow_float* roUzNN, flow_float* roeNN,
 flow_float* res_ro, flow_float* res_roUx, flow_float* res_roUy, flow_float* res_roUz, flow_float* res_roe,

 int include_scalar,
 flow_float* roK,   flow_float* roOmega,
 flow_float* roKN,  flow_float* roOmegaN,
 flow_float* roKNN, flow_float* roOmegaNN,
 flow_float* res_roK, flow_float* res_roOmega,
 int energyK   // sstEnergyIncludesK: エネルギー行の BDF を E_t = roe + roK に (plan turbulence-sst-energy-includes-k §4.4a)
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float coef = vol[ic] / max(dt_phys, static_cast<flow_float>(1.0e-30));
        res_ro[ic]   -= coef * (a*ro[ic]   - b*roN[ic]   + c*roNN[ic]);
        res_roUx[ic] -= coef * (a*roUx[ic] - b*roUxN[ic] + c*roUxNN[ic]);
        res_roUy[ic] -= coef * (a*roUy[ic] - b*roUyN[ic] + c*roUyNN[ic]);
        res_roUz[ic] -= coef * (a*roUz[ic] - b*roUzN[ic] + c*roUzNN[ic]);
        res_roe[ic]  -= coef * (a*roe[ic]  - b*roeN[ic]  + c*roeNN[ic]);
        if (energyK != 0 && include_scalar) {
            // E_t の時間項: 平均流 E_m の項に加えて ρk の BDF を足す (k 式側の時間項とは別勘定。E_t 式の保存は subiter 収束で成立)
            res_roe[ic] -= coef * (a*roK[ic] - b*roKN[ic] + c*roKNN[ic]);
        }
        if (include_scalar) {
            res_roK[ic]     -= coef * (a*roK[ic]     - b*roKN[ic]     + c*roKNN[ic]);
            res_roOmega[ic] -= coef * (a*roOmega[ic] - b*roOmegaN[ic] + c*roOmegaNN[ic]);
        }
    }
}

void addUnsteadyTimeTerm_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var,
                                   flow_float a, flow_float b, flow_float c, int include_scalar)
{
    addUnsteadyTimeTerm_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells,
        var.c_d["volume"],
        cfg.dt,
        a, b, c,
        var.c_d["ro"],   var.c_d["roUx"],   var.c_d["roUy"],   var.c_d["roUz"],   var.c_d["roe"],
        var.c_d["roN"],  var.c_d["roUxN"],  var.c_d["roUyN"],  var.c_d["roUzN"],  var.c_d["roeN"],
        var.c_d["roNN"], var.c_d["roUxNN"], var.c_d["roUyNN"], var.c_d["roUzNN"], var.c_d["roeNN"],
        var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roUz"], var.c_d["res_roe"],
        include_scalar,
        var.c_d["roK"],   var.c_d["roOmega"],
        var.c_d["roKN"],  var.c_d["roOmegaN"],
        var.c_d["roKNN"], var.c_d["roOmegaNN"],
        var.c_d["res_roK"], var.c_d["res_roOmega"],
        (cfg.sstEnergyIncludesK != 0 && cfg.LESorRANS == 2 && cfg.RANSmodel == 1) ? 1 : 0
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

// dual-time: 擬似時間サブ反復の commit。roN=Q^n は BDF 基準として固定するため、
// 定常の Q=roN+dq ではなく現在反復値への in-place 加算 Q += dq を行う。
__global__ void applyBlockImplicitCorrectionInPlace_d
(
 geom_int nCells,
 flow_float* ro, flow_float* roUx, flow_float* roUy, flow_float* roUz, flow_float* roe,
 flow_float* dq0, flow_float* dq1, flow_float* dq2, flow_float* dq3, flow_float* dq4,
 flow_float updateGuardAlpha
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        flow_float d0 = dq0[ic], d1 = dq1[ic], d2 = dq2[ic], d3 = dq3[ic], d4 = dq4[ic];
        if (updateGuardAlpha > (flow_float)0.0) {
            const flow_float s = updateGuardScale(updateGuardAlpha,
                ro[ic], roUx[ic], roUy[ic], roUz[ic], roe[ic], d0, d1, d2, d3, d4);
            d0 *= s; d1 *= s; d2 *= s; d3 *= s; d4 *= s;
        }
        ro[ic]   += d0;
        roUx[ic] += d1;
        roUy[ic] += d2;
        roUz[ic] += d3;
        roe[ic]  += d4;
    }
}

void applyBlockImplicitCorrectionInPlace_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    applyBlockImplicitCorrectionInPlace_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells,
        var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"],
        var.c_d["dq_block_old_0"], var.c_d["dq_block_old_1"], var.c_d["dq_block_old_2"],
        var.c_d["dq_block_old_3"], var.c_d["dq_block_old_4"],
        cfg.updateGuardAlpha
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

// dual-time の物理時間レベルシフト: roNN ← roN, roN ← ro（mean flow 5 + scalar 2）。
__global__ void shiftDualTimeLevels_d
(
 geom_int nCells_all,
 flow_float* ro, flow_float* roUx, flow_float* roUy, flow_float* roUz, flow_float* roe,
 flow_float* roN, flow_float* roUxN, flow_float* roUyN, flow_float* roUzN, flow_float* roeN,
 flow_float* roNN, flow_float* roUxNN, flow_float* roUyNN, flow_float* roUzNN, flow_float* roeNN,
 flow_float* roK, flow_float* roOmega,
 flow_float* roKN, flow_float* roOmegaN,
 flow_float* roKNN, flow_float* roOmegaNN
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells_all) {
        roNN[ic]=roN[ic]; roUxNN[ic]=roUxN[ic]; roUyNN[ic]=roUyN[ic]; roUzNN[ic]=roUzN[ic]; roeNN[ic]=roeN[ic];
        roN[ic]=ro[ic];   roUxN[ic]=roUx[ic];   roUyN[ic]=roUy[ic];   roUzN[ic]=roUz[ic];   roeN[ic]=roe[ic];
        roKNN[ic]=roKN[ic]; roOmegaNN[ic]=roOmegaN[ic];
        roKN[ic]=roK[ic];   roOmegaN[ic]=roOmega[ic];
    }
}

void shiftDualTimeLevels_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    shiftDualTimeLevels_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells_all,
        var.c_d["ro"], var.c_d["roUx"], var.c_d["roUy"], var.c_d["roUz"], var.c_d["roe"],
        var.c_d["roN"], var.c_d["roUxN"], var.c_d["roUyN"], var.c_d["roUzN"], var.c_d["roeN"],
        var.c_d["roNN"], var.c_d["roUxNN"], var.c_d["roUyNN"], var.c_d["roUzNN"], var.c_d["roeNN"],
        var.c_d["roK"], var.c_d["roOmega"],
        var.c_d["roKN"], var.c_d["roOmegaN"],
        var.c_d["roKNN"], var.c_d["roOmegaNN"]
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}