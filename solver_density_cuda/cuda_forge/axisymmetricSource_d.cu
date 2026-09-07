#include "axisymmetricSource_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"


// 軸対称 (B 流儀) — 半径方向運動量の圧力・粘性ソース項を residual に加算する。
//
// 円柱座標で書いた半径運動量保存則を r 倍してから (r,z) 平面で積分すると、
// 圧力勾配項のうち r 倍では落ちないテルムとして
//
//     R_{rho u_r} = ∫∫ (p / r) · r dr dz = p_cell · A_planar
//
// が残る。粘性は $-\tau_{\theta\theta}$ を同じ面積で積分した source として加える
// (詳しくは methods/axisymmetric/theory.md)。
//
// res_roUy には -∑(F·S) が累積されており、時間進行は
//   roUy^{n+1} = roUy^n + coef_DT * res_roUy * dt / V
// で書かれる (timeIntegration_d.cu)。よってソース項は res_roUy に
// "+= (p - tau_theta_theta) * A_planar" として加算する。
__global__ void axisymmetricSource_d
(
    geom_int nCells,
    flow_float* P,
    flow_float* A_planar,
    flow_float* vis_lam,
    flow_float* vis_turb,
    flow_float* axisym_divU,
    flow_float* axisym_uy_over_r,
    flow_float* res_roUy,
    geom_int* axis_flag,   // node-centered で軸上 CV を示す (nullptr 可)。軸上は半径方向ソースを課さない。
    flow_float* ccy, flow_float axisRFloor,   // r 床 (axisRFloor>0): A_planar は離散閉性 y 面積に置換済
    flow_float* res_roUx, flow_float* A_closure_x,  // 床の x 閉性欠損補正 (床なしでは 0)
    // free-stream 保存の基準静圧 (space.pRef)。対流流束は (p_tilde − pRef)·S で組まれるので、
    // hoop ソースも同じゲージ (P − pRef)·A にしないと一様圧 p=pRef で source だけが残り、
    // 偽の半径力 pRef·A が立つ (case/43 自由流テストで 1.25e6 m/s² を実測)。既定 pRef=0 で従来どおり。
    flow_float pRef,
    flow_float* roK_e, int energyK   // sstEnergyIncludesK: hoop 圧力に p* = p + (2/3)ρk (nullptr/0 で無効)
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        // 軸上 CV (R=0) は対称条件で半径方向の正味ソースが 0。node-centered では CV が軸に乗り、
        // P*A_planar/(ro*r_eff) が極端に stiff になって偽の u_r を駆動するため、ソースを課さない。
        if (axis_flag != nullptr && axis_flag[ic] == 1) return;
        const flow_float mu_total = vis_lam[ic] + vis_turb[ic];
        const flow_float tau_theta_theta =
            (flow_float)2.0 * mu_total * axisym_uy_over_r[ic]
            - (flow_float)(2.0 / 3.0) * mu_total * axisym_divU[ic];
        // axisRFloor>0: A_planar = 床適用後の Σ_f S_f,y (全面床の CV で 0)。τθθ は床帯で
        // uy_over_r=0 のため発散項のみだが面積も ~0 で実質不活性。x は圧力の閉性欠損補正のみ。
        const flow_float pk = (energyK != 0 && roK_e != nullptr) ? (flow_float)(2.0/3.0)*max(roK_e[ic], (flow_float)0.0) : (flow_float)0.0;
        res_roUy[ic] += (P[ic] + pk - pRef - tau_theta_theta) * A_planar[ic];
        if (axisRFloor > (flow_float)0.0) {
            res_roUx[ic] += (P[ic] + pk - pRef) * A_closure_x[ic];
        }
    }
}

// 軸対称ソース項の前計算。半径方向ひずみ u_r/r と発散 divU = ∂xUx + ∂yUy + u_r/r を
// セル毎に求める。これらは turbulent_viscosity (S_θθ) / ransSource / axisymmetricSource が
// 消費するため、それらより前に実行する必要がある。
// 有効半径 r_eff は revolved 体積 (per-radian, V = r·A_planar) を planar 面積で割って得る。
__global__ void axisymmetricGeomTerms_d
(
    geom_int nCells,
    flow_float* Uy,
    flow_float* volume,
    flow_float* A_planar,
    flow_float* dUxdx,
    flow_float* dUydy,
    flow_float* axisym_uy_over_r,
    flow_float* axisym_divU,
    // axisymMethod==1 (SU2 流): r_eff = V/A_planar は planar 幾何で成立しないため
    // セル重心 y を直接使う。軸ノード (axis_flag==1, node のみ) は SU2 同様 0 (y=0 の対称極限)。
    int axisymMethod, geom_float* ccy, geom_int* axis_flag, flow_float axisRFloor
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        flow_float r_eff;
        if (axisymMethod == 1) {
            const bool onAxis = (axis_flag != nullptr && axis_flag[ic] == 1);
            r_eff = onAxis ? (flow_float)0.0 : max((flow_float)ccy[ic], (flow_float)0.0);
        } else {
            const flow_float area = A_planar[ic];
            r_eff = (area > (flow_float)0.0) ? volume[ic] / area : (flow_float)0.0;
            // axisRFloor 帯 (r 床) は planar 扱い: hoop ひずみ u_r/r を課さない
            if (axisRFloor > (flow_float)0.0 && ccy[ic] < axisRFloor) r_eff = (flow_float)0.0;
        }

        if (r_eff > (flow_float)0.0) {
            axisym_uy_over_r[ic] = Uy[ic] / r_eff;
        } else {
            axisym_uy_over_r[ic] = (flow_float)0.0;
        }

        axisym_divU[ic] = dUxdx[ic] + dUydy[ic] + axisym_uy_over_r[ic];
    }
}


void axisymmetricSource_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.isAxisymmetric != 1 || cfg.axisymMethod != 0) return;   // hoop 源は r 重み方式専用

    // 軸 (R≤eps) でソース OFF。ソースは r 重みされない唯一の項で source/volume=P/r→∞ (r→0) が発散源。
    // node-centered 軸対称でのみ作用 (axis_flag = ノード R≈0)。SU2 の y<EPS ソース OFF に相当。
    geom_int* axis_flag = (cfg.discretization == "node" && cfg.isAxisymmetric == 1) ? msh.axis_flag_d : nullptr;
    axisymmetricSource_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells,
        var.c_d["P"],
        (cfg.axisRFloor > (flow_float)0.0 || cfg.hoopAreaFromClosure == 1)
            ? var.c_d["A_closure_y"] : var.c_d["A_planar"],
        var.c_d["vis_lam"],
        var.c_d["vis_turb"],
        var.c_d["axisym_divU"],
        var.c_d["axisym_uy_over_r"],
        var.c_d["res_roUy"],
        axis_flag,
        var.c_d["ccy"],
        (cfg.hoopAreaFromClosure == 1 && cfg.axisRFloor <= (flow_float)0.0) ? (flow_float)1.0e-30 : cfg.axisRFloor,
        var.c_d["res_roUx"], var.c_d["A_closure_x"],
        cfg.pRef,
        (cfg.sstEnergyIncludesK != 0 && cfg.LESorRANS == 2 && cfg.RANSmodel == 1) ? var.c_d["roK"] : nullptr,
        (cfg.sstEnergyIncludesK != 0 && cfg.LESorRANS == 2 && cfg.RANSmodel == 1) ? 1 : 0
    );

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

// ===================== axisymMethod==1: SU2 流 planar + 1/y ソース =====================
// SU2 `CSourceAxisymmetric_Flow` / `ResidualDiffusion` の移植 (plan axisymmetric-su2-source-formulation)。
// 幾何は planar のまま、軸対称の効果を全てセルソースで表す。軸ノード (axis_flag==1 / y≤eps) はソース 0。
// forge 符号: res += 物理ソース (SU2 は LinSysRes += flux-out 側なので符号反転して移植)。

// AuxVar 計算: A0 = μ_tot·v/y, A1 = A0·v, A2 = A0·u (SU2 ComputeAxisymmetricAuxGradients)。
__global__ void axisymAuxFields_d
(
    geom_int nCells,
    flow_float* Ux, flow_float* Uy,
    flow_float* vis_lam, flow_float* vis_turb,
    flow_float* ccy, geom_int* axis_flag,
    flow_float* aux0, flow_float* aux1, flow_float* aux2
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const bool onAxis = (axis_flag != nullptr && axis_flag[ic] == 1);
    const flow_float y = ccy[ic];
    if (onAxis || y <= (flow_float)1.0e-12) {
        aux0[ic] = (flow_float)0.0; aux1[ic] = (flow_float)0.0; aux2[ic] = (flow_float)0.0;
        return;
    }
    const flow_float mu_tot = vis_lam[ic] + max(vis_turb[ic], (flow_float)0.0);
    const flow_float a0 = mu_tot * Uy[ic] / y;
    aux0[ic] = a0;
    aux1[ic] = a0 * Uy[ic];
    aux2[ic] = a0 * Ux[ic];
}

// AuxVar の Green-Gauss 勾配 (planar, per-cell gather)。境界半割面/ghost は owner 値で閉じる。
// 必要成分のみ: ∂A0/∂x, ∂A0/∂y, ∂A1/∂y, ∂A2/∂x。
__global__ void axisymAuxGradGG_d
(
    geom_int nCells,
    geom_int* plane_cells, geom_int* cell_planes_index, geom_int* cell_planes,
    flow_float* sx, flow_float* sy, flow_float* volume,
    flow_float* aux0, flow_float* aux1, flow_float* aux2,
    flow_float* dA0dx, flow_float* dA0dy, flow_float* dA1dy, flow_float* dA2dx
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    flow_float g0x = 0.0, g0y = 0.0, g1y = 0.0, g2x = 0.0;
    const geom_int pb = cell_planes_index[ic];
    const geom_int pe = cell_planes_index[ic + 1];
    for (geom_int po = pb; po < pe; ++po) {
        const geom_int ip = cell_planes[po];
        const geom_int ic0 = plane_cells[2*ip + 0];
        const geom_int ic1 = plane_cells[2*ip + 1];
        const geom_int other = (ic0 == ic) ? ic1 : ic0;
        const flow_float sgn = (ic0 == ic) ? (flow_float)1.0 : (flow_float)-1.0;
        // ghost (境界) は owner 値で閉じる (SU2 GG の境界 vertex 値と同旨の 0 次閉包)
        const flow_float a0f = (other < nCells) ? (flow_float)0.5*(aux0[ic] + aux0[other]) : aux0[ic];
        const flow_float a1f = (other < nCells) ? (flow_float)0.5*(aux1[ic] + aux1[other]) : aux1[ic];
        const flow_float a2f = (other < nCells) ? (flow_float)0.5*(aux2[ic] + aux2[other]) : aux2[ic];
        g0x += a0f * sgn * sx[ip];
        g0y += a0f * sgn * sy[ip];
        g1y += a1f * sgn * sy[ip];
        g2x += a2f * sgn * sx[ip];
    }
    const flow_float iv = (flow_float)1.0 / max(volume[ic], (flow_float)1.0e-30);
    dA0dx[ic] = g0x * iv;
    dA0dy[ic] = g0y * iv;
    dA1dy[ic] = g1y * iv;
    dA2dx[ic] = g2x * iv;
}

// SU2 流 軸対称ソース本体 (非粘性 + 粘性)。res += S·V (forge 符号)。
__global__ void axisymmetricSourceSU2_d
(
    geom_int nCells,
    flow_float* vol,          // planar 体積
    flow_float* ccy, geom_int* axis_flag,
    flow_float* ro, flow_float* Ux, flow_float* Uy,
    flow_float* roe, flow_float* P,
    flow_float* vis_lam, flow_float* vis_turb, flow_float* thermCond, flow_float* cp_arr,
    flow_float Prt,
    flow_float* dUxdx, flow_float* dUxdy, flow_float* dUydx, flow_float* dUydy,
    flow_float* dTdy,
    flow_float* dA0dx, flow_float* dA0dy, flow_float* dA1dy, flow_float* dA2dx,
    int viscous,
    flow_float* res_ro, flow_float* res_roUx, flow_float* res_roUy, flow_float* res_roe,
    flow_float* roK_e, int energyK   // sstEnergyIncludesK: H → (E_t + p*)/ρ (nullptr/0 で無効)
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const bool onAxis = (axis_flag != nullptr && axis_flag[ic] == 1);
    const flow_float y = ccy[ic];
    if (onAxis || y <= (flow_float)1.0e-12) return;   // SU2: 軸上はソース 0

    const flow_float yinv = (flow_float)1.0 / y;
    const flow_float v = vol[ic];
    const flow_float rho = ro[ic];
    const flow_float u = Ux[ic];
    const flow_float w = Uy[ic];   // 半径方向速度
    // sstEnergyIncludesK: H* = (E_m + ρk + p + (2/3)ρk)/ρ
    const flow_float rk_e = (energyK != 0 && roK_e != nullptr) ? max(roK_e[ic], (flow_float)0.0) : (flow_float)0.0;
    const flow_float H = (roe[ic] + P[ic] + (flow_float)(5.0/3.0)*rk_e) / max(rho, (flow_float)1.0e-30);

    // 非粘性 (SU2 residual を符号反転): S = -(1/y)·[ρv, ρuv, ρv², ρvH]
    flow_float s0 = -yinv * rho * w;
    flow_float s1 = -yinv * rho * u * w;
    flow_float s2 = -yinv * rho * w * w;
    flow_float s3 = -yinv * rho * w * H;

    if (viscous) {
        const flow_float mu_tot = vis_lam[ic] + max(vis_turb[ic], (flow_float)0.0);
        const flow_float k_tot  = thermCond[ic] + cp_arr[ic] * max(vis_turb[ic], (flow_float)0.0) / Prt;
        const flow_float TWO3 = (flow_float)(2.0/3.0);
        // SU2 ResidualDiffusion (residual -= 粘性項 → forge S += 粘性項)。
        // ※ SU2 の +ρk (turb_ke) 項は次元不整合の疑いがあり移植しない (plan §2)。
        s1 += yinv * mu_tot * (dUxdy[ic] + dUydx[ic]) - TWO3 * dA0dx[ic];
        s2 += yinv * mu_tot * (flow_float)2.0 * (dUydy[ic] - w * yinv) - TWO3 * dA0dy[ic];
        s3 += yinv * (mu_tot * (u * (dUydx[ic] + dUxdy[ic])
                                + w * TWO3 * ((flow_float)2.0 * dUydy[ic] - dUxdx[ic] - w * yinv))
                      + k_tot * dTdy[ic])
              - TWO3 * (dA1dy[ic] + dA2dx[ic]);
    }

    res_ro[ic]   += s0 * v;
    res_roUx[ic] += s1 * v;
    res_roUy[ic] += s2 * v;
    res_roe[ic]  += s3 * v;
}

void axisymmetricSourceSU2_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.isAxisymmetric != 1 || cfg.axisymMethod != 1) return;
    geom_int* axis_flag = (cfg.discretization == "node") ? msh.axis_flag_d : nullptr;
    // 診断 (env): FORGE_DIAG_SU2VISC_OFF=1 で粘性軸対称ソースを落とす (非粘性のみ)。
    static const bool diagViscOff = (getenv("FORGE_DIAG_SU2VISC_OFF") != nullptr);
    const int viscous = (cfg.viscMethod != 0 && !diagViscOff) ? 1 : 0;
    if (viscous) {
        axisymAuxFields_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
            msh.nCells,
            var.c_d["Ux"], var.c_d["Uy"],
            var.c_d["vis_lam"], var.c_d["vis_turb"],
            var.c_d["ccy"], axis_flag,
            var.c_d["axisym_aux0"], var.c_d["axisym_aux1"], var.c_d["axisym_aux2"]);
        gpuErrchk( cudaPeekAtLastError() );
        axisymAuxGradGG_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
            msh.nCells,
            msh.map_plane_cells_d, msh.map_cell_planes_index_d, msh.map_cell_planes_d,
            var.p_d["sx"], var.p_d["sy"], var.c_d["volume"],
            var.c_d["axisym_aux0"], var.c_d["axisym_aux1"], var.c_d["axisym_aux2"],
            var.c_d["dAux0dx"], var.c_d["dAux0dy"], var.c_d["dAux1dy"], var.c_d["dAux2dx"]);
        gpuErrchk( cudaPeekAtLastError() );
    }
    axisymmetricSourceSU2_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells,
        (msh.volumePartial_d != nullptr) ? msh.volumePartial_d : var.c_d["volume"],
        var.c_d["ccy"], axis_flag,
        var.c_d["ro"], var.c_d["Ux"], var.c_d["Uy"],
        var.c_d["roe"], var.c_d["P"],
        var.c_d["vis_lam"], var.c_d["vis_turb"], var.c_d["thermCond"], var.c_d["cp"],
        cfg.turbulentPrandtl,
        var.c_d["dUxdx"], var.c_d["dUxdy"], var.c_d["dUydx"], var.c_d["dUydy"],
        var.c_d["dTdy"],
        var.c_d["dAux0dx"], var.c_d["dAux0dy"], var.c_d["dAux1dy"], var.c_d["dAux2dx"],
        viscous,
        var.c_d["res_ro"], var.c_d["res_roUx"], var.c_d["res_roUy"], var.c_d["res_roe"],
        (cfg.sstEnergyIncludesK != 0 && cfg.LESorRANS == 2 && cfg.RANSmodel == 1) ? var.c_d["roK"] : nullptr,
        (cfg.sstEnergyIncludesK != 0 && cfg.LESorRANS == 2 && cfg.RANSmodel == 1) ? 1 : 0);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

// node-centered 軸対称: 軸上 CV (R=0) で半径方向運動量を 0 に (軸対称の対称条件)。
// エネルギー整合のため、除去する半径方向運動エネルギー 0.5*roUy^2/ro を roe からも引く
// (これを怠ると内部エネルギー=圧力が過大になり発散する)。
__global__ void enforceAxisSymmetry_d
(
    geom_int nCells, geom_int* axis_flag,
    flow_float* ro, flow_float* roUy, flow_float* roe, flow_float* Uy,
    // 擬似時間ステップ開始値も同時にピン (implicit commit ro=roN+dq が stale roUyN を巻き戻さないため)。nullptr 可。
    flow_float* roUyN, flow_float* roeN, flow_float* roN
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells && axis_flag[ic] == 1) {
        const flow_float r = ro[ic];
        if (r > (flow_float)0.0) {
            roe[ic] -= (flow_float)0.5 * roUy[ic] * roUy[ic] / r; // 半径方向 KE を除去
        }
        roUy[ic] = (flow_float)0.0;
        Uy[ic]   = (flow_float)0.0;
        if (roUyN != nullptr) {
            const flow_float rN = roN[ic];
            if (rN > (flow_float)0.0) roeN[ic] -= (flow_float)0.5 * roUyN[ic] * roUyN[ic] / rN;
            roUyN[ic] = (flow_float)0.0;
        }
    }
}

// SU2 流の軸対称対称面処理 (MARKER_SYM): 軸上 CV (R=0) で半径方向運動量の「残差」を 0 にする。
// res_roUy=0 にすると roUy は更新されず初期値 (=0, 対称条件) を保つ。状態やエネルギーをいじらず、
// flux+source を全て積んだ後に残差を射影するだけなので陰解法とも整合し発散しない (状態強制版・source
// 抑制版は発散したのに対し、これは SU2 と同じ「残差から法線運動量成分を除去」する方式)。
__global__ void zeroAxisRadialResidual_d
(
    geom_int nCells, geom_int* axis_flag, flow_float* res_roUy
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells && axis_flag[ic] == 1) {
        res_roUy[ic] = (flow_float)0.0;
    }
}

void zeroAxisRadialResidual_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.isAxisymmetric != 1 || cfg.discretization != "node" || msh.axis_flag_d == nullptr) return;
    zeroAxisRadialResidual_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells, msh.axis_flag_d, var.c_d["res_roUy"]
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void enforceAxisSymmetry_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    // node-centered 軸対称で常時作用 (2026-08-16 常時 ON 化)。
    // 呼び出し位置は assembleResidual 冒頭 (壁 no-slip の enforceWallNoSlip と同相): 更新後状態を毎ステージ
    // u_r=0 へ射影し、以降の派生量・勾配・flux が全てピン後状態を見る。commit 後の外部手術ではない
    // (旧 open issue の発散は commit 後強制 + Jacobian 非整合の組合せ。今回は roUy 行 decouple と対)。
    if (cfg.isAxisymmetric != 1 || cfg.discretization != "node" || msh.axis_flag_d == nullptr) return;

    enforceAxisSymmetry_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells, msh.axis_flag_d,
        var.c_d["ro"], var.c_d["roUy"], var.c_d["roe"], var.c_d["Uy"],
        var.c_d.count("roUyN") ? var.c_d["roUyN"] : nullptr,
        var.c_d.count("roeN") ? var.c_d["roeN"] : nullptr,
        var.c_d.count("roN") ? var.c_d["roN"] : nullptr
    );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void axisymmetricGeomTerms_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (cfg.isAxisymmetric != 1) return;

    axisymmetricGeomTerms_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>>(
        msh.nCells,
        var.c_d["Uy"],
        var.c_d["volume"],
        var.c_d["A_planar"],
        var.c_d["dUxdx"],
        var.c_d["dUydy"],
        var.c_d["axisym_uy_over_r"],
        var.c_d["axisym_divU"],
        cfg.axisymMethod,
        var.c_d["ccy"],
        (cfg.discretization == "node") ? msh.axis_flag_d : nullptr,
        cfg.axisRFloor
    );

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}
