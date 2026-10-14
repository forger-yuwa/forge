#include "speciesTransport_d.cuh"

#include "scalarTransport_d.cuh"
#include "thermo_d.cuh"
#include "species_eos_coupling_d.cuh"   // 案C EOS クロス応答 (解析 JVP)
#include "passiveKernels_d.cuh"          // species_advection_faceY_d / passive_bounds_d / passive_diffusion_d (受動種と共用)
#include "passiveTransport_d.cuh"        // 受動種基盤 (本 TU で実装)
#include "periodicNode_d.cuh"            // node 周期の gather/mirror (化学種 DPLUR dq・EOS クロス項・受動種)
#include "calcGradient_d.cuh"             // スカラー勾配の LSQ 経路 (mesh.scalarGradient: lsq)
#include "passiveFct_d.cuh"              // dual-time 物理 step 末尾の保存的 FCT 補正 (§4.7)
#include "condensationTransport_d.cuh"   // condensationSource_d_wrapper (FCT の凍結ソース)
#include "condensationCorrReasons_d.cuh"   // COND_REASON_* / cond_atomic_max_double (再正規化の理由別監視)
#include "gasPhaseComposition_d.cuh"     // 凝縮 carrier の気相組成 (拡散係数の組成; plan condensation-two-phase-transport §4.1)
#include "twoPhaseDiffusion_d.cuh"      // 二相拡散の面流束 (plan condensation-two-phase-transport §4.2, #4e)
#include "twoPhaseFaceDiag_d.cuh"       // 診断 D1 の面の中間量 (FORGE_DIAG_TP_FACES; plan condensation-two-phase-default #4)
#include "renormGate_d.cuh"              // 再正規化の受入ゲートの計測 (#1b-pre)
#include "twoPhaseUpdateDiag_d.cuh"      // 診断 G3-b の記録スロット (FORGE_DIAG_TP_UPDATE; plan condensation-two-phase-default #4g3)
#include "twoPhaseOperatorDiag_d.cuh"    // 診断 G3-a の面の記録 (FORGE_DIAG_TP_OPERATOR; plan condensation-two-phase-default #4g3・#4pjg)
#include "commitLossDiag.hpp"            // 診断 FORGE_DIAG_COMMIT_LOSS (plan architecture-float-state-double-geometry §4.5、既定 off)

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <array>
#include <string>
#include <utility>
#include <vector>

// 再正規化の受入ゲート (#1b-pre): 更新ごとの集計バッファと履歴 (定義は本ファイル末尾)
static double* rngBegin();                 // 1 更新の集計を 0 にして返す (履歴の容量も確保)
static void    rngEnd(cudaConfig& cuda_cfg); // 1 更新の集計を履歴の 1 行にする
static flow_float** rngMomentPtrs();       // 液・Q の device ポインタ配列 (順序 g, Q2, Q1, Q0; 無ければ nullptr)

namespace {

// device の roY / res_roY / transport_diag ポインタ配列 (flow_float*[nSpecies])。
// speciesInit_d で 1 度だけ構築 (kinetic 拡散カーネルが化学種ループするため)。
flow_float** g_roY_dev = nullptr;
flow_float** g_roYN_dev = nullptr;   // roYN: species ステップ始点ベースライン (speciesUpdateOuter で roY に一致)
flow_float*  g_roPred_dev = nullptr;  // 案C 予測時点の ρ (dual-time では roN=ρ^n が固定なので commit の δρ 基準はこちら; chem e296f0d0)
geom_int     g_roPred_cap = 0;
flow_float** g_resroY_dev = nullptr;
flow_float** g_transdiag_dev = nullptr;
int          g_nSpecies = 0;

// face 整合再構成 (speciesFaceReconstruction==1) 用: Y{s} と セル勾配 ∇Y{s} のポインタ配列。
flow_float** g_Y_dev    = nullptr;
flow_float** g_dYdx_dev = nullptr;
flow_float** g_dYdy_dev = nullptr;
flow_float** g_dYdz_dev = nullptr;
flow_float** g_limiterY_dev = nullptr;   // ψ_Y[s] (Venkat on Y) ポインタ配列
// S3: convectiveFlux が書き出す upwind 再構成 face 組成 (layout [ip*nSpecies+s])。species 移流が同一面組成で読む。
flow_float*  g_Yface_dev = nullptr;
int          g_Yface_nPlanes = 0;

// 案C 用 scratch: dq_roY{s}_old のポインタ配列 (sweep の swap でポインタが変わるため毎回再構築) と
// セルごとの δp_Y [Pa]。speciesInit_d でなく案C 経路の初回に遅延確保。
flow_float** g_dqYold_dev = nullptr;   // device 上の flow_float*[nSpecies]
flow_float*  g_dpY_eos    = nullptr;   // device 上の flow_float[nCells_all]
int          g_dqYold_cap = 0;         // 確保済みポインタ数
geom_int     g_dpY_cap    = 0;         // 確保済みセル数

constexpr flow_float kSmall = static_cast<flow_float>(1.0e-30);

inline bool speciesEnabled(const variables& var)
{
    return var.nSpeciesRegistered >= 2;
}

// 化学種 s 用のスカラ輸送記述子を構築する (RANS buildScalarDescs と同形)。
// floor=0 (Y_s>=0)、sigma=0 (M2 は移流のみ。拡散は M4)。
ScalarTransportDesc buildSpeciesDesc(variables& var, int s)
{
    const std::string i = std::to_string(s);
    return ScalarTransportDesc{
        var.c_d["Y"+i], nullptr, nullptr, nullptr,  // dphidx/y/z: 化学種は汎用拡散(diffusion=0)未使用
        var.c_d["roY"+i], var.c_d["roY"+i+"N"], var.c_d["roY"+i+"M"],
        var.c_d["res_roY"+i], var.c_d["res_roY"+i+"_m"],
        var.c_d["src_jac_Y"+i], var.c_d["transport_diag_Y"+i],
        static_cast<flow_float>(0.0), static_cast<flow_float>(0.0),
        0  // 汎用拡散は使わない (化学種 Fick 拡散は species_diffusion_d で別途)
    };
}

__global__ void species_primitive_d(
    geom_int nCells_all,
    flow_float* ro,
    flow_float* roY,
    flow_float* Y)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells_all) {
        Y[ic] = roY[ic] / max(ro[ic], kSmall);
    }
}

// Neumann (zero-gradient) ghost 充填: roY[ig]=roY[ic], Y[ig]=Y[ic]。
__global__ void species_neumann_boundary_d(
    geom_int nb,
    geom_int* bplane_cell,
    geom_int* bplane_cell_ghst,
    flow_float* roY,
    flow_float* Y)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib < nb) {
        const geom_int ic = bplane_cell[ib];
        const geom_int ig = bplane_cell_ghst[ib];
        roY[ig] = roY[ic];
        Y[ig]   = Y[ic];
    }
}

// Dirichlet ghost 充填 (M5: 組成依存超音速入口): Y[ig]=Y_s^in, roY[ig]=ρ[ig]·Y_s^in。
// ρ[ig] は applyBconds (inlet カーネル) が設定済みの ghost 密度。Yb は per-face 入口組成 bvar。
__global__ void species_dirichlet_boundary_d(
    geom_int nb,
    geom_int* bplane_cell,
    geom_int* bplane_cell_ghst,
    flow_float* ro,
    flow_float* Yb,
    flow_float* roY,
    flow_float* Y,
    flow_float* scalarDirichletPin,
    int isNode)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib < nb) {
        const geom_int ig = bplane_cell_ghst[ib];
        const flow_float Yin = max(Yb[ib], static_cast<flow_float>(0.0));   // 入口カーネルの負値クリップと整合
        Y[ig]   = Yin;
        roY[ig] = ro[ig] * Yin;
        // node-centered: 境界半割面の化学種流束は ghost を読まず境界ノード自身の組成を使う
        // (species_advection_faceY_d / scalarTransport の nodeBnd 分岐)。ghost だけ書いても入口ノードの
        // Y_s は初期値のまま凍結し、入口分布 (inletProfile の Y{s} 列) も一様値の変更も場に入らなかった
        // (2026-09-08, case/16 run_0317)。k/ω (rans_dirichlet_scalar_boundary_d) と同じく境界ノードを
        // Dirichlet 値にピンし、scalarDirichletPin で残差を除外する (speciesPinResidual_d_wrapper)。
        if (isNode != 0) {
            const geom_int ic = bplane_cell[ib];
            Y[ic]   = Yin;
            roY[ic] = ro[ic] * Yin;
            scalarDirichletPin[ic] = static_cast<flow_float>(1.0);
        }
    }
}

// node 入口ピン (scalarDirichletPin==1) のノードで化学種残差・ソース Jacobian を 0 化する
// (ransSource の k/ω 残差除外と同形)。移流残差と化学ソースの集計後に呼ぶ。
__global__ void species_pin_residual_d(
    geom_int nCells,
    int nSpecies,
    flow_float** res_roY,
    flow_float** src_jac,
    flow_float* scalarDirichletPin)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells && scalarDirichletPin[ic] == static_cast<flow_float>(1.0)) {
        for (int s = 0; s < nSpecies; ++s) {
            res_roY[s][ic] = static_cast<flow_float>(0.0);
            src_jac[s][ic] = static_cast<flow_float>(0.0);
        }
    }
}

// 実現可能性 + 再正規化: 各 ρY_s>=0 にクランプ後、Σ_s ρY_s = ρ となるよう再スケール (ΣY_s=1)。
// reasons (nullptr 可): 凝縮 carrier の理由別監視 (plan condensation-two-phase-transport §4.3)。凝縮種 iw の ρY_w に掛かった補正
//   Σ|Δ(ρY_w)|V・max|係数−1|・変化ノード数を計上する (root のみ; 算術と書き込み値は変えない)。
__global__ void species_renormalize_d(
    geom_int nCells,
    int nSpecies,
    flow_float** roY,
    flow_float* ro,
    double* reasons, int iw, const geom_float* vol, const geom_int* root,
    double* rnAcc, flow_float** rophiRO,   // 受入ゲートの計測 (#1b-pre; nullptr で計測しない)。rophiRO: 液・Q (読むだけ; この経路は液に掛けない)
    double* tpu = nullptr, int iwTpu = -1)   // 診断 G3-b (既定 nullptr): 水 (iwTpu) の負値の 0 化と係数の乗算を before/after で記録
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const double qw_m = (rnAcc != nullptr && iw >= 0) ? (double)roY[iw][ic] : 0.0;   // 再正規化の直前 (負値の 0 化の前)
        const flow_float tpu_w0 = (tpu != nullptr && iwTpu >= 0) ? roY[iwTpu][ic] : (flow_float)0.0;
        double sum = 0.0;
        for (int s = 0; s < nSpecies; s++) {
            flow_float v = roY[s][ic];
            if (v < 0.0) v = 0.0;
            roY[s][ic] = v;
            sum += (double)v;
        }
        const double factor = (double)ro[ic] / (sum > (double)kSmall ? sum : (double)kSmall);
        const flow_float yw_in = (reasons != nullptr && iw >= 0) ? roY[iw][ic] : (flow_float)0.0;   // クランプ前の値ではなく ≥0 化後
        const flow_float tpu_w1 = (tpu != nullptr && iwTpu >= 0) ? roY[iwTpu][ic] : (flow_float)0.0;
        for (int s = 0; s < nSpecies; s++) {
            roY[s][ic] = (flow_float)((double)roY[s][ic] * factor);
        }
        if (tpu != nullptr && iwTpu >= 0) {
            tpu_ba(tpu, TPU_OP_RNF, 0, nCells, ic, (double)tpu_w0, (double)tpu_w1);
            tpu_ba(tpu, TPU_OP_RNS, 0, nCells, ic, (double)tpu_w1, (double)roY[iwTpu][ic]);
            tpu_put(tpu, TPU_RN_FACTOR, nCells, ic, factor);
        }
        if (reasons != nullptr && iw >= 0 && (root == nullptr || root[ic] == ic)) {
            const double d = fabs((double)roY[iw][ic] - (double)yw_in);
            cond_atomic_max_double(&reasons[COND_REASON_RN_MAX], fabs(factor - 1.0));
            cond_atomic_max_double(&reasons[COND_REASON_RN_MAXC], fabs(factor - 1.0));   // 累積 (#1b-pre)
            if (rnAcc != nullptr && rophiRO != nullptr) {
                double qm[RNG_NC], qp[RNG_NC];
                qm[0] = qw_m; qp[0] = (double)roY[iw][ic];
                for (int m = 0; m < RNG_NC - 1; ++m) { qm[1+m] = (double)rophiRO[m][ic]; qp[1+m] = qm[1+m]; }   // 液・Q は変えない
                rng_accumulate(rnAcc, qm, qp, (vol != nullptr) ? (double)vol[ic] : 1.0, fabs(factor - 1.0));
            }
            if (d > 0.0) {
                atomicAdd(&reasons[COND_REASON_RN_SUM], d*((vol != nullptr) ? (double)vol[ic] : 1.0));
                atomicAdd(&reasons[COND_REASON_RN_N], 1.0);
            }
        }
    }
}

// 化学種拡散係数に渡す組成 (plan condensation-two-phase-transport §4.1): 総組成 Yf (正規化済み) から、同じ正規化の単位の
//   液 gl (= g/ΣY) を凝縮種から除いた気相組成のモル分率 X (gas_phase_composition; X の正規化が 1−g の再正規化を兼ねる)。
//   液なし (iw < 0 / gl ≤ 0) は Yf の X そのもの (現行とビット一致)。流束の駆動勾配と補正の Yf は変えない (§4.2 は別項目)。
__device__ inline void species_transport_X_f(const SpeciesThermoF* sp, int n, const flow_float* Yf, int iw, flow_float gl, flow_float* X)
{
    if (iw >= 0 && gl > 0.0f) {
        flow_float Yg[THERMO_MAX_SPECIES];
        for (int s = 0; s < n; ++s) Yg[s] = Yf[s];
        gas_phase_composition(Yg, iw, gl);
        thermo_X_from_Y_f(sp, n, Yg, X);
    } else {
        thermo_X_from_Y_f(sp, n, Yf, X);
    }
}

// M4: 化学種 Fick 拡散 + 質量保存補正 (ΣJ=0) + エンタルピー拡散のエネルギー結合。
//   J_s = ρ D_s ∇Y_s (over-relaxed 法線), 補正 J_s* = J_s - Y_s Σ_k J_k (Σ J_s*=0),
//   res_roY_s += J_s*, res_roe += Σ_s h_s(T) J_s*。
//   D_s: diffMethod==1 で混合平均 (Chapman-Enskog), ==0 で定数 Schmidt D=μ/(ρ Sc)。
__global__ void species_diffusion_d(
    geom_int nCells,
    geom_int nNormalHaloPlanes,
    geom_int* normal_halo_planes,
    geom_int* plane_cells,
    // 面ごとの差 e = cc[ic1] − cc[ic0] (double の座標から 1 回だけ丸めた値、var.p_d["ge_*"]、
    // plans/active/architecture-float-state-double-geometry.md §4.2a、段 ③)
    const flow_float* ge_x, const flow_float* ge_y, const flow_float* ge_z,
    geom_float* fx, geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,
    const SpeciesThermoF* sp, int nSpecies,
    flow_float** roY, flow_float** res_roY, flow_float** transport_diag,
    flow_float* ro, flow_float* T, flow_float* P, flow_float* vis_lam, flow_float* vis_turb,
    flow_float* res_roe,
    int diffMethod, flow_float Sc, flow_float Sc_t,
    int isNode, flow_float** dYdx, flow_float** dYdy, flow_float** dYdz,
    GasPhaseLiquid liq,   // 凝縮 carrier の液 (拡散係数は気相組成で評価; 液なしは {nullptr, -1})
    TpoFacePtr tpo = TpoFacePtr{})   // 診断 G3-a (既定は無効): res_roY[s][ic0] に足した Jc を面スロットへ写す
{
    // 面ループは float32 で評価する (係数は SpeciesThermoF, 評価点は従来どおり面状態 T_f/P_f/Y_f。
    // 離散式は不変, plan performance-3d-node-sst-speedup §4.2-2)。旧 double 版は FP64 パイプ律速で
    // 3D 2.37 M 節点 13.6 ms/step を食っていた。
    geom_int ih = blockDim.x * blockIdx.x + threadIdx.x;
    if (ih >= nNormalHaloPlanes) return;

    const geom_int ip  = normal_halo_planes[ih];
    const geom_int ic0 = plane_cells[2 * ip + 0];
    const geom_int ic1 = plane_cells[2 * ip + 1];

    const flow_float f   = fx[ip];
    const flow_float sxx = sx[ip], syy = sy[ip], szz = sz[ip], sss = ss[ip];

    // node モードの境界半割面 (ghost を含む面): 拡散流束を加えない (skip)。
    // 根拠 (plan diffusion-node-boundary-real-distance.md §3 (c), scalar_diffusion_first_order_d と同方針):
    //  - Dirichlet (固定組成入口) では境界ノード Y はピンで上書きされ無意味、内部ノードへの拡散は
    //    内部双対面 W↔I が実距離で運ぶ (本ループの非境界 plane で計算済)。
    //  - Neumann (壁 zero-grad / slip) は ∂Y/∂n=0 ＝半割面フラックス 0。
    //  → ghost mirror の dcc≈0 退化も ∇Y·S 弱形式の境界閉包依存も不要。エネルギー結合 (Σh_sJ_s) も
    //    半割面では 0。cell は ghost で正しく閉じるので従来どおり。
    if (isNode != 0 && (ic0 >= nCells || ic1 >= nCells)) {
        if (tpo.val != nullptr) {   // 診断 G3-a: 足さない面も明示的に 0 を書く
            tpo_face_code(tpo, ih, 2);
            for (int s = 0; s < nSpecies; ++s) tpo_face_put(tpo, ih, tpo_comp(tpo, s), 0.0);
        }
        return;
    }

    // 座標の差 cc[ic1] − cc[ic0] の代わりに面ごとの差 e (同じ向き) を読む (§4.2a、段 ③)
    const flow_float dccx = ge_x[ip];
    const flow_float dccy = ge_y[ip];
    const flow_float dccz = ge_z[ip];
    const flow_float dcc  = sqrtf(dccx*dccx + dccy*dccy + dccz*dccz);
    const flow_float denom = dccx*sxx + dccy*syy + dccz*szz;
    const flow_float Dsafe = (fabsf(denom) < 1.0e-30f) ? ((denom>=0.0f)?1.0e-30f:-1.0e-30f) : denom;
    const flow_float delta = dcc * sss * sss / Dsafe;       // over-relaxed 法線

    const flow_float ro0 = max(ro[ic0], (flow_float)1.0e-30f);
    const flow_float ro1 = max(ro[ic1], (flow_float)1.0e-30f);
    const flow_float inv_ro0 = 1.0f/ro0, inv_ro1 = 1.0f/ro1;
    const flow_float g = 1.0f - f;
    const flow_float ro_face = f*ro0 + g*ro1;
    const flow_float T_face  = f*T[ic0] + g*T[ic1];
    const flow_float P_face  = f*P[ic0] + g*P[ic1];

    // 面の質量分率 Yf (正規化) -> モル分率 X
    flow_float Yf[THERMO_MAX_SPECIES], X[THERMO_MAX_SPECIES];
    flow_float Ys0[THERMO_MAX_SPECIES], Ys1[THERMO_MAX_SPECIES];
    flow_float ysum = 0.0f;
    for (int s=0;s<nSpecies;s++){
        Ys0[s] = roY[s][ic0]*inv_ro0;
        Ys1[s] = roY[s][ic1]*inv_ro1;
        flow_float y = f*Ys0[s] + g*Ys1[s];
        if (y<0.0f) y=0.0f; Yf[s]=y; ysum+=y;
    }
    const flow_float yinv = 1.0f/(ysum>1.0e-30f?ysum:1.0e-30f);
    for (int s=0;s<nSpecies;s++) Yf[s]*=yinv;
    // 拡散係数の組成は気相組成 (面の液 g_f = f g0 + (1−f) g1 を Yf と同じ正規化の単位で)
    const flow_float gl_f = (liq.iw >= 0 && liq.rog != nullptr)
        ? (f*liq.rog[ic0]*inv_ro0 + g*liq.rog[ic1]*inv_ro1)*yinv : 0.0f;
    species_transport_X_f(sp, nSpecies, Yf, liq.iw, gl_f, X);

    const flow_float mu_face  = f*vis_lam[ic0]  + g*vis_lam[ic1];
    const flow_float mut_face = f*vis_turb[ic0] + g*vis_turb[ic1];
    // 乱流化学種拡散 D_t = μ_t/(ρ Sc_t) は全種共通で加える。
    const flow_float Dt = (mut_face > 0.0f) ? mut_face/(ro_face*Sc_t) : 0.0f;
    const flow_float inv_dcc = 1.0f/dcc;
    const flow_float diag_geo = fabsf(delta) / max(dcc,(flow_float)1.0e-30f);

    // 各化学種の非補正 Fick flux J_s と Σ
    flow_float Js[THERMO_MAX_SPECIES];
    flow_float sumJ = 0.0f;
    for (int s=0;s<nSpecies;s++){
        flow_float D;
        if (diffMethod == 1) D = thermo_Dmix_species_f(sp, nSpecies, X, s, T_face, P_face);
        else                 D = mu_face/(ro_face*Sc);
        D += Dt;  // 層流 (Fick/Sc) + 乱流 (μ_t/Sc_t)
        const flow_float roD = ro_face * D;
        Js[s] = roD * ((Ys1[s] - Ys0[s])*inv_dcc) * delta;
        sumJ += Js[s];
        // point-implicit 拡散対角 (各セル ρ で正規化)
        const flow_float diag = roD * diag_geo;
        if (ic0 < nCells) atomicAdd(&transport_diag[s][ic0], diag*inv_ro0);
        if (ic1 < nCells) atomicAdd(&transport_diag[s][ic1], diag*inv_ro1);
    }

    // 補正 J_s* = J_s - Y_s Σ_k J_k (Σ J_s* = 0) と エネルギー結合 Σ h_s J_s*
    flow_float q = 0.0f;
    for (int s=0;s<nSpecies;s++){
        const flow_float Jc = Js[s] - Yf[s]*sumJ;
        if (tpo.val != nullptr) { tpo_face_code(tpo, ih, 1); tpo_face_put(tpo, ih, tpo_comp(tpo, s), (double)Jc); }
        if (ic0 < nCells) atomicAdd(&res_roY[s][ic0],  Jc);
        if (ic1 < nCells) atomicAdd(&res_roY[s][ic1], -Jc);
        const flow_float hs = thermo_h_mass_f(sp[s], T_face);   // NASA エンタルピー [J/kg] (datum 込み)
        q += hs * Jc;
    }
    if (ic0 < nCells) atomicAdd(&res_roe[ic0],  q);
    if (ic1 < nCells) atomicAdd(&res_roe[ic1], -q);
}

// TP 多成分気体の組成-エネルギー整合補正 (M-stagger 修正)。
// speciesTimeIntegration で roY が更新されたとき、roe を同じ温度が維持されるよう補正する:
//   roe[ic] += Σ_s (roY_s[ic] - roYN_s[ic]) * h_s(T[ic])
// 背景: roY が変化しても T 一定を保つには roe を Σ_s Δ(roY_s)*h_s(T) だけ変化させる必要がある
// (等温条件下でのエンタルピー差)。staggered update では roe が未補正のまま speciesPrimitive
// の Newton 反転に入るため、TP h_mix が大きく負の成分 (H2O: ~-13.4 MJ/kg) では T が不物理的
// にジャンプし発散を引き起こす。本カーネルはその inconsistency を 1 次で除去する。
// 適用条件: thermalMethod==2 (TP) かつ nSpecies>=2。単成分/CPG では no-op。
__global__ void species_energy_correction_kernel(
    geom_int nCells,
    int nSpecies,
    flow_float* const* roY_dev,
    flow_float* const* roYN_dev,
    const flow_float* T,
    flow_float* roe,
    const SpeciesThermo* sp)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    const double Tic = (double)T[ic];
    double droe = 0.0;
    for (int s = 0; s < nSpecies; s++) {
        const double droY_s = (double)roY_dev[s][ic] - (double)roYN_dev[s][ic];
        if (droY_s != 0.0) {
            droe += droY_s * thermo_h_mass(sp[s], Tic);
        }
    }
    roe[ic] += (flow_float)droe;
}

// 緩和整合 scalar-DPLUR の 1 Jacobi sweep (speciesImplicitCoupling==1)。
// 流れ scalar-DPLUR (timeIntegration_d.cu implicit_defect_correction_d) を化学種 1 本にミラーする。
//   D·δ(ρY_s) = res_roY_s + Σ_f offdiag_f·δ(ρY_s)_old[nbr]
//   D    = V/Δτ + transport_diag[ic]              (対角=流出質量流束/ρ + 粘性時 Fick 拡散対角)
//   offdiag_f = max(∓ṁ_f,0)/ρ_nbr                 (非対角=流入質量流束/ρ_nbr, 1次風上の凍結 Jacobian)
// massflux[ip]>0 は ic0→ic1。ic=ic0 の流入は ṁ<0、ic=ic1 の流入は ṁ>0。ρ_nbr は基準 roN(=ρⁿ) で
// transport_diag の ρ と整合させる (flow commit 後の ρ ではなく)。dq_old は呼び出し側で memset 済み。
__global__ void species_dplur_sweep_d(
    flow_float implicit_relax,
    flow_float dtScale,          // dt_local の倍率 (scalarCflMax; 1.0 で厳密に不変)
    flow_float* dt_local,
    geom_int nCells,
    geom_float* vol,
    geom_int* plane_cells,
    geom_int* cell_planes_index,
    geom_int* cell_planes,
    flow_float* massflux,
    flow_float* roN,
    flow_float* res_roY,
    flow_float* transport_diag,
    flow_float* src_jac,
    flow_float* dq_old,
    flow_float* dq_new,
    flow_float* scalarDirichletPin)   // node 入口ピン (==1 の行は δ(ρY)=0 に拘束。cell/非ピンは無効)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        // ピン行: 残差 0 だけでは隣接の δ(ρY) が Jacobi sweep で漏れ込み、第 2 sweep 以降に入口の補正が
        // 非ゼロになる (方式 2 では EOS クロス流束にも伝播)。行そのものを δ=0 に拘束する。
        if (scalarDirichletPin != nullptr && scalarDirichletPin[ic] == static_cast<flow_float>(1.0)) {
            dq_new[ic] = static_cast<flow_float>(0.0);
            return;
        }
        const flow_float dt_l = dt_local[ic] * dtScale;
        const geom_float v = vol[ic];

        flow_float neighbor = 0.0;
        const geom_int plane_begin = cell_planes_index[ic];
        const geom_int plane_end   = cell_planes_index[ic + 1];
        for (geom_int po = plane_begin; po < plane_end; ++po) {
            const geom_int ip  = cell_planes[po];
            const geom_int ic0 = plane_cells[2 * ip + 0];
            const geom_int ic1 = plane_cells[2 * ip + 1];
            const geom_int other_ic = (ic0 == ic) ? ic1 : ic0;
            const flow_float mdot = massflux[ip];
            // 流入 (other→ic) の質量流束。ic=ic0: ṁ<0 で流入、ic=ic1: ṁ>0 で流入。
            const flow_float inflow = (ic0 == ic) ? max(-mdot, static_cast<flow_float>(0.0))
                                                  : max( mdot, static_cast<flow_float>(0.0));
            if (other_ic < nCells) {
                const flow_float offdiag = inflow / max(roN[other_ic], static_cast<flow_float>(1.0e-30));
                neighbor += offdiag * dq_old[other_ic];
            }
        }

        // D = V/Δτ + V·src_jac + transport_diag (点陰的対角と同形。src_jac は化学種では 0)。
        const flow_float diag = max(
            static_cast<flow_float>(v / max(dt_l, static_cast<flow_float>(1.0e-30)))
                + static_cast<flow_float>(v) * src_jac[ic]
                + transport_diag[ic],
            static_cast<flow_float>(1.0e-30));

        dq_new[ic] = implicit_relax * (res_roY[ic] + neighbor) / diag;
    }
}

// node 周期用の 2 段 sweep (codex result M1): 近傍 (非対角) 寄与を独立バッファに組み、周期 group で合算 (root+=members → broadcast)
// してから解く。従来の 1 段 kernel は各ノード自身の cell_planes (部分 CV) だけで neighbor を組むため、残差・対角は合併済みなのに
// 非対角が片側だけになっていた。非周期では従来 kernel (ビット不変)。
__global__ void species_dplur_neighbor_d(
    geom_int nCells, geom_int* plane_cells, geom_int* cell_planes_index, geom_int* cell_planes,
    flow_float* massflux, flow_float* roN, flow_float* dq_old, flow_float* scalarDirichletPin, flow_float* nb_out)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        if (scalarDirichletPin != nullptr && scalarDirichletPin[ic] == static_cast<flow_float>(1.0)) { nb_out[ic] = 0.0; return; }
        flow_float neighbor = 0.0;
        const geom_int plane_begin = cell_planes_index[ic];
        const geom_int plane_end   = cell_planes_index[ic + 1];
        for (geom_int po = plane_begin; po < plane_end; ++po) {
            const geom_int ip  = cell_planes[po];
            const geom_int ic0 = plane_cells[2 * ip + 0];
            const geom_int ic1 = plane_cells[2 * ip + 1];
            const geom_int other_ic = (ic0 == ic) ? ic1 : ic0;
            const flow_float mdot = massflux[ip];
            const flow_float inflow = (ic0 == ic) ? max(-mdot, static_cast<flow_float>(0.0))
                                                  : max( mdot, static_cast<flow_float>(0.0));
            if (other_ic < nCells) {
                neighbor += (inflow / max(roN[other_ic], static_cast<flow_float>(1.0e-30))) * dq_old[other_ic];
            }
        }
        nb_out[ic] = neighbor;
    }
}
__global__ void species_dplur_solve_d(
    flow_float implicit_relax, flow_float dtScale, flow_float* dt_local, geom_int nCells, geom_float* vol,
    flow_float* res_roY, flow_float* transport_diag, flow_float* src_jac, flow_float* nb, flow_float* dq_new, flow_float* scalarDirichletPin)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        if (scalarDirichletPin != nullptr && scalarDirichletPin[ic] == static_cast<flow_float>(1.0)) { dq_new[ic] = 0.0; return; }
        const flow_float dt_l = dt_local[ic] * dtScale;
        const geom_float v = vol[ic];
        const flow_float diag = max(
            static_cast<flow_float>(v / max(dt_l, static_cast<flow_float>(1.0e-30)))
                + static_cast<flow_float>(v) * src_jac[ic] + transport_diag[ic],
            static_cast<flow_float>(1.0e-30));
        dq_new[ic] = implicit_relax * (res_roY[ic] + nb[ic]) / diag;
    }
}

// 緩和整合 scalar-DPLUR の commit: ρY_s = ρY_s^N + δ(ρY_s)。実現可能性フロア (ρY_s>=0)。
// Σ_s ρY_s = ρ の再正規化は呼び出し側 (speciesRenormalize_d_wrapper) が行う。
// tpu (診断 G3-b; 既定 nullptr): 水の commit の before・δ・格納前の候補 fl(ρY_N + δ)・after を記録する (演算は同じ式)。
__global__ void species_commit_correction_d(
    geom_int nCells,
    flow_float* roY,
    flow_float* roYN,
    flow_float* dq,
    double* tpu = nullptr)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float cand = roYN[ic] + dq[ic];
        roY[ic] = max(cand, static_cast<flow_float>(0.0));
        if (tpu != nullptr)
            tpu_commit(tpu, TPU_OP_SPC, 0, 0, nCells, ic, roYN[ic], (double)dq[ic], dq[ic], cand, roY[ic]);
    }
}

// ===== 案C: block-triangular roe↔roY coupling (EOS クロス応答) =====

// 組成接空間への射影: z_s = δz*_s - Y_s Σ_r δz*_r (Σ_s z_s = 0)。dq[s] を in-place で z_s に書換える。
// Y_s = ρY_s^N/ρ^N (現組成で補正を分配; 単純等配分でない)。
__global__ void species_eos_project_tangent_d(
    geom_int nCells, int nSpecies,
    flow_float** dq, flow_float** roYN, flow_float* roN)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const double roc = max((double)roN[ic], (double)kSmall);
        double sumdq = 0.0;
        for (int s = 0; s < nSpecies; ++s) sumdq += (double)dq[s][ic];
        for (int s = 0; s < nSpecies; ++s) {
            const double Ys = (double)roYN[s][ic] / roc;
            dq[s][ic] = (flow_float)((double)dq[s][ic] - Ys * sumdq);
        }
    }
}

// セルごとの解析 EOS-JVP δp_Y = Σ_s (∂p/∂(ρY_s)) z_s を評価して dpY[ic] に格納。
// Y は ρY_s^N/ρ^N の正規化, T は現状態, ρ は ρ^N, z は射影済み dq。
__global__ void species_eos_dp_cell_d(
    geom_int nCells, int nSpecies, const SpeciesThermo* sp,
    flow_float** dq, flow_float** roYN, flow_float* roN, flow_float* T,
    flow_float* dpY)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const double roc = max((double)roN[ic], (double)kSmall);
        double Y[THERMO_MAX_SPECIES];
        double z[THERMO_MAX_SPECIES];
        double ysum = 0.0;
        for (int s = 0; s < nSpecies; ++s) {
            double y = (double)roYN[s][ic] / roc;
            if (y < 0.0) y = 0.0;
            Y[s] = y; ysum += y;
            z[s] = (double)dq[s][ic];
        }
        const double inv = 1.0 / (ysum > (double)kSmall ? ysum : (double)kSmall);
        for (int s = 0; s < nSpecies; ++s) Y[s] *= inv;
        double dp = 0.0, dT = 0.0;
        species_eos_cross_response(sp, nSpecies, Y, (double)T[ic], roc, z, &dp, &dT);
        dpY[ic] = (flow_float)dp;
    }
}

// クロスエネルギー流束 A_QY δY を res_roe へ移項する。エネルギー流束 mdot·H に対し
// δH = δp_Y/ρ (保存 ρe・ρ 固定) を上流差分で載せる (対流流束 res_roe_temp と同符号規約)。
__global__ void species_eos_cross_flux_d(
    geom_int nNormalPlanes, geom_int* plane_cells,
    flow_float* massflux, flow_float* dpY, flow_float* roN,
    flow_float* res_roe)
{
    const geom_int ip = blockDim.x * blockIdx.x + threadIdx.x;
    if (ip < nNormalPlanes) {
        const geom_int ic0 = plane_cells[2*ip+0];
        const geom_int ic1 = plane_cells[2*ip+1];
        const flow_float mdot = massflux[ip];
        const flow_float dH0  = dpY[ic0] / max(roN[ic0], kSmall);
        const flow_float dH1  = dpY[ic1] / max(roN[ic1], kSmall);
        const flow_float cross = 0.5f*(mdot + fabsf(mdot))*dH0
                               + 0.5f*(mdot - fabsf(mdot))*dH1;
        atomicAdd(&res_roe[ic0], -cross);
        atomicAdd(&res_roe[ic1],  cross);
    }
}

// 案C 最終 commit: ρY_s = ρY_s^N + z_s + Y_s^N δρ (δρ=ρ-ρ^N)。Σ_s δ(ρY_s)=δρ を満たす。
__global__ void species_eos_final_commit_d(
    geom_int nCells, int nSpecies,
    flow_float** roY, flow_float** roYN, flow_float** dq,
    flow_float* ro, flow_float* roN)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float dro = ro[ic] - roN[ic];
        const flow_float roc = max(roN[ic], kSmall);
        for (int s = 0; s < nSpecies; ++s) {
            const flow_float Ys = roYN[s][ic] / roc;
            const flow_float v  = roYN[s][ic] + dq[s][ic] + Ys * dro;
            roY[s][ic] = max(v, static_cast<flow_float>(0.0));
        }
    }
}

}  // namespace

static flow_float** g_srcjac_dev = nullptr;   // src_jac_Y{s} の device ポインタ配列 (化学反応ソース用)

// 1 回の scalar-DPLUR sweep (化学種・受動種共通)。node 周期では近傍寄与を独立バッファで周期合算してから解く (M1)。
static void dplurSweepOnce(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                           flow_float relax, flow_float dtScale, flow_float* roRef,
                           flow_float* res, flow_float* td, flow_float* sj, flow_float* dq_old, flow_float* dq_new, flow_float* pin)
{
    if (periodicNodeActive(cfg, msh)) {
        static flow_float* s_nb = nullptr; static geom_int s_cap = 0;
        if (s_cap < msh.nCells_all) { if (s_nb) cudaFree(s_nb); gpuErrchk( cudaMalloc((void**)&s_nb, (size_t)msh.nCells_all*sizeof(flow_float)) ); s_cap = msh.nCells_all; }
        species_dplur_neighbor_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
            msh.nCells, msh.map_plane_cells_d, msh.map_cell_planes_index_d, msh.map_cell_planes_d,
            var.p_d["massflux"], roRef, dq_old, pin, s_nb);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
        periodicGatherArray_d_wrapper(cfg, cuda_cfg, msh, s_nb);   // 合併 CV の非対角 = 両側部分 CV の和
        species_dplur_solve_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
            relax, dtScale, var.c_d["dt_local"], msh.nCells, var.c_d["volume"], res, td, sj, s_nb, dq_new, pin);
    } else {
        species_dplur_sweep_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
            relax, dtScale, var.c_d["dt_local"], msh.nCells, var.c_d["volume"],
            msh.map_plane_cells_d, msh.map_cell_planes_index_d, msh.map_cell_planes_d,
            var.p_d["massflux"], roRef, res, td, sj, dq_old, dq_new, pin);
    }
}
flow_float** species_resroY_device_ptr() { return g_resroY_dev; }
flow_float** species_srcjac_device_ptr() { return g_srcjac_dev; }

void speciesInit_d(solverConfig& cfg, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var)) {
        g_roY_dev = nullptr;
        g_nSpecies = 0;
        return;
    }
    g_nSpecies = var.nSpeciesRegistered;

    std::vector<flow_float*> hroY(g_nSpecies), hroYN(g_nSpecies), hres(g_nSpecies), htd(g_nSpecies), hsj(g_nSpecies);
    for (int s = 0; s < g_nSpecies; s++) {
        const std::string i = std::to_string(s);
        hroY[s]  = var.c_d["roY"+i];
        hroYN[s] = var.c_d["roY"+i+"N"];
        hres[s]  = var.c_d["res_roY"+i];
        htd[s]   = var.c_d["transport_diag_Y"+i];
        hsj[s]   = var.c_d["src_jac_Y"+i];
    }
    const size_t pbytes = g_nSpecies*sizeof(flow_float*);
    gpuErrchk( cudaMalloc((void**)&g_roY_dev,      pbytes) );
    gpuErrchk( cudaMalloc((void**)&g_roYN_dev,     pbytes) );
    gpuErrchk( cudaMalloc((void**)&g_resroY_dev,   pbytes) );
    gpuErrchk( cudaMalloc((void**)&g_transdiag_dev, pbytes) );
    gpuErrchk( cudaMemcpy(g_roY_dev,       hroY.data(),  pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(g_roYN_dev,      hroYN.data(), pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(g_resroY_dev,    hres.data(),  pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(g_transdiag_dev, htd.data(),   pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMalloc((void**)&g_srcjac_dev, pbytes) );
    gpuErrchk( cudaMemcpy(g_srcjac_dev, hsj.data(), pbytes, cudaMemcpyHostToDevice) );

    // face 整合再構成用 Y / ∇Y ポインタ配列。
    std::vector<flow_float*> hY(g_nSpecies), hdx(g_nSpecies), hdy(g_nSpecies), hdz(g_nSpecies);
    for (int s = 0; s < g_nSpecies; s++) {
        const std::string i = std::to_string(s);
        hY[s]  = var.c_d["Y"+i];
        hdx[s] = var.c_d["dY"+i+"dx"]; hdy[s] = var.c_d["dY"+i+"dy"]; hdz[s] = var.c_d["dY"+i+"dz"];
    }
    gpuErrchk( cudaMalloc((void**)&g_Y_dev,    pbytes) ); gpuErrchk( cudaMemcpy(g_Y_dev,    hY.data(),  pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMalloc((void**)&g_dYdx_dev, pbytes) ); gpuErrchk( cudaMemcpy(g_dYdx_dev, hdx.data(), pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMalloc((void**)&g_dYdy_dev, pbytes) ); gpuErrchk( cudaMemcpy(g_dYdy_dev, hdy.data(), pbytes, cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMalloc((void**)&g_dYdz_dev, pbytes) ); gpuErrchk( cudaMemcpy(g_dYdz_dev, hdz.data(), pbytes, cudaMemcpyHostToDevice) );
    std::vector<flow_float*> hlim(g_nSpecies);
    for (int s = 0; s < g_nSpecies; s++) hlim[s] = var.c_d["limiter_Y"+std::to_string(s)];
    gpuErrchk( cudaMalloc((void**)&g_limiterY_dev, pbytes) ); gpuErrchk( cudaMemcpy(g_limiterY_dev, hlim.data(), pbytes, cudaMemcpyHostToDevice) );

    std::cout << "speciesInit_d: built device roY/roYN/res/diag/Y/dY/limiterY[] for nSpecies=" << g_nSpecies << "\n";
}

// face 整合再構成用アクセサ
flow_float** species_Y_device_ptr()    { return g_Y_dev; }
flow_float** species_dYdx_device_ptr() { return g_dYdx_dev; }
flow_float** species_dYdy_device_ptr() { return g_dYdy_dev; }
flow_float** species_dYdz_device_ptr() { return g_dYdz_dev; }
flow_float** species_limiterY_device_ptr() { return g_limiterY_dev; }

// S3: convectiveFlux 用 face 組成バッファを確保し device ポインタを返す。
flow_float* species_Yface_alloc(int nPlanes)
{
    if (g_nSpecies < 2) return nullptr;
    if (g_Yface_dev == nullptr || g_Yface_nPlanes < nPlanes) {
        if (g_Yface_dev) cudaFree(g_Yface_dev);
        gpuErrchk( cudaMalloc((void**)&g_Yface_dev, (size_t)nPlanes*g_nSpecies*sizeof(flow_float)) );
        g_Yface_nPlanes = nPlanes;
    }
    return g_Yface_dev;
}

// species_advection_faceY_d は passiveKernels_d.cuh (受動種と共用)。
void speciesAdvectionFaceY_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var) || g_Yface_dev == nullptr) return;
    dim3 dimGrid_nh = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));
    species_advection_faceY_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(
        msh.nCells, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
        var.c_d["ro"], var.p_d["massflux"], g_nSpecies, g_Yface_dev, g_resroY_dev, g_transdiag_dev,
        (cfg.discretization == "node") ? 1 : 0, g_roY_dev, g_nSpecies, farfieldFaceYDevice(),
        tpoFace(TPO_K_ADV_SP, TPO_KIND_ADV, -1));   // 診断 G3-a (既定は無効)
    gpuErrchk( cudaPeekAtLastError() );
}

// 化学種セル勾配 ∇Y{s} を Green-Gauss で計算する (calcGradient と同形)。speciesFaceReconstruction==1 のみ。
// 境界は Neumann ghost (applySpeciesBoundaries 済) を用い、内部面と同様に集計する。
// excludePeriodic (node 周期; plan species-passive-scalar-unification §4.1-5-1): 周期半割面 (面フラグ planePeriodic。
// 2026-09-26 まで「ip>=nNormalPlanes かつ相手が実 CV」で判定しており、周期にもゴーストが付くため一度も除外していなかった) を積算から除外し、勾配は内部双対面だけ (片側) にしておく。後段 periodicGradientGather が両側を合併体積で
// 厳密合併する (流れの calcGradient と同じ扱い)。0 のとき (cell / 非周期) は従来どおり全 plane。
__global__ void species_gradient_d(
    geom_int nCells, geom_int nPlanes, geom_int* plane_cells,
    geom_float* vol, geom_float* fx, geom_float* sx, geom_float* sy, geom_float* sz,
    int nSpecies, flow_float** Y, flow_float** dYdx, flow_float** dYdy, flow_float** dYdz,
    int excludePeriodic, const unsigned char* planePeriodic, flow_float* dumpFace)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;
    if (ip < nPlanes) {
        geom_int ic0 = plane_cells[2*ip+0];
        geom_int ic1 = plane_cells[2*ip+1];
        if (excludePeriodic != 0 && planePeriodic[ip] != 0) return;   // node 周期半割面 (面フラグで判定、gradient-fix §4.2a)
        geom_float f = fx[ip];
        const geom_float sxx = sx[ip], syy = sy[ip], szz = sz[ip];
        for (int s = 0; s < nSpecies; ++s) {
            const flow_float Yf = f*Y[s][ic0] + (1.0-f)*Y[s][ic1];
            atomicAdd(&dYdx[s][ic0],  sxx*Yf); atomicAdd(&dYdy[s][ic0],  syy*Yf); atomicAdd(&dYdz[s][ic0],  szz*Yf);
            if (dumpFace != nullptr) {   // 診断 (FORGE_DUMP_SCALARGRAD): atomicAdd に渡す同じ値を面ごとに非 atomic で書く
                const size_t o = 3*((size_t)s*nPlanes + ip);
                dumpFace[o+0] = sxx*Yf; dumpFace[o+1] = syy*Yf; dumpFace[o+2] = szz*Yf;
            }
            if (ic1 < nCells) { atomicAdd(&dYdx[s][ic1], -sxx*Yf); atomicAdd(&dYdy[s][ic1], -syy*Yf); atomicAdd(&dYdz[s][ic1], -szz*Yf); }
        }
    }
}

__global__ void species_gradient_normalize_d(
    geom_int nCells, geom_float* vol, int nSpecies,
    flow_float** dYdx, flow_float** dYdy, flow_float** dYdz)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float invv = 1.0/max(vol[ic], (flow_float)1.0e-30);
        for (int s = 0; s < nSpecies; ++s) {
            dYdx[s][ic] *= invv; dYdy[s][ic] *= invv; dYdz[s][ic] *= invv;
        }
    }
}

// 診断ダンプ (env `FORGE_DUMP_SCALARGRAD=<path>`、既定 off。**数値の振る舞いは変えない**)。
// species_gradient_d が atomicAdd に渡す面寄与 (sx,sy,sz)·φ_f を、化学種・受動種それぞれ**最初の呼び出しだけ**
// 面ごとの配列 [nVar][nPlanes][3] に非 atomic で書き、raw float で <path>.<tag> へ出す。除外した面は 0 のまま。
// 場 (res_*.h5) は atomicAdd の集積順序でビット再現しないが、この面寄与は 1 面 1 スレッドなので面レベルで比べられる
// (plan boundary-node-periodic-gradient-fix §5.1 #8a。先例は convectiveFlux_d.cu の FORGE_DUMP_MASSFLUX)。
static flow_float* scalarGradDumpBegin(const char* tag, bool& done, int nVar, geom_int nPlanes)
{
    if (done) return nullptr;
    const char* p = std::getenv("FORGE_DUMP_SCALARGRAD");
    if (!p || !*p || nVar <= 0) return nullptr;
    done = true;
    flow_float* d = nullptr;
    const size_t n = (size_t)3*nVar*nPlanes;
    gpuErrchk( cudaMalloc((void**)&d, n*sizeof(flow_float)) );
    gpuErrchk( cudaMemset(d, 0, n*sizeof(flow_float)) );
    (void)tag;
    return d;
}

static void scalarGradDumpEnd(const char* tag, flow_float* d, int nVar, geom_int nPlanes)
{
    if (d == nullptr) return;
    const size_t n = (size_t)3*nVar*nPlanes;
    std::vector<flow_float> h(n);
    gpuErrchk( cudaMemcpy(h.data(), d, n*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    cudaFree(d);
    const std::string path = std::string(std::getenv("FORGE_DUMP_SCALARGRAD")) + "." + tag;
    std::ofstream ofs(path, std::ios::binary);
    if (ofs) {
        ofs.write(reinterpret_cast<const char*>(h.data()), (std::streamsize)(n*sizeof(flow_float)));
        std::cout << "[FORGE_DUMP_SCALARGRAD] wrote " << nVar << " x " << nPlanes << " x 3 face contributions to " << path << '\n';
    } else {
        std::cout << "[FORGE_DUMP_SCALARGRAD] cannot open " << path << '\n';
    }
}

void speciesGradient_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;
    const int n = var.nSpeciesRegistered;
    for (int s = 0; s < n; ++s) {
        const std::string i = std::to_string(s);
        cudaMemset(var.c_d["dY"+i+"dx"], 0, msh.nCells*sizeof(flow_float));
        cudaMemset(var.c_d["dY"+i+"dy"], 0, msh.nCells*sizeof(flow_float));
        cudaMemset(var.c_d["dY"+i+"dz"], 0, msh.nCells*sizeof(flow_float));
    }
    // mesh.scalarGradient: lsq (node のみ) — NS と同じ事前計算 LSQ 係数の差分形 gather (plan gradient-scalar-lsq-unification §4.2)。
    // 体積除算 (normalize) はしない。周期の和→broadcast は periodicSeamMergeActive のときだけここで行い、
    // periodicGradientGather は lsq のとき dY を登録しない (二重合併の回避)。FORGE_DUMP_SCALARGRAD は GG 経路専用。
    if (scalarGradientLsqActive(cfg)) {
        lsqScalarGradient_d_wrapper(cuda_cfg, msh, n, g_Y_dev, g_dYdx_dev, g_dYdy_dev, g_dYdz_dev);
        if (preGatherDumpEnabled()) {   // 診断 (FORGE_DUMP_PREGATHER、既定 off・出力専用): 周期合併の直前の局所配列
            std::vector<std::string> nm; std::vector<std::array<const flow_float*, 3>> gp;
            for (int s = 0; s < n; ++s) {
                const std::string i = std::to_string(s);
                nm.push_back("Y" + i); gp.push_back({var.c_d["dY"+i+"dx"], var.c_d["dY"+i+"dy"], var.c_d["dY"+i+"dz"]});
            }
            preGatherDump("species_lsq.loop1", msh.nCells, nm, gp);
        }
        if (periodicSeamMergeActive(cfg, msh)) {
            for (int s = 0; s < n; ++s) {
                const std::string i = std::to_string(s);
                for (const char* c : {"x", "y", "z"}) periodicGatherArray_d_wrapper(cfg, cuda_cfg, msh, var.c_d["dY"+i+"d"+c]);
            }
        }
        return;
    }
    flow_float* gvol = (cfg.isAxisymmetric == 1) ? var.c_d["A_planar"] : var.c_d["volume"];
    flow_float* gsx = (cfg.isAxisymmetric == 1) ? var.p_d["sx_planar"] : var.p_d["sx"];
    flow_float* gsy = (cfg.isAxisymmetric == 1) ? var.p_d["sy_planar"] : var.p_d["sy"];
    flow_float* gsz = (cfg.isAxisymmetric == 1) ? var.p_d["sz_planar"] : var.p_d["sz"];
    static bool s_sgDumped = false;
    flow_float* sgDump = scalarGradDumpBegin("species", s_sgDumped, n, msh.nPlanes);
    species_gradient_d<<<cuda_cfg.dimGrid_plane, cuda_cfg.dimBlock>>>(
        msh.nCells, msh.nPlanes, msh.map_plane_cells_d, gvol, var.p_d["fx"], gsx, gsy, gsz,
        n, g_Y_dev, g_dYdx_dev, g_dYdy_dev, g_dYdz_dev,
        periodicSeamMergeActive(cfg, msh) ? 1 : 0, msh.planePeriodic_d, sgDump);
    scalarGradDumpEnd("species", sgDump, n, msh.nPlanes);
    species_gradient_normalize_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, gvol, n, g_dYdx_dev, g_dYdy_dev, g_dYdz_dev);
    gpuErrchk( cudaPeekAtLastError() );
}

flow_float** species_roY_device_ptr()
{
    return g_roY_dev;
}

void speciesEnergyCorrection_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (cfg.thermalMethod != 2) return;
    if (!speciesEnabled(var) || g_roY_dev == nullptr || g_roYN_dev == nullptr) return;

    dim3 dimGrid_c = dim3((geom_int)ceil(msh.nCells / (flow_float)cuda_cfg.blocksize));
    species_energy_correction_kernel<<<dimGrid_c, cuda_cfg.dimBlock>>>(
        msh.nCells,
        g_nSpecies,
        g_roY_dev,
        g_roYN_dev,
        var.c_d["T"],
        var.c_d["roe"],
        thermo_species_device_ptr());

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void speciesPrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var)) return;

    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        species_primitive_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
            msh.nCells_all,
            var.c_d["ro"],
            var.c_d["roY"+i],
            var.c_d["Y"+i]);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void speciesBoundary_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, bcond& bc, mesh& msh, variables& var)
{
    (void)msh;
    if (!speciesEnabled(var)) return;
    if (bc.iPlanes.empty()) return;

    const geom_int nb = static_cast<geom_int>(bc.iPlanes.size());
    // M5: 入口種別は組成依存 Dirichlet (roY[ig]=ρ[ig]·Y_s^in)。入口組成 bvar_d["Y{s}"] は
    // readBcondConfig が inlet_* に対して登録済み。他種別は従来通り Neumann (zero-gradient)。
    const bool isInlet = bc.bcondKind.rfind("inlet_", 0) == 0;
    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        const auto ybIt = bc.bvar_d.find("Y"+i);
        if (isInlet && ybIt != bc.bvar_d.end()) {
            species_dirichlet_boundary_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(
                nb,
                bc.map_bplane_cell_d,
                bc.map_bplane_cell_ghst_d,
                var.c_d["ro"],
                ybIt->second,
                var.c_d["roY"+i],
                var.c_d["Y"+i],
                var.c_d["scalarDirichletPin"],
                (cfg.discretization == "node") ? 1 : 0);
        } else {
            species_neumann_boundary_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(
                nb,
                bc.map_bplane_cell_d,
                bc.map_bplane_cell_ghst_d,
                var.c_d["roY"+i],
                var.c_d["Y"+i]);
        }
    }
}

void speciesPinResidual_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var) || cfg.discretization != "node" || g_resroY_dev == nullptr || g_srcjac_dev == nullptr) return;
    species_pin_residual_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, g_nSpecies, g_resroY_dev, g_srcjac_dev, var.c_d["scalarDirichletPin"]);
    gpuErrchk( cudaPeekAtLastError() );
}

void applySpeciesBoundaries(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;

    for (auto& bc : msh.bconds) {
        speciesBoundary_d_wrapper(cfg, cuda_cfg, bc, msh, var);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void speciesTransport_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;

    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["res_roY"+i], 0, msh.nCells * sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["transport_diag_Y"+i], 0, msh.nCells * sizeof(flow_float)));
        CHECK_CUDA_ERROR(cudaMemset(var.c_d["src_jac_Y"+i], 0, msh.nCells * sizeof(flow_float)));
    }
    tpoSnap("sp_zero");   // 診断 G3-a: 残差の写し (記録中のみ; 既定は何もしない)

    if (cfg.speciesFaceReconstruction >= 2 && g_Yface_dev != nullptr) {
        // S3: convectiveFlux が書いた同一 face 組成で移流 (energy 流束と整合)。diag は 1 次のまま。
        speciesAdvectionFaceY_d_wrapper(cfg, cuda_cfg, msh, var);
    } else {
        // 化学種の 1 次風上移流を最大 4 種ずつ 1 面ループで融合 (massflux・ρ の読みを共有)。
        std::vector<ScalarTransportDesc> descs;
        for (int s = 0; s < var.nSpeciesRegistered; s++) {
            descs.push_back(buildSpeciesDesc(var, s));
            descs.back().ext_face = farfieldFaceScalar("Y" + std::to_string(s));   // farfield 面の流入組成 (無ければ nullptr)
        }
        scalarTransportResidualMulti_d(cfg, cuda_cfg, msh, var, descs.data(), (int)descs.size());
    }
    tpoSnap("sp_adv");   // 診断 G3-a (記録中のみ)

    // M4: 粘性ケースのみ Fick 拡散 + ΣJ=0 補正 + エンタルピー拡散 (res_roe へ加算)。
    // 二相拡散 (condTwoPhaseDiffusion, #4e) が働く構成では水・気相種の分子拡散を twoPhaseDiffusion_d_wrapper (凝縮モーメントの
    // 残差初期化の後) が組むので、ここでは足さない (二重に足さない; plan condensation-two-phase-transport §4.2)。
    if (cfg.viscMethod != 0 && g_roY_dev != nullptr && !condTwoPhaseDiffusionActive(cfg)) {
        dim3 dimGrid_nh = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));
        species_diffusion_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(
            msh.nCells, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
            var.p_d["ge_x"], var.p_d["ge_y"], var.p_d["ge_z"],
            var.p_d["fx"], var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
            thermo_species_device_ptr_f(), g_nSpecies,
            g_roY_dev, g_resroY_dev, g_transdiag_dev,
            var.c_d["ro"], var.c_d["T"], var.c_d["P"], var.c_d["vis_lam"], var.c_d["vis_turb"],
            var.c_d["res_roe"],
            cfg.speciesDiffusionMethod, cfg.Sc, cfg.Sc_t,
            (cfg.discretization == "node") ? 1 : 0, g_dYdx_dev, g_dYdy_dev, g_dYdz_dev,
            gasPhaseLiquid(cfg, var),
            tpoFace(TPO_K_DIFF_OFF, TPO_KIND_DIFF, -1));   // 診断 G3-a (既定は無効)
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
    tpoSnap("sp_diff");   // 診断 G3-a (記録中のみ; 二相 ON では拡散を足さないので sp_adv と同じ)
}

// 試験用 (FORGE_TRANSPORT_PROBE; tests/unit/test_transport_gpu.py)。計算経路からは呼ばない。
//   セル ic を面の両側に置いた (f = 1) ときの species_diffusion_d と同じ組成 (clip → 正規化 → 気相組成 species_transport_X_f) と
//   同じ関数 thermo_Dmix_species_f で、分子拡散係数 D_s (乱流分を含まない) を全セル (ghost 込み) について書く。D は [s*nCells_all + ic]。
__global__ void species_Dmix_probe_d(geom_int nCells_all, const SpeciesThermoF* sp, int nSpecies, flow_float** roY,
                                     const flow_float* ro, const flow_float* T, const flow_float* P, GasPhaseLiquid liq, float* D)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells_all) return;
    const flow_float inv_ro = 1.0f/max(ro[ic], (flow_float)1.0e-30f);
    flow_float Yf[THERMO_MAX_SPECIES], X[THERMO_MAX_SPECIES];
    flow_float ysum = 0.0f;
    for (int s = 0; s < nSpecies; s++) { flow_float y = roY[s][ic]*inv_ro; if (y < 0.0f) y = 0.0f; Yf[s] = y; ysum += y; }
    const flow_float yinv = 1.0f/(ysum>1.0e-30f?ysum:1.0e-30f);
    for (int s = 0; s < nSpecies; s++) Yf[s] *= yinv;
    const flow_float gl = (liq.iw >= 0 && liq.rog != nullptr) ? liq.rog[ic]*inv_ro*yinv : 0.0f;
    species_transport_X_f(sp, nSpecies, Yf, liq.iw, gl, X);
    for (int s = 0; s < nSpecies; s++) D[(size_t)s*nCells_all + ic] = thermo_Dmix_species_f(sp, nSpecies, X, s, T[ic], P[ic]);
}

bool speciesDmixProbe_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, float* D_d)
{
    if (!speciesEnabled(var) || g_roY_dev == nullptr) return false;
    species_Dmix_probe_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells_all, thermo_species_device_ptr_f(), g_nSpecies, g_roY_dev,
        var.c_d["ro"], var.c_d["T"], var.c_d["P"], gasPhaseLiquid(cfg, var), D_d);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
    return true;
}

void speciesTimeIntegration_d_wrapper(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;

    // speciesImplicitRelax (既定 1.0 = 現行と同じ写像) と scalarCflMax は timeIntegration 11 の point-implicit 経路にだけ効く。
    const flow_float relax = static_cast<flow_float>(cfg.speciesImplicitRelax);
    const flow_float dts   = scalarDtScale(cfg);
    // 診断 FORGE_DIAG_COMMIT_LOSS (既定 off): 計る step だけ、point-implicit のカーネルに relax·δ の書き先を渡し、
    // commit の直後 (再正規化の前) に Q_before = roY{s}N・Q_after = roY{s} と比べる。off では nullptr のまま。
    const bool commitLoss = commitLossDiag::active() && cfg.timeIntegration == 11;
    double* diagDq = commitLoss ? commitLossDiag::speciesDqBuffer() : nullptr;
    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const ScalarTransportDesc desc = buildSpeciesDesc(var, s);
        scalarTimeIntegration_d(loop, cfg, cuda_cfg, msh, var, desc, relax, dts, diagDq);
        if (commitLoss)
            commitLossDiag::speciesCommitPointImplicit(msh, s, desc.rho_phi_N, desc.rho_phi, diagDq,
                                                       static_cast<double>(desc.floor));
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void speciesRenormalize_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var) || g_roY_dev == nullptr) return;

    // 診断 (FORGE_SPECIES_RAW_DIAG=1, registerSpecies が roYraw{s} を登録): 再正規化前の生更新値を退避する
    // (plan §6-1 の制御試験: 受動トレーサは再正規化を受けないので、比較対象は生更新値; §4.0)。
    for (int s = 0; s < g_nSpecies; ++s) {
        auto it = var.c_d.find("roYraw"+std::to_string(s));
        if (it == var.c_d.end() || it->second == nullptr) continue;
        gpuErrchk( cudaMemcpy(it->second, var.c_d["roY"+std::to_string(s)], (size_t)msh.nCells_all*sizeof(flow_float), cudaMemcpyDeviceToDevice) );
    }

    int iwRn = -1;
    double* reasonsRn = condCorrReasonsForRenormalize(cfg, var, &iwRn);   // 凝縮 carrier の理由別監視 (計上のみ)
    flow_float** mom = (reasonsRn != nullptr) ? rngMomentPtrs() : nullptr;
    double* rnAcc = (mom != nullptr) ? rngBegin() : nullptr;   // 受入ゲートの計測 (TP carrier 凝縮のみ; 計上だけ)
    species_renormalize_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells,
        g_nSpecies,
        g_roY_dev,
        var.c_d["ro"],
        reasonsRn, iwRn, var.c_d["volume"], passive_periodic_root(cfg, msh), rnAcc, mom,
        tpuSlots(), cfg.condGasSpecies);   // 診断 G3-b (既定 nullptr)
    if (rnAcc != nullptr) rngEnd(cuda_cfg);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

bool speciesImplicitCoupled(solverConfig& cfg, variables& var)
{
    return cfg.speciesImplicitCoupling == 1 && speciesEnabled(var);
}

// 緩和整合 scalar-DPLUR ソルバ (speciesImplicitCoupling==1)。
// 凍結 res_roY/transport_diag/massflux/dt_local (assembleResidual で確定) に対し、各化学種の補正
// δ(ρY_s) を dq=0 から nStepInner 回 Jacobi sweep で緩和し (流れ block と同一 implicitRelax/sweep 回数)、
// ρY_s = ρY_s^N + δ(ρY_s) を commit する。呼び出し前に speciesUpdateOuter で ρY_s^N=ρY_s を取ること。
// commit 後の Σ_s ρY_s=ρ 再正規化・Y=ρY/ρ 同期は呼び出し側 (main.cpp) で行う。
void speciesImplicitDPLURSolve_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;

    const int nSpecies = var.nSpeciesRegistered;
    const size_t bytes = static_cast<size_t>(msh.nCells_all) * sizeof(flow_float);

    // dq (new/old) を 0 初期化 (古典 DPLUR は dq=0 から開始)。
    for (int s = 0; s < nSpecies; s++) {
        const std::string i = std::to_string(s);
        cudaMemset(var.c_d["dq_roY"+i],        0, bytes);
        cudaMemset(var.c_d["dq_roY"+i+"_old"], 0, bytes);
    }

    const int nSweep = std::max(1, cfg.nStepInner);
    // 非対角の ρ_nbr: 定常は roN (=assemble 時の ρ, ビット不変)。dual-time では roN は物理レベル ρ^n で反復値と違うので現在の ρ。
    flow_float* roRef = (cfg.unsteady == 1 && cfg.dualTime == 1) ? var.c_d["ro"] : var.c_d["roN"];
    for (int iSweep = 0; iSweep < nSweep; ++iSweep) {
        for (int s = 0; s < nSpecies; s++) {
            const std::string i = std::to_string(s);
            dplurSweepOnce(cfg, cuda_cfg, msh, var, cfg.implicitRelax, scalarDtScale(cfg), roRef,
                var.c_d["res_roY"+i], var.c_d["transport_diag_Y"+i], var.c_d["src_jac_Y"+i],
                var.c_d["dq_roY"+i+"_old"], var.c_d["dq_roY"+i],
                (cfg.discretization == "node") ? var.c_d["scalarDirichletPin"] : nullptr);
        }
        // sweep 後に old↔new を swap (最終補正は dq_old 側に残る)。
        for (int s = 0; s < nSpecies; s++) {
            const std::string i = std::to_string(s);
            std::swap(var.c_d["dq_roY"+i], var.c_d["dq_roY"+i+"_old"]);
            // node 周期 DOF 同一視: 最新補正 dq_old を root→member でミラー (流れ側 periodicMirrorDq と同じ; §4.1-5)。
            periodicBroadcastArray_d_wrapper(cfg, cuda_cfg, msh, var.c_d["dq_roY"+i+"_old"]);
        }
    }

    // commit: ρY_s = ρY_s^N + δ(ρY_s) (最終補正は dq_roY{s}_old)。
    for (int s = 0; s < nSpecies; s++) {
        const std::string i = std::to_string(s);
        species_commit_correction_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
            msh.nCells,
            var.c_d["roY"+i],
            var.c_d["roY"+i+"N"],
            var.c_d["dq_roY"+i+"_old"],
            (s == cfg.condGasSpecies) ? tpuSlots() : nullptr);   // 診断 G3-b (既定 nullptr)
        // 診断 FORGE_DIAG_COMMIT_LOSS (既定 off): commit の直後に Q_before = roY{s}N・dq_req = dq_roY{s}_old・Q_after = roY{s}
        // (床 0 は species_commit_correction_d の max(·, 0))
        if (commitLossDiag::active())
            commitLossDiag::speciesCommitDPLUR(msh, s, var.c_d["roY"+i+"N"], var.c_d["roY"+i], var.c_d["dq_roY"+i+"_old"], 0.0);
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

// ===== 案C: block-triangular roe↔roY coupling のオーケストレーション =====

bool speciesEOSCoupled(solverConfig& cfg, variables& var)
{
    return cfg.speciesImplicitCoupling == 2 && speciesEnabled(var);
}

// dq_roY{s}_old のポインタ配列を (再)構築する。sweep の swap でポインタが変わるため毎回呼ぶ。
static void rebuildDqYoldPtrs(variables& var, int nSpecies)
{
    if (g_dqYold_cap < nSpecies) {
        if (g_dqYold_dev) cudaFree(g_dqYold_dev);
        gpuErrchk( cudaMalloc((void**)&g_dqYold_dev, nSpecies*sizeof(flow_float*)) );
        g_dqYold_cap = nSpecies;
    }
    std::vector<flow_float*> h(nSpecies);
    for (int s = 0; s < nSpecies; ++s) h[s] = var.c_d["dq_roY"+std::to_string(s)+"_old"];
    gpuErrchk( cudaMemcpy(g_dqYold_dev, h.data(), nSpecies*sizeof(flow_float*), cudaMemcpyHostToDevice) );
}

void speciesEOSCrossPredictInject_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var)) return;
    const int nSpecies = var.nSpeciesRegistered;
    const size_t bytes = static_cast<size_t>(msh.nCells_all) * sizeof(flow_float);

    // 予測時点の ρ を保存 (commit の δρ = ρ_after − ρ_predict)。定常陰解法では roN==ro なので従来と同一値 (ビット不変)、
    // dual-time では roN=ρ^n (物理レベル, サブ反復中は固定) を基準にすると δρ が累積/欠落して ΣρY≠ρ になる (chem e296f0d0)。
    if (g_roPred_cap < msh.nCells_all) {
        if (g_roPred_dev) cudaFree(g_roPred_dev);
        gpuErrchk( cudaMalloc((void**)&g_roPred_dev, bytes) );
        g_roPred_cap = msh.nCells_all;
    }
    gpuErrchk( cudaMemcpy(g_roPred_dev, var.c_d["ro"], bytes, cudaMemcpyDeviceToDevice) );

    // --- species scalar-DPLUR sweep で δ(ρY_s)* を予測 (dq=0 から, commit しない) ---
    for (int s = 0; s < nSpecies; ++s) {
        const std::string i = std::to_string(s);
        cudaMemset(var.c_d["dq_roY"+i],        0, bytes);
        cudaMemset(var.c_d["dq_roY"+i+"_old"], 0, bytes);
    }
    const int nSweep = std::max(1, cfg.nStepInner);
    for (int iSweep = 0; iSweep < nSweep; ++iSweep) {
        for (int s = 0; s < nSpecies; ++s) {
            const std::string i = std::to_string(s);
            dplurSweepOnce(cfg, cuda_cfg, msh, var, cfg.implicitRelax, scalarDtScale(cfg),
                (cfg.unsteady == 1 && cfg.dualTime == 1) ? var.c_d["ro"] : var.c_d["roN"],
                var.c_d["res_roY"+i], var.c_d["transport_diag_Y"+i], var.c_d["src_jac_Y"+i],
                var.c_d["dq_roY"+i+"_old"], var.c_d["dq_roY"+i],
                (cfg.discretization == "node") ? var.c_d["scalarDirichletPin"] : nullptr);
        }
        for (int s = 0; s < nSpecies; ++s) {
            const std::string i = std::to_string(s);
            std::swap(var.c_d["dq_roY"+i], var.c_d["dq_roY"+i+"_old"]);
            periodicBroadcastArray_d_wrapper(cfg, cuda_cfg, msh, var.c_d["dq_roY"+i+"_old"]);   // 周期 dq 整合 (§4.1-5)
        }
    }
    // 最終補正 δ(ρY_s)* は dq_roY{s}_old に残る。ポインタ配列を再構築。
    rebuildDqYoldPtrs(var, nSpecies);

    // --- 接空間射影 z_s = δz*_s - Y_s Σδz* (dq_old を in-place で z に) ---
    species_eos_project_tangent_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, nSpecies, g_dqYold_dev, g_roYN_dev, g_roPred_dev);

    // --- 解析 δp_Y をセルごとに評価 ---
    if (g_dpY_cap < msh.nCells_all) {
        if (g_dpY_eos) cudaFree(g_dpY_eos);
        gpuErrchk( cudaMalloc((void**)&g_dpY_eos, bytes) );
        g_dpY_cap = msh.nCells_all;
    }
    cudaMemset(g_dpY_eos, 0, bytes);
    species_eos_dp_cell_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, nSpecies, thermo_species_device_ptr(),
        g_dqYold_dev, g_roYN_dev, g_roPred_dev, var.c_d["T"], g_dpY_eos);

    // --- クロスエネルギー流束を res_roe へ移項 (内部 normal plane のみ) ---
    if (periodicNodeActive(cfg, msh)) {
        // node 周期 (codex plan-2 M3): res_roe は gather 済みなので直接足すと seam で部分 CV 分になる。クロス項は
        // 独立バッファに組み、その追加分だけを周期 gather (和→broadcast) してから res_roe に加える。
        static flow_float* s_cross = nullptr; static geom_int s_cap = 0;
        if (s_cap < msh.nCells_all) { if (s_cross) cudaFree(s_cross); gpuErrchk( cudaMalloc((void**)&s_cross, bytes) ); s_cap = msh.nCells_all; }
        cudaMemset(s_cross, 0, bytes);
        species_eos_cross_flux_d<<<cuda_cfg.dimGrid_plane, cuda_cfg.dimBlock>>>(
            msh.nNormalPlanes, msh.map_plane_cells_d,
            var.p_d["massflux"], g_dpY_eos, g_roPred_dev, s_cross);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
        periodicGatherArray_d_wrapper(cfg, cuda_cfg, msh, s_cross);
        passive_axpy1_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells, var.c_d["res_roe"], s_cross);
    } else {
        species_eos_cross_flux_d<<<cuda_cfg.dimGrid_plane, cuda_cfg.dimBlock>>>(
            msh.nNormalPlanes, msh.map_plane_cells_d,
            var.p_d["massflux"], g_dpY_eos, g_roPred_dev, var.c_d["res_roe"]);
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void speciesEOSFinalCommit_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var)) return;
    const int nSpecies = var.nSpeciesRegistered;
    rebuildDqYoldPtrs(var, nSpecies);   // PredictInject 後と同じ z (dq_old) を指す
    species_eos_final_commit_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, nSpecies, g_roY_dev, g_roYN_dev, g_dqYold_dev,
        var.c_d["ro"], g_roPred_dev ? g_roPred_dev : var.c_d["roN"]);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void speciesUpdateOuter_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!speciesEnabled(var)) return;

    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        gpuErrchk( cudaMemcpy(var.c_d["roY"+i+"N"], var.c_d["roY"+i], bytes, cudaMemcpyDeviceToDevice) );
        gpuErrchk( cudaMemcpy(var.c_d["roY"+i+"M"], var.c_d["roY"+i], bytes, cudaMemcpyDeviceToDevice) );
    }
}

void speciesUpdateInner_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg; (void)cuda_cfg;
    if (!speciesEnabled(var)) return;

    const size_t bytes = msh.nCells_all * sizeof(flow_float);
    for (int s = 0; s < var.nSpeciesRegistered; s++) {
        const std::string i = std::to_string(s);
        gpuErrchk( cudaMemcpy(var.c_d["roY"+i+"M"], var.c_d["roY"+i], bytes, cudaMemcpyDeviceToDevice) );
    }
}

// =============================================================================
// 受動種 (排気トレーサ・凝縮モーメント) 基盤 — 化学種カーネルの受動種ポインタ配列への適用
// (plans/active/species-passive-scalar-unification.md §4.1–4.3; passiveScalarScheme 1)。
// 化学種の anonymous-namespace カーネル (dirichlet/neumann/pin/gradient/dplur sweep) を共用するため本 TU に置く。
// =============================================================================
namespace {

int g_nPassive = 0, g_qTracer = -1, g_qMom0 = -1;
std::vector<std::string> g_pCons, g_pPrim;
// host 側ポインタ表 (q ごとの単一配列アクセス用) と device ポインタ配列 (カーネルの flow_float** 引数用)。
std::vector<flow_float*> h_p_rophi, h_p_rophiN, h_p_res, h_p_diag, h_p_sj, h_p_prim, h_p_gx, h_p_gy, h_p_gz, h_p_lim, h_p_corr, h_p_limc;
flow_float** g_p_rophi_dev = nullptr; flow_float** g_p_res_dev = nullptr; flow_float** g_p_diag_dev = nullptr;
flow_float** g_p_sj_dev = nullptr;    flow_float** g_p_prim_dev = nullptr;
flow_float** g_p_gx_dev = nullptr;    flow_float** g_p_gy_dev = nullptr;  flow_float** g_p_gz_dev = nullptr;
flow_float** g_p_lim_dev = nullptr;
flow_float*  g_Pface_dev = nullptr;   int g_Pface_nPlanes = 0;
flow_float*  g_p_roPre_dev = nullptr; bool g_p_roPre_valid = false;   // ρ_pre (流れ更新前の ρ; §5.1 #19)
double*      g_p_stats_dev = nullptr; // [nPassive*8]: 0 lo, 1 hi, 2 abs (floor; 全期間積算), 3 total (最新), 4 lim 量 Σ(1−θ)|δ|V (積算), 5 lim 作動セル数 (積算)
int*         g_p_thetaMin_dev = nullptr; // [nPassive]: θ_b の最小 (×1e9, atomicMin; ログ時にリセット)

std::vector<double> g_p_stats_last;   // 前回ログ時の積算 (per-step 平均用)
int          g_p_stats_last_step = 0;
flow_float*  g_zero_bvar = nullptr;   geom_int g_zero_cap = 0;

constexpr flow_float kNoFloor = static_cast<flow_float>(-1.0e30);
std::vector<double> g_fct_stats_total;   // FCT 補正の全期間積算 (log 用; 実体は下の FCT 節)
double g_fct_last_lin_res = 0.0; int g_fct_last_sweeps = 0;
std::vector<double> g_fct_budget_total, g_fct_max_relres, g_fct_max_rh, g_fct_closure_total;   // 全期間積算 / monitor 区間の最大 / 密度整合 (E_ρ rel 最大, 上限逸脱量, 個数)
std::vector<double> g_fct_run_max_relres, g_fct_run_max_rh, g_p_initial_total;   // 全期間の最大 / 計算開始前の総量 (収支の独立照合用; passiveRecordInitialTotals)
bool g_fct_nonfinite = false;   // 診断値 (残差ノルム等) に非有限が出たら true のまま (log に出す; ゲートは FAIL)   // 更新カーネルの floor を無効化 (floor は passive_bounds_d が担う)

flow_float** uploadPtrs(const std::vector<flow_float*>& h)
{
    flow_float** d = nullptr;
    const size_t b = h.size()*sizeof(flow_float*);
    gpuErrchk( cudaMalloc((void**)&d, b) );
    gpuErrchk( cudaMemcpy(d, h.data(), b, cudaMemcpyHostToDevice) );
    return d;
}

// ρφ = ρφ_N + δ (DPLUR 増分の commit; floor なし — 後段 passive_bounds_d)。
__global__ void passive_commit_d(geom_int nCells, flow_float* rophi, flow_float* rophiN, flow_float* dq)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < nCells) rophi[ic] = rophiN[ic] + dq[ic];
}

ScalarTransportDesc buildPassiveDesc(variables& var, int q)
{
    const std::string& c = g_pCons[q];
    const std::string& p = g_pPrim[q];
    return ScalarTransportDesc{
        var.c_d[p], nullptr, nullptr, nullptr,
        var.c_d[c], var.c_d[c+"N"], var.c_d[c+"M"],
        var.c_d["res_"+c], var.c_d["res_"+c+"_m"],
        var.c_d["src_jac_"+p], var.c_d["transport_diag_"+p],
        kNoFloor, static_cast<flow_float>(0.0),
        0
    };
}

}  // namespace

int  passive_count()          { return g_nPassive; }
int  passive_tracer_index()   { return g_qTracer; }
int  passive_moment_index0()  { return g_qMom0; }
const std::vector<std::string>& passive_cons_names() { return g_pCons; }
const std::vector<std::string>& passive_prim_names() { return g_pPrim; }
flow_float** passive_P_device_ptr()       { return g_p_prim_dev; }
flow_float** passive_dPdx_device_ptr()    { return g_p_gx_dev; }
flow_float** passive_dPdy_device_ptr()    { return g_p_gy_dev; }
flow_float** passive_dPdz_device_ptr()    { return g_p_gz_dev; }
flow_float** passive_limiter_device_ptr() { return g_p_lim_dev; }
flow_float*  passive_Pface_device_ptr()   { return g_Pface_dev; }

bool passiveSchemeEnabled(const solverConfig& cfg) { return cfg.passiveScalarScheme == 1 && g_nPassive > 0; }

flow_float scalarDtScale(const solverConfig& cfg)
{
    if (cfg.timeIntegration == 11 && cfg.scalarCflMax > 0.0 && cfg.cfl_pseudo > cfg.scalarCflMax) {
        return static_cast<flow_float>(cfg.scalarCflMax / cfg.cfl_pseudo);
    }
    return static_cast<flow_float>(1.0);
}

void passiveInit_d(solverConfig& cfg, variables& var)
{
    g_nPassive = 0; g_qTracer = -1; g_qMom0 = -1;
    g_pCons.clear(); g_pPrim.clear();
    if (var.tracerRegistered != 0) { g_qTracer = 0; g_pCons.push_back("roXi"); g_pPrim.push_back("Xi"); }
    if (!var.condMomentConsNames.empty()) {
        g_qMom0 = (int)g_pCons.size();
        for (const auto& c : var.condMomentConsNames) { g_pCons.push_back(c); g_pPrim.push_back(c.substr(2)); }
    }
    g_nPassive = (int)g_pCons.size();
    if (g_nPassive == 0) return;

    h_p_rophi.clear(); h_p_rophiN.clear(); h_p_res.clear(); h_p_diag.clear(); h_p_sj.clear(); h_p_prim.clear();
    h_p_gx.clear(); h_p_gy.clear(); h_p_gz.clear(); h_p_lim.clear(); h_p_corr.clear(); h_p_limc.clear();
    for (int q = 0; q < g_nPassive; ++q) {
        const std::string& c = g_pCons[q]; const std::string& p = g_pPrim[q];
        h_p_rophi.push_back(var.c_d.at(c));
        h_p_rophiN.push_back(var.c_d.at(c+"N"));
        h_p_res.push_back(var.c_d.at("res_"+c));
        h_p_diag.push_back(var.c_d.at("transport_diag_"+p));
        h_p_sj.push_back(var.c_d.at("src_jac_"+p));
        h_p_prim.push_back(var.c_d.at(p));
        h_p_gx.push_back(var.c_d.at("d"+p+"dx")); h_p_gy.push_back(var.c_d.at("d"+p+"dy")); h_p_gz.push_back(var.c_d.at("d"+p+"dz"));
        h_p_lim.push_back(var.c_d.at("limiter_"+p));
        h_p_corr.push_back(var.c_d.at("passiveFloorCorr_"+p));
        h_p_limc.push_back(var.c_d.at("passiveLimCorr_"+p));
    }
    g_p_rophi_dev = uploadPtrs(h_p_rophi);
    g_p_res_dev   = uploadPtrs(h_p_res);
    g_p_diag_dev  = uploadPtrs(h_p_diag);
    g_p_sj_dev    = uploadPtrs(h_p_sj);
    g_p_prim_dev  = uploadPtrs(h_p_prim);
    g_p_gx_dev = uploadPtrs(h_p_gx); g_p_gy_dev = uploadPtrs(h_p_gy); g_p_gz_dev = uploadPtrs(h_p_gz);
    g_p_lim_dev   = uploadPtrs(h_p_lim);
    gpuErrchk( cudaMalloc((void**)&g_p_stats_dev, (size_t)g_nPassive*8*sizeof(double)) );
    gpuErrchk( cudaMemset(g_p_stats_dev, 0, (size_t)g_nPassive*8*sizeof(double)) );
    g_p_stats_last.assign((size_t)g_nPassive*8, 0.0);
    gpuErrchk( cudaMalloc((void**)&g_p_thetaMin_dev, (size_t)g_nPassive*sizeof(int)) );
    { std::vector<int> one((size_t)g_nPassive, 1000000000); gpuErrchk( cudaMemcpy(g_p_thetaMin_dev, one.data(), one.size()*sizeof(int), cudaMemcpyHostToDevice) ); }
    g_p_stats_last_step = 0;
    if (cfg.passiveScalarScheme == 1 && cfg.speciesFaceReconstruction >= 2 && !(cfg.solver == "SLAU" || cfg.solver == "SLAU2"))
        std::cout << "passiveInit_d: WARNING speciesFaceReconstruction=" << cfg.speciesFaceReconstruction << " but solver=" << cfg.solver
                  << " is not SLAU/SLAU2: passive S3 face reconstruction and passiveFct are INACTIVE (passives use 1st-order upwind)\n";
    std::cout << "passiveInit_d: nPassive=" << g_nPassive << " (tracer q=" << g_qTracer << ", moments q0=" << g_qMom0
              << ") passiveScalarScheme=" << cfg.passiveScalarScheme
              << (cfg.passiveScalarScheme == 1 ? (std::string(" passiveImplicitCoupling=") + std::to_string(cfg.passiveImplicitCoupling)
                                                  + " passiveImplicitRelax=" + std::to_string(cfg.passiveImplicitRelax)) : std::string(" (legacy scalar path)"))
              << "\n";
}

flow_float* passive_Pface_alloc(int nPlanes)
{
    if (g_nPassive == 0) return nullptr;
    if (g_Pface_dev == nullptr || g_Pface_nPlanes < nPlanes) {
        if (g_Pface_dev) cudaFree(g_Pface_dev);
        gpuErrchk( cudaMalloc((void**)&g_Pface_dev, (size_t)nPlanes*g_nPassive*sizeof(flow_float)) );
        gpuErrchk( cudaMemset(g_Pface_dev, 0, (size_t)nPlanes*g_nPassive*sizeof(flow_float)) );
        g_Pface_nPlanes = nPlanes;
    }
    return g_Pface_dev;
}

void passivePrimitive_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int q0, int nq)
{
    (void)cfg;
    for (int q = q0; q < q0+nq; ++q) {
        species_primitive_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, var.c_d["ro"], h_p_rophi[q], h_p_prim[q]);
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void passiveGradient_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!passiveSchemeEnabled(cfg) || cfg.speciesFaceReconstruction < 1) return;
    const size_t bytes = (size_t)msh.nCells_all*sizeof(flow_float);
    for (int q = 0; q < g_nPassive; ++q) { cudaMemset(h_p_gx[q], 0, bytes); cudaMemset(h_p_gy[q], 0, bytes); cudaMemset(h_p_gz[q], 0, bytes); }
    // mesh.scalarGradient: lsq (node のみ)。化学種と同じ扱い (normalize なし、周期合併は periodicSeamMergeActive のときここで)。
    if (scalarGradientLsqActive(cfg)) {
        lsqScalarGradient_d_wrapper(cuda_cfg, msh, g_nPassive, g_p_prim_dev, g_p_gx_dev, g_p_gy_dev, g_p_gz_dev);
        if (preGatherDumpEnabled()) {   // 診断 (FORGE_DUMP_PREGATHER、既定 off・出力専用): 周期合併の直前の局所配列
            std::vector<std::array<const flow_float*, 3>> gp;
            for (int q = 0; q < g_nPassive; ++q) gp.push_back({h_p_gx[q], h_p_gy[q], h_p_gz[q]});
            preGatherDump("passive_lsq.loop1", msh.nCells, g_pPrim, gp);
        }
        if (periodicSeamMergeActive(cfg, msh)) {
            for (int q = 0; q < g_nPassive; ++q) {
                periodicGatherArray_d_wrapper(cfg, cuda_cfg, msh, h_p_gx[q]);
                periodicGatherArray_d_wrapper(cfg, cuda_cfg, msh, h_p_gy[q]);
                periodicGatherArray_d_wrapper(cfg, cuda_cfg, msh, h_p_gz[q]);
            }
        }
        return;
    }
    flow_float* gvol = (cfg.isAxisymmetric == 1) ? var.c_d["A_planar"] : var.c_d["volume"];
    flow_float* gsx = (cfg.isAxisymmetric == 1) ? var.p_d["sx_planar"] : var.p_d["sx"];
    flow_float* gsy = (cfg.isAxisymmetric == 1) ? var.p_d["sy_planar"] : var.p_d["sy"];
    flow_float* gsz = (cfg.isAxisymmetric == 1) ? var.p_d["sz_planar"] : var.p_d["sz"];
    static bool s_pgDumped = false;
    flow_float* pgDump = scalarGradDumpBegin("passive", s_pgDumped, g_nPassive, msh.nPlanes);
    species_gradient_d<<<cuda_cfg.dimGrid_plane, cuda_cfg.dimBlock>>>(
        msh.nCells, msh.nPlanes, msh.map_plane_cells_d, gvol, var.p_d["fx"], gsx, gsy, gsz,
        g_nPassive, g_p_prim_dev, g_p_gx_dev, g_p_gy_dev, g_p_gz_dev,
        periodicSeamMergeActive(cfg, msh) ? 1 : 0, msh.planePeriodic_d, pgDump);
    scalarGradDumpEnd("passive", pgDump, g_nPassive, msh.nPlanes);
    species_gradient_normalize_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, gvol, g_nPassive, g_p_gx_dev, g_p_gy_dev, g_p_gz_dev);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void passiveBoundary_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, bcond& bc, mesh& msh, variables& var,
                               int q, flow_float* bvar_dirichlet)
{
    (void)msh;
    if (bc.iPlanes.empty()) return;
    const geom_int nb = static_cast<geom_int>(bc.iPlanes.size());
    const bool isInlet = bc.bcondKind.rfind("inlet_", 0) == 0;
    if (isInlet) {
        flow_float* bv = bvar_dirichlet;
        if (bv == nullptr) {
            // 凝縮モーメントの dry 入口 (ρφ=0): ゼロの per-face 配列を共用する。
            if (g_zero_cap < nb) {
                if (g_zero_bvar) cudaFree(g_zero_bvar);
                gpuErrchk( cudaMalloc((void**)&g_zero_bvar, (size_t)nb*sizeof(flow_float)) );
                gpuErrchk( cudaMemset(g_zero_bvar, 0, (size_t)nb*sizeof(flow_float)) );
                g_zero_cap = nb;
            }
            bv = g_zero_bvar;
        }
        species_dirichlet_boundary_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(
            nb, bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d, var.c_d["ro"], bv,
            h_p_rophi[q], h_p_prim[q], var.c_d["scalarDirichletPin"], (cfg.discretization == "node") ? 1 : 0);
    } else {
        species_neumann_boundary_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(
            nb, bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d, h_p_rophi[q], h_p_prim[q]);
    }
}

void passiveAdvection_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int q0, int nq)
{
    const size_t bytes = (size_t)msh.nCells*sizeof(flow_float);
    for (int q = q0; q < q0+nq; ++q) {
        CHECK_CUDA_ERROR(cudaMemset(h_p_res[q],  0, bytes));
        CHECK_CUDA_ERROR(cudaMemset(h_p_diag[q], 0, bytes));
        CHECK_CUDA_ERROR(cudaMemset(h_p_sj[q],   0, bytes));
    }
    // 診断 G3-a: 凝縮モーメント (液・Q2・Q1・Q0 の 4 本を一括で呼ぶとき) だけ記録する (トレーサの呼び出しは記録しない)
    const bool tpoMom = (g_qMom0 >= 0 && q0 == g_qMom0 && nq == 4);
    if (tpoMom) tpoSnap("cm_zero");
    const bool s3 = (cfg.speciesFaceReconstruction >= 2 && g_Pface_dev != nullptr
                     && (cfg.solver == "SLAU" || cfg.solver == "SLAU2"));
    if (s3) {
        // S3: convectiveFlux (SLAU) が書いた受動種の upwind 面値 Pface[ip*nPassive+q] で移流 (対角は 1 次風上)。
        dim3 dimGrid_nh = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));
        species_advection_faceY_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(
            msh.nCells, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
            var.c_d["ro"], var.p_d["massflux"], nq, g_Pface_dev + q0, g_p_res_dev + q0, g_p_diag_dev + q0,
            (cfg.discretization == "node") ? 1 : 0, g_p_rophi_dev + q0, g_nPassive, nullptr,
            tpoMom ? tpoFace(TPO_K_ADV_PA, TPO_KIND_ADV, 1) : TpoFacePtr{});   // 診断 G3-a (既定は無効; 成分 1 + s = g, Q2, Q1, Q0)
    } else {
        // 1 次風上 (化学種の既定経路と同じ融合カーネル)。
        std::vector<ScalarTransportDesc> descs;
        for (int q = q0; q < q0+nq; ++q) descs.push_back(buildPassiveDesc(var, q));
        scalarTransportResidualMulti_d(cfg, cuda_cfg, msh, var, descs.data(), (int)descs.size());
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void passiveDiffusion_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int q)
{
    if (cfg.viscMethod == 0) return;
    dim3 dimGrid_nh = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));
    passive_diffusion_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(
        msh.nCells, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
        var.p_d["ge_x"], var.p_d["ge_y"], var.p_d["ge_z"],   // 面ごとの差 e (段 ③)
        var.p_d["fx"], var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
        h_p_prim[q], h_p_res[q], h_p_diag[q],
        var.c_d["ro"], var.c_d["vis_lam"], var.c_d["vis_turb"],
        static_cast<flow_float>(cfg.Sc), static_cast<flow_float>(cfg.Sc_t), (cfg.discretization == "node") ? 1 : 0);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void passivePinResidual_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!passiveSchemeEnabled(cfg) || cfg.discretization != "node") return;
    species_pin_residual_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, g_nPassive, g_p_res_dev, g_p_sj_dev, var.c_d["scalarDirichletPin"]);
    gpuErrchk( cudaPeekAtLastError() );
}

bool passiveDPLURIncrement_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int q0, int nq)
{
    if (cfg.timeIntegration != 11 || cfg.passiveImplicitCoupling != 1) return false;
    const size_t bytes = (size_t)msh.nCells_all*sizeof(flow_float);
    for (int q = q0; q < q0+nq; ++q) {
        cudaMemset(var.c_d["dq_"+g_pCons[q]], 0, bytes);
        cudaMemset(var.c_d["dq_"+g_pCons[q]+"_old"], 0, bytes);
    }
    const int nSweep = std::max(1, cfg.nStepInner);
    const flow_float relax = static_cast<flow_float>(cfg.passiveImplicitRelax);
    const flow_float dts   = scalarDtScale(cfg);
    flow_float* roRef = (cfg.unsteady == 1 && cfg.dualTime == 1) ? var.c_d["ro"] : var.c_d["roN"];   // dual-time は現在の ρ
    for (int iSweep = 0; iSweep < nSweep; ++iSweep) {
        for (int q = q0; q < q0+nq; ++q) {
            const std::string& c = g_pCons[q];
            dplurSweepOnce(cfg, cuda_cfg, msh, var, relax, dts, roRef,
                h_p_res[q], h_p_diag[q], h_p_sj[q],
                var.c_d["dq_"+c+"_old"], var.c_d["dq_"+c],
                (cfg.discretization == "node") ? var.c_d["scalarDirichletPin"] : nullptr);
        }
        for (int q = q0; q < q0+nq; ++q) {
            const std::string& c = g_pCons[q];
            std::swap(var.c_d["dq_"+c], var.c_d["dq_"+c+"_old"]);
            periodicBroadcastArray_d_wrapper(cfg, cuda_cfg, msh, var.c_d["dq_"+c+"_old"]);   // 周期 dq 整合
        }
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
    return true;
}

const geom_int* passive_periodic_root(solverConfig& cfg, mesh& msh)
{
    return periodicNodeActive(cfg, msh) ? msh.periodicRoot_d : nullptr;
}
flow_float* passive_limCorr_cell_ptr(int q) { return (q >= 0 && q < g_nPassive) ? h_p_limc[q] : nullptr; }
double*     passive_lim_stats_ptr(int q)    { return (q >= 0 && q < g_nPassive) ? g_p_stats_dev + (size_t)q*8 + 4 : nullptr; }
// 診断 G3-b: 受動種 q が凝縮種 0 のモーメント (g, Q2, Q1, Q0) なら成分番号 1..4、それ以外は −1 (記録しない)。
static int tpuMomentComp(int q) { return (g_qMom0 >= 0 && q >= g_qMom0 && q < g_qMom0 + 4) ? 1 + (q - g_qMom0) : -1; }

void passiveSaveRhoPre_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cuda_cfg;
    if (!passiveSchemeEnabled(cfg) || cfg.timeIntegration != 11) { g_p_roPre_valid = false; return; }
    const size_t bytes = (size_t)msh.nCells_all*sizeof(flow_float);
    if (g_p_roPre_dev == nullptr) gpuErrchk( cudaMalloc((void**)&g_p_roPre_dev, bytes) );
    gpuErrchk( cudaMemcpy(g_p_roPre_dev, var.c_d["ro"], bytes, cudaMemcpyDeviceToDevice) );
    g_p_roPre_valid = true;
}

void passiveAddRhoTerm_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int q0, int nq)
{
    if (!passiveSchemeEnabled(cfg) || cfg.timeIntegration != 11 || !g_p_roPre_valid) return;
    for (int q = q0; q < q0+nq; ++q) {
        passive_add_rho_term_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells, h_p_rophi[q], h_p_rophiN[q], g_p_roPre_dev, var.c_d["ro"],
            tpuSlots(), tpuMomentComp(q));   // 診断 G3-b (既定 nullptr)
    }
    gpuErrchk( cudaPeekAtLastError() );
}

bool passiveRecordStage(const solverConfig& cfg, int loop)
{
    if (cfg.timeIntegration == 11) return true;
    return loop == cfg.perStepIterationCount() - 1;
}

void passiveLimitIncrement_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int q0, int nq, bool record)
{
    // FCT 有効 (dual-time S3) では sub-iter 内の増分スケーリング θ_b は使わない: 非保存で、step 末尾の FCT が落とした分を戻せない
    // (Venkat のステップで −0.7 %; 2026-09-17 run_0134/0135)。有界化は FCT 補正 + 最後の floor が担う。
    if (passiveFctActive(cfg)) return;
    const geom_int* root = passive_periodic_root(cfg, msh);
    static double* s_scratch = nullptr; static int* s_scratchI = nullptr;
    if (!s_scratch) { gpuErrchk( cudaMalloc((void**)&s_scratch, 8*sizeof(double)) ); gpuErrchk( cudaMalloc((void**)&s_scratchI, sizeof(int)) ); }
    for (int q = q0; q < q0+nq; ++q) {
        passive_limit_increment_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
            msh.nCells, h_p_rophi[q], h_p_rophiN[q], (q == g_qTracer) ? 1 : 0, var.c_d["ro"], var.c_d["volume"],
            record ? h_p_limc[q] : nullptr, record ? g_p_stats_dev + (size_t)q*8 + 4 : s_scratch, record ? g_p_thetaMin_dev + q : s_scratchI, root,
            (cfg.timeIntegration == 11 && g_p_roPre_valid) ? g_p_roPre_dev : nullptr,
            tpuSlots(), tpuMomentComp(q));   // 診断 G3-b (既定 nullptr)
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void passiveBounds_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int q0, int nq, bool record, double* tcAcc)
{
    (void)cfg;
    const geom_int* root = passive_periodic_root(cfg, msh);
    static double* s_scratch = nullptr;
    if (!s_scratch) gpuErrchk( cudaMalloc((void**)&s_scratch, 8*sizeof(double)) );
    for (int q = q0; q < q0+nq; ++q) {
        // FCT 有効時、sub-iter 内 (record=true でも物理 step 末尾以外) のトレーサ floor は掛けない (熱力学に入らないので中間の逸脱は無害;
        // 掛けると FCT 前に非保存で切ってしまう)。モーメントは g<0 が EOS を壊すので floor を残す (収支に記録)。物理 step 末尾の floor は
        // main が passiveBounds_d_wrapper(..., record=true) で呼ぶ (cfg.dualTimeSubIter = -1 で区別)。
        if (passiveFctActive(cfg) && q == g_qTracer && cfg.dualTimeSubIter >= 0) continue;
        if (record) gpuErrchk( cudaMemset(g_p_stats_dev + (size_t)q*8 + 3, 0, sizeof(double)) );   // 総量は最新値
        passive_bounds_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
            msh.nCells, h_p_rophi[q], (q == g_qTracer) ? 1 : 0, var.c_d["ro"], var.c_d["volume"],
            record ? h_p_corr[q] : nullptr, record ? g_p_stats_dev + (size_t)q*8 : s_scratch, root,
            tcAcc, (tcAcc != nullptr && g_qMom0 >= 0 && q >= g_qMom0 && q < g_qMom0 + 4) ? 2 + (q - g_qMom0) : -1,
            tpuSlots(), tpuMomentComp(q));   // 診断 G3-b (既定 nullptr)
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void passiveCommitIncrement_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int q)
{
    (void)cfg;
    passive_commit_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, h_p_rophi[q], h_p_rophiN[q], var.c_d["dq_"+g_pCons[q]+"_old"]);
    gpuErrchk( cudaPeekAtLastError() );
}

void passiveTracerUpdate_d_wrapper(int loop, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    const int q = g_qTracer;
    if (q < 0) return;
    if (passiveDPLURIncrement_d_wrapper(cfg, cuda_cfg, msh, var, q, 1)) {
        passiveCommitIncrement_d_wrapper(cfg, cuda_cfg, msh, var, q);
    } else {
        const ScalarTransportDesc desc = buildPassiveDesc(var, q);
        const flow_float relax = (cfg.timeIntegration == 11) ? static_cast<flow_float>(cfg.passiveImplicitRelax) : static_cast<flow_float>(1.0);
        scalarTimeIntegration_d(loop, cfg, cuda_cfg, msh, var, desc, relax, scalarDtScale(cfg));
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
    const bool rec = passiveRecordStage(cfg, loop);
    passiveAddRhoTerm_d_wrapper(cfg, cuda_cfg, msh, var, q, 1);           // + φ_N δρ (流れの密度更新と整合; #19)
    passiveLimitIncrement_d_wrapper(cfg, cuda_cfg, msh, var, q, 1, rec);   // θ_b で輸送増分 z を縮め [0,ρ] を保つ (M5)
    passiveBounds_d_wrapper(cfg, cuda_cfg, msh, var, q, 1, rec);           // 最後の砦 (通常は無作用)
}

// 周期 root のみの総量 Σ ρφ V (double) を受動種ごとに (現在の確定状態)。
std::vector<double> passiveTotalsNow(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    std::vector<double> h((size_t)g_nPassive, 0.0);
    if (g_nPassive == 0) return h;
    static double* d = nullptr; if (!d) gpuErrchk( cudaMalloc((void**)&d, (size_t)g_nPassive*sizeof(double)) );
    gpuErrchk( cudaMemset(d, 0, (size_t)g_nPassive*sizeof(double)) );
    const geom_int* root = passive_periodic_root(cfg, msh);
    for (int q = 0; q < g_nPassive; ++q)
        passive_total_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(msh.nCells, h_p_rophi[q], var.c_d["volume"], root, d + q);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    gpuErrchk( cudaMemcpy(h.data(), d, h.size()*sizeof(double), cudaMemcpyDeviceToHost) );
    return h;
}
void passiveRecordInitialTotals_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!passiveSchemeEnabled(cfg)) return;
    g_p_initial_total = passiveTotalsNow(cfg, cuda_cfg, msh, var);
    for (int q = 0; q < g_nPassive; ++q) printf("[passive] initial total %-8s %.12e (root-only, before the first physical step)\n", g_pCons[q].c_str(), g_p_initial_total[(size_t)q]);
}

std::vector<double> passiveFloorCorrTotals()
{
    std::vector<double> h((size_t)g_nPassive*8, 0.0);
    if (g_nPassive > 0) gpuErrchk( cudaMemcpy(h.data(), g_p_stats_dev, h.size()*sizeof(double), cudaMemcpyDeviceToHost) );
    return h;
}

void passiveFloorCorrLog_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int iStep)
{
    if (!passiveSchemeEnabled(cfg)) return;
    const std::vector<double> h = passiveFloorCorrTotals();
    const std::vector<double> totNow = passiveTotalsNow(cfg, cuda_cfg, msh, var);   // 全後処理後の確定状態 (root のみ)
    const int nstep = std::max(1, (iStep + 1) - g_p_stats_last_step);
    std::vector<int> thmin((size_t)g_nPassive, 1000000000);
    gpuErrchk( cudaMemcpy(thmin.data(), g_p_thetaMin_dev, thmin.size()*sizeof(int), cudaMemcpyDeviceToHost) );
    for (int q = 0; q < g_nPassive; ++q) {
        const double* c = h.data() + (size_t)q*8;
        const double* l = g_p_stats_last.data() + (size_t)q*8;
        const double tot = totNow[(size_t)q];
        if (g_p_initial_total.size() != (size_t)g_nPassive) g_p_initial_total.assign((size_t)g_nPassive, -1.0);
        const double rel = (tot > 0.0) ? c[2]/tot : (c[2] > 0.0 ? 1.0 : 0.0);    // 総量 0 で補正が非ゼロなら 1 (= FAIL; plan-7 M1)
        const double rell = (tot > 0.0) ? c[4]/tot : (c[4] > 0.0 ? 1.0 : 0.0);
        // floorCorr = 硬い floor (最後の砦) の収支、limCorr = 増分スケーリング θ_b の縮小量 (M5)。総量・収支は周期 root のみ (M2)。
        printf("[passive] step %d floorCorr %-8s per-step(avg %d): lo %.3e hi %.3e abs %.3e | cumulative: lo %.9e hi %.9e abs %.9e | total %.12e rel(abs/total) %.6e"
               " | limCorr per-step %.3e cumulative abs %.9e signed %.9e rel %.6e cells %.0f thetaMin(interval) %.4f initialTotal %.12e\n",
               iStep + 1, g_pCons[q].c_str(), nstep,
               (c[0]-l[0])/nstep, (c[1]-l[1])/nstep, (c[2]-l[2])/nstep, c[0], c[1], c[2], tot, rel,
               (c[4]-l[4])/nstep, c[4], c[6], rell, c[5]-l[5], thmin[(size_t)q]*1.0e-9, g_p_initial_total[(size_t)q]);
        if (passiveFctActive(cfg) && (size_t)g_nPassive*8 == g_fct_stats_total.size()) {
            const double* f = g_fct_stats_total.data() + (size_t)q*8;
            const double* bg = g_fct_budget_total.data() + (size_t)q*4;
            printf("[passive]   fctCorr %-8s cumulative: dropped antidiffusion %.9e (rel %.6e) faces %.0f prelimited %.9e pinCorr %.9e (rel %.6e) baseViol %.9e (rel %.6e)"
                   " bndFluxSigned %.12e bndDropped %.9e upperViol %.9e (rel %.6e) | budget: srcHist %.12e remSigned %.12e remAbs %.9e (rel %.6e) increment %.12e"
                   " | qL rel-residual interval-max %.2e run-max %.2e (sweeps last %d) HO residual rel interval-max %.2e run-max %.2e nonfinite %d\n",
                   g_pCons[q].c_str(), f[0], (tot > 0.0) ? f[0]/tot : (f[0] > 0.0 ? 1.0 : 0.0), f[1], f[2], f[3], (tot > 0.0) ? f[3]/tot : (f[3] > 0.0 ? 1.0 : 0.0),
                   f[4], (tot > 0.0) ? f[4]/tot : (f[4] > 0.0 ? 1.0 : 0.0), f[5], f[6], f[7], (tot > 0.0) ? f[7]/tot : (f[7] > 0.0 ? 1.0 : 0.0),
                   bg[0], bg[1], bg[2], (tot > 0.0) ? bg[2]/tot : (bg[2] > 0.0 ? 1.0 : 0.0), bg[3],
                   g_fct_max_relres[(size_t)q], g_fct_run_max_relres[(size_t)q], g_fct_last_sweeps, g_fct_max_rh[(size_t)q], g_fct_run_max_rh[(size_t)q], g_fct_nonfinite ? 1 : 0);
            g_fct_max_relres[(size_t)q] = 0.0; g_fct_max_rh[(size_t)q] = 0.0;
        }
    }
    if (passiveFctActive(cfg) && g_fct_closure_total.size() == 4) {
        printf("[passive]   fctDensity: max rel E_rho (BE-form continuity remainder) %.2e | tracer upper-bound margin violation cumulative %.6e (cells %.0f)\n",
               g_fct_closure_total[0], g_fct_closure_total[1], g_fct_closure_total[2]);
        g_fct_closure_total[0] = 0.0;
    }
    if (g_qMom0 >= 0) {
        int ndeg = 0; const int nproj = condRealizViolReadReset(&ndeg);
        printf("[passive] step %d moment realizability corrections since last log: nearest-point %d, degenerate->monodisperse %d\n", iStep + 1, nproj, ndeg);
        const int nsp = (g_nPassive - g_qMom0)/4;
        const std::vector<double> cb = condClampBudgetTotals(nsp);
        for (int sp = 0; sp < nsp && (size_t)(sp+1)*8 <= cb.size(); ++sp) {
            const double* bq = cb.data() + (size_t)sp*8;
            const double* totg = h.data() + (size_t)(g_qMom0 + 4*sp)*8;   // 総量 [g,Q2,Q1,Q0] は受動種の stats[3]
            const double tg = totg[0*8+3], tQ2 = totg[1*8+3], tQ1 = totg[2*8+3], tQ0 = totg[3*8+3];
            auto relz = [](double ab, double t) { return (t > 0.0) ? ab/t : (ab > 0.0 ? 1.0 : 0.0); };   // 総量 0 で補正が非ゼロなら 1 (= FAIL)
            printf("[passive]   clampBudget species %d cumulative (signed/abs, rel to total): g %.6e/%.6e (%.6e) Q0 %.6e/%.6e (%.6e) Q1 %.6e/%.6e (%.6e) Q2 %.6e/%.6e (%.6e)\n", sp,
                   bq[0], bq[1], relz(bq[1], tg), bq[2], bq[3], relz(bq[3], tQ0), bq[4], bq[5], relz(bq[5], tQ1), bq[6], bq[7], relz(bq[7], tQ2));
        }
    }
    { std::vector<int> one((size_t)g_nPassive, 1000000000); gpuErrchk( cudaMemcpy(g_p_thetaMin_dev, one.data(), one.size()*sizeof(int), cudaMemcpyHostToDevice) ); }
    g_p_stats_last = h;
    g_p_stats_last_step = iStep + 1;
    fflush(stdout);
}

void passiveMirrorPeriodic_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)var;
    if (!passiveSchemeEnabled(cfg)) return;
    for (int q = 0; q < g_nPassive; ++q) periodicBroadcastArray_d_wrapper(cfg, cuda_cfg, msh, h_p_rophi[q]);
}

// =============================================================================
// dual-time (plan species-passive-scalar-unification §4.4; chem e296f0d0 の移植 + 受動種)
//   物理時間レベル Q^n (P), Q^{n-1} (PP) と BDF 項。履歴の有効数・係数は呼び出し側 (cfg.nHistoryValid, main.cpp) が
//   流れ・化学種・受動種で共有して決める (移植元のプロセス内カウンタは使わない)。
// =============================================================================
namespace {
__global__ void dt_shift_levels_d(geom_int nCells_all, const flow_float* q, flow_float* P, flow_float* PP)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells_all) { PP[ic] = P[ic]; P[ic] = q[ic]; }
}
__global__ void dt_init_levels_d(geom_int nCells_all, const flow_float* q, flow_float* P, flow_float* PP)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells_all) { P[ic] = q[ic]; PP[ic] = q[ic]; }
}
// res −= (V/Δt)(a q − b P + c PP), tdiag += V a/Δt (流れの addUnsteadyTimeTerm_d と同じ符号規約)。
__global__ void dt_add_unsteady_d(geom_int nCells, const geom_float* vol, flow_float dt, flow_float a, flow_float b, flow_float c,
                                  const flow_float* q, const flow_float* P, const flow_float* PP, flow_float* res, flow_float* tdiag)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic < nCells) {
        const flow_float V = static_cast<flow_float>(vol[ic]);
        res[ic]   -= (V / dt) * (a * q[ic] - b * P[ic] + c * PP[ic]);
        tdiag[ic] += V * a / dt;
    }
}
}  // namespace

void speciesInitDualTimeLevels_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var)) return;
    for (int s = 0; s < var.nSpeciesRegistered; ++s) {
        const std::string i = std::to_string(s);
        dt_init_levels_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, var.c_d["roY"+i], var.c_d["roY"+i+"P"], var.c_d["roY"+i+"PP"]);
    }
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
}

void speciesShiftDualTimeLevels_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg;
    if (!speciesEnabled(var)) return;
    for (int s = 0; s < var.nSpeciesRegistered; ++s) {
        const std::string i = std::to_string(s);
        dt_shift_levels_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, var.c_d["roY"+i], var.c_d["roY"+i+"P"], var.c_d["roY"+i+"PP"]);
    }
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
}

void speciesAddUnsteadyTimeTerm_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                                          flow_float a, flow_float b, flow_float c)
{
    if (!speciesEnabled(var)) return;
    for (int s = 0; s < var.nSpeciesRegistered; ++s) {
        const std::string i = std::to_string(s);
        dt_add_unsteady_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells, var.c_d["volume"], cfg.dt, a, b, c,
            var.c_d["roY"+i], var.c_d["roY"+i+"P"], var.c_d["roY"+i+"PP"], var.c_d["res_roY"+i], var.c_d["transport_diag_Y"+i]);
    }
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
}

void passiveInitDualTimeLevels_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    (void)cfg;
    for (int q = 0; q < g_nPassive; ++q) {
        dt_init_levels_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, h_p_rophi[q], var.c_d[g_pCons[q]+"P"], var.c_d[g_pCons[q]+"PP"]);
    }
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
}

void passiveShiftDualTimeLevels_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!passiveSchemeEnabled(cfg)) return;
    for (int q = 0; q < g_nPassive; ++q) {
        dt_shift_levels_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells_all, h_p_rophi[q], var.c_d[g_pCons[q]+"P"], var.c_d[g_pCons[q]+"PP"]);
    }
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
}

void passiveAddUnsteadyTimeTerm_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var,
                                          flow_float a, flow_float b, flow_float c)
{
    if (!passiveSchemeEnabled(cfg)) return;   // 旧経路 0 は物理時間項なし (従来どおり)
    for (int q = 0; q < g_nPassive; ++q) {
        dt_add_unsteady_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(msh.nCells, var.c_d["volume"], cfg.dt, a, b, c,
            h_p_rophi[q], var.c_d[g_pCons[q]+"P"], var.c_d[g_pCons[q]+"PP"], h_p_res[q], h_p_diag[q]);
    }
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
}

// =============================================================================
// dual-time の物理 step 末尾の保存的 FCT 補正 (plan §4.7 v4; passiveFct_d.cuh の kernel 群)。
//   呼び出し側 (main) が終了状態で assembleResidual を再評価してから passiveFctCorrect を呼び、後処理 (ピン・floor・実現可能性・EOS) の後に
//   passiveFctFinishHistory で確定状態から G^{n+1}, H^{n+1} を作る。BDF2 は流束形の BE (F^eff = (1/a)F + (c/a)G^n) として扱う。
// =============================================================================
namespace {
struct FctScratch {
    int nq = 0; geom_int nCells = 0; geom_int nPlanes = 0;
    std::vector<flow_float*> h_qL, h_rhs, h_nb, h_Pp, h_Pm, h_mx, h_mn, h_Rp, h_Rm, h_corr, h_base, h_qB, h_H, h_Hsrc, h_qP, h_qPP;
    flow_float **d_qL = nullptr, **d_rhs = nullptr, **d_nb = nullptr, **d_Pp = nullptr, **d_Pm = nullptr, **d_mx = nullptr, **d_mn = nullptr,
               **d_Rp = nullptr, **d_Rm = nullptr, **d_corr = nullptr, **d_base = nullptr, **d_qB = nullptr, **d_H = nullptr, **d_Hsrc = nullptr, **d_qP = nullptr, **d_qPP = nullptr;
    flow_float** d_roPtr = nullptr; flow_float** d_roNPtr = nullptr;   // 上限条件 Lρ − f の診断用 ({ro}, {roN})
    flow_float* divMeff = nullptr; flow_float* nbRho = nullptr; flow_float** d_nbRho = nullptr; flow_float* nbBnd = nullptr; flow_float** d_nbBnd = nullptr;
    flow_float** d_qLtr = nullptr; flow_float** d_qPtr = nullptr;   // トレーサの qL / q^n (単一エントリの pointer 配列)
    flow_float *diagAdv = nullptr, *diagDiff = nullptr, *cdiff = nullptr, *Aface = nullptr, *Apre = nullptr, *G = nullptr, *Gnew = nullptr, *mEff = nullptr;
    double* stats = nullptr;    // [nq*8]: 0 Σ|A^raw−αA^pre| (·Δt), 1 作動面数, 2 前制限量 (·Δt), 3 ピン行交換量, 4 基点の物理限界逸脱量 (·V), 5 境界の実現流束 (符号付き·Δt), 6 境界で落とした量 (·Δt), 7 予備
    double* acc = nullptr;      // [4*nq + 4]: [q] Jacobi 残差², [nq+q] rhs², [2nq+2q] r_H², [2nq+2q+1] (V a/Δt q)², 末尾 4: E_ρ², (V/Δt ρ)², 上限逸脱量, 個数
    double* budget = nullptr;   // [nq*4]: Σ H_src Δt, Σ H_rem Δt (符号付き), Σ|H_rem| Δt, Σ 総増分 (root のみ; 全期間積算は host)
    bool histValid = false;     // G/H が前 step (または checkpoint) から有効か
};
FctScratch g_fct;
double g_fct_last_rh_rel = 0.0;   // 最後の物理 step の max_q ||r_H||/||M q_H|| (sub-iter の反復誤差)

flow_float* allocCells(size_t n) { flow_float* p = nullptr; gpuErrchk( cudaMalloc((void**)&p, n*sizeof(flow_float)) ); gpuErrchk( cudaMemset(p, 0, n*sizeof(flow_float)) ); return p; }

void fctAlloc(mesh& msh, variables& var)
{
    if (g_fct.nq == g_nPassive && g_fct.nCells == msh.nCells_all) return;
    const size_t n = (size_t)msh.nCells_all;
    auto mk = [&](std::vector<flow_float*>& h) { h.clear(); for (int q = 0; q < g_nPassive; ++q) h.push_back(allocCells(n)); return uploadPtrs(h); };
    g_fct.d_qL = mk(g_fct.h_qL); g_fct.d_rhs = mk(g_fct.h_rhs); g_fct.d_nb = mk(g_fct.h_nb); g_fct.d_Pp = mk(g_fct.h_Pp); g_fct.d_Pm = mk(g_fct.h_Pm);
    g_fct.d_mx = mk(g_fct.h_mx); g_fct.d_mn = mk(g_fct.h_mn); g_fct.d_Rp = mk(g_fct.h_Rp); g_fct.d_Rm = mk(g_fct.h_Rm); g_fct.d_corr = mk(g_fct.h_corr);
    g_fct.d_base = mk(g_fct.h_base); g_fct.d_qB = mk(g_fct.h_qB); g_fct.d_H = mk(g_fct.h_H); g_fct.d_Hsrc = mk(g_fct.h_Hsrc);
    { std::vector<flow_float*> r{var.c_d.at("ro")}, rn{var.c_d.at("roN")}; g_fct.d_roPtr = uploadPtrs(r); g_fct.d_roNPtr = uploadPtrs(rn); }
    g_fct.divMeff = allocCells(n); g_fct.nbRho = allocCells(n); { std::vector<flow_float*> nb{g_fct.nbRho}; g_fct.d_nbRho = uploadPtrs(nb); }
    g_fct.nbBnd = allocCells(n); { std::vector<flow_float*> nb{g_fct.nbBnd}; g_fct.d_nbBnd = uploadPtrs(nb); }
    if (g_qTracer >= 0) { std::vector<flow_float*> a{g_fct.h_qL[(size_t)g_qTracer]}, b{var.c_d.at(g_pCons[(size_t)g_qTracer]+"P")}; g_fct.d_qLtr = uploadPtrs(a); g_fct.d_qPtr = uploadPtrs(b); }
    g_fct.h_qP.clear(); g_fct.h_qPP.clear();
    for (int q = 0; q < g_nPassive; ++q) { g_fct.h_qP.push_back(var.c_d.at(g_pCons[q]+"P")); g_fct.h_qPP.push_back(var.c_d.at(g_pCons[q]+"PP")); }
    g_fct.d_qP = uploadPtrs(g_fct.h_qP); g_fct.d_qPP = uploadPtrs(g_fct.h_qPP);
    g_fct.diagAdv = allocCells(n); g_fct.diagDiff = allocCells(n);
    const size_t fb = (size_t)msh.nPlanes*g_nPassive*sizeof(flow_float);
    gpuErrchk( cudaMalloc((void**)&g_fct.cdiff, (size_t)msh.nPlanes*sizeof(flow_float)) );
    gpuErrchk( cudaMalloc((void**)&g_fct.Aface, fb) ); gpuErrchk( cudaMemset(g_fct.Aface, 0, fb) );
    gpuErrchk( cudaMalloc((void**)&g_fct.Apre,  fb) ); gpuErrchk( cudaMemset(g_fct.Apre,  0, fb) );
    gpuErrchk( cudaMalloc((void**)&g_fct.G,     fb) ); gpuErrchk( cudaMemset(g_fct.G,     0, fb) );
    gpuErrchk( cudaMalloc((void**)&g_fct.Gnew,  fb) ); gpuErrchk( cudaMemset(g_fct.Gnew,  0, fb) );
    gpuErrchk( cudaMalloc((void**)&g_fct.mEff, (size_t)msh.nPlanes*sizeof(flow_float)) ); gpuErrchk( cudaMemset(g_fct.mEff, 0, (size_t)msh.nPlanes*sizeof(flow_float)) );
    gpuErrchk( cudaMalloc((void**)&g_fct.stats, (size_t)g_nPassive*8*sizeof(double)) );
    gpuErrchk( cudaMemset(g_fct.stats, 0, (size_t)g_nPassive*8*sizeof(double)) );
    gpuErrchk( cudaMalloc((void**)&g_fct.acc, (size_t)(4*g_nPassive + 4)*sizeof(double)) );
    gpuErrchk( cudaMalloc((void**)&g_fct.budget, (size_t)g_nPassive*4*sizeof(double)) ); gpuErrchk( cudaMemset(g_fct.budget, 0, (size_t)g_nPassive*4*sizeof(double)) );
    g_fct_stats_total.assign((size_t)g_nPassive*8, 0.0);
    g_fct_budget_total.assign((size_t)g_nPassive*4, 0.0); g_fct_max_relres.assign((size_t)g_nPassive, 0.0); g_fct_max_rh.assign((size_t)g_nPassive, 0.0);
    g_fct_run_max_relres.assign((size_t)g_nPassive, 0.0); g_fct_run_max_rh.assign((size_t)g_nPassive, 0.0);
    g_fct_closure_total.assign(4, 0.0);
    g_fct.nq = g_nPassive; g_fct.nCells = msh.nCells_all; g_fct.nPlanes = msh.nPlanes;
}
}  // namespace

// 設定だけで決まる判定 (checkpoint の layout と restart 契約に使う; 実行時の Pface 確保には依存しない — 2026-09-17 run_0131 の layout 不一致の修正)
bool passiveFctConfigured(const solverConfig& cfg)
{
    if (!passiveSchemeEnabled(cfg) || cfg.passiveFct == 0) return false;
    if (!(cfg.timeIntegration == 11 && cfg.unsteady == 1 && cfg.dualTime == 1)) return false;
    return cfg.speciesFaceReconstruction >= 2 && (cfg.solver == "SLAU" || cfg.solver == "SLAU2");
}
bool passiveFctActive(const solverConfig& cfg)
{
    return passiveFctConfigured(cfg) && g_Pface_dev != nullptr;
}

void passiveFctCorrect_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, flow_float a, flow_float b, flow_float c)
{
    if (!passiveFctActive(cfg)) return;
    fctAlloc(msh, var);
    { static bool once = false; if (!once) { printf("[passiveFct] active: post-step conservative FCT for %d passive scalars (prelimit %d, sweeps %d, tol %.1e)\n", g_nPassive, cfg.passiveFctPrelimit, cfg.passiveFctSweeps, cfg.passiveFctTol); once = true; } }
    const int nq = g_nPassive;
    const int nUnit = (g_qTracer >= 0) ? 1 : 0;
    const geom_int nC = msh.nCells;
    const size_t cb = (size_t)msh.nCells_all*sizeof(flow_float);
    const bool isNode = (cfg.discretization == "node");
    const bool perNode = periodicNodeActive(cfg, msh);
    const geom_int* root = perNode ? msh.periodicRoot_d : nullptr;
    flow_float* pin = isNode ? var.c_d["scalarDirichletPin"] : nullptr;
    const flow_float dt = cfg.dt;
    const flow_float oneOverDt = static_cast<flow_float>(1.0) / std::max(dt, static_cast<flow_float>(1.0e-30));
    const flow_float invA = static_cast<flow_float>(1.0) / a, cOverA = c / a;
    const int qDiff = (g_qTracer >= 0 && cfg.viscMethod != 0) ? g_qTracer : -1;
    dim3 dimGrid_nh = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));
    auto gatherAll = [&](std::vector<flow_float*>& h) { for (int q = 0; q < nq; ++q) periodicGatherArray_d_wrapper(cfg, cuda_cfg, msh, h[q]); };
    auto zeroAll = [&](std::vector<flow_float*>& h) { for (int q = 0; q < nq; ++q) gpuErrchk( cudaMemset(h[q], 0, cb) ); };
    const bool useHist = (c != static_cast<flow_float>(0.0));   // BDF2: 前 step の G/H を使う
    if (useHist && !g_fct.histValid) {
        // 契約上ここには来ない (G/H/ṁ^eff が無い restart は main が全系を BDF1 に揃える)。来たら局所形で続行し警告。
        static int warned = 0;
        if (warned < 3) { printf("[passiveFct] WARNING: flux-form history (G/H) missing for a BDF2 step: using the local form H=(V/dt)(q^n-q^{n-1}) (not limitable)\n"); ++warned; }
    }
    const flow_float* Gn = (useHist && g_fct.histValid) ? g_fct.G : nullptr;
    flow_float** Hn = (useHist && g_fct.histValid) ? g_fct.d_H : nullptr;
    // 低次作用素の質量流束: 流れの BE 形連続式と整合する ṁ^eff = invA·ṁ^{n+1} + cOverA·ṁ^eff,n (履歴が無ければ ṁ^{n+1})
    flow_float* meff = g_fct.mEff;
    {
        dim3 dimGrid_pl = dim3(ceil(msh.nPlanes / (flow_float)cuda_cfg.blocksize));
        passive_fct_meff_d<<<dimGrid_pl, cuda_cfg.dimBlock>>>(msh.nPlanes, (useHist && g_fct.histValid) ? invA : static_cast<flow_float>(1.0),
            (useHist && g_fct.histValid) ? cOverA : static_cast<flow_float>(0.0), var.p_d["massflux"], meff);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    }

    // (0) 診断: 終了状態の HO 残差から r_H のノルム (成分ごと)
    const size_t nacc = (size_t)4*nq + 4;
    gpuErrchk( cudaMemset(g_fct.acc, 0, nacc*sizeof(double)) );
    passive_fct_rh_norm_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(nC, var.c_d["volume"], dt, a, b, c, nq,
        g_p_res_dev, g_p_rophi_dev, g_fct.d_qP, g_fct.d_qPP, root, g_fct.acc + 2*nq);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    // (a) 凍結ソース S·V (モーメント: 終了状態で再評価; 周期は部分体積 → 和 gather)。トレーサは 0。
    for (int q = 0; q < nq; ++q) { gpuErrchk( cudaMemset(h_p_res[q], 0, cb) ); gpuErrchk( cudaMemset(h_p_sj[q], 0, cb) ); }
    if (g_qMom0 >= 0) {
        condensationSource_d_wrapper(cfg, cuda_cfg, msh, var);
        for (int q = g_qMom0; q < nq; ++q) periodicGatherArray_d_wrapper(cfg, cuda_cfg, msh, h_p_res[q]);
    }
    // (b) 低次作用素の対角と rhs
    gpuErrchk( cudaMemset(g_fct.diagAdv, 0, cb) ); gpuErrchk( cudaMemset(g_fct.diagDiff, 0, cb) );
    gpuErrchk( cudaMemset(g_fct.cdiff, 0, (size_t)msh.nPlanes*sizeof(flow_float)) );
    passive_fct_lo_diag_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(
        nC, msh.nNormalPlanes, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
        var.c_d["ro"], meff, isNode ? 1 : 0,
        (qDiff >= 0) ? 1 : 0, var.p_d["ge_x"], var.p_d["ge_y"], var.p_d["ge_z"],   // 面ごとの差 e (段 ③)
        var.p_d["fx"], var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
        var.c_d["vis_lam"], var.c_d["vis_turb"], static_cast<flow_float>(cfg.Sc), static_cast<flow_float>(cfg.Sc_t),
        g_fct.diagAdv, g_fct.diagDiff, g_fct.cdiff);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    periodicGatherArray_d_wrapper(cfg, cuda_cfg, msh, g_fct.diagAdv);
    periodicGatherArray_d_wrapper(cfg, cuda_cfg, msh, g_fct.diagDiff);
    passive_fct_lo_rhs_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(nC, var.c_d["volume"], dt, useHist ? invA : static_cast<flow_float>(1.0), useHist ? cOverA : static_cast<flow_float>(0.0),
        nq, g_fct.d_qP, g_fct.d_qPP, g_p_res_dev, Hn, root, g_fct.d_rhs, g_fct.acc + nq);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    // (c) Jacobi sweep (初期値 = q_H の物理限界クリップ; 受入は線形残差 ||r_L||/(||rhs|| + tolAbs) ≤ tol)
    passive_fct_lo_init_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(nC, var.c_d["ro"], nq, nUnit, g_p_rophi_dev, g_fct.d_qL);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    for (int q = 0; q < nq; ++q) periodicBroadcastArray_d_wrapper(cfg, cuda_cfg, msh, g_fct.h_qL[q]);
    std::vector<double> hacc(nacc, 0.0);
    gpuErrchk( cudaMemcpy(hacc.data(), g_fct.acc, nacc*sizeof(double), cudaMemcpyDeviceToHost) );
    std::vector<double> rhsNorm((size_t)nq), relres((size_t)nq, 0.0);
    g_fct_last_rh_rel = 0.0;
    for (int q = 0; q < nq; ++q) {
        rhsNorm[(size_t)q] = std::sqrt(hacc[(size_t)nq + q]) + (double)cfg.passiveFctTolAbs;
        const double rh = (hacc[(size_t)2*nq + 2*q + 1] > 0.0) ? std::sqrt(hacc[(size_t)2*nq + 2*q] / hacc[(size_t)2*nq + 2*q + 1]) : 0.0;
        if (!std::isfinite(rh)) g_fct_nonfinite = true;
        g_fct_last_rh_rel = std::max(g_fct_last_rh_rel, rh);
        g_fct_max_rh[(size_t)q] = std::max(g_fct_max_rh[(size_t)q], rh); g_fct_run_max_rh[(size_t)q] = std::max(g_fct_run_max_rh[(size_t)q], rh);
    }
    int sweeps = 0; double relmax = 0.0;
    for (sweeps = 1; sweeps <= cfg.passiveFctSweeps; ++sweeps) {
        zeroAll(g_fct.h_nb);
        passive_fct_lo_nb_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(
            nC, msh.nNormalPlanes, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
            var.c_d["ro"], var.c_d["roN"], meff, isNode ? 1 : 0, nq, qDiff, g_fct.cdiff, g_fct.d_qL, g_p_rophi_dev, g_fct.d_qP, g_fct.d_nb, 0);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
        gatherAll(g_fct.h_nb);
        gpuErrchk( cudaMemset(g_fct.acc, 0, (size_t)nq*sizeof(double)) );
        passive_fct_lo_solve_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
            nC, var.c_d["volume"], oneOverDt, nq, qDiff, g_fct.diagAdv, g_fct.diagDiff,
            g_fct.d_rhs, g_fct.d_nb, g_p_rophi_dev, pin, root, g_fct.d_qL, g_fct.acc);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
        for (int q = 0; q < nq; ++q) periodicBroadcastArray_d_wrapper(cfg, cuda_cfg, msh, g_fct.h_qL[q]);
        std::vector<double> r2((size_t)nq); gpuErrchk( cudaMemcpy(r2.data(), g_fct.acc, (size_t)nq*sizeof(double), cudaMemcpyDeviceToHost) );
        relmax = 0.0;
        for (int q = 0; q < nq; ++q) { relres[(size_t)q] = std::sqrt(r2[(size_t)q]) / rhsNorm[(size_t)q]; if (!std::isfinite(relres[(size_t)q])) g_fct_nonfinite = true; relmax = std::max(relmax, relres[(size_t)q]); }
        if (relmax <= cfg.passiveFctTol) break;   // 全成分が条件を満たす (plan-7 M5)
    }
    g_fct_last_sweeps = std::min(sweeps, cfg.passiveFctSweeps);
    g_fct_last_lin_res = relmax;
    for (int q = 0; q < nq; ++q) { g_fct_max_relres[(size_t)q] = std::max(g_fct_max_relres[(size_t)q], relres[(size_t)q]); g_fct_run_max_relres[(size_t)q] = std::max(g_fct_run_max_relres[(size_t)q], relres[(size_t)q]); }
    if (relmax > cfg.passiveFctTol) {
        static int warned = 0;
        if (warned < 5) { printf("[passiveFct] WARNING: low-order solve did not reach tol %.1e after %d sweeps (max component rel residual %.2e)\n", cfg.passiveFctTol, cfg.passiveFctSweeps, relmax); ++warned; }
    }
    // 密度整合の診断 (plan-7 M6): E_ρ = (V/Δt)(ρ^{n+1} − ρ^n) − Σ s ṁ^eff と、トレーサの上限条件 Lρ − f ≥ 0 の逸脱
    {
        gpuErrchk( cudaMemset(g_fct.divMeff, 0, cb) ); gpuErrchk( cudaMemset(g_fct.acc + 4*nq, 0, 4*sizeof(double)) );
        passive_fct_div_meff_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(nC, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d, meff, g_fct.divMeff);
        gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
        periodicGatherArray_d_wrapper(cfg, cuda_cfg, msh, g_fct.divMeff);
        passive_fct_density_closure_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(nC, var.c_d["volume"], dt, var.c_d["ro"], var.c_d["roN"], g_fct.divMeff, root, g_fct.acc + 4*nq);
        if (g_qTracer >= 0) {
            // Lρ の非対角: 内部近傍のみ (bndMode 1) を ρ で; f_full の境界定数項: bndMode 2 をトレーサ q で (plan-8 M3)
            gpuErrchk( cudaMemset(g_fct.nbRho, 0, cb) ); gpuErrchk( cudaMemset(g_fct.nbBnd, 0, cb) );
            passive_fct_lo_nb_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(
                nC, msh.nNormalPlanes, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
                var.c_d["ro"], var.c_d["roN"], meff, isNode ? 1 : 0, 1, (qDiff >= 0) ? 0 : -1, g_fct.cdiff, g_fct.d_roPtr, g_fct.d_roPtr, g_fct.d_roNPtr, g_fct.d_nbRho, 1);
            passive_fct_lo_nb_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(
                nC, msh.nNormalPlanes, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
                var.c_d["ro"], var.c_d["roN"], meff, isNode ? 1 : 0, 1, -1, nullptr, g_fct.d_qLtr, g_p_rophi_dev + g_qTracer, g_fct.d_qPtr, g_fct.d_nbBnd, 2);
            gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
            periodicGatherArray_d_wrapper(cfg, cuda_cfg, msh, g_fct.nbRho); periodicGatherArray_d_wrapper(cfg, cuda_cfg, msh, g_fct.nbBnd);
            passive_fct_upper_margin_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(nC, var.c_d["volume"], oneOverDt, var.c_d["ro"], g_fct.diagAdv, (qDiff >= 0) ? g_fct.diagDiff : nullptr,
                g_qTracer, g_fct.nbRho, g_fct.nbBnd, g_fct.d_rhs, pin, root, g_fct.acc + 4*nq + 2);
            gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
        }
        double hc[4] = {0.0, 0.0, 0.0, 0.0}; gpuErrchk( cudaMemcpy(hc, g_fct.acc + 4*nq, 4*sizeof(double), cudaMemcpyDeviceToHost) );
        const double erel = (hc[1] > 0.0) ? std::sqrt(hc[0]/hc[1]) : 0.0;
        g_fct_closure_total[0] = std::max(g_fct_closure_total[0], erel); g_fct_closure_total[1] += hc[2]*(double)dt; g_fct_closure_total[2] += hc[3];
        if (g_qTracer >= 0) g_fct_stats_total[(size_t)g_qTracer*8 + 7] += hc[2]*(double)dt;   // 上限条件の逸脱量 (·Δt) を FAIL 判定へ (plan-8 M3)
    }
    // (d) A^raw (全輸送面) と基点 q_B
    zeroAll(g_fct.h_base); zeroAll(g_fct.h_Pp); zeroAll(g_fct.h_Pm); zeroAll(g_fct.h_corr);
    gpuErrchk( cudaMemset(g_fct.stats, 0, (size_t)nq*8*sizeof(double)) );
    passive_fct_raw_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(
        nC, msh.nNormalPlanes, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
        var.c_d["ro"], var.c_d["roN"], var.p_d["massflux"], meff, isNode ? 1 : 0, useHist ? invA : static_cast<flow_float>(1.0), useHist ? cOverA : static_cast<flow_float>(0.0),
        nq, nq, g_Pface_dev, qDiff, (qDiff >= 0) ? g_fct.cdiff : nullptr, Gn,
        g_p_rophi_dev, g_fct.d_qL, g_fct.d_qP, g_fct.Aface, g_fct.d_base);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    gatherAll(g_fct.h_base);
    passive_fct_base_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(nC, var.c_d["volume"], dt, var.c_d["ro"], nq, nUnit,
        g_p_rophi_dev, g_fct.d_base, root, g_fct.d_qB, g_fct.stats);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    // (e) 前制限 → P± → 局所極値 → R±
    passive_fct_prelimit_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(nC, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d, var.c_d["ro"],
        nq, nq, g_fct.Aface, cfg.passiveFctPrelimit, g_fct.d_qB, g_fct.Apre, g_fct.d_Pp, g_fct.d_Pm, g_fct.stats);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    gatherAll(g_fct.h_Pp); gatherAll(g_fct.h_Pm);
    passive_fct_extrema_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        nC, msh.nNormalPlanes, msh.map_plane_cells_d, msh.map_cell_planes_index_d, msh.map_cell_planes_d,
        var.c_d["ro"], var.c_d["roN"], nq, g_fct.d_qP, g_fct.d_qB, g_fct.d_mx, g_fct.d_mn);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    for (int q = 0; q < nq; ++q) { periodicGatherMaxArray_d_wrapper(cfg, cuda_cfg, msh, g_fct.h_mx[q]); periodicGatherMinArray_d_wrapper(cfg, cuda_cfg, msh, g_fct.h_mn[q]); }
    passive_fct_ratio_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        nC, var.c_d["volume"], oneOverDt, var.c_d["ro"], nq, nUnit,
        g_fct.d_qB, g_fct.d_mx, g_fct.d_mn, g_fct.d_Pp, g_fct.d_Pm, pin, g_fct.d_Rp, g_fct.d_Rm);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    // (f) α_f, 補正 (独立バッファ → gather → commit → mirror), 実現流束 G^{n+1}
    passive_fct_apply_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(nC, msh.nNormalPlanes, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
        var.c_d["ro"], var.c_d["roN"], meff, isNode ? 1 : 0, nq, nq, g_fct.Aface, g_fct.Apre, g_fct.d_Rp, g_fct.d_Rm, 0,
        qDiff, (qDiff >= 0) ? g_fct.cdiff : nullptr, g_fct.d_qL, g_fct.d_qP, g_fct.d_corr, g_fct.Gnew, g_fct.stats);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    gatherAll(g_fct.h_corr);
    passive_fct_commit_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(nC, var.c_d["volume"], dt, nq, g_fct.d_corr, pin, root, g_p_rophi_dev, g_fct.stats);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    for (int q = 0; q < nq; ++q) periodicBroadcastArray_d_wrapper(cfg, cuda_cfg, msh, h_p_rophi[q]);
    std::vector<double> st((size_t)nq*8);
    gpuErrchk( cudaMemcpy(st.data(), g_fct.stats, st.size()*sizeof(double), cudaMemcpyDeviceToHost) );
    for (int q = 0; q < nq; ++q) {
        double* t = g_fct_stats_total.data() + (size_t)q*8; const double* f = st.data() + (size_t)q*8;
        t[0] += f[0]*(double)dt; t[1] += f[1]; t[2] += f[2]*(double)dt; t[3] += f[3]; t[4] += f[4]; t[5] += f[5]*(double)dt; t[6] += f[6]*(double)dt;
    }
}

void passiveFctFinishHistory_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, flow_float a, flow_float b, flow_float c)
{
    (void)b;
    if (!passiveFctActive(cfg)) return;
    fctAlloc(msh, var);
    const int nq = g_nPassive; const geom_int nC = msh.nCells;
    const size_t cb = (size_t)msh.nCells_all*sizeof(flow_float);
    const bool useHist = (c != static_cast<flow_float>(0.0)) && g_fct.histValid;
    const flow_float invA = useHist ? static_cast<flow_float>(1.0)/a : static_cast<flow_float>(1.0), cOverA = useHist ? c/a : static_cast<flow_float>(0.0);
    const geom_int* root = periodicNodeActive(cfg, msh) ? msh.periodicRoot_d : nullptr;
    dim3 dimGrid_nh = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));
    std::swap(g_fct.G, g_fct.Gnew);   // G ← 実現流束 G^{n+1}
    for (int q = 0; q < nq; ++q) gpuErrchk( cudaMemset(g_fct.h_base[q], 0, cb) );   // divG の一時
    passive_fct_divG_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(nC, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d, nq, nq, g_fct.G, g_fct.d_base);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    for (int q = 0; q < nq; ++q) periodicGatherArray_d_wrapper(cfg, cuda_cfg, msh, g_fct.h_base[q]);
    gpuErrchk( cudaMemset(g_fct.budget, 0, (size_t)nq*4*sizeof(double)) );
    // res_* は passiveFctCorrect (a) で組んだ凍結ソース (gather 済み) のまま → H_src の再帰に使う
    passive_fct_hist_local_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(nC, var.c_d["volume"], cfg.dt, invA, cOverA, nq,
        g_p_rophi_dev, g_fct.d_qP, g_fct.d_base, g_p_res_dev, g_fct.d_Hsrc, g_fct.d_H, root, g_fct.budget);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    std::vector<double> hb((size_t)nq*4); gpuErrchk( cudaMemcpy(hb.data(), g_fct.budget, hb.size()*sizeof(double), cudaMemcpyDeviceToHost) );
    for (size_t i = 0; i < hb.size(); ++i) { if (!std::isfinite(hb[i])) g_fct_nonfinite = true; g_fct_budget_total[i] += hb[i]; }
    g_fct.histValid = true;
}

bool passiveFctHistoryToHost(const mesh& msh, std::vector<std::vector<flow_float>>& G, std::vector<std::vector<flow_float>>& H, std::vector<flow_float>& mEff, std::vector<std::vector<flow_float>>& Hsrc)
{
    if (!g_fct.histValid || g_fct.nq == 0) return false;
    Hsrc.assign((size_t)g_fct.nq, std::vector<flow_float>((size_t)msh.nCells));
    for (int q = 0; q < g_fct.nq; ++q) gpuErrchk( cudaMemcpy(Hsrc[(size_t)q].data(), g_fct.h_Hsrc[(size_t)q], (size_t)msh.nCells*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    mEff.assign((size_t)msh.nPlanes, 0.0f);
    gpuErrchk( cudaMemcpy(mEff.data(), g_fct.mEff, mEff.size()*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    G.assign((size_t)g_fct.nq, std::vector<flow_float>((size_t)msh.nPlanes));
    H.assign((size_t)g_fct.nq, std::vector<flow_float>((size_t)msh.nCells));
    std::vector<flow_float> gall((size_t)msh.nPlanes*g_fct.nq);
    gpuErrchk( cudaMemcpy(gall.data(), g_fct.G, gall.size()*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    for (geom_int ip = 0; ip < msh.nPlanes; ++ip) for (int q = 0; q < g_fct.nq; ++q) G[(size_t)q][(size_t)ip] = gall[(size_t)ip*g_fct.nq + q];
    for (int q = 0; q < g_fct.nq; ++q) gpuErrchk( cudaMemcpy(H[(size_t)q].data(), g_fct.h_H[(size_t)q], (size_t)msh.nCells*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    return true;
}

bool passiveFctHistoryFromHost(solverConfig& cfg, mesh& msh, variables& var, const std::vector<std::vector<flow_float>>& G, const std::vector<std::vector<flow_float>>& H, const std::vector<flow_float>& mEff, const std::vector<std::vector<flow_float>>& Hsrc)
{
    // restart 時点では Pface が未確保なので実行時判定 (passiveFctActive) は使えない (検証 2 巡目 B: 履歴が一度も復元されなかった真因)。設定だけで判定する。
    if (!passiveFctConfigured(cfg)) return false;
    fctAlloc(msh, var);
    if ((int)G.size() != g_fct.nq || (int)H.size() != g_fct.nq || (int)Hsrc.size() != g_fct.nq || (geom_int)mEff.size() < msh.nPlanes) return false;
    for (int q = 0; q < g_fct.nq; ++q) { if ((geom_int)Hsrc[(size_t)q].size() < msh.nCells) return false; gpuErrchk( cudaMemcpy(g_fct.h_Hsrc[(size_t)q], Hsrc[(size_t)q].data(), (size_t)msh.nCells*sizeof(flow_float), cudaMemcpyHostToDevice) ); }
    gpuErrchk( cudaMemcpy(g_fct.mEff, mEff.data(), (size_t)msh.nPlanes*sizeof(flow_float), cudaMemcpyHostToDevice) );
    std::vector<flow_float> gall((size_t)msh.nPlanes*g_fct.nq, 0.0f);
    for (int q = 0; q < g_fct.nq; ++q) {
        if ((geom_int)G[(size_t)q].size() < msh.nPlanes || (geom_int)H[(size_t)q].size() < msh.nCells) return false;
        for (geom_int ip = 0; ip < msh.nPlanes; ++ip) gall[(size_t)ip*g_fct.nq + q] = G[(size_t)q][(size_t)ip];
        gpuErrchk( cudaMemcpy(g_fct.h_H[(size_t)q], H[(size_t)q].data(), (size_t)msh.nCells*sizeof(flow_float), cudaMemcpyHostToDevice) );
    }
    gpuErrchk( cudaMemcpy(g_fct.G, gall.data(), gall.size()*sizeof(flow_float), cudaMemcpyHostToDevice) );
    g_fct.histValid = true;
    return true;
}

std::vector<double> passiveFctStatsTotals() { return g_fct_stats_total; }
double passiveFctLastLinRes() { return g_fct_last_lin_res; }
int    passiveFctLastSweeps() { return g_fct_last_sweeps; }
double passiveFctLastRhRel()  { return g_fct_last_rh_rel; }

// =============================================================================
// 二相拡散 (plans/active/condensation-two-phase-transport.md §4.2, §5.1 #4e; methods/condensation.md §7c)。
//   TP carrier 凝縮で、気相内の分子拡散 (z 基準・風上補正) と全輸送量共通の乱流拡散を面ごとに 1 回だけ組み、
//   化学種 (水は総水分)・液 ρg・Q2/Q1/Q0・エネルギーの残差と点対角に同じ面流束を足す。面の代数は twoPhaseDiffusion_d.cuh の
//   tp_face_flux (GPU 単体試験 tests/unit/test_twophase_kernel.cu と共用)。呼び出しは凝縮モーメントの残差ゼロ化 (passiveAdvection) の後、
//   境界ピン・周期集約の前 (condensationTransport_d_wrapper の末尾)。
// =============================================================================
namespace {
// 面の入力 (格納値と面の係数) を組む。本番カーネル twophase_diffusion_d と収束受入の監査 twophase_audit_face_d (#4f) が共有する。
// node 境界半割面 (流束 0) は false。
__device__ inline bool tp_build_face_in(
    geom_int nCells, geom_int ip, geom_int* plane_cells,
    const flow_float* ge_x, const flow_float* ge_y, const flow_float* ge_z,   // 面ごとの差 e = cc[ic1] − cc[ic0] (§4.2a、段 ③)
    geom_float* fx, geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,
    const SpeciesThermoF* sp, int nSpecies, int iw,
    flow_float** roY, flow_float** rophi,
    flow_float* ro, flow_float* T, flow_float* P, flow_float* vis_lam, flow_float* vis_turb,
    int diffMethod, flow_float Sc, flow_float Sc_t, int isNode, const CondSpeciesProps& cprops,
    TpFaceIn& in, geom_int& ic0, geom_int& ic1)
{
    ic0 = plane_cells[2 * ip + 0];
    ic1 = plane_cells[2 * ip + 1];
    // node 境界半割面は流束 0 (species_diffusion_d と同方針: Dirichlet はピン、Neumann は流束 0)
    if (isNode != 0 && (ic0 >= nCells || ic1 >= nCells)) return false;

    const flow_float f   = fx[ip];
    const flow_float sxx = sx[ip], syy = sy[ip], szz = sz[ip], sss = ss[ip];
    // 座標の差 cc[ic1] − cc[ic0] の代わりに面ごとの差 e (同じ向き) を読む (§4.2a、段 ③)
    const flow_float dccx = ge_x[ip];
    const flow_float dccy = ge_y[ip];
    const flow_float dccz = ge_z[ip];
    const flow_float dcc  = sqrtf(dccx*dccx + dccy*dccy + dccz*dccz);
    const flow_float denom = dccx*sxx + dccy*syy + dccz*szz;
    const flow_float Dsafe = (fabsf(denom) < 1.0e-30f) ? ((denom>=0.0f)?1.0e-30f:-1.0e-30f) : denom;
    const flow_float delta = dcc * sss * sss / Dsafe;       // over-relaxed 法線 (species_diffusion_d と同じ)

    in.n = nSpecies; in.iw = iw; in.f = f;
    in.rho0 = max(ro[ic0], (flow_float)1.0e-30f);
    in.rho1 = max(ro[ic1], (flow_float)1.0e-30f);
    for (int s = 0; s < nSpecies; ++s) { in.rY0[s] = roY[s][ic0]; in.rY1[s] = roY[s][ic1]; }
    in.rg0 = rophi[0][ic0]; in.rg1 = rophi[0][ic1];
    for (int m = 0; m < TP_NQ; ++m) { in.rQ0[m] = rophi[1+m][ic0]; in.rQ1[m] = rophi[1+m][ic1]; }
    in.geo = delta/dcc;
    in.geo_abs = fabsf(delta)/max(dcc, (flow_float)1.0e-30f);

    const flow_float g = 1.0f - f;
    const flow_float T_face = f*T[ic0] + g*T[ic1];
    const flow_float P_face = f*P[ic0] + g*P[ic1];
    // 拡散係数の組成: species_diffusion_d と同じ (面の総組成を clip・正規化し、液 g_f を除いた気相組成 X; §4.1 / #3)
    const flow_float inv_ro0 = 1.0f/in.rho0, inv_ro1 = 1.0f/in.rho1;
    flow_float Yf[THERMO_MAX_SPECIES], X[THERMO_MAX_SPECIES];
    flow_float ysum = 0.0f;
    for (int s = 0; s < nSpecies; ++s) {
        flow_float y = f*in.rY0[s]*inv_ro0 + g*in.rY1[s]*inv_ro1;
        if (y < 0.0f) y = 0.0f; Yf[s] = y; ysum += y;
    }
    const flow_float yinv = 1.0f/(ysum > 1.0e-30f ? ysum : 1.0e-30f);
    for (int s = 0; s < nSpecies; ++s) Yf[s] *= yinv;
    const flow_float gl_f = (f*in.rg0*inv_ro0 + g*in.rg1*inv_ro1)*yinv;
    species_transport_X_f(sp, nSpecies, Yf, iw, gl_f, X);
    const flow_float mu_face  = f*vis_lam[ic0]  + g*vis_lam[ic1];
    const flow_float mut_face = f*vis_turb[ic0] + g*vis_turb[ic1];
    // 定数 Schmidt (diffMethod 0) は ρ_g,f D = μ_f/Sc (気相の分子拡散)
    const flow_float rgf = f*fmaxf(in.rho0 - in.rg0, 1.0e-30f) + g*fmaxf(in.rho1 - in.rg1, 1.0e-30f);
    for (int s = 0; s < nSpecies; ++s) {
        in.D[s] = (diffMethod == 1) ? thermo_Dmix_species_f(sp, nSpecies, X, s, T_face, P_face) : mu_face/(rgf*Sc);
        in.h[s] = thermo_h_mass_f(sp[s], T_face);
    }
    in.ct = (mut_face > 0.0f) ? mut_face/Sc_t : 0.0f;
    in.L  = (flow_float)cond_latent(cprops, (double)T_face);   // EOS と同じ潜熱 (double を 1 回評価)

    return true;
}

__global__ void twophase_diffusion_d(
    geom_int nCells, geom_int nNormalHaloPlanes, geom_int* normal_halo_planes, geom_int* plane_cells,
    const flow_float* ge_x, const flow_float* ge_y, const flow_float* ge_z,   // 面ごとの差 e (段 ③)
    geom_float* fx, geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,
    const SpeciesThermoF* sp, int nSpecies, int iw,
    flow_float** roY, flow_float** res_roY, flow_float** transport_diag,
    flow_float** rophi, flow_float** res_phi, flow_float** diag_phi,   // 受動種ポインタ配列のモーメント先頭 (順序 g, Q2, Q1, Q0)
    flow_float* ro, flow_float* T, flow_float* P, flow_float* vis_lam, flow_float* vis_turb,
    flow_float* res_roe, int diffMethod, flow_float Sc, flow_float Sc_t, int isNode, CondSpeciesProps cprops,
    TpoFacePtr tpo = TpoFacePtr{})   // 診断 G3-a (既定は無効): res[ic0] に足した J・Jl・JQ を面スロットへ写す
{
    const geom_int ih = blockDim.x * blockIdx.x + threadIdx.x;
    if (ih >= nNormalHaloPlanes) return;
    TpFaceIn in; geom_int ic0, ic1;
    if (!tp_build_face_in(nCells, normal_halo_planes[ih], plane_cells, ge_x, ge_y, ge_z, fx, sx, sy, sz, ss, sp, nSpecies, iw,
                          roY, rophi, ro, T, P, vis_lam, vis_turb, diffMethod, Sc, Sc_t, isNode, cprops, in, ic0, ic1)) {
        if (tpo.val != nullptr) {   // 診断 G3-a: node 境界半割面 (足さない) も明示的に 0 を書く
            tpo_face_code(tpo, ih, 2);
            for (int s = 0; s < nSpecies; ++s) tpo_face_put(tpo, ih, tpo_comp(tpo, s), 0.0);
            for (int m = 0; m < 1 + TP_NQ; ++m) tpo_face_put(tpo, ih, 1 + m, 0.0);
        }
        return;
    }
    TpFaceOut o;
    tp_face_flux(in, o);
    if (tpo.val != nullptr) {   // 診断 G3-a: atomic の前の値 (res[ic0] に足す値; res[ic1] には符号反転)
        tpo_face_code(tpo, ih, 1);
        for (int s = 0; s < nSpecies; ++s) tpo_face_put(tpo, ih, tpo_comp(tpo, s), (double)o.J[s]);
        tpo_face_put(tpo, ih, 1, (double)o.Jl);
        for (int m = 0; m < TP_NQ; ++m) tpo_face_put(tpo, ih, 2 + m, (double)o.JQ[m]);
    }

    if (ic0 < nCells) {
        for (int s = 0; s < nSpecies; ++s) { atomicAdd(&res_roY[s][ic0], o.J[s]); atomicAdd(&transport_diag[s][ic0], o.diag0[s]); }
        atomicAdd(&res_phi[0][ic0], o.Jl); atomicAdd(&diag_phi[0][ic0], o.diagt0);
        for (int m = 0; m < TP_NQ; ++m) { atomicAdd(&res_phi[1+m][ic0], o.JQ[m]); atomicAdd(&diag_phi[1+m][ic0], o.diagt0); }
        atomicAdd(&res_roe[ic0], o.q);
    }
    if (ic1 < nCells) {
        for (int s = 0; s < nSpecies; ++s) { atomicAdd(&res_roY[s][ic1], -o.J[s]); atomicAdd(&transport_diag[s][ic1], o.diag1[s]); }
        atomicAdd(&res_phi[0][ic1], -o.Jl); atomicAdd(&diag_phi[0][ic1], o.diagt1);
        for (int m = 0; m < TP_NQ; ++m) { atomicAdd(&res_phi[1+m][ic1], -o.JQ[m]); atomicAdd(&diag_phi[1+m][ic1], o.diagt1); }
        atomicAdd(&res_roe[ic1], -o.q);
    }
}

// 二相拡散の更新の再正規化: species_renormalize_d と同じ係数 ρ/ΣρY を ρY_s に掛け、同じ係数を液 ρg と Q にも掛ける
// (plan §4.2「化学種の再正規化で ρY_w に掛けた係数を ρg と Q にも」; 現行は液更新の前に水だけ再正規化される)。
__global__ void species_renormalize_liquid_d(
    geom_int nCells, int nSpecies, flow_float** roY, flow_float* ro,
    double* reasons, int iw, const geom_float* vol, const geom_int* root, flow_float** rophi, double* rnAcc,
    double* diagRn,   // 診断 (condTwoPhaseDiag; nullptr で書かない): [ic] f−1, [nCells+ic] Δ(ρY_w), [2nCells+ic] Δ(ρg)
    double* tpu = nullptr, int iwTpu = -1)   // 診断 G3-b (既定 nullptr): 水の負値の 0 化・係数の乗算 (水・液・Q) を before/after で記録
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const double qw_m = ((rnAcc != nullptr || diagRn != nullptr) && iw >= 0) ? (double)roY[iw][ic] : 0.0;   // 再正規化の直前 (負値の 0 化の前)
    const flow_float tpu_w0 = (tpu != nullptr && iwTpu >= 0) ? roY[iwTpu][ic] : (flow_float)0.0;
    double sum = 0.0;
    for (int s = 0; s < nSpecies; s++) {
        flow_float v = roY[s][ic];
        if (v < 0.0) v = 0.0;
        roY[s][ic] = v;
        sum += (double)v;
    }
    const double factor = (double)ro[ic] / (sum > (double)kSmall ? sum : (double)kSmall);
    const flow_float yw_in = (reasons != nullptr && iw >= 0) ? roY[iw][ic] : (flow_float)0.0;
    const flow_float tpu_w1 = (tpu != nullptr && iwTpu >= 0) ? roY[iwTpu][ic] : (flow_float)0.0;
    for (int s = 0; s < nSpecies; s++) roY[s][ic] = (flow_float)((double)roY[s][ic] * factor);
    double dl[1 + TP_NQ], qin[1 + TP_NQ];
    for (int m = 0; m < 1 + TP_NQ; ++m) { const flow_float in_ = rophi[m][ic]; qin[m] = (double)in_; rophi[m][ic] = (flow_float)((double)in_ * factor); dl[m] = fabs((double)rophi[m][ic] - (double)in_); }
    if (tpu != nullptr && iwTpu >= 0) {   // 成分 1..4 = 液・Q2・Q1・Q0 (rophi の順と同じ)
        tpu_ba(tpu, TPU_OP_RNF, 0, nCells, ic, (double)tpu_w0, (double)tpu_w1);
        tpu_ba(tpu, TPU_OP_RNS, 0, nCells, ic, (double)tpu_w1, (double)roY[iwTpu][ic]);
        for (int m = 0; m < 1 + TP_NQ; ++m) tpu_ba(tpu, TPU_OP_RNS, 1 + m, nCells, ic, qin[m], (double)rophi[m][ic]);
        tpu_put(tpu, TPU_RN_FACTOR, nCells, ic, factor);
    }
    if (diagRn != nullptr) {
        diagRn[ic] = factor - 1.0;
        diagRn[(size_t)nCells + ic] = (iw >= 0) ? ((double)roY[iw][ic] - qw_m) : 0.0;
        diagRn[2*(size_t)nCells + ic] = (double)rophi[0][ic] - qin[0];
    }
    if (reasons != nullptr && iw >= 0 && (root == nullptr || root[ic] == ic)) {
        const double d = fabs((double)roY[iw][ic] - (double)yw_in);
        const double Vc = (vol != nullptr) ? (double)vol[ic] : 1.0;
        // 成分別の補正量 (#1b-pre): 液と Q にも同じ係数を掛けた量
        if (dl[0] > 0.0) atomicAdd(&reasons[COND_REASON_RNG_SUM], dl[0]*Vc);
        if (dl[1] > 0.0) atomicAdd(&reasons[COND_REASON_RNQ2_SUM], dl[1]*Vc);
        if (dl[2] > 0.0) atomicAdd(&reasons[COND_REASON_RNQ1_SUM], dl[2]*Vc);
        if (dl[3] > 0.0) atomicAdd(&reasons[COND_REASON_RNQ0_SUM], dl[3]*Vc);
        cond_atomic_max_double(&reasons[COND_REASON_RN_MAXC], fabs(factor - 1.0));
        cond_atomic_max_double(&reasons[COND_REASON_RN_MAX], fabs(factor - 1.0));
        if (rnAcc != nullptr) {
            double qm[RNG_NC], qp[RNG_NC];
            qm[0] = qw_m; qp[0] = (double)roY[iw][ic];
            for (int m = 0; m < RNG_NC - 1; ++m) { qm[1+m] = qin[m]; qp[1+m] = (double)rophi[m][ic]; }
            rng_accumulate(rnAcc, qm, qp, Vc, fabs(factor - 1.0));
        }
        if (d > 0.0) {
            atomicAdd(&reasons[COND_REASON_RN_SUM], d*((vol != nullptr) ? (double)vol[ic] : 1.0));
            atomicAdd(&reasons[COND_REASON_RN_N], 1.0);
        }
    }
}
}  // namespace

void twoPhaseDiffusion_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!condTwoPhaseDiffusionActive(cfg) || !speciesEnabled(var) || g_roY_dev == nullptr || g_qMom0 < 0) return;
    const CondSpeciesProps cprops = condProps_make(cfg.condModel, cond_prop_opts(cfg));
    dim3 dimGrid_nh = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));
    twophase_diffusion_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(
        msh.nCells, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
        var.p_d["ge_x"], var.p_d["ge_y"], var.p_d["ge_z"],
        var.p_d["fx"], var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
        thermo_species_device_ptr_f(), g_nSpecies, cfg.condGasSpecies,
        g_roY_dev, g_resroY_dev, g_transdiag_dev,
        g_p_rophi_dev + g_qMom0, g_p_res_dev + g_qMom0, g_p_diag_dev + g_qMom0,
        var.c_d["ro"], var.c_d["T"], var.c_d["P"], var.c_d["vis_lam"], var.c_d["vis_turb"],
        var.c_d["res_roe"], cfg.speciesDiffusionMethod, cfg.Sc, cfg.Sc_t,
        (cfg.discretization == "node") ? 1 : 0, cprops,
        tpoFace(TPO_K_DIFF_ON, TPO_KIND_DIFF, -1));   // 診断 G3-a (既定は無効; 化学種は化学種の成分、液・Q は 1..4)
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

void speciesRenormalizeTwoPhase_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!speciesEnabled(var) || g_roY_dev == nullptr || g_qMom0 < 0) return;
    int iwRn = -1;
    double* reasonsRn = condCorrReasonsForRenormalize(cfg, var, &iwRn);
    double* rnAcc = (reasonsRn != nullptr) ? rngBegin() : nullptr;
    species_renormalize_liquid_d<<<cuda_cfg.dimGrid_cell, cuda_cfg.dimBlock>>>(
        msh.nCells, g_nSpecies, g_roY_dev, var.c_d["ro"],
        reasonsRn, iwRn, var.c_d["volume"], passive_periodic_root(cfg, msh), g_p_rophi_dev + g_qMom0, rnAcc, twoPhaseDiagRenormPtr(),
        tpuSlots(), cfg.condGasSpecies);   // 診断 G3-b (既定 nullptr)
    if (rnAcc != nullptr) rngEnd(cuda_cfg);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}

// =============================================================================
// 収束受入の独立残差監査 (plans/active/condensation-two-phase-transport.md §5.1 #4f (4); 設計メモ §15)。
//   格納状態から呼び出し側が assembleResidual を 1 回回して EOS (double 研磨の温度反転)・面値 (massflux と S3 の面組成)・
//   気相物性を作り直した後、二相系の各成分 (化学種 [水は総水分 w]・蒸気 v = w − l・液 l・Q2/Q1/Q0) の定常残差を
//   **面ごとの流束を double で組み直して double で集計**する (float の atomicAdd 残差や R_w − R_g の float 減算を使わない)。
//   移流は ṁ·φ_face (massflux・面組成は流れのソルバが作った float の値を入力として扱う)、拡散は tp_face_flux<double>
//   (格納値を double に上げる; D_k・h_k・L・μ_t の係数は本番と同じ float 値)、相変化ソースは同じソースカーネルの値 (セルごとの項)。
//   判定は成分ごと: max_i |r_q,i| ≤ max(1e-7·r0_q, 6ε₃₂·max_i A_q,i) (r0 = 計算開始時の監査、A = Σ|移流| + Σ|拡散| + |S|V)。
//   エネルギーは流れの全残差 (対流・粘性・熱伝導) を double で組み直していないので対象外 (rms_roe は check_convergence が見る)。
// =============================================================================
namespace {
__global__ void twophase_audit_face_d(
    geom_int nCells, geom_int nNormalHaloPlanes, geom_int* normal_halo_planes, geom_int* plane_cells,
    const flow_float* ge_x, const flow_float* ge_y, const flow_float* ge_z,   // 面ごとの差 e (本番と同じ、段 ③)
    geom_float* fx, geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,
    const SpeciesThermoF* sp, int nSpecies, int iw, flow_float** roY, flow_float** rophi,
    flow_float* ro, flow_float* T, flow_float* P, flow_float* vis_lam, flow_float* vis_turb,
    int diffMethod, flow_float Sc, flow_float Sc_t, int isNode, CondSpeciesProps cprops,
    const flow_float* massflux, const flow_float* Yface, const flow_float* Pface, int nPassive, int q0,
    int NQ, double* r, double* A, int oldOp, int visc, GasPhaseLiquid liq,
    double* brk = nullptr)   // #1b-r2: 液 R_g の内訳 [c*3+0] 移流, [c*3+1] 拡散 (セルへ入る向き; nullptr で書かない)
{
    const geom_int ih = blockDim.x * blockIdx.x + threadIdx.x;
    if (ih >= nNormalHaloPlanes) return;
    const geom_int ip  = normal_halo_planes[ih];
    const geom_int c0 = plane_cells[2 * ip + 0];
    const geom_int c1 = plane_cells[2 * ip + 1];
    const double mdot = (double)massflux[ip];
    const bool nodeBnd = (isNode != 0 && c1 >= nCells);
    const flow_float inv0 = 1.0f/max(ro[c0], (flow_float)1.0e-30);   // 境界ノード自身の組成 (species_advection_faceY_d と同じ float)
    auto add = [&](int q, double F) {   // F: セル 0 へ入る向き
        if (c0 < nCells) { atomicAdd(&r[(size_t)c0*NQ + q],  F); atomicAdd(&A[(size_t)c0*NQ + q], fabs(F)); }
        if (c1 < nCells) { atomicAdd(&r[(size_t)c1*NQ + q], -F); atomicAdd(&A[(size_t)c1*NQ + q], fabs(F)); }
    };
    // 移流 (S3 面値; node 境界半割面は境界ノード自身の値)
    double Yw_f = 0.0, g_f = 0.0;
    for (int s = 0; s < nSpecies; ++s) {
        const double Yf = nodeBnd ? (double)(roY[s][c0]*inv0) : (double)Yface[(size_t)ip*nSpecies + s];
        if (s == iw) Yw_f = Yf;
        add(s, -mdot*Yf);
    }
    for (int m = 0; m < 1 + TP_NQ; ++m) {
        const double Pf = nodeBnd ? (double)(rophi[m][c0]*inv0) : (double)Pface[(size_t)ip*nPassive + q0 + m];
        if (m == 0) g_f = Pf;
        add(nSpecies + 1 + m, -mdot*Pf);
        if (m == 0 && brk != nullptr) {
            if (c0 < nCells) atomicAdd(&brk[(size_t)c0*3 + 0], -mdot*Pf);
            if (c1 < nCells) atomicAdd(&brk[(size_t)c1*3 + 0],  mdot*Pf);
        }
    }
    add(nSpecies, -mdot*(Yw_f - g_f));   // 蒸気
    if (oldOp) {
        // 旧作用素 (二相拡散 OFF; #1b-pre (1)): species_diffusion_d の式 — 総組成 Y_s を ρ_f(D_s + D_t) で Fick 拡散、面の算術平均 Y_f で ΣJ=0 補正。
        //   液とモーメントは拡散しない。係数 (ρ_f, D_s, D_t) は本番と同じ float の評価順で作り、Y の差と流束・補正は double。
        if (visc == 0) return;
        if (isNode != 0 && (c0 >= nCells || c1 >= nCells)) return;
        const flow_float f = fx[ip];
        const flow_float sxx = sx[ip], syy = sy[ip], szz = sz[ip], sss = ss[ip];
        const flow_float dccx = ge_x[ip], dccy = ge_y[ip], dccz = ge_z[ip];   // = cc[c1] − cc[c0] の面ごとの差 e (段 ③)
        const flow_float dcc  = sqrtf(dccx*dccx + dccy*dccy + dccz*dccz);
        const flow_float denom = dccx*sxx + dccy*syy + dccz*szz;
        const flow_float Dsafe = (fabsf(denom) < 1.0e-30f) ? ((denom>=0.0f)?1.0e-30f:-1.0e-30f) : denom;
        const flow_float delta = dcc * sss * sss / Dsafe;
        const flow_float ro0 = max(ro[c0], (flow_float)1.0e-30f), ro1 = max(ro[c1], (flow_float)1.0e-30f);
        const flow_float inv_ro0 = 1.0f/ro0, inv_ro1 = 1.0f/ro1;
        const flow_float g = 1.0f - f;
        const flow_float ro_face = f*ro0 + g*ro1;
        const flow_float T_face  = f*T[c0] + g*T[c1];
        const flow_float P_face  = f*P[c0] + g*P[c1];
        flow_float Yf[THERMO_MAX_SPECIES], X[THERMO_MAX_SPECIES], Ys0[THERMO_MAX_SPECIES], Ys1[THERMO_MAX_SPECIES];
        flow_float ysum = 0.0f;
        for (int s = 0; s < nSpecies; s++) {
            Ys0[s] = roY[s][c0]*inv_ro0; Ys1[s] = roY[s][c1]*inv_ro1;
            flow_float y = f*Ys0[s] + g*Ys1[s];
            if (y < 0.0f) y = 0.0f; Yf[s] = y; ysum += y;
        }
        const flow_float yinv = 1.0f/(ysum>1.0e-30f?ysum:1.0e-30f);
        for (int s = 0; s < nSpecies; s++) Yf[s] *= yinv;
        const flow_float gl_f = (liq.iw >= 0 && liq.rog != nullptr) ? (f*liq.rog[c0]*inv_ro0 + g*liq.rog[c1]*inv_ro1)*yinv : 0.0f;
        species_transport_X_f(sp, nSpecies, Yf, liq.iw, gl_f, X);
        const flow_float mu_face  = f*vis_lam[c0]  + g*vis_lam[c1];
        const flow_float mut_face = f*vis_turb[c0] + g*vis_turb[c1];
        const flow_float Dt = (mut_face > 0.0f) ? mut_face/(ro_face*Sc_t) : 0.0f;
        const double geo = (double)delta/(double)dcc;
        // 駆動差と補正は double (格納 ρY を double に上げる)
        const double r0d = (double)ro0, r1d = (double)ro1;
        double Js[THERMO_MAX_SPECIES], sumJ = 0.0, Yd[THERMO_MAX_SPECIES], ys = 0.0;
        for (int s = 0; s < nSpecies; s++) {
            flow_float D = (diffMethod == 1) ? thermo_Dmix_species_f(sp, nSpecies, X, s, T_face, P_face) : mu_face/(ro_face*Sc);
            D += Dt;
            const flow_float roD = ro_face*D;   // 係数 (本番と同じ float)
            const double y0 = (double)roY[s][c0]/r0d, y1 = (double)roY[s][c1]/r1d;
            Js[s] = (double)roD*(y1 - y0)*geo;
            sumJ += Js[s];
            double y = (double)f*y0 + (1.0 - (double)f)*y1; if (y < 0.0) y = 0.0; Yd[s] = y; ys += y;
        }
        const double yi = 1.0/(ys > 1.0e-30 ? ys : 1.0e-30);
        for (int s = 0; s < nSpecies; s++) {
            const double Jc = Js[s] - Yd[s]*yi*sumJ;
            add(s, Jc);
            if (s == iw) add(nSpecies, Jc);   // 蒸気 = 総水分 (液は拡散しない; 内訳の拡散は 0)
        }
        return;
    }
    // 拡散 (格納値を double に上げて同じ式)
    TpFaceIn in; geom_int ic0, ic1;
    if (!tp_build_face_in(nCells, ip, plane_cells, ge_x, ge_y, ge_z, fx, sx, sy, sz, ss, sp, nSpecies, iw,
                          roY, rophi, ro, T, P, vis_lam, vis_turb, diffMethod, Sc, Sc_t, isNode, cprops, in, ic0, ic1)) return;
    TpFaceInT<double> d;
    d.n = in.n; d.iw = in.iw;
    d.rho0 = in.rho0; d.rho1 = in.rho1; d.rg0 = in.rg0; d.rg1 = in.rg1; d.f = in.f; d.geo = in.geo; d.geo_abs = in.geo_abs;
    d.ct = in.ct; d.L = in.L;
    for (int s = 0; s < nSpecies; ++s) { d.rY0[s] = in.rY0[s]; d.rY1[s] = in.rY1[s]; d.D[s] = in.D[s]; d.h[s] = in.h[s]; }
    for (int m = 0; m < TP_NQ; ++m) { d.rQ0[m] = in.rQ0[m]; d.rQ1[m] = in.rQ1[m]; }
    TpFaceOutT<double> o;
    tp_face_flux(d, o);
    for (int s = 0; s < nSpecies; ++s) add(s, o.J[s]);
    add(nSpecies, o.Jv);
    add(nSpecies + 1, o.Jl);
    for (int m = 0; m < TP_NQ; ++m) add(nSpecies + 2 + m, o.JQ[m]);
    if (brk != nullptr) {
        if (c0 < nCells) atomicAdd(&brk[(size_t)c0*3 + 1],  o.Jl);
        if (c1 < nCells) atomicAdd(&brk[(size_t)c1*3 + 1], -o.Jl);
    }
}

// mode 0: float の確定残差 (assembleResidual の結果; 比較用) を F に写す。mode 1: res_* に入れたソースだけの値を r・A に足す。
__global__ void twophase_audit_cell_d(geom_int nCells, int nSpecies, int iw, flow_float** res_roY, flow_float** res_phi,
                                      int NQ, int mode, double* r, double* A, double* F)
{
    const geom_int c = blockDim.x * blockIdx.x + threadIdx.x;
    if (c >= nCells) return;
    const size_t b = (size_t)c*NQ;
    if (mode == 0) {
        for (int s = 0; s < nSpecies; ++s) F[b + s] = (double)res_roY[s][c];
        F[b + nSpecies] = (double)res_roY[iw][c] - (double)res_phi[0][c];
        for (int m = 0; m < 1 + TP_NQ; ++m) F[b + nSpecies + 1 + m] = (double)res_phi[m][c];
    } else {
        for (int m = 0; m < 1 + TP_NQ; ++m) {
            const double S = (double)res_phi[m][c];
            r[b + nSpecies + 1 + m] += S; A[b + nSpecies + 1 + m] += fabs(S);
        }
        const double Sg = (double)res_phi[0][c];
        r[b + nSpecies] -= Sg; A[b + nSpecies] += fabs(Sg);   // 蒸気 = 総水分 − 液 (総水分にソースは無い)
    }
}
}  // namespace

static std::vector<double> g_tpAuditR0;   // 開始時の監査の max|r| (成分ごと; 未取得は −1)
const std::vector<double>& twoPhaseAuditR0() { return g_tpAuditR0; }

void twoPhaseAudit_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int iStep, bool final)
{
    if (!condResidualAuditActive(cfg) || !speciesEnabled(var) || g_roY_dev == nullptr || g_qMom0 < 0) return;
    const int oldOp = condTwoPhaseDiffusionActive(cfg) ? 0 : 1;   // 二相拡散 OFF の run は旧作用素を監査する (#1b-pre (1))
    const bool s3 = (cfg.speciesFaceReconstruction >= 2 && g_Yface_dev != nullptr && g_Pface_dev != nullptr
                     && (cfg.solver == "SLAU" || cfg.solver == "SLAU2"));
    if (!s3) { printf("[twophase-audit] step %d: skipped — the audit re-evaluates the S3 face-value advection only (speciesFaceReconstruction 2 with SLAU)\n", iStep); return; }
    if (periodicNodeActive(cfg, msh)) { printf("[twophase-audit] step %d: skipped — node periodic gather is not reproduced by the audit\n", iStep); return; }
    const int n = g_nSpecies, iw = cfg.condGasSpecies, NQ = n + 5;
    const size_t N = (size_t)msh.nCells*NQ;
    static double *r = nullptr, *A = nullptr, *F = nullptr; static size_t cap = 0;
    if (cap < N) { cudaFree(r); cudaFree(A); cudaFree(F);
        gpuErrchk( cudaMalloc((void**)&r, N*sizeof(double)) ); gpuErrchk( cudaMalloc((void**)&A, N*sizeof(double)) ); gpuErrchk( cudaMalloc((void**)&F, N*sizeof(double)) ); cap = N; }
    gpuErrchk( cudaMemset(r, 0, N*sizeof(double)) ); gpuErrchk( cudaMemset(A, 0, N*sizeof(double)) );
    const CondSpeciesProps cprops = condProps_make(cfg.condModel, cond_prop_opts(cfg));
    dim3 dimGrid_nh = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));
    twophase_audit_face_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(
        msh.nCells, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
        var.p_d["ge_x"], var.p_d["ge_y"], var.p_d["ge_z"], var.p_d["fx"], var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
        thermo_species_device_ptr_f(), n, iw, g_roY_dev, g_p_rophi_dev + g_qMom0,
        var.c_d["ro"], var.c_d["T"], var.c_d["P"], var.c_d["vis_lam"], var.c_d["vis_turb"],
        cfg.speciesDiffusionMethod, cfg.Sc, cfg.Sc_t, (cfg.discretization == "node") ? 1 : 0, cprops,
        var.p_d["massflux"], g_Yface_dev, g_Pface_dev, g_nPassive, g_qMom0, NQ, r, A,
        oldOp, (cfg.viscMethod != 0) ? 1 : 0, gasPhaseLiquid(cfg, var));
    twophase_audit_cell_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(msh.nCells, n, iw, g_resroY_dev, g_p_res_dev + g_qMom0, NQ, 0, r, A, F);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    // ソースだけを同じソースカーネルで組み直す (モーメント残差を 0 にしてから; 状態・温度は assembleResidual のまま)
    for (int m = 0; m < 1 + TP_NQ; ++m) gpuErrchk( cudaMemset(h_p_res[g_qMom0 + m], 0, (size_t)msh.nCells*sizeof(flow_float)) );
    condensationSource_d_wrapper(cfg, cuda_cfg, msh, var);
    twophase_audit_cell_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(msh.nCells, n, iw, g_resroY_dev, g_p_res_dev + g_qMom0, NQ, 1, r, A, F);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();

    std::vector<double> hr(N), hA(N), hF(N);
    gpuErrchk( cudaMemcpy(hr.data(), r, N*sizeof(double), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(hA.data(), A, N*sizeof(double), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(hF.data(), F, N*sizeof(double), cudaMemcpyDeviceToHost) );
    std::vector<flow_float> pin(msh.nCells, 0.0f);
    auto itp = var.c_d.find("scalarDirichletPin");
    if (cfg.discretization == "node" && itp != var.c_d.end() && itp->second != nullptr)
        gpuErrchk( cudaMemcpy(pin.data(), itp->second, msh.nCells*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    const double eps = 1.1920928955078125e-7;
    std::vector<double>& r0 = g_tpAuditR0;
    if (r0.size() != (size_t)NQ) r0.assign(NQ, -1.0);
    long npin = 0; for (geom_int c = 0; c < msh.nCells; ++c) if (pin[c] == 1.0f) ++npin;
    bool pass = true, nonfinite = false;
    printf("[twophase-audit] step %d %s: %s operator, double re-evaluation from the stored state (EOS via assembleResidual, face fluxes in double, source per cell; energy not covered); pinned nodes excluded %ld\n",
           iStep, final ? "(final state)" : "(initial state: r0)", oldOp ? "existing (key off: species Fick, liquid/moments advected only)" : "two-phase (key on)", npin);
    for (int q = 0; q < NQ; ++q) {
        double mr = 0.0, mA = 0.0, md = 0.0;
        for (geom_int c = 0; c < msh.nCells; ++c) {
            if (pin[c] == 1.0f) continue;
            const double rv = hr[(size_t)c*NQ + q], Av = hA[(size_t)c*NQ + q], Fv = hF[(size_t)c*NQ + q];
            if (!std::isfinite(rv) || !std::isfinite(Av) || !std::isfinite(Fv)) { nonfinite = true; continue; }
            mr = std::max(mr, fabs(rv)); mA = std::max(mA, Av);
            md = std::max(md, fabs(rv - Fv));
        }
        // 照合 (監査が同じ離散式を組んでいるか): float の確定残差との差を成分の尺度 ε·max A で。セル局所の |flux| で割ると、流束が小さく
        // 被演算子が大きい壁ノード (流束 1e-10、float の打ち消し誤差 1e-11) で 1e6 倍に見えるので使わない (2026-10-02 初回の試行)。
        const double cons = (mA > 0.0) ? md/(eps*mA) : (md == 0.0 ? 0.0 : INFINITY);
        if (!final || r0[q] < 0.0) r0[q] = mr;
        const double tol = std::max(1.0e-7*r0[q], 6.0*eps*mA);
        const double ratio = (tol > 0.0) ? mr/tol : (mr == 0.0 ? 0.0 : INFINITY);
        if (final && !(ratio <= 1.0)) pass = false;
        std::string nm = (q < n) ? ("roY" + std::to_string(q) + (q == iw ? "(w)" : "")) : (q == n ? "v" : (q == n + 1 ? "g" : (q == n + 2 ? "Q2" : (q == n + 3 ? "Q1" : "Q0"))));
        printf("[twophase-audit]   %-8s max|r| %.3e  max A %.3e  r0 %.3e  ratio %.3e (<= 1)  max|r_double - r_float|/(eps max A) %.3f\n",
               nm.c_str(), mr, mA, r0[q], ratio, cons);
    }
    if (final) printf("[twophase-audit] VERDICT: %s (component ratios max|r| / max(1e-7 r0, 6 eps max A); non-finite %s)\n",
                      (pass && !nonfinite) ? "PASS" : "NOT CONVERGED", nonfinite ? "yes" : "no");
    fflush(stdout);
}

// =============================================================================
// 再正規化の受入ゲートの計測 (#1b-pre; codex diagnose 2026-10-02 final)。更新ごとの C_q,n と F_n を device の履歴に積み、
// ログ (区間) と終了時 (末尾窓) は履歴から集計する (renormGate_d.cuh)。プロセス内の更新だけを数える (restart で 0 から)。
// =============================================================================
static double* g_rng_acc = nullptr;
static double* g_rng_hist = nullptr;
static size_t  g_rng_cap = 0, g_rng_n = 0, g_rng_logged = 0;

static flow_float** rngMomentPtrs() { return (g_qMom0 >= 0 && g_p_rophi_dev != nullptr) ? g_p_rophi_dev + g_qMom0 : nullptr; }
static double* rngBegin()
{
    if (g_rng_acc == nullptr) gpuErrchk( cudaMalloc((void**)&g_rng_acc, RNG_ACC*sizeof(double)) );
    if (g_rng_n + 1 > g_rng_cap) {
        const size_t cap = std::max<size_t>(1024, 2*g_rng_cap);
        double* h = nullptr;
        gpuErrchk( cudaMalloc((void**)&h, cap*RNG_ENTRY*sizeof(double)) );
        if (g_rng_hist != nullptr) {
            gpuErrchk( cudaMemcpy(h, g_rng_hist, g_rng_n*RNG_ENTRY*sizeof(double), cudaMemcpyDeviceToDevice) );
            cudaFree(g_rng_hist);
        }
        g_rng_hist = h; g_rng_cap = cap;
    }
    gpuErrchk( cudaMemset(g_rng_acc, 0, RNG_ACC*sizeof(double)) );
    return g_rng_acc;
}
static void rngEnd(cudaConfig& cuda_cfg)
{
    (void)cuda_cfg;
    rng_finalize_d<<<1, 1>>>(g_rng_acc, g_rng_hist + g_rng_n*RNG_ENTRY);
    gpuErrchk( cudaPeekAtLastError() );
    ++g_rng_n;
}

void renormGateLog(const solverConfig& cfg, int iStep, bool final)
{
    if (g_rng_hist == nullptr) return;
    const size_t n = g_rng_n;
    const size_t b = final ? rng_final_window_begin(n) : g_rng_logged;
    if (!final) g_rng_logged = n;
    std::vector<double> h(n*RNG_ENTRY, 0.0);
    if (n > b) gpuErrchk( cudaMemcpy(h.data() + b*RNG_ENTRY, g_rng_hist + b*RNG_ENTRY, (n - b)*RNG_ENTRY*sizeof(double), cudaMemcpyDeviceToHost) );
    const RngWindow w = rng_window_max(h, b, n);
    const double kappa = rng_kappa(cfg.nSpecies);
    char head[160];
    if (final) std::snprintf(head, sizeof(head), "[renorm-gate] FINAL step %d window: last %zu of %zu updates [%zu,%zu)", iStep, n - b, n, b, n);
    else       std::snprintf(head, sizeof(head), "[renorm-gate] step %d interval: %zu updates [%zu,%zu)", iStep + 1, n - b, b, n);
    printf("%s | max_n C_q: rhoYw %.6e rhog %.6e rhoQ2 %.6e rhoQ1 %.6e rhoQ0 %.6e | max_n max|f-1| %.6e | kappa %.7e | nonfinite %ld",
           head, w.C[0], w.C[1], w.C[2], w.C[3], w.C[4], w.F, kappa, w.nonfinite);
    if (final) printf(" | VERDICT: %s\n", rng_judge(w, kappa) ? "PASS" : "FAIL");
    else printf("\n");
    fflush(stdout);
}

// =============================================================================
// #1b-r2 診断の組立 B (condTwoPhaseDiag 3; 読むだけ)。assembleResidual の直後 (更新の前) に呼ぶ。
//   同じ格納状態・面値 (massflux・S3 面組成)・係数・ソース値から、二相系の各成分の残差を面ごとに double で組む (監査と同じカーネル)。
//   ソース値は同じソースカーネルの値を使うため、ソースが書く 14 配列 (モーメント残差 4・src_jac 4・診断 6) を退避 → 残差を 0 にして
//   ソースだけを組む → 退避した値をそのまま戻す。更新が使う組立 A (float の残差・src_jac・状態) はビット単位で変わらない。
//   この時点の B は終了時の独立残差監査 (EOS は直前の assembleResidual) と同じ定義の残差である。
// =============================================================================
static double *g_tpB_r = nullptr, *g_tpB_A = nullptr, *g_tpB_brk = nullptr;
static flow_float* g_tpB_save = nullptr;
static size_t g_tpB_cells = 0; static int g_tpB_NQ = 0;
namespace {
__global__ void tpB_source_add_d(geom_int nCells, int nSpecies, flow_float** res_phi, int NQ, double* r, double* A, double* brk)
{
    const geom_int c = blockDim.x * blockIdx.x + threadIdx.x;
    if (c >= nCells) return;
    const size_t b = (size_t)c*NQ;
    for (int m = 0; m < 1 + TP_NQ; ++m) { const double S = (double)res_phi[m][c]; r[b + nSpecies + 1 + m] += S; A[b + nSpecies + 1 + m] += fabs(S); }
    const double Sg = (double)res_phi[0][c];
    r[b + nSpecies] -= Sg; A[b + nSpecies] += fabs(Sg);
    brk[(size_t)c*3 + 2] = Sg;
}
}  // namespace
void twoPhaseDiagB_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!condTwoPhaseDiffusionActive(cfg) || cfg.condTwoPhaseDiag != 3 || g_roY_dev == nullptr || g_qMom0 < 0) return;
    if (!(cfg.speciesFaceReconstruction >= 2 && g_Yface_dev != nullptr && g_Pface_dev != nullptr && (cfg.solver == "SLAU" || cfg.solver == "SLAU2"))
        || periodicNodeActive(cfg, msh)) { g_tpB_NQ = 0; return; }   // 監査と同じ対象範囲
    const int n = g_nSpecies, iw = cfg.condGasSpecies, NQ = n + 5;
    const size_t nc = (size_t)msh.nCells;
    if (g_tpB_r == nullptr) {
        gpuErrchk( cudaMalloc((void**)&g_tpB_r, nc*NQ*sizeof(double)) ); gpuErrchk( cudaMalloc((void**)&g_tpB_A, nc*NQ*sizeof(double)) );
        gpuErrchk( cudaMalloc((void**)&g_tpB_brk, nc*3*sizeof(double)) ); gpuErrchk( cudaMalloc((void**)&g_tpB_save, nc*14*sizeof(flow_float)) );
        g_tpB_cells = nc; g_tpB_NQ = NQ;
    }
    gpuErrchk( cudaMemset(g_tpB_r, 0, nc*NQ*sizeof(double)) ); gpuErrchk( cudaMemset(g_tpB_A, 0, nc*NQ*sizeof(double)) );
    gpuErrchk( cudaMemset(g_tpB_brk, 0, nc*3*sizeof(double)) );
    const CondSpeciesProps cprops = condProps_make(cfg.condModel, cond_prop_opts(cfg));
    dim3 dimGrid_nh = dim3(ceil(msh.nNormal_halo_Planes / (flow_float)cuda_cfg.blocksize));
    twophase_audit_face_d<<<dimGrid_nh, cuda_cfg.dimBlock>>>(
        msh.nCells, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
        var.p_d["ge_x"], var.p_d["ge_y"], var.p_d["ge_z"], var.p_d["fx"], var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
        thermo_species_device_ptr_f(), n, iw, g_roY_dev, g_p_rophi_dev + g_qMom0,
        var.c_d["ro"], var.c_d["T"], var.c_d["P"], var.c_d["vis_lam"], var.c_d["vis_turb"],
        cfg.speciesDiffusionMethod, cfg.Sc, cfg.Sc_t, (cfg.discretization == "node") ? 1 : 0, cprops,
        var.p_d["massflux"], g_Yface_dev, g_Pface_dev, g_nPassive, g_qMom0, NQ, g_tpB_r, g_tpB_A,
        0, (cfg.viscMethod != 0) ? 1 : 0, gasPhaseLiquid(cfg, var), g_tpB_brk);
    gpuErrchk( cudaPeekAtLastError() );
    // ソース値: ソースが書く 14 配列を退避 → モーメント残差 0 → ソース → 足す → 戻す
    const char* names[14] = {"res_rog_0", "res_roQ0_0", "res_roQ1_0", "res_roQ2_0", "src_jac_g_0", "src_jac_Q0_0", "src_jac_Q1_0", "src_jac_Q2_0",
                             "condS_0", "condDrdt_0", "condR30_0", "condTsat_0", "condTheta_0", "condLim_0"};
    for (int k = 0; k < 14; ++k) gpuErrchk( cudaMemcpy(g_tpB_save + (size_t)k*nc, var.c_d[names[k]], nc*sizeof(flow_float), cudaMemcpyDeviceToDevice) );
    for (int k = 0; k < 4; ++k) gpuErrchk( cudaMemset(var.c_d[names[k]], 0, nc*sizeof(flow_float)) );
    condensationSource_d_wrapper(cfg, cuda_cfg, msh, var);
    tpB_source_add_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(msh.nCells, n, g_p_res_dev + g_qMom0, NQ, g_tpB_r, g_tpB_A, g_tpB_brk);
    gpuErrchk( cudaPeekAtLastError() );
    for (int k = 0; k < 14; ++k) gpuErrchk( cudaMemcpy(var.c_d[names[k]], g_tpB_save + (size_t)k*nc, nc*sizeof(flow_float), cudaMemcpyDeviceToDevice) );
    gpuErrchkKernelSync();
}
bool twoPhaseDiagBGet(const double** r, const double** A, const double** brk, int* NQ)
{
    if (g_tpB_NQ == 0 || g_tpB_r == nullptr) return false;
    *r = g_tpB_r; *A = g_tpB_A; *brk = g_tpB_brk; *NQ = g_tpB_NQ;
    return true;
}

// =============================================================================
// 診断 D1 (G1: 0 step の面作用素 A/B; plans/active/condensation-two-phase-default.md §5.1 #4、methods/condensation.md §7c)。
//   FORGE_DIAG_TP_FACES=<出力 h5> のときだけ main (runTpFacesDiag) が 1 回呼ぶ。計算経路からは呼ばない (既定 off は何もしない)。
//   呼び出し側が組立の前処理と後処理 (assembleResidualPre/Post) を 1 回ずつ通した状態から、全通常面 (normal_halo_planes; 本番の
//   species_diffusion_d / twophase_diffusion_d と同じ面の並び) について次を評価して host に返す:
//   - OFF (作用素 A): species_diffusion_d の面計算の**診断専用の写し** tpfd_off_face (本番カーネルは変更しない)。
//     写しの正しさは、写しの面流束を double で節点に組み直した値と、本番カーネル species_diffusion_d を同じ状態で
//     **退避用の 0 初期化した残差配列**に流した値 (= 拡散段の残差寄与そのもの) を節点ごとに照合して確かめる。
//   - ON (作用素 B): 本番の tp_build_face_in + tp_face_flux<float> (読むだけ) と、同じ float 入力を double に上げた
//     tp_face_flux<double> (監査 #4f と同じ経路)。中間量 (j_v⁰・補正項・乱流分) は tp_face_flux_diag (twoPhaseFaceDiag_d.cuh)。
//   符号: J は「セル 0 (ic0) へ入る向き」が正の面積分値 (本番と同じ; res[ic0] += J, res[ic1] −= J)。
//   node 境界半割面は両作用素とも拡散を足さないので skip = 1 (状態と幾何だけ記録)。
// =============================================================================
namespace {
// 出力列の配置 (列ごとに nFaces 個を並べる; host で [nFaces][width] に並べ替えて書く)
struct TpfdOutCols { int J, Jv, Jl, JQ, q, Sm, jv0, jv_corr, jv_mol, jv_turb, rgf, zv0, zv1, zupv; };
struct TpfdLayout {
    int n;
    // float 列 (F): 面の幾何・面値
    int f, sx, sy, sz, ss, pcx, pcy, pcz, dcc, delta, mu_f, mut_f, T_f, P_f;
    int x[2], y[2], z[2];
    // 両端の状態
    int ro[2], T[2], P[2], vl[2], vt[2], rg[2], wd[2], rY[2], rQ[2];
    // OFF の写し: 補正後 J、補正前 Js、分子拡散係数 D (乱流分を含まない)、面組成 Yf、乱流 D_t、ρ_f、ΣJs、q、拡散係数の組成の液 g_f
    int offJ, offJs, offD, offYf, offDt, offRof, offSumJ, offQ, offGlf;
    // ON の入力 TpFaceIn (float、全体)
    int inRho0, inRho1, inRY0, inRY1, inRg0, inRg1, inRQ0, inRQ1, inF, inGeo, inGeoAbs, inD, inCt, inH, inL;
    TpfdOutCols onF;   // F: ON float の出力
    TpfdOutCols onD;   // D: ON double 参照の出力
    int Al, Av;        // D: 誤差尺度
};
// int 列 (I) の固定番号
enum { TPFD_I_IP = 0, TPFD_I_IC0, TPFD_I_IC1, TPFD_I_SKIP, TPFD_I_UP0F, TPFD_I_UP0D, TPFD_I_SELF, TPFD_NI };

struct TpfdBlk { std::string name; int width; int kind; int off; };   // kind 0 float / 1 double / 2 int
struct TpfdBuilder {
    std::vector<TpfdBlk> blks; int nF = 0, nD = 0, nI = TPFD_NI;
    int add(const std::string& nm, int w, int kind) { int& c = (kind == 0) ? nF : nD; const int o = c; c += w; blks.push_back({nm, w, kind, o}); return o; }
    void outCols(const std::string& g, int n, int kind, TpfdOutCols& c) {
        c.J = add(g + "/J", n, kind); c.Jv = add(g + "/Jv", 1, kind); c.Jl = add(g + "/Jl", 1, kind); c.JQ = add(g + "/JQ", TP_NQ, kind);
        c.q = add(g + "/q", 1, kind); c.Sm = add(g + "/Sm", 1, kind); c.jv0 = add(g + "/jv0", 1, kind); c.jv_corr = add(g + "/jv_corr", 1, kind);
        c.jv_mol = add(g + "/jv_mol", 1, kind); c.jv_turb = add(g + "/jv_turb", 1, kind); c.rgf = add(g + "/rgf", 1, kind);
        c.zv0 = add(g + "/zv0", 1, kind); c.zv1 = add(g + "/zv1", 1, kind); c.zupv = add(g + "/zupv", 1, kind);
    }
};
static void tpfdBuildLayout(int n, TpfdLayout& L, TpfdBuilder& b)
{
    L.n = n;
    L.f = b.add("face/f", 1, 0); L.sx = b.add("face/sx", 1, 0); L.sy = b.add("face/sy", 1, 0); L.sz = b.add("face/sz", 1, 0);
    L.ss = b.add("face/ss", 1, 0); L.pcx = b.add("face/pcx", 1, 0); L.pcy = b.add("face/pcy", 1, 0); L.pcz = b.add("face/pcz", 1, 0);
    L.dcc = b.add("face/dcc", 1, 0); L.delta = b.add("face/delta", 1, 0);
    L.mu_f = b.add("face/mu_f", 1, 0); L.mut_f = b.add("face/mut_f", 1, 0); L.T_f = b.add("face/T_f", 1, 0); L.P_f = b.add("face/P_f", 1, 0);
    for (int e = 0; e < 2; ++e) {
        const std::string g = "state" + std::to_string(e);
        L.x[e] = b.add(g + "/x", 1, 0); L.y[e] = b.add(g + "/y", 1, 0); L.z[e] = b.add(g + "/z", 1, 0);
        L.ro[e] = b.add(g + "/ro", 1, 0); L.T[e] = b.add(g + "/T", 1, 0); L.P[e] = b.add(g + "/P", 1, 0);
        L.vl[e] = b.add(g + "/vis_lam", 1, 0); L.vt[e] = b.add(g + "/vis_turb", 1, 0); L.rg[e] = b.add(g + "/rog", 1, 0);
        L.wd[e] = b.add(g + "/wall_dist", 1, 0); L.rY[e] = b.add(g + "/roY", n, 0); L.rQ[e] = b.add(g + "/roQ", TP_NQ, 0);
    }
    L.offJ = b.add("off/J", n, 0); L.offJs = b.add("off/J_uncorrected", n, 0); L.offD = b.add("off/D_mol", n, 0);
    L.offYf = b.add("off/Yf", n, 0); L.offDt = b.add("off/D_t", 1, 0); L.offRof = b.add("off/ro_f", 1, 0);
    L.offSumJ = b.add("off/sumJ", 1, 0); L.offQ = b.add("off/q", 1, 0); L.offGlf = b.add("off/gl_f", 1, 0);
    L.inRho0 = b.add("on_in/rho0", 1, 0); L.inRho1 = b.add("on_in/rho1", 1, 0); L.inRY0 = b.add("on_in/rY0", n, 0); L.inRY1 = b.add("on_in/rY1", n, 0);
    L.inRg0 = b.add("on_in/rg0", 1, 0); L.inRg1 = b.add("on_in/rg1", 1, 0); L.inRQ0 = b.add("on_in/rQ0", TP_NQ, 0); L.inRQ1 = b.add("on_in/rQ1", TP_NQ, 0);
    L.inF = b.add("on_in/f", 1, 0); L.inGeo = b.add("on_in/geo", 1, 0); L.inGeoAbs = b.add("on_in/geo_abs", 1, 0);
    L.inD = b.add("on_in/D", n, 0); L.inCt = b.add("on_in/ct", 1, 0); L.inH = b.add("on_in/h", n, 0); L.inL = b.add("on_in/L", 1, 0);
    b.outCols("on_f", n, 0, L.onF);
    b.outCols("on_d", n, 1, L.onD);
    L.Al = b.add("scale/A_l", 1, 1); L.Av = b.add("scale/A_v", 1, 1);
}

// OFF の面計算 (species_diffusion_d の写し; 診断専用)。演算・評価順は species_diffusion_d (本ファイル) と同じに保つこと。
// 戻り値 false = node 境界半割面 (本番も流束を足さない)。
struct TpfdOff {
    flow_float dcc, delta, ro0, ro1, inv_ro0, inv_ro1, ro_face, T_face, P_face, mu_face, mut_face, Dt, gl_f, sumJ, q, diag_geo;
    flow_float Yf[THERMO_MAX_SPECIES], Dm[THERMO_MAX_SPECIES], Js[THERMO_MAX_SPECIES], Jc[THERMO_MAX_SPECIES], diag[THERMO_MAX_SPECIES];
};
__device__ inline bool tpfd_off_face(
    geom_int nCells, geom_int ip, geom_int ic0, geom_int ic1,
    const flow_float* ge_x, const flow_float* ge_y, const flow_float* ge_z,   // 面ごとの差 e (species_diffusion_d と同じ、段 ③)
    geom_float* fx, geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,
    const SpeciesThermoF* sp, int nSpecies, flow_float** roY,
    flow_float* ro, flow_float* T, flow_float* P, flow_float* vis_lam, flow_float* vis_turb,
    int diffMethod, flow_float Sc, flow_float Sc_t, int isNode, GasPhaseLiquid liq, TpfdOff& o)
{
    const flow_float f   = fx[ip];
    const flow_float sxx = sx[ip], syy = sy[ip], szz = sz[ip], sss = ss[ip];
    if (isNode != 0 && (ic0 >= nCells || ic1 >= nCells)) return false;

    const flow_float dccx = ge_x[ip];
    const flow_float dccy = ge_y[ip];
    const flow_float dccz = ge_z[ip];
    const flow_float dcc  = sqrtf(dccx*dccx + dccy*dccy + dccz*dccz);
    const flow_float denom = dccx*sxx + dccy*syy + dccz*szz;
    const flow_float Dsafe = (fabsf(denom) < 1.0e-30f) ? ((denom>=0.0f)?1.0e-30f:-1.0e-30f) : denom;
    const flow_float delta = dcc * sss * sss / Dsafe;

    const flow_float ro0 = max(ro[ic0], (flow_float)1.0e-30f);
    const flow_float ro1 = max(ro[ic1], (flow_float)1.0e-30f);
    const flow_float inv_ro0 = 1.0f/ro0, inv_ro1 = 1.0f/ro1;
    const flow_float g = 1.0f - f;
    const flow_float ro_face = f*ro0 + g*ro1;
    const flow_float T_face  = f*T[ic0] + g*T[ic1];
    const flow_float P_face  = f*P[ic0] + g*P[ic1];

    flow_float X[THERMO_MAX_SPECIES];
    flow_float Ys0[THERMO_MAX_SPECIES], Ys1[THERMO_MAX_SPECIES];
    flow_float ysum = 0.0f;
    for (int s=0;s<nSpecies;s++){
        Ys0[s] = roY[s][ic0]*inv_ro0;
        Ys1[s] = roY[s][ic1]*inv_ro1;
        flow_float y = f*Ys0[s] + g*Ys1[s];
        if (y<0.0f) y=0.0f; o.Yf[s]=y; ysum+=y;
    }
    const flow_float yinv = 1.0f/(ysum>1.0e-30f?ysum:1.0e-30f);
    for (int s=0;s<nSpecies;s++) o.Yf[s]*=yinv;
    const flow_float gl_f = (liq.iw >= 0 && liq.rog != nullptr)
        ? (f*liq.rog[ic0]*inv_ro0 + g*liq.rog[ic1]*inv_ro1)*yinv : 0.0f;
    species_transport_X_f(sp, nSpecies, o.Yf, liq.iw, gl_f, X);

    const flow_float mu_face  = f*vis_lam[ic0]  + g*vis_lam[ic1];
    const flow_float mut_face = f*vis_turb[ic0] + g*vis_turb[ic1];
    const flow_float Dt = (mut_face > 0.0f) ? mut_face/(ro_face*Sc_t) : 0.0f;
    const flow_float inv_dcc = 1.0f/dcc;
    const flow_float diag_geo = fabsf(delta) / max(dcc,(flow_float)1.0e-30f);

    flow_float sumJ = 0.0f;
    for (int s=0;s<nSpecies;s++){
        flow_float D;
        if (diffMethod == 1) D = thermo_Dmix_species_f(sp, nSpecies, X, s, T_face, P_face);
        else                 D = mu_face/(ro_face*Sc);
        o.Dm[s] = D;
        D += Dt;
        const flow_float roD = ro_face * D;
        o.Js[s] = roD * ((Ys1[s] - Ys0[s])*inv_dcc) * delta;
        sumJ += o.Js[s];
        o.diag[s] = roD * diag_geo;
    }
    flow_float q = 0.0f;
    for (int s=0;s<nSpecies;s++){
        const flow_float Jc = o.Js[s] - o.Yf[s]*sumJ;
        o.Jc[s] = Jc;
        const flow_float hs = thermo_h_mass_f(sp[s], T_face);
        q += hs * Jc;
    }
    o.dcc = dcc; o.delta = delta; o.ro0 = ro0; o.ro1 = ro1; o.inv_ro0 = inv_ro0; o.inv_ro1 = inv_ro1;
    o.ro_face = ro_face; o.T_face = T_face; o.P_face = P_face; o.mu_face = mu_face; o.mut_face = mut_face;
    o.Dt = Dt; o.gl_f = gl_f; o.sumJ = sumJ; o.q = q; o.diag_geo = diag_geo;
    return true;
}

// 面の記録 + OFF の写し + 写しの節点への組み直し (double)
__global__ void tpfd_off_face_d(
    geom_int nCells, geom_int nF, geom_int* normal_halo_planes, geom_int* plane_cells,
    geom_float* ccx, geom_float* ccy, geom_float* ccz,   // 記録 (L.x/L.y/L.z) だけに使う
    const flow_float* ge_x, const flow_float* ge_y, const flow_float* ge_z,   // 面ごとの差 e (OFF の写しが読む、段 ③)
    geom_float* fx, geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,
    geom_float* pcx, geom_float* pcy, geom_float* pcz,
    const SpeciesThermoF* sp, int nSpecies, flow_float** roY, flow_float** rophi,
    flow_float* ro, flow_float* T, flow_float* P, flow_float* vis_lam, flow_float* vis_turb, const flow_float* wall_dist,
    int diffMethod, flow_float Sc, flow_float Sc_t, int isNode, GasPhaseLiquid liq,
    TpfdLayout L, float* F, int* Iv,
    double* Rc, double* Ra, double* Tdc, double* Rq, double* Rqa, int* nface)   // 節点 [s*nCells + ic] (Rq/Rqa/nface は [ic])
{
    const geom_int ih = blockDim.x * blockIdx.x + threadIdx.x;
    if (ih >= nF) return;
    const size_t N = (size_t)nF;
    auto FW = [&](int c, flow_float v) { F[(size_t)c*N + ih] = v; };
    const geom_int ip  = normal_halo_planes[ih];
    const geom_int ic0 = plane_cells[2 * ip + 0];
    const geom_int ic1 = plane_cells[2 * ip + 1];
    Iv[TPFD_I_IP*N + ih] = ip; Iv[TPFD_I_IC0*N + ih] = ic0; Iv[TPFD_I_IC1*N + ih] = ic1;
    FW(L.f, fx[ip]); FW(L.sx, sx[ip]); FW(L.sy, sy[ip]); FW(L.sz, sz[ip]); FW(L.ss, ss[ip]);
    FW(L.pcx, pcx[ip]); FW(L.pcy, pcy[ip]); FW(L.pcz, pcz[ip]);
    const geom_int cc[2] = {ic0, ic1};
    for (int e = 0; e < 2; ++e) {
        const geom_int c = cc[e];
        FW(L.x[e], ccx[c]); FW(L.y[e], ccy[c]); FW(L.z[e], ccz[c]);
        FW(L.ro[e], ro[c]); FW(L.T[e], T[c]); FW(L.P[e], P[c]); FW(L.vl[e], vis_lam[c]); FW(L.vt[e], vis_turb[c]);
        FW(L.rg[e], rophi[0][c]);
        if (wall_dist != nullptr) FW(L.wd[e], wall_dist[c]);
        for (int s = 0; s < nSpecies; ++s) FW(L.rY[e] + s, roY[s][c]);
        for (int m = 0; m < TP_NQ; ++m) FW(L.rQ[e] + m, rophi[1+m][c]);
    }
    TpfdOff o;
    if (!tpfd_off_face(nCells, ip, ic0, ic1, ge_x, ge_y, ge_z, fx, sx, sy, sz, ss, sp, nSpecies, roY, ro, T, P, vis_lam, vis_turb,
                       diffMethod, Sc, Sc_t, isNode, liq, o)) {
        Iv[TPFD_I_SKIP*N + ih] = 1;   // node 境界半割面 (拡散は足さない)
        return;
    }
    Iv[TPFD_I_SKIP*N + ih] = 0;
    FW(L.dcc, o.dcc); FW(L.delta, o.delta); FW(L.mu_f, o.mu_face); FW(L.mut_f, o.mut_face); FW(L.T_f, o.T_face); FW(L.P_f, o.P_face);
    for (int s = 0; s < nSpecies; ++s) {
        FW(L.offJ + s, o.Jc[s]); FW(L.offJs + s, o.Js[s]); FW(L.offD + s, o.Dm[s]); FW(L.offYf + s, o.Yf[s]);
    }
    FW(L.offDt, o.Dt); FW(L.offRof, o.ro_face); FW(L.offSumJ, o.sumJ); FW(L.offQ, o.q); FW(L.offGlf, o.gl_f);
    // 写しの節点への組み直し (本番と同じ加算先・符号; 和は double)。対角は本番と同じ float 積 diag·(1/ρ) を足す。
    const size_t nc = (size_t)nCells;
    if (ic0 < nCells) {
        for (int s = 0; s < nSpecies; ++s) {
            atomicAdd(&Rc[s*nc + ic0], (double)o.Jc[s]); atomicAdd(&Ra[s*nc + ic0], fabs((double)o.Jc[s]));
            atomicAdd(&Tdc[s*nc + ic0], (double)(o.diag[s]*o.inv_ro0));
        }
        atomicAdd(&Rq[ic0], (double)o.q); atomicAdd(&Rqa[ic0], fabs((double)o.q)); atomicAdd(&nface[ic0], 1);
    }
    if (ic1 < nCells) {
        for (int s = 0; s < nSpecies; ++s) {
            atomicAdd(&Rc[s*nc + ic1], -(double)o.Jc[s]); atomicAdd(&Ra[s*nc + ic1], fabs((double)o.Jc[s]));
            atomicAdd(&Tdc[s*nc + ic1], (double)(o.diag[s]*o.inv_ro1));
        }
        atomicAdd(&Rq[ic1], -(double)o.q); atomicAdd(&Rqa[ic1], fabs((double)o.q)); atomicAdd(&nface[ic1], 1);
    }
}

__device__ inline void tpfd_write_out_f(float* F, size_t N, geom_int ih, const TpfdOutCols& c, int n, const TpFaceOut& o, const TpFaceDiagT<float>& d)
{
    for (int s = 0; s < n; ++s) F[(size_t)(c.J + s)*N + ih] = o.J[s];
    F[(size_t)c.Jv*N + ih] = o.Jv; F[(size_t)c.Jl*N + ih] = o.Jl;
    for (int m = 0; m < TP_NQ; ++m) F[(size_t)(c.JQ + m)*N + ih] = o.JQ[m];
    F[(size_t)c.q*N + ih] = o.q; F[(size_t)c.Sm*N + ih] = o.Sm;
    F[(size_t)c.jv0*N + ih] = d.jv0; F[(size_t)c.jv_corr*N + ih] = d.jv_corr; F[(size_t)c.jv_mol*N + ih] = d.jv_mol; F[(size_t)c.jv_turb*N + ih] = d.jv_turb;
    F[(size_t)c.rgf*N + ih] = d.rgf; F[(size_t)c.zv0*N + ih] = d.zv0; F[(size_t)c.zv1*N + ih] = d.zv1; F[(size_t)c.zupv*N + ih] = d.zupv;
}
__device__ inline void tpfd_write_out_d(double* D, size_t N, geom_int ih, const TpfdOutCols& c, int n, const TpFaceOutT<double>& o, const TpFaceDiagT<double>& d)
{
    for (int s = 0; s < n; ++s) D[(size_t)(c.J + s)*N + ih] = o.J[s];
    D[(size_t)c.Jv*N + ih] = o.Jv; D[(size_t)c.Jl*N + ih] = o.Jl;
    for (int m = 0; m < TP_NQ; ++m) D[(size_t)(c.JQ + m)*N + ih] = o.JQ[m];
    D[(size_t)c.q*N + ih] = o.q; D[(size_t)c.Sm*N + ih] = o.Sm;
    D[(size_t)c.jv0*N + ih] = d.jv0; D[(size_t)c.jv_corr*N + ih] = d.jv_corr; D[(size_t)c.jv_mol*N + ih] = d.jv_mol; D[(size_t)c.jv_turb*N + ih] = d.jv_turb;
    D[(size_t)c.rgf*N + ih] = d.rgf; D[(size_t)c.zv0*N + ih] = d.zv0; D[(size_t)c.zv1*N + ih] = d.zv1; D[(size_t)c.zupv*N + ih] = d.zupv;
}

// ON: 本番の tp_build_face_in + tp_face_flux<float> (読むだけ) と double 参照
__global__ void tpfd_on_face_d(
    geom_int nCells, geom_int nF, geom_int* normal_halo_planes, geom_int* plane_cells,
    const flow_float* ge_x, const flow_float* ge_y, const flow_float* ge_z,   // 面ごとの差 e (本番と同じ、段 ③)
    geom_float* fx, geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,
    const SpeciesThermoF* sp, int nSpecies, int iw, flow_float** roY, flow_float** rophi,
    flow_float* ro, flow_float* T, flow_float* P, flow_float* vis_lam, flow_float* vis_turb,
    int diffMethod, flow_float Sc, flow_float Sc_t, int isNode, CondSpeciesProps cprops,
    TpfdLayout L, float* F, double* D, int* Iv)
{
    const geom_int ih = blockDim.x * blockIdx.x + threadIdx.x;
    if (ih >= nF) return;
    const size_t N = (size_t)nF;
    TpFaceIn in; geom_int a0, a1;
    if (!tp_build_face_in(nCells, normal_halo_planes[ih], plane_cells, ge_x, ge_y, ge_z, fx, sx, sy, sz, ss, sp, nSpecies, iw,
                          roY, rophi, ro, T, P, vis_lam, vis_turb, diffMethod, Sc, Sc_t, isNode, cprops, in, a0, a1)) {
        if (Iv[TPFD_I_SKIP*N + ih] == 0) Iv[TPFD_I_SKIP*N + ih] = 2;   // OFF は評価したが ON は skip (起こらないはず; 記録)
        return;
    }
    if (Iv[TPFD_I_SKIP*N + ih] != 0) Iv[TPFD_I_SKIP*N + ih] = 3;       // ON は評価したが OFF は skip (同上)
    auto FW = [&](int c, flow_float v) { F[(size_t)c*N + ih] = v; };
    FW(L.inRho0, in.rho0); FW(L.inRho1, in.rho1); FW(L.inRg0, in.rg0); FW(L.inRg1, in.rg1);
    for (int s = 0; s < nSpecies; ++s) { FW(L.inRY0 + s, in.rY0[s]); FW(L.inRY1 + s, in.rY1[s]); FW(L.inD + s, in.D[s]); FW(L.inH + s, in.h[s]); }
    for (int m = 0; m < TP_NQ; ++m) { FW(L.inRQ0 + m, in.rQ0[m]); FW(L.inRQ1 + m, in.rQ1[m]); }
    FW(L.inF, in.f); FW(L.inGeo, in.geo); FW(L.inGeoAbs, in.geo_abs); FW(L.inCt, in.ct); FW(L.inL, in.L);
    // float (本番と同じ関数)
    TpFaceOut o; tp_face_flux(in, o);
    TpFaceDiagT<float> d; tp_face_flux_diag(in, d);
    tpfd_write_out_f(F, N, ih, L.onF, nSpecies, o, d);
    Iv[TPFD_I_UP0F*N + ih] = o.up0;
    Iv[TPFD_I_SELF*N + ih] = (o.Sm == d.Sm && o.up0 == d.up0) ? 0 : 1;   // 中間量の写しが本番と同じ Σj⁰・風上を作ったか
    // double 参照 (監査 #4f と同じ昇格)
    TpFaceInT<double> dd; tp_face_in_to_double(in, dd);
    TpFaceOutT<double> od; tp_face_flux(dd, od);
    TpFaceDiagT<double> d2; tp_face_flux_diag(dd, d2);
    tpfd_write_out_d(D, N, ih, L.onD, nSpecies, od, d2);
    Iv[TPFD_I_UP0D*N + ih] = od.up0;
    double Al, Av; tp_face_err_scales(dd, Al, Av);
    D[(size_t)L.Al*N + ih] = Al; D[(size_t)L.Av*N + ih] = Av;
}
}  // namespace

bool twoPhaseFaceDiag_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, TpFaceDiagData& out, std::string& why)
{
    const int n = g_nSpecies, iw = cfg.condGasSpecies;
    const GasPhaseLiquid liq = gasPhaseLiquid(cfg, var);
    if (!speciesEnabled(var) || g_roY_dev == nullptr) { why = "species transport is not active (nSpecies < 2)"; return false; }
    if (g_qMom0 < 0 || g_p_rophi_dev == nullptr) { why = "condensation moments are not registered (passive-scalar moment arrays missing)"; return false; }
    if (liq.iw < 0 || iw < 0) { why = "not a TP-carrier condensation run (condGasSpecies < 0, thermalMethod != 2 or nSpecies < 2)"; return false; }
    if (cfg.viscMethod == 0) { why = "viscMethod 0: neither operator diffuses"; return false; }
    if (n > THERMO_MAX_SPECIES) { why = "nSpecies exceeds THERMO_MAX_SPECIES"; return false; }

    const geom_int nF = msh.nNormal_halo_Planes;
    const size_t N = (size_t)nF, nc = (size_t)msh.nCells, nca = (size_t)msh.nCells_all;
    const int isNode = (cfg.discretization == "node") ? 1 : 0;
    TpfdLayout L; TpfdBuilder b; tpfdBuildLayout(n, L, b);

    float* F = nullptr; double* D = nullptr; int* Iv = nullptr;
    double *Rc = nullptr, *Ra = nullptr, *Tdc = nullptr, *Rq = nullptr, *Rqa = nullptr; int* nface = nullptr;
    gpuErrchk( cudaMalloc((void**)&F, (size_t)b.nF*N*sizeof(float)) );
    gpuErrchk( cudaMalloc((void**)&D, (size_t)b.nD*N*sizeof(double)) );
    gpuErrchk( cudaMalloc((void**)&Iv, (size_t)b.nI*N*sizeof(int)) );
    gpuErrchk( cudaMemset(F, 0xFF, (size_t)b.nF*N*sizeof(float)) );     // 未評価は NaN
    gpuErrchk( cudaMemset(D, 0xFF, (size_t)b.nD*N*sizeof(double)) );
    gpuErrchk( cudaMemset(Iv, 0, (size_t)b.nI*N*sizeof(int)) );
    gpuErrchk( cudaMalloc((void**)&Rc, n*nc*sizeof(double)) ); gpuErrchk( cudaMalloc((void**)&Ra, n*nc*sizeof(double)) );
    gpuErrchk( cudaMalloc((void**)&Tdc, n*nc*sizeof(double)) );
    gpuErrchk( cudaMalloc((void**)&Rq, nc*sizeof(double)) ); gpuErrchk( cudaMalloc((void**)&Rqa, nc*sizeof(double)) );
    gpuErrchk( cudaMalloc((void**)&nface, nc*sizeof(int)) );
    gpuErrchk( cudaMemset(Rc, 0, n*nc*sizeof(double)) ); gpuErrchk( cudaMemset(Ra, 0, n*nc*sizeof(double)) );
    gpuErrchk( cudaMemset(Tdc, 0, n*nc*sizeof(double)) );
    gpuErrchk( cudaMemset(Rq, 0, nc*sizeof(double)) ); gpuErrchk( cudaMemset(Rqa, 0, nc*sizeof(double)) );
    gpuErrchk( cudaMemset(nface, 0, nc*sizeof(int)) );

    auto itw = var.c_d.find("wall_dist");
    const flow_float* wd = (itw != var.c_d.end()) ? itw->second : nullptr;
    dim3 grid = dim3(ceil(nF / (flow_float)cuda_cfg.blocksize));
    tpfd_off_face_d<<<grid, cuda_cfg.dimBlock>>>(
        msh.nCells, nF, msh.normal_halo_planes_d, msh.map_plane_cells_d,
        var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["ge_x"], var.p_d["ge_y"], var.p_d["ge_z"],
        var.p_d["fx"], var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
        var.p_d["pcx"], var.p_d["pcy"], var.p_d["pcz"],
        thermo_species_device_ptr_f(), n, g_roY_dev, g_p_rophi_dev + g_qMom0,
        var.c_d["ro"], var.c_d["T"], var.c_d["P"], var.c_d["vis_lam"], var.c_d["vis_turb"], wd,
        cfg.speciesDiffusionMethod, cfg.Sc, cfg.Sc_t, isNode, liq, L, F, Iv, Rc, Ra, Tdc, Rq, Rqa, nface);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();
    const CondSpeciesProps cprops = condProps_make(cfg.condModel, cond_prop_opts(cfg));
    tpfd_on_face_d<<<grid, cuda_cfg.dimBlock>>>(
        msh.nCells, nF, msh.normal_halo_planes_d, msh.map_plane_cells_d,
        var.p_d["ge_x"], var.p_d["ge_y"], var.p_d["ge_z"],
        var.p_d["fx"], var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
        thermo_species_device_ptr_f(), n, iw, g_roY_dev, g_p_rophi_dev + g_qMom0,
        var.c_d["ro"], var.c_d["T"], var.c_d["P"], var.c_d["vis_lam"], var.c_d["vis_turb"],
        cfg.speciesDiffusionMethod, cfg.Sc, cfg.Sc_t, isNode, cprops, L, F, D, Iv);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();

    // 本番カーネル species_diffusion_d を同じ状態・同じ引数で、0 初期化した退避用の残差・対角・エネルギー配列に流す
    // (= 拡散段の残差寄与そのもの。本番の res_* には触れない)。
    std::vector<flow_float*> hRes(n), hTd(n);
    for (int s = 0; s < n; ++s) {
        gpuErrchk( cudaMalloc((void**)&hRes[s], nca*sizeof(flow_float)) ); gpuErrchk( cudaMemset(hRes[s], 0, nca*sizeof(flow_float)) );
        gpuErrchk( cudaMalloc((void**)&hTd[s], nca*sizeof(flow_float)) );  gpuErrchk( cudaMemset(hTd[s], 0, nca*sizeof(flow_float)) );
    }
    flow_float* roeP = nullptr;
    gpuErrchk( cudaMalloc((void**)&roeP, nca*sizeof(flow_float)) ); gpuErrchk( cudaMemset(roeP, 0, nca*sizeof(flow_float)) );
    flow_float** dRes = uploadPtrs(hRes);
    flow_float** dTd  = uploadPtrs(hTd);
    species_diffusion_d<<<grid, cuda_cfg.dimBlock>>>(
        msh.nCells, msh.nNormal_halo_Planes, msh.normal_halo_planes_d, msh.map_plane_cells_d,
        var.p_d["ge_x"], var.p_d["ge_y"], var.p_d["ge_z"],
        var.p_d["fx"], var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
        thermo_species_device_ptr_f(), g_nSpecies,
        g_roY_dev, dRes, dTd,
        var.c_d["ro"], var.c_d["T"], var.c_d["P"], var.c_d["vis_lam"], var.c_d["vis_turb"],
        roeP,
        cfg.speciesDiffusionMethod, cfg.Sc, cfg.Sc_t,
        isNode, g_dYdx_dev, g_dYdy_dev, g_dYdz_dev,
        liq);
    gpuErrchk( cudaPeekAtLastError() ); gpuErrchkKernelSync();

    // host へ
    out = TpFaceDiagData();
    out.nFaces = nF; out.nNodes = msh.nCells; out.nSpecies = n; out.iw = iw;
    std::vector<float> hF((size_t)b.nF*N); std::vector<double> hD((size_t)b.nD*N); std::vector<int> hI((size_t)b.nI*N);
    gpuErrchk( cudaMemcpy(hF.data(), F, hF.size()*sizeof(float), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(hD.data(), D, hD.size()*sizeof(double), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(hI.data(), Iv, hI.size()*sizeof(int), cudaMemcpyDeviceToHost) );
    // 列 → [nFaces][width] (行優先)
    for (const auto& bk : b.blks) {
        TpFaceDiagData::Block B; B.name = bk.name; B.width = bk.width; B.kind = bk.kind;
        if (bk.kind == 0) { B.f.resize(N*bk.width); for (int k = 0; k < bk.width; ++k) for (size_t i = 0; i < N; ++i) B.f[i*bk.width + k] = hF[(size_t)(bk.off + k)*N + i]; }
        else              { B.d.resize(N*bk.width); for (int k = 0; k < bk.width; ++k) for (size_t i = 0; i < N; ++i) B.d[i*bk.width + k] = hD[(size_t)(bk.off + k)*N + i]; }
        out.face.push_back(std::move(B));
    }
    const char* iname[TPFD_NI] = {"face/ip", "face/ic0", "face/ic1", "face/skip", "on_f/up0", "on_d/up0", "on_f/selfcheck_mismatch"};
    for (int c = 0; c < TPFD_NI; ++c) {
        TpFaceDiagData::Block B; B.name = iname[c]; B.width = 1; B.kind = 2;
        B.i.assign(hI.begin() + (size_t)c*N, hI.begin() + (size_t)(c + 1)*N);
        out.face.push_back(std::move(B));
    }

    // OFF の写しの照合 (節点): 写しの double 組み直し vs 本番カーネルの float 寄与
    std::vector<double> hRc(n*nc), hRa(n*nc), hTdc(n*nc), hRq(nc), hRqa(nc); std::vector<int> hNf(nc);
    gpuErrchk( cudaMemcpy(hRc.data(), Rc, hRc.size()*sizeof(double), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(hRa.data(), Ra, hRa.size()*sizeof(double), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(hTdc.data(), Tdc, hTdc.size()*sizeof(double), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(hRq.data(), Rq, nc*sizeof(double), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(hRqa.data(), Rqa, nc*sizeof(double), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(hNf.data(), nface, nc*sizeof(int), cudaMemcpyDeviceToHost) );
    std::vector<float> hRp(n*nc), hTdp(n*nc), hRoeP(nc);
    for (int s = 0; s < n; ++s) {
        gpuErrchk( cudaMemcpy(hRp.data() + s*nc, hRes[s], nc*sizeof(float), cudaMemcpyDeviceToHost) );
        gpuErrchk( cudaMemcpy(hTdp.data() + s*nc, hTd[s], nc*sizeof(float), cudaMemcpyDeviceToHost) );
    }
    gpuErrchk( cudaMemcpy(hRoeP.data(), roeP, nc*sizeof(float), cudaMemcpyDeviceToHost) );
    // 節点ブロック ([nNodes][width])
    auto nodeBlkD = [&](const std::string& nm, const std::vector<double>& v, int w) {
        TpFaceDiagData::Block B; B.name = nm; B.width = w; B.kind = 1; B.d.resize(nc*w);
        for (int k = 0; k < w; ++k) for (size_t i = 0; i < nc; ++i) B.d[i*w + k] = v[(size_t)k*nc + i];
        out.node.push_back(std::move(B)); };
    auto nodeBlkF = [&](const std::string& nm, const std::vector<float>& v, int w) {
        TpFaceDiagData::Block B; B.name = nm; B.width = w; B.kind = 0; B.f.resize(nc*w);
        for (int k = 0; k < w; ++k) for (size_t i = 0; i < nc; ++i) B.f[i*w + k] = v[(size_t)k*nc + i];
        out.node.push_back(std::move(B)); };
    nodeBlkF("off_check/R_prod", hRp, n);   nodeBlkD("off_check/R_copy", hRc, n);   nodeBlkD("off_check/R_absum", hRa, n);
    nodeBlkF("off_check/diag_prod", hTdp, n); nodeBlkD("off_check/diag_copy", hTdc, n);
    nodeBlkF("off_check/roe_prod", hRoeP, 1); nodeBlkD("off_check/roe_copy", hRq, 1); nodeBlkD("off_check/roe_absum", hRqa, 1);
    { TpFaceDiagData::Block B; B.name = "off_check/nface"; B.width = 1; B.kind = 2; B.i = hNf; out.node.push_back(std::move(B)); }

    // 要約 (ログ): 成分ごとに max|copy − prod|、max|prod|、相対差、float 和の厳密な丸め上界 nface·ε₃₂·Σ|J| に対する比の最大
    const double eps = 1.1920928955078125e-7;
    char line[512];
    auto summarize = [&](const std::string& nm, const double* cp, const float* pr, const double* ab) {
        double md = 0.0, mp = 0.0, mr = 0.0; long nbad = 0;
        for (size_t i = 0; i < nc; ++i) {
            const double d = fabs(cp[i] - (double)pr[i]);
            md = std::max(md, d); mp = std::max(mp, fabs((double)pr[i]));
            const double bound = (double)std::max(hNf[i], 1)*eps*(ab ? ab[i] : 0.0);
            const double r = (bound > 0.0) ? d/bound : (d == 0.0 ? 0.0 : INFINITY);
            mr = std::max(mr, r); if (!(r <= 1.0)) ++nbad;
        }
        snprintf(line, sizeof(line), "[tp-faces] OFF copy vs production species_diffusion_d  %-12s max|copy-prod| %.3e  max|prod| %.3e  rel %.3e  max ratio to nface*eps32*sum|J| %.3f  nodes over bound %ld",
                 nm.c_str(), md, mp, (mp > 0.0) ? md/mp : 0.0, mr, nbad);
        out.summary.push_back(line);
    };
    for (int s = 0; s < n; ++s) {
        const std::string nm = "roY" + std::to_string(s) + (s == iw ? "(w)" : "");
        summarize(nm, hRc.data() + s*nc, hRp.data() + s*nc, hRa.data() + s*nc);
    }
    summarize("roe", hRq.data(), hRoeP.data(), hRqa.data());
    for (int s = 0; s < n; ++s) {   // 対角は全項が正なので Σ|項| = 値そのもの
        const std::string nm = "diag" + std::to_string(s);
        summarize(nm, hTdc.data() + s*nc, hTdp.data() + s*nc, hTdc.data() + s*nc);
    }
    long nEval = 0, nSkip = 0, nOdd = 0, nSelf = 0;
    for (size_t i = 0; i < N; ++i) {
        const int sk = hI[(size_t)TPFD_I_SKIP*N + i];
        if (sk == 0) ++nEval; else if (sk == 1) ++nSkip; else ++nOdd;
        if (sk == 0 && hI[(size_t)TPFD_I_SELF*N + i] != 0) ++nSelf;
    }
    snprintf(line, sizeof(line), "[tp-faces] faces %ld: evaluated %ld, node boundary half-faces skipped %ld, OFF/ON skip disagreement %ld, intermediate self-check mismatches (Sm/up0) %ld",
             (long)N, nEval, nSkip, nOdd, nSelf);
    out.summary.push_back(line);

    for (int s = 0; s < n; ++s) { cudaFree(hRes[s]); cudaFree(hTd[s]); }
    cudaFree(roeP); cudaFree(dRes); cudaFree(dTd);
    cudaFree(F); cudaFree(D); cudaFree(Iv); cudaFree(Rc); cudaFree(Ra); cudaFree(Tdc); cudaFree(Rq); cudaFree(Rqa); cudaFree(nface);
    return true;
}
