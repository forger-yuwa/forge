#pragma once
// =============================================================================
// 対流フラックス スキーム共通の device ヘルパと診断グローバル。
//   - convectiveFlux_d.cu からのみ include される (単一 TU)。__device__/__constant__
//     の定義をここに置くため、他 TU から include しないこと (多重定義になる)。
//   - sign_sano は betaPls_slau/betaMns_slau より前に置き正順で参照させる。
// =============================================================================
#include "cuda_forge/reconIncrement_d.cuh"


// free-stream 保存: 対流流束の圧力項を (p_tilde - d_pRef)*s で組むための基準静圧。
// wrapper で cfg.pRef を cudaMemcpyToSymbol。既定 0.0 で従来挙動 (ビット不変)。
// 非直交メッシュで大きな p*s を float32 加算する際の桁落ち(metric closure 由来の
// 偽運動量源)を抑える。詳細: plans/active/convection-freestream-preserving-flux.md
__constant__ flow_float d_pRef;

// U∞≠0 の動く一様流保存 (KEEP CPG): 移流項を参照一様流 (d_roRef, d_uRef*) の流束 F∞(s) との
// 差分形因数分解で組むための参照状態。wrapper で cfg.roRef/uRef* を cudaMemcpyToSymbol。
// 既定 0.0 = off (ビット不変)。詳細: plans/active/convection-freestream-preserving-flux.md §8
__constant__ flow_float d_roRef;
__constant__ flow_float d_uRefX;
__constant__ flow_float d_uRefY;
__constant__ flow_float d_uRefZ;

// D2a/D4 診断 (env ゲート, 既定 off = ビット不変): 多成分 TP の contact 混合層 limit-cycle 切り分け。
//   FORGE_CONTACT_1ST=1 : 面組成センサ s_Y=max_s|Y_sR-Y_sL|/(Y_sR+Y_sL+ε) が g_contactThresh を超える
//                         内部面で flow(ρ,p,u) の MUSCL 再構成を無効化し1次 (セル値) にする (D2a)。species は不変。
//   FORGE_CONTACT_LOG=1 : s_Y>g_contactLogThresh の面で L/R 状態を printf (mixed-order/limiter chatter 観察, D4)。
__device__ int        g_contact1st = 0;
__device__ flow_float g_contactThresh = 0.05f;
__device__ int        g_contactLog = 0;
__device__ flow_float g_contactLogThresh = 0.3f;
// 連続ブレンド版 (chatter-free): ω_Y=min(1, s_Y/g_contactBlend) で flow 再構成を 2次→1次 へ滑らかに寄せる。
// 0=off。hard 切替 (FORGE_CONTACT_1ST) と違いセンサ閾値の bang-bang chatter を起こさない (D2a 清浄版/製品候補)。
__device__ flow_float g_contactBlend = 0.0f;
// S2 診断 (FORGE_FACE_THERMOY=1, 既定0): TP face thermo の組成 mixed-order を解消。
// R_mix/T/h を **face 補間組成** Y_face=f·Y_ic0+(1-f)·Y_ic1 (正規化) で評価し、ρ_L,P_L (2次再構成) と
// 同じ face 位置の組成に整合させる (現行は owner セル組成=1次 → ΔT_f^MO~30K)。species 流束は不変。
__device__ int g_faceThermoY = 0;
// node-centered: 2 次再構成の目標点を双対面重心 (pcx) でなく **エッジ中点 ½(x_A+x_B)** にする (SU2 流)。
// 値位置=ノードでは双対面重心が伸縮メッシュでエッジ中点から法線方向にずれ (Δ=(q−1)δ/4)、その分の
// 壁法線勾配 (∂φ/∂r Δ) が接線面の面値に混入する。高 AR (y+1) で不安定 (case/43 run_0042 等)。
// 1 で dc0p = ½(cc1−cc0), dc1p = −½(cc1−cc0)。cell/既定 0 でビット不変。内部面のみ (境界半割面は従来)。
__device__ int g_reconEdgeMid = 0;
// speciesFaceReconstruction==1: Y_s を ρ/Y 勾配 + min(ψ_ρ,ψ_Y) で face へ再構成し thermo/species 流束で
// 同一 face 組成を使う (proper S2/S3)。wrapper で cfg.speciesFaceReconstruction を設定。
__device__ int g_speciesFaceRecon = 0;
// 再構成 Y が [0,1] を外れた face 数のカウンタ (診断; clamp 前に atomicAdd)。
__device__ unsigned long long g_speciesOvershoot = 0;
// multispeciesRhoYCommonLimiter==1 (opt-in 診断): ρ と全 species に共通リミタ
//   ψ_ρY = min(ψ_ρ, min_s ψ_Y_s) を適用し、ρ_f=ρ(Y_f) の熱力学整合だけを切り分ける。
//   p・速度は従来どおり各自のリミタ。default 0 でビット不変。
__device__ int g_rhoYCommonLim = 0;
// rho-Y 共通リミタ診断カウンタ (すべて g_rhoYCommonLim 時のみ更新)。
__device__ int g_psiRhoY_min_scaled = 2000000;          // min(ψ_ρY)·1e6 (atomicMin)
__device__ unsigned long long g_psiRhoY_lt001 = 0;      // ψ_ρY<0.01 の face-side 数
__device__ unsigned long long g_psiRhoY_lt01  = 0;      // ψ_ρY<0.1  の face-side 数
__device__ unsigned long long g_rhoYMinByRho     = 0;   // min を ρ が決めた face-side 数
__device__ unsigned long long g_rhoYMinBySpecies = 0;   // min を species が決めた face-side 数
__device__ unsigned long long g_rhoYFallback = 0;       // 非実現可能で cell 値へ fallback した face-side 数
// W2 (plan convection-node-wall-reconstruction §4.23 / §6.4 V1): 面単位の非物理な再構成の**発火計測**。
// 既定 0 で atomicAdd は一切走らない。フォールバック本体は V1 で発火の実在を確かめてから入れる。
__device__ int   g_badReconDiag  = 0;      // 1 で計測 ON (space.badReconDiag)
__device__ flow_float g_badReconRoMin = 0.0f;   // ρ の床 (physProp.roMin)
__device__ flow_float g_badReconPMin  = 0.0f;   // P の床 (physProp.pMin)
__device__ unsigned long long g_badReconFaces = 0;   // ρ か P が床以下になった面の数
__device__ unsigned long long g_badReconRo    = 0;   // うち ρ が床以下
__device__ unsigned long long g_badReconP     = 0;   // うち P が床以下
__device__ unsigned long long g_badReconTotal = 0;   // 判定した面の総数 (分母)
// フォールバック本体 (§4.23): 発火した面を N 回の訪問だけ 1 次に落とす。SU2 `UpdateNonPhysicalEdgeCounter` と同一。
__device__ int g_badReconHyst = 0;                   // 0 = OFF、N = 1 次に落とす訪問回数 (SU2 は 20)
__device__ signed char* g_badReconCnt = nullptr;     // 面ごとのカウンタ (size = nPlanes)
__device__ unsigned long long g_badReconActive = 0;  // カウンタが生きていて 1 次化した面の数

__device__ int g_Yface_min_scaled =  2000000;           // min(Y_face)·1e6 (atomicMin)
__device__ int g_Yface_max_scaled = -2000000;           // max(Y_face)·1e6 (atomicMax)

// K7: pow(x,2.0) は exp(2*log(x)) に展開され重い。2乗は乗算 1 命令で済むため sq() に置換。
__device__ __forceinline__ flow_float sq(flow_float x) { return x * x; }

__device__ flow_float interp_MUSCL_2nd(int scheme, int limit_scheme,
                                       flow_float phiC, flow_float phiD, 
                                       flow_float dphidx, flow_float dphidy, flow_float dphidz,
                                       flow_float dphidxD, flow_float dphidyD, flow_float dphidzD,
                                       flow_float dx , flow_float dy , flow_float dz,
                                       flow_float cpdx, flow_float cpdy, flow_float cpdz,
                                       flow_float f, flow_float limiter
                                      )
{
    flow_float phif;
    flow_float k;
    flow_float r;
    flow_float psi_r;

    // 増分は recon_increment (reconIncrement_d.cuh) が唯一の定義。リミッタも同じ関数を呼ぶ。
    phif = phiC + limiter*recon_increment(1, phiC, phiD, dphidx, dphidy, dphidz, cpdx, cpdy, cpdz);

    return phif;
};

__device__ flow_float interp_MUSCL_3rd(int scheme, int limit_scheme,
                                       flow_float phiC, flow_float phiD, 
                                       flow_float dphidx, flow_float dphidy, flow_float dphidz,
                                       flow_float dphidxD, flow_float dphidyD, flow_float dphidzD,
                                       flow_float dx , flow_float dy , flow_float dz,
                                       flow_float cpdx, flow_float cpdy, flow_float cpdz,
                                       flow_float f, flow_float limiter
                                      )
{
    flow_float phif;
    flow_float k;
    flow_float r;
    flow_float psi_r;

    (void)k;
    phif = phiC + limiter*recon_increment(2, phiC, phiD, dphidx, dphidy, dphidz, cpdx, cpdy, cpdz);

    return phif;
};


__device__ flow_float interp_1stUp(int scheme, int limit_scheme,
                                 flow_float phiC, flow_float phiD, 
                                 flow_float dphidx, flow_float dphidy, flow_float dphidz,
                                 flow_float dphidxD, flow_float dphidyD, flow_float dphidzD,
                                 flow_float dx , flow_float dy , flow_float dz,
                                 flow_float cpdx, flow_float cpdy, flow_float cpdz,
                                 flow_float f, flow_float limiter
                                )
{
    return phiC;
};

__device__ flow_float interp_MINMOD(int scheme, int limit_scheme,
                                    flow_float phiC, flow_float phiD, 
                                    flow_float dphidxC, flow_float dphidyC, flow_float dphidzC,
                                    flow_float dphidxD, flow_float dphidyD, flow_float dphidzD,
                                    flow_float dx , flow_float dy , flow_float dz,
                                    flow_float cpdx, flow_float cpdy, flow_float cpdz,
                                    flow_float f, flow_float limiter
                                   )
{
    flow_float phif;
    flow_float k;
    flow_float r;
    flow_float psi_r;
    flow_float DD2dx;
    flow_float DD2dy;
    flow_float DD2dz;
    flow_float phiDD;
    flow_float phiU;
    flow_float limit;

    DD2dx = (2.0f*f-1.0f)*dx;
    DD2dy = (2.0f*f-1.0f)*dy;
    DD2dz = (2.0f*f-1.0f)*dz;

    phiDD = phiD - (DD2dx*dphidxD + DD2dy*dphidyD + DD2dz*dphidzD );
    phiU  = phiDD -4.0f*(1.0f-f)*(dx*dphidxC + dy*dphidyC + dz*dphidzC );

    r = (phiC - phiU)/(phiD - phiC);

    //limit = sign_sano(r)*max(0.0, (min(abs(r), sign_sano(r))));
    limit = max(0.0f, min(1.0f,r));

    phif = phiC + limit*(dphidxC*cpdx +dphidyC*cpdy +dphidzC*cpdz);

    return phif;
};

__device__ __forceinline__ flow_float interp_dispatch(int scheme, int limit_scheme,
                                                      flow_float phiC, flow_float phiD,
                                                      flow_float dphidxC, flow_float dphidyC, flow_float dphidzC,
                                                      flow_float dphidxD, flow_float dphidyD, flow_float dphidzD,
                                                      flow_float dx, flow_float dy, flow_float dz,
                                                      flow_float cpdx, flow_float cpdy, flow_float cpdz,
                                                      flow_float f, flow_float limiter)
{
    if (scheme == 0 || scheme == -1) {
        return interp_1stUp(scheme, limit_scheme,
                            phiC, phiD,
                            dphidxC, dphidyC, dphidzC,
                            dphidxD, dphidyD, dphidzD,
                            dx, dy, dz,
                            cpdx, cpdy, cpdz,
                            f, limiter);
    }

    if (scheme == 1 && limit_scheme >= 0) {
        return interp_MUSCL_2nd(scheme, limit_scheme,
                                phiC, phiD,
                                dphidxC, dphidyC, dphidzC,
                                dphidxD, dphidyD, dphidzD,
                                dx, dy, dz,
                                cpdx, cpdy, cpdz,
                                f, limiter);
    }

    if (scheme == 2 && limit_scheme >= 0) {
        return interp_MUSCL_3rd(scheme, limit_scheme,
                                phiC, phiD,
                                dphidxC, dphidyC, dphidzC,
                                dphidxD, dphidyD, dphidzD,
                                dx, dy, dz,
                                cpdx, cpdy, cpdz,
                                f, limiter);
    }

    return interp_MINMOD(scheme, limit_scheme,
                         phiC, phiD,
                         dphidxC, dphidyC, dphidzC,
                         dphidxD, dphidyD, dphidzD,
                         dx, dy, dz,
                         cpdx, cpdy, cpdz,
                         f, limiter);
}

__device__ __forceinline__ flow_float apply_ducros_limiter(flow_float limiter, flow_float duc)
{
    if (duc <= 0.8f) {
        return limiter;
    }

    return max(0.0f, (1.0f - duc) * limiter);
}




inline __device__ flow_float sign_sano(flow_float x)
{
    return x > 0 ? 1 : (x<0 ? -1 : 0);
}

__device__ flow_float betaPls_slau(flow_float M)
{
    if (abs(M) >= 1.0f) {
        return 0.25f*(2.0f-M)*sq(M+1.0f);
    } else {
        return 0.5f*(1.0f+sign_sano(+M));
    }
}

__device__ flow_float betaMns_slau(flow_float M)
{
    if (abs(M) >= 1.0f) {
        return 0.25f*(2.0f+M)*sq(M-1.0f);
    } else {
        return 0.5f*(1.0f+sign_sano(-M));
    }
}

// =============================================================================
// 対流フラックス カーネル共通の引数バンドル。
//   - いずれも device ポインタ (とスカラー) を束ねるだけの POD。配列実体はコピーせず、
//     アドレスのみを値渡しする (カーネル引数はもともと値渡し)。
//   - SLAU/HLLE/ROE で完全に同一な mesh幾何/保存量/残差/リミタ/勾配を 5 つに束ね、
//     カーネル冒頭でローカルポインタへ展開して本体を無改変に保つ。
//   - loop_planes は主ループ対象 plane 列 (cell: 全 halo, node: 内部+periodic)。
//     旧シグネチャの nNormal_ghst_Planes (SLAU) / nNormal_halo_Planes (HLLE/ROE) を統一。
// =============================================================================
struct FaceGeom {
    geom_int   nCells, nPlanes, nNormalPlanes;
    geom_int*  plane_cells;
    geom_int   nLoopPlanes;
    geom_int*  loop_planes;
    geom_float *vol, *ccx, *ccy, *ccz;
    geom_float *pcx, *pcy, *pcz, *fx;
    geom_float *sx, *sy, *sz, *ss;
    flow_float* massflux;
    // node の壁ノードフラグ [nCells] (cell 方式・未設定では nullptr)。space.slauWallNormalChi が読む。
    // ゴースト index (>= nCells) には無いので参照前に範囲で弾くこと。
    geom_int*   wall_flag = nullptr;
};
struct PrimState {
    flow_float *ro, *roUx, *roUy, *roUz, *roe;
    flow_float *Ux, *Uy, *Uz, *Ps, *Ht, *sonic;
    // space.reconT=1 で MUSCL 再構成を ρ でなく T に対して行うときの温度場 (0 では未使用)。
    flow_float *T = nullptr;
};
struct ResidualOut {
    flow_float *res_ro, *res_roUx, *res_roUy, *res_roUz, *res_roe;
};
struct LimiterFields {
    flow_float *limiter_ro, *limiter_Ux, *limiter_Uy, *limiter_Uz, *limiter_P;
    flow_float *ducros;
    flow_float *limiter_T = nullptr;   // space.reconT=1 のみ (nullptr なら limiter_P を流用)
};
struct GradFields {
    flow_float *drodx, *drody, *drodz;
    flow_float *dUxdx, *dUxdy, *dUxdz;
    flow_float *dUydx, *dUydy, *dUydz;
    flow_float *dUzdx, *dUzdy, *dUzdz;
    flow_float *dPdx, *dPdy, *dPdz;
    // space.reconT=1 用 (0 では未使用)。
    flow_float *dTdx = nullptr, *dTdy = nullptr, *dTdz = nullptr;
};
// 非平衡凝縮 (二相) のエネルギー流束補正。g_total==nullptr で従来挙動 (ビット不変)。全スキーム共通。
struct CondArgs {
    flow_float  cp_cpg;
    flow_float* g_total;
    flow_float* T_cell;
    int         condModel;
    // SST 全エネルギーに ρk を含める (sstEnergyIncludesK, plan turbulence-sst-energy-includes-k):
    //   kturb = セル k (nullptr で無効)。energyK!=0 のとき面エンタルピーに +(5/3)k、圧力流束に p*=p+(2/3)ρk。
    flow_float* kturb = nullptr;
    int         energyK = 0;
    // 凝縮種物性 (config オプション反映済み) と CPG carrier の Y_w (<=0 で pure)。面状態 T_f=p_f/(ρ_f R_eff) の再構成に使う。
    // (brace 初期化の位置引数順を壊さないよう末尾に置く; 構築後に代入)
    CondSpeciesProps cprops{};
    double      Yw = -1.0;
    // float 経路 (condFloat=1): 面の潜熱 L(T_cell) / L(T_f) を物性表で引く (plans/active/condensation-float-speedup.md §4.2-4)。g=0 の面は評価しない。
    CondTablesF tables{};
    int         condFloat = 0;
};
// thermally-perfect 多成分の化学種データと face 整合再構成用配列 (SLAU の種別再構成で使用)。
struct SpeciesArgs {
    int                  thermalMethod;
    const SpeciesThermo* sp;
    const SpeciesThermoF* spf;   // float32 ミラー (面エンタルピー h_mix(Y_f,T_f) の評価用)
    int                  nSpecies;
    flow_float**         roY;
    flow_float**         Yd_recon;
    flow_float**         dYdx_recon;
    flow_float**         dYdy_recon;
    flow_float**         dYdz_recon;
    flow_float**         limiterY_recon;
    flow_float*          Yface_out;
    flow_float*          Rmix_cell;
    // 受動種 (passiveScalarScheme 1; plan species-passive-scalar-unification §4.1): S3 の面再構成入力。
    //   nPassive 本の原始量 P / 勾配 / 受動種ごとの ψ_P、upwind 面値の出力 Pface_out[ip*nPassive+q]。
    //   nPassiveUnit: 先頭から何本が上限 1 (トレーサ ξ∈[0,1]) を持つか (残りは下限 0 のみ)。
    //   構築後に代入する (brace 初期化の位置引数順を壊さない)。
    int                  nPassive = 0;
    int                  nPassiveUnit = 0;
    flow_float**         P_recon = nullptr;
    flow_float**         dPdx_recon = nullptr;
    flow_float**         dPdy_recon = nullptr;
    flow_float**         dPdz_recon = nullptr;
    flow_float**         limiterP_recon = nullptr;
    flow_float*          Pface_out = nullptr;
};
