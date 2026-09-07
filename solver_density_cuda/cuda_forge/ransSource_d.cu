#include "ransSource_d.cuh"

#include "cuda_forge/cudaWrapper.cuh"

namespace {

constexpr flow_float kBetaStar  = static_cast<flow_float>(0.09);
constexpr flow_float kAlpha1    = static_cast<flow_float>(5.0 / 9.0);
constexpr flow_float kAlpha2    = static_cast<flow_float>(0.44);
constexpr flow_float kBeta1     = static_cast<flow_float>(0.075);
constexpr flow_float kBeta2     = static_cast<flow_float>(0.0828);
constexpr flow_float kSigmaK1   = static_cast<flow_float>(0.85);
constexpr flow_float kSigmaK2   = static_cast<flow_float>(1.0);
constexpr flow_float kSigmaW1   = static_cast<flow_float>(0.5);
constexpr flow_float kSigmaW2   = static_cast<flow_float>(0.856);
constexpr flow_float kA1        = static_cast<flow_float>(0.31);
constexpr flow_float kSmall     = static_cast<flow_float>(1.0e-12);

// SST source terms (production, destruction, cross-diffusion) for k and omega.
// vis_turb にはこの段階での渦粘性が入る（Stage 4 以前はゼロまたは簡易値）。
// 生産項の mu_t が 0 の場合は内部で rho*a1*k/omega で代替計算する。
__global__ void rans_sst_source_d(
    geom_int nCells,
    geom_float* vol,
    flow_float* ro,
    flow_float* k,
    flow_float* omega,
    flow_float* vis_lam,
    flow_float* vis_turb,
    flow_float* wall_dist,
    flow_float* dUxdx, flow_float* dUxdy, flow_float* dUxdz,
    flow_float* dUydx, flow_float* dUydy, flow_float* dUydz,
    flow_float* dUzdx, flow_float* dUzdy, flow_float* dUzdz,
    flow_float* dKdx,     flow_float* dKdy,     flow_float* dKdz,
    flow_float* dOmegadx, flow_float* dOmegady, flow_float* dOmegadz,
    int isAxisymmetric,
    flow_float* axisym_uy_over_r,
    // axisymMethod==1 (SU2 流): k/ω の 1/y 移流拡散ソース (ResidualAxisymmetricConvectionDiffusion)
    // に使う。method 0 では未使用。axis_flag==1 (軸ノード, node のみ) はソース 0。
    int axisymMethod, flow_float* Uy, flow_float* ccy, geom_int* axis_flag_axis,
    int dilatationCorrection,
    int katoLaunder,
    int wallTreatment,
    flow_float* wf_pk,
    flow_float* Pk_diag,   // 診断: 確定した最終 k 生産 (wall-function 置換後)
    // 診断 (§4.1): omega 方程式の項別収支。nullptr 可 (診断オフ)。単位は res_roOmega と同じ (体積込み)。
    flow_float* omg_prod, flow_float* omg_dest, flow_float* omg_cross, flow_float* omg_trans,
    flow_float* omg_axisym,   // 診断: 軸対称 1/y ソース·V (omg_trans はこれが加わる前の残差)
    // E3 (§5): 壁関数第一内層の壁法則整合ひずみ二乗 S_wf^2 (非対象 -1)。omegaWfSource!=0 のときだけ
    // omega 生産の入力に使う (k 生産は既に wf_pk で置換済み)。nullptr / 0 で従来ビット不変。
    flow_float* wf_sprod, int omegaWfSource,
    flow_float* res_roK,
    flow_float* res_roOmega,
    // SST 陰解法 (point-implicit) 用: 消散項ヤコビアン対角（β* ω, 2 β ω）
    flow_float* src_jac_k,
    flow_float* src_jac_omega,
    // SST-DES (methods/turbulence §8): DESmode>0 で k 消滅項を ρk^{3/2}/l_des に切替える。
    int DESmode,
    flow_float* l_des,
    // ω 交差拡散 Jacobian (turbulence-sst-omega-crossdiff-jacobian): 1 で CD_ω>0 のとき対角へ
    // +CD_ω/(ρω) を加える (∂CD_ω/∂(ρω) = -CD_ω/(ρω) < 0 の陰化)。0 = 従来 (ビット不変)。
    int crossDiffJac,
    // node: ω ピン位置 (=壁ノード) の識別用。wf_pk は第一内層にも付くため、ω 残差ゼロ化は壁ノード限定にする。
    // cell では nullptr (wf_pk>=0 で第一セルを識別)。
    geom_int* wall_flag, int isNode,
    // node k Dirichlet (SU2 SetTurbVars_WF 流): roK_wf>=0 の第一内層ノードで res_roK を 0 化 (k は固定なので
    // 残差不要・rms_roK 汚染回避)。nullptr/全-1 で無効 (cell 不変)。
    flow_float* roK_wf,
    flow_float* roOmega_wf,
    // node 入口 (Dirichlet スカラー境界) ピン: ==1 のノードで res_roK/res_roOmega/src_jac_* を 0 化
    // (保存量は ransBoundary で ρ·Dirichlet 値にピン済)。nullptr/全0 で無効 (cell 不変)。
    flow_float* scalarDirichletPin,
    int omegaProdFromPk,
    int energyKSource,
    flow_float* res_roe,
    flow_float* sstF1,
    int nodeWallKPin)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    const flow_float rho    = max(ro[ic], kSmall);
    const flow_float k_c    = max(k[ic],  static_cast<flow_float>(0.0));
    const flow_float w_c    = max(omega[ic], kSmall);
    const flow_float mu_lam = vis_lam[ic];
    const flow_float mu_t   = max(vis_turb[ic], static_cast<flow_float>(0.0));
    const flow_float y      = max(wall_dist[ic], kSmall);
    const geom_float v      = vol[ic];

    // 渦粘性がゼロ（SST eddy viscosity 未実装段階）では k/omega から代替
    const flow_float mu_t_eff = (mu_t > kSmall) ? mu_t : rho * kA1 * k_c / w_c;

    // ひずみ速度の大きさ squared: S^2 = 2 Sij Sij
    const flow_float S11 = dUxdx[ic];
    const flow_float S22 = dUydy[ic];
    const flow_float S33 = dUzdz[ic];
    const flow_float S12 = static_cast<flow_float>(0.5) * (dUxdy[ic] + dUydx[ic]);
    const flow_float S13 = static_cast<flow_float>(0.5) * (dUxdz[ic] + dUzdx[ic]);
    const flow_float S23 = static_cast<flow_float>(0.5) * (dUydz[ic] + dUzdy[ic]);
    flow_float S_sq = static_cast<flow_float>(2.0) *
        (S11*S11 + S22*S22 + S33*S33 +
         static_cast<flow_float>(2.0) * (S12*S12 + S13*S13 + S23*S23));
    // 軸対称: フープひずみ S_thetatheta = u_r/r を加算 (methods/turbulence/theory.md §7.2)
    if (isAxisymmetric) {
        const flow_float S_tt = axisym_uy_over_r[ic];
        S_sq += static_cast<flow_float>(2.0) * S_tt * S_tt;
    }

    // 速度の発散 divU = ∂xUx + ∂yUy + ∂zUz (+ 軸対称 u_r/r) = axisym_divU と一致
    flow_float divU = S11 + S22 + S33;
    if (isAxisymmetric) {
        divU += axisym_uy_over_r[ic];
    }

    // (A) deviatoric トレース除去 (dilatationCorrection >= 1, methods/turbulence/theory.md §7.3)
    //   S^2 = 2 Sij Sij から 2/3 (divU)^2 を引く (= 2 * 1/3 (divU)^2)
    if (dilatationCorrection >= 1) {
        S_sq -= static_cast<flow_float>(2.0 / 3.0) * divU * divU;
        S_sq = max(S_sq, static_cast<flow_float>(0.0));
    }

    // Kato-Launder 補正 (methods/turbulence/theory.md §7.5): 生産の片方の S を渦度 Omega に
    // 置換し、非回転の強加速 (よどみ点・ノズル喉中心線) での偽生産を抑える。
    //   S_prod = S^2 (標準, katoLaunder==0)  /  S*Omega (katoLaunder==1)
    // せん断層 S~Omega では従来同等、非回転 Omega->0 では生産->0。Omega は無旋回でも
    // フープ補正不要 (反対称部、ひずみ S^2 の +2 S_tt^2 と非対称)。
    flow_float S_prod = S_sq;
    if (katoLaunder == 1) {
        const flow_float wx = dUzdy[ic] - dUydz[ic];   // omega = rot u
        const flow_float wy = dUxdz[ic] - dUzdx[ic];
        const flow_float wz = dUydx[ic] - dUxdy[ic];
        const flow_float Om_sq = wx*wx + wy*wy + wz*wz; // = 2 Omega_ij Omega_ij = |omega|^2
        S_prod = sqrt(S_sq) * sqrt(Om_sq);              // S*Omega。積 S_sq*Om_sq は float32 で overflow し得る (相似試験 α=1e-3 で Inf, codex 2026-09-08)
    }

    // 交差拡散項（F1 ブレンドに使用）
    const flow_float grad_k_dot_w =
        dKdx[ic] * dOmegadx[ic] + dKdy[ic] * dOmegady[ic] + dKdz[ic] * dOmegadz[ic];
    const flow_float CD_kw = max(
        static_cast<flow_float>(2.0) * rho * kSigmaW2 / w_c * grad_k_dot_w,
        static_cast<flow_float>(1.0e-10));

    // F1 ブレンド関数
    const flow_float nu     = mu_lam / rho;
    const flow_float arg1_a = sqrt(k_c) / (kBetaStar * w_c * y);
    const flow_float arg1_b = static_cast<flow_float>(500.0) * nu / (w_c * y * y);
    const flow_float arg1_c = static_cast<flow_float>(4.0) * rho * kSigmaW2 * k_c / (CD_kw * y * y);
    const flow_float arg1   = min(max(arg1_a, arg1_b), arg1_c);
    const flow_float F1     = tanh(arg1 * arg1 * arg1 * arg1);
    // sstF1 は ransBlendF1_d_wrapper (拡散の前) が同式・同入力で書いており、ここでは同じ値を再書込するだけ
    // (σ ブレンドと生成の F1 が同一 step で一致することの保証)。
    if (sstF1 != nullptr) sstF1[ic] = F1;

    // ブレンドされた係数
    const flow_float alpha = F1 * kAlpha1 + (static_cast<flow_float>(1.0) - F1) * kAlpha2;
    const flow_float beta  = F1 * kBeta1  + (static_cast<flow_float>(1.0) - F1) * kBeta2;

    // k 生産項（10 beta* rho k omega でリミット）。S_prod は katoLaunder で S^2 / S*Omega。
    // (B) 等方項 -2/3 rho k divU (dilatationCorrection >= 2, theory.md §7.3) は**リミッタの前**に足す:
    //   Menter/Wilcox の圧縮性形は P_k = τ_ij ∂u_i/∂x_j (τ に -2/3ρk δ_ij を含む) をリミットする。旧順序
    //   (リミット後に加算) では最終 P_k が 10β*ρkω を超え得た (codex 2026-09-08, plan
    //   turbulence-sst-consistency-options §2.1)。膨張(divU>0)でシンク, 圧縮(divU<0)でソース。負生産は 0 でクリップ。
    flow_float Pk_raw = mu_t_eff * S_prod;
    if (dilatationCorrection >= 2) {
        Pk_raw -= static_cast<flow_float>(2.0 / 3.0) * rho * k_c * divU;
    }
    flow_float Pk = min(Pk_raw, static_cast<flow_float>(10.0) * kBetaStar * rho * k_c * w_c);
    Pk = max(Pk, static_cast<flow_float>(0.0));

    // automatic wall treatment (methods/turbulence §6.5(d)): wall-adjacent セル (wf_pk>=0) では
    // 解像勾配ベースの P_k を wall-function 生産 P_k=ρu_τ⁴/ν·g(1-g) に置換する。ω はピン留め済み
    // (ransBoundary) なので k は平衡値 u_τ²/√β* に自律収束する。非 wall セルは wf_pk<0 で不変。
    if (wallTreatment == 1 && wf_pk[ic] >= static_cast<flow_float>(0.0)) {
        Pk = wf_pk[ic];
    }
    Pk_diag[ic] = Pk;  // 診断: 確定 k 生産 (wall-function 置換後)

    // k 消滅項。標準 SST は D_k = β* ρ k ω = ρ k^{3/2}/l_RANS。
    // SST-DES (DESmode>0) では l_RANS を l_DDES に置換: D_k = ρ k^{3/2}/l_des (methods/turbulence §8)。
    // 剥離域 (l_des < l_RANS) で消滅が強まり渦粘性が下がる。ω 方程式・渦粘性式は不変。
    // DESmode==0 では従来式そのまま (係数・順序・キャスト一致) = ビット不変。
    flow_float Dk;
    flow_float jac_k;  // 陰解法ヤコビアン対角 ∂Dk/∂(ρk)
    if (DESmode > 0) {
        const flow_float l_hyb = max(l_des[ic], static_cast<flow_float>(1.0e-20));
        Dk    = rho * pow(k_c, static_cast<flow_float>(1.5)) / l_hyb;
        jac_k = static_cast<flow_float>(1.5) * sqrt(k_c) / l_hyb;
    } else {
        Dk    = kBetaStar * rho * k_c * w_c;
        jac_k = kBetaStar * w_c;
    }

    // omega 生産項。
    // E3 (§5): 壁関数の第一内層ノードでは入力ひずみを壁法則整合 S_wf^2 に差し替える
    // (k 側は既に wf_pk で壁関数値。omega だけ解像ひずみのままという非対称の解消)。
    flow_float S_prod_omega = S_prod;
    if (omegaWfSource != 0 && wf_sprod != nullptr
        && wf_sprod[ic] >= static_cast<flow_float>(0.0)) {
        S_prod_omega = wf_sprod[ic];
    }
    // sstOmegaProdFromPk=1: P_ω = α P_k/ν_t (k 側のリミッタ・dilatation・壁関数置換後の P_k と整合, TMR/SST-1994 形)。
    // 0 (既定): P_ω = α ρ S² (SST-2003 形, 現行)。リミッタ非発動時は両者一致。
    // ν_t は closure の正本 vis_turb (mu_t) を使う (mu_t_eff の 1e-12 切替は a1 分ずれる, codex 2026-09-08)。μt→0 の極限では
    // P_k→0 かつ P_k/ν_t→ρS² なので αρS_prod にフォールバック (相対閾値: μt < 1e-6 μ)。
    const flow_float Pw = (omegaProdFromPk != 0 && mu_t > static_cast<flow_float>(1.0e-6) * mu_lam)
        ? alpha * rho * Pk / mu_t
        : alpha * rho * S_prod_omega;

    // omega 消滅項
    const flow_float Dw = beta * rho * w_c * w_c;

    // omega 交差拡散ソース（k-epsilon 側でのみ有効: 1 - F1）
    const flow_float CDw = (static_cast<flow_float>(1.0) - F1) *
        static_cast<flow_float>(2.0) * rho * kSigmaW2 / w_c * grad_k_dot_w;

    // 診断 (§4.1): ソース加算 **前** の res_roOmega = 輸送 (対流+拡散) の寄与。
    // 各ソース項も体積込みで保存し、平衡がどの項で決まっているかを直接見られるようにする。
    if (omg_trans != nullptr) {
        omg_trans[ic] = res_roOmega[ic];
        omg_prod[ic]  = Pw  * v;
        omg_dest[ic]  = Dw  * v;
        omg_cross[ic] = CDw * v;
    }

    // 残差に加算（符号: ソース項は保存変数を増加させる正方向）
    atomicAdd(&res_roK[ic],     (Pk - Dk) * v);
    atomicAdd(&res_roOmega[ic], (Pw - Dw + CDw) * v);
    // sstEnergyKSource=1: E=e+u²/2 (k を含まない) の定式化で、k に溜まる (減る) 分 P_k−D_k を平均流エネルギーから
    // 引く (足す)。局所平衡 P_k≈D_k では 0。SU2 型 (E に k を含める) の代替 (plan turbulence-sst-node-corner-heating §3.5)。
    if (energyKSource != 0 && res_roe != nullptr) {
        atomicAdd(&res_roe[ic], -(Pk - Dk) * v);
    }

    // axisymMethod==1 (SU2 流): k/ω の 1/y 移流拡散ソース (SU2 ResidualAxisymmetricConvectionDiffusion
    // の符号反転移植)。S_φ = -(1/y)·(ρv·φ − (μ + σ_φ μ_t)·∂φ/∂y)。軸ノード (axis_flag==1) と y≤eps は 0。
    flow_float axisym_jac_diag = static_cast<flow_float>(0.0);
    if (isAxisymmetric == 1 && axisymMethod == 1) {
        const bool onAxis = (axis_flag_axis != nullptr && axis_flag_axis[ic] == 1);
        const flow_float yc = ccy[ic];
        if (!onAxis && yc > static_cast<flow_float>(1.0e-12)) {
            const flow_float yinv = static_cast<flow_float>(1.0) / yc;
            const flow_float sigK = F1 * kSigmaK1 + (static_cast<flow_float>(1.0) - F1) * kSigmaK2;
            const flow_float sigW = F1 * kSigmaW1 + (static_cast<flow_float>(1.0) - F1) * kSigmaW2;
            const flow_float rhov = rho * Uy[ic];
            const flow_float sk = -yinv * (rhov * k_c - (mu_lam + sigK * mu_t) * dKdy[ic]);
            const flow_float sw = -yinv * (rhov * w_c - (mu_lam + sigW * mu_t) * dOmegady[ic]);
            atomicAdd(&res_roK[ic],     sk * v);
            atomicAdd(&res_roOmega[ic], sw * v);
            if (omg_axisym != nullptr) omg_axisym[ic] = sw * v;
            // point-implicit 対角: -∂S/∂(ρφ) = +v_r/y (SU2 Jacobian と同じ)。対角正値性を守るため非負側のみ。
            axisym_jac_diag = max(yinv * Uy[ic], static_cast<flow_float>(0.0));
        }
    }

    // SST 陰解法用の消散ヤコビアン対角（point-implicit で D に加える正の量）。
    //   Dk = β* ρ k ω = β* (ρk) ω  →  ∂Dk/∂(ρk) = β* ω
    //   Dω = β ρ ω²  = β (ρω)²/ρ   →  ∂Dω/∂(ρω) = 2 β ω
    // 生産項は limiter 付きで bounded のため lagged（陰化しない）。生産振幅を正の対角に足す
    // under-relaxation も試したが、安定 cfl_pseudo を一切上げず（律速は k/ω 輸送項の陽的扱い）
    // 不採用。輸送項の point-implicit 化は scalarTransport_d.cu（transport_diag）を参照。
    src_jac_k[ic]     = jac_k + axisym_jac_diag;  // DESmode==0: β* ω (従来式と一致), DESmode>0: 1.5 √k / l_des
    src_jac_omega[ic] = static_cast<flow_float>(2.0) * beta * w_c + axisym_jac_diag;
    // ω 交差拡散 Jacobian (opt-in): CD_ω>0 のときのみ対角を強める (負側は対角を弱めるため除外, SU2 同様)。
    if (crossDiffJac != 0 && CDw > static_cast<flow_float>(0.0)) {
        src_jac_omega[ic] += CDw / max(rho * w_c, static_cast<flow_float>(1.0e-30));
    }

    // automatic wall treatment: ω ピン留め (Dirichlet, ransBoundary) セルの残差は 0 が正しい。
    // res_roOmega を 0 にし無駄な update と rms_roOmega 汚染を防ぐ。k は解くので残差は残す。
    // ω ピン位置は「壁ノード/壁隣接第一セル」: node は wall_flag==1 (wf_pk は第一内層にも付くので不可)、
    // cell は wf_pk>=0 (=第一セル ic)。第一内層ノード (node, wf_pk>=0 だが非壁) の ω は解くので 0 化しない。
    const bool omega_pinned = (wallTreatment == 1) &&
        (isNode != 0 ? (wall_flag != nullptr && wall_flag[ic] == 1)
                     : (wf_pk[ic] >= static_cast<flow_float>(0.0)));
    if (omega_pinned) {
        res_roOmega[ic]   = static_cast<flow_float>(0.0);
        src_jac_omega[ic] = static_cast<flow_float>(0.0);
    }
    // node 低 Re 壁ノード: k=0 / ω=ω_w とも ransBoundary でピン済 → 残差・対角を 0 化 (更新なし, rms 汚染回避)
    // (sstNodeWallKPin=0 で旧挙動: 壁ノードの ω は ransBoundary で毎 step 再ピンされるだけで残差は残る)
    if (nodeWallKPin != 0 && isNode != 0 && wallTreatment == 0 && wall_flag != nullptr && wall_flag[ic] == 1) {
        res_roOmega[ic]   = static_cast<flow_float>(0.0);
        src_jac_omega[ic] = static_cast<flow_float>(0.0);
        res_roK[ic]       = static_cast<flow_float>(0.0);
        src_jac_k[ic]     = static_cast<flow_float>(0.0);
    }
    // node k Dirichlet (第一内層ノード): k は roK_wf に固定するので残差・対角は不要 (rms_roK 汚染回避)。
    if (roK_wf != nullptr && roK_wf[ic] >= static_cast<flow_float>(0.0)) {
        res_roK[ic]   = static_cast<flow_float>(0.0);
        src_jac_k[ic] = static_cast<flow_float>(0.0);
    }
    // node ω Dirichlet (第一内層ノード, roOmega_wf>=0): ω は固定なので残差・対角不要 (k と対)。
    if (roOmega_wf != nullptr && roOmega_wf[ic] >= static_cast<flow_float>(0.0)) {
        res_roOmega[ic]   = static_cast<flow_float>(0.0);
        src_jac_omega[ic] = static_cast<flow_float>(0.0);
    }
    // node 入口 (Dirichlet スカラー境界) ノード: k/ω とも ρ·Dirichlet 値にピン済なので残差・対角を 0 化
    // (放置すると入口×壁コーナーで res_roOmega が突出し rms_roOmega プラトーを作る)。
    if (scalarDirichletPin != nullptr && scalarDirichletPin[ic] == static_cast<flow_float>(1.0)) {
        res_roK[ic]       = static_cast<flow_float>(0.0);
        res_roOmega[ic]   = static_cast<flow_float>(0.0);
        src_jac_k[ic]     = static_cast<flow_float>(0.0);
        src_jac_omega[ic] = static_cast<flow_float>(0.0);
    }
}

// F1 前処理 (σ ブレンド用)。rans_sst_source_d の F1 と**同じ式・同じ演算順**で書くこと (両者が一致する前提)。
__global__ void rans_sst_blend_f1_d(
    geom_int nCells,
    flow_float* ro, flow_float* k, flow_float* omega, flow_float* vis_lam, flow_float* wall_dist,
    flow_float* dKdx,     flow_float* dKdy,     flow_float* dKdz,
    flow_float* dOmegadx, flow_float* dOmegady, flow_float* dOmegadz,
    flow_float* sstF1)
{
    geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    const flow_float rho    = max(ro[ic], kSmall);
    const flow_float k_c    = max(k[ic],  static_cast<flow_float>(0.0));
    const flow_float w_c    = max(omega[ic], kSmall);
    const flow_float mu_lam = vis_lam[ic];
    const flow_float y      = max(wall_dist[ic], kSmall);
    const flow_float grad_k_dot_w =
        dKdx[ic] * dOmegadx[ic] + dKdy[ic] * dOmegady[ic] + dKdz[ic] * dOmegadz[ic];
    const flow_float CD_kw = max(
        static_cast<flow_float>(2.0) * rho * kSigmaW2 / w_c * grad_k_dot_w,
        static_cast<flow_float>(1.0e-10));
    const flow_float nu     = mu_lam / rho;
    const flow_float arg1_a = sqrt(k_c) / (kBetaStar * w_c * y);
    const flow_float arg1_b = static_cast<flow_float>(500.0) * nu / (w_c * y * y);
    const flow_float arg1_c = static_cast<flow_float>(4.0) * rho * kSigmaW2 * k_c / (CD_kw * y * y);
    const flow_float arg1   = min(max(arg1_a, arg1_b), arg1_c);
    sstF1[ic] = tanh(arg1 * arg1 * arg1 * arg1);
}

}

void ransBlendF1_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!(cfg.LESorRANS == 2 && cfg.RANSmodel == 1)) return;
    if (!var.c_d.count("sstF1")) return;
    rans_sst_blend_f1_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        msh.nCells,
        var.c_d["ro"], var.c_d["k"], var.c_d["omega"], var.c_d["vis_lam"], var.c_d["wall_dist"],
        var.c_d["dKdx"],     var.c_d["dKdy"],     var.c_d["dKdz"],
        var.c_d["dOmegadx"], var.c_d["dOmegady"], var.c_d["dOmegadz"],
        var.c_d["sstF1"]);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchkKernelSync();
}

void ransSource_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!(cfg.LESorRANS == 2 && cfg.RANSmodel == 1)) {
        return;
    }

    // v は res_roK/res_roOmega の体積ソース加算にのみ使われる → 周期 node では部分体積を渡す
    // (merged volume だと periodicNodeGather の合算で seam 2 倍/コーナー 4 倍。mesh.hpp 参照)。
    rans_sst_source_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        msh.nCells,
        (msh.volumePartial_d != nullptr) ? msh.volumePartial_d : var.c_d["volume"],
        var.c_d["ro"],
        var.c_d["k"],
        var.c_d["omega"],
        var.c_d["vis_lam"],
        var.c_d["vis_turb"],
        var.c_d["wall_dist"],
        var.c_d["dUxdx"], var.c_d["dUxdy"], var.c_d["dUxdz"],
        var.c_d["dUydx"], var.c_d["dUydy"], var.c_d["dUydz"],
        var.c_d["dUzdx"], var.c_d["dUzdy"], var.c_d["dUzdz"],
        var.c_d["dKdx"],     var.c_d["dKdy"],     var.c_d["dKdz"],
        var.c_d["dOmegadx"], var.c_d["dOmegady"], var.c_d["dOmegadz"],
        cfg.isAxisymmetric,
        var.c_d["axisym_uy_over_r"],
        // 診断 (env): FORGE_DIAG_SU2SST_OFF=1 で k/ω の 1/y 軸対称ソースを落とす。
        (getenv("FORGE_DIAG_SU2SST_OFF") != nullptr) ? 0 : cfg.axisymMethod,
        var.c_d["Uy"],
        var.c_d["ccy"],
        (cfg.discretization == "node") ? msh.axis_flag_d : nullptr,
        cfg.dilatationCorrection,
        cfg.katoLaunder,
        cfg.wallTreatmentSST,
        var.c_d["wf_pk"],
        var.c_d["Pk_diag"],
        var.c_d.count("omg_prod")  ? var.c_d["omg_prod"]  : nullptr,
        var.c_d.count("omg_dest")  ? var.c_d["omg_dest"]  : nullptr,
        var.c_d.count("omg_cross") ? var.c_d["omg_cross"] : nullptr,
        var.c_d.count("omg_trans") ? var.c_d["omg_trans"] : nullptr,
        var.c_d.count("omg_axisym") ? var.c_d["omg_axisym"] : nullptr,
        var.c_d.count("wf_sprod") ? var.c_d["wf_sprod"] : nullptr,
        // E3 ゲート: FORGE_WF_OMEGA_SOURCE=1 で有効 (既定 0 = 従来ビット不変)。
        // ★ node SST 壁関数が有効なときだけ許可する。そうでないと initWallFunctionPk_d_wrapper が
        //   早期 return して wf_sprod が確保時の 0 のまま残り、`wf_sprod>=0` が全 DOF で真になって
        //   **omega 生産を全域で消してしまう** (limiter env と同じ落とし穴)。
        (!(cfg.discretization == "node" && cfg.LESorRANS == 2 && cfg.RANSmodel == 1
           && cfg.wallTreatmentSST == 1))
        ? 0
        : [] {
            static const int v = [] {
                const char* s = std::getenv("FORGE_WF_OMEGA_SOURCE");
                const int x = s ? std::atoi(s) : 0;
                if (x) printf("[EXPERIMENT] E3: omega production uses wall-law S_wf^2 "
                              "at the first interior node (FORGE_WF_OMEGA_SOURCE=%d)\n", x);
                return x;
            }();
            return v;
        }(),
        var.c_d["res_roK"],
        var.c_d["res_roOmega"],
        var.c_d["src_jac_k"],
        var.c_d["src_jac_omega"],
        cfg.DESmode,
        var.c_d["l_des"],
        cfg.sstCrossDiffJac,
        (cfg.discretization == "node") ? msh.wall_flag_d : nullptr,
        (cfg.discretization == "node") ? 1 : 0,
        (cfg.discretization == "node" && cfg.wallTreatmentSST == 1 && cfg.nodeKwfDirichlet == 1) ? var.c_d["roK_wf"] : nullptr,
        (cfg.discretization == "node" && cfg.wallTreatmentSST == 1 && cfg.nodeOmegaWfDirichlet == 1) ? var.c_d["roOmega_wf"] : nullptr,
        (cfg.discretization == "node") ? var.c_d["scalarDirichletPin"] : nullptr,
        cfg.sstOmegaProdFromPk,
        cfg.sstEnergyKSource,
        var.c_d["res_roe"],
        var.c_d.count("sstF1") ? var.c_d["sstF1"] : nullptr,
        cfg.sstNodeWallKPin);

    gpuErrchk(cudaPeekAtLastError());
    gpuErrchkKernelSync();
}
