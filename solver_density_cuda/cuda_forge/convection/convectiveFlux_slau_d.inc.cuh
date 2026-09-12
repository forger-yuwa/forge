// =============================================================================
// 対流フラックス SLAU_d の実装断片。
//   - convectiveFlux_d.cu に textual include される断片であり **スタンドアロン TU ではない**。
//     単一 TU を維持するため (rdc 不要 / __device__ グローバル・cudaMemcpyToSymbol の共有を保つ)、
//     CMake の source には .cu/.cuh として足さないこと。
//   - 参照する共通ヘルパ (interp_dispatch 等) と診断グローバルは include 元の上方で定義済み。
// =============================================================================

__global__ void SLAU_d
(
 int conv_scheme, int limit_scheme,
 int slauVariant,   // 1: SLAU, 2: SLAU2 (圧力束第3項のみ低マッハ改良)
 int lowMachPrecond, flow_float precondEps,   // 低マッハ前処理 (1: 散逸スケールを c'、0: 従来 c_hat)
 int lowMachThornber,                         // Thornber 再構成補正 (1: L/R 速度ジャンプを z=min(M,1) で縮約)
 flow_float ga,
 SpeciesArgs   spA,    // thermalMethod / sp / nSpecies / roY / Y{,d}_recon / limiterY_recon / Yface_out / Rmix_cell
 CondArgs      cnd,    // cp_cpg / g_total / T_cell / condModel (非平衡凝縮の二相エネルギー補正)
 FaceGeom      geom,   // mesh 構造 (nCells..massflux)
 PrimState     st,     // 保存量/原始量 (ro..sonic)
 ResidualOut   reso,   // 残差出力 (res_*)
 LimiterFields lim,    // リミタ (limiter_* / ducros)
 GradFields    grd     // 勾配 (drod*..dPd*)
)
{
    // --- struct 引数をローカルへ展開 (以降の本体は旧シグネチャのまま無改変) ---
    const int            thermalMethod  = spA.thermalMethod;
    const SpeciesThermo* sp             = spA.sp;
    const int            nSpecies       = spA.nSpecies;
    flow_float**         roY            = spA.roY;
    flow_float**         Yd_recon       = spA.Yd_recon;
    flow_float**         dYdx_recon     = spA.dYdx_recon;
    flow_float**         dYdy_recon     = spA.dYdy_recon;
    flow_float**         dYdz_recon     = spA.dYdz_recon;
    flow_float**         limiterY_recon = spA.limiterY_recon;
    flow_float*          Yface_out      = spA.Yface_out;
    flow_float*          Rmix_cell      = spA.Rmix_cell;

    const flow_float cp_cpg   = cnd.cp_cpg;
    flow_float*      g_total  = cnd.g_total;
    flow_float*      T_cell   = cnd.T_cell;
    const int        condModel= cnd.condModel;
    flow_float*      kturb    = cnd.kturb;      // sstEnergyIncludesK: セル k (nullptr で無効)
    const int        energyK  = cnd.energyK;

    const geom_int nCells       = geom.nCells;
    const geom_int nPlanes      = geom.nPlanes;
    const geom_int nNormalPlanes= geom.nNormalPlanes;
    geom_int*      plane_cells  = geom.plane_cells;
    const geom_int nNormal_ghst_Planes = geom.nLoopPlanes;
    geom_int*      normal_ghst_planes_d = geom.loop_planes;
    geom_float *vol=geom.vol, *ccx=geom.ccx, *ccy=geom.ccy, *ccz=geom.ccz;
    geom_float *pcx=geom.pcx, *pcy=geom.pcy, *pcz=geom.pcz, *fx=geom.fx;
    geom_float *sx=geom.sx, *sy=geom.sy, *sz=geom.sz, *ss=geom.ss;
    flow_float* massflux=geom.massflux;

    flow_float *ro=st.ro, *roUx=st.roUx, *roUy=st.roUy, *roUz=st.roUz, *roe=st.roe;
    flow_float *Ux=st.Ux, *Uy=st.Uy, *Uz=st.Uz, *Ps=st.Ps, *Ht=st.Ht, *sonic=st.sonic;

    flow_float *res_ro=reso.res_ro, *res_roUx=reso.res_roUx, *res_roUy=reso.res_roUy, *res_roUz=reso.res_roUz, *res_roe=reso.res_roe;

    flow_float *limiter_ro=lim.limiter_ro, *limiter_Ux=lim.limiter_Ux, *limiter_Uy=lim.limiter_Uy, *limiter_Uz=lim.limiter_Uz, *limiter_P=lim.limiter_P;
    flow_float* ducros=lim.ducros;

    flow_float *drodx=grd.drodx, *drody=grd.drody, *drodz=grd.drodz;
    flow_float *dUxdx=grd.dUxdx, *dUxdy=grd.dUxdy, *dUxdz=grd.dUxdz;
    flow_float *dUydx=grd.dUydx, *dUydy=grd.dUydy, *dUydz=grd.dUydz;
    flow_float *dUzdx=grd.dUzdx, *dUzdy=grd.dUzdy, *dUzdz=grd.dUzdz;
    flow_float *dPdx=grd.dPdx, *dPdy=grd.dPdy, *dPdz=grd.dPdz;
    // --- ローカル展開ここまで ---

    geom_int ip_orig = blockDim.x*blockIdx.x + threadIdx.x;


    //if (ip < nPlanes) {
    if (ip_orig < nNormal_ghst_Planes) {

        geom_int ip = normal_ghst_planes_d[ip_orig];

        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        // At non-periodic boundary planes, ic1 is a ghost cell (index >= nCells)
        // whose state has been set by applyBconds(). Gradients/limiters are not
        // computed for ghost cells, so we must NOT reconstruct from ghost; force
        // 1st-order upwind here so that the face state equals the ghost state
        // exactly (which encodes the boundary condition).
        if (ic1 >= nCells) {
            conv_scheme = -1;
        }

        // D2a/D4 診断: 多成分 TP 内部面の組成センサ s_Y = max_s|Y_sR-Y_sL|/(Y_sR+Y_sL+ε)。
        flow_float sY_dbg = 0.0f;
        if ((g_contact1st || g_contactLog || g_contactBlend > 0.0f) && nSpecies > 1 && roY != nullptr && ic1 < nCells) {
            const flow_float ro0 = max(ro[ic0], (flow_float)1.0e-30);
            const flow_float ro1 = max(ro[ic1], (flow_float)1.0e-30);
            for (int s = 0; s < nSpecies; ++s) {
                const flow_float yl = roY[s][ic0] / ro0;
                const flow_float yr = roY[s][ic1] / ro1;
                sY_dbg = max(sY_dbg, fabsf(yr - yl) / (yr + yl + (flow_float)1.0e-12));
            }
            // D2a: 組成勾配が強い面で flow(ρ,p,u) を1次化 (既存 ghost 用 conv_scheme=-1 を流用)。species 流束は不変。
            if (g_contact1st && sY_dbg > g_contactThresh) conv_scheme = -1;
        }

        //__syncthreads();

        geom_float f = fx[ip];
        
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

        flow_float ccx_0 = ccx[ic0];
        flow_float ccy_0 = ccy[ic0];
        flow_float ccz_0 = ccz[ic0];

        flow_float ccx_1 = ccx[ic1];
        flow_float ccy_1 = ccy[ic1];
        flow_float ccz_1 = ccz[ic1];

        flow_float dcc_x = ccx_1 - ccx_0;
        flow_float dcc_y = ccy_1 - ccy_0;
        flow_float dcc_z = ccz_1 - ccz_0;

        flow_float dc0p_x = pcx[ip] - ccx_0;
        flow_float dc0p_y = pcy[ip] - ccy_0;
        flow_float dc0p_z = pcz[ip] - ccz_0;

        flow_float dc1p_x = pcx[ip] - ccx_1;
        flow_float dc1p_y = pcy[ip] - ccy_1;
        flow_float dc1p_z = pcz[ip] - ccz_1;
        if (g_reconEdgeMid == 1 && ip < nNormalPlanes) {   // node: 目標点 = エッジ中点 (SU2 流)
            dc0p_x = 0.5f*dcc_x; dc0p_y = 0.5f*dcc_y; dc0p_z = 0.5f*dcc_z;
            dc1p_x = -0.5f*dcc_x; dc1p_y = -0.5f*dcc_y; dc1p_z = -0.5f*dcc_z;
        }

        //flow_float lim_ro = min(limiter_ro[ic0], limiter_ro[ic1]);
        //flow_float lim_Ux = min(limiter_Ux[ic0], limiter_Ux[ic1]);
        //flow_float lim_Uy = min(limiter_Uy[ic0], limiter_Uy[ic1]);
        //flow_float lim_Uz = min(limiter_Uz[ic0], limiter_Uz[ic1]);
        //flow_float lim_P  = min(limiter_P [ic0], limiter_P[ic1]);

        // rho-Y 共通リミタ (opt-in 診断): ψ_ρY = min(ψ_ρ, min_s ψ_Y_s) を ρ と全 species へ共通適用。
        //   off (default) では lim_rho_* = limiter_ro[*] となりビット不変。p・速度は各自のリミタのまま。
        flow_float lim_rho_L = limiter_ro[ic0];
        flow_float lim_rho_R = limiter_ro[ic1];
        if (g_rhoYCommonLim && g_speciesFaceRecon && nSpecies > 1 && limiterY_recon != nullptr) {
            flow_float minY_L = limiter_ro[ic0], minY_R = limiter_ro[ic1];
            for (int s = 0; s < nSpecies; ++s) {
                minY_L = min(minY_L, limiterY_recon[s][ic0]);
                minY_R = min(minY_R, limiterY_recon[s][ic1]);
            }
            lim_rho_L = minY_L;
            lim_rho_R = minY_R;
            // 診断: min を決めたのが ρ か species か / バケット集計 / 最小値。
            if (lim_rho_L < limiter_ro[ic0]) atomicAdd(&g_rhoYMinBySpecies, 1ULL); else atomicAdd(&g_rhoYMinByRho, 1ULL);
            if (lim_rho_R < limiter_ro[ic1]) atomicAdd(&g_rhoYMinBySpecies, 1ULL); else atomicAdd(&g_rhoYMinByRho, 1ULL);
            if (lim_rho_L < 0.01f) atomicAdd(&g_psiRhoY_lt001, 1ULL);
            if (lim_rho_R < 0.01f) atomicAdd(&g_psiRhoY_lt001, 1ULL);
            if (lim_rho_L < 0.1f)  atomicAdd(&g_psiRhoY_lt01, 1ULL);
            if (lim_rho_R < 0.1f)  atomicAdd(&g_psiRhoY_lt01, 1ULL);
            atomicMin(&g_psiRhoY_min_scaled, (int)(min(lim_rho_L, lim_rho_R)*1.0e6f));
        }

        flow_float ro_L = interp_dispatch(conv_scheme, limit_scheme, ro[ic0] , ro[ic1], drodx[ic0], drody[ic0], drodz[ic0], drodx[ic1], drody[ic1], drodz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, lim_rho_L);
        flow_float Ux_L = interp_dispatch(conv_scheme, limit_scheme, Ux[ic0] , Ux[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Ux[ic0]);
        flow_float Uy_L = interp_dispatch(conv_scheme, limit_scheme, Uy[ic0] , Uy[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uy[ic0]);
        flow_float Uz_L = interp_dispatch(conv_scheme, limit_scheme, Uz[ic0] , Uz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uz[ic0]);
        flow_float P_L  = interp_dispatch(conv_scheme, limit_scheme, Ps[ic0] , Ps[ic1], dPdx[ic0] , dPdy[ic0] , dPdz[ic0] , dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_P[ic0]);
        //flow_float Ux_L = roUx_L/ro_L;
        //flow_float Uy_L = roUy_L/ro_L;
        //flow_float Uz_L = roUz_L/ro_L;
        // velocity2_L / h_p はブレンド後に算出するため後段へ移動 (lowMachThornber 対応)。

        flow_float ro_R  = interp_dispatch(conv_scheme, limit_scheme, ro[ic1], ro[ic0], drodx[ic1], drody[ic1], drodz[ic1], drodx[ic0], drody[ic0], drodz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0f-f, lim_rho_R);
        flow_float Ux_R  = interp_dispatch(conv_scheme, limit_scheme, Ux[ic1], Ux[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0f-f, limiter_Ux[ic1]);
        flow_float Uy_R  = interp_dispatch(conv_scheme, limit_scheme, Uy[ic1], Uy[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0f-f, limiter_Uy[ic1]);
        flow_float Uz_R  = interp_dispatch(conv_scheme, limit_scheme, Uz[ic1], Uz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0f-f, limiter_Uz[ic1]);
        flow_float P_R   = interp_dispatch(conv_scheme, limit_scheme, Ps[ic1], Ps[ic0], dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dPdx[ic0] , dPdy[ic0] , dPdz[ic0] ,-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0f-f, limiter_P[ic1]);

        // D2a 連続ブレンド (chatter-free): 強組成勾配ほど flow 再構成をセル値(1次)へ滑らかに寄せる。
        if (g_contactBlend > 0.0f && sY_dbg > 0.0f) {
            const flow_float w = min((flow_float)1.0, sY_dbg / g_contactBlend);
            if (w > (flow_float)0.0) {
                const flow_float w1 = (flow_float)1.0 - w;
                ro_L = w1*ro_L + w*ro[ic0]; Ux_L = w1*Ux_L + w*Ux[ic0]; Uy_L = w1*Uy_L + w*Uy[ic0]; Uz_L = w1*Uz_L + w*Uz[ic0]; P_L = w1*P_L + w*Ps[ic0];
                ro_R = w1*ro_R + w*ro[ic1]; Ux_R = w1*Ux_R + w*Ux[ic1]; Uy_R = w1*Uy_R + w*Uy[ic1]; Uz_R = w1*Uz_R + w*Uz[ic1]; P_R = w1*P_R + w*Ps[ic1];
            }
        }
        // 低マッハ Thornber 再構成補正: L/R 速度ジャンプを z=min(M,1) で縮約し、低マッハで
        // O(1/M) に増大する速度ジャンプ由来の散逸を抑える。lowMachThornber==0 で恒等 (ビット不変)、
        // M>=1 で z=1 (超音速域不変)。圧力 P_L/R・密度 ro_L/R は不変。理論は methods/convection/theory.md。
        if (lowMachThornber == 1) {
            flow_float c_hat_th = 0.5f*(sonic[ic0] + sonic[ic1]);
            flow_float v2L_th = Ux_L*Ux_L + Uy_L*Uy_L + Uz_L*Uz_L;
            flow_float v2R_th = Ux_R*Ux_R + Uy_R*Uy_R + Uz_R*Uz_R;
            flow_float z_th = min(static_cast<flow_float>(1.0), sqrt(0.5f*(v2L_th + v2R_th))/c_hat_th);
            flow_float um, du;
            um = 0.5f*(Ux_L+Ux_R); du = 0.5f*(Ux_L-Ux_R); Ux_L = um + z_th*du; Ux_R = um - z_th*du;
            um = 0.5f*(Uy_L+Uy_R); du = 0.5f*(Uy_L-Uy_R); Uy_L = um + z_th*du; Uy_R = um - z_th*du;
            um = 0.5f*(Uz_L+Uz_R); du = 0.5f*(Uz_L-Uz_R); Uz_L = um + z_th*du; Uz_R = um - z_th*du;
        }
        // velocity2_L / h_p はブレンド後に算出 (L 再構成直後から移動)。
        flow_float velocity2_L = Ux_L*Ux_L + Uy_L*Uy_L + Uz_L*Uz_L;
        flow_float velocity2_R = Ux_R*Ux_R + Uy_R*Uy_R + Uz_R*Uz_R;
        flow_float h_p, h_m;
        double YLf_s3[THERMO_MAX_SPECIES], YRf_s3[THERMO_MAX_SPECIES];   // S3: 再構成 face 組成の保存
        bool haveYf_s3 = false;
        if (thermalMethod == 2) {
            // TP gas: 被移流全エンタルピーを NASA 絶対エンタルピーで再構成する。
            //   ga/(ga-1)·P/ρ = cp·T であって h(T)=∫cp dT' ではないため誤り。
            //   面温度 T_face = P_face/(ρ_face·R_mix) → h_mix(T_face) を NASA で評価し ek を加える。
            //   多成分 (roY!=nullptr): L 側は ic0, R 側は ic1 のセル組成 Y を使う (組成は 1 次)。
            //   He/N2 の質量基準 cp は ~5 倍違うため、sp[0] 固定だと混合 contact で
            //   エネルギー束が不整合になり圧力発散する。
            if (nSpecies > 1 && roY != nullptr) {
                double YL[THERMO_MAX_SPECIES], YR[THERMO_MAX_SPECIES];
                const double roL_c = (double)max(ro[ic0], (flow_float)1.0e-30);
                const double roR_c = (double)max(ro[ic1], (flow_float)1.0e-30);
                double sL=0.0, sR=0.0;
                for (int s=0;s<nSpecies;s++){
                    double yl=(double)roY[s][ic0]/roL_c; if(yl<0.0)yl=0.0; YL[s]=yl; sL+=yl;
                    double yr=(double)roY[s][ic1]/roR_c; if(yr<0.0)yr=0.0; YR[s]=yr; sR+=yr;
                }
                const double iL=1.0/(sL>1.0e-30?sL:1.0e-30), iR=1.0/(sR>1.0e-30?sR:1.0e-30);
                for (int s=0;s<nSpecies;s++){ YL[s]*=iL; YR[s]*=iR; }
                double RgL, RgR;
                if (g_speciesFaceRecon && Yd_recon != nullptr && g_rhoYCommonLim) {
                    // rho-Y 共通リミタ (opt-in 診断): ρ と同一の ψ_ρY=min(ψ_ρ,min_s ψ_Y) (lim_rho_L/R) で Y を再構成。
                    //   → ρ_f と Y_f が常に同次数 = ρ_f=ρ(Y_f) の熱力学整合だけを切り分ける。
                    //   boundedness: Y∉[0,1] の face-side は Y だけ clamp せず、その side 全体 (ρ と全 Y) を
                    //   セル値 (1次) へ fallback して ρ-Y 整合を保つ (fallback 数を記録)。
                    double cellYL[THERMO_MAX_SPECIES], cellYR[THERMO_MAX_SPECIES];
                    for (int s=0;s<nSpecies;s++){ cellYL[s]=YL[s]; cellYR[s]=YR[s]; }
                    bool ovL=false, ovR=false;
                    for (int s=0;s<nSpecies;s++){
                        flow_float ylr = interp_dispatch(conv_scheme, limit_scheme, Yd_recon[s][ic0], Yd_recon[s][ic1],
                            dYdx_recon[s][ic0], dYdy_recon[s][ic0], dYdz_recon[s][ic0],
                            dYdx_recon[s][ic1], dYdy_recon[s][ic1], dYdz_recon[s][ic1],
                            dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, lim_rho_L);
                        flow_float yrr = interp_dispatch(conv_scheme, limit_scheme, Yd_recon[s][ic1], Yd_recon[s][ic0],
                            dYdx_recon[s][ic1], dYdy_recon[s][ic1], dYdz_recon[s][ic1],
                            dYdx_recon[s][ic0], dYdy_recon[s][ic0], dYdz_recon[s][ic0],
                            -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0f-f, lim_rho_R);
                        if (ylr < 0.0f || ylr > 1.0f) ovL = true;
                        if (yrr < 0.0f || yrr > 1.0f) ovR = true;
                        if (ylr < 0.0f || ylr > 1.0f || yrr < 0.0f || yrr > 1.0f) atomicAdd(&g_speciesOvershoot, 1ULL);
                        YL[s]=(double)ylr; YR[s]=(double)yrr;
                    }
                    if (ovL) { ro_L = ro[ic0]; for (int s=0;s<nSpecies;s++) YL[s]=cellYL[s]; atomicAdd(&g_rhoYFallback, 1ULL); }
                    else { double sLr=0.0; for(int s=0;s<nSpecies;s++) sLr+=YL[s]; const double iLr=1.0/(sLr>1.0e-30?sLr:1.0e-30); for(int s=0;s<nSpecies;s++) YL[s]*=iLr; }
                    if (ovR) { ro_R = ro[ic1]; for (int s=0;s<nSpecies;s++) YR[s]=cellYR[s]; atomicAdd(&g_rhoYFallback, 1ULL); }
                    else { double sRr=0.0; for(int s=0;s<nSpecies;s++) sRr+=YR[s]; const double iRr=1.0/(sRr>1.0e-30?sRr:1.0e-30); for(int s=0;s<nSpecies;s++) YR[s]*=iRr; }
                    for (int s=0;s<nSpecies;s++){
                        atomicMin(&g_Yface_min_scaled, (int)(min(YL[s],YR[s])*1.0e6f));
                        atomicMax(&g_Yface_max_scaled, (int)(max(YL[s],YR[s])*1.0e6f));
                    }
                    RgL = thermo_R_mix(sp, nSpecies, YL);
                    RgR = thermo_R_mix(sp, nSpecies, YR);
                    for (int s=0;s<nSpecies;s++){ YLf_s3[s]=YL[s]; YRf_s3[s]=YR[s]; }
                    haveYf_s3 = true;
                } else if (g_speciesFaceRecon && Yd_recon != nullptr) {
                    // proper S2/S3: Y_s を ρ と同一リミタ (limiter_ro) で face へ再構成 (ρ と整合 + Y の boundedness)。
                    // positivity は clamp([0,1])+正規化 (ΣY=1) で閉じる。clamp 前に overshoot をカウント。
                    double sLr=0.0, sRr=0.0;
                    for (int s=0;s<nSpecies;s++){
                        // ρ と同一リミタ (limiter_ro) を使う = ρ と Y が常に同次数 (consistent)。
                        // min(ψ_ρ,ψ_Y) は Y を ρ より強く落とし ρ-Y 不整合を生み case/28 で発散したため不採用。
                        // 共通 min を試すには multispeciesRhoYCommonLimiter=1 (上の分岐)。
                        const flow_float limL = limiter_ro[ic0];
                        const flow_float limR = limiter_ro[ic1];
                        flow_float ylr = interp_dispatch(conv_scheme, limit_scheme, Yd_recon[s][ic0], Yd_recon[s][ic1],
                            dYdx_recon[s][ic0], dYdy_recon[s][ic0], dYdz_recon[s][ic0],
                            dYdx_recon[s][ic1], dYdy_recon[s][ic1], dYdz_recon[s][ic1],
                            dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limL);
                        flow_float yrr = interp_dispatch(conv_scheme, limit_scheme, Yd_recon[s][ic1], Yd_recon[s][ic0],
                            dYdx_recon[s][ic1], dYdy_recon[s][ic1], dYdz_recon[s][ic1],
                            dYdx_recon[s][ic0], dYdy_recon[s][ic0], dYdz_recon[s][ic0],
                            -dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0f-f, limR);
                        if (ylr < 0.0f || ylr > 1.0f || yrr < 0.0f || yrr > 1.0f) atomicAdd(&g_speciesOvershoot, 1ULL);
                        ylr = min(max(ylr,(flow_float)0.0),(flow_float)1.0);
                        yrr = min(max(yrr,(flow_float)0.0),(flow_float)1.0);
                        YL[s]=(double)ylr; YR[s]=(double)yrr; sLr+=(double)ylr; sRr+=(double)yrr;
                    }
                    const double iLr=1.0/(sLr>1.0e-30?sLr:1.0e-30), iRr=1.0/(sRr>1.0e-30?sRr:1.0e-30);
                    for (int s=0;s<nSpecies;s++){ YL[s]*=iLr; YR[s]*=iRr; }
                    RgL = thermo_R_mix(sp, nSpecies, YL);
                    RgR = thermo_R_mix(sp, nSpecies, YR);
                    // S3: 正規化済み再構成 face 組成 (ΣY=1) を保存し、mdot 確定後に upwind を Yface_out へ。
                    for (int s=0;s<nSpecies;s++){ YLf_s3[s]=YL[s]; YRf_s3[s]=YR[s]; }
                    haveYf_s3 = true;
                } else if (g_faceThermoY) {
                    // S2: 整合 face 組成 Y_face=f·Y_ic0+(1-f)·Y_ic1 (正規化) を L/R thermo に使い、
                    // ρ_L,P_L (2次) と同じ face 位置の組成に揃える。R_mix も Y_face から作る。
                    double Yf[THERMO_MAX_SPECIES], sf=0.0;
                    for (int s=0;s<nSpecies;s++){ Yf[s] = (double)f*YL[s] + (1.0-(double)f)*YR[s]; sf+=Yf[s]; }
                    const double inv = 1.0/(sf>1.0e-30?sf:1.0e-30);
                    for (int s=0;s<nSpecies;s++){ Yf[s]*=inv; YL[s]=Yf[s]; YR[s]=Yf[s]; }
                    RgL = RgR = thermo_R_mix(sp, nSpecies, Yf);
                } else {
                    // R_mix はセル組成のみの定数 (現行 mixed-order)。dependentVariables で計算済の per-cell 値。
                    RgL = (double)Rmix_cell[ic0];
                    RgR = (double)Rmix_cell[ic1];
                }
                // 二相 (凝縮 g>0): 状態の圧力は蒸気のみ P=ρT(R_mix−gR_w) (dependentVariables) なので、
                // 面温度の再構成も同じ気体定数で行う。全水分を気相と数えた R_mix で T=P/(ρR_mix) と
                // すると T が g R_w/R_mix (~2 %) 低く出て面エンタルピーが状態より cp·ΔT (~6 kJ/kg) 小さく
                // なり、凝縮帯で全エンタルピーが +0.6 % 増える非保存が出た (case/44 va2 M4.75: 軸 h0
                // 836→841 kJ/kg、出口ノードだけ境界流束 [Ht 直読] で整合し T が 5 K 低い「跳ね」)。
                if (g_total != nullptr) {
                    const CondSpeciesProps cprR = (condModel == 1) ? condProps_H2O() : condProps_N2();
                    RgL -= (double)g_total[ic0]*cprR.R; if (RgL < 1.0) RgL = 1.0;
                    RgR -= (double)g_total[ic1]*cprR.R; if (RgR < 1.0) RgR = 1.0;
                }
                const double Tl = (double)P_L/((double)ro_L*RgL);
                const double Tr = (double)P_R/((double)ro_R*RgR);
                h_p = (flow_float)(thermo_h_mix(sp, nSpecies, YL, Tl) + 0.5*(double)velocity2_L);
                h_m = (flow_float)(thermo_h_mix(sp, nSpecies, YR, Tr) + 0.5*(double)velocity2_R);
            } else {
                double Rg = thermo_R_species(sp[0]);
                double RgL1 = Rg, RgR1 = Rg;
                if (g_total != nullptr) {   // pure-condensible TP: P=ρ(1−g)RT
                    RgL1 = Rg*(1.0 - (double)g_total[ic0]); if (RgL1 < 1.0) RgL1 = 1.0;
                    RgR1 = Rg*(1.0 - (double)g_total[ic1]); if (RgR1 < 1.0) RgR1 = 1.0;
                }
                const double Tl = (double)P_L/((double)ro_L*RgL1);
                const double Tr = (double)P_R/((double)ro_R*RgR1);
                h_p = (flow_float)(thermo_h_mass(sp[0], Tl) + 0.5*(double)velocity2_L);
                h_m = (flow_float)(thermo_h_mass(sp[0], Tr) + 0.5*(double)velocity2_R);
            }
            // carrier+condensible (H2O in N2) 二相補正: 凝縮した水の潜熱を引く
            //   h_2phase = h_gas^全蒸気 - g L_w + ek。気相混合 h は全蒸気なので -g L_w を足すだけ。
            if (g_total != nullptr) {
                const CondSpeciesProps cprL = (condModel == 1) ? condProps_H2O() : condProps_N2();
                h_p -= (flow_float)((double)g_total[ic0]*cond_latent(cprL, (double)T_cell[ic0]));
                h_m -= (flow_float)((double)g_total[ic1]*cond_latent(cprL, (double)T_cell[ic1]));
            }
        } else {
            h_p = ga*P_L/((ga-1.0f)*ro_L) + 0.5f*velocity2_L;
            h_m = ga*P_R/((ga-1.0f)*ro_R) + 0.5f*velocity2_R;
            // 非平衡凝縮 (二相): 単相 h=cp(1-g)T+ek を二相全エンタルピー h=cpT-gL+ek に補正。
            //   差 = g(cpT - L)。これを落とすとエネルギー流束が潜熱分を運ばず全エンタルピー非保存になる。
            //   c_hat (音速) は気相近似で単相のまま。g・T はセル値(1次, g は元来1次風上移流で整合)。
            //   g_total==nullptr (凝縮 off) で従来 (ビット不変)。
            if (g_total != nullptr) {
                const CondSpeciesProps cprC = (condModel == 1) ? condProps_H2O() : condProps_N2();
                const flow_float gL0 = g_total[ic0], gR0 = g_total[ic1];
                const flow_float TL0 = T_cell[ic0],  TR0 = T_cell[ic1];
                h_p += gL0*(cp_cpg*TL0 - (flow_float)cond_latent(cprC, (double)TL0));
                h_m += gR0*(cp_cpg*TR0 - (flow_float)cond_latent(cprC, (double)TR0));
            }
        }

        //flow_float Vn_p = ((Ux[ic0])*sxx +(Uy[ic0])*syy +(Uz[ic0])*szz)/sss;
        //flow_float Vn_m = ((Ux[ic1])*sxx +(Uy[ic1])*syy +(Uz[ic1])*szz)/sss;

        //flow_float roVn_p = (roUx[ic0]*sxx + roUy[ic0]*syy + roUz[ic0]*szz)/sss;
        //flow_float roVn_m = (roUx[ic1]*sxx + roUy[ic1]*syy + roUz[ic1]*szz)/sss;

        //flow_float VnL = ((Ux[ic0])*sxx +(Uy[ic0])*syy +(Uz[ic0])*szz)/sss;
        //flow_float VnR = ((Ux[ic1])*sxx +(Uy[ic1])*syy +(Uz[ic1])*szz)/sss;
        // SST 全エネルギー E_t = E_m + ρk (sstEnergyIncludesK, plan turbulence-sst-energy-includes-k §4):
        //   エネルギー流束は H* = h + u²/2 + k + (2/3)k (E_t + p*)/ρ、圧力流束は p* = p + (2/3)ρk。
        //   k はセル値 (1 次, k 式の 1 次風上移流と同じ upwind で運ぶ)。熱力学 (面温度・音速) は p のまま。
        flow_float Pf_L = P_L, Pf_R = P_R;   // 圧力流束用 (p*)
        if (energyK != 0 && kturb != nullptr) {
            const flow_float kL = max(kturb[ic0], (flow_float)0.0);
            // ghost (ic1>=nCells) の k は node では書かれない (ghostless) ので内部値で代用 (壁/出口は Neumann、入口は境界カーネル側の bvar k が担う)
            const flow_float kR = (ic1 < nCells) ? max(kturb[ic1], (flow_float)0.0) : kL;
            h_p  += (flow_float)(5.0f/3.0f)*kL;
            h_m  += (flow_float)(5.0f/3.0f)*kR;
            Pf_L += (flow_float)(2.0f/3.0f)*ro_L*kL;
            Pf_R += (flow_float)(2.0f/3.0f)*ro_R*kR;
        }

        flow_float Vn_p = ((Ux_L)*sxx +(Uy_L)*syy +(Uz_L)*szz)/sss;
        flow_float Vn_m = ((Ux_R)*sxx +(Uy_R)*syy +(Uz_R)*szz)/sss;

        flow_float Vn_p_abs = abs(Vn_p);
        flow_float Vn_m_abs = abs(Vn_m);

//TODO: change c
        flow_float c_hat = 0.5f*(sonic[ic0] + sonic[ic1]);

        // D4 診断: 強組成勾配面の L/R 状態を出力 (mixed-order face-state / limiter chatter 観察)。
        // 一面=一行/step。後処理で ip でフィルタすれば pseudo-step 時系列になる。
        if (g_contactLog && sY_dbg > g_contactLogThresh) {
            const flow_float RgL = (nSpecies > 1 && Rmix_cell != nullptr) ? Rmix_cell[ic0] : (flow_float)0.0;
            const flow_float RgR = (nSpecies > 1 && Rmix_cell != nullptr) ? Rmix_cell[ic1] : (flow_float)0.0;
            const flow_float TL  = (RgL > 0.0f) ? P_L / (ro_L * RgL) : (flow_float)0.0;
            const flow_float TR  = (RgR > 0.0f) ? P_R / (ro_R * RgR) : (flow_float)0.0;
            // mixed-order 誤差: 整合 face 組成 (面平均 Y) で R_mix を作った場合の T との差 ΔT_f^MO。
            flow_float Rmix_face = 0.0f;
            if (nSpecies > 1 && roY != nullptr) {
                const flow_float ro0c = max(ro[ic0],(flow_float)1e-30), ro1c = max(ro[ic1],(flow_float)1e-30);
                for (int s = 0; s < nSpecies; ++s) {
                    const flow_float yf = (flow_float)0.5*(roY[s][ic0]/ro0c + roY[s][ic1]/ro1c);
                    Rmix_face += yf * (flow_float)thermo_R_species(sp[s]);
                }
            }
            const flow_float TL_cons = (Rmix_face > 0.0f) ? P_L/(ro_L*Rmix_face) : (flow_float)0.0;
            const flow_float dT_MO   = TL - TL_cons;
            // limiter (ψ_ρ,ψ_p,ψ_u) と cell 値も出し、再構成が実際に効いていたか (sharp で ψ≈0 か) を判定可能に。
            printf("CLOG ip=%d ic0=%d ic1=%d sY=%.4f roL=%.5f roR=%.5f PL=%.1f PR=%.1f UxL=%.4f UxR=%.4f limP=%.4f limRo=%.4f limUx=%.4f roC=%.5f PC=%.1f UxC=%.4f RgL=%.2f RgR=%.2f TL=%.2f TR=%.2f hL=%.6e cL=%.2f\n",
                   (int)ip, (int)ic0, (int)ic1, (double)sY_dbg, (double)ro_L, (double)ro_R,
                   (double)P_L, (double)P_R, (double)Ux_L, (double)Ux_R, (double)limiter_P[ic0],
                   (double)limiter_ro[ic0], (double)limiter_Ux[ic0],
                   (double)ro[ic0], (double)Ps[ic0], (double)Ux[ic0],
                   (double)RgL, (double)RgR, (double)TL, (double)TR, (double)h_p, (double)sonic[ic0]);
            printf("CMO ip=%d sY=%.4f TL=%.3f TLcons=%.3f dT_MO=%.4f RgL=%.2f Rmix_face=%.2f\n",
                   (int)ip, (double)sY_dbg, (double)TL, (double)TL_cons, (double)dT_MO, (double)RgL, (double)Rmix_face);
        }

        flow_float M_p = Vn_p/c_hat;
        flow_float M_m = Vn_m/c_hat;

        flow_float P_del = Pf_R - Pf_L;   // p* 差 (sstEnergyIncludesK off では P_R−P_L)

        flow_float beta_p, beta_m;

        if (abs(M_p)>=1.0f){
            beta_p = 0.5f*(M_p + abs(M_p))/M_p;
        } else {
            beta_p = 0.25f*sq((M_p+1.0f))*(2.0f-M_p);
        }

        if (abs(M_m)>=1.0f){
            beta_m = 0.5f*(M_m - abs(M_m))/M_m;
        } else {
            beta_m = 0.25f*sq((M_m-1.0f))*(2.0f+M_m);
        }

        flow_float zero = 0.0f;
        flow_float one  = 1.0f;
        flow_float half = 0.5f;

        flow_float g = -max(min(M_p,zero),-one)*min(max(M_m,zero),one);
        flow_float Vn_hat_abs   = (ro_L*Vn_p_abs + ro_R*Vn_m_abs)/(ro_L + ro_R);
        flow_float Vn_hat_p_abs = (1.0f-g)*Vn_hat_abs + g*Vn_p_abs;
        flow_float Vn_hat_m_abs = (1.0f-g)*Vn_hat_abs + g*Vn_m_abs;

        flow_float M_hat = min(one, sqrt(half*(velocity2_R + velocity2_L))/c_hat);
        flow_float chi = (1.0f-M_hat)*(1.0f-M_hat);

        flow_float pressure_sum = Pf_L + Pf_R;
        // 圧力束の第3項のみ slauVariant で分岐 (mdot は SLAU/SLAU2 共通)。
        // SLAU : (1-chi)(beta_p+beta_m-1) * (P_L+P_R)/2   ... 低マッハで消失し圧力散逸が乏しい
        // SLAU2: (beta_p+beta_m-1) * sqrt((|u_L|^2+|u_R|^2)/2) * roBar * c_hat  ... M に比例した低マッハ散逸
        flow_float p_third;
        if (slauVariant == 2) {
            flow_float roBar = half*(ro_L + ro_R);
            p_third = (beta_p+beta_m-one)*sqrt(half*(velocity2_L+velocity2_R))*roBar*c_hat;
        } else {
            p_third = (one-chi)*(beta_p+beta_m-one)*half*pressure_sum;
        }
        flow_float p_tilde = half*pressure_sum + half*(beta_p-beta_m)*(Pf_L-Pf_R) + p_third;

        // 低マッハ前処理: 圧力散逸項のスケール c_hat を前処理音速 c_diss に置き換える。
        // lowMachPrecond==0 では c_diss==c_hat でビット不変。M>=1 でも c'=c_hat に復帰。
        // lowMachPrecond==3 (LHS のみの前処理) も RHS 散逸は触らず c_diss==c_hat (=0 とビット不変)。
        flow_float c_diss = c_hat;
        if (lowMachPrecond == 1 || lowMachPrecond == 2) {   // 1: フラックス散逸前処理のみ / 2: + LHS 完全前処理 (どちらも c' 散逸を使う)
            flow_float velMag_face = sqrt(half*(velocity2_L + velocity2_R));
            flow_float Un_face     = half*(Vn_p + Vn_m);
            c_diss = lowMachCprime(c_hat, velMag_face, Un_face, precondEps);
        }

        flow_float mdot = sss*0.5f*((ro_L*(Vn_p+Vn_hat_p_abs)+ro_R*(Vn_m-Vn_hat_m_abs)) -chi/(c_diss)*P_del);
        massflux[ip] = mdot;

        // S3: 同一の再構成 face 組成 (upwind) を Yface_out[ip*nSpecies+s] へ書き出す。species 移流がこれを読む。
        // → energy 流束と species 流束が同一 face 組成 (ΣY=1) を共有し thermo/advection mismatch を解消。
        if (g_speciesFaceRecon >= 2 && Yface_out != nullptr && haveYf_s3) {
            const bool up0 = (mdot >= (flow_float)0.0);
            for (int s = 0; s < nSpecies; ++s) Yface_out[(size_t)ip*nSpecies + s] = (flow_float)(up0 ? YLf_s3[s] : YRf_s3[s]);
        }

        flow_float res_ro_temp   = mdot;
        flow_float p_tilde_r = p_tilde - d_pRef;   // free-stream 保存: 基準静圧を差し引いて float32 桁落ちを抑制
        flow_float res_roUx_temp = 0.5f*(mdot+abs(mdot))*Ux_L +0.5f*(mdot-abs(mdot))*Ux_R +p_tilde_r*sxx;
        flow_float res_roUy_temp = 0.5f*(mdot+abs(mdot))*Uy_L +0.5f*(mdot-abs(mdot))*Uy_R +p_tilde_r*syy;
        flow_float res_roUz_temp = 0.5f*(mdot+abs(mdot))*Uz_L +0.5f*(mdot-abs(mdot))*Uz_R +p_tilde_r*szz;
        flow_float res_roe_temp  = 0.5f*(mdot+abs(mdot))*h_p +0.5f*(mdot-abs(mdot))*h_m ;

        atomicAdd(&res_ro[ic0]  , -res_ro_temp);
        atomicAdd(&res_roUx[ic0], -res_roUx_temp);
        atomicAdd(&res_roUy[ic0], -res_roUy_temp);
        atomicAdd(&res_roUz[ic0], -res_roUz_temp);
        atomicAdd(&res_roe[ic0] , -res_roe_temp);

        atomicAdd(&res_ro[ic1]  , res_ro_temp);
        atomicAdd(&res_roUx[ic1], res_roUx_temp);
        atomicAdd(&res_roUy[ic1], res_roUy_temp);
        atomicAdd(&res_roUz[ic1], res_roUz_temp);
        atomicAdd(&res_roe[ic1] , res_roe_temp);
    }

    __syncthreads();
}
