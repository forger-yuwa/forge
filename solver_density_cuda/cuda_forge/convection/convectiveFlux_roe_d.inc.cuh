// =============================================================================
// 対流フラックス ROE_d の実装断片。
//   - convectiveFlux_d.cu に textual include される断片であり **スタンドアロン TU ではない**。
//     単一 TU を維持するため (rdc 不要 / __device__ グローバル・cudaMemcpyToSymbol の共有を保つ)、
//     CMake の source には .cu/.cuh として足さないこと。
//   - 参照する共通ヘルパ (interp_dispatch 等) と診断グローバルは include 元の上方で定義済み。
// =============================================================================

__global__ void ROE_d
(
 int conv_scheme, int limit_scheme,
 flow_float ga,
 // SU2 の ENTROPY_FIX_COEFF と同形の固有値下限 (0 = 従来 = forge の Harten 補正のみ)。
 // lam[i] = max(lam[i], coeff*(|Ua|+ca))。SU2 既定は 0.001。
 // 現行 Harten 補正 `0.1*(|Ua|/ca+1)` は次元不整合で実効が ~500 倍弱い (下の該当行を参照)。
 flow_float entropyFixCoeff,
 int thermalMethod,                           // 0: calorically perfect, 2: thermally-perfect (NASA-9)
 const SpeciesThermo* sp, int nSpecies,       // thermally-perfect 用化学種データ
 CondArgs      cnd,    // cp_cpg / g_total / T_cell / condModel (非平衡凝縮の二相エネルギー補正)
 FaceGeom      geom,   // mesh 構造 (nCells..massflux)
 PrimState     st,     // 保存量/原始量 (ro..sonic)
 ResidualOut   reso,   // 残差出力 (res_*)
 LimiterFields lim,    // リミタ (limiter_* / ducros)
 GradFields    grd     // 勾配 (drod*..dPd*)
)
{
    // --- struct 引数をローカルへ展開 (以降の本体は旧シグネチャのまま無改変) ---
    const flow_float cp_cpg   = cnd.cp_cpg;
    flow_float*      g_total  = cnd.g_total;
    flow_float*      T_cell   = cnd.T_cell;
    const int        condModel= cnd.condModel;

    const geom_int nCells       = geom.nCells;
    const geom_int nPlanes      = geom.nPlanes;
    const geom_int nNormalPlanes= geom.nNormalPlanes;
    geom_int*      plane_cells  = geom.plane_cells;
    const geom_int nNormal_halo_Planes = geom.nLoopPlanes;
    geom_int*      normal_halo_planes_d = geom.loop_planes;
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
    if (ip_orig < nNormal_halo_Planes) {

        geom_int ip = normal_halo_planes_d[ip_orig];

        flow_float dro;
        flow_float dP;
        flow_float dUx;
        flow_float dUy;
        flow_float dUz;
        flow_float dU;

        flow_float lam[5];
        flow_float R[5][5];
        flow_float L[5][5];
        flow_float dW[5];
        flow_float difQ1[5];
        flow_float difQ2[5];
        flow_float delQ[5];

        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

        // Force 1st-order upwind at non-periodic boundary planes (ghost on R side).
        // Gradients/limiters are not maintained for ghost cells, so reconstructing
        // from them would corrupt the face state. The ghost value itself encodes
        // the boundary condition and must be used as-is.
        if (ic1 >= nCells) {
            conv_scheme = -1;
        }

        //__syncthreads();

        geom_float f = fx[ip];
        
        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

        geom_float nx = sxx/sss;
        geom_float ny = syy/sss;
        geom_float nz = szz/sss;

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
            dc0p_x = 0.5*dcc_x; dc0p_y = 0.5*dcc_y; dc0p_z = 0.5*dcc_z;
            dc1p_x = -0.5*dcc_x; dc1p_y = -0.5*dcc_y; dc1p_z = -0.5*dcc_z;
        }

        flow_float small_rho = 1.0e-4f;
        flow_float small_p   = 1.0f;
        flow_float small_a2  = 1.0f;

        flow_float duc = max(min(max(ducros[ic0], ducros[ic1]), 1.0), 0.0);
        flow_float limiter_ro_L = apply_ducros_limiter(limiter_ro[ic0], duc);
        flow_float limiter_Ux_L = apply_ducros_limiter(limiter_Ux[ic0], duc);
        flow_float limiter_Uy_L = apply_ducros_limiter(limiter_Uy[ic0], duc);
        flow_float limiter_Uz_L = apply_ducros_limiter(limiter_Uz[ic0], duc);
        flow_float limiter_P_L = apply_ducros_limiter(limiter_P[ic0], duc);
        flow_float limiter_ro_R = apply_ducros_limiter(limiter_ro[ic1], duc);
        flow_float limiter_Ux_R = apply_ducros_limiter(limiter_Ux[ic1], duc);
        flow_float limiter_Uy_R = apply_ducros_limiter(limiter_Uy[ic1], duc);
        flow_float limiter_Uz_R = apply_ducros_limiter(limiter_Uz[ic1], duc);
        flow_float limiter_P_R = apply_ducros_limiter(limiter_P[ic1], duc);

        flow_float ro_L = interp_dispatch(conv_scheme, limit_scheme, ro[ic0] , ro[ic1], drodx[ic0], drody[ic0], drodz[ic0], drodx[ic1], drody[ic1], drodz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_ro_L);
        flow_float Ux_L = interp_dispatch(conv_scheme, limit_scheme, Ux[ic0] , Ux[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Ux_L);
        flow_float Uy_L = interp_dispatch(conv_scheme, limit_scheme, Uy[ic0] , Uy[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uy_L);
        flow_float Uz_L = interp_dispatch(conv_scheme, limit_scheme, Uz[ic0] , Uz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_Uz_L);
        flow_float P_L  = interp_dispatch(conv_scheme, limit_scheme, Ps[ic0] , Ps[ic1], dPdx[ic0] , dPdy[ic0] , dPdz[ic0] , dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dcc_x, dcc_y, dcc_z, dc0p_x, dc0p_y, dc0p_z, f, limiter_P_L);

        ro_L = max(ro_L, small_rho);
        P_L = max(P_L, small_p);

        //flow_float Ux_L = roUx_L/ro_L;
        //flow_float Uy_L = roUy_L/ro_L;
        //flow_float Uz_L = roUz_L/ro_L;
        flow_float v2_L = Ux_L*Ux_L +Uy_L*Uy_L +Uz_L*Uz_L;
        flow_float roe_L, Ht_L, ca_L;
        if (thermalMethod == 2) {
            // TP gas: 内部エネルギー・全エンタルピー・音速を NASA で再構成 (sp[0], 単成分)
            const double Rg = thermo_R_species(sp[0]);
            const double Tl = (double)P_L/((double)ro_L*Rg);
            const double hl = thermo_h_mass(sp[0], Tl);
            const double cpl= thermo_cp_mass(sp[0], Tl);
            const double gl = cpl/((cpl-Rg) > 1.0e-6 ? (cpl-Rg) : 1.0e-6);
            roe_L = (flow_float)((double)ro_L*((hl - Rg*Tl) + 0.5*(double)v2_L));
            Ht_L  = (flow_float)(hl + 0.5*(double)v2_L);
            ca_L  = (flow_float)sqrt(max(gl*Rg*Tl, (double)small_a2));
        } else {
            roe_L = P_L/(ga-1.0) + 0.5*ro_L*v2_L;
            Ht_L  = roe_L/ro_L + P_L/ro_L;
            ca_L  = sqrt(max(ga*P_L/ro_L, small_a2));
            // 非平衡凝縮 (二相): roe/Ht を二相に補正 (差 g(cpT-L)/質量)。g_total==nullptr で従来。
            if (g_total != nullptr) {
                const CondSpeciesProps cprC = (condModel == 1) ? condProps_H2O() : condProps_N2();
                const flow_float dL0 = g_total[ic0]*(cp_cpg*T_cell[ic0] - (flow_float)cond_latent(cprC, (double)T_cell[ic0]));
                roe_L += ro_L*dL0;  Ht_L += dL0;
            }
        }

        flow_float ro_R  = interp_dispatch(conv_scheme, limit_scheme, ro[ic1], ro[ic0], drodx[ic1], drody[ic1], drodz[ic1], drodx[ic0], drody[ic0], drodz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_ro_R); 
        flow_float Ux_R  = interp_dispatch(conv_scheme, limit_scheme, Ux[ic1], Ux[ic0], dUxdx[ic1], dUxdy[ic1], dUxdz[ic1], dUxdx[ic0], dUxdy[ic0], dUxdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Ux_R);
        flow_float Uy_R  = interp_dispatch(conv_scheme, limit_scheme, Uy[ic1], Uy[ic0], dUydx[ic1], dUydy[ic1], dUydz[ic1], dUydx[ic0], dUydy[ic0], dUydz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uy_R);
        flow_float Uz_R  = interp_dispatch(conv_scheme, limit_scheme, Uz[ic1], Uz[ic0], dUzdx[ic1], dUzdy[ic1], dUzdz[ic1], dUzdx[ic0], dUzdy[ic0], dUzdz[ic0],-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_Uz_R);
        flow_float P_R   = interp_dispatch(conv_scheme, limit_scheme, Ps[ic1], Ps[ic0], dPdx[ic1] , dPdy[ic1] , dPdz[ic1] , dPdx[ic0] , dPdy[ic0] , dPdz[ic0] ,-dcc_x, -dcc_y, -dcc_z, dc1p_x, dc1p_y, dc1p_z, 1.0-f, limiter_P_R);

        ro_R = max(ro_R, small_rho);
        P_R = max(P_R, small_p);

        flow_float v2_R = Ux_R*Ux_R +Uy_R*Uy_R +Uz_R*Uz_R;
        flow_float roe_R, Ht_R, ca_R;
        if (thermalMethod == 2) {
            const double Rg = thermo_R_species(sp[0]);
            const double Tr = (double)P_R/((double)ro_R*Rg);
            const double hr = thermo_h_mass(sp[0], Tr);
            const double cpr= thermo_cp_mass(sp[0], Tr);
            const double gr = cpr/((cpr-Rg) > 1.0e-6 ? (cpr-Rg) : 1.0e-6);
            roe_R = (flow_float)((double)ro_R*((hr - Rg*Tr) + 0.5*(double)v2_R));
            Ht_R  = (flow_float)(hr + 0.5*(double)v2_R);
            ca_R  = (flow_float)sqrt(max(gr*Rg*Tr, (double)small_a2));
        } else {
            roe_R = P_R/(ga-1.0) + 0.5*ro_R*v2_R;
            Ht_R  = roe_R/ro_R + P_R/ro_R;
            ca_R  = sqrt(max(ga*P_R/ro_R, small_a2));
            if (g_total != nullptr) {
                const CondSpeciesProps cprC = (condModel == 1) ? condProps_H2O() : condProps_N2();
                const flow_float dR0 = g_total[ic1]*(cp_cpg*T_cell[ic1] - (flow_float)cond_latent(cprC, (double)T_cell[ic1]));
                roe_R += ro_R*dR0;  Ht_R += dR0;
            }
        }

        flow_float U_L = ((Ux_L)*sxx +(Uy_L)*syy +(Uz_L)*szz)/sss;
        flow_float U_R = ((Ux_R)*sxx +(Uy_R)*syy +(Uz_R)*szz)/sss;

        // calc delta value
        dro = ro_R - ro_L;
        dP  = P_R  - P_L;
        dUx = Ux_R - Ux_L;
        dUy = Uy_R - Uy_L;
        dUz = Uz_R - Uz_L;
        dU  = (dUx*sxx +dUy*syy +dUz*szz)/sss;

        // calc roe average //ok
        flow_float sqrt_ro_L = sqrt(ro_L);
        flow_float sqrt_ro_R = sqrt(ro_R);
        flow_float sqrt_ro_sum = max(sqrt_ro_R + sqrt_ro_L, small_rho);
        flow_float roa = sqrt_ro_R*sqrt_ro_L;
        flow_float ua = (sqrt_ro_L*Ux_L + sqrt_ro_R*Ux_R)/sqrt_ro_sum;
        flow_float va = (sqrt_ro_L*Uy_L + sqrt_ro_R*Uy_R)/sqrt_ro_sum;
        flow_float wa = (sqrt_ro_L*Uz_L + sqrt_ro_R*Uz_R)/sqrt_ro_sum;
        flow_float Ha = (sqrt_ro_L*Ht_L + sqrt_ro_R*Ht_R)/sqrt_ro_sum;
        flow_float Ua  = ua*nx +va*ny +wa*nz;

        // Roe 平均状態の有効比熱比 ga_eff と音速 ca。
        //   CPG: ga_eff=ga, ca=sqrt((ga-1)(Ha-ek))。
        //   TP : Roe 平均静エンタルピー h~=Ha-ek から T~ を反転し ga_eff=cp~/(cp~-R),
        //        ca=sqrt(ga_eff·R·T~) (= sqrt(ga_eff·P/ρ))。圧力微分ブロックは κ=ga_eff-1 を
        //        採用 (χ=∂P/∂ρ|_{ρe} 項は stage-A で無視; defect-correction で収束解は不変)。
        flow_float ga_eff = ga;
        flow_float ca;
        const flow_float ek_a = 0.5*(ua*ua + va*va + wa*wa);
        if (thermalMethod == 2) {
            const double Yone = 1.0;   // M1 単成分の質量分率
            const double Rg = thermo_R_species(sp[0]);
            const double ha = (double)Ha - (double)ek_a;
            const double Ta = thermo_T_from_h(sp, 1, &Yone, ha, 300.0, 50.0, 6000.0);
            const double cpa = thermo_cp_mass(sp[0], Ta);
            ga_eff = (flow_float)(cpa/((cpa-Rg) > 1.0e-6 ? (cpa-Rg) : 1.0e-6));
            ca = (flow_float)sqrt(max((double)ga_eff*Rg*Ta, (double)small_a2));
        } else {
            ca = sqrt(max((ga-1.0)*(Ha-ek_a), small_a2));
        }

        flow_float P_ro = 0.5*(ga_eff-1.0)*(ua*ua + va*va + wa*wa);
        flow_float P_roe = ga_eff-1.0;
        flow_float P_rou = -ua*(ga_eff-1.0);
        flow_float P_rov = -va*(ga_eff-1.0);
        flow_float P_row = -wa*(ga_eff-1.0);
        flow_float z = ua*ua+va*va+wa*wa - P_ro/P_roe;
        flow_float fai = P_ro-ca*ca;

        R[0][0] = 1.0;      R[0][1] = nx;               R[0][2] = ny;              R[0][3] = nz;              R[0][4] = 1.0;
        R[1][0] = Ha-ca*Ua; R[1][1] = z*nx+va*nz-wa*ny; R[1][2] = z*ny+wa*nx-ua*nz;R[1][3] = z*nz+ua*ny-va*nx;R[1][4] = Ha+ca*Ua;
        R[2][0] = ua-ca*nx; R[2][1] = ua*nx;            R[2][2] = ua*ny-nz        ;R[2][3] = ua*nz+ny        ;R[2][4] = ua+ca*nx;
        R[3][0] = va-ca*ny; R[3][1] = va*nx+nz;         R[3][2] = va*ny           ;R[3][3] = va*nz-nx        ;R[3][4] = va+ca*ny;
        R[4][0] = wa-ca*nz; R[4][1] = wa*nx-ny;         R[4][2] = wa*ny+nx        ;R[4][3] = wa*nz;           R[4][4] = wa+ca*nz;

        L[0][0] = 0.5*(P_ro+ca*Ua);            L[0][1] = P_roe*0.5; L[0][2] = (P_rou-ca*nx)*0.5 ; L[0][3] = (P_rov-ca*ny)*0.5 ; L[0][4] =  (P_row-ca*nz)*0.5;
        L[1][0] = -fai*nx-ca*ca*(va*nz-wa*ny); L[1][1] = -P_roe*nx; L[1][2] = -P_rou*nx         ; L[1][3] = -P_rov*nx+ca*ca*nz; L[1][4] = -P_row*nx-ca*ca*ny;
        L[2][0] = -fai*ny-ca*ca*(wa*nx-ua*nz); L[2][1] = -P_roe*ny; L[2][2] = -P_rou*ny-ca*ca*nz; L[2][3] = -P_rov*ny         ; L[2][4] = -P_row*ny+ca*ca*nx;
        L[3][0] = -fai*nz-ca*ca*(ua*ny-va*nx); L[3][1] = -P_roe*nz; L[3][2] = -P_rou*nz+ca*ca*ny; L[3][3] = -P_rov*nz-ca*ca*nx; L[3][4] = -P_row*nz;
        L[4][0] = 0.5*(P_ro-ca*Ua);
        L[4][1] = 0.5*P_roe;
        L[4][2] = 0.5*(P_rou+ca*nx);
        L[4][3] = 0.5*(P_rov+ca*ny);
        L[4][4] = 0.5*(P_row+ca*nz);

        delQ[0] = ro_R - ro_L;
        delQ[1] = (ro_R*Ht_R-P_R) - (ro_L*Ht_L-P_L);
        delQ[2] = ro_R*Ux_R - ro_L*Ux_L;
        delQ[3] = ro_R*Uy_R - ro_L*Uy_L;
        delQ[4] = ro_R*Uz_R - ro_L*Uz_L;

        for (int i=0; i<5; i++) {
            for (int j=0; j<5; j++) {
                L[i][j] = L[i][j]/(ca*ca);
            }
        }

        lam[0] = abs(Ua - ca);
        lam[1] = abs(Ua     );
        lam[2] = abs(Ua     );
        lam[3] = abs(Ua     );
        lam[4] = abs(Ua + ca);

        // van leer's entropy fix
        flow_float eta_vl ;
        if (false) {
            flow_float eta_vl = max(dU/sss-(ca_R-ca_L),0.0);
            if (lam[0] < 2*eta_vl) lam[0] = lam[0]*lam[0]/(4*eta_vl) + eta_vl;

            eta_vl = max(dU/sss,0.0);
            if (lam[1] < 2*eta_vl) {
                lam[1] = lam[1]*lam[1]/(4*eta_vl) + eta_vl;
                lam[2] = lam[1];
                lam[3] = lam[1];
            }
            eta_vl = max(dU/sss+(ca_R-ca_L),0.0);
            if (lam[4] < 2*eta_vl) lam[4] = lam[4]*lam[4]/(4*eta_vl) + eta_vl;
        }
        // harten entropy correction
        if (false) {
            eta_vl = 0.1; // 0.05~0.15
            if (lam[0] <= eta_vl) lam[0] = (lam[0]*lam[0] + eta_vl*eta_vl)/(2*eta_vl);

            if (lam[1] < eta_vl) {
                lam[1] = (lam[1]*lam[1] + eta_vl*eta_vl)/(2*eta_vl);
                lam[2] = lam[1];
                lam[3] = lam[1];
            }
            if (lam[4] <= eta_vl) lam[4] = (lam[4]*lam[4] + eta_vl*eta_vl)/(2*eta_vl);
        }
 

        // harten entropy correction
        if (true) {
            eta_vl = 0.1*(abs(Ua)/ca+1.0); 

            if (lam[0] <= eta_vl) lam[0] = (lam[0]*lam[0] + eta_vl*eta_vl)/(2*eta_vl);

            if (lam[1] < eta_vl) {
                lam[1] = (lam[1]*lam[1] + eta_vl*eta_vl)/(2*eta_vl);
                lam[2] = lam[1];
                lam[3] = lam[1];
            }
            if (lam[4] <= eta_vl) lam[4] = (lam[4]*lam[4] + eta_vl*eta_vl)/(2*eta_vl);
        }
 
        // SU2 同形の固有値下限 (opt-in)。**全 5 固有値**に coeff*(|Ua|+ca) の床を課す。
        // 壁近傍は |Ua|->0 なので、接触波 (lam[1..3] = |Ua|) がここで初めて有限の散逸を得る。
        // 既定 0 でビット不変。SU2: roe.cpp:200 `max(fabs(Lambda), coeff*MaxLambda)`。
        if (entropyFixCoeff > (flow_float)0.0) {
            const flow_float floorLam = entropyFixCoeff*(abs(Ua) + ca);
            for (int i = 0; i < 5; i++) lam[i] = max(lam[i], floorLam);
        }

        // entropy fix
        //flow_float h_fix = 0.05*(abs(Ua)+ca);

        //flow_float l_temp = lam[0];
        //lam[0] = min(l_temp, 0.5*(l_temp*l_temp/h_fix + h_fix));

        //l_temp = lam[1];
        //lam[1] = min(l_temp, 0.5*(l_temp*l_temp/h_fix + h_fix));

        //l_temp = lam[2];
        //lam[2] = min(l_temp, 0.5*(l_temp*l_temp/h_fix + h_fix));

        //l_temp = lam[3];
        //lam[3] = min(l_temp, 0.5*(l_temp*l_temp/h_fix + h_fix));

        //l_temp = lam[4];
        //lam[4] = min(l_temp, 0.5*(l_temp*l_temp/h_fix + h_fix));
        dW[0] = 0.0;
        dW[1] = 0.0;
        dW[2] = 0.0;
        dW[3] = 0.0;
        dW[4] = 0.0;

        for (int i=0; i<5; i++) {
            for (int j=0; j<5; j++) {
                dW[i] += L[i][j]*delQ[j];
            }
        }

        difQ1[0] = 0.0;
        difQ1[1] = 0.0;
        difQ1[2] = 0.0;
        difQ1[3] = 0.0;
        difQ1[4] = 0.0;

        for (int i=0; i<5; i++) {
            difQ1[i] = lam[i]*dW[i];
        }

        difQ2[0] = 0.0;
        difQ2[1] = 0.0;
        difQ2[2] = 0.0;
        difQ2[3] = 0.0;
        difQ2[4] = 0.0;

        for (int i=0; i<5; i++) {
            for (int j=0; j<5; j++) {
                difQ2[i] += R[i][j]*difQ1[j];
            }
        }

        flow_float mdot = 0.5*(ro_L*(Ux_L*nx+Uy_L*ny+Uz_L*nz) + ro_R*(Ux_R*nx+Uy_R*ny+Uz_R*nz))*sss ;
        massflux[ip] = mdot;

        flow_float res_ro_temp   = mdot                                        -0.5*difQ2[0]*sss;
        flow_float res_roe_temp  = mdot*0.5*(Ht_L + Ht_R)                      -0.5*difQ2[1]*sss;
        flow_float res_roUx_temp = mdot*0.5*(Ux_L + Ux_R) +0.5*(P_L + P_R)*sxx -0.5*difQ2[2]*sss;
        flow_float res_roUy_temp = mdot*0.5*(Uy_L + Uy_R) +0.5*(P_L + P_R)*syy -0.5*difQ2[3]*sss;
        flow_float res_roUz_temp = mdot*0.5*(Uz_L + Uz_R) +0.5*(P_L + P_R)*szz -0.5*difQ2[4]*sss;



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
