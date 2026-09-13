#include "convection/convectiveFlux_d.cuh"
#include "cuda_forge/wmlesWallModel_d.cuh"   // wmlesActiveForBcond / wmlesNodeActive (WMLES ゲート)
#include <cstdlib>
#include <cstdio>
#include <vector>

__global__ void viscousFlux_d
( 
 // mesh structure
 geom_int nCells,
 geom_int nPlanes, geom_int nNormalPlanes, geom_int* plane_cells,  
 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 flow_float mu ,  flow_float Prt , flow_float* cp , flow_float* thermCond,
 flow_float* vis_lam   , flow_float* vis_turb  ,

 // variables
//flow_float* convx , flow_float* convy , flow_float* convz,
// flow_float* diffx , flow_float* diffy , flow_float* diffz,
 flow_float* ro   ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* Ps  ,
 flow_float* Ht  ,
 flow_float* sonic,
 flow_float* Ts  ,

 flow_float* res_ro   ,
 flow_float* res_roUx  ,
 flow_float* res_roUy  ,
 flow_float* res_roUz  ,
 flow_float* res_roe   ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,

 // 軸対称: 完全な発散 ∇·u = ∂xUx+∂yUy+u_r/r (axisym_divU)。-2/3 μ(∇·u) の体積粘性項を
 // τθθ ソース (axisymmetricSource_d) と整合させるため planar 面でもフープ項込みの divu を使う。
 int isAxisymmetric, flow_float* axisym_divU,

 // SST automatic wall treatment (node, wallTreatmentSST==1) / WMLES (node): 壁ノードに格納した
 // τ_w=ρu_τ² (>0)。片端だけ壁ノードの内部双対面 (W-I) で接線粘性応力を τ_w に再スケールする
 // (SU2 AddTauWall)。cell/非壁関数では Tau_Wall≡-1 で再スケール無効 (ビット不変)。nullptr 可 (常に無効)。
 flow_float* Tau_Wall,
 // 診断 (§4.2): W-I 面で実際に残差へ加えた接線力/法線力と、再スケール前の解像接線力を壁ノードへ集計。
 // nullptr 可 (診断オフ)。
 flow_float* wi_ftan, flow_float* wi_fnrm, flow_float* wi_fnrm_abs, flow_float* wi_ftan_res,

 // WMLES 等温壁 (node): 壁ノードに格納した q_w [W/m²] (壁→流体正, 非対象は -1)。片端だけ壁ノードの
 // W-I 面で解像伝導熱流束を q_w·S に置換する (AddQWall, methods/turbulence §10.4)。nullptr で無効。
 flow_float* Qw_Wall,
 // SST 断熱壁 SU2 式熱結合 (experimental, sstThermalWallFunction==2, plan
 // turbulence-sst-su2-taw-coupling): 勾配計算後の primitive 温度 overlay (壁ノード=Taw, 非対象 -1)。
 // overlay 端点を持つ内部辺のみ、熱流束を SU2 corrected-gradient (算術平均 k_eff)、粘性仕事の
 // 面速度を算術平均へ切替える。nullptr (mode 0/1) では従来経路ビット不変。
 //
 // 【注意 (2026-08-11 凍結場収支診断)】旧 compact 端点置換 (over-relaxed 式へ Taw を単純代入,
 // 壁 μt 未修正) は壁半 CV に −12 W 級の無補償ドレインを作り EOS 床まで異常冷却して発散した
 // (run_0053/0057)。mode 2 は SU2 と同じく (i) corrected-gradient + 算術平均, (ii) 壁ノード k_wf
 // Dirichlet による壁法則整合 μt (ransWallFunction) を併せて初めてドレインが SU2 同等 (−2 W 級,
 // 過渡分) に閉じる。片方だけの適用は禁止。plans/active/turbulence-sst-su2-taw-coupling.md 参照。
 flow_float* Taw_Ov,
 // 【実験専用・恒久実装禁止】mode 2 overlay 辺の k_eff で「overlay 端点 (壁) 側 μt」に掛ける
 // スケール (env FORGE_TAW_WALL_MUT_SCALE, 既定 1.0=無効)。壁 μt 過大 (SU2 比 ~5x) が
 // ドレインの主因という §11 診断の動的確証用。恒久対策は壁 μt モデルの plan 決着後に実装。
 flow_float tawWallMutScale,
 // SST 断熱壁 defect-flux 閉包 (experimental, sstThermalWallFunction==3, theory §6.5(f) mode 3):
 // 壁ノードに格納した H⃗=H_T·m̂ (ransWallFunction, 非壁 0)。片端のみ非零の W-I 辺で、
 // エネルギー行の (伝導 + τ·u 仕事) を defect 流束 F=(S_out·H⃗)(Taw−T_W) の ± ペアに置換する
 // (Couette 恒等式 q+τu=0: T_W=Taw で全流束ゼロ、∂F/∂T_W<0 の負帰還で T_W→Taw)。
 // 運動量行は不変。nullptr (mode≠3) で従来経路ビット不変。Taw_diag は defect の目標温度。
 flow_float* Taw_HTnx, flow_float* Taw_HTny, flow_float* Taw_HTnz, flow_float* Taw_diag
,
 flow_float* kturb, int isoStress,   // sstIsotropicStress: -(2/3) rho k delta_ij (nullptr/0 で無効)
 // sstEnergyIncludesK (plan turbulence-sst-energy-includes-k §4.2): E_t のエネルギー流束に k 拡散
 //   (μ + σ_k μt)(∂k/∂n) S を足す。離散化は scalarTransport の k 拡散 (法線 over-relaxed 項のみ、相対ゼロ割ガード、
 //   σ_k は sstF1 ブレンド) と同形。nullptr/0 で無効。
 flow_float* dKdx_e, flow_float* dKdy_e, flow_float* dKdz_e, flow_float* sstF1_e, int sigmaBlend_e, int energyK
)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;
    (void)dKdx_e; (void)dKdy_e; (void)dKdz_e;   // 勾配クロス項は k 式と同じく使わない (法線項のみ)


    //if (ip < nPlanes) { 
    if (ip < nNormalPlanes) {

        geom_int  ic0 = plane_cells[2*ip+0];
        geom_int  ic1 = plane_cells[2*ip+1];

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
        flow_float dcc   = sqrt(dcc_x*dcc_x +dcc_y*dcc_y +dcc_z*dcc_z) ;

        flow_float Uxf = f*Ux[ic0] + (1.0f-f)*Ux[ic1];
        flow_float Uyf = f*Uy[ic0] + (1.0f-f)*Uy[ic1];
        flow_float Uzf = f*Uz[ic0] + (1.0f-f)*Uz[ic1];

        flow_float dUxdxf = f*dUxdx[ic0] + (1.0f-f)*dUxdx[ic1];
        flow_float dUxdyf = f*dUxdy[ic0] + (1.0f-f)*dUxdy[ic1];
        flow_float dUxdzf = f*dUxdz[ic0] + (1.0f-f)*dUxdz[ic1];

        flow_float dUydxf = f*dUydx[ic0] + (1.0f-f)*dUydx[ic1];
        flow_float dUydyf = f*dUydy[ic0] + (1.0f-f)*dUydy[ic1];
        flow_float dUydzf = f*dUydz[ic0] + (1.0f-f)*dUydz[ic1];

        flow_float dUzdxf = f*dUzdx[ic0] + (1.0f-f)*dUzdx[ic1];
        flow_float dUzdyf = f*dUzdy[ic0] + (1.0f-f)*dUzdy[ic1];
        flow_float dUzdzf = f*dUzdz[ic0] + (1.0f-f)*dUzdz[ic1];

        flow_float dTdxf = f*dTdx[ic0] + (1.0f-f)*dTdx[ic1];
        flow_float dTdyf = f*dTdy[ic0] + (1.0f-f)*dTdy[ic1];
        flow_float dTdzf = f*dTdz[ic0] + (1.0f-f)*dTdz[ic1];
#ifdef VISC_DIAG_NOGRAD
        // 診断: 面勾配 (Green-Gauss) 由来の項 (転置・非直交補正) を 0 にし、over-relaxed 法線 Laplacian のみ残す。
        dUxdxf=dUxdyf=dUxdzf=0; dUydxf=dUydyf=dUydzf=0; dUzdxf=dUzdyf=dUzdzf=0; dTdxf=dTdyf=dTdzf=0;
#endif

        // implicit 側（timeIntegration_d.cu）と同じ fabs()+下限で分母を保護する
        flow_float D_safe  = max(fabs(dcc_x*sxx +dcc_y*syy +dcc_z*szz), (flow_float)1.0e-30);
        flow_float delta   = dcc*sss*sss/D_safe; // over-relaxed
        flow_float delta_x = dcc_x*sss*sss/D_safe;
        flow_float delta_y = dcc_y*sss*sss/D_safe;
        flow_float delta_z = dcc_z*sss*sss/D_safe;
        flow_float k_x = sxx - delta_x;
        flow_float k_y = syy - delta_y;
        flow_float k_z = szz - delta_z;
        // 軸対称は完全発散 (u_r/r 込み) を面補間。デカルトは従来どおり ∂xUx+∂yUy+∂zUz。
        flow_float divu = (isAxisymmetric == 1)
            ? (f*axisym_divU[ic0] + (1.0f-f)*axisym_divU[ic1])
            : (dUxdxf+dUydyf+dUzdzf);

        flow_float v_lam  = f*vis_lam [ic0] + (1.0f-f)*vis_lam [ic1] ;
        flow_float v_turb = f*vis_turb[ic0] + (1.0f-f)*vis_turb[ic1] ;
        flow_float mu_total = v_lam + v_turb;

        // 完全な Newton 応力 tau_ij S_j = mu(du_i/dx_j + du_j/dx_i)S_j - (2/3)mu divu S_i。
        // 第1項 (Laplacian, mu grad(u_i).S) は over-relaxed: 法線スカラー delta + 同成分勾配.k。
        // 第2項 (転置, mu du_j/dx_i S_j) は面平均勾配にフル S を内積。第3項 (発散) は成分 s**。
        flow_float tau_x = mu_total*((Ux[ic1] -Ux[ic0])/dcc)*delta;
        tau_x += mu_total*(dUxdxf*k_x +dUxdyf*k_y +dUxdzf*k_z);
        tau_x += mu_total*(dUxdxf*sxx +dUydxf*syy +dUzdxf*szz);
        tau_x += -mu_total*2.0f/3.0f*(divu)*sxx;
        if (isoStress != 0 && kturb != nullptr) tau_x += -(2.0f/3.0f)*(f*ro[ic0]+(1.0f-f)*ro[ic1])*(f*kturb[ic0]+(1.0f-f)*kturb[ic1])*sxx;

        flow_float tau_y = mu_total*((Uy[ic1] -Uy[ic0])/dcc)*delta;
        tau_y += mu_total*(dUydxf*k_x +dUydyf*k_y +dUydzf*k_z);
        tau_y += mu_total*(dUxdyf*sxx +dUydyf*syy +dUzdyf*szz);
        tau_y += -mu_total*2.0f/3.0f*(divu)*syy;
        if (isoStress != 0 && kturb != nullptr) tau_y += -(2.0f/3.0f)*(f*ro[ic0]+(1.0f-f)*ro[ic1])*(f*kturb[ic0]+(1.0f-f)*kturb[ic1])*syy;

        flow_float tau_z = mu_total*((Uz[ic1] -Uz[ic0])/dcc)*delta;
        tau_z += mu_total*(dUzdxf*k_x +dUzdyf*k_y +dUzdzf*k_z);
        tau_z += mu_total*(dUxdzf*sxx +dUydzf*syy +dUzdzf*szz);
        tau_z += -mu_total*2.0f/3.0f*(divu)*szz;
        if (isoStress != 0 && kturb != nullptr) tau_z += -(2.0f/3.0f)*(f*ro[ic0]+(1.0f-f)*ro[ic1])*(f*kturb[ic0]+(1.0f-f)*kturb[ic1])*szz;

        // SST node 壁関数 (SU2 AddTauWall): 片端のみ壁ノードの内部双対面 (W-I) で、解像した粘性 traction
        // の接線成分をモデル τ_w に再スケールする。粗い y+ メッシュでは生の解像勾配が τ_w を過小評価する
        // ため、壁関数の τ_w を内部ノード I の運動量・エネルギーに正しく伝える (壁ノード残差は Dirichlet で
        // 別途ゼロ化)。cell/非壁関数では Tau_Wall≡-1 で無効 (ビット不変)。
        if (Tau_Wall != nullptr) {
            const flow_float tw0 = Tau_Wall[ic0];
            const flow_float tw1 = Tau_Wall[ic1];
            const bool w0 = (tw0 > (flow_float)0.0);
            const bool w1 = (tw1 > (flow_float)0.0);
            if (w0 != w1) {                                   // 片端のみ壁ノード (xor)
                const flow_float tauw  = w0 ? tw0 : tw1;      // 壁ノード側の τ_w [Pa]
                const flow_float inv_s = (flow_float)1.0 / max(sss, (flow_float)1.0e-30);
                const flow_float nhx = sxx*inv_s, nhy = syy*inv_s, nhz = szz*inv_s;
                const flow_float Tn  = tau_x*nhx + tau_y*nhy + tau_z*nhz;   // 法線 traction
                const flow_float ttx = tau_x - Tn*nhx;
                const flow_float tty = tau_y - Tn*nhy;
                const flow_float ttz = tau_z - Tn*nhz;
                const flow_float tmag = sqrt(ttx*ttx + tty*tty + ttz*ttz); // 接線 traction の大きさ [N]
                const flow_float scale = (tauw*sss) / max(tmag, (flow_float)1.0e-30); // |接線|→τ_w·面積
                tau_x *= scale; tau_y *= scale; tau_z *= scale;
                // 診断 (§4.2): 実際に残差へ入る力を壁ノードへ集計する。
                if (wi_ftan != nullptr) {
                    const geom_int icW = w0 ? ic0 : ic1;
                    const flow_float Tn2 = tau_x*nhx + tau_y*nhy + tau_z*nhz;   // 再スケール後の法線
                    const flow_float tx2 = tau_x - Tn2*nhx;
                    const flow_float ty2 = tau_y - Tn2*nhy;
                    const flow_float tz2 = tau_z - Tn2*nhz;
                    atomicAdd(&wi_ftan[icW],     sqrt(tx2*tx2 + ty2*ty2 + tz2*tz2));
                    atomicAdd(&wi_fnrm[icW],     Tn2);          // 符号付き (相殺しうる)
                    atomicAdd(&wi_fnrm_abs[icW], fabs(Tn2));    // 絶対値 (相殺なし)
                    atomicAdd(&wi_ftan_res[icW], tmag);
                }
            }
        }

        // 乱流熱伝導: 有効熱伝導率 k_eff = k_lam + cp*mu_turb/Pr_t。
        // 応力(摩擦発熱)が mu_total=v_lam+v_turb を使うのに対し従来は層流 k のみで
        // 熱を運んでいたため、乱流境界層で散逸熱が逃げ場を失い T が全温を超えて
        // overshoot していた。RANS のエネルギー保存には乱流熱伝導が必須。
        // cp は thermally-perfect の温度依存を反映するためセル配列を面平均で使う。
        flow_float cp_face = f*cp[ic0] + (1.0f-f)*cp[ic1];
        flow_float tc_face = f*thermCond[ic0] + (1.0f-f)*thermCond[ic1];
        tc_face += cp_face*v_turb/Prt;
        // W-I 内部熱拡散は既定で DOF 状態 (Ts) と DOF 勾配で評価する。モデル温度の単純 compact
        // 代入は禁止 (上の Taw_Ov コメント参照)。例外は mode 2 の SU2 corrected-gradient (下)。
        flow_float heatflux = tc_face*((Ts[ic1] -Ts[ic0])/dcc)*delta;
        heatflux += tc_face*(dTdxf*k_x +dTdyf*k_y +dTdzf*k_z);

        // SST 断熱壁 SU2 式熱結合 (mode 2): overlay 端点 (Taw) を持つ内部辺は SU2 corrected-gradient
        //   g_corr = ḡ + (ΔT_flux − ḡ·d)·d/|d|²,  q = k_eff·(g_corr·S)
        // で熱流束を置換する (ḡ=端点 DOF 勾配の算術平均, ΔT_flux=overlay 済み端点温度差,
        // k_eff は SU2 と同じ算術平均)。粘性仕事の面速度も算術平均へ (f 補間は壁側重み ~2/3 で
        // 内点速度を過小評価し、壁半 CV の仕事項が SU2 比 ~0.6 倍になる — 凍結場収支診断で確認)。
        if (Taw_Ov != nullptr) {
            const flow_float ov0 = Taw_Ov[ic0];
            const flow_float ov1 = Taw_Ov[ic1];
            if (ov0 > (flow_float)0.0 || ov1 > (flow_float)0.0) {
                const flow_float T0 = (ov0 > (flow_float)0.0) ? ov0 : Ts[ic0];
                const flow_float T1 = (ov1 > (flow_float)0.0) ? ov1 : Ts[ic1];
                const flow_float mut0 = (ov0 > (flow_float)0.0) ? tawWallMutScale*vis_turb[ic0] : vis_turb[ic0];
                const flow_float mut1 = (ov1 > (flow_float)0.0) ? tawWallMutScale*vis_turb[ic1] : vis_turb[ic1];
                const flow_float tc_arith = (flow_float)0.5*(thermCond[ic0]+thermCond[ic1])
                    + (flow_float)0.5*(cp[ic0]+cp[ic1])
                      * (flow_float)0.5*(mut0+mut1) / Prt;
                const flow_float gmx = (flow_float)0.5*(dTdx[ic0]+dTdx[ic1]);
                const flow_float gmy = (flow_float)0.5*(dTdy[ic0]+dTdy[ic1]);
                const flow_float gmz = (flow_float)0.5*(dTdz[ic0]+dTdz[ic1]);
                const flow_float dd = max(dcc_x*dcc_x + dcc_y*dcc_y + dcc_z*dcc_z, (flow_float)1.0e-30);
                const flow_float corr = ((T1 - T0) - (gmx*dcc_x + gmy*dcc_y + gmz*dcc_z)) / dd;
                heatflux = tc_arith*((gmx + corr*dcc_x)*sxx
                                    +(gmy + corr*dcc_y)*syy
                                    +(gmz + corr*dcc_z)*szz);
                Uxf = (flow_float)0.5*(Ux[ic0] + Ux[ic1]);
                Uyf = (flow_float)0.5*(Uy[ic0] + Uy[ic1]);
                Uzf = (flow_float)0.5*(Uz[ic0] + Uz[ic1]);
            }
        }

        // WMLES 等温壁 (node, AddQWall): 片端のみ壁ノードの W-I 面で、解像伝導熱流束を壁モデルの
        // q_w·S に置換する (AddTauWall の熱版)。粗い y+ メッシュでは解像 ∇T が q_w を過小評価する
        // ため、Kader 由来の q_w を内部ノード I のエネルギーへ正しく伝える (壁ノード側 res_roe は
        // Qw_Wall マーカで別途ゼロ化)。q_w>0 = 壁→流体。非対象面 (Qw_Wall=-1) は解像値のまま。
        if (Qw_Wall != nullptr) {
            const flow_float q0 = Qw_Wall[ic0];
            const flow_float q1 = Qw_Wall[ic1];
            const bool wq0 = (q0 > (flow_float)-0.5f);
            const bool wq1 = (q1 > (flow_float)-0.5f);
            if (wq0 != wq1) {                          // 片端のみ壁ノード (xor)
                const flow_float qw = wq0 ? q0 : q1;
                // 符号: res_roe[ic0]+=heatflux / res_roe[ic1]-=heatflux。内部側が +q_w·S を受ける向き。
                heatflux = (wq1 ? qw : -qw) * sss;
            }
        }

        flow_float res_ro_temp   = 0.0f;
        flow_float res_roUx_temp = tau_x;
        flow_float res_roUy_temp = tau_y;
        flow_float res_roUz_temp = tau_z;
        flow_float res_roe_temp  = tau_x*Uxf +tau_y*Uyf +tau_z*Uzf;
        res_roe_temp += heatflux;

        // sstEnergyIncludesK: k 拡散のエネルギー流束 (k 式の拡散と同形: 法線項のみ・相対ガード・σ_k ブレンド)
        if (energyK != 0 && kturb != nullptr) {
            const flow_float denom_k = dcc_x*sxx + dcc_y*syy + dcc_z*szz;
            const flow_float floor_k = (flow_float)1.0e-6 * dcc * sss;
            const flow_float safe_k  = (fabs(denom_k) < floor_k) ? ((denom_k >= (flow_float)0.0) ? floor_k : -floor_k) : denom_k;
            const flow_float delta_k = dcc*sss*sss/safe_k;
            const flow_float F1a = (sigmaBlend_e != 0 && sstF1_e != nullptr) ? sstF1_e[(ic0 < nCells) ? ic0 : ic1] : (flow_float)1.0;
            const flow_float F1b = (sigmaBlend_e != 0 && sstF1_e != nullptr) ? sstF1_e[(ic1 < nCells) ? ic1 : ic0] : (flow_float)1.0;
            const flow_float sig0 = F1a*(flow_float)0.85 + ((flow_float)1.0 - F1a)*(flow_float)1.0;
            const flow_float sig1 = F1b*(flow_float)0.85 + ((flow_float)1.0 - F1b)*(flow_float)1.0;
            const flow_float mu0 = vis_lam[ic0] + sig0*max(vis_turb[ic0], (flow_float)0.0);
            const flow_float mu1 = vis_lam[ic1] + sig1*max(vis_turb[ic1], (flow_float)0.0);
            const flow_float mu_k = f*mu0 + ((flow_float)1.0 - f)*mu1;
            res_roe_temp += mu_k * ((kturb[ic1] - kturb[ic0])/dcc) * delta_k;
        }

        // SST 断熱壁 defect-flux 閉包 (mode 3): 片端のみ壁 (|H⃗|>0) の W-I 辺で、エネルギー行の
        // (伝導+仕事) を F=(S_out·H⃗)(Taw−T_W) に置換 (±保存)。運動量 (AddTauWall 済) は不変。
        // S_out·m̂ の和は双対 CV 面閉性から壁面積に一致し、多重 W-I 辺へ自動分配される。
        if (Taw_HTnx != nullptr) {
            const flow_float h0 = Taw_HTnx[ic0]*Taw_HTnx[ic0] + Taw_HTny[ic0]*Taw_HTny[ic0]
                                + Taw_HTnz[ic0]*Taw_HTnz[ic0];
            const flow_float h1 = Taw_HTnx[ic1]*Taw_HTnx[ic1] + Taw_HTny[ic1]*Taw_HTny[ic1]
                                + Taw_HTnz[ic1]*Taw_HTnz[ic1];
            const bool w0 = (h0 > (flow_float)0.0);
            const bool w1 = (h1 > (flow_float)0.0);
            if (w0 != w1) {                                    // 片端のみ壁ノード (xor)
                const geom_int  icW  = w0 ? ic0 : ic1;
                const flow_float sgnW = w0 ? (flow_float)1.0 : (flow_float)-1.0f; // S を W から離れる向きへ
                const flow_float a = sgnW * (sxx*Taw_HTnx[icW] + syy*Taw_HTny[icW] + szz*Taw_HTnz[icW]);
                const flow_float Fin = a * (Taw_diag[icW] - Ts[icW]);  // 壁ノードへ入る向き正 [W]
                // res_roe[ic0] += temp / res_roe[ic1] -= temp の符号系で W が +Fin を受ける形に置換
                res_roe_temp = sgnW * Fin;
            }
        }

        atomicAdd(&res_ro[ic0]  , res_ro_temp);
        atomicAdd(&res_roUx[ic0], res_roUx_temp);
        atomicAdd(&res_roUy[ic0], res_roUy_temp);
        atomicAdd(&res_roUz[ic0], res_roUz_temp);
        atomicAdd(&res_roe[ic0] , res_roe_temp);

        atomicAdd(&res_ro[ic1]  , -res_ro_temp);
        atomicAdd(&res_roUx[ic1], -res_roUx_temp);
        atomicAdd(&res_roUy[ic1], -res_roUy_temp);
        atomicAdd(&res_roUz[ic1], -res_roUz_temp);
        atomicAdd(&res_roe[ic1] , -res_roe_temp);
    }

    __syncthreads();
}


__global__ void viscousFlux_wall_d
( 
  // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  

 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 flow_float mu ,  flow_float Prt , flow_float* cp , flow_float* thermCond,
 flow_float* vis_lam   , flow_float* vis_turb  ,

 // variables
//flow_float* convx , flow_float* convy , flow_float* convz,
// flow_float* diffx , flow_float* diffy , flow_float* diffz,
 flow_float* ro   ,
 flow_float* roUx  ,
 flow_float* roUy  ,
 flow_float* roUz  ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* Ps  ,
 flow_float* Ht  ,
 flow_float* sonic,
 flow_float* Ts  ,

 flow_float* res_ro   ,
 flow_float* res_roUx  ,
 flow_float* res_roUy  ,
 flow_float* res_roUz  ,
 flow_float* res_roe   ,

 flow_float* dUxdx  , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx  , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx  , flow_float* dUzdy , flow_float* dUzdz,
 flow_float* dTdx   , flow_float* dTdy  , flow_float* dTdz,

 // bvar
 flow_float* ypls_b ,
 flow_float* twall_x_b ,
 flow_float* twall_y_b ,
 flow_float* twall_z_b ,
 flow_float* Ux_b ,
 flow_float* Uy_b ,
 flow_float* Uz_b ,
 flow_float* Ts_b ,
//flow_float* sx_b ,
 //flow_float* sy_b ,
 //flow_float* sz_b

 // 軸対称: 完全発散 (u_r/r 込み) を体積粘性項に使う (内部面と同趣旨)
 int isAxisymmetric, flow_float* axisym_divU,

 // 壁処理: 0=なし (解像壁) / 1=SST automatic wall treatment (methods/turbulence §6.5 (c)):
 //   接線せん断を modeled τ_w=ρu_τ² に置換。utau_b は事前算出済み u_τ。
 // 2=WMLES 代数壁応力モデル (§10.4): wmlesWallModel が書いた twall_*_b (流束/面積) と
 //   qwall_b (q_w, 壁→流体正) で粘性流束を置換。
 int wallTreatment, flow_float* utau_b, flow_float* qwall_b,
 // node-centered (1): 壁法線 Laplacian/熱流束を ghost+dcc でなくセル勾配 ∇φ·S (bvar 壁閉包) で評価。
 // 壁ノードが壁面に乗り dcc≈0 に退化するのを回避する。cell (0) は従来の (φ[ig]-φ[ic])/dcc。
 int isNode,
 // 断熱壁 (kind: wall) = 1: 壁面伝導熱流束を厳密に 0 にする。cell はミラーゴーストで従来も
 // ビット厳密 0 だが、node は壁ノードのセル勾配 ∇T·S が法線成分を持ち断熱が破れていた
 // (壁エントロピー市松の直接容疑, plan boundary-node-nozzle-wall-outlet-stability §2.8)。
 // 等温壁 (wall_isothermal) = 0 で従来どおり解像熱流束を入れる。WMLES (==2) は qwall が担う。
 int adiabaticWall,
 // SST エネルギー壁関数 (§6.5(g), sstEnergyWallFunction==1 × wall_isothermal): 壁面熱流束を
 // ransWallFunction が書いた qwall_b (Kader q_w) に置換する (運動量は wallTreatment==1 のまま)。
 // node では壁ノード res_roe が Dirichlet で 0 化されるため実効は cell のみ (書いても無害)。
 int sstEnergyWf
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

        flow_float ccx_0 = ccx[ic];
        flow_float ccy_0 = ccy[ic];
        flow_float ccz_0 = ccz[ic];

        flow_float pcx_1 = pcx[ip];
        flow_float pcy_1 = pcy[ip];
        flow_float pcz_1 = pcz[ip];

        flow_float ccx_1 = ccx[ig];
        flow_float ccy_1 = ccy[ig];
        flow_float ccz_1 = ccz[ig];

        //flow_float dcc_x = pcx_1 - ccx_0;
        //flow_float dcc_y = pcy_1 - ccy_0;
        //flow_float dcc_z = pcz_1 - ccz_0;

        flow_float dcc_x = ccx_1 - ccx_0;
        flow_float dcc_y = ccy_1 - ccy_0;
        flow_float dcc_z = ccz_1 - ccz_0;
 
        flow_float dcc   = sqrt(dcc_x*dcc_x +dcc_y*dcc_y +dcc_z*dcc_z) ;

        flow_float dUxdxf = dUxdx[ic] ;
        flow_float dUxdyf = dUxdy[ic] ;
        flow_float dUxdzf = dUxdz[ic] ;

        flow_float dUydxf = dUydx[ic] ;
        flow_float dUydyf = dUydy[ic] ;
        flow_float dUydzf = dUydz[ic] ;

        flow_float dUzdxf = dUzdx[ic] ;
        flow_float dUzdyf = dUzdy[ic] ;
        flow_float dUzdzf = dUzdz[ic] ;

        flow_float v_turb = vis_turb[ic] ;
        flow_float mu_total = vis_lam[ic] + v_turb;

        // 軸対称は完全発散 (u_r/r 込み, セル値)。壁面はセル中心勾配を直接使うため cell の axisym_divU。
        flow_float divu = (isAxisymmetric == 1) ? axisym_divU[ic] : (dUxdxf+dUydyf+dUzdzf);

        // ミラーゴースト (u^g=-u^c) では d∥S なので over-relaxed の delta=sss、k=0。
        // 法線項: 法線勾配 (u^g-u^c)/dcc に面積 sss。転置項: セル中心勾配にフル S。発散項: 成分 s**。
        // 注: dcc フロア等で物理 flux をマスクしない方針 (退化 dcc はメッシュ/境界クロージャの欠陥であり、
        // 直すべきはメッシュとゴーストレス弱形式。診断は viscousWallDiag_d で別途行う)。
        // 法線 Laplacian 項: node は ∇U·S (セル勾配, bvar 壁閉包) で退化 dcc を回避、cell は従来の法線差分。
        flow_float tau_x = (isNode != 0) ? mu_total*(dUxdxf*sxx +dUxdyf*syy +dUxdzf*szz)
                                         : mu_total*((Ux[ig] - Ux[ic])/dcc)*sss;
        tau_x += mu_total*(dUxdxf*sxx +dUydxf*syy +dUzdxf*szz);
        tau_x += -mu_total*2.0f/3.0f*(divu)*sxx;

        flow_float tau_y = (isNode != 0) ? mu_total*(dUydxf*sxx +dUydyf*syy +dUydzf*szz)
                                         : mu_total*((Uy[ig] - Uy[ic])/dcc)*sss;
        tau_y += mu_total*(dUxdyf*sxx +dUydyf*syy +dUzdyf*szz);
        tau_y += -mu_total*2.0f/3.0f*(divu)*syy;

        flow_float tau_z = (isNode != 0) ? mu_total*(dUzdxf*sxx +dUzdyf*syy +dUzdzf*szz)
                                         : mu_total*((Uz[ig] - Uz[ic])/dcc)*sss;
        tau_z += mu_total*(dUxdzf*sxx +dUydzf*syy +dUzdzf*szz);
        tau_z += -mu_total*2.0f/3.0f*(divu)*szz;

        // SST automatic wall treatment (methods/turbulence §6.5 (c)): 粗メッシュで分子勾配が τ_w を
        // 過小評価するため、接線せん断を modeled τ_w=ρu_τ² (有効壁粘性) に置換する。u_τ は壁関数で
        // 事前算出済み (ransWallFunction)。法線は単位法線 n、接線方向 ê_t はセル接線速度から取る。
        // y⁺→0 では u_τ²=νU_t/y より分子勾配に縮退、対数層では壁関数化し、全層で連続。
        if (wallTreatment == 1) {
            const flow_float n_x = sxx / sss;
            const flow_float n_y = syy / sss;
            const flow_float n_z = szz / sss;
            const flow_float uc = Ux[ic];
            const flow_float vc = Uy[ic];
            const flow_float wc = Uz[ic];
            const flow_float un = uc*n_x + vc*n_y + wc*n_z;
            const flow_float utx = uc - un*n_x;
            const flow_float uty = vc - un*n_y;
            const flow_float utz = wc - un*n_z;
            const flow_float Ut = sqrt(utx*utx + uty*uty + utz*utz);
            const flow_float utau_w = utau_b[ib];
            const flow_float tauw = ro[ic]*utau_w*utau_w;   // modeled 壁せん断応力の大きさ
            flow_float etx = 0.0f, ety = 0.0f, etz = 0.0f;
            if (Ut > 1.0e-12f) { etx = utx/Ut; ety = uty/Ut; etz = utz/Ut; }
            // 流れを減速させる向き (-ê_t) に τ_w·面積 を課す。法線粘性・体積項は落とす (壁関数)。
            tau_x = -tauw*etx*sss;
            tau_y = -tauw*ety*sss;
            tau_z = -tauw*etz*sss;
        }
        else if (wallTreatment == 2) {
            // WMLES: 壁モデル (wmlesWallModel) が算出済みの τ_w ベクトル (流束/面積, 向き -ê_∥ は
            // マッチング点瞬時速度) で置換。法線粘性・体積項は落とす (壁モデルが全応力を与える)。
            // node では壁ノード残差は Dirichlet でゼロ化されるため実効は cell のみ (書いても無害)。
            tau_x = twall_x_b[ib]*sss;
            tau_y = twall_y_b[ib]*sss;
            tau_z = twall_z_b[ib]*sss;
        }

        // 乱流熱伝導 (内部面と同じ k_eff = k_lam + cp*mu_turb/Pr_t)。
        flow_float tc_w = thermCond[ic] + cp[ic]*v_turb/Prt;
        // 熱流束: WMLES (==2) は壁モデル q_w·S。断熱壁 (adiabaticWall) は**厳密 0** — cell は
        // ミラーゴーストで従来もビット 0、node は従来 ∇T[W]·S の法線成分が漏れていたのを修正。
        // 等温壁のみ解像熱流束 (node: ∇T·S / cell: ghost 差分)。
        flow_float heatflux;
        if (wallTreatment == 2 || sstEnergyWf != 0)  heatflux = qwall_b[ib]*sss;
        else if (adiabaticWall != 0)    heatflux = (flow_float)0.0;
        else                            heatflux = (isNode != 0)
                                            ? tc_w*(dTdx[ic]*sxx +dTdy[ic]*syy +dTdz[ic]*szz)
                                            : tc_w*((Ts[ig]- Ts[ic])/dcc)*sss;

        flow_float res_ro_temp   = 0.0f;
        flow_float res_roUx_temp = tau_x;
        flow_float res_roUy_temp = tau_y;
        flow_float res_roUz_temp = tau_z;
        flow_float res_roe_temp  = tau_x*Ux_b[ib] +tau_y*Uy_b[ib] +tau_z*Uz_b[ib]; 

        //printf("res_roe_temp=%e\n", res_roe_temp);
        //printf("Uyb[ib]=%e\n", Uy_b[ib]);
        //printf("Uzb[ib]=%e\n", Uz_b[ib]);

        res_roe_temp += heatflux;

        //res_roe_temp = 0.0;

        atomicAdd(&res_ro[ic]  , res_ro_temp);
        atomicAdd(&res_roUx[ic], res_roUx_temp);
        atomicAdd(&res_roUy[ic], res_roUy_temp);
        atomicAdd(&res_roUz[ic], res_roUz_temp);
        atomicAdd(&res_roe[ic] , res_roe_temp);

        twall_x_b[ib] = tau_x/sss;
        twall_y_b[ib] = tau_y/sss;
        twall_z_b[ib] = tau_z/sss;

        flow_float twall = sqrt(tau_x*tau_x + tau_y*tau_y + tau_z*tau_z)/sss;
        flow_float utau = sqrt(twall/ro[ic]);

        // mode 1/2 では壁関数/壁モデルカーネルが y⁺=u_τ y/ν を既に格納済み。
        // ここで dcc/mu_total ベースの値で上書きすると定義が不整合になるため mode 0 のみ更新する。
        if (wallTreatment == 0) {
            ypls_b[ib] = ro[ic]*utau*dcc/mu_total;
        }
        //if (ib == 1) {
        //    printf("ib = %d\n", ib);
        //    printf("ip = %d\n", ip);
        //    printf("tau_x = %e, tau_y = %e, tau_z = %e\n", tau_x, tau_y, tau_z);
        //    printf("c Ux = %e, Uy = %e, Uz = %e\n", Ux[ic], Uy[ic], Uz[ic]);
        //    printf("g Ux = %e, Uy = %e, Uz = %e\n", Ux[ig], Uy[ig], Uz[ig]);
        //    printf("b Ux = %e, Uy = %e, Uz = %e\n", Ux_b[ib], Uy_b[ib], Uz_b[ib]);
        //}
    }

    __syncthreads();
}

// ============================================================================
// node-centered 壁摩擦応力 (twall) の専用算出カーネル (cfg.nodeWallStressEdgeKernel==1)
// ----------------------------------------------------------------------------
// 動機: node の壁ノード W は壁面上に乗るため、viscousFlux_wall_d の twall は退化ミラーゴースト dcc /
//   壁ノード勾配 ∇U[W]·S で評価され近壁で特異スパイク・偶奇振動する (plan §2.1/§6.2)。
// 方針 (plan §11, ユーザー確認済 案 A): 壁ノード W に接続する「内部双対面 (W↔内部ノード I)」の粘性
//   運動量 flux を W の CV に集約し、壁半割面積で割って twall とする。W は Dirichlet (u≈0) なので、これは
//   流体が壁 CV に及ぼす粘性 traction (magnitude=τ_w)。**運動量残差・y+ には触れない** (内部ノード I の
//   運動量は viscousFlux_d の内部双対面が既に担うため、ここで res に足すと二重計上になる; twall は
//   出力専用で場は不変)。壁端の速度は no-slip 値 Uxb=0 を使い (ghost 不使用)、面勾配は 1/2(∇[W]+∇[I])。
//   over-relaxed 分解は viscousFlux_d と同形だが、距離は退化 dcc でなく W-I 間の物理距離 |cc[I]-cc[W]|。
// 1 スレッド = 1 壁半割面 ib。W=bplane_cell[ib] の入射内部面を CSR (cell_planes) で走査して集約する。
// **このカーネルは出力専用 (res/y+ に一切触れない; 名前の通り twall を後処理として算出するだけ)**。
// 壁関数 active (Tau_Wall[W]>0) のときは、解像 traction を modeled 壁せん断 τ_w=ρu_τ² (Tau_Wall, AddTauWall と
// 同値) に再スケールして向きだけ解像値・大きさをモデル値に揃える (壁関数メッシュで解像 μ_total·du/dy が
// 過大評価する問題の補正。utau/y+ と整合)。壁解像 (Tau_Wall<0) では解像値のまま (それが真の τ_w)。
__global__ void wallStressForOutput_node_d
(
 geom_int nb,
 geom_int* bplane_plane,
 geom_int* bplane_cell,            // node: 壁ノード W
 geom_int  nNormalPlanes,
 geom_int* cell_planes_index, geom_int* cell_planes, geom_int* plane_cells,
 geom_int* wall_flag,

 geom_float* ccx , geom_float* ccy , geom_float* ccz,
 geom_float* sx  , geom_float* sy  , geom_float* sz , geom_float* ss,

 flow_float* vis_lam , flow_float* vis_turb,
 flow_float* Ux , flow_float* Uy , flow_float* Uz,
 flow_float* Ux_b , flow_float* Uy_b , flow_float* Uz_b,   // 壁面値 (no-slip=0)

 flow_float* dUxdx , flow_float* dUxdy , flow_float* dUxdz,
 flow_float* dUydx , flow_float* dUydy , flow_float* dUydz,
 flow_float* dUzdx , flow_float* dUzdy , flow_float* dUzdz,

 int isAxisymmetric, flow_float* axisym_divU,

 flow_float* twall_x_b , flow_float* twall_y_b , flow_float* twall_z_b,
 flow_float* Tau_Wall              // 壁関数 modeled τ_w=ρu_τ² (>0 で active, ransWallFunction が算出)。再スケール用。
)
{
    geom_int ib = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        const geom_int W   = bplane_cell[ib];
        const geom_int ipw = bplane_plane[ib];
        const flow_float ss_wall = ss[ipw];

        // 壁面値 (no-slip では 0)。ghost は使わない。
        const flow_float Uwx = Ux_b[ib];
        const flow_float Uwy = Uy_b[ib];
        const flow_float Uwz = Uz_b[ib];

        // W の CV に内部双対面から入る粘性運動量 flux の総和 (= 流体が壁 CV に及ぼす粘性力)
        flow_float fX = 0.0f, fY = 0.0f, fZ = 0.0f;

        for (geom_int j = cell_planes_index[W]; j < cell_planes_index[W+1]; ++j) {
            const geom_int ip = cell_planes[j];
            if (ip >= nNormalPlanes) continue;          // 内部双対面のみ (境界半割面は skip)
            const geom_int a = plane_cells[2*ip+0];
            const geom_int b = plane_cells[2*ip+1];
            const geom_int I = (a == W) ? b : a;        // 相手ノード
            if (wall_flag[I] != 0) continue;            // 内部ノードのみ (壁-壁エッジ除外, 寄与≈0)

            // 面法線 S を W->I 向きへ揃える (格納は a->b)。
            const flow_float sgn = (a == W) ? (flow_float)1.0 : (flow_float)-1.0f;
            const flow_float sxx = sgn*sx[ip];
            const flow_float syy = sgn*sy[ip];
            const flow_float szz = sgn*sz[ip];
            const flow_float sss = ss[ip];

            // W->I の物理距離ベクトル (退化ミラーゴースト dcc は使わない)。
            const flow_float dcx = ccx[I] - ccx[W];
            const flow_float dcy = ccy[I] - ccy[W];
            const flow_float dcz = ccz[I] - ccz[W];
            const flow_float dcc = sqrt(dcx*dcx + dcy*dcy + dcz*dcz);

            // over-relaxed 分解 (viscousFlux_d と同形)。
            const flow_float D_safe  = max(fabs(dcx*sxx + dcy*syy + dcz*szz), (flow_float)1.0e-30);
            const flow_float delta   = dcc*sss*sss/D_safe;
            const flow_float delta_x = dcx*sss*sss/D_safe;
            const flow_float delta_y = dcy*sss*sss/D_safe;
            const flow_float delta_z = dcz*sss*sss/D_safe;
            const flow_float k_x = sxx - delta_x;
            const flow_float k_y = syy - delta_y;
            const flow_float k_z = szz - delta_z;

            // 面勾配 = 1/2(∇[W]+∇[I]) (dual 平均, plan §11.3)。
            const flow_float dUxdxf = (flow_float)0.5*(dUxdx[W]+dUxdx[I]);
            const flow_float dUxdyf = (flow_float)0.5*(dUxdy[W]+dUxdy[I]);
            const flow_float dUxdzf = (flow_float)0.5*(dUxdz[W]+dUxdz[I]);
            const flow_float dUydxf = (flow_float)0.5*(dUydx[W]+dUydx[I]);
            const flow_float dUydyf = (flow_float)0.5*(dUydy[W]+dUydy[I]);
            const flow_float dUydzf = (flow_float)0.5*(dUydz[W]+dUydz[I]);
            const flow_float dUzdxf = (flow_float)0.5*(dUzdx[W]+dUzdx[I]);
            const flow_float dUzdyf = (flow_float)0.5*(dUzdy[W]+dUzdy[I]);
            const flow_float dUzdzf = (flow_float)0.5*(dUzdz[W]+dUzdz[I]);

            const flow_float divu = (isAxisymmetric == 1)
                ? (flow_float)0.5*(axisym_divU[W]+axisym_divU[I])
                : (dUxdxf + dUydyf + dUzdzf);

            const flow_float mu_total = (flow_float)0.5*
                ((vis_lam[W]+vis_turb[W]) + (vis_lam[I]+vis_turb[I]));

            // 法線 Laplacian は壁端 no-slip 値 Uw を使う: (U[I]-Uw)/dcc。W=ic0(from) なので符号は
            // viscousFlux_d の「ic0 に +tau」と同じ (= W の CV に入る力)。
            flow_float tx = mu_total*((Ux[I]-Uwx)/dcc)*delta;
            tx += mu_total*(dUxdxf*k_x + dUxdyf*k_y + dUxdzf*k_z);
            tx += mu_total*(dUxdxf*sxx + dUydxf*syy + dUzdxf*szz);
            tx += -mu_total*(flow_float)(2.0f/3.0f)*divu*sxx;

            flow_float ty = mu_total*((Uy[I]-Uwy)/dcc)*delta;
            ty += mu_total*(dUydxf*k_x + dUydyf*k_y + dUydzf*k_z);
            ty += mu_total*(dUxdyf*sxx + dUydyf*syy + dUzdyf*szz);
            ty += -mu_total*(flow_float)(2.0f/3.0f)*divu*syy;

            flow_float tz = mu_total*((Uz[I]-Uwz)/dcc)*delta;
            tz += mu_total*(dUzdxf*k_x + dUzdyf*k_y + dUzdzf*k_z);
            tz += mu_total*(dUxdzf*sxx + dUydzf*syy + dUzdzf*szz);
            tz += -mu_total*(flow_float)(2.0f/3.0f)*divu*szz;

            fX += tx; fY += ty; fZ += tz;
        }

        // 壁半割面積で割って traction (応力) に。res には一切加算しない。
        twall_x_b[ib] = fX/ss_wall;
        twall_y_b[ib] = fY/ss_wall;
        twall_z_b[ib] = fZ/ss_wall;

        // 壁関数 active (Tau_Wall[W]>0) のとき、向きは解像値・大きさを modeled τ_w=ρu_τ² に揃える
        // (壁関数メッシュで解像 μ_total·du/dy が ρu_τ² の ~十数倍に過大評価する補正、utau/y+ と整合)。
        // 壁解像 (Tau_Wall<0) では解像値=真の τ_w なので何もしない。
        if (Tau_Wall != nullptr) {
            const geom_int icW = bplane_cell[ib];
            const flow_float tauM = Tau_Wall[icW];
            if (tauM > static_cast<flow_float>(0.0)) {
                const flow_float twMag = sqrt(twall_x_b[ib]*twall_x_b[ib]
                                            + twall_y_b[ib]*twall_y_b[ib]
                                            + twall_z_b[ib]*twall_z_b[ib]);
                if (twMag > static_cast<flow_float>(1.0e-30)) {
                    const flow_float sc = tauM / twMag;
                    twall_x_b[ib] *= sc; twall_y_b[ib] *= sc; twall_z_b[ib] *= sc;
                }
            }
        }
    }

    __syncthreads();
}

// ---- 壁粘性 flux の距離診断 (env FORGE_VISC_WALL_DIAG=1 で 1 回だけ実行) ----
// 壁半割面ごとに、signed wall distance dn=(pc-cc)·n、ミラーゴースト距離 dcc=|cc_ghost-cc|、
// 接線オフセット |(pc-cc)-dn n| を書き出す。node-centered で dcc≈0 (壁ノードが壁面上) に
// 退化しているか、cell-centered で dcc≈2dn が成立しているかを切り分けるための診断。
__global__ void viscousWallDiag_d
(
 geom_int nb,
 geom_int* bplane_plane, geom_int* bplane_cell, geom_int* bplane_cell_ghst,
 geom_float* ccx, geom_float* ccy, geom_float* ccz,
 geom_float* pcx, geom_float* pcy, geom_float* pcz,
 geom_float* sx , geom_float* sy , geom_float* sz, geom_float* ss,
 flow_float* dn_o, flow_float* dcc_o, flow_float* tang_o
)
{
    geom_int ib = blockDim.x*blockIdx.x + threadIdx.x;
    if (ib < nb) {
        geom_int ip = bplane_plane[ib];
        geom_int ic = bplane_cell[ib];
        geom_int ig = bplane_cell_ghst[ib];
        flow_float nx = sx[ip]/ss[ip], ny = sy[ip]/ss[ip], nz = sz[ip]/ss[ip];
        flow_float dx = pcx[ip]-ccx[ic], dy = pcy[ip]-ccy[ic], dz = pcz[ip]-ccz[ic];
        flow_float dn = dx*nx + dy*ny + dz*nz;
        flow_float tx = dx - dn*nx, ty = dy - dn*ny, tz = dz - dn*nz;
        flow_float gx = ccx[ig]-ccx[ic], gy = ccy[ig]-ccy[ic], gz = ccz[ig]-ccz[ic];
        dn_o[ib]   = dn;
        dcc_o[ib]  = sqrt(gx*gx + gy*gy + gz*gz);
        tang_o[ib] = sqrt(tx*tx + ty*ty + tz*tz);
    }
}

// SST エネルギー壁関数 (§6.5(g)) の node AddQWall ゲート: node × SST × wallTreatmentSST==1 ×
// sstEnergyWallFunction==1 で wall_isothermal bcond が存在するとき true。この条件下でのみ
// Qw_Wall が初期化 (init_wf_pk_d) + 書き込み (compute_wall_friction_sst_d) されている。
static bool sstEnergyWfNodeActive(const solverConfig& cfg, const mesh& msh)
{
    if (!(cfg.discretization == "node" && msh.wall_flag_d != nullptr
          && cfg.LESorRANS == 2 && cfg.RANSmodel == 1
          && cfg.wallTreatmentSST == 1 && cfg.sstEnergyWallFunction == 1)) return false;
    for (const auto& bc : msh.bconds) if (bc.bcondKind == "wall_isothermal") return true;
    return false;
}

void viscousFlux_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var , matrix& mat_ns)
{
    // W-I 実力診断 (§4.2): env FORGE_WI_FORCE_DIAG=1 のときだけ有効 (通常経路の性能を変えない)。
    // 有効時は毎ステップ 0 クリアしてから積む (出力時点の値 = その step の合計)。
    static const bool wiDiagOn = [] {
        const char* s = std::getenv("FORGE_WI_FORCE_DIAG");
        const bool on = s && std::atoi(s) != 0;
        if (on) printf("[DIAG] W-I force diagnostic ON (FORGE_WI_FORCE_DIAG)\n");
        return on;
    }();
    if (wiDiagOn && var.c_d.count("wi_ftan")) {
        gpuErrchk(cudaMemset(var.c_d["wi_ftan"],     0, sizeof(flow_float)*msh.nCells_all));
        gpuErrchk(cudaMemset(var.c_d["wi_fnrm"],     0, sizeof(flow_float)*msh.nCells_all));
        gpuErrchk(cudaMemset(var.c_d["wi_fnrm_abs"], 0, sizeof(flow_float)*msh.nCells_all));
        gpuErrchk(cudaMemset(var.c_d["wi_ftan_res"], 0, sizeof(flow_float)*msh.nCells_all));
    }

    // 距離診断 (1 回限り)。FORGE_VISC_WALL_DIAG=1 のとき壁半割面の dn/dcc/tangential を集計表示。
    static bool diag_done = false;
    if (!diag_done && getenv("FORGE_VISC_WALL_DIAG") && atoi(getenv("FORGE_VISC_WALL_DIAG")) != 0) {
        diag_done = true;
        for (auto& bc : msh.bconds) {
            if (!(bc.bcondKind == "wall" || bc.bcondKind == "wall_isothermal")) continue;
            geom_int nbf = (geom_int)bc.iPlanes.size();
            if (nbf == 0) { printf("[viscWallDiag] bcond '%s' (physID %d): node 壁 plane なし\n",
                                    bc.bcondKind.c_str(), bc.physID); continue; }
            flow_float *dn_d, *dcc_d, *tang_d;
            gpuErrchk(cudaMalloc(&dn_d , sizeof(flow_float)*nbf));
            gpuErrchk(cudaMalloc(&dcc_d, sizeof(flow_float)*nbf));
            gpuErrchk(cudaMalloc(&tang_d,sizeof(flow_float)*nbf));
            viscousWallDiag_d<<<cuda_cfg.dimGrid_bplane, cuda_cfg.dimBlock>>>(
                nbf, bc.map_bplane_plane_d, bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
                var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["pcx"], var.p_d["pcy"], var.p_d["pcz"],
                var.p_d["sx"], var.p_d["sy"], var.p_d["sz"], var.p_d["ss"],
                dn_d, dcc_d, tang_d);
            gpuErrchkKernelSync();
            std::vector<flow_float> dn(nbf), dcc(nbf), tang(nbf);
            gpuErrchk(cudaMemcpy(dn.data() , dn_d , sizeof(flow_float)*nbf, cudaMemcpyDeviceToHost));
            gpuErrchk(cudaMemcpy(dcc.data(), dcc_d, sizeof(flow_float)*nbf, cudaMemcpyDeviceToHost));
            gpuErrchk(cudaMemcpy(tang.data(),tang_d,sizeof(flow_float)*nbf, cudaMemcpyDeviceToHost));
            // 統計
            double dn_mn=1e30,dn_mx=-1e30,dn_amn=1e30, dcc_mn=1e30,dcc_mx=-1e30;
            double tang_mx=-1e30, ratio_mn=1e30, ratio_mx=-1e30; int nFlip=0;
            for (geom_int i=0;i<nbf;i++){
                double adn=fabs(dn[i]);
                dn_mn=fmin(dn_mn,dn[i]); dn_mx=fmax(dn_mx,dn[i]); dn_amn=fmin(dn_amn,adn);
                dcc_mn=fmin(dcc_mn,dcc[i]); dcc_mx=fmax(dcc_mx,dcc[i]);
                tang_mx=fmax(tang_mx,tang[i]);
                double r = dcc[i]/fmax(2.0*adn,1e-30);
                ratio_mn=fmin(ratio_mn,r); ratio_mx=fmax(ratio_mx,r);
                if (dn[i]<0) nFlip++;
            }
            printf("[viscWallDiag] bcond '%s' physID %d nFaces=%ld disc=%s\n",
                   bc.bcondKind.c_str(), bc.physID, (long)nbf, cfg.discretization.c_str());
            printf("  dn (signed): min=%.3e max=%.3e | |dn| min=%.3e | dn<0 count=%d\n",
                   dn_mn, dn_mx, dn_amn, nFlip);
            printf("  dcc=|cc_ghost-cc|: min=%.3e max=%.3e\n", dcc_mn, dcc_mx);
            printf("  dcc/(2|dn|): min=%.3e max=%.3e (cell では ~1, node 退化なら ~0/大ばらつき)\n",
                   ratio_mn, ratio_mx);
            printf("  tangential offset |(pc-cc)-dn n|: max=%.3e (node 壁ノードでは支配的のはず)\n", tang_mx);
            cudaFree(dn_d); cudaFree(dcc_d); cudaFree(tang_d);
        }
    }

    // ------------------------------
    // *** sum over normal planes ***
    // ------------------------------
    viscousFlux_d<<<cuda_cfg.dimGrid_plane , cuda_cfg.dimBlock>>> ( 
        // mesh structure
        msh.nCells,
        msh.nPlanes , msh.nNormalPlanes , msh.map_plane_cells_d,
        var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
        var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
        var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

        cfg.visc , cfg.turbulentPrandtl , var.c_d["cp"] , var.c_d["thermCond"],
        var.c_d["vis_lam"], var.c_d["vis_turb"],

        // basic variables
        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["Ux"]  , 
        var.c_d["Uy"]  , 
        var.c_d["Uz"]  , 
        var.c_d["P"]  , 
        var.c_d["Ht"]  , 
        var.c_d["sonic"]  , 
        var.c_d["T"]  , 

        var.c_d["res_ro"] ,
        var.c_d["res_roUx"] ,
        var.c_d["res_roUy"] ,
        var.c_d["res_roUz"] ,
        var.c_d["res_roe"]  ,
       
        // gradient
        var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
        var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
        var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
        var.c_d["dTdx"] , var.c_d["dTdy"] , var.c_d["dTdz"],
        cfg.isAxisymmetric, var.c_d["axisym_divU"],
        // node SST 壁関数 / node WMLES のとき Tau_Wall を渡し AddTauWall 再スケール。
        // それ以外は nullptr (Tau_Wall 未初期化のため)。
        ((cfg.discretization == "node" && cfg.LESorRANS == 2 && cfg.RANSmodel == 1 && cfg.wallTreatmentSST == 1)
         || wmlesNodeActive(cfg, msh))
            ? var.c_d["Tau_Wall"] : nullptr,
        // 診断 (§4.2): W-I 実力集計。毎ステップ 0 クリアしてから渡す。
        (wiDiagOn && var.c_d.count("wi_ftan"))     ? var.c_d["wi_ftan"]     : nullptr,
        (wiDiagOn && var.c_d.count("wi_fnrm"))     ? var.c_d["wi_fnrm"]     : nullptr,
        (wiDiagOn && var.c_d.count("wi_fnrm_abs")) ? var.c_d["wi_fnrm_abs"] : nullptr,
        (wiDiagOn && var.c_d.count("wi_ftan_res")) ? var.c_d["wi_ftan_res"] : nullptr,
        // node WMLES 等温壁 / node SST エネルギー壁関数 (§6.5(g)) のとき Qw_Wall を渡し
        // AddQWall (W-I 熱流束置換)。それ以外は nullptr (Qw_Wall 未初期化のため)。
        (wmlesNodeIsothermalActive(cfg, msh) || sstEnergyWfNodeActive(cfg, msh))
            ? var.c_d["Qw_Wall"] : nullptr,
        // SST 断熱壁 SU2 式熱結合 (experimental mode 2) のときのみ overlay を渡す。
        // mode 0/1 は nullptr (従来経路ビット不変)。
        (cfg.discretization == "node" && cfg.LESorRANS == 2 && cfg.RANSmodel == 1
         && cfg.wallTreatmentSST == 1 && cfg.sstThermalWallFunction == 2)
            ? var.c_d["Taw_Prim_Overlay"] : nullptr,
        // 【実験専用】overlay 辺の壁側 μt スケール (FORGE_TAW_WALL_MUT_SCALE, 既定 1.0)
        [] {
            static const flow_float v = [] {
                const char* s = std::getenv("FORGE_TAW_WALL_MUT_SCALE");
                const flow_float x = s ? static_cast<flow_float>(std::atof(s)) : static_cast<flow_float>(1.0);
                if (s) printf("[EXPERIMENT] taw wall mut scale = %g (FORGE_TAW_WALL_MUT_SCALE)\n", (double)x);
                return x;
            }();
            return v;
        }(),
        // SST 断熱壁 defect-flux 閉包 (experimental mode 3) のときのみ H⃗/Taw を渡す。
        // mode 0/1/2 は nullptr (従来経路ビット不変)。
        (cfg.discretization == "node" && cfg.LESorRANS == 2 && cfg.RANSmodel == 1
         && cfg.wallTreatmentSST == 1 && cfg.sstThermalWallFunction == 3)
            ? var.c_d["Taw_HTnx"] : nullptr,
        (cfg.discretization == "node" && cfg.LESorRANS == 2 && cfg.RANSmodel == 1
         && cfg.wallTreatmentSST == 1 && cfg.sstThermalWallFunction == 3)
            ? var.c_d["Taw_HTny"] : nullptr,
        (cfg.discretization == "node" && cfg.LESorRANS == 2 && cfg.RANSmodel == 1
         && cfg.wallTreatmentSST == 1 && cfg.sstThermalWallFunction == 3)
            ? var.c_d["Taw_HTnz"] : nullptr,
        (cfg.discretization == "node" && cfg.LESorRANS == 2 && cfg.RANSmodel == 1
         && cfg.wallTreatmentSST == 1 && cfg.sstThermalWallFunction == 3)
            ? var.c_d["Taw_diag"] : nullptr,
        (cfg.LESorRANS == 2 && cfg.RANSmodel == 1 && var.c_d.count("k")) ? var.c_d["k"] : nullptr,
        cfg.sstIsotropicStress,
        // sstEnergyIncludesK: k 拡散のエネルギー流束 (SST のときのみ)
        var.c_d["dKdx"], var.c_d["dKdy"], var.c_d["dKdz"],
        var.c_d.count("sstF1") ? var.c_d["sstF1"] : nullptr,
        cfg.sstSigmaBlend,
        (cfg.sstEnergyIncludesK != 0 && cfg.LESorRANS == 2 && cfg.RANSmodel == 1) ? 1 : 0
    ) ;

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

    for (auto& bc : msh.bconds)
    {
        if (bc.bcondKind == "wall" or bc.bcondKind == "wall_isothermal") {
            viscousFlux_wall_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> ( 
                // mesh structure
                bc.iPlanes.size(),
                bc.map_bplane_plane_d,  
                bc.map_bplane_cell_d,  
                bc.map_bplane_cell_ghst_d,

                var.c_d["volume"], var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["pcx"]   , var.p_d["pcy"], var.p_d["pcz"], var.p_d["fx"],
                var.p_d["sx"]    , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],  

                cfg.visc , cfg.turbulentPrandtl , var.c_d["cp"] , var.c_d["thermCond"],
                var.c_d["vis_lam"], var.c_d["vis_turb"],

                // basic variables
                //var.c_d["convx"] , var.c_d["convy"] , var.c_d["convz"] ,
                //var.c_d["diffx"] , var.c_d["diffy"] , var.c_d["diffz"] ,
                var.c_d["ro"] ,
                var.c_d["roUx"] ,
                var.c_d["roUy"] ,
                var.c_d["roUz"] ,
                var.c_d["roe"] ,
                var.c_d["Ux"]  , 
                var.c_d["Uy"]  , 
                var.c_d["Uz"]  , 
                var.c_d["P"]  , 
                var.c_d["Ht"]  , 
                var.c_d["sonic"]  , 
                var.c_d["T"]  , 

                var.c_d["res_ro"] ,
                var.c_d["res_roUx"] ,
                var.c_d["res_roUy"] ,
                var.c_d["res_roUz"] ,
                var.c_d["res_roe"]  ,
       
                // gradient
                var.c_d["dUxdx"] , var.c_d["dUxdy"] , var.c_d["dUxdz"],
                var.c_d["dUydx"] , var.c_d["dUydy"] , var.c_d["dUydz"],
                var.c_d["dUzdx"] , var.c_d["dUzdy"] , var.c_d["dUzdz"],
                var.c_d["dTdx"]  , var.c_d["dTdy"]  , var.c_d["dTdz"],

                bc.bvar_d["ypls"],
                bc.bvar_d["twall_x"],
                bc.bvar_d["twall_y"],
                bc.bvar_d["twall_z"],
                bc.bvar_d["Ux"],
                bc.bvar_d["Uy"],
                bc.bvar_d["Uz"],
                bc.bvar_d["Ts"],

                //bc.bvar_d["sx"],
                //bc.bvar_d["sy"],
                //bc.bvar_d["sz"]
                cfg.isAxisymmetric, var.c_d["axisym_divU"],

                // 壁処理: SST automatic wall treatment (1, ransWallFunction 算出の u_τ) /
                //         WMLES 壁応力モデル (2, wmlesWallModel 算出の twall_*/qwall bvar)
                (cfg.LESorRANS == 2 && cfg.RANSmodel == 1) ? cfg.wallTreatmentSST
                    : (wmlesActiveForBcond(cfg, bc) ? 2 : 0),
                bc.bvar_d["utau"], bc.bvar_d["qwall"],
                (cfg.discretization == "node") ? 1 : 0,   // node: 壁法線/熱流束を ∇φ·S で評価 (ghostless)
                (bc.bcondKind == "wall") ? 1 : 0,         // 断熱壁: 伝導熱流束を厳密 0 (等温壁は 0=従来)
                // SST エネルギー壁関数 (§6.5(g)): 等温壁の壁面熱流束を Kader q_w に置換
                (cfg.LESorRANS == 2 && cfg.RANSmodel == 1 && cfg.wallTreatmentSST == 1
                 && cfg.sstEnergyWallFunction == 1 && bc.bcondKind == "wall_isothermal") ? 1 : 0
            ) ;
        }
    }

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();

    // node-centered の壁摩擦応力 (twall) を、退化ミラーゴースト dcc/壁ノード勾配ベース (上の
    // viscousFlux_wall_d が書いた値) から「壁ノード↔内部ノードの内部双対面集約」ベースに上書きする
    // (plan §11, 案 A)。出力専用で res/y+ には触れないため場は不変。cell モードと flag=0 では何もしない。
    if (cfg.discretization == "node" && cfg.nodeWallStressEdgeKernel == 1 && msh.wall_flag_d != nullptr) {
        for (auto& bc : msh.bconds) {
            if (!(bc.bcondKind == "wall" || bc.bcondKind == "wall_isothermal")) continue;
            if (bc.iPlanes.empty()) continue;
            wallStressForOutput_node_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> (
                bc.iPlanes.size(),
                bc.map_bplane_plane_d,
                bc.map_bplane_cell_d,
                msh.nNormalPlanes,
                msh.map_cell_planes_index_d, msh.map_cell_planes_d, msh.map_plane_cells_d,
                msh.wall_flag_d,

                var.c_d["ccx"], var.c_d["ccy"], var.c_d["ccz"],
                var.p_d["sx"] , var.p_d["sy"] , var.p_d["sz"] , var.p_d["ss"],

                var.c_d["vis_lam"], var.c_d["vis_turb"],
                var.c_d["Ux"], var.c_d["Uy"], var.c_d["Uz"],
                bc.bvar_d["Ux"], bc.bvar_d["Uy"], bc.bvar_d["Uz"],

                var.c_d["dUxdx"], var.c_d["dUxdy"], var.c_d["dUxdz"],
                var.c_d["dUydx"], var.c_d["dUydy"], var.c_d["dUydz"],
                var.c_d["dUzdx"], var.c_d["dUzdy"], var.c_d["dUzdz"],

                cfg.isAxisymmetric, var.c_d["axisym_divU"],

                bc.bvar_d["twall_x"], bc.bvar_d["twall_y"], bc.bvar_d["twall_z"],
                // 壁関数 / WMLES active のときだけ Tau_Wall=ρu_τ² を渡し twall を再スケール。それ以外は
                // nullptr (Tau_Wall 未初期化のため; カーネルは nullptr で再スケール無し=解像値のまま)。
                ((cfg.LESorRANS == 2 && cfg.RANSmodel == 1 && cfg.wallTreatmentSST == 1)
                 || wmlesActiveForBcond(cfg, bc))
                    ? var.c_d["Tau_Wall"] : nullptr
            ) ;
        }
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchkKernelSync();
    }
}