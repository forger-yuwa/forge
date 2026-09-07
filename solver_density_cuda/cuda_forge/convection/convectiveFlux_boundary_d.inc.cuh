// =============================================================================
// 対流フラックス convectiveFlux_boundary_d の実装断片。
//   - convectiveFlux_d.cu に textual include される断片であり **スタンドアロン TU ではない**。
//     単一 TU を維持するため (rdc 不要 / __device__ グローバル・cudaMemcpyToSymbol の共有を保つ)、
//     CMake の source には .cu/.cuh として足さないこと。
//   - 参照する共通ヘルパ (interp_dispatch 等) と診断グローバルは include 元の上方で定義済み。
// =============================================================================

__global__ void convectiveFlux_boundary_d // slau
(
 flow_float ga,
 int advGauge,   // 移流基準差分 (plan §8.5): 主ループ KEEP_d の advGauge と厳密に同条件で on/off
                 // (host 判定 solver==KEEP && CPG && roRef>0)。ゲージは CV の全面に載らないと
                 // telescoping が破れるため、片側だけの有効化は禁止。
 // node × outlet_statPress のみ 1: 境界圧力流束 p_tilde を bvar Psb (規定背圧, 入力専用)
 // でなく内部値 Ps[ic] で評価する。超音速流出で背圧を読み出口列が沈む欠陥の根治
 // (plan boundary-node-nozzle-wall-outlet-stability §2.11)。亜音速の背圧アンカーは
 // 特性構成された bvar (rob/Uxb 等) の質量流束が担う。入口/壁/slip は 0 (従来どおり
 // bvar 圧力 = BC の効き筋 or P[ic] と同値)。
 int pTildeInterior,

  // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  

 geom_float* vol ,  geom_float* ccx ,  geom_float* ccy, geom_float* ccz,
 geom_float* pcx ,  geom_float* pcy ,  geom_float* pcz, geom_float* fx,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 flow_float* massflux,  // node 弱形式境界の質量流束を書き戻す (スカラ移流 k/ω/species が境界で読む)

 // variables
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

 // bvar
 flow_float* rob ,
 flow_float* roUxb ,
 flow_float* roUyb ,
 flow_float* roUzb ,
 flow_float* roeb ,
 flow_float* Uxb ,
 flow_float* Uyb ,
 flow_float* Uzb ,
 flow_float* Ttb ,
 flow_float* Ptb ,
 flow_float* Tsb ,
 flow_float* Psb ,
 
 flow_float* res_ro   ,
 flow_float* res_roUx  ,
 flow_float* res_roUy  ,
 flow_float* res_roUz  ,
 flow_float* res_roe   ,
 // sstEnergyIncludesK (plan turbulence-sst-energy-includes-k): 内部 k / 境界 k (nullptr なら内部値) / 有効フラグ
 flow_float* kturb, flow_float* kb, int energyK
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        flow_float dro;
        flow_float dP;
        flow_float dUx;
        flow_float dUy;
        flow_float dUz;
        flow_float dU;

        geom_float sxx = sx[ip];
        geom_float syy = sy[ip];
        geom_float szz = sz[ip];
        geom_float sss = ss[ip];

        geom_float nx = sxx/sss;
        geom_float ny = syy/sss;
        geom_float nz = szz/sss;

        flow_float ccx_0 = ccx[ic];
        flow_float ccy_0 = ccy[ic];
        flow_float ccz_0 = ccz[ic];

        flow_float pcx_1 = pcx[ip];
        flow_float pcy_1 = pcy[ip];
        flow_float pcz_1 = pcz[ip];

        flow_float dcc_x = pcx_1 - ccx_0;
        flow_float dcc_y = pcy_1 - ccy_0;
        flow_float dcc_z = pcz_1 - ccz_0;
        flow_float dcc   = sqrt(dcc_x*dcc_x +dcc_y*dcc_y +dcc_z*dcc_z) ;

        flow_float ro_L = ro[ic];
        flow_float Ux_L = Ux[ic];
        flow_float Uy_L = Uy[ic];
        flow_float Uz_L = Uz[ic];
        flow_float P_L  = Ps[ic];
        flow_float Ht_L = Ht[ic];

        flow_float ro_R= rob[ib]; 
        flow_float Ux_R= Uxb[ib];
        flow_float Uy_R= Uyb[ib];
        flow_float Uz_R= Uzb[ib];
        flow_float P_R = Psb[ib];
        //flow_float Ht_R= Htb[ib];
        flow_float Ht_R= (roeb[ib] + P_R)/ro_R;

        flow_float ro_p= ro_L;
        flow_float ro_m= ro_R;
        flow_float u_p = Ux_L;
        flow_float v_p = Uy_L;
        flow_float w_p = Uz_L;
        flow_float h_p = Ht_L;
        flow_float u_m = Ux_R;
        flow_float v_m = Uy_R;
        flow_float w_m = Uz_R;
        flow_float h_m = Ht_R;
        flow_float P_p = P_L;
        flow_float P_m = P_R;

        //flow_float Vn_p = ((Ux_L)*sxx +(Uy_L)*syy +(Uz_L)*szz)/sss;
        flow_float Vn_m = ((Ux_R)*sxx +(Uy_R)*syy +(Uz_R)*szz)/sss;

        //flow_float roVn_p = (ro_L*Ux_L*sxx + ro_L*Uy_L*syy + ro_L*Uz_L*szz)/sss;
        //flow_float roVn_m = (ro_R*Ux_R*sxx + ro_R*Uy_R*syy + ro_R*Uz_R*szz)/sss;

        //flow_float VnL = Vn_p;
        //flow_float VnR = Vn_m;

        //flow_float VnL_abs = abs(VnL);
        //flow_float VnR_abs = abs(VnR);

//TODO: change c
        //flow_float c_hat = 0.5*(sonic[ic] + sonic[ig]);

        //flow_float M_p = Vn_p/c_hat;
        //flow_float M_m = Vn_m/c_hat;

        //flow_float ro_del= ro_m - ro_p;
        //flow_float Vn_del= Vn_m - Vn_p;
        //flow_float P_del = P_m  - P_p;

        //flow_float beta_p, beta_m;

        //if (abs(M_p)>=1.0){
        //    beta_p = 0.5*(M_p + abs(M_p))/M_p;
        //} else {
        //    beta_p = 0.25*pow((M_p+1.0),2.0)*(2.0-M_p);
        //}

        //if (abs(M_m)>=1.0){
        //    beta_m = 0.5*(M_m - abs(M_m))/M_m;
        //} else {
        //    beta_m = 0.25*pow((M_m-1.0),2.0)*(2.0+M_m);
        //}

        //flow_float zero = 0.0;
        //flow_float one  = 1.0;
        //flow_float half = 0.5;

        //flow_float g = -max(min(M_p,zero),-one)*min(max(M_m,zero),one);
        //flow_float Vn_hat_abs   = (ro_p*abs(Vn_p) + ro_m*abs(Vn_m))/(ro_p + ro_m);
        //flow_float Vn_hat_p_abs = (1.0-g)*Vn_hat_abs + g*VnL_abs;
        //flow_float Vn_hat_m_abs = (1.0-g)*Vn_hat_abs + g*VnR_abs;

        //flow_float M_hat = min(one, sqrt(half*(u_m*u_m+v_m*v_m+w_m*w_m +u_p*u_p+v_p*v_p+w_p*w_p))/c_hat);
        //flow_float chi = (1.0-M_hat)*(1.0-M_hat);

        flow_float p_tilde = (pTildeInterior != 0) ? Ps[ic] : P_R;
        // SST 全エネルギー E_t (sstEnergyIncludesK): h に +(5/3)k、圧力流束に p* = p + (2/3)ρk (境界側 k は bvar kb / 入口 k / 内部値)。
        if (energyK != 0 && kturb != nullptr) {
            const flow_float kL = max(kturb[ic], (flow_float)0.0);
            const flow_float kR = (kb != nullptr) ? max(kb[ib], (flow_float)0.0) : kL;   // 入口は bvar k (指定値)、他は内部値 (Neumann)
            h_p += (flow_float)(5.0/3.0)*kL;
            h_m += (flow_float)(5.0/3.0)*kR;
            p_tilde += (pTildeInterior != 0) ? (flow_float)(2.0/3.0)*ro_L*kL : (flow_float)(2.0/3.0)*ro_R*kR;
        }

        flow_float mdot = sss*(ro_R*Vn_m);

        // 質量流束を書き戻す: node 弱形式では主ループが境界半割面を除外する (5a) ため massflux[ip]=0 のまま残り、
        // スカラ移流 (k/ω/species, scalar_advection_first_order_d) が境界で流出入を 0 と誤認して k 等が出口で
        // 流出できず蓄積→過剰乱流になる。境界の mdot をここで書き戻し正しく移流させる (符号は ic0→ghost 外向き)。
        // cell モードは主ループが境界 massflux を書くので二重書き回避のため nullptr (書かない)。
        if (massflux != nullptr) massflux[ip] = mdot;

        flow_float res_ro_temp, res_roUx_temp, res_roUy_temp, res_roUz_temp, res_roe_temp;
        flow_float p_tilde_r = p_tilde - d_pRef;   // free-stream 保存: 基準静圧を差し引く (境界面)
        if (advGauge == 0) {
        res_ro_temp   = mdot;
        res_roUx_temp = 0.5*(mdot+abs(mdot))*u_p +0.5*(mdot-abs(mdot))*u_m +p_tilde_r*sxx;
        res_roUy_temp = 0.5*(mdot+abs(mdot))*v_p +0.5*(mdot-abs(mdot))*v_m +p_tilde_r*syy;
        res_roUz_temp = 0.5*(mdot+abs(mdot))*w_p +0.5*(mdot-abs(mdot))*w_m +p_tilde_r*szz;
        res_roe_temp  = 0.5*(mdot+abs(mdot))*h_p +0.5*(mdot-abs(mdot))*h_m ;
        } else {
        // ---- advGauge (plan §8.5): F − F∞(s) を差分形因数分解で組む ----
        // mdot は R 状態に線形なので a+b=mdot の厳密恒等式で upwind のまま因数分解できる:
        //   質量:   dCS = S[(ρ_R−ρ∞)Un_R + ρ∞ (u_R−u∞)·n]            (= mdot − C∞S)
        //   運動量: a(u_p−u∞) + b(u_m−u∞) + dCS·u∞ (+ (p_R−pRef)s)
        //   エネルギー: a(h_p−H∞) + b(h_m−H∞) + dCS·H∞
        // 一様場 (bvar=自由流) で全因子がビット単位ゼロ。壁/slip (Un_R=0 強制) では −F∞(s) 級だが
        // 厳密な定数係数×s なので CV 内の他面 (主ループ KEEP_d のゲージ済流束) と telescoping で相殺。
        // 式は tools/verify_advective_gauge.py (B1)-(B4) で数値検証済。massflux は物理 mdot のまま。
        const flow_float a = 0.5*(mdot+abs(mdot));
        const flow_float b = 0.5*(mdot-abs(mdot));
        const flow_float dLx = u_p - d_uRefX, dLy = v_p - d_uRefY, dLz = w_p - d_uRefZ;
        const flow_float dRx = u_m - d_uRefX, dRy = v_m - d_uRefY, dRz = w_m - d_uRefZ;
        const flow_float dUnR = (dRx*sxx + dRy*syy + dRz*szz)/sss;
        const flow_float dCS  = sss*((ro_R - d_roRef)*Vn_m + d_roRef*dUnR);   // = mdot − C∞S
        const flow_float kRef = 0.5*(d_uRefX*d_uRefX + d_uRefY*d_uRefY + d_uRefZ*d_uRefZ);
        const flow_float Hinf = ga/(ga-1.0)*d_pRef/d_roRef + kRef;
        res_ro_temp   = dCS;
        res_roUx_temp = a*dLx + b*dRx + dCS*d_uRefX + p_tilde_r*sxx;
        res_roUy_temp = a*dLy + b*dRy + dCS*d_uRefY + p_tilde_r*syy;
        res_roUz_temp = a*dLz + b*dRz + dCS*d_uRefZ + p_tilde_r*szz;
        res_roe_temp  = a*(h_p - Hinf) + b*(h_m - Hinf) + dCS*Hinf;
        }

        atomicAdd(&res_ro[ic]  , -res_ro_temp);
        atomicAdd(&res_roUx[ic], -res_roUx_temp);
        atomicAdd(&res_roUy[ic], -res_roUy_temp);
        atomicAdd(&res_roUz[ic], -res_roUz_temp);
        atomicAdd(&res_roe[ic] , -res_roe_temp);

    }

    __syncthreads();
}
