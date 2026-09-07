#include "cuda_forge/boundaryCond_d.cuh"
#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/thermo_d.cuh"
#include "cuda_forge/speciesTransport_d.cuh"  // species_roY_device_ptr()

// 内部セル ic の正規化質量分率 Y を取り出す (多成分 ghost の熱再構成用)。
// 壁/slip/outflow/outlet の ghost 組成は隣接内部セルと同一 (zero-gradient) なので、
// TP の R/cp/h/e/γ をその混合組成で再構成する。roY==nullptr または n<2 (単成分) では
// false を返し、呼び出し側は従来の sp[0] 経路をビット不変で維持する。
// 共通 helper (thermo_d.cuh の thermo_cell_Y) への転送。旧名は BC カーネル多数の呼び出し互換用。
__device__ inline bool bc_cell_Y(flow_float* const* roY, int n, geom_int ic, double* Y)
{
    return thermo_cell_Y(roY, n, ic, Y);
}

__global__
void slip_d
(
 // gas properties
 flow_float ga,
 flow_float cp,
 int thermalMethod, const SpeciesThermo* sp,
 flow_float* const* roY, int nSpecies,   // 多成分: ghost 組成 = 内部セル組成

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic  ,
 flow_float* T   ,

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
 flow_float* Psb
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        geom_float Un =  (sx[ip]*Ux[ic] + sy[ip]*Uy[ic] + sz[ip]*Uz[ic])/ss[ip];

        const flow_float Ux_ig = Ux[ic] - 2 * Un * sx[ip] / ss[ip];
        const flow_float Uy_ig = Uy[ic] - 2 * Un * sy[ip] / ss[ip];
        const flow_float Uz_ig = Uz[ic] - 2 * Un * sz[ip] / ss[ip];
        const flow_float Ux_b = Ux[ic] - Un * sx[ip] / ss[ip];
        const flow_float Uy_b = Uy[ic] - Un * sy[ip] / ss[ip];
        const flow_float Uz_b = Uz[ic] - Un * sz[ip] / ss[ip];

        ro[ig]   = ro[ic];
        P[ig]    = P[ic];
        Ux[ig]   = Ux_ig;
        Uy[ig]   = Uy_ig;
        Uz[ig]   = Uz_ig;


        roUx[ig] = ro[ic] * Ux_ig;
        roUy[ig] = ro[ic] * Uy_ig;
        roUz[ig] = ro[ic] * Uz_ig;

        const flow_float ek_ig = 0.5*(Ux_ig*Ux_ig + Uy_ig*Uy_ig + Uz_ig*Uz_ig);
        const flow_float ek_b  = 0.5*(Ux_b*Ux_b   + Uy_b*Uy_b   + Uz_b*Uz_b);
        rob[ib]   = ro[ic];
        Psb[ib]   = P[ic];
        if (thermalMethod == 2) {
            // slip 反射は熱力学状態を保存。e は内部温度・内部組成の NASA 値 (共通 helper)。
            const GasStateAtT gs = thermo_state_at_T(thermalMethod, ga, cp, sp, roY, nSpecies, ic, T[ic]);
            roe[ig]   = ro[ic]*(gs.e + ek_ig);
            Ht[ig]    = (roe[ig] + P[ig])/ro[ig];
            sonic[ig] = sonic[ic];
            T[ig]     = T[ic];
            roeb[ib]  = ro[ic]*(gs.e + ek_b);
            Tsb[ib]   = T[ic];
        } else {
            roe[ig]  = P[ic] / (ga - 1.0) + ro[ic] * ek_ig;
            Ht[ig]   = (roe[ig] + P[ig])/ro[ig];
            sonic[ig]= sqrt(ga*P[ig]/ro[ig]);
            T[ig]    = P[ig]*ga/(ro[ig]*(ga-1.0)*cp);
            roeb[ib] = P[ic] / (ga - 1.0) + ro[ic] * ek_b;
            Tsb[ib]  = Psb[ib]*ga/(rob[ib]*(ga-1.0)*cp);
        }

        Uxb[ib]   = Ux_b;
        Uyb[ib]   = Uy_b;
        Uzb[ib]   = Uz_b;

        roUxb[ib] = ro[ic] * Ux_b;
        roUyb[ib] = ro[ic] * Uy_b;
        roUzb[ib] = ro[ic] * Uz_b;
    }
};

void slip_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    slip_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>(
        cfg.gamma,
        cfg.cp,
        cfg.thermalMethod, thermo_species_device_ptr(),
        species_roY_device_ptr(), cfg.nSpecies,

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,
        bc.map_bplane_cell_d,
        bc.map_bplane_cell_ghst_d,

        var.p_d["sx"],
        var.p_d["sy"],
        var.p_d["sz"],
        var.p_d["ss"],

        var.c_d["ro"],
        var.c_d["roUx"],
        var.c_d["roUy"],
        var.c_d["roUz"],
        var.c_d["roe"],
        var.c_d["Ux"],
        var.c_d["Uy"],
        var.c_d["Uz"],
        var.c_d["P"],
        var.c_d["Ht"],
        var.c_d["sonic"],
        var.c_d["T"],

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]
    );
}

__global__ void wall_d
(

 // gas properties
 flow_float ga,
 flow_float cp,
 int thermalMethod, const SpeciesThermo* sp,
 flow_float* const* roY, int nSpecies,   // 多成分: ghost 組成 = 内部セル組成

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,
 geom_int* bplane_cell,
 geom_int* bplane_cell_ghst,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 geom_float* fx  ,

 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
flow_float* sonic,
flow_float* T,

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
 flow_float* Psb
)
{
    geom_int ib = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];


        ro[ig]   = ro[ic];
        Ux[ig]   = -roUx[ic]/ro[ic];
        Uy[ig]   = -roUy[ic]/ro[ic];
        Uz[ig]   = -roUz[ic]/ro[ic];
        roUx[ig] = -roUx[ic];
        roUy[ig] = -roUy[ic];
        roUz[ig] = -roUz[ic];
        roe[ig]  = roe[ic];

        flow_float ek = 0.5*(Ux[ig]*Ux[ig] +Uy[ig]*Uy[ig] +Uz[ig]*Uz[ig]);
        flow_float R;
        flow_float Tc;
        if (thermalMethod == 2) {
            // 鏡像反射は熱力学状態 (ρ,e,P,T) を保存するため内部値をコピー (NASA 整合)。
            double Yc[THERMO_MAX_SPECIES];
            const bool mix = bc_cell_Y(roY, nSpecies, ic, Yc);
            R  = (flow_float)(mix ? thermo_R_mix(sp, nSpecies, Yc) : thermo_R_species(sp[0]));
            P[ig]     = P[ic];
            Ht[ig]    = Ht[ic];
            sonic[ig] = sonic[ic];
            Tc        = T[ic];
        } else {
            P[ig] =(ga-1.0)*(roe[ig]-ro[ig]*ek);
            Ht[ig]   = ga*roe[ig]/ro[ig] + (1.0-ga)*ek;
            sonic[ig]= sqrt(ga*P[ig]/ro[ig]);
            R  = (ga-1.0)/ga*cp;
            Tc = P[ic]/(ro[ic]*R);
        }

        // 境界値 (bvar): node-centered 弱形式が直接消費する壁面値。no-slip なので速度=0 を明示
        // (対流境界 flux=圧力のみ・勾配 φ_b 速度=0・粘性壁応力の no-work を保証)。
        T[ig]     = T[ic];
        Tsb[ib]   = Tc;
        Psb[ib]   = P[ic];
        rob[ib]   = Psb[ib]/(R*Tc);
        Uxb[ib]   = (flow_float)0.0;
        Uyb[ib]   = (flow_float)0.0;
        Uzb[ib]   = (flow_float)0.0;
        roUxb[ib] = (flow_float)0.0;
        roUyb[ib] = (flow_float)0.0;
        roUzb[ib] = (flow_float)0.0;
        // roeb: 壁面の全エネルギー (KE=0 なので内部エネルギーのみ)。内部セルの内部エネルギーを採る
        // (断熱壁 ~ 内部温度)。CPG/TP 双方で成立。
        roeb[ib]  = roe[ic] - (flow_float)0.5*(roUx[ic]*roUx[ic] + roUy[ic]*roUy[ic] + roUz[ic]*roUz[ic])/ro[ic];
    }
};

void wall_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    wall_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> (
        cfg.gamma,
        cfg.cp,
        cfg.thermalMethod, thermo_species_device_ptr(),
        species_roY_device_ptr(), cfg.nSpecies,

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,
        bc.map_bplane_cell_d,
        bc.map_bplane_cell_ghst_d,

        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  
        var.p_d["fx"],  

        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["P"], 
        var.c_d["Ht"], 
        var.c_d["sonic"],
        var.c_d["T"],

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]

    ) ;
}

__global__ void wall_isothermal_d
(

 // gas properties
 flow_float ga,
 flow_float cp,
 int thermalMethod, const SpeciesThermo* sp,
 flow_float* const* roY, int nSpecies,   // 多成分: ghost 組成 = 内部セル組成

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 geom_float* fx  , 

 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,
 flow_float* T,

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
 flow_float* Psb

)
{
    geom_int ib = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        // R は TP では内部セル組成の混合 (NASA), CPG では cp-cp/ga。
        double Yc[THERMO_MAX_SPECIES];
        const bool mix = bc_cell_Y(roY, nSpecies, ic, Yc);
        flow_float R = (thermalMethod == 2)
                       ? (flow_float)(mix ? thermo_R_mix(sp, nSpecies, Yc) : thermo_R_species(sp[0]))
                       : (cp - cp/ga);

        Psb[ib] = P[ic];
        rob[ib] = Psb[ib]/(R*Tsb[ib]);
        // ghost 温度は鏡像外挿 (セル-ghost 中点=壁面で T=Tw になる)。Tw 直置きは壁熱流束
        // (Ts[ig]-Ts[ic])/dcc が dcc=2y1 と組んで正しい値の 1/2 に縮む (純伝導検証 2026-07-20,
        // case/24 run_isoT_cond* で確定)。強冷却壁 (Tc>2Tw) での負温度は 0.2Tw でクランプ。
        const flow_float Tg = max(static_cast<flow_float>(2.0)*Tsb[ib] - T[ic],
                                  static_cast<flow_float>(0.2)*Tsb[ib]);
        T[ig]   = Tg;

        ro[ig]   = Psb[ib]/(R*Tg);
        // 移動壁対応: ghost 速度は壁速度 (bvar, YAML Ux/Uy/Uz) 回りの鏡像 2U_w−U_c
        // (面平均が U_w になる)。U_w=0 (既定) では従来の −U_c とビット同一。粘性仕事項
        // τ·U_b は viscousFlux_wall_d が bvar 面速度で既に加算する。
        Ux[ig]   = static_cast<flow_float>(2.0)*Uxb[ib] - roUx[ic]/ro[ic];
        Uy[ig]   = static_cast<flow_float>(2.0)*Uyb[ib] - roUy[ic]/ro[ic];
        Uz[ig]   = static_cast<flow_float>(2.0)*Uzb[ib] - roUz[ic]/ro[ic];
        // ghost 保存量は ghost の原始変数と整合させる: roU = ro[ig]·U[ig]。
        // Ux[ig] には既に鏡像 (2U_w−U_c) の符号が入っているため負号を付けると二重反転し、
        // ghost 運動量が内部セルと同じ向きになる。密度も ghost 値 ro[ig] でなく壁面値
        // rob[ib] を使っていた。断熱壁 wall_d は roUx[ig]=-roUx[ic] (=ro[ig]·Ux[ig]) で
        // 元から整合しており、等温壁だけが不整合だった。
        // 注: cell モードの ghost 構成であり、node (median-dual) の壁ノードは別経路
        // (case/38 run_0018 の y₁⁺=1 発散はこの修正では直らないことを実測済み)。
        roUx[ig] = ro[ig]*Ux[ig];
        roUy[ig] = ro[ig]*Uy[ig];
        roUz[ig] = ro[ig]*Uz[ig];
        P[ig]    = P[ic];
        if (thermalMethod == 2) {
            // 等温壁 ghost を T=Tg (鏡像) と整合: roe=ρ(e_NASA(Tg)+ek), sonic=√(γ_mix P/ρ)。
            const flow_float ekg = 0.5f*(Ux[ig]*Ux[ig] + Uy[ig]*Uy[ig] + Uz[ig]*Uz[ig]);
            const GasStateAtT gw = thermo_state_at_T(thermalMethod, ga, cp, sp, roY, nSpecies, ic, Tg);
            roe[ig]  = ro[ig]*(gw.e + ekg);
            sonic[ig]= sqrt(gw.gamma*P[ig]/ro[ig]);
        } else {
            const flow_float ekg = 0.5f*(Ux[ig]*Ux[ig] + Uy[ig]*Uy[ig] + Uz[ig]*Uz[ig]);
            roe[ig]  = P[ig]/(ga - static_cast<flow_float>(1.0)) + ro[ig]*ekg;
            sonic[ig]= sqrt(ga*P[ig]/ro[ig]);
        }
        Ht[ig]   = (roe[ig] + P[ig])/ro[ig];
//
//        flow_float ek = 0.5*(Ux[ig]*Ux[ig] +Uy[ig]*Uy[ig] +Uz[ig]*Uz[ig]);
//        P[ig] =(ga-1.0)*(roe[ig]-ro[ig]*ek);
//        Ht[ig]   = (roe[ig] + P[ig])/ro[ig];
//        sonic[ig]= sqrt(ga*P[ig]/ro[ig]);
//
//        rob[ib]   = ro[ic];
//        //Uxb[ib]   = roUx[ic]/ro[ic];
//        //Uyb[ib]   = roUy[ic]/ro[ic];
//        //Uzb[ib]   = roUz[ic]/ro[ic];
//        Uxb[ib]   = 0.0;
//        Uyb[ib]   = 0.0;
//        Uzb[ib]   = 0.0;
//        roUxb[ib] = 0.0;
//        roUyb[ib] = 0.0;
//        roUzb[ib] = 0.0;
//        roeb[ib]  = roe[ic];
//
//        //flow_float ek = 0.5*(Uxb[ib]*Uxb[ib] +Uyb[ib]*Uyb[ib] +Uzb[ib]*Uzb[ib]);
//        Psb[ib] =(ga-1.0)*(roeb[ib]-rob[ib]*ek);
//        //Tsb[ib] =Psb[ib]*ga/(rob[ib]*(ga-1.0)*cp);
//
        flow_float Ux_b = Uxb[ib];
        flow_float Uy_b = Uyb[ib];
        flow_float Uz_b = Uzb[ib];
        flow_float ek = 0.5*(Ux_b*Ux_b +Uy_b*Uy_b +Uz_b*Uz_b);

        roUxb[ib] = rob[ib]*Ux_b;
        roUyb[ib] = rob[ib]*Uy_b;
        roUzb[ib] = rob[ib]*Uz_b;
        Psb[ib]   = P[ic];
        rob[ib]   = Psb[ib]/(R*Tsb[ib]);
        if (thermalMethod == 2) {
            const GasStateAtT gw = thermo_state_at_T(thermalMethod, ga, cp, sp, roY, nSpecies, ic, Tsb[ib]);
            roeb[ib] = rob[ib]*(gw.e + ek);
        } else {
            roeb[ib]  = Psb[ib]/(ga-1.0) + rob[ib]*ek;
        }
    }
};

void wall_isothermal_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    wall_isothermal_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> (
        cfg.gamma,
        cfg.cp,
        cfg.thermalMethod, thermo_species_device_ptr(),
        species_roY_device_ptr(), cfg.nSpecies,

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,
        bc.map_bplane_cell_d,
        bc.map_bplane_cell_ghst_d,

        var.p_d["sx"],
        var.p_d["sy"],
        var.p_d["sz"],
        var.p_d["ss"],  
        var.p_d["fx"],  

        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["P"], 
        var.c_d["Ht"], 
        var.c_d["sonic"],
        var.c_d["T"],

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]

    ) ;
}



__global__ 
void outlet_statPress_d
(
 // gas properties
 flow_float ga,
 flow_float cp,
 int thermalMethod, const SpeciesThermo* sp,
 flow_float* const* roY, int nSpecies,   // 多成分: ghost/backflow 組成 = 内部セル組成
 const flow_float* gamma_cell,           // TP: per-cell γ_mix(T) (dependentVariables)。特性構成の γ に使う (下記)

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,
 geom_int* bplane_cell,
 geom_int* bplane_cell_ghst,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic ,
 flow_float* T   ,

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
flow_float* Psb


)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        // ghost / backflow 流入ガスの組成は隣接内部セルと同一 (zero-gradient)。
        double Yc[THERMO_MAX_SPECIES];
        const bool mix = bc_cell_Y(roY, nSpecies, ic, Yc);

        //if (ib == 0) {
        //    printf("outlet_statPress_d: ib = %d, ip = %d, ic = %d, ig = %d\n", ib, ip, ic, ig);
        //}

        flow_float sxx = sx[ip];
        flow_float syy = sy[ip];
        flow_float szz = sz[ip];
        flow_float sss = ss[ip];

        //if (ib == 0) {
        //    printf("outlet_statPress_d: sxx = %f, syy = %f, szz = %f, sss = %f\n", sxx, syy, szz, sss);
        //}

        // 規定背圧 P_exit は bvar Psb を**入力専用**として読む (面ごとの分布指定
        // [inletProfile CSV の Ps 列など] を保持するため、本カーネルは Psb を書き換えない)。
        // node 弱形式の境界圧力流束 p_tilde と勾配境界閉包は、outlet では内部値 P[ic]/T[ic] を
        // 参照する (流束/勾配カーネル側の outlet フラグ。超音速流出で背圧を読む旧欠陥の根治 —
        // plan boundary-node-nozzle-wall-outlet-stability §2.11。亜音速の背圧アンカーは
        // 本カーネルの特性構成 bvar [ronew/Vn_exit] の質量流束が担う)。
        flow_float Pnew = Psb[ib];
        flow_float ronew = ro[ic];
        flow_float Uxnew = Ux[ic];
        flow_float Uynew = Uy[ic];
        flow_float Uznew = Uz[ic];
        flow_float Umagc = sqrt(Uxnew*Uxnew+Uynew*Uynew+Uznew*Uznew);
        flow_float roec = roe[ic];
        flow_float Umag_new;

        flow_float Un = (Uxnew*sxx+Uynew*syy+Uznew*szz)/sss;

        // SU2 BC_Outlet (亜音速静圧出口) の特性ベース構築 (CEulerSolver::BC_Outlet)。
        // naive な「P=Ps 規定 + ρ・u を内部から独立外挿」は ρ/u/P が相互不整合で、毎ステップ
        // 流束が規定状態と争い、壁∩出口コーナーで成長する圧力振動を生む (case36 で実証)。
        // 代わりに: P=P_exit を課し、ρ は内部等エントロピー上、法線速度は外向き Riemann 不変量、
        // 接線速度は内部外挿、として自己整合な1-incoming-characteristic 状態を作る。
        //
        // 逆流 (Un<0) も同じ静圧 P_exit アンカーで扱う (SU2 は方向分岐しない; Vn_exit<0 を許容し
        // upwind flux が incoming/outgoing 特性を捌く)。旧実装は逆流時に ambient 全圧 (Pt/Tt) の
        // stagnation 流入へ切替えていたが、これは SU2 の *inlet* TOTAL_CONDITIONS の構築で、剥離域
        // (壁∩出口コーナー) へ高 stagnation エンタルピ/全圧を注入し過加圧→発散させた (case36 層流で実証)。
        // さらに旧速度構築 (Uxnew=-Ux[ic]*nx) は物理的に妥当な速度ベクトルでなく、forward↔backflow の
        // ばたつきと相まって不整合ゴーストを生んでいた。逆流出口でも静圧のみ規定するのが整合的。
        // 注: thermalMethod==2 (TP) では γ が温度・組成依存。特性構成 (等エントロピー・Riemann 不変量・c_exit) の
        //     γ は **セル局所 γ_mix(T)** (gamma_cell) を使い、c_int も同じ γ で sqrt(γP/ρ) と再構成する。
        //     旧実装は c_int=sonic[ic] (真の TP 音速) と c_exit=sqrt(ga_cfg·P_exit/ρ) (config 定数 γ) を混用しており、
        //     P_int=P_exit でも c_int≠c_exit → Vn_exit = Un + 2(c_int−c_exit)/(γ−1) が数十 m/s ずれる。
        //     低マッハ出口 (case/48 Cabra coflow 3.5 m/s, γ_mix 1.30 vs config 1.35) では −65 m/s の偽逆流を
        //     毎ステップ課し、出口列に質量が蓄積して圧力上昇→全出口列が逆流→発散した (2026-09-05 特定)。
        //     超音速流出 (全量外挿) と CPG (γ 定数) はビット不変。
        const flow_float ga_loc = (thermalMethod == 2 && gamma_cell != nullptr) ? max(gamma_cell[ic], (flow_float)1.0001) : ga;
        const flow_float c_int = (thermalMethod == 2 && gamma_cell != nullptr)
                                 ? max(sqrt(ga_loc*P[ic]/max(ro[ic], (flow_float)1.0e-30)), (flow_float)1.0e-20)
                                 : max(sonic[ic], (flow_float)1.0e-20);
        const bool supersonic_out = (Un > (flow_float)0.0 && Umagc/max(sonic[ic], (flow_float)1.0e-20) >= (flow_float)1.0);
        if (supersonic_out) {
            // 超音速流出: P_back を課さず全量外挿 (SU2: Mach>=1 で V_outlet=V_domain)。
            Pnew = P[ic]; ronew = ro[ic];   // Uxnew/Uynew/Uznew は既に内部値
        } else {
            // 亜音速流出 / 局所逆流: 静圧 P_exit + 内部エントロピー + 外向き Riemann で統一構築。
            const flow_float P_exit = Psb[ib];
            const flow_float s_int  = P[ic]/pow(ro[ic], ga_loc);         // 内部エントロピー P/ρ^γ (γ=局所)
            const flow_float Rplus  = Un + (flow_float)2.0*c_int/(ga_loc-(flow_float)1.0); // 外向き Riemann 不変量
            ronew = pow(P_exit/s_int, (flow_float)1.0/ga_loc);            // ρ on interior isentrope
            const flow_float c_exit = sqrt(ga_loc*P_exit/ronew);
            const flow_float Vn_exit= Rplus - (flow_float)2.0*c_exit/(ga_loc-(flow_float)1.0); // 逆流なら <0 を許容 (クランプ無)
            const flow_float nx=sxx/sss, ny=syy/sss, nz=szz/sss;
            Pnew  = P_exit;
            Uxnew = Ux[ic] + (Vn_exit-Un)*nx;   // 接線=内部外挿, 法線のみ Vn_exit へ補正
            Uynew = Uy[ic] + (Vn_exit-Un)*ny;
            Uznew = Uz[ic] + (Vn_exit-Un)*nz;
            // ghost エネルギーを (ronew,Pnew,構築速度) と整合 (CPG)。TP は後段で T から再計算。
            roec = Pnew/(ga-(flow_float)1.0) + (flow_float)0.5*ronew*(Uxnew*Uxnew+Uynew*Uynew+Uznew*Uznew);
        }
        Umag_new = Umagc;   // uscale=1 (ghost 速度 = 構築値そのもの)

        // 速度方向スケール。静止場 (Umagc=0) では方向が未定義なので 0 とし、0/0=NaN を防ぐ
        // (Uxnew 等も 0 のため ghost 速度は 0 になる)。u=0 初期化での発散を回避。
        const flow_float uscale = (Umagc > (flow_float)1.0e-30) ? (Umag_new/Umagc) : (flow_float)0.0;

        ro[ig]    = ronew;
        P[ig]     = Pnew;
        Ux[ig]    = Uxnew*uscale;
        Uy[ig]    = Uynew*uscale;
        Uz[ig]    = Uznew*uscale;
        roUx[ig]  = ronew*Uxnew*uscale;
        roUy[ig]  = ronew*Uynew*uscale;
        roUz[ig]  = ronew*Uznew*uscale;

        if (supersonic_out) {
            // 超音速流出は roe/T/sonic も内部値をそのまま外挿する (ρ,P,U と同じ「全量外挿」)。
            // T=P/(ρR_mix) から e を再構築すると二相 (凝縮 g>0) の内部状態と不整合になる
            // (R_mix は全水分を気相として数え、e に液相項 g(R_wT−L) が無い) → 出口ノードだけ T が
            // 5 K 低下し S 上昇・g 跳ね上がりが出た (case/44 va2 平衡凝縮 run, 2026-08-18)。
            // 内部値コピーなら EOS の種類 (CPG/TP/二相) に依らず整合。
            roe[ig]   = roe[ic];
            sonic[ig] = sonic[ic];
            T[ig]     = T[ic];
        } else if (thermalMethod == 2) {
            // TP: ghost を (ρ=ronew, P=Pnew) と整合させる。T=P/(ρ R_mix)、以降の e/γ は共通 helper。
            //     定数 ga/cp は使わない。
            //     注: 凝縮 (g>0) の亜音速流出では R_mix/e が液相を無視するため近似 (要フォロー)。
            const flow_float ekg = 0.5f*(Ux[ig]*Ux[ig] + Uy[ig]*Uy[ig] + Uz[ig]*Uz[ig]);
            const double R   = mix ? thermo_R_mix(sp, nSpecies, Yc) : thermo_R_species(sp[0]);
            const flow_float Tg = (flow_float)((double)Pnew/((double)ronew*R));
            const GasStateAtT gs = thermo_state_at_T(thermalMethod, ga, cp, sp, roY, nSpecies, ic, Tg);
            roe[ig]   = ronew*(gs.e + ekg);
            sonic[ig] = sqrt(gs.gamma*Pnew/ronew);
            T[ig]     = Tg;
        } else {
            roe[ig]   = roec;
            sonic[ig] = sqrt(ga*Pnew/ronew);
            T[ig]     = Pnew*ga/(ronew*(ga-1.0)*cp);
        }
        Ht[ig]    = roe[ig]/ronew + Pnew/ronew;

        rob[ib]    = ronew;
        // Psb は入力専用 (規定背圧の面分布を保持) — ここでは書き換えない。
        Uxb[ib]    = Uxnew*uscale;
        Uyb[ib]    = Uynew*uscale;
        Uzb[ib]    = Uznew*uscale;
        roUxb[ib]  = ronew*Uxnew*uscale;
        roUyb[ib]  = ronew*Uynew*uscale;
        roUzb[ib]  = ronew*Uznew*uscale;
        roeb[ib]   = roe[ig];
        // Tsb は境界面の動的圧力 Pnew と整合させる (旧: 常に規定背圧から計算しており、
        // 超音速流出で ro=内部・P=背圧 の不整合な温度が勾配境界閉包 dT を汚染していた)。
        if (supersonic_out) {
            Tsb[ib] = T[ic];    // 全量外挿 (二相でも整合)
        } else if (thermalMethod == 2) {
            const double R = mix ? thermo_R_mix(sp, nSpecies, Yc) : thermo_R_species(sp[0]);
            Tsb[ib] = (flow_float)((double)Pnew/((double)rob[ib]*R));
        } else {
            Tsb[ib] = Pnew*ga/(rob[ib]*(ga-1.0)*cp);
        }
    }
};

void outlet_statPress_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    outlet_statPress_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> (
        cfg.gamma,
        cfg.cp,
        cfg.thermalMethod, thermo_species_device_ptr(),
        species_roY_device_ptr(), cfg.nSpecies,
        (cfg.thermalMethod == 2) ? var.c_d["gamma"] : nullptr,

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,
        bc.map_bplane_cell_d,
        bc.map_bplane_cell_ghst_d,

        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["P"],
        var.c_d["Ht"],
        var.c_d["sonic"],
        var.c_d["T"],

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]

    ) ;
}

__global__
void inlet_uniformVelocity_d
(
 // gas
 flow_float ga,
 flow_float cp,
 int thermalMethod, const SpeciesThermo* sp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,
 geom_int* bplane_cell,
 geom_int* bplane_cell_ghst,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,
flow_float* T,

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

 // 多成分 TP (M5): 入口組成 Y_s^in (device per-face array, [nSpecies])。
 // nSpecies>=2 かつ Yb!=nullptr のとき混合則 thermo を使う。単成分では nullptr。
 flow_float** Yb ,
 int nSpecies


)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        //geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];        // 内面セル
        geom_int  ig = bplane_cell_ghst[ib];

        // 速度は常に指定値 (Uxb/Uyb/Uzb) を課す。
        // 超音速流入: 全量を境界指定 (ρ=rob, P=Psb) で固定 (従来挙動)。
        // 亜音速流入: 1 特性が内部から出るため ρ,P は内面セルから外挿し、
        //             速度のみ指定する (velocity inlet)。
        // 判定は指定速度の大きさ vs 内面音速 (Mach = |U_in| / c_interior)。
        flow_float ek_b = 0.5f*(Uxb[ib]*Uxb[ib] + Uyb[ib]*Uyb[ib] + Uzb[ib]*Uzb[ib]);
        flow_float Umag = sqrt(2.0f*ek_b);
        bool supersonic = (Umag >= sonic[ic]);

        // 超音速流入: 全量を境界指定 (ρ=rob, P=Psb)。
        // 亜音速流入: 非反射 (Riemann 不変量) 速度入口。純粋な速度固定+圧力外挿は
        //   音響波を反射して不安定化するため、外向き特性 R+=un_i+2c_i/(γ-1) を内面から取り、
        //   指定法線速度 un_b と合わせて境界音速 c_b を決める。熱力学レベルは config エントロピー
        //   s=Psb/rob^γ で閉じ、ρ,P を一意化する (一様流では基準状態を厳密保存)。
        flow_float ro_use, P_use;
        if (supersonic) {
            ro_use = rob[ib];
            P_use  = Psb[ib];
        } else if (thermalMethod != 2) {
            geom_int ip = bplane_plane[ib];
            flow_float inv_ss = 1.0f/ss[ip];
            flow_float nx = sx[ip]*inv_ss, ny = sy[ip]*inv_ss, nz = sz[ip]*inv_ss;
            flow_float un_i = Ux[ic]*nx + Uy[ic]*ny + Uz[ic]*nz;   // 内面法線速度
            flow_float c_i  = sonic[ic];
            flow_float un_b = Uxb[ib]*nx + Uyb[ib]*ny + Uzb[ib]*nz; // 指定法線速度
            flow_float Rplus = un_i + 2.0f*c_i/(ga-1.0f);          // 内面からの外向き特性
            flow_float c_b   = 0.5f*(ga-1.0f)*(Rplus - un_b);
            flow_float c_floor = 1.0e-3f*c_i;
            if (c_b < c_floor) c_b = c_floor;                      // 非物理 (c_b<=0) 回避
            flow_float s_cfg = Psb[ib]/pow(rob[ib], ga);           // config エントロピー
            ro_use = pow(c_b*c_b/(ga*s_cfg), 1.0f/(ga-1.0f));
            P_use  = ro_use*c_b*c_b/ga;
        } else {
            // TP (thermalMethod 2) は等エントロピー閉包が非定数 γ で複雑なため、当面は
            // 圧力のみ内面外挿・等温で ρ を整合 (簡易)。必要なら別途精緻化。
            P_use  = P[ic];
            ro_use = rob[ib]*P_use/Psb[ib];
        }

        flow_float roe_b, sonic_b, Ts_b;
        if (thermalMethod == 2 && Yb != nullptr && nSpecies >= 2) {
            // 多成分 TP: 入口組成 Y_s^in と整合。T=P/(ρ R_mix(Y)),
            // roe=ρ(e_mix(Y,T)+ek), sonic=√(γ_mix(Y,T) P/ρ)。
            double Yin[THERMO_MAX_SPECIES];
            double ysum = 0.0;
            for (int s = 0; s < nSpecies; s++) {
                double y = (double)Yb[s][ib];
                if (y < 0.0) y = 0.0;
                Yin[s] = y;
                ysum += y;
            }
            const double yinv = 1.0/(ysum > 1.0e-30 ? ysum : 1.0e-30);
            for (int s = 0; s < nSpecies; s++) Yin[s] *= yinv;

            const double Rg  = thermo_R_mix(sp, nSpecies, Yin);
            const double Tg  = (double)P_use/((double)ro_use*Rg);
            const double e   = thermo_e_mix(sp, nSpecies, Yin, Tg);
            const double gmx = thermo_gamma_mix(sp, nSpecies, Yin, Tg);
            roe_b   = (flow_float)((double)ro_use*(e + (double)ek_b));
            sonic_b = (flow_float)sqrt(gmx*(double)P_use/(double)ro_use);
            Ts_b    = (flow_float)Tg;
        } else if (thermalMethod == 2) {
            // 単成分 TP: 指定 (ρ=ro_use, P=P_use, u) と整合。T=P/(ρ R_s), roe=ρ(e_NASA(T)+ek), sonic=√(γ_s P/ρ)。
            const double Rg = thermo_R_species(sp[0]);
            const double Tg = (double)P_use/((double)ro_use*Rg);
            const double e  = thermo_h_mass(sp[0], Tg) - Rg*Tg;
            const double cpv= thermo_cp_mass(sp[0], Tg);
            const double gmx= cpv/((cpv-Rg) > 1.0e-6 ? (cpv-Rg) : 1.0e-6);
            roe_b   = (flow_float)((double)ro_use*(e + (double)ek_b));
            sonic_b = (flow_float)sqrt(gmx*(double)P_use/(double)ro_use);
            Ts_b    = (flow_float)Tg;
        } else {
            roe_b   = P_use/(ga-1.0) + ro_use*ek_b;
            sonic_b = sqrt(ga*P_use/ro_use);
            Ts_b    = P_use*ga/(ro_use*(ga-1.0)*cp);
        }

        ro[ig]    = ro_use;
        roUx[ig]  = ro_use*Uxb[ib];
        roUy[ig]  = ro_use*Uyb[ib];
        roUz[ig]  = ro_use*Uzb[ib];
        roe[ig]   = roe_b;
        P[ig]     = P_use;
        Ux[ig]    = Uxb[ib];
        Uy[ig]    = Uyb[ib];
        Uz[ig]    = Uzb[ib];
        Ht[ig]    = roe_b/ro_use + P_use/ro_use;
        sonic[ig] = sonic_b;
        T[ig]     = Ts_b;

        // 境界流束 (convectiveFlux_boundary_d) が読む R 状態 (bvar) を更新。
        // 速度 Uxb/Uyb/Uzb は指定値のまま保持。亜音速では rob/Psb を内面値で上書きする
        // (outflow_d と同じ作法; 指定 ρ,P は亜音速では使わない)。
        rob[ib]    = ro_use;
        roUxb[ib]  = ro_use*Uxb[ib];
        roUyb[ib]  = ro_use*Uyb[ib];
        roUzb[ib]  = ro_use*Uzb[ib];
        roeb[ib]   = roe_b;
        Psb[ib]    = P_use;
        Tsb[ib]    = Ts_b;
        //Uxb[ib]    = Uxb[ib];  // 速度は指定値を保持
        //Uyb[ib]    = Uyb[ib];
        //Uzb[ib]    = Uzb[ib];

    }
};

void inlet_uniformVelocity_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    // 多成分 TP (M5): 入口組成 Y_s^in の device ポインタ配列を遅延構築 (bc.Yb_d にキャッシュ)。
    // 単成分 (nSpecies<2) では nullptr のままで、カーネルは従来の sp[0] 経路を使う。
    if (cfg.thermalMethod == 2 && cfg.nSpecies >= 2 && bc.Yb_d == nullptr) {
        std::vector<flow_float*> hYb(cfg.nSpecies);
        for (int s = 0; s < cfg.nSpecies; s++) {
            hYb[s] = bc.bvar_d["Y" + std::to_string(s)];
        }
        gpuErrchk( cudaMalloc((void**)&bc.Yb_d, cfg.nSpecies*sizeof(flow_float*)) );
        gpuErrchk( cudaMemcpy(bc.Yb_d, hYb.data(), cfg.nSpecies*sizeof(flow_float*), cudaMemcpyHostToDevice) );
    }

    inlet_uniformVelocity_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> (
        cfg.gamma,
        cfg.cp,
        cfg.thermalMethod, thermo_species_device_ptr(),

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

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
        var.c_d["T"],

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"],

        bc.Yb_d,        // 多成分 TP: 入口組成 (単成分では nullptr)
        cfg.nSpecies
    ) ;
}


__global__
void inlet_Pressure_d
(
 // gas
 flow_float ga,
 flow_float cp,
 int thermalMethod, const SpeciesThermo* sp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,
 geom_int* bplane_cell,
 geom_int* bplane_cell_ghst,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,
 flow_float* T   ,

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

 // 多成分 TP (M5): 入口組成 Y_s^in (device per-face array, [nSpecies])。
 // nSpecies>=2 かつ Yb!=nullptr のとき混合則 thermo を使う。単成分では nullptr。
 flow_float** Yb ,
 int nSpecies
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        flow_float Ux_c;
        flow_float Uy_c;
        flow_float Uz_c;
        flow_float Umag_c;

        flow_float mach_c;
        flow_float mach_b;
        flow_float mach_new;

        //flow_float T_c;
        //flow_float P_c;

        flow_float Pt_b;
        flow_float Tt_b;

        flow_float Ps_new;
        flow_float Ts_new;
        flow_float ro_new;
        flow_float sonic_new;

        // TP 亜音速入口の mach ブレンド係数 (静圧参照 mach_b と速度参照 mach_c)。
        // [2026-08-16] rf=0.5 → 0.0 (速度参照のみ = CPG 分岐と同じ閉包)。静圧参照を混ぜると
        // 大径低速入口 (M≈0.01–0.03) で P[ic]↔mach_b の帰還が振動し、case/42 M4.2 (M_in 0.03) では
        // 残差床 rms_ro 3e-2・入口 P/Pt 1.03 振動、M6 (M_in 0.011) では 16k step 以降に指数成長して
        // P/Pt 6 まで発散した (run_0040)。CPG 分岐 (速度参照のみ) は同条件で rms_ro 1e-9 まで収束する。
        flow_float rf=0.0;

        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        Pt_b = Ptb[ib];
        Tt_b = Ttb[ib];

        Ux_c = Ux[ic];
        Uy_c = Uy[ic];
        Uz_c = Uz[ic];

        flow_float sxx = sx[ip];
        flow_float syy = sy[ip];
        flow_float szz = sz[ip];
        flow_float sss = ss[ip];


        //Umag_c = sqrt(Ux_c*Ux_c + Uy_c*Uy_c + Uz_c*Uz_c);
        flow_float Un_c = (Ux_c*sxx + Uy_c*syy + Uz_c*szz)/sss;
        mach_c = Un_c/sonic[ic];

        // 多成分 TP: 入口指定組成 Yin を正規化 (nSpecies>=2 のとき混合則 thermo を使う)。
        const bool mix = (thermalMethod == 2 && Yb != nullptr && nSpecies >= 2);
        double Yin[THERMO_MAX_SPECIES];
        if (mix) {
            double ysum = 0.0;
            for (int s = 0; s < nSpecies; s++) { double y=(double)Yb[s][ib]; if(y<0.0)y=0.0; Yin[s]=y; ysum+=y; }
            const double yinv = 1.0/(ysum > 1.0e-30 ? ysum : 1.0e-30);
            for (int s = 0; s < nSpecies; s++) Yin[s] *= yinv;
        }

        flow_float Ux_new;
        flow_float Uy_new;
        flow_float Uz_new;

        //P_c = P[ic];
        //T_c = T[ic];


        if ( Un_c < 0) { 

            //mach_b = sqrt((pow((P_c/Pt_b),-(ga-1.0)/ga) -1.0)*2.0/(ga-1.0));

            //flow_float mach_new = rf*mach_b + (1.0-rf)*mach_c;

            //Ps_new = P_c;
            ////Ts_new = Tt_b/(1.0+0.5*(ga-1.0)*mach_b*mach_b);
            //Ts_new = Tt_b/(1.0+0.5*(ga-1.0)*mach_new*mach_new);
            //sonic_new = sqrt((ga-1.0)*cp*Ts_new);
            //ro_new = ga*Ps_new/((ga-1.0)*cp*Ts_new);

            ////Ux_new = -mach_b*sonic_new*sxx/sss;
            ////Uy_new = -mach_b*sonic_new*syy/sss;
            ////Uz_new = -mach_b*sonic_new*szz/sss;
            //Ux_new = -mach_new*sonic_new*sxx/sss;
            //Uy_new = -mach_new*sonic_new*syy/sss;
            //Uz_new = -mach_new*sonic_new*szz/sss;

            Ux_new = Un_c*sxx/sss;
            Uy_new = Un_c*syy/sss;
            Uz_new = Un_c*szz/sss;

            if (thermalMethod == 2) {
                // TP 亜音速圧力入口 (inlet_Pressure_dir_d と同方式の緩和)。
                // 内部静圧 P[ic] から mach_b、内部速度から mach_c を作り、mach_new=rf*mach_b+(1-rf)*mach_c
                // で緩和する (rf=0.5)。純静圧参照 (rf=1) は不足減衰で振動、純速度参照 (rf=0) は静的状態を
                // 過拘束し P>Pt の反射を生むため、両者の中間で安定化。blend 後マッハで全状態から再構成。
                // 多成分は混合則、単成分は sp[0]。
                // P[ic] > Pt は等エントロピー計算で M²<0 → NaN になるため Pt*(1-ε) でクランプ。
                // 起動過渡期の圧力反射で静圧が全圧を超えるときに発生し、mach_bd→0 として処理する。
                const double Psc = fmin((double)P[ic], (double)Pt_b * (1.0 - 1.0e-6));
                double Ts0_d, ro0_d, um0_d, a0;
                if (mix) {
                    thermo_isentropic_from_total_Ps_mix(sp, nSpecies, Yin, Pt_b, Tt_b, Psc, &Ts0_d, &ro0_d, &um0_d);
                    a0 = sqrt(thermo_gamma_mix(sp, nSpecies, Yin, Ts0_d)*thermo_R_mix(sp, nSpecies, Yin)*Ts0_d);
                } else {
                    thermo_isentropic_from_total_Ps_single(sp, Pt_b, Tt_b, Psc, &Ts0_d, &ro0_d, &um0_d);
                    const double Rs = thermo_R_species(sp[0]), cps = thermo_cp_mass(sp[0], Ts0_d);
                    a0 = sqrt((cps/((cps-Rs) > 1.0e-6 ? (cps-Rs) : 1.0e-6))*Rs*Ts0_d);
                }
                const double mach_bd  = um0_d/(a0 > 1.0e-6 ? a0 : 1.0e-6);
                const double mach_new = (double)rf*mach_bd + (1.0-(double)rf)*fabs((double)mach_c);
                double Ts_d, Ps_d, ro_d, um_d;
                if (mix)
                    thermo_isentropic_from_total_mix(sp, nSpecies, Yin, Pt_b, Tt_b, mach_new, &Ts_d, &Ps_d, &ro_d, &um_d);
                else
                    thermo_isentropic_from_total_single(sp, Pt_b, Tt_b, mach_new, &Ts_d, &Ps_d, &ro_d, &um_d);
                Ps_new = (flow_float)Ps_d; Ts_new = (flow_float)Ts_d; ro_new = (flow_float)ro_d;
                const double gmx = mix ? thermo_gamma_mix(sp, nSpecies, Yin, Ts_d)
                                       : (thermo_cp_mass(sp[0],Ts_d)/((thermo_cp_mass(sp[0],Ts_d)-thermo_R_species(sp[0])) > 1.0e-6 ? (thermo_cp_mass(sp[0],Ts_d)-thermo_R_species(sp[0])) : 1.0e-6));
                sonic_new = (flow_float)sqrt(gmx*(double)Ps_new/(double)ro_new);
                // 流入方向は境界法線の内向き (-n)、大きさ |u| は blend 後マッハから。
                const flow_float invs = 1.0f/sss;
                Ux_new = -(flow_float)um_d * sxx*invs;
                Uy_new = -(flow_float)um_d * syy*invs;
                Uz_new = -(flow_float)um_d * szz*invs;
            } else {
                Ps_new = Pt_b/pow(1.0+0.5*(ga-1.0)*mach_c*mach_c, ga/(ga-1.0));
                Ts_new = Tt_b/(1.0+0.5*(ga-1.0)*mach_c*mach_c);
                sonic_new = sqrt((ga-1.0)*cp*Ts_new);
                ro_new = ga*Ps_new/((ga-1.0)*cp*Ts_new);
            }

        } else { // reverse-flow detected at pressure inlet: clamp ghost to stagnation
            Ux_new = 0.0;
            Uy_new = 0.0;
            Uz_new = 0.0;

            Ps_new = Pt_b;
            Ts_new = Tt_b;
            if (thermalMethod == 2) {
                flow_float R = (flow_float)(mix ? thermo_R_mix(sp, nSpecies, Yin) : thermo_R_species(sp[0]));
                ro_new = Ps_new/(R*Ts_new);
                sonic_new = sqrt(ga*Ps_new/ro_new);
            } else {
                sonic_new = sqrt((ga-1.0)*cp*Ts_new);
                ro_new = ga*Ps_new/((ga-1.0)*cp*Ts_new);
            }
        }

        flow_float ek = 0.5*(Ux_new*Ux_new +Uy_new*Uy_new +Uz_new*Uz_new);
        flow_float roe_new;
        if (thermalMethod == 2) {
            // TP: roe = ρ(e(Ts) + ek), sonic = √(γ Ps/ρ)。多成分は混合則、単成分は sp[0]。
            const double Tsd = (double)Ts_new;
            double e_in, gmx;
            if (mix) {
                e_in = thermo_e_mix(sp, nSpecies, Yin, Tsd);
                gmx  = thermo_gamma_mix(sp, nSpecies, Yin, Tsd);
            } else {
                const double R   = thermo_R_species(sp[0]);
                const double cpv = thermo_cp_mass(sp[0], Tsd);
                e_in = thermo_h_mass(sp[0], Tsd) - R*Tsd;
                gmx  = cpv/((cpv-R) > 1.0e-6 ? (cpv-R) : 1.0e-6);
            }
            roe_new   = (flow_float)((double)ro_new*(e_in + (double)ek));
            sonic_new = (flow_float)sqrt(gmx*(double)Ps_new/(double)ro_new);
        } else {
            roe_new   = Ps_new/(ga-1.0) + ro_new*ek;
        }

        ro[ig]    = ro_new;
        Ux[ig]    = Ux_new;
        Uy[ig]    = Uy_new;
        Uz[ig]    = Uz_new;
        roUx[ig]  = ro_new*Ux_new;
        roUy[ig]  = ro_new*Uy_new;
        roUz[ig]  = ro_new*Uz_new;
        roe[ig]   = roe_new;
        P[ig]     = Ps_new;
        Ht[ig]    = roe[ig]/ro_new + Ps_new/ro_new;
        sonic[ig] = sonic_new;
        T[ig]     = Ts_new;

        rob[ib]    = ro_new;
        Uxb[ib]    = Ux_new;
        Uyb[ib]    = Uy_new;
        Uzb[ib]    = Uz_new;
        roUxb[ib]  = ro_new*Ux_new;
        roUyb[ib]  = ro_new*Uy_new;
        roUzb[ib]  = ro_new*Uz_new;
        roeb[ib]   = roe_new;
        Psb[ib]    = Ps_new;
        Tsb[ib]    = Ts_new;

        //Htb[ib]    = roeb[ib]/ro_new + Ps_new/ro_new;
        //sonicb[ib] = sqrt(ga*Ps_new/ro_new);
 
    }
};

void inlet_Pressure_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    // 多成分 TP (M5): 入口組成 Y_s^in の device ポインタ配列を遅延構築 (bc.Yb_d にキャッシュ)。
    // 単成分 (nSpecies<2) では nullptr のままで、カーネルは従来の sp[0] 経路を使う。
    if (cfg.thermalMethod == 2 && cfg.nSpecies >= 2 && bc.Yb_d == nullptr) {
        std::vector<flow_float*> hYb(cfg.nSpecies);
        for (int s = 0; s < cfg.nSpecies; s++) {
            hYb[s] = bc.bvar_d["Y" + std::to_string(s)];
        }
        gpuErrchk( cudaMalloc((void**)&bc.Yb_d, cfg.nSpecies*sizeof(flow_float*)) );
        gpuErrchk( cudaMemcpy(bc.Yb_d, hYb.data(), cfg.nSpecies*sizeof(flow_float*), cudaMemcpyHostToDevice) );
    }

    inlet_Pressure_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> (
        cfg.gamma,
        cfg.cp,
        cfg.thermalMethod, thermo_species_device_ptr(),

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

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

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"],

        bc.Yb_d,        // 多成分 TP: 入口組成 (単成分では nullptr)
        cfg.nSpecies
    ) ;
}

__global__
void inlet_Pressure_dir_d
(
 // gas
 flow_float ga,
 flow_float cp,
 int thermalMethod, const SpeciesThermo* sp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,
 geom_int* bplane_cell,
 geom_int* bplane_cell_ghst,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,
 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,
flow_float* T,

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
 flow_float* Psb

)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        flow_float Ux_c;
        flow_float Uy_c;
        flow_float Uz_c;
        flow_float Umag_c;

        flow_float mach_c;
        flow_float mach_b;
        flow_float mach_new;

        flow_float T_c;
        flow_float P_c;

        flow_float Pt_b;
        flow_float Tt_b;

        flow_float Ps_new;
        flow_float Ts_new;
        flow_float ro_new;
        flow_float sonic_new;

        flow_float rf=0.5;

        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        Pt_b = Ptb[ib];
        Tt_b = Ttb[ib];

        flow_float Ubmag = sqrt(Uxb[ib]*Uxb[ib] + Uyb[ib]*Uyb[ib] + Uzb[ib]*Uzb[ib]);

        flow_float Unx_b = Uxb[ib]/Ubmag;
        flow_float Uny_b = Uyb[ib]/Ubmag;
        flow_float Unz_b = Uzb[ib]/Ubmag;

        Ux_c = Ux[ic];
        Uy_c = Uy[ic];
        Uz_c = Uz[ic];

        Umag_c = sqrt(Ux_c*Ux_c + Uy_c*Uy_c + Uz_c*Uz_c);
        mach_c = Umag_c/sonic[ic];

        T_c = T[ic];
        P_c = P[ic];

        flow_float Umag_new_b;
        if (thermalMethod == 2) {
            // TP: 全状態 (Pt,Tt) と外挿静圧 Ps=P_c から NASA 等エントロピーで静的状態を反転。
            //     CPG 同様 cell マッハとブレンドするため、blend 後の静圧で再評価する。
            double Ts0_d, ro0_d, um0_d;
            thermo_isentropic_from_total_Ps_single(sp, (double)Pt_b, (double)Tt_b, (double)P_c,
                                                   &Ts0_d, &ro0_d, &um0_d);
            const double a0 = sqrt(thermo_cp_mass(sp[0], Ts0_d)
                               /((thermo_cp_mass(sp[0], Ts0_d)-thermo_R_species(sp[0]))>1e-6
                                 ? (thermo_cp_mass(sp[0], Ts0_d)-thermo_R_species(sp[0])):1e-6)
                               * thermo_R_species(sp[0]) * Ts0_d);
            mach_b = (flow_float)(um0_d/(a0 > 1.0e-6 ? a0 : 1.0e-6));
            mach_new = rf*mach_b + (1.0-rf)*mach_c;
            // blend 後マッハに対応する静温・速度を全状態から再構成。
            double Ts_d, Ps_d, ro_d, um_d;
            thermo_isentropic_from_total_single(sp, (double)Pt_b, (double)Tt_b, (double)mach_new,
                                                &Ts_d, &Ps_d, &ro_d, &um_d);
            Ps_new = (flow_float)Ps_d; Ts_new = (flow_float)Ts_d; ro_new = (flow_float)ro_d;
            Umag_new_b = (flow_float)um_d;
            const double cpv = thermo_cp_mass(sp[0], Ts_d), Rg = thermo_R_species(sp[0]);
            const double gmx = cpv/((cpv-Rg) > 1.0e-6 ? (cpv-Rg) : 1.0e-6);
            sonic_new = (flow_float)sqrt(gmx*(double)Ps_new/(double)ro_new);
        } else {
            mach_b = sqrt((pow((P_c/Pt_b),-(ga-1.0)/ga) -1.0)*2.0/(ga-1.0));
            mach_new = rf*mach_b + (1.0-rf)*mach_c;
            Ps_new = P_c;
            Ts_new = Tt_b/(1.0+0.5*(ga-1.0)*mach_new*mach_new);
            sonic_new = sqrt((ga-1.0)*cp*Ts_new);
            ro_new = ga*Ps_new/((ga-1.0)*cp*Ts_new);
            Umag_new_b = mach_new*sonic_new;
        }

        flow_float Ux_new = Umag_new_b*Unx_b;
        flow_float Uy_new = Umag_new_b*Uny_b;
        flow_float Uz_new = Umag_new_b*Unz_b;

        flow_float ek_new = 0.5f*(Ux_new*Ux_new + Uy_new*Uy_new + Uz_new*Uz_new);
        flow_float roe_new;
        if (thermalMethod == 2) {
            const double Rg = thermo_R_species(sp[0]);
            const double e  = thermo_h_mass(sp[0], (double)Ts_new) - Rg*(double)Ts_new;
            roe_new = (flow_float)((double)ro_new*(e + (double)ek_new));
            // sonic_new は TP 分岐で γ_mix 済み。
        } else {
            roe_new = Ps_new/(ga-1.0) + ro_new*ek_new;
        }

        ro[ig]    = ro_new;
        Ux[ig]    = Ux_new;
        Uy[ig]    = Uy_new;
        Uz[ig]    = Uz_new;
        roUx[ig]  = ro_new*Ux_new;
        roUy[ig]  = ro_new*Uy_new;
        roUz[ig]  = ro_new*Uz_new;
        roe[ig]   = roe_new;
        P[ig]     = Ps_new;
        Ht[ig]    = roe[ig]/ro_new + Ps_new/ro_new;
        sonic[ig] = (thermalMethod == 2) ? sonic_new : sqrt(ga*Ps_new/ro_new);
        T[ig]     = Ts_new;

        rob[ib]    = ro_new;
        Uxb[ib]    = Ux_new;
        Uyb[ib]    = Uy_new;
        Uzb[ib]    = Uz_new;
        roUxb[ib]  = ro_new*Ux_new;
        roUyb[ib]  = ro_new*Uy_new;
        roUzb[ib]  = ro_new*Uz_new;
        roeb[ib]   = roe_new;
        Psb[ib]    = Ps_new;
        Tsb[ib]    = Ts_new;
    }
};

void inlet_Pressure_dir_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    inlet_Pressure_dir_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>> (
        cfg.gamma,
        cfg.cp,
        cfg.thermalMethod, thermo_species_device_ptr(),

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  bc.map_bplane_cell_d, bc.map_bplane_cell_ghst_d,
        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"]  , 
        var.c_d["Ux"]  , 
        var.c_d["Uy"]  , 
        var.c_d["Uz"]  , 
        var.c_d["P"]  , 
        var.c_d["Ht"]  , 
        var.c_d["sonic"]  , 
        var.c_d["T"]  , 

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]

    ) ;
}


__global__ 
void outflow_d
(
 // gas properties
 flow_float ga,
 flow_float cp,
 int thermalMethod, const SpeciesThermo* sp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,
 geom_int* bplane_cell,
 geom_int* bplane_cell_ghst,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic, 
flow_float* T,

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
 flow_float* Psb 
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        flow_float Pnew ;
        flow_float ronew ;
        flow_float Uxnew = Ux[ic];
        flow_float Uynew = Uy[ic];
        flow_float Uznew = Uz[ic];
        flow_float Umagc = sqrt(Uxnew*Uxnew+Uynew*Uynew+Uznew*Uznew);
        flow_float Umag_new;

        flow_float sxx = sx[ip];
        flow_float syy = sy[ip];
        flow_float szz = sz[ip];
        flow_float sss = ss[ip];

        flow_float Un = (Ux[ic]*sxx+Uy[ic]*syy+Uz[ic]*szz)/sss;

        if (Un <= 0) {
            Uxnew = -Uxnew*sxx/sss;
            Uynew = -Uynew*syy/sss;
            Uznew = -Uznew*szz/sss;

            Umag_new = sqrt(Uxnew*Uxnew +Uynew*Uynew +Uznew*Uznew);

            flow_float mach_new = Umag_new/sonic[ic];

            flow_float Pt_b = Ptb[ib];
            flow_float Tt_b = Ttb[ib];
            Pnew = Pt_b/pow(1.0+0.5*(ga-1.0)*mach_new*mach_new, ga/(ga-1.0));
            flow_float Tnew = Tt_b/(1.0+0.5*(ga-1.0)*mach_new*mach_new);
            ronew = ga*Pnew/((ga-1.0)*cp*Tnew);
        }


        ro[ig]   = ro[ic];
        P[ig]    = P[ic];
        Ux[ig]   = Ux[ic];
        Uy[ig]   = Uy[ic];
        Uz[ig]   = Uz[ic];
        roUx[ig] = roUx[ic];
        roUy[ig] = roUy[ic];
        roUz[ig] = roUz[ic];
        roe[ig]  = roe[ic];

        Ht[ig]   = (roe[ig] + P[ig])/ro[ig];   // 純代数 (TP でも整合)
        T[ig]    = T[ic];

        rob[ib]   = ro[ic];
        Psb[ib]   = P[ic];
        Uxb[ib]   = Ux[ic];
        Uyb[ib]   = Uy[ic];
        Uzb[ib]   = Uz[ic];
        roUxb[ib] = roUx[ic];
        roUyb[ib] = roUy[ic];
        roUzb[ib] = roUz[ic];
        roeb[ib]  = roe[ic];
        if (thermalMethod == 2) {
            // 外挿: 熱力学状態は内部コピー。sonic/Tsb を NASA 整合に。
            sonic[ig] = sonic[ic];
            Tsb[ib]   = T[ic];
        } else {
            sonic[ig] = sqrt(ga*P[ig]/ro[ig]);
            Tsb[ib]   = Psb[ib]*ga/(rob[ib]*(ga-1.0)*cp);
        }
    }
};

void outflow_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    outflow_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>(
        cfg.gamma,
        cfg.cp,
        cfg.thermalMethod, thermo_species_device_ptr(),

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  
        bc.map_bplane_cell_d,  
        bc.map_bplane_cell_ghst_d,

        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["P"], 
        var.c_d["Ht"], 
        var.c_d["sonic"],
        var.c_d["T"],

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]
    );
}

__global__ 
void periodic_d 
( 
 // gas properties
 flow_float ga,
 flow_float cp,

 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  

 geom_int* bplane_partnerCellID,  
 geom_float dtheta,

 geom_float* x   ,  geom_float* y   ,  geom_float* z ,
 geom_float* sx  ,  geom_float* sy  ,  geom_float* sz , geom_float* ss,

 // variables
 flow_float* ro   ,
 flow_float* roUx ,
 flow_float* roUy ,
 flow_float* roUz ,
 flow_float* roe ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,
 flow_float* P   ,
 flow_float* T   ,
 flow_float* Ht  ,
 flow_float* sonic,

 // 多成分: 化学種 ρY_s (flow_float*[nSpecies], 単成分時 nullptr)
 flow_float** roY,
 int nSpecies,

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
flow_float* Psb 
)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ig         = bplane_cell_ghst[ib];
        // periodic: ghost には自セルでなく **partner セル** の状態を写す (zero-gradient は誤り)。
        // 速度ベクトルのみ回転 periodic (dtheta) で座標変換。Cartesian は dtheta=0 で恒等。
        geom_int  ic         = bplane_partnerCellID[ib];

        ro[ig]   = ro[ic];
        P[ig]    = P[ic];
        Ux[ig]   = Ux[ic];
        Uy[ig]   = Uy[ic]*cos(-dtheta) -Uz[ic]*sin(-dtheta);
        Uz[ig]   = Uy[ic]*sin(-dtheta) +Uz[ic]*cos(-dtheta);
        roUx[ig] = roUx[ic];
        roUy[ig] = roUy[ic]*cos(-dtheta) -roUz[ic]*sin(-dtheta);
        roUz[ig] = roUy[ic]*sin(-dtheta) +roUz[ic]*cos(-dtheta);
        roe[ig]  = roe[ic];
        // T/Ht/sonic は partner から直接コピー (多成分 TP 整合; 旧 CPG 再計算 P·ga/(ρ(ga-1)cp) は誤り)。
        T[ig]    = T[ic];
        Ht[ig]   = Ht[ic];
        sonic[ig]= sonic[ic];
        // 多成分: 化学種 ρY_s も partner からコピー (未コピーだと seam で組成・R_mix が壊れる)。単成分は roY=nullptr。
        if (roY != nullptr) { for (int s = 0; s < nSpecies; ++s) roY[s][ig] = roY[s][ic]; }

        rob[ib]   = ro[ic];
        Psb[ib]   = P[ic];
        Uxb[ib]   = Ux[ic];
        Uyb[ib]   = Uy[ic]*cos(-dtheta) -Uz[ic]*sin(-dtheta);
        Uzb[ib]   = Uy[ic]*sin(-dtheta) +Uz[ic]*cos(-dtheta);
        roUxb[ib] = roUx[ic];
        roUyb[ib] = roUy[ic]*cos(-dtheta) -roUz[ic]*sin(-dtheta);
        roUzb[ib] = roUy[ic]*sin(-dtheta) +roUz[ic]*cos(-dtheta);
        roeb[ib]  = roe[ic];
        Tsb[ib]   = T[ic];
    }
};

void periodic_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    periodic_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>( 
        cfg.gamma,
        cfg.cp,

        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  
        bc.map_bplane_cell_d,  
        bc.map_bplane_cell_ghst_d,

        bc.bint_d["partnerCellID"],
        bc.inputFloats["dtheta"],   // 回転周期の回転角 [rad] (float)。inputInts から読むと常に 0 で回転周期が無効化される。
                                    // Cartesian (type=0) はキー無し → 0.0 で恒等 (従来どおり無害)。

        var.p_d["x"],  
        var.p_d["y"],  
        var.p_d["z"],  
        var.p_d["sx"],  
        var.p_d["sy"],  
        var.p_d["sz"],  
        var.p_d["ss"],  

        var.c_d["ro"] ,
        var.c_d["roUx"] ,
        var.c_d["roUy"] ,
        var.c_d["roUz"] ,
        var.c_d["roe"] ,
        var.c_d["Ux"] ,
        var.c_d["Uy"] ,
        var.c_d["Uz"] ,
        var.c_d["P"],
        var.c_d["T"],
        var.c_d["Ht"],
        var.c_d["sonic"],

        species_roY_device_ptr(),
        cfg.nSpecies,

        bc.bvar_d["ro"],
        bc.bvar_d["roUx"],
        bc.bvar_d["roUy"],
        bc.bvar_d["roUz"],
        bc.bvar_d["roe"],
        bc.bvar_d["Ux"],
        bc.bvar_d["Uy"],
        bc.bvar_d["Uz"],
        bc.bvar_d["Tt"],
        bc.bvar_d["Pt"],
        bc.bvar_d["Ts"],
        bc.bvar_d["Ps"]
    );
}



__global__ 
void copyBcondsGradient_d 
( 
 // mesh structure
 geom_int nb,
 geom_int* bplane_plane,  
 geom_int* bplane_cell,  
 geom_int* bplane_cell_ghst,  

 // variables
 flow_float* dTdx,
 flow_float* dTdy,
 flow_float* dTdz,

 flow_float* dHtdx,
 flow_float* dHtdy,
 flow_float* dHtdz,

 flow_float* drodx,
 flow_float* drody,
 flow_float* drodz,

 flow_float* dUxdx,
 flow_float* dUydx,
 flow_float* dUzdx,

 flow_float* dUxdy,
 flow_float* dUydy,
 flow_float* dUzdy,
 
 flow_float* dUxdz,
 flow_float* dUydz,
 flow_float* dUzdz,
 
 flow_float* dPdx,
 flow_float* dPdy,
 flow_float* dPdz

)
{
    geom_int ib  = blockDim.x*blockIdx.x + threadIdx.x;

    if (ib < nb) {
        geom_int  ip = bplane_plane[ib];
        geom_int  ic = bplane_cell[ib];
        geom_int  ig = bplane_cell_ghst[ib];

        dTdx[ig] = dTdx[ic];
        dTdy[ig] = dTdy[ic];
        dTdz[ig] = dTdz[ic];

        dHtdx[ig] = dHtdx[ic];
        dHtdy[ig] = dHtdy[ic];
        dHtdz[ig] = dHtdz[ic];

        drodx[ig] = drodx[ic];
        drody[ig] = drody[ic];
        drodz[ig] = drodz[ic];

        dUxdx[ig] = dUxdx[ic];
        dUydx[ig] = dUydx[ic];
        dUzdx[ig] = dUzdx[ic];

        dUxdy[ig] = dUxdy[ic];
        dUydy[ig] = dUydy[ic];
        dUzdy[ig] = dUzdy[ic];
 
        dUxdz[ig] = dUxdz[ic];
        dUydz[ig] = dUydz[ic];
        dUzdz[ig] = dUzdz[ic];

        dPdx[ig] = dPdx[ic];
        dPdy[ig] = dPdy[ic];
        dPdz[ig] = dPdz[ic];
    }
};

void copyBcondsGradient_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , bcond& bc , mesh& msh , variables& var , matrix& mat_p)
{
    copyBcondsGradient_d<<<cuda_cfg.dimGrid_bplane , cuda_cfg.dimBlock>>>( 
        bc.iPlanes.size(),
        bc.map_bplane_plane_d,  
        bc.map_bplane_cell_d,  
        bc.map_bplane_cell_ghst_d,

        var.c_d["dTdx"] ,
        var.c_d["dTdy"] ,
        var.c_d["dTdz"] ,
        var.c_d["dHtdx"] ,
        var.c_d["dHtdy"] ,
        var.c_d["dHtdz"] ,
        var.c_d["drodx"] ,
        var.c_d["drody"] ,
        var.c_d["drodz"] ,
        var.c_d["dUxdx"] ,
        var.c_d["dUydx"] ,
        var.c_d["dUzdx"] ,
        var.c_d["dUxdy"] ,
        var.c_d["dUydy"] ,
        var.c_d["dUzdy"] ,
        var.c_d["dUxdz"] ,
        var.c_d["dUydz"] ,
        var.c_d["dUzdz"] ,
        var.c_d["dPdx"] ,
        var.c_d["dPdy"] ,
        var.c_d["dPdz"] 
    );
}



