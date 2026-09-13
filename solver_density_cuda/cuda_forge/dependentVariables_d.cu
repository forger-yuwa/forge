#include "dependentVariables_d.cuh"
#include "thermo_d.cuh"
#include "speciesTransport_d.cuh"  // species_roY_device_ptr() (多成分 thermalMethod==2)
#include "condensationTransport_d.cuh"  // cond_rog_device_ptr() (二相 EOS)
#include "condensationEOS_d.cuh"        // cond_T_from_e_onetemp (一温度 二相 温度反転)

// 温度反転のクランプ範囲 (NASA-9 の有効域より広めに取り, 範囲外は外挿)
#define DEPVAR_TMIN 50.0
// CPG 二相の温度反転が収束しなかったセル数 (診断; plans/accepted/condensation-air.md §4.1)。wrapper が毎ステップ読み出して警告する。
__device__ unsigned int g_condTinvFail = 0u;
#define DEPVAR_TMAX 6000.0

__global__ void dependentVariables_d
(
 // gas properties
 int thermalMethod ,
 flow_float gamma , flow_float cp ,

 // EOS 正値化フロア (config 化。既定: pMin=1.0, roMin=1e-4, tMin=1e-4)
 flow_float pMin , flow_float roMin , flow_float tMin ,

 // thermally-perfect (thermalMethod==2) 用化学種データ
 const SpeciesThermo* sp , int nSpecies , flow_float** roY ,
 // ハイブリッド温度反転 (thermoFloat==1): float ミラー係数。0 のとき未使用。
 const SpeciesThermoF* spf , int thermoFloat ,

 // 非平衡凝縮 (一温度 二相 EOS)。condensation==0 のとき従来経路 (ビット不変)。
 int condensation , int nCondSpecies , flow_float** rog ,
 int condGasSpecies , int condModel ,   // carrier+condensible: 凝縮気相種 index / モデル(0:N2,1:H2O)
 int condEquilibrium ,                  // 2: EOS 拘束形平衡 ((T,g) 同時反転 → rog 射影)。0/1 は従来経路
 int condSonicModel ,                   // 1: 二相 frozen 音速/γ (TP 分岐, g>0 セル)。0: 旧 (全蒸気 √(γ_mix R_mix T))
 CondPropOpts condOpts ,
 int condFloat , CondTablesF condTb ,   // 凝縮 float 経路 (二相ハイブリッド反転; plans/active/condensation-float-speedup.md §4.2-5)                 // σ 倍率・N2 低温物性・CPG carrier の Y_w (plans/accepted/condensation-air.md)

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

 flow_float* P   ,
 flow_float* Ht  ,
 flow_float* sonic,
 flow_float* k   ,
 flow_float* omega,
 flow_float* T   ,
 flow_float* Ux  ,
 flow_float* Uy  ,
 flow_float* Uz  ,

 // thermally-perfect 用に毎ステップ更新する per-cell 物性
 flow_float* gam_array ,
 flow_float* cp_array ,
 flow_float* Rmix_array   // M6: 混合比気体定数 R[ic] (SLAU 面エンタルピーが毎面再計算するのを回避)
)
{
    geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;

    flow_float ek;
    flow_float intE;

    flow_float ro_temp;
    flow_float T_temp;
    flow_float P_temp;

    if (ic < nCells_all) {
        // 密度・圧力フロア: 膨張領域で非物理的な ro→0, P→0 が生じても速度爆発を防ぐ。
        // 既定 ro_min=1e-4 kg/m³, P_min=1.0 Pa は大気スケール想定。無次元・低圧ケースでは
        // config (physProp.pMin/roMin/tMin) で下げる (既定値は従来ハードコードと同一)。
        ro_temp = max(ro[ic], roMin);

        Ux[ic] = roUx[ic]/ro_temp;
        Uy[ic] = roUy[ic]/ro_temp;
        Uz[ic] = roUz[ic]/ro_temp;

        ek = 0.5f*(Ux[ic]*Ux[ic] +Uy[ic]*Uy[ic] +Uz[ic]*Uz[ic]);
        intE =(roe[ic]/ro_temp -ek);

        if (thermalMethod == 2) {
            // ---- 多成分 thermally-perfect gas (NASA-9) ----
            // 内部計算は全て double。組成 Y を構築 (nSpecies==1 は Y={1})。
            double Y[THERMO_MAX_SPECIES];
            float  Yf[THERMO_MAX_SPECIES];   // thermoFloat 用 (float で組んで double へ昇格: DP 除算を避ける)
            // 凝縮 ON でも condFloat=1 なら g≈0 セルはハイブリッド (組成 Y の float 構築も同条件で切替; plan condensation-float-speedup §4.2-6)
            const bool useHybrid = (thermoFloat != 0 && spf != nullptr && (condensation == 0 || condFloat != 0));
            if (nSpecies <= 1 || roY == nullptr) {
                Y[0] = 1.0; Yf[0] = 1.0f;
            } else if (useHybrid) {
                const float inv_ro = 1.0f/ro_temp;
                float ysum = 0.0f;
                for (int s=0;s<nSpecies;s++){
                    float y = roY[s][ic]*inv_ro;
                    if (y < 0.0f) y = 0.0f;
                    Yf[s] = y; ysum += y;
                }
                const float inv = 1.0f/(ysum > 1.0e-30f ? ysum : 1.0e-30f);
                for (int s=0;s<nSpecies;s++) { Yf[s] *= inv; Y[s] = (double)Yf[s]; }
            } else {
                double ysum = 0.0;
                for (int s=0;s<nSpecies;s++){
                    double y = (double)roY[s][ic]/(double)ro_temp;
                    if (y < 0.0f) y = 0.0f;
                    Y[s] = y; ysum += y;
                }
                double inv = 1.0/(ysum > 1.0e-30 ? ysum : 1.0e-30);
                for (int s=0;s<nSpecies;s++) Y[s] *= inv;
            }

            const double e_in = (double)intE;
            const double Tg   = ((double)T[ic] > DEPVAR_TMIN) ? (double)T[ic] : 300.0; // warm start

            // 非平衡凝縮 (一温度 二相 EOS): 総液相質量分率 g を集計。condensation==0 で g=0 (従来経路)。
            double g_liq = 0.0;
            const bool carrier = (condensation == 1 && condGasSpecies >= 0);
            if (condensation == 1 && rog != nullptr) {
                for (int s = 0; s < nCondSpecies; ++s) {
                    double gs = (double)rog[s][ic] / (double)ro_temp;
                    if (gs > 0.0f) g_liq += gs;
                }
                // realizability: carrier は g≤Y_凝縮種 (蒸気以上は凝縮しない)、pure は g≤0.99。
                if (carrier && roY != nullptr) {
                    double Yw = (double)roY[condGasSpecies][ic]/(double)ro_temp;
                    if (g_liq > Yw) g_liq = (Yw > 0.0f ? Yw : 0.0f);
                } else if (g_liq > 0.99f) g_liq = 0.99f;
                if (g_liq < 0.0f) g_liq = 0.0f;
            }

            const CondSpeciesProps cprops = condProps_make(condModel, condOpts);
            const double Rw = cprops.R;   // 凝縮種の比気体定数 (carrier: 蒸気分圧/EOS に使用)

            // 温度反転。g≈0 は従来 thermo_T_from_e で厳密縮約。
            double Tnew;
            bool hybrid = false; double hybrid_cp = 0.0, hybrid_h = 0.0;
            bool twophaseFail = false;   // 二相反転が収束しなかったセル: roe を上書きしない (保存量保護)
            if (condensation == 1 && condEquilibrium == 2 && rog != nullptr) {
                // EOS 拘束形平衡: g を状態量として (T,g) を同時反転し、rog[0] に射影する
                // (plans/accepted/condensation-equilibrium-eos.md)。輸送値 g_liq は初期値にだけ使う。
                double g_eq;
                if (carrier) {
                    const double Yw = (roY != nullptr) ? (double)roY[condGasSpecies][ic]/(double)ro_temp : 0.0;
                    g_eq = cond_equilibrium_Tg_carrier(sp, nSpecies, Y, e_in, (double)ro_temp, (Yw > 0.0 ? Yw : 0.0),
                                                       Rw, cprops, Tg, g_liq, DEPVAR_TMIN, DEPVAR_TMAX, &Tnew);
                } else {
                    g_eq = cond_equilibrium_Tg_pure_tp(sp, nSpecies, Y, e_in, (double)ro_temp, cprops,
                                                       Tg, g_liq, DEPVAR_TMIN, DEPVAR_TMAX, &Tnew);
                }
                g_liq = g_eq;
                rog[0][ic] = (flow_float)((double)ro_temp*g_eq);
            } else if (g_liq > 1.0e-12f) {
                // 一温度二相反転: condFloat=1 は float Newton + double 研磨 + 成功判定、それ以外は現行 double (成功判定のみ追加)。
                bool okc = true;
                if (useHybrid && condFloat != 0 && condTb.valid) {
                    Tnew = cond_T_from_e_twophase_hybrid(sp, spf, nSpecies, Y, Yf, condTb, e_in, g_liq, Rw, carrier ? 1 : 0, cprops,
                                                         Tg, DEPVAR_TMIN, DEPVAR_TMAX, &okc);
                } else {
                    // 旧 double Newton → 同じ残差条件 (1e-9|e|+0.05 J/kg) まで研磨 (全経路で成功条件を統一; codex result M2)。
                    Tnew = carrier ? cond_T_from_e_carrier(sp, nSpecies, Y, e_in, g_liq, Rw, cprops, Tg, DEPVAR_TMIN, DEPVAR_TMAX)
                                   : cond_T_from_e_onetemp(sp, nSpecies, Y, e_in, g_liq, Tg, DEPVAR_TMIN, DEPVAR_TMAX);
                    Tnew = cond_twophase_polish(sp, nSpecies, Y, Tnew, e_in, g_liq, Rw, carrier ? 1 : 0, cprops, DEPVAR_TMIN, DEPVAR_TMAX, &okc);
                }
                if (!okc) { twophaseFail = true; atomicAdd(&g_condTinvFail, 1u); }
            } else if (useHybrid) {
                // ハイブリッド: float Newton + double 1 段研磨。cp/h は研磨点 T_f の double 値から Taylor で組む
                // (double 評価 1 回で従来の反復数+1 回分を置換)。凝縮 off のときのみ (二相 EOS は従来経路)。
                double cpTf, hTf, Tf;
                Tnew = thermo_T_from_e_hybrid(sp, spf, nSpecies, Y, Yf, e_in, Tg, DEPVAR_TMIN, DEPVAR_TMAX, &cpTf, &hTf, &Tf, 12);
                hybrid_cp = cpTf; hybrid_h = hTf + cpTf*(Tnew - Tf); hybrid = true;
            } else {
                Tnew = thermo_T_from_e(sp, nSpecies, Y, e_in, Tg, DEPVAR_TMIN, DEPVAR_TMAX);
            }

            double Rmix;
            if (hybrid && sp[0].invMW > 0.0) { double s_ = 0.0; for (int s=0;s<nSpecies;s++) s_ += Y[s]*sp[s].invMW; Rmix = THERMO_RU * s_; }
            else Rmix = thermo_R_mix (sp, nSpecies, Y);
            double cpmix, hmix;
            if (hybrid) { cpmix = hybrid_cp; hmix = hybrid_h; }
            else thermo_cph_mix(sp, nSpecies, Y, Tnew, &cpmix, &hmix);  // cp,h を 1 スイープ (全蒸気混合)
            const double cvmix = cpmix - Rmix;
            // γ は出力 float なので、ハイブリッド経路は float の除算で十分 (double 除算を避ける)。
            const double gmix  = hybrid ? (double)((float)cpmix / (float)(cvmix > 1.0e-6 ? cvmix : 1.0e-6))
                                        : cpmix / (cvmix > 1.0e-6 ? cvmix : 1.0e-6);

            const double e_v   = hmix - Rmix*Tnew;
            const double Lcond = (g_liq > 1.0e-12) ? cond_latent(cprops, Tnew) : 0.0;
            double e_mix, Pnew, oneMg;
            if (carrier) {
                // carrier+condensible: e_mix=e_全蒸気+g(R_w T-L)、p=ρT(R_mix-g R_w)(凝縮で蒸気モル減)。
                e_mix = e_v + g_liq*(Rw*Tnew - Lcond);
                Pnew  = (double)ro_temp * Tnew * (Rmix - g_liq*Rw);
                oneMg = 1.0f;  // 気相質量は別途、p で表現済
            } else {
                // pure-condensible (気相=凝縮種): e_l=e_v+R_vT-L、p=(1-g)ρR T。
                e_mix = e_v + g_liq*Rmix*Tnew - g_liq*Lcond;
                oneMg = 1.0f - g_liq;
                Pnew  = (double)ro_temp * oneMg * Rmix * Tnew;
            }
            if (Pnew < (double)pMin) Pnew = (double)pMin;

            T[ic]         = (flow_float)Tnew;
            P[ic]         = (flow_float)Pnew;
            ro[ic]        = ro_temp;
            // roe を (floor 済 ro, 反転 T, 混合内部エネルギー) と整合させて再構成
            if (!twophaseFail) roe[ic] = (flow_float)((double)ro_temp * (e_mix + (double)ek));
            // 総エンタルピー Ht = e_mix + p/ρ + ek (g=0 → hmix+ek と一致)。反転失敗セルは保存量 roe から (整合を保つ)。
            Ht[ic]        = twophaseFail ? (flow_float)((double)roe[ic]/(double)ro_temp + (double)Pnew/(double)ro_temp)
                                         : (flow_float)(e_mix + (double)Pnew/(double)ro_temp + (double)ek);
            // 音速と γ: 既定 (condSonicModel 0 / g=0) は全蒸気気相 √(γ_mix R_mix T)。condSonicModel 1 かつ g>0 では
            // 一温度二相 EOS と整合する固定 g,Y の frozen 音速 c²=γ_2φ R_eff T (cond_twophase_sonic)。γ_2φ は block-DPLUR の
            // κ=γ−1 (固定 g,Y の frozen 近似) と TP 出口 BC が読む。g<1e-12 は式順序も従来と同一 (dry セル bit 同一)。
            double sonic2 = gmix * Rmix * Tnew, gam_out = gmix;
            if (condSonicModel == 1 && g_liq > 1.0e-12f) {
                // 二相音速の dL/dT は従来の double 両側差分のまま (表の片側微分は液相クランプ点 373.15 K で ±0.1 K の差分と 4 % 違い、
                // γ_2φ・音速 (陰解法の係数) が変わる: codex result-4 M2 で表微分案を撤回)。表微分は Newton の傾き (結果に効かない) にだけ使う。
                const double dL   = (cond_latent(cprops, Tnew + 0.1) - cond_latent(cprops, Tnew - 0.1)) / 0.2;
                const double Reff = carrier ? (Rmix - g_liq*Rw) : ((1.0 - g_liq)*Rmix);
                double g2, c2;
                if (cond_twophase_sonic(cpmix, Reff, g_liq, dL, Tnew, &g2, &c2)) { sonic2 = c2; gam_out = g2; }
            }
            sonic[ic]     = hybrid ? sqrtf((flow_float)sonic2) : (flow_float)sqrt(sonic2);
            gam_array[ic] = (flow_float)gam_out;
            cp_array[ic]  = (flow_float)cpmix;
            Rmix_array[ic]= (flow_float)Rmix;
        } else {
            // ---- calorically perfect gas ----
            // 非平衡凝縮 (一温度 二相 EOS, CPG): 総液相質量分率 g を集計。condensation==0 で g=0 (従来経路)。
            // CPG carrier (空気の N2 選択凝縮, condOpts.Yw>0): g<=Y_w、p=ρT(R_air−gR_w)、e=(c_v+gR_w)T−gL (plans/accepted/condensation-air.md §4.1)。
            const bool carrierCpg = (condensation == 1 && condOpts.Yw > 0.0);
            double g_liq = 0.0;
            if (condensation == 1 && rog != nullptr) {
                for (int s = 0; s < nCondSpecies; ++s) {
                    double gs = (double)rog[s][ic] / (double)ro_temp;
                    if (gs > 0.0f) g_liq += gs;
                }
                const double gcap = carrierCpg ? condOpts.Yw : 0.99;
                if (g_liq > gcap) g_liq = gcap;   // realizability
                if (g_liq < 0.0)  g_liq = 0.0;
            }

            const double cv = (double)cp/(double)gamma;            // cv = cp/γ
            const double Rgas = ((double)gamma - 1.0)*cv;          // R = cp - cv = (γ-1)cv (CPG carrier では空気の R)

            // EOS 拘束形平衡 (pure CPG): (T,g) 同時反転 → rog[0] 射影。g_eq>0 なら下の二相分岐へ流す。
            double Tn_eq = 0.0;
            const bool eq2 = (condensation == 1 && condEquilibrium == 2 && rog != nullptr);
            if (eq2) {
                const CondSpeciesProps cpropsEq = condProps_make(condModel, condOpts);
                const double Tguess = ((double)T[ic] > 1.0) ? (double)T[ic] : (double)max(intE/(cp/gamma), tMin);
                g_liq = cond_equilibrium_Tg_pure_cpg((double)intE, (double)ro_temp, cv, Rgas, cpropsEq, Tguess, g_liq, &Tn_eq);
                rog[0][ic] = (flow_float)((double)ro_temp*g_liq);
            }

            if (g_liq > 1.0e-12) {
                // 二相: e = (cv + g R_w) T - g L(T) = intE を括弧付き Newton で反転、p=ρ T R_eff (pure: R_w=R, R_eff=(1-g)R; carrier: R_eff=R_air−gR_w)。
                const double e_in = (double)intE;
                const CondSpeciesProps cpropsCpg = condProps_make(condModel, condOpts);
                const double Rw = carrierCpg ? cpropsCpg.R : Rgas;
                // 初期推定は前ステップの T (warm start; 無ければ単相値)
                const double Tguess = ((double)T[ic] > 1.0) ? (double)T[ic] : (double)max(intE/(cp/gamma), tMin);
                bool ok = true;
                const double Tn = eq2 ? Tn_eq : cond_T_from_e_cpg(e_in, g_liq, cv, Rw, Tguess, cpropsCpg, &ok);
                if (!ok) {
                    // 反転が収束しなかったセル: T,P,sonic,Ht は前ステップ値のまま残し、roe を反転結果で上書きしない (密度床・速度は上で更新済み)。診断カウンタに数える
                    // (codex 2026-09-12 M2 / 2026-09-13 M1: 失敗した温度で流束・核生成を評価しない)。密度床だけ反映。
                    atomicAdd(&g_condTinvFail, 1u);
                    ro[ic] = ro_temp;
                    Rmix_array[ic] = (flow_float)Rgas;
                } else {
                const double L = cond_latent(cpropsCpg, Tn);
                const double e_mix = (cv + g_liq*Rw)*Tn - g_liq*L;   // = e_in
                const double Reff = carrierCpg ? (Rgas - g_liq*Rw) : ((1.0 - g_liq)*Rgas);
                double Pn = (double)ro_temp*Reff*Tn;
                if (Pn < (double)pMin) Pn = (double)pMin;

                T[ic]   = (flow_float)Tn;
                P[ic]   = (flow_float)Pn;
                ro[ic]  = ro_temp;
                roe[ic] = (flow_float)((double)ro_temp*(e_mix + (double)ek));
                Ht[ic]  = (flow_float)(e_mix + Pn/(double)ro_temp + (double)ek);
                sonic[ic] = (flow_float)sqrt((double)gamma*Rgas*Tn); // 気相 frozen 音速 (loose coupling; CPG は旧式のまま)
                Rmix_array[ic] = (flow_float)Rgas;
                }
            } else {
                // 単相 CPG (従来経路, フロア未指定ならビット不変)
                T_temp = max(intE/(cp/gamma), tMin);
                P_temp = max((gamma-1.0f)*(roe[ic]-ro_temp*ek), pMin);

                T[ic] = T_temp;
                P[ic] = P_temp;

                ro[ic] = ro_temp;
                roe[ic] = P_temp/(gamma-1.0f) + ro_temp*ek;

                Ht[ic] = roe[ic]/ro_temp + P_temp/ro_temp;

                sonic[ic] = sqrt(gamma*P_temp/ro_temp);
                // CPG の混合比気体定数 R = cp - cv = (γ-1)cp/γ (定数。SLAU 単成分経路は未使用だが整合のため埋める)
                Rmix_array[ic] = (gamma-1.0f)*cp/gamma;
            }
        }

        k[ic] = max(roK[ic]/ro_temp, static_cast<flow_float>(0.0));
        omega[ic] = max(roOmega[ic]/ro_temp, static_cast<flow_float>(0.0));
    }
}


void dependentVariables_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    dependentVariables_d<<<cuda_cfg.dimGrid_cell , cuda_cfg.dimBlock>>> (
        // gas properties
        cfg.thermalMethod ,
        cfg.gamma , cfg.cp ,

        // EOS 正値化フロア (config 化)
        cfg.pMin , cfg.roMin , cfg.tMin ,

        // thermally-perfect 用化学種データ。多成分 (M2, nSpecies>=2) では device roY 配列を渡し、
        // 単成分のときは nullptr (混合則は Y={1} に縮退)。
        thermo_species_device_ptr() , cfg.nSpecies , species_roY_device_ptr() ,
        thermo_species_device_ptr_f() , cfg.thermoFloat ,

        // 非平衡凝縮 (二相 EOS)。condensation==0 で rog=nullptr/g=0 → 従来経路ビット不変。
        cfg.condensation , var.nCondSpeciesRegistered , cond_rog_device_ptr() ,
        cfg.condGasSpecies , cfg.condModel ,
        cfg.condEquilibrium ,
        cfg.condSonicModel ,
        cond_prop_opts(cfg) ,
        cfg.condFloat , cond_tables_device() ,

        // mesh structure
        msh.nCells_all , msh.nCells ,

        // basic variables
        var.c_d["ro"]  , var.c_d["roUx"], var.c_d["roUy"] , var.c_d["roUz"], var.c_d["roe"] ,
        var.c_d["roK"] , var.c_d["roOmega"],
        var.c_d["P"]   , var.c_d["Ht"]  , var.c_d["sonic"], var.c_d["k"], var.c_d["omega"], var.c_d["T"],
        var.c_d["Ux"]  , var.c_d["Uy"]  , var.c_d["Uz"] ,

        var.c_d["gamma"] , var.c_d["cp"] , var.c_d["Rmix"]
    ) ;
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
    // 二相の温度反転失敗セル数 (診断; CPG と TP 一温度二相の両方)。0 でなければ警告 (roe を保持したセルがある; plan condensation-float-speedup §5.1 #9)。
    if (cfg.condensation == 1) {
        unsigned int nfail = 0u;
        gpuErrchk( cudaMemcpyFromSymbol(&nfail, g_condTinvFail, sizeof(unsigned int)) );
        if (nfail > 0u) {
            std::cerr << "[condensation] WARNING: two-phase temperature inversion failed in " << nfail << " cells (roe kept; " << (cfg.thermalMethod == 0 ? "CPG: primitives kept from previous step" : "TP: T/P from the unconverged inversion") << ")\n";
            const unsigned int zero = 0u; gpuErrchk( cudaMemcpyToSymbol(g_condTinvFail, &zero, sizeof(unsigned int)) );
        }
    }
    // EOS 拘束形平衡 (condEquilibrium==2): kernel が rog[0] を g_eq に射影したので、SLAU 面温度・潜熱補正・出力が
    // 読む原始量 g_<s>=rog/ρ を同期する (realizability クランプも通る)。他モードでは呼ばない (従来経路不変)。
    if (cfg.condensation == 1 && cfg.condEquilibrium == 2 && var.nCondSpeciesRegistered >= 1) {
        condensationPrimitive_d_wrapper(cfg , cuda_cfg , msh , var);
    }
}
