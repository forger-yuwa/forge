// test_cond_float.cpp — 凝縮経路 float 化の host 単体検証 (plans/active/condensation-float-speedup.md §4.3)。
//   build: nvcc -x cu --expt-relaxed-constexpr -I. -o test_cond_float tests/unit/test_cond_float.cpp
// (1) 物性表 vs double 関数: 0.01 K 刻み全域 (接続点・床・臨界直下を含む)。
//     許容 (plan §4.3 v2): ln p_sat 絶対 ≤ 2e-6 + 1.5e-7·|ln p| (|ln p|~20–50 では float の表現限界 6e-8·|ln p| が支配),
//     L/ρ_l/k_gas/μ_gas 相対 ≤2e-6 (H2O の L は凝縮が起き得る T ≤ 400 K で厳密、それ以上は L 床 1.5 MJ/kg の折れ点 (~968 K) を含むので ≤1e-4),
//     σ は相対 ≤2e-6 または絶対 ≤1e-8 N/m (σ<1e-3 N/m の T_c 直下は絶対 ≤1e-5; べき 1.247/1.256 の微分特異点)。
//     微分: d ln p/dT 相対 ≤1e-4、dL/dT は |ΔL'| ≤ 2e-4|L'| + 0.1 J/kg/K (極値近傍の相対は意味が無い; Newton の傾き専用、音速の dL/dT は double のまま)。接続点 ±0.02 K は除外。
#include <cstdio>
#include <cmath>
#include <vector>
#include "cuda_forge/condensationTables_d.cuh"
#include "cuda_forge/condensationSourceF_d.cuh"
#include <cstdlib>

static int g_fail = 0;
static void check(bool ok, const char* what, double val, double tol)
{
    printf("  [%s] %-70s worst=%.3e tol=%.1e\n", ok ? "PASS" : "FAIL", what, val, tol);
    if (!ok) ++g_fail;
}

struct Worst { double v = 0.0, T = 0.0; void upd(double e, double T_) { if (e > v) { v = e; T = T_; } } };

static void test_tables(const char* name, const CondSpeciesProps& s, double Tlo, double Thi)
{
    CondTablesHost ht; cond_tables_build_host(s, ht);
    const CondTablesF tb = cond_tables_view_host(ht);
    printf("-- %s: T0=%.2f h=%.3f n=%d range [%.2f, %.2f] test [%.1f, %.1f]\n", name, ht.T0, ht.h, ht.n, tb.Tmin, tb.Tmax, Tlo, Thi);
    Worst wlp, wL, wLhi, ws, wr, wk, wm, wdlp, wdL;
    for (double T = Tlo; T <= Thi + 1e-9; T += 0.01) {
        const float Tf = (float)T;
        // 表の範囲外は端クランプ: double 側も同じ温度で評価 (範囲内の誤差だけを判定)
        const double Tc = (T < tb.Tmin) ? tb.Tmin : ((T > tb.Tmax) ? tb.Tmax : T);
        float dlp, dL;
        const double lp = cond_tab_lnpsat_f(tb, Tf, &dlp), lpd = log(cond_psat(s, Tc));
        wlp.upd(fabs(lp - lpd)/(1.0 + 0.15*fabs(lpd)), T);    // 許容 2e-6·(1 + 0.15|ln p|) = 2e-6 + 3e-7|ln p| (~5 ulp of |ln p|)
        const double L = cond_tab_latent_f(tb, Tf, &dL), Ld = cond_latent(s, Tc);
        if (s.model != COND_MODEL_H2O || T <= 400.0) wL.upd(fabs(L - Ld)/Ld, T); else wLhi.upd(fabs(L - Ld)/Ld, T);
        const double sg = cond_tab_sigma_f(tb, Tf), sgd = cond_sigma(s, Tc);
        { const double e = fabs(sg - sgd); ws.upd((e <= 1e-8 || (sgd < 1e-3 && e <= 1e-5)) ? 0.0 : e/sgd, T); }   // T_c 直下 (σ<1e-3) は絶対 1e-5
        const double rl = cond_tab_rhol_f(tb, Tf), rld = cond_rho_cond(s, Tc);       wr.upd(fabs(rl - rld)/rld, T);
        const double kg = cond_tab_kgas_f(tb, Tf), kgd = cond_kgas(s, Tc);           wk.upd(fabs(kg - kgd)/kgd, T);
        const double mu = cond_tab_mugas_f(tb, Tf), mud = n2_mu_gas(Tc);             wm.upd(fabs(mu - mud)/mud, T);
        // 微分 (中心差分の double と比較; 接続点 ±0.02 K は除外)
        const bool nearJ = (s.model == COND_MODEL_H2O) ? (fabs(T - 273.15) < 0.02 || fabs(T - 373.15) < 0.02 || fabs(T - 1000.0) < 0.02)
                                                       : (fabs(T - 45.0) < 0.02 || fabs(T - 50.0) < 0.02 || fabs(T - 70.0) < 0.02);
        if (T > tb.Tmin + 0.02 && T < tb.Tmax - 0.02 && !nearJ) {
            const double d = 1e-3;
            const double dlpd = (log(cond_psat(s, T + d)) - log(cond_psat(s, T - d)))/(2*d);
            const double dLd  = (cond_latent(s, T + d) - cond_latent(s, T - d))/(2*d);
            wdlp.upd(fabs(dlp - dlpd)/(fabs(dlpd) + 1e-4*fabs(lpd) + 1e-6), T);
            if (s.model != COND_MODEL_H2O || T <= 400.0) wdL.upd(0.5*fabs(dL - dLd)/(fabs(dLd) + 500.0), T);   // |ΔL'| ≤ 2e-4|L'| + 0.1 J/kg/K (Newton の傾きにしか使わない; H2O 液相 NASA-9 の高次項で 3 次内挿の微分誤差が 1.5e-4)
        }
    }
    char b[128];
    snprintf(b, sizeof b, "%s ln p_sat abs/(1+0.075|ln p|) (worst at %.2f K)", name, wlp.T); check(wlp.v <= 2e-6, b, wlp.v, 2e-6);
    snprintf(b, sizeof b, "%s L rel (worst at %.2f K)", name, wL.T);         check(wL.v <= 2e-6, b, wL.v, 2e-6);
    if (s.model == COND_MODEL_H2O) { snprintf(b, sizeof b, "%s L rel T>400K (L floor kink; worst at %.2f K)", name, wLhi.T); check(wLhi.v <= 1e-4, b, wLhi.v, 1e-4); }
    snprintf(b, sizeof b, "%s sigma rel|abs<=1e-8 (worst at %.2f K)", name, ws.T); check(ws.v <= 2e-6, b, ws.v, 2e-6);
    snprintf(b, sizeof b, "%s rho_l rel (worst at %.2f K)", name, wr.T);     check(wr.v <= 2e-6, b, wr.v, 2e-6);
    snprintf(b, sizeof b, "%s k_gas rel (worst at %.2f K)", name, wk.T);     check(wk.v <= 2e-6, b, wk.v, 2e-6);
    snprintf(b, sizeof b, "%s mu_gas rel (worst at %.2f K)", name, wm.T);    check(wm.v <= 2e-6, b, wm.v, 2e-6);
    snprintf(b, sizeof b, "%s d(ln p_sat)/dT rel (worst at %.2f K)", name, wdlp.T); check(wdlp.v <= 1e-4, b, wdlp.v, 1e-4);
    snprintf(b, sizeof b, "%s dL/dT rel (worst at %.2f K)", name, wdL.T);    check(wdL.v <= 1e-4, b, wdL.v, 1e-4);
}


// (2) 連鎖 (核生成 J・r*・成長 dr/dt・ソースベクトル・蒸発 λ・T_sat) を float 実体 vs double 関数で比較。
//     許容 (plan §4.3 v2): |Δ ln J| ≤ 2(ΔG*/k_BT)·2e-6 + 1e-4 (double の J が J_max 超のセルは float 上限で除外・計数),
//     r* 相対 ≤1e-5 + 1e-6/ln S, dr/dt ≤ 1e-4|dr/dt| + 1e-5×駆動流束尺度 (dr/dt>0 のみ), ソースベクトル相対 ≤1e-3 (物理的ゼロ [J<1e-20] は除外),
//     蒸発 λ ≤1e-4, T_sat ≤2e-3 K。N2 は表上端 125.5 K まで (125.6–125.69 K は double 側のクランプ点との差 = 文書化した近似)。
struct ChainWorst { double v = 0; double T = 0, S = 0; long n = 0, capped = 0; char info[200] = ""; void upd(double e, double T_, double S_, const char* i = nullptr) { ++n; if (e > v) { v = e; T = T_; S = S_; if (i) snprintf(info, sizeof info, "%s", i); } } };

static void test_chain(const char* name, const CondSpeciesProps& s, int carrier, double Tlo, double Thi, double dT)
{
    CondTablesHost ht; cond_tables_build_host(s, ht);
    const CondTablesF tb = cond_tables_view_host(ht);
    const CondSpeciesPropsF sf = condProps_to_f(s);
    const double Svals[] = {1.001, 1.01, 1.1, 1.5, 2.0, 5.0, 10.0, 30.0, 100.0, 300.0, 1000.0};
    const double Q0vals[] = {0.0, 1.0e14, 1.0e17, 1.0e19};
    const double rvals[] = {1.0e-9, 5.0e-9, 3.0e-8, 3.0e-7, 1.0e-6};
    ChainWorst wJ, wr, wg[2], wsv, wjac;
    long nJpos = 0;
    printf("-- chain %s (carrier=%d) T %.0f..%.0f\n", name, carrier, Tlo, Thi);
    for (int mode = 0; mode <= 3; ++mode) for (double T = Tlo; T <= Thi + 1e-9; T += dT) {
        const double psat = cond_psat(s, T);
        for (double S : Svals) {
            const double pv = S*psat;
            const double rho_v = carrier ? pv/(s.R*T) : 1.0;   // carrier: 蒸気密度=分圧/(R T), pure: 適当な気相密度
            CondNucCarrier car; car.a_v = carrier ? 0.02/s.M : 1.0; car.carrierSum = carrier ? 30.0 : 0.0; car.cvv_tilde = s.cv*s.M/COND_RU;
            CondNucCarrierF carf; carf.a_v = (float)car.a_v; carf.carrierSum = (float)car.carrierSum; carf.cvv_tilde = (float)car.cvv_tilde;
            const double gam = s.cp/s.cv;
            double Jd, rd; cond_nucleation(s, T, pv, rho_v, &Jd, &rd, mode, gam, &car);
            float  Jf, rf;  cond_nucleation_f(sf, tb, (float)T, (float)pv, (float)rho_v, &Jf, &rf, mode, (float)gam, &carf);
            if (rd > 0.0) wr.upd(fabs(rf - rd)/rd/(1.0 + 0.2/log(S)), T, S);   // 許容 1e-5·(1 + 0.2/ln S) = 1e-5 + 2e-6/ln S (ln S の float 絶対誤差 ~1e-6〜2e-6)
            if (Jd > 0.0 && Jf > 0.0) {
                ++nJpos;
                if (Jd > 1.0e35) { ++wJ.capped; }
                else {
                    const double sigma = cond_sigma(s, T), rho_l = cond_rho_cond(s, T);
                    const double dGkT = (4.0/3.0)*COND_PI*rd*rd*sigma/(COND_KB*T);
                    const double tol = 2.0*dGkT*(3.0e-6/log(S)) + 1.0e-4;   // δ(ΔG/kT) = 2(ΔG/kT)·δlnS/lnS, δlnS ≈ 3e-6 (表 2e-6 + ln p_v の float 丸め); S→1 で CNT 自体が発散的に敏感
                    wJ.upd(fabs(log(Jf) - log(Jd))/tol, T, S);   // 正規化: ≤1 で合格
                }
            } else if ((Jd > 0.0) != (Jf > 0.0)) {
                // 片方だけ 0: double が 1e-300 級 (float の -80 切り捨て) のときだけ許す
                if (Jd > 1.0e-30 || Jf > 1.0e-30) { printf("   J zero mismatch T=%.2f S=%.3f mode=%d Jd=%.3e Jf=%.3e\n", T, S, mode, Jd, (double)Jf); ++g_fail; }
            }
            if (mode == 0) for (int gm = 0; gm <= 1; ++gm) for (double Q0 : Q0vals) for (double rb : rvals) {
                const double Q1 = Q0*rb, Q2 = Q0*rb*rb;
                double a0,a1,a2,ag; cond_source_vector(s, T, pv, rho_v, Q0, Q1, Q2, &a0,&a1,&a2,&ag, 1, gm, gam, carrier ? 1.0e5 : -1.0, 3.18, 0, &car);
                float  b0,b1,b2,bg; cond_source_vector_f(sf, tb, (float)T, (float)pv, (float)rho_v, (float)Q0, (float)Q1, (float)Q2, &b0,&b1,&b2,&bg, 1, gm, (float)gam, carrier ? 1.0e5f : -1.0f, 3.18f, &carf);
                if (Jd > 1.0e35) continue;   // 上限セルは float 実体の設計変更 (対数上限) で不一致
                const double A[4] = {a0,a1,a2,ag}; const double B[4] = {b0,b1,b2,bg};
                // 許容: (a) 相対 1e-3 + ln J の許容 (J 駆動項は CNT の S→1 での敏感さを継承), または
                //       (b) 絶対 ≤ 1e-4 × 成長流束尺度 (Q_n·(p_v/ρ_l)/√(2πRT)): r̄≈r* の (1−r*/r̄) や S≈1 の (p_v−p_d) の相殺で相対は意味を失い、
                //           物理的には dr/dt≈0 (kernel も律速・切り上げで扱う)。J<1e-5 /m³/s 級は「物理的ゼロ」(1e-6 m³·1e-6 s で 1e-17 個) として除外。
                const double sigma_ = cond_sigma(s, T), rho_l_ = cond_rho_cond(s, T);
                const double dGkT_ = (rd > 0.0) ? (4.0/3.0)*COND_PI*rd*rd*sigma_/(COND_KB*T) : 0.0;
                const double tolJ_ = 2.0*dGkT_*(3.0e-6/log(S)) + 1.0e-4;
                const double dsc = (pv/rho_l_)/sqrt(2.0*COND_PI*s.R*T);
                const double scale_k[4] = {0.0, Q0*dsc, 2.0*Q1*dsc, 4.0*COND_PI*rho_l_*Q2*dsc};
                const double floor_k[4] = {1.0e-5, 1.0e-14, 1.0e-23, 1.0e-29};
                for (int k = 0; k < 4; ++k) {
                    if (fabs(A[k]) < floor_k[k] && fabs(B[k]) < floor_k[k]) continue;
                    const double den = fabs(A[k]) > fabs(B[k]) ? fabs(A[k]) : fabs(B[k]);
                    const double err = fabs(A[k] - B[k]);
                    if (err <= 1.0e-4*scale_k[k]) continue;
                    char inf[200]; snprintf(inf, sizeof inf, "k=%d gm=%d Q0=%.1e rb=%.1e A=%.4e B=%.4e Jd=%.3e rd=%.3e rf=%.3e tolJ=%.2e", k, gm, Q0, rb, A[k], B[k], Jd, rd, (double)rf, tolJ_);
                    wsv.upd(err/den/(1.0e-3 + tolJ_), T, S, inf);   // 正規化: ≤1 で合格
                }
                if (Q0 > 0.0 && rd > 0.0 && rb > rd) {
                    const double dd = cond_growth(s, T, pv, rb, rd, gm, carrier ? 1.0e5 : -1.0, 3.18, 0);
                    const float  df = cond_growth_f(sf, tb, (float)T, (float)pv, (float)rb, (float)rd, gm, carrier ? 1.0e5f : -1.0f, 3.18f);
                    if (dd > 0.0) {   // kernel は dr/dt<0 を 0 に切る
                        // 許容: 1e-4|dr/dt| + 1e-5 × 駆動流束尺度 (p_v/ρ_l)/√(2πRT) (S→1 の (p_v−p_d) 相殺で相対は意味を失う)
                        const double scale = (pv/cond_rho_cond(s, T))/sqrt(2.0*COND_PI*s.R*T);
                        wg[gm].upd(fabs(df - dd)/((1.0e-4 + 2.0e-6/log(S))*fabs(dd) + 1.0e-5*scale), T, S);   // 駆動力 ln S の float 絶対誤差 2e-6 を相対に含める
                    }
                }
            }
        }
        // 蒸発 (S<1) と T_sat
        for (double S : {0.2, 0.5, 0.9, 0.999}) for (double Q0 : {1.0e15, 1.0e18}) for (double g : {1.0e-6, 1.0e-4, 1.0e-2}) {
            const double pv = S*psat; const double rod = 1.0;
            double a0,a1,a2,ag,r30,dr; const double lamd = cond_evap_source(s, T, pv, rod, g, Q0, Q0*3e-8, Q0*9e-16, 1.0e-6, 1.0e-9, 0.5, 5.0e-3, 1.0, 742.0, 0, carrier ? 1.0e5 : -1.0, 3.18, 0, &a0,&a1,&a2,&ag,&r30,&dr);
            float  b0,b1,b2,bg,r30f,drf; const float lamf = cond_evap_source_f(sf, tb, (float)T, (float)pv, (float)rod, (float)g, (float)Q0, (float)(Q0*3e-8), (float)(Q0*9e-16), 1.0e-6f, 1.0e-9f, 0.5f, 5.0e-3f, 1.0f, 742.0f, 0, carrier ? 1.0e5f : -1.0f, 3.18f, 0, &b0,&b1,&b2,&bg,&r30f,&drf);
            wjac.upd(fabs(lamf - lamd), T, S);
        }
        {
            const double pv = 0.3*psat;
            const double Ts = cond_Tsat(s, pv, T*1.1); const float Tsf = cond_Tsat_f(tb, (float)pv, (float)(T*1.1));
            if (Ts > 0.0 && fabs(Tsf - Ts) > 2.0e-3) { printf("   Tsat mismatch T=%.2f: %.5f vs %.5f\n", T, Ts, (double)Tsf); ++g_fail; }
        }
    }
    char b[160];
    snprintf(b, sizeof b, "%s ln J (normalized by 2(dG/kT)(3e-6/lnS)+1e-4; %ld J>0, %ld capped) worst T=%.1f S=%.3g", name, nJpos, wJ.capped, wJ.T, wJ.S); check(wJ.v <= 1.0, b, wJ.v, 1.0);
    snprintf(b, sizeof b, "%s r* rel/(1+0.2/lnS) (worst T=%.1f S=%.3g)", name, wr.T, wr.S);           check(wr.v <= 1e-5, b, wr.v, 1e-5);
    snprintf(b, sizeof b, "%s dr/dt model0 /((1e-4+2e-6/lnS)|dd|+1e-5 scale) (worst T=%.1f S=%.3g)", name, wg[0].T, wg[0].S); check(wg[0].v <= 1.0, b, wg[0].v, 1.0);
    snprintf(b, sizeof b, "%s dr/dt model1 /((1e-4+2e-6/lnS)|dd|+1e-5 scale) (worst T=%.1f S=%.3g)", name, wg[1].T, wg[1].S); check(wg[1].v <= 1.0, b, wg[1].v, 1.0);
    snprintf(b, sizeof b, "%s source vector rel/(1e-3+tolJ) | abs 1e-4·scale (worst T=%.1f S=%.3g)", name, wsv.T, wsv.S); check(wsv.v <= 1.0, b, wsv.v, 1.0); if (wsv.v > 1.0) printf("      %s\n", wsv.info);
    snprintf(b, sizeof b, "%s evaporation lambda abs (worst T=%.1f S=%.3g)", name, wjac.T, wjac.S); check(wjac.v <= 5e-4, b, wjac.v, 5e-4);   // S→1 (0.999) の (p_v−p_d) 相殺で相対 2e-3 → λ 絶対 1.4e-4; 半径縮小比の誤差 <5e-4/step
}

int main()
{
    printf("== (1) property tables vs double ==\n");
    CondPropOpts o; o.latentLowT = 1; o.psatLowT = 1; o.liquidCp = 2000.0; o.gasKgasModel = 0; o.sigmaScale = 1.0; o.Yw = 0.0;
    test_tables("N2 (default)", condProps_make(COND_MODEL_N2, o), 20.0, 125.6);
    CondPropOpts o2 = o; o2.latentLowT = 0; o2.psatLowT = 0; o2.gasKgasModel = 1; o2.sigmaScale = 1.03;
    test_tables("N2 (old lowT, air kgas, sigma x1.03)", condProps_make(COND_MODEL_N2, o2), 20.0, 125.6);
    test_tables("H2O", condProps_make(COND_MODEL_H2O, o), 120.15, 1200.15);
    printf("== (2) source chain float vs double ==\n");
    test_chain("N2 pure", condProps_make(COND_MODEL_N2, o), 0, 39.0, 125.5, 1.0);
    test_chain("N2 (old lowT) pure", condProps_make(COND_MODEL_N2, o2), 0, 39.0, 125.5, 1.0);
    test_chain("H2O carrier", condProps_make(COND_MODEL_H2O, o), 1, 150.0, 400.0, 2.5);
    printf("%s (%d failures)\n", g_fail ? "FAILED" : "ALL PASS", g_fail);
    return g_fail ? 1 : 0;
}
