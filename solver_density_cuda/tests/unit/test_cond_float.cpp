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
#include "cuda_forge/condensationEOS_d.cuh"
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
        const double Tc = T;   // 参照は実温度 (codex result M1): 表範囲外は float 側が double 関数へ退避する (下の *_tab_or_d) ので一致するはず
        float dlp = 0.0f, dL = 0.0f;
        // ln p_sat は dry 診断 (S, T_sat) 用に表範囲 [Tmin, Tmax] で表を使う (外は端クランプ = 診断のみで判定対象外; 湿潤セルは double 実体へ退避)
        if (T >= tb.Tmin && T <= tb.Tmax) {
            const double lp = cond_tab_lnpsat_f(tb, Tf, &dlp), lpd = log(cond_psat(s, Tc));
            wlp.upd(fabs(lp - lpd)/(1.0 + 0.15*fabs(lpd)), T);    // 許容 2e-6·(1 + 0.15|ln p|) = 2e-6 + 3e-7|ln p| (~5 ulp of |ln p|)
        }
        // 湿潤経路の物性は表範囲 [Tmin, TwetMax] の外で double 関数へ退避: 退避値 = double の float 変換 (丸め 6e-8) を確認
        const double L = cond_tab_wet_ok(tb, Tf) ? cond_tab_latent_f(tb, Tf, &dL) : cond_latent_tab_or_d(tb, s, Tf), Ld = cond_latent(s, Tc);
        if (!cond_tab_wet_ok(tb, Tf)) { const double Lfb = cond_latent(s, (double)Tf); if (fabs(L - Lfb) > 1.5e-7*Lfb) { printf("   fallback L mismatch at %.2f K\n", T); ++g_fail; } continue; }   // 退避は同じ float 入力 Tf の double 関数値 (T_c 直下は L の傾きが急で T の float 丸めが効く)
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
            wdlp.upd(fabs(dlp - dlpd)/(fabs(dlpd) + 1e-4*fabs(log(cond_psat(s, T))) + 1e-6), T);
            // |ΔL'| ≤ 2e-4|L'| + 0.1 J/kg/K (Newton の傾きにしか使わない; H2O 液相 NASA-9 の高次項で 3 次内挿の微分誤差が 1.5e-4)。
            // H2O は 646.15 K の Watson 凍結点 (折れ点) をまたぐ区間だけ除く (2026-09-14)。
            if (s.model != COND_MODEL_H2O || T <= 645.8 || T >= 646.5) wdL.upd(0.5*fabs(dL - dLd)/(fabs(dLd) + 500.0), T);
        }
    }
    char b[128];
    snprintf(b, sizeof b, "%s ln p_sat abs/(1+0.075|ln p|) (worst at %.2f K)", name, wlp.T); check(wlp.v <= 2e-6, b, wlp.v, 2e-6);
    snprintf(b, sizeof b, "%s L rel (worst at %.2f K)", name, wL.T);         check(wL.v <= 2e-6, b, wL.v, 2e-6);
    if (s.model == COND_MODEL_H2O) { snprintf(b, sizeof b, "%s L rel T>400K (Watson branch + 646.15 K freeze; worst at %.2f K)", name, wLhi.T); check(wLhi.v <= 1e-4, b, wLhi.v, 1e-4); }
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

// (3) 一温度二相反転のハイブリッド (float Newton → double 研磨 → 残差判定) vs 厳密 double 参照。N2/H2O carrier (datum 298.15 K)。
//     許容: |ΔT| ≤ 1e-8·T (研磨後), float 格納 roe の 10 往復ドリフト ≤ 従来 double 反転のドリフト × 1.5 + 1e-7·T, 失敗 (ok=false) 0 件。
static SpeciesThermo mk(double MW,double sig,double eps,const double lo[9],const double hi[9]){
    SpeciesThermo s; s.MW=MW; s.sigma_LJ=sig; s.eps_kB=eps; s.Tlo=200.0; s.Tmid=1000.0; s.Thi=6000.0; s.h_datum=0.0; s.invMW=1.0/MW;
    for(int i=0;i<9;i++){ s.low[i]=lo[i]; s.high[i]=hi[i]; } return s; }
static SpeciesThermoF toF(const SpeciesThermo& s){
    SpeciesThermoF f; f.MW=(float)s.MW; f.invMW=(float)(1.0/s.MW); f.R=(float)(THERMO_RU/s.MW);
    f.sigma_LJ=(float)s.sigma_LJ; f.eps_kB=(float)s.eps_kB; f.Tlo=(float)s.Tlo; f.Tmid=(float)s.Tmid; f.Thi=(float)s.Thi;
    for(int k=0;k<9;k++){ f.low[k]=(float)s.low[k]; f.high[k]=(float)s.high[k]; } return f; }

static void test_inversion()
{
    const double N2lo[9]={2.210371497e+04,-3.818461820e+02,6.082738360e+00,-8.530914410e-03,1.384646189e-05,-9.625793620e-09,2.519705809e-12,7.108460860e+02,-1.076003744e+01};
    const double N2hi[9]={5.877124060e+05,-2.239249073e+03,6.066949220e+00,-6.139685500e-04,1.491806679e-07,-1.923105485e-11,1.061954386e-15,1.283210415e+04,-1.586640027e+01};
    const double H2Olo[9]={-3.947960830e+04,5.755731020e+02,9.317826530e-01,7.222712860e-03,-7.342557370e-06,4.955043490e-09,-1.336933246e-12,-3.303974310e+04,1.724205775e+01};
    const double H2Ohi[9]={1.034972096e+06,-2.412698562e+03,4.646110780e+00,2.291998307e-03,-6.836830480e-07,9.426468930e-11,-4.822380530e-15,-1.384286509e+04,-7.978148510e+00};
    std::vector<SpeciesThermo> sp = { mk(0.0280134,3.621,97.53,N2lo,N2hi), mk(0.0180153,2.605,572.4,H2Olo,H2Ohi) };
    for (auto& s : sp) { const double hr = thermo_h_molar(s, 298.15); s.low[7] += -hr/THERMO_RU; s.high[7] += -hr/THERMO_RU; }
    std::vector<SpeciesThermoF> spf = { toF(sp[0]), toF(sp[1]) };
    CondPropOpts o; o.latentLowT=1; o.psatLowT=1; o.liquidCp=2000.0; o.gasKgasModel=0; o.sigmaScale=1.0; o.Yw=0.0;
    const CondSpeciesProps cp = condProps_make(COND_MODEL_H2O, o);
    CondTablesHost ht; cond_tables_build_host(cp, ht); const CondTablesF tb = cond_tables_view_host(ht);
    const double Rw = cp.R;
    double wT = 0.0, wDriftH = 0.0, wDriftD = 0.0; int nfail = 0, n = 0;
    for (double Yw : {0.0113, 0.05, 0.2}) for (double gfrac : {0.001, 0.1, 0.5, 0.99}) for (double T : {150.0, 199.9, 200.1, 220.0, 250.0, 298.15, 350.0, 500.0, 900.0, 1500.0}) {
        const double g = gfrac*Yw; float Yf[2] = {(float)(1.0 - Yw), (float)Yw}; const double Y[2] = {(double)Yf[0], (double)Yf[1]};
        // 参照 e_in = e_mix(T) (double)
        double cpT, hT; thermo_cph_mix(sp.data(), 2, Y, T, &cpT, &hT); const double R = thermo_R_mix(sp.data(), 2, Y);
        const double e_in = (hT - R*T) + g*(Rw*T - cond_latent(cp, T));
        for (double Tg : {0.9*T, 1.1*T, 300.0}) {
            bool ok = true;
            const double Th = cond_T_from_e_twophase_hybrid(sp.data(), spf.data(), 2, Y, Yf, tb, e_in, g, Rw, 1, cp, Tg, 50.0, 6000.0, &ok);
            if (!ok) ++nfail;
            ++n;
            const double eT = fabs(Th - T)/T; if (eT > wT) wT = eT;
            if (eT > 1.0e-8 && n <= 400) { static int nprint = 0; if (nprint++ < 8) { double Gp; const double G = cond_twophase_resid(sp.data(), 2, Y, Th, g, Rw, 1, cp, e_in, &Gp);
                printf("      dbg: T=%.2f Yw=%.4f g=%.2e Tg=%.1f -> Th=%.6f (err %.2e) ok=%d G=%.3e Gp=%.1f tol=%.3e e_in=%.4e\n", T, Yw, g, Tg, Th, eT, (int)ok, G, Gp, 1.0e-9*fabs(e_in)+0.05, e_in); } }
            // 10 往復ドリフト: roe を float で格納 → 反転 → e_mix(T) を再構成 (double) → float 格納 ... (dependentVariables の roe 再構成相当)
            double Tc = Th, Td = cond_T_from_e_carrier(sp.data(), 2, Y, e_in, g, Rw, cp, Tg, 50.0, 6000.0);
            for (int k = 0; k < 10; ++k) {
                bool ok2 = true;
                const float ef = (float)cond_twophase_resid(sp.data(), 2, Y, Tc, g, Rw, 1, cp, 0.0);   // e_mix(Tc) を float 格納
                Tc = cond_T_from_e_twophase_hybrid(sp.data(), spf.data(), 2, Y, Yf, tb, (double)ef, g, Rw, 1, cp, Tc, 50.0, 6000.0, &ok2);
                const float ed = (float)cond_twophase_resid(sp.data(), 2, Y, Td, g, Rw, 1, cp, 0.0);
                Td = cond_T_from_e_carrier(sp.data(), 2, Y, (double)ed, g, Rw, cp, Td, 50.0, 6000.0);
            }
            const double dH = fabs(Tc - T)/T, dD = fabs(Td - T)/T;
            if (dH > wDriftH) wDriftH = dH; if (dD > wDriftD) wDriftD = dD;
        }
    }
    // 意図的な失敗: e_in を T_max=6000 K の外 (T≈9000 K 相当) に置く → 上限に張り付き残差が残る → ok=false (codex result M2)
    { float Yf[2] = {(float)(1.0 - 0.0113), 0.0113f}; const double Y[2] = {(double)Yf[0], (double)Yf[1]};
      double cpT, hT; thermo_cph_mix(sp.data(), 2, Y, 6000.0, &cpT, &hT); const double R = thermo_R_mix(sp.data(), 2, Y);
      const double e_big = (hT - R*6000.0)*1.5;
      bool ok = true; const double Tb = cond_T_from_e_twophase_hybrid(sp.data(), spf.data(), 2, Y, Yf, tb, e_big, 1.0e-3, Rw, 1, cp, 300.0, 50.0, 6000.0, &ok);
      bool ok2 = true; const double Tb2 = cond_twophase_polish(sp.data(), 2, Y, 300.0, e_big, 1.0e-3, Rw, 1, cp, 50.0, 6000.0, &ok2);
      printf("  [%s] forced failure (e beyond T_max): hybrid ok=%d T=%.0f, double-path polish ok=%d T=%.0f (both must be ok=0)\n", (!ok && !ok2) ? "PASS" : "FAIL", (int)ok, Tb, (int)ok2, Tb2);
      if (ok || ok2) ++g_fail; }
    char b[160];
    snprintf(b, sizeof b, "two-phase hybrid |dT|/T (%d cases, %d fail)", n, nfail); check(wT <= 1.0e-8 && nfail == 0, b, wT, 1.0e-8);
    snprintf(b, sizeof b, "two-phase 10-roundtrip drift hybrid %.2e vs double %.2e (limit 1.5x+1e-7)", wDriftH, wDriftD); check(wDriftH <= 1.5*wDriftD + 1.0e-7, b, wDriftH, 1.5*wDriftD + 1.0e-7);
}

int main()
{
    printf("== (1) property tables vs double ==\n");
    CondPropOpts o; o.latentLowT = 1; o.psatLowT = 1; o.liquidCp = 2000.0; o.gasKgasModel = 0; o.sigmaScale = 1.0; o.Yw = 0.0;
    test_tables("N2 (default)", condProps_make(COND_MODEL_N2, o), 15.0, 140.0);
    CondPropOpts o2 = o; o2.latentLowT = 0; o2.psatLowT = 0; o2.gasKgasModel = 1; o2.sigmaScale = 1.03;
    test_tables("N2 (old lowT, air kgas, sigma x1.03)", condProps_make(COND_MODEL_N2, o2), 15.0, 140.0);
    test_tables("H2O", condProps_make(COND_MODEL_H2O, o), 100.0, 1250.0);
    printf("== (2) source chain float vs double ==\n");
    test_chain("N2 pure", condProps_make(COND_MODEL_N2, o), 0, 39.0, 125.5, 1.0);
    test_chain("N2 (old lowT) pure", condProps_make(COND_MODEL_N2, o2), 0, 39.0, 125.5, 1.0);
    test_chain("H2O carrier", condProps_make(COND_MODEL_H2O, o), 1, 150.0, 400.0, 2.5);
    printf("== (3) two-phase hybrid inversion ==\n");
    test_inversion();
    printf("%s (%d failures)\n", g_fail ? "FAILED" : "ALL PASS", g_fail);
    return g_fail ? 1 : 0;
}
