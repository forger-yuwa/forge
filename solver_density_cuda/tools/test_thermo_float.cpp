// float32 熱力学 (SpeciesThermoF / thermo_T_from_e_f) の double 参照に対する単体検証。
// plan performance-3d-node-sst-speedup §4.2-3 (codex M2): 温度誤差・エネルギー残差・反復数・反転→再構成ドリフトを
// 使用 DB (MIXDRY=N2 係数 / H2O / AIR / HE)・組成端点・微量成分・50–6000 K・区間境界・外挿・datum 有無で測る。
// 区間可変 (plan thermophysics-solver-owned-species-db #13-1 G1-d): CEA の 3 区間 N2/O2 (200/1000/6000/20000) を含む組成を 20000 K まで
// (反転の上限 T_max も 20000 K; 6000 K の区切りの前後・20000 K 端を含む) 同じ判定で測る。
// 判定の再スコープ (plan §6 V3f, 2026-10-01 diagnostician; 事前固定):
//   (i)  冷間開始 (50.1 K) と T_max 開始は目標 T <= 6000 K に限る (1 呼び出しの到達上限 50.1·1.5^12 ≈ 6500 K の設計値)。
//        T > 6000 K の目標は warm start 3 通り (0.9T・1.1T・300 K) で判定。errHyb/T < 3e-8 は不変。
//   (ii) 区切り温度そのもの (組成のいずれかの種の内側の区切り; 1000 K・6000 K) では許容を 3e-8 + Δh_step/(c_v·T)
//        (Δh_step = 係数から計算した両側の混合 h の差 (上の区間 − 下の区間, 同じ T)、c_v はその点の混合 c_v)。区切り以外は 3e-8。
//   info (合否にしない): 6000 K 超での float 面経路 thermo_h_mix_f の h 相対誤差 (double 評価比) の最大 hF/h>6k。
// ビルド: g++ -O2 -I solver_density_cuda solver_density_cuda/tools/test_thermo_float.cpp -o /tmp/tthf && /tmp/tthf
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "../cuda_forge/thermo_d.cuh"

// 区間 k の係数で強制評価した h [J/mol] (区切りでの段差用)
static double h_molar_k(const SpeciesThermo& s, int k, double T){
    const double* a=s.coef[k]; const double Ti=1.0/T, lnT=log(T);
    return THERMO_RU*T*(-a[0]*Ti*Ti + a[1]*lnT*Ti + a[2] + a[3]*T/2.0 + a[4]*T*T/3.0 + a[5]*T*T*T/4.0 + a[6]*T*T*T*T/5.0 + a[7]*Ti);
}
// 組成 Y の混合 h の区切り段差 Δh_step [J/kg] (T がどの種の内側の区切りでもなければ 0)
static double h_step_mix(const std::vector<SpeciesThermo>& sp, const std::vector<double>& Y, double T, bool* isBrk){
    double d=0.0; *isBrk=false;
    for (size_t i=0;i<sp.size();++i) for (int k=0;k+1<sp[i].nInt;++k) if (sp[i].Tbrk[k]==T) { *isBrk=true; d+=Y[i]*(h_molar_k(sp[i],k+1,T)-h_molar_k(sp[i],k,T))/sp[i].MW; }
    return fabs(d);
}
static SpeciesThermo mk3(double MW,double sig,double eps,const double a[3][9]){
    SpeciesThermo s; s.MW=MW; s.sigma_LJ=sig; s.eps_kB=eps; s.h_datum=0.0; s.invMW=1.0/MW;
    const double Tb[4]={200.0,1000.0,6000.0,20000.0}; thermo_set_intervals(s,3,Tb,a); return s;
}
static SpeciesThermo mk(double MW,double sig,double eps,const double lo[9],const double hi[9]){
    SpeciesThermo s; s.MW=MW; s.sigma_LJ=sig; s.eps_kB=eps; s.h_datum=0.0; s.invMW=1.0/MW;
    thermo_set_nasa9_2(s, 200.0, 1000.0, 6000.0, lo, hi); return s;
}
static SpeciesThermoF toF(const SpeciesThermo& s){
    return thermo_to_float(s);
}
int main(){
    const double N2lo[9]={2.210371497e+04,-3.818461820e+02,6.082738360e+00,-8.530914410e-03,1.384646189e-05,-9.625793620e-09,2.519705809e-12,7.108460860e+02,-1.076003744e+01};
    const double N2hi[9]={5.877124060e+05,-2.239249073e+03,6.066949220e+00,-6.139685500e-04,1.491806679e-07,-1.923105485e-11,1.061954386e-15,1.283210415e+04,-1.586640027e+01};
    const double H2Olo[9]={-3.947960830e+04,5.755731020e+02,9.317826530e-01,7.222712860e-03,-7.342557370e-06,4.955043490e-09,-1.336933246e-12,-3.303974310e+04,1.724205775e+01};
    const double H2Ohi[9]={1.034972096e+06,-2.412698562e+03,4.646110780e+00,2.291998307e-03,-6.836830480e-07,9.426468930e-11,-4.822380530e-15,-1.384286509e+04,-7.978148510e+00};
    const double AIRc[9]={0.0,0.0,3.5,0.0,0.0,0.0,0.0,-1.0431373e+03,3.0};
    const double HElo[9]={0.0,0.0,2.5,0.0,0.0,0.0,0.0,-7.453750000e+02,9.287239740e-01};
    // CEA thermo.inp の 3 区間 (200/1000/6000/20000)
    const double N2c3[3][9]={{2.210371497e+04,-3.818461820e+02,6.082738360e+00,-8.530914410e-03,1.384646189e-05,-9.625793620e-09,2.519705809e-12,7.108460860e+02,-1.076003744e+01},
                             {5.877124060e+05,-2.239249073e+03,6.066949220e+00,-6.139685500e-04,1.491806679e-07,-1.923105485e-11,1.061954386e-15,1.283210415e+04,-1.586640027e+01},
                             {831013916.0,-642073.354,202.0264635,-0.03065092046,2.486903333e-06,-9.70595411e-11,1.437538881e-15,4938707.04,-1672.09974}};
    const double O2c3[3][9]={{-34255.6342,484.700097,1.119010961,0.00429388924,-6.83630052e-07,-2.0233727e-09,1.039040018e-12,-3391.45487,18.4969947},
                             {-1037939.022,2344.830282,1.819732036,0.001267847582,-2.188067988e-07,2.053719572e-11,-8.19346705e-16,-16890.10929,17.38716506},
                             {497529430.0,-286610.6874,66.9035225,-0.00616995902,3.016396027e-07,-7.4214166e-12,7.27817577e-17,2293554.027,-553.062161}};
    std::vector<std::pair<std::string,SpeciesThermo>> db = {
        {"N2",mk(0.0280134,3.621,97.53,N2lo,N2hi)}, {"H2O",mk(0.0180153,2.605,572.4,H2Olo,H2Ohi)},
        {"AIR",mk(0.0289647,3.711,78.6,AIRc,AIRc)}, {"HE",mk(0.0040026,2.551,10.22,HElo,HElo)},
        {"N2c3",mk3(0.0280134,3.621,97.53,N2c3)}, {"O2c3",mk3(0.0319988,3.458,107.4,O2c3)} };
    struct Mix { std::string name; std::vector<int> idx; std::vector<double> Y; double Tmax; };
    std::vector<Mix> mixes = {
        {"N2",{0},{1.0},6000.0}, {"H2O",{1},{1.0},6000.0}, {"AIR",{2},{1.0},6000.0}, {"HE",{3},{1.0},6000.0},
        {"N2/H2O 0.989/0.011 (case16)",{0,1},{0.98905,0.01095},6000.0}, {"N2/H2O 0.5/0.5",{0,1},{0.5,0.5},6000.0},
        {"N2/H2O trace 1e-6",{0,1},{1.0-1e-6,1e-6},6000.0}, {"N2/HE 0.5/0.5",{0,3},{0.5,0.5},6000.0},
        // 区間可変 (#13-1 G1-d): 3 区間種を含む組成を 20000 K まで
        {"N2c3 (3-int, to 20000 K)",{4},{1.0},20000.0}, {"N2c3/O2c3 0.767/0.233",{4,5},{0.767,0.233},20000.0},
        {"N2c3/H2O 0.989/0.011 (H2O extrap>6000)",{4,1},{0.98905,0.01095},20000.0}, {"N2c3/HE 0.5/0.5",{4,3},{0.5,0.5},20000.0} };
    // T_max を超える e は反転不能 (クランプ) なので範囲内のみ。<200 K は低温端の線形外挿 (反転可)。
    const double Tlist_edges[] = {50,60,100,150,199.9,200,200.1,250,298.15,300,400,600,800,950,999.9,1000,1000.1,1200,1500,2000,3000,4000,5000,5900,5999,6000,
                                  6000.1,6500,7000,8000,10000,12000,15000,18000,19000,19999,20000};
    printf("%-30s %-6s %9s %9s %9s %9s %9s %9s %9s %9s %6s\n","mix","datum","errF[K]","errF/T","errD[K]","errHyb/T","hyb/tol","driftH/T","driftD/T","driftF/T","itF");
    // errF/errD: 厳密参照 (double Newton, tol 1e-9, 60 反復) に対する float 版 / 生産 double 版 (tol 1e-3+1e-6T) の誤差
    int fails=0;
    for (int datum=0; datum<2; ++datum) {
        for (auto& m : mixes) {
            const int n=(int)m.idx.size();
            std::vector<SpeciesThermo> sp(n); std::vector<SpeciesThermoF> spf(n); std::vector<float> Yf(n);
            // 本番 (dependentVariables useHybrid) と同じく Y は float で組んで正規化し、double 側は (double)Yf を使う
            { float ys=0.f; for (int i=0;i<n;i++){ sp[i]=db[m.idx[i]].second; Yf[i]=(float)m.Y[i]; ys+=Yf[i]; }
              const float inv=1.0f/(ys>1e-30f?ys:1e-30f); for (int i=0;i<n;i++){ Yf[i]*=inv; m.Y[i]=(double)Yf[i]; } }
            if (datum) for (int i=0;i<n;i++){ const double h_ref=thermo_h_molar(sp[i],298.15); const double da7=-h_ref/THERMO_RU; thermo_add_a7(sp[i], da7); }
            for (int i=0;i<n;i++) spf[i]=toF(sp[i]);
            double maxdT=0, maxrel=0, maxres=0, maxdrift=0, maxdh=0, maxdTrelT=0, maxresRelT=0, maxHyb=0, maxdriftD=0, maxdriftF=0; int itFmax=0, itDmax=0;
            double hybT=0, hybG=0;   // errHyb/T の最悪点 (T, 初期値) — 20000 K までの 3 区間種で原因を切り分けるため (#13-1 G1-d)
            double maxRatio=0, ratT=0, ratTol=0;   // (errHyb/T)/許容 の最大 (許容は点ごと; 再スコープ (ii))
            double maxHF=-1;                       // info: 6000 K 超の thermo_h_mix_f 相対誤差 (-1 = 該当点なし)
            const double TM = m.Tmax;
            for (double T : Tlist_edges) {
                if (T > TM) continue;
                // 参照: double で e(T)
                double cpd, hd; thermo_cph_mix(sp.data(), n, m.Y.data(), T, &cpd, &hd);
                const double R=thermo_R_mix(sp.data(),n,m.Y.data()); const double e=hd-R*T;
                // 反転 (double, warm start は T の 0.9 倍・1.1 倍・300K の 3 通り)
                // 再スコープ (ii): 区切りそのものの許容
                bool isBrk=false; double tolT=3.0e-8;
                { const double dstep=h_step_mix(sp,m.Y,T,&isBrk);
                  if (isBrk) { double cpb,hb; thermo_cph_mix(sp.data(),n,m.Y.data(),T,&cpb,&hb); const double cvb=cpb-thermo_R_mix(sp.data(),n,m.Y.data()); tolT += dstep/(cvb*T); } }
                // info: float 面経路の h (6000 K 超)
                if (T > 6000.0) { const float Tf32=(float)T; double cpr,hr; thermo_cph_mix(sp.data(),n,m.Y.data(),(double)Tf32,&cpr,&hr);
                  const float hF=thermo_h_mix_f(spf.data(),n,Yf.data(),Tf32); maxHF=std::max(maxHF, fabs((double)hF-hr)/fabs(hr)); }
                // 再スコープ (i): 冷間開始 (50.1 K) と T_max 開始は目標 <= 6000 K のみ、6000 K 超は warm start 3 通り
                std::vector<double> starts = {0.9*T, 1.1*T, 300.0};
                if (T <= 6000.0) { starts.push_back(50.1); starts.push_back(TM); }
                for (double g : starts) {
                    const double Td = thermo_T_from_e(sp.data(), n, m.Y.data(), e, g, 50.0, TM);
                    int itF=0;
                    const float Tf = thermo_T_from_e_f(spf.data(), n, Yf.data(), (float)e, (float)g, 50.0f, (float)TM, &itF);
                    // double の反復数も概算 (同じロジックを再実行)
                    int itD=0; { double Tt=g; if(!(Tt>50.0))Tt=50.0; if(Tt>TM)Tt=TM; for(;itD<20;++itD){ double h_T,cp_T; thermo_cph_mix(sp.data(),n,m.Y.data(),Tt,&cp_T,&h_T); double e_T=h_T-R*Tt, cv=cp_T-R, cvf=(cv>1e-2*R?cv:1e-2*R); double dT=(e_T-e)/cvf; if(dT>0.5*Tt)dT=0.5*Tt; if(dT<-0.5*Tt)dT=-0.5*Tt; Tt-=dT; if(Tt<50.0)Tt=50.0; if(Tt>TM)Tt=TM; if(dT<0)dT=-dT; if(dT<1e-3+1e-6*Tt){++itD;break;} } }
                    // 厳密参照
                    double Tref=g; { for(int k=0;k<60;++k){ double h_T,cp_T; thermo_cph_mix(sp.data(),n,m.Y.data(),Tref,&cp_T,&h_T); double dT=((h_T-R*Tref)-e)/(cp_T-R); if(dT>0.5*Tref)dT=0.5*Tref; if(dT<-0.5*Tref)dT=-0.5*Tref; Tref-=dT; if(Tref<50.0)Tref=50.0; if(Tref>TM)Tref=TM; if(fabs(dT)<1e-9)break; } }
                    maxdT=std::max(maxdT,fabs((double)Tf-Tref)); maxrel=std::max(maxrel,fabs(Td-Tref)); maxdTrelT=std::max(maxdTrelT,fabs((double)Tf-Tref)/Tref);
                    // エネルギー残差 (float 解を double で評価)
                    double cp2,h2; thermo_cph_mix(sp.data(),n,m.Y.data(),(double)Tf,&cp2,&h2);
                    maxres=std::max(maxres, fabs((h2-R*(double)Tf)-e)/(cp2-R)); maxresRelT=std::max(maxresRelT, fabs((h2-R*(double)Tf)-e)/(cp2-R)/Tref);
                    // 反転→再構成ドリフト: float で e_v(Tf) を組み直し、もう一度反転して T がどれだけ動くか (dependentVariables の roe 再構成相当)
                    float cpf,hf; thermo_cph_mix_f(spf.data(),n,Yf.data(),Tf,&cpf,&hf);
                    const float Rf=thermo_R_mix_f(spf.data(),n,Yf.data()); const float e2=hf-Rf*Tf;
                    const float Tf2=thermo_T_from_e_f(spf.data(),n,Yf.data(),e2,Tf,50.0f,(float)TM,nullptr);
                    maxdriftF=std::max(maxdriftF,(double)fabsf(Tf2-Tf)/Tref);   // 純 float 反転の再反転ドリフト (相対, 参考値)
                    // ハイブリッド (float 8 反復 + double 研磨 1 段) の誤差
                    { double cpH,hH,TfH; const double Th=thermo_T_from_e_hybrid(sp.data(),spf.data(),n,m.Y.data(),Yf.data(),e,g,50.0,TM,&cpH,&hH,&TfH,12);
                      if (fabs(Th-Tref)/Tref > maxHyb) { maxHyb=fabs(Th-Tref)/Tref; hybT=T; hybG=g; }
                      if (fabs(Th-Tref)/Tref/tolT > maxRatio) { maxRatio=fabs(Th-Tref)/Tref/tolT; ratT=T; ratTol=tolT; }
                      // 再格納ドリフト (本番 dependentVariables と同じ経路): h(T)=h(T_f)+cp·(T−T_f) から e_mix を組み、
                      // roe=ρ(e_mix+ek) を **float に格納**して読み戻し、e=roe/ρ−ek を float→double で再反転する。
                      // 反復 (10 回) で T が漂わないこと・有限であることを合否に含める (codex result-2 m5)。
                      { const double R=thermo_R_mix(sp.data(),n,m.Y.data()); const float rho=1.0f; const float ek=0.5f*300.0f*300.0f;
                        double Tk=Th, cpk=cpH, hk=hH, Tfk=TfH; double Tfirst=Th;
                        double Td_k=Td;                                          // 同じ float 格納往復を従来 double 反転で
                        for (int it=0; it<10; ++it) {
                          const double e_mix=(hk+cpk*(Tk-Tfk))-R*Tk;
                          const float roe=(float)(rho*(e_mix+(double)ek));      // float 格納 (本番)
                          const float e_re=roe/rho-ek;                           // float 読み戻し (本番 intE)
                          double cp2,h2,Tf2; const double Tn=thermo_T_from_e_hybrid(sp.data(),spf.data(),n,m.Y.data(),Yf.data(),(double)e_re,Tk,50.0,TM,&cp2,&h2,&Tf2,12);
                          if (!std::isfinite(Tn)) { maxdrift=1e9; break; }
                          Tk=Tn; cpk=cp2; hk=h2; Tfk=Tf2;
                          double cpd,hd2; thermo_cph_mix(sp.data(),n,m.Y.data(),Td_k,&cpd,&hd2);
                          const float roed=(float)(rho*((hd2-R*Td_k)+(double)ek)); const float e_red=roed/rho-ek;
                          Td_k=thermo_T_from_e(sp.data(),n,m.Y.data(),(double)e_red,Td_k,50.0,TM);
                        }
                        maxdrift=std::max(maxdrift,fabs(Tk-Tfirst)/Tfirst);
                        maxdriftD=std::max(maxdriftD,fabs(Td_k-Td)/Td); } }
                    maxdh=std::max(maxdh, fabs((double)hf-hd)/(cpd*T));
                    itFmax=std::max(itFmax,itF); itDmax=std::max(itDmax,itD);
                }
            }
            // 合否: float の誤差が 2e-3 K 未満、かつ生産 double 版の誤差の 2 倍 + 1e-4 K 以内 (float 固有の劣化がない)
            // 合否: float の誤差 (厳密参照比) が T の 1e-6 (≈8 ulp) 未満、エネルギー残差/cv も同様、20 反復張り付き無し。
            // (生産 double 版は 1e-10 K まで落ちるので double との差ではなく float 分解能を基準にする。)
            // 合否 (採用経路はハイブリッド): errHyb/T < 3e-8 = float の T 格納分解能 (ulp 6e-8) の半分未満。
            // 1 段研磨後の残差は NASA 係数の温度域境界 (Tmid/Tlo) を float 段が跨ぐときの cp 段差由来で ~1e-8 まで残る。
            // float 単独 (errF/T ~1e-6, abs datum H2O は 20 反復張り付き) は参考値。
            // 再格納ドリフト (10 往復, float 格納込み) は float の e 分解能 (ulp/e ≈ 6e-8 → ΔT/T ≈ 3e-8·(e/(cv T)) ≲ 1e-7) の範囲: < 3e-7·T。
            // abs datum (thermoHrefTemp 無し) は本番でハイブリッドを使わない (自動で double 反転) ので、ドリフトは参考値扱い。
            // ドリフトは float 格納そのものに由来する (従来 double 反転でも同程度)。判定: ハイブリッドのドリフトが double 反転の 2 倍 + 1e-8 以内。
            const bool ok = (maxRatio < 1.0) && (datum == 0 || maxdrift <= 2.0*maxdriftD + 1.0e-8);   // errHyb/T < 許容 (区切り以外 3e-8)
            if (!ok) fails++;
            char hfs[32]; if (maxHF<0) snprintf(hfs,sizeof(hfs),"-"); else snprintf(hfs,sizeof(hfs),"%.2e",maxHF);
            printf("%-30s %-6s %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %6d %s  (errHyb worst at T=%g, start %g; worst ratio at T=%g tol %.3e; info hF/h>6k %s)\n", m.name.c_str(), datum?"298K":"abs", maxdT, maxdTrelT, maxrel, maxHyb, maxRatio, maxdrift, maxdriftD, maxdriftF, itFmax, ok?"OK":"FAIL", hybT, hybG, ratT, ratTol, hfs);
        }
    }
    printf("VERDICT: %s (fails=%d; 判定 errHyb/T < 3e-8 [float の T 格納分解能未満; 区切りそのものは 3e-8+Δh_step/(c_v T), 6000 K 超は warm start のみ] かつ float 格納 roe 10 往復のハイブリッド反転ドリフト driftH/T が従来 double 反転 driftD/T の 2 倍+1e-8 以内; driftF は純 float 反転の参考値)\n", fails?"FAIL":"PASS", fails);
    return fails?1:0;
}
