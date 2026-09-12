// float32 熱力学 (SpeciesThermoF / thermo_T_from_e_f) の double 参照に対する単体検証。
// plan performance-3d-node-sst-speedup §4.2-3 (codex M2): 温度誤差・エネルギー残差・反復数・反転→再構成ドリフトを
// 使用 DB (MIXDRY=N2 係数 / H2O / AIR / HE)・組成端点・微量成分・50–6000 K・区間境界・外挿・datum 有無で測る。
// ビルド: g++ -O2 -I solver_density_cuda solver_density_cuda/tools/test_thermo_float.cpp -o /tmp/tthf && /tmp/tthf
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include "../cuda_forge/thermo_d.cuh"

static SpeciesThermo mk(double MW,double sig,double eps,const double lo[9],const double hi[9]){
    SpeciesThermo s; s.MW=MW; s.sigma_LJ=sig; s.eps_kB=eps; s.Tlo=200.0; s.Tmid=1000.0; s.Thi=6000.0; s.h_datum=0.0; s.invMW=1.0/MW;
    for(int i=0;i<9;i++){ s.low[i]=lo[i]; s.high[i]=hi[i]; } return s;
}
static SpeciesThermoF toF(const SpeciesThermo& s){
    SpeciesThermoF f; f.MW=(float)s.MW; f.invMW=(float)(1.0/s.MW); f.R=(float)(THERMO_RU/s.MW);
    f.sigma_LJ=(float)s.sigma_LJ; f.eps_kB=(float)s.eps_kB; f.Tlo=(float)s.Tlo; f.Tmid=(float)s.Tmid; f.Thi=(float)s.Thi;
    for(int k=0;k<9;k++){ f.low[k]=(float)s.low[k]; f.high[k]=(float)s.high[k]; } return f;
}
int main(){
    const double N2lo[9]={2.210371497e+04,-3.818461820e+02,6.082738360e+00,-8.530914410e-03,1.384646189e-05,-9.625793620e-09,2.519705809e-12,7.108460860e+02,-1.076003744e+01};
    const double N2hi[9]={5.877124060e+05,-2.239249073e+03,6.066949220e+00,-6.139685500e-04,1.491806679e-07,-1.923105485e-11,1.061954386e-15,1.283210415e+04,-1.586640027e+01};
    const double H2Olo[9]={-3.947960830e+04,5.755731020e+02,9.317826530e-01,7.222712860e-03,-7.342557370e-06,4.955043490e-09,-1.336933246e-12,-3.303974310e+04,1.724205775e+01};
    const double H2Ohi[9]={1.034972096e+06,-2.412698562e+03,4.646110780e+00,2.291998307e-03,-6.836830480e-07,9.426468930e-11,-4.822380530e-15,-1.384286509e+04,-7.978148510e+00};
    const double AIRc[9]={0.0,0.0,3.5,0.0,0.0,0.0,0.0,-1.0431373e+03,3.0};
    const double HElo[9]={0.0,0.0,2.5,0.0,0.0,0.0,0.0,-7.453750000e+02,9.287239740e-01};
    std::vector<std::pair<std::string,SpeciesThermo>> db = {
        {"N2",mk(0.0280134,3.621,97.53,N2lo,N2hi)}, {"H2O",mk(0.0180153,2.605,572.4,H2Olo,H2Ohi)},
        {"AIR",mk(0.0289647,3.711,78.6,AIRc,AIRc)}, {"HE",mk(0.0040026,2.551,10.22,HElo,HElo)} };
    struct Mix { std::string name; std::vector<int> idx; std::vector<double> Y; };
    std::vector<Mix> mixes = {
        {"N2",{0},{1.0}}, {"H2O",{1},{1.0}}, {"AIR",{2},{1.0}}, {"HE",{3},{1.0}},
        {"N2/H2O 0.989/0.011 (case16)",{0,1},{0.98905,0.01095}}, {"N2/H2O 0.5/0.5",{0,1},{0.5,0.5}},
        {"N2/H2O trace 1e-6",{0,1},{1.0-1e-6,1e-6}}, {"N2/HE 0.5/0.5",{0,3},{0.5,0.5}} };
    // T_max=6000 を超える e は反転不能 (クランプ) なので範囲内のみ。<200 K は低温端の線形外挿 (反転可)。
    const double Tlist_edges[] = {50,60,100,150,199.9,200,200.1,250,298.15,300,400,600,800,950,999.9,1000,1000.1,1200,1500,2000,3000,4000,5000,5900,5999,6000};
    printf("%-30s %-6s %9s %9s %9s %9s %9s %9s %6s %6s\n","mix","datum","errF[K]","errF/T","errD[K]","errHyb/T","maxDrift","maxdh/cpT","itF","itD");
    // errF/errD: 厳密参照 (double Newton, tol 1e-9, 60 反復) に対する float 版 / 生産 double 版 (tol 1e-3+1e-6T) の誤差
    int fails=0;
    for (int datum=0; datum<2; ++datum) {
        for (auto& m : mixes) {
            const int n=(int)m.idx.size();
            std::vector<SpeciesThermo> sp(n); std::vector<SpeciesThermoF> spf(n); std::vector<float> Yf(n);
            // 本番 (dependentVariables useHybrid) と同じく Y は float で組んで正規化し、double 側は (double)Yf を使う
            { float ys=0.f; for (int i=0;i<n;i++){ sp[i]=db[m.idx[i]].second; Yf[i]=(float)m.Y[i]; ys+=Yf[i]; }
              const float inv=1.0f/(ys>1e-30f?ys:1e-30f); for (int i=0;i<n;i++){ Yf[i]*=inv; m.Y[i]=(double)Yf[i]; } }
            if (datum) for (int i=0;i<n;i++){ const double h_ref=thermo_h_molar(sp[i],298.15); const double da7=-h_ref/THERMO_RU; sp[i].low[7]+=da7; sp[i].high[7]+=da7; }
            for (int i=0;i<n;i++) spf[i]=toF(sp[i]);
            double maxdT=0, maxrel=0, maxres=0, maxdrift=0, maxdh=0, maxdTrelT=0, maxresRelT=0, maxHyb=0; int itFmax=0, itDmax=0;
            for (double T : Tlist_edges) {
                // 参照: double で e(T)
                double cpd, hd; thermo_cph_mix(sp.data(), n, m.Y.data(), T, &cpd, &hd);
                const double R=thermo_R_mix(sp.data(),n,m.Y.data()); const double e=hd-R*T;
                // 反転 (double, warm start は T の 0.9 倍・1.1 倍・300K の 3 通り)
                for (double g : {0.9*T, 1.1*T, 300.0, 50.1, 6000.0}) {
                    const double Td = thermo_T_from_e(sp.data(), n, m.Y.data(), e, g, 50.0, 6000.0);
                    int itF=0;
                    const float Tf = thermo_T_from_e_f(spf.data(), n, Yf.data(), (float)e, (float)g, 50.0f, 6000.0f, &itF);
                    // double の反復数も概算 (同じロジックを再実行)
                    int itD=0; { double Tt=g; if(!(Tt>50.0))Tt=50.0; if(Tt>6000.0)Tt=6000.0; for(;itD<20;++itD){ double h_T,cp_T; thermo_cph_mix(sp.data(),n,m.Y.data(),Tt,&cp_T,&h_T); double e_T=h_T-R*Tt, cv=cp_T-R, cvf=(cv>1e-2*R?cv:1e-2*R); double dT=(e_T-e)/cvf; if(dT>0.5*Tt)dT=0.5*Tt; if(dT<-0.5*Tt)dT=-0.5*Tt; Tt-=dT; if(Tt<50.0)Tt=50.0; if(Tt>6000.0)Tt=6000.0; if(dT<0)dT=-dT; if(dT<1e-3+1e-6*Tt){++itD;break;} } }
                    // 厳密参照
                    double Tref=g; { for(int k=0;k<60;++k){ double h_T,cp_T; thermo_cph_mix(sp.data(),n,m.Y.data(),Tref,&cp_T,&h_T); double dT=((h_T-R*Tref)-e)/(cp_T-R); if(dT>0.5*Tref)dT=0.5*Tref; if(dT<-0.5*Tref)dT=-0.5*Tref; Tref-=dT; if(Tref<50.0)Tref=50.0; if(Tref>6000.0)Tref=6000.0; if(fabs(dT)<1e-9)break; } }
                    maxdT=std::max(maxdT,fabs((double)Tf-Tref)); maxrel=std::max(maxrel,fabs(Td-Tref)); maxdTrelT=std::max(maxdTrelT,fabs((double)Tf-Tref)/Tref);
                    // エネルギー残差 (float 解を double で評価)
                    double cp2,h2; thermo_cph_mix(sp.data(),n,m.Y.data(),(double)Tf,&cp2,&h2);
                    maxres=std::max(maxres, fabs((h2-R*(double)Tf)-e)/(cp2-R)); maxresRelT=std::max(maxresRelT, fabs((h2-R*(double)Tf)-e)/(cp2-R)/Tref);
                    // 反転→再構成ドリフト: float で e_v(Tf) を組み直し、もう一度反転して T がどれだけ動くか (dependentVariables の roe 再構成相当)
                    float cpf,hf; thermo_cph_mix_f(spf.data(),n,Yf.data(),Tf,&cpf,&hf);
                    const float Rf=thermo_R_mix_f(spf.data(),n,Yf.data()); const float e2=hf-Rf*Tf;
                    const float Tf2=thermo_T_from_e_f(spf.data(),n,Yf.data(),e2,Tf,50.0f,6000.0f,nullptr);
                    maxdrift=std::max(maxdrift,(double)fabsf(Tf2-Tf));
                    // ハイブリッド (float 8 反復 + double 研磨 1 段) の誤差
                    { double cpH,hH,TfH; const double Th=thermo_T_from_e_hybrid(sp.data(),spf.data(),n,m.Y.data(),Yf.data(),e,g,50.0,6000.0,&cpH,&hH,&TfH,12);
                      maxHyb=std::max(maxHyb,fabs(Th-Tref)/Tref);
                      // 再格納ドリフト: h(T)=h(T_f)+cp·(T−T_f) (本番の Taylor 再構成) から e を組み直し、もう一度反転して T の動きを見る
                      const double e2=(hH+cpH*(Th-TfH))-thermo_R_mix(sp.data(),n,m.Y.data())*Th;
                      double cp2,h2,Tf2; const double Th2=thermo_T_from_e_hybrid(sp.data(),spf.data(),n,m.Y.data(),Yf.data(),e2,Th,50.0,6000.0,&cp2,&h2,&Tf2,12);
                      maxdrift=std::max(maxdrift,fabs(Th2-Th)/Th); }
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
            const bool ok = (maxHyb < 3.0e-8);
            if (!ok) fails++;
            printf("%-30s %-6s %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %6d %6d %s\n", m.name.c_str(), datum?"298K":"abs", maxdT, maxdTrelT, maxrel, maxHyb, maxdrift, maxdh, itFmax, itDmax, ok?"OK":"FAIL");
        }
    }
    printf("VERDICT: %s (fails=%d; 判定 ハイブリッド反転の誤差 errHyb/T < 3e-8 [float の T 格納分解能未満])\n", fails?"FAIL":"PASS", fails);
    return fails?1:0;
}
