// M6 検証: FP32 化した kinetic-theory 輸送係数 (thermo_mu/lambda_mix) を
// NIST 実測値および double 参照と照合する単体テスト。
// __CUDACC__ 未定義のため THERMO_HD は inline 展開され、純 host C++ でビルド可。
#include <cstdio>
#include <cmath>
#include "../cuda_forge/thermo_d.cuh"

// builtinDB() (thermo_d.cu) と同一の係数で N2 / AIR / He を再構成
static SpeciesThermo mk(double MW,double sig,double eps,const double lo[9],const double hi[9]){
    SpeciesThermo s; s.MW=MW; s.sigma_LJ=sig; s.eps_kB=eps;
    thermo_set_nasa9_2(s, 200.0, 1000.0, 6000.0, lo, hi);
    return s;
}
// double 参照: FP32 化前と同じ Chapman-Enskog/Eucken を全 double で評価
static double omega22_d(double Ts){
    if(Ts<0.3)Ts=0.3; if(Ts>100.0)Ts=100.0;
    return 1.16145*pow(Ts,-0.14874)+0.52487*exp(-0.77320*Ts)+2.16178*exp(-2.43787*Ts);
}
static double mu_d(const SpeciesThermo&s,double T){
    double Ts=T/s.eps_kB, Mg=s.MW*1000.0, sig2=s.sigma_LJ*s.sigma_LJ;
    return 2.6693e-6*sqrt(Mg*T)/(sig2*omega22_d(Ts));
}
static double lam_d(const SpeciesThermo&s,double T){
    return mu_d(s,T)*(thermo_cp_mass(s,T)+1.25*thermo_R_species(s));
}

int main(){
    const double N2lo[9]={2.210371497e+04,-3.818461820e+02,6.082738360e+00,-8.530914410e-03,1.384646189e-05,-9.625793620e-09,2.519705809e-12,7.108460860e+02,-1.076003744e+01};
    const double N2hi[9]={5.877124060e+05,-2.239249073e+03,6.066949220e+00,-6.139685500e-04,1.491806679e-07,-1.923105485e-11,1.061954386e-15,1.283210415e+04,-1.586640027e+01};
    const double AIRc[9]={0.0,0.0,3.5,0.0,0.0,0.0,0.0,-1.0431373e+03,3.0};
    const double Helo[9]={0.0,0.0,2.5,0.0,0.0,0.0,0.0,-7.453750000e+02,9.287239740e-01};
    SpeciesThermo N2=mk(0.0280134,3.621,97.53,N2lo,N2hi);
    SpeciesThermo AIR=mk(0.0289647,3.711,78.6,AIRc,AIRc);
    SpeciesThermo He=mk(0.0040026,2.551,10.22,Helo,Helo);

    const double X1[1]={1.0};
    struct Case{const char*name;SpeciesThermo*sp;double T;const char*q;double nist;};
    Case cs[]={
        {"N2 mu(300)",   &N2,300.0,"mu", 1.78e-5},
        {"N2 mu(1000)",  &N2,1000.0,"mu",4.15e-5},
        {"N2 lam(300)",  &N2,300.0,"lam",0.0260},
        {"AIR mu(300)",  &AIR,300.0,"mu",1.85e-5},
        {"AIR lam(300)", &AIR,300.0,"lam",0.0263},
        {"He mu(300)",   &He,300.0,"mu",1.99e-5},
    };
    printf("%-14s %12s %12s %12s %10s %10s\n","case","FP32","double","NIST","rel(d-32)","rel(NIST)");
    bool ok=true;
    for(auto&c:cs){
        bool ismu=(c.q[0]=='m');
        double f32 = ismu ? (double)thermo_mu_mix(c.sp,1,X1,c.T) : (double)thermo_lambda_mix(c.sp,1,X1,c.T);
        double dd  = ismu ? mu_d(*c.sp,c.T) : lam_d(*c.sp,c.T);
        double rd  = fabs(f32-dd)/dd;
        double rn  = fabs(f32-c.nist)/c.nist;
        printf("%-14s %12.4e %12.4e %12.4e %9.2e%% %9.1f%%\n",c.name,f32,dd,c.nist,rd*100,rn*100);
        if(rd>1e-4) ok=false;            // FP32 vs double: 1e-4 以内であること
        if(rn>0.12) ok=false;            // NIST: kinetic theory 近似として 12% 以内
    }
    // binary He/N2 混合 (Wilke) も NaN/正値チェック
    SpeciesThermo mix[2]={He,N2};
    double Xm[2]={0.5,0.5};
    double mmix=(double)thermo_mu_mix(mix,2,Xm,300.0);
    double lmix=(double)thermo_lambda_mix(mix,2,Xm,300.0);
    printf("He/N2 50/50 @300: mu=%.4e lam=%.4e\n",mmix,lmix);
    if(!(mmix>0&&lmix>0&&std::isfinite(mmix)&&std::isfinite(lmix))) ok=false;
    printf("\nRESULT: %s\n",ok?"PASS":"FAIL");
    return ok?0:1;
}
