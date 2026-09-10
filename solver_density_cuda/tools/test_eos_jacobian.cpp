// 一般EOS flux Jacobian 固有系の単体検証 (Level1/2)。
// plan time_integration-general-eos-jacobian.md §6。__CUDACC__ 未定義で純 host ビルド可。
//   Level1: ||LR-I||, ||RΛL - A_FD||/||A_FD||  (CPG/TP, 複数温度/Mach/法線)
//   Level2: A+ + A- = A,  法線反転  A_{-n}^+ = -A_n^- ,  A_{-n}^- = -A_n^+
// ビルド: g++ -O2 -I solver_density_cuda solver_density_cuda/tools/test_eos_jacobian.cpp -o /tmp/teij
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "../cuda_forge/thermo_d.cuh"
#include "../cuda_forge/eos_jacobian_d.cuh"
#include "../cuda_forge/block_dplur_jacobian_d.cuh"  // Level3: 実カーネル関数 accumulate_split_jacobian_cf
#include "../cuda_forge/condensationEOS_d.cuh"        // mode 2: 一温度二相 EOS (固定 g,Y) の frozen 音速/κ

static SpeciesThermo mkN2(){
    const double lo[9]={2.210371497e+04,-3.818461820e+02,6.082738360e+00,-8.530914410e-03,1.384646189e-05,-9.625793620e-09,2.519705809e-12,7.108460860e+02,-1.076003744e+01};
    const double hi[9]={5.877124060e+05,-2.239249073e+03,6.066949220e+00,-6.139685500e-04,1.491806679e-07,-1.923105485e-11,1.061954386e-15,1.283210415e+04,-1.586640027e+01};
    SpeciesThermo s; s.MW=0.0280134; s.sigma_LJ=3.621; s.eps_kB=97.53;
    s.Tlo=200.0; s.Tmid=1000.0; s.Thi=6000.0;
    for(int i=0;i<9;i++){ s.low[i]=lo[i]; s.high[i]=hi[i]; }
    return s;
}
static SpeciesThermo N2 = mkN2();
static const double Y1[1]={1.0};
// mode 2 (plans/active/condensation-kantrowitz-gamma-twophase-sonic.md §6 g): TP carrier N2+H2O, 液相分率 g 固定の二相 EOS。
static SpeciesThermo mkH2O(){
    const double lo[9]={-3.947960830e+04,5.755731020e+02,9.317826530e-01,7.222712860e-03,-7.342557370e-06,4.955043490e-09,-1.336933246e-12,-3.303974310e+04,1.724205775e+01};
    const double hi[9]={1.034972096e+06,-2.412698562e+03,4.646110780e+00,2.291998307e-03,-6.836830480e-07,9.426468930e-11,-4.821580530e-15,-1.384286509e+04,-7.978148510e+00};
    SpeciesThermo s; s.MW=0.0180153; s.sigma_LJ=2.605; s.eps_kB=572.4;
    s.Tlo=200.0; s.Tmid=1000.0; s.Thi=6000.0;
    for(int i=0;i<9;i++){ s.low[i]=lo[i]; s.high[i]=hi[i]; }
    return s;
}
static SpeciesThermo SP2[2] = {mkN2(), mkH2O()};
static const double Y2[2] = {0.98905, 0.01095};
static const double G_FIX = 0.0090;          // 液相質量分率 (固定)
static const CondSpeciesProps H2O = condProps_H2O();

struct Prim { double ro,ux,uy,uz,P,T,c,h,Ht,kappa,chi,e; };

// mode: 0=CPG (gamma,cp), 1=TP (NASA N2), 2=二相 TP carrier (N2+H2O, g=G_FIX 固定)
static Prim prim_from_Q(const double Q[5], int mode){
    Prim p; p.ro=Q[0]; p.ux=Q[1]/Q[0]; p.uy=Q[2]/Q[0]; p.uz=Q[3]/Q[0];
    double ek=0.5*(p.ux*p.ux+p.uy*p.uy+p.uz*p.uz);
    p.e = Q[4]/Q[0] - ek;
    if (mode==0){
        const double gam=1.4, cp=1039.0; double R=cp*(gam-1.0)/gam, cv=cp-R;
        p.T=p.e/cv; p.P=p.ro*R*p.T; p.h=cp*p.T; p.c=sqrt(gam*R*p.T); p.kappa=gam-1.0;
    } else if (mode==1) {
        p.T = thermo_T_from_e(&N2,1,Y1,p.e, 300.0, 50.0, 6000.0);
        ThermoDerivatives D; thermo_derivatives_mix(&N2,1,Y1,p.T,&D);
        p.P=p.ro*D.R*p.T; p.h=D.h; p.c=sqrt(D.a2); p.kappa=D.kappa;
    } else {
        const double g=G_FIX, Rw=H2O.R;
        p.T = cond_T_from_e_carrier(SP2,2,Y2,p.e,g,Rw,H2O,250.0,50.0,3000.0);
        double cpm,hm; thermo_cph_mix(SP2,2,Y2,p.T,&cpm,&hm);
        const double Rmix=thermo_R_mix(SP2,2,Y2), Reff=Rmix-g*Rw;
        p.P=p.ro*Reff*p.T; p.h=p.e+p.P/p.ro;   // 二相比エンタルピー h=e_mix+p/ρ
        const double dL=(cond_latent(H2O,p.T+0.1)-cond_latent(H2O,p.T-0.1))/0.2;
        double g2,c2; cond_twophase_sonic(cpm,Reff,g,dL,p.T,&g2,&c2);
        p.c=sqrt(c2); p.kappa=g2-1.0;
    }
    p.Ht = p.h + ek;
    p.chi = p.c*p.c - p.kappa*p.h;   // = ∂P/∂ρ|_{ρe} (CPG で 0)
    return p;
}
static void Fn(const double Q[5], const double n[3], int mode, double F[5]){
    Prim p=prim_from_Q(Q,mode); double Un=p.ux*n[0]+p.uy*n[1]+p.uz*n[2];
    F[0]=p.ro*Un;
    F[1]=p.ro*p.ux*Un+p.P*n[0];
    F[2]=p.ro*p.uy*Un+p.P*n[1];
    F[3]=p.ro*p.uz*Un+p.P*n[2];
    F[4]=(Q[4]+p.P)*Un;
}
static void A_fd(const double Q[5], const double n[3], int mode, double A[5][5]){
    for(int j=0;j<5;j++){
        double h=1e-6*std::max(std::fabs(Q[j]),1e-3);
        double Qp[5],Qm[5],Fp[5],Fm[5];
        for(int k=0;k<5;k++){Qp[k]=Q[k];Qm[k]=Q[k];} Qp[j]+=h; Qm[j]-=h;
        Fn(Qp,n,mode,Fp); Fn(Qm,n,mode,Fm);
        for(int i=0;i<5;i++) A[i][j]=(Fp[i]-Fm[i])/(2*h);
    }
}
static double fro(double A[5][5]){ double s=0; for(int i=0;i<5;i++)for(int j=0;j<5;j++)s+=A[i][j]*A[i][j]; return sqrt(s); }
static double frodiff(double A[5][5],double B[5][5]){ double s=0; for(int i=0;i<5;i++)for(int j=0;j<5;j++){double d=A[i][j]-B[i][j];s+=d*d;} return sqrt(s); }

// 全 Jacobian A = R Λ L (split せず) と L=R^-1 を作り、||LR-I|| を返す。
static double build_A_and_LRerr(const Prim&p,const double n[3],double A[5][5]){
    double R[5][5],L[5][5],lam[5];
    eos_eigvecs_general(p.ux,p.uy,p.uz,n[0],n[1],n[2],p.c,p.Ht,p.kappa,R,lam);
    eos_inv5(R,L);
    for(int i=0;i<5;i++)for(int j=0;j<5;j++){double s=0;for(int k=0;k<5;k++)s+=R[i][k]*lam[k]*L[k][j];A[i][j]=s;}
    // ||LR-I||
    double e=0; for(int i=0;i<5;i++)for(int j=0;j<5;j++){double s=0;for(int k=0;k<5;k++)s+=L[i][k]*R[k][j];double d=s-(i==j?1.0:0.0);e+=d*d;}
    return sqrt(e);
}

int main(){
    struct Case{const char*name;int mode;double T;double M;double dir;};
    // dir: 0=軸平行 (1,0,0)、1=一般方向
    Case cs[]={
        {"CPG T250 M1.5 axis", 0,250,1.5,0},
        {"CPG T250 M1.5 genN", 0,250,1.5,1},
        {"TP  T250 M0.0",      1,250,0.0,1},
        {"TP  T250 M0.5",      1,250,0.5,1},
        {"TP  T250 M0.98",     1,250,0.98,1},
        {"TP  T250 M1.5",      1,250,1.5,1},
        {"TP  T250 M3.0",      1,250,3.0,1},
        {"TP  T1200 M1.5",     1,1200,1.5,1},  // hi 係数側 (Tmid=1000 の seam は FD 基準が不連続なので避ける)
        {"TP  T2500 M1.5",     1,2500,1.5,1},  // 高温 NASA 区間 (hi 係数)
        {"2PH T210 M1.6 genN", 2,210,1.6,1},   // 二相 (g=0.009 固定) Wyslouzil 出口相当
        {"2PH T230 M0.8",      2,230,0.8,1},
        {"2PH T260 M2.5",      2,260,2.5,1},
    };
    printf("%-22s %12s %14s %14s %10s\n","case","||LR-I||","||RLL-Afd||rel","||clsd-num||","minpiv");
    bool ok=true;
    for(auto&c:cs){
        // build Q
        double ro=0.4;
        double T=c.T;
        double R,cc,e;
        if(c.mode==0){ const double gam=1.4,cp=1039.0; R=cp*(gam-1)/gam; cc=sqrt(gam*R*T); e=(cp-R)*T; }
        else if(c.mode==1){ ThermoDerivatives D; thermo_derivatives_mix(&N2,1,Y1,T,&D); R=D.R; cc=sqrt(D.a2); e=D.e; }
        else { double cpm,hm; thermo_cph_mix(SP2,2,Y2,T,&cpm,&hm); const double Rmix=thermo_R_mix(SP2,2,Y2);
               R=Rmix-G_FIX*H2O.R; e=(hm-Rmix*T)+G_FIX*(H2O.R*T-cond_latent(H2O,T));
               const double dL=(cond_latent(H2O,T+0.1)-cond_latent(H2O,T-0.1))/0.2; double g2,c2; cond_twophase_sonic(cpm,R,G_FIX,dL,T,&g2,&c2); cc=sqrt(c2); }
        double n[3]; if(c.dir==0){n[0]=1;n[1]=0;n[2]=0;} else {n[0]=0.6;n[1]=0.7;n[2]=0.39; double l=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]); n[0]/=l;n[1]/=l;n[2]/=l;}
        double sp=c.M*cc; // speed along x
        double Q[5]={ro, ro*sp, ro*0.05*cc, 0.0, ro*(e+0.5*((sp)*(sp)+(0.05*cc)*(0.05*cc)))};
        Prim p=prim_from_Q(Q,c.mode);
        double Aeig[5][5]; double lrerr=build_A_and_LRerr(p,n,Aeig);
        double Afd[5][5]; A_fd(Q,n,c.mode,Afd);
        double rel=frodiff(Aeig,Afd)/fro(Afd);
        // numerical split A± (verified vs FD)
        double Ap[5][5],Am[5][5]; double mp=eos_split_jacobian_general(p.ux,p.uy,p.uz,n[0],n[1],n[2],p.c,p.Ht,p.kappa,Ap,Am);
        // closed-form split A± (Method B, no LU)
        double Cp[5][5],Cm[5][5];
        eos_split_jacobian_general_closed(p.ux,p.uy,p.uz,n[0],n[1],n[2],p.c,p.Ht,p.kappa,p.chi,Cp,Cm);
        // closed A+ vs numerical A+, closed A- vs numerical A-
        double clsdErr = (frodiff(Cp,Ap)+frodiff(Cm,Am))/fmax(fro(Ap)+fro(Am),1e-30);
        bool pass = (rel<1e-5 && lrerr<1e-9 && clsdErr<1e-9);
        printf("%-22s %12.2e %14.2e %14.2e %10.2e %s\n",c.name,lrerr,rel,clsdErr,mp,pass?"":" <-- FAIL");
        if(!pass) ok=false;
        // Level2: A+ + A- = A,  flip (closed form)
        double Asum[5][5]; for(int i=0;i<5;i++)for(int j=0;j<5;j++)Asum[i][j]=Cp[i][j]+Cm[i][j];
        double splitErr=frodiff(Asum,Aeig)/fro(Aeig);
        double nn[3]={-n[0],-n[1],-n[2]}; double Bp[5][5],Bm[5][5];
        eos_split_jacobian_general_closed(p.ux,p.uy,p.uz,nn[0],nn[1],nn[2],p.c,p.Ht,p.kappa,p.chi,Bp,Bm);
        // A_{-n}^+ = -A_n^-
        double flipErr=0; { double M1[5][5]; for(int i=0;i<5;i++)for(int j=0;j<5;j++)M1[i][j]=Bp[i][j]+Am[i][j]; flipErr=fro(M1)/fmax(fro(Aeig),1e-30); }
        if(splitErr>1e-9||flipErr>1e-9){ printf("    Level2: split=%.1e flip=%.1e <-- FAIL\n",splitErr,flipErr); ok=false; }
    }
    // ===== Level3: 実カーネル関数 accumulate_split_jacobian_cf の組み立てが解析 A± と一致するか =====
    // diag += face_area·A⁺,  nbr += (−A⁻)·sdq を複数面で累積し、検証済み eos_split_jacobian_general_closed と照合。
    {
        double ro=0.4, T=350.0; ThermoDerivatives D; thermo_derivatives_mix(&N2,1,Y1,T,&D);
        double e=D.e, ux=480,uy=40,uz=10;
        double Q[5]={ro, ro*ux, ro*uy, ro*uz, ro*(e+0.5*(ux*ux+uy*uy+uz*uz))};
        Prim p=prim_from_Q(Q,1);
        const double gamma=p.kappa+1.0;   // accumulate は kappa=gamma-1 を再構成
        double faces[2][4]={{0.6,0.7,0.39,1.3},{-0.2,0.5,-0.84,0.8}}; // nx,ny,nz(非正規),area
        double sdqs[2][5]={{0.01,-0.02,0.03,0.0,5.0},{0.04,0.01,-0.01,0.02,-3.0}};
        double diag[5][5]={{0}}, nbr[5]={0}, ediag[5][5]={{0}}, enbr[5]={0};
        for(int f=0;f<2;f++){
            double n[3]={faces[f][0],faces[f][1],faces[f][2]}; double L=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]); n[0]/=L;n[1]/=L;n[2]/=L;
            double fa=faces[f][3]; double* sdq=sdqs[f];
            block_dplur::accumulate_split_jacobian_cf<double>(gamma,n[0],n[1],n[2],p.ux,p.uy,p.uz,p.c,p.Ht,true,fa,true,sdq,diag,nbr);
            double Cp[5][5],Cm[5][5]; eos_split_jacobian_general_closed(p.ux,p.uy,p.uz,n[0],n[1],n[2],p.c,p.Ht,p.kappa,p.chi,Cp,Cm);
            for(int i=0;i<5;i++)for(int j=0;j<5;j++) ediag[i][j]+=fa*Cp[i][j];       // A⁺·face_area
            for(int i=0;i<5;i++){double s=0;for(int j=0;j<5;j++)s+=-Cm[i][j]*sdq[j]; enbr[i]+=s;} // (−A⁻)·sdq
        }
        double derr=frodiff(diag,ediag)/fmax(fro(ediag),1e-30);
        double nr=0,nd=0; for(int i=0;i<5;i++){double d=nbr[i]-enbr[i];nr+=d*d;nd+=enbr[i]*enbr[i];}
        double nerr=sqrt(nr)/fmax(sqrt(nd),1e-30);
        bool l3=(derr<1e-12 && nerr<1e-12);
        printf("Level3 (assembly): ||diag-face*A+||=%.2e  ||nbr-(-A-)sdq||=%.2e  %s\n", derr, nerr, l3?"PASS":"FAIL");
        if(!l3) ok=false;
    }
    // ===== Level3b: 二相 (固定 g,Y) 状態で実カーネル関数 (double / float) を照合 =====
    //   double: 閉形式 A± と 1e-12。float: 同じ入力を float に落として組み立て、diag·dq (方向微分) を double と比較 (相対 1e-5)。
    {
        double ro=0.15, T=215.0; double cpm,hm; thermo_cph_mix(SP2,2,Y2,T,&cpm,&hm); const double Rmix=thermo_R_mix(SP2,2,Y2);
        double e=(hm-Rmix*T)+G_FIX*(H2O.R*T-cond_latent(H2O,T)), ux=470,uy=30,uz=5;
        double Q[5]={ro, ro*ux, ro*uy, ro*uz, ro*(e+0.5*(ux*ux+uy*uy+uz*uz))};
        Prim p=prim_from_Q(Q,2);
        const double gamma=p.kappa+1.0;
        double faces[2][4]={{0.6,0.7,0.39,1.3},{-0.2,0.5,-0.84,0.8}};
        double sdqs[2][5]={{0.01,-0.02,0.03,0.0,5.0},{0.04,0.01,-0.01,0.02,-3.0}};
        double diag[5][5]={{0}}, nbr[5]={0}, ediag[5][5]={{0}}, enbr[5]={0};
        float fdiag[5][5]={{0}}, fnbr[5]={0};
        for(int f=0;f<2;f++){
            double n[3]={faces[f][0],faces[f][1],faces[f][2]}; double L=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]); n[0]/=L;n[1]/=L;n[2]/=L;
            double fa=faces[f][3]; double* sdq=sdqs[f];
            block_dplur::accumulate_split_jacobian_cf<double>(gamma,n[0],n[1],n[2],p.ux,p.uy,p.uz,p.c,p.Ht,true,fa,true,sdq,diag,nbr);
            float fsdq[5]; for(int i=0;i<5;i++) fsdq[i]=(float)sdq[i];
            block_dplur::accumulate_split_jacobian_cf<float>((float)gamma,(float)n[0],(float)n[1],(float)n[2],(float)p.ux,(float)p.uy,(float)p.uz,(float)p.c,(float)p.Ht,true,(float)fa,true,fsdq,fdiag,fnbr);
            double Cp[5][5],Cm[5][5]; eos_split_jacobian_general_closed(p.ux,p.uy,p.uz,n[0],n[1],n[2],p.c,p.Ht,p.kappa,p.chi,Cp,Cm);
            for(int i=0;i<5;i++)for(int j=0;j<5;j++) ediag[i][j]+=fa*Cp[i][j];
            for(int i=0;i<5;i++){double s=0;for(int j=0;j<5;j++)s+=-Cm[i][j]*sdq[j]; enbr[i]+=s;}
        }
        double derr=frodiff(diag,ediag)/fmax(fro(ediag),1e-30);
        double nr=0,nd=0; for(int i=0;i<5;i++){double d=nbr[i]-enbr[i];nr+=d*d;nd+=enbr[i]*enbr[i];}
        double nerr=sqrt(nr)/fmax(sqrt(nd),1e-30);
        // float 方向微分: diag·dq
        double dq[5]={0.02,-0.01,0.015,0.004,8.0}; double vd[5]={0},vf[5]={0}; double er=0,en=0;
        for(int i=0;i<5;i++){ for(int j=0;j<5;j++){ vd[i]+=diag[i][j]*dq[j]; vf[i]+=(double)fdiag[i][j]*dq[j]; } er+=(vd[i]-vf[i])*(vd[i]-vf[i]); en+=vd[i]*vd[i]; }
        double ferr=sqrt(er)/fmax(sqrt(en),1e-30);
        double fn=0,fnn=0; for(int i=0;i<5;i++){ fn+=((double)fnbr[i]-nbr[i])*((double)fnbr[i]-nbr[i]); fnn+=nbr[i]*nbr[i]; }
        double fnerr=sqrt(fn)/fmax(sqrt(fnn),1e-30);
        bool l3b=(derr<1e-12 && nerr<1e-12 && ferr<1e-5 && fnerr<1e-5);
        printf("Level3b (two-phase g=%.4f, gamma_2ph=%.5f c=%.3f): ||diag-face*A+||=%.2e ||nbr-(-A-)sdq||=%.2e float diag.dq rel=%.2e float nbr rel=%.2e  %s\n",
               G_FIX, gamma, p.c, derr, nerr, ferr, fnerr, l3b?"PASS":"FAIL");
        if(!l3b) ok=false;
    }
    printf("\nLevel1/2/3 %s\n", ok?"PASS":"FAIL");
    return ok?0:1;
}
