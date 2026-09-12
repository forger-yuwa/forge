#pragma once
// block-DPLUR の閉形式 FVS 蓄積 (rank-2 外積で A⁺S を diag・−A⁻·sdq を nbr へ)。
// host/device 双方でビルド可 (Level3 単体テスト tools/test_eos_jacobian.cpp から呼ぶため共有ヘッダ化)。
// 数式・検証は methods/time_integration + plan time_integration-general-eos-jacobian.md。
#include <cmath>
#ifdef FORGE_JACOBIAN_SANITY
#include <cstdio>
#endif

#if defined(__CUDACC__)
  #define BLKDPLUR_HD __host__ __device__ inline
#else
  #define BLKDPLUR_HD inline
#endif

namespace block_dplur {

// 非 NaN 値で CUDA の max/fmax とビット一致する三元 max (host/device 共通、double 昇格なし)。
template<typename T>
BLKDPLUR_HD T blkdplur_max(T a, T b) { return a > b ? a : b; }

// 閉形式 FVS 蓄積。M(g)=g₂I+(g₁−g₂)r₁⊗l₁+(g₅−g₂)r₅⊗l₅, a_plus=M(λ⁺), k_off=M((−λ)⁺)。
// diag += face_area·A⁺,  nbr += (−A⁻)·sdq  (sdq は呼び出し側で face_area·ΔQ_nbr 済)。
// thermallyPerfect=false: CPG 閉形式 (従来式・ビット不変)。
// thermallyPerfect=true:  thermally-perfect ideal gas (固定組成)。音響右固有ベクトルのエネルギーに実全
//   エンタルピー Htot を、左固有ベクトルの密度成分に EOS の χ_eos=c²−κh (κ=γ−1, h=Htot−ek) を加える。
//   CPG (χ_eos=0, Htot=ek+c²/κ) に厳密に還元。**注意**: κ=γ−1 を仮定 (TP ideal gas 専用。多成分 TP は
//   κ=γ_mix−1 で成立)、任意の実在気体 EOS では κ を熱力学ルーチンから受け取る必要がある。
// 命名注意: ローカル `kappa_over_c=(γ−1)/c` は **EOS の χ ではない**。EOS の χ_eos=c²−κh は `chi_eos`。
// DoDiag=false: diag への A⁺S 蓄積を省く (block-DPLUR 対角キャッシュの sweep≥1 用。nbr 側のみ計算)。
template<typename T, bool DoDiag = true>
BLKDPLUR_HD void accumulate_split_jacobian_cf(
    T gamma, T nx, T ny, T nz, T u, T v, T w, T c, T Htot, bool thermallyPerfect,
    T face_area, bool has_nbr, const T sdq[5], T diag[5][5], T nbr[5])
{
    const T sonic = blkdplur_max(c, static_cast<T>(1.0e-8));  // 本番下限 (異常検出は FORGE_JACOBIAN_SANITY)
    const T ek = static_cast<T>(0.5)*(u*u+v*v+w*w);
    const T V = u*nx+v*ny+w*nz;
    const T kappa        = gamma - static_cast<T>(1.0);  // Grüneisen κ = γ−1 (TP)
    const T kappa_over_c = kappa/sonic;                  // κ/c (EOS の χ ではない)
    const T c_over_kappa = sonic/kappa;                  // c/κ
    const T inv_sqrt2    = static_cast<T>(0.7071067811865475244);  // 1/√2 (音響固有ベクトルの双直交正規化)
    const T inv_sonic    = static_cast<T>(1.0)/sonic;    // 1/c
#ifdef FORGE_JACOBIAN_SANITY
    if (!(c > static_cast<T>(0.0)) || !isfinite(c) || kappa < static_cast<T>(1.0e-3) || !isfinite(Htot)) {
        printf("[jacobian_sanity] c=%g gamma=%g kappa=%g Htot=%g\n",
               (double)c, (double)gamma, (double)kappa, (double)Htot);
    }
#endif
    // eE = r 音響エネルギー成分 (∓V を除く), dD = l 音響密度成分の TP 補正 (χ_eos·(1/c))。
    // CPG: eE=ek/c + c/κ (=Htot_cpg/c), dD=0 → 従来式とビット一致。
    T eE, dD;
    if (thermallyPerfect) {
        const T h_real  = Htot - ek;                       // 実エンタルピー h(T)
        const T chi_eos = sonic*sonic - kappa*h_real;       // χ_eos = c²−κh (EOS の真の χ。CPG で 0)
        eE = Htot * inv_sonic;                              // = (ek+h)/c
        dD = chi_eos * inv_sonic;                           // = χ_eos/c
    } else {
        eE = ek*inv_sonic + c_over_kappa;
        dD = static_cast<T>(0.0);
    }
    const T r1[5]={ inv_sqrt2*inv_sonic, (u*inv_sonic+nx)*inv_sqrt2, (v*inv_sonic+ny)*inv_sqrt2, (w*inv_sonic+nz)*inv_sqrt2, (eE+V)*inv_sqrt2 };
    const T r5[5]={ inv_sqrt2*inv_sonic, (u*inv_sonic-nx)*inv_sqrt2, (v*inv_sonic-ny)*inv_sqrt2, (w*inv_sonic-nz)*inv_sqrt2, (eE-V)*inv_sqrt2 };
    const T l1[5]={ (kappa_over_c*ek-V+dD)*inv_sqrt2, (-kappa_over_c*u+nx)*inv_sqrt2, (-kappa_over_c*v+ny)*inv_sqrt2, (-kappa_over_c*w+nz)*inv_sqrt2, kappa_over_c*inv_sqrt2 };
    const T l5[5]={ (kappa_over_c*ek+V+dD)*inv_sqrt2, (-kappa_over_c*u-nx)*inv_sqrt2, (-kappa_over_c*v-ny)*inv_sqrt2, (-kappa_over_c*w-nz)*inv_sqrt2, kappa_over_c*inv_sqrt2 };
    const T lam1=V+sonic, lam2=V, lam5=V-sonic;
    const T zero=static_cast<T>(0.0);
    const T p2=blkdplur_max(lam2,zero), pa1=blkdplur_max(lam1,zero)-p2, pa5=blkdplur_max(lam5,zero)-p2;
    if (DoDiag) {
    #pragma unroll
    for (int i=0;i<5;++i){
        const T c1=pa1*r1[i], c5=pa5*r5[i];
        #pragma unroll
        for (int j=0;j<5;++j){
            T m = c1*l1[j] + c5*l5[j];
            if (i==j) m += p2;
            diag[i][j] += face_area * m;
        }
    }
    }
    if (has_nbr){
        const T n2=blkdplur_max(-lam2,zero), na1=blkdplur_max(-lam1,zero)-n2, na5=blkdplur_max(-lam5,zero)-n2;
        T l1s=zero, l5s=zero;
        #pragma unroll
        for (int j=0;j<5;++j){ l1s+=l1[j]*sdq[j]; l5s+=l5[j]*sdq[j]; }
        const T cc1=na1*l1s, cc5=na5*l5s;
        #pragma unroll
        for (int i=0;i<5;++i) nbr[i] += n2*sdq[i] + cc1*r1[i] + cc5*r5[i];
    }
}

}  // namespace block_dplur
