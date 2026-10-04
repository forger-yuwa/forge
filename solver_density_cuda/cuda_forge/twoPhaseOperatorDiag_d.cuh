#pragma once
// =============================================================================
// twoPhaseOperatorDiag_d.cuh — 診断 G3-a (作用素の収支; FORGE_DIAG_TP_OPERATOR=<出力 h5>、既定 off) の記録の配置
//   plans/active/condensation-two-phase-default.md §5.1 #4g3 (G3-a)・#4pjg、設計の判断
//   notes/reviews/2026-10-04-twophase-g3-design-diagnose.md・notes/reviews/2026-10-04-twophase-g1-g3-diag-design-diagnose.md。
//
// 本番カーネルの中で、診断ポインタ (tpoFace() / tpoSrcSlots(); 既定は無効) が有効なときだけ、そのカーネルが**実際に残差へ足した値**
// (atomic の前の面流束・格納した S·V) を面・節点の専用スロットへ写す (積算は後処理で double)。本番の演算式・評価順は変えない。
//   成分 c: 0 = ρY_w (総水分), 1 = ρg (液), 2 = ρQ2, 3 = ρQ1, 4 = ρQ0 (twoPhaseUpdateDiag_d.cuh と同じ番号)、
//           5.. = 水以外の化学種 (化学種番号の昇順、水を飛ばす)。ρv = ρY_w − ρg は後処理で double で作る。
//
// 面の記録 (面 = normal_halo_planes の並び ih; 内部面と node 境界半割面を含む):
//   val[(kind·nComp + c)·nF + ih] = 残差 res[ic0] に足した値 a0 (res[ic1] には −a0 を足す; 面積分値)。kind 0 = 移流、1 = 拡散。
//     移流 (species_advection_faceY_d): a0 = −ṁ·Y_f (流出が正の ṁ のとき ic0 から出る)。node 境界半割面 (ic1 = ghost) も足す (ic0 だけ)。
//     拡散 (species_diffusion_d / twophase_diffusion_d): a0 = J (セル ic0 へ入る向きが正)。node 境界半割面は足さない (値 0・code 2)。
//   code[kernel·nF + ih]: −1 = その kernel は経路に無い (未実行)、1 = 内部面を評価、2 = node 境界半割面で skip (拡散)、
//                         3 = node 境界半割面を評価 (移流; ic0 だけに足す)。
// 節点の記録 (凝縮ソース; slot·n + ic、n = nCells): 実際に足した項 (成分 g, Q2, Q1, Q0)、分岐、成長分岐の Sg のクリップ、
//   クリップ前の Sg、θ、1 項あたりの丸め回数 (float 実体は積と和で 2、double 実体は (float)(S·V) を足すので 1、足さなければ 0)、核生成率 J。
// main.cpp (g++) からも include する (定数・名前表・POD だけ; __device__ 関数は __CUDACC__ のときだけ)。
// =============================================================================
#include "flowFormat.hpp"

#define TPO_NC_CORE 5   // 成分 (w, g, Q2, Q1, Q0)。水以外の化学種はこの後ろ

// 面の kernel (code の行)
enum TpoKernel {
    TPO_K_ADV_SP = 0,   // 化学種の移流 (species_advection_faceY_d; S3 面組成)
    TPO_K_ADV_PA,       // 凝縮モーメントの移流 (受動種経路の species_advection_faceY_d; S3 面値)
    TPO_K_DIFF_OFF,     // 化学種拡散 (species_diffusion_d; 二相拡散 OFF)
    TPO_K_DIFF_ON,      // 二相拡散 (twophase_diffusion_d; 化学種・液・Q)
    TPO_NK
};
#define TPO_KIND_ADV  0
#define TPO_KIND_DIFF 1

// 凝縮ソースの節点スロット
enum TpoSrcSlot {
    TPO_S_TERM_G = 0,   // res_rog に足した項 (double; float 実体は double(S)·double(V) = 積の厳密値)
    TPO_S_TERM_Q2,
    TPO_S_TERM_Q1,
    TPO_S_TERM_Q0,
    TPO_S_BRANCH,       // 0 = 足さない (早期退出: ρ≈0・乾燥・表範囲外で乾燥)、1 = 核生成・成長、2 = 蒸発、
                        // 10/11/12 = double 実体 (表範囲外の委譲) の 足さない/核生成・成長/蒸発
    TPO_S_CLIP,         // 成長分岐の if (Sg < 0) Sg = 0 が作動したら 1 (それ以外の分岐は 0)
    TPO_S_SG_RAW,       // クリップ前の Sg (成長分岐のみ; それ以外は NaN)
    TPO_S_THETA,        // ソースに掛けた θ (成長分岐; 蒸発分岐は 1、足さない節点は NaN)
    TPO_S_NROUND,       // 1 項あたりの丸め回数 (0 / 1 / 2)
    TPO_S_J,            // 核生成率 J (成長分岐のみ)
    TPO_NSRC
};

// 面の記録の device ポインタ (値渡し; val == nullptr で無効)。compBase ≥ 0: 成分 = compBase + s (モーメント)、
// compBase < 0: 化学種 s → 成分 (s == iw ? 0 : 5 + 水を飛ばした番号)。
struct TpoFacePtr {
    double* val = nullptr;
    int* code = nullptr;
    geom_int nF = 0;
    int nComp = 0;
    int kind = 0;       // TPO_KIND_ADV / TPO_KIND_DIFF
    int kernel = -1;    // TpoKernel
    int compBase = 0;
    int iw = -1;
};

static inline const char* tpoKernelName(int k)
{
    static const char* nm[TPO_NK] = {"adv_species", "adv_moments", "diff_species_off", "diff_twophase_on"};
    return (k >= 0 && k < TPO_NK) ? nm[k] : "?";
}
static inline const char* tpoSrcName(int s)
{
    static const char* nm[TPO_NSRC] = {"term_g", "term_Q2", "term_Q1", "term_Q0", "branch", "clip", "Sg_unclipped", "theta", "nround", "J"};
    return (s >= 0 && s < TPO_NSRC) ? nm[s] : "?";
}

#ifdef __CUDACC__
__device__ inline int tpo_comp(const TpoFacePtr& t, int s)
{
    if (t.compBase >= 0) return t.compBase + s;
    return (s == t.iw) ? 0 : TPO_NC_CORE + ((s < t.iw) ? s : s - 1);
}
__device__ inline void tpo_face_put(const TpoFacePtr& t, geom_int ih, int comp, double a0)
{
    t.val[((size_t)t.kind*(size_t)t.nComp + (size_t)comp)*(size_t)t.nF + (size_t)ih] = a0;
}
__device__ inline void tpo_face_code(const TpoFacePtr& t, geom_int ih, int code)
{
    t.code[(size_t)t.kernel*(size_t)t.nF + (size_t)ih] = code;
}
__device__ inline void tpo_src_put(double* s, geom_int n, geom_int ic, int slot, double v) { s[(size_t)slot*(size_t)n + (size_t)ic] = v; }
// 「足さない」記録 (早期退出も明示的に 0 を書く)
__device__ inline void tpo_src_none(double* s, geom_int n, geom_int ic, int branch)
{
    for (int k = TPO_S_TERM_G; k <= TPO_S_TERM_Q0; ++k) tpo_src_put(s, n, ic, k, 0.0);
    tpo_src_put(s, n, ic, TPO_S_BRANCH, (double)branch);
    tpo_src_put(s, n, ic, TPO_S_CLIP, 0.0);
    const double qnan = __longlong_as_double(0x7ff8000000000000LL);
    tpo_src_put(s, n, ic, TPO_S_SG_RAW, qnan);
    tpo_src_put(s, n, ic, TPO_S_THETA, qnan);
    tpo_src_put(s, n, ic, TPO_S_NROUND, 0.0);
    tpo_src_put(s, n, ic, TPO_S_J, qnan);
}
// 足した 4 項 (g, Q2, Q1, Q0) と分岐・丸め回数
__device__ inline void tpo_src_terms(double* s, geom_int n, geom_int ic, int branch, double tg, double tq2, double tq1, double tq0, int nround)
{
    tpo_src_put(s, n, ic, TPO_S_TERM_G, tg);
    tpo_src_put(s, n, ic, TPO_S_TERM_Q2, tq2);
    tpo_src_put(s, n, ic, TPO_S_TERM_Q1, tq1);
    tpo_src_put(s, n, ic, TPO_S_TERM_Q0, tq0);
    tpo_src_put(s, n, ic, TPO_S_BRANCH, (double)branch);
    tpo_src_put(s, n, ic, TPO_S_NROUND, (double)nround);
}
#endif
