#pragma once

// ============================================================================
// CHT Phase 2: 一般 2D 固体 (`fem2d`) の伝導 — **ソルバ内連成の固体側**
//
// 仕様は plans/accepted/boundary-conjugate-heat-transfer.md §4.6a、
// 界面熱量の符号は conjugateWall.hpp、外部ループの参照実装は
// solver_density_cuda/tools/solid_fem2d.py (**こちらが真値**。移植の同値試験 §6 V4b(a))。
//
// 【離散化】線形三角形 (1 次要素)。要素の熱伝導率は**節点値の平均**で評価する
//   (Python の `assemble_full` と同一規約。ここを変えると a1 が 1e-12 で合わなくなる)。
//     K_e = k_mean/(4A) (b b^T + c c^T),   A = |det|/2
//   冷却孔は Robin 辺で、**consistent** な辺行列を使う:
//     M_e = hL/6 [[2,1],[1,2]],   f_e += h T_c L/2  (各節点)
//
// 【解き方】Schur 補元は作らず**全節点系を直接解く** (§4.6a)。
//   未知数は全節点温度 u なので内部温度が状態になり、k_s(T) の自己整合が
//   「復元してから組み直す」操作を要さない (plan §5.1 #33 と同型の事故が起きない)。
//   行列は SPD (Robin 辺がある限り) で、変換時に RCM 済みなので**バンド Cholesky**で解く。
//   帯幅は C3X で 82 / Mark II で 69 (N b^2 = 5.0e7 / 2.9e7 flop)。
//
// 【格納】下三角バンド。ab_[k + j*(bw+1)] = A(j+k, j)  (k = 0..bw, 行 >= 列)。
// ============================================================================

#include <string>
#include <vector>

namespace conjugate {

// 固体メッシュ (tools/solid_mesh_to_h5.py が書く HDF5 をそのまま持つ)。
struct SolidMesh {
    int nNodes = 0;
    int nTris = 0;
    int bandwidth = 0;                 // 変換時に測った帯幅 (RCM 後)

    std::vector<double> x, y;          // nNodes
    std::vector<int>    tris;          // 3*nTris
    std::vector<int>    ifaceNodes;    // nIface
    std::vector<double> ifaceX, ifaceY;// nIface (流体の壁節点と 1 対 1 に突き合わせる)
    std::vector<int>    ifaceEdges;    // 2*nIfaceEdges
    std::vector<int>    robinEdges;    // 2*nRobin
    std::vector<double> robinH, robinTc;
    std::vector<double> kT, kV;        // k_s(T) のテーブル (1 点なら定数)
    std::string ifaceSha1;
    std::string contentSha1;   // 中身全体 (節点順・接続・物性・冷却条件) のハッシュ
    std::string sourceNpz;

    static SolidMesh read(const std::string& path);   // 読めなければ exit(1)
    double kOf(double T) const;                       // 線形内挿 (両端はクランプ)
    // 界面節点の**集中辺長** [m] (IFACE/NODES の順)。連成荷重は q_eff [W/m2] × これ。
    // **流体側の surfArea を使ってはいけない** (押し出し疑似 2D で奥行きが乗る。§4.6a)。
    std::vector<double> ifaceLumped() const;
    int nIface() const { return (int)ifaceNodes.size(); }
};

class SolidFem2D {
public:
    explicit SolidFem2D(const SolidMesh& m);

    // K(u) と b を組む。dfDiag が非空なら**界面節点の対角**に加える (連成の D_f)。
    // dfRhs が非空なら界面節点の右辺に加える (D_f T^k の項)。
    void assemble(const std::vector<double>& u,
                  const std::vector<double>& dfDiag = {},
                  const std::vector<double>& dfRhs  = {});

    // バンド Cholesky。**組んだ行列 ab_ は壊さず、複製 abF_ を分解する**ので、
    // 分解を再利用しながら同じ更新で残差 (matvec) も測れる。
    void factorize();
    void solveInPlace(std::vector<double>& v) const;  // 前進後退代入 (abF_ を使う)
    std::vector<double> matvec(const std::vector<double>& u) const;  // 組んだ行列 (分解前) の K u

    // 界面に節点荷重 Qf [W/m] を与えて解く (D_f = 0、k_s(T) は Picard)。
    // 戻り値は反復回数。u は入出力 (初期値として使う)。
    int solveLoad(const std::vector<double>& Qf, std::vector<double>& u,
                  int maxIter = 30, double tol = 1e-12);

    // 現在の u で組み直した残差 K(u)u - b。界面成分が「固体が受け取った熱」[W/m]。
    std::vector<double> residual(const std::vector<double>& u);

    const std::vector<double>& rhs() const { return b_; }
    int bandwidth() const { return bw_; }
    bool factored() const { return factored_; }

private:
    const SolidMesh& m_;
    int n_ = 0;
    int bw_ = 0;
    std::vector<double> ab_;     // 下三角バンド (組んだまま)
    std::vector<double> abF_;    // その複製を分解したもの
    std::vector<double> b_;
    bool factored_ = false;

    inline double& at(int i, int j) { return ab_[(i - j) + (size_t)j * (bw_ + 1)]; }
    inline double  at(int i, int j) const { return ab_[(i - j) + (size_t)j * (bw_ + 1)]; }
    inline double& atF(int i, int j) { return abF_[(i - j) + (size_t)j * (bw_ + 1)]; }
    inline double  atF(int i, int j) const { return abF_[(i - j) + (size_t)j * (bw_ + 1)]; }
};

// 固体場を forge と同じ XDMF 規約で書く (plan §4.6a。ParaView で流体と重ねられる)。
//   VALUE/T (節点温度) / k_s / q_iface (ガス側から受け取った熱 [W/m]) / q_hole (孔が持ち去った熱)
// `stem` は拡張子なしのパス (`res_solid_1000` → `.h5` と `.xmf` を書く)。
void writeSolidField(const std::string& stem, const SolidMesh& m,
                     const std::vector<double>& u, const std::vector<double>& qIface,
                     double timeValue);

} // namespace conjugate
