// 閉性の照合の診断 (plans/active/axisymmetric-freestream-hoop-gauge.md §4.6 の 3・5)。説明と書く量は hoopClosureDiag.hpp。
// 既定 off: main が起動の段取りの後に 1 回呼び、環境変数 FORGE_DIAG_HOOP_CLOSURE があるときだけ書く (計算は続ける)。

#include "cuda_forge/hoopClosureDiag.hpp"
#include "cuda_forge/cudaWrapper.cuh"

#include <highfive/H5File.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

namespace hoopClosureDiag {
namespace {

template <typename T>
std::vector<T> d2h(const T* d, size_t n)
{
    std::vector<T> h(n);
    if (n > 0 && d != nullptr) gpuErrchk( cudaMemcpy(h.data(), d, n*sizeof(T), cudaMemcpyDeviceToHost) );
    return h;
}

template <typename T>
void put(HighFive::File& h5, const std::string& name, const std::vector<T>& v)
{
    if (v.empty()) {   // 長さ 0 は書き込みを呼ばずに空のデータセットだけ作る
        h5.createDataSet<T>(name, HighFive::DataSpace(std::vector<size_t>{0}));
        return;
    }
    h5.createDataSet(name, v);
}

}  // namespace

void runIfRequested(const solverConfig& cfg, const mesh& msh, variables& var)
{
    const char* path = std::getenv("FORGE_DIAG_HOOP_CLOSURE");
    if (path == nullptr || *path == '\0') return;

    const size_t nP = (size_t)msh.nPlanes, nC = (size_t)msh.nCells, nCa = (size_t)msh.nCells_all;
    const bool axisym = (cfg.isAxisymmetric == 1);
    const bool rWeighted = axisym && cfg.axisymMethod == 0;   // Σ±S_y の目標が A_planar になる経路
    std::string reason;
    const bool segActive = axisSegmentRWeightApplies(cfg, msh, &reason);

    // デバイスに渡した最終の値 (壁・軸の射影とは無関係の幾何)
    const std::vector<flow_float> sx = d2h(var.p_d.at("sx"), nP), sy = d2h(var.p_d.at("sy"), nP);
    const std::vector<flow_float> sz = d2h(var.p_d.at("sz"), nP), ss = d2h(var.p_d.at("ss"), nP);
    const std::vector<flow_float> vol = d2h(var.c_d.at("volume"), nCa);
    const std::vector<flow_float> Ap = axisym ? d2h(var.c_d.at("A_planar"), nCa) : std::vector<flow_float>{};
    const std::vector<geom_int> pc = d2h(msh.map_plane_cells_d, 2*nP);

    // 面: ss の有限性・正値、法線の長さ ‖(sx,sy,sz)‖/ss − 1 (double)
    std::vector<double> normErr(nP, 0.0);
    long long nSsNonFinite = 0, nSsNonPos = 0;
    double maxNormErr = 0.0;
    for (size_t ip = 0; ip < nP; ++ip) {
        const double a = (double)ss[ip];
        if (!std::isfinite(a)) { ++nSsNonFinite; normErr[ip] = std::numeric_limits<double>::quiet_NaN(); continue; }
        if (!(a > 0.0))        { ++nSsNonPos;    normErr[ip] = std::numeric_limits<double>::quiet_NaN(); continue; }
        const double x = sx[ip], y = sy[ip], z = sz[ip];
        normErr[ip] = std::sqrt(x*x + y*y + z*z) / a - 1.0;
        maxNormErr = std::max(maxNormErr, std::fabs(normErr[ip]));
    }

    // CV: Σ±S、Σ|S| を double で (面の向きは plane_cells の ic0 に +、ic1 に −。ゴーストは数えない)
    std::vector<double> sumx(nC, 0.0), sumy(nC, 0.0), absx(nC, 0.0), absy(nC, 0.0);
    for (size_t ip = 0; ip < nP; ++ip) {
        const double x = sx[ip], y = sy[ip];
        for (int k = 0; k < 2; ++k) {
            const geom_int ic = pc[2*ip + k];
            if (ic < 0 || (size_t)ic >= nC) continue;
            const double sg = (k == 0) ? 1.0 : -1.0;
            sumx[ic] += sg*x; sumy[ic] += sg*y;
            absx[ic] += std::fabs(x); absy[ic] += std::fabs(y);
        }
    }
    const double eps64 = std::numeric_limits<double>::epsilon();
    std::vector<double> Ex(nC, 0.0), Ey(nC, 0.0);
    double maxEx = 0.0, maxEy = 0.0;
    long long nOverX = 0, nOverY = 0, nBadA = 0;
    for (size_t ic = 0; ic < nC; ++ic) {
        const double A = axisym ? (double)Ap[ic] : (double)vol[ic];
        const double Ty = rWeighted ? (double)Ap[ic] : 0.0;
        const double dx = std::fabs(sumx[ic]), dy = std::fabs(sumy[ic] - Ty);
        if (!(A > 0.0) || !std::isfinite(A)) { ++nBadA; Ex[ic] = Ey[ic] = std::numeric_limits<double>::quiet_NaN(); continue; }
        Ex[ic] = dx / A; Ey[ic] = dy / A;
        maxEx = std::max(maxEx, Ex[ic]); maxEy = std::max(maxEy, Ey[ic]);
        if (dx > 100.0*eps64*(A + absx[ic])) ++nOverX;
        if (dy > 100.0*eps64*(A + absy[ic])) ++nOverY;
    }

    HighFive::File h5(path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    h5.createAttribute("plan", std::string("plans/active/axisymmetric-freestream-hoop-gauge.md §4.6 (3, 5)"));
    h5.createAttribute("flow_float_bytes", (long long)sizeof(flow_float));
    h5.createAttribute("geom_float_bytes", (long long)sizeof(geom_float));
    h5.createAttribute("nCells", (long long)msh.nCells);
    h5.createAttribute("nCells_all", (long long)msh.nCells_all);
    h5.createAttribute("nPlanes", (long long)msh.nPlanes);
    h5.createAttribute("nNormalPlanes", (long long)msh.nNormalPlanes);
    h5.createAttribute("discretization", cfg.discretization);
    h5.createAttribute("mesh_file", cfg.meshFileName);
    h5.createAttribute("isAxisymmetric", (long long)cfg.isAxisymmetric);
    h5.createAttribute("axisymMethod", (long long)cfg.axisymMethod);
    h5.createAttribute("axisRFloor", (double)cfg.axisRFloor);
    h5.createAttribute("hoopAreaFromClosure", (long long)cfg.hoopAreaFromClosure);
    h5.createAttribute("axisSegmentRWeight", (long long)cfg.axisSegmentRWeight);
    h5.createAttribute("mesh_has_rSurfVect", (long long)(msh.rSurfVect64.empty() ? 0 : 1));
    h5.createAttribute("segment_rweight_active", (long long)(segActive ? 1 : 0));
    h5.createAttribute("segment_rweight_reason", segActive ? std::string("active") : reason);
    h5.createAttribute("closure_target_y", std::string(rWeighted ? "A_planar" : "0"));
    h5.createAttribute("closure_norm", std::string(axisym ? "A_planar" : "volume"));
    h5.createAttribute("n_ss_nonfinite", nSsNonFinite);
    h5.createAttribute("n_ss_nonpositive", nSsNonPos);
    h5.createAttribute("max_abs_normlen_err", maxNormErr);
    h5.createAttribute("max_E_x", maxEx);
    h5.createAttribute("max_E_y", maxEy);
    h5.createAttribute("n_over_100eps64_x", nOverX);
    h5.createAttribute("n_over_100eps64_y", nOverY);
    h5.createAttribute("n_bad_A", nBadA);

    put(h5, "/faces/sx", sx); put(h5, "/faces/sy", sy); put(h5, "/faces/sz", sz); put(h5, "/faces/ss", ss);
    put(h5, "/faces/plane_cells", pc);
    put(h5, "/faces/normlen_err", normErr);
    if (axisym) put(h5, "/cells/A_planar", Ap);
    put(h5, "/cells/volume", vol);
    put(h5, "/cells/sum64_x", sumx); put(h5, "/cells/sum64_y", sumy);
    put(h5, "/cells/abs64_x", absx); put(h5, "/cells/abs64_y", absy);
    put(h5, "/cells/E_x", Ex); put(h5, "/cells/E_y", Ey);
    if (msh.cc64.size() >= 3*nC) put(h5, "/cells/cc64", std::vector<double>(msh.cc64.begin(), msh.cc64.begin() + 3*nC));

    printf("[hoopClosure] FORGE_DIAG_HOOP_CLOSURE=%s: segment r-weight %s; ss non-finite %lld, non-positive %lld, "
           "max|normlen-1| %.3e; max E_x %.3e, max E_y %.3e (target_y %s); CVs over 100*eps64*(A+sum|S|): x %lld y %lld "
           "(continuing the run)\n",
           path, segActive ? "ON" : ("OFF: " + reason).c_str(), nSsNonFinite, nSsNonPos, maxNormErr, maxEx, maxEy,
           rWeighted ? "A_planar" : "0", nOverX, nOverY);
    fflush(stdout);
}

}  // namespace hoopClosureDiag
