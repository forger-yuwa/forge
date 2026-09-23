#include "solidFem2d.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <highfive/highfive.hpp>

namespace conjugate {

namespace {

template <typename T>
std::vector<T> readVec(const HighFive::File& f, const std::string& name)
{
    std::vector<T> v;
    f.getDataSet(name).read(v);
    return v;
}

// (n,2) / (n,3) の 2 次元データセットを平坦に読む。
template <typename T>
std::vector<std::vector<T>> readMat(const HighFive::File& f, const std::string& name)
{
    std::vector<std::vector<T>> v;
    f.getDataSet(name).read(v);
    return v;
}

} // namespace

SolidMesh SolidMesh::read(const std::string& path)
{
    SolidMesh m;
    try {
        HighFive::File f(path, HighFive::File::ReadOnly);

        const auto coord = readMat<double>(f, "MESH/COORD");
        m.nNodes = (int)coord.size();
        m.x.resize(m.nNodes);
        m.y.resize(m.nNodes);
        for (int i = 0; i < m.nNodes; i++) { m.x[i] = coord[i][0]; m.y[i] = coord[i][1]; }

        const auto tris = readMat<int>(f, "MESH/TRIS");
        m.nTris = (int)tris.size();
        m.tris.resize(3 * m.nTris);
        for (int e = 0; e < m.nTris; e++)
            for (int k = 0; k < 3; k++) m.tris[3*e + k] = tris[e][k];

        m.ifaceNodes = readVec<int>(f, "IFACE/NODES");
        const auto ic = readMat<double>(f, "IFACE/COORD");
        m.ifaceX.resize(ic.size());
        m.ifaceY.resize(ic.size());
        for (size_t i = 0; i < ic.size(); i++) { m.ifaceX[i] = ic[i][0]; m.ifaceY[i] = ic[i][1]; }

        const auto ie = readMat<int>(f, "IFACE/EDGES");
        m.ifaceEdges.resize(2 * ie.size());
        for (size_t i = 0; i < ie.size(); i++) { m.ifaceEdges[2*i] = ie[i][0]; m.ifaceEdges[2*i+1] = ie[i][1]; }

        const auto re = readMat<int>(f, "ROBIN/EDGES");
        m.robinEdges.resize(2 * re.size());
        for (size_t i = 0; i < re.size(); i++) { m.robinEdges[2*i] = re[i][0]; m.robinEdges[2*i+1] = re[i][1]; }
        m.robinH  = readVec<double>(f, "ROBIN/H");
        m.robinTc = readVec<double>(f, "ROBIN/TC");

        m.kT = readVec<double>(f, "SOLID/K_T");
        m.kV = readVec<double>(f, "SOLID/K_V");

        if (f.hasAttribute("bandwidth")) f.getAttribute("bandwidth").read(m.bandwidth);
        if (f.hasAttribute("iface_sha1")) f.getAttribute("iface_sha1").read(m.ifaceSha1);
        if (f.hasAttribute("source_npz")) f.getAttribute("source_npz").read(m.sourceNpz);
    } catch (const std::exception& e) {
        std::cerr << "[solidFem2d] ERROR: 固体メッシュ " << path << " を読めない: " << e.what() << "\n"
                  << "[solidFem2d]   tools/solid_mesh_to_h5.py で作ること "
                     "(MESH/COORD, MESH/TRIS, IFACE/*, ROBIN/*, SOLID/K_*)。\n";
        exit(EXIT_FAILURE);
    }

    if (m.nNodes <= 0 || m.nTris <= 0 || m.ifaceNodes.empty()) {
        std::cerr << "[solidFem2d] ERROR: 固体メッシュが空 (nodes=" << m.nNodes
                  << " tris=" << m.nTris << " iface=" << m.ifaceNodes.size() << ")\n";
        exit(EXIT_FAILURE);
    }
    if (m.robinEdges.empty()) {
        // 全断熱は定数零空間を持ち、正味入熱が非零なら定常解が無い (plan §4.4a)。
        std::cerr << "[solidFem2d] ERROR: Robin 辺が 1 本も無い (背面・孔がすべて断熱)。"
                     "定数零空間を持つので定常解が定まらない。\n";
        exit(EXIT_FAILURE);
    }
    return m;
}

double SolidMesh::kOf(double T) const
{
    if (kT.size() <= 1) return kV.empty() ? 0.0 : kV[0];
    if (T <= kT.front()) return kV.front();
    if (T >= kT.back())  return kV.back();
    const auto it = std::upper_bound(kT.begin(), kT.end(), T);
    const size_t j = (size_t)(it - kT.begin());
    const double t0 = kT[j-1], t1 = kT[j];
    const double w = (t1 > t0) ? (T - t0) / (t1 - t0) : 0.0;
    return kV[j-1] + w * (kV[j] - kV[j-1]);
}

std::vector<double> SolidMesh::ifaceLumped() const
{
    std::vector<int> local(nNodes, -1);
    for (size_t i = 0; i < ifaceNodes.size(); i++) local[ifaceNodes[i]] = (int)i;
    std::vector<double> a(ifaceNodes.size(), 0.0);
    for (size_t e = 0; e + 1 < ifaceEdges.size(); e += 2) {
        const int n0 = ifaceEdges[e], n1 = ifaceEdges[e+1];
        const double L = std::hypot(x[n1] - x[n0], y[n1] - y[n0]);
        if (local[n0] >= 0) a[local[n0]] += 0.5 * L;
        if (local[n1] >= 0) a[local[n1]] += 0.5 * L;
    }
    return a;
}

SolidFem2D::SolidFem2D(const SolidMesh& m) : m_(m), n_(m.nNodes)
{
    // 帯幅は**実際の非零構造から測り直す** (属性は三角形の辺だけで測っているので、
    // Robin 辺が別の節点対を結んでいれば足りなくなる)。
    int bw = 0;
    auto upd = [&bw](int a, int b) { bw = std::max(bw, std::abs(a - b)); };
    for (int e = 0; e < m_.nTris; e++) {
        const int* t = &m_.tris[3*e];
        upd(t[0], t[1]); upd(t[1], t[2]); upd(t[2], t[0]);
    }
    for (size_t i = 0; i < m_.robinEdges.size(); i += 2) upd(m_.robinEdges[i], m_.robinEdges[i+1]);
    bw_ = bw;
    if (m_.bandwidth > 0 && bw_ > m_.bandwidth)
        std::cout << "[solidFem2d] 注意: 実測帯幅 " << bw_ << " が h5 の属性 " << m_.bandwidth
                  << " より大きい (Robin 辺が三角形の辺でない)。実測値を使う。\n";
    ab_.assign((size_t)(bw_ + 1) * n_, 0.0);
    abF_.assign((size_t)(bw_ + 1) * n_, 0.0);
    b_.assign(n_, 0.0);
}

void SolidFem2D::assemble(const std::vector<double>& u,
                          const std::vector<double>& dfDiag,
                          const std::vector<double>& dfRhs)
{
    std::fill(ab_.begin(), ab_.end(), 0.0);
    std::fill(b_.begin(), b_.end(), 0.0);
    factored_ = false;

    // ---- 伝導 (線形三角形。要素の k は**節点値の平均** = Python と同一規約) ----
    for (int e = 0; e < m_.nTris; e++) {
        const int* t = &m_.tris[3*e];
        const double x0 = m_.x[t[0]], x1 = m_.x[t[1]], x2 = m_.x[t[2]];
        const double y0 = m_.y[t[0]], y1 = m_.y[t[1]], y2 = m_.y[t[2]];
        const double bb[3] = {y1 - y2, y2 - y0, y0 - y1};
        const double cc[3] = {x2 - x1, x0 - x2, x1 - x0};
        const double det = (x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0);
        const double area = 0.5 * std::fabs(det);
        if (!(area > 0.0)) continue;
        const double kmean = (m_.kOf(u[t[0]]) + m_.kOf(u[t[1]]) + m_.kOf(u[t[2]])) / 3.0;
        const double f = kmean / (4.0 * area);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                if (t[i] < t[j]) continue;                    // 下三角のみ
                at(t[i], t[j]) += f * (bb[i]*bb[j] + cc[i]*cc[j]);
            }
    }

    // ---- Robin 辺 (consistent: M = hL/6 [[2,1],[1,2]], f += h Tc L/2) ----
    for (size_t r = 0; r < m_.robinH.size(); r++) {
        const int n0 = m_.robinEdges[2*r], n1 = m_.robinEdges[2*r+1];
        const double L = std::hypot(m_.x[n1] - m_.x[n0], m_.y[n1] - m_.y[n0]);
        if (!(L > 0.0)) continue;
        const double h = m_.robinH[r], Tc = m_.robinTc[r];
        const double m2 = h * L / 6.0;
        const int idx[2] = {n0, n1};
        const double Me[2][2] = {{2.0*m2, 1.0*m2}, {1.0*m2, 2.0*m2}};
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                if (idx[i] < idx[j]) continue;
                at(idx[i], idx[j]) += Me[i][j];
            }
            b_[idx[i]] += h * Tc * L / 2.0;
        }
    }

    // ---- 連成項 D_f (界面対角のみ。§4.6a) ----
    if (!dfDiag.empty()) {
        for (int i = 0; i < m_.nIface(); i++) at(m_.ifaceNodes[i], m_.ifaceNodes[i]) += dfDiag[i];
    }
    if (!dfRhs.empty()) {
        for (int i = 0; i < m_.nIface(); i++) b_[m_.ifaceNodes[i]] += dfRhs[i];
    }
}

std::vector<double> SolidFem2D::matvec(const std::vector<double>& u) const
{
    std::vector<double> y(n_, 0.0);
    for (int j = 0; j < n_; j++) {
        const int kmax = std::min(bw_, n_ - 1 - j);
        y[j] += at(j, j) * u[j];
        for (int k = 1; k <= kmax; k++) {
            const double a = at(j + k, j);
            if (a == 0.0) continue;
            y[j + k] += a * u[j];
            y[j]     += a * u[j + k];
        }
    }
    return y;
}

void SolidFem2D::factorize()
{
    abF_ = ab_;                    // 組んだ行列は残す (同じ更新で残差も測るため)
    // 下三角バンド Cholesky (LAPACK dpbtf2 と同じ順序)。
    for (int j = 0; j < n_; j++) {
        double d = atF(j, j);
        if (!(d > 0.0)) {
            std::cerr << "[solidFem2d] ERROR: 固体行列が正定でない (節点 " << j
                      << " で対角 " << d << ")。孔・背面の Robin が付いていない"
                         "(= 定数零空間) か、メッシュが壊れている。\n";
            exit(EXIT_FAILURE);
        }
        d = std::sqrt(d);
        atF(j, j) = d;
        const int kmax = std::min(bw_, n_ - 1 - j);
        for (int i = 1; i <= kmax; i++) atF(j + i, j) /= d;
        for (int k = 1; k <= kmax; k++) {
            const double ajk = atF(j + k, j);
            if (ajk == 0.0) continue;
            for (int i = k; i <= kmax; i++) atF(j + i, j + k) -= atF(j + i, j) * ajk;
        }
    }
    factored_ = true;
}

void SolidFem2D::solveInPlace(std::vector<double>& v) const
{
    // L y = v
    for (int j = 0; j < n_; j++) {
        v[j] /= atF(j, j);
        const int kmax = std::min(bw_, n_ - 1 - j);
        for (int k = 1; k <= kmax; k++) v[j + k] -= atF(j + k, j) * v[j];
    }
    // L^T x = y
    for (int j = n_ - 1; j >= 0; j--) {
        const int kmax = std::min(bw_, n_ - 1 - j);
        for (int k = 1; k <= kmax; k++) v[j] -= atF(j + k, j) * v[j + k];
        v[j] /= atF(j, j);
    }
}

int SolidFem2D::solveLoad(const std::vector<double>& Qf, std::vector<double>& u,
                          int maxIter, double tol)
{
    if ((int)u.size() != n_) u.assign(n_, 300.0);
    std::vector<double> rhs(n_);
    for (int it = 1; it <= maxIter; it++) {
        assemble(u);
        rhs = b_;
        for (int i = 0; i < m_.nIface(); i++) rhs[m_.ifaceNodes[i]] += Qf[i];
        factorize();
        solveInPlace(rhs);
        double dmax = 0.0, umax = 1.0;
        for (int i = 0; i < n_; i++) {
            dmax = std::max(dmax, std::fabs(rhs[i] - u[i]));
            umax = std::max(umax, std::fabs(rhs[i]));
        }
        u = rhs;
        if (dmax < tol * umax) return it;
    }
    return maxIter;
}

std::vector<double> SolidFem2D::residual(const std::vector<double>& u)
{
    assemble(u);                       // 現在の u の k_s(T) で組み直す
    std::vector<double> r = matvec(u);
    for (int i = 0; i < n_; i++) r[i] -= b_[i];
    return r;
}

void writeSolidField(const std::string& stem, const SolidMesh& m,
                     const std::vector<double>& u, const std::vector<double>& qIface,
                     double timeValue)
{
    const int N = m.nNodes;
    std::vector<std::vector<double>> coord(N, std::vector<double>(3, 0.0));
    for (int i = 0; i < N; i++) { coord[i][0] = m.x[i]; coord[i][1] = m.y[i]; }

    // XDMF Mixed: 三角形は要素コード 4 (output/output.cpp と同じ規約)
    std::vector<int> conne;
    conne.reserve(4 * m.nTris);
    for (int e = 0; e < m.nTris; e++) {
        conne.push_back(4);
        for (int k = 0; k < 3; k++) conne.push_back(m.tris[3*e + k]);
    }

    std::vector<double> ks(N), qh(N, 0.0), qi(N, 0.0);
    for (int i = 0; i < N; i++) ks[i] = m.kOf(u[i]);
    for (size_t r = 0; r < m.robinH.size(); r++) {
        const int n0 = m.robinEdges[2*r], n1 = m.robinEdges[2*r+1];
        const double L = std::hypot(m.x[n1] - m.x[n0], m.y[n1] - m.y[n0]);
        const double h = m.robinH[r], Tc = m.robinTc[r];
        qh[n0] += 0.5 * L * h * (u[n0] - Tc);
        qh[n1] += 0.5 * L * h * (u[n1] - Tc);
    }
    for (int i = 0; i < m.nIface(); i++)
        if (i < (int)qIface.size()) qi[m.ifaceNodes[i]] = qIface[i];

    const std::string h5name = stem + ".h5";
    {
        HighFive::File f(h5name, HighFive::File::Overwrite);
        f.createDataSet("MESH/COORD", coord);
        f.createDataSet("MESH/CONNE", conne);
        f.createDataSet("VALUE/T", u);
        f.createDataSet("VALUE/k_s", ks);
        f.createDataSet("VALUE/q_iface", qi);
        f.createDataSet("VALUE/q_hole", qh);
    }

    const size_t slash = h5name.find_last_of('/');
    const std::string h5base = (slash == std::string::npos) ? h5name : h5name.substr(slash + 1);
    std::ofstream ofs(stem + ".xmf");
    ofs << "<?xml version='1.0' ?>\n<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>\n<Xdmf>\n  <Domain>\n";
    ofs << "    <Grid GridType='Collection' CollectionType='Spatial' Name='Mixed'>\n";
    ofs << "    <Time TimeType='Single' Value='" << timeValue << "'/> \n";
    ofs << "      <Grid Name='solid'>\n";
    ofs << "        <Topology Type='Mixed' NumberOfElements='" << m.nTris << "'>\n";
    ofs << "          <DataItem Format='HDF' DataType='Int' Dimensions='" << conne.size() << "'>\n";
    ofs << "            " << h5base << ":MESH/CONNE\n          </DataItem>\n        </Topology>\n";
    ofs << "        <Geometry Type='XYZ'>\n";
    ofs << "          <DataItem Format='HDF' DataType='Float' Precision='8' Dimensions='" << N*3 << "'>\n";
    ofs << "            " << h5base << ":MESH/COORD\n          </DataItem>\n        </Geometry>\n";
    for (const char* nm : {"T", "k_s", "q_iface", "q_hole"}) {
        ofs << "        <Attribute Name='" << nm << "' Center='Node' >\n";
        ofs << "          <DataItem Format='HDF' DataType='Float' Precision='8' Dimensions='" << N << "'>\n";
        ofs << "            " << h5base << ":VALUE/" << nm << "\n          </DataItem>\n        </Attribute>\n";
    }
    ofs << "      </Grid>\n    </Grid>\n  </Domain>\n</Xdmf>\n";
}

} // namespace conjugate
