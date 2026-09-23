// CHT Phase 2: 固体 FE (`fem2d`) の**同値試験用 CLI** (plan §6 V4b (a))。
//
// ソルバ本体に組み込む前に、参照実装 (tools/solid_fem2d.py) と**同じ問題を同じ精度で**
// 解けることをここで閉じる。突き合わせは tools/test_solid_fem2d_cpp.py が行う。
//
//   solid_fem2d_tool info   <solid.h5>
//   solid_fem2d_tool matvec <solid.h5> <u.txt>  <out.txt>   … K(u) u を書く (a1)
//   solid_fem2d_tool solve  <solid.h5> <Qf.txt> <out.txt>   … 界面荷重で解いて u を書く (a2)
//
// テキストは 1 行 1 値 (倍精度)。`solve` の Qf は **IFACE/NODES の順**、単位は W/m。

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "solidFem2d.hpp"

namespace {

std::vector<double> readText(const std::string& path)
{
    std::ifstream ifs(path);
    if (!ifs) { std::cerr << "[solid_fem2d_tool] cannot open " << path << "\n"; exit(EXIT_FAILURE); }
    std::vector<double> v;
    double t;
    while (ifs >> t) v.push_back(t);
    return v;
}

void writeText(const std::string& path, const std::vector<double>& v)
{
    std::ofstream ofs(path);
    ofs << std::setprecision(17) << std::scientific;
    for (const double t : v) ofs << t << "\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "usage: solid_fem2d_tool info|matvec|solve <solid.h5> [in.txt out.txt]\n";
        return 1;
    }
    const std::string mode = argv[1];
    const conjugate::SolidMesh m = conjugate::SolidMesh::read(argv[2]);
    conjugate::SolidFem2D fem(m);

    if (mode == "info") {
        std::cout << "nodes " << m.nNodes << "  tris " << m.nTris
                  << "  iface " << m.nIface() << "  robin " << m.robinH.size()
                  << "  bandwidth " << fem.bandwidth() << " (attr " << m.bandwidth << ")\n"
                  << "k_s(300) " << m.kOf(300.0) << "  k_s(900) " << m.kOf(900.0)
                  << "  iface_sha1 " << m.ifaceSha1 << "\n";
        return 0;
    }
    if (argc < 5) { std::cerr << "usage: ... <in.txt> <out.txt>\n"; return 1; }

    if (mode == "matvec") {
        const std::vector<double> u = readText(argv[3]);
        if ((int)u.size() != m.nNodes) {
            std::cerr << "[solid_fem2d_tool] u の長さ " << u.size() << " != nodes " << m.nNodes << "\n";
            return 1;
        }
        fem.assemble(u);
        writeText(argv[4], fem.matvec(u));
        return 0;
    }
    if (mode == "solve") {
        const std::vector<double> Qf = readText(argv[3]);
        if ((int)Qf.size() != m.nIface()) {
            std::cerr << "[solid_fem2d_tool] Qf の長さ " << Qf.size() << " != iface " << m.nIface() << "\n";
            return 1;
        }
        std::vector<double> u(m.nNodes, 300.0);
        const int it = fem.solveLoad(Qf, u);
        const std::vector<double> r = fem.residual(u);

        double rint = 0.0, qsum = 0.0;
        std::vector<char> isIface(m.nNodes, 0);
        for (const int i : m.ifaceNodes) isIface[i] = 1;
        for (int i = 0; i < m.nNodes; i++) {
            if (isIface[i]) qsum += r[i];
            else            rint = std::max(rint, std::fabs(r[i]));
        }
        double tmin = u[0], tmax = u[0];
        for (const double t : u) { tmin = std::min(tmin, t); tmax = std::max(tmax, t); }
        std::cout << std::setprecision(10)
                  << "picard " << it << "  T " << tmin << " .. " << tmax
                  << "  sum(q_iface) " << qsum << " W/m  max|r_interior| " << rint << "\n";
        writeText(argv[4], u);
        return 0;
    }
    std::cerr << "[solid_fem2d_tool] unknown mode " << mode << "\n";
    return 1;
}
