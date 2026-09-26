// CHT Phase 2: 固体 FE (`fem2d`) の**同値試験用 CLI** (plan §6 V4b (a))。
//
// ソルバ本体に組み込む前に、参照実装 (tools/solid_fem2d.py) と**同じ問題を同じ精度で**
// 解けることをここで閉じる。突き合わせは tools/test_solid_fem2d_cpp.py が行う。
//
//   solid_fem2d_tool [--axisym] info   <solid.h5>
//   solid_fem2d_tool [--axisym] matvec <solid.h5> <u.txt>  <out.txt>   … K(u) u を書く (a1)
//   solid_fem2d_tool [--axisym] solve  <solid.h5> <Qf.txt> <out.txt>   … 界面荷重で解いて u を書く (a2)
//   solid_fem2d_tool [--axisym] matrix <solid.h5> <u.txt>  <out.txt>   … K(u) の下三角 "i j A(i,j)" の後に
//                                                                        "b i b_i" (Robin 荷重) を書く
//   solid_fem2d_tool [--axisym] lumped <solid.h5> <out.txt>             … 界面の集中量 (IFACE/NODES 順)
//   solid_fem2d_tool [--axisym] field  <solid.h5> <u.txt>  <stem>       … writeSolidField で <stem>.h5/.xmf
//
// テキストは 1 行 1 値 (倍精度)。`solve` の Qf は **IFACE/NODES の順**、単位は W/m。
// `--axisym` は r = y の重みを入れる (単位は W/rad。plan boundary-cht-axisymmetric-fem2d §4.2)。
// 固体 h5 には持たない (content_sha1 を変えない) ので、ここで SolidMesh::axisym を立てる。

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

int main(int argc0, char** argv0)
{
    // `--axisym` はどこに置いてもよい (取り除いて残りを位置引数として読む)
    bool axisym = false;
    std::vector<char*> args;
    for (int i = 0; i < argc0; i++) {
        if (std::string(argv0[i]) == "--axisym") axisym = true;
        else args.push_back(argv0[i]);
    }
    const int argc = (int)args.size();
    char** argv = args.data();
    if (argc < 3) {
        std::cerr << "usage: solid_fem2d_tool [--axisym] info|matvec|solve|matrix|lumped|field "
                     "<solid.h5> [in.txt out.txt]\n";
        return 1;
    }
    const std::string mode = argv[1];
    conjugate::SolidMesh m = conjugate::SolidMesh::read(argv[2]);
    m.axisym = axisym;
    conjugate::SolidFem2D fem(m);

    if (mode == "info") {
        std::cout << (axisym ? "axisym  " : "planar  ") << "nodes " << m.nNodes << "  tris " << m.nTris
                  << "  iface " << m.nIface() << "  robin " << m.robinH.size()
                  << "  bandwidth " << fem.bandwidth() << " (attr " << m.bandwidth << ")\n"
                  << "k_s(300) " << m.kOf(300.0) << "  k_s(900) " << m.kOf(900.0)
                  << "  iface_sha1 " << m.ifaceSha1 << "\n";
        return 0;
    }
    if (mode == "lumped") {
        if (argc < 4) { std::cerr << "usage: ... lumped <solid.h5> <out.txt>\n"; return 1; }
        writeText(argv[3], m.ifaceLumped());
        return 0;
    }
    if (argc < 5) { std::cerr << "usage: ... <in.txt> <out.txt>\n"; return 1; }

    if (mode == "matrix" || mode == "field") {
        const std::vector<double> u = readText(argv[3]);
        if ((int)u.size() != m.nNodes) {
            std::cerr << "[solid_fem2d_tool] u の長さ " << u.size() << " != nodes " << m.nNodes << "\n";
            return 1;
        }
        if (mode == "field") {
            conjugate::writeSolidField(argv[4], m, u, {}, 0.0);
            return 0;
        }
        fem.assemble(u);
        std::ofstream ofs(argv[4]);
        ofs << std::setprecision(17) << std::scientific;
        const int bw = fem.bandwidth();
        for (int j = 0; j < m.nNodes; j++)
            for (int i = j; i <= std::min(m.nNodes - 1, j + bw); i++) {
                const double a = fem.entry(i, j);
                if (a != 0.0) ofs << i << " " << j << " " << a << "\n";
            }
        const std::vector<double>& b = fem.rhs();
        for (int i = 0; i < m.nNodes; i++) ofs << "b " << i << " " << b[i] << "\n";
        return 0;
    }

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
