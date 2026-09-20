#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <vector>
#include <cmath>

#include "weakIsothermalWall_d.cuh"
#include "../conjugateWall.hpp"
#include "cudaWrapper.cuh"

namespace weakIsoWall {

bool active(const solverConfig& cfg, const mesh& msh)
{
    if (cfg.nodeIsothermalEnergyBC == 0) return false;
    if (cfg.discretization != "node") return false;
    for (const auto& bc : msh.bconds) if (bc.bcondKind == "wall_isothermal") return true;
    return false;
}

static void die(const char* what)
{
    std::fprintf(stderr,
        "[nodeIsothermalEnergyBC] %s\n"
        "  mesh.nodeIsothermalEnergyBC: 1 (弱形式) は初版では構成を限定している。\n"
        "  仕様: methods/boundary.md「等温壁のエネルギー境界: 強制 と 弱形式」\n", what);
    std::exit(1);
}

void validate(const solverConfig& cfg, const mesh& msh)
{
    if (cfg.nodeIsothermalEnergyBC == 0) return;
    if (cfg.discretization != "node")
        die("cell 方式では使えない (node 限定)。");
    if (cfg.nodeWallDirichlet == 0)
        die("nodeWallDirichlet: 1 が必須 (0 では運動量まで弱くなる)。");
    if (cfg.thermalMethod != 0)
        die("初版は thermalMethod: 0 (CPG) 限定。TP は e(T) なので対角の CPG 式を一般化できない。");
    if (cfg.wallTreatmentSST == 1)
        die("wallTreatmentSST: 1 (壁関数) とは併用不可。");
    if (cfg.isAxisymmetric == 1)
        die("軸対称とは併用不可 (r 重みが未対応)。");

    bool any = false;
    for (const auto& bc : msh.bconds) {
        if (bc.bcondKind != "wall_isothermal") continue;
        any = true;
        auto itw = bc.inputInts.find("wallModelLES");
        if (itw != bc.inputInts.end() && itw->second != 0)
            die("wallModelLES: 1 の壁とは併用不可。");
        for (const char* k : {"Ux", "Uy", "Uz"}) {
            auto it = bc.inputFloats.find(k);
            if (it != bc.inputFloats.end() && std::fabs((double)it->second) > 0.0)
                die("移動壁 (wall_isothermal の Ux/Uy/Uz != 0) とは併用不可。");
        }
    }
    if (!any) die("wall_isothermal の bcond が無い。");
}

const Geom& geom(const solverConfig& cfg, const mesh& msh, const bcond& bc)
{
    static std::map<geom_int, Geom> cache;
    auto it = cache.find(bc.physID);
    if (it != cache.end()) return it->second;

    Geom g;
    const geom_int nbp = (geom_int)bc.iPlanes.size();
    g.n = nbp;
    if (nbp == 0) { cache[bc.physID] = g; return cache[bc.physID]; }

    // **診断 q_compact と同一規約**の第一内部点 (conjugateWall::firstInterior)。
    const conjugateWall::FirstInterior& fi = conjugateWall::firstInterior(cfg, msh, bc);

    // 壁 DOF 集合 (全 wall 種 bcond)。第一内部点が別の壁ノードなら評価不能にする。
    std::set<geom_int> wallDof;
    for (const auto& b2 : msh.bconds) {
        if (b2.bcondKind != "wall" && b2.bcondKind != "wall_isothermal") continue;
        for (const geom_int ic : b2.iCells) wallDof.insert(ic);
    }

    std::vector<geom_int>   jh(nbp);
    std::vector<flow_float> dh(nbp);
    geom_int nbad = 0;
    for (geom_int ib = 0; ib < nbp; ib++) {
        const geom_int j = fi.jdof[ib];
        const double   d = fi.d1[ib];
        const bool ok = (fi.ok[ib] != 0) && (j >= 0) && (j < msh.nCells)
                     && std::isfinite(d) && (d > 0.0)
                     && (wallDof.count(j) == 0);   // 別の壁ノードを選んでいない
        if (!ok) {
            if (nbad < 5)
                std::fprintf(stderr, "[nodeIsothermalEnergyBC] physID %d bplane %ld: "
                                     "第一内部点が評価不能 (j=%ld d1=%g align=%g)\n",
                             (int)bc.physID, (long)ib, (long)j, d, fi.align[ib]);
            nbad++;
            jh[ib] = -1; dh[ib] = (flow_float)0.0;
        } else {
            jh[ib] = j; dh[ib] = (flow_float)d;
        }
    }
    if (nbad > 0) {
        std::fprintf(stderr, "[nodeIsothermalEnergyBC] physID %d: 評価不能な壁半割面が %ld / %ld ある。\n"
                             "  弱形式は第一内部点が定まらない面を扱えない (合格に変換しない)。\n",
                     (int)bc.physID, (long)nbad, (long)nbp);
        std::exit(1);
    }

    gpuErrchk( cudaMalloc(&g.j_d , sizeof(geom_int)*nbp) );
    gpuErrchk( cudaMalloc(&g.d1_d, sizeof(flow_float)*nbp) );
    gpuErrchk( cudaMemcpy(g.j_d , jh.data(), sizeof(geom_int)*nbp  , cudaMemcpyHostToDevice) );
    gpuErrchk( cudaMemcpy(g.d1_d, dh.data(), sizeof(flow_float)*nbp, cudaMemcpyHostToDevice) );

    double dmin = 1e30, dmax = -1e30;
    for (geom_int ib = 0; ib < nbp; ib++) { dmin = std::fmin(dmin, dh[ib]); dmax = std::fmax(dmax, dh[ib]); }
    std::printf("[nodeIsothermalEnergyBC] physID %d: %ld 壁半割面, d1 (法線投影) %.3e .. %.3e m\n",
                (int)bc.physID, (long)nbp, dmin, dmax);

    cache[bc.physID] = g;
    return cache[bc.physID];
}


// --- 陰解法の近似対角項の素材 -------------------------------------------------
namespace {
    flow_float* g_diag = nullptr;
    geom_int    g_n    = 0;
}

flow_float* diagBuf(const mesh& msh)
{
    if (g_diag == nullptr) {
        g_n = msh.nCells;
        gpuErrchk( cudaMalloc(&g_diag, sizeof(flow_float)*g_n) );
        gpuErrchk( cudaMemset(g_diag, 0, sizeof(flow_float)*g_n) );
    }
    return g_diag;
}

void diagReset(const mesh& msh)
{
    flow_float* p = diagBuf(msh);
    gpuErrchk( cudaMemset(p, 0, sizeof(flow_float)*g_n) );
}

}  // namespace weakIsoWall
