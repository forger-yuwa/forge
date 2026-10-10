// 段 ② の診断 (plans/active/architecture-float-state-double-geometry.md §6.3)。説明と書く量は geomStage2Dump.hpp。
// 既定 off: main が環境変数 FORGE_DIAG_GEOM_STAGE2_DUMP を見たときだけ呼ぶ。

#include "cuda_forge/geomStage2Dump.hpp"
#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/calcGradient_d.cuh"
#include "cuda_forge/wallRepPoint_d.cuh"
#include "conjugateWall.hpp"

#include <highfive/H5File.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

namespace geomStage2Dump {
namespace {

template <typename T>
std::vector<T> d2h(const T* d, size_t n)
{
    std::vector<T> h(n);
    if (n > 0 && d != nullptr) gpuErrchk( cudaMemcpy(h.data(), d, n*sizeof(T), cudaMemcpyDeviceToHost) );
    return h;
}

template <typename T>
void put(HighFive::File& h5, const std::string& name, const std::vector<T>& v, const std::string& layout)
{
    if (v.empty()) {   // 長さ 0 は書き込みを呼ばずに空のデータセットだけ作る
        auto ds = h5.createDataSet<T>(name, HighFive::DataSpace(std::vector<size_t>{0}));
        if (!layout.empty()) ds.createAttribute("layout", layout);
        return;
    }
    auto ds = h5.createDataSet(name, v);
    if (!layout.empty()) ds.createAttribute("layout", layout);
}

template <typename T>
long long countBitDiff(const std::vector<T>& a, const std::vector<T>& b)
{
    if (a.size() != b.size()) return -1;
    long long n = 0;
    for (size_t i = 0; i < a.size(); ++i) if (std::memcmp(&a[i], &b[i], sizeof(T)) != 0) ++n;
    return n;
}

bool isWallKind(const bcond& bc) { return bc.bcondKind == "wall" || bc.bcondKind == "wall_isothermal"; }

}  // namespace

int run(const char* path, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    const size_t nP = (size_t)msh.nPlanes, nC = (size_t)msh.nCells, nCa = (size_t)msh.nCells_all;
    const bool node = (cfg.discretization == "node");
    printf("[geomStage2] FORGE_DIAG_GEOM_STAGE2_DUMP=%s: writing load-time geometry and connectivity (new = from double copies, "
           "legacy = rebuilt from geom_float/flow_float coordinates), then exit without update (flow_float %zu bytes, geom64 %s)\n",
           path, sizeof(flow_float), msh.hasGeom64() ? "yes" : "NO");

    HighFive::File h5(path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    h5.createAttribute("plan", std::string("plans/active/architecture-float-state-double-geometry.md §6.3 (stage 2)"));
    h5.createAttribute("flow_float_bytes", (long long)sizeof(flow_float));
    h5.createAttribute("has_geom64", (long long)(msh.hasGeom64() ? 1 : 0));
    h5.createAttribute("nCells", (long long)msh.nCells);
    h5.createAttribute("nCells_all", (long long)msh.nCells_all);
    h5.createAttribute("nPlanes", (long long)msh.nPlanes);
    h5.createAttribute("nNormalPlanes", (long long)msh.nNormalPlanes);
    h5.createAttribute("discretization", cfg.discretization);
    h5.createAttribute("gradLSQ", (long long)cfg.gradLSQ);
    h5.createAttribute("isAxisymmetric", (long long)cfg.isAxisymmetric);
    h5.createAttribute("axisymMethod", (long long)cfg.axisymMethod);
    h5.createAttribute("axisRFloor", (double)cfg.axisRFloor);
    h5.createAttribute("hoopAreaFromClosure", (long long)cfg.hoopAreaFromClosure);
    h5.createAttribute("lineImplicit", (long long)cfg.lineImplicit);
    h5.createAttribute("interfaceDiagAlignMin", (double)cfg.interfaceDiagAlignMin);
    h5.createAttribute("mesh_file", cfg.meshFileName);
    h5.createAttribute("input_value_file", cfg.valueFileName);

    // ---- 接続の表 (面・incidence の番号で突き合わせるため) ----
    std::vector<int> planeCells(2*nP, -1);
    for (size_t ip = 0; ip < nP; ++ip) {
        const auto& pc = msh.planes[ip].iCells;
        for (size_t k = 0; k < 2 && k < pc.size(); ++k) planeCells[2*ip + k] = pc[k];
    }
    put(h5, "/mesh/plane_cells", planeCells, "[2*ip+k] = planes[ip].iCells[k]");
    const std::vector<geom_int> cpi = d2h(msh.map_cell_planes_index_d, nC + 1);
    const size_t nInc = (size_t)cpi[nC];
    put(h5, "/mesh/cell_planes_index", cpi, "CSR [nCells+1]");
    put(h5, "/mesh/cell_planes", d2h(msh.map_cell_planes_d, nInc), "CSR [nInc]: plane id of incidence ilp");

    // ---- ラインの接続 (lineImplicit の有無によらず、本番と同じ関数で作る) ----
    {
        // 座標の退避先 (節点座標を持たない CV だけが読む) は本番 (main の buildImplicitLines) と同じホストの ccx..ccz。
        // ホストに無ければデバイスの写しを使う。
        std::vector<flow_float> cxh, cyh, czh;
        const flow_float* hx = var.c.at("ccx").data();
        const flow_float* hy = var.c.at("ccy").data();
        const flow_float* hz = var.c.at("ccz").data();
        if (var.c.at("ccx").size() < nC || var.c.at("ccy").size() < nC || var.c.at("ccz").size() < nC) {
            cxh = d2h(var.c_d.at("ccx"), nC); cyh = d2h(var.c_d.at("ccy"), nC); czh = d2h(var.c_d.at("ccz"), nC);
            hx = cxh.data(); hy = cyh.data(); hz = czh.data();
        }
        std::vector<geom_int> off[2], cel[2], prv[2], nxt[2];
        geom_int lenMax[2] = {0, 0};
        for (int k = 0; k < 2; ++k) msh.implicitLineTopology(k == 0, hx, hy, hz, off[k], cel[k], prv[k], nxt[k], lenMax[k]);
        const char* tag[2] = {"new", "legacy"};
        for (int k = 0; k < 2; ++k) {
            const std::string g = std::string("/lines/") + tag[k] + "/";
            put(h5, g + "offsets", off[k], "CSR [nLines+1]");
            put(h5, g + "cells", cel[k], "[sum of line lengths] CV id from the wall outwards");
            put(h5, g + "prev", prv[k], "[nCells] wall-side neighbour on the line (-1 = none)");
            put(h5, g + "next", nxt[k], "[nCells] outer neighbour on the line (-1 = none)");
        }
        long long devEq = -1;
        if (msh.line_offsets_d != nullptr) {
            const std::vector<geom_int> dOff = d2h(msh.line_offsets_d, (size_t)msh.nImplicitLines + 1);
            const std::vector<geom_int> dCel = d2h(msh.line_cells_d, (size_t)dOff.back());
            devEq = (countBitDiff(dOff, off[0]) == 0 && countBitDiff(dCel, cel[0]) == 0) ? 1 : 0;
        }
        const long long nLinesNew = (long long)off[0].size() - 1, nLinesOld = (long long)off[1].size() - 1;
        const bool same = (countBitDiff(off[0], off[1]) == 0 && countBitDiff(cel[0], cel[1]) == 0);
        h5.createAttribute("lines_built_in_run", (long long)(msh.line_offsets_d != nullptr ? 1 : 0));
        h5.createAttribute("lines_device_equal_new", devEq);
        h5.createAttribute("lines_new_equal_legacy", (long long)(same ? 1 : 0));
        printf("[geomStage2] lines: new %lld lines / %zu CVs, legacy %lld lines / %zu CVs, new == legacy: %s, device == new: %lld\n",
               nLinesNew, cel[0].size(), nLinesOld, cel[1].size(), same ? "yes" : "NO", devEq);
    }

    // ---- 周期の相手 ----
    {
        const std::map<int, std::vector<geom_int>> legacy = msh.periodicPartnerPlanes(false);
        long long nMis = 0, nTot = 0;
        for (const bcond& bc : msh.bconds) {
            if (bc.bcondKind != "periodic") continue;
            const std::string g = "/periodic/" + std::to_string(bc.physID) + "/";
            const auto itP = bc.bint.find("partnerPlnID");
            const auto itC = bc.bint.find("partnerCellID");
            if (itP != bc.bint.end()) put(h5, g + "partnerPlnID", itP->second, "[boundary face of this bcond] partner plane id (production)");
            if (itC != bc.bint.end()) put(h5, g + "partnerCellID", itC->second, "[boundary face of this bcond] partner cell id (production)");
            put(h5, g + "iPlanes", bc.iPlanes, "[boundary face of this bcond] plane id");
            const auto itL = legacy.find(bc.physID);
            if (itL != legacy.end()) {
                put(h5, g + "partnerPlnID_legacy", itL->second, "rebuilt from geom_float face centroids");
                if (itP != bc.bint.end()) {
                    for (size_t i = 0; i < itL->second.size() && i < itP->second.size(); ++i) {
                        ++nTot;
                        if (itL->second[i] != itP->second[i]) ++nMis;
                    }
                }
            }
        }
        if (!msh.periodicRoot.empty()) put(h5, "/periodic/periodicRoot", msh.periodicRoot, "[nCells] node periodic group root (production)");
        h5.createAttribute("periodic_faces_compared", nTot);
        h5.createAttribute("periodic_new_vs_legacy_mismatch", nMis);
        if (nTot > 0) printf("[geomStage2] periodic partners: %lld faces, new vs legacy mismatch %lld\n", nTot, nMis);
    }

    // ---- LSQ の係数 (gradLSQ 1/2、node) ----
    if (node && (cfg.gradLSQ == 1 || cfg.gradLSQ == 2)) {
        LsqDiagCoef cNew, cOld;
        lsqCoefForDiag(cfg, cuda_cfg, msh, var, false, cNew);
        lsqCoefForDiag(cfg, cuda_cfg, msh, var, true,  cOld);
        if (cfg.gradLSQ == 2) {
            const LsqCoefView v = lsq_coef_view();
            std::vector<flow_float> prod;
            if (v.cInt != nullptr && (size_t)v.nInc == nInc) prod = d2h(v.cInt, 3*nInc);
            const long long rebuildMis = prod.empty() ? -1 : countBitDiff(prod, cNew.cInt);
            put(h5, "/lsq/new/cInt", prod.empty() ? cNew.cInt : prod, "[3*ilp+k] production coefficients (lsq_coef_view)");
            put(h5, "/lsq/legacy/cInt", cOld.cInt, "[3*ilp+k] rebuilt with the stage-1 displacement (flow_float coordinate difference)");
            put(h5, "/lsq/new/seam_group", cNew.seamGroup, "[nInc] periodic seam group root (-1 = not on a seam group)");
            put(h5, "/lsq/new/seam_class", cNew.seamClass, "[nInc] equivalence class within the group (-1 = not on a seam group)");
            put(h5, "/lsq/legacy/seam_group", cOld.seamGroup, "[nInc]");
            put(h5, "/lsq/legacy/seam_class", cOld.seamClass, "[nInc]");
            h5.createAttribute("lsq_production_vs_rebuild_mismatch", rebuildMis);
            h5.createAttribute("lsq_nDegen_new", (long long)cNew.nDegen);
            h5.createAttribute("lsq_nDegen_legacy", (long long)cOld.nDegen);
            long long seamMis = 0;
            for (size_t i = 0; i < cNew.seamClass.size() && i < cOld.seamClass.size(); ++i)
                if (cNew.seamClass[i] != cOld.seamClass[i] || cNew.seamGroup[i] != cOld.seamGroup[i]) ++seamMis;
            h5.createAttribute("lsq_seam_new_vs_legacy_mismatch", seamMis);
            printf("[geomStage2] LSQ (gradLSQ=2): %zu incidences, production vs rebuild (new) bit mismatch %lld, "
                   "truncated new/legacy %d/%d, seam class mismatch new vs legacy %lld\n",
                   nInc, rebuildMis, cNew.nDegen, cOld.nDegen, seamMis);
        } else {
            put(h5, "/lsq/new/M6", cNew.M6, "[k*nCells+ic], k = Mxx,Mxy,Mxz,Myy,Myz,Mzz (flow_float as stored)");
            put(h5, "/lsq/legacy/M6", cOld.M6, "same layout, stage-1 displacement");
            printf("[geomStage2] LSQ (gradLSQ=1): M written for %zu CVs\n", nC);
        }
    }

    // ---- 壁の代表内点 (node、全壁 bcond。本番の配列 + 旧版の作り直し) ----
    if (node && msh.wall_flag_d != nullptr) {
        long long nTot = 0, nIrepMis = 0, nNoCandNew = 0, nNoCandOld = 0;
        for (const bcond& bc : msh.bconds) {
            if (!isWallKind(bc) || bc.iPlanes.empty()) continue;
            const WallRepPoints& w = wallRepPoints(cuda_cfg, msh, var, bc);
            const size_t nb = (size_t)w.nb;
            const std::string g = "/wall/" + std::to_string(bc.physID) + "/";
            const std::vector<geom_int> irep = d2h(w.irep_d, nb);
            put(h5, g + "iPlanes", bc.iPlanes, "[ib] boundary plane id");
            put(h5, g + "iCells", bc.iCells, "[ib] wall node (CV) id");
            put(h5, g + "irep", irep, "[ib] representative interior node (-1 = no candidate), production");
            put(h5, g + "dn", d2h(w.dn_d, nb), "[ib] inward wall distance; kernels use y = max(dn, 1e-12)");
            put(h5, g + "dist", d2h(w.dist_d, nb), "[ib] |x_I - x_W|");
            put(h5, g + "cos", d2h(w.cos_d, nb), "[ib] dn/|d| (-2 = no candidate)");
            std::vector<geom_int> li; std::vector<flow_float> ldn, ldist, lcos;
            wallRepPointsHost(cuda_cfg, msh, var, bc, true, li, ldn, ldist, lcos);
            put(h5, g + "irep_legacy", li, "stage-1 selection (flow_float coordinate difference)");
            put(h5, g + "dn_legacy", ldn, "");
            put(h5, g + "dist_legacy", ldist, "");
            put(h5, g + "cos_legacy", lcos, "");
            for (size_t i = 0; i < nb && i < li.size(); ++i) {
                ++nTot;
                if (irep[i] != li[i]) ++nIrepMis;
                if (irep[i] < 0) ++nNoCandNew;
                if (li[i] < 0) ++nNoCandOld;
            }
        }
        h5.createAttribute("wall_faces", nTot);
        h5.createAttribute("wall_irep_new_vs_legacy_mismatch", nIrepMis);
        printf("[geomStage2] wall representative points: %lld faces, irep new vs legacy mismatch %lld, "
               "no candidate new/legacy %lld/%lld\n", nTot, nIrepMis, nNoCandNew, nNoCandOld);
    }

    // ---- 軸対称の closure (method 0) ----
    if (cfg.isAxisymmetric == 1 && cfg.axisymMethod == 0) {
        const bool active = (cfg.axisRFloor > (flow_float)0.0 || cfg.hoopAreaFromClosure == 1);
        h5.createAttribute("closure_active", (long long)(active ? 1 : 0));
        put(h5, "/axisym/A_closure_x", d2h(var.c_d.at("A_closure_x"), nCa), "[nCells_all] device (used by the hoop source)");
        put(h5, "/axisym/A_closure_y", d2h(var.c_d.at("A_closure_y"), nCa), "[nCells_all] device");
        put(h5, "/axisym/A_planar", d2h(var.c_d.at("A_planar"), nCa), "[nCells_all] device");
        // 診断: デバイスの最終の面ベクトル (半径の重みの後) の double の和 (丸める前)・絶対値の和、生の double の幾何からの closure
        const std::vector<flow_float> sx = d2h(var.p_d.at("sx"), nP), sy = d2h(var.p_d.at("sy"), nP);
        std::vector<double> sumx(nCa, 0.0), sumy(nCa, 0.0), absx(nCa, 0.0), absy(nCa, 0.0), rawx(nCa, 0.0), rawy(nCa, 0.0);
        const double rFloor = (cfg.axisRFloor > (flow_float)0.0) ? (double)cfg.axisRFloor : 1.0e-20;
        const bool has64 = msh.hasGeom64();
        for (size_t ip = 0; ip < nP; ++ip) {
            const auto& pc = msh.planes[ip].iCells;
            if (pc.empty()) continue;
            double rx = 0.0, ry = 0.0;
            if (has64) {
                const double r = std::max(msh.planeCent64[3*ip + 1], rFloor);
                rx = msh.surfVect64[3*ip + 0]*r;
                ry = msh.surfVect64[3*ip + 1]*r;
            }
            for (size_t k = 0; k < 2 && k < pc.size(); ++k) {
                const geom_int ic = pc[k];
                if (ic < 0 || ic >= msh.nCells) continue;
                const double sg = (k == 0) ? 1.0 : -1.0;
                sumx[ic] += sg*(double)sx[ip]; sumy[ic] += sg*(double)sy[ip];
                absx[ic] += std::fabs((double)sx[ip]); absy[ic] += std::fabs((double)sy[ip]);
                rawx[ic] += sg*rx; rawy[ic] += sg*ry;
            }
        }
        put(h5, "/axisym/closure_sum64_x", sumx, "sum of +-S_device in double before the final rounding");
        put(h5, "/axisym/closure_sum64_y", sumy, "");
        put(h5, "/axisym/closure_abs64_x", absx, "sum of |S_device| (normalisation before cancellation)");
        put(h5, "/axisym/closure_abs64_y", absy, "");
        if (has64) {
            put(h5, "/axisym/closure_raw64_x", rawx, "diagnostic: closure from raw double surfVect and max(planeCent_y, rFloor)");
            put(h5, "/axisym/closure_raw64_y", rawy, "");
        }
        printf("[geomStage2] axisymmetric closure written (active %d)\n", active ? 1 : 0);
    }

    // ---- 第一・第二内部点 d1/d2 (全壁 bcond) ----
    {
        int nb = 0;
        for (const bcond& bc : msh.bconds) {
            if (!isWallKind(bc) || bc.iPlanes.empty()) continue;
            const std::string g = "/d1d2/" + std::to_string(bc.physID) + "/";
            for (int k = 0; k < 2; ++k) {
                const conjugateWall::FirstInterior fi = conjugateWall::firstInteriorForDiag(cfg, msh, bc, k == 1);
                const std::string sfx = (k == 0) ? "" : "_legacy";
                std::vector<int> ok(fi.ok.begin(), fi.ok.end());
                put(h5, g + "d1" + sfx, fi.d1, k == 0 ? "[ib] normal distance to the first interior DOF (NaN = not resolved)" : "");
                put(h5, g + "d2" + sfx, fi.d2, k == 0 ? "[ib] normal distance to the second interior DOF" : "");
                put(h5, g + "jdof" + sfx, fi.jdof, "");
                put(h5, g + "jdof2" + sfx, fi.jdof2, "");
                put(h5, g + "ok" + sfx, ok, "");
            }
            ++nb;
        }
        printf("[geomStage2] first/second interior points written for %d wall bconds\n", nb);
    }

    // ---- delta_les (本番のデバイスの値と、段 ① までの式の作り直し) ----
    {
        put(h5, "/delta_les", d2h(var.c_d.at("delta_les"), nCa), "[nCells_all] device");
        const std::vector<flow_float> cx = d2h(var.c_d.at("ccx"), nCa), cy = d2h(var.c_d.at("ccy"), nCa), cz = d2h(var.c_d.at("ccz"), nCa);
        std::vector<flow_float> dl(nCa, (flow_float)0.0);
        for (size_t ip = 0; ip < nP; ++ip) {
            const geom_int ic1 = msh.planes[ip].iCells[0];
            const geom_int ic2 = msh.planes[ip].iCells[1];
            const geom_float dx = cx[ic1] - cx[ic2];
            const geom_float dy = cy[ic1] - cy[ic2];
            const geom_float dz = cz[ic1] - cz[ic2];
            const geom_float d  = sqrt(dx*dx + dy*dy + dz*dz);
            if (ic1 < msh.nCells && d > dl[ic1]) dl[ic1] = d;
            if (ic2 < msh.nCells && d > dl[ic2]) dl[ic2] = d;
        }
        put(h5, "/delta_les_legacy", dl, "stage-1 formula on the flow_float coordinates");
    }

    printf("[geomStage2] wrote %s; exiting without any update\n", path);
    fflush(stdout);
    return 0;
}

}  // namespace geomStage2Dump
