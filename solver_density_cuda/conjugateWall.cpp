#include "conjugateWall.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <memory>
#include <set>
#include <unordered_map>
#include <utility>

#include <highfive/highfive.hpp>

#include "conjugate/solidFem2d.hpp"
#include "cuda_forge/cudaWrapper.cuh"

namespace conjugateWall {

namespace {

// DOF -> 隣接 DOF (内部面が繋ぐ 2 つの CV)。メッシュ全体で 1 度だけ作る。
const std::vector<std::vector<geom_int>>& dofNeighbors(const mesh& msh)
{
    static std::vector<std::vector<geom_int>> nb;
    static const mesh* built_for = nullptr;
    if (built_for == &msh && !nb.empty()) return nb;

    nb.assign(msh.nCells, {});
    for (geom_int ip = 0; ip < (geom_int)msh.planes.size(); ip++) {
        const auto& ic = msh.planes[ip].iCells;
        if (ic.size() < 2) continue;
        const geom_int a = ic[0], b = ic[1];
        if (a < 0 || b < 0 || a >= msh.nCells || b >= msh.nCells) continue;  // ゴーストは除く
        nb[a].push_back(b);
        nb[b].push_back(a);
    }
    built_for = &msh;
    return nb;
}

// 値の位置 (node モードはノード座標 = T の DOF 位置、cell モードは CV 重心)。
inline void dofCoords(const mesh& msh, bool nodeMode, geom_int ic, double xyz[3])
{
    if (nodeMode && (geom_int)msh.nodes.size() > ic && msh.nodes[ic].coords.size() >= 3) {
        xyz[0] = msh.nodes[ic].coords[0];
        xyz[1] = msh.nodes[ic].coords[1];
        xyz[2] = msh.nodes[ic].coords[2];
    } else {
        xyz[0] = msh.cells[ic].centCoords[0];
        xyz[1] = msh.cells[ic].centCoords[1];
        xyz[2] = msh.cells[ic].centCoords[2];
    }
}

std::vector<flow_float> pullField(const solverConfig& cfg, variables& var, const std::string& name, geom_int n)
{
    std::vector<flow_float> h(n, (flow_float)0.0);
    const auto itd = var.c_d.find(name);
    if (cfg.gpu == 1 && itd != var.c_d.end() && itd->second != nullptr) {
        cudaWrapper::cudaMemcpy_D2H_wrapper(itd->second, h.data(), n);
        return h;
    }
    const auto ith = var.c.find(name);
    if (ith != var.c.end() && (geom_int)ith->second.size() >= n) {
        std::copy(ith->second.begin(), ith->second.begin() + n, h.begin());
    }
    return h;
}

} // namespace

const FirstInterior& firstInterior(const solverConfig& cfg, const mesh& msh, const bcond& bc)
{
    const bool nodeMode  = (cfg.discretization == "node");
    const double alignMin = cfg.interfaceDiagAlignMin;
    static std::map<geom_int, FirstInterior> cache;
    const auto it = cache.find(bc.physID);
    if (it != cache.end()) return it->second;

    const auto& nb = dofNeighbors(msh);
    const geom_int nbp = (geom_int)bc.iPlanes.size();

    // 壁 DOF ごとの法線 = 属する境界面ベクトルの合算 (角では平均法線)。
    std::unordered_map<geom_int, std::array<double,3>> nrm;
    for (geom_int ib = 0; ib < nbp; ib++) {
        const geom_int ic = bc.iCells[ib];
        const geom_int ip = bc.iPlanes[ib];
        auto& v = nrm[ic];
        v[0] += msh.planes[ip].surfVect[0];
        v[1] += msh.planes[ip].surfVect[1];
        v[2] += msh.planes[ip].surfVect[2];
    }

    FirstInterior fi;
    fi.jdof.assign(nbp, -1);
    fi.jdof2.assign(nbp, -1);
    fi.d1.assign(nbp, std::numeric_limits<double>::quiet_NaN());
    fi.d2.assign(nbp, std::numeric_limits<double>::quiet_NaN());
    fi.align.assign(nbp, 0.0);
    fi.nx.assign(nbp, 0.0); fi.ny.assign(nbp, 0.0); fi.nz.assign(nbp, 0.0);
    fi.ok.assign(nbp, 0);

    for (geom_int ib = 0; ib < nbp; ib++) {
        const geom_int ic = bc.iCells[ib];
        if (ic < 0 || ic >= msh.nCells) continue;
        const auto& v = nrm[ic];
        const double ln = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        if (!(ln > 0.0)) continue;
        const double nh[3] = { v[0]/ln, v[1]/ln, v[2]/ln };
        fi.nx[ib] = nh[0]; fi.ny[ib] = nh[1]; fi.nz[ib] = nh[2];

        double xw[3]; dofCoords(msh, nodeMode, ic, xw);
        double bestAl = -1.0, bestD = 0.0; geom_int bestJ = -1;
        for (const geom_int j : nb[ic]) {
            double xj[3]; dofCoords(msh, nodeMode, j, xj);
            const double d[3] = { xj[0]-xw[0], xj[1]-xw[1], xj[2]-xw[2] };
            const double dn = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
            if (!(dn > 0.0)) continue;
            const double proj = std::fabs(d[0]*nh[0] + d[1]*nh[1] + d[2]*nh[2]);
            const double al = proj/dn;
            if (al > bestAl) { bestAl = al; bestD = proj; bestJ = j; }
        }
        fi.align[ib] = (bestAl > 0.0) ? bestAl : 0.0;
        if (bestJ >= 0 && bestAl >= alignMin && bestD > 0.0) {
            fi.jdof[ib] = bestJ;
            fi.d1[ib]   = bestD;
            fi.ok[ib]   = 1;
            // 第二内部点: 第一内部点の隣で、同じ法線に沿ってさらに遠いもの (2 次片側差分用)。
            double best2Al = -1.0, best2D = 0.0; geom_int best2J = -1;
            for (const geom_int j2 : nb[bestJ]) {
                if (j2 == ic) continue;
                double x2[3]; dofCoords(msh, nodeMode, j2, x2);
                const double d2v[3] = { x2[0]-xw[0], x2[1]-xw[1], x2[2]-xw[2] };
                const double dn2 = std::sqrt(d2v[0]*d2v[0] + d2v[1]*d2v[1] + d2v[2]*d2v[2]);
                if (!(dn2 > 0.0)) continue;
                const double proj2 = std::fabs(d2v[0]*nh[0] + d2v[1]*nh[1] + d2v[2]*nh[2]);
                if (proj2 <= bestD) continue;              // 壁から遠い側だけ
                const double al2 = proj2/dn2;
                if (al2 > best2Al) { best2Al = al2; best2D = proj2; best2J = j2; }
            }
            if (best2J >= 0 && best2Al >= alignMin) { fi.jdof2[ib] = best2J; fi.d2[ib] = best2D; }
        }
    }

    geom_int nok = 0;
    for (geom_int ib = 0; ib < nbp; ib++) nok += (fi.ok[ib] ? 1 : 0);
    std::cout << "[interfaceDiag] physID=" << bc.physID << " (" << bc.physName << "): "
              << "first interior point resolved for " << nok << " / " << nbp
              << " wall DOFs (alignMin=" << alignMin << ")" << std::endl;

    return cache.emplace(bc.physID, std::move(fi)).first->second;
}

void fillInterfaceDiagnostics(const solverConfig& cfg, const mesh& msh, variables& var, bcond& bc)
{
    if (cfg.interfaceDiag != 1) return;
    if (!(bc.bcondKind == "wall" || bc.bcondKind == "wall_isothermal")) return;

    const geom_int nbp = (geom_int)bc.iPlanes.size();
    if (nbp <= 0) return;

    const FirstInterior& fi = firstInterior(cfg, msh, bc);

    const std::vector<flow_float> T        = pullField(cfg, var, "T",         msh.nCells);
    const std::vector<flow_float> thermCond= pullField(cfg, var, "thermCond", msh.nCells);
    const std::vector<flow_float> cp       = pullField(cfg, var, "cp",        msh.nCells);
    const std::vector<flow_float> visTurb  = pullField(cfg, var, "vis_turb",  msh.nCells);

    std::vector<flow_float> d_T1(nbp, 0.0), d_d1(nbp, 0.0), d_keff(nbp, 0.0);
    std::vector<flow_float> d_qc(nbp, 0.0), d_qr(nbp, 0.0), d_ok(nbp, 0.0), d_al(nbp, 0.0);
    std::vector<flow_float> d_q2(nbp, 0.0);
    // M4: 指定壁温と保存壁温を**別々に**記録する (弱形式では一致しない)。
    std::vector<flow_float> d_Twbc(nbp, 0.0), d_Tnode(nbp, 0.0), d_qcn(nbp, 0.0);
    const auto itTs = bc.bvar.find("Ts");
    const bool hasTsb = (itTs != bc.bvar.end() && (geom_int)itTs->second.size() >= nbp);

    const auto itq = bc.bvar.find("qwall");
    const bool hasQwall = (itq != bc.bvar.end() && (geom_int)itq->second.size() >= nbp);

    // ---- 保存的な実効界面熱量 $Q_f = \sum F^E - C$ (plan boundary-conjugate-heat-transfer §4.3) ----
    // 定常の Dirichlet 行では $C=-R^{raw}$ なので、壁ノード 1 点あたり
    //   $Q_f = \sum F^E_{\partial w} + R^{raw}$,  $\sum F^E_{\partial w} = -(\text{壁面が res\_roe に入れた寄与})$
    // である。`ifaceFw` が後者、`ifaceRraw` が前者の素材。面積で割って [W/m²] にする。
    // **適用範囲** (超えたら NaN を出す。散文で保証しない):
    //   - 定常 (`unsteady: 0`) のみ。dual-time は $C=D_t(VE)-R^{raw}$ で式が違う (§4.3)。
    //   - node の等温壁 Dirichlet ピンがある壁のみ (断熱壁には $R^{raw}$ が無い)。
    //   - 周期境界に属する壁ノードは除く (root 単位の 1 回集計が要る。依存 plan の解除待ち)。
    std::vector<flow_float> d_qeff(nbp, std::numeric_limits<flow_float>::quiet_NaN());
    // 未収束の擬似時間では壁 CV の質量残差 $R_\rho$ が 0 でなく、壁 CV (固定体積・$u=0$・$T=T_w$) のエネルギーが
    // $d(V\rho e_w)/d\tau = e_w R_\rho$ だけ変わる。これは蓄積であって固体へ渡る熱ではないので `iface_q_eff` から引く
    // (dual-time の $D_t(VE)$ と同じ位置づけ。plan boundary-conjugate-heat-transfer §5.1 #61)。
    // **係数は $e_w=\rho E/\rho$ であって $H_w$ ではない** (codex 2026-09-21): 対流が運ぶのは $H_wR_\rho$ だが、
    // 差 $(P/\rho)R_\rho$ は等温のまま質量を押し込む流動仕事で、壁が実際に受け取る熱である。
    // 引く前の値は `iface_q_eff_raw`。定常 ($R_\rho=0$) では一致する。block-DPLUR の実際の更新量とは
    // 一般に一致しないので「半離散式に基づく推定」である。
    std::vector<flow_float> d_qeff_raw(nbp, std::numeric_limits<flow_float>::quiet_NaN());
    const auto itRo = bc.bvar.find("ifaceRro");
    const bool hasRro = (itRo != bc.bvar.end() && (geom_int)itRo->second.size() >= nbp);
    std::vector<flow_float> f_ro, f_roe;
    if (hasRro) {
        f_ro  = pullField(cfg, var, "ro",  msh.nCells);
        f_roe = pullField(cfg, var, "roe", msh.nCells);
    }
    const auto itFw = bc.bvar.find("ifaceFw");
    const auto itRr = bc.bvar.find("ifaceRraw");
    // 弱形式 (mesh.nodeIsothermalEnergyBC=1) では壁エネルギーを拘束しないので **C = 0**。
    // したがって $Q_f = \sum F^E$ がそのまま収支に一致し、`ifaceRraw` は不要 (読んでも 0)。
    // ただし $C=0$ は流体側の定常性も固体との収支も保証しない。G-cons と収束判定を別途通すこと
    // (plan boundary-weak-isothermal-wall §4.5)。
    const bool weakIso = (cfg.nodeIsothermalEnergyBC != 0) && (cfg.discretization == "node");
    const bool hasEff = (cfg.unsteady == 0)
                     && (bc.bcondKind == "wall_isothermal")
                     && (itFw != bc.bvar.end() && (geom_int)itFw->second.size() >= nbp)
                     && (weakIso || (itRr != bc.bvar.end() && (geom_int)itRr->second.size() >= nbp));
    if (hasEff) {
        const bool hasPer = !msh.periodicRoot.empty();
        for (geom_int ib = 0; ib < nbp; ib++) {
            const geom_int ic = bc.iCells[ib];
            if (hasPer && ic < (geom_int)msh.periodicRoot.size() && msh.periodicRoot[ic] != ic) continue;
            const geom_int ip = bc.iPlanes[ib];
            const double area = (double)msh.planes[ip].surfArea;
            if (!(area > 0.0)) continue;
            const double Fw   = (double)itFw->second[ib];    // 壁面が res_roe に入れた寄与
            const double Rraw = weakIso ? 0.0                // 弱形式: 拘束が無いので C=0
                                        : (double)itRr->second[ib];    // 射影前の res_roe
            d_qeff_raw[ib] = (flow_float)((Rraw - Fw) / area);   // 固体向き正
            double store = 0.0;                                   // 壁 CV の蓄積 $e_w R_\rho$ (弱形式は拘束なしなので対象外)
            if (hasRro && !weakIso && f_ro[ic] > (flow_float)0.0)
                store = (double)f_roe[ic] / (double)f_ro[ic] * (double)itRo->second[ib];
            d_qeff[ib] = (flow_float)((Rraw - Fw - store) / area);
        }
    }

    for (geom_int ib = 0; ib < nbp; ib++) {
        const geom_int ic = bc.iCells[ib];
        d_al[ib] = (flow_float)fi.align[ib];
        // k_eff は viscousFlux の壁経路と同じ定義 (k_lam + cp mu_t / Pr_t)。
        const flow_float keff = thermCond[ic] + cp[ic]*visTurb[ic]/cfg.turbulentPrandtl;
        d_keff[ib] = keff;
        // q_recon: qwall は「壁→流体が正」なので固体向き正へ反転する。
        d_qr[ib] = hasQwall ? (flow_float)(-itq->second[ib]) : (flow_float)0.0;

        if (!fi.ok[ib]) { d_ok[ib] = 0.0; continue; }
        const geom_int j = fi.jdof[ib];
        d_ok[ib] = 1.0;
        d_T1[ib] = T[j];
        d_d1[ib] = (flow_float)fi.d1[ib];
        // q_compact: 固体向き正 (T_1 > T_w なら流体から固体へ)。
        // **壁温は指定値 (bvar Ts) を使う** — 強制側は T[ic]=Tw なので従来とビット同一だが、
        // 弱形式では T[ic] が自由 DOF なので保存温度を使うと境界流束と別量になる
        // (codex plan レビュー M4)。保存温度からの勾配は `iface_q_compact_Tnode` に分ける。
        const double Tw_bc = hasTsb ? (double)itTs->second[ib] : (double)T[ic];
        d_Twbc[ib] = (flow_float)Tw_bc;
        d_Tnode[ib] = T[ic];
        d_qc[ib]  = (flow_float)(keff * ((double)T[j] - Tw_bc)   / fi.d1[ib]);
        d_qcn[ib] = (flow_float)(keff * (double)(T[j] - T[ic])   / fi.d1[ib]);
        // q_2nd: 3 点非等間隔の 2 次片側差分 (壁 0 / 第一内部点 a / 第二内部点 b)。
        //   f'(0) = -((a+b)/(a b)) f0 + (b/(a(b-a))) f1 - (a/(b(b-a))) f2
        if (fi.jdof2[ib] >= 0) {
            const double a = fi.d1[ib], b = fi.d2[ib];
            if (b > a && a > 0.0) {
                const double dTdn = -((a+b)/(a*b))*Tw_bc
                                  + (b/(a*(b-a)))*(double)T[j]
                                  - (a/(b*(b-a)))*(double)T[fi.jdof2[ib]];
                d_q2[ib] = (flow_float)(keff * dTdn);
            }
        }
    }

    bc.diagVar["iface_T1"]        = std::move(d_T1);
    bc.diagVar["iface_d1"]        = std::move(d_d1);
    bc.diagVar["iface_keff"]      = std::move(d_keff);
    bc.diagVar["iface_Tw_bc"]     = std::move(d_Twbc);
    bc.diagVar["iface_T_node"]    = std::move(d_Tnode);
    bc.diagVar["iface_q_compact"] = std::move(d_qc);
    bc.diagVar["iface_q_compact_Tnode"] = std::move(d_qcn);
    bc.diagVar["iface_q_recon"]   = std::move(d_qr);
    bc.diagVar["iface_q_2nd"]     = std::move(d_q2);
    bc.diagVar["iface_ok"]        = std::move(d_ok);
    bc.diagVar["iface_align"]     = std::move(d_al);
    bc.diagVar["iface_q_eff"]     = std::move(d_qeff);
    bc.diagVar["iface_q_eff_raw"] = std::move(d_qeff_raw);   // 蓄積項 $e_w R_\rho$ を引く前 (旧定義)
}

void checkWallTemperatureSharing(const solverConfig& cfg, const mesh& msh)
{
    // 壁温を陽に扱っている run でのみ検査する (既定 run のログと挙動を変えない)。
    bool anyProfile = false;
    for (const bcond& bc : msh.bconds) {
        const auto it = bc.inputInts.find("wallProfile");
        if (it != bc.inputInts.end() && it->second == 1) anyProfile = true;
    }
    if (!anyProfile && cfg.interfaceDiag != 1) return;

    // DOF -> (bcond index, その bplane の Ts)
    std::unordered_map<geom_int, std::vector<std::pair<int,double>>> owners;
    for (int ib = 0; ib < (int)msh.bconds.size(); ib++) {
        const bcond& bc = msh.bconds[ib];
        if (bc.bcondKind != "wall_isothermal") continue;   // 温度を拘束するのはこの種別
        const auto itv = bc.bvar.find("Ts");
        if (itv == bc.bvar.end()) continue;
        for (size_t k = 0; k < bc.iCells.size() && k < itv->second.size(); k++) {
            owners[bc.iCells[k]].push_back({ib, (double)itv->second[k]});
        }
    }

    int nConflict = 0;
    const double tol = 1.0e-6;   // 同一壁温とみなす差 [K]
    std::string first;
    for (const auto& kv : owners) {
        const auto& lst = kv.second;
        if (lst.size() < 2) continue;
        double tmin = lst[0].second, tmax = lst[0].second;
        for (const auto& e : lst) { tmin = std::min(tmin, e.second); tmax = std::max(tmax, e.second); }
        if (tmax - tmin <= tol) continue;                 // 同じ壁温なら順序に依らない
        // 同一 bcond 内の複数 bplane (node では起きない) は競合としない
        bool crossBcond = false;
        for (size_t a = 0; a + 1 < lst.size(); a++)
            for (size_t b = a + 1; b < lst.size(); b++)
                if (lst[a].first != lst[b].first && std::fabs(lst[a].second - lst[b].second) > tol) crossBcond = true;
        if (!crossBcond) continue;
        if (nConflict == 0) {
            std::ostringstream oss;
            oss << "CV " << kv.first << " shared by";
            for (const auto& e : lst)
                oss << " {physID " << msh.bconds[e.first].physID << " (" << msh.bconds[e.first].physName
                    << ") Ts=" << e.second << "}";
            first = oss.str();
        }
        ++nConflict;
    }

    if (nConflict > 0) {
        std::cerr << "[conjugateWall] ERROR: " << nConflict
                  << " control volume(s) are shared by isothermal walls that prescribe different wall temperatures.\n"
                  << "[conjugateWall]   " << first << "\n"
                  << "[conjugateWall]   node の温度ピンは bcond 順に適用されるため、共有ノードの壁温は\n"
                  << "[conjugateWall]   設定順で決まってしまう (後勝ち)。形状を分割して共有 CV に矛盾する Ts を与えないこと。\n"
                  << "[conjugateWall]   詳細: methods/boundary.md「共役熱伝達 (CHT)」の「起動時に拒否する構成」。\n";
        exit(EXIT_FAILURE);
    }
    std::cout << "[conjugateWall] wall temperature sharing check: OK ("
              << owners.size() << " wall CVs scanned)" << std::endl;
}


// ---------------------------------------------------------------- 内部 CHT (Phase 2a)
namespace {

bool bcondIsConjugate(const bcond& bc)
{
    const auto it = bc.inputInts.find("conjugate");
    return (it != bc.inputInts.end() && it->second == 1);
}

// ---- mode: fem2d の固体状態 (physID ごとに 1 つ持つ) ----
struct SolidState {
    conjugate::SolidMesh mesh;
    std::unique_ptr<conjugate::SolidFem2D> fem;
    std::vector<int>    wallOfIface;   // 界面節点 -> 壁 bplane index (座標一致で 1 対 1)
    std::vector<double> lumped;        // 界面節点の集中辺長 [m]
    std::vector<double> u;             // 固体全節点温度 [K] (これが状態)
    std::vector<double> qIface;        // 直近の界面荷重 [W/m] (出力用。平均後の値)
    std::vector<std::vector<double>> qBuf;  // 界面熱量の環状バッファ (flux_avg 更新分)
    size_t qPos = 0, qFilled = 0;      // 書き込み位置と充填数
    double dTprev = -1.0;
    int    nGrow  = 0;                 // dTw が増え続けた回数 (発散検知)
    bool   ready  = false;
};

std::map<int, SolidState>& solidStates()
{
    static std::map<int, SolidState> s;
    return s;
}

double backResistance(const solverConfig& cfg)
{
    if (cfg.conjugateBackKind == "coolant") return 1.0 / cfg.conjugateHc;
    return 0.0;                                   // isothermal
}

} // namespace

// mode: fem2d の起動時準備 — 固体を読み、流体の壁節点と **座標一致で 1 対 1** に結ぶ。
// 合わないものは黙って内挿せず落とす (plan §4.4d「初版は流体と固体の外周節点を一致させる」)。
void initSolidFem2d(const solverConfig& cfg, const mesh& msh, const bcond& bc)
{
    SolidState& st = solidStates()[bc.physID];
    st.mesh = conjugate::SolidMesh::read(cfg.conjugateSolidFile);
    st.fem.reset(new conjugate::SolidFem2D(st.mesh));
    st.lumped = st.mesh.ifaceLumped();

    const bool nodeMode = (cfg.discretization == "node");
    const int nb = (int)bc.iCells.size();
    const int ni = st.mesh.nIface();

    // **平面 2D ($z\equiv0$) 以外は拒否する** (§4.6a)。押し出し疑似 2D では流体側の量に
    // 奥行きが乗るので、固体の集中辺長 [m] と単位が合わない。
    double zmax = 0.0;
    for (int ib = 0; ib < nb; ib++) {
        double xyz[3];
        dofCoords(msh, nodeMode, bc.iCells[ib], xyz);
        zmax = std::max(zmax, std::fabs(xyz[2]));
    }
    if (zmax > 1.0e-9) {
        std::cerr << "[conjugateWall] ERROR: mode fem2d は平面 2D (z=0) のみ。physID " << bc.physID
                  << " の壁節点に |z| = " << zmax << " m がある。\n"
                  << "[conjugateWall]   押し出し疑似 2D では流体側の熱量に奥行きが乗り、"
                     "固体 (W/m) と単位が合わない。\n";
        exit(EXIT_FAILURE);
    }
    if (nb != ni) {
        std::cerr << "[conjugateWall] ERROR: 壁節点 " << nb << " と固体界面節点 " << ni
                  << " の数が違う (gen_solid_mesh.py --outer-from で合わせること)。\n";
        exit(EXIT_FAILURE);
    }

    st.wallOfIface.assign(ni, -1);
    std::vector<char> used(nb, 0);
    double dworst = 0.0;
    for (int i = 0; i < ni; i++) {
        int best = -1;
        double dbest = 1e30;
        for (int ib = 0; ib < nb; ib++) {
            double xyz[3];
            dofCoords(msh, nodeMode, bc.iCells[ib], xyz);
            const double d = std::hypot(xyz[0] - st.mesh.ifaceX[i], xyz[1] - st.mesh.ifaceY[i]);
            if (d < dbest) { dbest = d; best = ib; }
        }
        if (dbest > 1.0e-7 || best < 0 || used[best]) {
            std::cerr << "[conjugateWall] ERROR: 固体界面節点 " << i << " ("
                      << st.mesh.ifaceX[i] << ", " << st.mesh.ifaceY[i] << ") に対応する壁節点が無い"
                      << " (最近傍 " << dbest << " m" << (best >= 0 && used[best] ? ", 既に使用済み" : "")
                      << ")。**内挿はしない**。\n";
            exit(EXIT_FAILURE);
        }
        used[best] = 1;
        st.wallOfIface[i] = best;
        dworst = std::max(dworst, dbest);
    }

    // 初期の固体温度: 界面は現在の壁温、内部はその平均から始める (最初の更新で解き直す)。
    const auto itTs = bc.bvar.find("Ts");
    double tsum = 0.0;
    for (int i = 0; i < ni; i++)
        tsum += (itTs != bc.bvar.end() && st.wallOfIface[i] < (int)itTs->second.size())
                ? (double)itTs->second[st.wallOfIface[i]] : 300.0;
    st.u.assign(st.mesh.nNodes, tsum / std::max(1, ni));
    for (int i = 0; i < ni; i++)
        if (itTs != bc.bvar.end() && st.wallOfIface[i] < (int)itTs->second.size())
            st.u[st.mesh.ifaceNodes[i]] = (double)itTs->second[st.wallOfIface[i]];
    st.ready = true;

    std::cout << "[conjugateWall] physID " << bc.physID << ": solid fem2d "
              << st.mesh.nNodes << " nodes / " << st.mesh.nTris << " tris / "
              << ni << " iface (最大ずれ " << dworst << " m, 帯幅 " << st.fem->bandwidth()
              << ", iface_sha1 " << st.mesh.ifaceSha1 << ")" << std::endl;
}

bool conjugateActive(const solverConfig& cfg, const mesh& msh)
{
    if (cfg.conjugateEnabled != 1) return false;
    for (const bcond& bc : msh.bconds) if (bcondIsConjugate(bc)) return true;
    return false;
}

void initConjugateWalls(const solverConfig& cfg, const mesh& msh)
{
    if (!conjugateActive(cfg, msh)) {
        // ブロックだけ書いて bcond で有効化していない = 黙って無効になるので警告する
        if (cfg.conjugateEnabled == 1)
            std::cerr << "[conjugateWall] WARNING: 'conjugate:' block is present but no bcond has "
                         "ints: {conjugate: 1} — the coupling is inactive.\n";
        return;
    }
    if (cfg.discretization != "node") {
        std::cerr << "[conjugateWall] ERROR: in-solver CHT is implemented for discretization: node only "
                     "(cell は wallProfile までが対象)。\n";
        exit(EXIT_FAILURE);
    }
    if (cfg.unsteady == 1) {
        std::cerr << "[conjugateWall] ERROR: in-solver CHT does not support dual-time (unsteady: 1). "
                     "過渡では拘束反力が C = D_t(VE) - R^raw になり、本実装の定常仮定が崩れる。\n";
        exit(EXIT_FAILURE);
    }
    if (cfg.conjugateFlux == "q_eff" && cfg.interfaceDiag != 1) {
        std::cerr << "[conjugateWall] ERROR: conjugate flux=q_eff (既定) は保存形の界面熱量 iface_q_eff を使うので "
                     "output: {interfaceDiag: 1} が要る。\n"
                  << "[conjugateWall]   旧実装の抵抗加重平均に戻すなら conjugate: {flux: q_compact} "
                     "(ただし §4.3 の保存形と 1〜2 % ずれる。plan §5.1 #66)。\n";
        exit(EXIT_FAILURE);
    }
    for (const bcond& bc : msh.bconds) {
        if (!bcondIsConjugate(bc)) continue;
        if (bc.bcondKind != "wall_isothermal") {
            std::cerr << "[conjugateWall] ERROR: physID " << bc.physID << " has ints: {conjugate: 1} but kind="
                      << bc.bcondKind << " (must be wall_isothermal).\n";
            exit(EXIT_FAILURE);
        }
        const FirstInterior& fi = firstInterior(cfg, msh, bc);
        if (cfg.conjugateMode == "fem2d") initSolidFem2d(cfg, msh, bc);
        geom_int nbad = 0;
        for (size_t ib = 0; ib < fi.ok.size(); ++ib) if (!fi.ok[ib]) ++nbad;
        if (nbad > 0) {
            std::cerr << "[conjugateWall] ERROR: physID " << bc.physID << ": " << nbad << " / "
                      << fi.ok.size() << " wall CVs have no usable first interior point "
                      << "(alignment < " << cfg.interfaceDiagAlignMin << "). 角・斜交で界面抵抗が定義できない。\n";
            exit(EXIT_FAILURE);
        }
    }
    // 界面ゲートの許容を run に**事前登録**する (判定は tools/check_cht_interface.py)。
    if (cfg.conjugateGateSet == 1) {
        std::ofstream ofs("conjugate_gate.json");
        ofs.precision(10);
        ofs << "{\n  \"eps_rel\": " << cfg.conjugateGateEpsRel
            << ",\n  \"eps_abs_Wm2\": " << cfg.conjugateGateEpsAbs
            << ",\n  \"dT_K\": " << cfg.conjugateGateDtK
            << ",\n  \"n_consec\": " << cfg.conjugateGateNConsec;
        if (cfg.conjugateGateTolSolid >= 0.0) ofs << ",\n  \"tol_solid\": " << cfg.conjugateGateTolSolid;
        ofs << "\n}\n";
        std::cout << "[conjugateWall] G-if の許容を conjugate_gate.json に登録した "
                     "(eps_rel " << cfg.conjugateGateEpsRel << ", eps_abs " << cfg.conjugateGateEpsAbs
                  << " W/m2, dT " << cfg.conjugateGateDtK << " K, n_consec "
                  << cfg.conjugateGateNConsec << ")" << std::endl;
    } else {
        std::cout << "[conjugateWall] WARNING: conjugate.gate が無い。界面の合否は判定できない "
                     "(結果を見てから許容を決めないこと。plan §6 G-if)" << std::endl;
    }

    std::cout << "[conjugateWall] in-solver CHT active (mode=" << cfg.conjugateMode
              << ", flux=" << cfg.conjugateFlux
              << ", interval=" << cfg.conjugateInterval << ", warmup=" << cfg.conjugateWarmup
              << ", relax=" << cfg.conjugateRelax << ")" << std::endl;
}

namespace {

// G-if の履歴 1 行 (local1d と共通の書式にする)。
void appendInterfaceHistory(int iStep, int physID, int n, double tmean, double tmin, double tmax,
                            double dTmax, double resAbs, double resMax, double resRel,
                            double qTotal, double resSolid)
{
    static bool headerDone = false;
    const std::string fn = "conjugate_history.csv";
    const bool needHeader = !headerDone && !std::ifstream(fn).good();
    std::ofstream ofs(fn, std::ios::app);
    if (needHeader)
        ofs << "step,physID,n,Tw_mean,Tw_min,Tw_max,dTw_max,res_abs_Wm2,res_max_W,res_rel,"
               "q_total,res_solid\n";
    headerDone = true;
    ofs.precision(10);
    ofs << iStep << "," << physID << "," << n << ","
        << std::scientific << tmean << "," << tmin << "," << tmax << "," << dTmax << ","
        << resAbs << "," << resMax << "," << resRel << "," << qTotal << "," << resSolid << "\n";
}

} // namespace

// mode: fem2d の 1 回の更新 (plan §4.6a)。
//   (K_s(u^k) + E^T D_f E) u^{k+1} = b_s + E^T [ Q_f + D_f E u^k ]
// Q_f は**固体側の集中辺長** × `iface_q_eff` [W/m2] で作る (流体の surfArea は使わない)。
// D_f は界面対角のみ (g_f A_i)。受理判定・line search・Anderson は持ち込まない (§4.6a)。
void updateFem2dWall(const solverConfig& cfg, const mesh& msh, bcond& bc,
                     const FirstInterior& fi, const std::vector<flow_float>* qeff,
                     const std::vector<flow_float>& T, const std::vector<flow_float>& thermCond,
                     const std::vector<flow_float>& cp, const std::vector<flow_float>& visTurb,
                     int iStep)
{
    SolidState& st = solidStates()[bc.physID];
    if (!st.ready) return;
    auto& Ts = bc.bvar["Ts"];
    const int ni = st.mesh.nIface();

    std::vector<double> Qf(ni, 0.0), Df(ni, 0.0), Twk(ni, 0.0);
    for (int i = 0; i < ni; i++) {
        const int ib = st.wallOfIface[i];
        const geom_int ic = bc.iCells[ib];
        double q = 0.0;
        if (qeff) {
            q = (double)(*qeff)[ib];
            if (!std::isfinite(q)) {
                std::cerr << "[conjugateWall] ERROR: physID " << bc.physID
                          << ": iface_q_eff が界面節点 " << i << " で非有限。\n";
                exit(EXIT_FAILURE);
            }
        } else {
            // flux: q_compact (A/B 用)。壁温は bvar の指定値を使う。
            const geom_int j = fi.jdof[ib];
            const double keff = (double)thermCond[ic] + (double)cp[ic]*(double)visTurb[ic]/cfg.turbulentPrandtl;
            q = (j >= 0) ? keff * ((double)T[j] - (double)Ts[ib]) / fi.d1[ib] : 0.0;
        }
        const double keff = (double)thermCond[ic] + (double)cp[ic]*(double)visTurb[ic]/cfg.turbulentPrandtl;
        Qf[i]  = q * st.lumped[i];                                    // [W/m]
        Df[i]  = keff / fi.d1[ib] * st.lumped[i] * cfg.conjugateDfScale;
        Twk[i] = (double)Ts[ib];
    }

    // ---- 流束の時間平均 (plan §5.1 #70 / §6 V4b(g)) ----
    // 流体側に局所振動があると**瞬時の $Q_f$ では界面ゲートが床に当たる** (C3X 吸込面の k オンセット
    // 前線が周期 1012 step で揺れ、最悪節点の $Q_f$ が中央値の 78 倍ばらつく)。N 更新の後方移動平均を
    // 使う。$D_f$ は平均しない (固定点に効かず、反復経路だけを決めるため)。
    const int navg = std::max(1, cfg.conjugateFluxAvg);
    if (navg > 1) {
        if (st.qBuf.size() != (size_t)navg) { st.qBuf.assign(navg, std::vector<double>(ni, 0.0)); st.qPos = 0; st.qFilled = 0; }
        st.qBuf[st.qPos] = Qf;
        st.qPos = (st.qPos + 1) % (size_t)navg;
        st.qFilled = std::min(st.qFilled + 1, (size_t)navg);
        for (int i = 0; i < ni; i++) {
            double sum = 0.0;
            for (size_t k = 0; k < st.qFilled; k++) sum += st.qBuf[k][i];
            Qf[i] = sum / (double)st.qFilled;
        }
    }

    // ---- G-if: **更新前 (未緩和)** の界面残差 r_i = (K u^k - b)_i - Q_f,i ----
    // (平均を使っているときは**平均後の $Q_f$ で**測る = 実際に解いている方程式の残差にする)
    const std::vector<double> r0 = st.fem->residual(st.u);
    double resAbs = 0.0, resMax = 0.0, qfMax = 0.0, qTotal = 0.0, resSolid = 0.0;
    {
        std::vector<char> isIface(st.mesh.nNodes, 0);
        for (int i = 0; i < ni; i++) isIface[st.mesh.ifaceNodes[i]] = 1;
        for (int i = 0; i < ni; i++) {
            const double r = r0[st.mesh.ifaceNodes[i]] - Qf[i];
            resMax = std::max(resMax, std::fabs(r));
            if (st.lumped[i] > 0.0) resAbs = std::max(resAbs, std::fabs(r) / st.lumped[i]);
            qfMax  = std::max(qfMax, std::fabs(Qf[i]));
            qTotal += Qf[i];
        }
        for (int i = 0; i < st.mesh.nNodes; i++)
            if (!isIface[i]) resSolid = std::max(resSolid, std::fabs(r0[i]));
    }

    // ---- 更新 (固体全節点系をバンド Cholesky で 1 回解く) ----
    st.fem->assemble(st.u, Df);
    st.fem->factorize();
    std::vector<double> rhs = st.fem->rhs();
    for (int i = 0; i < ni; i++) rhs[st.mesh.ifaceNodes[i]] += Qf[i] + Df[i] * Twk[i];
    st.fem->solveInPlace(rhs);

    // ---- 安全装置 (自動調整はしない。止めて報告する。§4.6a) ----
    double tcMin = 1e30;
    for (const double t : st.mesh.robinTc) tcMin = std::min(tcMin, t);
    for (int i = 0; i < st.mesh.nNodes; i++) {
        if (!std::isfinite(rhs[i]) || rhs[i] < tcMin - 20.0) {
            std::cerr << "[conjugateWall] ERROR: step " << iStep << " physID " << bc.physID
                      << ": 固体温度が範囲外 (節点 " << i << " = " << rhs[i] << " K, 下限 "
                      << tcMin - 20.0 << ")。連成が発散している。\n"
                      << "[conjugateWall]   conjugate.Df_scale を上げる / interval を短くする / "
                         "流体を先に落ち着かせること。\n";
            exit(EXIT_FAILURE);
        }
    }

    // ---- 反映 ----
    double dTmax = 0.0, tsum = 0.0, tmin = 1e30, tmax = -1e30;
    for (int i = 0; i < ni; i++) {
        const int ib = st.wallOfIface[i];
        const double Tnew = rhs[st.mesh.ifaceNodes[i]];
        const double Tw = (1.0 - cfg.conjugateRelax) * Twk[i] + cfg.conjugateRelax * Tnew;
        dTmax = std::max(dTmax, std::fabs(Tw - Twk[i]));
        Ts[ib] = (flow_float)Tw;
        rhs[st.mesh.ifaceNodes[i]] = Tw;              // 状態と壁温を一致させる
        tsum += Tw; tmin = std::min(tmin, Tw); tmax = std::max(tmax, Tw);
    }
    st.u = rhs;
    st.qIface.assign(ni, 0.0);
    for (int i = 0; i < ni; i++) st.qIface[i] = Qf[i];

    if (st.dTprev >= 0.0 && dTmax > st.dTprev) st.nGrow++; else st.nGrow = 0;
    if (st.nGrow >= 10 && st.dTprev > 0.0 && dTmax > 2.0 * st.dTprev) {
        std::cerr << "[conjugateWall] ERROR: step " << iStep << " physID " << bc.physID
                  << ": max|dTw| が 10 更新連続で増え、かつ 2 倍を超えた ("
                  << st.dTprev << " -> " << dTmax << " K)。連成が発散している。\n";
        exit(EXIT_FAILURE);
    }
    st.dTprev = dTmax;

    if (cfg.gpu == 1 && bc.bvar_d.count("Ts"))
        cudaWrapper::cudaMemcpy_H2D_wrapper(Ts.data(), bc.bvar_d["Ts"], (geom_int)Ts.size());

    appendInterfaceHistory(iStep, bc.physID, ni, tsum / std::max(1, ni), tmin, tmax, dTmax,
                           resAbs, resMax, (qfMax > 0.0 ? resMax / qfMax : 0.0), qTotal, resSolid);

    if (iStep % (cfg.conjugateInterval * 20) == 0)
        std::cout << "[conjugateWall] step " << iStep << " physID " << bc.physID
                  << (navg > 1 ? " (flux_avg " + std::to_string(st.qFilled) + "/" + std::to_string(navg) + ")" : "")
                  << ": Tw " << tmin << " .. " << tmax << " K, max|dTw| " << dTmax
                  << " K, if-res " << resAbs << " W/m2 (rel "
                  << (qfMax > 0.0 ? resMax / qfMax : 0.0) << "), solid-res " << resSolid
                  << " W/m" << std::endl;
}

void updateConjugateWalls(const solverConfig& cfg, mesh& msh, variables& var, int iStep)
{
    if (!conjugateActive(cfg, msh)) return;
    if (iStep < cfg.conjugateWarmup) return;
    if ((iStep - cfg.conjugateWarmup) % cfg.conjugateInterval != 0) return;

    const std::vector<flow_float> T        = pullField(cfg, var, "T",         msh.nCells);
    const std::vector<flow_float> thermCond= pullField(cfg, var, "thermCond", msh.nCells);
    const std::vector<flow_float> cp       = pullField(cfg, var, "cp",        msh.nCells);
    const std::vector<flow_float> visTurb  = pullField(cfg, var, "vis_turb",  msh.nCells);
    const double Rback = backResistance(cfg);
    const bool useEff  = (cfg.conjugateFlux == "q_eff");

    const bool femMode = (cfg.conjugateMode == "fem2d");

    for (bcond& bc : msh.bconds) {
        if (!bcondIsConjugate(bc)) continue;
        const FirstInterior& fi = firstInterior(cfg, msh, bc);

        // ---- 保存形 $Q_f$ を使う場合は、この step の壁面素材を host に降ろして作る ----
        // `iface_q_eff` = $(R^{raw}-F_w-e_wR_\rho)/A$ (定義と符号は本ファイル冒頭と methods/boundary.md)。
        // 素材 (`ifaceFw`/`ifaceRraw`/`ifaceRro`) はデバイス側で毎 step 採れているので D2H だけで足りる。
        const std::vector<flow_float>* qeff = nullptr;
        if (useEff) {
            bc.copyVariables_bplane_D2H();
            fillInterfaceDiagnostics(cfg, msh, var, bc);
            const auto itq = bc.diagVar.find("iface_q_eff");
            if (itq == bc.diagVar.end() || itq->second.size() < bc.iCells.size()) {
                std::cerr << "[conjugateWall] ERROR: physID " << bc.physID
                          << ": conjugate flux=q_eff だが iface_q_eff が作れない。"
                             "output.interfaceDiag: 1 と定常 (unsteady: 0) の node 等温壁が要る。\n";
                exit(EXIT_FAILURE);
            }
            qeff = &itq->second;
        }

        if (femMode) { updateFem2dWall(cfg, msh, bc, fi, qeff, T, thermCond, cp, visTurb, iStep); continue; }

        auto& Ts = bc.bvar["Ts"];
        double dTmax = 0.0;
        // G-if の素材 (plan boundary-conjugate-heat-transfer §6 G-if)。
        //   r_i = Q_{f,i} - g_s (T_{w,i}-T_b) A_i   [W]  … **更新前 (未緩和) の界面残差**
        double resAbsW = 0.0;      // max_i |r_i| / A_i  [W/m2]  (局所面積で規格化)
        double resMaxW = 0.0;      // max_i |r_i|        [W or W/m]
        double qfMaxW  = 0.0;      // max_i |Q_{f,i}|
        double qTotal  = 0.0;      // sum_i Q_{f,i}
        double Tsum = 0.0, Tmin = 1e30, Tmax = -1e30;
        geom_int nUsed = 0, nBadQ = 0;

        for (size_t ib = 0; ib < bc.iCells.size() && ib < Ts.size(); ++ib) {
            const geom_int ic = bc.iCells[ib];
            const geom_int j  = fi.jdof[ib];
            if (j < 0) continue;
            const double keff = (double)thermCond[ic] + (double)cp[ic]*(double)visTurb[ic]/cfg.turbulentPrandtl;
            const double gf   = keff / fi.d1[ib];                    // [W/m2K] 流体側コンダクタンス
            const double Rtot = cfg.conjugateThickness / cfg.conjugateKsolid + Rback;
            const double gs   = 1.0 / Rtot;
            const double Twk  = (double)Ts[ib];

            double Tnew;
            if (useEff) {
                // §4.2 の更新式 (固定点を保存する形): $(g_s+D_f)T^{k+1}=g_sT_b+q_{\rm eff}+D_fT^k$。
                // 収束すると $g_s(T_w-T_b)=q_{\rm eff}$ = **保存形の界面熱量と固体の 1 次元法則の釣り合い**。
                // $D_f$ は界面抵抗の初期推定 $g_f$ (plan §5.1 #32: `hA` は初期推定に留める)。収束速度だけを決める。
                const double q = (double)(*qeff)[ib];
                if (!std::isfinite(q)) { ++nBadQ; continue; }
                const double Df = gf;
                Tnew = (gs*cfg.conjugateTb + q + Df*Twk) / (gs + Df);

                const geom_int ip = bc.iPlanes[ib];
                const double area = (double)msh.planes[ip].surfArea;
                const double Qf   = q * area;
                const double r    = Qf - gs*(Twk - cfg.conjugateTb)*area;   // 未緩和の界面残差 [W]
                if (area > 0.0) resAbsW = std::max(resAbsW, std::fabs(r)/area);
                resMaxW = std::max(resMaxW, std::fabs(r));
                qfMaxW  = std::max(qfMaxW,  std::fabs(Qf));
                qTotal += Qf;
            } else {
                // 旧実装 (`flux: q_compact`) — 抵抗加重平均。**ビット不変**にするため式をそのまま残す。
                Tnew = (gf * (double)T[j] + gs * cfg.conjugateTb) / (gf + gs);
            }

            const double Tw = (1.0 - cfg.conjugateRelax) * Twk + cfg.conjugateRelax * Tnew;
            dTmax = std::max(dTmax, std::fabs(Tw - Twk));
            Ts[ib] = (flow_float)Tw;
            Tsum += Tw; Tmin = std::min(Tmin, Tw); Tmax = std::max(Tmax, Tw); ++nUsed;
        }

        if (nBadQ > 0) {
            std::cerr << "[conjugateWall] ERROR: physID " << bc.physID << ": iface_q_eff が " << nBadQ
                      << " 節点で非有限。適用範囲外の構成 (dual-time / 周期に属する壁ノード / 断熱壁) では\n"
                      << "[conjugateWall]   保存形の界面熱量が定義できない (conjugateWall.cpp の適用範囲を参照)。\n"
                      << "[conjugateWall]   外部ループ (tools/cht_loop.py) を使うか、flux: q_compact で旧定義に戻すこと。\n";
            exit(EXIT_FAILURE);
        }

        // **注意**: H2D ラッパの引数順は (host, device, n) で、D2H の (device, host, n) と逆。
        if (cfg.gpu == 1 && bc.bvar_d.count("Ts"))
            cudaWrapper::cudaMemcpy_H2D_wrapper(Ts.data(), bc.bvar_d["Ts"], (geom_int)Ts.size());

        // ---- G-if: 更新ごとに界面の残差と温度更新量を残す (判定は check_* 側で行う) ----
        if (useEff && nUsed > 0) {
            static std::set<int> headerWritten;
            const std::string fn = "conjugate_history.csv";
            const bool needHeader = headerWritten.empty() && !std::ifstream(fn).good();
            std::ofstream ofs(fn, std::ios::app);
            if (needHeader)
                ofs << "step,physID,n,Tw_mean,Tw_min,Tw_max,dTw_max,res_abs_Wm2,res_max_W,res_rel,q_total\n";
            headerWritten.insert(bc.physID);
            ofs.precision(10);
            ofs << iStep << "," << bc.physID << "," << nUsed << ","
                << std::scientific << Tsum/nUsed << "," << Tmin << "," << Tmax << ","
                << dTmax << "," << resAbsW << "," << resMaxW << ","
                << (qfMaxW > 0.0 ? resMaxW/qfMaxW : 0.0) << "," << qTotal << "\n";
        }

        if (iStep % (cfg.conjugateInterval * 20) == 0) {
            std::cout << "[conjugateWall] step " << iStep << " physID " << bc.physID
                      << ": Tw " << Tmin << " .. " << Tmax << " K, max|dTw| " << dTmax << " K";
            if (useEff) std::cout << ", if-res " << resAbsW << " W/m2 (rel "
                                  << (qfMaxW > 0.0 ? resMaxW/qfMaxW : 0.0) << ")";
            std::cout << std::endl;
        }
    }
}

void writeConjugateState(const solverConfig& cfg, const mesh& msh, int iStep)
{
    if (!conjugateActive(cfg, msh)) return;
    const bool nodeMode = (cfg.discretization == "node");
    for (const bcond& bc : msh.bconds) {
        if (!bcondIsConjugate(bc)) continue;

        // mode: fem2d は**固体場そのもの**を流体と同じ間隔で出す (plan §4.6a)。
        // 可視化 (`res_solid_<step>`) と再開 (`conjugate_state.h5`) は**別ファイル**にする。
        // **本関数は毎 step 呼ばれる** (壁温 CSV は上書きなので無害だった)。固体場は step ごとに
        // 別ファイルになるので、流体の出力間隔で間引かないと run が数万ファイル・10 GB 級になる
        // (2026-09-23 に 19751 ファイル/13 GB を作った)。
        const bool outStep = (cfg.outStepInterval > 0)
                          && (iStep % cfg.outStepInterval == 0) && (iStep >= cfg.outStepStart);
        if (cfg.conjugateMode == "fem2d" && outStep) {
            const auto it = solidStates().find(bc.physID);
            if (it != solidStates().end() && it->second.ready) {
                SolidState& st = it->second;
                std::ostringstream oss;
                oss << "res_solid_" << bc.physID << "_" << iStep;
                conjugate::writeSolidField(oss.str(), st.mesh, st.u, st.qIface,
                                           (double)iStep);
                // 再開用の状態 (可視化の間引きと混ぜない)。
                try {
                    HighFive::File f("conjugate_state_" + std::to_string(bc.physID) + ".h5",
                                     HighFive::File::Overwrite);
                    f.createDataSet("SOLID/T", st.u);
                    f.createAttribute("step", iStep);
                    f.createAttribute("iface_sha1", st.mesh.ifaceSha1);
                } catch (const std::exception& e) {
                    std::cerr << "[conjugateWall] WARNING: conjugate_state を書けない: "
                              << e.what() << "\n";
                }
            }
        }
        const auto it = bc.bvar.find("Ts");
        if (it == bc.bvar.end()) continue;
        std::ofstream ofs("conjugate_Tw_" + std::to_string(bc.physID) + ".csv");
        ofs << "x y z Ts\n";
        ofs.precision(10);
        for (size_t ib = 0; ib < bc.iCells.size() && ib < it->second.size(); ++ib) {
            double xyz[3]; dofCoords(msh, nodeMode, bc.iCells[ib], xyz);
            ofs << std::scientific << xyz[0] << " " << xyz[1] << " " << xyz[2] << " "
                << (double)it->second[ib] << "\n";
        }
    }
}

} // namespace conjugateWall
