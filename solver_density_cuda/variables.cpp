#include <cstdlib>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <array>
#include <algorithm>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/calcStructualVariables_d.cuh"



variables::variables()
{
    for (const auto& cellValName : cellValNames)
    {
        this->c.emplace(cellValName, std::vector<flow_float>{});
        this->c_d.emplace(cellValName, nullptr);
    }

    for (const auto& planeValName : planeValNames)
    {
        this->p.emplace(planeValName, std::vector<flow_float>{});
        this->p_d.emplace(planeValName, nullptr);
    }
};

// 化学種 (M2): 1 化学種ごとに必要なセル変数名を生成する。
//   roY{s}   : 保存質量分率 ρY_s
//   Y{s}     : 原始質量分率 Y_s = ρY_s/ρ
//   roY{s}N  : RK ステップ始点, roY{s}M : RK ステージ始点
//   res_roY{s}, res_roY{s}_m : 残差 / 4thRunge 累積
//   transport_diag_Y{s} : 輸送 (移流+拡散) ヤコビアン対角 [m³/s]
//   src_jac_Y{s}        : 源項ヤコビアン対角 (M2 では化学反応なしで 0)
static std::list<std::string> speciesCellVarNames(int s)
{
    const std::string i = std::to_string(s);
    return {
        "roY"+i, "Y"+i, "roY"+i+"N", "roY"+i+"M",
        "res_roY"+i, "res_roY"+i+"_m",
        "transport_diag_Y"+i, "src_jac_Y"+i,
        // 緩和整合 scalar-DPLUR (speciesImplicitCoupling==1) の Jacobi 補正 new/old バッファ。
        // 既定経路 (=0) では未使用だが確保コストは僅少。
        "dq_roY"+i, "dq_roY"+i+"_old",
        // 質量分率セル勾配 + Venkat リミタ (speciesFaceReconstruction==1 の face 整合再構成用)。
        "dY"+i+"dx", "dY"+i+"dy", "dY"+i+"dz", "limiter_Y"+i
    };
}

void variables::registerSpecies(int nSpecies, int chemistry)
{
    // 単成分 (<=1) は M1 と同一経路を保つため化学種変数を登録しない。
    if (nSpecies <= 1) {
        this->nSpeciesRegistered = 0;
        return;
    }

    this->nSpeciesRegistered = nSpecies;
    this->speciesVarNames.clear();

    for (int s = 0; s < nSpecies; s++) {
        this->speciesVarNames.push_back("roY"+std::to_string(s));
        for (const auto& name : speciesCellVarNames(s)) {
            // cellValNames へ追加し、host/device マップにも空エントリを作る
            // (allocVariables がこのリストを走査して確保する)。
            this->cellValNames.push_back(name);
            this->c.emplace(name, std::vector<flow_float>{});
            this->c_d.emplace(name, nullptr);
        }
        // 保存質量分率 roY{s} と原始 Y{s} を HDF5 出力対象に加える (roY はリスタートに必須)。
        this->output_cellValNames.push_back("roY"+std::to_string(s));
        this->output_cellValNames.push_back("Y"+std::to_string(s));
    }

    // 有限速度化学 (chemistry.enabled): 反応熱 Q̇ [W/m3] と化学時間 τ_c [s] の診断出力。
    if (chemistry != 0) {
        for (const auto& name : {std::string("chemQdot"), std::string("chemTau")}) {
            this->cellValNames.push_back(name);
            this->c.emplace(name, std::vector<flow_float>{});
            this->c_d.emplace(name, nullptr);
            this->output_cellValNames.push_back(name);
        }
    }

    std::cout << "registerSpecies: nSpecies=" << nSpecies
              << " -> registered " << nSpecies*10 << " cell variables" << (chemistry ? " (+chemistry diagnostics)" : "") << "\n";
}

// 非平衡凝縮 (Phase 1): 1 モーメント (保存量名 consName 例 "rog_0") ごとに必要なセル変数名を生成する。
//   <consName>        : 保存量 ρφ (例 rog_0, roQ0_0)
//   <prim>            : 原始量 φ = ρφ/ρ。consName の先頭 "ro" を外した名前 (g_0, Q0_0)
//   <consName>N / M   : RK ステップ始点 / ステージ始点
//   res_<consName>, res_<consName>_m : 残差 / 4thRunge 累積
//   src_jac_<prim>        : 源項ヤコビアン対角 (Phase 1 は 0)
//   transport_diag_<prim> : 輸送 (移流) ヤコビアン対角 [m³/s]
static std::list<std::string> condMomentCellVarNames(const std::string& consName)
{
    const std::string prim = consName.substr(2); // 先頭 "ro" を除去
    return {
        consName, prim, consName+"N", consName+"M",
        "res_"+consName, "res_"+consName+"_m",
        "src_jac_"+prim, "transport_diag_"+prim
    };
}

void variables::registerCondensation(int nCondSpecies)
{
    if (nCondSpecies <= 0) {
        this->nCondSpeciesRegistered = 0;
        return;
    }

    this->nCondSpeciesRegistered = nCondSpecies;
    this->condMomentConsNames.clear();

    // 1 凝縮種あたり 4 モーメント。順序は論文の保存ベクトル後半 (ρg,ρQ2,ρQ1,ρQ0)。
    const std::array<std::string, 4> bases = {"g", "Q2", "Q1", "Q0"};

    for (int s = 0; s < nCondSpecies; s++) {
        for (const auto& b : bases) {
            const std::string consName = "ro" + b + "_" + std::to_string(s);
            this->condMomentConsNames.push_back(consName);
            for (const auto& name : condMomentCellVarNames(consName)) {
                this->cellValNames.push_back(name);
                this->c.emplace(name, std::vector<flow_float>{});
                this->c_d.emplace(name, nullptr);
            }
            // 原始量・保存量を HDF5 出力対象に加える (可視化)。
            this->output_cellValNames.push_back(consName);
            this->output_cellValNames.push_back(consName.substr(2));
        }
        // 診断 (source kernel が毎ステップ書く): 過飽和 S=p_v/p_sat, 成長率 dr/dt [m/s] (負=蒸発),
        // 体積平均半径 r30 [m] (蒸発分岐で評価; 0=未評価)。確保時に 0 初期化 (prefix "cond")。
        // condTsat_: 飽和温度 T_sat(p_v) [K] (過冷却度 T_sat-T の評価用, plans/accepted/condensation-equilibrium.md)
        for (const auto& d : {std::string("condS_"), std::string("condDrdt_"), std::string("condR30_"), std::string("condTsat_")}) {
            const std::string name = d + std::to_string(s);
            this->cellValNames.push_back(name);
            this->c.emplace(name, std::vector<flow_float>{});
            this->c_d.emplace(name, nullptr);
            this->output_cellValNames.push_back(name);
        }
    }

    std::cout << "registerCondensation: nCondSpecies=" << nCondSpecies
              << " -> registered " << nCondSpecies*4*8 << " cell variables\n";
}

variables::~variables() {
    for (auto& cellValName : cellValNames)
    {
        cudaWrapper::cudaFree_wrapper(this->c_d.at(cellValName));
    }

    for (auto& planeValName : planeValNames)
    {
        cudaWrapper::cudaFree_wrapper(this->p_d.at(planeValName));
    }
}

void variables::allocVariables(const int &useGPU , mesh& msh)
{
    // W-I 実力診断 (§4.2) は FORGE_WI_FORCE_DIAG=1 のときだけ有効。OFF なら
    // 確保も出力もしない (通常経路のメモリ・D2H・HDF5 を増やさない)。
    {
        const char* e = std::getenv("FORGE_WI_FORCE_DIAG");
        if (!(e && std::atoi(e) != 0)) {
            for (const char* n : {"wi_ftan", "wi_fnrm", "wi_fnrm_abs", "wi_ftan_res"}) {
                cellValNames.remove(n);
                output_cellValNames.remove(n);
                c.erase(n);
                c_d.erase(n);
            }
        }
    }
    // E3 (§5) の wf_sprod も env ゲート (OFF ならメモリ・D2H・HDF5 を増やさない)
    {
        const char* e = std::getenv("FORGE_WF_OMEGA_SOURCE");
        if (!(e && std::atoi(e) != 0)) {
            cellValNames.remove("wf_sprod");
            output_cellValNames.remove("wf_sprod");
            c.erase("wf_sprod");
            c_d.erase("wf_sprod");
        }
    }
    // 閉包則診断 (§5.2 ③) の wf_g も env ゲート
    {
        const char* e = std::getenv("FORGE_WF_CLOSURE_DIAG");
        if (!(e && std::atoi(e) != 0)) {
            cellValNames.remove("wf_g");
            output_cellValNames.remove("wf_g");
            c.erase("wf_g");
            c_d.erase("wf_g");
        }
    }
    // 代表点幾何診断 (§3.1) も env ゲート
    {
        const char* e = std::getenv("FORGE_WF_REP_DIAG");
        if (!(e && std::atoi(e) != 0)) {
            for (const char* n : {"rep_id", "rep_y", "rep_dist", "rep_cos", "rep_toff", "rep_wdratio",
                                  "rep_nx", "rep_ny", "rep_nz"}) {
                cellValNames.remove(n);
                output_cellValNames.remove(n);
                c.erase(n);
                c_d.erase(n);
            }
        }
    }
    // omega 項別収支 (§4.1) も env ゲート
    {
        const char* e = std::getenv("FORGE_OMEGA_BUDGET");
        if (!(e && std::atoi(e) != 0)) {
            for (const char* n : {"omg_prod", "omg_dest", "omg_cross", "omg_trans", "omg_axisym"}) {
                cellValNames.remove(n);
                output_cellValNames.remove(n);
                c.erase(n);
                c_d.erase(n);
            }
        }
    }
    for (auto& cellValName : cellValNames)
    {
        //this->c[cellValName].resize(msh.nCells);
        std::vector<flow_float>& cellValues = this->c.at(cellValName);
        cellValues.resize(msh.nCells_all); // including ghost cells

        if (useGPU == 1) 
        {
            gpuErrchk( cudaMalloc((void**) &(this->c_d.at(cellValName)), (msh.nCells_all)*sizeof(flow_float)) );
            // 診断配列は毎ステップ書かれるとは限らないので確保直後に 0 初期化する
            // (未初期化 device メモリを出力しないため)。
            if (cellValName.rfind("wi_", 0) == 0 || cellValName.rfind("omg_", 0) == 0
                || cellValName.rfind("rep_", 0) == 0 || cellValName.rfind("cond", 0) == 0
                || cellValName == "wf_irep_flag" || cellValName == "wf_sprod"
                || cellValName == "wf_g") {
                gpuErrchk( cudaMemset(this->c_d.at(cellValName), 0, (msh.nCells_all)*sizeof(flow_float)) );
            }
        }

    }
    for (auto& planeValName : planeValNames)
    {
        std::vector<flow_float>& planeValues = this->p.at(planeValName);
        planeValues.resize(msh.nPlanes);

        if (useGPU == 1) 
        {
            gpuErrchk( cudaMalloc((void**) &(this->p_d.at(planeValName)), msh.nPlanes*sizeof(flow_float)) );
        }
    }
}



void variables::copyVariables_cell_plane_H2D_all()
{
    for (auto& name : this->cellValNames)
    {
        std::vector<flow_float>& cellValues = this->c.at(name);
        cudaWrapper::cudaMemcpy_H2D_wrapper(cellValues.data() , this->c_d.at(name), cellValues.size());
    }
    for (auto& name : this->planeValNames)
    {
        std::vector<flow_float>& planeValues = this->p.at(name);
        cudaWrapper::cudaMemcpy_H2D_wrapper(planeValues.data() , this->p_d.at(name), planeValues.size());
    }
}

//void variables::copyVariables_cell_H2D(std::string name)
void variables::copyVariables_cell_H2D(std::list<std::string> names)
{
    for (auto& name : names) {
        std::vector<flow_float>& cellValues = this->c.at(name);
        cudaWrapper::cudaMemcpy_H2D_wrapper(cellValues.data() , this->c_d.at(name), cellValues.size());
    }
}
void variables::copyVariables_plane_H2D(std::list<std::string> names)
{
    for (auto& name : names) {
        std::vector<flow_float>& planeValues = this->p.at(name);
        cudaWrapper::cudaMemcpy_H2D_wrapper(planeValues.data() , this->p_d.at(name), planeValues.size());
    }
}

void variables::copyVariables_cell_plane_D2H_all()
{
    for (auto& name : this->cellValNames)
    {
        std::vector<flow_float>& cellValues = this->c.at(name);
        cudaWrapper::cudaMemcpy_D2H_wrapper(this->c_d.at(name), cellValues.data() , cellValues.size());
    }
    for (auto& name : this->planeValNames)
    {
        std::vector<flow_float>& planeValues = this->p.at(name);
        cudaWrapper::cudaMemcpy_D2H_wrapper(this->p_d.at(name), planeValues.data() , planeValues.size());
    }
}

void variables::copyVariables_cell_D2H(std::list<std::string> names)
{
    for (auto& name : names) {
        std::vector<flow_float>& cellValues = this->c.at(name);
        cudaWrapper::cudaMemcpy_D2H_wrapper(this->c_d.at(name), cellValues.data(), cellValues.size());
    }
}

void variables::copyVariables_plane_D2H(std::list<std::string> names)
{
    for (auto& name : names) {
        std::vector<flow_float>& planeValues = this->p.at(name);
        cudaWrapper::cudaMemcpy_D2H_wrapper(this->p_d.at(name), planeValues.data(), planeValues.size());
    }
}

void variables::setStructuralVariables(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh)
{
    if (cfg.gpu==1) {
        variables::setStructuralVariables_d(cfg, cuda_cfg, msh);
        return;
    }

    //geom_float ss;
    std::vector<geom_float> sv(3);

    std::vector<geom_float> dccv(3);
    geom_float dcc;

    std::vector<geom_float> dc1pv(3);
    std::vector<geom_float> dc2pv(3);
    geom_float dc2p;

    std::vector<geom_float> pcent(3);
    std::vector<geom_float> c1cent(3);
    std::vector<geom_float> c2cent(3);

    geom_int ic1;
    geom_int ic2;
    std::vector<flow_float>& fxp  = this->p.at("fx");
    std::vector<flow_float>& dccp = this->p.at("dcc");
    std::vector<flow_float>& vvol = this->c.at("volume");
    std::vector<flow_float>& vccx = this->c.at("ccx");
    std::vector<flow_float>& vccy = this->c.at("ccy");
    std::vector<flow_float>& vccz = this->c.at("ccz");

    for (geom_int ip=0 ; ip<msh.nPlanes ; ip++)
    {
        ic1     = msh.planes[ip].iCells[0];
        ic2     = msh.planes[ip].iCells[1];
        sv      = msh.planes[ip].surfVect;

        pcent   = msh.planes[ip].centCoords;

        c1cent  = msh.cells[ic1].centCoords;
        c2cent  = msh.cells[ic2].centCoords;

        dccv[0] = c2cent[0] - c1cent[0];
        dccv[1] = c2cent[1] - c1cent[1];
        dccv[2] = c2cent[2] - c1cent[2];
        dcc     = sqrt( pow(dccv[0], 2.0) + pow(dccv[1], 2.0) + pow(dccv[2], 2.0));

        dc2pv[0] = pcent[0] - c2cent[0];
        dc2pv[1] = pcent[1] - c2cent[1];
        dc2pv[2] = pcent[2] - c2cent[2];
        dc2p     = sqrt( pow(dc2pv[0], 2.0) + pow(dc2pv[1], 2.0) + pow(dc2pv[2], 2.0));

        fxp[ip]  = dc2p/dcc;
        dccp[ip] = dcc;
    }

    // cell
    for (geom_int ic=0 ; ic<msh.nCells ; ic++)
    {
        c1cent   = msh.cells[ic].centCoords;
        vvol[ic] = msh.cells[ic].volume;
        vccx[ic] = c1cent[0];
        vccy[ic] = c1cent[1];
        vccz[ic] = c1cent[2];
    }

}

void variables::setStructuralVariables_d(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh )
{
    geom_float* sx;
    geom_float* sy;
    geom_float* sz;
    geom_float* ss;
    geom_float* sx_planar;
    geom_float* sy_planar;
    geom_float* sz_planar;
    geom_float* ss_planar;
    geom_float* pcx;
    geom_float* pcy;
    geom_float* pcz;
    geom_float* ccx;
    geom_float* ccy;
    geom_float* ccz;
    geom_float* volume;

    sx = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    sy = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    sz = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    ss = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    sx_planar = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    sy_planar = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    sz_planar = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    ss_planar = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    pcx = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    pcy = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);
    pcz = (geom_float*)malloc(sizeof(geom_float)*msh.nPlanes);

    ccx = (geom_float*)malloc(sizeof(geom_float)*(msh.nCells_all));
    ccy = (geom_float*)malloc(sizeof(geom_float)*(msh.nCells_all));
    ccz = (geom_float*)malloc(sizeof(geom_float)*(msh.nCells_all));
    volume = (geom_float*)malloc(sizeof(geom_float)*msh.nCells_all);
    geom_float* A_planar_h = (geom_float*)malloc(sizeof(geom_float)*msh.nCells_all);


    for (geom_int ip=0; ip<msh.nPlanes; ip++)
    {
        sx[ip] = msh.planes[ip].surfVect[0];
        sy[ip] = msh.planes[ip].surfVect[1];
        sz[ip] = msh.planes[ip].surfVect[2];
        ss[ip] = msh.planes[ip].surfArea;
        sx_planar[ip] = sx[ip];
        sy_planar[ip] = sy[ip];
        sz_planar[ip] = sz[ip];
        ss_planar[ip] = ss[ip];
        pcx[ip] = msh.planes[ip].centCoords[0];
        pcy[ip] = msh.planes[ip].centCoords[1];
        pcz[ip] = msh.planes[ip].centCoords[2];
    }

    for (geom_int ic=0; ic<msh.nCells_all; ic++)
    {
        ccx[ic] = msh.cells[ic].centCoords[0];
        ccy[ic] = msh.cells[ic].centCoords[1];
        ccz[ic] = msh.cells[ic].centCoords[2];
        volume[ic] = msh.cells[ic].volume;
        A_planar_h[ic] = 0.0;
    }

    if (cfg.isAxisymmetric == 1 && cfg.axisymMethod == 0) {
        // B 流儀: 幾何量に r 重み付け、半径方向の圧力ソース用に planar 面積を保存。
        // 軸 (r=0) 上の face で S を厳密に 0 にすると、下流の flux/BC カーネルで
        // n = S/|S| = 0/0 = NaN になるため、極小の r フロアを入れて n の方向を
        // 保ちつつ寄与を実質ゼロにする。
        // axisRFloor > 0 (SU2 y<EPS ガードの r 重み版): 床を物理値に引き上げ、軸帯の
        // 面積・体積を消さない (軸半 CV の真空化対策)。hoop ソース/Jacobian/uy_over_r は
        // 別途 ccy < axisRFloor で skip する (ソース・ヤコビアンも入れない)。
        const geom_float r_floor = (cfg.axisRFloor > (flow_float)0.0)
            ? (geom_float)cfg.axisRFloor : (geom_float)1.0e-20;
        for (geom_int ip=0; ip<msh.nPlanes; ip++) {
            const geom_float r_face = (pcy[ip] > r_floor) ? pcy[ip] : r_floor;
            sx[ip] *= r_face;
            sy[ip] *= r_face;
            sz[ip] *= r_face;
            ss[ip] *= r_face;
        }
        // nodeValueAtNode: 実 CV (ic<nCells) の回転半径は双対重心 r̄ (mesh::rEff)。ccy はノード座標 (軸で 0)。
        const bool useREff = (msh.nodeValueAtNode == 1 && (geom_int)msh.rEff.size() == msh.nCells);
        // ゴースト CV も所有 CV の r̄ を使う。値位置=ノードでは境界ノードが境界面上に乗り鏡映距離が 0 =
        // ゴーストがノードと同位置になるため、軸∩境界コーナー (r=0) でゴーストの回転体積が r 床 (1e-20) に
        // 潰れ、setDT の dx=vol_ghost/|S| が ~1e-20 → 局所 CFL ~1e13 → dt_local ~1e-22 でその CV が
        // 完全に凍結する (case/43 run_0003 の入口軸/出口軸ノードで実測)。
        std::vector<geom_float> rEffGhost;
        if (useREff) {
            rEffGhost.assign(msh.nCells_all, (geom_float)-1.0);
            for (const bcond& bc : msh.bconds)
                for (size_t k = 0; k < bc.iCells_ghst.size() && k < bc.iCells.size(); ++k) {
                    const geom_int ig = bc.iCells_ghst[k], io = bc.iCells[k];
                    if (ig >= msh.nCells && ig < msh.nCells_all && io >= 0 && io < msh.nCells)
                        rEffGhost[ig] = msh.rEff[io];
                }
        }
        for (geom_int ic=0; ic<msh.nCells_all; ic++) {
            const geom_float r_src = useREff ? ((ic < msh.nCells) ? msh.rEff[ic]
                                                : (rEffGhost[ic] >= (geom_float)0.0 ? rEffGhost[ic] : ccy[ic]))
                                             : ccy[ic];
            const geom_float r_cell = (r_src > r_floor) ? r_src : r_floor;
            A_planar_h[ic] = volume[ic];
            volume[ic]     = volume[ic] * r_cell;
        }
        // axisRFloor>0: 床適用後の面ベクトルで各セルの離散閉性 Σ_f S_f (outward) を計算し、
        // A_planar を y 成分 (=一般化 hoop 面積)、A_closure_x を x 成分に置換する。
        // 床なし領域では y 成分は解析 A_planar に一致し x 成分は 0 (厳密幾何)。
        // 全面が床の CV では両成分とも 0 → ソース・Jacobian が自然に消える。
        if (cfg.axisRFloor > (flow_float)0.0 || cfg.hoopAreaFromClosure == 1) {
            // 面の向き規約: planes[ip].iCells[0] にとって外向き (+S)、iCells[1] にとって内向き (-S)。
            // (solver 側 mesh は iPlanesDir を持たないため plane 走査で集計する。)
            std::vector<geom_float> aclx(msh.nCells_all, 0.0), acly(msh.nCells_all, 0.0);
            for (geom_int ip = 0; ip < msh.nPlanes; ++ip) {
                const auto& pc = msh.planes[ip].iCells;
                if (pc.empty()) continue;
                const geom_int ic0 = pc[0];
                if (ic0 >= 0 && ic0 < msh.nCells) { aclx[ic0] += sx[ip]; acly[ic0] += sy[ip]; }
                if (pc.size() > 1) {
                    const geom_int ic1 = pc[1];
                    if (ic1 >= 0 && ic1 < msh.nCells) { aclx[ic1] -= sx[ip]; acly[ic1] -= sy[ip]; }
                }
            }
            // A_planar (勾配分母) は不変のまま、閉性面積は専用配列へ (hoop ソース/Jacobian が参照)。
            cudaMemcpy(this->c_d.at("A_closure_x"), aclx.data(), msh.nCells_all*sizeof(geom_float), cudaMemcpyHostToDevice);
            cudaMemcpy(this->c_d.at("A_closure_y"), acly.data(), msh.nCells_all*sizeof(geom_float), cudaMemcpyHostToDevice);
        }
    } else if (cfg.isAxisymmetric == 1 && cfg.axisymMethod == 1) {
        // SU2 流: 幾何は planar のまま (r 重みなし)。1/y はソース項側 (axisymmetricSourceSU2) で扱う。
        // A_planar は「planar 体積」の意味で保持 (勾配計算の grad_volume がこれを参照するため)。
        for (geom_int ic=0; ic<msh.nCells_all; ic++) {
            A_planar_h[ic] = volume[ic];
        }
    }

    flow_float* p_sx = this->p_d.at("sx");
    flow_float* p_sy = this->p_d.at("sy");
    flow_float* p_sz = this->p_d.at("sz");
    flow_float* p_ss = this->p_d.at("ss");
    flow_float* p_sx_planar = this->p_d.at("sx_planar");
    flow_float* p_sy_planar = this->p_d.at("sy_planar");
    flow_float* p_sz_planar = this->p_d.at("sz_planar");
    flow_float* p_ss_planar = this->p_d.at("ss_planar");
    flow_float* p_pcx = this->p_d.at("pcx");
    flow_float* p_pcy = this->p_d.at("pcy");
    flow_float* p_pcz = this->p_d.at("pcz");
    flow_float* c_ccx = this->c_d.at("ccx");
    flow_float* c_ccy = this->c_d.at("ccy");
    flow_float* c_ccz = this->c_d.at("ccz");
    flow_float* c_volume = this->c_d.at("volume");

    cudaMemcpy(p_sx , sx , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_sy , sy , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_sz , sz , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_ss , ss , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_sx_planar , sx_planar , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_sy_planar , sy_planar , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_sz_planar , sz_planar , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_ss_planar , ss_planar , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);

    cudaMemcpy(p_pcx , pcx , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_pcy , pcy , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(p_pcz , pcz , msh.nPlanes*sizeof(geom_float) , cudaMemcpyHostToDevice);

    cudaMemcpy(c_ccx , ccx , msh.nCells_all*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(c_ccy , ccy , msh.nCells_all*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(c_ccz , ccz , msh.nCells_all*sizeof(geom_float) , cudaMemcpyHostToDevice);
    cudaMemcpy(c_volume , volume , msh.nCells_all*sizeof(geom_float) , cudaMemcpyHostToDevice);

    if (cfg.isAxisymmetric == 1) {
        flow_float* c_A_planar = this->c_d.at("A_planar");
        cudaMemcpy(c_A_planar , A_planar_h , msh.nCells_all*sizeof(geom_float) , cudaMemcpyHostToDevice);
    }

    // SST-DES グリッドスケール Δmax (methods/turbulence §8, plan §4.3)。
    // 各 CV について隣接面を介した重心間距離の最大を取る。BL の高アスペクト比セルで V^{1/3} が
    // 接線スケールを過小評価し DES リミッタが誤発火するのを避けるため、Δmax を採用する。
    // 実行時の重心 (ccx/ccy/ccz = CV 中心) と面接続のみから計算するので、node-centered (median-dual)
    // でも双対 CV の Δ が自動的に得られる (plan §5.6: primal mesh の volume を直接参照しない)。
    // 静的量ゆえ幾何セットアップ時に host で 1 回計算し H2D 転送する。
    {
        geom_float* delta_les_h = (geom_float*)malloc(sizeof(geom_float)*msh.nCells_all);
        for (geom_int ic=0; ic<msh.nCells_all; ic++) delta_les_h[ic] = 0.0;
        for (geom_int ip=0; ip<msh.nPlanes; ip++) {
            const geom_int ic1 = msh.planes[ip].iCells[0];
            const geom_int ic2 = msh.planes[ip].iCells[1];
            const geom_float dx = ccx[ic1] - ccx[ic2];
            const geom_float dy = ccy[ic1] - ccy[ic2];
            const geom_float dz = ccz[ic1] - ccz[ic2];
            const geom_float d  = sqrt(dx*dx + dy*dy + dz*dz);
            if (ic1 < msh.nCells && d > delta_les_h[ic1]) delta_les_h[ic1] = d;
            if (ic2 < msh.nCells && d > delta_les_h[ic2]) delta_les_h[ic2] = d;
        }
        cudaMemcpy(this->c_d.at("delta_les"), delta_les_h, msh.nCells_all*sizeof(geom_float), cudaMemcpyHostToDevice);
        free(delta_les_h);
    }

    calcStructualVariables_d_wrapper(cfg , cuda_cfg , msh , *this);

    free(sx) ; free(sy) ; free(sz) ; free(ss);
    free(sx_planar) ; free(sy_planar) ; free(sz_planar) ; free(ss_planar);
    free(pcx); free(pcy); free(pcz);
    free(ccx); free(ccy); free(ccz);
    free(volume);
    free(A_planar_h);

}

void variables::readValueHDF5(std::string fname , mesh& msh,
                              flow_float kInit, flow_float omegaInit)
{
    HighFive::File file(fname, HighFive::File::ReadOnly);

    // read basic 
    HighFive::Group group = file.getGroup("/VALUE");

    // nodes
    std::vector<geom_float> ro;
    file.getDataSet("/VALUE/ro").read(ro);
    std::vector<geom_float> roUx;
    file.getDataSet("/VALUE/roUx").read(roUx);
    std::vector<geom_float> roUy;
    file.getDataSet("/VALUE/roUy").read(roUy);
    std::vector<geom_float> roUz;
    file.getDataSet("/VALUE/roUz").read(roUz);
    std::vector<geom_float> roe;
    file.getDataSet("/VALUE/roe").read(roe);
    std::vector<geom_float> wall_dist;
    file.getDataSet("/VALUE/wall_dist").read(wall_dist);
    std::vector<geom_float> roK;
    std::vector<geom_float> roOmega;
    const bool has_roK = file.exist("/VALUE/roK");
    const bool has_roOmega = file.exist("/VALUE/roOmega");
    if (has_roK) {
        file.getDataSet("/VALUE/roK").read(roK);
    }
    if (has_roOmega) {
        file.getDataSet("/VALUE/roOmega").read(roOmega);
    }
 
   
    std::vector<flow_float>& v_ro = this->c.at("ro");
    std::vector<flow_float>& v_roUx = this->c.at("roUx");
    std::vector<flow_float>& v_roUy = this->c.at("roUy");
    std::vector<flow_float>& v_roUz = this->c.at("roUz");
    std::vector<flow_float>& v_roe = this->c.at("roe");
    std::vector<flow_float>& v_wall_dist = this->c.at("wall_dist");
    std::vector<flow_float>& v_roK = this->c.at("roK");
    std::vector<flow_float>& v_roOmega = this->c.at("roOmega");

    for (geom_int i=0; i<msh.nCells; i++)
    {
        v_ro[i] = ro[i];
        v_roUx[i] = roUx[i];
        v_roUy[i] = roUy[i];
        v_roUz[i] = roUz[i];
        v_roe[i] = roe[i];
        v_wall_dist[i] = wall_dist[i];
        // IC に roK/roOmega が無い場合は freestream 初期値 ro*kInit / ro*omegaInit を使う。
        // kInit=omegaInit=0 (既定) では従来どおり 0 (ビット不変)。ただし ω=0 は mu_t=k/ω が
        // ill-posed で SST cold start から発散しやすい (特に node wt=1) ため、非 SST IC から
        // SST を始める場合は config で freestream 値を設定すること (variables.hpp / solverConfig)。
        v_roK[i] = has_roK ? roK[i] : ro[i] * kInit;
        v_roOmega[i] = has_roOmega ? roOmega[i] : ro[i] * omegaInit;
    }

    // 壁なしメッシュのガード: 壁 bcond の無いメッシュは converter が wall_dist≡0 のまま出力する。
    // 0 のままだと WALE の Ls=min(κd, CwΔ) が 0 に潰れ vis_turb≡0 (SGS 不活性) になるため、
    // 「全 CV で ≤0 = 壁情報なし」のときは実質∞ (1e30) で充填する (plans/active/turbulence-wale-fix.md)。
    // 壁ありメッシュは壁面上の CV (d=0) があっても max>0 なので不変。
    {
        flow_float wd_max = 0.0;
        for (geom_int i=0; i<msh.nCells; i++) wd_max = std::max(wd_max, v_wall_dist[i]);
        if (wd_max <= 0.0) {
            for (geom_int i=0; i<msh.nCells; i++) v_wall_dist[i] = 1.0e30;
            std::cout << "[variables] wall_dist is zero everywhere (no wall bcond): "
                         "filled with 1e30 so WALE/DES length scales use Cw*Delta.\n";
        }
    }

    std::list<std::string> names = {"ro", "roUx", "roUy", "roUz", "roe", "wall_dist", "roK", "roOmega"};
    this->copyVariables_cell_H2D(names);

    // --- 化学種 (M2): 保存質量分率 ρY_s を読み込む ---
    // 優先順: (1) VALUE/roY{s} (保存量) → (2) VALUE/Y{s} (原始量) から roY=ro*Y を復元
    // → (3) デフォルト: s==0 を ρ (Y0=1)、他を 0。
    // 旧出力 (roY 未保存・Y のみ) のリスタートも (2) で吸収する。
    if (this->nSpeciesRegistered >= 2) {
        std::list<std::string> sp_names;
        for (int s = 0; s < this->nSpeciesRegistered; s++) {
            const std::string si = std::to_string(s);
            const std::string roYname = "roY"+si;
            const std::string Yname    = "Y"+si;
            std::vector<flow_float>& v_roY = this->c.at(roYname);
            std::vector<flow_float>& v_Y   = this->c.at(Yname);

            std::vector<geom_float> roY_in;
            bool has_roY = file.exist("/VALUE/"+roYname);
            bool has_Y   = file.exist("/VALUE/"+Yname);
            if (has_roY) {
                file.getDataSet("/VALUE/"+roYname).read(roY_in);
            } else if (has_Y) {
                // 旧ファイル互換: Y (原始) → roY = ro * Y
                std::vector<geom_float> Y_in;
                file.getDataSet("/VALUE/"+Yname).read(Y_in);
                roY_in.resize(Y_in.size());
                for (std::size_t k = 0; k < Y_in.size(); k++) {
                    roY_in[k] = static_cast<geom_float>(this->c.at("ro")[k]) * Y_in[k];
                }
                has_roY = true;  // フォールバック成功
            }

            for (geom_int i=0; i<msh.nCells; i++) {
                const flow_float roi = this->c.at("ro")[i];
                flow_float roYi;
                if (has_roY)        roYi = roY_in[i];
                else if (s == 0)    roYi = roi;        // 既定: 第 1 化学種のみ
                else                roYi = 0.0;
                v_roY[i] = roYi;
                v_Y[i]   = roYi / std::max(roi, static_cast<flow_float>(1.0e-30));
            }
            sp_names.push_back(roYname);
            sp_names.push_back(Yname);
        }
        this->copyVariables_cell_H2D(sp_names);
    }

    // --- 非平衡凝縮 (Phase 1): 液相モーメント ρφ を読み込む ---
    // 入力 HDF5 に VALUE/<consName> があれば読む。無ければ 0 (dry リスタート)。
    // 原始量 φ=ρφ/ρ も同時に設定する。
    if (this->nCondSpeciesRegistered >= 1) {
        std::list<std::string> cond_names;
        for (const auto& consName : this->condMomentConsNames) {
            const std::string primName = consName.substr(2);
            std::vector<flow_float>& v_cons = this->c.at(consName);
            std::vector<flow_float>& v_prim = this->c.at(primName);

            std::vector<geom_float> cons_in;
            const bool has = file.exist("/VALUE/"+consName);
            if (has) file.getDataSet("/VALUE/"+consName).read(cons_in);

            for (geom_int i=0; i<msh.nCells; i++) {
                const flow_float roi = this->c.at("ro")[i];
                const flow_float consi = has ? cons_in[i] : static_cast<flow_float>(0.0);
                v_cons[i] = consi;
                v_prim[i] = consi / std::max(roi, static_cast<flow_float>(1.0e-30));
            }
            cond_names.push_back(consName);
            cond_names.push_back(primName);
        }
        this->copyVariables_cell_H2D(cond_names);
    }
}