#include <cstdlib>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <array>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"
#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/calcStructualVariables_d.cuh"
#include "cuda_forge/qAccumulator.hpp"



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
        "dY"+i+"dx", "dY"+i+"dy", "dY"+i+"dz", "limiter_Y"+i,
        // dual-time の物理時間レベル Q^n, Q^{n-1} (BDF 項用; N/M は擬似時間の始点。plan species-passive-scalar-unification §4.4)
        "roY"+i+"P", "roY"+i+"PP"
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
        for (const char* c : {"x", "y", "z"}) this->extraOnly_cellValNames.push_back("dY"+std::to_string(s)+"d"+c);   // extraFields 専用
    }

    // 診断 (FORGE_SPECIES_RAW_DIAG=1): 再正規化前の生更新値 roYraw{s} を出力する (plan species-passive-scalar-unification §6-1 の
    // 制御試験: 受動トレーサは ΣρY=ρ 再正規化を受けないので比較対象は生更新値)。既定では登録しない。
    {
        const char* e = std::getenv("FORGE_SPECIES_RAW_DIAG");
        if (e && std::atoi(e) != 0) {
            for (int s = 0; s < nSpecies; s++) {
                const std::string name = "roYraw"+std::to_string(s);
                this->cellValNames.push_back(name);
                this->c.emplace(name, std::vector<flow_float>{});
                this->c_d.emplace(name, nullptr);
                this->output_cellValNames.push_back(name);
            }
            std::cout << "registerSpecies: FORGE_SPECIES_RAW_DIAG=1 -> roYraw{s} (pre-renormalization update) registered\n";
        }
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

// 受動トレーサ (排気率 ξ; physProp.tracer: exhaust)。凝縮モーメントと同じ 8 本構成。
//   roXi : 保存量 ρξ, Xi : 原始量 ξ=ρξ/ρ, roXiN/roXiM : RK ステップ/ステージ始点,
//   res_roXi, res_roXi_m : 残差 / 4thRunge 累積, src_jac_Xi : 源項ヤコビアン (0), transport_diag_Xi : 輸送対角 [m³/s]
void variables::registerTransition(int enabled, int diag)
{
    this->transitionRegistered = (enabled != 0) ? 1 : 0;
    if (enabled == 0) return;
    std::vector<std::string> names = {"roGamma", "roReth", "gammaTr", "reTheta", "gammaEff", "res_roGamma", "res_roReth",
                                      "src_jac_gamma", "src_jac_reth", "transport_diag_gamma", "transport_diag_reth"};
    std::vector<std::string> outs  = {"roGamma", "roReth", "gammaTr", "reTheta", "gammaEff"};
    if (diag != 0) {
        for (const char* nm : {"lmFonset", "lmFlength", "lmFtheta", "lmRethCorr", "lmGammaSep", "lmPgamma", "lmEgamma", "lmPtheta", "lmCorrIter"}) {
            names.emplace_back(nm); outs.emplace_back(nm);
        }
        outs.emplace_back("src_jac_gamma"); outs.emplace_back("src_jac_reth");   // 陰的対角も単体検査の対象 (codex result M3)
    }
    for (const auto& name : names) {
        this->cellValNames.push_back(name);
        this->c.emplace(name, std::vector<flow_float>{});
        this->c_d.emplace(name, nullptr);
    }
    for (const auto& name : outs) this->output_cellValNames.push_back(name);
    std::cout << "registerTransition: LM2009 gamma-Re_theta_t registered (" << names.size() << " cell variables)\n";
}

void variables::registerTracer(int enabled)
{
    if (enabled == 0) {
        this->tracerRegistered = 0;
        return;
    }
    this->tracerRegistered = 1;
    for (const auto& name : {"roXi", "Xi", "roXiN", "roXiM", "res_roXi", "res_roXi_m", "src_jac_Xi", "transport_diag_Xi"}) {
        this->cellValNames.push_back(name);
        this->c.emplace(name, std::vector<flow_float>{});
        this->c_d.emplace(name, nullptr);
    }
    this->output_cellValNames.push_back("roXi");
    this->output_cellValNames.push_back("Xi");
    // 受動種経路 (passiveScalarScheme 1; plan species-passive-scalar-unification §4.1): 勾配 ∇ξ・リミッタ ψ_ξ (S3 面再構成)、
    // scalar-DPLUR 増分 dq、floor 補正の累積 |Δ(ρξ)| (passiveFloorCorr_Xi; 確保時 0 初期化)。勾配/リミッタ/補正は level 2 出力のみ。
    for (const auto& name : {"dXidx", "dXidy", "dXidz", "limiter_Xi", "passiveFloorCorr_Xi", "passiveLimCorr_Xi", "dq_roXi", "dq_roXi_old", "roXiP", "roXiPP"}) {
        this->cellValNames.push_back(name);
        this->c.emplace(name, std::vector<flow_float>{});
        this->c_d.emplace(name, nullptr);
    }
    for (const auto& name : {"dXidx", "dXidy", "dXidz", "limiter_Xi", "passiveFloorCorr_Xi", "passiveLimCorr_Xi"}) this->output_cellValNames.push_back(name);
    std::cout << "registerTracer: exhaust tracer roXi registered (17 cell variables)\n";
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
        "src_jac_"+prim, "transport_diag_"+prim,
        // 受動種経路 (passiveScalarScheme 1): 勾配・リミッタ (S3)、scalar-DPLUR 増分、floor 補正の累積 |Δ(ρφ)|
        "d"+prim+"dx", "d"+prim+"dy", "d"+prim+"dz", "limiter_"+prim, "passiveFloorCorr_"+prim, "passiveLimCorr_"+prim,
        "dq_"+consName, "dq_"+consName+"_old",
        consName+"P", consName+"PP"   // dual-time の物理時間レベル (受動種 BDF 項用)
    };
}

void variables::registerTwoPhaseVaporResidual(int enabled)
{
    if (enabled == 0) return;
    // 蒸気の残差 R_v = R_w − R_g (監視のみ; HDF5 には出さない)。plans/active/condensation-two-phase-transport.md §5.1 #4e
    const std::string name = "res_roYv";
    this->cellValNames.push_back(name);
    this->c.emplace(name, std::vector<flow_float>{});
    this->c_d.emplace(name, nullptr);
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
            // 受動種経路の診断 (level 2 のみ): 勾配・リミッタ・floor 補正の累積。
            const std::string prim = consName.substr(2);
            for (const auto& name : {"d"+prim+"dx", "d"+prim+"dy", "d"+prim+"dz", "limiter_"+prim, "passiveFloorCorr_"+prim, "passiveLimCorr_"+prim}) {
                this->output_cellValNames.push_back(name);
            }
        }
        // 診断 (source kernel が毎ステップ書く): 過飽和 S=p_v/p_sat, 成長率 dr/dt [m/s] (負=蒸発),
        // 体積平均半径 r30 [m] (蒸発分岐で評価; 0=未評価)。確保時に 0 初期化 (prefix "cond")。
        // condTsat_: 飽和温度 T_sat(p_v) [K] (過冷却度 T_sat-T の評価用, plans/accepted/condensation-equilibrium.md)
        // condTheta_: 非等温核生成補正 θ (condKantrowitz 1–3; S>1 のセル), condLim_: 更新クランプ係数 θ_u (condLimiterMode 1; 収束時 ≈1)
        //   [旧 mode 0 では残差ソース律速係数], condClampCorr_: 更新後の硬クランプ (ρg<0 → 0) 補正量 [質量分率] (収束時 0)。
        for (const auto& d : {std::string("condS_"), std::string("condDrdt_"), std::string("condR30_"), std::string("condTsat_"),
                              std::string("condTheta_"), std::string("condLim_"), std::string("condClampCorr_"), std::string("condClampCorrQ_")}) {
            const std::string name = d + std::to_string(s);
            this->cellValNames.push_back(name);
            this->c.emplace(name, std::vector<flow_float>{});
            this->c_d.emplace(name, nullptr);
            this->output_cellValNames.push_back(name);
        }
    }

    std::cout << "registerCondensation: nCondSpecies=" << nCondSpecies
              << " -> registered " << nCondSpecies*(4*17+8) << " cell variables\n";
}

// --- FP64 影アキュムレータ (plans/active/time_integration-fp64-accumulator.md §4.3) ---
// **内点 CV だけ**確保する (nCells。ゴースト nCells_all-nCells は境界条件が毎 step 書くので正本を持たない)。
// 5 保存量 x 8 B = 40 B/CV。1000 万 CV で +400 MB。
static const char* const s_qaccNames[5] = {"ro", "roUx", "roUy", "roUz", "roe"};

void variables::allocQAccumulator(geom_int nCells)
{
    if (this->qacc_d[0] != nullptr) return;   // 二重確保を防ぐ
    for (int i = 0; i < 5; i++) {
        gpuErrchk( cudaMalloc((void**) &(this->qacc_d[i]), nCells*sizeof(double)) );
    }
    gpuErrchk( cudaMalloc((void**) &(this->qaccAdopt_d), sizeof(int)) );
    qaccResetAdoptCounter(this->qaccAdopt_d);
    std::cout << "allocQAccumulator: FP64 影アキュムレータ " << nCells << " CV x 5 変数 ("
              << (double)nCells*5.0*8.0/1024.0/1024.0 << " MB)\n";
}

void variables::initQAccumulatorFromQ(geom_int nCells)
{
    // **黙って no-op にしない** (2026-09-23): 以前はここで return していたため、確保より先に
    // 呼ばれていたことに気づけず、Qacc=0 のまま commit されて ro≈0 → 発散した。
    if (this->qacc_d[0] == nullptr) {
        std::cerr << "initQAccumulatorFromQ: 正本が未確保のまま呼ばれた "
                     "(allocQAccumulator より前に呼んでいる)\n";
        std::exit(1);
    }
    flow_float* q[5];
    for (int i = 0; i < 5; i++) q[i] = this->c_d.at(s_qaccNames[i]);
    qaccInitFromQ(this->qacc_d, q, nCells);
}

void variables::freeQAccumulator()
{
    for (int i = 0; i < 5; i++) {
        if (this->qacc_d[i] != nullptr) { cudaWrapper::cudaFree_wrapper(this->qacc_d[i]); this->qacc_d[i] = nullptr; }
    }
    if (this->qaccAdopt_d != nullptr) { cudaWrapper::cudaFree_wrapper(this->qaccAdopt_d); this->qaccAdopt_d = nullptr; }
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

    this->freeQAccumulator();
}

void variables::allocVariables(const int &useGPU , mesh& msh)
{
    // W-I 実力診断 (§4.2) は FORGE_WI_FORCE_DIAG=1 のときだけ有効。OFF なら
    // 確保も出力もしない (通常経路のメモリ・D2H・HDF5 を増やさない)。
    {
        const char* e = std::getenv("FORGE_WI_FORCE_DIAG");
        if (!(e && std::atoi(e) != 0)) {
            for (const char* n : {"wi_ftan", "wi_fnrm", "wi_fnrm_abs", "wi_ftan_res", "wi_eheat", "wi_ework"}) {
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
                || cellValName.rfind("passiveFloorCorr_", 0) == 0 || cellValName.rfind("passiveLimCorr_", 0) == 0 || cellValName.rfind("roYraw", 0) == 0
                || cellValName == "wf_irep_flag" || cellValName == "wf_sprod"
                || cellValName == "wf_g"
                || cellValName.rfind("lm", 0) == 0 || cellValName == "gammaEff" || cellValName == "gammaTr" || cellValName == "reTheta"
                || cellValName == "roGamma" || cellValName == "roReth" || cellValName.find("_gamma") != std::string::npos
                || cellValName.find("_reth") != std::string::npos || cellValName == "res_roGamma" || cellValName == "res_roReth") {
                gpuErrchk( cudaMemset(this->c_d.at(cellValName), 0, (msh.nCells_all)*sizeof(flow_float)) );
            }
        }
        // SST F1 の初期値は 1 (配列の確保時に入れる。以前は buildScalarDescs が初回に 1 で埋めており、直前に計算した F1 を
        // 上書きしていた。plan boundary-node-periodic-gradient-fix §4.2 / codex m5・実装レビュー m3)
        if (this->c_d.count("sstF1") && this->c_d.at("sstF1") != nullptr) {
            std::vector<flow_float> ones(msh.nCells_all, (flow_float)1.0);
            gpuErrchk( cudaMemcpy(this->c_d.at("sstF1"), ones.data(), sizeof(flow_float)*msh.nCells_all, cudaMemcpyHostToDevice) );
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

// 面ごとの差 e = cc[ic1] − cc[ic0] を double の値の位置 (mesh::cc64) で引き、1 回だけ flow_float に丸めて
// ホストの p["ge_x"/"ge_y"/"ge_z"] に入れる (plans/active/architecture-float-state-double-geometry.md §4.2a)。
// 全面 (内部面・周期面・境界面のゴースト側) が対象。段 ① では本番のカーネルは読まない。
// cc64 が無い mesh (readMesh を通らない) では NaN を入れる (誤って使えば見えるように)。
static void fillGeomDiffE(mesh& msh, variables& v)
{
    std::vector<flow_float>& ex = v.p.at("ge_x");
    std::vector<flow_float>& ey = v.p.at("ge_y");
    std::vector<flow_float>& ez = v.p.at("ge_z");
    ex.resize(msh.nPlanes); ey.resize(msh.nPlanes); ez.resize(msh.nPlanes);
    const bool has64 = (msh.cc64.size() == 3*(size_t)msh.nCells_all);
    if (!has64) {
        std::cout << "[variables] ge_x/ge_y/ge_z: mesh has no double centre copy (cc64); filled with NaN" << std::endl;
    }
    const flow_float nan = std::numeric_limits<flow_float>::quiet_NaN();
    for (geom_int ip = 0; ip < msh.nPlanes; ip++) {
        const auto& pc = msh.planes[ip].iCells;
        if (!has64 || pc.size() < 2) { ex[ip] = nan; ey[ip] = nan; ez[ip] = nan; continue; }
        const size_t c0 = 3*(size_t)pc[0], c1 = 3*(size_t)pc[1];
        const double dx = msh.cc64[c1+0] - msh.cc64[c0+0];
        const double dy = msh.cc64[c1+1] - msh.cc64[c0+1];
        const double dz = msh.cc64[c1+2] - msh.cc64[c0+2];
        ex[ip] = (flow_float)dx;
        ey[ip] = (flow_float)dy;
        ez[ip] = (flow_float)dz;
    }
}

// 面の両側の pc − cc を double の値の位置 (面重心 mesh::planeCent64、値の位置 mesh::cc64) で引き、1 回だけ flow_float に丸めて
// ホストの p["gr0_*"] (= pc[ip] − cc[ic0]) と p["gr1_*"] (= pc[ip] − cc[ic1]) に入れる
// (plans/active/architecture-float-state-double-geometry.md §4.2a、段 ④)。ic0/ic1 は fillGeomDiffE と同じ planes[ip].iCells。
// 全面 (内部面・周期面・境界面のゴースト側) が対象。デバイスの pcx..pcz は planes[].centCoords (= planeCent64 を geom_float に
// 丸めた値) なので、FP64 のビルドでは今の pcx − ccx と同じ値になる。double の写しが無い mesh では NaN を入れる。
static void fillGeomDiffR(mesh& msh, variables& v)
{
    std::vector<flow_float>* r[2][3] = {
        {&v.p.at("gr0_x"), &v.p.at("gr0_y"), &v.p.at("gr0_z")},
        {&v.p.at("gr1_x"), &v.p.at("gr1_y"), &v.p.at("gr1_z")} };
    for (int s = 0; s < 2; ++s) for (int k = 0; k < 3; ++k) r[s][k]->resize(msh.nPlanes);
    const bool has64 = (msh.cc64.size() == 3*(size_t)msh.nCells_all) && (msh.planeCent64.size() == 3*(size_t)msh.nPlanes);
    if (!has64) {
        std::cout << "[variables] gr0_*/gr1_*: mesh has no double copy (cc64/planeCent64); filled with NaN" << std::endl;
    }
    const flow_float nan = std::numeric_limits<flow_float>::quiet_NaN();
    for (geom_int ip = 0; ip < msh.nPlanes; ip++) {
        const auto& pc = msh.planes[ip].iCells;
        for (int s = 0; s < 2; ++s) {
            if (!has64 || pc.size() < 2) { for (int k = 0; k < 3; ++k) (*r[s][k])[ip] = nan; continue; }
            const size_t c = 3*(size_t)pc[s], p = 3*(size_t)ip;
            for (int k = 0; k < 3; ++k) {
                const double d = msh.planeCent64[p+k] - msh.cc64[c+k];
                (*r[s][k])[ip] = (flow_float)d;
            }
        }
    }
}

// 診断 (環境変数 FORGE_DIAG_GEOM_STAGE4_CHECK=1 のときだけ、既定 off で何もしない。plan §6.6 の 2):
// 段 ④ で再構成とリミタが読む値と、今のカーネルが作っていた座標の差のビット不一致の件数を数えて印字する。
// 今の差はデバイスの配列 (ccx..ccz・pcx..pcz・plane_cells) を写して flow_float で引く (カーネルと同じ 1 回の引き算)。
//   e      : ge  対 ccx[ic1] − ccx[ic0]
//   e_side1: 0 − ge 対 ccx[ic0] − ccx[ic1] (リミタが ic1 側のセルから見る辺中点。符号付きゼロも比べる)
//   r0 / r1: gr0 対 pcx − ccx[ic0]、gr1 対 pcx − ccx[ic1]
//   bnd_r0 : 境界カーネル (convectiveFlux_boundary_d) の ic = bc.iCells[ib] で作る pcx − ccx[ic] 対 gr0
//            (と、ic ≠ plane_cells[2ip+0] の件数)
// 内部面 (ip < nNormalPlanes) とそれ以外を分けて数える。FP64 のビルドでは全項目 0 になるはず。
// float のビルドでは、double の写しで引いた値が float の座標の差と違う面が数えられる (それが段 ④ の目的)。
static void checkGeomStage4(mesh& msh, variables& v)
{
    const char* env = std::getenv("FORGE_DIAG_GEOM_STAGE4_CHECK");
    if (env == nullptr || *env == '\0' || std::string(env) == "0") return;

    const size_t nP = (size_t)msh.nPlanes, nCa = (size_t)msh.nCells_all;
    auto d2h = [](const flow_float* d, size_t n) {
        std::vector<flow_float> h(n);
        cudaMemcpy(h.data(), d, n*sizeof(flow_float), cudaMemcpyDeviceToHost);
        return h;
    };
    std::vector<geom_int> pc(2*nP);
    cudaMemcpy(pc.data(), msh.map_plane_cells_d, 2*nP*sizeof(geom_int), cudaMemcpyDeviceToHost);
    const std::vector<flow_float> cc[3] = {d2h(v.c_d.at("ccx"), nCa), d2h(v.c_d.at("ccy"), nCa), d2h(v.c_d.at("ccz"), nCa)};
    const std::vector<flow_float> pcc[3] = {d2h(v.p_d.at("pcx"), nP), d2h(v.p_d.at("pcy"), nP), d2h(v.p_d.at("pcz"), nP)};
    const std::vector<flow_float> ge[3] = {d2h(v.p_d.at("ge_x"), nP), d2h(v.p_d.at("ge_y"), nP), d2h(v.p_d.at("ge_z"), nP)};
    const std::vector<flow_float> r0[3] = {d2h(v.p_d.at("gr0_x"), nP), d2h(v.p_d.at("gr0_y"), nP), d2h(v.p_d.at("gr0_z"), nP)};
    const std::vector<flow_float> r1[3] = {d2h(v.p_d.at("gr1_x"), nP), d2h(v.p_d.at("gr1_y"), nP), d2h(v.p_d.at("gr1_z"), nP)};

    auto sameBits = [](flow_float a, flow_float b) { return std::memcmp(&a, &b, sizeof(flow_float)) == 0; };
    // 項目ごと・面の種類ごと (0 = 内部面、1 = それ以外) の件数と、ベクトルの相対差 ‖new − old‖/‖new‖ の最大
    enum { kE = 0, kE1, kR0, kR1, kN };
    const char* names[kN] = {"e", "e_side1", "r0", "r1"};
    unsigned long long mis[kN][2] = {}, nonfin[kN][2] = {}, nface[2] = {0, 0};
    double maxRel[kN][2] = {};
    for (size_t ip = 0; ip < nP; ++ip) {
        const int cls = (ip < (size_t)msh.nNormalPlanes) ? 0 : 1;
        const geom_int i0 = pc[2*ip+0], i1 = pc[2*ip+1];
        if (i0 < 0 || i1 < 0 || (size_t)i0 >= nCa || (size_t)i1 >= nCa) continue;
        ++nface[cls];
        for (int it = 0; it < kN; ++it) {
            flow_float nw[3], od[3];
            for (int k = 0; k < 3; ++k) {
                switch (it) {
                case kE:  nw[k] = ge[k][ip];                   od[k] = cc[k][i1] - cc[k][i0];   break;
                case kE1: nw[k] = (flow_float)0.0 - ge[k][ip]; od[k] = cc[k][i0] - cc[k][i1];   break;
                case kR0: nw[k] = r0[k][ip];                   od[k] = pcc[k][ip] - cc[k][i0];  break;
                default:  nw[k] = r1[k][ip];                   od[k] = pcc[k][ip] - cc[k][i1];  break;
                }
            }
            bool bad = false, fin = true;
            double dn = 0.0, nn = 0.0;
            for (int k = 0; k < 3; ++k) {
                if (!sameBits(nw[k], od[k])) bad = true;
                if (!std::isfinite((double)nw[k])) fin = false;
                dn += ((double)nw[k] - (double)od[k])*((double)nw[k] - (double)od[k]);
                nn += (double)nw[k]*(double)nw[k];
            }
            if (bad) ++mis[it][cls];
            if (!fin) { ++nonfin[it][cls]; continue; }
            if (nn > 0.0) maxRel[it][cls] = std::max(maxRel[it][cls], std::sqrt(dn/nn));
        }
    }
    // 境界カーネルの ic (= bc.iCells[ib]) で作る差と gr0 の照合
    unsigned long long bndN = 0, bndCellMis = 0, bndR0Mis = 0;
    for (const bcond& bc : msh.bconds) {
        for (size_t ib = 0; ib < bc.iPlanes.size() && ib < bc.iCells.size(); ++ib) {
            const geom_int ip = bc.iPlanes[ib], ic = bc.iCells[ib];
            if (ip < 0 || (size_t)ip >= nP || ic < 0 || (size_t)ic >= nCa) continue;
            ++bndN;
            if (ic != pc[2*(size_t)ip+0]) ++bndCellMis;
            bool bad = false;
            for (int k = 0; k < 3; ++k) if (!sameBits(r0[k][ip], (flow_float)(pcc[k][ip] - cc[k][ic]))) bad = true;
            if (bad) ++bndR0Mis;
        }
    }
    printf("[geomStage4] FORGE_DIAG_GEOM_STAGE4_CHECK: sizeof(flow_float)=%zu, has_geom64=%s, faces internal=%llu other=%llu\n",
           sizeof(flow_float), msh.hasGeom64() ? "yes" : "NO", nface[0], nface[1]);
    for (int it = 0; it < kN; ++it) {
        printf("[geomStage4]   %-8s bit-mismatch internal=%llu other=%llu  nonfinite internal=%llu other=%llu"
               "  max|new-old|/|new| internal=%.3e other=%.3e\n",
               names[it], mis[it][0], mis[it][1], nonfin[it][0], nonfin[it][1], maxRel[it][0], maxRel[it][1]);
    }
    printf("[geomStage4]   bnd_r0   boundary faces=%llu  bc.iCells != plane_cells[2ip+0]: %llu  bit-mismatch vs pcx-ccx[bc.iCells]=%llu\n",
           bndN, bndCellMis, bndR0Mis);
    fflush(stdout);
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

    // 面ごとの差 e (double の座標から作る、§4.2a)。段 ① では読む経路なし。
    fillGeomDiffE(msh, *this);
    // 面の両側の pc − cc (同、段 ④)
    fillGeomDiffR(msh, *this);
}

// 区間ごとの r 重みの面ベクトルを使うか (宣言の説明は variables.hpp)。plans/active/axisymmetric-freestream-hoop-gauge.md §4.5:
// データセットの有無だけでは有効にせず、設定と格子の条件をすべて満たすときだけ使う。
bool axisSegmentRWeightApplies(const solverConfig& cfg, const mesh& msh, std::string* reason)
{
    auto no = [reason](const char* r) { if (reason) *reason = r; return false; };
    if (cfg.axisSegmentRWeight != 1)        return no("mesh.axisSegmentRWeight 0");
    if (cfg.isAxisymmetric != 1)            return no("isAxisymmetric 0");
    if (cfg.axisymMethod != 0)              return no("axisymMethod 1");
    if (cfg.axisRFloor > (flow_float)0.0)   return no("axisRFloor > 0");
    if (cfg.discretization != "node")       return no("not node");
    if (msh.rSurfVect64.empty())            return no("no /PLANES/rSurfVect in the mesh (reconvert with isAxisymmetric 1 to use it)");
    if (msh.rSurfVect64.size() != 3*(size_t)msh.nPlanes || msh.surfVect64.size() != 3*(size_t)msh.nPlanes
     || msh.planeCent64.size() != 3*(size_t)msh.nPlanes || msh.coord64.size() < 3*msh.nodes.size() || msh.nodes.empty())
        return no("double geometry copies missing or of the wrong size");
    // 平面の 2D (true 2D): 全節点の z が同じで、全面ベクトルの z 成分が 0 (押し出しの疑似 2D や 3D を除く)
    const double z0 = msh.coord64[2];
    for (size_t i = 0; i < msh.nodes.size(); ++i)
        if (msh.coord64[3*i + 2] != z0) return no("not planar 2D (node z differs)");
    for (geom_int ip = 0; ip < msh.nPlanes; ++ip)
        if (msh.surfVect64[3*(size_t)ip + 2] != 0.0) return no("not planar 2D (face vector has a z component)");
    if (reason) reason->clear();
    return true;
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

    // 区間ごとの r 重みの面ベクトル (plans/active/axisymmetric-freestream-hoop-gauge.md §4.5) を使うか。起動時に 1 行出す。
    std::string segRWReason;
    const bool segRW = axisSegmentRWeightApplies(cfg, msh, &segRWReason);
    if (!segRW && (cfg.isAxisymmetric == 1 || !msh.rSurfVect64.empty()))
        printf("[axisym] axisSegmentRWeight: OFF (%s) -> face vectors r̄_f·S_f as before\n", segRWReason.c_str());

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
        if (segRW) {
            // 区間ごとの r 重み W_f = Σ_k r_k S_k (変換器が double で作った /PLANES/rSurfVect)。折れた双対面を 1 本にまとめてから
            // 重心の半径 r̄_f を掛けると Σ r_k S_k と食い違い、FP64 でも軸の近くに偽の半径力が立つ (§4.1 #2)。
            // sx..sz = W_f を丸めたもの、ss = ‖W_f‖ (double で取ってから丸める。法線 sx/ss の整合)。
            // 軸の上の面 (double の面重心の半径 ≤ 1e-20、幾何の条件) は今と同じく S·r_floor で向きを保つ (床の後のベクトルの
            // ノルムを double で取る)。axisRFloor 0 が条件なので床の値は 1e-20。
            const double rFloor64 = 1.0e-20;
            geom_int nAxisFace = 0;
            for (geom_int ip=0; ip<msh.nPlanes; ip++) {
                const size_t p = 3*(size_t)ip;
                double wx, wy, wz;
                if (msh.planeCent64[p + 1] <= rFloor64) {
                    wx = msh.surfVect64[p + 0]*rFloor64;
                    wy = msh.surfVect64[p + 1]*rFloor64;
                    wz = msh.surfVect64[p + 2]*rFloor64;
                    ++nAxisFace;
                } else {
                    wx = msh.rSurfVect64[p + 0];
                    wy = msh.rSurfVect64[p + 1];
                    wz = msh.rSurfVect64[p + 2];
                }
                const double wn = std::sqrt(wx*wx + wy*wy + wz*wz);
                sx[ip] = (geom_float)wx;
                sy[ip] = (geom_float)wy;
                sz[ip] = (geom_float)wz;
                ss[ip] = (geom_float)wn;
                if (!(std::isfinite((double)ss[ip]) && ss[ip] > (geom_float)0.0 && std::isfinite((double)sx[ip])
                      && std::isfinite((double)sy[ip]) && std::isfinite((double)sz[ip]))) {
                    fprintf(stderr, "[axisym] ERROR: axisSegmentRWeight: face %lld has a non-finite or non-positive area after "
                            "rounding (W = %.17g %.17g %.17g, |W| = %.17g, ss = %.9g, face centroid r = %.17g)\n",
                            (long long)ip, wx, wy, wz, wn, (double)ss[ip], msh.planeCent64[p + 1]);
                    exit(EXIT_FAILURE);
                }
            }
            printf("[axisym] axisSegmentRWeight: ON -> face vectors W_f = sum_k r_k S_k from /PLANES/rSurfVect "
                   "(%lld faces, %lld axis faces at r = 1e-20)\n", (long long)msh.nPlanes, (long long)nAxisFace);
        } else {
            for (geom_int ip=0; ip<msh.nPlanes; ip++) {
                const geom_float r_face = (pcy[ip] > r_floor) ? pcy[ip] : r_floor;
                sx[ip] *= r_face;
                sy[ip] *= r_face;
                sz[ip] *= r_face;
                ss[ip] *= r_face;
            }
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
            // 足し上げるのはデバイスに渡す最終の面ベクトル (半径の重みを掛けて geom_float に丸めた sx/sy) のまま。
            // 和は double に貯め、最後に 1 回だけ丸める (plans/active/architecture-float-state-double-geometry.md §4.2b、段 ②。
            // 薄い半径方向のセルでは S·r_top − S·r_bot が打ち消し合うので、float の和では丸めが残る。FP64 のビルドでは従来と同じ)。
            std::vector<double> aclx64(msh.nCells_all, 0.0), acly64(msh.nCells_all, 0.0);
            for (geom_int ip = 0; ip < msh.nPlanes; ++ip) {
                const auto& pc = msh.planes[ip].iCells;
                if (pc.empty()) continue;
                const geom_int ic0 = pc[0];
                if (ic0 >= 0 && ic0 < msh.nCells) { aclx64[ic0] += (double)sx[ip]; acly64[ic0] += (double)sy[ip]; }
                if (pc.size() > 1) {
                    const geom_int ic1 = pc[1];
                    if (ic1 >= 0 && ic1 < msh.nCells) { aclx64[ic1] -= (double)sx[ip]; acly64[ic1] -= (double)sy[ip]; }
                }
            }
            std::vector<geom_float> aclx(msh.nCells_all), acly(msh.nCells_all);
            for (geom_int ic = 0; ic < msh.nCells_all; ++ic) {
                aclx[ic] = (geom_float)aclx64[ic];
                acly[ic] = (geom_float)acly64[ic];
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
    // 段 ② (plans/active/architecture-float-state-double-geometry.md §4.2 3.): 重心間の差と距離・最大は double の値の位置
    // (mesh::cc64) で取り、最後に 1 回だけ丸める (丸めは単調なので、丸めた値の最大と同じ順序になる。FP64 のビルドでは従来と同じ)。
    // double の写しを持たない mesh は従来の geom_float の差。
    {
        const bool use64 = msh.hasGeom64();
        std::vector<double> delta_les64(msh.nCells_all, 0.0);
        geom_float* delta_les_h = (geom_float*)malloc(sizeof(geom_float)*msh.nCells_all);
        for (geom_int ic=0; ic<msh.nCells_all; ic++) delta_les_h[ic] = 0.0;
        for (geom_int ip=0; ip<msh.nPlanes; ip++) {
            const geom_int ic1 = msh.planes[ip].iCells[0];
            const geom_int ic2 = msh.planes[ip].iCells[1];
            if (use64) {
                const double dx = msh.cc64[3*(size_t)ic1 + 0] - msh.cc64[3*(size_t)ic2 + 0];
                const double dy = msh.cc64[3*(size_t)ic1 + 1] - msh.cc64[3*(size_t)ic2 + 1];
                const double dz = msh.cc64[3*(size_t)ic1 + 2] - msh.cc64[3*(size_t)ic2 + 2];
                const double d  = sqrt(dx*dx + dy*dy + dz*dz);
                if (ic1 < msh.nCells && d > delta_les64[ic1]) delta_les64[ic1] = d;
                if (ic2 < msh.nCells && d > delta_les64[ic2]) delta_les64[ic2] = d;
            } else {
                const geom_float dx = ccx[ic1] - ccx[ic2];
                const geom_float dy = ccy[ic1] - ccy[ic2];
                const geom_float dz = ccz[ic1] - ccz[ic2];
                const geom_float d  = sqrt(dx*dx + dy*dy + dz*dz);
                if (ic1 < msh.nCells && d > delta_les_h[ic1]) delta_les_h[ic1] = d;
                if (ic2 < msh.nCells && d > delta_les_h[ic2]) delta_les_h[ic2] = d;
            }
        }
        if (use64) {
            for (geom_int ic=0; ic<msh.nCells_all; ic++) delta_les_h[ic] = (geom_float)delta_les64[ic];
        }
        cudaMemcpy(this->c_d.at("delta_les"), delta_les_h, msh.nCells_all*sizeof(geom_float), cudaMemcpyHostToDevice);
        free(delta_les_h);
    }

    // 面ごとの差 e = fl(cc64[ic1] − cc64[ic0]) (plans/active/architecture-float-state-double-geometry.md §4.2a、段 ①)。
    // 段 ③ から粘性・拡散・陰解法の対角、段 ④ から再構成の辺中点 (±0.5·e) とリミタが読む。
    fillGeomDiffE(msh, *this);
    this->copyVariables_plane_H2D({"ge_x", "ge_y", "ge_z"});
    // 面の両側の pc − cc: r0 = fl(pc64 − cc64[ic0])、r1 = fl(pc64 − cc64[ic1]) (同 plan §4.2a、段 ④)。
    // 再構成 (辺中点でない面・cell の双対面重心・境界面) とリミタが読む。
    fillGeomDiffR(msh, *this);
    this->copyVariables_plane_H2D({"gr0_x", "gr0_y", "gr0_z", "gr1_x", "gr1_y", "gr1_z"});

    calcStructualVariables_d_wrapper(cfg , cuda_cfg , msh , *this);

    // 診断 (FORGE_DIAG_GEOM_STAGE4_CHECK=1、既定 off): e・r0・r1 と今の座標の差のビット不一致の件数を印字する (plan §6.6 の 2)。
    checkGeomStage4(msh, *this);

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

    // --- 遷移モデル: ργ, ρRe_θt を読み込む。無ければ初回の transitionPrimitive が γ=1 / 自由流相関で初期化する ---
    if (this->transitionRegistered != 0) {
        const bool has = file.exist("/VALUE/roGamma") && file.exist("/VALUE/roReth");
        if (has) {
            std::vector<geom_float> g, r;
            file.getDataSet("/VALUE/roGamma").read(g);
            file.getDataSet("/VALUE/roReth").read(r);
            std::vector<flow_float>& vg = this->c.at("roGamma");
            std::vector<flow_float>& vr = this->c.at("roReth");
            for (geom_int i=0; i<msh.nCells; i++) { vg[i] = g[i]; vr[i] = r[i]; }
            this->copyVariables_cell_H2D({"roGamma", "roReth"});
        }
        this->transitionNeedsInit = has ? 0 : 1;
        std::cout << "[variables] transition roGamma/roReth " << (has ? "read from input" : "not in input: will be initialised on the first step") << "\n";
    }

    // --- 受動トレーサ: ρξ を読み込む (VALUE/roXi → VALUE/Xi×ρ → 0 の優先順) ---
    if (this->tracerRegistered != 0) {
        std::vector<flow_float>& v_roXi = this->c.at("roXi");
        std::vector<flow_float>& v_Xi   = this->c.at("Xi");
        std::vector<geom_float> in;
        bool has = false;
        if (file.exist("/VALUE/roXi")) {
            file.getDataSet("/VALUE/roXi").read(in); has = true;
        } else if (file.exist("/VALUE/Xi")) {
            std::vector<geom_float> xi_in;
            file.getDataSet("/VALUE/Xi").read(xi_in);
            in.resize(xi_in.size());
            for (std::size_t k = 0; k < xi_in.size(); k++) in[k] = static_cast<geom_float>(this->c.at("ro")[k]) * xi_in[k];
            has = true;
        }
        for (geom_int i=0; i<msh.nCells; i++) {
            const flow_float roi = this->c.at("ro")[i];
            flow_float v = has ? static_cast<flow_float>(in[i]) : static_cast<flow_float>(0.0);
            if (v < 0.0) v = 0.0;
            if (v > roi) v = roi;
            v_roXi[i] = v;
            v_Xi[i]   = v / std::max(roi, static_cast<flow_float>(1.0e-30));
        }
        std::cout << "[variables] tracer roXi " << (has ? "read from input" : "not in input: initialized to 0") << "\n";
        this->copyVariables_cell_H2D({"roXi", "Xi"});
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