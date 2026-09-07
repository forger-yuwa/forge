#include "input/solverConfig.hpp"


solverConfig::solverConfig(){};

// getValidatedValue 関数の定義
template <typename T>
T getValidatedValue(const YAML::Node& node, const std::string& key, const std::string& parent = "") {
    try {
        if (!node[key].IsDefined()) {
            throw std::runtime_error("Key '" + key + "' is missing" + (parent.empty() ? "" : " in '" + parent + "'."));
        }
        T value = node[key].as<T>();
        std::cout << "'" << key << "'" << (parent.empty() ? "" : " in '" + parent + "'") << ": " << value << std::endl;
        return value;
    } catch (const YAML::BadConversion& e) {
        throw std::runtime_error("Key '" + key + "' in '" + parent + "' has an incorrect type.");
    }
}

template <typename T>
T getOptionalValidatedValue(const YAML::Node& node, const std::string& key, const T& default_value, const std::string& parent = "") {
    if (!node[key].IsDefined()) {
        return default_value;
    }

    try {
        T value = node[key].as<T>();
        std::cout << "'" << key << "'" << (parent.empty() ? "" : " in '" + parent + "'") << ": " << value << std::endl;
        return value;
    } catch (const YAML::BadConversion&) {
        throw std::runtime_error("Key '" + key + "' in '" + parent + "' has an incorrect type.");
    }
}


void solverConfig::read(std::string fname)
{
//    try {
//
//
//        this->solConfigFileName = fname;
//        YAML::Node config = YAML::LoadFile(this->solConfigFileName);
//
//        
//        std::string meshFormat    = config["mesh"]["meshFormat"].as<std::string>();
//        std::string meshFileName  = config["mesh"]["meshFileName"].as<std::string>();
//        std::string valueFileName = config["mesh"]["valueFileName"].as<std::string>();
//        std::string solver = config["solver"].as<std::string>();
//
//        this->gpu = config["gpu"].as<int>();
//
//        int unsteady = config["time"]["unsteady"].as<int>();
//        int dualTime = config["time"]["dualTime"].as<int>();
//
//        auto last = config["time"]["last"];
//            int endTimeControl = last["control"].as<int>();
//            int nStep = last["nStep"].as<int>();
//            int endTime = last["time"].as<flow_float>();
//
//        auto deltaT = config["time"]["deltaT"];
//            int dtControl = deltaT["control"].as<int>();
//            flow_float dt = deltaT["dt"].as<flow_float>();
//            flow_float cfl = deltaT["cfl"].as<flow_float>();
//            flow_float cfl_pseudo = deltaT["cfl_pseudo"].as<flow_float>();
//            flow_float dt_max= deltaT["dt_max"].as<flow_float>();
//            flow_float dt_min= deltaT["dt_min"].as<flow_float>();
//
//        int outStepInterval = config["time"]["outStepInterval"].as<int>();
//        int outStepStart    = config["time"]["outStepStart"].as<int>();
//        int timeIntegration = config["time"]["timeIntegration"].as<int>();
//        int nInnerLoop      = config["time"]["nInnerLoop"].as<int>();
//
//        auto space = config["space"];
//            int convMethod = space["convMethod"].as<flow_float>();
//            int limiter    = space["limiter"].as<int>();
//
//        auto turb = config["turbulence"];
//            int LESorRANS = turb["LESorRANS"].as<int>();
//            int LESmodel  = turb["LESmodel"].as<int>();
//
//        auto physProp = config["physProp"];
//            int isCompressible = physProp["isCompressible"].as<int>();
//            int thermalMethod = physProp["thermalMethod"].as<int>();
//            int viscMethod = physProp["viscMethod"].as<int>();
//            printf("viscMethod=%d\n", viscMethod);
//            flow_float ro = physProp["ro"].as<flow_float>();
//            flow_float visc = physProp["visc"].as<flow_float>();
//            flow_float thermCond = physProp["thermCond"].as<flow_float>();
//            flow_float cp = physProp["cp"].as<flow_float>();
//            flow_float gamma = physProp["gamma"].as<flow_float>();
//
//
//        this->initial = config["initial"].as<std::string>();
//
//        std::cout << "Mesh Name : " << meshFileName << '\n';
//        std::cout << "Mesh Name size : " << meshFileName.size() << '\n';
//        std::cout << "Value File Name : " << valueFileName << '\n';
//
//        std::cout << "Solver Name : " << solver << '\n';
//        if        (endTimeControl == 0) { std::cout << "End Step : " << nStep << '\n';
//        } else if (endTimeControl == 1) { std::cout << "End Time : " << endTime << '\n';
//        } else {
//            std::cerr << "something wrong in solver yaml" << std::endl;
//            exit(EXIT_FAILURE);
//        }
//
//        if        (dtControl == 0) { std::cout << "delta Time : " << dt << '\n'; 
//        } else if (dtControl == 1) { std::cout << "CFL : " << cfl << '\n';
//        } else {
//            std::cerr << "something wrong in solver yaml" << std::endl;
//            exit(EXIT_FAILURE);
//        }
//
//        std::cout << "convection method : " << convMethod<< '\n';
//        std::cout << "limiter : " << limiter<< '\n';
//
//        this->meshFormat   = meshFormat;
//        this->meshFileName = meshFileName;
//        this->valueFileName = valueFileName;
//        this->solver = solver;
//
//        this->endTimeControl = endTimeControl;
//        this->nStep = nStep;
//        this->outStepInterval = outStepInterval;
//        this->outStepStart    = outStepStart;
//
//        this->dtControl = dtControl;
//        this->dt = dt;
//        this->cfl = cfl;
//        this->cfl_pseudo = cfl_pseudo;
//        this->dt_max = dt_max;
//        this->dt_min = dt_min;
//        this->unsteady = unsteady;
//        this->dualTime = dualTime;
//        this->timeIntegration = timeIntegration;
//        this->nInnerLoop = nInnerLoop;
//        this->cfl_pseudo = cfl_pseudo;
//
//        initTimeIntegrationScheme(timeIntegration);
//
//        this->convMethod = convMethod;
//        this->limiter = limiter;
//
//
//        this->LESorRANS = LESorRANS;
//        this->LESmodel  = LESmodel;
//
//
//        this->isCompressible = isCompressible ;
//
//        this->thermalMethod = thermalMethod;
//        this->viscMethod = viscMethod;
//        this->ro = ro;
//        this->visc = visc;
//        this->thermCond = thermCond;
//        this->cp = cp;
//        this->gamma = gamma;
//
//        this->initial = initial;
//
//    } catch (const YAML::BadFile& e) {
//        std::cerr << e.msg << std::endl;
//        exit(EXIT_FAILURE);
//
//    } catch (const YAML::ParserException& e) {
//        std::cerr << e.msg << std::endl;
//        exit(EXIT_FAILURE);
//    }
//
//
    try {
        //YAML::Node config = YAML::LoadFile("config.yaml");
        this->solConfigFileName = fname;
        YAML::Node config = YAML::LoadFile(this->solConfigFileName);

        // mesh関連
        this->meshFormat = getValidatedValue<std::string>(config["mesh"], "meshFormat", "mesh");
        this->meshFileName = getValidatedValue<std::string>(config["mesh"], "meshFileName", "mesh");
        this->valueFileName = getValidatedValue<std::string>(config["mesh"], "valueFileName", "mesh");

        // 離散化レイアウト (任意, 既定 "cell")。"cell": cell-centered、"node": node-centered
        // (中点双対 median-dual)。methods/discretization/ 参照。
        if (config["mesh"]["discretization"]) {
            this->discretization = config["mesh"]["discretization"].as<std::string>();
            if (this->discretization != "cell" && this->discretization != "node") {
                throw std::runtime_error("Key 'discretization' in 'mesh' must be 'cell' or 'node'.");
            }
        } else {
            this->discretization = "cell";
        }
        if (config["mesh"]["bndFirstOrder"]) {
            this->bndFirstOrder = config["mesh"]["bndFirstOrder"].as<int>();
        }
        if (config["mesh"]["axisCentroidShift"]) {
            this->axisCentroidShift = config["mesh"]["axisCentroidShift"].as<int>();
        }
        if (config["mesh"]["nodeWallDirichlet"]) {
            this->nodeWallDirichlet = config["mesh"]["nodeWallDirichlet"].as<int>();
        }
        if (config["mesh"]["hoopAreaFromClosure"]) {
            this->hoopAreaFromClosure = config["mesh"]["hoopAreaFromClosure"].as<int>();
        }
        for (const char* k : {"nodeAxisDirichlet", "nodeMidpointFx", "nodeValueAtNode", "nodeReconEdgeMidpoint", "nodeAxisUrDirichlet"}) {
            if (config["mesh"][k]) {
                std::cerr << "[config] mesh." << k << " は撤去されました (node の整合セットは常時 ON, "
                          << "plan architecture-node-option-consolidation)。キーを削除してください。" << std::endl;
                std::exit(1);
            }
        }
        if (config["mesh"]["gradLSQ"]) {
            this->gradLSQ = config["mesh"]["gradLSQ"].as<int>();
            if (this->gradLSQ < 0 || this->gradLSQ > 2) {
                throw std::runtime_error("'gradLSQ' in 'mesh' must be 0 (GG), 1 (LSQ per-step) or 2 (LSQ precomputed).");
            }
        }
        if (this->discretization == "node") {
            // node は LSQ 固定 (内部隣接ノードのみ・点値整合)。GG は伸縮した壁行で市松を作る (case/43)。
            if (config["mesh"]["gradLSQ"] && this->gradLSQ != 2) {
                std::cerr << "[config] node では mesh.gradLSQ は 2 (LSQ 事前計算) 固定です (GG は壁行市松の原因, case/43)。" << std::endl;
                std::exit(1);
            }
            this->gradLSQ = 2;
        }
        if (config["mesh"]["gradLSQDegenThresh"]) {
            this->gradLSQDegenThresh = config["mesh"]["gradLSQDegenThresh"].as<double>();
        }
        if (config["mesh"]["nodeWallStressEdgeKernel"]) {
            this->nodeWallStressEdgeKernel = config["mesh"]["nodeWallStressEdgeKernel"].as<int>();
        }

        // solver関連
        this->solver = getValidatedValue<std::string>(config, "solver");

        // KEEP 用 opt-in 散逸レイヤ (solver=="KEEP" のみ有効。plans/active/convection-keep-es-dissipation.md)
        // 一様体積力 (top-level: bodyForce: [fx, fy, fz])。既定 [0,0,0] = off。
        if (config["bodyForce"]) {
            std::vector<double> f;
            for (const auto& v : config["bodyForce"]) f.push_back(v.as<double>());
            if (f.size() != 3) {
                throw std::runtime_error("Key 'bodyForce' must be a list of 3 floats [fx, fy, fz].");
            }
            this->bodyForceX = f[0]; this->bodyForceY = f[1]; this->bodyForceZ = f[2];
            std::cout << "'bodyForce': [" << f[0] << ", " << f[1] << ", " << f[2] << "]\n";
        }

        // 質量流量一定制御 (plans/active/timeint-bodyforce-massflow-control.md)。
        // bodyForceCtrl=1 で毎物理ステップ bodyForceX を deadbeat 更新し
        // 体積平均 ρu_x を bodyForceCtrlTarget に保つ。既定 0 = 固定値駆動 (ビット不変)。
        this->bodyForceCtrl = getOptionalValidatedValue<int>(config, "bodyForceCtrl", 0, "");
        this->bodyForceCtrlTarget = getOptionalValidatedValue<double>(config, "bodyForceCtrlTarget", 0.0, "");
        this->bodyForceCtrlRelax = getOptionalValidatedValue<double>(config, "bodyForceCtrlRelax", 1.0, "");

        this->keepDissType  = getOptionalValidatedValue<int>(config, "keepDissType", 0, "");
        this->keepDissCoeff = getOptionalValidatedValue<double>(config, "keepDissCoeff", 0.05, "");
        this->keepDissCprime = getOptionalValidatedValue<int>(config, "keepDissCprime", 1, "");
        this->keepDissJump = getOptionalValidatedValue<int>(config, "keepDissJump", 0, "");
        this->keepDissPrecond = getOptionalValidatedValue<int>(config, "keepDissPrecond", 0, "");
        this->keepDissFdBlend = getOptionalValidatedValue<int>(config, "keepDissFdBlend", 0, "");
        this->keepDissCoeffMax = getOptionalValidatedValue<double>(config, "keepDissCoeffMax", 1.0, "");
        this->keepDissCbCoeff = getOptionalValidatedValue<double>(config, "keepDissCbCoeff", 0.0, "");
        this->keepDissCbEps   = getOptionalValidatedValue<double>(config, "keepDissCbEps", 0.10, "");
        this->keepDissOpBlendRaw = getOptionalValidatedValue<int>(config, "keepDissOpBlendRaw", 0, "");
        if (this->keepDissOpBlendRaw == 1 && (this->keepDissFdBlend != 1 || this->keepDissPrecond != 1 || this->keepDissJump != 2)) {
            throw std::runtime_error("'keepDissOpBlendRaw: 1' requires the operator blend (keepDissFdBlend: 1 + keepDissPrecond: 1) with keepDissJump: 2.");
        }
        if (this->keepDissCbCoeff < 0.0) {
            throw std::runtime_error("Key 'keepDissCbCoeff' must be >= 0 (0 = off).");
        }
        if (this->keepDissCbCoeff > 0.0) {
            // 高周波圧力欠陥 δp^HF は面再構成 (keepDissJump>=1) の不整合として抽出するため、
            // 再構成なしでは市松狙い撃ちにならず生 Δp 拡散になる → ハードエラー。
            if (this->keepDissJump < 1) {
                throw std::runtime_error("'keepDissCbCoeff > 0' requires keepDissJump >= 1 "
                                         "(the HF pressure defect is the face-reconstruction mismatch).");
            }
            if (this->keepDissType != 2) {
                throw std::runtime_error("'keepDissCbCoeff > 0' requires keepDissType: 2 (matrix ES dissipation, CPG).");
            }
            if (this->keepDissCbEps <= 0.0 || this->keepDissCbEps > 1.0) {
                throw std::runtime_error("Key 'keepDissCbEps' must be in (0, 1].");
            }
        }
        if (this->keepDissPrecond < 0 || this->keepDissPrecond > 1) {
            throw std::runtime_error("Key 'keepDissPrecond' must be 0 (off) or 1 (Turkel-preconditioned acoustic dissipation, matrix CPG only).");
        }
        if (this->keepDissJump < 0 || this->keepDissJump > 2) {
            throw std::runtime_error("Key 'keepDissJump' must be 0 (raw jump), 1 (reconstructed jump), "
                                     "or 2 (reconstructed + sign-property clip, provably entropy-stable).");
        }
        if (this->keepDissType < 0 || this->keepDissType > 2) {
            std::cerr << "Error: keepDissType must be 0 (off), 1 (scalar ES) or 2 (matrix ES), got "
                      << this->keepDissType << std::endl;
            std::exit(EXIT_FAILURE);
        }

        // GPU設定
        this->gpu = getValidatedValue<int>(config, "gpu");

        // time関連
        this->unsteady = getValidatedValue<int>(config["time"], "unsteady", "time");
        this->dualTime = getValidatedValue<int>(config["time"], "dualTime", "time");

        auto last = config["time"]["last"];
        this->endTimeControl = getValidatedValue<int>(last, "control", "time.last");
        if (last["nStep"].IsDefined()) {
            throw std::runtime_error(
                "Key 'nStep' in 'time.last' is no longer supported. Rename it to 'nStepOuter'."
            );
        }
        this->nStepOuter = getValidatedValue<int>(last, "nStepOuter", "time.last");

        auto deltaT = config["time"]["deltaT"];
        this->dtControl = getValidatedValue<int>(deltaT, "control", "time.deltaT");
        this->dt = getValidatedValue<double>(deltaT, "dt", "time.deltaT");
        this->cfl = getValidatedValue<double>(deltaT, "cfl", "time.deltaT");
        this->cfl_pseudo = getValidatedValue<double>(deltaT, "cfl_pseudo", "time.deltaT");
        this->implicitRelax = getOptionalValidatedValue<double>(deltaT, "implicitRelax", 1.0, "time.deltaT");
        this->updateGuardAlpha = getOptionalValidatedValue<double>(deltaT, "updateGuardAlpha", 0.0, "time.deltaT");
        this->lineImplicit = getOptionalValidatedValue<int>(deltaT, "lineImplicit", 0, "time.deltaT");
        // line-implicit v2 試作 (plans/active/time_integration-line-implicit-viscous-v2.md):
        //   lineKFreeze: dual-time のサブ反復間で K/diag/LU 分解を凍結 (subiter 0 のみ抽出・分解)。
        //   lineViscCoupling: line 面にスカラー粘性結合 (K += α·I, 対角は 2α→α で真の [−α,2α,−α] 化)。
        //   lineViscousDtRelief: on-line セルの擬似 dt 粘性スペクトル半径を (1−θ) 倍に割引 (θ∈[0,1])。
        this->lineKFreeze = getOptionalValidatedValue<int>(deltaT, "lineKFreeze", 0, "time.deltaT");
        this->lineViscCoupling = getOptionalValidatedValue<int>(deltaT, "lineViscCoupling", 0, "time.deltaT");
        this->lineViscousDtRelief = getOptionalValidatedValue<double>(deltaT, "lineViscousDtRelief", 0.0, "time.deltaT");
        this->lineDtDirectional = getOptionalValidatedValue<int>(deltaT, "lineDtDirectional", 0, "time.deltaT");
        this->lineDtWallRelief = getOptionalValidatedValue<int>(deltaT, "lineDtWallRelief", 0, "time.deltaT");
        {
            double raw = getOptionalValidatedValue<double>(deltaT, "implicitRelaxSST", -1.0, "time.deltaT");
            this->implicitRelaxSST = (raw < 0.0) ? this->implicitRelax : (flow_float)raw;
        }
        this->blockDPLUR = getOptionalValidatedValue<int>(deltaT, "blockDPLUR", 0, "time.deltaT");
        // 多成分 TP 陰解法の化学種更新方式: 既定 0 (従来 segregated 点陰的・ビット不変)。
        // 1 で緩和整合 scalar-DPLUR (流れ block と同一緩和。plan thermophysics-species-implicit-coupling.md)。
        this->speciesImplicitCoupling = getOptionalValidatedValue<int>(deltaT, "speciesImplicitCoupling", 0, "time.deltaT");
        // 多成分 face 整合再構成: 既定 0 (mixed-order・ビット不変)。1 で Y を ρ と同じ再構成し thermo/species 流束整合。
        this->speciesFaceReconstruction = getOptionalValidatedValue<int>(deltaT, "speciesFaceReconstruction", 0, "time.deltaT");
        // multispeciesRhoYCommonLimiter: opt-in 診断 (既定 0・ビット不変)。1 で ρ と全 species に共通 min リミタ。
        this->multispeciesRhoYCommonLimiter = getOptionalValidatedValue<int>(deltaT, "multispeciesRhoYCommonLimiter", 0, "time.deltaT");
        // 軸対称 near-axis 安定化係数 β_axis: 既定 0 (不変)。擬似時間スペクトル半径に λ_axis=β(|u_r|+c)A_planar を加える。
        this->axisTimestepBeta = getOptionalValidatedValue<flow_float>(deltaT, "axisTimestepBeta", 0.0, "time.deltaT");
        // block-DPLUR 線形 solve の内部精度: 既定 0 (float・従来高速)。1 で double 化 (軸対称近軸の根治用)。
        this->implicitSolvePrecision = getOptionalValidatedValue<int>(deltaT, "implicitSolvePrecision", 0, "time.deltaT");
        // Ducros リミタ 1 次化: 既定 0 (off・使わない; MUSCL 2 次のまま)。1 で衝撃近傍の強制 1 次化 (従来挙動)。
        this->ducrosLimiter = getOptionalValidatedValue<int>(deltaT, "ducrosLimiter", 0, "time.deltaT");
        // NaN 検知診断モード: 既定 0 で従来挙動 (検査なし・ビット不変)。1 で検査 (fused device フラグ)。
        // 歴史的に time.deltaT 配下だが、既存 run config の多く (case/38 等) がトップレベルに書いて
        // 黙って off になる事故があったため (2026-07-21 case/39 で発覚)、トップレベルも受理する
        // (トップレベル優先)。
        this->detectNaN = getOptionalValidatedValue<int>(deltaT, "detectNaN", 0, "time.deltaT");
        if (config["detectNaN"]) this->detectNaN = config["detectNaN"].as<int>();
        // detectNaN フラグの host 読み出し間隔 [step]: 既定 1 (毎ステップ)。大で per-step 同期を間引く。
        this->detectNaNInterval = getOptionalValidatedValue<int>(deltaT, "detectNaNInterval", 1, "time.deltaT");
        if (config["detectNaNInterval"]) this->detectNaNInterval = config["detectNaNInterval"].as<int>();
        if (this->detectNaNInterval < 1) this->detectNaNInterval = 1;
        // 残差 RMS の device バッファ flush 間隔 [step]: 既定 1 (毎ステップ書き出し=従来挙動)。
        // モニタリング出力 (残差 CSV flush + max cfl/dt console 出力) の共通間隔 [step]: 既定 1 (毎ステップ=従来)。
        // 大で per-step の host 同期 (残差 D2H・max cfl 読み出し) をまとめて間引く。dt 適応とは独立。
        this->monitorInterval = getOptionalValidatedValue<int>(deltaT, "monitorInterval", 1, "time.deltaT");
        if (this->monitorInterval < 1) this->monitorInterval = 1;
        // 低マッハ前処理 (Weiss-Smith): 既定 0 で従来挙動 (ビット不変)。
        // 0: off / 1: フラックス散逸前処理のみ (RHS, c') / 2: RHS(c')+LHS 完全前処理 /
        // 3: LHS 完全前処理のみ (RHS 散逸は c_hat のまま=収束解は 0 と一致・保存性不変)。
        this->lowMachPrecond = getOptionalValidatedValue<int>(deltaT, "lowMachPrecond", 0, "time.deltaT");
        if (this->lowMachPrecond < 0 || this->lowMachPrecond > 3) {
            throw std::runtime_error(
                "lowMachPrecond must be 0 (off), 1 (flux dissipation only), "
                "2 (RHS+LHS full precond), or 3 (LHS-only full precond).");
        }
        this->precondEps = getOptionalValidatedValue<double>(deltaT, "precondEps", 0.15, "time.deltaT");
        // Thornber 型低マッハ再構成補正: 既定 0 で従来挙動 (ビット不変)。lowMachPrecond と直交・併用可。
        this->lowMachThornber = getOptionalValidatedValue<int>(deltaT, "lowMachThornber", 0, "time.deltaT");
        // 旧 space.keepDissipation は廃止 (KEEP_d は常に純粋 KEEP)。既存 config に残っていても無視される。
        this->dt_max = getValidatedValue<double>(deltaT, "dt_max", "time.deltaT");
        this->dt_min = getValidatedValue<double>(deltaT, "dt_min", "time.deltaT");

        this->outStepInterval = getValidatedValue<int>(config["time"], "outStepInterval", "time");
        this->outStepStart = getValidatedValue<int>(config["time"], "outStepStart", "time");
        this->timeIntegration = getValidatedValue<int>(config["time"], "timeIntegration", "time");
        if (config["time"]["nInnerLoop"].IsDefined()) {
            throw std::runtime_error(
                "Key 'nInnerLoop' in 'time' is no longer supported. Rename it to 'nStepInner'."
            );
        }
        if (config["time"]["dualTime_InnerLoop"].IsDefined()) {
            throw std::runtime_error(
                "Key 'dualTime_InnerLoop' in 'time' is no longer supported. Rename it to 'nStepInner'."
            );
        }
        this->nStepInner = getOptionalValidatedValue<int>(config["time"], "nStepInner", 1, "time");
        // dual-time (非定常陰解法): 1 物理ステップあたりの擬似時間サブ反復数と BDF 次数。
        this->nSubIterDualTime = getOptionalValidatedValue<int>(config["time"], "nSubIterDualTime", 20, "time");
        this->bdfOrder = getOptionalValidatedValue<int>(config["time"], "bdfOrder", 2, "time");

        // 空間設定
        auto space = config["space"];
        this->convMethod = getValidatedValue<int>(space, "convMethod", "space");
        this->limiter = getValidatedValue<int>(space, "limiter", "space");
        if (this->limiter != 0 && this->limiter != 1 && this->limiter != 2 && this->limiter != -1) {
            throw std::runtime_error("Key 'limiter' in 'space' must be one of 0, 1, 2, or -1.");
        }
        // free-stream 保存用の基準静圧 (既定 0.0 = 従来挙動・ビット不変)
        this->pRef = getOptionalValidatedValue<double>(space, "pRef", 0.0, "space");

        // 移流項の基準差分 (U∞≠0 の動く一様流保存, KEEP CPG)。roRef>0 で有効 (既定 0.0 = off・ビット不変)。
        this->roRef = getOptionalValidatedValue<double>(space, "roRef", 0.0, "space");
        if (space["uRef"]) {
            std::vector<double> u;
            for (const auto& v : space["uRef"]) u.push_back(v.as<double>());
            if (u.size() != 3) {
                throw std::runtime_error("Key 'uRef' in 'space' must be a list of 3 floats [x, y, z].");
            }
            this->uRefX = u[0]; this->uRefY = u[1]; this->uRefZ = u[2];
            std::cout << "'uRef' in 'space': [" << u[0] << ", " << u[1] << ", " << u[2] << "]\n";
        }
        if (this->roRef > 0.0) {
            if (this->pRef <= 0.0) {
                throw std::runtime_error("'roRef' in 'space' requires 'pRef' > 0 (advective gauge uses e_inf = pRef/roRef).");
            }
            // node は境界半割面 (convectiveFlux_boundary_d) にも同じゲージを載せて対応済 (plan §8.5)。
            // 有効化条件 (KEEP && CPG) は wrapper が主ループ/境界で厳密に揃える。
        }

        // turbulence model
        auto turb = config["turbulence"];
        // 乱流 closure は turbulence.model (文字列) 1 キーで指定する
        // (plans/active/turbulence-config-model-key.md)。旧 4 キー (LESorRANS/LESmodel/RANSmodel/DESmode)
        // は廃止・即エラー (nStepInner 改名と同じ方針)。内部 int 群へのマッピングのみ行い、カーネルは不変。
        for (const char* old_key : {"LESorRANS", "LESmodel", "RANSmodel", "DESmode"}) {
            if (turb[old_key]) {
                throw std::runtime_error(std::string("Key '") + old_key +
                    "' in 'turbulence' is no longer supported. Use 'model' instead: "
                    "none | wale | sigma | sst | sst-ddes | sst-iddes.");
            }
        }
        const std::string turbModel = getValidatedValue<std::string>(turb, "model", "turbulence");
        if      (turbModel == "none")      { this->LESorRANS = 0; this->LESmodel = 0; this->RANSmodel = 0; this->DESmode = 0; }
        else if (turbModel == "wale")      { this->LESorRANS = 1; this->LESmodel = 1; this->RANSmodel = 0; this->DESmode = 0; }
        else if (turbModel == "sigma")     { this->LESorRANS = 1; this->LESmodel = 2; this->RANSmodel = 0; this->DESmode = 0; }
        else if (turbModel == "sst")       { this->LESorRANS = 2; this->LESmodel = 0; this->RANSmodel = 1; this->DESmode = 0; }
        else if (turbModel == "sst-ddes")  { this->LESorRANS = 2; this->LESmodel = 0; this->RANSmodel = 1; this->DESmode = 1; }
        else if (turbModel == "sst-iddes") { this->LESorRANS = 2; this->LESmodel = 0; this->RANSmodel = 1; this->DESmode = 2; }
        else {
            throw std::runtime_error("Key 'model' in 'turbulence' must be one of: "
                "none | wale | sigma | sst | sst-ddes | sst-iddes (got '" + turbModel + "').");
        }
        std::cout << "'model' in 'turbulence': " << turbModel << "\n";
        this->scalarDiffusion = getOptionalValidatedValue<int>(turb, "scalarDiffusion", 1, "turbulence");
        this->dilatationCorrection = getOptionalValidatedValue<int>(turb, "dilatationCorrection", 2, "turbulence");

        if (this->dilatationCorrection < 0 || this->dilatationCorrection > 2) {
            throw std::runtime_error("Key 'dilatationCorrection' in 'turbulence' must be one of 0, 1, or 2.");
        }
        // ω 交差拡散 Jacobian (点陰対角強化)。既定 0 = ビット不変。
        this->sstCrossDiffJac = getOptionalValidatedValue<int>(turb, "sstCrossDiffJac", 0, "turbulence");

        this->katoLaunder = getOptionalValidatedValue<int>(turb, "katoLaunder", 0, "turbulence");

        if (this->katoLaunder < 0 || this->katoLaunder > 1) {
            throw std::runtime_error("Key 'katoLaunder' in 'turbulence' must be 0 or 1.");
        }

        // SST 壁処理 (methods/turbulence §6.5): 0=low-Re 壁解像 (60ν/β₁y², 既定), 1=automatic (y⁺ 非依存)
        this->wallTreatmentSST = getOptionalValidatedValue<int>(turb, "wallTreatmentSST", 1, "turbulence"); // 既定 1 (automatic)

        if (this->wallTreatmentSST < 0 || this->wallTreatmentSST > 1) {
            throw std::runtime_error("Key 'wallTreatmentSST' in 'turbulence' must be 0 or 1.");
        }

        // SST 壁関数の熱的閉包 (methods/turbulence §6.5(f)): 断熱壁の壁面温度を Crocco 型回復温度
        // T_aw = T_rep + Pr^{1/3}U_t²/(2cp) (Taw_diag) にする。wallTreatmentSST==1 のときのみ有効。
        this->sstThermalWallFunction = getOptionalValidatedValue<int>(turb, "sstThermalWallFunction", 0, "turbulence");
        if (this->sstThermalWallFunction < 0 || this->sstThermalWallFunction > 3) {
            throw std::runtime_error("Key 'sstThermalWallFunction' in 'turbulence' must be 0 (off), 1 (output-only), 2 (experimental SU2 coupled), or 3 (experimental defect-flux).");
        }
        this->sstEnergyWallFunction = getOptionalValidatedValue<int>(turb, "sstEnergyWallFunction", 0, "turbulence");
        if (this->sstEnergyWallFunction < 0 || this->sstEnergyWallFunction > 1) {
            throw std::runtime_error("Key 'sstEnergyWallFunction' in 'turbulence' must be 0 or 1.");
        }

        this->nodeKwfDirichlet = getOptionalValidatedValue<int>(turb, "nodeKwfDirichlet", 1, "turbulence"); // 既定 1 (node k Dirichlet)
        this->nodeOmegaWfDirichlet = getOptionalValidatedValue<int>(turb, "nodeOmegaWfDirichlet", 0, "turbulence");
        if (this->nodeKwfDirichlet < 0 || this->nodeKwfDirichlet > 1) {
            throw std::runtime_error("Key 'nodeKwfDirichlet' in 'turbulence' must be 0 or 1.");
        }

        // SST 初期乱流 (IC に roK/roOmega が無いときの初期値, 既定 0 = 従来動作・ビット不変)
        this->kInit     = getOptionalValidatedValue<flow_float>(turb, "kInit", 0.0, "turbulence");
        this->omegaInit = getOptionalValidatedValue<flow_float>(turb, "omegaInit", 0.0, "turbulence");
        if (this->kInit < 0.0 || this->omegaInit < 0.0) {
            throw std::runtime_error("Keys 'kInit'/'omegaInit' in 'turbulence' must be >= 0.");
        }

        // SST-DES length scale 修正 (methods/turbulence §8): DESmode は turbulence.model
        // (sst-ddes/sst-iddes) から設定済み。旧 'DESmode' キーは上のガードで拒否される。
        this->C_DES_kw = getOptionalValidatedValue<flow_float>(turb, "C_DES_kw", 0.78, "turbulence");
        this->C_DES_ke = getOptionalValidatedValue<flow_float>(turb, "C_DES_ke", 0.61, "turbulence");
        if (this->C_DES_kw <= 0.0 || this->C_DES_ke <= 0.0) {
            throw std::runtime_error("Keys 'C_DES_kw'/'C_DES_ke' in 'turbulence' must be positive.");
        }

        // 乱流プラントル数 Pr_t (乱流熱伝導 k_t = cp*mu_t/Pr_t)。既定 0.9。
        this->turbulentPrandtl = getOptionalValidatedValue<flow_float>(turb, "turbulentPrandtl", 0.85, "turbulence");
        if (this->turbulentPrandtl <= 0.0) {
            throw std::runtime_error("Key 'turbulentPrandtl' in 'turbulence' must be positive.");
        }

        // WMLES 代数壁応力モデル (methods/turbulence §10)。有効化は bcondConfig の壁単位 wallModelLES。
        this->wmlesNewtonTol   = getOptionalValidatedValue<flow_float>(turb, "wmlesNewtonTol", 1.0e-6, "turbulence");
        this->wmlesNewtonMaxIt = getOptionalValidatedValue<int>(turb, "wmlesNewtonMaxIt", 20, "turbulence");
        this->wmlesPrt         = getOptionalValidatedValue<flow_float>(turb, "wmlesPrt", 0.9, "turbulence");
        if (this->wmlesNewtonTol <= 0.0 || this->wmlesNewtonMaxIt < 1 || this->wmlesPrt <= 0.0) {
            throw std::runtime_error("Keys 'wmlesNewtonTol'/'wmlesNewtonMaxIt'/'wmlesPrt' in 'turbulence' must be positive.");
        }

        if (this->LESorRANS < 0 || this->LESorRANS > 2) {
            throw std::runtime_error("Key 'LESorRANS' in 'turbulence' must be one of 0, 1, or 2.");
        }
        if (this->LESorRANS == 1 && this->LESmodel <= 0) {
            throw std::runtime_error("Key 'LESmodel' in 'turbulence' must be positive when LESorRANS == 1.");
        }
        if (this->LESorRANS == 2 && this->RANSmodel <= 0) {
            throw std::runtime_error("Key 'RANSmodel' in 'turbulence' must be set when LESorRANS == 2.");
        }
        if (this->RANSmodel < 0 || this->RANSmodel > 1) {
            throw std::runtime_error("Key 'RANSmodel' in 'turbulence' must be one of 0 or 1.");
        }

        // 物理プロパティ
        auto physProp = config["physProp"];
        this->isCompressible = getValidatedValue<int>(physProp, "isCompressible", "physProp");
        // 軸対称スイッチの正本は mesh ブロック (幾何/離散化の設定であり物性ではない)。
        // physProp 配下は後方互換の deprecated 読み (既存 case config 用, 警告のみ)。
        this->isAxisymmetric = 0;
        if (physProp["isAxisymmetric"]) {
            this->isAxisymmetric = physProp["isAxisymmetric"].as<int>();
            std::cout << "[config] deprecated: physProp.isAxisymmetric -> use mesh.isAxisymmetric" << std::endl;
        }
        if (physProp["axisymMethod"]) {
            this->axisymMethod = physProp["axisymMethod"].as<int>();
            std::cout << "[config] deprecated: physProp.axisymMethod -> use mesh.axisymMethod" << std::endl;
        }
        if (config["mesh"] && config["mesh"]["isAxisymmetric"]) {
            this->isAxisymmetric = config["mesh"]["isAxisymmetric"].as<int>();
        }
        if (config["mesh"] && config["mesh"]["axisymMethod"]) {
            this->axisymMethod = config["mesh"]["axisymMethod"].as<int>();
        }
        if (config["mesh"] && config["mesh"]["axisRFloor"]) {
            this->axisRFloor = config["mesh"]["axisRFloor"].as<double>();
        }
        if (this->axisymMethod < 0 || this->axisymMethod > 1) {
            std::cerr << "axisymMethod must be 0 (r-weight) or 1 (SU2 source)" << std::endl;
            exit(EXIT_FAILURE);
        }
        if (this->discretization == "node") {
            std::cout << "[config] node scheme: values at nodes, edge-midpoint reconstruction, LSQ gradient"
                      << (this->isAxisymmetric == 1 ? ", axis u_r Dirichlet triple" : "") << std::endl;
        }
        this->thermalMethod = getValidatedValue<int>(physProp, "thermalMethod", "physProp");
        this->viscMethod = getValidatedValue<int>(physProp, "viscMethod", "physProp");
        this->ro = getValidatedValue<double>(physProp, "ro", "physProp");
        this->visc = getValidatedValue<double>(physProp, "visc", "physProp");
        this->thermCond = getValidatedValue<double>(physProp, "thermCond", "physProp");
        this->thermCondMethod = getOptionalValidatedValue<int>(physProp, "thermCondMethod", 0, "physProp");
        this->prandtlLam = getOptionalValidatedValue<double>(physProp, "prandtlLam", 0.72, "physProp");
        this->cp = getValidatedValue<double>(physProp, "cp", "physProp");
        //double gamma = getValidatedValue<double>(physProp, "gamma", "physProp");
        this->gamma = getValidatedValue<double>(physProp, "gamma", "physProp");

        // EOS 正値化フロア (任意)。未指定なら従来ハードコード値 (pMin=1.0, roMin=1e-4, tMin=1e-4)。
        // 無次元・低圧ケースでは小さくしてフロアによる場の破壊を防ぐ。
        this->pMin  = getOptionalValidatedValue<double>(physProp, "pMin",  1.0,    "physProp");
        this->roMin = getOptionalValidatedValue<double>(physProp, "roMin", 1.0e-4, "physProp");
        this->tMin  = getOptionalValidatedValue<double>(physProp, "tMin",  1.0e-4, "physProp");

        // 多成分 thermally-perfect gas 設定 (任意, thermalMethod==2 で使用)
        if (physProp["species"]) {
            this->speciesNames.clear();
            for (const auto& sn : physProp["species"]) this->speciesNames.push_back(sn.as<std::string>());
            this->nSpecies = static_cast<int>(this->speciesNames.size());
        }

        if (this->isAxisymmetric == 1 &&
            (this->bodyForceX != 0.0 || this->bodyForceY != 0.0 || this->bodyForceZ != 0.0)) {
            throw std::runtime_error("'bodyForce' is not supported with isAxisymmetric=1.");
        }
        if (this->keepDissFdBlend == 1) {
            if (this->DESmode <= 0) {
                throw std::runtime_error("'keepDissFdBlend: 1' requires a DES model ('model: sst-ddes' or 'sst-iddes').");
            }
            if (this->solver != "KEEP" || this->keepDissType <= 0) {
                throw std::runtime_error("'keepDissFdBlend: 1' requires solver: KEEP with keepDissType > 0.");
            }
            if (this->keepDissCoeffMax < this->keepDissCoeff || this->keepDissCoeffMax > 1.0) {
                throw std::runtime_error("'keepDissCoeffMax' must be in [keepDissCoeff, 1.0].");
            }
            // fdblend×precond は演算子ブレンドで併用可 (2026-07-22): 実効散逸を
            // σ_min·D_precond + (σ_f−σ_min)·D_std に内分する (kernel の opBlend)。
            // 素朴な併用 (σ_f を precond にそのまま掛ける) は Δp 増強 (∝c²/Ur)×σ→1 が壁近傍で
            // 発散するため不可 (case/39 run_diag_fdb_*: σmax0.3 でも step2 NaN)。precond を外した
            // c' 退避構成は大波長圧力モードの成長不安定を許す (run_0013 の壁 ΔP~100kPa 汚染) ため、
            // fdblend 時も precond (演算子ブレンド) が推奨。
        }
        if (this->bodyForceCtrl == 1) {
            if (this->unsteady != 1) {
                throw std::runtime_error("'bodyForceCtrl: 1' requires 'time.unsteady: 1' (physical dt).");
            }
            if (this->bodyForceCtrlTarget <= 0.0) {
                throw std::runtime_error("'bodyForceCtrl: 1' requires positive 'bodyForceCtrlTarget' (volume-averaged rho*Ux).");
            }
        }
        if (physProp["speciesDBFile"])          this->speciesDBFile = physProp["speciesDBFile"].as<std::string>();
        if (physProp["speciesDiffusionMethod"]) this->speciesDiffusionMethod = physProp["speciesDiffusionMethod"].as<int>();
        if (physProp["thermoHrefTemp"])         this->thermoHrefTemp = physProp["thermoHrefTemp"].as<double>();
        if (physProp["chemistry"]) {
            const YAML::Node ch = physProp["chemistry"];
            this->chemEnabled       = getOptionalValidatedValue<int>(ch, "enabled", 0, "physProp.chemistry");
            if (ch["mechanismFile"]) this->chemMechanismFile = ch["mechanismFile"].as<std::string>();
            this->chemTmaxReaction  = getOptionalValidatedValue<double>(ch, "tMaxReaction", 6000.0, "physProp.chemistry");
            this->chemFreezeBelowT  = getOptionalValidatedValue<double>(ch, "freezeBelowT", 0.0, "physProp.chemistry");
            this->chemJacobianMode  = getOptionalValidatedValue<int>(ch, "jacobianMode", 1, "physProp.chemistry");
            if (this->chemEnabled != 0) {
                if (this->thermalMethod != 2) throw std::runtime_error("'physProp.chemistry.enabled: 1' requires thermalMethod: 2.");
                if (this->nSpecies < 2)      throw std::runtime_error("'physProp.chemistry.enabled: 1' requires >=2 species.");
                if (this->chemMechanismFile.empty()) throw std::runtime_error("'physProp.chemistry.mechanismFile' is required.");
            }
        }
        if (physProp["Sc"])                     this->Sc = physProp["Sc"].as<double>();
        if (physProp["Sc_t"])                   this->Sc_t = physProp["Sc_t"].as<double>();
        // 乱流シュミット数は turbulence.turbulentSchmidt でも設定可 (turbulentPrandtl と同じ場所)。
        // 後方互換: physProp.Sc_t も有効。両方あれば turbulence 側を優先する。
        if (config["turbulence"] && config["turbulence"]["turbulentSchmidt"]) {
            this->Sc_t = config["turbulence"]["turbulentSchmidt"].as<flow_float>();
            std::cout << "'turbulentSchmidt' in 'turbulence': " << this->Sc_t << std::endl;
        }
        if (this->Sc_t <= 0.0) {
            throw std::runtime_error("Turbulent Schmidt number (Sc_t / turbulentSchmidt) must be positive.");
        }
        if (this->thermalMethod == 2 && this->speciesNames.empty()) {
            // 既定: 単成分 N2 (CEA 検証用)
            this->speciesNames.push_back("N2");
            this->nSpecies = 1;
        }

        // 非平衡凝縮 (任意セクション)。methods/condensation/ 参照。
        // Phase 1 は受動スカラー輸送のみ (核生成/成長係数はまだ読まない)。
        if (config["condensation"]) {
            auto cond = config["condensation"];
            this->condensation = getOptionalValidatedValue<int>(cond, "condensation", 0, "condensation");
            this->nCondSpecies = getOptionalValidatedValue<int>(cond, "nCondSpecies", 0, "condensation");
            this->condModel = getOptionalValidatedValue<int>(cond, "condModel", 0, "condensation");
            this->condGasSpecies = getOptionalValidatedValue<int>(cond, "condGasSpecies", -1, "condensation");
            this->condKantrowitz = getOptionalValidatedValue<int>(cond, "condKantrowitz", 0, "condensation");
            this->condGrowthModel = getOptionalValidatedValue<int>(cond, "condGrowthModel", 0, "condensation");
            this->condGyarmathyC = getOptionalValidatedValue<double>(cond, "condGyarmathyC", 3.18, "condensation");
            this->condTwoTemp = getOptionalValidatedValue<int>(cond, "condTwoTemp", 0, "condensation");
            this->condEvaporation = getOptionalValidatedValue<int>(cond, "condEvaporation", 1, "condensation");
            this->condEvapRmin    = getOptionalValidatedValue<double>(cond, "condEvapRmin", 1.0e-9, "condensation");
            this->condEvapKelvin  = getOptionalValidatedValue<int>(cond, "condEvapKelvin", 0, "condensation");
            this->condEquilibrium = getOptionalValidatedValue<int>(cond, "condEquilibrium", 0, "condensation");
            this->condEqRelax     = getOptionalValidatedValue<double>(cond, "condEqRelax", 1.0, "condensation");
            this->condEqDTmax     = getOptionalValidatedValue<double>(cond, "condEqDTmax", 10.0, "condensation");
            this->condEqDgMax     = getOptionalValidatedValue<double>(cond, "condEqDgMax", 0.05, "condensation");
        }
        if (this->condensation != 0 && this->condensation != 1) {
            throw std::runtime_error("Key 'condensation' in 'condensation' must be 0 or 1.");
        }
        if (this->condensation == 1 && this->nCondSpecies < 1) {
            throw std::runtime_error("Key 'nCondSpecies' in 'condensation' must be >= 1 when condensation == 1.");
        }
        if (this->condEquilibrium < 0 || this->condEquilibrium > 2) {
            throw std::runtime_error("Key 'condEquilibrium' in 'condensation' must be 0, 1 or 2.");
        }
        if (this->condensation == 1 && this->condEquilibrium == 2 && this->nCondSpecies != 1) {
            throw std::runtime_error("condEquilibrium == 2 (EOS-constrained equilibrium) supports nCondSpecies == 1 only.");
        }
        if (this->condensation == 0) {
            this->nCondSpecies = 0; // off のときは登録しない
        }

        // 初期条件
        this->initial = getValidatedValue<std::string>(config, "initial");



        // デバッグ用出力
        std::cout << "Mesh Format: " << meshFormat << '\n';
        std::cout << "Mesh File Name: " << meshFileName << '\n';
        std::cout << "Solver: " << solver << '\n';

    } catch (const std::runtime_error& e) {
        std::cerr << "Configuration Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    } catch (const YAML::Exception& e) {
        std::cerr << "YAML Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    initTimeIntegrationScheme(timeIntegration);
}

void solverConfig::initTimeIntegrationScheme(int timeIntegration){

    int ns;
    switch (timeIntegration) {
        case 1: // 1st Euler Explicit
            ns = 1;
            this->isImplicit = 0;
            this->nStage = ns;
            this->coef_N.resize(ns);
            this->coef_M.resize(ns);
            this->coef_Res.resize(ns);

            coef_N[0]   = 1.0;
            coef_M[0]   = 0.0;
            coef_Res[0] = 1.0;
            break;
        case 3: // 3rd TVD Runge-Kutta Explicit
            ns = 3;
            this->isImplicit = 0;
            this->nStage = ns;
            this->coef_N.resize(ns);
            this->coef_M.resize(ns);
            this->coef_Res.resize(ns);

            coef_N[0]   = 1.0 ; coef_N[1]   = 0.75 ; coef_N[2]   = 1.0/3.0;
            coef_M[0]   = 0.0 ; coef_M[1]   = 0.25 ; coef_M[2]   = 2.0/3.0;
            coef_Res[0] = 1.0 ; coef_Res[1] = 0.25 ; coef_Res[2] = 2.0/3.0;
            break;

        case 4: // 4th Runge-Kutta Explicit
            ns = 4;
            this->isImplicit = 0;
            this->nStage = ns;
            this->coef_DT_4thRunge.resize(ns);
            this->coef_Res_4thRunge.resize(ns);

            this->coef_DT_4thRunge[0] = 0.5;
            this->coef_DT_4thRunge[1] = 0.5;
            this->coef_DT_4thRunge[2] = 1.0;
            this->coef_DT_4thRunge[3] = -1e+30;

            this->coef_Res_4thRunge[0] = 1.0/6.0;
            this->coef_Res_4thRunge[1] = 1.0/3.0;
            this->coef_Res_4thRunge[2] = 1.0/3.0;
            this->coef_Res_4thRunge[3] = 1.0/6.0;

            break;

        case 11: // implicit pseudo-time defect-correction (block DPLUR)
            this->isImplicit = 1;
            this->nStage = 1;
            this->coef_N.clear();
            this->coef_M.clear();
            this->coef_Res.clear();
            this->coef_DT_4thRunge.clear();
            this->coef_Res_4thRunge.clear();
            // blockDPLUR==1: 5×5 block (LU-SGS A±)、blockDPLUR==0: scalar 対角 (スペクトル半径) 版。
            // どちらも古典 DPLUR 制御フロー (残差固定 + nStepInner sweep + 単一 commit) に対応。
            if (this->blockDPLUR != 0 && this->blockDPLUR != 1) {
                throw std::runtime_error(
                    "timeIntegration==11: blockDPLUR must be 0 (scalar-diagonal) or 1 (5x5 block).");
            }
            break;

        default:
            throw std::runtime_error("Unsupported timeIntegration: " + std::to_string(timeIntegration));

    }

    // 低マッハ完全前処理 (lowMachPrecond>=2: 2=RHS+LHS / 3=LHS-only) は block-DPLUR 陰解法
    // (timeIntegration==11 && blockDPLUR==1) でのみ LHS 前処理カーネルが起動する。これ以外の経路では
    // setDT の擬似時間刻み拡大だけが効いて LHS 前処理が黙って抜け落ちる事故になるため明示エラーで停止する。
    // (lowMachPrecond==1 はフラックス散逸のみなので経路非依存・本チェック対象外。)
    if (this->lowMachPrecond >= 2) {
        if (timeIntegration != 11 || this->blockDPLUR != 1) {
            throw std::runtime_error(
                "lowMachPrecond>=2 (full preconditioned block DPLUR: 2=RHS+LHS, 3=LHS-only) "
                "requires timeIntegration==11 (implicit) and blockDPLUR==1.");
        }
    }

    // implicitSolvePrecision (block-DPLUR 線形 solve の double 化) のサポート範囲チェック。
    // 未対応の組み合わせはフラグが黙って無視され「効いていない」事故になるため、明示エラーで停止する。
    // 実装済みは block 経路 (blockDPLUR==1) かつ非 precond (lowMachPrecond 0/1) の implicit (timeIntegration==11) のみ。
    // 詳細: plans/archived/precision-mixed-axisym.md §13。
    if (this->implicitSolvePrecision != 0 && this->implicitSolvePrecision != 1) {
        throw std::runtime_error(
            "implicitSolvePrecision must be 0 (float) or 1 (double).");
    }
    if (this->implicitSolvePrecision == 1) {
        if (timeIntegration != 11) {
            throw std::runtime_error(
                "implicitSolvePrecision==1 (double solve) requires timeIntegration==11 "
                "(implicit block DPLUR). It has no effect for explicit schemes.");
        }
        if (this->blockDPLUR != 1) {
            throw std::runtime_error(
                "implicitSolvePrecision==1 (double solve) is only supported with blockDPLUR==1 "
                "(5x5 block DPLUR); scalar-diagonal DPLUR (blockDPLUR==0) is not supported.");
        }
        if (this->lowMachPrecond >= 2) {
            throw std::runtime_error(
                "implicitSolvePrecision==1 (double solve) is not supported with lowMachPrecond>=2 "
                "(preconditioned block DPLUR: 2=RHS+LHS, 3=LHS-only). Use lowMachPrecond 0 or 1.");
        }
    }

}

int solverConfig::mainLoopCount() const
{
    return std::max(1, this->nStepOuter);
}

int solverConfig::perStepIterationCount() const
{
    if (this->isImplicit == 1) {
        return std::max(1, this->nStepInner);
    }

    return std::max(1, this->nStage);
}

const char* solverConfig::perStepIterationLabel() const
{
    if (this->isImplicit == 1) {
        return "Inner";
    }

    return "Stage";
}
