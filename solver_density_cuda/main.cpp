#include <algorithm>
#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <stdio.h>                                                                                       
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <limits>
#include <cstdio>
#include <map>
#include <tuple>
#include <cstring>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

#include "input/solverConfig.hpp"
#include "input/condSonicResolve.hpp"
#include "input/setInitial.hpp"

#include "mesh/gmshReader.hpp"
#include "mesh/memlog.hpp"
#include "output/output.hpp"
#include "output/outputFieldNames.hpp"
#include "conjugateWall.hpp"
#include "boundaryCond.hpp"

#include "setStructualVariables.hpp"

#include "gradient.hpp"

#include "dependentVariables.hpp"
//#include "solvePoisson_amgx.hpp"

#include "convectiveFlux.hpp"
//#include "timeIntegration.hpp"
#include "update.hpp"

#include "common/stringUtil.hpp"
#include "common/vectorUtil.hpp"

//#include "setDT.hpp"

// cuda
#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/cudaConfig.cuh"

#include "cuda_forge/calcGradient_d.cuh"
#include "cuda_forge/convection/convectiveFlux_d.cuh"
#include "cuda_forge/ransTransport_d.cuh"
#include "cuda_forge/ransSource_d.cuh"
#include "cuda_forge/transition_d.cuh"
#include "cuda_forge/speciesTransport_d.cuh"
#include "cuda_forge/chemistrySource_d.cuh"
#include "cuda_forge/condensationTransport_d.cuh"
#include "cuda_forge/twoPhaseUpdateDiag_d.cuh"   // 診断 G3-b の記録スロットの配置 (FORGE_DIAG_TP_UPDATE; 名前表と定数だけ)
#include "cuda_forge/twoPhaseOperatorDiag_d.cuh" // 診断 G3-a の記録の配置 (FORGE_DIAG_TP_OPERATOR; 名前表と定数だけ)
#include "cuda_forge/tracerTransport_d.cuh"
#include "cuda_forge/passiveTransport_d.cuh"
#include "cuda_forge/gasPhaseComposition_d.cuh"   // gasPhaseLiquid (transport probe の液)
#include <highfive/H5File.hpp>
#include "input/speciesDB.hpp"
#include "input/zeroThicknessEdge.hpp"   // space.zeroThicknessEdgeVelocity の起動時検査・読込
#include "cuda_forge/viscousFlux_d.cuh"
#include "cuda_forge/updateCenterVelocity_d.cuh"
#include "cuda_forge/interpVelocity_c2p_d.cuh"
#include "cuda_forge/timeIntegration_d.cuh"
#include "cuda_forge/limiter_d.cuh"
#include "cuda_forge/ducrosSensor_d.cuh"
#include "cuda_forge/axisymmetricSource_d.cuh"
#include "cuda_forge/wmlesWallModel_d.cuh"
#include "cuda_forge/nodeWallDirichlet_d.cuh"
#include "cuda_forge/weakIsothermalWall_d.cuh"
#include "cuda_forge/bodyForce_d.cuh"
#include "cuda_forge/turbulent_viscosity_d.cuh"
#include "cuda_forge/residualMonitor_d.cuh"

#include "cuda_forge/fluct_variables_d.cuh"
#include "cuda_forge/gasProperties_d.cuh"
#include "cuda_forge/transportTables_d.cuh"   // 表引きの試験ハーネス (FORGE_TRANSPORT_TABLE_PROBE; #5t2-3)
#include "cuda_forge/thermo_d.cuh"

#include "probe/point_probes.cuh"
#include "cuda_forge/setDT_d.cuh"
#include "cuda_forge/periodicNode_d.cuh"
#include "cuda_forge/qAccumulator.hpp"
#include "cuda_forge/floorEvents_d.cuh"   // 毎更新の EOS 床事象のカウンタ (output.floorEvents、出力専用)

#include <cuda_runtime.h>
#include <sys/stat.h>
#include <ctime>
#include <fstream>
#include <sstream>

// 起動ごとの実効設定の記録 (plan convection-slau-wall-normal-chi-default §4.3、codex plan M1)。
// run_case.sh は forge_run.log / RUN_PROVENANCE.txt を起動ごとに上書きし、段階起動の runner は残差しか退避しないので、
// 段ごとの実効値は forge 自身が **追記** で残す。stage_manifest.py が cfg_fnv (solverConfig.yaml の FNV-1a 64) で段と結び付ける。
static unsigned long long fnv1a64File(const std::string& path)
{
    std::ifstream f(path, std::ios::binary);
    unsigned long long h = 14695981039346656037ULL;
    char buf[1 << 16];
    while (f) {
        f.read(buf, sizeof(buf));
        for (std::streamsize i = 0; i < f.gcount(); ++i) { h ^= (unsigned char)buf[i]; h *= 1099511628211ULL; }
    }
    return h;
}

// 起動記録の本体 (閉じ括弧の前まで)。厚さ 0 の板の自由端の処置の記録 (appendLaunchRecordZte) が同じ値で 2 行目を書くために残す。
static std::string g_launchRecordBody;

static void appendLaunchRecord(const solverConfig& cfg)
{
    struct stat st{};
    long long exeSize = 0, exeMtime = 0;
    if (stat("/proc/self/exe", &st) == 0) { exeSize = (long long)st.st_size; exeMtime = (long long)st.st_mtime; }
    std::ostringstream os;
    os << "{\"time\": " << (long long)std::time(nullptr)
       << ", \"cfg_fnv\": \"" << std::hex << fnv1a64File("solverConfig.yaml") << "\""
       << ", \"bcond_fnv\": \"" << fnv1a64File("bcondConfig.yaml") << std::dec << "\""
       << ", \"exe_size\": " << exeSize << ", \"exe_mtime\": " << exeMtime
       << ", \"slauWallNormalChi\": " << cfg.slauWallNormalChi
       << ", \"slauWallNormalChi_source\": \"" << cfg.slauWallNormalChiReason << "\""
       << ", \"scalarGradient\": \"" << cfg.scalarGradient << "\""
       << ", \"scalarGradient_source\": \"" << cfg.scalarGradientReason << "\"";
    g_launchRecordBody = os.str();
    std::ofstream out("forge_launches.jsonl", std::ios::app);
    if (out) out << g_launchRecordBody << "}" << "\n";
}

// 厚さ 0 の板の自由端の速度再構成の処置 (space.zeroThicknessEdgeVelocity) の起動記録。appendLaunchRecord は格子を読む前に
// 走るので、格子署名と w の照合が済んだ後に**同じ起動の 2 行目**として追記する (1 行目の形式は変えない)。2 行目は 1 行目と
// 同じキー・同じ値 (最後の行を読む側が slauWallNormalChi 等を失わないため) に "record" と入れ子の処置の記録を足したもの。
// 値は属性の転記でなく、ソルバが読み込んだ格子と w から再計算して照合した値。無効は {"enabled": 0} だけ。
static void appendLaunchRecordZte(const ZteLaunchInfo& z)
{
    if (g_launchRecordBody.empty()) return;
    std::ostringstream os;
    os << g_launchRecordBody << ", \"record\": \"zeroThicknessEdgeVelocity\", \"zeroThicknessEdgeVelocity\": {\"enabled\": " << z.enabled;
    if (z.enabled) {
        os << ", \"field\": \"" << z.field << "\", \"mesh_signature\": \"" << z.meshSignature
           << "\", \"mesh_signature_version\": \"" << z.meshSignatureVersion << "\", \"field_hash\": \"" << z.fieldSha256
           << "\", \"n_w_zero\": " << z.nZero << ", \"n_w_lt1\": " << z.nLt1 << ", \"n_nodes\": " << z.nNodes;
    }
    os << "}}";
    std::ofstream out("forge_launches.jsonl", std::ios::app);
    if (out) out << os.str() << "\n";
}

// FORGE_MEMLOG=1 の工程別メモリ計測 (mesh/memlog.hpp) で forge 本体の主要コンテナと GPU 使用量を 1 行にまとめる。
// 計測専用で挙動・出力には影響しない (MEMLOG マクロの extra 引数としてだけ評価されるので、既定では呼ばれもしない)。
// ホスト側は size/capacity × sizeof と glibc のチャンク概算 (memlog::heapChunk)。GPU 側は cudaMemGetInfo の
// (total − free) と、最初の呼び出し時点 (= CUDA context 作成直後) からの増分、c_d/p_d の確保済み本数からの推定。
static std::string memSolverSummary(const mesh& msh, const variables& var, const matrix* mat)
{
    using namespace memlog;
    std::ostringstream os;
    // var.c / var.p (ホストの cell / plane 変数): 確保済み (capacity > 0) の本数と合計
    size_t nc = 0, bc = 0, np = 0, bp = 0;
    for (const auto& kv : var.c) if (kv.second.capacity() > 0) { ++nc; bc += flatBytes(kv.second); }
    for (const auto& kv : var.p) if (kv.second.capacity() > 0) { ++np; bp += flatBytes(kv.second); }
    os << item("var.c(host)", nc, bc) << " " << item("var.p(host)", np, bp);
    // msh の格子構造 (vector-of-structs、要素ごとに内側 vector のヒープ確保を持つ)
    auto nodeIn  = [](const node& n)  { return innerVec(n.iCells) + innerVec(n.iPlanes) + innerVec(n.coords); };
    auto planeIn = [](const plane& p) { return innerVec(p.iNodes) + innerVec(p.iCells) + innerVec(p.surfVect) + innerVec(p.centCoords); };
    auto cellIn  = [](const cell& c)  { return innerVec(c.iNodes) + innerVec(c.iPlanes) + innerVec(c.iPlanesDir) + innerVec(c.centCoords); };
    os << " " << item("msh.nodes", msh.nodes.size(), nestedBytes(msh.nodes, nodeIn))
       << " " << item("msh.planes", msh.planes.size(), nestedBytes(msh.planes, planeIn))
       << " " << item("msh.cells", msh.cells.size(), nestedBytes(msh.cells, cellIn));
    if (!msh.vizCONNE.empty()) os << " " << item("msh.vizCONNE", msh.vizCONNE.size(), flatBytes(msh.vizCONNE));
    // 境界: 面 index 列・可視化面・bvar/bint (host)・出力用の局所格子
    size_t bIdx = 0, bViz = 0, bVar = 0, bLoc = 0, nBvar = 0;
    for (const auto& b : msh.bconds) {
        bIdx += flatBytes(b.iPlanes) + flatBytes(b.iBPlanes) + flatBytes(b.iCells) + flatBytes(b.iCells_ghst);
        bViz += nestedBytes(b.vizBfaceNodes, [](const std::vector<geom_int>& f) { return innerVec(f); });
        for (const auto& kv : b.bvar) { bVar += flatBytes(kv.second); ++nBvar; }
        for (const auto& kv : b.bint) bVar += flatBytes(kv.second);
        bLoc += nestedBytes(b.nodes_local, nodeIn) + nestedBytes(b.planes_local, planeIn)
              + flatBytes(b.inodes_l2g) + flatBytes(b.inodes_g2l);
    }
    os << " " << item("bconds.idx", msh.bconds.size(), bIdx) << " " << item("bconds.vizBface", msh.bconds.size(), bViz)
       << " " << item("bconds.bvar+bint(host)", nBvar, bVar) << " " << item("bconds.local(output)", msh.bconds.size(), bLoc);
    // 計算に使っていない行列 (initializeSimulation の "Init Matrix (but not used now)")
    if (mat != nullptr) {
        const size_t bMat = nestedBytes(mat->structure, [](const std::vector<geom_int>& v) { return innerVec(v); })
                          + nestedBytes(mat->lhs, [](const std::vector<flow_float>& v) { return innerVec(v); })
                          + flatBytes(mat->rhs)
                          + nestedBytes(mat->localPlnOfCell, [](const std::vector<geom_int>& v) { return innerVec(v); });
        os << " " << item("mat_ns", mat->structure.size(), bMat);
    }
    if (!msh.periodicRoot.empty() || !msh.rEff.empty())
        os << " " << item("msh.periodicRoot+rEff", msh.periodicRoot.size(), flatBytes(msh.periodicRoot) + flatBytes(msh.rEff));
    // GPU: c_d / p_d の確保済み本数からの推定と、cudaMemGetInfo の実使用量 (デバイス全体。他プロセス分を含む)
    size_t ncd = 0, npd = 0;
    for (const auto& kv : var.c_d) if (kv.second != nullptr) ++ncd;
    for (const auto& kv : var.p_d) if (kv.second != nullptr) ++npd;
    const size_t bcd = ncd * (size_t)msh.nCells_all * sizeof(flow_float), bpd = npd * (size_t)msh.nPlanes * sizeof(flow_float);
    os << " | " << item("var.c_d(est)", ncd, bcd) << " " << item("var.p_d(est)", npd, bpd);
    size_t fr = 0, tot = 0;
    if (cudaMemGetInfo(&fr, &tot) == cudaSuccess) {
        static const size_t used0 = tot - fr;   // 最初の呼び出し (context 作成直後) の使用量
        const size_t used = tot - fr;
        os << " GPU used=" << used / 1048576 << "MB (since first memlog +" << ((long long)used - (long long)used0) / 1048576 << "MB)";
    }
    return os.str();
}
// 工程の境目に置く 1 行 (msh / var / mat_ns がスコープにある所で使う)
#define MEMLOG_SOLVER(label) MEMLOG((label), memSolverSummary(msh, var, &mat_ns))

#define CHECK_LAST_CUDA_ERROR() checkLast(__FILE__, __LINE__)
void checkLast(const char* const file, const int line)
{
    cudaError_t err{cudaGetLastError()};
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA Runtime Error at: " << file << ":" << line
                  << std::endl;
        std::cerr << cudaGetErrorString(err) << std::endl;
        // We don't exit when we encounter CUDA errors in this example.
        std::exit(EXIT_FAILURE);
    }
}

namespace {

enum class ProfileSection {
    UpdateInner,
    DependentVariables,
    GasProperties,
    ApplyBconds,
    CalcGradient,
    Limiter,
    DucrosSensor,
    TurbulenceModel,
    ConvectiveFlux,
    AxisymmetricSource,
    ViscousFlux,
    ImplicitCorrection,
    TimeIntegration,
    UpdateOuter,
    WriteOutputs,
    SetDt,
    StepTotal,
    Count
};

constexpr std::size_t kProfileSectionCount = static_cast<std::size_t>(ProfileSection::Count);
enum class ResidualPhase {
    OuterBegin,
    InnerBegin,
    InnerIter,
    OuterEnd,
};

constexpr std::array<const char*, 5> kResidualEquationNames = {
    "ro", "roUx", "roUy", "roUz", "roe"
};

constexpr std::array<const char*, 2> kScalarResidualEquationNames = {
    "roK", "roOmega"
};

bool scalarResidualEnabled(const solverConfig& cfg)
{
    return cfg.LESorRANS == 2 && cfg.RANSmodel == 1;
}

// 非平衡凝縮 (Phase 1): 4 モーメントの輸送が有効か。
bool condensationEnabled(const solverConfig& cfg)
{
    return cfg.condensation == 1 && cfg.nCondSpecies >= 1;
}

// 凝縮モーメントの保存量名 (登録順)。variables::registerCondensation と同じ命名規約で cfg のみから導出。
// 1 凝縮種あたり 4 本: rog_{s}, roQ2_{s}, roQ1_{s}, roQ0_{s}。
std::vector<std::string> condMomentConsNames(const solverConfig& cfg)
{
    std::vector<std::string> names;
    if (!condensationEnabled(cfg)) return names;
    const std::array<const char*, 4> bases = {"g", "Q2", "Q1", "Q0"};
    for (int s = 0; s < cfg.nCondSpecies; ++s) {
        for (const auto* b : bases) {
            names.emplace_back(std::string("ro") + b + "_" + std::to_string(s));
        }
    }
    return names;
}

// 化学種の保存量名 (nSpecies>=2 のとき roY0..roY{n-1}; variables::registerSpecies と同じ命名規約)。
// 残差配列は res_roY{s} (speciesTransport_d.cu)。単成分では空。
std::vector<std::string> speciesResidualNames(const solverConfig& cfg)
{
    std::vector<std::string> names;
    if (cfg.nSpecies < 2) return names;
    for (int s = 0; s < cfg.nSpecies; ++s) names.emplace_back("roY" + std::to_string(s));
    return names;
}

std::vector<std::string> residualEquationNames(const solverConfig& cfg)
{
    std::vector<std::string> names;
    names.reserve(kResidualEquationNames.size() + kScalarResidualEquationNames.size());

    for (const auto* name : kResidualEquationNames) {
        names.emplace_back(name);
    }

    if (scalarResidualEnabled(cfg)) {
        for (const auto* name : kScalarResidualEquationNames) {
            names.emplace_back(name);
        }
    }
    // 遷移モデル rms_roGamma / rms_roReth (SST の直後)。
    if (cfg.transitionEnabled()) {
        names.emplace_back("roGamma");
        names.emplace_back("roReth");
    }

    // 化学種 rms_roY{s} (流れ/RANS の後、凝縮モーメントの前)。受動トレーサ rms_roXi はその次。
    for (const auto& name : speciesResidualNames(cfg)) {
        names.emplace_back(name);
    }
    if (cfg.tracerEnabled()) {
        names.emplace_back("roXi");
    }

    if (condensationEnabled(cfg)) {
        for (const auto& name : condMomentConsNames(cfg)) {
            names.emplace_back(name);
        }
    }
    // 二相拡散 (#4e) の蒸気残差 rms_roYv = rms(res_roY_w − res_rog_0) (新キー ON の構成だけ; 既定の列構成は不変)。
    if (condTwoPhaseDiffusionActive(cfg)) {
        names.emplace_back("roYv");
    }

    return names;
}

struct ResidualSnapshot {
    std::vector<flow_float> rms;
};

struct CorrectionSnapshot {
    std::array<flow_float, kResidualEquationNames.size()> rms{};
};

struct ImplicitDiagSnapshot {
    double mean_pseudo_diag = 0.0;
    double mean_face_diag = 0.0;
    double mean_pseudo_ratio = 0.0;
    double min_pseudo_ratio = 0.0;
    double max_pseudo_ratio = 0.0;
    double mean_face_to_pseudo = 0.0;
    double mean_dt_local = 0.0;
};

// 残差 RMS を CSV (residual_history.csv) へ書き出すロガー。
// GPU 経路 (cfg.gpu==1) では残差二乗和を device 常駐バッファに async 集約し、monitorInterval
// ステップごとに 1 回だけ D2H 転送してまとめて書き出す (per-step host 同期を避ける)。CPU 経路は
// 従来通り即時計算・即時書き出し。CSV の列・行構成・値 (毎ステップ・3 行/step・rms_dq_*=0) は不変。
// 設計: plans/accepted/architecture-residual-monitor-async.md
class ResidualCsvLogger {
public:
    ResidualCsvLogger(const std::string& file_name, const solverConfig& cfg, mesh& msh, variables& var)
        : residual_names_(residualEquationNames(cfg)),
          stream_(file_name, std::ios::out | std::ios::trunc),
          gpu_(cfg.gpu == 1)
    {
        if (!stream_) {
            throw std::runtime_error("Failed to open residual log file: " + file_name);
        }

        stream_ << "step,inner,phase";
        for (const auto& name : residual_names_) {
            stream_ << ",rms_" << name;
        }
        for (const auto* name : kResidualEquationNames) {
            stream_ << ",rms_dq_" << name;
        }
        stream_ << "\n";

        nVar_ = static_cast<int>(residual_names_.size());
        nDq_ = static_cast<int>(kResidualEquationNames.size());
        last_snapshot_.rms.assign(residual_names_.size(), static_cast<flow_float>(0.0));

        if (gpu_) {
            // 残差配列名は "res_" + 方程式名 (gatherResidualSnapshot と同じ)。
            std::vector<std::string> dev_names;
            dev_names.reserve(residual_names_.size());
            for (const auto& name : residual_names_) {
                dev_names.emplace_back("res_" + name);
            }
            std::vector<std::string> resolved;
            reducer_ = makeDeviceResidualReducer(msh, var, dev_names, resolved);
            if (reducer_.nVar != nVar_) {
                throw std::runtime_error(
                    "ResidualCsvLogger: device residual arrays (res_*) do not match active equations; "
                    "found " + std::to_string(reducer_.nVar) + " of " + std::to_string(nVar_));
            }
            capacity_ = std::max(1, cfg.monitorInterval);
            // 1 step に複数回 gather する経路 (RK/dual-time) の余裕として margin を足す。
            buf_capacity_ = capacity_ + kIntraStepMargin;
            d_rms_buf_ = allocDeviceRmsBuffer(buf_capacity_, reducer_.nVar);
            host_rms_.resize(static_cast<size_t>(buf_capacity_) * reducer_.nVar);
            rows_.reserve(static_cast<size_t>(buf_capacity_) * 3);
        }
    }

    ~ResidualCsvLogger()
    {
        if (gpu_) {
            flushDevice();
            freeDeviceRmsBuffer(d_rms_buf_);
            freeDeviceResidualReducer(reducer_);
        }
    }

    // ---- GPU 経路: device へ async reduce し行を buffer する (同期しない) ----
    void recordGpu(int step, int inner_index, mesh& /*msh*/, variables& /*var*/)
    {
        const int slot = device_count_++;
        reduceResidualToSlot(reducer_, d_rms_buf_, slot);
        last_slot_ = slot;
        if (inner_index == 0) {
            rows_.push_back({step, -1, ResidualPhase::OuterBegin, slot});
            rows_.push_back({step, inner_index, ResidualPhase::InnerBegin, slot});
        } else {
            rows_.push_back({step, inner_index, ResidualPhase::InnerIter, slot});
        }
    }

    // ---- CPU 経路: 即時書き出し (従来挙動) ----
    void logOuterBegin(
        int step,
        int inner,
        const ResidualSnapshot& snapshot,
        const CorrectionSnapshot& correction_snapshot = {})
    {
        last_snapshot_ = snapshot;
        last_correction_snapshot_ = correction_snapshot;

        writeRowImmediate(step, -1, ResidualPhase::OuterBegin, snapshot, correction_snapshot);
        writeRowImmediate(step, inner, ResidualPhase::InnerBegin, snapshot, correction_snapshot);
    }

    void logInnerIter(
        int step,
        int inner_index,
        const ResidualSnapshot& snapshot,
        const CorrectionSnapshot& correction_snapshot = {})
    {
        last_snapshot_ = snapshot;
        last_correction_snapshot_ = correction_snapshot;
        writeRowImmediate(step, inner_index, ResidualPhase::InnerIter, snapshot, correction_snapshot);
    }

    // outer_end 行。GPU 経路は last_slot_ を参照する行を buffer し、batch が満ちたら flush。
    void logOuterEnd(
        int step,
        const ResidualSnapshot* snapshot = nullptr,
        const CorrectionSnapshot* correction_snapshot = nullptr)
    {
        if (gpu_) {
            rows_.push_back({step, -1, ResidualPhase::OuterEnd, last_slot_});
            if (device_count_ >= capacity_) {
                flushDevice();
            }
            return;
        }
        if (snapshot != nullptr) {
            last_snapshot_ = *snapshot;
        }
        if (correction_snapshot != nullptr) {
            last_correction_snapshot_ = *correction_snapshot;
        }
        writeRowImmediate(step, -1, ResidualPhase::OuterEnd, last_snapshot_, last_correction_snapshot_);
    }

    // buffer に残った残差を強制 flush する (停止前・console モニタ行の直前などに使用)。
    void flush()
    {
        if (gpu_) {
            flushDevice();
        }
    }

    // console モニタ行用の要約: 最新 outer_end 行の rms と、step 0 outer_begin の rms (基準)。
    // flush 済みの値 (= CSV に書いた値そのもの) を返すので追加 reduction は無い。
    struct Summary {
        bool valid = false;
        int step = -1;
        std::vector<std::string> names;
        std::vector<double> rms;    // 最新 outer_end
        std::vector<double> rms0;   // step 0 outer_begin (基準)。未取得なら空
    };
    Summary latestSummary() const
    {
        Summary sm;
        sm.valid = last_end_valid_;
        sm.step = last_end_step_;
        sm.names = residual_names_;
        sm.rms = last_end_rms_;
        if (has_rms0_) sm.rms0 = rms0_;
        return sm;
    }

private:
    // 要約状態の更新 (CPU/GPU 両経路の行書き出しから呼ぶ)。
    void noteRow(int step, ResidualPhase phase, const flow_float* v, int n)
    {
        if (phase == ResidualPhase::OuterBegin && !has_rms0_) {
            rms0_.assign(v, v + n);
            has_rms0_ = true;
        }
        if (phase == ResidualPhase::OuterEnd) {
            last_end_rms_.assign(v, v + n);
            last_end_step_ = step;
            last_end_valid_ = true;
        }
    }

    struct RowDesc {
        int step;
        int inner;
        ResidualPhase phase;
        int slot;
    };
    static constexpr int kIntraStepMargin = 1024;

    static const char* phaseName(ResidualPhase phase)
    {
        switch (phase) {
            case ResidualPhase::OuterBegin: return "outer_begin";
            case ResidualPhase::InnerBegin: return "inner_begin";
            case ResidualPhase::InnerIter: return "inner_iter";
            case ResidualPhase::OuterEnd: return "outer_end";
        }

        return "unknown";
    }

    void writeRowImmediate(
        int step,
        int inner,
        ResidualPhase phase,
        const ResidualSnapshot& snapshot,
        const CorrectionSnapshot& correction_snapshot)
    {
        stream_ << step << ',' << inner << ',' << phaseName(phase);
        for (flow_float value : snapshot.rms) {
            stream_ << ',' << std::setprecision(16) << value;
        }
        for (flow_float value : correction_snapshot.rms) {
            stream_ << ',' << std::setprecision(16) << value;
        }
        stream_ << "\n";
        stream_.flush();
        noteRow(step, phase, snapshot.rms.data(), static_cast<int>(snapshot.rms.size()));
    }

    // device buffer を host へ一括転送し、buffer 済みの全行を書き出して reset する (同期点はここだけ)。
    void flushDevice()
    {
        if (device_count_ <= 0) {
            return;
        }
        downloadRmsBuffer(d_rms_buf_, host_rms_.data(), device_count_, reducer_.nVar);
        for (const RowDesc& rd : rows_) {
            stream_ << rd.step << ',' << rd.inner << ',' << phaseName(rd.phase);
            const flow_float* v = &host_rms_[static_cast<size_t>(rd.slot) * reducer_.nVar];
            for (int i = 0; i < reducer_.nVar; ++i) {
                stream_ << ',' << std::setprecision(16) << v[i];
            }
            for (int i = 0; i < nDq_; ++i) {
                stream_ << ',' << std::setprecision(16) << static_cast<flow_float>(0.0);  // rms_dq_* は常に 0
            }
            stream_ << "\n";
            noteRow(rd.step, rd.phase, v, reducer_.nVar);
        }
        stream_.flush();
        rows_.clear();
        device_count_ = 0;
    }

    std::ofstream stream_;
    std::vector<std::string> residual_names_;
    bool gpu_ = false;
    int nVar_ = 0;
    int nDq_ = 0;
    ResidualSnapshot last_snapshot_{};
    CorrectionSnapshot last_correction_snapshot_{};

    // GPU buffering 状態
    DeviceResidualReducer reducer_{};
    flow_float* d_rms_buf_ = nullptr;
    std::vector<flow_float> host_rms_;
    std::vector<RowDesc> rows_;
    int capacity_ = 1;
    int buf_capacity_ = 1;
    int device_count_ = 0;
    int last_slot_ = 0;

    // console モニタ行用の要約状態
    std::vector<double> rms0_;
    bool has_rms0_ = false;
    std::vector<double> last_end_rms_;
    int last_end_step_ = -1;
    bool last_end_valid_ = false;
};

// console モニタ行 (methods/architecture/overview.md §8.5, plans/accepted/architecture-runtime-monitor-line.md)。
// monitorInterval ステップごとに 1 行: [unsteady のみ t/dt/maxCFL] | 壁時計 ms/step・経過・ETA | 残差要約。
// 定常 (unsteady==0) では cfg.dt/max cfl は dt_local から打ち消されて無意味なので出さない (陽・陰とも)。
class StepMonitor {
public:
    using clock = std::chrono::steady_clock;

    StepMonitor(const solverConfig& cfg, ResidualCsvLogger& logger)
        : cfg_(cfg), logger_(logger), t_start_(clock::now()), t_last_(t_start_) {}

    // 起動時のモード要約 (1 回)。
    void printHeader() const
    {
        const bool unsteady = (cfg_.unsteady == 1);
        const bool implicit = (cfg_.isImplicit == 1);
        std::string mode;
        if (!unsteady && implicit)      mode = "steady implicit (block-DPLUR, local pseudo-dt)";
        else if (!unsteady && !implicit) mode = "steady explicit (local pseudo-dt)";
        else if (unsteady && implicit)  mode = "unsteady dual-time implicit";
        else                            mode = "unsteady explicit RK";
        std::cout << "[monitor] mode: " << mode
                  << "  nStepOuter=" << cfg_.mainLoopCount()
                  << "  monitorInterval=" << cfg_.monitorInterval;
        if (implicit) {
            std::cout << "  cfl_pseudo=" << cfg_.cfl_pseudo
                      << "  implicitRelax=" << cfg_.implicitRelax
                      << "  nStepInner=" << cfg_.nStepInner;
            if (unsteady) std::cout << "  nSubIterDualTime=" << cfg_.nSubIterDualTime;
        } else if (!unsteady) {
            std::cout << "  cfl_pseudo=" << cfg_.cfl_pseudo;
        }
        if (unsteady) {
            std::cout << "  dt=" << cfg_.dt
                      << (cfg_.dtControl == 1 ? " (CFL-adaptive)" : " (fixed)");
        }
        std::cout << "\n";
        std::cout << "[monitor] columns: step"
                  << (unsteady ? ((cfg_.dtControl == 1 && !implicit) ? " | t dt maxCFL dt_next" : " | t dt maxCFL") : "")
                  << " | ms/step elapsed eta | rms_ro (log10 change vs step 0) worst column (log10 change) [from0 column]\n";
    }

    // 各ステップ末尾で呼ぶ。monitor step のみ出力 (それ以外は何もしない = 同期なし)。
    void report(int iStep)
    {
        steps_since_++;
        if (iStep % cfg_.monitorInterval != 0) return;

        logger_.flush();   // 残差を host へ (monitor step ではいずれ起きる D2H を前倒し)
        const auto now = clock::now();
        const double elapsed = std::chrono::duration<double>(now - t_start_).count();
        const double interval = std::chrono::duration<double>(now - t_last_).count();
        const double ms_per_step = (steps_since_ > 0) ? interval * 1.0e3 / steps_since_ : 0.0;
        const int remaining = cfg_.mainLoopCount() - (iStep + 1);
        const double eta = (remaining > 0) ? remaining * ms_per_step * 1.0e-3 : 0.0;
        t_last_ = now;
        steps_since_ = 0;

        std::ostringstream os;
        os << "step " << std::setw(8) << (iStep + 1);
        if (cfg_.unsteady == 1) {
            // dt と maxCFL は同じ刻みで対にする: monitorCflDt = この step の前進に使った dt で評価した CFL。
            // dtControl==1 (適応) では setDT が既に cfg.dt を次 step 用に更新しているので dt_next として別表示。
            const bool haveCfl = (cfg_.monitorCflMax >= 0.0 && cfg_.monitorCflDt > 0.0);
            const double dtUsed = haveCfl ? cfg_.monitorCflDt : cfg_.dt;
            os << " | t " << std::scientific << std::setprecision(4) << cfg_.totalTime
               << " dt " << std::setprecision(2) << dtUsed;
            if (haveCfl) {
                os << " maxCFL " << std::fixed << std::setprecision(2) << cfg_.monitorCflMax;
            }
            if (cfg_.dtControl == 1 && cfg_.isImplicit == 0) {
                os << " dt_next " << std::scientific << std::setprecision(2) << cfg_.dt;
            }
        }
        os << " | " << std::fixed << std::setprecision(2) << ms_per_step << " ms/step"
           << " elapsed " << std::setprecision(1) << elapsed << " s"
           << " eta " << std::setprecision(0) << eta << " s";

        const auto sm = logger_.latestSummary();
        if (sm.valid && !sm.rms.empty()) {
            bool nonfinite = false;
            int worst = -1;
            double worst_dec = std::numeric_limits<double>::infinity();
            int from0 = -1;          // 初期残差 0 → 現在非 0 の列 (低下桁数が定義できないので別表示)
            double from0_val = 0.0;
            std::vector<double> dec(sm.rms.size(), std::numeric_limits<double>::quiet_NaN());
            for (size_t i = 0; i < sm.rms.size(); ++i) {
                if (!std::isfinite(sm.rms[i])) { nonfinite = true; continue; }
                if (i < sm.rms0.size() && sm.rms0[i] > 0.0 && sm.rms[i] > 0.0) {
                    dec[i] = std::log10(sm.rms0[i] / sm.rms[i]);
                    if (dec[i] < worst_dec) { worst_dec = dec[i]; worst = static_cast<int>(i); }
                } else if (i < sm.rms0.size() && sm.rms0[i] == 0.0 && sm.rms[i] > from0_val) {
                    from0 = static_cast<int>(i); from0_val = sm.rms[i];
                }
            }
            os << " | rms_" << sm.names[0] << " " << std::scientific << std::setprecision(2) << sm.rms[0];
            if (std::isfinite(dec[0])) os << " (" << std::fixed << std::showpos << std::setprecision(1) << -dec[0] << std::noshowpos << ")";
            if (worst >= 0 && worst != 0) {
                os << " worst rms_" << sm.names[worst] << " " << std::scientific << std::setprecision(2) << sm.rms[worst]
                   << " (" << std::fixed << std::showpos << std::setprecision(1) << -dec[worst] << std::noshowpos << ")";
            }
            if (from0 >= 0) {
                os << " from0 rms_" << sm.names[from0] << " " << std::scientific << std::setprecision(2) << sm.rms[from0];
            }
            if (nonfinite) os << " *** NaN ***";
        }
        std::cout << os.str() << "\n";
        std::cout.flush();
    }

    double elapsedSeconds() const
    {
        return std::chrono::duration<double>(clock::now() - t_start_).count();
    }

private:
    const solverConfig& cfg_;
    ResidualCsvLogger& logger_;
    clock::time_point t_start_;
    clock::time_point t_last_;
    int steps_since_ = 0;
};

class ImplicitDiagLogger {
public:
    ImplicitDiagLogger()
    {
        const char* path = std::getenv("FORGE_IMPLICIT_DIAG_CSV");
        if (path == nullptr || path[0] == '\0') {
            return;
        }

        stream_.open(path, std::ios::out | std::ios::trunc);
        if (!stream_) {
            throw std::runtime_error("Failed to open implicit diag log file: " + std::string(path));
        }

        enabled_ = true;
        stream_ << "step,inner,mean_pseudo_diag,mean_face_diag,mean_pseudo_ratio,min_pseudo_ratio,max_pseudo_ratio,mean_face_to_pseudo,mean_dt_local\n";
    }

    bool enabled() const
    {
        return enabled_;
    }

    void log(int step, int inner, const ImplicitDiagSnapshot& snapshot)
    {
        if (!enabled_) {
            return;
        }

        stream_ << step << ',' << inner
                << ',' << std::setprecision(16) << snapshot.mean_pseudo_diag
                << ',' << snapshot.mean_face_diag
                << ',' << snapshot.mean_pseudo_ratio
                << ',' << snapshot.min_pseudo_ratio
                << ',' << snapshot.max_pseudo_ratio
                << ',' << snapshot.mean_face_to_pseudo
                << ',' << snapshot.mean_dt_local
                << "\n";
        stream_.flush();
    }

private:
    bool enabled_ = false;
    std::ofstream stream_;
};

ResidualSnapshot gatherResidualSnapshot(solverConfig& cfg, mesh& msh, variables& var)
{
    std::vector<std::string> variable_names = {
        "res_ro", "res_roUx", "res_roUy", "res_roUz", "res_roe"
    };
    if (scalarResidualEnabled(cfg)) {
        variable_names.emplace_back("res_roK");
        variable_names.emplace_back("res_roOmega");
    }
    for (const auto& name : speciesResidualNames(cfg)) {
        variable_names.emplace_back("res_" + name);
    }
    if (cfg.tracerEnabled()) {
        variable_names.emplace_back("res_roXi");
    }
    if (condensationEnabled(cfg)) {
        for (const auto& name : condMomentConsNames(cfg)) {
            variable_names.emplace_back("res_" + name);
        }
    }

    ResidualSnapshot snapshot;
    if (cfg.gpu == 1) {
        snapshot.rms = gatherVariableRms_d(msh, var, variable_names);
        return snapshot;
    }

    snapshot.rms.assign(variable_names.size(), static_cast<flow_float>(0.0));

    for (std::size_t i = 0; i < variable_names.size(); ++i) {
        const std::vector<flow_float>& values = var.c.at(variable_names[i]);
        flow_float sum_sq = 0.0;
        for (geom_int ic = 0; ic < msh.nCells; ++ic) {
            const flow_float value = values[ic];
            sum_sq += value * value;
        }
        snapshot.rms[i] = std::sqrt(sum_sq / static_cast<flow_float>(std::max<geom_int>(msh.nCells, 1)));
    }

    return snapshot;
}

// frozen scalar 陰解法（次フェーズ）の inner 補正ログ用。現状の経路では未使用だが温存する。
[[maybe_unused]] CorrectionSnapshot gatherCorrectionSnapshot(solverConfig& cfg, mesh& msh, variables& var)
{
    CorrectionSnapshot snapshot;

    if (cfg.timeIntegration != 11) {
        return snapshot;
    }

    const std::array<const char*, kResidualEquationNames.size()> variable_names = {
        cfg.blockDPLUR == 1 ? "dq_block_old_0" : "dq_ro_old",
        cfg.blockDPLUR == 1 ? "dq_block_old_1" : "dq_roUx_old",
        cfg.blockDPLUR == 1 ? "dq_block_old_2" : "dq_roUy_old",
        cfg.blockDPLUR == 1 ? "dq_block_old_3" : "dq_roUz_old",
        cfg.blockDPLUR == 1 ? "dq_block_old_4" : "dq_roe_old"
    };

    if (cfg.gpu == 1) {
        snapshot.rms = gatherVariableRms_d(msh, var, variable_names);
        return snapshot;
    }

    for (std::size_t i = 0; i < variable_names.size(); ++i) {
        const std::vector<flow_float>& values = var.c.at(variable_names[i]);
        flow_float sum_sq = 0.0;
        for (geom_int ic = 0; ic < msh.nCells; ++ic) {
            const flow_float value = values[ic];
            sum_sq += value * value;
        }
        snapshot.rms[i] = std::sqrt(sum_sq / static_cast<flow_float>(std::max<geom_int>(msh.nCells, 1)));
    }

    return snapshot;
}

ImplicitDiagSnapshot gatherImplicitDiagSnapshot(solverConfig& cfg, mesh& msh, variables& var)
{
    // 名前は確保側 (hostCellSet) と共用 (output/outputFieldNames.cpp)。ホスト側は hostCell で長さを検査してから読む。
    const std::list<std::string> variable_names = implicitDiagCellNames();

    if (cfg.gpu == 1) {
        var.copyVariables_cell_D2H(variable_names);
    }

    const std::vector<flow_float>& ro = var.hostCell("ro");
    const std::vector<flow_float>& ux = var.hostCell("Ux");
    const std::vector<flow_float>& uy = var.hostCell("Uy");
    const std::vector<flow_float>& uz = var.hostCell("Uz");
    const std::vector<flow_float>& sonic = var.hostCell("sonic");
    const std::vector<flow_float>& vis_turb = var.hostCell("vis_turb");
    const std::vector<flow_float>& dt_local = var.hostCell("dt_local");

    ImplicitDiagSnapshot snapshot;
    snapshot.min_pseudo_ratio = std::numeric_limits<double>::infinity();

    for (geom_int ic = 0; ic < msh.nCells; ++ic) {
        const flow_float volume = static_cast<flow_float>(msh.cells[ic].volume);
        const flow_float dt_l = std::max(dt_local[ic], static_cast<flow_float>(1.0e-30));
        const flow_float density = std::max(ro[ic], static_cast<flow_float>(1.0e-30));
        const flow_float nu_eff = (cfg.visc + std::max(vis_turb[ic], static_cast<flow_float>(0.0))) / density;
        const flow_float local_sonic = std::max(sonic[ic], static_cast<flow_float>(0.0));

        double face_diag_sum = 0.0;
        for (const geom_int ip : msh.cells[ic].iPlanes) {
            const plane& pln = msh.planes[ip];
            const flow_float sx = static_cast<flow_float>(pln.surfVect[0]);
            const flow_float sy = static_cast<flow_float>(pln.surfVect[1]);
            const flow_float sz = static_cast<flow_float>(pln.surfVect[2]);
            const flow_float face_area = std::max(static_cast<flow_float>(pln.surfArea), static_cast<flow_float>(1.0e-30));

            const flow_float advective_radius = std::fabs(
                ux[ic] * sx + uy[ic] * sy + uz[ic] * sz
            ) / face_area + local_sonic;

            const geom_int ic0 = pln.iCells[0];
            const geom_int ic1 = pln.iCells[1];
            const geom_int other_ic = (ic0 == ic) ? ic1 : ic0;

            const geom_float dcc_x = msh.cells[other_ic].centCoords[0] - msh.cells[ic].centCoords[0];
            const geom_float dcc_y = msh.cells[other_ic].centCoords[1] - msh.cells[ic].centCoords[1];
            const geom_float dcc_z = msh.cells[other_ic].centCoords[2] - msh.cells[ic].centCoords[2];
            const flow_float dcc = std::max(
                static_cast<flow_float>(std::sqrt(dcc_x * dcc_x + dcc_y * dcc_y + dcc_z * dcc_z)),
                static_cast<flow_float>(1.0e-30)
            );
            const flow_float dcc_dot_s = std::max(
                static_cast<flow_float>(std::fabs(dcc_x * sx + dcc_y * sy + dcc_z * sz)),
                static_cast<flow_float>(1.0e-30)
            );
            const flow_float delta = std::max(
                dcc * face_area * face_area / dcc_dot_s,
                static_cast<flow_float>(1.0e-30)
            );
            const flow_float viscous_radius = static_cast<flow_float>(2.0) * nu_eff / delta;
            const flow_float face_coeff = face_area * (advective_radius + viscous_radius);
            face_diag_sum += static_cast<double>(face_coeff);
        }

        const double pseudo_diag = static_cast<double>(volume / dt_l);
        const double total_diag = std::max(pseudo_diag + face_diag_sum, 1.0e-30);
        const double pseudo_ratio = pseudo_diag / total_diag;
        const double face_to_pseudo = face_diag_sum / std::max(pseudo_diag, 1.0e-30);

        snapshot.mean_pseudo_diag += pseudo_diag;
        snapshot.mean_face_diag += face_diag_sum;
        snapshot.mean_pseudo_ratio += pseudo_ratio;
        snapshot.min_pseudo_ratio = std::min(snapshot.min_pseudo_ratio, pseudo_ratio);
        snapshot.max_pseudo_ratio = std::max(snapshot.max_pseudo_ratio, pseudo_ratio);
        snapshot.mean_face_to_pseudo += face_to_pseudo;
        snapshot.mean_dt_local += static_cast<double>(dt_l);
    }

    const double inv_n_cells = 1.0 / static_cast<double>(std::max<geom_int>(msh.nCells, 1));
    snapshot.mean_pseudo_diag *= inv_n_cells;
    snapshot.mean_face_diag *= inv_n_cells;
    snapshot.mean_pseudo_ratio *= inv_n_cells;
    snapshot.mean_face_to_pseudo *= inv_n_cells;
    snapshot.mean_dt_local *= inv_n_cells;
    if (!std::isfinite(snapshot.min_pseudo_ratio)) {
        snapshot.min_pseudo_ratio = 0.0;
    }

    return snapshot;
}

struct ProfileStat {
    double total_ms = 0.0;
    std::size_t calls = 0;
};

class RuntimeProfiler {
public:
    RuntimeProfiler()
        : enabled_(readFlag("FORGE_PROFILE")),
          verbose_(readFlag("FORGE_PROFILE_VERBOSE"))
    {
        if (enabled_) {
            CHECK_CUDA_ERROR(cudaEventCreate(&start_event_));
            CHECK_CUDA_ERROR(cudaEventCreate(&stop_event_));
        }
    }

    ~RuntimeProfiler()
    {
        if (!enabled_) {
            return;
        }

        cudaEventDestroy(start_event_);
        cudaEventDestroy(stop_event_);
    }

    bool enabled() const
    {
        return enabled_;
    }

    template <typename Func>
    void measureWall(ProfileSection section, Func&& func)
    {
        if (!enabled_) {
            func();
            return;
        }

        const auto start = std::chrono::steady_clock::now();
        func();
        const auto end = std::chrono::steady_clock::now();
        const double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
        add(section, elapsed_ms);
    }

    template <typename Func>
    void measureCuda(ProfileSection section, Func&& func)
    {
        if (!enabled_) {
            func();
            return;
        }

        CHECK_CUDA_ERROR(cudaEventRecord(start_event_));
        func();
        CHECK_CUDA_ERROR(cudaEventRecord(stop_event_));
        CHECK_CUDA_ERROR(cudaEventSynchronize(stop_event_));

        float elapsed_ms = 0.0f;
        CHECK_CUDA_ERROR(cudaEventElapsedTime(&elapsed_ms, start_event_, stop_event_));
        add(section, static_cast<double>(elapsed_ms));
    }

    void printSummary() const
    {
        if (!enabled_) {
            return;
        }

        const ProfileStat& total = stats_[index(ProfileSection::StepTotal)];
        std::cout << "\n=== Runtime Profile Summary ===\n";
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "Profiled steps: " << total.calls << "\n";
        if (total.calls == 0 || total.total_ms <= 0.0) {
            std::cout << "No profiled steps were recorded.\n";
            return;
        }

        for (std::size_t i = 0; i < kProfileSectionCount; ++i) {
            const auto section = static_cast<ProfileSection>(i);
            if (section == ProfileSection::StepTotal) {
                continue;
            }

            const ProfileStat& stat = stats_[i];
            if (stat.calls == 0) {
                continue;
            }

            const double avg_ms = stat.total_ms / static_cast<double>(stat.calls);
            const double share_pct = (stat.total_ms / total.total_ms) * 100.0;
            std::cout << std::setw(20) << std::left << sectionName(section)
                      << " total_ms=" << std::setw(12) << stat.total_ms
                      << " avg_ms=" << std::setw(10) << avg_ms
                      << " share=" << std::setw(8) << share_pct << "%"
                      << " calls=" << stat.calls << "\n";
        }

        if (verbose_) {
            std::cout << "Profiling source: CPU sections use steady_clock, GPU wrapper sections use CUDA events on the default stream.\n";
        }
    }

private:
    static bool readFlag(const char* name)
    {
        const char* value = std::getenv(name);
        if (value == nullptr) {
            return false;
        }

        return std::string(value) != "0";
    }

    static constexpr std::size_t index(ProfileSection section)
    {
        return static_cast<std::size_t>(section);
    }

    void add(ProfileSection section, double elapsed_ms)
    {
        ProfileStat& stat = stats_[index(section)];
        stat.total_ms += elapsed_ms;
        stat.calls += 1;
    }

    static const char* sectionName(ProfileSection section)
    {
        switch (section) {
            case ProfileSection::UpdateInner: return "update_inner";
            case ProfileSection::DependentVariables: return "dependent_vars";
            case ProfileSection::GasProperties: return "gas_properties";
            case ProfileSection::ApplyBconds: return "apply_bconds";
            case ProfileSection::CalcGradient: return "calc_gradient";
            case ProfileSection::Limiter: return "limiter";
            case ProfileSection::DucrosSensor: return "ducros_sensor";
            case ProfileSection::TurbulenceModel: return "turbulence_model";
            case ProfileSection::ConvectiveFlux: return "convective_flux";
            case ProfileSection::AxisymmetricSource: return "axisym_source";
            case ProfileSection::ViscousFlux: return "viscous_flux";
            case ProfileSection::ImplicitCorrection: return "implicit_corr";
            case ProfileSection::TimeIntegration: return "time_integr";
            case ProfileSection::UpdateOuter: return "update_outer";
            case ProfileSection::WriteOutputs: return "write_outputs";
            case ProfileSection::SetDt: return "set_dt";
            case ProfileSection::StepTotal: return "step_total";
            case ProfileSection::Count: return "count";
        }

        return "unknown";
    }

    bool enabled_ = false;
    bool verbose_ = false;
    cudaEvent_t start_event_{};
    cudaEvent_t stop_event_{};
    std::array<ProfileStat, kProfileSectionCount> stats_{};
};


// dual-time の履歴契約 (plan species-passive-scalar-unification §4.4, codex plan-2 M6): valueFileName の /CHECKPOINT が
// 流れ・化学種・受動種の前物理レベル (流れ *N, 化学種 roY{s}P, 受動種 <cons>P) と totalTime/dt/nHistoryValid/layout を
// 全部持ち layout が一致するときだけ復元して BDF2 を継続する。1 つでも欠ければ (旧形式 res_*.h5 を含む) 全系を
// P = PP = 現在値に揃え nHistoryValid=0 (最初の物理 step は BDF1) で再開する。updateVariablesOuter (roN=ro) の後に呼ぶ。
static void initDualTimeHistory(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!(cfg.unsteady == 1 && cfg.dualTime == 1)) return;
    // まず全系を P = PP = 現在値に (復元に失敗しても一貫した状態)。
    speciesInitDualTimeLevels_d_wrapper(cfg, cuda_cfg, msh, var);
    passiveInitDualTimeLevels_d_wrapper(cfg, cuda_cfg, msh, var);
    cfg.nHistoryValid = 0;

    // 履歴の名前は書出し側 (output.cpp の /CHECKPOINT) と確保側 (hostCellSet) と共用 (output/outputFieldNames.cpp)。
    // scheme 0 の受動種は物理時間項を持たないので履歴は要らない (書きもしない)。
    const std::list<std::string> hist = dualTimeHistoryNames(cfg, var);
    // 復元条件 (codex result M3): 配列構成に加え、履歴を生成した物理 dt (相対 1e-12)・bdfOrder・passiveScalarScheme・
    // speciesImplicitCoupling が一致すること。1 つでも違えば全系を BDF1 から再開する (刻み変更 restart は最初の step の時間微分が狂う)。
    char dtbuf[64]; std::snprintf(dtbuf, sizeof(dtbuf), "%.17g", (double)cfg.dt);
    const std::string layout = "nSpecies=" + std::to_string(var.nSpeciesRegistered) + ";tracer=" + std::to_string(var.tracerRegistered)
                             + ";nCond=" + std::to_string(var.nCondSpeciesRegistered)
                             + ";dt=" + std::string(dtbuf) + ";bdfOrder=" + std::to_string(cfg.bdfOrder)
                             + ";passiveScalarScheme=" + std::to_string(cfg.passiveScalarScheme)
                             + ";speciesImplicitCoupling=" + std::to_string(cfg.speciesImplicitCoupling)
                             + ";passiveFct=" + std::to_string(passiveFctConfigured(cfg) ? 1 : 0);
    std::string why;
    try {
        HighFive::File file(cfg.valueFileName, HighFive::File::ReadOnly);
        if (!file.exist("/CHECKPOINT")) {
            why = "no /CHECKPOINT group (old-format input)";
        } else {
            HighFive::Group ck = file.getGroup("/CHECKPOINT");
            std::string lay; double tt = 0.0, dtp = 0.0; int nh = 0;
            if (!ck.hasAttribute("layout") || !ck.hasAttribute("totalTime") || !ck.hasAttribute("dt") || !ck.hasAttribute("nHistoryValid")) {
                why = "missing checkpoint attributes";
            } else {
                ck.getAttribute("layout").read(lay); ck.getAttribute("totalTime").read(tt);
                ck.getAttribute("dt").read(dtp);     ck.getAttribute("nHistoryValid").read(nh);
                if (std::abs(dtp - (double)cfg.dt) > 1.0e-12 * std::max(std::abs(dtp), std::abs((double)cfg.dt)))
                    why = "physical dt differs from the checkpoint (file " + std::to_string(dtp) + " vs run " + std::to_string((double)cfg.dt) + ")";
                else if (lay != layout) why = (lay.find("passiveFct=0") != std::string::npos && layout.find("passiveFct=1") != std::string::npos)
                                            ? "passive FCT flux-form history missing (checkpoint written without FCT: layout '" + lay + "')"
                                            : "layout mismatch (file '" + lay + "' vs run '" + layout + "')";
                else if (nh < 1) why = "checkpoint has no valid history (nHistoryValid=0)";
                else {
                    for (const auto& nm : hist) if (!file.exist("/CHECKPOINT/" + nm)) { why = "missing dataset /CHECKPOINT/" + nm; break; }
                }
            }
            if (why.empty()) {
                std::list<std::string> names;
                for (const auto& nm : hist) {
                    std::vector<geom_float> in; file.getDataSet("/CHECKPOINT/" + nm).read(in);
                    if ((geom_int)in.size() < msh.nCells) { why = "dataset /CHECKPOINT/" + nm + " too short"; break; }
                    std::vector<flow_float>& v = var.hostCell(nm);   // 書く前に長さを検査する (H に無ければ名前つきで停止)
                    for (geom_int i = 0; i < msh.nCells; ++i) v[i] = static_cast<flow_float>(in[i]);
                    names.push_back(nm);
                }
                if (why.empty()) {
                    var.copyVariables_cell_H2D(names);
                    // 受動種 P は host 名で H2D 済み; 化学種 roY{s}P も同様。PP はシフトで上書きされるので不要。
                    cfg.totalTime = static_cast<flow_float>(tt); cfg.totalTimeD = tt;
                    cfg.nHistoryValid = std::min(nh, 2);
                    // 受動種 FCT の流束形履歴 G/H/ṁ^eff (§4.7 v4, plan-6 M5): FCT 有効なら**必須** (無ければ全系を BDF1 に揃えて再開)。
                    if (passiveFctConfigured(cfg)) {
                        std::vector<std::vector<flow_float>> G, H, Hs; std::vector<flow_float> mEff; std::string miss;
                        const auto& cons = passive_cons_names();
                        for (const auto& cn : cons) {
                            if (!file.exist("/CHECKPOINT/" + cn + "_fctG")) { miss = cn + "_fctG"; break; }
                            if (!file.exist("/CHECKPOINT/" + cn + "_fctH")) { miss = cn + "_fctH"; break; }
                            if (!file.exist("/CHECKPOINT/" + cn + "_fctHsrc")) { miss = cn + "_fctHsrc"; break; }
                            std::vector<geom_float> g, h, hs; file.getDataSet("/CHECKPOINT/" + cn + "_fctG").read(g); file.getDataSet("/CHECKPOINT/" + cn + "_fctH").read(h); file.getDataSet("/CHECKPOINT/" + cn + "_fctHsrc").read(hs);
                            G.emplace_back(g.begin(), g.end()); H.emplace_back(h.begin(), h.end()); Hs.emplace_back(hs.begin(), hs.end());
                        }
                        if (miss.empty() && !file.exist("/CHECKPOINT/passive_fctMeff")) miss = "passive_fctMeff";
                        if (miss.empty()) { std::vector<geom_float> m; file.getDataSet("/CHECKPOINT/passive_fctMeff").read(m); mEff.assign(m.begin(), m.end()); }
                        if (!miss.empty() || !passiveFctHistoryFromHost(cfg, msh, var, G, H, mEff, Hs)) {
                            cfg.nHistoryValid = 0;
                            std::cout << "[dual-time] passive FCT flux-form history missing or inconsistent (" << (miss.empty() ? std::string("size mismatch") : miss)
                                      << "): all systems start from P=PP=current, first physical step is BDF1\n";
                            passiveInitDualTimeLevels_d_wrapper(cfg, cuda_cfg, msh, var);
                            return;
                        }
                        std::cout << "[dual-time] passive FCT flux-form history (G/H/mEff) restored\n";
                    }
                    std::cout << "[dual-time] history restored from " << cfg.valueFileName << " (/CHECKPOINT: " << names.size()
                              << " levels, totalTime=" << tt << ", dt_file=" << dtp << ", nHistoryValid=" << cfg.nHistoryValid
                              << ")\n";
                    return;
                }
            }
        }
    } catch (const std::exception& e) {
        why = std::string("read error: ") + e.what();
    }
    std::cout << "[dual-time] history NOT restored (" << why << "): all systems start from P=PP=current, first physical step is BDF1\n";
}

// 化学種 datum の有効温度 (thermo_init_db と同じ条件: thermoHrefTemp>0 のときだけ datum を適用)
static double speciesRecordTref(const solverConfig& cfg)
{
    return (cfg.thermoHrefTemp > 0.0) ? cfg.thermoHrefTemp : 0.0;
}

// valueFileName の species_hash 属性を照合し、run ディレクトリ (cwd) に解決済み記録を書く。
//   属性あり: 互換性ハッシュが一致 → 通す / 不一致 → 差を示して終了 (許可手段なし)。
//   属性なし: 照合不能として終了。その実行だけ FORGE_ALLOW_UNVERIFIED_SPECIES=1 で通し、出力に未検証の印を付ける。
//   CPG (thermalMethod != 2) は対象外。
static void checkInputSpeciesAndWriteRecord(const solverConfig& cfg)
{
    if (cfg.thermalMethod != 2) return;
    const ResolvedSpeciesDB* db = speciesDB_current();
    if (db == nullptr) db = &speciesDB_init(cfg);
    const double Tref = speciesRecordTref(cfg);

    std::string fieldHash, fieldRecordSha;
    int fieldUnverified = -1;
    try {
        HighFive::File f(cfg.valueFileName, HighFive::File::ReadOnly);
        if (f.hasAttribute("species_hash"))          f.getAttribute("species_hash").read(fieldHash);
        if (f.hasAttribute("species_record_sha256")) f.getAttribute("species_record_sha256").read(fieldRecordSha);
        if (f.hasAttribute("species_input_unverified")) f.getAttribute("species_input_unverified").read(fieldUnverified);
    } catch (const std::exception& e) {
        std::cerr << "[species] cannot read attributes of " << cfg.valueFileName << ": " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    // 属性なしの場は既定で停止する (plan thermophysics-solver-owned-species-db #3c)。許可はその実行だけの
    // FORGE_ALLOW_UNVERIFIED_SPECIES=1 のみ (出力に未検証の印)。係数不一致は常に停止 (allow は不一致を通さない)。
    const char* env = std::getenv("FORGE_ALLOW_UNVERIFIED_SPECIES");
    const bool allow = (env != nullptr && std::string(env) == "1");
    std::vector<std::string> dirs;
    {
        std::filesystem::path d = std::filesystem::path(cfg.valueFileName).parent_path();
        dirs.push_back(d.empty() ? std::string(".") : d.string());
        dirs.push_back(".");
    }
    std::string status, msg;
    int unverified = 0;
    const bool ok = speciesDB_checkInputField(*db, Tref, fieldHash, fieldRecordSha, fieldUnverified, cfg.valueFileName,
                                              dirs, allow, status, unverified, msg);
    if (!ok) {
        std::cerr << "[species] ERROR: " << msg << std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::cout << "[species] " << (unverified ? "WARNING: " : "") << msg << "\n";
    try {
        const SpeciesRecordInfo rec = speciesDB_writeRecord(*db, Tref, cfg.speciesDBFile, cfg.valueFileName, status, unverified, ".");
        speciesDB_setCurrentRecord(rec);
        std::cout << "[species] species_hash " << rec.compatHash << " (record " << rec.recordFile
                  << ", sha256 " << rec.recordSha256.substr(0, 16) << ", input " << status << ")" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

// 試験用 (FORGE_TRANSPORT_PROBE=<states.txt>; plan thermophysics-solver-owned-species-db #5t2-2、tests/unit/test_transport_gpu.py)。
//   物性だけを評価する CFD 0 step のハーネス: 初期化後、時間更新をせずに
//     (1) 状態表の (T, ρ, Y) を全セル (ghost 込み nCells_all) と全境界面の解点へ割り当てて float で device に書き、
//     (2) 実際のセル経路 gasProperties_d_wrapper (vis_lam・thermCond を float で格納) を 1 回、
//     (3) 同じ組成・評価関数の double 値 (gasPropertiesTransportProbe_d_wrapper) を 1 回、
//     (4) 壁モデルの物性関数 (wmlesTransportProbe; 全 bcond の境界面) を 1 回、
//     (5) 状態表の double の Y を float の roY を経由せずに同じ評価関数へ渡す (gasPropertiesTransportProbeStates_d_wrapper)、
//     (6) 表引き (#5t2-3) があれば、同じ Y・T を float に丸めて表引きの経路へ渡す (gasPropertiesTransportProbeStatesTab_d_wrapper)
//   を行い、device から読み戻した実際の入力 (T・ρ・roY) と出力を <states.txt>.bin に、配置を <states.txt>.json に書いて終了する。
//   状態の割り当ては pass ごとにずらし、どの状態もセル・ghost・境界面の解点に少なくとも 1 回は乗るようにする。
//   states.txt: 1 行目 "K nSpecies"、以降 K 行 "T ro Y0 .. Y{n-1}" (Y は輸送種の質量分率)。
//   凝縮 carrier の気相組成の試験 (plan condensation-two-phase-transport §4.1, #3): 1 行目を "K nSpecies 1" にすると各行の末尾に
//   液の質量分率 g を読み、rog_0 = ρg を全セルに書く (TP carrier の凝縮 run でなければエラー)。このとき pass ごとに
//   読み戻した rog・P と、化学種拡散と同じ組成・関数の分子拡散係数 D_s (speciesDmixProbe_d_wrapper) も書き、
//   (5)(6) は同じ g で気相組成にしてから評価する。
static int runTransportProbe(const char* statesPath, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    std::ifstream in(statesPath);
    int K = 0, nS = 0, nLiq = 0;
    {
        std::string head;
        std::getline(in, head);
        std::istringstream hs(head);
        hs >> K >> nS;
        if (!(hs >> nLiq)) nLiq = 0;
    }
    if (K < 1 || nS != cfg.nSpecies || nLiq < 0 || nLiq > 1) {
        std::cerr << "[transport-probe] bad states file " << statesPath << " (K=" << K << ", nSpecies=" << nS
                  << ", config nSpecies=" << cfg.nSpecies << ", nLiq=" << nLiq << ")" << std::endl;
        return 2;
    }
    const GasPhaseLiquid liq = gasPhaseLiquid(cfg, var);
    const bool hasLiq = (nLiq == 1);
    if (hasLiq && liq.iw < 0) {
        std::cerr << "[transport-probe] liquid column given but the run is not a TP-carrier condensation run" << std::endl;
        return 2;
    }
    std::vector<double> sT(K), sRo(K), sY(static_cast<size_t>(K)*nS), sG(hasLiq ? K : 0);
    for (int k = 0; k < K; ++k) {
        in >> sT[k] >> sRo[k];
        for (int s = 0; s < nS; ++s) in >> sY[static_cast<size_t>(k)*nS + s];
        if (hasLiq) in >> sG[k];
    }
    if (!in) { std::cerr << "[transport-probe] states file truncated" << std::endl; return 2; }
    const TransportTableD* tt = thermo_transport_table();
    if (tt == nullptr || cfg.viscMethod != 2) {
        std::cerr << "[transport-probe] requires physProp.transport and viscMethod 2" << std::endl;
        return 2;
    }
    const int nR = tt->nReal;
    const geom_int nC = msh.nCells, nA = msh.nCells_all, nG = nA - nC;
    flow_float** roYdev = species_roY_device_ptr();
    const bool hasRoY = (roYdev != nullptr && nS >= 2);
    std::vector<flow_float*> roYptr(nS, nullptr);
    if (hasRoY) for (int s = 0; s < nS; ++s) roYptr[s] = var.c_d["roY" + std::to_string(s)];

    geom_int nWall = 0;
    for (const auto& bc : msh.bconds) nWall += static_cast<geom_int>(bc.iCells.size());
    // 境界面の解点は bcond をまたいで重複しうる (角のセル) ので、重複を除いた列に状態を割り当てる
    // (重複を数えると後の割り当てが先を上書きし、状態点が境界に乗らないことがある)。
    std::vector<geom_int> wallCellsU;
    {
        std::vector<char> seen(nA, 0);
        for (const auto& bc : msh.bconds)
            for (geom_int ic : bc.iCells) if (!seen[ic]) { seen[ic] = 1; wallCellsU.push_back(ic); }
    }
    const geom_int nWallU = static_cast<geom_int>(wallCellsU.size());
    const geom_int M = std::max<geom_int>(1, std::min<geom_int>(nWallU, std::max<geom_int>(nG, 1)));
    // FORGE_TRANSPORT_PROBE_NOCELLS=1: セル・ghost・壁の pass を省き、状態表の (5)(6) だけ評価する (大量の状態点の照合用)
    const char* noCellsEnv = getenv("FORGE_TRANSPORT_PROBE_NOCELLS");
    const bool noCells = (noCellsEnv != nullptr && std::string(noCellsEnv) == "1");
    const int P = noCells ? 0 : static_cast<int>((K + M - 1)/M);

    const std::string base(statesPath);
    std::ofstream bin(base + ".bin", std::ios::binary);
    auto w = [&](const void* p, size_t n) { bin.write(static_cast<const char*>(p), static_cast<std::streamsize>(n)); };

    double *mu_d = nullptr, *lam_d = nullptr, *X_d = nullptr;
    gpuErrchk(cudaMalloc(&mu_d, nA*sizeof(double)));
    gpuErrchk(cudaMalloc(&lam_d, nA*sizeof(double)));
    gpuErrchk(cudaMalloc(&X_d, static_cast<size_t>(nA)*nR*sizeof(double)));
    std::vector<size_t> wallCounts;
    for (int p = 0; p < P; ++p) {
        const geom_int o = static_cast<geom_int>(p)*M;
        std::vector<int> st(nA);
        for (geom_int i = 0; i < nC; ++i) st[i] = static_cast<int>((i + o) % K);
        for (geom_int g = 0; g < nG; ++g) st[nC + g] = static_cast<int>((g + o) % K);
        for (geom_int q = 0; q < nWallU; ++q) st[wallCellsU[q]] = static_cast<int>((q + o) % K);
        std::vector<flow_float> hT(nA), hRo(nA), hY(static_cast<size_t>(nA)*(hasRoY ? nS : 0)), hG(hasLiq ? nA : 0);
        for (geom_int i = 0; i < nA; ++i) {
            const int k = st[i];
            hT[i]  = static_cast<flow_float>(sT[k]);
            hRo[i] = static_cast<flow_float>(sRo[k]);
            if (hasLiq) hG[i] = static_cast<flow_float>(sRo[k]*sG[k]);
            if (hasRoY)
                for (int s = 0; s < nS; ++s)
                    hY[static_cast<size_t>(s)*nA + i] = static_cast<flow_float>(sRo[k]*sY[static_cast<size_t>(k)*nS + s]);
        }
        gpuErrchk(cudaMemcpy(var.c_d["T"], hT.data(), nA*sizeof(flow_float), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(var.c_d["ro"], hRo.data(), nA*sizeof(flow_float), cudaMemcpyHostToDevice));
        if (hasRoY)
            for (int s = 0; s < nS; ++s)
                gpuErrchk(cudaMemcpy(roYptr[s], hY.data() + static_cast<size_t>(s)*nA, nA*sizeof(flow_float), cudaMemcpyHostToDevice));
        if (hasLiq) gpuErrchk(cudaMemcpy(var.c_d["rog_0"], hG.data(), nA*sizeof(flow_float), cudaMemcpyHostToDevice));
        // 出力を NaN で埋めてから評価 (未更新の要素を検出できるように)
        {
            std::vector<flow_float> nanf(nA, std::numeric_limits<flow_float>::quiet_NaN());
            gpuErrchk(cudaMemcpy(var.c_d["vis_lam"], nanf.data(), nA*sizeof(flow_float), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemcpy(var.c_d["thermCond"], nanf.data(), nA*sizeof(flow_float), cudaMemcpyHostToDevice));
            gpuErrchk(cudaMemset(mu_d, 0xff, nA*sizeof(double)));
            gpuErrchk(cudaMemset(lam_d, 0xff, nA*sizeof(double)));
            gpuErrchk(cudaMemset(X_d, 0xff, static_cast<size_t>(nA)*nR*sizeof(double)));
        }
        gasProperties_d_wrapper(cfg, cuda_cfg, msh, var);                     // 実際のセル経路 (ghost 込み)
        gasPropertiesTransportProbe_d_wrapper(cfg, cuda_cfg, msh, var, mu_d, lam_d, X_d);
        WmlesTransportProbeOut wo;
        wmlesTransportProbe(cfg, cuda_cfg, msh, var, wo);                      // 壁モデルの物性関数

        // 実際に device にある入力と出力を読み戻して書く
        std::vector<flow_float> rT(nA), rRo(nA), rY(static_cast<size_t>(nA)*(hasRoY ? nS : 0)), rV(nA), rK(nA);
        std::vector<double> rMu(nA), rLam(nA), rX(static_cast<size_t>(nA)*nR);
        gpuErrchk(cudaMemcpy(rT.data(), var.c_d["T"], nA*sizeof(flow_float), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(rRo.data(), var.c_d["ro"], nA*sizeof(flow_float), cudaMemcpyDeviceToHost));
        if (hasRoY)
            for (int s = 0; s < nS; ++s)
                gpuErrchk(cudaMemcpy(rY.data() + static_cast<size_t>(s)*nA, roYptr[s], nA*sizeof(flow_float), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(rV.data(), var.c_d["vis_lam"], nA*sizeof(flow_float), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(rK.data(), var.c_d["thermCond"], nA*sizeof(flow_float), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(rMu.data(), mu_d, nA*sizeof(double), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(rLam.data(), lam_d, nA*sizeof(double), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(rX.data(), X_d, rX.size()*sizeof(double), cudaMemcpyDeviceToHost));
        w(st.data(), st.size()*sizeof(int));
        w(rT.data(), rT.size()*sizeof(flow_float));
        w(rRo.data(), rRo.size()*sizeof(flow_float));
        w(rY.data(), rY.size()*sizeof(flow_float));
        w(rV.data(), rV.size()*sizeof(flow_float));
        w(rK.data(), rK.size()*sizeof(flow_float));
        w(rMu.data(), rMu.size()*sizeof(double));
        w(rLam.data(), rLam.size()*sizeof(double));
        w(rX.data(), rX.size()*sizeof(double));
        const size_t nw = wo.cell.size();
        wallCounts.push_back(nw);
        w(wo.bcond.data(), nw*sizeof(int));
        w(wo.cell.data(), nw*sizeof(long long));
        w(wo.mu_w.data(), nw*sizeof(flow_float));
        w(wo.lam_w.data(), nw*sizeof(flow_float));
        w(wo.mu.data(), nw*sizeof(double));
        w(wo.lam.data(), nw*sizeof(double));
        if (hasLiq) {
            // 液 (読み戻し)・P (拡散係数の入力)・分子拡散係数 D_s (Fick 拡散と同じ組成・関数)
            std::vector<flow_float> rG(nA), rP(nA);
            std::vector<float> rD(static_cast<size_t>(nA)*nS, std::numeric_limits<float>::quiet_NaN());
            float* D_d = nullptr;
            gpuErrchk(cudaMalloc(&D_d, rD.size()*sizeof(float)));
            gpuErrchk(cudaMemcpy(D_d, rD.data(), rD.size()*sizeof(float), cudaMemcpyHostToDevice));
            speciesDmixProbe_d_wrapper(cfg, cuda_cfg, msh, var, D_d);
            gpuErrchk(cudaMemcpy(rD.data(), D_d, rD.size()*sizeof(float), cudaMemcpyDeviceToHost));
            cudaFree(D_d);
            gpuErrchk(cudaMemcpy(rG.data(), var.c_d["rog_0"], nA*sizeof(flow_float), cudaMemcpyDeviceToHost));
            gpuErrchk(cudaMemcpy(rP.data(), var.c_d["P"], nA*sizeof(flow_float), cudaMemcpyDeviceToHost));
            w(rG.data(), rG.size()*sizeof(flow_float));
            w(rP.data(), rP.size()*sizeof(flow_float));
            w(rD.data(), rD.size()*sizeof(float));
        }
    }
    cudaFree(mu_d); cudaFree(lam_d); cudaFree(X_d);

    // (5) double の Y を直接
    {
        double *Yd = nullptr, *Td = nullptr, *md = nullptr, *ld = nullptr, *Xd = nullptr, *Gd = nullptr;
        if (hasLiq) {
            gpuErrchk(cudaMalloc(&Gd, K*sizeof(double)));
            gpuErrchk(cudaMemcpy(Gd, sG.data(), K*sizeof(double), cudaMemcpyHostToDevice));
        }
        gpuErrchk(cudaMalloc(&Yd, sY.size()*sizeof(double)));
        gpuErrchk(cudaMalloc(&Td, K*sizeof(double)));
        gpuErrchk(cudaMalloc(&md, K*sizeof(double)));
        gpuErrchk(cudaMalloc(&ld, K*sizeof(double)));
        gpuErrchk(cudaMalloc(&Xd, static_cast<size_t>(K)*nR*sizeof(double)));
        gpuErrchk(cudaMemcpy(Yd, sY.data(), sY.size()*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(Td, sT.data(), K*sizeof(double), cudaMemcpyHostToDevice));
        gasPropertiesTransportProbeStates_d_wrapper(cfg, K, Yd, Gd, liq.iw, Td, md, ld, Xd);
        std::vector<double> hm(K), hl(K), hx(static_cast<size_t>(K)*nR);
        gpuErrchk(cudaMemcpy(hm.data(), md, K*sizeof(double), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(hl.data(), ld, K*sizeof(double), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(hx.data(), Xd, hx.size()*sizeof(double), cudaMemcpyDeviceToHost));
        w(hm.data(), hm.size()*sizeof(double));
        w(hl.data(), hl.size()*sizeof(double));
        w(hx.data(), hx.size()*sizeof(double));
        cudaFree(Yd); cudaFree(Td); cudaFree(md); cudaFree(ld); cudaFree(Xd); cudaFree(Gd);
    }
    // (6) 表引き (#5t2-3) があれば、同じ状態表の Y・T を float に丸めて表引きの経路を呼ぶ (float の μ・λ)
    const bool hasTab = (tt->tab.valid != 0);
    if (hasTab) {
        double *Yd = nullptr, *Td = nullptr, *Gd = nullptr;
        float *md = nullptr, *ld = nullptr;
        if (hasLiq) {
            gpuErrchk(cudaMalloc(&Gd, K*sizeof(double)));
            gpuErrchk(cudaMemcpy(Gd, sG.data(), K*sizeof(double), cudaMemcpyHostToDevice));
        }
        gpuErrchk(cudaMalloc(&Yd, sY.size()*sizeof(double)));
        gpuErrchk(cudaMalloc(&Td, K*sizeof(double)));
        gpuErrchk(cudaMalloc(&md, K*sizeof(float)));
        gpuErrchk(cudaMalloc(&ld, K*sizeof(float)));
        gpuErrchk(cudaMemcpy(Yd, sY.data(), sY.size()*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(Td, sT.data(), K*sizeof(double), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemset(md, 0xff, K*sizeof(float)));
        gpuErrchk(cudaMemset(ld, 0xff, K*sizeof(float)));
        gasPropertiesTransportProbeStatesTab_d_wrapper(cfg, K, Yd, Gd, liq.iw, Td, md, ld);
        std::vector<float> hm(K), hl(K);
        gpuErrchk(cudaMemcpy(hm.data(), md, K*sizeof(float), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(hl.data(), ld, K*sizeof(float), cudaMemcpyDeviceToHost));
        w(hm.data(), hm.size()*sizeof(float));
        w(hl.data(), hl.size()*sizeof(float));
        cudaFree(Yd); cudaFree(Td); cudaFree(md); cudaFree(ld); cudaFree(Gd);
    }
    bin.close();

    std::ofstream js(base + ".json");
    js << "{\"K\": " << K << ", \"nSpecies\": " << nS << ", \"nReal\": " << nR << ", \"nCells\": " << nC
       << ", \"nCells_all\": " << nA << ", \"hasRoY\": " << (hasRoY ? 1 : 0) << ", \"passes\": " << P
       << ", \"sizeof_flow_float\": " << sizeof(flow_float) << ", \"table\": " << (hasTab ? 1 : 0)
       << ", \"hasLiq\": " << (hasLiq ? 1 : 0) << ", \"iw\": " << liq.iw << ", \"wall_counts\": [";
    for (size_t i = 0; i < wallCounts.size(); ++i) js << (i ? ", " : "") << wallCounts[i];
    js << "], \"bconds\": [";
    for (size_t b = 0; b < msh.bconds.size(); ++b)
        js << (b ? ", " : "") << "{\"physID\": " << msh.bconds[b].physID << ", \"kind\": \"" << msh.bconds[b].bcondKind << "\"}";
    js << "]}" << std::endl;
    std::cout << "[transport-probe] wrote " << base << ".bin/.json (K=" << K << ", passes=" << P << ", nCells=" << nC
              << ", nCells_all=" << nA << ", boundary planes=" << nWall << ", nReal=" << nR << ")" << std::endl;
    return 0;
}

// 試験用 (FORGE_TRANSPORT_TABLE_PROBE=<points.txt>; plan thermophysics-solver-owned-species-db #5t2-3、tests/unit/test_transport_gpu.py)。
//   輸送表の単体値を GPU で評価して終了する (時間更新なし)。points.txt: 1 行目 N、以降 N 行 "kind idx T"
//   (kind 0: 実種 idx の μ・λ、kind 1: 組 idx の η; T は float に丸めて使う)。出力: <points.txt>.bin (float v0[N], v1[N]) と
//   <points.txt>.json (表の範囲・各表の分割区間 [Ta, Tb, 小区間数, 所属の閾値 Tupper]・組の種類・メモリ量)。
static int runTransportTableProbe(const char* path)
{
    const TransportTableD* tt = thermo_transport_table();
    const TransportTablesHost* H = thermo_transport_tables_host();
    if (tt == nullptr || H == nullptr) {
        std::cerr << "[transport-table-probe] requires physProp.transport with tables enabled" << std::endl;
        return 2;
    }
    std::ifstream in(path);
    long long N = 0;
    if (!(in >> N) || N < 0) { std::cerr << "[transport-table-probe] bad points file " << path << std::endl; return 2; }
    std::vector<int> kind(N), idx(N);
    std::vector<float> T(N);
    for (long long k = 0; k < N; ++k) { double t; in >> kind[k] >> idx[k] >> t; T[k] = static_cast<float>(t); }
    if (!in) { std::cerr << "[transport-table-probe] points file truncated" << std::endl; return 2; }
    std::vector<float> v0(N, 0.0f), v1(N, 0.0f);
    if (N > 0) {
        int *kd = nullptr, *id = nullptr;
        float *Td = nullptr, *a = nullptr, *b = nullptr;
        gpuErrchk(cudaMalloc(&kd, N*sizeof(int)));
        gpuErrchk(cudaMalloc(&id, N*sizeof(int)));
        gpuErrchk(cudaMalloc(&Td, N*sizeof(float)));
        gpuErrchk(cudaMalloc(&a, N*sizeof(float)));
        gpuErrchk(cudaMalloc(&b, N*sizeof(float)));
        gpuErrchk(cudaMemcpy(kd, kind.data(), N*sizeof(int), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(id, idx.data(), N*sizeof(int), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(Td, T.data(), N*sizeof(float), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemset(a, 0xff, N*sizeof(float)));
        gpuErrchk(cudaMemset(b, 0xff, N*sizeof(float)));
        transportTableSingles_d_wrapper(static_cast<int>(N), kd, id, Td, a, b);
        gpuErrchk(cudaMemcpy(v0.data(), a, N*sizeof(float), cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(v1.data(), b, N*sizeof(float), cudaMemcpyDeviceToHost));
        cudaFree(kd); cudaFree(id); cudaFree(Td); cudaFree(a); cudaFree(b);
    }
    const std::string base(path);
    {
        std::ofstream bin(base + ".bin", std::ios::binary);
        bin.write(reinterpret_cast<const char*>(v0.data()), static_cast<std::streamsize>(N*sizeof(float)));
        bin.write(reinterpret_cast<const char*>(v1.data()), static_cast<std::streamsize>(N*sizeof(float)));
    }
    std::ofstream js(base + ".json");
    js << std::setprecision(17);
    const int nR = tt->nReal;
    js << "{\"N\": " << N << ", \"nReal\": " << nR << ", \"Tmin\": " << static_cast<double>(H->Tmin)
       << ", \"Tmax\": " << static_cast<double>(H->Tmax) << ", \"bytes\": " << H->bytes() << ", \"tables\": [";
    const int nTab = static_cast<int>(H->segLo.size());
    for (int t = 0; t < nTab; ++t) {
        const TransportTabRefF& r = (t < nR) ? H->spTab[t] : H->pairTab[t - nR];
        js << (t ? ", " : "") << "{\"kind\": " << (t < nR ? 0 : 1) << ", \"idx\": " << (t < nR ? t : t - nR);
        if (t >= nR) {
            const TransportPairD& p = speciesDB_current()->transport.pairs[t - nR];
            js << ", \"a\": " << p.a << ", \"b\": " << p.b << ", \"pair_kind\": " << p.kind;
        }
        js << ", \"segs\": [";
        for (int k = 0; k < r.nseg; ++k) {
            const TransportSegF& g = H->seg[r.seg0 + k];
            js << (k ? ", " : "") << "[" << H->segLo[t][k] << ", " << H->segHi[t][k] << ", " << g.m << ", "
               << static_cast<double>(g.Tupper) << "]";
        }
        js << "]}";
    }
    js << "]}" << std::endl;
    std::cout << "[transport-table-probe] wrote " << base << ".bin/.json (N=" << N << ", tables=" << nTab
              << ", bytes=" << H->bytes() << ")" << std::endl;
    return 0;
}

// forge --resolve-species: GPU を使わず solverConfig.yaml を読み、解決済み記録を cwd に書いて互換性ハッシュを標準出力の最終行に出す。
// IC 生成・runner が宛先の物性を得るため (#3b)。記録を既存場へ貼るだけで検証済みにはしない。CPG は終了コード 2。
static int resolveSpeciesOnly()
{
    std::streambuf* orig = std::cout.rdbuf(std::cerr.rdbuf());   // ログは stderr へ、stdout はハッシュだけ
    int rc = 0;
    std::string hash;
    try {
        solverConfig cfg;
        cfg.read("solverConfig.yaml");
        if (cfg.thermalMethod != 2) {
            std::cerr << "[species] thermalMethod != 2 (calorically perfect): no species record / hash" << std::endl;
            rc = 2;
        } else {
            const ResolvedSpeciesDB& db = speciesDB_init(cfg);
            speciesDB_printTable(cfg, db);
            const SpeciesRecordInfo rec = speciesDB_writeRecord(db, speciesRecordTref(cfg), cfg.speciesDBFile, "",
                                                                "not_checked_resolve_only", 0, ".");
            std::cerr << "[species] record " << rec.recordFile << " (sha256 " << rec.recordSha256 << ")" << std::endl;
            hash = rec.compatHash;
        }
    } catch (const std::exception& e) {
        std::cerr << "[species] resolve failed: " << e.what() << std::endl;
        rc = 1;
    }
    std::cout.rdbuf(orig);
    if (rc == 0) std::cout << hash << std::endl;
    return rc;
}

cudaConfig initializeSimulation(
    solverConfig& cfg,
    mesh& msh,
    matrix& mat_ns,
    variables& var,
    fluct_variables& fluct,
    point_probes& pprobes)
{
    MEMLOG_SOLVER("start");
    cout << "Read Solver Config \n";
    cfg.read("solverConfig.yaml");
    appendLaunchRecord(cfg);
    MEMLOG_SOLVER("after cfg.read");

    // 化学種 DB の host 側解決 (GPU 非依存; 未知種名はここで exit)。bcond の X{s}→Y{s} 換算と
    // 起動ログ (種表) が使う。thermo_init_db は同じ結果を device へ上げる。
    speciesDB_printTable(cfg, speciesDB_init(cfg));
    // physProp.transport (種ごとの輸送物性の出所): 段 2 (plan thermophysics-solver-owned-species-db #5t2-2) で GPU に接続した。
    // セル (gasProperties_d) と壁 (wmlesWallModel_d) の viscMethod 2 経路だけが新しい μ・λ を使う。viscMethod 0/1 では
    // μ・λ は定数・Sutherland のままで記録 (輸送ブロック) と実際の計算が食い違うので、その組み合わせは起動を止める。
    if (!cfg.speciesTransport.empty() && cfg.viscMethod != 2) {
        cerr << "[transport] ERROR: physProp.transport is used only with viscMethod: 2 (the GPU transport path replaces the "
                "kinetic-theory mixture). With viscMethod " << cfg.viscMethod << " the recorded transport would not be the one "
                "computed. Set viscMethod: 2 or remove physProp.transport." << endl;
        std::exit(EXIT_FAILURE);
    }
    // viscMethod 2 は種ごとの輸送物性 (physProp.transport 必須、CEA 形 frozen 混合則) に置き換えた (plan §4.3c 案 C)。
    // 旧 kinetic 経路 (LJ + Chapman-Enskog、Wilke の φ を μ と λ で共用) は計算から外したので、指定なしは起動を止める。
    // 旧結果の再現は旧バイナリで行う (procedures/solver-settings.md)。
    if (cfg.viscMethod == 2 && cfg.speciesTransport.empty()) {
        cerr << "[transport] ERROR: viscMethod: 2 requires physProp.transport (the transport model of every real species, "
                "e.g. physProp: {thermalMethod: 2, viscMethod: 2, transport: {N2: cea, O2: cea, H2O: custom:h2o_iapws_cea_v1}}; "
                "models: cea, kinetic, fit, custom:<name>_v<version>; thermalMethod: 2 only). The former kinetic-theory "
                "mixture (Wilke phi shared by mu and lambda) has been removed. If air Sutherland viscosity is enough, use "
                "viscMethod: 1. To reproduce old viscMethod 2 results, run the old binary." << endl;
        std::exit(EXIT_FAILURE);
    }

    cout << "Init Thermo DB \n";
    thermo_init_db(cfg);   // NASA-9/LJ 化学種 DB を構築し device へアップロード (thermalMethod==2 用)
    chemistry_init(cfg);   // 有限速度化学: 反応機構を読み device へ (chemistry.enabled==1 のみ)
    MEMLOG_SOLVER("after thermo/chemistry init");

    cout << "Read Mesh \n";
    if (cfg.meshFormat == "hdf5") {
        msh.nodeValueAtNode = (cfg.discretization == "node") ? 1 : 0;   // node は常に値=ノード座標
        msh.readMesh(cfg.meshFileName);
    } else {
        cerr << "Error unknown mesh format: " << cfg.meshFormat << endl;
        std::exit(EXIT_FAILURE);
    }

    MEMLOG_SOLVER("after readMesh");
    // mat_ns のメンバ (structure/lhs/rhs/localPlnOfCell) を読むのはビルド対象外の solvePoisson_amgx.cpp・
    // solvePoisson_amgcl.cpp だけで、GPU 経路は各 wrapper へ参照を渡すだけ。gpu: 1 では作らない
    // (約 312 B/節点。plan architecture-solver-host-memory §4.1)。オブジェクト自体は引数として渡るので残す。
    if (cfg.gpu == 0) {
        cout << "Init Matrix (but not used now) \n";
        mat_ns.initMatrix(msh);
    } else {
        cout << "Init Matrix: skipped (gpu: 1 does not use it) \n";
    }
    MEMLOG_SOLVER("after initMatrix");

    cout << "Read Boundary Conditions \n";
    readBcondConfig(cfg , msh.bconds);
    MEMLOG_SOLVER("after readBcondConfig");

    // 厚さ 0 の板の自由端の近傍の速度再構成 (space.zeroThicknessEdgeVelocity、既定 無効)。bcond 種別 (周期の検査) が
    // 揃った直後に、契約とメッシュ h5 の /AUX の重み w を検査して device へ上げる (無効なら起動エコー 1 行だけで何もしない)。
    // plans/active/convection-zero-thickness-edge-reconstruction.md §4
    {
        ZteLaunchInfo zinfo;
        zteSetVelocityWeightDevice(zteLoadVelocityWeight(cfg , msh , &zinfo));
        appendLaunchRecordZte(zinfo);   // 照合済みの値を起動記録 forge_launches.jsonl に追記 (無効は enabled 0)
    }

    // 凝縮セルの二相 frozen 音速 (condSonicModel) の自動解決: bcond 種別が揃った後に確定し理由をログに出す
    // (plans/active/condensation-kantrowitz-gamma-twophase-sonic.md §4.2)。
    {
        std::vector<std::string> kinds;
        for (const auto& bc : msh.bconds) kinds.push_back(bc.bcondKind);
        std::string reason, warn;
        const int resolved = resolveCondSonicModel(cfg.condSonicModel, cfg.condensation, cfg.thermalMethod, cfg.condGasSpecies,
                                                   cfg.condModel, cfg.condEquilibrium, kinds, reason, warn);
        if (cfg.condensation == 1) {
            cout << "[condensation] condSonicModel=" << resolved << " (" << reason << ")\n";
            if (!warn.empty()) cout << "[condensation] WARNING: " << warn << "\n";
            // 読込後の実効値 (省略時既定の確認用; codex 2026-09-13 carrier result M3): condKantrowitz 0=補正なし/1=Kantrowitz/2,3=Feder carrier
            // dual-time は受動種経路 (passiveScalarScheme 1) ならモーメントに BDF 物理時間項が入る (plan species-passive-scalar-unification §4.4)
            // ので降格しない。旧経路 0 は従来どおり降格 (物理時間項なし)。
            if (cfg.condLimiterMode == 1 && (cfg.timeIntegration != 11 || (cfg.dualTime != 0 && cfg.passiveScalarScheme == 0))) {
                // 更新クランプ経路は定常 point-implicit (timeIntegration 11, dual-time なし) だけで検証済み。RK 陽解法は未制限残差の
                // 累積バッファを持ち、dual-time は凝縮モーメントに物理時間項が無い (followups F-cf8) ので旧経路に降格する。
                cout << "[condensation] condLimiterMode 1 is verified for steady timeIntegration 11 only; falling back to 0 (legacy residual theta) for this run\n";
                cfg.condLimiterMode = 0;
            }
            cout << "[condensation] condLimiterMode=" << cfg.condLimiterMode
                 << (cfg.condLimiterMode == 1 ? " (theta clamps the moment update; residual source is pseudo-time-step independent)" : " (legacy: theta multiplies the residual source)")
                 << " condDgMaxStep=" << cfg.condDgMaxStep << " condDTmaxStep=" << cfg.condDTmaxStep << "\n";
            cout << "[condensation] condKantrowitz=" << cfg.condKantrowitz << " condKantrowitzGammaMode=" << cfg.condKantrowitzGammaMode
                 << " condSigmaScale=" << cfg.condSigmaScale << " condVaporMassFraction=" << cfg.condVaporMassFraction
                 << " condN2LatentLowT=" << cfg.condN2LatentLowT << " condN2PsatLowT=" << cfg.condN2PsatLowT << " condN2LiquidCp=" << cfg.condN2LiquidCp << "\n";
        }
        cfg.condSonicModel = resolved;
        // 二相拡散 (condTwoPhaseDiffusion, #4e): 実効状態 (active / inactive-a / inactive-b / unsupported-c) を判定して cfg に記録し
        // [twophase] 行を出す。指定 ON + 未対応 (c) は理由を出して終了 (plan condensation-two-phase-default §4-2)。周期は bcond 種別で判定。
        condTwoPhaseDiffusionValidate(cfg, std::find(kinds.begin(), kinds.end(), std::string("periodic")) != kinds.end());
        // CPG carrier 形 (空気の N2 選択凝縮) の境界受付範囲 (plans/accepted/condensation-air.md §4.1; codex 2026-09-13 M3):
        //   ghost/ピンを単相 EOS で再構成する境界 (wall, wall_isothermal, inlet_Pressure, outflow, periodic 等) は未対応。
        //   出口 outlet_statPress は超音速全量外挿のときだけ整合 (実行時条件) → 後処理 onset_analysis.py --series の u_n/c>1 で確認する。
        if (cfg.condensation == 1 && cfg.condVaporMassFraction > 0.0) {
            for (const auto& k : kinds) {
                if (!(k == "inlet_uniformVelocity" || k == "outlet_statPress" || k == "slip")) {
                    std::cerr << "Configuration Error: condVaporMassFraction (CPG carrier) supports boundary kinds inlet_uniformVelocity / outlet_statPress (supersonic) / slip only; found '" << k << "'.\n";
                    std::exit(EXIT_FAILURE);
                }
            }
            cout << "[condensation] CPG carrier (air): Y_w=" << cfg.condVaporMassFraction << ", boundaries verified kinds only; outlet must stay supersonic (checked in post-processing)\n";
        }
    }

    // 入口分布プロファイル (ints:{inletProfile:1} の inlet): per-face bvar を CSV から face 重心で補間。
    // mesh.planes (重心) と bvar が揃った後・最初の applyBconds より前に適用。未指定 inlet は一様のまま。
    applyInletProfiles(cfg , msh);
    // 壁温分布 (ints: {wallProfile: 1})。inlet と同じ位相で、最初の applyBconds より前に入れる。
    applyWallProfiles(cfg , msh);
    // 壁温を陽に扱う run では、共有 CV に矛盾する壁温を与えていないかを検査する (後勝ち事故の防止)。
    conjugateWall::checkWallTemperatureSharing(cfg , msh);
    // ソルバ内 CHT (conjugate: ブロック + bcond ints: {conjugate: 1}) の対応範囲を検査する。
    conjugateWall::initConjugateWalls(cfg , msh);
    MEMLOG_SOLVER("after inlet/wall profiles, conjugate init");

    // 化学種変数を登録 (allocVariables より前)。nSpecies<=1 では no-op。
    var.registerSpecies(cfg.nSpecies, cfg.chemEnabled);

    // 非平衡凝縮モーメント変数を登録 (allocVariables より前)。condensation==0 では no-op。
    var.registerCondensation(cfg.nCondSpecies);
    var.registerTwoPhaseVaporResidual(condTwoPhaseDiffusionActive(cfg) ? 1 : 0);   // 二相拡散の蒸気残差 res_roYv (#4e)

    // 受動トレーサ roXi を登録 (allocVariables より前)。physProp.tracer 未指定では no-op。
    var.registerTracer(cfg.tracerEnabled() ? 1 : 0);

    // 遷移モデル (turbulence.transition: lm2009) の変数を登録。none では no-op。診断場は output.level>=2 のときだけ確保。
    transitionValidateConfig(cfg);
    var.registerTransition(cfg.transitionEnabled() ? 1 : 0, (cfg.outputLevel >= 2) ? 1 : 0);

    // 変数・診断の登録 → 出力と checkpoint が要る名前の確定 → H の構築 → 確保 (plan architecture-solver-host-memory §4・§4.3)。
    // 環境変数による出力名の追加 (旧: main() の確保の後) と無効な診断変数の除去 (旧: allocVariables の先頭) を確保の前に済ませ、
    // 出力側と同じ関数 (output/outputFieldNames.cpp) でホスト確保集合 H を作る。gpu: 1 ではホストのセル配列を H の名前だけ確保する。
    registerOutputDiagnostics(var);
    applyEnvGatedRemovals(var);
    const std::set<std::string> hostCells = hostCellSet(cfg, var);
    if (cfg.gpu == 1) {
        std::ostringstream os;
        os << "[variables] host cell arrays (gpu: 1): " << hostCells.size() << " of " << var.cellValNames.size() << " registered:";
        for (const auto& n : hostCells) os << " " << n;
        cout << os.str() << "\n";
    }
    var.allocVariables(cfg.gpu , msh, hostCells);
    MEMLOG_SOLVER("after allocVariables");

    // device roY[] ポインタ配列を構築 (c_d 確保後, dependentVariables より前)。
    speciesInit_d(cfg , var);

    // device rog[] ポインタ配列を構築 (二相 EOS が液相質量分率を読む)。condensation==0 で no-op。
    condensationInit_d(cfg , var);
    // 受動種 (トレーサ + 凝縮モーメント) の device ポインタ配列 (passiveScalarScheme 1 の化学種経路用; 0 では表だけ)。
    passiveInit_d(cfg , var);

    // 入力場の化学種照合と解決済み記録 (plans/active/thermophysics-solver-owned-species-db.md §4.3, #3a)。
    // TP (thermalMethod 2) だけが対象。数値には触れない (記録と照合のみ)。
    checkInputSpeciesAndWriteRecord(cfg);

    MEMLOG_SOLVER("after species/condensation/passive init_d");
    cout << "Read Initial Values \n";
    var.readValueHDF5(cfg.valueFileName , msh, cfg.kInit, cfg.omegaInit);
    MEMLOG_SOLVER("after readValueHDF5");

    cout << "Set mesh connection map for cuda \n";
    cudaConfig cuda_cfg(msh);
    msh.setMeshMap_d();
    MEMLOG_SOLVER("after setMeshMap_d (mesh maps H2D)");
    // node × 回転周期 (bcond periodic の type != 0) は**起動エラー** (plan boundary-node-rotational-periodic §5.1 #0a)。
    // node の seam 経路は速度・勾配・陰解法の dq を回さず、NS の periodicGradientGather も並進を前提に部分和を足すので、
    // 走らせると黙って誤った解になる。回転の実装 (同 plan #2–#4) が入るまで止める。周期対の構築 (失敗すると例外) より前に判定する。
    if (cfg.discretization == "node") {
        for (const auto& bc : msh.bconds) {
            if (bc.bcondKind != "periodic") continue;
            auto it = bc.inputInts.find("type");
            if (it != bc.inputInts.end() && it->second != 0) {
                std::cerr << "[config] node の回転周期 (bcond periodic, physID " << bc.physID << ", type " << it->second
                          << ") は未対応です (plan boundary-node-rotational-periodic)。回転の実装が入るまで起動しません。"
                          << "並進周期 (type 0) か cell を使ってください。" << std::endl;
                std::exit(1);
            }
        }
    }
    msh.setPeriodicPartner();
    // setPeriodicPartner は host の bint["partnerCellID"/"partnerPlnID"] を埋めるが device (bint_d) へは
    // 転送されない (bcondInitVariables の H2D は yaml uniform=type==1 の bint のみ)。これを怠ると periodic_d が
    // 未初期化の device partnerCellID を読み、ghost が幾何的に誤ったセル値を取得して cell 周期境界が非保存になる
    // (node は host の partnerCellID を直接使うため無傷)。明示的に H2D 転送する。
    for (auto& bc : msh.bconds) {
        if (bc.bcondKind != "periodic") continue;
        for (const char* key : {"partnerCellID","partnerPlnID"}) {
            if (bc.bint.count(key)==0) continue;
            if (bc.bint_d.count(key)==0 || bc.bint_d[key]==nullptr)
                gpuErrchk( cudaMalloc(&(bc.bint_d[key]), bc.bint[key].size()*sizeof(geom_int)) );
            gpuErrchk( cudaMemcpy(bc.bint_d[key], bc.bint[key].data(),
                       bc.bint[key].size()*sizeof(geom_int), cudaMemcpyHostToDevice) );
        }
    }

    MEMLOG_SOLVER("after setPeriodicPartner");
    var.setStructuralVariables(cfg , cuda_cfg , msh);
    MEMLOG_SOLVER("after setStructuralVariables");
    // node-centered 周期境界 DOF 同一視 (median-dual M4, §4.5): partnerCellID から周期ノード group(union-find)
    // を構築し、各 group の合併体積を var.c_d["volume"] へ書き戻す。setDT より前 (volume を使うため)。
    // cell モード / 非周期では no-op。
    msh.buildPeriodicNodeGroups(cfg.discretization == "node", var.c_d["volume"]);
    // node 周期境界 DOF 同一視 (§4.5.9): 読み込んだ保存量を root→member ミラーで slave=master に初期同期する。
    // seed (2D複製+非周期 Uz 摂動 / restart / 丸め) で周期ペアの保存量が desync していることがあり、残差 gather だけ
    // では desync が永続して継ぎ目フラックス不整合 (seam 圧力欠陥) を生む。dependentVariables より前に揃える。
    periodicMirrorNSState_d_wrapper(cfg , cuda_cfg , msh , var);
    // SST の roK/roOmega も初期同期する (plan gradient-scalar-lsq-unification §4.2、codex plan M2)。従来は更新後
    // (point-implicit / 陰解法 sweep の直後) にしかミラーしておらず、読み込んだ場の desync が最初の勾配に入っていた。
    // 非 SST / cell / 非周期では no-op。
    periodicMirrorScalarState_d_wrapper(cfg , cuda_cfg , msh , var);
    speciesPrimitive_d_wrapper(cfg , cuda_cfg , msh , var);  // Y_s = ρY_s/ρ (roY を読込済)
    condensationPrimitive_d_wrapper(cfg , cuda_cfg , msh , var);  // φ = ρφ/ρ (液相モーメント読込済)
    tracerPrimitive_d_wrapper(cfg , cuda_cfg , msh , var);  // ξ = ρξ/ρ (トレーサ読込済)
    // 毎更新の EOS 床事象のカウンタ (output.floorEvents、既定 off): 記録を開き、初期化の EOS は更新事象と別枠 (init) で数える
    // (初期化では EOS の後に updateVariablesOuter が走るので、ここでの床は基準状態へ入る。plan tooling-sern-te-wake-grid §4)
    {
        std::ostringstream fnv; fnv << std::hex << fnv1a64File("solverConfig.yaml");
        floorEventsOpen(cfg, msh, fnv.str());
    }
    floorEventsSetPhase(FloorEventPhase::Init, 0);
    dependentVariables(cfg , cuda_cfg , msh , var, mat_ns);
    floorEventsSetPhase(FloorEventPhase::Aux, 0);
    // node-centered 壁 Dirichlet: IC の壁ノード速度を厳密 0 に初期化 (KE を roe から除去)。
    // この後 gasProperties が補正 roe から P/T を再計算する。cell/非 node では no-op。
    enforceWallNoSlip_d_wrapper(cfg , cuda_cfg , msh , var);
    gasProperties_d_wrapper(cfg , cuda_cfg , msh , var);

    MEMLOG_SOLVER("after initial EOS/gasProperties");
    fluct.allocVariables();
    fluct.set_fluctVelocity(cfg , cuda_cfg , msh , var);

    applyBconds(cfg , cuda_cfg , msh , var, mat_ns , fluct);
    applyRansScalarBoundaries(cfg , cuda_cfg , msh , var);
    // 遷移モデル: 初期出力 (res_0) の前に原始量の生成 (入力に無ければ初期化)・入口ピン・周期同期を済ませる
    // (codex result M1: 未初期化の roGamma/roReth/gammaTr/reTheta/gammaEff が初期出力に出ていた)。none では no-op。
    periodicMirrorTransitionState_d_wrapper(cfg , cuda_cfg , msh , var);
    transitionPrimitive_d_wrapper(cfg , cuda_cfg , msh , var);
    applyTransitionBoundaries(cfg , cuda_cfg , msh , var);
    applyWmlesWallModel(cfg , cuda_cfg , msh , var);   // WMLES 壁応力モデル (wallModelLES 壁のみ, §10)
    weakIsoWall::validate(cfg, msh);   // 弱形式 (nodeIsothermalEnergyBC=1) の構成検査 (併用不可を拒否)
    applyNodeIsothermalWallPin(cfg , cuda_cfg , msh , var);   // 素の node 等温壁の壁ノード T ピン
    applySstThermalWallFunction(cfg , cuda_cfg , msh , var);  // SST 熱的壁関数: 断熱壁 T_aw (§6.5(f))
    applySpeciesBoundaries(cfg , cuda_cfg , msh , var);
    applyCondensationBoundaries(cfg , cuda_cfg , msh , var);
    applyTracerBoundaries(cfg , cuda_cfg , msh , var);
    calcGradient_d_wrapper(cfg , cuda_cfg , msh , var);
    // 初期 setup でも周期勾配 gather を適用 (assembleResidual と整合; res_0 出力と初期診断を正しい合併勾配にする)。
    preGatherDumpMain(cfg , msh , var , "init");   // 診断 (FORGE_DUMP_PREGATHER、既定 off・出力専用)
    periodicGradientGather_d_wrapper(cfg , cuda_cfg , msh , var);
    axisymmetricGeomTerms_d_wrapper(cfg , cuda_cfg , msh , var);
    updateVariablesOuter(cfg , cuda_cfg , msh , var , mat_ns);
    // FP64 正本は **updateVariablesOuter の後**に「確保 → 初期化」を**隣接させて**置く (§5.1 S1a)。
    // ここより前 (周期ミラー :1221・初期ピン) に置くと、それらの初期射影が Qacc に入らない。
    // **離すと壊れる**: 以前は確保を main() 側に置いていたため初期化が先に走って no-op になり、
    // Qacc=0 のまま commit されて ro≈0 → step 4 で発散した (2026-09-23, run_0023_g1_on)。
    if (cfg.qAccumulatorFP64 == 1) {
        var.allocQAccumulator(msh.nCells);
        var.initQAccumulatorFromQ(msh.nCells);
    }
    speciesUpdateOuter_d_wrapper(cfg , cuda_cfg , msh , var);  // roY{s}N/M ベースライン
    condensationUpdateOuter_d_wrapper(cfg , cuda_cfg , msh , var);  // 液相モーメント N/M ベースライン
    tracerUpdateOuter_d_wrapper(cfg , cuda_cfg , msh , var);  // トレーサ N/M ベースライン
    initDualTimeHistory(cfg , cuda_cfg , msh , var);   // dual-time: 前物理レベルの復元 (checkpoint) または BDF1 再開 (§4.4)
    setDT_d_wrapper(cfg , cuda_cfg , msh , var);
    MEMLOG_SOLVER("after initial BC/gradient/updateOuter/setDT");

    pprobes.init(cfg , cuda_cfg , msh);
    MEMLOG_SOLVER("after pprobes.init");


    // リミッタの無次元化基準 (plan convection-node-wall-reconstruction §4.13)。**run 中固定**。
    // 初期場から体積加重平均で ro_ref / p_ref (絶対圧) / a_ref (音速) を取り、
    // 基準長は未指定ならメッシュ境界箱の対角 (格子細分で変わらない長さ)。
    // h_i は 2D / 軸対称では半径重み前の面積の平方根、3D では体積の立方根 (ieleType で判別)。
    // limiterDiag>0 のときも基準値を計算する: G1 の許容幅に**構成によらない絶対床**を与えるため
    // (Uy/Uz が恒等 0 の場では近傍レンジも 0 になり、float ノイズを逸脱と数えてしまう)。
    if (cfg.limiterScaled == 1 || cfg.limiterDiag > 0) {
        std::vector<flow_float> h_ro(msh.nCells), h_P(msh.nCells), h_a(msh.nCells);
        gpuErrchk( cudaMemcpy(h_ro.data(), var.c_d["ro"],    msh.nCells*sizeof(flow_float), cudaMemcpyDeviceToHost) );
        gpuErrchk( cudaMemcpy(h_P.data(),  var.c_d["P"],     msh.nCells*sizeof(flow_float), cudaMemcpyDeviceToHost) );
        gpuErrchk( cudaMemcpy(h_a.data(),  var.c_d["sonic"], msh.nCells*sizeof(flow_float), cudaMemcpyDeviceToHost) );
        double wro=0.0, wP=0.0, wa=0.0, wv=0.0;
        double xmin=1e300,xmax=-1e300,ymin=1e300,ymax=-1e300,zmin=1e300,zmax=-1e300;
        bool all2D = true;
        for (geom_int ic=0; ic<msh.nCells; ic++) {
            const double v = msh.cells[ic].volume;
            wv += v; wro += v*std::fabs((double)h_ro[ic]); wP += v*std::fabs((double)h_P[ic]); wa += v*std::fabs((double)h_a[ic]);
        }
        for (const auto& nd : msh.nodes) {
            if (nd.coords.size() < 3) continue;
            xmin=std::min(xmin,(double)nd.coords[0]); xmax=std::max(xmax,(double)nd.coords[0]);
            ymin=std::min(ymin,(double)nd.coords[1]); ymax=std::max(ymax,(double)nd.coords[1]);
            zmin=std::min(zmin,(double)nd.coords[2]); zmax=std::max(zmax,(double)nd.coords[2]);
        }
        // 明示指定 (>0) は上書きしない。自動だと restart ごとに値が変わり、分割実行が
        // 連続実行と同じ数値作用素にならない (codex plan-3 Major 6)。
        const bool roAuto = !(cfg.limiterRoRef > 0.0), pAuto = !(cfg.limiterPRef > 0.0), aAuto = !(cfg.limiterARef > 0.0);
        if (roAuto) cfg.limiterRoRef = (wv>0.0 && wro>0.0) ? wro/wv : 1.0;
        if (pAuto)  cfg.limiterPRef  = (wv>0.0 && wP >0.0) ? wP /wv : 1.0;
        if (aAuto)  cfg.limiterARef  = (wv>0.0 && wa >0.0) ? wa /wv : 1.0;
        if (!(cfg.limiterRefLength > 0.0)) {
            const double dx=xmax-xmin, dy=ymax-ymin, dz=zmax-zmin;
            const double diag = std::sqrt(dx*dx+dy*dy+dz*dz);
            cfg.limiterRefLength = (diag > 0.0) ? diag : 1.0;
        }
        // 2D 判定は**境界箱が 1 方向に潰れているか**で行う (平面 2D は z 幅 0)。
        // ieleType は node の双対 CV では primal の型を持たないので使えない。
        {
            const double dx=xmax-xmin, dy=ymax-ymin, dz=zmax-zmin;
            const double dmax=std::max(dx,std::max(dy,dz));
            all2D = (dmax > 0.0) && (std::min(dx,std::min(dy,dz)) < 1.0e-6*dmax);
        }
        cfg.limiterLengthFromArea = (all2D || cfg.isAxisymmetric == 1) ? 1 : 0;
        std::cout << "[limiter] scaled: ro_ref=" << cfg.limiterRoRef << (roAuto ? "(auto)" : "(fixed)")
                  << " p_ref=" << cfg.limiterPRef << (pAuto ? "(auto)" : "(fixed)")
                  << " a_ref=" << cfg.limiterARef << (aAuto ? "(auto)" : "(fixed)")
                  << " L_ref=" << cfg.limiterRefLength
                  << " K=" << cfg.venkatK << " h_i=" << (cfg.limiterLengthFromArea ? "sqrt(A_planar)" : "cbrt(volume)")
                  << std::endl;
        if (roAuto || pAuto || aAuto) {
            std::cout << "[limiter] 警告: 基準値が自動決定なので、この run は**開始場に依存する作用素**である"
                      << " (分割実行が連続実行と一致しない)。固定するには solverConfig.yaml の space へ次をそのまま貼ること"
                      << " (plan convection-node-wall-reconstruction §4.25):" << std::endl;
            std::cout << "[limiter]   limiterRoRef: " << std::setprecision(10) << cfg.limiterRoRef << std::endl;
            std::cout << "[limiter]   limiterPRef: "  << cfg.limiterPRef  << std::endl;
            std::cout << "[limiter]   limiterARef: "  << cfg.limiterARef  << std::endl;
        }
    }

    MEMLOG_SOLVER("end of initializeSimulation");
    return cuda_cfg;
}

void writeStepOutputs(
    solverConfig& cfg,
    cudaConfig& cuda_cfg,
    mesh& msh,
    variables& var,
    point_probes& pprobes,
    int iStep)
{
    MEMLOG("before writeStepOutputs", memSolverSummary(msh, var, nullptr));
    outputH5_XDMF(cfg , msh, var, iStep);
    MEMLOG("after outputH5_XDMF", memSolverSummary(msh, var, nullptr));
    outputBconds_H5_XDMF(cfg , msh, var, iStep);
    MEMLOG("after outputBconds_H5_XDMF", memSolverSummary(msh, var, nullptr));
    pprobes.outputProbes(cfg , cuda_cfg , msh , var , iStep);
    // ソルバ内 CHT の壁温を残す (再開時は wall_profile_<physID>.csv にコピーして使う)。
    conjugateWall::writeConjugateState(cfg , msh , iStep);
}

void writeInitialOutputs(
    solverConfig& cfg,
    mesh& msh,
    variables& var)
{
    outputH5_XDMF(cfg , msh, var, 0);
}

// 1 ステップ分の状態を束ねる軽量コンテキスト。旧 advanceOneStep の [&] capture を置換し、
// 自由関数間でデータフローを明示する。
struct StepContext {
    solverConfig&       cfg;
    cudaConfig&         cuda_cfg;
    mesh&               msh;
    matrix&             mat_ns;   // 陰解法では未使用だが非陰解法のシグネチャに必要なため保持
    variables&          var;
    fluct_variables&    fluct;
    point_probes&       pprobes;
    RuntimeProfiler&    profiler;
    ResidualCsvLogger&  residual_logger;
    ImplicitDiagLogger& implicit_diag_logger;
    int                 iStep;
};

// 残差組み立ての単一情報源（旧 assembleCurrentState）。保存量から派生量・境界・勾配・各フラックス・
// ソース項を計算し res_* を確定する。explicit / implicit 双方が呼ぶ。
// 組立の前半 (状態射影・EOS・物性・BC・勾配まで) と後半 (リミッタ以降の流束・ソース) に分ける。
// 通常は assembleResidual が両方を続けて呼ぶだけで、挙動は分割前と同一。後半だけを複数回呼ぶのは
// 診断 (limiter-inlet-column-oscillation §5.1 #5e) が共通入力から枝分かれして評価するため。
static void assembleResidualPre(StepContext& s);
static void assembleResidualPost(StepContext& s);
void assembleResidual(StepContext& s, int stage_index)
{
    (void)stage_index;  // 現状カーネルは stage_index を使わない（dual-time 拡張用に interface 保持）
    assembleResidualPre(s);
    assembleResidualPost(s);
}

static void assembleResidualPre(StepContext& s)
{
    ledgerBeginAssemble(s.msh);   // 診断 (FORGE_DUMP_LEDGER、既定 off・出力専用)
    ledgerCapture(s.msh , s.var , "entry" , false);
    s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
        updateVariablesInner(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
    });
    // node-centered 壁 Dirichlet (SU2 SetVelocity_Old 相当): 更新後の壁ノード保存量を毎ステージ u=0 へ
    // 射影 (roe から KE 除去 + ρu=0)。運動量残差ゼロ (zeroWallMomentumResidual) だけでは ρ 変化で
    // u=ρu/ρ がドリフトするため状態再設定が必須。マルチマーカー emit (コーナー出口流出) と併用。
    // nodeWallDirichlet=0 / cell / 非 node では no-op。
    enforceWallNoSlip_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    // node × 軸対称: 軸ノード u_r=0 の状態ピン (roUy=0, roe から半径 KE 除去)。壁 no-slip と同相で、
    // 残差 0 化 (zeroAxisRadialResidual) と block-DPLUR の roUy 行 decouple (axis_ur_flag) と三点セット。
    enforceAxisSymmetry_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    // node 等温壁の温度ピンは**状態更新の直後・EOS/物性/RANS 壁境界より前**に置く (2026-09-20 修正)。
    // 以前は applyBconds の後に置いていたため、EOS は「更新で巻き戻った roe」と「更新後の ρ」から
    // 温度を作り、gasProperties と壁 ω がその誤った温度を使っていた。壁で ρ が動くケース
    // (翼列の前縁よどみ点) では EOS が ρ を床に張り付かせ、T が 9.2e6 K、μ が 151 倍、
    // 壁 ω が 8e15 s⁻¹ になって発散した (case/53 で実測。codex 診断 2026-09-20)。
    applyNodeIsothermalWallPin(s.cfg , s.cuda_cfg , s.msh , s.var);
    s.profiler.measureCuda(ProfileSection::DependentVariables, [&]() {
        speciesPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // Y_s = ρY_s/ρ (混合則 thermo の前)
        condensationPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // φ = ρφ/ρ (スカラ移流の上流値)
        tracerPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // ξ = ρξ/ρ (スカラ移流の上流値)
    });
    s.profiler.measureWall(ProfileSection::DependentVariables, [&]() {
        dependentVariables(s.cfg , s.cuda_cfg , s.msh , s.var, s.mat_ns);
        transitionPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // γ, Re_θt (初回初期化は Ux/k が要るのでここ)
    });
    s.profiler.measureCuda(ProfileSection::GasProperties, [&]() {
        gasProperties_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    tpuSnap("pre_state");   // 診断 G3-b: 状態の射影・クランプの後、境界の上書きの前 (記録中のみ; 既定は何もしない)
    s.profiler.measureWall(ProfileSection::ApplyBconds, [&]() {
        applyBconds(s.cfg , s.cuda_cfg , s.msh , s.var, s.mat_ns , s.fluct);
    });
    s.profiler.measureWall(ProfileSection::ApplyBconds, [&]() {
        applyRansScalarBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);
        applyTransitionBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);   // 遷移モデル: node 入口ピン (k のピンの後)
        applyWmlesWallModel(s.cfg , s.cuda_cfg , s.msh , s.var);   // WMLES 壁応力モデル (§10)
        applySstThermalWallFunction(s.cfg , s.cuda_cfg , s.msh , s.var);  // SST 熱的壁関数: 断熱壁 T_aw (§6.5(f))
        applySpeciesBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);
        applyCondensationBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);
        applyTracerBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    ledgerCapture(s.msh , s.var , "after_eos_bc" , false);   // 診断 (既定 no-op)
    ledgerCapture(s.msh , s.var , "res_before_conv" , true);
    s.profiler.measureCuda(ProfileSection::CalcGradient, [&]() {
        calcGradient_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        // 多成分 face 整合再構成 (speciesFaceReconstruction==1): ∇Y_s を Green-Gauss で計算。
        // species ghost は直前の applySpeciesBoundaries で Neumann 充填済み。既定 0 で no-op。
        // 旧: node + 粘性多成分でも ∇Y を計算していた (境界半割面の ghostless 弱形式 J_s=ρD∇Y·S 用) が、
        // 現在の species_diffusion_d は node 境界半割面を skip する (plan diffusion-node-boundary-real-distance §3(c))
        // ので ∇Y の読者は面整合再構成 (speciesFaceReconstruction≥1: SLAU Yd_recon / limiter_Y) だけ。
        // 3D 2.37 M 節点で毎ステップ 1.2 ms の無駄だった (plan performance-3d-node-sst-speedup)。
        if (s.cfg.speciesFaceReconstruction >= 1) {
            speciesGradient_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            passiveGradient_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // 受動種 ∇φ (passiveScalarScheme 1 のみ)
        }
        // node 周期境界 DOF 同一視 (§4.5 拡張): boundary periodic node の Green-Gauss 勾配を「和→broadcast」で
        // 厳密合併に直す (calcGradient_b_d で periodic 半割面は除外済み)。2次再構成・粘性の seam 精度向上。
        preGatherDumpMain(s.cfg , s.msh , s.var , "loop1");   // 診断 (FORGE_DUMP_PREGATHER、初回のみ・出力専用)
        periodicGradientGather_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::AxisymmetricSource, [&]() {
        axisymmetricGeomTerms_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
}

static void assembleResidualPost(StepContext& s)
{
    s.profiler.measureCuda(ProfileSection::Limiter, [&]() {
        limiter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::DucrosSensor, [&]() {
        ducrosSensor_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        turbulent_viscosity_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::ConvectiveFlux, [&]() {
        convectiveFlux_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var, s.mat_ns);
    });
    captureNodeIsothermalEnergyResidual(s.cfg , s.cuda_cfg , s.msh , s.var , "ifaceRconv");   // 診断 (既定 no-op)
    ledgerFlushFaces();                                            // 診断 (FORGE_DUMP_LEDGER、既定 no-op)
    ledgerCapture(s.msh , s.var , "res_after_conv" , true);
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        // k/ω 勾配と F1 を拡散の**前**に評価する (2026-09-08, plan turbulence-sst-consistency-options §2.1):
        // 旧順序 (transport → gradient → source) では拡散の非直交補正と σ ブレンドの F1 が前回評価の値
        // (Picard ラグ) だった。ransSource の F1 は同式・同入力なので sstF1 と一致する。
        ransGradient_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        ransBlendF1_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        ransTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);
        transitionTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);   // γ / Re_θt の移流拡散 (none で no-op)
    });
    ledgerCapture(s.msh , s.var , "res_after_rans_transport" , true);   // 診断 (既定 no-op)
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        speciesTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);  // 化学種移流残差
        chemistrySource_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);   // 有限速度化学ソース (ω_s, Q̇, 対角 Jacobian)
        tpoSnap("chem");      // 診断 G3-a: 残差の写し (記録中のみ; 既定は何もしない)
        speciesPinResidual_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var); // node 入口ピンノードの化学種残差除外 (cell は no-op)
        tpoSnap("sp_pin1");   // 診断 G3-a (記録中のみ)
    });
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        condensationTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);  // 液相モーメント移流残差 (Phase 1)
        condensationSource_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);     // 核生成+成長ソース (Phase 2)
        tpoSnap("cond_src");  // 診断 G3-a (記録中のみ)
        condThetaScan_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var, 0);       // θ_src の全評価を覆う集計 (#1b-pre; 計上のみ)
        tracerTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);        // 受動トレーサ移流残差 (node 入口ピン込み)
        tpoSnap("tracer");    // 診断 G3-a (記録中のみ; この 2 つは記録する成分に触れない)
        passivePinResidual_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);     // 受動種経路: node 入口ピンノードの残差除外 (ソース集計の後)
        tpoSnap("pas_pin");   // 診断 G3-a (記録中のみ)
        // 二相拡散 (#4e) は化学種の残差へ condensationTransport の中で足すので、化学種のピン除去をもう一度掛ける (周期集約の前)。
        if (condTwoPhaseDiffusionActive(s.cfg)) speciesPinResidual_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);
        tpoSnap("sp_pin2");   // 診断 G3-a (記録中のみ; 二相 OFF では pas_pin と同じ)
    });
    ledgerCapture(s.msh , s.var , "res_after_species" , true);   // 診断 (既定 no-op)
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        transitionSource_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // γ_eff を先に確定 (SST の k 式が同じ反復の値を読む)
        ransSource_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // k/ω 勾配は上 (ransTransport の前) で評価済み
    });
    s.profiler.measureCuda(ProfileSection::AxisymmetricSource, [&]() {
        axisymmetricSource_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);      // method 0 (r 重み): hoop 源
        axisymmetricSourceSU2_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // method 1 (SU2 流): 1/y 全ソース
        bodyForce_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // 一様体積力 (bodyForce, off なら no-op)
    });
    captureNodeIsothermalEnergyResidual(s.cfg , s.cuda_cfg , s.msh , s.var , "ifaceRpre");    // 診断 (既定 no-op)
    ledgerCapture(s.msh , s.var , "res_after_sources" , true);   // 診断 (既定 no-op)
    s.profiler.measureCuda(ProfileSection::ViscousFlux, [&]() {
        viscousFlux_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var, s.mat_ns);
    });
    // SU2 流の軸対称対称面 (MARKER_SYM): 軸上 CV の半径方向運動量残差を 0 に射影し roUy=0 を保つ。
    // explicit では ΔroUy=0 になり直接効く (implicit/block-DPLUR では連成 solve が補正を漏らすため Jacobian
    // 整合が別途必要・open issue, docs §7.1)。cell/非軸対称/平面では no-op。
    ledgerCapture(s.msh , s.var , "res_after_viscous" , true);   // 診断 (既定 no-op)
    zeroAxisRadialResidual_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    // node-centered 壁 Dirichlet: 壁ノードの運動量残差を 0 に射影し u=0 を保つ (壁ゴースト撤廃の代替)。
    // 軸射影の後に置き、コーナー (壁∩軸はまれだが) でも壁 no-slip を最終確定する。cell/非 node では no-op。
    captureNodeIsothermalEnergyResidual(s.cfg , s.cuda_cfg , s.msh , s.var , "ifaceRro" , "res_ro");   // 診断 (既定 no-op)
    zeroWallDirichletResiduals_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    // node-centered 周期境界 DOF 同一視 (median-dual M4, §4.5): 全 flux/source 積算 + 壁/軸射影の後に、
    // 周期 group の保存量残差を全員で足し合わせ全員へ書き戻す。合併体積と合わせ両側部分 CV を 1 CV として
    // 同期更新する (継ぎ目に双対面を作らず、両側内部双対面が res を組む)。cell/非周期では no-op。
    periodicNodeGather_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    ledgerCapture(s.msh , s.var , "res_final" , true);   // 診断 (既定 no-op): 壁射影・周期合併の後
    // 二相拡散 (#4e) の蒸気残差 res_roYv = res_roY_w − res_rog_0 (監視; 周期集約の後の確定残差から)。無効構成は no-op。
    twoPhaseVaporResidual_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    tpoSnap("final");   // 診断 G3-a: 組立後の確定残差 (記録中のみ; 既定は何もしない)
    // TODO(dual-time): unsteady のとき addUnsteadyTimeTerm(s) で BDF 物理時間項を res_* と
    // 対角に加える。定常では no-op。本体は次フェーズ。
}

// 残差スナップショットを CSV へ記録（inner_index==0 で outer_begin、それ以外で inner_iter）。
void logResidualSnapshot(StepContext& s, int inner_index)
{
    if (s.cfg.gpu == 1) {
        // GPU: device へ async reduce し行を buffer する (同期しない)。flush は logOuterEnd で行う。
        s.residual_logger.recordGpu(s.iStep, inner_index, s.msh, s.var);
    } else {
        const ResidualSnapshot residual_snapshot = gatherResidualSnapshot(s.cfg, s.msh, s.var);
        if (inner_index == 0) {
            s.residual_logger.logOuterBegin(s.iStep, inner_index, residual_snapshot);
        } else {
            s.residual_logger.logInnerIter(s.iStep, inner_index, residual_snapshot);
        }
    }
    if (s.implicit_diag_logger.enabled() && s.cfg.timeIntegration == 11) {
        s.implicit_diag_logger.log(s.iStep, inner_index, gatherImplicitDiagSnapshot(s.cfg, s.msh, s.var));
    }
}

// 古典 DPLUR 線形ソルバ。固定残差 res_* に対し Q を更新せず dq_block を nStepInner 回 Jacobi 緩和する。
// 各 sweep 後にバッファを swap し、最終補正は dq_block_old に残る（commit は呼び出し側）。
void blockDPLURSolve(StepContext& s, int subiter = 0)
{
    // 古典 DPLUR は dq=0 から開始する。前ステップの残留値による近傍参照を避けるため明示ゼロ化。
    // blockDPLUR==1: 5×5 block 版 (dq_block_*)、blockDPLUR==0: scalar 対角版 (dq_ro_* 等)。
    const size_t bytes = static_cast<size_t>(s.msh.nCells_all) * sizeof(flow_float);
    const bool useBlock = (s.cfg.blockDPLUR == 1);
    if (useBlock) {
        cudaMemset(s.var.c_d["dq_block_old_0"], 0, bytes);
        cudaMemset(s.var.c_d["dq_block_old_1"], 0, bytes);
        cudaMemset(s.var.c_d["dq_block_old_2"], 0, bytes);
        cudaMemset(s.var.c_d["dq_block_old_3"], 0, bytes);
        cudaMemset(s.var.c_d["dq_block_old_4"], 0, bytes);
    } else {
        cudaMemset(s.var.c_d["dq_ro_old"],   0, bytes);
        cudaMemset(s.var.c_d["dq_roUx_old"], 0, bytes);
        cudaMemset(s.var.c_d["dq_roUy_old"], 0, bytes);
        cudaMemset(s.var.c_d["dq_roUz_old"], 0, bytes);
        cudaMemset(s.var.c_d["dq_roe_old"],  0, bytes);
    }

    // line-implicit v2: lineKFreeze==1 の dual-time では K/diag 抽出と LU 分解をサブ反復 0 に
    // 限定し、以後のサブ反復は保存因子での solve だけにする (LHS 凍結 = defect-correction の
    // 近似強化。収束経路のみ変わり収束解は不変)。定常経路は subiter=0 固定で従来どおり毎回構築。
    const int lineStoreK =
        (s.cfg.lineKFreeze == 1) ? ((subiter == 0) ? 1 : 0) : 1;
    const int nSweep = std::max(1, s.cfg.nStepInner);
    for (int iSweep = 0; iSweep < nSweep; ++iSweep) {
        s.profiler.measureCuda(ProfileSection::TimeIntegration, [&]() {
            timeIntegration_d_wrapper(iSweep, s.cfg , s.cuda_cfg , s.msh , s.var, lineStoreK);
        });
        // line-implicit: ライン CV の dq_new を block-Thomas で上書き (swap 前)
        if (useBlock && s.cfg.lineImplicit == 1) {
            s.profiler.measureCuda(ProfileSection::TimeIntegration, [&]() {
                if (iSweep == 0 && lineStoreK == 1) {
                    lineThomasFactor_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
                }
                lineThomas_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            });
        }
        if (useBlock) {
            swapBlockImplicitCorrectionBuffers(s.var);
        } else {
            swapScalarImplicitCorrectionBuffers(s.var);
        }
        // node-centered 周期境界 DOF 同一視 (§4.5.7): swap 後の最新補正 dq_*_old を root から member へミラー。
        // 周期同一視ノードが master/slave で別 dq になり drift→発散するのを防ぐ。次 sweep の隣接 dq 読みと
        // 最終 commit を整合させる。cell/非周期では no-op。
        periodicMirrorDq_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    }
}

// 1 回の非線形（擬似時間）更新。定常・dual-time 共有の核。
// 診断用 freeze フラグ (環境変数ゲート, 既定 off, 構造体非変更)。CFL 律速の切り分け専用。
//   FORGE_FREEZE_SPECIES=1 : 化学種 ρY_s を更新しない (組成を凍結, flow だけ TP で解く → H1 切り分け)
//   FORGE_FREEZE_TURB=1    : SST k-ω を更新しない (μ_t は凍結 k,ω から再計算され実質固定 → H3 切り分け)
static bool freezeSpeciesEnabled() {
    static const bool v = [](){ const char* e = getenv("FORGE_FREEZE_SPECIES"); return e && atoi(e) != 0; }();
    return v;
}
static bool freezeTurbEnabled() {
    static const bool v = [](){ const char* e = getenv("FORGE_FREEZE_TURB"); return e && atoi(e) != 0; }();
    return v;
}


// 診断 (FORGE_PIN_DIAG=1, plan species-passive-scalar-unification §6-6-iii): dual-time の各サブ反復で
//   (a) BDF 追加 + ピン再除去の直後に、入口ピンノード (scalarDirichletPin==1) の化学種/受動種残差の max|res| (厳密 0 を期待)
//   (b) 更新 + 入口 Dirichlet 再適用の直後に、入口 bcond の境界ノード値と bvar 入口値の差 max|Y−Y_in| / |Xi−Xi_in| / |φ_moment|
// を host 読み戻しで印字する。既定 off (性能・出力に影響なし)。
static bool pinDiagEnabled() {
    static const bool v = [](){ const char* e = getenv("FORGE_PIN_DIAG"); return e && atoi(e) != 0; }();
    return v;
}
static void pinRowDiagnosticResidual(StepContext& s, int m)
{
    if (!pinDiagEnabled() || s.cfg.discretization != "node") return;
    const std::list<std::string> names = pinDiagResidualNames(s.var);   // 確保側 (hostCellSet) と共用
    s.var.copyVariables_cell_D2H(names);
    const auto& pin = s.var.hostCell("scalarDirichletPin");
    geom_int nPin = 0; for (geom_int ic = 0; ic < s.msh.nCells; ++ic) if (pin[ic] == (flow_float)1.0) ++nPin;
    std::ostringstream os; os << "[pin-diag] step " << s.iStep + 1 << " subiter " << m << " pinned nodes " << nPin << " max|res| on pinned:";
    for (const auto& nm : names) {
        if (nm == "scalarDirichletPin") continue;
        const auto& r = s.var.hostCell(nm); double mx = 0.0;
        for (geom_int ic = 0; ic < s.msh.nCells; ++ic) if (pin[ic] == (flow_float)1.0) mx = std::max(mx, (double)std::abs(r[ic]));
        os << " " << nm << " " << std::scientific << std::setprecision(2) << mx;
    }
    std::cout << os.str() << "\n";
}
static void pinRowDiagnosticState(StepContext& s, int m)
{
    if (!pinDiagEnabled() || s.cfg.discretization != "node") return;
    const std::list<std::string> names = pinDiagStateNames(s.var);   // 確保側 (hostCellSet) と共用
    if (names.empty()) return;
    s.var.copyVariables_cell_D2H(names);
    std::ostringstream os; os << "[pin-diag] step " << s.iStep + 1 << " subiter " << m << " inlet-node values vs bvar:";
    for (auto& bc : s.msh.bconds) {
        if (bc.bcondKind.rfind("inlet_", 0) != 0 || bc.iPlanes.empty()) continue;
        const geom_int nb = (geom_int)bc.iPlanes.size();
        std::vector<geom_int> cell(nb);
        gpuErrchk( cudaMemcpy(cell.data(), bc.map_bplane_cell_d, nb*sizeof(geom_int), cudaMemcpyDeviceToHost) );
        for (const auto& nm : names) {
            const auto& v = s.var.hostCell(nm);
            std::vector<flow_float> bv;
            const bool isMoment = std::find(s.var.condMomentConsNames.begin(), s.var.condMomentConsNames.end(), "ro" + nm) != s.var.condMomentConsNames.end();
            if (!isMoment) {
                auto it = bc.bvar_d.find(nm);
                if (it == bc.bvar_d.end() || it->second == nullptr) continue;
                bv.resize(nb); gpuErrchk( cudaMemcpy(bv.data(), it->second, nb*sizeof(flow_float), cudaMemcpyDeviceToHost) );
            }
            double mx = 0.0;
            for (geom_int ib = 0; ib < nb; ++ib) {
                const double ref = isMoment ? 0.0 : (double)std::max(bv[ib], (flow_float)0.0);
                mx = std::max(mx, std::abs((double)v[cell[ib]] - ref));
            }
            os << " " << bc.bcondKind << "/" << nm << " " << std::scientific << std::setprecision(2) << mx;
        }
    }
    std::cout << os.str() << "\n";
}

// ---- limiter-inlet-column-oscillation §5.1 #5e (診断専用・既定 off) ------------------------------------------------
// FORGE_DIAG_PSI_DUALEVAL="N0,N1": 外反復 N0 (0 始まりの iStep) の通常組立の後に流れ 5 変数の ψ を保存し、
// N0 < iStep <= N1 の各外反復で、組立の前半 (状態射影・EOS・BC・勾配) を 1 回だけ行ってから、
// 全 device 配列 (var.c_d・var.p_d・各境界の bvar_d) を退避し、後半 (リミッタ以降) を 3 回評価する:
//   B  = 入口集合 (ccx < FORGE_DIAG_PSI_XMAX [m]、省略時は全節点) だけ保存 ψ に差し替え
//   A' = 差し替えなし (退避から復元した同じ入力で。再評価の誤差の基準)
//   A  = 差し替えなし (退避から復元した同じ入力で。**これを時間更新に使う**)
// 各評価の前に退避した配列を書き戻すので 3 枝の入力は同一 (codex diagnose 2026-10-03 dualeval-result の指定)。
// 退避の外にある状態 (static な作業配列など) は検査できないので、A'−A を再評価誤差として記録する。
struct PsiDualEval {
    bool on = false; int n0 = -1, n1 = -1; flow_float xmax = (flow_float)1e30; bool saved = false;
    std::vector<char> inS; std::ofstream csv;
};
static PsiDualEval g_pde;

static void pdeInit(StepContext& s)
{
    static bool done = false; if (done) return; done = true;
    const char* e = std::getenv("FORGE_DIAG_PSI_DUALEVAL"); if (!e) return;
    if (std::sscanf(e, "%d,%d", &g_pde.n0, &g_pde.n1) != 2 || g_pde.n0 < 0 || g_pde.n1 <= g_pde.n0) {
        std::cerr << "[psi-dualeval] FORGE_DIAG_PSI_DUALEVAL must be \"N0,N1\" with 0 <= N0 < N1\n"; exit(EXIT_FAILURE);
    }
    if (const char* x = std::getenv("FORGE_DIAG_PSI_XMAX")) g_pde.xmax = (flow_float)std::atof(x);
    g_pde.on = true;
    std::vector<flow_float> cx(s.msh.nCells_all);
    gpuErrchk( cudaMemcpy(cx.data(), s.var.c_d["ccx"], cx.size()*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    g_pde.inS.assign(s.msh.nCells_all, 0); long nS = 0;
    for (geom_int i = 0; i < s.msh.nCells; ++i) if (cx[i] < g_pde.xmax) { g_pde.inS[i] = 1; ++nS; }
    g_pde.csv.open("psi_dualeval.csv");
    g_pde.csv << "step,var,nS,S_A,S_B,S_BmA,S_ApmA,S_AdotBmA,all_A,all_B,all_BmA,all_ApmA,all_AdotBmA,psi_maxdiff_S,restore_mismatch\n";
    printf("[psi-dualeval] ON (branching): save psi at iStep %d, dual evaluation for iStep %d..%d, xmax %.6g m, %ld nodes in S\n",
           g_pde.n0, g_pde.n0 + 1, g_pde.n1, (double)g_pde.xmax, nS);
}

// 退避: 名前 → (device ポインタ, 要素数, host 複製)
struct PdeSnap { std::vector<std::tuple<flow_float*, size_t, std::vector<flow_float>>> arr; };
static size_t pdeSize(const std::map<std::string, std::vector<flow_float>>& host, const std::string& k, size_t fallback)
{
    auto it = host.find(k); return (it != host.end() && !it->second.empty()) ? it->second.size() : fallback;
}
static PdeSnap pdeSnapshot(StepContext& s)
{
    PdeSnap sn;
    auto add = [&](flow_float* p, size_t n) { if (!p || n == 0) return; std::vector<flow_float> h(n);
        gpuErrchk( cudaMemcpy(h.data(), p, n*sizeof(flow_float), cudaMemcpyDeviceToHost) ); sn.arr.emplace_back(p, n, std::move(h)); };
    for (auto& kv : s.var.c_d) add(kv.second, pdeSize(s.var.c, kv.first, s.msh.nCells_all));
    // 面配列: gpu: 1 ではホスト側を確保しない (長さ 0) ので、fallback は allocVariables の cudaMalloc と同じ nPlanes
    // (0 にすると面配列が黙って退避から外れる。plan architecture-solver-host-memory §4.2)
    for (auto& kv : s.var.p_d) add(kv.second, pdeSize(s.var.p, kv.first, s.msh.nPlanes));
    for (auto& bc : s.msh.bconds) for (auto& kv : bc.bvar_d) add(kv.second, pdeSize(bc.bvar, kv.first, 0));
    static bool logged = false;
    if (!logged) { logged = true; size_t nb = 0; for (auto& a : sn.arr) nb += std::get<1>(a);
        size_t nTot = s.var.c_d.size() + s.var.p_d.size(); for (auto& bc : s.msh.bconds) nTot += bc.bvar_d.size();
        printf("[psi-dualeval] snapshot: %zu of %zu device arrays (%.1f MB)\n", sn.arr.size(), nTot, nb*sizeof(flow_float)/1048576.0); }
    return sn;
}
static void pdeRestore(const PdeSnap& sn)
{
    for (auto& a : sn.arr) gpuErrchk( cudaMemcpy(std::get<0>(a), std::get<2>(a).data(), std::get<1>(a)*sizeof(flow_float), cudaMemcpyHostToDevice) );
}
static long pdeRestoreMismatch(const PdeSnap& sn)   // 書き戻しが効いたかの検査 (配列数)
{
    long n = 0;
    for (auto& a : sn.arr) { std::vector<flow_float> h(std::get<1>(a));
        gpuErrchk( cudaMemcpy(h.data(), std::get<0>(a), h.size()*sizeof(flow_float), cudaMemcpyDeviceToHost) );
        n += (std::memcmp(h.data(), std::get<2>(a).data(), h.size()*sizeof(flow_float)) != 0); }
    return n;
}
static std::vector<std::string> pdeResKeys(StepContext& s)
{
    std::vector<std::string> k;
    for (auto& kv : s.var.c_d) if (kv.second && kv.first.rfind("res_", 0) == 0) k.push_back(kv.first);
    return k;
}
static std::map<std::string, std::vector<flow_float>> pdeCopy(StepContext& s, const std::vector<std::string>& keys)
{
    std::map<std::string, std::vector<flow_float>> m;
    for (auto& k : keys) { auto& v = m[k]; v.resize(s.msh.nCells_all);
        gpuErrchk( cudaMemcpy(v.data(), s.var.c_d[k], v.size()*sizeof(flow_float), cudaMemcpyDeviceToHost) ); }
    return m;
}

// 窓内なら枝分かれ評価を行って true を返す (このとき通常の assembleResidual は呼ばない)。
static bool pdeAssemble(StepContext& s)
{
    pdeInit(s);
    if (!g_pde.on || !g_pde.saved || s.iStep <= g_pde.n0 || s.iStep > g_pde.n1) return false;
    const auto resKeys = pdeResKeys(s);
    assembleResidualPre(s);
    const PdeSnap sn = pdeSnapshot(s);
    limiterPsiOverride(true);
    assembleResidualPost(s);                                  // B
    limiterPsiOverride(false);
    const auto B = pdeCopy(s, resKeys);
    double psiMax = 0.0;
    {
        const char* ln[5] = {"limiter_ro","limiter_Ux","limiter_Uy","limiter_Uz","limiter_P"};
        for (int k = 0; k < 5; ++k) {
            std::vector<flow_float> a(s.msh.nCells_all), b(s.msh.nCells_all);
            gpuErrchk( cudaMemcpy(a.data(), s.var.c_d[ln[k]], a.size()*sizeof(flow_float), cudaMemcpyDeviceToHost) );
            gpuErrchk( cudaMemcpy(b.data(), limiterPsiSaved(k), b.size()*sizeof(flow_float), cudaMemcpyDeviceToHost) );
            for (geom_int i = 0; i < s.msh.nCells; ++i) if (g_pde.inS[i]) psiMax = std::max(psiMax, (double)std::fabs(a[i]-b[i]));
        }
    }
    pdeRestore(sn); long mm = pdeRestoreMismatch(sn);
    assembleResidualPost(s);                                  // A'
    const auto Ap = pdeCopy(s, resKeys);
    pdeRestore(sn); mm += pdeRestoreMismatch(sn);
    assembleResidualPost(s);                                  // A (時間更新に使う)
    const auto A = pdeCopy(s, resKeys);
    long nS = 0; for (geom_int i = 0; i < s.msh.nCells; ++i) nS += g_pde.inS[i];
    for (auto& k : resKeys) {
        const auto& a = A.at(k); const auto& b = B.at(k); const auto& ap = Ap.at(k);
        double sA=0, sB=0, sD=0, sP=0, sX=0, gA=0, gB=0, gD=0, gP=0, gX=0;
        for (geom_int i = 0; i < s.msh.nCells; ++i) {
            const double va = a[i], vb = b[i], d = vb - va, dp = (double)ap[i] - va;
            gA += va*va; gB += vb*vb; gD += d*d; gP += dp*dp; gX += va*d;
            if (g_pde.inS[i]) { sA += va*va; sB += vb*vb; sD += d*d; sP += dp*dp; sX += va*d; }
        }
        g_pde.csv << s.iStep << "," << k << "," << nS << "," << std::sqrt(sA) << "," << std::sqrt(sB) << "," << std::sqrt(sD) << ","
                  << std::sqrt(sP) << "," << sX << "," << std::sqrt(gA) << "," << std::sqrt(gB) << "," << std::sqrt(gD) << ","
                  << std::sqrt(gP) << "," << gX << "," << psiMax << "," << mm << "\n";
    }
    g_pde.csv.flush();
    return true;
}
static void pdeAfterNormal(StepContext& s)
{
    if (!g_pde.on) return;
    if (s.iStep == g_pde.n0 && !g_pde.saved) {               // 通常組立の ψ を保存
        limiterPsiSave(s.cfg, s.msh, s.var, g_pde.xmax); g_pde.saved = true;
        printf("[psi-dualeval] psi saved at iStep %d\n", s.iStep);
    }
}

// 診断 G3-b (FORGE_DIAG_TP_UPDATE): runTpUpdateDiag が組立を自分で (前処理・後処理を分けて) 通したあと、更新だけを 1 回行うための印。
// 既定 false (通常の計算では立たない)。立っていたら 1 回だけ組立を飛ばして下ろす。
static bool g_tpuSkipAssemblyOnce = false;

// 残差 1 回構築 → 局所擬似時間 dτ → 古典 DPLUR 線形解 → Q への commit。
void implicitNonlinearUpdate(StepContext& s, int inner_index)
{
    // limiter-inlet-column-oscillation §5.1 #5e 診断 (既定 off): 窓内は共通入力から枝分かれして組み、A を残す
    if (g_tpuSkipAssemblyOnce) g_tpuSkipAssemblyOnce = false;   // 診断 G3-b: 組立は呼び出し側で済んでいる
    else if (!pdeAssemble(s)) { assembleResidual(s, 1); pdeAfterNormal(s); }
    logResidualSnapshot(s, inner_index);
    // #1b-r2 診断 (condTwoPhaseDiag 3, 窓内だけ): 同じ状態・面値・係数・ソース値の double 組立 B (状態・組立 A は不変; 読むだけ)
    if (s.cfg.condTwoPhaseDiag == 3 && twoPhaseDiagInWindow(s.cfg)) {
        // 試験用 FORGE_TPD3_VERIFY=1: 組立 B の前後で全セル配列 (ghost 込み) をバイト比較し、B が状態を書かないことを確かめる
        const bool verify = (std::getenv("FORGE_TPD3_VERIFY") != nullptr);
        std::map<std::string, std::vector<flow_float>> before;
        if (verify) for (auto& kv : s.var.c_d) if (kv.second) { auto& v = before[kv.first]; v.resize(s.msh.nCells_all);
            gpuErrchk( cudaMemcpy(v.data(), kv.second, v.size()*sizeof(flow_float), cudaMemcpyDeviceToHost) ); }
        twoPhaseDiagB_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        if (verify) {
            long nArr = 0, nDiff = 0; std::string first;
            for (auto& kv : before) { std::vector<flow_float> a(kv.second.size());
                gpuErrchk( cudaMemcpy(a.data(), s.var.c_d[kv.first], a.size()*sizeof(flow_float), cudaMemcpyDeviceToHost) ); ++nArr;
                if (std::memcmp(a.data(), kv.second.data(), a.size()*sizeof(flow_float)) != 0) { ++nDiff; if (first.empty()) first = kv.first; } }
            printf("[twophase-diag] verify step %d: %ld cell arrays compared byte-wise before/after assembly B, %ld changed%s%s\n",
                   s.iStep + 1, nArr, nDiff, nDiff ? " (first: " : "", nDiff ? (first + ")").c_str() : "");
        }
    }
    // 定常 (unsteady==0) implicit では dt_local=cfl_pseudo·dx/λ で cfg.dt が打ち消され、dt 適応も表示も
    // monitorInterval ごとで足りる (per-step host 同期を回避)。dt 適応と表示は同一 (monitor 時のみ host 読み)。
    // unsteady でここに来る経路は無い (implicit unsteady は dual-time) が、防御的に毎ステップ adapt にする。
    const bool onMonitor = (s.iStep % s.cfg.monitorInterval == 0);
    const bool adaptDt  = (s.cfg.unsteady != 0) || onMonitor;
    // max cfl の host 読みは unsteady のモニタ行にしか要らない (定常では cfg.dt が無意味なので読まない)。
    const bool printCfl = onMonitor && (s.cfg.unsteady == 1);
    s.profiler.measureCuda(ProfileSection::SetDt, [&]() {
        setDT_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var, adaptDt, printCfl);
    });
    const bool freezeSpecies = freezeSpeciesEnabled();
    const bool freezeTurb    = freezeTurbEnabled();
    // 案C (speciesImplicitCoupling==2): block-triangular roe↔roY coupling。
    // flow block 解の前に species 仮更新 δ(ρY)* を予測 → 接空間射影 z → 解析 EOS-JVP δp_Y を
    // res_roe へ移項し、flow が組成変化 (T,p,h への影響) を同一 block 解の中で見るようにする。
    // freezeSpecies 時は予測/移項もしない (組成完全凍結)。
    const bool eosCoupled = speciesEOSCoupled(s.cfg, s.var) && !freezeSpecies;
    if (eosCoupled) {
        s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
            speciesUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // roY_N = roY
            speciesEOSCrossPredictInject_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        });
    }

    passiveSaveRhoPre_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // 受動種の φ_N δρ 項用に更新前 ρ を退避 (#19)
    blockDPLURSolve(s);
    s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
        if (s.cfg.blockDPLUR == 1) {
            applyBlockImplicitCorrection(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
        } else {
            applyScalarImplicitCorrection(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
        }
    });
    // 注: 軸上 roUy=0 の commit 後強制は block-DPLUR と非整合で Mach~1000 に発散 (外部状態手術不可)。
    // SU2 流の対称面は Jacobian 内で対称化する必要がある (open issue, docs §7.1)。暫定で無効。
    // enforceAxisSymmetry_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    // SST (k-ω) を segregated point-implicit で更新（凍結解除）。残差・消散ヤコビアンは
    // 直前の assembleResidual (ransSource) で確定済み、dt_local は setDT 済み。
    if (scalarResidualEnabled(s.cfg) && !freezeTurb) {
        s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
            sstEnergyKCorrection_begin_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // sstEnergyIncludesK: roK 退避
            applySSTPointImplicit(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
            // node 周期 DOF 同一視 (§4.5): point-implicit SST 更新後に k/ω 状態を root→member ミラーし drift を防ぐ。
            periodicMirrorScalarState_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            applyTransitionPointImplicit_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // 遷移モデル γ / Re_θt (周期ミラー込み)
            sstEnergyKCorrection_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var, 0);   // E_t 保存: roe -= (roK − roK_prev) (増分更新)
        });
    }

    // 化学種 (多成分 TP) の陰解法更新（凍結解除）。res_roY/transport_diag/src_jac は assembleResidual の
    // speciesTransport で確定済み。baseline roY_N=roY を取り、実現可能性再正規化で閉じる。
    // 各 wrapper は nSpecies<2 で no-op (単成分/CPG は不変)。
    // speciesImplicitCoupling==2: 案C block-triangular。予測 δ(ρY)* は上で済んでいるので、ここでは
    //   flow 密度更新 δρ を含めた最終 commit ρY_s=ρY_s^N+z_s+Y_s^N δρ を行う。
    // speciesImplicitCoupling==1: 緩和整合 scalar-DPLUR (流れ block と同一 dt/implicitRelax/nStepInner sweep)。
    //                          =0: 従来 segregated 点陰的 forward-Euler (既定・ビット不変)。
    // freezeSpecies 時は化学種更新を完全にスキップ (ρY_s 凍結)。EOS は凍結 ρY/ρ で評価される。
    // 二相拡散 (condTwoPhaseDiffusion, #4e): 水 ρY_w は化学種の更新で commit せず、液・Q と一緒に非分割更新 (蒸気/液の増分) で commit し、
    // その後に再正規化 (係数を液・Q にも) する (plan condensation-two-phase-transport §4.2、設計メモ §6.1・§14)。化学種凍結時は現行経路。
    const bool twoPhase = condTwoPhaseDiffusionActive(s.cfg) && !freezeSpecies;
    if (!freezeSpecies) {
        s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
            if (eosCoupled) {
                speciesEOSFinalCommit_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            } else {
                speciesUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // roY_N = roY_M = roY
                if (speciesImplicitCoupled(s.cfg, s.var)) {
                    speciesImplicitDPLURSolve_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
                } else {
                    speciesTimeIntegration_d_wrapper(0, s.cfg , s.cuda_cfg , s.msh , s.var);
                }
            }
            tpuSnap("sp_commit");   // 診断 G3-b (記録中のみ; 既定は何もしない)
            if (twoPhase) {
                twoPhaseHoldWater_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // 水は更新前に戻す (液の後で commit)
                tpuSnap("hold");   // 診断 G3-b: 水の戻し (化学種の commit を取り消す格納差; 記録中のみ)
            } else {
                speciesRenormalize_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
                tpuSnap("sp_renorm");   // 診断 G3-b (記録中のみ)
                periodicMirrorSpeciesState_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // node 周期: 化学種状態を root→member (§4.1-5)
                speciesPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);     // Y=roY/ρ (出力/次残差用に同期)
            }
        });
    }

    // 液相モーメント (非平衡凝縮) を segregated point-implicit で更新 (Phase 1 はソース=0 の純移流)。
    // res_/transport_diag は assembleResidual の condensationTransport で確定済み。各 wrapper は
    // condensation==0 で no-op。
    s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
        condensationUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // ro*_N = ro*_M = ro*
        if (twoPhase) {
            // 蒸気・液・Q の非分割更新 → 再正規化 (係数を液・Q にも) → 受動種の砦・周期ミラー (wrapper 内)。化学種の周期ミラーと Y の同期はここ。
            twoPhaseUpdate_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            periodicMirrorSpeciesState_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            speciesPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        } else {
            condensationTimeIntegration_d_wrapper(0, s.cfg , s.cuda_cfg , s.msh , s.var);
        }
        // 更新の θ (二相の非分割更新 / 更新クランプ θ_u) の全更新を覆う集計 (#1b-pre; 計上のみ)。θ を書かない経路では呼ばない。
        if (twoPhase || (s.cfg.timeIntegration == 11 && s.cfg.condLimiterMode == 1 && s.cfg.condEquilibrium == 0))
            condThetaScan_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var, 1);
        condensationPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);     // φ=ρφ/ρ (出力/次残差用に同期)
        tpuSnap("realizability");   // 診断 G3-b: 実現可能性クランプの後 (記録中のみ; 既定は何もしない)
        if (twoPhase) twoPhaseCorrGateEnd();   // #4h: この更新の補正計測を閉じる (実現可能性クランプの後)
        // 受動トレーサ (segregated point-implicit)。tracer 無効で no-op。
        tracerUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        tracerTimeIntegration_d_wrapper(0, s.cfg , s.cuda_cfg , s.msh , s.var);
        tracerPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
}

// 陽解法 Runge-Kutta（tI 1/3/4）。steady/unsteady 両対応。
void advanceExplicitRK(StepContext& s)
{
    const int iteration_count = s.cfg.perStepIterationCount();
    const char* iteration_label = s.cfg.perStepIterationLabel();

    for (int iloop = 0 ; iloop < iteration_count ; iloop++) {
        s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
            updateVariablesInner(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
            speciesUpdateInner_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // roY{s}M ステージ始点
            condensationUpdateInner_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // 液相モーメント M ステージ始点
            tracerUpdateInner_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // トレーサ M ステージ始点
        });

        (void)iteration_label;   // 旧 "Stage : n" 行は廃止 (console モニタ行に集約)
        assembleResidual(s, iloop + 1);
        logResidualSnapshot(s, iloop);
        s.profiler.measureCuda(ProfileSection::TimeIntegration, [&]() {
            timeIntegration_d_wrapper(iloop, s.cfg , s.cuda_cfg , s.msh , s.var);
            // node 周期境界 DOF 同一視 (§4.5.9): NS 保存量更新直後に root→member ミラーで slave=master を強制。
            // 残差 gather だけでは初期 desync (非周期 seed 摂動) が残り継ぎ目フラックス不整合を生むため。cell/非周期で no-op。
            periodicMirrorNSState_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            ransTimeIntegration_d_wrapper(iloop, s.cfg , s.cuda_cfg , s.msh , s.var);
            sstEnergyKCorrection_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var, 1);   // E_t 保存: roe -= (roK − roKN) (RK stage は N から組み直す)
            speciesTimeIntegration_d_wrapper(iloop, s.cfg , s.cuda_cfg , s.msh , s.var);
            speciesRenormalize_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // ρY_s>=0, ΣρY_s=ρ
            periodicMirrorSpeciesState_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // node 周期: 化学種状態ミラー (§4.1-5)
            condensationTimeIntegration_d_wrapper(iloop, s.cfg , s.cuda_cfg , s.msh , s.var);  // 液相モーメント (Phase 1 ソース=0)
            tracerTimeIntegration_d_wrapper(iloop, s.cfg , s.cuda_cfg , s.msh , s.var);  // 受動トレーサ
        });
        // 注: explicit 軸対称は軸 CV が step1 で発散するため (recipe 併用でも不変)、enforce は呼んでも
        // 検証できない。explicit の near-axis 安定化は別途要 (open issue)。暫定で無効。
        // enforceAxisSymmetry_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    }

    s.profiler.measureWall(ProfileSection::UpdateOuter, [&]() {
        updateVariablesOuter(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
        speciesUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // roY{s}N/M 次ステップ用
        speciesPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);    // 出力 Y_s を最終 roY_s と同期
        condensationUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // 液相モーメント N/M 次ステップ用
        condensationPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);    // 出力 φ を最終 ρφ と同期
        tracerUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // トレーサ N/M 次ステップ用
        tracerPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);    // 出力 ξ を最終 ρξ と同期
    });
    s.profiler.measureWall(ProfileSection::WriteOutputs, [&]() {
        writeStepOutputs(s.cfg , s.cuda_cfg , s.msh , s.var , s.pprobes , s.iStep+1);
    });
    // explicit (陽解法) は cfg.dt を時間前進に使うため dt 適応は毎ステップ行う。max cfl の host 読みは
    // unsteady のモニタ行用に monitorInterval で間引く (定常局所 dt では cfg.dt は無意味なので読まない)。
    const bool printCflExp = (s.iStep % s.cfg.monitorInterval == 0) && (s.cfg.unsteady == 1);
    // 物理時間はこの step の前進に使った dt で進める。setDT (dtControl==1 の適応) の後に足すと次 step 用の dt が
    // 加算され t がずれる (旧実装のバグ。適応 sod で step 1 の t=1.087e-6 ≠ 使用 dt 1.097e-6 として露見, 2026-09-09)。
    if (s.cfg.unsteady == 1) {
        s.cfg.totalTimeD += (double)s.cfg.dt; s.cfg.totalTime = static_cast<flow_float>(s.cfg.totalTimeD);
    }
    s.profiler.measureCuda(ProfileSection::SetDt, [&]() {
        setDT_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var, /*adaptDt=*/true, /*printCfl=*/printCflExp);
    });
    s.residual_logger.logOuterEnd(s.iStep);
}

// 定常 block DPLUR 陰解法。メインループ（nStepOuter）を擬似時間とし、1 ステップ = 1 非線形更新の縮退形。
void advanceImplicitSteady(StepContext& s)
{
    // baseline (roN) は前ステップ末尾 / 初期化の updateVariablesOuter で設定済み（ro == roN）。
    implicitNonlinearUpdate(s, 0);

    // **最後の更新のあとにもピンを当てる**: dq_roe=0 なので更新は roe を step 冒頭の値へ戻す。
    // ここで当てないと、出力される保存量と次ステップの基準 (roN) が等温条件を満たさず、
    // ρ が動くほど roe/(ρ c_v T_w) がずれていく (実測: 壁ノードの roe が 1000 step ビット不変)。
    applyNodeIsothermalWallPin(s.cfg , s.cuda_cfg , s.msh , s.var);

    s.profiler.measureWall(ProfileSection::UpdateOuter, [&]() {
        updateVariablesOuter(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
    });
    // FP64 影アキュムレータの**計器**: reconcile が採用したセル数 (§4.4)。
    // 「上書き型の writer に書き換えられた保存量の数」であり、run ごとに**期待値を事前登録して
    // 突き合わせる** (case/56 なら等温壁ノード数、case/09 や case/44 なら 0)。
    // 期待を超える = 棚卸しできていない writer がいる、という意味なので計器として出す。
    if (s.cfg.qAccumulatorFP64 == 1 && s.var.qaccAdopt_d != nullptr
        && s.iStep % s.cfg.monitorInterval == 0) {
        const int nAdopt = qaccReadAdoptCounter(s.var.qaccAdopt_d);
        printf("[qAccumulatorFP64] step %d: reconcile 採用 %d (= 上書き型 writer が触った保存量の数)\n",
               s.iStep + 1, nAdopt);
        qaccResetAdoptCounter(s.var.qaccAdopt_d);
    }
    s.profiler.measureWall(ProfileSection::WriteOutputs, [&]() {
        writeStepOutputs(s.cfg , s.cuda_cfg , s.msh , s.var , s.pprobes , s.iStep+1);
    });
    // 末尾の setDT は撤去 (冗長): 定常では dt_local は次ステップ冒頭の implicitNonlinearUpdate→setDT で
    // 再計算され、ここで計算した cfl/dt_local は使われる前に上書きされる純粋な無駄 (~80µs/step の setCFL
    // カーネル×3)。cfg.dt も dt_local に効かず cosmetic。max cfl/dt 表示は冒頭 setDT が monitorInterval ごとに行う。
    s.residual_logger.logOuterEnd(s.iStep);
}

// 非定常 dual-time 陰解法。1 物理ステップ = 時間レベルシフト → 擬似時間サブ反復（BDF 物理時間項つき
// 非線形更新を nSubIterDualTime 回）→ 物理時間前進。BDF1（初回）/ BDF2（以降, bdfOrder==2 のとき）。
void advanceImplicitDualTime(StepContext& s)
{
    if (s.cfg.dualTime != 1) {
        throw std::runtime_error(
            "Unsteady implicit (timeIntegration==11 && unsteady==1) requires dualTime=1.");
    }
    if (s.cfg.blockDPLUR != 1) {
        throw std::runtime_error(
            "Dual-time implicit currently supports blockDPLUR=1 (5x5 block) only.");
    }
    if (s.cfg.dtControl != 0) {
        // 物理 Δt は固定でなければならない（setDT が cfg.dt を CFL 適応すると BDF 項が壊れる）。
        throw std::runtime_error(
            "Dual-time implicit requires time.deltaT.control=0 (fixed physical dt).");
    }

    // 物理時間レベルシフト: roNN ← roN, roN ← ro（現在の ro = Q^n）。化学種 roY_PP ← roY_P ← roY、受動種 (scheme 1) も同時に
    // (plan species-passive-scalar-unification §4.4: 履歴は流れ・化学種・受動種で 1 つの状態)。
    s.profiler.measureWall(ProfileSection::UpdateOuter, [&]() {
        shiftDualTimeLevels_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        speciesShiftDualTimeLevels_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        passiveShiftDualTimeLevels_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });

    // BDF 係数: 履歴が無い (fresh start / 旧形式 restart) 最初の物理 step または bdfOrder==1 は BDF1 (1,1,0)、以降 BDF2 (3/2,2,1/2)。
    // 履歴有効数 cfg.nHistoryValid は全系共有 (checkpoint から復元、物理 step 完了ごとに +1)。
    const bool useBDF2 = (s.cfg.bdfOrder >= 2) && (s.cfg.nHistoryValid >= 1);
    const flow_float a = useBDF2 ? static_cast<flow_float>(1.5) : static_cast<flow_float>(1.0);
    const flow_float b = useBDF2 ? static_cast<flow_float>(2.0) : static_cast<flow_float>(1.0);
    const flow_float c = useBDF2 ? static_cast<flow_float>(0.5) : static_cast<flow_float>(0.0);
    const int include_scalar = scalarResidualEnabled(s.cfg) ? 1 : 0;
    // 対角へ加える物理時間項係数 a/Δt（block/scalar/SST カーネルが cfg 経由で参照）。
    s.cfg.unsteadyDiagCoef = a / std::max(s.cfg.dt, static_cast<flow_float>(1.0e-30));

    // 化学種 (多成分): 2026-09-13 まで dual-time は化学種を一切更新していなかった (ρY_s が初期場のまま凍結、ρ だけ動いて
    // ΣY_s=ρ⁰/ρ になる; chem e296f0d0 で修正)。定常陰解法と同じ更新 (coupling 0 point-implicit / 1 scalar-DPLUR / 2 案C 予測→block→commit)
    // を各サブ反復で回し、残差には BDF 項を入れる。案C の commit の δρ 基準は予測時点の ρ (speciesEOSCrossPredictInject が保存)。
    const bool freezeSpecies = freezeSpeciesEnabled();
    const bool eosCoupled = speciesEOSCoupled(s.cfg, s.var) && !freezeSpecies;
    const bool haveSpecies = (s.var.nSpeciesRegistered > 1) && !freezeSpecies;

    const int nSub = std::max(1, s.cfg.nSubIterDualTime);
    for (int m = 0; m < nSub; ++m) {
        // (1) 空間残差 (+ 化学種/受動種のピン残差除去 + 周期 gather) → (2) BDF 項 (合併体積で一度だけ; 流れ・k/ω・化学種・受動種)
        //  → (3) ピン残差の再除去 (BDF がピン行に残差を戻すため)。
        assembleResidual(s, 1);
        addUnsteadyTimeTerm_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var, a, b, c, include_scalar);
        speciesAddUnsteadyTimeTerm_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var, a, b, c);
        passiveAddUnsteadyTimeTerm_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var, a, b, c);
        speciesPinResidual_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        passivePinResidual_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        pinRowDiagnosticResidual(s, m);   // FORGE_PIN_DIAG=1: ピン行の残差が 0 か
        logResidualSnapshot(s, m);   // BDF 込みのサブ反復残差 (inner_iter 行; 化学種・受動種列を含む)
        // dual-time は dtControl==0 を強制している (上の検査) ので adaptDt=true でも cfg.dt は変わらない
        // (host 読みは printCflDt のときだけ発生)。max cfl (物理 CFL) の格納は monitorInterval で間引く (モニタ行が表示)。
        const bool printCflDt = (s.iStep % s.cfg.monitorInterval == 0);
        s.profiler.measureCuda(ProfileSection::SetDt, [&]() {
            setDT_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var, /*adaptDt=*/true, /*printCfl=*/printCflDt);
        });
        // (4) 案C (coupling 2): 擬似刻み確定後、BDF 込み・ピン除去済み残差で予測し EOS クロス項を流れ RHS へ (周期は独立バッファで gather)。
        if (eosCoupled) {
            s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
                speciesUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // 擬似時間の始点 roY_N = 現在の反復値
                speciesEOSCrossPredictInject_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            });
        }
        // (5) 流れ block 解 → in-place commit（roN=Q^n は BDF 基準で固定のため roN+dq は使えない）。
        s.cfg.dualTimeSubIter = m;
        passiveSaveRhoPre_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // 受動種の φ_N δρ 項用 (#19)
        blockDPLURSolve(s, m);
        // 診断 (FORGE_RESID_SNAP=1): subiter 0 の res (BDF 込み R*) と最終 dq を退避 (局所収縮率 g の分母)。
        {
            // FORGE_RESID_SNAP=<m>: 退避する subiter 番号 (既定 0 相当は "0"。未設定なら無効)。
            static const int snapM = [](){ const char* e = getenv("FORGE_RESID_SNAP"); return e ? atoi(e) : -1; }();
            if (snapM >= 0 && m == snapM) {
                const size_t nb = (size_t)s.msh.nCells_all * sizeof(flow_float);
                const char* src[] = {"res_ro","res_roUx","res_roUy","res_roUz","res_roe",
                                     "dq_block_old_0","dq_block_old_1","dq_block_old_2","dq_block_old_3","dq_block_old_4"};
                const char* dst[] = {"res_ro_m","res_roUx_m","res_roUy_m","res_roUz_m","res_roe_m",
                                     "dq_ro_new","dq_roUx_new","dq_roUy_new","dq_roUz_new","dq_roe_new"};
                for (int q = 0; q < 10; ++q)
                    gpuErrchk(cudaMemcpy(s.var.c_d[dst[q]], s.var.c_d[src[q]], nb, cudaMemcpyDeviceToDevice));
            }
        }
        s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
            applyBlockImplicitCorrectionInPlace_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        });
        // FORGE_FREEZE_TURB=1 は dual-time でも SST 状態更新を凍結する (2026-09-03 修正: 従来この
        // 経路は無条件更新で、freeze 診断が定常専用だった — dual-time A/B は無効だった)。
        if (include_scalar && !freezeTurbEnabled()) {
            s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
                // sstEnergyIncludesK: dual-time は addUnsteadyTimeTerm でエネルギー行の BDF に ρk を含めるので、ここでの roe 補正は不要
                applySSTPointImplicit(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
                periodicMirrorScalarState_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var); // §4.5 k/ω 周期ミラー
            });
        } else if (include_scalar) {
            static bool logged = false;
            if (!logged) { printf("[FREEZE_TURB] dual-time: SST state update frozen\n"); logged = true; }
        }
        // (6) 化学種更新 (BDF 込み res/transport_diag) → 再正規化 → 周期ミラー → primitive → 入口 Dirichlet の再適用。
        if (haveSpecies) {
            s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
                if (eosCoupled) {
                    speciesEOSFinalCommit_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // ρY = ρY_N + z + Y_N (ρ − ρ_pred)
                } else {
                    speciesUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);      // roY_N = 現在の反復値
                    if (speciesImplicitCoupled(s.cfg, s.var)) {
                        speciesImplicitDPLURSolve_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // coupling 1 (定常経路と同じ sweep)
                    } else {
                        speciesTimeIntegration_d_wrapper(0, s.cfg , s.cuda_cfg , s.msh , s.var);   // coupling 0 (speciesImplicitRelax)
                    }
                }
                speciesRenormalize_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);            // ΣρY_s = ρ
                periodicMirrorSpeciesState_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
                speciesPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);              // Y = ρY/ρ (次の残差・出力用)
                applySpeciesBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);                  // 入口 Dirichlet (node ピン値) の再適用
            });
        }
        // (7) 受動種 (凝縮モーメント・トレーサ): scheme 1 は BDF 込みの res/transport_diag で更新クランプ / point-implicit / DPLUR
        //     → 上下限と補正収支 → 周期ミラー (各 wrapper 内)。scheme 0 (旧経路) は従来どおり物理時間項なし。
        s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
            condensationUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            condensationTimeIntegration_d_wrapper(0, s.cfg , s.cuda_cfg , s.msh , s.var);
            condensationPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            tracerUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            tracerTimeIntegration_d_wrapper(0, s.cfg , s.cuda_cfg , s.msh , s.var);
            tracerPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            applyCondensationBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);   // 入口 Dirichlet の再適用 (dry=0 / Xi)
            applyTracerBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);
        });
        pinRowDiagnosticState(s, m);   // FORGE_PIN_DIAG=1: 入口ノード値が bvar 入口値のままか
    }

    // (8) 受動種の物理 step 末尾の保存的 FCT 補正 (plan species-passive-scalar-unification §4.7; SLAU S3 のみ作動):
    //     低次陰解 q_L を限界に、収束した HO 解から制限した反拡散を面共有 α で落とす → 入口 Dirichlet 再適用 → floor (収支) → 実現可能性 → primitive。
    if (passiveFctActive(s.cfg)) {
        // 終了状態で残差を再評価 (ṁ, P_face, ソース, 依存変数を q_H で固定; 受動種の res_* = 空間 HO 残差 → r_H 診断) してから補正。
        assembleResidual(s, 1);
        s.cfg.dualTimeSubIter = -1;   // 物理 step 末尾の印 (トレーサ floor は sub-iter 内では掛けない)
        s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
            passiveFctCorrect_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var, a, b, c);
            applyCondensationBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);
            applyTracerBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);
            passiveBounds_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var, 0, passive_count(), true);
            passiveMirrorPeriodic_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            condensationPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // g の上限・消滅・非負 (射影は下で EOS 更新後)
            tracerPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        });
    }
    // 物理 step 末尾: 二相 EOS (T, P) を更新 → その T で Q1/Q2 の実現可能性射影 (dual-time では sub-iter 内で射影しない; plan §4.7 v7) → FCT の流束形履歴。
    // 受動種経路 (passiveScalarScheme 1) 限定 (codex result-3 M3): 旧経路では**この後処理ブロック自体を回さない**
    // (射影だけでなく、ここで追加で呼ぶ dependentVariables [EOS 再更新] も step 末の状態を変えてしまうため)。
    if (passiveFctActive(s.cfg) || (condensationEnabled(s.cfg) && s.cfg.condRealizProject != 0 && s.cfg.passiveScalarScheme == 1)) {
        s.profiler.measureWall(ProfileSection::DependentVariables, [&]() {
            dependentVariables(s.cfg , s.cuda_cfg , s.msh , s.var, s.mat_ns);
            condensationRealizabilityProject_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        });
        passiveFctFinishHistory_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var, a, b, c);
    }
    s.cfg.unsteadyDiagCoef = 0.0; // 定常側へ影響しないようリセット
    s.cfg.totalTimeD += (double)s.cfg.dt; s.cfg.totalTime = static_cast<flow_float>(s.cfg.totalTimeD);
    s.cfg.nHistoryValid = std::min(s.cfg.nHistoryValid + 1, 2);   // 物理 step 完了: 次 step から Q^{n-1} が有効 (全系共有)

    s.profiler.measureWall(ProfileSection::WriteOutputs, [&]() {
        writeStepOutputs(s.cfg , s.cuda_cfg , s.msh , s.var , s.pprobes , s.iStep+1);
    });
    s.residual_logger.logOuterEnd(s.iStep);
}

// detectNaN==1 のときのみ呼ぶ診断ルーチン。保存量 (ro,roUx,roUy,roUz,roe) と圧力 P
// (RANS 時は roK,roOmega も) を内部セルにわたって検査し、NaN/Inf があれば解を res_nan_<step>.h5 に
// 強制ダンプして即停止する。off のときは一切呼ばれないので通常実行の性能・結果には影響しない。
// GPU 経路は fused 1 カーネルで device int フラグへ集約し、detectNaNInterval ステップごとにのみ
// フラグを host 読み出しする (per-step 同期を避ける)。検知時のみ重い per-var 特定経路に入る。
void checkNonFiniteAndHalt(StepContext& s)
{
    std::vector<std::string> names = {"ro", "roUx", "roUy", "roUz", "roe", "P"};
    if (scalarResidualEnabled(s.cfg)) {
        names.emplace_back("roK");
        names.emplace_back("roOmega");
    }
    if (s.cfg.transitionEnabled()) {
        names.emplace_back("roGamma");
        names.emplace_back("roReth");
    }

    std::string offending;
    bool bad = false;
    if (s.cfg.gpu == 1) {
        // detectNaNInterval ステップごとにのみ device フラグを host 読み出し (それ以外は同期しない)。
        if (((s.iStep + 1) % s.cfg.detectNaNInterval) != 0) {
            return;
        }
        static DeviceResidualReducer nanReducer;
        static bool nanReducerInited = false;
        if (!nanReducerInited) {
            std::vector<std::string> resolved;
            nanReducer = makeDeviceResidualReducer(s.msh, s.var, names, resolved);
            nanReducerInited = true;
        }
        scanNonFiniteToFlag(nanReducer);
        if (downloadNonFiniteFlag(nanReducer) == 0) {
            return;   // 非有限なし: 高頻度で通る軽量経路 (kernel 1 + memcpy 1 のみ)
        }
        // 非有限あり (低頻度): どの変数かを per-var 特定してダンプ・停止する。
        bad = hasNonFiniteCellValue_d(s.msh, s.var, names, offending);
        if (!bad) {
            offending = "(unknown)";   // フラグは立ったが特定できない稀ケースも停止扱い
            bad = true;
        }
    } else {
        for (const std::string& name : names) {
            auto it = s.var.c.find(name);
            if (it == s.var.c.end()) continue;
            for (geom_int ic = 0; ic < s.msh.nCells; ic++) {
                if (!std::isfinite(it->second[ic])) { bad = true; offending = name; break; }
            }
            if (bad) break;
        }
    }

    if (!bad) return;

    s.residual_logger.flush();   // buffer 済みの残差を失わないよう停止前に書き出す
    const int dumpStep = s.iStep + 1;
    std::cerr << "[detectNaN] Non-finite value detected in '" << offending
              << "' at step " << dumpStep
              << ". Dumping solution to res_nan_" << dumpStep
              << ".h5 and halting." << std::endl;
    dumpSolutionH5_force(s.cfg, s.msh, s.var, dumpStep, "res_nan_");
    std::exit(EXIT_FAILURE);
}

void advanceOneStep(
    solverConfig& cfg,
    cudaConfig& cuda_cfg,
    mesh& msh,
    matrix& mat_ns,
    variables& var,
    fluct_variables& fluct,
    point_probes& pprobes,
    RuntimeProfiler& profiler,
    ResidualCsvLogger& residual_logger,
    ImplicitDiagLogger& implicit_diag_logger,
    int iStep)
{
    // 旧 "Step : N / Time : t" ヘッダは廃止。console 出力は StepMonitor (main ループ) の 1 行に集約。
    StepContext s{cfg, cuda_cfg, msh, mat_ns, var, fluct, pprobes,
                  profiler, residual_logger, implicit_diag_logger, iStep};

    // 質量流量一定制御 (bodyForceCtrl=1 のみ)。物理ステップ境界で bodyForceX を更新し、
    // ステップ内 (RK 段・dual-time サブ反復) は固定値として扱う。
    bodyForceCtrlUpdate(cfg , msh , var , iStep);

    // 床事象のカウンタ: この step の中の EOS を step 番号 (1 起点) で数える (無効なら何もしない)
    floorEventsSetPhase(FloorEventPhase::Step, iStep + 1);

    profiler.measureWall(ProfileSection::StepTotal, [&]() {
        if (cfg.isImplicit == 1) {
            if (cfg.unsteady == 1) {
                advanceImplicitDualTime(s);
            } else {
                advanceImplicitSteady(s);
            }
        } else {
            advanceExplicitRK(s);
        }
    });

    floorEventsSetPhase(FloorEventPhase::Aux, iStep + 1);   // step の外で呼ばれる EOS は更新事象と別枠

    // ソルバ内 CHT: ステップ完了後に壁温を更新する (次ステップの残差組立ての前に効く)。
    // 行番号でなく「完了した定常ステップ数」で interval を数える (advanceImplicitSteady 経路が主対象)。
    conjugateWall::updateConjugateWalls(cfg , msh , var , iStep);

    // detectNaN 診断モード (既定 off): 保存量+P を検査し NaN/Inf があればダンプして停止する。
    if (cfg.detectNaN == 1) {
        checkNonFiniteAndHalt(s);
    }
}

}

// ---- 診断 D1 (plans/active/condensation-two-phase-default.md §5.1 #4; FORGE_DIAG_TP_FACES=<出力 h5>, 既定 off) -------------------
// 組立の前処理 (assembleResidualPre: 状態射影・EOS・BC・勾配) を 1 回だけ通し、その前後の保存量の差 (クランプ・境界上書き) を記録する。
// 続けて後処理 (assembleResidualPost) を 1 回通す: 面の拡散係数が読む μ_t (vis_turb) は後処理の turbulent_viscosity が作るので、
// 本番の拡散カーネルが見るのと同じ状態にするため。後処理が面計算の入力 (ρ・T・P・μ・ρY・液・Q) を変えないことはバイト比較で記録する。
// その状態から全通常面の OFF/ON 面作用素を評価して h5 に書き (twoPhaseFaceDiag_d_wrapper)、**時間更新へ進まず終了する**
// (nStepOuter に頼らない)。周期境界は二相拡散の包絡外 (plan §4-1) で、面診断は周期集約を再現しないので拒否する。
static int runTpFacesDiag(const char* path, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, matrix& mat_ns, variables& var,
                          fluct_variables& fluct, point_probes& pprobes, RuntimeProfiler& profiler,
                          ResidualCsvLogger& residual_logger, ImplicitDiagLogger& implicit_diag_logger)
{
    bool hasPeriodic = (msh.nPeriodicMembers > 0);
    for (auto& bc : msh.bconds) if (bc.bcondKind == "periodic") hasPeriodic = true;
    if (hasPeriodic) {
        fprintf(stderr, "[tp-faces] refused: periodic boundaries are outside the two-phase diffusion envelope (plan condensation-two-phase-default §4-1) "
                        "and the face diagnostic does not reproduce the node periodic gather\n");
        return EXIT_FAILURE;
    }
    printf("[tp-faces] FORGE_DIAG_TP_FACES=%s: one pre-part + one post-part of the assembly, face operators A (OFF copy) / B (ON float + double), then exit without update\n", path);
    StepContext s{cfg, cuda_cfg, msh, mat_ns, var, fluct, pprobes, profiler, residual_logger, implicit_diag_logger, 0};
    const size_t nca = (size_t)msh.nCells_all, nc = (size_t)msh.nCells;
    auto present = [&](const std::string& k) { auto it = var.c_d.find(k); return it != var.c_d.end() && it->second != nullptr; };
    auto d2h = [&](const std::vector<std::string>& keys) {
        std::map<std::string, std::vector<flow_float>> m;
        for (auto& k : keys) { auto& v = m[k]; v.resize(nca); gpuErrchk( cudaMemcpy(v.data(), var.c_d[k], nca*sizeof(flow_float), cudaMemcpyDeviceToHost) ); }
        return m;
    };
    // 前後差: 内点 [0, nCells) と ghost [nCells, nCells_all) を分けて、変わった要素数と最大絶対差 (NaN の出入りも 1 件と数える)
    struct Diff { std::vector<std::string> name; std::vector<long long> nReal, nGhost; std::vector<double> mReal, mGhost; };
    auto diffOf = [&](const std::vector<std::string>& keys, const std::map<std::string, std::vector<flow_float>>& a,
                      const std::map<std::string, std::vector<flow_float>>& b) {
        Diff D;
        for (auto& k : keys) {
            const auto& x = a.at(k); const auto& y = b.at(k);
            long long nr = 0, ng = 0; double mr = 0.0, mg = 0.0;
            for (size_t i = 0; i < nca; ++i) {
                if (std::memcmp(&x[i], &y[i], sizeof(flow_float)) == 0) continue;
                const double d = std::fabs((double)y[i] - (double)x[i]);
                if (i < nc) { ++nr; if (std::isnan(d) || d > mr) mr = d; } else { ++ng; if (std::isnan(d) || d > mg) mg = d; }
            }
            D.name.push_back(k); D.nReal.push_back(nr); D.nGhost.push_back(ng); D.mReal.push_back(mr); D.mGhost.push_back(mg);
        }
        return D;
    };
    // 保存量 (存在するものだけ)
    std::vector<std::string> cons;
    for (const char* k : {"ro","roUx","roUy","roUz","roe","roK","roOmega","roXi","roGamma","roReth"}) if (present(k)) cons.push_back(k);
    for (int k = 0; k < var.nSpeciesRegistered; ++k) if (present("roY" + std::to_string(k))) cons.push_back("roY" + std::to_string(k));
    for (const auto& k : var.condMomentConsNames) if (present(k)) cons.push_back(k);
    const auto consBefore = d2h(cons);
    assembleResidualPre(s);
    const auto consAfter = d2h(cons);
    const Diff preDiff = diffOf(cons, consBefore, consAfter);

    // 面計算の入力 (後処理の前後で変わらないこと) と μ_t (後処理が作る; 変わるのが正常)
    std::vector<std::string> inp;
    for (const char* k : {"ro","T","P","vis_lam","vis_turb","ccx","ccy","ccz"}) if (present(k)) inp.push_back(k);
    for (int k = 0; k < var.nSpeciesRegistered; ++k) if (present("roY" + std::to_string(k))) inp.push_back("roY" + std::to_string(k));
    for (const auto& k : var.condMomentConsNames) if (present(k)) inp.push_back(k);
    const auto inpBefore = d2h(inp);
    assembleResidualPost(s);
    const auto inpAfter = d2h(inp);
    const Diff postDiff = diffOf(inp, inpBefore, inpAfter);

    TpFaceDiagData d; std::string why;
    if (!twoPhaseFaceDiag_d_wrapper(cfg, cuda_cfg, msh, var, d, why)) {
        fprintf(stderr, "[tp-faces] refused: %s\n", why.c_str());
        return EXIT_FAILURE;
    }

    // h5
    {
        HighFive::File h5(path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
        auto writeBlk = [&](const TpFaceDiagData::Block& B, size_t rows) {
            if (rows == 0) return;
            const std::vector<size_t> dims = (B.width == 1) ? std::vector<size_t>{rows} : std::vector<size_t>{rows, (size_t)B.width};
            if (B.kind == 0)      { auto ds = h5.createDataSet<float>("/" + B.name, HighFive::DataSpace(dims));  ds.write_raw(B.f.data()); }
            else if (B.kind == 1) { auto ds = h5.createDataSet<double>("/" + B.name, HighFive::DataSpace(dims)); ds.write_raw(B.d.data()); }
            else                  { auto ds = h5.createDataSet<int>("/" + B.name, HighFive::DataSpace(dims));    ds.write_raw(B.i.data()); }
        };
        for (const auto& B : d.face) writeBlk(B, (size_t)d.nFaces);
        for (const auto& B : d.node) writeBlk(B, (size_t)d.nNodes);
        auto writeDiff = [&](const std::string& g, const Diff& D) {
            h5.createDataSet("/" + g + "/names", D.name);
            h5.createDataSet("/" + g + "/count_real", D.nReal); h5.createDataSet("/" + g + "/maxabs_real", D.mReal);
            h5.createDataSet("/" + g + "/count_ghost", D.nGhost); h5.createDataSet("/" + g + "/maxabs_ghost", D.mGhost);
        };
        writeDiff("pre", preDiff);
        writeDiff("post_inputs", postDiff);
        const int iw = d.iw, ns = d.nSpecies, nq = 3 /* Q2, Q1, Q0 (TP_NQ) */, dm = cfg.speciesDiffusionMethod, req = cfg.condTwoPhaseDiffusion;
        const long long nF = d.nFaces, nN = d.nNodes;
        const double Sc = cfg.Sc, Sct = cfg.Sc_t, eps32 = 1.1920928955078125e-7;
        const int hasWd = var.c_d.count("wall_dist") ? 1 : 0;
        const std::string sign = "J > 0 = into cell ic0 (face-integrated; res[ic0] += J, res[ic1] -= J)";
        const std::string skip = "face/skip: 0 evaluated, 1 node boundary half-face (diffusion not added by either operator), 2/3 OFF/ON skip disagreement";
        h5.createAttribute("nSpecies", ns); h5.createAttribute("iw", iw); h5.createAttribute("nMoments", nq);
        h5.createAttribute("nFaces", nF); h5.createAttribute("nNodes", nN);
        h5.createAttribute("speciesDiffusionMethod", dm); h5.createAttribute("Sc", Sc); h5.createAttribute("Sc_t", Sct);
        h5.createAttribute("condTwoPhaseDiffusion_requested", req);
        h5.createAttribute("twophase_diffusion_state", cfg.condTwoPhaseDiffusionState);
        h5.createAttribute("discretization", cfg.discretization);
        h5.createAttribute("input_value_file", cfg.valueFileName);
        h5.createAttribute("species_names", cfg.speciesNames);
        h5.createAttribute("eps32", eps32); h5.createAttribute("has_wall_dist", hasWd);
        h5.createAttribute("sign_convention", sign); h5.createAttribute("skip_codes", skip);
        h5.createAttribute("summary", d.summary);
    }
    for (size_t i = 0; i < preDiff.name.size(); ++i)
        printf("[tp-faces] pre-part change %-10s real %lld (max|d| %.3e)  ghost %lld (max|d| %.3e)\n", preDiff.name[i].c_str(),
               preDiff.nReal[i], preDiff.mReal[i], preDiff.nGhost[i], preDiff.mGhost[i]);
    for (size_t i = 0; i < postDiff.name.size(); ++i)
        if (postDiff.nReal[i] + postDiff.nGhost[i] > 0)
            printf("[tp-faces] post-part changed face input %-10s real %lld  ghost %lld%s\n", postDiff.name[i].c_str(),
                   postDiff.nReal[i], postDiff.nGhost[i], (postDiff.name[i] == "vis_turb") ? " (expected: mu_t is built in the post-part)" : "  ** unexpected **");
    for (const auto& l : d.summary) printf("%s\n", l.c_str());
    printf("[tp-faces] wrote %s (%ld faces, %ld nodes); exiting without any update\n", path, d.nFaces, d.nNodes);
    fflush(stdout);
    return 0;
}

// ---- 診断 G3-b (plans/active/condensation-two-phase-default.md §5.1 #4g3; FORGE_DIAG_TP_UPDATE=<出力 h5>, 既定 off) ---------------
// 更新写像の収支: 再開状態から組立を 1 回 (前処理・後処理) 通し、前処理の格納差 (射影・クランプ・境界の上書き) を「前処理収支」として記録する。
// 続けて本番と同じ設定で**外側の陰的更新を 1 回だけ**行い (implicitNonlinearUpdate; 組立は上で済んでいるので飛ばす)、
// 本番カーネルが実際に読んだ値・書いた値を成分 (ρY_w, ρg, ρQ2, ρQ1, ρQ0) ごとに節点スロットへ書かせる
// (commit の before・δq_limited・格納前の候補・after、各補正操作の before/after; 配置は twoPhaseUpdateDiag_d.cuh)。
// 操作の境界ごとに格納値の写しも取る (記録点の欠落を後処理の E_bookkeeping で見るため)。h5 を書いて終了する (2 回目の更新・res 出力はしない)。
// 積算・判定は notes/investigations/2026-10-04-twophase-g3/g3b_judge.py (double)。周期境界は包絡外なので拒否する。
static int runTpUpdateDiag(const char* path, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, matrix& mat_ns, variables& var,
                           fluct_variables& fluct, point_probes& pprobes, RuntimeProfiler& profiler,
                           ResidualCsvLogger& residual_logger, ImplicitDiagLogger& implicit_diag_logger)
{
    bool hasPeriodic = (msh.nPeriodicMembers > 0);
    for (auto& bc : msh.bconds) if (bc.bcondKind == "periodic") hasPeriodic = true;
    if (hasPeriodic) {
        fprintf(stderr, "[tp-update] refused: periodic boundaries are outside the two-phase diffusion envelope (plan condensation-two-phase-default §4-1)\n");
        return EXIT_FAILURE;
    }
    {
        std::string why;
        if (!tpuBegin(cfg, msh, var, why)) { fprintf(stderr, "[tp-update] refused: %s\n", why.c_str()); return EXIT_FAILURE; }
    }
    const bool twoPhase = condTwoPhaseDiffusionActive(cfg);
    printf("[tp-update] FORGE_DIAG_TP_UPDATE=%s: one pre-part + one post-part of the assembly, then exactly one outer implicit update (%s path), then exit\n",
           path, twoPhase ? "two-phase ON" : "two-phase OFF");
    StepContext s{cfg, cuda_cfg, msh, mat_ns, var, fluct, pprobes, profiler, residual_logger, implicit_diag_logger, 0};

    // 前処理収支 (前処理のバッファ): init → pre_state (状態の射影・クランプ) → pre_end (境界の上書き) → post_end (後処理は状態を変えない)
    tpuArm(1);
    tpuSnap("init");
    assembleResidualPre(s);
    tpuSnap("pre_end");
    assembleResidualPost(s);
    tpuSnap("post_end");
    // 更新写像 (更新のバッファ): upd_start (= q_before) → 本番の 1 更新 (各操作の境界で写し) → final (= q_after)
    tpuArm(2);
    tpuSnap("upd_start");
    g_tpuSkipAssemblyOnce = true;
    implicitNonlinearUpdate(s, 0);
    applyNodeIsothermalWallPin(cfg , cuda_cfg , msh , var);   // advanceImplicitSteady と同じ (この 5 成分には触れない)
    tpuSnap("final");
    tpuArm(0);

    TpuDiagData d;
    if (!tpuCollect(d)) { fprintf(stderr, "[tp-update] internal error: no diagnostic state\n"); return EXIT_FAILURE; }
    const size_t n = (size_t)d.n;
    auto d2hReal = [&](const std::string& k) {
        std::vector<float> v(n);
        gpuErrchk( cudaMemcpy(v.data(), var.c_d[k], n*sizeof(flow_float), cudaMemcpyDeviceToHost) );
        return v;
    };
    auto present = [&](const std::string& k) { auto it = var.c_d.find(k); return it != var.c_d.end() && it->second != nullptr; };

    // スロット名 (h5 の /slot_names; 行 = スロット番号)
    std::vector<std::string> slotNames(TPU_NSLOT);
    for (int op = 0; op < TPU_NOP; ++op)
        for (int c = 0; c < TPU_NC; ++c)
            for (int k = 0; k < 2; ++k) slotNames[TPU_BA(op, c, k)] = std::string(tpuOpName(op)) + "/" + tpuCompName(c) + (k == 0 ? "/before" : "/after");
    const int commitOp[3] = {TPU_OP_SPC, TPU_OP_TPC, TPU_OP_CMC};
    for (int ci = 0; ci < 3; ++ci)
        for (int c = 0; c < TPU_NC; ++c)
            for (int x = 0; x < TPU_NX; ++x) slotNames[TPU_XS(ci, c, x)] = std::string(tpuOpName(commitOp[ci])) + "/" + tpuCompName(c) + "/" + tpuXName(x);
    slotNames[TPU_RN_FACTOR] = "RN_factor";
    {   // 診断 #4pj: 二相 DPLUR の分母 (本番 / sweep が使った値) と射影が使う T・ρ_l・Q3・分岐
        const char* rowName[5] = {"v", "g", "Q2", "Q1", "Q0"};
        for (int q = 0; q < 5; ++q) {
            slotNames[TPU_DIAG(q, 0)] = std::string("DPLUR_denom/") + rowName[q] + "/production";
            slotNames[TPU_DIAG(q, 1)] = std::string("DPLUR_denom/") + rowName[q] + "/used";
        }
        slotNames[TPU_RZ_T] = "RZP_input/T";
        slotNames[TPU_RZ_RHOL] = "RZP_input/rho_l";
        slotNames[TPU_RZ_Q3] = "RZP_input/Q3";
        slotNames[TPU_RZ_KIND] = "RZP_input/kind";
    }

    // 写しの照合 (要約): post_end と upd_start は一致するはず、更新全体の Σ V Δq と Σ V C_round (commit ごと) を出す
    std::vector<float> vol = d2hReal("volume");
    auto snapIdx = [&](const char* lab) { for (size_t i = 0; i < d.snapLabel.size(); ++i) if (d.snapLabel[i] == lab) return (int)i; return -1; };
    const int iS = snapIdx("upd_start"), iF = snapIdx("final"), iP = snapIdx("post_end");
    std::vector<std::string> summary;
    {
        char buf[512];
        long nPostDiff = 0;
        if (iP >= 0 && iS >= 0) for (size_t j = 0; j < (size_t)TPU_NC*n; ++j) if (std::memcmp(&d.snap[iP][j], &d.snap[iS][j], sizeof(float)) != 0) ++nPostDiff;
        snprintf(buf, sizeof(buf), "[tp-update] post-part changed %ld stored values of the 5 components (expected 0)", nPostDiff);
        summary.push_back(buf);
        for (int c = 0; c < TPU_NC; ++c) {
            double lhs = 0.0, dsum = 0.0, cr = 0.0; long nLost = 0;
            for (size_t i = 0; i < n; ++i) {
                const double V = (double)vol[i];
                if (iS >= 0 && iF >= 0) lhs += V*((double)d.snap[iF][(size_t)c*n + i] - (double)d.snap[iS][(size_t)c*n + i]);
                for (int ci = 0; ci < 3; ++ci) {
                    const double dd = d.upd[(size_t)TPU_XS(ci, c, TPU_X_DELTA_D)*n + i];
                    if (!std::isfinite(dd)) continue;
                    if (twoPhase && c == 0 && ci == 0) continue;   // ON の水: 化学種の commit は戻されるので δ に数えない (hold と対で別表示)
                    dsum += V*dd;
                    cr += V*(d.upd[(size_t)TPU_XS(ci, c, TPU_X_CAND_F)*n + i] - d.upd[(size_t)TPU_XS(ci, c, TPU_X_CAND_D)*n + i]);
                    if (d.upd[(size_t)TPU_XS(ci, c, TPU_X_LOST)*n + i] != 0.0) ++nLost;
                }
            }
            snprintf(buf, sizeof(buf), "[tp-update] %-2s  sum V dq = %+.6e  sum V delta_limited = %+.6e  sum V C_round = %+.6e  increments lost to storage: %ld nodes",
                     tpuCompName(c), lhs, dsum, cr, nLost);
            summary.push_back(buf);
        }
    }

    {
        HighFive::File h5(path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
        { auto ds = h5.createDataSet<double>("/pre/slots", HighFive::DataSpace({(size_t)TPU_NSLOT, n}));    ds.write_raw(d.pre.data()); }
        { auto ds = h5.createDataSet<double>("/update/slots", HighFive::DataSpace({(size_t)TPU_NSLOT, n})); ds.write_raw(d.upd.data()); }
        h5.createDataSet("/slot_names", slotNames);
        {
            const size_t ns = d.snap.size();
            std::vector<float> all(ns*(size_t)TPU_NC*n);
            for (size_t k = 0; k < ns; ++k) std::memcpy(all.data() + k*(size_t)TPU_NC*n, d.snap[k].data(), (size_t)TPU_NC*n*sizeof(float));
            auto ds = h5.createDataSet<float>("/snap/data", HighFive::DataSpace({ns, (size_t)TPU_NC, n})); ds.write_raw(all.data());
            h5.createDataSet("/snap/labels", d.snapLabel);
            h5.createDataSet("/snap/phase", d.snapPhase);
        }
        h5.createDataSet("/node/volume", vol);
        for (const char* k : {"ccx", "ccy", "ccz", "wall_dist", "scalarDirichletPin"})
            if (present(k)) h5.createDataSet(std::string("/node/") + k, d2hReal(k));
        std::vector<std::string> comps, ops;
        for (int c = 0; c < TPU_NC; ++c) comps.push_back(tpuCompName(c));
        for (int op = 0; op < TPU_NOP; ++op) ops.push_back(tpuOpName(op));
        h5.createAttribute("components", comps);
        h5.createAttribute("ops", ops);
        h5.createAttribute("nNodes", (long)n);
        h5.createAttribute("nSlots", (int)TPU_NSLOT);
        h5.createAttribute("iw", d.iw);
        h5.createAttribute("twophase_active", twoPhase ? 1 : 0);
        h5.createAttribute("twophase_diffusion_state", cfg.condTwoPhaseDiffusionState);
        h5.createAttribute("condTwoPhaseSolver", cfg.condTwoPhaseSolver);
        h5.createAttribute("condTwoPhaseNonnegLimit", cfg.condTwoPhaseNonnegLimit);
        h5.createAttribute("condTwoPhaseRelax", cfg.condTwoPhaseRelax);
        h5.createAttribute("cfl_pseudo", (double)cfg.cfl_pseudo);
        h5.createAttribute("implicitRelax", (double)cfg.implicitRelax);
        h5.createAttribute("nStepInner", cfg.nStepInner);
        h5.createAttribute("passiveImplicitRelax", (double)cfg.passiveImplicitRelax);
        h5.createAttribute("input_value_file", cfg.valueFileName);
        h5.createAttribute("species_names", cfg.speciesNames);
        h5.createAttribute("eps32", 1.1920928955078125e-7);
        h5.createAttribute("unit_roundoff_u", 5.9604644775390625e-8);
        h5.createAttribute("note_identity",
            std::string("per op: after - before (double of stored floats); commit ops: after - before = delta_d + C_round + C_floor, ")
            + "C_round = cand_f - cand_d, cand_d = before + delta_d, C_floor = after - cand_f; NaN slot = op not on this path");
        h5.createAttribute("note_fma",
            std::string("delta_f_unfused is a non-fused float evaluation for reference only; the production expression (e.g. rg + Th*dg) may be ")
            + "contracted to FMA, so no separately rounded float increment exists; delta_d (exact products of the float operands in double) is the budget term");
        h5.createAttribute("summary", summary);
        // 診断 #4pj: 分母の共通化 (B) の有無と、追加スロットの意味
        h5.createAttribute("common_diag", d.commonDiag);
        h5.createAttribute("note_common_diag",
            std::string("common_diag 1 (FORGE_DIAG_TP_COMMON_DIAG=1): in the two-phase DPLUR (tp_dplur_prep_d) the denominators of g, Q2, Q1, Q0 ")
            + "were replaced by the per-node max(D_g, D_Q2, D_Q1, D_Q0); DPLUR_denom/<row>/production = tp_denoms, /used = what the sweep used "
            + "(rows: v vapour, g, Q2, Q1, Q0; equal when common_diag 0). NaN = DPLUR path not taken");
        h5.createAttribute("note_rzp_input",
            std::string("RZP_input/*: written by the float realizability clamp for every node: T, rho_l = cond_rho_cond(T) (same as the projection), ")
            + "Q3 = rhog/((4/3) pi rho_l) with rhog after the lower/upper clamp; kind: -1 projection branch not entered (rhoQ0 <= 0 or rhog <= 0; "
            + "g = 0 dust is left to droplet removal), 0 unchanged, 1 nearest-point projection, 2 degenerate monodisperse re-init, "
            + "3 Q3 not representable -> Q1 = Q2 = 0. The update buffer holds the clamp inside the one update; the pre buffer the assembly pre-part");
    }
    for (const auto& l : summary) printf("%s\n", l.c_str());
    printf("[tp-update] snapshots:");
    for (size_t i = 0; i < d.snapLabel.size(); ++i) printf(" %s", d.snapLabel[i].c_str());
    printf("\n[tp-update] wrote %s (%ld nodes, %d slots); exiting after one update\n", path, d.n, TPU_NSLOT);
    fflush(stdout);
    return 0;
}

// ---- 診断 G3-a (plans/active/condensation-two-phase-default.md §5.1 #4g3・#4pjg; FORGE_DIAG_TP_OPERATOR=<出力 h5>, 既定 off) -------
// 作用素の収支 (0 step): 再開状態から組立を 1 回 (前処理・後処理) 通し、後処理の中で本番カーネルが実際に残差へ足した値を記録する:
//   面: 移流 (化学種・凝縮モーメント; S3 面値) と拡散 (OFF: species_diffusion_d / ON: twophase_diffusion_d) の atomic の前の値と所属 (ic0, ic1)、
//   節点: 凝縮ソースの実際に足した S·V (分岐・成長分岐の Sg クリップ・クリップ前の Sg・θ)、
//   記録点ごとの残差の写し (speciesTransport のゼロ化後・移流後・OFF 拡散後、化学ソース後、1 回目のピン後、モーメントのゼロ化後・移流後・
//   二相拡散後、凝縮ソース後、トレーサ後、受動種のピン後、2 回目のピン後、最終)。
// 更新へ進まず h5 を書いて終了する (nStepOuter に頼らない)。積算・判定は notes/investigations/2026-10-04-twophase-g3/g3a_judge.py (double)。
// 周期境界は包絡外 (集約後の残差を全 member で二重に数えてしまう) なので拒否する。化学ソースは記録していないので有効なら拒否する。
static int runTpOperatorDiag(const char* path, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, matrix& mat_ns, variables& var,
                             fluct_variables& fluct, point_probes& pprobes, RuntimeProfiler& profiler,
                             ResidualCsvLogger& residual_logger, ImplicitDiagLogger& implicit_diag_logger)
{
    bool hasPeriodic = (msh.nPeriodicMembers > 0);
    for (auto& bc : msh.bconds) if (bc.bcondKind == "periodic") hasPeriodic = true;
    if (hasPeriodic) {
        fprintf(stderr, "[tp-operator] refused: periodic boundaries are outside the two-phase diffusion envelope (plan condensation-two-phase-default §4-1) "
                        "and the node periodic gather would double-count the summed residual\n");
        return EXIT_FAILURE;
    }
    if (chemistryEnabled(cfg)) {
        fprintf(stderr, "[tp-operator] refused: finite-rate chemistry is on; the chemistry source is not instrumented (its residual increment could not be closed)\n");
        return EXIT_FAILURE;
    }
    {
        std::string why;
        if (!tpoBegin(cfg, msh, var, why)) { fprintf(stderr, "[tp-operator] refused: %s\n", why.c_str()); return EXIT_FAILURE; }
    }
    const bool twoPhase = condTwoPhaseDiffusionActive(cfg);
    printf("[tp-operator] FORGE_DIAG_TP_OPERATOR=%s: one pre-part + one post-part of the assembly (%s path), recording the residual terms; exit without update\n",
           path, twoPhase ? "two-phase ON" : "two-phase OFF");
    StepContext s{cfg, cuda_cfg, msh, mat_ns, var, fluct, pprobes, profiler, residual_logger, implicit_diag_logger, 0};
    assembleResidualPre(s);
    tpoArm(true);
    assembleResidualPost(s);
    tpoArm(false);

    TpoDiagData d;
    if (!tpoCollect(msh, d)) { fprintf(stderr, "[tp-operator] internal error: no diagnostic state\n"); return EXIT_FAILURE; }
    const size_t n = (size_t)d.n, nF = (size_t)d.nF, nC = (size_t)d.nComp;
    auto present = [&](const std::string& k) { auto it = var.c_d.find(k); return it != var.c_d.end() && it->second != nullptr; };
    auto d2hReal = [&](const std::string& k) {
        std::vector<float> v(n);
        gpuErrchk( cudaMemcpy(v.data(), var.c_d[k], n*sizeof(flow_float), cudaMemcpyDeviceToHost) );
        return v;
    };

    // 記録の網羅 (経路にある kernel が全面を書いたか、ソースが全節点を書いたか)
    std::vector<std::string> summary;
    bool coverageOk = true;
    {
        char buf[512];
        for (int k = 0; k < TPO_NK; ++k) {
            long cnt[5] = {0, 0, 0, 0, 0};   // −1, 1, 2, 3, その他
            for (size_t ih = 0; ih < nF; ++ih) {
                const int c = d.faceCode[(size_t)k*nF + ih];
                if (c == -1) ++cnt[0]; else if (c == 1) ++cnt[1]; else if (c == 2) ++cnt[2]; else if (c == 3) ++cnt[3]; else ++cnt[4];
            }
            const bool expected = (k == TPO_K_ADV_SP || k == TPO_K_ADV_PA)
                                  || (cfg.viscMethod != 0 && ((k == TPO_K_DIFF_ON) == twoPhase));
            const bool ok = expected ? (cnt[0] == 0 && cnt[4] == 0) : (cnt[0] == (long)nF);
            if (!ok) coverageOk = false;
            snprintf(buf, sizeof(buf), "[tp-operator] face kernel %-17s: %s; not run %ld, internal %ld, half-face skipped %ld, half-face evaluated %ld, other %ld%s",
                     tpoKernelName(k), expected ? "on path" : "off path", cnt[0], cnt[1], cnt[2], cnt[3], cnt[4], ok ? "" : "  ** coverage mismatch **");
            summary.push_back(buf);
        }
        long nNaN = 0, nBr[4] = {0, 0, 0, 0}, nClip = 0;
        for (size_t i = 0; i < n; ++i) {
            const double b = d.src[(size_t)TPO_S_BRANCH*n + i];
            if (!std::isfinite(b)) { ++nNaN; continue; }
            const int bi = (int)b % 10;
            if (bi >= 0 && bi < 4) ++nBr[bi];
            if (d.src[(size_t)TPO_S_CLIP*n + i] == 1.0) ++nClip;
        }
        if (nNaN != 0) coverageOk = false;
        snprintf(buf, sizeof(buf), "[tp-operator] source: nodes not written %ld%s; branch none %ld, growth %ld, evaporation %ld; growth-branch Sg clip fired %ld",
                 nNaN, nNaN ? "  ** coverage mismatch **" : "", nBr[0], nBr[1], nBr[2], nClip);
        summary.push_back(buf);
        std::string lab = "[tp-operator] residual snapshots:";
        for (const auto& l : d.snapLabel) lab += " " + l;
        summary.push_back(lab);
    }

    // 面の境界種別 (−1 = 内部面; それ以外は msh.bconds の番号)
    std::vector<int> faceBc(nF, -1);
    {
        std::vector<int> ipBc((size_t)msh.nPlanes, -1);
        for (size_t b = 0; b < msh.bconds.size(); ++b)
            for (geom_int ip : msh.bconds[b].iPlanes) if (ip >= 0 && ip < msh.nPlanes) ipBc[(size_t)ip] = (int)b;
        for (size_t ih = 0; ih < nF; ++ih) faceBc[ih] = ipBc[(size_t)d.faceIp[ih]];
    }
    std::vector<float> mdot(nF);
    {
        std::vector<float> mp((size_t)msh.nPlanes);
        gpuErrchk( cudaMemcpy(mp.data(), var.p_d["massflux"], mp.size()*sizeof(flow_float), cudaMemcpyDeviceToHost) );
        for (size_t ih = 0; ih < nF; ++ih) mdot[ih] = mp[(size_t)d.faceIp[ih]];
    }

    {
        HighFive::File h5(path, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
        h5.createDataSet("/face/ip", d.faceIp);
        h5.createDataSet("/face/ic0", d.faceIc0);
        h5.createDataSet("/face/ic1", d.faceIc1);
        h5.createDataSet("/face/bcond", faceBc);
        h5.createDataSet("/face/massflux", mdot);
        { auto ds = h5.createDataSet<double>("/face/adv", HighFive::DataSpace({nC, nF}));  ds.write_raw(d.faceVal.data()); }
        { auto ds = h5.createDataSet<double>("/face/diff", HighFive::DataSpace({nC, nF})); ds.write_raw(d.faceVal.data() + nC*nF); }
        { auto ds = h5.createDataSet<int>("/face/code", HighFive::DataSpace({(size_t)TPO_NK, nF})); ds.write_raw(d.faceCode.data()); }
        { auto ds = h5.createDataSet<double>("/source/slots", HighFive::DataSpace({(size_t)TPO_NSRC, n})); ds.write_raw(d.src.data()); }
        {
            std::vector<std::string> sn; for (int k = 0; k < TPO_NSRC; ++k) sn.push_back(tpoSrcName(k));
            h5.createDataSet("/source/slot_names", sn);
            std::vector<std::string> kn; for (int k = 0; k < TPO_NK; ++k) kn.push_back(tpoKernelName(k));
            h5.createDataSet("/face/code_kernels", kn);
        }
        {
            const size_t ns = d.snap.size();
            std::vector<float> all(ns*nC*n);
            for (size_t k = 0; k < ns; ++k) std::memcpy(all.data() + k*nC*n, d.snap[k].data(), nC*n*sizeof(float));
            auto ds = h5.createDataSet<float>("/res/snap", HighFive::DataSpace({ns, nC, n})); ds.write_raw(all.data());
            h5.createDataSet("/res/labels", d.snapLabel);
        }
        {
            std::vector<float> st(nC*n);
            for (size_t c = 0; c < nC; ++c) { const auto v = d2hReal(d.compCons[c]); std::memcpy(st.data() + c*n, v.data(), n*sizeof(float)); }
            auto ds = h5.createDataSet<float>("/node/state", HighFive::DataSpace({nC, n})); ds.write_raw(st.data());
        }
        h5.createDataSet("/node/volume", d2hReal("volume"));
        if (msh.volumePartial_d != nullptr) {
            std::vector<float> v(n);
            gpuErrchk( cudaMemcpy(v.data(), msh.volumePartial_d, n*sizeof(geom_float), cudaMemcpyDeviceToHost) );
            h5.createDataSet("/node/volume_partial", v);
        }
        for (const char* k : {"ccx", "ccy", "ccz", "wall_dist", "scalarDirichletPin", "ro", "T", "P",
                              "condS_0", "condLim_0", "condDrdt_0", "condR30_0"})
            if (present(k)) h5.createDataSet(std::string("/node/") + k, d2hReal(k));
        std::vector<std::string> comps = {"w", "g", "Q2", "Q1", "Q0"};
        for (int sp = 0; sp < var.nSpeciesRegistered; ++sp)
            if (sp != d.iw) comps.push_back((sp < (int)cfg.speciesNames.size()) ? cfg.speciesNames[sp] : ("Y" + std::to_string(sp)));
        std::vector<std::string> bcNames, bcKinds;
        for (auto& bc : msh.bconds) { bcNames.push_back(bc.physName); bcKinds.push_back(bc.bcondKind); }
        h5.createAttribute("components", comps);
        h5.createAttribute("comp_cons", d.compCons);
        h5.createAttribute("comp_res", d.compRes);
        h5.createAttribute("bcond_names", bcNames);
        h5.createAttribute("bcond_kinds", bcKinds);
        h5.createAttribute("nNodes", (long)n);
        h5.createAttribute("nFaces", (long)nF);
        h5.createAttribute("nComp", (int)nC);
        h5.createAttribute("iw", d.iw);
        h5.createAttribute("twophase_active", twoPhase ? 1 : 0);
        h5.createAttribute("twophase_diffusion_state", cfg.condTwoPhaseDiffusionState);
        h5.createAttribute("viscMethod", cfg.viscMethod);
        h5.createAttribute("speciesDiffusionMethod", cfg.speciesDiffusionMethod);
        h5.createAttribute("Sc", (double)cfg.Sc);
        h5.createAttribute("Sc_t", (double)cfg.Sc_t);
        h5.createAttribute("condFloat", cfg.condFloat);
        h5.createAttribute("condLimiterMode", cfg.condLimiterMode);
        h5.createAttribute("discretization", cfg.discretization);
        h5.createAttribute("input_value_file", cfg.valueFileName);
        h5.createAttribute("species_names", cfg.speciesNames);
        h5.createAttribute("coverage_ok", coverageOk ? 1 : 0);
        h5.createAttribute("unit_roundoff_u", 5.9604644775390625e-8);
        h5.createAttribute("flt_min", 1.1754943508222875e-38);
        h5.createAttribute("note_face",
            std::string("/face/adv and /face/diff [component, ih]: value a0 added to res[ic0] (res[ic1] receives -a0; face-integrated). ")
            + "adv: a0 = -mdot*Y_f (node boundary half-face, ic1 >= nNodes: added to ic0 only, code 3). diff: a0 = J (into ic0); "
            + "node boundary half-faces are not added (value 0, code 2). NaN = no kernel wrote this component on this face. "
            + "/face/code [kernel, ih]: -1 kernel not on path, 1 internal face evaluated, 2 half-face skipped, 3 half-face evaluated. "
            + "/face/bcond: -1 internal, else index into bcond_names/bcond_kinds. /face/massflux = massflux of the plane (positive ic0 -> ic1)");
        h5.createAttribute("note_source",
            std::string("/source/slots [slot, node] (names /source/slot_names): term_* = value added to res (float path: exact double(S)*double(V), ")
            + "production rounds the product and the sum, nround 2; double delegate: (float)(S*V), nround 1; nothing added: 0 and nround 0). "
            + "branch: 0 none, 1 nucleation/growth, 2 evaporation (float); 10/11/12 the same in the double delegate. clip: growth-branch "
            + "'if (Sg < 0) Sg = 0' fired. Sg_unclipped: Sg before that clip (growth branch only, else NaN). theta: source scale, J: nucleation rate");
        h5.createAttribute("note_snap",
            std::string("/res/snap [label, component, node] (float stored residuals at the recording points, labels /res/labels). ")
            + "Moment residuals are zeroed in passiveAdvection (label cm_zero); species residuals in speciesTransport (sp_zero). "
            + "Pins overwrite (no rounding); every other change between labels is a sum of the recorded terms plus assembly rounding");
        h5.createAttribute("summary", summary);
    }
    for (const auto& l : summary) printf("%s\n", l.c_str());
    printf("[tp-operator] wrote %s (%ld nodes, %ld faces, %d components); exiting without any update%s\n", path, (long)n, (long)nF, (int)nC,
           coverageOk ? "" : "  ** recording coverage mismatch: do not use for the budget **");
    fflush(stdout);
    return coverageOk ? 0 : EXIT_FAILURE;
}

int main(int argc, char** argv) {
    // --resolve-species: 化学種の解決済み記録だけ書いて終了 (GPU 不使用; plan thermophysics-solver-owned-species-db §4.3)
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--resolve-species") return resolveSpeciesOnly();
    }
    // 診断 #4pj (FORGE_DIAG_TP_COMMON_DIAG=1; 既定 off): 分母の共通化は FORGE_DIAG_TP_UPDATE の 1 更新の中だけで使う診断。
    // 本番の計算に混ざらないよう、単独で立っていたら計算を始めずに拒否する。
    if (const char* e = getenv("FORGE_DIAG_TP_COMMON_DIAG"); e != nullptr && *e != '\0' && std::strcmp(e, "0") != 0) {
        const char* u = getenv("FORGE_DIAG_TP_UPDATE");
        if (u == nullptr || *u == '\0') {
            fprintf(stderr, "[tp-update] refused: FORGE_DIAG_TP_COMMON_DIAG=%s is a diagnostic for FORGE_DIAG_TP_UPDATE only "
                            "(plan condensation-two-phase-default #4pj); unset it for production runs\n", e);
            return EXIT_FAILURE;
        }
    }
    RuntimeProfiler profiler;

    solverConfig cfg;
    mesh msh;
    matrix mat_ns;
    variables var;
    fluct_variables fluct;
    point_probes pprobes;
    ImplicitDiagLogger implicit_diag_logger;

    cudaConfig cuda_cfg = initializeSimulation(cfg, msh, mat_ns, var, fluct, pprobes);
    // 試験用: 物性だけを評価して終了 (時間更新なし; runTransportProbe の説明)
    if (const char* e = getenv("FORGE_TRANSPORT_PROBE"); e != nullptr && *e != '\0') {
        return runTransportProbe(e, cfg, cuda_cfg, msh, var);
    }
    if (const char* e = getenv("FORGE_TRANSPORT_TABLE_PROBE"); e != nullptr && *e != '\0') {
        return runTransportTableProbe(e);
    }
    // FORGE_OUT_RESIDUALS / FORGE_RESID_SNAP の出力名の登録は確保の前 (initializeSimulation の registerOutputDiagnostics) へ移した。
    // 保存量の FP64 影アキュムレータ (plans/active/time_integration-fp64-accumulator.md §4.4)。
    // **非対応の経路で明示 ON されたら黙って劣化させず拒否する** (累積が消える経路があるため)。
    if (cfg.qAccumulatorFP64 == 1) {
        // **v1a の対応範囲** (§5.1 S1a)。どれも原理的な制限ではなく「まだ Qacc を扱っていない」だけで、
        // S1b-①〜④ で順に外す。**黙って劣化させるくらいなら拒否する** (codex plan M5)。
        const char* why = nullptr;
        // **検証済みの経路だけ通す** (2026-09-24, codex result-2 M5)。
        // 以前は node / GPU の検査が無く、**cell も CPU 経路も素通り**していた。
        // CPU の commit (`update.cpp` の applyScalarImplicitCorrection 等) は FP32 のままなので、
        // `gpu != 1` で ON にすると device 正本だけ確保されて一切使われない。
        // cell はユーザ方針で使わない (AGENTS.md / [[user-prefers-node-base]]) ため**未検証**。
        if (cfg.gpu != 1)
            why = "GPU 経路 (gpu=1) のみ対応 (CPU の commit は FP32 のまま)";
        else if (cfg.discretization != "node")
            why = "node のみ対応 (cell は未検証。ユーザ方針で cell は使わない)";
        else if (cfg.timeIntegration != 11)
            why = "v1a は timeIntegration=11 のみ (陽解法 tI 1/4 は S1b-④、tI 3 は凸結合なので Qacc_N/Qacc_M が要る)";
        else if (cfg.unsteady != 0)
            why = "v1a は unsteady=0 のみ (dual-time は QaccN/QaccNN の shift が要る。陽解法 unsteady は S1b-④)";
        else if (cfg.isAxisymmetric != 0)
            // enforceAxisSymmetry は commit の**基準** roeN/roUyN を射影する (axisymmetricSource_d.cu:312-323)。
            // Qacc は Q_N を読まないので、その射影を Qacc に当てるまでは対応できない。
            why = "node 軸対称は S1b-① 待ち (軸ピンが commit の基準 roeN/roUyN を射影するため)";
        else if (cfg.sstEnergyIncludesK != 0)
            // ransTransport_d.cu:175 は roe への**増分**なので、Qacc にも同じ増分を当てる必要がある。
            why = "sstEnergyIncludesK=1 は S1b-② 待ち (roe への増分を Qacc にも当てる必要がある)";
        else if (msh.nPeriodicMembers > 0)
            // periodicNode_d.cu:119 は FP32 値だけを root→member に配るので、member の Qacc の
            // 下位ビットが root と食い違ったまま残る (値が一致すると reconcile も発火しない)。
            why = "node 周期は S1b-③ 待ち (root→member ミラーが FP32 値だけを配るため)";
        if (why != nullptr) {
            fprintf(stderr, "[qAccumulatorFP64] 拒否: %s\n", why);
            fprintf(stderr, "[qAccumulatorFP64] v1a の対応: GPU・node・timeIntegration=11 && unsteady=0、block/scalar DPLUR、"
                            "軸対称なし、周期なし、sstEnergyIncludesK=0\n");
            exit(1);
        }
        // 確保と初期化は initializeSimulation の中 (updateVariablesOuter の直後) で隣接して行う。
        // ここでは**それが済んでいることを確かめるだけ**にする (黙って未初期化のまま進ませない)。
        if (var.qacc_d[0] == nullptr) {
            fprintf(stderr, "[qAccumulatorFP64] 内部エラー: 正本が確保されていない "
                            "(initializeSimulation での確保・初期化が走っていない)\n");
            exit(1);
        }
        printf("[qAccumulatorFP64] 有効: 保存量 5 本の正本を FP64 に置く (内点 %ld CV)\n", (long)msh.nCells);
        // **対象は流れの 5 本だけ**。乱流・化学種・凝縮モーメント・受動種・遷移モデルの保存量は
        // **float32 のまま**で、commit の丸めを受ける経路に残る。拒否はしない (case/56 では
        // 全域 FP64 ビルドと物理量が 0.00-3.77 % で一致しており実害が出ていない) が、
        // **混在していることを黙らせない** (2026-09-24)。横展開は
        // plans/active/time_integration-fp64-accumulator-rollout.md。
        {
            std::vector<std::string> notAcc;
            if (cfg.LESorRANS == 2) notAcc.push_back("乱流 (roK, roOmega)");
            if (cfg.nSpecies > 1)   notAcc.push_back("化学種 (roY*)");
            if (var.condMomentConsNames.size() > 0) notAcc.push_back("凝縮モーメント");
            if (var.tracerRegistered != 0)     notAcc.push_back("受動トレーサ (roXi)");
            if (var.transitionRegistered != 0) notAcc.push_back("遷移モデル (roGamma, roReth)");
            if (!notAcc.empty()) {
                printf("[qAccumulatorFP64] **注意**: 次は FP64 正本を持たない (float32 のまま):");
                for (size_t i = 0; i < notAcc.size(); ++i) printf("%s %s", i ? " /" : "", notAcc[i].c_str());
                printf("\n[qAccumulatorFP64]   流れの保存量だけが累積される混在状態になる。\n");
            }
        }
    }

    // line-implicit (plans/active/time_integration-line-implicit.md): 壁法線ラインを構築。
    // blockDPLUR==1 専用・完全前処理 (lowMachPrecond>=2) とは併用不可。
    if (cfg.lineImplicit == 1) {
        if (cfg.blockDPLUR != 1 || cfg.lowMachPrecond >= 2 || cfg.timeIntegration != 11) {
            fprintf(stderr, "[lineImplicit] requires timeIntegration=11, blockDPLUR=1, lowMachPrecond<2\n");
            exit(1);
        }
        // ホストの ccx/ccy/ccz は GPU 経路では書かれない (0 のまま。node はノード座標を使う。監査 §5 危険 5、今の挙動のまま)
        msh.buildImplicitLines(var.hostCell("ccx").data(), var.hostCell("ccy").data(), var.hostCell("ccz").data());
    } else if (cfg.lineKFreeze != 0 || cfg.lineViscCoupling != 0 ||
               cfg.lineViscousDtRelief != (flow_float)0.0 || cfg.lineDtDirectional != 0) {
        fprintf(stderr, "[lineImplicit] lineKFreeze/lineViscCoupling/lineViscousDtRelief require lineImplicit=1\n");
        exit(1);
    }
    ResidualCsvLogger residual_logger("residual_history.csv", cfg, msh, var);
    MEMLOG_SOLVER("after lineImplicit/ResidualCsvLogger");

    // 診断 D1 (FORGE_DIAG_TP_FACES=<出力 h5>; 既定 off): 組立を 1 回だけ通した状態の面作用素 A/B を書いて終了する (時間更新・res_0 出力なし)。
    if (const char* e = getenv("FORGE_DIAG_TP_FACES"); e != nullptr && *e != '\0') {
        return runTpFacesDiag(e, cfg, cuda_cfg, msh, mat_ns, var, fluct, pprobes, profiler, residual_logger, implicit_diag_logger);
    }
    // 診断 G3-a (FORGE_DIAG_TP_OPERATOR=<出力 h5>; 既定 off): 組立 1 回の作用素の記録 (面流束・ソース・残差の写し) を書いて終了する (更新・res_0 出力なし)。
    if (const char* e = getenv("FORGE_DIAG_TP_OPERATOR"); e != nullptr && *e != '\0') {
        return runTpOperatorDiag(e, cfg, cuda_cfg, msh, mat_ns, var, fluct, pprobes, profiler, residual_logger, implicit_diag_logger);
    }
    // 診断 G3-b (FORGE_DIAG_TP_UPDATE=<出力 h5>; 既定 off): 組立 1 回 + 外側の陰的更新 1 回の更新写像の記録を書いて終了する (res_0 出力なし)。
    if (const char* e = getenv("FORGE_DIAG_TP_UPDATE"); e != nullptr && *e != '\0') {
        return runTpUpdateDiag(e, cfg, cuda_cfg, msh, mat_ns, var, fluct, pprobes, profiler, residual_logger, implicit_diag_logger);
    }

    MEMLOG_SOLVER("before writeInitialOutputs");
    writeInitialOutputs(cfg , msh , var);
    MEMLOG_SOLVER("after writeInitialOutputs");

    StepMonitor monitor(cfg, residual_logger);
    monitor.printHeader();
    passiveRecordInitialTotals_d_wrapper(cfg, cuda_cfg, msh, var);   // 収支の独立照合の始点 (計算開始前の総量)
    // 二相拡散 (#4f (4)): 収束受入の独立残差監査。開始時に r0 を記録し、終了時に格納状態から組み直して判定する (無効構成は no-op)。
    auto twoPhaseAudit = [&](int iStepAudit, bool final) {
        if (!condResidualAuditActive(cfg)) return;
        StepContext sa{cfg, cuda_cfg, msh, mat_ns, var, fluct, pprobes, profiler, residual_logger, implicit_diag_logger, iStepAudit};
        assembleResidual(sa, 1);
        twoPhaseAudit_d_wrapper(cfg, cuda_cfg, msh, var, iStepAudit, final);
    };
    twoPhaseAudit(0, false);
    cout << "Start Calculation \n";
    for (int iStep = 0 ; iStep < cfg.mainLoopCount() ; iStep++) {
        advanceOneStep(cfg , cuda_cfg , msh , mat_ns , var , fluct , pprobes , profiler , residual_logger , implicit_diag_logger , iStep);
        floorEventsTestInject(cfg , cuda_cfg , msh , var , iStep + 1);   // 試験用 (FORGE_FLOOR_TEST_INJECT、既定 off)
        monitor.report(iStep);
        // 受動種経路の補正収支 (floor による保存量補正の体積積分; monitorInterval ごと)。scheme 0 / 受動種なしでは no-op。
        if (iStep % cfg.monitorInterval == 0) passiveFloorCorrLog_d_wrapper(cfg, cuda_cfg, msh, var, iStep);
        // 凝縮の理由別の補正量 (plan condensation-two-phase-transport §4.3; 凝縮なしは no-op)
        if (iStep % cfg.monitorInterval == 0) condCorrectionLog_d_wrapper(cfg, cuda_cfg, msh, var, iStep);
        // 1 step 目の後 (陰解法・SST・化学種などの遅延確保はここまでに入る) と 2 step 目の後 (定常に達したか)
        if (iStep == 0) MEMLOG_SOLVER("after step 1");
        if (iStep == 1) MEMLOG_SOLVER("after step 2");
    }
    MEMLOG_SOLVER("after main loop");
    // 床事象のカウンタ: 最後の更新の結果 Q_N は次の EOS が無いので、コピー上で同じ前処理・同じ EOS にかけて監査する
    // (通常の配列には触れない。二相の残差監査など通常の配列で EOS を回す処理より前に置く)
    floorEventsAudit(cfg , cuda_cfg , msh , var , cfg.mainLoopCount());
    // 終了時に受動種の収支を必ず出す (最終 step が monitorInterval に乗らないと末尾の補正が記録されない; plan-8 M1)
    if (cfg.mainLoopCount() > 0 && ((cfg.mainLoopCount() - 1) % cfg.monitorInterval) != 0) {
        passiveFloorCorrLog_d_wrapper(cfg, cuda_cfg, msh, var, cfg.mainLoopCount() - 1);
        condCorrectionLog_d_wrapper(cfg, cuda_cfg, msh, var, cfg.mainLoopCount() - 1);
    }

    twoPhaseDiagWrite(cfg, cuda_cfg, msh, var);   // 二相更新の診断 CSV (condTwoPhaseDiag; 無効なら no-op)
    renormGateLog(cfg, cfg.mainLoopCount(), true);   // 再正規化の受入ゲート: 末尾 ceil(0.1N) 更新の max と VERDICT (#1b-pre)
    twoPhaseCorrGateLog(cfg, cfg.mainLoopCount(), true);   // 二相の補正ゲート: 末尾 ceil(0.1N) 更新の max と VERDICT (#4h)
    twoPhaseAudit(cfg.mainLoopCount(), true);   // 最終の格納状態 (出力は書き終えている)
    limiterDiag_finalize(cfg);   // 有界性診断の末尾取りこぼしを回収して累計を確定 (plan §4.35)

    // 壁時計 (旧実装は clock() = CPU 時間で、GPU 待ちを含まなかった)。書式 "Time = %.3f s" は grep 互換のため維持。
    printf("Time = %.3f s (wall, %d steps, %.2f ms/step)\n", monitor.elapsedSeconds(), cfg.mainLoopCount(),
           monitor.elapsedSeconds() * 1.0e3 / std::max(1, cfg.mainLoopCount())); 
    profiler.printSummary();
    floorEventsClose(cfg.mainLoopCount());   // 正常終了の印 (session_end) を書いて閉じる
    MEMLOG_SOLVER("end");

	return 0;
}