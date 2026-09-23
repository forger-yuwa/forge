#include <algorithm>
#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <stdio.h>                                                                                       
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <limits>
#include <cstdio>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

#include "input/solverConfig.hpp"
#include "input/condSonicResolve.hpp"
#include "input/setInitial.hpp"

#include "mesh/gmshReader.hpp"
#include "output/output.hpp"
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
#include "cuda_forge/tracerTransport_d.cuh"
#include "cuda_forge/passiveTransport_d.cuh"
#include <highfive/H5File.hpp>
#include "input/speciesDB.hpp"
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
#include "cuda_forge/thermo_d.cuh"

#include "probe/point_probes.cuh"
#include "cuda_forge/setDT_d.cuh"
#include "cuda_forge/periodicNode_d.cuh"
#include "cuda_forge/qAccumulator.hpp"

#include <cuda_runtime.h>


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
    const std::list<std::string> variable_names = {
        "ro", "Ux", "Uy", "Uz", "sonic", "vis_turb", "dt_local"
    };

    if (cfg.gpu == 1) {
        var.copyVariables_cell_D2H(variable_names);
    }

    const std::vector<flow_float>& ro = var.c.at("ro");
    const std::vector<flow_float>& ux = var.c.at("Ux");
    const std::vector<flow_float>& uy = var.c.at("Uy");
    const std::vector<flow_float>& uz = var.c.at("Uz");
    const std::vector<flow_float>& sonic = var.c.at("sonic");
    const std::vector<flow_float>& vis_turb = var.c.at("vis_turb");
    const std::vector<flow_float>& dt_local = var.c.at("dt_local");

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

    std::list<std::string> hist = {"roN", "roUxN", "roUyN", "roUzN", "roeN", "roKN", "roOmegaN"};
    for (const auto& nm : var.speciesVarNames) hist.push_back(nm + "P");
    if (cfg.passiveScalarScheme == 1) {   // scheme 0 の受動種は物理時間項を持たないので履歴は要らない (書きもしない)
        if (var.tracerRegistered != 0) hist.push_back("roXiP");
        for (const auto& nm : var.condMomentConsNames) hist.push_back(nm + "P");
    }
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
                    std::vector<flow_float>& v = var.c.at(nm);
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

cudaConfig initializeSimulation(
    solverConfig& cfg,
    mesh& msh,
    matrix& mat_ns,
    variables& var,
    fluct_variables& fluct,
    point_probes& pprobes)
{
    cout << "Read Solver Config \n";
    cfg.read("solverConfig.yaml");

    // 化学種 DB の host 側解決 (GPU 非依存; 未知種名はここで exit)。bcond の X{s}→Y{s} 換算と
    // 起動ログ (種表) が使う。thermo_init_db は同じ結果を device へ上げる。
    speciesDB_printTable(cfg, speciesDB_init(cfg));

    cout << "Init Thermo DB \n";
    thermo_init_db(cfg);   // NASA-9/LJ 化学種 DB を構築し device へアップロード (thermalMethod==2 用)
    chemistry_init(cfg);   // 有限速度化学: 反応機構を読み device へ (chemistry.enabled==1 のみ)

    cout << "Read Mesh \n";
    if (cfg.meshFormat == "hdf5") {
        msh.nodeValueAtNode = (cfg.discretization == "node") ? 1 : 0;   // node は常に値=ノード座標
        msh.readMesh(cfg.meshFileName);
    } else {
        cerr << "Error unknown mesh format: " << cfg.meshFormat << endl;
        std::exit(EXIT_FAILURE);
    }

    cout << "Init Matrix (but not used now) \n";
    mat_ns.initMatrix(msh);

    cout << "Read Boundary Conditions \n";
    readBcondConfig(cfg , msh.bconds);

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

    // 化学種変数を登録 (allocVariables より前)。nSpecies<=1 では no-op。
    var.registerSpecies(cfg.nSpecies, cfg.chemEnabled);

    // 非平衡凝縮モーメント変数を登録 (allocVariables より前)。condensation==0 では no-op。
    var.registerCondensation(cfg.nCondSpecies);

    // 受動トレーサ roXi を登録 (allocVariables より前)。physProp.tracer 未指定では no-op。
    var.registerTracer(cfg.tracerEnabled() ? 1 : 0);

    // 遷移モデル (turbulence.transition: lm2009) の変数を登録。none では no-op。診断場は output.level>=2 のときだけ確保。
    transitionValidateConfig(cfg);
    var.registerTransition(cfg.transitionEnabled() ? 1 : 0, (cfg.outputLevel >= 2) ? 1 : 0);

    var.allocVariables(cfg.gpu , msh);

    // device roY[] ポインタ配列を構築 (c_d 確保後, dependentVariables より前)。
    speciesInit_d(cfg , var);

    // device rog[] ポインタ配列を構築 (二相 EOS が液相質量分率を読む)。condensation==0 で no-op。
    condensationInit_d(cfg , var);
    // 受動種 (トレーサ + 凝縮モーメント) の device ポインタ配列 (passiveScalarScheme 1 の化学種経路用; 0 では表だけ)。
    passiveInit_d(cfg , var);

    cout << "Read Initial Values \n";
    var.readValueHDF5(cfg.valueFileName , msh, cfg.kInit, cfg.omegaInit);

    cout << "Set mesh connection map for cuda \n";
    cudaConfig cuda_cfg(msh);
    msh.setMeshMap_d();
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

    var.setStructuralVariables(cfg , cuda_cfg , msh);
    // node-centered 周期境界 DOF 同一視 (median-dual M4, §4.5): partnerCellID から周期ノード group(union-find)
    // を構築し、各 group の合併体積を var.c_d["volume"] へ書き戻す。setDT より前 (volume を使うため)。
    // cell モード / 非周期では no-op。
    msh.buildPeriodicNodeGroups(cfg.discretization == "node", var.c_d["volume"]);
    // node 周期境界 DOF 同一視 (§4.5.9): 読み込んだ保存量を root→member ミラーで slave=master に初期同期する。
    // seed (2D複製+非周期 Uz 摂動 / restart / 丸め) で周期ペアの保存量が desync していることがあり、残差 gather だけ
    // では desync が永続して継ぎ目フラックス不整合 (seam 圧力欠陥) を生む。dependentVariables より前に揃える。
    periodicMirrorNSState_d_wrapper(cfg , cuda_cfg , msh , var);
    speciesPrimitive_d_wrapper(cfg , cuda_cfg , msh , var);  // Y_s = ρY_s/ρ (roY を読込済)
    condensationPrimitive_d_wrapper(cfg , cuda_cfg , msh , var);  // φ = ρφ/ρ (液相モーメント読込済)
    tracerPrimitive_d_wrapper(cfg , cuda_cfg , msh , var);  // ξ = ρξ/ρ (トレーサ読込済)
    dependentVariables(cfg , cuda_cfg , msh , var, mat_ns);
    // node-centered 壁 Dirichlet: IC の壁ノード速度を厳密 0 に初期化 (KE を roe から除去)。
    // この後 gasProperties が補正 roe から P/T を再計算する。cell/非 node では no-op。
    enforceWallNoSlip_d_wrapper(cfg , cuda_cfg , msh , var);
    gasProperties_d_wrapper(cfg , cuda_cfg , msh , var);

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

    pprobes.init(cfg , cuda_cfg , msh);


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
    outputH5_XDMF(cfg , msh, var, iStep);
    outputBconds_H5_XDMF(cfg , msh, var, iStep);
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
void assembleResidual(StepContext& s, int stage_index)
{
    (void)stage_index;  // 現状カーネルは stage_index を使わない（dual-time 拡張用に interface 保持）
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
        periodicGradientGather_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::AxisymmetricSource, [&]() {
        axisymmetricGeomTerms_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
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
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        // k/ω 勾配と F1 を拡散の**前**に評価する (2026-09-08, plan turbulence-sst-consistency-options §2.1):
        // 旧順序 (transport → gradient → source) では拡散の非直交補正と σ ブレンドの F1 が前回評価の値
        // (Picard ラグ) だった。ransSource の F1 は同式・同入力なので sstF1 と一致する。
        ransGradient_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        ransBlendF1_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        ransTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);
        transitionTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);   // γ / Re_θt の移流拡散 (none で no-op)
    });
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        speciesTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);  // 化学種移流残差
        chemistrySource_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);   // 有限速度化学ソース (ω_s, Q̇, 対角 Jacobian)
        speciesPinResidual_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var); // node 入口ピンノードの化学種残差除外 (cell は no-op)
    });
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        condensationTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);  // 液相モーメント移流残差 (Phase 1)
        condensationSource_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);     // 核生成+成長ソース (Phase 2)
        tracerTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);        // 受動トレーサ移流残差 (node 入口ピン込み)
        passivePinResidual_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);     // 受動種経路: node 入口ピンノードの残差除外 (ソース集計の後)
    });
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
    s.profiler.measureCuda(ProfileSection::ViscousFlux, [&]() {
        viscousFlux_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var, s.mat_ns);
    });
    // SU2 流の軸対称対称面 (MARKER_SYM): 軸上 CV の半径方向運動量残差を 0 に射影し roUy=0 を保つ。
    // explicit では ΔroUy=0 になり直接効く (implicit/block-DPLUR では連成 solve が補正を漏らすため Jacobian
    // 整合が別途必要・open issue, docs §7.1)。cell/非軸対称/平面では no-op。
    zeroAxisRadialResidual_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    // node-centered 壁 Dirichlet: 壁ノードの運動量残差を 0 に射影し u=0 を保つ (壁ゴースト撤廃の代替)。
    // 軸射影の後に置き、コーナー (壁∩軸はまれだが) でも壁 no-slip を最終確定する。cell/非 node では no-op。
    captureNodeIsothermalEnergyResidual(s.cfg , s.cuda_cfg , s.msh , s.var , "ifaceRro" , "res_ro");   // 診断 (既定 no-op)
    zeroWallDirichletResiduals_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    // node-centered 周期境界 DOF 同一視 (median-dual M4, §4.5): 全 flux/source 積算 + 壁/軸射影の後に、
    // 周期 group の保存量残差を全員で足し合わせ全員へ書き戻す。合併体積と合わせ両側部分 CV を 1 CV として
    // 同期更新する (継ぎ目に双対面を作らず、両側内部双対面が res を組む)。cell/非周期では no-op。
    periodicNodeGather_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
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
    std::list<std::string> names = {"scalarDirichletPin"};
    for (int k = 0; k < s.var.nSpeciesRegistered; ++k) names.push_back("res_roY" + std::to_string(k));
    if (s.var.tracerRegistered != 0) names.push_back("res_roXi");
    for (const auto& nm : s.var.condMomentConsNames) names.push_back("res_" + nm);
    s.var.copyVariables_cell_D2H(names);
    const auto& pin = s.var.c.at("scalarDirichletPin");
    geom_int nPin = 0; for (geom_int ic = 0; ic < s.msh.nCells; ++ic) if (pin[ic] == (flow_float)1.0) ++nPin;
    std::ostringstream os; os << "[pin-diag] step " << s.iStep + 1 << " subiter " << m << " pinned nodes " << nPin << " max|res| on pinned:";
    for (const auto& nm : names) {
        if (nm == "scalarDirichletPin") continue;
        const auto& r = s.var.c.at(nm); double mx = 0.0;
        for (geom_int ic = 0; ic < s.msh.nCells; ++ic) if (pin[ic] == (flow_float)1.0) mx = std::max(mx, (double)std::abs(r[ic]));
        os << " " << nm << " " << std::scientific << std::setprecision(2) << mx;
    }
    std::cout << os.str() << "\n";
}
static void pinRowDiagnosticState(StepContext& s, int m)
{
    if (!pinDiagEnabled() || s.cfg.discretization != "node") return;
    std::list<std::string> names;
    for (int k = 0; k < s.var.nSpeciesRegistered; ++k) names.push_back("Y" + std::to_string(k));
    if (s.var.tracerRegistered != 0) names.push_back("Xi");
    for (const auto& nm : s.var.condMomentConsNames) names.push_back(nm.substr(2));
    if (names.empty()) return;
    s.var.copyVariables_cell_D2H(names);
    std::ostringstream os; os << "[pin-diag] step " << s.iStep + 1 << " subiter " << m << " inlet-node values vs bvar:";
    for (auto& bc : s.msh.bconds) {
        if (bc.bcondKind.rfind("inlet_", 0) != 0 || bc.iPlanes.empty()) continue;
        const geom_int nb = (geom_int)bc.iPlanes.size();
        std::vector<geom_int> cell(nb);
        gpuErrchk( cudaMemcpy(cell.data(), bc.map_bplane_cell_d, nb*sizeof(geom_int), cudaMemcpyDeviceToHost) );
        for (const auto& nm : names) {
            const auto& v = s.var.c.at(nm);
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

// 残差 1 回構築 → 局所擬似時間 dτ → 古典 DPLUR 線形解 → Q への commit。
void implicitNonlinearUpdate(StepContext& s, int inner_index)
{
    assembleResidual(s, 1);
    logResidualSnapshot(s, inner_index);
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
            speciesRenormalize_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            periodicMirrorSpeciesState_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // node 周期: 化学種状態を root→member (§4.1-5)
            speciesPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);     // Y=roY/ρ (出力/次残差用に同期)
        });
    }

    // 液相モーメント (非平衡凝縮) を segregated point-implicit で更新 (Phase 1 はソース=0 の純移流)。
    // res_/transport_diag は assembleResidual の condensationTransport で確定済み。各 wrapper は
    // condensation==0 で no-op。
    s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
        condensationUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // ro*_N = ro*_M = ro*
        condensationTimeIntegration_d_wrapper(0, s.cfg , s.cuda_cfg , s.msh , s.var);
        condensationPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);     // φ=ρφ/ρ (出力/次残差用に同期)
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

    // ソルバ内 CHT: ステップ完了後に壁温を更新する (次ステップの残差組立ての前に効く)。
    // 行番号でなく「完了した定常ステップ数」で interval を数える (advanceImplicitSteady 経路が主対象)。
    conjugateWall::updateConjugateWalls(cfg , msh , var , iStep);

    // detectNaN 診断モード (既定 off): 保存量+P を検査し NaN/Inf があればダンプして停止する。
    if (cfg.detectNaN == 1) {
        checkNonFiniteAndHalt(s);
    }
}

}

int main(void) {
    RuntimeProfiler profiler;

    solverConfig cfg;
    mesh msh;
    matrix mat_ns;
    variables var;
    fluct_variables fluct;
    point_probes pprobes;
    ImplicitDiagLogger implicit_diag_logger;

    cudaConfig cuda_cfg = initializeSimulation(cfg, msh, mat_ns, var, fluct, pprobes);
    // 診断 (FORGE_OUT_RESIDUALS=1): 流れ残差場と陰的補正 dq を h5 出力へ追加する
    // (サブ反復収縮の空間局在の測定用。既定 off = 出力不変)。書かれる値は「最終サブ反復・
    // 最終 sweep 時点」の res_* (BDF 項込み R*) と dq_block_new_* (implicitRelax 適用後)。
    // 注意: blockDPLURSolve は sweep 毎に new/old を swap するため、最終補正は dq_block_old_* に
    // 残る (dq_block_new_* は 1 sweep 前 — 2026-09-03 Codex 指摘で修正)。
    if (const char* e = getenv("FORGE_OUT_RESIDUALS"); e && atoi(e) != 0) {
        for (const char* n : {"res_ro","res_roUx","res_roUy","res_roUz","res_roe",
                              "dq_block_old_0","dq_block_old_1","dq_block_old_2",
                              "dq_block_old_3","dq_block_old_4"})
            var.output_cellValNames.push_back(n);
        printf("[FORGE_OUT_RESIDUALS] residual/dq fields added to h5 outputs\n");
    }
    // FORGE_RESID_SNAP=1: dual-time の subiter 0 直後の res/dq を未使用スロット (res_*_m / dq_*_new
    // スカラー枠) へ退避して出力に含める → 局所収縮率 g=|dq_final|/|dq_sub0| を場で測れる。
    if (const char* e = getenv("FORGE_RESID_SNAP"); e != nullptr) {  // 値は退避 subiter 番号 ("0" も有効)
        for (const char* n : {"res_ro_m","res_roUx_m","res_roUy_m","res_roUz_m","res_roe_m",
                              "dq_ro_new","dq_roUx_new","dq_roUy_new","dq_roUz_new","dq_roe_new"})
            var.output_cellValNames.push_back(n);
        printf("[FORGE_RESID_SNAP] subiter-0 residual/dq snapshots added to h5 outputs\n");
    }
    // 保存量の FP64 影アキュムレータ (plans/active/time_integration-fp64-accumulator.md §4.4)。
    // **非対応の経路で明示 ON されたら黙って劣化させず拒否する** (累積が消える経路があるため)。
    if (cfg.qAccumulatorFP64 == 1) {
        // **v1a の対応範囲** (§5.1 S1a)。どれも原理的な制限ではなく「まだ Qacc を扱っていない」だけで、
        // S1b-①〜④ で順に外す。**黙って劣化させるくらいなら拒否する** (codex plan M5)。
        const char* why = nullptr;
        if (cfg.timeIntegration != 11)
            why = "v1a は timeIntegration=11 のみ (陽解法 tI 1/4 は S1b-④、tI 3 は凸結合なので Qacc_N/Qacc_M が要る)";
        else if (cfg.unsteady != 0)
            why = "v1a は unsteady=0 のみ (dual-time は QaccN/QaccNN の shift が要る。陽解法 unsteady は S1b-④)";
        else if (cfg.isAxisymmetric != 0 && cfg.discretization == "node")
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
            fprintf(stderr, "[qAccumulatorFP64] v1a の対応: timeIntegration=11 && unsteady=0、block/scalar DPLUR、"
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
        msh.buildImplicitLines(var.c.at("ccx").data(), var.c.at("ccy").data(), var.c.at("ccz").data());
    } else if (cfg.lineKFreeze != 0 || cfg.lineViscCoupling != 0 ||
               cfg.lineViscousDtRelief != (flow_float)0.0 || cfg.lineDtDirectional != 0) {
        fprintf(stderr, "[lineImplicit] lineKFreeze/lineViscCoupling/lineViscousDtRelief require lineImplicit=1\n");
        exit(1);
    }
    ResidualCsvLogger residual_logger("residual_history.csv", cfg, msh, var);

    writeInitialOutputs(cfg , msh , var);

    StepMonitor monitor(cfg, residual_logger);
    monitor.printHeader();
    passiveRecordInitialTotals_d_wrapper(cfg, cuda_cfg, msh, var);   // 収支の独立照合の始点 (計算開始前の総量)
    cout << "Start Calculation \n";
    for (int iStep = 0 ; iStep < cfg.mainLoopCount() ; iStep++) {
        advanceOneStep(cfg , cuda_cfg , msh , mat_ns , var , fluct , pprobes , profiler , residual_logger , implicit_diag_logger , iStep);
        monitor.report(iStep);
        // 受動種経路の補正収支 (floor による保存量補正の体積積分; monitorInterval ごと)。scheme 0 / 受動種なしでは no-op。
        if (iStep % cfg.monitorInterval == 0) passiveFloorCorrLog_d_wrapper(cfg, cuda_cfg, msh, var, iStep);
    }
    // 終了時に受動種の収支を必ず出す (最終 step が monitorInterval に乗らないと末尾の補正が記録されない; plan-8 M1)
    if (cfg.mainLoopCount() > 0 && ((cfg.mainLoopCount() - 1) % cfg.monitorInterval) != 0) passiveFloorCorrLog_d_wrapper(cfg, cuda_cfg, msh, var, cfg.mainLoopCount() - 1);

    limiterDiag_finalize(cfg);   // 有界性診断の末尾取りこぼしを回収して累計を確定 (plan §4.35)

    // 壁時計 (旧実装は clock() = CPU 時間で、GPU 待ちを含まなかった)。書式 "Time = %.3f s" は grep 互換のため維持。
    printf("Time = %.3f s (wall, %d steps, %.2f ms/step)\n", monitor.elapsedSeconds(), cfg.mainLoopCount(),
           monitor.elapsedSeconds() * 1.0e3 / std::max(1, cfg.mainLoopCount())); 
    profiler.printSummary();

	return 0;
}