#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <stdio.h>                                                                                       
#include <fstream>
#include <string>
#include <time.h>
#include <limits>
#include <cstdio>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

#include "input/solverConfig.hpp"
#include "input/setInitial.hpp"

#include "mesh/gmshReader.hpp"
#include "output/output.hpp"
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
#include "cuda_forge/speciesTransport_d.cuh"
#include "cuda_forge/chemistrySource_d.cuh"
#include "cuda_forge/condensationTransport_d.cuh"
#include "cuda_forge/viscousFlux_d.cuh"
#include "cuda_forge/updateCenterVelocity_d.cuh"
#include "cuda_forge/interpVelocity_c2p_d.cuh"
#include "cuda_forge/timeIntegration_d.cuh"
#include "cuda_forge/limiter_d.cuh"
#include "cuda_forge/ducrosSensor_d.cuh"
#include "cuda_forge/axisymmetricSource_d.cuh"
#include "cuda_forge/wmlesWallModel_d.cuh"
#include "cuda_forge/nodeWallDirichlet_d.cuh"
#include "cuda_forge/bodyForce_d.cuh"
#include "cuda_forge/turbulent_viscosity_d.cuh"
#include "cuda_forge/residualMonitor_d.cuh"

#include "cuda_forge/fluct_variables_d.cuh"
#include "cuda_forge/gasProperties_d.cuh"
#include "cuda_forge/thermo_d.cuh"

#include "probe/point_probes.cuh"
#include "cuda_forge/setDT_d.cuh"
#include "cuda_forge/periodicNode_d.cuh"

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

    // buffer に残った残差を強制 flush する (停止前などに使用)。
    void flush()
    {
        if (gpu_) {
            flushDevice();
        }
    }

private:
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

    // 入口分布プロファイル (ints:{inletProfile:1} の inlet): per-face bvar を CSV から face 重心で補間。
    // mesh.planes (重心) と bvar が揃った後・最初の applyBconds より前に適用。未指定 inlet は一様のまま。
    applyInletProfiles(cfg , msh);

    // 化学種変数を登録 (allocVariables より前)。nSpecies<=1 では no-op。
    var.registerSpecies(cfg.nSpecies, cfg.chemEnabled);

    // 非平衡凝縮モーメント変数を登録 (allocVariables より前)。condensation==0 では no-op。
    var.registerCondensation(cfg.nCondSpecies);

    var.allocVariables(cfg.gpu , msh);

    // device roY[] ポインタ配列を構築 (c_d 確保後, dependentVariables より前)。
    speciesInit_d(cfg , var);

    // device rog[] ポインタ配列を構築 (二相 EOS が液相質量分率を読む)。condensation==0 で no-op。
    condensationInit_d(cfg , var);

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
    dependentVariables(cfg , cuda_cfg , msh , var, mat_ns);
    // node-centered 壁 Dirichlet: IC の壁ノード速度を厳密 0 に初期化 (KE を roe から除去)。
    // この後 gasProperties が補正 roe から P/T を再計算する。cell/非 node では no-op。
    enforceWallNoSlip_d_wrapper(cfg , cuda_cfg , msh , var);
    gasProperties_d_wrapper(cfg , cuda_cfg , msh , var);

    fluct.allocVariables();
    fluct.set_fluctVelocity(cfg , cuda_cfg , msh , var);

    applyBconds(cfg , cuda_cfg , msh , var, mat_ns , fluct);
    applyRansScalarBoundaries(cfg , cuda_cfg , msh , var);
    applyWmlesWallModel(cfg , cuda_cfg , msh , var);   // WMLES 壁応力モデル (wallModelLES 壁のみ, §10)
    applyNodeIsothermalWallPin(cfg , cuda_cfg , msh , var);   // 素の node 等温壁の壁ノード T ピン
    applySstThermalWallFunction(cfg , cuda_cfg , msh , var);  // SST 熱的壁関数: 断熱壁 T_aw (§6.5(f))
    applySpeciesBoundaries(cfg , cuda_cfg , msh , var);
    applyCondensationBoundaries(cfg , cuda_cfg , msh , var);
    calcGradient_d_wrapper(cfg , cuda_cfg , msh , var);
    // 初期 setup でも周期勾配 gather を適用 (assembleResidual と整合; res_0 出力と初期診断を正しい合併勾配にする)。
    periodicGradientGather_d_wrapper(cfg , cuda_cfg , msh , var);
    axisymmetricGeomTerms_d_wrapper(cfg , cuda_cfg , msh , var);
    updateVariablesOuter(cfg , cuda_cfg , msh , var , mat_ns);
    speciesUpdateOuter_d_wrapper(cfg , cuda_cfg , msh , var);  // roY{s}N/M ベースライン
    condensationUpdateOuter_d_wrapper(cfg , cuda_cfg , msh , var);  // 液相モーメント N/M ベースライン
    setDT_d_wrapper(cfg , cuda_cfg , msh , var);

    pprobes.init(cfg , cuda_cfg , msh);

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
    s.profiler.measureCuda(ProfileSection::DependentVariables, [&]() {
        speciesPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // Y_s = ρY_s/ρ (混合則 thermo の前)
        condensationPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // φ = ρφ/ρ (スカラ移流の上流値)
    });
    s.profiler.measureWall(ProfileSection::DependentVariables, [&]() {
        dependentVariables(s.cfg , s.cuda_cfg , s.msh , s.var, s.mat_ns);
    });
    s.profiler.measureCuda(ProfileSection::GasProperties, [&]() {
        gasProperties_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureWall(ProfileSection::ApplyBconds, [&]() {
        applyBconds(s.cfg , s.cuda_cfg , s.msh , s.var, s.mat_ns , s.fluct);
    });
    s.profiler.measureWall(ProfileSection::ApplyBconds, [&]() {
        applyRansScalarBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);
        applyWmlesWallModel(s.cfg , s.cuda_cfg , s.msh , s.var);   // WMLES 壁応力モデル (§10)
        applyNodeIsothermalWallPin(s.cfg , s.cuda_cfg , s.msh , s.var);   // 素の node 等温壁 T ピン
        applySstThermalWallFunction(s.cfg , s.cuda_cfg , s.msh , s.var);  // SST 熱的壁関数: 断熱壁 T_aw (§6.5(f))
        applySpeciesBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);
        applyCondensationBoundaries(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::CalcGradient, [&]() {
        calcGradient_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        // 多成分 face 整合再構成 (speciesFaceReconstruction==1): ∇Y_s を Green-Gauss で計算。
        // species ghost は直前の applySpeciesBoundaries で Neumann 充填済み。既定 0 で no-op。
        // node + 粘性多成分でも ∇Y が要る: species_diffusion_d の境界半割面 ghostless 弱形式
        // (J_s=ρD∇Y·S) が ∇Y を参照するため、その場合も計算しておく。
        if (s.cfg.speciesFaceReconstruction >= 1 ||
            (s.cfg.discretization == "node" && s.cfg.viscMethod != 0)) {
            speciesGradient_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
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
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        ransTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        speciesTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);  // 化学種移流残差
        chemistrySource_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);   // 有限速度化学ソース (ω_s, Q̇, 対角 Jacobian)
    });
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        condensationTransport_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);  // 液相モーメント移流残差 (Phase 1)
        condensationSource_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var);     // 核生成+成長ソース (Phase 2)
    });
    s.profiler.measureCuda(ProfileSection::TurbulenceModel, [&]() {
        ransGradient_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        ransSource_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });
    s.profiler.measureCuda(ProfileSection::AxisymmetricSource, [&]() {
        axisymmetricSource_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);      // method 0 (r 重み): hoop 源
        axisymmetricSourceSU2_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // method 1 (SU2 流): 1/y 全ソース
        bodyForce_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);   // 一様体積力 (bodyForce, off なら no-op)
    });
    s.profiler.measureCuda(ProfileSection::ViscousFlux, [&]() {
        viscousFlux_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var, s.mat_ns);
    });
    // SU2 流の軸対称対称面 (MARKER_SYM): 軸上 CV の半径方向運動量残差を 0 に射影し roUy=0 を保つ。
    // explicit では ΔroUy=0 になり直接効く (implicit/block-DPLUR では連成 solve が補正を漏らすため Jacobian
    // 整合が別途必要・open issue, docs §7.1)。cell/非軸対称/平面では no-op。
    zeroAxisRadialResidual_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    // node-centered 壁 Dirichlet: 壁ノードの運動量残差を 0 に射影し u=0 を保つ (壁ゴースト撤廃の代替)。
    // 軸射影の後に置き、コーナー (壁∩軸はまれだが) でも壁 no-slip を最終確定する。cell/非 node では no-op。
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
    const bool printCfl = onMonitor;
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
            applySSTPointImplicit(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
            // node 周期 DOF 同一視 (§4.5): point-implicit SST 更新後に k/ω 状態を root→member ミラーし drift を防ぐ。
            periodicMirrorScalarState_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
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
        });

        cout << "       " << iteration_label << " : " << iloop+1 << "\n";
        assembleResidual(s, iloop + 1);
        logResidualSnapshot(s, iloop);
        s.profiler.measureCuda(ProfileSection::TimeIntegration, [&]() {
            timeIntegration_d_wrapper(iloop, s.cfg , s.cuda_cfg , s.msh , s.var);
            // node 周期境界 DOF 同一視 (§4.5.9): NS 保存量更新直後に root→member ミラーで slave=master を強制。
            // 残差 gather だけでは初期 desync (非周期 seed 摂動) が残り継ぎ目フラックス不整合を生むため。cell/非周期で no-op。
            periodicMirrorNSState_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            ransTimeIntegration_d_wrapper(iloop, s.cfg , s.cuda_cfg , s.msh , s.var);
            speciesTimeIntegration_d_wrapper(iloop, s.cfg , s.cuda_cfg , s.msh , s.var);
            speciesRenormalize_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);  // ρY_s>=0, ΣρY_s=ρ
            condensationTimeIntegration_d_wrapper(iloop, s.cfg , s.cuda_cfg , s.msh , s.var);  // 液相モーメント (Phase 1 ソース=0)
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
    });
    s.profiler.measureWall(ProfileSection::WriteOutputs, [&]() {
        writeStepOutputs(s.cfg , s.cuda_cfg , s.msh , s.var , s.pprobes , s.iStep+1);
    });
    // explicit (陽解法) は cfg.dt を時間前進に使うため dt 適応は毎ステップ行う。表示のみ monitorInterval で間引く。
    const bool printCflExp = (s.iStep % s.cfg.monitorInterval == 0);
    s.profiler.measureCuda(ProfileSection::SetDt, [&]() {
        setDT_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var, /*adaptDt=*/true, /*printCfl=*/printCflExp);
    });
    s.residual_logger.logOuterEnd(s.iStep);
    if (s.cfg.unsteady == 1) {
        s.cfg.totalTime += s.cfg.dt;
    }
}

// 定常 block DPLUR 陰解法。メインループ（nStepOuter）を擬似時間とし、1 ステップ = 1 非線形更新の縮退形。
void advanceImplicitSteady(StepContext& s)
{
    // baseline (roN) は前ステップ末尾 / 初期化の updateVariablesOuter で設定済み（ro == roN）。
    implicitNonlinearUpdate(s, 0);

    s.profiler.measureWall(ProfileSection::UpdateOuter, [&]() {
        updateVariablesOuter(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
    });
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

    // 物理時間レベルシフト: roNN ← roN, roN ← ro（現在の ro = Q^n）。
    s.profiler.measureWall(ProfileSection::UpdateOuter, [&]() {
        shiftDualTimeLevels_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
    });

    // BDF 係数: 初回ステップ or bdfOrder==1 は BDF1 (1,1,0)、以降 BDF2 (3/2,2,1/2)。
    const bool useBDF2 = (s.cfg.bdfOrder >= 2) && (s.iStep > 0);
    const flow_float a = useBDF2 ? static_cast<flow_float>(1.5) : static_cast<flow_float>(1.0);
    const flow_float b = useBDF2 ? static_cast<flow_float>(2.0) : static_cast<flow_float>(1.0);
    const flow_float c = useBDF2 ? static_cast<flow_float>(0.5) : static_cast<flow_float>(0.0);
    const int include_scalar = scalarResidualEnabled(s.cfg) ? 1 : 0;
    // 対角へ加える物理時間項係数 a/Δt（block/scalar/SST カーネルが cfg 経由で参照）。
    s.cfg.unsteadyDiagCoef = a / std::max(s.cfg.dt, static_cast<flow_float>(1.0e-30));

    const int nSub = std::max(1, s.cfg.nSubIterDualTime);
    for (int m = 0; m < nSub; ++m) {
        assembleResidual(s, 1);
        // 残差に物理時間 BDF 項を加える: res* = res - (V/Δt)(a Q - b Q^n + c Q^{n-1})。
        addUnsteadyTimeTerm_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var, a, b, c, include_scalar);
        logResidualSnapshot(s, m);
        // dual-time も pseudo/physical の時間を CFL に基づき変えうるため、dt 適応は毎サブ反復で行う
        // (adaptDt=true → dtControl==1 のとき cfg.dt を適応; host 読み出しもそのとき発生)。表示のみ
        // monitorInterval で間引く。
        const bool printCflDt = (s.iStep % s.cfg.monitorInterval == 0);
        s.profiler.measureCuda(ProfileSection::SetDt, [&]() {
            setDT_d_wrapper(s.cfg , s.cuda_cfg, s.msh , s.var, /*adaptDt=*/true, /*printCfl=*/printCflDt);
        });
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
        // dual-time の commit は in-place（roN=Q^n は BDF 基準で固定のため roN+dq は使えない）。
        s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
            applyBlockImplicitCorrectionInPlace_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        });
        // FORGE_FREEZE_TURB=1 は dual-time でも SST 状態更新を凍結する (2026-09-03 修正: 従来この
        // 経路は無条件更新で、freeze 診断が定常専用だった — dual-time A/B は無効だった)。
        if (include_scalar && !freezeTurbEnabled()) {
            s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
                applySSTPointImplicit(s.cfg , s.cuda_cfg , s.msh , s.var , s.mat_ns);
                periodicMirrorScalarState_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var); // §4.5 k/ω 周期ミラー
            });
        } else if (include_scalar) {
            static bool logged = false;
            if (!logged) { printf("[FREEZE_TURB] dual-time: SST state update frozen\n"); logged = true; }
        }
        // 液相モーメント (非平衡凝縮) を segregated point-implicit で更新。condensation==0 で no-op。
        s.profiler.measureWall(ProfileSection::UpdateInner, [&]() {
            condensationUpdateOuter_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
            condensationTimeIntegration_d_wrapper(0, s.cfg , s.cuda_cfg , s.msh , s.var);
            condensationPrimitive_d_wrapper(s.cfg , s.cuda_cfg , s.msh , s.var);
        });
    }

    s.cfg.unsteadyDiagCoef = 0.0; // 定常側へ影響しないようリセット
    s.cfg.totalTime += s.cfg.dt;

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
    cout << "----------------------------\n";
    cout << "Step : " << iStep;
    if (cfg.unsteady == 1) {
        cout << "  Time : " << cfg.totalTime;
    }
    cout << "\n";

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

    // detectNaN 診断モード (既定 off): 保存量+P を検査し NaN/Inf があればダンプして停止する。
    if (cfg.detectNaN == 1) {
        checkNonFiniteAndHalt(s);
    }
}

}

int main(void) {
    clock_t start = clock();
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
    // line-implicit (plans/active/time_integration-line-implicit.md): 壁法線ラインを構築。
    // blockDPLUR==1 専用・完全前処理 (lowMachPrecond>=2) とは併用不可。
    if (cfg.lineImplicit == 1) {
        if (cfg.blockDPLUR != 1 || cfg.lowMachPrecond >= 2 || cfg.timeIntegration != 11) {
            fprintf(stderr, "[lineImplicit] requires timeIntegration=11, blockDPLUR=1, lowMachPrecond<2\n");
            exit(1);
        }
        msh.buildImplicitLines(var.c.at("ccx").data(), var.c.at("ccy").data(), var.c.at("ccz").data());
        // lineDtWallRelief 診断用: wall 種 bcond の境界面フラグ (kind に "wall" を含む bcond のみ。
        // inlet/outlet/periodic は除外しない — 壁境界 Λ の律速仮説の切り分け専用)。
        if (cfg.lineDtWallRelief == 1) {
            std::vector<unsigned char> pw(msh.nPlanes, 0);
            geom_int nw = 0;
            for (auto& bc : msh.bconds) {
                if (bc.bcondKind.find("wall") == std::string::npos) continue;
                for (auto ip : bc.iPlanes) { pw[ip] = 1; ++nw; }
            }
            gpuErrchk(cudaMalloc((void**)&msh.plane_wall_flag_d, sizeof(unsigned char)*msh.nPlanes));
            gpuErrchk(cudaMemcpy(msh.plane_wall_flag_d, pw.data(), sizeof(unsigned char)*msh.nPlanes, cudaMemcpyHostToDevice));
            printf("[lineDtWallRelief] wall boundary planes flagged: %d\n", (int)nw);
        }
    } else if (cfg.lineKFreeze != 0 || cfg.lineViscCoupling != 0 ||
               cfg.lineViscousDtRelief != (flow_float)0.0 || cfg.lineDtDirectional != 0 ||
               cfg.lineDtWallRelief != 0) {
        fprintf(stderr, "[lineImplicit] lineKFreeze/lineViscCoupling/lineViscousDtRelief require lineImplicit=1\n");
        exit(1);
    }
    ResidualCsvLogger residual_logger("residual_history.csv", cfg, msh, var);

    writeInitialOutputs(cfg , msh , var);

    cout << "Start Calculation \n";
    for (int iStep = 0 ; iStep < cfg.mainLoopCount() ; iStep++) {
        advanceOneStep(cfg , cuda_cfg , msh , mat_ns , var , fluct , pprobes , profiler , residual_logger , implicit_diag_logger , iStep);
    }

    clock_t end = clock();
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time = %.3f s\n", time); 
    profiler.printSummary();

	return 0;
}