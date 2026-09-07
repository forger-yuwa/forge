// =============================================================================
// chemistry_d.cu — 有限速度化学ソース項 (Phase 1: 陽ソース + 対角 point-implicit Jacobian)
//   methods/chemistry.md 実装 §2。反応表は chemistry_mech_io.hpp で読み device へ 1 度アップロード。
// =============================================================================
#include <iostream>
#include <vector>
#include <cstring>
#include "chemistry_d.cuh"
#include "chemistry_mech_io.hpp"
#include "chemistrySource_d.cuh"
#include "speciesTransport_d.cuh"

static ReactionTable* g_rt_dev = nullptr;
static ReactionTable  g_rt_host;
static bool g_chem_ready = false;

bool chemistryEnabled(const solverConfig& cfg)
{
    return cfg.chemEnabled != 0 && cfg.thermalMethod == 2 && cfg.nSpecies >= 2;
}

void chemistry_init(solverConfig& cfg)
{
    g_chem_ready = false;
    if (cfg.chemEnabled == 0) return;
    if (cfg.thermalMethod != 2 || cfg.nSpecies < 2) {
        std::cerr << "[chemistry] chemistry.enabled=1 requires thermalMethod=2 and >=2 species" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::vector<std::string> eqs;
    try {
        chem_io::loadMechanism(cfg.chemMechanismFile, cfg.speciesNames, g_rt_host, &eqs);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl; std::exit(EXIT_FAILURE);
    }
    if (g_rt_dev) { cudaFree(g_rt_dev); g_rt_dev = nullptr; }
    gpuErrchk( cudaMalloc((void**)&g_rt_dev, sizeof(ReactionTable)) );
    gpuErrchk( cudaMemcpy(g_rt_dev, &g_rt_host, sizeof(ReactionTable), cudaMemcpyHostToDevice) );
    g_chem_ready = true;

    std::cout << "[chemistry] mechanism '" << cfg.chemMechanismFile << "': " << g_rt_host.nReac
              << " reactions, " << g_rt_host.nSpecies << " species (flow order)."
              << " Tmax=" << cfg.chemTmaxReaction << "K freezeBelowT=" << cfg.chemFreezeBelowT
              << "K jacobianMode=" << cfg.chemJacobianMode << std::endl;
    for (size_t r = 0; r < eqs.size(); ++r) {
        std::cout << "  (" << r+1 << ") " << eqs[r] << "  A=" << g_rt_host.A[r] << " b=" << g_rt_host.b[r]
                  << " Ea=" << g_rt_host.Ea[r] << " J/mol" << (g_rt_host.thirdBody[r] ? " [M]" : "")
                  << (g_rt_host.reversible[r] ? "" : " (irrev)") << "\n";
    }
    // datum の確認 (反応熱は h_datum から出る)
    const SpeciesThermo* sp = thermo_species_host();
    std::cout << "[chemistry] species datum c_s=h_abs(Tref) [J/kg]:";
    for (int s = 0; s < g_rt_host.nSpecies; ++s) std::cout << " " << cfg.speciesNames[s] << "=" << sp[s].h_datum/sp[s].MW;
    std::cout << std::endl;
    if (cfg.thermoHrefTemp <= 0.0) {
        std::cout << "[chemistry] WARNING: thermoHrefTemp=0 (absolute datum). Reaction heat enters through absolute enthalpies;"
                     " block-DPLUR is known to be unstable with absolute datum (methods/thermophysics.md)." << std::endl;
    }
}

const ReactionTable* chemistry_table_host() { return g_chem_ready ? &g_rt_host : nullptr; }

// -----------------------------------------------------------------------------
// セルごとに ω_s, Q̇, 対角 Jacobian を評価し残差へ加算する。
//   res_roY[s] += V ω_s,  src_jac_Y[s] += max(0, −∂ω_s/∂(ρY_s)),  res_roe += V Q̇
//   chemQdot = Q̇ [W/m3], chemTau = 1/max_s|∂ω_s/∂(ρY_s)| [s] (診断)
// -----------------------------------------------------------------------------
__global__ void chemistry_source_d(
    geom_int nCells, int nSpecies, int jacMode,
    double Tmax, double Tfreeze,
    const SpeciesThermo* sp, const ReactionTable* rt,
    const flow_float* vol,      // 周期 node では部分体積 (seam 二重計上防止, bodyForce_d と同じ)
    const flow_float* ro, const flow_float* T,
    flow_float* const* roY_dev,
    flow_float* const* res_roY_dev,
    flow_float* const* src_jac_dev,
    flow_float* res_roe,
    flow_float* chemQdot, flow_float* chemTau)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;

    const double rho = (double)ro[ic];
    double Tc = (double)T[ic];
    chemQdot[ic] = 0.0; chemTau[ic] = 1.0e30;
    if (!(rho > 0.0) || !(Tc > 0.0) || Tc < Tfreeze) return;
    if (Tc > Tmax) Tc = Tmax;
    if (Tc < 200.0) Tc = 200.0;

    double Y[THERMO_MAX_SPECIES], omega[THERMO_MAX_SPECIES], J[THERMO_MAX_SPECIES*THERMO_MAX_SPECIES];
    double ysum = 0.0;
    for (int s = 0; s < nSpecies; ++s) { double y = (double)roY_dev[s][ic] / rho; if (y < 0.0) y = 0.0; Y[s] = y; ysum += y; }
    if (ysum > 0.0) for (int s = 0; s < nSpecies; ++s) Y[s] /= ysum;

    double Qdot = 0.0;
    chem_source(sp, rt, rho, Tc, Y, omega, &Qdot, jacMode, J, nullptr);

    const double V = (double)vol[ic];
    double jmax = 0.0;
    for (int s = 0; s < nSpecies; ++s) {
        res_roY_dev[s][ic] += (flow_float)(V * omega[s]);
        if (jacMode > 0) {
            const double d = J[s*nSpecies+s];
            if (d < 0.0) src_jac_dev[s][ic] += (flow_float)(-d);
            if (fabs(d) > jmax) jmax = fabs(d);
        }
    }
    res_roe[ic] += (flow_float)(V * Qdot);
    chemQdot[ic] = (flow_float)Qdot;
    chemTau[ic]  = (flow_float)((jmax > 0.0) ? 1.0/jmax : 1.0e30);
}

void chemistrySource_d_wrapper(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var)
{
    if (!g_chem_ready || !chemistryEnabled(cfg)) return;
    flow_float** roY = species_roY_device_ptr();
    flow_float** res = species_resroY_device_ptr();
    flow_float** sj  = species_srcjac_device_ptr();
    if (roY == nullptr || res == nullptr || sj == nullptr) return;

    chemistry_source_d<<<cuda_cfg.dimGrid_normalcell, cuda_cfg.dimBlock>>>(
        msh.nCells, var.nSpeciesRegistered, cfg.chemJacobianMode,
        cfg.chemTmaxReaction, cfg.chemFreezeBelowT,
        thermo_species_device_ptr(), g_rt_dev,
        // 周期 node の volume は合併体積なので体積ソースには部分体積を使う (periodicNodeGather で seam 2 倍/コーナー 4 倍になる)。
        (msh.volumePartial_d != nullptr) ? msh.volumePartial_d : var.c_d["volume"],
        var.c_d["ro"], var.c_d["T"],
        roY, res, sj,
        var.c_d["res_roe"], var.c_d["chemQdot"], var.c_d["chemTau"]);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}
