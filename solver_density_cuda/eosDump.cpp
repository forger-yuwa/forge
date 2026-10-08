#include "eosDump.hpp"

// 診断 (FORGE_DUMP_EOS_STEP / FORGE_DUMP_EOS_FILE、既定 off・出力専用)。仕様は eosDump.hpp。
// 書くのは読み出し (D2H) だけで、device の配列・設定には触れない。終了は書き終えた直後 (exit 0)。

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <cuda_runtime.h>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>

#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/thermo_d.cuh"

namespace {

struct EosDumpState {
    bool parsed = false;
    bool on = false;
    int step = -1;           // 1 起点
    std::string path;
    bool beforeDone = false;
    std::unique_ptr<HighFive::File> f;
};
EosDumpState g;

void parseOnce()
{
    if (g.parsed) return;
    g.parsed = true;
    const char* s = std::getenv("FORGE_DUMP_EOS_STEP");
    if (s == nullptr || *s == '\0') return;
    const char* p = std::getenv("FORGE_DUMP_EOS_FILE");
    char* end = nullptr;
    const long k = std::strtol(s, &end, 10);
    if (end == s || *end != '\0' || k < 1 || p == nullptr || *p == '\0') {
        std::fprintf(stderr, "[eos-dump] 拒否: FORGE_DUMP_EOS_STEP は 1 以上の整数、FORGE_DUMP_EOS_FILE は出力先 (step '%s', file '%s')\n",
                     s, p ? p : "(未設定)");
        std::exit(2);
    }
    g.on = true;
    g.step = (int)k;
    g.path = p;
    std::printf("[eos-dump] armed: step %d の EOS の前後の全 cell 配列を %s に書いて終了する\n", g.step, g.path.c_str());
    std::fflush(stdout);
}

// var.c_d の全配列 (std::map なので名前順) を group/<名前> に nCells_all 長の float32 で書く。null の配列は名前だけ記録する。
void writeCellArrays(HighFive::File& f, const std::string& group, mesh& msh, variables& var)
{
    gpuErrchk( cudaDeviceSynchronize() );
    const size_t n = (size_t)msh.nCells_all;
    std::vector<flow_float> h(n);
    std::vector<std::string> nulls;
    for (auto& kv : var.c_d) {
        if (kv.second == nullptr) { nulls.push_back(kv.first); continue; }
        gpuErrchk( cudaMemcpy(h.data(), kv.second, n*sizeof(flow_float), cudaMemcpyDeviceToHost) );
        f.createDataSet("/" + group + "/" + kv.first, h);
    }
    std::string joined;
    for (auto& s : nulls) joined += (joined.empty() ? "" : ";") + s;
    f.getGroup("/" + group).createAttribute<std::string>("null_arrays", HighFive::DataSpace::From(joined)).write(joined);
    f.getGroup("/" + group).createAttribute<long long>("n_arrays", (long long)(var.c_d.size() - nulls.size()));
}

template <class T>
void pack(std::vector<unsigned char>& b, const T& v)
{
    const size_t o = b.size();
    b.resize(o + sizeof(T));
    std::memcpy(b.data() + o, &v, sizeof(T));
}

// 物性 DB (EOS が読む device の係数) を詰め物 (padding) なしでバイト列にする。フィールドの並びは構造体の宣言順。
void writeThermoDb(HighFive::File& f)
{
    const int n = thermo_num_species();
    f.createAttribute<int>("db_n_species", n);
    if (n <= 0) return;
    std::vector<unsigned char> b, bf;
    if (const SpeciesThermo* d = thermo_species_device_ptr()) {
        std::vector<SpeciesThermo> h(n);
        gpuErrchk( cudaMemcpy(h.data(), d, n*sizeof(SpeciesThermo), cudaMemcpyDeviceToHost) );
        for (const auto& s : h) {
            pack(b, s.MW); pack(b, s.sigma_LJ); pack(b, s.eps_kB); pack(b, s.Tlo); pack(b, s.Thi); pack(b, s.nInt);
            for (int k = 0; k < THERMO_MAX_INTERVALS - 1; ++k) pack(b, s.Tbrk[k]);
            for (int k = 0; k < THERMO_MAX_INTERVALS; ++k) for (int c = 0; c < 9; ++c) pack(b, s.coef[k][c]);
            pack(b, s.h_datum); pack(b, s.invMW);
        }
        f.createDataSet("/db/species_thermo", b);
    }
    if (const SpeciesThermoF* d = thermo_species_device_ptr_f()) {
        std::vector<SpeciesThermoF> h(n);
        gpuErrchk( cudaMemcpy(h.data(), d, n*sizeof(SpeciesThermoF), cudaMemcpyDeviceToHost) );
        for (const auto& s : h) {
            pack(bf, s.MW); pack(bf, s.invMW); pack(bf, s.R); pack(bf, s.sigma_LJ); pack(bf, s.eps_kB);
            pack(bf, s.Tlo); pack(bf, s.Thi); pack(bf, s.nInt);
            for (int k = 0; k < THERMO_MAX_INTERVALS - 1; ++k) pack(bf, s.Tbrk[k]);
            for (int k = 0; k < THERMO_MAX_INTERVALS; ++k) for (int c = 0; c < 9; ++c) pack(bf, s.coef[k][c]);
        }
        f.createDataSet("/db/species_thermo_f", bf);
    }
}

// EOS (dependentVariables_d) が読む設定のスカラーと経路。値は格納型のまま書く。
void writeMeta(HighFive::File& f, solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int step)
{
    auto a = [&](const char* name, auto v) { f.createAttribute(name, v); };
    a("step", step);
    a("nCells", (long long)msh.nCells);
    a("nCells_all", (long long)msh.nCells_all);
    a("sizeof_flow_float", (int)sizeof(flow_float));
    a("thermalMethod", cfg.thermalMethod);
    a("gamma", cfg.gamma);
    a("cp", cfg.cp);
    a("pMin", cfg.pMin);
    a("roMin", cfg.roMin);
    a("tMin", cfg.tMin);
    a("nSpecies", cfg.nSpecies);
    a("nSpeciesRegistered", var.nSpeciesRegistered);
    a("thermoFloat", cfg.thermoFloat);
    a("condensation", cfg.condensation);
    a("nCondSpeciesRegistered", var.nCondSpeciesRegistered);
    a("condGasSpecies", cfg.condGasSpecies);
    a("condModel", cfg.condModel);
    a("condEquilibrium", cfg.condEquilibrium);
    a("condSonicModel", cfg.condSonicModel);
    a("condFloat", cfg.condFloat);
    a("isImplicit", cfg.isImplicit);
    a("unsteady", cfg.unsteady);
    a("dualTime", cfg.dualTime);
    a("blocksize", cuda_cfg.blocksize);
    f.createAttribute<std::string>("discretization", HighFive::DataSpace::From(cfg.discretization)).write(cfg.discretization);
    std::ifstream in("solverConfig.yaml", std::ios::binary);
    std::ostringstream ss; ss << in.rdbuf();
    const std::string txt = ss.str();
    f.createAttribute<std::string>("solverConfig_yaml", HighFive::DataSpace::From(txt)).write(txt);
}

}  // namespace

void eosDumpBefore(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int step)
{
    parseOnce();
    if (!g.on || g.beforeDone || step != g.step) return;
    g.f = std::make_unique<HighFive::File>(g.path, HighFive::File::Truncate);
    writeMeta(*g.f, cfg, cuda_cfg, msh, var, step);
    writeThermoDb(*g.f);
    writeCellArrays(*g.f, "pre", msh, var);
    g.beforeDone = true;
}

void eosDumpAfterAndExit(solverConfig& cfg, cudaConfig& cuda_cfg, mesh& msh, variables& var, int step)
{
    (void)cfg; (void)cuda_cfg;
    if (!g.on || !g.beforeDone || step != g.step) return;
    writeCellArrays(*g.f, "post", msh, var);
    g.f->flush();
    g.f.reset();
    std::printf("[eos-dump] wrote %s (step %d, %lld cells incl. ghost, %zu cell arrays x {pre, post}); exiting\n",
                g.path.c_str(), step, (long long)msh.nCells_all, var.c_d.size());
    std::fflush(stdout);
    std::exit(0);
}
