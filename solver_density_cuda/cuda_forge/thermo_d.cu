#include "thermo_d.cuh"
#include "input/solverConfig.hpp"

#include <cuda_runtime.h>
#include <yaml-cpp/yaml.h>

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <cstdlib>

// =============================================================================
// thermo_d.cu
//   化学種熱物性 DB の構築と device へのアップロードを担当する。
//   - 内蔵 DB: 代表的化学種の NASA-9 係数 (CEA / NASA Glenn, McBride 2002) と
//     Lennard-Jones パラメータ (Svehla 1962 等) を保持する。
//   - cfg.speciesDBFile が与えられた場合は yaml で上書き/追加できる。
//   - host 配列と device 配列を 1 本ずつ保持する (run 中は不変なので singleton)。
// =============================================================================

namespace {

// 簡易 CUDA エラーチェック
#define THERMO_CUDA_CHECK(call)                                                  \
    do {                                                                         \
        cudaError_t err__ = (call);                                              \
        if (err__ != cudaSuccess) {                                              \
            std::cerr << "[thermo_d] CUDA error " << cudaGetErrorString(err__)   \
                      << " at " << __FILE__ << ":" << __LINE__ << std::endl;     \
            std::exit(EXIT_FAILURE);                                             \
        }                                                                        \
    } while (0)

std::vector<SpeciesThermo> g_host;   // host 側化学種データ (length = g_n)
SpeciesThermo*             g_dev = nullptr; // device 側コピー
int                        g_n   = 0;

// ---- 内蔵 DB --------------------------------------------------------------
// NASA-9 係数 (a0..a8): cp/R = a0/T^2 + a1/T + a2 + a3 T + a4 T^2 + a5 T^3 + a6 T^4
// MW [kg/mol], sigma_LJ [Angstrom], eps_kB [K]。
// 各 species・各定数の出典 (CEA / Svehla / forge 擬似種 AIR の区別) は
// methods/thermophysics.md 「内蔵 species DB の一覧と出典」を参照。
SpeciesThermo makeSpecies(double MW, double sigma, double eps_kB,
                          double Tlo, double Tmid, double Thi,
                          const double low[9], const double high[9])
{
    SpeciesThermo s;
    s.MW = MW; s.sigma_LJ = sigma; s.eps_kB = eps_kB;
    s.Tlo = Tlo; s.Tmid = Tmid; s.Thi = Thi;
    for (int i=0;i<9;i++){ s.low[i]=low[i]; s.high[i]=high[i]; }
    s.h_datum = 0.0;
    return s;
}

std::map<std::string, SpeciesThermo> builtinDB()
{
    std::map<std::string, SpeciesThermo> db;

    // ---- N2 (CEA) ----
    {
        const double low[9] = { 2.210371497e+04,-3.818461820e+02, 6.082738360e+00,
                               -8.530914410e-03, 1.384646189e-05,-9.625793620e-09,
                                2.519705809e-12, 7.108460860e+02,-1.076003744e+01 };
        const double high[9]= { 5.877124060e+05,-2.239249073e+03, 6.066949220e+00,
                               -6.139685500e-04, 1.491806679e-07,-1.923105485e-11,
                                1.061954386e-15, 1.283210415e+04,-1.586640027e+01 };
        db["N2"] = makeSpecies(0.0280134, 3.621, 97.53, 200.0, 1000.0, 6000.0, low, high);
    }
    // ---- O2 (CEA) ----
    {
        const double low[9] = {-3.425563420e+04, 4.847000970e+02, 1.119010961e+00,
                                4.293889240e-03,-6.836300520e-07,-2.023372700e-09,
                                1.039040018e-12,-3.391454870e+03, 1.849699470e+01 };
        const double high[9]= {-1.037939022e+06, 2.344830282e+03, 1.819732036e+00,
                                1.267847582e-03,-2.188067988e-07, 2.053719572e-11,
                               -8.193467050e-16,-1.689010929e+04, 1.738716506e+01 };
        db["O2"] = makeSpecies(0.0319988, 3.458, 107.4, 200.0, 1000.0, 6000.0, low, high);
    }
    // ---- Ar (単原子, cp/R = 2.5) ----
    {
        const double low[9] = { 0.0,0.0,2.5,0.0,0.0,0.0,0.0,-7.453750000e+02, 4.379674910e+00 };
        const double high[9]= { 0.0,0.0,2.5,0.0,0.0,0.0,0.0,-7.453750000e+02, 4.379674910e+00 };
        db["AR"] = makeSpecies(0.039948, 3.330, 136.5, 200.0, 1000.0, 6000.0, low, high);
        db["Ar"] = db["AR"];
    }
    // ---- CO2 (CEA) ----
    {
        const double low[9] = { 4.943650540e+04,-6.264116010e+02, 5.301725240e+00,
                                2.503813816e-03,-2.127308728e-07,-7.689988780e-10,
                                2.849677801e-13,-4.528198460e+04,-7.048279440e+00 };
        const double high[9]= { 1.176962419e+05,-1.788791477e+03, 8.291523190e+00,
                               -9.223156780e-05, 4.863676880e-09,-1.891053312e-12,
                                6.330036590e-16,-3.908350590e+04,-2.652669281e+01 };
        db["CO2"] = makeSpecies(0.0440095, 3.763, 244.0, 200.0, 1000.0, 6000.0, low, high);
    }
    // ---- He (単原子) ----
    {
        const double low[9] = { 0.0,0.0,2.5,0.0,0.0,0.0,0.0,-7.453750000e+02, 9.287239740e-01 };
        const double high[9]= { 0.0,0.0,2.5,0.0,0.0,0.0,0.0,-7.453750000e+02, 9.287239740e-01 };
        db["HE"] = makeSpecies(0.0040026, 2.551, 10.22, 200.0, 1000.0, 6000.0, low, high);
        db["He"] = db["HE"];
    }
    // ---- H2O (CEA, McBride-Gordon 2002)。非平衡凝縮 (Wyslouzil) のキャリア+凝縮種で使用 ----
    {
        const double low[9] = {-3.947960830e+04, 5.755731020e+02, 9.317826530e-01,
                                7.222712860e-03,-7.342557370e-06, 4.955043490e-09,
                               -1.336933246e-12,-3.303974310e+04, 1.724205775e+01 };
        const double high[9]= { 1.034972096e+06,-2.412698562e+03, 4.646110780e+00,
                                2.291998307e-03,-6.836830480e-07, 9.426468930e-11,
                               -4.822380530e-15,-1.384286509e+04,-7.978148510e+00 };
        db["H2O"] = makeSpecies(0.0180153, 2.605, 572.4, 200.0, 1000.0, 6000.0, low, high);
        db["h2o"] = db["H2O"];
        db["WATER"] = db["H2O"];
    }
    // ---- AIR (単成分擬似空気: 二原子 cp/R≈3.5 + 空気の MW) ----
    // CEA の「Air」混合に近い熱力学を単成分で近似したい場合の簡便種。
    // 低温で cp≈1004.5 J/kgK となるよう定数 cp/R=3.5 (γ=1.4) の擬似多項式。
    {
        const double c[9] = { 0.0,0.0,3.5,0.0,0.0,0.0,0.0,-1.0431373e+03, 3.0 };
        db["AIR"] = makeSpecies(0.0289647, 3.711, 78.6, 200.0, 1000.0, 6000.0, c, c);
        db["Air"] = db["AIR"];
        db["air"] = db["AIR"];
    }

    return db;
}

// yaml ノードから 9 係数を読む
bool read9(const YAML::Node& n, double out[9])
{
    if (!n || !n.IsSequence() || n.size() != 9) return false;
    for (int i=0;i<9;i++) out[i] = n[i].as<double>();
    return true;
}

} // anonymous namespace

void thermo_init_db(solverConfig& cfg)
{
    // calorically-perfect (thermalMethod!=2) でも N=1 のダミーを 1 つ用意しておく
    // (TP 経路以外からは参照されないが、デバイスポインタを null にしないため)。
    std::map<std::string, SpeciesThermo> db = builtinDB();

    // 外部 DB ファイルがあれば内蔵 DB を上書き/追加
    if (!cfg.speciesDBFile.empty()) {
        try {
            YAML::Node root = YAML::LoadFile(cfg.speciesDBFile);
            for (auto it = root.begin(); it != root.end(); ++it) {
                std::string name = it->first.as<std::string>();
                YAML::Node s = it->second;
                SpeciesThermo sp{};
                sp.MW       = s["MW"]       ? s["MW"].as<double>()       : 0.0;
                sp.sigma_LJ = s["LJ_sigma"] ? s["LJ_sigma"].as<double>() : 3.6;
                sp.eps_kB   = s["LJ_eps_kB"]? s["LJ_eps_kB"].as<double>(): 97.0;
                sp.Tlo      = s["Tlo"]      ? s["Tlo"].as<double>()      : 200.0;
                sp.Tmid     = s["Tmid"]     ? s["Tmid"].as<double>()     : 1000.0;
                sp.Thi      = s["Thi"]      ? s["Thi"].as<double>()      : 6000.0;
                double lo[9], hi[9];
                if (read9(s["nasa9_low"], lo) && read9(s["nasa9_high"], hi)) {
                    for (int i=0;i<9;i++){ sp.low[i]=lo[i]; sp.high[i]=hi[i]; }
                    if (sp.MW <= 0.0) { std::cerr << "[thermo_d] species '" << name
                        << "' has invalid MW in DB file" << std::endl; std::exit(EXIT_FAILURE); }
                    db[name] = sp;
                } else {
                    std::cerr << "[thermo_d] species '" << name
                              << "' in DB file lacks valid nasa9_low/nasa9_high (need 9 coeffs)" << std::endl;
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "[thermo_d] failed to read speciesDBFile '" << cfg.speciesDBFile
                      << "': " << e.what() << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    // cfg.speciesNames の順で host 配列を構築
    std::vector<std::string> names = cfg.speciesNames;
    if (names.empty()) names.push_back("N2");  // ダミー (calorically-perfect 用)

    if (static_cast<int>(names.size()) > THERMO_MAX_SPECIES) {
        std::cerr << "[thermo_d] nSpecies=" << names.size()
                  << " exceeds THERMO_MAX_SPECIES=" << THERMO_MAX_SPECIES << std::endl;
        std::exit(EXIT_FAILURE);
    }

    g_host.clear();
    for (const auto& nm : names) {
        auto it = db.find(nm);
        if (it == db.end()) {
            std::cerr << "[thermo_d] species '" << nm
                      << "' not found in built-in DB nor speciesDBFile" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        g_host.push_back(it->second);
    }
    g_n = static_cast<int>(g_host.size());

    // -------------------------------------------------------------------------
    // エンタルピー基準オフセット (sensible-enthalpy datum, thermoHrefTemp>0)
    //   各化学種の絶対 (生成込み) エンタルピー h_s(T) を h_s(Tref) だけ平行移動し
    //   h_s(Tref)=0 に揃える。NASA-9 では h_molar = Ru*T*(...) + Ru*a7 と a7 が
    //   定数項 Ru*a7 を与えるので、両温度域の a[7] に Δa7 = -h_ref/Ru を加算すれば
    //   全 T で一定オフセット c_s=-h_ref を与え Tmid 連続性も保たれる。係数に焼き込む
    //   ため thermo_h_*/cph/T_from_e/混合則/SLAU 面エンタルピー/BC が自動で同一基準に
    //   なり (一点改変で全経路整合)、エントロピー a[8] は不変 (基準はエントロピーに無関係)。
    //   非反応流では支配方程式が不変 (Σh_s J_s* と e_mix の基準移動が連続式で相殺)。
    //   生成エンタルピーの桁違い (H2O≈-13.4MJ/kg) を除くので、多成分 implicit の
    //   roe(block-DPLUR)/roY(point-implicit) 緩和ミスマッチによる Newton 温度ジャンプを抑制。
    if (cfg.thermalMethod == 2 && cfg.thermoHrefTemp > 0.0) {
        const double Tref = cfg.thermoHrefTemp;
        for (int i=0;i<g_n;i++) {
            const double h_ref = thermo_h_molar(g_host[i], Tref);   // 移動前の絶対 h [J/mol]
            const double da7   = -h_ref / THERMO_RU;
            g_host[i].low[7]  += da7;
            g_host[i].high[7] += da7;
            g_host[i].h_datum  = h_ref;   // 反応熱・K_c 用に除いた分を保持 (chemistry_d.cuh)
        }
        std::cout << "[thermo_d] enthalpy datum offset applied: h_s(Tref="
                  << Tref << "K)=0 for all species" << std::endl;
    }

    // device へアップロード
    if (g_dev) { cudaFree(g_dev); g_dev = nullptr; }
    THERMO_CUDA_CHECK(cudaMalloc((void**)&g_dev, g_n*sizeof(SpeciesThermo)));
    THERMO_CUDA_CHECK(cudaMemcpy(g_dev, g_host.data(), g_n*sizeof(SpeciesThermo),
                                 cudaMemcpyHostToDevice));

    std::cout << "[thermo_d] initialized " << g_n << " species:";
    for (const auto& nm : names) std::cout << " " << nm;
    std::cout << std::endl;
    if (cfg.thermalMethod == 2) {
        for (int i=0;i<g_n;i++) {
            std::cout << "  [" << i << "] MW=" << g_host[i].MW
                      << " kg/mol, R=" << (THERMO_RU/g_host[i].MW) << " J/kgK"
                      << ", cp(300K)=" << thermo_cp_mass(g_host[i], 300.0) << " J/kgK"
                      << ", cp(1500K)=" << thermo_cp_mass(g_host[i], 1500.0) << " J/kgK"
                      << ", h(300K)=" << thermo_h_mass(g_host[i], 300.0) << " J/kg"
                      << std::endl;
        }
    }
}

const SpeciesThermo* thermo_species_device_ptr() { return g_dev; }
int                  thermo_num_species()        { return g_n; }
const SpeciesThermo* thermo_species_host()        { return g_host.data(); }
