#include "thermo_d.cuh"
#include "input/solverConfig.hpp"
#include "input/speciesDB.hpp"

#include <cuda_runtime.h>

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
SpeciesThermo*             g_dev = nullptr;
SpeciesThermoF* g_dev_f = nullptr; // device 側コピー
int                        g_n   = 0;

// 内蔵 DB・yaml 上書き・名前解決は host 側 input/speciesDB.cpp (speciesDB_resolve) に集約した
// (convertGmshToForge も GPU 無しで同じ解決を使う)。ここは device アップロードと datum オフセットだけ。

} // anonymous namespace

void thermo_init_db(solverConfig& cfg)
{
    // calorically-perfect (thermalMethod!=2) でも N=1 のダミーを 1 つ用意しておく
    // (TP 経路以外からは参照されないが、デバイスポインタを null にしないため)。
    // 名前→係数の解決は host 関数 speciesDB_resolve (cfg.read() 直後に speciesDB_init 済みならそれを再利用)。
    const ResolvedSpeciesDB* pre = speciesDB_current();
    const ResolvedSpeciesDB  db  = pre ? *pre : speciesDB_init(cfg);   // 未知種名は speciesDB_init が exit
    const std::vector<std::string>& names = db.names;

    g_host = db.species;
    g_n = static_cast<int>(g_host.size());
    for (auto& s : g_host) s.invMW = 1.0/s.MW;   // 研磨段の乗算用 (yaml 由来の種も含め全種)

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
    // float32 ミラー (面ループ用)。datum オフセット焼き込み後の係数をそのまま float へ。
    {
        std::vector<SpeciesThermoF> hf(g_n);
        for (int i=0;i<g_n;i++) {
            const SpeciesThermo& s = g_host[i];
            SpeciesThermoF& f = hf[i];
            f.MW = (float)s.MW; f.invMW = (float)(1.0/s.MW); f.R = (float)(THERMO_RU/s.MW);
            f.sigma_LJ = (float)s.sigma_LJ; f.eps_kB = (float)s.eps_kB;
            f.Tlo = (float)s.Tlo; f.Tmid = (float)s.Tmid; f.Thi = (float)s.Thi;
            for (int k=0;k<9;k++) { f.low[k] = (float)s.low[k]; f.high[k] = (float)s.high[k]; }
        }
        if (g_dev_f) { cudaFree(g_dev_f); g_dev_f = nullptr; }
        THERMO_CUDA_CHECK(cudaMalloc((void**)&g_dev_f, g_n*sizeof(SpeciesThermoF)));
        THERMO_CUDA_CHECK(cudaMemcpy(g_dev_f, hf.data(), g_n*sizeof(SpeciesThermoF), cudaMemcpyHostToDevice));
    }

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
const SpeciesThermoF* thermo_species_device_ptr_f() { return g_dev_f; }
int                  thermo_num_species()        { return g_n; }
const SpeciesThermo* thermo_species_host()        { return g_host.data(); }
