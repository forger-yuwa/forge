#include "thermo_d.cuh"
#include "transportMix_d.cuh"
#include "transportTables_d.cuh"
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

// 種ごとの輸送物性 (physProp.transport; plan thermophysics-solver-owned-species-db #5t2-2)。
//   g_trans は device ポインタを持つ host 側の表 (カーネルへ値渡し)。書かれていない run では g_transOn = false。
SpeciesTransportD* g_trans_sp     = nullptr;
TransportPairD*    g_trans_pairs  = nullptr;
double*            g_trans_expand = nullptr;
TransportTableD    g_trans        = {0, 0, nullptr, nullptr, nullptr};
bool               g_transOn      = false;
// 表引き (#5t2-3)。host 側の表 (記録・試験用に保持) と device 側の配列。
TransportTablesHost g_tabHost;
std::vector<void*>  g_tabDev;

template <class T>
const T* thermo_upload_vec(const std::vector<T>& v)
{
    if (v.empty()) return nullptr;
    void* d = nullptr;
    THERMO_CUDA_CHECK(cudaMalloc(&d, v.size()*sizeof(T)));
    THERMO_CUDA_CHECK(cudaMemcpy(d, v.data(), v.size()*sizeof(T), cudaMemcpyHostToDevice));
    g_tabDev.push_back(d);
    return static_cast<const T*>(d);
}

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
    // 区間数 (区間可変 #13-1)。speciesDB が全入口で検査済みだが、カーネルの固定長配列を越えないことをここでも確かめる。
    for (int i = 0; i < g_n; ++i) {
        if (g_host[i].nInt < 1 || g_host[i].nInt > THERMO_MAX_INTERVALS) {
            std::cerr << "[thermo_d] species '" << names[i] << "' has " << g_host[i].nInt << " temperature intervals (supported 1 to "
                      << THERMO_MAX_INTERVALS << ")" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    // -------------------------------------------------------------------------
    // エンタルピー基準オフセット (sensible-enthalpy datum, thermoHrefTemp>0)
    //   各化学種の絶対 (生成込み) エンタルピー h_s(T) を h_s(Tref) だけ平行移動し
    //   h_s(Tref)=0 に揃える。NASA-9 では h_molar = Ru*T*(...) + Ru*a7 と a7 が
    //   定数項 Ru*a7 を与えるので、全温度区間の a[7] に Δa7 = -h_ref/Ru を加算すれば
    //   全 T で一定オフセット c_s=-h_ref を与え区切りでの連続性 (段差) も変わらない。係数に焼き込む
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
            thermo_add_a7(g_host[i], da7);   // 全区間の a7 に同じ Δa7 (区切りの段差は不変)
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
        for (int i=0;i<g_n;i++) hf[i] = thermo_to_float(g_host[i]);
        if (g_dev_f) { cudaFree(g_dev_f); g_dev_f = nullptr; }
        THERMO_CUDA_CHECK(cudaMalloc((void**)&g_dev_f, g_n*sizeof(SpeciesThermoF)));
        THERMO_CUDA_CHECK(cudaMemcpy(g_dev_f, hf.data(), g_n*sizeof(SpeciesThermoF), cudaMemcpyHostToDevice));
    }

    // 種ごとの輸送物性 (physProp.transport があるときだけ; #5t2-2)。実種表・組の表・展開行列 (密, 行 = 輸送種) を上げる。
    //   上限 (実種 TRANSPORT_MAX_REAL_SPECIES・輸送種 THERMO_MAX_SPECIES) は resolver でも検査済みだが、
    //   カーネルの固定長配列を越えないことをここでもう一度確かめる。
    if (g_trans_sp)     { cudaFree(g_trans_sp);     g_trans_sp = nullptr; }
    if (g_trans_pairs)  { cudaFree(g_trans_pairs);  g_trans_pairs = nullptr; }
    if (g_trans_expand) { cudaFree(g_trans_expand); g_trans_expand = nullptr; }
    for (void* p : g_tabDev) cudaFree(p);
    g_tabDev.clear();
    g_tabHost = TransportTablesHost{};
    g_trans   = TransportTableD{0, 0, nullptr, nullptr, nullptr};
    g_transOn = false;
    if (db.transport.enabled) {
        const ResolvedTransport& tr = db.transport;
        const int nR = tr.nReal();
        const int nT = static_cast<int>(tr.expand.size());
        if (nR < 1 || nR > TRANSPORT_MAX_REAL_SPECIES || nT != g_n || nT < 1 || nT > THERMO_MAX_SPECIES
            || static_cast<int>(tr.pairs.size()) != nR*(nR - 1)/2) {
            std::cerr << "[thermo_d] transport table out of range: nReal=" << nR << " (max " << TRANSPORT_MAX_REAL_SPECIES
                      << "), transported=" << nT << " (species " << g_n << ", max " << THERMO_MAX_SPECIES
                      << "), pairs=" << tr.pairs.size() << std::endl;
            std::exit(EXIT_FAILURE);
        }
        std::vector<double> E(static_cast<size_t>(nT)*nR, 0.0);
        for (int s = 0; s < nT; ++s)
            for (const auto& e : tr.expand[s]) E[static_cast<size_t>(s)*nR + e.first] += e.second;
        THERMO_CUDA_CHECK(cudaMalloc((void**)&g_trans_sp, nR*sizeof(SpeciesTransportD)));
        THERMO_CUDA_CHECK(cudaMemcpy(g_trans_sp, tr.sp.data(), nR*sizeof(SpeciesTransportD), cudaMemcpyHostToDevice));
        if (!tr.pairs.empty()) {
            THERMO_CUDA_CHECK(cudaMalloc((void**)&g_trans_pairs, tr.pairs.size()*sizeof(TransportPairD)));
            THERMO_CUDA_CHECK(cudaMemcpy(g_trans_pairs, tr.pairs.data(), tr.pairs.size()*sizeof(TransportPairD), cudaMemcpyHostToDevice));
        }
        THERMO_CUDA_CHECK(cudaMalloc((void**)&g_trans_expand, E.size()*sizeof(double)));
        THERMO_CUDA_CHECK(cudaMemcpy(g_trans_expand, E.data(), E.size()*sizeof(double), cudaMemcpyHostToDevice));
        g_trans   = TransportTableD{nR, nT, g_trans_sp, g_trans_pairs, g_trans_expand};
        g_transOn = true;
        std::cout << "[thermo_d] transport (physProp.transport) uploaded: " << nR << " real species, "
                  << tr.pairs.size() << " pairs, expand " << nT << "x" << nR << " (mole basis)" << std::endl;

        // 表引き (#5t2-3): 既定で使う。FORGE_TRANSPORT_TABLE=0 のときだけ段 2 の double 評価のまま (性能比較・切り分け用)。
        const char* tabEnv = std::getenv("FORGE_TRANSPORT_TABLE");
        if (tabEnv != nullptr && std::string(tabEnv) == "0") {
            std::cout << "[thermo_d] transport tables disabled (FORGE_TRANSPORT_TABLE=0): double evaluation" << std::endl;
        } else {
            std::vector<double> MWs(nT);
            for (int s = 0; s < nT; ++s) MWs[s] = g_host[s].MW;
            if (!transport_tables_build_host(tr.sp, tr.pairs, MWs, E, g_tabHost)) {
                std::cerr << "[thermo_d] ERROR: transport table build failed: " << g_tabHost.error << std::endl;
                std::exit(EXIT_FAILURE);
            }
            TransportTablesF& t = g_trans.tab;
            t.valid   = 1;
            t.Tmin    = g_tabHost.Tmin;
            t.Tmax    = g_tabHost.Tmax;
            t.seg     = thermo_upload_vec(g_tabHost.seg);
            t.spc     = thermo_upload_vec(g_tabHost.spc);
            t.pairc   = thermo_upload_vec(g_tabHost.pairc);
            t.spTab   = thermo_upload_vec(g_tabHost.spTab);
            t.pairTab = thermo_upload_vec(g_tabHost.pairTab);
            t.pm      = thermo_upload_vec(g_tabHost.pm);
            t.invMWs  = thermo_upload_vec(g_tabHost.invMWs);
            t.expandF = thermo_upload_vec(g_tabHost.expandF);
            int nTabPairs = 0;
            for (const auto& r : g_tabHost.pairTab) nTabPairs += (r.nseg > 0) ? 1 : 0;
            std::cout << "[thermo_d] transport tables (float, ln T Hermite): T in [" << t.Tmin << ", " << t.Tmax
                      << "] K (outside: double), " << nR << " species + " << nTabPairs << " tabulated pairs ("
                      << (tr.pairs.size() - nTabPairs) << " rigid-sphere at run time), " << g_tabHost.seg.size()
                      << " segments, " << (g_tabHost.spc.size()/2 + g_tabHost.pairc.size()) << " sub-intervals, "
                      << g_tabHost.bytes() << " bytes" << std::endl;
        }
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
const TransportTableD* thermo_transport_table()     { return g_transOn ? &g_trans : nullptr; }
const TransportTablesHost* thermo_transport_tables_host() { return (g_transOn && g_trans.tab.valid) ? &g_tabHost : nullptr; }
