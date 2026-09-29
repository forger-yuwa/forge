#pragma once
// =============================================================================
// cond_latent_test_helper.cuh — 単体試験用: H2O 潜熱の気液ペアを共通データ (埋め込み) から作る
//   (plans/active/thermophysics-solver-owned-species-db.md §4.8, #10)。
//   ソルバは起動時に種 DB の解決結果から作る (condensationTransport_d.cu cond_latent_pair_for)。試験はここで同じ経路
//   (speciesDB_resolve → speciesDB_attachCondensed → cond_latent_pair_make) を通す。気相係数を試験側に書き写さない。
//
// これを include する試験のビルドには、共通データの埋め込みヘッダと種 DB の TU が要る:
//   cmake -DIN=solver_density_cuda/data/species/forge_species_v1.yaml -DOUT=/tmp/forge_species_gen/forge_species_data.hpp
//       -P solver_density_cuda/cmake/embed_species_data.cmake
//   nvcc ... -I solver_density_cuda -I /tmp/forge_species_gen <試験>.cu -x none
//       solver_density_cuda/input/speciesDB.cpp solver_density_cuda/input/speciesTransportDB.cpp -lyaml-cpp
// =============================================================================
#include <string>
#include <vector>
#include <cuda_runtime.h>
#include "input/speciesDB.hpp"
#include "cuda_forge/condensationProperties_d.cuh"

// datum 温度 Tref (0 = datum なし) の気液ペア (H2O + H2O(L))。
inline CondLatentPair cond_test_latent_pair(double Tref)
{
    ResolvedSpeciesDB db = speciesDB_resolve(std::vector<std::string>{"H2O"}, "");
    speciesDB_attachCondensed(db, "H2O(L)", "H2O", true);
    const ResolvedCondensed& c = db.condensed;
    return cond_latent_pair_make(c.gas, c.coeffs, c.Tlo, c.Thi, Tref);
}

// 参照 (host は常に、device は with_device のときだけ複製)。ペアはプロセス終了まで保持する (試験用)。
inline CondLatentRef cond_test_latent_ref(double Tref, bool with_device)
{
    CondLatentPair* h = new CondLatentPair(cond_test_latent_pair(Tref));
    CondLatentRef r;
    r.h = h;
#ifdef __CUDACC__
    if (with_device) {
        CondLatentPair* d = nullptr;
        cudaMalloc((void**)&d, sizeof(CondLatentPair));
        cudaMemcpy(d, h, sizeof(CondLatentPair), cudaMemcpyHostToDevice);
        r.d = d;
    }
#else
    (void)with_device;   // 純 host ビルド (g++) では device 複製を作らない
#endif
    return r;
}

// condProps_H2O() に気液ペアを入れたもの (ソルバでは condProps_make(COND_MODEL_H2O, cond_prop_opts(cfg)) に相当)。
inline CondSpeciesProps cond_test_props_H2O(bool with_device, double Tref = 0.0)
{
    CondSpeciesProps s = condProps_H2O();
    s.lat = cond_test_latent_ref(Tref, with_device);
    return s;
}
