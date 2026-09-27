#pragma once

// =============================================================================
// speciesLump.hpp
//   physProp.species の lump (擬似種) 指定 (plans/active/thermophysics-solver-owned-species-db.md §4.2, #6a;
//   仕様 methods/thermophysics.md §1b.2)。solverConfig (読み込み) と speciesDB (起動時合成) の両方が使うので
//   依存の無い小さなヘッダに分けた (solverConfig.hpp から thermo_d.cuh を引き込まないため)。
//
//   physProp:
//     species:
//       - {name: MIXDRY, lump: {N2: 0.708873, O2: 0.230376, AR: 0.00850387, CO2: 0.0522474}, basis: mole}
//       - H2O
// =============================================================================

#include <string>
#include <vector>

struct SpeciesLumpSpec {
    std::string              name;       // 輸送する擬似種の名前 (physProp.species の位置 = index s)
    std::string              basis;      // "mole" | "mass" (members の分率の基準)
    std::vector<std::string> members;    // 構成種名 (config に書いた順・書いた綴り)
    std::vector<double>      fractions;  // 同順。config に書いた値そのまま (正規化・検査は speciesDB_resolve)
};
