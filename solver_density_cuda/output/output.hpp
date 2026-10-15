#pragma once

#include "mesh/mesh.hpp"
#include "input/solverConfig.hpp"
#include "variables.hpp"

//void outputH5_XDMF(const mesh &msh , const variables &var)
// 全場 res_<step>.h5。間隔は time.fieldOutStepInterval (省略時 = outStepInterval)、最後の step は必ず書く。
// force = true (協調停止の step) は間隔に関係なく書く (plans/active/tooling-nozzle-card-local-mesh.md §4.4 1′)。
void outputH5_XDMF(const solverConfig& , const mesh& , variables& , const int& , bool force = false);

// detectNaN 診断モード用: 出力間隔ガードを無視して現在の解を prefix 付きで強制ダンプする。
void dumpSolutionH5_force(const solverConfig& , const mesh& , variables& , const int& , const std::string& );


// 境界の h5 (outputHDFflg: 1 の境界)。間隔は outStepInterval、最後の step は必ず書く。force は上と同じ。
void outputBconds_H5_XDMF(const solverConfig& cfg , mesh& msh , variables& var , const int& iStep , bool force = false);
