#pragma once
// 遠方境界 farfield の面の値配列 (plan boundary-node-farfield-characteristic §4.3)。実体は convectiveFlux_d.cu。
// 化学種・k・ω の輸送カーネルが、farfield の境界半割面で流入 (質量流束 < 0) のときに運ぶ外側状態の値。
// farfield が無ければ nullptr。name = "Y{s}" / "k" / "omega"。farfield 以外の面は NaN (使わない)。
#include <string>
#include "flowFormat.hpp"
flow_float* farfieldFaceScalar(const std::string& name);
flow_float** farfieldFaceYDevice();
