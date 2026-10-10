#pragma once

#include <cstddef>
#include <vector>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"

// 変換器 (convertGmshToForge / gmshReader) 専用の幾何の正本 (double)。
// plans/active/architecture-float-state-double-geometry.md §4.6 (案 (b)、codex diagnose 2026-10-10)。
//
// ビルドの型 (geom_float) に依らず、座標の読み取り (stod) から HDF5 の出力までを、この double の配列で持つ。
// 共用の mesh (node/plane/cell の geom_float) へは roundTo* で一方向に丸めて渡すだけで、
// 丸めた値から幾何を計算し直したり、丸めた値を HDF5 に書いたりはしない。
// FP64 のビルド (geom_float = double) では丸めは恒等なので、共用の mesh の値も正本と同じになる。
//
// 並びは HDF5 と同じ [3*i+k]。node モードでは replacePrimalWithDual が双対の量で置き換える
// (CV = 節点、面 = 内部双対面 + 境界半割面)。
struct convGeom64
{
    std::vector<double> nodeCoord;     // [3*nNodes]  節点座標 (/MESH/COORD)
    std::vector<double> planeSurfVect; // [3*nPlanes] 面ベクトル (整向後, /PLANES/surfVect)
    std::vector<double> planeSurfArea; // [nPlanes]   面積 (/PLANES/surfArea)
    std::vector<double> planeCent;     // [3*nPlanes] 面重心 (/PLANES/centCoords)
    std::vector<double> cellVolume;    // [nCells]    体積 (/CELLS/volume)
    std::vector<double> cellCent;      // [3*nCells]  セル (CV) 重心 (/CELLS/centCoords)
    // 区間ごとの r 重みの面ベクトル W_f = Σ_k r_k S_k (/PLANES/rSurfVect、向きは planeSurfVect と同じ)。
    // 2D の node の格子で、変換の solverConfig が isAxisymmetric 1 のときだけ持つ (それ以外は空で、書かない)。
    // plans/active/axisymmetric-freestream-hoop-gauge.md §4.5。共用の plane へは渡さない (ソルバが HDF5 から読む)。
    std::vector<double> planeRSurfVect; // [3*nPlanes] または空

    // 共用の plane へ一方向に丸めて渡す (面ベクトル・面積・面重心)。
    void roundToShared(std::vector<plane>& planes) const
    {
        for (size_t ip = 0; ip < planes.size(); ++ip) {
            plane& p = planes[ip];
            p.surfVect.resize(3);
            p.centCoords.resize(3);
            for (int k = 0; k < 3; ++k) {
                p.surfVect[k]   = static_cast<geom_float>(planeSurfVect[3*ip + k]);
                p.centCoords[k] = static_cast<geom_float>(planeCent[3*ip + k]);
            }
            p.surfArea = static_cast<geom_float>(planeSurfArea[ip]);
        }
    }

    // 共用の cell へ一方向に丸めて渡す (体積・重心)。
    void roundToShared(std::vector<cell>& cells) const
    {
        for (size_t ic = 0; ic < cells.size(); ++ic) {
            cell& c = cells[ic];
            c.centCoords.resize(3);
            for (int k = 0; k < 3; ++k) c.centCoords[k] = static_cast<geom_float>(cellCent[3*ic + k]);
            c.volume = static_cast<geom_float>(cellVolume[ic]);
        }
    }
};
