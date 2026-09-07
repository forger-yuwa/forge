#include "input/calcWallDistance_kdtree.hpp"

namespace {
inline geom_float dist(const Point &a, const Point &b) {
    const geom_float dx = a.x - b.x;
    const geom_float dy = a.y - b.y;
    const geom_float dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

#ifdef HAVE_KDTREE
void compute_with_kdtree(const std::vector<Point> &walls,
                                                 const std::vector<Point> &cells,
                                                 std::vector<geom_float> &distance) {
    kdtree *tree = kd_create(3);
    std::unique_ptr<int[]> indices(new int[walls.size()]);
    for (int i = 0; i < static_cast<int>(walls.size()); ++i) {
        indices[i] = i;
        kd_insert3(tree, walls[i].x, walls[i].y, walls[i].z, &indices[i]);
    }
    distance.reserve(cells.size());
    for (int j = 0; j < static_cast<int>(cells.size()); ++j) {
        kdres *set = kd_nearest3(tree, cells[j].x, cells[j].y, cells[j].z);
        int i = *(int *)kd_res_item_data(set);
        kd_res_free(set);
        distance.push_back(dist(walls[i], cells[j]));
    }
    kd_free(tree);
}
#endif

void compute_bruteforce(const std::vector<Point> &walls,
                                                const std::vector<Point> &cells,
                                                std::vector<geom_float> &distance) {
    distance.reserve(cells.size());
    for (const auto &c : cells) {
        geom_float dmin = std::numeric_limits<geom_float>::max();
        for (const auto &w : walls) {
            dmin = std::min(dmin, dist(w, c));
        }
        distance.push_back(dmin == std::numeric_limits<geom_float>::max() ? 0.0 : dmin);
    }
}
} // namespace

void calcWallDistance_kdtree(solverConfig &cfg, mesh &msh, variables &var) {
    // node-centered (median-dual): 壁ノードが CV 中心 (=壁面上) で、off-wall ノードはその直上に並ぶ。
    // 壁点集合に「壁半割面の重心」を使うと重心が壁ノードから x 方向に dx/8 ずれており、近壁ノードの
    // 最近接壁点距離が法線距離 y でなく x ずれ (≈dx/8, 下流で増大) に支配され wall_dist が大きく誤る
    // (壁ノード自身も 0 にならない)。→ ω_w=60ν/(β·wall_dist²) 過小→過剰乱流。node では壁点に「壁ノード
    // 座標 (bc.iCells の CV 中心)」を使う。直下の壁ノードが最近接になり wall_dist=y を正しく返す。
    // cell モードは従来どおり壁面の plane 重心 (セル中心の直下に整列) を使う (挙動不変)。
    const bool nodeMode = (cfg.discretization == "node");

    std::vector<Point> wall_points;
    std::vector<Point> cell_points;
    std::vector<geom_float> distance;

    int wall_count = 0;
    for (auto &bc : msh.bconds) {
        if (bc.bcondKind == "wall_isothermal" || bc.bcondKind == "wall") {
            if (nodeMode) {
                for (auto &ic : bc.iCells) {
                    if (ic < 0 || ic >= static_cast<geom_int>(msh.cells.size())) continue;
                    geom_float x = msh.cells[ic].centCoords[0];
                    geom_float y = msh.cells[ic].centCoords[1];
                    geom_float z = msh.cells[ic].centCoords[2];
                    wall_points.emplace_back(x, y, z);
                }
            } else {
                for (auto &ip : bc.iPlanes) {
                    geom_float x = msh.planes[ip].centCoords[0];
                    geom_float y = msh.planes[ip].centCoords[1];
                    geom_float z = msh.planes[ip].centCoords[2];
                    wall_points.emplace_back(x, y, z);
                }
            }
            ++wall_count;
        }
    }

    if (wall_count > 0) {
        cell_points.reserve(msh.cells.size());
        for (auto &icell : msh.cells) {
            geom_float x = icell.centCoords[0];
            geom_float y = icell.centCoords[1];
            geom_float z = icell.centCoords[2];
            cell_points.emplace_back(x, y, z);
        }

#ifdef HAVE_KDTREE
        compute_with_kdtree(wall_points, cell_points, distance);
#else
        compute_bruteforce(wall_points, cell_points, distance);
#endif

        // Copy into variable array (assumes var.c["wall_dist"] sized nCells_all)
        if (var.c.count("wall_dist") && var.c["wall_dist"].size() >= distance.size()) {
            std::copy(distance.begin(), distance.end(), var.c["wall_dist"].begin());
        }
    } else {
        distance.assign(msh.nCells_all, geom_float(0));
        if (var.c.count("wall_dist") && var.c["wall_dist"].size() >= distance.size()) {
            std::copy(distance.begin(), distance.end(), var.c["wall_dist"].begin());
        }
    }

    std::list<std::string> names = {"wall_dist"};
    var.copyVariables_cell_H2D(names);
}