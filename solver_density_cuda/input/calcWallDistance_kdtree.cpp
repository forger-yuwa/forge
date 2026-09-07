#include "input/calcWallDistance_kdtree.hpp"
#include <cstdlib>

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

// 内蔵 k-d tree (外部ライブラリ不要)。壁点集合を再帰的に中央値分割し、最近接探索は枝刈り付き。
// 2026-09-08: 外部 kdtree が見つからないビルド (WSL native / AWS) はブルートフォース O(N·M) に落ちていて、
// 2M 節点 × 13 万壁点で数分〜十数分かかっていた (変換時間の主因)。本実装で数秒。
// 最近接距離は一意 (同距離のタイは距離値が同じ) なので結果はブルートフォースと一致する (float 丸めの範囲)。
struct BuiltinKdTree {
    struct Node { int lo, hi; int axis; geom_float split; int left, right; };
    const std::vector<Point> &pts;
    std::vector<int> idx;
    std::vector<Node> nodes;
    static constexpr int LEAF = 8;
    explicit BuiltinKdTree(const std::vector<Point> &p) : pts(p), idx(p.size()) {
        for (size_t i = 0; i < p.size(); ++i) idx[i] = (int)i;
        if (!p.empty()) build(0, (int)p.size(), 0);
    }
    static geom_float coord(const Point &q, int ax) { return ax == 0 ? q.x : (ax == 1 ? q.y : q.z); }
    int build(int lo, int hi, int depth) {
        Node nd; nd.lo = lo; nd.hi = hi; nd.left = nd.right = -1; nd.axis = -1; nd.split = 0;
        const int me = (int)nodes.size(); nodes.push_back(nd);
        if (hi - lo > LEAF) {
            // 分割軸: 範囲が最大の軸
            geom_float mn[3] = { std::numeric_limits<geom_float>::max(), std::numeric_limits<geom_float>::max(), std::numeric_limits<geom_float>::max() };
            geom_float mx[3] = { -mn[0], -mn[1], -mn[2] };
            for (int i = lo; i < hi; ++i) for (int a = 0; a < 3; ++a) { const geom_float v = coord(pts[idx[i]], a); mn[a] = std::min(mn[a], v); mx[a] = std::max(mx[a], v); }
            int ax = 0; for (int a = 1; a < 3; ++a) if (mx[a] - mn[a] > mx[ax] - mn[ax]) ax = a;
            const int mid = (lo + hi) / 2;
            std::nth_element(idx.begin() + lo, idx.begin() + mid, idx.begin() + hi,
                             [&](int a, int b) { return coord(pts[a], ax) < coord(pts[b], ax); });
            nodes[me].axis = ax; nodes[me].split = coord(pts[idx[mid]], ax);
            const int l = build(lo, mid, depth + 1);
            const int r = build(mid, hi, depth + 1);
            nodes[me].left = l; nodes[me].right = r;
        }
        return me;
    }
    void nearest(int n, const Point &q, geom_float &best2) const {
        const Node &nd = nodes[n];
        if (nd.axis < 0) {
            for (int i = nd.lo; i < nd.hi; ++i) {
                const Point &w = pts[idx[i]];
                const geom_float dx = w.x - q.x, dy = w.y - q.y, dz = w.z - q.z;
                best2 = std::min(best2, dx*dx + dy*dy + dz*dz);
            }
            return;
        }
        const geom_float d = coord(q, nd.axis) - nd.split;
        const int first = (d < 0) ? nd.left : nd.right, second = (d < 0) ? nd.right : nd.left;
        nearest(first, q, best2);
        if (d * d < best2) nearest(second, q, best2);
    }
};

void compute_with_builtin_kdtree(const std::vector<Point> &walls,
                                 const std::vector<Point> &cells,
                                 std::vector<geom_float> &distance) {
    if (walls.empty()) { distance.assign(cells.size(), geom_float(0)); return; }
    BuiltinKdTree tree(walls);
    distance.resize(cells.size());
    #pragma omp parallel for schedule(static)
    for (long j = 0; j < (long)cells.size(); ++j) {
        geom_float best2 = std::numeric_limits<geom_float>::max();
        tree.nearest(0, cells[j], best2);
        distance[j] = std::sqrt(best2);
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
        (void)compute_bruteforce;   // 参照用に残す (検証時は FORGE_WALLDIST_BRUTE=1 で切替)
        if (std::getenv("FORGE_WALLDIST_BRUTE") != nullptr)
            compute_bruteforce(wall_points, cell_points, distance);
        else
            compute_with_builtin_kdtree(wall_points, cell_points, distance);
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