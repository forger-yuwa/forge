// Optional KD-Tree based wall distance calculation.
// If k-d tree library is not available (HAVE_KDTREE undefined), a brute force
// fallback will be used. This removes the previous hard-coded absolute path.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#ifdef HAVE_KDTREE
#include <kdtree.h>
#include <memory>
#endif

#include "flowFormat.hpp"
#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

struct Point {
  geom_float x, y, z;
  Point(geom_float x_, geom_float y_, geom_float z_) : x(x_), y(y_), z(z_) {}
};

void calcWallDistance_kdtree(solverConfig &cfg, mesh &msh, variables &var);