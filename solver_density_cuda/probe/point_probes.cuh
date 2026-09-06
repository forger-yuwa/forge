#pragma once

#include <iostream>
#include <fstream>

#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

#include <vector>
#include <list>

#include "flowFormat.hpp"
#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include <memory>
#ifdef HAVE_KDTREE
#include <kdtree.h>
#endif


struct Point_r1 {
  double x, y, z;
  Point_r1(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

class point_probes {
public:
    int nProbes;
    std::vector<geom_int > cell_id;
    geom_int* cell_id_d;
    std::map<std::string, std::vector<flow_float>> var; // host probe variables
    std::map<std::string, flow_float*> var_d; // device probe variables

    std::vector<Point_r1> inputProbeXYZ;
    std::vector<Point_r1> nearestCellXYZ;

    int outStepStart; 
    int outStepInterval; 
    std::string outFileName;

    dim3 dimGrid_pprobe;

    const std::list<std::string> valNames = 
    {
        "T" , "P" , "Ux" , "Uy" , "Uz" 
    };

    point_probes();

    //variables(const int& ,  mesh&);
    ~point_probes();

    void init(solverConfig& cfg, cudaConfig& cuda_cfg , mesh& msh);

    void allocVariables();

    void copyVariables_probes_H2D_all();

    void copyVariables_probes_D2H_all();

    void readYAML();

    void setNearestCell(solverConfig& cfg, cudaConfig& cuda_cfg , mesh& msh);

    void setProbeVariables_d_wrapper(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var);

    void outputProbes(solverConfig& cfg , cudaConfig& cudaCfg , mesh& msh , variables& var , int iStep);
};

