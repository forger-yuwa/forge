#pragma once

#include <iostream>
#include <fstream>

#include <vector>
#include <list>

#include "flowFormat.hpp"
#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"
#include <memory>
#include </home/kumpei/app/kdtree/install/include/kdtree.h>


struct Point {
  double x, y, z;
  Point(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

class point_probes {
public:
    int nProbes;
    std::vector<geom_int > cell_id;
    geom_int* cell_id_d;
    std::map<std::string, std::vector<flow_float>> var; // host probe variables
    std::map<std::string, flow_float*> var_d; // device probe variables

    std::vector<Point> inputProbeXYZ;
    std::vector<Point> nearestCellXYZ;

    int outStepStart; 
    int outStepInterval; 
    std::string outFileName;

    const std::list<std::string> valNames = 
    {
        "T" , "P" , "Ux" , "Uy" , "Uz" 
    };

    point_probes();

    //variables(const int& ,  mesh&);
    ~point_probes();

    void allocVariables(const int &useGPU , mesh& msh);


    void copyVariables_probes_H2D_all();

    void copyVariables_probes_D2H_all();

    void setNearestCell(solverConfig& cfg, cudaConfig& cuda_cfg , mesh& msh);

    void outputProbes(solverConfig& cfg,mesh& msh , variables& var);
};
