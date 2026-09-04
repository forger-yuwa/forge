#include <iostream>
#include <vector>

#include "input/solverConfig.hpp"
#include "input/setInitial.hpp"

#include "mesh/mesh.hpp"
#include "mesh/gmshReader.hpp"
#include "boundaryCond.hpp"

using namespace std;

int main(int argc , char *argv[]) 
{
    if (argc != 3) 
    {
        cerr << "usage: convertGmshToNagare gmshFileName inputMeshName \n"; 
    }

    cout << "-------------------------- \n";
    cout << "*** Read Solver Config *** \n";
    cout << "-------------------------- \n";
    solverConfig cfg = solverConfig();
    string fname = "solverConfig.yaml";
    cfg.read(fname);

    cout << "----------------- \n";
    cout << "*** Read Mesh *** \n";
    cout << "----------------- \n";
    gmshReader gmsh = gmshReader(argv[1]);

    cout << "-------------------------------- \n";
    cout << "*** Read Boundary Conditions *** \n";
    cout << "-------------------------------- \n";
    readBcondConfig(cfg , gmsh.bconds);

    // node-centered (median-dual) モードでは、双対メッシュを構築して primal を置き換え、
    // 「双対メッシュを primary mesh」として書き出す (solver 側は無変更で node-centered を扱える)。
    // 初期値は置き換え後の CV (=ノード) 上で設定する必要があるため、setInitial より前に置換する。
    if (cfg.discretization == "node") {
        cout << "-------------------------------- \n";
        cout << "*** Build Median-Dual Mesh   *** \n";
        cout << "-------------------------------- \n";
        gmsh.axisCentroidShift = (cfg.axisCentroidShift != 0);
        gmsh.inletCornerWall   = (cfg.nodeInletCornerWall != 0);
        gmsh.buildMedianDual();
        gmsh.replacePrimalWithDual();
    }

    cout << "-------------------------- \n";
    cout << "*** Set Initial Values *** \n";
    cout << "-------------------------- \n";
    variables var = variables();
    var.allocVariables(cfg.gpu , gmsh);
    setInitial(cfg , gmsh , var);

    cout << "------------------------ \n";
    cout << "*** Write Input HDF5 *** \n";
    cout << "------------------------ \n";
    gmsh.writeInputH5(argv[2] , var);

    return 0;
}