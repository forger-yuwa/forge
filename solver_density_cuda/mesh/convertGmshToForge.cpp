#include <iostream>
#include <vector>

#include "input/solverConfig.hpp"
#include "input/setInitial.hpp"
#include "input/speciesDB.hpp"

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
    // 化学種 DB の host 側解決 (GPU 無し)。bcond の X{s}→Y{s} 換算が MW を使う。未知種名はここで exit。
    speciesDB_printTable(cfg, speciesDB_init(cfg));

    cout << "----------------- \n";
    cout << "*** Read Mesh *** \n";
    cout << "----------------- \n";
    gmshReader::renumberRCM = (cfg.meshRenumber == "rcm");
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
        gmsh.axisymRWeight     = (cfg.isAxisymmetric == 1);   // 2D は /PLANES/rSurfVect を書く (plan axisymmetric-freestream-hoop-gauge §4.5)
        gmsh.buildMedianDual();
        gmsh.replacePrimalWithDual();
    }

    cout << "-------------------------- \n";
    cout << "*** Set Initial Values *** \n";
    cout << "-------------------------- \n";
    variables var = variables();
    var.allocVariables(cfg.gpu , gmsh);
    setInitial(cfg , gmsh , var);

    // 壁距離の正本 (double)。幾何の正本 gmsh.geo64 の位置から double で計算し、/VALUE/wall_dist はこれから書く
    // (plan architecture-float-state-double-geometry §4.6)。setInitial が var.c["wall_dist"] に入れた値は
    // 共用の mesh (geom_float) の位置から作ったソルバ用の写しで、HDF5 には書かない。
    WallDistPositions64 wallPos;
    wallPos.nodeCoord = gmsh.geo64.nodeCoord.data();
    wallPos.cellCent  = gmsh.geo64.cellCent.data();
    wallPos.planeCent = gmsh.geo64.planeCent.data();
    const std::vector<double> wallDist64 = calcWallDistance64(cfg , gmsh , wallPos);

    cout << "------------------------ \n";
    cout << "*** Write Input HDF5 *** \n";
    cout << "------------------------ \n";
    gmsh.writeInputH5(argv[2] , var , &wallDist64);

    return 0;
}