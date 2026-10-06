#include <iostream>
#include <vector>

#include "input/solverConfig.hpp"
#include "input/setInitial.hpp"
#include "input/speciesDB.hpp"

#include "mesh/mesh.hpp"
#include "mesh/gmshReader.hpp"
#include "mesh/memlog.hpp"
#include "boundaryCond.hpp"

using namespace std;

int main(int argc , char *argv[]) 
{
    if (argc != 3) 
    {
        cerr << "usage: convertGmshToNagare gmshFileName inputMeshName \n"; 
    }

    // FORGE_MEMLOG=1 で工程別の VmRSS/VmHWM を出す (mesh/memlog.hpp; plan tooling-sern-mesh-blocking B4-5 (1))
    MEMLOG("開始", std::string());

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
    MEMLOG("gmshReader 構築後 (読込み+makeMesh)", gmsh.memSummary());

    cout << "-------------------------------- \n";
    cout << "*** Read Boundary Conditions *** \n";
    cout << "-------------------------------- \n";
    readBcondConfig(cfg , gmsh.bconds);
    MEMLOG("readBcondConfig 後", std::string());

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
        MEMLOG("buildMedianDual 後", gmsh.memSummary());
        gmsh.replacePrimalWithDual();
        MEMLOG("replacePrimalWithDual 後", gmsh.memSummary());
    }

    cout << "-------------------------- \n";
    cout << "*** Set Initial Values *** \n";
    cout << "-------------------------- \n";
    variables var = variables();
    var.allocVariables(cfg.gpu , gmsh);
    MEMLOG("allocVariables 後", [&]{ size_t b = 0; for (const auto& kv : var.c) b += memlog::flatBytes(kv.second);
                                     size_t bp = 0; for (const auto& kv : var.p) bp += memlog::flatBytes(kv.second);
                                     return memlog::item("var.c(cell 変数)", var.c.size(), b) + " " + memlog::item("var.p(plane 変数)", var.p.size(), bp); }());
    setInitial(cfg , gmsh , var);
    MEMLOG("setInitial 後 (壁距離 kd-tree を含む)", std::string());

    cout << "------------------------ \n";
    cout << "*** Write Input HDF5 *** \n";
    cout << "------------------------ \n";
    gmsh.writeInputH5(argv[2] , var);
    MEMLOG("writeInputH5 後", std::string());

    return 0;
}