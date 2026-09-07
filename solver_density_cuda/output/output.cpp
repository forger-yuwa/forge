#include "output.hpp"

#include <iostream>
#include <fstream>

#include <vector>
#include <string>
#include <list>
#include <algorithm>

#include "mesh/mesh.hpp"
#include "mesh/elementType.hpp"
#include "common/vectorUtil.hpp"

#include <highfive/H5File.hpp>

using HighFive::File;

namespace {

flow_float outputTimeValue(const solverConfig& cfg, int iStep)
{
    if (cfg.unsteady == 1) {
        return cfg.totalTime;
    }

    return static_cast<flow_float>(iStep);
}

}

// 出力する場の量を config output.level で絞る (procedures/solver-settings.md「output」)。
//   level 2: output_cellValNames 全部 (従来)。level 0/1: 下の基本集合 + extraFields を output_cellValNames の順で。
//   h0 (全エンタルピー) は level>=1 で合成出力 (Ht [+k]) し、属性 h0_includes_k を付ける。
static std::list<std::string> effectiveOutputNames(const solverConfig& cfg, const variables& var)
{
    if (cfg.outputLevel >= 2) return var.output_cellValNames;
    std::vector<std::string> base = {"ro","roUx","roUy","roUz","roe","roK","roOmega"};
    for (const auto& n : var.speciesVarNames) base.push_back(n);            // roY{s}
    for (const auto& n : var.condMomentConsNames) base.push_back(n);        // 凝縮モーメント保存量
    if (cfg.outputLevel >= 1) {
        for (const char* n : {"P","T","Ux","Uy","Uz","k","omega","sonic","vis_lam","vis_turb","wall_dist"}) base.push_back(n);
        for (const auto& n : var.speciesVarNames) base.push_back(n.substr(2));   // Y{s}
        for (const auto& n : var.condMomentConsNames) base.push_back(n.substr(2));
    }
    for (const auto& n : cfg.outputExtraFields) base.push_back(n);
    std::list<std::string> out;
    for (const auto& n : var.output_cellValNames) {
        if (std::find(base.begin(), base.end(), n) != base.end()) out.push_back(n);
    }
    for (const auto& n : cfg.outputExtraFields) {
        if (std::find(var.output_cellValNames.begin(), var.output_cellValNames.end(), n) == var.output_cellValNames.end()) {
            static bool warned = false;
            if (!warned) { std::cerr << "[output] extraFields: '" << n << "' is not an output variable (ignored)\n"; warned = true; }
        }
    }
    return out;
}

static void writeSolutionH5_XDMF(const solverConfig& cfg , const mesh& msh , variables& var , const int& iStep , const std::string& prefix)
{
    const std::list<std::string> outNames = effectiveOutputNames(cfg, var);
    const bool writeH0 = (cfg.outputLevel >= 1) && var.c.count("Ht") && var.c.count("k");
    const bool h0IncludesK = (cfg.sstEnergyIncludesK != 0 && cfg.LESorRANS == 2 && cfg.RANSmodel == 1);
    if (cfg.gpu == 1) {
        std::list<std::string> d2h = outNames;
        if (writeH0) { d2h.push_back("Ht"); d2h.push_back("k"); }
        var.copyVariables_cell_D2H(d2h);
    }

    elementTypeMap eleTypeMap;
    
    ostringstream oss;

    // ------------
    // *** HDF5 *** 
    // ------------
    oss << iStep;
    string fnameH5 = prefix+oss.str()+".h5";
    ofstream ofsH5(fnameH5);

    File file(fnameH5, File::ReadWrite | File::Truncate);

    // write mesh structure
    vector<geom_float> COORD;
    for (auto& nod : msh.nodes)
    {
        COORD.push_back(nod.coords[0]);
        COORD.push_back(nod.coords[1]);
        COORD.push_back(nod.coords[2]);
    }

    file.createDataSet("/MESH/COORD",COORD);

    // node-centered (median-dual) モードでは msh.cells が双対 CV なので primal トポロジを
    // 持たない。msh.vizCONNE (primal セル接続) があればそれを使い、属性を Center='Node' で書く
    // (CV index == primal node index)。無ければ従来の cell-centered (Center='Cell')。
    const bool nodeViz = !msh.vizCONNE.empty();

    vector<geom_int> CONNE;
    geom_int CONNE_dim = 0;
    geom_int nElements = 0;
    if (nodeViz)
    {
        CONNE = msh.vizCONNE;
        CONNE_dim = msh.vizCONNE_dim;
        nElements = msh.nVizCells;
    }
    else
    {
        geom_int CONNE0;
        for (geom_int ic = 0 ; ic<msh.nCells ; ic++)
        {
            auto cell = msh.cells[ic];
            geom_int nn = eleTypeMap.mapElementFromGmshID[cell.ieleType].nNodes;
            string name = eleTypeMap.mapElementFromGmshID[cell.ieleType].name;

            if (name == "hex") CONNE0 = 9;
            if (name == "prism") CONNE0 = 8;
            if (name == "pyramid") CONNE0 = 7;
            if (name == "tetra") CONNE0 = 6;
            if (name == "quad") CONNE0 = 5;
            if (name == "triangle") CONNE0 = 4;

            CONNE.push_back(CONNE0);
            CONNE_dim += nn + 1;

            for (auto& nod : cell.iNodes)
            {
                CONNE.push_back(nod);
            }
        }
        nElements = msh.nCells;
    }
    file.createDataSet("/MESH/CONNE",CONNE);

    // write variables
    //for (string name : var.output_cellValNames)
    //TODO: ややこしくしているのでシンプルにoutput_cellValNamesで回したい。が、なぜかエラーになる
    for (auto& v : var.c)
    {
        string name = v.first;

        auto itr = std::find(outNames.begin(), outNames.end(), name);
        if (itr == outNames.end()) {
            continue; // notfound
        }

        std::vector<flow_float> vtemp;
        vtemp.resize(msh.nCells);
        copy(v.second.begin(), v.second.begin()+msh.nCells, vtemp.begin());

        //file.createDataSet("/VALUE/"+name , v.second);
        file.createDataSet("/VALUE/"+name , vtemp);
    }
    // h0: 全エンタルピー (単位質量) = Ht = e + p/ρ + u²/2 (+ k は sstEnergyIncludesK のときだけ)。
    // 全温・全圧の後処理はこれを逆算して作る (自前で T + u²/2c_p を組まない: plan output-level-and-h0)。
    if (writeH0) {
        std::vector<flow_float> h0(msh.nCells);
        const auto& Ht = var.c.at("Ht"); const auto& kk = var.c.at("k");
        for (geom_int i = 0; i < msh.nCells; ++i) h0[i] = Ht[i] + (h0IncludesK ? std::max(kk[i], (flow_float)0.0) : (flow_float)0.0);
        auto ds = file.createDataSet("/VALUE/h0", h0);
        const int flag = h0IncludesK ? 1 : 0;
        ds.createAttribute<int>("h0_includes_k", HighFive::DataSpace::From(flag)).write(flag);
    }

    // ------------
    // *** XDMF ***
    // ------------
    string fnameXDMF = prefix+oss.str()+".xmf";
    ofstream ofs(fnameXDMF);

    ofs << "<?xml version='1.0' ?>\n";
    ofs << "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>\n";
    ofs << "<Xdmf>\n";
    ofs << "  <Domain>\n";
    ofs << "    <Grid GridType='Collection' CollectionType='Spatial' Name='Mixed'>\n";
    ofs << "    <Time TimeType='Single' Value='" << outputTimeValue(cfg, iStep) <<"'/> \n";
    ofs << "      <Grid Name='gridooo'>\n";
    ofs << "        <Topology Type='Mixed' NumberOfElements='" << nElements << "'>\n";
    ofs << "          <DataItem Format='HDF' DataType='Int' Dimensions='" << CONNE_dim << "'>\n";
    ofs << "            " << prefix << oss.str() <<".h5:MESH/CONNE\n";
    ofs << "          </DataItem>\n";
    ofs << "        </Topology>\n";

    ofs << "        <Geometry Type='XYZ'>\n";
    ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << msh.nNodes*3 << "'>\n";
    ofs << "            " << prefix << oss.str() <<".h5:MESH/COORD\n";
    ofs << "          </DataItem>\n";
    ofs << "        </Geometry>\n";

    const char* centerAttr = nodeViz ? "Node" : "Cell";
    std::list<std::string> xmfNames = outNames;
    if (writeH0) xmfNames.push_back("h0");
    for (string name : xmfNames)
    {
    //for (auto& v : var.c) {
        //string name = v.first;
        ofs << "        <Attribute Name='"  << name << "' Center='" << centerAttr << "' >\n";
        ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << msh.nCells << "'>\n";
        ofs << "            " << prefix << oss.str() << ".h5:VALUE/" << name << "\n";
        ofs << "          </DataItem>\n";
        ofs << "        </Attribute>\n";
    }

    ofs << "      </Grid>\n";
    ofs << "    </Grid>\n";
    ofs << "  </Domain>\n";
    ofs << "</Xdmf>\n";
}

void outputH5_XDMF(const solverConfig& cfg , const mesh& msh , variables& var , const int& iStep)
{
    if (iStep%cfg.outStepInterval != 0 or iStep < cfg.outStepStart) return;

    writeSolutionH5_XDMF(cfg , msh , var , iStep , "res_");
}

// detectNaN 診断モード用: 出力間隔ガードを無視して現在の解を強制ダンプする。
// prefix で通常出力 (res_) と区別する (例: res_nan_)。
void dumpSolutionH5_force(const solverConfig& cfg , const mesh& msh , variables& var , const int& iStep , const std::string& prefix)
{
    writeSolutionH5_XDMF(cfg , msh , var , iStep , prefix);
}

void outputBconds_H5_XDMF(const solverConfig& cfg , mesh& msh , variables& var , const int& iStep)
{
    if (iStep%cfg.outStepInterval != 0 or iStep < cfg.outStepStart) return;


    for (auto& bc : msh.bconds) {
        //if (bc.bcondKind  != "wall") continue;
        if (bc.outputHDFflg != 1) continue;

        bc.copyVariables_bplane_D2H();

        elementTypeMap eleTypeMap;
        ostringstream oss;
        ostringstream oss_id;

        // ------------
        // *** HDF5 *** 
        // ------------
        bc.output_preparation(msh.nodes, msh.planes);

        oss << iStep;
        oss_id << bc.physID;
        string fnameH5 = "res_"+bc.physName+"_"+oss_id.str()+"_"+oss.str()+".h5";
        ofstream ofsH5(fnameH5);

        File file(fnameH5, File::ReadWrite | File::Truncate);

        // write boundary
        vector<geom_float> COORD;
        geom_int inl;

        for (geom_int inl=0; inl<bc.inodes_l2g.size(); inl++) {
            geom_int ing = bc.inodes_l2g[inl];
            COORD.push_back(msh.nodes[ing].coords[0]);
            COORD.push_back(msh.nodes[ing].coords[1]);
            COORD.push_back(msh.nodes[ing].coords[2]);
        }

        file.createDataSet("/MESH/COORD",COORD);

        // --- 壁サーフェスのトポロジ (XDMF Mixed) ---
        // cell モード: planes_local の面 (2D=line, 3D=tri/quad)、値は面単位 (Center='Cell')。
        // node モード: 境界も半割双対面(1ノード)になるため、replacePrimalWithDual で退避した primal 境界面
        //   (bc.vizBfaceNodes) で壁ノードを接続し、値はノード単位に並べ替えて Center='Node' で書く。
        // XDMF Mixed の可変ノード型 (Polyvertex=1/Polyline=2) は [type, nNodes, node...] とノード数が要る。
        //   固定型 (Triangle=4/Quad=5) は [type, node...]。(旧版は nn 分岐漏れで未初期化型を書き ParaView 不可視だった)
        const bool nodeWallViz = !bc.vizBfaceNodes.empty();

        vector<geom_int> CONNE;
        geom_int CONNE_dim = 0;
        auto pushElem = [&](const vector<geom_int>& gnodes) {
            geom_int nn = (geom_int)gnodes.size();
            if      (nn == 1) { CONNE.push_back(1); CONNE.push_back(1);  CONNE_dim += 3; }      // Polyvertex
            else if (nn == 2) { CONNE.push_back(2); CONNE.push_back(2);  CONNE_dim += 4; }      // Polyline
            else if (nn == 3) { CONNE.push_back(4);                      CONNE_dim += nn + 1; } // Triangle
            else if (nn == 4) { CONNE.push_back(5);                      CONNE_dim += nn + 1; } // Quadrilateral
            else              { CONNE.push_back(1); CONNE.push_back(nn); CONNE_dim += nn + 2; } // fallback Polyvertex
            for (const geom_int g : gnodes) CONNE.push_back(bc.inodes_g2l[g]);
        };

        geom_int nElements;
        geom_int attrDim;
        string attrCenter;
        if (nodeWallViz) {
            for (auto& face : bc.vizBfaceNodes) pushElem(face);
            nElements  = (geom_int)bc.vizBfaceNodes.size();
            attrDim    = (geom_int)bc.inodes_l2g.size();
            attrCenter = "Node";
        } else {
            for (auto& pln : bc.planes_local) pushElem(pln.iNodes);
            nElements  = (geom_int)bc.planes_local.size();
            attrDim    = (geom_int)bc.planes_local.size();
            attrCenter = "Cell";
        }
        file.createDataSet("/MESH/CONNE",CONNE);

        // write variables
        for (auto& vt : bc.valueTypes)
        {
            string name = vt.first;
            std::vector<flow_float> vtemp;
            if (nodeWallViz) {
                // bvar は bplane(=半割面=ノード) 順。Center='Node' 用に局所ノード順 (inodes_l2g) へ並べ替える。
                vtemp.assign(bc.inodes_l2g.size(), (flow_float)0.0);
                for (geom_int j = 0; j < (geom_int)bc.planes_local.size(); j++) {
                    if (bc.planes_local[j].iNodes.empty()) continue;
                    geom_int L = bc.inodes_g2l[bc.planes_local[j].iNodes[0]];
                    if (L >= 0) vtemp[L] = bc.bvar[name][j];
                }
            } else {
                vtemp.resize(bc.planes_local.size());
                copy(bc.bvar[name].begin(), bc.bvar[name].begin()+bc.planes_local.size(), vtemp.begin());
            }
            file.createDataSet("/VALUE/"+name , vtemp);
        }

//
//        for (auto& v : bc.bvar)
//        {
//            string name = v.first;
//
//
//            //auto itr = std::find(var.output_cellValNames.begin(), var.output_cellValNames.end(), name);
//            //if (itr == var.output_cellValNames.end()) {
//            //    continue; // notfound
//            //}
//
//            std::vector<flow_float> vtemp;
//            vtemp.resize(bc.planes_local.size());
//            copy(v.second.begin(), v.second.begin()+bc.planes_local.size(), vtemp.begin());
//
//            file.createDataSet("/VALUE/"+name , vtemp);
//        }
//
        //ghost // ------------------
        //ghost // *** Ghost cell ***
        //ghost // ------------------

        //ghost vector<geom_float> COORD_ghst;

        //ghost // out ghost cell 
        //ghost for (auto& ighst : bc.iCells_ghst) {
        //ghost     // write boundary
        //ghost     COORD_ghst.push_back(msh.cells[ighst].centCoords[0]);
        //ghost     COORD_ghst.push_back(msh.cells[ighst].centCoords[1]);
        //ghost     COORD_ghst.push_back(msh.cells[ighst].centCoords[2]);
        //ghost }

        //ghost file.createDataSet("/MESH_ghst/COORD",COORD_ghst);

        //ghost vector<geom_int> CONNE_ghst;
        //ghost vector<geom_int> Indexes_ghst;
        //ghost geom_int CONNE_dim_ghst = 0;
        //ghost geom_int CONNE0_ghst;    

        //ghost geom_int ighst_l = 0;
        //ghost for (auto& ighst : bc.iCells_ghst) {
        //ghost     geom_int nn = 1;
        //ghost     string name = "point";
        //ghost     CONNE0_ghst = 1;

        //ghost     CONNE_ghst.push_back(CONNE0_ghst);
        //ghost     CONNE_dim_ghst += nn + 1;

        //ghost     CONNE_ghst.push_back(ighst_l);
        //ghost     Indexes_ghst.push_back(ighst_l);
        //ghost     ighst_l++;
        //ghost }
        //ghost file.createDataSet("/MESH_ghst/CONNE",CONNE_ghst);
        //ghost file.createDataSet("/MESH_ghst/Indexes",Indexes_ghst);

        //ghost // write variables
        //ghost //for (string name : var.output_cellValNames)
        //ghost //TODO: ややこしくしているのでシンプルにoutput_cellValNamesで回したい。が、なぜかエラーになる
        //ghost for (auto& v : var.c)
        //ghost {
        //ghost     string name = v.first;
    
        //ghost     auto itr = std::find(var.output_cellValNames.begin(), var.output_cellValNames.end(), name);
        //ghost     if (itr == var.output_cellValNames.end()) {
        //ghost         continue; // notfound
        //ghost     }
    
        //ghost     std::vector<flow_float> vtemp;
        //ghost     vtemp.resize(bc.iCells_ghst.size());
        //ghost     for (geom_int igl = 0 ; igl<bc.iCells_ghst.size() ; igl++) {
        //ghost         geom_int ig = bc.iCells_ghst[igl];
        //ghost         vtemp[igl] = v.second[ig];
        //ghost     }
    
        //ghost     //file.createDataSet("/VALUE/"+name , v.second);
        //ghost     file.createDataSet("/VALUE_ghst/"+name , vtemp);
        //ghost }

        // ------------
        // *** XDMF ***
        // ------------
        //string fnameXDMF = "res_wall_"+oss.str()+".xmf";

        string fnameXDMF = "res_"+bc.physName+"_"+oss_id.str()+"_"+oss.str()+".xmf";
        ofstream ofs(fnameXDMF);

        ofs << "<?xml version='1.0' ?>\n";
        ofs << "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>\n";
        ofs << "<Xdmf>\n";
        ofs << "  <Domain>\n";
        ofs << "    <Grid  GridType='Collection' CollectionType='Spatial' Name='Mixed'>\n";
        ofs << "    <Time TimeType='Single' Value='" << outputTimeValue(cfg, iStep) << "'/> \n";
        ofs << "      <Grid Name='wall'>\n";
        ofs << "        <Topology Type='Mixed' NumberOfElements='" << nElements << "'>\n";
        ofs << "          <DataItem Format='HDF' DataType='Int' Dimensions='" << CONNE_dim << "'>\n";
        ofs << "            res_" << bc.physName << "_" << oss_id.str() <<"_" << oss.str() <<".h5:MESH/CONNE\n";
        ofs << "          </DataItem>\n";
        ofs << "        </Topology>\n";

        ofs << "        <Geometry Type='XYZ'>\n";
        ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << bc.inodes_l2g.size()*3 << "'>\n";
        //ofs << "            res_wall_"<< oss.str() <<".h5:MESH/COORD\n";
        ofs << "            res_" << bc.physName << "_" << oss_id.str() <<"_" << oss.str() <<".h5:MESH/COORD\n";
        ofs << "          </DataItem>\n";
        ofs << "        </Geometry>\n";


        for (auto& vt : bc.valueTypes)
        {
            string name = vt.first;
            ofs << "        <Attribute Name='"  << name << "' Center='" << attrCenter << "' >\n";
            ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << attrDim << "'>\n";
            ofs << "            res_" << bc.physName << "_" << oss_id.str() <<"_" << oss.str() <<".h5:VALUE/" << name << "\n";
            ofs << "          </DataItem>\n";
            ofs << "        </Attribute>\n";
        }
        ofs << "      </Grid>\n";
        //ghost // ghost cell
        //ghost ofs << "      <Grid Name='wall_ghst'>\n";
        //ghost ofs << "        <Topology Type='Polyvertex' Dimensions='" << bc.iCells_ghst.size() <<"' NodesPerElement='1'>\n";
        //ghost //ofs << "          <DataItem Format='HDF' DataType='Int' Dimensions='" << CONNE_dim_ghst << "'>\n";
        //ghost ofs << "          <DataItem Format='HDF' DataType='Int' Dimensions='" << bc.iCells_ghst.size() << "'>\n";
        //ghost //ofs << "            res_wall_" << oss.str() <<".h5:MESH_ghst/CONNE\n";
        //ghost ofs << "            res_wall_" << oss.str() <<".h5:MESH_ghst/Indexes\n";
        //ghost ofs << "          </DataItem>\n";
        //ghost ofs << "        </Topology>\n";

        //ghost ofs << "        <Geometry Type='XYZ'>\n";
        //ghost ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << bc.iCells_ghst.size()*3 << "'>\n";
        //ghost ofs << "            res_wall_"<< oss.str() <<".h5:MESH_ghst/COORD\n";
        //ghost ofs << "          </DataItem>\n";
        //ghost ofs << "        </Geometry>\n";

        //ghost for (string name : var.output_cellValNames)
        //ghost {
        //ghost     ofs << "        <Attribute Name='"  << name << "' Center='Node' >\n";
        //ghost     ofs << "          <DataItem Format='HDF' DataType='Float' Dimensions='" << bc.iCells_ghst.size() << "'>\n";
        //ghost     ofs << "            res_wall_"<< oss.str() << ".h5:VALUE_ghst/" << name << "\n";
        //ghost     ofs << "          </DataItem>\n";
        //ghost     ofs << "        </Attribute>\n";
        //ghost }
        //ghost ofs << "      </Grid>\n";

        ofs << "    </Grid>\n";
        ofs << "  </Domain>\n";
        ofs << "</Xdmf>\n";
    }
}
