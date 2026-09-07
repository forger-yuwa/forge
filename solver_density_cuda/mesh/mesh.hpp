#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <string>
//#include <Eigen/Dense>
#include "flowFormat.hpp"
#include "elementType.hpp"

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>

struct node 
{
    std::vector<geom_int> iCells;
    std::vector<geom_int> iPlanes;
    std::vector<geom_float> coords; 

    node();
    node(geom_float& , geom_float& , geom_float& );
};

struct plane 
{
    std::vector<geom_int> iNodes;
    std::vector<geom_int> iCells; 

    std::vector<geom_float> surfVect;
    geom_float surfArea;
    std::vector<geom_float> centCoords;
};

struct cell 
{
    std::vector<geom_int> iNodes;

    std::vector<geom_int> iPlanes;
    std::vector<geom_int> iPlanesDir;

    geom_float volume; 

    geom_int ieleType;

    std::vector<geom_float> centCoords;

    geom_int regionId;

    cell();
    cell(std::vector<geom_int>&) ;
};

struct bcond
{
public:
    std::list<std::string> bplaneValNames =  
    {
        "ro"   , "roUx" , "roUy" , "roUz" , "roUz" , "roe",
        "Ux"   , "Uy"   , "Uz"   , "Tt"   , "Pt"   , "Ts" , "Ps",
        "k"    , "omega", "kb", "omegab",
        "Ux0"  , "Uy0"  , "Uz0"  , 
        "ypls" , "twall_x" , "twall_y" , "twall_z" , "utau" , "qwall"
    };

    std::list<std::string> bplaneIntNames =  
    {
        "partnerPlnID",  // for periodic 
        "partnerCellID" // for periodic 
    };

    geom_int physID; 
    std::string physName; 

    std::vector<geom_int> iPlanes;
    std::vector<geom_int> iBPlanes;
    std::vector<geom_int> iCells; //chk
    std::vector<geom_int> iCells_ghst; //ghost cell

    int outputHDFflg; 
    int output_preparation_flg = 0;
    std::vector<node> nodes_local;
    std::vector<plane> planes_local;
    std::vector<geom_int> inodes_l2g;
    std::vector<geom_int> inodes_g2l;

    // node モード壁可視化用: primal 境界面の global ノードid列を一度だけ退避 (replacePrimalWithDual)。
    // node モードでは境界も半割双対面(1ノード)に置換され面を作れないため、msh.vizCONNE の境界版として
    // primal 面接続を保存し h5 (/BCONDS/<id>/vizBface*) 経由で solver の壁出力へ渡す (Center='Node')。
    // 空なら cell モード (planes_local の面を使い Center='Cell')。
    std::vector<std::vector<geom_int>> vizBfaceNodes;

    std::map<std::string, int> inputInts;
    std::map<std::string, flow_float> inputFloats;

    std::string bcondKind; 

    std::map<std::string,int> valueTypes; // {P, 0}, {T, 1}, etc.

    std::map<std::string, std::vector<flow_float>> bvar;
    std::map<std::string, flow_float* > bvar_d; // cuda

    // 多成分 TP (M5): 入口組成 Y_s^in の device ポインタ配列 (flow_float*[nSpecies])。
    // inlet カーネルが混合則 thermo を計算するため。inlet_uniformVelocity_d_wrapper で
    // bvar_d["Y{s}"] から遅延構築する。単成分では nullptr のまま。
    flow_float** Yb_d = nullptr; // cuda (device array of nSpecies bvar_d["Y{s}"] pointers)

    // interger values (only for periodic now)
    std::map<std::string, std::vector<geom_int>> bint;
    std::map<std::string, geom_int* > bint_d; // cuda

    //cuda
    geom_int* map_bplane_plane_d;
    geom_int* map_bplane_cell_d;
    geom_int* map_bplane_cell_ghst_d;

    bcond();
    bcond(const geom_int& , const std::vector<geom_int>& , 
          const std::vector<geom_int>& , const std::vector<geom_int>& );
    ~bcond();

    void bcondInitVariables(const int &useGPU);
    void copyVariables_bplane_D2H();
    //void set_nodes_local(mesh& msh);
    void output_preparation(std::vector<node>& nodes, std::vector<plane>& planes);
};

struct mesh 
{
public:
    geom_int nNodes , nNodes_all, nNodes_halo;
    geom_int nPlanes, nNormalPlanes, nNormal_halo_Planes, nPlanes_all; 
    geom_int nCells , nCells_halo, nCells_ghst , nCells_all;
    geom_int nBPlanes;
    geom_int nWallHaloPlanes = 0; // normal_halo_planes_d 末尾に並ぶ壁境界 plane 数 (node 弱形式で主ループから除外)
    geom_int nBoundaryHaloPlanes = 0; // 非 periodic 境界 plane 総数 (wall+inlet+outlet+slip...)。node 弱形式 (Phase 2) で
                                      // 主対流ループから全境界を外し convPlaneBound=内部+periodic にするのに使う
    geom_int nBconds;

    std::vector<node> nodes;
    std::vector<plane> planes;
    std::vector<cell> cells;
    std::vector<bcond> bconds;

    // node-centered (median-dual) 用の可視化トポロジ。空でなければ「CV=ノード」モードで、
    // 出力は primal セルのトポロジ + Center='Node' で書く (CV index == primal node index)。
    // cell-centered モードでは空のまま (出力は従来の Center='Cell')。
    // vizCONNE は XDMF mixed-topology 形式 (各セル: [typeCode, node0, node1, ...] の連結)。
    std::vector<geom_int> vizCONNE;
    geom_int vizCONNE_dim = 0; // vizCONNE の総長 (XDMF DataItem Dimensions)
    geom_int nVizCells    = 0; // primal セル数 (XDMF NumberOfElements)

    // cuda
    //geom_int* map_nplane_cells_d; // normal plane
    geom_int* map_plane_cells_d; // 

    geom_int* map_cell_planes_index_d; 
    geom_int* map_cell_planes_d; 

    geom_int* normal_halo_planes_d; //

    // 境界隣接 CV フラグ [nCells] (境界面を持つ CV=1)。node-centered の壁近傍ロバスト化
    // (bndFirstOrder) の診断・対策で、境界ノードの 2 次再構成を 1 次に落とすのに使う。
    geom_int* bnode_flag_d = nullptr;

    // 軸上 CV フラグ [nCells] (ノード R=0 の CV=1)。node-centered 軸対称で、軸上ノードが CV 中心=
    // 特異点になり半径方向運動量が偽駆動されるのを防ぐため、軸上で roUy=0 (対称条件) を課すのに使う。
    geom_int* axis_flag_d = nullptr;


    // nodeValueAtNode (main.cpp が readMesh 前に node なら 1 をセット): readMesh で cells[ic].centCoords ← nodes[ic].coords
    // (ic<nCells) に置換し、置換前の双対 CV 重心 y (=面積加重半径 r̄) を rEff[ic] に退避する。
    // 軸対称 r 重み (variables.cpp) は rEff を使う (centCoords.y=0 の軸ノードで回転体積が消えないように)。
    int nodeValueAtNode = 0;
    std::vector<geom_float> rEff;

    // 壁 CV フラグ [nCells] (wall 種別 bcond の CV=1)。node-centered 壁 Dirichlet で、壁ノード速度を
    // 厳密に 0 に固定する (state 初期化 + 運動量残差射影) のに使う。壁ゴーストを撤廃する代替。
    geom_int* wall_flag_d = nullptr;

    // --- line-implicit (plans/active/time_integration-line-implicit.md) ---
    // 壁法線ライン (積層方向の CV 鎖)。buildImplicitLines() が greedy に構築して上げる。
    // 載らない CV は point-DPLUR fallback (line_prev/next = -1)。
    geom_int nImplicitLines = 0;
    geom_int* line_prev_d = nullptr;      // [nCells] ライン上の前隣 (壁側)。-1 = なし
    geom_int* line_next_d = nullptr;      // [nCells] ライン上の次隣 (外側)。-1 = なし
    geom_int* line_offsets_d = nullptr;   // [nLines+1] CSR
    geom_int* line_cells_d = nullptr;     // [Σ長] ライン順 (壁→終端) の CV id
    flow_float* line_Kprev_d = nullptr;   // [25*nCells] 前隣への近傍行列 K (rhs += K·dq_prev)
    flow_float* line_Knext_d = nullptr;   // [25*nCells] 次隣への近傍行列
    double* line_W_d = nullptr;           // [25*nCells] Thomas scratch (D̃⁻¹·Knext)
    double* line_y_d = nullptr;           // [5*nCells]  Thomas scratch (前進 rhs)
    double* line_LU_d = nullptr;          // [25*nCells] v2: D̃ の LU 因子 (factor/solve 分離)
    signed char* line_piv_d = nullptr;    // [5*nCells]  v2: ピボット
    unsigned char* line_fail_d = nullptr; // [nLines]    v2: factor 失敗フラグ (solve は dq 据え置き)
    unsigned char* plane_wall_flag_d = nullptr; // [nPlanes] wall 種 bcond の境界面フラグ (lineDtWallRelief 診断用)
    void buildImplicitLines(const flow_float* ccx, const flow_float* ccy, const flow_float* ccz);

    // 等温壁 CV フラグ [nCells] (wall_isothermal bcond の CV=1)。node-centered 等温壁の壁ノード
    // T ピン (applyNodeIsothermalWallPin / WMLES pin) と対で、block-DPLUR のエネルギー行 (row4) を
    // decouple するのに使う (運動量 wall_flag と同型。state ピンと Jacobian の不整合による発散防止)。
    geom_int* iso_wall_flag_d = nullptr;

    // node-centered 周期境界 DOF 同一視 (median-dual M4, §4.5)。周期 partner ノードを union-find で
    // グループ化し、各 CV の root(=master) index を持つ。periodicRoot[c]==c なら root/非周期、それ以外は slave。
    // periodicNodeGather (res 和を group 全員に書く) と合併体積で「両側部分 CV を 1 CV」として扱う。
    std::vector<geom_int> periodicRoot;       // [nCells] host: 各 CV の group root
    geom_int* periodicRoot_d = nullptr;       // [nCells] device
    // 合併前の部分 CV 体積 [nCells] (node periodic 時のみ確保、他は nullptr)。
    // 体積ソース (bodyForce, ransSource の k/ω 源) は periodicNodeGather で group 合算されるため、
    // merged volume (var volume) を使うと seam で 2 倍 (コーナー group は 4 倍) になる。
    // ソース側はこの部分体積を使うこと (plans/active/turbulence-iddes-sst.md 変更ログ 2026-07-22)。
    geom_float* volumePartial_d = nullptr;
    geom_int  nPeriodicMembers = 0;           // root != self の CV 数 (>0 なら周期同一視を実行)

    mesh();
    ~mesh();
    mesh(geom_int& , geom_int& ,geom_int& , geom_int& ,
         geom_int& , geom_int& ,
         std::vector<node>& , std::vector<plane>& , std::vector<cell>& , std::vector<bcond>& );

    void readMesh(std::string);

    void setPeriodicPartner();
    void setMeshMap_d();
    // node モード: setPeriodicPartner の partnerCellID から周期ノード group(union-find) を構築し
    // periodicRoot/periodicRoot_d を埋め、各 group の合併体積を var_volume へ書き戻す (§4.5.3)。
    void buildPeriodicNodeGroups(bool nodeMode, geom_float* var_volume_d);
};

struct matrix
{
private:
    geom_int itemp;
    geom_int ic0;
    geom_int ic1;
    geom_int ip;

    std::vector<geom_int> cellPlnCounter;

public:
    std::vector<std::vector<geom_int>> structure;
    std::vector<std::vector<flow_float>> lhs;
    std::vector<flow_float> rhs;

    std::vector<std::vector<geom_int>> localPlnOfCell;

    std::vector<flow_float> row_index;
    std::vector<flow_float> col_index;
    std::vector<flow_float> value;

    matrix();

    void initMatrix(mesh& );
};
