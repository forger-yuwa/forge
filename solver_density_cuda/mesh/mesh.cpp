#include "mesh.hpp"
#include "cuda_forge/cudaWrapper.cuh"

using namespace std;
using namespace HighFive;

node::node() {};
node::node(geom_float &x, geom_float &y, geom_float &z)
{
    this->coords.push_back(x);
    this->coords.push_back(y);
    this->coords.push_back(z);
}


cell::cell() : regionId(0) {}
cell::cell(vector<geom_int> &iNodes) 
{
    this->iNodes = iNodes;
}

bcond::bcond() {}

bcond::~bcond() 
{
    for (auto& item : this->bvar_d)
    {
        //gpuErrchk( cudaFree(item.second) );
        cudaWrapper::cudaFree_wrapper(item.second);
    }
    cudaWrapper::cudaFree_wrapper(this->map_bplane_plane_d);
    cudaWrapper::cudaFree_wrapper(this->map_bplane_cell_d);
    if (this->Yb_d != nullptr) cudaFree(this->Yb_d);
}
bcond::bcond(const geom_int& physID, const vector<geom_int>& iPlanes, 
             const vector<geom_int>& iCells , const vector<geom_int>& iBPlanes)
{
    this->physID   = physID;
    this->iPlanes  = iPlanes;
    this->iCells   = iCells;
    this->iBPlanes = iBPlanes;
}

void bcond::bcondInitVariables(const int &useGPU)
{
    // ------------------------------------------------------
    // *** allocate and init variables for boundary cells ***
    // ------------------------------------------------------
    // float
    for (auto& bValName : this->bplaneValNames)
    {
        if (valueTypes.find(bValName) == valueTypes.end()) { // not neccesally variable for this bc
            continue;

        } else { // needed
            this->bvar[bValName].resize(this->iPlanes.size());

            int type = valueTypes[bValName];

            if (type == 1) { // read uniform float from yaml
                //cout << "read " << bValName << " of " << bcondKind  << " from config\n";
                for (flow_float& var : this->bvar[bValName])
                {
                    var = this->inputFloats[bValName];
                }
            }

            //cout << "alloc " << bValName << " of " << bcondKind  << " from config\n";
            gpuErrchk( cudaMalloc(&(this->bvar_d[bValName]) , this->iPlanes.size()*sizeof(flow_float)) );
            type = valueTypes[bValName];

            if (type == 1) { // read uniform float from yaml
                gpuErrchk( cudaMemcpy(this->bvar_d[bValName] , &(this->bvar[bValName][0]) ,
                                             (size_t)(this->iPlanes.size()*sizeof(flow_float)), 
                                             cudaMemcpyHostToDevice) );
            }
            
       }
    }

    // ----------------------------------------------------------
    // *** allocate and init INT variables for boundary cells ***
    // ----------------------------------------------------------
    for (auto& bIntName : this->bplaneIntNames)
    {

        if (valueTypes.find(bIntName) == valueTypes.end()) { // not found
            continue;

        } else { // found
            this->bint[bIntName].resize(this->iPlanes.size());
            int type = valueTypes[bIntName];

            if (type == 1) { // read uniform int from yaml
                cout << "read uniform" << bIntName << " of " << bcondKind  << " from config\n";
                for (geom_int& var : this->bint[bIntName])
                {
                    var = this->inputFloats[bIntName];
                }
            }

            if (useGPU == 1){
                gpuErrchk( cudaMalloc(&(this->bint_d[bIntName]) , this->iPlanes.size()*sizeof(flow_float)) );

                int type = valueTypes[bIntName];

                if (type == 1) { // read uniform int from yaml
                    gpuErrchk( cudaMemcpy(this->bint_d[bIntName] , &(this->bint[bIntName][0]) ,
                                                 (size_t)(this->iPlanes.size()*sizeof(flow_float)), 
                                                 cudaMemcpyHostToDevice) );
                } // TODO: set special treatment
            }
        }

    }
}

void bcond::copyVariables_bplane_D2H()
{
    // float
    for (auto& bValName : this->bplaneValNames)
    {
        if (valueTypes.find(bValName) == valueTypes.end()) { // not neccesally variable for this bc
            continue;

        } else { // needed
            cudaWrapper::cudaMemcpy_D2H_wrapper(this->bvar_d[bValName], this->bvar[bValName].data(), this->bvar[bValName].size());
        }
    }
}

void bcond::output_preparation(std::vector<node>& nodes, std::vector<plane>& planes)
{
    if (this->output_preparation_flg == 1) return;

    vector<geom_int> nodes_bc_flg;

    nodes_bc_flg.resize(nodes.size());
    this->inodes_g2l.resize(nodes.size());
    std::fill(nodes_bc_flg.begin(), nodes_bc_flg.end(), 0);
    std::fill(inodes_g2l.begin(), inodes_g2l.end(), -1);

    geom_int ipl = 0;
    for (auto& ip : this->iPlanes) {
        this->planes_local.push_back(planes[ip]);

        for (auto& in : planes[ip].iNodes) {
            //cout << "inl=" << inl << endl;
            nodes_bc_flg[in] = 1;
            //this->inodes_g2l[in] = inl;
            //this->inodes_l2g.push_back(in);
            //inl++;
        }
        ipl++;
    }

    geom_int inl = -1;
    for (geom_int ing=0 ; ing<nodes.size(); ing++) {
        if (nodes_bc_flg[ing] == 1){
            inl++;
            this->inodes_g2l[ing] = inl;
            this->inodes_l2g.push_back(ing);
        }
    }

    this->output_preparation_flg = 1;

}

mesh::mesh(){}
mesh::~mesh()
{
    cudaWrapper::cudaFree_wrapper(this->map_plane_cells_d);
}

mesh::mesh(geom_int& nNodes,geom_int& nPlanes,geom_int& nCells, geom_int& nNormalPlanes, 
    geom_int& nBPlanes, geom_int& nBconds,
    vector<node> &nodes , vector<plane> &planes , vector<cell>& cells , vector<bcond>& bconds)
{
    this->nNodes = nNodes;
    this->nPlanes = nPlanes;
    this->nCells = nCells;
    this->nNormalPlanes= nNormalPlanes;
    this->nBPlanes = nBPlanes;
    this->nBconds = nBconds;
    this->nodes = nodes;
    this->planes = planes;
    this->cells = cells;
    this->bconds = bconds;
}

void mesh::readMesh(string fname)
{
    File file(fname, File::ReadOnly);

    // read basic 
    Group group = file.getGroup("/MESH");

    Attribute a = group.getAttribute("nNodes");
    this->nNodes = a.read<geom_int>();

    a = group.getAttribute("nCells");
    this->nCells = a.read<geom_int>();
    cout << "Number of Cells: " << this->nCells << endl;

    a = group.getAttribute("nPlanes");
    this->nPlanes = a.read<geom_int>();
    cout << "Number of Planes: " << this->nPlanes << endl;

    a = group.getAttribute("nNormalPlanes");
    this->nNormalPlanes = a.read<geom_int>();
    cout << "Number of Normal Planes: " << this->nNormalPlanes << endl;

    a = group.getAttribute("nBPlanes");
    this->nBPlanes= a.read<geom_int>();
    cout << "Number of Boundary Planes: " << this->nBPlanes << endl;

    a = group.getAttribute("nBconds");
    this->nBconds = a.read<geom_int>();
    cout << "Number of Boundary Conditions: " << this->nBconds << endl;

//ghst>
    this->nCells_ghst = this->nBPlanes;
    this->nCells_all  = this->nCells + this->nCells_ghst;

    cout << "Number of Ghost Cells: " << this->nCells_ghst << endl;
    cout << "Number of All   Cells: " << this->nCells_all << endl;
    //this->nCells_all  = this->nCells ;
//ghst<

    // nodes
    this->nodes.resize(this->nNodes);
    std::vector<geom_float> coord;
    file.getDataSet("/MESH/COORD").read(coord);

    geom_int ii = 0;
    for (geom_int i=0; i<(this->nNodes); i++)
    {
        nodes[i].coords.resize(3);
        nodes[i].coords[0] = coord[ii];
        ii += 1;
        nodes[i].coords[1] = coord[ii];
        ii += 1;
        nodes[i].coords[2] = coord[ii];
        ii += 1;
    }
    
    // planes
    this->planes.resize(this->nPlanes);
    std::vector<geom_int> strct;
    std::vector<geom_float> surfVect;
    std::vector<geom_float> surfArea;
    std::vector<geom_float> centCoords;
    file.getDataSet("/PLANES/STRUCT").read(strct);
    file.getDataSet("/PLANES/surfVect").read(surfVect);
    file.getDataSet("/PLANES/surfArea").read(surfArea);
    file.getDataSet("/PLANES/centCoords").read(centCoords);

    geom_int ipp = 0;
    for (geom_int ip=0; ip<this->nPlanes; ip++)
    {
        geom_int nn = strct[ipp];
        this->planes[ip].iNodes.resize(nn);
        ipp += 1;
        for (geom_int in=0; in<nn; in++)
        {
            this->planes[ip].iNodes[in] = strct[ipp];
            ipp += 1;
        }

        nn = strct[ipp];
        this->planes[ip].iCells.resize(nn);
        ipp += 1;
        for (geom_int in=0; in<nn; in++)
        {
            this->planes[ip].iCells[in] = strct[ipp];
            ipp += 1;
        }

        this->planes[ip].surfVect.resize(3);
        this->planes[ip].surfVect[0] = surfVect[3*ip + 0];
        this->planes[ip].surfVect[1] = surfVect[3*ip + 1];
        this->planes[ip].surfVect[2] = surfVect[3*ip + 2];

        this->planes[ip].surfArea    = surfArea[ip];

        this->planes[ip].centCoords.resize(3);
        this->planes[ip].centCoords[0] = centCoords[3*ip + 0];
        this->planes[ip].centCoords[1] = centCoords[3*ip + 1];
        this->planes[ip].centCoords[2] = centCoords[3*ip + 2];
    }
 
//ghst>
    // Cells including ghost cells 
//    this->cells.resize(this->nCells);
    this->cells.resize(this->nCells_all);
//ghst<
    std::vector<geom_int> strct2;
    file.getDataSet("/CELLS/STRUCT").read(strct2);

    std::vector<geom_float> volume;
    std::vector<geom_float> centCoords2;
    file.getDataSet("/CELLS/volume").read(volume);
    file.getDataSet("/CELLS/centCoords").read(centCoords2);

    std::vector<geom_int> regionIds;
    if (file.exist("/CELLS/regionId")) {
        file.getDataSet("/CELLS/regionId").read(regionIds);
    }

    geom_int icc = 0;
    for (geom_int ic=0; ic<this->nCells; ic++)
    {
        geom_int nn = strct2[icc];
        this->cells[ic].iNodes.resize(nn);
        icc += 1;
        for (geom_int in=0; in<nn; in++)
        {
            this->cells[ic].iNodes[in] = strct2[icc];
            icc += 1;
        }

        nn = strct2[icc];
        this->cells[ic].iPlanes.resize(nn);
        icc += 1;
        for (geom_int in=0; in<nn; in++)
        {
            this->cells[ic].iPlanes[in] = strct2[icc];
            icc += 1;
        }

        nn = strct2[icc];
        this->cells[ic].iPlanesDir.resize(nn);
        icc += 1;
        for (geom_int in=0; in<nn; in++)
        {
            this->cells[ic].iPlanesDir[in] = strct2[icc];
            icc += 1;
        }
        this->cells[ic].ieleType = strct2[icc];
        icc += 1;

        this->cells[ic].volume = volume[ic];

        this->cells[ic].centCoords.resize(3);
        this->cells[ic].centCoords[0] = centCoords2[3*ic + 0];
        this->cells[ic].centCoords[1] = centCoords2[3*ic + 1];
        this->cells[ic].centCoords[2] = centCoords2[3*ic + 2];

        this->cells[ic].regionId = regionIds.empty() ? 0 : regionIds[ic];
    }

    // nodeValueAtNode: 値の位置をノード座標に統一 (node モード, CV ic == ノード ic)。双対重心の半径は rEff へ。
    // ゴースト生成 (下) より前に行い、ゴースト鏡映もノード基準にする。
    if (this->nodeValueAtNode == 1 && (geom_int)this->nodes.size() >= this->nCells) {
        this->rEff.assign(this->nCells, (geom_float)0.0);
        geom_int nswap = 0; geom_float maxShift = 0.0;
        for (geom_int ic = 0; ic < this->nCells; ic++) {
            if (this->nodes[ic].coords.size() < 3) continue;
            this->rEff[ic] = this->cells[ic].centCoords[1];
            const geom_float dx = this->cells[ic].centCoords[0] - this->nodes[ic].coords[0];
            const geom_float dy = this->cells[ic].centCoords[1] - this->nodes[ic].coords[1];
            const geom_float dz = this->cells[ic].centCoords[2] - this->nodes[ic].coords[2];
            maxShift = std::max(maxShift, std::sqrt(dx*dx + dy*dy + dz*dz));
            this->cells[ic].centCoords[0] = this->nodes[ic].coords[0];
            this->cells[ic].centCoords[1] = this->nodes[ic].coords[1];
            this->cells[ic].centCoords[2] = this->nodes[ic].coords[2];
            ++nswap;
        }
        cout << "[mesh] nodeValueAtNode: centCoords <- node coords for " << nswap
             << " CVs (max centroid shift " << maxShift << "), dual-centroid r kept in rEff" << endl;
    }

    // boundary conditions
    Group grp = file.getGroup("/BCONDS");
    geom_int nb = grp.getNumberObjects();
    bconds.resize(nb);

    geom_int ib = 0;
    geom_int nGhost = 0;
    for (string oname : grp.listObjectNames())
    {
        this->bconds[ib].physID = stoi(oname);
        grp = file.getGroup("/BCONDS/"+oname);

        a = grp.getAttribute("bcondKind");
        this->bconds[ib].bcondKind = a.read<std::string>();
        cout << "in mesh.cpp  physID=" << this->bconds[ib].physID << endl;
        cout << "in mesh.cpp  bcondKind=" << this->bconds[ib].bcondKind << endl;

        // iCells
        std::vector<geom_int> iCells;
        grp.getDataSet("iCells").read(iCells);

        this->bconds[ib].iCells.resize(iCells.size());
        for (geom_int ic = 0 ; ic<iCells.size() ; ic++)
        {
            this->bconds[ib].iCells[ic] = iCells[ic];
        }

        // iPlanes
        std::vector<geom_int> iPlanes;
        grp.getDataSet("iPlanes").read(iPlanes);

        this->bconds[ib].iPlanes.resize(iPlanes.size());
        for (geom_int ip = 0 ; ip<iPlanes.size() ; ip++)
        {
            this->bconds[ib].iPlanes[ip] = iPlanes[ip];
        }
        // node モードの壁 Dirichlet 等で境界 plane を持たない bcond (iPlanes 空) を許容する。
        if (!iPlanes.empty())
            cout << "               ip min=" << iPlanes[0] << ", ip max=" << iPlanes[iPlanes.size()-1]<< endl;
        else
            cout << "               (no boundary planes; e.g. node-mode wall Dirichlet)" << endl;

        // node モード壁可視化用 primal 境界面接続 (あれば復元)。
        this->bconds[ib].vizBfaceNodes.clear();
        if (grp.exist("vizBfaceSizes")) {
            std::vector<geom_int> bfNodes, bfSizes;
            grp.getDataSet("vizBfaceNodes").read(bfNodes);
            grp.getDataSet("vizBfaceSizes").read(bfSizes);
            geom_int off = 0;
            for (geom_int s : bfSizes) {
                this->bconds[ib].vizBfaceNodes.emplace_back(bfNodes.begin()+off, bfNodes.begin()+off+s);
                off += s;
            }
        }

        // iBPlanes
        std::vector<geom_int> iBPlanes;
        grp.getDataSet("iBPlanes").read(iBPlanes);

        this->bconds[ib].iBPlanes.resize(iBPlanes.size());
        for (geom_int ibp = 0 ; ibp<iBPlanes.size() ; ibp++)
        {
            this->bconds[ib].iBPlanes[ibp] = iBPlanes[ibp];
        }

//ghst>
        // add ghost cells
        for (geom_int ibl = 0 ; ibl<iBPlanes.size() ; ibl++)
        {
            geom_int ip = iPlanes[ibl];
            geom_int ic;

            // change plane information
            this->bconds[ib].iCells_ghst.resize(iBPlanes.size());
            this->bconds[ib].iCells_ghst[ibl] = nCells+nGhost;
            if (this->planes[ip].iCells.size() == 1) {
                this->planes[ip].iCells.resize(2);
                this->planes[ip].iCells[1] =nCells+nGhost;

                ic = this->planes[ip].iCells[0];
            } else {
                std::cerr << "Error: something wrong to add ghost cells" << endl;
                exit(3);
            }

            // make ghost cell
            this->cells[nCells+nGhost].iPlanes.resize(1);
            this->cells[nCells+nGhost].iPlanes[0] = ip;

            this->cells[nCells+nGhost].ieleType = this->cells[ic].ieleType;
            this->cells[nCells+nGhost].volume   = this->cells[ic].volume;
            this->cells[nCells+nGhost].centCoords.resize(3);

            geom_float xc = this->cells[ic].centCoords[0];
            geom_float yc = this->cells[ic].centCoords[1];
            geom_float zc = this->cells[ic].centCoords[2];

            geom_float xp = this->planes[ip].centCoords[0];
            geom_float yp = this->planes[ip].centCoords[1];
            geom_float zp = this->planes[ip].centCoords[2];

            geom_float dx = xp - xc;
            geom_float dy = yp - yc;
            geom_float dz = zp - zc;

            geom_float ss = this->planes[ip].surfArea;
            geom_float nx = this->planes[ip].surfVect[0]/ss;
            geom_float ny = this->planes[ip].surfVect[1]/ss;
            geom_float nz = this->planes[ip].surfVect[2]/ss;

            geom_float dnx = (dx*nx +dy*ny + dz*nz)*nx;
            geom_float dny = (dx*nx +dy*ny + dz*nz)*ny;
            geom_float dnz = (dx*nx +dy*ny + dz*nz)*nz;
            
            this->cells[nCells+nGhost].centCoords[0] = xc + 2*dnx;
            this->cells[nCells+nGhost].centCoords[1] = yc + 2*dny;
            this->cells[nCells+nGhost].centCoords[2] = zc + 2*dnz;
            nGhost++;
        }
//ghst<
        ib++;
    }

    // node-centered 可視化トポロジ (/VIZMESH)。存在すれば node モードとして出力側が使う。
    if (file.exist("/VIZMESH"))
    {
        Group vgrp = file.getGroup("/VIZMESH");
        this->nVizCells    = vgrp.getAttribute("nVizCells").read<geom_int>();
        this->vizCONNE_dim = vgrp.getAttribute("vizCONNE_dim").read<geom_int>();
        file.getDataSet("/VIZMESH/CONNE").read(this->vizCONNE);
        cout << "readMesh: loaded /VIZMESH (node-centered viz, nVizCells=" << this->nVizCells << ")" << endl;
    }
}


void mesh::setPeriodicPartner()
{
    std::list<int> checked_bcIDs;

    for (auto& bc0 : this->bconds) {
        if (bc0.bcondKind == "periodic") {
            int bcID         = bc0.physID;
            int bcID_partner = bc0.inputInts["partnerBCID"];

            flow_float dx;
            flow_float dy;
            flow_float dz;
            flow_float dtheta;

            if (bc0.inputInts["type"] == 0) { // Cartesian
                dx = bc0.inputFloats["dx"];
                dy = bc0.inputFloats["dy"];
                dz = bc0.inputFloats["dz"];

            } else if (bc0.inputInts["type"] == 1) { // rotation
                dtheta = bc0.inputFloats["dtheta"];
            }

            if (std::find(checked_bcIDs.begin(), checked_bcIDs.end(), bcID) == checked_bcIDs.end()){ // not saved BC
                checked_bcIDs.push_back(bcID);
                checked_bcIDs.push_back(bcID_partner);

                std::vector<geom_int> map_ib1_iplane;

                // find nearest planes of partner
                for (auto& bc1 : this->bconds) {
                    if (bc1.physID != bcID_partner) continue;

                    //Eigen::VectorXd XYZ(3);
                    //Eigen::MatrixXd partnerXYZ(3, bc1.iPlanes.size());
                    //Eigen::MatrixXd::Index index;

                    vector<geom_float> XYZ(3);
                    vector<vector<geom_float>> partnerXYZ;
                    geom_float index;

                    partnerXYZ.resize(3);
                    partnerXYZ[0].resize(bc1.iPlanes.size());
                    partnerXYZ[1].resize(bc1.iPlanes.size());
                    partnerXYZ[2].resize(bc1.iPlanes.size());

                    geom_float x1;
                    geom_float y1;
                    geom_float z1;

                    geom_float x1_r1;
                    geom_float y1_r1;
                    geom_float z1_r1;

                    geom_int ib1_local = 0;
                    for (geom_int& ip1 : bc1.iPlanes) {
                        x1 = this->planes[ip1].centCoords[0];
                        y1 = this->planes[ip1].centCoords[1];
                        z1 = this->planes[ip1].centCoords[2];

                        if (bc0.inputInts["type"] == 0) {
                            x1_r1 = x1 - dx;
                            y1_r1 = y1 - dy;
                            z1_r1 = z1 - dz;

                        } else if (bc0.inputInts["type"] == 1) {
                            x1_r1 = x1;
                            y1_r1 = cos(-dtheta)*y1 -sin(-dtheta)*z1;
                            z1_r1 = sin(-dtheta)*y1 +cos(-dtheta)*z1;
                        }

                        partnerXYZ[0][ib1_local] = x1_r1;
                        partnerXYZ[1][ib1_local] = y1_r1;
                        partnerXYZ[2][ib1_local] = z1_r1;

                        map_ib1_iplane.push_back(ip1);

                        ib1_local++;
                    }

                    geom_int ib0_local = 0;
                    for (geom_int& ip0 : bc0.iPlanes) {
                        geom_float x0 = this->planes[ip0].centCoords[0];
                        geom_float y0 = this->planes[ip0].centCoords[1];
                        geom_float z0 = this->planes[ip0].centCoords[2];

                        XYZ[0] = x0;
                        XYZ[1] = y0;
                        XYZ[2] = z0;

                        geom_float dist2 = 1e+30;
                        geom_float dist2_temp;
                        geom_int ib1_local = 0;
                        for (geom_int& ip1 : bc1.iPlanes) {
                            x1 = partnerXYZ[0][ib1_local] ;
                            y1 = partnerXYZ[1][ib1_local] ;
                            z1 = partnerXYZ[2][ib1_local] ;

                            dist2_temp = pow((x1-x0),2) +pow((y1-y0),2) +pow((z1-z0),2);

                            if (dist2>dist2_temp) {
                                dist2 = dist2_temp;
                                index = ib1_local;
                            }

                            ib1_local++;
                        }

                        // don't use eigen because too many warnings
                        //(partnerXYZ.colwise() - XYZ).colwise().squaredNorm().minCoeff(&index);

                        geom_int ip1 = map_ib1_iplane[index];

                        geom_int ic0 = this->planes[ip0].iCells[0];
                        geom_int ic1 = this->planes[ip1].iCells[0];

                        bc0.bint["partnerPlnID"][ib0_local] = ip1;
                        bc1.bint["partnerPlnID"][index]     = ip0;

                        bc0.bint["partnerCellID"][ib0_local] = ic1;
                        bc1.bint["partnerCellID"][index]     = ic0;


                        ib0_local++;
                    }
                } 
            }else {
                continue; // already added bc
            }
        }
    }
};


// node-centered 周期境界 DOF 同一視 (median-dual M4, §4.5)。
// setPeriodicPartner が埋めた各 periodic bcond の partnerCellID (= partner ノード CV) を使って、
// 周期 partner ノードを union-find でグループ化する。各 group は「継ぎ目で割れた 1 つの CV の集合」で、
// res を group 全員へ足し合わせ (periodicNodeGather) + 合併体積で更新すると 1 CV として振る舞う。
// cell モード / 非周期では何もしない。
void mesh::buildPeriodicNodeGroups(bool nodeMode, geom_float* var_volume_d)
{
    this->nPeriodicMembers = 0;
    this->periodicRoot.assign(this->nCells, 0);
    for (geom_int c = 0; c < this->nCells; ++c) this->periodicRoot[c] = c;

    if (!nodeMode) return;

    // 周期 bcond が無ければ何もしない
    bool hasPeriodic = false;
    for (auto& bc : this->bconds) if (bc.bcondKind == "periodic") hasPeriodic = true;
    if (!hasPeriodic) return;

    // union-find (path halving, iterative なので std::function 不要)
    auto findRoot = [&](geom_int x) {
        while (this->periodicRoot[x] != x) {
            this->periodicRoot[x] = this->periodicRoot[this->periodicRoot[x]];
            x = this->periodicRoot[x];
        }
        return x;
    };
    auto unite = [&](geom_int a, geom_int b) {
        geom_int ra = findRoot(a), rb = findRoot(b);
        if (ra == rb) return;
        // 小さい index を root に固定 (決定的)
        if (rb < ra) { geom_int t = ra; ra = rb; rb = t; }
        this->periodicRoot[rb] = ra;
    };

    for (auto& bc : this->bconds) {
        if (bc.bcondKind != "periodic") continue;
        if (bc.bint.count("partnerCellID") == 0) {
            std::cerr << "[buildPeriodicNodeGroups] periodic bcond physID=" << bc.physID
                      << " has no partnerCellID (setPeriodicPartner not run?). skip.\n";
            continue;
        }
        const auto& partner = bc.bint["partnerCellID"];
        for (size_t l = 0; l < bc.iPlanes.size(); ++l) {
            const geom_int ip = bc.iPlanes[l];
            const geom_int node = this->planes[ip].iCells[0]; // node モードでは CV index
            if (l >= partner.size()) break;
            const geom_int part = partner[l];
            if (node < 0 || node >= this->nCells || part < 0 || part >= this->nCells) continue;
            unite(node, part);
        }
    }

    // 全 root を平坦化 + member 数カウント
    for (geom_int c = 0; c < this->nCells; ++c) this->periodicRoot[c] = findRoot(c);
    for (geom_int c = 0; c < this->nCells; ++c) if (this->periodicRoot[c] != c) this->nPeriodicMembers++;

    // ★ 合併前の部分体積を退避 (体積ソースの二重計上防止用)。gather が group 合算するため、
    // ソース項は部分体積で加算しないと seam で 2 倍 (x∩z コーナー group は 4 倍) になる
    // (case/39 z-seam の +75% 速度異常で発覚)。
    {
        std::vector<geom_float> partVol(this->nCells);
        for (geom_int c = 0; c < this->nCells; ++c) partVol[c] = this->cells[c].volume;
        gpuErrchk(cudaMalloc((void **)&(this->volumePartial_d), sizeof(geom_float)*this->nCells));
        gpuErrchk(cudaMemcpy(this->volumePartial_d, partVol.data(),
                             sizeof(geom_float)*this->nCells, cudaMemcpyHostToDevice));
    }

    // 合併体積: group ごとに体積を合算し、group 全員の volume をその合算値にする
    // (両者が同 res・同 vol で更新され bit 一致同期するため。詳細は plans §4.5.3)。
    std::vector<geom_float> groupVol(this->nCells, 0.0);
    for (geom_int c = 0; c < this->nCells; ++c) groupVol[this->periodicRoot[c]] += this->cells[c].volume;
    std::vector<geom_float> newVol(this->nCells, 0.0);
    for (geom_int c = 0; c < this->nCells; ++c) {
        newVol[c] = groupVol[this->periodicRoot[c]];
        this->cells[c].volume = newVol[c]; // host も整合させる
    }
    if (var_volume_d != nullptr) {
        gpuErrchk(cudaMemcpy(var_volume_d, newVol.data(), sizeof(geom_float)*this->nCells, cudaMemcpyHostToDevice));
    }

    // device へ root をアップロード
    gpuErrchk(cudaMalloc((void **)&(this->periodicRoot_d), sizeof(geom_int)*this->nCells));
    gpuErrchk(cudaMemcpy(this->periodicRoot_d, this->periodicRoot.data(),
                         sizeof(geom_int)*this->nCells, cudaMemcpyHostToDevice));

    std::cout << "[buildPeriodicNodeGroups] node periodic DOF identification: "
              << this->nPeriodicMembers << " slave CVs merged into masters (union-find).\n";
}


void mesh::setMeshMap_d()
{
    gpuErrchk(cudaMalloc((void **)&(this->map_plane_cells_d), sizeof(geom_int)*this->nPlanes*2));

    // Build a list of plane indices to be processed by convective-flux kernels.
    // Layout:
    //   [0, nNormalPlanes)              : interior (normal) planes
    //   [nNormalPlanes, ...) periodic   : periodic boundary planes
    //   [..., n_normal_halo_planes)     : non-periodic boundary planes (slip/wall/inlet/outlet/...)
    // For non-periodic boundary planes, plane_cells[2*ip+1] is the ghost-cell index,
    // so the convective-flux kernel can use the ghost-cell state directly without
    // a dedicated boundary flux kernel.
    geom_int n_normal_halo_planes =  this->nNormalPlanes;

    for (auto& bc : this->bconds)
    {
        n_normal_halo_planes += bc.iPlanes.size();
    }

    this->nNormal_halo_Planes = n_normal_halo_planes;


    gpuErrchk(cudaMalloc((void **)&(this->normal_halo_planes_d), sizeof(geom_int)*n_normal_halo_planes));

    geom_int* normal_halo_planes;
    normal_halo_planes = (geom_int *)malloc(sizeof(geom_int)*n_normal_halo_planes);


    geom_int ip_sum = 0;
    for (geom_int ip=0; ip<this->nNormalPlanes; ip++ ) {
        normal_halo_planes[ip] = ip_sum;
        ip_sum += 1;

    }

     for (auto& bc : this->bconds)
    {
        if (bc.bcondKind == "periodic") {
            for (auto& ip : bc.iPlanes){
                normal_halo_planes[ip_sum] = ip;
                ip_sum += 1;
            }
        }
    }

    // 非 periodic 境界 plane: 非壁を先に、壁 (wall/wall_isothermal) を末尾に並べる。
    // これにより node 弱形式では主対流ループを「内部+非壁境界」に限定 (末尾の壁 plane を除外) でき、
    // 壁だけ別途 pressure-only の境界 flux (convectiveFlux_boundary_d) で扱える。
    // cell モードは従来どおり全 plane (壁含む) を主ループでゴースト処理する (順序は無影響)。
    auto isWallKind = [](const std::string& k){ return k == "wall" || k == "wall_isothermal"; };
    this->nBoundaryHaloPlanes = 0;  // 非 periodic 境界 plane 総数 (非壁 + 壁)
    for (auto& bc : this->bconds)
    {
        if (bc.bcondKind != "periodic" && !isWallKind(bc.bcondKind)) {
            for (auto& ip : bc.iPlanes){
                normal_halo_planes[ip_sum] = ip;
                ip_sum += 1;
                this->nBoundaryHaloPlanes += 1;
            }
        }
    }
    this->nWallHaloPlanes = 0;
    for (auto& bc : this->bconds)
    {
        if (isWallKind(bc.bcondKind)) {
            for (auto& ip : bc.iPlanes){
                normal_halo_planes[ip_sum] = ip;
                ip_sum += 1;
                this->nWallHaloPlanes += 1;
                this->nBoundaryHaloPlanes += 1;
            }
        }
    }

    gpuErrchk(cudaMemcpy(this->normal_halo_planes_d  , normal_halo_planes ,
                         sizeof(geom_int)*n_normal_halo_planes , cudaMemcpyHostToDevice));

    free(normal_halo_planes); 


    geom_int* pc_h;
    geom_int* bp_h;
    geom_int* bc_h;
    geom_int* bcg_h;
    pc_h = (geom_int *)malloc(sizeof(geom_int)*this->nPlanes*2);

    for (geom_int ip=0; ip<this->nPlanes; ip++)
    {
        pc_h[2*ip + 0] = this->planes[ip].iCells[0]; 
        pc_h[2*ip + 1] = this->planes[ip].iCells[1]; 
        //printf("ip=%d, ic1=%d, ic2=%d\n", ip, pc_h[2*ip + 0], pc_h[2*ip + 1]);
    }

    gpuErrchk(cudaMemcpy(this->map_plane_cells_d  , pc_h , 
                     sizeof(geom_int)*(this->nPlanes*2) , cudaMemcpyHostToDevice));

    free(pc_h); 


    // cell -> planes
    gpuErrchk(cudaMalloc((void **)&(this->map_cell_planes_index_d), sizeof(geom_int)*(this->nCells+1)));

    geom_int* cpi_h;
    cpi_h = (geom_int *)malloc(sizeof(geom_int)*(this->nCells+1));

    cpi_h[0] = 0; 
    for (geom_int ic=0; ic<this->nCells; ic++)
    {
        geom_int nplane_local = this->cells[ic].iPlanes.size();
        cpi_h[ic+1] = cpi_h[ic] + nplane_local; 
    }

    gpuErrchk(cudaMalloc((void **)&(this->map_cell_planes_d), sizeof(geom_int)*cpi_h[this->nCells]));
    geom_int* cp_h;
    cp_h = (geom_int *)malloc(sizeof(geom_int)*cpi_h[this->nCells]);

    geom_int ipln = 0;
    for (geom_int ic=0; ic<this->nCells; ic++)
    {
        for (auto ip :this->cells[ic].iPlanes) {
            cp_h[ipln] = ip; 
            ipln += 1;
        }
    }

    gpuErrchk(cudaMemcpy(this->map_cell_planes_index_d  , cpi_h , 
                         sizeof(geom_int)*(this->nCells+1) , cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(this->map_cell_planes_d  , cp_h , 
                     sizeof(geom_int)*cpi_h[this->nCells] , cudaMemcpyHostToDevice));

    free(cpi_h); 
    free(cp_h); 


    for (bcond& bc : this->bconds)
    {
        bp_h = (geom_int *)malloc(sizeof(geom_int)*bc.iPlanes.size());
        bc_h = (geom_int *)malloc(sizeof(geom_int)*bc.iPlanes.size());
        bcg_h = (geom_int *)malloc(sizeof(geom_int)*bc.iPlanes.size()); // ghost cell
        gpuErrchk(cudaMalloc((void **)&(bc.map_bplane_plane_d), sizeof(geom_int)*bc.iPlanes.size()));
        gpuErrchk(cudaMalloc((void **)&(bc.map_bplane_cell_d) , sizeof(geom_int)*bc.iPlanes.size()));
        gpuErrchk(cudaMalloc((void **)&(bc.map_bplane_cell_ghst_d) , sizeof(geom_int)*bc.iPlanes.size()));

        for (geom_int ibl=0 ; ibl<bc.iPlanes.size() ; ibl++)
        {
            bp_h[ibl]  = bc.iPlanes[ibl];
            bc_h[ibl]  = bc.iCells[ibl];
            bcg_h[ibl] = bc.iCells_ghst[ibl];
        }

        gpuErrchk(cudaMemcpy(bc.map_bplane_plane_d , bp_h , 
                             sizeof(geom_int)*(bc.iPlanes.size()) , cudaMemcpyHostToDevice));

        gpuErrchk(cudaMemcpy(bc.map_bplane_cell_d , bc_h , 
                             sizeof(geom_int)*(bc.iPlanes.size()) , cudaMemcpyHostToDevice));

        gpuErrchk(cudaMemcpy(bc.map_bplane_cell_ghst_d , bcg_h , 
                             sizeof(geom_int)*(bc.iPlanes.size()) , cudaMemcpyHostToDevice));
 
        free(bp_h);
        free(bc_h);
        free(bcg_h);
    }

    // 境界隣接 CV フラグ: いずれかの bcond の境界 CV (iCells) を 1 にする。
    {
        std::vector<geom_int> bnf(this->nCells, 0);
        for (bcond& bc : this->bconds)
            for (geom_int ic : bc.iCells)
                if (ic >= 0 && ic < this->nCells) bnf[ic] = 1;
        gpuErrchk(cudaMalloc((void **)&(this->bnode_flag_d), sizeof(geom_int)*this->nCells));
        gpuErrchk(cudaMemcpy(this->bnode_flag_d, bnf.data(),
                             sizeof(geom_int)*this->nCells, cudaMemcpyHostToDevice));
    }

    // 軸上 CV フラグ: ノード R(=coords[1])≈0 の CV を 1 にする (node-centered 軸対称用)。
    // node モードでは CV ic == ノード ic。cell モードでは無意味だが使用側 (axisym+node) で限定する。
    {
        std::vector<geom_int> axf(this->nCells, 0);
        const geom_int nlim = std::min((geom_int)this->nodes.size(), this->nCells);
        for (geom_int ic = 0; ic < nlim; ++ic)
            if (this->nodes[ic].coords.size() >= 2 &&
                std::abs(this->nodes[ic].coords[1]) < (geom_float)1.0e-9) axf[ic] = 1;
        gpuErrchk(cudaMalloc((void **)&(this->axis_flag_d), sizeof(geom_int)*this->nCells));
        gpuErrchk(cudaMemcpy(this->axis_flag_d, axf.data(),
                             sizeof(geom_int)*this->nCells, cudaMemcpyHostToDevice));

    }

    // 壁 CV フラグ: wall 種別 bcond の境界 CV (iCells) を 1 にする (node-centered 壁 Dirichlet 用)。
    // node モードはマルチマーカ emit なので、コーナー (壁∩出口 等) は壁・他境界の双方の iCells に出現する
    // (gmshReader buildMedianDual, ow=ib)。ここで wall bcond の iCells を走査すれば、そのコーナーも
    // 壁ノードとして確実にフラグされる (壁優先で no-slip を効かせる)。壁 bcond も自分の半割面 ghost を持つ。
    {
        std::vector<geom_int> wf(this->nCells, 0);
        for (bcond& bc : this->bconds)
            if (bc.bcondKind == "wall" || bc.bcondKind == "wall_isothermal")
                for (geom_int ic : bc.iCells)
                    if (ic >= 0 && ic < this->nCells) wf[ic] = 1;
        gpuErrchk(cudaMalloc((void **)&(this->wall_flag_d), sizeof(geom_int)*this->nCells));
        gpuErrchk(cudaMemcpy(this->wall_flag_d, wf.data(),
                             sizeof(geom_int)*this->nCells, cudaMemcpyHostToDevice));

        // 等温壁ノードフラグ (エネルギー行 decouple 用)。壁∩等温壁コーナーは等温側を優先 (T ピン対象)。
        std::vector<geom_int> iwf(this->nCells, 0);
        for (bcond& bc : this->bconds)
            if (bc.bcondKind == "wall_isothermal")
                for (geom_int ic : bc.iCells)
                    if (ic >= 0 && ic < this->nCells) iwf[ic] = 1;
        gpuErrchk(cudaMalloc((void **)&(this->iso_wall_flag_d), sizeof(geom_int)*this->nCells));
        gpuErrchk(cudaMemcpy(this->iso_wall_flag_d, iwf.data(),
                             sizeof(geom_int)*this->nCells, cudaMemcpyHostToDevice));
    }
};

matrix::matrix(){}
void matrix::initMatrix(mesh& msh)
{
    structure.resize(msh.nCells);
    lhs.resize(msh.nCells);
    rhs.resize(msh.nCells);

    cellPlnCounter.resize(msh.nCells);
    localPlnOfCell.resize(msh.nNormalPlanes);

    for (ic0=0; ic0<msh.nCells; ic0++)
    {
        structure[ic0].push_back(ic0);
    }

    for (auto& i : cellPlnCounter)
    {
        i = 1;
    }

    for (ip=0; ip<msh.nNormalPlanes; ip++)
    {
        ic0 = msh.planes[ip].iCells[0];
        ic1 = msh.planes[ip].iCells[1];
        structure[ic0].push_back(ic1);
        structure[ic1].push_back(ic0);

        localPlnOfCell[ip].push_back(cellPlnCounter[ic0]);
        localPlnOfCell[ip].push_back(cellPlnCounter[ic1]);

        cellPlnCounter[ic0] += 1;
        cellPlnCounter[ic1] += 1;
    }

    for (geom_int ist=0; ist<structure.size(); ist++)
    {
        lhs[ist].resize(structure[ist].size());
    }

}


// =============================================================================
// line-implicit のライン構築 (plans/active/time_integration-line-implicit.md)。
// 壁 CV を種に「最初は最近接隣 (壁法線 = 最短エッジ)、以降は前進方向と最も揃う隣」へ
// greedy に鎖を伸ばす。構造化 TFI/積層メッシュでは j 列がそのまま再構成される。
// 載らない CV は line_prev/next = -1 のまま point-DPLUR fallback。
void mesh::buildImplicitLines(const flow_float* ccx, const flow_float* ccy, const flow_float* ccz)
{
    const geom_int n = this->nCells;
    // 内部面の隣接リスト (両側とも実 CV の面のみ)
    std::vector<std::vector<geom_int>> adj(n);
    for (geom_int ip = 0; ip < this->nNormalPlanes; ++ip) {
        const geom_int c0 = this->planes[ip].iCells[0];
        const geom_int c1 = this->planes[ip].iCells[1];
        if (c0 >= 0 && c0 < n && c1 >= 0 && c1 < n) {
            adj[c0].push_back(c1);
            adj[c1].push_back(c0);
        }
    }
    // 壁種 bcond の CV を種にする
    std::vector<geom_int> seeds;
    {
        std::vector<char> isWall(n, 0);
        for (bcond& bc : this->bconds)
            if (bc.bcondKind == "wall" || bc.bcondKind == "wall_isothermal")
                for (geom_int ic : bc.iCells)
                    if (ic >= 0 && ic < n && !isWall[ic]) { isWall[ic] = 1; seeds.push_back(ic); }
    }
    // 座標: node モードは nodes[ic].coords (CV=節点) を正とする。host の ccx 配列は
    // 初期化順によって未充填 (全ゼロ) のことがある (2026-09-02 実測) ため信用しない。
    std::vector<double> px(n), py(n), pz(n);
    const bool useNodes = ((geom_int)this->nodes.size() >= n);
    for (geom_int ic = 0; ic < n; ++ic) {
        if (useNodes && this->nodes[ic].coords.size() >= 2) {
            px[ic] = (double)this->nodes[ic].coords[0];
            py[ic] = (double)this->nodes[ic].coords[1];
            pz[ic] = (this->nodes[ic].coords.size() >= 3) ? (double)this->nodes[ic].coords[2] : 0.0;
        } else {
            px[ic] = (double)ccx[ic]; py[ic] = (double)ccy[ic]; pz[ic] = (double)ccz[ic];
        }
    }
    auto dist2 = [&](geom_int a, geom_int b) {
        const double dx = px[a] - px[b];
        const double dy = py[a] - py[b];
        const double dz = pz[a] - pz[b];
        return dx*dx + dy*dy + dz*dz;
    };
    std::vector<char> visited(n, 0);
    std::vector<geom_int> prevArr(n, -1), nextArr(n, -1);
    std::vector<geom_int> offsets(1, 0), cells;
    const double cosMin = 0.7;   // 前進方向との整列条件 (45°)
    // 実験用: ライン長上限 (FORGE_LINE_MAXLEN, 既定 0 = 無制限)。壁側から数えて打ち切る。
    long maxLen_env = 0;
    if (const char* e = getenv("FORGE_LINE_MAXLEN")) maxLen_env = atol(e);
    geom_int nLines = 0, lenMax = 0;
    for (geom_int seed : seeds) {
        if (visited[seed]) continue;
        std::vector<geom_int> line;
        line.push_back(seed);
        visited[seed] = 1;
        geom_int cur = seed;
        double dx = 0.0, dy = 0.0, dz = 0.0;   // 前進方向 (未定 = 最初)
        bool haveDir = false;
        while (true) {
            geom_int best = -1; double bestScore = -1.0e300;
            for (geom_int nb : adj[cur]) {
                if (visited[nb]) continue;
                if (!haveDir) {
                    const double s = -dist2(cur, nb);       // 最近接 (壁法線 = 最短エッジ)
                    if (s > bestScore) { bestScore = s; best = nb; }
                } else {
                    const double ex = px[nb] - px[cur];
                    const double ey = py[nb] - py[cur];
                    const double ez = pz[nb] - pz[cur];
                    const double el = std::sqrt(ex*ex + ey*ey + ez*ez) + 1.0e-300;
                    const double c  = (ex*dx + ey*dy + ez*dz) / el;
                    if (c > cosMin && c > bestScore) { bestScore = c; best = nb; }
                }
            }
            if (best < 0) break;
            const double ex = px[best] - px[cur];
            const double ey = py[best] - py[cur];
            const double ez = pz[best] - pz[cur];
            const double el = std::sqrt(ex*ex + ey*ey + ez*ez) + 1.0e-300;
            dx = ex/el; dy = ey/el; dz = ez/el; haveDir = true;
            line.push_back(best);
            visited[best] = 1;
            cur = best;
            if (maxLen_env > 0 && (long)line.size() >= maxLen_env) break;
        }
        if ((geom_int)line.size() < 3) {          // 短すぎる鎖はライン化しない (point fallback)
            for (geom_int ic : line) visited[ic] = 0;
            visited[seed] = 1;                    // 種は再訪しない
            continue;
        }
        for (size_t k = 0; k < line.size(); ++k) {
            if (k > 0)               prevArr[line[k]] = line[k-1];
            if (k + 1 < line.size()) nextArr[line[k]] = line[k+1];
            cells.push_back(line[k]);
        }
        offsets.push_back((geom_int)cells.size());
        lenMax = std::max(lenMax, (geom_int)line.size());
        ++nLines;
    }
    this->nImplicitLines = nLines;
    const geom_int nOn = (geom_int)cells.size();
    printf("[lineImplicit] lines=%d  covered CVs=%d/%d (%.1f%%)  maxLen=%d\n",
           (int)nLines, (int)nOn, (int)n, 100.0*(double)nOn/(double)std::max(n,(geom_int)1), (int)lenMax);
    gpuErrchk(cudaMalloc((void**)&this->line_prev_d, sizeof(geom_int)*n));
    gpuErrchk(cudaMalloc((void**)&this->line_next_d, sizeof(geom_int)*n));
    gpuErrchk(cudaMemcpy(this->line_prev_d, prevArr.data(), sizeof(geom_int)*n, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(this->line_next_d, nextArr.data(), sizeof(geom_int)*n, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMalloc((void**)&this->line_offsets_d, sizeof(geom_int)*offsets.size()));
    gpuErrchk(cudaMemcpy(this->line_offsets_d, offsets.data(), sizeof(geom_int)*offsets.size(), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMalloc((void**)&this->line_cells_d, sizeof(geom_int)*std::max(nOn,(geom_int)1)));
    gpuErrchk(cudaMemcpy(this->line_cells_d, cells.data(), sizeof(geom_int)*nOn, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMalloc((void**)&this->line_Kprev_d, sizeof(flow_float)*25*this->nCells_all));
    gpuErrchk(cudaMalloc((void**)&this->line_Knext_d, sizeof(flow_float)*25*this->nCells_all));
    gpuErrchk(cudaMemset(this->line_Kprev_d, 0, sizeof(flow_float)*25*this->nCells_all));
    gpuErrchk(cudaMemset(this->line_Knext_d, 0, sizeof(flow_float)*25*this->nCells_all));
    gpuErrchk(cudaMalloc((void**)&this->line_W_d, sizeof(double)*25*this->nCells_all));
    gpuErrchk(cudaMalloc((void**)&this->line_y_d, sizeof(double)*5*this->nCells_all));
    // v2: factor/solve 分離用 (D̃ の LU 因子・ピボット・ライン失敗フラグ)。
    gpuErrchk(cudaMalloc((void**)&this->line_LU_d, sizeof(double)*25*this->nCells_all));
    gpuErrchk(cudaMalloc((void**)&this->line_piv_d, sizeof(signed char)*5*this->nCells_all));
    gpuErrchk(cudaMalloc((void**)&this->line_fail_d, sizeof(unsigned char)*std::max(nLines,(geom_int)1)));
    gpuErrchk(cudaMemset(this->line_fail_d, 0, sizeof(unsigned char)*std::max(nLines,(geom_int)1)));
}
