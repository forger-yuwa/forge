#include <vector>
#include <list>

#include "flowFormat.hpp"
#include "mesh/mesh.hpp"
#include "point_probes.hpp"


point_probes::point_probes() {};

point_probes::~point_probes() {
    //for (auto& valName : valNames)
    //{
    //    cudaWrapper::cudaFree_wrapper(this->var[valName]);
    //}
}

void point_probes::allocVariables(const int &useGPU , mesh& msh)
{
    for (auto& valName : valNames)
    {
        //this->c[cellValName].resize(msh.nCells);
        this->var[valName].resize(this->nProbes); // including ghost cells

        gpuErrchk( cudaMalloc((void**) &(this->var_d[valName]), (this->nProbes)*sizeof(flow_float)) );
    }
}

void point_probes::copyVariables_probes_H2D_all()
{
    for (auto& name : this->valNames)
    {
        cudaWrapper::cudaMemcpy_H2D_wrapper(this->var[name].data() , this->var_d[name], this->var[name].size());
    }
}

void point_probes::copyVariables_probes_D2H_all()
{
    for (auto& name : this->valNames)
    {
        cudaWrapper::cudaMemcpy_D2H_wrapper(this->var_d[name], this->var[name].data() , this->var[name].size());
    }
}

geom_float norm2_(const Point &l, const Point &r) {
  const geom_float dx = l.x - r.x;
  const geom_float dy = l.y - r.y;
  const geom_float dz = l.z - r.z;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

static void useKDTree(const std::vector<Point> &left, const std::vector<Point> &right, std::vector<geom_float> &distance, std::vector<geom_int> &index) {
    kdtree *tree = kd_create(3);

    std::unique_ptr<int[]> indexes(new int[left.size()]);
    for (int i = 0; i < (int)left.size(); ++i) {
        indexes[i] = i;
        kd_insert3(tree, left[i].x, left[i].y, left[i].z, &indexes[i]);
    }

    int l, r;
    geom_float minimum = std::numeric_limits<geom_float>::max();
    for (int j = 0; j < (int)right.size(); ++j) {
        kdres *set = kd_nearest3(tree, right[j].x, right[j].y, right[j].z);

        int i = *(int *)kd_res_item_data(set);
        kd_res_free(set);
        geom_float d = norm2_(left[i], right[j]);

        distance.push_back(d);
        index.push_back(d);
    }

    kd_free(tree);
    //return minimum;
}

void point_probes::setNearestCell(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh)
{
    if (this->nProbes < 1) return ;

    std::vector<Point> cellPoints;
    std::vector<geom_float> distance; 
    std::vector<geom_int> index; 

    for (auto& cell : msh.cells) {
        geom_float x = cell.centCoords[0];
        geom_float y = cell.centCoords[1];
        geom_float z = cell.centCoords[2];
        cellPoints.emplace_back(x,y,z);
    } 

    if (this->nProbes > 0) {

        useKDTree(cellPoints, inputProbeXYZ, distance, index);
        std::copy(index.begin(), index.end(), this->cell_id.begin());
    } else {
        distance.resize(this->nProbes);
        index.resize(this->nProbes);
        geom_float zero = 0.0;
        geom_int mone = -1;
        std::fill(distance.begin(), distance.end(), zero);
        std::fill(index.begin(), index.end(), mone);
    }

    cudaWrapper::cudaMemcpy_H2D_wrapper(this->cell_id.data() , this->cell_id_d, this->cell_id.size());
}

void point_probes::outputProbes(solverConfig& cfg, mesh& msh , variables& var)
{
    if (this->nProbes == 0) return;

    for (int i = 0; i<this->nProbes; i++) {
        std::fstream outfile; // read write
        std::ostringstream oss;
        oss << i;
        std::string fname = "probe_" + oss.str() + ".out";
        outfile.open(this->outFileName);

        outfile << cfg.nStepOuter << " , " << cfg.totalTime;

        for (auto varName : this->valNames) {
            outfile << " , " << this->var[varName][i] ;
        }
        outfile << "\n";

        outfile.close();
    }
}