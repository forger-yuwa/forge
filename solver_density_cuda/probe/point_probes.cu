#include <vector>
#include <list>

#include "flowFormat.hpp"
#include "point_probes.cuh"


point_probes::point_probes()
    : nProbes(0), cell_id_d(nullptr), outStepStart(0), outStepInterval(1), dimGrid_pprobe(0)
{};

point_probes::~point_probes() {
}

void point_probes::init(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh)
{
    point_probes::readYAML();
    this->dimGrid_pprobe = dim3(ceil(this->nProbes/ (flow_float)cuda_cfg.blocksize));
    point_probes::allocVariables();
    point_probes::setNearestCell(cfg , cuda_cfg , msh);
}

void point_probes::allocVariables()
{
    if (this->nProbes == 0) return;

    for (auto& valName : valNames)
    {
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

void point_probes::readYAML()
{
    std::cout << "Read probe setting form probe.yaml" << std::endl;


    try {
        std::string fname = "probe.yaml";

        YAML::Node config = YAML::LoadFile(fname);
        
        this-> outStepInterval = config["outStepInterval"].as<int>();
        this-> outStepStart    = config["outStepStart"].as<int>();

        auto probes = config["points"];

        this->nProbes = 0;
        this->inputProbeXYZ.clear();

        if (!probes || !probes.IsMap()) {
            return;
        }

        std::vector<Point_r1> Points;

        for (const auto& point : probes) {
            std::string point_name = point.first.as<std::string>();
            flow_float x = point.second["x"].as<flow_float>();
            flow_float y = point.second["y"].as<flow_float>();
            flow_float z = point.second["z"].as<flow_float>();

            std::cout << point_name << ": x = " << x << ", y = " << y << ", z = " << z << std::endl;

            Points.emplace_back(x,y,z);
        }
        this->nProbes = Points.size();
        this->inputProbeXYZ = Points;


    } catch (const YAML::BadFile& e) {
        std::cerr << e.msg << std::endl;
        exit(EXIT_FAILURE);

    } catch (const YAML::ParserException& e) {
        std::cerr << e.msg << std::endl;
        exit(EXIT_FAILURE);
    }

}



geom_float norm2_(const Point_r1 &l, const Point_r1 &r) {
  const geom_float dx = l.x - r.x;
  const geom_float dy = l.y - r.y;
  const geom_float dz = l.z - r.z;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

#ifdef HAVE_KDTREE
static void useKDTree(const std::vector<Point_r1> &left, const std::vector<Point_r1> &right, std::vector<geom_float> &distance, std::vector<geom_int> &index) {
    kdtree *tree = kd_create(3);
    std::unique_ptr<int[]> indexes(new int[left.size()]);
    for (int i = 0; i < (int)left.size(); ++i) {
        indexes[i] = i;
        kd_insert3(tree, left[i].x, left[i].y, left[i].z, &indexes[i]);
    }
    for (int j = 0; j < (int)right.size(); ++j) {
        kdres *set = kd_nearest3(tree, right[j].x, right[j].y, right[j].z);
        int i = *(int *)kd_res_item_data(set);
        kd_res_free(set);
        geom_float d = norm2_(left[i], right[j]);
        distance.push_back(d);
        index.push_back(i);
    }
    kd_free(tree);
}
#endif

flow_float norm2(const Point_r1 &l, const Point_r1 &r) {
  const flow_float dx = l.x - r.x;
  const flow_float dy = l.y - r.y;
  const flow_float dz = l.z - r.z;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void bruteForce(const std::vector<Point_r1> &left, const Point_r1 &right, geom_float &distance, geom_int &index) {
    int l, r;

    flow_float minimum = std::numeric_limits<flow_float>::max();

    for (int i = 0; i < (int)left.size(); ++i) {
        flow_float d = norm2(left[i], right);
        if (minimum > d) {
            minimum = d;
            l = i;
        }
    }
    distance = minimum;
    index = l;
}

void point_probes::setNearestCell(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh)
{
    if (this->nProbes < 1) return ;

    std::vector<Point_r1> cellPoints;
    std::vector<geom_float> distance; 
    std::vector<geom_int> index; 

    for (auto& cell : msh.cells) {
        geom_float x = cell.centCoords[0];
        geom_float y = cell.centCoords[1];
        geom_float z = cell.centCoords[2];
        cellPoints.emplace_back(x,y,z);
    } 

    flow_float dist;
    geom_int ind;

    for (int j=0 ; j<this->nProbes ; j++) {
#ifdef HAVE_KDTREE
        // If kdtree is available we could batch via useKDTree; fallback keeps bruteForce for now.
        bruteForce(cellPoints, inputProbeXYZ[j], dist, ind);
#else
        bruteForce(cellPoints, inputProbeXYZ[j], dist, ind);
#endif
        // (optional) verbose output removed (cfg.verbose not defined)
        // std::cout << "probe " << j << ": dist=" << dist << " cell=" << ind << std::endl;
        index.push_back(ind);
    }
    this->cell_id = index;

    // Detailed debug print disabled (cfg.verbose not present in solverConfig)
    // for (auto d : distance) std::cout << "dist " << d << "\n";
    // for (auto id : index) std::cout << "ind=" << id << "\n";

    cudaWrapper::cudaMalloc_wrapper(&(this->cell_id_d) , this->cell_id.size());
    cudaWrapper::cudaMemcpy_H2D_wrapper(this->cell_id.data() , this->cell_id_d, this->cell_id.size());
}

void point_probes::outputProbes(solverConfig& cfg , cudaConfig& cudaCfg ,  mesh& msh , variables& var , int iStep)
{

    if (this->nProbes == 0) return;

    if (iStep%this->outStepInterval != 0 or iStep < this->outStepStart) return;

    setProbeVariables_d_wrapper(cfg , cudaCfg , msh , var);
    copyVariables_probes_D2H_all();

    for (int i = 0; i<this->nProbes; i++) {
        std::ofstream outfile; // read write
        std::ostringstream oss;
        oss << i;
        std::string fname = "point_probe_" + oss.str() + ".out";

        outfile.open(fname, std::ios::app);
        if (outfile.tellp() == 0) { 
            outfile << "Step , TotalTime";
            for (auto varName : this->valNames) {
                outfile << " , " << varName;
            }
            outfile << "\n";
        }

        outfile << iStep << " , " << cfg.totalTime;

        for (auto varName : this->valNames) {
            outfile << " , " << this->var[varName][i] ;
        }
        outfile << "\n";

        outfile.close();
    }
}

__global__ void setProbeVariables_d
( 
 int nProbes, 
 geom_int* probe_nearest_cell_id,

 // variables
 flow_float* var_c , flow_float* var_p 
)
{
    geom_int ip = blockDim.x*blockIdx.x + threadIdx.x;


    if (ip < nProbes) {
        geom_int  ic = probe_nearest_cell_id[ip];
        var_p[ip] = var_c[ic];
    }
}


void point_probes::setProbeVariables_d_wrapper
(solverConfig& cfg , cudaConfig& cuda_cfg , mesh& msh , variables& var)
{
    if (this->nProbes == 0) return;

    for (auto& varName : this->valNames) {
        std::cout << varName << std::endl;
        setProbeVariables_d<<<this->dimGrid_pprobe , cuda_cfg.blocksize>>>
        (
            this->nProbes,
            this->cell_id_d,

            var.c_d[varName] , this->var_d[varName]
        );
    }
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

