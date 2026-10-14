#include "wallRepPoint_d.cuh"

#include <cstdio>
#include <map>
#include <vector>

#include "cuda_forge/cudaWrapper.cuh"
#include "cuda_forge/wallLaw_d.cuh"

// 代表内点と壁法線距離の作り方は wallRepPoint_d.cuh の説明のとおり
// (plans/active/architecture-float-state-double-geometry.md §4.2 3.、段 ②)。
namespace {

constexpr flow_float kSmall = kWallLawSmall;   // ransWallFunction_d / wmlesWallModel_d の kSmall (1e-12) と同じ値

// 境界面 ib ごとの代表内点。T は座標の型: 本番は double (mesh::cc64 の写し)、診断の旧版は flow_float (デバイスの ccx..ccz)。
// 式・比較・候補の除外・同率の扱いは段 ① までの compute_wall_friction_sst_d / wmles_wall_model_d の選択と同じ
// (T = flow_float ならそれと同じ式・同じ型)。n̂ はデバイスの面ベクトル (flow_float) を T に広げてから作る。
template <typename T>
__global__ void wall_rep_point_d(
    geom_int nb, const geom_int* bplane_plane, const geom_int* bplane_cell,
    const flow_float* sx, const flow_float* sy, const flow_float* sz, const flow_float* ss,
    geom_int nNormalPlanes, const geom_int* cell_planes_index, const geom_int* cell_planes, const geom_int* plane_cells,
    const geom_int* wall_flag,
    const T* ccx, const T* ccy, const T* ccz,
    geom_int* irep_out, flow_float* dn_out, flow_float* dist_out, flow_float* cos_out)
{
    const geom_int ib = blockDim.x * blockIdx.x + threadIdx.x;
    if (ib >= nb) return;

    const geom_int ip = bplane_plane[ib];
    const geom_int ic = bplane_cell[ib];

    const T kS = static_cast<T>(kSmall);

    // 壁単位法線 (外向き)
    const T sss = max(static_cast<T>(ss[ip]), kS);
    const T nx = static_cast<T>(sx[ip]) / sss;
    const T ny = static_cast<T>(sy[ip]) / sss;
    const T nz = static_cast<T>(sz[ip]) / sss;

    T best_cos = static_cast<T>(-2.0);
    geom_int bestI = -1;
    T bestDn = kS;
    T bestDist = kS;
    const T xw = ccx[ic], yw = ccy[ic], zw = ccz[ic];
    for (geom_int j = cell_planes_index[ic]; j < cell_planes_index[ic + 1]; ++j) {
        const geom_int ipn = cell_planes[j];
        if (ipn >= nNormalPlanes) continue;                 // 内部双対面のみ
        const geom_int a = plane_cells[2 * ipn + 0];
        const geom_int b = plane_cells[2 * ipn + 1];
        const geom_int cand = (a == ic) ? b : a;
        if (wall_flag[cand] != 0) continue;                 // 内部ノードのみ
        const T dx = ccx[cand] - xw;
        const T dy = ccy[cand] - yw;
        const T dz = ccz[cand] - zw;
        const T dist = sqrt(dx * dx + dy * dy + dz * dz);
        if (dist <= kS) continue;
        const T dn_in = -(dx * nx + dy * ny + dz * nz);     // 壁内向き距離 (>0 が内部側)
        if (dn_in <= static_cast<T>(0.0)) continue;
        const T cosv = dn_in / dist;                        // 内向き法線との cos
        if (cosv > best_cos) { best_cos = cosv; bestI = cand; bestDn = dn_in; bestDist = dist; }
    }
    irep_out[ib] = bestI;
    dn_out[ib]   = static_cast<flow_float>(bestDn);
    dist_out[ib] = static_cast<flow_float>(bestDist);
    cos_out[ib]  = static_cast<flow_float>(best_cos);
}

bool isWallKind(const bcond& bc) { return bc.bcondKind == "wall" || bc.bcondKind == "wall_isothermal"; }

template <typename T>
void launchRep(mesh& msh, variables& var, const bcond& bc, const T* cx, const T* cy, const T* cz,
               geom_int* irep, flow_float* dn, flow_float* dist, flow_float* cosv)
{
    const geom_int nb = static_cast<geom_int>(bc.iPlanes.size());
    if (nb == 0) return;
    const int bs = 256;
    const int ng = (int)((nb + bs - 1) / bs);
    wall_rep_point_d<T><<<ng, bs>>>(
        nb, bc.map_bplane_plane_d, bc.map_bplane_cell_d,
        var.p_d.at("sx"), var.p_d.at("sy"), var.p_d.at("sz"), var.p_d.at("ss"),
        msh.nNormalPlanes, msh.map_cell_planes_index_d, msh.map_cell_planes_d, msh.map_plane_cells_d,
        msh.wall_flag_d, cx, cy, cz, irep, dn, dist, cosv);
    gpuErrchk( cudaPeekAtLastError() );
}

// double の値の位置 (mesh::cc64、[3*i+k]) を x・y・z の 3 本にしてデバイスへ上げる (作り終えたら解放する)。
struct Dev64 { double* x = nullptr; double* y = nullptr; double* z = nullptr; };
Dev64 uploadCc64(const mesh& msh)
{
    const size_t n = (size_t)msh.nCells_all;
    std::vector<double> h[3];
    for (int k = 0; k < 3; ++k) {
        h[k].resize(n);
        for (size_t i = 0; i < n; ++i) h[k][i] = msh.cc64[3*i + k];
    }
    Dev64 d;
    double** p[3] = {&d.x, &d.y, &d.z};
    for (int k = 0; k < 3; ++k) {
        gpuErrchk( cudaMalloc((void**)p[k], n*sizeof(double)) );
        gpuErrchk( cudaMemcpy(*p[k], h[k].data(), n*sizeof(double), cudaMemcpyHostToDevice) );
    }
    return d;
}
void freeDev64(Dev64& d)
{
    cudaFree(d.x); cudaFree(d.y); cudaFree(d.z);
    d = Dev64{};
}

// 1 つの bcond の代表点を、座標の出どころを選んでデバイスの書き先へ作る。
void buildOne(mesh& msh, variables& var, const bcond& bc, bool use64, const Dev64& c64,
              geom_int* irep, flow_float* dn, flow_float* dist, flow_float* cosv)
{
    if (use64) launchRep<double>(msh, var, bc, c64.x, c64.y, c64.z, irep, dn, dist, cosv);
    else       launchRep<flow_float>(msh, var, bc, var.c_d.at("ccx"), var.c_d.at("ccy"), var.c_d.at("ccz"), irep, dn, dist, cosv);
}

std::map<const bcond*, WallRepPoints> g_rep;
bool g_built = false;

}  // namespace

const WallRepPoints& wallRepPoints(cudaConfig& cuda_cfg, mesh& msh, variables& var, const bcond& bc)
{
    (void)cuda_cfg;
    static const WallRepPoints empty{};
    if (!g_built) {
        g_built = true;
        if (msh.wall_flag_d != nullptr) {
            const bool use64 = msh.hasGeom64();
            Dev64 c64;
            if (use64) c64 = uploadCc64(msh);
            int nBc = 0;
            long long nFaces = 0, nNoCand = 0;
            for (const bcond& b : msh.bconds) {
                if (!isWallKind(b) || b.iPlanes.empty()) continue;
                WallRepPoints w;
                w.nb = static_cast<geom_int>(b.iPlanes.size());
                gpuErrchk( cudaMalloc((void**)&w.irep_d, w.nb*sizeof(geom_int)) );
                gpuErrchk( cudaMalloc((void**)&w.dn_d,   w.nb*sizeof(flow_float)) );
                gpuErrchk( cudaMalloc((void**)&w.dist_d, w.nb*sizeof(flow_float)) );
                gpuErrchk( cudaMalloc((void**)&w.cos_d,  w.nb*sizeof(flow_float)) );
                buildOne(msh, var, b, use64, c64, w.irep_d, w.dn_d, w.dist_d, w.cos_d);
                gpuErrchkKernelSync();
                std::vector<geom_int> h(w.nb);
                gpuErrchk( cudaMemcpy(h.data(), w.irep_d, w.nb*sizeof(geom_int), cudaMemcpyDeviceToHost) );
                for (geom_int v : h) if (v < 0) ++nNoCand;
                nFaces += w.nb;
                ++nBc;
                g_rep[&b] = w;
            }
            if (use64) freeDev64(c64);
            printf("[wallRepPoint] wall representative points built once for %d wall bconds, %lld faces "
                   "(%lld without an interior candidate), from %s positions\n",
                   nBc, nFaces, nNoCand, use64 ? "double" : "flow_float");
        }
    }
    const auto it = g_rep.find(&bc);
    return (it != g_rep.end()) ? it->second : empty;
}

void wallRepPointsHost(cudaConfig& cuda_cfg, mesh& msh, variables& var, const bcond& bc, bool legacy,
                       std::vector<geom_int>& irep, std::vector<flow_float>& dn,
                       std::vector<flow_float>& dist, std::vector<flow_float>& cosv)
{
    (void)cuda_cfg;
    irep.clear(); dn.clear(); dist.clear(); cosv.clear();
    const geom_int nb = static_cast<geom_int>(bc.iPlanes.size());
    if (nb == 0 || msh.wall_flag_d == nullptr) return;
    const bool use64 = !legacy && msh.hasGeom64();
    Dev64 c64;
    if (use64) c64 = uploadCc64(msh);
    geom_int* i_d = nullptr;
    flow_float *dn_d = nullptr, *dist_d = nullptr, *cos_d = nullptr;
    gpuErrchk( cudaMalloc((void**)&i_d,    nb*sizeof(geom_int)) );
    gpuErrchk( cudaMalloc((void**)&dn_d,   nb*sizeof(flow_float)) );
    gpuErrchk( cudaMalloc((void**)&dist_d, nb*sizeof(flow_float)) );
    gpuErrchk( cudaMalloc((void**)&cos_d,  nb*sizeof(flow_float)) );
    buildOne(msh, var, bc, use64, c64, i_d, dn_d, dist_d, cos_d);
    gpuErrchkKernelSync();
    irep.resize(nb); dn.resize(nb); dist.resize(nb); cosv.resize(nb);
    gpuErrchk( cudaMemcpy(irep.data(), i_d,    nb*sizeof(geom_int),   cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(dn.data(),   dn_d,   nb*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(dist.data(), dist_d, nb*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    gpuErrchk( cudaMemcpy(cosv.data(), cos_d,  nb*sizeof(flow_float), cudaMemcpyDeviceToHost) );
    cudaFree(i_d); cudaFree(dn_d); cudaFree(dist_d); cudaFree(cos_d);
    if (use64) freeDev64(c64);
}
