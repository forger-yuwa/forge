// 二相拡散 既定化 plan §5.1 #4s: G1(ii) の判別 A/B (指数スケール)。
// 保存した面入力 (tp_faces.h5 → export_faces.py の binary) から、本番の tp_face_flux<float> と double 参照を GPU で評価する。
// A = 液 (rg0, rg1) をそのまま、B = 2^k 倍 (k はコマンド引数、事前登録 64)。出力 Jl (float/double) を binary で書く。
// ビルドは本番と同じ CUDA フラグで (build/CMakeFiles/.../flags.make の CUDA_FLAGS を写す)。
//   nvcc <本番の CUDA_FLAGS> -std=c++17 -I . -o tp_jl_scale_ab tests/unit/tp_jl_scale_ab.cu
#include "cuda_forge/twoPhaseDiffusion_d.cuh"
#include "cuda_forge/twoPhaseFaceDiag_d.cuh"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

// 入力: 面ごとに float の列 (順序は export_faces.py と同じ)
struct Cols { int n; int iw; };
enum { C_RHO0, C_RHO1, C_RG0, C_RG1, C_F, C_GEO, C_GEOABS, C_CT, C_L, C_FIXED };   // + 種ごと rY0,rY1,D,h (4n) + Q (6)

__global__ void eval_d(int N, int n, int iw, int ncol, const float* X, double scale, float* Jf, double* Jd)
{
    const int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= N) return;
    const float* x = X + (size_t)i * ncol;
    TpFaceIn in; in.n = n; in.iw = iw;
    in.rho0 = x[C_RHO0]; in.rho1 = x[C_RHO1];
    in.rg0 = (float)((double)x[C_RG0] * scale); in.rg1 = (float)((double)x[C_RG1] * scale);   // 2^k 倍は正確 (溢れない範囲)
    in.f = x[C_F]; in.geo = x[C_GEO]; in.geo_abs = x[C_GEOABS]; in.ct = x[C_CT]; in.L = x[C_L];
    int c = C_FIXED;
    for (int s = 0; s < n; ++s) { in.rY0[s] = x[c++]; }
    for (int s = 0; s < n; ++s) { in.rY1[s] = x[c++]; }
    for (int s = 0; s < n; ++s) { in.D[s] = x[c++]; }
    for (int s = 0; s < n; ++s) { in.h[s] = x[c++]; }
    for (int m = 0; m < TP_NQ; ++m) in.rQ0[m] = x[c++];
    for (int m = 0; m < TP_NQ; ++m) in.rQ1[m] = x[c++];
    TpFaceOut o; tp_face_flux(in, o);
    TpFaceInT<double> dd; tp_face_in_to_double(in, dd);
    TpFaceOutT<double> od; tp_face_flux(dd, od);
    Jf[i] = o.Jl; Jd[i] = od.Jl;
}

int main(int argc, char** argv)
{
    if (argc < 4) { std::fprintf(stderr, "usage: %s faces.bin k out.bin\n", argv[0]); return 2; }
    FILE* fp = std::fopen(argv[1], "rb"); if (!fp) return 2;
    int hdr[4]; if (std::fread(hdr, sizeof(int), 4, fp) != 4) return 2;   // N, n, iw, ncol
    const int N = hdr[0], n = hdr[1], iw = hdr[2], ncol = hdr[3];
    std::vector<float> X((size_t)N * ncol);
    if (std::fread(X.data(), sizeof(float), X.size(), fp) != X.size()) return 2;
    std::fclose(fp);
    const int k = std::atoi(argv[2]);
    const double scale = std::ldexp(1.0, k);
    float* dX; float* dJf; double* dJd;
    cudaMalloc(&dX, X.size() * sizeof(float)); cudaMalloc(&dJf, N * sizeof(float)); cudaMalloc(&dJd, N * sizeof(double));
    cudaMemcpy(dX, X.data(), X.size() * sizeof(float), cudaMemcpyHostToDevice);
    eval_d<<<(N + 255) / 256, 256>>>(N, n, iw, ncol, dX, scale, dJf, dJd);
    cudaDeviceSynchronize();
    std::vector<float> Jf(N); std::vector<double> Jd(N);
    cudaMemcpy(Jf.data(), dJf, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(Jd.data(), dJd, N * sizeof(double), cudaMemcpyDeviceToHost);
    FILE* fo = std::fopen(argv[3], "wb"); if (!fo) return 2;
    std::fwrite(Jf.data(), sizeof(float), N, fo); std::fwrite(Jd.data(), sizeof(double), N, fo); std::fclose(fo);
    std::printf("evaluated %d faces, k %d, err %s\n", N, k, cudaGetErrorString(cudaGetLastError()));
    return 0;
}
