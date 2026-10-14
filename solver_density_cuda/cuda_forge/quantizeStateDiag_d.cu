// 診断 FORGE_DIAG_QUANTIZE_STATE のカーネル (plans/active/architecture-float-state-double-geometry.md §5.1 #19・§6.36、既定 off)。
// 1 本の保存量の配列 q[0..n) を、CV ごとに float32 に丸めて flow_float に戻す (quantizeStateValue)。読み書きはその配列だけ。
// 呼ぶのは quantizeStateDiag::quantize (step の末尾、量ごとに 1 回) だけで、off では起動しない。

#include "cuda_forge/quantizeStateDiag_d.cuh"
#include "cuda_forge/cudaConfig.cuh"
#include "cuda_forge/cudaWrapper.cuh"

__global__ void quantizeStateToFloat32_d(geom_int n, flow_float* q)
{
    const geom_int ic = blockDim.x*blockIdx.x + threadIdx.x;
    if (ic < n) q[ic] = quantizeStateValue(q[ic]);
}

void quantizeStateToFloat32_d_wrapper(cudaConfig& cuda_cfg, geom_int n, flow_float* q)
{
    if (n <= 0) return;
    const unsigned int block = cuda_cfg.dimBlock.x;
    const dim3 grid((unsigned int)((n + (geom_int)block - 1) / (geom_int)block));
    quantizeStateToFloat32_d<<<grid, cuda_cfg.dimBlock>>>(n, q);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchkKernelSync();
}
