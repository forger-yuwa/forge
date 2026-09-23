// qAccumulator_d.cu — FP64 影アキュムレータの初期化・カウンタ操作
//   (plans/active/time_integration-fp64-accumulator.md §4.3 / §5.1 S1)
#include <cuda_runtime.h>

#include "cudaWrapper.cuh"
#include "qAccumulator.hpp"

__global__ void qaccInitFromQ_d(double* q0, double* q1, double* q2, double* q3, double* q4,
                                const flow_float* p0, const flow_float* p1, const flow_float* p2,
                                const flow_float* p3, const flow_float* p4, geom_int nCells)
{
    const geom_int ic = blockDim.x * blockIdx.x + threadIdx.x;
    if (ic >= nCells) return;
    q0[ic] = (double)p0[ic];
    q1[ic] = (double)p1[ic];
    q2[ic] = (double)p2[ic];
    q3[ic] = (double)p3[ic];
    q4[ic] = (double)p4[ic];
}

void qaccInitFromQ(double* const qacc[5], flow_float* const q[5], geom_int nCells)
{
    const int block = 256;
    const int grid  = (int)((nCells + block - 1) / block);
    if (grid <= 0) return;
    qaccInitFromQ_d<<<grid, block>>>(qacc[0], qacc[1], qacc[2], qacc[3], qacc[4],
                                     q[0], q[1], q[2], q[3], q[4], nCells);
    gpuErrchk( cudaPeekAtLastError() );
}

void qaccResetAdoptCounter(int* adopt_d)
{
    if (adopt_d != nullptr) gpuErrchk( cudaMemset(adopt_d, 0, sizeof(int)) );
}

int qaccReadAdoptCounter(const int* adopt_d)
{
    if (adopt_d == nullptr) return 0;
    int n = 0;
    gpuErrchk( cudaMemcpy(&n, adopt_d, sizeof(int), cudaMemcpyDeviceToHost) );
    return n;
}
