// test_qacc_axis.cu — 軸ピンを FP64 の正本 (Qacc) に当てる契約の単体試験
//   (plans/active/architecture-float-state-double-geometry.md §5.1 #18 の (iv)・§6.34 の前提 1)
//   build: nvcc --expt-relaxed-constexpr -I. -o test_qacc_axis tests/unit/test_qacc_axis.cu
//   (solver_density_cuda で。float のビルド (flowFormat.hpp の既定) 専用。生のビットを uint32 で比べるため)
//
// **製品の __device__ 関数を直接駆動する**: 軸ピン axisSymmetryPinCell (cuda_forge/axisSymmetryPin_d.cuh) と
// commit / reconcile の qaccCommitCell / qaccReconcileCell (cuda_forge/qAccumulator_d.cuh)。
// 1 step の順序は製品と同じ: (壁 no-slip) → 軸ピン → commit → reconcile → Q_N = Q (updateVariablesOuter_d)。
//   - 壁 no-slip は nodeWallDirichlet_d.cu の enforceWallNoSlip_d (匿名名前空間で呼べない) と同じ式を下の wallPinRepro に写した。
//   - commit は applyBlockImplicitCorrection_d (update_d.cu) の経路と同じ: 軸の分岐は nullptr (無効)、updateGuardAlpha 0。
//     ON は 5 量とも qaccCommitCell(Qacc, Q, dq, 1)、OFF は Q = Q_N + dq。
// 腕:
//   OFF : qAccumulatorFP64 0 (軸ピンに Qacc を渡さない、commit は Q_N + dq)。
//   ON  : 契約 (軸ピンが Q_N を射影する前の roUyN から ke_N を作り Qacc_roe −= ke_N、Qacc_roUy = 0)。
//   OLD : 旧い作り (軸ピンは Qacc に当てない、commit は Qacc から)。軸対称の拒否を外しただけの状態。
//   REJ : 却下した案 (codex diagnose 2026-10-11 Major 1): 作用用の Q から運動エネルギーを取って Qacc に当て、
//         作用用の roe を (float)Qacc で上書きする。試験だけにある。
// 項目 (plan #18 (iv)): (a) rollout の反例 (b) 1 ULP 未満の残余の保持 (c) 残余ゼロで (float)Qacc が OFF の Q_N と生のビットで
// 一致 (乱数 20 万) (d) 壁→軸→commit で壁と軸に重なる節点 (e) 軸でない節点は不変 (f) dq ≠ 0 (g) 残余を持つ複数 step。
#include <cstdio>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <cuda_runtime.h>
#include "flowFormat.hpp"
#include "cuda_forge/qAccumulator_d.cuh"
#include "cuda_forge/axisSymmetryPin_d.cuh"

static_assert(sizeof(flow_float) == 4, "test_qacc_axis は float のビルド専用 (生のビットを uint32 で比べる)");

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++g_fail; printf("  FAIL: " __VA_ARGS__); printf("\n"); } } while (0)

enum Mode { OFF = 0, ON = 1, OLD = 2, REJ = 3 };
static const char* const kModeName[4] = {"OFF", "ON", "OLD", "REJ"};

// 1 節点の状態。作用用の Q・commit の基準 Q_N・FP64 の正本 Qacc。並びは ro, roUx, roUy, roUz, roe。
struct Node {
    flow_float q[5];      // 作用用の Q
    flow_float U[3];      // 速度 (ピンが 0 にする)
    flow_float qN[5];     // commit の基準 Q_N
    double     qa[5];     // FP64 の正本
    geom_int   axis, wall;
};

// 最後の step の途中の値 (観測用)
struct Snap {
    flow_float roeAfterPin;   // 軸ピンの後の作用用の roe (= 残差の評価に渡る値)
    flow_float q[5];          // commit の直後の Q
    double     qa[5];         // commit の直後の Qacc
};

// nodeWallDirichlet_d.cu の enforceWallNoSlip_d の 1 節点分と同じ式 (演算の順も同じ)。
__device__ void wallPinRepro(Node& x)
{
    if (x.wall != 1) return;
    const flow_float r = x.q[0];
    if (r > (flow_float)0.0) {
        x.q[4] -= (flow_float)0.5 * (x.q[1]*x.q[1] + x.q[2]*x.q[2] + x.q[3]*x.q[3]) / r;
    }
    x.q[1] = (flow_float)0.0; x.q[2] = (flow_float)0.0; x.q[3] = (flow_float)0.0;
    x.U[0] = (flow_float)0.0; x.U[1] = (flow_float)0.0; x.U[2] = (flow_float)0.0;
}

// 製品の軸ピン (ic = 0 として 1 節点を渡す)。
__device__ void axisPinProduct(Node& x, bool acc)
{
    axisSymmetryPinCell(0, &x.axis, &x.q[0], &x.q[2], &x.q[4], &x.U[1], &x.qN[2], &x.qN[4], &x.qN[0],
                        acc ? &x.qa[2] : nullptr, acc ? &x.qa[4] : nullptr);
}

// 却下した案 (試験だけ): 作用用の Q から ke を取り Qacc に当て、作用用の roe を (float)Qacc で上書きする。
__device__ void axisPinRejected(Node& x)
{
    if (x.axis != 1) return;
    const flow_float r = x.q[0];
    const flow_float ke = (r > (flow_float)0.0) ? (flow_float)0.5 * x.q[2] * x.q[2] / r : (flow_float)0.0;
    x.qa[4] -= (double)ke;
    x.qa[2] = 0.0;
    axisPinProduct(x, false);              // 作用用の Q と Q_N は製品と同じ射影
    x.q[4] = (flow_float)x.qa[4];          // 作用用の roe を Qacc から上書き (← 壁が除いた KE を戻す)
}

__device__ void applyPin(Node& x, int mode)
{
    if (mode == ON) axisPinProduct(x, true);
    else if (mode == REJ) axisPinRejected(x);
    else axisPinProduct(x, false);          // OFF・OLD
}

// update_d.cu の applyBlockImplicitCorrection_d と同じ経路 (軸の分岐 nullptr、updateGuardAlpha 0)。
__device__ void commit(Node& x, const flow_float* d, int mode)
{
    if (mode == OFF) {
        for (int k = 0; k < 5; ++k) x.q[k] = x.qN[k] + d[k];
    } else {
        for (int k = 0; k < 5; ++k) qaccCommitCell(x.qa[k], x.q[k], d[k], (flow_float)1.0);
    }
}

// updateVariablesOuter_d: reconcile (ON 系だけ) → Q_N = Q。戻り値は採用数。
__device__ int reconcileShift(Node& x, int mode)
{
    int n = 0;
    if (mode != OFF) for (int k = 0; k < 5; ++k) n += qaccReconcileCell(x.qa[k], x.q[k]) ? 1 : 0;
    for (int k = 0; k < 5; ++k) x.qN[k] = x.q[k];
    return n;
}

// nStep 回の外側の step。d は節点ごとに 5 つ (毎 step 同じ値)。最後の step の途中の値を snap に残す。
__global__ void k_steps(Node* nodes, Snap* snap, const flow_float* d, int n, int mode, int wallFirst, int nStep, int* nAdopt)
{
    const int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= n) return;
    Node x = nodes[i];
    int c = 0;
    for (int s = 0; s < nStep; ++s) {
        if (wallFirst) wallPinRepro(x);
        applyPin(x, mode);
        if (s == nStep - 1) snap[i].roeAfterPin = x.q[4];
        commit(x, d + 5*i, mode);
        if (s == nStep - 1) for (int k = 0; k < 5; ++k) { snap[i].q[k] = x.q[k]; snap[i].qa[k] = x.qa[k]; }
        c += reconcileShift(x, mode);
    }
    nodes[i] = x;
    if (c > 0) atomicAdd(nAdopt, c);
}

// 軸ピンだけ (c・e)
__global__ void k_pinOnly(Node* nodes, int n, int acc)
{
    const int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= n) return;
    axisPinProduct(nodes[i], acc != 0);
}

// ---- host ----
static uint32_t bits32(flow_float v) { uint32_t u; std::memcpy(&u, &v, 4); return u; }
static uint64_t bits64(double v) { uint64_t u; std::memcpy(&u, &v, 8); return u; }
static double ulpf(flow_float v) { return (double)std::nextafter(std::fabs(v), INFINITY) - (double)std::fabs(v); }

struct RunOut { std::vector<Node> nodes; std::vector<Snap> snap; int nAdopt; };

static RunOut run(const std::vector<Node>& in, const std::vector<flow_float>& d, int mode, int wallFirst, int nStep)
{
    const int n = (int)in.size();
    Node* dn; Snap* ds; flow_float* dd; int* da;
    cudaMalloc(&dn, n * sizeof(Node)); cudaMalloc(&ds, n * sizeof(Snap));
    cudaMalloc(&dd, 5 * n * sizeof(flow_float)); cudaMalloc(&da, sizeof(int));
    cudaMemcpy(dn, in.data(), n * sizeof(Node), cudaMemcpyHostToDevice);
    cudaMemcpy(dd, d.data(), 5 * n * sizeof(flow_float), cudaMemcpyHostToDevice);
    cudaMemset(ds, 0, n * sizeof(Snap)); cudaMemset(da, 0, sizeof(int));
    k_steps<<<(n + 127) / 128, 128>>>(dn, ds, dd, n, mode, wallFirst, nStep, da);
    cudaError_t e = cudaDeviceSynchronize();
    if (e != cudaSuccess) { printf("CUDA error: %s\n", cudaGetErrorString(e)); std::exit(2); }
    RunOut o; o.nodes.resize(n); o.snap.resize(n);
    cudaMemcpy(o.nodes.data(), dn, n * sizeof(Node), cudaMemcpyDeviceToHost);
    cudaMemcpy(o.snap.data(), ds, n * sizeof(Snap), cudaMemcpyDeviceToHost);
    cudaMemcpy(&o.nAdopt, da, sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(dn); cudaFree(ds); cudaFree(dd); cudaFree(da);
    return o;
}

static std::vector<Node> pinOnly(const std::vector<Node>& in, int acc)
{
    const int n = (int)in.size();
    Node* dn; cudaMalloc(&dn, n * sizeof(Node));
    cudaMemcpy(dn, in.data(), n * sizeof(Node), cudaMemcpyHostToDevice);
    k_pinOnly<<<(n + 255) / 256, 256>>>(dn, n, acc);
    cudaError_t e = cudaDeviceSynchronize();
    if (e != cudaSuccess) { printf("CUDA error: %s\n", cudaGetErrorString(e)); std::exit(2); }
    std::vector<Node> out(n);
    cudaMemcpy(out.data(), dn, n * sizeof(Node), cudaMemcpyDeviceToHost);
    cudaFree(dn);
    return out;
}

// 残余ゼロの節点 (Q = Q_N = (float)Qacc、Qacc = (double)Q)
static Node makeNode(flow_float ro, flow_float rux, flow_float ruy, flow_float ruz, flow_float roe, int axis, int wall)
{
    Node x{};
    const flow_float v[5] = {ro, rux, ruy, ruz, roe};
    for (int k = 0; k < 5; ++k) { x.q[k] = v[k]; x.qN[k] = v[k]; x.qa[k] = (double)v[k]; }
    x.U[0] = rux / ro; x.U[1] = ruy / ro; x.U[2] = ruz / ro;
    x.axis = axis; x.wall = wall;
    return x;
}

static bool sameNode(const Node& a, const Node& b)
{
    for (int k = 0; k < 5; ++k) {
        if (bits32(a.q[k]) != bits32(b.q[k]) || bits32(a.qN[k]) != bits32(b.qN[k]) || bits64(a.qa[k]) != bits64(b.qa[k]))
            return false;
    }
    for (int k = 0; k < 3; ++k) if (bits32(a.U[k]) != bits32(b.U[k])) return false;
    return a.axis == b.axis && a.wall == b.wall;
}

struct Lcg {
    unsigned s;
    explicit Lcg(unsigned seed) : s(seed) {}
    double u() { s = s * 1664525u + 1013904223u; return (double)(s >> 8) / 16777216.0; }   // [0,1)
    double logU(double lo, double hi) { return std::pow(10.0, lo + (hi - lo) * u()); }
    double sgn() { return u() < 0.5 ? -1.0 : 1.0; }
};

int main()
{
    const std::vector<flow_float> d0(5, (flow_float)0.0);

    // ---- (a) rollout の反例: ρ=1, ρu_r=0.125, ρE=2, dq=0 → 軸ピン → commit ----
    {
        printf("(a) rollout の反例: ρ=1, ρu_r=0.125, ρE=2, dq=0, 軸の節点 (壁でない)\n");
        const std::vector<Node> in{makeNode(1.0f, 0.0f, 0.125f, 0.0f, 2.0f, 1, 0)};
        for (int m : {OFF, ON, OLD}) {
            const RunOut o = run(in, d0, m, 0, 1);
            const Snap& s = o.snap[0];
            printf("    %-3s commit 後 Q: ρu_r=%.9g ρE=%.9g   Qacc: ρu_r=%.9g ρE=%.9g   採用 %d\n",
                   kModeName[m], (double)s.q[2], (double)s.q[4], m == OFF ? NAN : s.qa[2], m == OFF ? NAN : s.qa[4], o.nAdopt);
            if (m == OFF || m == ON) {
                CHECK(s.q[2] == 0.0f && s.q[4] == 1.9921875f, "(a) %s: commit 後の Q が (0, 1.9921875) でない", kModeName[m]);
            }
            if (m == ON) {
                CHECK(s.qa[2] == 0.0 && s.qa[4] == 1.9921875, "(a) ON: commit 後の Qacc が (0, 1.9921875) でない");
                CHECK(o.nAdopt == 0, "(a) ON: reconcile が採用した (%d)", o.nAdopt);
            }
            if (m == OLD) {
                // 旧い作りでは Qacc に当たっていないので commit がピンを巻き戻す (反例の再現)
                CHECK(s.q[2] == 0.125f && s.q[4] == 2.0f, "(a) OLD: 反例が再現しない (Q が 0.125, 2 に戻らない)");
                CHECK(s.qa[2] == 0.125 && s.qa[4] == 2.0, "(a) OLD: Qacc が 0.125, 2 でない");
            }
        }
    }

    // ---- (b) 1 ULP 未満の残余はピンの後も Qacc に残る ----
    {
        printf("\n(b) 1 ULP 未満の残余: Qacc_ρE = 2 + 1e-9 (Q_ρE = 2.0f, 1 ULP = %.3e)、ρ=1, ρu_r=0.125, dq=0\n", ulpf(2.0f));
        Node x = makeNode(1.0f, 0.0f, 0.125f, 0.0f, 2.0f, 1, 0);
        x.qa[4] = 2.0 + 1e-9;
        const double resid0 = x.qa[4] - 2.0;
        const std::vector<Node> pinned = pinOnly({x}, 1);
        const double expectPin = (2.0 + 1e-9) - (double)0.0078125f;
        const double residPin = pinned[0].qa[4] - (double)(flow_float)pinned[0].qa[4];
        printf("    ピン後: Qacc_ρE = %.17g (期待 %.17g)、残余 Qacc − (double)(float)Qacc = %.6e (ピン前 %.6e)\n",
               pinned[0].qa[4], expectPin, residPin, resid0);
        CHECK(bits64(pinned[0].qa[4]) == bits64(expectPin), "(b) ピン後の Qacc_ρE が (2+1e-9) − ke_N と一致しない");
        CHECK(residPin == resid0, "(b) ピンで残余が変わった (%.6e → %.6e)", resid0, residPin);
        CHECK(pinned[0].qa[2] == 0.0, "(b) ピン後の Qacc_ρu_r が 0 でない");
        const RunOut o = run({x}, d0, ON, 0, 1);
        const double residEnd = o.nodes[0].qa[4] - (double)o.nodes[0].q[4];
        printf("    step の末尾 (commit・reconcile の後): Q_ρE = %.9g、残余 %.6e、採用 %d\n",
               (double)o.nodes[0].q[4], residEnd, o.nAdopt);
        CHECK(o.nodes[0].q[4] == 1.9921875f, "(b) commit 後の Q_ρE が 1.9921875 でない");
        CHECK(residEnd == resid0 && o.nAdopt == 0, "(b) step の末尾で残余が残っていない (%.6e、採用 %d)", residEnd, o.nAdopt);
    }

    // ---- (c) 残余ゼロ・乱数 20 万: ピン後の (float)Qacc が OFF のピン後の Q_N と生のビットで一致 ----
    {
        const int M = 200000;
        Lcg r(20261011u);
        std::vector<Node> in(M);
        int nGuard = 0;
        for (int i = 0; i < M; ++i) {
            Node& x = in[i];
            x = Node{};
            flow_float rN = (flow_float)r.logU(-3.0, 3.0);
            if (i % 100 == 0) { rN = (i % 200 == 0) ? 0.0f : -rN; ++nGuard; }   // 密度のガード (rN ≤ 0) も混ぜる
            const flow_float uyN = (flow_float)(r.sgn() * r.logU(-8.0, 3.0));
            const double keApprox = (rN > 0.0f) ? 0.5 * (double)uyN * uyN / rN : 0.0;
            const flow_float eN = (flow_float)(keApprox * (0.5 + 2.5 * r.u()) + r.logU(-3.0, 6.0));
            const flow_float qN[5] = {rN, (flow_float)(r.sgn() * r.logU(-3.0, 3.0)), uyN,
                                      (flow_float)(r.sgn() * r.logU(-3.0, 3.0)), eN};
            for (int k = 0; k < 5; ++k) { x.qN[k] = qN[k]; x.qa[k] = (double)qN[k]; }
            // 作用用の Q は Q_N と別の値 (ピンの作用用の側が Qacc に関係しないことも見る)
            x.q[0] = (flow_float)r.logU(-3.0, 3.0);
            x.q[1] = (flow_float)(r.sgn() * r.logU(-3.0, 3.0));
            x.q[2] = (flow_float)(r.sgn() * r.logU(-8.0, 3.0));
            x.q[3] = (flow_float)(r.sgn() * r.logU(-3.0, 3.0));
            x.q[4] = (flow_float)r.logU(-3.0, 6.0);
            x.U[0] = 1.0f; x.U[1] = 2.0f; x.U[2] = 3.0f;
            x.axis = 1; x.wall = 0;
        }
        const std::vector<Node> on = pinOnly(in, 1), off = pinOnly(in, 0);
        long long nRoe = 0, nRuy = 0, nFloat = 0, nOther = 0;
        for (int i = 0; i < M; ++i) {
            if (bits32((flow_float)on[i].qa[4]) != bits32(off[i].qN[4])) ++nRoe;
            if (bits32((flow_float)on[i].qa[2]) != bits32(off[i].qN[2])) ++nRuy;
            for (int k = 0; k < 5; ++k) {
                if (bits32(on[i].q[k]) != bits32(off[i].q[k]) || bits32(on[i].qN[k]) != bits32(off[i].qN[k])) { ++nFloat; break; }
            }
            if (bits64(on[i].qa[0]) != bits64(in[i].qa[0]) || bits64(on[i].qa[1]) != bits64(in[i].qa[1])
                || bits64(on[i].qa[3]) != bits64(in[i].qa[3])) ++nOther;
        }
        printf("\n(c) 残余ゼロ・乱数 %d 節点 (うち密度のガード rN ≤ 0 が %d): ピン後の生のビットの不一致\n", M, nGuard);
        printf("    (float)Qacc_ρE vs OFF の roeN: %lld、(float)Qacc_ρu_r vs OFF の roUyN: %lld、"
               "float の側 (作用用の Q・Q_N) の ON/OFF: %lld、Qacc の ρ・ρu_x・ρu_z の変化: %lld\n", nRoe, nRuy, nFloat, nOther);
        CHECK(nRoe == 0, "(c) (float)Qacc_ρE が OFF の roeN とビット一致しない (%lld)", nRoe);
        CHECK(nRuy == 0, "(c) (float)Qacc_ρu_r が OFF の roUyN とビット一致しない (%lld)", nRuy);
        CHECK(nFloat == 0, "(c) ON で float の側が OFF と違う (%lld)", nFloat);
        CHECK(nOther == 0, "(c) ピンが Qacc の ρ・ρu_x・ρu_z を変えた (%lld)", nOther);
    }

    // ---- (d) 壁→軸→commit: 壁と軸に重なる節点 (codex の CPU 反例) ----
    {
        printf("\n(d) 壁と軸に重なる節点: ρ=1, ρu_x=1, ρu_r=0.125, ρE=2, dq=0、壁 no-slip → 軸ピン → commit\n");
        const std::vector<Node> in{makeNode(1.0f, 1.0f, 0.125f, 0.0f, 2.0f, 1, 1)};
        for (int m : {OFF, ON, REJ}) {
            const RunOut o = run(in, d0, m, 1, 1);
            const Snap& s = o.snap[0];
            printf("    %-3s 軸ピン後の作用用の ρE = %.9g、commit 後 Q: ρu_x=%.9g ρu_r=%.9g ρE=%.9g\n",
                   kModeName[m], (double)s.roeAfterPin, (double)s.q[1], (double)s.q[2], (double)s.q[4]);
            if (m == OFF || m == ON) {
                CHECK(s.roeAfterPin == 1.4921875f, "(d) %s: 軸ピン後の作用用の ρE が 1.4921875 でない", kModeName[m]);
                CHECK(s.q[4] == 1.9921875f, "(d) %s: commit 後の ρE が 1.9921875 でない", kModeName[m]);
                CHECK(s.q[1] == 1.0f && s.q[2] == 0.0f, "(d) %s: commit 後の ρu_x・ρu_r が (1, 0) でない", kModeName[m]);
            }
            if (m == ON) CHECK(s.qa[4] == 1.9921875 && s.qa[2] == 0.0, "(d) ON: commit 後の Qacc が (0, 1.9921875) でない");
            if (m == REJ) {
                // 却下した案は壁が除いた KE を戻し、軸の KE も引かない (反例の再現)
                CHECK(s.roeAfterPin == 2.0f && s.q[4] == 2.0f, "(d) REJ: 反例が再現しない (2.0 にならない)");
            }
        }
    }

    // ---- (e) 軸でない節点は ON・OFF とも不変 ----
    {
        const int M = 200000;
        Lcg r(777u);
        std::vector<Node> in(M);
        for (int i = 0; i < M; ++i) {
            Node& x = in[i];
            x = Node{};
            for (int k = 0; k < 5; ++k) {
                x.q[k]  = (flow_float)(r.sgn() * r.logU(-3.0, 3.0));
                x.qN[k] = (flow_float)(r.sgn() * r.logU(-3.0, 3.0));
                x.qa[k] = (double)x.qN[k] + r.u() * 1e-12;   // 残余つき
            }
            x.q[0] = std::fabs(x.q[0]); x.qN[0] = std::fabs(x.qN[0]);
            x.U[0] = 1.0f; x.U[1] = 2.0f; x.U[2] = 3.0f;
            x.axis = (i % 3 == 0) ? 0 : ((i % 3 == 1) ? -1 : 2);   // 1 以外は軸でない
            x.wall = i % 2;
        }
        const std::vector<Node> on = pinOnly(in, 1), off = pinOnly(in, 0);
        long long nOn = 0, nOff = 0;
        for (int i = 0; i < M; ++i) { nOn += !sameNode(on[i], in[i]); nOff += !sameNode(off[i], in[i]); }
        printf("\n(e) 軸でない節点 %d (axis_flag ∈ {0, -1, 2}): ピンで変わった節点 ON %lld、OFF %lld\n", M, nOn, nOff);
        CHECK(nOn == 0 && nOff == 0, "(e) 軸でない節点が変わった (ON %lld、OFF %lld)", nOn, nOff);
    }

    // ---- (f) dq ≠ 0 ----
    {
        printf("\n(f) dq ≠ 0\n");
        // 決まった 1 節点
        {
            const std::vector<Node> in{makeNode(1.0f, 0.0f, 0.125f, 0.0f, 2.0f, 1, 0)};
            const std::vector<flow_float> d{1e-3f, 2e-3f, 3e-3f, -1e-3f, -4e-3f};
            const RunOut on = run(in, d, ON, 0, 1), off = run(in, d, OFF, 0, 1);
            printf("    1 節点 (ρ=1, ρu_r=0.125, ρE=2, dq=(1e-3,2e-3,3e-3,-1e-3,-4e-3)): ON ρu_r=%.9g ρE=%.9g、OFF ρu_r=%.9g ρE=%.9g\n",
                   (double)on.snap[0].q[2], (double)on.snap[0].q[4], (double)off.snap[0].q[2], (double)off.snap[0].q[4]);
            bool same = true;
            for (int k = 0; k < 5; ++k) same = same && bits32(on.snap[0].q[k]) == bits32(off.snap[0].q[k]);
            CHECK(same, "(f) 1 節点: ON と OFF の commit 後の Q がビット一致しない (ke_N が float で正確な場合)");
            CHECK(on.snap[0].q[2] == 3e-3f && on.snap[0].qa[2] == (double)3e-3f, "(f) 1 節点: commit 後の ρu_r が dq_r でない");
            CHECK(on.snap[0].qa[4] == (1.9921875 + (double)(-4e-3f)), "(f) 1 節点: Qacc_ρE が (2 − ke_N) + dq でない");
        }
        // 乱数 20 万
        {
            const int M = 200000;
            Lcg r(4242u);
            std::vector<Node> in(M);
            std::vector<flow_float> d(5 * M);
            for (int i = 0; i < M; ++i) {
                const flow_float rN = (flow_float)r.logU(-3.0, 3.0);
                const flow_float uy = (flow_float)(r.sgn() * r.logU(-6.0, 2.0));
                const flow_float e = (flow_float)(0.5 * (double)uy * uy / rN * (1.0 + r.u()) + r.logU(-3.0, 6.0));
                in[i] = makeNode(rN, (flow_float)(r.sgn() * r.logU(-3.0, 3.0)), uy, (flow_float)(r.sgn() * r.logU(-3.0, 3.0)), e, 1, 0);
                for (int k = 0; k < 5; ++k) d[5*i + k] = (flow_float)(r.sgn() * r.logU(-9.0, -1.0) * std::fabs((double)in[i].qN[k] + 1e-3));
            }
            const RunOut on = run(in, d, ON, 0, 1), off = run(in, d, OFF, 0, 1);
            long long nExact = 0, nMirror = 0, nOther = 0, nRoeDiff = 0, nWorse = 0;
            double maxResidUlp = 0.0;
            for (int i = 0; i < M; ++i) {
                const Node& x0 = in[i];
                const flow_float keN = (flow_float)0.5 * x0.qN[2] * x0.qN[2] / x0.qN[0];   // 製品と同じ式 (host の IEEE float)
                const double wantQa4 = ((double)x0.qN[4] - (double)keN) + (double)d[5*i + 4];
                const double wantQa2 = 0.0 + (double)d[5*i + 2];
                const Snap& s = on.snap[i];
                if (bits64(s.qa[4]) != bits64(wantQa4) || bits64(s.qa[2]) != bits64(wantQa2)) ++nExact;
                bool mir = true;
                for (int k = 0; k < 5; ++k) mir = mir && bits32(s.q[k]) == bits32((flow_float)s.qa[k]);
                if (!mir) ++nMirror;
                for (int k = 0; k < 4; ++k) if (bits32(s.q[k]) != bits32(off.snap[i].q[k])) { ++nOther; break; }
                if (bits32(s.q[4]) != bits32(off.snap[i].q[4])) {
                    ++nRoeDiff;
                    // ON は正本の最近接の float (OFF は ke_N を引いた時点で 1 回丸めるので、それより近くはならない)
                    if (std::fabs((double)s.q[4] - s.qa[4]) > std::fabs((double)off.snap[i].q[4] - s.qa[4])) ++nWorse;
                }
                maxResidUlp = std::fmax(maxResidUlp, std::fabs(s.qa[4] - (double)s.q[4]) / ulpf(s.q[4]));
            }
            printf("    乱数 %d 節点 (残余ゼロから 1 step): Qacc が (Q_N − ke_N) + dq と不一致 %lld、Q ≠ (float)Qacc %lld、"
                   "ρ・ρu_x・ρu_r・ρu_z の ON/OFF 不一致 %lld\n", M, nExact, nMirror, nOther);
            printf("    ρE が ON と OFF で違う節点 %lld (軸ピンの ke_N を FP64 で引くので 1 回分の丸めが残余に残る。想定どおり)、"
                   "うち ON の方が正本から遠い %lld、ON の残余の最大 %.3f ULP\n", nRoeDiff, nWorse, maxResidUlp);
            CHECK(nExact == 0, "(f) Qacc が (Q_N − ke_N) + dq と一致しない (%lld)", nExact);
            CHECK(nMirror == 0, "(f) commit 後の Q が (float)Qacc でない (%lld)", nMirror);
            CHECK(nOther == 0, "(f) ρ・ρu_x・ρu_r・ρu_z が ON と OFF でビット一致しない (%lld)", nOther);
            CHECK(nWorse == 0, "(f) ON の ρE が OFF より正本から遠い (%lld)", nWorse);
            CHECK(maxResidUlp <= 0.5, "(f) ON の残余が 1/2 ULP を超えた (%.3f)", maxResidUlp);
            CHECK(on.nAdopt == 0, "(f) reconcile が採用した (%d)", on.nAdopt);
        }
    }

    // ---- (g) 残余を持つ複数 step ----
    {
        // 値はすべて 2 進の有限桁にして、double の期待値を正確に書けるようにする。
        // 初期: ρ=1, ρu_r=0.125, ρE=1.5 (Q)、Qacc_ρE = 1.5 + 2^-25 (= 0.25 ULP の残余)。毎 step dq_r = 2^-10、dq_E = 2^-26 (1/8 ULP)。
        const int N = 100;
        Node x = makeNode(1.0f, 0.0f, 0.125f, 0.0f, 1.5f, 1, 0);
        x.qa[4] = 1.5 + std::ldexp(1.0, -25);
        const std::vector<flow_float> d{0.0f, 0.0f, (flow_float)std::ldexp(1.0, -10), 0.0f, (flow_float)std::ldexp(1.0, -26)};
        const RunOut on = run({x}, d, ON, 0, N), off = run({x}, d, OFF, 0, N), old = run({x}, d, OLD, 0, N);
        // 期待値: step 1 の ke_N = 2^-7、以降は ke_N = 0.5·(2^-10)^2 = 2^-21。
        const double wantQa4 = 1.5 + std::ldexp(1.0, -25) - std::ldexp(1.0, -7) + N * std::ldexp(1.0, -26)
                             - (N - 1) * std::ldexp(1.0, -21);
        const flow_float wantOff = (flow_float)(1.4921875 - (N - 1) * std::ldexp(1.0, -21));   // dq_E は毎回丸めで消える
        const double ulp = std::ldexp(1.0, -23);
        printf("\n(g) 残余を持つ %d step (残余 0.25 ULP、毎 step dq_r = 2^-10、dq_E = 1/8 ULP):\n", N);
        printf("    ON : Qacc_ρE = %.17g (期待 %.17g)、Q_ρE = %.9g、ρu_r = %.9g、採用 %d\n",
               on.nodes[0].qa[4], wantQa4, (double)on.nodes[0].q[4], (double)on.nodes[0].q[2], on.nAdopt);
        printf("    OFF: Q_ρE = %.9g (期待 %.9g)、ρu_r = %.9g。ON − OFF = %.2f ULP (dq_E の累積 12.5 ULP + 初期残余 0.25 ULP を反映)\n",
               (double)off.nodes[0].q[4], (double)wantOff, (double)off.nodes[0].q[2],
               ((double)on.nodes[0].q[4] - (double)off.nodes[0].q[4]) / ulp);
        printf("    OLD: Q_ρE = %.9g、ρu_r = %.9g (ピンが巻き戻り ρu_r が step ごとに積もる)\n",
               (double)old.nodes[0].q[4], (double)old.nodes[0].q[2]);
        CHECK(bits64(on.nodes[0].qa[4]) == bits64(wantQa4), "(g) ON: Qacc_ρE が期待値と一致しない");
        CHECK(on.nodes[0].q[4] == (flow_float)wantQa4, "(g) ON: Q_ρE が (float)Qacc でない");
        CHECK(on.nodes[0].q[2] == (flow_float)std::ldexp(1.0, -10) && on.nodes[0].qa[2] == std::ldexp(1.0, -10),
              "(g) ON: ρu_r が 1 step 分の dq_r (2^-10) でない");
        CHECK(on.nAdopt == 0, "(g) ON: reconcile が採用した (%d)", on.nAdopt);
        CHECK(off.nodes[0].q[4] == wantOff, "(g) OFF: Q_ρE が期待値でない");
        CHECK(off.nodes[0].q[2] == on.nodes[0].q[2], "(g) ρu_r が ON と OFF で違う");
        CHECK(std::fabs(((double)on.nodes[0].q[4] - (double)off.nodes[0].q[4]) / ulp - 13.0) < 1e-12,
              "(g) ON − OFF が 13 ULP (= round(12.75)) でない");
        CHECK(old.nodes[0].q[2] == (flow_float)(0.125 + N * std::ldexp(1.0, -10)),
              "(g) OLD: ρu_r が 0.125 + N·2^-10 に積もらない (旧い作りの欠陥が再現しない)");
    }

    printf("\n%s (失敗 %d)\n", g_fail ? "VERDICT: FAIL" : "VERDICT: PASS", g_fail);
    return g_fail ? 1 : 0;
}
