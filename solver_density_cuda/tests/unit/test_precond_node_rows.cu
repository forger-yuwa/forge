// test_precond_node_rows.cu — 前処理の block カーネル (lowMachPrecond>=2) の node の行の処理の単体試験
//   (plans/active/architecture-float-state-double-geometry.md §5.1 #26 の R2b: β < 1 で拘束の行が行列のレベルで正しいこと)
//
// **製品のカーネル implicit_defect_correction_block_precond_d (cuda_forge/timeIntegration_d.cu) をそのまま起動する**。
// カーネルは別の翻訳単位 (ビルド済みの timeIntegration_d.cu.o) にあり、ここでは宣言だけを書いてリンクする
// (rdc なしでも __global__ の起動は翻訳単位をまたげる)。timeIntegration_d.cu.o の外部参照は weakIsoWall::active /
// diagBuf だけで、前処理のカーネルは呼ばないので下で止まる実体 (abort) を置く。
//
// build (solver_density_cuda で。<obj> は試験するビルドの cuda_forge/CMakeFiles/cuda_forge.dir/timeIntegration_d.cu.o):
//   nvcc -std=c++17 --expt-relaxed-constexpr -arch=sm_86 -I. -o test_precond_node_rows tests/unit/test_precond_node_rows.cu <obj>
//   FP64 のビルドの object を試すときは、typedef を double にした flowFormat.hpp の写しを置いた <fp64> を先に探させる
//   (-I<fp64> -I.。この試験は flowFormat.hpp のほかに forge のヘッダを読まない)。flow_float の型がカーネルの宣言の型と
//   揃わないとリンクで落ちる。
//   移植前のカーネル (引数 4 つ少ない) と比べるときは -DBASE_KERNEL で移植前の object にリンクする (項目 (5) の dump だけ使う)。
// run: ./test_precond_node_rows [dump.bin]   (dump.bin を渡すと項目 (5) の出力を書く。基準との比較は cmp で行う)
//   環境変数 PNR_VERBOSE を立てると項目 (2) の成分ごとの差を出す。
//
// 試験の場:
//   内部節点 4 個 (0〜3) とゴースト 3 個 (4〜6)。面 8 枚 (内部 5・境界半割面 3)。節点 3 だけ境界半割面を持たない。
//   低マッハ (|u| ≈ 3 m/s, c ≈ 347 m/s, M ≈ 0.01)、precondEps 0.15 → β = 0.0225 (Γ_c のランク 1 項が (1−β)/β ≈ 43 倍)。
//   軸対称 (isAxisymmetric 1) なので軸の行 (roUy) にはソース Jacobian が乗っており、拘束はそれを上書きする必要がある。
//   印: 節点 0 = 軸 (roUy 行)、節点 1 = 壁 (運動量 3 行) + 等温壁 (roe 行)、節点 2・3 = 拘束なし。
// 項目:
//   (1) 拘束の行の更新は 0: 乱数の残差と乱数の近傍 dq_old (sweep の 2 回目以降) で、節点 0 の dq[2] と節点 1 の dq[1..4] が
//       全部厳密に 0。前提として D0 の拘束の列で部分ピボットの行交換が起きないこと (|D0[i][r]| < 0.5, i≠r) を確かめる
//       (交換が起きると標準の block カーネルと同じく丸め誤差ぶんだけ 0 から外れる)。ランク 1 項 γ g rᵀ を拘束の行から
//       外していなければ、その行の解は −γ g[r] (rᵀx) で 0 にならない (この項目が g[row] = 0 の検査を兼ねる)。
//   (2) 拘束の無い成分は、拘束の無い同じ行列 D (= D0 + γ g rᵀ、isNode 1・印なしで単位ベクトルへの応答から D⁻¹ を取る) の
//       部分行列の解 x_U = (D_UU)⁻¹ b_U に一致する。(D_UU)⁻¹ = M_UU − M_UC M_CC⁻¹ M_CU (M = D⁻¹) で、D を逆算しない。
//       比較の感度: 同じ式に β = 1 の M (ランク 1 項なし) を入れると、許容差を大きく超えて外れることも確かめる。
//       許容差は float のビルドで 2e-2 (丸めの拡大ぶん。粗い確認)、FP64 のビルドで 1e-8 (厳密な確認)。両方のビルドで回す。
//   (3) 拘束の無い節点 (2・3) の出力は、印の有無で生のビットが一致する。
//   (4) isNode 1 は境界半割面の粘性対角だけを外す: 境界半割面を持たない節点 3 は isNode 0/1 で生のビットが一致し、
//       持つ節点は D(isNode 0) − D(isNode 1) が 2ν·δ/dcc·I (host で同じ式) に一致する (β = 1 で比較)。
//       さらに node→ghost が退化した (dcc を 1e-7 倍にした) 境界の節点は isNode 0 で応答がほぼ 0 に凍結し、isNode 1 では凍結しない。
//   (5) dump: isNode 0・印なし (cell と同じ経路) の出力を書く。移植前の object (-DBASE_KERNEL) の出力と cmp で比べ、
//       前処理のカーネルの cell の経路が移植の前後でビット一致することを確かめる。
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <cuda_runtime.h>
#include "flowFormat.hpp"

class solverConfig;
class mesh;
namespace weakIsoWall {
bool active(const solverConfig&, const mesh&) { fprintf(stderr, "weakIsoWall::active は試験では呼ばれない\n"); abort(); }
flow_float* diagBuf(const mesh&) { fprintf(stderr, "weakIsoWall::diagBuf は試験では呼ばれない\n"); abort(); }
}

// 製品のカーネルの宣言 (cuda_forge/timeIntegration_d.cu と同じ並び・同じ型)。
__global__ void __launch_bounds__(128) implicit_defect_correction_block_precond_d
(
 int loop, flow_float dt, flow_float* dt_local, flow_float implicit_relax, flow_float* gamma_arr, flow_float precondEps,
 geom_int nCells_all, geom_int nCells, geom_float* vol,
 geom_int* plane_cells, geom_int* cell_planes_index, geom_int* cell_planes,
 geom_float* ccx, geom_float* ccy, geom_float* ccz,
 const flow_float* ge_x, const flow_float* ge_y, const flow_float* ge_z,
 geom_float* sx, geom_float* sy, geom_float* sz, geom_float* ss,
 flow_float* ro, flow_float* roUx, flow_float* roUy, flow_float* roUz, flow_float* roe,
 flow_float laminar_visc, flow_float* vis_turb, flow_float* sonic,
 flow_float* Ux, flow_float* Uy, flow_float* Uz, flow_float* Ht,
 flow_float* res_ro, flow_float* res_roUx, flow_float* res_roUy, flow_float* res_roUz, flow_float* res_roe,
 flow_float* dq_old_0, flow_float* dq_old_1, flow_float* dq_old_2, flow_float* dq_old_3, flow_float* dq_old_4,
 flow_float* dq_new_0, flow_float* dq_new_1, flow_float* dq_new_2, flow_float* dq_new_3, flow_float* dq_new_4,
 int isAxisymmetric, flow_float* A_planar, flow_float axisRFloor, flow_float unsteady_diag
#ifndef BASE_KERNEL
 , const geom_int* axis_ur_flag, const geom_int* wall_flag, const geom_int* iso_wall_flag, int isNode
#endif
);

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { ++g_fail; printf("  FAIL: " __VA_ARGS__); printf("\n"); } } while (0)
#define CU(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { fprintf(stderr, "CUDA %s @%d: %s\n", #x, __LINE__, cudaGetErrorString(e_)); exit(2); } } while (0)

constexpr int NC = 4, NALL = 7, NP = 8;
// 4・5・6 はゴースト (境界半割面。節点 0・1・2 にそれぞれ 1 枚。節点 3 は内部の面だけ)
const int kPlanes[NP][2] = {{0,1},{1,2},{2,3},{0,3},{0,4},{1,3},{1,5},{2,6}};

template<typename T> static T* dev(const std::vector<T>& h) {
    T* d = nullptr; CU(cudaMalloc((void**)&d, h.size() * sizeof(T)));
    CU(cudaMemcpy(d, h.data(), h.size() * sizeof(T), cudaMemcpyHostToDevice)); return d;
}

struct Field {
    // 幾何
    std::vector<geom_float> vol, ccx, ccy, ccz, sx, sy, sz, ss;
    std::vector<flow_float> gex, gey, gez, Apl;
    std::vector<geom_int> pc, cpi, cp;
    // 状態
    std::vector<flow_float> dt, gam, ro, roUx, roUy, roUz, roe, vt, c, Ux, Uy, Uz, Ht;
    flow_float mu = 1.8e-5;
};

// ghostScale: ゴーストの位置を境界面の外へ置く距離の倍率 (1 で内部の間隔の半分、1e-7 で退化)
static Field makeField(double ghostScale)
{
    Field f;
    const double h = 2.0e-4;   // 間隔 [m]。軸対称の近軸の小さい面 (S ~ 1e-8 m²) を想定し、D0 の成分を 1 未満に保つ
    const double g = 0.5 * h * ghostScale;   // 節点→ゴーストの距離
    const double X[NALL][3] = {{0,0,0},{h,0.1*h,0},{2*h,-0.05*h,0},{1.1*h,h,0},
                               {-g, -0.04*g, 0}, {h + 0.05*g, 0.1*h - g, 0}, {2*h + g, -0.05*h + 0.06*g, 0}};
    const double depth = 5.0e-5;   // 奥行き (軸対称の r 重み相当)
    f.ccx.resize(NALL); f.ccy.resize(NALL); f.ccz.resize(NALL); f.vol.resize(NALL);
    for (int i = 0; i < NALL; ++i) { f.ccx[i] = (geom_float)X[i][0]; f.ccy[i] = (geom_float)(X[i][1] + 1.0e-3); f.ccz[i] = (geom_float)X[i][2]; }
    const double V[NC] = {0.6, 1.0, 0.8, 0.9};
    for (int i = 0; i < NC; ++i) f.vol[i] = (geom_float)(V[i] * h * h * depth);
    for (int i = NC; i < NALL; ++i) f.vol[i] = (geom_float)(0.5 * h * h * depth);
    f.sx.resize(NP); f.sy.resize(NP); f.sz.resize(NP); f.ss.resize(NP); f.gex.resize(NP); f.gey.resize(NP); f.gez.resize(NP);
    f.pc.resize(2 * NP);
    for (int p = 0; p < NP; ++p) {
        const int a = kPlanes[p][0], b = kPlanes[p][1];
        f.pc[2*p] = a; f.pc[2*p+1] = b;
        // 面ベクトル: 節点間の向きから少し傾けた法線 × 面積 (非直交の成分を入れる)
        const double ex = X[b][0] - X[a][0], ey = X[b][1] - X[a][1];
        double nx = ex + 0.15 * ey, ny = ey - 0.1 * ex;
        const double nn = std::sqrt(nx*nx + ny*ny); nx /= nn; ny /= nn;
        const double area = (p < 4 || p == 5) ? 0.9 * h * depth : 0.7 * h * depth;
        f.sx[p] = (geom_float)(nx * area); f.sy[p] = (geom_float)(ny * area); f.sz[p] = 0;
        f.ss[p] = (geom_float)area;
        f.gex[p] = (flow_float)ex; f.gey[p] = (flow_float)ey; f.gez[p] = 0;
    }
    // CSR (内部節点ごとの面)
    f.cpi.assign(NALL + 1, 0);
    std::vector<std::vector<int>> lst(NALL);
    for (int p = 0; p < NP; ++p) { lst[kPlanes[p][0]].push_back(p); lst[kPlanes[p][1]].push_back(p); }
    for (int i = 0; i < NALL; ++i) { f.cpi[i+1] = f.cpi[i] + (int)lst[i].size(); for (int p : lst[i]) f.cp.push_back(p); }
    f.Apl.resize(NALL); for (int i = 0; i < NALL; ++i) f.Apl[i] = (flow_float)(f.vol[i] / (f.ccy[i]));   // A_planar ~ V/r

    // 低マッハの状態 (CPG, γ = 1.4, R = 287)
    const double gm = 1.4, R = 287.0, T0 = 300.0, P0 = 1.0e5;
    const double U[NALL][3] = {{3.0,0.0,0.0},{2.5,1.2,0.0},{3.2,-0.8,0.1},{2.8,0.5,-0.2},{3.0,0.0,0.0},{3.0,0.0,0.0},{3.0,0.0,0.0}};
    f.dt.resize(NALL); f.gam.resize(NALL); f.ro.resize(NALL); f.roUx.resize(NALL); f.roUy.resize(NALL); f.roUz.resize(NALL);
    f.roe.resize(NALL); f.vt.resize(NALL); f.c.resize(NALL); f.Ux.resize(NALL); f.Uy.resize(NALL); f.Uz.resize(NALL); f.Ht.resize(NALL);
    for (int i = 0; i < NALL; ++i) {
        const double T = T0 * (1.0 + 0.01 * i), rho = P0 / (R * T), cs = std::sqrt(gm * R * T);
        const double ek = 0.5 * (U[i][0]*U[i][0] + U[i][1]*U[i][1] + U[i][2]*U[i][2]);
        f.gam[i] = (flow_float)gm; f.ro[i] = (flow_float)rho; f.c[i] = (flow_float)cs;
        f.Ux[i] = (flow_float)U[i][0]; f.Uy[i] = (flow_float)U[i][1]; f.Uz[i] = (flow_float)U[i][2];
        f.roUx[i] = (flow_float)(rho*U[i][0]); f.roUy[i] = (flow_float)(rho*U[i][1]); f.roUz[i] = (flow_float)(rho*U[i][2]);
        f.roe[i] = (flow_float)(P0/(gm-1.0) + rho*ek);
        f.Ht[i] = (flow_float)(cs*cs/(gm-1.0) + ek);
        f.vt[i] = (flow_float)(5.0e-2 * (1 + i));   // 乱流粘性 (粘性対角を V/Δτ と同程度にし、(4b) の差を D の丸めより十分大きくする)
        f.dt[i] = (flow_float)(5.0e-8 * (1.0 + 0.2 * i));   // 局所の擬似時間刻み (前処理で拡大済みの値として渡す)
    }
    return f;
}

struct Run {
    // 入力
    int loop = 1; flow_float eps = 0.15; int isNode = 1; bool flags = false; bool axisym = true;
    std::vector<flow_float> res[5];   // NALL
    std::vector<flow_float> dqo[5];   // NALL
    // 出力
    std::vector<flow_float> dq[5];
};

static void launch(const Field& f, Run& r)
{
    std::vector<geom_int> axis(NALL, 0), wall(NALL, 0), iso(NALL, 0);
    axis[0] = 1; wall[1] = 1; iso[1] = 1;
    geom_float *vol = dev(f.vol), *ccx = dev(f.ccx), *ccy = dev(f.ccy), *ccz = dev(f.ccz);
    geom_float *sx = dev(f.sx), *sy = dev(f.sy), *sz = dev(f.sz), *ss = dev(f.ss);
    flow_float *gex = dev(f.gex), *gey = dev(f.gey), *gez = dev(f.gez), *Apl = dev(f.Apl);
    geom_int *pc = dev(f.pc), *cpi = dev(f.cpi), *cp = dev(f.cp);
    flow_float *dt = dev(f.dt), *gam = dev(f.gam), *ro = dev(f.ro), *rux = dev(f.roUx), *ruy = dev(f.roUy), *ruz = dev(f.roUz), *roe = dev(f.roe);
    flow_float *vt = dev(f.vt), *c = dev(f.c), *ux = dev(f.Ux), *uy = dev(f.Uy), *uz = dev(f.Uz), *Ht = dev(f.Ht);
    flow_float *res[5], *dqo[5], *dqn[5];
    for (int k = 0; k < 5; ++k) { res[k] = dev(r.res[k]); dqo[k] = dev(r.dqo[k]); dqn[k] = dev(std::vector<flow_float>(NALL, (flow_float)-777.0)); }
    geom_int *axd = dev(axis), *wad = dev(wall), *isd = dev(iso);
    implicit_defect_correction_block_precond_d<<<1, 128>>>(
        r.loop, (flow_float)1.0e-6, dt, (flow_float)1.0, gam, r.eps,
        NALL, NC, vol, pc, cpi, cp, ccx, ccy, ccz, gex, gey, gez, sx, sy, sz, ss,
        ro, rux, ruy, ruz, roe, f.mu, vt, c, ux, uy, uz, Ht,
        res[0], res[1], res[2], res[3], res[4], dqo[0], dqo[1], dqo[2], dqo[3], dqo[4],
        dqn[0], dqn[1], dqn[2], dqn[3], dqn[4],
        r.axisym ? 1 : 0, Apl, (flow_float)0.0, (flow_float)0.0
#ifndef BASE_KERNEL
        , r.flags ? axd : nullptr, r.flags ? wad : nullptr, r.flags ? isd : nullptr, r.isNode
#endif
    );
    CU(cudaGetLastError()); CU(cudaDeviceSynchronize());
    for (int k = 0; k < 5; ++k) {
        r.dq[k].resize(NALL);
        CU(cudaMemcpy(r.dq[k].data(), dqn[k], NALL * sizeof(flow_float), cudaMemcpyDeviceToHost));
        cudaFree(res[k]); cudaFree(dqo[k]); cudaFree(dqn[k]);
    }
    for (void* p : {(void*)vol,(void*)ccx,(void*)ccy,(void*)ccz,(void*)sx,(void*)sy,(void*)sz,(void*)ss,(void*)gex,(void*)gey,(void*)gez,(void*)Apl,
                    (void*)pc,(void*)cpi,(void*)cp,(void*)dt,(void*)gam,(void*)ro,(void*)rux,(void*)ruy,(void*)ruz,(void*)roe,(void*)vt,(void*)c,
                    (void*)ux,(void*)uy,(void*)uz,(void*)Ht,(void*)axd,(void*)wad,(void*)isd}) cudaFree(p);
}

static Run blank() { Run r; for (int k = 0; k < 5; ++k) { r.res[k].assign(NALL, 0); r.dqo[k].assign(NALL, 0); } return r; }

typedef double Mat5[5][5];
// 全節点の M = D⁻¹ を単位ベクトルへの応答から取る (dq_old = 0、relax = 1)。M[node][i][j] = 応答 i / 入力 j。
static void inverseOperator(const Field& f, flow_float eps, int isNode, Mat5* M)
{
    for (int j = 0; j < 5; ++j) {
        Run r = blank(); r.eps = eps; r.isNode = isNode; r.flags = false;
        for (int n = 0; n < NC; ++n) r.res[j][n] = (flow_float)1.0;
        launch(f, r);
        for (int n = 0; n < NC; ++n) for (int i = 0; i < 5; ++i) M[n][i][j] = (double)r.dq[i][n];
    }
}

// 小さい密行列の逆 (部分ピボット Gauss-Jordan、double)
static bool invert(int n, const double* A, double* Ai)
{
    std::vector<double> a(A, A + n*n), b(n*n, 0.0);
    for (int i = 0; i < n; ++i) b[i*n+i] = 1.0;
    for (int col = 0; col < n; ++col) {
        int p = col; for (int r = col+1; r < n; ++r) if (std::fabs(a[r*n+col]) > std::fabs(a[p*n+col])) p = r;
        if (std::fabs(a[p*n+col]) < 1e-300) return false;
        for (int k = 0; k < n; ++k) { std::swap(a[col*n+k], a[p*n+k]); std::swap(b[col*n+k], b[p*n+k]); }
        const double iv = 1.0 / a[col*n+col];
        for (int k = 0; k < n; ++k) { a[col*n+k] *= iv; b[col*n+k] *= iv; }
        for (int r = 0; r < n; ++r) if (r != col) { const double fct = a[r*n+col]; if (fct != 0.0) for (int k = 0; k < n; ++k) { a[r*n+k] -= fct*a[col*n+k]; b[r*n+k] -= fct*b[col*n+k]; } }
    }
    std::copy(b.begin(), b.end(), Ai); return true;
}

// x_U = (M_UU − M_UC M_CC⁻¹ M_CU) b_U、x_C = 0
static void expectedConstrained(const Mat5& M, const bool C[5], const double b[5], double x[5])
{
    int iu[5], ic[5], nu = 0, nc = 0;
    for (int i = 0; i < 5; ++i) { if (C[i]) ic[nc++] = i; else iu[nu++] = i; }
    double Mcc[25], Mcci[25];
    for (int a = 0; a < nc; ++a) for (int bb = 0; bb < nc; ++bb) Mcc[a*nc+bb] = M[ic[a]][ic[bb]];
    if (nc > 0 && !invert(nc, Mcc, Mcci)) { printf("  M_CC が特異\n"); exit(2); }
    for (int i = 0; i < 5; ++i) x[i] = 0.0;
    for (int a = 0; a < nu; ++a) {
        double s = 0.0;
        for (int bb = 0; bb < nu; ++bb) {
            double S = M[iu[a]][iu[bb]];
            for (int p = 0; p < nc; ++p) for (int q = 0; q < nc; ++q) S -= M[iu[a]][ic[p]] * Mcci[p*nc+q] * M[ic[q]][iu[bb]];
            s += S * b[iu[bb]];
        }
        x[iu[a]] = s;
    }
}

static bool bitsEqual(flow_float a, flow_float b) { return std::memcmp(&a, &b, sizeof(flow_float)) == 0; }

int main(int argc, char** argv)
{
    const Field f = makeField(1.0);
    std::mt19937_64 rng(20261012);
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    // 残差の成分ごとの大きさ (ρ, ρu, ρv, ρw, ρE の次元の違いを模す)
    const double rs[5] = {1.0e-6, 3.0e-6, 3.0e-6, 3.0e-6, 0.3};
    auto randRun = [&](bool withNbr) {
        Run r = blank();
        for (int k = 0; k < 5; ++k) for (int n = 0; n < NALL; ++n) {
            r.res[k][n] = (flow_float)(rs[k] * U(rng));
            if (withNbr) r.dqo[k][n] = (flow_float)(rs[k] * 1.0e2 * U(rng));
        }
        return r;
    };

#ifdef BASE_KERNEL
    // 移植前の object: 項目 (5) の dump だけ
    if (argc < 2) { printf("BASE_KERNEL では dump の出力先を渡す\n"); return 2; }
#endif
    // ---------------------------------------------------------------- (5) cell と同じ経路の dump
    {
        FILE* fp = (argc >= 2) ? fopen(argv[1], "wb") : nullptr;
        for (int t = 0; t < 8; ++t) {
            Run r = randRun(true); r.loop = 1 + (t % 3); r.eps = (t < 6) ? (flow_float)0.15 : (flow_float)1.0;
            r.isNode = 0; r.flags = false; r.axisym = (t % 2 == 0);
            launch(f, r);
            if (fp) for (int k = 0; k < 5; ++k) fwrite(r.dq[k].data(), sizeof(flow_float), NC, fp);
        }
        if (fp) { fclose(fp); printf("(5) isNode 0・印なしの出力 8 組を %s に書いた (移植前の object の出力と cmp で比べる)\n", argv[1]); }
    }
#ifdef BASE_KERNEL
    return 0;
#else
    // ---------------------------------------------------------------- 準備: M = D⁻¹ (isNode 1、印なし)
    static Mat5 Mb[NC], M1[NC];   // β < 1 (eps 0.15) と β = 1 (eps 1.0)
    inverseOperator(f, (flow_float)0.15, 1, Mb);
    inverseOperator(f, (flow_float)1.0, 1, M1);
    // β < 1 の効き: M(β) と M(1) の差が M(1) に比べて十分大きい
    {
        double dmax = 0, ref = 0;
        for (int n = 0; n < NC; ++n) for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) {
            const double sc = std::sqrt(std::fabs(M1[n][i][i] * M1[n][j][j])) + 1e-300;
            dmax = std::max(dmax, std::fabs(Mb[n][i][j] - M1[n][i][j]) / sc); ref = std::max(ref, std::fabs(M1[n][i][j]) / sc);
        }
        printf("準備: β < 1 と β = 1 の D⁻¹ の差 (対角で正規化した最大) %.3e (D⁻¹(β=1) 自体 %.3e)\n", dmax, ref);
        CHECK(dmax > 1.0e-2, "β < 1 のランク 1 項が D⁻¹ にほとんど効いていない (%.3e)。試験の場が弱い", dmax);
    }
    // 前提 (項目 1): D0 = M(β=1)⁻¹ の拘束の列 r で、拘束の行以外の成分が 1 未満 (部分ピボットの行交換が起きない)
    {
        bool ok = true;
        const int cons[2][5] = {{0,0,1,0,0},{0,1,1,1,1}};
        for (int n = 0; n < 2; ++n) {
            double A[25], Ai[25];
            for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) A[i*5+j] = M1[n][i][j];
            invert(5, A, Ai);   // Ai = D0 (β = 1 では γ = 0)
            for (int r = 0; r < 5; ++r) if (cons[n][r]) for (int i = 0; i < 5; ++i) if (i != r && std::fabs(Ai[i*5+r]) >= 0.5) ok = false;
            double mx = 0; for (int i = 0; i < 25; ++i) mx = std::max(mx, std::fabs(Ai[i]));
            printf("準備: 節点 %d の D0 の成分の最大 %.3e\n", n, mx);
        }
        CHECK(ok, "D0 の拘束の列に 0.5 以上の成分がある (部分ピボットの行交換で 0 が丸めに崩れうる。場の大きさを見直す)");
    }

    // ---------------------------------------------------------------- (1)・(2)・(3)
    {
        const bool C[NC][5] = {{false,false,true,false,false},{false,true,true,true,true},{false},{false}};
        int nZero = 0, nCons = 0, nBitSame = 0, nFree = 0;
        // (2) の許容差: 期待値もカーネルの出力も flow_float の解から作るので、Sherman-Morrison の打ち消し ((1−β)/β ≈ 43 倍) と
        // D の条件数 (~1e3) のぶん丸めが拡大する。小さい成分 (ρ・ρu、エネルギー成分の 1e-4 倍) の相対差で ε·43·1e3 ≈ 3e-3 程度。
        // float は 2e-2 (粗い整合の確認、感度は下の β = 1 の対照で示す)、double は 1e-8 (行列のレベルの厳密な確認)。
        const double tol2 = (sizeof(flow_float) == 4) ? 2.0e-2 : 1.0e-8;
        double errMax = 0.0, errNoRank1 = 0.0;
        for (int t = 0; t < 40; ++t) {
            const bool withNbr = (t >= 20);
            Run r = randRun(withNbr); r.loop = withNbr ? 2 : 1; r.eps = (flow_float)0.15; r.isNode = 1; r.flags = true;
            Run q = r; q.flags = false;
            launch(f, r); launch(f, q);
            for (int n = 0; n < NC; ++n) {
                // (1) 拘束の行は厳密に 0
                for (int i = 0; i < 5; ++i) if (C[n][i]) { ++nCons; if (r.dq[i][n] == (flow_float)0.0) ++nZero; else printf("  節点 %d 行 %d: %.6e\n", n, i, (double)r.dq[i][n]); }
                // (3) 拘束の無い節点は印の有無でビット一致
                if (n >= 2) for (int i = 0; i < 5; ++i) { ++nFree; if (bitsEqual(r.dq[i][n], q.dq[i][n])) ++nBitSame; }
                // (2) 拘束の無い成分 (近傍なしの試行だけ: b がそのまま残差)
                if (!withNbr && n < 2) {
                    double b[5], x[5], x1[5];
                    for (int i = 0; i < 5; ++i) b[i] = (double)r.res[i][n];
                    expectedConstrained(Mb[n], C[n], b, x);
                    expectedConstrained(M1[n], C[n], b, x1);   // 感度: ランク 1 項を無視した期待値
                    for (int i = 0; i < 5; ++i) if (!C[n][i]) {
                        const double sc = std::fabs(x[i]) + 1.0e-3 * rs[i] * std::fabs(Mb[n][i][i]);   // 偶然 0 に近い成分の下限
                        errMax = std::max(errMax, std::fabs((double)r.dq[i][n] - x[i]) / sc);
                        if (getenv("PNR_VERBOSE")) printf("    t%d n%d i%d x_on %.9e x_exp %.9e rel %.3e\n", t, n, i, (double)r.dq[i][n], x[i], std::fabs((double)r.dq[i][n] - x[i]) / sc);
                        errNoRank1 = std::max(errNoRank1, std::fabs(x1[i] - x[i]) / sc);
                    }
                }
            }
        }
        printf("(1) 拘束の行の dq が厳密に 0: %d / %d (近傍 dq_old あり・なしの各 20 試行)\n", nZero, nCons);
        CHECK(nZero == nCons, "拘束の行の dq が 0 でない (%d / %d)", nCons - nZero, nCons);
        printf("(2) 拘束の無い成分と (D_UU)⁻¹ b_U の相対差の最大 %.3e (許容 %.0e)。β = 1 の M で作った期待値との相対差 %.3e\n", errMax, tol2, errNoRank1);
        CHECK(errMax < tol2, "拘束の無い成分が期待値から外れる (%.3e)", errMax);
        CHECK(errNoRank1 > 1.0e2 * tol2, "比較の感度が足りない: ランク 1 項を無視しても許容差の 100 倍以内 (%.3e)", errNoRank1);
        printf("(3) 拘束の無い節点 (2・3) の出力が印の有無で生のビット一致: %d / %d\n", nBitSame, nFree);
        CHECK(nBitSame == nFree, "拘束の無い節点の出力が印で変わった");
    }

    // ---------------------------------------------------------------- (4) isNode: 境界半割面の粘性対角
    {
        // 節点 3 (境界半割面なし) は isNode 0/1 でビット一致
        int same = 0, tot = 0;
        for (int t = 0; t < 10; ++t) {
            Run a = randRun(t >= 5); a.loop = (t >= 5) ? 2 : 1; a.eps = (flow_float)0.15; a.isNode = 1;
            Run b2 = a; b2.isNode = 0;
            launch(f, a); launch(f, b2);
            for (int i = 0; i < 5; ++i) { ++tot; if (bitsEqual(a.dq[i][3], b2.dq[i][3])) ++same; }
        }
        printf("(4a) 境界半割面を持たない節点 3 の出力が isNode 0/1 で生のビット一致: %d / %d\n", same, tot);
        CHECK(same == tot, "境界半割面を持たない節点の出力が isNode で変わった");

        // 境界半割面を持つ節点 0・1・2: D(isNode 0) − D(isNode 1) = Σ_ghost 2ν·δ/dcc·I (β = 1 で比較)
        static Mat5 M0n[NC], M1n[NC];
        inverseOperator(f, (flow_float)1.0, 0, M0n);
        inverseOperator(f, (flow_float)1.0, 1, M1n);
        double worst = 0.0, offd = 0.0;
        for (int n = 0; n < 3; ++n) {
            double A0[25], A1[25], D0n[25], D1n[25];
            for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) { A0[i*5+j] = M0n[n][i][j]; A1[i*5+j] = M1n[n][i][j]; }
            invert(5, A0, D0n); invert(5, A1, D1n);
            // host の期待値: 製品と同じ式 (dcc = |e|, delta = dcc·S²/|e·S|, 2ν·delta/dcc, ν = (μ + μ_t)/ρ)
            double v = 0.0;
            for (int p = 0; p < NP; ++p) {
                const int a = kPlanes[p][0], b = kPlanes[p][1];
                if (!((a == n && b >= NC) || (b == n && a >= NC))) continue;
                const double ex = f.gex[p], ey = f.gey[p], ez = f.gez[p];
                const double dcc = std::sqrt(ex*ex + ey*ey + ez*ez);
                const double ds = std::fabs(ex*f.sx[p] + ey*f.sy[p] + ez*f.sz[p]);
                const double S = f.ss[p];
                const double delta = dcc * S * S / ds;
                const double nu = ((double)f.mu + (double)f.vt[n]) / (double)f.ro[n];
                v += 2.0 * nu * delta / dcc;
            }
            double scale = 0.0; for (int i = 0; i < 25; ++i) scale = std::max(scale, std::fabs(D1n[i]));
            for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) {
                const double d = D0n[i*5+j] - D1n[i*5+j];
                if (i == j) worst = std::max(worst, std::fabs(d - v) / v);
                else offd = std::max(offd, std::fabs(d) / scale);
            }
            printf("(4b) 節点 %d: 境界半割面の粘性対角 2ν·δ/dcc = %.4e (D の成分の最大 %.3e)\n", n, v, scale);
        }
        printf("(4b) D(isNode 0) − D(isNode 1) の対角と 2ν·δ/dcc の相対差の最大 %.3e、非対角 (D の最大で正規化) %.3e (許容 1e-3)\n", worst, offd);
        CHECK(worst < 1.0e-3 && offd < 1.0e-3, "isNode が外す量が境界半割面の粘性対角と違う");

        // 退化した node→ghost: isNode 0 は境界の節点が凍結し、isNode 1 は凍結しない
        const Field fd = makeField(1.0e-7);
        Run a = randRun(false); a.eps = (flow_float)0.15; a.isNode = 0;
        Run b2 = a; b2.isNode = 1;
        launch(fd, a); launch(fd, b2);
        double ratioMax = 0.0, ratioMin = 1e300;
        for (int n = 0; n < 3; ++n) for (int i = 0; i < 5; ++i) {
            const double rr = std::fabs((double)a.dq[i][n]) / (std::fabs((double)b2.dq[i][n]) + 1e-300);
            ratioMax = std::max(ratioMax, rr); ratioMin = std::min(ratioMin, rr);
        }
        printf("(4c) node→ghost を 1e-7 倍に退化: 境界の節点の |dq(isNode 0)| / |dq(isNode 1)| = %.3e〜%.3e (旧い前処理のカーネルの凍結の再現)\n", ratioMin, ratioMax);
        CHECK(ratioMax < 1.0e-3, "退化した境界半割面で isNode 0 が凍結しない (試験の場が凍結を再現していない)");
    }

    printf(g_fail == 0 ? "PASS\n" : "FAIL (%d)\n", g_fail);
    return g_fail == 0 ? 0 : 1;
#endif
}
