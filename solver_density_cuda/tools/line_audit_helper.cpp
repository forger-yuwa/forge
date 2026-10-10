// 製品の経路の照合 (plan time_integration-line-viscous-jacobian-faceh §6.7) の補助: 監査用のビルド (-DFORGE_LINE_AUDIT) が書いた
// audit_face.f64・audit_node.f64 を読み、記録された ST (float) の入力から、製品と同じ共通関数・同じ加算の順序を host の float で再現する。
//   面ごと <dir>/audit_host_face.f64: [節点][AUDIT_NF][HF]
//     0-24   この面の後の D (記録された「この面の前の D」から再現)
//     25-49  この面の後の D (時間項の後の D から面の順に連ねて再現)
//     50-74  丸めの尺度 (前の D・後の D・対流の増分・粘性の増分の絶対値の最大、要素ごと)
//     75-99  対流の K (列の抽出)、100-124 薄層の D (零から)、125-149 薄層の K
//   節点ごと <dir>/audit_host_node.f64: [節点][HN]
//     0-24   D (時間項の後、生の体積と dt_local から)、25-49 D (軸対称の後、記録された面のループの後の D から)、50-74 軸対称の丸めの尺度
// 製品の加算の順序 (timeIntegration_d.cu の implicit_defect_correction_block_d、factor の sweep は sdq = 0):
//   面ごとに 対流 (accumulate_split_jacobian_cf) → 粘性 (has_nbr のとき): 枝 1 = [値 3 ならスカラーを単位行列で] + 薄層 / 枝 3 = スカラーを
//   行 0..nScalarRows−1 の対角 + 行 4 の温度の項 / 枝 4 = スカラーを単位行列で。軸対称 (r 重み) は面のループの後。
// 使い方: line_audit_helper <dir> <nNodes> <thermallyPerfect 0|1>
// ビルド: g++ -O2 -std=c++17 -Wno-unknown-pragmas -I solver_density_cuda solver_density_cuda/tools/line_audit_helper.cpp -o /tmp/lah
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include "../cuda_forge/block_dplur_jacobian_d.cuh"

static const int NF = 12, FREC = 200, NREC = 160, HF = 150, HN = 75;

static std::vector<double> load(const std::string& path, size_t n) {
    std::vector<double> v(n);
    std::ifstream f(path, std::ios::binary);
    if (!f.read((char*)v.data(), sizeof(double) * n)) { std::fprintf(stderr, "%s を読めない (%zu 個)\n", path.c_str(), n); std::exit(3); }
    return v;
}
static void add_id(float m[5][5], float s) { for (int r = 0; r < 5; ++r) m[r][r] += s; }   // add_identity_scaled と同じ
static void put25(double* o, const float m[5][5]) { for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c) o[r * 5 + c] = m[r][c]; }
static void get25(float m[5][5], const double* in) { for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c) m[r][c] = (float)in[r * 5 + c]; }

struct NodeS { float rho, u, v, w, c, Ht, gamma, roe; };

// 面 1 枚の寄与を D に足す (製品と同じ順序)。conv_inc・visc_inc には増分 (零から) を返す
static void apply_face(float D[5][5], const double* fr, const NodeS& s, bool tp, float conv_inc[5][5], float visc_inc[5][5], float Dthin[5][5], float Kthin[5][5]) {
    const float fa = (float)fr[18], nx = (float)fr[19], ny = (float)fr[20], nz = (float)fr[21];
    const bool has_nbr = fr[4] != 0.0;
    const int br = (int)fr[45];
    const float sdq0[5] = {0, 0, 0, 0, 0};
    float nbr[5] = {};
    for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c) { conv_inc[r][c] = 0; visc_inc[r][c] = 0; Dthin[r][c] = 0; Kthin[r][c] = 0; }
    block_dplur::accumulate_split_jacobian_cf<float>(s.gamma, nx, ny, nz, s.u, s.v, s.w, s.c, s.Ht, tp, fa, has_nbr, sdq0, D, nbr);
    float nb2[5] = {};
    block_dplur::accumulate_split_jacobian_cf<float>(s.gamma, nx, ny, nz, s.u, s.v, s.w, s.c, s.Ht, tp, fa, has_nbr, sdq0, conv_inc, nb2);
    if (!has_nbr) return;                       // node: 境界半割面は粘性の対角を課さない
    const float vd = (float)fr[28];
    if (br == 1) {
        const float cp_i = std::max((float)fr[35], 1.0e-30f), cp_j = std::max((float)fr[36], 1.0e-30f);
        auto thin = [&](float M[5][5], float (*K)[5]) {
            block_dplur::accumulate_thinlayer_visc_jacobian<float>(
                (float)fr[41], (float)fr[42], nx, ny, nz, (float)fr[30],
                s.rho, s.u, s.v, s.w, s.roe, s.gamma, cp_i,
                std::max((float)fr[52], 1.0e-30f), (float)fr[53], (float)fr[54], (float)fr[55], (float)fr[56], (float)fr[51], cp_j,
                fr[43] != 0.0, fr[44] != 0.0, M, K, (int)fr[48]);
        };
        if (fr[49] != 0.0) { add_id(D, vd); add_id(visc_inc, vd); }
        thin(D, nullptr);
        thin(visc_inc, nullptr);
        thin(Dthin, Kthin);
    } else if (br == 3) {
        const int nsr = (int)fr[190];
        const float cfac = (float)fr[47];
        for (int r = 0; r < nsr; ++r) { D[r][r] += vd; visc_inc[r][r] += vd; }
        const float q2 = s.u * s.u + s.v * s.v + s.w * s.w;
        const float e_int = s.roe / s.rho - 0.5f * q2;
        const float t0 = -cfac * (e_int - 0.5f * q2), t1 = -cfac * s.u, t2 = -cfac * s.v, t3 = -cfac * s.w;
        D[4][0] += t0; D[4][1] += t1; D[4][2] += t2; D[4][3] += t3; D[4][4] += cfac;
        visc_inc[4][0] += t0; visc_inc[4][1] += t1; visc_inc[4][2] += t2; visc_inc[4][3] += t3; visc_inc[4][4] += cfac;
    } else if (br == 4) {
        add_id(D, vd); add_id(visc_inc, vd);
    }
}

int main(int argc, char** argv) {
    if (argc != 4) { std::fprintf(stderr, "使い方: %s <dir> <nNodes> <thermallyPerfect>\n", argv[0]); return 3; }
    const std::string dir = argv[1];
    const size_t nn = std::strtoul(argv[2], nullptr, 10);
    const bool tp = std::atoi(argv[3]) != 0;
    if (nn == 0) { std::fprintf(stderr, "節点が 0\n"); return 3; }
    const auto F = load(dir + "/audit_face.f64", nn * NF * FREC);
    const auto N = load(dir + "/audit_node.f64", nn * NREC);
    std::vector<double> of(nn * NF * HF, 0.0), on(nn * HN, 0.0);
    long nface = 0;
    for (size_t k = 0; k < nn; ++k) {
        const double* nd = &N[k * NREC];
        const int nfaces = std::min((int)nd[4], NF);
        const NodeS s{(float)nd[9], (float)nd[10], (float)nd[11], (float)nd[12], (float)nd[13], (float)nd[14], (float)nd[15], (float)nd[16]};
        // 時間項 (生の体積と dt_local から)
        float Dt[5][5] = {};
        const float v = (float)nd[133], dt_l = (float)nd[134];
        add_id(Dt, v / std::max(dt_l, 1.0e-30f));
        add_id(Dt, v * (float)nd[135]);
        put25(&on[k * HN + 0], Dt);
        // 面: 記録の前の D から (局所) と、時間項の後の D から連ねて (全体)
        float Dchain[5][5]; get25(Dchain, &nd[24]);
        for (int sl = 0; sl < nfaces; ++sl) {
            const double* fr = &F[(k * NF + sl) * FREC];
            double* o = &of[(k * NF + sl) * HF];
            ++nface;
            float Dloc[5][5], ci[5][5], vi[5][5], dth[5][5], kth[5][5];
            get25(Dloc, &fr[57]);
            apply_face(Dloc, fr, s, tp, ci, vi, dth, kth);
            float c2[5][5], v2[5][5], d2[5][5], k2[5][5];
            apply_face(Dchain, fr, s, tp, c2, v2, d2, k2);
            put25(o + 0, Dloc); put25(o + 25, Dchain);
            for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c)
                o[50 + r * 5 + c] = std::max({std::fabs(fr[57 + r * 5 + c]), std::fabs(fr[157 + r * 5 + c]), (double)std::fabs(ci[r][c]), (double)std::fabs(vi[r][c])});
            if (fr[5] != 0.0) {   // ライン面: 対流の K を列の抽出で
                const float fa = (float)fr[18], nx = (float)fr[19], ny = (float)fr[20], nz = (float)fr[21];
                for (int j = 0; j < 5; ++j) {
                    float dd[5][5] = {}, kc[5] = {}, ev[5] = {0, 0, 0, 0, 0}; ev[j] = fa;
                    block_dplur::accumulate_split_jacobian_cf<float>(s.gamma, nx, ny, nz, s.u, s.v, s.w, s.c, s.Ht, tp, fa, true, ev, dd, kc);
                    for (int i = 0; i < 5; ++i) o[75 + i * 5 + j] = kc[i];
                }
            }
            put25(o + 100, dth); put25(o + 125, kth);
        }
        // 軸対称 (r 重み、枝 1): 記録された面のループの後の D から
        float Da[5][5]; get25(Da, &nd[49]);
        float sc[5][5]; for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c) sc[r][c] = std::fabs(Da[r][c]);
        if ((int)nd[124] == 1) {
            const float A = (float)nd[125], hoop = (float)nd[126], al = (float)nd[127];
            const float g1 = s.gamma - 1.0f, q2 = s.u * s.u + s.v * s.v + s.w * s.w;
            const float t0 = -A * (0.5f * g1 * q2 + hoop * s.v), t1 = A * (g1 * s.u), t2 = A * (g1 * s.v + hoop), t3 = A * (g1 * s.w), t4 = -A * g1, t5 = al * A * s.c;
            Da[2][0] += t0; Da[2][1] += t1; Da[2][2] += t2; Da[2][3] += t3; Da[2][4] += t4; Da[2][2] += t5;
            sc[2][0] = std::max(sc[2][0], std::fabs(t0)); sc[2][1] = std::max(sc[2][1], std::fabs(t1)); sc[2][2] = std::max({sc[2][2], std::fabs(t2), std::fabs(t5)});
            sc[2][3] = std::max(sc[2][3], std::fabs(t3)); sc[2][4] = std::max(sc[2][4], std::fabs(t4));
        }
        for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c) sc[r][c] = std::max(sc[r][c], std::fabs(Da[r][c]));
        put25(&on[k * HN + 25], Da);
        for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c) on[k * HN + 50 + r * 5 + c] = sc[r][c];
    }
    std::ofstream f1(dir + "/audit_host_face.f64", std::ios::binary); f1.write((const char*)of.data(), sizeof(double) * of.size());
    std::ofstream f2(dir + "/audit_host_node.f64", std::ios::binary); f2.write((const char*)on.data(), sizeof(double) * on.size());
    if (!f1 || !f2) { std::fprintf(stderr, "出力を書けない\n"); return 3; }
    std::printf("line_audit_helper: 節点 %zu、面 %ld を再現した\n", nn, nface);
    return 0;
}
