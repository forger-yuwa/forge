// =============================================================================
// transport_table_ab_host.cpp — 輸送表の区間選択の判別 A/B (host; plan thermophysics-solver-owned-species-db #5t2-3)
//   codex diagnose 2026-09-27 (notes/reviews/2026-09-27-transport-tables-diagnose.md) の判別:
//   同じ分割表・係数・float 演算で、区間の選び方だけを変える。
//     A = float の ln T で選ぶ (ln T を float に丸め、境界の ln b も float に丸めて比べる) — この試験ハーネスにだけある
//     B = 元の T で選ぶ (ソルバの transport_tab_species = cuda_forge/transportTables_d.cuh)
//   対象は単成分 fit: T ≤ 1000 K で μ 1e-5 Pa s・λ 0.02 W/(m K)、1000 K 超で μ 1.01e-5・λ 0.0202 (1 % の段差)。
//   T = 1000 K と両隣の float (と少し離れた点) を評価し、1 行の JSON を出す。判定は test_transport_gpu.py が独立参照で行う。
//
//   g++ -O1 -std=c++17 -I solver_density_cuda solver_density_cuda/tests/unit/transport_table_ab_host.cpp -o transport_table_ab_host
// =============================================================================
#include "cuda_forge/transportTables_d.cuh"

#include <cmath>
#include <cstdio>
#include <vector>

// A: float の ln T で区間を選ぶ (試験専用。ソルバには置かない)
static void species_select_by_float_lnT(const TransportTablesHost& H, const TransportTablesF& tb, int r, float T,
                                        float* mu, float* lam)
{
    const TransportTabRefF& ref = tb.spTab[r];
    const float lnTf = logf(T);
    int k = ref.seg0;
    const int e = ref.seg0 + ref.nseg - 1;
    while (k < e && lnTf > (float)log(H.segHi[r][k - ref.seg0])) ++k;
    float u;
    const int j = transport_tab_locate(tb.seg[k], (double)lnTf, &u);
    *mu  = transport_tab_hermite(tb.spc[2*j],     u);
    *lam = transport_tab_hermite(tb.spc[2*j + 1], u);
}

int main()
{
    SpeciesTransportD s{};
    s.model = TRANSPORT_MODEL_FIT;
    s.MW = 0.028;
    for (TransportFitD* f : {&s.V, &s.C}) {
        f->n = 2;
        f->Tlo[0] = 150.0;  f->Thi[0] = 1000.0;
        f->Tlo[1] = 1000.0; f->Thi[1] = 20000.0;
        for (int k = 0; k < 2; ++k) { f->A[k] = 0.0; f->B[k] = 0.0; f->C[k] = 0.0; }
    }
    s.V.D[0] = log(100.0); s.V.D[1] = log(101.0);   // μP: 1e-5 / 1.01e-5 Pa s
    s.C.D[0] = log(200.0); s.C.D[1] = log(202.0);   // μW/(cm K): 0.02 / 0.0202 W/(m K)

    TransportTablesHost H;
    if (!transport_tables_build_host({s}, {}, {s.MW}, {1.0}, H)) {
        std::printf("{\"ok\": false, \"error\": \"%s\"}\n", H.error.c_str());
        return 2;
    }
    const TransportTablesF tb = transport_tables_view_host(H);
    std::vector<float> Ts;
    const float b = 1000.0f;
    Ts.push_back(std::nextafter(std::nextafter(b, 0.0f), 0.0f));
    Ts.push_back(std::nextafter(b, 0.0f));
    Ts.push_back(b);
    Ts.push_back(std::nextafter(b, 2000.0f));
    Ts.push_back(std::nextafter(std::nextafter(b, 2000.0f), 2000.0f));
    Ts.push_back(999.5f);
    Ts.push_back(1000.5f);
    std::printf("{\"ok\": true, \"segments\": %d, \"points\": [", tb.spTab[0].nseg);
    for (size_t i = 0; i < Ts.size(); ++i) {
        const float T = Ts[i];
        float muA, laA, muB, laB;
        species_select_by_float_lnT(H, tb, 0, T, &muA, &laA);
        transport_tab_species(tb, 0, T, log((double)T), &muB, &laB);
        std::printf("%s{\"T\": %.9g, \"lnT_float\": %.9g, \"muA\": %.9g, \"lamA\": %.9g, \"muB\": %.9g, \"lamB\": %.9g}",
                    i ? ", " : "", (double)T, (double)logf(T), (double)muA, (double)laA, (double)muB, (double)laB);
    }
    std::printf("]}\n");
    return 0;
}
