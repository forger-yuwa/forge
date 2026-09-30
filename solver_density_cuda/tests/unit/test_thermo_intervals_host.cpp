// =============================================================================
// test_thermo_intervals_host.cpp — 区間可変 (plan thermophysics-solver-owned-species-db §5.1 #13-1) の host 試験 (GPU 不要)
//
//   (C) G1-c = §6 V3: 区切りの違う種を畳んだ lump (区切りの和集合) が、構成種ごとの質量分率加重和と 100–20000 K の全点
//       (各区切り温度そのものと ±1e-9 K、構成種の外挿域 (< Tlo・> Thi) を含む) で相対 1e-12 かつ絶対 cp 1e-9 J/(kg K)・
//       h 1e-6 J/kg・s° 1e-9 J/(kg K) 以内 (h の相対は範囲の max|h| で規格化; test_species_lump_solver.py と同じ)。datum 0 / 298.15 K。
//       6000 K 超 (第 3 区間を持つ種) の絶対許容は同じ係数の long double 評価との差 (合成だけの誤差) で判定する (plan §6 V3 の再スコープ 2026-10-01)。
//       区切りでの h・cp の段差 (上の区間 − 下の区間を同じ T で評価) は、構成種の段差の加重和からの増分が同じ絶対許容差以内。
//       試験種: CEA thermo.inp の N2 3 区間 (200/1000/6000/20000; 外部 DB)、Tmid 1500 の 2 区間種、Tlo 298.15 の種
//       (N2 の係数を 298.15/1000/6000/20000 に置いた 3 区間型と、298.15/1000/6000 の 2 区間型)、CEA の e- (298.15–20000, 3 区間)。
//   (F) G1-f: 3 区間種を含む記録を書き → speciesDB_diffRecord で読み直し (差なし・互換性ハッシュ再計算一致)、
//       schema と外挿規約が *_nint、2 区間だけの記録は従来の schema・外挿規約で、本文に Tbounds / coef[ が現れない。
//       lump の合成規約は区切りが揃えば従来の文字列 (SPECIES_LUMP_SYNTHESIS)、違えば SPECIES_LUMP_SYNTHESIS_UNION。
//   (R) 上限超過の拒否: 外部 DB の 4 区間種、和集合が 4 区間になる lump (N2 + e-)、境界の非単調・書式の混在。
//
// ビルド/実行 (共通データの埋め込みヘッダを先に生成する):
//   cmake -DIN=solver_density_cuda/data/species/forge_species_v1.yaml -DOUT=/tmp/forge_species_gen/forge_species_data.hpp
//       -P solver_density_cuda/cmake/embed_species_data.cmake
//   g++ -O1 -std=c++17 -I solver_density_cuda -I /tmp/forge_species_gen solver_density_cuda/tests/unit/test_thermo_intervals_host.cpp
//       solver_density_cuda/input/speciesDB.cpp solver_density_cuda/input/speciesTransportDB.cpp -lyaml-cpp -o /tmp/test_thermo_intervals_host
//   /tmp/test_thermo_intervals_host [出力ディレクトリ]   (出力ディレクトリに試験用の外部 DB と記録を書く; 既定 /tmp/test_thermo_intervals_host.d)
//   規約: [PASS]/[FAIL]、失敗があれば非ゼロ終了。
// =============================================================================
#include "input/speciesDB.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <set>
#include <string>
#include <vector>

namespace fs = std::filesystem;
static int g_fail = 0;
static void check(bool ok, const std::string& what)
{
    std::printf("%s %s\n", ok ? "[PASS]" : "[FAIL]", what.c_str());
    if (!ok) ++g_fail;
}

// CEA thermo.inp (/home/sano/work/forge/.venv-cea/nasa_cea/thermo.inp; SHA-256 b1bc0707...) の N2 と e-
static const double N2a[3][9] = {
    {22103.71497, -381.846182, 6.08273836, -0.00853091441, 1.384646189e-05, -9.62579362e-09, 2.519705809e-12, 710.846086, -10.76003744},
    {587712.406, -2239.249073, 6.06694922, -0.00061396855, 1.491806679e-07, -1.923105485e-11, 1.061954386e-15, 12832.10415, -15.86640027},
    {831013916.0, -642073.354, 202.0264635, -0.03065092046, 2.486903333e-06, -9.70595411e-11, 1.437538881e-15, 4938707.04, -1672.09974}};
static const double Ea[9] = {0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, -745.375, -11.72081224};

static std::string g17(double v) { char b[64]; std::snprintf(b, sizeof(b), "%.17g", v); return b; }
static std::string arr(const double* a, int n) { std::string s = "["; for (int i = 0; i < n; ++i) s += (i ? ", " : "") + g17(a[i]); return s + "]"; }

static void writeDb(const fs::path& p, bool withBad = false)
{
    std::ofstream f(p);
    // 3 区間 (CEA N2 そのもの; 区間可変の書式)
    f << "N2CEA3:\n  MW: 0.0280134\n  LJ_sigma: 3.621\n  LJ_eps_kB: 97.53\n  Tbounds: [200.0, 1000.0, 6000.0, 20000.0]\n  nasa9_intervals:\n";
    for (int k = 0; k < 3; ++k) f << "    - " << arr(N2a[k], 9) << "\n";
    // Tlo 298.15 の 3 区間型 (N2 の係数を 298.15 から)
    f << "N2T298:\n  MW: 0.0280134\n  LJ_sigma: 3.621\n  LJ_eps_kB: 97.53\n  Tbounds: [298.15, 1000.0, 6000.0, 20000.0]\n  nasa9_intervals:\n";
    for (int k = 0; k < 3; ++k) f << "    - " << arr(N2a[k], 9) << "\n";
    // Tlo 298.15 の 2 区間型 (従来の書式)
    f << "N2T298B:\n  MW: 0.0280134\n  LJ_sigma: 3.621\n  LJ_eps_kB: 97.53\n  Tlo: 298.15\n  Tmid: 1000.0\n  Thi: 6000.0\n"
      << "  nasa9_low: " << arr(N2a[0], 9) << "\n  nasa9_high: " << arr(N2a[1], 9) << "\n";
    // Tmid 1500 の 2 区間型 (区切り以外は N2 の係数; 1500 K に段差がある)
    f << "TMID1500:\n  MW: 0.0280134\n  LJ_sigma: 3.621\n  LJ_eps_kB: 97.53\n  Tlo: 200.0\n  Tmid: 1500.0\n  Thi: 6000.0\n"
      << "  nasa9_low: " << arr(N2a[0], 9) << "\n  nasa9_high: " << arr(N2a[1], 9) << "\n";
    // 2 区間 200/1000/6000 の N2・O2 (内蔵の先頭 2 区間の写し)。段 3 (#13-3) から内蔵 N2/O2/AR は CEA そのもの (3 区間) なので、
    // 「区切りの違う 2 区間種と畳む」試験には 2 区間の写しを外部 DB で与える (内蔵そのものの lump は V2-type の case で見る)
    for (const char* nm : {"N2", "O2"}) {
        const SpeciesThermo b = speciesDB_builtin().at(nm);
        f << nm << "B2:\n  MW: " << g17(b.MW) << "\n  LJ_sigma: " << g17(b.sigma_LJ) << "\n  LJ_eps_kB: " << g17(b.eps_kB)
          << "\n  Tlo: 200.0\n  Tmid: 1000.0\n  Thi: 6000.0\n  nasa9_low: " << arr(b.coef[0], 9) << "\n  nasa9_high: " << arr(b.coef[1], 9) << "\n";
    }
    // CEA の e- (298.15/1000/6000/20000; 3 区間とも同じ係数)
    f << "Eminus:\n  MW: 5.48579903e-07\n  LJ_sigma: 3.0\n  LJ_eps_kB: 10.0\n  Tbounds: [298.15, 1000.0, 6000.0, 20000.0]\n  nasa9_intervals:\n";
    for (int k = 0; k < 3; ++k) f << "    - " << arr(Ea, 9) << "\n";
    if (withBad) {
        f << "FOURINT:\n  MW: 0.0280134\n  Tbounds: [200.0, 1000.0, 3000.0, 6000.0, 20000.0]\n  nasa9_intervals:\n";
        for (int k = 0; k < 4; ++k) f << "    - " << arr(N2a[k % 3], 9) << "\n";
        f << "NONMONO:\n  MW: 0.0280134\n  Tbounds: [200.0, 6000.0, 1000.0, 20000.0]\n  nasa9_intervals:\n";
        for (int k = 0; k < 3; ++k) f << "    - " << arr(N2a[k], 9) << "\n";
        f << "MIXEDFMT:\n  MW: 0.0280134\n  Tmid: 1000.0\n  Tbounds: [200.0, 1000.0, 6000.0]\n  nasa9_intervals:\n";
        for (int k = 0; k < 2; ++k) f << "    - " << arr(N2a[k], 9) << "\n";
    }
}

// 区間 k の係数で強制評価 (段差の検査用; 範囲の内側だけで使う)
static double hk(const SpeciesThermo& s, int k, double T)
{
    const double* a = s.coef[k];
    const double Ti = 1.0/T, lnT = std::log(T);
    return THERMO_RU*T*(-a[0]*Ti*Ti + a[1]*lnT*Ti + a[2] + a[3]*T/2.0 + a[4]*T*T/3.0 + a[5]*T*T*T/4.0 + a[6]*T*T*T*T/5.0 + a[7]*Ti)/s.MW;
}
static double cpk(const SpeciesThermo& s, int k, double T)
{
    const double* a = s.coef[k];
    const double Ti = 1.0/T;
    return THERMO_RU*(a[0]*Ti*Ti + a[1]*Ti + a[2] + a[3]*T + a[4]*T*T + a[5]*T*T*T + a[6]*T*T*T*T)/s.MW;
}
// 種 s の T での段差 (上の区間 − 下の区間)。T が s の内側の区切りでなければ 0。
static void jumpAt(const SpeciesThermo& s, double T, double* dh, double* dcp)
{
    *dh = 0.0; *dcp = 0.0;
    for (int k = 0; k + 1 < s.nInt; ++k)
        if (s.Tbrk[k] == T) { *dh = hk(s, k+1, T) - hk(s, k, T); *dcp = cpk(s, k+1, T) - cpk(s, k, T); }
}

// long double で同じ規約 (区間選択・端でのクランプ・線形/対数外挿) を評価する (係数は double のまま)。
// 合成 (係数の和) そのものの誤差を double の評価の丸め (CEA 6000–20000 K 区間の N2 単体で ~1.2e-6 J/kg) から切り分ける参考値用。
static void evalLD(const SpeciesThermo& s, double Td, long double* cp, long double* h, long double* s0)
{
    auto poly = [&](long double T, long double* c, long double* hh, long double* ss) {
        const double* a = s.coef[thermo_interval(s, (double)T)];
        const long double Ti = 1.0L/T, lnT = std::log(T);
        *c  = (long double)THERMO_RU*(a[0]*Ti*Ti + a[1]*Ti + a[2] + a[3]*T + a[4]*T*T + a[5]*T*T*T + a[6]*T*T*T*T);
        *hh = (long double)THERMO_RU*T*(-a[0]*Ti*Ti + a[1]*lnT*Ti + a[2] + a[3]*T/2 + a[4]*T*T/3 + a[5]*T*T*T/4 + a[6]*T*T*T*T/5 + a[7]*Ti);
        *ss = (long double)THERMO_RU*(-a[0]*Ti*Ti/2 - a[1]*Ti + a[2]*lnT + a[3]*T + a[4]*T*T/2 + a[5]*T*T*T/3 + a[6]*T*T*T*T/4 + a[8]);
    };
    const long double T = Td;
    const long double Tb = (Td < s.Tlo) ? s.Tlo : ((Td > s.Thi) ? s.Thi : Td);
    long double c, hh, ss;
    poly(Tb, &c, &hh, &ss);
    if (T != Tb) { hh += c*(T - Tb); ss += c*std::log(T/Tb); }
    *cp = c/s.MW; *h = hh/s.MW; *s0 = ss/s.MW;
}

static void applyDatum(SpeciesThermo& s, double Tref)
{
    if (Tref > 0.0) thermo_add_a7(s, -thermo_h_molar(s, Tref)/THERMO_RU);
}

// lump db.species[li] (構成は db.lumps[li]) を構成種の質量分率加重和と比べる
static void checkLump(const ResolvedSpeciesDB& db, int li, const std::string& tag, const std::string& expectSynth, int expectN)
{
    const ResolvedLump& L = db.lumps[li];
    check(L.synthesis == expectSynth, tag + ": synthesis rule '" + L.synthesis.substr(0, 60) + "...'");
    check(db.species[li].nInt == expectN, tag + ": synthesized nInt " + std::to_string(db.species[li].nInt) + " (expected " + std::to_string(expectN) + ")");
    for (double Tref : {0.0, 298.15}) {
        SpeciesThermo lump = db.species[li];
        applyDatum(lump, Tref);
        std::vector<SpeciesThermo> m = L.memberSpecies;
        for (auto& x : m) applyDatum(x, Tref);
        std::vector<double> Y(m.size());
        for (size_t k = 0; k < m.size(); ++k) Y[k] = L.x[k]*m[k].MW/lump.MW;
        // 評価点: 100–20000 K の 1 K 刻み (< 1000 K) と 10 K 刻み、全区切り (和集合と構成種) とその ±1e-9 K (20000 K の +1e-9 は lump の外挿)
        std::set<double> Ts;
        for (double T = 100.0; T < 1000.0; T += 1.0) Ts.insert(T);
        for (double T = 1000.0; T <= 20000.0; T += 10.0) Ts.insert(T);
        std::vector<double> brk;
        for (int k = 0; k <= lump.nInt; ++k) brk.push_back(thermo_bound(lump, k));
        for (const auto& x : m) for (int k = 0; k <= x.nInt; ++k) brk.push_back(thermo_bound(x, k));
        for (double b : brk) { Ts.insert(b); Ts.insert(b - 1e-9); Ts.insert(b + 1e-9); }
        // 判定の分け方 (plan §6 V3 の再スコープ 2026-10-01, diagnostician; 事前固定):
        //   (a) T <= 6000 K: double 同士 (lump の double 評価 vs 構成種の double 評価の加重和) で相対 1e-12・絶対 cp 1e-9 / h 1e-6 / s° 1e-9。
        //   (b) T > 6000 K で第 3 区間を持つ種 (lump か構成種の Thi > 6000 K) のとき: 相対 1e-12 (h は範囲 max|h| 規格化) は double 同士のまま、
        //       絶対は**同じ係数の long double 評価**どうしの差 (合成だけの誤差) で |dh| <= 1e-6・|dcp| <= 1e-9・|ds°| <= 1e-9。
        //       double 同士の |dh| と単種の double 丸め (|double - long double|) は info (合否に使わない)。
        //   第 3 区間を持たない lump (Thi <= 6000 K) は 6000 K 超も外挿域として (a) で判定する (再スコープの対象外)。
        bool hi3 = lump.Thi > 6000.0;
        for (const auto& x : m) hi3 = hi3 || x.Thi > 6000.0;
        double hmax = 0.0;
        for (double T : Ts) hmax = std::max(hmax, std::fabs(thermo_h_mass(lump, T)));
        double wcp = 0, wh = 0, ws = 0;                    // 相対 (全点)
        double acp = 0, ah = 0, as = 0, Tah = 0;           // 絶対 double 同士 (a の点)
        double lcp = 0, lh = 0, ls = 0, Tlh = 0;           // 絶対 long double 同士 = 合成だけの誤差 (b の点)
        double ihDD = 0, iRound = 0, TiR = 0;              // info (b の点): double 同士の |dh|、単種の double 丸め |h_double - h_ld|
        size_t na = 0, nb = 0;
        for (double T : Ts) {
            double cp = 0, h = 0, s = 0;
            for (size_t k = 0; k < m.size(); ++k) { cp += Y[k]*thermo_cp_mass(m[k], T); h += Y[k]*thermo_h_mass(m[k], T); s += Y[k]*thermo_s0_mass(m[k], T); }
            const double cpl = thermo_cp_mass(lump, T), hl = thermo_h_mass(lump, T), sl = thermo_s0_mass(lump, T);
            wcp = std::max(wcp, std::fabs(cpl - cp)/std::fabs(cp));
            wh  = std::max(wh,  std::fabs(hl - h)/hmax);
            ws  = std::max(ws,  std::fabs(sl - s)/std::fabs(s));
            if (!(hi3 && T > 6000.0)) {
                ++na;
                acp = std::max(acp, std::fabs(cpl - cp));
                if (std::fabs(hl - h) > ah) { ah = std::fabs(hl - h); Tah = T; }
                as  = std::max(as,  std::fabs(sl - s));
            } else {
                ++nb;
                long double c, hh, ss, cS = 0, hS = 0, sS = 0;
                for (size_t k = 0; k < m.size(); ++k) {
                    evalLD(m[k], T, &c, &hh, &ss); cS += Y[k]*c; hS += Y[k]*hh; sS += Y[k]*ss;
                    const double r = (double)std::fabs((long double)thermo_h_mass(m[k], T) - hh);
                    if (r > iRound) { iRound = r; TiR = T; }
                }
                evalLD(lump, T, &c, &hh, &ss);
                lcp = std::max(lcp, (double)std::fabs(c - cS));
                if ((double)std::fabs(hh - hS) > lh) { lh = (double)std::fabs(hh - hS); Tlh = T; }
                ls  = std::max(ls,  (double)std::fabs(ss - sS));
                ihDD = std::max(ihDD, std::fabs(hl - h));
                const double r = (double)std::fabs((long double)hl - hh);
                if (r > iRound) { iRound = r; TiR = T; }
            }
        }
        char b[640];
        std::snprintf(b, sizeof(b), "%s datum %g: lump = sum of constituents, %zu points 100-20000 K (breakpoints +-1e-9 K): "
                      "rel cp %.2e h %.2e s %.2e (<=1e-12); (a) %zu pts%s double vs double abs cp %.2e h %.2e (at %.6g K) s %.2e (<=1e-9/1e-6/1e-9)",
                      tag.c_str(), Tref, Ts.size(), wcp, wh, ws, na, hi3 ? " T<=6000 K" : " (no 3rd interval: all T)", acp, ah, Tah, as);
        check(wcp <= 1e-12 && wh <= 1e-12 && ws <= 1e-12 && acp <= 1e-9 && ah <= 1e-6 && as <= 1e-9, b);
        if (hi3) {
            std::snprintf(b, sizeof(b), "%s datum %g: (b) %zu pts T>6000 K synthesis-only error (long double of the same coefficients): "
                          "|dcp| %.2e |dh| %.2e (at %.6g K) |ds| %.2e (<=1e-9/1e-6/1e-9)", tag.c_str(), Tref, nb, lcp, lh, Tlh, ls);
            check(lcp <= 1e-9 && lh <= 1e-6 && ls <= 1e-9, b);
            std::printf("       (info, T>6000 K) double vs double |dh| %.2e J/kg; single-species double rounding |h_double - h_longdouble| up to %.2e J/kg (at %.6g K)\n",
                        ihDD, iRound, TiR);
        }
        // 区切りでの段差の増分
        double dJh = 0, dJcp = 0, maxJ = 0;
        for (int k = 0; k + 1 < lump.nInt; ++k) {
            const double T = lump.Tbrk[k];
            const double jl = hk(lump, k+1, T) - hk(lump, k, T), jcl = cpk(lump, k+1, T) - cpk(lump, k, T);
            double jm = 0, jcm = 0;
            for (size_t i = 0; i < m.size(); ++i) { double a, c; jumpAt(m[i], T, &a, &c); jm += Y[i]*a; jcm += Y[i]*c; }
            dJh = std::max(dJh, std::fabs(jl - jm)); dJcp = std::max(dJcp, std::fabs(jcl - jcm)); maxJ = std::max(maxJ, std::fabs(jm));
        }
        std::snprintf(b, sizeof(b), "%s datum %g: breakpoint step increment vs weighted constituent steps: |dh| %.2e J/kg (<=1e-6), "
                      "|dcp| %.2e (<=1e-9); constituent step max %.3e J/kg", tag.c_str(), Tref, dJh, dJcp, maxJ);
        check(dJh <= 1e-6 && dJcp <= 1e-9, b);
    }
}

static bool has(const std::string& s, const std::string& k) { return s.find(k) != std::string::npos; }

int main(int argc, char** argv)
{
    const fs::path dir = (argc > 1) ? fs::path(argv[1]) : fs::path("/tmp/test_thermo_intervals_host.d");
    fs::create_directories(dir);
    const fs::path db = dir / "species_db_intervals.yaml";
    writeDb(db);
    const std::vector<SpeciesLumpSpec> none;

    // ---- 外部 DB の区間可変の書式を読む ----
    {
        const ResolvedSpeciesDB r = speciesDB_resolve({"N2CEA3", "N2T298", "N2T298B", "TMID1500", "Eminus"}, db.string());
        check(r.species[0].nInt == 3 && r.species[0].Tlo == 200.0 && r.species[0].Tbrk[0] == 1000.0 && r.species[0].Tbrk[1] == 6000.0
              && r.species[0].Thi == 20000.0 && r.species[0].coef[2][0] == N2a[2][0] && r.species[0].coef[2][8] == N2a[2][8],
              "external DB Tbounds/nasa9_intervals: N2CEA3 has 3 intervals 200/1000/6000/20000 with the CEA coefficients");
        check(r.species[2].nInt == 2 && r.species[2].Tlo == 298.15 && r.species[3].Tbrk[0] == 1500.0, "external DB 2-interval format (Tlo 298.15, Tmid 1500)");
        // 3 区間種の評価: 区間ちょうどは上の区間、Thi (20000) の外は外挿
        const SpeciesThermo& s = r.species[0];
        check(thermo_interval(s, 999.9999999) == 0 && thermo_interval(s, 1000.0) == 1 && thermo_interval(s, 6000.0) == 2
              && thermo_interval(s, 20000.0) == 2 && thermo_interval(s, std::nan("")) == 2,
              "interval selection: breakpoint belongs to the upper interval, NaN -> last interval");
    }

    // ---- (C) V3: 区切りの違う種を畳む ----
    struct LumpCase { std::string tag; std::vector<std::string> mem; std::vector<double> fr; std::string synth; int n; };
    const std::vector<LumpCase> cases = {
        {"V3 N2B2 (2-int) + N2CEA3 (3-int 200/1000/6000/20000)", {"N2B2", "N2CEA3"}, {0.6, 0.4}, SPECIES_LUMP_SYNTHESIS_UNION, 3},
        {"V3 O2B2 + TMID1500 (Tmid 1000 vs 1500)", {"O2B2", "TMID1500"}, {0.3, 0.7}, SPECIES_LUMP_SYNTHESIS_UNION, 3},
        {"V3 O2B2 + N2T298B (Tlo 200 vs 298.15)", {"O2B2", "N2T298B"}, {0.5, 0.5}, SPECIES_LUMP_SYNTHESIS_UNION, 3},
        {"V3 N2T298 + Eminus (equal breakpoints 298.15/1000/6000/20000)", {"N2T298", "Eminus"}, {0.999, 0.001}, SPECIES_LUMP_SYNTHESIS, 3},
        // 内蔵そのもの: 段 3 (#13-3) から 3 区間 200/1000/6000/20000 で揃う (以前は 2 区間)
        {"V2-type N2 + O2 + AR (builtin, equal breakpoints 200/1000/6000/20000)", {"N2", "O2", "AR"}, {0.78, 0.21, 0.01}, SPECIES_LUMP_SYNTHESIS, 3},
    };
    for (const auto& c : cases) {
        SpeciesLumpSpec lp; lp.name = "LMP"; lp.basis = "mole"; lp.members = c.mem; lp.fractions = c.fr;
        ResolvedSpeciesDB r;
        try {
            r = speciesDB_resolve({"LMP", "H2O"}, db.string(), {lp});
        } catch (const std::exception& e) {
            check(false, c.tag + ": resolves (" + e.what() + ")");
            continue;
        }
        checkLump(r, 0, c.tag, c.synth, c.n);
    }

    // ---- (F) G1-f: 記録の書き → 読み → 再ハッシュ ----
    {
        SpeciesLumpSpec lp; lp.name = "LMP"; lp.basis = "mole"; lp.members = {"O2B2", "TMID1500"}; lp.fractions = {0.3, 0.7};
        const ResolvedSpeciesDB r = speciesDB_resolve({"LMP", "N2CEA3", "H2O"}, db.string(), {lp});
        for (double Tref : {0.0, 298.15}) {
            const fs::path d = dir / ("rec_nint_" + std::to_string(static_cast<int>(Tref)));
            fs::remove_all(d); fs::create_directories(d);
            const SpeciesRecordInfo info = speciesDB_writeRecord(r, Tref, db.string(), "", "not_checked_resolve_only", 0, d.string());
            const std::string p = (d / info.recordFile).string();
            std::ifstream f(p); std::string text((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
            const std::vector<std::string> diff = speciesDB_diffRecord(p, r, Tref);
            check(diff.empty(), "G1-f C++: 3-interval record written and read back without differences (Tref " + g17(Tref) + ")"
                  + (diff.empty() ? std::string() : ": " + diff.front()));
            check(has(text, "schema: \"" SPECIES_RECORD_SCHEMA_NINT "\"") && has(text, SPECIES_RECORD_EXTRAPOLATION_NINT)
                  && has(text, "Tbounds: [200, 1000, 6000, 20000]") && has(text, "nasa9_intervals:") && has(text, SPECIES_LUMP_SYNTHESIS_UNION),
                  "G1-f C++: record carries the *_nint schema, the n-interval convention, Tbounds/nasa9_intervals and the union synthesis");
            const std::string ct = speciesDB_compatText(r, Tref);
            check(has(ct, "species[1].T: 200 1000 6000 20000\n") && has(ct, "species[1].coef[2]: 831013916") && has(ct, "species[0].T: 200 1000 1500 6000\n")
                  && has(ct, "species[2].low:"), "G1-f C++: compat text has all bounds and coef[k] lines for nInt != 2, low/high for H2O");
        }
        // 係数 1 つを書き換えた記録は差 (区間可変のキー名) と自己整合の破れを示す
        {
            const fs::path d = dir / "rec_nint_0";
            std::string p;
            for (const auto& e : fs::directory_iterator(d)) p = e.path().string();
            std::ifstream f(p); std::string text((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
            const std::string key = "831013916";
            const size_t at = text.find(key);
            text.replace(at, key.size(), "831013917");
            const std::string q = (d / "edited.yaml").string();
            std::ofstream(q) << text;
            const std::vector<std::string> diff = speciesDB_diffRecord(q, r, 0.0);
            bool self = false, coef = false;
            for (const auto& x : diff) { self = self || has(x, "compat_hash"); coef = coef || has(x, "N2CEA3.nasa9_intervals[2][0]"); }
            check(self && coef, "G1-f C++: edited 3rd-interval coefficient -> self-consistency broken and N2CEA3.nasa9_intervals[2][0] shown");
        }
        // 2 区間だけの記録: 従来の schema・外挿規約、本文に区間可変の行が無い
        const ResolvedSpeciesDB r2 = speciesDB_resolve({"N2B2", "H2O", "TMID1500"}, db.string());   // 内蔵 N2 は段 3 から 3 区間
        const std::string ct2 = speciesDB_compatText(r2, 298.15);
        check(has(ct2, "schema: " SPECIES_RECORD_SCHEMA "\n") && has(ct2, "extrapolation: " SPECIES_RECORD_EXTRAPOLATION "\n")
              && !has(ct2, "coef[") && !has(ct2, "_nint"), "G1-f C++: 2-interval-only record keeps schema v1 and the 2-interval convention");
    }

    // ---- (R) 上限超過・不正の拒否 ----
    {
        const fs::path dbBad = dir / "species_db_bad.yaml";
        writeDb(dbBad, true);
        auto refused = [&](const std::vector<std::string>& names, const std::vector<SpeciesLumpSpec>& lumps, const std::string& key, const std::string& tag) {
            std::string what;
            try { speciesDB_resolve(names, dbBad.string(), lumps); } catch (const std::exception& e) { what = e.what(); }
            check(!what.empty() && has(what, key), tag + " -> refused: " + what.substr(0, 200));
        };
        refused({"N2"}, none, "THERMO_MAX_INTERVALS", "external DB species with 4 intervals (MAXI+1)");
        SpeciesLumpSpec lp; lp.name = "LMP"; lp.basis = "mole"; lp.members = {"N2", "Eminus"}; lp.fractions = {0.99, 0.01};
        std::string what;
        try { speciesDB_resolve({"LMP"}, db.string(), {lp}); } catch (const std::exception& e) { what = e.what(); }
        check(has(what, "makes 4 intervals") && has(what, "THERMO_MAX_INTERVALS"), "lump N2 + e- (union 200/298.15/1000/6000/20000 = 4 intervals) -> refused: " + what.substr(0, 160));
    }
    {
        // 非単調・書式の混在 (各 1 種だけの DB で)
        for (const std::string nm : {"NONMONO", "MIXEDFMT"}) {
            const fs::path p = dir / ("species_db_" + nm + ".yaml");
            std::ofstream f(p);
            if (nm == "NONMONO") {
                f << "NONMONO:\n  MW: 0.0280134\n  Tbounds: [200.0, 6000.0, 1000.0, 20000.0]\n  nasa9_intervals:\n";
                for (int k = 0; k < 3; ++k) f << "    - " << arr(N2a[k], 9) << "\n";
            } else {
                f << "MIXEDFMT:\n  MW: 0.0280134\n  Tmid: 1000.0\n  Tbounds: [200.0, 1000.0, 6000.0]\n  nasa9_intervals:\n";
                for (int k = 0; k < 2; ++k) f << "    - " << arr(N2a[k], 9) << "\n";
            }
            f.close();
            std::string what;
            try { speciesDB_resolve({"N2"}, p.string()); } catch (const std::exception& e) { what = e.what(); }
            check(!what.empty(), "external DB " + nm + " -> refused: " + what.substr(0, 160));
        }
    }
    std::printf(g_fail ? "FAILED (%d failures)\n" : "ALL PASSED (0 failures)\n", g_fail);
    return g_fail ? 1 : 0;
}
