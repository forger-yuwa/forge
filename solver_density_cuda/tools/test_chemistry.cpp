// =============================================================================
// test_chemistry.cpp — 有限速度化学 (chemistry_d.cuh) のホスト単体検証
//   (1) 解析 Jacobian ∂ω/∂(ρY) と ∂ω/∂T を有限差分と照合
//   (2) 0-D 定積断熱反応器 (BDF1 + Newton) で量論 H2-air の着火遅れを Cantera 参照 CSV と比較
//   (3) 絶対 datum と sensible datum (thermoHrefTemp) で結果が一致することを確認 (反応熱 Q̇ 経路の検証)
//
// ビルド: g++ -O2 -std=c++17 -I solver_density_cuda solver_density_cuda/tools/test_chemistry.cpp -lyaml-cpp -o /tmp/tchem
// 実行:   /tmp/tchem <mechanism.yaml> <species_db.yaml> [ref.csv] [T0=1200] [p_atm=1] [t_end=1e-3]
//   ref.csv は tools/chem_reference_cantera.py が出力する (time,T,Y_...)。
// =============================================================================
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <yaml-cpp/yaml.h>
#include "../cuda_forge/chemistry_d.cuh"
#include "../cuda_forge/chemistry_mech_io.hpp"

static std::vector<SpeciesThermo> loadSpeciesDB(const std::string& path, const std::vector<std::string>& names)
{
    YAML::Node root = YAML::LoadFile(path);
    std::vector<SpeciesThermo> out;
    for (const auto& nm : names) {
        YAML::Node s = root[nm];
        if (!s) { fprintf(stderr, "species %s not in %s\n", nm.c_str(), path.c_str()); exit(1); }
        SpeciesThermo sp{};
        sp.MW = s["MW"].as<double>(); sp.sigma_LJ = s["LJ_sigma"].as<double>(); sp.eps_kB = s["LJ_eps_kB"].as<double>();
        sp.Tlo = s["Tlo"].as<double>(); sp.Tmid = s["Tmid"].as<double>(); sp.Thi = s["Thi"].as<double>();
        for (int i = 0; i < 9; ++i) { sp.low[i] = s["nasa9_low"][i].as<double>(); sp.high[i] = s["nasa9_high"][i].as<double>(); }
        sp.h_datum = 0.0;
        out.push_back(sp);
    }
    return out;
}

static void applyDatum(std::vector<SpeciesThermo>& sp, double Tref)
{
    for (auto& s : sp) {
        const double h_ref = thermo_h_molar(s, Tref);
        s.low[7] += -h_ref / THERMO_RU; s.high[7] += -h_ref / THERMO_RU; s.h_datum = h_ref;
    }
}

// 小さな密行列 LU (部分ピボット) で A x = b を解く
static void solveDense(int n, std::vector<double> A, std::vector<double>& b)
{
    std::vector<int> piv(n);
    for (int k = 0; k < n; ++k) {
        int p = k; for (int i = k + 1; i < n; ++i) if (fabs(A[i*n+k]) > fabs(A[p*n+k])) p = i;
        if (p != k) { for (int j = 0; j < n; ++j) std::swap(A[k*n+j], A[p*n+j]); std::swap(b[k], b[p]); }
        const double d = A[k*n+k];
        for (int i = k + 1; i < n; ++i) {
            const double f = A[i*n+k] / d;
            if (f == 0.0) continue;
            for (int j = k; j < n; ++j) A[i*n+j] -= f * A[k*n+j];
            b[i] -= f * b[k];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        double s = b[i]; for (int j = i + 1; j < n; ++j) s -= A[i*n+j] * b[j];
        b[i] = s / A[i*n+i];
    }
}

struct Reactor {
    const SpeciesThermo* sp; const ReactionTable* rt; int n; double rho; double e0;   // e0 は sensible datum では Q̇ で更新される
    double T_of(const std::vector<double>& U, double Tg) const {
        std::vector<double> Y(n); for (int s = 0; s < n; ++s) Y[s] = std::max(U[s], 0.0) / rho;
        return thermo_T_from_e(sp, n, Y.data(), e0, Tg, 200.0, 6000.0);
    }
    // ω(U) と (I − dt J_total) を返す
    void rhs(const std::vector<double>& U, double Tg, std::vector<double>& om, double& T, std::vector<double>* Jtot, double dt) const {
        std::vector<double> Y(n); for (int s = 0; s < n; ++s) Y[s] = std::max(U[s], 0.0) / rho;
        T = thermo_T_from_e(sp, n, Y.data(), e0, Tg, 200.0, 6000.0);
        std::vector<double> J(n*n), dOdT(n); double Q;
        chem_source(sp, rt, rho, T, Y.data(), om.data(), &Q, Jtot ? 2 : 0, J.data(), Jtot ? dOdT.data() : nullptr);
        if (Jtot) {
            const double cv = thermo_cp_mix(sp, n, Y.data(), T) - thermo_R_mix(sp, n, Y.data());
            for (int s = 0; s < n; ++s) for (int k = 0; k < n; ++k) {
                const double e_k = thermo_h_mass(sp[k], T) - thermo_R_species(sp[k]) * T;
                const double dTdU = -e_k / (rho * cv);
                (*Jtot)[s*n+k] = ((s == k) ? 1.0 : 0.0) - dt * (J[s*n+k] + dOdT[s] * dTdU);
            }
        }
    }
};

int main(int argc, char** argv)
{
    if (argc < 3) { fprintf(stderr, "usage: %s mech.yaml species_db.yaml [ref.csv] [T0] [p_atm] [t_end]\n", argv[0]); return 1; }
    const std::string mech = argv[1], dbpath = argv[2];
    const std::string ref = (argc > 3) ? argv[3] : "";
    const double T0 = (argc > 4) ? atof(argv[4]) : 1200.0;
    const double p0 = ((argc > 5) ? atof(argv[5]) : 1.0) * 101325.0;
    const double tEnd = (argc > 6) ? atof(argv[6]) : 1.0e-3;
    // 刻み制御 (環境変数で厳しくできる: TCHEM_RELMAX, TCHEM_DTMAX, TCHEM_DTCAP)
    const double RELMAX = getenv("TCHEM_RELMAX") ? atof(getenv("TCHEM_RELMAX")) : 0.3;
    const double DTMAX  = getenv("TCHEM_DTMAX")  ? atof(getenv("TCHEM_DTMAX"))  : 40.0;
    const double DTCAP  = getenv("TCHEM_DTCAP")  ? atof(getenv("TCHEM_DTCAP"))  : 2e-6;

    // 機構ファイルの phases[0].species を流れ側の順序として使う
    YAML::Node root = YAML::LoadFile(mech);
    std::vector<std::string> names; for (const auto& s : root["phases"][0]["species"]) names.push_back(s.as<std::string>());
    const int n = (int)names.size();
    std::vector<SpeciesThermo> spAbs = loadSpeciesDB(dbpath, names);
    ReactionTable rt; chem_io::loadMechanism(mech, names, rt);
    printf("mechanism %s: %d reactions, %d species\n", mech.c_str(), rt.nReac, n);

    // 量論 H2-air: H2:O2:N2 = 2:1:3.76 (mol)
    std::vector<double> X0(n, 0.0);
    auto idx = [&](const char* nm){ return chem_io::findSpecies(names, nm); };
    X0[idx("H2")] = 2.0; X0[idx("O2")] = 1.0; if (idx("N2") >= 0) X0[idx("N2")] = 3.76;
    double xs = 0; for (double x : X0) xs += x; for (double& x : X0) x /= xs;
    double Wmix = 0; for (int s = 0; s < n; ++s) Wmix += X0[s] * spAbs[s].MW;
    std::vector<double> Y0(n); for (int s = 0; s < n; ++s) Y0[s] = X0[s] * spAbs[s].MW / Wmix;
    const double rho = p0 * Wmix / (THERMO_RU * T0);

    // ---- (1) Jacobian FD 照合 (途中状態: 少しラジカルを混ぜる) ----
    {
        std::vector<double> Y = Y0; for (double& y : Y) y += 2e-4;   // 全種を正にして片側差分 (クランプ) を避ける
        Y[idx("H")] += 1e-4; Y[idx("OH")] += 1e-3; Y[idx("O")] += 1e-4;
        double ys = 0; for (double y : Y) ys += y; for (double& y : Y) y /= ys;
        const double T = 1500.0;
        std::vector<double> om(n), J(n*n), dOdT(n); double Q;
        chem_source(spAbs.data(), &rt, rho, T, Y.data(), om.data(), &Q, 2, J.data(), dOdT.data());
        double maxRel = 0, maxRelT = 0;
        for (int k = 0; k < n; ++k) {
            const double h = 1e-5 * rho * Y[k];
            std::vector<double> Yp = Y, Ym = Y; Yp[k] += h / rho; Ym[k] -= h / rho;
            std::vector<double> op(n), omn(n), Jd(n*n); double Qp, Qm;
            chem_source(spAbs.data(), &rt, rho, T, Yp.data(), op.data(), &Qp, 0, Jd.data(), nullptr);
            chem_source(spAbs.data(), &rt, rho, T, Ym.data(), omn.data(), &Qm, 0, Jd.data(), nullptr);
            for (int s = 0; s < n; ++s) {
                const double fd = (op[s] - omn[s]) / (2 * h);
                const double scale = std::max(fabs(J[s*n+k]), 1e-3 * fabs(om[s]) / std::max(rho*Y[k], 1e-12));
                if (scale > 1e-30) maxRel = std::max(maxRel, fabs(fd - J[s*n+k]) / scale);
            }
        }
        {
            const double h = 1e-3;
            std::vector<double> op(n), omn(n), Jd(n*n); double Qp, Qm;
            chem_source(spAbs.data(), &rt, rho, T + h, Y.data(), op.data(), &Qp, 0, Jd.data(), nullptr);
            chem_source(spAbs.data(), &rt, rho, T - h, Y.data(), omn.data(), &Qm, 0, Jd.data(), nullptr);
            for (int s = 0; s < n; ++s) {
                const double fd = (op[s] - omn[s]) / (2 * h);
                if (fabs(dOdT[s]) > 1e-30) maxRelT = std::max(maxRelT, fabs(fd - dOdT[s]) / fabs(dOdT[s]));
            }
        }
        printf("[jacobian] max rel diff vs central FD: dω/d(ρY) %.2e, dω/dT %.2e  (%s)\n", maxRel, maxRelT,
               (maxRel < 1e-5 && maxRelT < 1e-5) ? "PASS" : "FAIL");
        double sum = 0; for (int s = 0; s < n; ++s) sum += om[s];
        printf("[mass]     Σω_s = %.3e kg/m3/s (|ω|max %.3e)\n", sum, *std::max_element(om.begin(), om.end(), [](double a, double b){ return fabs(a) < fabs(b); }));
    }

    // ---- (2)(3) 0-D 定積反応器: 絶対 datum と sensible datum で 2 回 ----
    std::vector<std::vector<std::pair<double,double>>> hist(2);   // (t, T)
    for (int pass = 0; pass < 2; ++pass) {
        std::vector<SpeciesThermo> sp = spAbs;
        if (pass == 1) applyDatum(sp, 298.15);
        Reactor R{sp.data(), &rt, n, rho, thermo_e_mix(sp.data(), n, Y0.data(), T0)};
        std::vector<double> U(n); for (int s = 0; s < n; ++s) U[s] = rho * Y0[s];
        double t = 0, dt = 1e-9, T = T0;
        std::vector<double> om(n), Jt(n*n), rhsv(n);
        int nStep = 0;
        hist[pass].push_back({0.0, T0});
        while (t < tEnd) {
            if (t + dt > tEnd) dt = tEnd - t;
            // BDF1: U^{n+1} = U^n + dt ω(U^{n+1}), Newton
            std::vector<double> Un = U; double Tn = T; bool ok = false; int it;
            for (it = 0; it < 12; ++it) {
                double Tk; R.rhs(U, Tn, om, Tk, &Jt, dt);
                double res = 0;
                for (int s = 0; s < n; ++s) { rhsv[s] = -(U[s] - Un[s] - dt * om[s]); res = std::max(res, fabs(rhsv[s]) / (rho * 1e-10 + fabs(U[s]) * 1e-8 + 1e-30)); }
                if (res < 1.0) { ok = true; T = Tk; break; }
                solveDense(n, Jt, rhsv);
                for (int s = 0; s < n; ++s) U[s] = std::max(U[s] + rhsv[s], 0.0);
                Tn = Tk;
            }
            if (!ok) { U = Un; T = Tn; dt *= 0.5; continue; }
            // 受理判定: 主要種の相対変化 <30%, |ΔT|<40K。超えたら刻みを半分にしてやり直し (BDF1 の誘導期増幅を防ぐ)
            double relmax = 0; for (int s = 0; s < n; ++s) if (Un[s] > 1e-7 * rho) relmax = std::max(relmax, fabs(U[s] - Un[s]) / Un[s]);
            const double dT = fabs(T - Tn);
            if ((relmax > RELMAX || dT > DTMAX) && dt > 1e-12) { U = Un; T = Tn; dt *= 0.5; continue; }
            // sensible datum: 保存エネルギーは反応熱 Δe = −Σ_s (h_datum_s/W_s) Δ(ρY_s)/ρ で更新 (CFD の res_roe += VQ̇ と同じ離散恒等式)
            double de = 0; for (int s = 0; s < n; ++s) de -= (sp[s].h_datum / sp[s].MW) * (U[s] - Un[s]) / rho;
            R.e0 += de;
            // 質量分率の再正規化
            double us = 0; for (double u : U) us += u; for (double& u : U) u *= rho / us;
            T = R.T_of(U, T);
            t += dt; ++nStep;
            hist[pass].push_back({t, T});
            if (it <= 2 && relmax < 0.3*RELMAX && dT < 0.25*DTMAX) dt *= 1.5;
            dt = std::min(dt, DTCAP);
        }
        // 着火遅れ = max dT/dt
        double tauIgn = 0, dTmax = 0;
        for (size_t i = 1; i < hist[pass].size(); ++i) {
            const double d = (hist[pass][i].second - hist[pass][i-1].second) / (hist[pass][i].first - hist[pass][i-1].first);
            if (d > dTmax) { dTmax = d; tauIgn = hist[pass][i].first; }
        }
        printf("[reactor%s] T0=%.0fK p=%.2fatm rho=%.4f: steps=%d  T_end=%.2fK  tau_ign=%.3f us\n",
               pass == 0 ? " abs-datum" : " sens-datum", T0, p0/101325.0, rho, nStep, T, tauIgn * 1e6);
        if (pass == 1) {
            // 同一時刻 T の最大差 (線形補間)
            double maxd = 0;
            for (const auto& pt : hist[1]) {
                const auto& h0 = hist[0];
                auto it2 = std::lower_bound(h0.begin(), h0.end(), pt, [](const std::pair<double,double>& a, const std::pair<double,double>& b){ return a.first < b.first; });
                if (it2 == h0.begin() || it2 == h0.end()) continue;
                const auto a = *(it2 - 1), b = *it2;
                const double Ti = a.second + (b.second - a.second) * (pt.first - a.first) / (b.first - a.first);
                maxd = std::max(maxd, fabs(Ti - pt.second));
            }
            printf("[datum]    max |T_abs − T_sens| at same t = %.1f K (刻み列の差で着火前線がずれる分を含む)\n", maxd);
        }
        if (pass == 0) {
            std::ofstream f("/tmp/tchem_forge.csv"); f << "time,T\n";
            for (const auto& pt : hist[0]) f << pt.first << "," << pt.second << "\n";
        }
    }

    // ---- Cantera 参照との比較 ----
    if (!ref.empty()) {
        std::ifstream f(ref); std::string line; std::getline(f, line);
        std::vector<std::pair<double,double>> rh;
        while (std::getline(f, line)) { std::stringstream ss(line); std::string a, b; std::getline(ss, a, ','); std::getline(ss, b, ','); rh.push_back({atof(a.c_str()), atof(b.c_str())}); }
        double tauRef = 0, dTmax = 0;
        for (size_t i = 1; i < rh.size(); ++i) { const double d = (rh[i].second - rh[i-1].second) / (rh[i].first - rh[i-1].first); if (d > dTmax) { dTmax = d; tauRef = rh[i].first; } }
        double maxd = 0;
        for (const auto& pt : rh) {
            const auto& h0 = hist[0];
            auto it2 = std::lower_bound(h0.begin(), h0.end(), pt, [](const std::pair<double,double>& a, const std::pair<double,double>& b){ return a.first < b.first; });
            if (it2 == h0.begin() || it2 == h0.end()) continue;
            const auto a = *(it2 - 1), b = *it2;
            const double Ti = a.second + (b.second - a.second) * (pt.first - a.first) / (b.first - a.first);
            maxd = std::max(maxd, fabs(Ti - pt.second));
        }
        printf("[cantera]  tau_ign ref=%.3f us, T_end ref=%.2fK, max|ΔT| at same t=%.2fK\n", tauRef * 1e6, rh.back().second, maxd);
    }
    return 0;
}
