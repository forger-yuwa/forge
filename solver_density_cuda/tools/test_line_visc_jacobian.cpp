// U-J (plan time_integration-line-viscous-jacobian §6): 薄層の粘性・熱伝導の Jacobian
// (block_dplur::accumulate_thinlayer_visc_jacobian、`lineViscCoupling: 2`) の単体照合。
//   (1) D = −∂R_i/∂Q_i・K = ∂R_i/∂Q_j を、同じ薄層の流束モデル R の中心差分と照合 (列ごとの誤差を行列の最大値で正規化、≤ 1e-6)
//   (2) 零空間: δQ = δρ(1, u, v, w, E) に対する D・K の作用 (≤ 1e-12、相対)
//   (3) 等温壁の拘束 Δ(ρE)_w = e_w Δρ_w (u_w = 0) の下で ΔT_w = 0、壁のフラグで K の該当行が 0
//   (4) 8 節点のラインを組み、block-Thomas (lineThomas_d と同じ前進消去・後退代入の host の写し) と密行列の解の差 (≤ 1e-10、相対)
//   (5) float 版と double 版の差 (≤ 1e-4、相対)
// ビルド: g++ -O2 -std=c++17 -I solver_density_cuda solver_density_cuda/tools/test_line_visc_jacobian.cpp -o /tmp/tlvj && /tmp/tlvj
#include <cstdio>
#include <cmath>
#include <random>
#include <vector>
#include <array>
#include <algorithm>
#include "../cuda_forge/block_dplur_jacobian_d.cuh"

struct Face { double beta, kappa, n[3], fi; };
struct Node { double rho, u[3], rhoE, gam, cp; };

static Node from_q(const double q[5], double gam, double cp) {
    Node s; s.rho = q[0]; for (int a = 0; a < 3; ++a) s.u[a] = q[1 + a] / q[0]; s.rhoE = q[4]; s.gam = gam; s.cp = cp; return s;
}
static void to_q(const Node& s, double q[5]) { q[0] = s.rho; for (int a = 0; a < 3; ++a) q[1 + a] = s.rho * s.u[a]; q[4] = s.rhoE; }
// 物性凍結の温度 T = (γ/c_p) e (定数の基準のずれは微分に効かない)
static double temp(const Node& s) { const double q2 = s.u[0]*s.u[0] + s.u[1]*s.u[1] + s.u[2]*s.u[2]; return (s.gam / s.cp) * (s.rhoE / s.rho - 0.5 * q2); }

// 薄層の流束モデル (自節点 i の残差への寄与): [0, βPΔu, κΔT + (βPΔu)·ū] (別実装)
static void flux(const Face& f, const Node& i, const Node& j, double R[5]) {
    double P[3][3];
    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) P[a][b] = (a == b ? 1.0 : 0.0) + f.n[a] * f.n[b] / 3.0;
    double tau[3], ub[3];
    for (int a = 0; a < 3; ++a) { tau[a] = 0; for (int b = 0; b < 3; ++b) tau[a] += f.beta * P[a][b] * (j.u[b] - i.u[b]); ub[a] = f.fi * i.u[a] + (1 - f.fi) * j.u[a]; }
    R[0] = 0; for (int a = 0; a < 3; ++a) R[1 + a] = tau[a];
    R[4] = f.kappa * (temp(j) - temp(i)) + tau[0] * ub[0] + tau[1] * ub[1] + tau[2] * ub[2];
}

template<typename T>
static void jac(const Face& f, const Node& i, const Node& j, bool vf, bool tf, double D[5][5], double K[5][5]) {
    T Dt[5][5] = {}, Kt[5][5] = {};
    block_dplur::accumulate_thinlayer_visc_jacobian<T>((T)f.beta, (T)f.kappa, (T)f.n[0], (T)f.n[1], (T)f.n[2], (T)f.fi,
        (T)i.rho, (T)i.u[0], (T)i.u[1], (T)i.u[2], (T)i.rhoE, (T)i.gam, (T)i.cp,
        (T)j.rho, (T)j.u[0], (T)j.u[1], (T)j.u[2], (T)j.rhoE, (T)j.gam, (T)j.cp, vf, tf, Dt, Kt);
    for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c) { D[r][c] = Dt[r][c]; K[r][c] = Kt[r][c]; }
}

static double maxabs(const double A[5][5]) { double m = 0; for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c) m = std::max(m, std::fabs(A[r][c])); return m; }

int main() {
    std::mt19937_64 g(20261009);
    auto U = [&](double a, double b) { return std::uniform_real_distribution<double>(a, b)(g); };
    auto rnd_node = [&](bool wall) {
        Node s; s.rho = U(0.01, 2.0); s.gam = U(1.1, 1.4); s.cp = U(1000, 2500);
        for (int a = 0; a < 3; ++a) s.u[a] = wall ? 0.0 : U(-1700, 1700);
        const double T = U(250, 2000), e = T * s.cp / s.gam, q2 = s.u[0]*s.u[0] + s.u[1]*s.u[1] + s.u[2]*s.u[2];
        s.rhoE = s.rho * (e + 0.5 * q2);
        return s;
    };
    auto rnd_face = [&]() {
        Face f; f.beta = std::pow(10.0, U(-6, -1)); f.kappa = std::pow(10.0, U(-4, 1));
        double n[3] = {U(-1, 1), U(-1, 1), U(-1, 1)}; const double L = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        for (int a = 0; a < 3; ++a) f.n[a] = n[a] / L;
        const double f0 = U(0.2, 0.8); f.fi = (g() & 1) ? f0 : 1.0 - f0; return f;
    };
    double e_fd = 0, e_null = 0, e_wallT = 0, e_flag = 0, e_f32 = 0;
    for (int t = 0; t < 200; ++t) {
        const Face f = rnd_face(); const Node i = rnd_node(false), j = rnd_node(false);
        double D[5][5], K[5][5]; jac<double>(f, i, j, false, false, D, K);
        double qi[5], qj[5]; to_q(i, qi); to_q(j, qj);
        const double sD = std::max(maxabs(D), 1e-300), sK = std::max(maxabs(K), 1e-300);
        for (int c = 0; c < 5; ++c) {
            const double h = 1e-6 * std::max(std::fabs(qi[c]), 1e-3), hj = 1e-6 * std::max(std::fabs(qj[c]), 1e-3);
            double qp[5], qm[5], Rp[5], Rm[5];
            std::copy(qi, qi + 5, qp); std::copy(qi, qi + 5, qm); qp[c] += h; qm[c] -= h;
            flux(f, from_q(qp, i.gam, i.cp), j, Rp); flux(f, from_q(qm, i.gam, i.cp), j, Rm);
            for (int r = 0; r < 5; ++r) e_fd = std::max(e_fd, std::fabs(-(Rp[r] - Rm[r]) / (2 * h) - D[r][c]) / sD);
            std::copy(qj, qj + 5, qp); std::copy(qj, qj + 5, qm); qp[c] += hj; qm[c] -= hj;
            flux(f, i, from_q(qp, j.gam, j.cp), Rp); flux(f, i, from_q(qm, j.gam, j.cp), Rm);
            for (int r = 0; r < 5; ++r) e_fd = std::max(e_fd, std::fabs((Rp[r] - Rm[r]) / (2 * hj) - K[r][c]) / sK);
        }
        // 零空間
        const double dqi[5] = {1, i.u[0], i.u[1], i.u[2], i.rhoE / i.rho}, dqj[5] = {1, j.u[0], j.u[1], j.u[2], j.rhoE / j.rho};
        for (int r = 0; r < 5; ++r) {
            double a = 0, b = 0, na = 0, nb = 0;
            for (int c = 0; c < 5; ++c) { a += D[r][c] * dqi[c]; b += K[r][c] * dqj[c]; na += std::fabs(D[r][c] * dqi[c]); nb += std::fabs(K[r][c] * dqj[c]); }
            e_null = std::max({e_null, std::fabs(a) / std::max(na, 1e-300), std::fabs(b) / std::max(nb, 1e-300)});
        }
        // 等温壁の隣: u_w = 0、拘束 δQ_w = δρ(1, 0, 0, 0, e_w) で ΔT_w = 0 → フラグ無しの K の作用も 0、フラグ付きは該当行が 0
        const Node w = rnd_node(true);
        double Dw[5][5], Kw[5][5], Kf[5][5]; jac<double>(f, i, w, false, false, Dw, Kw); jac<double>(f, i, w, true, true, Dw, Kf);
        const double dqw[5] = {1, 0, 0, 0, w.rhoE / w.rho};
        for (int r = 0; r < 5; ++r) {
            double a = 0, na = 0; for (int c = 0; c < 5; ++c) { a += Kw[r][c] * dqw[c]; na += std::fabs(Kw[r][c] * dqw[c]); }
            e_wallT = std::max(e_wallT, std::fabs(a) / std::max(na, 1e-300));
            for (int c = 0; c < 5; ++c) e_flag = std::max(e_flag, std::fabs(Kf[r][c]));
        }
        // float と double
        double Df[5][5], Kff[5][5]; jac<float>(f, i, j, false, false, Df, Kff);
        for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c)
            e_f32 = std::max({e_f32, std::fabs(Df[r][c] - D[r][c]) / sD, std::fabs(Kff[r][c] - K[r][c]) / sK});
    }
    // (4) 8 節点のライン: D_k = (V/Δτ) I + (両面の薄層の D)、−K_prev/−K_next。Thomas (host の写し) と密行列の解
    const int N = 8; std::vector<Node> s(N); for (auto& x : s) x = rnd_node(false);
    std::vector<Face> fc(N - 1); for (auto& x : fc) x = rnd_face();
    std::vector<std::array<std::array<double, 5>, 5>> Dk(N), Kp(N), Kn(N);
    for (int k = 0; k < N; ++k) for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c) { Dk[k][r][c] = (r == c) ? std::pow(10.0, U(-2, 1)) : 0.0; Kp[k][r][c] = Kn[k][r][c] = 0; }
    for (int k = 0; k < N; ++k) {
        for (int side = 0; side < 2; ++side) {
            const int jn = side ? k + 1 : k - 1; if (jn < 0 || jn >= N) continue;
            Face f = fc[side ? k : k - 1]; if (side == 0) f.fi = 1.0 - f.fi;   // 面の向きで自節点の重みが入れ替わる
            double D[5][5], K[5][5]; jac<double>(f, s[k], s[jn], false, false, D, K);
            for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c) { Dk[k][r][c] += D[r][c]; (side ? Kn : Kp)[k][r][c] += K[r][c]; }
        }
    }
    std::vector<double> rhs(5 * N); for (auto& x : rhs) x = U(-1, 1);
    // 密行列: 行 k は D_k x_k − Kp_k x_{k−1} − Kn_k x_{k+1} = rhs_k
    const int M = 5 * N; std::vector<double> A(M * M, 0.0), b = rhs;
    for (int k = 0; k < N; ++k) for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c) {
        A[(5 * k + r) * M + 5 * k + c] = Dk[k][r][c];
        if (k > 0) A[(5 * k + r) * M + 5 * (k - 1) + c] = -Kp[k][r][c];
        if (k < N - 1) A[(5 * k + r) * M + 5 * (k + 1) + c] = -Kn[k][r][c];
    }
    for (int p = 0; p < M; ++p) {                                          // 部分ピボットの Gauss 消去
        int pv = p; for (int r = p + 1; r < M; ++r) if (std::fabs(A[r * M + p]) > std::fabs(A[pv * M + p])) pv = r;
        if (pv != p) { for (int c = 0; c < M; ++c) std::swap(A[p * M + c], A[pv * M + c]); std::swap(b[p], b[pv]); }
        for (int r = p + 1; r < M; ++r) { const double m = A[r * M + p] / A[p * M + p]; for (int c = p; c < M; ++c) A[r * M + c] -= m * A[p * M + c]; b[r] -= m * b[p]; }
    }
    std::vector<double> xd(M); for (int r = M - 1; r >= 0; --r) { double acc = b[r]; for (int c = r + 1; c < M; ++c) acc -= A[r * M + c] * xd[c]; xd[r] = acc / A[r * M + r]; }
    // Thomas (lineThomas_d と同じ: M̃_k = D_k − Kp_k W_{k−1}、b̃_k = rhs_k + Kp_k y_{k−1}、y_k = M̃⁻¹ b̃_k、W_k = M̃⁻¹ Kn_k、x_k = y_k + W_k x_{k+1})
    auto solve5 = [](std::array<std::array<double, 5>, 5> Mm, double x[5]) {
        for (int p = 0; p < 5; ++p) { int pv = p; for (int r = p + 1; r < 5; ++r) if (std::fabs(Mm[r][p]) > std::fabs(Mm[pv][p])) pv = r;
            std::swap(Mm[p], Mm[pv]); std::swap(x[p], x[pv]);
            for (int r = p + 1; r < 5; ++r) { const double m = Mm[r][p] / Mm[p][p]; for (int c = p; c < 5; ++c) Mm[r][c] -= m * Mm[p][c]; x[r] -= m * x[p]; } }
        for (int r = 4; r >= 0; --r) { double acc = x[r]; for (int c = r + 1; c < 5; ++c) acc -= Mm[r][c] * x[c]; x[r] = acc / Mm[r][r]; }
    };
    std::vector<std::array<std::array<double, 5>, 5>> W(N); std::vector<std::array<double, 5>> y(N);
    for (int k = 0; k < N; ++k) {
        auto Mk = Dk[k]; double bk[5]; for (int r = 0; r < 5; ++r) bk[r] = rhs[5 * k + r];
        if (k > 0) for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c) { double acc = 0; for (int m = 0; m < 5; ++m) acc += Kp[k][r][m] * W[k - 1][m][c]; Mk[r][c] -= acc; }
        if (k > 0) for (int r = 0; r < 5; ++r) for (int m = 0; m < 5; ++m) bk[r] += Kp[k][r][m] * y[k - 1][m];
        solve5(Mk, bk); for (int r = 0; r < 5; ++r) y[k][r] = bk[r];
        for (int c = 0; c < 5; ++c) { double col[5]; for (int r = 0; r < 5; ++r) col[r] = Kn[k][r][c]; solve5(Mk, col); for (int r = 0; r < 5; ++r) W[k][r][c] = col[r]; }
    }
    std::vector<double> xt(M); for (int r = 0; r < 5; ++r) xt[5 * (N - 1) + r] = y[N - 1][r];
    for (int k = N - 2; k >= 0; --k) for (int r = 0; r < 5; ++r) { double acc = y[k][r]; for (int m = 0; m < 5; ++m) acc += W[k][r][m] * xt[5 * (k + 1) + m]; xt[5 * k + r] = acc; }
    double nx = 0, nd = 0; for (int r = 0; r < M; ++r) { nd = std::max(nd, std::fabs(xt[r] - xd[r])); nx = std::max(nx, std::fabs(xd[r])); }
    const double e_line = nd / std::max(nx, 1e-300);

    const bool ok = e_fd <= 1e-6 && e_null <= 1e-12 && e_wallT <= 1e-12 && e_flag == 0.0 && e_line <= 1e-10 && e_f32 <= 1e-4;
    std::printf("(1) 中心差分との差 (行列の最大で正規化) %.3e (≤ 1e-6)\n", e_fd);
    std::printf("(2) 零空間 δρ(1,u,v,w,E) の作用 %.3e (≤ 1e-12)\n", e_null);
    std::printf("(3) 等温壁の拘束で ΔT_w の作用 %.3e (≤ 1e-12)、壁のフラグで K の該当行の最大 %.3e (= 0)\n", e_wallT, e_flag);
    std::printf("(4) 8 節点のライン: Thomas と密行列の解の差 %.3e (≤ 1e-10)\n", e_line);
    std::printf("(5) float と double の差 %.3e (≤ 1e-4)\n", e_f32);
    std::printf("VERDICT: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
