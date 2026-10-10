// U-J (plan time_integration-line-viscous-jacobian §6、拡張は plan time_integration-line-viscous-jacobian-faceh §6.3): 薄層の粘性・熱伝導の Jacobian
// (block_dplur::accumulate_thinlayer_visc_jacobian、`lineViscCoupling: 2`) の単体照合。A = 共通関数の D・K、B = 試験の中の独立な流束 `flux` の中心差分。
//   (1)  D = −∂R_i/∂Q_i・K = ∂R_i/∂Q_j と中心差分の差を行列の最大値で正規化 (記録だけ。次元の違う小さい列の誤差を保証しない)
//   (1′) 同じ照合を列ごとに分類する。列 c の大きさ s_c = max_r |J_rc| で割った誤差 e(h)・e(h/2) と再現 |FD_h − FD_{h/2}| / s_c:
//        再現 > 1e-7 → 解像しない、両方 ≤ 1e-6 → 合格、両方 > 1e-6 → 不一致 (FAIL)、片方だけ → 閾値をまたぐ (判別不能)。
//        零列 (s_c = 0) は事前に固定した尺度で無次元化した |FD_rc|·q_ref,c / R_ref,r に同じ規則を許容 1e-10・再現 1e-11 で使う。
//        J に非有限 → FAIL、差分に非有限 → 判別不能。
//        標本: 乱数の 200 組 (rand)、端の場合 (edge: Δu = 0、両側 |u| 1700 m/s、ρ 0.01 と 2、軸に沿う法線、f_i = 0.5)、零列を持つ標本 (zero: β = 0・静止)。
//        edge・zero は全列の解像を要求し、rand は解像しない・またぐ列を rand の列数の 1 % まで許す (標本・側・列・誤差を全部書き出す)。
//   (1″) 面の向きの整合: (j, i, 1 − f_i) の D・K が (i, j, f_i) の K・D と一致 (列ごとに s_c で割って ≤ 1e-12)
//   (2)  零空間: δQ = δρ(1, u, v, w, E) に対する D・K の作用 (≤ 1e-12、相対)
//   (3)  等温壁の拘束 Δ(ρE)_w = e_w Δρ_w (u_w = 0) の下で ΔT_w = 0、壁のフラグで K の該当行が 0
//   (3′) 片側だけ固定の隣を、自由度を消去した系で列ごとに照合 (規則は (1′) と同じ、全列の解像を要求):
//        速度だけ固定 (u_j = 0): 自由な列 (ρ, ρE) の K(jVelFixed) と差分。温度だけ固定: 自由な列 (ρ, ρu, ρv, ρw) を δT_j = 0 の写像 B で写した K(jTempFixed)·B と方向微分
//   (4)  8 節点のラインを組み、block-Thomas (lineThomas_d と同じ前進消去・後退代入の host の写し) と密行列の解の差 (≤ 1e-10、相対)
//   (5)  float 版と double 版の差 (≤ 1e-4、行列の最大で正規化)。列ごとの値は記録だけ
//   判定に使う量に非有限が 1 件でもあれば PASS にしない。FAIL は失敗の項目を分けて書く (微分の不一致だけが微分の実装への反証)。
// ビルド: g++ -O2 -std=c++17 -Wno-unknown-pragmas -I solver_density_cuda solver_density_cuda/tools/test_line_visc_jacobian.cpp -o /tmp/tlvj && /tmp/tlvj
#include <cstdio>
#include <cmath>
#include <random>
#include <vector>
#include <array>
#include <string>
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
static double speed(const Node& s) { return std::sqrt(s.u[0]*s.u[0] + s.u[1]*s.u[1] + s.u[2]*s.u[2]); }

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

// 判定に使う量の非有限の件数 (1 件でも PASS にしない)。std::max は NaN を落とすので最大の更新はこれで行う
static long g_nonfinite = 0;
static void upd(double& m, double v) { if (!std::isfinite(v)) { ++g_nonfinite; return; } if (v > m) m = v; }
static double maxabs(const double A[5][5]) { double m = 0; for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c) upd(m, std::fabs(A[r][c])); return m; }

// 中心差分の Jacobian: side 0 = −∂R_i/∂Q_i (D と同じ符号)、side 1 = ∂R_i/∂Q_j (K と同じ符号)。差分幅は hs·max(|q_c|, 1e-3)
static void fd_jac(const Face& f, const Node& i, const Node& j, int side, double hs, double J[5][5]) {
    double q[5]; to_q(side ? j : i, q);
    for (int c = 0; c < 5; ++c) {
        const double h = hs * std::max(std::fabs(q[c]), 1e-3);
        double qp[5], qm[5], Rp[5], Rm[5];
        std::copy(q, q + 5, qp); std::copy(q, q + 5, qm); qp[c] += h; qm[c] -= h;
        if (side) { flux(f, i, from_q(qp, j.gam, j.cp), Rp); flux(f, i, from_q(qm, j.gam, j.cp), Rm); }
        else      { flux(f, from_q(qp, i.gam, i.cp), j, Rp); flux(f, from_q(qm, i.gam, i.cp), j, Rm); }
        for (int r = 0; r < 5; ++r) J[r][c] = (side ? 1.0 : -1.0) * (Rp[r] - Rm[r]) / (2 * h);
    }
}
// 隣 j の状態の方向 dq への方向微分 ∂R_i/∂Q_j · dq (差分幅 hs)
static void dir_fd(const Face& f, const Node& i, const Node& j, const double dq[5], double hs, double out[5]) {
    double q[5], qp[5], qm[5], Rp[5], Rm[5]; to_q(j, q);
    for (int c = 0; c < 5; ++c) { qp[c] = q[c] + hs * dq[c]; qm[c] = q[c] - hs * dq[c]; }
    flux(f, i, from_q(qp, j.gam, j.cp), Rp); flux(f, i, from_q(qm, j.gam, j.cp), Rm);
    for (int r = 0; r < 5; ++r) out[r] = (Rp[r] - Rm[r]) / (2 * hs);
}

// 列の分類 (§6.3 (1′))
enum ColClass { C_PASS = 0, C_FAIL, C_UNRES, C_STRADDLE, C_NF_J, C_NF_FD, C_NCLS };
static const char* CLS_NAME[C_NCLS] = {"合格", "不一致", "解像しない", "閾値をまたぐ", "J が非有限", "差分が非有限"};
struct ColRec { int sample; std::string set; char side; int col; ColClass cls; double e1, e2, rep; bool zero; };
// Jc・F1c (幅 h)・F2c (幅 h/2) は列ベクトル。zs[r] は零列のときの無次元化の係数 q_ref,c / R_ref,r
static ColClass classify(const double Jc[5], const double F1c[5], const double F2c[5], const double zs[5], double& e1, double& e2, double& rep, bool& zero) {
    e1 = e2 = rep = 0; zero = false;
    for (int r = 0; r < 5; ++r) if (!std::isfinite(Jc[r])) return C_NF_J;
    for (int r = 0; r < 5; ++r) if (!std::isfinite(F1c[r]) || !std::isfinite(F2c[r])) return C_NF_FD;
    double s = 0; for (int r = 0; r < 5; ++r) s = std::max(s, std::fabs(Jc[r]));
    double tol = 1e-6, reptol = 1e-7;
    if (s == 0.0) {
        zero = true; tol = 1e-10; reptol = 1e-11;
        for (int r = 0; r < 5; ++r) { e1 = std::max(e1, std::fabs(F1c[r]) * zs[r]); e2 = std::max(e2, std::fabs(F2c[r]) * zs[r]); rep = std::max(rep, std::fabs(F1c[r] - F2c[r]) * zs[r]); }
    } else {
        for (int r = 0; r < 5; ++r) { e1 = std::max(e1, std::fabs(F1c[r] - Jc[r]) / s); e2 = std::max(e2, std::fabs(F2c[r] - Jc[r]) / s); rep = std::max(rep, std::fabs(F1c[r] - F2c[r]) / s); }
    }
    if (!std::isfinite(e1) || !std::isfinite(e2) || !std::isfinite(rep)) return C_NF_FD;
    if (rep > reptol) return C_UNRES;
    if (e1 <= tol && e2 <= tol) return C_PASS;
    if (e1 > tol && e2 > tol) return C_FAIL;
    return C_STRADDLE;
}
// 零列の無次元化: q_ref は節点の状態 (ρ・ρ max(|u|, 100)・ρE)、R_ref は行ごと (運動量 β U、エネルギー κ T + β U²、U = max(|u_i|, |u_j|, 100)、T = max(T_i, T_j))
static void zero_scale(const Face& f, const Node& i, const Node& j, const Node& col_node, int c, double zs[5]) {
    const double U = std::max({speed(i), speed(j), 100.0}), Tm = std::max(temp(i), temp(j));
    const double qref = (c == 0) ? col_node.rho : (c == 4) ? col_node.rhoE : col_node.rho * std::max(speed(col_node), 100.0);
    const double Rm = std::max(f.beta, 1e-12) * U, Re = f.kappa * Tm + std::max(f.beta, 1e-12) * U * U;
    for (int r = 0; r < 4; ++r) zs[r] = qref / Rm;
    zs[4] = qref / Re;
}

int main() {
    std::mt19937_64 g(20261009);
    auto U = [&](double a, double b) { return std::uniform_real_distribution<double>(a, b)(g); };
    auto rnd_node = [&](bool wall) {
        Node s; s.rho = U(0.01, 2.0); s.gam = U(1.1, 1.4); s.cp = U(1000, 2500);
        for (int a = 0; a < 3; ++a) s.u[a] = wall ? 0.0 : U(-1700, 1700);   // 速度の各成分が ±1700 m/s (ノルムは最大約 2944 m/s)
        const double T = U(250, 2000), e = T * s.cp / s.gam, q2 = s.u[0]*s.u[0] + s.u[1]*s.u[1] + s.u[2]*s.u[2];
        s.rhoE = s.rho * (e + 0.5 * q2);
        return s;
    };
    auto rnd_face = [&]() {   // 法線は立方体の一様乱数を正規化したもの (球面一様ではない)
        Face f; f.beta = std::pow(10.0, U(-6, -1)); f.kappa = std::pow(10.0, U(-4, 1));
        double n[3] = {U(-1, 1), U(-1, 1), U(-1, 1)}; const double L = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        for (int a = 0; a < 3; ++a) f.n[a] = n[a] / L;
        const double f0 = U(0.2, 0.8); f.fi = (g() & 1) ? f0 : 1.0 - f0; return f;
    };
    // (1)・(2)・(3)・(5): 既存の 200 組
    double e_fd = 0, e_null = 0, e_wallT = 0, e_flag = 0, e_f32 = 0, e_f32_col = 0;
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
            for (int r = 0; r < 5; ++r) upd(e_fd, std::fabs(-(Rp[r] - Rm[r]) / (2 * h) - D[r][c]) / sD);
            std::copy(qj, qj + 5, qp); std::copy(qj, qj + 5, qm); qp[c] += hj; qm[c] -= hj;
            flux(f, i, from_q(qp, j.gam, j.cp), Rp); flux(f, i, from_q(qm, j.gam, j.cp), Rm);
            for (int r = 0; r < 5; ++r) upd(e_fd, std::fabs((Rp[r] - Rm[r]) / (2 * hj) - K[r][c]) / sK);
        }
        // 零空間
        const double dqi[5] = {1, i.u[0], i.u[1], i.u[2], i.rhoE / i.rho}, dqj[5] = {1, j.u[0], j.u[1], j.u[2], j.rhoE / j.rho};
        for (int r = 0; r < 5; ++r) {
            double a = 0, b = 0, na = 0, nb = 0;
            for (int c = 0; c < 5; ++c) { a += D[r][c] * dqi[c]; b += K[r][c] * dqj[c]; na += std::fabs(D[r][c] * dqi[c]); nb += std::fabs(K[r][c] * dqj[c]); }
            upd(e_null, std::fabs(a) / std::max(na, 1e-300)); upd(e_null, std::fabs(b) / std::max(nb, 1e-300));
        }
        // 等温壁の隣: u_w = 0、拘束 δQ_w = δρ(1, 0, 0, 0, e_w) で ΔT_w = 0 → フラグ無しの K の作用も 0、フラグ付きは該当行が 0
        const Node w = rnd_node(true);
        double Dw[5][5], Kw[5][5], Kf[5][5]; jac<double>(f, i, w, false, false, Dw, Kw); jac<double>(f, i, w, true, true, Dw, Kf);
        const double dqw[5] = {1, 0, 0, 0, w.rhoE / w.rho};
        for (int r = 0; r < 5; ++r) {
            double a = 0, na = 0; for (int c = 0; c < 5; ++c) { a += Kw[r][c] * dqw[c]; na += std::fabs(Kw[r][c] * dqw[c]); }
            upd(e_wallT, std::fabs(a) / std::max(na, 1e-300));
            for (int c = 0; c < 5; ++c) upd(e_flag, std::fabs(Kf[r][c]));
        }
        // float と double (行列の最大で正規化。列ごとの値は記録だけ)
        double Df[5][5], Kff[5][5]; jac<float>(f, i, j, false, false, Df, Kff);
        for (int r = 0; r < 5; ++r) for (int c = 0; c < 5; ++c) { upd(e_f32, std::fabs(Df[r][c] - D[r][c]) / sD); upd(e_f32, std::fabs(Kff[r][c] - K[r][c]) / sK); }
        for (int c = 0; c < 5; ++c) {
            double cD = 0, cK = 0; for (int r = 0; r < 5; ++r) { cD = std::max(cD, std::fabs(D[r][c])); cK = std::max(cK, std::fabs(K[r][c])); }
            for (int r = 0; r < 5; ++r) { if (cD > 0) upd(e_f32_col, std::fabs(Df[r][c] - D[r][c]) / cD); if (cK > 0) upd(e_f32_col, std::fabs(Kff[r][c] - K[r][c]) / cK); }
        }
    }

    // (1′)・(1″): 乱数の 200 組 (rand)、端の場合 (edge)、零列を持つ標本 (zero)
    struct Trial { Face f; Node i, j; std::string set; };
    std::vector<Trial> trials;
    for (int t = 0; t < 200; ++t) { Trial x{rnd_face(), rnd_node(false), rnd_node(false), "rand"}; trials.push_back(x); }
    auto unit_dir = [&](double mag, double u[3]) { double n[3] = {U(-1, 1), U(-1, 1), U(-1, 1)}; const double L = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]); for (int a = 0; a < 3; ++a) u[a] = mag * n[a] / L; };
    auto set_u = [](Node& s, const double u[3]) {   // 内部エネルギーを保って速度を置き換える
        const double q2o = s.u[0]*s.u[0] + s.u[1]*s.u[1] + s.u[2]*s.u[2], e = s.rhoE / s.rho - 0.5 * q2o;
        for (int a = 0; a < 3; ++a) s.u[a] = u[a];
        const double q2 = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
        s.rhoE = s.rho * (e + 0.5 * q2);
    };
    for (int t = 0; t < 6; ++t) {
        Trial x{rnd_face(), rnd_node(false), rnd_node(false), "edge"}; set_u(x.j, x.i.u); trials.push_back(x);                          // Δu = 0
        Trial y{rnd_face(), rnd_node(false), rnd_node(false), "edge"}; double ua[3], ub[3]; unit_dir(1700, ua); unit_dir(1700, ub);
        set_u(y.i, ua); set_u(y.j, ub); trials.push_back(y);                                                                               // 両側 |u| 1700 m/s
        Trial z{rnd_face(), rnd_node(false), rnd_node(false), "edge"};
        { const double ei = z.i.rhoE / z.i.rho, ej = z.j.rhoE / z.j.rho; z.i.rho = (t & 1) ? 0.01 : 2.0; z.j.rho = (t & 1) ? 2.0 : 0.01; z.i.rhoE = z.i.rho * ei; z.j.rhoE = z.j.rho * ej; }
        trials.push_back(z);                                                                                                               // ρ 0.01 と 2
        Trial w{rnd_face(), rnd_node(false), rnd_node(false), "edge"}; for (int a = 0; a < 3; ++a) w.f.n[a] = (a == t % 3) ? ((t & 1) ? -1.0 : 1.0) : 0.0;
        trials.push_back(w);                                                                                                               // 軸に沿う法線
        Trial v{rnd_face(), rnd_node(false), rnd_node(false), "edge"}; v.f.fi = 0.5; trials.push_back(v);                                  // f_i = 0.5
        Trial s{rnd_face(), rnd_node(true), rnd_node(true), "zero"}; s.f.beta = 0.0; trials.push_back(s);                                  // β = 0・静止 → D・K の列 1〜3 が零
    }
    std::vector<ColRec> recs;
    double e_orient = 0;
    auto add_cols = [&](int id, const std::string& set, char side, const double J[5][5], const double F1[5][5], const double F2[5][5], const Face& f, const Node& i, const Node& j) {
        for (int c = 0; c < 5; ++c) {
            double Jc[5], F1c[5], F2c[5], zs[5];
            for (int r = 0; r < 5; ++r) { Jc[r] = J[r][c]; F1c[r] = F1[r][c]; F2c[r] = F2[r][c]; }
            zero_scale(f, i, j, side == 'D' ? i : j, c, zs);
            ColRec rc{id, set, side, c, C_PASS, 0, 0, 0, false};
            rc.cls = classify(Jc, F1c, F2c, zs, rc.e1, rc.e2, rc.rep, rc.zero);
            recs.push_back(rc);
        }
    };
    for (size_t k = 0; k < trials.size(); ++k) {
        const Trial& x = trials[k];
        double D[5][5], K[5][5], F1[5][5], F2[5][5]; jac<double>(x.f, x.i, x.j, false, false, D, K);
        fd_jac(x.f, x.i, x.j, 0, 1e-6, F1); fd_jac(x.f, x.i, x.j, 0, 5e-7, F2); add_cols((int)k, x.set, 'D', D, F1, F2, x.f, x.i, x.j);
        fd_jac(x.f, x.i, x.j, 1, 1e-6, F1); fd_jac(x.f, x.i, x.j, 1, 5e-7, F2); add_cols((int)k, x.set, 'K', K, F1, F2, x.f, x.i, x.j);
        // (1″) 面の向き: (j, i, 1 − f_i) の D・K は (i, j, f_i) の K・D
        Face fr = x.f; fr.fi = 1.0 - x.f.fi;
        double Dr[5][5], Kr[5][5]; jac<double>(fr, x.j, x.i, false, false, Dr, Kr);
        for (int c = 0; c < 5; ++c) {
            double sK = 0, sD = 0; for (int r = 0; r < 5; ++r) { sK = std::max(sK, std::fabs(K[r][c])); sD = std::max(sD, std::fabs(D[r][c])); }
            for (int r = 0; r < 5; ++r) {
                upd(e_orient, sK > 0 ? std::fabs(Dr[r][c] - K[r][c]) / sK : (Dr[r][c] != 0.0 ? 1.0 : 0.0));
                upd(e_orient, sD > 0 ? std::fabs(Kr[r][c] - D[r][c]) / sD : (Kr[r][c] != 0.0 ? 1.0 : 0.0));
            }
        }
    }
    // (3′) 片側だけ固定の隣を、自由度を消去した系で列ごとに
    for (int t = 0; t < 200; ++t) {
        const Face f = rnd_face(); const Node i = rnd_node(false);
        {   // 速度だけ固定 (u_j = 0): 自由な列 (ρ, ρE)。ρu_j = 0 を保つので u_j = 0 のまま
            const Node j = rnd_node(true);
            double D[5][5], K[5][5], F1[5][5], F2[5][5]; jac<double>(f, i, j, true, false, D, K);
            fd_jac(f, i, j, 1, 1e-6, F1); fd_jac(f, i, j, 1, 5e-7, F2);
            for (int c : {0, 4}) {
                double Jc[5], F1c[5], F2c[5], zs[5];
                for (int r = 0; r < 5; ++r) { Jc[r] = K[r][c]; F1c[r] = F1[r][c]; F2c[r] = F2[r][c]; }
                zero_scale(f, i, j, j, c, zs);
                ColRec rc{t, "vfix", 'K', c, C_PASS, 0, 0, 0, false};
                rc.cls = classify(Jc, F1c, F2c, zs, rc.e1, rc.e2, rc.rep, rc.zero); recs.push_back(rc);
            }
        }
        {   // 温度だけ固定: 自由な列 (ρ, ρu, ρv, ρw)、δ(ρE) = (e_j − ½|u_j|²)δρ + u_j·δ(ρu) (e_j は内部エネルギー) で δT_j = 0
            const Node j = rnd_node(false);
            const double q2 = j.u[0]*j.u[0] + j.u[1]*j.u[1] + j.u[2]*j.u[2], ein = j.rhoE / j.rho - 0.5 * q2;
            double D[5][5], K[5][5]; jac<double>(f, i, j, false, true, D, K);
            for (int b = 0; b < 4; ++b) {
                double dq[5] = {0, 0, 0, 0, 0};
                if (b == 0) { dq[0] = j.rho; dq[4] = (ein - 0.5 * q2) * j.rho; }
                else { const double m = j.rho * std::max(speed(j), 100.0); dq[b] = m; dq[4] = j.u[b - 1] * m; }
                double Jc[5], F1c[5], F2c[5], zs[5];
                for (int r = 0; r < 5; ++r) { Jc[r] = 0; for (int c = 0; c < 5; ++c) Jc[r] += K[r][c] * dq[c]; }
                dir_fd(f, i, j, dq, 1e-6, F1c); dir_fd(f, i, j, dq, 5e-7, F2c);
                // 零列の尺度: dq はすでに状態の大きさなので q_ref = 1
                const double Um = std::max({speed(i), speed(j), 100.0}), Tm = std::max(temp(i), temp(j));
                for (int r = 0; r < 4; ++r) zs[r] = 1.0 / (std::max(f.beta, 1e-12) * Um);
                zs[4] = 1.0 / (f.kappa * Tm + std::max(f.beta, 1e-12) * Um * Um);
                ColRec rc{t, "tfix", 'K', b, C_PASS, 0, 0, 0, false};
                rc.cls = classify(Jc, F1c, F2c, zs, rc.e1, rc.e2, rc.rep, rc.zero); recs.push_back(rc);
            }
        }
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
    double nx = 0, nd = 0; for (int r = 0; r < M; ++r) { upd(nd, std::fabs(xt[r] - xd[r])); upd(nx, std::fabs(xd[r])); }
    const double e_line = nd / std::max(nx, 1e-300);

    // 集計と判定
    long cnt[4][C_NCLS] = {}; long ncols[4] = {};      // 集合: 0 rand、1 edge、2 zero、3 拘束 (vfix・tfix)
    auto set_ix = [](const std::string& s) { return s == "rand" ? 0 : s == "edge" ? 1 : s == "zero" ? 2 : 3; };
    double e_res_max = 0; long nzero = 0;
    for (const ColRec& r : recs) {
        const int k = set_ix(r.set); ++ncols[k]; ++cnt[k][r.cls]; if (r.zero) ++nzero;
        if (r.cls == C_PASS && !r.zero) e_res_max = std::max({e_res_max, r.e1, r.e2});
    }
    const char* SETN[4] = {"rand", "edge", "zero", "拘束 (vfix・tfix)"};
    bool fail_deriv = false, nf_j = false;
    for (int k = 0; k < 4; ++k) { if (cnt[k][C_FAIL]) fail_deriv = true; if (cnt[k][C_NF_J]) nf_j = true; }
    const long bad_rand = cnt[0][C_UNRES] + cnt[0][C_STRADDLE] + cnt[0][C_NF_FD];
    bool indet = bad_rand > ncols[0] / 100 || g_nonfinite > 0;
    for (int k = 1; k < 4; ++k) if (cnt[k][C_UNRES] + cnt[k][C_STRADDLE] + cnt[k][C_NF_FD]) indet = true;
    if (nzero == 0) indet = true;                         // 零列の判定が一度も通っていない
    std::vector<std::string> fails;
    if (fail_deriv) fails.push_back("微分の不一致 (1′)/(3′)");
    if (nf_j) fails.push_back("J が非有限");
    if (e_orient > 1e-12) fails.push_back("面の向き (1″)");
    if (e_null > 1e-12 || e_wallT > 1e-12 || e_flag != 0.0) fails.push_back("零空間・等温壁 (2)/(3)");
    if (e_line > 1e-10) fails.push_back("短いラインの解 (4)");
    if (e_f32 > 1e-4) fails.push_back("float と double (5)");

    std::printf("(1)  中心差分との差 (行列の最大で正規化、記録だけ) %.3e\n", e_fd);
    std::printf("(1′) 列ごと (標本 %zu、拘束の標本 200×2):\n", trials.size());
    for (int k = 0; k < 4; ++k) {
        std::printf("     %-18s 列 %5ld:", SETN[k], ncols[k]);
        for (int c = 0; c < C_NCLS; ++c) std::printf(" %s %ld", CLS_NAME[c], cnt[k][c]);
        std::printf("\n");
    }
    std::printf("     合格した非零列の誤差の最大 %.3e (≤ 1e-6)、零列の数 %ld\n", e_res_max, nzero);
    long shown = 0;
    for (const ColRec& r : recs) if (r.cls != C_PASS && shown < 40) {
        ++shown; std::printf("     [%s] 標本 %d・%c・列 %d: %s (e(h) %.3e、e(h/2) %.3e、再現 %.3e%s)\n", r.set.c_str(), r.sample, r.side, r.col, CLS_NAME[r.cls], r.e1, r.e2, r.rep, r.zero ? "、零列" : "");
    }
    std::printf("(1″) 面の向きの整合 %.3e (≤ 1e-12)\n", e_orient);
    std::printf("(2)  零空間 δρ(1,u,v,w,E) の作用 %.3e (≤ 1e-12)\n", e_null);
    std::printf("(3)  等温壁の拘束で ΔT_w の作用 %.3e (≤ 1e-12)、壁のフラグで K の該当行の最大 %.3e (= 0)\n", e_wallT, e_flag);
    std::printf("(4)  8 節点のライン: Thomas と密行列の解の差 %.3e (≤ 1e-10)\n", e_line);
    std::printf("(5)  float と double の差 %.3e (≤ 1e-4)、列ごと (記録だけ) %.3e\n", e_f32, e_f32_col);
    std::printf("     判定に使う量の非有限 %ld 件\n", g_nonfinite);
    if (!fails.empty()) {
        std::printf("VERDICT: FAIL (");
        for (size_t k = 0; k < fails.size(); ++k) std::printf("%s%s", k ? "、" : "", fails[k].c_str());
        std::printf(")\n");
        return 1;
    }
    if (indet) { std::printf("VERDICT: INDETERMINATE (edge・zero・拘束の未解像、rand の未解像 %ld > 1 %%、非有限、または零列の判定が無い)\n", bad_rand); return 2; }
    std::printf("VERDICT: PASS (解像した対象で不一致を検出しなかった。rand の未解像 %ld 列は保留)\n", bad_rand);
    return 0;
}
