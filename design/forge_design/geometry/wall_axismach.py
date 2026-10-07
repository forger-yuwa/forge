r"""axis-Mach チェーンの CFD ドメイン壁と壁 QA (①風洞, 無次元 r* = 1)。

計画: plans/accepted/tooling-nozzle-axismach-chain.md §5.3 +
plans/active/tooling-nozzle-axismach-throat-characteristic.md (A8)。壁構成 (上流→下流):

  入口直管 (r_U) → U→T 5次 Hermite (`UpstreamThroatPoly` 流用, r''(T)=1/R)
  → 逆 MOC 壁流線 [x_T, x_F] (端条件クランプ 5 次 B-spline = `ModeFWall` と同じ流儀)

**スロート T (x=0) から下流は全て MOC の出力**で、円弧も放物線も挟まない
(ユーザ指定 2026-08-15)。初期値線を**スロート特性線**にしたことで壁流線が
スロート壁点そのものから始まるため、旧実装にあった [T, x0] の骨接放物線区間は
不要になった (2026-08-15)。T での接続は spline 左端クランプ (r'=0, r''=1/R) と
Hermite の端点条件 (同じ 0 と 1/R) の一致で C²。

旧・縦 starting line 構成では壁流線が x0>0 から始まるため、[T, x0] を Hall 模型が
仮定する骨接放物線 r = 1 + x²/(2R) で埋めていた (`throat_start=False` で当時の
挙動を再現できる)。

`wall_qa` は逆 MOC 壁テーブルの品質指標 (単調性・最大壁角・出口角・x_F/r_F の
理論比較・曲率) を返す (原方針 §8.3, §28)。

物理壁を入口から出口まで 1 本の 5 次 B-spline に作り直す `SingleBSplinePhysicalWall` と、run に保存した壁
(`wall_repr.json`) の書き出し・復元 (`save_wall_file` / `load_wall_file`) は末尾
(plans/accepted/tooling-nozzle-wall-single-bspline.md)。

joint 壁の物理壁の上流の作り方 (`geometry.pw_upstream`: `poly` = 配管〜設計スロートの 5 次多項式 [既定] / `ramp` = 旧来の
δ_r ランプ) と、物理スロートを物理壁の大域最小として求める `wall_global_min` も末尾
(plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md §4.1)。
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline, make_interp_spline

from .wall_walldriven import UpstreamThroatPoly


def area_ratio_isentropic(M_d: float, gamma=1.4) -> float:
    """A_e/A_t (1D 等エントロピー)。r_F/r_t = sqrt(この値)。
    `gamma` にガスモデル (`GasSemiPerfect`) を渡せば thermally-perfect のテーブル値。"""
    if hasattr(gamma, "area_ratio"):
        return float(gamma.area_ratio(M_d))
    g = gamma
    return float((1.0 / M_d) * ((2.0 + (g - 1.0) * M_d * M_d) / (g + 1.0))
                 ** ((g + 1.0) / (2.0 * (g - 1.0))))


class AxisMachCFDWall:
    """直管 + U→T Hermite + 逆 MOC 壁 (クランプ 5 次 B-spline)。
    旧・縦 starting line 構成では U→T Hermite と設計壁の間に骨接放物線が入る。

    `mesh2d.generate_axisym_mesh` / `paste_isentropic_ic` 互換
    (`x_in` / `x_e` / `r(x, deriv)`)。
    """

    def __init__(self, wall_pts, R: float, r_U: float = 2.5, L_U: float = 3.5,
                 L_pipe: float = 0.5) -> None:
        """wall_pts: (n,>=2) [x, r] — 逆 MOC 壁テーブル。先頭点は
        スロート (x=0, r=1) か、旧構成では x0>0 (骨接放物線上)。"""
        wall_pts = np.asarray(wall_pts, dtype=float)
        if wall_pts.ndim != 2 or len(wall_pts) < 10:
            raise ValueError("wall_pts は (n>=10, >=2) のテーブル")
        # θ 列があれば保持 (validate のリンギング検査に使う)。壁の構築自体は (x,r) のみ
        self._th_tbl = wall_pts[:, 2].copy() if wall_pts.shape[1] >= 3 else None
        self._x_tbl = wall_pts[:, 0].copy()
        self.R = float(R)
        self.up = UpstreamThroatPoly(r_U=float(r_U), R_t=self.R, L_U=float(L_U))
        self.L_pipe = float(L_pipe)
        self.x_in = -self.up.L_U - self.L_pipe
        self.x0 = float(wall_pts[0, 0])
        if self.x0 < -1e-9:
            raise ValueError(f"設計壁始点 x0 = {self.x0:.4g} < 0 (スロート上流)")
        self.throat_start = self.x0 < 1e-6
        if self.throat_start:
            if abs(float(wall_pts[0, 1]) - 1.0) > 5e-3:
                raise ValueError(f"スロート始点なのに r = {wall_pts[0, 1]:.5f} ≠ 1")
        else:
            # [旧] 縦 starting line 構成: 始点が骨接放物線に乗っているか
            r0_par = 1.0 + self.x0 ** 2 / (2.0 * self.R)
            if abs(float(wall_pts[0, 1]) - r0_par) > 5e-3:
                raise ValueError(f"設計壁始点 r = {wall_pts[0, 1]:.5f} が放物線 "
                                 f"{r0_par:.5f} から乖離")
        self.x_e = float(wall_pts[-1, 0])
        # 左端 = スロート (r'=0, r''=1/R) / 旧構成では放物線の解析微分にクランプ
        d0 = self.x0 / self.R
        s0 = 1.0 / self.R
        _cs = CubicSpline(wall_pts[:, 0], wall_pts[:, 1])
        self._spl = make_interp_spline(
            wall_pts[:, 0], wall_pts[:, 1], k=5,
            bc_type=([(1, d0), (2, s0)],
                     [(1, float(_cs(self.x_e, 1))), (2, float(_cs(self.x_e, 2)))]))

    def r(self, x, deriv: int = 0):
        x = np.asarray(x, dtype=float)
        xU = -self.up.L_U
        out = np.empty_like(x)
        m_pipe = x < xU
        m_up = (x >= xU) & (x < 0.0)
        m_par = (x >= 0.0) & (x < self.x0)
        m_dsg = x >= self.x0
        if deriv == 0:
            out[m_pipe] = self.up.r_U
            out[m_par] = 1.0 + x[m_par] ** 2 / (2.0 * self.R)
        elif deriv == 1:
            out[m_pipe] = 0.0
            out[m_par] = x[m_par] / self.R
        elif deriv == 2:
            out[m_pipe] = 0.0
            out[m_par] = 1.0 / self.R
        elif deriv == 3:
            # 直管・骨接放物線は deriv>=3 で恒等的に 0 (壁曲率診断用、2026-08-16 追加。
            # 設計壁区間 [x0, x_F] は throat_char 構成で x0=0 なので m_par は通常空)。
            out[m_pipe] = 0.0
            out[m_par] = 0.0
        else:
            raise ValueError("deriv は 0..3")
        if m_up.any():
            out[m_up] = self.up.r(x[m_up], deriv)
        if m_dsg.any():
            out[m_dsg] = self._spl(np.minimum(x[m_dsg], self.x_e), deriv) \
                if deriv else self._spl(np.minimum(x[m_dsg], self.x_e))
        return out

    def theta(self, x):
        return np.arctan(self.r(x, 1))

    def kappa(self, x):
        rp, rpp = self.r(x, 1), self.r(x, 2)
        return rpp / (1.0 + rp * rp) ** 1.5

    def validate(self, n: int = 4000) -> list:
        msgs = ["U→T: " + m for m in self.up.validate()]
        xs = np.linspace(self.x_in, self.x_e, n)
        rv = self.r(xs)
        if np.any(rv <= 0.0):
            msgs.append("壁半径が非正")
        if abs(float(rv.min()) - 1.0) > 5e-3:
            msgs.append(f"最小半径 {float(rv.min()):.4f} != 1")
        m_dsg = xs >= self.x0
        if np.any(np.diff(rv[m_dsg]) < -1e-9):
            msgs.append("設計壁区間で半径が非単調")
        # spline リンギング検査 (2026-08-16 追加): 補間スプラインはテーブル点を通るが、
        # 端点クランプとデータが矛盾すると**点の間**で振動する。テーブル点上の
        # θ = arctan(spline') とテーブル θ (MOC 出力) の乖離を上限 0.2° で検査
        # (実測: 健全な設計は ≤0.05°、破綻した L_c=4 は 12.7°)。
        if self._th_tbl is not None:
            th_spl = np.arctan(self._spl(self._x_tbl, 1))
            dev = float(np.max(np.abs(np.degrees(th_spl - self._th_tbl))))
            if dev > 0.2:
                msgs.append(f"壁 spline リンギング (テーブル点上の θ 乖離 {dev:.2f}° > 0.2°)")
        # 接合の C1/C2 (構成的に成り立つはずだが実測で保証 — ModeFWall と同じ流儀)
        h = 1e-6
        joints = [("直管/U→T", -self.up.L_U)]
        joints += ([("U→T/設計壁 (スロート)", 0.0)] if self.throat_start
                   else [("U→T/放物線", 0.0), ("放物線/設計壁", self.x0)])
        for name, xc in joints:
            dl = float(self.r(np.array([xc - h]), 1)[0])
            dr_ = float(self.r(np.array([xc + h]), 1)[0])
            if abs(dl - dr_) > 5e-3:
                msgs.append(f"{name} の接線不連続 ({dl:.4f} vs {dr_:.4f})")
            cl = float(self.r(np.array([xc - h]), 2)[0])
            cr = float(self.r(np.array([xc + h]), 2)[0])
            if abs(cl - cr) > 0.05 * max(abs(cl), abs(cr), 1.0):
                msgs.append(f"{name} の曲率不連続 ({cl:.4f} vs {cr:.4f})")
        return msgs


def throat_curvature_fit(wall, x_max: float = 0.25) -> tuple:
    r"""MOC 壁テーブル自身が**スロートで要求する曲率** $\kappa_0$ のロバスト推定。

    $r'(0)=0$ を使い $r-1=\frac{\kappa_0}{2}x^2$ を $0<x<x_{max}$ で最小二乗
    (1 パラメータ。生の 2 階差分は非一様間隔で雑音支配になるため使わない)。
    戻り: (κ₀, フィット最大残差, 使用点数)。点不足なら (nan, nan, n)。
    """
    w = np.asarray(wall, dtype=float)
    m = (w[:, 0] > 0.0) & (w[:, 0] < x_max)
    if int(m.sum()) < 5:
        return float("nan"), float("nan"), int(m.sum())
    x, r = w[m, 0], w[m, 1]
    A = (x ** 2 / 2.0)[:, None]
    k0, *_ = np.linalg.lstsq(A, r - 1.0, rcond=None)
    resid = float(np.max(np.abs(A @ k0 - (r - 1.0))))
    return float(k0[0]), resid, int(m.sum())


def wall_qa(wall, M_d: float, x_E: float, gamma: float = 1.4,
            R: float | None = None) -> dict:
    r"""逆 MOC 壁テーブル (n,4)[x,r,θ,M] の品質指標 (原方針 §8.3, §28)。

    - 単調性 / 最大壁角 / 出口壁角 (リップで θ→0 に戻っているか)
    - $x_F, r_F$ と理論予測 ($r_F=\sqrt{A_e/A_t}$、$x_F-x_E\approx r_F\sqrt{M_d^2-1}$)
      の比較 — terminal Mach line 近似は一様域でのみ厳密なので比率は情報値
    - 壁 M の終端値 (設計 $M_d$ との差)
    - **スロート曲率の整合** (`R` を渡した時のみ、2026-08-16 追加): MOC 壁が要求する
      $\kappa_0$ (`throat_curvature_fit`) が遷音速モデルの $1/R$ と ±15% で一致し、
      2 次フィット残差が小さいこと。破れは「壁の問題」ではなく**設計の自己不整合**
      (遷音速場の仮定と超音速壁の要求が矛盾 — $L_c$ が短すぎる兆候。実測で
      $L_c\le5$ を設計段階で検出し、CFD の内部衝撃波の崖 [$L_c\le5.5$] とほぼ一致)。
      壁表現は補間スプライン + スロート端の上流クランプを維持する (ユーザ判断
      2026-08-16: 近似スプライン化は矛盾を壁の変形に隠すだけなので不採用)。
    violations が空なら合格。
    """
    w = np.asarray(wall, dtype=float)
    x, r, th = w[:, 0], w[:, 1], w[:, 2]
    r_F, x_F = float(r[-1]), float(x[-1])
    r_F_pred = float(np.sqrt(area_ratio_isentropic(M_d, gamma)))
    dxEF_pred = r_F_pred * float(np.sqrt(M_d * M_d - 1.0))
    v = []
    if np.any(np.diff(r) < -1e-9):
        v.append("壁半径が非単調")
    th_exit_deg = float(np.rad2deg(th[-1]))
    if abs(th_exit_deg) > 0.2:
        v.append(f"出口壁角 {th_exit_deg:.3f}° (|θ|>0.2° — F 未到達)")
    if abs(r_F / r_F_pred - 1.0) > 0.03:
        v.append(f"r_F = {r_F:.4f} が 1D 理論 {r_F_pred:.4f} から 3% 超乖離")
    k0 = k0_resid = k0_ratio = None
    if R is not None:
        k0, k0_resid, _n = throat_curvature_fit(w)
        if np.isfinite(k0):
            k0_ratio = k0 * float(R)
            if abs(k0_ratio - 1.0) > 0.15:
                v.append(f"スロート曲率不整合: MOC 壁の要求 κ₀={k0:.3f} が 1/R={1/R:.3f} "
                         f"から {100*abs(k0_ratio-1):.0f}% 乖離 (設計の自己不整合 — "
                         "L_c が短すぎる兆候)")
            if k0_resid > 2e-3:
                v.append(f"スロート近傍が放物線から逸脱 (fit 残差 {k0_resid:.1e} > 2e-3)")
    out = {
        "x_F": x_F, "r_F": r_F, "r_F_pred_1d": r_F_pred,
        "r_F_err_rel": float(r_F / r_F_pred - 1.0),
        "xF_minus_xE": x_F - float(x_E),
        "xF_minus_xE_pred": dxEF_pred,
        "theta_max_deg": float(np.rad2deg(th.max())),
        "theta_exit_deg": th_exit_deg,
        "M_wall_exit": float(w[-1, 3]),
        "kappa0_fit": k0, "kappa0_R_ratio": k0_ratio, "kappa0_fit_resid": k0_resid,
        "violations": v,
    }
    return out


class PhysicalNozzleWall:
    r"""**物理壁** (A13): 非粘性設計壁 $C_I$ + 上流履歴込み $\delta^*(s)$ の法線オフセット。

    計画: plans/active/tooling-nozzle-axismach-physical-throat.md。

    A12 との違い (外部レビュー採択 2026-08-17):

    1. **δ\* の弧長は入口起点** (縁 Mach は上流 = $A/A^*$ の亜音速枝、下流 = MOC 壁の M)。
       スロートで $\delta^*_t > 0$ となり、収縮部からの境界層発達履歴が入る。
    2. **真の幾何スロートを探索**: オフセット後の輪郭は $d\delta^*/ds > 0$ のため
       最小半径点が設計スロートの**やや上流**へ動き、そこでの壁角が 0 でなくなる。
       $r_W' = 0$ の点 $(x_{t,W}, r_{t,W})$ を放物線フィットで求め、その曲率
       $r''_{t,W}$ を下流スプラインの左端クランプと上流 Hermite の端条件に渡す
       (**下流がマスター、上流は作り直し** — スロート C² はユーザ要求どおり両側一致)。
    3. 遷音速の設計基準 (Hall・軸 law) は据え置き。物理壁で Hall を解き直すと
       二重計上になる。実際の遷音速場の変化は CFD で測り、必要ならアンカー帰還。

    インターフェースは `AxisMachCFDWall` 互換 (`x_in`/`x_e`/`r(x, deriv)`/`validate`)。
    追加属性 `x_throat`/`r_throat` (IC の 1D 等エントロピーが参照)。
    """

    def __init__(self, design_wall, wall_tbl, rt_m: float, Pt: float, Tt: float,
                 gamma: float = 1.4, cp: float = 1004.5, dstar_x=None,
                 x_offset_lo: float = -0.8, offset: str = "normal", delta_r_x=None,
                 analytic: bool | None = None, ramp: tuple | None = None,
                 ramp_gate_d2: float = 5e-3, upstream: str | None = None,
                 poly_gate_d2: float = 5e-3) -> None:
        """design_wall: AxisMachCFDWall (非粘性)。wall_tbl: (n,4) MOC 壁 [x,r,θ,M]。
        dstar_x: 任意の δ*(x) [r_t 単位] (None = 上流履歴込み相関)。
        offset: "normal" (旧: 壁法線オフセット + x シフト) / "radial" (半径方向 r_W = r + δ_r(x)、
        plans/active/tooling-nozzle-deltastar-core-matched-euler.md §4.5 の生産経路)。
        delta_r_x: 半径方向補正 δ_r(x) [r_t 単位] を直接与える (与えると offset="radial" 固定)。
        analytic: 解析経路 (下記)。None = 設計壁が `JointFitCFDWall` (wall_repr: joint) のとき True。
        ramp: 解析経路の δ_r を入れる区間 [x_lo, x_hi] (x ≤ x_lo で 0、x ≥ x_hi で全量、間は 5 次 smoothstep)。
        problem の `geometry.pw_ramp`。None = 既定「直管接合の直後 (設計縮流部で r′ < −0.05 になる最初の x) から −0.5·L_U」
        (`default_pw_ramp`)。
        ramp_gate_d2: ランプ区間のゲート。|r″ − r″_design| ≤ ramp_gate_d2 かつ r′ < 0 でなければ例外で止める
        (短い縮流部ではランプが δ_r·s″ で縮流部の曲率を壊すため; codex result M5)。

        **解析経路** (plan tooling-nozzle-cfd-pinned-initial-line §5.1 #6b, 仕様 A″ diagnostician 2026-10-05):
        κ_t の窓 LSQ 再推定・上流 Hermite の作り直し・オフセット点群の補間スプラインでの作り直しをせず、
        r_W(x) = r_design(x) + s(x)·δ_r(x) をそのまま使う (r′・r″・r‴ も設計壁の解析微分 + (s·δ_r) の
        解析微分)。δ_r は導関数を返せる callable (`delta_r_x(x, deriv)`, `supports_deriv`) に限る — 差分では代用しない。
        物理スロート (x_t, r_t, κ_t) は r_W′ = 0 の根。

        upstream: 解析経路の上流 (x < 0) の作り方 (problem の `geometry.pw_upstream`、plan tooling-nozzle-upstream-poly-and-throat-sizing
        §4.1)。None = 既定 `poly` (ユーザ決定 2026-10-07)。
        - `poly`: 直管 r_U ([x_in, −L_U))・5 次多項式 Q ([−L_U, 0)、端条件は x = −L_U で (r_U, 0, 0)、x = 0 で下流の物理壁 S + δ_r の
          (r, r′, r″))・S + δ_r ([0, x_e]、`ramp` とビット同一)。上流に δ_r を足さない。`ramp` (pw_ramp) を併せて渡したら例外。
          物理スロートは物理壁の大域最小 (`wall_global_min`、固定の囲い込みを使わない)。ゲート (`upstream_gate`): 正の半径・
          一意な最小・その前で r′ ≤ 0・後で r′ ≥ 0、[−L_U, 0) で |Q″ − H″| ≤ poly_gate_d2、継ぎ目 (−L_U・0) の値・1 階・2 階微分の
          跳び ≤ `POLY_SEAM_TOL`。不合格は例外。
        - `ramp`: 上の s(x)·δ_r (旧来の作り方。`ramp` を明示したときだけ)。変更前とビット同一。
        解析経路でない (joint でない) 壁に `poly` を渡したら例外。"""
        from ..evaluate.ic import invert_area_ratio
        from ..feedback.deltastar import dstar_flatplate
        if upstream is not None and not (isinstance(upstream, str) and upstream in PW_UPSTREAMS):
            raise ValueError(f"PhysicalNozzleWall: upstream は {' / '.join(repr(c) for c in PW_UPSTREAMS)} のどれか (受け取った値: {upstream!r})")
        if not (analytic if analytic is not None else isinstance(design_wall, JointFitCFDWall)) and upstream == "poly":
            raise ValueError("PhysicalNozzleWall: upstream 'poly' は joint 壁の物理壁の解析経路だけ "
                             f"(設計壁 {type(design_wall).__name__}, analytic={analytic})")
        self.design = design_wall
        self.r_U = float(design_wall.up.r_U)
        self.L_U = float(design_wall.up.L_U)
        self.x_in = float(design_wall.x_in)
        wall_tbl = np.asarray(wall_tbl, dtype=float)
        x_F = float(wall_tbl[-1, 0])

        # --- δ*(x): 入口からの弧長 + 区分縁 Mach ------------------------------
        xg = np.linspace(self.x_in + 1e-6, x_F, 3000)
        rg = design_wall.r(xg)
        rpg = design_wall.r(xg, 1)
        s_g = np.concatenate([[0.0], np.cumsum(
            np.hypot(np.diff(xg), np.diff(rg)))]) * rt_m
        M_g = np.empty_like(xg)
        m_up = xg <= 0.0
        AR_up = np.maximum(rg[m_up], 1.0 + 1e-12) ** 2
        gas_obj = gamma if hasattr(gamma, "area_ratio") else None
        if gas_obj is not None:
            from ..evaluate.ic import _invert_area_ratio_gas
            M_g[m_up] = _invert_area_ratio_gas(AR_up, np.zeros(int(m_up.sum()), bool), gas_obj)
            g_corr = float(gas_obj.gamma_throat(Tt))       # 相関 (Eckert) は代表 γ で十分
        else:
            M_g[m_up] = invert_area_ratio(AR_up, np.zeros(int(m_up.sum()), bool), gamma)
            g_corr = float(gamma)
        M_g[~m_up] = np.interp(xg[~m_up], wall_tbl[:, 0], wall_tbl[:, 3])
        ds_g = dstar_flatplate(s_g, M_g, Pt, Tt, g_corr, cp) / rt_m
        self._dstar_hist = lambda x: np.interp(x, xg, ds_g)
        dstar = self._dstar_hist if dstar_x is None else dstar_x

        if analytic is None:
            analytic = isinstance(design_wall, JointFitCFDWall)
        self.analytic = bool(analytic)
        if self.analytic:
            self._init_analytic(design_wall, delta_r_x, ramp, xg, x_offset_lo, ramp_gate_d2,
                                upstream=upstream, poly_gate_d2=poly_gate_d2)
            return
        # 解析経路でない壁は上流を作り直す従来経路 (pw_upstream の概念が無い; `ramp` の明示は効かない)
        self.pw_upstream = None
        self.pw_upstream_source = None if upstream is None else f"explicit {upstream!r} (解析経路でないので無効)"

        # --- 法線オフセット (窓 [x_offset_lo, x_F] — 真のスロート探索を含む) ---
        mw = xg >= x_offset_lo
        th_w = np.arctan(rpg[mw])
        if delta_r_x is not None:
            offset = "radial"
        self.offset_mode = offset
        if offset == "radial":
            drv = np.asarray((delta_r_x if delta_r_x is not None else dstar)(xg[mw]), dtype=float)
            xw = xg[mw].copy()
            rw = rg[mw] + drv
            self._delta_r_applied = lambda x, _x=xw.copy(), _d=drv.copy(): np.interp(x, _x, _d)
        elif offset == "normal":
            dsv = np.asarray(dstar(xg[mw]), dtype=float)
            xw = xg[mw] - dsv * np.sin(th_w)
            rw = rg[mw] + dsv * np.cos(th_w)
            o = np.argsort(xw)
            xw, rw = xw[o], rw[o]
            self._delta_r_applied = lambda x, _x=xw.copy(), _r=rw.copy(), _dw=design_wall: \
                np.interp(x, _x, _r) - _dw.r(np.asarray(x, dtype=float))
        else:
            raise ValueError("offset は 'normal' か 'radial'")
        self._xw_dbg, self._rw_dbg = xw.copy(), rw.copy()   # 感度診断用に保持

        # --- 真の幾何スロート: オフセット点群を spline 化 → r'(x)=0 を root solve ---
        # (Codex 指摘 2026-08-17: 放物線フィットは近似探索。spline の導関数の零点を
        #  直接解き、その点で r, r'' を評価する。感度は sensitivity_report() で別途)
        self.x_throat, self.r_throat, self.kappa_throat, self._throat_diag = \
            self._locate_throat(xw, rw)
        if not (-0.5 < self.x_throat < 0.2):
            raise RuntimeError(f"PhysicalNozzleWall: 真のスロート x={self.x_throat:.3f} "
                               "が想定窓 (-0.5, 0.2) の外")

        # --- 下流: スロートクランプ付き補間 5 次 B-spline ----------------------
        md = xw > self.x_throat + 0.02
        xd = np.concatenate([[self.x_throat], xw[md]])
        rd = np.concatenate([[self.r_throat], rw[md]])
        _cs = CubicSpline(xd, rd)
        self.x_e = float(xd[-1])
        self._spl = make_interp_spline(
            xd, rd, k=5,
            bc_type=([(1, 0.0), (2, self.kappa_throat)],
                     [(1, float(_cs(self.x_e, 1))), (2, float(_cs(self.x_e, 2)))]))

        # --- 上流: 入口 → 新スロートの quintic Hermite (作り直し) -------------
        from .wall import _hermite_quintic, _poly_eval
        self._herm_x0 = -self.L_U
        self._herm_c = _hermite_quintic(self._herm_x0, self.x_throat,
                                        self.r_U, 0.0, 0.0,
                                        self.r_throat, 0.0, self.kappa_throat)
        self._poly_eval = _poly_eval

    def _init_analytic(self, design_wall, delta_r_x, ramp, xg, x_offset_lo, ramp_gate_d2: float = 5e-3,
                       upstream: str | None = None, poly_gate_d2: float = 5e-3) -> None:
        """解析経路の構築 (`__init__` の docstring 参照)。"""
        from scipy.optimize import brentq
        mode = PW_UPSTREAM_DEFAULT if upstream is None else upstream
        if mode == "poly" and ramp is not None:
            raise ValueError(f"PhysicalNozzleWall (解析経路): pw_upstream 'poly' ({'既定' if upstream is None else '明示'}) と "
                             f"pw_ramp {ramp} の併記は不可 — 旧来のランプを使うなら pw_upstream: ramp を明示する "
                             "(plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1)")
        if delta_r_x is None:
            raise ValueError("PhysicalNozzleWall (解析経路 / joint 壁): delta_r_x (半径方向 δ_r) が必須 "
                             "— 相関 δ* の法線オフセットは解析経路に無い")
        if not getattr(delta_r_x, "supports_deriv", False):
            raise ValueError("PhysicalNozzleWall (解析経路 / joint 壁): δ_r が導関数を返せない "
                             "(delta_r_from_table の関数を使う。積分法の smooth_delta_quintic は未対応 — 差分で代用しない)")
        self.pw_upstream = mode
        self.pw_upstream_source = "default" if upstream is None else "explicit"
        if mode == "poly":
            self._init_poly(design_wall, delta_r_x, xg, x_offset_lo, poly_gate_d2)
            return
        if ramp is None:
            ramp = default_pw_ramp(design_wall)
            self._ramp_source = "default"
        else:
            self._ramp_source = "pw_ramp"
        if len(ramp) != 2:
            raise ValueError(f"pw_ramp {ramp} は 2 要素 [x_lo, x_hi]")
        lo, hi = float(ramp[0]), float(ramp[1])
        if not (-self.L_U <= lo < hi < 0.0):
            raise ValueError(f"pw_ramp {ramp} は [−L_U, 0) = [{-self.L_U:g}, 0) 内の増加区間であること")
        self._ramp = (lo, hi)
        self._dr = delta_r_x
        self.offset_mode = "radial"
        self.x_e = float(design_wall.x_e)
        self._herm_x0 = -self.L_U                              # validate の上流区間 (設計の U→T Hermite 始点)
        self._delta_r_applied = lambda x: self._sdr(np.asarray(x, dtype=float), 0)
        mw = xg >= x_offset_lo
        self._xw_dbg, self._rw_dbg = xg[mw].copy(), self.r(xg[mw])
        d1 = lambda x: float(self.r(np.array([x]), 1)[0])  # noqa: E731
        xl, xr = -0.3, 0.2
        if not (d1(xl) < 0.0 < d1(xr)):
            raise RuntimeError("PhysicalNozzleWall (解析経路): r′=0 の囲い込みに失敗 (−0.3, 0.2)")
        self.x_throat = float(brentq(d1, xl, xr, xtol=1e-14))
        self.r_throat = float(self.r(np.array([self.x_throat]))[0])
        self.kappa_throat = float(self.r(np.array([self.x_throat]), 2)[0])
        # ランプ区間のゲート: 縮流部の曲率を δ_r·s″ で壊していないか・収縮が単調か
        xr = np.linspace(lo, hi, 6001)
        d2 = float(np.max(np.abs(self.r(xr, 2) - design_wall.r(xr, 2))))
        r1 = float(np.max(self.r(xr, 1)))
        self.ramp_gate = {"ramp": [lo, hi], "source": self._ramp_source, "max_abs_d2_change": d2,
                          "max_r1": r1, "limit_d2": float(ramp_gate_d2),
                          "pass": bool(d2 <= ramp_gate_d2 and r1 < 0.0)}
        if not self.ramp_gate["pass"]:
            raise ValueError(f"PhysicalNozzleWall (解析経路): pw_ramp [{lo:g}, {hi:g}] のゲート不合格 — "
                             f"|r″ − r″_design| max {d2:.2e} (≤ {ramp_gate_d2:g}), max r′ {r1:.2e} (< 0) "
                             "(縮流部が短くランプが曲率を壊す: geometry.pw_ramp を広げるか L_U を見直す)")
        self._throat_diag = {"method": "analytic", "ramp": [lo, hi], "ramp_gate": self.ramp_gate}

    def _init_poly(self, design_wall, delta_r_x, xg, x_offset_lo, poly_gate_d2: float) -> None:
        """解析経路の上流 `poly` の構築 (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1)。
        Q は x = −L_U で (r_U, 0, 0)、x = 0 で下流の物理壁 (S + δ_r) の x = 0⁺ の (r, r′, r″) を満たす 5 次 Hermite
        (6 条件で一意、ξ = (x + L_U)/L_U 基底)。x ≥ 0 は `ramp` と同じ式 (設計壁 + δ_r) でビット同一。"""
        from .wall import _hermite_quintic, _poly_eval
        L = self.L_U
        self._dr = delta_r_x
        self.offset_mode = "radial"
        self.x_e = float(design_wall.x_e)
        self._herm_x0 = -L                                     # validate の上流区間 (多項式 Q の始点)
        x0 = np.array([0.0])
        e0 = [float(design_wall.r(x0, n)[0] + delta_r_x(x0, n)[0]) for n in range(3)]    # x = 0⁺ (S + δ_r)
        self._q_c = _hermite_quintic(-L, 0.0, self.r_U, 0.0, 0.0, e0[0], e0[1], e0[2])
        self._q_eval = _poly_eval
        self.upstream_poly = {"coef": [float(v) for v in self._q_c], "basis": "Σ coef[j] ξ^j, ξ = (x + L_U)/L_U",
                              "domain": [-L, 0.0],
                              "end_conditions": {"x_minus_L_U": [self.r_U, 0.0, 0.0], "x_0": e0,
                                                 "note": "x = −L_U で直管 (r_U, 0, 0)、x = 0 で下流の物理壁 S + δ_r の x = 0⁺ の (r, r′, r″)"}}
        self._delta_r_applied = lambda x: self.r(np.asarray(x, dtype=float)) - self.design.r(np.asarray(x, dtype=float))
        mw = xg >= x_offset_lo
        self._xw_dbg, self._rw_dbg = xg[mw].copy(), self.r(xg[mw])
        # --- 物理スロート = 物理壁の大域最小 (Q と下流 S + δ_r の全ノット区間) ---
        srch = wall_global_min(self.r, self.poly_breaks())
        self.x_throat, self.r_throat, self.kappa_throat = srch["x_t"], srch["r_t"], srch["kappa_t"]
        # --- 上流の曲率の変化 |Q″ − H″| ([−L_U, 0)、差は ξ の 3 次式なので区間端と極値で厳密) ---
        d2max, x_d2 = _quintic_d2_diff_absmax(self._q_c, np.asarray(design_wall.up._c, dtype=float), L)
        # --- 継ぎ目の値・1 階・2 階微分の跳び (左右の極限) ---
        qL = [float(_poly_eval(self._q_c, -L, 0.0, np.array([-L]), n)[0]) for n in range(3)]
        q0 = [float(_poly_eval(self._q_c, -L, 0.0, np.array([0.0]), n)[0]) for n in range(3)]
        seams = [{"x": -L, "left": [self.r_U, 0.0, 0.0], "right": qL,
                  **{f"jump_d{n}": float(abs(qL[n] - (self.r_U if n == 0 else 0.0))) for n in range(3)}},
                 {"x": 0.0, "left": q0, "right": e0, **{f"jump_d{n}": float(abs(e0[n] - q0[n])) for n in range(3)}}]
        seam_max = max(s_[f"jump_d{n}"] for s_ in seams for n in range(3))
        self.upstream_gate = _poly_gate_record(srch, d2max, x_d2, poly_gate_d2, seams, seam_max,
                                               source=self.pw_upstream_source, repr_="legacy")
        self.upstream_gate["vs_ramp_default"] = self._diff_vs_default_ramp(design_wall)
        self._throat_diag = {"method": "global_min", "pw_upstream": "poly", "search": srch}
        if not self.upstream_gate["pass"]:
            raise ValueError("PhysicalNozzleWall (解析経路、pw_upstream poly): 上流のゲート不合格 — "
                             + _poly_gate_message(self.upstream_gate))

    def poly_breaks(self) -> np.ndarray:
        """`poly` の物理壁が区間ごとに 5 次以下の多項式になる区切り: 端・継ぎ目 (−L_U・0)・設計壁 S のノット・δ_r の補間スプラインの
        ノットと表の端 (表の外は端値一定)。`wall_global_min` に渡す。"""
        S = self.design._spl
        br = [self.x_in, -self.L_U, 0.0, float(self.design.x0), self.x_e]      # x0: 旧・縦 starting line の放物線と S の継ぎ目
        br += list(np.asarray(S.t, dtype=float))
        sp = getattr(self._dr, "spline", None)
        if sp is None or not hasattr(self._dr, "x_range"):
            raise ValueError("PhysicalNozzleWall (pw_upstream poly): δ_r の関数に spline・x_range が無い "
                             "(delta_r_from_table の関数を使う — 区間の区切りが決まらない)")
        br += list(np.asarray(sp.t, dtype=float)) + [float(v) for v in self._dr.x_range]
        br = np.unique(np.asarray(br, dtype=float))
        return br[((br >= 0.0) & (br <= self.x_e)) | np.isin(br, [self.x_in, -self.L_U])]

    def _diff_vs_default_ramp(self, design_wall) -> dict:
        """診断: 旧来のランプ (既定の `default_pw_ramp`) の物理壁との差 (x < 0 だけ。x ≥ 0 は同じ式)。ゲートには使わない。"""
        try:
            lo, hi = default_pw_ramp(design_wall)
        except ValueError as e:
            return {"status": f"既定ランプが決まらない ({e})"}
        import types
        ns = types.SimpleNamespace(_ramp=(lo, hi))
        xs = np.linspace(self.x_in, 0.0, 24001)[:-1]
        d = self.r(xs) - (design_wall.r(xs) + PhysicalNozzleWall._s(ns, xs, 0) * self._dr(xs, 0))
        i = int(np.argmax(np.abs(d)))
        return {"ramp_ref": [lo, hi], "ramp_ref_source": "default_pw_ramp", "n_points": int(len(xs)),
                "max_abs_diff_rt": float(abs(d[i])), "x_at_max_abs_diff": float(xs[i]),
                "min_signed_rt": float(d.min()), "max_signed_rt": float(d.max()),
                "note": "poly − ramp (既定ランプ、均等点)。生産の pw_ramp と比べるときは検証スクリプトで"}

    def _s(self, x, n: int):
        """δ_r の入れ方 s(x): 5 次 smoothstep (端で 1・2 階微分 0)。n 階微分。"""
        lo, hi = self._ramp
        L = hi - lo
        u = np.clip((x - lo) / L, 0.0, 1.0)
        m = (x > lo) & (x < hi)
        if n == 0:
            return u ** 3 * (10.0 - 15.0 * u + 6.0 * u * u)
        if n == 1:
            return np.where(m, 30.0 * u * u * (1.0 - u) ** 2 / L, 0.0)
        if n == 2:
            return np.where(m, 60.0 * u * (1.0 - u) * (1.0 - 2.0 * u) / L ** 2, 0.0)
        if n == 3:
            return np.where(m, 60.0 * (1.0 - 6.0 * u + 6.0 * u * u) / L ** 3, 0.0)
        raise ValueError("deriv は 0..3")

    def _sdr(self, x, n: int):
        """(s·δ_r) の n 階微分 (Leibniz)。s = 0 の区間 (直管・縮流部入口側) は厳密に 0。"""
        out = np.zeros_like(x)
        on = x > self._ramp[0]
        if on.any():
            xo = x[on]
            binom = (1, 1, 1, 1) if n == 0 else ((1, 1), (1, 2, 1), (1, 3, 3, 1))[n - 1]
            out[on] = sum(binom[k] * self._dr(xo, k) * self._s(xo, n - k) for k in range(n + 1))
        return out

    @staticmethod
    def _locate_throat(xw, rw, window: float = 0.6, smooth_lam: float = 1e-9,
                       kappa_window: float = 0.15):
        r"""オフセット点群から真の幾何スロート $(x_t, r_t, r''_t)$ を求める。

        **位置と曲率で推定器を分ける** (2026-08-17 感度測定に基づく):

        - **位置 $x_t, r_t$**: 点群を微小 λ の平滑化スプラインにし $r'(x)=0$ を Brent 法で
          root solve。λ を 1e-9〜1e-3 で振っても $x_t$ は ±0.004、$r_t$ は 5 桁不変。
          放物線フィットは窓 0.08〜0.40 で $x_t$ が −0.013〜+0.006 と動くので位置には使わない。
        - **曲率 $r''_t$**: 法線オフセットした点群は $x$ 間隔が非一様で $r''$ に ±10% の
          ジッタを持つ (生 2 階差分の IQR [0.46, 0.61])。スプラインの $r''$ は λ に
          4.36→0.53 と極端に依存し信頼できない。$r'(x_t)=0$ を既知として
          $r - r_t = \frac{\kappa}{2}(x-x_t)^2$ の 1 パラメータ LSQ (窓 ±kappa_window) で
          robust 推定する — 窓 0.08〜0.40 で 0.526〜0.552 と安定、オフセット理論値
          $\kappa_I/(1-\kappa_I\delta^*)=0.503$ + $\delta^{*\prime\prime}$ 項と整合。

        戻り: (x_t, r_t, r''_t, diag)。diag に両推定器の感度を含める。"""
        from scipy.interpolate import make_smoothing_spline
        from scipy.optimize import brentq
        i0 = int(np.argmin(rw))
        m = np.abs(xw - xw[i0]) < window
        if int(m.sum()) < 12:
            raise RuntimeError("PhysicalNozzleWall: スロート近傍の点が不足")
        spl = make_smoothing_spline(xw[m], rw[m], lam=smooth_lam)
        d1 = spl.derivative(1)
        xl, xr = float(xw[max(i0 - 3, 0)]), float(xw[min(i0 + 3, len(xw) - 1)])
        for _ in range(6):
            if d1(xl) < 0.0 < d1(xr):
                break
            xl -= 0.02; xr += 0.02
        else:
            raise RuntimeError("PhysicalNozzleWall: r'=0 の囲い込みに失敗")
        x_t = float(brentq(lambda x: float(d1(x)), xl, xr, xtol=1e-12))
        r_t = float(spl(x_t))
        # 曲率: r'(x_t)=0 拘束の 1 パラメータ LSQ
        mk = np.abs(xw - x_t) < kappa_window
        A = (0.5 * (xw[mk] - x_t) ** 2)[:, None]
        k_t = float(np.linalg.lstsq(A, rw[mk] - r_t, rcond=None)[0][0])
        # 感度診断
        ks = {}
        for kw in (0.08, 0.15, 0.30):
            mm = np.abs(xw - x_t) < kw
            AA = (0.5 * (xw[mm] - x_t) ** 2)[:, None]
            ks[kw] = float(np.linalg.lstsq(AA, rw[mm] - r_t, rcond=None)[0][0])
        xs_ = {}
        for lam in (1e-9, 1e-5, 1e-3):
            sp2 = make_smoothing_spline(xw[m], rw[m], lam=lam)
            try:
                xs_[lam] = float(brentq(lambda x: float(sp2.derivative(1)(x)), xl, xr))
            except ValueError:
                xs_[lam] = float("nan")
        diag = {"kappa_by_window": ks, "x_throat_by_lambda": xs_,
                "kappa_spread": float(max(ks.values()) - min(ks.values())),
                "x_throat_spread": float(np.nanmax(list(xs_.values()))
                                         - np.nanmin(list(xs_.values()))),
                "n_pts": int(m.sum()), "kappa_window": kappa_window}
        return x_t, r_t, k_t, diag

    def r(self, x, deriv: int = 0):
        x = np.asarray(x, dtype=float)
        if self.analytic:
            if self.pw_upstream == "poly":
                # 直管 r_U / 多項式 Q / S + δ_r (x ≥ 0 は ramp の design.r + s·δ_r と同じ値 — x ≥ x_hi で s ≡ 1)
                if deriv not in (0, 1, 2, 3):
                    raise ValueError("deriv は 0..3")
                out = np.empty_like(x)
                m_pipe = x < -self.L_U
                m_q = (x >= -self.L_U) & (x < 0.0)
                m_dn = x >= 0.0
                out[m_pipe] = self.r_U if deriv == 0 else 0.0
                if m_q.any():
                    out[m_q] = self._q_eval(self._q_c, -self.L_U, 0.0, x[m_q], deriv)
                if m_dn.any():
                    xd = x[m_dn]
                    out[m_dn] = self.design.r(xd, deriv) + self._dr(xd, deriv)
                return out
            return self.design.r(x, deriv) + self._sdr(x, deriv)
        out = np.empty_like(x)
        m_pipe = x < self._herm_x0
        m_up = (x >= self._herm_x0) & (x < self.x_throat)
        m_dn = x >= self.x_throat
        out[m_pipe] = self.r_U if deriv == 0 else 0.0
        if m_up.any():
            out[m_up] = self._poly_eval(self._herm_c, self._herm_x0, self.x_throat,
                                        x[m_up], deriv)
        if m_dn.any():
            xq = np.minimum(x[m_dn], self.x_e)
            out[m_dn] = self._spl(xq, deriv) if deriv else self._spl(xq)
        return out

    def theta(self, x):
        return np.arctan(self.r(x, 1))

    def validate(self, n: int = 4000) -> list:
        msgs = []
        xs = np.linspace(self.x_in, self.x_e, n)
        rv = self.r(xs)
        if np.any(rv <= 0.0):
            msgs.append("壁半径が非正")
        if abs(float(rv.min()) - self.r_throat) > 5e-3:
            msgs.append(f"最小半径 {float(rv.min()):.4f} != r_throat {self.r_throat:.4f}")
        m_dn = xs >= self.x_throat
        if np.any(np.diff(rv[m_dn]) < -1e-9):
            msgs.append("スロート下流で半径が非単調")
        m_up = (xs >= self._herm_x0) & (xs <= self.x_throat)
        if np.any(np.diff(rv[m_up]) > 1e-9):
            msgs.append("上流 Hermite が非単調収縮")
        # C1/C2 接合 (構成的に成り立つはず — 実測で保証)
        h = 1e-6
        joints = ((("直管/縮流部", self._herm_x0), ("縮流部/設計壁 (x=0)", 0.0)) if self.analytic else
                  (("直管/Hermite", self._herm_x0), ("Hermite/下流 (物理スロート)", self.x_throat)))
        for name, xc in joints:
            dl = float(self.r(np.array([xc - h]), 1)[0])
            dr_ = float(self.r(np.array([xc + h]), 1)[0])
            if abs(dl - dr_) > 5e-3:
                msgs.append(f"{name} の接線不連続 ({dl:.5f} vs {dr_:.5f})")
            cl = float(self.r(np.array([xc - h]), 2)[0])
            cr = float(self.r(np.array([xc + h]), 2)[0])
            if abs(cl - cr) > 0.05 * max(abs(cl), abs(cr), 1.0):
                msgs.append(f"{name} の曲率不連続 ({cl:.4f} vs {cr:.4f})")
        return msgs


# --- A14: 制約付き最小二乗 B-spline (本流の形状表現候補、補間壁と A/B) -----------
def lsq_bspline_wall(wall_tbl, n_cp: int, kappa_t: float, k: int = 5,
                     theta_e: float | None = None):
    r"""MOC 壁テーブル (n,4)[x,r,θ,M] を **制約付き最小二乗 B-spline** $r(x)$ で近似する。

    計画: plans/active/tooling-nozzle-axismach-physical-throat.md (A14)。

    - 未知: 制御点 $c_i$ ($i=1..n_{cp}$)、ノットは弧長で等分配 (端は $k+1$ 重)。
    - **ハード拘束** (等式): $r(x_t)=r_t$, $r'(x_t)=0$, $r''(x_t)=\kappa_t$ (上流 Hermite と
      C² — ユーザ要求), $r(x_e)=r_e$, 任意で $r'(x_e)=\tan\theta_e$。
    - **重み** $w_j \sim \Delta s_j$ (MOC の点密度が形状を支配しないよう弧長重み)。
    - 解法: KKT 系 (正規方程式 + ラグランジュ乗数)。

    戻り: `scipy.interpolate.BSpline` と診断 dict (max|Δr|・max|Δθ| [deg]・
    曲率振動 $\int\kappa'^2 ds$・拘束残差)。raw MOC 点は呼び出し側で保持し、
    診断を**別ゲート**として監視する (LSQ が不整合を隠さないため — Codex 指摘)。
    """
    from scipy.interpolate import BSpline
    w = np.asarray(wall_tbl, dtype=float)
    x, r, th = w[:, 0], w[:, 1], w[:, 2]
    x_t, r_t, x_e, r_e = float(x[0]), float(r[0]), float(x[-1]), float(r[-1])
    # 弧長重み
    ds = np.hypot(np.diff(x), np.diff(r))
    wt = np.concatenate([[ds[0]], 0.5 * (ds[1:] + ds[:-1]), [ds[-1]]])
    wt = wt / wt.mean()
    # 弧長等分配ノット (内部ノット n_cp - k - 1 個)
    s_cum = np.concatenate([[0.0], np.cumsum(ds)])
    n_int = n_cp - k - 1
    if n_int < 0:
        raise ValueError(f"n_cp={n_cp} は次数 k={k} に対して少なすぎる (最低 {k+1})")
    s_int = np.linspace(0.0, s_cum[-1], n_int + 2)[1:-1]
    x_int = np.interp(s_int, s_cum, x)
    t = np.concatenate([[x_t] * (k + 1), x_int, [x_e] * (k + 1)])
    # 設計行列
    B = BSpline.design_matrix(x, t, k).toarray()            # (n, n_cp)
    def _row(xq, d):
        # 基底関数の d 階微分行 (design_matrix は nu 非対応 → 単位係数で評価)
        n_cp_ = len(t) - k - 1
        return np.array([BSpline(t, np.eye(n_cp_)[i], k)(xq, d) for i in range(n_cp_)])
    C = [_row(x_t, 0), _row(x_t, 1), _row(x_t, 2), _row(x_e, 0)]
    dvec = [r_t, 0.0, float(kappa_t), r_e]
    if theta_e is not None:
        C.append(_row(x_e, 1)); dvec.append(float(np.tan(theta_e)))
    C = np.asarray(C); dvec = np.asarray(dvec)
    # KKT
    W = wt[:, None]
    A = B.T @ (W * B)
    b = B.T @ (wt * r)
    m = len(dvec)
    K = np.block([[A, C.T], [C, np.zeros((m, m))]])
    rhs = np.concatenate([b, dvec])
    sol = np.linalg.lstsq(K, rhs, rcond=None)[0]
    c = sol[:n_cp]
    spl = BSpline(t, c, k)
    # 診断
    r_fit = spl(x)
    th_fit = np.arctan(spl(x, 1))
    xs = np.linspace(x_t, x_e, 4000)
    rp, rpp, rppp = spl(xs, 1), spl(xs, 2), spl(xs, 3)
    kap = rpp / (1 + rp ** 2) ** 1.5
    dk = np.gradient(kap, xs)
    ds_x = np.sqrt(1 + rp ** 2)
    J_fair = float(np.trapezoid(dk ** 2 * ds_x, xs))
    # MOC 生テーブルの同じ指標 (比較基準): 差分曲率
    rp_t = np.gradient(r, x); rpp_t = np.gradient(rp_t, x)
    kap_t = rpp_t / (1 + rp_t ** 2) ** 1.5
    dk_t = np.gradient(kap_t, x)
    J_fair_tbl = float(np.trapezoid(dk_t ** 2 * np.sqrt(1 + rp_t ** 2), x))
    diag = {"n_cp": int(n_cp), "k": int(k),
            "max_dr": float(np.max(np.abs(r_fit - r))),
            "rms_dr": float(np.sqrt(np.mean((r_fit - r) ** 2))),
            "max_dtheta_deg": float(np.degrees(np.max(np.abs(th_fit - th)))),
            "J_fair": J_fair, "J_fair_tbl": J_fair_tbl,
            "constraint_resid": float(np.max(np.abs(C @ c - dvec))),
            "monotone": bool(np.all(np.diff(spl(xs)) > -1e-9))}
    return spl, diag


def select_lsq_ncp(wall_tbl, kappa_t: float, candidates=(12, 16, 20, 24, 32),
                   tol_dr: float = 5e-4, tol_dtheta_deg: float = 0.05, k: int = 5,
                   theta_e: float | None = None) -> tuple:
    """max|Δr|・max|Δθ| のゲートを満たす**最小の制御点数**を選ぶ (固定 20 にしない)。
    戻り: (spline, diag, sweep_table)。どれも満たさなければ最大 n_cp を返し
    diag["gate_ok"]=False。"""
    table = []
    best = None
    for n in candidates:
        try:
            spl, dg = lsq_bspline_wall(wall_tbl, n, kappa_t, k=k, theta_e=theta_e)
        except Exception as e:  # noqa: BLE001
            table.append({"n_cp": n, "error": str(e)}); continue
        dg["gate_ok"] = (dg["max_dr"] <= tol_dr and dg["max_dtheta_deg"] <= tol_dtheta_deg
                         and dg["monotone"] and dg["constraint_resid"] < 1e-8)
        table.append(dg)
        best = (spl, dg)
        if dg["gate_ok"]:
            break
    return best[0], best[1], table


class LSQBsplineCFDWall(AxisMachCFDWall):
    r"""**A14 の CFD 壁**: 設計区間を `lsq_bspline_wall` の近似 B-spline で表現する
    `AxisMachCFDWall` 亜種。上流 (直管 + U→T Hermite) は同一、スロート端は同じ
    $(r'=0,\ r''=1/R)$ ハード拘束で C² 接続。**メッシュ生成・CFD・CAD 出力が同じ
    スプラインを使う** (Codex: エクスポート層のみは CAD≠CFD 形状になるため不採用)。
    `fit_diag` に近似誤差ゲート (max|Δr|・max|Δθ|・J_fair) を保持し、raw MOC 点との
    乖離は別ゲートとして監視する。"""

    def __init__(self, wall_tbl, R: float, n_cp=None, r_U: float = 2.5,
                 L_U: float = 3.5, L_pipe: float = 0.5, k: int = 5,
                 tol_dr: float = 5e-4, tol_dtheta_deg: float = 0.05) -> None:
        wall_tbl = np.asarray(wall_tbl, dtype=float)
        # 親の検査 (throat_start・接合) を通してから設計区間の表現だけ差し替える
        super().__init__(wall_tbl, R=R, r_U=r_U, L_U=L_U, L_pipe=L_pipe)
        if n_cp is None:
            spl, dg, sweep = select_lsq_ncp(wall_tbl, 1.0 / R, k=k, tol_dr=tol_dr,
                                            tol_dtheta_deg=tol_dtheta_deg,
                                            theta_e=float(wall_tbl[-1, 2]))
            self.fit_sweep = sweep
        else:
            spl, dg = lsq_bspline_wall(wall_tbl, int(n_cp), 1.0 / R, k=k,
                                       theta_e=float(wall_tbl[-1, 2]))
            self.fit_sweep = [dg]
        self._spl = spl
        self.fit_diag = dg
        self.n_cp = int(dg["n_cp"])
        # 親の 0.2° リンギング検査は「補間スプラインの点間振動」用。LSQ 壁では
        # テーブル点上の θ 乖離は意図した近似残差そのものなので、その基準は
        # fit_diag (max_dtheta_deg / max_dr) 側で別ゲートとして扱う (Codex: 別々のゲート)。
        self._th_tbl = None
        self.tol_dtheta_deg = float(tol_dtheta_deg)
        self.tol_dr = float(tol_dr)

    def validate(self, n: int = 4000) -> list:
        msgs = super().validate(n)
        dg = self.fit_diag
        if dg["max_dr"] > self.tol_dr:
            msgs.append(f"LSQ 近似誤差 max|Δr|={dg['max_dr']:.2e} > {self.tol_dr:.1e}")
        if dg["max_dtheta_deg"] > self.tol_dtheta_deg:
            msgs.append(f"LSQ 近似誤差 max|Δθ|={dg['max_dtheta_deg']:.3f}° > {self.tol_dtheta_deg}°")
        return msgs


def default_pw_ramp(design_wall, r1_on: float = -0.05) -> tuple:
    """`pw_ramp` の既定: 直管接合の直後 (設計縮流部で r′ < r1_on になる最初の x) から −0.5·L_U。"""
    L_U = float(design_wall.up.L_U)
    xs = np.linspace(-L_U, 0.0, 200001)
    r1 = design_wall.r(xs, 1)
    i = int(np.argmax(r1 < r1_on))
    if not r1[i] < r1_on:
        raise ValueError(f"default_pw_ramp: 設計縮流部で r′ < {r1_on} の点が無い")
    return (float(xs[i]), -0.5 * L_U)


# --- 上流の作り方と物理スロートの大域最小 (plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md §4.1) -----------
PW_UPSTREAMS = ("ramp", "poly")      # problem の geometry.pw_upstream
PW_UPSTREAM_DEFAULT = "poly"         # joint 壁の物理壁の解析経路の既定 (ユーザ決定 2026-10-07)。`ramp` は明示したときだけ
POLY_SEAM_TOL = 1e-8                 # 継ぎ目 (−L_U・0) の値・1 階・2 階微分の跳びの上限 (左右の極限)
THROAT_R1_TOL = 1e-12                # 単調性の判定で r′ の符号を見るときの丸めの許容 (r′ は無次元)


def wall_global_min(f, breaks, deg: int = 4, r1_tol: float = THROAT_R1_TOL, xtol: float = 1e-14) -> dict:
    r"""物理壁 r(x) の**大域最小** (物理スロート) を区間ごとの多項式から求め、一意性と前後の単調性を判定する
    (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1、plan 段レビュー M2)。

    f: `f(x, deriv)` (deriv 0..2 を使う)。breaks: r が各区間で次数 ≤ deg + 1 の多項式になる区切り (端・継ぎ目・ノット)。
    固定の囲い込みや「最小は Q の中にある」の仮定は使わない。

    1. 各区間で r′ を Chebyshev 点 deg + 1 個から次数 deg の多項式として復元する (多項式なので厳密)。区間内の別の 2 点で照合し、
       合わなければ区切りが足りないとして例外。r′ ≡ 0 の区間 (直管) は平らな区間として扱う。
    2. r′ の実根 (区間内、符号が変われば Brent 法で詰める)・区間の端・継ぎ目を候補にし、弱い極小 (下の 4.) のうち r が最小の
       候補を x_t とする (極小でない候補の方が丸めの許容を超えて低ければ不合格)。
    3. 単調性: 区間を r′ の実根で切った小区間の中点で r′ の符号を見る (多項式は根の間で符号が変わらない)。x_t より前で
       r′ ≤ r1_tol、後で r′ ≥ −r1_tol。
    4. 一意性: 弱い極小 (左で r′ ≤ r1_tol かつ右で r′ ≥ −r1_tol の候補) のうち、r が最小値 + 1e-12·max(1, |r_t|) 以内のものが
       1 か所だけ (1e-12 以内の重複は 1 か所とみなす)。
    戻り: x_t・r_t・kappa_t (= r″(x_t)) と判定 (positive・unique・monotone_before・monotone_after・pass) と記録。"""
    from numpy.polynomial import chebyshev as Ch
    from scipy.optimize import brentq
    br = np.unique(np.asarray(breaks, dtype=float))
    if len(br) < 2 or not np.all(np.isfinite(br)):
        raise ValueError("wall_global_min: 区切りが 2 点未満・非有限")
    a, b = br[:-1], br[1:]
    n = len(a)
    mid, half = 0.5 * (a + b), 0.5 * (b - a)
    m = deg + 1
    uc = np.cos(np.pi * (np.arange(m) + 0.5) / m)
    xc = mid[:, None] + half[:, None] * uc[None, :]
    v = np.asarray(f(xc.ravel(), 1), dtype=float).reshape(xc.shape)
    coef = np.linalg.solve(Ch.chebvander(uc, deg), v.T)                 # (deg+1, n)
    uk = np.array([-0.77, 0.41])                                        # 照合点 (Chebyshev 点と別)
    xk = mid[:, None] + half[:, None] * uk[None, :]
    vk = np.asarray(f(xk.ravel(), 1), dtype=float).reshape(xk.shape)
    pk = Ch.chebval(uk, coef)                                           # (n, 2)
    scl = np.maximum(np.abs(v).max(axis=1), np.abs(vk).max(axis=1))
    resid = np.abs(vk - pk).max(axis=1)
    rel = resid / np.maximum(scl, 1e-300)
    bad = resid > 1e-8 * scl + 1e-12                                    # 絶対の床 1e-12: r ~ 10 の B-spline の微分の丸め
    if bad.any():
        j = int(np.argmax(np.where(bad, rel, -1.0)))
        raise ValueError(f"wall_global_min: 区間 [{a[j]:.9g}, {b[j]:.9g}] で r′ が次数 {deg} の多項式でない "
                         f"(照合の差 {resid[j]:.2e}, r′ の大きさ {scl[j]:.2e}) — 区切り (ノット・継ぎ目) が足りない")

    def d1(x):
        return float(np.asarray(f(np.array([x]), 1), dtype=float)[0])

    flat = np.zeros(n, dtype=bool)
    roots = []                                                          # (x, 区間 j)
    for j in range(n):
        c = coef[:, j]
        cm = float(np.abs(c).max())
        if cm <= 1e-15:
            flat[j] = True
            continue
        ct = Ch.chebtrim(c, tol=1e-13 * cm)
        if len(ct) < 2:
            continue
        rr = Ch.chebroots(ct)
        rr = rr[np.abs(rr.imag) <= 1e-6].real                          # 重根・接する根は小さい虚部で出る
        rr = np.clip(rr[(rr >= -1.0 - 1e-9) & (rr <= 1.0 + 1e-9)], -1.0, 1.0)
        for u in np.unique(rr):
            x = float(mid[j] + half[j] * u)
            dl, dh = max(float(a[j]), x - 1e-6 * float(half[j])), min(float(b[j]), x + 1e-6 * float(half[j]))
            fl, fh = d1(dl), d1(dh)
            if fl * fh < 0.0:
                x = float(brentq(d1, dl, dh, xtol=xtol))
            roots.append((x, j))
    # 候補 = 区切り + 根 (1e-12 以内は 1 点に)
    cand = np.unique(np.r_[br, [x for x, _ in roots]])
    keep = np.r_[True, np.diff(cand) > 1e-12 * np.maximum(1.0, np.abs(cand[1:]))]
    cand = cand[keep]
    rc = np.asarray(f(cand, 0), dtype=float)
    # 小区間 (候補で切った区間) の中点での r′ の符号と、弱い極小 (左で r′ ≤ tol かつ右で r′ ≥ −tol; 端は片側だけ)
    sub_a, sub_b = cand[:-1], cand[1:]
    sm = 0.5 * (sub_a + sub_b)
    r1m = np.asarray(f(sm, 1), dtype=float)
    left = np.r_[-np.inf, r1m]                                          # 候補 i の左の小区間の r′ (i = 0 は無し)
    right = np.r_[r1m, np.inf]
    is_min = (left <= r1_tol) & (right >= -r1_tol)
    # x_t = 弱い極小のうち r が最小のもの。全候補の最小で選ぶと、丸めで r が並ぶ点 (根のすぐ隣の継ぎ目など) を拾って
    # 単調性の判定を誤る (2026-10-07 U1 の負例: δ_r′(0) = 1.6e-11 で根 −3.1e-11 と継ぎ目 0 の r が丸めで同じ)
    idx = np.where(is_min)[0]
    i_min = int(idx[int(np.argmin(rc[idx]))]) if len(idx) else int(np.argmin(rc))
    x_t, r_t = float(cand[i_min]), float(rc[i_min])
    global_below = float(rc.min()) < r_t - 1e-12 * max(1.0, abs(r_t))      # 極小でない候補の方が低い (起こらないはず — 不合格にする)
    before, after = sm < x_t, sm > x_t
    viol_b = np.where(before, r1m, -np.inf)
    viol_a = np.where(after, -r1m, -np.inf)
    max_r1_before = float(viol_b.max()) if before.any() else None
    min_r1_after = float(-viol_a.max()) if after.any() else None
    mono_b = (max_r1_before is None) or max_r1_before <= r1_tol
    mono_a = (min_r1_after is None) or min_r1_after >= -r1_tol
    tol_r = 1e-12 * max(1.0, abs(r_t))
    at_min = np.where(is_min & (rc <= r_t + tol_r))[0]
    out = {"x_t": x_t, "r_t": r_t, "kappa_t": float(np.asarray(f(np.array([x_t]), 2), dtype=float)[0]),
           "r1_at_x_t": d1(x_t),
           "n_intervals": int(n), "n_flat_intervals": int(flat.sum()), "n_roots": int(len(roots)),
           "roots_x": [float(x) for x, _ in sorted(roots)][:64],
           "local_minima": [[float(cand[i]), float(rc[i])] for i in np.where(is_min)[0]][:64],
           "n_at_min": int(len(at_min)), "unique": bool(len(at_min) == 1),
           "positive": bool(r_t > 0.0),
           "monotone_before": bool(mono_b), "max_r1_before": max_r1_before,
           "x_max_r1_before": (float(sm[int(np.argmax(viol_b))]) if before.any() else None),
           "monotone_after": bool(mono_a), "min_r1_after": min_r1_after,
           "x_min_r1_after": (float(sm[int(np.argmax(viol_a))]) if after.any() else None),
           "r1_tol": float(r1_tol), "poly_check_max_rel": float(rel.max()), "deg_r1": int(deg),
           "x_t_is_local_min": bool(len(idx) > 0), "global_min_below_x_t": bool(global_below)}
    out["pass"] = bool(out["positive"] and out["unique"] and mono_b and mono_a and len(idx) > 0 and not global_below)
    return out


def _quintic_d2_diff_absmax(cq, ch, L: float) -> tuple:
    """[−L, 0) で |Q″ − H″| の最大 (Q・H は ξ = (x + L)/L 基底の 5 次の係数)。差の 2 階微分は ξ の 3 次式なので、区間端と
    その導関数 (2 次式) の実根で厳密に求める。戻り: (最大, その x)。"""
    P = np.polynomial.Polynomial
    d2 = P(np.asarray(cq, dtype=float) - np.asarray(ch, dtype=float)).deriv(2) / (L * L)
    crit = [0.0, 1.0]
    if d2.degree() >= 2:
        rr = d2.deriv().roots()
        crit += [float(z.real) for z in np.atleast_1d(rr) if abs(z.imag) <= 1e-12 and 0.0 < z.real < 1.0]
    vals = np.abs(d2(np.asarray(crit)))
    i = int(np.argmax(vals))
    return float(vals[i]), float(-L + L * crit[i])


def _poly_absmax(fun, a: float, b: float, deg: int) -> tuple:
    """[a, b] で次数 ≤ deg の多項式 `fun(x)` の |·| の最大 (Chebyshev 点で復元し、区間端と導関数の実根で評価する)。戻り: (最大, x)。"""
    from numpy.polynomial import chebyshev as Ch
    m = deg + 1
    uc = np.cos(np.pi * (np.arange(m) + 0.5) / m)
    mid, half = 0.5 * (a + b), 0.5 * (b - a)
    c = Ch.chebfit(uc, np.asarray(fun(mid + half * uc), dtype=float), deg)
    crit = [-1.0, 1.0]
    if deg >= 2:
        rr = Ch.chebroots(Ch.chebder(c))
        crit += [float(z.real) for z in np.atleast_1d(rr) if abs(z.imag) <= 1e-9 and -1.0 < z.real < 1.0]
    xs = mid + half * np.asarray(crit)
    vals = np.abs(np.asarray(fun(xs), dtype=float))
    i = int(np.argmax(vals))
    return float(vals[i]), float(xs[i])


def _poly_gate_record(srch: dict, d2max: float, x_d2: float, limit_d2: float, seams: list, seam_max: float,
                      source, repr_: str) -> dict:
    """`poly` の上流のゲートの記録 (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1)。"""
    rec = {"pw_upstream": "poly", "source": source, "repr": repr_,
           "throat": {"x": srch["x_t"], "r": srch["r_t"], "kappa": srch["kappa_t"]},
           "positive_radius": srch["positive"], "unique_min": srch["unique"], "n_at_min": srch["n_at_min"],
           "monotone_before": srch["monotone_before"], "max_r1_before": srch["max_r1_before"],
           "monotone_after": srch["monotone_after"], "min_r1_after": srch["min_r1_after"], "r1_tol": srch["r1_tol"],
           "max_abs_d2_change_vs_H": float(d2max), "x_at_max_d2_change": float(x_d2), "limit_d2": float(limit_d2),
           "limit_d2_note": "幾何の変化の上限 (今のランプのゲートと同じ値)。流れの品質の保証ではない",
           "seam_jumps": seams, "max_seam_jump": float(seam_max), "limit_seam": POLY_SEAM_TOL,
           "throat_search": {k: srch[k] for k in ("n_intervals", "n_flat_intervals", "n_roots", "roots_x", "local_minima",
                                                   "x_max_r1_before", "x_min_r1_after", "poly_check_max_rel", "r1_at_x_t")}}
    rec["pass"] = bool(srch["pass"] and d2max <= limit_d2 and seam_max <= POLY_SEAM_TOL)
    return rec


def _poly_gate_message(g: dict) -> str:
    return (f"正の半径 {g['positive_radius']}, 一意な最小 {g['unique_min']} ({g['n_at_min']} か所), "
            f"前で r′ ≤ 0 {g['monotone_before']} (max {g['max_r1_before']}), 後で r′ ≥ 0 {g['monotone_after']} "
            f"(min {g['min_r1_after']}), |Q″ − H″| max {g['max_abs_d2_change_vs_H']:.3e} (≤ {g['limit_d2']:g}), "
            f"継ぎ目の跳び max {g['max_seam_jump']:.2e} (≤ {g['limit_seam']:g}); 物理スロート x = {g['throat']['x']:.6g}")


# --- 位置 + 壁角の同時当てはめ (wall_repr: joint, V0 型壁) ---------------------------
def joint_fit_wall(wall_tbl, R: float, lam: float = 1e-9, k: int = 5, sig_r: float = 1e-6,
                   sig_th: float = 1e-4, h0: float = 0.0125, h1: float = 0.5, x_g: float = 6.0,
                   mono_r2=None, diag: dict | None = None):
    r"""MOC 壁点 (n,>=3)[x,r,θ] に 5 次 B-spline $r(x)$ を**位置と壁角の両方**で当てはめる (V0 型壁)。

    計画: plans/active/tooling-nozzle-cfd-pinned-initial-line.md §4.5・§5.1 #6 (試作は
    case/45.isobutane_m6_d155/moc_wall_fit_ab.py::joint_fit、同じ式・同じ既定値)。

    目的関数 $\sum_j w_j[(r(x_j)-r_j)^2/\sigma_r^2 + (r'(x_j)-\tan\theta_j)^2/\sigma_\theta^2]
    + (\lambda/\sigma_r^2)\int (r''')^2 dx$ ($w_j\propto\Delta x_j$)。ハード拘束 (KKT):
    $r(x_0)=r_0$, $r'(x_0)=x_0/R$, $r''(x_0)=1/R$, $r(x_e)=r_e$, $r'(x_e)=\tan\theta_e$。
    ノット間隔は始点 h0 → x_g で h1 へ smoothstep、以降 h1 (最後の内部ノットは x_e − h1/2 まで)。
    戻り: (BSpline, 制御点数)。

    `mono_r2=(a, b)` (plans/accepted/tooling-nozzle-throat-monotone-r2.md §4.1): r‴ (2 次スプライン) の B-spline 係数の
    うち台 [t₃ⱼ, t₃ⱼ₊₃] が (a, b) にかかるものを ≤ 0 に拘束する (凸包性により r″ が [a, b] で単調非増加の**十分条件**)。
    2 次形式を対角最大で、不等式行列を行ごとの最大絶対値で正規化し、有効制約法 (違反最大を追加・負の乗数を除去、
    上限 200 回で例外) で解く。停止条件は違反 max(Gc) ≤ 1e-10·max(1, max|Gc|)、乗数 ≥ −1e-12·max(1, max|μ|)
    (解法用。形状の単調性の保証は区間多項式で別に検査する)。`mono_r2=None` は従来の連立方程式 1 回 (正規化なし) で
    ビット同一。`diag` (dict) を渡すと、等式残差・ノット・係数・次数 (と mono 時は KKT の記録) を書き込む。"""
    from scipy.interpolate import BSpline
    tb = np.asarray(wall_tbl, dtype=float)
    x, r, th = tb[:, 0], tb[:, 1], tb[:, 2]
    x0, xe = x[0], x[-1]
    if mono_r2 is not None:
        # 不正入力は当てはめ前に拒否 (黙って拘束なし・空拘束で解かない)
        try:
            mono = [float(v) for v in mono_r2]
        except (TypeError, ValueError):
            raise ValueError(f"mono_r2 は 2 要素の数値列 (a, b): {mono_r2!r}") from None
        if isinstance(mono_r2, (str, bytes)) or len(mono) != 2 or any(isinstance(v, (bool, np.bool_)) for v in mono_r2):
            raise ValueError(f"mono_r2 は 2 要素の数値列 (a, b): {mono_r2!r}")
        a_m, b_m = mono
        if not (np.isfinite(a_m) and np.isfinite(b_m)):
            raise ValueError(f"mono_r2 に非有限値: {mono_r2!r}")
        if not a_m < b_m:
            raise ValueError(f"mono_r2 は a < b: {mono_r2!r}")
        if a_m < x0 or b_m > xe:
            raise ValueError(f"mono_r2 = {mono_r2!r} が設計壁の x 範囲 [{x0:.6g}, {xe:.6g}] の外")
    xs = [x0]
    while xs[-1] < xe:
        u = min((xs[-1] - x0) / (x_g - x0), 1.0)
        xs.append(xs[-1] + h0 + (h1 - h0) * u * u * (3 - 2 * u))
    xi = np.array(xs[1:-1])
    xi = xi[xi < xe - 0.5 * h1]
    t = np.r_[[x0] * (k + 1), xi, [xe] * (k + 1)]
    nc = len(t) - k - 1
    E = np.eye(nc)

    def D(xq, d):
        return np.array([BSpline(t, E[i], k)(xq, d) for i in range(nc)]).T
    B0, B1 = D(x, 0), D(x, 1)
    w = np.gradient(x)
    w = w / w.mean()
    A = (B0.T * (w / sig_r ** 2)) @ B0 + (B1.T * (w / sig_th ** 2)) @ B1
    b = B0.T @ (w * r / sig_r ** 2) + B1.T @ (w * np.tan(th) / sig_th ** 2)
    if lam > 0:
        xq = np.linspace(x0, xe, 8000)
        B3 = D(xq, 3)
        A = A + lam * (B3.T * np.gradient(xq)) @ B3 / sig_r ** 2
    Cm = np.array([D(np.r_[x0], 0)[0], D(np.r_[x0], 1)[0], D(np.r_[x0], 2)[0],
                   D(np.r_[xe], 0)[0], D(np.r_[xe], 1)[0]])
    dv = np.array([r[0], x0 / R, 1.0 / R, r[-1], np.tan(th[-1])])
    kkt = None
    if mono_r2 is None:
        K = np.block([[A, Cm.T], [Cm, np.zeros((5, 5))]])
        c = np.linalg.solve(K, np.r_[b, dv])[:nc]
    else:
        c, kkt = _joint_fit_mono_qp(A, b, Cm, dv, t, k, nc, (a_m, b_m))
    if diag is not None:
        diag.update({"eq_resid": float(np.max(np.abs(Cm @ c - dv))),
                     "spline": {"k": int(k), "t": [float(v) for v in t], "c": [float(v) for v in c]}})
        if kkt is not None:
            diag.update(kkt)
    return BSpline(t, c, k), nc


def _joint_fit_mono_qp(A, b, Ce, de, t, k: int, nc: int, mono, max_iter: int = 200):
    r"""`joint_fit_wall(mono_r2=…)` の有効制約法 (plan tooling-nozzle-throat-monotone-r2 §4.1)。
    等式 Ce c = de、不等式 G c ≤ 0 (G = 台が mono にかかる r‴ の B-spline 係数の行、行ごとに最大絶対値で正規化)。
    試作 case/45.isobutane_m6_d155/throat_r2_explainer.py::fit_variant と同じ式・同じ停止条件。戻り: (c, KKT 記録)。"""
    from scipy.interpolate import BSpline
    sc = float(np.abs(np.diag(A)).max())
    A = A / sc
    b = b / sc
    E = np.eye(nc)
    # r‴ (2 次スプライン, ノット t[3:-3]) の係数を c の一次式として組む: r‴ の係数 = M3 @ c
    M3 = np.array([BSpline(t, E[i], k).derivative(3).c[:nc - 3] for i in range(nc)]).T
    t3 = t[3:len(t) - 3]
    sel = [j for j in range(nc - 3) if t3[j + 3] > mono[0] + 1e-12 and t3[j] < mono[1]]
    if not sel:
        raise ValueError(f"mono_r2 = {mono!r} にかかる r‴ の係数が無い")
    G0 = M3[sel]
    gs = np.abs(G0).max(axis=1, keepdims=True)
    G = G0 / gs
    ne = len(de)
    act: list = []
    for it in range(max_iter):
        Cm = np.vstack([Ce, G[act]]) if act else Ce
        dd = np.r_[de, np.zeros(len(act))]
        m = len(dd)
        sol = np.linalg.solve(np.block([[A, Cm.T], [Cm, np.zeros((m, m))]]), np.r_[b, dd])
        c, mu = sol[:nc], sol[nc + ne:]
        if len(act) and mu.min() < -1e-12 * max(1.0, float(np.abs(mu).max())):
            act.pop(int(np.argmin(mu)))
            continue
        viol = G @ c
        if viol.max() <= 1e-10 * max(1.0, float(np.abs(viol).max())):
            break
        j = int(np.argmax(viol))
        if j in act:
            raise RuntimeError("joint_fit_wall(mono_r2): 有効制約法が循環した (有効な制約が違反最大)")
        act.append(j)
    else:
        raise RuntimeError(f"joint_fit_wall(mono_r2): 有効制約法が {max_iter} 回で収束しない")
    viol = G @ c
    kkt = {"mono_r2": [float(mono[0]), float(mono[1])], "n_constraints": int(len(sel)),
           "n_active": int(len(act)), "iters": int(it + 1),
           "ineq_max_normalized": float(viol.max()),
           "mu_min": (float(mu.min()) if len(act) else None),
           "r3_coef_max_unnormalized": float((G0 @ c).max()),
           "active_support": [[float(t3[sel[a]]), float(t3[sel[a] + 3])] for a in act],
           "qp_scale_diag_max": sc}
    return c, kkt


def r3_piecewise_exact(f, breaks, a: float, b: float, tol_sign: float = 1e-6) -> dict:
    r"""[a, b] で r‴ が区間ごとの 2 次式 (5 次スプライン、またはノットの違う 5 次スプラインの和) である壁の形状量を
    **区間多項式から厳密に**求める (plan tooling-nozzle-throat-monotone-r2 §6 S1・S6・S8 の形状用検査。均等点ではない)。

    f: `f(x, deriv)` (deriv = 2, 3 を使う) か `BSpline` (breaks 省略可)。breaks: 区分の境界 (ノット; [a, b] 外は無視)。
    各区間で r‴ を内部の Gauss 3 点から 2 次式として復元し (2 次式なので厳密)、次を返す
    (区間端の r‴ は f を直接評価する。5 次スプラインは単純ノットで C⁴ なので r‴ は区間端で連続):

    - `r3_max` / `x_r3_max`: r‴ の最大 (区間端・頂点で評価)
    - `r2_max_increase`: max_{x<y} [r″(y) − r″(x)] (r″ は r‴ の零点の間で単調なので、区間端と零点の r″ 列で厳密)
    - `r2_max` / `x_r2_max`: r″ の最大 (同じ候補点)
    - `r4_absmax`: max|r⁗| (r⁗ は区間ごとの 1 次式なので区間端の片側極限)
    - `int_r3sq`: ∫_a^b (r‴)² dx (区間 3 点 Gauss、被積分は 4 次なので厳密)
    - `n_r2_extrema`: r″ の内部極値の数 = r‴ の符号反転の数 (|r‴| < tol_sign の点は符号なしとして飛ばす)"""
    if hasattr(f, "t") and hasattr(f, "c") and hasattr(f, "k"):
        spl = f
        if breaks is None:
            breaks = spl.t
        f = lambda x, d, _s=spl: _s(x, d)  # noqa: E731
    a, b = float(a), float(b)
    if not (np.isfinite(a) and np.isfinite(b) and a < b):
        raise ValueError(f"r3_piecewise_exact: 区間 [{a}, {b}] が不正")
    br = np.unique(np.r_[a, b, np.asarray(breaks, dtype=float)])
    br = br[(br >= a) & (br <= b)]
    lo_, hi_ = br[:-1], br[1:]
    keep = hi_ - lo_ > 0.0
    lo_, hi_ = lo_[keep], hi_[keep]
    gp, gw = np.polynomial.legendre.leggauss(3)
    h = hi_ - lo_
    xq = 0.5 * h[:, None] * gp[None, :] + 0.5 * (lo_ + hi_)[:, None]          # (n, 3)
    r3q = np.asarray(f(xq.ravel(), 3), dtype=float).reshape(xq.shape)
    int_r3sq = float((r3q ** 2 * (0.5 * h[:, None] * gw[None, :])).sum())
    # 区間ごとの 2 次式 q(u) = c2 u² + c1 u + c0 (u = x − lo) を 3 点から復元
    cand_x, cand_r3, r4_end = [], [], []
    sign_seq = []
    roots_all = []
    for i in range(len(lo_)):
        u = xq[i] - lo_[i]
        c2, c1, c0 = np.linalg.solve(np.vander(u, 3), r3q[i])
        q = lambda uu, c2=c2, c1=c1, c0=c0: c2 * uu * uu + c1 * uu + c0  # noqa: E731
        pts = [0.0]
        if c2 != 0.0:
            uv = -c1 / (2.0 * c2)
            if 0.0 < uv < h[i]:
                pts.append(uv)
        pts.append(h[i])
        for j, uu in enumerate(pts):
            # 区間端は f を直接評価 (2 次式の外挿は丸めで ~1e-9 相対ずれる)、頂点だけ復元した 2 次式で
            v = q(uu) if 0 < j < len(pts) - 1 else float(np.asarray(f(np.array([lo_[i] + uu]), 3), dtype=float)[0])
            cand_x.append(lo_[i] + uu)
            cand_r3.append(v)
            sign_seq.append(v)
        r4_end += [abs(c1), abs(2.0 * c2 * h[i] + c1)]
        # r‴ の零点 (r″ の極値候補)
        if c2 != 0.0:
            disc = c1 * c1 - 4.0 * c2 * c0
            rts = [] if disc < 0 else [(-c1 - np.sqrt(disc)) / (2 * c2), (-c1 + np.sqrt(disc)) / (2 * c2)]
        elif c1 != 0.0:
            rts = [-c0 / c1]
        else:
            rts = []
        roots_all += [lo_[i] + rr for rr in rts if 0.0 < rr < h[i]]
    cand_x = np.asarray(cand_x)
    cand_r3 = np.asarray(cand_r3)
    i3 = int(np.argmax(cand_r3))
    xs2 = np.unique(np.r_[br, roots_all])
    r2 = np.asarray(f(xs2, 2), dtype=float)
    inc = float(np.max(r2 - np.minimum.accumulate(r2)))
    i2 = int(np.argmax(r2))
    sg = np.sign([v if abs(v) >= tol_sign else 0.0 for v in sign_seq])
    sg = sg[sg != 0]
    return {"a": a, "b": b, "r3_max": float(cand_r3[i3]), "x_r3_max": float(cand_x[i3]),
            "r2_max_increase": inc, "r2_max": float(r2[i2]), "x_r2_max": float(xs2[i2]),
            "r4_absmax": float(max(r4_end)), "int_r3sq": int_r3sq,
            "n_r2_extrema": int(np.count_nonzero(np.diff(sg) != 0)), "n_intervals": int(len(lo_))}


class JointFitCFDWall(AxisMachCFDWall):
    r"""**`wall_repr: joint` の CFD 壁**: 設計区間を `joint_fit_wall` (位置 + 壁角の同時当てはめ) で表現する
    `AxisMachCFDWall` 亜種。上流 (直管 + U→T Hermite) と検査 (接合・リンギング 0.2°) は親と同じ。
    `PhysicalNozzleWall` はこの壁を受けると、スロートの再推定をせず設計壁の解析 $r(x)$ に $\delta_r(x)$ を足す
    (`analytic=True` 経路)。`fit_diag` に当てはめの設定と点上の乖離 (max|Δr|・max|Δθ|) を持つ。"""

    def __init__(self, wall_tbl, R: float, r_U: float = 2.5, L_U: float = 3.5,
                 L_pipe: float = 0.5, lam: float = 1e-9, h0: float = 0.0125, h1: float = 0.5,
                 x_g: float = 6.0, mono_r2=None) -> None:
        wall_tbl = np.asarray(wall_tbl, dtype=float)
        super().__init__(wall_tbl, R=R, r_U=r_U, L_U=L_U, L_pipe=L_pipe)
        fd: dict = {}
        spl, nc = joint_fit_wall(wall_tbl, self.R, lam=lam, h0=h0, h1=h1, x_g=x_g, mono_r2=mono_r2, diag=fd)
        self._spl = spl
        xt = wall_tbl[:, 0]
        # mono_r2 (plan tooling-nozzle-throat-monotone-r2 §4.1): r″ の単調拘束区間 (None = 拘束なし、従来と同一)。
        # spline (ノット・係数・次数) は run の壁の証拠 (変換後メッシュの壁節点と解析壁の照合) に使う
        self.fit_diag = {"kind": "joint_fit", "lam": float(lam), "h0": float(h0), "h1": float(h1),
                         "x_g": float(x_g), "n_cp": int(nc),
                         "max_dr_pts": float(np.max(np.abs(spl(xt) - wall_tbl[:, 1]))),
                         "max_dtheta_pts_deg": float(np.degrees(np.max(np.abs(
                             np.arctan(spl(xt, 1)) - wall_tbl[:, 2])))),
                         "mono_r2": (None if mono_r2 is None else [float(v) for v in mono_r2]),
                         **fd}


# --- 物理壁の全域 1 本の 5 次 B-spline と保存した壁の復元 (plans/accepted/tooling-nozzle-wall-single-bspline.md) -----------
PHYSICAL_WALL_REPRS = ("legacy", "single_bspline")   # problem の geometry.physical_wall_repr (キー無し = legacy・壁ファイルなし)
WALL_FILE = "wall_repr.json"                          # run_dir 直下の壁ファイル (保存した壁の復元の取り決め §4.2)
WALL_FILE_FORMAT = "forge_design.nozzle_wall"
# 版 2 (2026-10-07, plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1): 上流の方式 `pw_upstream` と多項式 Q (係数・基底・区間)
# を保存する。版 1 (方式の記録が無い) は ramp として読む
WALL_FILE_VERSION = 2
WALL_FILE_VERSIONS_READ = (1, 2)
SINGLE_BSPLINE_TOL = {"r": 1.3e-7, "r1": 1e-7, "r2": 1e-5}   # 元の物理壁との差の許容 (要求の上限、plan wall-single-bspline §4.1)。超えたら例外
KNOT_REMOVAL_TOL = 1e-12           # 継ぎ目のノット除去の誤差の上限 [r_t] (plan upstream-poly §4.1b)。超えたら例外
THROAT_BRACKET = (-0.3, 0.2)       # `ramp` の解析経路の r′ = 0 の囲い込み (`PhysicalNozzleWall._init_analytic`)。`poly` は `wall_global_min`

def bspline_piece_limits(spl, xq, nmax: int = 4):
    """B-spline `spl` の 0..nmax 階微分の、点 xq (ノット) での左右の極限を区間多項式から厳密に求める。
    戻り: (left, right) — 形 (nmax+1, len(xq))。端 (左の区間が無い・右の区間が無い) は nan。"""
    from scipy.interpolate import PPoly
    pp = PPoly.from_spline(spl)
    xb = pp.x
    idx = np.where(np.diff(xb) > 0.0)[0]                    # 長さのある区間
    a, b = xb[idx], xb[idx + 1]
    pos_r = {float(v): i for i, v in enumerate(a)}            # 左端 = xq の区間 (右極限)
    pos_l = {float(v): i for i, v in enumerate(b)}            # 右端 = xq の区間 (左極限)
    xq = np.asarray(xq, dtype=float)
    left = np.full((nmax + 1, len(xq)), np.nan)
    right = np.full((nmax + 1, len(xq)), np.nan)
    for n in range(nmax + 1):
        c = (pp.derivative(n) if n else pp).c[:, idx]          # (deg+1, m): 区間 j の多項式 Σ c[i, j] (x − a_j)^(deg−i)
        for q, x in enumerate(xq):
            j = pos_l.get(float(x))
            if j is not None:
                left[n, q] = np.polyval(c[:, j], b[j] - a[j])
            j = pos_r.get(float(x))
            if j is not None:
                right[n, q] = c[-1, j]
    return left, right


def bspline_wall_errors(spl, f, intervals, n_dense: int = 41, n_cheb: int = 16, deg: int = 10, nmax: int = 3) -> dict:
    """1 本の B-spline `spl` と元の壁 `f(x, deriv)` の差を、各ノット区間 [a, b] の (1) 内部の密な点 (n_dense 点の内側) と
    (2) 区間ごとの差の多項式の極値 (差は区間ごとに次数 ≤ deg の多項式: ランプは 5 次 × 5 次 + 5 次 = 10 次。n_cheb 点の
    Chebyshev 当てはめで復元し、導関数の根を候補点にして差そのもの (当てはめではない) を評価) の両方で測る。
    intervals: (m, 2) の [a, b]。戻り: {n: (各区間の最大 |Δ| (m,), 各区間の最大の位置 (m,))} (n = 0..nmax)。"""
    from numpy.polynomial import chebyshev as Ch
    iv = np.asarray(intervals, dtype=float)
    a, b = iv[:, 0], iv[:, 1]
    m = len(a)
    mid, half = 0.5 * (a + b), 0.5 * (b - a)
    ud = np.linspace(-1.0, 1.0, n_dense)[1:-1]
    uc = np.cos(np.pi * (np.arange(n_cheb) + 0.5) / n_cheb)
    xd = (mid[:, None] + half[:, None] * ud[None, :])
    xc = (mid[:, None] + half[:, None] * uc[None, :])
    out = {}
    for n in range(nmax + 1):
        dd = spl(xd.ravel(), n).reshape(xd.shape) - f(xd.ravel(), n).reshape(xd.shape)
        dc = spl(xc.ravel(), n).reshape(xc.shape) - f(xc.ravel(), n).reshape(xc.shape)
        dg = max(deg - n, 1)
        V = Ch.chebvander(uc, dg)                              # (n_cheb, dg+1) — 全区間で同じ
        coef = np.linalg.lstsq(V, dc.T, rcond=None)[0]         # (dg+1, m)
        cand_u, cand_i = [], []
        for j in range(m):
            r = Ch.chebroots(Ch.chebder(coef[:, j])) if dg >= 2 else np.array([])
            r = r[np.isreal(r)].real if len(r) else r
            r = r[(r > -1.0) & (r < 1.0)] if len(r) else r
            cand_u.append(r)
            cand_i.append(np.full(len(r), j))
        cu = np.concatenate(cand_u) if cand_u else np.array([])
        ci = np.concatenate(cand_i).astype(int) if cand_i else np.array([], int)
        best = np.abs(dd).max(axis=1)
        xbest = xd[np.arange(m), np.argmax(np.abs(dd), axis=1)]
        if len(cu):
            xe = mid[ci] + half[ci] * cu
            de = np.abs(spl(xe, n) - f(xe, n))
            for j, v, x in zip(ci, de, xe):
                if v > best[j]:
                    best[j], xbest[j] = v, x
        out[n] = (best, xbest)
    return out


def _power_to_bernstein(a) -> np.ndarray:
    """ξ ∈ [0, 1] の多項式 Σ aᵢ ξ^i (次数 k = len(a) − 1) の Bernstein 係数 b_j = Σ_{i≤j} C(j, i)/C(k, i)·aᵢ (Bézier = 端の重複度
    k + 1 の B-spline の係数)。"""
    from math import comb
    a = np.asarray(a, dtype=float)
    k = len(a) - 1
    return np.array([sum(comb(j, i) / comb(k, i) * a[i] for i in range(j + 1)) for j in range(k + 1)])


def _insert_knot(t, c, k: int, u: float) -> tuple:
    """Boehm のノット挿入 (1 個)。新しい係数は隣り合う係数の凸結合だけ。u は基本区間の内部。戻り: (t, c)。"""
    n = len(c)
    l_ = min(int(np.searchsorted(t, u, side="right")) - 1, n - 1)
    if l_ < k:
        raise ValueError(f"_insert_knot: u = {u!r} が基本区間 [{t[k]!r}, {t[n]!r}] の内部でない")
    cn = np.empty(n + 1)
    cn[:l_ - k + 1] = c[:l_ - k + 1]
    i = np.arange(l_ - k + 1, l_ + 1)
    al = (u - t[i]) / (t[i + k] - t[i])
    cn[i] = (1.0 - al) * c[i - 1] + al * c[i]
    cn[l_ + 1:] = c[l_:]
    return np.insert(t, l_ + 1, u), cn


def _restrict_bspline(t, c, k: int, lo: float, hi: float) -> tuple:
    """B-spline を [lo, hi] に切り出す: lo・hi を重複度 k + 1 まで挿入し (基本区間の端ならそのまま)、その区間の係数とノットを返す。"""
    t, c = np.asarray(t, dtype=float), np.asarray(c, dtype=float)
    n_ins = 0
    for u in (lo, hi):
        if t[k] < u < t[len(t) - k - 1]:
            for _ in range(k + 1 - int(np.sum(t == u))):
                t, c = _insert_knot(t, c, k, u)
                n_ins += 1
        elif not (u == t[k] or u == t[len(t) - k - 1]):
            raise ValueError(f"_restrict_bspline: {u!r} が基本区間 [{t[k]!r}, {t[len(t) - k - 1]!r}] の外")
    i0 = int(np.searchsorted(t, lo, side="left"))
    i1 = int(np.searchsorted(t, hi, side="left"))
    return t[i0:i1 + k + 1].copy(), c[i0:i1].copy(), n_ins


def _bspline_union_add(tA, cA, tB, cB, k: int) -> tuple:
    """同じ基本区間 (両端の重複度 k + 1) の 2 本の B-spline を、内部ノットの和集合 (重複度は大きい方) にノット挿入で揃えて
    係数を足す。戻り: (t, c, A への挿入数, B への挿入数)。"""
    if not (tA[0] == tB[0] and tA[-1] == tB[-1]):
        raise ValueError("_bspline_union_add: 基本区間が違う")
    va, ma = np.unique(tA[k + 1:len(tA) - k - 1], return_counts=True)
    vb, mb = np.unique(tB[k + 1:len(tB) - k - 1], return_counts=True)
    da, db = dict(zip(va.tolist(), ma.tolist())), dict(zip(vb.tolist(), mb.tolist()))
    nA = nB = 0
    for u in np.union1d(va, vb).tolist():
        m = max(da.get(u, 0), db.get(u, 0))
        for _ in range(m - da.get(u, 0)):
            tA, cA = _insert_knot(tA, cA, k, u)
            nA += 1
        for _ in range(m - db.get(u, 0)):
            tB, cB = _insert_knot(tB, cB, k, u)
            nB += 1
    if not np.array_equal(tA, tB):
        raise RuntimeError("_bspline_union_add: 挿入後のノットが一致しない")
    return tA, cA + cB, nA, nB


def _remove_knot(t, c, k: int, u: float) -> tuple:
    """ノット u を 1 個除去する (NURBS Book A5.8 の 1 回分)。除去後の係数を左右から解き、余った 1 本の式の食い違い (重複度 k + 1
    のときは隣り合う 2 係数の差 = 値の跳び) を誤差とする。戻り: (t, c, 誤差)。"""
    t, c = np.asarray(t, dtype=float), np.asarray(c, dtype=float)
    r = int(np.searchsorted(t, u, side="right")) - 1                   # 最後の u の位置
    s = int(np.sum(t == u))
    if s < 1 or t[r] != u:
        raise ValueError(f"_remove_knot: {u!r} はノットでない")
    tau = np.delete(t, r)
    if s == k + 1:                                                     # 値の連続 (C⁰) だけを確かめて 2 係数を 1 つに
        j = r - k - 1
        err = abs(c[j + 1] - c[j])
        cn = np.r_[c[:j], 0.5 * (c[j] + c[j + 1]), c[j + 2:]]
        return tau, cn, float(err)
    lo, hi = r - k, r - s                                              # 式の番号 i ∈ [lo, hi]: c_i = αᵢ Q_i + (1 − αᵢ) Q_{i−1}
    Q = np.empty(len(c) - 1)
    Q[:lo] = c[:lo]                                                    # Q_j = c_j (j < lo)
    Q[hi:] = c[hi + 1:]                                                # Q_j = c_{j+1} (j ≥ hi)

    def al(i):
        return (u - tau[i]) / (tau[i + k] - tau[i])
    a_ = (hi - lo) // 2
    for i in range(lo, lo + a_):                                       # 左から
        Q[i] = (c[i] - (1.0 - al(i)) * Q[i - 1]) / al(i)
    for i in range(hi, lo + a_, -1):                                   # 右から
        Q[i - 1] = (c[i] - al(i) * Q[i]) / (1.0 - al(i))
    i = lo + a_
    err = abs(c[i] - (al(i) * Q[i] + (1.0 - al(i)) * Q[i - 1]))
    return tau, Q, float(err)


_REQ_COMMON = ("x_in", "x_e", "x_throat", "r_throat", "kappa_throat", "offset_mode", "_dstar_hist", "design", "r_U", "L_U",
               "_herm_x0", "_throat_diag", "_delta_r_applied", "pw_upstream")
REQUIRED_ATTRS_BY_UPSTREAM = {"ramp": ("_ramp", "ramp_gate"), "poly": ("upstream_gate", "upstream_poly")}
_REQ_SINGLE = ("spline", "joints", "fit_diag")


class SingleBSplinePhysicalWall:
    r"""物理壁を入口から出口まで **1 本の x の 5 次 B-spline** で表した壁 (plan tooling-nozzle-wall-single-bspline §4.2、
    係数の求め方は plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1b)。

    元の物理壁は `PhysicalNozzleWall` の解析経路 (joint 設計壁、offset radial) の **`pw_upstream: poly`** に限る: 直管 r_U・
    5 次多項式 Q・S + $\delta_r$ (S も $\delta_r$ も x の 5 次 B-spline) で、全区間が x の 5 次の区分多項式なので、当てはめでなく
    **ノット挿入で代数的に**組む (最小二乗の版は外した — ユーザ決定 2026-10-07)。`ramp` の壁 (ランプ区間の $s\,\delta_r$ は区間ごとに
    10 次で厳密に表せない) は例外。

    1. 直管は定数の係数 (Bézier)、Q は 5 次多項式の Bernstein 係数 ([−L_U, 0] の Bézier)。
    2. [0, x_e] は S と $\delta_r$ (`delta_r_from_table` の補間スプライン) を、$\delta_r$ は x = 0 と x_e で重複度 6 まで挿入して切り出し、
       内部ノットの和集合にノット挿入で揃えて係数を足す。$\delta_r$ の表が [0, x_e] を覆わなければ例外。
    3. 区間ごとの表現を重複度 6 で並べたものから、継ぎ目 (−L_U・0) のノットを 3 回ずつ除去して重複度 3 ($C^2$) にする
       (除去の誤差を記録し、`KNOT_REMOVAL_TOL` を超えたら例外)。
    4. 元の壁との差を各ノット区間の密な点と区間多項式の極値で測り (`bspline_wall_errors`)、`SINGLE_BSPLINE_TOL` を超えたら例外。
    5. スロート量は 1 本の B-spline から大域最小 (`wall_global_min`、ノット区間ごと) で求め直す。上流のゲート (一意な最小・前後の
       単調性・|r″ − H″|・継ぎ目の跳び) も 1 本で評価し直し、不合格なら例外。

    下流 (メッシュ・初期値・NS の準備・報告) には `PhysicalNozzleWall` と同じ属性を出す (`REQUIRED_ATTRS`、欠けたら
    `check_required_attrs` が例外)。`r(x, deriv)` は定義域の外を外挿しない: 座標の float32 丸め分 (`_dom_tol`) だけ外れた点は
    端で評価し、それより外は例外。"""

    REQUIRED_ATTRS = _REQ_COMMON + REQUIRED_ATTRS_BY_UPSTREAM["poly"] + _REQ_SINGLE

    def __init__(self, src, tol: dict | None = None) -> None:
        import time
        from scipy.interpolate import BSpline
        t_start = time.perf_counter()
        tol = dict(SINGLE_BSPLINE_TOL if tol is None else tol)
        if not isinstance(src, PhysicalNozzleWall) or not src.analytic or src.offset_mode != "radial":
            raise ValueError("SingleBSplinePhysicalWall: joint 壁 + 物理壁の解析経路 (offset: radial) の PhysicalNozzleWall 専用 "
                             f"(受け取った: {type(src).__name__}, analytic={getattr(src, 'analytic', '?')}, "
                             f"offset={getattr(src, 'offset_mode', '?')})")
        if src.pw_upstream != "poly":
            raise ValueError(f"SingleBSplinePhysicalWall: pw_upstream {src.pw_upstream!r} は不可 — ランプ区間の s·δ_r は区間ごとに 10 次で "
                             "5 次の B-spline で厳密に表せない (pw_upstream: poly の壁だけ; plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1b)")
        design = src.design
        if not isinstance(design, JointFitCFDWall):
            raise ValueError(f"SingleBSplinePhysicalWall: 設計壁が JointFitCFDWall でない ({type(design).__name__})")
        if float(design.x0) != 0.0:
            raise ValueError(f"SingleBSplinePhysicalWall: 設計壁の始点 x0 = {design.x0!r} ≠ 0 (旧・縦 starting line 構成) は未対応")
        dr = src._dr
        for a_ in ("spline", "x_range"):
            if not hasattr(dr, a_):
                raise ValueError(f"SingleBSplinePhysicalWall: δ_r の関数に `{a_}` が無い (delta_r_from_table の関数を使う)")
        k = 5
        S, D = design._spl, dr.spline
        if int(S.k) != k or int(D.k) != k:
            raise ValueError(f"SingleBSplinePhysicalWall: 設計壁 S (k={S.k})・δ_r (k={D.k}) は 5 次の B-spline であること")
        L = float(src.L_U)
        x_in, x_e = float(src.x_in), float(src.x_e)
        tab_lo, tab_hi = (float(v) for v in dr.x_range)
        if not (tab_lo <= 0.0 and tab_hi >= x_e):
            raise ValueError(f"SingleBSplinePhysicalWall: δ_r の表の範囲 [{tab_lo:.9g}, {tab_hi:.9g}] が設計スロート 0 から出口 "
                             f"{x_e:.9g} までを覆わない (表の外は端値延長で導関数 0 — 1 本の 5 次 B-spline で表せない)")
        if not (x_in < -L < 0.0 < x_e):
            raise ValueError(f"SingleBSplinePhysicalWall: 区間の並び x_in {x_in!r} < −L_U {-L!r} < 0 < x_e {x_e!r} でない")
        tS = np.asarray(S.t, dtype=float)
        cS = np.asarray(S.c, dtype=float)[:len(tS) - k - 1]
        if not (np.all(tS[:k + 1] == 0.0) and np.all(tS[len(tS) - k - 1:] == x_e)):
            raise ValueError("SingleBSplinePhysicalWall: 設計壁 S のノットが [0, x_e] で閉じていない (端の重複度 6)")
        tD0 = np.asarray(D.t, dtype=float)
        # 1) [0, x_e]: δ_r を切り出し、S と内部ノットの和集合に揃えて係数を足す
        tD, cD, n_cut = _restrict_bspline(tD0, np.asarray(D.c, dtype=float)[:len(tD0) - k - 1], k, 0.0, x_e)
        tR, cR, nS_ins, nD_ins = _bspline_union_add(tS, cS, tD, cD, k)
        # 2) 直管 (定数)・Q (Bézier)・[0, x_e] を重複度 6 で並べる
        bq = _power_to_bernstein(src._q_c)
        t = np.r_[[x_in] * (k + 1), [-L] * (k + 1), [0.0] * (k + 1), tR[k + 1:]]
        c = np.r_[[float(src.r_U)] * (k + 1), bq, cR]
        # 3) 継ぎ目のノットを除去して重複度 3 (C²) に
        removal = []
        for u in (-L, 0.0):
            errs = []
            for _ in range(k + 1 - 3):
                t, c, e = _remove_knot(t, c, k, u)
                errs.append(e)
            removal.append({"x": u, "errors_rt": errs, "max_rt": float(max(errs)), "note": "1 回目 = 値、2 回目 = 1 階、3 回目 = 2 階の連続"})
        rem_max = max(r_["max_rt"] for r_ in removal)
        if not rem_max <= KNOT_REMOVAL_TOL:
            raise ValueError(f"SingleBSplinePhysicalWall: 継ぎ目のノット除去の誤差 {rem_max:.3e} r_t > {KNOT_REMOVAL_TOL:g} "
                             "(継ぎ目で 2 階微分まで連続でない)")
        self.spline = BSpline(t, c, k)
        self.design = design
        self.r_U, self.L_U = float(src.r_U), L
        self.x_in, self.x_e = x_in, x_e
        self._herm_x0 = -L
        self.offset_mode = "radial"
        self.pw_upstream = "poly"
        self.pw_upstream_source = src.pw_upstream_source
        self.upstream_poly = src.upstream_poly
        self._dstar_hist = src._dstar_hist                    # 相関 δ* (prepare_info の診断。壁には使わない)
        self.joints = [-L, 0.0]
        self._dom_tol = 4.0 * float(np.finfo(np.float32).eps) * max(abs(x_in), abs(x_e))
        self._delta_r_applied = lambda x: self.r(np.asarray(x, dtype=float)) - self.design.r(np.asarray(x, dtype=float))
        dist = np.unique(t)
        iv = np.c_[dist[:-1], dist[1:]]
        # 4) 元の壁 (区分表現) との差 (密な点と区間多項式の極値)
        err = bspline_wall_errors(self.spline, src.r, iv)
        brk = [x_in, -L, 0.0, x_e]
        seg_rows = []
        for a, b, name in zip(brk[:-1], brk[1:], ("pipe", "Q", "S + delta_r")):
            m = (iv[:, 0] >= a) & (iv[:, 1] <= b)
            row = {"range": [a, b], "parts": name, "n_intervals": int(m.sum())}
            for n in range(4):
                e, xe_ = err[n]
                j = int(np.argmax(np.where(m, e, -1.0)))
                row[f"max_d{n}"], row[f"x_at_max_d{n}"] = float(e[j]), float(xe_[j])
            seg_rows.append(row)
        emax = {n: float(err[n][0].max()) for n in range(4)}
        lim = {0: float(tol["r"]), 1: float(tol["r1"]), 2: float(tol["r2"])}
        bad = [f"{('r', 'r′', 'r″')[n]} {emax[n]:.3e} > {lim[n]:g}" for n in (0, 1, 2) if not emax[n] <= lim[n]]
        if bad:
            raise ValueError("SingleBSplinePhysicalWall: 元の物理壁との差が許容誤差を超えた — " + "; ".join(bad))
        # 継ぎ目の値・1〜3 階微分の跳びと、元の壁の左右の極限との差
        jl, jr = bspline_piece_limits(self.spline, self.joints, nmax=3)
        src_lim = {-L: ([self.r_U, 0.0, 0.0, 0.0], [float(src._q_eval(src._q_c, -L, 0.0, np.array([-L]), n)[0]) for n in range(4)]),
                   0.0: ([float(src._q_eval(src._q_c, -L, 0.0, np.array([0.0]), n)[0]) for n in range(4)],
                         [float(design.r(np.array([0.0]), n)[0] + dr(np.array([0.0]), n)[0]) for n in range(4)])}
        joint_jumps = []
        for q, xj in enumerate(self.joints):
            sl, sr = src_lim[xj]
            joint_jumps.append({"x": xj, **{f"jump_d{n}": float(abs(jr[n, q] - jl[n, q])) for n in range(4)},
                                **{f"left_vs_source_d{n}": float(abs(jl[n, q] - sl[n])) for n in range(3)},
                                **{f"right_vs_source_d{n}": float(abs(jr[n, q] - sr[n])) for n in range(3)}})
        # 5) スロート (1 本の B-spline の大域最小) と上流のゲート
        srch = wall_global_min(lambda x, n: self.spline(np.asarray(x, dtype=float), n), dist)
        self.x_throat, self.r_throat, self.kappa_throat = srch["x_t"], srch["r_t"], srch["kappa_t"]
        throat_diff = {"x": self.x_throat - float(src.x_throat), "r": self.r_throat - float(src.r_throat),
                       "kappa": self.kappa_throat - float(src.kappa_throat)}
        d2max, x_d2 = _poly_absmax(lambda x: self.spline(x, 2) - design.up.r(x, 2), -L, 0.0, 3)
        seams = [{"x": xj, **{f"jump_d{n}": float(abs(jr[n, q] - jl[n, q])) for n in range(3)}} for q, xj in enumerate(self.joints)]
        seam_max = max(s_[f"jump_d{n}"] for s_ in seams for n in range(3))
        self.upstream_gate = _poly_gate_record(srch, d2max, x_d2, float(src.upstream_gate["limit_d2"]), seams, seam_max,
                                               source=self.pw_upstream_source, repr_="single_bspline")
        self.upstream_gate["source_wall"] = {k_: src.upstream_gate[k_] for k_ in
                                             ("pass", "max_abs_d2_change_vs_H", "max_seam_jump", "throat", "max_r1_before", "min_r1_after")}
        if not self.upstream_gate["pass"]:
            raise ValueError("SingleBSplinePhysicalWall: 1 本の B-spline で上流のゲートが不合格 — " + _poly_gate_message(self.upstream_gate))
        self._throat_diag = {"method": "single_bspline_global_min", "search": srch, "diff_vs_source": throat_diff}
        dk = np.diff(dist)
        self.fit_diag = {"kind": "single_bspline", "construction": "knot_insertion", "k": k, "n_coef": int(len(c)),
                         "n_knots": int(len(t)), "n_distinct_knots": int(len(dist)), "min_knot_gap_rt": float(dk.min()),
                         "x_min_knot_gap": float(dist[int(np.argmin(dk))]), "joints": self.joints,
                         "domain": [x_in, x_e], "delta_r_table_range": [tab_lo, tab_hi],
                         "insertion": {"delta_r_cut_at": [0.0, x_e], "n_inserted_to_cut_delta_r": int(n_cut),
                                       "n_inserted_into_S": int(nS_ins), "n_inserted_into_delta_r": int(nD_ins)},
                         "knot_removal": removal, "knot_removal_tol": KNOT_REMOVAL_TOL,
                         "tol": {"r": lim[0], "r1": lim[1], "r2": lim[2]},
                         "max_err": {f"d{n}": emax[n] for n in range(4)}, "segments": seg_rows,
                         "joint_jumps": joint_jumps, "throat": {"x": self.x_throat, "r": self.r_throat, "kappa": self.kappa_throat},
                         "throat_source": {"x": float(src.x_throat), "r": float(src.r_throat), "kappa": float(src.kappa_throat)},
                         "throat_diff": throat_diff, "build_seconds": float(time.perf_counter() - t_start)}
        check_required_attrs(self)

    def r(self, x, deriv: int = 0):
        x = np.asarray(x, dtype=float)
        if deriv not in (0, 1, 2, 3, 4):
            raise ValueError("deriv は 0..4")
        if x.size and (float(np.min(x)) < self.x_in - self._dom_tol or float(np.max(x)) > self.x_e + self._dom_tol):
            raise ValueError(f"SingleBSplinePhysicalWall.r: x ∈ [{float(np.min(x)):.9g}, {float(np.max(x)):.9g}] が定義域 "
                             f"[{self.x_in:.9g}, {self.x_e:.9g}] の外 (外挿しない)")
        xc = np.clip(x, self.x_in, self.x_e)
        return self.spline(xc, deriv) if deriv else self.spline(xc)

    def theta(self, x):
        return np.arctan(self.r(x, 1))

    def kappa(self, x):
        rp, rpp = self.r(x, 1), self.r(x, 2)
        return rpp / (1.0 + rp * rp) ** 1.5

    def validate(self, n: int = 4000) -> list:
        """`PhysicalNozzleWall.validate` と同じ検査 (半径 > 0・最小半径 = r_throat・スロート下流の単調・上流の単調収縮・
        継ぎ目の C¹/C²) + 必須属性。"""
        msgs = []
        try:
            check_required_attrs(self)
        except ValueError as e:
            msgs.append(str(e))
        xs = np.linspace(self.x_in, self.x_e, n)
        rv = self.r(xs)
        if np.any(rv <= 0.0):
            msgs.append("壁半径が非正")
        if abs(float(rv.min()) - self.r_throat) > 5e-3:
            msgs.append(f"最小半径 {float(rv.min()):.4f} != r_throat {self.r_throat:.4f}")
        if np.any(np.diff(rv[xs >= self.x_throat]) < -1e-9):
            msgs.append("スロート下流で半径が非単調")
        if np.any(np.diff(rv[(xs >= self._herm_x0) & (xs <= self.x_throat)]) > 1e-9):
            msgs.append("上流 (多項式 Q) が非単調収縮")
        h = 1e-6
        for xc in self.joints:
            dl = float(self.r(np.array([xc - h]), 1)[0])
            dr_ = float(self.r(np.array([xc + h]), 1)[0])
            if abs(dl - dr_) > 5e-3:
                msgs.append(f"継ぎ目 x={xc:g} の接線不連続 ({dl:.5f} vs {dr_:.5f})")
            cl = float(self.r(np.array([xc - h]), 2)[0])
            cr = float(self.r(np.array([xc + h]), 2)[0])
            if abs(cl - cr) > 0.05 * max(abs(cl), abs(cr), 1.0):
                msgs.append(f"継ぎ目 x={xc:g} の曲率不連続 ({cl:.4f} vs {cr:.4f})")
        return msgs


def required_attrs(wall) -> tuple:
    """物理壁の必須属性 (上流の方式別、plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1): 共通 + `ramp` なら `_ramp`・
    `ramp_gate`、`poly` なら `upstream_gate`・`upstream_poly` + 1 本の B-spline なら `spline`・`joints`・`fit_diag`。
    `pw_upstream` が無い・不正なら ValueError (方式が分からないまま既定に落とさない)。"""
    mode = wall.__dict__.get("pw_upstream") if hasattr(wall, "__dict__") else None
    if mode not in PW_UPSTREAMS:
        raise ValueError(f"物理壁 ({type(wall).__name__}) の pw_upstream が無い・不正 ({mode!r})")
    req = _REQ_COMMON + REQUIRED_ATTRS_BY_UPSTREAM[mode]
    if isinstance(wall, SingleBSplinePhysicalWall):
        req = req + _REQ_SINGLE
    return req


def check_required_attrs(wall) -> None:
    """`required_attrs` (下流が読む属性と診断属性、方式別) が全部あり None でないこと。欠けたら ValueError
    (`getattr(..., None)` で黙って既定に落ちる読み出しがあるため、構築時と validate で必須として検査する)。"""
    req = required_attrs(wall)
    miss = [a for a in req if getattr(wall, a, None) is None]
    if miss:
        raise ValueError(f"物理壁 ({type(wall).__name__}, pw_upstream {wall.pw_upstream}) に必須属性が無い・None: {miss}")


def _delta_r_eval(spl, lo: float, hi: float):
    """δ_r の評価関数 (`runner_axismach.delta_r_from_table` と同じ規約: 表の範囲外は値を端値クリップ・導関数 0)。
    保存した δ_r の B-spline から legacy の物理壁を復元するのに使う。"""
    def f(xq, deriv: int = 0, _s=spl, _lo=lo, _hi=hi):
        xq = np.asarray(xq, dtype=float)
        if deriv == 0:
            return _s(np.clip(xq, _lo, _hi))
        if deriv not in (1, 2, 3):
            raise ValueError("deriv は 0..3")
        return np.where((xq >= _lo) & (xq <= _hi), _s(np.clip(xq, _lo, _hi), deriv), 0.0)
    f.supports_deriv = True
    f.spline = spl
    f.x_range = (float(lo), float(hi))
    return f


def _spline_rec(spl, lo: float, hi: float) -> dict:
    return {"k": int(spl.k), "t": [float(v) for v in spl.t], "c": [float(v) for v in np.asarray(spl.c, dtype=float).ravel()],
            "domain": [float(lo), float(hi)]}


def _gate_summary(g: dict) -> dict:
    """上流のゲートの記録から壁ファイルに写す要約 (探索の根の一覧などの大きいものを除く)。"""
    return {k: g[k] for k in ("pw_upstream", "source", "repr", "pass", "throat", "positive_radius", "unique_min", "monotone_before",
                              "max_r1_before", "monotone_after", "min_r1_after", "max_abs_d2_change_vs_H", "x_at_max_d2_change",
                              "limit_d2", "max_seam_jump", "limit_seam") if k in g}


def wall_file_record(wall, scale_m: float, kind: str) -> dict:
    """保存する壁ファイルの中身 (plan tooling-nozzle-wall-single-bspline §4.2「保存した壁の復元の取り決め」、版 2 は plan
    tooling-nozzle-upstream-poly-and-throat-sizing §4.1)。
    kind = 'single_bspline' (wall は `SingleBSplinePhysicalWall`、上流は poly) / 'legacy' (wall は解析経路の `PhysicalNozzleWall` の区分表現)。
    長さは r_t 単位 (x_m = scale_m·x)。設計壁は直管 (r_U・x_in)・上流 Hermite (係数・区間 [−L_U, 0])・S (ノット・係数・
    有効域 [0, x_e]) を持つ。上流の方式 `pw_upstream` を最上位と物理壁に書く。物理壁は 1 本の B-spline か、legacy なら復元に要る全要素
    (ramp: ランプと δ_r の B-spline と表の範囲 / poly: 多項式 Q の係数・基底・区間と δ_r の B-spline と表の範囲)。"""
    if kind not in PHYSICAL_WALL_REPRS:
        raise ValueError(f"wall_file_record: kind={kind!r} は {PHYSICAL_WALL_REPRS} のどれか")
    if kind == "single_bspline" and not isinstance(wall, SingleBSplinePhysicalWall):
        raise ValueError("wall_file_record: single_bspline には SingleBSplinePhysicalWall が要る")
    if kind == "legacy" and not (isinstance(wall, PhysicalNozzleWall) and wall.analytic and wall.offset_mode == "radial"):
        raise ValueError("wall_file_record: legacy の壁ファイルは joint 壁 + 解析経路 (offset radial) の PhysicalNozzleWall だけ")
    mode = wall.pw_upstream
    if mode not in PW_UPSTREAMS:
        raise ValueError(f"wall_file_record: 物理壁の pw_upstream {mode!r} が {PW_UPSTREAMS} のどれでもない")
    check_required_attrs(wall)
    dw = wall.design
    if not isinstance(dw, JointFitCFDWall) or float(dw.x0) != 0.0:
        raise ValueError("wall_file_record: 設計壁は始点 x0 = 0 の JointFitCFDWall に限る")
    S = float(scale_m)
    design = {"kind": "joint_axismach", "R": float(dw.R), "x_in": float(dw.x_in), "x_e": float(dw.x_e), "x0": float(dw.x0),
              "pipe": {"r_U": float(dw.up.r_U), "domain": [float(dw.x_in), -float(dw.up.L_U)]},
              "upstream_hermite": {"L_U": float(dw.up.L_U), "R_t": float(dw.up.R_t), "r_t": float(dw.up.r_t),
                                   "coef": [float(v) for v in dw.up._c], "basis": "Σ coef[j] ξ^j, ξ = (x + L_U)/L_U",
                                   "domain": [-float(dw.up.L_U), 0.0]},
              "S": _spline_rec(dw._spl, float(dw.x0), float(dw.x_e))}
    throat = {"x": float(wall.x_throat), "r": float(wall.r_throat), "kappa": float(wall.kappa_throat)}
    if kind == "single_bspline":
        phys = {"kind": "single_bspline", "pw_upstream": mode, **_spline_rec(wall.spline, wall.x_in, wall.x_e),
                "joints": list(wall.joints), "throat": throat, "n_coef": wall.fit_diag["n_coef"],
                "n_distinct_knots": wall.fit_diag["n_distinct_knots"], "min_knot_gap_rt": wall.fit_diag["min_knot_gap_rt"],
                "max_err_vs_legacy": wall.fit_diag["max_err"], "tol": wall.fit_diag["tol"],
                "construction": wall.fit_diag["construction"], "knot_removal": wall.fit_diag["knot_removal"],
                "upstream_poly": wall.upstream_poly, "upstream_gate": _gate_summary(wall.upstream_gate)}
    else:
        lo_t, hi_t = (float(v) for v in wall._dr.x_range)
        drec = {**_spline_rec(wall._dr.spline, lo_t, hi_t), "x_range": [lo_t, hi_t],
                "outside": "値は端値クリップ・導関数 0 (delta_r_from_table と同じ)"}
        if mode == "ramp":
            phys = {"kind": "legacy", "pw_upstream": "ramp", "ramp": [float(v) for v in wall._ramp], "delta_r": drec,
                    "throat": throat, "formula": "r = r_design + s(x)·δ_r(x), s = 5 次 smoothstep on ramp"}
        else:
            phys = {"kind": "legacy", "pw_upstream": "poly", "upstream_poly": wall.upstream_poly, "delta_r": drec, "throat": throat,
                    "formula": "x < −L_U: r_U (直管) / −L_U ≤ x < 0: Q(x) (upstream_poly) / 0 ≤ x ≤ x_e: S(x) + δ_r(x)",
                    "upstream_gate": _gate_summary(wall.upstream_gate)}
    return {"format": WALL_FILE_FORMAT, "version": WALL_FILE_VERSION, "physical_wall_repr": kind, "pw_upstream": mode,
            "units": {"length": "r_t", "scale_m": S, "to_m": "x_m = scale_m · x, r_m = scale_m · r (n 階微分は scale_m^(1−n) 倍)"},
            "origin": "設計スロート (設計壁 S の始点 x = 0)。物理スロートは physical_wall.throat.x",
            "domain": [float(wall.x_in), float(wall.x_e)], "domain_m": [float(wall.x_in) * S, float(wall.x_e) * S],
            "design_wall": design, "physical_wall": phys}


def save_wall_file(run_dir, wall, scale_m: float, kind: str) -> tuple:
    """壁ファイル `WALL_FILE` を run_dir に書き、書いたものを読み直して復元できることを確かめる。戻り: (path, record, sha256)。"""
    import hashlib
    import json
    from pathlib import Path
    rec = wall_file_record(wall, scale_m, kind)
    path = Path(run_dir) / WALL_FILE
    txt = json.dumps(rec, indent=1, ensure_ascii=False)
    path.write_text(txt)
    back = load_wall_file(path)
    xs = np.linspace(rec["domain"][0], rec["domain"][1], 2001)
    if not np.array_equal(back["physical"].r(xs), wall.r(xs)):
        raise RuntimeError(f"save_wall_file: 書いた {path} から復元した物理壁が元と一致しない")
    return path, rec, hashlib.sha256(txt.encode()).hexdigest()


class _SavedBSpline:
    """保存した B-spline (k, t, c) を有効域 [lo, hi] でだけ評価する (域外は例外; scipy の BSpline の既定の外挿をしない)。"""

    def __init__(self, rec: dict, name: str) -> None:
        from scipy.interpolate import BSpline
        k, t, c = int(rec["k"]), np.asarray(rec["t"], dtype=float), np.asarray(rec["c"], dtype=float)
        lo, hi = (float(v) for v in rec["domain"])
        if len(t) != len(c) + k + 1:
            raise ValueError(f"壁ファイル: {name} のノット数 {len(t)} ≠ 係数 {len(c)} + 次数 {k} + 1")
        if np.any(np.diff(t) < 0.0) or not (np.all(np.isfinite(t)) and np.all(np.isfinite(c))):
            raise ValueError(f"壁ファイル: {name} のノットが非減少でない・非有限値がある")
        if not (t[k] <= lo < hi <= t[len(t) - k - 1]):
            raise ValueError(f"壁ファイル: {name} の有効域 [{lo}, {hi}] がノットの基本区間 [{t[k]}, {t[len(t) - k - 1]}] の外")
        self.spl, self.lo, self.hi, self.name = BSpline(t, c, k), lo, hi, name

    def __call__(self, x, deriv: int = 0):
        x = np.asarray(x, dtype=float)
        if x.size and (float(np.min(x)) < self.lo or float(np.max(x)) > self.hi):
            raise ValueError(f"{self.name}: x ∈ [{float(np.min(x)):.9g}, {float(np.max(x)):.9g}] が有効域 "
                             f"[{self.lo:.9g}, {self.hi:.9g}] の外 (外挿しない)")
        return self.spl(x, deriv) if deriv else self.spl(x)


class RestoredDesignWall:
    """壁ファイルから復元した設計壁 (直管 + 上流 Hermite + S)。`AxisMachCFDWall.r` と同じ式で評価し、定義域の外は例外。"""

    def __init__(self, rec: dict) -> None:
        from .wall_walldriven import UpstreamThroatPoly
        self.x_in, self.x_e, self.x0 = float(rec["x_in"]), float(rec["x_e"]), float(rec["x0"])
        if self.x0 != 0.0:
            raise ValueError(f"壁ファイル: 設計壁の始点 x0 = {self.x0} ≠ 0 は未対応")
        uh = rec["upstream_hermite"]
        self.up = UpstreamThroatPoly(r_U=float(rec["pipe"]["r_U"]), R_t=float(uh["R_t"]), L_U=float(uh["L_U"]), r_t=float(uh["r_t"]))
        if not np.array_equal(np.asarray(self.up._c, dtype=float), np.asarray(uh["coef"], dtype=float)):
            raise ValueError("壁ファイル: 上流 Hermite の係数が (r_U, R_t, L_U) から作り直した係数と一致しない")
        self.S = _SavedBSpline(rec["S"], "設計壁 S")
        if (self.S.lo, self.S.hi) != (self.x0, self.x_e):
            raise ValueError(f"壁ファイル: 設計壁 S の有効域 {[self.S.lo, self.S.hi]} ≠ [x0, x_e] = {[self.x0, self.x_e]}")

    def r(self, x, deriv: int = 0):
        x = np.asarray(x, dtype=float)
        if deriv not in (0, 1, 2, 3):
            raise ValueError("deriv は 0..3")
        if x.size and (float(np.min(x)) < self.x_in or float(np.max(x)) > self.x_e):
            raise ValueError(f"復元した設計壁: x ∈ [{float(np.min(x)):.9g}, {float(np.max(x)):.9g}] が定義域 "
                             f"[{self.x_in:.9g}, {self.x_e:.9g}] の外 (外挿しない)")
        out = np.empty_like(x)
        xU = -self.up.L_U
        m_pipe, m_up, m_S = x < xU, (x >= xU) & (x < 0.0), x >= 0.0
        out[m_pipe] = self.up.r_U if deriv == 0 else 0.0
        if m_up.any():
            out[m_up] = self.up.r(x[m_up], deriv)
        if m_S.any():
            out[m_S] = self.S(x[m_S], deriv)
        return out


def _restored_delta_r(d: dict):
    """壁ファイルの δ_r の記録 → `delta_r_from_table` と同じ規約の評価関数。"""
    from scipy.interpolate import BSpline
    lo, hi = (float(v) for v in d["x_range"])
    t, c, k = np.asarray(d["t"], dtype=float), np.asarray(d["c"], dtype=float), int(d["k"])
    if len(t) != len(c) + k + 1:
        raise ValueError("壁ファイル: δ_r のノット数と係数の数が合わない")
    return _delta_r_eval(BSpline(t, c, k), lo, hi)


class _RestoredLegacyPhysicalWall:
    """壁ファイル (legacy、pw_upstream ramp・版 1 を含む) から復元した物理壁 r = r_design + s·δ_r
    (`PhysicalNozzleWall` の解析経路の ramp と同じ式・同じ関数)。"""
    _s = PhysicalNozzleWall._s
    _sdr = PhysicalNozzleWall._sdr

    def __init__(self, design: RestoredDesignWall, rec: dict) -> None:
        self.design = design
        self.pw_upstream = "ramp"
        self._ramp = tuple(float(v) for v in rec["ramp"])
        self._dr = _restored_delta_r(rec["delta_r"])
        self.x_in, self.x_e = design.x_in, design.x_e

    def r(self, x, deriv: int = 0):
        x = np.asarray(x, dtype=float)
        if x.size and (float(np.min(x)) < self.x_in or float(np.max(x)) > self.x_e):
            raise ValueError(f"復元した物理壁 (legacy): x が定義域 [{self.x_in:.9g}, {self.x_e:.9g}] の外 (外挿しない)")
        return self.design.r(x, deriv) + self._sdr(x, deriv)


class _RestoredPolyPhysicalWall:
    """壁ファイル (legacy、pw_upstream poly) から復元した物理壁: 直管 r_U / 多項式 Q / S + δ_r (`PhysicalNozzleWall` の poly と
    同じ式・同じ関数)。Q の端条件 (x = −L_U で (r_U, 0, 0)、x = 0 で S + δ_r) が `POLY_SEAM_TOL` で合わなければ例外 (壊れた壁ファイル)。"""

    def __init__(self, design: RestoredDesignWall, rec: dict) -> None:
        from .wall import _poly_eval
        self.design = design
        self.pw_upstream = "poly"
        up = rec["upstream_poly"]
        miss = [k_ for k_ in ("coef", "domain") if not isinstance(up, dict) or up.get(k_) is None]
        if miss:
            raise ValueError(f"壁ファイル: physical_wall.upstream_poly に必要な要素が無い: {miss}")
        self.L_U, self.r_U = float(design.up.L_U), float(design.up.r_U)
        self._q_c = np.asarray(up["coef"], dtype=float)
        if self._q_c.shape != (6,) or not np.all(np.isfinite(self._q_c)):
            raise ValueError(f"壁ファイル: upstream_poly.coef は有限の 6 個 (受け取った形 {self._q_c.shape})")
        if [float(v) for v in up["domain"]] != [-self.L_U, 0.0]:
            raise ValueError(f"壁ファイル: upstream_poly.domain {up['domain']} ≠ [−L_U, 0] = {[-self.L_U, 0.0]}")
        self._q_eval = _poly_eval
        self._dr = _restored_delta_r(rec["delta_r"])
        self.x_in, self.x_e = design.x_in, design.x_e
        L = self.L_U
        qL = [float(_poly_eval(self._q_c, -L, 0.0, np.array([-L]), n)[0]) for n in range(3)]
        q0 = [float(_poly_eval(self._q_c, -L, 0.0, np.array([0.0]), n)[0]) for n in range(3)]
        e0 = [float(design.r(np.array([0.0]), n)[0] + self._dr(np.array([0.0]), n)[0]) for n in range(3)]
        jump = max([abs(qL[0] - self.r_U), abs(qL[1]), abs(qL[2])] + [abs(q0[n] - e0[n]) for n in range(3)])
        if not jump <= POLY_SEAM_TOL:
            raise ValueError(f"壁ファイル: 多項式 Q の端条件が直管・下流の物理壁と合わない (継ぎ目の跳び {jump:.2e} > {POLY_SEAM_TOL:g})")

    def r(self, x, deriv: int = 0):
        x = np.asarray(x, dtype=float)
        if x.size and (float(np.min(x)) < self.x_in or float(np.max(x)) > self.x_e):
            raise ValueError(f"復元した物理壁 (legacy, poly): x が定義域 [{self.x_in:.9g}, {self.x_e:.9g}] の外 (外挿しない)")
        if deriv not in (0, 1, 2, 3):
            raise ValueError("deriv は 0..3")
        out = np.empty_like(x)
        m_pipe = x < -self.L_U
        m_q = (x >= -self.L_U) & (x < 0.0)
        m_dn = x >= 0.0
        out[m_pipe] = self.r_U if deriv == 0 else 0.0
        if m_q.any():
            out[m_q] = self._q_eval(self._q_c, -self.L_U, 0.0, x[m_q], deriv)
        if m_dn.any():
            xd = x[m_dn]
            out[m_dn] = self.design.r(xd, deriv) + self._dr(xd, deriv)
        return out


class _RestoredSingleBSplineWall:
    """壁ファイル (single_bspline) から復元した物理壁 (1 本の B-spline。定義域の外は例外)。"""

    def __init__(self, b: _SavedBSpline, pw_upstream: str) -> None:
        self._b, self.spline = b, b.spl
        self.x_in, self.x_e = b.lo, b.hi
        self.pw_upstream = pw_upstream

    def r(self, x, deriv: int = 0):
        return self._b(x, deriv)


_WALL_FILE_KEYS = {
    "top": ("format", "version", "physical_wall_repr", "units", "domain", "design_wall", "physical_wall"),
    "design_wall": ("R", "x_in", "x_e", "x0", "pipe", "upstream_hermite", "S"),
    "spline": ("k", "t", "c", "domain"),
    "single_bspline": ("kind", "k", "t", "c", "domain", "joints", "throat"),
    "legacy": ("kind", "ramp", "delta_r", "throat"),
    # 版 2 (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1): 上流の方式別の必須要素
    ("legacy", "ramp"): ("kind", "pw_upstream", "ramp", "delta_r", "throat"),
    ("legacy", "poly"): ("kind", "pw_upstream", "upstream_poly", "delta_r", "throat"),
    ("single_bspline", "poly"): ("kind", "pw_upstream", "k", "t", "c", "domain", "joints", "throat", "upstream_poly"),
}


def load_wall_file(path) -> dict:
    """壁ファイル (`WALL_FILE`、run_dir を渡してもよい) を読み、物理壁と設計壁の評価関数を返す (共通の読み込み関数)。
    形式・版・表現の種類・上流の方式・必要な要素を検査し、欠けていたら ValueError。評価関数は有効域の外で例外 (外挿しない)。
    版 1 (上流の方式の記録が無い) は ramp として読む。版 2 の single_bspline は poly だけ (ramp + single_bspline は作れない)。
    戻り: {"record", "physical" (r(x, deriv)), "design" (RestoredDesignWall), "kind", "pw_upstream", "version", "scale_m",
    "domain", "throat"}。"""
    import json
    from pathlib import Path
    p = Path(path)
    if p.is_dir():
        p = p / WALL_FILE
    rec = json.loads(p.read_text())

    def need(d, keys, where):
        miss = [k_ for k_ in keys if not isinstance(d, dict) or k_ not in d or d[k_] is None]
        if miss:
            raise ValueError(f"壁ファイル {p}: {where} に必要な要素が無い: {miss}")

    need(rec, _WALL_FILE_KEYS["top"], "最上位")
    if rec["format"] != WALL_FILE_FORMAT:
        raise ValueError(f"壁ファイル {p}: format {rec['format']!r} ≠ {WALL_FILE_FORMAT!r}")
    version = rec["version"]
    if isinstance(version, bool) or version not in WALL_FILE_VERSIONS_READ:
        raise ValueError(f"壁ファイル {p}: 版 {version!r} は未対応 (読める版 {WALL_FILE_VERSIONS_READ})")
    kind = rec["physical_wall_repr"]
    if kind not in PHYSICAL_WALL_REPRS:
        raise ValueError(f"壁ファイル {p}: physical_wall_repr {kind!r} は {PHYSICAL_WALL_REPRS} のどれでもない")
    if version == 1:
        mode = "ramp"                                                  # 旧形式: 方式の記録が無い → ランプ
        if "pw_upstream" in rec or "pw_upstream" in (rec.get("physical_wall") or {}):
            raise ValueError(f"壁ファイル {p}: 版 1 に pw_upstream がある (版と中身が食い違う)")
    else:
        need(rec, ("pw_upstream",), "最上位 (版 2)")
        mode = rec["pw_upstream"]
        if not isinstance(mode, str) or mode not in PW_UPSTREAMS:
            raise ValueError(f"壁ファイル {p}: pw_upstream {mode!r} は {PW_UPSTREAMS} のどれでもない")
        if kind == "single_bspline" and mode != "poly":
            raise ValueError(f"壁ファイル {p}: single_bspline は pw_upstream poly だけ (受け取った {mode!r})")
    need(rec["units"], ("length", "scale_m"), "units")
    if rec["units"]["length"] != "r_t":
        raise ValueError(f"壁ファイル {p}: 長さの単位 {rec['units']['length']!r} は未対応 (r_t)")
    dw = rec["design_wall"]
    need(dw, _WALL_FILE_KEYS["design_wall"], "design_wall")
    need(dw["S"], _WALL_FILE_KEYS["spline"], "design_wall.S")
    need(dw["pipe"], ("r_U",), "design_wall.pipe")
    need(dw["upstream_hermite"], ("L_U", "R_t", "r_t", "coef"), "design_wall.upstream_hermite")
    design = RestoredDesignWall(dw)
    ph = rec["physical_wall"]
    need(ph, _WALL_FILE_KEYS[kind] if version == 1 else _WALL_FILE_KEYS[(kind, mode)], "physical_wall")
    if ph["kind"] != kind:
        raise ValueError(f"壁ファイル {p}: physical_wall.kind {ph['kind']!r} ≠ physical_wall_repr {kind!r}")
    if version >= 2 and ph["pw_upstream"] != mode:
        raise ValueError(f"壁ファイル {p}: physical_wall.pw_upstream {ph['pw_upstream']!r} ≠ 最上位の pw_upstream {mode!r}")
    dom = [float(v) for v in rec["domain"]]
    if dom != [design.x_in, design.x_e]:
        raise ValueError(f"壁ファイル {p}: domain {dom} ≠ 設計壁の [x_in, x_e] {[design.x_in, design.x_e]}")
    if kind == "single_bspline":
        phys = _SavedBSpline(ph, "物理壁 (single_bspline)")
        if [phys.lo, phys.hi] != dom:
            raise ValueError(f"壁ファイル {p}: 物理壁の有効域 {[phys.lo, phys.hi]} ≠ domain {dom}")
        physical = _RestoredSingleBSplineWall(phys, mode)
    elif mode == "ramp":
        need(ph["delta_r"], ("k", "t", "c", "x_range"), "physical_wall.delta_r")
        physical = _RestoredLegacyPhysicalWall(design, ph)
    else:
        need(ph["delta_r"], ("k", "t", "c", "x_range"), "physical_wall.delta_r")
        physical = _RestoredPolyPhysicalWall(design, ph)
    need(ph["throat"], ("x", "r", "kappa"), "physical_wall.throat")
    return {"record": rec, "path": str(p), "kind": kind, "pw_upstream": mode, "version": int(version),
            "physical": physical, "design": design,
            "scale_m": float(rec["units"]["scale_m"]), "domain": dom, "throat": dict(ph["throat"])}
