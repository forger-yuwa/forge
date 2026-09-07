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
                 x_offset_lo: float = -0.8, offset: str = "normal", delta_r_x=None) -> None:
        """design_wall: AxisMachCFDWall (非粘性)。wall_tbl: (n,4) MOC 壁 [x,r,θ,M]。
        dstar_x: 任意の δ*(x) [r_t 単位] (None = 上流履歴込み相関)。
        offset: "normal" (旧: 壁法線オフセット + x シフト) / "radial" (半径方向 r_W = r + δ_r(x)、
        plans/active/tooling-nozzle-deltastar-core-matched-euler.md §4.5 の生産経路)。
        delta_r_x: 半径方向補正 δ_r(x) [r_t 単位] を直接与える (与えると offset="radial" 固定)。"""
        from ..evaluate.ic import invert_area_ratio
        from ..feedback.deltastar import dstar_flatplate
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
        for name, xc in (("直管/Hermite", self._herm_x0),
                         ("Hermite/下流 (物理スロート)", self.x_throat)):
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
