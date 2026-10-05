r"""CFD ピン: node Euler 場からスロート特性線 (MOC の初期線) と軸アンカーを取る provider。

計画: plans/active/tooling-nozzle-cfd-pinned-initial-line.md §4.1・§4.2 (diagnostician 2026-10-05 採用版)。
仕様: methods/design/overview.md「初期線の出所: Hall / CFD ピン」。case/45 の試作
(`case/45.isobutane_m6_d155/cfd_initial_line.py`) を生産コードへ移したもので、挙動・入力の拒否は試作と同じ。

`CFDPinnedThroat` は `HallThroat` を継承し `throat_characteristic` / `axis_anchor` / `mach` / `theta` だけを
上書きする (`cd_series` 等は Hall のまま)。`design_chain` は `geometry.initial_line: cfd` でこれを使う。

- 場: 構造格子 (断面は x 一定) の M, θ を (x, η = r/r_w(x)) 平面の 3 次スプラインで補間。壁 r_w(x) は断面の壁節点の 3 次補間。
- 線: 壁足は (0, 1, θ=0) に厳密。追跡は壁足 1e-6 内側から dr/dx = tan(θ−μ) を RK4 (ds 2e-4) で軸を跨ぐまで。
  壁足の M は内側 10 点の線形外挿、軸端 (r=0) の x・M・θ は追跡点の r ≤ 0.05 の 2 次多項式 (偶関数当てはめは使わない —
  C⁻ は軸を斜めに横切る)。
- 軸アンカー: M_A = 線の軸端の M (線と同じ出所)、M′_A = 軸の evenfit (各断面で r>0 の 4 点に M=a₀+a₂r²) の 4 次窓フィット
  (窓 ±0.25)、M″_A は Hall の式を x₀ で評価 (CFD の局所 M″ は law を実測軸から離すので使わない)。
- 凍結源: 同じ形・ガスで Hall 初期線の V0 型壁を Euler で解いた場 (ピン壁自身の Euler は循環定義になるので使わない)。
"""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline

from ..geometry.transonic import HallThroat


def _sha16(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()[:16]


def _check_run_config(rd: Path, info: dict) -> dict:
    """凍結源の実効設定を YAML の構造として照合する (codex result M2: 文字列一致では SST・viscMethod 2 を通していた)。
    node・軸対称・Euler (viscMethod 0 かつ visc 0・乱流なし・物理壁はすべて slip) でなければ拒否。戻り: 照合した値。"""
    import yaml
    cfg = yaml.safe_load((rd / "solverConfig.yaml").read_text()) or {}
    bc = yaml.safe_load((rd / "bcondConfig.yaml").read_text()) or {}
    mesh = cfg.get("mesh") or {}
    pp = cfg.get("physProp") or {}
    turb = cfg.get("turbulence") or {}
    errs = []
    if str(mesh.get("discretization")) != "node":
        errs.append(f"mesh.discretization={mesh.get('discretization')!r} (node でない)")
    if int(mesh.get("isAxisymmetric", 0) or 0) != 1:
        errs.append(f"mesh.isAxisymmetric={mesh.get('isAxisymmetric')!r} (軸対称でない)")
    if int(pp.get("viscMethod", -1)) != 0 or float(pp.get("visc", -1.0)) != 0.0:
        errs.append(f"physProp.viscMethod={pp.get('viscMethod')!r}, visc={pp.get('visc')!r} (Euler = viscMethod 0 かつ visc 0 でない)")
    if str(turb.get("model", "")).lower() != "none":
        errs.append(f"turbulence.model={turb.get('model')!r} (乱流なしでない)")
    walls = {k: v for k, v in bc.items() if isinstance(v, dict) and str(v.get("kind", "")).lower()
             not in ("axis", "symmetry") and not str(v.get("kind", "")).lower().startswith(("inlet", "outlet"))}
    if not walls:
        errs.append("bcondConfig に物理壁が無い")
    for k, v in walls.items():
        if str(v.get("kind")) != "slip":
            errs.append(f"壁 {k}: kind={v.get('kind')!r} (slip でない)")
    if info.get("viscous") not in (None, False, "none", "euler"):
        errs.append(f"prepare_info.viscous={info.get('viscous')!r}")
    if errs:
        raise ValueError("CFD ピン: 凍結源が node・軸対称・Euler (slip) の run でない — " + "; ".join(errs))
    return {"discretization": "node", "isAxisymmetric": 1, "viscMethod": 0, "visc": 0.0, "turbulence": "none",
            "walls": {k: "slip" for k in walls}}


def load_run_field(run_dir, res=None, n_tail: int = 1):
    """run の node 場 (r_t 単位の格子 X, R と M, θ) を読む。res: ファイル名 (str) かその列 (平均する) — **必須**
    (省略すると最新場を黙って選ぶので拒否する; codex result M2)。n_tail は互換のため残すが使わない。
    対応範囲: 構造化 node・軸対称・Euler (slip) の run だけ (`_check_run_config`)。不適合は拒否する。
    戻り: X, R, M, θ, source (run・res・照合した設定・凍結源の形と熱力学条件・ハッシュ)。"""
    import h5py
    rd = Path(run_dir)
    if not res:
        raise ValueError("CFD ピン: 凍結源の snapshot (initial_line_res) の明示が必須 (省略すると最新場を選んでしまう)")
    info = json.loads((rd / "prepare_info.json").read_text())
    cfg_checked = _check_run_config(rd, info)
    il = info.get("initial_line") or {}
    if il.get("source", "hall") != "hall":
        raise ValueError(f"CFD ピン: 凍結源は Hall 初期線の場に限る (この run の initial_line.source={il.get('source')!r} — 循環定義)")
    if info.get("start_line", "throat_char") != "throat_char":
        raise ValueError(f"CFD ピン: 凍結源の start_line={info.get('start_line')!r} (throat_char でない)")
    S = float(info["scale_m"])
    ni = int(info["mesh"]["ni"])
    with h5py.File(rd / "nozzle.h5") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    X = (nc[:, 0] / S).reshape(ni, -1)
    R = (nc[:, 1] / S).reshape(ni, -1)
    if isinstance(res, (str, Path)):
        res = [res]
    files = [rd / r for r in res]
    Ms, Ts, snap = [], [], []
    for fp in files:
        snap.append(_sha16(fp.read_bytes()))
        with h5py.File(fp) as f:
            Ux, Uy, son = (f["/VALUE/" + k][:] for k in ("Ux", "Uy", "sonic"))
        Ms.append(np.hypot(Ux, Uy) / son)
        Ts.append(np.arctan2(Uy, Ux))
    thermo = {"gas": info.get("gas"), "gamma_hall": info.get("gamma_hall"),
              "physProp": {k: v for k, v in (yaml_physprop(rd)).items()}}
    src = {"run": str(rd), "res": [f.name for f in files],
           "config": cfg_checked,
           "shape": {"R": info.get("R"), "L_U": info.get("L_U"), "scale_m": S,
                     "x_in_mesh": float(X[0, 0]), "r_U_mesh": float(R[0, -1])},
           "gamma_hall": info.get("gamma_hall"), "gas": info.get("gas"),
           "sha256_16": {"snapshot": snap[0] if len(snap) == 1 else snap,
                         "thermo": _sha16(json.dumps(thermo, sort_keys=True, default=str).encode())},
           "_wall_x": X[:, -1].copy(), "_wall_r": R[:, -1].copy()}
    return (X, R, np.mean(Ms, 0).reshape(ni, -1), np.mean(Ts, 0).reshape(ni, -1), src)


def yaml_physprop(rd: Path) -> dict:
    import yaml
    return dict((yaml.safe_load((Path(rd) / "solverConfig.yaml").read_text()) or {}).get("physProp") or {})


def check_source_matches(src: dict, R: float, gamma: float, expect: dict | None = None,
                         tol_rel: float = 1e-9, tol_wall: float = 5e-6) -> dict:
    """凍結源の形・熱力学条件を使う側と照合し、食い違えば拒否する (codex result M2: R=3・γ=1.4 を受理していた)。

    - R・γ_Hall (prepare_info の値) と使う側の R・γ_Hall
    - expect (design_chain から): gas (kind・Tt・組成 Y)、r_U・L_U・L_pipe と U→T Hermite —
      凍結源のメッシュの壁節点 (入口直管と縮流部) を使う側の設計縮流部 `UpstreamThroatPoly` と比べる
    tol_wall [r_t]: メッシュ座標は float32 (case/45 の r_U 6.485 で 1 ulp ≈ 5e-7 r_t、実測の残差 5.1e-7) なので 5e-6。
    戻り: 照合結果 (照合した項目と照合できなかった項目)。"""
    errs, checked, unchecked = [], [], []
    for key, mine in (("R", R), ("gamma_hall", gamma)):
        v = src.get("shape", {}).get(key) if key == "R" else src.get(key)
        if v is None:
            unchecked.append(key)
        elif abs(float(v) - float(mine)) > tol_rel * max(1.0, abs(float(mine))):
            errs.append(f"{key}: 凍結源 {v} vs 使う側 {mine}")
        else:
            checked.append(key)
    expect = expect or {}
    g_src, g_me = src.get("gas") or {}, expect.get("gas")
    if g_me is not None:
        if not g_src:
            unchecked.append("gas")
        else:
            if str(g_src.get("kind")) != str(g_me.get("kind")):
                errs.append(f"gas.kind: {g_src.get('kind')} vs {g_me.get('kind')}")
            for k in ("Tt",):
                if g_src.get(k) is not None and g_me.get(k) is not None and abs(float(g_src[k]) - float(g_me[k])) > 1e-9:
                    errs.append(f"gas.{k}: {g_src[k]} vs {g_me[k]}")
            ys, ym = g_src.get("Y") or {}, g_me.get("Y") or {}
            if set(ys) != set(ym) or any(abs(float(ys[k]) - float(ym[k])) > 1e-12 for k in ys):
                errs.append(f"gas.Y: {ys} vs {ym}")
            checked.append("gas")
    if "L_U" in expect:
        sh = src.get("shape", {})
        L_U, r_U, L_pipe = float(expect["L_U"]), float(expect["r_U"]), float(expect.get("L_pipe", 0.5))
        if sh.get("L_U") is not None and abs(float(sh["L_U"]) - L_U) > 1e-12:
            errs.append(f"L_U: 凍結源 {sh['L_U']} vs 使う側 {L_U}")
        if abs(sh["x_in_mesh"] - (-L_U - L_pipe)) > tol_wall:
            errs.append(f"入口位置 (L_U+L_pipe): 凍結源のメッシュ x_in {sh['x_in_mesh']:.6f} vs 使う側 {-L_U - L_pipe:.6f}")
        if abs(sh["r_U_mesh"] - r_U) > tol_wall:
            errs.append(f"r_U: 凍結源のメッシュ {sh['r_U_mesh']:.6f} vs 使う側 {r_U:.6f}")
        from ..geometry.wall_walldriven import UpstreamThroatPoly
        up = UpstreamThroatPoly(r_U=r_U, R_t=float(R), L_U=L_U)
        xw, rw = src["_wall_x"], src["_wall_r"]
        m = (xw > -L_U) & (xw < 0.0)
        if m.sum() < 10:
            errs.append("凍結源のメッシュに縮流部の壁節点が少ない")
        else:
            dmax = float(np.max(np.abs(up.r(xw[m]) - rw[m])))
            if dmax > tol_wall:
                errs.append(f"縮流部の形: 凍結源のメッシュ壁と使う側の U→T Hermite の差 {dmax:.2e} > {tol_wall:g}")
        checked += ["L_U", "r_U", "L_pipe", "contraction_shape"]
    if errs:
        raise ValueError("CFD ピン: 凍結源が使う側と同じ形・ガスでない — " + "; ".join(errs))
    return {"checked": checked, "unchecked": unchecked}


class CFDPinnedThroat(HallThroat):
    """`HallThroat` 互換の CFD ピン provider (初期線・軸アンカー・場の M/θ を CFD 場から)。"""

    def __init__(self, R, gamma, X, Rr, M, TH, source=None, ds=2e-4, x_window=(-2.0, 2.5),
                 anchor_halfwidth=0.25):
        super().__init__(R=R, gamma=gamma)
        self.source = dict(source or {})
        if self.source.get("run"):
            # 凍結源の形・熱力学条件と使う側の照合 (source が run 由来のときだけ。合成場の単体試験は source 無し)
            self.source["match"] = check_source_matches(self.source, R, gamma, self.source.pop("_expect", None))
        for k in ("_wall_x", "_wall_r", "_expect"):
            self.source.pop(k, None)
        xs = X[:, 0]
        eta = Rr[0] / Rr[0, -1]
        if np.abs(X - xs[:, None]).max() > 1e-12 or np.abs(Rr / Rr[:, -1:] - eta[None, :]).max() > 1e-6:
            raise ValueError("断面が x 一定でない / η 分布が断面ごとに違う")
        m = (xs >= x_window[0]) & (xs <= x_window[1])
        self._xs, self._eta = xs[m], eta
        self._rw = CubicSpline(xs, Rr[:, -1])
        self._Ms = RectBivariateSpline(self._xs, eta, M[m], kx=3, ky=3)
        self._Ts = RectBivariateSpline(self._xs, eta, TH[m], kx=3, ky=3)
        # 軸の evenfit: 各断面で r>0 の 4 点に M = a0 + a2 r²
        rr = Rr[m, 1:5]
        mm = M[m, 1:5]
        self._ax_x = self._xs
        self._ax_M = np.array([np.polyfit(rr[i] ** 2, mm[i], 1)[1] for i in range(len(self._xs))])
        self._build_line(ds)
        self._anchor_hw = anchor_halfwidth
        self.source.setdefault("sha256_16", {})["line"] = _sha16(
            np.concatenate([self._line_x, self._line_r, self._line_M, self._line_T]).tobytes())

    # -- 場 --
    def _field(self, x, r):
        x = float(x)
        e = min(max(float(r) / float(self._rw(x)), 0.0), 1.0)
        return float(np.ravel(self._Ms.ev(x, e))[0]), float(np.ravel(self._Ts.ev(x, e))[0])

    def _build_line(self, ds):
        def slope(p):
            if not (self._xs[0] <= p[0] <= self._xs[-1]):
                raise ValueError(f"CFD ピン: 追跡点 x={p[0]:.4f} が場の窓 "
                                 f"[{self._xs[0]:.3f}, {self._xs[-1]:.3f}] の外")
            Mv, tv = self._field(p[0], p[1])
            if not (np.isfinite(Mv) and np.isfinite(tv)):
                raise ValueError(f"CFD ピン: 非有限の場 (x={p[0]:.4f}, r={p[1]:.4f})")
            if Mv <= 1.0:
                raise ValueError(f"CFD ピン: 亜音速 M={Mv:.5f} (x={p[0]:.4f}, r={p[1]:.4f}) — "
                                 "スロート特性線が引けない")
            a = tv - np.arcsin(1.0 / Mv)
            return np.r_[np.cos(a), np.sin(a)]
        p = np.r_[0.0, 1.0 - 1e-6]
        P = [p.copy()]
        while p[1] > 0 and len(P) < 400000:
            k1 = slope(p)
            k2 = slope(p + 0.5 * ds * k1)
            k3 = slope(p + 0.5 * ds * k2)
            k4 = slope(p + ds * k3)
            p = p + ds * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            P.append(p.copy())
        if p[1] > 0:
            raise ValueError("CFD ピン: 軸に到達しない")
        P = np.array(P[:-1])
        P = P[P[:, 1] > 0]                                                 # 軸を跨いだ最後の点は捨てる
        F = np.array([self._field(a, b) for a, b in P])
        nr = P[:, 1] <= 0.05
        ax = [np.polyval(np.polyfit(P[nr, 1], v, 2), 0.0) for v in (P[nr, 0], F[nr, 0], F[nr, 1])]
        top = np.argsort(-P[:, 1])[:10]                                    # 壁足の M: 内側 10 点の線形外挿
        M_foot = float(np.polyval(np.polyfit(P[top, 1], F[top, 0], 1), 1.0))
        o = np.argsort(P[:, 1])
        self._line_r = np.r_[0.0, P[o, 1], 1.0]
        self._line_x = np.r_[ax[0], P[o, 0], 0.0]
        self._line_M = np.r_[ax[1], F[o, 0], M_foot]
        self._line_T = np.r_[0.0, F[o, 1], 0.0]
        self.x0_cfd, self.M_axis_line = float(ax[0]), float(ax[1])

    # -- HallThroat 互換 API の上書き --
    def throat_characteristic(self, n: int = 61, **_):
        r = np.linspace(0.0, 1.0, n)
        return (np.interp(r, self._line_r, self._line_x), r,
                np.interp(r, self._line_r, self._line_M), np.interp(r, self._line_r, self._line_T))

    def axis_anchor(self, x: float) -> tuple:
        x = float(x)
        if abs(x - self.x0_cfd) > 1e-9:
            raise ValueError(f"CFD ピンの軸アンカーは線の軸着地 x0={self.x0_cfd:.6f} でだけ定義 (要求 x={x:.6f})")
        w = np.abs(self._ax_x - x) <= self._anchor_hw
        c = np.polyfit(self._ax_x[w] - x, self._ax_M[w], 4)
        Mpp_hall = super().axis_anchor(x)[2]
        return self.M_axis_line, float(c[-2]), float(Mpp_hall)

    def mach(self, x, r):
        x = np.atleast_1d(np.asarray(x, float))
        r = np.broadcast_to(np.asarray(r, float), x.shape)
        out = np.array([self._field(a, b)[0] for a, b in zip(x, r)])
        return out if out.size > 1 else float(out[0])

    def theta(self, x, r):
        x = np.atleast_1d(np.asarray(x, float))
        r = np.broadcast_to(np.asarray(r, float), x.shape)
        out = np.array([self._field(a, b)[1] for a, b in zip(x, r)])
        return out if out.size > 1 else float(out[0])


def pinned_factory(run_dir, res=None, expect: dict | None = None):
    """`HallThroat(R=, gamma=)` と同じ呼び方で CFD ピンの throat を返す factory (場は 1 回だけ読む)。
    res (snapshot) は必須。expect (gas・r_U・L_U・L_pipe) を渡すと凍結源の形・ガスを照合する (`check_source_matches`)。"""
    X, Rr, M, TH, src = load_run_field(run_dir, res)

    def make(R, gamma):
        s = copy.deepcopy(src)
        s["_expect"] = expect
        return CFDPinnedThroat(R, gamma, X, Rr, M, TH, source=s)
    return make
