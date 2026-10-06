"""逆 MOC マーチ (モード F 逆設計): 軸目標 M(x) → 場充填 → 壁流線抽出。

親計画 §4.7 v1 の実体。数理 (2026-08-13 の設計議論で確定):

- 軸は条件 2 個 (M 指定 + θ=0) を持つ Cauchy データ線。初期データ =
  starting line (遷音速パッチ出力, 軸→壁) + 軸目標点列 (下流へ延長) の L 字。
- **三角充填**: データ列の隣接ペア (右=C⁻ 担体, 左=C⁺ 担体) から interior
  単位プロセスで 1 段上の点を作る。レベルごとに 1 点ずつ減り、決定領域 =
  「左端 C⁺ (starting line 壁点から) と右端 C⁻ (最下流軸点へ降りる) に挟まれた
  三角形」を過不足なく埋める。壁点の依存域は軸の両側に足を持つ (C⁺ 足=上流,
  C⁻ 足=下流) ため、**軸目標は壁端の C⁻ 足まで下流延長が必要** (一様出口では
  M_d 一定でタダ)。
- **壁抽出**: 場を線形補間し、starting line 壁端から流線 dr/dx = tanθ を積分。
  質量流量の整合 (∫ρu dA vs スロート) は診断として併記。

単位プロセス (interior/軸対称源項/予測子修正子) は moc_kernel.KernelMOC を継承。
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import LinearNDInterpolator

from . import moc_kernel as _mk
from .moc_kernel import (AXIS_LIMIT_MODES, AXIS_R_EPS, CORR_MAX, CORR_TOL,
                         CORRECTOR_MODES, PAIR_BELOW_AXIS, PAIR_CONVERGED, PAIR_DONE,
                         PAIR_GEOM, PAIR_MAXITER, PAIR_MISSING, PAIR_NONFINITE,
                         RESID_TOL, KernelMOC, _Pt, axis_theta_r, interior_vec,
                         pm_mach, pm_mach_vec, pm_nu, source_branches)


def field_interpolator(pts):
    r"""全計算点の Delaunay を **1 回だけ**作り、$(\theta, \nu)$ をまとめて返す補間器。

    充填をベクトル化した後 (2026-08-15) は**三角形分割が支配コスト**になった
    (実測 n_axis=1000 で 1 個 8.4 s、n_axis=2000 で ~35 s)。従来は θ 用と ν 用に
    別々の `LinearNDInterpolator` を作り、さらに壁流線・質量流束診断・終端特性線が
    それぞれ作り直していたため同じ分割を 6 回構築していた。値を (n,2) にまとめ、
    生成した補間器を呼び出し側で使い回すことで 1 回に減らす (数値は完全に同一 —
    同じ点集合の同じ線形補間)。"""
    xy = np.array([[p.x, p.r] for p in pts])
    val = np.array([[p.th, p.nu] for p in pts])
    return LinearNDInterpolator(xy, val)


def _mass_flux_density(M, g):
    """等エントロピー流の ρ|V|/(ρ0 a0) (よどみ量規格化)。ガスモデルなら ρV/(ρ*a*)
    (規格化定数が違うが、本モジュールでは**比**しか使わないので等価)。"""
    if hasattr(g, "mass_flux_density"):
        return np.asarray(g.mass_flux_density(M), dtype=float)
    M = np.asarray(M, dtype=float)
    t = 1.0 + 0.5 * (g - 1.0) * M * M
    return M * t ** (-0.5 * (g + 1.0) / (g - 1.0))


class MocFillDiag:
    r"""充填 1 回分の単位過程の集計 (plans/active/discretization-moc-axis-limit-and-corrector.md §4.2)。

    対の 5 分類: 入力時点の対象 (`input`) = もともとの欠損 (`missing`) + 幾何的棄却 (`geom_parallel`・
    `geom_below_axis`) + 反復の失敗 (`iter_nonfinite`・`iter_maxiter`) + 収束 (`converged`、fixed2 は
    固定回数を回し終えた `done_fixed`)。反復回数の分布、最終残差の最大、源項の分岐 (`source_branches`) の数と
    `AXIS_LIMIT_FRAC` が発火した位置、幾何的棄却・反復の失敗の位置も持つ。値は変えない (記録専用)。"""

    CAP = 200          # 位置の記録の上限 (数は全数)

    def __init__(self, axis_limit: str, corrector: str, n_corr: int, tol: float, max_corr: int):
        self.cfg = {"axis_limit": axis_limit, "corrector": corrector, "n_corr": int(n_corr),
                    "tol": float(tol), "max_corr": int(max_corr), "resid_tol": RESID_TOL}
        self.counts = dict.fromkeys(("input", "missing", "geom_parallel", "geom_below_axis",
                                     "iter_nonfinite", "iter_maxiter", "converged", "done_fixed"), 0)
        self.iter_hist: dict = {}
        self.resid = {"geom_max": 0.0, "comp_max": 0.0, "at_geom": None, "at_comp": None}
        self.branch = dict.fromkeys(("axis_analytic", "axis_borrow", "axis_zero", "frac_A", "frac_B",
                                     "frac_PA", "frac_PB", "P_on_axis"), 0)
        self.frac_pos: list = []           # [level, i, 評価, x_p, r_p, r_相手]
        self.frac_rng = [np.inf, -np.inf, np.inf, -np.inf]     # x_min, x_max, r_min, r_max
        self.geom = []                     # 幾何的棄却 (全数): (level, i, kind, xA, rA, xB, rB)
        self.fail: list = []               # 反復の失敗: [level, i, kind, n_iter, xA, rA, xB, rB]
        self.n_levels = 0

    def add(self, level: int, Ax, Ar, Bx, Br, xP, rP, st: dict, below):
        status = np.array(st["status"], dtype=np.int8)
        status[(status == PAIR_CONVERGED) & below] = PAIR_BELOW_AXIS
        status[(status == PAIR_DONE) & below] = PAIR_BELOW_AXIS
        n_it = st["n_iter"]
        c = self.counts
        c["input"] += int(status.size)
        for key, code in (("missing", PAIR_MISSING), ("geom_parallel", PAIR_GEOM),
                          ("geom_below_axis", PAIR_BELOW_AXIS), ("iter_nonfinite", PAIR_NONFINITE),
                          ("iter_maxiter", PAIR_MAXITER), ("converged", PAIR_CONVERGED),
                          ("done_fixed", PAIR_DONE)):
            c[key] += int(np.sum(status == code))
        self.n_levels = max(self.n_levels, int(level))
        done = (status == PAIR_CONVERGED) | (status == PAIR_DONE)
        if np.any(done):
            v, cnt = np.unique(n_it[done], return_counts=True)
            for a, b in zip(v.tolist(), cnt.tolist()):
                self.iter_hist[int(a)] = self.iter_hist.get(int(a), 0) + int(b)
        for key, arr in (("geom_max", st["resid_geom"]), ("comp_max", st["resid_comp"])):
            fa = np.where(np.isfinite(arr), arr, -1.0)
            if fa.size and fa.max() > self.resid[key]:
                j = int(np.argmax(fa))
                self.resid[key] = float(fa[j])
                self.resid["at_" + key.split("_")[0]] = [int(level), j, float(xP[j]), float(rP[j])]
        b = st["branch"]
        valid = status != PAIR_MISSING
        for key in ("axis_analytic", "axis_borrow", "axis_zero"):
            self.branch[key] += int(np.sum(np.where(valid, b[key], 0)))
        for key, (xp, rp, ro) in (("frac_A", (Ax, Ar, Br)), ("frac_B", (Bx, Br, Ar)),
                                  ("frac_PA", (xP, np.maximum(rP, 0.0), Ar)),
                                  ("frac_PB", (xP, np.maximum(rP, 0.0), Br))):
            m = b[key] & valid
            nm = int(np.sum(m))
            self.branch[key] += nm
            if nm:
                idx = np.flatnonzero(m)
                self.frac_rng = [min(self.frac_rng[0], float(xp[idx].min())), max(self.frac_rng[1], float(xp[idx].max())),
                                 min(self.frac_rng[2], float(rp[idx].min())), max(self.frac_rng[3], float(rp[idx].max()))]
                for j in idx[: max(self.CAP - len(self.frac_pos), 0)]:
                    self.frac_pos.append([int(level), int(j), key, float(xp[j]), float(rp[j]), float(ro[j])])
        self.branch["P_on_axis"] += int(np.sum(b["P_on_axis"] & ((status == PAIR_CONVERGED) | (status == PAIR_DONE))))
        for code, kind in ((PAIR_GEOM, "parallel"), (PAIR_BELOW_AXIS, "below_axis")):
            for j in np.flatnonzero(status == code):
                self.geom.append((int(level), int(j), kind, float(Ax[j]), float(Ar[j]), float(Bx[j]), float(Br[j])))
        for code, kind in ((PAIR_NONFINITE, "nonfinite"), (PAIR_MAXITER, "maxiter")):
            for j in np.flatnonzero(status == code)[: max(self.CAP - len(self.fail), 0)]:
                self.fail.append([int(level), int(j), kind, int(n_it[j]), float(Ax[j]), float(Ar[j]),
                                  float(Bx[j]), float(Br[j])])

    def summary(self) -> dict:
        h = self.iter_hist
        n = sum(h.values())
        it = ({"max": int(max(h)), "mean": float(sum(k * v for k, v in h.items()) / n), "n_pairs": int(n),
               "hist": {str(k): int(h[k]) for k in sorted(h)}} if n else
              {"max": None, "mean": None, "n_pairs": 0, "hist": {}})
        fr = None if not np.isfinite(self.frac_rng[0]) else {
            "x_min": self.frac_rng[0], "x_max": self.frac_rng[1], "r_min": self.frac_rng[2], "r_max": self.frac_rng[3]}
        return {**self.cfg, "pairs": dict(self.counts), "iters": it,
                "resid": (dict(self.resid) if self.cfg["corrector"] == "converge" else None),
                "branch": dict(self.branch),
                "axis_limit_frac": {"frac": float(_mk.AXIS_LIMIT_FRAC), "n_evals": int(sum(self.branch[k] for k in
                                                                     ("frac_A", "frac_B", "frac_PA", "frac_PB"))),
                                    "range": fr, "positions": self.frac_pos,
                                    "positions_note": f"[level, i, 評価, x_p, r_p, r_相手] (先頭 {self.CAP} 件)"},
                "geom_rejects": [list(g) for g in self.geom[: self.CAP]],
                "_geom_all": list(self.geom),         # moc_gate の判定用 (記録からは外す)
                "iter_failures": self.fail, "n_levels": int(self.n_levels)}


def _posthoc_stats(A, B, out, n_corr: int) -> dict:
    """既定 (legacy + fixed2) の対の分類を、従来の呼び出しの入出力から作る (値は変えない)。
    `interior_vec` を位置引数だけで呼ぶ (試験スクリプトが単位過程を差し替える経路との互換) ため、
    幾何的棄却と反復中の非有限の区別は出力の有限性で付ける (棄却対は den=1 で続行されるので通常は有限)。"""
    xP, rP, thP, nuP, ok = out
    fin_in = np.ones(np.shape(xP), dtype=bool)
    for a in (*A, *B):
        fin_in &= np.isfinite(a)
    fin_out = np.isfinite(xP) & np.isfinite(rP) & np.isfinite(thP) & np.isfinite(nuP)
    st = np.full(np.shape(xP), PAIR_DONE, dtype=np.int8)
    st[fin_in & ~ok & fin_out] = PAIR_GEOM
    st[fin_in & ~ok & ~fin_out] = PAIR_NONFINITE
    st[~fin_in] = PAIR_MISSING
    return {"status": st, "n_iter": np.full(np.shape(xP), int(n_corr), dtype=np.int16),
            "resid_geom": np.full(np.shape(xP), np.nan), "resid_comp": np.full(np.shape(xP), np.nan),
            "branch": source_branches(A[1], B[1], rP)}


class InverseMOC(KernelMOC):
    """軸 Cauchy データからの三角充填と壁流線抽出。

    `axis_limit` / `corrector` は `moc_kernel.interior_vec` の選択肢 (既定 legacy + fixed2 = 従来とビット同一。
    plans/active/discretization-moc-axis-limit-and-corrector.md)。充填のたびに対の集計を `last_diag` に置く。"""

    def __init__(self, gamma=1.4, delta=1.0, n_corr=2, axis_limit: str = "legacy",
                 corrector: str = "fixed2", tol: float = CORR_TOL, max_corr: int = CORR_MAX):
        super().__init__(gamma=gamma, delta=delta, n_corr=n_corr)
        if axis_limit not in AXIS_LIMIT_MODES:
            raise ValueError(f"axis_limit は {AXIS_LIMIT_MODES} のどれか (受け取った値: {axis_limit!r})")
        if corrector not in CORRECTOR_MODES:
            raise ValueError(f"corrector は {CORRECTOR_MODES} のどれか (受け取った値: {corrector!r})")
        _mk._check_corr_params(tol, max_corr)
        self.axis_limit, self.corrector = axis_limit, corrector
        self.tol, self.max_corr = float(tol), int(max_corr)
        self.last_diag: dict | None = None

    def _new_diag(self) -> MocFillDiag:
        return MocFillDiag(self.axis_limit, self.corrector, self.n_corr, self.tol, self.max_corr)

    def _front_thr(self, n: int, axis_thr):
        """初期前線の θ_r (analytic のときだけ。軸端点以外は NaN)。legacy では使わない (None)。"""
        if self.axis_limit != "analytic":
            return None
        if axis_thr is None:
            raise ValueError("axis_limit='analytic' には初期前線の軸端点の θ_r (axis_thr) が必要")
        thr = np.asarray(axis_thr, dtype=float)
        if thr.shape != (n,):
            raise ValueError(f"axis_thr の長さ {thr.shape} が初期前線 ({n}) と違う")
        return thr

    def _pairs(self, x, r, th, nu, M, thr):
        """隣接対 (B=front[:-1] が C⁻ 担体、A=front[1:] が C⁺ 担体) を単位過程に通す。戻り: 5 値 + 集計用 stats。"""
        A = (x[1:], r[1:], th[1:], nu[1:], M[1:])
        B = (x[:-1], r[:-1], th[:-1], nu[:-1], M[:-1])
        if self.axis_limit == "legacy" and self.corrector == "fixed2":
            # 既定: 従来と同じ呼び出し (位置引数だけ — 単位過程を差し替える試験スクリプトとの互換)
            out = interior_vec(*A, *B, self.g, self.delta, self.n_corr)
            return out + (_posthoc_stats(A, B, out, self.n_corr),)
        st: dict = {}
        out = interior_vec(*A, *B, self.g, self.delta, self.n_corr, axis_limit=self.axis_limit,
                           thrA=None if thr is None else thr[1:], thrB=None if thr is None else thr[:-1],
                           corrector=self.corrector, tol=self.tol, max_corr=self.max_corr, stats=st)
        return out + (st,)

    # -- 場充填 ---------------------------------------------------------------
    def fill(self, init_front, axis_thr=None) -> list:
        """init_front: _Pt 列。**最下流軸点 → 上流軸点 → starting line を壁へ**
        の順 (L 字に沿って単調)。戻り値: 全計算点 (init 含む)。

        **既知の限界と、棄却された「修正」の記録 (2026-08-14)**: 単位プロセスの
        役割 (A=C⁺ 担体 / B=C⁻ 担体) は本来前線の向きに依存し、縦の starting
        line 区間では下の点が C⁺ 担体であるべき。この固定割当てでは縦線区間で
        逆転するため交点が棄却され、**スロート直後の楔領域が無計算のまま残る**
        (壁流線はそこを凸包補間で跨ぐ)。そこで「両割当てを試し前方交点を採る」
        向き非依存の充填を試したが、**実設計は悪化した** — 楔が埋まる代わりに
        誤ったデータで埋まり、壁流線がそれを拾うため:

        | 設計 | 出口 ε_M | コア M (目標 4.0) |
        | --- | --- | --- |
        | 本実装 (楔は空・凸包補間) | 0.173% | 3.9995 |
        | 向き非依存 fill | 0.932% | 4.0353 |

        よって**本実装 (楔を空のまま残す) を維持**する。なお初期値線をスロート
        特性線にした (A8) 時点で縦線区間が無くなり、この楔自体が構造的に消えた。
        壁を古典的に閉じる経路は `cplus_flux_wall` (A10) を使う。

        `axis_thr`: `axis_limit='analytic'` のときの初期前線の θ_r (`fill_arrays` と同じ)。"""
        arr = self.fill_arrays(init_front, axis_thr=axis_thr)
        g = self.g
        return [_Pt(float(a[0]), float(a[1]), float(a[2]), float(a[3]), g, float(a[4]))
                for a in arr]

    def fill_arrays(self, init_front, axis_thr=None) -> np.ndarray:
        r"""`fill` の**フロント一括ベクトル版**。戻り: (n,5) [x, r, th, nu, M]。

        レベルごとに隣接ペアを**まとめて**単位過程に通す (計算内容はスカラー版と
        同一で、ループの順序だけが変わる)。スカラー版は 1 点ずつ Python で回すため
        点数 $\sim n^2/2$ に比例した Python オーバヘッドが支配していた
        (実測: $n_{axis}$=2000 で 369 秒)。ベクトル版は**レベル数 $n$ 回**の
        numpy 呼び出しで済む (同 2.4 秒, 150 倍)。numba/C++ を持ち込まずに済むのは、
        1 レベル内のペアが互いに独立だから (レベル間の依存だけが逐次)。

        `axis_thr`: `axis_limit='analytic'` のときの初期前線の各点の θ_r (長さ = 初期前線、
        軸端点以外は NaN)。生成した点 (軸外) は θ_r を持たない。対の集計は `self.last_diag`。
        """
        x = np.array([p.x for p in init_front], dtype=float)
        r = np.array([p.r for p in init_front], dtype=float)
        th = np.array([p.th for p in init_front], dtype=float)
        nu = np.array([p.nu for p in init_front], dtype=float)
        M = np.array([p.M for p in init_front], dtype=float)
        thr = self._front_thr(len(x), axis_thr)
        diag = self._new_diag()
        out = [np.column_stack([x, r, th, nu, M])]
        k = 0
        while len(x) >= 2:
            k += 1
            # B=右(下流)側 C⁻ 担体 = front[:-1], A=左側 C⁺ 担体 = front[1:]
            xP, rP, thP, nuP, ok, st = self._pairs(x, r, th, nu, M, thr)
            below = ~(rP >= -1e-12)
            ok &= rP >= -1e-12
            diag.add(k, x[1:], r[1:], x[:-1], r[:-1], xP, rP, st, below & np.isfinite(rP))
            if not np.any(ok):
                break
            x, r, th, nu = xP[ok], rP[ok], thP[ok], nuP[ok]
            M = pm_mach_vec(nu, self.g)
            thr = None if thr is None else np.full(len(x), np.nan)
            out.append(np.column_stack([x, r, th, nu, M]))
        self.last_diag = diag.summary()
        return np.vstack(out)

    def fill_levels(self, init_front, axis_thr=None) -> np.ndarray:
        r"""`fill_arrays` の**レベル構造を保ったまま**返す版。戻り: (n_lev, n_pt, 5)。

        `[k, i]` = レベル $k$ の位置 $i$ の点 $[x, r, \theta, \nu, M]$、
        欠損 (棄却された対・前線の縮み) は NaN。**添字を詰めない**のが要点で、
        これにより特性線が添字だけで読み出せる (`cplus_lines`)。

        $L_k[i]$ は $L_{k-1}[i]$ (C⁻ 担体) と $L_{k-1}[i+1]$ (C⁺ 担体) から作るので:

        - **C⁻ 線** = 添字固定の列 $L_k[i],\ k=0,1,\dots$
        - **C⁺ 線** = 反対角線 $L_k[m-k],\ k=0,1,\dots,m$ (起点 $L_0[m]$)

        つまり三角充填の網と特性線網は同じもので、走査順が違うだけ。

        `axis_thr`: `fill_arrays` と同じ。対の集計 (欠損 = 死んだ対を含む) は `self.last_diag`。
        """
        x = np.array([p.x for p in init_front], dtype=float)
        r = np.array([p.r for p in init_front], dtype=float)
        th = np.array([p.th for p in init_front], dtype=float)
        nu = np.array([p.nu for p in init_front], dtype=float)
        M = np.array([p.M for p in init_front], dtype=float)
        n = len(x)
        thr = self._front_thr(n, axis_thr)
        diag = self._new_diag()
        out = np.full((n, n, 5), np.nan)
        out[0] = np.column_stack([x, r, th, nu, M])
        live = np.ones(n, dtype=bool)          # レベル k で有効な添字
        for k in range(1, n):
            # 対 (i, i+1) がともに有効なときだけ新点を作る (添字は i を継承)
            pair = live[:-1] & live[1:]
            if not np.any(pair):
                break
            xP, rP, thP, nuP, ok, st = self._pairs(x, r, th, nu, M, thr)
            diag.add(k, x[1:], r[1:], x[:-1], r[:-1], xP, rP, st, pair & np.isfinite(rP) & ~(rP >= -1e-12))
            ok &= pair & (rP >= -1e-12)
            thr = None if thr is None else np.full(len(x) - 1, np.nan)    # 生成した点は軸外 (θ_r なし)
            if not np.any(ok):
                break
            MP = pm_mach_vec(np.where(ok, nuP, 0.0), self.g)
            x = np.where(ok, xP, np.nan)
            r = np.where(ok, rP, np.nan)
            th = np.where(ok, thP, np.nan)
            nu = np.where(ok, nuP, np.nan)
            M = np.where(ok, MP, np.nan)
            out[k, :len(x)] = np.column_stack([x, r, th, nu, M])
            live = ok
        self.last_diag = diag.summary()
        return out

    # -- 壁流線 ---------------------------------------------------------------
    def wall_streamline(self, pts, x_start: float, r_start: float,
                        dx: float = 0.02, itp=None) -> np.ndarray:
        """(x_start, r_start) から dr/dx = tanθ を RK2 積分。場外に出たら終了。

        `itp`: `field_interpolator(pts)` を使い回す場合に渡す (省略時は自前で構築)。
        戻り値: (n,4) [x, r, theta, M]。"""
        itp = field_interpolator(pts) if itp is None else itp
        out = []
        x, r = float(x_start), float(r_start)
        t0, n0 = (float(v) for v in itp(x, r))
        while np.isfinite(t0):
            out.append((x, r, t0, pm_mach(max(n0, 1e-9), self.g)))
            xm, rm = x + 0.5 * dx, r + 0.5 * dx * np.tan(t0)
            tm = float(itp(xm, rm)[0])
            if not np.isfinite(tm):
                break
            x, r = x + dx, r + dx * np.tan(tm)
            t0, n0 = (float(v) for v in itp(x, r))
        return np.array(out)

    # -- 質量流量診断 -----------------------------------------------------------
    def massflux_across(self, pts, x_plane: float, r_top: float, n: int = 120,
                        itp=None) -> float:
        """x=x_plane の縦断面 (0..r_top) の ∫ρu·2πr dr (よどみ量規格化)。"""
        itp = field_interpolator(pts) if itp is None else itp
        rr = np.linspace(1e-6, r_top, n)
        val = itp(np.full(n, x_plane), rr)
        tt, vv = val[:, 0], val[:, 1]
        ok = np.isfinite(tt) & np.isfinite(vv)
        M = np.array([pm_mach(max(v, 1e-9), self.g) for v in vv[ok]])
        f = _mass_flux_density(M, self.g) * np.cos(tt[ok])
        return float(np.trapezoid(f * 2.0 * np.pi * rr[ok], rr[ok]))


def flux_closure_wall(itp, mdot_star: float, xs, r_top, gamma: float = 1.4,
                      n_r: int = 400) -> np.ndarray:
    r"""**壁 = 累積質量流束が $\dot m^*$ に達する半径** (各 $x$ 断面, ベクトル化)。

    計画: plans/archived/tooling-nozzle-moc-flux-closure-wall.md (A9)。

    $$\int_0^{r_w}\rho(\mathbf u\cdot\mathbf n)\,dA
      = \int_0^{r_w} \rho V\cos\theta\,2\pi r\,dr = \dot m^*$$

    **⚠ 既定では使わない — CFD で棄却済み (2026-08-16)**。定常流では流線面と質量流束
    一定面は同一の対象なので「保存則を構成的に満たす閉包の方が良い壁が出る」と期待して
    実装したが、**実測は逆だった**。同一場 ($n_{\rm axis}$=2000) で `wall_mode` だけを
    振った CFD A/B (case/41 run_0050 [流線] vs run_0052 [本関数]):

    | 指標 | 流線 | 流束閉包 |
    | --- | --- | --- |
    | ‖ΔM‖∞ [% M_d] | 0.353 | 0.507 |
    | overshoot | +0.068% | +0.217% |
    | 出口 ε_M rms | 0.020% | 0.057% |

    厳密解を持つ放射源流でも流線 6.1e-4 に対し本閉包は +1.43e-3 (2.3 倍悪い、しかも
    $x$ によらずほぼ一定の系統バイアス)。実装は正しい (一様流で $r_w\equiv1$ を 1.7e-7)。

    **理由**: MOC 場は特性線適合関係 (運動量+等エントロピー) から作られており、
    離散レベルでは質量保存を満たさない。断面積分で壁を決めてもその不整合は消えず、
    **壁半径の誤差から壁形状の誤差へ移るだけ**。しかも断面積分は各断面の半径方向
    プロファイル全体 (軸対称ソース項 $1/r$ が効く近軸域を含む) の精度に依存するのに対し、
    流線は $\theta$ を 1 本の曲線上で積分するだけなので誤差が小さい。

    **紛らわしい指標に注意**: 本閉包を使うと `mdot_ratio_moc` は構成的に 1 になる
    (循環指標で品質の証拠にならない) し、$r_F$ も $C_D$ 整合値に近づく。どちらも
    改善に見えるが CFD は悪化する。**判定は必ず CFD の軸 M / 出口一様性で行うこと**。

    **適用範囲** (使う場合): 断面がまるごと計算領域に入る $x$ (初期値線の軸着地より
    下流) のみ。上流側は呼び出し側で流線壁と合成する。

    `itp`: `field_interpolator` の補間器。`xs`: 断面位置。`r_top`: 各断面の探索上限
    (スカラーまたは `xs` と同形)。戻り: (n,) の $r_w$。
    """
    xs = np.asarray(xs, dtype=float)
    r_top = np.broadcast_to(np.asarray(r_top, dtype=float), xs.shape)
    s = np.linspace(0.0, 1.0, n_r)                       # 断面ごとの正規化半径
    rr = r_top[:, None] * s[None, :]                     # (n_x, n_r)
    val = itp(np.repeat(xs, n_r), rr.ravel())            # 補間器呼び出しは 1 回だけ
    th = val[:, 0].reshape(rr.shape)
    nu = val[:, 1].reshape(rr.shape)
    M = pm_mach_vec(np.where(np.isfinite(nu), nu, 0.0), gamma)
    f = _mass_flux_density(M, gamma) * np.cos(th) * 2.0 * np.pi * rr
    good = np.isfinite(th) & np.isfinite(nu)
    # 場の外 (凸包外) は NaN。軸から連続する有限区間だけを積分に使う
    good &= np.cumprod(good, axis=1).astype(bool)
    f = np.where(good, f, 0.0)
    dr = np.diff(rr, axis=1)
    cum = np.concatenate([np.zeros((len(xs), 1)),
                          np.cumsum(0.5 * (f[:, 1:] + f[:, :-1]) * dr, axis=1)], axis=1)
    cum = np.where(good, cum, np.nan)
    out = np.empty(len(xs))
    for i in range(len(xs)):
        c, r_i, ok = cum[i], rr[i], good[i]
        n_ok = int(ok.sum())
        if n_ok < 3:
            out[i] = np.nan
            continue
        c, r_i = c[:n_ok], r_i[:n_ok]
        if c[-1] >= mdot_star:
            out[i] = float(np.interp(mdot_star, c, r_i))
        else:
            # 場の縁まで使っても不足 — 縁の状態を凍結して外挿 (薄い帯の古典的処理)。
            # dṁ/dr = k·r (k = ρV cosθ·2π 一定) から r_w を陽に解く。
            k = f[i, n_ok - 1] / max(r_i[-1], 1e-30)
            out[i] = float(np.sqrt(max(r_i[-1] ** 2
                                       + 2.0 * (mdot_star - c[-1]) / max(k, 1e-30), 0.0)))
    return out


def cplus_lines(levels: np.ndarray, m: int) -> np.ndarray:
    r"""レベル配列から**起点 $L_0[m]$ の C⁺ 線**を取り出す (`fill_levels` の反対角線)。

    $L_k[i]$ は $A=L_{k-1}[i+1]$ を C⁺ 担体として作られるので、$L_0[m]$ から出る
    C⁺ 線は $L_k[m-k]$ ($k=0,\dots,m$)。**計算は一切せず添字を読むだけ**。
    戻り: (L,5)、NaN の手前まで。"""
    n = levels.shape[1]
    d = np.diagonal(levels[:, ::-1], offset=n - 1 - m, axis1=0, axis2=1).T
    ok = np.isfinite(d[:, 0]) & np.isfinite(d[:, 1])
    cut = int(np.argmin(ok)) if not ok.all() else len(ok)
    return d[:cut]


def cplus_flux_wall(levels: np.ndarray, init_cum: np.ndarray, mdot_star: float,
                    gamma: float = 1.4) -> np.ndarray:
    r"""**壁点 = C⁺ 線上で累積質量流束が $\dot m^*$ に達する点** (古典的な逆設計閉包)。

    計画: plans/active/tooling-nozzle-moc-wall-unit-process.md。

    定常流では壁は流線であり、流線は質量流束一定面なので、「軸から数えて
    $\dot m^*$ 分の流れが通った高さ」がその位置の壁になる。各 C⁺ 線に沿って

    $$\Delta\dot m = \rho V(\cos\theta\,\Delta r - \sin\theta\,\Delta x)\cdot 2\pi\bar r$$

    を足し上げ、$\dot m^*$ を跨いだ線分を線形補間して切る。

    **`flux_closure_wall` (A9, 棄却済) との決定的な違いは積分する線と使う値**:
    あちらは鉛直断面 × Delaunay 補間場だったので半径方向プロファイル全体の
    補間誤差を拾った。本関数は **C⁺ 線 × 網の点そのもの**で、補間を一切通さない。
    壁が網を作る過程から出てくるので、Delaunay も流線 ODE 積分も不要になる
    (`wall_streamline` は補間場の中で $dr/dx=\tan\theta$ を積分するため、
    誤差が下流へ累積し、それを抑えるために網を極端に細かくする必要があった)。

    `init_cum[m]`: 起点 $L_0[m]$ までに**すでに軸から流れている**流束。初期値線上の
    点なら軸からその点までの累積、軸上の点なら 0。
    戻り: (n,4) [x, r, theta, M] を x 昇順で。"""
    n = levels.shape[1]
    out = []
    for m in range(n - 1, -1, -1):
        line = cplus_lines(levels, m)
        if len(line) < 2:
            continue
        x, r, th, nu, M = (line[:, i] for i in range(5))
        Mm = 0.5 * (M[1:] + M[:-1])
        thm = 0.5 * (th[1:] + th[:-1])
        rm = 0.5 * (r[1:] + r[:-1])
        d = (_mass_flux_density(Mm, gamma)
             * (np.cos(thm) * np.diff(r) - np.sin(thm) * np.diff(x))
             * 2.0 * np.pi * rm)
        cum = init_cum[m] + np.concatenate([[0.0], np.cumsum(d)])
        if cum[0] >= mdot_star:                 # 起点が既に壁 (初期値線の壁足)
            out.append((x[0], r[0], th[0], M[0]))
            continue
        j = int(np.argmax(cum >= mdot_star))
        if cum[j] < mdot_star:                  # 線が尽きた = 壁に届かない
            continue
        s = (mdot_star - cum[j - 1]) / max(cum[j] - cum[j - 1], 1e-30)
        out.append((x[j - 1] + s * (x[j] - x[j - 1]),
                    r[j - 1] + s * (r[j] - r[j - 1]),
                    th[j - 1] + s * (th[j] - th[j - 1]),
                    M[j - 1] + s * (M[j] - M[j - 1])))
    w = np.asarray(out, dtype=float)
    return w[np.argsort(w[:, 0])] if len(w) else w


def _smoothstep(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def _flux_along(line, g: float) -> np.ndarray:
    """折れ線 (軸→上) を横切る累積質量流束 (よどみ量規格化, 2πr 込み)。

    セグメント法線 n dℓ = (dr, −dx) (下流向き) で
    Δflux = ρV (cosθ·dr − sinθ·dx) · 2π r_mid。戻り値は各点までの累積 (len(line))。"""
    acc = [0.0]
    for P0, P1 in zip(line[:-1], line[1:]):
        Mm = 0.5 * (P0.M + P1.M)
        thm = 0.5 * (P0.th + P1.th)
        rm = 0.5 * (P0.r + P1.r)
        f = float(_mass_flux_density(Mm, g))
        acc.append(acc[-1] + f * (np.cos(thm) * (P1.r - P0.r)
                                  - np.sin(thm) * (P1.x - P0.x)) * 2.0 * np.pi * rm)
    return np.asarray(acc)



def terminal_exit(pts, wall, x_E: float, M_d: float, gamma: float = 1.4,
                  ds: float = 2e-3, n_max: int = 400000, itp=None) -> dict:
    r"""**終端特性線による物理出口 F の決定** (2026-08-15, ユーザ指摘)。

    軸上で $M=M_d$ に達する点 $E=(x_E,0)$ から出る C⁺ 特性線
    ($dr/dx=\tan(\theta+\mu)$) を MOC 場の中で追跡し、**壁流線との交点**を物理出口
    $F$ とする。この特性線は一様域の上流境界であり、$F$ より下流の断面は
    $(M_d,\theta=0)$ で埋まる — これが「一様出口」の定義そのもの。

    **なぜ壁角しきい値ではだめか**: 従来は壁流線の $\theta$ が 0.05° を切った点で
    切っていた。壁角は $F$ に向けて**漸近的に**ゼロへ近づくので、しきい値は $F$ の
    数 $r_t$ 上流で発火する (実測 M_d=4, R=2 で 2.5$r_t$ 手前)。半径差は
    $10^{-3}r_t$ と無視できるが、**その分だけ一様コアが出口面に届かない**
    (実測: コア半径が出口半径の 80%)。長さの定義として不正確。

    場の外へ出た区間は一様出口状態 $(\theta=0, M=M_d)$ で凍結して直線延長する
    (終端特性線の下流側は定義上一様域なので近似ではない)。戻り値に凍結ステップ数を
    含めるので、追跡の大半が凍結なら「場が短すぎる」と分かる。
    """
    itp = field_interpolator(pts) if itp is None else itp
    wx, wr = wall[:, 0], wall[:, 1]
    mu_d = float(np.arcsin(1.0 / max(M_d, 1.0 + 1e-12)))
    x, r, n_frozen = float(x_E), 0.0, 0
    path = [(x, r)]
    for _ in range(n_max):
        th, nu = (float(v) for v in itp(x, r))
        if np.isfinite(th) and np.isfinite(nu):
            M = pm_mach(max(nu, 1e-9), gamma)
            slope = np.tan(th + np.arcsin(1.0 / max(M, 1.0 + 1e-12)))
        else:
            slope = np.tan(mu_d)
            n_frozen += 1
        x_n, r_n = x + ds, r + ds * slope
        rw, rw_n = float(np.interp(x, wx, wr)), float(np.interp(x_n, wx, wr))
        if r_n >= rw_n:                       # 壁を越えた — 線形補間で交点
            t = (rw - r) / max((r_n - r) - (rw_n - rw), 1e-30)
            x_F = x + t * ds
            return {"x_F": float(x_F), "r_F": float(np.interp(x_F, wx, wr)),
                    "ok": True, "n_frozen": n_frozen, "path": np.asarray(path)}
        x, r = x_n, r_n
        path.append((x, r))
        if x > wx[-1]:
            break
    return {"x_F": float("nan"), "r_F": float("nan"), "ok": False,
            "n_frozen": n_frozen, "path": np.asarray(path)}


def _design_cplus(inv, init, n_ax, ax, mdot_star, g, x0, x_wall0,
                  exit_mode, x_E, M_d, start_line, axis_thr=None) -> dict:
    r"""`wall_mode='cplus'` の設計本体 — **補間構造を一切作らない**経路。

    壁は各 C⁺ 線上の流束閉包 (`cplus_flux_wall`)、物理出口 $F$ は $E=(x_E,0)$ を
    起点とする**網の C⁺ 線**と壁の交点。どちらも網の点だけで決まるので Delaunay も
    流線 ODE も不要 (設計コストの 82% を占めていた三角形分割が丸ごと消える)。

    診断は `mdot_ratio_moc` を返さない — 本閉包では構成的に 1 になる**循環指標**に
    なるため (A9 の教訓)。代わりに**流線整合残差** $\max|dr/dx-\tan\theta_{\rm net}|$
    を返す: 壁は流線でもあるはずなので、独立な 2 つの閉包の食い違いを測っている。

    `axis_thr`: `axis_limit='analytic'` のときの初期前線の θ_r (`axis_theta_r_init`)。
    """
    lev = inv.fill_levels(init, axis_thr=axis_thr)
    cum0 = np.zeros(len(init))
    cum0[n_ax:] = _flux_along(init[n_ax:], g)
    wall = cplus_flux_wall(lev, cum0, mdot_star, g)
    if len(wall) < 10:
        raise RuntimeError(f"cplus: 壁点が {len(wall)} 個しか取れない (場が不足)")
    wall_full = wall.copy()
    exit_info: dict = {"mode": exit_mode}
    if exit_mode == "characteristic":
        if x_E is None or M_d is None:
            raise ValueError("exit_mode='characteristic' には x_E と M_d が必要")
        m_E = int(np.argmin(np.abs(np.asarray(ax) - float(x_E))))
        term = cplus_lines(lev, m_E)
        if len(term) < 3:
            raise RuntimeError("cplus: 終端 C⁺ が短すぎる")
        # 終端 C⁺ と壁の交点 (どちらも網由来の折れ線)
        rw_on_term = np.interp(term[:, 0], wall[:, 0], wall[:, 1],
                               left=np.nan, right=np.nan)
        d = term[:, 1] - rw_on_term
        j = int(np.argmax(d >= 0.0)) if np.any(d >= 0.0) else -1
        if j <= 0:
            raise RuntimeError("終端特性線が壁に到達しない (軸 target の延長不足)")
        s = -d[j - 1] / max(d[j] - d[j - 1], 1e-30)
        x_F = float(term[j - 1, 0] + s * (term[j, 0] - term[j - 1, 0]))
        exit_info.update({"x_F": x_F,
                          "r_F": float(np.interp(x_F, wall[:, 0], wall[:, 1])),
                          "x_E_index": m_E, "n_frozen": 0})
        wall = wall[wall[:, 0] <= x_F]
    elif exit_mode == "lip":
        i_pk = int(np.argmax(wall[:, 2]))
        after = wall[i_pk:, 2] <= np.deg2rad(0.05)
        if np.any(after):
            wall = wall[: i_pk + int(np.argmax(after)) + 1]
        exit_info.update({"x_F": float(wall[-1, 0]), "r_F": float(wall[-1, 1])})
    else:
        raise ValueError("exit_mode は 'characteristic' か 'lip'")
    # 独立診断: 壁は流線でもあるはず (循環しない整合チェック)
    m = wall[:, 0] > wall[0, 0] + 0.5
    res = (np.gradient(wall[m, 1], wall[m, 0]) - np.tan(wall[m, 2])) if m.sum() > 5 \
        else np.array([np.nan])
    return {"wall": wall, "wall_full": wall_full, "pts": None, "x0": float(x0),
            "start_line": start_line, "x_wall0": float(x_wall0),
            "wall_mode": {"mode": "cplus", "n_wall": int(len(wall)),
                          "streamline_residual_max": float(np.max(np.abs(res))),
                          "streamline_residual_rms": float(np.sqrt(np.mean(res ** 2)))},
            "levels": lev, "exit": exit_info,
            "mdot_start": mdot_star, "mdot_exit": float("nan")}


def _axis_grid(x0: float, x_end: float, n: int, dx0: float | None,
               q: float = 1.05) -> np.ndarray:
    r"""軸点列。`dx0=None` は等間隔 (従来)。`dx0` を与えると**スロート側だけ**
    初項 dx0・公比 q の等比で細分し、間隔が dx_max に達したら等間隔に切り替える
    (dx_max は総長が x_end−x0 になるよう解く。全体を等比にすると下流が粗くなり
    出口平坦部の壁に 1e-5 のジッタが出るので、細分はスロート側に限る)。

    背景 (2026-08-16, M6): 初期値線 (スロート特性線) と壁の間の場は、初期値線点が
    C⁻ 担体として退化する (線自体が C⁻) ため、**軸点から後退する C⁻ だけで**埋まる。
    スロート直後の壁点は最初の数本の C⁻ でしか決まらないので、軸間隔が粗い
    (M6 は場が x_end≈180 $r_t$ に伸び、n_axis=1200 でも dx=0.15) と壁の最初の
    ~0.3 $r_t$ が ±2e-4 $r_t$・±0.2° でジッタし (円弧からの乖離)、壁 QA
    (単調 / spline リンギング) を落とす。dx₀≈0.03 で M4.2 (dx 0.1) 相当の滑らかさ。"""
    if dx0 is None or dx0 <= 0.0:
        return np.linspace(x0, x_end, n)
    L = float(x_end - x0)
    if dx0 * (n - 1) >= L:                    # 等間隔でも dx0 より細かい
        return np.linspace(x0, x_end, n)

    def spacing(dx_max):
        k = int(np.floor(np.log(dx_max / dx0) / np.log(q))) + 1   # 等比区間の点数
        k = max(1, min(k, n - 1))
        d = dx0 * q ** np.arange(k)
        d = np.minimum(d, dx_max)
        rest = n - 1 - k
        return np.concatenate([d, np.full(rest, dx_max)])
    a, b = dx0, L                             # 総長は dx_max について単調増加
    for _ in range(200):
        m = 0.5 * (a + b)
        if spacing(m).sum() < L:
            a = m
        else:
            b = m
        if b - a < 1e-13 * L:
            break
    d = spacing(0.5 * (a + b))
    x = np.concatenate([[x0], x0 + np.cumsum(d)])
    x *= 1.0
    x[-1] = x_end
    return x


# 軸端の接続検査 (target の軸節点と throat の初期線の軸端) の許容差。どちらも同じ量を別経路で作っているので
# 構成上は丸め誤差で一致するはずの量 (plan §4.1)。M・ν は絶対値、M′ は max(1, |M′|) に対する相対値
CONN_TOL = 1e-9


def axis_theta_r_init(init, n_ax: int, target, g, target_dM=None, axis_anchor=None,
                      M_line_axis: float | None = None):
    r"""初期前線の軸端点の解析極限 $\theta_r$ (plans/active/discretization-moc-axis-limit-and-corrector.md §4.1)。

    - 軸節点 `init[:n_ax]` ($r=0$): $M$ = `target(x)`、$M'$ = `target_dM(x)` (軸則の解析微分) →
      $\theta_r=\frac12\sqrt{M^2-1}\,\nu_M(M)M'$ (`moc_kernel.axis_theta_r`)。
    - 初期線の軸端 `init[n_ax]` ($r\le$ `AXIS_R_EPS` のとき): アンカー `axis_anchor` = (x_A, M_A, M′_A) が
      その x にあれば (M_A, M′_A)。一般経路 ($x_A\ne x_0$、`runner_axismach` の CFD 反復アンカー) では
      初期線の軸端を $x_A$ と扱わず、$M$ = 初期線の軸端の M (`M_line_axis`)、$M'$ = `target_dM(x_0)`。
    - `target_dM` が無い (軸則に微分が無い) ときの代わり: 軸端点の ν を x の 3 次スプラインにして
      $\theta_r=\frac12\sqrt{M^2-1}\,d\nu/dx$ ($\nu_M$ を掛けない)。
    戻り: (thr [len(init)、軸端点以外 NaN], info)。info: `source` (軸節点・初期線の軸端それぞれの出所)、
    `connection` (target と throat の軸端の x・M・ν・M′ の差と `ok`)。"""
    xs = np.array([p.x for p in init[:n_ax]], dtype=float)
    M_ax = np.array([float(target(float(x))) for x in xs], dtype=float)
    p_end = init[n_ax] if len(init) > n_ax else None
    on_end = p_end is not None and p_end.r <= AXIS_R_EPS
    M_l = (float(M_line_axis) if M_line_axis is not None else float(p_end.M)) if p_end is not None else None
    thr = np.full(len(init), np.nan)
    src: dict = {}
    spl = None
    if target_dM is None:
        from scipy.interpolate import CubicSpline
        xx = np.r_[[p_end.x] if on_end else [], xs[::-1]]
        nn = np.r_[[float(pm_nu(M_l, g))] if on_end else [], np.array([p.nu for p in init[:n_ax]])[::-1]]
        spl = CubicSpline(xx, nn).derivative()
        thr[:n_ax] = 0.5 * np.sqrt(np.maximum(M_ax ** 2 - 1.0, 0.0)) * spl(xs)
        src["axis_nodes"] = "nu_spline"
    else:
        thr[:n_ax] = axis_theta_r(M_ax, np.array([float(target_dM(float(x))) for x in xs]), g)
        src["axis_nodes"] = "law_derivative"
    conn: dict = {"x_line_axis": None if p_end is None else float(p_end.x),
                  "r_line_axis": None if p_end is None else float(p_end.r), "tol": CONN_TOL}
    if on_end:
        x_l = float(p_end.x)
        M_t = float(target(x_l))
        Mp_t = None if target_dM is None else float(target_dM(x_l))
        conn.update(x_axis_node_first=(float(xs[-1]) if n_ax else None), M_line=M_l, M_target=M_t,
                    dM=M_l - M_t, dnu=float(pm_nu(M_l, g)) - float(pm_nu(M_t, g)), Mp_target=Mp_t)
        okc = abs(conn["dM"]) <= CONN_TOL and abs(conn["dnu"]) <= CONN_TOL
        same_x = False
        if axis_anchor is not None:
            x_A, M_A, Mp_A = (float(v) for v in axis_anchor)
            same_x = abs(x_l - x_A) <= AXIS_R_EPS
            conn.update(x_A=x_A, dx_line_minus_xA=x_l - x_A, line_axis_is_x_A=bool(same_x))
            if same_x:
                conn.update(dM_anchor=M_A - M_l)
                okc = okc and abs(M_A - M_l) <= CONN_TOL
                if Mp_t is not None:
                    conn.update(dMp_anchor_minus_target=Mp_A - Mp_t)
                    okc = okc and abs(Mp_A - Mp_t) <= CONN_TOL * max(1.0, abs(Mp_A))
        conn["ok"] = bool(okc)
        if same_x:
            thr[n_ax] = float(axis_theta_r(M_A, Mp_A, g))
            src["line_axis_end"] = "anchor"
            conn["Mp_used"] = Mp_A
        elif target_dM is not None:
            thr[n_ax] = float(axis_theta_r(M_l, Mp_t, g))
            src["line_axis_end"] = "target_derivative@x_line"
            conn["Mp_used"] = Mp_t
        else:
            thr[n_ax] = 0.5 * np.sqrt(max(M_l * M_l - 1.0, 0.0)) * float(spl(x_l))
            src["line_axis_end"] = "nu_spline"
            conn["Mp_used"] = None
    else:
        conn["ok"] = True
        src["line_axis_end"] = None            # 初期線の始点が軸上にない (θ_r 不要)
    if not np.all(np.isfinite(thr[:n_ax])) or (on_end and not np.isfinite(thr[n_ax])):
        raise ValueError("axis_theta_r_init: 軸端点の θ_r が非有限 (軸則・アンカーの微分を確認)")
    info = {"source": src, "connection": conn,
            "first": {"x": float(xs[-1]) if n_ax else None, "theta_r": float(thr[n_ax - 1]) if n_ax else None},
            "line_axis_end": {"x": conn["x_line_axis"], "theta_r": float(thr[n_ax]) if on_end else None}}
    return thr, info


def moc_gate(diag: dict, thr_info: dict | None, wall_full=None) -> dict:
    r"""単位過程のゲート (plan §4.2)。**設計は止めない** (診断に記録し、検証・生産のゲートが読む)。

    `corrector='converge'` のときだけ合否を出す (fixed2 は収束を判定しないので `pass` = None):
    - 反復の失敗 (反復中の NaN・Inf、上限到達) が 1 対でもあれば不合格
    - 最終状態の残差 (幾何の交点式・適合式) の最大が `RESID_TOL` を超えたら不合格
    - 幾何的棄却 (平行な特性線・軸より下) が許す領域の外で起きたら不合格。許す領域は従来どおり
      「壁の外」(対の両端が設計壁 `wall_full` [x, r] 以上の半径にある、または壁の x 範囲の外 = 網の端)
    - `axis_limit='analytic'` で軸端の接続 (target と throat の x・M・ν・M′) が一致しなければ不合格"""
    reasons = []
    applicable = diag["corrector"] == "converge"
    n_fail = diag["pairs"]["iter_nonfinite"] + diag["pairs"]["iter_maxiter"]
    rmax = None
    if applicable:
        if n_fail:
            reasons.append(f"反復の失敗 {n_fail} 対 (非有限 {diag['pairs']['iter_nonfinite']}・"
                           f"上限到達 {diag['pairs']['iter_maxiter']})")
        rmax = max(diag["resid"]["geom_max"], diag["resid"]["comp_max"])
        if not rmax <= RESID_TOL:
            reasons.append(f"最終残差 {rmax:.3e} > {RESID_TOL:g}")
    n_geom = diag["pairs"]["geom_parallel"] + diag["pairs"]["geom_below_axis"]
    inside = None
    if n_geom:
        inside = 0
        for (_lev, _i, _kind, xa, ra, xb, rb) in diag["_geom_all"]:
            if wall_full is None or len(wall_full) < 2:
                inside += 1
                continue
            wx, wr = wall_full[:, 0], wall_full[:, 1]
            out_x = (min(xa, xb) > wx[-1]) or (max(xa, xb) < wx[0])
            above = (ra >= np.interp(xa, wx, wr)) and (rb >= np.interp(xb, wx, wr))
            if not (out_x or above):
                inside += 1
        if inside:
            reasons.append(f"幾何的棄却 {inside} 対が壁の内側 (許す領域 = 壁の外・網の端の外)")
    conn_ok = None
    if diag["axis_limit"] == "analytic" and thr_info is not None:
        conn_ok = bool(thr_info["connection"]["ok"])
        if not conn_ok:
            reasons.append("軸端の接続不一致 (target と throat の x・M・ν・M′)")
    return {"applicable": applicable, "pass": (not reasons) if applicable else None, "reasons": reasons,
            "n_iter_fail": int(n_fail), "resid_max": rmax, "n_geom_reject": int(n_geom),
            "n_geom_reject_inside_wall": inside, "connection_ok": conn_ok}


def inverse_design(throat, target, x_axis_end: float, n_axis: int = 260,
                   n_start: int = 41, gamma: float = 1.4, dx_wall: float = 0.02,
                   th_wall0: float | None = None, M_start: float = 1.05,
                   exit_mode: str = "lip", x_E: float | None = None,
                   M_d: float | None = None, start_line: str = "vertical",
                   wall_mode: str = "streamline", blend_width: float = 1.0,
                   axis_dx0: float | None = None, axis_limit: str = "legacy",
                   corrector: str = "fixed2", target_dM=None, axis_anchor=None):
    """starting line (throat: SauerThroat) + 軸目標 target(x) から壁を逆設計する。

    target: 呼び出し可能 M(x) — starting line の軸点 x0 で場に C1 整合していること
    (モード F: MachBezier.from_constraints の start に Sauer 軸微分を渡す)。
    x_axis_end: 軸目標の下流端 (壁端の C⁻ 足より下流まで — 一様出口なら
    M_d 一定を伸ばすだけ)。
    axis_limit / corrector: 単位過程の選択肢 (`moc_kernel.interior_vec`。既定 legacy + fixed2 は従来とビット同一。
    plans/active/discretization-moc-axis-limit-and-corrector.md)。`analytic` では軸端点の θ_r を
    `axis_theta_r_init` で作る — `target_dM` (軸則の dM/dx。無ければ ν の 3 次スプライン微分で代用) と
    `axis_anchor` = (x_A, M_A, M′_A) (初期線の軸端が x_A にあるとき、その θ_r に使う)。
    戻り値 dict: wall (n,4 [x,r,θ,M]), pts, mdot_start, mdot_exit, moc (単位過程の集計・θ_r の出所・ゲート)。
    """
    # wall_mode: 'cplus' = 壁点を C⁺ 線上の流束閉包で決める古典法 (A10, `_design_cplus`。
    # Delaunay も流線 ODE も通らない) / 'streamline' = 三角充填 + 補間場の流線積分 (旧既定。
    # 向き非依存 fill への「修正」は出口指標を悪化させ撤回 — plan §9.1)。
    # streamline の既知の限界: 壁足の曲率が円弧 1/R でなく ~0 に寝る。
    inv = InverseMOC(gamma=gamma, delta=1.0, axis_limit=axis_limit, corrector=corrector)
    g = gamma
    if start_line == "throat_char":
        # **スロート特性線を初期値線にする** (CONTUR 流, 2026-08-15 ユーザ指摘)。
        # 線自体が C⁻ なので担体割当てが区間内で反転せず、縦線構成で残っていた
        # スロート直後の未計算楔が構造的に消える。壁流線はスロート壁点から始まる。
        xs_c, rr, MM, tt = throat.throat_characteristic(n=n_start)
        x0 = float(xs_c[0])            # 軸着地 = 設計が軸に効き始める位置
        x_line = np.asarray(xs_c, dtype=float)
    elif start_line == "vertical":
        x0, rr, MM, tt = throat.starting_line(M_start=M_start, n=n_start)
        x_line = np.full(n_start, float(x0))
    else:
        raise ValueError("start_line は 'throat_char' か 'vertical'")
    ax = _axis_grid(x0, x_axis_end, n_axis, axis_dx0)[1:][::-1]
    init = [_Pt(float(x), 0.0, 0.0, float(pm_nu(float(target(x)), g)), g) for x in ax]
    init += [_Pt(float(x_line[i]), float(rr[i]), float(tt[i]),
                 float(pm_nu(float(MM[i]), g)), g) for i in range(n_start)]
    if th_wall0 is not None and start_line == "vertical":
        # 壁足の θ を厳密壁接線で上書き (Sauer 線形化 vbar は過小 — kernel と同じ補正)
        # throat_char では壁足 = 幾何スロートで θ=0 が厳密に成り立つので不要。
        init[-1].th = float(th_wall0)
    mdot_star = float(_flux_along(init[len(ax):], g)[-1])
    axis_thr, thr_info = None, None
    if axis_limit == "analytic":
        axis_thr, thr_info = axis_theta_r_init(init, len(ax), target, g, target_dM=target_dM,
                                               axis_anchor=axis_anchor, M_line_axis=float(MM[0]))
    if wall_mode == "cplus":
        res = _design_cplus(inv, init, len(ax), ax, mdot_star, g, x0,
                            float(x_line[-1]), exit_mode, x_E, M_d, start_line, axis_thr=axis_thr)
        res["moc"] = _moc_record(inv, thr_info, res["wall_full"])
        return res
    pts = inv.fill(init, axis_thr=axis_thr)
    # Delaunay は充填ベクトル化後の支配コスト → 1 個だけ作って全診断で使い回す
    itp = field_interpolator(pts)
    wall = inv.wall_streamline(pts, float(x_line[-1]), float(rr[-1]), dx=dx_wall,
                               itp=itp)
    wall_diag: dict = {"mode": wall_mode}
    if wall_mode == "flux":
        # **壁を流束閉包で置き換える** (A9)。断面が計算領域に入る x0 より下流だけ
        # 置換し、[x0, x0+blend_width] で流線壁と smoothstep 合成する
        # (r(x_wall0)=r_t を厳密に保つため — 合成帯では両者の差はまだ小さい)。
        x_lo = float(x0) + 0.05
        m = wall[:, 0] >= x_lo
        if int(m.sum()) < 10:
            raise RuntimeError("wall_mode='flux': 置換区間が短すぎる")
        xs_f = wall[m, 0]
        r_f = flux_closure_wall(itp, mdot_star, xs_f, 1.08 * wall[m, 1], gamma=g)
        if not np.all(np.isfinite(r_f)):
            raise RuntimeError("wall_mode='flux': 流束閉包が解けない断面がある")
        r_s = wall[m, 1].copy()
        w = _smoothstep((xs_f - x_lo) / max(blend_width, 1e-12))
        wall = wall.copy()
        wall[m, 1] = (1.0 - w) * r_s + w * r_f
        v = itp(wall[:, 0], wall[:, 1])          # θ, M は新しい壁上で取り直す
        okv = np.isfinite(v[:, 0]) & np.isfinite(v[:, 1])
        wall[okv, 2] = v[okv, 0]
        wall[okv, 3] = pm_mach_vec(v[okv, 1], g)
        wall_diag.update({
            "x_lo": x_lo, "blend_width": float(blend_width),
            # 合成帯で両者がどれだけ違うか (小さいほど合成の恣意性が小さい)
            "d_blend_max": float(np.max(np.abs((r_f - r_s)[w < 1.0]))) if np.any(w < 1.0) else 0.0,
            # 置換で壁がどれだけ動いたか (= 流線積分が溜め込んでいた誤差)
            "dr_max": float(np.max(np.abs(r_f - r_s))),
            "dr_exit": float((r_f - r_s)[-1])})
    elif wall_mode != "streamline":
        raise ValueError("wall_mode は 'streamline' か 'flux'")
    wall_full = wall.copy()
    exit_info: dict = {"mode": exit_mode}
    if exit_mode == "characteristic":
        # **物理出口 = 終端 C⁺ (E 発) と壁流線の交点** (正しい定義, `terminal_exit`)
        if x_E is None or M_d is None:
            raise ValueError("exit_mode='characteristic' には x_E と M_d が必要")
        te = terminal_exit(pts, wall, float(x_E), float(M_d), gamma=gamma, itp=itp)
        if not te["ok"]:
            raise RuntimeError("終端特性線が壁に到達しない (軸 target の延長不足の疑い"
                               f" — x_axis_end={x_axis_end:.3g}, wall_end={wall[-1,0]:.3g})")
        i_cut = int(np.searchsorted(wall[:, 0], te["x_F"]))
        wall = wall[:max(i_cut, 10)]
        exit_info.update({k: te[k] for k in ("x_F", "r_F", "n_frozen")})
        exit_info["term_path"] = te["path"]
    elif exit_mode == "lip":
        # [旧] リップ (θ が最大値を経て 0.05° を切る点) で切る。壁角は F へ漸近的に
        # 落ちるため、この閾値は F の数 r_t 上流で発火する (`terminal_exit` docstring)。
        if len(wall) > 10:
            i_pk = int(np.argmax(wall[:, 2]))
            after = wall[i_pk:, 2] <= np.deg2rad(0.05)
            if np.any(after):
                wall = wall[: i_pk + int(np.argmax(after)) + 1]
        exit_info.update({"x_F": float(wall[-1, 0]), "r_F": float(wall[-1, 1])})
    else:
        raise ValueError("exit_mode は 'characteristic' か 'lip'")
    # 質量流量診断: 壁 (=流管であるべき曲線) までの断面積分を ṁ* と比較
    mdot1 = (inv.massflux_across(pts, wall[-1, 0] - 0.3,
                                 float(np.interp(wall[-1, 0] - 0.3, wall[:, 0], wall[:, 1])),
                                 itp=itp)
             if len(wall) > 10 else float("nan"))
    return {"wall": wall, "wall_full": wall_full, "pts": pts, "x0": float(x0),
            "start_line": start_line, "x_wall0": float(x_line[-1]),
            "wall_mode": wall_diag,
            "exit": exit_info, "mdot_start": mdot_star, "mdot_exit": mdot1,
            "moc": _moc_record(inv, thr_info, wall_full)}


def _moc_record(inv, thr_info: dict | None, wall_full) -> dict:
    """`inverse_design` の戻り値 `moc`: 単位過程の集計 (`InverseMOC.last_diag`)・θ_r の出所と接続検査・ゲート。"""
    diag = dict(inv.last_diag)
    gate = moc_gate(diag, thr_info, wall_full)
    diag.pop("_geom_all", None)
    return {**diag, "theta_r": thr_info, "gate": gate}
