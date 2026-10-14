"""plan tooling-nozzle-core-grid §4.98 P-hi の準備: 入口半径を半分にした形 H の格子・変換・初期値 (AWS で回す)。

形: 配管 r = R_U_NEW (x ∈ [x_in, X_U_NEW])、縮流部 = (X_U_NEW, R_U_NEW, 0, 0) → x = 0 の今の物理の壁の (r, r′, r″) の 5 次 Hermite、x ≥ 0 は今の物理の壁。
`prepare_ns` (ns_n012._prepare_ns、Gref の準備と同じ引数) の中で次の 2 つの関数だけを差し替える:
  - `runner_axismach.build_physical_wall` → 元の壁を組んでから縮流部だけを差し替えた壁 (HalfInletWall) を返す。
  - `mesh2d._first_frac_profile` → 元の比 × r_今(x)/r_新(x) (旧格子の実寸の第一層厚を再現する)。
問題 YAML は Gref のものから name と physical_wall_repr だけを変えたもの (差し替えた壁を single_bspline の壁ファイルに書かせないため)。
初期値は run_0608 (Gref の B0、最後の res は res_60000) から interp_field (prepare_ns の ic_from)。差し替えた壁の記録は <out>/half_inlet_wall.json。
usage (AWS の case/45 で): FORGE_BIN=... REAL_CONVERTER=<fgeom7 の FP64 の変換器> FORGE_CONVERTER=$(pwd)/conv_tolerant.sh python3 cg_hi_prep.py <out_dir>"""
import hashlib, json, os, sys
from pathlib import Path
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_pair as CP   # noqa: E402
import ns_n012 as NS     # noqa: E402
RA = NS.runner()
from forge_design.meshing import mesh2d  # noqa: E402
from forge_design.geometry.wall_walldriven import _hermite_quintic, _poly_eval  # noqa: E402

R_U_OLD, R_U_NEW, X_U_NEW = 6.485, 6.485 / 2, -9.0
SRC_PROBLEM = HERE / "problem_d155_ns_prod_coldmesh_tw300_cgref.yaml"
H_PROBLEM = HERE / "problem_d155_ns_prod_coldmesh_tw300_cgref_hi.yaml"
IC_RUN, IC_RES = "run_0608_cg_gref_cont_k1", "res_60000.h5"
STATE = {}


class HalfInletWall:
    """元の物理の壁の x < 0 だけを差し替えた壁。他の属性は元の壁に任せる。"""

    def __init__(self, base):
        self._b = base
        # required_attrs は wall.__dict__ から pw_upstream を直接読む (wall_axismach.required_attrs) ので、インスタンスの属性として写す
        self.pw_upstream = base.__dict__.get("pw_upstream")
        j = [float(np.asarray(base.r(np.array([0.0]), d) if d else base.r(np.array([0.0])))[0]) for d in (0, 1, 2)]
        self.join = j
        self._L = -X_U_NEW
        self._c = _hermite_quintic(X_U_NEW, 0.0, R_U_NEW, 0.0, 0.0, j[0], j[1], j[2])

    def __getattr__(self, k):
        if k.startswith("__"):
            raise AttributeError(k)
        return getattr(self._b, k)

    def r(self, x, deriv: int = 0):
        xa = np.asarray(x, dtype=float); scalar = xa.ndim == 0; xa = np.atleast_1d(xa)
        out = np.array(self._b.r(xa, deriv) if deriv else self._b.r(xa), dtype=float, copy=True)
        m = xa < 0.0
        if np.any(m):
            xi = (np.clip(xa[m], X_U_NEW, 0.0) - X_U_NEW) / self._L
            v = _poly_eval(self._c, xi, deriv) / self._L ** deriv
            v = np.where(xa[m] < X_U_NEW, (R_U_NEW if deriv == 0 else 0.0), v)
            out[m] = v
        return float(out[0]) if scalar else out

    def own_checks(self):
        x = np.linspace(float(self._b.x_in), 0.0, 125001)
        r, rp = self.r(x), self.r(x, 1)
        msgs = []
        if np.any(rp[x < -0.01] > 0.0):
            msgs.append(f"縮流部で dr/dx > 0 (最大 {rp[x < -0.01].max():.3e})")
        if not (abs(float(r.max()) - R_U_NEW) < 1e-12 and float(r.min()) > 1.0):
            msgs.append(f"r の範囲が想定外 ({r.min():.6f}..{r.max():.6f})")
        e = 1e-9
        for d in (0, 1, 2):
            a, b = float(self.r(np.array([-e]), d)[0]), float(self.r(np.array([0.0]), d)[0])
            if abs(a - b) > 1e-6 * max(1.0, abs(b)):
                msgs.append(f"x = 0 での {d} 階の不連続 ({a} 対 {b})")
        return msgs

    def validate(self, n: int = 4000):
        return list(self._b.validate()) + self.own_checks()


_orig_bpw = RA.build_physical_wall
def _bpw(*a, **k):
    base = _orig_bpw(*a, **k)
    w = HalfInletWall(base); STATE["base"], STATE["wall"] = base, w
    return w


_orig_ff = mesh2d._first_frac_profile
def _ff(xs, prm):
    f = _orig_ff(xs, prm)
    b, w = STATE["base"], STATE["wall"]
    xs = np.asarray(xs, dtype=float)
    return f * np.asarray(b.r(xs), dtype=float) / np.asarray(w.r(xs), dtype=float)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    out = Path(sys.argv[1]).resolve()
    if out.exists():
        raise SystemExit(f"{out} が既にある — 止める")
    NS.check_dry_env(False)
    conv = Path(os.environ.get("REAL_CONVERTER", ""))
    if not conv.is_file():
        raise SystemExit("REAL_CONVERTER が無い — 止める")
    # 問題 YAML: Gref のものから name と physical_wall_repr だけを変える
    txt = SRC_PROBLEM.read_text()
    lines = []
    for ln in txt.splitlines():
        if ln.startswith("name: "):
            ln = "name: problem_d155_ns_prod_coldmesh_tw300_cgref_hi"
        elif ln.strip().startswith("physical_wall_repr:"):
            continue
        lines.append(ln)
    hdr = ("# plan tooling-nozzle-core-grid §4.98 P-hi: cg_hi_prep.py が Gref の問題 YAML から作る (手で編集しない)。違いは name と physical_wall_repr の行だけ。\n"
           "# 入口半径の半減 (配管 r 3.2425・縮流部 x ∈ [−9, 0] の 5 次 Hermite) は cg_hi_prep.py が壁を組む関数を差し替えて入れる (この YAML の geometry は今の形のまま)。\n")
    H_PROBLEM.write_text(hdr + "\n".join(lines) + "\n")
    p0 = yaml.safe_load(SRC_PROBLEM.read_text()); p1 = yaml.safe_load(H_PROBLEM.read_text())
    p0.pop("name"); p1.pop("name"); p0["geometry"].pop("physical_wall_repr")
    if p0 != p1:
        raise SystemExit("問題 YAML の差が name と physical_wall_repr 以外にもある — 止める")
    src = HERE / IC_RUN
    rs = NS.res_files(src)
    if not rs or rs[-1].name != IC_RES:
        raise SystemExit(f"IC のドナー {src.name} の最後の res が {rs[-1].name if rs else None} ({IC_RES} であること) — 止める")
    dr_csv = HERE / CP.WALL_REF / "delta_r_initial.csv"
    RA.build_physical_wall = _bpw
    mesh2d._first_frac_profile = _ff
    info = NS._prepare_ns(H_PROBLEM, out, nsteps=12000, ic_from=src, delta_r_csv=str(dr_csv), offset="radial",
                          euler_ref=str(HERE / CP.EULER_REF), cfl_main=5.0, implicit_relax=CP.RELAX)
    NS.jdump(out / "prepare_info.json", info)
    w = STATE["wall"]
    rec = {"plan": "plans/active/tooling-nozzle-core-grid.md §4.98", "tool": "cg_hi_prep.py", "created": NS.now(),
           "shape": {"r_U_old": R_U_OLD, "r_U_new": R_U_NEW, "x_U_new": X_U_NEW, "x_in": float(w._b.x_in), "x_e": float(w._b.x_e),
                     "pipe_length": float(X_U_NEW - w._b.x_in), "join_at_x0": {"r": w.join[0], "r1": w.join[1], "r2": w.join[2]},
                     "mu": X_U_NEW ** 2 / (2.0 * (R_U_NEW - w.join[0])), "own_checks": w.own_checks()},
           "first_layer": "元の wall_first_frac_table の比 × r_今(x)/r_新(x) (旧格子の実寸)",
           "problem": H_PROBLEM.name, "problem_sha256": sha(H_PROBLEM), "src_problem_sha256": sha(SRC_PROBLEM),
           "converter": str(conv), "converter_sha256": sha(conv), "forge_bin": os.environ.get("FORGE_BIN"),
           "ic": {"src_run": IC_RUN, "src_res": IC_RES, "src_res_sha256": sha(rs[-1])},
           "delta_r_csv_sha256": sha(dr_csv), "nozzle_h5_sha256": sha(out / "nozzle.h5")}
    (out / "half_inlet_wall.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str))
    print(json.dumps(rec["shape"], ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
