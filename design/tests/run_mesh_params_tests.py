"""runner_axismach.mesh_params の試験 (plan tooling-nozzle-cfd-pinned-initial-line §5.1 #11f 手順 1, 2026-10-06)。
(a) 既存の Euler 問題 (格子キーが ni/nj/wall_first_frac/throat_refine だけ) で修正前の Euler 経路と同じ Mesh2DParams になる。
(b) NS の細分格子のキー (throat_width・wall_first_frac_throat・前後ブレンド) を Euler 経路でも読む。
(c) 同じ壁で Euler 経路と NS 経路の格子パラメータから生成した座標が一致する。
2026-10-07 (plan verification-case45-euler-total-enthalpy §6 E3) から Euler の prepare は mesh_euler を読む (`mesh_params_euler`;
試験は run_mesh_euler_tests.py)。ここの (a)〜(d) は mesh を読む `mesh_params` (NS の経路。E3 以前は Euler も同じ関数) の試験として残す。"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.evaluate.runner_axismach import mesh_params  # noqa: E402
from forge_design.meshing.mesh2d import Mesh2DParams  # noqa: E402
fails = 0


def check(name, ok, info=""):
    global fails
    print(("ok   " if ok else "FAIL ") + name + (f" ({info})" if info else "")); fails += (not ok)


class P:
    def __init__(self, mesh): self.mesh = mesh


old_euler = P({"ni": 1100, "nj": 65, "wall_first_frac": 5.0e-3, "throat_refine": 3.0})
ref = Mesh2DParams(ni=1100, nj=65, wall_first_frac=5.0e-3, throat_refine=3.0, scale=0.0768075)
check("(a) 既存 Euler 問題で修正前と同じ", mesh_params(old_euler, 0.0768075, 321, 65, 5.0e-3) == ref)
fine = {"ni": 2000, "nj": 97, "wall_first_frac": 1.3e-5, "wall_first_frac_throat": 4.5e-6, "wall_first_up_x0": -9.0, "wall_first_up_x1": -4.0,
        "wall_first_blend_x0": 1.0, "wall_first_blend_x1": 17.0, "throat_refine": 4.0, "throat_width": 3.0}
me = mesh_params(P(fine), 0.0767531, 321, 65, 5.0e-3); mn = mesh_params(P(fine), 0.0767531, 561, 97, 4.5e-5)
check("(b) Euler 経路が細分格子のキーを読む", (me.throat_width, me.wall_first_frac_throat, me.wall_first_up_x0, me.wall_first_blend_x1) == (3.0, 4.5e-6, -9.0, 17.0))
check("(c) 同じキーなら Euler 経路と NS 経路のパラメータが一致", me == mn)
check("(c') 既定値は dataclass と同じ (キー無しでビット不変)", mesh_params(P({}), 1.0, 241, 81, 0.002) == Mesh2DParams())
m = mesh_params(P(dict(fine, axis_cap_frac=0.02)), 0.0767531, 321, 65, 5.0e-3)
check("(d) axis_cap_frac が Euler/NS 両経路に渡る", m.axis_cap_frac == 0.02 and mesh_params(P(dict(fine, axis_cap_frac=0.02)), 0.0767531, 561, 97, 4.5e-5) == m)
from forge_design.meshing.mesh2d import _radial_fracs_capfixed  # noqa: E402
import numpy as np  # noqa: E402
s, ratio = _radial_fracs_capfixed(257, 1.3e-5, 0.02); d = np.diff(s)
check("(e) capfixed: 第一セル・上限・単調・和", abs(d[-1] - 1.3e-5) < 1e-12 and abs(d[0] - 0.02) < 1e-9 and np.all(d > 0) and abs(s[-1] - 1) < 1e-15 and ratio < 1.05, (d[-1], d[0], ratio))
for args, why in (((65, 1.3e-5, 0.0133), "和が 1 にならない"), ((257, 0.03, 0.02), "第一セル > 上限")):
    try:
        _radial_fracs_capfixed(*args); check(f"(f) 実現不能 {args} を拒否", False)
    except ValueError:
        check(f"(f) 実現不能 {args} を拒否 ({why})", True)
import ast  # noqa: E402
for mod in ("evaluate/runner_axismach.py", "feedback/deltastar_loop.py"):
    body = ast.parse((Path(__file__).resolve().parents[1] / "forge_design" / mod).read_text()).body
    mi = [i for i, n in enumerate(body) if isinstance(n, ast.If) and "__main__" in ast.unparse(n.test)]
    ld = max(i for i, n in enumerate(body) if isinstance(n, (ast.FunctionDef, ast.ClassDef)))
    check(f"(g) {mod}: `if __name__ == \"__main__\"` が全関数定義の後 (CLI から後方の関数に届く; codex result-2 2026-10-06)", bool(mi) and mi[0] > ld)
print(f"FAIL 件数: {fails}"); sys.exit(1 if fails else 0)
