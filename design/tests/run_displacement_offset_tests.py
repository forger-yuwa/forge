#!/usr/bin/env python3
"""③ベル (M) の δ* 補正のずれ `geometry.wall.DisplacementOffset` と `RaoBellWall(offset=...)` の試験
(plan tooling-nozzle-bell-two-mode §5.1 #7 (7b-4))。

合格条件: (1) 立ち上げの始点 a で δ・δ′・δ″ が 0 (a = 0 ならスロートの r・r′・r″ を保つ)、a 以下は 0、立ち上げの終点 b 以降は
平滑化スプラインそのもの、(2) 1・2 階が数値微分と相対 1e-5 で一致 (積の微分)、(3) 終点の値が与えた値に 1e-10 で一致 (終点の重み)、
(4) 保存・読み込みで値がビット一致、(5) x_end を 1e-9 越えた点で 0 に落ちない (終点の値で頭打ち)、(6) `RaoBellWall` にずれを
足すと拡大部の r が非粘性壁 + δ、収縮部は不変、スロートの接合 (r・r′・r″) が非粘性壁と同じ。
実行: design/.venv-opt の Python で `python design/tests/run_displacement_offset_tests.py`。
"""
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from forge_design.geometry.wall import DisplacementOffset  # noqa: E402

FAIL = 0


def check(name, ok, detail=""):
    global FAIL
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAIL += 1


def synth(ramp):
    x = np.r_[0.0, np.sort(np.random.default_rng(0).uniform(0.0, 20.0, 600)), 20.0]
    x = np.unique(x)
    d = 1e-3 + 6e-3 * x + 1e-4 * x * x + 2e-7 * np.random.default_rng(1).standard_normal(len(x))
    return x, d, DisplacementOffset.fit(x, d, ramp=ramp, meta={"eps_inv": 48.0})


def main():
    for ramp in ((0.0, 0.41), (1.17, 2.0)):
        x, d, off = synth(ramp)
        a, b = ramp
        z = np.array([a])
        v = [float(off(z, k)[0]) for k in range(3)]
        check(f"ramp {ramp}: a で δ・δ′・δ″ = 0", max(abs(t) for t in v) < 1e-14, str(v))
        xl = np.linspace(0.0, a, 50) if a > 0 else np.array([0.0])
        check(f"ramp {ramp}: x ≤ a で 0", float(np.max(np.abs(off(xl)))) == 0.0)
        xr = np.linspace(b, 20.0, 200)
        check(f"ramp {ramp}: x ≥ b でスプライン", float(np.max(np.abs(off(xr) - off.spl(xr)))) < 1e-15)
        xs = np.linspace(0.05, 19.95, 400); h = 1e-5
        d1n = (off(xs + h) - off(xs - h)) / (2 * h); d2n = (off(xs + h) - 2 * off(xs) + off(xs - h)) / (h * h)
        e1 = float(np.max(np.abs(off(xs, 1) - d1n)) / np.max(np.abs(d1n)))
        e2 = float(np.max(np.abs(off(xs, 2) - d2n)) / np.max(np.abs(off(xs, 2))))
        check(f"ramp {ramp}: 1 階・2 階が数値微分と一致", e1 < 1e-5 and e2 < 1e-3, f"(1 階 {e1:.1e}・2 階 {e2:.1e})")
        check(f"ramp {ramp}: 終点の値", abs(float(off(np.array([20.0]))[0]) - d[-1]) < 1e-10,
              f"({float(off(np.array([20.0]))[0]) - d[-1]:.1e})")
        check(f"ramp {ramp}: 終点の 1e-9 先で頭打ち", abs(float(off(np.array([20.0 + 1.5e-9]))[0]) - d[-1]) < 1e-10)
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "o.npz"; off.save(p); o2 = DisplacementOffset.load(p)
            xx = np.linspace(0, 20, 1001)
            check(f"ramp {ramp}: 保存・読み込みでビット一致",
                  all(np.array_equal(off(xx, k), o2(xx, k)) for k in range(3)) and o2.meta == off.meta and o2.ramp == off.ramp)
    # RaoBellWall に足す (L100 のピン版の壁があれば)
    root = Path(__file__).resolve().parents[2]
    wf = root / "case/68.bell_two_mode/rao_m/rao_pin_parts/L100_pin_n321_wall.npz"
    if wf.exists():
        from forge_design.geometry.rao_axisym import load_rao_wall
        from forge_design.geometry.wall import RaoBellWall
        from forge_design.geometry.wall_bell_c2 import BellC2Spec
        rw, _ = load_rao_wall(wf)
        spec = BellC2Spec(eps=float(rw.r_end) ** 2, R_t=1.0, contraction_ratio=3.0, L_U=3.0, L_pipe=1.0)
        xg = np.r_[0.0, np.linspace(1e-3, float(rw.x_end), 300)]
        xg[-1] = float(rw.x_end)
        off = DisplacementOffset.fit(xg, 1e-3 + 5e-3 * xg, ramp=(0.0, 0.41))
        w0, w1 = RaoBellWall(rw, spec), RaoBellWall(rw, spec, offset=off)
        xe = np.linspace(0.0, w0.x_e, 2001)
        check("RaoBellWall: 拡大部 = 非粘性壁 + δ", all(float(np.max(np.abs(w1.r(xe, k) - w0.r(xe, k) - off(xe, k)))) < 1e-14 for k in range(3)))
        xc = np.linspace(w0.x_in, -1e-6, 500)
        check("RaoBellWall: 収縮部・直管は不変", all(np.array_equal(w1.r(xc, k), w0.r(xc, k)) for k in range(3)))
        tj0, tj1 = w0.throat_jump(), w1.throat_jump()
        check("RaoBellWall: スロートの接合が同じ", all(abs(tj1[k] - tj0[k]) < 1e-14 for k in (0, 1, 2)), str(tj1))
        check("RaoBellWall: validate", w1.validate() == [], str(w1.validate()))
    else:
        print(f"[SKIP] RaoBellWall の試験 ({wf} が無い)")
    print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
