import sys, numpy as np
from forge_design.report.nozzle_report import eta_line, load_field
run, res = sys.argv[1], sys.argv[2]
F = load_field(run, res); S = float(F["info"]["scale_m"])
for x0, x1 in ((4.0, 7.95),):
    xq = np.linspace(x0 / S, x1 / S, 2001)
    for eta in (0.1, 0.0):
        d = 100 * (eta_line(F, "M", eta, xq) / 5.5 - 1)
        print(f"eta {eta}: x {x0}-{x1} m  mean {d.mean():+.4f} %  min {d.min():+.4f} (x {xq[d.argmin()]*S:.2f})  max {d.max():+.4f} (x {xq[d.argmax()]*S:.2f})  halfamp {(d.max()-d.min())/2:.4f}")
xs = np.linspace(3.0 / S, 7.95 / S, 34)
d = 100 * (eta_line(F, "M", 0.1, xs) / 5.5 - 1)
print("profile eta0.1:", " ".join(f"{x*S:.2f}:{v:+.3f}" for x, v in zip(xs, d)))
