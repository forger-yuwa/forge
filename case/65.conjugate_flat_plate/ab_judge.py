#!/usr/bin/env python3
r"""case/65 B1 本段の判定 (plan `boundary-cht-conjugate-flat-plate.md` §4.1・§4.1.2、事前登録どおり)。

    python3 ab_judge.py <run_A (連成継続)> <run_B (壁温固定)> <元 run (run_0007)>

1. 残差: `residual_history.csv` の outer_end 行の `rms_ro/rms_roUy/rms_roe` を 20000 step × 3 区間で算術平均。
   横ばい = 3 区間平均の max/min ≤ 1.1。元 run は最後の 60000 step (540000–600000) を同じ 3 区間で。
2. 振幅: `ab_series.npz` (毎 step・FP64) の量ごと・節点群ごとに、区間内の規格化系列の標準偏差の節点最大。
   尺度 P q∞ 709.275 Pa、T 10 K、Uy U∞、界面 q は元 run の評価窓平均 |q|。
   「減衰」= B の最終区間の振幅 ≤ A の最終区間の 1/10、かつ B の 3 区間の振幅が厳密減少。
   測定床: A・B とも最終区間の振幅が 1e-12 未満の系列は比で評価しない。
3. A の再現: A の区間平均残差が元 run の区間平均の 1/2〜2 倍 (3 列 × 3 区間)。かつ、界面残差
   (`conjugate_iface_log_5.csv` の |r_W|/A、`conjugate_history.csv` の res_abs と同じ定義) の最大を取る節点が、
   元 run と A の同じ長さの区間で同じ帯 (前縁 x/L≤0.02 / 後縁 x/L≥0.98 / その間) にあること。
4. 判別:
   - B の最終 2 区間で 3 列がすべて A の 1/10 以下、かつ局所変動も減衰 → 「流体側だけで同程度の停滞を維持する」を棄却
   - B でも 3 列が A の 1/2 以上に残り、B の 3 区間が横ばい → 「動的な連成だけが停滞を維持する」を棄却
   - それ以外 → 中間 (120000 step まで延長)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import plate_common as pc  # noqa: E402
from ab_series import q_scale, refuse  # noqa: E402

COLS = ("rms_ro", "rms_roUy", "rms_roe")
NI, LI = 3, 20000
SCALE = {"P": 709.275, "T": 10.0, "Uy": pc.U_INF}


def resid_means(run, start):
    import csv
    vals = {c: [[] for _ in range(NI)] for c in COLS}
    with open(Path(run) / "residual_history.csv") as f:
        for row in csv.DictReader(f):
            if row["phase"] != "outer_end":
                continue
            s = int(row["step"]) - start
            if 0 <= s < NI * LI:
                for c in COLS:
                    v = float(row[c])
                    if not np.isfinite(v):
                        refuse(f"{run}: {c} に非有限値 (step {row['step']})")
                    vals[c][s // LI].append(v)
    out = {}
    for c in COLS:
        for k in range(NI):
            if len(vals[c][k]) != LI:
                refuse(f"{run}: {c} 区間 {k} の行数 {len(vals[c][k])} (期待 {LI})")
        out[c] = np.array([np.mean(v) for v in vals[c]])
    return out


def iface_loc(run, start):
    """区間ごとに、界面残差 |r_W|/A の最大を取る節点の x/L (更新ごとの最大の中の最大)。"""
    nodes = np.genfromtxt(Path(run) / "conjugate_iface_nodes_5.csv", delimiter=",", names=True)
    A = nodes["A"]; x = nodes["x"] / pc.L
    d = np.genfromtxt(Path(run) / "conjugate_iface_log_5.csv", delimiter=",", names=True)
    st = d["step"].astype(int) - start; i = d["i"].astype(int); r = np.abs(d["r_W"]) / A[i]
    res = []
    for k in range(NI):
        m = (st >= k * LI) & (st < (k + 1) * LI)
        if not m.any():
            refuse(f"{run}: 界面ログに区間 {k} の行が無い")
        j = np.argmax(np.where(m, r, -1))
        res.append((x[i[j]], r[j]))
    return res


def band(xl):
    return "前縁" if xl <= 0.02 + 1e-9 else ("後縁" if xl >= 0.98 - 1e-9 else "中間")


def amps(run, qs):
    d = np.load(Path(run) / "ab_series.npz", allow_pickle=True)
    step = d["step"]
    if len(step) != NI * LI or (np.diff(step) != 1).any():
        refuse(f"{run}: 系列の step が {len(step)} 点・不連続")
    grp = d["group"]; wx = d["wall_x"] / pc.L
    out = {}
    for qn in ("P", "T", "Uy"):
        a = d[qn] / SCALE[qn]
        for g in np.unique(grp):
            out[(g, qn)] = np.array([np.std(a[k * LI:(k + 1) * LI][:, grp == g], axis=0).max() for k in range(NI)])
    q = d["q"] / qs
    for g, m in (("iface_le", wx <= 0.02 + 1e-9), ("iface_te", wx >= 0.98 - 1e-9), ("iface_mid", (wx > 0.02 + 1e-9) & (wx < 0.98 - 1e-9))):
        out[(g, "q")] = np.array([np.std(q[k * LI:(k + 1) * LI][:, m], axis=0).max() for k in range(NI)])
    return out


def main():
    rA, rB, r0 = (Path(p) for p in sys.argv[1:4])
    qs = q_scale()
    print(f"=== B1 判定  A={rA.name}  B={rB.name}  元={r0.name}  (20000 step × 3 区間)")
    RA, RB, R0 = resid_means(rA, 0), resid_means(rB, 0), resid_means(r0, 540000)
    print("  残差の区間平均:")
    for c in COLS:
        print(f"    {c:9s} 元 {' '.join(f'{v:.3e}' for v in R0[c])} | A {' '.join(f'{v:.3e}' for v in RA[c])} | B {' '.join(f'{v:.3e}' for v in RB[c])}")
    # A の再現
    rep = all(0.5 <= RA[c][k] / R0[c][k] <= 2.0 for c in COLS for k in range(NI))
    LA, L0 = iface_loc(rA, 0), iface_loc(r0, 540000)
    same = all(band(a[0]) == band(b[0]) for a, b in zip(LA, L0))
    print(f"  A の再現: 残差比 A/元 ∈ [0.5,2] {'OK' if rep else 'NG'} ("
          + ", ".join(f"{c} {min(RA[c]/R0[c]):.2f}–{max(RA[c]/R0[c]):.2f}" for c in COLS) + ")")
    print(f"            界面残差の最大の位置 {'OK' if same else 'NG'}: 元 " + ", ".join(f"x/L {x:.4f} ({r:.0f} W/m²)" for x, r in L0)
          + " | A " + ", ".join(f"x/L {x:.4f} ({r:.0f} W/m²)" for x, r in LA))
    # 振幅
    AA, AB = amps(rA, qs), amps(rB, qs)
    print("  局所変動の振幅 (規格化標準偏差の節点最大、区間 1/2/3):")
    decay_all = True
    for key in AA:
        a, b = AA[key], AB[key]
        floor = a[-1] < 1e-12 and b[-1] < 1e-12
        dec = floor or (b[-1] <= a[-1] / 10 and b[0] > b[1] > b[2])
        decay_all &= dec
        print(f"    {key[0]:9s} {key[1]:2s} A {' '.join(f'{v:.3e}' for v in a)} | B {' '.join(f'{v:.3e}' for v in b)}"
              f"  B/A {b[-1]/a[-1] if a[-1] > 0 else float('nan'):.3f}{'  測定床' if floor else ('  減衰' if dec else '')}")
    flatB = all(RB[c].max() / RB[c].min() <= 1.1 for c in COLS)
    flatA = all(RA[c].max() / RA[c].min() <= 1.1 for c in COLS)
    low = all((RB[c][1:] <= RA[c][1:] / 10).all() for c in COLS)
    high = all((RB[c] >= RA[c] / 2).all() for c in COLS)
    print(f"  横ばい (max/min ≤ 1.1): A {'yes' if flatA else 'no'} ("
          + ", ".join(f"{c} {RA[c].max()/RA[c].min():.3f}" for c in COLS) + "), B " + ('yes' if flatB else 'no') + " ("
          + ", ".join(f"{c} {RB[c].max()/RB[c].min():.3f}" for c in COLS) + ")")
    print("  B/A 残差比 (区間 1/2/3): " + ", ".join(f"{c} " + "/".join(f"{v:.3f}" for v in RB[c] / RA[c]) for c in COLS))
    if not (rep and same):
        print("VERDICT: 判定不能 (A が元の停滞を再現していない)")
        return 3
    if low and decay_all:
        print("VERDICT: 「流体側だけで同程度の停滞を維持する」を棄却 (連成が停滞を維持)")
    elif high and flatB:
        print("VERDICT: 「動的な連成だけが停滞を維持する」を棄却 (流体側の要因。slip とは断定しない)")
    else:
        print("VERDICT: 中間 (120000 step まで延長)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
