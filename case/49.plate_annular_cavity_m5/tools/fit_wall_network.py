#!/usr/bin/env python3
"""壁間熱回路 (gas-mediated conductance network) の同定 (plan §4.7.6 / §5.1 #22)。

§4.7.6 で「等温壁 3 本の $h_{aw}$ 包絡は上限にならない」ことが確定した。原因は熱回路

    Q_i = G_i0 (T_aw - T_i) + Σ_{j≠i} G_ij (T_j - T_i)          [W]

の第 2 項 (他の壁からガス越しに受ける熱) で、**全壁同温の run では恒等的に 0** になり
同定できない。非一様壁温 run を足すとこの項が測れる。

  G_i0 [W/K] … 外部流からその壁への実効コンダクタンス (= h_aw,i * A_i)
  G_ij [W/K] … 壁 i と壁 j がガスを介して交わすコンダクタンス (対称, G_ij = G_ji)

3 壁 (外筒内壁 / 円柱側面 / 底面) なら未知数は G_i0 が 3 個 + G_ij が 3 個 = 6 個。
一様 run 1 本で 3 式 (G_ij は寄与しない)、非一様 run 1 本で 3 式。**非一様 2 本で
G_ij が過決定**になり、整合性まで確認できる。

**これは定数 G の最小モデル**である。実測では G_i0 自体が壁温で ±17 % 動く (§4.7.2) ので、
残差はその分残る。残差を各 run ごとに出すので、定数近似で足りるかを見て判断すること。

usage: python3 tools/fit_wall_network.py RUN [RUN ...] [--predict RUN]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(CASE))
import cavity_eval as ce  # noqa: E402

WALLS = ["cav_outer", "cyl_side", "cav_floor"]
JP = {"cav_outer": "外筒内壁", "cyl_side": "円柱側面", "cav_floor": "底面"}
PAIRS = [(0, 1), (0, 2), (1, 2)]


def read_run(run):
    run = Path(run)
    ev = json.loads((run / "cavity_eval.json").read_text())
    D = ce.run_conditions(run)
    by = D.get("wall_T_by_group", {}) or {}
    T = np.array([float(by.get(w, D["wall_T"])) for w in WALLS])
    Q = np.array([float(ev["wall"][w]["Q_W"]) for w in WALLS])      # 半割 [W]
    return dict(name=run.name, T=T, Q=Q, Taw=float(ce.taw_of(D)),
                uniform=bool(len(set(np.round(T, 6))) == 1))


def build(rows):
    """未知数 x = [G_00, G_10, G_20, G_01, G_02, G_12] に対する A x = b を組む。"""
    A, b = [], []
    for r in rows:
        for i in range(3):
            row = np.zeros(6)
            row[i] = r["Taw"] - r["T"][i]
            for k, (p, q) in enumerate(PAIRS):
                if i == p:
                    row[3 + k] = r["T"][q] - r["T"][i]
                elif i == q:
                    row[3 + k] = r["T"][p] - r["T"][i]
            A.append(row)
            b.append(r["Q"][i])
    return np.array(A), np.array(b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--predict", default=None, help="同定に使わず、予測誤差だけ見る run")
    a = ap.parse_args()
    rows = [read_run(r) for r in a.runs]
    if not any(not r["uniform"] for r in rows):
        raise SystemExit("非一様壁温の run が 1 本も無い -> G_ij は同定できない "
                         "(gen_runs.py --wall-temps で作る)")

    A, b = build(rows)
    x, *_ = np.linalg.lstsq(A, b, rcond=None)
    G0, Gij = x[:3], x[3:]
    print("同定に使った run (半割 Q [W]):")
    for r in rows:
        print("  %-44s %s  Tw %s" % (r["name"], "一様" if r["uniform"] else "非一様",
                                     " ".join("%.1f" % t for t in r["T"])))
    print("\n外部流からのコンダクタンス G_i0 [W/K] (= h_aw,i x A_i)")
    for i, w in enumerate(WALLS):
        print("   %-10s %10.5f" % (JP[w], G0[i]))
    print("\n壁どうし (ガス越し) のコンダクタンス G_ij [W/K]")
    for k, (p, q) in enumerate(PAIRS):
        print("   %-10s - %-10s %10.5f" % (JP[WALLS[p]], JP[WALLS[q]], Gij[k]))
        if Gij[k] < 0:
            print("      ** 負 = このモデルで説明できていない (定数 G の限界か run 不足) **")

    print("\n各 run の再現性 (実測 Q vs モデル Q, 半割 [W])")
    for r in rows + ([read_run(a.predict)] if a.predict else []):
        Ar, br = build([r])
        pred = Ar @ x
        tag = "  [予測のみ]" if a.predict and r["name"] == Path(a.predict).name else ""
        print("  %s%s" % (r["name"], tag))
        for i, w in enumerate(WALLS):
            e = (pred[i] - br[i]) / max(abs(br[i]), 1e-30) * 100
            print("     %-10s 実測 %9.4g   モデル %9.4g   誤差 %+7.1f %%"
                  % (JP[w], br[i], pred[i], e))

    # 包絡法との比較 (§4.7.6 の主張を数値で残す)
    uni = [r for r in rows if r["uniform"]]
    if uni:
        print("\n参考: 一様 run だけから作る「包絡」予測 (G_ij = 0 と置いたもの)")
        for r in rows:
            if r["uniform"]:
                continue
            for i, w in enumerate(WALLS):
                env = max(u["Q"][i] / max(u["Taw"] - u["T"][i], 1e-30) for u in uni)
                q_env = env * (r["Taw"] - r["T"][i])
                print("  %-30s %-10s 実測 %9.4g   包絡 %9.4g   **実測/包絡 %.2f 倍**"
                      % (r["name"], JP[w], r["Q"][i], q_env,
                         r["Q"][i] / max(q_env, 1e-30)))


if __name__ == "__main__":
    main()
