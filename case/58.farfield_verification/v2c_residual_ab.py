#!/usr/bin/env python3
"""V2c 角の履歴差の切り分け (plan boundary-node-farfield-characteristic §5.1 #3a、codex diagnose 2 回目 2026-09-29)。
同じ低領域 B の格子・BC・設定に、2 つの最終場 (一様 IC から = run_0064_v2c_B/res_6000、C から = run_0094_v2c_B_fromC/res_6000) を
index コピーで与え、1 step の初回評価の**絶対残差** (帳簿 res_final = 壁射影等の後、更新に渡る値) を全節点・全 5 保存量で見る。
診断閾値 (事前固定): |R_q| / Σ面|F_q| ≤ 1e-5。分母 = その節点に接する内部面流束の絶対和 (帳簿 .faces) + 境界半割面の規模
(質量 |ρ u·S|、運動量 |ρ u (u·S)| + |P − pRef| |S|、エネルギー |ρ H u·S|; 節点の状態から)。規模 0 は自由流基準 × 節点の境界面積。
これは収束基準の代用ではない。
  python3 v2c_residual_ab.py setup | compare
一般化 (2026-10-03、#3c の次数 A/B): 状態を環境変数で渡せる。V2C_STATES="DST=SNAP@CFGRUN;..." (CFGRUN の solverConfig/bcond を使う、省略時は SRC)、
V2C_BIN=<forge 実行ファイル>、V2C_SRC=<格子を読む run>。
"""
import csv, os, shutil, subprocess, sys
from collections import defaultdict
import h5py, numpy as np

SRC = os.environ.get("V2C_SRC", "run_0064_v2c_B")
STATES = {"run_0097_v2c_resab_fromU": "run_0064_v2c_B/res_6000.h5", "run_0098_v2c_resab_fromC": "run_0094_v2c_B_fromC/res_6000.h5"}
CFG = {}
if os.environ.get("V2C_STATES"):
    STATES = {}
    for item in os.environ["V2C_STATES"].split(";"):
        dst, rest = item.split("=", 1); snap, _, cfg = rest.partition("@")
        STATES[dst] = snap
        if cfg:
            CFG[dst] = cfg
BIN = os.environ.get("V2C_BIN", os.path.expanduser("~/forge-pgrad-new/solver_density_cuda/build-ff/forge"))
CONS = ("ro", "roUx", "roUy", "roUz", "roe")
P_REF, GAM, CP = 101325.0, 1.4, 1004.5
RI = P_REF / ((CP - CP / GAM) * 300.0); UI = 2.5 * (GAM * (CP - CP / GAM) * 300.0) ** 0.5
HI = CP / GAM * 300.0 + P_REF / RI + 0.5 * UI * UI   # e = cv T (CPG、setup_v2c と同じ) + P/ρ + U²/2


def setup():
    with h5py.File(os.path.join(SRC, "ramp.h5")) as f:
        n = f["VALUE/ro"].shape[0]
    for dst, snap in STATES.items():
        os.makedirs(dst)
        for f in ("solverConfig.yaml", "bcondConfig.yaml", "probe.yaml"):
            shutil.copy(os.path.join(CFG.get(dst, SRC), f), dst)
        shutil.copy(os.path.join(SRC, "ramp.h5"), dst)
        import re as _re
        s = _re.sub(r"nStepOuter: \d+", "nStepOuter: 1", open(os.path.join(dst, "solverConfig.yaml")).read())
        s = _re.sub(r"outStepInterval: \d+", "outStepInterval: 1", s)
        open(os.path.join(dst, "solverConfig.yaml"), "w").write(s)
        with h5py.File(snap) as s5, h5py.File(os.path.join(dst, "ramp.h5"), "r+") as f:
            assert s5["VALUE/ro"].shape[0] == n, "同一格子でない"
            for q in CONS:
                f["VALUE"][q][...] = s5["VALUE"][q][:]
        open(os.path.join(dst, "IC_FROM.txt"), "w").write(f"{snap} の保存量を index コピー (同一格子 {n} 節点)\n")
        env = dict(os.environ, FORGE_CUDA_BLOCKSIZE="128", FORGE_DUMP_LEDGER="ledger.csv", FORGE_DUMP_LEDGER_CALLS="1",
                   FORGE_DUMP_LEDGER_NODES="all", FORGE_DUMP_FARFIELD="ffdump")
        r = subprocess.run([BIN], cwd=dst, env=env, stdin=subprocess.DEVNULL, stdout=open(os.path.join(dst, "forge_run.log"), "w"), stderr=subprocess.STDOUT)
        print(dst, "exit", r.returncode)


def compare():
    with h5py.File(os.path.join(SRC, "ramp.h5")) as f:
        xyz = np.array(f["MESH/COORD"]).reshape(-1, 3)
        sv = np.array(f["PLANES/surfVect"]).reshape(-1, 3)
        bpl, slipn = {}, {}
        SLIP = {3, 5, 6} | ({4} if "kind: slip" in open(os.path.join(SRC, "bcondConfig.yaml")).read().split("ymax")[1].split("\n")[0] else set())
        for b in f["BCONDS"].keys():
            ic, ip = f[f"BCONDS/{b}/iCells"][:], f[f"BCONDS/{b}/iPlanes"][:]
            for c, p in zip(ic, ip):
                bpl.setdefault(int(c), []).append(sv[int(p)])
                if int(b) in SLIP:
                    slipn.setdefault(int(c), {}).setdefault(int(b), np.zeros(3))
                    slipn[int(c)][int(b)] += sv[int(p)]
    n = len(xyz)
    ok_all = True
    for run in STATES:
        st, res = defaultdict(dict), defaultdict(dict)
        with open(os.path.join(run, "ledger.csv")) as fh:
            rd = csv.reader(fh); next(rd)
            for row in rd:
                if row[0] != "1":
                    continue
                if row[1] == "after_eos_bc":
                    st[row[3]][int(row[2])] = float(row[4])
                elif row[1] == "res_final":
                    res[row[3]][int(row[2])] = float(row[4])
        sc = defaultdict(lambda: np.zeros(n))
        area = np.zeros(n)
        seen = set()
        with open(os.path.join(run, "ledger.csv.faces")) as fh:
            for r in csv.DictReader(fh):
                if r["call"] != "1" or r["ip"] in seen:
                    continue
                seen.add(r["ip"])
                # 帳簿の F_roU* は圧力項を含まない (接線流の y・z 運動量で分母が 0 に近くなる) ので、面の圧力 × 面積成分を足す
                pf = abs(0.5 * (float(r["P_L"]) + float(r["P_R"])) - P_REF)
                sd = {"roUx": abs(float(r["sx"])), "roUy": abs(float(r["sy"])), "roUz": abs(float(r["sz"]))}
                for c, q in (("F_ro", "ro"), ("F_roUx", "roUx"), ("F_roUy", "roUy"), ("F_roUz", "roUz"), ("F_roe", "roe")):
                    v = abs(float(r[c])) + (pf * sd[q] if q in sd else 0.0)
                    sc[q][int(r["ic0"])] += v; sc[q][int(r["ic1"])] += v
                a = float(r["ss"]); area[int(r["ic0"])] += a; area[int(r["ic1"])] += a
        for c, planes in bpl.items():
            ro = st["ro"][c]; u = np.array([st["Ux"][c], st["Uy"][c], st["Uz"][c]]); P = st["P"][c]
            H = (st["roe"][c] + P) / ro
            for S in planes:
                un = float(u @ S); A = float(np.linalg.norm(S))
                sc["ro"][c] += abs(ro * un); sc["roe"][c] += abs(ro * H * un)
                for d, q in enumerate(("roUx", "roUy", "roUz")):
                    sc[q][c] += abs(ro * u[d] * un) + abs(P - P_REF) * abs(S[d])
        print(f"{run}: 内部面 {len(seen)}、境界半割面を持つ節点 {len(bpl)} (うち slip {len(slipn)}、運動量は{'接線成分 (--proj)' if PROJ else '射影なしの全成分'}で判定)")
        # slip 面を持つ節点は運動量残差から壁法線成分 (面ごとの法線、Gram–Schmidt) を除く: 更新時に射影で捨てられる成分なので収束の指標でない
        Rm = np.stack([np.array([res["res_" + q].get(i, np.nan) for i in range(n)]) for q in ("roUx", "roUy", "roUz")], 1)
        for c, dct in (slipn.items() if PROJ else []):
            basis = []
            for v in dct.values():
                w = v / np.linalg.norm(v)
                for b in basis:
                    w = w - (w @ b) * b
                if np.linalg.norm(w) > 1e-6:
                    basis.append(w / np.linalg.norm(w))
            for b in basis:
                Rm[c] = Rm[c] - (Rm[c] @ b) * b
        Rproj = {"roUx": Rm[:, 0], "roUy": Rm[:, 1], "roUz": Rm[:, 2]}
        ok = True
        for q in CONS:
            R = Rproj[q] if q in Rproj else np.array([res["res_" + q].get(i, np.nan) for i in range(n)])
            # 規模が数値的に 0 の成分 (一様流域の y・z 運動量など) は自由流基準 × 節点の面積 (単位を合わせる)
            ref = {"ro": RI * UI, "roe": RI * UI * HI}.get(q, RI * UI * UI) * area
            # 運動量 3 成分は共通の規模 (面流束 3 成分の絶対和、V1d と同じ): 成分別では流れに垂直な成分の分母が痩せる
            s = (sc["roUx"] + sc["roUy"] + sc["roUz"]).copy() if q in ("roUx", "roUy", "roUz") else sc[q].copy()
            zero = s < 1e-6 * ref; s[zero] = ref[zero]
            r = np.abs(R) / s
            w = int(np.nanargmax(r))
            corner = (xyz[:, 0] >= 0.195) & (xyz[:, 0] <= 0.23) & (xyz[:, 1] <= 0.01)
            exitb = (xyz[:, 0] >= 1.195) & (xyz[:, 1] <= 0.09)
            convex = (xyz[:, 0] >= 0.69) & (xyz[:, 0] <= 0.76) & (xyz[:, 1] <= 0.1)
            nbad = int(np.sum(r > 1e-5))
            ok &= nbad == 0 and np.all(np.isfinite(R))
            top = np.argsort(-np.nan_to_num(r))[:3]
            w = int(np.nanargmax(r))
            print(f"  {q:5s}: max |R|/Σ|F| {r[w]:.3e} (x {xyz[w,0]:.3f} y {xyz[w,1]:.4f})、> 1e-5 の節点 {nbad}、圧縮角 {np.nanmax(r[corner]):.2e}、凸角近傍 {np.nanmax(r[convex]):.2e}、出口下端 {np.nanmax(r[exitb]):.2e}、非有限 {int(np.sum(~np.isfinite(R)))}、上位 {[(round(float(xyz[t,0]),3), round(float(xyz[t,1]),4), round(float(xyz[t,2]),3), f'{r[t]:.1e}') for t in top]}")
        print(f"  → {'全項目 ≤ 1e-5' if ok else '閾値超えあり'}")
        ok_all &= ok
    print("VERDICT: " + ("両状態とも全節点で ≤ 1e-5 → 残差を残した停止の仮説を棄却 (異なる離散平衡の候補として扱う)" if ok_all
                         else "少なくとも片方が閾値超え → 「二つとも離散定常解」を棄却し、非零残差が更新されない箇所を追う"))


PROJ = "--proj" in sys.argv   # 既定は射影しない (codex diagnose 3 回目: slip は弱形式でソルバは法線残差を捨てない。旧版は射影していた)

if __name__ == "__main__":
    {"setup": setup, "compare": compare}[sys.argv[1]]()
