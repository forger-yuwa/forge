#!/usr/bin/env python3
"""`slauWallNormalChi` の適用診断が**できる設定か**を検査し、できるなら面ごとの表を出す
(plan convection-slau-wall-normal-chi-usage-rule §4.1 / §4.1.1、§5.1 #3b)。

`diag_wall_cv_budget.py` は 1 次 SLAU の再計算で、実カーネルの `slauContactFloor` の質量流束追加
(`convectiveFlux_slau_d.inc.cuh:620`) や `sstEnergyIncludesK` の p* 差 (:516) を持たない。設定によっては
流束の符号まで食い違う (codex plan-3 M3 の再計算: +0.099 対 −0.191)。そこで

  1. **キー検査**: 診断区間の実効値を `solverConfig.yaml` と `forge_run.log` の起動エコーの両方から取り
     (省略 ≡ 既定値に正規化、両者が食い違えば診断不能)、§4.1.1 の表を要求する。
  2. **面流束の照合**: 同じ run を 1 step、`FORGE_DUMP_MASSFLUX=<path>` で回した `massflux` (カーネル値) と
     `<path>.state` (カーネルが読んだ状態: ro,Ux,Uy,Uz,P,sonic) から、対象 CV の**全接続面**で
     ツール値 (`slau_mdot`) との差が前 plan V5 の許容 2 τ_b (τ_b = 1e-5 A ρ̄ ĉ) 以内か。
  3. 1・2 とも通ったときだけ「診断可能」とし、面ごとに chi, chi_n, ΔP, mdot(chi), mdot(chi_n) の表と、
     CV への正味補充の変化 Σ[流出(chi_n) − 流出(chi)] (< 0 で補充が増える、§4.1 条件 iii) を出す。
     ρ_i = 対象壁 CV の**壁でない隣接ノード全部の平均**。

**判定は「診断可能 / 診断不能 (理由)」**。診断不能は「flag 1 が不要」の意味ではない (既知構成以外で 1 にしない、の根拠になる)。

usage:
  python3 diag_applicability.py RUN_DIR --mesh RUN_DIR/sern.h5 --massflux mf.bin --ids 153797,153880 \
      --wall-phys-ids 1,2,3,4,10,11,12,13,15 [--config-only]
"""
import argparse, re, sys
from pathlib import Path
import h5py, numpy as np, yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from diag_wall_cv_budget import slau_mdot, parse_struct

# (表示名, yaml パス, 要求 (集合 or 関数), 省略時の既定値, 理由)  — plan §4.1.1 の表
REQ = [
    ("mesh.discretization", ("mesh", "discretization"), {"node"}, None, "対象構成"),
    ("mesh.nodeWallDirichlet", ("mesh", "nodeWallDirichlet"), {"1"}, "0", "対象構成"),
    ("solver", ("solver",), {"SLAU", "SLAU2"}, None, "質量流束式が同じ (SLAU2 は圧力束のみ差)"),
    ("space.convMethod", ("space", "convMethod"), {"0"}, None, "ツールは再構成なし"),
    ("space.slauContactFloor", ("space", "slauContactFloor"), {"0", "0.0"}, "0", ":620 の追加項がツールに無い"),
    ("turbulence.sstEnergyIncludesK", ("turbulence", "sstEnergyIncludesK"), {"0"}, "0", ":516 の p* 差"),
    ("space.lowMachPrecond", ("space", "lowMachPrecond"), {"0"}, "0", "散逸スケール c'"),
    ("time.deltaT.lowMachPrecond", ("time", "deltaT", "lowMachPrecond"), {"0"}, "0", "散逸スケール c' (deltaT 側の綴り)"),
    ("space.slauWallNormalChi", ("space", "slauWallNormalChi"), {"0"}, "0", "--wall-normal-chi は置換予測用"),
    ("space.badReconFallback", ("space", "badReconFallback"), {"0"}, "0", "面状態の差し替え"),
    ("mesh.isAxisymmetric", ("mesh", "isAxisymmetric"), {"0"}, "0", "適用範囲外"),
]
ABSENT = [("condensation", "凝縮・二相 (面エンタルピー・P=ρ(1−g)RT)")]


def yget(doc, path):
    cur = doc
    for k in path:
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
    return None if isinstance(cur, (dict, list)) else str(cur).strip().strip('"\'')


def log_echo(log_text):
    """`'key' in 'section': value` / `'key': value` の起動エコーを {section.key: value} で返す。"""
    out = {}
    for key, sec, val in re.findall(r"^'(\w+)' in '([\w.]+)':\s*(\S+)", log_text, re.M):
        out[f"{sec}.{key}"] = val
    for k, v in re.findall(r"^'(\w+)':\s*(\S+)", log_text, re.M):
        out.setdefault(k, v)
    return out


def check_config(run_dir):
    cfg = yaml.safe_load((Path(run_dir) / "solverConfig.yaml").read_text()) or {}
    logp = Path(run_dir) / "forge_run.log"
    echo = log_echo(logp.read_text(errors="replace")) if logp.exists() else {}
    rows, bad = [], []
    for name, path, want, dflt, why in REQ:
        yv = yget(cfg, path)
        v = yv if yv is not None else dflt
        ev = echo.get(name) if name in echo else echo.get(path[-1]) if len(path) == 1 else None
        vn = v.upper() if name == "solver" and v else v
        evn = ev.upper() if name == "solver" and ev else ev
        ok = vn in want if vn is not None else False
        if ev is not None and vn is not None and evn != vn:
            try:
                same = float(evn) == float(vn)
            except (TypeError, ValueError):
                same = False
            if not same:
                ok = False
                why = f"yaml と起動エコーが食い違う ({vn} / {evn})"
        rows.append((name, vn, evn, sorted(want), ok, why))
        if not ok:
            bad.append(f"{name}={vn} (要求 {sorted(want)}; {why})")
    for key, why in ABSENT:
        if key in cfg and cfg[key]:
            bad.append(f"{key} が有効 ({why})")
            rows.append((key, "有効", None, ["なし"], False, why))
    tm = yget(cfg, ("physProp", "thermalMethod"))
    rows.append(("physProp.thermalMethod", tm, echo.get("physProp.thermalMethod"), ["記録のみ (CPG/TP をツールの EOS 経路と合わせる)"], True, "ĉ・T"))
    return rows, bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--mesh", default=None)
    ap.add_argument("--massflux", default=None, help="FORGE_DUMP_MASSFLUX の出力 (同名 .state も読む)")
    ap.add_argument("--ids", default="", help="対象壁 CV (節点) をカンマ区切り")
    ap.add_argument("--wall-phys-ids", default="", help="壁 physID (bcondConfig の wall 系) をカンマ区切り")
    ap.add_argument("--config-only", action="store_true", help="キー検査だけ行う")
    a = ap.parse_args()

    rows, bad = check_config(a.run_dir)
    print("=== キー検査 (plan §4.1.1) ===")
    for name, v, ev, want, ok, why in rows:
        print(f"  {'OK ' if ok else 'NG '} {name:32s} yaml/既定={v!s:8s} 起動エコー={ev!s:8s} 要求={want}  ({why})")
    if bad:
        print("VERDICT: 診断不能 (キー検査)\n  " + "\n  ".join(bad))
        return 2
    if a.config_only:
        print("VERDICT: キー検査のみ OK (面流束の照合は未実施 — これだけでは「診断可能」と言わない)")
        return 0
    if not (a.mesh and a.massflux and a.ids and a.wall_phys_ids):
        sys.exit("--mesh --massflux --ids --wall-phys-ids が要る (または --config-only)")

    ids = [int(s) for s in a.ids.split(",") if s]
    with h5py.File(a.mesh) as f:
        S = np.asarray(f["PLANES/surfVect"], np.float64).reshape(-1, 3)
        own, nei, _ = parse_struct(np.asarray(f["PLANES/STRUCT"]), len(S))
        wall = set()
        for pid in [s for s in a.wall_phys_ids.split(",") if s]:
            if pid in f["BCONDS"]:
                c = np.asarray(f["BCONDS"][pid]["iCells"]).ravel(); wall.update(c[c >= 0].tolist())
    own = own.astype(int); nei = nei.astype(int)
    mf = np.fromfile(a.massflux, dtype=np.float32).astype(np.float64)
    if mf.size != len(S):
        sys.exit(f"massflux {mf.size} 面 != メッシュ {len(S)} 面")
    st = np.fromfile(a.massflux + ".state", dtype=np.float32).astype(np.float64)
    nC = st.size // 6
    ro, Ux, Uy, Uz, P, c = st.reshape(6, nC)
    state = lambda i: {"ro": ro[i], "Ux": Ux[i], "Uy": Uy[i], "Uz": Uz[i], "P": P[i], "sonic": c[i]}

    worst, nf, fails = 0.0, 0, []
    all_ok = True
    for cv in ids:
        faces = np.where(((own == cv) | (nei == cv)) & (nei >= 0))[0]
        nbr = [int(nei[k]) if own[k] == cv else int(own[k]) for k in faces]
        nonwall = [j for j in nbr if j not in wall]
        rho_i = float(np.mean(ro[nonwall])) if nonwall else float("nan")
        print(f"\n=== CV {cv}  ρ_w {ro[cv]:.5e}  ρ_i (非壁隣接 {len(nonwall)} 点の平均) {rho_i:.5e}  ρ_w/ρ_i {ro[cv]/rho_i:.4f} ===")
        print(f"  {'face':>9} {'nbr':>9} {'chi':>7} {'chi_n':>7} {'dP':>10} {'mdot_k':>12} {'mdot_tool':>12} {'|d|/2τb':>8} {'mdot(chi_n)':>12}")
        d_out = 0.0
        for k, j in zip(faces, nbr):
            A = float(np.linalg.norm(S[k])); n = S[k] / A
            L, R = state(own[k]), state(nei[k])
            m0, _, _, chi, ch, _, _ = slau_mdot(A, *n, L, R, wall_face=False)
            is_wf = (own[k] in wall) or (nei[k] in wall)
            m1, _, _, chin, _, _, _ = slau_mdot(A, *n, L, R, wall_face=is_wf)
            tau = 1e-5 * A * 0.5 * (L["ro"] + R["ro"]) * ch
            r = abs(mf[k] - m0) / (2 * tau)
            worst = max(worst, r); nf += 1
            if r > 1.0:
                all_ok = False; fails.append((cv, int(k), r))
            sgn = 1.0 if own[k] == cv else -1.0            # CV から見た流出を正に
            d_out += sgn * (m1 - m0)
            print(f"  {k:9d} {j:9d} {chi:7.4f} {chin:7.4f} {R['P']-L['P']:10.3e} {mf[k]:12.4e} {m0:12.4e} {r:8.3f} {m1:12.4e}")
        print(f"  Σ[流出(chi_n) − 流出(chi)] = {d_out:+.4e} kg/s  ({'補充が増える' if d_out < 0 else '補充は増えない'})")
    print(f"\n=== 面流束の照合: {nf} 面、max |mdot_kernel − mdot_tool|/(2τ_b) = {worst:.3f} ===")
    if not all_ok:
        print("VERDICT: 診断不能 (ツールとカーネルの面流束が一致しない)\n  " +
              "\n  ".join(f"CV {cv} face {k}: {r:.2f}×2τ_b" for cv, k, r in fails[:20]))
        return 3
    print("VERDICT: 診断可能 (キー検査 OK・面流束がツールと 2τ_b 以内で一致)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
