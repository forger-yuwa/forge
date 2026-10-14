"""作用素の切り分け: 近似 Jacobian の作用と実残差の方向微分の比較 (plan time_integration-line-viscous-jacobian §6.13、2026-10-09)。

forge のライン行列の書き出し (FORGE_LINE_DUMP_DIR: D・Kprev・Knext・各 sweep の dq・dt_local と体積・拘束のフラグ) から補正の方向 p を作る。
p = 最後の sweep の dq ÷ implicitRelax (= 緩和前のライン解)、書き出したライン上だけ非零。
実残差の方向微分 J_t p ≈ −{R(Q+εp) − R(Q−εp)}/(2ε) と、近似の作用 J_a p = (D − V/Δτ I) p − Kprev p_{k−1} − Knext p_{k+1} をライン上の節点で行ごとに比べる。
forge は LHS ΔQ = R (R = res_*) を解くので LHS ≈ V/Δτ − ∂R/∂Q、比べるのは −∂R/∂Q の作用。拘束の行 (壁の運動量・等温壁のエネルギー・軸の半径運動量) は除く。

摂動の作り方 (§6.13 の事前登録): ρ・ρu・ρE に s·p を足す。ρk・ρω は保存量のまま固定 (ソルバの 1 step でも乱流は流れの更新と別に解く)。
ρY_k は Y_k を固定 (ρY_k += Y_k s p_ρ; ライン行列は組成を凍結した 5×5)。

usage (AWS の case dir):
  python3 jprobe.py fields <state_res.h5> <dump_dir> <tag> --eps E [--sweep 4]   → _jprobe/<tag>_{pe,me,ph,mh,pp}/res_0.h5
  python3 jprobe.py extract <run>                     → _jprobe/npz/<run>.npz (全節点の res_*・wi_*、非有限・非正の検査)
  python3 jprobe.py locate <run>                      → 最後の res で |res_roUy| が最大の節点
  python3 jprobe.py compare <dir_dump> <tag> --eps E --ops 7=<dump>,5=<dump> --runs q0=..,q0b=..,pe=..,me=..,ph=..,mh=..,pp=..
"""
import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "_jprobe"
CONS = ("ro", "roUx", "roUy", "roUz", "roe")
RES = ("res_ro", "res_roUx", "res_roUy", "res_roUz", "res_roe")
ROWS = ("mass", "momx", "momy", "momz", "energy")


def load_dump(d: Path):
    meta, relax = {}, None
    for line in (d / "meta.txt").read_text().splitlines():
        if line.startswith("#"):
            relax = float(line.split("implicitRelax")[1].split()[0]); continue
        n, r, c = line.split(); meta[n] = (int(r), int(c))
    a = {k: np.fromfile(d / f"{k}.f64").reshape(v) for k, v in meta.items()}
    a["_relax"] = relax
    a["_nodes"] = a["node_line"][:, 0].astype(int)
    a["_lines"] = a["node_line"][:, 1].astype(int)
    return a


def direction(dump, sweep):
    return dump[f"dqnew_s{sweep}"] / dump["_relax"]          # 緩和前のライン解 (節点, 5)


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def cmd_fields(state: Path, dumpdir: Path, tag: str, eps: float, sweep: int, no_pp: bool = False):
    dump = load_dump(dumpdir)
    p = direction(dump, sweep)
    nodes = dump["_nodes"]
    with h5py.File(state, "r") as h:
        ro0 = np.asarray(h["VALUE"]["ro"][:])[nodes]
        if not np.allclose(ro0, dump["state_ro_roU_roe_cp_gamma"][:, 0], rtol=0, atol=0):
            raise SystemExit("状態の ρ が書き出しの状態と一致しない — 止める")
        ys = sorted(k for k in h["VALUE"] if re.fullmatch(r"roY\d+", k))
        Y = {k: np.asarray(h["VALUE"][k][:])[nodes] / ro0 for k in ys}
    OUT.mkdir(exist_ok=True)
    cases = {"pe": eps, "me": -eps, "ph": 0.5 * eps, "mh": -0.5 * eps, "pp": dump["_relax"]}
    if no_pp:                                   # ε の再試行・診断では実補正の場 (pp) を作り直さない
        cases.pop("pp")
    rec = {"state": str(state), "state_sha256": sha256(state), "dump": str(dumpdir), "eps": eps, "sweep": sweep,
           "relax": dump["_relax"], "species_fixed_Y": ys, "fields": {}}
    for name, s in cases.items():
        d = OUT / f"{tag}_{name}"
        if d.exists():
            raise SystemExit(f"{d} が既にある — 止める")
        d.mkdir()
        shutil.copy(state, d / "res_0.h5")
        for y in state.parent.glob("resolved_species_*.yaml"):     # restart_field が場の化学種の記録を照合する
            shutil.copy2(y, d / y.name)
        with h5py.File(d / "res_0.h5", "a") as h:
            V = h["VALUE"]
            for j, k in enumerate(CONS):
                v = np.asarray(V[k][:]); v[nodes] = v[nodes] + s * p[:, j]; V[k][...] = v
            for k in ys:
                v = np.asarray(V[k][:]); v[nodes] = v[nodes] + Y[k] * s * p[:, 0]; V[k][...] = v
            sumY = sum(np.asarray(V[k][:])[nodes] for k in ys) if ys else None
            ro1 = np.asarray(V["ro"][:])[nodes]
        rec["fields"][name] = {"coef": s, "max_rel_drho": float(np.max(np.abs(s * p[:, 0]) / ro0)),
                               "max_sumY_minus_rho_rel": float(np.max(np.abs(sumY - ro1) / ro1)) if ys else None,
                               "sha256": sha256(d / "res_0.h5")}
        print(f"[jprobe] {d.name}: 係数 {s:+.4g}、max |δρ/ρ| {rec['fields'][name]['max_rel_drho']:.3e}、"
              f"Σρ_Y と ρ の相対差 {rec['fields'][name]['max_sumY_minus_rho_rel']}")
    (OUT / f"{tag}_fields.json").write_text(json.dumps(rec, indent=1, ensure_ascii=False))


def _last_res(run: Path) -> Path:
    fs = sorted((f for f in run.glob("res_*.h5") if re.fullmatch(r"res_\d+\.h5", f.name)), key=lambda f: int(f.name[4:-3]))
    if not fs:
        raise SystemExit(f"{run} に res が無い")
    return fs[-1]


def cmd_extract(run: Path):
    f = _last_res(run)
    (OUT / "npz").mkdir(parents=True, exist_ok=True)
    with h5py.File(f, "r") as h:
        V = h["VALUE"]
        rec = {k: np.asarray(V[k][:], float) for k in RES + ("wi_eheat", "wi_ework", "volume") if k in V}
        bad = {k: int(np.count_nonzero(~np.isfinite(np.asarray(V[k][:])))) for k in V if np.asarray(V[k][:]).dtype.kind == "f"}
        nonpos = {k: int(np.count_nonzero(np.asarray(V[k][:]) <= 0)) for k in ("ro", "P", "T") if k in V}
    missing = [k for k in RES if k not in rec]
    np.savez(OUT / "npz" / f"{run.name}.npz", **rec)
    info = {"run": run.name, "res": f.name, "nonfinite": {k: v for k, v in bad.items() if v}, "nonpositive": {k: v for k, v in nonpos.items() if v},
            "missing": missing, "fields": sorted(rec)}
    (OUT / "npz" / f"{run.name}.json").write_text(json.dumps(info, indent=1, ensure_ascii=False))
    print(f"[jprobe] {run.name}/{f.name}: {sorted(rec)}; 非有限 {info['nonfinite'] or 'なし'}・非正 {info['nonpositive'] or 'なし'}・欠け {missing or 'なし'}")
    if missing or info["nonfinite"] or info["nonpositive"]:
        raise SystemExit(1)


def cmd_locate(run: Path):
    f = _last_res(run)
    with h5py.File(f, "r") as h:
        r = np.abs(np.asarray(h["VALUE"]["res_roUy"][:]))
    i = int(np.argmax(r))
    print(f"[jprobe] {run.name}/{f.name}: |res_roUy| 最大の節点 {i} ({r[i]:.4e}、二乗平均 {np.sqrt(np.mean(r ** 2)):.4e})")
    (OUT / f"{run.name}_locate.json").write_text(json.dumps({"run": run.name, "res": f.name, "node": i, "value": float(r[i])}))


def apply_op(dump, p, lines):
    """J_a p = (D − V/Δτ I) p − Kprev p_{k−1} − Knext p_{k+1} (書き出したライン上)。"""
    n = len(lines)
    D = dump["D"].reshape(n, 5, 5); Kp = dump["Kprev"].reshape(n, 5, 5); Kn = dump["Knext"].reshape(n, 5, 5)
    vdt = dump["dt_vol"][:, 1] / dump["dt_vol"][:, 0]
    Ja = np.einsum("nij,nj->ni", D, p) - vdt[:, None] * p
    for i in range(n):
        if i > 0 and lines[i - 1] == lines[i]:
            Ja[i] -= Kp[i] @ p[i - 1]
        if i < n - 1 and lines[i + 1] == lines[i]:
            Ja[i] -= Kn[i] @ p[i + 1]
    return Ja


def _rel(a, b):
    nb = np.linalg.norm(b)
    return float(np.linalg.norm(a - b) / nb) if nb > 0 else float("nan")


def cmd_compare(dirdump: Path, tag: str, eps: float, ops: dict, runs: dict):
    meta = json.loads((OUT / f"{tag}_fields.json").read_text())
    if abs(meta["eps"] - eps) > 0 or Path(meta["dump"]).resolve() != dirdump.resolve():
        raise SystemExit("fields の記録と ε・方向の書き出しが一致しない — 止める")
    dd = load_dump(dirdump)
    nodes, lines = dd["_nodes"], dd["_lines"]
    n = len(nodes)
    p = direction(dd, meta["sweep"])
    opd = {}
    for k, v in ops.items():
        o = load_dump(v)
        if not (np.array_equal(o["_nodes"], nodes) and np.array_equal(o["_lines"], lines)):
            raise SystemExit(f"作用素 {k} の書き出しの節点が方向の書き出しと違う — 止める")
        opd[k] = o
    R = {k: np.load(OUT / "npz" / f"{v}.npz") for k, v in runs.items()}
    Rm = {k: np.stack([R[k][c][nodes] for c in RES], 1) for k in R}
    fl = dd["flags_wall_iso_axis"]
    dec = np.zeros((n, 5), bool)
    dec[:, 1:4] |= (fl[:, 0] == 1)[:, None]; dec[:, 4] |= (fl[:, 1] == 1); dec[:, 2] |= (fl[:, 2] == 1)
    Jt = -(Rm["pe"] - Rm["me"]) / (2 * eps)
    Jth = -(Rm["ph"] - Rm["mh"]) / eps
    noise = -(Rm["q0b"] - Rm["q0"]) / (2 * eps)            # 同じ Q の再評価の差を方向微分の尺度に換算
    relax = meta["relax"]
    fin = -(Rm["pp"] - Rm["q0"])                          # 実際に掛けた補正 relax·p の有限振幅の応答
    Ja = {k: apply_op(o, p, lines) for k, o in opd.items()}
    # ライン外の節点の応答 (書き出していない K_offline の作用、比べないが大きさを記録)
    allR = {k: np.stack([R[k][c] for c in RES], 1) for k in ("pe", "me")}
    Jt_all = -(allR["pe"] - allR["me"]) / (2 * eps)
    on = np.zeros(Jt_all.shape[0], bool); on[nodes] = True
    off_ratio = [float(np.linalg.norm(Jt_all[~on, r]) / max(np.linalg.norm(Jt_all[on, r]), 1e-300)) for r in range(5)]
    heat = work = None
    if "wi_eheat" in R["pe"]:
        heat = -(R["pe"]["wi_eheat"][nodes] - R["me"]["wi_eheat"][nodes]) / (2 * eps)
        work = -(R["pe"]["wi_ework"][nodes] - R["me"]["wi_ework"][nodes]) / (2 * eps)
    # 出力の残差が書き出しの rhs (sweep 0、ライン外の補正 0) と同じ R(Q) か (rhs は float で持つので相対 1e-6 程度の差は丸め)
    q0_vs_rhs = [_rel(Rm["q0"][~dec[:, r], r], dd["rhs_s0"][~dec[:, r], r]) for r in range(5)]
    ul = np.unique(lines)
    adjacent = [(int(a), int(b)) for a, b in zip(ul[:-1], ul[1:]) if b - a <= 1]
    out = {"tag": tag, "eps": eps, "direction_dump": str(dirdump), "ops": {k: str(v) for k, v in ops.items()}, "runs": runs,
           "q0_vs_rhs_s0": dict(zip(ROWS, q0_vs_rhs)), "adjacent_line_indices": adjacent,
           "offline_response_norm_over_online": dict(zip(ROWS, off_ratio)), "lines": []}
    print("  q0 の残差と書き出しの rhs_s0 の相対差: " + "・".join(f"{r} {v:.1e}" for r, v in zip(ROWS, q0_vs_rhs))
          + (f"; 隣り合う番号のライン {adjacent} (互いの摂動が混ざりうる)" if adjacent else ""))
    print(f"== {tag}: ε {eps}、方向 {dirdump.name}、作用素 {list(ops)}; ライン外の応答のノルム ÷ ライン上: "
          + "・".join(f"{r} {v:.2f}" for r, v in zip(ROWS, off_ratio)))
    for L in np.unique(lines):
        m = lines == L
        rec = {"line": int(L), "node_first": int(nodes[m][0])}
        print(f" ライン {L} (節点 {rec['node_first']}〜):")
        for r, nm in enumerate(ROWS):
            ok = m & ~dec[:, r]
            nt = np.linalg.norm(Jt[ok, r])
            if not ok.any() or nt == 0:
                continue
            x = {"norm_true": float(nt), "fd_eps_vs_half": _rel(Jth[ok, r], Jt[ok, r]), "rerun_noise_over_true": float(np.linalg.norm(noise[ok, r]) / nt),
                 "finite_vs_linear": _rel(fin[ok, r], relax * Jt[ok, r])}
            for k in Ja:
                a = Ja[k][ok, r]
                x[f"op{k}_secant_rel"] = _rel(relax * a, fin[ok, r])      # 実際に掛けた補正での残差の変化と近似の予測 (relax·J_a p) の差
                x[f"op{k}_rel"] = _rel(a, Jt[ok, r])
                x[f"op{k}_norm_ratio"] = float(np.linalg.norm(a) / nt)
                x[f"op{k}_cos"] = float(a @ Jt[ok, r] / max(np.linalg.norm(a) * nt, 1e-300))
                x[f"op{k}_worst_k"] = int(np.flatnonzero(ok)[np.argmax(np.abs(a - Jt[ok, r]))] - np.flatnonzero(m)[0])
            rec[nm] = x
            print(f"  {nm:6s}: 差分 ε/ε/2 {x['fd_eps_vs_half']:.1e}・再評価のノイズ {x['rerun_noise_over_true']:.1e}・有限振幅 {x['finite_vs_linear']:.3f} | "
                  + " | ".join(f"作用素 {k}: 相対差 {x[f'op{k}_rel']:.3f}・大きさ {x[f'op{k}_norm_ratio']:.3f}・cos {x[f'op{k}_cos']:+.3f}・最悪 k={x[f'op{k}_worst_k']}・実補正の予測差 {x[f'op{k}_secant_rel']:.3f}" for k in Ja))
        if heat is not None:
            okE = m & ~dec[:, 4]
            nE = max(np.linalg.norm(Jt[okE, 4]), 1e-300)
            parts = {"heat": heat[okE], "work": work[okE], "rest": Jt[okE, 4] - heat[okE] - work[okE]}
            rec["energy_true_parts_norm_over_total"] = {k: float(np.linalg.norm(v) / nE) for k, v in parts.items()}
            print("  エネルギーの行の実の内訳 (ノルム ÷ 全体): " + "・".join(f"{k} {v:.3f}" for k, v in rec["energy_true_parts_norm_over_total"].items()))
        out["lines"].append(rec)
    (OUT / f"{tag}_compare.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))


def cmd_fhjudge():
    """§6.15 の判定: ライン 2183 のエネルギーの行 (拘束の行を除く) で A (float の面エンタルピー) と B (double) の方向微分を比べる。"""
    dump = load_dump(HERE / "run_0313_e1bdump_m7_linedump")
    nodes, lines = dump["_nodes"], dump["_lines"]
    eps = 1e-6
    m = (lines == 2183) & (dump["flags_wall_iso_axis"][:, 1] != 1)
    x = nodes[m]
    p = direction(dump, 4)
    L = lambda r: np.load(OUT / "npz" / f"{r}.npz")
    Jt, ok = {}, {}
    out = {"line": 2183, "row": "energy", "eps": eps}
    for side in ("a", "b"):
        c = json.loads((OUT / f"s0p7h{side}_compare.json").read_text())
        e = next(l for l in c["lines"] if l["line"] == 2183)["energy"]
        R = {k: L(f"run_0323_jph_{side}_{k}") for k in ("q0", "q0b", "pe", "me")}
        Jt[side] = -(R["pe"]["res_roe"][x] - R["me"]["res_roe"][x]) / (2 * eps)
        ok[side] = e["fd_eps_vs_half"] <= 0.01 and e["rerun_noise_over_true"] <= 1e-3
        Ja = apply_op(dump, p, lines)[m, 4]
        Ja_ns = Ja - dump["scalar_visc_line"][m, 0] * p[m, 4]          # 探索: 値 3 が足したライン面のスカラーを除いた作用
        out[side] = {"fd_eps_vs_half": e["fd_eps_vs_half"], "rerun_noise_over_true": e["rerun_noise_over_true"], "valid": ok[side],
                     "op7_rel": e["op7_rel"], "op7_norm_ratio": e["op7_norm_ratio"], "op7_cos": e["op7_cos"],
                     "op7_noscalar_rel_exploratory": _rel(Ja_ns, Jt[side]), "norm_Jt": float(np.linalg.norm(Jt[side]))}
    qa, qb = L("run_0323_jph_a_q0"), L("run_0323_jph_b_q0")
    g = L("run_0317_jp_s0_q0")
    out["a_q0_vs_lineG_q0_max_abs"] = {c: float(np.max(np.abs(qa[c] - g[c]))) for c in RES}
    out["b_q0_minus_a_q0_rel"] = {c: _rel(qb[c], qa[c]) for c in RES}
    if not (ok["a"] and ok["b"]):
        verdict = "判別不能 (有効の条件を満たさない側がある)"
    else:
        ch = _rel(Jt["b"], Jt["a"])
        out["Jt_change_B_vs_A"] = ch
        if ch >= 0.1:
            verdict = "面の熱力学の精度に依存することを支持"
            verdict += ("; 近似作用素との相対差が半分以下になったので H-c の説明としても支持" if out["b"]["op7_rel"] <= 0.5 * out["a"]["op7_rel"]
                        else "; 近似作用素との相対差は半分以下にならない (H-c の説明としては支持しない)")
        else:
            verdict = "棄却 (J_t p の変化が 10 % 未満)"
    out["verdict"] = verdict
    (OUT / "fh_judge.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    a = sp.add_parser("fields"); a.add_argument("state"); a.add_argument("dump"); a.add_argument("tag")
    a.add_argument("--eps", type=float, required=True); a.add_argument("--sweep", type=int, default=4); a.add_argument("--no-pp", action="store_true")
    a = sp.add_parser("extract"); a.add_argument("run")
    a = sp.add_parser("locate"); a.add_argument("run")
    a = sp.add_parser("fhjudge")
    a = sp.add_parser("compare"); a.add_argument("dump"); a.add_argument("tag"); a.add_argument("--eps", type=float, required=True)
    a.add_argument("--ops", required=True, help="7=<dump>,5=<dump>")
    a.add_argument("--runs", required=True, help="q0=run,q0b=run,pe=run,me=run,ph=run,mh=run,pp=run")
    g = ap.parse_args()
    kv = lambda s: dict(t.split("=", 1) for t in s.split(","))
    if g.cmd == "fields":
        cmd_fields(HERE / g.state, HERE / g.dump, g.tag, g.eps, g.sweep, g.no_pp)
    elif g.cmd == "extract":
        cmd_extract(HERE / g.run)
    elif g.cmd == "locate":
        cmd_locate(HERE / g.run)
    elif g.cmd == "fhjudge":
        cmd_fhjudge()
    else:
        cmd_compare(HERE / g.dump, g.tag, g.eps, {k: HERE / v for k, v in kv(g.ops).items()}, kv(g.runs))
