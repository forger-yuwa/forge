"""⑤ SERN の閉じた制御体積の運動量収支 (plan §5.1 R2, codex M5 採用 2026-09-13)。

力の帳簿 (`sern_forces` / `runner_sern3d.forces3d`: 壁出力 res_<tag>_<id>_<step>.h5 の積分) を、**同じ run の
メッシュ境界面 (BCONDS) 全体**の運動量流束で検算する。定常なら

    Σ_open ∮ [ρ u (u·n) + (p − p_a) n] dA  +  Σ_wall ∮ [ρ u (u·n) + (p − p_a) n] dA  =  0     (n: 流体から外向き)

で、壁項 = 流体が壁に及ぼす力 (帳簿の F)。閉曲面なので p_a ゲージは自由 (Σ n A = 0)。残差を F_ideal で無次元化して
`closure` に、境界群ごとの寄与を `groups` に残す。節点値 (node) は面の節点平均、cell 値は owner cell 値で面に置く
(cell は owner 値、node は境界ノード値をその双対面に置く = ソルバの境界流束の 1 次評価。厳密な離散流束ではないので、帳簿の**範囲・符号・面積**の誤りを検出する目的。数 % の閉じ残差は離散化差)。
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import h5py
import numpy as np

WALL_KINDS = {"slip", "wall", "wall_isothermal", "wall_heatflux"}
OPEN_KINDS = {"inlet_uniformVelocity", "inlet_Pressure", "inlet_Pressure_dir", "outlet_statPress", "outflow", "inlet_profile"}


def parse_bcond(bcond_yaml) -> dict:
    """bcondConfig.yaml → {physID: (name, kind)}。1 行 1 境界の書式 (`name: {physID: n, kind: k, ...}`) 専用。"""
    out = {}
    for line in Path(bcond_yaml).read_text().splitlines():
        m = re.match(r"\s*(\w+):\s*\{physID:\s*(\d+),\s*kind:\s*(\w+)", line)
        if m:
            out[int(m.group(2))] = (m.group(1), m.group(3))
    return out


def momentum_balance(mesh_h5, res_h5, bcond_yaml, p_a: float, F_ideal: float, vehicle_groups=("vehicle",), z_half_w=None,
                     ramp_group: str = "ramp") -> dict:
    """戻り値: {"groups": {name: {kind, Fx_p, Fy_p, Fz_p, Fx_m, Fy_m, Fz_m, area, mdot}}, "closure": {...}, ...}。
    z_half_w を与えると `ramp_group` の面を z ≤ z_half_w (ノズル幅内) と外 (機体) に分けて別群にする (旧 3D メッシュ用)。"""
    bc = parse_bcond(bcond_yaml)
    with h5py.File(mesh_h5, "r") as m, h5py.File(res_h5, "r") as r:
        coord = m["MESH/COORD"][:].reshape(-1, 3).astype(float)
        sv = m["PLANES/surfVect"][:].reshape(-1, 3).astype(float)
        V = r["VALUE"]
        ro = V["ro"][:].astype(float); n_val = len(ro)
        node_mode = (n_val == len(coord))
        u = np.column_stack([V[k][:].astype(float) for k in ("Ux", "Uy", "Uz")]) if "Ux" in V else \
            np.column_stack([V[k][:].astype(float) for k in ("roUx", "roUy", "roUz")]) / ro[:, None]
        P = V["P"][:].astype(float) if "P" in V else V["Ps"][:].astype(float)
        groups = {}
        for b in m["BCONDS"]:
            pid = int(b)
            name, kind = bc.get(pid, (f"phys{pid}", "unknown"))
            g = m[f"BCONDS/{b}"]
            # BCONDS/<id>/iPlanes = 境界 CV 面 (cell: primal 面, node: 境界ノードの双対面 = 境界ノード数)、
            # iCells = その面の owner (cell: セル, node: ノード)。どちらも VALUE[iCells] が面の状態になる
            # (vizBface* は可視化用の primal 面で、node では面数が 1 つ少ない — 収支には使わない)
            ip = g["iPlanes"][:]; nA = sv[ip]                       # 外向き n·A
            ic = g["iCells"][:]; rof, uf, Pf = ro[ic], u[ic], P[ic]
            cen = m["PLANES/centCoords"][:].reshape(-1, 3)[ip]
            un = np.einsum("ij,ij->i", uf, nA)                     # u·n A
            Fm = (rof * un)[:, None] * uf                          # ρ u (u·n) A
            Fp = (Pf - p_a)[:, None] * nA
            parts = [(name, np.ones(len(ip), bool))]
            if z_half_w is not None and name == ramp_group:
                inside = cen[:, 2] <= z_half_w * (1 + 1e-9)
                parts = [(name, inside), (name + "_outside", ~inside)]
            for nm, sel in parts:
                groups[nm] = {"physID": pid, "kind": kind, "wall": kind in WALL_KINDS, "n_faces": int(sel.sum()),
                              "area": float(np.linalg.norm(nA[sel], axis=1).sum()), "mdot": float((rof * un)[sel].sum()),
                              "Fx_p": float(Fp[sel, 0].sum()), "Fy_p": float(Fp[sel, 1].sum()), "Fz_p": float(Fp[sel, 2].sum()),
                              "Fx_m": float(Fm[sel, 0].sum()), "Fy_m": float(Fm[sel, 1].sum()), "Fz_m": float(Fm[sel, 2].sum())}
    tot = np.zeros(3); wall = np.zeros(3); openf = np.zeros(3); mdot_in = mdot_out = 0.0
    for nm, g in groups.items():
        f = np.array([g["Fx_p"] + g["Fx_m"], g["Fy_p"] + g["Fy_m"], g["Fz_p"] + g["Fz_m"]])
        tot += f
        (wall if g["wall"] else openf).__iadd__(f)
        if not g["wall"]:
            mdot_in += -min(g["mdot"], 0.0); mdot_out += max(g["mdot"], 0.0)
    # 帳簿と同じ規約: 推力 = 壁力の −x。壁力 F_wall = Σ_wall (p−p_a) n A (+ 運動量項は slip/no-slip で ≈0)
    T_wall_p = -sum(g["Fx_p"] for g in groups.values() if g["wall"])
    T_wall_nozzle_p = -sum(g["Fx_p"] for nm, g in groups.items() if g["wall"] and nm not in vehicle_groups and not nm.endswith("_outside"))
    return {"node_mode": bool(node_mode), "p_a": p_a, "F_ideal": F_ideal, "groups": groups,
            "closure": {"sum_all": tot.tolist(), "sum_wall": wall.tolist(), "sum_open": openf.tolist(),
                        "residual_over_F_ideal": (np.abs(tot) / F_ideal).tolist(),
                        "mass_in": mdot_in, "mass_out": mdot_out, "mass_imbalance_frac": (mdot_out - mdot_in) / max(mdot_in, 1e-30)},
            "T_wall_all_p": T_wall_p, "T_wall_nozzle_p": T_wall_nozzle_p,
            "C_T_wall_all": T_wall_p / F_ideal, "C_T_wall_nozzle": T_wall_nozzle_p / F_ideal}


def check_run(run_dir, mesh_name: str = "sern.h5", step=None, vehicle_groups=("vehicle",), z_half_w=None, out_path=None) -> dict:
    """run ディレクトリの prepare_info.json / metrics.json と突き合わせて momentum_balance.json を書く。"""
    run_dir = Path(run_dir)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    if step is not None:
        res = [f for f in res if int("".join(c for c in f.stem if c.isdigit())) == int(step)]
    if not res:
        raise FileNotFoundError(f"{run_dir}: res_*.h5 が無い")
    H = info["H_m"]; F_ideal = info["F_ideal_N_per_m"] * (info.get("half_W_m", 1.0) if info.get("dim") == 3 else 1.0)
    mb = momentum_balance(run_dir / mesh_name, res[-1], run_dir / "bcondConfig.yaml", info["states"]["ext"]["P"], F_ideal,
                          vehicle_groups=vehicle_groups, z_half_w=z_half_w)
    mb["res"] = res[-1].name
    mp = run_dir / "metrics.json"
    if mp.exists():
        met = json.loads(mp.read_text())
        mb["ledger"] = {"C_T": met.get("C_T"), "C_T_wall": met.get("C_T_wall"), "C_T_vehicle": met.get("C_T_vehicle"),
                        "C_T_total_with_vehicle": met.get("C_T_total_with_vehicle")}
        if met.get("C_T_wall") is not None:
            mb["ledger"]["C_T_wall_nozzle_diff"] = mb["C_T_wall_nozzle"] - met["C_T_wall"]
    (Path(out_path) if out_path else run_dir / "momentum_balance.json").write_text(json.dumps(mb, indent=1))
    return mb


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="SERN 制御体積の運動量収支 (R2)")
    ap.add_argument("run_dir"); ap.add_argument("--step", type=int, default=None); ap.add_argument("--z-half-w", type=float, default=None,
                                                                                                    help="旧 3D メッシュ用: ramp を z ≤ W/2 と外に分ける")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    mb = check_run(a.run_dir, step=a.step, z_half_w=a.z_half_w, out_path=a.out)
    c = mb["closure"]
    print(f"{a.run_dir} [{mb['res']}] node_mode={mb['node_mode']}")
    print(f"  closure residual / F_ideal = ({c['residual_over_F_ideal'][0]:.4f}, {c['residual_over_F_ideal'][1]:.4f}, {c['residual_over_F_ideal'][2]:.4f})"
          f"  mass imbalance {c['mass_imbalance_frac']*100:+.3f} %")
    for nm, g in mb["groups"].items():
        print(f"  {nm:16s} {g['kind']:22s} Fx_p {g['Fx_p']:+.4e} Fx_m {g['Fx_m']:+.4e} | Fy_p {g['Fy_p']:+.4e} Fy_m {g['Fy_m']:+.4e} | mdot {g['mdot']:+.4e}")
    print(f"  C_T_wall (all walls) {mb['C_T_wall_all']:.5f} / nozzle-only {mb['C_T_wall_nozzle']:.5f}" + (f" / ledger C_T_wall {mb['ledger']['C_T_wall']}" if "ledger" in mb else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
