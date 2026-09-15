#!/usr/bin/env python3
"""化学種の集合/順序が違う run へ場を移す (種変換 restart)。同一メッシュ・index コピー。

旧 `[MIXDRY, H2O]` (擬似種) の収束場を新順序 (例 `[H2O, N2, O2, AR, CO2]` の full) へ移すとき、
擬似種の保存量 ρY_MIXDRY を構成種へ `species_meta.yaml` の `expansion` (lump 内質量分率) で分配し、
H2O (凝縮種) と凝縮モーメント `rog_*/roQ*_*`・トレーサ `roXi` を**名前で**移す。逆 (full → lumped) は
destination の lump に構成種を合算する。DB (`species_db.yaml`) やエンタルピー基準 (`thermoHrefTemp`) が
変わる場合は `roe` を再構成する: `roe += ρ [e_gas,dst(Y_dst, T) − e_gas,src(Y_src, T)]` (両 DB とも NASA-9,
`total_quantities._TPGas` と同式; 差分形なので凝縮セルの液相エネルギー (二相 EOS, g·L ≈ 2.4 MJ/kg) はそのまま保たれる。
source の DB が読めないときだけ乾き気相の完全再構成 ρ(e_gas,dst(T)+u²/2) に落ちる = 凝縮セルでは不正確)。
移した後に ΣρY=ρ (1e-6)・総水量 (ρY_w + ρg の総和)・T の保存を検査して要約を出す。
plans/active/thermophysics-cea-mole-fraction-species.md §2 (forge 本体) / §4.5 M4。

  convert_species_field.py SRC_res.h5 DST_input.h5 --meta DST/species_meta.yaml [--src-meta SRC/species_meta.yaml]
      [--src-run SRC_DIR] [--dst-run DST_DIR] [--reconstruct-roe auto|always|never] [--dry-run]

- SRC: res_*.h5 (原始量 P,T,Ux,.. + Y{s}) か input h5 (保存量 roY{s}; T は SRC DB で反転)。
- DST: convertGmshToForge 直後 (同一メッシュ・同一 CV 数) の input h5。ro/roU/roe/roK/roOmega も index コピーする。
- 種名は `species_meta.yaml` (`species`, `expansion`, `condensing_species`, `MW`) を正とし、無ければ
  run dir の `solverConfig.yaml` (`physProp.species`) から取る (`--src-run/--dst-run` 省略時は h5 の隣)。
- 行き先の無い実種 (destination のどの種にも展開されない) があれば拒否する。
"""
import argparse, os, sys
import numpy as np, h5py, yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from total_quantities import _TPGas, RU  # noqa: E402
from forge_species import species_info  # noqa: E402


def _up(s):
    return str(s).upper()


def load_layout(meta_path, run_dir, label):
    """{names, expansion {transported: {real: w}}, condensing, tracer, MW, db, Tref} を返す。"""
    meta = yaml.safe_load(open(meta_path)) if meta_path else None
    info = species_info(run_dir) if run_dir and os.path.exists(os.path.join(run_dir, "solverConfig.yaml")) else None
    if meta is None and info is None:
        raise SystemExit(f"{label}: species_meta.yaml も solverConfig.yaml も無い (--meta / --src-meta / --src-run / --dst-run)")
    names = [_up(s) for s in (meta["species"] if meta else info["names"])]
    if info and [_up(s) for s in info["names"]] != names:
        raise SystemExit(f"{label}: species_meta.yaml の species {names} と solverConfig の {info['names']} が違う")
    exp = {}
    for s in names:
        row = (meta or {}).get("expansion", {}).get(s) if meta else None
        if row is None:
            row = {s: 1.0}
        exp[s] = {_up(k): float(v) for k, v in row.items()}
        tot = sum(exp[s].values())
        if abs(tot - 1.0) > 1e-9:
            raise SystemExit(f"{label}: expansion[{s}] の重み和 {tot} が 1 でない")
    cond = None
    if meta and meta.get("condensing_species"):
        cond = _up(meta["condensing_species"])
    elif info and info["condensing"]:
        cond = _up(info["condensing"])
    tracer = bool((meta or {}).get("tracer", {}).get("enabled")) if meta else bool(info and info["tracer"])
    MW = {}
    if meta and meta.get("MW"):
        MW = {_up(k): float(v) for k, v in meta["MW"].items()}
    elif info:
        MW = {_up(k): float(v) for k, v in info["MW"].items()}
    db = None; Tref = 0.0; db_file = None
    if run_dir and os.path.exists(os.path.join(run_dir, "solverConfig.yaml")):
        cfg = yaml.safe_load(open(os.path.join(run_dir, "solverConfig.yaml")))
        pp = cfg.get("physProp") or {}
        Tref = float(pp.get("thermoHrefTemp", 0.0))
        db_file = pp.get("speciesDBFile")
        if db_file:
            p = db_file if os.path.isabs(db_file) else os.path.join(run_dir, db_file)
            if os.path.exists(p):
                db = {_up(k): v for k, v in (yaml.safe_load(open(p)) or {}).items()}
    return {"names": names, "expansion": exp, "condensing": cond, "tracer": tracer, "MW": MW,
            "db": db, "Tref": Tref, "db_file": db_file, "run_dir": run_dir}


def build_transfer(src, dst):
    """T[j, s]: 旧輸送種 s の質量分率のうち新輸送種 j へ行く割合。実種ごとに行き先を 1 つ決める。"""
    # 実種 r → 行き先 j: dst 種 j の expansion に r があれば j (恒等または lump)。
    home = {}
    for j, dj in enumerate(dst["names"]):
        for r in dst["expansion"][dj]:
            if r in home and home[r] != j:
                raise SystemExit(f"destination: 実種 {r} が {dst['names'][home[r]]} と {dj} の両方に入っている")
            home[r] = j
    ns, nd = len(src["names"]), len(dst["names"])
    T = np.zeros((nd, ns))
    missing = []
    for s, ss in enumerate(src["names"]):
        for r, w in src["expansion"][ss].items():
            if r not in home:
                missing.append((ss, r)); continue
            T[home[r], s] += w
    if missing:
        raise SystemExit("REFUSED: source の実種に destination の行き先が無い: "
                         + ", ".join(f"{r} (from {ss})" for ss, r in missing)
                         + f". destination species {dst['names']} に入れるか lump に含めること")
    return T, home


def db_differs(src, dst):
    """DB (係数・MW) か Tref が違えば True (roe 再構成が要る)。片方の DB が読めなければ True。"""
    if src["Tref"] != dst["Tref"]:
        return True, f"thermoHrefTemp {src['Tref']} -> {dst['Tref']}"
    if src["db"] is None or dst["db"] is None:
        return True, "species_db.yaml が片方で読めない"
    if src["names"] != dst["names"]:
        return True, "species set/order changed"
    for n in dst["names"]:
        a, b = src["db"].get(n), dst["db"].get(n)
        if a is None or b is None:
            return True, f"{n} not in one DB"
        for k in ("MW", "nasa9_low", "nasa9_high", "Tmid"):
            if not np.allclose(np.asarray(a[k], float), np.asarray(b[k], float), rtol=0, atol=0):
                return True, f"{n}.{k} differs"
    return False, "same DB and datum"


def T_from_e(gas, Y, e, T0):
    """e = h(T) − R T を Newton で反転 (ベクトル)。"""
    T = np.asarray(T0, float).copy(); R = gas.Rmix(Y)
    for _ in range(60):
        f = gas.h(Y, T) - R * T - e
        dfdT = gas.cp(Y, T) - R
        dT = np.clip(f / np.maximum(dfdT, 1.0), -0.3 * T, 0.3 * T)
        T = np.maximum(T - dT, 10.0)
        if np.max(np.abs(dT)) < 1e-9 * np.max(T):
            break
    return T


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src"); ap.add_argument("dst")
    ap.add_argument("--meta", required=True, help="destination の species_meta.yaml")
    ap.add_argument("--src-meta", help="source の species_meta.yaml (無ければ SRC run dir の solverConfig.yaml)")
    ap.add_argument("--src-run", help="source run dir (既定: SRC h5 の隣)")
    ap.add_argument("--dst-run", help="destination run dir (既定: DST h5 の隣)")
    ap.add_argument("--reconstruct-roe", choices=["auto", "always", "never"], default="auto",
                    help="roe を T から destination DB で再構成する (auto: DB/datum/種集合が変わるときだけ)")
    ap.add_argument("--dry-run", action="store_true", help="書き込まず検査だけ")
    a = ap.parse_args()

    src_run = a.src_run or os.path.dirname(os.path.abspath(a.src))
    dst_run = a.dst_run or os.path.dirname(os.path.abspath(a.dst))
    src = load_layout(a.src_meta, src_run, "source")
    dst = load_layout(a.meta, dst_run, "destination")
    T, home = build_transfer(src, dst)
    print(f"[convert] source species      : {src['names']}")
    print(f"[convert] destination species : {dst['names']}")
    for j, dj in enumerate(dst["names"]):
        parts = [f"{T[j, s]:.6g}*{src['names'][s]}" for s in range(len(src["names"])) if T[j, s] > 0]
        print(f"[convert]   {dj:8s} <- " + (" + ".join(parts) if parts else "0 (not present in source)"))
    if src["condensing"] and dst["condensing"] and src["condensing"] != dst["condensing"]:
        raise SystemExit(f"REFUSED: 凝縮種が違う (source {src['condensing']} / destination {dst['condensing']}); 凝縮モーメントを移せない")
    cond_name = dst["condensing"] or src["condensing"]

    # ---- source 読込 ----
    with h5py.File(a.src, "r") as f:
        V = f["VALUE"]; ro = np.array(V["ro"], np.float64); n = ro.shape[0]
        is_res = "P" in V and "Ux" in V
        if is_res:
            Ux, Uy, Uz = (np.array(V[k], np.float64) for k in ("Ux", "Uy", "Uz"))
            Tsrc = np.array(V["T"], np.float64)
            Ysrc = [np.array(V[f"Y{s}"], np.float64) if f"Y{s}" in V else None for s in range(len(src["names"]))]
            roe = np.array(V["roe"], np.float64) if "roe" in V else None
            roK = ro * np.array(V["k"], np.float64) if "k" in V else None
            roOm = ro * np.array(V["omega"], np.float64) if "omega" in V else None
            moments = {k: ro * np.array(V[k], np.float64) for k in V if k.startswith(("g_", "Q0_", "Q1_", "Q2_"))}
            moments = {"ro" + k: v for k, v in moments.items()}
            roXi = np.array(V["roXi"], np.float64) if "roXi" in V else (ro * np.array(V["Xi"], np.float64) if "Xi" in V else None)
        else:
            roUx, roUy, roUz = (np.array(V[k], np.float64) for k in ("roUx", "roUy", "roUz"))
            Ux, Uy, Uz = roUx / ro, roUy / ro, roUz / ro
            roe = np.array(V["roe"], np.float64)
            Ysrc = [np.array(V[f"roY{s}"], np.float64) / ro if f"roY{s}" in V else None for s in range(len(src["names"]))]
            roK = np.array(V["roK"], np.float64) if "roK" in V else None
            roOm = np.array(V["roOmega"], np.float64) if "roOmega" in V else None
            moments = {k: np.array(V[k], np.float64) for k in V if k.startswith(("rog_", "roQ0_", "roQ1_", "roQ2_"))}
            roXi = np.array(V["roXi"], np.float64) if "roXi" in V else None
            Tsrc = None
        vol = np.array(V["volume"], np.float64) if "volume" in V else None
    if len(src["names"]) == 1 and Ysrc[0] is None:
        Ysrc[0] = np.ones(n)
    for s, y in enumerate(Ysrc):
        if y is None:
            raise SystemExit(f"source に Y{s}/roY{s} ({src['names'][s]}) が無い")
    Ysrc = np.array(Ysrc)                    # [ns, n]
    ssum = Ysrc.sum(axis=0)
    print(f"[convert] source ΣY: min {ssum.min():.9f} max {ssum.max():.9f} (正規化して使う)")
    Ysrc = Ysrc / np.maximum(ssum, 1e-30)

    # source T (input h5 のときは SRC DB で反転)
    ke = 0.5 * (Ux**2 + Uy**2 + Uz**2)
    if Tsrc is None:
        if src["db"] is None:
            raise SystemExit("source が input h5 で T が無く、source の species_db.yaml も読めない")
        gs = _TPGas({k: v for k, v in src["db"].items()}, src["names"], src["Tref"])
        Tsrc = T_from_e(gs, list(Ysrc), roe / ro - ke, np.full(n, 300.0))
        print(f"[convert] source T を SRC DB で反転: {Tsrc.min():.2f}..{Tsrc.max():.2f} K")

    # ---- 種の移送 ----
    Ydst = T @ Ysrc                           # [nd, n]
    dsum = Ydst.sum(axis=0)
    err_sum = np.max(np.abs(dsum - 1.0))
    print(f"[convert] destination ΣY−1: max |{err_sum:.3e}|  {'OK' if err_sum < 1e-6 else 'FAIL (>1e-6)'}")
    if err_sum >= 1e-6:
        raise SystemExit("REFUSED: ΣρY=ρ が 1e-6 で成り立たない (expansion 行列を確認)")
    # 実種ごとの保存 (source 展開 = destination 展開)
    real_src = {}
    for s, ss in enumerate(src["names"]):
        for r, w in src["expansion"][ss].items():
            real_src[r] = real_src.get(r, 0.0) + w * Ysrc[s]
    real_dst = {}
    for j, dj in enumerate(dst["names"]):
        for r, w in dst["expansion"][dj].items():
            real_dst[r] = real_dst.get(r, 0.0) + w * Ydst[j]
    worst = 0.0
    for r in real_src:
        d = np.max(np.abs(real_src[r] - real_dst.get(r, 0.0)) * ro)
        worst = max(worst, d)
    print(f"[convert] 実種ごとの ρY 保存 (per-cell max |Δ|): {worst:.3e} kg/m³")

    # 総水量 (ρY_w + ρg): volume があれば体積重み、無ければ CV 単純和 (同一メッシュなので比較には十分)
    if cond_name:
        w_src = real_src.get(cond_name, np.zeros(n)) * ro
        w_dst = real_dst.get(cond_name, np.zeros(n)) * ro
        g = sum(v for k, v in moments.items() if k.startswith("rog_")) if moments else 0.0
        wt = vol if vol is not None else np.ones(n)
        tot_s = np.sum((w_src + g) * wt); tot_d = np.sum((w_dst + g) * wt)
        print(f"[convert] 総水量 ({cond_name} 気相+液相, {'体積重み' if vol is not None else 'CV 単純和'}): "
              f"source {tot_s:.9e} destination {tot_d:.9e} rel diff {abs(tot_d - tot_s) / max(abs(tot_s), 1e-300):.3e}")

    # ---- roe ----
    differs, why = db_differs(src, dst)
    do_rec = (a.reconstruct_roe == "always") or (a.reconstruct_roe == "auto" and differs)
    gliq = sum(v for k, v in moments.items() if k.startswith("rog_")) / ro if moments else np.zeros(n)
    dry = gliq <= 0.0
    Tchk = None
    if do_rec:
        if dst["db"] is None:
            raise SystemExit("roe 再構成に destination の species_db.yaml が要る (--dst-run)")
        gd = _TPGas({k: v for k, v in dst["db"].items()}, dst["names"], dst["Tref"])
        Yl = list(Ydst)
        e_dst = gd.h(Yl, Tsrc) - gd.Rmix(Yl) * Tsrc
        if src["db"] is not None and roe is not None:
            # 差分形: 液相エネルギー (二相 EOS) を含む source の roe に気相 e の DB 差だけを足す。
            gs = _TPGas({k: v for k, v in src["db"].items()}, src["names"], src["Tref"])
            Ys = list(Ysrc)
            e_src = gs.h(Ys, Tsrc) - gs.Rmix(Ys) * Tsrc
            roe_new = roe + ro * (e_dst - e_src)
            how = "差分形 roe += ρ[e_dst(T) − e_src(T)]"
        else:
            roe_new = ro * (e_dst + ke)
            how = "乾き気相の完全再構成 ρ(e_dst(T)+u²/2) (source DB 無し; 凝縮セルでは液相分が落ちる)"
        de = (roe_new - roe) / ro if roe is not None else np.full(n, np.nan)
        print(f"[convert] roe 再構成 ({why}; {how}): Δe max {np.nanmax(np.abs(de)):.3e} J/kg, mean {np.nanmean(np.abs(de)):.3e} J/kg")
        roe_out = roe_new
    else:
        if roe is None:
            raise SystemExit("source に roe が無く再構成も指定されていない (--reconstruct-roe always)")
        roe_out = roe
        print(f"[convert] roe はそのまま ({why})")
    if dst["db"] is not None:
        # T 保存の検査: 乾きセル (g=0) で destination DB により roe を反転し source の T と比較。湿潤セルは二相 EOS
        # (潜熱) を本ツールは持たないので件数だけ報告する。
        gd = _TPGas({k: v for k, v in dst["db"].items()}, dst["names"], dst["Tref"])
        Tchk = T_from_e(gd, list(Ydst), roe_out / ro - ke, Tsrc)
        dTdry = np.max(np.abs(Tchk - Tsrc)[dry]) if dry.any() else 0.0
        print(f"[convert] T 保存 (destination DB で roe を反転, 乾きセル {int(dry.sum())}): max |ΔT| {dTdry:.3e} K"
              f"{'' if dTdry < 0.05 else '  <-- WARNING: > 0.05 K'}; 湿潤セル {int((~dry).sum())} は二相 EOS のため未検査")

    if a.dry_run:
        print("[convert] --dry-run: 書き込みなし"); return

    # ---- 書き込み (同一メッシュ index コピー) ----
    with h5py.File(a.dst, "r+") as d:
        nd_ = d["VALUE/ro"].shape[0]
        if nd_ != n:
            raise SystemExit(f"REFUSED: CV 数が違う (source {n}, destination {nd_}); 同一メッシュの input h5 を指定する"
                             " (別メッシュは interp_field.py --force-species → 本ツール の順)")
        dt = d["VALUE/ro"].dtype

        def put(name, arr):
            ds = "VALUE/" + name
            if ds in d:
                d[ds][...] = arr.astype(dt)
            else:
                d.create_dataset(ds, data=arr.astype(dt))

        put("ro", ro); put("roUx", ro * Ux); put("roUy", ro * Uy); put("roUz", ro * Uz); put("roe", roe_out)
        if roK is not None: put("roK", roK)
        if roOm is not None: put("roOmega", roOm)
        # 旧 roY を消してから新順序で書く (destination に多い/少ない index が残らないように)
        for k in list(d["VALUE"].keys()):
            if k.startswith("roY") and k[3:].isdigit():
                del d["VALUE/" + k]
        for j in range(len(dst["names"])):
            put(f"roY{j}", ro * Ydst[j])
        for k, v in moments.items():
            put(k, v)
        if roXi is not None:
            put("roXi", roXi)
        moved = ["ro", "roUx", "roUy", "roUz", "roe"] + [f"roY{j}" for j in range(len(dst["names"]))] + list(moments) \
            + (["roK", "roOmega"] if roK is not None else []) + (["roXi"] if roXi is not None else [])
    print(f"[convert] wrote {a.dst}: {moved}")
    print("[convert] SUMMARY: ΣY OK; 実種 ρY 保存 max |Δ| {:.2e}; T 保存 (乾きセル) max |ΔT| {:.2e} K".format(
        worst, float(np.max(np.abs(Tchk - Tsrc)[dry])) if (Tchk is not None and dry.any()) else float("nan")))


if __name__ == "__main__":
    main()
