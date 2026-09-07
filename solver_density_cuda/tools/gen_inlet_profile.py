#!/usr/bin/env python3
"""入口分布 (全温 Tt / 全圧 Pt / 組成 Y_s / k / ω / 超音速入口の ρ,U,Ps) の `inlet_profile_<physID>.csv` を生成・照合する。

forge の入口分布機能 (bcondConfig の対象 inlet に `ints: {inletProfile: 1}`) は、run ディレクトリの
`inlet_profile_<physID>.csv` を face 重心座標で補間して per-face 境界値 (bvar) に書く (kernel 無改修)。
本ツールはその CSV を「座標の式」または「測定表」から作り、run 後に境界ノード/セルの場と照合する。
手順書: procedures/inlet-profile.md。

  gen:    式/表 → CSV
    gen_inlet_profile.py gen --run RUN_DIR --physID 1 --axis y --range -0.0127 0.0127 [--n 201] \
        --Tt "286.65 + 30*exp(-(y/0.004)**2)" --Y H2O="0.005 + 0.012*exp(-(y/0.006)**2)" [--plot]
    gen_inlet_profile.py gen --run RUN_DIR --physID 1 --table measured.csv   (ヘッダ: y Tt Pt Y_H2O k omega ...)
    - 与えなかった量は bcondConfig の一様値のまま (CSV に列を書かない = 変更しない)。
    - 式の変数: 座標 x,y,z (使わない軸は 0), r=sqrt(y²+z²), numpy 関数 (exp, tanh, where, ...),
      bcond 一様値 cfg_Tt, cfg_Pt, cfg_k, cfg_omega, cfg_Y_<name>。
    - 化学種は名前指定 (`--Y H2O=...`)。指定しなかった種が残り (1−ΣY指定) を bcond 比率で受け持つ。
    - inlet_Pressure (亜音速): 列 Tt Pt (+Y_s, k, omega)。Pt の代わりに --Ps と --M でも可 (等エントロピーで Pt に換算)。
    - inlet_uniformVelocity / inlet_fluctVelocity (超音速・全量固定): --Tt と --M と (--Ps | --Pt) から
      ρ, |U|, Ps を換算して列 ro Ux Uy Uz Ps を書く (方向は --dir、既定は bcond の速度方向)。
      CPG は閉形式、TP (thermalMethod 2) は NASA-9 (species_db.yaml) の h/s° で解く (凍結組成)。
    - --set NAME=EXPR で任意の bvar 列 (Ux, ro, Ps, Ts ...) を直接書ける (換算より優先)。
    - 2D 分布は --axis "y z" --range ylo yhi zlo zhi --n ny nz (forge 側は最近傍補間)。

  verify: run 結果と CSV の照合 (境界ノード/owner セルの T0 (VALUE/h0 由来), T, Y_s, k, ω を目標と比較)
    gen_inlet_profile.py verify --run RUN_DIR --physID 1 [--res res_N.h5] [--plot]
"""
import argparse, os, sys, glob
import numpy as np, yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RU = 8.31446261815324
AXIS = {"x": 0, "y": 1, "z": 2}


# ----------------------------------------------------------------------------- config
def load_run(run_dir):
    cfg = yaml.safe_load(open(os.path.join(run_dir, "solverConfig.yaml")))
    bcs = yaml.safe_load(open(os.path.join(run_dir, "bcondConfig.yaml")))
    return cfg, bcs


def find_bc(bcs, physID):
    for name, b in bcs.items():
        if isinstance(b, dict) and int(b.get("physID", -1)) == physID:
            return name, b
    raise SystemExit(f"bcondConfig.yaml に physID={physID} が無い")


def species_names(cfg):
    """physProp.species (TP)。省略時はソルバ既定と同じ単成分 N2 (solverConfig.cpp)。CPG は []。"""
    pp = cfg["physProp"]
    if int(pp.get("thermalMethod", 0)) == 2:
        return list(pp["species"]) if pp.get("species") else ["N2"]
    return []


class Gas:
    """CPG (thermalMethod 0) / TP NASA-9 凍結組成 (thermalMethod 2) の h, cp, s°, R, γ (質量基準)。"""
    def __init__(self, run_dir, cfg):
        pp = cfg["physProp"]; self.tm = int(pp.get("thermalMethod", 0))
        if self.tm == 2:
            from total_quantities import _TPGas
            self.names = species_names(cfg)
            db = yaml.safe_load(open(os.path.join(run_dir, pp.get("speciesDBFile", "species_db.yaml"))))
            for n in self.names:
                if n not in db:
                    raise SystemExit(f"species_db に {n} が無い (physProp.species / 既定 N2)")
            self.tp = _TPGas(db, self.names, float(pp.get("thermoHrefTemp", 0.0)))
        elif self.tm == 0:
            self.cp0 = float(pp["cp"]); self.ga0 = float(pp["gamma"]); self.R0 = self.cp0 * (self.ga0 - 1.0) / self.ga0
            self.names = []
        else:
            raise SystemExit(f"thermalMethod {self.tm} は未対応")

    # Y: list of arrays (nSpecies) or None
    def R(self, Y):
        return self.tp.Rmix(Y) if self.tm == 2 else self.R0

    def h(self, Y, T):
        return self.tp.h(Y, T) if self.tm == 2 else self.cp0 * T

    def cp(self, Y, T):
        return self.tp.cp(Y, T) if self.tm == 2 else np.full_like(np.asarray(T, float), self.cp0)

    def s0(self, Y, T):
        return self.tp.s0(Y, T) if self.tm == 2 else self.cp0 * np.log(T)

    def gamma(self, Y, T):
        cp = self.cp(Y, T); R = self.R(Y)
        return cp / (cp - R)

    def a2(self, Y, T):  # 音速²
        return self.gamma(Y, T) * self.R(Y) * T

    def static_from_total(self, Y, Tt, M):
        """h(Tt) - h(Ts) = M² a(Ts)²/2 を Ts について解く (単調・二分法, ベクトル)。"""
        Tt = np.asarray(Tt, float); M = np.asarray(M, float)
        if self.tm == 0:
            return Tt / (1.0 + 0.5 * (self.ga0 - 1.0) * M * M)
        lo = Tt / (1.0 + 0.5 * 0.7 * M * M) * 0.5; hi = Tt.copy()
        hT = self.h(Y, Tt)
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            f = hT - self.h(Y, mid) - 0.5 * M * M * self.a2(Y, mid)   # Ts↑ で単調減少
            lo = np.where(f > 0, mid, lo); hi = np.where(f > 0, hi, mid)
        return 0.5 * (lo + hi)

    def p_ratio_t_over_s(self, Y, Tt, Ts):
        """Pt/Ps = exp((s°(Tt) − s°(Ts))/R) (等エントロピー)。"""
        return np.exp((self.s0(Y, Tt) - self.s0(Y, Ts)) / self.R(Y))


# ----------------------------------------------------------------------------- expressions
def make_namespace(coords, axes, bc_floats, names):
    ns = {k: getattr(np, k) for k in ("exp", "log", "sqrt", "tanh", "sin", "cos", "abs", "where",
                                        "minimum", "maximum", "clip", "pi", "heaviside")}
    ns["np"] = np
    xyz = [np.zeros(coords.shape[0]) for _ in range(3)]
    for j, ax in enumerate(axes):
        xyz[AXIS[ax]] = coords[:, j]
    ns.update(x=xyz[0], y=xyz[1], z=xyz[2], r=np.sqrt(xyz[1] ** 2 + xyz[2] ** 2))
    for k, v in (bc_floats or {}).items():
        ns["cfg_" + k] = float(v)
    for s, n in enumerate(names):
        ns["cfg_Y_" + n] = float((bc_floats or {}).get(f"Y{s}", 1.0 if s == 0 else 0.0))
    return ns


def evaluate(expr, ns, npts):
    v = eval(expr, {"__builtins__": {}}, ns)
    v = np.asarray(v, dtype=float)
    if v.ndim == 0:
        v = np.full(npts, float(v))
    if v.shape[0] != npts:
        raise SystemExit(f"式 '{expr}' の長さ {v.shape[0]} が点数 {npts} と合わない")
    return v


def parse_kv(items):
    out = {}
    for it in items or []:
        if "=" not in it:
            raise SystemExit(f"NAME=EXPR 形式で指定: {it}")
        k, e = it.split("=", 1); out[k.strip()] = e.strip()
    return out


# ----------------------------------------------------------------------------- gen
def cmd_gen(a):
    run_dir = os.path.abspath(a.run)
    cfg, bcs = load_run(run_dir)
    bcname, bc = find_bc(bcs, a.physID)
    kind = bc["kind"]; bcf = bc.get("floats") or {}
    names = species_names(cfg)
    gas = Gas(run_dir, cfg)

    # ---- sample points ----
    table = {}
    if a.table:
        with open(a.table) as f:
            hdr = None; rows = []
            for line in f:
                t = line.split()
                if not t or t[0].startswith("#"):
                    continue
                if hdr is None:
                    hdr = t; continue
                rows.append([float(v) for v in t])
        rows = np.array(rows)
        axes = [h for h in hdr if h in AXIS]
        if not axes or hdr[:len(axes)] != axes:
            raise SystemExit("--table のヘッダは先頭に座標列 (x/y/z) が要る")
        coords = rows[:, :len(axes)]
        for j, h in enumerate(hdr[len(axes):]):
            table[h] = rows[:, len(axes) + j]
    else:
        axes = a.axis.split()
        if len(a.range) != 2 * len(axes):
            raise SystemExit("--range は軸ごとに lo hi")
        ns_ = a.n if len(a.n) == len(axes) else [a.n[0]] * len(axes)
        grids = [np.linspace(a.range[2 * j], a.range[2 * j + 1], ns_[j]) for j in range(len(axes))]
        mesh = np.meshgrid(*grids, indexing="ij")
        coords = np.stack([m.ravel() for m in mesh], axis=1)
    npts = coords.shape[0]
    ns = make_namespace(coords, axes, bcf, names)

    def get(name, expr, default=None):
        """列の値: --式 > 表の列 > default (None なら未指定)。"""
        if expr is not None:
            return evaluate(expr, ns, npts)
        if name in table:
            return table[name]
        return None

    Tt = get("Tt", a.Tt); Pt = get("Pt", a.Pt); Ps = get("Ps", a.Ps); M = get("M", a.M)
    k = get("k", a.k); om = get("omega", a.omega)
    direct = {n: evaluate(e, ns, npts) for n, e in parse_kv(a.set).items()}
    for n in list(table.keys()):
        if n not in ("Tt", "Pt", "Ps", "M", "k", "omega") and not n.startswith("Y_") and n not in direct:
            direct[n] = table[n]

    # ---- species ----
    Ycols = None
    yspec = parse_kv(a.Y)
    for n in table:
        if n.startswith("Y_") and n[2:] not in yspec:
            yspec[n[2:]] = None
    if yspec:
        if not names:
            raise SystemExit("--Y は thermalMethod 2 (physProp.species) のときだけ使える")
        Ycols = [None] * len(names)
        for nm, e in yspec.items():
            if nm not in names:
                raise SystemExit(f"化学種 {nm} は physProp.species {names} に無い")
            Ycols[names.index(nm)] = evaluate(e, ns, npts) if e is not None else table["Y_" + nm]
        for nm, c in zip(names, Ycols):
            if c is not None and (not np.all(np.isfinite(c)) or np.any(c < -1e-12) or np.any(c > 1.0 + 1e-9)):
                raise SystemExit(f"Y_{nm} が [0,1] を外れるか非有限: min {np.nanmin(c):.4g} max {np.nanmax(c):.4g}")
        given = np.sum([c for c in Ycols if c is not None], axis=0)
        if np.any(given > 1.0 + 1e-9) or np.any(given < -1e-12):
            raise SystemExit(f"指定した Y の和が [0,1] を外れる: min {given.min():.4g} max {given.max():.4g}")
        rest_idx = [s for s, c in enumerate(Ycols) if c is None]
        if rest_idx:
            w = np.array([float(bcf.get(f"Y{s}", 1.0 if s == 0 else 0.0)) for s in rest_idx])
            w = w / w.sum() if w.sum() > 0 else np.full(len(rest_idx), 1.0 / len(rest_idx))
            for wi, s in zip(w, rest_idx):
                Ycols[s] = wi * (1.0 - given)
        elif np.any(np.abs(given - 1.0) > 1e-6):
            raise SystemExit("全化学種を指定したが和が 1 でない")
    # 換算に使う組成 (指定が無ければ bcond 一様値)
    Yuse = Ycols if Ycols is not None else (
        [np.full(npts, float(bcf.get(f"Y{s}", 1.0 if s == 0 else 0.0))) for s in range(len(names))] or None)

    cols = {}
    if kind in ("inlet_Pressure", "inlet_Pressure_dir"):
        if Tt is None and Pt is None and Ps is None and Ycols is None and k is None and om is None and not direct:
            raise SystemExit("何も分布が指定されていない (--Tt/--Pt/--Y/--k/--omega/--set)")
        Tt_use = Tt if Tt is not None else np.full(npts, float(bcf["Tt"]))
        if Pt is None and Ps is not None:
            if M is None:
                raise SystemExit("--Ps から Pt を作るには --M が要る")
            Ts = gas.static_from_total(Yuse, Tt_use, M); Pt = Ps * gas.p_ratio_t_over_s(Yuse, Tt_use, Ts)
        if Tt is not None: cols["Tt"] = Tt
        if Pt is not None: cols["Pt"] = Pt
    elif kind in ("inlet_uniformVelocity", "inlet_fluctVelocity"):
        if Tt is not None or M is not None:
            if Tt is None or M is None or (Ps is None and Pt is None):
                raise SystemExit("超音速入口の換算には --Tt, --M と (--Ps | --Pt) の 3 つが要る")
            Ts = gas.static_from_total(Yuse, Tt, M)
            if Ps is None:
                Ps = Pt / gas.p_ratio_t_over_s(Yuse, Tt, Ts)
            R = gas.R(Yuse); ro = Ps / (R * Ts); U = M * np.sqrt(gas.a2(Yuse, Ts))
            if a.dir is not None:
                d = np.array(a.dir, float)
            else:
                d = np.array([float(bcf.get("Ux", 1.0)), float(bcf.get("Uy", 0.0)), float(bcf.get("Uz", 0.0))])
            d = d / np.linalg.norm(d)
            pre = "Ux0" if kind == "inlet_fluctVelocity" else "Ux"
            cols["ro"] = ro; cols["Ps"] = Ps
            cols[pre] = U * d[0]; cols[pre.replace("x", "y")] = U * d[1]; cols[pre.replace("x", "z")] = U * d[2]
            # 換算の自己検証: 出力から Tt を戻す
            Ts_chk = Ps / (ro * R); M_chk = U / np.sqrt(gas.a2(Yuse, Ts_chk))
            Tt_chk = np.array([0.0])
            hs = gas.h(Yuse, Ts_chk) + 0.5 * U * U
            T0 = Ts_chk.copy()
            for _ in range(50):
                r = gas.h(Yuse, T0) - hs; T0 = np.maximum(T0 - r / np.maximum(gas.cp(Yuse, T0), 1.0), 1.0)
            print(f"[gen] 換算チェック: Tt 復元誤差 max {np.max(np.abs(T0 - Tt)):.3e} K, M 復元誤差 max {np.max(np.abs(M_chk - M)):.3e}")
        elif not direct and Ycols is None and k is None and om is None:
            raise SystemExit("何も分布が指定されていない")
    else:
        raise SystemExit(f"kind {kind} は入口ではない (inletProfile は inlet_* のみ)")

    if Ycols is not None:
        for s, c in enumerate(Ycols):
            cols[f"Y{s}"] = c
    if k is not None: cols["k"] = k
    if om is not None: cols["omega"] = om
    cols.update(direct)

    out = a.out or os.path.join(run_dir, f"inlet_profile_{a.physID}.csv")
    with open(out, "w") as f:
        f.write(" ".join(axes + list(cols.keys())) + "\n")
        for i in range(npts):
            f.write(" ".join(f"{v:.9g}" for v in coords[i]) + " " + " ".join(f"{cols[c][i]:.9g}" for c in cols) + "\n")
    print(f"[gen] {out}: kind={kind} ({bcname}), {npts} 点, 補間 {'1D 線形' if len(axes) == 1 else f'{len(axes)}D 最近傍'} ({' '.join(axes)})")
    for c, v in cols.items():
        tag = ""
        if c in bcf:
            tag = f"  (bcond 一様値 {float(bcf[c]):.6g})"
        print(f"  {c:>6s}: min {v.min():.6g}  max {v.max():.6g}{tag}")
    ints = bc.get("ints") or {}
    if int(ints.get("inletProfile", 0)) != 1:
        print(f"[gen] 注意: bcondConfig.yaml の {bcname} に `ints: {{inletProfile: 1}}` が無い (CSV は読まれない)")
    if a.plot:
        plot_profiles(coords, axes, cols, os.path.splitext(out)[0] + ".png", title=f"{bcname} (physID {a.physID})")
    return 0


def plot_profiles(coords, axes, cols, png, title="", overlay=None):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    if len(axes) != 1:
        print("[plot] 2D/3D 分布の図化は未対応 (1D のみ)"); return
    n = len(cols); fig, axs = plt.subplots(1, n, figsize=(3.6 * n, 3.4), squeeze=False)
    s = coords[:, 0]; o = np.argsort(s)
    for ax, (c, v) in zip(axs[0], cols.items()):
        ax.plot(v[o], s[o], "-", label="target (csv)")
        if overlay and c in overlay:
            ss, vv = overlay[c]; oo = np.argsort(ss); ax.plot(vv[oo], ss[oo], ".", ms=3, label="field")
        ax.set_xlabel(c); ax.set_ylabel(axes[0]); ax.grid(alpha=.3)
        if overlay and c in overlay: ax.legend(fontsize=7)
    fig.suptitle(title); fig.tight_layout(); fig.savefig(png, dpi=130); print(f"[plot] {png}")


# ----------------------------------------------------------------------------- verify
def interp_like_forge(coords_tab, vals_tab, axes, fc):
    """forge applyInletProfiles と同じ規則: 1D は昇順ソート線形 (端はクランプ)、2D/3D は最近傍。"""
    if len(axes) == 1:
        s = coords_tab[:, 0]; o = np.argsort(s)
        return np.interp(fc[:, AXIS[axes[0]]], s[o], vals_tab[o])
    idx = [AXIS[ax] for ax in axes]
    d = ((fc[:, None, idx] - coords_tab[None, :, :]) ** 2).sum(axis=2)
    return vals_tab[np.argmin(d, axis=1)]


def cmd_verify(a):
    import h5py
    run_dir = os.path.abspath(a.run)
    cfg, bcs = load_run(run_dir)
    bcname, bc = find_bc(bcs, a.physID)
    names = species_names(cfg)
    csv = a.csv or os.path.join(run_dir, f"inlet_profile_{a.physID}.csv")
    with open(csv) as f:
        hdr = f.readline().split(); rows = np.loadtxt(f, ndmin=2)
    axes = [h for h in hdr if h in AXIS]; nax = len(axes)
    ctab = rows[:, :nax]; tab = {h: rows[:, nax + j] for j, h in enumerate(hdr[nax:])}
    mesh_path = os.path.join(run_dir, cfg["mesh"]["meshFileName"])
    with h5py.File(mesh_path, "r") as m:
        pc = m["PLANES/centCoords"][:].reshape(-1, 3)
        ip = m[f"BCONDS/{a.physID}/iPlanes"][:]; ic = m[f"BCONDS/{a.physID}/iCells"][:]
        fc = pc[ip]
    res = a.res or sorted(glob.glob(os.path.join(run_dir, "res_[0-9]*.h5")), key=lambda s: int(s.split("_")[-1][:-3]))[-1]
    node = cfg["mesh"].get("discretization", "cell") == "node"
    print(f"[verify] {bcname} (physID {a.physID}, kind {bc['kind']}), {len(ip)} faces, res={os.path.basename(res)}, "
          f"{'node (境界ノード値)' if node else 'cell (第 1 セル値: 境界値そのものではない)'}")
    field = {}
    with h5py.File(res, "r") as r:
        V = r["VALUE"]
        for key in ("T", "P", "k", "omega", "Ux", "Uy", "Uz", "ro"):
            if key in V: field[key] = V[key][:][ic].astype(float)
        for s in range(len(names)):
            if f"Y{s}" in V: field[f"Y{s}"] = V[f"Y{s}"][:][ic].astype(float)
        has_h0 = "h0" in V
    if has_h0:
        from total_quantities import total_state
        st = total_state(run_dir, res); field["Tt"] = st["T0"][ic]; field["Pt"] = st["P0"][ic]
    # 超音速入口 (ro/U/Ps 列) は換算元の Tt, M を目標として復元し、h0 由来 T0 と比較する
    targets = {c: interp_like_forge(ctab, v, axes, fc) for c, v in tab.items()}
    ukeys = ("Ux", "Uy", "Uz") if "Ux" in targets else ("Ux0", "Uy0", "Uz0")
    if bc["kind"] in ("inlet_uniformVelocity", "inlet_fluctVelocity") and "ro" in targets and "Ps" in targets and has_h0:
        gas = Gas(run_dir, cfg)
        Yt = [targets[f"Y{s}"] for s in range(len(names))] if names and all(f"Y{s}" in targets for s in range(len(names))) else (
            [np.full(len(fc), float((bc.get("floats") or {}).get(f"Y{s}", 1.0 if s == 0 else 0.0))) for s in range(len(names))] or None)
        U2 = sum(targets.get(k, np.zeros(len(fc))) ** 2 for k in ukeys)
        Ts = targets["Ps"] / (targets["ro"] * gas.R(Yt)); hs = gas.h(Yt, Ts) + 0.5 * U2; T0 = Ts.copy()
        for _ in range(50):
            T0 = np.maximum(T0 - (gas.h(Yt, T0) - hs) / np.maximum(gas.cp(Yt, T0), 1.0), 1.0)
        targets["Tt"] = T0; targets["M"] = np.sqrt(U2 / gas.a2(Yt, Ts))
        if "Ux" in field:
            field["M"] = np.sqrt(sum(field.get(k, 0.0) ** 2 for k in ("Ux", "Uy", "Uz"))) / np.sqrt(
                gas.a2(Yt, field["T"])) if "T" in field else None
    # CSV 列名 → res の場名 (静圧 Ps は P, 静温 Ts は T)。超音速入口の Ps は境界で固定されるので第 1 ノードと比較可。
    alias = {"Ps": "P", "Ts": "T", "Ux0": "Ux", "Uy0": "Uy", "Uz0": "Uz"}
    # 目標が全点 0 の列 (Uy=0 等) は同じ単位群 (速度成分) の最大値でスケールし、実測の非ゼロを隠さない
    groups = [("Ux", "Uy", "Uz", "Ux0", "Uy0", "Uz0")]
    def unit_scale(c):
        for g in groups:
            if c in g:
                return max([abs(targets[k]).max() for k in g if k in targets] + [1e-30])
        return 1e-30
    overlay = {}; worst = 0.0
    print(f"  {'col':>6s} {'target min':>12s} {'target max':>12s} {'field min':>12s} {'field max':>12s} {'max|Δ|':>10s} {'max|Δ|/scale':>12s}")
    for c, tgt in targets.items():
        fv = field.get(c, field.get(alias.get(c, "")))
        if fv is None:
            print(f"  {c:>6s} {tgt.min():12.6g} {tgt.max():12.6g}   (res に対応する場が無い)"); continue
        d = np.abs(fv - tgt); rng = tgt.max() - tgt.min(); amax = abs(tgt).max()
        if rng > 1e-3 * amax:
            scale = rng                      # 分布幅基準
        elif amax > 0:
            scale = amax                     # ほぼ一様な列は絶対値基準
        else:
            scale = unit_scale(c)            # 全点 0 の列は同じ単位群の大きさ基準
        rel = d.max() / scale
        worst = max(worst, rel)
        print(f"  {c:>6s} {tgt.min():12.6g} {tgt.max():12.6g} {fv.min():12.6g} {fv.max():12.6g} {d.max():10.3g} {rel:12.3g}")
        overlay[c] = (fc[:, AXIS[axes[0]]] if nax == 1 else np.arange(len(fc)), fv)
    note = "(Tt は VALUE/h0 由来の T0; 亜音速 inlet_Pressure の境界 T は Tt より u²/2cp だけ低いのが正常)" if has_h0 else \
           "(VALUE/h0 が無いので Tt は比較不可: output.level>=1 で出力する)"
    print(f"[verify] 目標との最大相対差 (分布幅基準, 一様列は絶対値基準) = {worst:.3g}  {note}")
    if a.plot:
        cols = {c: t for c, t in targets.items() if field.get(c) is not None}
        plot_profiles(fc[:, [AXIS[axes[0]]]] if nax == 1 else fc, axes if nax == 1 else ["x"], cols,
                      os.path.join(run_dir, f"inlet_profile_{a.physID}_verify.png"),
                      title=f"{bcname} target vs {os.path.basename(res)}", overlay=overlay)
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gen"); g.add_argument("--run", required=True); g.add_argument("--physID", type=int, required=True)
    g.add_argument("--axis", default="y", help='補間軸: "y" / "x" / "z" / "y z" (2D は最近傍)')
    g.add_argument("--range", type=float, nargs="+", default=[], help="軸ごとの lo hi")
    g.add_argument("--n", type=int, nargs="+", default=[201]); g.add_argument("--table", help="測定表 CSV (空白区切り, ヘッダ先頭に座標列)")
    g.add_argument("--Tt"); g.add_argument("--Pt"); g.add_argument("--Ps"); g.add_argument("--M")
    g.add_argument("--Y", action="append", help="NAME=EXPR (例 H2O=0.01+0.005*exp(-(y/0.004)**2))")
    g.add_argument("--k"); g.add_argument("--omega"); g.add_argument("--set", action="append", help="NAME=EXPR (任意 bvar 列)")
    g.add_argument("--dir", type=float, nargs=3, help="超音速入口の速度方向 (既定: bcond の Ux,Uy,Uz)")
    g.add_argument("--out"); g.add_argument("--plot", action="store_true")
    v = sub.add_parser("verify"); v.add_argument("--run", required=True); v.add_argument("--physID", type=int, required=True)
    v.add_argument("--res"); v.add_argument("--csv"); v.add_argument("--plot", action="store_true")
    a = ap.parse_args()
    if a.cmd == "gen" and not a.table and not a.range:
        ap.error("--range (または --table) が要る")
    return cmd_gen(a) if a.cmd == "gen" else cmd_verify(a)


if __name__ == "__main__":
    sys.exit(main())
