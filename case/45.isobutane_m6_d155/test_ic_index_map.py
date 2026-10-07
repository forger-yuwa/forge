"""ic_index_map (plan tooling-nozzle-throat-monotone-r2 §6 E1、§5.1 #5b) の試験。合成した小さな構造格子 (ni 6 × nj 6、node 方式の
nozzle.h5 / res_*.h5 と同じデータセット構成、CPG) で、正例 (番号写像・最近傍・再適用) と負例 (接続違い・j ずれ・境界の節点ずれ・
境界種別違い・移動超過・反転・量欠落・SRC の保存量が DST に無い・非有限・physProp 違い・単位 (scale_m) 違い・保存場の座標違い) を調べる。
負例では宛先 h5 が 1 バイトも変わらず、記録の VERDICT が REFUSED であることも確かめる。
格子は壁際の層を 0.3 µm・0.7 µm にし、2 列で壁際 3 層を内向きに 0.5 µm 動かすと最近傍が 1 層内側 (j − 1) を拾う (実格子の 37 節点と同じ型)。
usage: python3 test_ic_index_map.py   (FAIL 0 で exit 0)
"""
import json
import shutil
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ic_index_map as M  # noqa: E402

FAIL = 0
NI, NJ = 6, 6
SCALE = 0.003
R_COL = np.array([0.0, 1e-3, 2e-3, 3e-3 - 1.0e-6, 3e-3 - 0.3e-6, 3e-3])   # 壁際: 0.3 µm と 0.7 µm の層
CFG = """mesh: {discretization: "node", isAxisymmetric: 1, axisCentroidShift: 1, meshFileName: "nozzle.h5", valueFileName: "nozzle.h5"}
physProp: {thermalMethod: 0, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: 1004.5, gamma: 1.4}
turbulence: {model: "none"}
"""
KINDS = {"1": "inlet_Pressure", "2": "outlet_statPress", "3": "slip", "4": "axis"}
CONS = ("ro", "roUx", "roUy", "roUz", "roe")


def check(name, cond):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name)
    if not cond:
        FAIL += 1


def grid(dR=None):
    X = np.repeat((np.arange(NI) * 1e-3)[:, None], NJ, axis=1)
    R = np.repeat(R_COL[None, :], NI, axis=0)
    if dR is not None:
        R = R + dR
    return X, R


def viz_conne(ni, nj):
    q = []
    for i in range(ni - 1):
        for j in range(nj - 1):
            b = i * nj + j
            q += [5, b, b + nj, b + nj + 1, b + 1]
    return np.array(q, np.int32)


def edge_nodes(edge):
    idx = np.arange(NI * NJ).reshape(NI, NJ)
    e = {"1": idx[0, :], "2": idx[-1, :], "3": idx[:, -1], "4": idx[:, 0]}[edge]
    return np.stack([e[:-1], e[1:]], axis=1).ravel().astype(np.int32)


def write_mesh(path, X, R, values, kinds=KINDS, bnodes=None):
    n = NI * NJ
    with h5py.File(path, "w") as f:
        f.create_dataset("MESH/COORD", data=np.stack([X.ravel(), R.ravel(), np.zeros(n)], 1).astype(np.float32).ravel())
        f.create_dataset("MESH/CONNE", data=np.stack([np.full(n, 4), np.arange(n)], 1).astype(np.int32).ravel())
        f["MESH"].attrs["nNodes"] = np.int32(n)
        f["MESH"].attrs["nCells"] = np.int32(n)
        f.create_dataset("VIZMESH/CONNE", data=viz_conne(NI, NJ))
        f.create_dataset("CELLS/STRUCT", data=np.arange(3 * n, dtype=np.int32))
        f.create_dataset("CELLS/regionId", data=np.zeros(n, np.int32))
        f.create_dataset("CELLS/volume", data=(1.0 + np.arange(n)).astype(np.float32) * 1e-9)     # 幾何量 (保持されること)
        f.create_dataset("PLANES/STRUCT", data=np.arange(5 * n, dtype=np.int32))
        f.create_dataset("PLANES/surfArea", data=np.linspace(1e-7, 2e-7, 2 * n).astype(np.float32))
        for b, k in kinds.items():
            g = f.create_group(f"BCONDS/{b}")
            g.attrs["bcondKind"] = k
            nodes = (bnodes or {}).get(b, edge_nodes(b))
            g.create_dataset("vizBfaceNodes", data=nodes)
            g.create_dataset("vizBfaceSizes", data=np.full(nodes.size // 2, 2, np.int32))
            for t in ("iBPlanes", "iCells", "iPlanes"):
                g.create_dataset(t, data=np.arange(nodes.size // 2 + 1, dtype=np.int32))
            g.create_dataset("VALUE/Ps", data=np.full(nodes.size // 2, 101325.0, np.float32))
        for k, v in values.items():
            f.create_dataset("VALUE/" + k, data=np.asarray(v, np.float32))


def src_values(seed=0):
    rng = np.random.default_rng(seed)
    n = NI * NJ
    v = {k: (1.0 + rng.random(n)).astype(np.float32) * s for k, s in zip(CONS, (1.0, 300.0, 50.0, 0.0, 2e5))}
    v["roUz"] = np.zeros(n, np.float32)
    return v


def setup(tmp, dR_dst=None, mutate=None):
    """tmp/src (IC run: nozzle.h5 + res_100.h5) と tmp/dst (宛先 prep: nozzle.h5) を作る。mutate(paths) で負例の改変。"""
    src, dst = tmp / "src", tmp / "dst"
    for d in (src, dst):
        d.mkdir(parents=True)
        (d / "solverConfig.yaml").write_text(CFG)
        (d / "prepare_info.json").write_text(json.dumps({"scale_m": SCALE, "mesh": {"ni": NI, "nj": NJ}}))
    X, R = grid()
    sv = src_values()
    wd = np.linspace(0.0, 1e-3, NI * NJ).astype(np.float32)
    write_mesh(src / "nozzle.h5", X, R, {**{k: v * 0 + 1 for k, v in sv.items()}, "wall_dist": wd})
    with h5py.File(src / "res_100.h5", "w") as f:
        f.create_dataset("MESH/COORD", data=np.stack([X.ravel(), R.ravel(), np.zeros(NI * NJ)], 1).astype(np.float32).ravel())
        f.create_dataset("MESH/CONNE", data=viz_conne(NI, NJ))
        for k, v in sv.items():
            f.create_dataset("VALUE/" + k, data=v)
        f.create_dataset("VALUE/P", data=sv["roe"] * 0.4)
        f.create_dataset("VALUE/Ux", data=sv["roUx"] / sv["ro"])
        f.create_dataset("VALUE/wall_dist", data=wd)
        f.create_dataset("VALUE/h0", data=sv["roe"] / sv["ro"]).attrs["h0_includes_k"] = np.int32(0)
    Xd, Rd = grid(dR_dst)
    write_mesh(dst / "nozzle.h5", Xd, Rd, {**{k: np.full(NI * NJ, 7.0) for k in CONS}, "wall_dist": wd * 1.01})
    paths = {"src": src, "dst": dst, "res": src / "res_100.h5", "src_mesh": src / "nozzle.h5", "dst_h5": dst / "nozzle.h5", "sv": sv}
    if mutate:
        mutate(paths)
    return paths


def shift_dR(cols=(2, 3), layers=(3, 4, 5), d=-0.5e-6):
    dR = np.zeros((NI, NJ))
    for i in cols:
        for j in layers:
            dR[i, j] = d
    return dR


def run(p, **kw):
    return M.apply_map(p["res"], p["dst_h5"], p["src_mesh"], p["src"], p["dst"], record=p["dst"] / "IC_MAP.json", **kw)


def expect_refused(name, dR=None, mutate=None, want=(), **kw):
    with tempfile.TemporaryDirectory() as t:
        p = setup(Path(t), dR, mutate)
        sha0 = M._sha_file(p["dst_h5"])
        try:
            run(p, **kw)
        except M.MapRefused as e:
            failed = {k for k, v in e.record["checks"].items() if not v["ok"]}
            rec = json.loads((p["dst"] / "IC_MAP.json").read_text())
            check(f"{name}: 拒否 (不成立 {sorted(failed)}) / 期待 {sorted(want)} を含む", set(want) <= failed)
            check(f"{name}: 宛先 h5 は 1 バイトも変わらない・記録は REFUSED",
                  M._sha_file(p["dst_h5"]) == sha0 and rec["VERDICT"] == "REFUSED" and rec["failures"])
            return failed
        check(f"{name}: 拒否されるはずが通った", False)


# --- 正例: 番号写像 --------------------------------------------------------------------------------------------------
with tempfile.TemporaryDirectory() as t:
    p = setup(Path(t), shift_dR())
    with h5py.File(p["dst_h5"]) as f:
        geo0 = {k: f[k][:].tobytes() for k in ("CELLS/volume", "PLANES/surfArea", "VALUE/wall_dist", "MESH/COORD")}
    out = run(p, resolve_species=True)       # CPG: 化学種の属性は対象外 (forge を起動しない)
    check(f"番号写像: VERDICT OK ({out['VERDICT']})", out["VERDICT"] == "OK")
    with h5py.File(p["dst_h5"]) as f:
        same = all(f["VALUE/" + k][:].tobytes() == p["sv"][k].tobytes() for k in CONS)
        geo1 = {k: f[k][:].tobytes() for k in geo0}
    check("番号写像: 転送量が SRC とビット一致 (保存量を直接コピー)", same)
    check("番号写像: wall_dist・幾何量・座標は宛先のまま", geo0 == geo1)
    check(f"番号写像: 転送した量 = 宛先の /VALUE − wall_dist ({out['transferred']})", sorted(out["transferred"]) == sorted(CONS))
    dr = out["checks"]["displacement"]["detail"]
    check(f"番号写像: 移動の最大 {dr['max_m']:.3e} m ≈ 0.5 µm、動いた節点 6 ({dr['n_moved']})",
          abs(dr["max_m"] - 0.5e-6) < 2e-9 and dr["n_moved"] == 6 and dr["argmax"]["i"] in (2, 3))
    nv = out["nearest_vs_index"]
    check(f"番号写像の記録: 最近傍は 4 節点で 1 層内側 (j − 1) を拾う ({nv['n_mismatch']}, {nv['dj_counts']})",
          nv["n_mismatch"] == 4 and nv["dj_counts"] == {"-1": 4} and {r["layer_from_wall"] for r in nv["nodes"]} == {0, 1})
    check("記録: 全検査が ok", all(v["ok"] for v in out["checks"].values()))
    out2 = run(p, resolve_species=True)
    check("再適用 (同じ宛先にもう一度) も OK でビット一致", out2["VERDICT"] == "OK" and out2["post_write"]["bit_exact_vs_src"])
    rc = M.main([str(p["res"]), str(p["dst_h5"]), "--src-mesh", str(p["src_mesh"]), "--dst-run", str(p["dst"])])
    check(f"CLI: 成立なら exit 0 ({rc})", rc == 0)

# --- 正例: 最近傍 (予備 A/B の α) ------------------------------------------------------------------------------------
with tempfile.TemporaryDirectory() as t:
    p = setup(Path(t), shift_dR())
    out = run(p, mode="nearest", resolve_species=False)
    from scipy.spatial import cKDTree
    X, R = grid(); Xd, Rd = grid(shift_dR())
    cs = np.stack([X.ravel(), R.ravel()], 1).astype(np.float32).astype(float)
    cd = np.stack([Xd.ravel(), Rd.ravel()], 1).astype(np.float32).astype(float)
    idx = cKDTree(cs).query(cd)[1]
    with h5py.File(p["dst_h5"]) as f:
        same = all(f["VALUE/" + k][:].tobytes() == p["sv"][k][idx].tobytes() for k in CONS)
        nd = sum(int(np.count_nonzero(f["VALUE/" + k][:] != p["sv"][k])) for k in ("ro",))
    check(f"最近傍: VERDICT OK、転送量が SRC[最近傍] とビット一致 ({out['VERDICT']})", out["VERDICT"] == "OK" and same)
    check(f"最近傍: 番号写像と違う値の節点は 4 ({nd})", nd == 4)
    d = out["nearest_vs_index"]["conserved_diff_nearest_minus_index"]["ro"]
    want = float(np.abs(p["sv"]["ro"][idx].astype(float) - p["sv"]["ro"].astype(float)).max())
    check(f"最近傍: 保存量の差 (最大) を記録 ({d['max_abs']:.6g} = {want:.6g})、RMS > 0", d["max_abs"] == want and d["rms"] > 0 and d["n_differs"] == 4)
    check("記録: 乾式確認では属性を付けない (species_resolved False)", out["species_resolved"] is False)

# --- 負例 ----------------------------------------------------------------------------------------------------------------
def _cells_struct(p):
    with h5py.File(p["dst_h5"], "r+") as f:
        a = f["CELLS/STRUCT"][:]
        a[3] += 1
        f["CELLS/STRUCT"][...] = a


expect_refused("接続違い (CELLS/STRUCT)", shift_dR(), _cells_struct, want={"connectivity"})


def _viz(p):
    with h5py.File(p["dst_h5"], "r+") as f:
        a = f["VIZMESH/CONNE"][:]
        a[1:5] = a[[4, 3, 2, 1]]                     # 最初の四角形の節点順を逆に
        f["VIZMESH/CONNE"][...] = a


expect_refused("接続違い (VIZMESH/CONNE の節点順)", shift_dR(), _viz, want={"connectivity", "logical_ij"})


def _roll_j(p):
    with h5py.File(p["dst_h5"], "r+") as f:
        c = f["MESH/COORD"][:].reshape(NI, NJ, 3)
        f["MESH/COORD"][...] = np.roll(c, -1, axis=1).ravel()   # 各列の節点番号に 1 層上の座標 (j ずれ)


expect_refused("j ずれ (各列の座標が 1 層ずれ)", None, _roll_j, want={"logical_ij", "displacement"})


def _bnodes_shift(p):
    with h5py.File(p["dst_h5"], "r+") as f:
        a = f["BCONDS/3/vizBfaceNodes"][:]
        f["BCONDS/3/vizBfaceNodes"][...] = a - 1      # 壁の節点集合を j − 1 に
    with h5py.File(p["src_mesh"], "r+") as f:
        a = f["BCONDS/3/vizBfaceNodes"][:]
        f["BCONDS/3/vizBfaceNodes"][...] = a - 1      # 両側でずらしても論理辺 (j = nj − 1) と合わない


expect_refused("j ずれ (壁の節点集合が 1 層内側)", shift_dR(), _bnodes_shift, want={"boundaries"})


def _kind(p):
    with h5py.File(p["dst_h5"], "r+") as f:
        f["BCONDS/3"].attrs["bcondKind"] = "wall_noSlip"


expect_refused("境界種別違い (slip → wall_noSlip)", shift_dR(), _kind, want={"boundaries"})
f_ = expect_refused("移動超過 (壁を外向きに 1.5 µm)", shift_dR(cols=(2,), layers=(5,), d=+1.5e-6), want={"displacement"})
check(f"移動超過: ほかの検査は通る ({sorted(f_ or [])})", f_ == {"displacement"})
f_ = expect_refused("要素のねじれ (壁 1 層だけ内向きに 0.5 µm、0.3 µm の層を越える; 面積の符号は変わらない)",
                    shift_dR(cols=(2,), layers=(5,)), want={"no_inversion"})
check(f"ねじれ: 移動量 (0.5 µm) は上限内 ({sorted(f_ or [])})", f_ is not None and "displacement" not in f_)
f_ = expect_refused("要素の反転 (壁際 2 層を内向きに 0.9 µm、0.7 µm の層を越える; 面積の符号が変わる)",
                    shift_dR(cols=(2, 3), layers=(4, 5), d=-0.9e-6), want={"no_inversion"})


def _drop_roe(p):
    with h5py.File(p["res"], "r+") as f:
        del f["VALUE/roe"]


expect_refused("量欠落 (SRC に roe が無い)", shift_dR(), _drop_roe, want={"species_energy"})


def _extra(p):
    with h5py.File(p["res"], "r+") as f:
        f.create_dataset("VALUE/roXi", data=np.zeros(NI * NJ, np.float32))


expect_refused("SRC の保存量が DST に無い (roXi を黙って落とさない)", shift_dR(), _extra, want={"species_energy"})


def _nan(p):
    with h5py.File(p["res"], "r+") as f:
        a = f["VALUE/roUx"][:]
        a[5] = np.nan
        f["VALUE/roUx"][...] = a


expect_refused("SRC に非有限値", shift_dR(), _nan, want={"species_energy"})


def _phys(p):
    (p["dst"] / "solverConfig.yaml").write_text(CFG.replace("gamma: 1.4", "gamma: 1.3"))


expect_refused("physProp 違い (エネルギー基準)", shift_dR(), _phys, want={"species_energy"})


CFG_BLOCK = """# 同じ設定をブロック形式・コメント・引用符付きキー・余分な空白で書く (解釈は同じ)
mesh:
  discretization : node
  "isAxisymmetric": 1
  axisCentroidShift: 1
  meshFileName: nozzle.h5
  valueFileName: nozzle.h5
base: &pp {thermalMethod: 0, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: 1004.5, gamma: 1.4}
physProp:
  <<: *pp
turbulence: {model: none}
"""


def _cfg_block(p):
    (p["dst"] / "solverConfig.yaml").write_text(CFG_BLOCK.replace("base: &pp", "_base: &pp"))
    (p["src"] / "solverConfig.yaml").write_text(CFG.replace("physProp:", "_base: {thermalMethod: 0, viscMethod: 0, visc: 0.0, thermCond: 0.0, cp: 1004.5, gamma: 1.4}\nphysProp:"))


with tempfile.TemporaryDirectory() as t:
    p = setup(Path(t), shift_dR(), _cfg_block)
    try:
        out = run(p)
        check("YAML: ブロック形式・コメント・引用符付きキー・merge キー (<<) でも同じ解釈なら通る", out["VERDICT"] == "OK")
    except M.MapRefused as e:
        check(f"YAML: ブロック形式・merge キーでも同じ解釈なら通る (拒否された: {e.failures})", False)


def _dup(p):
    (p["dst"] / "solverConfig.yaml").write_text(CFG.replace("gamma: 1.4}", "gamma: 1.4, gamma: 1.4}"))


expect_refused("YAML: 重複キー (同じ値でも) を拒否", shift_dR(), _dup, want={"coordinate_system_units"})


def _dup_block(p):
    (p["dst"] / "solverConfig.yaml").write_text(CFG + "turbulence: {model: none}\n")


expect_refused("YAML: トップレベルの重複キーを拒否", shift_dR(), _dup_block, want={"coordinate_system_units"})


def _sci(p):
    (p["dst"] / "solverConfig.yaml").write_text(CFG.replace("cp: 1004.5", "cp: 1.0045e3"))


f_ = expect_refused("YAML: 指数表記で値が同じでも型が違えば (文字列 1.0045e3 と 1004.5) 拒否", shift_dR(), _sci, want={"species_energy"})


def _nocfg(p):
    (p["dst"] / "solverConfig.yaml").unlink()


expect_refused("solverConfig.yaml が無い", shift_dR(), _nocfg, want={"coordinate_system_units"})


def _scale(p):
    (p["dst"] / "prepare_info.json").write_text(json.dumps({"scale_m": SCALE * 1000, "mesh": {"ni": NI, "nj": NJ}}))


expect_refused("単位違い (scale_m が mm)", shift_dR(), _scale, want={"coordinate_system_units"})


def _axisym(p):
    (p["dst"] / "solverConfig.yaml").write_text(CFG.replace("isAxisymmetric: 1", "isAxisymmetric: 0"))


expect_refused("座標系違い (軸対称の有無)", shift_dR(), _axisym, want={"coordinate_system_units"})


def _res_coord(p):
    with h5py.File(p["res"], "r+") as f:
        c = f["MESH/COORD"][:]
        c[4] += np.float32(1e-6)
        f["MESH/COORD"][...] = c


expect_refused("保存場の座標が IC 格子と違う", shift_dR(), _res_coord, want={"coordinate_system_units"})
with tempfile.TemporaryDirectory() as t:
    p = setup(Path(t), shift_dR(cols=(2,), layers=(5,), d=+1.5e-6))
    rc = M.main([str(p["res"]), str(p["dst_h5"]), "--src-mesh", str(p["src_mesh"]), "--dst-run", str(p["dst"])])
    check(f"CLI: 不成立なら exit 2 ({rc})", rc == 2)
    rc = M.main([str(p["res"]), str(p["dst_h5"]), "--src-mesh", str(p["src_mesh"]), "--dst-run", str(p["dst"]), "--max-disp-m", "2e-6"])
    check(f"CLI: --max-disp-m で上限を変えると通る ({rc})", rc == 0)
    try:
        M.apply_map(p["res"], p["dst_h5"], p["src_mesh"], p["src"], p["dst"], max_disp_m=float("nan"))
        check("max_disp_m = NaN を拒否 (通ってしまった)", False)
    except ValueError:
        check("max_disp_m = NaN を拒否", True)
    try:
        M.apply_map(p["res"], p["dst_h5"], p["src_mesh"], p["src"], p["dst"], mode="cubic")
        check("未知のモードを拒否 (通ってしまった)", False)
    except ValueError:
        check("未知のモードを拒否", True)

# --- 補助関数 ------------------------------------------------------------------------------------------------------------
try:
    M.structured_shape(np.arange(7), 36)
    check("structured_shape: 5 の倍数でない長さを拒否 (通ってしまった)", False)
except ValueError:
    check("structured_shape: 5 の倍数でない長さを拒否", True)
check("structured_shape: 正しい並びから (ni, nj) を復元", M.structured_shape(viz_conne(NI, NJ), NI * NJ) == (NI, NJ))
with tempfile.TemporaryDirectory() as t:
    p = setup(Path(t), shift_dR())
    shutil.copy(p["dst_h5"], Path(t) / "copy.h5")
    a, b = M.mesh_digest(p["dst_h5"]), M.mesh_digest(Path(t) / "copy.h5")
    with h5py.File(Path(t) / "copy.h5", "r+") as f:
        f["VALUE/ro"][...] = 0.0
    c = M.mesh_digest(Path(t) / "copy.h5")
    d = M.mesh_digest(p["src_mesh"])
    check("mesh_digest: 同じ格子で一致・/VALUE に依存しない・座標が違えば不一致", a == b == c and a != d)
print("FAIL 件数:", FAIL)
sys.exit(1 if FAIL else 0)
