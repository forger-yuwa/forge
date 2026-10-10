#!/usr/bin/env python3
"""g3 → g4 の初期場の内外ラベルの 0 step A/B の道具 (`case/46.sern_design/diag/twin_ic_check.py`) の試験
(plan tooling-sern-te-wake-grid §5.1 #4、codex diagnose 2026-10-10 g4-initial-field)。forge は回さない。

合格 (測る前に固定):
  (a) 合成の小さな双子の格子 (mesh_sern3d.generate_sern_mesh3d に仮の設計を渡し、変換器を通さず node の h5 を書く。
      板厚 0 = カウルの上下の面も座標一致):
      (a1) 同じ格子を g3・g4 に渡す: 双子 = カウル・側壁・角の 3 分類とも > 0、正解なし 0、構造上の対 (dup1・dup2) がすべて
           座標一致・座標一致の組がすべて構造上の対。A は全組が同一 donor で誤割当 = 組数 (座標が同じ = 照会も同じ)。
           B (paste_region_ic3d) は誤割当 0、B2 (h5 に prepare と同じく書いた初期場) も一致 → `B_OK (A n 件誤り)`、n = 組数、終了コード 0
      (a2) g3 粗・g4 細 (station・層の数が違う): B_OK、A の誤割当 > 0
      (a3) g4 の節点番号を並べ替えた h5 (RENUMBER_PERM 付き、変換器の `mesh.renumber: rcm` と同じ形): paste の index 分類が
           崩れる → B_NG、終了コード 1 (構造上の対は RENUMBER_PERM で写して再現できる)
      (a4) g4 の sern.h5 の VALUE で双子の 1 節点だけ反対側の状態に書き換える: B は 0 件だが B2 が不一致 → B_NG
      (a5) g4 の sern.h5 の VALUE をどちらの状態でもない値にする: B2 は照合なし、判定は B だけで B_OK
      (a6) 正解が決まらない: 双子の 1 節点に反対側のタグの面を足す (両側のタグ) → UNDECIDABLE、終了コード 2
      (a7) 双子が 0 組 (双子の片方の座標をずらす) → UNDECIDABLE
      (a8) g4 にタグ (sidewall_out) の境界面の接続が無い → UNDECIDABLE
      (a9) prepare_info.json に格子署名があるのに h5 から再計算できない → UNDECIDABLE。discretization が node でない
           (paste が座標の分類に落ちる) → UNDECIDABLE
      (a10) 座標一致の組は数値で作る: −0.0 と 0.0 の双子も 1 組
      (a11) g3 の内部節点のラベル (格子の辺に沿った距離): カウルの上の第 1 層 = 排気側・下 = 外気側、側壁の内側の第 1 層 = 排気側・外側 = 外気側
            (自由端・後縁・上流端から離れた所)
  (b) runner_sern3d.prepare で変換した粗い格子 (生産 YAML problem_3d_prod_3op_wallres_lswx08_tewake_B10.yaml の構造で
      z・ノズル区間・層を減らした g3 相当と、それより細かい g4 相当、m10_on。板厚 0.005 H):
      双子が側壁・角とも > 0、正解なし 0 (両格子)、座標一致の組がすべて構造上の対、A の誤割当 > 0 で双子の全組が同一 donor、
      B の誤割当 0・B2 一致 → `TWIN_IC VERDICT: B_OK (A n 件誤り)`、終了コード 0。
      A の donor は `interp_field.py` の CLI を実際に回した結果とビット一致 (SRC の VALUE/ro に節点番号を入れて転送し、DST の値を読む)。
      prepare_info.json の格子署名を書き換えた g4 (別の格子の情報) → UNDECIDABLE

  /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_sern_twin_ic_tests.py
"""
import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOOLS = ROOT / "solver_density_cuda" / "tools"
CHECK = ROOT / "case" / "46.sern_design" / "diag" / "twin_ic_check.py"
sys.path.insert(0, str(HERE.parent))
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402
from forge_design.meshing.mesh_sern3d import PHYS_SERN3D, SernMesh3DParams, generate_sern_mesh3d  # noqa: E402

FAIL = 0
tmp = Path(tempfile.mkdtemp(prefix="twin_ic_tests_"))
H = 0.1


def check(name, cond, detail="", only_on_fail=False):
    global FAIL
    show = detail and (not only_on_fail or not cond)
    print(("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if show else ""))
    if not cond:
        FAIL += 1


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, str(path))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


TC = load_module(CHECK, "twin_ic_check_tests")


def cli(*args):
    return subprocess.run([sys.executable, str(CHECK), *map(str, args)], capture_output=True, text=True)


def verdict_of(r):
    lines = [l for l in r.stdout.splitlines() if l.startswith("TWIN_IC VERDICT:")]
    return lines[-1].split(":", 1)[1].strip() if lines else None


# ------------------------------------------------------------------------------------------------ (a) 合成の格子
ST = {"gas_model": "cpg",
      "exhaust": {"ro": 2.0, "u": 300.0, "P": 2.0e5, "k": 10.0, "omega": 1.0e4},
      "ext": {"ro": 1.0, "u": 600.0, "P": 1.0e5, "k": 1.0, "omega": 1.0e3}}
DESIGN = SimpleNamespace(cowl_xy=np.array([[0.0, 0.0], [0.6, -0.03], [1.2, -0.1]]),
                         ramp_xy=np.array([[0.0, 1.0], [1.0, 0.7], [2.0, 0.4]]), L_ramp=2.0)
XDMF = lambda h: np.column_stack([np.full(len(h), 9, np.int64), h]).ravel().astype(np.int32)


def synth_prm(**kw):
    p = SernMesh3DParams(ni_up=3, ni_noz=9, ni_plume=8, nj_top=6, nj_bot=5, nz_in=4, nz_out=3, L_sw=0.8, cowl_thickness=0.0,
                         first_wall_frac=0.02, first_z_frac=0.05, t_base=0.0, scale=H)
    return replace(p, **kw)


def write_synth(run, prm, perm_seed=None, disc="node"):
    """合成の run: sern.h5 (node: CV 数 = 節点数、MESH/COORD float32・VIZMESH/CONNE・BCONDS/*/vizBface*・CELLS/centCoords・VALUE)、
    prepare と同じく paste_region_ic3d で初期場を書く。prepare_info.json (mesh・states・H_m・half_W_m・discretization)。"""
    run.mkdir(parents=True)
    coords, hexes, B, info, _ = generate_sern_mesh3d(DESIGN, prm)
    n = len(coords)
    perm = None
    if perm_seed is not None:                    # perm[新] = 旧 (変換器の RENUMBER_PERM と同じ向き)
        perm = np.random.default_rng(perm_seed).permutation(n)
        inv = np.empty_like(perm); inv[perm] = np.arange(n)
        coords = coords[perm]; hexes = inv[hexes]; B = {k: [tuple(inv[list(q)]) for q in v] for k, v in B.items()}
    with h5py.File(run / "sern.h5", "w") as f:
        f["MESH/COORD"] = coords.astype(np.float32).ravel()
        f["MESH"].attrs["nCells"] = np.int32(n); f["MESH"].attrs["nNodes"] = np.int32(n)
        if perm is not None:
            f["MESH/RENUMBER_PERM"] = perm.astype(np.int64)
        f["VIZMESH/CONNE"] = XDMF(hexes)
        f["CELLS/centCoords"] = coords.astype(np.float32).ravel()
        for nm, q in B.items():
            if len(q):
                q = np.asarray(q, np.int32)
                f[f"BCONDS/{PHYS_SERN3D[nm]}/vizBfaceSizes"] = np.full(len(q), 4, np.int32)
                f[f"BCONDS/{PHYS_SERN3D[nm]}/vizBfaceNodes"] = q.ravel()
        for k in ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "wall_dist"):
            f["VALUE/" + k] = np.zeros(n, np.float32)
    info_j = json.loads(json.dumps(info))
    R3.paste_region_ic3d(run / "sern.h5", None, H, 0.5 * prm.W * H, ST, 1.4, minfo=info_j)
    (run / "prepare_info.json").write_text(json.dumps({"mesh": info_j, "states": ST, "H_m": H, "half_W_m": 0.5 * prm.W * H,
                                                       "discretization": disc}))
    return info, coords, B


def run_tool(g3, g4, *extra):
    r = cli("--g3", g3, "--g4", g4, "--json", tmp / "last.json", *extra)
    rep = json.loads((tmp / "last.json").read_text()) if (tmp / "last.json").exists() else {}
    (tmp / "last.json").unlink(missing_ok=True)
    return r, rep


def synth_tests():
    print("--- (a) 合成の小さな双子の格子")
    a3, a4 = tmp / "s_g3", tmp / "s_g4"
    info3, _, _ = write_synth(a3, synth_prm())
    info4, c4, B4 = write_synth(a4, synth_prm(ni_noz=13, nj_top=8, nj_bot=6, nz_in=5, nz_out=4))
    # (a1) 同じ格子
    r, rep = run_tool(a3, a3)
    t = rep.get("g4_truth", {})
    np_ = t.get("pairs", 0)
    check("(a1) 同じ格子: 双子がカウル・側壁・角とも > 0、正解なし 0",
          all(t.get("by_category", {}).get(c, 0) > 0 for c in TC.CATS) and t.get("n_undecided") == 0, json.dumps(t, ensure_ascii=False)[:400], only_on_fail=True)
    s = rep.get("g4_struct", {})
    check("(a1) 構造上の対 (dup1・dup2) はすべて座標一致、座標一致の組はすべて構造上の対",
          s.get("dup1") == s.get("dup1_coincident") == info3["n_dup_cowl"] and s.get("dup2") == s.get("dup2_coincident") == info3["n_dup_side"]
          and s.get("coincident_not_structural") == 0 and np_ == info3["n_dup_cowl"] + info3["n_dup_side"], str(s))
    A = rep.get("A", {}).get("twins_total", {})
    check("(a1) A: 双子の全組が同一 donor、誤割当 = 組数 (各組の片側)", A.get("same_donor_pairs") == np_ and A.get("wrong") == np_, str(A))
    Bt = rep.get("B", {})
    check("(a1) B: 双子・壁節点とも誤割当 0、B2 も一致", Bt.get("双子 (計)", {}).get("wrong") == 0 and Bt.get("壁節点", {}).get("wrong") == 0
          and Bt.get("B2_available") and Bt.get("B2_mismatch") == 0, json.dumps(Bt, ensure_ascii=False)[:300], only_on_fail=True)
    check(f"(a1) VERDICT: B_OK (A {np_} 件誤り)、終了コード 0", r.returncode == 0 and verdict_of(r) == f"B_OK (A {np_} 件誤り)",
          r.stdout[-600:] + r.stderr[-600:], only_on_fail=True)
    # (a2) 粗 → 細
    r, rep = run_tool(a3, a4)
    A = rep.get("A", {}).get("twins_total", {})
    check("(a2) g3 粗・g4 細: B_OK、A の誤割当 > 0", r.returncode == 0 and verdict_of(r).startswith("B_OK") and A.get("wrong", 0) > 0,
          f"{verdict_of(r)} {A}")
    # (a3) 番号の並べ替え
    a4p = tmp / "s_g4_perm"
    write_synth(a4p, synth_prm(ni_noz=13, nj_top=8, nj_bot=6, nz_in=5, nz_out=4), perm_seed=7)
    r, rep = run_tool(a3, a4p)
    s = rep.get("g4_struct", {})
    check("(a3) RENUMBER_PERM 付き: paste の index 分類が崩れて B_NG、終了コード 1",
          r.returncode == 1 and verdict_of(r) == "B_NG" and rep.get("B", {}).get("双子 (計)", {}).get("wrong", 0) > 0, r.stdout[-500:], only_on_fail=True)
    check("(a3) 構造上の対は RENUMBER_PERM で写して全組が座標一致", s.get("dup1") == s.get("dup1_coincident") and s.get("dup2") == s.get("dup2_coincident")
          and s.get("coincident_not_structural") == 0, str(s))

    def mutated(name, fn):
        d = tmp / name
        shutil.copytree(a4, d)
        with h5py.File(d / "sern.h5", "r+") as f:
            fn(f, d)
        return d
    lab4, _, _ = TC.wall_labels(TC.load_grid(str(a4 / "sern.h5"), str(a4), PHYS_SERN3D))
    pairs4, _ = TC.coincident_groups(np.asarray(h5py.File(a4 / "sern.h5", "r")["MESH/COORD"][...]).reshape(-1, 3) + np.float32(0))
    tin = int(pairs4[0][0] if lab4[pairs4[0][0]] == TC.IN else pairs4[0][1])

    # (a4) B2 だけ不一致
    def flip(f, d):
        from forge_design.evaluate import runner_sern as R2
        ext = R2.region_ic_arrays(np.array([False]), ST, 1.4)      # 外気の状態 (排気側の双子に書く)
        for k, a in ext.items():
            v = f["VALUE/" + k][...]; v[tin] = a[0]; f["VALUE/" + k][...] = v
    r, rep = run_tool(a3, mutated("s_g4_b2flip", flip))
    Bt = rep.get("B", {})
    check("(a4) h5 の初期場で双子の 1 節点を外気の状態に: B は 0 件、B2 不一致 1 → B_NG",
          verdict_of(r) == "B_NG" and r.returncode == 1 and Bt.get("双子 (計)", {}).get("wrong") == 0 and Bt.get("B2_mismatch") == 1,
          json.dumps(Bt, ensure_ascii=False)[:300], only_on_fail=True)

    # (a5) B2 照合なし
    def scramble(f, d):
        f["VALUE/ro"][...] = np.full(f["VALUE/ro"].shape, 1.5, np.float32)
    r, rep = run_tool(a3, mutated("s_g4_b2none", scramble))
    check("(a5) h5 の初期場が読めない: B2 照合なし、判定は B だけで B_OK", r.returncode == 0 and verdict_of(r).startswith("B_OK")
          and rep.get("B", {}).get("B2_available") is False and "照合なし" in r.stdout, r.stdout[-400:], only_on_fail=True)

    # (a6) 正解が決まらない双子
    def both_tags(f, d):
        pid = PHYS_SERN3D["cowl_out"] if lab4[tin] == TC.IN else PHYS_SERN3D["cowl_in"]
        g = f[f"BCONDS/{pid}"]
        sz = np.concatenate([g["vizBfaceSizes"][...], [4]]); nd = np.concatenate([g["vizBfaceNodes"][...], [tin, tin, tin, tin]])
        del g["vizBfaceSizes"], g["vizBfaceNodes"]
        g["vizBfaceSizes"] = sz.astype(np.int32); g["vizBfaceNodes"] = nd.astype(np.int32)
    r, rep = run_tool(a3, mutated("s_g4_both", both_tags))
    check("(a6) 双子の 1 節点が両側のタグ → 正解なし 1 組、UNDECIDABLE、終了コード 2",
          r.returncode == 2 and verdict_of(r) == "UNDECIDABLE" and rep.get("g4_truth", {}).get("n_undecided") == 1, r.stdout[-500:], only_on_fail=True)

    # (a7) 双子 0 組
    def unpair(f, d):
        c = f["MESH/COORD"][...].reshape(-1, 3)
        for a, b in pairs4:
            c[b, 1] += np.float32(1e-3)
        f["MESH/COORD"][...] = c.ravel()
    r, rep = run_tool(a3, mutated("s_g4_nopair", unpair))
    check("(a7) 双子 0 組 → UNDECIDABLE", r.returncode == 2 and verdict_of(r) == "UNDECIDABLE" and "0 組" in r.stdout, r.stdout[-400:], only_on_fail=True)

    # (a8) タグの境界面の接続が無い
    def drop_tag(f, d):
        del f[f"BCONDS/{PHYS_SERN3D['sidewall_out']}"]
    r, _ = run_tool(a3, mutated("s_g4_notag", drop_tag))
    check("(a8) sidewall_out の境界面の接続が無い → UNDECIDABLE", r.returncode == 2 and verdict_of(r) == "UNDECIDABLE"
          and "sidewall_out" in r.stdout, r.stdout[-300:], only_on_fail=True)

    # (a9) 格子署名の不一致・node でない
    def bad_sig(f, d):
        j = json.loads((d / "prepare_info.json").read_text()); j["mesh_provenance"] = {"mesh_signature": "0" * 64}
        (d / "prepare_info.json").write_text(json.dumps(j))
    r, _ = run_tool(a3, mutated("s_g4_badsig", bad_sig))
    check("(a9) prepare_info に格子署名があるのに h5 から再計算できない (合成の h5 は PLANES が無い) → UNDECIDABLE",
          r.returncode == 2 and "格子署名を再計算できない" in r.stdout, r.stdout[-300:] + r.stderr[-300:], only_on_fail=True)

    def cell_disc(f, d):
        j = json.loads((d / "prepare_info.json").read_text()); j["discretization"] = "cell"
        (d / "prepare_info.json").write_text(json.dumps(j))
    r, _ = run_tool(a3, mutated("s_g4_cell", cell_disc))
    check("(a9) discretization cell (paste が座標の分類に落ちる) → UNDECIDABLE", r.returncode == 2 and "座標の分類" in r.stdout, r.stdout[-300:], only_on_fail=True)

    # (a10) −0.0 と 0.0
    c = np.array([[0.0, 0.0, 0.1], [0.0, -0.0, 0.1], [1.0, 2.0, 3.0]], np.float32) + np.float32(0.0)
    p, big = TC.coincident_groups(c)
    c2 = np.array([[0.0, 0.0, 0.1], [0.0, -0.0, 0.1]], np.float32)
    p2, _ = TC.coincident_groups(c2 + np.float32(0.0))
    check("(a10) −0.0 と 0.0 の双子も 1 組 (load_grid と同じく +0.0 で揃える)", p.tolist() == [[0, 1]] and not big and p2.tolist() == [[0, 1]])

    # (a11) g3 の内部節点のラベル
    g3 = TC.load_grid(str(a3 / "sern.h5"), str(a3), PHYS_SERN3D, need_hex=True)
    lab3, _, _ = TC.wall_labels(g3)
    side, _ = TC.side_labels(g3, lab3, 0.5 * H)
    m = info3
    NJ, nz, jm, ksw, isw = m["NJ"], m["nz"], m["jm"], m["k_sw"], m["i_sw"]
    base = lambda i, j, k: (i * NJ + j) * nz + k
    above = [base(i, jm + 1, k) for i in range(1, isw - 1) for k in range(ksw)]
    below = [base(i, jm - 1, k) for i in range(1, isw - 1) for k in range(ksw)]
    s_in = [base(i, j, ksw - 1) for i in range(1, isw - 1) for j in range(jm + 1, NJ - 1)]
    s_out = [base(i, j, ksw + 1) for i in range(1, isw - 1) for j in range(jm + 1, NJ - 1)]
    check("(a11) g3 の内部節点: カウルの上の第 1 層 = 排気側・下 = 外気側、側壁の内側の第 1 層 = 排気側・外側 = 外気側",
          np.all(side[above] == TC.IN) and np.all(side[below] == TC.OUT) and np.all(side[s_in] == TC.IN) and np.all(side[s_out] == TC.OUT),
          f"上 {np.unique(side[above])} 下 {np.unique(side[below])} 内 {np.unique(side[s_in])} 外 {np.unique(side[s_out])}")


# ------------------------------------------------------------------------------------------------ (b) 変換した粗い格子
def converted_tests():
    import yaml
    from forge_design.evaluate import runner_sern as R2
    from forge_design.probdef import load_problem
    print("--- (b) runner_sern3d.prepare で変換した粗い格子 (g3 相当・g4 相当)")
    prob = ROOT / "case" / "46.sern_design" / "problem_3d_prod_3op_wallres_lswx08_tewake_B10.yaml"
    p = load_problem(prob)
    cached = R2.design_from_problem(p)                   # 生産の設計 MOC (~10 s)。2 つの prepare で使い回す
    saved = R2.design_from_problem
    R2.design_from_problem = lambda p_, design=None: cached
    try:
        runs = {}
        common = {"first_wake_frac": 0.004, "first_z_frac": 0.02, "nj_vside": 5}
        for lab, upd in (("g3c", {"first_wall_frac": 2.56e-3, "nz_in": 4, "nz_out": 3, "ni_noz": 30}),
                         ("g4c", {"first_wall_frac": 1.6e-3, "nz_in": 5, "nz_out": 4, "ni_noz": 37, "nj_top": 81, "nj_bot": 61})):
            y = yaml.safe_load(open(prob))
            y["mesh3d"].update({**common, **upd})
            y["mesh"].update({"nj_ext_top": 9, "first_top_frac": 0.004})
            pth = tmp / f"prob_{lab}.yaml"
            pth.write_text(yaml.safe_dump(y, allow_unicode=True))
            R3.prepare(pth, tmp / f"run_{lab}", op="m10_on")
            runs[lab] = tmp / f"run_{lab}"
    finally:
        R2.design_from_problem = saved
    g3, g4 = runs["g3c"], runs["g4c"]
    n3 = json.loads((g3 / "prepare_info.json").read_text())["mesh"]["nodes"]
    n4 = json.loads((g4 / "prepare_info.json").read_text())["mesh"]["nodes"]
    print(f"    g3 相当 {n3} 節点・g4 相当 {n4} 節点")
    r, rep = run_tool(g3, g4)
    print("\n".join("    " + l for l in r.stdout.splitlines()))
    t3, t4 = rep.get("g3_truth", {}), rep.get("g4_truth", {})
    check("(b) 双子が側壁・角とも > 0 (両格子)", all(t.get("by_category", {}).get(c, 0) > 0 for t in (t3, t4) for c in TC.CATS[1:]),
          f"{t3.get('by_category')} {t4.get('by_category')}")
    check("(b) 正解の決まらない双子 0 (両格子)", t3.get("n_undecided") == 0 and t4.get("n_undecided") == 0)
    check("(b) 座標一致の組はすべて構造上の対 (両格子)", rep.get("g3_struct", {}).get("coincident_not_structural") == 0
          and rep.get("g4_struct", {}).get("coincident_not_structural") == 0)
    A = rep.get("A", {}).get("twins_total", {})
    check("(b) A: 誤割当 > 0、双子の全組が同一 donor", A.get("wrong", 0) > 0 and A.get("same_donor_pairs") == t4.get("pairs"), str(A))
    Bt = rep.get("B", {})
    check("(b) B: 双子・壁節点とも誤割当 0、B2 (prepare が書いた初期場) も一致",
          Bt.get("双子 (計)", {}).get("wrong") == 0 and Bt.get("壁節点", {}).get("wrong") == 0 and Bt.get("B2_available") and Bt.get("B2_mismatch") == 0)
    check("(b) VERDICT: B_OK (A n 件誤り)、終了コード 0", r.returncode == 0 and verdict_of(r) == f"B_OK (A {A.get('wrong')} 件誤り)", r.stderr[-500:], only_on_fail=True)

    # A の donor = interp_field.py の CLI の結果 (SRC の VALUE/ro に節点番号を入れて転送)
    d = tmp / "if_eq"; d.mkdir()
    shutil.copy(g3 / "sern.h5", d / "src.h5"); shutil.copy(g4 / "sern.h5", d / "dst.h5")
    with h5py.File(d / "src.h5", "r+") as f:
        n = f["VALUE/ro"].shape[0]
        f["VALUE/ro"][...] = np.arange(n, dtype=np.float32)       # 節点数 < 2^24 なので float32 で正確
    rr = subprocess.run([sys.executable, str(TOOLS / "interp_field.py"), str(d / "src.h5"), str(d / "dst.h5"), "--force-species"],
                        capture_output=True, text=True)
    with h5py.File(d / "dst.h5", "r") as f:
        donor_cli = np.asarray(f["VALUE/ro"][...]).astype(np.int64)
    R = TC.run_check(str(g3), str(g4))[2]
    ev, dn = R["_A_eval"], R["_A_donor"]
    check("(b) A の donor が interp_field.py の CLI の転送結果とビット一致 (評価節点すべて)",
          rr.returncode == 0 and np.array_equal(donor_cli[ev], dn), (rr.stdout + rr.stderr)[-600:], only_on_fail=True)

    # prepare_info.json の格子署名が g4 の h5 と違う (別の格子の情報) → UNDECIDABLE
    gx = tmp / "run_g4c_badsig"
    shutil.copytree(g4, gx)
    j = json.loads((gx / "prepare_info.json").read_text())
    check("(b) prepare が格子署名を prepare_info.json に書いている", bool((j.get("mesh_provenance") or {}).get("mesh_signature")))
    j["mesh_provenance"]["mesh_signature"] = "0" * 64
    (gx / "prepare_info.json").write_text(json.dumps(j))
    r, _ = run_tool(g3, gx)
    check("(b) prepare_info の格子署名が g4 の h5 と違う → UNDECIDABLE、終了コード 2",
          r.returncode == 2 and verdict_of(r) == "UNDECIDABLE" and "格子署名" in r.stdout and "違う" in r.stdout, r.stdout[-300:], only_on_fail=True)


if __name__ == "__main__":
    try:
        synth_tests()
        converted_tests()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
    sys.exit(1 if FAIL else 0)
