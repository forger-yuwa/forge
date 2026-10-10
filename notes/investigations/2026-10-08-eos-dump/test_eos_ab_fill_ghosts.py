#!/usr/bin/env python3
"""eos_ab_dir.py の make-replay-fill-ghosts と fill-ghosts の再生判定の単体試験 (合成の格子・ダンプ、GPU 不要)。

plan §5.1 #2 の追加の受入れ試験 (5) の入力を作る道具の試験:
  (a) ghost と owner の対応: BCONDS を名前のバイト列の昇順 ("1" "10" "2") に回った通し番号で ghost を組む
      (数値順・iBPlanes の値から組むと別の境界を指す配置にしてある)。owner は iCells と PLANES/STRUCT の 2 経路で一致。
  (b) 正常: 非有限の index 集合 = 対象 ghost の集合 → 再生ファイルを書く。対象 ghost の読む入力一式 = owner (ビット)、
      それ以外の配列・要素はビット不変、読む入力は全節点で有限、記録の kind = fill_ghosts。
  (c)-(k) 止まるべき入力: 件数は同じで集合が違う・実節点の非有限・iCells と STRUCT の不一致・ghost 数の不一致・
      同じ種類の境界が 2 つ (--pid で解決)・--pid の種類違い・forge_run.log の順の不一致・再生のダンプを BASE に渡す・
      cell が 2 つの境界面・既存の出力。
  (l) 再生の判定 fill_replay_verdict: PASS / FAIL (旧新の差) / 試験不成立 (v2 INVALID・腕内の非再現・/pre 不一致・
      post の非有限・新版の eos 行の欠落・組の不足)。
  (m) check_fill_dump: /pre の一致・post の非有限の数。
  (n) コマンド行: fill-ghosts の再生に --node を付ける → エラー、make-replay の再生に --node が無い → エラー。

    python3 test_eos_ab_fill_ghosts.py [--tools solver_density_cuda/tools]
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import eos_ab_dir as E   # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--tools", default=os.path.join(HERE, "..", "..", "..", "solver_density_cuda", "tools"))
TOOLS = os.path.abspath(ap.parse_args().tools)
ced = E.load_ced(TOOLS)

ok_all = True
NC = 12                                  # 実節点
# 境界条件 (名前, 種類, owner 列)。面の index は変換器と同じく physID の数値順に並べる: "1" → "2" → "10"
BC = [("1", "inlet_uniformVelocity", [0, 1]), ("2", "slip", [2, 3, 4]), ("10", "farfield", [5, 6, 7, 8])]
NNORMAL = 20
NB = sum(len(o) for _, _, o in BC)       # 9
NALL = NC + NB
# forge の順 ("1", "10", "2") の通し番号: "1" → 12,13 / "10" → 14..17 / "2" → 18..20
GHOST_FF = np.arange(14, 18)
OWNER_FF = np.array([5, 6, 7, 8])
READS = ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1", "T"]


def chk(name, cond, detail=""):
    global ok_all
    ok_all = ok_all and bool(cond)
    print("  [%s] %-78s %s" % ("OK " if cond else "NG ", name, detail))


def write_mesh(path, bc=BC, mutate=None, nbplanes=None):
    """合成の格子: 内部面 NNORMAL 枚 (cell 2 つ) + 境界面 (cell 1 つ = owner)。BCONDS は数値順の面 index。"""
    strct, ip0, bcd = [], NNORMAL, {}
    for i in range(NNORMAL):
        strct += [2, i, i + 1, 2, i % NC, (i + 1) % NC]
    for nm, kind, own in bc:
        ips = list(range(ip0, ip0 + len(own)))
        for p, o in zip(ips, own):
            strct += [3, p, p + 1, p + 2, 1, o]
        bcd[nm] = {"kind": kind, "iPlanes": np.array(ips, np.int32), "iBPlanes": np.array(ips, np.int32) - NNORMAL,
                   "iCells": np.array(own, np.int32)}
        ip0 += len(own)
    if mutate:
        mutate(bcd, strct)
    with h5py.File(path, "w") as f:
        g = f.create_group("MESH")
        g.attrs["nCells"] = np.int32(NC)
        g.attrs["nNodes"] = np.int32(NC)
        g.attrs["nPlanes"] = np.int32(ip0)
        g.attrs["nNormalPlanes"] = np.int32(NNORMAL)
        g.attrs["nBPlanes"] = np.int32(nbplanes if nbplanes is not None else ip0 - NNORMAL)
        f.create_dataset("PLANES/STRUCT", data=np.array(strct, np.int32))
        for nm, d in bcd.items():
            gb = f.create_group(f"BCONDS/{nm}")
            gb.attrs["bcondKind"] = d["kind"]
            for k in ("iPlanes", "iBPlanes", "iCells"):
                gb.create_dataset(k, data=d[k])


def base_pre(rng):
    pre = {}
    for nm in READS + ["P", "Ht", "sonic", "gamma", "cp", "Rmix", "Ux", "Uy", "Uz", "k", "omega", "Y0", "Y1", "dUxdx", "volume"]:
        pre[nm] = (rng.random(NALL, dtype=np.float32) + np.float32(0.5)).astype(np.float32)
    for nm in ("roe", "T"):             # g3 と同じ形: farfield の ghost の roe・T が NaN
        pre[nm][GHOST_FF] = np.nan
    for nm in ("gamma", "cp"):          # 書き込み専用の配列は全 ghost が NaN (判定外)
        pre[nm][NC:] = np.nan
    return pre


def write_base(path, pre, meta=None, post=True):
    m = {"step": 1, "nCells": NC, "nCells_all": NALL, "sizeof_flow_float": 4, "thermalMethod": 2, "nSpeciesRegistered": 2,
         "nSpecies": 2, "thermoFloat": 1, "condensation": 0, "isImplicit": 1, "unsteady": 0, "dualTime": 0,
         "discretization": "node", "solverConfig_yaml": "a: 1\n", "db_n_species": 2}
    m.update(meta or {})
    with h5py.File(path, "w") as f:
        for k, v in m.items():
            f.attrs[k] = v
        f.create_dataset("db/species_thermo", data=np.frombuffer(b"\x01\x02\x03\x04", dtype=np.uint8))
        g = f.create_group("pre")
        g.attrs["null_arrays"] = "foo;bar"
        g.attrs["n_arrays"] = np.int64(len(pre))
        for k, v in pre.items():
            g.create_dataset(k, data=v)
        if post:
            for k, v in pre.items():
                f.create_dataset(f"post/{k}", data=np.nan_to_num(v, nan=1.0))


def make_input(d, mesh_mutate=None, nbplanes=None, bc=BC):
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "solverConfig.yaml"), "w").write('mesh: {meshFileName: "mesh.h5", valueFileName: "mesh.h5"}\n')
    write_mesh(os.path.join(d, "mesh.h5"), bc=bc, mutate=mesh_mutate, nbplanes=nbplanes)
    return d


def args(tmp, name, inp, base, **kw):
    a = SimpleNamespace(base=base, input=inp, tools=TOOLS, out=os.path.join(tmp, f"{name}.h5"), bcond_kind="farfield",
                        pid=None, forge_log=None)
    for k, v in kw.items():
        setattr(a, k, v)
    return a


def expect_stop(name, fn, want, a=None):
    try:
        fn()
    except SystemExit as e:
        msg = str(e)
        chk(name, want in msg and (a is None or not os.path.exists(a.out)), msg[:110])
        return msg
    chk(name, False, "止まらなかった")
    return ""


def main():
    tmp = tempfile.mkdtemp(prefix="fill_ghosts_tests_")
    rng = np.random.default_rng(0)
    pre = base_pre(rng)
    inp = make_input(os.path.join(tmp, "in"))
    base = os.path.join(tmp, "base.h5")
    write_base(base, pre)

    print("--- (a) ghost と owner の対応")
    gm = E.ghost_owner_map(os.path.join(inp, "mesh.h5"))
    chk("(a) 回る順 = 名前のバイト列の昇順 (1, 10, 2)", [r["name"] for r in gm["order"]] == ["1", "10", "2"], str([r["name"] for r in gm["order"]]))
    chk("(a) farfield の ghost = nCells + len(BCONDS/1) から 4 個 (14..17)", np.array_equal(gm["ghost"], GHOST_FF), str(gm["ghost"]))
    with h5py.File(os.path.join(inp, "mesh.h5"), "r") as f:
        naive_val = NC + int(f["BCONDS/10/iBPlanes"][0])                                     # iBPlanes の値から組む
    naive_num = NC + sum(len(o) for nm, _, o in BC if int(nm) < 10)                          # 数値順の通し番号
    naive = set(range(naive_val, naive_val + 4))
    chk("(a) iBPlanes の値・数値順から組むと 17..20 (うち 3 個は BCONDS/2 の ghost 18..20) を指す配置になっている",
        naive_val == naive_num == 17 and naive != set(gm["ghost"].tolist()) and len(naive & {18, 19, 20}) == 3,
        f"値から {naive_val}・数値順 {naive_num}")
    chk("(a) owner = iCells (= PLANES/STRUCT の cell)", np.array_equal(gm["owner"], OWNER_FF) and gm["owner_sources_agree"], str(gm["owner"]))
    chk("(a) bcond_order は数値順ではない", E.bcond_order(["2", "10", "1", "9"]) == ["1", "10", "2", "9"], "")

    print("--- (b) 正常")
    a = args(tmp, "ok", inp, base)
    E.make_replay_fill_ghosts(a)
    rec = json.load(open(a.out + ".json"))
    with h5py.File(a.out, "r") as o, h5py.File(base, "r") as b:
        names_same = sorted(o["pre"]) == sorted(b["pre"])
        ghost_eq = all(np.array_equal(o["pre"][k][()][GHOST_FF].view(np.uint32), b["pre"][k][()][OWNER_FF].view(np.uint32)) for k in READS)
        mask = np.ones(NALL, bool)
        mask[GHOST_FF] = False
        outside = all(np.array_equal(o["pre"][k][()][mask].view(np.uint32), b["pre"][k][()][mask].view(np.uint32)) for k in READS)
        others = all(np.array_equal(o["pre"][k][()].view(np.uint32), b["pre"][k][()].view(np.uint32)) for k in b["pre"] if k not in READS)
        finite = all(np.isfinite(o["pre"][k][()]).all() for k in READS)
        gamma_nan_kept = int((~np.isfinite(o["pre"]["gamma"][()])).sum()) == NB
        attrs_ok = o["pre"].attrs["null_arrays"] == "foo;bar" and int(o.attrs["nCells_all"]) == NALL and bytes(o["db/species_thermo"][()]) == b"\x01\x02\x03\x04"
        no_post = "post" not in o
    chk("(b) /pre の名前の集合・属性・物性 DB が BASE と同じ、/post は持たない", names_same and attrs_ok and no_post, "")
    chk("(b) 対象 ghost の読む入力一式 = owner (ビット)", ghost_eq, "")
    chk("(b) 読む配列の対象 ghost の外はビット不変", outside, "")
    chk("(b) 読む集合以外の配列はビット不変 (書き込み専用 gamma の ghost の NaN もそのまま)", others and gamma_nan_kept, "")
    chk("(b) 充填後の読む入力は全節点で有限", finite, "")
    chk("(b) 記録: kind fill_ghosts・pid 10・ghost [14, 18)・全検査 True・充填の数", rec["kind"] == "fill_ghosts" and rec["pid"] == "10"
        and rec["ghost_first"] == 14 and rec["ghost_end"] == 18 and all(rec["checks"].values()) and rec["changed_elements"]["roe"] == 4,
        str(rec["changed_elements"]))
    chk("(b) replay_kind() = fill_ghosts", E.replay_kind(a.out) == "fill_ghosts", "")

    print("--- (c)-(k) 止まるべき入力")
    p2 = {k: v.copy() for k, v in pre.items()}
    p2["roe"][17] = np.float32(1.0)
    p2["T"][17] = np.float32(1.0)
    p2["roe"][18] = np.nan            # 件数は 4 のまま、1 個を BCONDS/2 の ghost へ
    b2 = os.path.join(tmp, "base_c.h5")
    write_base(b2, p2)
    a = args(tmp, "c", inp, b2)
    expect_stop("(c) 件数は同じで集合が違う → 試験不成立 (書かない)", lambda: E.make_replay_fill_ghosts(a), "試験不成立", a)
    inv = json.load(open(a.out + ".invalid.json"))
    chk("(c) 不成立の記録に内訳 (対象外 1 個 = BCONDS/2 の ghost、対象のうち有限 1 個)", inv["nonfinite_extra"] == 1 and inv["target_ghost_finite"] == 1
        and inv["nonfinite_extra_where"].get("ghost BCONDS/2 (slip)") == 1, str(inv["nonfinite_extra_where"]))

    p3 = {k: v.copy() for k, v in pre.items()}
    p3["roUx"][3] = np.inf
    b3 = os.path.join(tmp, "base_d.h5")
    write_base(b3, p3)
    a = args(tmp, "d", inp, b3)
    expect_stop("(d) 実節点に非有限 → 試験不成立", lambda: E.make_replay_fill_ghosts(a), "試験不成立", a)
    p3b = {k: v.copy() for k, v in pre.items()}
    p3b["roOmega"][OWNER_FF[1]] = np.nan   # owner の非有限は (d) と同じく実節点の非有限なので集合の照合で止まる
    b3b = os.path.join(tmp, "base_d2.h5")
    write_base(b3b, p3b)
    a = args(tmp, "d2", inp, b3b)
    expect_stop("(d) owner の非有限 → 試験不成立 (集合の照合で止まる)", lambda: E.make_replay_fill_ghosts(a), "試験不成立", a)

    def bad_icells(bcd, strct):
        bcd["10"]["iCells"] = bcd["10"]["iCells"].copy()
        bcd["10"]["iCells"][1] = 9
    inp_e = make_input(os.path.join(tmp, "in_e"), mesh_mutate=bad_icells)
    a = args(tmp, "e", inp_e, base)
    expect_stop("(e) iCells と PLANES/STRUCT の owner が食い違う → 止まる", lambda: E.make_replay_fill_ghosts(a), "食い違う", a)

    inp_f = make_input(os.path.join(tmp, "in_f"), nbplanes=NB + 1)
    a = args(tmp, "f", inp_f, base)
    expect_stop("(f) len(iBPlanes) の和 ≠ nBPlanes → 止まる", lambda: E.make_replay_fill_ghosts(a), "nBPlanes", a)

    bc2 = [("1", "inlet_uniformVelocity", [0, 1]), ("2", "farfield", [2, 3, 4]), ("10", "farfield", [5, 6, 7, 8])]
    inp_g = make_input(os.path.join(tmp, "in_g"), bc=bc2)
    a = args(tmp, "g", inp_g, base)
    expect_stop("(g) farfield が 2 つで --pid 無し → 止まる", lambda: E.make_replay_fill_ghosts(a), "1 つに決まらない", a)
    gm2 = E.ghost_owner_map(os.path.join(inp_g, "mesh.h5"), pid="10")
    chk("(g) --pid 10 で選べる (ghost 14..17)", np.array_equal(gm2["ghost"], GHOST_FF), str(gm2["ghost"]))

    a = args(tmp, "h", inp, base, pid="2")
    expect_stop("(h) --pid の境界の種類が --bcond-kind と違う → 止まる", lambda: E.make_replay_fill_ghosts(a), "違う", a)

    def log(path, names):
        kinds = {nm: k for nm, k, _ in BC}
        ip = {"1": (20, 21), "2": (22, 24), "10": (25, 28)}
        with open(path, "w") as fh:
            fh.write(f"Number of Ghost Cells: {NB}\nNumber of All   Cells: {NALL}\n")
            for nm in names:
                fh.write(f"in mesh.cpp  physID={nm}\nin mesh.cpp  bcondKind={kinds[nm]}\n               ip min={ip[nm][0]}, ip max={ip[nm][1]}\n")
    lg_ok, lg_bad = os.path.join(tmp, "ok.log"), os.path.join(tmp, "bad.log")
    log(lg_ok, ["1", "10", "2"])
    log(lg_bad, ["1", "2", "10"])
    c_ok = E.check_forge_log_order(lg_ok, gm)
    chk("(i) forge_run.log の順 (1, 10, 2)・面の範囲・ghost 数が一致 → ok", c_ok["ok"], str(c_ok["mismatch"]))
    a = args(tmp, "i", inp, base, forge_log=lg_bad)
    expect_stop("(i) ログが数値順 (1, 2, 10) → 止まる", lambda: E.make_replay_fill_ghosts(a), "forge_run.log", a)

    b4 = os.path.join(tmp, "base_j.h5")
    write_base(b4, pre, meta={"replay": 1})
    a = args(tmp, "j", inp, b4)
    expect_stop("(j) BASE が再生のダンプ → 止まる", lambda: E.make_replay_fill_ghosts(a), "再生のダンプ", a)

    def two_cells(bcd, strct):
        # 最初の farfield の面 (index 25) の cell を 2 つにする: [3, p, p+1, p+2, 1, o] → [3, p, p+1, p+2, 2, o, 0]
        pos = 6 * NNORMAL + 6 * (2 + 3)
        assert strct[pos + 4] == 1
        strct[pos + 4:pos + 6] = [2, strct[pos + 5], 0]
    inp_k = make_input(os.path.join(tmp, "in_k"), mesh_mutate=two_cells)
    a = args(tmp, "k", inp_k, base)
    expect_stop("(k) cell が 2 つの境界面 → 止まる", lambda: E.make_replay_fill_ghosts(a), "cell が 1 つでない", a)
    a = args(tmp, "ok", inp, base)
    try:
        E.make_replay_fill_ghosts(a)
        chk("(k) 既存の出力は上書きしない", False, "止まらなかった")
    except SystemExit as e:
        chk("(k) 既存の出力は上書きしない", "上書きしない" in str(e), str(e)[:80])

    print("--- (l) 再生の判定 fill_replay_verdict")

    def S_of(mut=None):
        cmp = {k: {"verdict_v2": "IDENTICAL", "different": [], "invalid_v2": [], "eos_args_post_diff_bytes": 0, "eos_args_bytes": 100,
                   "eos_args_pre_diff_bytes": 0}
               for k in ("old_r1_vs_new_r1", "old_r1_vs_old_r2", "old_r2_vs_new_r1", "new_r1_vs_new_r2", "old_r1_vs_new_r2")}
        eos_row = {"kind": "eos", "step": "1", "inner": "0", "q_index": "0", "nT_real": "0", "ids_T": "", "nRho_real": "0", "nP_real": "0",
                   "n_near_real": "0", "overflow": "0"}
        runs = {n: {"floor_rows": ([dict(eos_row)] if n.startswith("new") else [])} for n in ("old_r1", "new_r1", "old_r2", "new_r2")}
        rc = {n: {"pre_equals_replay": True, "pre_mismatch_arrays": [], "reads_pre_nonfinite_total": 0, "reads_pre_nonfinite": {},
                  "writes_post_nonfinite_total": 0, "writes_post_nonfinite": {}, "writes_post_nonfinite_target": {}, "missing_arrays": [],
                  "unsupported": [], "target_post": {}} for n in runs}
        S = {"compare": cmp, "runs": runs, "replay_checks": rc}
        if mut:
            mut(S)
        return E.fill_replay_verdict(SimpleNamespace(step=1), S)

    chk("(l) 5 組とも v2 IDENTICAL・/pre 一致・新版の eos 行あり → PASS", S_of()[0] == "PASS", "")

    def m_fail(S):
        S["compare"]["old_r1_vs_new_r1"].update(verdict_v2="DIFFERENT", different=["T"])
    chk("(l) 旧版と新版の間の v2 DIFFERENT → FAIL", S_of(m_fail)[0] == "FAIL", "")

    def m_arm(S):
        S["compare"]["new_r1_vs_new_r2"].update(verdict_v2="DIFFERENT", different=["T"])
    chk("(l) 同じ版どうしの v2 DIFFERENT (腕内の非再現) → 試験不成立", S_of(m_arm)[0] == "試験不成立", "")

    def m_inv(S):
        S["compare"]["old_r2_vs_new_r1"].update(verdict_v2="INVALID", invalid_v2=["EOS の書き込み先 T の post に非有限"])
    chk("(l) v2 INVALID → 試験不成立", S_of(m_inv)[0] == "試験不成立", "")

    def m_pre(S):
        S["replay_checks"]["old_r2"].update(pre_equals_replay=False, pre_mismatch_arrays=["roe"])
    chk("(l) ダンプの /pre が再生ファイルと違う → 試験不成立", S_of(m_pre)[0] == "試験不成立", "")

    def m_post(S):
        S["replay_checks"]["new_r1"].update(writes_post_nonfinite_total=2, writes_post_nonfinite={"T": 2})
    chk("(l) 書き込み先の post に非有限 → 試験不成立", S_of(m_post)[0] == "試験不成立", "")

    def m_row(S):
        S["runs"]["new_r2"]["floor_rows"] = []
    chk("(l) 新版の step の eos 行が無い (カウンタ有効の腕でない) → 試験不成立", S_of(m_row)[0] == "試験不成立", "")

    def m_n(S):
        del S["compare"]["old_r1_vs_new_r2"]
    chk("(l) 比較の組が 5 でない → 試験不成立", S_of(m_n)[0] == "試験不成立", "")

    print("--- (m) check_fill_dump")
    rp = os.path.join(tmp, "ok.h5")
    with h5py.File(rp, "r") as r:
        rpre = {k: r["pre"][k][()] for k in r["pre"]}
    dump = os.path.join(tmp, "dump_ok.h5")
    post = {k: np.nan_to_num(v, nan=1.0) for k, v in rpre.items()}
    with h5py.File(dump, "w") as f:
        for k, v in h5py.File(base, "r").attrs.items():
            f.attrs[k] = v
        f.attrs["replay"] = 1
        for k, v in rpre.items():
            f.create_dataset(f"pre/{k}", data=v)
        for k, v in post.items():
            f.create_dataset(f"post/{k}", data=v)
    c = E.check_fill_dump(dump, rp, json.load(open(rp + ".json")), ced)
    chk("(m) /pre = 再生ファイル・読む入力の pre と書き込み先の post に非有限なし", c["pre_equals_replay"] and c["reads_pre_nonfinite_total"] == 0
        and c["writes_post_nonfinite_total"] == 0, "")
    with h5py.File(dump, "r+") as f:
        x = f["pre/roe"][()]
        x[0] = np.float32(x[0] * 2)
        f["pre/roe"][...] = x
        y = f["post/T"][()]
        y[15] = np.nan
        f["post/T"][...] = y
    c = E.check_fill_dump(dump, rp, json.load(open(rp + ".json")), ced)
    chk("(m) /pre の 1 要素の差・対象 ghost の post T の NaN を検出", (not c["pre_equals_replay"]) and c["pre_mismatch_arrays"] == ["roe"]
        and c["writes_post_nonfinite_target"]["T"] == 1, str(c["pre_mismatch_arrays"]))

    print("--- (n) コマンド行")
    me = os.path.join(HERE, "eos_ab_dir.py")
    r = subprocess.run([sys.executable, me, "run", "--input", inp, "--work", os.path.join(tmp, "w"), "--old", "/bin/true", "--new", "/bin/true",
                        "--tools", TOOLS, "--replay", rp, "--node", "3"], capture_output=True, text=True)
    chk("(n) fill-ghosts の再生に --node → エラー (終了コード 2)", r.returncode == 2 and "--node" in r.stderr, r.stderr.strip().splitlines()[-1][:90])
    rp_old = os.path.join(tmp, "floor_kind.h5")
    with h5py.File(rp_old, "w") as f:
        f.create_group("replay_info").attrs["record_json"] = json.dumps({"node": 3})
    r = subprocess.run([sys.executable, me, "run", "--input", inp, "--work", os.path.join(tmp, "w"), "--old", "/bin/true", "--new", "/bin/true",
                        "--tools", TOOLS, "--replay", rp_old], capture_output=True, text=True)
    chk("(n) make-replay の再生に --node・--expect-floor が無い → エラー (終了コード 2)", r.returncode == 2 and "--node" in r.stderr,
        r.stderr.strip().splitlines()[-1][:90])

    if ok_all:
        shutil.rmtree(tmp)
    else:
        print(f"作業ディレクトリ (失敗したので残す) {tmp}")
    print("ALL OK" if ok_all else "FAILED")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
