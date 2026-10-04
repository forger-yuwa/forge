#!/usr/bin/env python3
"""#4b (plan gradient-scalar-lsq-unification §5.1 #4b、判定 A/B/B′ は測る前に固定) と #4a の「ダンプ有効でも res が変わらない」確認。

tgv (三重周期 SST + ξ、`s1_nointerference.py` の `s1_prep/tgv`) を gg 3 本・lsq 3 本、`FORGE_DUMP_PREGATHER` を付けて 1 step。
  (i)  NS 18 配列 (d{ro,Ux,Uy,Uz,P,T}d[xyz]) の gather 前の局所配列 (tag init = res_0 の元、loop1 = res_1 の元) が 6 本すべてでビット同一か。
  (ii) gather 後 (res_0 / res_1) の NS 18 配列で 6 本の間に不一致がある節点について、各 run の値が member 部分和の float32 順列和
       (root の値から member を任意順に足す = `periodicGather1ToRoot_d` の atomicAdd がとりうる値) の集合に含まれるか。
       参考に、不一致の無い節点も含めて 2 member 以上の全 group で同じ検査をする。
  判定: A = (i) 同一 かつ (ii) 全点で含まれる / B = (i) 不一致 / B′ = (i) 同一だが (ii) が外れる。
#4a 確認: 同じ入力でダンプ無しの gg 3 本 (S1 の s1_tgv_gg_{a,b,c}) とダンプ有りの gg 3 本の res_0・res_1 全配列を noise_rule で比べる。

使い方 (AWS):
  python3 pregather_check.py run <scratch> --prep <s1_prep/tgv> --bin BIN
  python3 pregather_check.py compare <scratch> --nodump-runs DIR_a DIR_b DIR_c [--out PREGATHER_4b.txt]
"""
import argparse
import os
import shutil
import sys

import h5py
import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gharness as G  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = [("gg", r) for r in "abc"] + [("lsq", r) for r in "abc"]
NS6 = ("ro", "Ux", "Uy", "Uz", "P", "T")


def rdir(scr, kind, r):
    return os.path.join(scr, f"pg_tgv_{kind}_{r}")


def cmd_run(scr, prep, fbin):
    for kind, r in RUNS:
        d = rdir(scr, kind, r)
        if os.path.exists(os.path.join(d, "res_1.h5")):
            print("exists", d); continue
        if os.path.exists(d):
            shutil.rmtree(d)
        shutil.copytree(prep, d)
        if kind == "lsq":
            p = os.path.join(d, "solverConfig.yaml")
            c = yaml.safe_load(open(p)); c["mesh"]["scalarGradient"] = "lsq"
            yaml.safe_dump(c, open(p, "w"), sort_keys=False, default_flow_style=None)
        os.environ["FORGE_DUMP_PREGATHER"] = os.path.join(d, "pregather")
        G.run_forge(d, expect_merge=True, forge_bin=fbin)
        os.environ.pop("FORGE_DUMP_PREGATHER", None)
        print("ran", d, sorted(f for f in os.listdir(d) if f.startswith("pregather")))


def load(path):
    with h5py.File(path, "r") as f:
        return {k: np.array(f["VALUE"][k]) for k in f["VALUE"].keys()}


def cmd_compare(scr, nodump, outp):
    out = ["# #4b pre-gather ダンプによる tgv S1(2) の判定 (plan gradient-scalar-lsq-unification §5.1 #4b) と #4a のダンプ非干渉",
           f"harness revision: {G.git_rev()}"]
    P = out.append
    d0 = rdir(scr, "gg", "a")
    msh = G.Mesh(os.path.join(d0, "mesh.h5"), os.path.join(d0, "bcondConfig.yaml"))
    P(f"run: {scr}/pg_tgv_{{gg,lsq}}_{{a,b,c}}")
    P(G.provenance(d0).rstrip())
    P(f"CV {msh.nCells}、周期 group {len(msh.groups)} (member 数の内訳 " +
      ", ".join(f"{n}:{int((msh.nmember == n).sum())}" for n in (1, 2, 4, 8)) + ")")
    for kind, r in RUNS:
        fs = sorted(f for f in os.listdir(rdir(scr, kind, r)) if f.startswith("pregather") and not f.endswith(".names"))
        P(f"- {kind}_{r}: ダンプ {fs}")
    # (i)
    P("\n## (i) NS 18 配列の gather 前 (tag init → res_0、loop1 → res_1) が 6 本でビット同一か")
    P("| tag | 変数 | 6 本のうち a (gg_a) と異なる run の数 | 不一致節点数 (最大) |")
    P("| --- | --- | --- | --- |")
    ok_i = True
    pre = {}
    for tag in ("init", "loop1"):
        pre[tag] = {(k, r): G.read_pregather(os.path.join(rdir(scr, k, r), "pregather." + tag)) for k, r in RUNS}
        ref = pre[tag][("gg", "a")]
        for v in NS6:
            nd, nn = 0, 0
            for key in RUNS:
                a = pre[tag][key][v]
                m = np.any(a.view(np.uint32) != ref[v].view(np.uint32), axis=1)
                if m.any():
                    nd += 1; nn = max(nn, int(m.sum()))
            ok_i &= nd == 0
            P(f"| {tag} | d{v} | {nd} | {nn} |")
    P(f"(i): {'ビット同一' if ok_i else '不一致あり'}")
    # (ii)
    P("\n## (ii) gather 後の値が member 部分和の float32 順列和の集合に含まれるか")
    P("gather 前は gg_a のダンプを使う ((i) が同一ならどれでも同じ)。対象 = 6 本の間で gather 後の値に不一致がある節点、"
      "参考 = 2 member 以上の全 group の全節点")
    P("| res | 変数 | 6 本で不一致の節点 (その member 数) | 不一致節点で集合外の (run, 節点, 成分) | 全 group で集合外 | 判定 |")
    P("| --- | --- | --- | --- | --- | --- |")
    ok_ii = True
    for st, tag in ((0, "init"), (1, "loop1")):
        R = {key: load(os.path.join(rdir(scr, *key), f"res_{st}.h5")) for key in RUNS}
        for v in NS6:
            p0 = pre[tag][("gg", "a")][v]
            sets = G.gather_perm_sums(msh, p0)
            post = {key: np.stack([R[key][f"d{v}d{c}"] for c in "xyz"], 1)[:msh.nCells] for key in RUNS}
            ref = post[("gg", "a")]
            mism = np.zeros(msh.nCells, bool)
            for key in RUNS:
                mism |= np.any(post[key].view(np.uint32) != ref.view(np.uint32), axis=1)
            bad_m = bad_all = 0
            for key in RUNS:
                u = post[key].view(np.uint32)
                for rt, s3 in sets.items():
                    for c in range(3):
                        for node in msh.groups[rt]:
                            inset = u[node, c].item() in s3[c]
                            if not inset:
                                bad_all += 1
                                if mism[node]:
                                    bad_m += 1
            ok = bad_m == 0 and bad_all == 0
            ok_ii &= ok
            mm = sorted(set(msh.nmember[mism].tolist()))
            P(f"| {st} | d{v} | {int(mism.sum())} ({mm}) | {bad_m} | {bad_all} | {'ok' if ok else 'NG'} |")
    verdict = "A" if (ok_i and ok_ii) else ("B" if not ok_i else "B′")
    P(f"\nVERDICT #4b: {verdict} ((i) {'同一' if ok_i else '不一致'}、(ii) {'全点で含まれる' if ok_ii else '外れあり'})")
    # #4a: ダンプ有効でも res が変わらない
    if nodump:
        out.extend(dump_noninterference(nodump, [rdir(scr, "gg", r) for r in "abc"], "tgv"))
    txt = "\n".join(out) + "\n"
    open(outp, "w").write(txt)
    print(txt)


def dump_noninterference(nodump, withdump, label):
    out = []
    P = out.append
    if True:
        P("\n## #4a ダンプ非干渉: ダンプ無し gg 3 本 vs ダンプ有り gg 3 本 (res_0・res_1 の全配列、noise_rule)")
        P("無し: " + ", ".join(nodump))
        okd = True
        rows = []
        nbit = ntot = 0
        for st in (0, 1):
            A = [load(os.path.join(d, f"res_{st}.h5")) for d in nodump]
            B = [load(os.path.join(d, f"res_{st}.h5")) for d in withdump]
            for k in sorted(A[0]):
                if A[0][k].dtype.kind != "f":
                    continue
                if any(k not in x for x in B):
                    rows.append(f"| res_{st} | {k} | ダンプ有りに無い | FAIL |"); okd = False; continue
                ntot += 1
                r = G.noise_rule([x[k] for x in A], [x[k] for x in B])
                if all(t[0] == 0 for t in r["aa"] + r["bb"] + r["ab"]):
                    nbit += 1; continue
                okd &= r["verdict"] == "PASS"
                rows.append(f"| res_{st} | {k} | 無し同士 n≤{max(t[0] for t in r['aa'])}、有り同士 n≤{max(t[0] for t in r['bb'])}、"
                            f"無し–有り n≤{max(t[0] for t in r['ab'])} / 最大差 {max(t[1] for t in r['ab']):.3g} | {r['verdict']} |")
        P(f"全 {ntot} 配列中 {nbit} 配列は全 15 対でビット一致")
        P("| res | 配列 | 不一致 | 判定 |")
        P("| --- | --- | --- | --- |")
        out.extend(rows)
        P(f"VERDICT #4a ダンプ非干渉 ({label}): {'PASS' if okd else 'FAIL'}")
    return out


def cmd_dumpcheck(scr, prep, fbin, nodump, label, outp):
    """任意ケースで gg をダンプ有りで 3 本回し、ダンプ無し 3 本と res_0・res_1 の全配列を noise_rule で比べる。"""
    runs = []
    for r in "abc":
        d = os.path.join(scr, f"dc_{label}_{r}")
        runs.append(d)
        if os.path.exists(os.path.join(d, "res_1.h5")):
            continue
        if os.path.exists(d):
            shutil.rmtree(d)
        shutil.copytree(prep, d)
        os.environ["FORGE_DUMP_PREGATHER"] = os.path.join(d, "pregather")
        G.run_forge(d, expect_merge=False, forge_bin=fbin)
        os.environ.pop("FORGE_DUMP_PREGATHER", None)
    out = [f"# #4a ダンプ非干渉 ({label})", f"harness revision: {G.git_rev()}", G.provenance(runs[0]).rstrip(),
           "ダンプ有り: " + ", ".join(runs) + " (ダンプ " + str(sorted(f for f in os.listdir(runs[0]) if f.startswith("pregather"))) + ")"]
    out.extend(dump_noninterference(nodump, runs, label))
    txt = "\n".join(out) + "\n"
    open(outp, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("run", "compare", "dumpcheck"))
    ap.add_argument("--label", default="case48")
    ap.add_argument("scratch")
    ap.add_argument("--prep")
    ap.add_argument("--bin")
    ap.add_argument("--nodump-runs", nargs="*", default=[])
    ap.add_argument("--out", default=os.path.join(HERE, "PREGATHER_4b.txt"))
    a = ap.parse_args()
    if a.cmd == "run":
        cmd_run(a.scratch, a.prep, a.bin)
    elif a.cmd == "dumpcheck":
        cmd_dumpcheck(a.scratch, a.prep, a.bin, a.nodump_runs, a.label, a.out)
    else:
        cmd_compare(a.scratch, a.nodump_runs, a.out)
