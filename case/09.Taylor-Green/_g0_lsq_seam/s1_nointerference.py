#!/usr/bin/env python3
"""S1 非干渉 (plan gradient-scalar-lsq-unification §6 S1、判定規則は 2026-09-26 の訂正版 = gharness.noise_rule)。

(1) `scalarGradient` 省略 (= gg) の新バイナリが実装前バイナリと全配列で一致する (旧経路の保存、roK/roOmega の初期ミラー追加を含む)。
    旧 3 本・新 gg 3 本の同一設定反復。res_0 と res_1 の全 VALUE 配列を noise_rule で判定する。
(2) `scalarGradient: lsq` で NS の勾配・リミタ配列が gg と同じ step でビット一致 (1 step。nStepInner は線形 sweep 数なので
    1 step の勾配評価は 1 回 = res_1 の勾配配列が初回の評価)。新 gg 3 本・新 lsq 3 本で同じ規則。
    対象 = d{ro,Ux,Uy,Uz,P,T}d[xyz] と limiter_{ro,Ux,Uy,Uz,P,T} (スカラーの勾配・リミタは作用素の差で変わるので除く)。

ケース (入力は同一のものを 9 run へ複製):
  case48 / case16 / axi : 前提 plan R3 の準備物 `r3_prep/<case>` (`r3_prepare.py` 製。case48 = SST 壁あり、case16 = 2 成分 SFR 1、
                          axi = 軸対称 2 成分)
  tgv                   : TGV 32³ 三重周期に SST + ξ (SFR 1)、k・ω・ξ を sin 場で焼いたもの (周期 group の初期ミラーを通す)

使い方 (AWS):
  python3 s1_nointerference.py prepare <scratch> --r3prep DIR --tgv-src DIR/Taylor-Green.h5
  python3 s1_nointerference.py run <scratch> --old BIN --new BIN
  python3 s1_nointerference.py compare <scratch> [--out S1.txt]
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
CASES = ("case48", "case16", "axi", "tgv")
REPS = ("a", "b", "c")
KINDS = ("old", "gg", "lsq")     # old = 実装前バイナリ (既定 gg)、gg / lsq = 新バイナリ
NS_ARR = [f"d{v}d{c}" for v in ("ro", "Ux", "Uy", "Uz", "P", "T") for c in "xyz"] + \
         [f"limiter_{v}" for v in ("ro", "Ux", "Uy", "Uz", "P", "T")]


def rd(scr, case, kind, rep):
    return os.path.join(scr, f"s1_{case}_{kind}_{rep}")


def prep_tgv(d, tgv_src):
    os.makedirs(d, exist_ok=False)
    cfg = G.base_solver_cfg(sst=True, tracer=True)
    h5 = G.make_run(d + "_tmp", tgv_src, cfg, os.path.join(HERE, "bcondConfig.yaml"))
    msh = G.Mesh(h5, os.path.join(d + "_tmp", "bcondConfig.yaml"))
    w, _, _ = G.wrapped_coords(msh)
    x = msh.xyz.astype(np.float64); ext = x.max(0) - x.min(0)
    k = G.periodic_sin(w, ext, 2.0, 0.5, (0.3, 1.1, -0.4))
    om = G.periodic_sin(w, ext, 100.0, 20.0, (1.7, -0.6, 0.9))
    xi = G.periodic_sin(w, ext, 0.5, 0.2, (-0.8, 0.4, 2.1))
    ux = G.periodic_sin(w, ext, 10.0, 2.0, (0.1, 0.2, 0.3))
    G.bake(h5, G.prim_state_fields(msh, 1.0, ux, 0.0 * ux, 0.0 * ux, 101325.0, k=k, om=om, xi=xi))
    for f in os.listdir(d + "_tmp"):
        shutil.move(os.path.join(d + "_tmp", f), d)
    os.rmdir(d + "_tmp")


def cmd_prepare(scr, r3prep, tgv_src):
    pdir = os.path.join(scr, "s1_prep")
    os.makedirs(pdir, exist_ok=False)
    for c in ("case48", "case16", "axi"):
        d = os.path.join(pdir, c); os.makedirs(d)
        for f in os.listdir(os.path.join(r3prep, c)):
            if not f.endswith((".msh", ".geo", ".xmf", ".log", "quality.txt")):
                shutil.copy(os.path.join(r3prep, c, f), d)
    prep_tgv(os.path.join(pdir, "tgv"), tgv_src)
    for c in CASES:
        for kind in KINDS:
            for r in REPS:
                d = rd(scr, c, kind, r)
                shutil.copytree(os.path.join(pdir, c), d)
                if kind == "lsq":
                    p = os.path.join(d, "solverConfig.yaml")
                    cfg = yaml.safe_load(open(p))
                    cfg["mesh"]["scalarGradient"] = "lsq"
                    yaml.safe_dump(cfg, open(p, "w"), sort_keys=False, default_flow_style=None)
    print("prepared", scr)


def cmd_run(scr, old, new):
    for c in CASES:
        for kind in KINDS:
            for r in REPS:
                d = rd(scr, c, kind, r)
                if os.path.exists(os.path.join(d, "res_1.h5")):
                    print("exists", d); continue
                G.run_forge(d, expect_merge=(c == "tgv"), forge_bin=(old if kind == "old" else new))
                print("ran", d)


def load(path):
    with h5py.File(path, "r") as f:
        return {k: np.array(f["VALUE"][k]) for k in f["VALUE"].keys()}


def cmd_compare(scr, outp):
    out = ["# S1 非干渉 (plan gradient-scalar-lsq-unification §6 S1、2026-09-26 訂正規則)",
           f"harness revision: {G.git_rev()}",
           "規則: 旧同士 (3 対) がビット一致の配列は旧新 (9 対) もビット一致。それ以外は旧新の最大差 ≤ ノイズ対 (旧旧 3 + 新新 3) の最大差 × 2 "
           "かつ 旧新の不一致数の最大 ≤ ノイズ対の不一致数の最大 × 2。不一致位置が 3 member 以上の group に限られるかは記録 (判定外、tgv のみ)"]
    P = out.append
    allok = True
    for c in CASES:
        P(f"\n## {c}")
        for kind in KINDS:
            prov = [l.strip() for l in open(os.path.join(rd(scr, c, kind, "a"), "RUN_PROVENANCE.txt"))
                    if l.startswith(("forge_sha256", "'scalarGradient'"))]
            P(f"- {kind}: {rd(scr, c, kind, 'a')} (.._b, .._c) {prov}")
        nan = []
        for kind in KINDS:
            for r in REPS:
                for st in (0, 1):
                    for k, a in load(os.path.join(rd(scr, c, kind, r), f"res_{st}.h5")).items():
                        if a.dtype.kind == "f" and not np.isfinite(a).all():
                            nan.append(f"{kind}_{r}:res_{st}:{k}")
        P(f"NaN/Inf: {nan if nan else 'なし'}")
        allok &= not nan
        msh = None
        if c == "tgv":
            d0 = rd(scr, c, "gg", "a")
            msh = G.Mesh(os.path.join(d0, "mesh.h5"), os.path.join(d0, "bcondConfig.yaml"))
        for title, A_kind, B_kind, names in (("(1) gg: 旧 vs 新", "old", "gg", None),
                                             ("(2) NS 勾配・リミタ: 新 gg vs 新 lsq", "gg", "lsq", NS_ARR)):
            P(f"\n### {title}")
            ok_all = True
            nbit = 0; ntot = 0
            rows = []
            for st in (0, 1):
                As = [load(os.path.join(rd(scr, c, A_kind, r), f"res_{st}.h5")) for r in REPS]
                Bs = [load(os.path.join(rd(scr, c, B_kind, r), f"res_{st}.h5")) for r in REPS]
                keys = sorted(set(As[0]) | set(Bs[0])) if names is None else [k for k in names if k in As[0]]
                for k in keys:
                    if any(k not in x for x in As + Bs):
                        rows.append(f"| res_{st} | {k} | 配列が片方に無い | - | - | FAIL |"); ok_all = False; continue
                    if As[0][k].dtype.kind != "f":
                        continue
                    ntot += 1
                    r = G.noise_rule([x[k] for x in As], [x[k] for x in Bs])
                    ok = r["verdict"] == "PASS"
                    ok_all &= ok
                    if all(t[0] == 0 for t in r["aa"] + r["bb"] + r["ab"]):
                        nbit += 1; continue
                    loc = "-"
                    if msh is not None and r["mask"] is not None and r["mask"].shape[0] >= msh.nCells:
                        mm = r["mask"][:msh.nCells]
                        if mm.any():
                            loc = "≥3 member のみ" if (msh.nmember[mm] >= 3).all() else \
                                f"member {sorted(set(msh.nmember[mm].tolist()))}"
                    f = lambda L: f"n {max(t[0] for t in L)} / max {max(t[1] for t in L):.3g}"
                    rows.append(f"| res_{st} | {k} | {f(r['aa'])} | {f(r['bb'])} | {f(r['ab'])} | {r['verdict']} | {loc} |")
            P(f"全 {ntot} 配列中 {nbit} 配列は全 15 対でビット一致 (表に出さない)")
            P("| res | 配列 | A 同士 (不一致数の最大 / 最大差) | B 同士 | A–B (9 対) | 判定 | 不一致位置 |")
            P("| --- | --- | --- | --- | --- | --- | --- |")
            out.extend(rows)
            P(f"VERDICT S1 {c} {title}: {'PASS' if ok_all else 'FAIL'}")
            allok &= ok_all
    P(f"\nVERDICT S1 (全ケース): {'PASS' if allok else 'FAIL'}")
    txt = "\n".join(out) + "\n"
    open(outp, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("prepare", "run", "compare"))
    ap.add_argument("scratch")
    ap.add_argument("--r3prep")
    ap.add_argument("--tgv-src")
    ap.add_argument("--old")
    ap.add_argument("--new")
    ap.add_argument("--out", default=os.path.join(HERE, "S1.txt"))
    a = ap.parse_args()
    if a.cmd == "prepare":
        cmd_prepare(a.scratch, a.r3prep, a.tgv_src)
    elif a.cmd == "run":
        cmd_run(a.scratch, a.old, a.new)
    else:
        cmd_compare(a.scratch, a.out)
