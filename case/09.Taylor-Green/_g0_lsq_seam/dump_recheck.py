#!/usr/bin/env python3
"""#4d (plan gradient-scalar-lsq-unification §5.1 #4d、§6「ダンプ非干渉」): 既存の 1 step 出力を E/P/N 区分で再判定する (事後規則)。

区分はコード上の依存関係で先に決める (観測値で決めない。2026-09-26 codex (diagnose)):
  E (厳密一致): res_0 の状態量 (初期値 + BC、atomicAdd を通らない)。NS 18 勾配 (+ divU) の周期 group に属さない節点
                (LSQ gather は atomic なし)。周期の無いケースのリミタ配列 (入力が E、limiter の atomicAdd は診断用で既定 off)。
  P (順列集合): NS 18 勾配の 2 member 以上の節点。基準 run の pre-gather 部分和から作る float32 順列和の集合に含まれること
                (2 member は唯一の和)。基準 run の pre-gather は全 run で同一であることも E として確認する。
  N (両側ノイズ): 上以外。ノイズ対 (無し同士 3 + 有り同士 3) の最大差・不一致数の 2 倍以内。片側条項は使わない。
                不一致数と最大差を 3 組 (無し同士・有り同士・無し–有り) とも出す。
判定: E/P と N を別々に出す (2026-09-26 codex (diagnose) 2)。E/P: A = E が全て厳密一致 かつ P が全点で集合内 / B = それ以外。
  N: 全配列で 2 倍規則を満たせば PASS。E/P の A をダンプ非干渉の合格と読まない。

使い方 (AWS):
  python3 dump_recheck.py tgv   --off RUN... --on RUN... [--pre-ref ON_RUN] --out DUMP_RECHECK_tgv.txt
  python3 dump_recheck.py case48 --off RUN... --on RUN...                 --out DUMP_RECHECK_case48.txt
  --off / --on は同じ入力・同じバイナリの 1 step run (ダンプ無し / 有り)。tgv は gg・lsq を混ぜて渡してよい
  (E・P は全 run を 1 本の基準と比べる。N は gg 同士・lsq 同士で別に組む)。
"""
import argparse
import itertools
import os
import re
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gharness as G  # noqa: E402

NS6 = ("ro", "Ux", "Uy", "Uz", "P", "T")
NS_GRAD = {f"d{v}d{c}" for v in NS6 for c in "xyz"}
# res_0 の状態量 (初期値 + applyBconds。勾配・リミタ・乱流粘性など gather の下流は含めない)
STATE_RE = re.compile(r"^(ro|roUx|roUy|roUz|roe|roK|roOmega|roXi|roY\d+|Ux|Uy|Uz|P|T|Xi|Y\d+|k|omega|h0|Ht|roGamma|roReth)$")
LIM_RE = re.compile(r"lim", re.IGNORECASE)


def load(path):
    with h5py.File(path, "r") as f:
        return {k: np.array(f["VALUE"][k]) for k in f["VALUE"].keys()}


def kind_of(run):
    """run の scalarGradient (起動エコー)。無ければ gg。"""
    for fn in ("forge_run.log", "log.txt", "forge.log"):
        p = os.path.join(run, fn)
        if os.path.exists(p):
            t = open(p, errors="replace").read()
            m = re.search(r"'scalarGradient' effective: (\w+)", t)
            if m:
                return m.group(1)
    return "gg"


def mesh_path(run):
    for fn in ("mesh.h5", "axi.h5"):
        if os.path.exists(os.path.join(run, fn)):
            return os.path.join(run, fn)
    hs = [f for f in os.listdir(run) if f.endswith(".h5") and not f.startswith("res_")]
    return os.path.join(run, sorted(hs)[0])


def pairs_stats(X, Y=None):
    if Y is None:
        return [G.bitcmp(a, b)[:2] for a, b in itertools.combinations(X, 2)]
    return [G.bitcmp(a, b)[:2] for a in X for b in Y]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("label")
    ap.add_argument("--off", nargs="+", required=True)
    ap.add_argument("--on", nargs="+", required=True)
    ap.add_argument("--pre-ref", default=None, help="pre-gather 部分和の基準 run (ダンプ有り)。既定 = --on の先頭")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    out = [f"# #4d ダンプ非干渉と S1 (2) の再判定 ({a.label}) — E/P/N 区分 (事後規則、plan §6「ダンプ非干渉」)",
           f"harness revision: {G.git_rev()}"]
    P = out.append
    ref_run = a.pre_ref or a.on[0]
    P(G.provenance(ref_run).rstrip())
    runs = [("off", r) for r in a.off] + [("on", r) for r in a.on]
    kinds = {r: kind_of(r) for _, r in runs}
    import hashlib
    def sha(path):
        h = hashlib.sha256()
        with open(path, "rb") as fp:
            for b in iter(lambda: fp.read(1 << 20), b""):
                h.update(b)
        return h.hexdigest()[:16]
    for side, r in runs:
        fsha = next((l.split(":", 1)[1].strip()[:16] for l in G.provenance(r).splitlines() if l.startswith("forge_sha256")), "?")
        P(f"- {side} {kinds[r]}: {r} (forge {fsha}、mesh {sha(mesh_path(r))})")
    msh = G.Mesh(mesh_path(ref_run), os.path.join(ref_run, "bcondConfig.yaml"))
    periodic = len(msh.groups) > 0
    P(f"CV {msh.nCells}、周期 group {len(msh.groups)}" + (" (member 数の内訳 " + ", ".join(
        f"{n}:{int((msh.nmember == n).sum())}" for n in sorted(set(msh.nmember.tolist()))) + ")" if periodic else ""))
    ingrp = msh.nmember >= 2 if periodic else np.zeros(msh.nCells, bool)
    okE = okP = okN = True

    # --- pre-gather の同一性 (ダンプ有りの全 run、E) ---
    pre = {}
    if periodic:
        P("\n## E0: pre-gather NS 18 配列がダンプ有りの全 run でビット同一か (P の部分和の前提)")
        P("| tag | 変数 | 基準と異なる run 数 | 不一致節点数 (最大) |")
        P("| --- | --- | --- | --- |")
        for tag in ("init", "loop1"):
            pre[tag] = {r: G.read_pregather(os.path.join(r, "pregather." + tag)) for r in a.on}
            for v in NS6:
                nd = nn = 0
                for r in a.on:
                    m = np.any(pre[tag][r][v].view(np.uint32) != pre[tag][ref_run][v].view(np.uint32), axis=1)
                    if m.any():
                        nd += 1; nn = max(nn, int(m.sum()))
                okE &= nd == 0
                P(f"| {tag} | d{v} | {nd} | {nn} |")

    for st, tag in ((0, "init"), (1, "loop1")):
        R = {r: load(os.path.join(r, f"res_{st}.h5")) for _, r in runs}
        names = sorted(k for k in R[ref_run] if R[ref_run][k].dtype.kind == "f")
        E_rows, P_rows, N_rows = [], [], []
        nbitN = 0
        for k in names:
            if any(k not in R[r] for _, r in runs):
                E_rows.append(f"| res_{st} | {k} | 一部の run に無い | FAIL |"); okE = False; continue
            base = R[ref_run][k][:msh.nCells]
            is_nsg = k in NS_GRAD or k == "divU"
            if (st == 0 and STATE_RE.match(k)) or (is_nsg) or (not periodic and LIM_RE.search(k)):
                # E (NS 勾配は周期 group 外の節点だけ)
                sel = ~ingrp if is_nsg else np.ones(msh.nCells, bool)
                nrun = nn = 0
                for _, r in runs:
                    x = R[r][k][:msh.nCells]
                    m = x.view(np.uint32) != base.view(np.uint32) if x.dtype == np.float32 else x != base
                    m = (np.any(m, axis=1) if m.ndim > 1 else m) & sel
                    if m.any():
                        nrun += 1; nn = max(nn, int(m.sum()))
                okE &= nrun == 0
                E_rows.append(f"| res_{st} | {k}{' (group 外)' if is_nsg and periodic else ''} | {nrun} | {nn} | {'ok' if nrun == 0 else 'NG'} |")
                if is_nsg and periodic and k != "divU":
                    # P: 2 member 以上の節点が基準 pre-gather の順列和の集合内か
                    v, c = k[1:-2], "xyz".index(k[-1])
                    sets = G.gather_perm_sums(msh, pre[tag][ref_run][v])
                    bad = 0; mism = 0
                    for _, r in runs:
                        u = R[r][k][:msh.nCells].view(np.uint32)
                        mism += int(np.count_nonzero((u != base.view(np.uint32)) & ingrp))
                        for rt, s3 in sets.items():
                            for node in msh.groups[rt]:
                                if u[node].item() not in s3[c]:
                                    bad += 1
                    okP &= bad == 0
                    P_rows.append(f"| res_{st} | {k} | {mism} | {bad} | {'ok' if bad == 0 else 'NG'} |")
                continue
            # N: gg 同士・lsq 同士で別に組む
            for kd in sorted(set(kinds.values())):
                A = [R[r][k] for s, r in runs if s == "off" and kinds[r] == kd]
                B = [R[r][k] for s, r in runs if s == "on" and kinds[r] == kd]
                if len(A) < 2 or len(B) < 2:
                    continue
                aa, bb, ab = pairs_stats(A), pairs_stats(B), pairs_stats(A, B)
                if all(t[0] == 0 for t in aa + bb + ab):
                    nbitN += 1; continue
                noise = aa + bb
                n_no, m_no = max(t[0] for t in noise), max(t[1] for t in noise)
                n_ab, m_ab = max(t[0] for t in ab), max(t[1] for t in ab)
                ok = (m_ab <= 2.0 * m_no) and (n_ab <= 2.0 * n_no)
                f = lambda L: f"n≤{max(t[0] for t in L)} / max {max(t[1] for t in L):.3g}"
                N_rows.append(f"| res_{st} | {k} | {kd} | {f(aa)} | {f(bb)} | {f(ab)} | {'PASS' if ok else 'FAIL'} |")
        P(f"\n## res_{st}: E (厳密一致、基準 = {os.path.basename(ref_run)})")
        P("| res | 配列 | 基準と異なる run 数 | 不一致節点数 (最大) | 判定 |")
        P("| --- | --- | --- | --- | --- |")
        out.extend(E_rows)
        if P_rows:
            P(f"\n## res_{st}: P (2 member 以上の節点、pre-gather {tag} の順列和の集合)")
            P("| res | 配列 | 基準と値が違う (run, 節点) 数 | 集合外の (run, 節点) 数 | 判定 |")
            P("| --- | --- | --- | --- | --- |")
            out.extend(P_rows)
        P(f"\n## res_{st}: N (両側ノイズ、記録。全対ビット一致 {nbitN} 件は省略)")
        P("| res | 配列 | kind | 無し同士 | 有り同士 | 無し–有り | 2 倍規則 |")
        P("| --- | --- | --- | --- | --- | --- | --- |")
        out.extend(N_rows)
        okN_st = all(not r.endswith("FAIL |") for r in N_rows)
        P(f"N 区分 (res_{st}): {'全て 2 倍以内' if okN_st else '2 倍を超える配列あり'}")
        okN &= okN_st
    verdict = "A" if (okE and okP) else "B"
    P(f"\nVERDICT E/P ({a.label}): {verdict} (E {'全て厳密一致' if okE else '不一致あり'}、P {'全点集合内' if okP else '集合外あり' if periodic else '該当なし'})")
    P(f"VERDICT N ({a.label}): {'PASS (全配列 2 倍以内)' if okN else 'FAIL (2 倍を超える配列あり)'}")
    txt = "\n".join(out) + "\n"
    open(a.out, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
