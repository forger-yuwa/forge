#!/usr/bin/env python3
"""#8a (plan boundary-node-periodic-gradient-fix §5.1 #8a、codex result-2 M1): スカラー GG の面寄与ダンプを 3 バイナリで比較する。

`FORGE_DUMP_SCALARGRAD=<path>` を付けると `species_gradient_d` が atomicAdd に渡す面寄与 (sx,sy,sz)·φ_f を、化学種・受動種
それぞれ最初の呼び出しだけ raw float32 [nVar][nPlanes][3] で `<path>.species` / `<path>.passive` に書く (除外した面は 0)。
1 面 1 スレッドの非 atomic 書き込みなので、場 (res_*.h5) と違い面レベルでビット比較できる。

バイナリ (FORGE_BIN で渡す。タグは run 名に入る):
  new     : f67fe877 (新 + ダンプ)。除外条件 excludePeriodic = periodicSeamMergeActive (軸対称では false)、面フラグ planePeriodic
  olddump : 1266aba1 + 同じダンプ (notes/sessions/boundary-node-periodic-gradient-fix-8a-olddump.patch)。除外条件は死に条件のまま
  wrong   : f67fe877 の 2 か所を periodicNodeActive に戻した誤り版 (…-8a-wrong.patch)。軸対称 × 周期でも半割面を除外する

ケース (各 1 step、同一入力を 3 run へ cp):
  axi   : 軸対称 × 並進周期 (x)。#6a の準備物 `r6a_prep` (61×41、MIXDRY/H2O + 受動トレーサ ξ、speciesFaceReconstruction 1) をそのまま使う
  trans : 並進周期 (三重周期 TGV 32³、`lsqseam_m1/Taylor-Green.h5`、bcond はハーネスの `bcondConfig.yaml`)。設定は axi と同じ
          (化学種 2 + ξ、SFR 1) で isAxisymmetric 0。状態は周期関数で焼く (ρ・u・Y1・ξ、e0 は入力の一様値)

合格 (plan §5.1 #8a、測る前に固定。並進周期の行は 2026-09-26 に親が測定前に訂正):
  axi   : olddump = new が全面ビット同一 (化学種・受動種とも)。
  trans : 周期半割面以外の全面で olddump = new がビット同一。周期半割面は new で 0、olddump で非 0 (面の 3 成分が全て 0 でない)。
  検出力: wrong は axi で周期半割面の寄与が 0 になり、new と差が出る。
周期半割面 = h5 の BCONDS/<physID>/iPlanes のうち bcondConfig の kind が periodic のもの (PLANES の index)。

使い方 (AWS):
  python3 r8a_facedump.py prepare <scratch>
  FORGE_BIN=... python3 r8a_facedump.py run <scratch> <axi|trans> <new|olddump|wrong>
  python3 r8a_facedump.py compare <scratch> [--out R8a_facedump.txt]
"""
import argparse
import os
import shutil
import subprocess
import sys

import h5py
import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
RUN_CASE = os.path.join(REPO, "solver_density_cuda", "tools", "run_case.sh")
CASES = ("axi", "trans")
BINS = ("new", "olddump", "wrong")
FILES = {"axi": ("axi.h5", "solverConfig.yaml", "bcondConfig.yaml", "species_db.yaml", "probe.yaml"),
         "trans": ("tgv.h5", "solverConfig.yaml", "bcondConfig.yaml", "species_db.yaml", "probe.yaml")}


def prep_dir(scratch, case):
    return os.path.join(scratch, f"r8a_prep_{case}")


def run_dir(scratch, case, b):
    return os.path.join(scratch, f"r8a_{case}_{b}")


def bake_trans(h5):
    """三重周期 [0,2π]³ の周期関数で ρ・u・Y1・ξ を焼く。e0 = 入力 (一様) の比内部エネルギーの最小値。"""
    with h5py.File(h5, "r+") as f:
        V = f["VALUE"]
        xyz = np.array(f["MESH/COORD"]).reshape(-1, 3).astype(np.float64)
        n = V["ro"].shape[0]
        if xyz.shape[0] != n:
            raise SystemExit(f"MESH/COORD と VALUE の長さが合わない ({xyz.shape}, {n})")
        x, y, z = xyz.T
        ro_in = np.array(V["ro"], dtype=np.float64)
        e0 = float((np.array(V["roe"], dtype=np.float64) / ro_in
                    - 0.5 * (np.array(V["roUx"], dtype=np.float64) / ro_in) ** 2).min())
        ro = ro_in[0] * (1.0 + 0.05 * np.cos(x) * np.sin(y))
        ux = 10.0 + 2.0 * np.sin(y) * np.cos(z)
        uy = 2.0 * np.sin(z) * np.cos(x)
        uz = 2.0 * np.sin(x) * np.cos(y)
        y1 = 0.02 + 0.01 * np.sin(x + 0.3) * np.cos(y) * np.sin(z + 0.5)
        xi = 0.5 + 0.3 * np.sin(x + 0.7) * np.cos(y + 0.2) * np.cos(z)
        out = {"ro": ro, "roUx": ro * ux, "roUy": ro * uy, "roUz": ro * uz,
               "roe": ro * (e0 + 0.5 * (ux ** 2 + uy ** 2 + uz ** 2)),
               "roY0": ro * (1.0 - y1), "roY1": ro * y1, "roXi": ro * xi}
        for k in ("roK", "roOmega"):
            if k in V:
                del V[k]
        dt = V["ro"].dtype
        for k, a in out.items():
            if k in V:
                V[k][...] = a.astype(V[k].dtype)
            else:
                V.create_dataset(k, data=a.astype(dt))


def cmd_prepare(scratch):
    # axi: #6a の準備物をそのまま使う
    src = os.path.join(scratch, "r6a_prep")
    d = prep_dir(scratch, "axi")
    os.makedirs(d, exist_ok=False)
    for fn in FILES["axi"]:
        shutil.copyfile(os.path.join(src, fn), os.path.join(d, fn))
    # trans: 同じ設定を非軸対称にし、TGV 32³ の三重周期メッシュへ
    d = prep_dir(scratch, "trans")
    os.makedirs(d, exist_ok=False)
    with open(os.path.join(src, "solverConfig.yaml")) as fp:
        cfg = yaml.safe_load(fp)
    cfg["mesh"] = {"discretization": "node", "isAxisymmetric": 0, "meshFileName": "tgv.h5", "valueFileName": "tgv.h5"}
    with open(os.path.join(d, "solverConfig.yaml"), "w") as fp:
        yaml.safe_dump(cfg, fp, sort_keys=False, default_flow_style=None)
    shutil.copyfile(os.path.join(scratch, "lsqseam_m1", "Taylor-Green.h5"), os.path.join(d, "tgv.h5"))
    shutil.copyfile(os.path.join(HERE, "bcondConfig.yaml"), os.path.join(d, "bcondConfig.yaml"))
    for fn in ("species_db.yaml", "probe.yaml"):
        shutil.copyfile(os.path.join(src, fn), os.path.join(d, fn))
    bake_trans(os.path.join(d, "tgv.h5"))
    for case in CASES:
        for b in BINS:
            rd = run_dir(scratch, case, b)
            os.makedirs(rd, exist_ok=False)
            for fn in FILES[case]:
                shutil.copyfile(os.path.join(prep_dir(scratch, case), fn), os.path.join(rd, fn))
    print("prepared", scratch)


def cmd_run(scratch, case, b):
    rd = run_dir(scratch, case, b)
    if os.path.exists(os.path.join(rd, "res_1.h5")):
        raise SystemExit(f"既に回した跡がある: {rd}")
    env = dict(os.environ, FORGE_CUDA_BLOCKSIZE="128", FORGE_CUDA_BLOCKSIZE_SMALL="128",
               FORGE_DUMP_SCALARGRAD=os.path.join(rd, "sgdump"))
    r = subprocess.run(["bash", RUN_CASE, rd], env=env, capture_output=True, text=True)
    print(rd, "rc", r.returncode)


def periodic_planes(rd, h5name):
    with open(os.path.join(rd, "bcondConfig.yaml")) as fp:
        bc = yaml.safe_load(fp)
    per_ids = {int(v["physID"]) for v in bc.values() if isinstance(v, dict) and v.get("kind") == "periodic"}
    with h5py.File(os.path.join(rd, h5name), "r") as f:
        nP = int(f["MESH"].attrs["nPlanes"]); nN = int(f["MESH"].attrs["nNormalPlanes"])
        ip = [np.asarray(f[f"BCONDS/{k}/iPlanes"], dtype=np.int64) for k in f["BCONDS"] if int(k) in per_ids]
    mask = np.zeros(nP, bool)
    if ip:
        mask[np.concatenate(ip)] = True
    return mask, nP, nN


def load_dump(rd, tag, nP):
    p = os.path.join(rd, f"sgdump.{tag}")
    if not os.path.exists(p):
        return None
    a = np.fromfile(p, dtype=np.float32)
    if a.size % (3 * nP) != 0:
        raise SystemExit(f"{p}: 大きさ {a.size} が 3·nPlanes={3 * nP} の倍数でない")
    return a.reshape(-1, nP, 3)


def cmd_compare(scratch, out):
    L = []
    P = L.append
    P("#8a: スカラー GG の面寄与ダンプ (FORGE_DUMP_SCALARGRAD) を 3 バイナリで比較 (plan boundary-node-periodic-gradient-fix §5.1 #8a、codex result-2 M1)")
    P(f"run root: {os.path.abspath(scratch)}/r8a_<axi|trans>_<new|olddump|wrong>")
    P("比較は float32 のビット比較 (uint32 view)。面 = PLANES の index、周期半割面 = BCONDS/<physID>/iPlanes (kind periodic) の和集合")
    P("注: RUN_PROVENANCE の git_head は run_case.sh を置いたリポジトリ (forge-pgrad-new = f67fe877) の値で、バイナリの版ではない。"
      "バイナリは forge_sha256 で識別する (new = f67fe877、olddump = 1266aba1 + olddump.patch、wrong = f67fe877 + wrong.patch)")
    verdict = {}
    for case in CASES:
        h5name = FILES[case][0]
        P("")
        P(f"## {case}")
        for b in BINS:
            prov = os.path.join(run_dir(scratch, case, b), "RUN_PROVENANCE.txt")
            if os.path.exists(prov):
                P(f"{b}: " + " ".join(l.strip() for l in open(prov) if l.startswith(("forge_sha256", "git_head"))))
            log = os.path.join(run_dir(scratch, case, b), "forge_run.log")
            if os.path.exists(log):
                for l in open(log):
                    if "FORGE_DUMP_SCALARGRAD" in l or "buildPeriodicNodeGroups" in l or "periodic seam" in l:
                        P(f"  {b} log: {l.rstrip()}")
        for b in BINS:
            r1 = os.path.join(run_dir(scratch, case, b), "res_1.h5")
            if os.path.exists(r1):
                with h5py.File(r1, "r") as f:
                    nbad = sum(int((~np.isfinite(np.asarray(f["VALUE"][k]))).sum()) for k in f["VALUE"])
                P(f"  {b} res_1.h5 の VALUE 非有限: {nbad}")
        mask, nP, nN = periodic_planes(run_dir(scratch, case, "new"), h5name)
        P(f"nPlanes {nP} (内部 {nN})、周期半割面 {int(mask.sum())}、それ以外 {int((~mask).sum())}")
        ok_case = True
        det = None
        for tag in ("species", "passive"):
            D = {b: load_dump(run_dir(scratch, case, b), tag, nP) for b in BINS}
            miss = [b for b in BINS if D[b] is None]
            if miss:
                P(f"- {tag}: ダンプが無い ({','.join(miss)}) → 判定不能")
                ok_case = False
                continue
            nv = D["new"].shape[0]
            if any(D[b].shape != D["new"].shape for b in BINS):
                P(f"- {tag}: 形が合わない {[D[b].shape for b in BINS]} → 判定不能")
                ok_case = False
                continue
            U = {b: D[b].view(np.uint32) for b in BINS}
            face_diff = lambda a, c: np.any(U[a] != U[c], axis=2)        # (nVar, nPlanes)
            zero_face = {b: np.all(D[b] == 0, axis=2) for b in BINS}      # 3 成分とも 0 (±0)
            nonfin = {b: int((~np.isfinite(D[b])).sum()) for b in BINS}
            P(f"- {tag}: nVar {nv}。非有限 " + ", ".join(f"{b} {nonfin[b]}" for b in BINS))
            P("  | 比較 | 周期半割面で差のある (変数,面) | それ以外で差のある (変数,面) |")
            P("  | --- | --- | --- |")
            for a, c in (("olddump", "new"), ("wrong", "new"), ("olddump", "wrong")):
                fd = face_diff(a, c)
                P(f"  | {a} vs {c} | {int(fd[:, mask].sum())} / {nv * int(mask.sum())} | {int(fd[:, ~mask].sum())} / {nv * int((~mask).sum())} |")
            P("  | バイナリ | 周期半割面で 3 成分とも 0 の (変数,面) | それ以外で 3 成分とも 0 の (変数,面) |")
            P("  | --- | --- | --- |")
            for b in BINS:
                P(f"  | {b} | {int(zero_face[b][:, mask].sum())} / {nv * int(mask.sum())} | {int(zero_face[b][:, ~mask].sum())} / {nv * int((~mask).sum())} |")
            if case == "axi":
                ok = not face_diff("olddump", "new").any()
                P(f"  判定 (axi: olddump = new 全面ビット同一): {'ok' if ok else 'NG'}")
                d_ok = bool(zero_face["wrong"][:, mask].all()) and bool(face_diff("wrong", "new")[:, mask].all()) \
                    and not face_diff("wrong", "new")[:, ~mask].any()
                P(f"  検出力 (wrong は周期半割面が全て 0・その全てで new と差、他の面は new と同一): {'ok' if d_ok else 'NG'}")
                det = d_ok if det is None else (det and d_ok)
            else:
                ok_other = not face_diff("olddump", "new")[:, ~mask].any()
                ok_new0 = bool(zero_face["new"][:, mask].all())
                ok_old = not zero_face["olddump"][:, mask].any()
                ok = ok_other and ok_new0 and ok_old
                P(f"  判定 (trans: 周期半割面以外で olddump = new ビット同一 {'ok' if ok_other else 'NG'}、"
                  f"周期半割面 new = 0 {'ok' if ok_new0 else 'NG'}、olddump ≠ 0 {'ok' if ok_old else 'NG'}): {'ok' if ok else 'NG'}")
            ok_case &= ok
        verdict[case] = ok_case
        if case == "axi":
            verdict["detect"] = bool(det)
    P("")
    P(f"VERDICT #8a axi (軸対称 × 周期、olddump = new 全面ビット同一): {'PASS' if verdict.get('axi') else 'FAIL'}")
    P(f"VERDICT #8a trans (並進周期、半割面以外ビット同一・半割面 new 0 / olddump 非 0): {'PASS' if verdict.get('trans') else 'FAIL'}")
    P(f"VERDICT #8a 検出力 (wrong は axi で周期半割面の寄与 0 → new と差): {'PASS' if verdict.get('detect') else 'FAIL'}")
    allok = verdict.get("axi") and verdict.get("trans") and verdict.get("detect")
    P(f"VERDICT #8a (3 条件): {'PASS' if allok else 'FAIL'}")
    txt = "\n".join(L) + "\n"
    with open(out, "w") as fp:
        fp.write(txt)
    sys.stdout.write(txt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("prepare", "run", "compare"))
    ap.add_argument("scratch")
    ap.add_argument("args", nargs="*")
    ap.add_argument("--out", default=os.path.join(HERE, "R8a_facedump.txt"))
    a = ap.parse_args()
    s = os.path.abspath(a.scratch)
    if a.cmd == "prepare":
        cmd_prepare(s)
    elif a.cmd == "run":
        cmd_run(s, a.args[0], a.args[1])
    else:
        cmd_compare(s, a.out)


if __name__ == "__main__":
    main()
