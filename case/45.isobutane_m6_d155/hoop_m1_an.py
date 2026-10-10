#!/usr/bin/env python3
"""plan axisymmetric-freestream-hoop-gauge §4.11 (m1、事前登録) の判定: hoop_m1.sh のログ (hoop_verify_conv.py check --rw の出力) を読む。
合格: |ΔW|/|W| ≤ 1e-12 (内部の面・境界の半割面)、対応の無い半割面 0、r 重みの閉性の超過 0 (x・y)、W·S ≤ 0 の面 0、平面の閉性の超過 0 (x・y)。
行が欠ける・変換が失敗したら判定不能。"""
import re, sys
log = open(sys.argv[1]).read()
blocks = re.split(r"^== (T\d_\w+) 変換 rc=(\d+)\n", log, flags=re.M)
res = {}
for i in range(1, len(blocks), 3):
    T, rc, body = blocks[i], int(blocks[i + 1]), blocks[i + 2]
    g = lambda pat: re.findall(pat, body)
    pl = g(r"planar closure Σ±S_[xy]: .*?CVs over 100·ε64·Σ\|S\|: (\d+)")
    rw = g(r"r-weighted closure \(rSurfVect\): max E_x (\S+)\s+max E_y (\S+)\s+CVs over 100·ε64·\(A\+Σ\|W\|\): x (\d+) y (\d+)")
    ref = g(r"\(reference\) r̄·S on the same file: max E_x (\S+)\s+max E_y (\S+)")
    ori = g(r"faces with W·S <= 0 among W != 0: (\d+)")
    wi = g(r"interior faces: rebuilt W vs rSurfVect max \|ΔW\|/\|W\| (\S+)")
    wb = g(r"boundary half faces \((\d+), unmapped (\d+)\): rebuilt W vs rSurfVect max \|ΔW\|/\|W\| (\S+?);")
    nn = g(r"nNodes (\d+) nPlanes (\d+) \(interior (\d+)\)")
    if rc != 0 or len(pl) != 2 or not (rw and ref and ori and wi and wb and nn):
        res[T] = ("判定不能", f"変換 rc={rc} または照合の行が欠けた"); continue
    ex, ey, ox, oy = rw[0]; nb, um, dwb = wb[0]
    ok = (float(wi[0]) <= 1e-12 and float(dwb) <= 1e-12 and um == "0" and ox == "0" and oy == "0" and ori[0] == "0" and pl == ["0", "0"])
    res[T] = ("PASS" if ok else "FAIL",
              f"節点 {nn[0][0]}・面 {nn[0][1]} (内部 {nn[0][2]}、半割 {nb})、|ΔW|/|W| 内部 {wi[0]}・半割 {dwb} (対応なし {um})、"
              f"r 重みの閉性 max E_x {ex}・E_y {ey} (超過 x {ox}・y {oy})、W·S≤0 {ori[0]}、平面の閉性の超過 {pl}、参考 r̄·S max E_x {ref[0][0]}・E_y {ref[0][1]}")
for T in ("T1_tri_curved", "T2_mixed", "T3_step_multicorner"):
    v, d = res.get(T, ("判定不能", "ログに無い"))
    print(f"[{v}] {T}: {d}")
bad = [T for T in ("T1_tri_curved", "T2_mixed", "T3_step_multicorner") if res.get(T, ("",))[0] != "PASS"]
print("== VERDICT §4.11 (m1):", "ALL PASS" if not bad else f"PASS でないもの {bad}")
