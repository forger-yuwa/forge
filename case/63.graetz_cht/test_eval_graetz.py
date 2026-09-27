#!/usr/bin/env python3
r"""`eval_graetz.py` の自己試験 (流体計算なし)。plan `boundary-cht-axisymmetric-graetz.md` §5.1 #5。

合成 run (N_r=16 のメッシュ、乾式の壁ダンプ座標) を**作業ディレクトリに実ファイルで**作る (元の run には触らない。
2026-09-26 の事故: h5 以外もリンクして成果物が元 run を上書きした)。

- 加熱 run: ρ 一定・放物速度・Poiseuille 圧力、内部の T は各 x で基準解の混合平均 T_b(x⁺) (半径方向一様)、
  壁は T_w = T_c、q = −k Nu_ref (T_c − T_b)/D (Nu_ref は評価器と同じ双対面区間平均)
- 対照 run: T = T_in 一様、q0 = 0、T_w0 = T_in

合否は**判定語と終了コードの両方**で見る (「(no VERDICT)」を合格と取り違えない)。

    python3 case/63.graetz_cht/test_eval_graetz.py [--work DIR]
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

import eval_graetz as ev
import graetz_common as gc

HERE = Path(__file__).resolve().parent
MESH = HERE / "mesh" / "graetz_r16.h5"
DRY = HERE / "run_0001_dry_r16" / "res_wall_heat_4_1.h5"
DT = 10.0


def write_run(d: Path, T, wx, q, tw, ok, steps=(1000,), drift=None):
    """run ディレクトリ d に mesh.h5 (実コピー) と res_<st>.h5 / res_wall_heat_4_<st>.h5 を書く。
    drift = (節点 index, 相対変化の最終値) なら、スナップショットごとにその節点の q を線形に動かす。"""
    d.mkdir(parents=True)
    shutil.copy(MESH, d / "mesh.h5")
    with h5py.File(MESH, "r") as m:
        xyz = np.asarray(m["MESH/COORD"][:], float).reshape(-1, 3)
    x, r = xyz[:, 0], xyz[:, 1]
    u = 2 * gc.U_M * np.clip(1 - (r / gc.R) ** 2, 0, None)
    P = gc.P_OUT + 8 * gc.MU * gc.U_M / gc.R ** 2 * (gc.L_HEAT + gc.L_DOWN - x)
    ro = np.full_like(x, gc.RHO)
    with h5py.File(DRY, "r") as w:
        wc = np.asarray(w["MESH/COORD"][:], float)
        wconne = np.asarray(w["MESH/CONNE"][:])
    wxyz = wc.reshape(-1, 3)
    order = np.argsort(wxyz[:, 0])                       # 壁ダンプの並び → x 昇順の対応
    inv = np.empty_like(order); inv[order] = np.arange(len(order))
    for k, st in enumerate(steps):
        with h5py.File(d / f"res_{st}.h5", "w") as h:
            for name, v in (("ro", ro), ("Ux", u), ("Uy", 0 * u), ("T", T), ("P", P)):
                h.create_dataset(f"VALUE/{name}", data=v)
        qk = q.copy()
        if drift is not None:
            i, amp = drift
            qk[i] *= 1 + amp * k / max(1, len(steps) - 1)
        with h5py.File(d / f"res_wall_heat_4_{st}.h5", "w") as h:
            h.create_dataset("MESH/COORD", data=wc)
            h.create_dataset("MESH/CONNE", data=wconne)
            # 壁ダンプは元の並び (x 昇順でない) で書く: 評価器が並べ替えることも試す
            h.create_dataset("VALUE/iface_q_eff", data=qk[inv])
            h.create_dataset("VALUE/iface_Tw_bc", data=tw[inv])
            h.create_dataset("VALUE/iface_ok", data=ok[inv].astype(float))


def build(work: Path):
    run = ev.Run.__new__(ev.Run)
    with h5py.File(MESH, "r") as m:
        xyz = np.asarray(m["MESH/COORD"][:], float).reshape(-1, 3)
    wx = None
    with h5py.File(DRY, "r") as w:
        wx = np.sort(np.asarray(w["MESH/COORD"][:], float).reshape(-1, 3)[:, 0])
    # 評価器と同じ Pe (Simpson の ṁ) を得るため、放物速度・ρ 一定の列積分から Pe を作る
    r = np.linspace(0, gc.R, 17)
    m = ev.rint(gc.RHO * 2 * gc.U_M * (1 - (r / gc.R) ** 2) * r, r)
    Pe = (2 * m / gc.R ** 2) * gc.D / gc.MU * gc.PR
    xs, nu_m, tb_m = ev.reference()
    Tc = gc.T_IN + DT

    def fields(pe_true, origin_shift=0.0):
        xp_nodes = (xyz[:, 0] + origin_shift) / (gc.D * pe_true)
        thb = np.where(xp_nodes > 0, np.interp(np.log(np.clip(xp_nodes, 1e-9, None)), np.log(xs), tb_m), 1.0)
        T = Tc - (Tc - gc.T_IN) * thb
        xpw = (wx + origin_shift) / (gc.D * pe_true)
        nu_ref = ev.ref_for_nodes(xpw)
        thw = np.interp(np.log(np.clip(xpw, 1e-9, None)), np.log(xs), tb_m)
        Tb_w = Tc - (Tc - gc.T_IN) * thw
        q = -gc.K_F * nu_ref * (Tc - Tb_w) / gc.D
        q[~np.isfinite(q)] = 0.0
        return T, q

    ok = np.ones(len(wx), bool)
    T, q = fields(Pe)
    Tctl = np.full(len(xyz), gc.T_IN)
    zero = np.zeros(len(wx))
    tw_h, tw_c = np.full(len(wx), Tc), np.full(len(wx), gc.T_IN)
    steps6 = tuple(range(1000, 7000, 1000))
    write_run(work / "ctl", Tctl, wx, zero, tw_c, ok, steps=steps6)
    write_run(work / "ok", T, wx, q, tw_h, ok, steps=steps6)
    write_run(work / "signflip", T, wx, -q, tw_h, ok)
    T3, q3 = fields(Pe, origin_shift=gc.L_UP)            # x⁺ の原点を入口にした取り違え
    write_run(work / "origin", T3, wx, q3, tw_h, ok)
    T4, q4 = fields(2 * Pe)                              # Pe の取り違え
    write_run(work / "pe2", T4, wx, q4, tw_h, ok)
    ok5 = ok.copy()
    xp = wx / (gc.D * Pe)
    iw = np.where((xp >= 3e-3) & (xp <= 0.1))[0]
    ok5[iw[len(iw) // 2]] = False
    write_run(work / "ifaceok0", T, wx, q, tw_h, ok5)
    write_run(work / "wrongctl", Tctl, wx, zero, tw_c, ok)          # 加熱のつもりが対照と同じ (分母 0)
    # 4 観測点 (3e-3, 1e-2, 3e-2, 0.1) から最も遠い窓内節点を 1 % ドリフトさせる
    lp = np.log(np.array(ev.PTS))
    far = iw[np.argmax([np.min(np.abs(np.log(xp[i]) - lp)) for i in iw])]
    write_run(work / "drift", T, wx, q, tw_h, ok, steps=steps6, drift=(far, 0.01))
    # 最大誤差節点を壁ダンプから落とす: ok run の壁ダンプから 1 節点を削除した run
    d = work / "missing"
    write_run(d, T, wx, q, tw_h, ok)
    for f in d.glob("res_wall_heat_4_*.h5"):
        with h5py.File(f, "a") as h:
            for k in ("iface_q_eff", "iface_Tw_bc", "iface_ok"):
                v = h[f"VALUE/{k}"][:]; del h[f"VALUE/{k}"]; h.create_dataset(f"VALUE/{k}", data=v[:-1])
            c = h["MESH/COORD"][:].reshape(-1, 3)[:-1].ravel(); del h["MESH/COORD"]; h.create_dataset("MESH/COORD", data=c)
    return Pe, float(xp[far])


def run_eval(mode, heated, ctl):
    p = subprocess.run([sys.executable, str(HERE / "eval_graetz.py"), mode, str(heated), str(ctl)],
                       capture_output=True, text=True)
    out = p.stdout + p.stderr
    v = [l for l in out.splitlines() if l.startswith("VERDICT")]
    word = v[-1].split(":")[1].strip().split()[0] if v else "(no VERDICT)"
    return word, p.returncode, out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", default=None)
    a = ap.parse_args()
    work = Path(a.work) if a.work else Path(tempfile.mkdtemp(prefix="graetz_evaltest_"))
    if any(work.iterdir()) if work.exists() else False:
        raise SystemExit(f"{work} が空でない")
    work.mkdir(parents=True, exist_ok=True)
    Pe, xfar = build(work)
    print(f"合成 run: {work}  (Pe {Pe:.6f}、ドリフト節点 x⁺ {xfar:.4g})")
    cases = [  # (名前, mode, heated, 期待する判定語, 期待する終了コード, 説明)
        ("正常 (snap)", "snap", "ok", "PASS", 0, "基準解そのものの合成データ"),
        ("正常 (series)", "series", "ok", "PASS", 0, "6 スナップショット同一"),
        ("q の符号反転", "snap", "signflip", "FAIL", 1, "Nu が負"),
        ("x⁺ の原点を入口に", "snap", "origin", "FAIL", 1, "10 mm のずれ = x⁺ 0.0069"),
        ("Pe を 2 倍", "snap", "pe2", "FAIL", 1, ""),
        ("窓内 iface_ok=0", "snap", "ifaceok0", "REFUSED", 2, ""),
        ("壁ダンプの節点欠損", "snap", "missing", "REFUSED", 2, ""),
        ("対照と同じ run を加熱として渡す", "snap", "wrongctl", "FAIL", 1, "分母 0 → 非有限"),
        ("観測 4 点を避けた局所ドリフト 1 %", "series", "drift", "FAIL", 1, ""),
    ]
    bad = 0
    for name, mode, h, want, rc_want, note in cases:
        word, rc, out = run_eval(mode, work / h, work / "ctl")
        good = (word == want and rc == rc_want)
        bad += 0 if good else 1
        print(f"  {'ok ' if good else 'NG '} {name:<34} 判定 {word:<9} rc {rc}  (期待 {want}/{rc_want}) {note}")
        if not good:
            print("      " + "\n      ".join(out.splitlines()[-12:]))
    # 正常 run の誤差が丸め程度であること (評価器が基準解を再現する)
    csv = sorted((work / "ok").glob("graetz_nu_*.csv"))[-1]
    d = np.genfromtxt(csv, delimiter=",", names=True)
    xp = d["xplus"]; W = (xp >= 3e-3) & (xp <= 0.1)
    e = np.abs(d["err"][W]).max()
    # 許容 1e-6: 合成側 (T_b を対数補間) と評価側 (Simpson) の経路差の丸め。初版の 1e-9 は根拠なく置いた値で、
    # 実測 2.5e-9 (2026-09-27) で落ちたので、合否許容 2 % より 4 桁下の 1e-6 に直した
    good = e < 1e-6
    bad += 0 if good else 1
    print(f"  {'ok ' if good else 'NG '} 正常 run の窓内 max|err| = {e:.2e} (期待 < 1e-6)")
    print(f"\nVERDICT: {'PASS' if bad == 0 else 'FAIL'}  ({len(cases) + 1 - bad}/{len(cases) + 1})")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
