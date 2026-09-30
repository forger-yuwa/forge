#!/usr/bin/env python3
r"""case/65 B1 本段の毎 step・FP64 系列を、forge の通常出力から走行中に抜き出す (plan `boundary-cht-conjugate-flat-plate.md` §4.1.2)。

forge は毎 step `res_<N>.h5` (流体) と `res_plate_5_<N>.h5` (壁、`iface_q_eff`・`iface_Tw_bc`) を書く。
全場を 60000 step 分は残せないので、本ツールが並走して登録節点の値を FP64 のまま抜き出し、処理済みの h5 を消す。
**ソルバには手を入れない** — 値は既存の壁ダンプと同じ診断式・同じ出力時点のもの (probe の 6 桁出力を使わない。
2026-10-01 codex diagnose `sampling-redesign`)。

    python3 ab_extract.py <run> --until <最終出力 step> [--keep-every 2000]

- 処理するのは「次の step のファイルが出来た」か「forge が終わった」step だけ (書きかけを読まない)。
- `--keep-every` の倍数の step と最初・最後の step は全場を残す (後で場を確かめるため)。
- 出力: `<run>/ab_series.npz` (step、流体の `--fields` [step, 124]、界面 q/Tw [step, 161]、節点 ID・座標)。
  1000 step ごとに書き直す (途中で落ちても残る)。同じ step の欠落・重複は REFUSED。
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ab_series import node_sets, refuse  # noqa: E402
import plate_common as pc  # noqa: E402


def forge_running(run):
    p = subprocess.run(["pgrep", "-x", "forge"], capture_output=True, text=True)
    for pid in p.stdout.split():
        try:
            if Path(f"/proc/{pid}/cwd").resolve() == Path(run).resolve():
                return True
        except OSError:
            pass
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run"); ap.add_argument("--until", type=int, required=True); ap.add_argument("--keep-every", type=int, default=2000)
    ap.add_argument("--fields", default="P,T,Uy", help="流体節点で採る VALUE 名 (カンマ区切り)。B2′ は P,T,Uy,Ux,ro,dt_local,limiter_ro,limiter_Ux,limiter_Uy,limiter_P,limiter_T")
    ap.add_argument("--region-stats", action="store_true",
                    help="B3: 全域の res_ro/res_roUx/res_roUy/res_roe の領域別二乗和を毎 step、limiter_ro/Ux/Uy/P と res² の節点ごとの統計 (最終 1/3 区間) を保存")
    a = ap.parse_args()
    run = Path(a.run)
    while not (run / "res_0.h5").exists():
        time.sleep(2)
    time.sleep(2)
    sets, _, _ = node_sets(run)
    order = [k for k in ("le", "te", "up", "down", "le_in", "te_in", "up_in", "down_in")]
    ids = np.concatenate([sets[k] for k in order])
    groups = np.concatenate([[k] * len(sets[k]) for k in order])
    fields = a.fields.split(",")
    steps, Q, TW = [], [], []
    FV = {k: [] for k in fields}
    wx = None
    n = 0
    out = run / "ab_series.npz"
    RES = ("res_ro", "res_roUx", "res_roUy", "res_roe")
    LIM = ("limiter_ro", "limiter_Ux", "limiter_Uy", "limiter_P")
    RS = []
    if a.region_stats:
        with h5py.File(run / "res_0.h5", "r") as r:
            cc = np.asarray(r["MESH/COORD"][:], float).reshape(-1, 3)
        xl, yl = cc[:, 0] / pc.L, cc[:, 1] / pc.L
        near = yl <= 0.05 + 1e-12
        reg_le = near & (xl >= -0.1 - 1e-9) & (xl <= 0.02 + 1e-9)
        reg_te = near & (xl >= 0.98 - 1e-9) & (xl <= 1.1 + 1e-9)
        regions = np.stack([reg_le, reg_te, ~(reg_le | reg_te)])
        nn = len(xl)
        t0 = a.until - a.until // 3          # 最終 1/3 区間の開始
        acc = {"res2": np.zeros((len(RES), nn)), "lim_s": np.zeros((len(LIM), nn)), "lim_ss": np.zeros((len(LIM), nn)), "cnt": 0}

    def save():
        extra = {}
        if a.region_stats:
            extra = {"region_res2": np.array(RS), "region_names": np.array(["le", "te", "other"]), "region_count": regions.sum(axis=1),
                     "res_names": np.array(RES), "lim_names": np.array(LIM), "n_nodes": nn,
                     "tail_res2_sum": acc["res2"], "tail_lim_sum": acc["lim_s"], "tail_lim_sumsq": acc["lim_ss"], "tail_count": acc["cnt"],
                     "coord": cc}
        np.savez(out, step=np.array(steps), q=np.array(Q), Tw=np.array(TW), node_id=ids, group=groups, wall_x=wx,
                 **{k: np.array(v) for k, v in FV.items()}, **extra)

    while n <= a.until:
        f, w = run / f"res_{n}.h5", run / f"res_plate_5_{n}.h5"
        nxt = run / f"res_{n + 1}.h5"
        done = not forge_running(run)
        if not (f.exists() and (n == 0 or w.exists()) and (nxt.exists() or done)):
            if done and not f.exists():
                refuse(f"forge が終わったのに res_{n}.h5 が無い")
            time.sleep(0.5)
            continue
        with h5py.File(f, "r") as r:
            miss = [k for k in fields if f"VALUE/{k}" not in r]
            if miss:
                refuse(f"step {n}: 出力に {miss} が無い (extraFields を確認)")
            v = {q: np.asarray(r[f"VALUE/{q}"][:], float)[ids] for q in fields}
            if a.region_stats and n > 0:
                miss = [k for k in RES + LIM if f"VALUE/{k}" not in r]
                if miss:
                    refuse(f"step {n}: 出力に {miss} が無い (extraFields を確認)")
                r2 = np.array([np.asarray(r[f"VALUE/{k}"][:], float) ** 2 for k in RES])
                if not np.isfinite(r2).all():
                    refuse(f"step {n}: 残差場に非有限値")
                RS.append(np.array([[r2[j][m].sum() for m in regions] for j in range(len(RES))]))
                if n > t0:
                    lm = np.array([np.asarray(r[f"VALUE/{k}"][:], float) for k in LIM])
                    acc["res2"] += r2; acc["lim_s"] += lm; acc["lim_ss"] += lm ** 2; acc["cnt"] += 1
        if n > 0:
            with h5py.File(w, "r") as h:
                c = np.asarray(h["MESH/COORD"][:], float).reshape(-1, 3)
                o = np.argsort(c[:, 0])
                if wx is None:
                    wx = c[o, 0]
                q = -np.asarray(h["VALUE/iface_q_eff"][:], float)[o]
                tw = np.asarray(h["VALUE/iface_Tw_bc"][:], float)[o]
        else:
            q = tw = None
        if not all(np.isfinite(x).all() for x in v.values()) or (q is not None and not (np.isfinite(q).all() and np.isfinite(tw).all())):
            refuse(f"step {n}: 非有限値")
        if n > 0:
            steps.append(n); Q.append(q); TW.append(tw)
            for k in fields:
                FV[k].append(v[k])
        keep = n == 0 or n == a.until or n % a.keep_every == 0
        if not keep:
            for p in (f, w, f.with_suffix(".xmf"), w.with_suffix(".xmf"), run / f"res_solid_5_{n}.h5", run / f"res_solid_5_{n}.xmf"):
                p.unlink(missing_ok=True)
        if n % 1000 == 0:
            save()
        n += 1
    save()
    s = np.array(steps)
    if len(s) != a.until or (np.diff(s) != 1).any():
        refuse(f"step の欠落・重複 ({len(s)} 点)")
    print(f"VERDICT: OK ({len(s)} step、流体 {len(ids)} 節点、界面 {len(wx)} 節点 → {out.name})")


if __name__ == "__main__":
    main()
