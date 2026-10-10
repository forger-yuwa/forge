#!/usr/bin/env python3
"""方向別 dt の上限の A/B (plan time_integration-line-viscous-jacobian-dt §6) の事前のゲート。本試験 (run_0580〜0583、条件付きで 0588・0589) の前に回し、合格したときだけ本試験へ進む。

同じ状態 (run_0183 の res_100000) から段 ③ の FP64・値 3・マスク 7 で 1 step の書き出しを 4 本読む:
  a  = run_0584_lvcdt_a_dump  (方向別 dt、上限なし)
  a2 = run_0586_lvcdt_a_dump2 (a の再実行、atomicAdd による揺れを測る)
  b  = run_0585_lvcdt_b_dump  (方向別 dt、上限 50: lineDtDirectionalCap 50)
  c  = run_0587_lvcdt_c_dump  (方向別なし = point の dt: lineDtDirectional を書かない)
と、faceh §6.7 の監査の記録 (run_0550_lvcaudit、同じ状態・同じ 605 節点、上限なしの方向別 dt) を使い、
  I  入力のゲート (faceh §6.12 の lvc75_pregate.py と同じ項目。設定の期待値は腕ごと、4 本の設定は dt の 2 キー以外で同じ)
  R  残差の不変: 出発の状態の残差 (全節点、倍精度) が、a の再実行どうしでビット一致なら b・c もビット一致、そうでなければ ≤ 3 × 再実行の差
  S  構造の照合: 各書き出しで rhs_s0 = 拘束の処理をした float(出力 step 1 の残差) が全要素で完全に一致
  V  介入の成立 (b・c のそれぞれを a に対して): (a) Kprev・Knext が不変、(b) r = Δτ_a/Δτ_x が全節点で ≥ 1 − 4 ulp、縮流部の第一内部節点 (1571・4112・7984) で ≥ 10、
     全節点で r = 1 なら INVALID、r_c ≥ r_b、(c) D の対角以外と拘束の行はビット一致、拘束の無い行の対角は時間項 V/Δτ の差と 64 ulp 以内、(c') 第一内部節点で時間項の差が、拘束の無い行のどれかで許容の 100 倍超
を判定する。出力は _band_ab/cold_pair/lvcdt_pregate.json。終了コード: 0 = 合格 (本試験へ)、1 = INVALID (入力・証拠の不備)、2 = 判別不能 (残差の変化・構造の外れ・介入の不成立)。
使い方: lvcdt_pregate.py [--verify]
"""
import hashlib
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
SRC, SRC_RES, SRC_SHA16 = "run_0183_ns_coldmesh_tw300_ext", "res_100000.h5", "207d39f0e7f4aa03"
SHA_NEW = "129de3f4e7f67aa3a80dbd30d5f5cb75df1582d3974e702d76b61c8998598cec"
DUMPS = {"a": ("run_0584_lvcdt_a_dump", SHA_NEW), "a2": ("run_0586_lvcdt_a_dump2", SHA_NEW), "b": ("run_0585_lvcdt_b_dump", SHA_NEW), "c": ("run_0587_lvcdt_c_dump", SHA_NEW)}
MASK = {k: 7 for k in DUMPS}
# 腕ごとの設定の差 (EXPECT を上書きする。None は「書いていない」を要求)
ARM_CFG = {"a": {}, "a2": {}, "b": {"time/deltaT/lineDtDirectionalCap": 50.0},
           "c": {"time/deltaT/lineDtDirectional": None}}
DT_KEYS = {"time/deltaT/lineDtDirectional", "time/deltaT/lineDtDirectionalCap"}
WALL_FIRST = [1571, 4112, 7984]          # 縮流部の 3 本の壁法線のライン (1572・4113・7985 の列) の第一内部節点
AUDIT = "run_0550_lvcaudit"
WANT_NODES = [1572, 4113, 7985, 198560, 264263]
RES_FIELDS = ["res_ro", "res_roUx", "res_roUy", "res_roUz", "res_roe", "res_roK", "res_roOmega"]   # res_roUz は rhs_s0 の行 3 の照合に使う (§6.12)
RES_OPT = ["res_roY0", "res_roY1"]          # 化学種の残差 (3 本ともにあれば R に含める)
EXTRA = RES_FIELDS + RES_OPT + ["volume"]
SAME_FILES = ["bcondConfig.yaml", "probe.yaml", "species_meta.yaml", "wall_design.csv", "wall_physical.csv", "target_axis_M.csv", "wall_repr.json", "MESH_QUALITY.txt"]
EXPECT = {"time/deltaT/cfl": 4.0, "time/deltaT/cfl_pseudo": 4.0, "time/deltaT/implicitRelax": 0.7, "time/deltaT/implicitThermalJacobian": 5,
          "time/deltaT/lineDtDirectional": 1, "time/deltaT/lineImplicit": 1, "time/deltaT/lineViscCoupling": 3, "time/deltaT/blockDPLUR": 1,
          "time/deltaT/detectNaN": 1, "time/timeIntegration": 11, "time/nStepInner": 5, "time/last/nStepOuter": 1,
          "time/deltaT/lineDtDirectionalCap": None, "time/deltaT/implicitSolvePrecision": None, "time/outStepInterval": 1,
          "output/extraFields": EXTRA}
NF, FREC, NREC, HF = 12, 200, 160, 150
f32 = np.float32


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def meshsha(p):
    h = hashlib.sha256()

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            h.update(name.encode()); h.update(obj[...].tobytes())
    with h5py.File(p, "r") as f:
        f["MESH"].visititems(visit)
    return h.hexdigest()


def flat(d, p=""):
    out = {}
    for k, v in d.items():
        kk = f"{p}/{k}" if p else str(k)
        out.update(flat(v, kk) if isinstance(v, dict) else {kk: v})
    return out


def dump(run, name):
    d = HERE / run / "linedump"
    meta = {l.split()[0]: (int(l.split()[1]), int(l.split()[2])) for l in (d / "meta.txt").read_text().splitlines() if l.strip() and not l.startswith("#")}
    r, c = meta[name]
    a = np.fromfile(d / f"{name}.f64", dtype=np.float64)
    if a.size != r * c:
        raise ValueError(f"{run}/{name}: 大きさが違う")
    return a.reshape(r, c)


EVIDENCE_FILES = ["linedump/meta.txt", "linedump/D.f64", "linedump/Kprev.f64", "linedump/Knext.f64", "linedump/state_ro_roU_roe_cp_gamma.f64",
                  "linedump/dt_vol.f64", "linedump/flags_wall_iso_axis.f64", "linedump/node_line.f64", "linedump/rhs_s0.f64",
                  "res_0.h5", "res_1.h5", "solverConfig.yaml", "LVC_TERMS.txt", "forge_run.log", "RUN_PROVENANCE.txt"]
AUDIT_FILES = ["linedump/audit_face.f64", "linedump/audit_node.f64", "linedump/meta.txt", "linedump/node_line.f64",
               "linedump/state_ro_roU_roe_cp_gamma.f64", "linedump/dt_vol.f64", "linedump/flags_wall_iso_axis.f64"]


def evidence():
    """事前のゲートが読んだ証拠と、このスクリプト自身・書き出しの設定のハッシュ (本試験の起動と本判定で照合する)"""
    ev = {f"{run}/{fn}": sha(HERE / run / fn) for run, _ in DUMPS.values() for fn in EVIDENCE_FILES}
    ev.update({f"{AUDIT}/{fn}": sha(HERE / AUDIT / fn) for fn in AUDIT_FILES})
    ev["lvcdt_pregate.py"] = sha(Path(__file__).resolve())
    return ev


def verify():
    """--verify: 記録が PASS で、記録したハッシュが今のファイルと同じなら 0、そうでなければ 1"""
    rec = json.loads((HERE / "_band_ab" / "cold_pair" / "lvcdt_pregate.json").read_text())
    if not str(rec.get("VERDICT", "")).startswith("PASS"):
        print("事前のゲートの記録が PASS でない:", rec.get("VERDICT")); return 1
    ev = rec.get("evidence_sha256") or {}
    if not ev:
        print("証拠のハッシュが記録に無い"); return 1
    now = evidence()
    bad = sorted(k for k in set(ev) | set(now) if ev.get(k) != now.get(k))
    if bad:
        print("記録と違う証拠:", bad[:10]); return 1
    print(f"照合 OK: {len(ev)} ファイル"); return 0


def ulp32(x):
    return np.spacing(np.abs(np.asarray(x, dtype=np.float32))).astype(np.float64)


def main():
    out = HERE / "_band_ab" / "cold_pair" / "lvcdt_pregate.json"
    rec = {"plan": "time_integration-line-viscous-jacobian-dt §6 (事前のゲート)"}
    checks = []                 # (段, 名前, ok, 値)

    def ok(stage, name, cond, val=None):
        checks.append((stage, name, bool(cond), val))
    verdict = None
    try:
        # ---- I: 入力 ----
        src_full = sha(HERE / SRC / SRC_RES)
        ok("I", "出発の場の sha256", src_full.startswith(SRC_SHA16), src_full[:16])
        mref = meshsha(HERE / SRC / "nozzle.h5")
        for key, (run, bsha) in DUMPS.items():
            d = HERE / run
            prov = (d / "RUN_PROVENANCE.txt").read_text(errors="replace")
            ok("I", f"{run}: バイナリの sha256", bsha in prov, bsha[:16])
            cp = json.loads((d / "COLD_PAIR.json").read_text())
            ok("I", f"{run}: 出発の場", cp["parent"] == SRC and cp["parent_res"] == SRC_RES and cp["parent_res_sha256"] == src_full, cp["parent_res_sha256"][:16])
            m = meshsha(d / "nozzle.h5")
            ok("I", f"{run}: 格子の実体", m == mref and (d / "MESH_SHA.txt").read_text().strip() == m, m[:16])
            ok("I", f"{run}: 入力ファイル", all((d / fn).is_file() and sha(d / fn) == sha(HERE / SRC / fn) for fn in SAME_FILES), None)
            c = flat(yaml.safe_load((d / "solverConfig.yaml").read_text()))
            exp_ = {**EXPECT, **ARM_CFG[key]}
            bad = {k: c.get(k) for k, v in exp_.items() if (k in c if v is None else c.get(k) != v)}
            ok("I", f"{run}: 設定 (腕ごとの期待どおり。本試験との差は step 数・出力の間隔・extraFields だけ)", bad == {}, bad)
            log = (d / "forge_run.log").read_text(errors="replace")
            mk = MASK[key]
            shown = ("FORGE_LVC_TERMS=" not in log) if mk == 7 else (f"診断のマスク FORGE_LVC_TERMS={mk} " in log)
            ok("I", f"{run}: 起動ログ (値 3・マスク {mk} の表示・LAYOUT2・面エンタルピーの切替なし・書き出しあり)",
               "'lineViscCoupling' in 'time.deltaT': 3" in log and shown and "Thomas の配列の並び: LAYOUT2" in log
               and "FORGE_DIAG_FACE_H_DOUBLE" not in log and "[lineDump] factor の直前を書いた" in log, None)
            lt = (d / "LVC_TERMS.txt").read_text().strip() if (d / "LVC_TERMS.txt").is_file() else None
            ok("I", f"{run}: 台本が渡した FORGE_LVC_TERMS = {mk}", lt == str(mk), lt)
            ok("I", f"{run}: RUN_RC = 0", (d / "RUN_RC").read_text().strip() == "0", None)
            ok("I", f"{run}: extraFields の無視の警告なし", "は確保されていない変数なので無視する" not in log, None)
        cfgs = {k: flat(yaml.safe_load((HERE / DUMPS[k][0] / "solverConfig.yaml").read_text())) for k in DUMPS}
        def strip(c):
            return {k: v for k, v in c.items() if k not in DT_KEYS}
        ok("I", "4 本の書き出しの設定が lineDtDirectional・lineDtDirectionalCap 以外で完全に同じ", all(strip(cfgs[k]) == strip(cfgs["a"]) for k in DUMPS), None)
        ok("I", "A と A の再実行の設定が完全に同じ", cfgs["a2"] == cfgs["a"], None)
        rec["dump_config"] = {k: cfgs[k] for k in DUMPS}
        # 書き出しの形・節点
        arrays = {}
        for key, (run, _) in DUMPS.items():
            d = HERE / run
            a = {n: dump(run, n) for n in ("node_line", "D", "Kprev", "Knext", "state_ro_roU_roe_cp_gamma", "dt_vol", "flags_wall_iso_axis", "rhs_s0")}
            nn = a["node_line"].shape[0]
            shapes = {"node_line": 2, "D": 25, "Kprev": 25, "Knext": 25, "state_ro_roU_roe_cp_gamma": 7, "dt_vol": 2, "flags_wall_iso_axis": 3, "rhs_s0": 5}
            ok("I", f"{run}: 配列の形", nn > 0 and all(a[n].shape == (nn, w) for n, w in shapes.items()), {n: a[n].shape for n in shapes})
            nodes = [int(x) for x in a["node_line"][:, 0]]
            ok("I", f"{run}: 要求した 5 節点・5 本のライン・重複なし", all(w in nodes for w in WANT_NODES) and len(set(a["node_line"][:, 1])) == 5 and len(set(nodes)) == nn, None)
            ok("I", f"{run}: 有限", all(bool(np.all(np.isfinite(a[n]))) for n in shapes), None)
            st = a["state_ro_roU_roe_cp_gamma"]
            ok("I", f"{run}: ρ・c_p・dt・体積 > 0、γ > 1", bool(np.all(st[:, 0] > 0) and np.all(st[:, 5] > 0) and np.all(st[:, 6] > 1) and np.all(a["dt_vol"] > 0)), None)
            with h5py.File(d / "res_1.h5", "r") as h:
                have = set(h["VALUE"].keys())
                r0 = h["VALUE"]["res_ro"][...] if "res_ro" in have else None
            ok("I", f"{run}: res_1.h5 に残差の場がそろう", all(f in have for f in RES_FIELDS), sorted(set(RES_FIELDS) - have))
            # 1 step の出力の残差が出発の状態の残差であること (rhs_s0 の行 0 = float(res_ro)、loop 0 は隣の寄与 0)
            ok("I", f"{run}: res_1 の res_ro が書き出しの rhs_s0 の行 0 と float で一致 (出発の状態の残差)",
               r0 is not None and np.array_equal(r0[[int(x) for x in a["node_line"][:, 0]]].astype(np.float32).astype(np.float64), a["rhs_s0"][:, 0]), None)
            arrays[key] = a
            arrays[key]["_have"] = have
        for n in ("node_line", "state_ro_roU_roe_cp_gamma", "flags_wall_iso_axis"):
            ok("I", f"{n}: 4 本でビット一致 (同じ入力)",
               all(arrays[k][n].shape == arrays["a"][n].shape and np.array_equal(arrays[k][n].view(np.int64), arrays["a"][n].view(np.int64)) for k in DUMPS), None)
        ok("I", "dt_vol: A と A の再実行でビット一致、体積の列は 4 本でビット一致",
           np.array_equal(arrays["a2"]["dt_vol"].view(np.int64), arrays["a"]["dt_vol"].view(np.int64))
           and all(np.array_equal(arrays[k]["dt_vol"][:, 1].view(np.int64), arrays["a"]["dt_vol"][:, 1].view(np.int64)) for k in DUMPS), None)
        # 場全体の入力 (保存量・乱流・組成・物性・壁距離): 出力 step 0 の全データセットが 3 本でビット一致
        def h5vals(run):
            with h5py.File(HERE / run / "res_0.h5", "r") as h:
                return {k: h["VALUE"][k][...] for k in h["VALUE"].keys()}
        v0 = {k: h5vals(DUMPS[k][0]) for k in DUMPS}
        keys0 = sorted(v0["a"].keys())
        ok("I", "res_0.h5 の全データセット (保存量・k・ω・組成・μ・壁距離など) が 4 本でビット一致",
           len(keys0) > 0 and all(sorted(v0[k].keys()) == keys0 and all(np.array_equal(v0[k][n], v0["a"][n]) for n in keys0) for k in DUMPS), keys0)
        ok("I", "res_0.h5 に保存量・乱流・組成がある", all(n in keys0 for n in ("ro", "roUx", "roUy", "roe", "roK", "roOmega")), None)
        del v0
        aud = {n: dump(AUDIT, n) for n in ("node_line", "state_ro_roU_roe_cp_gamma", "dt_vol", "flags_wall_iso_axis", "D", "Kprev", "Knext")}
        for key in DUMPS:
            for n in ("node_line", "state_ro_roU_roe_cp_gamma", "flags_wall_iso_axis") + (("dt_vol",) if key in ("a", "a2") else ()):
                ok("I", f"{DUMPS[key][0]}: {n} が監査の記録 (上限なしの方向別 dt) とビット一致 (同じ入力・同じ節点)",
                   arrays[key][n].shape == aud[n].shape and np.array_equal(arrays[key][n].view(np.int64), aud[n].view(np.int64)), None)
        Nchk = np.fromfile(HERE / AUDIT / "linedump" / "audit_node.f64").reshape(aud["node_line"].shape[0], NREC)
        ok("I", "監査の記録: 記録の節点番号 = 書き出しの順、面の数が枠に入る、skipDiag なし",
           np.array_equal(Nchk[:, 0], aud["node_line"][:, 0]) and bool(np.all(Nchk[:, 5] == 0)) and bool(np.all(Nchk[:, 4] <= NF)) and bool(np.all(Nchk[:, 22] == 0)), None)
        if not all(c[2] for c in checks):
            raise ValueError("入力のゲートが不合格")

        # ---- R: 残差の不変 ----
        def res_field(run, f):
            with h5py.File(HERE / run / "res_1.h5", "r") as h:
                return h["VALUE"][f][...].astype(np.float64)
        resv = {}
        opt = [f for f in RES_OPT if all(f in arrays[k]["_have"] for k in DUMPS)]
        rec["residual_optional_used"] = opt
        for f in RES_FIELDS + opt:
            fa, fa2 = res_field(DUMPS["a"][0], f), res_field(DUMPS["a2"][0], f)
            fx = {k: res_field(DUMPS[k][0], f) for k in ("b", "c")}
            if not (np.all(np.isfinite(fa)) and np.all(np.isfinite(fa2)) and all(np.all(np.isfinite(v)) for v in fx.values())):
                raise ValueError(f"{f} に非有限がある (出発の状態の残差)")
            ok("I", f"{f}: 4 本で有限", True, None)
            rep_bit = bool(np.array_equal(fa.view(np.int64), fa2.view(np.int64)))
            d_oo = float(np.max(np.abs(fa2 - fa)))
            resv[f] = {"rerun_bit": rep_bit, "rerun": d_oo}
            for k, v in fx.items():
                d_no = float(np.max(np.abs(v - fa)))
                cond = bool(np.array_equal(v.view(np.int64), fa.view(np.int64))) if rep_bit else d_no <= 3 * d_oo
                ok("R", f"{f}: {k.upper()} と A の差 ({'ビット一致を要求' if rep_bit else '≤ 3 × 再実行の差'})", cond, [d_no, d_oo])
                resv[f][f"{k}_vs_a"] = d_no
        # ---- S: 構造の照合 (§6.12): rhs_s0 = 拘束の処理をした float(残差) が全要素で完全に一致 ----
        rows = [(0, "res_ro"), (1, "res_roUx"), (2, "res_roUy"), (3, "res_roUz"), (4, "res_roe")]
        struct = {}
        for key, (run, _) in DUMPS.items():
            a = arrays[key]
            idx = [int(x) for x in a["node_line"][:, 0]]
            fl = a["flags_wall_iso_axis"]
            with h5py.File(HERE / run / "res_1.h5", "r") as h:
                R = {f: h["VALUE"][f][...] for _, f in rows}
            mism = {}
            for row, f in rows:
                raw = R[f][idx]
                if not np.all(np.isfinite(raw)):
                    raise ValueError(f"{run}: {f} の書き出しの節点に非有限")
                exp = raw.astype(np.float32).astype(np.float64)
                if row in (1, 2, 3):
                    exp = np.where(fl[:, 0] == 1, 0.0, exp)       # 壁: 運動量の 3 行
                if row == 2:
                    exp = np.where(fl[:, 2] == 1, 0.0, exp)       # 軸: 半径の運動量の行
                if row == 4:
                    exp = np.where(fl[:, 1] == 1, 0.0, exp)       # 等温壁: エネルギーの行
                mism[row] = int(np.sum(a["rhs_s0"][:, row] != exp))
            struct[run] = mism
            ok("S", f"{run}: rhs_s0 = 拘束の処理をした float(残差) が 605 節点 × 5 行で完全に一致", sum(mism.values()) == 0, mism)
        rec["struct_rhs"] = struct
        # 記録だけ: rhs_s0 の差の ulp の最大と差のある要素の数
        def ulp_stats(x, ref):
            d = np.abs(x - ref)
            u = np.maximum(np.spacing(np.abs(ref).astype(np.float32)).astype(np.float64), float(np.finfo(np.float32).tiny))
            return {"n_diff": int(np.sum(d > 0)), "max_ulp": float(np.max(d / u)), "max_abs": float(np.max(d))}
        resv["rhs_s0_record_only"] = {f"{k}_vs_a": ulp_stats(arrays[k]["rhs_s0"], arrays["a"]["rhs_s0"]) for k in ("a2", "b", "c")}
        rec["residual"] = resv

        # ---- V: 介入の成立 (dt の介入: K は不変、D は時間項 V/Δτ の差だけ、Δτ は縮む向きにだけ変わる) ----
        nodes = [int(x) for x in arrays["a"]["node_line"][:, 0]]
        pos = {n: k for k, n in enumerate(nodes)}
        fl = arrays["a"]["flags_wall_iso_axis"]
        nn = len(nodes)
        def same_or_noise(name, x, sl=slice(None)):
            a, a2, b = (arrays[k][name][:, sl] for k in ("a", "a2", x))
            rbit = bool(np.array_equal(a.view(np.int64), a2.view(np.int64)))
            d_ba, d_rr = float(np.max(np.abs(b - a))), float(np.max(np.abs(a2 - a)))
            return (bool(np.array_equal(b.view(np.int64), a.view(np.int64))) if rbit else d_ba <= 3 * d_rr), {"rerun_bit": rbit, "x_vs_a": d_ba, "rerun": d_rr}
        dtA = arrays["a"]["dt_vol"][:, 0]
        vol = arrays["a"]["dt_vol"][:, 1]
        vinfo = {}
        dec_rows = np.zeros((nn, 5), dtype=bool)                          # 拘束の行 (壁 → 1〜3、軸 → 2、等温壁 → 4)
        dec_rows[fl[:, 0] == 1, 1:4] = True
        dec_rows[fl[:, 2] == 1, 2] = True
        dec_rows[fl[:, 1] == 1, 4] = True
        for x in ("b", "c"):
            info = {}
            for nm in ("Kprev", "Knext"):
                okK, info[nm] = same_or_noise(nm, x)
                ok("V", f"(a) {x.upper()}: {nm} が A と不変 (K は dt に依らない)", okK, info[nm])
            dtX = arrays[x]["dt_vol"][:, 0]
            r = dtA / dtX
            if not (np.all(np.isfinite(r)) and np.all(r > 0)):
                raise ValueError(f"{x}: Δτ の比が非有限か非正")
            tol_r = 4 * np.spacing(np.float64(1.0))
            ok("V", f"(b) {x.upper()}: 全節点で r = Δτ_A/Δτ_{x.upper()} ≥ 1 − 4 ulp (縮む向きだけ)", bool(np.all(r >= 1 - tol_r)), float(np.min(r)))
            if bool(np.all(r == 1.0)):
                raise ValueError(f"{x}: 全節点で Δτ が A と同じ (上限・方向別のキーが効いていない)")
            rf = {n: float(r[pos[n]]) for n in WALL_FIRST}
            ok("V", f"(b) {x.upper()}: 縮流部の第一内部節点 {WALL_FIRST} で r ≥ 10", all(v >= 10 for v in rf.values()), rf)
            info["r"] = {"min": float(np.min(r)), "median": float(np.median(r)), "max": float(np.max(r)), "wall_first": rf}
            # (c) D: 対角以外と拘束の行はビット一致、拘束の無い行の対角は時間項の差で再現
            DA, DX = arrays["a"]["D"].reshape(nn, 5, 5), arrays[x]["D"].reshape(nn, 5, 5)
            offdiag = ~np.eye(5, dtype=bool)
            same_off = all(np.array_equal(DA[k][offdiag | dec_rows[k][:, None]].view(np.int64), DX[k][offdiag | dec_rows[k][:, None]].view(np.int64)) for k in range(nn))
            ok("V", f"(c) {x.upper()}: D の対角以外と拘束の行が A とビット一致", same_off, None)
            tA = (np.float32(vol) / np.maximum(np.float32(dtA), np.float32(1e-30))).astype(np.float64)
            tX = (np.float32(vol) / np.maximum(np.float32(dtX), np.float32(1e-30))).astype(np.float64)
            worst, resolved = 0.0, {}
            for k in range(nn):
                for i in range(5):
                    if dec_rows[k, i]:
                        continue
                    scale = max(abs(DA[k, i, i]), abs(DX[k, i, i]), tA[k], tX[k])
                    tol = 64 * max(float(np.spacing(np.float32(scale))), float(np.finfo(np.float32).tiny))
                    e = abs((DX[k, i, i] - DA[k, i, i]) - (tX[k] - tA[k])) / tol
                    if not np.isfinite(e):
                        raise ValueError(f"{x}: D の比較に非有限 (節点 {nodes[k]} 行 {i})")
                    worst = max(worst, e)
                    if nodes[k] in WALL_FIRST:                    # 拘束の無い行のうち最もよく解像される行
                        resolved[nodes[k]] = max(resolved.get(nodes[k], 0.0), abs(tX[k] - tA[k]) / tol)
            ok("V", f"(c) {x.upper()}: 拘束の無い行の対角で D_X − D_A = V/Δτ_X − V/Δτ_A が 64 ulp 以内", worst <= 1.0, worst)
            ok("V", f"(c') {x.upper()}: 第一内部節点 {WALL_FIRST} のそれぞれで、拘束の無い行のどれかの時間項の差が許容の 100 倍超", len(resolved) == len(WALL_FIRST) and all(v > 100 for v in resolved.values()), resolved)
            info["D_max_err_over_tol"] = worst
            info["resolved_over_tol"] = resolved
            vinfo[x] = info
        # C は B より縮む (point の dt は上限 50 の dt 以下)
        rB, rC = dtA / arrays["b"]["dt_vol"][:, 0], dtA / arrays["c"]["dt_vol"][:, 0]
        ok("V", "(b) C: 全節点で r_C ≥ r_B − 4 ulp", bool(np.all(rC >= rB * (1 - 4 * np.spacing(np.float64(1.0))))), float(np.min(rC / rB)))
        rec["intervention"] = vinfo
        i_ok = all(c[2] for c in checks if c[0] == "I")
        r_ok = all(c[2] for c in checks if c[0] in ("R", "S"))
        v_ok = all(c[2] for c in checks if c[0] == "V")
        if not i_ok:
            raise ValueError("入力のゲートに外れがある")
        verdict = "PASS (本試験へ進む)" if (r_ok and v_ok) else "INDETERMINATE (残差の変化・構造の照合の外れ・介入の不成立のどれか: 本試験へ進まない)"
        if verdict.startswith("PASS"):
            rec["evidence_sha256"] = evidence()       # 本試験・本判定がこの記録と照合する (M3)
    except Exception as e:
        verdict = f"INVALID ({type(e).__name__}: {e})"
    rec["checks"] = [{"stage": a, "check": b, "ok": c, "value": d} for a, b, c, d in checks]
    rec["VERDICT"] = verdict
    out.write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str))
    for st in ("I", "R", "S", "V"):
        lst = [c for c in checks if c[0] == st]
        bad = [c for c in lst if not c[2]]
        print(f"{st}: {len(lst)} 項目、外れ {len(bad)}  {[(c[1], c[3]) for c in bad[:3]]}")
    for k in ("residual", "struct_rhs", "intervention"):
        if k in rec:
            print(k, rec[k])
    print("VERDICT:", verdict)
    return 0 if verdict.startswith("PASS") else (1 if verdict.startswith("INVALID") else 2)


if __name__ == "__main__":
    sys.exit(verify() if sys.argv[1:] == ["--verify"] else main())
