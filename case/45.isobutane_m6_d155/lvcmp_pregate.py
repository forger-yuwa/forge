#!/usr/bin/env python3
"""上限 50 での全部入り・粘性入り・生産の比較 (plan time_integration-line-viscous-jacobian-dt §6.3) の事前のゲート。本試験 (run_0553〜0555) の前に回し、合格したときだけ本試験へ進む。

同じ状態 (run_0183 の res_100000) から段 ③ の FP64・方向別 dt 上限 50・キー 5・cfl 4 で 1 step の書き出しを 4 本読む:
  m7  = run_0556_lvcmp_l7_dump  (値 3・FORGE_LVC_TERMS=7、全部入り)
  m7b = run_0559_lvcmp_l7_dump2 (m7 の再実行、atomicAdd による揺れを測る)
  m5  = run_0558_lvcmp_l5_dump  (値 3・FORGE_LVC_TERMS=5、粘性入り = 熱伝導の近傍 K なし)
  m0  = run_0557_lvcmp_l0_dump  (値 0、生産: ライン面の粘性はスカラーの対角だけ)
と、faceh §6.7 の監査の記録 (run_0550_lvcaudit、同じ状態・同じ 605 節点) の面ごとの生の入力を使い、
  I  入力のゲート (faceh §6.12 と同じ項目。設定は 4 本で lineViscCoupling 以外同じ、上限 50。dt_vol は 4 本と dt plan の run_0585 でビット一致)
  R  残差の不変: 出発の状態の残差 (全節点、倍精度) が、m7 の再実行どうしでビット一致なら m5・m0 もビット一致、そうでなければ ≤ 3 × 再実行の差
  S  構造の照合: 各書き出しで rhs_s0 = 拘束の処理をした float(出力 step 1 の残差) が全要素で完全に一致
  V  介入の成立 (期待の接続の全部、端は 0):
     (a) m7 − m5: D と K の行 0〜3 が不変、K の行 4 の変化が独立に計算した熱伝導の近傍 K の項と許容内 (faceh §6.12 と同じ)
     (b) m5 − m0: K の行 0 が不変、K の行 1〜3 の変化が独立に計算した運動量の薄層の項 β P ∂u_j/∂Q_j と許容内 (要素ごとに 1e-5 × |h| + 64 ulp)
     (c') それぞれ、変化が解像される面が 100 面以上
を判定する。出力は _band_ab/cold_pair/lvcmp_pregate.json。終了コード: 0 = 合格、1 = INVALID、2 = 判別不能。使い方: lvcmp_pregate.py [--verify]
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
DUMPS = {"m7": ("run_0556_lvcmp_l7_dump", SHA_NEW), "m7b": ("run_0559_lvcmp_l7_dump2", SHA_NEW),
         "m5": ("run_0558_lvcmp_l5_dump", SHA_NEW), "m0": ("run_0557_lvcmp_l0_dump", SHA_NEW)}
MASK = {"m7": 7, "m7b": 7, "m5": 5, "m0": None}          # None = 値 0 (FORGE_LVC_TERMS を渡さない、LVC_TERMS.txt は "-")
ARM_CFG = {"m7": {}, "m7b": {}, "m5": {}, "m0": {"time/deltaT/lineViscCoupling": None}}
CAP_REF = "run_0585_lvcdt_b_dump"                         # dt plan §6.2 で用量を確かめた上限 50 の書き出し (同じ状態)
AUDIT = "run_0550_lvcaudit"
WANT_NODES = [1572, 4113, 7985, 198560, 264263]
RES_FIELDS = ["res_ro", "res_roUx", "res_roUy", "res_roUz", "res_roe", "res_roK", "res_roOmega"]   # res_roUz は rhs_s0 の行 3 の照合に使う (§6.12)
RES_OPT = ["res_roY0", "res_roY1"]          # 化学種の残差 (3 本ともにあれば R に含める)
EXTRA = RES_FIELDS + RES_OPT + ["volume"]
SAME_FILES = ["bcondConfig.yaml", "probe.yaml", "species_meta.yaml", "wall_design.csv", "wall_physical.csv", "target_axis_M.csv", "wall_repr.json", "MESH_QUALITY.txt"]
EXPECT = {"time/deltaT/cfl": 4.0, "time/deltaT/cfl_pseudo": 4.0, "time/deltaT/implicitRelax": 0.7, "time/deltaT/implicitThermalJacobian": 5,
          "time/deltaT/lineDtDirectional": 1, "time/deltaT/lineImplicit": 1, "time/deltaT/lineViscCoupling": 3, "time/deltaT/blockDPLUR": 1,
          "time/deltaT/detectNaN": 1, "time/timeIntegration": 11, "time/nStepInner": 5, "time/last/nStepOuter": 1,
          "time/deltaT/lineDtDirectionalCap": 50.0, "time/deltaT/implicitSolvePrecision": None, "time/outStepInterval": 1,
          "output/extraFields": EXTRA, "time/deltaT/axisTimestepBeta": None, "time/deltaT/lineViscousDtRelief": None}
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
    ev[f"{CAP_REF}/linedump/dt_vol.f64"] = sha(HERE / CAP_REF / "linedump" / "dt_vol.f64")
    ev["lvcmp_pregate.py"] = sha(Path(__file__).resolve())
    return ev


def verify():
    """--verify: 記録が PASS で、記録したハッシュが今のファイルと同じなら 0、そうでなければ 1"""
    rec = json.loads((HERE / "_band_ab" / "cold_pair" / "lvcmp_pregate.json").read_text())
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
    out = HERE / "_band_ab" / "cold_pair" / "lvcmp_pregate.json"
    rec = {"plan": "time_integration-line-viscous-jacobian-dt §6.3 (事前のゲート)"}
    checks = []

    def ok(stage, name, cond, val=None):
        checks.append((stage, name, bool(cond), val))
    verdict = None
    try:
        # ---- I: 入力 ----
        src_full = sha(HERE / SRC / SRC_RES)
        ok("I", "出発の場の sha256", src_full.startswith(SRC_SHA16), src_full[:16])
        mref = meshsha(HERE / SRC / "nozzle.h5")
        cfgs = {}
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
            cfgs[key] = c
            exp_ = {**EXPECT, **ARM_CFG[key]}
            bad = {k: c.get(k) for k, v in exp_.items() if (k in c if v is None else c.get(k) != v)}
            ok("I", f"{run}: 設定 (腕ごとの期待どおり)", bad == {}, bad)
            log = (d / "forge_run.log").read_text(errors="replace")
            mk = MASK[key]
            shown = (f"診断のマスク FORGE_LVC_TERMS={mk} " in log) if mk == 5 else ("FORGE_LVC_TERMS=" not in log)
            vline = ("'lineViscCoupling' in 'time.deltaT': 3" in log) if mk is not None else ("'lineViscCoupling' in 'time.deltaT': 3" not in log)
            ok("I", f"{run}: 起動ログ (値・マスクの表示・LAYOUT2・面エンタルピーの切替なし・書き出しあり)",
               vline and shown and "Thomas の配列の並び: LAYOUT2" in log and "FORGE_DIAG_FACE_H_DOUBLE" not in log and "[lineDump] factor の直前を書いた" in log, None)
            lt = (d / "LVC_TERMS.txt").read_text().strip() if (d / "LVC_TERMS.txt").is_file() else None
            ok("I", f"{run}: 台本が渡した FORGE_LVC_TERMS = {mk if mk is not None else '-'}", lt == (str(mk) if mk is not None else "-"), lt)
            ok("I", f"{run}: RUN_RC = 0", (d / "RUN_RC").read_text().strip() == "0", None)
            ok("I", f"{run}: extraFields の無視の警告なし", "は確保されていない変数なので無視する" not in log, None)
        strip = lambda c: {k: v for k, v in c.items() if k != "time/deltaT/lineViscCoupling"}
        ok("I", "4 本の設定が lineViscCoupling 以外で完全に同じ", all(strip(cfgs[k]) == strip(cfgs["m7"]) for k in DUMPS), None)
        rec["dump_config"] = cfgs
        arrays = {}
        for key, (run, _) in DUMPS.items():
            d = HERE / run
            a = {n: dump(run, n) for n in ("node_line", "D", "Kprev", "Knext", "state_ro_roU_roe_cp_gamma", "dt_vol", "flags_wall_iso_axis", "rhs_s0")}
            nn = a["node_line"].shape[0]
            shapes = {"node_line": 2, "D": 25, "Kprev": 25, "Knext": 25, "state_ro_roU_roe_cp_gamma": 7, "dt_vol": 2, "flags_wall_iso_axis": 3, "rhs_s0": 5}
            ok("I", f"{run}: 配列の形", nn > 0 and all(a[n].shape == (nn, w) for n, w in shapes.items()), {n: a[n].shape for n in shapes})
            nodes_ = [int(x) for x in a["node_line"][:, 0]]
            ok("I", f"{run}: 要求した 5 節点・5 本のライン・重複なし", all(w in nodes_ for w in WANT_NODES) and len(set(a["node_line"][:, 1])) == 5 and len(set(nodes_)) == nn, None)
            ok("I", f"{run}: 有限", all(bool(np.all(np.isfinite(a[n]))) for n in shapes), None)
            st_ = a["state_ro_roU_roe_cp_gamma"]
            ok("I", f"{run}: ρ・c_p・dt・体積 > 0、γ > 1", bool(np.all(st_[:, 0] > 0) and np.all(st_[:, 5] > 0) and np.all(st_[:, 6] > 1) and np.all(a["dt_vol"] > 0)), None)
            with h5py.File(d / "res_1.h5", "r") as h:
                have = set(h["VALUE"].keys())
            ok("I", f"{run}: res_1.h5 に残差の場がそろう", all(f in have for f in RES_FIELDS), sorted(set(RES_FIELDS) - have))
            arrays[key] = a
            arrays[key]["_have"] = have
        for n in ("node_line", "state_ro_roU_roe_cp_gamma", "flags_wall_iso_axis", "dt_vol"):
            ok("I", f"{n}: 4 本でビット一致 (同じ入力・同じ上限 50 の dt)",
               all(arrays[k][n].shape == arrays["m7"][n].shape and np.array_equal(arrays[k][n].view(np.int64), arrays["m7"][n].view(np.int64)) for k in DUMPS), None)
        capref = dump(CAP_REF, "dt_vol")
        ok("I", f"dt_vol が {CAP_REF} (dt plan §6.2 で用量を確かめた上限 50) とビット一致", capref.shape == arrays["m7"]["dt_vol"].shape and np.array_equal(capref.view(np.int64), arrays["m7"]["dt_vol"].view(np.int64)), None)
        def h5vals(run):
            with h5py.File(HERE / run / "res_0.h5", "r") as h:
                return {k: h["VALUE"][k][...] for k in h["VALUE"].keys()}
        v0 = {k: h5vals(DUMPS[k][0]) for k in DUMPS}
        keys0 = sorted(v0["m7"].keys())
        ok("I", "res_0.h5 の全データセットが 4 本でビット一致",
           len(keys0) > 0 and all(sorted(v0[k].keys()) == keys0 and all(np.array_equal(v0[k][n], v0["m7"][n]) for n in keys0) for k in DUMPS), keys0)
        del v0
        aud = {n: dump(AUDIT, n) for n in ("node_line", "state_ro_roU_roe_cp_gamma", "flags_wall_iso_axis")}
        for key in DUMPS:
            for n in ("node_line", "state_ro_roU_roe_cp_gamma", "flags_wall_iso_axis"):
                ok("I", f"{DUMPS[key][0]}: {n} が監査の記録とビット一致", arrays[key][n].shape == aud[n].shape and np.array_equal(arrays[key][n].view(np.int64), aud[n].view(np.int64)), None)
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
        for f in RES_FIELDS + opt:
            fa, fa2 = res_field(DUMPS["m7"][0], f), res_field(DUMPS["m7b"][0], f)
            fx = {k: res_field(DUMPS[k][0], f) for k in ("m5", "m0")}
            if not (np.all(np.isfinite(fa)) and np.all(np.isfinite(fa2)) and all(np.all(np.isfinite(v)) for v in fx.values())):
                raise ValueError(f"{f} に非有限がある (出発の状態の残差)")
            rbit = bool(np.array_equal(fa.view(np.int64), fa2.view(np.int64)))
            d_oo = float(np.max(np.abs(fa2 - fa)))
            resv[f] = {"rerun_bit": rbit, "rerun": d_oo}
            for k, v in fx.items():
                d_no = float(np.max(np.abs(v - fa)))
                ok("R", f"{f}: {k} と m7 の差 ({'ビット一致を要求' if rbit else '≤ 3 × 再実行の差'})",
                   bool(np.array_equal(v.view(np.int64), fa.view(np.int64))) if rbit else d_no <= 3 * d_oo, [d_no, d_oo])
                resv[f][f"{k}_vs_m7"] = d_no
        rec["residual"] = resv
        # ---- S: 構造の照合 ----
        rows = [(0, "res_ro"), (1, "res_roUx"), (2, "res_roUy"), (3, "res_roUz"), (4, "res_roe")]
        struct = {}
        for key, (run, _) in DUMPS.items():
            a = arrays[key]
            idx = [int(x) for x in a["node_line"][:, 0]]
            fl_ = a["flags_wall_iso_axis"]
            with h5py.File(HERE / run / "res_1.h5", "r") as h:
                R = {f: h["VALUE"][f][...] for _, f in rows}
            mism = {}
            for row, f in rows:
                raw = R[f][idx]
                if not np.all(np.isfinite(raw)):
                    raise ValueError(f"{run}: {f} の書き出しの節点に非有限")
                exp = raw.astype(np.float32).astype(np.float64)
                if row in (1, 2, 3):
                    exp = np.where(fl_[:, 0] == 1, 0.0, exp)
                if row == 2:
                    exp = np.where(fl_[:, 2] == 1, 0.0, exp)
                if row == 4:
                    exp = np.where(fl_[:, 1] == 1, 0.0, exp)
                mism[row] = int(np.sum(a["rhs_s0"][:, row] != exp))
            struct[run] = mism
            ok("S", f"{run}: rhs_s0 = 拘束の処理をした float(残差) が全要素で完全に一致", sum(mism.values()) == 0, mism)
        rec["struct_rhs"] = struct

        # ---- V: 介入の成立 ----
        ad = HERE / AUDIT / "linedump"
        nn = aud["node_line"].shape[0]
        F = np.fromfile(ad / "audit_face.f64").reshape(nn, NF, FREC)
        N = np.fromfile(ad / "audit_node.f64").reshape(nn, NREC)
        st = arrays["m7"]["state_ro_roU_roe_cp_gamma"]
        nodes = [int(x) for x in arrays["m7"]["node_line"][:, 0]]
        pos = {n: k for k, n in enumerate(nodes)}
        flags = arrays["m7"]["flags_wall_iso_axis"]
        nl = arrays["m7"]["node_line"]

        def same_or_noise(name, x, ref, sl):
            a, a2, b = arrays[ref][name][:, sl], arrays["m7b"][name][:, sl], arrays[x][name][:, sl]
            a7, a7b = arrays["m7"][name][:, sl], arrays["m7b"][name][:, sl]
            rbit = bool(np.array_equal(a7.view(np.int64), a7b.view(np.int64)))
            d_ba, d_rr = float(np.max(np.abs(b - a))), float(np.max(np.abs(a7b - a7)))
            return (bool(np.array_equal(b.view(np.int64), a.view(np.int64))) if rbit else d_ba <= 3 * d_rr), {"rerun_bit": rbit, "diff": d_ba, "rerun": d_rr}
        same = {}
        okD, same["D_m5_vs_m7"] = same_or_noise("D", "m5", "m7", slice(0, 25))
        ok("V", "(a) D が m7 と m5 で不変 (熱伝導の D は常に入る)", okD, same["D_m5_vs_m7"])
        for nm in ("Kprev", "Knext"):
            okK, same[f"{nm}_rows0to3_m5_vs_m7"] = same_or_noise(nm, "m5", "m7", slice(0, 20))
            ok("V", f"(a) {nm} の行 0〜3 が m7 と m5 で不変", okK, same[f"{nm}_rows0to3_m5_vs_m7"])
            okK0, same[f"{nm}_row0_m0_vs_m5"] = same_or_noise(nm, "m0", "m5", slice(0, 5))
            ok("V", f"(b) {nm} の行 0 が m5 と m0 で不変 (連続の行に粘性の流束は無い)", okK0, same[f"{nm}_row0_m0_vs_m5"])
        # (b2) 等温壁の節点の D の行 4 (codex plan-2 m4): 値 ≥ 2 は拘束の行 [−ρE/ρ, 0, 0, 0, 1] (timeIntegration_d.cu 1263 行付近)、値 0・キー 5 は [0, 0, 0, 0, 1]
        isoN = [k for k in range(arrays["m7"]["node_line"].shape[0]) if int(arrays["m7"]["flags_wall_iso_axis"][k, 1]) == 1]
        stv = arrays["m7"]["state_ro_roU_roe_cp_gamma"]
        worst_iso, bad_iso = 0.0, []
        for k in isoN:
            rho = max(np.float32(stv[k, 0]), np.float32(1e-30)); roe = np.float32(stv[k, 4])
            e40 = float(-roe / rho)
            for key, exp in (("m7", [e40, 0, 0, 0, 1]), ("m5", [e40, 0, 0, 0, 1]), ("m0", [0, 0, 0, 0, 1])):
                row4 = arrays[key]["D"][k].reshape(5, 5)[4]
                tol = 4 * float(np.spacing(np.float32(abs(e40)))) if key != "m0" else 0.0
                err = abs(row4[0] - exp[0])
                if not (err <= tol and all(row4[j] == exp[j] for j in range(1, 5))):
                    bad_iso.append((nodes_k := int(arrays["m7"]["node_line"][k, 0]), key, row4.tolist()))
                worst_iso = max(worst_iso, err / tol if tol > 0 else err)
        ok("V", f"(b2) 等温壁の節点 ({len(isoN)} 個) の D の行 4: 値 3 は [−ρE/ρ, 0, 0, 0, 1] (4 ulp)、値 0 は [0, 0, 0, 0, 1]", len(isoN) > 0 and not bad_iso, bad_iso[:3] or worst_iso)
        rec["same"] = same
        exp_conn = []
        for k in range(nn):
            same_prev = k > 0 and nl[k - 1, 1] == nl[k, 1]
            same_next = k + 1 < nn and nl[k + 1, 1] == nl[k, 1]
            exp_conn.append((k, 0, nodes[k - 1] if same_prev else -1))
            exp_conn.append((k, 1, nodes[k + 1] if same_next else -1))
        conn_bad = []
        res = {"heat": [0.0, 0, 0], "mom": [0.0, 0, 0]}     # [比の最大, 比べた面, 解像された面]
        for k, side, nb in exp_conn:
            nd = N[k]
            if int(nd[0]) != nodes[k] or int(nd[1 + side]) != nb:
                conn_bad.append((nodes[k], side, "記録の line_prev/next が並びと違う")); continue
            wall, iso, ax = int(flags[k, 0]), int(flags[k, 1]), int(flags[k, 2])
            exp_dec = (14 if wall == 1 else 0) | (16 if iso == 1 else 0) | (4 if ax == 1 else 0)
            if int(nd[18]) != exp_dec:
                conn_bad.append((nodes[k], side, "拘束の行のビットがフラグと違う")); continue
            nm = "Kprev" if side == 0 else "Knext"
            K7, K5, K0 = (arrays[x][nm][k].reshape(5, 5) for x in ("m7", "m5", "m0"))
            if nb < 0:
                if not (np.all(K7 == 0) and np.all(K5 == 0) and np.all(K0 == 0)):
                    conn_bad.append((nodes[k], side, "ラインの端なのに K が 0 でない"))
                continue
            fs = [s_ for s_ in range(int(min(nd[4], NF))) if F[k, s_, 5] != 0 and int(F[k, s_, 3]) == nb]
            if len(fs) != 1:
                conn_bad.append((nodes[k], side, "期待の接続にライン面がちょうど 1 つでない")); continue
            fr = F[k, fs[0]]
            kb = pos.get(nb)
            back = [s_ for s_ in range(int(min(N[kb, 4], NF))) if F[kb, s_, 5] != 0 and int(F[kb, s_, 3]) == nodes[k]] if kb is not None else []
            if not (int(fr[6]) == side and int(fr[45]) == 1 and {int(fr[1]), int(fr[2])} == {nodes[k], nb} and kb is not None
                    and len(back) == 1 and int(F[kb, back[0], 0]) == int(fr[0]) and int(F[kb, back[0], 6]) == 1 - side):
                conn_bad.append((nodes[k], side, "面の向き・枝・両端・往復の面が合わない")); continue
            if bool(fr[43] != 0) != (int(flags[kb, 0]) == 1) or bool(fr[44] != 0) != (int(flags[kb, 1]) == 1):
                conn_bad.append((nodes[k], side, "隣の速度・温度固定のフラグが書き出しのフラグと違う")); continue
            f0 = fr[7]
            ic0_is_k = int(fr[1]) == nodes[k]
            e = np.array([fr[15] - fr[12], fr[16] - fr[13], fr[17] - fr[14]])
            dod = fr[11] ** 2 / abs(float(e @ fr[8:11]))
            rj, uj, vj, wj, rEj, gj = fr[52], fr[53], fr[54], fr[55], fr[56], fr[51]
            # 熱伝導の近傍 K (m7 − m5、行 4)
            kf = f0 * fr[186] + (1 - f0) * fr[187] + (f0 * fr[184] + (1 - f0) * fr[185]) * (f0 * fr[182] + (1 - f0) * fr[183]) / fr[50]
            kap = max(kf, 0.0) * dod
            cpj = max(st[kb, 5], 1e-30)
            q2 = uj * uj + vj * vj + wj * wj
            ej = rEj / rj - 0.5 * q2
            hh = kap * (gj / cpj) / rj * np.array([-(ej - 0.5 * q2), -uj, -vj, -wj, 1.0])
            if iso == 1 or int(flags[kb, 1]) == 1:
                hh = np.zeros(5)
            # 運動量の薄層 (m5 − m0、行 1〜3): K[1+a][0] = −β Σ_b P_ab u_j,b / ρ_j、K[1+a][1+b] = β P_ab / ρ_j
            mlam0, mlam1 = (fr[188], fr[189]) if ic0_is_k else (fr[189], fr[188])
            muf = f0 * (mlam0 + fr[182]) + (1 - f0) * (mlam1 + fr[183])
            beta = max(muf, 0.0) * dod
            nvec = np.array(fr[8:11]) / fr[11]
            P = np.eye(3) + np.outer(nvec, nvec) / 3.0
            ujv = np.array([uj, vj, wj])
            hm = np.zeros((3, 5)); hs = np.zeros((3, 5))
            for a_ in range(3):
                hm[a_, 0] = -beta * float(P[a_] @ ujv) / rj
                hs[a_, 0] = beta * float(np.abs(P[a_]) @ np.abs(ujv)) / rj
                hm[a_, 1:4] = beta * P[a_] / rj
                hs[a_, 1:4] = np.abs(hm[a_, 1:4])
            for a_ in range(3):
                if (wall == 1) or (ax == 1 and a_ == 1) or int(flags[kb, 0]) == 1:
                    hm[a_] = 0.0; hs[a_] = 0.0
            for tag, dK, h, hsc in (("heat", K7[4] - K5[4], hh, np.abs(hh)), ("mom", (K5[1:4] - K0[1:4]).ravel(), hm.ravel(), hs.ravel())):
                scale = np.maximum.reduce([np.abs(dK), hsc, np.abs(h)])
                if tag == "heat":
                    scale = np.maximum.reduce([scale, np.abs(K7[4]), np.abs(K5[4])])
                else:
                    scale = np.maximum.reduce([scale, np.abs(K5[1:4]).ravel(), np.abs(K0[1:4]).ravel()])
                floor = 64 * np.maximum(ulp32(scale), np.finfo(np.float32).tiny)
                tol = 1e-5 * hsc + floor
                ratio = np.abs(dK - h) / tol
                if not (np.isfinite(kap) and np.isfinite(beta) and np.all(np.isfinite(h)) and np.all(np.isfinite(ratio))):
                    raise ValueError(f"節点 {nodes[k]} 向き {side}: {tag} の比に非有限 (入力の不備)")
                res[tag][0] = max(res[tag][0], float(np.max(ratio))); res[tag][1] += 1
                if float(np.max(np.abs(h))) > 1e3 * float(np.max(floor)):
                    res[tag][2] += 1
        rec["intervention"] = {"heat_m7_minus_m5": res["heat"], "mom_m5_minus_m0": res["mom"], "connection_problems": conn_bad[:10], "expected_connections": len(exp_conn)}
        if conn_bad:
            raise ValueError(f"期待の接続と監査の記録が合わない ({len(conn_bad)} 件、最初 {conn_bad[0]})")
        ok("V", "(a) K の行 4 の変化 (m7 − m5) が独立に計算した熱伝導の項と許容内 (要素ごとに 1e-5 × |h| + 64 ulp)、端は 0", res["heat"][1] > 0 and res["heat"][0] <= 1.0, res["heat"])
        ok("V", "(a') 熱伝導の項が解像される面が 100 面以上", res["heat"][2] >= 100, res["heat"][2])
        ok("V", "(b) K の行 1〜3 の変化 (m5 − m0) が独立に計算した運動量の薄層の項と許容内 (要素ごとに 1e-5 × 尺度 + 64 ulp)、端は 0", res["mom"][1] > 0 and res["mom"][0] <= 1.0, res["mom"])
        ok("V", "(b') 運動量の項が解像される面が 100 面以上", res["mom"][2] >= 100, res["mom"][2])
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
    for k in ("residual", "struct_rhs", "same", "intervention"):
        if k in rec:
            print(k, rec[k])
    print("VERDICT:", verdict)
    return 0 if verdict.startswith("PASS") else (1 if verdict.startswith("INVALID") else 2)


if __name__ == "__main__":
    sys.exit(verify() if sys.argv[1:] == ["--verify"] else main())
