#!/usr/bin/env python3
"""LHS の座標の差の A/B (plan time_integration-line-viscous-jacobian-faceh §6.9) の事前のゲート。本試験 (run_0560〜0563) の前に回し、合格したときだけ本試験へ進む。

同じ状態 (run_0183 の res_100000) から 1 step の書き出し 3 本を読む:
  old  = run_0564_lvcgeom_old_dump  (段 ② の FP64、LHS は ST(ccx[o]) − ST(ccx[ic]))
  old2 = run_0566_lvcgeom_old_dump2 (段 ② の再実行、atomicAdd による揺れを測る)
  new  = run_0565_lvcgeom_new_dump  (段 ③ の FP64、LHS は ST(±ge_x))
と、§6.7 の監査の記録 (run_0550_lvcaudit、同じ状態・同じ 605 節点、LHS は旧の方式) を使い、
  I  入力のゲート: バイナリ・出発の場・格子・入力ファイル・設定 (本試験との差は step 数・出力の間隔・extraFields だけ)・起動ログ、書き出しの形・節点・有限性・正値、
     3 本と監査の記録の状態・フラグ・dt がビット一致 (同じ入力)
  R  残差の不変: 1 step の出力の更新前の残差 (res_ro・res_roUx・res_roUy・res_roe・res_roK・res_roOmega、全節点) と書き出しの rhs_s0 (ライン上の節点) が、
     旧の再実行どうしでビット一致なら新もビット一致、そうでなければ max|新 − 旧| ≤ 3 × max|旧 − 旧の再実行|
  V  介入の成立: 監査の記録の面ごとの入力から、新の方式 (double の座標の差を float に 1 回丸める) の dcc・δ・β・κ・スカラーの対角・cfac を float32 で作り、
     (a) ライン面 (薄層) の β・κ が double の参照と相対 ≤ 1e-5、
     (b) 薄層の K の変化 (新 − 旧の書き出し) が、共通関数で組んだ Kv(新) − Kv(旧) と要素ごとに ≤ 64 ulp (Kv の列の大きさ)、
     (c) D の変化 (新 − 旧の書き出し) が、line_audit_helper で同じ順序の float の加算を再現した D(新) − D(旧) と要素ごとに ≤ 64 ulp (D の大きさ)
を判定する。出力は _band_ab/cold_pair/lvcgeom_pregate.json。終了コード: 0 = 合格 (本試験へ)、1 = INVALID (入力・証拠の不備)、2 = 判別不能 (介入の不成立・残差の変化)。
使い方: lvcgeom_pregate.py <line_audit_helper>
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
SHA_OLD = "1b8590e84ece63ae50a05cd2d6934440be5f90cbf7f7723a40cd6776ee7285a1"
SHA_NEW = "129de3f4e7f67aa3a80dbd30d5f5cb75df1582d3974e702d76b61c8998598cec"
DUMPS = {"old": ("run_0564_lvcgeom_old_dump", SHA_OLD), "old2": ("run_0566_lvcgeom_old_dump2", SHA_OLD), "new": ("run_0565_lvcgeom_new_dump", SHA_NEW)}
AUDIT = "run_0550_lvcaudit"
WANT_NODES = [1572, 4113, 7985, 198560, 264263]
RES_FIELDS = ["res_ro", "res_roUx", "res_roUy", "res_roe", "res_roK", "res_roOmega"]
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


def ulp32(x):
    return np.spacing(np.abs(np.asarray(x, dtype=np.float32))).astype(np.float64)


def main():
    helper = sys.argv[1] if len(sys.argv) > 1 else None
    out = HERE / "_band_ab" / "cold_pair" / "lvcgeom_pregate.json"
    rec = {"plan": "time_integration-line-viscous-jacobian-faceh §6.9 (事前のゲート)"}
    checks = []                 # (段, 名前, ok, 値)

    def ok(stage, name, cond, val=None):
        checks.append((stage, name, bool(cond), val))
    verdict = None
    try:
        # ---- I: 入力 ----
        ok("I", "line_audit_helper の引数", helper is not None and Path(helper).is_file(), helper)
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
            bad = {k: c.get(k) for k, v in EXPECT.items() if (k in c if v is None else c.get(k) != v)}
            ok("I", f"{run}: 設定 (本試験との差は step 数・出力の間隔・extraFields だけ)", bad == {}, bad)
            log = (d / "forge_run.log").read_text(errors="replace")
            ok("I", f"{run}: 起動ログ (値 3・マスクの表示なし・LAYOUT2・面エンタルピーの切替なし・書き出しあり)",
               "'lineViscCoupling' in 'time.deltaT': 3" in log and "FORGE_LVC_TERMS=" not in log and "Thomas の配列の並び: LAYOUT2" in log
               and "FORGE_DIAG_FACE_H_DOUBLE" not in log and "[lineDump] factor の直前を書いた" in log, None)
            ok("I", f"{run}: RUN_RC = 0", (d / "RUN_RC").read_text().strip() == "0", None)
            ok("I", f"{run}: extraFields の無視の警告なし", "は確保されていない変数なので無視する" not in log, None)
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
        for n in ("node_line", "state_ro_roU_roe_cp_gamma", "dt_vol", "flags_wall_iso_axis"):
            ok("I", f"{n}: 旧・旧の再実行・新の 3 本でビット一致 (同じ入力)",
               all(arrays[k][n].shape == arrays["old"][n].shape and np.array_equal(arrays[k][n].view(np.int64), arrays["old"][n].view(np.int64)) for k in DUMPS), None)
        # 場全体の入力 (保存量・乱流・組成・物性・壁距離): 出力 step 0 の全データセットが 3 本でビット一致
        def h5vals(run):
            with h5py.File(HERE / run / "res_0.h5", "r") as h:
                return {k: h["VALUE"][k][...] for k in h["VALUE"].keys()}
        v0 = {k: h5vals(DUMPS[k][0]) for k in DUMPS}
        keys0 = sorted(v0["old"].keys())
        ok("I", "res_0.h5 の全データセット (保存量・k・ω・組成・μ・壁距離など) が 3 本でビット一致",
           len(keys0) > 0 and all(sorted(v0[k].keys()) == keys0 and all(np.array_equal(v0[k][n], v0["old"][n]) for n in keys0) for k in DUMPS), keys0)
        ok("I", "res_0.h5 に保存量・乱流・組成がある", all(n in keys0 for n in ("ro", "roUx", "roUy", "roe", "roK", "roOmega")), None)
        del v0
        aud = {n: dump(AUDIT, n) for n in ("node_line", "state_ro_roU_roe_cp_gamma", "dt_vol", "flags_wall_iso_axis", "D", "Kprev", "Knext")}
        for key in DUMPS:
            for n in ("node_line", "state_ro_roU_roe_cp_gamma", "dt_vol", "flags_wall_iso_axis"):
                ok("I", f"{DUMPS[key][0]}: {n} が監査の記録とビット一致 (同じ入力・同じ節点)",
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
            o, o2, n = (res_field(DUMPS[k][0], f) for k in ("old", "old2", "new"))
            if not (np.all(np.isfinite(o)) and np.all(np.isfinite(o2)) and np.all(np.isfinite(n))):
                raise ValueError(f"{f} に非有限がある (出発の状態の残差)")   # 2026-10-10 plan-6 レビュー M1 の同型の修正
            ok("I", f"{f}: 3 本で有限", True, None)
            rep_bit = bool(np.array_equal(o.view(np.int64), o2.view(np.int64)))
            d_no, d_oo = float(np.max(np.abs(n - o))), float(np.max(np.abs(o2 - o)))
            cond = bool(np.array_equal(n.view(np.int64), o.view(np.int64))) if rep_bit else d_no <= 3 * d_oo
            ok("R", f"{f}: 新と旧の差 ({'ビット一致を要求' if rep_bit else '≤ 3 × 再実行の差'})", cond, [d_no, d_oo])
            resv[f] = {"rerun_bit": rep_bit, "new_vs_old": d_no, "rerun": d_oo}
        o, o2, n = arrays["old"]["rhs_s0"], arrays["old2"]["rhs_s0"], arrays["new"]["rhs_s0"]
        rep_bit = bool(np.array_equal(o.view(np.int64), o2.view(np.int64)))
        d_no, d_oo = float(np.max(np.abs(n - o))), float(np.max(np.abs(o2 - o)))
        ok("R", f"rhs_s0: 新と旧の差 ({'ビット一致' if rep_bit else '≤ 3 × 再実行の差'})", bool(np.array_equal(n.view(np.int64), o.view(np.int64))) if rep_bit else d_no <= 3 * d_oo, [d_no, d_oo])
        resv["rhs_s0"] = {"rerun_bit": rep_bit, "new_vs_old": d_no, "rerun": d_oo}
        rec["residual"] = resv

        # ---- V: 介入の成立 (監査の記録から新の係数を作る) ----
        ad = HERE / AUDIT / "linedump"
        nn = aud["node_line"].shape[0]
        F = np.fromfile(ad / "audit_face.f64").reshape(nn, NF, FREC)
        N = np.fromfile(ad / "audit_node.f64").reshape(nn, NREC)
        stA = aud["state_ro_roU_roe_cp_gamma"]
        def coef(fr, nd, cp_i, dx, dy, dz):
            """製品の式 (timeIntegration_d.cu の粘性の対角・薄層・枝 3) を float32 で。dx,dy,dz は ST の座標の差"""
            fa = f32(fr[18])
            dcc = max(np.sqrt(dx * dx + dy * dy + dz * dz), f32(1e-30))
            dds = max(abs(dx * f32(fr[8]) + dy * f32(fr[9]) + dz * f32(fr[10])), f32(1e-30))
            delta = max(dcc * fa * fa / dds, f32(1e-30))
            c = {25: dcc, 26: dds, 27: delta, 28: f32(2) * f32(fr[29]) * delta / dcc}
            if int(fr[45]) == 1:
                c[41] = max(f32(fr[39]), f32(0)) * delta / dcc
                c[42] = max(f32(fr[40]), f32(0)) * delta / dcc
            elif int(fr[45]) == 3:
                c[47] = (f32(fr[46]) * delta / dcc) * (f32(nd[15]) / max(f32(cp_i), f32(1e-30))) / f32(nd[9])
            return c
        Fn = F.copy()
        beta_err, old_ulp = [], 0.0
        for k in range(nn):
            nd = N[k]
            for s in range(int(min(nd[4], NF))):
                fr = Fn[k, s]
                if fr[4] == 0 or int(fr[45]) == 0:
                    continue
                # 旧の式 (座標を float にしてから差) を再計算し、記録された GPU の値と合うこと (この再計算の式の検証)
                co = coef(fr, nd, stA[k, 5], f32(fr[15]) - f32(fr[12]), f32(fr[16]) - f32(fr[13]), f32(fr[17]) - f32(fr[14]))
                for q, v in co.items():
                    ru = abs(float(v) - fr[q]) / max(float(np.spacing(np.float32(abs(fr[q])))), float(np.finfo(np.float32).tiny))
                    if not np.isfinite(ru):
                        raise ValueError(f"旧の式の再計算に非有限 (節点 {k} 面 {s})")
                    old_ulp = max(old_ulp, ru)
                # 新の式 (double の差 other − ic = ±ge を float に 1 回丸める)
                e = np.array([fr[15] - fr[12], fr[16] - fr[13], fr[17] - fr[14]])
                dx, dy, dz = f32(e[0]), f32(e[1]), f32(e[2])
                cn = coef(fr, nd, stA[k, 5], dx, dy, dz)
                fr[22], fr[23], fr[24] = dx, dy, dz
                for q, v in cn.items():
                    fr[q] = float(v)
                if int(fr[45]) == 1:
                    dod = fr[11] ** 2 / abs(float(e @ fr[8:11]))                         # double の参照の δ/dcc
                    beta_err.append(max(abs(fr[41] / (max(fr[39], 0) * dod) - 1), abs(fr[42] / (max(fr[40], 0) * dod) - 1)))
        rec["recompute_old_formula_max_ulp"] = old_ulp
        ok("V", "(0) 旧の式の再計算が記録された GPU の dcc・δ・粘性の対角・β・κ・cfac と ≤ 8 ulp (式の検証)", old_ulp <= 8, old_ulp)
        rec["intervention_new_coef_max_rel_vs_double"] = max(beta_err) if beta_err else None
        if not all(np.isfinite(x) for x in beta_err):
            raise ValueError("β・κ の参照との比に非有限")
        ok("V", "(a) 新の方式の薄層の β・κ が double の参照と ≤ 1e-5", bool(beta_err) and max(beta_err) <= 1e-5, max(beta_err) if beta_err else None)
        # helper で旧・新の面の加算を再現 (時間項の後から連ねる)
        with tempfile.TemporaryDirectory() as td:
            outs = {}
            for tag, arr in (("old", F), ("new", Fn)):
                dd = Path(td) / tag; dd.mkdir()
                arr.tofile(dd / "audit_face.f64"); N.tofile(dd / "audit_node.f64")
                r = subprocess.run([helper, str(dd), str(nn), "1"], capture_output=True, text=True)
                if r.returncode != 0:
                    raise ValueError(f"line_audit_helper ({tag}) が失敗: {r.stderr}")
                outs[tag] = np.fromfile(dd / "audit_host_face.f64").reshape(nn, NF, HF)
        # (b) 薄層の K の変化
        Ko, Kn = {"Kprev": arrays["old"]["Kprev"], "Knext": arrays["old"]["Knext"]}, {"Kprev": arrays["new"]["Kprev"], "Knext": arrays["new"]["Knext"]}
        worst_k = 0.0; nk = 0
        for k in range(nn):
            nd = N[k]; rd = int(nd[18])
            for s in range(int(min(nd[4], NF))):
                fr = F[k, s]
                if fr[5] == 0 or int(fr[45]) != 1:
                    continue
                nm = "Kprev" if int(fr[6]) == 0 else "Knext"
                dK_dump = (Kn[nm][k] - Ko[nm][k]).reshape(5, 5)
                kv_new, kv_old = outs["new"][k, s, 125:150].reshape(5, 5), outs["old"][k, s, 125:150].reshape(5, 5)
                pred = np.where(np.array([[(rd >> i) & 1 for _ in range(5)] for i in range(5)]) == 1, 0.0, kv_new - kv_old)
                scale = np.broadcast_to(np.maximum(np.max(np.abs(kv_new), axis=0), np.max(np.abs(kv_old), axis=0)), (5, 5))
                u = np.maximum(ulp32(scale), np.finfo(np.float32).tiny)
                rk = float(np.max(np.abs(dK_dump - pred) / u))
                if not np.isfinite(rk):
                    raise ValueError(f"K の比較に非有限 (節点 {k} 面 {s})")
                worst_k = max(worst_k, rk); nk += 1
        rec["intervention_K"] = {"faces": nk, "max_ulp": worst_k}
        ok("V", "(b) 薄層の K の変化が Kv(新) − Kv(旧) と ≤ 64 ulp", nk > 0 and worst_k <= 64, [nk, worst_k])
        # (c) D の変化 (拘束の前の面の和の再現の差を、書き出しの D の差と比べる。拘束の行は変化 0 を要求)
        Do, Dn = arrays["old"]["D"], arrays["new"]["D"]
        worst_d = 0.0
        for k in range(nn):
            nd = N[k]; nfc = int(min(nd[4], NF))
            if nfc == 0:
                continue
            dpred = (outs["new"][k, nfc - 1, 25:50] - outs["old"][k, nfc - 1, 25:50]).reshape(5, 5)
            rows = [i for i in range(5) if (int(nd[18]) >> i) & 1]
            if int(nd[21]) == 1:
                rows.append(2)
            for i in set(rows):
                dpred[i, :] = 0.0
            ddump = (Dn[k] - Do[k]).reshape(5, 5)
            hs = np.maximum(np.max(outs["old"][k, :nfc, 50:75], axis=0), np.max(outs["new"][k, :nfc, 50:75], axis=0))
            scale = np.maximum(np.maximum(np.abs(Do[k]), np.abs(Dn[k])), hs).reshape(5, 5)   # 面の途中の和の大きさも含める
            u = np.maximum(ulp32(scale), np.finfo(np.float32).tiny)
            rd_ = float(np.max(np.abs(ddump - dpred) / u))
            if not np.isfinite(rd_):
                raise ValueError(f"D の比較に非有限 (節点 {k})")
            worst_d = max(worst_d, rd_)
        rec["intervention_D"] = {"max_ulp": worst_d}
        ok("V", "(c) D の変化が係数の変更の再現と ≤ 64 ulp", worst_d <= 64, worst_d)
        # 記録: 旧の書き出しと監査の記録の D・K が一致するか (別のバイナリ)
        rec["old_dump_vs_audit_bit"] = {n: bool(np.array_equal(arrays["old"][n].view(np.int64), aud[n].view(np.int64))) for n in ("D", "Kprev", "Knext")}
        rec["near_wall_K_change_max_rel"] = float(np.max(np.abs(Kn["Kprev"] - Ko["Kprev"]) / np.maximum(np.abs(Ko["Kprev"]), 1e-300)))
        if not all(c[2] for c in checks if c[0] == "I"):
            raise ValueError("入力のゲートに外れがある")
        r_ok = all(c[2] for c in checks if c[0] == "R")
        v_ok = all(c[2] for c in checks if c[0] == "V")
        verdict = "PASS (本試験へ進む)" if (r_ok and v_ok) else "INDETERMINATE (残差の変化か介入の不成立: 本試験へ進まない)"
    except Exception as e:
        verdict = f"INVALID ({type(e).__name__}: {e})"
    rec["checks"] = [{"stage": a, "check": b, "ok": c, "value": d} for a, b, c, d in checks]
    rec["VERDICT"] = verdict
    out.write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str))
    for st in ("I", "R", "V"):
        lst = [c for c in checks if c[0] == st]
        bad = [c for c in lst if not c[2]]
        print(f"{st}: {len(lst)} 項目、外れ {len(bad)}  {[(c[1], c[3]) for c in bad[:3]]}")
    for k in ("residual", "intervention_new_coef_max_rel_vs_double", "intervention_K", "intervention_D", "old_dump_vs_audit_bit"):
        if k in rec:
            print(k, rec[k])
    print("VERDICT:", verdict)
    return 0 if verdict.startswith("PASS") else (1 if verdict.startswith("INVALID") else 2)


if __name__ == "__main__":
    sys.exit(main())
