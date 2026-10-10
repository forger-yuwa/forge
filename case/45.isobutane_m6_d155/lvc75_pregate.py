#!/usr/bin/env python3
"""値 3 のマスク 7/5 の比較 (plan time_integration-line-viscous-jacobian-faceh §6.11) の事前のゲート。本試験 (run_0570〜0573) の前に回し、合格したときだけ本試験へ進む。

同じ状態 (run_0183 の res_100000) から段 ③ の FP64 で 1 step の書き出しを 3 本読む:
  m7  = run_0574_lvc75_m7_dump  (FORGE_LVC_TERMS=7、既定の全部入り)
  m5  = run_0575_lvc75_m5_dump  (FORGE_LVC_TERMS=5、熱伝導の近傍 K を外す)
  m7b = run_0576_lvc75_m7_dump2 (7 の再実行、atomicAdd による揺れを測る)
と、§6.7 の監査の記録 (run_0550_lvcaudit、同じ状態・同じ 605 節点) の面ごとの生の入力を使い、
  I  入力のゲート (§6.9 の lvcgeom_pregate.py と同じ項目に、マスクの表示と台本が渡した値を足したもの)
  R  残差の不変: 出発の状態の残差 (全節点、倍精度) が、7 の再実行どうしでビット一致なら 5 もビット一致、そうでなければ ≤ 3 × 再実行の差
  S  構造の照合 (2026-10-10 改訂、§6.12): 各書き出しで rhs_s0 (605 節点 × 5 行) = 拘束の処理をした float(出力 step 1 の残差) が全要素で完全に一致。
     rhs_s0 の 5 − 7・再実行 − 7 の差 (ulp の最大・差のある要素の数) は記録だけ (1 回目の規則 max|5 − 7| ≤ 3 × 再実行は、丸めの反転の位置で決まるので外した)
  V  介入の成立: (a) D と (b) K の行 0〜3 が 7 と 5 で不変 (R と同じ規則)、
     (c) K の行 4 の変化 (7 − 5) が、生の double の入力から独立に計算した熱伝導の近傍 K の項 κ ∂T_j/∂Q_j と許容内 (要素ごとに 1e-5 × |h| + 64 ulp)、
     (c') 変化が解像される面が 100 面以上
を判定する。出力は _band_ab/cold_pair/lvc75_pregate.json。終了コード: 0 = 合格 (本試験へ)、1 = INVALID (入力・証拠の不備)、2 = 判別不能 (残差の変化・介入の不成立)。
使い方: lvc75_pregate.py
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
DUMPS = {"m7": ("run_0577_lvc75_m7_dump_r2", SHA_NEW), "m7b": ("run_0579_lvc75_m7_dump2_r2", SHA_NEW), "m5": ("run_0578_lvc75_m5_dump_r2", SHA_NEW)}   # 2 回目の登録 (§6.12)
MASK = {"m7": 7, "m7b": 7, "m5": 5}
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
    ev["lvc75_pregate.py"] = sha(Path(__file__).resolve())
    return ev


def verify():
    """--verify: 記録が PASS で、記録したハッシュが今のファイルと同じなら 0、そうでなければ 1"""
    rec = json.loads((HERE / "_band_ab" / "cold_pair" / "lvc75_pregate.json").read_text())
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
    out = HERE / "_band_ab" / "cold_pair" / "lvc75_pregate.json"
    rec = {"plan": "time_integration-line-viscous-jacobian-faceh §6.11 (事前のゲート)"}
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
            bad = {k: c.get(k) for k, v in EXPECT.items() if (k in c if v is None else c.get(k) != v)}
            ok("I", f"{run}: 設定 (本試験との差は step 数・出力の間隔・extraFields だけ)", bad == {}, bad)
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
        ok("I", "3 本の書き出しの設定が完全に同じ (マスクは環境変数)", all(cfgs[k] == cfgs["m7"] for k in DUMPS), None)
        rec["dump_config"] = cfgs["m7"]
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
            ok("I", f"{n}: 7・7 の再実行・5 の 3 本でビット一致 (同じ入力)",
               all(arrays[k][n].shape == arrays["m7"][n].shape and np.array_equal(arrays[k][n].view(np.int64), arrays["m7"][n].view(np.int64)) for k in DUMPS), None)
        # 場全体の入力 (保存量・乱流・組成・物性・壁距離): 出力 step 0 の全データセットが 3 本でビット一致
        def h5vals(run):
            with h5py.File(HERE / run / "res_0.h5", "r") as h:
                return {k: h["VALUE"][k][...] for k in h["VALUE"].keys()}
        v0 = {k: h5vals(DUMPS[k][0]) for k in DUMPS}
        keys0 = sorted(v0["m7"].keys())
        ok("I", "res_0.h5 の全データセット (保存量・k・ω・組成・μ・壁距離など) が 3 本でビット一致",
           len(keys0) > 0 and all(sorted(v0[k].keys()) == keys0 and all(np.array_equal(v0[k][n], v0["m7"][n]) for n in keys0) for k in DUMPS), keys0)
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
            o, o2, n = (res_field(DUMPS[k][0], f) for k in ("m7", "m7b", "m5"))
            if not (np.all(np.isfinite(o)) and np.all(np.isfinite(o2)) and np.all(np.isfinite(n))):
                raise ValueError(f"{f} に非有限がある (出発の状態の残差)")
            ok("I", f"{f}: 3 本で有限", True, None)
            rep_bit = bool(np.array_equal(o.view(np.int64), o2.view(np.int64)))
            d_no, d_oo = float(np.max(np.abs(n - o))), float(np.max(np.abs(o2 - o)))
            cond = bool(np.array_equal(n.view(np.int64), o.view(np.int64))) if rep_bit else d_no <= 3 * d_oo
            ok("R", f"{f}: マスク 5 と 7 の差 ({'ビット一致を要求' if rep_bit else '≤ 3 × 再実行の差'})", cond, [d_no, d_oo])
            resv[f] = {"rerun_bit": rep_bit, "new_vs_old": d_no, "rerun": d_oo}
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
        resv["rhs_s0_record_only"] = {"m5_vs_m7": ulp_stats(arrays["m5"]["rhs_s0"], arrays["m7"]["rhs_s0"]),
                                      "rerun_vs_m7": ulp_stats(arrays["m7b"]["rhs_s0"], arrays["m7"]["rhs_s0"])}
        rec["residual"] = resv

        # ---- V: 介入の成立 (マスクのビット 2 は熱伝導の近傍 K だけを切り替える: block_dplur_jacobian_d.cuh) ----
        ad = HERE / AUDIT / "linedump"
        nn = aud["node_line"].shape[0]
        F = np.fromfile(ad / "audit_face.f64").reshape(nn, NF, FREC)
        N = np.fromfile(ad / "audit_node.f64").reshape(nn, NREC)
        st = arrays["m7"]["state_ro_roU_roe_cp_gamma"]
        nodes = [int(x) for x in arrays["m7"]["node_line"][:, 0]]
        pos = {n: k for k, n in enumerate(nodes)}
        flags = arrays["m7"]["flags_wall_iso_axis"]

        def same_or_noise(name, sl):
            """7 と 7 の再実行がビット一致ならビット一致を、そうでなければ ≤ 3 × 再実行の差を要求する (sl: 列の範囲)"""
            a, a2, b = (arrays[k][name][:, sl] for k in ("m7", "m7b", "m5"))
            rbit = bool(np.array_equal(a.view(np.int64), a2.view(np.int64)))
            d_ba, d_rr = float(np.max(np.abs(b - a))), float(np.max(np.abs(a2 - a)))
            return (bool(np.array_equal(b.view(np.int64), a.view(np.int64))) if rbit else d_ba <= 3 * d_rr), {"rerun_bit": rbit, "m5_vs_m7": d_ba, "rerun": d_rr}
        same = {}
        okD, same["D"] = same_or_noise("D", slice(0, 25))
        ok("V", "(a) D が 7 と 5 で不変 (熱伝導の D は常に入る)", okD, same["D"])
        for nm in ("Kprev", "Knext"):
            okK, same[nm + "_rows0to3"] = same_or_noise(nm, slice(0, 20))
            ok("V", f"(b) {nm} の行 0〜3 が 7 と 5 で不変", okK, same[nm + "_rows0to3"])
        rec["same_D_K03"] = same
        # (c) K の行 4 の変化 (7 − 5) = 熱伝導の近傍 K を、生の double の入力から独立に計算した h と比べる
        #     h = κ (γ_j / c_p,j) / ρ_j · [−(e_j − ½|u_j|²), −u_j, −v_j, −w_j, 1]、κ = max(k_f, 0)·|S|²/|e·S|、
        #     k_f = f0·k0 + (1 − f0)·k1 + (f0·c_p0 + (1 − f0)·c_p1)(f0·μt0 + (1 − f0)·μt1)/Pr_t (面の両端 ic0・ic1 の生の値、f0 = fx)
        #     自節点の行 4 が拘束 (等温壁) か、隣が等温壁 (温度固定) なら h = 0。
        #     比べる対象は記録の側でなく、書き出しのラインの並びから作った期待の接続 (節点 k と向き 0 = prev / 1 = next) の全部 (plan-6 レビュー M2)。
        nl = arrays["m7"]["node_line"]
        exp_conn = []                                  # (k, 向き, 隣の節点番号 or −1)
        for k in range(nn):
            same_prev = k > 0 and nl[k - 1, 1] == nl[k, 1]
            same_next = k + 1 < nn and nl[k + 1, 1] == nl[k, 1]
            exp_conn.append((k, 0, nodes[k - 1] if same_prev else -1))
            exp_conn.append((k, 1, nodes[k + 1] if same_next else -1))
        conn_bad, worst, nres, ncmp, detail = [], 0.0, 0, 0, []
        for k, side, nb in exp_conn:
            nd = N[k]
            if int(nd[0]) != nodes[k] or int(nd[1 + side]) != nb:
                conn_bad.append((nodes[k], side, "記録の line_prev/next が並びと違う", int(nd[1 + side]), nb))
                continue
            # 拘束の行は書き出しのフラグから独立に (壁 → 行 1〜3、等温壁 → 行 4、軸 → 行 2)
            wall, iso, ax = int(flags[k, 0]), int(flags[k, 1]), int(flags[k, 2])
            exp_dec = (14 if wall == 1 else 0) | (16 if iso == 1 else 0) | (4 if ax == 1 else 0)
            if int(nd[18]) != exp_dec:
                conn_bad.append((nodes[k], side, "拘束の行のビットがフラグと違う", int(nd[18]), exp_dec))
                continue
            nm = "Kprev" if side == 0 else "Knext"
            d7, d5 = arrays["m7"][nm][k, 20:25], arrays["m5"][nm][k, 20:25]
            dK = d7 - d5
            if nb < 0:                                 # ラインの端: K は 7・5 とも 0 のはず
                if not (np.all(d7 == 0) and np.all(d5 == 0)):
                    conn_bad.append((nodes[k], side, "ラインの端なのに K の行 4 が 0 でない", d7.tolist(), d5.tolist()))
                continue
            fs = [s_ for s_ in range(int(min(nd[4], NF))) if F[k, s_, 5] != 0 and int(F[k, s_, 3]) == nb]
            if len(fs) != 1:
                conn_bad.append((nodes[k], side, "期待の接続にライン面がちょうど 1 つでない", len(fs), nb))
                continue
            fr = F[k, fs[0]]
            kb = pos.get(nb)
            back = [s_ for s_ in range(int(min(N[kb, 4], NF))) if F[kb, s_, 5] != 0 and int(F[kb, s_, 3]) == nodes[k]] if kb is not None else []
            if not (int(fr[6]) == side and int(fr[45]) == 1 and {int(fr[1]), int(fr[2])} == {nodes[k], nb} and kb is not None
                    and len(back) == 1 and int(F[kb, back[0], 0]) == int(fr[0]) and int(F[kb, back[0], 6]) == 1 - side):
                conn_bad.append((nodes[k], side, "面の向き・枝・両端・往復の面が合わない", None, nb))
                continue
            if bool(fr[44] != 0) != (int(flags[kb, 1]) == 1):
                conn_bad.append((nodes[k], side, "隣の温度固定のフラグが書き出しのフラグと違う", fr[44], flags[kb, 1]))
                continue
            f0 = fr[7]
            kf = f0 * fr[186] + (1 - f0) * fr[187] + (f0 * fr[184] + (1 - f0) * fr[185]) * (f0 * fr[182] + (1 - f0) * fr[183]) / fr[50]
            e = np.array([fr[15] - fr[12], fr[16] - fr[13], fr[17] - fr[14]])
            kap = max(kf, 0.0) * fr[11] ** 2 / abs(float(e @ fr[8:11]))
            rj, uj, vj, wj, rEj, gj = fr[52], fr[53], fr[54], fr[55], fr[56], fr[51]
            cpj = max(st[kb, 5], 1e-30)
            q2 = uj * uj + vj * vj + wj * wj
            ej = rEj / rj - 0.5 * q2
            c = kap * (gj / cpj) / rj
            h = c * np.array([-(ej - 0.5 * q2), -uj, -vj, -wj, 1.0])
            if iso == 1 or int(flags[kb, 1]) == 1:
                h = np.zeros(5)
            scale = np.maximum(np.maximum(np.abs(d7), np.abs(d5)), np.abs(h))
            floor = 64 * np.maximum(ulp32(scale), np.finfo(np.float32).tiny)
            tol = 1e-5 * np.abs(h) + floor   # 要素ごと
            ratio = np.abs(dK - h) / tol
            if not (np.isfinite(kf) and np.isfinite(kap) and np.all(np.isfinite(h)) and np.all(np.isfinite(tol)) and np.all(np.isfinite(ratio))):
                raise ValueError(f"節点 {nodes[k]} 向き {side}: κ・h・許容・比に非有限 (入力の不備)")
            ncmp += 1
            r = float(np.max(ratio))
            worst = max(worst, r)
            if float(np.max(np.abs(h))) > 1e3 * float(np.max(floor)):
                nres += 1
            if r > 1 and len(detail) < 5:
                detail.append({"node": nodes[k], "side": nm, "dK": dK.tolist(), "h": h.tolist()})
        rec["heatK"] = {"expected_connections": len(exp_conn), "compared_line_faces": ncmp, "resolved_faces": nres, "max_err_over_tol": worst,
                        "connection_problems": conn_bad[:10], "n_connection_problems": len(conn_bad), "first_bad": detail}
        if conn_bad:
            raise ValueError(f"期待の接続と監査の記録が合わない ({len(conn_bad)} 件、最初 {conn_bad[0]})")
        ok("V", "(c) 期待の接続の全部で、K の行 4 の変化 (7 − 5) が独立に計算した熱伝導の項と許容内 (要素ごとに 1e-5 × |h| + 64 ulp)、端は 0",
           ncmp > 0 and worst <= 1.0, [ncmp, worst])
        ok("V", "(c') 変化が解像される面 (熱伝導の項が許容の床の 1000 倍超) が 100 面以上", nres >= 100, nres)
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
    for k in ("residual", "struct_rhs", "same_D_K03", "heatK"):
        if k in rec:
            print(k, rec[k])
    print("VERDICT:", verdict)
    return 0 if verdict.startswith("PASS") else (1 if verdict.startswith("INVALID") else 2)


if __name__ == "__main__":
    sys.exit(verify() if sys.argv[1:] == ["--verify"] else main())
