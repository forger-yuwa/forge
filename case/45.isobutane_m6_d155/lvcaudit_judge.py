#!/usr/bin/env python3
"""製品の経路の照合 (plan time_integration-line-viscous-jacobian-faceh §6.7、2026-10-10 事前登録、codex plan-4 の採否を反映) の判定。

監査用のビルド (-DFORGE_LINE_AUDIT) が FORGE_LINE_DUMP_DIR に書いた記録と、同じソース・同じコンパイル条件で監査なしの通常のビルドの書き出し (対照) を読み、
  S  採取のゲート: 型の幅 (FP64)・実効の設定・バイナリと出発の場の sha256・5 本のラインの全節点・sweep 0〜4・必須の配列の形と有限性・正値
  G  構造: 面の接続 (lp/ln の向きごとにちょうど 1 面、往復の対応)・ライン上の並び・拘束のフラグ・隣の生の状態と物性・速度/温度固定の隣のフラグを、
     書き出した元の配列から独立に照合 (完全一致)
  B  係数の float の再現: 生の入力から製品と同じ ST = float の式を numpy の float32 で再現 (全面の面積・法線・dcc・δ・ν_eff・スカラーの粘性の対角、
     薄層の f_i・μ_f・k_f・β・κ、熱伝導の Jacobian の k_face・cfac、軸対称の A_pl・r_eff・μ_total・hoop、時間項の v・dt)。相対 ≤ 1e-6 合格、≤ 1e-4 判別不能、それより大きい = 不一致
  C  面ごとの寄与: line_audit_helper が記録の「この面の前の D」から製品と同じ順序の float の加算で再現した「この面の後の D」、対流の K、薄層の D・K と CUDA の値。
     要素ごとに ulp (丸めの尺度の float の ulp) で ≤ 8 合格、≤ 64 判別不能、それより大きい = 不一致
  D  組立: 時間項 (生の体積・dt_local から)、面の連なり (前の D = 一つ前の面の後の D、最初 = 時間項の後、最後 = 面のループの後: 完全一致)、
     時間項から連ねた再現と面のループの後の D (ulp)、軸対称 (ulp)、拘束の後の D・格納の D・Kprev/Knext の値と置き場所 (完全一致)
  A  元の入力から double で独立に計算した係数と、CUDA が使った ST の値の相対差 (特性の記録): 薄層の面の β・κ、ほかの面のスカラーの粘性の対角・熱伝導の k_face·δ/dcc。
     別に dcc・δ/dcc・法線の差も記録する。差 ≤ 1e-5 / ≤ 1e-3 / > 1e-3 に分ける
  T  対照: 監査用のビルドと通常のビルドの D・Kprev・Knext・状態・rhs・dq (sweep 0〜4) などのビット一致 (整数の表現で比べる)
を判定し、<dump>/lvcaudit_judge.json と標準出力に書く。規則を変えるときは plan §6.7 を先に改訂する。
使い方: lvcaudit_judge.py <dump> <line_audit_helper> --run <run dir> --expect-sha <監査用のビルドの forge の sha256> --ctrl <対照の dump> --ctrl-run <対照の run dir> --ctrl-sha <通常のビルドの sha256>
"""
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

NF, FREC, NREC, HF, HN = 12, 200, 160, 150, 75
SRC_SHA16 = "207d39f0e7f4aa03"
WANT_NODES = [1572, 4113, 7985, 198560, 264263]
EXPECT_CFG = {"implicitSolvePrecision": "0", "lineViscCoupling": "3", "implicitThermalJacobian": "5", "lineViscTerms": "7", "factorCall": "1",
              "targetCall": "1", "sizeof_flow_float": "8", "sizeof_geom_float": "8", "lineImplicit": "1", "lineDtDirectional": "1",
              "lineDtDirectionalCap": "0", "isAxisymmetric": "1", "axisymMethod": "0", "thermalMethod": "2"}
B_PASS, B_IND = 1e-6, 1e-4
ULP_PASS, ULP_IND = 8, 64
A1, A2 = 1e-5, 1e-3
CTRL_ARRAYS = ["D", "Kprev", "Knext", "state_ro_roU_roe_cp_gamma", "dt_vol", "scalar_visc_line", "flags_wall_iso_axis"] + \
              [f"{p}_s{k}" for p in ("rhs", "dqnew", "dqold") for k in range(5)]
f32 = np.float32


class Checks:
    def __init__(self):
        self.items = {}

    def add(self, stage, name, cls, val=None):
        """cls: 'pass' / 'ind' / 'fail'"""
        self.items.setdefault(stage, []).append((name, cls, val))

    def ok(self, stage, name, cond, val=None):
        self.add(stage, name, "pass" if cond else "fail", val)

    def summary(self, stage):
        lst = self.items.get(stage, [])
        out = {"n": len(lst), "fail": sum(c == "fail" for _, c, _ in lst), "ind": sum(c == "ind" for _, c, _ in lst)}
        out["first_fail"] = [[n, v] for n, c, v in lst if c == "fail"][:10]
        out["first_ind"] = [[n, v] for n, c, v in lst if c == "ind"][:10]
        return out


def load_meta(d):
    meta, cfg = {}, {}
    for line in (d / "meta.txt").read_text().splitlines():
        if line.startswith("# audit_cfg"):
            toks = line.split()[2:]
            cfg = {toks[i]: toks[i + 1] for i in range(0, len(toks) - 1, 2)}
            continue
        if line.startswith("#") or not line.strip():
            continue
        name, r, c = line.split()
        meta[name] = (int(r), int(c))
    return meta, cfg


def arr(d, meta, name):
    r, c = meta[name]
    a = np.fromfile(d / f"{name}.f64", dtype=np.float64)
    if a.size != r * c:
        raise ValueError(f"{name}: 大きさが違う ({a.size} != {r}×{c})")
    return a.reshape(r, c)


def ulp32(x):
    return np.spacing(np.abs(np.asarray(x, dtype=np.float32))).astype(np.float64)


def ulp_class(dev, host, scale):
    """要素ごとに |dev − host| / ulp_float(scale) の最大で分類する。"""
    dev, host, scale = (np.asarray(x, dtype=np.float64) for x in (dev, host, scale))
    u = np.maximum(ulp32(np.maximum(np.abs(scale), np.abs(dev))), np.finfo(np.float32).tiny)
    m = float(np.max(np.abs(dev - host) / u)) if dev.size else 0.0
    return ("pass" if m <= ULP_PASS else "ind" if m <= ULP_IND else "fail"), m


def rel_class(dev, ref):
    dev, ref = float(dev), float(ref)
    e = abs(dev - ref) / abs(ref) if ref != 0 else abs(dev)
    return ("pass" if e <= B_PASS else "ind" if e <= B_IND else "fail"), e


def sha_of(run):
    p = Path(run) / "RUN_PROVENANCE.txt"
    m = re.search(r"forge_sha256\s*:\s*([0-9a-f]{64})", p.read_text()) if p.is_file() else None
    return m.group(1) if m else None


def main():
    a = sys.argv[1:]
    if len(a) < 2:
        raise SystemExit(__doc__)
    opt = lambda k: a[a.index(k) + 1] if k in a else None
    d, helper = Path(a[0]), a[1]
    run, ctrl, ctrl_run = opt("--run"), opt("--ctrl"), opt("--ctrl-run")
    exp_sha, ctrl_sha = opt("--expect-sha"), opt("--ctrl-sha")
    out = d / "lvcaudit_judge.json"
    rec = {"plan": "time_integration-line-viscous-jacobian-faceh §6.7", "dump": str(d)}
    ch = Checks()
    try:
        meta, cfg = load_meta(d)
        rec["audit_cfg"] = cfg
        # ---- S: 採取のゲート ----
        for k, v in EXPECT_CFG.items():
            ch.ok("S", f"設定 {k} = {v}", cfg.get(k) == v, cfg.get(k))
        ch.ok("S", "引数 --run・--expect-sha・--ctrl・--ctrl-run・--ctrl-sha がある", all([run, exp_sha, ctrl, ctrl_run, ctrl_sha]), None)
        if run:
            ch.ok("S", "監査用のビルドの sha256", sha_of(run) == exp_sha, sha_of(run))
            cp = json.loads((Path(run) / "COLD_PAIR.json").read_text())
            ch.ok("S", "出発の場が run_0183 の res_100000 (sha256 固定)", cp.get("parent_res") == "res_100000.h5" and cp.get("parent_res_sha256", "").startswith(SRC_SHA16),
                  [cp.get("parent_res"), cp.get("parent_res_sha256", "")[:16]])
        if ctrl_run:
            ch.ok("S", "通常のビルドの sha256", sha_of(ctrl_run) == ctrl_sha, sha_of(ctrl_run))
            cpc = json.loads((Path(ctrl_run) / "COLD_PAIR.json").read_text())
            ch.ok("S", "対照の出発の場が同じ", cpc.get("parent_res_sha256") == cp.get("parent_res_sha256") if run else False, cpc.get("parent_res_sha256", "")[:16])
        need = ["node_line", "D", "Kprev", "Knext", "flags_wall_iso_axis", "state_ro_roU_roe_cp_gamma", "dt_vol", "node_vel_props", "node_cc",
                "audit_face", "audit_node"] + [f"{p}_s{k}" for p in ("rhs", "dqnew", "dqold") for k in range(5)]
        miss = [n for n in need if n not in meta]
        ch.ok("S", "必須の配列がそろう (sweep 0〜4 を含む)", not miss, miss)
        if miss:
            raise ValueError(f"必須の配列が無い {miss}")
        nl = arr(d, meta, "node_line"); nn = nl.shape[0]
        ch.ok("S", "節点が 0 でない", nn > 0, nn)
        nodes = [int(x) for x in nl[:, 0]]
        pos = {n: k for k, n in enumerate(nodes)}
        ch.ok("S", "節点の重複なし", len(pos) == nn, nn - len(pos))
        lines = sorted(set(int(x) for x in nl[:, 1]))
        ch.ok("S", "要求した 5 節点がすべて含まれる", all(w in pos for w in WANT_NODES), [w for w in WANT_NODES if w not in pos])
        ch.ok("S", "ラインが 5 本", len(lines) == 5, len(lines))
        D, Kp, Kn = arr(d, meta, "D"), arr(d, meta, "Kprev"), arr(d, meta, "Knext")
        flags, st, dtv = arr(d, meta, "flags_wall_iso_axis"), arr(d, meta, "state_ro_roU_roe_cp_gamma"), arr(d, meta, "dt_vol")
        props, cc = arr(d, meta, "node_vel_props"), arr(d, meta, "node_cc")
        F = arr(d, meta, "audit_face").reshape(nn, NF, FREC)
        N = arr(d, meta, "audit_node")
        for name, x in (("D", D), ("Kprev", Kp), ("Knext", Kn), ("状態", st), ("dt_vol", dtv), ("物性", props), ("座標", cc), ("節点の記録", N)):
            ch.ok("S", f"{name} が有限", bool(np.all(np.isfinite(x))), int(np.sum(~np.isfinite(x))))
        ch.ok("S", "ρ・dt・体積・c_p > 0、γ > 1", bool(np.all(st[:, 0] > 0) and np.all(dtv > 0) and np.all(st[:, 5] > 0) and np.all(st[:, 6] > 1)), None)
        for k in range(nn):
            nd = N[k]
            ch.ok("S", f"節点 {k}: 記録の節点番号・storeLU・loop 0・ST = float・型の幅 8", int(nd[0]) == nodes[k] and nd[136] == 1 and nd[137] == 0 and nd[138] == 4
                  and nd[147] == 8 and nd[148] == 8, [nd[0], nd[136], nd[137], nd[138], nd[147], nd[148]])
            ch.ok("S", f"節点 {k}: 面の数が枠に入る (あふれ 0、数 = plane の数)", nd[5] == 0 and nd[4] == nd[146] and nd[4] <= NF, [nd[4], nd[5], nd[146]])
            for s in range(int(min(nd[4], NF))):
                fr = F[k, s]
                ch.ok("S", f"節点 {k} 面 {s}: 枠の番号・既知の粘性の枝・有限", int(fr[191]) == s and int(fr[45]) in (0, 1, 3, 4) and bool(np.all(np.isfinite(fr))),
                      [fr[191], fr[45]])
        if any(c == "fail" for _, c, _ in ch.items["S"]):
            raise ValueError("採取のゲートが不合格")

        # ---- G: 構造 (元の配列から独立に) ----
        def face_of(k, other):
            return [s for s in range(int(min(N[k, 4], NF))) if F[k, s, 5] != 0 and int(F[k, s, 3]) == other]
        for k in range(nn):
            nd = N[k]
            ic, lp, ln = nodes[k], int(nd[1]), int(nd[2])
            for side, nb in ((0, lp), (1, ln)):
                if nb < 0:
                    continue
                fs = face_of(k, nb)
                ch.ok("G", f"節点 {k}: 向き {side} の隣 {nb} にちょうど 1 面", len(fs) == 1, len(fs))
                if len(fs) != 1:
                    continue
                fr = F[k, fs[0]]
                ch.ok("G", f"節点 {k}: 向き {side} のフラグ", int(fr[6]) == side, fr[6])
                ch.ok("G", f"節点 {k}: 面の両端が自分と隣", {int(fr[1]), int(fr[2])} == {ic, nb}, [fr[1], fr[2]])
                ch.ok("G", f"節点 {k}: 隣 {nb} が書き出しの中", nb in pos, None)
                if nb in pos:
                    kb = pos[nb]
                    back = face_of(kb, ic)
                    ch.ok("G", f"節点 {k}: 隣から見た往復の面 (同じ面の番号・逆向き)", len(back) == 1 and int(F[kb, back[0], 0]) == int(fr[0])
                          and int(N[kb, 2 if side == 0 else 1]) == ic, [len(back)])
                    # 隣の生の状態・物性は元の配列と完全一致
                    ch.ok("G", f"節点 {k}: 隣の生の状態 (γ・ρ・u・ρE)", fr[51] == st[kb, 6] and fr[52] == st[kb, 0] and fr[53] == props[kb, 0] and fr[54] == props[kb, 1]
                          and fr[55] == props[kb, 2] and fr[56] == st[kb, 4], None)
                    if int(fr[45]) == 1:
                        ch.ok("G", f"節点 {k}: 隣の速度固定のフラグ = 隣の壁のフラグ", (fr[43] != 0) == (flags[kb, 0] == 1), [fr[43], flags[kb, 0]])
                        ch.ok("G", f"節点 {k}: 隣の温度固定のフラグ = 隣の等温壁のフラグ", (fr[44] != 0) == (flags[kb, 1] == 1), [fr[44], flags[kb, 1]])
                        ch.ok("G", f"節点 {k}: マスク・値 3", fr[48] == 7 and fr[49] == 1, [fr[48], fr[49]])
                        ch.ok("G", f"節点 {k}: 両端の生の層流 μ が元の配列", fr[188] == props[k, 3] and fr[189] == props[kb, 3], None)
            nlines = sum(1 for s in range(int(min(nd[4], NF))) if F[k, s, 5] != 0)
            ch.ok("G", f"節点 {k}: ライン面の数 = 隣の数", nlines == (lp >= 0) + (ln >= 0), [nlines, lp, ln])
            if k + 1 < nn and nl[k + 1, 1] == nl[k, 1]:
                ch.ok("G", f"節点 {k}: 書き出しの次の節点がライン上の次 (ln)", nodes[k + 1] == ln, [nodes[k + 1], lp, ln])
            wall, iso, ax = int(flags[k, 0]), int(flags[k, 1]), int(flags[k, 2])
            exp = (14 if wall == 1 else 0) | (16 if iso == 1 else 0) | (4 if ax == 1 else 0)
            ch.ok("G", f"節点 {k}: 拘束の行のビット = 元のフラグ", int(nd[18]) == exp and nd[19] == wall and nd[20] == iso and nd[21] == ax, [nd[18], exp])
            ch.ok("G", f"節点 {k}: 生の ρ・ρE・体積・dt が元の配列", nd[17] == st[k, 0] and nd[16] == st[k, 4] and nd[133] == dtv[k, 1] and nd[134] == dtv[k, 0], None)
            ch.ok("G", f"節点 {k}: 生の μ_t が元の配列", nd[132] == props[k, 4], None)
        if any(c == "fail" for _, c, _ in ch.items.get("G", [])):
            raise ValueError("構造の照合が不合格")

        # ---- helper ----
        r = subprocess.run([helper, str(d), str(nn), "1"], capture_output=True, text=True)
        if r.returncode != 0:
            raise ValueError(f"line_audit_helper が失敗: {r.stderr}")
        HFc = np.fromfile(d / "audit_host_face.f64", dtype=np.float64).reshape(nn, NF, HF)
        HNn = np.fromfile(d / "audit_host_node.f64", dtype=np.float64).reshape(nn, HN)

        Arows = []
        visc = float(cfg["visc"])
        for k in range(nn):
            nd = N[k]
            ic = nodes[k]
            rho = max(f32(st[k, 0]), f32(1e-30)); g_i = f32(st[k, 6])
            ux, uy, uz = f32(props[k, 0]), f32(props[k, 1]), f32(props[k, 2])
            # B: 節点の量
            for name, dev, ref in (("ρ", nd[9], rho), ("γ", nd[15], g_i), ("u", nd[10], ux), ("v", nd[11], uy), ("w", nd[12], uz),
                                   ("音速", nd[13], max(f32(props[k, 6]), f32(1e-8))), ("H_t", nd[14], f32(props[k, 7])), ("体積", nd[6], f32(dtv[k, 1])), ("dt", nd[7], f32(dtv[k, 0]))):
                ch.add("B", f"節点 {k}: {name}", *rel_class(dev, ref))
            nu_eff = (f32(visc) + max(f32(props[k, 4]), f32(0))) / rho
            ch.add("B", f"節点 {k}: ν_eff", *rel_class(nd[23], nu_eff))
            if nd[124] == 1:
                A = f32(nd[128]); v = f32(dtv[k, 1])
                r_eff = max(v / max(A, f32(1e-30)), f32(1e-30))
                mu_t = f32(visc) + max(f32(props[k, 4]), f32(0))
                hoop = f32(2) * mu_t / (rho * r_eff)
                for name, dev, ref in (("A_pl", nd[125], A), ("r_eff", nd[129], r_eff), ("μ_total", nd[130], mu_t), ("hoop", nd[126], hoop)):
                    ch.add("B", f"節点 {k}: 軸対称 {name}", *rel_class(dev, ref))
            ch.ok("B", f"節点 {k}: 軸対称の枝が r 重み (1) か、軸対称の対象外 (0)", nd[124] in (0, 1), nd[124])
            for s in range(int(min(nd[4], NF))):
                fr = F[k, s]; br = int(fr[45]); has_nbr = fr[4] != 0
                ic0, ic1, o = int(fr[1]), int(fr[2]), int(fr[3])
                sgn = f32(1) if ic0 == ic else f32(-1)
                fa = max(f32(fr[11]), f32(1e-30))
                nf = (sgn * f32(fr[8]) / fa, sgn * f32(fr[9]) / fa, sgn * f32(fr[10]) / fa)
                ch.add("B", f"節点 {k} 面 {s}: 面積", *rel_class(fr[18], fa))
                for q in range(3):
                    cls, e = rel_class(fr[19 + q], nf[q]) if nf[q] != 0 else (("pass" if fr[19 + q] == 0 else "fail"), abs(fr[19 + q]))
                    ch.add("B", f"節点 {k} 面 {s}: 法線 {q}", cls, e)
                if not has_nbr:
                    ch.ok("B", f"節点 {k} 面 {s}: 境界の半割面は粘性の枝なし", br == 0, br)
                    continue
                dx, dy, dz = f32(fr[15]) - f32(fr[12]), f32(fr[16]) - f32(fr[13]), f32(fr[17]) - f32(fr[14])
                dcc = max(np.sqrt(dx * dx + dy * dy + dz * dz), f32(1e-30))
                dds = max(abs(dx * f32(fr[8]) + dy * f32(fr[9]) + dz * f32(fr[10])), f32(1e-30))
                delta = max(dcc * fa * fa / dds, f32(1e-30))
                vd = f32(2) * nu_eff * delta / dcc
                for name, dev, ref in (("dcc", fr[25], dcc), ("dcc·s", fr[26], dds), ("δ", fr[27], delta), ("スカラーの粘性の対角", fr[28], vd), ("ν_eff", fr[29], nu_eff)):
                    ch.add("B", f"節点 {k} 面 {s}: {name}", *rel_class(dev, ref))
                # A: 元の入力から double で
                dxd = np.array([fr[15] - fr[12], fr[16] - fr[13], fr[17] - fr[14]])
                dccd = float(np.linalg.norm(dxd)); ddsd = abs(float(dxd @ fr[8:11])); dod = fr[11] ** 2 / ddsd
                nud = (visc + max(props[k, 4], 0.0)) / st[k, 0]
                arow = {"node": k, "ic": ic, "other": o, "slot": s, "line": int(nl[k, 1]), "branch": br, "isLine": bool(fr[5] != 0),
                        "rel_dcc": abs(fr[25] / dccd - 1), "rel_delta_over_dcc": abs((fr[27] / fr[25]) / dod - 1),
                        "rel_normal": float(np.max(np.abs(fr[19:22] - (1 if ic0 == ic else -1) * fr[8:11] / fr[11]))),
                        "y_i": fr[13], "y_o": fr[16], "dcc_double": dccd, "dcc_float": fr[25], "ulp_y": float(np.spacing(np.float32(max(abs(fr[13]), abs(fr[16])))))}
                if br == 1:
                    kb = pos.get(o)
                    fi = f32(fr[7]) if ic0 == ic else f32(1) - f32(fr[7]); omf = f32(1) - fi
                    mli, mlj = f32(fr[188]), f32(fr[189]); mti, mtj = f32(props[k, 4]), f32(props[kb, 4])
                    cpi, cpj = f32(st[k, 5]), f32(st[kb, 5]); tci, tcj = f32(props[k, 5]), f32(props[kb, 5])
                    muf = fi * (mli + mti) + omf * (mlj + mtj)
                    kf = fi * tci + omf * tcj + (fi * cpi + omf * cpj) * (fi * mti + omf * mtj) / f32(fr[50])
                    beta, kappa = max(muf, f32(0)) * delta / dcc, max(kf, f32(0)) * delta / dcc
                    for name, dev, ref in (("f_i", fr[30], fi), ("μ_f", fr[39], muf), ("k_f", fr[40], kf), ("β", fr[41], beta), ("κ", fr[42], kappa),
                                           ("層流 μ_i", fr[31], mli), ("層流 μ_j", fr[32], mlj), ("μ_t,i", fr[33], mti), ("μ_t,j", fr[34], mtj), ("c_p,i", fr[35], cpi), ("c_p,j", fr[36], cpj)):
                        ch.add("B", f"節点 {k} 面 {s}: {name}", *rel_class(dev, ref))
                    fid = fr[7] if ic0 == ic else 1 - fr[7]
                    mufd = fid * (props[k, 3] + props[k, 4]) + (1 - fid) * (props[kb, 3] + props[kb, 4])
                    kfd = fid * props[k, 5] + (1 - fid) * props[kb, 5] + (fid * st[k, 5] + (1 - fid) * st[kb, 5]) * (fid * props[k, 4] + (1 - fid) * props[kb, 4]) / fr[50]
                    arow["rel_beta"] = abs(fr[41] / (max(mufd, 0) * dod) - 1); arow["rel_kappa"] = abs(fr[42] / (max(kfd, 0) * dod) - 1)
                    arow["coef"] = max(arow["rel_beta"], arow["rel_kappa"])
                elif br == 3:
                    f = f32(fr[7]); omf = f32(1) - f
                    cpf = f * f32(fr[184]) + omf * f32(fr[185]); mtf = f * f32(fr[182]) + omf * f32(fr[183])
                    kface = f * f32(fr[186]) + omf * f32(fr[187]) + cpf * mtf / f32(fr[50])
                    cfac = (kface * delta / dcc) * (g_i / max(f32(st[k, 5]), f32(1e-30))) / rho
                    ch.add("B", f"節点 {k} 面 {s}: k_face", *rel_class(fr[46], kface))
                    ch.add("B", f"節点 {k} 面 {s}: cfac", *rel_class(fr[47], cfac))
                    ch.ok("B", f"節点 {k} 面 {s}: スカラーの行の数 = 5 (キー 5)", fr[190] == 5, fr[190])
                    kfd = fr[7] * fr[186] + (1 - fr[7]) * fr[187] + (fr[7] * fr[184] + (1 - fr[7]) * fr[185]) * (fr[7] * fr[182] + (1 - fr[7]) * fr[183]) / fr[50]
                    arow["rel_vd"] = abs(fr[28] / (2 * nud * dod) - 1); arow["rel_kface_dod"] = abs((fr[46] * fr[27] / fr[25]) / (kfd * dod) - 1)
                    arow["coef"] = max(arow["rel_vd"], arow["rel_kface_dod"])
                else:
                    arow["rel_vd"] = abs(fr[28] / (2 * nud * dod) - 1); arow["coef"] = arow["rel_vd"]
                Arows.append(arow)
                # C: 面ごとの寄与
                hh = HFc[k, s]
                cls, m = ulp_class(fr[157:182], hh[0:25], hh[50:75])
                ch.add("C", f"節点 {k} 面 {s} (枝 {br}): この面の後の D (前の D から再現)", cls, m)
                if fr[5] != 0:
                    kc_dev, kc_h = fr[82:107].reshape(5, 5), hh[75:100].reshape(5, 5)
                    cls, m = ulp_class(kc_dev, kc_h, np.broadcast_to(np.max(np.abs(kc_h), axis=0), (5, 5)))
                    ch.add("C", f"節点 {k} 面 {s}: 対流の K", cls, m)
                    if br == 1:
                        for nm, dv, hv in (("薄層の D", fr[132:157], hh[100:125]), ("薄層の K", fr[107:132], hh[125:150])):
                            dv5, hv5 = dv.reshape(5, 5), hv.reshape(5, 5)
                            cls, m = ulp_class(dv5, hv5, np.broadcast_to(np.max(np.abs(hv5), axis=0), (5, 5)))
                            ch.add("C", f"節点 {k} 面 {s}: {nm}", cls, m)

        # ---- D: 組立 ----
        for k in range(nn):
            nd = N[k]; nfc = int(min(nd[4], NF)); hn = HNn[k]
            cls, m = ulp_class(nd[24:49], hn[0:25], nd[24:49]); ch.add("D", f"節点 {k}: D (時間項)", cls, m)
            if nfc:
                ch.ok("D", f"節点 {k}: 最初の面の前の D = 時間項の後の D", bool(np.array_equal(F[k, 0, 57:82], nd[24:49])), None)
                for s in range(1, nfc):
                    ch.ok("D", f"節点 {k}: 面 {s} の前の D = 面 {s - 1} の後の D", bool(np.array_equal(F[k, s, 57:82], F[k, s - 1, 157:182])), None)
                ch.ok("D", f"節点 {k}: 最後の面の後の D = 面のループの後の D", bool(np.array_equal(F[k, nfc - 1, 157:182], nd[49:74])), None)
                scale = np.max(np.stack([HFc[k, s, 50:75] for s in range(nfc)]), axis=0)
                cls, m = ulp_class(nd[49:74], HFc[k, nfc - 1, 25:50], scale)
                ch.add("D", f"節点 {k}: 面のループの後の D (時間項から連ねて再現)", cls, m)
            cls, m = ulp_class(nd[74:99], hn[25:50], hn[50:75]); ch.add("D", f"節点 {k}: D (軸対称の後)", cls, m)
            Dc = nd[74:99].reshape(5, 5).copy()
            if nd[21] == 1:
                Dc[2, :] = 0; Dc[2, 2] = 1
            if nd[19] == 1:
                for row in (1, 2, 3):
                    Dc[row, :] = 0; Dc[row, row] = 1
            if nd[20] == 1:
                Dc[4, :] = 0; Dc[4, 4] = 1; Dc[4, 0] = float(-f32(nd[16]) / max(f32(nd[17]), f32(1e-30)))
            ch.ok("D", f"節点 {k}: 拘束の後の D (完全一致)", bool(np.array_equal(Dc.ravel(), nd[99:124])), float(np.max(np.abs(Dc.ravel() - nd[99:124]))))
            ch.ok("D", f"節点 {k}: 格納の D = 拘束の後の D (完全一致)", bool(np.array_equal(D[k], nd[99:124])), None)
            rd = int(nd[18])
            for side, nb, K in ((0, int(nd[1]), Kp), (1, int(nd[2]), Kn)):
                if nb < 0:
                    continue
                s = [x for x in range(nfc) if F[k, x, 5] != 0 and int(F[k, x, 3]) == nb][0]
                fr = F[k, s]
                exp = np.zeros(25)
                for i in range(5):
                    if rd & (1 << i):
                        continue
                    for j in range(5):
                        v = float(np.float64(np.float32(fr[82 + i * 5 + j])))
                        if fr[45] == 1:
                            v = v + float(np.float64(np.float32(fr[107 + i * 5 + j])))
                        exp[i * 5 + j] = v
                ch.ok("D", f"節点 {k}: {'Kprev' if side == 0 else 'Knext'} の値と置き場所 (完全一致)", bool(np.array_equal(K[k].view(np.int64), exp.view(np.int64))),
                      float(np.max(np.abs(K[k] - exp))))

        # ---- T: 対照 (同じソース・監査なしの通常のビルド) とのビット一致 ----
        bit = {}
        if ctrl:
            mc, _ = load_meta(Path(ctrl))
            ncl = arr(Path(ctrl), mc, "node_line")
            same_nodes = bool(np.array_equal(ncl, nl))
            ch.ok("T", "対照の書き出しの節点と順が同じ", same_nodes, None)
            for name in CTRL_ARRAYS:
                if name not in meta or name not in mc:
                    ch.ok("T", f"{name} が両方にある", False, None); continue
                x, y = arr(d, meta, name), arr(Path(ctrl), mc, name)
                eq = x.shape == y.shape and bool(np.array_equal(x.view(np.int64), y.view(np.int64)))
                bit[name] = eq
                ch.ok("T", f"{name} がビット一致", eq, None)
        rec["bit_identical_vs_ctrl"] = bit

        # ---- A の分類 ----
        classes = {"≤1e-5": 0, "≤1e-3": 0, ">1e-3": 0}
        for row in Arows:
            c = "≤1e-5" if row["coef"] <= A1 else ("≤1e-3" if row["coef"] <= A2 else ">1e-3")
            row["class"] = c; classes[c] += 1
        Arows.sort(key=lambda r: -r["coef"])
        rec["A"] = {"n_faces": len(Arows), "classes": classes, "worst": Arows[:40],
                    "line_faces": {"n": sum(r["isLine"] for r in Arows), "gt_1e-3": sum(r["isLine"] and r["coef"] > A2 for r in Arows),
                                   "max_rel_beta": max((r.get("rel_beta", 0) for r in Arows), default=0), "max_rel_kappa": max((r.get("rel_kappa", 0) for r in Arows), default=0)},
                    "max_rel_dcc": max((r["rel_dcc"] for r in Arows), default=0), "max_rel_delta_over_dcc": max((r["rel_delta_over_dcc"] for r in Arows), default=0),
                    "max_abs_normal": max((r["rel_normal"] for r in Arows), default=0)}
        stages = {s: ch.summary(s) for s in ("S", "G", "B", "C", "D", "T")}
        rec.update(stages)
        if any(stages[s]["n"] == 0 for s in ("S", "G", "B", "C", "D")):
            verdict = "INVALID (照合の項目が 0 の段がある)"
        elif any(stages[s]["fail"] for s in ("B", "C", "D")):
            first = next(s for s in ("B", "C", "D") if stages[s]["fail"])
            verdict = f"FAIL (製品の組立が記録の入力からの再現と許容外: 最初に外れる段 {first})"
        elif any(stages[s]["ind"] for s in ("B", "C", "D")):
            verdict = "INDETERMINATE (丸めの差で説明できる範囲の外れがあり、合格とも不一致とも言えない)"
        else:
            verdict = "PASS (製品の組立は、記録の入力からの再現と許容内)"
        if stages["T"]["n"] == 0 or stages["T"]["fail"]:
            verdict += "; 対照とのビット一致が成り立たないので、通常の経路への結論は保留"
        rec["VERDICT"] = verdict
        rec["A_VERDICT"] = (f"元の入力から double で計算した係数との差 (面 {len(Arows)}): ≤1e-5 {classes['≤1e-5']}・≤1e-3 {classes['≤1e-3']}・>1e-3 {classes['>1e-3']}"
                            f" (ライン面の β・κ で >1e-3 は {rec['A']['line_faces']['gt_1e-3']})。再現の PASS と両立しうる。破綻の原因の断定には使わない")
    except Exception as e:   # 採取・構造の不備は理由付きの INVALID
        rec["VERDICT"] = f"INVALID ({type(e).__name__}: {e})"
        for s in ("S", "G"):
            if s in ch.items:
                rec[s] = ch.summary(s)
    out.write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=float))
    for s in ("S", "G", "B", "C", "D", "T"):
        if s in rec:
            x = rec[s]
            print(f"{s}: {x['n']} 項目、不一致 {x['fail']}、判別不能 {x['ind']}  {x['first_fail'][:3]} {x['first_ind'][:2]}")
    if "A" in rec:
        print("A:", rec["A_VERDICT"])
        for row in rec["A"]["worst"][:8]:
            print(f"  ライン {row['line']} 節点 {row['node']} ({row['ic']}→{row['other']}、枝 {row['branch']}{'・ライン面' if row['isLine'] else ''}): 係数 {row['coef']:.2e}、"
                  f"dcc {row['rel_dcc']:.2e}、δ/dcc {row['rel_delta_over_dcc']:.2e} (dcc double {row['dcc_double']:.3e} / float {row['dcc_float']:.3e}、y の ulp {row['ulp_y']:.1e})")
    print("VERDICT:", rec["VERDICT"])
    return 0 if rec["VERDICT"].startswith("PASS") else 1


if __name__ == "__main__":
    sys.exit(main())
