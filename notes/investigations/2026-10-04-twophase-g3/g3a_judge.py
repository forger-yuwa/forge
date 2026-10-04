"""plan condensation-two-phase-default §5.1 #4g3 (G3-a 作用素の収支) の判定と #4pjg (射影節点の残差の段ごとの分解)。
  python3 g3a_judge.py <ON tp_operator.h5> <OFF tp_operator.h5> [--length-scale 1000] [--dF-err w=..,g=..] [--dF w=..,g=..]
                       [--mask projection_nodes.npy] [--mask-arm on|off|both] [--csv out.csv]

入力は FORGE_DIAG_TP_OPERATOR=<h5> で forge が書いた記録 (main.cpp runTpOperatorDiag; 配置は
solver_density_cuda/cuda_forge/twoPhaseOperatorDiag_d.cuh)。ON = run_0561 の設定 + その res_48000、OFF = run_0567 の設定 + その res_48000。
組立 1 回 (前処理 + 後処理、更新なし) の中で本番カーネルが実際に残差へ足した値 (面: atomic の前の値と所属 ic0/ic1、節点: 凝縮ソースの S·V)
と、記録点ごとの残差の写しを読む。すべて double で積算する。

成分: w (ρY_w 総水分), g (ρg 液), Q2, Q1, Q0, 水以外の化学種 (h5 の components)、v = w − g は各項を double で差し引いて作る (記録)。

Ω (実節点 CV の集合; ghost を含めない): all、x 帯 [入口, 17)・[17, 35)・[35, 70)・[70, 出口] mm、x 35–70 mm の壁距離 < 0.4 mm と ≥ 0.4 mm、
  参考に入口ピン節点。長さは ccx・wall_dist に --length-scale (既定 1000 = m → mm) を掛けた値で切る (x の範囲を最初に表示する)。
  境界面は面の所属から決める (片端だけが Ω にある面; 両端が Ω の内部面は相殺)。x 帯の境界は CV 集合を切る双対面 (補間断面ではない)。
  node 境界半割面 (ic1 = ghost) の移流は残り (ic0 だけに足す)、拡散は足さない (code 2)。

恒等式 (Ω ごと・成分ごと):
  Σ_Ω R_final = −F_adv − F_diff + S⁺ + S⁻ + B_pin + R_other + E_assembly
  - R_final = 記録点 final の残差 (格納値; 体積を再乗算しない)。
  - F_adv・F_diff = Ω の境界面の外向き流束 (面の記録 a0 から: ic0 ∈ Ω, ic1 ∉ Ω なら −a0、ic1 ∈ Ω, ic0 ∉ Ω なら +a0)。
  - S⁺・S⁻ = 凝縮ソースが実際に足した項の正・負の部分 (液・Q のみ; 早期退出は 0 が明示的に書かれている)。
  - B_pin = ピン (残差の 0 上書き) の前後差 (写しから; 上書きなので丸めなし)。化学種は 1 回目・2 回目の speciesPinResidual、液・Q は passivePinResidual。
  - R_other = 命名した「この成分に触れないはずの段」の写しの差の和 (化学ソース・トレーサ・後段など)。0 でなければ未記録の書き手 (記録の欠落)
    として (1) を FAIL にする (E_assembly に混ぜない)。
  - E_assembly = 左辺 − 右辺 (float の加算の丸め + double の積算)。
組立誤差の上界 (節点 i への加算項 a_ij、加算回数 m_i; 内部面は両端で 1 回ずつ数える):
  B_assembly = Σ_i γ_{m_i} Σ_j |a_ij| + B_underflow + B_double、γ_m = m u/(1 − m u)、u = 2⁻²⁴。
  凝縮ソースの float 実体は積と和で 2 回丸める (FMA なら 1 回) ので 1 項を 2 回と数える。B_underflow = Σ_i m_i · 2⁻¹²⁶
  (float atomic は非正規化数を 0 にする; 非正規化数の項の数は別に表示)。B_double = judge 自身の double の積算の上界 (γ⁽⁵³⁾)。
合格 (1): |E_assembly| ≤ B_assembly (全 Ω・全成分)。記録の網羅 (経路にある kernel が全面・全節点を書いた) が崩れていれば INVALID。

物理の収支 (2): F_c,exit = 出口マーカー (bcond の種別が outlet で始まる) の node 境界半割面の外向き移流流束 (−a0) の和。出口半割面の拡散は 0
  (code 2 を確認する)。ΔF_c = |F_c,exit(ON) − F_c,exit(OFF)|。合格 (2): |Σ_Ω R_final| + B_assembly < 0.1·ΔF_c (成分 w, g, Q2, Q1, Q0; v は記録)。
  ΔF_c がその測定誤差以下なら UNDETERMINED。**注意 (plan との差)**: plan §5.1 #4g3 は ΔF を時間平均と新しい流束系列の check_quasisteady から
  作ると定めているが、この 0 step の記録が与えるのは再開スナップショット 1 枚の出口流束だけである。--dF で系列の平均差を、--dF-err で
  その測定誤差 (例: 窓平均の標準誤差) を外から与えられる。与えなければスナップショットの差と丸めの誤差だけを使い、判定に (snapshot) と付ける
  (確定の判定には使わない)。

#4pjg (--mask): 射影集合の節点番号 (extract_projection_mask.py が run_0592 の tp_update.h5 から作る npy) で、液・Q2・Q1・Q0 の残差を段
  (移流・拡散・凝縮ソース・ピン・その他・最終) ごとに分解し、和・中央値・正 (Q2 を増やす) / 負 (Q1 を減らす) の節点の割合、ソースの分岐と
  成長分岐の Sg クリップの作動を出す。段の寄与は記録した項 (移流・拡散・ソース) と写しの差 (ピン・その他) から作る。
"""
import argparse, csv, sys
import h5py
import numpy as np

U = 2.0 ** -24
U53 = 2.0 ** -53
FLT_MIN = 2.0 ** -126
LABELS = ["sp_zero", "sp_adv", "sp_diff", "chem", "sp_pin1", "cm_zero", "cm_adv", "tp_diff", "cond_src", "tracer", "pas_pin", "sp_pin2", "final"]
K_ADV_SP, K_ADV_PA, K_DIFF_OFF, K_DIFF_ON = 0, 1, 2, 3
S_TERM = {1: 0, 2: 1, 3: 2, 4: 3}   # 成分 g, Q2, Q1, Q0 → ソースの項スロット
S_BRANCH, S_CLIP, S_SGRAW, S_THETA, S_NROUND, S_J = 4, 5, 6, 7, 8, 9
CHECK2 = ["w", "g", "Q2", "Q1", "Q0"]


def gam(m, u):
    mu = m * u
    return np.where(mu < 1.0, mu / np.maximum(1.0 - mu, 1e-300), np.inf)


def dec(a):
    return [s.decode() if isinstance(s, bytes) else str(s) for s in a]


def is_species(c):
    return c == 0 or c >= 5


def stage_kinds(c):
    """成分 c の基点の写しと、各記録点までの段の種別 (adv / diff / src / pin / other)。"""
    if is_species(c):
        return "sp_zero", {"sp_adv": "adv", "sp_diff": "diff", "chem": "other", "sp_pin1": "pin", "cm_zero": "other", "cm_adv": "other",
                           "tp_diff": "diff", "cond_src": "other", "tracer": "other", "pas_pin": "other", "sp_pin2": "pin", "final": "other"}
    return "cm_zero", {"cm_adv": "adv", "tp_diff": "diff", "cond_src": "src", "tracer": "other", "pas_pin": "pin", "sp_pin2": "other", "final": "other"}


class Rec:
    def __init__(self, path, scale):
        f = h5py.File(path, "r")
        self.path = path
        self.attrs = {k: f.attrs[k] for k in f.attrs.keys()}
        self.comps = dec(self.attrs["components"])
        self.nC = len(self.comps)
        self.n = int(self.attrs["nNodes"]); self.nF = int(self.attrs["nFaces"])
        self.on = int(self.attrs["twophase_active"]) == 1
        self.coverage_ok = int(self.attrs.get("coverage_ok", 0)) == 1
        self.ic0 = f["/face/ic0"][...].astype(np.int64); self.ic1 = f["/face/ic1"][...].astype(np.int64)
        self.bc = f["/face/bcond"][...].astype(np.int64)
        self.adv = f["/face/adv"][...]; self.diff = f["/face/diff"][...]
        self.code = f["/face/code"][...]
        self.src = f["/source/slots"][...]
        self.labels = dec(f["/res/labels"][...])
        self.snap = f["/res/snap"][...]                       # [ns, nC, n] float
        self.state = f["/node/state"][...].astype(np.float64)  # [nC, n]
        self.V = f["/node/volume"][...].astype(np.float64)
        self.scale = scale
        self.x = f["/node/ccx"][...].astype(np.float64) * scale
        self.wd = f["/node/wall_dist"][...].astype(np.float64) * scale if "/node/wall_dist" in f else None
        self.pin = (f["/node/scalarDirichletPin"][...] == 1.0) if "/node/scalarDirichletPin" in f else np.zeros(self.n, bool)
        self.S = f["/node/condS_0"][...].astype(np.float64) if "/node/condS_0" in f else None
        self.bcKinds = dec(self.attrs["bcond_kinds"]); self.bcNames = dec(self.attrs["bcond_names"])
        f.close()
        self.problems = []
        if not self.coverage_ok:
            self.problems.append("forge reported a recording coverage mismatch (coverage_ok 0)")
        for l in LABELS:
            if self.labels.count(l) != 1:
                self.problems.append(f"snapshot '{l}' appears {self.labels.count(l)} times (expected once)")
        self._terms = {}

    def snapshot(self, label, c):
        return self.snap[self.labels.index(label), c].astype(np.float64)

    def diff_kernel(self):
        on = [k for k in (K_DIFF_OFF, K_DIFF_ON) if np.any(self.code[k] != -1)]
        if len(on) > 1:
            self.problems.append("both diffusion kernels recorded faces")
        return on[0] if on else None

    def end_masks(self):
        n = self.n
        r0 = (self.ic0 >= 0) & (self.ic0 < n); r1 = (self.ic1 >= 0) & (self.ic1 < n)
        return r0, r1

    def terms(self, c):
        """成分 c の節点ごとの項と加算の勘定 (キャッシュ)。v (c = 'v') は w − g。"""
        if c in self._terms:
            return self._terms[c]
        if c == "v":
            w, g = self.terms(0), self.terms(1)
            T = {k: w[k] - g[k] for k in ("adv", "diff", "src", "pin", "other", "R_final")}
            T["other_stages"] = {k: w["other_stages"].get(k, 0.0) - g["other_stages"].get(k, 0.0)
                                 for k in set(w["other_stages"]) | set(g["other_stages"])}
            T["m"] = w["m"] + g["m"]; T["A"] = w["A"] + g["A"]; T["nsub"] = w["nsub"] + g["nsub"]
            T["B"] = w["B"] + g["B"]; T["Bd"] = w["Bd"] + g["Bd"]
            T["stage_check"] = {}
            T["adv_face"] = w["adv_face"] - g["adv_face"]; T["diff_face"] = w["diff_face"] - g["diff_face"]
            T["adv_ok"] = w["adv_ok"] | g["adv_ok"]; T["diff_ok"] = w["diff_ok"] | g["diff_ok"]
            self._terms[c] = T
            return T
        n = self.n
        r0, r1 = self.end_masks()
        T = {k: np.zeros(n) for k in ("adv", "diff", "src")}
        m = np.zeros(n); A = np.zeros(n); nsub = np.zeros(n)

        def scatter(vals, ok, key):
            nonlocal m, A, nsub
            sub = (np.abs(vals) > 0) & (np.abs(vals) < FLT_MIN)
            for end, rr, sgn in ((self.ic0, r0, 1.0), (self.ic1, r1, -1.0)):
                sel = ok & rr
                idx = end[sel]
                T[key] += np.bincount(idx, weights=sgn * vals[sel], minlength=n)
                m += np.bincount(idx, minlength=n)
                A += np.bincount(idx, weights=np.abs(vals[sel]), minlength=n)
                nsub += np.bincount(idx, weights=sub[sel].astype(float), minlength=n)

        ka = K_ADV_PA if 1 <= c <= 4 else K_ADV_SP
        a = self.adv[c]
        adv_ok = np.isfinite(a) & (self.code[ka] != -1)
        if not adv_ok.any():
            self.problems.append(f"component {self.comps[c]}: no advection face recorded")
        scatter(np.where(adv_ok, a, 0.0), adv_ok, "adv")
        kd = self.diff_kernel()
        d = self.diff[c]
        diff_ok = (np.isfinite(d) & (self.code[kd] == 1)) if kd is not None else np.zeros(self.nF, bool)
        scatter(np.where(diff_ok, d, 0.0), diff_ok, "diff")
        if c in S_TERM:
            br = self.src[S_BRANCH]
            if not np.isfinite(br).all():
                self.problems.append("source branch slot has unwritten nodes")
            added = np.isin(np.nan_to_num(br, nan=-1) % 10, (1, 2))
            t = np.where(added, np.nan_to_num(self.src[S_TERM[c]]), 0.0)
            T["src"] = t
            nr = np.where(added, np.nan_to_num(self.src[S_NROUND]), 0.0)
            m += nr; A += np.abs(t)
            nsub += ((np.abs(t) > 0) & (np.abs(t) < FLT_MIN)).astype(float)
        # 写しの段
        base, kinds = stage_kinds(c)
        q0 = self.snapshot(base, c)
        if np.any(q0 != 0.0):
            self.problems.append(f"component {self.comps[c]}: residual at '{base}' is not zero ({int(np.count_nonzero(q0))} nodes)")
        pin = np.zeros(n); other = np.zeros(n); other_st = {}; stage_check = {}
        add_stage = {"adv": np.zeros(n), "diff": np.zeros(n), "src": np.zeros(n)}
        prev = q0
        for l in LABELS[LABELS.index(base) + 1:]:
            cur = self.snapshot(l, c); dlt = cur - prev
            kd_ = kinds[l]
            if kd_ == "pin":
                pin += dlt
            elif kd_ == "other":
                other += dlt; other_st[l] = dlt
            else:
                add_stage[kd_] += dlt
            prev = cur
        R = prev
        for k in ("adv", "diff", "src"):
            stage_check[k] = add_stage[k] - T[k]   # 段ごとの写しの差 − 記録した項 (段内の丸め)
        T.update(pin=pin, other=other, other_stages=other_st, R_final=R, m=m, A=A, nsub=nsub, stage_check=stage_check,
                 stage_snap=add_stage)
        # 上界 (節点)
        T["B"] = gam(m, U) * A + m * FLT_MIN
        T["Bd"] = gam(m + 6.0, U53) * (A + np.abs(R) + np.abs(pin) + np.abs(other))
        T["adv_face"] = np.where(adv_ok, a, 0.0); T["diff_face"] = np.where(diff_ok, d, 0.0)
        T["adv_ok"] = adv_ok; T["diff_ok"] = diff_ok
        self._terms[c] = T
        return T

    def outward(self, vals, ok, mask):
        """Ω の境界面の外向き流束 (和) と、境界種別ごとの内訳。"""
        n = self.n
        r0, r1 = self.end_masks()
        in0 = r0 & mask[np.clip(self.ic0, 0, n - 1)]
        in1 = r1 & mask[np.clip(self.ic1, 0, n - 1)]
        out = np.where(ok & in0 & ~in1, -vals, 0.0) + np.where(ok & in1 & ~in0, vals, 0.0)
        cats = {}
        internal = r0 & r1
        cats["internal (Omega cut)"] = float(np.sum(out[internal]))
        for b in np.unique(self.bc[~internal]):
            sel = (~internal) & (self.bc == b)
            nm = "boundary ?" if b < 0 else f"{self.bcKinds[b]}:{self.bcNames[b]}"
            cats[nm] = cats.get(nm, 0.0) + float(np.sum(out[sel]))
        nb = int(np.count_nonzero(ok & (in0 ^ in1)))
        return float(np.sum(out)), cats, float(np.sum(np.abs(out))), nb

    def exit_flux(self, c):
        """出口マーカーの外向き移流流束 (−a0 の和) と拡散 (0 のはず) の確認。"""
        T = self.terms(c)
        outlet = np.array([k.startswith("outlet") for k in self.bcKinds] + [False])
        bsel = np.where(self.bc >= 0, self.bc, len(self.bcKinds))
        sel = outlet[bsel]
        F = float(np.sum(-T["adv_face"][sel & T["adv_ok"]]))
        Fd = float(np.sum(T["diff_face"][sel & T["diff_ok"]]))
        err = float(gam(np.count_nonzero(sel) + 1.0, U53) * np.sum(np.abs(T["adv_face"][sel])))
        nf = int(np.count_nonzero(sel))
        return F, Fd, err, nf


def regions(rec):
    x = rec.x
    R = [("all", np.ones(rec.n, bool)),
         ("x[in,17)", x < 17.0), ("x[17,35)", (x >= 17.0) & (x < 35.0)),
         ("x[35,70)", (x >= 35.0) & (x < 70.0)), ("x[70,ex]", x >= 70.0)]
    if rec.wd is not None:
        band = (x >= 35.0) & (x < 70.0)
        R += [("x35-70 wd<0.4", band & (rec.wd < 0.4)), ("x35-70 wd>=0.4", band & (rec.wd >= 0.4))]
    R += [("pinned (ref)", rec.pin)]
    return R


def comp_list(rec):
    return list(range(rec.nC)) + ["v"]


def cname(rec, c):
    return "v" if c == "v" else rec.comps[c]


def region_budget(rec, c, mask):
    T = rec.terms(c)
    Fa, catA, absA, nbA = rec.outward(T["adv_face"], T["adv_ok"], mask)
    Fd, catD, absD, nbD = rec.outward(T["diff_face"], T["diff_ok"], mask)
    src = T["src"][mask]
    Sp = float(np.sum(np.maximum(src, 0.0))); Sm = float(np.sum(np.minimum(src, 0.0)))
    Bpin = float(np.sum(T["pin"][mask])); Ro = float(np.sum(T["other"][mask]))
    R = float(np.sum(T["R_final"][mask]))
    E = R - (-Fa - Fd + Sp + Sm + Bpin + Ro)
    Enode = float(np.sum((T["R_final"] - T["adv"] - T["diff"] - T["src"] - T["pin"] - T["other"])[mask]))
    Nn = int(mask.sum())
    Bas = float(np.sum(T["B"][mask]))
    Bdbl = float(np.sum(T["Bd"][mask])) + float(gam(Nn + nbA + nbD + 10.0, U53)) * (
        float(np.sum(T["A"][mask])) + absA + absD + abs(Sp) + abs(Sm) + abs(Bpin) + abs(Ro) + float(np.sum(np.abs(T["R_final"][mask]))))
    # 面の積分と節点の和の照合 (内部面の相殺; double の丸めだけのはず)
    cancel = float(np.sum((T["adv"] + T["diff"])[mask])) - (-Fa - Fd)
    return dict(R=R, Fadv=Fa, Fdiff=Fd, Sp=Sp, Sm=Sm, Bpin=Bpin, Rother=Ro, E=E, Enode=Enode, Bas=Bas, Bdbl=Bdbl,
                B=Bas + Bdbl, Bunder=float(np.sum((T["m"] * FLT_MIN)[mask])), nsub=int(np.sum(T["nsub"][mask])),
                madd=int(np.sum(T["m"][mask])), cancel=cancel, catA=catA, catD=catD, N=Nn,
                other_st={k: float(np.sum(v[mask])) for k, v in T["other_stages"].items()})


def fmt(v):
    return f"{v:+.4e}"


def operator_budget(rec, rows):
    tag = "ON" if rec.on else "OFF"
    print(f"\n=== {rec.path} [{tag}] G3-a operator budget (one assembly, no update) ===")
    all_ok = True
    worst = []
    other_gap = []
    for rname, mask in regions(rec):
        if not mask.any():
            print(f"\n  -- Omega = {rname}: no nodes"); continue
        cols = comp_list(rec)
        res = {c: region_budget(rec, c, mask) for c in cols}
        print(f"\n  -- Omega = {rname} ({int(mask.sum())} nodes) --")
        print("  " + f"{'term':34s}" + "".join(f"{cname(rec, c):>14s}" for c in cols))
        for key, lab in [("R", "sum R_final"), ("Fadv", "F_adv (outward)"), ("Fdiff", "F_diff (outward)"), ("Sp", "S+"), ("Sm", "S-"),
                         ("Bpin", "B_pin"), ("Rother", "R_other"), ("E", "E_assembly"), ("Enode", "  (node-wise sum of E_i)"),
                         ("Bas", "  sum gamma_m*sum|a| + underflow"), ("Bdbl", "  B_double"), ("B", "B_assembly")]:
            print("  " + f"{lab:34s}" + "".join(f"{fmt(res[c][key]):>14s}" for c in cols))
            for c in cols:
                rows.append([rec.path, tag, rname, cname(rec, c), key, res[c][key]])
        print("  " + f"{'  adds / subnormal terms':34s}" + "".join(f"{str(res[c]['madd'])+'/'+str(res[c]['nsub']):>14s}" for c in cols))
        print("  " + f"{'  internal-face cancel (double)':34s}" + "".join(f"{fmt(res[c]['cancel']):>14s}" for c in cols))
        okr = []
        for c in cols:
            r = res[c]; ok = abs(r["E"]) <= r["B"]
            okr.append(ok); all_ok &= ok
            worst.append((abs(r["E"]) / r["B"] if r["B"] > 0 else (np.inf if r["E"] != 0 else 0.0), rname, cname(rec, c)))
            rows.append([rec.path, tag, rname, cname(rec, c), "check1", "PASS" if ok else "FAIL"])
        print("  " + f"{'(1) |E| <= B_assembly':34s}" + "".join(f"{('PASS' if o else 'FAIL'):>14s}" for o in okr))
        nz = {k: [res[c]["other_st"].get(k, 0.0) for c in cols] for k in sorted({k for c in cols for k in res[c]["other_st"]})}
        for k, vals in nz.items():
            if any(v != 0.0 for v in vals):
                print("  " + f"{'  ** R_other stage ' + k:34s}" + "".join(f"{fmt(v):>14s}" for v in vals)
                      + "   (unrecorded writer between recording points)")
                # 命名した段は「この成分に触れない」段だけなので、0 でなければ記録の欠落 (E に入れずに FAIL とする)
                all_ok = False; other_gap.append((rname, k))
        if rname == "all":
            print("  outward advective flux by boundary (Omega = all):")
            keys = sorted({k for c in cols for k in res[c]["catA"]})
            for k in keys:
                print("  " + f"    {k:30s}" + "".join(f"{fmt(res[c]['catA'].get(k, 0.0)):>14s}" for c in cols))
            keys = sorted({k for c in cols for k in res[c]["catD"]})
            print("  outward diffusive flux by boundary (Omega = all; node half-faces are not added):")
            for k in keys:
                print("  " + f"    {k:30s}" + "".join(f"{fmt(res[c]['catD'].get(k, 0.0)):>14s}" for c in cols))
    # 段ごとの照合 (写しの差 − 記録した項; 段内の float の丸めだけのはず)
    print("\n  per-stage check over all nodes: sum|snapshot stage increment - recorded terms| (float rounding inside that stage only; "
          "a value near sum|terms| means a stage wrote without being recorded)")
    for c in range(rec.nC):
        T = rec.terms(c)
        parts = []
        for k in ("adv", "diff", "src"):
            dv = float(np.sum(np.abs(T["stage_check"][k])))
            parts.append(f"{k} {dv:.3e}")
        print(f"    {rec.comps[c]:>8s}: " + ", ".join(parts))
    w = max(worst) if worst else (0.0, "", "")
    print(f"\n  G3-a (1) [{tag}]: {'PASS' if all_ok else 'FAIL'}  (worst |E|/B = {w[0]:.3e} at {w[1]} / {w[2]})"
          + (f"; nonzero R_other stages (unrecorded writers): {sorted(set(k for _r, k in other_gap))}" if other_gap else ""))
    return all_ok


def physical_budget(on, off, dF_user, dFerr_user, rows):
    print("\n=== G3-a (2) physical budget: |sum R_final| + B_assembly < 0.1 dF_c (dF from the outlet-marker outward advective flux) ===")
    dF = {}; dFerr = {}; src = {}
    for c in CHECK2 + ["v"]:
        ci = "v" if c == "v" else on.comps.index(c)
        cj = "v" if c == "v" else off.comps.index(c)
        if ci == "v":
            Fw1, Fdw1, e1w, n1 = on.exit_flux(0); Fg1, Fdg1, e1g, _ = on.exit_flux(1)
            Fw0, Fdw0, e0w, n0 = off.exit_flux(0); Fg0, Fdg0, e0g, _ = off.exit_flux(1)
            F1, Fd1, e1 = Fw1 - Fg1, Fdw1 - Fdg1, e1w + e1g
            F0, Fd0, e0 = Fw0 - Fg0, Fdw0 - Fdg0, e0w + e0g
        else:
            F1, Fd1, e1, n1 = on.exit_flux(ci)
            F0, Fd0, e0, n0 = off.exit_flux(cj)
        dF[c] = abs(F1 - F0); dFerr[c] = e1 + e0; src[c] = "snapshot"
        if c in dF_user:
            dF[c] = dF_user[c]; src[c] = "user (--dF)"
        if c in dFerr_user:
            dFerr[c] = dFerr_user[c]; src[c] += " + user error"
        print(f"  {c:3s}: F_exit ON {F1:+.6e} (diffusion on outlet half-faces {Fd1:+.1e}), OFF {F0:+.6e} ({Fd0:+.1e}); "
              f"dF = {dF[c]:.6e} [{src[c]}], measurement error {dFerr[c]:.3e}; outlet half-faces {n1}/{n0}")
        rows.append(["ON-OFF", "", "outlet", c, "dF", dF[c]])
        if Fd1 != 0.0 or Fd0 != 0.0:
            print(f"       ** nonzero diffusion on outlet half-faces (expected 0 on this node path) **")
    verdict = {}
    for rec in (on, off):
        tag = "ON" if rec.on else "OFF"
        print(f"\n  [{tag}] {rec.path}")
        print("  " + f"{'Omega':18s}" + "".join(f"{c:>34s}" for c in CHECK2 + ["v"]))
        for rname, mask in regions(rec):
            if not mask.any():
                continue
            cells = []
            for c in CHECK2 + ["v"]:
                ci = "v" if c == "v" else rec.comps.index(c)
                r = region_budget(rec, ci, mask)
                lhs = abs(r["R"]) + r["B"]; rhs = 0.1 * dF[c]
                if dF[c] <= dFerr[c]:
                    v = "UNDETERMINED"
                else:
                    v = "PASS" if lhs < rhs else "FAIL"
                if src[c] == "snapshot" and v != "UNDETERMINED":
                    v += "(snapshot)"
                cells.append(f"{lhs:.2e}/{rhs:.2e} {v}")
                if c != "v":
                    verdict.setdefault((tag, rname), []).append(v)
                rows.append([rec.path, tag, rname, c, "check2", f"{lhs}/{rhs} {v}"])
            print("  " + f"{rname:18s}" + "".join(f"{s:>34s}" for s in cells))
    print("  (cells: |sum R_final| + B_assembly / 0.1 dF; v is recorded, not judged)")
    return verdict


def pjg(rec, idx, rows):
    tag = "ON" if rec.on else "OFF"
    if idx.size == 0:
        print(f"\n=== #4pjg [{tag}]: empty mask ==="); return
    if idx.max() >= rec.n or idx.min() < 0:
        raise SystemExit(f"mask indices out of range for {rec.path} ({rec.n} nodes)")
    print(f"\n=== #4pjg [{tag}] {rec.path}: residual decomposition on the projection set ({idx.size} nodes; "
          f"x [{rec.x[idx].min():.3f}, {rec.x[idx].max():.3f}] scaled; node-coordinate checksum sum(ccx) = {float(np.sum(rec.x)/rec.scale):.10e}, "
          f"compare with extract_projection_mask.py) ===")
    stages = ["adv", "diff", "src", "pin", "other", "R_final"]
    V = rec.V[idx]
    for c in (1, 2, 3, 4):
        T = rec.terms(c); nm = rec.comps[c]
        q = rec.state[c][idx]
        print(f"\n  {nm}: stage contribution to the residual R (R > 0 increases {nm}); relative = R/(V q) [1/s] where q > 0")
        print("  " + f"{'stage':10s}{'sum':>14s}{'median':>14s}{'frac>0':>9s}{'frac<0':>9s}{'med rel':>13s}{'sum|.|':>12s}")
        for st in stages:
            v = T[st][idx]
            with np.errstate(all="ignore"):
                rel = np.where(q > 0, v / (V * q), np.nan)
            mr = np.nanmedian(rel) if np.isfinite(rel).any() else np.nan
            print("  " + f"{st:10s}{np.sum(v):>+14.4e}{np.median(v):>+14.4e}{np.mean(v > 0):>9.3f}{np.mean(v < 0):>9.3f}{mr:>+13.4e}{np.sum(np.abs(v)):>12.3e}")
            rows.append([rec.path, tag, "mask", nm, f"pjg:{st}:sum", float(np.sum(v))])
            rows.append([rec.path, tag, "mask", nm, f"pjg:{st}:median", float(np.median(v))])
        # 段の写しと記録した項の照合
        chk = {k: float(np.sum(np.abs(T["stage_check"][k][idx]))) for k in ("adv", "diff", "src")}
        print("  " + "  snapshot-vs-recorded per stage sum|d|: " + ", ".join(f"{k} {v:.2e}" for k, v in chk.items()))
    # 向き: Q2 を増やす段 / Q1 を減らす段
    TQ2, TQ1 = rec.terms(2), rec.terms(3)
    up = TQ2["R_final"][idx] > 0; dn = TQ1["R_final"][idx] < 0
    print(f"\n  final residual: Q2 up (R>0) {int(up.sum())}/{idx.size}, Q1 down (R<0) {int(dn.sum())}/{idx.size}, both {int((up & dn).sum())}")
    print("  " + f"{'stage':10s}{'Q2 up | among R_Q2>0':>24s}{'Q1 down | among R_Q1<0':>26s}{'Q2 up (all)':>14s}{'Q1 down (all)':>15s}")
    for st in ["adv", "diff", "src", "pin", "other"]:
        a2 = TQ2[st][idx]; a1 = TQ1[st][idx]
        f2 = np.mean(a2[up] > 0) if up.any() else np.nan; f1 = np.mean(a1[dn] < 0) if dn.any() else np.nan
        print("  " + f"{st:10s}{f2:>24.3f}{f1:>26.3f}{np.mean(a2 > 0):>14.3f}{np.mean(a1 < 0):>15.3f}")
        rows.append([rec.path, tag, "mask", "Q2", f"pjg:{st}:frac_up", float(np.mean(a2 > 0))])
        rows.append([rec.path, tag, "mask", "Q1", f"pjg:{st}:frac_down", float(np.mean(a1 < 0))])
    # ソースの分岐と成長分岐のクリップ
    br = rec.src[S_BRANCH][idx]
    names = {0: "none", 1: "growth", 2: "evaporation", 10: "none (double)", 11: "growth (double)", 12: "evaporation (double)"}
    cnt = {names.get(int(b), str(b)): int(np.count_nonzero(br == b)) for b in np.unique(br[np.isfinite(br)])}
    nnan = int(np.count_nonzero(~np.isfinite(br)))
    clip = rec.src[S_CLIP][idx]; sgr = rec.src[S_SGRAW][idx]; th = rec.src[S_THETA][idx]
    print(f"\n  source branch on the set: {cnt}" + (f", not written {nnan}" if nnan else ""))
    print(f"  growth-branch clip 'if (Sg < 0) Sg = 0' fired: {int(np.nansum(clip == 1.0))} nodes; "
          f"unclipped Sg < 0: {int(np.count_nonzero(np.isfinite(sgr) & (sgr < 0)))}; theta < 1: {int(np.count_nonzero(np.isfinite(th) & (th < 1)))}")
    if np.isfinite(sgr).any():
        print(f"  unclipped Sg (growth branch): min {np.nanmin(sgr):+.3e}, median {np.nanmedian(sgr):+.3e}")
    if rec.S is not None:
        S = rec.S[idx]
        print(f"  saturation S = p_v/p_sat on the set: median {np.median(S):.4g}, S <= 1 (unsaturated) {np.mean(S <= 1.0):.3f}")
    for b, k in cnt.items():
        rows.append([rec.path, tag, "mask", "", f"pjg:branch:{b}", k])


def parse_kv(s):
    out = {}
    if not s:
        return out
    for kv in s.split(","):
        k, v = kv.split("=")
        out[k.strip()] = float(v)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("on_h5"); ap.add_argument("off_h5")
    ap.add_argument("--length-scale", type=float, default=1000.0, help="ccx・wall_dist に掛ける倍率 (既定 1000 = m → mm)")
    ap.add_argument("--dF", default="", help="成分=値,... 時間平均の出口流束差 |ΔF_c| を外から与える (既定はスナップショットの差)")
    ap.add_argument("--dF-err", default="", help="成分=値,... ΔF_c の測定誤差 (例: 窓平均の標準誤差)")
    ap.add_argument("--mask", default=None, help="#4pjg: 射影集合の節点番号 (npy; extract_projection_mask.py)")
    ap.add_argument("--mask-arm", choices=["on", "off", "both"], default="on")
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()
    on = Rec(a.on_h5, a.length_scale); off = Rec(a.off_h5, a.length_scale)
    rows = []
    for r in (on, off):
        print(f"# {r.path}: {r.n} nodes, {r.nF} faces, components {r.comps}, twophase_active={int(r.on)}, "
              f"x range [{r.x.min():.3f}, {r.x.max():.3f}] (scaled; mm expected)" + ("" if r.wd is None else f", wall_dist max {r.wd.max():.3f}"))
        print(f"#   input_value_file = {r.attrs.get('input_value_file')}, coverage_ok = {int(r.coverage_ok)}")
        for l in r.attrs.get("summary", []):
            print("#   " + (l.decode() if isinstance(l, bytes) else str(l)))
    invalid = []
    if not on.on or off.on:
        invalid.append("first file must be the two-phase ON record and the second the OFF record")
    if on.n != off.n or on.nF != off.nF or not np.array_equal(on.x, off.x) or not np.array_equal(on.ic0, off.ic0):
        invalid.append("ON and OFF records are not on the same mesh / face ordering")
    for r in (on, off):
        for c in range(r.nC):
            r.terms(c)
        invalid += [f"{r.path}: {p}" for p in r.problems]
    if invalid:
        print("\nG3-a INVALID:")
        for p in invalid:
            print("  - " + p)
        return 2
    ok1 = [operator_budget(r, rows) for r in (on, off)]
    v2 = physical_budget(on, off, parse_kv(a.dF), parse_kv(a.dF_err), rows)
    print(f"\nG3-a (1) assembly closure: ON {'PASS' if ok1[0] else 'FAIL'}, OFF {'PASS' if ok1[1] else 'FAIL'}")
    for (tag, rname), vs in v2.items():
        if rname == "all" or rname.startswith("x35-70"):
            print(f"G3-a (2) [{tag}] {rname}: " + ", ".join(f"{c} {v}" for c, v in zip(CHECK2, vs)))
    if a.mask:
        idx = np.load(a.mask).astype(np.int64)
        arms = {"on": [on], "off": [off], "both": [on, off]}[a.mask_arm]
        for r in arms:
            pjg(r, idx, rows)
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["file", "arm", "omega", "component", "quantity", "value"]); w.writerows(rows)
        print(f"\nwrote {a.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
