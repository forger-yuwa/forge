#!/usr/bin/env python3
"""`slauWallNormalChi` の適用診断が**できる設定か**を検査し、できるなら面ごとの表を出す
(plan convection-slau-wall-normal-chi-usage-rule §4.1 / §4.1.1、§5.1 #3b)。

`diag_wall_cv_budget.py` は 1 次 SLAU の再計算で、実カーネルの `slauContactFloor` の質量流束追加
(`convectiveFlux_slau_d.inc.cuh:620`) や `sstEnergyIncludesK` の p* 差 (:516) を持たない。設定によっては
流束の符号まで食い違う (codex plan-3 M3 の再計算: +0.099 対 −0.191)。そこで

  1. **キー検査**: 実効値は**起動ログのエコー > YAML** の順で解決する (ログ欠落は診断不能、両者の食い違いも診断不能)。
     エコーを出さないキーは「未確認 (YAML/既定値のみ)」として VERDICT に列挙し、「両方から確認」とは書かない。
     軸対称は新旧両配置 (`mesh.isAxisymmetric` と旧 `physProp.isAxisymmetric`。`solverConfig.cpp:881` が受理) を、
     周期は `bcondConfig.yaml` の kind (`periodic` / `cyclic` を含む) を拒否する。
  2. **面流束の照合**: 同じ run を 1 step、`FORGE_DUMP_MASSFLUX=<path>` で回した `massflux` (カーネル値) と
     `<path>.state` (カーネルが読んだ状態: ro,Ux,Uy,Uz,P,sonic) から、対象 CV の**全接続面**を照合する:
       - 内部面: |mdot_kernel − mdot_tool| ≤ 2 τ_b (τ_b = 1e-5 A ρ̄ ĉ、前 plan V5)
       - 境界半割面 (nei < 0): 壁種別 (bcond の kind) を記録し、|mdot_kernel| ≤ 2 τ_b (Dirichlet 壁の期待値は 0)
     状態・面積・許容差・誤差の有限性、正の密度・音速、CV 番号の範囲、照合面数 > 0 を必須とする。
  3. 1・2 とも通ったときだけ「診断可能」。面ごとに chi, chi_n, ΔP, mdot(chi), mdot(chi_n) と、
     CV への正味補充の変化 Σ[流出(chi_n) − 流出(chi)] (< 0 で補充が増える、§4.1 条件 iii) を出す。
     ρ_i = 対象壁 CV の**壁でない隣接ノード全部の平均**。

判定は「診断可能 / 診断不能 (理由)」。診断不能は「flag 1 が不要」の意味ではない。
**これが示すのは 1 dump での診断可能性と条件 (iii)**。条件 (i)(ii) (3 dump 以上の持続) は `diag_wall_cv_budget.py` で別に見る。

usage:
  python3 diag_applicability.py RUN_DIR --mesh RUN_DIR/sern.h5 --massflux mf.bin --ids 153797,153880 \
      --wall-phys-ids 1,2,3,4,10,11,12,13,15 [--config-only]
試験: test_diag_applicability.py
"""
import argparse, re, sys
from pathlib import Path
import numpy as np, yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from diag_wall_cv_budget import slau_mdot, parse_struct

# (表示名, yaml パス, 要求集合, 省略時の既定値, 理由)  — plan §4.1.1 の表
REQ = [
    ("mesh.discretization", ("mesh", "discretization"), {"node"}, None, "対象構成"),
    ("mesh.nodeWallDirichlet", ("mesh", "nodeWallDirichlet"), {"1"}, "0", "対象構成"),
    ("solver", ("solver",), {"SLAU", "SLAU2"}, None, "質量流束式が同じ (SLAU2 は圧力束のみ差)"),
    ("space.convMethod", ("space", "convMethod"), {"0"}, None, "ツールは再構成なし"),
    ("space.slauContactFloor", ("space", "slauContactFloor"), {"0"}, "0", ":620 の追加項がツールに無い"),
    ("turbulence.sstEnergyIncludesK", ("turbulence", "sstEnergyIncludesK"), {"0"}, "0", ":516 の p* 差"),
    ("space.lowMachPrecond", ("space", "lowMachPrecond"), {"0"}, "0", "散逸スケール c'"),
    ("time.deltaT.lowMachPrecond", ("time", "deltaT", "lowMachPrecond"), {"0"}, "0", "散逸スケール c' (deltaT 側の綴り)"),
    ("space.slauWallNormalChi", ("space", "slauWallNormalChi"), {"0"}, "0", "--wall-normal-chi は置換予測用"),
    ("space.badReconFallback", ("space", "badReconFallback"), {"0"}, "0", "面状態の差し替え"),
    ("mesh.isAxisymmetric", ("mesh", "isAxisymmetric"), {"0"}, "0", "適用範囲外 (軸対称)"),
    ("physProp.isAxisymmetric", ("physProp", "isAxisymmetric"), {"0"}, "0", "適用範囲外 (軸対称、旧配置)"),
]
RECORD = [("physProp.thermalMethod", ("physProp", "thermalMethod"), "ĉ・T (CPG/TP をツールの EOS 経路と合わせる)")]
ABSENT = [("condensation", "凝縮・二相 (面エンタルピー・P=ρ(1−g)RT)")]


def yget(doc, path):
    cur = doc
    for k in path:
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
    return None if isinstance(cur, (dict, list)) else str(cur).strip().strip('"\'')


def log_echo(log_text):
    """`'key' in 'section': value` / `'key': value` の起動エコーを {section.key: value} で返す。"""
    out = {}
    for key, sec, val in re.findall(r"^'(\w+)' in '([\w.]+)':\s*(\S+)", log_text, re.M):
        out[f"{sec}.{key}"] = val
    for k, v in re.findall(r"^'(\w+)':\s*(\S+)", log_text, re.M):
        out.setdefault(k, v)
    return out


def _norm(name, v):
    if v is None:
        return None
    if name == "solver":
        return v.upper()
    try:
        f = float(v)
        return str(int(f)) if f == int(f) else repr(f)
    except ValueError:
        return v


def check_config(cfg, echo, bcond, log_present):
    """純関数: (rows, bad, unconfirmed)。cfg/bcond は dict、echo は log_echo の結果。"""
    rows, bad, unconf = [], [], []
    if not log_present:
        bad.append("forge_run.log が無い (起動エコーで実効値を確認できない)")
    for name, path, want, dflt, why in REQ:
        yv = _norm(name, yget(cfg, path))
        ev = _norm(name, echo.get(name) if name in echo else (echo.get(path[0]) if len(path) == 1 else None))
        eff = ev if ev is not None else (yv if yv is not None else _norm(name, dflt))
        ok = eff in want
        note = why
        if ev is not None and yv is not None and ev != yv:
            ok, note = False, f"yaml と起動エコーが食い違う ({yv} / {ev})"
        if ev is None:
            unconf.append(name)
        rows.append((name, yv, ev, eff, sorted(want), ok, note))
        if not ok:
            bad.append(f"{name}={eff} (要求 {sorted(want)}; {note})")
    for name, path, why in RECORD:
        yv = _norm(name, yget(cfg, path)); ev = _norm(name, echo.get(name))
        ok = not (ev is not None and yv is not None and ev != yv)
        if ev is None:
            unconf.append(name)
        rows.append((name, yv, ev, ev or yv, ["記録のみ"], ok, why if ok else f"yaml と起動エコーが食い違う ({yv} / {ev})"))
        if not ok:
            bad.append(f"{name}: yaml と起動エコーが食い違う ({yv} / {ev})")
    for key, why in ABSENT:
        if isinstance(cfg.get(key), dict) and cfg.get(key):
            bad.append(f"{key} が有効 ({why})")
    for bname, b in (bcond or {}).items():
        kind = str((b or {}).get("kind", "")).lower() if isinstance(b, dict) else ""
        if "periodic" in kind or "cyclic" in kind:
            bad.append(f"bcond '{bname}' が周期 (kind {kind}; 適用範囲外)")
    return rows, bad, unconf


def check_faces(own, nei, S, mf, st, ids, wall, bface_kind):
    """純関数: 対象 CV の全接続面を照合する。

    own/nei: 面の owner/neighbor (nei<0 は境界半割面)、S: 面ベクトル (nF,3)、mf: カーネル massflux (nF)、
    st: dict ro,Ux,Uy,Uz,P,sonic (nC)、ids: 対象 CV、wall: 壁ノード集合、bface_kind: {面 id: bcond kind}。
    戻り: (ok, bad, per_cv, summary)"""
    bad, per_cv = [], []
    nF, nC = len(own), len(st["ro"])
    if mf.size != nF:
        return False, [f"massflux {mf.size} 面 != メッシュ {nF} 面 (欠落面)"], [], {}
    n_int = n_bnd = n_nonfin = 0; worst = 0.0
    for cv in ids:
        if not (0 <= cv < nC):
            bad.append(f"CV {cv} が範囲外 (0..{nC - 1})"); continue
        faces = np.where((own == cv) | (nei == cv))[0]
        rows, d_out = [], 0.0
        nbr_nonwall = []
        for k in faces:
            A = float(np.linalg.norm(S[k]))
            if not (np.isfinite(A) and A > 0):
                bad.append(f"CV {cv} 面 {k}: 面積が不正 ({A})"); n_nonfin += 1; continue
            n = S[k] / A
            if nei[k] < 0:                                       # 境界半割面
                kind = bface_kind.get(int(k), "?")
                r = abs(mf[k]) / (2 * 1e-5 * A * st["ro"][cv] * st["sonic"][cv])
                if not np.isfinite(r):
                    bad.append(f"CV {cv} 境界面 {k}: 非有限"); n_nonfin += 1; continue
                n_bnd += 1; worst = max(worst, r)
                if r > 1.0:
                    bad.append(f"CV {cv} 境界面 {k} ({kind}): |mdot_kernel| {mf[k]:.3e} が 2τ_b を超える ({r:.2f}×)")
                rows.append(("bnd", int(k), -1, kind, float(mf[k]), 0.0, r, None, None, None))
                continue
            j = int(nei[k]) if own[k] == cv else int(own[k])
            L = {q: st[q][own[k]] for q in st}; R = {q: st[q][nei[k]] for q in st}
            vals = [L[q] for q in ("ro", "sonic")] + [R[q] for q in ("ro", "sonic")]
            if not all(np.isfinite(list(L.values()) + list(R.values()))) or min(vals) <= 0 or not np.isfinite(mf[k]):
                bad.append(f"CV {cv} 面 {k}: 状態が非有限または ρ/c ≤ 0"); n_nonfin += 1; continue
            m0, _, _, chi, ch, _, _ = slau_mdot(A, *n, L, R, wall_face=False)
            is_wf = (own[k] in wall) or (nei[k] in wall)
            m1, _, _, chin, _, _, _ = slau_mdot(A, *n, L, R, wall_face=is_wf)
            tau = 1e-5 * A * 0.5 * (L["ro"] + R["ro"]) * ch
            r = abs(mf[k] - m0) / (2 * tau)
            if not np.isfinite(r):
                bad.append(f"CV {cv} 面 {k}: 誤差が非有限"); n_nonfin += 1; continue
            n_int += 1; worst = max(worst, r)
            if r > 1.0:
                bad.append(f"CV {cv} 面 {k}: |mdot_kernel − mdot_tool| が 2τ_b を超える ({r:.2f}×)")
            if j not in wall:
                nbr_nonwall.append(j)
            sgn = 1.0 if own[k] == cv else -1.0
            d_out += sgn * (m1 - m0)
            rows.append(("int", int(k), j, None, float(mf[k]), float(m0), r, chi, chin, (R["P"] - L["P"], m1)))
        rho_i = float(np.mean(st["ro"][nbr_nonwall])) if nbr_nonwall else float("nan")
        per_cv.append({"cv": cv, "faces": len(faces), "rows": rows, "rho_w": float(st["ro"][cv]),
                       "rho_i": rho_i, "d_out": d_out})
    if n_int + n_bnd == 0:
        bad.append("照合面が 0 面")
    summary = {"n_int": n_int, "n_bnd": n_bnd, "n_nonfinite": n_nonfin, "worst": worst}
    return (not bad), bad, per_cv, summary


def main():
    import h5py
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--mesh", default=None)
    ap.add_argument("--massflux", default=None, help="FORGE_DUMP_MASSFLUX の出力 (同名 .state も読む)")
    ap.add_argument("--ids", default="", help="対象壁 CV (節点) をカンマ区切り")
    ap.add_argument("--wall-phys-ids", default="", help="壁 physID (bcondConfig の wall 系) をカンマ区切り")
    ap.add_argument("--config-only", action="store_true", help="キー検査だけ行う")
    a = ap.parse_args()
    rd = Path(a.run_dir)
    cfg = yaml.safe_load((rd / "solverConfig.yaml").read_text()) or {}
    bcp = rd / "bcondConfig.yaml"
    bcond = yaml.safe_load(bcp.read_text()) if bcp.exists() else {}
    logp = rd / "forge_run.log"
    echo = log_echo(logp.read_text(errors="replace")) if logp.exists() else {}
    rows, bad, unconf = check_config(cfg, echo, bcond, logp.exists())
    print("=== キー検査 (plan §4.1.1; 実効値 = 起動エコー > YAML > 既定値) ===")
    for name, yv, ev, eff, want, ok, why in rows:
        print(f"  {'OK ' if ok else 'NG '} {name:30s} yaml={yv!s:6s} echo={ev!s:6s} 実効={eff!s:6s} 要求={want}  ({why})")
    if unconf:
        print("  未確認 (起動エコーに無く YAML/既定値のみ): " + ", ".join(unconf))
    if bad:
        print("VERDICT: 診断不能 (キー検査)\n  " + "\n  ".join(bad)); return 2
    if a.config_only:
        print("VERDICT: キー検査のみ OK (面流束の照合は未実施 — これだけでは「診断可能」と言わない)"); return 0
    if not (a.mesh and a.massflux and a.ids and a.wall_phys_ids):
        sys.exit("--mesh --massflux --ids --wall-phys-ids が要る (または --config-only)")
    ids = [int(s) for s in a.ids.split(",") if s]
    kind_of = {}
    for bname, b in (bcond or {}).items():
        if isinstance(b, dict) and "physID" in b:
            kind_of[str(b["physID"])] = f"{bname}:{b.get('kind')}"
    with h5py.File(a.mesh) as f:
        S = np.asarray(f["PLANES/surfVect"], np.float64).reshape(-1, 3)
        own, nei, _ = parse_struct(np.asarray(f["PLANES/STRUCT"]), len(S))
        wall, bface_kind = set(), {}
        for pid in f["BCONDS"]:
            ip = np.asarray(f["BCONDS"][pid]["iPlanes"]).ravel()
            for k in ip[ip >= 0]:
                bface_kind[int(k)] = kind_of.get(pid, f"physID {pid}")
        for pid in [s for s in a.wall_phys_ids.split(",") if s]:
            if pid in f["BCONDS"]:
                c = np.asarray(f["BCONDS"][pid]["iCells"]).ravel(); wall.update(c[c >= 0].tolist())
    mf = np.fromfile(a.massflux, dtype=np.float32).astype(np.float64)
    raw = np.fromfile(a.massflux + ".state", dtype=np.float32).astype(np.float64)
    if raw.size % 6:
        print("VERDICT: 診断不能 (.state の長さが 6 の倍数でない)"); return 3
    nC = raw.size // 6
    st = dict(zip(("ro", "Ux", "Uy", "Uz", "P", "sonic"), raw.reshape(6, nC)))
    ok, fbad, per_cv, sm = check_faces(own.astype(int), nei.astype(int), S, mf, st, ids, wall, bface_kind)
    for c in per_cv:
        print(f"\n=== CV {c['cv']}  接続面 {c['faces']}  ρ_w {c['rho_w']:.5e}  ρ_i (非壁隣接の平均) {c['rho_i']:.5e}  ρ_w/ρ_i {c['rho_w']/c['rho_i']:.4f} ===")
        print(f"  {'種別':>4} {'face':>9} {'nbr':>9} {'chi':>7} {'chi_n':>7} {'dP':>10} {'mdot_k':>12} {'mdot_tool':>12} {'|d|/2τb':>8} {'mdot(chi_n)':>12}")
        for kind, k, j, bk, mk, mt, r, chi, chin, extra in c["rows"]:
            if kind == "bnd":
                print(f"  {'境界':>4} {k:9d} {'-':>9} {'':>7} {'':>7} {'':>10} {mk:12.4e} {0.0:12.4e} {r:8.3f}   ({bk})")
            else:
                print(f"  {'内部':>4} {k:9d} {j:9d} {chi:7.4f} {chin:7.4f} {extra[0]:10.3e} {mk:12.4e} {mt:12.4e} {r:8.3f} {extra[1]:12.4e}")
        print(f"  Σ[流出(chi_n) − 流出(chi)] = {c['d_out']:+.4e} kg/s  ({'補充が増える' if c['d_out'] < 0 else '補充は増えない'})")
    print(f"\n=== 面流束の照合: 内部面 {sm.get('n_int')}・境界半割面 {sm.get('n_bnd')}・非有限 {sm.get('n_nonfinite')}、"
          f"max |誤差|/(2τ_b) = {sm.get('worst', float('nan')):.3f} ===")
    if unconf:
        print("  未確認のキー (起動エコーに無い): " + ", ".join(unconf))
    if not ok:
        print("VERDICT: 診断不能 (面流束の照合)\n  " + "\n  ".join(fbad[:20])); return 3
    print("VERDICT: 診断可能 (キー検査 OK・全接続面が 2τ_b 以内)。示したのは 1 dump での診断可能性と条件 (iii)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
