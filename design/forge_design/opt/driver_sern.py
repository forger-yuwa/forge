"""⑤ SERN の多作動点 MOO ドライバ (S6)。

dv = (M_c, f, theta_r0_deg, theta_c0_deg, L_cowl) [YAML の dv.min/max]。1 評価 = 作動点数分の forge run
(`spec.operating_points[]`、重み w_k)。目的 (最小化):
    f1 = −Σ_k w_k C_T^(k)      (Kriging。RANS では摩擦込み C_T)
    f2 = L_ramp / H            (逆設計の決定的結果。候補点では Kriging で近似)
ゲート (plan §5.1 R1, `metrics/sern_gates.py`): 逆設計成立 (kernel 内に key point、C⁻ が直線ランプに着地、
`L_ramp_max` は最終輪郭で再検査)、メッシュ品質 PASS、各作動点で **rc == 0 かつ 保存場が有限・正値、全残差に NaN/rising 無し、
実目的量 (RANS は C_T_with_shear) と C_T/C_L/C_M が STEADY**。発散した評価の採用 (旧 §4.13-3) は撤回した:
ゲート不合格はサロゲート学習と Pareto から外す (accepted `tooling-nozzle-moo-loop.md` と同じ)。
status: PASS / INFEASIBLE (物理: 設計不成立 `DESIGN`, `L_RAMP_MAX`, `CM_WINDOW`) / FAIL (数値: `DIVERGED`,
`RESIDUAL_RISING`, `NOT_CONVERGED`, `UNSTEADY`, `NO_FORCES`, `ERROR`)。`degraded` = 緩レシピ (再試行) で通った評価
(台帳と pareto.json の両方に残す)。C_M (機体 CG 基準) は台帳に記録し制約は任意 (`opt.cm_min/cm_max` = 加重平均、
`opt.cm_window: {op: [lo, hi]}` = 作動点別の窓 [R6(d)])。
既存 `opt/` (lhs / KrigingSet / propose_infill / hypervolume2d) を再利用。使い方:

  design/.venv-opt/bin/python -m forge_design.opt.driver_sern <problem.yaml> <campaign_dir> [--n-doe 12 --n-iter 3 --batch 2]
  design/.venv-opt/bin/python -m forge_design.opt.driver_sern <problem.yaml> <campaign_dir> --rejudge <out_dir>
      (CFD を回さず既存 run を現行ゲートで再判定: out_dir に ledger_rejudged.jsonl / pareto_rejudged.json / rejudge_summary.md)
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import yaml

from ..evaluate import runner_sern as R
from ..geometry.moc_sern import PlanarMOC, SernKernelSpec
from ..probdef import load_problem
from .doe import lhs
from .ehvi import hypervolume2d, nondominated_mask
from .moo import propose_infill
from .surrogate import KrigingSet

DV_ORDER = ("M_c", "f", "theta_r0_deg", "theta_c0_deg", "L_cowl")
GATE_KEYS = ("C_T", "C_L", "C_M")


# slauWallNormalChi の既定変更 (runner_sern.FLAG_POLICY と同じ値)
from ..evaluate.runner_sern import FLAG_POLICY  # noqa: E402

class DesignInfeasible(ValueError):
    """物理的に成立しない候補 (逆設計不成立 / L_ramp_max 超過)。数値失敗と区別する (R1)。"""


class EvalFailure(RuntimeError):
    """CFD 評価の数値失敗 (ゲート不合格)。fail_class を持つ。"""

    def __init__(self, fail_class: str, msg: str) -> None:
        super().__init__(msg); self.fail_class = fail_class


class _KrgBoth:
    def __init__(self, krg: KrigingSet) -> None:
        self.krg = krg

    def predict(self, X):
        return self.krg.predict(np.atleast_2d(np.asarray(X, dtype=float)))


REQUIRED_SCALAR_GRADIENT = "lsq"


def _learnable(r: dict) -> bool:
    """学習に使ってよい行: PASS・現行の flag_policy・全作動点の実効 scalarGradient が lsq。"""
    if r.get("status") != "PASS" or r.get("flag_policy") != FLAG_POLICY:
        return False
    ops = r.get("ops") or {}
    return bool(ops) and all(o.get("scalar_gradient_effective") == REQUIRED_SCALAR_GRADIENT for o in ops.values())


class SernCampaign:
    def __init__(self, base_yaml, campaign_dir, ref=(-0.90, 20.0), seed: int = 0) -> None:
        self.base_yaml = Path(base_yaml); self.dir = Path(campaign_dir); self.dir.mkdir(parents=True, exist_ok=True)
        self.ref = tuple(ref); self.seed = int(seed)
        self.base_raw = yaml.safe_load(self.base_yaml.read_text())
        self.bounds = np.array([[float(self.base_raw["dv"][k]["min"]), float(self.base_raw["dv"][k]["max"])] for k in DV_ORDER])
        self.ops = self.base_raw["spec"].get("operating_points") or [{"name": "default", "weight": 1.0, "external": self.base_raw["spec"]["external"]}]
        self.optcfg = self.base_raw.get("opt", {})
        self.ledger = self.dir / "ledger.jsonl"
        self.rows = [json.loads(l) for l in self.ledger.read_text().splitlines() if l.strip()] if self.ledger.exists() else []
        _old = sum(1 for r in self.rows if r.get("status") == "PASS" and r.get("flag_policy") != FLAG_POLICY)
        if _old:
            print(f"[campaign] flag_policy が {FLAG_POLICY} でない PASS 行 {_old} 件を学習から除外 "
                  f"(既定変更前の評価: slauWallNormalChi 2026-09-26 / scalarGradient 2026-09-27)", flush=True)
        _nolsq = sum(1 for r in self.rows if r.get("status") == "PASS" and r.get("flag_policy") == FLAG_POLICY and not _learnable(r))
        if _nolsq:
            print(f"[campaign] 実効 scalarGradient が全作動点で lsq と確認できない PASS 行 {_nolsq} 件を学習から除外", flush=True)

    def _write_problem(self, x, path: Path) -> Path:
        raw = json.loads(json.dumps(self.base_raw))
        for k, v in zip(DV_ORDER, np.asarray(x, dtype=float)):
            raw["dv"][k]["value"] = float(v)
        path.write_text(yaml.safe_dump(raw, allow_unicode=True, sort_keys=False)); return path

    # -- 逆設計の成立判定 (粗い kernel で高速化) ----------------------------------
    def feasible(self, x) -> bool:
        try:
            self._design_probe(x); return True
        except Exception:
            return False

    def _design_probe(self, x):
        raw = self.base_raw; geo = raw["geometry"]; st_in = raw["spec"]["inflow"]
        ex0 = raw["spec"]["external"]      # 設計点 (作動点ではない)
        g = float(raw["gas"].get("gamma", 1.4)); cp = float(raw["gas"].get("cp", 1004.5)); Rg = cp * (g - 1) / g
        p_ratio = float(ex0["p_inf"]) / float(st_in["p_in"])
        M_c, f, tr, tc, Lc = [float(v) for v in x]
        spec = SernKernelSpec(M_in=float(st_in["M_in"]), theta_r0=np.deg2rad(tr), theta_c0=np.deg2rad(tc), L_cowl=Lc, gamma=g,
                              p_ext_over_p_in=p_ratio, x_max=float(geo.get("x_max_kernel", 10.0)),
                              nj=int(geo.get("nj_moc_probe", 151)), dx=float(geo.get("dx_moc_probe", 4e-3)))
        k = PlanarMOC(spec).march(stop_at=(f, M_c)); d = k.design_ramp(M_c=M_c, f=f, ds=2e-3)
        if d.info["warnings"]:
            raise ValueError(str(d.info["warnings"]))
        Lmax = float(geo.get("L_ramp_max", 1e9))
        if d.L_ramp > Lmax:
            raise ValueError(f"L_ramp {d.L_ramp:.2f} > {Lmax}")
        return d

    def _eval_op(self, prob, tag: str, op: str, row: dict):
        """1 作動点を標準レシピで回し、ゲート不合格なら緩レシピ (step 倍・cfl 半分) で 1 回だけ再試行する (plan §4.13-1,2)。
        戻り値 (run_dir, rc, metrics, degraded)。**採用条件は rc == 0 かつ gates.verdict == PASS のみ** (旧 §4.13-3 の
        「rc != 0 でも力係数 STEADY なら採用」は撤回、R1)。degraded = 再試行 (緩レシピ) で通った場合 True。"""
        oc = self.optcfg
        # YAML の operating_points[].warm_from で「どの作動点の収束場から立ち上げるか」を指定する
        # (同一メッシュ・熱力学整合リマップ。NPR が遠い作動点には付けない — plan §5.1-1c)
        wf = next((o.get("warm_from") for o in self.ops if o["name"] == op), None)
        warm_src = self.dir / f"{tag}_{wf}" if wf else None
        if warm_src is not None and not (warm_src / "prepare_info.json").exists():
            print(f"   [{tag}/{op}] warm_from={wf} の run が無いので cold start", flush=True)
            warm_src = None
        ladder = [(dict(soft_steps=int(oc.get("soft_steps", 1500)), soft_cfl=float(oc.get("soft_cfl", 0.5)),
                        warm_lam_steps=int(oc.get("warm_lam_steps", 0)), warm_lam_cfl=float(oc.get("warm_lam_cfl", 0.2)),
                        mid_steps=int(oc.get("mid_steps", 0)), warm_src=warm_src,
                        warm_adapt_steps=int(oc.get("warm_adapt_steps", 500))), ""),
                   (dict(soft_steps=2 * int(oc.get("soft_steps", 1500)), soft_cfl=0.5 * float(oc.get("soft_cfl", 0.5)),
                         warm_lam_steps=2 * int(oc.get("warm_lam_steps", 0)), warm_lam_cfl=0.5 * float(oc.get("warm_lam_cfl", 0.2)),
                         mid_steps=2 * int(oc.get("mid_steps", 0)), warm_src=None), "_retry")]   # 再試行は cold start に落とす
        last = None
        for kw, suffix in ladder:
            rd = self.dir / f"{tag}_{op}{suffix}"
            if rd.exists():
                continue
            info = R.prepare(prob, rd, op=op)
            self._check_design(info)
            try:
                rc = R.run_staged(rd, "full", **kw)
            except RuntimeError as e:      # 起動段で落ちた: 壁出力が無いので次の梯子へ
                last = (rd, 1, {"gates": {"verdict": "FAIL", "fail_class": "DIVERGED", "reasons": [str(e)[:200]]}}, True)
                print(f"   [{tag}/{op}{suffix}] {e}", flush=True)
                continue
            out = R.collect(prob, rd, rc=rc, require_residual_pass=bool(oc.get("require_residual_pass", False)))
            g = out["gates"]
            if rc == 0 and g["verdict"] == "PASS":
                return rd, rc, out, suffix == "_retry"
            print(f"   [{tag}/{op}{suffix}] gate {g['fail_class']}: {'; '.join(g['reasons'])[:200]}", flush=True)
            last = (rd, rc, out, True)
        if last is None:                    # 梯子の全段が既存ディレクトリでスキップされた (resume 時)
            raise RuntimeError(f"{tag}/{op}: 梯子の全段が既存ディレクトリ。再開するなら該当 run を消すこと")
        return last

    def _check_design(self, info: dict) -> None:
        """設計段の物理的成立性 (INFEASIBLE)。L_ramp_max は probe (粗い kernel) でなく**最終輪郭**で再検査する (R6(d)/§5.1-8b)。"""
        if info["design"]["warnings"]:
            raise DesignInfeasible("design warnings: " + "; ".join(info["design"]["warnings"]))
        Lmax = float(self.base_raw["geometry"].get("L_ramp_max", 1e9))
        if float(info["design"]["L_ramp"]) > Lmax * (1 + float(self.optcfg.get("l_ramp_tol", 0.0))):
            raise DesignInfeasible(f"L_ramp {info['design']['L_ramp']:.3f} > L_ramp_max {Lmax} (final contour)")

    @staticmethod
    def _op_summary(out: dict, rd) -> dict:
        g = out.get("gates", {})
        return {"C_T": out.get(g.get("objective", "C_T"), out.get("C_T")), "C_T_p": out.get("C_T"), "C_L": out.get("C_L"), "C_M": out.get("C_M"),
                "step": out.get("step"), "run_dir": str(rd), "sep_frac_ramp": out.get("sep_frac_ramp"), "sep_x_min_ramp": out.get("sep_x_min_ramp"),
                "forge_rc": out.get("forge_rc"), "gate": g.get("verdict"), "gate_fail_class": g.get("fail_class"),
                "slau_wall_normal_chi_effective": out.get("slau_wall_normal_chi_effective"), "flag_policy": out.get("flag_policy"),
                "scalar_gradient_effective": out.get("scalar_gradient_effective"),
                "residual": g.get("residual", {}).get("verdict"), "objective": g.get("objective"),
                "steadiness": {k: v.get("verdict") for k, v in g.get("steadiness", {}).get("series", {}).items()}}

    def _cm_check(self, row: dict) -> None:
        """C_M 制約: 加重平均の窓 (`cm_min/cm_max`) と作動点別の窓 (`cm_window: {op: [lo, hi]}`, R6(d))。物理的 INFEASIBLE。"""
        lo, hi = self.optcfg.get("cm_min"), self.optcfg.get("cm_max")
        if (lo is not None and row["C_M_w"] < lo) or (hi is not None and row["C_M_w"] > hi):
            raise DesignInfeasible(f"C_M_w {row['C_M_w']:.3f} outside [{lo}, {hi}]")
        for op, win in (self.optcfg.get("cm_window") or {}).items():
            cm = row["ops"].get(op, {}).get("C_M")
            if cm is None:
                continue
            wlo, whi = (win + [None, None])[:2] if isinstance(win, list) else (win.get("min"), win.get("max"))
            if (wlo is not None and cm < wlo) or (whi is not None and cm > whi):
                raise DesignInfeasible(f"C_M[{op}] {cm:.3f} outside [{wlo}, {whi}]")

    # -- 1 点評価 -----------------------------------------------------------------
    def evaluate(self, x, tag: str) -> dict:
        x = [float(v) for v in np.asarray(x, dtype=float)]
        t0 = time.time(); prob = self._write_problem(x, self.dir / f"{tag}.yaml")
        row = {"tag": tag, "x": x, "status": "FAIL", "fail_class": None, "ops": {}, "note": "", "degraded": False, "degraded_ops": [],
               "flag_policy": FLAG_POLICY}
        try:
            ct_w, L_ramp, cm_w, wsum = 0.0, None, 0.0, 0.0
            for o in self.ops:
                rd, rc, out, degraded = self._eval_op(prob, tag, o["name"], row)
                info = json.loads((rd / "prepare_info.json").read_text())
                L_ramp = info["design"]["L_ramp"]
                g = out.get("gates", {"verdict": "FAIL", "fail_class": "ERROR", "reasons": ["no gates"]})
                row["ops"][o["name"]] = self._op_summary(out, rd)
                if rc != 0 or g["verdict"] != "PASS":
                    raise EvalFailure(g.get("fail_class") or "DIVERGED", f"{o['name']}: rc={rc} {'; '.join(g.get('reasons', []))[:200]}")
                ct = out[g["objective"]]
                if ct is None or not np.isfinite(ct):
                    raise EvalFailure("DIVERGED", f"{o['name']}: objective {g['objective']}={ct}")
                if degraded:
                    row["degraded_ops"].append(o["name"])
                w = float(o.get("weight", 1.0)); wsum += w
                ct_w += w * ct; cm_w += w * out["C_M"]
                # 容量節約: 場は最終ステップのみ残す。**.xmf も一緒に消す** (h5 だけ消すと XDMF が消えた
                # h5 を指し続けて「中身の無い xmf」が残る — 2026-09-05 修正)
                vol = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
                for f_ in vol[:-1]:
                    f_.unlink()
                    f_.with_suffix(".xmf").unlink(missing_ok=True)
            row.update({"status": "PASS", "C_T_w": ct_w / wsum, "C_M_w": cm_w / wsum, "L_ramp": float(L_ramp),
                        "degraded": bool(row["degraded_ops"])})
            self._cm_check(row)
        except DesignInfeasible as e:       # 物理的に不成立 (学習対象外だが「失敗」ではない)
            row["note"] = str(e)[:300]
            row["status"] = "INFEASIBLE"
            row["fail_class"] = "CM_WINDOW" if "C_M" in str(e) else ("L_RAMP_MAX" if "L_ramp_max" in str(e) else "DESIGN")
        except EvalFailure as e:            # 数値失敗 (ゲート不合格)
            row["note"] = str(e)[:300]; row["status"] = "FAIL"; row["fail_class"] = e.fail_class
        except Exception as e:              # 想定外 (I/O 等)。INFEASIBLE と混同しない
            row["note"] = f"{type(e).__name__}: {str(e)[:280]}"; row["status"] = "FAIL"; row["fail_class"] = "ERROR"
        row["elapsed_s"] = time.time() - t0
        with open(self.ledger, "a") as fh:
            fh.write(json.dumps(row, ensure_ascii=False) + "\n")
        self.rows.append(row)
        print(f"[{tag}] {row['status']}{'/' + row['fail_class'] if row['fail_class'] else ''}{' (degraded)' if row['degraded'] else ''} "
              f"x={np.round(x, 3).tolist()} C_T_w={row.get('C_T_w')} L={row.get('L_ramp')} C_M_w={row.get('C_M_w')} ({row['elapsed_s']:.0f}s) {row['note']}", flush=True)
        return row

    def _XF(self):
        # slauWallNormalChi の既定変更 (2026-09-26) 前後の評価を同じ応答関数として学習しない (plan
        # convection-slau-wall-normal-chi-default §4.4、codex plan M3)。flag_policy の無い旧行と不一致の行は除外する。
        # さらに mesh.scalarGradient の node 既定 lsq 化 (2026-09-27) 以降は、**全作動点の実効値が lsq と確認できた行だけ**を使う
        # (日付の一致だけでは gg 評価・不明が混ざる。codex diagnose 2026-09-27、plan gradient-scalar-lsq-unification #6)。
        ok = [r for r in self.rows if _learnable(r)]
        X = np.array([r["x"] for r in ok]); F = np.array([[-r["C_T_w"], r["L_ramp"]] for r in ok])
        return X, F

    def run(self, n_doe: int, n_iter: int, batch: int) -> None:
        done = {r["tag"] for r in self.rows}
        cand = lhs(self.bounds, 3 * n_doe, seed=self.seed)
        doe = []
        for x in cand:
            if len(doe) >= n_doe:
                break
            if self.feasible(x):
                doe.append(x)
        print(f"DOE {len(doe)} 点", flush=True)
        for i, x in enumerate(doe):
            tag = f"doe_{i:03d}"
            if tag not in done:
                self.evaluate(x, tag)
        for it in range(n_iter):
            X, F = self._XF()
            if len(X) < 4:
                raise RuntimeError("PASS 評価が少なすぎる (<4)")
            print(f"--- iter {it}: PASS {len(X)} 点, HV={hypervolume2d(F, self.ref):.4f} ---", flush=True)
            krg = KrigingSet(self.bounds).fit(X, F)
            xs = propose_infill(_KrgBoth(krg), X, F, self.ref, self.bounds, n_infill=batch, seed=self.seed + 1000 + it, feasible=self.feasible)
            for j, x in enumerate(xs):
                tag = f"inf_{it:02d}_{j}"
                if tag not in done:
                    self.evaluate(x, tag)
        summary = self.summary()
        (self.dir / "pareto.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
        print(json.dumps(summary, indent=1, ensure_ascii=False), flush=True)

    def summary(self, rows=None) -> dict:
        """Pareto 要約。**degraded / tag / 作動点ごとのゲート要約を落とさない** (R1: pareto.json でも追える)。"""
        rows = self.rows if rows is None else rows
        ok = [r for r in rows if r["status"] == "PASS"]
        X = np.array([r["x"] for r in ok]); F = np.array([[-r["C_T_w"], r["L_ramp"]] for r in ok])
        pareto = []
        if len(ok):
            for i in np.where(nondominated_mask(F))[0]:
                r = ok[i]
                pareto.append({"tag": r["tag"], "x": dict(zip(DV_ORDER, X[i].tolist())), "C_T_w": float(-F[i, 0]), "L_ramp": float(F[i, 1]),
                               "C_M_w": r["C_M_w"], "degraded": bool(r.get("degraded")), "degraded_ops": r.get("degraded_ops", []),
                               "ops": {op: {k: v.get(k) for k in ("C_T", "C_M", "gate", "residual", "steadiness")} for op, v in r.get("ops", {}).items()}})
        pareto.sort(key=lambda r: r["L_ramp"])
        classes = {}
        for r in rows:
            k = r["status"] + ("/" + r["fail_class"] if r.get("fail_class") else "")
            classes[k] = classes.get(k, 0) + 1
        return {"n_eval": len(rows), "n_pass": int(len(ok)), "n_degraded": int(sum(1 for r in ok if r.get("degraded"))),
                "hv": (hypervolume2d(F, self.ref) if len(ok) else 0.0), "ref": self.ref, "status_counts": classes,
                "gate_policy": "R1: rc==0 + finite field + residual no NaN/rising + objective & C_T/C_L/C_M STEADY (no divergent adoption)",
                "operating_points": self.ops, "pareto": pareto}

    # -- 既存キャンペーンの再判定 (CFD なし) ----------------------------------------
    def rejudge(self, out_dir) -> dict:
        """台帳の各評価を現行ゲートで判定し直す。元 run は読むだけ (metrics/force_history は out_dir/<tag>_<op>/ に書く)。
        rc は台帳に無いので run_case_stdout.log から復元する。"""
        out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
        new_rows = []; lines = ["| tag | 旧 status | 新 status | fail_class | degraded | C_T_w (旧→新) | 理由 |", "| --- | --- | --- | --- | --- | --- | --- |"]
        for r0 in self.rows:
            tag = r0["tag"]; prob = self.dir / f"{tag}.yaml"
            row = {"tag": tag, "x": r0["x"], "status": "FAIL", "fail_class": None, "ops": {}, "note": "", "degraded": False, "degraded_ops": [],
                   "old_status": r0["status"], "old_fail_class": r0.get("fail_class"), "old_C_T_w": r0.get("C_T_w")}
            try:
                if not prob.exists():
                    raise EvalFailure("ERROR", "problem yaml missing")
                ct_w, cm_w, wsum, L_ramp = 0.0, 0.0, 0.0, None
                for o in self.ops:
                    rd_std, rd_retry = self.dir / f"{tag}_{o['name']}", self.dir / f"{tag}_{o['name']}_retry"
                    cands = [d for d in (rd_std, rd_retry) if (d / "prepare_info.json").exists()]
                    if not cands:
                        raise EvalFailure("NO_FORCES", f"{o['name']}: no run dir")
                    chosen = None
                    for rd in cands:       # 標準 → 再試行の順に、現行ゲートを通る最初の run を採る
                        out = R.collect(prob, rd, out_dir=out_dir / rd.name, require_residual_pass=bool(self.optcfg.get("require_residual_pass", False)))
                        row["ops"][o["name"]] = self._op_summary(out, rd)
                        if out["forge_rc"] == 0 and out["gates"]["verdict"] == "PASS":
                            chosen = (rd, out); break
                    if chosen is None:
                        g = out["gates"]
                        raise EvalFailure(g.get("fail_class") or "DIVERGED", f"{o['name']}: rc={out['forge_rc']} {'; '.join(g['reasons'])[:200]}")
                    rd, out = chosen
                    info = json.loads((rd / "prepare_info.json").read_text()); self._check_design(info); L_ramp = info["design"]["L_ramp"]
                    if rd.name.endswith("_retry"):
                        row["degraded_ops"].append(o["name"])
                    w = float(o.get("weight", 1.0)); wsum += w; ct_w += w * out[out["gates"]["objective"]]; cm_w += w * out["C_M"]
                row.update({"status": "PASS", "C_T_w": ct_w / wsum, "C_M_w": cm_w / wsum, "L_ramp": float(L_ramp), "degraded": bool(row["degraded_ops"])})
                self._cm_check(row)
            except DesignInfeasible as e:
                row["note"] = str(e)[:300]; row["status"] = "INFEASIBLE"
                row["fail_class"] = "CM_WINDOW" if "C_M" in str(e) else ("L_RAMP_MAX" if "L_ramp_max" in str(e) else "DESIGN")
            except EvalFailure as e:
                row["note"] = str(e)[:300]; row["status"] = "FAIL"; row["fail_class"] = e.fail_class
            except Exception as e:
                row["note"] = f"{type(e).__name__}: {str(e)[:280]}"; row["status"] = "FAIL"; row["fail_class"] = "ERROR"
            new_rows.append(row)
            ctn = row.get("C_T_w"); cto = r0.get("C_T_w")
            lines.append(f"| {tag} | {r0['status']}{'/' + r0['fail_class'] if r0.get('fail_class') else ''} | {row['status']} | {row['fail_class'] or ''} | "
                         f"{'yes' if row['degraded'] else ''} | {'' if cto is None else f'{cto:.4f}'} → {'' if ctn is None else f'{ctn:.4f}'} | {row['note'][:120]} |")
        with open(out_dir / "ledger_rejudged.jsonl", "w") as fh:
            for r in new_rows:
                fh.write(json.dumps(r, ensure_ascii=False) + "\n")
        summ = self.summary(new_rows); summ["source_campaign"] = str(self.dir)
        old_ok = [r for r in self.rows if r["status"] == "PASS"]
        summ["old"] = {"n_pass": len(old_ok), "hv": (hypervolume2d(np.array([[-r["C_T_w"], r["L_ramp"]] for r in old_ok]), self.ref) if old_ok else 0.0)}
        (out_dir / "pareto_rejudged.json").write_text(json.dumps(summ, indent=1, ensure_ascii=False))
        md = [f"# 再判定 (R1 ゲート): {self.dir}", "", f"旧: PASS {summ['old']['n_pass']} / HV {summ['old']['hv']:.4f} → "
              f"新: PASS {summ['n_pass']} (degraded {summ['n_degraded']}) / HV {summ['hv']:.4f}", "", f"status 内訳: {summ['status_counts']}", ""] + lines
        (out_dir / "rejudge_summary.md").write_text("\n".join(md) + "\n")
        print("\n".join(md), flush=True)
        return summ


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="⑤ SERN 多作動点 MOO キャンペーン")
    ap.add_argument("problem"); ap.add_argument("campaign_dir")
    ap.add_argument("--n-doe", type=int, default=12); ap.add_argument("--n-iter", type=int, default=3)
    ap.add_argument("--batch", type=int, default=2); ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--ref", type=float, nargs=2, default=(-0.90, 20.0))
    ap.add_argument("--rejudge", default=None, metavar="OUT_DIR", help="CFD を回さず既存キャンペーンを現行ゲートで再判定して OUT_DIR に書く")
    a = ap.parse_args(argv)
    c = SernCampaign(a.problem, a.campaign_dir, ref=tuple(a.ref), seed=a.seed)
    if a.rejudge:
        c.rejudge(a.rejudge); return 0
    c.run(a.n_doe, a.n_iter, a.batch)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
