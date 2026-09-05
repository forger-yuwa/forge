"""⑤ SERN の多作動点 MOO ドライバ (S6)。

dv = (M_c, f, theta_r0_deg, theta_c0_deg, L_cowl) [YAML の dv.min/max]。1 評価 = 作動点数分の forge run
(`spec.operating_points[]`、重み w_k)。目的 (最小化):
    f1 = −Σ_k w_k C_T^(k)      (Kriging。RANS では摩擦込み C_T)
    f2 = L_ramp / H            (逆設計の決定的結果。候補点では Kriging で近似)
ゲート: 逆設計成立 (kernel 内に key point、C⁻ が直線ランプに着地)、メッシュ品質 PASS、各作動点で
NaN なし・力係数 STEADY。C_M (機体 CG 基準) は台帳に記録し制約は任意 (`opt.cm_min/cm_max`)。
既存 `opt/` (lhs / KrigingSet / propose_infill / hypervolume2d) を再利用。使い方:

  design/.venv-opt/bin/python -m forge_design.opt.driver_sern <problem.yaml> <campaign_dir> [--n-doe 12 --n-iter 3 --batch 2]
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


class _KrgBoth:
    def __init__(self, krg: KrigingSet) -> None:
        self.krg = krg

    def predict(self, X):
        return self.krg.predict(np.atleast_2d(np.asarray(X, dtype=float)))


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
        """1 作動点を標準レシピで回し、落ちたら緩レシピ (step 倍・cfl 半分) で 1 回だけ再試行する (plan §4.13)。
        戻り値 (run_dir, rc, metrics, degraded)。degraded = rc != 0 でも力係数が使えた場合に True。"""
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
            if info["design"]["warnings"]:
                raise ValueError("design warnings: " + "; ".join(info["design"]["warnings"]))
            try:
                rc = R.run_staged(rd, "full", **kw)
            except RuntimeError as e:      # 起動段で落ちた: 壁出力が無いので次の梯子へ
                last = (rd, 1, {}, True); print(f"   [{tag}/{op}{suffix}] {e}", flush=True)
                continue
            out = R.collect(prob, rd)
            st = out.get("steadiness", {})
            ok = all(st.get(k, {}).get("verdict") == "STEADY" for k in ("C_T", "C_L", "C_M"))
            if rc == 0 and ok:
                return rd, rc, out, False
            last = (rd, rc, out, True)
            if ok:                          # 発散はしたが力積分は収束済み → 再試行せず採用 (§4.13-3)
                return last
        if last is None:                    # 梯子の全段が既存ディレクトリでスキップされた (resume 時)
            raise RuntimeError(f"{tag}/{op}: 梯子の全段が既存ディレクトリ。再開するなら該当 run を消すこと")
        return last

    # -- 1 点評価 -----------------------------------------------------------------
    def evaluate(self, x, tag: str) -> dict:
        x = [float(v) for v in np.asarray(x, dtype=float)]
        t0 = time.time(); prob = self._write_problem(x, self.dir / f"{tag}.yaml")
        row = {"tag": tag, "x": x, "status": "FAIL", "fail_class": None, "ops": {}, "note": ""}
        try:
            ct_w, L_ramp, cm_w, wsum = 0.0, None, 0.0, 0.0
            for o in self.ops:
                rd, rc, out, degraded = self._eval_op(prob, tag, o["name"], row)
                info = json.loads((rd / "prepare_info.json").read_text())
                L_ramp = info["design"]["L_ramp"]
                st = out.get("steadiness", {}); ct = out.get("C_T_with_shear", out.get("C_T"))   # RANS は摩擦込み
                if ct is None or not np.isfinite(ct):
                    row["fail_class"] = "DIVERGED"; raise RuntimeError(f"{o['name']}: forge rc={rc} / C_T={ct}")
                # plan §4.13-3: 力係数が 4 点以上で 3 成分とも STEADY なら、rc != 0 でも degraded として採用する
                if any(st.get(k, {}).get("verdict") != "STEADY" for k in ("C_T", "C_L", "C_M")):
                    row["fail_class"] = "UNSTEADY"; raise RuntimeError(f"{o['name']}: C_T {st.get('C_T')}")
                if degraded:
                    row.setdefault("degraded_ops", []).append(o["name"])
                w = float(o.get("weight", 1.0)); wsum += w
                ct_w += w * ct; cm_w += w * out["C_M"]
                row["ops"][o["name"]] = {"C_T": ct, "C_T_p": out["C_T"], "C_L": out["C_L"], "C_M": out["C_M"], "step": out["step"], "run_dir": str(rd),
                                         "sep_frac_ramp": out.get("sep_frac_ramp"), "sep_x_min_ramp": out.get("sep_x_min_ramp")}
                # 容量節約: 場は最終ステップのみ残す。**.xmf も一緒に消す** (h5 だけ消すと XDMF が消えた
                # h5 を指し続けて「中身の無い xmf」が残る — 2026-09-05 修正)
                vol = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
                for f_ in vol[:-1]:
                    f_.unlink()
                    f_.with_suffix(".xmf").unlink(missing_ok=True)
            row.update({"status": "PASS", "C_T_w": ct_w / wsum, "C_M_w": cm_w / wsum, "L_ramp": float(L_ramp),
                        "degraded": bool(row.get("degraded_ops"))})
            lo, hi = self.optcfg.get("cm_min"), self.optcfg.get("cm_max")
            if (lo is not None and row["C_M_w"] < lo) or (hi is not None and row["C_M_w"] > hi):
                row["status"] = "INFEASIBLE"; row["fail_class"] = "CM_WINDOW"
        except Exception as e:
            row["note"] = str(e)[:300]
            if row["fail_class"] is None:
                row["fail_class"] = "INFEASIBLE" if isinstance(e, ValueError) else "RETRYABLE"
        row["elapsed_s"] = time.time() - t0
        with open(self.ledger, "a") as fh:
            fh.write(json.dumps(row, ensure_ascii=False) + "\n")
        self.rows.append(row)
        print(f"[{tag}] {row['status']} x={np.round(x, 3).tolist()} C_T_w={row.get('C_T_w')} L={row.get('L_ramp')} C_M_w={row.get('C_M_w')} ({row['elapsed_s']:.0f}s) {row['note']}", flush=True)
        return row

    def _XF(self):
        ok = [r for r in self.rows if r["status"] == "PASS"]
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
        X, F = self._XF(); mask = nondominated_mask(F)
        pareto = [{"x": dict(zip(DV_ORDER, X[i].tolist())), "C_T_w": float(-F[i, 0]), "L_ramp": float(F[i, 1]),
                   "C_M_w": [r for r in self.rows if r["status"] == "PASS"][i]["C_M_w"]} for i in np.where(mask)[0]]
        pareto.sort(key=lambda r: r["L_ramp"])
        summary = {"n_eval": len(self.rows), "n_pass": int(len(X)), "hv": hypervolume2d(F, self.ref), "ref": self.ref,
                   "operating_points": self.ops, "pareto": pareto}
        (self.dir / "pareto.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
        print(json.dumps(summary, indent=1, ensure_ascii=False), flush=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="⑤ SERN 多作動点 MOO キャンペーン")
    ap.add_argument("problem"); ap.add_argument("campaign_dir")
    ap.add_argument("--n-doe", type=int, default=12); ap.add_argument("--n-iter", type=int, default=3)
    ap.add_argument("--batch", type=int, default=2); ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--ref", type=float, nargs=2, default=(-0.90, 20.0))
    a = ap.parse_args(argv)
    SernCampaign(a.problem, a.campaign_dir, ref=tuple(a.ref), seed=a.seed).run(a.n_doe, a.n_iter, a.batch)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
