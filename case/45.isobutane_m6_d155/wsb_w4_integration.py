"""plan tooling-nozzle-wall-single-bspline §6.0 W4 の統合検査 (2026-10-07 改訂: 対象は run_0166_ns_n012_N1)。CFD 0 step。

1 本の B-spline はノット挿入の版で `pw_upstream: poly` の壁だけを受け付けるので、ランプの壁の run_0147 ではなく、
上流の多項式化の plan の U4 の N1 (run_0166) の結果で検査する。AWS 上で run_0166 と同じコード・環境・変換器で回す。

  prepare : N1 の問題で prepare_ns を 2 腕 (A キー無し / B physical_wall_repr: single_bspline) → _band_ab/wsb_w4i/{A_key_absent,B_single_bspline}
            (引数は ns_n012.py prep-dry と同じ。IC は貼らない)
  judge   : (1) A/B のソルバ入力のビット同一 (W3 の判定) / (2) 腕 A が run_0166 の壁を再現 /
            (3) run_0166 の結果に腕 B の壁ファイルを添えた写しで nozzle_report を回し、保存した係数から評価していること・
                壁形状以外の評価量が run_0166 自身の報告とビット同一 → _band_ab/wsb/W4_integration.json

usage (AWS の case dir で):
  FORGE_CONVERTER=<conv_tolerant.sh> REAL_CONVERTER=<変換器> FORGE_ALLOW_UNVERIFIED_SPECIES=1 python3 wsb_w4_integration.py prepare
  python3 wsb_w4_integration.py judge
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C))
sys.path.insert(0, str(C.parents[1] / "design"))
import wsb_prepare as WP  # noqa: E402
import wsb_verify as WV  # noqa: E402

OUT = C / "_band_ab" / "wsb_w4i"
RUN = C / "run_0166_ns_n012_N1"
EULER = C / "run_0164_euler_e4_recal_d1"
PN = "problem_d155_ns_n012_N1.yaml"
ARMS = {"A_key_absent": None, "B_single_bspline": "single_bspline"}
IL = "initial_line_run: run_0062_euler_wallfit_fit_r1_ext6k"
UP = "  pw_upstream: poly"


def problem_copy(arm: str, key) -> Path:
    txt = (C / PN).read_text()
    if txt.count(IL) != 1 or txt.count(UP) != 1:
        raise SystemExit(f"{PN} に {IL!r}・{UP!r} が 1 回ずつ無い")
    txt = txt.replace(IL, f"initial_line_run: {C / 'run_0062_euler_wallfit_fit_r1_ext6k'}")
    if key:
        i = txt.index(UP)
        j = txt.index("\n", i) + 1
        txt = txt[:j] + f"  physical_wall_repr: {key}\n" + txt[j:]
    inp = OUT / "inputs"
    inp.mkdir(parents=True, exist_ok=True)
    dst = inp / f"{arm}.problem.yaml"
    if dst.exists():
        raise SystemExit(f"{dst} が既にある (上書きしない)")
    dst.write_text(txt)
    return dst


def prepare():
    import ns_n012 as NS
    conv = os.environ.get("FORGE_CONVERTER")
    if not conv or not Path(conv).is_file():
        raise SystemExit("FORGE_CONVERTER (変換器) を指定する")
    kf = NS.k_f_checked(C / PN)
    for arm, key in ARMS.items():
        d = OUT / arm
        if d.exists():
            raise SystemExit(f"{d} が既にある (上書きしない)")
        prob = problem_copy(arm, key)
        info = NS._prepare_ns(prob, d, nsteps=12000, ic_from=None, initializer=NS.initializer(kf), cfl_main=5.0, implicit_relax=NS.RELAX)
        rec = {"arm": arm, "physical_wall_repr": key, "problem_copy": str(prob), "problem_copy_sha256": WP.sha_file(prob),
               "problem_source": str(C / PN), "problem_source_sha256": WP.sha_file(C / PN), "k_f": kf,
               "converter": conv, "real_converter": os.environ.get("REAL_CONVERTER"),
               "git_head": subprocess.run(["git", "-C", str(C), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
               "git_dirty_design": bool(subprocess.run(["git", "-C", str(C.parents[1]), "status", "--porcelain", "--", "design"],
                                                       capture_output=True, text=True).stdout.strip()),
               "prepare_info_keys": sorted(info)}
        (OUT / f"{arm}.record.json").write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=str))
        print(f"[w4i] {arm}: {d}", flush=True)


def _mesh_identical(pa: Path, pb: Path) -> dict:
    """nozzle.h5 の /MESH・/CELLS・/PLANES・/BCONDS の幾何 (IC を含まない部分) をデータセット単位で比べる。"""
    cmp = WV.h5_compare(pa, pb)
    geo = {k: v["identical"] for k, v in cmp["datasets"].items()
           if k.startswith(("MESH/", "CELLS/", "PLANES/")) or (k.startswith("BCONDS/") and "/VALUE/" not in k)}
    return {"n": len(geo), "all_identical": bool(geo) and all(geo.values()), "diff": [k for k, v in geo.items() if not v]}


def report_copy() -> Path:
    """run_0166 の結果 (report/ 以外) をシンボリックリンクで並べ、腕 B の壁ファイルを添えた写し。"""
    from forge_design.geometry.wall_axismach import WALL_FILE
    dst = OUT / "run_0166_sb"
    if dst.exists():
        raise SystemExit(f"{dst} が既にある (上書きしない)")
    dst.mkdir()
    for p in sorted(RUN.iterdir()):
        if p.name in ("report", WALL_FILE) or p.name.startswith("report_"):
            continue
        (dst / p.name).symlink_to(p.resolve())
    (dst / WALL_FILE).write_bytes((OUT / "B_single_bspline" / WALL_FILE).read_bytes())
    return dst


def judge():
    from forge_design.geometry.wall_axismach import WALL_FILE, load_wall_file
    A, B = OUT / "A_key_absent", OUT / "B_single_bspline"
    ck = {}
    # (1) A/B のソルバ入力 (W3 の判定)
    cmp = WV.compare_dirs(A, B)
    j = WV.judge_w3(cmp)
    ck["solver_input_A_vs_B"] = {"pass": j["pass"], "items": j["items"], "nozzle_msh_text": cmp["nozzle_msh"],
                                 "prepare_info_diff": cmp["prepare_info_diff"], "files": cmp["files"]}
    # (2) 腕 A が run_0166 の壁を再現
    same = {n: WV.sha_file(A / n) == WV.sha_file(RUN / n) for n in ("wall_physical.csv", "wall_design.csv", "delta_r_initial.csv")}
    ia, ir = (json.loads((d / "prepare_info.json").read_text()) for d in (A, RUN))
    mesh = _mesh_identical(A / "nozzle.h5", RUN / "nozzle.h5")
    ck["arm_A_reproduces_run_0166"] = {"pass": bool(all(same.values()) and mesh["all_identical"] and ia["throat_physical"] == ir["throat_physical"]),
                                       "files_identical": same, "nozzle_h5_geometry": mesh,
                                       "throat_physical": {"A": ia["throat_physical"], "run_0166": ir["throat_physical"]}}
    # (3) 報告: run_0166 の結果 + 腕 B の壁ファイル
    R = OUT / "run_0166_sb"
    if not R.exists():
        report_copy()
    rep_dir = R / "report"
    if not (rep_dir / "report.json").is_file():
        r = subprocess.run([sys.executable, "-m", "forge_design.report.nozzle_report", str(R), "--euler", str(EULER),
                            "--no-pptx", "--wall-over-frac", "5"], cwd=str(C.parents[1] / "design"), capture_output=True, text=True)
        (OUT / "report_stdout.log").write_text(r.stdout + "\n--- stderr ---\n" + r.stderr)
        if r.returncode != 0:
            raise SystemExit(f"nozzle_report が失敗 (rc {r.returncode}) — {OUT / 'report_stdout.log'}")
    new, ref = (json.loads((d / "report" / "report.json").read_text()) for d in (R, RUN))
    W = load_wall_file(R / WALL_FILE)
    S = float(ir["scale_m"])
    ws = new["metrics"]["wall_shape"]
    exit_eval = float(W["physical"].r(np.r_[W["domain"][1]])[0] * S)
    # 場所のラベル (run のパス・親ディレクトリ名) は写しで必ず変わるので比べない: パスの文字列は元の run に置き換え、
    # conditions の run・run_name・case は除く (2026-10-07 修正。初版は case と wall_resolution.cmd のパスまで比べて FAIL、
    # 初版の出力は _band_ab/wsb/W4_integration_v1_pathlabels.json に残した)
    def _norm(o):
        if isinstance(o, str):
            return o.replace(str(R), str(RUN))
        if isinstance(o, dict):
            return {k: _norm(v) for k, v in o.items()}
        if isinstance(o, list):
            return [_norm(v) for v in o]
        return o
    LOC = ("run", "run_name", "case")
    m_new = {k: _norm(v) for k, v in new["metrics"].items() if k != "wall_shape"}
    m_ref = {k: v for k, v in ref["metrics"].items() if k != "wall_shape"}
    c_new = {k: _norm(v) for k, v in new["conditions"].items() if k not in LOC}
    c_ref = {k: v for k, v in ref["conditions"].items() if k not in LOC}
    ck["report"] = {"pass": bool(ws.get("wall_source") == "saved_coefficients" and ws.get("wall_repr") == "single_bspline"
                                 and ws.get("exit_radius_m") == exit_eval and m_new == m_ref and c_new == c_ref),
                    "wall_source": ws.get("wall_source"), "wall_repr": ws.get("wall_repr"), "wall_file": ws.get("wall_file"),
                    "exit_radius_m": {"report": ws.get("exit_radius_m"), "saved_coefficients": exit_eval},
                    "metrics_except_wall_shape_identical": m_new == m_ref, "metrics_diff": WV.dict_diff(m_ref, m_new)[:30],
                    "conditions_except_run_identical": c_new == c_ref, "conditions_diff": WV.dict_diff(c_ref, c_new)[:30],
                    "wall_shape_ref_legacy_csv": {k: ref["metrics"]["wall_shape"].get(k) for k in ("wall_source", "exit_radius_m", "r2_highfreq_max_x_gt2")},
                    "wall_shape_new": {k: ws.get(k) for k in ("wall_source", "wall_repr", "exit_radius_m", "r2_highfreq_max_x_gt2")},
                    "figure": str(rep_dir / "fig_wall_shape.png")}
    ok = all(v["pass"] for v in ck.values())
    out = {"item": "W4 統合検査 (plan §6.0 W4、2026-10-07 改訂: run_0166_ns_n012_N1)", "run": str(RUN), "euler": str(EULER),
           "verdict": "PASS" if ok else "FAIL (差の経路を特定して止める)", "checks": ck,
           "records": {a: json.loads((OUT / f"{a}.record.json").read_text()) for a in ARMS}}
    dst = C / "_band_ab" / "wsb" / "W4_integration.json"
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str))
    print(json.dumps({"verdict": out["verdict"], **{k: v["pass"] for k, v in ck.items()}}, ensure_ascii=False))


if __name__ == "__main__":
    {"prepare": prepare, "judge": judge}[sys.argv[1]]()
