"""case/70: 要求カードへの M_knot の組み込み (plan tooling-nozzle-mknot-auto §5.1 #3) の §6 V3・V4 (設計計算だけ、forge 0 step)。

  v3    後方互換: 数値の overrides.geometry.M_knot のカード 2 枚 (M5 デモの例 + 2.25・c2b + 2.069270614) を、#3 の前の生成器
        (--old-rev の wt_card.py・mknot_select.py を一時ディレクトリの写しのパッケージに入れたもの) と作業ツリーの生成器で作り、
        生成物 (問題 YAML・card_log・Euler の写し) を比べる。正規化するのは生成日・所要時間 (time_s)・out_sha256・雛形の絶対パス・
        出力先だけ。作業ツリーの生成は探索の関数 (search_range・select_mknot・resolution_gate・evaluate_mknot・make_eval_fn) を
        呼んだら失敗する差し替えの中で回す。合格: 全部同じ・探索の関数を 1 度も呼ばない。
  v4    統合: card_m5_air_demo_c2b.yaml と同じ条件で M_knot を書かないカード (_mknot_card_v4/card_c2b_auto_{final,explore}.yaml) を
        生成し (CLI と同じ build_problem)、`check` と同じ判定をする。所要時間を記録する。
  check 生成した問題 YAML の判定 (§6 V4): 寸法ゲート PASS (長さ 1 mm・径 0.02 mm)、全部の制約 (記録の N・2N と、--reeval なら出力 YAML を
        n_axis_inv と 2 倍で独立に評価し直した値)、全部の境目からの距離 ≥ 0.05、解像度のゲート PASS、記録の選んだ値 = YAML の
        geometry.M_knot = card_log、r_t の一致。

usage (forge-deliv の直下、design/.venv-opt の python):
  python case/70.m5_air_wt_demo/mknot_card_v34.py v3 [--old-rev d9f554f4] [--out DIR]
  python case/70.m5_air_wt_demo/mknot_card_v34.py v4 final|explore
  python case/70.m5_air_wt_demo/mknot_card_v34.py check case/70.m5_air_wt_demo/_mknot_card_v4/problem_c2b_auto_final.yaml [--reeval]
所要 (2026-10-14 手元の実測・目安): v3 約 3 分 (生成 4 回)、v4 final 約 45 分 (並列 2/1 のとき; 1/1 ではもっと長い)、check --reeval 約 1 分
(n_axis_inv 4800 の評価 1 回、メモリ 約 3.9 GB)。M_knot の評価の並列は wt_card.MKNOT_WORKERS (N・2N とも 1)。
"""
import json
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import yaml

C = Path(__file__).resolve().parent
REPO = C.parents[1]
DESIGN = REPO / "design"
V4DIR = C / "_mknot_card_v4"
SEARCH_FUNCS = ("search_range", "select_mknot", "resolution_gate", "evaluate_mknot", "make_eval_fn")
TIMEKEYS = {"time_s", "elapsed_s", "generated", "out_sha256"}


def _arg(argv, name, default):
    if name in argv:
        k = argv.index(name)
        v = argv[k + 1]
        del argv[k:k + 2]
        return v
    return default


# --- V3 ---------------------------------------------------------------------------------------------------------------
BUILD = r'''
import sys, time
root, card, out, forbid = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] == "1"
sys.path.insert(0, root)
from forge_design.intake import wt_card as W
calls = []
if forbid:
    from forge_design.evaluate import mknot_select as MS
    def f_(name):
        def f(*a, **k):
            calls.append(name)
            raise AssertionError(f"探索の関数 {name} が呼ばれた")
        return f
    for n in %r:
        setattr(MS, n, f_(n))
t0 = time.time()
logd = W.build_problem(card, out, log=lambda *a, **k: None)
print("BUILD", W.__file__, round(time.time() - t0, 1), calls, logd["dimensions"]["pass"], logd["design_chain"]["M_knot"], flush=True)
''' % (SEARCH_FUNCS,)


def _norm_str(s, roots):
    for r in roots:
        s = s.replace(r, "<ROOT>")
    s = re.sub(r"/tmp/[^/\s'\"]+/", "<TMP>/", s)
    s = re.sub(r"(^|/)(old_out|new_out)/", r"\1<OUT>/", s)
    s = re.sub(r"\(\d{4}-\d{2}-\d{2}\)", "(DATE)", s)
    s = re.sub(r"generated: '\d{4}-\d{2}-\d{2}'", "generated: DATE", s)
    return re.sub(r"time_s: [0-9.e+-]+", "time_s: T", s)


def _norm(v, roots, key=None):
    if key in TIMEKEYS:
        return "<T>"
    if isinstance(v, dict):
        return {k: _norm(x, roots, k) for k, x in v.items()}
    if isinstance(v, list):
        return [_norm(x, roots) for x in v]
    return _norm_str(v, roots) if isinstance(v, str) else v


def _diff(a, b, path="", out=None):
    out = [] if out is None else out
    if type(a) is not type(b):
        out.append((path, a, b))
    elif isinstance(a, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append((f"{path}.{k}", a.get(k, "<無し>"), b.get(k, "<無し>")))
            else:
                _diff(a[k], b[k], f"{path}.{k}", out)
    elif isinstance(a, list):
        if len(a) != len(b):
            out.append((f"{path}[len]", len(a), len(b)))
        for i, (x, y) in enumerate(zip(a, b)):
            _diff(x, y, f"{path}[{i}]", out)
    elif a != b:
        out.append((path, a, b))
    return out


def run_v3(old_rev: str, outdir: Path) -> bool:
    py = sys.executable
    outdir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="mknot_v3_", dir=outdir))
    # #3 の前の生成器: 作業ツリーの forge_design の写しに old_rev の wt_card.py・mknot_select.py を入れる (他のファイルは同じ)。
    # パッケージは REPO 直下の solver_density_cuda を parents[3] で読むので、写しの横に symlink を置く
    old_root = work / "old"
    shutil.copytree(DESIGN / "forge_design", old_root / "forge_design", ignore=shutil.ignore_patterns("__pycache__"))
    for rel in ("design/forge_design/intake/wt_card.py", "design/forge_design/evaluate/mknot_select.py"):
        txt = subprocess.run(["git", "-C", str(REPO), "show", f"{old_rev}:{rel}"], check=True, capture_output=True, text=True).stdout
        (old_root / rel.removeprefix("design/")).write_text(txt)
    (work / "solver_density_cuda").symlink_to(REPO / "solver_density_cuda")
    m5 = yaml.safe_load((DESIGN / "forge_design" / "intake" / "examples" / "m5_air_demo.yaml").read_text())
    m5["overrides"] = {"geometry.M_knot": 2.25}
    c2b = yaml.safe_load((C / "card_m5_air_demo_c2b.yaml").read_text())
    c2b["overrides"] = dict(c2b["overrides"], **{"geometry.M_knot": 2.069270614})
    (work / "build.py").write_text(BUILD)
    ok_all = True
    for name, card in (("m5_num", m5), ("c2b_num", c2b)):
        (work / f"card_{name}.yaml").write_text(yaml.safe_dump(card, sort_keys=False, allow_unicode=True))
        for tag, root, forbid in (("old_out", old_root, "0"), ("new_out", DESIGN, "1")):
            (work / tag).mkdir(exist_ok=True)
            r = subprocess.run([py, "build.py", str(root), f"card_{name}.yaml", f"{tag}/problem_{name}.yaml", forbid], cwd=work,
                               capture_output=True, text=True)
            line = [ln for ln in r.stdout.splitlines() if ln.startswith("BUILD")]
            print(f"  {tag} {name}: rc {r.returncode} {line[-1] if line else r.stderr[-400:]}", flush=True)
            if r.returncode != 0:
                ok_all = False
        roots = [str(old_root), str(DESIGN)]
        o, n = work / "old_out" / f"problem_{name}.yaml", work / "new_out" / f"problem_{name}.yaml"
        if not (o.exists() and n.exists()):
            continue
        yo = [_norm_str(x, roots) for x in o.read_text().splitlines()]
        yn = [_norm_str(x, roots) for x in n.read_text().splitlines()]
        yd = [(i, a, b) for i, (a, b) in enumerate(zip(yo, yn)) if a != b] + ([("len", len(yo), len(yn))] if len(yo) != len(yn) else [])
        jd = _diff(_norm(json.loads(Path(str(o) + ".card_log.json").read_text()), roots),
                   _norm(json.loads(Path(str(n) + ".card_log.json").read_text()), roots))
        eo, en = o.with_name(o.stem + "_euler.yaml"), n.with_name(n.stem + "_euler.yaml")
        ed = None if not (eo.exists() or en.exists()) else (
            [_norm_str(x, roots) for x in eo.read_text().splitlines()] == [_norm_str(x, roots) for x in en.read_text().splitlines()])
        log_new = json.loads(Path(str(n) + ".card_log.json").read_text())
        same = not yd and not jd and ed in (None, True) and "M_knot_rule" not in log_new
        ok_all &= same
        print(f"V3 {name}: {'同じ' if same else '違いあり'} (YAML {len(yo)}/{len(yn)} 行・違う行 {len(yd)}、card_log の違い {len(jd)}、"
              f"Euler の写し {ed}、M_knot_rule の記録 {'あり' if 'M_knot_rule' in log_new else 'なし'})", flush=True)
        for d in (yd[:5] + jd[:10]):
            print("   ", str(d)[:300])
    print(f"V3 判定: {'PASS' if ok_all else 'FAIL'} (出力 {work})")
    return ok_all


# --- V4 ---------------------------------------------------------------------------------------------------------------
def check(out: Path, reeval: bool) -> bool:
    sys.path.insert(0, str(DESIGN))
    from forge_design.evaluate import mknot_select as MS
    from forge_design.probdef import load_problem
    raw = yaml.safe_load(out.read_text())
    logd = json.loads(Path(str(out) + ".card_log.json").read_text())
    mkr = raw["_card"]["M_knot_rule"]
    M = raw["geometry"]["M_knot"]
    N = int(raw["geometry"]["n_axis_inv"])
    dims = logd["dimensions"]
    load_problem(out)
    checks = {
        "寸法ゲート PASS (長さ 1 mm・径 0.02 mm)": dims["pass"] and raw["_card"]["tolerances"] == {"length_m": 1e-3, "diameter_m": 2e-5}
        and all(r["pass"] for r in dims["rows"]),
        "記録の選んだ値 = YAML の geometry.M_knot = card_log": mkr["chosen"] == M == logd["M_knot_rule"]["chosen"] == logd["design_chain"]["M_knot"],
        "全部の制約 (記録の N・2N)": mkr["constraints_at_output"]["N"]["ok"] and mkr["constraints_at_output"]["2N"]["ok"],
        "θ_max (記録 = 最終の設計チェーン)": mkr["theta_max_deg"] == mkr["theta_max_deg_final_chain"],
        "全部の境目からの距離 ≥ 0.05": mkr["distance_to_boundary"] is None or mkr["distance_to_boundary"] >= 0.05 - 1e-9,
        "解像度のゲート PASS": mkr["gate"]["pass"] and all(mkr["gate"]["checks"].values()),
        "規則の名前": mkr["rule"] == "min_theta_max_v1",
        "r_t (YAML = 記録 = card_log)": raw["spec"]["r_throat"] == mkr["r_t_m"] == logd["r_t"]["r_t_m"],
    }
    if reeval:
        base = {k: v for k, v in raw.items() if k != "_card"}
        rN = MS.evaluate_mknot(base, M, N, base_dir=out.parent)
        r2 = MS.evaluate_mknot(base, M, 2 * N, base_dir=out.parent)
        checks["全部の制約 (出力 YAML を独立に評価し直す N・2N)"] = rN["ok"] and r2["ok"]
        checks["θ_max (評価し直した N = 記録)"] = rN["theta_max_deg"] == mkr["theta_max_deg"]
        print(f"  評価し直し: N {N} θ_max {rN['theta_max_deg']!r} 成立 {rN['ok']}・2N {2 * N} θ_max {r2['theta_max_deg']!r} 成立 {r2['ok']}")
    print(f"V4 {out.name}: {'PASS' if all(checks.values()) else 'FAIL'}")
    for k, v in checks.items():
        print(f"  {'ok  ' if v else 'FAIL'} {k}")
    g = mkr["gate"]
    print(f"  M_knot {M!r}、r_t {raw['spec']['r_throat']!r} m、θ_max {mkr['theta_max_deg']:.6f}°、余裕 {mkr['margin_deg']:.4f}°、"
          f"最小 sin 角 {mkr['min_sin_angle']:.4f}、距離 {mkr['distance_to_boundary']}、境目 "
          f"{[(b['kind'], round(b['feasible_end'], 6), b['source']) for b in mkr['boundaries']]}")
    print(f"  ゲート: N {g['chosen_N']!r}・2N {g['chosen_2N']!r}、入力 {g['chosen']!r}、窓 {g['window']}、最初の r_t のゲートを再利用 "
          f"{mkr['gate_first_reused']}、やり直し {mkr['redo']}")
    for e in mkr["events"]:
        print(f"  [{e['stage']}] r_t {e['r_t_m']!r} 入力 {e['M_K_in']!r} → chosen_N {e['chosen_N']!r} (差 {e['diff_vs_used']}): "
              f"{e['action'][:160]} ({e['elapsed_s']:.0f} s)")
    for r in dims["rows"]:
        print(f"  寸法 {r['name']:<26} 差 {r['diff_m']:+.3e} m (公差 {r['tol_m']:g}) {'ok' if r['pass'] else 'FAIL'}")
    print(f"  所要: 全体 {logd['time_s']:.0f} s、うち M_knot {mkr['elapsed_s']:.0f} s (並列 {mkr['workers']})")
    return all(checks.values())


def run_v4(which: str) -> bool:
    card = V4DIR / f"card_c2b_auto_{which}.yaml"
    out = V4DIR / f"problem_c2b_auto_{which}.yaml"
    for p in (out, Path(str(out) + ".card_log.json"), out.with_name(out.stem + "_euler.yaml")):
        p.unlink(missing_ok=True)
    t0 = time.time()
    r = subprocess.run([sys.executable, "-m", "forge_design.intake.wt_card", str(card), "--out", str(out)], cwd=DESIGN,
                       capture_output=True, text=True)
    wall = time.time() - t0
    (V4DIR / f"problem_c2b_auto_{which}.stdout.txt").write_text(r.stdout + r.stderr + f"\nrc={r.returncode} wall {wall:.0f} s\n")
    print(f"V4 {which}: 生成 rc {r.returncode}、{wall:.0f} s ({V4DIR / f'problem_c2b_auto_{which}.stdout.txt'})", flush=True)
    return r.returncode == 0 and check(out, reeval=False)


if __name__ == "__main__":
    argv = list(sys.argv[1:])
    old_rev = _arg(argv, "--old-rev", "d9f554f4")
    outdir = Path(_arg(argv, "--out", str(C / "_mknot_card_v3")))
    reeval = "--reeval" in argv
    if reeval:
        argv.remove("--reeval")
    what = argv[0] if argv else ""
    if what == "v3":
        ok = run_v3(old_rev, outdir)
    elif what == "v4" and len(argv) == 2 and argv[1] in ("final", "explore"):
        ok = run_v4(argv[1])
    elif what == "check" and len(argv) == 2:
        ok = check(Path(argv[1]).resolve(), reeval)
    else:
        sys.exit(__doc__)
    sys.exit(0 if ok else 1)
