#!/usr/bin/env python3
"""判定器の不正な入力の試験 (plan time_integration-line-viscous-jacobian-dt・faceh の result 段 codex M3・m4・faceh M1、2026-10-11)。

模擬の入力で、直した判定器が不正な入力を合格にしないこと、正常な入力を通すことを確かめる。
  1. lvc_segcheck.recheck: 正常 → ok、残差 CSV なし・CSV の末尾を削る (last step 不一致)・保存の VERDICT を書き換える → NG
  2. lvcmp_judge.qs: 正常な 4 量 → ok、列を 1 つ落とす → ok=False (ツールは rc 1 + OVERALL を出すので、分類の行で見分ける)
  3. lvcmp_judge.assess: 見張りの途中の行に NaN・Inf・量の欠け → INVALID (到達の行は正常のまま)
  4. lvcmp_rerun_judge.floor_median: 残差の step 0〜149999 の CSV で、窓の終わり 149999 は通り、150000 は落ちる
  5. lvcaudit_judge (引数で保存済みの証拠の写しを渡したときだけ): 監査用の rhs_s1 に NaN → INVALID、
     対照 2 本の rhs_s1 に同じビット列の NaN → INVALID、scalar_visc_line の形の不一致 → INVALID
使い方: lvc_judge_badinput_test.py [<run_0550・0551・0552 の写しがある作業用ディレクトリ> <line_audit_helper>]
     作業用ディレクトリの中身は試験で書き換えるので、証拠そのものを渡さないこと。
"""
import csv
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
argv = sys.argv[1:]
sys.argv = sys.argv[:1]
import lvc_segcheck as SC  # noqa: E402
import lvcmp_judge as J  # noqa: E402
import lvcmp_rerun_judge as R  # noqa: E402

fails = []


def expect(name, cond):
    print(("ok   " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


def write_hist(d, n_steps, keys=("rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe", "rms_roK", "rms_roOmega")):
    with open(d / "residual_history_main.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["step", "phase"] + list(keys))
        for s in range(n_steps):
            v = [1e-3 * (0.999 ** s)] * len(keys)
            v[3] = 0.0
            w.writerow([s, "outer_begin"] + v); w.writerow([s, "outer_end"] + v)


def seg_run(t, n_steps):
    d = t / "segrun"; d.mkdir()
    (d / "stage_manifest.json").write_text(json.dumps({"manifest_version": 2, "stages": [
        {"order": 0, "tag": "main", "history": "residual_history_main.csv", "restart_from": None, "key": {"solver": "SLAU"}, "cfg_fnv": "x"}]}))
    write_hist(d, n_steps)
    o = subprocess.run([sys.executable, str(SC.TOOL), str(d), "--segment"], capture_output=True, text=True, cwd=d)
    (d / "CONVERGENCE_SEGMENT.txt").write_text(o.stdout + o.stderr)
    (d / "residual_history_segment.csv").unlink(missing_ok=True)
    return d


with tempfile.TemporaryDirectory() as t_:
    t = Path(t_)
    # ---- 1. lvc_segcheck ----
    d = seg_run(t, 3000)
    ok, rec = SC.recheck(d)
    expect(f"segcheck: 正常な run は通る ({rec.get('why')})", ok)
    d2 = t / "seg_nocsv"; shutil.copytree(d, d2); (d2 / "residual_history_main.csv").unlink()
    expect("segcheck: 残差 CSV なしは NG", not SC.recheck(d2)[0])
    d3 = t / "seg_trunc"; shutil.copytree(d, d3)
    rows = (d3 / "residual_history_main.csv").read_text().splitlines()
    (d3 / "residual_history_main.csv").write_text("\n".join(rows[:-200]) + "\n")
    expect("segcheck: CSV の末尾を削った (last step 不一致) は NG", not SC.recheck(d3)[0])
    d4 = t / "seg_text"; shutil.copytree(d, d4)
    txt = (d4 / "CONVERGENCE_SEGMENT.txt").read_text()
    stored = SC.parse(txt)[2]
    other = "CONVERGED" if stored != "CONVERGED" else "NOT CONVERGED"
    (d4 / "CONVERGENCE_SEGMENT.txt").write_text(txt.replace(f"-> {stored}", f"-> {other}"))
    expect("segcheck: 保存の VERDICT を書き換えたものは NG", not SC.recheck(d4)[0])

    # ---- 2. lvcmp_judge.qs ----
    J.D = t
    win = [{"step": s, "theta_r_40": 0.0778, "theta_r_70": 0.1149, "theta_r_94": 0.1277, "Q_w": 9.29e6} for s in range(90000, 110001, 2500)]
    q = J.qs(win, "test_ok")
    expect(f"qs: 正常な 4 量は ok ({q.get('classes')})", q.get("ok") is True)
    real_keys = J.KEYS
    J.KEYS = real_keys[:3] + ("Q_w_missing",)
    try:
        q2 = J.qs([{**x, "Q_w_missing": x["Q_w"]} for x in win], "test_bad")
    finally:
        J.KEYS = real_keys
    # 列名を変えた CSV を書いて、元の KEYS で照合するとツールは「列の欠け」で終わる
    p = t / "lvcmp_qs_test_bad.csv"
    o = subprocess.run([sys.executable, str(HERE / "../../solver_density_cuda/tools/check_quasisteady.py"), "--series-csv", str(p), "--series-cols", ",".join(real_keys),
                        "--tail", "1.0", "--drift", "0.0005", "--min-snaps", "5"], capture_output=True, text=True)
    txt = o.stdout + o.stderr
    expect(f"qs: ツールは列の欠けでも rc {o.returncode}・OVERALL を出す (前提の確認)", o.returncode == 1 and "OVERALL" in txt)
    J.KEYS = real_keys

    def qs_on(csvp):
        o = subprocess.run([sys.executable, str(HERE / "../../solver_density_cuda/tools/check_quasisteady.py"), "--series-csv", str(csvp), "--series-cols", ",".join(J.KEYS),
                            "--tail", "1.0", "--drift", "0.0005", "--min-snaps", "5"], capture_output=True, text=True)
        return o
    # qs() 本体の判定の部分を同じ出力で確かめるため、列の欠けた CSV を qs の出力先に置いて呼ぶ
    orig_run = J.subprocess.run

    def fake_run(cmd, **kw):
        return orig_run(cmd[:3] + [str(p)] + cmd[4:], **kw)   # --series-csv の値だけを差し替える
    J.subprocess.run = fake_run
    try:
        q3 = J.qs(win, "test_redirect")
    finally:
        J.subprocess.run = orig_run
    expect(f"qs: 列の欠けた入力は ok=False ({q3.get('classes')}, rc {q3.get('rc')})", q3.get("ok") is False)

    # ---- 3. lvcmp_judge.assess の途中の行 ----
    J.HERE = t
    for tag, bad in (("nan", float("nan")), ("inf", float("inf")), ("none", None)):
        run = f"run_test_{tag}"
        (t / run).mkdir()
        rows = [{"step": s, "deficit": 0.05, "theta_r_40": 0.0778, "theta_r_70": 0.1149, "theta_r_94": 0.1277, "Q_w": 9.29e6} for s in range(2500, 110001, 2500)]
        rows[10]["Q_w" if tag != "none" else "theta_r_40"] = bad
        if tag == "none":
            del rows[10]["theta_r_40"]
        (t / run / "m9_watch.json").write_text(json.dumps({"run": run, "phase": "line", "budget": J.BUDGET, "consec": 1, "status": "REACHED",
                                                           "reach_step": 110000, "rows": rows}, allow_nan=True))
        J.ARMS["T_" + tag] = (run, "7", "m7")
        cls, info = J.assess("T_" + tag)
        expect(f"assess: 途中の行の {tag} は INVALID ({info.get('why')})", cls == "INVALID" and "見張りの行の量" in str(info.get("why")))

    # ---- 4. lvcmp_rerun_judge.floor_median の窓 ----
    R.HERE = t
    run = "run_test_floor"; (t / run).mkdir()
    with open(t / run / "residual_history.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["step", "phase"] + list(R.RES_COLS))
        for s in range(J.BUDGET - 25000, J.BUDGET):
            w.writerow([s, "outer_begin"] + [1e-5] * len(R.RES_COLS))
    try:
        fm = R.floor_median(run, J.BUDGET - 1)
        expect(f"floor_median: 窓の終わり {J.BUDGET - 1} は通る ({fm['rows']} 行)", fm["rows"] == 20001)
    except Exception as e:
        expect(f"floor_median: 窓の終わり {J.BUDGET - 1} は通る ({e})", False)
    try:
        R.floor_median(run, J.BUDGET)
        expect(f"floor_median: 窓の終わり {J.BUDGET} は落ちる", False)
    except ValueError:
        expect(f"floor_median: 窓の終わり {J.BUDGET} は落ちる", True)

# ---- 5. lvcaudit_judge (保存済みの証拠の写しがあるときだけ) ----
if len(argv) >= 2:
    base, helper = Path(argv[0]).resolve(), str(Path(argv[1]).resolve())
    SHA_A = "eb45557958c6547f9ab7b29e17bdec7f82dc94466736b7c02e708efb5ae6452a"
    SHA_C = "4910b74d0e7c6ae62644f2a7262578f8aab7388cdacfbb842d66be919291391b"

    def run_judge(root):
        subprocess.run([sys.executable, str(HERE / "lvcaudit_judge.py"), str(root / "run_0550_lvcaudit/linedump"), helper,
                        "--run", str(root / "run_0550_lvcaudit"), "--expect-sha", SHA_A, "--ctrl", str(root / "run_0551_lvcaudit_ctrl/linedump"),
                        "--ctrl-run", str(root / "run_0551_lvcaudit_ctrl"), "--ctrl-sha", SHA_C, "--ctrl2", str(root / "run_0552_lvcaudit_ctrl2/linedump")],
                       capture_output=True, text=True)
        return json.loads((root / "run_0550_lvcaudit/linedump/lvcaudit_judge.json").read_text())["VERDICT"]

    def case(tag, mutate):
        root = Path(tempfile.mkdtemp(prefix=f"lvcaudit_{tag}_"))
        try:
            for r in ("run_0550_lvcaudit", "run_0551_lvcaudit_ctrl", "run_0552_lvcaudit_ctrl2"):
                shutil.copytree(base / r, root / r, ignore=shutil.ignore_patterns("res_*", "*.xmf", "nozzle.h5", "mesh.h5"))
            mutate(root)
            return run_judge(root)
        finally:
            shutil.rmtree(root, ignore_errors=True)

    def nan_at(path, k=7):
        a = np.fromfile(path, dtype=np.float64); a[k] = np.nan; a.tofile(path)

    v = case("normal", lambda r: None)
    expect(f"lvcaudit: 正常な写しは前と同じ FAIL ({v[:40]})", v.startswith("FAIL"))
    v = case("auditnan", lambda r: nan_at(r / "run_0550_lvcaudit/linedump/rhs_s1.f64"))
    expect(f"lvcaudit: 監査用の rhs_s1 の NaN は INVALID ({v[:60]})", v.startswith("INVALID"))
    v = case("ctrlnan", lambda r: [nan_at(r / x / "linedump/rhs_s1.f64") for x in ("run_0551_lvcaudit_ctrl", "run_0552_lvcaudit_ctrl2")])
    expect(f"lvcaudit: 対照 2 本の同じビット列の NaN は INVALID ({v[:60]})", v.startswith("INVALID"))

    def bad_shape(r):
        m = r / "run_0550_lvcaudit/linedump/meta.txt"
        m.write_text(m.read_text().replace("scalar_visc_line 605 1", "scalar_visc_line 121 5"))
    v = case("shape", bad_shape)
    expect(f"lvcaudit: scalar_visc_line の形の不一致は INVALID ({v[:60]})", v.startswith("INVALID"))
else:
    print("(lvcaudit_judge の試験は飛ばした: 証拠の写しのディレクトリと line_audit_helper を渡していない)")

print(f"\n{'ALL OK' if not fails else f'{len(fails)} 件 FAIL'}")
sys.exit(1 if fails else 0)
