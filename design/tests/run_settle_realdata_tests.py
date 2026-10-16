#!/usr/bin/env python3
r"""result 段レビュー 2 回目 (plan tooling-settle-convergence-and-monitor §5.1 #9・#10・#11; notes/reviews/
2026-10-16-tooling-settle-convergence-and-monitor-result-2.md の 1〜3) の再現と修正の確かめを、手元の V4・V6 の写し
(`case/71.m55_combustion_wt/run_0011_ns_v2`・`run_0014_chain_v2`) で行う。合成の run の試験は `solver_density_cuda/tools/test_settle.py`
(S11b〜S11f・F1〜F10) と `design/tests/run_settle_migration_tests.py` ((e′)) にある。GPU・forge は使わない。

run は一時ディレクトリへ写して使う (h5・xmf は symlink、それ以外は実ファイルの写し — 書き換える負例が元の run を上書きしないように;
元の run には何も書かない)。run が無ければその run の項目を飛ばす (飛ばしたことを出す)。

  (#9) 基準の 1 項目を差し替えた写し (grid.segment.stage main → S1) と基準の sha256 の不一致 → `check_settled.py` (既定の照合) が
       INDETERMINATE・終了コード 3。差し替えない写しは今どおり SETTLED・0。修正前 (PRE_REVIEW2_COMMIT の check_settled.py) は S1 の写しで
       SETTLED・0 だった (検出力)。
  (#10) 生成器 (`settle_presets.resolve`) が run の入力から解決した wt_nozzle_v3 の基準で、実データは SETTLED (見張りは毎 step の行・
       既定の 4 列を満たす)、残差を最終 1 行だけ・rms_ro の列だけ・残差の列の無い step,inner,phase の 1 行にした写しは INDETERMINATE
       (見張りが判定不能)。修正前の settle.py はどれも SETTLED だった (検出力)。
  (#11) 成果物の収束 (`wall_tables.thermal_use`) が、同じ run を相対・絶対のパスで渡しても、別の numpy (system python3 があれば) で
       判定しても、停止規則 v2 の再評価と一致して PASS (修正前は「未評価」; 検出力)。

実行 (作業ツリーの直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_settle_realdata_tests.py
所要: 約 1 分。
"""
from __future__ import annotations

import importlib.util
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
from contextlib import redirect_stdout
from pathlib import Path

import yaml

DESIGN = Path(__file__).resolve().parents[1]
ROOT = DESIGN.parent
TOOLS = ROOT / "solver_density_cuda" / "tools"
sys.path.insert(0, str(DESIGN))
sys.path.insert(0, str(TOOLS))
os.environ.setdefault("FORGE_NO_DESIGN_REPORT", "1")

from forge_design.evaluate import settle_presets as SP  # noqa: E402
from forge_design.evaluate import stop_rule as SR  # noqa: E402
from forge_design.export import wall_tables as WT  # noqa: E402
import check_settled  # noqa: E402

ST = SR._settle()
CASE = ROOT / "case" / "71.m55_combustion_wt"
RUNS = {"V4 run_0011_ns_v2": CASE / "run_0011_ns_v2", "V6 run_0014_chain_v2": CASE / "run_0014_chain_v2"}
PRE_REVIEW2_COMMIT = "b64dbfdd"      # result 段レビュー 2 回目の修正の前 (検出力の確かめ)
SYS_PY = Path("/usr/bin/python3")    # design/.venv-opt と numpy の版が違う (1.26 と 2.x; 再計算の丸めが違う)
TMP = Path(tempfile.mkdtemp(prefix="settle_realdata_"))
FAIL = 0
NCHK = 0
SKIPPED = []


def check(name, cond, detail=""):
    global FAIL, NCHK
    NCHK += 1
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


def copy_run(src: Path, name: str) -> Path:
    """run の写し: h5・xmf は symlink、それ以外のファイルは実ファイルの写し、ディレクトリは写さない。"""
    dst = TMP / name
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir()
    for f in src.iterdir():
        if f.is_dir():
            continue
        if f.suffix in (".h5", ".xmf"):
            (dst / f.name).symlink_to(f.resolve())
        else:
            shutil.copy2(f, dst / f.name)
    return dst


def load_at(commit: str, relpath: str, modname: str, package: str | None = None):
    """commit の relpath を git から読んでモジュールにする。package を渡すと同じパッケージの中の名前で読む (相対 import を解く)。"""
    src = subprocess.run(["git", "-C", str(ROOT), "show", f"{commit}:{relpath}"], capture_output=True, text=True, check=True).stdout
    d = TMP / "ref" / Path(relpath).parent
    d.mkdir(parents=True, exist_ok=True)
    f = d / (modname.rsplit(".", 1)[-1] + ".py")
    f.write_text(src)
    spec = importlib.util.spec_from_file_location(modname, f)
    m = importlib.util.module_from_spec(spec)
    if package:
        sys.modules[modname] = m
    spec.loader.exec_module(m)
    return m


def cli(mod, run: Path, *args):
    buf = io.StringIO()
    with redirect_stdout(buf):
        try:
            rc = mod.main([str(run), *args])
        except SystemExit as e:
            rc = e.code
    out = buf.getvalue()
    return rc, out, (out.strip().splitlines() or [""])[-1]


def rewrite_residual(run: Path, header: list, rows: list):
    p = run / "residual_history.csv"
    p.unlink()                          # 写しの実ファイルを消してから書く (symlink を通して元を書き換えない)
    with open(p, "w") as f:
        f.write(",".join(header) + "\n")
        for r in rows:
            f.write(r + "\n")


def v3_criteria(run: Path) -> Path:
    """生成器が run の入力 (prepare_info・solverConfig・bcondConfig) から解決した wt_nozzle_v3 の基準 (x 範囲は run の基準と同じ)。"""
    pi = json.loads((run / "prepare_info.json").read_text())
    sc = yaml.safe_load((run / "solverConfig.yaml").read_text())
    bc = yaml.safe_load((run / "bcondConfig.yaml").read_text())
    cur = yaml.safe_load((run / ST.CRITERIA_FILE).read_text())
    es = {"preset": "wt_nozzle_v3"}
    xr = (cur.get("source") or {}).get("evaluate_settle", {}).get("x_range_m")
    if xr:
        es["x_range_m"] = xr
    c = SP.resolve(es, pi, sc, bc, segment=SP.main_stage_segment())
    p = TMP / f"criteria_v3_{run.name}.yaml"
    p.write_text(yaml.safe_dump(c, allow_unicode=True, sort_keys=False))
    return p


def main() -> int:
    print(f"出力: {TMP}", flush=True)
    try:
        old_cli = load_at(PRE_REVIEW2_COMMIT, "solver_density_cuda/tools/check_settled.py", "check_settled_pre_review2")
        old_settle = load_at(PRE_REVIEW2_COMMIT, "solver_density_cuda/tools/settle.py", "settle_pre_review2")
        old_wt = load_at(PRE_REVIEW2_COMMIT, "design/forge_design/export/wall_tables.py", "forge_design.export._wall_tables_pre_review2",
                         package="forge_design.export")
        for k_ in ("REPO", "TOOLS"):
            setattr(old_wt, k_, getattr(WT, k_))
    except Exception as e:  # noqa: BLE001
        old_cli = old_settle = old_wt = None
        check(f"修正前 ({PRE_REVIEW2_COMMIT}) のモジュールを git から読めた (読めなければ検出力は判定不能 — 合格にしない)", False, str(e))
    for label, src in RUNS.items():
        if not (src / ST.CRITERIA_FILE).is_file() or not (src / SR.STATE).is_file():
            SKIPPED.append(f"{label} ({src.relative_to(ROOT)} が無い)")
            print(f"skip {label}: {src.relative_to(ROOT)} が無い", flush=True)
            continue
        print(f"== {label} ({src.relative_to(ROOT)})", flush=True)
        tag = src.name
        # ---------------------------------------------------------------- #9
        base = copy_run(src, f"{tag}_base")
        rc, out, last = cli(check_settled, base, "--no-write")
        check(f"(#9) {label}: 差し替えない写し → 今どおり SETTLED・終了コード 0", rc == 0 and last.startswith("VERDICT: SETTLED"), last)
        r = copy_run(src, f"{tag}_stage_s1")
        txt = (r / ST.CRITERIA_FILE).read_text()
        check(f"(#9) {label}: 基準に grid.segment.stage: main の行が 1 つある (差し替えの前提)", txt.count("\n    stage: main\n") == 1)
        (r / ST.CRITERIA_FILE).write_text(txt.replace("\n    stage: main\n", "\n    stage: S1\n"))
        rc, out, last = cli(check_settled, r, "--no-write")
        check(f"(#9) {label}: 基準の 1 項目を差し替え (stage main → S1) → INDETERMINATE (基準が本段の判定に使えない)・終了コード 3",
              rc == 3 and last == "VERDICT: INDETERMINATE (基準が本段の判定に使えない)", last)
        if old_cli is not None:
            rc_o, _, last_o = cli(old_cli, r, "--no-write")
            check(f"(#9) {label}: 検出力 — 修正前の check_settled.py は同じ写しで SETTLED・終了コード 0 (レビューの再現)",
                  rc_o == 0 and last_o.startswith("VERDICT: SETTLED"), last_o)
        rc, out, last = cli(check_settled, r, "--no-write", "--no-bind")
        check(f"(#9) {label}: --no-bind を明示 → 今どおり照合を外して判定 (振る舞いを変えない)", rc == 0 and "照合を外した" in last, last)
        r = copy_run(src, f"{tag}_sha")
        with open(r / ST.CRITERIA_FILE, "a") as f:
            f.write("# 本段の起動の後に書き換えた\n")
        rc, out, last = cli(check_settled, r, "--no-write")
        check(f"(#9) {label}: 基準の sha256 が本段の記録 (MAIN_STAGE_STARTED) と違う → INDETERMINATE (照合の不一致)・終了コード 3",
              rc == 3 and last == "VERDICT: INDETERMINATE (本段の起動との照合の不一致)" and "settle_criteria.sha256" in out, last)
        # ---------------------------------------------------------------- #10
        cp = v3_criteria(base)
        lc = ST.load_criteria(base, cp)
        res = ST.judge(base, None, loaded=lc)
        q = next(x for x in res["quantities"] if x["kind"] == "residual_finite")
        check(f"(#10) {label}: 生成器の wt_nozzle_v3 の基準で実データ → SETTLED (見張りは毎 step の行 {q['values']['rows_in_segment']} 行・"
              f"必要な列 {q['values']['required']})", res["verdict"] == "SETTLED" and q["verdict"] == "SETTLED"
              and q["values"]["rows_in_segment"] == res["end_step"], f"{res['verdict']} {[(x['name'], x['verdict']) for x in res['quantities']]}")
        lines = (base / "residual_history.csv").read_text().splitlines()
        hdr = lines[0].split(",")
        outer = [ln for ln in lines[1:] if ln.split(",")[2] == "outer_end"]
        i_ro = hdr.index("rms_ro")
        negs = (("最終 1 行だけ", hdr, outer[-1:]),
                ("rms_ro の列だけ", ["step", "inner", "phase", "rms_ro"], [",".join(ln.split(",")[:3] + [ln.split(",")[i_ro]]) for ln in outer]),
                ("残差の列の無い step,inner,phase の 1 行", ["step", "inner", "phase"], [",".join(outer[-1].split(",")[:3])]))
        for k, (what, h, rows) in enumerate(negs):
            r = copy_run(src, f"{tag}_m2_{k}")
            rewrite_residual(r, h, rows)
            res = ST.judge(r, None, loaded=ST.load_criteria(r, cp))
            q = next(x for x in res["quantities"] if x["kind"] == "residual_finite")
            check(f"(#10) {label}: 残差を{what}にした写し → 見張りが判定不能・全体 INDETERMINATE", q["verdict"] == "INDETERMINATE"
                  and res["verdict"] == "INDETERMINATE", str([t["text"][:120] for t in q["reasons"]]))
            if old_settle is not None:
                # 修正前の settle.py は residual_finite の `required` を知らない (知らないキーで例外) — 生成器が書くようになった
                # (plan §5.1 #10 の続き) ので、修正前のコードに渡す基準からだけ外す (検出力の確かめは列の要求ではなく行・列の欠けを見る)
                cp_old = TMP / f"{cp.stem}_pre_review2.yaml"
                c_old = yaml.safe_load(cp.read_text())
                for q_ in c_old["quantities"]:
                    if q_.get("kind") == "residual_finite":
                        q_.pop("required", None)
                cp_old.write_text(yaml.safe_dump(c_old, allow_unicode=True, sort_keys=False))
                ro = old_settle.judge(r, None, loaded=old_settle.load_criteria(r, cp_old))
                check(f"(#10) {label}: 検出力 — 修正前の settle.py は残差を{what}にした写しで SETTLED (レビューの再現)", ro["verdict"] == "SETTLED",
                      ro["verdict"])
        # ---------------------------------------------------------------- #11
        cwd = os.getcwd()
        try:
            os.chdir(base.parent)
            th_rel = WT.thermal_use(Path(base.name))
            th_rel_old = old_wt.thermal_use(Path(base.name)) if old_wt is not None else None
        finally:
            os.chdir(cwd)
        th_abs = WT.thermal_use(base.resolve())
        for how, th in (("相対のパス", th_rel), ("絶対のパス", th_abs)):
            c = th["convergence"]
            check(f"(#11) {label}: 成果物の収束 ({how}) → 停止規則 v2 の再評価と一致・収束 PASS (総合 {th['verdict']})",
                  c["result"] == "PASS" and c["stop_rule_reeval"].get("matches") is True, str(c.get("reason"))[:300])
        if th_rel_old is not None:
            c = th_rel_old["convergence"]
            check(f"(#11) {label}: 検出力 — 修正前の wall_tables は同じ写しで収束「未評価」(レビューの再現)",
                  c["result"] == "未評価" and c["stop_rule_reeval"].get("matches") is False, str(c.get("reason"))[:200])
        if SYS_PY.is_file():
            th = WT.thermal_use(base.resolve(), python=str(SYS_PY))
            c = th["convergence"]
            m = c["stop_rule_reeval"].get("compare") or {}
            check(f"(#11) {label}: 別の numpy ({SYS_PY}) で判定しても一致・収束 PASS", c["result"] == "PASS" and m.get("matches") is True,
                  str(c.get("reason"))[:300])
        else:
            SKIPPED.append(f"(#11) {label}: {SYS_PY} が無い (別の numpy の照合)")
    if SKIPPED:
        print("飛ばした項目: " + "; ".join(SKIPPED), flush=True)
    print(f"\n{NCHK - FAIL}/{NCHK} ok — {'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'} (一時ディレクトリ {TMP})", flush=True)
    if FAIL == 0:
        shutil.rmtree(TMP, ignore_errors=True)
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
