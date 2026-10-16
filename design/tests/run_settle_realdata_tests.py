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
  (#16) result 段レビュー 3 回目 M1 (notes/reviews/2026-10-16-tooling-settle-convergence-and-monitor-result-3.md の 1): 写しで
       「通常の CLI (`check_settled.py`、相対のパス・既定で SETTLE_VERDICT.json を上書き) → 成果物の収束 (`thermal_use`、相対・絶対の
       パス) → 受け渡しの照合 (`deltastar_loop._v2_binding`)」を続けて回し、成果物が「未評価」にならず照合も通る。古い形 (写しのまま;
       停止の時の専用ファイルが無い V4・V6 の形 → CLI が上書きの直前に停止の時の原本を同じ名前の専用ファイルへ移す) と新しい形 (写しの
       SETTLE_VERDICT.json を停止の時の専用ファイル `stop_rule_v2/SETTLE_VERDICT_stop_<step>.json` に置き、状態の記録に相対パスと
       sha256 を足した — 修正後の見張りが書く形) の両方。専用ファイル (新しい形・移した古い形) を書き換えると未評価・照合の失敗。
       修正前 (PRE_M1_COMMIT の stop_rule・deltastar_loop・wall_tables) は古い形で CLI の後に「未評価」・照合の失敗 (レビューの再現; 検出力)。
  (#16 続き) result 段レビュー 4 回目 M1 (notes/reviews/2026-10-16-tooling-settle-convergence-and-monitor-result-4.md の 1): 写しで
       流量の CSV の流入・流出の絶対値を各 1 kg/s 増やして通常の CLI (終了コード 0・SETTLED) → 成果物の収束は未評価・受け渡しの照合は
       失敗 (古い形・新しい形の両方; 入力の sha256 が停止の時の記録と違う)。古い形で移行の前に SETTLE_VERDICT.json が上書きされていた写し
       (修正前の CLI で上書き; 原本なし) は照合不能で拒否。修正前 (PRE_M1B_COMMIT = 61fc3d59 の stop_rule・deltastar_loop・wall_tables・
       check_settled) は、入力を変えて通常の CLI を回すと古い形の成果物 (記録の形 overwritten)・受け渡し、新しい形の受け渡しが通った
       (レビューの再現; 検出力)。

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
from forge_design.feedback import deltastar_loop as DL  # noqa: E402
import check_settled  # noqa: E402

ST = SR._settle()
CASE = ROOT / "case" / "71.m55_combustion_wt"
RUNS = {"V4 run_0011_ns_v2": CASE / "run_0011_ns_v2", "V6 run_0014_chain_v2": CASE / "run_0014_chain_v2"}
PRE_REVIEW2_COMMIT = "b64dbfdd"      # result 段レビュー 2 回目の修正の前 (検出力の確かめ)
PRE_M1_COMMIT = "8cf2a284"           # result 段レビュー 3 回目 M1 (§5.1 #16) の修正の前 (検出力の確かめ)
PRE_M1B_COMMIT = "61fc3d59"          # result 段レビュー 4 回目 M1 (§5.1 #16 の続き) の修正の前 (検出力の確かめ)
PRE_M5_COMMIT = "b283aa48"           # result 段レビュー 5 回目 M1 (旧形式の residual_finite の正規化) の修正の前 (検出力の確かめ)
# 手順書の試走 A (#10 より前に書いた v3 の停止の時の記録; residual_finite の values に required・required_source が無い)。手元の写し
TRIAL = {"試走 A p0 run_0003_ns_p0": ROOT / "case/90.qs_trial_a/run_0003_ns_p0", "試走 A p1 run_0005_ns_p1": ROOT / "case/90.qs_trial_a/run_0005_ns_p1"}
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


def rejudge_cli(run: Path, mod=None):
    """通常の CLI (check_settled.py; 既定で SETTLE_VERDICT.json に書く) を run の親から相対のパスで回す。mod: 別の版の check_settled
    (検出力の確かめ; 既定は今の作業ツリー)。戻り (終了コード, 最後の行, 出力の全文)。"""
    cwd = os.getcwd()
    try:
        os.chdir(run.parent)
        rc, out, last = cli(mod or check_settled, Path(run.name))
    finally:
        os.chdir(cwd)
    return rc, last, out


def bump_massflux(run: Path, dm: float = 1.0):
    """流量の CSV の流入・流出の絶対値を各 dm [kg/s] 増やす (result 段レビュー 4 回目 M1 の再現: 流入は負・流出は正の符号; 実寸の
    mdot_ と生の raw_ (= mdot/2π) の両方)。写しの実ファイルを消してから書く。"""
    import math
    p = run / "boundary_massflux.csv"
    lines = p.read_text().splitlines()
    hi = next(i for i, ln in enumerate(lines) if ln.startswith("step,"))
    cols = lines[hi].split(",")
    step = {}
    for c in cols:
        if c.startswith(("mdot_inlet", "mdot_outlet")):
            step[c] = (-1.0 if "inlet" in c else 1.0) * dm
        elif c.startswith(("raw_inlet", "raw_outlet")):
            step[c] = (-1.0 if "inlet" in c else 1.0) * dm / (2.0 * math.pi)
    out = lines[:hi + 1]
    for ln in lines[hi + 1:]:
        if not ln or ln.startswith("#"):
            out.append(ln)
            continue
        v = ln.split(",")
        for c, d in step.items():
            j = cols.index(c)
            v[j] = repr(float(v[j]) + d)
        out.append(",".join(v))
    p.unlink()
    p.write_text("\n".join(out) + "\n")


def thermal_both(run: Path, mod=None):
    """成果物の収束 (`thermal_use`) を相対のパス・絶対のパスで。戻り (相対の収束, 絶対の収束)。"""
    mod = mod or WT
    cwd = os.getcwd()
    try:
        os.chdir(run.parent)
        c_rel = mod.thermal_use(Path(run.name))["convergence"]
    finally:
        os.chdir(cwd)
    return c_rel, mod.thermal_use(run.resolve())["convergence"]


def to_new_form(run: Path) -> str:
    """写しを修正後の見張りが書く形にする: 停止の時のままの SETTLE_VERDICT.json を停止の時の専用ファイルに置き、stop_rule.json と
    STOP_RULE_STATE.json の reeval に相対パスと sha256 を足す (STOP_RULE_STATE の stop_rule.json の sha256 も合わせる)。戻り 相対パス。"""
    st = SR.read_state(run)
    ve = st["reeval"]
    assert SR._sha256(run / ST.VERDICT_FILE) == ve["sha256"], "SETTLE_VERDICT.json が停止の時のままでない"
    rel = SR.stop_verdict_relpath(ve["end_step"])
    (run / rel).parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(run / ST.VERDICT_FILE, run / rel)
    add = {"stop_file": rel, "stop_file_sha256": SR._sha256(run / rel)}
    rec = json.loads((run / SR.RECORD).read_text())
    rec["reeval"].update(add)
    SR.save_record(run, rec)
    st["reeval"].update(add)
    st["record_sha256"] = SR._sha256(run / SR.RECORD)
    SR.write_state(run, st)
    return rel


def with_stop_rule(mod, fn):
    """`from ..evaluate import stop_rule` を mod (修正前の stop_rule) に差し替えて fn() を呼ぶ (検出力の確かめ; 呼んだ後に戻す)。"""
    import forge_design.evaluate as EV
    key = "forge_design.evaluate.stop_rule"
    saved_m, saved_a = sys.modules[key], EV.stop_rule
    sys.modules[key], EV.stop_rule = mod, mod
    try:
        return fn()
    finally:
        sys.modules[key], EV.stop_rule = saved_m, saved_a


def test_legacy_v3(label: str, src: Path, old_sr):
    """(5 回目 M1) #10 より前の v3 の記録 (residual_finite の values に required・required_source が無い) を、成果物の収束が
    旧形式に限った正規化で照合できる (PASS)・入力を変えると拒否 (未評価)・修正前 (PRE_M5_COMMIT) の stop_rule では未評価 (検出力)。"""
    r = copy_run(src, f"{src.name}_legacy_v3")
    ve = SR.read_state(r)["reeval"]
    rec_q = json.loads((r / ST.VERDICT_FILE).read_text())["quantities"]
    rf = [q for q in rec_q if q.get("kind") == ST.KIND_RESIDUAL_FINITE]
    check(f"(5 回目 M1) {label}: 停止の時の記録は旧形式 (residual_finite の values に required・required_source が無い、v3)",
          len(rf) == 1 and not any(k in rf[0].get("values", {}) for k in SR.LEGACY_RF_KEYS), str([list(q.get("values", {})) for q in rf])[:200])
    c_rel, c_abs = thermal_both(r)
    check(f"(5 回目 M1) {label}: 成果物の収束 (相対・絶対) → PASS (旧形式の residual_finite だけ正規化して停止の時の記録と一致)",
          all(c["result"] == "PASS" and c["stop_rule_reeval"].get("matches") is True for c in (c_rel, c_abs)),
          f"{str(c_rel.get('reason'))[:200]}")
    if old_sr is not None:
        c_old = with_stop_rule(old_sr, lambda: WT.thermal_use(r.resolve())["convergence"])
        check(f"(5 回目 M1) {label}: 検出力 — 修正前 ({PRE_M5_COMMIT}) の stop_rule では同じ写しで収束「未評価」(レビューの再現: キーが違う)",
              c_old["result"] == "未評価" and "キーが違う" in str((c_old.get("stop_rule_reeval") or {}).get("problems") or c_old.get("reason")),
              str(c_old.get("reason"))[:200])
    r2 = copy_run(src, f"{src.name}_legacy_v3_input")
    bump_massflux(r2)
    c2_rel, c2_abs = thermal_both(r2)
    check(f"(5 回目 M1) {label}: 流量の CSV を各 +1 kg/s → 成果物の収束は未評価 (相対・絶対; 入力の sha256 が停止の時と違う)",
          all(c["result"] == "未評価" and "boundary_massflux.csv の中身 (sha256) が記録と違う" in (c.get("reason") or "") for c in (c2_rel, c2_abs)),
          str(c2_rel.get("reason"))[:200])
    # 正規化は旧形式に限る: 記録に required が有る (新形式) のに再計算と値が違えば今どおり食い違い
    a = [{"kind": ST.KIND_RESIDUAL_FINITE, "values": {"required": ["rms_ro"], "required_source": "既定 (…)"}}]
    b = [{"kind": ST.KIND_RESIDUAL_FINITE, "values": {"required": ["rms_ro", "rms_roK"]}}]
    check(f"(5 回目 M1) {label}: 正規化は旧形式 (2 つとも無い記録) に限る — 記録に required が有れば今どおり比べる",
          bool(SR.quantities_mismatch(a, b)), str(SR.quantities_mismatch(a, b))[:160])


def test_rejudge_chain(label: str, src: Path, old: dict | None, old_b: dict | None):
    """(#16) 通常の CLI → 成果物の収束 → 受け渡しの照合 (古い形・新しい形・検出力)。old・old_b は修正前 (PRE_M1_COMMIT・PRE_M1B_COMMIT)
    のモジュール {sr, dl, wt (, cs)} (読めなければ None)。"""
    tag = src.name

    def conv_ok(c, form):
        return (c["result"] == "PASS" and c["stop_rule_reeval"].get("matches") is True
                and (c["stop_rule_reeval"].get("compare") or {}).get("form") == form)

    def binding_both(r):
        cwd = os.getcwd()
        try:
            os.chdir(r.parent)
            b_rel = DL._v2_binding(Path(r.name), SR.read_state(r), ST)
        finally:
            os.chdir(cwd)
        return b_rel, DL._v2_binding(r.resolve(), SR.read_state(r), ST)

    def prepare(form, name):
        """写し (form "stop_file" は新しい形にする)。戻り (写し, 専用ファイルの相対パス, 状態の reeval)。"""
        r = copy_run(src, name)
        rel = to_new_form(r) if form == "stop_file" else SR.legacy_stop_verdict_relpath(SR.read_state(r)["reeval"])
        return r, rel, SR.read_state(r)["reeval"]

    for form in ("migrated", "stop_file"):
        what = "古い形 (写しのまま; 専用ファイルなし)" if form == "migrated" else "新しい形 (停止の時の専用ファイルあり)"
        r, rel, ve = prepare(form, f"{tag}_m1_{form}")
        orig = (r / ST.VERDICT_FILE).read_bytes()
        state_bytes = (r / SR.STATE).read_bytes()
        stop_bytes = (r / rel).read_bytes() if form == "stop_file" else None
        before_missing = form == "stop_file" or not (r / rel).exists()
        rc, last, out = rejudge_cli(r)
        if form == "stop_file":
            ok = (r / rel).read_bytes() == stop_bytes and "移行:" not in out
            tail = "、専用ファイルは変わらない"
        else:
            ok = (before_missing and "移行:" in out and (r / rel).read_bytes() == orig and SR._sha256(r / rel) == ve["sha256"]
                  and (r / SR.STATE).read_bytes() == state_bytes)
            tail = f"、上書きの前に停止の時の原本を {rel} へ移した (バイト同一・状態の記録の sha256)、STOP_RULE_STATE.json は変えない"
        check(f"(#16) {label} {what}: 通常の CLI (相対のパス・既定で書く) → SETTLED・終了コード 0、SETTLE_VERDICT.json は停止の時と違うバイト"
              + tail, rc == 0 and last.startswith("VERDICT: SETTLED") and SR._sha256(r / ST.VERDICT_FILE) != ve["sha256"] and ok,
              f"rc {rc} {last}")
        c_rel, c_abs = thermal_both(r)
        check(f"(#16) {label} {what}: CLI の後の成果物の収束 (相対・絶対のパス) → 停止規則 v2 の再評価と一致・PASS (記録の形 {form})",
              conv_ok(c_rel, form) and conv_ok(c_abs, form), f"{str(c_rel.get('reason'))[:200]} / {str(c_abs.get('reason'))[:200]}")
        b_rel, b = binding_both(r)
        check(f"(#16) {label} {what}: CLI の後の受け渡しの照合 (_v2_binding、相対・絶対のパス) が通る (記録の形 {form}、停止の時の原本 {rel})",
              all(x["pass"] and x["items"]["reeval"]["form"] == form and x["items"]["reeval"]["exact"] is True and x["verdict_file"] == rel
                  and "rejudge_match" not in x["items"]["reeval"] for x in (b_rel, b)),
              str((b_rel["problems"] + b["problems"])[:3])[:300])
        (r / rel).write_text((r / rel).read_text().replace('"verdict": "SETTLED"', '"verdict": "SETTLED" ', 1))
        c_rel, _ = thermal_both(r)
        b = DL._v2_binding(r, SR.read_state(r), ST)
        key = "停止の後に書き換えられた" if form == "stop_file" else "停止の時の原本でない"
        check(f"(#16) {label} {what}: 専用ファイル{'' if form == 'stop_file' else ' (移したもの)'}を停止の後に書き換えた → 成果物の収束は未評価・"
              f"受け渡しの照合は失敗 (理由に「{key}」)",
              c_rel["result"] == "未評価" and key in (c_rel.get("reason") or "") and not b["pass"] and any(key in p_ for p_ in b["problems"]),
              f"{str(c_rel.get('reason'))[:200]} / {b['problems'][:1]}")

    # (4 回目 M1) 入力を変える → 通常の CLI → 成果物・受け渡しが拒否される (古い形・新しい形)
    for form in ("migrated", "stop_file"):
        what = "古い形" if form == "migrated" else "新しい形"
        r, rel, ve = prepare(form, f"{tag}_m1b_input_{form}")
        orig = (r / ST.VERDICT_FILE).read_bytes()
        bump_massflux(r)
        rc, last, out = rejudge_cli(r)
        c_rel, c_abs = thermal_both(r)
        b_rel, b = binding_both(r)
        key_c, key_b = "boundary_massflux.csv の中身 (sha256) が記録と違う", "入力が今の run と違う"
        check(f"(#16 続き) {label} {what}: 流量の CSV の流入・流出を各 +1 kg/s → 通常の CLI (終了コード {rc}・{last[:24]}) → 成果物の収束"
              "は未評価 (相対・絶対)、受け渡しの照合は失敗 (相対・絶対) — 入力の sha256 が停止の時と違う"
              + ("; 原本は移した (拒否は入力の照合による)" if form == "migrated" else ""),
              rc == 0 and last.startswith("VERDICT: SETTLED")
              and all(c["result"] == "未評価" and key_c in (c.get("reason") or "") for c in (c_rel, c_abs))
              and all(not x["pass"] and any(key_b in p_ for p_ in x["problems"]) and x["items"]["reeval"]["form"] == form
                      for x in (b_rel, b))
              and (form == "stop_file" or (r / rel).read_bytes() == orig),
              f"{str(c_rel.get('reason'))[:200]} / {str(b['problems'][:1])[:200]}")

    # (4 回目 M1) 古い形で移行の前に SETTLE_VERDICT.json が上書きされていた (修正前の CLI で上書き; 原本なし) → 照合不能で拒否
    if old_b is not None:
        r, rel, ve = prepare("migrated", f"{tag}_m1b_lost")
        rejudge_cli(r, old_b["cs"])
        rc, last, out = rejudge_cli(r)
        c_rel, c_abs = thermal_both(r)
        b_rel, b = binding_both(r)
        key = "原本が失われた"
        check(f"(#16 続き) {label} 古い形: 移行の前に SETTLE_VERDICT.json が上書きされていた (修正前 {PRE_M1B_COMMIT} の CLI) → 今の CLI は"
              f"移行しない、成果物の収束は未評価・受け渡しの照合は失敗 (理由に「{key}」)",
              rc == 0 and "移行:" not in out and not (r / rel).exists()
              and all(c["result"] == "未評価" and key in (c.get("reason") or "") for c in (c_rel, c_abs))
              and all(not x["pass"] and any(key in p_ for p_ in x["problems"]) for x in (b_rel, b)),
              f"{str(c_rel.get('reason'))[:200]} / {str(b['problems'][:1])[:200]}")

    if old is not None:
        r = copy_run(src, f"{tag}_m1_detect")
        sr_o, dl_o, wt_o = old["sr"], old["dl"], old["wt"]
        pre = with_stop_rule(sr_o, lambda: (thermal_both(r, wt_o)[0], dl_o._v2_binding(r, sr_o.read_state(r), ST)))
        rc, last, _ = rejudge_cli(r)
        post = with_stop_rule(sr_o, lambda: (thermal_both(r, wt_o)[0], dl_o._v2_binding(r, sr_o.read_state(r), ST)))
        check(f"(#16) {label}: 検出力 — 修正前は CLI の前は成果物の収束 PASS・照合が通り、通常の CLI (終了コード {rc}) の後は成果物の収束"
              "「未評価」・受け渡しの照合が失敗 (レビューの再現)",
              rc == 0 and pre[0]["result"] == "PASS" and pre[1]["pass"] and post[0]["result"] == "未評価" and not post[1]["pass"],
              f"前 {pre[0]['result']}・{pre[1]['pass']} / 後 {post[0]['result']}・{str(post[1]['problems'][:1])[:200]}")
    if old_b is not None:
        sr_b, dl_b, wt_b, cs_b = old_b["sr"], old_b["dl"], old_b["wt"], old_b["cs"]
        res = {}
        for form in ("migrated", "stop_file"):
            r, _, _ = prepare(form, f"{tag}_m1b_detect_{form}")
            bump_massflux(r)
            rc, last, _ = rejudge_cli(r, cs_b)
            res[form] = (rc, with_stop_rule(sr_b, lambda: (thermal_both(r, wt_b)[0], dl_b._v2_binding(r, sr_b.read_state(r), ST))))
        (rc_l, (c_l, b_l)), (rc_n, (c_n, b_n)) = res["migrated"], res["stop_file"]
        check(f"(#16 続き) {label}: 検出力 — 修正前 ({PRE_M1B_COMMIT}) は流量の CSV を変えて通常の CLI (その版) を回すと、古い形で成果物の収束 "
              "PASS (記録の形 overwritten)・受け渡しの照合が通り、新しい形でも受け渡しの照合が通った (レビューの再現; 今は両方とも拒否)",
              rc_l == 0 and rc_n == 0 and c_l["result"] == "PASS" and (c_l["stop_rule_reeval"].get("compare") or {}).get("form") == "overwritten"
              and b_l["pass"] and b_n["pass"],
              f"古い形 {c_l['result']}・{b_l['pass']} / 新しい形 {c_n['result']}・{b_n['pass']}")


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
    try:
        old_m1 = {"sr": load_at(PRE_M1_COMMIT, "design/forge_design/evaluate/stop_rule.py", "forge_design.evaluate._stop_rule_pre_m1",
                                package="forge_design.evaluate"),
                  "dl": load_at(PRE_M1_COMMIT, "design/forge_design/feedback/deltastar_loop.py", "forge_design.feedback._deltastar_loop_pre_m1",
                                package="forge_design.feedback"),
                  "wt": load_at(PRE_M1_COMMIT, "design/forge_design/export/wall_tables.py", "forge_design.export._wall_tables_pre_m1",
                                package="forge_design.export")}
        for k_ in ("REPO", "TOOLS"):          # 道具 (check_settled.py) は今の作業ツリーのもの (比べるのは照合の中身だけ)
            setattr(old_m1["wt"], k_, getattr(WT, k_))
    except Exception as e:  # noqa: BLE001
        old_m1 = None
        check(f"修正前 ({PRE_M1_COMMIT}) のモジュールを git から読めた (読めなければ検出力は判定不能 — 合格にしない)", False, str(e))
    try:
        old_m1b = {"sr": load_at(PRE_M1B_COMMIT, "design/forge_design/evaluate/stop_rule.py", "forge_design.evaluate._stop_rule_pre_m1b",
                                 package="forge_design.evaluate"),
                   "dl": load_at(PRE_M1B_COMMIT, "design/forge_design/feedback/deltastar_loop.py",
                                 "forge_design.feedback._deltastar_loop_pre_m1b", package="forge_design.feedback"),
                   "wt": load_at(PRE_M1B_COMMIT, "design/forge_design/export/wall_tables.py", "forge_design.export._wall_tables_pre_m1b",
                                 package="forge_design.export"),
                   "cs": load_at(PRE_M1B_COMMIT, "solver_density_cuda/tools/check_settled.py", "check_settled_pre_m1b")}
        for k_ in ("REPO", "TOOLS"):          # 成果物の道具 (check_settled.py --bind --no-write) は今の作業ツリーのもの (書かないので移行は無関係)
            setattr(old_m1b["wt"], k_, getattr(WT, k_))
    except Exception as e:  # noqa: BLE001
        old_m1b = None
        check(f"修正前 ({PRE_M1B_COMMIT}) のモジュールを git から読めた (読めなければ検出力は判定不能 — 合格にしない)", False, str(e))
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
        # ---------------------------------------------------------------- #16
        test_rejudge_chain(label, src, old_m1, old_m1b)
    # ---------------------------------------------------------------- 5 回目 M1 (試走 A の旧形式の v3 の記録)
    try:
        old_sr5 = load_at(PRE_M5_COMMIT, "design/forge_design/evaluate/stop_rule.py", "forge_design.evaluate._stop_rule_pre_m5",
                          package="forge_design.evaluate")
    except Exception as e:  # noqa: BLE001
        old_sr5 = None
        check(f"修正前 ({PRE_M5_COMMIT}) の stop_rule を git から読めた (読めなければ検出力は判定不能 — 合格にしない)", False, str(e))
    for label, src in TRIAL.items():
        if not (src / ST.CRITERIA_FILE).is_file() or not (src / SR.STATE).is_file():
            SKIPPED.append(f"{label} ({src.relative_to(ROOT)} が無い)")
            continue
        print(f"== {label} ({src.relative_to(ROOT)})", flush=True)
        test_legacy_v3(label, src, old_sr5)
    if SKIPPED:
        print("飛ばした項目: " + "; ".join(SKIPPED), flush=True)
    print(f"\n{NCHK - FAIL}/{NCHK} ok — {'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'} (一時ディレクトリ {TMP})", flush=True)
    if FAIL == 0:
        shutil.rmtree(TMP, ignore_errors=True)
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
