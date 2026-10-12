#!/usr/bin/env python3
"""plan tooling-nozzle-sern-chain §5.1 R7b「④ の改訂 (2026-10-13)」の ④-0 を実行する台本 (判断: codex diagnose
notes/reviews/2026-10-13-sern-r7b-lsw-decision-diagnose.md — 全件採用)。forge は回さない (0 step)。

既存の判定区間 7 本 (g3 の L_sw 0.8 / 1.0 × m10_on・m6_on・m4_off と、g4 の L_sw 0.8・m10_on) を、元データ・元設定に結び付けて
現行の受理検査 (runner_sern3d.collect、受理検査の版 sern_gates.ACCEPTANCE_VERSION) で判定し直す。元記録は変えない:
再判定の出力は run の中の新しいサブディレクトリ (既定 REJUDGE_2026-10-13/、既にあれば失敗) に置く。受理検査だけを
旧版 (run に保存された metrics.json) / 現行版で変える A/B を残す。1 本ずつ処理する (h5 は読むだけ)。

run ごとの手順 (どれかが合わなければ「判定不能」として理由を記録し、そこで止める。判定はしない):
  (0) 問題 YAML: その run を作った te_wake の検証用 YAML を基本とし、判定区間の段の判定 (#1〜#3) を有効にするキー
      (evaluate.judgment_steps・judgment_extend_once・floor_events・floor_events_window_steps・outStepInterval) だけを生産 YAML の
      値にした YAML (problem_rejudge.yaml) を作る。使った値・両 YAML の sha256 と最後のコミット・他に違うキーを記録する
  (1) 来歴: run の prepare_info.json の problem (prepare に使った YAML) の名前と作動点が対応表と一致すること
  (2) 入力: prepare_info.json・solverConfig.yaml・bcondConfig.yaml・格子 (meshFileName)・residual_history.csv・最終場 res_N.h5
      (N = 判定区間の段の長さ)・壁の部品の出力 res_<部品>_<タグ>_<step>.h5 (力の部品の要求 × 予定 step 列の全部)・
      floor_events.csv・TE_MONITOR.csv・cowl_contour.csv・run_case_stdout.log の forge の終了コード。欠けていれば列挙する
  (3) 設定: 現行 runner の prepare を同じ YAML・作動点で prepare_current/ に回し (格子生成・変換・品質。段階起動はしない)、
      run の solverConfig.yaml を現行の判定区間の config (本段の config の nStepOuter を N にしたもの = runner_sern.judgment_config) と、
      run の solverConfig_main.yaml (あれば) を現行の本段の config と比べる。許す差は nStepOuter・出力系だけ
      (diag/judgment_r2c.py の ALLOWED_DIFF と同じ規則)。bcondConfig.yaml・species_meta.yaml・species_db_external.yaml は一致すること
      (同じく SAME_FILES)。格子署名・双対幾何のハッシュは、run の格子 (いま再計算) と prepare_info.json の記録 (あれば) が
      現行の prepare と一致すること。collect が読む prepare_info.json の量 (H・半幅・理想推力・入口と外気の状態・離散化・作動点・
      L_sw・te_wake の L_b・壁の部品の面の数) も現行の prepare と一致すること (台本の追加の照合)。格子・中間ファイルは照合の後に消す
  (4) 温度の監視: run の TE_MONITOR.csv の列が、現行の te_monitor.py を最終場に回し直した出力の列と一致し、最終 step の行の値も一致すること
  (5) 判定: runner_sern3d.collect(problem_rejudge.yaml, run, out_dir=<再判定のディレクトリ>) → metrics.json・force_history.csv
      (#2 の温度は run の既存の TE_MONITOR.csv、forge の終了コードは run_case_stdout.log から)
  (6) A/B: 保存された旧 metrics.json の gates と現行の gates を項目ごとに並べる (AB_GATES.md)
  最後に run の既存ファイル (再判定のサブディレクトリを除く) の大きさ・mtime・sha256 が処理の前後で同じことを確かめる。
7 本のまとめ (<case>/R7B_<tag>.md・.json、既にあれば失敗) は結果 A (7 本すべて必要条件を満たす) / 結果 B (欠落・不合格あり) を
機械的に出し、現行の verdict が PASS の行から g3 の設計差 Δᵢ・ΔJ と m10_on の格子感度の参考値を表示する (解釈はしない)。

使い方 (AWS-B、~/forge-r8 = feature/sern-design の checkout。solver_density_cuda/build/convertGmshToForge が要る。長いので nohup で):
  nohup python3 case/46.sern_design/diag/rejudge_r7b.py > case/46.sern_design/rejudge_r7b.log 2>&1 &
  python3 case/46.sern_design/diag/rejudge_r7b.py --only run_1084_tewake_g3_L10_m10      # 1 本だけ (まとめは記録のある run から)
  python3 case/46.sern_design/diag/rejudge_r7b.py --summary-only --summary-out <path>      # 既存の再判定の記録からまとめだけ
終了コード: 0 = 対象の全 run の記録を書いた (判定不能を含む)、2 = 再判定の出力先が既にある run があった、
1 = その他の失敗 (台本の例外・空き容量の不足・元ファイルの変化・まとめの出力先が既にある)。
"""
from __future__ import annotations

import argparse
import copy
import csv
import dataclasses
import difflib
import hashlib
import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DESIGN = ROOT / "design"
TOOLS = ROOT / "solver_density_cuda" / "tools"
for p_ in (str(DESIGN), str(TOOLS)):
    if p_ not in sys.path:
        sys.path.insert(0, p_)

import yaml  # noqa: E402

from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402
from forge_design.metrics import sern_gates as G  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

C46 = ROOT / "case" / "46.sern_design"
TAG = "REJUDGE_2026-10-13"
VERSION = "rejudge-r7b-v1"
# 判定区間の段の判定 (#1〜#3) を有効にするキー (生産 YAML の値で与える)
JUDGMENT_KEYS = ("judgment_steps", "judgment_extend_once", "floor_events", "floor_events_window_steps", "outStepInterval")
# plan §5.1 R7b ④ の加重 (ΔJ = 0.5Δ₆ + 0.3Δ₁₀ + 0.2Δ₄)、q = C_T_with_shear、u = max(a, d, 1e-6)
WEIGHTS = {"m6_on": 0.5, "m10_on": 0.3, "m4_off": 0.2}
OBJ = "C_T_with_shear"
QUANTS = ("C_T", "C_T_with_shear", "C_L", "C_M")
U_FLOOR = 1.0e-6
INFO_JSON, AB_MD, SUMMARY_TXT, CONFIG_DIFF_TXT = "REJUDGE_INFO.json", "AB_GATES.md", "REJUDGE_SUMMARY.txt", "CONFIG_DIFF.txt"
PREP_DIR, MERGED_YAML, TE_MON_CUR = "prepare_current", "problem_rejudge.yaml", "TE_MONITOR_current_final.csv"
_MISSING = "(無い)"
Y08, Y10 = "problem_3d_prod_3op_wallres_lswx08", "problem_3d_prod_3op_wallres_lswx10"


@dataclasses.dataclass(frozen=True)
class Spec:
    run: str
    L_sw: float
    op: str
    grid: str
    te_yaml: str          # その run を作った YAML (te_wake の検証用)
    prod_yaml: str        # 判定区間のキーを取る生産 YAML
    note: str = ""        # 来歴 (case README の run 一覧・diag/tewake_*.sh)


# 来歴は case/46.sern_design/README.md の run 一覧と diag/tewake_*.sh (prepare の YAML は各 run の prepare_info.json の problem で照合する)
RUNS = (
    Spec("run_1084_tewake_g3_L10_m10", 0.8, "m10_on", "g3", f"{Y08}_tewake_g3_L10.yaml", f"{Y08}.yaml",
         "prepare (YAML lswx08_tewake_g3_L10) → run_1055 の最終場から restart_field_deformed、20000 step (diag/tewake_g3_local.sh)"),
    Spec("run_1091_tewake_m6_L10_cont40k", 0.8, "m6_on", "g3", f"{Y08}_tewake_g3_L10.yaml", f"{Y08}.yaml",
         "run_1088 (prepare → run_1050 の最終場から restart_field_deformed、diag/tewake_c6.sh) の +20000 (diag/tewake_ext.sh)"),
    Spec("run_1090_tewake_m4_L10", 0.8, "m4_off", "g3", f"{Y08}_tewake_g3_L10.yaml", f"{Y08}.yaml",
         "prepare → run_1058 の最終場から restart_field_deformed、20000 step (diag/tewake_c6.sh)"),
    Spec("run_1096_tewake_lsw10_m10_cont20k", 1.0, "m10_on", "g3", f"{Y10}_tewake_g3_L10.yaml", f"{Y10}.yaml",
         "run_1093 (段階起動) の本段の最終場から 20000 step (diag/tewake_lsw10.sh → tewake_ext.sh)"),
    Spec("run_1099_tewake_lsw10_m6_cont40k", 1.0, "m6_on", "g3", f"{Y10}_tewake_g3_L10.yaml", f"{Y10}.yaml",
         "run_1097 (run_1094 の段階起動の本段の最終場から 20000 step) の +20000 (diag/tewake_ext.sh)"),
    Spec("run_1098_tewake_lsw10_m4_cont20k", 1.0, "m4_off", "g3", f"{Y10}_tewake_g3_L10.yaml", f"{Y10}.yaml",
         "run_1095 (段階起動) の本段の最終場から 20000 step (diag/tewake_lsw10.sh → tewake_ext.sh)"),
    Spec("run_1086_tewake_g4_L10_cont40k", 0.8, "m10_on", "g4", f"{Y08}_tewake_g4_L10.yaml", f"{Y08}.yaml",
         "run_1085 (run_1082 の段階起動の本段の最終場から 20000 step、diag/tewake_g4_local.sh) の +20000 (diag/tewake_g4_ext.sh)"),
)

# collect が読む prepare_info.json の量 (力の正規化・帳簿・監視の形状)。現行の prepare と一致すること
INFO_KEYS = (("H_m",), ("half_W_m",), ("half_W_vehicle_m",), ("F_ideal_N_per_m",), ("discretization",),
             ("states", "exhaust", "ro"), ("states", "exhaust", "u"), ("states", "exhaust", "P"), ("states", "ext", "P"),
             ("operating_point", "name"), ("mesh", "L_sw"), ("mesh", "te_wake_blend_H"))


class Stop(Exception):
    """判定不能 (入力の欠け・照合の不一致)。kind = yaml / provenance / inputs / prepare / config / mesh / prepare_info / te_monitor。"""

    def __init__(self, kind: str, reasons):
        super().__init__("; ".join(reasons))
        self.kind = kind
        self.reasons = list(reasons)


def _r2c():
    """diag/judgment_r2c.py (R2c の台本) の照合規則 (ALLOWED_DIFF・SAME_FILES) を同じものとして使う。"""
    spec = importlib.util.spec_from_file_location("_judgment_r2c", HERE / "judgment_r2c.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R2C = _r2c()


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def _git(*args):
    r = subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else None


def git_head() -> dict:
    st = _git("status", "--porcelain", "--untracked-files=no")
    return {"head": _git("rev-parse", "HEAD"), "dirty_tracked": None if st is None else len([x for x in st.splitlines() if x.strip()])}


def yaml_record(path: Path) -> dict:
    return {"path": str(path), "sha256": sha256(path), "git_last_commit": _git("log", "-1", "--format=%H", "--", str(path)),
            "git_modified": bool(_git("status", "--porcelain", "--", str(path)))}


def snapshot(rd: Path, exclude: Path) -> dict:
    """run の既存ファイル (exclude の下を除く) の {相対パス: {size, mtime_ns, sha256}}。読むだけ。"""
    out = {}
    ex = exclude.resolve()
    for dp, dns, fns in os.walk(rd):
        dpp = Path(dp)
        dns[:] = sorted(d for d in dns if (dpp / d).resolve() != ex)
        for fn in sorted(fns):
            f = dpp / fn
            if f.is_symlink() or not f.is_file():
                out[str(f.relative_to(rd))] = {"symlink_or_special": True, "target": os.readlink(f) if f.is_symlink() else None}
                continue
            st = f.stat()
            out[str(f.relative_to(rd))] = {"size": st.st_size, "mtime_ns": st.st_mtime_ns, "sha256": sha256(f)}
    return out


def snapshot_changes(a: dict, b: dict) -> list:
    ch = [f"消えた: {k}" for k in sorted(set(a) - set(b))] + [f"増えた: {k}" for k in sorted(set(b) - set(a))]
    ch += [f"変わった: {k} ({', '.join(x for x in ('size', 'mtime_ns', 'sha256') if a[k].get(x) != b[k].get(x))})"
           for k in sorted(set(a) & set(b)) if a[k] != b[k]]
    return ch


def _get(d, path):
    for k in path:
        if not isinstance(d, dict) or k not in d:
            return _MISSING
        d = d[k]
    return d


def _same(x, y) -> bool:
    if isinstance(x, bool) or isinstance(y, bool):
        return x == y
    if isinstance(x, (int, float)) and isinstance(y, (int, float)):
        return math.isclose(float(x), float(y), rel_tol=1e-12, abs_tol=0.0) or (math.isnan(float(x)) and math.isnan(float(y)))
    return x == y


def _flat(d, pre="") -> dict:
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(_flat(v, f"{pre}.{k}" if pre else str(k)))
    else:
        out[pre] = d
    return out


def _yaml_path(name: str, case_dir: Path) -> Path:
    p = Path(name)
    return p if p.is_absolute() else case_dir / p


# ----------------------------------------------------------------------------------------------- (0) 問題 YAML
def merged_problem(spec: Spec, case_dir: Path, out: Path) -> tuple:
    te, pr = _yaml_path(spec.te_yaml, case_dir), _yaml_path(spec.prod_yaml, case_dir)
    lack = [str(p) for p in (te, pr) if not p.is_file()]
    if lack:
        raise Stop("yaml", [f"問題 YAML が無い: {', '.join(lack)}"])
    yt, yp = yaml.safe_load(te.read_text()), yaml.safe_load(pr.read_text())
    evt, evp = (yt.get("evaluate") or {}), (yp.get("evaluate") or {})
    lk = [k for k in JUDGMENT_KEYS if k not in evp]
    if lk:
        raise Stop("yaml", [f"生産 YAML {pr.name} の evaluate に {', '.join(lk)} が無い"])
    L_sw = (yt.get("mesh3d") or {}).get("L_sw")
    if L_sw is None or not _same(float(L_sw), float(spec.L_sw)):
        raise Stop("yaml", [f"te_wake の YAML {te.name} の mesh3d.L_sw {L_sw!r} が対応表の {spec.L_sw} と違う"])
    names = [o.get("name") for o in ((yt.get("spec") or {}).get("operating_points") or [])]
    if spec.op not in names:
        raise Stop("yaml", [f"te_wake の YAML {te.name} に作動点 {spec.op} が無い ({names})"])
    merged = copy.deepcopy(yt)
    merged.setdefault("evaluate", {})
    for k in JUDGMENT_KEYS:
        merged["evaluate"][k] = copy.deepcopy(evp[k])
    # gas.species_db の相対パスは YAML の所在基準 (probdef.Problem.species_db) なので、置き場所を変える前に絶対パスにする
    sdb = (merged.get("gas") or {}).get("species_db")
    if sdb and not Path(str(sdb)).is_absolute():
        merged["gas"]["species_db"] = str((te.parent / str(sdb)).resolve())
    path = out / MERGED_YAML
    head = (f"# {VERSION}: plan tooling-nozzle-sern-chain §5.1 R7b ④-0 の再判定用 ({now()})。\n"
            f"# 基本 = {te} (run {spec.run} を作った YAML)、evaluate の {', '.join(JUDGMENT_KEYS)} だけ {pr} の値。\n")
    path.write_text(head + yaml.safe_dump(merged, allow_unicode=True, sort_keys=False))
    ft, fp = _flat(yt), _flat(yp)
    skip = {"name"} | {f"evaluate.{k}" for k in JUDGMENT_KEYS}
    other = [{"key": k, "te_wake": ft.get(k, _MISSING), "production": fp.get(k, _MISSING)}
             for k in sorted(set(ft) | set(fp)) if k not in skip and not _same(ft.get(k, _MISSING), fp.get(k, _MISSING))]
    rec = {"te_wake": yaml_record(te), "production": yaml_record(pr), "merged": {"path": str(path), "sha256": sha256(path)},
           "keys_from_production": {k: {"production": evp[k], "te_wake": evt.get(k, _MISSING), "used": evp[k]} for k in JUDGMENT_KEYS},
           "species_db_abs": merged.get("gas", {}).get("species_db") if sdb else None,
           "other_differences_te_vs_production": other}
    return path, rec


# ----------------------------------------------------------------------------------------------- (1)(2) 来歴・入力
def check_provenance(spec: Spec, rd: Path) -> dict:
    p = rd / "prepare_info.json"
    if not p.is_file():
        return {"checked": False, "why": "prepare_info.json が無い (入力の欠けとして列挙)"}
    info = json.loads(p.read_text())
    used = info.get("problem")
    op = (info.get("operating_point") or {}).get("name")
    rec = {"checked": True, "prepare_problem": used, "prepare_op": op, "expected_yaml": Path(spec.te_yaml).name, "expected_op": spec.op}
    bad = []
    if not used or Path(str(used)).name != Path(spec.te_yaml).name:
        bad.append(f"run の prepare_info.json の problem {used!r} が対応表の YAML {Path(spec.te_yaml).name} と違う")
    if op != spec.op:
        bad.append(f"run の prepare_info.json の作動点 {op!r} が対応表の {spec.op} と違う")
    if bad:
        raise Stop("provenance", bad)
    return rec


def _mesh_name(rd: Path):
    try:
        m = re.search(r'meshFileName:\s*"?([^",}\s]+)"?', (rd / "solverConfig.yaml").read_text())
    except OSError:
        return None
    return m.group(1) if m else None


def _fmt_steps(xs) -> str:
    xs = list(xs)
    return ", ".join(str(x) for x in xs[:4]) + (" …" if len(xs) > 4 else "")


def check_inputs(rd: Path, jp: dict) -> dict:
    """必要な入力の有無 (欠けを列挙)。jp = runner_sern.judgment_plan (予定 step 列)。"""
    N = int(jp["planned_steps"][-1])
    missing, notes = [], []
    for name in ("prepare_info.json", "solverConfig.yaml", "bcondConfig.yaml", "residual_history.csv", "floor_events.csv",
                 "TE_MONITOR.csv", "cowl_contour.csv", "run_case_stdout.log"):
        if not (rd / name).is_file():
            missing.append(name)
    mesh = _mesh_name(rd)
    if (rd / "solverConfig.yaml").is_file():
        if mesh is None:
            missing.append("solverConfig.yaml の meshFileName")
        elif not (rd / mesh).is_file():
            missing.append(f"格子 {mesh} (meshFileName)")
    if not (rd / f"res_{N}.h5").is_file():
        snaps = [f.name for _, f in R2._volume_snapshots(rd)]
        missing.append(f"最終場 res_{N}.h5 (判定区間の段の終了 step。ある全場: {snaps or 'なし'})")
    rc = None
    if (rd / "run_case_stdout.log").is_file():
        rc = G.forge_rc_from_log(rd)
        if rc is None:
            missing.append("run_case_stdout.log の forge の終了コード ([run_case] forge exit=)")
    parts = {}
    if (rd / "prepare_info.json").is_file() and (rd / "bcondConfig.yaml").is_file():
        fp = R3.force_parts_plan(rd)
        parts = {"required": fp["required"], "friction": fp["friction"], "faces_source": fp["faces_source"], "reasons": fp["reasons"]}
        if fp["reasons"]:
            notes += [f"力の部品の要求: {r} (collect が判定不能にする)" for r in fp["reasons"]]
        for n in fp["required"]:
            pid = R3._FORCE_PARTS[n][0]
            lack = [s for s in jp["planned_steps"] if not (rd / f"res_{n}_{pid}_{int(s)}.h5").is_file()]
            if lack:
                missing.append(f"壁の部品 {n} の出力 res_{n}_{pid}_<step>.h5 が {len(lack)}/{len(jp['planned_steps'])} 標本で無い "
                               f"(step {_fmt_steps(lack)})")
    return {"missing": missing, "notes": notes, "final_step": N, "mesh_file": mesh, "forge_rc_from_log": rc, "force_parts": parts,
            "planned": {"first": int(jp["planned_steps"][0]), "last": N, "n": len(jp["planned_steps"]),
                        "window_steps": jp["window_steps"], "out_interval": jp["out_interval"]}}


# ----------------------------------------------------------------------------------------------- (3) 設定・格子の照合
def run_prepare(merged: Path, op: str, prep: Path) -> dict:
    """現行 runner の prepare (runner_sern3d --prepare-only、別プロセス = 格子のメモリを run ごとに返す)。"""
    t0 = time.time()
    r = subprocess.run([sys.executable, "-m", "forge_design.evaluate.runner_sern3d", str(merged), str(prep), "--op", op, "--prepare-only"],
                       cwd=DESIGN, capture_output=True, text=True)
    (prep.parent / f"{PREP_DIR}.log").write_text(f"$ runner_sern3d {merged} {prep} --op {op} --prepare-only\n{r.stdout}{r.stderr}"
                                                 f"EXIT: {r.returncode}\n")
    if r.returncode != 0 or not (prep / "prepare_info.json").is_file():
        tail = (r.stderr or r.stdout or "").strip().splitlines()[-4:]
        raise Stop("prepare", [f"現行 runner の prepare が失敗 (rc {r.returncode}): " + " / ".join(tail)])
    info = json.loads((prep / "prepare_info.json").read_text())
    info["_elapsed_s"] = round(time.time() - t0, 1)
    return info


def _diff(old_text: str, new_text: str, old_label: str, new_label: str) -> tuple:
    """(unified diff, 許す差の行, 許さない差の行)。規則は judgment_r2c.ALLOWED_DIFF (nStepOuter・出力系)。"""
    ud = list(difflib.unified_diff(old_text.splitlines(), new_text.splitlines(), old_label, new_label, lineterm=""))
    changed = [ln[1:] for ln in ud if ln[:1] in "+-" and not ln.startswith(("+++", "---"))]
    return "\n".join(ud), [ln for ln in changed if R2C.ALLOWED_DIFF.match(ln)], [ln for ln in changed if not R2C.ALLOWED_DIFF.match(ln)]


def compare_configs(rd: Path, prep: Path, jo: dict, out: Path) -> dict:
    cur_main = (prep / "solverConfig_main.yaml").read_text()
    try:
        cur_judg = R2.judgment_config(cur_main, jo)
    except (RuntimeError, ValueError) as e:
        raise Stop("config", [f"現行の判定区間の config を作れない ({type(e).__name__}: {e})"])
    rec, lines = {"compared": {}, "files": {}}, []
    pairs = [("solverConfig.yaml", cur_judg, f"現行の判定区間の config (runner_sern.judgment_config({PREP_DIR}/solverConfig_main.yaml))")]
    if (rd / "solverConfig_main.yaml").is_file():
        pairs.append(("solverConfig_main.yaml", cur_main, f"現行の本段の config ({PREP_DIR}/solverConfig_main.yaml)"))
    else:
        rec["compared"]["solverConfig_main.yaml"] = "run に無い (比べない。run の判定区間は solverConfig.yaml が読んだ)"
    for name, cur, what in pairs:
        ud, ok, bad = _diff((rd / name).read_text(), cur, f"run/{name}", f"current/{name}")
        rec["compared"][name] = {"against": what, "allowed": ok, "disallowed": bad}
        lines += [f"# run/{name} と {what}", f"# 許す差 (nStepOuter・出力系) {len(ok)} 行、許さない差 {len(bad)} 行", ud or "(差なし)", ""]
    for name in R2C.SAME_FILES:
        a, b = rd / name, prep / name
        if not (a.exists() or b.exists()):
            continue
        ta, tb = (a.read_text() if a.exists() else None), (b.read_text() if b.exists() else None)
        same = ta == tb
        rec["files"][name] = same if (a.exists() and b.exists()) else ("run にだけある" if a.exists() else "現行の prepare にだけある")
        d = "" if same else ("\n".join(difflib.unified_diff((ta or "").splitlines(), (tb or "").splitlines(), f"run/{name}",
                                                            f"current/{name}", lineterm="")) or "片方にだけある")
        lines += [f"# {name}: {'差なし' if same else '差あり'}"] + ([d] if d else []) + [""]
    (out / CONFIG_DIFF_TXT).write_text("\n".join(lines) + "\n")
    bad = [f"{n}: {ln}" for n, v in rec["compared"].items() if isinstance(v, dict) for ln in v["disallowed"]]
    bad += [f"{n} が違う ({'差あり' if v is False else v})" for n, v in rec["files"].items() if v is not True]
    rec["disallowed"] = bad
    return rec


def compare_mesh(rd: Path, new_info: dict) -> dict:
    mp = R2.mesh_provenance(rd)       # run の格子 (meshFileName) をいま再計算 + prepare の記録
    new = new_info.get("mesh_provenance") or {}
    now_, pre = mp.get("now") or {}, mp.get("prepare") or {}
    rec = {"run_now": now_, "run_record": mp.get("prepare"), "current": new, "items": {}}
    bad = []
    if "error" in now_:
        bad.append(f"run の格子の来歴を再計算できない ({now_['error']})")
    for k in ("mesh_signature", "dual_hash", "discretization"):
        c, n_, r_ = new.get(k), now_.get(k), (pre.get(k) if isinstance(pre, dict) else None)
        rec["items"][k] = {"current": c, "run_now": n_, "run_record": r_}
        if c is None:
            bad.append(f"現行の prepare の {k} が無い")
            continue
        if "error" not in now_ and n_ != c:
            bad.append(f"{k}: run の格子 (いま再計算) {n_} ≠ 現行の prepare {c}")
        if r_ is not None and r_ != c:
            bad.append(f"{k}: run の prepare_info.json の記録 {r_} ≠ 現行の prepare {c}")
        if r_ is None:
            rec["items"][k]["note"] = "run の prepare_info.json に記録が無い (いま再計算した値だけで照合)"
    rec["disallowed"] = bad
    return rec


def compare_info(rd: Path, new_info: dict) -> dict:
    old = json.loads((rd / "prepare_info.json").read_text())
    rec, bad = {"items": {}}, []
    for path in INFO_KEYS:
        a, b = _get(old, path), _get(new_info, path)
        key = ".".join(path)
        rec["items"][key] = {"run": a, "current": b, "same": _same(a, b)}
        if not _same(a, b):
            bad.append(f"prepare_info.{key}: run {a!r} ≠ 現行の prepare {b!r}")
    wa, wb = old.get(R3.WALL_FACES_KEY), new_info.get(R3.WALL_FACES_KEY)
    if isinstance(wa, dict) and isinstance(wb, dict):
        rec["items"][R3.WALL_FACES_KEY] = {"run": wa, "current": wb, "same": wa == wb}
        if wa != wb:
            bad.append(f"prepare_info.{R3.WALL_FACES_KEY} (壁の部品の面の数) が違う: run {wa} / 現行 {wb}")
    else:
        rec["items"][R3.WALL_FACES_KEY] = {"run": wa, "current": wb, "same": None,
                                           "note": "片方に記録が無い (collect はメッシャの記録 prepare_info.mesh から面の有無を決める)"}
    em_old = (old.get("eval_method") or {}).get(R2.EVAL_METHOD_ID)
    em_new = (new_info.get("eval_method") or {}).get(R2.EVAL_METHOD_ID)
    rec["eval_method_id"] = {"run": em_old, "current": em_new, "note": "記録だけ (照合は格子署名・双対・config で)"}
    rec["disallowed"] = bad
    return rec


def _cleanup_prepare(prep: Path) -> list:
    """現行の prepare の格子・中間ファイル (再生成できる。容量のため) を消し、消した名前と大きさを返す。"""
    gone = []
    for f in sorted(prep.iterdir()):
        if f.is_file() and f.suffix in (".h5", ".msh", ".xmf"):
            gone.append({"name": f.name, "bytes": f.stat().st_size})
            f.unlink()
    return gone


# ----------------------------------------------------------------------------------------------- (4) 温度の監視
def _read_csv(path: Path) -> tuple:
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    return (rows[0] if rows else []), [r for r in rows[1:] if r]


def _cell_same(a: str, b: str) -> bool:
    try:
        return _same(float(a), float(b))
    except ValueError:
        return a == b


def check_te_monitor(rd: Path, N: int, out: Path) -> dict:
    cur = out / TE_MON_CUR
    lb = R3._te_monitor_lb(rd)
    cmd = [sys.executable, str(R3.TE_MONITOR_TOOL), str(rd), str(rd / f"res_{N}.h5"), "--csv", str(cur), "--lb-h", repr(lb)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    rec = {"tool": str(R3.TE_MONITOR_TOOL), "tool_sha256": sha256(R3.TE_MONITOR_TOOL), "lb_h": lb, "rc": r.returncode,
           "stdout": (r.stdout or "").strip()[-400:], "stderr": (r.stderr or "").strip()[-400:]}
    if r.returncode != 0 or not cur.is_file():
        raise Stop("te_monitor", [f"現行の te_monitor.py を最終場 res_{N}.h5 に回せない (rc {r.returncode}): {rec['stderr'][-200:]}"])
    hc, rc_rows = _read_csv(cur)
    ho, ro_rows = _read_csv(rd / R2.MONITOR_CSV)
    bad = []
    rec["columns"] = {"current": len(hc), "run": len(ho), "same": hc == ho}
    if hc != ho:
        lack, extra = [c for c in hc if c not in ho], [c for c in ho if c not in hc]
        bad.append(f"run の TE_MONITOR.csv の列が現行の te_monitor.py の出力と違う (無い列 {lack[:6]}・余分な列 {extra[:6]}"
                   + ("・順序が違う" if not lack and not extra else "") + ")")
    ragged = [i + 2 for i, row in enumerate(ro_rows) if len(row) != len(ho)]
    if ragged:
        bad.append(f"run の TE_MONITOR.csv に見出しと欄数の違う行が {len(ragged)} 行ある (行 {_fmt_steps(ragged)})")
    rec["n_rows_run"] = len(ro_rows)
    if hc == ho and not ragged and rc_rows:
        si = hc.index("step")
        last = [row for row in ro_rows if _cell_same(row[si], str(N))]
        if len(last) != 1:
            bad.append(f"run の TE_MONITOR.csv に最終 step {N} の行が {len(last)} 行 (1 行であること)")
        else:
            diff = [c for c, a, b in zip(hc, last[0], rc_rows[0]) if not _cell_same(a, b)]
            rec["final_row_same"] = not diff
            if diff:
                bad.append(f"最終 step {N} の行が現行の te_monitor.py の値と違う (列 {diff[:6]})")
    if bad:
        raise Stop("te_monitor", bad)
    return rec


# ----------------------------------------------------------------------------------------------- (6) A/B
def _ok(x) -> str:
    return "—" if x is None else ("OK" if x else "NG")


def gates_items(g: dict | None) -> list:
    """gates を項目ごとの (名前, 表示) にする (旧 metrics.json の gates にも使う。無い欄は "(なし)")。"""
    if not isinstance(g, dict):
        return [("gates", "(なし)")]
    fld, res, st = g.get("field") or {}, g.get("residual") or {}, g.get("steadiness") or {}
    fl, fe, rs = g.get("floors") or {}, g.get("floor_events") or {}, g.get("residual_scale") or {}
    jg = g.get("judgment")
    ev = G.residual_evidence_summary(res)
    ser = st.get("series") or {}
    fw, tp = ((jg or {}).get("force_window") or {}), ((jg or {}).get("temperature") or {})
    evs = "(なし)" if ev is None else (f"input_ok {ev.get('input_ok')}・欠けた列 {ev.get('missing_columns')}・"
                                      f"時系列 {(ev.get('series') or {}).get('n_outer_end') if ev.get('series') else None} 行"
                                      f" (step {(ev.get('series') or {}).get('first')}..{(ev.get('series') or {}).get('last')}、"
                                      f"最終場 {(ev.get('series') or {}).get('final_field')})")
    return [
        ("verdict", str(g.get("verdict"))),
        ("fail_class", str(g.get("fail_class"))),
        ("acceptance_version", str(g.get(G.ACCEPTANCE_VERSION_KEY) or "(なし)")),
        ("rc", str(g.get("rc"))),
        ("保存場 field", f"{_ok(fld.get('ok'))}" + (f" 欠け {fld.get('missing')}" if fld.get("missing") else "")),
        ("残差 residual", f"{_ok(res.get('ok'))} {res.get('verdict')}"),
        ("残差の入力 residual_evidence", evs),
        ("4 量の頭打ち steadiness", f"{_ok(st.get('ok'))} " + "・".join(f"{q} {v.get('verdict')}" for q, v in ser.items())
         + (f" 欠けた量 {st.get('missing')}" if st.get("missing") else "")),
        ("最終場の床 floors", _ok(fl.get("ok"))),
        ("#1 床事象 floor_events", f"{_ok(fe.get('ok'))} {fe.get('verdict')} 区間 {fe.get('window')} 更新 {fe.get('n_updates_checked')}"
         f" 必須 {fe.get('required')}" if fe else "(なし)"),
        ("力の部品 force_parts", "(なし)" if "force_parts" not in g else
         f"{_ok((g.get('force_parts') or {}).get('ok'))} 不完全 {(g.get('force_parts') or {}).get('n_incomplete')}"),
        ("判定区間 judgment (#2・#3)", "(なし)" if jg is None else f"{jg.get('verdict')} {jg.get('fail_class') or ''}".strip()),
        ("#2 温度 temperature", "(なし)" if jg is None else f"{tp.get('verdict')} 窓 {tp.get('window')} 標本 {tp.get('n_window')}"),
        ("#3 力の窓 force_window", "(なし)" if jg is None else
         f"{fw.get('verdict')} 窓 {fw.get('window')} 標本 {fw.get('n_tail')}/{fw.get('n_prev')}"),
        ("残差の桁 residual_scale", _ok(rs.get("ok")) if rs else "(なし)"),
        ("理由 reasons", "; ".join((g.get("reasons") or [])[:4]) or "(なし)"),
    ]


def write_ab(out: Path, run: str, old_m: dict | None, new_m: dict, old_rec: dict) -> dict:
    a, b = dict(gates_items((old_m or {}).get("gates"))), dict(gates_items(new_m.get("gates")))
    keys = list(dict.fromkeys(list(b) + list(a)))
    rows = [{"item": k, "old": a.get(k, "(なし)"), "current": b.get(k, "(なし)"), "same": a.get(k) == b.get(k)} for k in keys]
    md = [f"# A/B: 受理検査だけを旧版 (保存された metrics.json) / 現行版 ({G.ACCEPTANCE_VERSION}) で変える — {run}", "",
          f"- 旧: `{old_rec.get('path')}` (sha256 {old_rec.get('sha256')}, mtime {old_rec.get('mtime')}) — run の時点の collect の出力",
          f"- 現行: `{out / 'metrics.json'}` — `runner_sern3d.collect({MERGED_YAML}, run, out_dir=...)`", "",
          "| 項目 | 旧 | 現行 | 同じ |", "|---|---|---|---|"]
    md += [f"| {r['item']} | {str(r['old']).replace('|', '/')} | {str(r['current']).replace('|', '/')} | {'=' if r['same'] else '≠'} |"
           for r in rows]
    jg = (new_m.get("gates") or {}).get("judgment") or {}
    fwq = (jg.get("force_window") or {}).get("quantities") or {}
    if fwq:
        md += ["", "判定区間 #3 の窓 (現行、末尾 W step の平均 m・a = max|q − m|・前窓との差 d):", "",
               "| 量 | 平均 | a | d | 許容 | 定常 | 合否 |", "|---|---|---|---|---|---|---|"]
        md += [f"| {q} | {v.get('mean'):.7g} | {v.get('a'):.2e} | {v.get('d'):.2e} | {v.get('tol'):.0e} | {v.get('steady')} | "
               f"{'OK' if v.get('ok') else 'NG'} |" for q, v in fwq.items()]
    (out / AB_MD).write_text("\n".join(md) + "\n")
    return {"rows": rows}


# ----------------------------------------------------------------------------------------------- run 1 本
def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, ensure_ascii=False, default=str))


def _summary_of(rec: dict) -> dict:
    """まとめ・結果 A/B が読む欄 (現行の gates から)。"""
    g = rec.get("current_gates") or {}
    jg = g.get("judgment") or {}
    fe = g.get("floor_events") or {}
    ev = G.residual_evidence_summary(g.get("residual")) if g else None
    fwq = (jg.get("force_window") or {}).get("quantities") or {}
    return {"status": rec.get("status"), "verdict": g.get("verdict"), "fail_class": g.get("fail_class"),
            "acceptance_version": g.get(G.ACCEPTANCE_VERSION_KEY), "reasons": (g.get("reasons") or [])[:4],
            "judgment": {"verdict": jg.get("verdict"), "fail_class": jg.get("fail_class"),
                         "temperature": (jg.get("temperature") or {}).get("verdict"),
                         "force_window": (jg.get("force_window") or {}).get("verdict"),
                         "floor_events": fe.get("verdict"), "floor_window": fe.get("window")} if g else None,
            "force_window": {q: {k: v.get(k) for k in ("mean", "a", "d", "tol", "steady", "ok")} for q, v in fwq.items()},
            "residual_evidence_input_ok": (ev or {}).get("input_ok") if ev else None,
            "missing_inputs": (rec.get("inputs") or {}).get("missing") or [],
            "stop": rec.get("stop"), "run_files_unchanged": rec.get("run_files_unchanged"), "error": rec.get("error")}


def rejudge_one(spec: Spec, a) -> dict:
    case_dir = Path(a.case_dir).resolve() if a.case_dir else C46
    rd = case_dir / spec.run
    out = rd / a.tag
    t0 = time.time()
    rec = {"version": VERSION, "start": now(), "spec": dataclasses.asdict(spec), "run_dir": str(rd), "out_dir": str(out),
           "acceptance_version_current": G.ACCEPTANCE_VERSION, "git": git_head(), "status": None, "elapsed_s": {}}
    if not rd.is_dir():
        rec.update({"status": "ERROR", "error": f"run が無い: {rd}"})
        return rec
    if out.exists():
        rec.update({"status": "ERROR", "error": f"再判定の出力先 {out} が既にある (上書きしない。--tag で別の名前にする)", "exists": True})
        return rec
    ts = time.time()
    before = snapshot(rd, out)
    rec["elapsed_s"]["snapshot_before"] = round(time.time() - ts, 1)
    out.mkdir()
    rec["run_files"] = {"n": len(before), "bytes": sum(v.get("size", 0) for v in before.values())}
    for name in ("IC_FROM.txt", "RESTART_FROM.txt", "IC_TRANSFER.txt", "RUN_PROVENANCE.txt"):
        if (rd / name).is_file():
            rec.setdefault("run_records", {})[name] = (rd / name).read_text(errors="replace")[-2000:]
    try:
        # (0) 問題 YAML
        merged, rec["yaml"] = merged_problem(spec, case_dir, out)
        p = load_problem(merged)
        try:
            jo, jp = R2.judgment_options(p.evaluate), R2.judgment_plan(p.evaluate)
        except ValueError as e:
            raise Stop("yaml", [f"判定区間の段の設定を作れない ({e})"])
        if jp is None:
            raise Stop("yaml", ["問題 YAML に evaluate.judgment_steps が無い (判定区間の段の判定を有効にできない)"])
        rec["judgment_plan"] = {k: jp[k] for k in ("steps", "window_steps", "extend_once", "out_interval", "required")}
        # (1) 来歴
        rec["provenance"] = check_provenance(spec, rd)
        # (2) 入力
        rec["inputs"] = check_inputs(rd, jp)
        if rec["inputs"]["missing"]:
            raise Stop("inputs", [f"入力が欠けている: {m}" for m in rec["inputs"]["missing"]])
        # (3) 設定・格子・prepare_info の照合 (現行 runner の prepare)
        free = shutil.disk_usage(rd).free / 2 ** 30
        rec["free_gib_before_prepare"] = round(free, 1)
        if free < a.min_free_gb:
            rec.update({"status": "ERROR", "error": f"空き {free:.1f} GiB < {a.min_free_gb} GiB (現行の prepare の格子に要る)"})
            return rec
        prep = out / PREP_DIR
        rec["current_prepare"] = {"dir": str(prep), "converter": R2.converter_identity(), "mesher": R2.mesher_version(3)}
        try:
            ts = time.time()
            new_info = run_prepare(merged, spec.op, prep)
            rec["elapsed_s"]["prepare"] = new_info.pop("_elapsed_s", round(time.time() - ts, 1))
            mq = prep / "MESH_QUALITY.txt"
            rec["current_prepare"]["mesh_quality"] = (mq.read_text().strip().splitlines() or [""])[-1] if mq.is_file() else None
            ts = time.time()
            rec["mesh"] = compare_mesh(rd, new_info)
            rec["elapsed_s"]["mesh_identity"] = round(time.time() - ts, 1)
            rec["config"] = compare_configs(rd, prep, jo, out)
            rec["prepare_info"] = compare_info(rd, new_info)
        finally:      # 失敗・中断でも現行の prepare の格子・中間 (再生成できる) は消す
            rec["current_prepare"]["removed"] = _cleanup_prepare(prep) if prep.is_dir() else []
        for kind in ("mesh", "config", "prepare_info"):
            if rec[kind]["disallowed"]:
                raise Stop(kind, rec[kind]["disallowed"])
        # (4) 温度の監視の列
        ts = time.time()
        rec["te_monitor"] = check_te_monitor(rd, rec["inputs"]["final_step"], out)
        rec["elapsed_s"]["te_monitor"] = round(time.time() - ts, 1)
        # (5) 判定 (現行の受理検査)
        ts = time.time()
        o = p.raw.get("opt") or {}
        new_m = R3.collect(str(merged), str(rd), out_dir=str(out), require_residual_pass=bool(o.get("require_residual_pass", False)))
        rec["elapsed_s"]["collect"] = round(time.time() - ts, 1)
        rec["current_gates"] = new_m["gates"]
        rec["current_metrics"] = {"path": str(out / "metrics.json"), "sha256": sha256(out / "metrics.json"), "forge_rc": new_m.get("forge_rc"),
                                  "n_snapshots": new_m.get("n_snapshots"), "eval_method_id": new_m.get(R2.EVAL_METHOD_ID)}
        # (6) A/B
        oldp = rd / "metrics.json"
        old_m = json.loads(oldp.read_text()) if oldp.is_file() else None
        rec["old_metrics"] = ({"path": str(oldp), "sha256": sha256(oldp), "mtime": time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(oldp.stat().st_mtime)),
                               "forge_rc": old_m.get("forge_rc"), "n_snapshots": old_m.get("n_snapshots")}
                              if old_m is not None else {"path": None, "note": "run に metrics.json が無い (旧版の列は「なし」)"})
        rec["ab"] = write_ab(out, spec.run, old_m, new_m, rec["old_metrics"])
        rec["old_gates_items"] = dict(gates_items((old_m or {}).get("gates")))
        rec["status"] = "JUDGED"
    except Stop as e:
        rec.update({"status": "INDETERMINATE", "stop": {"kind": e.kind, "reasons": e.reasons}})
    except Exception as e:  # noqa: BLE001 (台本の失敗は記録して次の run へ)
        rec.update({"status": "ERROR", "error": f"{type(e).__name__}: {e}", "traceback": traceback.format_exc()[-3000:]})
    finally:
        ts = time.time()
        after = snapshot(rd, out)
        rec["elapsed_s"]["snapshot_after"] = round(time.time() - ts, 1)
        ch = snapshot_changes(before, after)
        rec["run_files_unchanged"] = not ch
        rec["run_files_changes"] = ch[:50]
        if ch and rec["status"] != "ERROR":
            rec.update({"status": "ERROR", "error": f"run の既存ファイルが処理の前後で変わった ({len(ch)} 件): {ch[:3]}"})
        rec["elapsed_s"]["total"] = round(time.time() - t0, 1)
        rec["end"] = now()
        rec["summary"] = _summary_of(rec)
        if out.is_dir():
            _write_json(out / INFO_JSON, rec)
            _write_run_summary(out, rec)
    return rec


def _write_run_summary(out: Path, rec: dict) -> None:
    s = rec["summary"]
    lines = [f"REJUDGE {rec['spec']['run']} ({VERSION}, {rec.get('end')}) status {rec['status']}",
             f"現行の受理検査 ({G.ACCEPTANCE_VERSION}): verdict {s.get('verdict')} fail_class {s.get('fail_class')}",
             f"判定区間: {s.get('judgment')}"]
    if rec.get("stop"):
        lines += [f"判定不能 ({rec['stop']['kind']}):"] + [f"  - {r}" for r in rec["stop"]["reasons"]]
    if rec.get("error"):
        lines.append(f"ERROR: {rec['error']}")
    for q, v in (s.get("force_window") or {}).items():
        lines.append(f"  {q:15s} 平均 {v.get('mean'):.7g}  a {v.get('a'):.2e}  d {v.get('d'):.2e}  {v.get('steady')} {'OK' if v.get('ok') else 'NG'}")
    lines.append(f"run の既存ファイル: {'処理の前後で同じ' if rec.get('run_files_unchanged') else '変わった ' + str(rec.get('run_files_changes'))}")
    lines.append(f"所要 [s]: {rec.get('elapsed_s')}")
    (out / SUMMARY_TXT).write_text("\n".join(lines) + "\n")


# ----------------------------------------------------------------------------------------------- まとめ
def necessary_failures(rec: dict | None) -> list:
    """結果 A の必要条件 (現行の受理検査で PASS・版・判定区間 PASS・残差の入力・元ファイル不変) に満たないものの理由。"""
    if rec is None:
        return ["再判定の記録が無い (未実施)"]
    s = rec.get("summary") or _summary_of(rec)
    if s.get("status") == "INDETERMINATE":
        return [f"判定不能 ({(s.get('stop') or {}).get('kind')}): " + "; ".join((s.get("stop") or {}).get("reasons", [])[:6])]
    if s.get("status") != "JUDGED":
        return [f"{s.get('status')}: {s.get('error')}"]
    out = []
    if s.get("verdict") != "PASS":
        out.append(f"現行の受理検査で {s.get('verdict')} ({s.get('fail_class')}): " + "; ".join(s.get("reasons") or []))
    if s.get("acceptance_version") != G.ACCEPTANCE_VERSION:
        out.append(f"受理検査の版 {s.get('acceptance_version')} が現行 {G.ACCEPTANCE_VERSION} でない")
    if (s.get("judgment") or {}).get("verdict") != "PASS":
        out.append(f"判定区間 (#2・#3) が PASS でない ({(s.get('judgment') or {}).get('verdict')})")
    if (s.get("judgment") or {}).get("floor_events") != "PASS":
        out.append(f"#1 床事象が PASS でない ({(s.get('judgment') or {}).get('floor_events')})")
    if s.get("residual_evidence_input_ok") is not True:
        out.append("残差の入力の検査済み要約が合格でない")
    if s.get("run_files_unchanged") is not True:
        out.append("run の既存ファイルが処理の前後で変わった")
    return out


def _u(v: dict) -> float:
    return max(float(v["a"]), float(v["d"]), U_FLOOR)


def design_diff(specs, recs: dict) -> dict:
    """現行の verdict が PASS の行から g3 の設計差 Δᵢ (L_sw 1.0 − 0.8)・ΔJ と m10_on の格子感度の参考値 (解釈しない)。"""
    def fw(spec):
        r = recs.get(spec.run) if spec else None
        s = (r or {}).get("summary") or {}
        ok = s.get("status") == "JUDGED" and s.get("verdict") == "PASS" and all(q in (s.get("force_window") or {}) for q in QUANTS)
        return (s.get("force_window") if ok else None), (s.get("verdict") if r else "未実施")
    by = {(s.grid, s.op, round(float(s.L_sw), 6)): s for s in specs}
    out = {"q": OBJ, "weights": WEIGHTS, "per_op": {}, "dJ": None, "U_time_weighted": None, "missing": [], "G10": None}
    dJ, uJ, n = 0.0, 0.0, 0
    for op in ("m10_on", "m6_on", "m4_off"):
        s0, s1 = by.get(("g3", op, 0.8)), by.get(("g3", op, 1.0))
        (f0, v0), (f1, v1) = fw(s0), fw(s1)
        if s0 is None or s1 is None or f0 is None or f1 is None:
            out["missing"].append(f"{op}: g3 の L_sw 0.8 ({s0.run if s0 else '対応表に無い'}: {v0}) / 1.0 ({s1.run if s1 else '対応表に無い'}: {v1})"
                                  " の両方が現行で PASS でない")
            continue
        row = {"run_0.8": s0.run, "run_1.0": s1.run, "quantities": {}}
        for q in QUANTS:
            row["quantities"][q] = {"q_0.8": f0[q]["mean"], "q_1.0": f1[q]["mean"], "delta": f1[q]["mean"] - f0[q]["mean"],
                                    "U": _u(f0[q]) + _u(f1[q])}
        out["per_op"][op] = row
        dJ += WEIGHTS[op] * row["quantities"][OBJ]["delta"]
        uJ += WEIGHTS[op] * row["quantities"][OBJ]["U"]
        n += 1
    if n == len(WEIGHTS):
        out["dJ"], out["U_time_weighted"] = dJ, uJ
    g4 = [s for s in specs if s.grid == "g4" and s.op == "m10_on"]
    g4_10 = by.get(("g4", "m10_on", 1.0))
    s3 = by.get(("g3", "m10_on", 0.8))
    G10 = {"status": "未確定", "why": ("Δ₁₀,g4 に要る g4・L_sw 1.0・m10_on の判定区間が無い (g4 は L_sw 0.8 だけ)" if g4_10 is None
                                      else "g4 の L_sw 1.0・m10_on の判定区間は対応表にあるが、式の値はここでは出さない"),
           "formula": "G₁₀ = |Δ₁₀,g4 − Δ₁₀,g3| + U₁₀,g4 + U₁₀,g3"}
    for s4 in g4:
        if round(float(s4.L_sw), 6) != 0.8 or s3 is None:
            continue
        (f4, v4), (f3, v3) = fw(s4), fw(s3)
        if f4 is None or f3 is None:
            G10["reference_0.8"] = f"未算出 ({s4.run}: {v4} / {s3.run}: {v3})"
        else:
            G10["reference_0.8"] = {"run_g4": s4.run, "run_g3": s3.run, "q_g4": f4[OBJ]["mean"], "q_g3": f3[OBJ]["mean"],
                                    "q_g4_minus_g3": f4[OBJ]["mean"] - f3[OBJ]["mean"], "u_g4": _u(f4[OBJ]), "u_g3": _u(f3[OBJ])}
    out["G10"] = G10
    return out


def _yaml_brief(y: dict | None, spec: Spec) -> dict:
    """まとめに載せる問題 YAML の対応 (どの YAML のどの値を使ったか)。記録が無ければ対応表の名前だけ。"""
    if not isinstance(y, dict):
        return {"te_wake": Path(spec.te_yaml).name, "production": Path(spec.prod_yaml).name, "keys": None}
    return {"te_wake": Path(y["te_wake"]["path"]).name, "te_wake_sha256": y["te_wake"]["sha256"][:12],
            "te_wake_commit": (y["te_wake"].get("git_last_commit") or "")[:8],
            "production": Path(y["production"]["path"]).name, "production_sha256": y["production"]["sha256"][:12],
            "production_commit": (y["production"].get("git_last_commit") or "")[:8],
            "keys": {k: v.get("used") for k, v in (y.get("keys_from_production") or {}).items()},
            "keys_te_wake": {k: v.get("te_wake") for k, v in (y.get("keys_from_production") or {}).items()},
            "n_other_differences": len(y.get("other_differences_te_vs_production") or [])}


def summarize(specs, recs: dict) -> dict:
    rows, fails = [], []
    for s in specs:
        r = recs.get(s.run)
        sm = ((r or {}).get("summary") or (_summary_of(r) if r else {}))
        nf = necessary_failures(r)
        if nf:
            fails.append({"run": s.run, "reasons": nf})
        rows.append({"run": s.run, "L_sw": s.L_sw, "op": s.op, "grid": s.grid, "status": (r or {}).get("status", "未実施"),
                     "verdict": sm.get("verdict"), "fail_class": sm.get("fail_class"), "judgment": sm.get("judgment"),
                     "force_window": sm.get("force_window") or {}, "missing_inputs": sm.get("missing_inputs") or [],
                     "stop": sm.get("stop"), "necessary_ok": not nf, "record": (r or {}).get("out_dir"),
                     "record_end": (r or {}).get("end"), "git": ((r or {}).get("git") or {}).get("head"),
                     "yaml": _yaml_brief((r or {}).get("yaml"), s)})
    return {"version": VERSION, "time": now(), "git": git_head(), "acceptance_version_current": G.ACCEPTANCE_VERSION,
            "result": "A" if not fails else "B",
            "result_text": ("結果 A: 対象の全 run が現行の受理検査の必要条件を満たす" if not fails else
                            "結果 B: 欠落・判定不能・不合格の run がある (版の付け替えだけで救済しない)"),
            "failures": fails, "rows": rows, "design_diff": design_diff(specs, recs)}


def _fmt_q(v) -> str:
    return "—" if not v else f"{v.get('mean'):.7g} ({v.get('a'):.1e} / {v.get('d'):.1e})"


def write_summary(summary: dict, md: Path, js: Path) -> None:
    _write_json(js, summary)
    L = [f"# R7b ④-0 既存 {len(summary['rows'])} 本の再判定 ({VERSION}) — 機械的な表示 (解釈は呼び出し側)", "",
         f"- 生成 {summary['time']}、git {summary['git'].get('head')}、現行の受理検査の版 `{summary['acceptance_version_current']}`",
         f"- **{summary['result_text']}**"]
    for f in summary["failures"]:
        L.append(f"  - `{f['run']}`: " + " / ".join(f["reasons"]))
    L += ["", "## 対象の表", "",
          "判定区間の列は #1 床事象 / #2 温度 / #3 力の窓。力は判定区間の末尾 W step の「平均 (a / d)」(a = max|q − 平均|、d = 前窓との差)。", "",
          "| run | L_sw | 作動点 | 格子 | 状態 | 現行 verdict | fail_class | 判定区間 #1/#2/#3 | C_T | C_T_with_shear | C_L | C_M | 欠けた入力・判定不能の理由 |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in summary["rows"]:
        j = r.get("judgment") or {}
        jtxt = f"{j.get('floor_events')}/{j.get('temperature')}/{j.get('force_window')}" if j else "—"
        why = r["missing_inputs"] or ((r.get("stop") or {}).get("reasons") or [])[:4]
        L.append(f"| `{r['run']}` | {r['L_sw']} | {r['op']} | {r['grid']} | {r['status']} | {r['verdict'] or '—'} | {r['fail_class'] or '—'} | "
                 f"{jtxt} | " + " | ".join(_fmt_q(r["force_window"].get(q)) for q in QUANTS)
                 + f" | {'; '.join(str(x) for x in why).replace('|', '/') or '—'} |")
    dd = summary["design_diff"]
    L += ["", f"## 参考: g3 の設計差 (L_sw 1.0 − 0.8、現行の verdict が PASS の行だけ、q = {dd['q']})", "",
          "| 作動点 | run 0.8 | run 1.0 | 量 | q(0.8) | q(1.0) | Δ | 時間幅 U = u₀.₈ + u₁.₀ |", "|---|---|---|---|---|---|---|---|"]
    for op, row in dd["per_op"].items():
        for q, v in row["quantities"].items():
            L.append(f"| {op} | `{row['run_0.8']}` | `{row['run_1.0']}` | {q} | {v['q_0.8']:+.7f} | {v['q_1.0']:+.7f} | {v['delta']:+.4e} | {v['U']:.2e} |")
    for m in dd["missing"]:
        L.append(f"- 未算出: {m}")
    if dd["dJ"] is not None:
        L.append(f"- ΔJ = 0.5Δ₆ + 0.3Δ₁₀ + 0.2Δ₄ ({dd['q']}) = **{dd['dJ']:+.4e}**、うち時間変動幅の加重 Σ wᵢUᵢ,g3 = {dd['U_time_weighted']:.2e} "
                 "(E_J = Σ wᵢ(Uᵢ,g3 + Gᵢ + Dᵢ) の G・D は未測定なので E_J は出さない)")
    else:
        L.append("- ΔJ: 未算出 (3 作動点の対がそろわない)")
    g = dd["G10"] or {}
    L.append(f"- m10_on の格子感度 {g.get('formula')}: **{g.get('status')}** — {g.get('why')}")
    ref = g.get("reference_0.8")
    if isinstance(ref, dict):
        L.append(f"  - 参考 (L_sw 0.8 の絶対値): q(g4) − q(g3) = {ref['q_g4_minus_g3']:+.4e} (`{ref['run_g4']}` {ref['q_g4']:+.7f} / "
                 f"`{ref['run_g3']}` {ref['q_g3']:+.7f})、u_g4 {ref['u_g4']:.2e}・u_g3 {ref['u_g3']:.2e}")
    elif ref:
        L.append(f"  - 参考 (L_sw 0.8): {ref}")
    L += ["", "## 使った問題 YAML (te_wake の YAML を基本に、evaluate の判定区間のキーだけ生産 YAML の値)", "",
          "| run | te_wake の YAML (sha256・最終コミット) | 生産 YAML (sha256・最終コミット) | 生産 YAML から取った値 (te_wake の YAML の値) | 他に違うキーの数 |",
          "|---|---|---|---|---|"]
    for r in summary["rows"]:
        y = r.get("yaml") or {}
        kv = "・".join(f"{k} {v} ({(y.get('keys_te_wake') or {}).get(k)})" for k, v in (y.get("keys") or {}).items()) or "(記録なし)"
        L.append(f"| `{r['run']}` | {y.get('te_wake')} ({y.get('te_wake_sha256', '—')}・{y.get('te_wake_commit', '—')}) | "
                 f"{y.get('production')} ({y.get('production_sha256', '—')}・{y.get('production_commit', '—')}) | {kv} | "
                 f"{y.get('n_other_differences', '—')} |")
    L += ["", "## 各 run の記録", ""]
    for r in summary["rows"]:
        L.append(f"- `{r['run']}`: `{r['record']}` ({INFO_JSON}・{AB_MD}・{SUMMARY_TXT}・metrics.json・{CONFIG_DIFF_TXT}) "
                 f"終了 {r['record_end']}・git {r['git']}")
    md.write_text("\n".join(L) + "\n")


def load_records(specs, case_dir: Path, tag: str) -> dict:
    out = {}
    for s in specs:
        p = case_dir / s.run / tag / INFO_JSON
        if p.is_file():
            out[s.run] = json.loads(p.read_text())
    return out


# ----------------------------------------------------------------------------------------------- main
def parse_spec(text: str) -> Spec:
    f = [x.strip() for x in text.split(",")]
    if len(f) != 6:
        raise argparse.ArgumentTypeError("--spec は RUN,L_SW,OP,GRID,TE_YAML,PROD_YAML")
    return Spec(f[0], float(f[1]), f[2], f[3], f[4], f[5], "(--spec)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="plan tooling-nozzle-sern-chain §5.1 R7b ④-0: 既存 7 本を現行の受理検査で再判定 (forge 0 step)")
    ap.add_argument("--only", nargs="*", default=None, help="この run だけ処理する (まとめは記録のある全対象から)")
    ap.add_argument("--spec", action="append", type=parse_spec, default=None,
                    help="対応表を差し替える RUN,L_SW,OP,GRID,TE_YAML,PROD_YAML (複数可。台本の試験用)")
    ap.add_argument("--case-dir", default=None, help="run と YAML の置き場 (既定 case/46.sern_design。台本の試験用)")
    ap.add_argument("--tag", default=TAG, help=f"再判定のサブディレクトリ名 (既定 {TAG}、既にあれば失敗)")
    ap.add_argument("--summary-out", default=None, help="まとめの出力 (拡張子なしのパス。既定 <case>/R7B_<tag>、.md と .json。既にあれば失敗)")
    ap.add_argument("--summary-only", action="store_true", help="run を処理せず、既存の再判定の記録からまとめだけ書く")
    ap.add_argument("--no-summary", action="store_true", help="まとめを書かない")
    ap.add_argument("--min-free-gb", type=float, default=8.0, help="現行の prepare の前に要る空き容量 [GiB] (格子・変換の中間)")
    a = ap.parse_args(argv)
    case_dir = Path(a.case_dir).resolve() if a.case_dir else C46
    specs = tuple(a.spec) if a.spec else RUNS
    base = Path(a.summary_out) if a.summary_out else case_dir / f"R7B_{a.tag}"
    md, js = base.with_name(base.name + ".md"), base.with_name(base.name + ".json")
    if not a.no_summary and (md.exists() or js.exists()):
        print(f"REJUDGE STOP: まとめの出力先 {md} / {js} が既にある (上書きしない。--summary-out で別の名前にする)", file=sys.stderr)
        return 1
    rc, exists, errors = 0, [], []
    if not a.summary_only:
        todo = [s for s in specs if a.only is None or s.run in a.only]
        if a.only is not None and len(todo) != len(a.only):
            print(f"REJUDGE STOP: --only の run が対応表に無い: {sorted(set(a.only) - {s.run for s in specs})}", file=sys.stderr)
            return 1
        print(f"REJUDGE_START {now()} git {git_head()} 対象 {len(todo)} 本", flush=True)
        for s in todo:
            print(f"== {s.run} (L_sw {s.L_sw}・{s.op}・{s.grid}) {now()}", flush=True)
            rec = rejudge_one(s, a)
            sm = rec.get("summary") or {}
            print(f"   status {rec['status']} verdict {sm.get('verdict')} fail_class {sm.get('fail_class')} "
                  f"判定区間 {sm.get('judgment')} 所要 {rec.get('elapsed_s', {}).get('total')} s", flush=True)
            if rec.get("stop"):
                print("   判定不能: " + " / ".join(rec["stop"]["reasons"][:4]), flush=True)
            if rec["status"] == "ERROR":
                print(f"   ERROR: {rec.get('error')}", flush=True)
                (exists if rec.get("exists") else errors).append(s.run)
        rc = 1 if errors else (2 if exists else 0)
    if not a.no_summary:
        summary = summarize(specs, load_records(specs, case_dir, a.tag))
        write_summary(summary, md, js)
        print(f"REJUDGE {summary['result_text']}")
        for f in summary["failures"]:
            print(f"   {f['run']}: {' / '.join(f['reasons'])}")
        print(f"まとめ: {md} / {js}")
    print(f"REJUDGE_DONE {now()} rc {rc}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
