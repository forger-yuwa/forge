"""plan tooling-nozzle-isothermal-wall-chain §5.1 #31・§6「V-g0c-chain の事前登録」(commit 68d56723、§4.8、codex diagnose
`notes/reviews/2026-10-13-g0c-chain-design-diagnose.md`) の実行: 設計チェーンの準備 (prepare_ns) → 前提の検査 (1)〜(6) → 段階起動 full
(run_staged_ns) → 本段だけに停止規則の見張り (m9_watch.py) → 判定 (vgc_judge.py)。B で、vgc_run.sh (環境変数・B の確認) から起動する。

  コード : ~/forge-g0cchain (vgc_prep.sh が作る、この台本を commit した後の commit の作業ツリー) の design/ を sys.path の先頭にして forge_design を読む。
           FORGE_ROOT は forge_design の場所から決まるので、run_case.sh・interp_field.py・check_mesh_quality.py・stage_manifest.py も同じ作業ツリーのもの。
  forge  : ~/forge-fgeom7-fp64 (sha256 7e95e5bbd8f50c9e…)、変換器は同じビルドの convertGmshToForge (711f4438…)。
  準備   : prepare_ns(問題, run, ic_from=run_0179_ns_n012_N2_ext (最後の res が res_20000.h5)、delta_r_csv=_band_ab/prod_confirm/prep/delta_r_initial.csv、
           offset="radial"、euler_ref=run_0174_euler_v5d_M_r1)。cfl_main・implicit_relax・nsteps は渡さない (問題 YAML の evaluate が決める)。
  前提   : (1) nozzle.msh の sha256 が G0c (285f0fa4ecc0ed1c)。(2) nozzle.h5 の格子 (MESH・CELLS・PLANES・VIZMESH・BCONDS の番号) と VALUE/wall_dist が
           run_0804_sp_g0c の入力 (= _g0c/q32_init.h5、6f8fcda6d97482e3) とビット一致。(3) check_mesh_quality.py --ar-max 10000 が PASS、AR > 5000 のセルの
           数・位置・スキューを記録。(4) ic_from の res・準備の後の nozzle.h5・本段の初期の場の sha256 を記録 (res は登録の fb65c9b37d9bb475 と照合)。
           (5) 本段の solverConfig.yaml が B の run_0717_c50_d1/solverConfig.yaml と nStepOuter・outStepInterval・mesh.axisSegmentRWeight の明示 (既定 1)
           以外でキーごとに一致、S1・S2 にライン陰解法・ITJ・上限・リミッタの基準値・extraFields が無く implicitRelax 0.7。(6) bcondConfig.yaml・
           species_meta.yaml・resolved_species_*.yaml が run_0804_sp_g0c と一致 (forge が起動時に書く記録で準備の時点に無いものは、本段の直前に照合し直す)。
           どれかが外れたら forge を回さない (判定不能)。
  実行   : START_EPOCH (準備の前の時刻) と 30 秒ごとの timing.tsv (sp_arm.sh と同じ列: 時刻・残差の履歴の最終 step・forge の数・ほかの forge の cwd)。
           run_staged_ns(run, stages="full")。本段の通知 MAIN_STAGE_STARTED を書く関数 (runner_axismach._write_main_stage_started) の直後に、通知の中身・
           前の段の出力が残っていないこと・(6) を照合してから m9_watch.py を §6.11 の引数で起動する (照合が外れたら本段の forge を起動させない)。
           観測だけの包み: run_forge (段ごとの rc・時刻 → vgc_stages.jsonl)、stage_gate (S1・S2 の最終場の有限性・ρ・P・T → vgc_stage_gate.jsonl)。
           どちらも元の関数の戻り値をそのまま返す (チェーンの手順は変えない)。
  見張り : ディスクの空きが 2 GiB を切ったら、cwd がこの run の forge だけを止める (DISK_STOP)。
  終わり : RUN_RC (run_staged_ns の戻り、例外なら EXC)、m9_watch.json の状態を RUN_END.json と vgc_run.log に書き、vgc_judge.py を回す。

usage (B、vgc_run.sh 経由): bash vgc_run.sh --commit <作業ツリーの commit> [--wait-max 3600]
       回す前の確かめだけ (run_0810 を使わない、forge を回さない): bash vgc_run.sh --commit <commit> --dry-prep _vgc_dry_1
"""
import argparse
import hashlib
import importlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent            # B: ~/forge-wallfit/case/45.isobutane_m6_d155 (vgc_prep.sh が作業ツリーから写す)
WT = Path.home() / "forge-g0cchain"
CASE = "case/45.isobutane_m6_d155"
PLAN = "plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #31・§6 V-g0c-chain (事前登録 commit 68d56723)"
PROBLEM = "problem_d155_ns_prod_coldmesh_tw300_g0c_chain.yaml"
RUN = "run_0810_vgc_chain"
IC_RUN, IC_RES, IC_SHA = "run_0179_ns_n012_N2_ext", "res_20000.h5", "fb65c9b37d9bb475"
PREP = "_band_ab/prod_confirm/prep"
PREP_SHA = {"delta_r_initial.csv": "0bc9da5b6e8f0197", "wall_physical.csv": "31b4fd80818b4a0c"}
EULER_REF = "run_0174_euler_v5d_M_r1"
B0_REF = "run_0717_c50_d1"
CTRL, CTRL_INPUT, CTRL_SHA = "run_0804_sp_g0c", "_g0c/q32_init.h5", "6f8fcda6d97482e3"
MSH_SHA = "285f0fa4ecc0ed1c"
BIN = Path.home() / "forge-fgeom7-fp64/solver_density_cuda/build/forge"
BIN_SHA, CONV_SHA = "7e95e5bbd8f50c9e", "711f4438"
ENV_WANT = {"FORGE_CUDA_BLOCKSIZE": "128", "FORGE_OMEGA_BUDGET": "1", "FORGE_DIAG_COMMIT_LOSS": "500", "FORGE_BIN": str(BIN)}
NUMPY_WANT = "2.5.2"
M9_ARGS = ["--phase", "line", "--budget", "200000", "--consec", "2", "--interval", "2500", "--min-points", "9", "--keep-tail", "9"]
MAIN_STEPS, OUT_INT, RELAX, AR_MAX = 200000, 2500, 0.7, 10000
DISK_MIN = 2 * 1024 ** 3
PHYS_FILES = ("bcondConfig.yaml", "species_meta.yaml")       # + resolved_species_*.yaml
RS_RE = re.compile(r"^resolved_species_[0-9a-f]{16}(?:_[0-9a-f]{16})?\.yaml$")
# 格子の照合 (2): ここに入るデータセットは全部ビット一致が必須。BCONDS/*/VALUE/* (境界の値) と VALUE/* (wall_dist 以外) は記録だけ
MESH_GROUPS = ("MESH", "CELLS", "PLANES", "VIZMESH")


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def log(msg: str) -> None:
    print(f"[vgc_run {now()}] {msg}", flush=True)


def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(t: str) -> str:
    return hashlib.sha256(t.encode()).hexdigest()


def jdump(p, d) -> None:
    tmp = Path(str(p) + ".tmp")
    tmp.write_text(json.dumps(d, indent=1, ensure_ascii=False, default=str))
    os.replace(tmp, p)


def jappend(p, d) -> None:
    with open(p, "a") as f:
        f.write(json.dumps(d, ensure_ascii=False, default=str) + "\n")


def forge_pids():
    """(pid, cwd の realpath) の一覧 (pgrep -x forge)。"""
    out = subprocess.run(["pgrep", "-x", "forge"], capture_output=True, text=True).stdout.split()
    r = []
    for p in out:
        try:
            r.append((int(p), os.path.realpath(f"/proc/{p}/cwd")))
        except OSError:
            r.append((int(p), "?"))
    return r


def stop_own_forge(run: Path) -> bool:
    """cwd がこの run の forge だけを止める (SIGTERM → 180 s → SIGKILL)。pkill は使わない。"""
    rd = os.path.realpath(run)
    own = lambda: [p for p, c in forge_pids() if c == rd]  # noqa: E731
    for sig, wait in ((signal.SIGTERM, 180), (signal.SIGKILL, 60)):
        for p in own():
            try:
                os.kill(p, sig)
            except OSError:
                pass
        t0 = time.time()
        while own() and time.time() - t0 < wait:
            time.sleep(3)
    return not own()


# --- 生成器が足した明示すべきキー (plan tooling-generator-defaults-and-launch-gate §5.1 #11・#26; h370_run.py と同じ関数) -----------
# #11 以降の生成器は段の config に「明示すべきキー」を旧い実効値で足す (config_builder.finish_config → yaml_strict.insert_scalars:
# 既にある写像へのキーの挿入だけ)。prepare_ns は足した後の設定を書き、足したキーを RECIPE_CONTEXT.json の explicit_added
# {family, prestage, main} に残す。(5) の照合は足す前の設定で行い (照合の意味は足さない版と同じ)、足したのが記録どおりのキーの挿入だけで
# あることは、除いた設定に生成器の仕上げを当て直して元の文字列に戻ることで確かめる。
EXPLICIT_CONTEXT = "RECIPE_CONTEXT.json"
# 前の段の元の設定・S1・S2 に B0 を入れない検査 (5) は、キーの有無でなく実効値で見る (同 §5.1 #26・§6.1 result (段取り 1) M1; h370_run.py と
# 同じ表)。値 = B0 だけのキーを書いていないときの実効値 (ソルバの既定、solverConfig.cpp:372・389・390・453・713〜715; どれも B0 を切った値)。
# この値の明示 (#11 以降の生成器が前の段に足す lineImplicit・lineDtDirectional・implicitThermalJacobian 0) は書いていないのと同じ、
# それ以外の値 (型の違う値・真偽値・null を含む) と表に無い B0 のキーの明示は B0 が残っている。
B0_OFF = {"lineImplicit": 0, "lineDtDirectional": 0, "lineDtDirectionalCap": 0.0, "implicitThermalJacobian": 0,
          "limiterRoRef": 0.0, "limiterPRef": 0.0, "limiterARef": 0.0}


def b0_left(RA, doc, paths) -> list:
    """読み込んだ設定 doc の B0 だけのパス paths のうち、実効値が書いていないとき (B0_OFF) と違うもの (点区切り)。整数のキーは整数だけ、
    実数のキーは整数か実数を書いていないときの値と同じとみなす。"""
    out = []
    for p in paths:
        v = RA._cfg_path_get(doc, p)
        if v is RA._CFG_MISSING:
            continue
        off = B0_OFF.get(p[-1])
        same = off is not None and not isinstance(v, bool) and (
            (type(v) is int and v == off) if isinstance(off, int) else (type(v) in (int, float) and float(v) == off))
        if not same:
            out.append(".".join(p))
    return out


def generator_finish(RA):
    """生成器の段の config の仕上げ (config_builder.finish_config) があればその関数、無い版 (足さない旧い生成器) なら None。"""
    name = RA.__package__ + ".config_builder"
    try:
        return importlib.import_module(name).finish_config
    except ModuleNotFoundError as e:
        if e.name != name:
            raise
        return None


def explicit_record(run: Path) -> dict:
    """prepare_ns が RECIPE_CONTEXT.json に書いた、足した明示すべきキー {family, prestage, main}。無ければ {} (足さない旧い生成器。
    記録の無い準備に足したキーがあれば、除かずに比べるので照合で不成立になる)。"""
    try:
        rec = json.loads((Path(run) / EXPLICIT_CONTEXT).read_text())
    except (OSError, ValueError):
        return {}
    return (rec.get("explicit_added") or {}) if isinstance(rec, dict) else {}


def explicit_remove(text: str, keys) -> str:
    """keys (点区切りのキー) の項目を text から除く (insert_scalars の逆: flow 形式は前の値の直後から自分の値の終わりまで、block 形式は行ごと)。
    キーが無い・flow 形式の先頭の項目 (除くと書式が変わる) は ValueError。"""
    import yaml
    root = yaml.compose(text, Loader=yaml.SafeLoader)
    spans = []
    for key in keys or ():
        node, parts = root, key.split(".")
        for i, part in enumerate(parts):
            items = node.value if isinstance(node, yaml.MappingNode) else []
            j = next((j for j, (k, _) in enumerate(items) if isinstance(k, yaml.ScalarNode) and k.value == part), None)
            if j is None:
                raise ValueError(f"足したキー {key} が設定に無い")
            if i < len(parts) - 1:
                node = items[j][1]
                continue
            k, v = items[j]
            if node.flow_style and j > 0:
                spans.append((items[j - 1][1].end_mark.index, v.end_mark.index))
            elif not node.flow_style:
                e = text.find("\n", v.end_mark.index)
                spans.append((text.rfind("\n", 0, k.start_mark.index) + 1, len(text) if e < 0 else e + 1))
            else:
                raise ValueError(f"足したキー {key} が flow 形式の先頭の項目 (挿入でできる位置でない)")
    for s, e in sorted(spans, reverse=True):
        text = text[:s] + text[e:]
    return text


def explicit_strip(RA, text: str, added, bcond, family, tag=None) -> str:
    """足す前の設定 = text から added (足したキー {キー: 値}) を除いたもの。family が None (記録が無い) なら text のまま。
    除いた設定に生成器の仕上げ (finish_config) を当て直すと text に戻り、足すキーが added と同じことを確かめる (違えば ValueError)。"""
    if family is None:
        return text
    fin = generator_finish(RA)
    if fin is None:
        raise ValueError(f"足したキーの記録 ({family}) があるのに、読み込んだ生成器に config_builder が無い")
    base = explicit_remove(text, list(added or {}))
    again, info = fin(base, bcond, family, tag=tag)
    if again != text or set(info.get("explicit_added") or {}) != set(added or {}):
        raise ValueError(f"足したキーを除いた設定に生成器の仕上げを当て直しても元に戻らない ({family}; 記録 {sorted(added or {})}、"
                         f"仕上げ {sorted(info.get('explicit_added') or {})})")
    return base


# --- 前提の検査 (5): 設定 ------------------------------------------------------------------------------------------------
def check_configs(RA, main: str, pre: str, ref: str, explicit: dict | None = None, bcond: str | None = None) -> dict:
    """本段の設定 main・前の段の元の設定 pre (PRESTAGE_CONFIG)・B0 の参照 ref (run_0717 の solverConfig.yaml) を照合する。
    本段と ref の差として許すのは nStepOuter (こちら 200000)・outStepInterval (こちら 2500)・mesh.axisSegmentRWeight (片側だけが明示の 1)。
    S1・S2 は RA.ns_stage_configs(main, "full", cfg_pre=pre) で導き、B0 だけのキー (ライン陰解法・上限・ITJ・リミッタの基準値; 起動する設定の
    実効値で見る: 書いていないときの値 B0_OFF の明示は無いのと同じ、b0_left) と本段の extraFields が無く、implicitRelax 0.7 であること。
    戻り {ok, failures, diff_vs_b0_ref, stages: {tag: {sha256, cfl, ...}}}。
    explicit (生成器が足した明示すべきキーの記録 {family, prestage, main}; explicit_record) があれば、main・pre から足したキーを除いた
    足す前の設定で上の照合をし (段の値も足す前のもの)、起動する段の設定から足したキーを除くと足す前の段と同じことを確かめる。段の sha は
    起動する (足した後の) 文字列のもの。bcond は除いた設定に生成器の仕上げを当て直すときの bcondConfig.yaml の文字列。"""
    M = RA._CFG_MISSING
    fails = []
    ex = explicit or {}
    try:
        main_b = explicit_strip(RA, main, ex.get("main"), bcond, ex.get("family"))
        pre_b = explicit_strip(RA, pre, ex.get("prestage"), bcond, ex.get("family"))
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "failures": [f"(5) 生成器が足した明示すべきキー ({EXPLICIT_CONTEXT}) を除いて比べられない: {e!r}"],
                "diff_vs_b0_ref": {}, "stages": {}, "explicit_added": ex}
    dm, dr = RA._cfg_load(main_b), RA._cfg_load(ref)
    diff = RA._cfg_diff_paths(dr, dm)
    shown = {".".join(map(str, k)): ["<無し>" if a is M else a, "<無し>" if b is M else b] for k, (a, b) in diff.items()}
    for k, (a, b) in diff.items():
        if k in (("time", "last", "nStepOuter"), ("time", "outStepInterval")):
            continue
        if k == ("mesh", "axisSegmentRWeight") and {type(a), type(b)} == {int, type(M)} and 1 in (a, b):
            continue                                   # ソルバの既定 1 を片側だけが明示 (solverConfig.hpp:454)
        fails.append(f"(5) 本段の設定が {B0_REF} と違う: {'.'.join(map(str, k))} = {shown['.'.join(map(str, k))]} (参照, こちら)")
    for pth, want in ((("time", "last", "nStepOuter"), MAIN_STEPS), (("time", "outStepInterval"), OUT_INT)):
        if RA._cfg_path_get(dm, pth) != want:
            fails.append(f"(5) 本段の {'.'.join(pth)} = {RA._cfg_path_get(dm, pth)!r} ({want} であること)")
    b0_only = sorted(p for p in RA._MAIN_LINE_ONLY_PATHS if p != ("output",))
    stages = {}
    try:
        pre_st, main_f, main_want = RA.ns_stage_configs(main, "full", cfg_pre=pre)
        pre_st_b = RA.ns_stage_configs(main_b, "full", cfg_pre=pre_b)[0] if ex else pre_st     # 足す前の段 (照合はこちら)
        differ = [t for (t, c, *_), (_, cb, *_) in zip(pre_st, pre_st_b) if explicit_remove(c, ex.get("prestage")) != cb]
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "failures": fails + [f"(5) ns_stage_configs が失敗: {e!r}"], "diff_vs_b0_ref": shown, "stages": {},
                "explicit_added": ex}
    if main_f != main:
        fails.append("(5) ns_stage_configs が本段の設定を変えた")
    if [t for t, *_ in pre_st] != ["S1_soft", "S2_mid"] or [t for t, *_ in pre_st_b] != ["S1_soft", "S2_mid"]:
        fails.append(f"(5) 前の段が S1_soft・S2_mid でない: {[t for t, *_ in pre_st]}")
    if differ:
        fails.append(f"(5) 起動する段の設定から足したキーを除いても足す前の段と同じにならない: {differ}")
    out_pre = RA._cfg_path_get(RA._cfg_load(pre_b), ("output",))
    for tag, cfg, cfg_b in [("prestage_config", pre, pre_b)] + [(t, c, cb) for (t, c, *_), (_, cb, *_) in zip(pre_st, pre_st_b)]:
        d = RA._cfg_load(cfg_b)
        # B0 は起動する (足した後の) 設定の実効値で見る (書いていないときの値 B0_OFF の明示は B0 なし)。足す前の設定はキーを除いただけ
        left = b0_left(RA, RA._cfg_load(cfg), b0_only)
        if left:
            fails.append(f"(5) {tag} に B0 のキー {left} がある (書いていないときと違う実効値)")
        o = RA._cfg_path_get(d, ("output",))
        if o != out_pre or (isinstance(o, dict) and "extraFields" in o):
            fails.append(f"(5) {tag} の output が前の段の元の設定と違うか extraFields がある ({'<無し>' if o is M else o!r})")
        if RA._cfg_path_get(d, ("time", "deltaT", "implicitRelax")) != RELAX:
            fails.append(f"(5) {tag} の implicitRelax = {RA._cfg_path_get(d, ('time', 'deltaT', 'implicitRelax'))!r} ({RELAX} であること)")
        stages[tag] = {"sha256": sha256_text(cfg), **{k: (None if RA._cfg_path_get(d, p) is M else RA._cfg_path_get(d, p)) for k, p in (
            ("cfl", ("time", "deltaT", "cfl")), ("cfl_pseudo", ("time", "deltaT", "cfl_pseudo")), ("convMethod", ("space", "convMethod")),
            ("nStepInner", ("time", "nStepInner")), ("nStepOuter", ("time", "last", "nStepOuter")),
            ("outStepInterval", ("time", "outStepInterval")), ("implicitRelax", ("time", "deltaT", "implicitRelax")))}}
    stages["main"] = {"sha256": sha256_text(main), "want": main_want}
    return {"ok": not fails, "failures": fails, "diff_vs_b0_ref": shown, "stages": stages,
            "stage_tags_by_sha": {v["sha256"]: k for k, v in stages.items() if k != "prestage_config"}, "explicit_added": ex}


# --- 前提の検査 (6): bcond・物性のファイル --------------------------------------------------------------------------------
def check_phys_files(run: Path, ctrl: Path, late: bool) -> dict:
    """bcondConfig.yaml・species_meta.yaml・resolved_species_*.yaml を run_0804 と照合する。late=False (準備の直後) は、forge が起動時に
    書く記録で run にまだ無いものを「未作成」として残す (後で照合し直す)。late=True (本段の直前・判定) は全部が同じであること。"""
    names = sorted(set(PHYS_FILES) | {f.name for f in ctrl.iterdir() if RS_RE.match(f.name)})
    ours_rs = {f.name for f in run.iterdir() if RS_RE.match(f.name)}
    fails, rec, pending = [], {}, []
    for n in names:
        a, b = run / n, ctrl / n
        if not b.is_file():
            fails.append(f"(6) {CTRL}/{n} が無い"); rec[n] = "参照に無い"; continue
        if not a.is_file():
            if late or n in PHYS_FILES:
                fails.append(f"(6) {n} がこの run に無い"); rec[n] = "無い"
            else:
                pending.append(n); rec[n] = "未作成 (forge が起動時に書く記録; 本段の直前に照合し直す)"
            continue
        same = a.read_bytes() == b.read_bytes()
        rec[n] = {"same": same, "sha256": sha256_file(a)[:16], "ref_sha256": sha256_file(b)[:16]}
        if not same:
            fails.append(f"(6) {n} が {CTRL} と違う")
    extra = sorted(ours_rs - set(names))
    if extra:
        fails.append(f"(6) {CTRL} に無い記録がこの run にある: {extra}")
    probe = (run / "probe.yaml", ctrl / "probe.yaml")
    rec["probe.yaml (記録だけ)"] = (probe[0].read_bytes() == probe[1].read_bytes()) if all(p.is_file() for p in probe) else "片側に無い"
    return {"ok": not fails, "failures": fails, "files": rec, "pending": pending, "late": late}


# --- 前提の検査 (2): 格子 --------------------------------------------------------------------------------------------------
def compare_mesh_h5(ours: Path, ref: Path) -> dict:
    """MESH・CELLS・PLANES・VIZMESH の全データセットと属性、BCONDS/*/ の番号のデータセット (VALUE 以外)、VALUE/wall_dist のビット一致。"""
    import h5py
    import numpy as np

    def collect(f):
        req, rec = {}, {}

        def visit(name, obj):
            if not isinstance(obj, h5py.Dataset):
                return
            top = name.split("/")[0]
            if top in MESH_GROUPS or (top == "BCONDS" and "/VALUE/" not in name) or name == "VALUE/wall_dist":
                req[name] = obj
            elif top == "BCONDS":
                rec[name] = obj
        f.visititems(visit)
        return req, rec
    fails, n_eq, rec_diff = [], 0, []
    with h5py.File(ours, "r") as fa, h5py.File(ref, "r") as fb:
        ra, ca = collect(fa)
        rb, cb = collect(fb)
        if set(ra) != set(rb):
            fails.append(f"(2) 格子のデータセットの集合が違う: こちらだけ {sorted(set(ra) - set(rb))[:5]}、参照だけ {sorted(set(rb) - set(ra))[:5]}")
        for k in sorted(set(ra) & set(rb)):
            a, b = ra[k], rb[k]
            if a.dtype != b.dtype or a.shape != b.shape:
                fails.append(f"(2) {k}: 型・形が違う ({a.dtype}{a.shape} / {b.dtype}{b.shape})"); continue
            if np.asarray(a[()]).tobytes() != np.asarray(b[()]).tobytes():
                fails.append(f"(2) {k}: ビット一致しない"); continue
            n_eq += 1
        for g in MESH_GROUPS:
            if g in fa and g in fb:
                ga = {k: np.asarray(v).tolist() for k, v in fa[g].attrs.items()}; gb = {k: np.asarray(v).tolist() for k, v in fb[g].attrs.items()}
                if ga != gb:
                    fails.append(f"(2) {g} の属性が違う ({ga} / {gb})")
        for k in sorted(set(ca) & set(cb)):
            if ca[k].shape != cb[k].shape or np.asarray(ca[k][()]).tobytes() != np.asarray(cb[k][()]).tobytes():
                rec_diff.append(k)
        root = {k: str(v) for k, v in fa.attrs.items()}, {k: str(v) for k, v in fb.attrs.items()}
        coord_dtype = str(fa["MESH/COORD"].dtype)
    return {"ok": not fails, "failures": fails, "n_bit_equal": n_eq, "coord_dtype": coord_dtype,
            "bconds_values_differ (記録だけ)": rec_diff, "root_attrs (こちら, 参照; 記録だけ)": root}


# --- 前提の検査 (3): 格子の品質 --------------------------------------------------------------------------------------------
def quality_record(RA, run: Path, nj: int, scale: float) -> dict:
    """check_mesh_quality.py --ar-max 10000 を変換後の nozzle.h5 (node: VIZMESH の primal セル) に回し、prepare_ns の MESH_QUALITY.txt
    (cell 変換の写し、--ar-max は問題の ar_max) と合わせて PASS を確かめる。AR > 5000 のセルの数・位置 (x/r_t・壁からの層)・スキューは
    格子の番号 (i 方向の主、j = nj−1 が壁) から作った四角形で数える (#29 の quality.py と同じ数え方)。"""
    import h5py
    import numpy as np
    fails = []
    mq = (run / "MESH_QUALITY.txt").read_text() if (run / "MESH_QUALITY.txt").is_file() else ""
    mq_line = next((l.strip() for l in mq.splitlines()[::-1] if "VERDICT" in l), "")
    if not (mq_line.startswith("VERDICT: PASS") and f"AR<= {AR_MAX}" in mq_line):
        fails.append(f"(3) prepare_ns の MESH_QUALITY.txt が --ar-max {AR_MAX} の PASS でない ({mq_line!r})")
    q = subprocess.run([sys.executable, str(RA.FORGE_TOOLS / "check_mesh_quality.py"), str(run / "nozzle.h5"), "--ar-max", str(AR_MAX)],
                       capture_output=True, text=True)
    (run / "vgc_quality.txt").write_text(q.stdout + q.stderr)
    q_line = next((l.strip() for l in q.stdout.splitlines()[::-1] if "VERDICT" in l), "")
    if q.returncode != 0 or not q_line.startswith("VERDICT: PASS"):
        fails.append(f"(3) check_mesh_quality.py nozzle.h5 --ar-max {AR_MAX} が PASS でない (rc {q.returncode}, {q_line!r})")
    out = {"mesh_quality_txt": mq_line, "check_mesh_quality_nozzle_h5": q_line, "rc": q.returncode}
    try:
        if str(RA.FORGE_TOOLS) not in sys.path:
            sys.path.insert(0, str(RA.FORGE_TOOLS))
        import check_mesh_quality as CMQ
        with h5py.File(run / "nozzle.h5", "r") as f:
            c = np.asarray(f["MESH/COORD"], dtype=float).reshape(-1, 3)[:, :2]
            wd = np.asarray(f["VALUE/wall_dist"], dtype=float)
        ni = c.shape[0] // nj
        if ni * nj != c.shape[0]:
            raise ValueError(f"節点数 {c.shape[0]} が nj {nj} で割り切れない")
        idx = np.arange(ni * nj).reshape(ni, nj)
        qd = np.stack([idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]], axis=-1).reshape(-1, 4)
        ars, sk = CMQ.metrics_vectorized(c, qd.ravel(), np.arange(1, qd.shape[0] + 1) * 4, np.full(qd.shape[0], 9), False)
        A = ars.reshape(ni - 1, nj - 1); S = sk.reshape(ni - 1, nj - 1)
        X = c[:, 0].reshape(ni, nj) / scale
        wall_ok = bool(np.all(wd.reshape(ni, nj)[:, -1] <= wd.reshape(ni, nj)[:, -2]))   # 番号の向きの確認 (j = nj−1 が壁の側)
        hi = A > 5000
        ii, jj = np.nonzero(hi)
        out["ar_gt5000"] = {"n": int(hi.sum()), "wall_layers": [int((nj - 2) - jj.max()), int((nj - 2) - jj.min())] if hi.any() else None,
                            "x_rt": [float(X[ii, -1].min()), float(X[ii, -1].max())] if hi.any() else None,
                            "skew_max": float(S[hi].max()) if hi.any() else None, "wall_index_is_last_j": wall_ok}
        out["structured"] = {"ni": ni, "nj": nj, "ar_max": float(ars.max()), "skew_max": float(sk.max())}
        out["ref_29"] = {"n": 1325, "wall_layers": [0, 5], "x_rt": [-7.8084353921111935, -2.2321787081815407], "skew_max": 0.002435322543597667,
                         "ar_max": 8704.636122634876}
    except Exception as e:  # noqa: BLE001
        out["ar_gt5000_error"] = repr(e)                    # 記録の失敗は検査の不成立にしない (PASS は CLI で判定済み)
    out["ok"] = not fails
    out["failures"] = fails
    return out


# --- 観測だけの包み ---------------------------------------------------------------------------------------------------------
class Observer:
    """run_forge・stage_gate・_write_main_stage_started を包む (元の戻り値をそのまま返す)。"""

    def __init__(self, RA, run: Path, tags_by_sha: dict, prepared: dict, ctrl: Path):
        self.RA, self.run, self.tags, self.prep, self.ctrl = RA, run, tags_by_sha, prepared, ctrl
        self.orig = {"run_forge": RA.run_forge, "stage_gate": RA.stage_gate, "write_main": RA._write_main_stage_started}
        self.m9 = None

    def install(self):
        self.RA.run_forge, self.RA.stage_gate, self.RA._write_main_stage_started = self.run_forge, self.stage_gate, self.write_main

    def uninstall(self):
        self.RA.run_forge, self.RA.stage_gate, self.RA._write_main_stage_started = (self.orig["run_forge"], self.orig["stage_gate"],
                                                                                     self.orig["write_main"])

    def run_forge(self, run_dir):
        cfg_sha = sha256_file(Path(run_dir) / "solverConfig.yaml")
        tag = self.tags.get(cfg_sha, "?")
        t0 = time.time()
        jappend(self.run / "vgc_stages.jsonl", {"event": "start", "tag": tag, "t": t0, "config_sha256": cfg_sha})
        log(f"forge 起動 段 {tag} (設定 {cfg_sha[:16]})")
        rc = self.orig["run_forge"](run_dir)
        t1 = time.time()
        jappend(self.run / "vgc_stages.jsonl", {"event": "end", "tag": tag, "t0": t0, "t1": t1, "rc": rc, "config_sha256": cfg_sha})
        log(f"forge 終了 段 {tag} rc={rc} ({t1 - t0:.0f} s)")
        return rc

    def stage_gate(self, res_h5, cfg_text):
        probs = self.orig["stage_gate"](res_h5, cfg_text)
        rec = {"tag": self.tags.get(sha256_text(cfg_text), "?"), "res": Path(res_h5).name, "problems": list(probs), "t": time.time()}
        try:
            import h5py
            import numpy as np
            st = {}
            with h5py.File(res_h5, "r") as h:
                V = h["VALUE"]
                nf = {k: int(np.count_nonzero(~np.isfinite(np.asarray(V[k][()])))) for k in V if V[k].dtype.kind == "f"}
                st["nonfinite_total"] = int(sum(nf.values()))
                st["nonfinite_by_field"] = {k: v for k, v in nf.items() if v}
                for k in ("ro", "P", "T"):
                    if k in V:
                        a = np.asarray(V[k][()], dtype=float)
                        st[f"{k}_min"], st[f"{k}_max"] = float(np.nanmin(a)), float(np.nanmax(a))
                        st[f"{k}_n_nonpositive"] = int(np.count_nonzero(a <= 0))
                    else:
                        st[f"{k}_missing"] = True
            rec["stats"] = st
            rec["res_sha256"] = sha256_file(res_h5)
        except Exception as e:  # noqa: BLE001
            rec["stats_error"] = repr(e)
        jappend(self.run / "vgc_stage_gate.jsonl", rec)
        log(f"段終了ゲート {rec['tag']} {rec['res']}: 問題 {len(probs)} 件、{rec.get('stats', rec.get('stats_error'))}")
        return probs

    def write_main(self, run_dir, cfg_main, stages, done, cfg_pre):
        rec = self.orig["write_main"](run_dir, cfg_main, stages, done, cfg_pre)
        run_dir = Path(run_dir)
        fails = []
        if rec.get("res_files_present"):
            fails.append(f"通知の時点で前の段の出力が残っている: {rec['res_files_present']}")
        left = sorted(f.name for f in run_dir.glob("res_*"))
        if left:
            fails.append(f"通知の後に res_* がある: {left}")
        if rec["main_config"]["sha256"] != self.prep["main_sha256"]:
            fails.append("本段の設定の sha256 が準備の時点の solverConfig.yaml と違う")
        want_prev = [("S1_soft", 3000, self.prep["s1_sha256"]), ("S2_mid", 3000, self.prep["s2_sha256"])]
        got_prev = [(s.get("tag"), s.get("nsteps"), s.get("config_sha256")) for s in rec.get("previous_stages", [])]
        if got_prev != want_prev:
            fails.append(f"前の段の記録が導いた S1・S2 と違う: {got_prev}")
        if (rec.get("prestage_config") or {}).get("sha256") != self.prep["prestage_sha256"]:
            fails.append("前の段の元の設定の sha256 が準備の時点と違う")
        late = check_phys_files(run_dir, self.ctrl, late=True)
        fails += late["failures"]
        out = {"t": time.time(), "main_stage_started": rec, "files_6_late": late, "failures": fails, "ok": not fails,
               "initial_field_sha256": rec["initial_field"]["sha256"]}
        if fails:
            jdump(run_dir / "VGC_MAIN_START.json", out)
            (run_dir / "WATCH_ABORT").write_text(now() + "\n" + "\n".join(fails) + "\n")
            log("本段の直前の照合が不成立 — 本段の forge を起動させない: " + "; ".join(fails))
            raise RuntimeError("V-g0c-chain: 本段の直前の照合が不成立 (WATCH_ABORT)")
        cmd = [sys.executable, str(HERE / "m9_watch.py"), run_dir.name] + M9_ARGS
        with open(run_dir / "m9_watch.log", "a") as lf:
            self.m9 = subprocess.Popen(cmd, cwd=str(HERE), stdout=lf, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                       start_new_session=True)
        out.update(m9_cmd=cmd, m9_pid=self.m9.pid, m9_started=time.time(), m9_watch_sha256=sha256_file(HERE / "m9_watch.py"))
        jdump(run_dir / "VGC_MAIN_START.json", out)
        log(f"本段の通知を確かめ、見張りを起動 (pid {self.m9.pid}): {' '.join(cmd[1:])}")
        return rec


# --- 30 秒ごとの時刻・ディスク -------------------------------------------------------------------------------------------
class Timing(threading.Thread):
    """sp_arm.sh と同じ列の timing.tsv。run が無い間 (準備中) は溜めておき、run ができたら書く。ディスクの空き < 2 GiB で自分の forge を止める。"""

    def __init__(self, run: Path):
        super().__init__(daemon=True)
        self.rd, self.buf, self.stop_evt, self.disk_stopped = run, [], threading.Event(), False   # self.run は Thread.run と衝突するので使わない

    def last_step(self):
        p = self.rd / "residual_history.csv"
        try:
            with open(p, "rb") as f:
                f.seek(0, 2); n = f.tell(); f.seek(max(0, n - 8192))
                lines = [l for l in f.read().decode(errors="replace").splitlines() if l.strip()]
            s = lines[-1].split(",")[0].strip() if lines else ""
            return s if s.isdigit() else ""
        except OSError:
            return ""

    def run_once(self):
        t = time.time()
        rd = os.path.realpath(self.rd)
        ps = forge_pids()
        others = " ".join(c for _, c in ps if c != rd)
        row = f"{t:.6f}\t{self.last_step()}\t{len(ps)}\t{others}\n"
        if self.rd.is_dir():
            with open(self.rd / "timing.tsv", "a") as f:
                f.write("".join(self.buf) + row)
            self.buf = []
        else:
            self.buf.append(row)
        free = shutil.disk_usage(HERE).free
        if free < DISK_MIN and not self.disk_stopped:
            self.disk_stopped = True
            msg = f"{now()} ディスクの空き {free / 2**30:.2f} GiB < 2 GiB — cwd がこの run の forge を止める"
            log(msg)
            if self.rd.is_dir():
                (self.rd / "DISK_STOP").write_text(msg + "\n")
            ok = stop_own_forge(self.rd)
            log(f"自分の forge を止めた: {ok}")

    def run(self):  # noqa: D401
        while not self.stop_evt.is_set() and not (self.rd / "RUN_RC").exists():
            try:
                self.run_once()
            except Exception as e:  # noqa: BLE001
                log(f"timing の記録の失敗: {e!r}")
            self.stop_evt.wait(30)


def sh(cmd) -> str:
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout.strip()
    except (OSError, subprocess.SubprocessError) as e:
        return f"({cmd[0]} を実行できない: {e!r})"


def die(msg: str, code: int = 1):
    log("中止: " + msg)
    sys.exit(code)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--commit", required=True, help="~/forge-g0cchain の HEAD であるべき commit (vgc_prep.sh に渡したもの)")
    ap.add_argument("--wait-max", type=int, default=3600, help="ほかの forge が無くなるまで待つ上限 [s] (sp_arm.sh と同じ 1 時間)")
    ap.add_argument("--shared-gpu", action="store_true",
                    help="占有を待たずにほかのセッションの forge と GPU を共有して回す (2026-10-13 の登録の変更: 判定は時間を使わないので、"
                         "時間は共有の参考値として記録するだけ。使えるメモリ ≥ 4.5 GiB を確かめる)")
    ap.add_argument("--no-judge", action="store_true")
    ap.add_argument("--dry-prep", default=None, metavar="_vgc_dry_<名前>",
                    help="準備と前提の検査 (1)〜(6) だけを、run_0810 を使わずに別のディレクトリで行う (forge・見張り・判定は回さない。占有も待たない)")
    a = ap.parse_args(argv)
    if a.dry_prep is not None and not re.fullmatch(r"_vgc_dry_[A-Za-z0-9_.-]+", a.dry_prep):
        die("--dry-prep の名前は _vgc_dry_<英数字> (run_* と取り違えない)")
    run = HERE / (a.dry_prep or RUN)
    log(f"開始 {PLAN} run {run}")
    # --- 環境・コード・バイナリ ---
    env = {k: v for k, v in os.environ.items() if k.startswith("FORGE_")}
    if env != ENV_WANT:
        die(f"FORGE_* の環境変数が登録 (sp_arm.sh と同じ) と違う: {env} (vgc_run.sh から起動すること)")
    prep_done = HERE / "vgc_prep.done"
    if not prep_done.is_file() or a.commit not in prep_done.read_text():
        die(f"vgc_prep.done が無いか commit {a.commit} と合わない")
    if run.exists():
        die(f"{run} が既にある (既存の run は消さない)")
    sys.path.insert(0, str(WT / "design"))
    import numpy
    import forge_design
    from forge_design.evaluate import runner_axismach as RA
    fd = Path(forge_design.__file__).resolve()
    if not str(fd).startswith(str((WT / "design").resolve()) + os.sep):
        die(f"forge_design が作業ツリーの外から読まれた: {fd}")
    head = subprocess.run(["git", "-C", str(WT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(WT), "status", "--porcelain", "--untracked-files=no"], capture_output=True, text=True).stdout.strip()
    if len(a.commit) < 7 or not re.fullmatch(r"[0-9a-f]+", a.commit) or not head.startswith(a.commit):
        die(f"作業ツリーの HEAD {head} が --commit {a.commit} と違う")
    if dirty:
        die(f"作業ツリーに commit されていない変更がある: {dirty[:300]}")
    if numpy.__version__ != NUMPY_WANT:
        die(f"numpy {numpy.__version__} が登録の {NUMPY_WANT} と違う")
    bin_sha = sha256_file(BIN)
    conv = RA.converter_path()
    conv_sha = sha256_file(conv) if conv.is_file() else None
    if not bin_sha.startswith(BIN_SHA) or not (conv_sha or "").startswith(CONV_SHA):
        die(f"バイナリが登録と違う (forge {bin_sha[:16]}、変換器 {conv} {str(conv_sha)[:16]})")
    if (HERE / PROBLEM).read_bytes() != (WT / CASE / PROBLEM).read_bytes():
        die(f"{PROBLEM} が作業ツリーの commit の版と違う")
    code = {"forge_design": str(fd), "worktree": str(WT), "worktree_head": head, "forge_root": str(RA.FORGE_TOOLS.parents[1]),
            "python": sys.executable, "python_version": sys.version.split()[0], "numpy": numpy.__version__,
            "forge_bin": str(BIN), "forge_sha256": bin_sha, "converter": str(conv), "converter_sha256": conv_sha, "env": env,
            "tools_sha256": {n: sha256_file(RA.FORGE_TOOLS / n)[:16] for n in ("run_case.sh", "interp_field.py", "check_mesh_quality.py",
                                                                              "stage_manifest.py", "restart_field.py", "rerun_conditions.py",
                                                                              "check_convergence.py") if (RA.FORGE_TOOLS / n).is_file()},
            "case_scripts_sha256": {n: sha256_file(HERE / n)[:16] for n in ("vgc_run.py", "vgc_judge.py", "m9_watch.py", "cold_series.py",
                                                                           "cold_xcheck.py", "sp_judge.py", PROBLEM) if (HERE / n).is_file()}}
    log(f"コード {json.dumps(code, ensure_ascii=False)}")
    # --- 入力 ---
    ic = HERE / IC_RUN
    rs = sorted((f for f in ic.glob("res_*.h5") if re.fullmatch(r"res_\d+\.h5", f.name)), key=lambda f: int(re.findall(r"\d+", f.name)[0]))
    if not rs or rs[-1].name != IC_RES:
        die(f"ic_from {IC_RUN} の最後の res が {rs[-1].name if rs else None} ({IC_RES} であること)")
    ic_sha = sha256_file(rs[-1])
    if not ic_sha.startswith(IC_SHA):
        die(f"{IC_RUN}/{IC_RES} の sha256 {ic_sha[:16]} が登録 {IC_SHA} と違う")
    prep_sha = {n: sha256_file(HERE / PREP / n) for n in PREP_SHA}
    bad = {n: s[:16] for n, s in prep_sha.items() if not s.startswith(PREP_SHA[n])}
    if bad:
        die(f"{PREP} の sha256 が登録と違う: {bad}")
    ctrl = HERE / CTRL
    if not (HERE / B0_REF / "solverConfig.yaml").is_file() or not ctrl.is_dir():
        die(f"{B0_REF}/solverConfig.yaml か {CTRL} が無い")
    # --- GPU の占有: ほかの forge が無くなるまで待つ ---
    t_wait = time.time()
    if a.shared_gpu:
        others = [str(Path(f"/proc/{p}/cwd").resolve()) for p in forge_pids()]
        if others:
            log(f"ほかの forge と GPU を共有して回す (時間は共有の参考値): {others}")
        avail = int(next(l.split()[1] for l in Path("/proc/meminfo").read_text().splitlines() if l.startswith("MemAvailable"))) // 1024
        if avail < 4608:
            die(f"使えるメモリが 4.5 GiB 未満 ({avail} MiB)")
    else:
        while forge_pids() and time.time() - t_wait < a.wait_max and a.dry_prep is None:
            time.sleep(30)
        if forge_pids() and a.dry_prep is None:
            die(f"ほかの forge が {a.wait_max} s 空かない (占有できない): {forge_pids()}")
    free = shutil.disk_usage(HERE).free
    if free < DISK_MIN:
        die(f"ディスクの空き {free / 2**30:.2f} GiB < 2 GiB")
    gpu = sh(["nvidia-smi", "--query-gpu=name,memory.used,memory.total,utilization.gpu", "--format=csv,noheader"])
    mem = sh(["free", "-m"]).splitlines()
    log(f"占有の待ち {time.time() - t_wait:.0f} s、ディスクの空き {free / 2**30:.1f} GiB、GPU {gpu}、メモリ {mem[1] if len(mem) > 1 else mem}")
    # --- 準備 (総時間の起点) ---
    t_start = time.time()
    timing = Timing(run)
    if a.dry_prep is None:
        timing.start()
    pc = {"plan": PLAN, "created": now(), "run": RUN, "code": code, "t_start": t_start, "wait_s": t_start - t_wait,
          "inputs": {"ic_from": {"run": IC_RUN, "res": IC_RES, "sha256": ic_sha}, "prep_sha256": prep_sha,
                     "euler_ref": {"run": EULER_REF, "exists": (HERE / EULER_REF).is_dir(),
                                   "note": "prepare_ns は euler_ref のパスを prepare_info.json に書くだけで、中のファイルは読まない"}},
          "failures": [], "ok": False}
    log("prepare_ns を呼ぶ")
    try:
        info = RA.prepare_ns(HERE / PROBLEM, run, ic_from=ic, delta_r_csv=str(HERE / PREP / "delta_r_initial.csv"), offset="radial",
                             euler_ref=str(HERE / EULER_REF))
    except Exception as e:  # noqa: BLE001
        tb = traceback.format_exc()
        if isinstance(e, subprocess.CalledProcessError):
            tb += f"\n--- stdout ---\n{(e.stdout or '')[-3000:]}\n--- stderr ---\n{(e.stderr or '')[-3000:]}"
        log("prepare_ns が失敗:\n" + tb)
        pc["failures"].append(f"prepare_ns が失敗: {e!r}")
        if run.is_dir():
            (run / "START_EPOCH").write_text(f"{t_start:.6f}\n")
            jdump(run / "VGC_PRECHECK.json", pc)
            (run / "PRECHECK_FAILED").write_text(now() + "\nprepare_ns が失敗\n")
        timing.stop_evt.set()
        if not a.no_judge and run.is_dir() and a.dry_prep is None:
            subprocess.run([sys.executable, str(HERE / "vgc_judge.py")], cwd=str(HERE))
        return 2
    t_prep = time.time()
    (run / "START_EPOCH").write_text(f"{t_start:.6f}\n")
    pc["t_prepared"] = t_prep
    log(f"prepare_ns 終了 ({t_prep - t_start:.0f} s)")
    # --- 前提の検査 (1)〜(6) ---
    items = {}
    msh_sha = sha256_file(run / "nozzle.msh")
    items["1_msh"] = {"sha256": msh_sha, "want": MSH_SHA, "ok": msh_sha.startswith(MSH_SHA)}
    if not items["1_msh"]["ok"]:
        pc["failures"].append(f"(1) nozzle.msh の sha256 {msh_sha[:16]} が G0c {MSH_SHA} と違う")
    ref_h5 = ctrl / "nozzle.h5"
    if not (ref_h5.is_file() and sha256_file(ref_h5).startswith(CTRL_SHA)):
        ref_h5 = HERE / CTRL_INPUT
    if ref_h5.is_file() and sha256_file(ref_h5).startswith(CTRL_SHA):
        try:
            items["2_mesh"] = {"ref": str(ref_h5), **compare_mesh_h5(run / "nozzle.h5", ref_h5)}
        except Exception as e:  # noqa: BLE001
            items["2_mesh"] = {"ref": str(ref_h5), "ok": False, "failures": [f"(2) 照合の失敗: {e!r}"]}
    else:
        items["2_mesh"] = {"ok": False, "failures": [f"(2) {CTRL}/nozzle.h5 も {CTRL_INPUT} も sha256 {CTRL_SHA} で見つからない"]}
    pc["failures"] += items["2_mesh"]["failures"]
    items["3_quality"] = quality_record(RA, run, int(RA.load_problem(HERE / PROBLEM).mesh["nj"]), float(info["scale_m"]))
    pc["failures"] += items["3_quality"]["failures"]
    items["4_sha"] = {"ic_res": {"path": str(rs[-1]), "sha256": ic_sha, "registered": IC_SHA},
                      "nozzle_h5_after_prepare": sha256_file(run / "nozzle.h5"), "prepare_info_ic_from": info.get("ic_from"),
                      "main_initial_field": "本段の通知 (VGC_MAIN_START.json の initial_field_sha256) に記録"}
    main_txt = (run / "solverConfig.yaml").read_text()
    pre_p = run / RA.PRESTAGE_CONFIG
    if not pre_p.is_file():
        items["5_config"] = {"ok": False, "failures": [f"(5) {RA.PRESTAGE_CONFIG} が無い (main_line が効いていない)"]}
    else:
        bc_p = run / "bcondConfig.yaml"
        items["5_config"] = check_configs(RA, main_txt, pre_p.read_text(), (HERE / B0_REF / "solverConfig.yaml").read_text(),
                                          explicit_record(run), bc_p.read_text() if bc_p.is_file() else None)
        nsc = info.get("ns_stage_configs") or {}
        st = items["5_config"].get("stages", {})
        rec_vs = {k: nsc.get(k) == st.get(k, {}).get("sha256") for k in ("S1_soft", "S2_mid", "main")}
        rec_vs["prestage_config"] = nsc.get("prestage_config") == st.get("prestage_config", {}).get("sha256")
        items["5_config"]["prepare_info_sha_match"] = rec_vs
        if not all(rec_vs.values()):
            items["5_config"]["failures"].append(f"(5) prepare_info の ns_stage_configs が導いた段の設定と違う: {rec_vs}")
        ml = info.get("main_line") or {}
        if not (ml.get("cap") == 50.0 and ml.get("itj") == 5 and ml.get("cfl") == 4.0 and ml.get("relax") == RELAX
                and info.get("implicit_relax") == RELAX):
            items["5_config"]["failures"].append(f"(5) prepare_info の main_line が登録と違う: {ml}")
        items["5_config"]["ok"] = not items["5_config"]["failures"]
    pc["failures"] += items["5_config"]["failures"]
    items["6_files"] = check_phys_files(run, ctrl, late=False)
    pc["failures"] += items["6_files"]["failures"]
    # 記録だけ: 物理壁の表現・壁の表が生産の準備 (prod_confirm/prep) と同じか
    items["wall_record"] = {n: ((run / n).read_bytes() == (HERE / PREP / n).read_bytes()) if (run / n).is_file() and (HERE / PREP / n).is_file()
                            else "片側に無い" for n in ("wall_repr.json", "wall_physical.csv", "bcondConfig.yaml")}
    pc["items"] = items
    pc["ok"] = not pc["failures"]
    pc["t_checked"] = time.time()
    pc["dry_prep"] = a.dry_prep
    jdump(run / "VGC_PRECHECK.json", pc)
    if a.dry_prep is not None:
        log(f"--dry-prep {a.dry_prep}: 前提の検査 {'合格' if pc['ok'] else '不成立'} (forge は回さない): " + ("; ".join(pc["failures"]) or "なし"))
        return 0 if pc["ok"] else 2
    if not pc["ok"]:
        log("前提の検査が不成立 — forge を回さない (判定不能):\n  " + "\n  ".join(pc["failures"]))
        (run / "PRECHECK_FAILED").write_text(now() + "\n" + "\n".join(pc["failures"]) + "\n")
        timing.stop_evt.set()
        if not a.no_judge:
            subprocess.run([sys.executable, str(HERE / "vgc_judge.py")], cwd=str(HERE))
        return 2
    log(f"前提の検査 (1)〜(6) 合格: msh {msh_sha[:16]}、格子のビット一致 {items['2_mesh']['n_bit_equal']} データセット、"
        f"品質 {items['3_quality']['check_mesh_quality_nozzle_h5']}、AR>5000 {items['3_quality'].get('ar_gt5000')}、"
        f"(6) 未作成 {items['6_files']['pending']}")
    # --- 段階起動 full ---
    st5 = items["5_config"]["stages"]
    obs = Observer(RA, run, items["5_config"]["stage_tags_by_sha"],
                   {"main_sha256": st5["main"]["sha256"], "s1_sha256": st5["S1_soft"]["sha256"], "s2_sha256": st5["S2_mid"]["sha256"],
                    "prestage_sha256": st5["prestage_config"]["sha256"]}, ctrl)
    obs.install()
    run_end = {"t_start": t_start}
    try:
        rc = RA.run_staged_ns(run, stages="full")
        rc_txt = str(rc)
    except Exception as e:  # noqa: BLE001
        run_end["exception"] = traceback.format_exc()
        log("run_staged_ns が例外:\n" + run_end["exception"])
        rc_txt = "EXC"
    finally:
        obs.uninstall()
    t_end = time.time()
    (run / "RUN_RC").write_text(rc_txt + "\n")
    run_end.update(t_end=t_end, run_rc=rc_txt, main_started=(run / RA.MAIN_STAGE_STARTED).is_file())
    log(f"run_staged_ns 終了 RUN_RC={rc_txt} ({t_end - t_start:.0f} s)")
    if obs.m9 is not None:
        try:
            obs.m9.wait(timeout=3600)
        except subprocess.TimeoutExpired:
            run_end["m9_note"] = "見張りが RUN_RC の後 1 時間で終わらない"
            log(run_end["m9_note"])
    timing.stop_evt.set()
    timing.join(timeout=60)
    m9 = {}
    if (run / "m9_watch.json").is_file():
        d = json.loads((run / "m9_watch.json").read_text())
        m9 = {k: d.get(k) for k in ("status", "reach_step", "reach_sha256", "fail_step", "last_step", "note", "rc")}
    run_end["m9_watch"] = m9
    run_end["m9_returncode"] = None if obs.m9 is None else obs.m9.returncode
    jdump(run / "RUN_END.json", run_end)
    log(f"{RUN} 終了 RUN_RC={rc_txt} 見張り {m9.get('status')} {m9.get('reach_step')} (総 {t_end - t_start:.0f} s)")
    if not a.no_judge:
        j = subprocess.run([sys.executable, str(HERE / "vgc_judge.py")], cwd=str(HERE), capture_output=True, text=True)
        print(j.stdout + j.stderr, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
