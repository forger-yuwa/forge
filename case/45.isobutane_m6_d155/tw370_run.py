"""plan tooling-nozzle-isothermal-wall-chain §5.1 #32・§6「V-tw370 の事前登録」(commit b5df7416) と「登録の訂正」(commit b6aeb263、codex diagnose
`notes/reviews/2026-10-13-tw370-registration-diagnose.md`) の実行: 壁温 370 K・固定壁で G0c・B0 (上限 50) の停止の時点の出力を G0・上限 1 の参照と比べる。
腕を 1 つずつ別のプロセスで回す (Observer はモジュールの関数を差し替えるので、同じプロセスのスレッドで並列にしない; 訂正 6)。B で tw370_run.sh から起動する。

  腕     : R0  = G0・上限 1 (参照)   run_0811_tw370_g0_c1   problem_d155_ns_prod_coldmesh_tw370_g0_c1.yaml
           C1  = G0c・上限 1         run_0812_tw370_g0c_c1  problem_d155_ns_prod_coldmesh_tw370_g0c_c1.yaml
           C50 = G0c・上限 50 (B0)   run_0813_tw370_g0c_c50 problem_d155_ns_prod_coldmesh_tw370_g0c_c50.yaml
  コード : ~/forge-g0cchain (tw370_prep.sh がこの台本を commit した後の commit に合わせる) の design/ を sys.path の先頭にして forge_design を読む。
  forge  : ~/forge-fgeom7-fp64 (sha256 7e95e5bbd8f50c9e…)、変換器は同じビルドの convertGmshToForge (711f4438…)。環境変数は sp_arm.sh と同じ。
  ドナー : run_0754_t2_d2 の res_142500.h5 (G0・300 K・上限 1 の停止の出力)。最後の res がこれで、m9_watch.json の到達の step・sha と一致すること。
  準備   : R0 = prepare_ns(問題, run, ic_from なし, delta_r_csv=_band_ab/prod_confirm/prep/delta_r_initial.csv, offset="radial",
           euler_ref=run_0174_euler_v5d_M_r1) の後に、restart_field.py でドナーを nozzle.h5 へ同じ格子のビット一致で写す (訂正 3: 同じ格子への
           interp_field は使わない)。移した後の補正はしない (ω の床は prepare_ns の ic_from の経路だけが掛けるので R0 には掛からない; 掛けたら
           変わる節点の数を記録だけする)。
           C1 = 同じ prepare_ns に ic_from=run_0754_t2_d2 (チェーンの格子をまたぐ移し interp_field + ω の床)。移しの最近傍の距離と ω の床が効いた
           節点の数を記録する。
           R0・C1 とも run_staged_ns(stages="full") で S1 → S2 → 本段。
           C50 = C1 が本段を始める直前 (MAIN_STAGE_STARTED を書いた後・本段の forge の前) に C1 の台本が固定して保存した本段の入力
           (C1/C50_HANDOFF/: S2 の後の nozzle.h5 と準備のファイル、HANDOFF.json に sha) を待ち、その nozzle.h5 をバイトの写しで受け取る
           (訂正 4)。設定は C50 の問題 YAML から `_ns_configs` で作り、C1 の本段の設定との差が lineDtDirectionalCap (1 → 50) だけであること、
           前の段の元の設定が C1 と文字列で同じで、導いた S1・S2 の設定の sha が C1 が実際に回した S1・S2 と同じことを確かめ、
           run_staged_ns(stages="none") で本段だけを回す。C1 の S1・S2 の記録は C1_HISTORY/ に写す (C50 の履歴として引き継ぐ)。
  前提   : (1) nozzle.msh の sha256 (G0 9a891968fd9d0e36 / G0c 285f0fa4ecc0ed1c)。(2) nozzle.h5 の格子 (MESH・CELLS・PLANES・VIZMESH・BCONDS の番号)
           と VALUE/wall_dist が参照 (G0: run_0754_t2_d2/nozzle.h5 → _tr/q32_init.h5 ba859e02、G0c: run_0804_sp_g0c の入力 6f8fcda6d97482e3) と
           ビット一致。(3) check_mesh_quality.py --ar-max (G0 5000・G0c 10000) が PASS、AR > 5000 のセルを記録。(4) ドナー・移した後・本段の初期の場の
           sha256。(5) 本段の solverConfig.yaml が run_0752_t2_d1 (上限 1.0) と nStepOuter・outStepInterval・mesh.axisSegmentRWeight の明示
           (と C50 の lineDtDirectionalCap 50) 以外でキーごとに一致、S1・S2 に B0 のキーが無く (実効値で) implicitRelax 0.7。壁温は solverConfig に出ない。
           (6) bcondConfig.yaml が run_0804_sp_g0c と壁の行の Ts (300.0 → 370.0) だけ違い、species_meta.yaml・resolved_species_*.yaml が同じ
           (forge が起動時に書く記録は本段の直前に照合し直す)。(7) R0 の移しが保存量・乱流量・組成でドナーとビット一致。
           (8) C50 の本段の入力が C1 の本段の入力と格子・保存量・乱流量・組成でビット一致 (ファイルの sha256 も)。
           q300: tw370_judge.py --fix-q300 が起動の前に固定した _band_ab/cold_pair/tw370_q300.json があり、ドナーの sha が同じで、300 K の帯の感度が
           登録の値と合うこと (訂正 2 の a_q の出典)。
           資源 (訂正 6): nvidia-smi の GPU のメモリの空き (既存のジョブを含む) ≥ 2 GiB、ホストの MemAvailable ≥ 4.5 GiB、ディスクの空き ≥ 8 GiB を
           起動の前と本段の forge の前 (C50 は待った後) に確かめ、既存の GPU のジョブ・forge を記録する。どれかが外れたら forge を回さない (判定不能)。
  実行   : START_EPOCH と 30 秒ごとの timing.tsv (sp_arm.sh と同じ列)。MAIN_STAGE_STARTED を書く関数の直後に、通知の中身・前の段の出力が
           残っていないこと・(6) を照合し (C1 はさらに C50 への受け渡しを作り)、m9_watch.py を §6.11 の引数で起動する。
  見張り : ディスクの空きが 2 GiB を切ったら、cwd がこの run の forge だけを止める (DISK_STOP)。
  終わり : RUN_RC・RUN_END.json。3 本とも終わっていたら (最後に終わった腕が) tw370_judge.py を回す。

usage (B、tw370_run.sh 経由、3 本を別のプロセスで同時に):
  bash tw370_run.sh --arm R0 --commit <commit>;  bash tw370_run.sh --arm C1 --commit <commit>;  bash tw370_run.sh --arm C50 --commit <commit>
回す前の確かめ (run_08xx を使わない、forge を回さない):
  bash tw370_run.sh --arm R0 --commit <c> --dry-prep _tw370_dry_1_R0
  bash tw370_run.sh --arm C1 --commit <c> --dry-prep _tw370_dry_1_C1      (準備の後の nozzle.h5 で C50 への受け渡しの写しも作る)
  bash tw370_run.sh --arm C50 --commit <c> --dry-prep _tw370_dry_1_C50 --c1-dir _tw370_dry_1_C1
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
import tempfile
import threading
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent            # B: ~/forge-wallfit/case/45.isobutane_m6_d155 (tw370_prep.sh が作業ツリーから写す)
WT = Path.home() / "forge-g0cchain"
CASE = "case/45.isobutane_m6_d155"
PLAN = ("plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #32・§6 V-tw370 (事前登録 commit b5df7416、登録の訂正 commit b6aeb263)")
ARMS = {
    "R0": {"problem": "problem_d155_ns_prod_coldmesh_tw370_g0_c1.yaml", "run": "run_0811_tw370_g0_c1", "grid": "G0", "cap": 1.0,
           "msh_sha": "9a891968fd9d0e36", "ar_max": 5000, "mesh_ref": ("run_0754_t2_d2/nozzle.h5", "_tr/q32_init.h5", "ba859e02"),
           "transfer": "restart_field", "stages": "full"},
    "C1": {"problem": "problem_d155_ns_prod_coldmesh_tw370_g0c_c1.yaml", "run": "run_0812_tw370_g0c_c1", "grid": "G0c", "cap": 1.0,
           "msh_sha": "285f0fa4ecc0ed1c", "ar_max": 10000, "mesh_ref": ("run_0804_sp_g0c/nozzle.h5", "_g0c/q32_init.h5", "6f8fcda6d97482e3"),
           "transfer": "interp_field", "stages": "full"},
    "C50": {"problem": "problem_d155_ns_prod_coldmesh_tw370_g0c_c50.yaml", "run": "run_0813_tw370_g0c_c50", "grid": "G0c", "cap": 50.0,
            "msh_sha": "285f0fa4ecc0ed1c", "ar_max": 10000, "mesh_ref": ("run_0804_sp_g0c/nozzle.h5", "_g0c/q32_init.h5", "6f8fcda6d97482e3"),
            "transfer": "handoff", "stages": "none"},
}
DONOR, DONOR_RES, DONOR_STEP = "run_0754_t2_d2", "res_142500.h5", 142500
PREP = "_band_ab/prod_confirm/prep"
PREP_SHA = {"delta_r_initial.csv": "0bc9da5b6e8f0197", "wall_physical.csv": "31b4fd80818b4a0c"}
EULER_REF = "run_0174_euler_v5d_M_r1"
CFG_REF = "run_0752_t2_d1"            # 本段の設定の参照 (T2a、G0・上限 1.0・300 K)。壁温は solverConfig に出ない
PHYS_REF = "run_0804_sp_g0c"          # (6) の参照 (V-g0c-chain でチェーンの bcond・物性のファイルと一致を確かめた run)
PHYS_REC = DONOR                      # (6) の記録だけの照合先
TW_OLD, TW_NEW = "Ts: 300.0", "Ts: 370.0"
Q300 = "_band_ab/cold_pair/tw370_q300.json"
BIN = Path.home() / "forge-fgeom7-fp64/solver_density_cuda/build/forge"
BIN_SHA, CONV_SHA = "7e95e5bbd8f50c9e", "711f4438"
ENV_WANT = {"FORGE_CUDA_BLOCKSIZE": "128", "FORGE_OMEGA_BUDGET": "1", "FORGE_DIAG_COMMIT_LOSS": "500", "FORGE_BIN": str(BIN)}
NUMPY_WANT = "2.5.2"
M9_ARGS = ["--phase", "line", "--budget", "200000", "--consec", "2", "--interval", "2500", "--min-points", "9", "--keep-tail", "9"]
MAIN_STEPS, OUT_INT, RELAX = 200000, 2500, 0.7
DISK_MIN = 2 * 1024 ** 3                  # 走行中: これを切ったら自分の forge を止める
DISK_START_MIN = 8 * 1024 ** 3            # 起動の前: 3 本の出力 (各 res 10 枚前後) と受け渡しの写しの分
DISK_MAIN_MIN = 4 * 1024 ** 3             # 本段の forge の前: この腕の出力の分
GPU_FREE_MIN_MIB = 2048                   # 起動の前: GPU のメモリの空き (既存のジョブを含めた nvidia-smi の used から)
HOST_MIN_MIB = 4608                       # 起動の前: ホストの MemAvailable
PHYS_FILES = ("bcondConfig.yaml", "species_meta.yaml")       # + resolved_species_*.yaml
RS_RE = re.compile(r"^resolved_species_[0-9a-f]{16}(?:_[0-9a-f]{16})?\.yaml$")
MESH_GROUPS = ("MESH", "CELLS", "PLANES", "VIZMESH")
CONSERVED_RE = re.compile(r"^(ro|roUx|roUy|roUz|roe|roK|roOmega|roY\d+|roGamma|roReth|roXi|rog_.+|roQ\d_.+)$")
CONSERVED_REQ = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega")
# C1 → C50 の受け渡し (訂正 4): 本段の入力 (C50 の run に同じ名前で置く) と、C1 の S1・S2 の履歴 (C50 の C1_HISTORY/ に置く)
HANDOFF_DIR, HANDOFF_JSON = "C50_HANDOFF", "HANDOFF.json"
HANDOFF_REQ = ("nozzle.h5", "nozzle.msh", "bcondConfig.yaml", "probe.yaml", "species_meta.yaml", "prepare_info.json", "solverConfig_prestage.yaml")
HANDOFF_OPT = ("MESH_QUALITY.txt", "wall_repr.json", "wall_physical.csv", "wall_design.csv", "target_axis_M.csv")
HANDOFF_HIST = (("solverConfig.yaml", "solverConfig_C1_main.yaml"), ("stage_manifest.json", "stage_manifest_C1.json"),
                ("residual_history_S1_soft.csv", None), ("residual_history_S2_mid.csv", None), ("MAIN_STAGE_STARTED", "MAIN_STAGE_STARTED_C1"),
                ("TW370_PRECHECK.json", "TW370_PRECHECK_C1.json"), ("tw370_stages.jsonl", None), ("tw370_stage_gate.jsonl", None),
                ("restart_field.log", None))


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


LOG_TAG = "tw370_run"


def log(msg: str) -> None:
    print(f"[{LOG_TAG} {now()}] {msg}", flush=True)


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


def jload(p):
    try:
        return json.loads(Path(p).read_text())
    except (OSError, ValueError):
        return None


def res_steps(d: Path) -> list:
    """d の res_<N>.h5 を N の順に (prepare_ns の ic_from と同じ glob と並べ方)。"""
    return sorted((f for f in d.glob("res_[0-9]*.h5") if re.fullmatch(r"res_\d+\.h5", f.name)), key=lambda f: int(f.stem.split("_")[1]))


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


def sh(cmd) -> str:
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout.strip()
    except (OSError, subprocess.SubprocessError) as e:
        return f"({cmd[0]} を実行できない: {e!r})"


# --- 資源 (訂正 6) ---------------------------------------------------------------------------------------------------------
def resources(where: str, disk_min: int = DISK_START_MIN) -> dict:
    """GPU のメモリ (nvidia-smi、既存のジョブを含む)・ホストの MemAvailable・ディスクの空きを確かめる。読めないのも不成立。
    disk_min は起動の前 8 GiB (3 本の出力の分)、本段の forge の前 4 GiB (この腕の出力の分; ほかの腕が既に書いた分を二重に数えない)。"""
    fails, rec = [], {"where": where, "t": time.time(), "disk_min_gib": disk_min / 2 ** 30}
    try:
        q = subprocess.run(["nvidia-smi", "--query-gpu=index,name,memory.used,memory.total", "--format=csv,noheader,nounits"],
                           capture_output=True, text=True, timeout=60)
        rows = [[c.strip() for c in l.split(",")] for l in q.stdout.strip().splitlines() if l.strip()]
        if q.returncode != 0 or not rows:
            raise RuntimeError(f"rc {q.returncode}: {q.stderr.strip()[:200]}")
        used, total = float(rows[0][2]), float(rows[0][3])
        rec["gpu"] = {"index": rows[0][0], "name": rows[0][1], "used_mib": used, "total_mib": total, "free_mib": total - used, "n_gpu": len(rows)}
        apps = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,process_name,used_memory", "--format=csv,noheader,nounits"],
                              capture_output=True, text=True, timeout=60).stdout.strip().splitlines()
        rec["gpu_apps"] = []
        for l in apps:
            c = [x.strip() for x in l.split(",")]
            if len(c) >= 3:
                try:
                    cwd = os.path.realpath(f"/proc/{c[0]}/cwd")
                except OSError:
                    cwd = "?"
                rec["gpu_apps"].append({"pid": c[0], "name": c[1], "used_mib": c[2], "cwd": cwd})
        if total - used < GPU_FREE_MIN_MIB:
            fails.append(f"GPU のメモリの空き {total - used:.0f} MiB < {GPU_FREE_MIN_MIB} MiB (既存のジョブ {len(rec['gpu_apps'])} 本)")
    except Exception as e:  # noqa: BLE001
        fails.append(f"nvidia-smi で GPU のメモリを読めない ({e!r})")
    try:
        avail = int(next(l.split()[1] for l in Path("/proc/meminfo").read_text().splitlines() if l.startswith("MemAvailable"))) // 1024
        rec["host_mem_available_mib"] = avail
        if avail < HOST_MIN_MIB:
            fails.append(f"ホストの MemAvailable {avail} MiB < {HOST_MIN_MIB} MiB")
    except Exception as e:  # noqa: BLE001
        fails.append(f"/proc/meminfo を読めない ({e!r})")
    free = shutil.disk_usage(HERE).free
    rec["disk_free_gib"] = free / 2 ** 30
    if free < disk_min:
        fails.append(f"ディスクの空き {free / 2**30:.2f} GiB < {disk_min / 2**30:.0f} GiB")
    rec["forge"] = [{"pid": p, "cwd": c} for p, c in forge_pids()]
    rec["ok"], rec["failures"] = not fails, fails
    return rec


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
def check_configs(RA, main: str, pre: str, ref: str, cap: float, explicit: dict | None = None, bcond: str | None = None) -> dict:
    """本段の設定 main・前の段の元の設定 pre・参照 ref (run_0752 の solverConfig.yaml、上限 1.0) を照合する。
    本段と ref の差として許すのは nStepOuter (こちら 200000)・outStepInterval (同じ 2500 のはず)・mesh.axisSegmentRWeight (片側だけが明示の 1)
    と、cap ≠ 1 の腕の lineDtDirectionalCap (1.0 → cap)。S1・S2 は RA.ns_stage_configs(main, "full", cfg_pre=pre) で導き (C50 は回さないが
    C1 が回した S1・S2 との照合に使う)、B0 だけのキー (起動する設定の実効値で見る: 書いていないときの値 B0_OFF の明示は無いのと同じ、b0_left) と本段の
    extraFields が無く implicitRelax 0.7 であること。
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
                "diff_vs_ref": {}, "stages": {}, "explicit_added": ex}
    dm, dr = RA._cfg_load(main_b), RA._cfg_load(ref)
    diff = RA._cfg_diff_paths(dr, dm)
    shown = {".".join(map(str, k)): ["<無し>" if a is M else a, "<無し>" if b is M else b] for k, (a, b) in diff.items()}
    for k, (a, b) in diff.items():
        if k in (("time", "last", "nStepOuter"), ("time", "outStepInterval")):
            continue
        if k == ("mesh", "axisSegmentRWeight") and {type(a), type(b)} == {int, type(M)} and 1 in (a, b):
            continue                                   # ソルバの既定 1 を片側だけが明示 (solverConfig.hpp:454)
        if k == ("time", "deltaT", "lineDtDirectionalCap") and cap != 1.0 and a == 1.0 and b == cap:
            continue
        fails.append(f"(5) 本段の設定が {CFG_REF} と違う: {'.'.join(map(str, k))} = {shown['.'.join(map(str, k))]} (参照, こちら)")
    for pth, want in ((("time", "last", "nStepOuter"), MAIN_STEPS), (("time", "outStepInterval"), OUT_INT),
                      (("time", "deltaT", "lineDtDirectionalCap"), cap), (("time", "deltaT", "implicitRelax"), RELAX)):
        if RA._cfg_path_get(dm, pth) != want:
            fails.append(f"(5) 本段の {'.'.join(pth)} = {RA._cfg_path_get(dm, pth)!r} ({want} であること)")
    # 壁温が solverConfig に効く所があれば 300 K の参照 (run_0752) との差に出る (上の許した差だけなら壁温は bcond にだけ出ている)
    b0_only = sorted(p for p in RA._MAIN_LINE_ONLY_PATHS if p != ("output",))
    stages = {}
    try:
        pre_st, main_f, main_want = RA.ns_stage_configs(main, "full", cfg_pre=pre)
        pre_st_b = RA.ns_stage_configs(main_b, "full", cfg_pre=pre_b)[0] if ex else pre_st     # 足す前の段 (照合はこちら)
        differ = [t for (t, c, *_), (_, cb, *_) in zip(pre_st, pre_st_b) if explicit_remove(c, ex.get("prestage")) != cb]
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "failures": fails + [f"(5) ns_stage_configs が失敗: {e!r}"], "diff_vs_ref": shown, "stages": {}, "explicit_added": ex}
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
    return {"ok": not fails, "failures": fails, "diff_vs_ref": shown, "stages": stages,
            "stage_tags_by_sha": {v["sha256"]: k for k, v in stages.items() if k != "prestage_config"}, "explicit_added": ex}


# --- 前提の検査 (6): bcond・物性のファイル --------------------------------------------------------------------------------
def bcond_vs_ref(ours: str, ref: str) -> tuple:
    """bcond が参照 (300 K) と壁の行の Ts (300.0 → 370.0) だけ違うか。戻り (ok, 説明)。"""
    if ours.count(TW_NEW) != 1 or ref.count(TW_OLD) != 1:
        return False, f"{TW_NEW!r} がこちらに {ours.count(TW_NEW)} 個、{TW_OLD!r} が参照に {ref.count(TW_OLD)} 個 (1 個ずつであること)"
    line = next(l for l in ours.splitlines() if TW_NEW in l)
    if not line.startswith("wall:") or "kind: wall_isothermal" not in line:
        return False, f"{TW_NEW!r} が壁の等温の行に無い ({line!r})"
    if ours.replace(TW_NEW, TW_OLD) != ref:
        return False, "壁の行の Ts を戻しても参照と一致しない"
    return True, line


def check_phys_files(run: Path, ref: Path, late: bool) -> dict:
    """bcondConfig.yaml (壁の行の Ts だけ違う)・species_meta.yaml・resolved_species_*.yaml (同じ) を参照 run と照合する。late=False (準備の直後)
    は、forge が起動時に書く記録で run にまだ無いものを「未作成」として残す (後で照合し直す)。late=True (本段の直前・判定) は全部がそろうこと。"""
    names = sorted(set(PHYS_FILES) | {f.name for f in ref.iterdir() if RS_RE.match(f.name)})
    ours_rs = {f.name for f in run.iterdir() if RS_RE.match(f.name)}
    fails, rec, pending = [], {}, []
    for n in names:
        a, b = run / n, ref / n
        if not b.is_file():
            fails.append(f"(6) {ref.name}/{n} が無い"); rec[n] = "参照に無い"; continue
        if not a.is_file():
            if late or n in PHYS_FILES:
                fails.append(f"(6) {n} がこの run に無い"); rec[n] = "無い"
            else:
                pending.append(n); rec[n] = "未作成 (forge が起動時に書く記録; 本段の直前に照合し直す)"
            continue
        if n == "bcondConfig.yaml":
            ok, why = bcond_vs_ref(a.read_text(), b.read_text())
            rec[n] = {"same_except_wall_Ts": ok, "detail": why, "sha256": sha256_file(a)[:16], "ref_sha256": sha256_file(b)[:16]}
            if not ok:
                fails.append(f"(6) bcondConfig.yaml が {ref.name} と壁の行の Ts 以外でも違う: {why}")
            continue
        same = a.read_bytes() == b.read_bytes()
        rec[n] = {"same": same, "sha256": sha256_file(a)[:16], "ref_sha256": sha256_file(b)[:16]}
        if not same:
            fails.append(f"(6) {n} が {ref.name} と違う")
    extra = sorted(ours_rs - set(names))
    if extra:
        fails.append(f"(6) {ref.name} に無い記録がこの run にある: {extra}")
    probe = (run / "probe.yaml", ref / "probe.yaml")
    rec["probe.yaml (記録だけ)"] = (probe[0].read_bytes() == probe[1].read_bytes()) if all(p.is_file() for p in probe) else "片側に無い"
    return {"ok": not fails, "failures": fails, "files": rec, "pending": pending, "late": late, "ref": ref.name}


def phys_record_vs_donor(run: Path, donor: Path) -> dict:
    """記録だけ: ドナー (run_0754) の bcond・物性のファイルとの照合。"""
    out = {}
    for n in sorted(set(PHYS_FILES) | {f.name for f in donor.iterdir() if RS_RE.match(f.name)} if donor.is_dir() else set(PHYS_FILES)):
        a, b = run / n, donor / n
        if not (a.is_file() and b.is_file()):
            out[n] = "片側に無い"
        elif n == "bcondConfig.yaml":
            out[n] = bcond_vs_ref(a.read_text(), b.read_text())[1]
        else:
            out[n] = a.read_bytes() == b.read_bytes()
    return out


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


def mesh_ref_path(arm: dict):
    """(2) の参照のファイル (参照 run の nozzle.h5、無いか sha が違えば共通の入力) と、見つかったかの説明。"""
    a, b, want = arm["mesh_ref"]
    for c in (HERE / a, HERE / b):
        if c.is_file() and sha256_file(c).startswith(want):
            return c, f"{c.relative_to(HERE)} → {c.resolve()} (sha256 {want}…)"
    return None, f"(2) {a} も {b} も sha256 {want} で見つからない"


# --- 前提の検査 (3): 格子の品質 --------------------------------------------------------------------------------------------
def quality_record(RA, run: Path, nj: int, scale: float, ar_max: int) -> dict:
    """check_mesh_quality.py --ar-max <ar_max> を変換後の nozzle.h5 に回し、prepare_ns の MESH_QUALITY.txt (cell 変換の写し) と合わせて PASS を
    確かめる。AR > 5000 のセルの数・位置 (x/r_t・壁からの層)・スキューは格子の番号 (i 方向の主、j = nj−1 が壁) から作った四角形で数える。"""
    import h5py
    import numpy as np
    fails = []
    mq = (run / "MESH_QUALITY.txt").read_text() if (run / "MESH_QUALITY.txt").is_file() else ""
    mq_line = next((l.strip() for l in mq.splitlines()[::-1] if "VERDICT" in l), "")
    if not (mq_line.startswith("VERDICT: PASS") and f"AR<= {ar_max}" in mq_line):
        fails.append(f"(3) prepare_ns の MESH_QUALITY.txt が --ar-max {ar_max} の PASS でない ({mq_line!r})")
    q = subprocess.run([sys.executable, str(RA.FORGE_TOOLS / "check_mesh_quality.py"), str(run / "nozzle.h5"), "--ar-max", str(ar_max)],
                       capture_output=True, text=True)
    (run / "tw370_quality.txt").write_text(q.stdout + q.stderr)
    q_line = next((l.strip() for l in q.stdout.splitlines()[::-1] if "VERDICT" in l), "")
    if q.returncode != 0 or not q_line.startswith("VERDICT: PASS"):
        fails.append(f"(3) check_mesh_quality.py nozzle.h5 --ar-max {ar_max} が PASS でない (rc {q.returncode}, {q_line!r})")
    out = {"ar_max": ar_max, "mesh_quality_txt": mq_line, "check_mesh_quality_nozzle_h5": q_line, "rc": q.returncode}
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
        wall_ok = bool(np.all(wd.reshape(ni, nj)[:, -1] <= wd.reshape(ni, nj)[:, -2]))
        hi = A > 5000
        ii, jj = np.nonzero(hi)
        out["ar_gt5000"] = {"n": int(hi.sum()), "wall_layers": [int((nj - 2) - jj.max()), int((nj - 2) - jj.min())] if hi.any() else None,
                            "x_rt": [float(X[ii, -1].min()), float(X[ii, -1].max())] if hi.any() else None,
                            "skew_max": float(S[hi].max()) if hi.any() else None, "wall_index_is_last_j": wall_ok}
        out["structured"] = {"ni": ni, "nj": nj, "ar_max": float(ars.max()), "skew_max": float(sk.max())}
    except Exception as e:  # noqa: BLE001
        out["ar_gt5000_error"] = repr(e)                    # 記録の失敗は検査の不成立にしない (PASS は CLI で判定済み)
    out["ok"] = not fails
    out["failures"] = fails
    return out


# --- 移し (訂正 3) -----------------------------------------------------------------------------------------------------------
def omega_floor_count(p, h5path: Path, mode: str) -> dict:
    """prepare_ns の ω の床 (6 ν/(0.075 d²)、T は CPG の式) を同じ式で作り、mode="would_raise" はそれより小さい節点の数 (R0: 掛けていない床が
    掛かったら変わる数)、mode="active" は床の値と一致する節点の数 (C1: 床が効いた数) を返す。記録だけ。"""
    import h5py
    import numpy as np
    sys.path.insert(0, str(WT / "design"))
    from forge_design.feedback.deltastar import _sutherland
    g = p.gamma
    R_gas = p.cp * (g - 1.0) / g
    with h5py.File(h5path, "r") as f:
        ro = f["/VALUE/ro"][:]; roUx, roUy = f["/VALUE/roUx"][:], f["/VALUE/roUy"][:]
        roe = f["/VALUE/roe"][:]; wd = f["/VALUE/wall_dist"][:]; rw = f["/VALUE/roOmega"][:]
    T = np.maximum((roe - 0.5 * (roUx ** 2 + roUy ** 2) / np.maximum(ro, 1e-12)) * (g - 1.0) / (np.maximum(ro, 1e-12) * R_gas), 50.0)
    nu = _sutherland(T) / np.maximum(ro, 1e-12)
    fl = ro * (6.0 * nu / (0.075 * np.maximum(wd, 1e-9) ** 2))
    if mode == "would_raise":
        n = int(np.count_nonzero(rw < fl))
        return {"mode": "掛けていない床が掛かったら変わる節点の数 (記録だけ)", "n": n, "n_nodes": int(ro.size)}
    n = int(np.count_nonzero(np.isclose(rw, fl, rtol=1e-13, atol=0.0)))
    return {"mode": "床の値と一致する節点の数 (床が効いた節点、記録だけ)", "n": n, "n_nodes": int(ro.size)}


def transfer_restart(RA, run: Path, donor_res: Path) -> dict:
    """R0: restart_field.py --keep-src-dtype でドナーを nozzle.h5 へ (同じ格子・保存量の index コピー)。写した後に、宛先の VALUE (wall_dist 以外)
    がドナーと型も含めてビット一致し、ドナーの保存量・乱流量・組成の全部が宛先にあることを独立に確かめる。
    --keep-src-dtype の理由: prepare_ns の初期の場 (paste_isentropic_ic) は roY0・roY1 を float32 のデータセットで作る (変換器の出力に無いので
    新しく作る、ic.py)。付けないと倍精度のドナーの roY* が float32 に丸められ、ビット一致にならない (restart_field の「型の縮小で丸めあり」)。
    付けると宛先の roY* は倍精度のデータセットに作り直される (変わった型は dtype_before/after に記録)。"""
    import h5py
    import numpy as np
    fails = []
    before = sha256_file(run / "nozzle.h5")
    with h5py.File(run / "nozzle.h5", "r") as d:
        dt_before = {k: str(d["VALUE"][k].dtype) for k in d["VALUE"]}
    cmd = [sys.executable, str(RA.FORGE_TOOLS / "restart_field.py"), str(donor_res), str(run / "nozzle.h5"), "--dst-run", str(run),
           "--forge", str(BIN), "--keep-src-dtype"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    (run / "tw370_transfer.log").write_text(" ".join(cmd) + "\n" + r.stdout + r.stderr)
    vline = next((l.strip() for l in r.stdout.splitlines()[::-1] if l.startswith("VERDICT")), "")
    if r.returncode != 0 or not vline.startswith("VERDICT: OK") or "SRC とビット一致" not in vline:
        fails.append(f"(7) restart_field.py が失敗かビット一致でない (rc {r.returncode}, {vline!r}, {(r.stdout + r.stderr).strip()[-300:]!r})")
    rec = {"tool": "restart_field.py", "cmd": cmd, "rc": r.returncode, "verdict": vline, "nozzle_sha256_before": before}
    if not fails:
        with h5py.File(donor_res, "r") as s, h5py.File(run / "nozzle.h5", "r") as d:
            sv, dv = s["VALUE"], d["VALUE"]
            src_cons = sorted(k for k in sv if CONSERVED_RE.match(k))
            moved, bad = [], []
            for k in sorted(dv):
                if k == "wall_dist":
                    continue
                if k not in sv:
                    bad.append(f"{k} (ドナーに無い)"); continue
                a, b = np.asarray(sv[k][()]), np.asarray(dv[k][()])
                if a.dtype != b.dtype or a.tobytes() != b.tobytes():
                    bad.append(f"{k} ({a.dtype}/{b.dtype})")
                moved.append(k)
            lack = [k for k in src_cons if k not in dv] + [k for k in CONSERVED_REQ if k not in dv]
            if bad:
                fails.append(f"(7) 移した後の nozzle.h5 がドナーとビット一致しない: {bad}")
            if lack:
                fails.append(f"(7) ドナーの保存量・乱流量・組成が宛先に無い (写っていない): {sorted(set(lack))}")
            attrs = {k: str(v) for k, v in d.attrs.items()}
            dt_after = {k: str(dv[k].dtype) for k in dv}
        rec.update(moved=moved, donor_conserved=src_cons, nozzle_attrs=attrs,
                   dtype_changed={k: [dt_before.get(k), v] for k, v in dt_after.items() if dt_before.get(k) != v})
    rec["nozzle_sha256_after"] = sha256_file(run / "nozzle.h5")
    rec["correction_after_transfer"] = "なし (ω の床は掛けない; prepare_ns は ic_from の経路だけで床を掛ける)"
    rec["ok"], rec["failures"] = not fails, fails
    return rec


def nearest_record(donor_res: Path, nozzle: Path, scale: float) -> dict:
    """C1: interp_field の最近傍の対応を独立に作り直した記録 (interp_field の出力は prepare_ns が捨てる)。interp_field はドナーの格子を
    res_*.h5 の MESH/COORD から読む (run_0754 の nozzle.h5 → _tr/q32_init.h5 は読まない)。"""
    import h5py
    import numpy as np
    from scipy.spatial import cKDTree
    with h5py.File(donor_res, "r") as s, h5py.File(nozzle, "r") as d:
        cs = np.asarray(s["MESH/COORD"]).reshape(-1, 3)[:, :2]
        cd = np.asarray(d["MESH/COORD"]).reshape(-1, 3)[:, :2]
        dts = (str(s["MESH/COORD"].dtype), str(d["MESH/COORD"].dtype))
    dist, _ = cKDTree(cs).query(cd)
    x = cd[:, 0] / scale
    out_reg = (x < -8.5) | (x > -1.5)
    return {"donor_coord_source": "res_142500.h5 の MESH/COORD (interp_field.centroids)", "coord_dtype (ドナー, 宛先)": dts,
            "n_dst": int(cd.shape[0]), "n_src": int(cs.shape[0]), "max_dist_rt": float(dist.max() / scale),
            "n_zero_dist": int(np.count_nonzero(dist == 0)),
            "outside_respace_n": int(out_reg.sum()), "outside_respace_n_zero_dist": int(np.count_nonzero(dist[out_reg] == 0)),
            "outside_respace_max_dist_rt": float(dist[out_reg].max() / scale) if out_reg.any() else None}


def start_identity(a: Path, b: Path) -> dict:
    """C50 の本段の入力 a と C1 の本段の入力 b の照合 (訂正 4): ファイルの sha256 と、全データセット (格子・保存量・乱流量・組成・wall_dist) と
    属性のビット一致。"""
    import h5py
    import numpy as np
    fails = []
    sa, sb = sha256_file(a), sha256_file(b)
    if sa != sb:
        fails.append(f"(8) 本段の入力の sha256 が C1 と違う ({sa[:16]} / {sb[:16]})")
    names = []
    with h5py.File(a, "r") as fa, h5py.File(b, "r") as fb:
        da, db = {}, {}
        fa.visititems(lambda n, o: da.__setitem__(n, o) if isinstance(o, h5py.Dataset) else None)
        fb.visititems(lambda n, o: db.__setitem__(n, o) if isinstance(o, h5py.Dataset) else None)
        if set(da) != set(db):
            fails.append(f"(8) データセットの集合が違う ({sorted(set(da) ^ set(db))[:5]})")
        for k in sorted(set(da) & set(db)):
            x, y = np.asarray(da[k][()]), np.asarray(db[k][()])
            if x.dtype != y.dtype or x.shape != y.shape or x.tobytes() != y.tobytes():
                fails.append(f"(8) {k} がビット一致しない")
            names.append(k)
        for g in ("/", "VALUE"):
            if {k: str(v) for k, v in fa[g].attrs.items()} != {k: str(v) for k, v in fb[g].attrs.items()}:
                fails.append(f"(8) {g} の属性が違う")
        value = sorted(k for k in fa["VALUE"])
    req = [k for k in CONSERVED_REQ if f"VALUE/{k}" not in names] + ([] if any(re.fullmatch(r"VALUE/roY\d+", n) for n in names) else ["roY*"])
    if req:
        fails.append(f"(8) 保存量・乱流量・組成が入力に無い: {req}")
    return {"ok": not fails, "failures": fails, "sha256": sa, "c1_sha256": sb, "n_datasets": len(names), "value_datasets": value}


# --- C1 → C50 の受け渡し (訂正 4) -----------------------------------------------------------------------------------------
def make_handoff(run_dir: Path, mss: dict, dry: bool = False) -> dict:
    """C1 の本段の入力 (MAIN_STAGE_STARTED の後・本段の forge の前の nozzle.h5 と準備のファイル) と S1・S2 の履歴を run_dir/C50_HANDOFF/ に写し、
    sha を確かめてから HANDOFF.json を (最後に・原子的に) 書く。C50 はこの JSON が出るまで待つ。"""
    hd = run_dir / HANDOFF_DIR
    hd.mkdir(exist_ok=False)
    inputs, hist, miss = {}, {}, []
    rs = sorted(f.name for f in run_dir.iterdir() if RS_RE.match(f.name))
    for n in HANDOFF_REQ + HANDOFF_OPT + tuple(rs):
        if (run_dir / n).is_file():
            shutil.copy2(run_dir / n, hd / n)
            inputs[n] = sha256_file(hd / n)
            if sha256_file(run_dir / n) != inputs[n]:
                raise RuntimeError(f"受け渡しの写し {n} の sha256 が元と違う")
        elif n in HANDOFF_REQ:
            miss.append(n)
    if miss:
        raise RuntimeError(f"受け渡しに要るファイルが無い: {miss}")
    for n, dst in HANDOFF_HIST:
        if (run_dir / n).is_file():
            shutil.copy2(run_dir / n, hd / (dst or n))
            hist[dst or n] = sha256_file(hd / (dst or n))
    init_sha = (mss.get("initial_field") or {}).get("sha256") if mss else None
    if not dry and inputs["nozzle.h5"] != init_sha:
        raise RuntimeError(f"受け渡しの nozzle.h5 の sha256 {inputs['nozzle.h5'][:16]} が MAIN_STAGE_STARTED の本段の初期の場 {str(init_sha)[:16]} と違う")
    rec = {"plan": PLAN, "created": now(), "t": time.time(), "dry": dry, "from_run": run_dir.name, "inputs": inputs, "history": hist,
           "resolved_species": rs, "initial_field_sha256": inputs["nozzle.h5"], "c1_main_stage_started": mss,
           "c1_main_config_sha256": sha256_file(run_dir / "solverConfig.yaml"),
           "c1_prestage_config_sha256": sha256_file(run_dir / "solverConfig_prestage.yaml")}
    jdump(hd / HANDOFF_JSON, rec)
    return rec


# --- 観測だけの包み ---------------------------------------------------------------------------------------------------------
class Observer:
    """run_forge・stage_gate・_write_main_stage_started を包む (元の戻り値をそのまま返す)。C1 は本段の直前に C50 への受け渡しを作る。"""

    def __init__(self, RA, run: Path, tags_by_sha: dict, prepared: dict, ref: Path, arm: str):
        self.RA, self.run, self.tags, self.prep, self.ref, self.arm = RA, run, tags_by_sha, prepared, ref, arm
        self.orig = {"run_forge": RA.run_forge, "stage_gate": RA.stage_gate, "write_main": RA._write_main_stage_started}
        self.m9 = None
        self.handoff = None

    def install(self):
        self.RA.run_forge, self.RA.stage_gate, self.RA._write_main_stage_started = self.run_forge, self.stage_gate, self.write_main

    def uninstall(self):
        self.RA.run_forge, self.RA.stage_gate, self.RA._write_main_stage_started = (self.orig["run_forge"], self.orig["stage_gate"],
                                                                                     self.orig["write_main"])

    def run_forge(self, run_dir):
        cfg_sha = sha256_file(Path(run_dir) / "solverConfig.yaml")
        tag = self.tags.get(cfg_sha, "?")
        t0 = time.time()
        jappend(self.run / "tw370_stages.jsonl", {"event": "start", "tag": tag, "t": t0, "config_sha256": cfg_sha})
        log(f"forge 起動 段 {tag} (設定 {cfg_sha[:16]})")
        rc = self.orig["run_forge"](run_dir)
        t1 = time.time()
        jappend(self.run / "tw370_stages.jsonl", {"event": "end", "tag": tag, "t0": t0, "t1": t1, "rc": rc, "config_sha256": cfg_sha})
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
        jappend(self.run / "tw370_stage_gate.jsonl", rec)
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
        got_prev = [(s.get("tag"), s.get("nsteps"), s.get("config_sha256")) for s in rec.get("previous_stages", [])]
        if got_prev != self.prep["want_prev"]:
            fails.append(f"前の段の記録が導いた S1・S2 と違う: {got_prev}")
        if (rec.get("prestage_config") or {}).get("sha256") != self.prep["prestage_sha256"]:
            fails.append("前の段の元の設定の sha256 が準備の時点と違う")
        if self.prep.get("initial_field_sha256") and rec["initial_field"]["sha256"] != self.prep["initial_field_sha256"]:
            fails.append("本段の初期の場の sha256 が C1 の本段の入力と違う (C50)")
        late = check_phys_files(run_dir, self.ref, late=True)
        fails += late["failures"]
        res = resources("本段の forge の前", DISK_MAIN_MIN)
        if not res["ok"]:
            fails += ["資源: " + f for f in res["failures"]]
        out = {"t": time.time(), "arm": self.arm, "main_stage_started": rec, "files_6_late": late, "resources": res, "failures": fails,
               "ok": not fails, "initial_field_sha256": rec["initial_field"]["sha256"]}
        if fails:
            jdump(run_dir / "TW370_MAIN_START.json", out)
            (run_dir / "WATCH_ABORT").write_text(now() + "\n" + "\n".join(fails) + "\n")
            if any(f.startswith("資源: ") for f in fails):
                (run_dir / "RESOURCE_STOP").write_text(now() + "\n" + "\n".join(res["failures"]) + "\n")
            log("本段の直前の照合が不成立 — 本段の forge を起動させない: " + "; ".join(fails))
            raise RuntimeError("V-tw370: 本段の直前の照合が不成立 (WATCH_ABORT)")
        if self.arm == "C1":
            try:
                self.handoff = make_handoff(run_dir, rec)
                out["handoff"] = {"dir": HANDOFF_DIR, "initial_field_sha256": self.handoff["initial_field_sha256"],
                                  "n_inputs": len(self.handoff["inputs"]), "n_history": len(self.handoff["history"])}
                log(f"C50 への受け渡しを作った ({HANDOFF_DIR}/、本段の入力 {self.handoff['initial_field_sha256'][:16]})")
            except Exception as e:  # noqa: BLE001
                out["handoff_error"] = traceback.format_exc()
                (run_dir / "HANDOFF_FAILED").write_text(now() + "\n" + repr(e) + "\n")
                log(f"C50 への受け渡しが失敗 (C1 は続ける; C50 は判定不能になる): {e!r}")
        cmd = [sys.executable, str(HERE / "m9_watch.py"), run_dir.name] + M9_ARGS
        with open(run_dir / "m9_watch.log", "a") as lf:
            self.m9 = subprocess.Popen(cmd, cwd=str(HERE), stdout=lf, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                       start_new_session=True)
        out.update(m9_cmd=cmd, m9_pid=self.m9.pid, m9_started=time.time(), m9_watch_sha256=sha256_file(HERE / "m9_watch.py"))
        jdump(run_dir / "TW370_MAIN_START.json", out)
        log(f"本段の通知を確かめ、見張りを起動 (pid {self.m9.pid}): {' '.join(cmd[1:])}")
        return rec


# --- 30 秒ごとの時刻・ディスク -------------------------------------------------------------------------------------------
class Timing(threading.Thread):
    """sp_arm.sh と同じ列の timing.tsv。run が無い間 (準備中) は溜めておき、run ができたら書く。ディスクの空き < 2 GiB で自分の forge を止める。"""

    def __init__(self, run: Path):
        super().__init__(daemon=True)
        self.rd, self.buf, self.stop_evt, self.disk_stopped = run, [], threading.Event(), False

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


def die(msg: str, code: int = 1):
    log("中止: " + msg)
    sys.exit(code)


def fail_record(run: Path, pc: dict, failures: list, t_start: float):
    """前提の不成立を run に残す (forge は回さない; 判定不能)。"""
    run.mkdir(parents=True, exist_ok=True)
    pc["failures"] = list(pc.get("failures", [])) + list(failures)
    pc["ok"] = False
    (run / "START_EPOCH").write_text(f"{t_start:.6f}\n")
    jdump(run / "TW370_PRECHECK.json", pc)
    (run / "PRECHECK_FAILED").write_text(now() + "\n" + "\n".join(pc["failures"]) + "\n")


def maybe_judge(no_judge: bool):
    """3 本とも終わっていたら (RUN_RC か PRECHECK_FAILED)、最初にここへ来た腕が 1 回だけ tw370_judge.py を回す。"""
    if no_judge:
        return
    done = all((HERE / a["run"] / "RUN_RC").is_file() or (HERE / a["run"] / "PRECHECK_FAILED").is_file() for a in ARMS.values())
    if not done:
        log("ほかの腕がまだ終わっていない — 判定は最後に終わった腕が回す (手で回すなら python3 tw370_judge.py)")
        return
    try:
        fd = os.open(str(HERE / "_band_ab/cold_pair/tw370_judge.started"), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.write(fd, f"{now()} {os.getpid()}\n".encode()); os.close(fd)
    except FileExistsError:
        log("判定はほかの腕が回している (tw370_judge.started がある)")
        return
    j = subprocess.run([sys.executable, str(HERE / "tw370_judge.py")], cwd=str(HERE), capture_output=True, text=True)
    print(j.stdout[-6000:] + j.stderr[-3000:], flush=True)


def main(argv=None) -> int:
    global LOG_TAG
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True, choices=sorted(ARMS))
    ap.add_argument("--commit", required=True, help="~/forge-g0cchain の HEAD であるべき commit (tw370_prep.sh に渡したもの)")
    ap.add_argument("--wait-c1-max", type=int, default=6 * 3600, help="C50: C1 の受け渡しを待つ上限 [s]")
    ap.add_argument("--c1-dir", default=None, help="C50 の --dry-prep だけ: 受け渡しを読む C1 の dry の dir (_tw370_dry_<名前>)")
    ap.add_argument("--no-judge", action="store_true")
    ap.add_argument("--dry-prep", default=None, metavar="_tw370_dry_<名前>",
                    help="準備と前提の検査だけを、run_08xx を使わずに別のディレクトリで行う (forge・見張り・判定は回さない)。"
                         "C1 は準備の後の nozzle.h5 で受け渡しの写しも作る (C50 の --dry-prep の --c1-dir に渡す)")
    a = ap.parse_args(argv)
    arm = ARMS[a.arm]
    LOG_TAG = f"tw370_run {a.arm}"
    if a.dry_prep is not None and not re.fullmatch(r"_tw370_dry_[A-Za-z0-9_.-]+", a.dry_prep):
        die("--dry-prep の名前は _tw370_dry_<英数字> (run_* と取り違えない)")
    if a.c1_dir is not None and (a.arm != "C50" or a.dry_prep is None or not re.fullmatch(r"_tw370_dry_[A-Za-z0-9_.-]+", a.c1_dir)):
        die("--c1-dir は C50 の --dry-prep のときだけ、_tw370_dry_<名前> を渡す")
    run = HERE / (a.dry_prep or arm["run"])
    log(f"開始 {PLAN} 腕 {a.arm} run {run}")
    # --- 環境・コード・バイナリ ---
    env = {k: v for k, v in os.environ.items() if k.startswith("FORGE_")}
    if env != ENV_WANT:
        die(f"FORGE_* の環境変数が登録 (sp_arm.sh と同じ) と違う: {env} (tw370_run.sh から起動すること)")
    prep_done = HERE / "tw370_prep.done"
    if not prep_done.is_file() or a.commit not in prep_done.read_text():
        die(f"tw370_prep.done が無いか commit {a.commit} と合わない")
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
    for n in (arm["problem"], "tw370_run.py", "tw370_judge.py", "m9_watch.py", "sp_judge.py", "cold_xcheck.py", "cold_series.py"):
        if not (HERE / n).is_file() or not (WT / CASE / n).is_file() or (HERE / n).read_bytes() != (WT / CASE / n).read_bytes():
            die(f"{n} が無いか作業ツリーの commit の版と違う")
    code = {"forge_design": str(fd), "worktree": str(WT), "worktree_head": head, "forge_root": str(RA.FORGE_TOOLS.parents[1]),
            "python": sys.executable, "python_version": sys.version.split()[0], "numpy": numpy.__version__,
            "forge_bin": str(BIN), "forge_sha256": bin_sha, "converter": str(conv), "converter_sha256": conv_sha, "env": env,
            "tools_sha256": {n: sha256_file(RA.FORGE_TOOLS / n)[:16] for n in ("run_case.sh", "interp_field.py", "check_mesh_quality.py",
                                                                              "stage_manifest.py", "restart_field.py", "rerun_conditions.py",
                                                                              "check_convergence.py") if (RA.FORGE_TOOLS / n).is_file()},
            "case_scripts_sha256": {n: sha256_file(HERE / n)[:16] for n in ("tw370_run.py", "tw370_judge.py", "m9_watch.py", "cold_series.py",
                                                                           "cold_xcheck.py", "sp_judge.py", arm["problem"]) if (HERE / n).is_file()}}
    log(f"コード {json.dumps(code, ensure_ascii=False)}")
    # --- 入力 (ドナー・壁・参照) ---
    donor = HERE / DONOR
    rs = res_steps(donor)
    pn = sorted(donor.glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))   # prepare_ns の ic_from の選び方そのまま
    if not rs or rs[-1].name != DONOR_RES or not pn or pn[-1].name != DONOR_RES:
        die(f"ドナー {DONOR} の最後の res が {rs[-1].name if rs else None} / prepare_ns の選び方で {pn[-1].name if pn else None} "
            f"({DONOR_RES} であること; prepare_ns の ic_from は最後の res を移す)")
    donor_sha = sha256_file(rs[-1])
    dm9 = jload(donor / "m9_watch.json") or {}
    if dm9.get("status") != "REACHED" or dm9.get("reach_step") != DONOR_STEP or dm9.get("reach_sha256") != donor_sha:
        die(f"ドナーが停止の出力でない (m9_watch {dm9.get('status')} {dm9.get('reach_step')}、sha の一致 {dm9.get('reach_sha256') == donor_sha})")
    q300 = jload(HERE / Q300)
    if not q300 or q300.get("donor", {}).get("sha256") != donor_sha or not q300.get("registered_band_values_match"):
        die(f"{Q300} が無いか、ドナーの sha が違うか、300 K の帯の感度が登録の値と合わない (tw370_judge.py --fix-q300 を tw370_prep.sh で)")
    prep_sha = {n: sha256_file(HERE / PREP / n) for n in PREP_SHA}
    bad = {n: s[:16] for n, s in prep_sha.items() if not s.startswith(PREP_SHA[n])}
    if bad:
        die(f"{PREP} の sha256 が登録と違う: {bad}")
    ref6 = HERE / PHYS_REF
    if not (HERE / CFG_REF / "solverConfig.yaml").is_file() or not ref6.is_dir():
        die(f"{CFG_REF}/solverConfig.yaml か {PHYS_REF} が無い")
    mref, mref_note = mesh_ref_path(arm)
    if mref is None:
        die(mref_note)
    # --- 資源 (訂正 6) ---
    res0 = resources("起動の前")
    log(f"資源 {json.dumps(res0, ensure_ascii=False)}")
    if not res0["ok"] and a.dry_prep is None:
        die("資源が足りない (forge を回さない): " + "; ".join(res0["failures"]))
    t_start = time.time()
    pc = {"plan": PLAN, "arm": a.arm, "created": now(), "run": arm["run"], "code": code, "t_start": t_start, "resources_start": res0,
          "inputs": {"donor": {"run": DONOR, "res": DONOR_RES, "sha256": donor_sha, "m9_reach_step": dm9.get("reach_step"),
                               "nozzle_h5": str((donor / "nozzle.h5").resolve()) if (donor / "nozzle.h5").exists() else None},
                     "q300": {"file": Q300, "sha256": sha256_file(HERE / Q300)}, "prep_sha256": prep_sha,
                     "euler_ref": {"run": EULER_REF, "exists": (HERE / EULER_REF).is_dir(),
                                   "note": "prepare_ns は euler_ref のパスを prepare_info.json に書くだけ (判定の δ_E・δ_target の抽出が読む)"},
                     "mesh_ref": mref_note, "phys_ref": PHYS_REF, "cfg_ref": CFG_REF},
          "failures": [], "ok": False}
    if a.arm == "C50":
        return main_c50(a, RA, arm, run, pc, t_start, mref, ref6, donor)
    # --- R0・C1: 準備 ---
    timing = Timing(run)
    if a.dry_prep is None:
        timing.start()
    log("prepare_ns を呼ぶ" + (" (ic_from なし、後で restart_field)" if arm["transfer"] == "restart_field" else f" (ic_from={DONOR})"))
    try:
        info = RA.prepare_ns(HERE / arm["problem"], run, ic_from=(donor if arm["transfer"] == "interp_field" else None),
                             delta_r_csv=str(HERE / PREP / "delta_r_initial.csv"), offset="radial", euler_ref=str(HERE / EULER_REF))
    except Exception as e:  # noqa: BLE001
        tb = traceback.format_exc()
        if isinstance(e, subprocess.CalledProcessError):
            tb += f"\n--- stdout ---\n{(e.stdout or '')[-3000:]}\n--- stderr ---\n{(e.stderr or '')[-3000:]}"
        log("prepare_ns が失敗:\n" + tb)
        if run.is_dir():
            fail_record(run, pc, [f"prepare_ns が失敗: {e!r}"], t_start)
        timing.stop_evt.set()
        if a.dry_prep is None:
            maybe_judge(a.no_judge)
        return 2
    t_prep = time.time()
    (run / "START_EPOCH").write_text(f"{t_start:.6f}\n")
    pc["t_prepared"] = t_prep
    pc["prepared_files"] = sorted(f.name for f in run.iterdir())
    log(f"prepare_ns 終了 ({t_prep - t_start:.0f} s)")
    p = RA.load_problem(HERE / arm["problem"])
    items = {}
    # --- (7) 移し ---
    if arm["transfer"] == "restart_field":
        try:
            items["7_transfer"] = transfer_restart(RA, run, rs[-1])
            items["7_transfer"]["omega_floor"] = omega_floor_count(p, run / "nozzle.h5", "would_raise")
        except Exception as e:  # noqa: BLE001
            items["7_transfer"] = {"ok": False, "failures": [f"(7) restart_field の移しの失敗: {e!r}"]}
    else:
        rec7 = {"tool": "interp_field.py (prepare_ns の ic_from)", "donor_res": str(rs[-1]), "ok": True, "failures": [],
                "correction_after_transfer": "ω の粘性底層の床 6 ν/(0.075 d²) (prepare_ns)", "nozzle_sha256_after": sha256_file(run / "nozzle.h5")}
        try:
            rec7["nearest"] = nearest_record(rs[-1], run / "nozzle.h5", float(info["scale_m"]))
            rec7["omega_floor"] = omega_floor_count(p, run / "nozzle.h5", "active")
        except Exception as e:  # noqa: BLE001
            rec7["record_error"] = repr(e)                  # 記録の失敗は検査の不成立にしない
        items["7_transfer"] = rec7
    pc["failures"] += items["7_transfer"]["failures"]
    # prepare_info.json に移し方を残す (R0 の prepare_info の ic_from は None のまま = prepare_ns は移していない)
    info["tw370_transfer"] = {k: items["7_transfer"].get(k) for k in ("tool", "verdict", "correction_after_transfer", "nozzle_sha256_after")}
    info["tw370_transfer"].update(arm=a.arm, donor=f"{DONOR}/{DONOR_RES}", donor_sha256=donor_sha)
    (run / "prepare_info.json").write_text(json.dumps(info, indent=1))
    # --- 前提の検査 (1)〜(6) ---
    msh_sha = sha256_file(run / "nozzle.msh")
    items["1_msh"] = {"sha256": msh_sha, "want": arm["msh_sha"], "grid": arm["grid"], "ok": msh_sha.startswith(arm["msh_sha"])}
    if not items["1_msh"]["ok"]:
        pc["failures"].append(f"(1) nozzle.msh の sha256 {msh_sha[:16]} が {arm['grid']} {arm['msh_sha']} と違う")
    try:
        items["2_mesh"] = {"ref": str(mref), **compare_mesh_h5(run / "nozzle.h5", mref)}
    except Exception as e:  # noqa: BLE001
        items["2_mesh"] = {"ref": str(mref), "ok": False, "failures": [f"(2) 照合の失敗: {e!r}"]}
    pc["failures"] += items["2_mesh"]["failures"]
    items["3_quality"] = quality_record(RA, run, int(p.mesh["nj"]), float(info["scale_m"]), arm["ar_max"])
    pc["failures"] += items["3_quality"]["failures"]
    items["4_sha"] = {"donor_res": {"path": str(rs[-1]), "sha256": donor_sha}, "nozzle_h5_after_transfer": sha256_file(run / "nozzle.h5"),
                      "prepare_info_ic_from": info.get("ic_from"), "main_initial_field": "本段の通知 (TW370_MAIN_START.json の initial_field_sha256)"}
    main_txt = (run / "solverConfig.yaml").read_text()
    pre_p = run / RA.PRESTAGE_CONFIG
    if not pre_p.is_file():
        items["5_config"] = {"ok": False, "failures": [f"(5) {RA.PRESTAGE_CONFIG} が無い (main_line が効いていない)"]}
    else:
        bc_p = run / "bcondConfig.yaml"
        items["5_config"] = check_configs(RA, main_txt, pre_p.read_text(), (HERE / CFG_REF / "solverConfig.yaml").read_text(), arm["cap"],
                                          explicit_record(run), bc_p.read_text() if bc_p.is_file() else None)
        nsc = info.get("ns_stage_configs") or {}
        st = items["5_config"].get("stages", {})
        rec_vs = {k: nsc.get(k) == st.get(k, {}).get("sha256") for k in ("S1_soft", "S2_mid", "main")}
        rec_vs["prestage_config"] = nsc.get("prestage_config") == st.get("prestage_config", {}).get("sha256")
        items["5_config"]["prepare_info_sha_match"] = rec_vs
        if not all(rec_vs.values()):
            items["5_config"]["failures"].append(f"(5) prepare_info の ns_stage_configs が導いた段の設定と違う: {rec_vs}")
        ml = info.get("main_line") or {}
        if not (ml.get("cap") == arm["cap"] and ml.get("itj") == 5 and ml.get("cfl") == 4.0 and ml.get("relax") == RELAX
                and info.get("implicit_relax") == RELAX and (info.get("wall_thermal") or {}).get("Tw") == 370.0):
            items["5_config"]["failures"].append(f"(5) prepare_info の main_line・壁温が登録と違う: {ml}、{info.get('wall_thermal')}")
        items["5_config"]["ok"] = not items["5_config"]["failures"]
    pc["failures"] += items["5_config"]["failures"]
    items["6_files"] = check_phys_files(run, ref6, late=False)
    items["6_files"]["vs_donor (記録だけ)"] = phys_record_vs_donor(run, donor)
    pc["failures"] += items["6_files"]["failures"]
    items["wall_record"] = {n: ((run / n).read_bytes() == (HERE / PREP / n).read_bytes()) if (run / n).is_file() and (HERE / PREP / n).is_file()
                            else "片側に無い" for n in ("wall_repr.json", "wall_physical.csv")}
    pc["items"] = items
    pc["ok"] = not pc["failures"]
    pc["t_checked"] = time.time()
    pc["dry_prep"] = a.dry_prep
    jdump(run / "TW370_PRECHECK.json", pc)
    if a.dry_prep is not None:
        if a.arm == "C1" and pc["ok"]:
            try:
                h = make_handoff(run, {}, dry=True)
                log(f"--dry-prep: C50 の確かめ用の受け渡し ({HANDOFF_DIR}/、準備の後の nozzle.h5 {h['initial_field_sha256'][:16]}) を作った")
            except Exception as e:  # noqa: BLE001
                log(f"--dry-prep: 受け渡しの写しが失敗: {e!r}")
                return 2
        log(f"--dry-prep {a.dry_prep}: 前提の検査 {'合格' if pc['ok'] else '不成立'} (forge は回さない): " + ("; ".join(pc["failures"]) or "なし"))
        return 0 if pc["ok"] else 2
    if not pc["ok"]:
        log("前提の検査が不成立 — forge を回さない (判定不能):\n  " + "\n  ".join(pc["failures"]))
        (run / "PRECHECK_FAILED").write_text(now() + "\n" + "\n".join(pc["failures"]) + "\n")
        timing.stop_evt.set()
        maybe_judge(a.no_judge)
        return 2
    log(f"前提の検査合格: msh {msh_sha[:16]}、格子のビット一致 {items['2_mesh']['n_bit_equal']} データセット、品質 "
        f"{items['3_quality']['check_mesh_quality_nozzle_h5']}、AR>5000 {items['3_quality'].get('ar_gt5000')}、移し {items['7_transfer'].get('tool')}、"
        f"(6) 未作成 {items['6_files']['pending']}")
    st5 = items["5_config"]["stages"]
    prepared = {"main_sha256": st5["main"]["sha256"], "prestage_sha256": st5["prestage_config"]["sha256"],
                "want_prev": [("S1_soft", 3000, st5["S1_soft"]["sha256"]), ("S2_mid", 3000, st5["S2_mid"]["sha256"])]}
    return run_arm(a, RA, run, items["5_config"]["stage_tags_by_sha"], prepared, ref6, "full", timing, t_start)


def run_arm(a, RA, run: Path, tags: dict, prepared: dict, ref6: Path, stages: str, timing, t_start: float) -> int:
    """段階起動 (R0・C1 は full、C50 は none) を Observer の下で回し、終わりの記録を残す。"""
    obs = Observer(RA, run, tags, prepared, ref6, a.arm)
    obs.install()
    run_end = {"t_start": t_start, "arm": a.arm, "stages": stages}
    try:
        rc = RA.run_staged_ns(run, stages=stages)
        rc_txt = str(rc)
    except Exception:  # noqa: BLE001
        run_end["exception"] = traceback.format_exc()
        log("run_staged_ns が例外:\n" + run_end["exception"])
        rc_txt = "EXC"
    finally:
        obs.uninstall()
    t_end = time.time()
    (run / "RUN_RC").write_text(rc_txt + "\n")
    run_end.update(t_end=t_end, run_rc=rc_txt, main_started=(run / RA.MAIN_STAGE_STARTED).is_file(),
                   handoff=(obs.handoff or {}).get("initial_field_sha256"))
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
    log(f"{run.name} 終了 RUN_RC={rc_txt} 見張り {m9.get('status')} {m9.get('reach_step')} (総 {t_end - t_start:.0f} s)")
    maybe_judge(a.no_judge)
    return 0


def main_c50(a, RA, arm: dict, run: Path, pc: dict, t_start: float, mref: Path, ref6: Path, donor: Path) -> int:
    """C50: C1 の受け渡しを待ち、本段の入力をバイトの写しで受け取り、C50 の設定で段なしの本段だけを回す (訂正 4)。"""
    c1 = HERE / (a.c1_dir or ARMS["C1"]["run"])
    hd = c1 / HANDOFF_DIR
    t_wait = time.time()
    log(f"C1 ({c1.name}) の受け渡し {HANDOFF_DIR}/{HANDOFF_JSON} を待つ (上限 {a.wait_c1_max} s)")
    why = None
    while not (hd / HANDOFF_JSON).is_file():
        for f, msg in (("PRECHECK_FAILED", "C1 の前提の検査が不成立"), ("HANDOFF_FAILED", "C1 の受け渡しの作成が失敗"),
                       ("WATCH_ABORT", "C1 の本段の直前の照合が不成立"), ("RUN_RC", "C1 が本段に入らずに終わった")):
            if (c1 / f).is_file():
                why = f"{msg} ({c1.name}/{f})"
        if why or time.time() - t_wait > a.wait_c1_max or (a.dry_prep is not None and not c1.is_dir()):
            break
        time.sleep(30)
    if not (hd / HANDOFF_JSON).is_file():
        why = why or (f"C1 の受け渡しが {a.wait_c1_max} s 来ない" if c1.is_dir() else f"{c1} が無い")
        log("C50 を始めない: " + why)
        pc["c50_not_started"] = {"reason": why, "c1_main_started": (c1 / RA.MAIN_STAGE_STARTED).is_file(), "waited_s": time.time() - t_wait}
        fail_record(run, pc, [f"C1 から分岐できない: {why}"], t_start)
        jdump(run / "C50_NOT_STARTED.json", pc["c50_not_started"])
        if a.dry_prep is None:
            maybe_judge(a.no_judge)
        return 2
    t_got = time.time()
    hj = json.loads((hd / HANDOFF_JSON).read_text())
    log(f"受け渡しを受け取った ({t_got - t_wait:.0f} s 待った、本段の入力 {hj['initial_field_sha256'][:16]}、dry {hj.get('dry')})")
    res_w = resources("C1 を待った後", DISK_MAIN_MIN)       # ほかの腕が既に書いた出力を二重に数えない
    t_start = time.time()
    pc.update(t_start=t_start, wait_c1_s=t_got - t_wait, handoff=hj, resources_after_wait=res_w)
    if not res_w["ok"] and a.dry_prep is None:
        fail_record(run, pc, ["資源: " + f for f in res_w["failures"]], t_start)
        (run / "RESOURCE_STOP").write_text(now() + "\n" + "\n".join(res_w["failures"]) + "\n")
        maybe_judge(a.no_judge)
        return 2
    fails = []
    if bool(hj.get("dry")) != (a.dry_prep is not None):
        fails.append(f"受け渡しの dry ({hj.get('dry')}) がこの実行と合わない")
    for n, s in list(hj["inputs"].items()) + list(hj["history"].items()):
        if not (hd / n).is_file() or sha256_file(hd / n) != s:
            fails.append(f"受け渡しの {n} が無いか sha256 が HANDOFF.json と違う")
    c1_mss = jload(c1 / RA.MAIN_STAGE_STARTED)
    if a.dry_prep is None:
        if not c1_mss or (c1_mss.get("initial_field") or {}).get("sha256") != hj["initial_field_sha256"]:
            fails.append("受け渡しの nozzle.h5 が C1 の MAIN_STAGE_STARTED の本段の初期の場と違う")
    if fails:
        fail_record(run, pc, fails, t_start)
        if a.dry_prep is None:
            maybe_judge(a.no_judge)
        return 2
    # --- run を作り、受け渡しを写す ---
    run.mkdir(parents=False, exist_ok=False)
    timing = Timing(run)
    if a.dry_prep is None:
        timing.start()
    (run / "START_EPOCH").write_text(f"{t_start:.6f}\n")
    for n in hj["inputs"]:
        shutil.copy2(hd / n, run / n)
    (run / "C1_HISTORY").mkdir()
    for n in hj["history"]:
        shutil.copy2(hd / n, run / "C1_HISTORY" / n)
    shutil.copy2(hd / HANDOFF_JSON, run / "C1_HISTORY" / HANDOFF_JSON)
    # --- C50 の設定 (問題 YAML から; C1 の本段の設定との差は上限だけ) ---
    items = {}
    p = RA.load_problem(HERE / arm["problem"])
    ml = RA._main_line_spec(p)
    n = int(p.evaluate.get("nStepOuter", 48000))
    out_int = int(p.evaluate.get("outStepInterval", max(n // 6, 1)))
    f5 = []
    if n % out_int:
        f5.append(f"(5) nStepOuter {n} が outStepInterval {out_int} で割り切れない (prepare_ns は間隔を直すので、ここでは止める)")
    # 生成器が明示すべきキーを足す版 (plan tooling-generator-defaults-and-launch-gate §5.1 #11) は、足すキーを run の bcond (C1 の準備の
    # もの) で決め、足したキーの記録 ex50 を (5) の照合に渡す (prepare_ns の RECIPE_CONTEXT.json と同じ中身)
    import inspect
    ex50 = {}
    kw50 = {"record": ex50} if "record" in inspect.signature(RA._ns_configs).parameters else {}
    with tempfile.TemporaryDirectory(prefix="tw370_c50_cfg_") as td:
        shutil.copy2(run / "bcondConfig.yaml", Path(td) / "bcondConfig.yaml")
        cfg_pre, cfg_main = RA._ns_configs(p, Path(td), n, out_int, None, None, ml, **kw50)
        gas_files = {f.name: (Path(td) / f.name).read_bytes() for f in Path(td).iterdir() if f.is_file() and f.name != "bcondConfig.yaml"}
    for gn, gb in gas_files.items():
        if not (run / gn).is_file() or (run / gn).read_bytes() != gb:
            f5.append(f"(5) _apply_gas_to_config が書く {gn} が C1 の準備のものと違う")
    c1_main = (hd / "solverConfig_C1_main.yaml").read_text() if (hd / "solverConfig_C1_main.yaml").is_file() else ""
    if cfg_pre != (run / RA.PRESTAGE_CONFIG).read_text():
        f5.append("(5) C50 の前の段の元の設定が C1 のものと文字列で違う")
    d50 = RA._cfg_diff_paths(RA._cfg_load(c1_main), RA._cfg_load(cfg_main)) if c1_main else {("C1 の本段の設定が無い",): (None, None)}
    if d50 != {("time", "deltaT", "lineDtDirectionalCap"): (1.0, 50.0)}:
        f5.append(f"(5) C1 の本段の設定との差が lineDtDirectionalCap (1.0 → 50.0) だけでない: {{{', '.join('.'.join(map(str, k)) for k in d50)}}}")
    if cfg_main != c1_main.replace("lineDtDirectionalCap: 1.0,", "lineDtDirectionalCap: 50.0,", 1):
        f5.append("(5) C1 の本段の設定の上限だけを 50.0 にした文字列と一致しない")
    (run / "solverConfig.yaml").write_text(cfg_main)
    (run / RA.PRESTAGE_CONFIG).write_text(cfg_pre)
    info = json.loads((run / "prepare_info.json").read_text())
    (run / "C1_HISTORY" / "prepare_info_C1.json").write_text(json.dumps(info, indent=1))
    info.update(main_line={**ml, "source": "evaluate.main_line", "prestage_config": RA.PRESTAGE_CONFIG},
                cfl_main=ml["cfl"], implicit_relax=ml["relax"], nStepOuter=n, ns_stage_configs=RA.ns_stage_sha256(cfg_main, cfg_pre),
                ic_from=f"{c1.name}/{HANDOFF_DIR}/nozzle.h5 (C1 の本段の入力 = S2 の後の場、登録の訂正 4)",
                tw370_branch={"from_run": c1.name, "handoff_initial_field_sha256": hj["initial_field_sha256"], "problem": arm["problem"],
                              "stages": "none (S1・S2 は C1 から引き継いだ履歴; C1_HISTORY/)", "dry": bool(hj.get("dry"))})
    (run / "prepare_info.json").write_text(json.dumps(info, indent=1))
    pc["prepared_files"] = sorted(f.name for f in run.iterdir())
    # (5) 本段の設定 (参照 run_0752、上限 50 を許す) と S1・S2 (C50 は回さないが、C1 が実際に回した設定と sha が同じこと)
    items["5_config"] = check_configs(RA, cfg_main, cfg_pre, (HERE / CFG_REF / "solverConfig.yaml").read_text(), arm["cap"],
                                      ex50, (run / "bcondConfig.yaml").read_text())
    items["5_config"]["failures"] = f5 + items["5_config"]["failures"]
    st = items["5_config"].get("stages", {})
    c1_prev = [(s.get("tag"), s.get("nsteps"), s.get("config_sha256")) for s in (hj.get("c1_main_stage_started") or {}).get("previous_stages", [])]
    want_prev = [("S1_soft", 3000, st.get("S1_soft", {}).get("sha256")), ("S2_mid", 3000, st.get("S2_mid", {}).get("sha256"))]
    items["5_config"]["c1_previous_stages"] = c1_prev
    if a.dry_prep is None and c1_prev != want_prev:
        items["5_config"]["failures"].append(f"(5) C1 が実際に回した S1・S2 が C50 の設定から導いた S1・S2 と違う: {c1_prev} / {want_prev}")
    items["5_config"]["ok"] = not items["5_config"]["failures"]
    pc["failures"] += items["5_config"]["failures"]
    # (1)・(2)・(3)・(6)・(8)
    msh_sha = sha256_file(run / "nozzle.msh")
    items["1_msh"] = {"sha256": msh_sha, "want": arm["msh_sha"], "grid": arm["grid"], "ok": msh_sha.startswith(arm["msh_sha"])}
    if not items["1_msh"]["ok"]:
        pc["failures"].append(f"(1) nozzle.msh の sha256 {msh_sha[:16]} が {arm['grid']} {arm['msh_sha']} と違う")
    try:
        items["2_mesh"] = {"ref": str(mref), **compare_mesh_h5(run / "nozzle.h5", mref)}
    except Exception as e:  # noqa: BLE001
        items["2_mesh"] = {"ref": str(mref), "ok": False, "failures": [f"(2) 照合の失敗: {e!r}"]}
    pc["failures"] += items["2_mesh"]["failures"]
    items["3_quality"] = quality_record(RA, run, int(p.mesh["nj"]), float(info["scale_m"]), arm["ar_max"])
    pc["failures"] += items["3_quality"]["failures"]
    items["4_sha"] = {"handoff_nozzle_h5": hj["initial_field_sha256"], "nozzle_h5": sha256_file(run / "nozzle.h5"),
                      "c1_main_initial_field": (c1_mss or {}).get("initial_field"), "donor_res_of_c1": {"run": DONOR, "res": DONOR_RES,
                                                                                                       "sha256": pc["inputs"]["donor"]["sha256"]}}
    items["6_files"] = check_phys_files(run, ref6, late=False)
    pc["failures"] += items["6_files"]["failures"]
    items["8_start_identity"] = start_identity(run / "nozzle.h5", hd / "nozzle.h5")
    pc["failures"] += items["8_start_identity"]["failures"]
    pc["items"] = items
    pc["ok"] = not pc["failures"]
    pc["t_checked"] = time.time()
    pc["dry_prep"] = a.dry_prep
    jdump(run / "TW370_PRECHECK.json", pc)
    if a.dry_prep is not None:
        log(f"--dry-prep {a.dry_prep}: 前提の検査 {'合格' if pc['ok'] else '不成立'} (forge は回さない): " + ("; ".join(pc["failures"]) or "なし"))
        return 0 if pc["ok"] else 2
    if not pc["ok"]:
        log("前提の検査が不成立 — forge を回さない (判定不能):\n  " + "\n  ".join(pc["failures"]))
        (run / "PRECHECK_FAILED").write_text(now() + "\n" + "\n".join(pc["failures"]) + "\n")
        timing.stop_evt.set()
        maybe_judge(a.no_judge)
        return 2
    log(f"前提の検査合格 (C50): 本段の入力 = C1 の本段の入力 ({items['8_start_identity']['sha256'][:16]}、{items['8_start_identity']['n_datasets']} "
        f"データセット)、C1 の本段の設定との差は上限だけ、格子のビット一致 {items['2_mesh']['n_bit_equal']}")
    prepared = {"main_sha256": st["main"]["sha256"], "prestage_sha256": st["prestage_config"]["sha256"], "want_prev": [],
                "initial_field_sha256": hj["initial_field_sha256"]}
    return run_arm(a, RA, run, {st["main"]["sha256"]: "main"}, prepared, ref6, "none", timing, t_start)


if __name__ == "__main__":
    sys.exit(main())
