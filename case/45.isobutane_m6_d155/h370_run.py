"""plan tooling-nozzle-isothermal-wall-chain §5.1 #36f・§6「V-h370 の事前登録」(§4.9 codex diagnose 3、
`notes/reviews/2026-10-14-h370-registration-diagnose.md`) の実行: 「H・B0 による暫定の設計」の最初の固定壁の pass。壁は更新しない。
B で h370_run.sh から起動する (h370_prep.sh <commit> を先に通す)。

  順番 (§6 V-h370): (i) Euler 参照 → 資格の判定 → 固定 → (ii) NS。
  (i)  --phase euler: run_0950_h370_euler_ref = problem_d155_euler_h_ref.yaml を runner_axismach の Euler の経路
       (prepare → run_staged: soft 3,000 → mid 3,000 → 本段 60,000 (出力 1,000・cfl_main 2.0 は YAML から)) で回す。
       forge の前の検査 (§6「投入前」): 正式の Euler の壁の保存・復元の照合 (壁ファイル → 格子の生成 → nozzle.msh とバイト一致)、
       Euler の格子の check_mesh_quality.py (--ar-max は mesh_euler の ar_max、無ければ 1000) の VERDICT、さらに NS の格子を
       停止規則を使わない準備 (prepare_ns(stop_rule=False)、_h370_prelaunch_ns) で作り、NS の check_mesh_quality.py --ar-max 5000、
       変換後の座標の delta_in_check (≤ 0.001 mm・全列の被覆; metrics.delta_in_check.check_runs)、NS が停止の出力で壁の出力 (局所 y₁⁺ の
       判定 check_wall_resolution.py に要る量) を書くか (prepare の時点の config で)、実効設定の sha。どれかが不成立なら forge を回さない。
       NS の正式の準備は Euler の参照の固定が要る (停止規則の帯・δ の記録が読む) ので、ここでは停止規則を使わない準備で格子と設定を確かめ、
       (ii) で正式の準備の格子・設定がこれとビット同一であることを確かめる。
       終了後に h370_judge.py --euler-prequalify で予備の資格を判定し (全残差の check_convergence.py --segment PASS・NaN/Inf なし・出口 M・
       流量・参照場の時系列の check_quasisteady.py の VERDICT; §5.1 #36g (1))、満たせば本段の最後の出力と sha256 を
       _band_ab/half_inlet/h370_euler_ref_pin.json に固定し、σ_ref の代わりの参照 (末尾の窓の全 20 出力と sha、§5.1 #36i) を記録する。満たさなければ NS に
       進まない。最終の資格は NS の判定の σ_ref (h370_judge.py) で決める。このスクリプトと runner_axismach.run_staged は Euler の本段の出力を
       消さない (run_staged が消すのは soft・mid の段の res_* だけで、本段の前)。
  (ii) --phase ns: 固定した {run, res, sha256} を euler_ref にして run_0951_h370_ns_b0 = problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml を
       prepare_ns(delta_r_csv=_band_ab/half_inlet/delta_r_initial_C.csv (sha 24251cb6…), offset radial) で準備し (CFD ピンは問題 YAML の
       run_0062_euler_wallfit_fit_r1_ext6k)、同じ検査 (品質・δ_in・壁の出力・設定・停止規則の設定・(i) の準備とのビット同一) の後、
       run_staged_ns(stages="full") で S1 → S2 → 本段 (停止規則は YAML のキーどおり: 上限 200,000・保持 40,000・窓 20,000・断面 40・70・94 r_t)。
       本段の見張りは run_staged_ns の中の stop_rule.watch_main_stage。本段の直前 (MAIN_STAGE_STARTED を書いた後・forge の前) に通知の中身・
       資源を照合し、合わなければ本段の forge を起動させない。終わったら h370_judge.py (§6 V-h370 の判定) を回す。
  --phase all: (i) で資格を満たしたら同じプロセスで (ii) へ。
  --phase ns-diag (§5.1 #36l・§6「V-h370-diag の事前登録」、codex diagnose `notes/reviews/2026-10-14-h370-diag-ns-diagnose.md`): **診断用**の
       run_0951_h370_ns_b0。Euler の参照は run_0950 の res_60000.h5 を**未適格 (予備の資格 FAIL) のまま**診断専用に sha で固定し、正式の記録
       (h370_euler_prequalify.json・h370_euler_ref_pin.json) とは別の _band_ab/half_inlet/h370_diag_euler_ref.json に残す (正式の記録は読むだけ;
       read_pin の資格の検査を迂回するのはこの経路の中だけ)。感度の判定の参照 (run_0952・run_0953 の res_40000.h5) の sha も同じ記録に固定する。
       準備・投入前の検査 (品質・δ_in・壁の出力・設定・停止規則・(i) の準備とのビット同一)・本段の直前の照合・NaN/Inf の見張り (停止規則)・資源は
       (ii) と同じで、足すのは保存の容量の照合 (上限 200,000 / 出力 2,500 = 80 出力 × 1 出力の大きさ (投入前の準備の nozzle.h5) + 壁の出力 + 余裕
       をディスクの空きと比べ、足りなければ forge を回さない) と、抽出器・見張りの版 (design の該当ファイルの sha) の記録。受け渡しの関門
       (deltastar_loop.handoff_gate) は**呼ばない**、壁の更新・次の pass の自動実行もしない、正式の判定 (h370_judge.py の既定) も呼ばない
       (update_conditions_met は常に偽; run に印 H370_DIAG_ONLY.json)。自動の延長はしない。このスクリプト・run_staged_ns・stop_rule は本段の
       出力を消さない (run_staged_ns が消すのは S1・S2 の出力だけで本段の前、stop_rule.delete_after が消すのは停止の出力 s より後の出力だけ)。
       終わったら終了の状態 (REACHED・CENSORED・実行の失敗) と停止の時刻 s・判定の窓 (CENSORED は「打切りの時点の診断の窓」= 末尾の 17 出力)・
       本段の出力の完全性・run の全ファイルの sha の一覧 (H370_DIAG_SHA256.json)・check_convergence.py --segment --drop 3 の VERDICT (記録)・正式の記録が
       書き換わっていないことを H370_DIAG_END.json に残す。感度の判定は別に `python3 h370_judge.py --diag-sens` (CFD 0 step)。

  コード : ~/forge-g0cchain (h370_prep.sh が commit に合わせる) の design/ を sys.path の先頭にして forge_design を読む。
  forge  : ~/forge-fgeom7-fp64 (sha256 7e95e5bbd8f50c9e…)、変換器は同じビルドの convertGmshToForge (711f4438…)。環境変数は V-tw370 と同じ
           (FORGE_CUDA_BLOCKSIZE=128・FORGE_OMEGA_BUDGET=1・FORGE_DIAG_COMMIT_LOSS=500・FORGE_BIN、ほかの FORGE_* は無し; h370_run.sh が整える)。
  資源   : 起動の前と本段の forge の前に、ディスクの空き ≥ 20 GiB (§5.1 #36f の指示)・GPU のメモリの空き ≥ 2 GiB・ホストの MemAvailable ≥ 4.5 GiB。
           走行中は 30 秒ごとに timing.tsv を書き、ディスクの空きが 2 GiB を切ったら cwd がこの run の forge だけを止める (DISK_STOP)、
           forge_run.log の `Write failed` を見つけたら WRITE_FAILED に残す。forge のプロセスは `pgrep -x forge` + /proc/PID/cwd で絞る。
  記録   : h370_run.log (h370_run.sh)、各 run の H370_PRECHECK.json (・PRECHECK_FAILED)・h370_stages.jsonl・timing.tsv・RUN_RC・RUN_END.json、
           NS は H370_MAIN_START.json (・WATCH_ABORT・RESOURCE_STOP)。

usage (B、h370_run.sh 経由):
  bash h370_run.sh --phase euler --commit <c> --dry-run r1   (回す前の確かめ: _h370_dry_r1_euler・_h370_dry_r1_ns に準備と検査だけ、forge は回さない)
  bash h370_run.sh --phase all --commit <c>                  ((i) → 資格 → (ii) → 判定)
  bash h370_run.sh --phase ns --commit <c>                   ((i) が済んでいて、(ii) だけを回すとき)
  bash h370_run.sh --phase ns-diag --commit <c> --dry-run d1 (診断の経路の回す前の確かめ: _h370_dry_d1_nsdiag に準備と検査だけ)
  bash h370_run.sh --phase ns-diag --commit <c>              (診断用の run_0951 (§6 V-h370-diag)。終わったら h370_judge.py --diag-sens)
  python3 h370_run.py --phase check-problems --commit <c>    (h370_prep.sh が呼ぶ: h370_problems.py check を作業ツリーの版で)
終了コード: 0 正常 (dry-run の検査の合格を含む) / 1 中止 (環境・入力) / 2 前提の検査の不成立・資格の未達・保存の容量の不足 (forge を回さない) /
           4 NS が停止の出力で壁の出力を書かない (書く設定を足すのは方針の変更なので親に返す)。
"""
import argparse
import hashlib
import importlib.util
import json
import math
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

import numpy as np

HERE = Path(__file__).resolve().parent            # B: ~/forge-wallfit/case/45.isobutane_m6_d155 (h370_prep.sh が作業ツリーから写す)
WT = Path.home() / "forge-g0cchain"
CASE = "case/45.isobutane_m6_d155"
PLAN = "plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36f・§6 V-h370 (§4.9 codex diagnose 3)"
EU_PROBLEM = "problem_d155_euler_h_ref.yaml"
NS_PROBLEM = "problem_d155_ns_prod_coldmesh_tw370_cgref_h_b0.yaml"
EU_RUN = "run_0950_h370_euler_ref"
NS_RUN = "run_0951_h370_ns_b0"
NS_PRELAUNCH = "_h370_prelaunch_ns"               # (i) の前に NS の格子・設定を確かめる準備 (停止規則なし、forge は回さない)
CSV = "_band_ab/half_inlet/delta_r_initial_C.csv"
CSV_SHA = "24251cb647226d6523109655c08c05210fdb8ba4a669aeaecb4407202945b918"
PIN = "_band_ab/half_inlet/h370_euler_ref_pin.json"
QUAL = "_band_ab/half_inlet/h370_euler_prequalify.json"   # h370_judge.QUAL と同じ (§5.1 #36g で h370_euler_qualify.json から改名)
PREP_CHECK_REC = "_band_ab/half_inlet/h370_prep_check_tr5.json"      # 手元の prep-check (生成の座標) の記録 (照合は記録だけ)
CFD_PIN_RUN = "run_0062_euler_wallfit_fit_r1_ext6k"
CFD_PIN_FILES = ("prepare_info.json", "nozzle.h5", "res_6000.h5", "solverConfig.yaml", "bcondConfig.yaml")
OLD_EULER = "run_0174_euler_v5d_M_r1"             # 抽出の参照の A/B の腕 A (判定で使う; ここでは有無と sha を記録)
BIN = Path.home() / "forge-fgeom7-fp64/solver_density_cuda/build/forge"
BIN_SHA, CONV_SHA = "7e95e5bbd8f50c9e", "711f4438"
ENV_WANT = {"FORGE_CUDA_BLOCKSIZE": "128", "FORGE_OMEGA_BUDGET": "1", "FORGE_DIAG_COMMIT_LOSS": "500", "FORGE_BIN": str(BIN)}
NUMPY_WANT = "2.5.2"
SCRIPTS = ("h370_run.py", "h370_judge.py", "h370_problems.py", EU_PROBLEM, NS_PROBLEM, CSV)
# 登録の値 (§6 V-h370; h370_problems.py の STOP_RULE_* と同じ値。試験 test_h370_launch.py が照合する)
EU_STEPS, EU_OUT, EU_CFL, EU_SOFT = 60000, 1000, 2.0, 3000
EU_TAIL = 20000                                    # σ_ref の Euler の末尾の窓 (h370_judge.EU_TAIL と同じ)
SIGMA_REF_STEPS = list(range(EU_STEPS - EU_TAIL, EU_STEPS, EU_OUT))   # σ_ref の代わりの参照 = 末尾の窓の全出力 (固定の 60,000 を除く 20 出力;
                                                   # h370_judge.sigma_ref_required() と同じ、§5.1 #36i M1)
EU_MESH = {"ni": 2000, "nj": 97, "throat_refine": 5.0}
NS_CAP, NS_OUT, NS_HOLD, NS_PRE = 200000, 2500, 40000, 3000
NS_MESH = {"ni": 4719, "nj": 160}
X_FRACS = [0.419738616669543, 0.734542579171700, 0.986385749173426]
THRESHOLDS = {"window": 20000, "min_points": 9, "r_mean_max": 5.0e-4, "r_amp_max": 5.0e-3, "theta_drift_max": 5.0e-4,
              "theta_amp_max": 5.0e-3, "x_fracs": X_FRACS}
PHYS_UP = {"r_inlet": 3.2425, "L_U": 9.0, "L_pipe": 3.5}
TW = 370.0
MAIN_LINE = {"line": "dir", "cap": 50.0, "itj": 5, "cfl": 4.0, "relax": 0.7, "prestage_relax": 0.7,
             "limiter_ref": {"limiterRoRef": 0.8739869154, "limiterPRef": 387491.0573, "limiterARef": 359.7712185}}
# NS の段の実効設定の要求 (h370_problems.py prep-check の want と同じ)
NS_STAGE_WANT = {"S1_soft": {"cfl": 0.5, "convMethod": 0, "nStepOuter": NS_PRE}, "S2_mid": {"cfl": 1.0, "convMethod": 0, "nStepOuter": NS_PRE},
                 "main": {"cfl": 4.0, "nStepOuter": NS_CAP, "outStepInterval": NS_OUT, "lineImplicit": 1, "lineDtDirectional": 1,
                          "lineDtDirectionalCap": 50.0, "implicitThermalJacobian": 5, "limiterRoRef": 0.8739869154,
                          "limiterPRef": 387491.0573, "limiterARef": 359.7712185}}
CFG_KEYS = {"cfl": ("time", "deltaT", "cfl"), "cfl_pseudo": ("time", "deltaT", "cfl_pseudo"), "convMethod": ("space", "convMethod"),
            "nStepOuter": ("time", "last", "nStepOuter"), "outStepInterval": ("time", "outStepInterval"),
            "implicitRelax": ("time", "deltaT", "implicitRelax"), "lineImplicit": ("time", "deltaT", "lineImplicit"),
            "lineDtDirectional": ("time", "deltaT", "lineDtDirectional"), "lineDtDirectionalCap": ("time", "deltaT", "lineDtDirectionalCap"),
            "implicitThermalJacobian": ("time", "deltaT", "implicitThermalJacobian"), "limiterRoRef": ("space", "limiterRoRef"),
            "limiterPRef": ("space", "limiterPRef"), "limiterARef": ("space", "limiterARef"), "output": ("output",),
            "blockDPLUR": ("time", "deltaT", "blockDPLUR"), "timeIntegration": ("time", "timeIntegration"),
            "nStepInner": ("time", "nStepInner"), "limiter": ("space", "limiter"), "wallTreatmentSST": ("turbulence", "wallTreatmentSST")}
B0_ONLY = ("lineImplicit", "lineDtDirectional", "lineDtDirectionalCap", "implicitThermalJacobian", "limiterRoRef", "limiterPRef", "limiterARef")
# S1・S2 に B0 を入れない検査 (§6 V-h370) はキーの有無でなく実効値で見る (plan tooling-generator-defaults-and-launch-gate §5.1 #26・
# §6.1 result (段取り 1) M1)。B0 だけのキーを書いていないときの実効値 = ソルバの既定 (solverConfig.cpp:372・389・390・453・713〜715;
# どれも B0 を切った値) を明示したものは書いていないのと同じ (#11 以降の生成器は前の段に lineImplicit・lineDtDirectional・
# implicitThermalJacobian 0 を明示すべきキーとして足す)。それ以外の値 (型の違う値・真偽値を含む) は B0 が残っているとみなす。
# h370_problems.B0_OFF・tw370_run.B0_OFF・vgc_run.B0_OFF と同じ (test_h370_launch.py A で照合)
B0_OFF = {"lineImplicit": 0, "lineDtDirectional": 0, "lineDtDirectionalCap": 0.0, "implicitThermalJacobian": 0,
          "limiterRoRef": 0.0, "limiterPRef": 0.0, "limiterARef": 0.0}
DISK_START_MIN = 20 * 1024 ** 3                   # 起動の前・本段の forge の前 (NS の本段は出力 1 枚 ~0.2 GB × 最大 80 枚 + 壁の出力)
DISK_MIN = 2 * 1024 ** 3                          # 走行中: これを切ったら自分の forge を止める
GPU_FREE_MIN_MIB = 2048
HOST_MIN_MIB = 4608
EXIT_WALL_OUTPUT = 4
# --- 診断の経路 (--phase ns-diag; §5.1 #36l・§6 V-h370-diag) -------------------------------------------------------------------
PLAN_DIAG = ("plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36l・§6 V-h370-diag (codex diagnose "
             "notes/reviews/2026-10-14-h370-diag-ns-diagnose.md)")
DIAG_REF = "_band_ab/half_inlet/h370_diag_euler_ref.json"   # 診断専用の Euler の参照の固定 (正式の PIN・QUAL とは別; h370_judge.DIAG_REF と同じ)
DIAG_SENS_REFS = {"A": ("run_0952_h370_euler_cfl2", "res_40000.h5"), "B": ("run_0953_h370_euler_cfl05", "res_40000.h5")}   # 感度の判定の腕
DIAG_MARK = "H370_DIAG_ONLY.json"                 # 診断専用の run の印 (正式の判定・関門に回さない)
DIAG_END = "H370_DIAG_END.json"                   # 終了の状態・判定の窓・保存の完全性の記録
DIAG_SHA = "H370_DIAG_SHA256.json"                # 終了の時点の run の全ファイルの sha256 の一覧
DIAG_CONV = "H370_DIAG_CONVERGENCE.txt"           # check_convergence.py --segment --drop 3 の出力 (記録)
DIAG_STATUS = "未適格 (予備の資格 FAIL) — 診断専用。正式の資格・固定ではなく、受け渡し・壁の更新に使わない"
DIAG_N_OUT = NS_CAP // NS_OUT                     # 本段の上限までの出力の数 (80)
DIAG_EXTRA_OUT = 2                                # S1・S2 の出力 (本段の前に消す) の分
DIAG_WALL_FRAC = 0.01                             # 1 出力の境界の出力 (res_<名前>_<id>_<step>.h5) の見積もり = nozzle.h5 × これ (case/45 run_0117
                                                  # の実測: res_wall_3 232 KB / nozzle.h5 50 MB = 0.46 %、outlet・inlet を足しても 1 % 未満)
DIAG_MARGIN = 1.10                                # 余裕 (×)、さらに走行中の停止の閾値 DISK_MIN を残す
VERSION_FILES = ("evaluate/stop_rule.py", "metrics/deltastar.py", "feedback/deltastar_loop.py", "evaluate/runner_axismach.py",
                 "evaluate/runner.py", "metrics/delta_in_check.py")   # 抽出器・見張り・準備の版 (forge_design の中の相対パス; h370_judge.VERSION_FILES と同じ)


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


LOG_TAG = "h370_run"


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
    return sorted((f for f in Path(d).glob("res_[0-9]*.h5") if re.fullmatch(r"res_\d+\.h5", f.name)), key=lambda f: int(f.stem.split("_")[1]))


def die(msg: str, code: int = 1):
    log("中止: " + msg)
    sys.exit(code)


# --- forge のプロセス (pgrep -x forge + /proc/PID/cwd; pgrep -f・pkill -f は使わない) --------------------------------------
def forge_pids():
    """(pid, cwd の realpath) の一覧。"""
    out = subprocess.run(["pgrep", "-x", "forge"], capture_output=True, text=True).stdout.split()
    r = []
    for p in out:
        try:
            r.append((int(p), os.path.realpath(f"/proc/{p}/cwd")))
        except OSError:
            r.append((int(p), "?"))
    return r


def stop_own_forge(run: Path) -> bool:
    """cwd がこの run の forge だけを止める (SIGTERM → 180 s → SIGKILL)。"""
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


# --- 資源 -----------------------------------------------------------------------------------------------------------------------
def resources(where: str, disk_min: int = DISK_START_MIN) -> dict:
    """GPU のメモリ (nvidia-smi、既存のジョブを含む)・ホストの MemAvailable・ディスクの空きを確かめる。読めないのも不成立。"""
    fails, rec = [], {"where": where, "t": time.time(), "disk_min_gib": disk_min / 2 ** 30}
    try:
        q = subprocess.run(["nvidia-smi", "--query-gpu=index,name,memory.used,memory.total", "--format=csv,noheader,nounits"],
                           capture_output=True, text=True, timeout=60)
        rows = [[c.strip() for c in ln.split(",")] for ln in q.stdout.strip().splitlines() if ln.strip()]
        if q.returncode != 0 or not rows:
            raise RuntimeError(f"rc {q.returncode}: {q.stderr.strip()[:200]}")
        used, total = float(rows[0][2]), float(rows[0][3])
        rec["gpu"] = {"index": rows[0][0], "name": rows[0][1], "used_mib": used, "total_mib": total, "free_mib": total - used}
        if total - used < GPU_FREE_MIN_MIB:
            fails.append(f"GPU のメモリの空き {total - used:.0f} MiB < {GPU_FREE_MIN_MIB} MiB")
    except Exception as e:  # noqa: BLE001
        fails.append(f"nvidia-smi で GPU のメモリを読めない ({e!r})")
    try:
        avail = int(next(ln.split()[1] for ln in Path("/proc/meminfo").read_text().splitlines() if ln.startswith("MemAvailable"))) // 1024
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


# --- 30 秒ごとの時刻・ディスク・Write failed ---------------------------------------------------------------------------------
class Timing(threading.Thread):
    """timing.tsv (時刻・最後の step・forge の本数・ほかの run の forge の cwd)。ディスクの空き < 2 GiB で自分の forge を止める (DISK_STOP)。
    forge_run.log に `Write failed` (HighFive の書き込みの失敗; ディスク満杯で forge が落ちる) が出たら WRITE_FAILED に残す。"""

    def __init__(self, run: Path):
        super().__init__(daemon=True)
        self.rd, self.buf, self.stop_evt, self.disk_stopped, self.write_failed = run, [], threading.Event(), False, False

    def last_step(self):
        try:
            with open(self.rd / "residual_history.csv", "rb") as f:
                f.seek(0, 2); n = f.tell(); f.seek(max(0, n - 8192))
                lines = [ln for ln in f.read().decode(errors="replace").splitlines() if ln.strip()]
            s = lines[-1].split(",")[0].strip() if lines else ""
            return s if s.isdigit() else ""
        except OSError:
            return ""

    def check_write_failed(self):
        if self.write_failed:
            return
        for name in ("forge_run.log", "run_case_stdout.log"):
            p = self.rd / name
            try:
                with open(p, "rb") as f:
                    f.seek(0, 2); n = f.tell(); f.seek(max(0, n - 65536))
                    tail = f.read().decode(errors="replace")
            except OSError:
                continue
            if "Write failed" in tail:
                self.write_failed = True
                line = next(ln for ln in tail.splitlines() if "Write failed" in ln)
                msg = f"{now()} {name} に Write failed (ディスクの空き {shutil.disk_usage(HERE).free / 2**30:.2f} GiB): {line.strip()[:300]}"
                log(msg)
                (self.rd / "WRITE_FAILED").write_text(msg + "\n")
                return

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
            self.check_write_failed()
        else:
            self.buf.append(row)
        free = shutil.disk_usage(HERE).free
        if free < DISK_MIN and not self.disk_stopped:
            self.disk_stopped = True
            msg = f"{now()} ディスクの空き {free / 2**30:.2f} GiB < 2 GiB — cwd がこの run の forge を止める"
            log(msg)
            if self.rd.is_dir():
                (self.rd / "DISK_STOP").write_text(msg + "\n")
            log(f"自分の forge を止めた: {stop_own_forge(self.rd)}")

    def run(self):  # noqa: D401
        while not self.stop_evt.is_set() and not (self.rd / "RUN_RC").exists():
            try:
                self.run_once()
            except Exception as e:  # noqa: BLE001
                log(f"timing の記録の失敗: {e!r}")
            self.stop_evt.wait(30)


# --- コード・環境 ---------------------------------------------------------------------------------------------------------------
def load_forge_design(commit: str) -> tuple:
    """作業ツリー WT の design/ から forge_design を読み、HEAD が commit であること・commit されていない変更が無いことを確かめる。
    戻り (runner_axismach, 記録)。"""
    sys.path.insert(0, str(WT / "design"))
    import forge_design
    from forge_design.evaluate import runner_axismach as RA
    fd = Path(forge_design.__file__).resolve()
    if not str(fd).startswith(str((WT / "design").resolve()) + os.sep):
        die(f"forge_design が作業ツリーの外から読まれた: {fd}")
    head = subprocess.run(["git", "-C", str(WT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(WT), "status", "--porcelain", "--untracked-files=no"], capture_output=True, text=True).stdout.strip()
    if len(commit) < 7 or not re.fullmatch(r"[0-9a-f]+", commit) or not head.startswith(commit):
        die(f"作業ツリーの HEAD {head} が --commit {commit} と違う")
    if dirty:
        die(f"作業ツリーに commit されていない変更がある: {dirty[:300]}")
    return RA, {"forge_design": str(fd), "worktree": str(WT), "worktree_head": head}


def code_checks(a) -> tuple:
    """環境変数・h370_prep.done・作業ツリー・numpy・バイナリ・台本と入力の版 (case dir と作業ツリーの commit が同じ) を確かめる。
    戻り (runner_axismach, 記録)。不成立は die (forge を回さない)。"""
    env = {k: v for k, v in os.environ.items() if k.startswith("FORGE_")}
    if env != ENV_WANT:
        die(f"FORGE_* の環境変数が登録 (V-tw370 と同じ) と違う: {env} (h370_run.sh から起動すること)")
    prep_done = HERE / "h370_prep.done"
    if not prep_done.is_file() or a.commit not in prep_done.read_text():
        die(f"h370_prep.done が無いか commit {a.commit} と合わない (先に h370_prep.sh <commit>)")
    RA, rec = load_forge_design(a.commit)
    import numpy
    import scipy
    if numpy.__version__ != NUMPY_WANT:
        die(f"numpy {numpy.__version__} が登録の {NUMPY_WANT} と違う")
    bin_sha = sha256_file(BIN)
    conv = RA.converter_path()
    conv_sha = sha256_file(conv) if conv.is_file() else None
    if not bin_sha.startswith(BIN_SHA) or not (conv_sha or "").startswith(CONV_SHA):
        die(f"バイナリが登録と違う (forge {bin_sha[:16]}、変換器 {conv} {str(conv_sha)[:16]})")
    for n in SCRIPTS:
        if not (HERE / n).is_file() or not (WT / CASE / n).is_file() or (HERE / n).read_bytes() != (WT / CASE / n).read_bytes():
            die(f"{n} が無いか作業ツリーの commit の版と違う")
    rec.update(python=sys.executable, python_version=sys.version.split()[0], numpy=numpy.__version__, scipy=scipy.__version__,
               forge_bin=str(BIN), forge_sha256=bin_sha, converter=str(conv), converter_sha256=conv_sha, env=env,
               tools_sha256={n: sha256_file(RA.FORGE_TOOLS / n)[:16] for n in (
                   "run_case.sh", "check_convergence.py", "check_quasisteady.py", "check_wall_resolution.py", "check_mesh_quality.py",
                   "stage_manifest.py", "restart_field.py", "rerun_conditions.py") if (RA.FORGE_TOOLS / n).is_file()},
               case_sha256={n: sha256_file(HERE / n)[:16] for n in SCRIPTS})
    return RA, rec


def input_checks() -> dict:
    """B にある入力: δ_r の CSV (sha)・CFD ピンの凍結元 run_0062 (問題 YAML の initial_line_run)・腕 A の旧い Euler run_0174 (記録)。"""
    fails, rec = [], {}
    s = sha256_file(HERE / CSV) if (HERE / CSV).is_file() else None
    rec["delta_r_csv"] = {"path": CSV, "sha256": s}
    if s != CSV_SHA:
        fails.append(f"{CSV} の sha256 {str(s)[:16]} が登録 {CSV_SHA[:16]}… と違う")
    pin = HERE / CFD_PIN_RUN
    rec["cfd_pin"] = {"run": CFD_PIN_RUN, "files": {f: (sha256_file(pin / f)[:16] if (pin / f).is_file() else None) for f in CFD_PIN_FILES}}
    miss = [f for f, v in rec["cfd_pin"]["files"].items() if v is None]
    if miss:
        fails.append(f"CFD ピンの凍結元 {CFD_PIN_RUN} に {miss} が無い")
    old = HERE / OLD_EULER
    rs = res_steps(old) if old.is_dir() else []
    rec["old_euler (腕 A、判定で使う)"] = {"run": OLD_EULER, "prepare_info": (old / "prepare_info.json").is_file(), "nozzle_h5": (old / "nozzle.h5").is_file(),
                                         "last_res": rs[-1].name if rs else None}
    rec["ok"], rec["failures"] = not fails, fails
    return rec


# --- 前提の検査 ---------------------------------------------------------------------------------------------------------------
def quality_check(RA, run: Path, ar_max: int) -> dict:
    """prepare の MESH_QUALITY.txt (cell 変換の写し) と、変換後の nozzle.h5 (node、VIZMESH の primal) への check_mesh_quality.py --ar-max の
    両方が VERDICT: PASS であること (SOFT-PASS は通さない)。"""
    fails = []
    mq = (run / "MESH_QUALITY.txt").read_text() if (run / "MESH_QUALITY.txt").is_file() else ""
    mq_line = next((ln.strip() for ln in mq.splitlines()[::-1] if "VERDICT" in ln), "")
    if not (mq_line.startswith("VERDICT: PASS") and f"AR<= {ar_max}" in mq_line):
        fails.append(f"prepare の MESH_QUALITY.txt が --ar-max {ar_max} の PASS でない ({mq_line!r})")
    q = subprocess.run([sys.executable, str(RA.FORGE_TOOLS / "check_mesh_quality.py"), str(run / "nozzle.h5"), "--ar-max", str(ar_max)],
                       capture_output=True, text=True)
    (run / "H370_QUALITY.txt").write_text(q.stdout + q.stderr)
    q_line = next((ln.strip() for ln in q.stdout.splitlines()[::-1] if "VERDICT" in ln), "")
    if q.returncode != 0 or not q_line.startswith("VERDICT: PASS"):
        fails.append(f"check_mesh_quality.py nozzle.h5 --ar-max {ar_max} が PASS でない (rc {q.returncode}, {q_line!r})")
    return {"ar_max": ar_max, "mesh_quality_txt": mq_line, "check_mesh_quality_nozzle_h5": q_line, "rc": q.returncode,
            "ok": not fails, "failures": fails}


def generator_sha256(coords, quads, bedges) -> str:
    """runner_axismach._mesh_hashes の generator_sha256 と同じ式 (h370_problems.generator_sha256 と同じ)。"""
    h = hashlib.sha256()

    def _upd(name, a):
        a = np.ascontiguousarray(a)
        h.update(name.encode() + str(a.dtype).encode() + str(a.shape).encode() + a.tobytes())
    _upd("coords", np.asarray(coords, dtype=np.float64))
    _upd("quads", np.asarray(quads, dtype=np.int64))
    for g in ("inlet", "outlet", "wall", "axis"):
        _upd(g, np.asarray(bedges[g], dtype=np.int64))
    return h.hexdigest()


def euler_wall_restore_check(RA, run: Path) -> dict:
    """正式の Euler の壁の保存・復元の照合 (§6「投入前」): run の壁ファイル (wall_repr.json、種類 design_upstream) の sha が prepare_info の
    euler_wall.wall_file と同じで、読み直した壁 (load_wall_file) から mesh_euler で生成した格子を msh に書くと run の nozzle.msh とバイト一致し、
    生成器の sha が prepare_info の mesh.hashes.generator_sha256 と同じ。手元の prep-check の記録 (生成の座標) との照合は記録だけ。"""
    from forge_design.geometry.wall_axismach import DESIGN_UPSTREAM_KIND, WALL_FILE, load_wall_file
    from forge_design.meshing.mesh2d import generate_axisym_mesh, write_msh41_2d
    fails, rec = [], {}
    info = jload(run / "prepare_info.json") or {}
    ew = info.get("euler_wall") or {}
    wf = ew.get("wall_file") or {}
    wsha = sha256_file(run / WALL_FILE) if (run / WALL_FILE).is_file() else None
    rec["wall_file"] = {"file": WALL_FILE, "sha256": wsha, "prepare_info_sha256": wf.get("sha256"), "repr": wf.get("repr")}
    if wf.get("repr") != DESIGN_UPSTREAM_KIND or wsha is None or wsha != wf.get("sha256"):
        fails.append(f"壁ファイルが prepare_info の記録と違う (種類 {wf.get('repr')!r}、sha {str(wsha)[:16]} / {str(wf.get('sha256'))[:16]})")
        return {"ok": False, "failures": fails, **rec}
    back = load_wall_file(run / WALL_FILE)
    rec["kind"] = back["kind"]
    if back["kind"] != DESIGN_UPSTREAM_KIND:
        fails.append(f"復元した壁の種類 {back['kind']!r} (design_upstream であること)")
    prob = RA.load_problem(HERE / EU_PROBLEM)
    mp = RA.mesh_params_euler(prob, float(info["scale_m"]))
    c, q, b = generate_axisym_mesh(back["physical"], mp)
    with tempfile.TemporaryDirectory(prefix="h370_wall_") as td:
        write_msh41_2d(Path(td) / "nozzle.msh", c, q, b)
        same_msh = (Path(td) / "nozzle.msh").read_bytes() == (run / "nozzle.msh").read_bytes()
    gsha = generator_sha256(c, q, b)
    want_g = ((info.get("mesh") or {}).get("hashes") or {}).get("generator_sha256")
    rec.update(msh_bytes_identical=same_msh, generator_sha256=gsha, prepare_info_generator_sha256=want_g)
    if not same_msh:
        fails.append("壁ファイルから復元した壁で作った格子の msh が run の nozzle.msh とバイト一致しない")
    if gsha != want_g:
        fails.append(f"復元した壁の格子の生成器の sha {gsha[:16]} が prepare_info の {str(want_g)[:16]} と違う")
    loc = jload(HERE / PREP_CHECK_REC) or {}
    rec["local_prep_check (記録だけ)"] = {"wall_file_sha256_same": (loc.get("euler_wall_file") or {}).get("sha256") == wsha,
                                         "generator_sha256_same": (loc.get("euler_mesh") or {}).get("generator_sha256") == gsha,
                                         "note": "手元 (numpy 2.5.2・scipy 1.18.0) の生成の座標の記録。違えば環境の差として記録する (合否にしない)"}
    rec["ok"], rec["failures"] = not fails, fails
    return rec


def wall_output_check(run: Path) -> dict:
    """NS が停止の出力で壁の出力 (局所 y₁⁺ の判定 check_wall_resolution.py --step s に要る量) を書くかを prepare の時点の設定で確かめる:
    - bcondConfig.yaml の壁 (kind wall_isothermal・Ts 370.0) が outputHDFflg: 1 → forge は出力の間隔ごとと最終 step に
      res_<名前>_<physID>_<step>.h5 (twall_x・twall_y・twall_z・ro・Ts) を書く (output.cpp outputBconds_H5_XDMF、boundaryCond.hpp)。
    - 本段の設定の output.level ≥ 1 (無ければ既定 1) → res_<step>.h5 に vis_lam (check_wall_resolution.py の μ_w。viscMethod 2 の多成分では
      vis_lam が無いと判定不能; output.cpp の level 1 の量)。
    - wallTreatmentSST 0 (低 Re、判定は局所 y₁⁺ ≤ 1)、mesh.meshFileName が run にある (check_wall_resolution.py が config から引く)。
    - 本段の outStepInterval が上限を割り切る (停止の出力 s は出力の間隔の倍数で、その step に壁の出力も書かれる)。
    書かない設定なら止めて親に返す (書く設定を足すのは方針の変更; 終了コード 4)。"""
    import yaml
    tools = WT / "solver_density_cuda" / "tools"
    sys.path.insert(0, str(tools if tools.is_dir() else HERE.parents[1] / "solver_density_cuda" / "tools"))
    fails, rec = [], {}
    try:
        bc = yaml.safe_load((run / "bcondConfig.yaml").read_text()) or {}
    except (OSError, yaml.YAMLError) as e:
        return {"ok": False, "failures": [f"bcondConfig.yaml を読めない ({e!r})"]}
    walls = {k: v for k, v in bc.items() if isinstance(v, dict) and str(v.get("kind", "")).startswith("wall")}
    rec["walls"] = {k: {"kind": v.get("kind"), "physID": v.get("physID"), "outputHDFflg": v.get("outputHDFflg"),
                        "Ts": (v.get("floats") or {}).get("Ts")} for k, v in walls.items()}
    iso = [k for k, v in walls.items() if v.get("kind") == "wall_isothermal"]
    if len(walls) != 1 or len(iso) != 1:
        fails.append(f"壁の境界が wall_isothermal 1 つでない ({rec['walls']})")
    else:
        w = walls[iso[0]]
        if int(w.get("outputHDFflg") or 0) != 1:
            fails.append(f"壁 {iso[0]} の outputHDFflg が 1 でない ({w.get('outputHDFflg')!r}) — 壁の出力 res_{iso[0]}_{w.get('physID')}_<step>.h5 が書かれない")
        if float((w.get("floats") or {}).get("Ts", float("nan"))) != TW:
            fails.append(f"壁の Ts が {TW} でない ({(w.get('floats') or {}).get('Ts')!r})")
        rec["wall_dump_name"] = f"res_{iso[0]}_{w.get('physID')}_<step>.h5"
    try:
        import yaml_strict
        cfg = yaml_strict.load((run / "solverConfig.yaml").read_text())
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "failures": fails + [f"solverConfig.yaml を読めない ({e!r})"], **rec}
    out = cfg.get("output") if isinstance(cfg.get("output"), dict) else {}
    level = int(out.get("level", 1)) if out is not None else 1
    tm = cfg.get("time") or {}
    out_int = int(tm.get("outStepInterval", 0) or 0)
    wt = int(((cfg.get("turbulence") or {}).get("wallTreatmentSST", 0)) or 0)
    mf = (cfg.get("mesh") or {}).get("meshFileName")
    nso = int(((tm.get("last") or {}).get("nStepOuter", 0)) or 0)
    rec.update(output_level=level, extra_fields=out.get("extraFields") if out else None, outStepInterval=out_int, wallTreatmentSST=wt,
               meshFileName=mf, nStepOuter=nso)
    if level < 1:
        fails.append(f"本段の output.level {level} < 1 — res_<step>.h5 に vis_lam が出ない (μ_w を引けない)")
    if wt != 0:
        fails.append(f"wallTreatmentSST {wt} (低 Re の 0 であること)")
    if not mf or not (run / str(mf)).is_file():
        fails.append(f"mesh.meshFileName {mf!r} が run に無い")
    if out_int <= 0 or NS_CAP % out_int or nso % out_int:
        fails.append(f"outStepInterval {out_int} が上限 {NS_CAP}・nStepOuter {nso} を割り切らない (停止の出力に壁の出力が揃わない)")
    rec["ok"], rec["failures"] = not fails, fails
    return rec


def eff_settings(RA, cfg: str) -> dict:
    doc = RA._cfg_load(cfg)
    return {k: (lambda v: None if v is RA._CFG_MISSING else v)(RA._cfg_path_get(doc, pth)) for k, pth in CFG_KEYS.items()}


def b0_off_value(v, off) -> bool:
    """B0 だけのキーの値 v (書いていなければ None) の実効値が、書いていないとき (off = B0_OFF の値) と同じか。整数のキーは整数だけ、
    実数のキーは整数か実数を認める (真偽値は認めない)。"""
    if v is None:
        return True
    if isinstance(v, bool):
        return False
    return (type(v) is int and v == off) if isinstance(off, int) else (type(v) in (int, float) and float(v) == off)


def b0_left(st: dict) -> list:
    """段の設定の値 st ({キー: 値}、eff_settings) のうち、B0 だけのキーで実効値が書いていないときと違うもの (B0 が残っている)。"""
    return [k for k in B0_ONLY if not b0_off_value(st.get(k), B0_OFF[k])]


# --- 生成器が足した明示すべきキー (plan tooling-generator-defaults-and-launch-gate §5.1 #11・#26) ------------------------------------
# #11 以降の生成器は段の config に「明示すべきキー」を旧い実効値で足す (config_builder.finish_config → yaml_strict.insert_scalars:
# 既にある写像へのキーの挿入だけ)。NS は prepare_ns が足した後の設定を書き、足したキーを RECIPE_CONTEXT.json の explicit_added
# {family, prestage, main} に残す。Euler は run_staged が段ごとに起動の直前に足す。登録の検査は足す前の設定で行い (検査の意味は足さない版と
# 同じ)、足したのが記録どおりのキーの挿入だけであることは、除いた設定に生成器の仕上げを当て直して元の文字列に戻ることで確かめる。
EXPLICIT_CONTEXT = "RECIPE_CONTEXT.json"


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
    記録の無い準備に足したキーがあれば、除かずに比べるので検査で不成立になる)。"""
    rec = jload(Path(run) / EXPLICIT_CONTEXT)
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


def explicit_strip(RA, text: str, added, bcond: str, family, tag=None) -> str:
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


def ns_config_checks(RA, run: Path) -> dict:
    """NS の段の実効設定 (ns_stage_configs で導く S1・S2・本段) が登録どおり (h370_problems.py prep-check の want と同じ)、S1・S2 に B0 の
    キー (起動する段の実効値で: 書いていないときの値 B0_OFF の明示は無いのと同じ; b0_left)・本段の output が無い、3 段の implicitRelax 0.7、
    sha が prepare_info の ns_stage_configs と同じ。手元の prep-check の sha は記録だけ。
    生成器が明示すべきキーを足した準備 (RECIPE_CONTEXT.json の explicit_added) は、足す前の設定から導いた段で上の検査をし (stages の値も
    足す前のもの)、起動する段の設定から足したキーを除くと足す前の段と同じことを確かめる。sha は起動する (足した後の) 文字列のもの。"""
    fails = []
    main = (run / "solverConfig.yaml").read_text()
    pre_p = run / RA.PRESTAGE_CONFIG
    if not pre_p.is_file():
        return {"ok": False, "failures": [f"{RA.PRESTAGE_CONFIG} が無い (main_line が効いていない)"]}
    pre = pre_p.read_text()
    try:
        pre_st, main_eff, main_want = RA.ns_stage_configs(main, "full", cfg_pre=pre)
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "failures": [f"ns_stage_configs が例外: {e!r}"]}
    ex = explicit_record(run)
    pre_st_b, main_eff_b = pre_st, main_eff
    if ex:
        bc_p = run / "bcondConfig.yaml"
        try:
            bc = bc_p.read_text() if bc_p.is_file() else None
            pre_b = explicit_strip(RA, pre, ex.get("prestage"), bc, ex.get("family"))
            main_b = explicit_strip(RA, main, ex.get("main"), bc, ex.get("family"))
            pre_st_b, main_eff_b, _ = RA.ns_stage_configs(main_b, "full", cfg_pre=pre_b)
            differ = [t for (t, c, _, _), (_, cb, _, _) in zip(pre_st, pre_st_b) if explicit_remove(c, ex.get("prestage")) != cb]
            differ += ["main"] if explicit_remove(main_eff, ex.get("main")) != main_eff_b else []
        except Exception as e:  # noqa: BLE001
            return {"ok": False, "failures": [f"生成器が足した明示すべきキー ({EXPLICIT_CONTEXT}) を除いて比べられない: {e!r}"],
                    "explicit_added": ex}
        if differ or [t for t, *_ in pre_st] != [t for t, *_ in pre_st_b]:
            fails.append(f"起動する段の設定から足したキーを除いても足す前の段と同じにならない: {differ}")
    st = {tag: eff_settings(RA, cfg) for tag, cfg, _, _ in pre_st_b}
    st["main"] = eff_settings(RA, main_eff_b)
    for tag, w in NS_STAGE_WANT.items():
        for k, v in w.items():
            if st.get(tag, {}).get(k) != v:
                fails.append(f"{tag} の {k} = {st.get(tag, {}).get(k)!r} ({v!r} であること)")
    # S1・S2 の B0 は起動する (足した後の) 段の実効値で見る: 書いていないときの値 B0_OFF の明示は B0 なし、それ以外は B0 あり。
    # 足す前の段は起動する段からキーを除いただけなので、起動する段に B0 が無ければ足す前の段にも無い
    st_run = {tag: eff_settings(RA, cfg) for tag, cfg, _, _ in pre_st}
    for tag in ("S1_soft", "S2_mid"):
        left = b0_left(st_run[tag])
        if left or st_run[tag]["output"] is not None or st[tag]["output"] is not None:
            fails.append(f"{tag} に B0 のキー {left} (書いていないときと違う実効値; output {st_run[tag]['output']}) が残っている")
    if not all(st[t]["implicitRelax"] == 0.7 for t in st):
        fails.append(f"implicitRelax が 3 段で 0.7 でない: {[st[t]['implicitRelax'] for t in st]}")
    sha = {tag: sha256_text(cfg) for tag, cfg, _, _ in pre_st}
    sha["main"] = sha256_text(main_eff)
    sha["prestage_config"] = sha256_text(pre)
    info = jload(run / "prepare_info.json") or {}
    nsc = info.get("ns_stage_configs") or {}
    match = {k: nsc.get(k) == v for k, v in sha.items()}
    if not all(match.values()):
        fails.append(f"prepare_info の ns_stage_configs が導いた段の設定と違う: {match}")
    loc = ((jload(HERE / PREP_CHECK_REC) or {}).get("stage_sha256") or {})
    return {"ok": not fails, "failures": fails, "stages": st, "sha256": sha, "prepare_info_match": match,
            "local_prep_check_same (記録だけ)": {k: loc.get(k) == v for k, v in sha.items()},
            "stage_tags_by_sha": {v: k for k, v in sha.items() if k != "prestage_config"},
            "previous_want": [("S1_soft", NS_PRE, sha["S1_soft"]), ("S2_mid", NS_PRE, sha["S2_mid"])], "explicit_added": ex}


def phys_up_ok(info: dict) -> bool:
    """prepare_info の physical_upstream.requested (`_physical_upstream` の解決済みの値: r_inlet・L_U・L_pipe に x_in・mu・frozen_design・
    source が付く) の 3 値が H (PHYS_UP) と同じか。"""
    req = (info.get("physical_upstream") or {}).get("requested") or {}
    return isinstance(req, dict) and {k: req.get(k) for k in PHYS_UP} == PHYS_UP


def ns_prepare_info_checks(info: dict, pin_rec, prelaunch: bool) -> dict:
    """NS の prepare_info.json: δ_r の CSV・物理壁の上流 (H)・壁温・B0・格子・CFD ピンの照合・停止規則の設定 (正式の準備だけ)。"""
    fails = []
    if info.get("delta_r_csv") != str(HERE / CSV):
        fails.append(f"prepare_info の delta_r_csv {info.get('delta_r_csv')!r} が {HERE / CSV} でない")
    if not phys_up_ok(info):
        fails.append(f"physical_upstream {(info.get('physical_upstream') or {}).get('requested')} が H {PHYS_UP} でない")
    if info.get("wall_thermal") != {"mode": "isothermal", "Tw": TW}:
        fails.append(f"wall_thermal {info.get('wall_thermal')} が 370 K の等温でない")
    ml = {k: (info.get("main_line") or {}).get(k) for k in MAIN_LINE}
    if ml != MAIN_LINE:
        fails.append(f"main_line {ml} が登録 {MAIN_LINE} と違う")
    m = info.get("mesh") or {}
    if (m.get("ni"), m.get("nj")) != (NS_MESH["ni"], NS_MESH["nj"]):
        fails.append(f"NS の格子 {m.get('ni')} × {m.get('nj')} が Gref {NS_MESH} でない")
    match = ((info.get("initial_line") or {}).get("match") or {})
    if not match.get("checked") or match.get("unchecked"):
        fails.append(f"CFD ピンの照合 (initial_line.match) が不完全: {match}")
    sr = info.get("stop_rule") or {}
    rec = {"stop_rule": {k: sr.get(k) for k in ("requested", "use", "cap", "hold_steps", "out_step_interval", "fallback_steps", "euler_ref",
                                                 "window", "min_points", "x_fracs", "reason")},
           "physical_wall_domain": (info.get("physical_wall") or {}).get("domain"), "euler_ref": info.get("euler_ref")}
    if prelaunch:
        if sr.get("requested") or sr.get("use"):
            fails.append("投入前の準備 (停止規則なし) の prepare_info が停止規則を要求している")
    else:
        from forge_design.feedback.deltastar_loop import handoff_extract_settings
        th = {k: sr.get(k) for k in THRESHOLDS}
        if not (sr.get("requested") and sr.get("use") and sr.get("cap") == NS_CAP and sr.get("hold_steps") == NS_HOLD
                and sr.get("out_step_interval") == NS_OUT and th == THRESHOLDS):
            fails.append(f"停止規則の設定が登録と違う: {rec['stop_rule']}")
        ex = {k: v for k, v in (sr.get("extract") or {}).items() if k != "source"}
        if ex != handoff_extract_settings():
            fails.append(f"停止規則の抽出の条件 {ex} が handoff_extract_settings {handoff_extract_settings()} と違う")
        if sr.get("euler_ref") != pin_rec or info.get("euler_ref") != pin_rec:
            fails.append(f"Euler の参照の記録 (stop_rule {sr.get('euler_ref')}、prepare_info {info.get('euler_ref')}) が固定 {pin_rec} と違う")
    return {"ok": not fails, "failures": fails, **rec}


def delta_in_check(ns_run: Path, eu_run: Path) -> dict:
    """変換後の座標の δ_in の照合 (metrics.delta_in_check.check_runs、x ≥ 0 で ≤ 0.001 mm・全列の被覆)。"""
    from forge_design.metrics.delta_in_check import check_runs
    try:
        r = check_runs(ns_run, eu_run, str(HERE / CSV))
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "failures": [f"delta_in_check が判定不能 ({type(e).__name__}: {e})"]}
    jdump(ns_run / "H370_DELTA_IN.json", r)
    j = r["judged"]
    out = {"pass": r["pass"], "covered": r["covered"], "n_uncovered": r["n_uncovered"], "judged": j, "by_region": r["by_region"],
           "tol_mm": r["tol_mm"], "scale_same": r["scale_same"], "coord_dtype": r["coord_dtype"], "file": "H370_DELTA_IN.json",
           "csv_sha256": r["delta_r_csv_sha256"]}
    fails = [] if r["pass"] else [f"δ_in の照合が不合格 (x ≥ 0 の max {j.get('max_abs_mm')} mm > {r['tol_mm']} mm、被覆 {r['covered']})"]
    if r["delta_r_csv_sha256"] != CSV_SHA:
        fails.append("δ_in の照合が読んだ CSV の sha が登録と違う")
    return {"ok": not fails, "failures": fails, **out}


def mesh_identity(a: Path, b: Path) -> dict:
    """2 つの準備の nozzle.h5 の照合 (変換器の HDF5 はバイトで再現しないのでデータセットで比べる): 格子 (/MESH・/CELLS・/PLANES・/VIZMESH・
    /BCONDS の VALUE 以外) と /VALUE/wall_dist はビット一致を要る (投入前の検査の品質・δ_in・壁の距離がこの run にも当てはまること)。
    初期の場 (ほかの VALUE) の違いは記録だけ (この run の初期の場の sha は別に記録する)。"""
    import h5py
    diff, rec_diff, n = [], [], 0

    def need(k):
        top = k.split("/")[0]
        return top in ("MESH", "CELLS", "PLANES", "VIZMESH") or (top == "BCONDS" and "/VALUE/" not in k) or k == "VALUE/wall_dist"
    with h5py.File(a, "r") as fa, h5py.File(b, "r") as fb:
        da, db = {}, {}
        fa.visititems(lambda k, o: da.__setitem__(k, o) if isinstance(o, h5py.Dataset) else None)
        fb.visititems(lambda k, o: db.__setitem__(k, o) if isinstance(o, h5py.Dataset) else None)
        ka, kb = {k for k in da if need(k)}, {k for k in db if need(k)}
        if ka != kb:
            diff.append(f"格子のデータセットの集合が違う ({sorted(ka ^ kb)[:5]})")
        for k in sorted(set(da) & set(db)):
            x, y = np.asarray(da[k][()]), np.asarray(db[k][()])
            same = x.dtype == y.dtype and x.shape == y.shape and x.tobytes() == y.tobytes()
            if not same:
                (diff if need(k) else rec_diff).append(k)
            n += 1
    return {"ok": not diff, "n_datasets": n, "differ": diff[:20], "initial_field_differ (記録だけ)": rec_diff[:20]}


# --- Euler の段の設定 (run_staged と同じ導出) ------------------------------------------------------------------------------------
def euler_stage_configs(RA, cfg_main: str, cfl_main: float, mid_stage: bool, bcond: str | None = None) -> dict:
    """runner_axismach.run_staged と同じ導出で、soft・mid・本段の起動直前の設定の文字列を返す (forge を回さずに sha を決める)。
    生成器が段の config に明示すべきキーを足す版 (run_staged が起動の直前に config_builder.finish_config を当てる; plan
    tooling-generator-defaults-and-launch-gate §5.1 #11) では、同じ仕上げ (機種 axismach_euler、段の tag、run の bcond) を当てた文字列。
    その版で bcond が無ければ ValueError (足すキーは bcond の適用条件で決まる)。"""
    main = RA._cfg_set(cfg_main, {"cfl": cfl_main, "cfl_pseudo": cfl_main})
    out = {"soft": RA._cfg_set(RA._first_order(RA._cfg_set(main, {"cfl": 0.5, "cfl_pseudo": 0.5})), {"nStepOuter": EU_SOFT, "outStepInterval": EU_SOFT})}
    if mid_stage:
        out["mid"] = RA._cfg_set(RA._cfg_set(main, {"cfl": 1.0, "cfl_pseudo": 1.0}), {"nStepOuter": EU_SOFT, "outStepInterval": EU_SOFT})
    out["main"] = main
    fin = generator_finish(RA)
    if fin is not None:
        if bcond is None:
            raise ValueError("生成器が段の config に明示すべきキーを足す版では、段の設定の導出に run の bcondConfig.yaml が要る")
        out = {t: fin(c, bcond, "axismach_euler", tag=t)[0] for t, c in out.items()}
    return out


def euler_prepare_checks(RA, run: Path, info: dict) -> dict:
    """Euler の prepare_info・段の設定: 格子 2000 × 97・throat_refine 5、H の上流、壁の種類、CFD ピンの照合、本段 60,000・出力 1,000・cfl 2.0。"""
    fails = []
    m = info.get("mesh") or {}
    tr = ((m.get("params") or {}).get("throat_refine"))
    if (m.get("ni"), m.get("nj"), tr) != (EU_MESH["ni"], EU_MESH["nj"], EU_MESH["throat_refine"]):
        fails.append(f"Euler の格子 {m.get('ni')} × {m.get('nj')}・throat_refine {tr} が {EU_MESH} でない")
    if not phys_up_ok(info):
        fails.append(f"physical_upstream {(info.get('physical_upstream') or {}).get('requested')} が H {PHYS_UP} でない")
    ew = info.get("euler_wall") or {}
    if ew.get("kind") != "design_with_physical_upstream" or ew.get("delta_r") != 0.0 or not (ew.get("gate") or {}).get("pass"):
        fails.append(f"Euler の参照の壁の記録 (種類 {ew.get('kind')!r}・δ_r {ew.get('delta_r')!r}・ゲート {(ew.get('gate') or {}).get('pass')})")
    match = ((info.get("initial_line") or {}).get("match") or {})
    if not match.get("checked") or match.get("unchecked"):
        fails.append(f"CFD ピンの照合 (initial_line.match) が不完全: {match}")
    p = RA.load_problem(HERE / EU_PROBLEM)
    cfl = float(p.evaluate.get("cfl_main"))
    mid = bool(p.evaluate.get("mid_stage", p.is_semiperfect))
    bc_p = run / "bcondConfig.yaml"
    try:
        cfgs = euler_stage_configs(RA, (run / "solverConfig.yaml").read_text(), cfl, mid, bc_p.read_text() if bc_p.is_file() else None)
    except Exception as e:  # noqa: BLE001
        return {"ok": False, "failures": fails + [f"Euler の段の設定を導けない: {e!r}"], "stages": {}, "sha256": {}, "stage_tags_by_sha": {},
                "cfl_main": cfl, "mid_stage": mid}
    eff ={t: {k: RA._cfg_value(c, k) for k in ("convMethod", "cfl", "cfl_pseudo", "nStepOuter", "outStepInterval")} for t, c in cfgs.items()}
    if (cfl, mid, eff["main"]["nStepOuter"], eff["main"]["outStepInterval"], eff["main"]["cfl"]) != (EU_CFL, True, EU_STEPS, EU_OUT, EU_CFL):
        fails.append(f"Euler の本段の設定が登録と違う (cfl_main {cfl}・mid {mid}・{eff['main']})")
    if eff["soft"]["convMethod"] != 0 or eff["soft"]["nStepOuter"] != EU_SOFT or eff.get("mid", {}).get("nStepOuter") != EU_SOFT:
        fails.append(f"Euler の前の段の設定が想定と違う: {eff}")
    sha = {t: sha256_text(c) for t, c in cfgs.items()}
    return {"ok": not fails, "failures": fails, "stages": eff, "sha256": sha, "stage_tags_by_sha": {v: k for k, v in sha.items()},
            "cfl_main": cfl, "mid_stage": mid, "physical_upstream": info.get("physical_upstream"), "euler_wall_throat": ew.get("throat")}


# --- 観測だけの包み ---------------------------------------------------------------------------------------------------------------
class Observer:
    """run_forge (前の段と Euler の全段) と _write_main_stage_started (NS の本段の直前) を包む。段の設定の sha が導いたものと違えば forge を
    起動させない (例外)。NS の本段の直前は通知の中身・資源を照合し、不成立なら WATCH_ABORT を書いて例外 (本段の forge を起動させない)。
    disk_min: 本段の前の資源の照合のディスクの空きの下限 (既定 DISK_START_MIN; 診断の経路は保存の容量の見積もりとの大きい方)。"""

    def __init__(self, RA, run: Path, tags: dict, main_tag: str, prepared: dict | None = None, disk_min: int = DISK_START_MIN):
        self.RA, self.run, self.tags, self.main_tag, self.prep = RA, run, tags, main_tag, prepared or {}
        self.disk_min = int(disk_min)
        self.orig = {"run_forge": RA.run_forge, "write_main": RA._write_main_stage_started}

    def install(self):
        self.RA.run_forge, self.RA._write_main_stage_started = self.run_forge, self.write_main

    def uninstall(self):
        self.RA.run_forge, self.RA._write_main_stage_started = self.orig["run_forge"], self.orig["write_main"]

    def run_forge(self, run_dir):
        cfg_sha = sha256_file(Path(run_dir) / "solverConfig.yaml")
        tag = self.tags.get(cfg_sha)
        if tag is None:
            jappend(self.run / "h370_stages.jsonl", {"event": "refused", "t": time.time(), "config_sha256": cfg_sha})
            raise RuntimeError(f"段の設定 {cfg_sha[:16]} が準備の時点で導いた段のどれとも違う — forge を起動しない")
        if tag == self.main_tag:
            res = resources("本段の forge の前", self.disk_min)
            jappend(self.run / "h370_stages.jsonl", {"event": "resources", "tag": tag, **res})
            if not res["ok"]:
                (self.run / "RESOURCE_STOP").write_text(now() + "\n" + "\n".join(res["failures"]) + "\n")
                raise RuntimeError("V-h370: 本段の forge の前の資源が足りない: " + "; ".join(res["failures"]))
        t0 = time.time()
        jappend(self.run / "h370_stages.jsonl", {"event": "start", "tag": tag, "t": t0, "config_sha256": cfg_sha})
        log(f"forge 起動 段 {tag} (設定 {cfg_sha[:16]})")
        rc = self.orig["run_forge"](run_dir)
        t1 = time.time()
        jappend(self.run / "h370_stages.jsonl", {"event": "end", "tag": tag, "t0": t0, "t1": t1, "rc": rc, "config_sha256": cfg_sha})
        log(f"forge 終了 段 {tag} rc={rc} ({t1 - t0:.0f} s)")
        return rc

    def write_main(self, run_dir, cfg_main, stages, done, cfg_pre):
        rec = self.orig["write_main"](run_dir, cfg_main, stages, done, cfg_pre)
        run_dir = Path(run_dir)
        fails = []
        if rec.get("res_files_present"):
            fails.append(f"通知の時点で前の段の出力が残っている: {rec['res_files_present']}")
        if rec["main_config"]["sha256"] != self.prep.get("main_sha256"):
            fails.append(f"本段の設定の sha256 {str(rec['main_config']['sha256'])[:16]} が準備の時点で導いた本段 (上限 {NS_CAP}) と違う")
        if rec["main_config"].get("nStepOuter") != NS_CAP:
            fails.append(f"本段の nStepOuter {rec['main_config'].get('nStepOuter')} が上限 {NS_CAP} でない")
        got_prev = [(s.get("tag"), s.get("nsteps"), s.get("config_sha256")) for s in rec.get("previous_stages", [])]
        if got_prev != [tuple(x) for x in self.prep.get("want_prev", [])]:
            fails.append(f"前の段の記録が導いた S1・S2 と違う: {got_prev}")
        if (rec.get("prestage_config") or {}).get("sha256") != self.prep.get("prestage_sha256"):
            fails.append("前の段の元の設定の sha256 が準備の時点と違う")
        res = resources("本段の forge の前", self.disk_min)
        if not res["ok"]:
            fails += ["資源: " + f for f in res["failures"]]
        out = {"t": time.time(), "main_stage_started": rec, "resources": res, "failures": fails, "ok": not fails,
               "initial_field_sha256": rec["initial_field"]["sha256"]}
        jdump(run_dir / "H370_MAIN_START.json", out)
        if fails:
            (run_dir / "WATCH_ABORT").write_text(now() + "\n" + "\n".join(fails) + "\n")
            if any(f.startswith("資源: ") for f in fails):
                (run_dir / "RESOURCE_STOP").write_text(now() + "\n" + "\n".join(res["failures"]) + "\n")
            log("本段の直前の照合が不成立 — 本段の forge を起動させない: " + "; ".join(fails))
            raise RuntimeError("V-h370: 本段の直前の照合が不成立 (WATCH_ABORT)")
        log(f"本段の通知を確かめた (初期の場 {str(rec['initial_field']['sha256'])[:16]}、上限 {NS_CAP}) — 見張りは stop_rule.watch_main_stage")
        return rec


def fail_record(run: Path, pc: dict, failures: list) -> None:
    """前提の不成立を run に残す (forge は回さない)。"""
    run.mkdir(parents=True, exist_ok=True)
    pc["failures"] = list(pc.get("failures", [])) + list(failures)
    pc["ok"] = False
    jdump(run / "H370_PRECHECK.json", pc)
    (run / "PRECHECK_FAILED").write_text(now() + "\n" + "\n".join(pc["failures"]) + "\n")


def collect_failures(items: dict) -> list:
    return [f for it in items.values() if isinstance(it, dict) for f in it.get("failures", [])]


# --- (i) Euler -------------------------------------------------------------------------------------------------------------------
def phase_euler(a, RA, code: dict) -> int:
    dry = a.dry_run
    run = HERE / (f"_h370_dry_{dry}_euler" if dry else EU_RUN)
    ns_pre = HERE / (f"_h370_dry_{dry}_ns" if dry else NS_PRELAUNCH)
    for d in (run, ns_pre):
        if d.exists():
            die(f"{d} が既にある (既存の run・準備は消さない。中身を確かめてから退ける)")
    pc = {"plan": PLAN, "phase": "euler", "created": now(), "run": run.name, "dry_run": dry, "code": code, "failures": [], "ok": False}
    pc["inputs"] = input_checks()
    res0 = resources("起動の前")
    pc["resources_start"] = res0
    log(f"資源 {json.dumps(res0, ensure_ascii=False)}")
    if not pc["inputs"]["ok"]:
        die("入力が欠けているか登録と違う: " + "; ".join(pc["inputs"]["failures"]))
    if not res0["ok"] and not dry:
        die("資源が足りない (forge を回さない): " + "; ".join(res0["failures"]))
    timing = Timing(run)
    if not dry:
        timing.start()
    t_start = time.time()
    log(f"Euler の prepare ({EU_PROBLEM} → {run.name})")
    try:
        info = RA.prepare(HERE / EU_PROBLEM, run)
        info["stages"] = "full"                       # runner_axismach.main と同じ記録
        (run / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    except Exception as e:  # noqa: BLE001
        log("Euler の prepare が失敗:\n" + traceback.format_exc())
        if run.is_dir():
            fail_record(run, pc, [f"Euler の prepare が失敗: {e!r}"])
        timing.stop_evt.set()
        return 2
    (run / "START_EPOCH").write_text(f"{t_start:.6f}\n")
    items = {}
    items["E1_wall_restore"] = euler_wall_restore_check(RA, run)
    m_eu = RA.mesh_euler_block(RA.load_problem(HERE / EU_PROBLEM))
    items["E2_quality"] = quality_check(RA, run, int(m_eu.get("ar_max", 1000)))
    items["E3_prepare"] = euler_prepare_checks(RA, run, info)
    items["E4_sha"] = {"nozzle_h5_initial": sha256_file(run / "nozzle.h5"), "nozzle_msh": sha256_file(run / "nozzle.msh"),
                       "bcondConfig": sha256_file(run / "bcondConfig.yaml"), "solverConfig_main": sha256_file(run / "solverConfig.yaml"),
                       "problem": sha256_file(HERE / EU_PROBLEM), "mesh_hashes": (info.get("mesh") or {}).get("hashes")}
    log(f"Euler の検査: 壁 {items['E1_wall_restore']['ok']}、品質 {items['E2_quality'].get('check_mesh_quality_nozzle_h5')}、"
        f"設定 {items['E3_prepare']['ok']}")
    # NS の投入前の準備 (停止規則なし): 格子・δ_in・壁の出力・設定
    log(f"NS の投入前の準備 ({NS_PROBLEM} → {ns_pre.name}、停止規則なし・Euler の参照なし; forge は回さない)")
    try:
        ns_info = RA.prepare_ns(HERE / NS_PROBLEM, ns_pre, delta_r_csv=str(HERE / CSV), offset="radial", euler_ref=None, stop_rule=False)
    except Exception as e:  # noqa: BLE001
        tb = traceback.format_exc()
        if isinstance(e, subprocess.CalledProcessError):
            tb += f"\n--- stdout ---\n{(e.stdout or '')[-3000:]}\n--- stderr ---\n{(e.stderr or '')[-3000:]}"
        log("NS の投入前の準備が失敗:\n" + tb)
        pc["items"] = items
        fail_record(run, pc, [f"NS の投入前の準備 (prepare_ns) が失敗: {e!r}"])
        timing.stop_evt.set()
        return 2
    items["N1_quality"] = quality_check(RA, ns_pre, 5000)
    items["N2_delta_in"] = delta_in_check(ns_pre, run)
    items["N3_wall_output"] = wall_output_check(ns_pre)
    items["N4_config"] = ns_config_checks(RA, ns_pre)
    items["N5_prepare"] = ns_prepare_info_checks(ns_info, None, prelaunch=True)
    loc = jload(HERE / PREP_CHECK_REC) or {}
    gsha = ((ns_info.get("mesh") or {}).get("hashes") or {}).get("generator_sha256")
    items["N6_generator (記録だけ)"] = {"generator_sha256": gsha, "local_prep_check_same": (loc.get("ns_mesh") or {}).get("generator_sha256") == gsha}
    pc["items"] = items
    pc["t_checked"] = time.time()
    fails = collect_failures(items)
    pc["failures"], pc["ok"] = fails, not fails
    jdump(run / "H370_PRECHECK.json", pc)
    log(f"Euler の投入前の検査 {'合格' if pc['ok'] else '不成立'}: 品質 Euler {items['E2_quality'].get('check_mesh_quality_nozzle_h5')} / NS "
        f"{items['N1_quality'].get('check_mesh_quality_nozzle_h5')}、δ_in max {items['N2_delta_in'].get('judged', {}).get('max_abs_mm')} mm、"
        f"壁の出力 {items['N3_wall_output']['ok']}" + ("" if pc["ok"] else " — " + "; ".join(fails)))
    if not items["N3_wall_output"]["ok"]:
        (run / "PRECHECK_FAILED").write_text(now() + "\n" + "\n".join(fails) + "\n")
        log("NS が停止の出力で壁の出力を書かない設定 — 書く設定を足すのは方針の変更なので止めて親に返す (終了コード 4)")
        timing.stop_evt.set()
        return EXIT_WALL_OUTPUT
    if dry:
        log(f"--dry-run {dry}: Euler と NS の投入前の検査 {'合格' if pc['ok'] else '不成立'} (forge は回さない)")
        return 0 if pc["ok"] else 2
    if not pc["ok"]:
        (run / "PRECHECK_FAILED").write_text(now() + "\n" + "\n".join(fails) + "\n")
        log("前提の検査が不成立 — Euler の forge を回さない")
        timing.stop_evt.set()
        return 2
    # --- 実行 (runner_axismach.main と同じ: run_staged → collect) ---
    obs = Observer(RA, run, items["E3_prepare"]["stage_tags_by_sha"], "main")
    obs.install()
    run_end = {"t_start": t_start, "phase": "euler"}
    try:
        rc = RA.run_staged(run, cfl_main=items["E3_prepare"]["cfl_main"], mid_stage=items["E3_prepare"]["mid_stage"], stages="full")
        rc_txt = str(int(rc))
    except Exception:  # noqa: BLE001
        run_end["exception"] = traceback.format_exc()
        log("run_staged が例外:\n" + run_end["exception"])
        rc_txt = "EXC"
    finally:
        obs.uninstall()
    (run / "RUN_RC").write_text(rc_txt + "\n")
    timing.stop_evt.set()
    timing.join(timeout=60)
    if rc_txt == "0":
        try:
            run_end["collect"] = RA.collect(HERE / EU_PROBLEM, run)
        except Exception as e:  # noqa: BLE001 — 記録 (metrics.json) の失敗は資格の判定を止めない
            run_end["collect_error"] = repr(e)
    run_end.update(t_end=time.time(), run_rc=rc_txt, last_res=(res_steps(run)[-1].name if res_steps(run) else None))
    jdump(run / "RUN_END.json", run_end)
    log(f"{run.name} 終了 RUN_RC={rc_txt} 最後の出力 {run_end['last_res']} ({run_end['t_end'] - t_start:.0f} s)")
    # --- 予備の資格の判定 (h370_judge.py --euler-prequalify) ---
    q = subprocess.run([sys.executable, str(HERE / "h370_judge.py"), "--euler-prequalify", "--design", str(WT / "design"),
                        "--tools", str(RA.FORGE_TOOLS)], cwd=str(HERE), capture_output=True, text=True)
    print(q.stdout[-8000:] + q.stderr[-3000:], flush=True)
    log(f"Euler の予備の資格の判定 rc={q.returncode} ({'満たす → 固定 ' + PIN if q.returncode == 0 else '満たさない・判定不能 → NS に進まない'})")
    return 0 if q.returncode == 0 else 2


# --- (ii) NS ----------------------------------------------------------------------------------------------------------------------
def read_pin() -> dict:
    """固定した Euler の参照 (h370_judge.py --euler-prequalify が書く) を読み、sha を照合する。σ_ref の代わりの参照 (予備の資格の
    sigma_ref_snapshots、§5.1 #36g (2)・#36i (M1)) が末尾の窓の全 20 出力 (SIGMA_REF_STEPS、登録の定数から列挙) ちょうどで、各出力が在り
    sha が記録と同じことを NS を回す前に確かめる (欠ければ σ_ref を測れず、NS の後で更新に進めない)。戻り {spec, record, pin}。不成立は die。"""
    pin, qual = jload(HERE / PIN), jload(HERE / QUAL)
    if not pin or not qual:
        die(f"{PIN} か {QUAL} が無い (Euler の予備の資格の判定が済んでいない)")
    if not qual.get("prequalified") or (qual.get("pin") or {}).get("sha256") != pin.get("sha256") or pin.get("run") != EU_RUN:
        die(f"Euler の予備の資格の記録と固定が合わない (prequalified {qual.get('prequalified')}、run {pin.get('run')})")
    f = HERE / pin["run"] / pin["res"]
    if not f.is_file() or sha256_file(f) != pin["sha256"]:
        die(f"固定した Euler の参照 {f} が無いか sha256 が {pin['sha256'][:16]}… と違う")
    if pin.get("res") != f"res_{EU_STEPS}.h5":
        die(f"固定した Euler の参照 {pin.get('res')} が本段の上限の出力 res_{EU_STEPS}.h5 でない")
    alts = qual.get("sigma_ref_snapshots") or []
    if not alts:
        die(f"{QUAL} に σ_ref の代わりの参照 (sigma_ref_snapshots) が無い")
    got = sorted(a["step"] for a in alts if isinstance(a, dict) and isinstance(a.get("step"), int))
    if got != SIGMA_REF_STEPS or len(got) != len(alts) or any(a.get("res") != f"res_{a.get('step')}.h5" for a in alts):
        miss = [k for k in SIGMA_REF_STEPS if k not in got]
        die(f"{QUAL} の σ_ref の代わりの参照が末尾の窓の全 {len(SIGMA_REF_STEPS)} 出力 ({SIGMA_REF_STEPS[0]}〜{SIGMA_REF_STEPS[-1]}) でない "
            f"(無い {miss[:8]}、記録 {len(got)} 出力)")
    for a in alts:
        g = HERE / pin["run"] / str(a.get("res"))
        if not g.is_file() or sha256_file(g) != a.get("sha256"):
            die(f"σ_ref の代わりの参照 {g} が無いか sha256 が予備の資格の記録と違う (Euler の出力が消された・書き換えられた)")
    spec = {"run": str(HERE / pin["run"]), "res": pin["res"], "sha256": pin["sha256"]}
    return {"spec": spec, "record": spec, "pin": pin}


def phase_ns(a, RA, code: dict) -> int:
    dry = a.dry_run
    pinr = read_pin()
    run = HERE / (f"_h370_dry_{dry}_ns2" if dry else NS_RUN)
    if run.exists():
        die(f"{run} が既にある (既存の run は消さない)")
    pc = {"plan": PLAN, "phase": "ns", "created": now(), "run": run.name, "dry_run": dry, "code": code, "euler_ref": pinr["spec"],
          "failures": [], "ok": False}
    pc["inputs"] = input_checks()
    res0 = resources("起動の前")
    pc["resources_start"] = res0
    log(f"資源 {json.dumps(res0, ensure_ascii=False)}")
    if not pc["inputs"]["ok"]:
        die("入力が欠けているか登録と違う: " + "; ".join(pc["inputs"]["failures"]))
    if not res0["ok"] and not dry:
        die("資源が足りない (forge を回さない): " + "; ".join(res0["failures"]))
    timing = Timing(run)
    if not dry:
        timing.start()
    t_start = time.time()
    log(f"NS の prepare_ns ({NS_PROBLEM} → {run.name}、euler_ref {pinr['spec']['res']} {pinr['spec']['sha256'][:16]}…)")
    try:
        info = RA.prepare_ns(HERE / NS_PROBLEM, run, delta_r_csv=str(HERE / CSV), offset="radial", euler_ref=pinr["spec"])
    except Exception as e:  # noqa: BLE001
        tb = traceback.format_exc()
        if isinstance(e, subprocess.CalledProcessError):
            tb += f"\n--- stdout ---\n{(e.stdout or '')[-3000:]}\n--- stderr ---\n{(e.stderr or '')[-3000:]}"
        log("prepare_ns が失敗:\n" + tb)
        if run.is_dir():
            fail_record(run, pc, [f"prepare_ns が失敗: {e!r}"])
        timing.stop_evt.set()
        return 2
    (run / "START_EPOCH").write_text(f"{t_start:.6f}\n")
    items = {}
    items["N1_quality"] = quality_check(RA, run, 5000)
    items["N2_delta_in"] = delta_in_check(run, HERE / EU_RUN)
    items["N3_wall_output"] = wall_output_check(run)
    items["N4_config"] = ns_config_checks(RA, run)
    items["N5_prepare"] = ns_prepare_info_checks(info, pinr["record"], prelaunch=False)
    pre = HERE / NS_PRELAUNCH
    if pre.is_dir():
        same = {n: (run / n).read_bytes() == (pre / n).read_bytes() for n in ("solverConfig.yaml", RA.PRESTAGE_CONFIG, "bcondConfig.yaml",
                                                                             "nozzle.msh", "probe.yaml", "species_meta.yaml")
                if (run / n).is_file() and (pre / n).is_file()}
        mi = mesh_identity(run / "nozzle.h5", pre / "nozzle.h5")
        fails7 = [f"(i) の投入前の準備と {n} が違う" for n, v in same.items() if not v] + ([] if mi["ok"] else [f"(i) の投入前の準備と nozzle.h5 が違う: {mi['differ']}"])
        items["N7_same_as_prelaunch"] = {"ok": not fails7, "failures": fails7, "files": same, "nozzle_h5": mi, "prelaunch": NS_PRELAUNCH}
    else:
        items["N7_same_as_prelaunch"] = {"ok": False, "failures": [f"(i) の投入前の準備 {NS_PRELAUNCH} が無い (照合できない)"]}
    items["N8_sha"] = {"nozzle_h5_initial": sha256_file(run / "nozzle.h5"), "nozzle_msh": sha256_file(run / "nozzle.msh"),
                       "bcondConfig": sha256_file(run / "bcondConfig.yaml"), "solverConfig_main": sha256_file(run / "solverConfig.yaml"),
                       "solverConfig_prestage": sha256_file(run / RA.PRESTAGE_CONFIG), "problem": sha256_file(HERE / NS_PROBLEM),
                       "euler_ref": pinr["spec"], "mesh_hashes": (info.get("mesh") or {}).get("hashes")}
    pc["items"] = items
    pc["t_checked"] = time.time()
    fails = collect_failures(items)
    pc["failures"], pc["ok"] = fails, not fails
    jdump(run / "H370_PRECHECK.json", pc)
    log(f"NS の投入前の検査 {'合格' if pc['ok'] else '不成立'}: 品質 {items['N1_quality'].get('check_mesh_quality_nozzle_h5')}、"
        f"δ_in max {items['N2_delta_in'].get('judged', {}).get('max_abs_mm')} mm、壁の出力 {items['N3_wall_output']['ok']}、"
        f"(i) の準備と同一 {items['N7_same_as_prelaunch']['ok']}" + ("" if pc["ok"] else " — " + "; ".join(fails)))
    if not items["N3_wall_output"]["ok"]:
        (run / "PRECHECK_FAILED").write_text(now() + "\n" + "\n".join(fails) + "\n")
        log("NS が停止の出力で壁の出力を書かない設定 — 止めて親に返す (終了コード 4)")
        timing.stop_evt.set()
        return EXIT_WALL_OUTPUT
    if dry:
        log(f"--dry-run {dry}: NS の投入前の検査 {'合格' if pc['ok'] else '不成立'} (forge は回さない)")
        return 0 if pc["ok"] else 2
    if not pc["ok"]:
        (run / "PRECHECK_FAILED").write_text(now() + "\n" + "\n".join(fails) + "\n")
        log("前提の検査が不成立 — NS の forge を回さない")
        timing.stop_evt.set()
        return 2
    # --- 実行: run_staged_ns (S1 → S2 → 本段; 本段は停止規則の見張り) ---
    st4 = items["N4_config"]
    main_cap = RA._cfg_set(RA.ns_stage_configs((run / "solverConfig.yaml").read_text(), "full",
                                                cfg_pre=(run / RA.PRESTAGE_CONFIG).read_text())[1], {"nStepOuter": NS_CAP})
    prepared = {"main_sha256": sha256_text(main_cap), "prestage_sha256": st4["sha256"]["prestage_config"], "want_prev": st4["previous_want"]}
    obs = Observer(RA, run, {k: v for k, v in st4["stage_tags_by_sha"].items() if v != "main"}, "main", prepared)
    obs.install()
    run_end = {"t_start": t_start, "phase": "ns"}
    try:
        rc = RA.run_staged_ns(run, stages="full")
        rc_txt = str(int(rc))
        run_end["stop_rule"] = getattr(rc, "stop_rule", None)
    except Exception:  # noqa: BLE001
        run_end["exception"] = traceback.format_exc()
        log("run_staged_ns が例外:\n" + run_end["exception"])
        rc_txt = "EXC"
    finally:
        obs.uninstall()
    (run / "RUN_RC").write_text(rc_txt + "\n")
    timing.stop_evt.set()
    timing.join(timeout=60)
    run_end.update(t_end=time.time(), run_rc=rc_txt, main_started=(run / RA.MAIN_STAGE_STARTED).is_file())
    jdump(run / "RUN_END.json", run_end)
    log(f"{run.name} 終了 RUN_RC={rc_txt} 停止の状態 {(run_end.get('stop_rule') or {}).get('status')} "
        f"({run_end['t_end'] - t_start:.0f} s)")
    if not a.no_judge:
        j = subprocess.run([sys.executable, str(HERE / "h370_judge.py"), "--design", str(WT / "design"), "--tools", str(RA.FORGE_TOOLS)],
                           cwd=str(HERE), capture_output=True, text=True)
        print(j.stdout[-8000:] + j.stderr[-3000:], flush=True)
        log(f"判定 h370_judge.py rc={j.returncode}")
    return 0


# --- 診断の経路 (--phase ns-diag; §5.1 #36l・§6 V-h370-diag) ---------------------------------------------------------------------
def formal_records_state() -> dict:
    """正式の記録 (予備の資格 QUAL・固定 PIN) の有無と sha256 (読むだけ。診断の経路はこれを書かない)。"""
    out = {}
    for name, rel in (("prequalify", QUAL), ("pin", PIN)):
        p = HERE / rel
        out[name] = {"file": rel, "exists": p.is_file(), "sha256": sha256_file(p) if p.is_file() else None}
    return out


def module_versions() -> dict:
    """抽出器・見張り・準備の版: 読み込んだ forge_design の VERSION_FILES の sha256 (投入の時点に記録し、感度の判定が照合する)。"""
    import forge_design
    root = Path(forge_design.__file__).resolve().parent
    return {n: (sha256_file(root / n) if (root / n).is_file() else None) for n in VERSION_FILES}


def diag_euler_ref() -> dict:
    """診断専用の Euler の参照の固定 (§6 V-h370-diag「状態と入力の固定」): run_0950 の res_60000.h5 を**未適格 (予備の資格 FAIL) のまま**
    sha で固定し、感度の判定の腕 A・B (run_0952・run_0953 の res_40000.h5) の sha も固定して DIAG_REF に書く (無ければ O_EXCL で作る、
    あれば今の sha と照合して違えば止める; 消さない・書き換えない)。正式の記録 (QUAL・PIN) は**読むだけ**: 予備の資格の記録があって
    prequalified が偽であること (資格があるなら正式の --phase ns を使う)、正式の固定 PIN が無いこと、予備の資格の記録の末尾の出力の sha
    (tail_outputs) に res_60000 があれば一致することを確かめる。read_pin の資格の検査を迂回するのはこの関数を呼ぶ診断の経路の中だけ。
    戻り {"spec" (prepare_ns の euler_ref = {run (絶対パス), res, sha256}), "record", "formal_before"}。不成立は die。"""
    formal = formal_records_state()
    qual = jload(HERE / QUAL)
    if qual is None:
        die(f"{QUAL} が無い (予備の資格の判定が済んでいない) — 診断の経路は予備の資格 FAIL の参照を診断専用に使う経路")
    if qual.get("prequalified") is not False:
        die(f"{QUAL} の prequalified が偽でない ({qual.get('prequalified')!r}) — 資格があるなら正式の --phase ns を使う (診断の経路は未適格の参照のため)")
    if formal["pin"]["exists"]:
        die(f"正式の固定 {PIN} がある (予備の資格 FAIL なら書かれない) — 記録が食い違う。中身を確かめてから進める")
    f = HERE / EU_RUN / f"res_{EU_STEPS}.h5"
    if not f.is_file():
        die(f"{f} が無い")
    sha = sha256_file(f)
    tail = {o.get("step"): o.get("sha256") for o in (qual.get("tail_outputs") or []) if isinstance(o, dict)}
    cross = tail.get(EU_STEPS)
    if cross is not None and cross != sha:
        die(f"{EU_RUN}/res_{EU_STEPS}.h5 の sha256 {sha[:16]} が予備の資格の記録 (tail_outputs) の {str(cross)[:16]} と違う")
    sens = {}
    for arm, (rn, res) in DIAG_SENS_REFS.items():
        g = HERE / rn / res
        if not g.is_file():
            die(f"感度の判定の腕 {arm} の参照 {rn}/{res} が無い (§6 V-h370-diag の入力)")
        sens[arm] = {"run": rn, "res": res, "sha256": sha256_file(g)}
    eref = {"run": EU_RUN, "res": f"res_{EU_STEPS}.h5", "step": EU_STEPS, "sha256": sha}
    rec = {"plan": PLAN_DIAG, "created": now(), "kind": DIAG_STATUS, "qualified": False, "diagnostic_only": True,
           "update_conditions_met": False, "euler_ref": eref, "sens_refs": sens,
           "prequalify_record": {"file": QUAL, "sha256": formal["prequalify"]["sha256"], "prequalified": qual.get("prequalified"),
                                 "rc": qual.get("rc"), "failures": (qual.get("failures") or [])[:4], "undecidable": (qual.get("undecidable") or [])[:4],
                                 "tail_res_sha256": cross, "cross_check": ("一致" if cross is not None else "予備の資格の記録に無い (照合なし)")},
           "formal_pin": formal["pin"],
           "note": ("正式の記録 (予備の資格・固定) は読むだけで書き換えない。この記録は診断の NS (run_0951) の停止規則の帯・δ の記録・投入前の照合と、"
                    "感度の判定 (h370_judge.py --diag-sens) が使う")}
    p = HERE / DIAG_REF
    old = jload(p)
    if old is not None:
        if (old.get("euler_ref"), old.get("sens_refs")) != (eref, sens):
            die(f"{DIAG_REF} が既にあり、今の出力の sha と違う (記録 {old.get('euler_ref')} {old.get('sens_refs')}) — 消さない")
        used = old
    else:
        p.parent.mkdir(parents=True, exist_ok=True)
        fd = os.open(str(p), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.write(fd, json.dumps(rec, indent=1, ensure_ascii=False).encode())
        os.close(fd)
        used = rec
    return {"spec": {"run": str(HERE / EU_RUN), "res": f"res_{EU_STEPS}.h5", "sha256": sha}, "record": used, "formal_before": formal}


def disk_free() -> int:
    return shutil.disk_usage(HERE).free


def capacity_check(run: Path, pre: Path) -> dict:
    """保存の容量の照合 (§6 V-h370-diag「保存」): 本段の全期間の場と境界の出力を残すので、上限まで回したときの出力
    (DIAG_N_OUT = 80 + S1・S2 の DIAG_EXTRA_OUT) × (1 出力の大きさ + 境界の出力) × 余裕 DIAG_MARGIN + 走行中の停止の閾値 DISK_MIN を、
    今のディスクの空きと比べる。1 出力の大きさは投入前の準備 (pre、停止規則なし) の nozzle.h5 と、この run の nozzle.h5 の大きい方
    (nozzle.h5 は格子と初期の場を持つので res_<step>.h5 より大きい = 多めの見積もり)。足りなければ不成立 (forge を回さない)。"""
    sizes = {d.name: (d / "nozzle.h5").stat().st_size for d in (pre, run) if (d / "nozzle.h5").is_file()}
    if pre.name not in sizes:
        return {"ok": False, "failures": [f"1 出力の大きさを見積もる投入前の準備の nozzle.h5 ({pre.name}) が無い"], "nozzle_h5_bytes": sizes}
    per = max(sizes.values())
    wall = per * DIAG_WALL_FRAC
    n = DIAG_N_OUT + DIAG_EXTRA_OUT
    need = int(math.ceil(n * (per + wall) * DIAG_MARGIN)) + DISK_MIN
    free = disk_free()
    rec = {"nozzle_h5_bytes": sizes, "per_output_bytes": per, "boundary_per_output_bytes": int(wall), "n_outputs": n,
           "margin": DIAG_MARGIN, "keep_free_bytes": DISK_MIN, "need_bytes": need, "free_bytes": free,
           "need_gib": round(need / 2 ** 30, 2), "free_gib": round(free / 2 ** 30, 2),
           "rule": (f"({DIAG_N_OUT} + {DIAG_EXTRA_OUT}) 出力 × (nozzle.h5 + 境界の出力 {DIAG_WALL_FRAC:g} 倍) × {DIAG_MARGIN} + 走行中の停止の閾値 "
                    f"{DISK_MIN / 2 ** 30:.0f} GiB ≤ ディスクの空き")}
    fails = [] if free >= need else [f"保存の容量が足りない: 要る {need / 2**30:.2f} GiB > 空き {free / 2**30:.2f} GiB ({rec['rule']})"]
    return {"ok": not fails, "failures": fails, **rec}


def files_sha256(run: Path, skip=()) -> dict:
    """run の全ファイル (下のディレクトリを含む、symlink は辿らない) の {相対パス: sha256} と大きさの合計。"""
    out, total = {}, 0
    for p in sorted(run.rglob("*")):
        rel = str(p.relative_to(run))
        if p.is_symlink() or not p.is_file() or rel in skip:
            continue
        out[rel] = sha256_file(p)
        total += p.stat().st_size
    return {"n": len(out), "bytes": total, "files": out}


def write_diag_mark(run: Path, dref: dict) -> None:
    """診断専用の run の印 (h370_judge.py の正式の判定はこの印のある run を判定せず、受け渡しの関門を呼ばない)。"""
    if run.is_dir() and not (run / DIAG_MARK).exists():
        jdump(run / DIAG_MARK, {"plan": PLAN_DIAG, "created": now(), "diagnostic_only": True, "update_conditions_met": False,
                                "euler_ref": dref["spec"], "euler_ref_status": DIAG_STATUS, "diag_euler_ref_record": DIAG_REF,
                                "forbidden": ["受け渡しの関門 (deltastar_loop.handoff_gate)", "壁の更新", "次の pass の自動実行",
                                              "診断の成果物を正式の受け渡しの表として公開すること", "#36h の上書き"],
                                "sensitivity_judge": "python3 h370_judge.py --diag-sens"})


def diag_end_record(RA, run: Path, rc_txt: str, run_end: dict, dref: dict, wall_dump_name) -> dict:
    """診断の NS の終了の記録 (§6 V-h370-diag「実行・終了」「保存」): 終了の状態 REACHED・CENSORED・実行の失敗、停止の時刻 s (CENSORED は
    最後の出力)、判定の窓 [s − 40,000, s] の 17 出力 (CENSORED は「打切りの時点の診断の窓」)、本段の出力 (場・境界の出力) の完全性、
    停止規則の到達・中断・保持の記録、check_convergence.py --segment --drop 3 の VERDICT (記録)、run の全ファイルの sha の一覧
    (DIAG_SHA)、正式の記録が書き換わっていないこと。何も消さない。"""
    from forge_design.evaluate import stop_rule as SR
    state = SR.read_state(run) or {}
    rec = jload(run / SR.RECORD) or {}
    st = state.get("status")
    if rc_txt == "EXC" or not state:
        kind = "実行の失敗"
    elif st in ("REACHED", "CENSORED"):
        kind = st
    else:
        kind = "実行の失敗"
    out = {"plan": PLAN_DIAG, "created": now(), "diagnostic_only": True, "update_conditions_met": False, "run_rc": rc_txt,
           "end_state": kind, "stop_rule_status": st, "stop_rule_state": state, "exception": (run_end.get("exception") or "")[-2000:] or None}
    s = state.get("stop_step") if kind in ("REACHED", "CENSORED") else None
    out["stop_step"] = s
    if s is not None:
        req = SR.required_window_steps(int(s), int(THRESHOLDS["window"]), NS_OUT)
        out["window"] = {"steps": req, "range": [int(s) - 2 * int(THRESHOLDS["window"]), int(s)], "n": len(req),
                         "kind": ("停止の時刻の窓 [s − 40,000, s] (2 窓 × 20,000、出力 2,500 ごと)" if kind == "REACHED" else
                                  "打切りの時点の診断の窓 (上限の打切り CENSORED: 最後の出力 s で終わる末尾の 17 出力; §6 V-h370-diag で先に登録)"),
                         "missing_res": [k for k in req if not (run / f"res_{k}.h5").is_file()],
                         "missing_wall": ([k for k in req if not (run / wall_dump_name.replace("<step>", str(k))).is_file()] if wall_dump_name else None)}
    steps = sorted(int(m.group(1)) for f in run.iterdir() if (m := re.fullmatch(r"res_(\d+)\.h5", f.name)) and int(m.group(1)) > 0)
    last = steps[-1] if steps else None
    want = list(range(NS_OUT, (last or 0) + 1, NS_OUT))
    out["main_outputs"] = {"n": len(steps), "first": steps[0] if steps else None, "last": last,
                           "gaps": [k for k in want if k not in steps], "off_grid": [k for k in steps if k % NS_OUT],
                           "boundary_dumps": len([f for f in run.iterdir() if re.fullmatch(r"res_(?!nan_).+_\d+_\d+\.h5", f.name)]),
                           "complete": bool(steps) and steps == want}
    out["stop_rule"] = {"hold": rec.get("hold"), "deleted_after_stop": rec.get("deleted_after_stop"), "stop_sha256": rec.get("stop_sha256"),
                        "band": {k: (rec.get("band") or {}).get(k) for k in ("step", "res_sha256", "euler_ref", "out_dir")},
                        "n_rows": len(rec.get("rows") or []), "record": SR.RECORD,
                        "note": "停止規則の全系列 (行)・到達・中断・保持の記録・参照から作った帯 (stop_rule_band/) は run に残す"}
    try:
        cmd = [sys.executable, str(RA.FORGE_TOOLS / "check_convergence.py"), "--segment", "--drop", "3", str(run)]
        cv = subprocess.run(cmd, capture_output=True, text=True)
        (run / DIAG_CONV).write_text(cv.stdout + cv.stderr)
        out["convergence (記録)"] = {"cmd": " ".join(cmd[1:]), "rc": cv.returncode, "file": DIAG_CONV,
                                     "verdict": [ln.strip() for ln in cv.stdout.splitlines() if ln.startswith("===") or ln.startswith("OVERALL")]}
    except Exception as e:  # noqa: BLE001 — 記録だけ
        out["convergence (記録)"] = {"error": repr(e)}
    out["residual_png"] = (run / "residual_history.png").is_file()
    after = formal_records_state()
    out["formal_records"] = {"before": dref["formal_before"], "after": after, "unchanged": after == dref["formal_before"]}
    shal = files_sha256(run, skip=(DIAG_SHA, DIAG_END))
    jdump(run / DIAG_SHA, {"plan": PLAN_DIAG, "created": now(), **shal})
    out["sha_list"] = {"file": DIAG_SHA, "n": shal["n"], "bytes": shal["bytes"]}
    return out


def phase_ns_diag(a, RA, code: dict) -> int:
    """--phase ns-diag (§6 V-h370-diag): 上の docstring。正式の (ii) と同じ準備・検査・本段の直前の照合・見張りで、Euler の参照だけが
    診断専用の固定 (diag_euler_ref)。関門・正式の判定・壁の更新・次の pass は呼ばない。"""
    dry = a.dry_run
    dref = diag_euler_ref()
    run = HERE / (f"_h370_dry_{dry}_nsdiag" if dry else NS_RUN)
    if run.exists():
        die(f"{run} が既にある (既存の run は消さない; 診断の run 名は登録の {NS_RUN} のまま、未使用であること)")
    pc = {"plan": PLAN_DIAG, "phase": "ns-diag", "diagnostic_only": True, "update_conditions_met": False, "created": now(), "run": run.name,
          "dry_run": dry, "code": code, "euler_ref": dref["spec"], "euler_ref_status": DIAG_STATUS, "diag_euler_ref_record": DIAG_REF,
          "sens_refs": dref["record"].get("sens_refs"), "formal_records_before": dref["formal_before"], "versions": module_versions(),
          "failures": [], "ok": False}
    pc["inputs"] = input_checks()
    res0 = resources("起動の前")
    pc["resources_start"] = res0
    log(f"資源 {json.dumps(res0, ensure_ascii=False)}")
    if not pc["inputs"]["ok"]:
        die("入力が欠けているか登録と違う: " + "; ".join(pc["inputs"]["failures"]))
    if not res0["ok"] and not dry:
        die("資源が足りない (forge を回さない): " + "; ".join(res0["failures"]))
    timing = Timing(run)
    if not dry:
        timing.start()
    t_start = time.time()
    log(f"診断の NS の prepare_ns ({NS_PROBLEM} → {run.name}、euler_ref {dref['spec']['res']} {dref['spec']['sha256'][:16]}… = {DIAG_STATUS})")
    try:
        info = RA.prepare_ns(HERE / NS_PROBLEM, run, delta_r_csv=str(HERE / CSV), offset="radial", euler_ref=dref["spec"])
    except Exception as e:  # noqa: BLE001
        tb = traceback.format_exc()
        if isinstance(e, subprocess.CalledProcessError):
            tb += f"\n--- stdout ---\n{(e.stdout or '')[-3000:]}\n--- stderr ---\n{(e.stderr or '')[-3000:]}"
        log("prepare_ns が失敗:\n" + tb)
        if run.is_dir():
            write_diag_mark(run, dref)
            fail_record(run, pc, [f"prepare_ns が失敗: {e!r}"])
        timing.stop_evt.set()
        return 2
    write_diag_mark(run, dref)
    (run / "START_EPOCH").write_text(f"{t_start:.6f}\n")
    items = {}
    items["N1_quality"] = quality_check(RA, run, 5000)
    items["N2_delta_in"] = delta_in_check(run, HERE / EU_RUN)
    items["N3_wall_output"] = wall_output_check(run)
    items["N4_config"] = ns_config_checks(RA, run)
    items["N5_prepare"] = ns_prepare_info_checks(info, dref["spec"], prelaunch=False)
    pre = HERE / NS_PRELAUNCH
    if pre.is_dir():
        same = {n: (run / n).read_bytes() == (pre / n).read_bytes() for n in ("solverConfig.yaml", RA.PRESTAGE_CONFIG, "bcondConfig.yaml",
                                                                             "nozzle.msh", "probe.yaml", "species_meta.yaml")
                if (run / n).is_file() and (pre / n).is_file()}
        mi = mesh_identity(run / "nozzle.h5", pre / "nozzle.h5")
        fails7 = [f"(i) の投入前の準備と {n} が違う" for n, v in same.items() if not v] + ([] if mi["ok"] else [f"(i) の投入前の準備と nozzle.h5 が違う: {mi['differ']}"])
        items["N7_same_as_prelaunch"] = {"ok": not fails7, "failures": fails7, "files": same, "nozzle_h5": mi, "prelaunch": NS_PRELAUNCH}
    else:
        items["N7_same_as_prelaunch"] = {"ok": False, "failures": [f"(i) の投入前の準備 {NS_PRELAUNCH} が無い (照合できない)"]}
    items["N8_sha"] = {"nozzle_h5_initial": sha256_file(run / "nozzle.h5"), "nozzle_msh": sha256_file(run / "nozzle.msh"),
                       "bcondConfig": sha256_file(run / "bcondConfig.yaml"), "solverConfig_main": sha256_file(run / "solverConfig.yaml"),
                       "solverConfig_prestage": sha256_file(run / RA.PRESTAGE_CONFIG), "problem": sha256_file(HERE / NS_PROBLEM),
                       "euler_ref": dref["spec"], "mesh_hashes": (info.get("mesh") or {}).get("hashes")}
    items["D1_capacity"] = capacity_check(run, pre)
    items["D2_prepared_files_sha256 (記録)"] = files_sha256(run, skip=("H370_PRECHECK.json", "timing.tsv"))
    pc["items"] = items
    pc["t_checked"] = time.time()
    fails = collect_failures(items)
    pc["failures"], pc["ok"] = fails, not fails
    jdump(run / "H370_PRECHECK.json", pc)
    cap = items["D1_capacity"]
    log(f"診断の NS の投入前の検査 {'合格' if pc['ok'] else '不成立'}: 品質 {items['N1_quality'].get('check_mesh_quality_nozzle_h5')}、"
        f"δ_in max {items['N2_delta_in'].get('judged', {}).get('max_abs_mm')} mm、壁の出力 {items['N3_wall_output']['ok']}、"
        f"(i) の準備と同一 {items['N7_same_as_prelaunch']['ok']}、保存の容量 要る {cap.get('need_gib')} GiB / 空き {cap.get('free_gib')} GiB"
        + ("" if pc["ok"] else " — " + "; ".join(fails)))
    if not items["N3_wall_output"]["ok"]:
        (run / "PRECHECK_FAILED").write_text(now() + "\n" + "\n".join(fails) + "\n")
        log("NS が停止の出力で壁の出力を書かない設定 — 止めて親に返す (終了コード 4)")
        timing.stop_evt.set()
        return EXIT_WALL_OUTPUT
    if dry:
        log(f"--dry-run {dry}: 診断の NS の投入前の検査 {'合格' if pc['ok'] else '不成立'} (forge は回さない)")
        return 0 if pc["ok"] else 2
    if not pc["ok"]:
        (run / "PRECHECK_FAILED").write_text(now() + "\n" + "\n".join(fails) + "\n")
        log("前提の検査が不成立 — 診断の NS の forge を回さない")
        timing.stop_evt.set()
        return 2
    # --- 実行: run_staged_ns (S1 → S2 → 本段; 本段は停止規則の見張り)。関門・正式の判定は呼ばない ---
    st4 = items["N4_config"]
    main_cap = RA._cfg_set(RA.ns_stage_configs((run / "solverConfig.yaml").read_text(), "full",
                                                cfg_pre=(run / RA.PRESTAGE_CONFIG).read_text())[1], {"nStepOuter": NS_CAP})
    prepared = {"main_sha256": sha256_text(main_cap), "prestage_sha256": st4["sha256"]["prestage_config"], "want_prev": st4["previous_want"]}
    obs = Observer(RA, run, {k: v for k, v in st4["stage_tags_by_sha"].items() if v != "main"}, "main", prepared,
                   disk_min=max(DISK_START_MIN, int(cap["need_bytes"])))
    obs.install()
    run_end = {"t_start": t_start, "phase": "ns-diag"}
    try:
        rc = RA.run_staged_ns(run, stages="full")
        rc_txt = str(int(rc))
        run_end["stop_rule"] = getattr(rc, "stop_rule", None)
    except Exception:  # noqa: BLE001
        run_end["exception"] = traceback.format_exc()
        log("run_staged_ns が例外:\n" + run_end["exception"])
        rc_txt = "EXC"
    finally:
        obs.uninstall()
    (run / "RUN_RC").write_text(rc_txt + "\n")
    timing.stop_evt.set()
    timing.join(timeout=60)
    run_end.update(t_end=time.time(), run_rc=rc_txt, main_started=(run / RA.MAIN_STAGE_STARTED).is_file())
    jdump(run / "RUN_END.json", run_end)
    try:
        end = diag_end_record(RA, run, rc_txt, run_end, dref, items["N3_wall_output"].get("wall_dump_name"))
    except Exception:  # noqa: BLE001 — 記録の失敗も残す (run は消さない)
        end = {"plan": PLAN_DIAG, "diagnostic_only": True, "update_conditions_met": False, "run_rc": rc_txt, "end_state": "記録の失敗",
               "exception": traceback.format_exc()[-3000:]}
    jdump(run / DIAG_END, end)
    log(f"{run.name} (診断) 終了 RUN_RC={rc_txt} 終了の状態 {end.get('end_state')} 停止の時刻 {end.get('stop_step')} "
        f"本段の出力 {(end.get('main_outputs') or {}).get('n')} (欠け {(end.get('main_outputs') or {}).get('gaps')}) 正式の記録の書き換え "
        f"{'なし' if (end.get('formal_records') or {}).get('unchanged') else 'あり・不明'} ({run_end['t_end'] - t_start:.0f} s) — "
        "感度の判定は python3 h370_judge.py --diag-sens (関門・壁の更新は呼ばない)")
    return 0


# --- h370_problems.py check (作業ツリーの版で; h370_prep.sh が呼ぶ) ------------------------------------------------------------
def check_problems(a) -> int:
    """作業ツリーの h370_problems.py を読み込み (HERE = 作業ツリーの case → design も作業ツリー)、V5d の起動の記録 (moc_v5d_run.log、
    commit されていない) だけを B の case dir のものに差し替えて cmd_check を回す。出力は _band_ab/half_inlet/h370_problems_B.json
    (commit した h370_problems.json を上書きしない)。"""
    load_forge_design(a.commit)
    path = WT / CASE / "h370_problems.py"
    spec = importlib.util.spec_from_file_location("h370_problems_wt", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    log_b = HERE / "_band_ab" / "moc_v5d_run.log"
    if not log_b.is_file():
        die(f"{log_b} が無い (h370_problems.py check が V5d の起動の記録を読む) — 手元の case/45 の _band_ab/moc_v5d_run.log を写す")
    mod.V5D_LOG = log_b
    out = HERE / "_band_ab" / "half_inlet" / "h370_problems_B.json"
    rc = mod.cmd_check(out_json=out)
    log(f"h370_problems.py check (作業ツリー {path}、記録 {log_b.name}) rc={rc} → {out}")
    return rc


def main(argv=None) -> int:
    global LOG_TAG
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", required=True, choices=("euler", "ns", "all", "check-problems", "ns-diag"))
    ap.add_argument("--commit", required=True, help="~/forge-g0cchain の HEAD であるべき commit (h370_prep.sh に渡したもの)")
    ap.add_argument("--dry-run", default=None, metavar="NAME",
                    help="準備と前提の検査だけを run_09xx を使わずに _h370_dry_<NAME>_* で行う (forge・資格の判定・判定は回さない)")
    ap.add_argument("--no-judge", action="store_true", help="NS の後に h370_judge.py を回さない")
    a = ap.parse_args(argv)
    LOG_TAG = f"h370_run {a.phase}"
    if a.dry_run is not None and not re.fullmatch(r"[A-Za-z0-9_.-]+", a.dry_run):
        die("--dry-run の名前は英数字 (_h370_dry_<NAME>_* を作る)")
    log(f"開始 {PLAN} phase {a.phase} commit {a.commit} dry-run {a.dry_run}")
    if a.phase == "check-problems":
        return check_problems(a)
    RA, code = code_checks(a)
    log(f"コード {json.dumps(code, ensure_ascii=False)}")
    if a.phase == "ns-diag":
        return phase_ns_diag(a, RA, code)
    if a.phase in ("euler", "all"):
        rc = phase_euler(a, RA, code)
        if rc != 0 or a.phase == "euler" or a.dry_run:
            return rc
    return phase_ns(a, RA, code)


if __name__ == "__main__":
    sys.exit(main())
