#!/usr/bin/env python3
r"""計算の後の判定を 1 本にまとめる (plan tooling-convergence-and-wall-resolution-gates §4.6・§5.1 #12; ユーザ 2026-10-15
「計算後の収束判定は完全に機械化してほしい」「あほなAIでもね」)。

`run_case.sh` が forge の後に自動で呼ぶ。手で回すときも引数は run だけでよい (**選ぶ場面を作らない**):

    python3 solver_density_cuda/tools/judge_run.py <run_dir>                # 判定して RUN_VERDICT.{json,txt} を書く
    python3 solver_density_cuda/tools/judge_run.py <run_dir> --full-audit   # 非有限の監査を全スナップショット・境界の出力にも
    python3 solver_density_cuda/tools/judge_run.py <run_dir> --check-fresh  # 判定が今の入力に対応するか (書かない)

既存の判定ツール (`check_convergence.py`・`check_quasisteady.py`・`check_wall_resolution.py`・停止規則の記録) の**共通の集約層**であり、
数値の判定基準は各ツールのまま。設計の受理の条件 (deltastar_loop の受け渡しの関門・SERN のゲートなど) は各ゲートに残す。

判定の項目 (`RUN_VERDICT.json` の `checks`、各項目の状態は PASS / FAIL / INDETERMINATE / NOT_APPLICABLE):

| 項目 | 何を見るか | INDETERMINATE (判定不能) になるとき |
| --- | --- | --- |
| `solver` | forge の終了コード (`--solver-rc`、無ければ前の判定の記録 (同じ `forge_run.log`)、無ければ `forge_run.log` の正常終了の行 `Time = … s (wall`)。0 は PASS、停止規則の見張りが止めた (143/137 かつ `STOP_RULE_STATE.json` が REACHED/DIVERGED) は PASS (正常な監視の停止)、他の非 0 は FAIL | 終了コードが分からない |
| `nonfinite` | 残差の全履歴 (`residual_history.csv` と manifest の段の履歴の全行・全 `rms_*` 列) と最終の場 `res_<N>.h5` の `VALUE/*` の非有限。`--full-audit` で全 `res_<N>.h5` と境界の出力も。run の開始より後の `res_nan_*.h5` は FAIL | 履歴・最終の場が無い・読めない・`rms_*` 列が無い |
| `physical` | 最終の場の `ro`>0・`P`>0・`T`>0 (`P`・`T` はあれば) と**静圧 ≤ 全圧** (`total_quantities.total_state`: `VALUE/h0` から全温・全圧; `P > P0·(1 + PT_REL_TOL)` の点が 1 つでもあれば FAIL、P・T が有限で正なのに全圧が非有限の点も FAIL) | 最終の場・`ro` が無い・静圧 ≤ 全圧を確かめられない (`h0` が無い・`thermalMethod` が未対応・物性 (`cp`・`gamma`・種の記録) が読めない) |
| `convergence` | `check_convergence.py` (仕様の mode: auto は manifest があれば `--segment`、無ければ**単段と確かめられたとき**だけ履歴全体)。PASS は終了 0、列の DIVERGED・RISING・STALLED・still converging・ABOVE FLOOR は FAIL | 判定不能の行・参照の拒否・単段と確かめられない・非定常で仕様が無い |
| `stop_rule` | `STOP_RULE_STATE.json` (無ければ `stop_rule.json`) の状態を**そのまま** (`raw`) 記録。REACHED→PASS、CENSORED・DIVERGED・EXEC_ERROR→FAIL、UNAVAILABLE・NOT_USED→NOT_APPLICABLE | DATA_ERROR・MISSING (prepare_info.json が要求したのに状態が無い)・RUNNING・記録が読めない |
| `quasisteady` | 仕様の量だけ (`check_quasisteady.classify_series`)。要求した量と評価した量の完全一致を照合。必須の量の FAIL (NONFINITE・DRIFTING・TRANSIENT-UNSETTLED・受け入れない OSCILLATING) は FAIL | 仕様が無い (対象未設定)・系列・列・量が無い・step の重複・逆行・点数や窓の不足・抽出が適用不能 |
| `wall_resolution` | 実効の粘性と境界条件で適用を決め (`wall_resolution_applies`: 定数粘性 `viscMethod 0` も粘性)、全 no-slip 壁 (kind wall*) の最後の step の出力で `check_wall_resolution.py --step N --groups <全壁> --weight area` | 壁の出力が欠ける・要る量が無い・ツールの行が全壁・step N で揃わない・仕様の applies が実効と違う |

**総合 (真理値表)**: 必須の項目のどれかが FAIL → FAIL (判定不能が併存しても**既知の FAIL を隠さない**); FAIL が無く INDETERMINATE が
ある → INDETERMINATE; 残りに PASS がある → PASS; 全部 NOT_APPLICABLE → NOT_APPLICABLE。必須でない (`required: false`) 量の
FAIL・判定不能は参考として出すが総合には入れない。

**終了コード** (判定モード): 0 = 総合 PASS / 1 = FAIL / 2 = INDETERMINATE (総合 NOT_APPLICABLE も 2: 合格ではない) /
3 = 判定器の例外・引数の誤り (それでも RUN_VERDICT に `state: error`・INDETERMINATE を書こうとする) /
4 = 段階起動の段の結果だけを記録した (`state: deferred`、総合は確定しない)。`--check-fresh` は 0 = 今の入力に対応 /
1 = 入力が変わった (古い判定) / 2 = 判定が無い・読めない。

**確定の時点** (§4.6-2): 単段の run は `run_case.sh` が forge の後に `--from-run-case` で呼んで確定する。段階起動の run
(run に `stage_manifest.json` がある、または親が `FORGE_JUDGE_DEFER=1` を渡した) では、`--from-run-case` の呼び出しは
判定を確定せず、段の結果 (solver・残差の履歴とその時点の最新の場の非有限・物理性・`CONVERGENCE_VERDICT.txt` の副判定) を
`RUN_VERDICT.json` の `stages` に足して `state: deferred` で書く (総合は既知の FAIL があれば FAIL、無ければ INDETERMINATE)。
親 (runner) が段の記録・manifest・停止の記録の確定・出力の整理を終えた後に引数なしで呼んで一度だけ確定する。設計チェーンの
runner は全部の段 (S1・stages=none の本段を含む) で `FORGE_JUDGE_DEFER=1` を渡し、最後に `design/forge_design/evaluate/judge.finalize`
で確定する (段の途中の例外でも、forge を起動していれば確定する)。手で段階起動した run は最後に引数なしで回す (Stop フックが
`state: deferred` のまま終えるのを止める)。`--defer` は段の記録を強制する。

**判定の仕様** (§4.6-1): run の `JUDGE_SPEC.json` (生成器が実行の前に書く) → 無ければ run の `judge.yaml` → case の
`judge.yaml` (手書き) を読む。どれも無ければ準定常は「対象未設定」= INDETERMINATE (黙って省かない)、収束は auto、壁解像は
実効の設定から決める。`judge.yaml` が run の開始より後に書かれていたら `spec.written_after_start: true` を記録する (判定は変えない)。
仕様の未知のキー・型の誤り・重複キーは仕様の誤り (仕様に依る項目はすべて INDETERMINATE)。形式 (版 1):

    {
      "judge_spec_version": 1,                 # JSON では必須 (judge.yaml では省略可 = 1)
      "source": "generator:runner_axismach",   # 任意: 書いたもの
      "convergence": {                         # 任意 (無ければ mode auto)
        "mode": "auto",                        # auto | segment | whole | from_floor | not_applicable
        "drop": 3.0, "tail": 0.2,              # check_convergence の --drop・--tail (既定はツールと同じ)
        "reference": "../run_0101_converged",  # from_floor のとき必須 (run からの相対か絶対)
        "floor_factor": 1.5,                   # from_floor の --floor-factor
        "reason": "..."                        # not_applicable のとき必須 (例: 非定常 LES)
      },
      "quasisteady": {                         # 無ければ「対象未設定」= INDETERMINATE
        "not_applicable": "<理由>",            # 報告する派生量が無い run はこれ (quantities とどちらか一方)
        "quantities": [{
          "name": "C_T",                       # 必須・一意
          "definition": "推力係数 (壁の圧力 + せん断)",  # 必須: 量の定義
          "unit": "-",                         # 必須: 単位
          "region": "ノズル内壁の全体",         # 必須: 対象の領域
          "extractor": {"kind": "csv_column", "file": "force_history.csv", "column": "C_T_with_shear",
                        "step_column": "step"},          # 系列の CSV の列 (run からの相対)
          #  または {"kind": "snapshot", "quantity": "machmax"}  # check_quasisteady の量 (shock・asym・machmax・pmax・
          #          theta・cf_retheta・cf_momentum; 任意で cf_x・cf_ytop・wall_law・mesh) を全 res_<N>.h5 から
          "segment": {"step_min": null, "step_max": null},  # 任意: 判定に使う step の区間 (両端を含む)
          "window": {"tail_frac": 0.4, "min_points": 4, "min_tail_span": 0},   # 末尾窓の割合・区間の最少点数・末尾窓の最少 step 幅
          "criterion": "relative",             # 任意: relative (既定、classify_series の相対の drift・osc) | absolute (下)
          "thresholds": {"drift": 0.05, "osc": 0.10},       # relative: classify の --drift・--osc
          #  absolute (0 を横切る量、例 質量の収支 r): 末尾窓の |平均| ≤ abs_mean かつ 振幅 (最大 − 最小)/2 ≤ abs_amp で STEADY、
          #  振幅の超過は OSCILLATING、|平均| の超過は OFFSET (受け入れられない)。"thresholds": {"abs_mean": 5e-4, "abs_amp": 5e-3}
          "accept": ["STEADY"],                # 合格にする判定 (OSCILLATING を入れたら平均±振幅で報告すること)
          "required": true                     # false なら参考 (総合に入れない)
        }]
      },
      "wall_resolution": {                     # 任意 (無ければ実効の粘性と BC から決める)
        "applies": true,                       # 実効の設定と違えば INDETERMINATE (仕様と実効の食い違い)
        "reason": "...", "target": 1.0, "over_frac": 2.0  # check_wall_resolution の --target・--over-frac
      },
      "notes": ["..."]                         # 任意: 仕様を書いたものの注記 (RUN_VERDICT.txt に写す; 判定は変えない)
    }

適用になる run で no-slip 壁の `outputHDFflg` が 1 でなければ、投入の関門 (`launch_gate.py`) が forge を起動しない (計算の後では壁解像を
判定できず、計算し直しになるため; §4.6-1)。設計チェーンの生成器は機種の既定の仕様を `design/forge_design/evaluate/judge.py` で書く。

手書きの run の `judge.yaml` の例 (平板の Cf を報告する run):

    convergence: {mode: auto}
    quasisteady:
      quantities:
        - name: cf_momentum
          definition: 運動量積分の Cf / Cf_KS (x = 0.6 m)
          unit: "-"
          region: 平板 x ∈ [0.1, 0.95] m
          extractor: {kind: snapshot, quantity: cf_momentum, cf_x: 0.6}
          window: {tail_frac: 0.4, min_points: 6}
    wall_resolution: {applies: true, reason: 低 Re SST の平板}

**出力**: `RUN_VERDICT.json` (各判定の原文・道具の sha256・入力の識別・実行の識別子・確定の状態) と `RUN_VERDICT.txt` (先頭の
1 行に総合)。どちらも一時ファイル → rename で原子的に書く。入力の識別は、テキストの入力 (config・残差の履歴・manifest・
停止の記録・仕様・系列の CSV・`forge_run.log`) は sha256、h5 (最終の場・境界の出力・全出力の一覧) は名前・大きさ・mtime_ns。
`--check-fresh` はこれを今の run と比べる (Stop フックの読み替えは後半)。

**監査の範囲** (§4.6-6、AGENTS.md「NaN / 発散チェック」): 通常は残差の全履歴 + 最終の場。途中の全場は `--full-audit` のときだけ。
"""
import argparse
import csv
import datetime
import glob
import hashlib
import json
import math
import os
import re
import socket
import subprocess
import sys
import tempfile
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

JUDGE_VERSION = 1
SPEC_VERSION = 1
VERDICT_JSON = "RUN_VERDICT.json"
VERDICT_TXT = "RUN_VERDICT.txt"
SPEC_JSON = "JUDGE_SPEC.json"
SPEC_YAML = "judge.yaml"
STATUSES = ("PASS", "FAIL", "INDETERMINATE", "NOT_APPLICABLE")
EXIT_OF = {"PASS": 0, "FAIL": 1, "INDETERMINATE": 2, "NOT_APPLICABLE": 2}
EXIT_ERROR = 3
EXIT_DEFERRED = 4
CHECK_ORDER = ("solver", "nonfinite", "physical", "convergence", "stop_rule", "quasisteady", "wall_resolution")
DEFER_CHECKS = ("solver", "nonfinite", "physical")
TOOLS = ("judge_run.py", "check_convergence.py", "check_quasisteady.py", "check_wall_resolution.py", "stage_manifest.py",
         "total_quantities.py")
# 停止規則の状態 (design/forge_design/evaluate/stop_rule.py の STATES + NOT_USED、ここで足した MISSING・RUNNING) → 項目の状態
STOP_STATUS = {"REACHED": "PASS", "CENSORED": "FAIL", "DIVERGED": "FAIL", "EXEC_ERROR": "FAIL", "DATA_ERROR": "INDETERMINATE",
               "UNAVAILABLE": "NOT_APPLICABLE", "NOT_USED": "NOT_APPLICABLE", "MISSING": "INDETERMINATE", "RUNNING": "INDETERMINATE"}
MONITOR_STOP_RC = (143, 137, -15, -9)          # SIGTERM / SIGKILL (run_case.sh は 128+sig、Python の subprocess は −sig)
QS_SNAPSHOT = ("shock", "asym", "machmax", "pmax", "theta", "cf_retheta", "cf_momentum")
QS_VERDICTS = ("STEADY", "OSCILLATING", "TRANSIENT-UNSETTLED", "DRIFTING")      # accept に書けるもの (absolute の OFFSET は書けない)
QS_CRITERIA = ("relative", "absolute")
# 静圧 ≤ 全圧の許容 (P > P0·(1 + PT_REL_TOL) を違反とする)。float32 の出力 (h0・T・P) の丸めで、よどみ点の近くの P/P0 − 1 は
# 正になる (2026-10-15 実測、最終の場: case/57・case/56 の 25 run で最大 6.4e-7 (CPG 5.2e-7〜6.4e-7、TP 1.6e-7〜2.7e-7)、
# case/45 run_0138_rerun_fullpath_blk4 (イソブタン TP) 1.0e-6; 倍精度の case/57 run_0016_t3a_lm_double は 1.1e-15)。許容はその約 10 倍
PT_REL_TOL = 1.0e-5
WALL_NEED = ("twall_x", "twall_y", "twall_z", "ro", "Ts")                       # check_wall_resolution.py の main と同じ
CONV_FAIL_MARKS = ("<-- DIVERGED", "<-- RISING", "<-- STALLED", "<-- still converging", "<-- ABOVE FLOOR")
_RES_RE = re.compile(r"^res_(\d+)\.h5$")
_BND_RE = re.compile(r"^res_(.+)_(\d+)_(\d+)\.h5$")
HASH_MAX = 512 << 20                           # これより大きいテキストは名前・大きさ・mtime で識別する


class SpecError(Exception):
    """判定の仕様 (JUDGE_SPEC.json・judge.yaml) の誤り。仕様に依る項目は INDETERMINATE にする。"""


class _Parser(argparse.ArgumentParser):
    def error(self, message):                  # 引数の誤りは判定器の誤り (終了 2 = INDETERMINATE と取り違えない)
        self.print_usage(sys.stderr)
        print("%s: error: %s" % (self.prog, message), file=sys.stderr)
        sys.exit(EXIT_ERROR)


# ------------------------------------------------------------------------------------------------ 小道具

def _now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _digest(path):
    """入力の識別: h5 と大きいファイルは名前・大きさ・mtime_ns、他は sha256。無ければ None。"""
    try:
        st = os.stat(path)
    except OSError:
        return None
    if path.endswith(".h5") or st.st_size > HASH_MAX:
        return "stat:%d:%d" % (st.st_size, st.st_mtime_ns)
    return "sha256:" + _sha256_file(path)


def _atomic_write(path, text):
    d = os.path.dirname(os.path.abspath(path))
    fd, tmp = tempfile.mkstemp(dir=d, prefix="." + os.path.basename(path) + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(text)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def _rel(path):
    a = os.path.abspath(path)
    return os.path.relpath(a, ROOT) if a.startswith(ROOT + os.sep) else a


def _yaml_load(text):
    import yaml_strict
    return yaml_strict.load(text)


def _read_yaml(path):
    """(dict か None, 誤りの文字列か None)。重複キー・merge key は誤り (solver と読みが割れる; yaml_strict)。"""
    if not os.path.exists(path):
        return None, "%s が無い" % os.path.basename(path)
    try:
        with open(path, encoding="utf-8") as f:
            doc = _yaml_load(f.read())
    except Exception as e:                                                       # noqa: BLE001
        return None, "%s を読めない (%s: %s)" % (os.path.basename(path), type(e).__name__, str(e)[:200])
    if doc is None:
        doc = {}
    if not isinstance(doc, dict):
        return None, "%s が mapping でない" % os.path.basename(path)
    return doc, None


def _check(status, reason, **kw):
    assert status in STATUSES, status
    out = {"status": status, "reason": reason}
    out.update(kw)
    return out


def res_steps(run):
    """主スナップショット res_<N>.h5 の step の昇順 (境界の出力・res_nan_* は入れない)。"""
    out = []
    for f in os.listdir(run):
        m = _RES_RE.match(f)
        if m:
            out.append(int(m.group(1)))
    return sorted(out)


def boundary_dumps(run):
    """境界の出力 res_<名前>_<physID>_<N>.h5 → [(名前, physID, N, ファイル名)] (res_nan_* は入れない)。"""
    out = []
    for f in os.listdir(run):
        m = _BND_RE.match(f)
        if m and not f.startswith("res_nan_"):
            out.append((m.group(1), int(m.group(2)), int(m.group(3)), f))
    return sorted(out)


def provenance(run):
    """RUN_PROVENANCE.txt (run_case.sh が起動のたびに書く) の主な行 → dict。無ければ {}。"""
    p = os.path.join(run, "RUN_PROVENANCE.txt")
    out = {}
    if not os.path.exists(p):
        return out
    for ln in open(p, encoding="utf-8", errors="replace"):
        m = re.match(r"^(date|forge_bin|forge_sha256|git_head|host)\s*:\s*(.*)$", ln.rstrip("\n"))
        if m and m.group(1) not in out:
            out[m.group(1)] = m.group(2).strip()
    try:
        out["mtime_ns"] = os.stat(p).st_mtime_ns
    except OSError:
        pass
    return out


def run_start_ns(run):
    """run (今の起動) の開始の時刻 [ns]: RUN_PROVENANCE.txt の date、読めなければ None。"""
    d = provenance(run).get("date")
    if not d:
        return None
    try:
        return int(datetime.datetime.fromisoformat(d).timestamp() * 1e9)
    except ValueError:
        return None


# ------------------------------------------------------------------------------------------------ 仕様

def _need(d, key, typ, where, default=None, required=False):
    if key not in d or d[key] is None:
        if required:
            raise SpecError("%s.%s が必須" % (where, key))
        return default
    v = d[key]
    if typ is float:
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(float(v)):
            raise SpecError("%s.%s は有限の数 (%r)" % (where, key, v))
        return float(v)
    if typ is int:
        if isinstance(v, bool) or not isinstance(v, int):
            raise SpecError("%s.%s は整数 (%r)" % (where, key, v))
        return v
    if typ is bool:
        if not isinstance(v, bool):
            raise SpecError("%s.%s は true / false (%r)" % (where, key, v))
        return v
    if typ is str:
        if not isinstance(v, str) or not v.strip():
            raise SpecError("%s.%s は空でない文字列 (%r)" % (where, key, v))
        return v.strip()
    if not isinstance(v, typ):
        raise SpecError("%s.%s の型が違う (%r)" % (where, key, v))
    return v


def _only(d, allowed, where):
    if not isinstance(d, dict):
        raise SpecError("%s は mapping" % where)
    extra = sorted(set(d) - set(allowed))
    if extra:
        raise SpecError("%s に未知のキー %s (使えるのは %s)" % (where, extra, sorted(allowed)))


def normalize_spec(raw, from_json):
    """仕様を検査し、既定を埋めた dict を返す。誤りは SpecError。"""
    _only(raw, ("judge_spec_version", "source", "notes", "convergence", "quasisteady", "wall_resolution"), "spec")
    ver = raw.get("judge_spec_version")
    if ver is None and from_json:
        raise SpecError("JUDGE_SPEC.json に judge_spec_version が無い")
    if ver is not None and ver != SPEC_VERSION:
        raise SpecError("judge_spec_version %r は未対応 (対応は %d)" % (ver, SPEC_VERSION))
    out = {"judge_spec_version": SPEC_VERSION, "source": raw.get("source")}
    notes = raw.get("notes")
    if notes is not None:
        if isinstance(notes, str):
            notes = [notes]
        if not isinstance(notes, list) or any(not isinstance(n, str) for n in notes):
            raise SpecError("notes は文字列か文字列のリスト (%r)" % (notes,))
        out["notes"] = list(notes)
    c = raw.get("convergence") or {}
    _only(c, ("mode", "drop", "tail", "reference", "floor_factor", "reason"), "convergence")
    mode = c.get("mode", "auto")
    if mode not in ("auto", "segment", "whole", "from_floor", "not_applicable"):
        raise SpecError("convergence.mode %r は auto / segment / whole / from_floor / not_applicable" % (mode,))
    conv = {"mode": mode, "drop": _need(c, "drop", float, "convergence", 3.0),
            "tail": _need(c, "tail", float, "convergence", 0.2)}
    if not (0.0 < conv["tail"] < 0.5) or conv["drop"] <= 0.0:
        raise SpecError("convergence.tail は (0, 0.5)、drop は正 (%r, %r)" % (conv["tail"], conv["drop"]))
    if mode == "from_floor":
        conv["reference"] = _need(c, "reference", str, "convergence", required=True)
        conv["floor_factor"] = _need(c, "floor_factor", float, "convergence", 1.5)
    elif "reference" in c or "floor_factor" in c:
        raise SpecError("convergence.reference・floor_factor は mode from_floor のときだけ")
    if mode == "not_applicable":
        conv["reason"] = _need(c, "reason", str, "convergence", required=True)
    elif c.get("reason") is not None:
        conv["reason"] = _need(c, "reason", str, "convergence")
    out["convergence"] = conv
    q = raw.get("quasisteady")
    if q is None:
        out["quasisteady"] = None
    else:
        _only(q, ("not_applicable", "quantities"), "quasisteady")
        if ("not_applicable" in q) == ("quantities" in q):
            raise SpecError("quasisteady は not_applicable (理由) と quantities のどちらか一方")
        if "not_applicable" in q:
            out["quasisteady"] = {"not_applicable": _need(q, "not_applicable", str, "quasisteady", required=True)}
        else:
            qs = q["quantities"]
            if not isinstance(qs, list) or not qs:
                raise SpecError("quasisteady.quantities は空でないリスト")
            items, names = [], set()
            for i, it in enumerate(qs):
                w = "quasisteady.quantities[%d]" % i
                _only(it, ("name", "definition", "unit", "region", "extractor", "segment", "window", "criterion", "thresholds",
                           "accept", "required"), w)
                nm = _need(it, "name", str, w, required=True)
                if nm in names:
                    raise SpecError("%s.name %r が重複" % (w, nm))
                names.add(nm)
                ex = it.get("extractor")
                _only(ex, ("kind", "file", "column", "step_column", "quantity", "cf_x", "cf_ytop", "wall_law", "mesh"),
                      w + ".extractor")
                kind = ex.get("kind")
                if kind == "csv_column":
                    _only(ex, ("kind", "file", "column", "step_column"), w + ".extractor")
                    exn = {"kind": kind, "file": _need(ex, "file", str, w + ".extractor", required=True),
                           "column": _need(ex, "column", str, w + ".extractor", required=True),
                           "step_column": _need(ex, "step_column", str, w + ".extractor", "step")}
                elif kind == "snapshot":
                    _only(ex, ("kind", "quantity", "cf_x", "cf_ytop", "wall_law", "mesh"), w + ".extractor")
                    qn = _need(ex, "quantity", str, w + ".extractor", required=True)
                    if qn not in QS_SNAPSHOT:
                        raise SpecError("%s.extractor.quantity %r は %s のどれか" % (w, qn, list(QS_SNAPSHOT)))
                    law = ex.get("wall_law", "reichardt")
                    if law not in ("reichardt", "spalding"):
                        raise SpecError("%s.extractor.wall_law %r は reichardt / spalding" % (w, law))
                    exn = {"kind": kind, "quantity": qn, "cf_x": _need(ex, "cf_x", float, w + ".extractor", 0.6),
                           "cf_ytop": _need(ex, "cf_ytop", float, w + ".extractor", None), "wall_law": law,
                           "mesh": _need(ex, "mesh", str, w + ".extractor", None)}
                else:
                    raise SpecError("%s.extractor.kind %r は csv_column / snapshot" % (w, kind))
                sg = it.get("segment") or {}
                _only(sg, ("step_min", "step_max"), w + ".segment")
                wn = it.get("window") or {}
                _only(wn, ("tail_frac", "min_points", "min_tail_span"), w + ".window")
                crit = it.get("criterion", "relative")
                if crit not in QS_CRITERIA:
                    raise SpecError("%s.criterion %r は %s のどれか" % (w, crit, list(QS_CRITERIA)))
                th = it.get("thresholds") or {}
                if crit == "relative":
                    _only(th, ("drift", "osc"), w + ".thresholds")
                    thn = {"drift": _need(th, "drift", float, w + ".thresholds", 0.05),
                           "osc": _need(th, "osc", float, w + ".thresholds", 0.10)}
                else:
                    _only(th, ("abs_mean", "abs_amp"), w + ".thresholds")
                    thn = {"abs_mean": _need(th, "abs_mean", float, w + ".thresholds", required=True),
                           "abs_amp": _need(th, "abs_amp", float, w + ".thresholds", required=True)}
                if any(v <= 0.0 for v in thn.values()):
                    raise SpecError("%s.thresholds は正の数 (%r)" % (w, thn))
                acc = it.get("accept", ["STEADY"])
                if not isinstance(acc, list) or not acc or any(a not in QS_VERDICTS for a in acc):
                    raise SpecError("%s.accept は %s の空でないリスト (%r)" % (w, list(QS_VERDICTS), acc))
                win = {"tail_frac": _need(wn, "tail_frac", float, w + ".window", 0.4),
                       "min_points": _need(wn, "min_points", int, w + ".window", 4),
                       "min_tail_span": _need(wn, "min_tail_span", float, w + ".window", 0.0)}
                if not (0.0 < win["tail_frac"] <= 1.0) or win["min_points"] < 2 or win["min_tail_span"] < 0:
                    raise SpecError("%s.window の値の範囲 (tail_frac (0,1]・min_points ≥ 2・min_tail_span ≥ 0)" % w)
                items.append({
                    "name": nm, "definition": _need(it, "definition", str, w, required=True),
                    "unit": _need(it, "unit", str, w, required=True), "region": _need(it, "region", str, w, required=True),
                    "extractor": exn,
                    "segment": {"step_min": _need(sg, "step_min", float, w + ".segment", None),
                                "step_max": _need(sg, "step_max", float, w + ".segment", None)},
                    "window": win, "criterion": crit, "thresholds": thn,
                    "accept": list(acc), "required": _need(it, "required", bool, w, True)})
            out["quasisteady"] = {"quantities": items}
    wr = raw.get("wall_resolution")
    if wr is None:
        out["wall_resolution"] = None
    else:
        _only(wr, ("applies", "reason", "target", "over_frac"), "wall_resolution")
        out["wall_resolution"] = {"applies": _need(wr, "applies", bool, "wall_resolution", required=True),
                                  "reason": _need(wr, "reason", str, "wall_resolution", None),
                                  "target": _need(wr, "target", float, "wall_resolution", 1.0),
                                  "over_frac": _need(wr, "over_frac", float, "wall_resolution", 2.0)}
    return out


def spec_candidates(run):
    """仕様の置き場所 (優先順): run の JUDGE_SPEC.json → run の judge.yaml → case (run の親) の judge.yaml。"""
    run = os.path.abspath(run)
    return [(SPEC_JSON, os.path.join(run, SPEC_JSON)), ("judge.yaml (run)", os.path.join(run, SPEC_YAML)),
            ("judge.yaml (case)", os.path.join(os.path.dirname(run), SPEC_YAML))]


def load_spec(run):
    """{source, path, sha256, written_after_start, spec (normalize 済み) か None, error か None}。"""
    for src, p in spec_candidates(run):
        if not os.path.exists(p):
            continue
        info = {"source": src, "path": _rel(p), "sha256": _sha256_file(p), "spec": None, "error": None}
        st = run_start_ns(run)
        info["written_after_start"] = (st is not None and os.stat(p).st_mtime_ns > st)
        try:
            text = open(p, encoding="utf-8").read()
            if p.endswith(".json"):
                raw = json.loads(text)
            else:
                raw = _yaml_load(text)
            if not isinstance(raw, dict):
                raise SpecError("仕様が mapping でない")
            info["spec"] = normalize_spec(raw, p.endswith(".json"))
        except SpecError as e:
            info["error"] = "仕様の誤り: %s" % e
        except Exception as e:                                                    # noqa: BLE001
            info["error"] = "仕様を読めない (%s: %s)" % (type(e).__name__, str(e)[:300])
        return info
    return {"source": None, "path": None, "sha256": None, "written_after_start": False, "spec": None, "error": None}


# ------------------------------------------------------------------------------------------------ 適用の判定 (仕様側)

def no_slip_walls(bcond):
    """bcondConfig の kind が wall* の境界 (no-slip 壁) → {名前: physID} (slip・出入口・軸は入れない;
    deltastar_loop.no_slip_walls と同じ選び方)。physID が整数に読めなければ ValueError。"""
    return {str(k): int(v["physID"]) for k, v in (bcond or {}).items()
            if isinstance(v, dict) and str(v.get("kind", "")).startswith("wall")}


def wall_resolution_applies(cfg, bcond):
    """壁解像の判定が適用になるか (**実効の粘性と境界条件**で決める、§4.6-1・codex plan 段 judge M4)。
    戻り (applies: True / False / None (決められない), 理由, {壁の名前: physID})。

    粘性: `physProp.viscMethod` 0 = 定数粘性 `visc` (**0 より大きければ粘性**; `gasProperties_d.cu` の viscMethod 0 は
    vis_lam = visc)、1 = Sutherland、2 = 多成分の輸送物性 → 粘性。visc 0 の viscMethod 0 は非粘性 (Euler)。
    境界: kind が wall* の境界が 1 つ以上 (no-slip 壁)。どちらかが読めなければ None。"""
    if not isinstance(bcond, dict):
        return None, "bcondConfig.yaml を読めない", {}
    try:
        walls = no_slip_walls(bcond)
    except (KeyError, TypeError, ValueError) as e:
        return None, "bcondConfig.yaml の壁の physID を読めない (%s)" % e, {}
    if not walls:
        return False, "no-slip 壁 (kind wall*) が無い", {}
    pp = (cfg or {}).get("physProp") if isinstance(cfg, dict) else None
    if not isinstance(pp, dict) or pp.get("viscMethod") is None:
        return None, "physProp.viscMethod を読めない", walls
    try:
        vm = int(pp["viscMethod"])
    except (TypeError, ValueError):
        return None, "physProp.viscMethod %r が整数でない" % (pp["viscMethod"],), walls
    if vm == 0:
        try:
            visc = float(pp.get("visc"))
        except (TypeError, ValueError):
            return None, "viscMethod 0 で physProp.visc を読めない", walls
        if not math.isfinite(visc) or visc < 0:
            return None, "physProp.visc %r が有限の非負でない" % (visc,), walls
        if visc == 0.0:
            return False, "非粘性 (viscMethod 0・visc 0)", walls
        desc = "定数粘性 viscMethod 0・visc %g Pa·s" % visc
    elif vm == 1:
        desc = "Sutherland (viscMethod 1)"
    elif vm == 2:
        desc = "多成分の輸送物性 (viscMethod 2)"
    else:
        return None, "未知の viscMethod %d" % vm, walls
    return True, "粘性 (%s) と no-slip 壁 %s" % (desc, ",".join(sorted(walls))), walls


def is_unsteady(cfg):
    t = (cfg or {}).get("time") if isinstance(cfg, dict) else None
    try:
        return isinstance(t, dict) and int(t.get("unsteady", 0) or 0) == 1
    except (TypeError, ValueError):
        return False


# ------------------------------------------------------------------------------------------------ 各項目

def check_solver(run, solver_rc, prior, stop):
    log = os.path.join(run, "forge_run.log")
    log_sha = _sha256_file(log) if os.path.exists(log) else None
    rc, src = None, None
    if solver_rc is not None:
        rc, src = int(solver_rc), "run_case.sh (--solver-rc)"
    elif prior and log_sha and (prior.get("execution") or {}).get("forge_log_sha256") == log_sha \
            and (prior.get("execution") or {}).get("solver_rc") is not None:
        rc, src = int(prior["execution"]["solver_rc"]), "前の判定の記録 (同じ forge_run.log)"
    elif log_sha:
        with open(log, "rb") as f:
            f.seek(0, 2)
            f.seek(max(0, f.tell() - (256 << 10)))
            tail = f.read().decode("utf-8", "replace")
        if re.search(r"^Time = [0-9.eE+-]+ s \(wall", tail, re.M):
            rc, src = 0, "forge_run.log の正常終了の行 (Time = … s (wall …)"
    ex = {"rc": rc, "rc_source": src, "forge_log_sha256": log_sha}
    if rc is None:
        return _check("INDETERMINATE", "forge の終了コードが分からない (--solver-rc が無く、forge_run.log に正常終了の行が無い)", **ex)
    if rc == 0:
        return _check("PASS", "forge exit 0 (%s)" % src, **ex)
    sst = (stop or {}).get("raw")
    if rc in MONITOR_STOP_RC and sst in ("REACHED", "DIVERGED"):
        return _check("PASS", "forge exit %d = 停止規則の見張りが止めた (正常な監視の停止、停止規則 %s)" % (rc, sst),
                      monitor_stop=True, **ex)
    return _check("FAIL", "forge が異常終了 (終了コード %d、%s)" % (rc, src), **ex)


def _residual_files(run):
    """監査する残差の履歴 (run の residual_history.csv と manifest の段の履歴)。戻り ([パス], [欠けの理由])。"""
    files, missing = [], []
    main = os.path.join(run, "residual_history.csv")
    man = os.path.join(run, "stage_manifest.json")
    stages = []
    if os.path.exists(man):
        try:
            stages = json.load(open(man)).get("stages") or []
        except (OSError, ValueError, AttributeError) as e:
            missing.append("stage_manifest.json を読めない (%s)" % e)
    for st in stages:
        p = os.path.join(run, str(st.get("history", "")))
        if os.path.isfile(p):
            files.append(p)
        elif len(stages) == 1 and os.path.isfile(main):
            pass                                  # 段が 1 つの run は素の履歴がその段 (check_convergence と同じ)
        else:
            missing.append("段 %s の履歴 %s が無い" % (st.get("tag"), st.get("history")))
    if os.path.isfile(main):
        files.append(main)
    elif not stages:
        missing.append("residual_history.csv が無い")
    return files, missing


def audit_residual(path):
    """残差の履歴の全行・全 rms_* 列の非有限 → {rows, cols, nonfinite, first_bad, unparsable}。"""
    out = {"file": os.path.basename(path), "rows": 0, "cols": 0, "nonfinite": 0, "first_bad": None, "unparsable": 0}
    with open(path, newline="") as f:
        rd = csv.reader(f)
        hdr = next(rd, None)
        if not hdr:
            return out
        hdr = [h.strip() for h in hdr]
        idx = [i for i, h in enumerate(hdr) if h.startswith("rms_")]
        out["cols"] = len(idx)
        for row in rd:
            if not row:
                continue
            out["rows"] += 1
            for i in idx:
                try:
                    v = float(row[i])
                except (IndexError, ValueError):
                    out["unparsable"] += 1
                    continue
                if not math.isfinite(v):
                    out["nonfinite"] += 1
                    if out["first_bad"] is None:
                        out["first_bad"] = {"step": row[0] if row else None, "column": hdr[i], "value": row[i]}
    return out


def audit_h5_values(path):
    """h5 の VALUE/* の浮動小数点の量の非有限 → {file, n_fields, bad: {量: 個数}}。"""
    import h5py
    import numpy as np
    out = {"file": os.path.basename(path), "n_fields": 0, "bad": {}}
    with h5py.File(path, "r") as f:
        if "VALUE" not in f:
            out["error"] = "VALUE が無い"
            return out
        for k in f["VALUE"]:
            ds = f["VALUE"][k]
            if not hasattr(ds, "dtype") or ds.dtype.kind not in "fc":
                continue
            out["n_fields"] += 1
            nb = int(np.count_nonzero(~np.isfinite(ds[()])))
            if nb:
                out["bad"][k] = nb
    return out


def check_nonfinite(run, full_audit=False, residual_only_main=False):
    reasons_i, reasons_f, detail = [], [], {"residual": [], "fields": []}
    if residual_only_main:
        p = os.path.join(run, "residual_history.csv")
        files, missing = ([p], []) if os.path.isfile(p) else ([], ["residual_history.csv が無い"])
    else:
        files, missing = _residual_files(run)
    reasons_i += missing
    for p in files:
        try:
            a = audit_residual(p)
        except OSError as e:
            reasons_i.append("%s を読めない (%s)" % (os.path.basename(p), e))
            continue
        detail["residual"].append(a)
        if a["rows"] == 0 or a["cols"] == 0:
            reasons_i.append("%s に行か rms_* 列が無い" % a["file"])
        if a["unparsable"]:
            reasons_i.append("%s に数値に読めない値 %d 個" % (a["file"], a["unparsable"]))
        if a["nonfinite"]:
            reasons_f.append("%s に非有限 %d 個 (最初: step %s・%s)" % (a["file"], a["nonfinite"], a["first_bad"]["step"],
                                                                 a["first_bad"]["column"]))
    steps = res_steps(run)
    if full_audit:
        targets = ["res_%d.h5" % s for s in steps] + [b[3] for b in boundary_dumps(run)]
    else:
        targets = ["res_%d.h5" % steps[-1]] if steps else []
    if not steps:
        reasons_i.append("最終の場 res_<N>.h5 が無い")
    for fn in targets:
        try:
            a = audit_h5_values(os.path.join(run, fn))
        except Exception as e:                                                    # noqa: BLE001
            reasons_i.append("%s を読めない (%s)" % (fn, type(e).__name__))
            continue
        detail["fields"].append(a)
        if a.get("error"):
            reasons_i.append("%s: %s" % (fn, a["error"]))
        if a["bad"]:
            reasons_f.append("%s に非有限 %s" % (fn, a["bad"]))
    st = run_start_ns(run)
    nan_dumps = sorted(f for f in os.listdir(run) if f.startswith("res_nan_") and f.endswith(".h5")
                       and (st is None or os.stat(os.path.join(run, f)).st_mtime_ns >= st))
    if nan_dumps:
        reasons_f.append("forge が NaN を検出したダンプ %s がある" % nan_dumps[:3])
    detail["res_nan"] = nan_dumps
    scope = "残差の履歴 %d 本 + %s" % (len(files), ("全出力 %d 個" % len(targets)) if full_audit else
                                    ("最終の場 %s" % (targets[0] if targets else "なし")))
    if reasons_f:
        return _check("FAIL", "; ".join(reasons_f + reasons_i), scope=scope, full_audit=full_audit, **detail)
    if reasons_i:
        return _check("INDETERMINATE", "; ".join(reasons_i), scope=scope, full_audit=full_audit, **detail)
    nrow = sum(a["rows"] for a in detail["residual"])
    return _check("PASS", "%s に非有限なし (残差 %d 行)" % (scope, nrow), scope=scope, full_audit=full_audit, **detail)


def check_physical(run):
    import h5py
    import numpy as np
    steps = res_steps(run)
    if not steps:
        return _check("INDETERMINATE", "最終の場 res_<N>.h5 が無い")
    fn = "res_%d.h5" % steps[-1]
    out = {"file": fn, "min": {}, "n_nonpositive": {}}
    try:
        with h5py.File(os.path.join(run, fn), "r") as f:
            V = f.get("VALUE")
            if V is None or "ro" not in V:
                return _check("INDETERMINATE", "%s に VALUE/ro が無い" % fn, **out)
            for k in ("ro", "P", "T"):
                if k in V:
                    a = np.asarray(V[k][()], float)
                    fin = a[np.isfinite(a)]
                    out["min"][k] = float(fin.min()) if fin.size else None
                    out["n_nonpositive"][k] = int(np.count_nonzero(fin <= 0.0))
    except Exception as e:                                                        # noqa: BLE001
        return _check("INDETERMINATE", "%s を読めない (%s)" % (fn, type(e).__name__), **out)
    bad = {k: n for k, n in out["n_nonpositive"].items() if n}
    names = "・".join("%s>0" % k for k in out["min"])
    pt_status, pt_why, out["static_le_total"] = check_static_le_total(run, fn)
    fails = []
    if bad:
        fails.append("%s で 0 以下: %s (最小 %s)" % (fn, bad, {k: out["min"][k] for k in bad}))
    if pt_status == "FAIL":
        fails.append(pt_why)
    if fails:
        return _check("FAIL", "; ".join(fails + ([pt_why] if pt_status == "INDETERMINATE" else [])), **out)
    if pt_status == "INDETERMINATE":
        return _check("INDETERMINATE", "%s (%s) は PASS、%s" % (names, fn, pt_why), **out)
    return _check("PASS", "%s・%s (%s)" % (names, pt_why, fn), **out)


def check_static_le_total(run, fn):
    """静圧 ≤ 全圧 (plan §5.1 #12a (6)): `total_quantities.total_state` (ソルバの `VALUE/h0` から全温・全圧; スクリプトで
    T + u²/2c_p を組まない) で最終の場の P0 を作り、P > P0·(1 + PT_REL_TOL) の点を数える。P・T が有限で正なのに P0 が非有限の点も
    違反 (h0 が負など)。戻り (PASS / FAIL / INDETERMINATE, 理由, 詳細)。h0 が無い・thermalMethod が未対応・物性が読めないときは
    INDETERMINATE (total_state は SystemExit で止まるので BaseException で受ける)。"""
    import numpy as np
    import total_quantities as TQ
    try:
        st = TQ.total_state(run, os.path.join(run, fn))
    except KeyboardInterrupt:
        raise
    except BaseException as e:                                                    # noqa: BLE001
        return ("INDETERMINATE", "静圧 ≤ 全圧を確かめられない (total_quantities: %s: %s)" % (type(e).__name__, str(e)[:240]),
                {"error": "%s: %s" % (type(e).__name__, str(e)[:400])})
    P = np.asarray(st["P"], float)
    T = np.asarray(st["T"], float)
    P0 = np.asarray(st["P0"], float)
    ok_in = np.isfinite(P) & np.isfinite(T) & (P > 0.0) & (T > 0.0)
    nonfin = int(np.count_nonzero(ok_in & ~np.isfinite(P0)))
    sel = ok_in & np.isfinite(P0)
    rel = P[sel] / P0[sel] - 1.0
    nbad = int(np.count_nonzero(rel > PT_REL_TOL))
    det = {"method": st.get("method"), "h0_includes_k": st.get("includes_k"), "tol": PT_REL_TOL,
           "n_points": int(np.count_nonzero(sel)), "n_violations": nbad, "n_nonfinite_P0": nonfin,
           "max_P_over_P0_minus_1": (float(rel.max()) if rel.size else None), "warn": st.get("warn") or []}
    if st.get("species_source"):
        det["species_source"] = st["species_source"]
    if not rel.size:
        return "INDETERMINATE", "静圧 ≤ 全圧: P・T が有限で正の点が無い", det
    if nbad or nonfin:
        return ("FAIL", "静圧 > 全圧: %d 点 (P/P0 − 1 の最大 %.3e > 許容 %.0e)%s (%s)"
                % (nbad, rel.max(), PT_REL_TOL, ("・全圧が非有限 %d 点" % nonfin) if nonfin else "", st.get("method")), det)
    return "PASS", "P ≤ P0 (P/P0 − 1 の最大 %.2e ≤ 許容 %.0e; %s)" % (rel.max(), PT_REL_TOL, st.get("method")), det


def confirm_single_stage(run):
    """manifest の無い run が単段と確かめられるか → (True/False, 理由)。段の履歴 (residual_history_<段>.csv) が無く、
    残差の履歴の step が単調に増える (継ぎ足しが無い) ときだけ True。"""
    p = os.path.join(run, "residual_history.csv")
    if not os.path.isfile(p):
        return False, "residual_history.csv が無い"
    extra = sorted(os.path.basename(f) for f in glob.glob(os.path.join(run, "residual_history_*.csv"))
                   if os.path.basename(f) != "residual_history_segment.csv")
    if extra:
        return False, "段の履歴 %s があるのに stage_manifest.json が無い (段階起動の区間を決められない)" % extra[:3]
    import check_convergence as cc
    try:
        rows, _ = cc.load_series(p)
        st = [int(float(r["step"])) for r in rows]
    except (KeyError, ValueError, OSError) as e:
        return False, "残差の履歴の step を読めない (%s)" % e
    if not st:
        return False, "残差の履歴に行が無い"
    back = sum(1 for a, b in zip(st, st[1:]) if b <= a)
    if back:
        return False, "残差の履歴の step が単調に増えない (%d か所; 継ぎ足しの履歴か)" % back
    return True, "stage_manifest.json 無し・段の履歴無し・step が単調に増える → 単段とみなした"


def _parse_convergence(rc, out):
    lines = out.splitlines()
    verdict = [ln.strip() for ln in lines if ln.startswith("=== ") and "->" in ln]
    seg = [ln.strip() for ln in lines if "[segment] 判定区間" in ln]
    fails = sorted({m for m in CONV_FAIL_MARKS for ln in lines if m in ln})
    undet = [ln.strip() for ln in lines if "判定不能" in ln or "COLUMN MISMATCH" in ln or "REFUSED" in ln
             or "NO residual_history" in ln or "empty residual" in ln]
    info = {"verdict_lines": verdict, "segment_lines": seg, "fail_marks": fails, "undetermined_lines": undet[:10]}
    short = None
    if verdict:                                  # "=== <run>  [last step N]  -> VERDICT ===" → "VERDICT [last step N]"
        m = re.match(r"^=== .*?(\[last step [^\]]*\])\s*->\s*(.*?)\s*===$", verdict[-1])
        short = "%s %s" % (m.group(2), m.group(1)) if m else verdict[-1]
    if rc == 0 and "OVERALL: ALL PASS" in out:
        return "PASS", (short or "ALL PASS"), info
    if fails:
        return "FAIL", "%s (%s)" % (short or "NOT CONVERGED", ", ".join(fails)), info
    return "INDETERMINATE", (undet[0] if undet else "check_convergence の終了コード %d (判定の行が無い)" % rc), info


def check_convergence(run, spec, cfg):
    conv = (spec or {}).get("convergence") or {"mode": "auto", "drop": 3.0, "tail": 0.2}
    mode = conv["mode"]
    has_man = os.path.exists(os.path.join(run, "stage_manifest.json"))
    base = {"mode": mode}
    if mode == "not_applicable":
        return _check("NOT_APPLICABLE", "仕様: %s" % conv["reason"], **base)
    use_segment = False
    if mode == "segment" or (mode in ("auto", "from_floor") and has_man):
        use_segment = True
    elif mode in ("auto", "from_floor"):
        if mode == "auto" and is_unsteady(cfg):
            return _check("INDETERMINATE", "非定常の run (time.unsteady 1): 収束の扱いが仕様に無い "
                                           "(JUDGE_SPEC.json・judge.yaml の convergence.mode で決める)", **base)
        ok, why = confirm_single_stage(run)
        base["single_stage"] = why
        if not ok:
            return _check("INDETERMINATE", "単段と確かめられない: %s" % why, **base)
    tool = os.path.join(HERE, "check_convergence.py")
    cmd = [sys.executable, tool, run] + (["--segment"] if use_segment else []) + \
          ["--drop", "%g" % conv["drop"], "--tail", "%g" % conv["tail"]]
    if mode == "from_floor":
        ref = conv["reference"] if os.path.isabs(conv["reference"]) else os.path.normpath(os.path.join(run, conv["reference"]))
        cmd += ["--from-floor", ref, "--floor-factor", "%g" % conv["floor_factor"]]
        base["reference"] = _rel(ref)
    base["interval"] = "segment (stage_manifest.json の最後の区間)" if use_segment else "whole (residual_history.csv 全体)"
    r = subprocess.run(cmd, capture_output=True, text=True)
    out = r.stdout + r.stderr
    status, why, info = _parse_convergence(r.returncode, out)
    return _check(status, "[%s] %s" % ("segment" if use_segment else "whole", why), cmd=[_rel(c) if os.path.exists(c) else c
                                                                                          for c in cmd],
                  rc=r.returncode, output=out, **base, **info)


def read_stop_state(run):
    """停止規則の状態を**そのまま** → {raw, source, final, record}。"""
    p = os.path.join(run, "STOP_RULE_STATE.json")
    if os.path.exists(p):
        try:
            st = json.load(open(p))
            return {"raw": st.get("status"), "source": "STOP_RULE_STATE.json", "final": True, "record": st}
        except (OSError, ValueError) as e:
            return {"raw": None, "source": "STOP_RULE_STATE.json", "final": True, "error": "読めない (%s)" % e}
    p = os.path.join(run, "stop_rule.json")
    if os.path.exists(p):
        try:
            st = json.load(open(p))
            return {"raw": st.get("status"), "source": "stop_rule.json (本段の終了の状態はまだ無い)", "final": False,
                    "record": {k: st.get(k) for k in ("status", "stop_step", "reason", "started", "finished")}}
        except (OSError, ValueError) as e:
            return {"raw": None, "source": "stop_rule.json", "final": False, "error": "読めない (%s)" % e}
    try:
        sr = json.load(open(os.path.join(run, "prepare_info.json"))).get("stop_rule")
    except (OSError, ValueError, AttributeError):
        sr = None
    if isinstance(sr, dict) and sr.get("requested"):
        return {"raw": "MISSING", "source": "prepare_info.json は停止規則を要求しているのに STOP_RULE_STATE.json が無い",
                "final": False}
    return {"raw": "NOT_USED", "source": "停止規則の記録が無い (停止規則を使わない run)", "final": True}


def check_stop_rule(stop):
    raw = stop.get("raw")
    if stop.get("error"):
        return _check("INDETERMINATE", "%s %s" % (stop["source"], stop["error"]), **stop)
    if raw not in STOP_STATUS:
        return _check("INDETERMINATE", "停止規則の状態 %r は既知でない (%s)" % (raw, stop["source"]), **stop)
    st = STOP_STATUS[raw]
    if raw in ("REACHED", "CENSORED", "DIVERGED", "EXEC_ERROR", "DATA_ERROR", "UNAVAILABLE") and not stop.get("final"):
        st = "INDETERMINATE"                     # 見張りの途中の記録 (STOP_RULE_STATE.json がまだ無い) では確定しない
    return _check(st, "%s (%s)" % (raw, stop["source"]), **stop)


def _series_csv(run, ex):
    """csv_column の系列 → (steps, values, 誤り)。"""
    p = ex["file"] if os.path.isabs(ex["file"]) else os.path.join(run, ex["file"])
    if not os.path.isfile(p):
        return None, None, "系列のファイル %s が無い" % ex["file"]
    with open(p, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return None, None, "%s に行が無い" % ex["file"]
    miss = [c for c in (ex["step_column"], ex["column"]) if c not in rows[0]]
    if miss:
        return None, None, "%s に列 %s が無い" % (ex["file"], miss)
    steps, vals = [], []
    for r in rows:
        try:
            steps.append(float(r[ex["step_column"]]))
        except (TypeError, ValueError):
            return None, None, "%s の %s に数値でない値 %r" % (ex["file"], ex["step_column"], r[ex["step_column"]])
        v = r[ex["column"]]
        try:
            vals.append(float(v) if v not in ("", None, "None") else float("nan"))
        except ValueError:
            vals.append(float("nan"))
    return steps, vals, None


def classify_absolute(steps, vals, tail_frac, abs_mean, abs_amp):
    """0 を横切る量 (質量の収支 r など; 相対の drift・osc は平均が 0 に近いと意味を持たない) の判定: 末尾窓 (classify と同じ
    k = max(3, ⌈tail_frac·n⌉) 点) の平均 m と振幅 a = (最大 − 最小)/2 で、|m| > abs_mean は OFFSET (受け入れられない)、
    a > abs_amp は OSCILLATING (平均±振幅で報告)、どちらも満たせば STEADY。戻り (verdict, detail, (平均, 振幅))。
    呼ぶ前に非有限・step の単調性・点数は検査済み (check_quasisteady の内側)。"""
    k = max(3, int(math.ceil(tail_frac * len(vals))))
    st, vt = list(steps)[-k:], [float(x) for x in vals][-k:]
    m = sum(vt) / len(vt)
    a = (max(vt) - min(vt)) / 2.0
    detail = ("tail mean=%.4g (|平均| ≤ %g)  amp=%.4g (≤ %g)  [末尾 %d 点 step %g–%g]"
              % (m, abs_mean, a, abs_amp, len(vt), st[0], st[-1]))
    if abs(m) > abs_mean:
        return "OFFSET", detail + "  -> |平均| が許容を超える", (m, a)
    if a > abs_amp:
        return "OSCILLATING", detail + "  -> report %.4g +/- %.2g" % (m, a), (m, a)
    return "STEADY", detail, (m, a)


def check_quasisteady(run, spec_info):
    spec = spec_info.get("spec")
    if spec_info.get("error"):
        return _check("INDETERMINATE", spec_info["error"])
    q = (spec or {}).get("quasisteady")
    if q is None:
        return _check("INDETERMINATE", "対象未設定 (%s)" % ("仕様に quasisteady が無い" if spec else
                                                       "JUDGE_SPEC.json も judge.yaml も無い"))
    if "not_applicable" in q:
        return _check("NOT_APPLICABLE", "仕様: %s" % q["not_applicable"])
    import numpy as np
    import check_quasisteady as cq
    items = q["quantities"]
    results = {}
    series = {}
    # 抽出: csv_column は 1 本ずつ、snapshot は抽出の設定ごとにまとめて res_<N>.h5 を 1 回ずつ読む
    snap_groups = {}
    for it in items:
        ex = it["extractor"]
        if ex["kind"] == "csv_column":
            s, v, err = _series_csv(run, ex)
            series[it["name"]] = (s, v, err)
        else:
            key = json.dumps({k: ex[k] for k in ("cf_x", "cf_ytop", "wall_law", "mesh")}, sort_keys=True)
            snap_groups.setdefault(key, []).append(it)
    for key, its in snap_groups.items():
        opt = json.loads(key)
        mesh = cq.find_mesh(run, opt["mesh"])
        fs, steps = cq.res_files(run)
        if mesh is None or not fs:
            for it in its:
                series[it["name"]] = (None, None, "メッシュ (/CELLS/centCoords) か res_<N>.h5 が無い")
            continue
        cc_ = cq.centroids(mesh)
        want = [it["extractor"]["quantity"] for it in its]
        exs = cq.build_extractors(cc_, want, opt["cf_x"], opt["cf_ytop"], opt["wall_law"])
        ok_q = [w for w in want if w in exs]
        ser = cq.extract_series(fs, cc_, exs, ok_q)
        for it in its:
            w = it["extractor"]["quantity"]
            if w not in exs:
                series[it["name"]] = (None, None, "抽出が適用不能 (%s: 例 上下対称でないメッシュの asym)" % w)
            else:
                series[it["name"]] = ([float(s_) for s_ in steps], ser[w], None)
    for it in items:
        nm = it["name"]
        s, v, err = series.get(nm, (None, None, "評価していない"))
        r = {"required": it["required"], "definition": it["definition"], "unit": it["unit"], "region": it["region"],
             "extractor": it["extractor"], "accept": it["accept"]}
        if err:
            results[nm] = _check("INDETERMINATE", err, **r)
            continue
        lo, hi = it["segment"]["step_min"], it["segment"]["step_max"]
        keep = [i for i, x in enumerate(s) if (lo is None or x >= lo) and (hi is None or x <= hi)]
        s = [s[i] for i in keep]
        v = [v[i] for i in keep]
        r.update(n=len(v), steps=[s[0], s[-1]] if s else None)
        if not v:
            results[nm] = _check("INDETERMINATE", "区間 [%s, %s] に点が無い (空集合)" % (lo, hi), **r)
            continue
        if any(not math.isfinite(x) for x in v):
            nb = sum(1 for x in v if not math.isfinite(x))
            results[nm] = _check("FAIL", "NONFINITE (%d/%d 点が非有限)" % (nb, len(v)), verdict="NONFINITE", **r)
            continue
        ds = np.diff(np.asarray(s, float))
        if ds.size and (np.any(ds <= 0)):
            results[nm] = _check("INDETERMINATE", "step が単調に増えない (重複 %d・逆行 %d)" %
                                 (int(np.count_nonzero(ds == 0)), int(np.count_nonzero(ds < 0))), **r)
            continue
        win = it["window"]
        if len(v) < win["min_points"]:
            results[nm] = _check("INDETERMINATE", "窓の不足: 点数 %d < min_points %d" % (len(v), win["min_points"]), **r)
            continue
        k = max(3, int(math.ceil(win["tail_frac"] * len(v))))
        span = s[-1] - s[-min(k, len(s))]
        r["tail_points"], r["tail_span"] = min(k, len(s)), span
        if span < win["min_tail_span"]:
            results[nm] = _check("INDETERMINATE", "窓の不足: 末尾窓の step 幅 %g < min_tail_span %g" % (span, win["min_tail_span"]), **r)
            continue
        r["criterion"] = it["criterion"]
        if it["criterion"] == "absolute":
            verdict, detail, ma = classify_absolute(s, v, win["tail_frac"], it["thresholds"]["abs_mean"],
                                                    it["thresholds"]["abs_amp"])
        else:
            verdict, detail, ma = cq.classify_series(s, v, win["tail_frac"], it["thresholds"]["drift"], it["thresholds"]["osc"],
                                                     win["min_points"])
        r.update(verdict=verdict, detail=detail)
        if ma is not None:
            r.update(mean=float(ma[0]), amp=float(ma[1]))
        if verdict in it["accept"]:
            results[nm] = _check("PASS", "%s  %s" % (verdict, detail), **r)
        elif verdict == "INDETERMINATE":
            results[nm] = _check("INDETERMINATE", detail, **r)
        else:
            results[nm] = _check("FAIL", "%s  %s" % (verdict, detail), **r)
    want_names = [it["name"] for it in items]
    if sorted(results) != sorted(want_names):            # 要求と評価の完全一致 (防御; ここに来ないはず)
        return _check("INDETERMINATE", "要求した量 %s と評価した量 %s が一致しない" % (want_names, sorted(results)),
                      quantities=results)
    req = {n: r for n, r in results.items() if r["required"]}
    adv = {n: r for n, r in results.items() if not r["required"]}
    st = _combine([r["status"] for r in req.values()]) if req else "NOT_APPLICABLE"
    bad = ["%s %s (%s)" % (n, r["status"], r["reason"][:200]) for n, r in req.items() if r["status"] != "PASS"]
    why = ("必須の量 %d 個すべて受け入れ" % len(req)) if st == "PASS" else ("; ".join(bad) if bad else "必須の量が無い")
    if adv:
        why += " (参考の量: %s)" % ", ".join("%s %s" % (n, r["status"]) for n, r in adv.items())
    return _check(st, why, quantities=results)


def check_wall_resolution(run, cfg, cfg_err, bcond, bcond_err, spec_info):
    spec = spec_info.get("spec")
    if spec_info.get("error"):
        return _check("INDETERMINATE", spec_info["error"])
    if cfg_err:
        return _check("INDETERMINATE", "solverConfig.yaml: %s" % cfg_err)
    if bcond_err:
        return _check("INDETERMINATE", "bcondConfig.yaml: %s" % bcond_err)
    applies, why, walls = wall_resolution_applies(cfg, bcond)
    wr = (spec or {}).get("wall_resolution")
    base = {"applies": applies, "applies_reason": why, "walls": walls, "spec": wr}
    if applies is None:
        return _check("INDETERMINATE", "適用を決められない: %s" % why, **base)
    if wr is not None and wr["applies"] != applies:
        return _check("INDETERMINATE", "仕様の applies %s が実効の設定 (%s → %s) と違う" % (wr["applies"], why, applies), **base)
    if not applies:
        return _check("NOT_APPLICABLE", why, **base)
    steps = res_steps(run)
    if not steps:
        return _check("INDETERMINATE", "最後の step を決められない (res_<N>.h5 が無い)", **base)
    n = steps[-1]
    base["step"] = n
    miss = []
    import h5py
    for name, pid in sorted(walls.items()):
        fn = "res_%s_%d_%d.h5" % (name, pid, n)
        p = os.path.join(run, fn)
        if not os.path.isfile(p):
            miss.append("%s (%s が無い; 壁の出力 outputHDFflg を確かめる)" % (name, fn))
            continue
        try:
            with h5py.File(p, "r") as f:
                lack = [k for k in WALL_NEED if "VALUE/" + k not in f] + ([] if "MESH/COORD" in f else ["MESH/COORD"])
        except Exception as e:                                                    # noqa: BLE001
            miss.append("%s (%s を読めない: %s)" % (name, fn, type(e).__name__))
            continue
        if lack:
            miss.append("%s (%s に %s が無い)" % (name, fn, ",".join(lack)))
    if miss:
        return _check("INDETERMINATE", "最後の step %d の壁の出力が揃っていない: %s" % (n, "; ".join(miss)), **base)
    target = wr["target"] if wr else 1.0
    over = wr["over_frac"] if wr else 2.0
    tool = os.path.join(HERE, "check_wall_resolution.py")
    cmd = [sys.executable, tool, run, "--step", str(n), "--groups", ",".join(sorted(walls)), "--weight", "area",
           "--target", "%g" % target, "--over-frac", "%g" % over]
    r = subprocess.run(cmd, capture_output=True, text=True)
    out = r.stdout + r.stderr
    rows = {m.group(1): int(m.group(2)) for m in re.finditer(r"^  (\S+)\s+step\s+(\d+)\s+y1 = ", out, re.M)}
    vline = [ln.strip() for ln in out.splitlines() if ln.startswith("VERDICT")]
    base.update(cmd=[_rel(c) if os.path.exists(c) else c for c in cmd], rc=r.returncode, output=out, rows=rows,
                verdict_line=vline[-1] if vline else None)
    status = {0: "PASS", 1: "FAIL"}.get(r.returncode, "INDETERMINATE")
    if status != "INDETERMINATE":
        off = sorted(nm for nm in walls if rows.get(nm) != n) + sorted(nm for nm in rows if nm not in walls)
        if off:
            return _check("INDETERMINATE", "ツールの壁ごとの行が全部の壁・step %d で揃っていない (%s)" % (n, ", ".join(off)), **base)
    m = re.search(r"を超える面積割合 最大 ([0-9.]+) % \(許容 ([0-9.eE+-]+) %\)", out)
    if m:
        base.update(over_area_pct=float(m.group(1)), over_area_allow_pct=float(m.group(2)))
    return _check(status, "%s (全 no-slip 壁 %d 枚・step %d・面積の割合%s)" % (
        vline[-1] if vline else "終了コード %d" % r.returncode, len(walls), n,
        (" y₁⁺ > %g の面積 最大 %s %% / 許容 %s %%" % (target, m.group(1), m.group(2))) if m else ""), **base)


# ------------------------------------------------------------------------------------------------ 総合・記録

def _combine(statuses):
    """真理値表: FAIL > INDETERMINATE > PASS > NOT_APPLICABLE (既知の FAIL を判定不能で隠さない)。"""
    s = list(statuses)
    if "FAIL" in s:
        return "FAIL"
    if "INDETERMINATE" in s:
        return "INDETERMINATE"
    if "PASS" in s:
        return "PASS"
    return "NOT_APPLICABLE"


def overall_of(checks):
    st = _combine(c["status"] for c in checks.values())
    fails = [k for k, c in checks.items() if c["status"] == "FAIL"]
    undet = [k for k, c in checks.items() if c["status"] == "INDETERMINATE"]
    if st == "FAIL":
        why = "FAIL: %s" % ", ".join(fails) + ((" (判定不能も: %s)" % ", ".join(undet)) if undet else "")
    elif st == "INDETERMINATE":
        why = "判定不能: %s" % ", ".join(undet)
    elif st == "PASS":
        why = "全項目 PASS か適用外"
    else:
        why = "全項目が適用外"
    return st, why


def tool_hashes():
    return {t: _sha256_file(os.path.join(HERE, t)) for t in TOOLS if os.path.exists(os.path.join(HERE, t))}


def collect_inputs(run, spec_info=None):
    """判定が依る入力の識別 {相対名: digest か None (無い)}。鮮度の照合 (--check-fresh) はこれを比べる。"""
    run = os.path.abspath(run)
    if spec_info is None:
        spec_info = load_spec(run)
    names = ["solverConfig.yaml", "bcondConfig.yaml", "residual_history.csv", "stage_manifest.json", "forge_launches.jsonl",
             "STOP_RULE_STATE.json", "stop_rule.json", "prepare_info.json", "forge_run.log", SPEC_JSON, SPEC_YAML]
    out = {nm: _digest(os.path.join(run, nm)) for nm in names}
    out["../" + SPEC_YAML] = _digest(os.path.join(os.path.dirname(run), SPEC_YAML))
    try:
        for st in json.load(open(os.path.join(run, "stage_manifest.json"))).get("stages") or []:
            h = str(st.get("history", ""))
            out[h] = _digest(os.path.join(run, h))
    except (OSError, ValueError, AttributeError):
        pass
    spec = spec_info.get("spec") or {}
    for it in ((spec.get("quasisteady") or {}).get("quantities") or []):
        ex = it["extractor"]
        if ex["kind"] == "csv_column":
            out[ex["file"]] = _digest(ex["file"] if os.path.isabs(ex["file"]) else os.path.join(run, ex["file"]))
    conv = spec.get("convergence") or {}
    if conv.get("mode") == "from_floor":
        ref = conv["reference"] if os.path.isabs(conv["reference"]) else os.path.normpath(os.path.join(run, conv["reference"]))
        out["reference:" + conv["reference"]] = _digest(os.path.join(ref, "residual_history.csv"))
    outs = sorted(f for f in os.listdir(run) if f.startswith("res_") and f.endswith(".h5"))
    h = hashlib.sha256()
    for f in outs:
        h.update(("%s:%s\n" % (f, _digest(os.path.join(run, f)))).encode())
    out["res_*.h5 (名前・大きさ・mtime の一覧)"] = "sha256:%s (%d 個)" % (h.hexdigest(), len(outs))
    return out


def execution_info(run, solver):
    pv = provenance(run)
    ex = {"provenance_date": pv.get("date"), "forge_sha256": pv.get("forge_sha256"), "git_head": pv.get("git_head"),
          "host": pv.get("host"), "forge_log_sha256": solver.get("forge_log_sha256"), "solver_rc": solver.get("rc"),
          "solver_rc_source": solver.get("rc_source")}
    ex["run_id"] = hashlib.sha256(("%s|%s" % (ex["provenance_date"], ex["forge_log_sha256"])).encode()).hexdigest()[:16]
    return ex


def _short(c):
    return c["reason"] if len(c["reason"]) < 400 else c["reason"][:400] + " …"


def render_txt(v):
    lines = []
    head = "OVERALL: %s" % v["overall"]
    if v["state"] != "final":
        head += "  (state %s: %s)" % (v["state"], v.get("state_note", ""))
    lines.append(head)
    lines.append("  run      : %s  (判定 %s、run_id %s)" % (v["run_rel"], v["judged_at"], (v.get("execution") or {}).get("run_id")))
    lines.append("  総合     : %s" % v["overall_reason"])
    for k in CHECK_ORDER:
        c = v["checks"].get(k)
        if c is None:
            continue
        lines.append("  %-16s: %-15s %s" % (k, c["status"], _short(c)))
        if k == "quasisteady":
            for nm, r in (c.get("quantities") or {}).items():
                lines.append("      %-24s %-15s %s%s" % (nm, r["status"], r["reason"][:300], "" if r["required"] else "  [参考]"))
        if k == "convergence" and c.get("segment_lines"):
            lines.append("      %s" % c["segment_lines"][-1])
    for s in v.get("stages") or []:
        lines.append("  段 %-3s %s: %s" % (s["seq"], s.get("recorded_at"), "; ".join(
            "%s %s" % (k, s["checks"][k]["status"]) for k in DEFER_CHECKS if k in s["checks"])))
    sp = v.get("spec") or {}
    lines.append("  仕様     : %s%s%s" % (sp.get("source") or "なし (JUDGE_SPEC.json も judge.yaml も無い)",
                                       ("  sha256 %s" % sp["sha256"][:16]) if sp.get("sha256") else "",
                                       "  (run の開始より後に書かれた)" if sp.get("written_after_start") else ""))
    if sp.get("error"):
        lines.append("  仕様の誤り: %s" % sp["error"])
    for n in ((sp.get("spec") or {}).get("notes") or []):
        lines.append("  仕様の注記: %s" % n)
    lines.append("  詳細     : %s (各判定の原文・道具と入力の識別)" % VERDICT_JSON)
    return "\n".join(lines) + "\n"


def write_verdict(run, v):
    _atomic_write(os.path.join(run, VERDICT_JSON), json.dumps(v, indent=1, ensure_ascii=False, default=str) + "\n")
    txt = render_txt(v)
    _atomic_write(os.path.join(run, VERDICT_TXT), txt)
    return txt


def _load_prior(run):
    try:
        with open(os.path.join(run, VERDICT_JSON), encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def _base_record(run, args, state):
    return {"judge_run_version": JUDGE_VERSION, "state": state, "run": os.path.abspath(run), "run_rel": _rel(run),
            "judged_at": _now(),
            "invocation": {"argv": sys.argv[1:], "host": socket.gethostname(), "pid": os.getpid(),
                           "from_run_case": bool(getattr(args, "from_run_case", False))},
            "plan": "plans/active/tooling-convergence-and-wall-resolution-gates.md §4.6・§5.1 #12"}


def deferred_reason(run, args):
    if args.defer:
        return "--defer"
    if args.from_run_case:
        if os.path.exists(os.path.join(run, "stage_manifest.json")):
            return "段階起動の run (stage_manifest.json がある)"
        if os.environ.get("FORGE_JUDGE_DEFER", "0") not in ("", "0"):
            return "親が FORGE_JUDGE_DEFER=%s を渡した" % os.environ["FORGE_JUDGE_DEFER"]
    return None


def judge(run, args):
    """判定して RUN_VERDICT.{json,txt} を書く。戻り (終了コード, RUN_VERDICT.txt の本文)。"""
    run = os.path.abspath(run)
    prior = _load_prior(run)
    spec_info = load_spec(run)
    stop = read_stop_state(run)
    solver = check_solver(run, args.solver_rc, prior, stop)
    why_defer = deferred_reason(run, args)
    if why_defer:
        stage_checks = {"solver": solver, "nonfinite": check_nonfinite(run, residual_only_main=True), "physical": check_physical(run)}
        cv = os.path.join(run, "CONVERGENCE_VERDICT.txt")
        sub = None
        if os.path.exists(cv):
            txt = open(cv, encoding="utf-8", errors="replace").read()
            m = re.findall(r"^CONVERGENCE_EXIT:\s*(-?\d+)", txt, re.M)
            sub = {"exit": int(m[-1]) if m else None,
                   "verdict_lines": [ln.strip() for ln in txt.splitlines() if ln.startswith("=== ") and "->" in ln]}
        stages = list((prior or {}).get("stages") or []) if (prior or {}).get("state") == "deferred" else []
        stages.append({"seq": len(stages) + 1, "recorded_at": _now(), "checks": stage_checks, "convergence_sub": sub,
                       "config_sha256": _digest(os.path.join(run, "solverConfig.yaml")),
                       "residual_history": _digest(os.path.join(run, "residual_history.csv"))})
        st = _combine(c["status"] for s in stages for c in s["checks"].values())
        st = "FAIL" if st == "FAIL" else "INDETERMINATE"
        v = _base_record(run, args, "deferred")
        v.update(overall=st, state_note="段の結果だけを記録した (%s); 総合の判定は親が段の記録・停止の記録の確定の後に確定する" % why_defer,
                 overall_reason=("既知の FAIL がある段: %s" % [s["seq"] for s in stages
                                                       if any(c["status"] == "FAIL" for c in s["checks"].values())])
                 if st == "FAIL" else "確定待ち (段階起動)",
                 exit_code=EXIT_DEFERRED, execution=execution_info(run, solver), checks={}, stages=stages,
                 spec={k: spec_info.get(k) for k in ("source", "path", "sha256", "written_after_start", "error")},
                 tools=tool_hashes(), inputs=collect_inputs(run, spec_info))
        return EXIT_DEFERRED, write_verdict(run, v)
    cfg, cfg_err = _read_yaml(os.path.join(run, "solverConfig.yaml"))
    bcond, bcond_err = _read_yaml(os.path.join(run, "bcondConfig.yaml"))
    spec = spec_info.get("spec")
    checks = {"solver": solver,
              "nonfinite": check_nonfinite(run, full_audit=args.full_audit),
              "physical": check_physical(run)}
    if spec_info.get("error"):
        checks["convergence"] = _check("INDETERMINATE", spec_info["error"])
    else:
        checks["convergence"] = check_convergence(run, spec, cfg)
    checks["stop_rule"] = check_stop_rule(stop)
    checks["quasisteady"] = check_quasisteady(run, spec_info)
    checks["wall_resolution"] = check_wall_resolution(run, cfg, cfg_err, bcond, bcond_err, spec_info)
    st, why = overall_of(checks)
    v = _base_record(run, args, "final")
    v.update(overall=st, overall_reason=why, exit_code=EXIT_OF[st], execution=execution_info(run, solver), checks=checks,
             spec={k: spec_info.get(k) for k in ("source", "path", "sha256", "written_after_start", "error", "spec")},
             tools=tool_hashes(), inputs=collect_inputs(run, spec_info))
    if (prior or {}).get("stages"):
        # 段階起動の段の記録 (state deferred のときに run_case.sh が足したもの) は参考として引き継ぐ。総合への入れ方は親の確定 (後半) で決める
        v["stages"] = prior["stages"]
    return EXIT_OF[st], write_verdict(run, v)


def check_fresh(run):
    """判定 (RUN_VERDICT.json) が今の入力に対応するか → (終了コード, 本文)。0 対応 / 1 入力が変わった / 2 判定が無い。"""
    prior = _load_prior(run)
    if not prior or "inputs" not in prior:
        return 2, "NO VERDICT: %s が無いか読めない (%s)" % (VERDICT_JSON, run)
    now = collect_inputs(run)
    old = prior["inputs"]
    changed = sorted(k for k in set(now) | set(old) if now.get(k) != old.get(k))
    tools_now = tool_hashes()
    tool_changed = sorted(k for k in tools_now if (prior.get("tools") or {}).get(k) != tools_now[k])
    note = ("  (判定の後に道具が変わった: %s)" % ", ".join(tool_changed)) if tool_changed else ""
    if changed:
        return 1, "STALE: 判定 (%s, state %s, overall %s) の後に入力が変わった: %s%s" % (
            prior.get("judged_at"), prior.get("state"), prior.get("overall"), ", ".join(changed), note)
    return 0, "FRESH: 判定 (%s, state %s, overall %s) は今の入力に対応する%s" % (
        prior.get("judged_at"), prior.get("state"), prior.get("overall"), note)


def main(argv=None):
    ap = _Parser(description="計算の後の判定を 1 本に (RUN_VERDICT.json / RUN_VERDICT.txt)。詳しくは docstring。")
    ap.add_argument("run")
    ap.add_argument("--full-audit", action="store_true", help="非有限の監査を全 res_<N>.h5 と境界の出力にも (既定は最終の場だけ)")
    ap.add_argument("--solver-rc", type=int, default=None, help="forge の終了コード (run_case.sh が渡す)")
    ap.add_argument("--from-run-case", action="store_true",
                    help="run_case.sh からの呼び出し: 段階起動の run (stage_manifest.json がある / FORGE_JUDGE_DEFER=1) は段の結果だけ記録")
    ap.add_argument("--defer", action="store_true", help="段の結果だけを記録する (state deferred、総合は確定しない)")
    ap.add_argument("--check-fresh", action="store_true", help="判定が今の入力に対応するかを見る (書かない)")
    args = ap.parse_args(argv)
    run = os.path.abspath(args.run)
    if not os.path.isdir(run):
        print("judge_run: run が無い: %s" % args.run, file=sys.stderr)
        return EXIT_ERROR
    if args.check_fresh:
        rc, msg = check_fresh(run)
        print(msg)
        return rc
    try:
        rc, txt = judge(run, args)
        sys.stdout.write(txt)
        return rc
    except Exception as e:                                                        # noqa: BLE001
        tb = traceback.format_exc()
        print(tb, file=sys.stderr)
        try:
            v = _base_record(run, args, "error")
            v.update(overall="INDETERMINATE", overall_reason="判定器の例外 (%s: %s)" % (type(e).__name__, str(e)[:300]),
                     state_note="判定器の例外", exit_code=EXIT_ERROR, checks={}, error=tb,
                     tools=tool_hashes())
            sys.stdout.write(write_verdict(run, v))
        except Exception:                                                         # noqa: BLE001
            print("judge_run: RUN_VERDICT も書けなかった", file=sys.stderr)
        return EXIT_ERROR


if __name__ == "__main__":
    sys.exit(main())
