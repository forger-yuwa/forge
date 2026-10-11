"""⑤ SERN 評価の受理ゲート (plan §5.1 R1, codex C1 採用 2026-09-09)。

1 作動点 run を「サロゲート学習・Pareto に入れてよいか」を判定する。ゲートは互いに独立で、全て必須:
  (1) forge の終了コード rc == 0 (発散 run の採用は撤回。旧 §4.13-3 は廃止)
  (2) 保存場の有限値・正値 (最終 res_*.h5 の ro/roU/roe/P/T が有限、ro,P,T > 0、res_nan_*.h5 が無い)。量は乱流モデルから
      必須の集合を決め (SST は roK・roOmega も)、欠ければ判定不能 (plan tooling-sern-te-wake-grid §5.1 R2g。h5 にある量だけを
      検査して欠落を素通りさせない)。最終場の床の roOmega も SST では必須
  (3) 全残差 (check_convergence.analyze): NaN/Inf 無し、末尾で rising な列が無い。PASS (3 桁低下) は
      `require_residual_pass` のときだけ必須 (本ケースは残差が 1–2.5 桁でプラトーする性格 — その事実は台帳に残す)
  (4) **実際に最適化する量** (SST は C_T_with_shear、Euler は C_T) と C_T / C_L / C_M の頭打ち
      (check_quasisteady.classify_series による STEADY。NONFINITE / DRIFTING / TRANSIENT-UNSETTLED / OSCILLATING は不合格)。
      量と目的量は問題設定の乱流モデル (evaluate.model) から決め、量の列が力の履歴に欠ければ判定不能
      (plan tooling-sern-te-wake-grid §5.1 R2f。列の有無で量を減らさない)
  (5) 問題 YAML が要求するときだけ: 毎更新の床事象 (floor_event_gate)、判定区間の段の §6 #2・#3 (judgment_gate:
      後縁まわりの最低温度と 4 力量の窓。plan tooling-sern-te-wake-grid §5.1 R2)
判定結果は metrics.json の `gates` に残し、数値失敗 (fail_class) と物理的 INFEASIBLE (設計不成立・C_M 窓) を混同しない。
"""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import h5py
import numpy as np

from .sern_forces import steadiness

_FORGE_ROOT = Path(os.environ.get("FORGE_ROOT", Path(__file__).resolve().parents[3]))
_TOOLS = str(_FORGE_ROOT / "solver_density_cuda" / "tools")
if _TOOLS not in sys.path:
    sys.path.insert(0, _TOOLS)
import check_convergence as _cc  # noqa: E402
import check_floor_events as _cfe  # noqa: E402
from check_quasisteady import classify_series as _classify_series  # noqa: E402

FIELD_VARS = ("ro", "roUx", "roUy", "roUz", "roe", "P", "T", "roK", "roOmega")
POSITIVE_VARS = ("ro", "P", "T")
FORCE_KEYS = ("C_T", "C_L", "C_M")
# 乱流モデル → 最終場に必須の量 (plan tooling-sern-te-wake-grid §5.1 R2g、R2f の洗い出し 2026-10-12)。旧実装は FIELD_VARS のうち
# h5 にある量だけを検査したので、SST で roK・roOmega が無くても素通りした。必須の量が無ければ判定不能 (evaluate_gates で
# JUDGMENT_UNVERIFIABLE)。必須でない FIELD_VARS (Euler の roK・roOmega) はあれば従来どおり検査する
FIELD_REQUIRED = {"euler": ("ro", "roUx", "roUy", "roUz", "roe", "P", "T"),
                  "sst": ("ro", "roUx", "roUy", "roUz", "roe", "P", "T", "roK", "roOmega")}


def _volume_res(run_dir):
    return sorted(Path(run_dir).glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))


def _state_model(model) -> str:
    """場の必須量の表 (FIELD_REQUIRED) を引く乱流モデル。表の引き方は §6 #3 の必須量 (judgment_required_quantities) と同じ。
    None・未知は ValueError (collect が問題設定から解決して渡す。既定に読み替えない: R2f の _gate_quantities と同じ扱い)。"""
    if model is None:
        raise ValueError("場の必須量の乱流モデル (model) が None: collect は evaluate.model を解決して渡す")
    m = str(model).strip().lower()
    if m not in FIELD_REQUIRED:
        raise ValueError(f"乱流モデル {model!r} から場の必須量を決められない ({sorted(FIELD_REQUIRED)} のどれか)")
    return m


def field_health(run_dir, *, model) -> dict:
    """最終体積出力の保存量・原始量が有限で正値かを見る。res_nan_*.h5 (発散ダンプ) があれば不合格。
    model = 問題設定の乱流モデル (省略不可)。その必須の量 (FIELD_REQUIRED) が最終 h5 に無ければ判定不能 (missing に挙げて不合格)。"""
    m = _state_model(model)          # モデルの検査は run の読み込みより先に
    run_dir = Path(run_dir)
    nan_files = sorted(f.name for f in run_dir.glob("res_nan*.h5"))
    res = _volume_res(run_dir)
    out = {"ok": True, "file": None, "nonfinite": {}, "nonpositive": {}, "min": {}, "nan_files": nan_files, "reasons": [],
           "required": list(FIELD_REQUIRED[m]), "missing": []}
    if nan_files:
        out["ok"] = False; out["reasons"].append(f"res_nan dump present: {nan_files[:3]}")
    if not res:
        out["ok"] = False; out["reasons"].append("no volume res_*.h5"); return out
    out["file"] = res[-1].name
    with h5py.File(res[-1], "r") as f:
        # 固定の保存量に加え、**ファイルにある化学種・受動スカラーを全部**見る (codex plan-2 M3)。
        # FIELD_VARS だけだと roY* の NaN が field_health も floor_gate も素通りする実例があった。
        extra = tuple(sorted(k for k in f["VALUE"].keys() if re.fullmatch(r"(roY|Y|roXi|roQ)\d*[\d_]*", k)))
        out["species_checked"] = list(extra)
        out["missing"] = [v for v in FIELD_REQUIRED[m] if f"VALUE/{v}" not in f]
        for v in FIELD_VARS + extra:
            if f"VALUE/{v}" not in f:
                continue
            a = f[f"VALUE/{v}"][:]
            nbad = int(np.count_nonzero(~np.isfinite(a)))
            if nbad:
                out["nonfinite"][v] = nbad
            if v in POSITIVE_VARS:
                fin = a[np.isfinite(a)]
                out["min"][v] = float(fin.min()) if fin.size else None
                npos = int(np.count_nonzero(fin <= 0.0))
                if npos:
                    out["nonpositive"][v] = npos
    if out["nonfinite"]:
        out["ok"] = False; out["reasons"].append(f"non-finite field values: {out['nonfinite']}")
    if out["nonpositive"]:
        out["ok"] = False; out["reasons"].append(f"non-positive field values: {out['nonpositive']}")
    if out["missing"]:
        out["ok"] = False
        out["reasons"].append(f"判定不能: 必須の場の量 {', '.join(out['missing'])} が最終 h5 ({out['file']}) に無い (乱流モデル {m})")
    return out


def residual_health(run_dir, min_drop: float = 3.0, tail: float = 0.2) -> dict:
    """check_convergence.analyze で全残差列を見る (VERDICT の語彙も同じ)。
    ok = NaN/Inf 無し かつ rising 列無し。converged = 全列 PASS (3 桁低下)。"""
    path = Path(run_dir) / "residual_history.csv"
    out = {"ok": False, "converged": False, "verdict": "MISSING", "nan": False, "rising": [], "stalled": [], "converging": [],
           "columns": {}, "last_step": None}
    if not path.exists():
        out["reasons"] = ["no residual_history.csv"]; return out
    res = _cc.analyze(str(path), min_drop, tail)
    if res is None:
        out["reasons"] = ["empty residual_history.csv"]; return out
    laststep, report, ok, any_nan, any_stalled, any_converging = res
    out["last_step"] = laststep; out["nan"] = bool(any_nan); out["converged"] = bool(ok)
    # **trend は msg から直接読む** (2026-09-19, codex plan レビュー M1)。
    # 旧実装は `<-- RISING` 等の**注記文字列**で分類していたが、`check_convergence.py:107` は
    # 3 桁低下した列を col_ok = True として注記を付けない。そのため
    # **3 桁落ちてなお下降中 (falling) の列がどのリストにも入らず**、
    # 「全列プラトー」を要求したはずのゲートを素通りしていた。
    # msg は必ず `drop=X.Xdec <trend>` を含むので、そこから分類する。
    import re as _re
    for c, (msg, col_ok) in report.items():
        m = _re.search(r"drop=\s*([-\d.]+)dec\s+(rising|falling|flat)", msg)
        trend = m.group(2) if m else ""
        drop = float(m.group(1)) if m else float("nan")
        out["columns"][c] = {"msg": msg, "ok": bool(col_ok), "trend": trend, "drop": drop}
        if trend == "rising":
            out["rising"].append(c)
        elif trend == "flat":
            out["stalled"].append(c)
        elif trend == "falling":
            out["converging"].append(c)
    if ok:
        out["verdict"] = "PASS (converged)"
    elif any_nan:
        out["verdict"] = "DIVERGED (NaN/Inf)"
    elif out["rising"]:
        out["verdict"] = "NOT CONVERGED (rising)"
    elif any_stalled:
        out["verdict"] = "NOT CONVERGED (stalled/plateau)"
    elif any_converging:
        out["verdict"] = "NOT CONVERGED (still converging)"
    else:
        out["verdict"] = "NOT CONVERGED"
    out["ok"] = (not any_nan) and (not out["rising"])
    out["reasons"] = [] if out["ok"] else [f"residual {out['verdict']}" + (f" columns {out['rising']}" if out["rising"] else "")]
    return out


# --- 受理ゲート (4) の量 (plan tooling-sern-te-wake-grid §5.1 R2f、R2a の横展開 2026-10-12) --------------------------------
# 旧実装は「最後の履歴に C_T_with_shear の列があるか」で目的量と量の数を決めていたので、判定区間の段を持たない SST キャンペーンでは
# 列が欠けると C_T の 3 量で通った。量は問題設定の乱流モデルから決め (表は §6 #3 の必須量 JUDGMENT_QUANTITIES を共用)、
# 量の列が欠けた履歴は判定不能 (evaluate_gates で JUDGMENT_UNVERIFIABLE) にする。model は collect が問題設定 (evaluate.model、
# runner の _solver_config と同じ既定 euler) から解決して渡す。None は渡し忘れと区別できないので拒否する (既定に読み替えない)
def _gate_quantities(model) -> tuple:
    """(目的量, 量の並び)。並びは目的量が先頭 (旧実装の metrics.json と同じ)。None・未知のモデルは ValueError。"""
    if model is None:
        raise ValueError("受理ゲートの乱流モデル (model) が None: collect は evaluate.model を runner の既定 euler で解決して渡す")
    req = judgment_required_quantities(model)
    obj = "C_T_with_shear" if "C_T_with_shear" in req else "C_T"
    return obj, (obj,) + tuple(q for q in req if q != obj)


def objective_key(model) -> str:
    """実際に最適化する量: 乱流モデルが SST なら C_T_with_shear、Euler なら C_T (列の有無では決めない)。"""
    return _gate_quantities(model)[0]


def steadiness_gate(hist, model) -> dict:
    """目的量と力係数の頭打ち (受理ゲート (4))。量は乱流モデルから (SST 4 量・Euler 3 量)。量の列が履歴のどれかの標本に
    無ければ missing に挙げて不合格 (判定不能。evaluate_gates が JUDGMENT_UNVERIFIABLE にする)。値が非有限なら従来どおり NONFINITE。"""
    obj, keys = _gate_quantities(model)
    steps = [h["step"] for h in hist]
    missing = [k for k in keys if any(k not in h for h in hist)]
    series = {k: steadiness([h.get(k, float("nan")) for h in hist], steps=steps) for k in keys if k not in missing}
    unsteady = [k for k, v in series.items() if v["verdict"] != "STEADY"]
    reasons = ([f"判定不能: 必須量 {k} の列が力の履歴に無い ({sum(k not in h for h in hist)}/{len(hist)} 標本、乱流モデル {model})"
                for k in missing]
               + [f"{k}: {series[k]['verdict']} ({series[k]['detail']})" for k in unsteady])
    return {"ok": not unsteady and not missing, "objective": obj, "model": str(model).strip().lower(), "missing": missing,
            "series": series, "unsteady": unsteady, "reasons": reasons}



# --- 床・下限への張り付き (codex plan レビュー 2 の C1, 2026-09-19) ------------------------------
# 有限で非上昇なら通る、というだけのゲートは**床に張り付いた解**を受理してしまう。
# run_0122 では `sym` 上の 1 ノードで k=0・roOmega=1e-20 に張り付き、交差拡散が ω の下限 1e-12 で
# 割られて、その 1 点だけで rms_roOmega = 8.95e15 を作っていた (保存場は有限・力係数は ALL STEADY)。
FLOORS = {"P": ("pMin", 1.0), "T": ("tMin", 50.0), "ro": ("roMin", 1.0e-4)}
OMEGA_FLOOR = 1.0e-20          # update_d.cu の roOmega 下限
# ~~K_FLOOR~~ **削除 (2026-09-19, codex plan レビュー M6)**: `k <= 0` を一律に異常とするのは誤り。
# node 低 Re 壁は `sstNodeWallKPin` (既定 1) が**壁ノードの k を 0 にピンするのが正しい仕様**
# (`procedures/solver-settings.md` の SST 表)。run_0122 で「k≤0 が 28 ノード」と報告したのは
# その壁ピンを拾っていた可能性が高い。壁ピンを除いた領域で評価する仕組みが要る (残作業 R-f)。


def _solver_floors(run_dir) -> dict:
    """run の `solverConfig.yaml` から実効の床を読む (2026-09-19, codex plan レビュー M6)。
    書かれていなければソルバ既定 (`solverConfig.hpp`: pMin 1.0 / roMin 1e-4 / tMin 1e-4) を使う。

    **温度の実効下限は経路で違う** (2026-09-20, codex plan レビュー M2 で誤りが判明):

    - **単相 CPG** (`thermalMethod: 0` かつ凝縮なし) は `dependentVariables_d.cu` L290 の
      `T = max(intE/(cp/gamma), tMin)` で、床は **config の `tMin` (既定 1e-4 K)**。
    - **TP / 多成分 / 凝縮**の温度反転だけが `DEPVAR_TMIN` = **50 K** でクランプされる
      (同 L148/151/160)。

    旧実装は経路を見ずに一律 50 K を下限にしていたため、CPG run の 37.6 K を
    「床に張り付いた」と誤判定していた (`case/46` `run_0415`)。"""
    out = {"P": 1.0, "ro": 1.0e-4, "T": 1.0e-4}
    f = Path(run_dir) / "solverConfig.yaml"
    if f.exists():
        txt = f.read_text()
        m = re.search(r"thermalMethod\s*:\s*(\d+)", txt)
        tp = (m is not None and int(m.group(1)) != 0)
        mc = re.search(r"condensation\s*:\s*(\d+)", txt)
        if tp or (mc is not None and int(mc.group(1)) != 0):
            out["T"] = 50.0          # DEPVAR_TMIN: 温度反転を通る経路だけ
        for key, name in (("pMin", "P"), ("roMin", "ro"), ("tMin", "T")):
            m = re.search(rf"{key}\s*:\s*([-\d.eE+]+)", txt)
            if m:
                v = float(m.group(1))
                out[name] = max(v, out[name]) if name == "T" else v
    return out


def floor_gate(run_dir, p_min: float | None = None, tol: float = 1.0e-6, *, model) -> dict:
    """最終保存場で EOS 床・乱流下限に張り付いたノードを数える。1 個でも在れば NG。
    床は run の `solverConfig.yaml` から読む。**読み取れない・保存場が無いときは判定不能で不合格**
    (旧実装は例外を握り潰して ok=True を返していた: codex M6)。
    `k <= 0` は見ない — node 低 Re 壁は `sstNodeWallKPin` (既定 1) が正当に k=0 をピンする。
    model = 問題設定の乱流モデル (省略不可)。SST では roOmega の床を必須にし、h5 に無ければ P・T・ro が無いときと同じく
    「不在」(-1) で不合格 (plan tooling-sern-te-wake-grid §5.1 R2g。旧実装は roOmega があるときだけ数えた)。Euler は従来どおり
    roOmega があれば数える。"""
    import glob
    import os
    m = _state_model(model)          # モデルの検査は run の読み込みより先に (例外を判定不能に紛れさせない)
    run_dir = Path(run_dir)
    fs = sorted((q for q in glob.glob(str(run_dir / "res_[0-9]*.h5"))
                 if os.path.basename(q)[4:-3].isdigit()),      # 鏡像 (_full) など派生物を除く
                key=lambda q: int(os.path.basename(q)[4:-3]))
    if not fs:
        return {"ok": False, "counts": {}, "reasons": ["判定不能: res_*.h5 が無い"]}
    fl = _solver_floors(run_dir)
    if p_min is not None:
        fl["P"] = float(p_min)
    counts = {}
    try:
        import h5py
        with h5py.File(fs[-1], "r") as f:
            V = f["VALUE"]
            for name, lo in fl.items():
                if name not in V:
                    counts[f"{name} 不在"] = -1
                    continue
                a = np.asarray(V[name][:])
                counts[f"{name}<={lo:g}"] = int(np.sum(np.isfinite(a) & (a <= lo * (1.0 + tol))))
            # 化学種: 非負性と ΣρY = ρ (codex M6)
            ys = sorted(k for k in V.keys() if re.fullmatch(r"roY\d+", k))
            if ys and "ro" in V:
                roa = np.asarray(V["ro"][:]); tot = np.zeros_like(roa)
                neg = 0; nonfin = 0
                for k in ys:
                    a = np.asarray(V[k][:]); tot += a
                    # **非有限を先に数える** (codex plan-2 M3)。isfinite で絞ってから負値を見ると
                    # NaN が「負でない」として合格側に落ちる
                    nonfin += int(np.sum(~np.isfinite(a)))
                    neg += int(np.sum(np.isfinite(a) & (a < -tol * np.maximum(roa, 1e-30))))
                counts["roY 非有限"] = nonfin
                counts["roY<0"] = neg
                good = np.isfinite(tot) & np.isfinite(roa) & (roa > 0)
                # 全点判定不能なら -1 を返し、下の bad 判定で不合格にする
                counts["|sum(roY)/ro-1| max"] = float(np.max(np.abs(tot[good] / roa[good] - 1.0))) if good.any() else -1.0
            if "roOmega" in V:
                a = np.asarray(V["roOmega"][:])
                counts[f"roOmega<={OMEGA_FLOOR:g}"] = int(np.sum(np.isfinite(a) & (a <= OMEGA_FLOOR * (1.0 + tol))))
            elif m == "sst":
                counts["roOmega 不在"] = -1          # SST では必須 (無ければ床を数えられない = 不合格)
            # k の床判定は壁ピン (正当) と区別できないので**当面外す** (R-f)
    except Exception as e:                                 # **読めないときは判定不能で不合格**
        return {"ok": False, "counts": {}, "reasons": [f"判定不能: {type(e).__name__}: {e}"]}
    # `|sum...|` は許容超過 **と 判定不能 (-1)** の両方を不合格にする (codex plan-2 M3)
    bad = {k: v for k, v in counts.items()
           if (k.startswith("|sum") and (v > 1.0e-4 or v < 0.0)) or (not k.startswith("|sum") and v != 0)}
    return {"ok": not bad, "counts": counts, "file": os.path.basename(fs[-1]),
            "reasons": [] if not bad else ["床/下限に張り付き: " + ", ".join(f"{k} {v} ノード" for k, v in bad.items())]}


# --- 毎更新の EOS 床事象 (plan tooling-sern-te-wake-grid §4「床の判定」、codex plan レビュー M1・M2) -------------
# `floor_gate()` は最後の HDF5 だけを見るので、判定区間中に床を使い最終場で回復した評価が通る。生産系列は forge の
# `output: {floorEvents: 1}` が毎更新に書く floor_events.csv を必須入力にし、判定区間の全更新で床事象 0 を要求する。
# 判定の正本は solver_density_cuda/tools/check_floor_events.py (記録の欠落・区間の欠落・未検証の経路は判定不能 = 不合格)。
# 有効にする経路: 問題 YAML の `evaluate.floor_events: 1` (runner が config に書き、collect がゲートを必須にする)。
# 判定区間は `evaluate.floor_events_window_steps` (末尾の step 数)。省略時は最後の session の末尾 50 % (plan §6)。
FLOOR_EVENTS_KEY = "floor_events"
FLOOR_EVENTS_WINDOW_KEY = "floor_events_window_steps"


def floor_events_options(evaluate) -> dict:
    """問題 YAML の evaluate 節から evaluate_gates の床事象ゲートの引数を作る。キーが無ければ必須にしない (既存の挙動)。"""
    ev = evaluate or {}
    req = ev.get(FLOOR_EVENTS_KEY, 0)
    if req not in (0, 1, True, False, None):
        raise ValueError(f"evaluate.{FLOOR_EVENTS_KEY} は 0 か 1: {req!r}")
    w = ev.get(FLOOR_EVENTS_WINDOW_KEY)
    if w is not None and (not isinstance(w, int) or isinstance(w, bool) or w <= 0):
        raise ValueError(f"evaluate.{FLOOR_EVENTS_WINDOW_KEY} は正の整数 [step]: {w!r}")
    return {"require_floor_events": bool(req), "floor_window_steps": w}


def floor_events_config_line(evaluate) -> str:
    """solverConfig.yaml に足す行 (evaluate.floor_events: 1 のとき `output: {floorEvents: 1}`、それ以外は空)。"""
    return "output: {floorEvents: 1}\n" if floor_events_options(evaluate)["require_floor_events"] else ""


def floor_event_gate(run_dir, required: bool = False, window_steps: int | None = None, tail: float = 0.5) -> dict:
    """floor_events.csv の判定 (check_floor_events.judge)。`required=False` (既定) では結果を載せるだけで合否に使わない。
    必須のときは PASS 以外を不合格にする: 区間内の床事象は FLOOR_EVENT、判定不能 (記録の欠落・区間の欠落・restart 後の
    区間の不足・未検証の経路・監査なし・食い違い) は FLOOR_UNVERIFIABLE。"""
    res = _cfe.judge(run_dir, window_steps=window_steps, tail=tail)
    fc = None if res["ok"] else ("FLOOR_EVENT" if res["verdict"] == "FAIL" else "FLOOR_UNVERIFIABLE")
    return {"ok": bool(res["ok"] or not required), "required": bool(required), "verdict": res["verdict"],
            "fail_class": fc if required else None, "window": res["window"], "path": res["path"],
            "last_step": res["last_step"], "session": res["session"], "n_updates_checked": res["n_updates_checked"],
            "events": res["events"], "first_events": res["first_events"], "entry": res["entry"], "init": res["init"],
            "aux": res["aux"], "near_window": res["near_window"], "ghost_window": res["ghost_window"],
            "reasons": list(res["reasons"]) if (required and not res["ok"]) else [],
            "notes": [] if (required and not res["ok"]) else list(res["reasons"])}


def floor_events_summary(fev: dict | None) -> dict | None:
    """metrics.json の上位と台帳に載せる要約 (必須か・判定・区間・種類別の件数・理由)。"""
    if not fev:
        return None
    return {"required": fev.get("required"), "verdict": fev.get("verdict"), "fail_class": fev.get("fail_class"),
            "window": fev.get("window"), "n_updates_checked": fev.get("n_updates_checked"),
            "n_events": {k: v.get("n_events") for k, v in (fev.get("events") or {}).items()},
            "first_events": (fev.get("first_events") or [])[:2],
            "reasons": (fev.get("reasons") or fev.get("notes") or [])[:3]}


# --- 判定区間の §6 #2・#3 (plan tooling-sern-te-wake-grid §5.1 R2、codex result レビュー 2026-10-11 M2) -------------
# 生産の標準経路 (runner の段階起動 → 判定区間の段) で、検証の手作業 (case/46.sern_design/diag/tewake_*.sh) が使った判定を
# 保存データから再計算できる形で持つ (collect と --rejudge が同じ関数を通る)。#1 (床事象) は floor_event_gate。
#   #3 窓: collect の力の履歴の末尾 W step の平均 m・a = max|q − m|・前窓 (その前の W step) の平均との差 d
#          (case/46.sern_design/v3_farfield_eval.py と同じ切り方)。max(a, d) ≤ 許容 (C_T・C_T_with_shear・C_L 5e-5、C_M 5e-4)、
#          かつ必須量 (SST 4 量・Euler 3 量) が check_quasisteady --tail 0.5 (drift 0.05・osc 0.10・min 4 の既定) で STEADY。
#          力の履歴の step は予定 step 列 (両窓とも 20 標本以上) と一致すること
#   #2 温度: 監視 CSV (diag/te_monitor.py が全場のスナップショットごとに 1 行) の R_TE・R_SE・RET・ALL で、末尾 W step に
#          T < 150 K の節点が無い・各領域の最低温度が --tail 0.5 で STEADY。T < 200 K の数と体積は記録だけ
# 判定: PASS / UNSETTLED (窓・定常性が未達 = 同一設定の延長の対象) / FAIL (低温・非有限の温度) / INDETERMINATE (標本の欠け・
# 空の領域・履歴が 2 窓に足りない・非有限の力 = 判定不能、延長の対象外)。**許容と窓は事前登録の値で固定** (結果を見て緩めない)
# **必須量と標本の格子は問題設定から決める** (codex diagnose 2026-10-12 M1、plan §5.1 R2a): 必須量は乱流モデル (evaluate.model) で
# 決め (SST は 4 量、Euler は 3 量。列の有無で減らさない)、標本は判定区間の段の予定 step 列 (N・W・outStepInterval) と照合する
# (欠測・重複・予定外・終了 step の不一致は判定不能。標本数だけ見ると、各窓の後半だけに 20 標本を置いた入力が通っていた)
JUDGMENT_FORCE_TOL = {"C_T": 5.0e-5, "C_T_with_shear": 5.0e-5, "C_L": 5.0e-5, "C_M": 5.0e-4}
JUDGMENT_MIN_SAMPLES = 20
JUDGMENT_STEADY = {"tail": 0.5, "drift": 0.05, "osc": 0.10, "min_snaps": 4}    # check_quasisteady.py --tail 0.5 (他は既定)
# 乱流モデル → §6 #3 の必須量。runner の _solver_config と同じ既定 (evaluate.model が無ければ euler)。受理ゲート (4) の量も
# この表から (_gate_quantities、§5.1 R2f)
JUDGMENT_QUANTITIES = {"sst": ("C_T", "C_T_with_shear", "C_L", "C_M"), "euler": ("C_T", "C_L", "C_M")}
TE_REGIONS = ("R_TE", "R_SE", "RET", "ALL")
T_COLD = 150.0


def _steady_verdict(steps, vals) -> tuple:
    s = JUDGMENT_STEADY
    v, detail, _ = _classify_series(steps, vals, s["tail"], s["drift"], s["osc"], s["min_snaps"])
    return v, detail


def judgment_required_quantities(model) -> tuple:
    """問題設定の乱流モデル (evaluate.model。None は runner と同じ既定 euler) から §6 #3 の必須量。未知のモデルは ValueError
    (必須量を決められないまま判定しない)。"""
    m = "euler" if model is None else str(model).strip().lower()
    if m not in JUDGMENT_QUANTITIES:
        raise ValueError(f"evaluate.model {model!r} から §6 #3 の必須量を決められない ({sorted(JUDGMENT_QUANTITIES)} のどれか)")
    return JUDGMENT_QUANTITIES[m]


def judgment_planned_steps(n_steps, window_steps, out_interval) -> list:
    """判定区間の段の予定の標本 step 列 = outStepInterval の倍数 (0 < s ≤ N)。forge は段ごとに step 0 から数え、壁面の出力
    (力の履歴の元) は outStepInterval の倍数の step だけ (step 0 は書かない、output.cpp の outputBconds_H5_XDMF)。
    N ≥ 2W・outStepInterval が W と N を割り切り・窓あたり JUDGMENT_MIN_SAMPLES 以上でなければ ValueError
    (runner_sern.judgment_config が config に課す条件と同じ。窓を短くして通さない)。"""
    vals = (n_steps, window_steps, out_interval)
    if not all(isinstance(v, (int, np.integer)) and not isinstance(v, bool) and v > 0 for v in vals):
        raise ValueError(f"判定区間の予定: N・W・outStepInterval は正の整数 ({vals!r})")
    N, W, oi = (int(v) for v in vals)
    if N < 2 * W or W % oi or N % oi or W // oi < JUDGMENT_MIN_SAMPLES:
        raise ValueError(f"判定区間の予定: N {N}・W {W}・outStepInterval {oi} は N ≥ 2W・outStepInterval が W と N を割り切り・"
                         f"窓あたり {JUDGMENT_MIN_SAMPLES} 標本以上を満たさない")
    return list(range(oi, N + 1, oi))


def _planned_of(planned_steps) -> list | None:
    """予定 step 列を整数のリストにする (空・整数でない・昇順でない・重複は None = 予定として使えない)。"""
    try:
        p = [int(s) for s in planned_steps]
    except (TypeError, ValueError):
        return None
    if not p or any(float(s) != float(q) for s, q in zip(planned_steps, p)) or any(b <= a for a, b in zip(p, p[1:])):
        return None
    return p


def _step_mismatch(steps, planned, what: str) -> list:
    """観測 step 列と予定 step 列の照合: 欠測・重複・予定外・終了 step の不一致・昇順でないことを理由の文字列で返す (一致なら空)。"""
    from collections import Counter
    out = []
    obs = [float(s) for s in steps]
    cnt = Counter(obs)
    dup = sorted(s for s, c in cnt.items() if c > 1)
    plan = set(float(s) for s in planned)
    missing = sorted(s for s in plan if s not in cnt)
    extra = sorted(s for s in cnt if s not in plan)

    def _fmt(xs):
        return ", ".join(f"{x:g}" for x in xs[:5]) + (" …" if len(xs) > 5 else "")
    if obs and max(obs) != float(planned[-1]):
        out.append(f"判定不能: {what}の最終 step {max(obs):g} が予定の終了 step {planned[-1]} と違う")
    if dup:
        out.append(f"判定不能: {what}の step が重複している ({_fmt(dup)})")
    if missing:
        out.append(f"判定不能: {what}に予定の標本が欠けている ({len(missing)} 個: {_fmt(missing)})")
    if extra:
        out.append(f"判定不能: {what}に予定外の step がある ({len(extra)} 個: {_fmt(extra)})")
    if any(b < a for a, b in zip(obs, obs[1:])):
        out.append(f"判定不能: {what}の step が昇順でない")
    return out


def judgment_window_gate(hist, window_steps: int, required, planned_steps) -> dict:
    """§6 #3 (窓と必須量の STEADY)。hist は collect の力の履歴 (step 昇順)。required = 必須量 (judgment_required_quantities)、
    planned_steps = 判定区間の段の予定 step 列 (judgment_planned_steps)。必須量の列の欠落・予定 step 列との不一致は判定不能。"""
    W = int(window_steps)
    req = tuple(required or ())
    planned = _planned_of(planned_steps or ())
    out = {"verdict": "INDETERMINATE", "ok": False, "window": None, "n_tail": 0, "n_prev": 0, "quantities": {}, "reasons": [],
           "required": list(req), "planned": {"end": planned[-1] if planned else None, "n": len(planned) if planned else 0}}
    if not req or any(q not in JUDGMENT_FORCE_TOL for q in req) or not {"C_T", "C_L", "C_M"} <= set(req):
        out["reasons"] = [f"判定不能: 必須量 {list(req)} が §6 #3 の量 ({list(JUDGMENT_FORCE_TOL)}) の C_T・C_L・C_M を含む部分集合でない"]
        return out
    if planned is None:
        out["reasons"] = ["判定不能: 予定 step 列が無い・整数の昇順でない"]
        return out
    N = planned[-1]
    out["window"] = [N - W, N]
    if N < 2 * W:
        out["reasons"] = [f"判定不能: 予定の終了 step {N} が 2 窓 ({2 * W} step) に足りない (前窓との差 d が取れない)"]
        return out
    if not hist:
        out["reasons"] = ["判定不能: 力の履歴が無い"]
        return out
    steps = np.array([float(h["step"]) if h.get("step") is not None else np.nan for h in hist])
    if not np.all(np.isfinite(steps)):
        out["reasons"] = ["判定不能: 力の履歴の step が読めない"]
        return out
    mism = _step_mismatch(steps, planned, "力の履歴")
    if mism:
        out["reasons"] = mism
        return out
    tail = steps > N - W
    prev = (steps > N - 2 * W) & ~tail
    out["n_tail"], out["n_prev"] = int(tail.sum()), int(prev.sum())
    if out["n_tail"] < JUDGMENT_MIN_SAMPLES or out["n_prev"] < JUDGMENT_MIN_SAMPLES:
        out["reasons"] = [f"判定不能: 窓の標本が足りない (末尾 {out['n_tail']}・前窓 {out['n_prev']} < {JUDGMENT_MIN_SAMPLES})"]
        return out
    # 必須量は問題設定から (列の有無で減らさない)。列が無い量は判定不能
    lack = [q for q in req if any(q not in h for h in hist)]
    if lack:
        out["reasons"] = [f"判定不能: 必須量 {q} の列が力の履歴に無い ({sum(q not in h for h in hist)}/{len(hist)} 標本)" for q in lack]
        return out
    qs = list(req)
    bad, unsettled = [], []
    for q in qs:
        v = np.array([float(h[q]) if h.get(q) is not None else np.nan for h in hist])
        if not np.all(np.isfinite(v)):
            bad.append(f"{q} に非有限・欠けた値 {int(np.count_nonzero(~np.isfinite(v)))} 個")
            continue
        m = float(v[tail].mean())
        a = float(np.max(np.abs(v[tail] - m)))
        d = float(abs(m - v[prev].mean()))
        st, detail = _steady_verdict(steps, v)
        tol = JUDGMENT_FORCE_TOL[q]
        ok = max(a, d) <= tol and st == "STEADY"
        out["quantities"][q] = {"mean": m, "a": a, "d": d, "tol": tol, "steady": st, "detail": detail, "ok": bool(ok)}
        if not ok:
            unsettled.append(f"{q}: a {a:.2e}・d {d:.2e} (許容 {tol:.0e})・{st}")
    if bad:
        out["reasons"] = ["判定不能: " + b for b in bad]
        return out
    out["verdict"] = "UNSETTLED" if unsettled else "PASS"
    out["ok"] = not unsettled
    out["reasons"] = ["窓・定常性が未達: " + u for u in unsettled]
    return out


def read_monitor_csv(path) -> list | None:
    """監視 CSV (te_monitor.py の出力) を読む。無ければ None。数値にできない欄は NaN。step 昇順に並べる。"""
    import csv
    p = Path(path)
    if not p.exists():
        return None
    rows = []
    with open(p, newline="") as fh:
        for r in csv.DictReader(fh):
            row = {}
            for k, v in r.items():
                try:
                    row[k] = float(v)
                except (TypeError, ValueError):
                    row[k] = float("nan")
            rows.append(row)
    return sorted(rows, key=lambda r: r.get("step", float("nan")))


def judgment_temperature_gate(rows, window_steps: int, last_step: int | None, expect_steps=None) -> dict:
    """§6 #2 (端の近くの低温)。rows = read_monitor_csv の戻り。last_step = 判定区間の段の予定の終了 step (監視も同じ場で終わること)、
    expect_steps = 判定区間の段の予定 step 列 (末尾窓で欠測・重複・予定外があれば判定不能)。"""
    W = int(window_steps)
    out = {"verdict": "INDETERMINATE", "ok": False, "window": None, "n_window": 0, "regions": {}, "reasons": []}
    if rows is None:
        out["reasons"] = ["判定不能: 監視の記録 (TE_MONITOR.csv) が無い"]
        return out
    steps = np.array([r.get("step", np.nan) for r in rows])
    if not rows or not np.all(np.isfinite(steps)):
        out["reasons"] = ["判定不能: 監視の記録が空・step が読めない"]
        return out
    N = int(steps.max())
    out["window"] = [N - W, N]
    if last_step is not None and N != int(last_step):
        out["reasons"] = [f"判定不能: 監視の最終 step {N} が予定の終了 step {last_step} と違う (最後の場を監視していない)"]
        return out
    win = steps > N - W
    out["n_window"] = int(win.sum())
    reasons = []
    if expect_steps is not None:          # 末尾窓の予定 step 列と照合 (欠測・重複・予定外)
        plan_w = sorted({int(s) for s in expect_steps if N - W < int(s) <= N})
        reasons += _step_mismatch(steps[win], plan_w, "末尾窓の監視") if plan_w else ["判定不能: 末尾窓に予定の標本が無い"]
    if out["n_window"] < JUDGMENT_MIN_SAMPLES:
        reasons.append(f"判定不能: 末尾窓の監視の標本 {out['n_window']} < {JUDGMENT_MIN_SAMPLES}")
    nonfin = float(np.sum([r.get("nonfinite_T", np.nan) for r, w in zip(rows, win) if w]))
    cold, unsettled = [], []
    for R in TE_REGIONS:
        n = np.array([r.get(f"{R}_n", np.nan) for r in rows])
        if not np.all(np.isfinite(n[win])) or np.any(n[win] <= 0):
            reasons.append(f"判定不能: 領域 {R} が空・欠けている")
            continue
        tmin = np.array([r.get(f"{R}_Tmin", np.nan) for r in rows])
        lt150 = np.array([r.get(f"{R}_n_lt150", np.nan) for r in rows])
        lt200 = np.array([r.get(f"{R}_n_lt200", np.nan) for r in rows])
        vol200 = np.array([r.get(f"{R}_vol_lt200", np.nan) for r in rows])
        if not (np.all(np.isfinite(tmin[win])) and np.all(np.isfinite(lt150[win]))):
            reasons.append(f"判定不能: 領域 {R} の最低温度・低温の数が読めない")
            continue
        st, detail = _steady_verdict(steps, tmin)
        jmin = int(np.flatnonzero(win)[np.argmin(tmin[win])])
        reg = {"Tmin_window": float(tmin[jmin]), "Tmin_step": int(steps[jmin]), "Tmin_node": rows[jmin].get(f"{R}_Tmin_node"),
               "Tmin_xyz": [rows[jmin].get(f"{R}_Tmin_{c}") for c in "xyz"], "n_lt150_window_sum": int(np.sum(lt150[win])),
               "n_lt200_window_max": int(np.nanmax(lt200[win])) if np.any(np.isfinite(lt200[win])) else None,
               "vol_lt200_window_max": float(np.nanmax(vol200[win])) if np.any(np.isfinite(vol200[win])) else None,
               "steady": st, "detail": detail}
        out["regions"][R] = reg
        if reg["n_lt150_window_sum"] > 0:
            cold.append(f"{R}: T < {T_COLD:g} K が延べ {reg['n_lt150_window_sum']} 節点 (最低 {reg['Tmin_window']:.1f} K、step {reg['Tmin_step']})")
        if st != "STEADY":
            unsettled.append(f"{R} の最低温度 {st}")
    if not (np.isfinite(nonfin) and nonfin == 0):
        cold.append(f"末尾窓に非有限の温度 (延べ {nonfin:g})")
    if cold:
        out["verdict"] = "FAIL"
        out["reasons"] = ["低温: " + c for c in cold] + reasons
        return out
    if reasons:
        out["reasons"] = reasons
        return out
    out["verdict"] = "UNSETTLED" if unsettled else "PASS"
    out["ok"] = not unsettled
    out["reasons"] = ["最低温度が未定常: " + u for u in unsettled]
    return out


def judgment_gate(hist, window_steps: int, required, planned_steps, monitor_rows=None, monitor_unavailable: str | None = None) -> dict:
    """§6 #2・#3 をまとめる。fail_class: COLD_SPOT (低温・非有限の温度。延長の対象外) / JUDGMENT_UNVERIFIABLE (判定不能。対象外) /
    JUDGMENT_UNSETTLED (窓・定常性の未達。同一設定の延長の対象)。monitor_unavailable を渡すと #2 は判定不能 (監視が使えない理由)。
    required・planned_steps は問題設定から作った必須量と予定 step 列 (judgment_plan)。力の窓・温度の監視とも予定 step 列と照合する。"""
    win = judgment_window_gate(hist, window_steps, required, planned_steps)
    planned = _planned_of(planned_steps or ())
    last = planned[-1] if planned else None          # 予定の終了 step (予定が使えなければ None。力の窓が判定不能になる)
    if monitor_unavailable:
        temp = {"verdict": "INDETERMINATE", "ok": False, "window": None, "n_window": 0, "regions": {},
                "reasons": [f"判定不能: 温度の監視が使えない ({monitor_unavailable})"]}
    else:
        temp = judgment_temperature_gate(monitor_rows, window_steps, last, expect_steps=planned or [])
    vs = (win["verdict"], temp["verdict"])
    if temp["verdict"] == "FAIL":
        fc, ext = "COLD_SPOT", False
    elif "INDETERMINATE" in vs:
        fc, ext = "JUDGMENT_UNVERIFIABLE", False
    elif "UNSETTLED" in vs:
        fc, ext = "JUDGMENT_UNSETTLED", True
    else:
        fc, ext = None, False
    verdict = {None: "PASS", "COLD_SPOT": "FAIL", "JUDGMENT_UNVERIFIABLE": "INDETERMINATE", "JUDGMENT_UNSETTLED": "UNSETTLED"}[fc]
    reasons = [] if fc is None else [r for g in (win, temp) for r in g["reasons"]]
    return {"ok": fc is None, "verdict": verdict, "fail_class": fc, "extendable": ext, "window_steps": int(window_steps),
            "force_window": win, "temperature": temp, "reasons": reasons}


def judgment_summary(jg: dict | None) -> dict | None:
    """台帳の作動点要約に載せる判定区間の §6 #2・#3 の要約 (plan §5.1 R2b: driver_sern の学習の採否が判定・窓の長さ・
    必須量の充足・温度の領域を照合する)。judgment_gate を通っていない評価は None。"""
    if not jg:
        return None
    fw, tp = jg.get("force_window") or {}, jg.get("temperature") or {}
    return {"verdict": jg.get("verdict"), "fail_class": jg.get("fail_class"), "window_steps": jg.get("window_steps"),
            "force_window": {"verdict": fw.get("verdict"), "window": fw.get("window"), "n_tail": fw.get("n_tail"), "n_prev": fw.get("n_prev"),
                             "required": fw.get("required"), "planned": fw.get("planned"),
                             "ok": {q: bool(v.get("ok")) for q, v in (fw.get("quantities") or {}).items()}},
            "temperature": {"verdict": tp.get("verdict"), "window": tp.get("window"), "n_window": tp.get("n_window"),
                            "regions": sorted(tp.get("regions") or {})},
            "reasons": (jg.get("reasons") or [])[:3]}


def judgment_extendable(gates: dict) -> tuple:
    """同一設定の延長 (+N を 1 回) の対象か。**未達が「未定常」だけ**のとき (既存ゲートの UNSTEADY と判定区間の
    JUDGMENT_UNSETTLED) に限る。rc ≠ 0・場・残差・床 (最終場・毎更新)・低温・判定不能があれば延長しない (戻り = (可否, 理由))。"""
    if gates.get("verdict") == "PASS":
        return False, "合格"
    blockers = []
    rc = gates.get("rc")
    if rc != 0:
        blockers.append(f"rc={rc}")
    for key, label in (("field", "保存場"), ("floors", "最終場の床"), ("floor_events", "床事象"), ("residual_scale", "残差の桁")):
        if not (gates.get(key) or {}).get("ok", False):
            blockers.append(label)
    res = gates.get("residual") or {}
    if not res.get("ok", False) or (gates.get("require_residual_pass") and not res.get("converged")):
        blockers.append("残差")
    if gates.get("require_residual_plateau") and (gates.get("plateau") or {}).get("falling"):
        blockers.append("残差のプラトー")
    stead = gates.get("steadiness") or {}
    if not stead.get("series"):
        blockers.append("力の履歴")
    if stead.get("missing"):             # 量の列の欠落は判定不能 (延長しても列は増えない。§5.1 R2f)
        blockers.append("必須量の列 " + "・".join(stead["missing"]))
    jg = gates.get("judgment")
    if jg is not None and not jg.get("ok") and not jg.get("extendable"):
        blockers.append(f"判定区間 {jg.get('fail_class')}")
    if blockers:
        return False, "未定常以外の不合格: " + "・".join(blockers)
    unsteady = (not stead.get("ok", False)) or (jg is not None and jg.get("fail_class") == "JUDGMENT_UNSETTLED")
    return (True, "未定常 (窓・定常性の未達) だけ") if unsteady else (False, "不合格の理由が分類できない")


def residual_scale_gate(run_dir, ratio: float = 1.0e6, gate: bool = False) -> dict:
    """残差列の**桁の揃い**を見る **補助警報** (2026-09-19, codex plan レビュー M3 で受理条件から降格)。

    次元の違う量 (質量・運動量・エネルギー・乱流・化学種) を同じ集合に入れて中央値と比べているので、
    **収束や方程式間の釣り合いを意味しない** (`residualMonitor_d.cu` の RMS に方程式間の正規化は無い)。
    極端な異常 (run_0122 の `rms_roOmega` 8.95e15 / 他は 1e-4〜1e-1) の検知には有効なので、
    `gate=True` のときだけ不合格にする。既定は警告のみ。"""
    import csv
    path = Path(run_dir) / "residual_history.csv"
    if not path.exists():
        return {"ok": True, "skipped": "no residual_history.csv"}
    rows = [r for r in csv.DictReader(path.open()) if r.get("phase") == "outer_end"]
    if not rows:
        return {"ok": True, "skipped": "no outer_end rows"}
    last = {}
    for k, v in rows[-1].items():
        if not k.startswith("rms_"):
            continue
        try:
            fv = float(v)
        except (TypeError, ValueError):
            continue
        if np.isfinite(fv) and fv > 0.0:
            last[k] = fv
    if len(last) < 2:
        return {"ok": True, "skipped": "columns < 2"}
    med = float(np.median(list(last.values())))
    bad = {k: v for k, v in last.items() if v > med * ratio}
    return {"ok": (not bad) or (not gate), "gate": bool(gate), "median": med, "outliers": bad,
            "warnings": [] if not bad else [f"残差の桁が揃わない (中央値 {med:.2e}): "
                                            + ", ".join(f"{k} {v:.2e}" for k, v in bad.items())],
            "reasons": [] if (not bad or not gate) else
            [f"残差の桁が揃わない (中央値 {med:.2e}): " + ", ".join(f"{k} {v:.2e}" for k, v in bad.items())]}


def evaluate_gates(run_dir, hist, rc, require_residual_pass: bool = False,
                   p_min: float | None = None, require_residual_plateau: bool = False,
                   require_floor_events: bool = False, floor_window_steps: int | None = None, floor_tail: float = 0.5,
                   judgment: dict | None = None, *, model) -> dict:
    """全ゲートを評価して verdict / fail_class を返す。fail_class は数値失敗の種別:
    DIVERGED (rc≠0 / 発散ダンプ / 非有限・非正の場 / 残差 NaN), RESIDUAL_RISING, NOT_CONVERGED (require 時のみ),
    NO_FORCES (壁出力が無い), UNSTEADY (目的量・力係数が頭打ちしていない),
    FLOOR_EVENT / FLOOR_UNVERIFIABLE (`require_floor_events` のときだけ: 判定区間の毎更新の床事象 / その記録が判定不能),
    COLD_SPOT / JUDGMENT_UNVERIFIABLE / JUDGMENT_UNSETTLED (`judgment` [judgment_gate の戻り] を渡したときだけ: 判定区間の §6 #2・#3)。
    JUDGMENT_UNVERIFIABLE は受理ゲート (4) の量の列が力の履歴に欠けたときにも出す (判定区間の有無によらない、§5.1 R2f)。
    最終場に乱流モデルの必須の量 (SST の roK・roOmega など) が無く、場の不合格がその欠落だけのときも JUDGMENT_UNVERIFIABLE (§5.1 R2g)。
    model = 問題設定の乱流モデル (collect が evaluate.model を runner の既定 euler で解決して渡す。省略不可・None 不可):
    目的量と (4) の量、場と最終場の床の必須の量はここから決める (旧引数 obj と、列・量の有無から推測する経路は廃止)。
    物理的 INFEASIBLE はここでは出さない。"""
    obj = objective_key(model)          # モデルの検査 (None・未知は ValueError) を履歴・run の読み込みより先に
    field = field_health(run_dir, model=model)
    resid = residual_health(run_dir)
    floors = floor_gate(run_dir, p_min, model=model)
    fevents = floor_event_gate(run_dir, required=require_floor_events, window_steps=floor_window_steps, tail=floor_tail)
    rscale = residual_scale_gate(run_dir)
    stead = steadiness_gate(hist, model) if hist else {"ok": False, "objective": obj, "model": str(model).strip().lower(), "missing": [],
                                                      "series": {}, "unsteady": [], "reasons": ["no force history (no wall output)"]}
    reasons = []; fail = None
    if rc is None:
        reasons.append("forge rc unknown"); fail = fail or "DIVERGED"
    elif rc != 0:
        reasons.append(f"forge rc={rc}"); fail = fail or "DIVERGED"
    if not field["ok"]:
        # 必須の量の欠落**だけ**なら判定不能 (延長しても量は増えない)。非有限・非正・発散ダンプ・場が無いは従来どおり DIVERGED
        only_missing = bool(field.get("missing")) and not (field["nonfinite"] or field["nonpositive"] or field["nan_files"])
        reasons += field["reasons"]; fail = fail or ("JUDGMENT_UNVERIFIABLE" if only_missing else "DIVERGED")
    if not resid["ok"]:
        reasons += resid.get("reasons", []); fail = fail or ("DIVERGED" if resid["nan"] else "RESIDUAL_RISING")
    elif require_residual_pass and not resid["converged"]:
        reasons.append(f"residual {resid['verdict']} (require_residual_pass)"); fail = fail or "NOT_CONVERGED"
    elif require_residual_plateau and resid.get("converging"):
        # **プラトー要求は既定 OFF に降格** (2026-09-19, codex plan レビュー M1)。
        # プラトーは収束の**十分条件ではない**: 残差の大きさを問わないので、float32 の更新消失や
        # クランプで動かなくなった状態と、方程式を満たして止まった状態を区別できない。
        # 受理には方程式別の無次元残差上限・保存収支・場/目的量の定常性が要る (残作業 R-a)。
        # 停滞の**診断情報**としては有用なので opt-in で残す。
        reasons.append(f"residual まだ低下中 (プラトー未達) columns {resid['converging']}")
        fail = fail or "NOT_PLATEAU"
    if not floors["ok"]:
        reasons += floors["reasons"]; fail = fail or "FLOOR_STUCK"
    if not fevents["ok"]:
        reasons += ["床事象: " + r for r in fevents["reasons"]]; fail = fail or fevents["fail_class"]
    if not rscale["ok"]:
        reasons += rscale["reasons"]; fail = fail or "RESIDUAL_UNBALANCED"
    if not hist:
        reasons += stead["reasons"]; fail = fail or "NO_FORCES"
    elif stead["missing"]:
        # 量の列の欠落は判定不能 (UNSTEADY にすると「未定常だけ」として延長の対象になる)
        reasons += stead["reasons"]; fail = fail or "JUDGMENT_UNVERIFIABLE"
    elif not stead["ok"]:
        reasons += stead["reasons"]; fail = fail or "UNSTEADY"
    if judgment is not None and not judgment.get("ok"):
        reasons += ["判定区間: " + r for r in judgment.get("reasons", [])]; fail = fail or judgment.get("fail_class") or "JUDGMENT_UNVERIFIABLE"
    return {"verdict": "PASS" if fail is None else "FAIL", "fail_class": fail, "reasons": reasons, "rc": rc,
            "objective": stead["objective"], "field": field, "residual": resid, "steadiness": stead,
            "floors": floors, "residual_scale": rscale, "floor_events": fevents,
            **({"judgment": judgment} if judgment is not None else {}),
            # プラトー到達は**診断**として常に載せる (受理条件ではない)
            "plateau": {"all_flat": not resid.get("converging") and not resid.get("rising"),
                        "falling": resid.get("converging", []), "rising": resid.get("rising", []),
                        "flat": resid.get("stalled", [])},
            "require_residual_pass": bool(require_residual_pass),
            "require_residual_plateau": bool(require_residual_plateau)}


def forge_rc_from_log(run_dir) -> int | None:
    """run_case.sh の stdout ログから最後の forge 終了コードを読む (rc を保存していない旧 run の再判定用)。"""
    p = Path(run_dir) / "run_case_stdout.log"
    if not p.exists():
        return None
    m = re.findall(r"\[run_case\] forge exit=(\d+)", p.read_text())
    return int(m[-1]) if m else None
