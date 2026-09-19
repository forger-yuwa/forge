"""⑤ SERN 評価の受理ゲート (plan §5.1 R1, codex C1 採用 2026-09-09)。

1 作動点 run を「サロゲート学習・Pareto に入れてよいか」を判定する。ゲートは互いに独立で、全て必須:
  (1) forge の終了コード rc == 0 (発散 run の採用は撤回。旧 §4.13-3 は廃止)
  (2) 保存場の有限値・正値 (最終 res_*.h5 の ro/roU/roe/P/T が有限、ro,P,T > 0、res_nan_*.h5 が無い)
  (3) 全残差 (check_convergence.analyze): NaN/Inf 無し、末尾で rising な列が無い。PASS (3 桁低下) は
      `require_residual_pass` のときだけ必須 (本ケースは残差が 1–2.5 桁でプラトーする性格 — その事実は台帳に残す)
  (4) **実際に最適化する量** (RANS は C_T_with_shear、Euler は C_T) と C_T / C_L / C_M の頭打ち
      (check_quasisteady.classify_series による STEADY。NONFINITE / DRIFTING / TRANSIENT-UNSETTLED / OSCILLATING は不合格)
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

FIELD_VARS = ("ro", "roUx", "roUy", "roUz", "roe", "P", "T", "roK", "roOmega")
POSITIVE_VARS = ("ro", "P", "T")
FORCE_KEYS = ("C_T", "C_L", "C_M")


def _volume_res(run_dir):
    return sorted(Path(run_dir).glob("res_[0-9]*.h5"), key=lambda f: int("".join(c for c in f.stem if c.isdigit())))


def field_health(run_dir) -> dict:
    """最終体積出力の保存量・原始量が有限で正値かを見る。res_nan_*.h5 (発散ダンプ) があれば不合格。"""
    run_dir = Path(run_dir)
    nan_files = sorted(f.name for f in run_dir.glob("res_nan*.h5"))
    res = _volume_res(run_dir)
    out = {"ok": True, "file": None, "nonfinite": {}, "nonpositive": {}, "min": {}, "nan_files": nan_files, "reasons": []}
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


def objective_key(hist) -> str:
    """実際に最適化する量: RANS (壁摩擦出力あり) は C_T_with_shear、Euler は C_T。"""
    return "C_T_with_shear" if hist and "C_T_with_shear" in hist[-1] else "C_T"


def steadiness_gate(hist, obj: str | None = None) -> dict:
    obj = obj or objective_key(hist)
    keys = [obj] + [k for k in FORCE_KEYS if k != obj]
    steps = [h["step"] for h in hist]
    series = {k: steadiness([h.get(k, float("nan")) for h in hist], steps=steps) for k in keys}
    unsteady = [k for k, v in series.items() if v["verdict"] != "STEADY"]
    return {"ok": not unsteady, "objective": obj, "series": series, "unsteady": unsteady,
            "reasons": [f"{k}: {series[k]['verdict']} ({series[k]['detail']})" for k in unsteady]}



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
    ただし温度は `dependentVariables_d.cu` の反転クランプ `DEPVAR_TMIN` 50 K が実効下限。"""
    out = {"P": 1.0, "ro": 1.0e-4, "T": 50.0}
    f = Path(run_dir) / "solverConfig.yaml"
    if f.exists():
        txt = f.read_text()
        for key, name in (("pMin", "P"), ("roMin", "ro"), ("tMin", "T")):
            m = re.search(rf"{key}\s*:\s*([-\d.eE+]+)", txt)
            if m:
                v = float(m.group(1))
                out[name] = max(v, out[name]) if name == "T" else v
    return out


def floor_gate(run_dir, p_min: float | None = None, tol: float = 1.0e-6) -> dict:
    """最終保存場で EOS 床・乱流下限に張り付いたノードを数える。1 個でも在れば NG。
    床は run の `solverConfig.yaml` から読む。**読み取れない・保存場が無いときは判定不能で不合格**
    (旧実装は例外を握り潰して ok=True を返していた: codex M6)。
    `k <= 0` は見ない — node 低 Re 壁は `sstNodeWallKPin` (既定 1) が正当に k=0 をピンする。"""
    import glob
    import os
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
            # k の床判定は壁ピン (正当) と区別できないので**当面外す** (R-f)
    except Exception as e:                                 # **読めないときは判定不能で不合格**
        return {"ok": False, "counts": {}, "reasons": [f"判定不能: {type(e).__name__}: {e}"]}
    # `|sum...|` は許容超過 **と 判定不能 (-1)** の両方を不合格にする (codex plan-2 M3)
    bad = {k: v for k, v in counts.items()
           if (k.startswith("|sum") and (v > 1.0e-4 or v < 0.0)) or (not k.startswith("|sum") and v != 0)}
    return {"ok": not bad, "counts": counts, "file": os.path.basename(fs[-1]),
            "reasons": [] if not bad else ["床/下限に張り付き: " + ", ".join(f"{k} {v} ノード" for k, v in bad.items())]}


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


def evaluate_gates(run_dir, hist, rc, require_residual_pass: bool = False, obj: str | None = None,
                   p_min: float | None = None, require_residual_plateau: bool = False) -> dict:
    """全ゲートを評価して verdict / fail_class を返す。fail_class は数値失敗の種別:
    DIVERGED (rc≠0 / 発散ダンプ / 非有限・非正の場 / 残差 NaN), RESIDUAL_RISING, NOT_CONVERGED (require 時のみ),
    NO_FORCES (壁出力が無い), UNSTEADY (目的量・力係数が頭打ちしていない)。物理的 INFEASIBLE はここでは出さない。"""
    field = field_health(run_dir)
    resid = residual_health(run_dir)
    floors = floor_gate(run_dir, p_min)
    rscale = residual_scale_gate(run_dir)
    stead = steadiness_gate(hist, obj) if hist else {"ok": False, "objective": obj or objective_key(hist), "series": {}, "unsteady": [],
                                                    "reasons": ["no force history (no wall output)"]}
    reasons = []; fail = None
    if rc is None:
        reasons.append("forge rc unknown"); fail = fail or "DIVERGED"
    elif rc != 0:
        reasons.append(f"forge rc={rc}"); fail = fail or "DIVERGED"
    if not field["ok"]:
        reasons += field["reasons"]; fail = fail or "DIVERGED"
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
    if not rscale["ok"]:
        reasons += rscale["reasons"]; fail = fail or "RESIDUAL_UNBALANCED"
    if not hist:
        reasons += stead["reasons"]; fail = fail or "NO_FORCES"
    elif not stead["ok"]:
        reasons += stead["reasons"]; fail = fail or "UNSTEADY"
    return {"verdict": "PASS" if fail is None else "FAIL", "fail_class": fail, "reasons": reasons, "rc": rc,
            "objective": stead["objective"], "field": field, "residual": resid, "steadiness": stead,
            "floors": floors, "residual_scale": rscale,
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
