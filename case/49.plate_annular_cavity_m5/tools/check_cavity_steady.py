#!/usr/bin/env python3
"""キャビティの結論量 (Q, q'', h_ref, T_mid, T_floor, 侵入深さ …) の準定常判定。

AGENTS.md の「準定常確認 (必須)」の実体化。汎用の `check_quasisteady.py` は衝撃位置や
CL/CD 用の量しか持たないので、**判定ロジック (`classify`) だけを借りて**、対象量は
`cavity_eval.py --series` が吐く `cavity_series.csv` (全 `res_*.h5` スナップショットの時系列)
から取る。

usage:
  python3 tools/cavity_eval.py <run> --series          # 先に時系列を作る
  python3 tools/check_cavity_steady.py <run> [...]     # 判定 (複数 run 可)
    [--tail 0.4] [--drift 0.02] [--osc 0.05] [--min-snaps 4] [--quantity Q_cav,h_ref]

終了コード: 全量 STEADY/OSCILLATING なら 0、DRIFTING か UNSETTLED があれば 1。
"""
import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
# **`classify` を直呼びしない**。非有限を黙って落とすため、[1,1,1,1,1,NaN] が STEADY になる。
# 正本の `classify_series` (非有限を NONFINITE にする) を使う (2026-09-19 codex Major 4)。
from check_quasisteady import classify_series as classify, SEV      # noqa: E402

# **必須列**: 1 つでも欠けていたら / 非有限が混じっていたら FAIL にする
# (2026-09-19 codex Major: 旧実装は欠けた列を無言で飛ばし、NaN を除去する classify を呼ぶため
#  h_ref が全点 NaN でも他が STEADY なら PASS になっていた)
DEFAULT = ["q_wall_sum", "q_outer", "q_cylside", "q_floor", "dT_mouth", "dT_mid", "dT_floor",
           "zpen_25", "mdot_in", "mdot_out", "mdot_imbalance", "h_ref"]

# **非一様壁温では判定に使えない列** (2026-09-20)。`h_ref = <q''/(T0_mid - Tw)>_A` は
# 局所係数の面積平均なので、壁ごとに壁温が違うと**すきま中央のガス温度が壁温の中間**に来て
# ΔT ≈ 0 の場所が必ず生じ、`q''/ΔT` が発散する (mixB 実測: 平均 1.71e4 W/m2K・変動 252 %、
# 一方で壁入熱は 5-6 桁一致の完全な定常)。§4.7.3 の「(A) 局所 T0 基準は非一様壁温では
# 定義できない」と同じ理由。**壁温が一様でない run では除外する**。
NONUNIFORM_EXCLUDE = ["h_ref"]
LABEL = {"q_wall_sum": "3壁 総入熱 Q [W]", "q_outer": "外筒壁 Q [W]", "q_cylside": "円柱側面 Q [W]",
         "q_floor": "底面 Q [W]", "dT_mouth": "開口 dT [K]", "dT_mid": "中央 dT [K]",
         "dT_floor": "底 dT [K]", "zpen_25": "侵入深さ(25K) [m]", "mdot_in": "開口流入 [kg/s]",
         "mdot_out": "開口流出 [kg/s]", "mdot_imbalance": "正味/片道", "h_ref": "h_ref [W/m2K]"}
# **絶対許容** (plan §4.7 / case.json の eval.tol_*)。平均で正規化する相対 drift では
# ゼロ近傍の量が判定できないので、固定尺度でも見る。
# 許容は **max(絶対下限, 値に対する割合)**。固定の絶対値だけだと作動点が変わったときに
# 非現実的な厳しさになる (2026-09-19: M5 で dT_mouth ~70 K を想定して 2 K としたが、
# M9 では dT_mouth が 449 K になり 2 K = 0.45 % を要求してしまい、振れ 1.1 % で不合格になった)。
# **plan §6.4 の絶対許容をそのまま使う** (2026-09-21 codex result M8)。
# 従来は `max(2 K, 平均の 2 %)` としており、ΔT が数百 K のケースで許容が 9 K まで
# 緩んでいた (末尾窓で 6 K 動く単調増加でも `STEADY` を返す反例を codex が再現)。
# plan は 2 K・0.5 mm を要求しているので、**相対緩和を外して plan と一致させる**。
# 緩めたいときは plan §6.4 を先に直し、根拠と対象量を書くこと。
ABS_TOL = {"dT_mouth": ("K", 2.0, 0.0), "dT_mid": ("K", 2.0, 0.0), "dT_floor": ("K", 2.0, 0.0),
           "zpen_25": ("m", 0.5e-3, 0.0), "mdot_imbalance": ("-", 0.005, 0.0)}
# 壁ごとの入熱は **総入熱に対する割合**で許容する (絶対値が小さい壁ほど相対 drift は当てにならない)
Q_REL_TOL = 0.005      # 3 壁合計の 0.5 %


def is_nonuniform_wall(run):
    """`conditions.json` に `wall_T_by_group` があり、壁温が実際に割れているか。"""
    import json as _json
    f = Path(run) / "conditions.json"
    if not f.exists():
        return False
    d = _json.loads(f.read_text())
    by = d.get("wall_T_by_group") or {}
    if not by:
        return False
    vals = set(round(float(v), 6) for v in by.values()) | {round(float(d["wall_T"]), 6)}
    return len(vals) > 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--tail", type=float, default=0.4)
    ap.add_argument("--drift", type=float, default=0.02)
    ap.add_argument("--osc", type=float, default=0.05)
    ap.add_argument("--min-snaps", type=int, default=4)
    ap.add_argument("--quantity", default=None, help="カンマ区切り (既定は主要量すべて)")
    a = ap.parse_args()
    want = a.quantity.split(",") if a.quantity else DEFAULT

    worst = 0
    for run in a.runs:
        cols = list(want)
        if not a.quantity and is_nonuniform_wall(run):
            drop = [c for c in NONUNIFORM_EXCLUDE if c in cols]
            for c in drop:
                cols.remove(c)
            if drop:
                print("  [非一様壁温] 判定から除外: %s  "
                      "(壁温が割れていると ΔT≈0 の点が生じ q''/ΔT が発散する。§4.7.3)"
                      % ", ".join(drop))
        csvp = Path(run) / "cavity_series.csv"
        print("\n=== %s ===" % run)
        if not csvp.exists():
            print("  cavity_series.csv が無い -> `cavity_eval.py %s --series` を先に回す" % run)
            worst = max(worst, SEV["TRANSIENT-UNSETTLED"])
            continue
        with open(csvp) as f:
            rows = list(csv.DictReader(f))
        if not rows:
            print("  空の時系列")
            worst = max(worst, SEV["TRANSIENT-UNSETTLED"])
            continue
        steps = [float(r["step"]) for r in rows]
        print("  スナップショット %d 点 (step %d..%d),  tail=%.0f%%  drift<%.3g  osc<%.3g"
              % (len(steps), steps[0], steps[-1], 100 * a.tail, a.drift, a.osc))
        qsum_scale = 0.0
        if "q_wall_sum" in rows[0]:
            qv = [float(r["q_wall_sum"]) for r in rows if r["q_wall_sum"] not in ("", None)]
            qv = [x for x in qv if x == x]
            if qv:
                n0 = max(2, int(len(qv) * a.tail))
                qsum_scale = abs(sum(qv[-n0:]) / n0)
        missing = [k for k in cols if k not in rows[0]]
        if missing:
            print("  **必須列が無い**: %s  -> FAIL (系列の作り直しが要る)" % ", ".join(missing))
            worst = max(worst, SEV["NONFINITE"])
        for k in cols:
            if k not in rows[0]:
                continue
            raw = [r[k] for r in rows]
            vals = [float(x) if x not in ("", None) else float("nan") for x in raw]
            nbad = sum(1 for x in vals if not (x == x and abs(x) != float("inf")))
            if nbad:
                print("  %-22s **非有限 %d/%d 点** -> FAIL" % (LABEL.get(k, k), nbad, len(vals)))
                worst = max(worst, SEV["NONFINITE"])
                continue
            verd, detail, ma = classify(steps, vals, a.tail, a.drift, a.osc, a.min_snaps)
            tag = LABEL.get(k, k)
            ext = ""
            if ma and verd == "OSCILLATING":
                ext = "   -> 平均 %.4g +/- %.3g で報告" % ma
            tol = None
            if k in ABS_TOL:
                unit, floor_, rel_ = ABS_TOL[k]
                tol = max(floor_, rel_ * abs(ma[0]) if ma else floor_)
            elif k in ("q_wall_sum", "q_outer", "q_cylside", "q_floor") and qsum_scale:
                unit, tol = "W", Q_REL_TOL * qsum_scale
            if tol is not None and ma:
                n = max(2, int(len(vals) * a.tail))
                swing = max(vals[-n:]) - min(vals[-n:])
                within = swing <= tol
                ext += "   [絶対 %s %.3g / 許容 %.3g %s]" % (
                    "OK" if within else "**超過**", swing, tol, unit)
                if not within:
                    worst = max(worst, SEV["DRIFTING"])
                elif verd in ("DRIFTING", "TRANSIENT-UNSETTLED"):
                    # **ゼロ近傍の量を相対 drift で落とさない**。絶対許容は「工学的に意味のある
                    # 変化量」なので、そこに収まっていれば相対がいくら大きくても定常扱いにする
                    # (実例: 偏心時の底面 Q は 1.1e-4 W = 総入熱の 0.0004 % で相対 drift 14.6 %、
                    #  絶対の振れは 2e-5 W で許容 0.13 W の 1/6000。2026-09-19)。
                    verd = "STEADY"
                    ext += "  ← 相対は大きいが**絶対許容内**なので定常扱い"
                    worst = max(worst, SEV["STEADY"])
                    print("  %-22s %-20s %s%s" % (tag, verd, detail, ext))
                    continue
            worst = max(worst, SEV.get(verd, 4))
            print("  %-22s %-20s %s%s" % (tag, verd, detail, ext))
    name = [k for k, v in SEV.items() if v == worst]
    print("\nVERDICT: %s" % ("STEADY (全量)" if worst == 0 else (name[0] if name else "?")))
    return 0 if worst <= SEV["OSCILLATING"] else 1


if __name__ == "__main__":
    sys.exit(main())
