"""plan tooling-nozzle-core-grid §4.13 の起動の A/B の判定 (事前登録どおり)。usage: python3 cg_startab_judge.py <A_run> <B_run>
- 見る量: 全保存量の残差 rms_* (rms_dq_* は除く)、非有限の数、最後の場の ρ・T・P の最小値。
- 残差の爆発: 各列の末尾 200 step の中央値が step 100〜300 の最大値の 100 倍を超える (基準値が 0 の列は判定不能)。
- 判定: 両腕が通過 → 本段から直接始める / A だけ失敗 → 起動の手順を登録し直す / 両腕が失敗 → 長い run に進まない。
結果は標準出力と _band_ab/core_grid/startab.json。"""
import csv, json, math, sys
from pathlib import Path
import h5py, numpy as np
HERE = Path(__file__).resolve().parent
def arm(run):
    run = Path(run); rows = {}
    with open(run / "residual_history.csv") as fh:
        rd = csv.DictReader(fh)
        cols = [c for c in rd.fieldnames if c.startswith("rms_") and not c.startswith("rms_dq_")]
        for r in rd:
            try: st = int(float(r["step"]))
            except (KeyError, ValueError): continue
            rows[st] = r                                     # 同じ step の行は最後を残す (内反復の行)
    steps = sorted(rows); out = {"run": run.name, "last_step": steps[-1] if steps else None, "cols": {}}
    nonfin = 0; fails = []; undecided = []
    for c in cols:
        v = np.array([float(rows[s][c]) for s in steps]); st = np.array(steps)
        nf = int(np.count_nonzero(~np.isfinite(v))); nonfin += nf
        base_w = (st >= 100) & (st <= 300); tail_w = st > st[-1] - 200
        base = float(np.nanmax(v[base_w])) if base_w.any() else float("nan"); tail = float(np.nanmedian(v[tail_w])) if tail_w.any() else float("nan")
        ratio = tail / base if base > 0 else None
        out["cols"][c] = {"nonfinite": nf, "base_max_100_300": base, "tail_median_200": tail, "ratio": ratio}
        if ratio is None: undecided.append(c)
        elif ratio > 100: fails.append(c)
    res = sorted(run.glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))
    with h5py.File(res[-1], "r") as h:
        V = {k: np.asarray(h["VALUE/" + k][:], np.float64) for k in ("ro", "T", "P")}
        nf_field = int(sum(np.count_nonzero(~np.isfinite(np.asarray(h["VALUE/" + k][:], np.float64))) for k in h["VALUE"] if h["VALUE/" + k].ndim == 1))
    out.update(res=res[-1].name, field_nonfinite=nf_field, min={k: float(np.nanmin(v)) for k, v in V.items()}, residual_nonfinite=nonfin,
               exploded=fails, undecidable_cols=undecided)
    out["pass"] = bool(out["last_step"] is not None and out["last_step"] >= 1999 and nonfin == 0 and nf_field == 0 and not fails
                       and all(m > 0 for m in out["min"].values()))
    return out
A, B = arm(sys.argv[1]), arm(sys.argv[2])
verdict = ("両腕が通過 → 本段から直接始める" if A["pass"] and B["pass"] else
           "A だけ失敗・B は通過 → 起動の手順を登録し直す" if (not A["pass"]) and B["pass"] else
           "両腕が失敗 → 長い run に進まない (上位に諮る)" if not (A["pass"] or B["pass"]) else "A は通過・B は失敗 (想定外 → 上位に諮る)")
res = {"A_cfl4": A, "B_cfl05": B, "VERDICT": verdict}
(HERE / "_band_ab" / "core_grid").mkdir(parents=True, exist_ok=True)
(HERE / "_band_ab" / "core_grid" / "startab.json").write_text(json.dumps(res, ensure_ascii=False, indent=1))
for k, a in (("A (cfl_pseudo 4)", A), ("B (cfl_pseudo 0.5)", B)):
    print(f"{k} {a['run']}: 最後の step {a['last_step']}、残差の非有限 {a['residual_nonfinite']}、場の非有限 {a['field_nonfinite']}、最小 {a['min']}、爆発 {a['exploded']}、判定不能 {a['undecidable_cols']} → {'通過' if a['pass'] else '失敗'}")
print("VERDICT:", verdict)
