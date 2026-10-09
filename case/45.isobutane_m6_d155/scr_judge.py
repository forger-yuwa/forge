"""§6.16 (2) のふるいの判定 (plan time_integration-line-implicit-speed、codex plan-6 の反映)。
各巡 i の並び B0 v B0 v … B0 について、案 v の 1 step の時間を前後の B0 の平均と比べる: 差 d_i = v − (B0_前 + B0_後)/2、相対 d_i / B0。
雑音の尺度 = 全巡の隣り合う B0 どうしの差の最大値 (同じ設定の繰り返しの揺れ)。
判定: 3 巡の d がすべて負で中央値の大きさが雑音を超える → 速い、すべて正で同じ → 遅い、ほか → 判別不能 (再測定の対象)。
総時間へ進む = 速い かつ 相対の中央値 ≤ −5 % かつ (M 系は) 部分被覆の照合 run_0372_ml64_cmp/cmp_judge_LAYOUT2.json が PASS。
5 % は計算予算の優先の基準で、通らない案は「総時間は未評価」(効果なしとは判定しない)。どの run かの確認 (scr_check.json) が欠けた・ok でない → 判定不能。
usage: python3 scr_judge.py → _band_ab/cold_pair/scr_judge.json"""
import json, statistics
from pathlib import Path
HERE = Path(__file__).resolve().parent; D = HERE / "_band_ab" / "cold_pair"
VARS = ("M100", "M64", "M48", "S3", "S2")
def chk(run):
    p = HERE / run / "scr_check.json"
    if not p.exists(): return None
    j = json.loads(p.read_text()); return j["ms"] if j.get("ok") else None
out = {"rounds": {}, "variants": {}}
cmp = HERE / "run_0372_ml64_cmp" / "cmp_judge_LAYOUT2.json"
cmp_v = json.loads(cmp.read_text())["verdict"] if cmp.exists() else "記録なし"
out["cmp_partial_coverage"] = cmp_v
diffs = {v: [] for v in VARS}; noise = []; missing = []
for i in (1, 2, 3):
    runs = sorted(p.name for p in HERE.glob(f"run_037{2 + i}_scr_r{i}_*"))
    seq = [(r, r.rsplit("_", 1)[1], chk(r)) for r in runs]
    out["rounds"][i] = [(r, v, ms) for r, v, ms in seq]
    if len(seq) != 11 or [v for _, v, _ in seq][0::2] != ["B0"] * 6:
        missing.append(f"巡 {i}: 並びが B0 v … B0 の 11 本でない"); continue
    b0 = [ms for _, _, ms in seq[0::2]]
    noise += [abs(b - a) for a, b in zip(b0, b0[1:]) if a is not None and b is not None]
    for k in range(1, 11, 2):
        r, v, ms = seq[k]; a, b = seq[k - 1][2], seq[k + 1][2]
        if None in (ms, a, b): missing.append(f"{r} かその前後の B0 の確認が欠けた"); continue
        diffs[v].append({"run": r, "d_ms": ms - 0.5 * (a + b), "rel": (ms - 0.5 * (a + b)) / (0.5 * (a + b)), "b0_ms": 0.5 * (a + b), "v_ms": ms})
nz = max(noise) if noise else None
out["noise_ms"] = nz; out["missing"] = missing
for v in VARS:
    d = diffs[v]
    if len(d) != 3 or nz is None:
        out["variants"][v] = {"verdict": "判定不能 (3 巡そろわない)", "diffs": d}; continue
    med = statistics.median(x["d_ms"] for x in d); rel = statistics.median(x["rel"] for x in d)
    pv = "速い" if all(x["d_ms"] < 0 for x in d) and abs(med) > nz else "遅い" if all(x["d_ms"] > 0 for x in d) and abs(med) > nz else "判別不能 (再測定の対象)"
    if pv == "速い" and rel <= -0.05 and (not v.startswith("M") or cmp_v == "PASS"): verdict = "総時間へ進む"
    elif pv == "速い" and rel <= -0.05: verdict = f"進めない (部分被覆の照合が {cmp_v})"
    elif pv == "速い": verdict = "5 % 未満 (総時間は未評価)"
    else: verdict = pv
    out["variants"][v] = {"diffs": d, "median_d_ms": med, "median_rel": rel, "pair_verdict": pv, "verdict": verdict}
(D / "scr_judge.json").write_text(json.dumps(out, indent=1, ensure_ascii=False)); print(json.dumps(out["variants"], indent=1, ensure_ascii=False)); print("noise_ms", nz, "cmp", cmp_v, "missing", missing)
