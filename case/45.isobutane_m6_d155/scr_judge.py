"""§6.16 (2) のふるいの判定 (plan time_integration-line-implicit-speed)。各案について time_pairs_judge の結果 (scr_time_<案>.json) と
実効の設定の確認 (各 run の scr_modes.txt が OK、B0 を含む) を読み、組の判定が「速い」かつ組の差の中央値が B0 の中央値の −5 % 以下 → 「総時間へ進む」、
それ以外 → 「捨てる」。確認が欠けた・外れた案は「判定不能」(総時間へ進めない)。usage: python3 scr_judge.py → _band_ab/cold_pair/scr_judge.json"""
import json, statistics
from pathlib import Path
HERE = Path(__file__).resolve().parent; D = HERE / "_band_ab" / "cold_pair"
def modes_ok(v):
    for i in (1, 2, 3):
        f = HERE / f"run_037{2 + i}_scr_{v}_{i}" / "scr_modes.txt"
        if not f.exists() or f.read_text().strip() != "OK": return False
    return True
out = {}
b0ok = modes_ok("B0")
for v in ("M100", "M64", "M48", "S3", "S2"):
    p = D / f"scr_time_{v}.json"
    if not (b0ok and modes_ok(v) and p.exists()):
        out[v] = {"verdict": "判定不能 (設定の確認か組の判定が欠けた)"}; continue
    j = json.loads(p.read_text()); b0 = statistics.median(j["A_ms"]); rel = j["median_diff_ms"] / b0
    go = j["verdict"] == "速い" and rel <= -0.05
    out[v] = {"B0_median_ms": b0, "variant_ms": j["B_ms"], "median_diff_ms": j["median_diff_ms"], "rel": rel, "pair_verdict": j["verdict"],
              "verdict": "総時間へ進む" if go else "捨てる"}
(D / "scr_judge.json").write_text(json.dumps(out, indent=1, ensure_ascii=False)); print(json.dumps(out, indent=1, ensure_ascii=False))
