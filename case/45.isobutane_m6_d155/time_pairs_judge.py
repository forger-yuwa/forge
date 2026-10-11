"""性能の組の判定 (plan time_integration-line-implicit-speed §6.2 (3) と同じ規則、codex 2026-10-10 m5 の反映):
各 run について (i) ログの "[line] factor 1 回目: モード X" と "[line] solve 1 回目: モード X" が期待のモード (A か B) で比較・診断なし、
(ii) step 101..999 がちょうど 1 回ずつで ms が有限かつ正、(iii) 投入前の GPU の他の計算プロセス (run の gpu_procs.txt) が 0 本、を確かめる。
外れた run があれば判定不能 (終了コード 2)。組ごとの差 (B − A) が 3 組とも負でその中央値の大きさが A の 3 本の幅を超える → 速い、3 組とも正で同じ条件 → 遅い、ほかは判別不能。
usage: python3 time_pairs_judge.py <out.json> <A のモード> <B のモード> A1 B1 A2 B2 A3 B3"""
import json, math, re, statistics, sys
from collections import Counter
from pathlib import Path
HERE = Path(__file__).resolve().parent
def check(run, mode):
    log = (HERE / run / "forge_run.log").read_text(errors="replace")
    for ph in ("factor", "solve"):
        m = re.search(rf"^\[line\] {ph} 1 回目: モード (\S+?)( \(比較あり\))?$", log, re.M)
        if not m or m.group(1) != mode or m.group(2):
            print(f"{run}: {ph} の実効のモードが {mode} (比較なし) でない ({m.group(0) if m else '行なし'}) — 判定不能"); sys.exit(2)
    g = HERE / run / "gpu_procs.txt"
    if not g.exists() or g.read_text().strip() != "0":
        print(f"{run}: 投入前の GPU の他の計算プロセスが 0 本でない ({g.read_text().strip() if g.exists() else '記録なし'}) — 判定不能"); sys.exit(2)
    rows = [(int(k), float(v)) for k, v in re.findall(r"^step\s+(\d+) \| ([0-9.eE+-]+|nan|inf) ms/step", log, re.M)]
    sel = [(k, v) for k, v in rows if 101 <= k <= 999]
    if Counter(k for k, _ in sel) != Counter(range(101, 1000)) or not all(math.isfinite(v) and v > 0 for _, v in sel):
        print(f"{run}: step 101..999 がちょうど 1 回ずつでないか ms が有限・正でない — 判定不能"); sys.exit(2)
    return sum(v for _, v in sel) / len(sel)
out_path = HERE / sys.argv[1]; ma, mb = sys.argv[2], sys.argv[3]; runs = sys.argv[4:]; assert len(runs) == 6
out_path.unlink(missing_ok=True)
A = [check(r, ma) for r in runs[0::2]]; Bv = [check(r, mb) for r in runs[1::2]]
d = [b - a for a, b in zip(A, Bv)]; spread = max(A) - min(A); med = statistics.median(d)
verdict = "速い" if all(x < 0 for x in d) and abs(med) > spread else "遅い" if all(x > 0 for x in d) and abs(med) > spread else "判別不能"
out = {"runs": runs, "modes": [ma, mb], "A_ms": A, "B_ms": Bv, "pair_diff_ms": d, "median_diff_ms": med, "A_spread_ms": spread, "verdict": verdict}
out_path.write_text(json.dumps(out, indent=1, ensure_ascii=False)); print(json.dumps(out, indent=1, ensure_ascii=False))
