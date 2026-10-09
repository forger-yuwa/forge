"""§6.16 (2) の 1 run の確認 (plan time_integration-line-implicit-speed、codex plan-6 の反映)。出力を消す前に呼ぶ。
確かめること: (a) 投入前の GPU 照会が成功して計算プロセス 0 本 (gpu_pre.txt = "ok 0")、計測中の見張り (gpu_watch.txt) で照会の失敗 0 回・最大本数 ≤ 1 (自分だけ)、
(b) ログの "[line] factor/solve 1 回目: モード LAYOUT2" (比較なし)、被覆と最長が期待どおり、設定の nStepInner が期待どおり、
(c) 表示の step 101..999 がちょうど 1 回ずつで ms が有限・正、(d) residual_history.csv の全 rms_* 列が全行で有限、(e) 最終場 res_1000.h5 の VALUE が全量有限。
結果は <run>/scr_check.json ({"ok": bool, "ms": 101..999 の平均, "why": [...]})。usage: python3 scr_check.py <run> <被覆 CV 数> <最長> <nStepInner>"""
import csv, json, math, re, sys
from collections import Counter
from pathlib import Path
import h5py, numpy as np
run, cov, mlen, inner = Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
why = []
pre = (run / "gpu_pre.txt").read_text().split() if (run / "gpu_pre.txt").exists() else []
if pre != ["ok", "0"]: why.append(f"投入前の GPU 照会が成功・0 本でない ({pre})")
w = (run / "gpu_watch.txt").read_text().split() if (run / "gpu_watch.txt").exists() else []
if len(w) != 3 or w[0] != "done" or int(w[1]) != 0 or int(w[2]) > 1: why.append(f"計測中の GPU の見張りが失敗 0・最大 1 本以下でない ({w})")
log = (run / "forge_run.log").read_text(errors="replace")
for ph in ("factor", "solve"):
    m = re.search(rf"^\[line\] {ph} 1 回目: モード (\S+?)( \(比較あり\))?$", log, re.M)
    if not m or m.group(1) != "LAYOUT2" or m.group(2): why.append(f"{ph} のモードが LAYOUT2 (比較なし) でない")
m = re.search(r"^\[lineImplicit\] lines=(\d+)  covered CVs=(\d+)/(\d+) .*maxLen=(\d+)", log, re.M)
if not m or (int(m.group(1)), int(m.group(2)), int(m.group(4))) != (4719, cov, mlen): why.append(f"被覆・最長が期待 (4719, {cov}, {mlen}) と違う ({m.group(0) if m else '行なし'})")
if not re.search(rf"nStepInner: *{inner}([,}} ]|$)", (run / "solverConfig.yaml").read_text(), re.M): why.append(f"nStepInner が {inner} でない")
rows = [(int(k), float(v)) for k, v in re.findall(r"^step\s+(\d+) \| ([0-9.eE+-]+|nan|inf) ms/step", log, re.M)]
sel = [(k, v) for k, v in rows if 101 <= k <= 999]
ms = None
if Counter(k for k, _ in sel) != Counter(range(101, 1000)) or not all(math.isfinite(v) and v > 0 for _, v in sel): why.append("step 101..999 がそろわないか ms が有限・正でない")
else: ms = sum(v for _, v in sel) / len(sel)
with open(run / "residual_history.csv") as f:
    rd = csv.DictReader(f); cols = [c for c in rd.fieldnames if c.startswith("rms_")]; n = bad = 0
    for r in rd:
        n += 1; bad += sum(1 for c in cols if not math.isfinite(float(r[c])))
if n == 0 or bad: why.append(f"残差に非有限 {bad} 個 (行 {n})")
with h5py.File(run / "res_1000.h5", "r") as h:
    nf = sum(int(np.count_nonzero(~np.isfinite(np.asarray(d[:])))) for d in h["VALUE"].values() if np.asarray(d[:]).dtype.kind == "f")
if nf: why.append(f"最終場の非有限 {nf} 個")
out = {"ok": not why, "ms": ms, "why": why}
(run / "scr_check.json").write_text(json.dumps(out, ensure_ascii=False)); print(run.name, json.dumps(out, ensure_ascii=False))
