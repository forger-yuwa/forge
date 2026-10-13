"""`nonpos_ref` 検査の有無だけを変える抽出の A/B (plan tooling-nozzle-deltastar-band-root-zone の result 段の諮問 2026-10-13 の推奨、forge 0 step)。

腕 A = 探索域の q_ref の有限・正の検査 (`nonpos_ref`) だけを無効にした `band_local_deficit` (探索域の限定・帯選択・平滑化は今のまま)、腕 B = 今の実装。
同じ NS/Euler の対を `extract_and_merge(.., band_select="edge")` に通し、最終の `delta_r_equiv.csv`・`delta_r_next.csv` のバイト同一、
最終の断面の `nonpos_ref` の件数、使用域 [0, x_e] の hard 不合格 (held) の数、呼び出しの件数 (梯子・感度を含む) を比べる。
入力の run は読むだけ (出力は OUT_DIR)。

usage (design/.venv-opt):  python nonpos_ref_ab.py NS_RUN EULER_RUN X_E OUT_DIR ARM(A|B)
"""
import inspect
import json
import sys
from pathlib import Path

import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C.parents[1] / "design"))
import forge_design.metrics.deltastar as DS  # noqa: E402
from forge_design.feedback.deltastar_loop import extract_and_merge, read_delta_r_next  # noqa: E402

ns, eu, x_e, out, arm = Path(sys.argv[1]), Path(sys.argv[2]), float(sys.argv[3]), Path(sys.argv[4]), sys.argv[5]
src = inspect.getsource(DS.band_local_deficit)
chk = "    if not (np.all(np.isfinite(qz)) and np.all(qz > 0.0)):\n"
assert src.count(chk) == 1, "nonpos_ref の検査の行が見つからない"
cnt = {"calls": 0, "nonpos": 0}
inst = ("    _CNT['calls'] += 1\n"
        "    if not (np.all(np.isfinite(qz)) and np.all(qz > 0.0)):\n"
        "        _CNT['nonpos'] += 1\n")
new = src.replace(chk, inst + ("    if False:\n" if arm == "A" else chk))
g = dict(DS.__dict__)
g["_CNT"] = cnt
exec(compile(new, "patched", "exec"), g)
DS.band_local_deficit = g["band_local_deficit"]
g["band_local_deficit"] = g["band_local_deficit"]
od = out / arm
od.mkdir(parents=True, exist_ok=True)
extract_and_merge(ns, eu, band_select="edge", out_dir=od)
d = json.loads((od / "delta_r_equiv_diag.json").read_text())
x = np.array(d["x"], float); rs = list(d["reason"]); hard = np.array(d["hard_ok"], bool)
use = (x >= 0.0) & (x <= x_e)
nx = read_delta_r_next(od / "delta_r_next.csv")
held = np.asarray(nx.get("held", np.zeros(len(nx["x_rt"]))), float)
xn = np.asarray(nx["x_rt"], float)
res = {"ns_run": str(ns), "euler_run": str(eu), "arm": arm, "x_e": x_e,
       "calls_reaching_check": cnt["calls"], "calls_nonpos_ref": cnt["nonpos"],
       "final_rows": int(len(x)), "final_nonpos_ref": int(sum("nonpos_ref" in r for r in rs)),
       "final_nonpos_ref_in_use": int(sum(("nonpos_ref" in r) and u for r, u in zip(rs, use))),
       "hard_fail_in_use": int((~hard & use).sum()),
       "held_in_use": int(((held > 0) & (xn >= 0.0) & (xn <= x_e)).sum())}
(od / "nonpos_ref_ab.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
print(json.dumps(res, ensure_ascii=False))
