# 本段の残差履歴 (residual_history.csv の outer_end 行) を 200 step ごとに間引き、全場の出力の時刻 (mtime) も書く
import csv, os, re, sys, json
run, out = sys.argv[1], sys.argv[2]
with open(os.path.join(run, "residual_history.csv")) as f:
    r = csv.reader(f); h = next(r)
    idx = [i for i, k in enumerate(h) if k.startswith("rms_")]
    rows = []
    for row in r:
        if len(row) < 3 or row[2] != "outer_end":
            continue
        s = int(row[0])
        if s % 200 == 0:
            rows.append([s] + [row[i] for i in idx])
with open(out, "w") as f:
    f.write(",".join(["step"] + [h[i] for i in idx]) + "\n")
    for row in rows:
        f.write(",".join(map(str, row)) + "\n")
t = {}
for fn in os.listdir(run):
    m = re.fullmatch(r"res_(\d+)\.h5", fn)
    if m:
        t[int(m.group(1))] = os.path.getmtime(os.path.join(run, fn))
ms = os.path.getmtime(os.path.join(run, "MAIN_STAGE_STARTED")) if os.path.exists(os.path.join(run, "MAIN_STAGE_STARTED")) else None
json.dump({"res_mtime": t, "main_started": ms, "now": __import__("time").time()}, open(out + ".times.json", "w"))
print(len(rows), rows[-1][0] if rows else None)
