#!/usr/bin/env python3
"""S3 性能 (plan gradient-scalar-lsq-unification §6 S3): lsq/gg の step 時間比 ≤ 1.05。

同一 GPU・同一 BLOCKSIZE (128/128)・同一バイナリで、S2 双子の最終場から 2500 step を回し、
ログの 1 step ごとの壁時計 (`step N | X ms/step`) の step 501–2499 の平均を 1 本の測定値にする
(先頭 500 step はウォームアップ、最終 step は出力を含むので除く)。gg と lsq を交互に 3 反復し、
各 kind の中央値の比を出す。各 run の前に他の forge プロセスが無いことを確かめ、あれば待つ (GPU 占有の測定)。

  python3 s3_perf.py --case case/48.flat_plate_cooled_m4 --gg-src RUN --lsq-src RUN --res res_24000.h5 --out S3_case48.txt
"""
import argparse
import os
import re
import shutil
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
FBIN = os.path.expanduser("~/sglsq/forge_f99f236d")


def gpu_pids():
    """GPU を使っている計算プロセスの PID (`nvidia-smi --query-compute-apps`)。プロセス名に依らない
    (codex result-1 M4: 測定バイナリ `forge_f99f236d` は `pgrep -x forge` に数えられず、競合を見逃していた)。"""
    r = subprocess.run(["nvidia-smi", "--query-compute-apps=pid", "--format=csv,noheader"], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("nvidia-smi が使えない (GPU 占有を確認できない)")
    return {int(x) for x in r.stdout.split() if x.strip().isdigit()}


def solver_pids(root_pid):
    """root_pid (run_case.sh の bash) の子孫のうち、実行ファイルが FBIN のプロセス。"""
    target = os.path.realpath(FBIN)
    out, stack = set(), [root_pid]
    while stack:
        p = stack.pop()
        try:
            kids = open(f"/proc/{p}/task/{p}/children").read().split()
        except OSError:
            kids = []
        for k in map(int, kids):
            stack.append(k)
            try:
                if os.path.realpath(f"/proc/{k}/exe") == target:
                    out.add(k)
            except OSError:
                pass
    return out


def run_one(case, src, res, grad, dst):
    if os.path.exists(dst):
        shutil.rmtree(dst)
    subprocess.run([sys.executable, os.path.join(HERE, "s2_setup.py"), "--start", src, "--res", os.path.join(src, res),
                    "--dst", dst, "--grad", grad, "--nstep", "2500", "--out-int", "2500"], check=True,
                   capture_output=True, text=True)
    env = dict(os.environ, FORGE_BIN=FBIN, FORGE_CUDA_BLOCKSIZE="128", FORGE_CUDA_BLOCKSIZE_SMALL="128")
    for _ in range(5):
        while gpu_pids():
            time.sleep(20)
        # 走行中に他の forge が現れたら測定を捨てて取り直す (2026-09-26: 開始前だけ見ていて case/48 の 2 反復が重なった)
        pr = subprocess.Popen(["bash", os.path.join(REPO, "solver_density_cuda", "tools", "run_case.sh"), dst], env=env,
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # 競合判定 (codex result-1 M4): GPU 計算プロセスの PID を、run_case.sh の子孫で実行ファイルが FBIN の
        # solver PID と照合する。solver 以外の GPU プロセスが 1 回でも見えたら汚染、solver が GPU 上に 1 回も
        # 観測されない・異常終了も無効。0.5 s ごとのポーリングなので、それより短い競合は検出できない。
        clean, observed, mine = True, False, set()
        while pr.poll() is None:
            time.sleep(0.5)
            mine |= solver_pids(pr.pid)
            now = gpu_pids()
            if now & mine:
                observed = True
            if now - mine:
                clean = False
                print(f"{dst}: solver 以外の GPU プロセス {sorted(now - mine)}", flush=True)
        if clean and observed and pr.returncode == 0:
            break
        if not observed:
            print(f"{dst}: solver を GPU 上で観測できなかった (rc={pr.returncode})", flush=True)
        print(f"{dst}: 走行中に他の forge が重なったので取り直す", flush=True)
    else:
        raise SystemExit(f"{dst}: 5 回とも他の forge と重なった")
    ms = []
    for line in open(os.path.join(dst, "forge_run.log"), errors="replace"):
        m = re.match(r"step\s+(\d+) \|\s+([\d.]+) ms/step", line)
        if m and 501 <= int(m.group(1)) <= 2499:
            ms.append(float(m.group(2)))
    if len(ms) < 1900:
        raise SystemExit(f"{dst}: 測定 step が足りない ({len(ms)})")
    return sum(ms) / len(ms), len(ms)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True)
    ap.add_argument("--gg-src", required=True)
    ap.add_argument("--lsq-src", required=True)
    ap.add_argument("--res", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.chdir(REPO)
    got = {"gg": [], "lsq": []}
    lines = [f"# S3 性能 ({a.case})", f"binary: {FBIN}", "block 128/128、2500 step、step 501–2499 の 1 step 壁時計の平均、gg/lsq 交互 3 反復"]
    for rep in "abc":
        for g in (("gg", "lsq") if rep != "b" else ("lsq", "gg")):
            src = a.gg_src if g == "gg" else a.lsq_src
            dst = os.path.join(a.case, f"_s3_{g}_{rep}")
            v, n = run_one(a.case, src, a.res, g, dst)
            got[g].append(v)
            lines.append(f"- {g}_{rep}: {v:.3f} ms/step (n={n}) {dst}")
    mg, ml = statistics.median(got["gg"]), statistics.median(got["lsq"])
    lines.append(f"中央値: gg {mg:.3f} / lsq {ml:.3f} ms/step → 比 lsq/gg = {ml / mg:.4f} (上限 1.05)")
    lines.append(f"VERDICT S3 ({a.case}): {'PASS' if ml / mg <= 1.05 else 'FAIL'}")
    txt = "\n".join(lines) + "\n"
    open(a.out, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
