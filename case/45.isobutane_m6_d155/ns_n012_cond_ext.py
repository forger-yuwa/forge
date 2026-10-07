"""plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U4「N0 の凝縮の延長」(2026-10-07 登録、諮問
notes/reviews/2026-10-07-u4-v5prime-result-interpretation-diagnose.md): 凝縮の run を 1 回だけ延長して判定する。

  prep SRC_COND EXT   : SRC_COND (凝縮、res_18000) の入力を写し、restart_field.py で保存量をビット一致で継ぐ。solverConfig の差は
                        nStepOuter (18000 → 20000) だけ。凝縮は再初期化しない。記録 NS_N012.json (role cond_ext)
  run EXT             : ns_n012.py run と同じ経路で forge を回す (run_staged_ns、段なし)
  judge SRC_COND EXT  : 凝縮 4 量 (cond_series.py の量) を親と延長を連結した通算 34000〜38000 の 5 枚で check_quasisteady
                        (末尾 5 枚の引数、ns_n012_eval.quasisteady と同じ呼び方) → STEADY、残差は親と延長を連結した全列で
                        check_convergence → RISING・DIVERGED なし。両方なら K の合格。→ _band_ab/ns_n012_cond_ext_<EXT>.json

usage (AWS の case dir で): FORGE_BIN=... FORGE_CUDA_BLOCKSIZE=128 python3 ns_n012_cond_ext.py prep run_0168_ns_n012_N0_cond run_0180_ns_n012_N0_cond_ext
"""
import csv
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ns_n012 as NS  # noqa: E402
import ns_n012_eval as EV  # noqa: E402

PARENT_STEPS = NS.COND_STEPS          # 18000
EXT_STEPS = 20000
OUT = NS.COND_OUT                      # 1000
QS_COLS = ["onset_x_axis", "S_max", "exit_g_core", "exit_core_M"]
COND_COPY = NS.EXT_COPY + ("species_db.yaml",)


def prep(src: Path, run: Path) -> dict:
    NS.check_dry_env(False)
    if not NS.RUN_RE.match(run.name):
        raise SystemExit(f"run 名 {run.name!r} が run_NNNN_<slug> でない")
    if run.exists():
        raise SystemExit(f"{run} が既にある — 止める")
    srec = NS.jload(src / NS.RECORD)
    if srec.get("role") != "cond" or srec.get("dry"):
        raise SystemExit(f"src {src.name} は凝縮の run でない (role {srec.get('role')}) — 止める")
    rs = NS.res_files(src)
    if not rs or rs[-1].name != f"res_{PARENT_STEPS}.h5":
        raise SystemExit(f"src の最後の res が {rs[-1].name if rs else None} (res_{PARENT_STEPS}.h5 であること) — 止める")
    if (src / "RUN_RC").read_text().strip() != "0":
        raise SystemExit("src の RUN_RC が 0 でない — 止める")
    ys = NS.yaml_strict()
    ptext = (src / "solverConfig.yaml").read_text()
    P = ys.load(ptext)
    if NS.cfg_get(P, NS.NSTEP) != PARENT_STEPS or NS.cfg_get(P, NS.OUTINT) != OUT:
        raise SystemExit("src の solverConfig が凝縮の本段 (18000・1000 ごと) でない — 止める")
    ctext = ys.replace_scalars(ptext, {NS.NSTEP: str(EXT_STEPS)})
    dd = NS.MK.diff_paths(P, ys.load(ctext))
    if set(dd) != {NS.NSTEP}:
        raise SystemExit(f"延長の solverConfig の差が nStepOuter だけでない ({sorted(dd)}) — 止める")
    run.mkdir(parents=True)
    copied = []
    for f in COND_COPY:
        if (src / f).is_file():
            shutil.copy2(src / f, run / f)
            copied.append(f)
    for p in sorted(src.glob("resolved_species_*.yaml")):
        shutil.copy2(p, run / p.name)
        copied.append(p.name)
    (run / NS.RECORD).unlink(missing_ok=True)
    (run / "solverConfig.yaml").write_text(ctext)
    # restart_field は DST にあるデータセットにしか書かない。凝縮の run の入力 h5 にはモーメント (rog_*・roQ*_*) が無いので
    # (forge は無ければ 0 で始める)、src の保存量のうち DST に無いものを src と同じ型・形で先に作る (値は restart_field が書いて検査する)
    import h5py
    with h5py.File(rs[-1], "r") as h:
        want = sorted(k for k in h["VALUE"] if NS.CONS_RE.match(k))
        created = []
        with h5py.File(run / "nozzle.h5", "r+") as d:
            for k in want:
                if k not in d["VALUE"]:
                    d["VALUE"].create_dataset(k, data=h["VALUE"][k][()])
                    created.append(k)
    log = run / "restart_field.log"
    cmd = [sys.executable, str(NS.TOOLS / "restart_field.py"), str(rs[-1]), str(run / "nozzle.h5"), "--dst-run", str(run)]
    r = subprocess.run(cmd, capture_output=True, text=True, env=NS.runner()._ENV)
    out = r.stdout + r.stderr
    log.write_text(out)
    if r.returncode != 0 or "ビット一致" not in out:
        print(out[-3000:])
        raise SystemExit(f"restart_field がビット一致を確認していない (rc {r.returncode}) — 止める ({log})")
    # 移した量の数: 凝縮の run の /VALUE の保存量 (CONS_RE) の全部であること (化学種・モーメントを落とさない)
    tail = [l for l in out.splitlines() if "VERDICT" in l]
    import re
    m = re.search(r"(\d+) 量を移した", " ".join(tail))
    if not m or int(m.group(1)) != len(want):
        raise SystemExit(f"移した量の数 {m.group(1) if m else None} が src の保存量 {len(want)} ({want}) と違う — 止める")
    info = NS.jload(run / "prepare_info.json")
    info.update(stages={"stages": "none", "ramp": None, "ramp_steps": 1000}, extends=src.name,
                restart_from=f"{src.name}/{rs[-1].name} (restart_field)")
    NS.jdump(run / "prepare_info.json", info)
    rec = {**{k: srec.get(k) for k in ("plan", "condition", "md_offset", "md_offset_token", "problem", "problem_sha256", "k_f",
                                       "r_throat", "r_throat_requested", "out_interval")},
           "tool": "ns_n012_cond_ext.py prep", "created": NS.now(), "git_head": NS.git_head(), "dry": False, "role": "cond_ext",
           "parent": src.name, "parent_end": PARENT_STEPS, "ext_steps": EXT_STEPS, "stages": "none",
           "src": rs[-1].name, "src_sha256": NS.sha256_file(rs[-1]), "copied": copied, "conserved_moved": want,
           "datasets_created_in_dst_before_restart": created,
           "restart_field_tail": tail[-1:], "nozzle_sha256_after_prep": NS.sha256_file(run / "nozzle.h5")}
    NS.jdump(run / NS.RECORD, rec)
    print(f"[prep-cond-ext] {run.name} ← {src.name}/{rs[-1].name}: {EXT_STEPS} step、{len(want)} 量、{rec['restart_field_tail']}")
    return rec


def concat_rows(a: list, b: list, offset: int) -> list:
    out = [dict(r) for r in a]
    for r in b:
        r = dict(r)
        r["step"] = str(int(r["step"]) + offset)
        out.append(r)
    return out


def judge(src: Path, ext: Path) -> dict:
    out = {"item": "N0 の凝縮の延長の判定 (plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U4、2026-10-07 登録)",
           "parent": src.name, "ext": ext.name, "window_steps": list(range(PARENT_STEPS + EXT_STEPS - 4 * OUT, PARENT_STEPS + EXT_STEPS + 1, OUT))}
    # 凝縮 4 量: 親 (cond_series.csv、1000〜18000) + 延長 (1000〜20000 → +18000)
    pa, pe = EV.read_rows(src / "cond_series.csv"), EV.read_rows(ext / "cond_series.csv")
    EV.expect_steps(pa, OUT, PARENT_STEPS, src.name)
    EV.expect_steps(pe, OUT, EXT_STEPS, ext.name)
    rows = concat_rows(pa, pe, PARENT_STEPS)
    cat = ext / "cond_series_concat.csv"
    EV.write_rows(cat, rows)
    qs = EV.quasisteady(cat, QS_COLS, len(rows))
    win = {c: EV.window_stats(rows, c, out["window_steps"], "連結") for c in QS_COLS}
    out["quasisteady"] = qs
    out["window"] = {c: {k: win[c][k] for k in ("mean", "min", "max", "range", "finite")} for c in QS_COLS}
    qs_ok = all(qs[c]["verdict"] == "STEADY" for c in QS_COLS)
    # 残差: 親 (residual_history.csv) + 延長 (step を +18000) を連結した一時 dir で check_convergence (全列)
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        with open(src / "residual_history.csv") as f1, open(ext / "residual_history.csv") as f2, open(td / "residual_history.csv", "w", newline="") as fo:
            r1, r2 = csv.reader(f1), csv.reader(f2)
            h1, h2 = next(r1), next(r2)
            if h1 != h2:
                raise SystemExit("親と延長の残差履歴の列が違う — 止める")
            w = csv.writer(fo)
            w.writerow(h1)
            si = h1.index("step") if "step" in h1 else 0
            n1 = n2 = 0
            for row in r1:
                w.writerow(row); n1 += 1
            for row in r2:
                row[si] = str(int(row[si]) + PARENT_STEPS); w.writerow(row); n2 += 1
        q = subprocess.run([sys.executable, str(NS.TOOLS / "check_convergence.py"), str(td)], capture_output=True, text=True)
        text = q.stdout + q.stderr
    (ext / "CONVERGENCE_VERDICT_concat.txt").write_text(f"# 親 {src.name} ({n1} 行) + 延長 {ext.name} ({n2} 行、step +{PARENT_STEPS}) を連結\n" + text)
    head = next((l for l in text.splitlines() if l.startswith("===")), "")
    rising = [l.strip() for l in text.splitlines() if "RISING" in l or "DIVERGED" in l]
    out["convergence_concat"] = {"head": head, "rising_or_diverged": rising, "rows": [n1, n2]}
    res_ok = bool(head) and not rising and "DIVERGED" not in head
    out["checks"] = {"quasisteady_4_STEADY": qs_ok, "residual_no_rising": res_ok}
    out["verdict"] = ("合格 (K: 4 量 STEADY・残差 RISING なし" + ("、plateau は収束ではない)" if "NOT CONVERGED" in head else ")")) if (qs_ok and res_ok) \
        else "未達 — 再延長せず、未達のままユーザ判断 (登録)"
    dst = HERE / "_band_ab" / f"ns_n012_cond_ext_{ext.name}.json"
    NS.jdump(dst, out)
    print(json.dumps({"verdict": out["verdict"], "qs": {c: qs[c]["verdict"] for c in QS_COLS}, "conv": head, "rising": rising}, ensure_ascii=False, indent=1))
    return out


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "prep":
        prep(HERE / sys.argv[2], HERE / sys.argv[3])
    elif mode == "run":
        sys.exit(NS.run_one(HERE / sys.argv[2]))
    elif mode == "judge":
        judge(HERE / sys.argv[2], HERE / sys.argv[3])
    else:
        raise SystemExit(__doc__)
