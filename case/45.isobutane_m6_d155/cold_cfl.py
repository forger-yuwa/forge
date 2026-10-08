"""冷却壁の腕の擬似 CFL の試行と、質量の収支の監査の準備 (plan tooling-nozzle-isothermal-wall-chain §5.1 #27)。
cold_pair.py の prep_ext (延長) と同じ手順 (restart_field でビット一致・FP64 の型のまま、段なし) で、変えてよい設定を
nStepOuter・cfl・cfl_pseudo・outStepInterval と output.extraFields に限る (implicitRelax は変えない — codex 2026-10-08)。
走行中の run が cold_pair.py を使っているので、cold_pair.py は書き換えずにここで包む。

usage (AWS の case dir):
  python3 cold_cfl.py prep <src_run> <run> --steps N --cfl C [--out 5000] [--extra res_ro,volume] [--limiter-ref-from <run>]
    --limiter-ref-from: リミッタの基準値 (limiterRoRef・limiterPRef・limiterARef) を指定した run の forge_run.log の値に固定する。
      既定 (自動) では開始場から決まるので、restart した run は親と別の作用素になる (forge の警告; run_0182 → run_0183 で a_ref が 1.6 % 違った)。
  python3 cold_pair.py run <run>       (投入は既存の run_one; FORGE_DUMP_MASSFLUX は投入側の環境変数で渡す)
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_pair as CP  # noqa: E402
NS = CP.NS


def limiter_refs(run: Path) -> dict:
    """forge_run.log の最後の limiterRoRef・limiterPRef・limiterARef (forge が貼るよう促す行) を読む。"""
    import re
    ref = {}
    for line in (run / "forge_run.log").read_text(errors="replace").splitlines():
        m = re.match(r"^\[limiter\]\s+(limiterRoRef|limiterPRef|limiterARef):\s*([0-9.eE+-]+)\s*$", line)
        if m:
            ref[m.group(1)] = m.group(2)
    if set(ref) != {"limiterRoRef", "limiterPRef", "limiterARef"}:
        raise SystemExit(f"{run.name}/forge_run.log にリミッタの基準値の 3 行がそろっていない: {ref}")
    return ref


def prep(src: Path, run: Path, steps: int, cfl: float, out_int: int, extra: list[str], ref_from: Path | None = None) -> dict:
    NS.check_dry_env(False)
    binrec = CP.binary_record()
    if not NS.RUN_RE.match(run.name) or run.exists():
        raise SystemExit(f"{run} の名前が不正か既にある — 止める")
    srec = NS.jload(src / CP.RECORD)
    rs = NS.res_files(src)
    if not rs:
        raise SystemExit(f"{src} に res が無い")
    src_h5 = rs[-1]
    ys = NS.yaml_strict()
    ptext = (src / "solverConfig.yaml").read_text()
    ctext = ys.replace_scalars(ptext, {NS.NSTEP: str(int(steps)), NS.CFL: repr(float(cfl)), NS.CFLP: repr(float(cfl)),
                                       NS.OUTINT: str(int(out_int))})
    pcfg = ys.load(ptext)
    if extra:
        if "output" in pcfg:
            raise SystemExit("親の設定に output がある — extraFields の足し方を決めていないので止める")
        ctext = ctext.rstrip("\n") + "\noutput: {level: 1, extraFields: [" + ", ".join(extra) + "]}\n"
    refs = limiter_refs(ref_from) if ref_from is not None else {}
    if refs:
        if any(k in pcfg.get("space", {}) for k in refs):
            raise SystemExit("親の設定に既にリミッタの基準値がある — 止める")
        if ctext.count("space: {") != 1:
            raise SystemExit("space が 1 行のフロー形式でない — 基準値の足し方を決めていないので止める")
        ctext = ctext.replace("space: {", "space: {" + ", ".join(f"{k}: {v}" for k, v in refs.items()) + ", ")
    allowed = ({NS.NSTEP, NS.CFL, NS.CFLP, NS.OUTINT} | ({("output",)} if extra else set())
               | {("space", k) for k in refs})
    diff = set(NS.MK.diff_paths(pcfg, ys.load(ctext)))
    if not diff <= allowed:
        raise SystemExit(f"許していない設定の差がある: {sorted(diff - allowed)} — 止める")
    if ys.load(ctext)["time"]["deltaT"]["implicitRelax"] != pcfg["time"]["deltaT"]["implicitRelax"]:
        raise SystemExit("implicitRelax が変わった — 止める")
    run.mkdir(parents=True)
    for fn in NS.EXT_COPY + ("wall_repr.json", "bcondConfig.yaml", "species_meta.yaml"):
        if (src / fn).is_file():
            shutil.copy2(src / fn, run / fn)
    for p in sorted(src.glob("resolved_species_*.yaml")):
        shutil.copy2(p, run / p.name)
    (run / "solverConfig.yaml").write_text(ctext)
    cmd = [sys.executable, str(NS.TOOLS / "restart_field.py"), str(src_h5), str(run / "nozzle.h5"), "--dst-run", str(run), "--keep-src-dtype"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=NS.runner()._ENV)
    (run / "restart_field.log").write_text(r.stdout + r.stderr)
    if r.returncode != 0 or "ビット一致" not in (r.stdout + r.stderr):
        print((r.stdout + r.stderr)[-3000:])
        raise SystemExit(f"restart_field がビット一致を確認していない (rc {r.returncode}) — 止める")
    info = NS.jload(run / "prepare_info.json")
    info.update(stages={"stages": "none", "ramp": None, "ramp_steps": 1000}, extends=src.name, restart_from=f"{src.name}/{src_h5.name}")
    NS.jdump(run / "prepare_info.json", info)
    keep = ("plan", "kind", "problem", "problem_sha256", "delta_r_csv", "euler_ref", "implicit_relax", "mesh_checks",
            "geometry_vs_production", "wall_thermal", "bcond_wall")
    rec = {**{k: srec[k] for k in keep if k in srec},
           "tool": "cold_cfl.py prep", "plan_item": "§5.1 #27", "created": NS.now(), "git_head": NS.git_head(), "binary": binrec,
           "stages": "none", "parent": src.name, "parent_res": src_h5.name, "parent_res_sha256": NS.sha256_file(src_h5),
           "ext_steps": int(steps), "cfl_main": float(cfl), "cfl_parent": srec.get("cfl_main"), "out_interval": int(out_int),
           "extra_fields": extra, "limiter_ref_from": ref_from.name if ref_from is not None else None, "limiter_refs": refs,
           "config_diff": sorted("/".join(p) for p in diff),
           "restart_field_tail": (r.stdout + r.stderr).strip().splitlines()[-1:], "nozzle_sha256_after_prep": NS.sha256_file(run / "nozzle.h5")}
    NS.jdump(run / CP.RECORD, rec)
    print(f"[cold_cfl prep] {run.name} ← {src.name}/{src_h5.name}: cfl {cfl}・{steps} step・出力 {out_int} ごと・extra {extra}; "
          f"{rec['restart_field_tail']}")
    return rec


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("prep"); p.add_argument("src"); p.add_argument("run")
    p.add_argument("--steps", type=int, required=True); p.add_argument("--cfl", type=float, required=True)
    p.add_argument("--out", type=int, default=5000); p.add_argument("--extra", default="")
    p.add_argument("--limiter-ref-from", default=None)
    a = ap.parse_args()
    prep(HERE / a.src, HERE / a.run, a.steps, a.cfl, a.out, [s for s in a.extra.split(",") if s],
         HERE / a.limiter_ref_from if a.limiter_ref_from else None)
