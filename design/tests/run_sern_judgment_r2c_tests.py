#!/usr/bin/env python3
"""R2c の台本 (case/46.sern_design/diag/judgment_r2c.sh・judgment_r2c.py) の分岐の試験 (plan tooling-sern-te-wake-grid §5.1 R2c・§6 #11)。
生産 YAML (lswx10) の格子を粗くした問題 YAML で、偽の起点 run (prepare だけ、場は格子の初期場) に対して台本の手順 1〜4 の分岐を確かめる。
変換器が無ければ省略。--gpu では起点 run を実際に段階起動 (粗い格子・短い段) し、台本を最後 (judge・post) まで通す (GPU・forge が要る)。

合格 (測る前に固定):
  (1) 通る入力: prepare が 0 で終わり、R2C_PREPARE.json に格子署名・双対の一致、config の差 (差なし / nStepOuter・出力系だけ)、
      化学種の状態 (属性なし → 許可が要る / 印付き → 要らない)、起点の場の複写 (sha256 一致) が残る
  (2) 既存ディレクトリ: シェルの台本も Python も 2 で止まり、既存の run を変えない
  (3) 格子署名の不一致・双対幾何のハッシュの不一致・節点数の不一致: 3 で止まる (R2C_STOP.txt)
  (4) 本段の config の差 (cfl)・bcondConfig.yaml の差・作動点の違い: 4 で止まり、差を R2C_CONFIG_DIFF.txt に残す
  (5) 起点の場の化学種の記録の破損 (species_hash はあるが記録が無い): 5 で止まる。起点 run の記録の欠け・空き容量の不足: 6
  (6) judge の環境の検査: FORGE_CUDA_BLOCKSIZE が 128 でない・FORGE_ALLOW_UNVERIFIED_SPECIES が prepare の判断と合わない → 6 で止まり、
      時間積分に進まない (forge を起動しない)
  (7) --gpu: 台本を通しで回すと JUDGMENT.json・R2C_JUDGE.json (環境の実際の値)・R2C_SUMMARY.txt・判定ツールの出力が run に残る。
      段階起動の出力 (印付きの場) を起点にすると許可なしで継続でき、段の移行は生産 config の必須量で検査される

  /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_sern_judgment_r2c_tests.py [--gpu] [--out DIR] [--keep]
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))

from forge_design.evaluate import runner_sern as R2  # noqa: E402
from forge_design.evaluate import runner_sern3d as R3  # noqa: E402

SH = ROOT / "case" / "46.sern_design" / "diag" / "judgment_r2c.sh"
PYT = ROOT / "case" / "46.sern_design" / "diag" / "judgment_r2c.py"
PROD = ROOT / "case" / "46.sern_design" / "problem_3d_prod_3op_wallres_lswx10.yaml"
COARSE3D = {"ni_up": 6, "ni_noz": 30, "ni_plume": 40, "nj_top": 21, "nj_bot": 15, "nz_in": 4, "nz_out": 3, "nj_vside": 5,
            "first_wall_frac": 2.0e-3, "first_z_frac": 0.02, "first_wake_frac": 0.004}
COARSE2 = {"nj_ext_top": 9, "first_top_frac": 0.004}
fails = 0


def check(name, ok, detail=""):
    global fails
    fails += not ok
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail and not ok else ""), flush=True)


ap = argparse.ArgumentParser()
ap.add_argument("--gpu", action="store_true", help="起点 run を段階起動し、台本を judge・post まで通す (GPU・forge が要る)")
ap.add_argument("--out", default=None)
ap.add_argument("--keep", action="store_true")
a = ap.parse_args()
if not R2.CONVERTER.exists():
    print(f"SKIP: 変換器 {R2.CONVERTER} が無い (prepare を回せない)")
    print("VERDICT: PASS (省略)")
    sys.exit(0)
S = Path(a.out) if a.out else Path(tempfile.mkdtemp(prefix="sern_r2c_"))
S.mkdir(parents=True, exist_ok=True)
y = yaml.safe_load(PROD.read_text())
y["name"] += "_r2c_test"
y["mesh3d"].update(COARSE3D); y["mesh"].update(COARSE2)
y["evaluate"].update({"nStepOuter": 20, "outStepInterval": 2, "floor_events_window_steps": 40, "judgment_steps": 80})
y["opt"].update({"warm_lam_steps": 20, "soft_steps": 20, "mid_steps": 20})
prob = S / "coarse.yaml"
prob.write_text(yaml.safe_dump(y, allow_unicode=True, sort_keys=False))
SRC = S / "src_run"
if not SRC.exists():
    R3.prepare(prob, SRC, op="m6_on")
    shutil.copyfile(SRC / "sern.h5", SRC / "res_4000.h5")        # 偽の起点の場 (格子の初期場。分岐の試験には中身は要らない)
ENV = {k: v for k, v in os.environ.items() if k != "FORGE_ALLOW_UNVERIFIED_SPECIES"}
COMMON = ["--case-dir", str(S), "--problem", str(prob), "--min-free-gb", "0.1"]


def run_py(stage, new, *extra, env=None):
    r = subprocess.run([sys.executable, str(PYT), stage, new, *COMMON, *extra], capture_output=True, text=True, env=env or ENV)
    return r.returncode, r.stdout + r.stderr


def variant(name, edit=None, files=("prepare_info.json", "solverConfig_main.yaml", "bcondConfig.yaml", "species_meta.yaml", "res_4000.h5")):
    """起点 run の複製 (記録だけ書き換える)。"""
    d = S / name
    if d.exists():
        shutil.rmtree(d)
    d.mkdir()
    for f in files:
        if (SRC / f).exists():
            shutil.copyfile(SRC / f, d / f)
    if edit:
        edit(d)
    return d


def edit_json(key_path, value):
    def _e(d):
        info = json.loads((d / "prepare_info.json").read_text())
        x = info
        for k in key_path[:-1]:
            x = x[k]
        x[key_path[-1]] = value
        (d / "prepare_info.json").write_text(json.dumps(info))
    return _e


# (1) 通る入力 (同じ config・属性なし → 許可が要る)
rc, log = run_py("prepare", "new_ok", "--src-run", str(SRC))
pj = json.loads((S / "new_ok" / "R2C_PREPARE.json").read_text()) if (S / "new_ok" / "R2C_PREPARE.json").exists() else {}
check("(1) 通る入力: prepare 0・格子署名と双対が一致・config の差なし・属性なしで許可が要る・起点の場の複写が一致・中間の sern.msh は消す",
      rc == 0 and pj.get("ok") and all(v["new"] == v["src"] for v in pj["mesh_identity"].values())
      and pj["config_diff"] == {"allowed": [], "disallowed": [], "files": {"bcondConfig.yaml": True, "species_meta.yaml": True}}
      and pj["species"]["state"] == "none" and pj["allow_unverified_species"] is True
      and pj["src_res_sha256"] == pj["copied_res_sha256"] and (S / "new_ok" / "res_4000.h5").exists()
      and not (S / "new_ok" / "sern.msh").exists(), log[-1500:])
# (1) nStepOuter・出力間隔だけ違う起点 run は通る (差を記録)
v = variant("src_nsteps", lambda d: (d / "solverConfig_main.yaml").write_text(
    (d / "solverConfig_main.yaml").read_text().replace("nStepOuter: 20", "nStepOuter: 4000").replace("outStepInterval: 2", "outStepInterval: 500")))
rc, log = run_py("prepare", "new_nsteps", "--src-run", str(v))
pj = json.loads((S / "new_nsteps" / "R2C_PREPARE.json").read_text()) if rc == 0 else {}
check("(1) nStepOuter・outStepInterval だけ違う: 通り、許す差として記録",
      rc == 0 and len(pj["config_diff"]["allowed"]) == 4 and not pj["config_diff"]["disallowed"]
      and "nStepOuter: 4000" in (S / "new_nsteps" / "R2C_CONFIG_DIFF.txt").read_text(), log[-800:])
# (1) 印付きの場 (species_hash と species_input_unverified=1) は許可が要らない
v = variant("src_marked")
with h5py.File(v / "res_4000.h5", "r+") as f:
    f.attrs["species_hash"] = "a" * 64; f.attrs["species_record_sha256"] = "b" * 64; f.attrs["species_record_file"] = "x.yaml"
    f.attrs.create("species_input_unverified", 1, dtype="int32")
rc, log = run_py("prepare", "new_marked", "--src-run", str(v))
pj = json.loads((S / "new_marked" / "R2C_PREPARE.json").read_text()) if rc == 0 else {}
check("(1) 印付きの場: 通り、許可は要らない (ハッシュの照合は restart_by_index とソルバに任せる)",
      rc == 0 and pj["species"]["state"] == "unverified" and pj["allow_unverified_species"] is False, log[-800:])

# (2) 既存ディレクトリ
before = sorted(p.name for p in (S / "new_ok").iterdir())
mt = (S / "new_ok" / "R2C_PREPARE.json").stat().st_mtime_ns
rc, log = run_py("prepare", "new_ok", "--src-run", str(SRC))
r = subprocess.run(["bash", str(SH), "new_ok", *COMMON, "--src-run", str(SRC)], capture_output=True, text=True,
                   env={**ENV, "PYTHON": sys.executable})
check("(2) 既存ディレクトリ: Python・シェルの台本とも 2 で止まり、既存の run を変えない",
      rc == 2 and r.returncode == 2 and "既にある" in log and "prepare で停止 (rc=2)" in r.stdout
      and sorted(p.name for p in (S / "new_ok").iterdir()) == before and (S / "new_ok" / "R2C_PREPARE.json").stat().st_mtime_ns == mt,
      (rc, r.returncode, r.stdout[-300:], r.stderr[-300:]))

# (3) 格子の不一致
for name, edit, label in (("src_sig", edit_json(["mesh_provenance", "mesh_signature"], "0" * 64), "格子署名"),
                          ("src_dual", edit_json(["mesh_provenance", "dual_hash"], "0" * 64), "双対幾何のハッシュ")):
    v = variant(name, edit)
    rc, log = run_py("prepare", "new_" + name, "--src-run", str(v))
    stopf = S / ("new_" + name) / "R2C_STOP.txt"
    check(f"(3) {label}の不一致: 3 で止まり R2C_STOP.txt に理由 (起点の場は置かない)", rc == 3 and stopf.exists() and label != ""
          and not (S / ("new_" + name) / "res_4000.h5").exists(), log[-600:])
v = variant("src_nodes")
with h5py.File(v / "res_4000.h5", "r+") as f:
    ro = f["VALUE/ro"][:]; del f["VALUE/ro"]; f.create_dataset("VALUE/ro", data=ro[:-1])
rc, log = run_py("prepare", "new_nodes", "--src-run", str(v))
check("(3) 起点の場の節点数の不一致: 3 で止まる", rc == 3 and "節点数" in log, log[-600:])

# (4) config・BC・作動点の差
v = variant("src_cfl", lambda d: (d / "solverConfig_main.yaml").write_text(
    (d / "solverConfig_main.yaml").read_text().replace("cfl: 0.25, cfl_pseudo: 0.25", "cfl: 0.5, cfl_pseudo: 0.5")))
rc, log = run_py("prepare", "new_cfl", "--src-run", str(v))
diff = (S / "new_cfl" / "R2C_CONFIG_DIFF.txt").read_text() if (S / "new_cfl" / "R2C_CONFIG_DIFF.txt").exists() else ""
check("(4) 本段の config の差 (cfl): 4 で止まり、差を R2C_CONFIG_DIFF.txt に残す", rc == 4 and "cfl: 0.5" in diff and "許さない差 2 行" in diff,
      (rc, log[-600:]))
v = variant("src_bc", lambda d: (d / "bcondConfig.yaml").write_text((d / "bcondConfig.yaml").read_text() + "# changed\n"))
rc, log = run_py("prepare", "new_bc", "--src-run", str(v))
check("(4) bcondConfig.yaml の差: 4 で止まる", rc == 4 and "bcondConfig.yaml が違う" in log, log[-600:])
v = variant("src_op", edit_json(["operating_point", "name"], "m10_on"))
rc, log = run_py("prepare", "new_op", "--src-run", str(v))
check("(4) 作動点の違い: 4 で止まる", rc == 4 and "作動点" in log, log[-600:])

# (5) 化学種の記録の破損・入力の欠け
v = variant("src_broken")
with h5py.File(v / "res_4000.h5", "r+") as f:
    f.attrs["species_hash"] = "a" * 64; f.attrs["species_record_sha256"] = "b" * 64; f.attrs["species_record_file"] = "x.yaml"
    f.attrs.create("species_input_unverified", 0, dtype="int32")
rc, log = run_py("prepare", "new_broken", "--src-run", str(v))
check("(5) 化学種の記録の破損 (species_hash はあるが記録が無い): 5 で止まる (起点の場は置かない)",
      rc == 5 and not (S / "new_broken" / "res_4000.h5").exists(), log[-600:])
v = variant("src_noprov", edit_json(["mesh_provenance"], None))
rc, log = run_py("prepare", "new_noprov", "--src-run", str(v))
check("(5) 起点 run に格子署名の記録が無い: 6 で止まる", rc == 6 and "記録が無い" in log, log[-600:])
rc, log = run_py("prepare", "new_nores", "--src-run", str(SRC), "--src-res", "res_9999.h5")
rc2, log2 = run_py("prepare", "new_disk", "--src-run", str(SRC), "--min-free-gb", "1e9")
check("(5) 起点の場が無い・空き容量の不足: 6 で止まり、run を作らない",
      rc == 6 and rc2 == 6 and not (S / "new_nores").exists() and not (S / "new_disk").exists(), (log[-300:], log2[-300:]))

# (6) judge の環境の検査 (forge を起動しない)
for label, env in (("FORGE_CUDA_BLOCKSIZE が無い", {**ENV, "FORGE_ALLOW_UNVERIFIED_SPECIES": "1"}),
                   ("FORGE_CUDA_BLOCKSIZE が 512", {**ENV, "FORGE_CUDA_BLOCKSIZE": "512", "FORGE_ALLOW_UNVERIFIED_SPECIES": "1"}),
                   ("許可が要るのに付いていない", {**ENV, "FORGE_CUDA_BLOCKSIZE": "128"})):
    env = {k: v_ for k, v_ in env.items() if not (label.startswith("FORGE_CUDA_BLOCKSIZE が無い") and k == "FORGE_CUDA_BLOCKSIZE")}
    rc, log = run_py("judge", "new_ok", env=env)
    check(f"(6) judge: {label} → 6 で止まり forge を起動しない", rc == 6 and not (S / "new_ok" / "runner_launches.jsonl").exists()
          and not (S / "new_ok" / R2.JUDGMENT_FILE).exists(), log[-400:])
rc, log = run_py("judge", "new_marked", env={**ENV, "FORGE_CUDA_BLOCKSIZE": "128", "FORGE_ALLOW_UNVERIFIED_SPECIES": "1"})
check("(6) judge: 許可が要らないのに付いている → 6 で止まる (その実行だけの許可を余計に付けない)", rc == 6, log[-400:])
rc, log = run_py("judge", "new_cfl", env={**ENV, "FORGE_CUDA_BLOCKSIZE": "128"})
check("(6) judge: prepare で止めた run は進めない (6)", rc == 6 and "prepare が通っていない" in log, log[-400:])

# (7) GPU: 起点 run を実際に段階起動して台本を通す
if a.gpu:
    G1 = S / "src_gpu"
    if not G1.exists():
        envg = {**ENV, "FORGE_CUDA_BLOCKSIZE": "128", "FORGE_ALLOW_UNVERIFIED_SPECIES": "1"}   # 起点を作る段階起動 (台本の外、検証の手順と同じ)
        code = ("import sys; sys.path.insert(0, %r)\n"
                "from forge_design.evaluate import runner_sern3d as R3\n"
                "R3.prepare(%r, %r, op='m6_on')\n"
                "rc = R3.run_staged(%r, 'full', 20, soft_cfl=0.5, soft_conv=0, warm_lam_steps=20, warm_lam_cfl=0.2, mid_steps=20)\n"
                "print('STAGED rc', rc)\n") % (str(HERE.parent), str(prob), str(G1), str(G1))
        rs = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=envg, timeout=3600)
        (S / "src_gpu.staged.log").write_text(rs.stdout + rs.stderr)
        print(rs.stdout[-300:])
    res = [p.name for _, p in R2._volume_snapshots(G1)]          # step の昇順 (最後 = 本段の最終場)
    if not res:
        check("(7) 起点 run の段階起動が場を残した", False, (S / "src_gpu.staged.log").read_text()[-800:])
    else:
        r = subprocess.run(["bash", str(SH), "new_gpu", *COMMON, "--src-run", str(G1), "--src-res", res[-1], "--poll-s", "1"],
                           capture_output=True, text=True, env={**ENV, "PYTHON": sys.executable}, timeout=3600)
        (S / "new_gpu.r2c.log").write_text(r.stdout + r.stderr)
        rd = S / "new_gpu"
        j = json.loads((rd / R2.JUDGMENT_FILE).read_text()) if (rd / R2.JUDGMENT_FILE).exists() else {}
        jj = json.loads((rd / "R2C_JUDGE.json").read_text()) if (rd / "R2C_JUDGE.json").exists() else {}
        print(r.stdout[-1500:])
        check("(7) 台本を通しで回すと JUDGMENT.json・R2C_JUDGE.json・R2C_SUMMARY.txt・判定ツールの出力が残る",
              j.get("verdict") in ("PASS", "FAIL", "INDETERMINATE") and (rd / "R2C_SUMMARY.txt").exists()
              and (rd / "R2C_FLOOR_EVENTS.txt").exists() and (rd / "R2C_CONVERGENCE_judgment1.txt").exists()
              and (rd / "R2C_QUASISTEADY_forces.txt").exists(), r.stdout[-600:] + r.stderr[-600:])
        en = (j.get("env") or {}).get("named") or {}
        pj = json.loads((rd / "R2C_PREPARE.json").read_text())
        # 段階起動 (許可つき) の出力は印付き (species_hash + species_input_unverified=1) なので、judge に許可は要らない
        check("(7) 起点は段階起動の印付きの場: prepare は許可を要らないと判断し、judge は許可なしで継続できた (継承の照合を通った)",
              pj["species"]["state"] == "unverified" and pj["allow_unverified_species"] is False
              and j["attempts"][0].get("rc") is not None, pj["species"])
        check("(7) 環境の実際の値: FORGE_CUDA_BLOCKSIZE 128、FORGE_ALLOW_UNVERIFIED_SPECIES は prepare の判断どおり (ここでは未設定 = null)、"
              "R2C_JUDGE.json・runner_launches.jsonl にも同じ値",
              en.get("FORGE_CUDA_BLOCKSIZE") == "128"
              and en.get("FORGE_ALLOW_UNVERIFIED_SPECIES") == ("1" if pj["allow_unverified_species"] else None)
              and jj.get("env", {}).get("named") == en
              and all(json.loads(x)["env"]["named"] == en for x in (rd / R2.RUNNER_LAUNCHES).read_text().splitlines() if x.strip()), en)
        check("(7) 段の移行の検査は生産 config の必須量 (流体 5・roK・roOmega・roY0・roY1) で通った",
              j.get("required_state") == ["ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega", "roY0", "roY1"]
              and j["attempts"][0]["restart_from"] == res[-1], j.get("required_state"))
        print(f"run (GPU): {rd}")
if not a.keep and not a.out:
    shutil.rmtree(S, ignore_errors=True)
print("VERDICT:", "PASS" if fails == 0 else f"FAIL ({fails})")
sys.exit(1 if fails else 0)
