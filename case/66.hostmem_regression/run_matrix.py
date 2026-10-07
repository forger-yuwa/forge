#!/usr/bin/env python3
"""case/66 回帰ハーネス: 構成 × ビルド × 反復の run を作って AWS で順に回す (plan architecture-solver-host-memory §5.1 #3)。

AWS の ~/forge-b4/case/66.hostmem_regression/ に置いて使う (ローカルでは forge を回さない)。

    # 1) 種の適用 (同一メッシュ restart の入力だけ。restart_field.py を掛ける。1 回だけ)
    python3 run_matrix.py seed --bin ~/bin-hostmem/forge_9c9f623c [入力名 ...] [--force-species]
    # 2) 投入 (構成ごとに run_NNNN_<構成>_<ビルド>_r<k> を作り、1 本ずつ順に回すワーカーを裏で起動して即戻る)
    python3 run_matrix.py launch --build base --bin ~/bin-hostmem/forge_9c9f623c \
        --convert-bin ~/bin-hostmem/convert_9c9f623c --reps 1 2 3 c36node c09ckpt100 ...
    # 3) dual-time 再開の checkpoint を置く (base r1 の ckpt100 が終わってから。全ビルド・全反復で同じファイルを使う)
    python3 run_matrix.py set-ckpt c09restart100 run_0004_c09ckpt100_base_r1
    # 4) 状態 / 待ち合わせ / README の表
    python3 run_matrix.py status [--build base]
    python3 run_matrix.py table

- バイナリは --bin で明示し、run_case.sh に FORGE_BIN で渡す。RUN_PROVENANCE.txt の forge_bin と sha256 を投入後に照合する (verify)。
- FORGE_CUDA_BLOCKSIZE=128 を常に付ける (skill forge-aws-run)。stdin は /dev/null、ワーカーは nohup 相当 (新しいセッション)。
- 他セッションの forge には触らない。生存確認・メモリ採取は自分が起動したプロセスの子を /proc でたどる (pgrep -f を使わない)。
- 各 run の終了後に NaN/Inf を検査して NANCHECK.txt を残す (plan §6 (c))。
- launch-dir は既に入力をそろえた run ディレクトリ (case/46 の SERN 等) を同じワーカーで回す (--memwatch でホスト RSS/HWM と GPU を 1 s ごとに採る)。
"""
import argparse
import csv
import fcntl
import glob
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
RUN_CASE = os.path.join(ROOT, "solver_density_cuda", "tools", "run_case.sh")
RESTART_FIELD = os.path.join(ROOT, "solver_density_cuda", "tools", "restart_field.py")
REGISTRY = os.path.join(HERE, "registry.tsv")
REG_COLS = ["run", "cfg", "input", "build", "rep", "kind", "bin", "bin_sha256", "env", "created"]
sys.path.insert(0, HERE)
import matrix_spec as ms  # noqa: E402


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_registry():
    if not os.path.exists(REGISTRY):
        return []
    with open(REGISTRY) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def run_state(rundir):
    p = os.path.join(rundir, ".state")
    return open(p).read().strip() if os.path.exists(p) else "?"


# ---------------------------------------------------------------- 種と checkpoint
def cmd_seed(a):
    names = a.names or [n for n, s in ms.INPUTS.items() if "seed" in s]
    for n in names:
        spec = ms.INPUTS[n]
        d = os.path.join(HERE, "inputs", n)
        if "seed" not in spec:
            print(f"[seed] {n}: 種なし (skip)")
            continue
        if os.path.exists(os.path.join(d, "SEEDED.txt")) and not a.redo:
            print(f"[seed] {n}: 適用済み (SEEDED.txt)。やり直すなら --redo (元の h5 を入れ直してから)")
            continue
        pre = []
        if spec["seed"].get("precreate"):
            # 宛先 (変換直後の h5) に無い保存量は restart_field が写さないので、先に 0 で作っておく (skill forge-aws-run §3)
            import h5py
            import numpy as np
            with h5py.File(os.path.join(d, spec["seed"]["dst"]), "r+") as f:
                ref = f["VALUE/ro"]
                for nm in spec["seed"]["precreate"]:
                    if "VALUE/" + nm not in f:
                        f.create_dataset("VALUE/" + nm, data=np.zeros(ref.shape, dtype=ref.dtype))
                        pre.append(nm)
        cmd = [sys.executable, RESTART_FIELD, os.path.join(d, "seed_src.h5"), os.path.join(d, spec["seed"]["dst"]),
               "--dst-run", d, "--forge", os.path.expanduser(a.bin)]
        if a.force_species or spec["seed"].get("force_species"):
            cmd.append("--force-species")
        env = dict(os.environ, FORGE_CUDA_BLOCKSIZE="128")
        r = subprocess.run(cmd, capture_output=True, text=True, env=env, stdin=subprocess.DEVNULL)
        out = r.stdout + r.stderr
        print(f"[seed] {n}: rc={r.returncode}\n" + "\n".join("    " + x for x in out.strip().splitlines()[-6:]))
        if r.returncode == 0 and "VERDICT: OK" in out:
            with open(os.path.join(d, "SEEDED.txt"), "w") as f:
                f.write(f"date: {time.strftime('%Y-%m-%dT%H:%M:%S%z')}\ncmd : {' '.join(cmd)}\n"
                        f"precreate (宛先に 0 で先に作った量): {pre}\n\n{out}")


def cmd_set_ckpt(a):
    spec = ms.INPUTS[a.input]
    c = spec["ckpt"]
    src_run = os.path.join(HERE, a.run)
    src = os.path.join(src_run, c["file"])
    dst = os.path.join(HERE, "inputs", a.input, c["dst"])
    if os.path.exists(dst) and not a.redo:
        sys.exit(f"[set-ckpt] {dst} は既にある (全ビルドで同じ checkpoint を使うため上書きしない。やり直すなら --redo)")
    if c["cfg"] not in os.path.basename(src_run):
        sys.exit(f"[set-ckpt] {a.run} は {c['cfg']} の run ではない")
    shutil.copy2(src, dst)
    with open(os.path.join(HERE, "inputs", a.input, "CKPT_FROM.txt"), "w") as f:
        f.write(f"{c['dst']} <- {a.run}/{c['file']} (sha256 {sha256(dst)})\ndate: {time.strftime('%Y-%m-%dT%H:%M:%S%z')}\n")
    print(f"[set-ckpt] {dst} <- {src}")


# ---------------------------------------------------------------- 投入
def forge_alias(binp, label):
    """バイナリを `forge` という名前のシンボリックリンク経由で起動する (.bin/<label>/forge)。

    AWS の idle 自動停止 (tools/cloud/idle_autostop.sh) と他セッションの確認は `pgrep -x forge` で forge の有無を見る。
    `forge_9c9f623c` のような名前のまま起動すると「forge 無し」と判定され、短い run が続く間 (GPU 使用率の瞬間値が 0) に
    インスタンスが止まる (2026-10-07 に本ハーネスの base 投入中に停止した。原因の推定は README)。comm は exec したパスの basename なので、
    リンク名を forge にすれば pgrep -x forge に一致する。来歴 (RUN_PROVENANCE の sha256・verify) は実体で照合する。
    """
    real = os.path.realpath(os.path.expanduser(binp))
    if os.path.basename(binp) == "forge":
        return os.path.expanduser(binp)
    d = os.path.join(HERE, ".bin", label)
    os.makedirs(d, exist_ok=True)
    link = os.path.join(d, "forge")
    if os.path.islink(link):
        if os.path.realpath(link) != real:
            sys.exit(f"[launch] {link} は別のバイナリ ({os.path.realpath(link)}) を指している。ビルド名 (--build) を変えること")
    else:
        os.symlink(real, link)
    return link


def next_number():
    nums = [int(m.group(1)) for r in read_registry() for m in [re.match(r"run_(\d+)_", r["run"])] if m]
    nums += [int(m.group(1)) for p in glob.glob(os.path.join(HERE, "run_*"))
             for m in [re.match(r"run_(\d+)_", os.path.basename(p))] if m]
    return max(nums, default=0) + 1


def populate(rundir, input_name):
    d = os.path.join(HERE, "inputs", input_name)
    spec = ms.INPUTS[input_name]
    if "seed" in spec and not os.path.exists(os.path.join(d, "SEEDED.txt")):
        sys.exit(f"[launch] 入力 {input_name} に種が未適用 (run_matrix.py seed を先に)")
    files = [x for x in open(os.path.join(d, "FILES")).read().split() if x]
    files += [os.path.basename(p) for p in glob.glob(os.path.join(d, "resolved_species_*.yaml"))]
    for fn in files:
        p = os.path.join(d, fn)
        if not os.path.exists(p):
            sys.exit(f"[launch] 入力 {input_name} に {fn} が無い (checkpoint なら run_matrix.py set-ckpt を先に)")
    os.makedirs(rundir)
    for fn in files:
        shutil.copy2(os.path.join(d, fn), os.path.join(rundir, fn))
    with open(os.path.join(rundir, "INPUT_FILES"), "w") as f:   # NaN 検査・比較で「入力」として除外するファイル
        f.write("\n".join(files) + "\n")
    for fn in ("SOURCE.txt", "SEEDED.txt", "CKPT_FROM.txt"):
        if os.path.exists(os.path.join(d, fn)):
            shutil.copy2(os.path.join(d, fn), os.path.join(rundir, "INPUT_" + fn))


def spawn_worker(jobs, tag):
    jf = os.path.join(HERE, f".jobs_{tag}.json")
    with open(jf, "w") as f:
        json.dump(jobs, f, indent=1)
    log = open(os.path.join(HERE, f"worker_{tag}.log"), "a")
    p = subprocess.Popen([sys.executable, os.path.abspath(__file__), "worker", jf], cwd=HERE, stdin=subprocess.DEVNULL,
                         stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
    print(f"[launch] worker pid {p.pid}: {len(jobs)} 本を順に回す (ログ worker_{tag}.log)")


def cmd_launch(a):
    binp = forge_alias(a.bin, a.build)
    cbin = os.path.expanduser(a.convert_bin) if a.convert_bin else None
    for c in a.cfgs:
        if c not in ms.CONFIGS:
            sys.exit(f"[launch] 未知の構成 {c}")
        if ms.CONFIGS[c].get("kind", "forge") == "convert" and not cbin:
            sys.exit(f"[launch] {c} は変換器の構成: --convert-bin が要る")
    if a.dry_run:
        n = next_number()
        for rep in a.reps:
            for c in a.cfgs:
                print(f"[launch] (dry-run) run_{n:04d}_{c}_{a.build}_r{rep}  env {ms.CONFIGS[c]['env'] or '-'}")
                n += 1
        print("[launch] --dry-run: 何も作らない")
        return
    jobs = []
    lockf = open(os.path.join(HERE, ".registry.lock"), "w")
    fcntl.flock(lockf, fcntl.LOCK_EX)
    new_reg = not os.path.exists(REGISTRY)
    with open(REGISTRY, "a", newline="") as rf:
        w = csv.DictWriter(rf, fieldnames=REG_COLS, delimiter="\t")
        if new_reg:
            w.writeheader()
        n = next_number()
        for rep in a.reps:
            for c in a.cfgs:
                cs = ms.CONFIGS[c]
                kind = cs.get("kind", "forge")
                b = cbin if kind == "convert" else binp
                name = f"run_{n:04d}_{c}_{a.build}_r{rep}"
                rundir = os.path.join(HERE, name)
                populate(rundir, cs["input"])
                env = dict(cs["env"])
                env.update(dict(kv.split("=", 1) for kv in a.env))
                rec = dict(run=name, cfg=c, input=cs["input"], build=a.build, rep=rep, kind=kind, bin=b,
                           bin_sha256=sha256(b), env=" ".join(f"{k}={v}" for k, v in env.items()),
                           created=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
                with open(os.path.join(rundir, "HARNESS.txt"), "w") as f:
                    f.write("\n".join(f"{k:10s}: {v}" for k, v in rec.items()) + "\n")
                with open(os.path.join(rundir, ".state"), "w") as f:
                    f.write("queued\n")
                w.writerow(rec)
                jobs.append(dict(rundir=rundir, kind=kind, bin=b, env=env, memwatch=a.memwatch))
                print(f"[launch] {name}  ({kind}, env {rec['env'] or '-'})")
                n += 1
    fcntl.flock(lockf, fcntl.LOCK_UN)
    spawn_worker(jobs, time.strftime("%Y%m%d_%H%M%S"))


def cmd_resume(a):
    """registry にあって .state が queued の run を (登録どおりの env・バイナリで) 回す。ワーカーが途中で止まったとき用。"""
    jobs = []
    for r in read_registry():
        d = os.path.join(HERE, r["run"])
        if run_state(d) != "queued" or (a.runs and r["run"] not in a.runs):
            continue
        env = dict(kv.split("=", 1) for kv in r["env"].split()) if r["env"] else {}
        # forge は名前 `forge` のリンク経由で起動する (forge_alias。停止前に登録した run は実ファイルのパスを持っている)
        b = forge_alias(r["bin"], r["build"]) if r["kind"] == "forge" else r["bin"]
        jobs.append(dict(rundir=d, kind=r["kind"], bin=b, env=env, memwatch=False))
        print(f"[resume] {r['run']}  ({r['kind']}, env {r['env'] or '-'})")
    if jobs:
        spawn_worker(jobs, "resume_" + time.strftime("%Y%m%d_%H%M%S"))


def cmd_launch_dir(a):
    binp = forge_alias(a.bin, a.label)
    env = dict(kv.split("=", 1) for kv in a.env)
    jobs = []
    for d in a.dirs:
        d = os.path.abspath(os.path.expanduser(d))
        if not os.path.exists(os.path.join(d, "solverConfig.yaml")):
            sys.exit(f"[launch-dir] {d} に solverConfig.yaml が無い")
        with open(os.path.join(d, "INPUT_FILES"), "w") as f:   # 投入前にあったファイル = 入力
            f.write("\n".join(sorted(x for x in os.listdir(d) if not x.startswith("."))) + "\n")
        with open(os.path.join(d, ".state"), "w") as f:
            f.write("queued\n")
        with open(os.path.join(d, "HARNESS.txt"), "w") as f:
            f.write(f"bin       : {binp}\nbin_sha256: {sha256(binp)}\nenv       : {' '.join(f'{k}={v}' for k, v in env.items())}\n"
                    f"memwatch  : {a.memwatch}\ncreated   : {time.strftime('%Y-%m-%dT%H:%M:%S%z')}\n"
                    f"harness   : {os.path.abspath(__file__)}\n")
        jobs.append(dict(rundir=d, kind="forge", bin=binp, env=env, memwatch=a.memwatch))
    spawn_worker(jobs, "dir_" + time.strftime("%Y%m%d_%H%M%S"))


# ---------------------------------------------------------------- ワーカー
def children(pid):
    out = []
    try:
        for t in os.listdir(f"/proc/{pid}/task"):
            with open(f"/proc/{pid}/task/{t}/children") as f:
                for c in f.read().split():
                    out.append(int(c))
                    out += children(int(c))
    except OSError:
        pass
    return out


def find_forge(root_pid, binp):
    real = os.path.realpath(binp)
    for c in children(root_pid):
        try:
            if os.path.realpath(f"/proc/{c}/exe") == real:
                return c
        except OSError:
            pass
    return None


def proc_mem(pid):
    rss = hwm = -1
    try:
        for line in open(f"/proc/{pid}/status"):
            if line.startswith("VmRSS:"):
                rss = int(line.split()[1])
            elif line.startswith("VmHWM:"):
                hwm = int(line.split()[1])
    except OSError:
        pass
    return rss, hwm


def gpu_mem(pid):
    proc_mib = total_mib = -1
    try:
        r = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,used_memory", "--format=csv,noheader,nounits"],
                           capture_output=True, text=True, timeout=20)
        for line in r.stdout.splitlines():
            p, m = [x.strip() for x in line.split(",")]
            if int(p) == pid:
                proc_mib = int(m)
        r = subprocess.run(["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"],
                           capture_output=True, text=True, timeout=20)
        total_mib = int(r.stdout.split()[0])
    except Exception:
        pass
    return proc_mib, total_mib


def nan_check(rundir):
    """plan §6 (c): 残差 CSV の数値列と、run が書いた全 h5 の浮動小数データセットに NaN/Inf が無いか。"""
    import numpy as np
    msgs, bad = [], 0
    p = os.path.join(rundir, "residual_history.csv")
    if os.path.exists(p):
        with open(p) as f:
            rd = csv.reader(f)
            hdr = next(rd, [])
            first_bad = None
            nrow = 0
            for row in rd:
                nrow += 1
                for k, v in zip(hdr, row):
                    if k in ("step", "inner", "phase"):
                        continue
                    try:
                        x = float(v)
                    except ValueError:
                        continue
                    if x != x or x in (float("inf"), float("-inf")):
                        bad += 1
                        if first_bad is None:
                            first_bad = (row[0], k)
            msgs.append(f"residual_history.csv: {nrow} 行, 非有限 {bad} 個" + (f" (最初 step {first_bad[0]} 列 {first_bad[1]})" if first_bad else ""))
    import h5py
    inputs = set(open(os.path.join(rundir, "INPUT_FILES")).read().split()) if os.path.exists(os.path.join(rundir, "INPUT_FILES")) else set()
    nds = 0
    for h in sorted(glob.glob(os.path.join(rundir, "*.h5"))):
        if os.path.basename(h) in inputs:   # 入力 (複製した mesh/初期場) は検査しない
            continue
        nb = 0
        with h5py.File(h, "r") as f:
            def visit(name, obj):
                nonlocal nb, nds
                if isinstance(obj, h5py.Dataset) and obj.dtype.kind == "f":
                    nds += 1
                    a = obj[()]
                    k = int(np.size(a) - np.count_nonzero(np.isfinite(a)))
                    if k:
                        nb += k
                        msgs.append(f"{os.path.basename(h)}:{name} 非有限 {k}")
            f.visititems(visit)
        bad += nb
    msgs.append(f"h5: 浮動小数データセット {nds} 個を検査")
    verdict = "PASS" if bad == 0 else "FAIL"
    with open(os.path.join(rundir, "NANCHECK.txt"), "w") as f:
        f.write("\n".join(msgs) + f"\nNANCHECK: {verdict}\n")
    return verdict


def cmd_worker(a):
    jobs = json.load(open(a.jobfile))
    for j in jobs:
        d = j["rundir"]
        with open(os.path.join(d, ".state"), "w") as f:
            f.write("running\n")
        t0 = time.time()
        env = dict(os.environ)
        env.update(FORGE_CUDA_BLOCKSIZE="128")
        env.update(j["env"])
        print(f"[worker] start {d}", flush=True)
        if j["kind"] == "convert":
            msh = "mesh.msh"
            with open(os.path.join(d, "convert.log"), "w") as lf:
                p = subprocess.Popen([j["bin"], msh, "converted.h5"], cwd=d, stdin=subprocess.DEVNULL, stdout=lf,
                                     stderr=subprocess.STDOUT, env=env)
        else:
            env["FORGE_BIN"] = j["bin"]
            with open(os.path.join(d, "run_case_stdout.log"), "w") as lf:
                p = subprocess.Popen(["bash", RUN_CASE, d], cwd=d, stdin=subprocess.DEVNULL, stdout=lf,
                                     stderr=subprocess.STDOUT, env=env)
        if j.get("memwatch"):
            with open(os.path.join(d, "mem_samples.csv"), "w") as mf:
                mf.write("t_s,forge_pid,VmRSS_kB,VmHWM_kB,gpu_proc_MiB,gpu_total_MiB\n")
                last = (-1, -1)
                while p.poll() is None:
                    fp = find_forge(p.pid, j["bin"])
                    if fp:
                        rss, hwm = proc_mem(fp)
                        g, gt = gpu_mem(fp)
                        if rss > 0:
                            last = (rss, hwm)
                            mf.write(f"{time.time() - t0:.1f},{fp},{rss},{hwm},{g},{gt}\n")
                            mf.flush()
                    time.sleep(1.0)
        rc = p.wait()
        try:
            v = nan_check(d)
        except Exception as e:  # 検査の失敗は合格にしない
            v = f"ERROR ({e})"
            with open(os.path.join(d, "NANCHECK.txt"), "w") as f:
                f.write(f"NANCHECK: ERROR {e}\n")
        extra = ""
        if j["kind"] == "convert":
            # 変換器は終了時の cudaFree で GPUassert (invalid argument) を出して exit 1 になる既知の罠がある (出力 h5 は完全。
            # memory aws-p1-instance-state)。完了は rc でなく出力 h5 が開けて /VALUE を持つことで判定する
            # (cell モードは /VIZMESH を書かないので書込みの行では判定できない)
            ok = 0
            try:
                import h5py
                with h5py.File(os.path.join(d, "converted.h5"), "r") as f:
                    ok = int("VALUE" in f and "MESH" in f)
            except Exception:
                ok = 0
            extra = f" written={ok}"
        with open(os.path.join(d, ".state"), "w") as f:
            f.write(f"done rc={rc} nan={v}{extra} wall={time.time() - t0:.0f}s\n")
        print(f"[worker] done {d} rc={rc} nan={v} ({time.time() - t0:.0f} s)", flush=True)


# ---------------------------------------------------------------- 状態・確認・表
def provenance_ok(rundir, binp):
    p = os.path.join(rundir, "RUN_PROVENANCE.txt")
    if not os.path.exists(p):
        return "no-provenance"
    txt = open(p).read()
    m = re.search(r"forge_bin\s*:\s*(\S+)", txt)
    s = re.search(r"forge_sha256:\s*(\S+)", txt)
    if not m or os.path.realpath(m.group(1)) != os.path.realpath(binp):
        return f"bin-mismatch ({m.group(1) if m else '?'})"
    if not s or s.group(1) != sha256(binp):
        return "sha-mismatch"
    return "ok"


def short_conv(rundir):
    p = os.path.join(rundir, "CONVERGENCE_VERDICT.txt")
    if not os.path.exists(p):
        return "-"
    t = open(p).read()
    m = re.search(r"->\s*([A-Z][A-Z \-]+?)\s*(\(|=|$)", t, re.M)
    return m.group(1).strip() if m else "?"


def cmd_status(a):
    for r in read_registry():
        if a.build and r["build"] != a.build:
            continue
        if a.cfg and r["cfg"] not in a.cfg:
            continue
        d = os.path.join(HERE, r["run"])
        st = run_state(d)
        prov = provenance_ok(d, r["bin"]) if r["kind"] == "forge" and st.startswith("done") else ""
        if st == "running" and not any(os.path.realpath(f"/proc/{p}/cwd") == os.path.realpath(d)
                                       for p in os.listdir("/proc") if p.isdigit() and os.path.exists(f"/proc/{p}/cwd")):
            st = "running? (プロセス無し = 中断。note --exclude して再投入)"
        print(f"{r['run']:44s} {st:48s} {prov} {short_conv(d) if st.startswith('done') else ''}")


def cmd_verify(a):
    """投入後の確認: RUN_PROVENANCE の forge_bin/sha256 が登録どおりか (skill forge-aws-run §2)。"""
    bad = 0
    for r in read_registry():
        if r["kind"] != "forge":
            continue
        d = os.path.join(HERE, r["run"])
        if not run_state(d).startswith("done"):
            continue
        v = provenance_ok(d, r["bin"])
        if v != "ok":
            bad += 1
            print(f"{r['run']}: {v}")
    print(f"VERIFY: {'PASS' if bad == 0 else 'FAIL'} ({bad} 件不一致)")


def cmd_table(a):
    """README の「計算 run 一覧」表の行 (run 名・目的と設定差分・元の入力・主要結果・状態)。"""
    print("| run | 目的・設定差分 | 元の入力 | 主要結果 | 状態 |")
    print("| --- | --- | --- | --- | --- |")
    for r in read_registry():
        d = os.path.join(HERE, r["run"])
        cs = ms.CONFIGS.get(r["cfg"], {})
        src = ms.INPUTS.get(r["input"], {}).get("src", "?")
        st = run_state(d)
        if st.startswith("done"):
            m = re.search(r"rc=(-?\d+) nan=(\S+)", st)
            res = f"rc {m.group(1)}, NaN {m.group(2)}" if m else st
            if r["kind"] == "forge":
                res += f", {short_conv(d)}"
        else:
            res = st
        env = f" + `{r['env']}`" if r["env"] else ""
        state = a.state
        np_ = os.path.join(d, ".note")
        if os.path.exists(np_):
            state, text = open(np_).read().rstrip("\n").split("\t", 1)
            if text:
                res += f"。{text}"
        print(f"| `{r['run']}` | {cs.get('group', r['cfg'])} ({r['build']} r{r['rep']}){env} | "
              f"`inputs/{r['input']}` ← `{src}` | {res} | {state} |")


def cmd_note(a):
    """README の状態列と備考を run に残す。--exclude で比較 (compare_runs.py) から外す (起動失敗など)。"""
    for r in a.runs:
        d = os.path.join(HERE, r)
        with open(os.path.join(d, ".note"), "w") as f:
            f.write(f"{a.state}\t{a.text}\n")
        if a.exclude:
            open(os.path.join(d, ".exclude"), "w").close()
        print(f"[note] {r}: {a.state} {a.text}{' (比較から除外)' if a.exclude else ''}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("seed")
    s.add_argument("names", nargs="*")
    s.add_argument("--bin", required=True)
    s.add_argument("--force-species", action="store_true")
    s.add_argument("--redo", action="store_true")
    s = sub.add_parser("set-ckpt")
    s.add_argument("input")
    s.add_argument("run")
    s.add_argument("--redo", action="store_true")
    s = sub.add_parser("launch")
    s.add_argument("cfgs", nargs="+")
    s.add_argument("--build", required=True)
    s.add_argument("--bin", required=True)
    s.add_argument("--convert-bin")
    s.add_argument("--reps", nargs="+", type=int, default=[1, 2, 3])
    s.add_argument("--env", nargs="*", default=[], help="全構成に足す環境変数 K=V (例 FORGE_MEMLOG=1)")
    s.add_argument("--memwatch", action="store_true")
    s.add_argument("--dry-run", action="store_true")
    s = sub.add_parser("launch-dir")
    s.add_argument("dirs", nargs="+")
    s.add_argument("--bin", required=True)
    s.add_argument("--label", default="dir", help="リンク .bin/<label>/forge の名前 (ビルドごとに変える)")
    s.add_argument("--env", nargs="*", default=[])
    s.add_argument("--memwatch", action="store_true")
    s = sub.add_parser("resume")
    s.add_argument("runs", nargs="*")
    s = sub.add_parser("worker")
    s.add_argument("jobfile")
    s = sub.add_parser("status")
    s.add_argument("--build")
    s.add_argument("--cfg", nargs="*")
    sub.add_parser("verify")
    s = sub.add_parser("table")
    s.add_argument("--state", default="active")
    s = sub.add_parser("note")
    s.add_argument("runs", nargs="+")
    s.add_argument("--state", required=True)
    s.add_argument("--text", default="")
    s.add_argument("--exclude", action="store_true")
    a = ap.parse_args()
    {"seed": cmd_seed, "set-ckpt": cmd_set_ckpt, "launch": cmd_launch, "launch-dir": cmd_launch_dir,
     "worker": cmd_worker, "status": cmd_status, "verify": cmd_verify, "table": cmd_table, "note": cmd_note,
     "resume": cmd_resume}[a.cmd](a)


if __name__ == "__main__":
    main()
