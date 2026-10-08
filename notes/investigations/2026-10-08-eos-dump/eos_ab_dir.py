#!/usr/bin/env python3
"""凍結入力に対する EOS 1 回の全出力のビット比較を、準備済みの入力ディレクトリ 1 つについて回す (AWS の g3 用、ローカルでも可)。

比較器 compare_eos_dump.py は改訂規則 (v2) と旧規則 (v1) の両方を出す。ここでの合否は v2 (plan §5.1 #2 の規則の改訂 (1))。

  prepare: 既存 run から入力を作る (res_*・ログ・判定・時系列を除いて複製 → restart_field で同じ res から保存量を移す →
           nStepOuter 1・outStepInterval 1・output.floorEvents 1 にする。両版で同じファイルを使う: 旧版は floorEvents を読まない)
      python3 eos_ab_dir.py prepare --src-run RUN --res res_20000.h5 --input IN --tools TOOLS --forge NEWBIN
  run: 旧版 2 回・新版 2 回 (入力ディレクトリを毎回複製し、入力ファイルの sha256 が複製元と同じことを確かめてから回す)、
       5 組 (旧同士・新同士・旧対新 3 組) を compare_eos_dump.py で比べ、ダンプは比較が済んだ順に消す (--keep-dumps で残す)。同時に残るダンプは最大 3 本。
      python3 eos_ab_dir.py run --input IN --work W --old OLDBIN --new NEWBIN --tools TOOLS [--env FORGE_CUDA_BLOCKSIZE=128]

床の下の入力の直接比較 (plan §5.1 #2 の追加の試験 (3)、診断専用の再生 FORGE_EOS_REPLAY_FILE):
  base: 旧版 1 回・新版 1 回を再生なしでダンプし、/pre の全配列が両版でビット一致することを確かめて新版のダンプを BASE に残す。
      python3 eos_ab_dir.py base --input IN --work W --old OLDBIN --new NEWBIN --out BASE.h5
  make-replay: BASE の /pre を写し、指定内部節点の roe だけを e_in = e_mix(T_min) + DE [J/kg] になる値にした再生ファイルを作る
       (DE < 0 = 床の下、DE > 0 = 床の上の対照)。余裕 (|DE| ≥ c_v(T_min)·1 K かつ ρ|DE| ≥ 64 ULP(格納 ρE)) と内部節点を確かめる。
      python3 eos_ab_dir.py make-replay --base BASE.h5 --input IN --node 18259 --de -10000 --out R.h5
  run --replay R.h5 --node N --expect-floor 1|0: 再生ありで旧版 2 回・新版 2 回。各ダンプの /pre が R の /pre とビット一致すること、
       実際の pre (ダンプの /pre) で床の条件と余裕、新版の floor_events.csv の step の eos 行で指定節点の温度床の件数を確かめ、
       REPLAY_VERDICT (PASS / FAIL / 試験不成立) を SUMMARY.json と標準出力に出す。
"""
import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import struct
import subprocess
import sys

import numpy as np

SKIP_PREFIX = ("res_",)
SKIP_NAMES = {"forge_run.log", "residual_history.csv", "residual_history.png", "floor_events.csv", "RUN_PROVENANCE.txt",
              "CONVERGENCE_VERDICT.txt", "forge_launches.jsonl", "eos_dump.h5"}
HDF5_LIB = "/usr/lib/x86_64-linux-gnu/hdf5/serial"
RU = 8.314462618               # thermo_d.cuh THERMO_RU
TP_TMIN = 50.0                 # dependentVariables_d.cu DEPVAR_TMIN (TP の温度床)
MAX_INTERVALS = 3              # thermo_d.cuh THERMO_MAX_INTERVALS
SP_FMT = "<5di" + "d" * (MAX_INTERVALS - 1) + "d" * (MAX_INTERVALS * 9) + "2d"   # eosDump.cpp writeThermoDb の詰め方 (詰め物なし)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def value_file(d):
    c = open(os.path.join(d, "solverConfig.yaml")).read()
    return re.search(r'valueFileName:\s*"?([^",}\s]+)', c).group(1)


def mesh_file(d):
    c = open(os.path.join(d, "solverConfig.yaml")).read()
    return re.search(r'meshFileName:\s*"?([^",}\s]+)', c).group(1)


def set_output_on(cfg):
    m = re.search(r"(?m)^output:\s*\{([^}]*)\}", cfg)
    if m:
        inner = m.group(1).strip()
        if "floorEvents" in inner:
            return re.sub(r"floorEvents:\s*\d+", "floorEvents: 1", cfg)
        return cfg[:m.start()] + "output: {" + (inner + ", " if inner else "") + "floorEvents: 1}" + cfg[m.end():]
    if re.search(r"(?m)^output:", cfg):
        raise SystemExit("block 形式の output 節は未対応 (手で output.floorEvents: 1 を入れて --no-edit で回す)")
    return cfg.rstrip("\n") + "\noutput: {floorEvents: 1}\n"


def prepare(a):
    if os.path.exists(a.input):
        raise SystemExit(f"{a.input} は既にある (上書きしない)")
    os.makedirs(a.input)
    copied = []
    for f in sorted(os.listdir(a.src_run)):
        p = os.path.join(a.src_run, f)
        if os.path.isdir(p) or f.startswith(SKIP_PREFIX) or f.endswith(".xmf") or f in SKIP_NAMES or f.endswith(".msh"):
            continue
        shutil.copy2(p, a.input)
        copied.append(f)
    val = value_file(a.input)
    if val not in copied:
        raise SystemExit(f"入力 h5 {val} が {a.src_run} に無い")
    src = os.path.join(a.src_run, a.res)
    env = dict(os.environ, FORGE_BIN=a.forge, LD_LIBRARY_PATH=HDF5_LIB + ":" + os.environ.get("LD_LIBRARY_PATH", ""))
    r = subprocess.run([sys.executable, os.path.join(a.tools, "restart_field.py"), src, os.path.join(a.input, val)],
                       env=env, capture_output=True, text=True)
    out = (r.stdout + r.stderr).strip()
    print(out)
    if r.returncode != 0 or "VERDICT: OK" not in out:
        raise SystemExit("restart_field が OK でない")
    c = open(os.path.join(a.input, "solverConfig.yaml")).read()
    c2 = re.sub(r"nStepOuter:\s*\d+", "nStepOuter: 1", c)
    c2 = re.sub(r"(?m)^(\s*)outStepInterval:\s*\d+", r"\g<1>outStepInterval: 1", c2)
    c2 = set_output_on(c2)
    open(os.path.join(a.input, "solverConfig.yaml"), "w").write(c2)
    try:
        import yaml
        y = yaml.safe_load(c2)
        print("config: nStepOuter", y["time"]["last"]["nStepOuter"], "output", y.get("output"))
    except Exception as e:  # noqa: BLE001
        print("config の YAML 検査をしなかった:", e)
    files = {f: sha(os.path.join(a.input, f)) for f in sorted(os.listdir(a.input)) if os.path.isfile(os.path.join(a.input, f))}
    json.dump({"src_run": a.src_run, "res": a.res, "restart_field": out.splitlines()[-1] if out else "", "files": files},
              open(os.path.join(a.input, "INPUT.json"), "w"), indent=1)
    print(f"入力 {a.input}: {len(files)} ファイル (INPUT.json に sha256)")


def input_files(inp):
    p = os.path.join(inp, "INPUT.json")
    if os.path.exists(p):
        return json.load(open(p))["files"]
    return {f: sha(os.path.join(inp, f)) for f in sorted(os.listdir(inp)) if os.path.isfile(os.path.join(inp, f))}


def floor_rows(rd):
    fe = os.path.join(rd, "floor_events.csv")
    if not os.path.exists(fe):
        return []
    with open(fe) as fh:
        return [r for r in csv.DictReader(fh) if r.get("kind") in ("init", "eos")]


def run_one(a, inp, files, name, binp, extra_env, replay=None):
    rd = os.path.join(a.work, name)
    if os.path.exists(rd):
        shutil.rmtree(rd)
    shutil.copytree(inp, rd)
    bad = [f for f, h in files.items() if f != "INPUT.json" and sha(os.path.join(rd, f)) != h]
    if bad:
        raise SystemExit(f"{name}: 複製した入力が複製元と違う {bad}")
    dump = os.path.join(rd, "eos_dump.h5")
    env = dict(os.environ, FORGE_DUMP_EOS_STEP=str(a.step), FORGE_DUMP_EOS_FILE=dump,
               LD_LIBRARY_PATH=HDF5_LIB + ":" + os.environ.get("LD_LIBRARY_PATH", ""))
    env.pop("FORGE_EOS_REPLAY_FILE", None)
    if replay:
        env["FORGE_EOS_REPLAY_FILE"] = os.path.abspath(replay)
    env.update(extra_env)
    with open(os.path.join(rd, "forge_run.log"), "w") as lf:
        r = subprocess.run([binp], cwd=rd, env=env, stdout=lf, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
    log = open(os.path.join(rd, "forge_run.log")).read()
    ok = r.returncode == 0 and "[eos-dump] wrote" in log and os.path.exists(dump)
    if replay:
        ok = ok and "[eos-dump] replay: loaded" in log and "[eos-dump] replay file checked" in log
    rows = floor_rows(rd)
    # 入力の大きいファイルと初期出力は run から消す (正本は入力ディレクトリ)
    for f in list(files) + [x for x in os.listdir(rd) if x.startswith("res_")]:
        p = os.path.join(rd, f)
        if f not in ("solverConfig.yaml", "INPUT.json") and os.path.isfile(p) and os.path.getsize(p) > (1 << 20):
            os.remove(p)
    info = {"rc": r.returncode, "ok": ok, "dump_sha256": sha(dump) if os.path.exists(dump) else None,
            "dump_bytes": os.path.getsize(dump) if os.path.exists(dump) else 0, "floor_rows": rows,
            "replay_log": [l for l in log.splitlines() if "[eos-dump]" in l]}
    print(f"  {name}: rc={r.returncode} dump={'あり' if os.path.exists(dump) else '無し'} ({info['dump_bytes'] / 2**30:.3f} GiB) "
          f"floor_events {len(rows)} 行", flush=True)
    for row in rows:
        print("     ", ",".join(f"{k}={row[k]}" for k in ("kind", "step", "inner", "q_index", "nT_real", "nRho_real", "nP_real",
                                                         "dRhoE_T_sum", "ids_T", "ids_P", "ids_near")))
    if not ok:
        raise SystemExit(f"{name}: forge が EOS ダンプを書いて終わらなかった (forge_run.log を見る)")
    return dump, info


def compare(a, x, y, dx, dy, res):
    out = os.path.join(a.work, f"cmp_{x}_vs_{y}")
    r = subprocess.run([sys.executable, os.path.join(a.tools, "compare_eos_dump.py"), dx, dy, "--json", out + ".json"],
                       capture_output=True, text=True)
    open(out + ".txt", "w").write(r.stdout + r.stderr)
    j = json.load(open(out + ".json"))
    t = j.get("totals", {})
    res[f"{x}_vs_{y}"] = {"verdict_v1": j.get("verdict_v1"), "verdict_v2": j.get("verdict_v2"),
                          "eos_args_post_diff_bytes": t.get("eos_args_post_diff_bytes"),
                          "eos_args_pre_diff_bytes": t.get("eos_args_pre_diff_bytes"), "eos_args_bytes": t.get("eos_args_bytes"),
                          "other_pre_diff_bytes": t.get("other_pre_diff_bytes"), "other_post_diff_bytes": t.get("other_post_diff_bytes"),
                          "invalid_v1": j.get("invalid_v1"), "invalid_v2": j.get("invalid_v2"), "info_v2": j.get("info_v2"),
                          "different": [d["array"] for d in j.get("different", [])]}
    print(f"  {x} vs {y}: VERDICT_v1 {j.get('verdict_v1')} / VERDICT_v2 {j.get('verdict_v2')}  EOS の配列 post 差分バイト "
          f"{t.get('eos_args_post_diff_bytes')} / {t.get('eos_args_bytes')}  pre 差分 {t.get('eos_args_pre_diff_bytes')}  "
          f"その他 pre/post 差分 {t.get('other_pre_diff_bytes')}/{t.get('other_post_diff_bytes')}", flush=True)
    for s in (j.get("invalid_v2") or [])[:4]:
        print("      不成立 (v2):", s[:160])


# ---------------- 再生 (床の下の入力の直接比較) ----------------

def thermo_db(raw, n):
    """/db/species_thermo (詰め物なしのバイト列) → 種ごとの dict。"""
    size = struct.calcsize(SP_FMT)
    if len(raw) != n * size:
        raise SystemExit(f"物性 DB の長さ {len(raw)} が {n} 種 × {size} バイトでない")
    out = []
    for i in range(n):
        v = struct.unpack_from(SP_FMT, raw, i * size)
        k = 0
        MW, sig, eps, Tlo, Thi, nInt = v[0:6]; k = 6
        Tbrk = list(v[k:k + MAX_INTERVALS - 1]); k += MAX_INTERVALS - 1
        coef = [list(v[k + 9 * j:k + 9 * (j + 1)]) for j in range(MAX_INTERVALS)]; k += 9 * MAX_INTERVALS
        h_datum, invMW = v[k:k + 2]
        out.append({"MW": MW, "Tlo": Tlo, "Thi": Thi, "nInt": nInt, "Tbrk": Tbrk, "coef": coef, "h_datum": h_datum, "invMW": invMW})
    return out


def _pick(sp, Tc):
    k = 0
    for j in range(1, MAX_INTERVALS):
        if not (Tc < sp["Tbrk"][j - 1]):
            k = j
    return sp["coef"][k]


def _cp_molar_c(sp, Tc):
    a = _pick(sp, Tc)
    Ti = 1.0 / Tc
    return RU * (a[0] * Ti * Ti + a[1] * Ti + a[2] + a[3] * Tc + a[4] * Tc ** 2 + a[5] * Tc ** 3 + a[6] * Tc ** 4)


def _h_molar_c(sp, Tc):
    a = _pick(sp, Tc)
    Ti = 1.0 / Tc
    hRT = (-a[0] * Ti * Ti + a[1] * math.log(Tc) * Ti + a[2] + a[3] * Tc / 2.0 + a[4] * Tc ** 2 / 3.0
           + a[5] * Tc ** 3 / 4.0 + a[6] * Tc ** 4 / 5.0 + a[7] * Ti)
    return RU * Tc * hRT


def h_molar(sp, T):
    """thermo_h_molar と同じ (範囲外は端の cp で線形外挿)。独立に書いた Python 版 (比較のための検算)。"""
    if T < sp["Tlo"]:
        return _h_molar_c(sp, sp["Tlo"]) + _cp_molar_c(sp, sp["Tlo"]) * (T - sp["Tlo"])
    if T > sp["Thi"]:
        return _h_molar_c(sp, sp["Thi"]) + _cp_molar_c(sp, sp["Thi"]) * (T - sp["Thi"])
    return _h_molar_c(sp, T)


def cp_molar(sp, T):
    return _cp_molar_c(sp, min(max(T, sp["Tlo"]), sp["Thi"]))


def node_state(pre, attrs, db, node):
    """EOS が節点 node で使う量を float32 の演算順で作り (dependentVariables_d.cu :86-93、TP のハイブリッド組成 :103-114)、
    床の述語の量 (TP: e_in と e_mix(T_min)、CPG: intE/c_v と tMin) を返す。GPU の FMA 縮約とは最後の数 ulp で違いうる (余裕は桁で大きい)。"""
    f32 = np.float32
    roMin = f32(attrs["roMin"])
    ro = f32(pre["ro"][node]); ro_t = max(ro, roMin)
    ux, uy, uz = (f32(pre[k][node]) / ro_t for k in ("roUx", "roUy", "roUz"))
    ek = f32(0.5) * (ux * ux + uy * uy + uz * uz)
    roe = f32(pre["roe"][node])
    intE = roe / ro_t - ek
    st = {"ro": float(ro), "ro_temp": float(ro_t), "ek": float(ek), "roe": float(roe), "e_in": float(intE),
          "ulp_roe": float(np.spacing(np.abs(roe)))}
    if attrs["thermalMethod"] == 2:
        nsr = int(attrs.get("nSpeciesRegistered", 1) or 1)
        n = len(db)
        if nsr >= 2:
            inv_ro = f32(1.0) / ro_t
            yf = [max(f32(pre[f"roY{s}"][node]) * inv_ro, f32(0.0)) for s in range(n)]
            ysum = f32(0.0)
            for y in yf:
                ysum = ysum + y
            inv = f32(1.0) / (ysum if ysum > f32(1e-30) else f32(1e-30))
            Y = [float(y * inv) for y in yf]
        else:
            Y = [1.0]
        R = RU * sum(Y[s] * db[s]["invMW"] for s in range(len(Y))) if db[0]["invMW"] > 0 else RU * sum(Y[s] / db[s]["MW"] for s in range(len(Y)))
        hs = sum(Y[s] * (h_molar(db[s], TP_TMIN) / db[s]["MW"]) for s in range(len(Y)))
        cpm = sum(Y[s] * (cp_molar(db[s], TP_TMIN) / db[s]["MW"]) for s in range(len(Y)))
        st.update(T_min=TP_TMIN, e_floor=hs - R * TP_TMIN, cv_Tmin=cpm - R, R=R, Y=Y)
    else:
        cv32 = f32(attrs["cp"]) / f32(attrs["gamma"])
        tMin = f32(attrs["tMin"])
        st.update(T_min=float(tMin), e_floor=float(cv32) * float(tMin), cv_Tmin=float(cv32), feTq=float(intE / cv32))
    st["margin"] = st["e_in"] - st["e_floor"]     # 負 = 床の下
    return st


def make_replay(a):
    import h5py
    if os.path.exists(a.out):
        raise SystemExit(f"{a.out} は既にある (上書きしない)")
    with h5py.File(os.path.join(a.input, mesh_file(a.input)), "r") as m:
        bnd = set()
        for k in m["BCONDS"]:
            bnd.update(int(i) for i in m[f"BCONDS/{k}/iCells"][:])
    with h5py.File(a.base, "r") as b:
        attrs = {k: (v.decode() if isinstance(v, bytes) else (v.item() if isinstance(v, np.generic) else v)) for k, v in b.attrs.items()}
        pre = {k: b["pre"][k][()] for k in b["pre"]}
        pre_attrs = dict(b["pre"].attrs)
        dbraw = bytes(b["db/species_thermo"][()]) if "db/species_thermo" in b else b""
    if not (0 <= a.node < attrs["nCells"]) or a.node in bnd:
        raise SystemExit(f"節点 {a.node} は内部の実節点でない (nCells {attrs['nCells']}, 境界 {a.node in bnd})")
    db = thermo_db(dbraw, int(attrs.get("db_n_species", 0))) if attrs["thermalMethod"] == 2 else None
    s0 = node_state(pre, attrs, db, a.node)
    # roe を ρ (e_mix(T_min) + ek + DE) にする (密度・運動量・組成・他の全配列はそのまま)
    roe_new = np.float32(s0["ro_temp"] * (s0["e_floor"] + s0["ek"] + a.de))
    pre2 = {k: v.copy() for k, v in pre.items()}
    pre2["roe"][a.node] = roe_new
    s1 = node_state(pre2, attrs, db, a.node)
    chg = [k for k in pre if np.count_nonzero(pre[k].view(np.uint8) != pre2[k].view(np.uint8))]
    nb = int(np.count_nonzero(pre["roe"].view(np.uint32) != pre2["roe"].view(np.uint32)))
    checks = {
        "only_roe_at_node_changed": chg == ["roe"] and nb == 1,
        "abs_margin_ge_cv_1K": abs(s1["margin"]) >= s1["cv_Tmin"] * 1.0,
        "rho_abs_margin_ge_64ulp": s1["ro_temp"] * abs(s1["margin"]) >= 64.0 * s1["ulp_roe"],
        "side": (s1["margin"] < 0) if a.de < 0 else (s1["margin"] > 0),
    }
    print(f"節点 {a.node}: ρ {s1['ro_temp']:.9g}  ek {s1['ek']:.6g}  roe {s0['roe']:.9g} -> {s1['roe']:.9g}  "
          f"e_in {s1['e_in']:.6f}  e_mix(T_min={s1['T_min']:.6g}) {s1['e_floor']:.6f}  e_in − e_mix(T_min) {s1['margin']:.4f} J/kg  "
          f"c_v(T_min) {s1['cv_Tmin']:.4f}  ρ|差| {s1['ro_temp'] * abs(s1['margin']):.6g}  64 ULP {64 * s1['ulp_roe']:.6g}")
    print("検査:", checks)
    if not all(checks.values()):
        raise SystemExit("再生ファイルの条件 (床の側・余裕・1 要素だけの変更) が成り立たない")
    with h5py.File(a.out, "w") as o:
        for k, v in attrs.items():
            o.attrs[k] = v
        g = o.create_group("pre")
        for k, v in pre2.items():
            g.create_dataset(k, data=v)
        for k, v in pre_attrs.items():
            g.attrs[k] = v
        if dbraw:
            o.create_dataset("db/species_thermo", data=np.frombuffer(dbraw, dtype=np.uint8))
        ri = o.create_group("replay_info")
        rec = {"base": os.path.abspath(a.base), "base_sha256": sha(a.base), "node": a.node, "de_J_per_kg": a.de,
               "roe_old": s0["roe"], "roe_new": float(roe_new), "e_in": s1["e_in"], "e_floor": s1["e_floor"], "margin": s1["margin"],
               "cv_Tmin": s1["cv_Tmin"], "rho": s1["ro_temp"], "ek": s1["ek"], "ulp_roe": s1["ulp_roe"], "T_min": s1["T_min"],
               "checks": checks}
        ri.attrs["record_json"] = json.dumps(rec)
    rec["out_sha256"] = sha(a.out)
    json.dump(rec, open(a.out + ".json", "w"), indent=1)
    print(f"再生ファイル {a.out} ({os.path.getsize(a.out) / 2**20:.1f} MiB, sha256 {rec['out_sha256'][:16]}…)")


def base(a):
    import h5py
    inp = os.path.abspath(a.input)
    files = input_files(inp)
    d_old, i_old = run_one(a, inp, files, "base_old", a.old, {})
    d_new, i_new = run_one(a, inp, files, "base_new", a.new, {})
    with h5py.File(d_old, "r") as fo, h5py.File(d_new, "r") as fn:
        no, nn = set(fo["pre"]), set(fn["pre"])
        diff = {k: int(np.count_nonzero(fo["pre"][k][()].view(np.uint8) != fn["pre"][k][()].view(np.uint8))) for k in no & nn}
    bad = sorted(k for k, v in diff.items() if v) + sorted(no ^ nn)
    print(f"base: /pre の配列 {len(no & nn)} 本、旧版と新版で差のある配列 {bad or 'なし'}")
    if bad:
        raise SystemExit("base: 再生なしの /pre が両版で一致しない")
    shutil.move(d_new, a.out)
    os.remove(d_old)
    json.dump({"input": inp, "old": {"bin": a.old, "sha256": sha(a.old), "dump_sha256": i_old["dump_sha256"]},
               "new": {"bin": a.new, "sha256": sha(a.new), "dump_sha256": i_new["dump_sha256"]},
               "pre_arrays": len(no), "pre_diff_arrays": bad, "out": os.path.abspath(a.out)},
              open(a.out + ".json", "w"), indent=1)
    print(f"base: 新版のダンプを {a.out} に残した (旧版のダンプは消した)")


def check_replay_dump(dump, replay, node, attrs_db):
    """ダンプの /pre が再生ファイルの /pre とビット一致するか、実際の pre での床の条件、post の節点の T。"""
    import h5py
    with h5py.File(dump, "r") as d, h5py.File(replay, "r") as r:
        nd, nr = set(d["pre"]), set(r["pre"])
        mism = sorted(nd ^ nr)
        nbytes = 0
        for k in nd & nr:
            x, y = d["pre"][k][()], r["pre"][k][()]
            if x.shape != y.shape or x.dtype != y.dtype:
                mism.append(k)
                continue
            c = int(np.count_nonzero(x.view(np.uint8) != y.view(np.uint8)))
            nbytes += c
            if c:
                mism.append(k)
        attrs = {k: (v.decode() if isinstance(v, bytes) else (v.item() if isinstance(v, np.generic) else v)) for k, v in d.attrs.items()}
        pre = {k: d["pre"][k][()] for k in ("ro", "roUx", "roUy", "roUz", "roe") + tuple(f"roY{s}" for s in range(8)) if k in d["pre"]}
        st = node_state(pre, attrs, attrs_db, node)
        post_T = float(d["post"]["T"][node])
        post_T_bits = int(d["post"]["T"][()].view(np.uint32)[node])
        rep_attr = attrs.get("replay"), attrs.get("replay_file")
    return {"pre_equals_replay": not mism, "pre_mismatch_arrays": mism, "pre_diff_bytes_vs_replay": nbytes,
            "state": {k: v for k, v in st.items() if k != "Y"}, "post_T_node": post_T, "post_T_node_hex": f"0x{post_T_bits:08x}",
            "replay_attrs": rep_attr}


def run(a):
    os.makedirs(a.work, exist_ok=True)
    inp = os.path.abspath(a.input)
    files = input_files(inp)
    extra = dict(kv.split("=", 1) for kv in (a.env or []))
    replay = os.path.abspath(a.replay) if a.replay else None
    S = {"input": inp, "old": {"bin": a.old, "sha256": sha(a.old)}, "new": {"bin": a.new, "sha256": sha(a.new)},
         "env": extra, "step": a.step, "replay": replay, "replay_sha256": sha(replay) if replay else None,
         "node": a.node, "expect_floor": a.expect_floor, "runs": {}, "compare": {}, "replay_checks": {}}
    print(f"old {S['old']['sha256'][:16]}… new {S['new']['sha256'][:16]}… env {extra} replay {replay}")
    db = None
    if replay:
        import h5py
        with h5py.File(replay, "r") as r:
            ra = {k: (v.decode() if isinstance(v, bytes) else (v.item() if isinstance(v, np.generic) else v)) for k, v in r.attrs.items()}
            if ra.get("thermalMethod") == 2:
                db = thermo_db(bytes(r["db/species_thermo"][()]), int(ra.get("db_n_species", 0)))
            S["replay_record"] = json.loads(r["replay_info"].attrs["record_json"]) if "replay_info" in r else None

    def one(name, binp):
        dump, info = run_one(a, inp, files, name, binp, extra, replay)
        S["runs"][name] = info
        if replay:
            S["replay_checks"][name] = c = check_replay_dump(dump, replay, a.node, db)
            st = c["state"]
            print(f"     /pre = 再生ファイル: {c['pre_equals_replay']} (差分バイト {c['pre_diff_bytes_vs_replay']})  実際の pre の節点 {a.node}: "
                  f"e_in {st['e_in']:.6f}  e_mix(T_min) {st['e_floor']:.6f}  差 {st['margin']:.4f} J/kg  post T {c['post_T_node']:.9g} K")
        return dump

    d = {}
    d["old_r1"] = one("old_r1", a.old)
    d["new_r1"] = one("new_r1", a.new)
    compare(a, "old_r1", "new_r1", d["old_r1"], d["new_r1"], S["compare"])
    d["old_r2"] = one("old_r2", a.old)
    compare(a, "old_r1", "old_r2", d["old_r1"], d["old_r2"], S["compare"])
    compare(a, "old_r2", "new_r1", d["old_r2"], d["new_r1"], S["compare"])
    if not a.keep_dumps:
        os.remove(d["old_r2"])
    d["new_r2"] = one("new_r2", a.new)
    compare(a, "new_r1", "new_r2", d["new_r1"], d["new_r2"], S["compare"])
    compare(a, "old_r1", "new_r2", d["old_r1"], d["new_r2"], S["compare"])
    if not a.keep_dumps:
        for k in ("old_r1", "new_r1", "new_r2"):
            os.remove(d[k])
    v = {k: (c["verdict_v1"], c["verdict_v2"]) for k, c in S["compare"].items()}
    print("SUMMARY (v1, v2):", v)
    if replay:
        S["replay_verdict"], S["replay_reasons"] = replay_verdict(a, S)
        print(f"REPLAY_VERDICT: {S['replay_verdict']}")
        for s in S["replay_reasons"]:
            print("   ", s)
    json.dump(S, open(os.path.join(a.work, "SUMMARY.json"), "w"), indent=1, ensure_ascii=False)
    return 0


def replay_verdict(a, S):
    """合否 (plan §5.1 #2 の追加の試験 (3)): 試験不成立 = 入力不一致・非有限の出力 (v2 INVALID)・床の条件の不成立・腕内の非再現。
    FAIL = 旧版と新版の間で再現する出力差 (v2 DIFFERENT)、または新版の温度床の件数が期待と違う。PASS = それ以外。"""
    inval, fail, notes = [], [], []
    for name, c in S["replay_checks"].items():
        if not c["pre_equals_replay"]:
            inval.append(f"{name}: ダンプの /pre が再生ファイルと違う {c['pre_mismatch_arrays'][:6]} (入力不一致)")
        st = c["state"]
        cv1 = st["cv_Tmin"] * 1.0
        ulp64 = 64.0 * st["ulp_roe"] / st["ro_temp"]
        if a.expect_floor == 1 and not (st["margin"] < 0 and -st["margin"] >= cv1 and -st["margin"] >= ulp64):
            inval.append(f"{name}: 実際の pre で床下の条件または余裕が成り立たない (差 {st['margin']:.4f} J/kg, c_v·1K {cv1:.4f}, 64 ULP/ρ {ulp64:.4g})")
        if a.expect_floor == 0 and not (st["margin"] > 0 and st["margin"] >= cv1 and st["margin"] >= ulp64):
            inval.append(f"{name}: 実際の pre で床の上の条件または余裕が成り立たない (差 {st['margin']:.4f} J/kg)")
    for k, c in S["compare"].items():
        same_arm = k.split("_")[0] == k.split("_vs_")[1].split("_")[0]
        if c["verdict_v2"] == "INVALID":
            inval.append(f"{k}: v2 INVALID {(c['invalid_v2'] or [''])[0][:120]}")
        elif c["verdict_v2"] == "DIFFERENT":
            (inval if same_arm else fail).append(f"{k}: v2 DIFFERENT {c['different']} ({'腕内の非再現' if same_arm else '旧版と新版の出力差'})")
    for name in ("new_r1", "new_r2"):
        rows = [r for r in S["runs"][name]["floor_rows"] if r["kind"] == "eos" and int(r["step"]) == a.step]
        if len(rows) != 1:
            fail.append(f"{name}: floor_events.csv に step {a.step} の eos 行が 1 行でない ({len(rows)} 行)")
            continue
        r = rows[0]
        nT, ids = int(r["nT_real"]), [x for x in r["ids_T"].split(";") if x]
        want_ids = [str(a.node)] if a.expect_floor == 1 else []
        if nT != a.expect_floor or ids != want_ids:
            fail.append(f"{name}: 対象 EOS 行の温度床 {nT} 件 ids {ids} (期待 {a.expect_floor} 件 {want_ids})")
        notes.append(f"{name}: eos 行 step {r['step']} inner {r['inner']} q_index {r['q_index']}: nT_real {r['nT_real']} ids_T '{r['ids_T']}' "
                     f"nRho_real {r['nRho_real']} nP_real {r['nP_real']} ids_P '{r['ids_P']}' dRhoE_T_sum {r['dRhoE_T_sum']} "
                     f"n_near_real {r['n_near_real']} ids_near '{r['ids_near']}' nT_mismatch_real {r['nT_mismatch_real']} overflow {r['overflow']}")
    for name in ("old_r1", "old_r2"):
        if S["runs"][name]["floor_rows"]:
            notes.append(f"{name}: 旧版に floor_events.csv の行がある (想定外: 旧版はカウンタを持たない)")
    verdict = "試験不成立" if inval else ("FAIL" if fail else "PASS")
    return verdict, inval + fail + notes


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--src-run", required=True)
    p.add_argument("--res", required=True)
    p.add_argument("--input", required=True)
    p.add_argument("--tools", required=True, help="restart_field.py のある tools ディレクトリ")
    p.add_argument("--forge", required=True, help="restart_field の化学種解決に使う forge (新版)")
    for name in ("run", "base"):
        r = sub.add_parser(name)
        r.add_argument("--input", required=True)
        r.add_argument("--work", required=True)
        r.add_argument("--old", required=True)
        r.add_argument("--new", required=True)
        r.add_argument("--step", type=int, default=1)
        if name == "run":
            r.add_argument("--tools", required=True, help="compare_eos_dump.py のある tools ディレクトリ")
            r.add_argument("--env", action="append", help="KEY=VALUE (両版に同じ値を渡す。例 FORGE_CUDA_BLOCKSIZE=128)")
            r.add_argument("--keep-dumps", action="store_true")
            r.add_argument("--replay", default=None, help="再生ファイル (make-replay の出力)。両版に FORGE_EOS_REPLAY_FILE で渡す")
            r.add_argument("--node", type=int, default=None, help="再生で roe を変えた節点")
            r.add_argument("--expect-floor", type=int, choices=(0, 1), default=None, help="新版の対象 EOS 行で期待する温度床の件数")
        else:
            r.add_argument("--out", required=True)
    m = sub.add_parser("make-replay")
    m.add_argument("--base", required=True)
    m.add_argument("--input", required=True, help="格子 (境界節点の判定) を読む入力ディレクトリ")
    m.add_argument("--node", type=int, required=True)
    m.add_argument("--de", type=float, required=True, help="e_in − e_mix(T_min) [J/kg] (負 = 床の下)")
    m.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.cmd == "run" and a.replay and (a.node is None or a.expect_floor is None):
        ap.error("--replay には --node と --expect-floor が要る")
    return {"prepare": prepare, "run": run, "base": base, "make-replay": make_replay}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main() or 0)
