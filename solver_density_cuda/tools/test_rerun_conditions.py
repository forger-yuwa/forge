#!/usr/bin/env python3
"""rerun_conditions.py と run_staged_ns の段終了ゲートの単体試験 (plan tooling-rerun-conditions §6「単体」(a)〜(m))。

    python3 solver_density_cuda/tools/test_rerun_conditions.py

fixture は case/45 run_0094 (粗格子 NS) の入力を一時ディレクトリへ複製して使う (元は書き換えない)。
`RERUN_TEST_REF` で別の参照 run を渡せる。3 種・トレーサ・等温壁・凝縮・inletProfile 等は run_0094 の config を書き換えて合成する。

`restart_field.py` まで通す試験 ((a) の作成経路と (j) の作成経路) は `forge --resolve-species` を要する。
--resolve-species を持つ forge が無い (既定の build/forge が旧版で FORGE_BIN も無い) ときは **理由つきで SKIP** し、
`--force-species` で通すことはしない。SKIP は FAIL 件数に数えないが、末尾に件数と理由を出す。
"""
import contextlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace

import h5py
import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "design"))
import forge_species as fsp  # noqa: E402
import rerun_conditions as rc  # noqa: E402

REF_SRC = os.environ.get("RERUN_TEST_REF",
                         "/home/sano/work/forge/case/45.isobutane_m6_d155/run_0094_ns_c2pin_pass2_ext6k")
RES = "res_6000.h5"
FIXTURE_GLOBS = ("bcondConfig.yaml", "solverConfig.yaml", "species_meta.yaml", "resolved_species_*.yaml",
                 "prepare_info.json", "probe.yaml", "nozzle.h5", "nozzle.xmf", RES, "res_outlet_2_6000.h5",
                 "wall_*.csv", "target_axis_M.csv", "delta_r_initial.*",
                 # 持ち込まれてはいけないもの (許可リスト外) も 1 つずつ置く
                 "CONVERGENCE_VERDICT.txt", "RUN_PROVENANCE.txt")

fails = 0
skips = []


def check(name, ok, info=""):
    global fails
    print(("ok   " if ok else "FAIL ") + name + (f" ({info})" if info else ""))
    fails += (not ok)


def skip(name, why):
    skips.append((name, why))
    print(f"SKIP {name} — {why}")


TMP = tempfile.mkdtemp(prefix="rerun_cond_test_")
BASE = os.path.join(TMP, "base_ref")


def build_base():
    import fnmatch
    os.makedirs(BASE)
    for fn in sorted(os.listdir(REF_SRC)):
        if any(fnmatch.fnmatch(fn, g) for g in FIXTURE_GLOBS):
            shutil.copy2(os.path.join(REF_SRC, fn), os.path.join(BASE, fn))


_n = [0]


def make_ref(edit_cfg=None, edit_bc=None):
    """BASE の複製 (実ファイル) を作り、config を文字列置換で書き換える。戻り (ref, new のパス)。"""
    _n[0] += 1
    ref = os.path.join(TMP, f"ref{_n[0]:02d}")
    shutil.copytree(BASE, ref)
    for fn, fx in (("solverConfig.yaml", edit_cfg), ("bcondConfig.yaml", edit_bc)):
        if fx:
            p = os.path.join(ref, fn)
            t = open(p).read()
            t2 = fx(t)
            assert t2 != t, f"{fn} の書き換えが当たらない"
            open(p, "w").write(t2)
    return ref, os.path.join(TMP, f"new{_n[0]:02d}")


def done(ref, new=None):
    shutil.rmtree(ref, ignore_errors=True)
    if new:
        shutil.rmtree(new, ignore_errors=True)


def run_main(argv):
    """rc.main を同一プロセスで回し (戻り値, stdout, stderr) を返す。"""
    o, e = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(o), contextlib.redirect_stderr(e):
        code = rc.main(argv)
    return code, o.getvalue(), e.getvalue()


def plan_of(argv):
    return rc.build_plan(rc.make_parser().parse_args(argv))


def forge_for_resolve():
    try:
        exe = fsp.find_forge(None)          # FORGE_BIN > solver_density_cuda/build/forge
    except fsp.SpeciesResolveUnavailable as e:
        return None, str(e)
    if exe is None:
        return None, "既定の solver_density_cuda/build/forge が --resolve-species を持たない旧版 (FORGE_BIN も未指定)"
    return exe, None


def lines_diff(a, b):
    la, lb = a.splitlines(), b.splitlines()
    return [i for i in range(max(len(la), len(lb))) if (la[i] if i < len(la) else None) != (lb[i] if i < len(lb) else None)]


# --- 合成 fixture の書き換え ---
def cfg_three_species_tracer(t):
    t = t.replace('"H2O"],', '"H2O", "CO2"],', 1)
    return t.replace("thermoHrefTemp: 298.15}", "thermoHrefTemp: 298.15, tracer: exhaust}", 1)


def bc_three_species(t):
    return t.replace("Y0: 0.91420000, Y1: 0.08580000,", "Y0: 0.91420000, Y1: 0.08580000, Y2: 0.0,", 1)


def meta_three_species(ref):
    """species_meta.yaml を 3 種 (MIXDRY, H2O, CO2) + トレーサに合わせる。"""
    p = os.path.join(ref, "species_meta.yaml")
    m = yaml.safe_load(open(p))
    m["species"] = ["MIXDRY", "H2O", "CO2"]
    m["expansion"]["CO2"] = {"CO2": 1.0}
    m["streams"]["inflow"]["Y_transport"] = list(m["streams"]["inflow"]["Y_transport"]) + [0.0]
    m["tracer"]["enabled"] = True
    open(p, "w").write(yaml.safe_dump(m, sort_keys=False, allow_unicode=True))


def add_three_species_fields(ref):
    """res と nozzle.h5 に roY2 (= 0、一部だけ正) と roXi を足す (ΣρY = ρ を保つよう roY0 から差し引く)。"""
    for fn in (RES, "nozzle.h5"):
        with h5py.File(os.path.join(ref, fn), "r+") as f:
            V = f["VALUE"]
            ro = np.asarray(V["ro"], dtype=np.float64)
            y2 = np.zeros_like(ro)
            y2[: ro.size // 3] = 0.01 * ro[: ro.size // 3]          # 1/3 は正、残りはゼロ
            V["roY0"][...] = (np.asarray(V["roY0"], dtype=np.float64) - y2).astype(np.float32)
            V.create_dataset("roY2", data=y2.astype(np.float32))
            xi = 0.3 * ro
            xi[ro.size // 2:] = 0.0                                   # 半分はゼロ
            V.create_dataset("roXi", data=xi.astype(np.float32))


def main():
    if not os.path.isdir(REF_SRC):
        print(f"参照 run {REF_SRC} が無い (RERUN_TEST_REF で指定)")
        return 2
    build_base()
    forge, forge_why = forge_for_resolve()
    print(f"fixture: {REF_SRC} -> {BASE}")
    print(f"forge (--resolve-species): {forge or 'なし — ' + forge_why}")

    # ------------------------------------------------------------------ (a)
    ref, new = make_ref()
    p = plan_of([ref, new])
    check("(a) 無変更: 計画の bcond・solverConfig が参照とバイト一致",
          p["new_bc_text"] == p["bc_text"] and p["new_cfg_text"] == p["cfg_text"])
    p2 = plan_of([ref, new, "--steps", "12000", "--out-interval", "500"])
    exp = p2["cfg_text"].replace("nStepOuter: 6000", "nStepOuter: 12000").replace("outStepInterval: 1000", "outStepInterval: 500")
    check("(a) --steps/--out-interval: solverConfig の差はその 2 トークンだけ・bcond はバイト一致",
          p2["new_cfg_text"] == exp and p2["new_bc_text"] == p2["bc_text"])
    code, _, err = run_main([ref, new, "--steps", "12000", "--out-interval", "5000"])
    check("(a) nStepOuter % outStepInterval != 0 は停止", code == 2 and not os.path.exists(new) and "倍数でない" in err)
    if forge is None:
        skip("(a) 作成経路 (restart_field VERDICT OK・必要保存量 array_equal・species_hash 一致)", forge_why)
    else:
        code, out, err = run_main([ref, new, "--steps", "12000", "--out-interval", "1000"])
        check("(a) 作成: 終了コード 0", code == 0, err[-400:])
        if code == 0:
            rec = json.load(open(os.path.join(new, "RERUN_CONDITIONS.json")))
            check("(a) restart_field の VERDICT OK 行を記録", rec["restart_field_verdict"].startswith("VERDICT: OK"))
            check("(a) bcond がバイト一致", open(os.path.join(new, "bcondConfig.yaml"), "rb").read()
                  == open(os.path.join(ref, "bcondConfig.yaml"), "rb").read())
            with h5py.File(os.path.join(ref, RES), "r") as s, h5py.File(os.path.join(new, "nozzle.h5"), "r") as d:
                eq = all(np.array_equal(np.asarray(s["VALUE"][k]), np.asarray(d["VALUE"][k])) for k in rec["required_conserved"])
                hs, hd = s.attrs.get("species_hash"), d.attrs.get("species_hash")
            check("(a) 必要保存量が参照 res と array_equal", eq)
            check("(a) species_hash 一致", hs is not None and hs == hd, f"{hs} / {hd}")
            copied = set(os.listdir(new))
            check("(a) 許可リスト外 (VERDICT・PROVENANCE・res) を持ち込まない",
                  not ({"CONVERGENCE_VERDICT.txt", "RUN_PROVENANCE.txt", RES} & copied), sorted(copied))
            pi = json.load(open(os.path.join(new, "prepare_info.json")))
            check("(a) prepare_info に rerun_of・ic_from", pi.get("rerun_of") and pi.get("ic_from", "").endswith(RES))
    done(ref, new)

    # ------------------------------------------------------------------ (b)
    ref, new = make_ref()
    p = plan_of([ref, new, "--Pt", "4.4e6", "--Ps", "1789.6"])
    a, b = p["bc_text"], p["new_bc_text"]
    dl = lines_diff(a, b)
    la, lb = a.splitlines(), b.splitlines()
    inl = [i for i, s in enumerate(la) if s.startswith("inlet:")][0]
    outl = [i for i, s in enumerate(la) if s.startswith("outlet:")][0]
    check("(b) 差分のある行は inlet と outlet だけ", sorted(dl) == sorted([inl, outl]), dl)
    check("(b) inlet 行の差は Pt トークンだけ",
          lb[inl] == la[inl].replace("Pt: 5500000.0", "Pt: 4400000.0"), lb[inl])
    check("(b) outlet 行の差は Ps・Pt トークンだけ (Tt 据え置き)",
          lb[outl] == la[outl].replace("Ps: 2237.0", "Ps: 1789.6").replace("Pt: 2237.0", "Pt: 1789.6"), lb[outl])
    y = yaml.safe_load(b)
    check("(b) YAML 再読込で要求値", y["inlet"]["floats"]["Pt"] == 4.4e6 and y["outlet"]["floats"]["Ps"] == 1789.6
          and y["outlet"]["floats"]["Pt"] == 1789.6 and y["outlet"]["floats"]["Tt"] == 300.0)
    check("(b) 他の行はバイト一致", all(la[i] == lb[i] for i in range(len(la)) if i not in (inl, outl)))
    check("(b) Ps/(f·P_exit_ref) を記録 (= 1789.6/(0.8·P_exit_ref))", p["Ps_over_fPexit"] is not None
          and abs(p["Ps_over_fPexit"] - 1789.6 / (0.8 * p["P_exit_ref"])) < 1e-12, p["Ps_over_fPexit"])
    check("(b) P_exit_ref は出口断面の内部節点の P (課した Ps 2237 ちょうどでない; 2238〜2246)",
          p["P_exit_ref"] != 2237.0 and 2238.0 < p["P_exit_ref"] < 2246.0 and "内部節点" in p["P_exit_ref_source"],
          (p["P_exit_ref"], p["P_exit_ref_source"]))
    done(ref, new)
    # 壁・軸の BC 節点を除いていること: 出口列の端 (壁・軸) の P を極端な値にしても P_exit_ref は動かない
    ref, new = make_ref()
    with h5py.File(os.path.join(ref, "nozzle.h5"), "r") as f:
        o = set(np.asarray(f["BCONDS/2/iCells"]).tolist())
        ends = sorted(o & (set(np.asarray(f["BCONDS/3/iCells"]).tolist()) | set(np.asarray(f["BCONDS/4/iCells"]).tolist())))
    with h5py.File(os.path.join(ref, RES), "r+") as f:
        a = np.asarray(f["VALUE/P"]); a[ends] = 1e9; f["VALUE/P"][...] = a
    p_end = plan_of([ref, new, "--Pt", "4.4e6", "--Ps", "1789.6"])
    check("(b) 壁・軸の節点を除外 (端の P を 1e9 にしても P_exit_ref 不変)",
          len(ends) == 2 and p_end["P_exit_ref"] == p["P_exit_ref"], (ends, p_end["P_exit_ref"]))
    done(ref, new)
    # 出口の BC 節点が単一 x に並ばないとき: 出口 BC の節点と同じ x を持つ節点で代替
    ref, new = make_ref()
    with h5py.File(os.path.join(ref, "nozzle.h5"), "r+") as f:
        o = np.asarray(f["BCONDS/2/iCells"])
        C = np.asarray(f["MESH/COORD"]).reshape(-1, 3)
        C[o[50], 0] -= 1e-3                      # 出口の内部節点 1 つだけ x をずらす
        f["MESH/COORD"][...] = C.reshape(-1)
    p_alt = plan_of([ref, new, "--Pt", "4.4e6", "--Ps", "1789.6"])
    check("(b) 出口が単一 x に並ばない → 同じ x の節点で代替 (値は出口列とほぼ同じ)",
          p_alt["P_exit_ref"] is not None and "出口 BC の節点と同じ x" in p_alt["P_exit_ref_source"] and abs(p_alt["P_exit_ref"] - p["P_exit_ref"]) < 5.0,
          (p_alt["P_exit_ref"], p_alt["P_exit_ref_source"]))
    done(ref, new)
    # 出口の BC 節点が取れないとき: 停止せず P_exit_ref = null と警告を記録
    ref, new = make_ref()
    with h5py.File(os.path.join(ref, "nozzle.h5"), "r+") as f:
        del f["BCONDS/2/iCells"]
    p_null = plan_of([ref, new, "--Pt", "4.4e6", "--Ps", "1789.6"])
    check("(b) P_exit_ref が取れない → 停止せず null・警告を記録",
          p_null["P_exit_ref"] is None and p_null["Ps_over_fPexit"] is None
          and any("P_exit_ref: null" in w for w in p_null["warnings"]), p_null["warnings"])
    done(ref, new)

    # ------------------------------------------------------------------ (c)
    ref, new = make_ref()
    r = subprocess.run([sys.executable, os.path.join(HERE, "rerun_conditions.py"), ref, new, "--lump", "CO2=0.2"],
                       capture_output=True, text=True)
    check("(c) --lump → 非ゼロ終了・NEW_RUN なし", r.returncode != 0 and not os.path.exists(new), r.stderr.strip()[-200:])
    done(ref, new)

    # ------------------------------------------------------------------ (d)
    ref, new = make_ref(edit_bc=lambda t: t.replace("kind: inlet_Pressure", "kind: inlet_uniformVelocity", 1))
    code, _, err = run_main([ref, new])
    check("(d) inlet_uniformVelocity → 拒否", code == 2 and not os.path.exists(new) and "対応外の境界種別" in err)
    done(ref, new)

    # ------------------------------------------------------------------ (e)
    ref, new = make_ref()
    p = plan_of([ref, new, "--Y", "H2O=0.09"])
    yb = yaml.safe_load(p["new_bc_text"])["inlet"]["floats"]
    check("(e) --Y H2O=0.09 → Y0 = 0.91・Y1 = 0.09", yb["Y0"] == 0.91 and yb["Y1"] == 0.09, (yb["Y0"], yb["Y1"]))
    check("(e) ΣY = 1 (1e-12)", abs(yb["Y0"] + yb["Y1"] - 1.0) <= 1e-12)
    meta = yaml.safe_load(p["meta_text_new"])
    check("(e) species_meta Y_transport = [0.91, 0.09]", meta["streams"]["inflow"]["Y_transport"] == [0.91, 0.09],
          meta["streams"]["inflow"]["Y_transport"])
    yr = meta["streams"]["inflow"]["Y"]
    check("(e) species_meta 実種 Y を expansion で同期 (H2O = 0.09, Σ = 1)",
          abs(yr["H2O"] - 0.09) < 1e-15 and abs(sum(yr.values()) - 1.0) < 1e-12, yr)
    p1 = plan_of([ref, new, "--Y1", "0.09"])
    check("(e) --Y1 0.09 は --Y H2O=0.09 と同じ bcond", p1["new_bc_text"] == p["new_bc_text"])
    check("(e) 組成変更 → recommended_stages full", p["recommended_stages"]["stages"] == "full")
    done(ref, new)
    ref, new = make_ref(edit_cfg=cfg_three_species_tracer, edit_bc=bc_three_species)
    add_three_species_fields(ref)
    meta_three_species(ref)
    code, _, err = run_main([ref, new, "--Y", "H2O=0.09", "--dry-run"])
    check("(e) 3 種で --balance なし → 停止", code == 2 and "--balance" in err and not os.path.exists(new), err[-200:])
    p3 = plan_of([ref, new, "--Y", "H2O=0.09", "--balance", "MIXDRY"])
    y3 = yaml.safe_load(p3["new_bc_text"])["inlet"]["floats"]
    check("(e) 3 種 + --balance MIXDRY → Y0 = 1 − 0.09 − Y2、Σ = 1 (1e-12)",
          abs(y3["Y0"] + y3["Y1"] + y3["Y2"] - 1.0) <= 1e-12 and y3["Y1"] == 0.09 and y3["Y2"] == 0.0, y3)
    done(ref, new)

    # ------------------------------------------------------------------ (f)
    iso = lambda t: re.sub(r"wall:\s*\{physID: 3, kind: wall, +outputHDFflg: 1, ints: , floats: \}",  # noqa: E731
                           "wall:   {physID: 3, kind: wall_isothermal,  outputHDFflg: 1, ints: , "
                           "floats: {Ux: 0.0, Uy: 0.0, Uz: 0.0, Ts: 500.0}}", t)
    ref, new = make_ref(edit_bc=iso)
    code, _, err = run_main([ref, new, "--Tt", "1500"])
    check("(f) 等温壁 + --Tt のみ → 停止", code == 2 and "--Tw" in err and not os.path.exists(new), err[-200:])
    p = plan_of([ref, new, "--Tt", "1500", "--Tw", "450"])
    yb = yaml.safe_load(p["new_bc_text"])
    check("(f) --Tw 450 → 等温壁の Ts = 450 (他の float は据え置き)",
          yb["wall"]["floats"] == {"Ux": 0.0, "Uy": 0.0, "Uz": 0.0, "Ts": 450.0} and yb["inlet"]["floats"]["Tt"] == 1500.0)
    p = plan_of([ref, new, "--Tt", "1500", "--keep-Tw"])
    check("(f) --keep-Tw → 壁行はバイト一致", lines_diff(p["bc_text"], p["new_bc_text"]) == [0])
    done(ref, new)

    # ------------------------------------------------------------------ (f′) Euler の滑り壁 (2026-10-06 追加)
    slip = lambda t: re.sub(r"wall:\s*\{physID: 3, kind: wall, +outputHDFflg: 1, ints: , floats: \}",  # noqa: E731
                            "wall:   {physID: 3, kind: slip,             outputHDFflg: 1, ints: , floats: }", t)
    ref, new = make_ref(edit_bc=slip)
    p = plan_of([ref, new, "--Pt", "4.4e6", "--Ps", "1789.6"])
    yb = yaml.safe_load(p["new_bc_text"])
    check("(f′) 滑り壁 (slip) の run を受理し、壁行はバイト一致", yb["wall"]["kind"] == "slip" and lines_diff(p["bc_text"], p["new_bc_text"]) == [0, 1],
          lines_diff(p["bc_text"], p["new_bc_text"]))
    check("(f′) slip 壁の節点は P_exit_ref から除かれる (内部節点 95 点)", "95/97" in str(p.get("P_exit_ref_source", "")), p.get("P_exit_ref_source"))
    done(ref, new)

    # ------------------------------------------------------------------ (f″) 乱流モデルなしの run に残る roK/roOmega (2026-10-06 追加)
    noturb = lambda t: re.sub(r"^turbulence:.*$", "turbulence: {model: \"none\"}", t, flags=re.M)  # noqa: E731
    ref, new = make_ref(edit_cfg=noturb, edit_bc=slip)
    p = plan_of([ref, new, "--Pt", "4.4e6", "--Ps", "1789.6", "--scale-ic", "pt"])
    check("(f″) 乱流なし + roK/roOmega の入れ物 → 受理し警告に記録、必要保存量に roK/roOmega を含めない",
          "roK" not in p["required"] and any("未使用量" in w for w in p["warnings"]), (p["required"], p["warnings"]))
    done(ref, new)

    # ------------------------------------------------------------------ (n) Pt 変更の推奨 (2026-10-06、§6 (ii′)・A3)
    ref, new = make_ref()
    p = plan_of([ref, new, "--Pt", "4.4e6", "--Ps", "1789.6", "--scale-ic", "pt"])
    check("(n) Pt 変更 → recommended_stages full・本段 cfl 1", p["recommended_stages"].get("stages") == "full" and p["recommended_stages"].get("cfl") == 1.0, p["recommended_stages"])
    check("(n) scale-ic pt なら none 警告なし", not any("--scale-ic none" in w for w in p["warnings"]), p["warnings"])
    p = plan_of([ref, new, "--Pt", "4.4e6", "--Ps", "1789.6", "--scale-ic", "none"])
    check("(n) Pt 変更 + scale-ic none → 警告", any("--scale-ic none" in w for w in p["warnings"]), p["warnings"])
    done(ref, new)

    # ------------------------------------------------------------------ (g)
    ref, new = make_ref(edit_cfg=lambda t: t + "condensation: {condensation: 1, nCondSpecies: 1, condensationSpecies: H2O}\n")
    code, _, err = run_main([ref, new, "--Pt", "4.4e6", "--Ps", "1789.6", "--scale-ic", "pt"])
    check("(g) 凝縮 block + --scale-ic pt → 停止", code == 2 and "--scale-ic pt は使えない" in err and "condensation" in err
          and not os.path.exists(new), err[-200:])
    done(ref, new)
    ref, new = make_ref()
    code, _, err = run_main([ref, new, "--Pt", "4.4e6", "--Ps", "1789.6", "--Tt", "1500", "--scale-ic", "pt"])
    check("(g') Tt 変更 + --scale-ic pt → 停止", code == 2 and "Tt を変える" in err)
    code, _, err = run_main([ref, new, "--scale-ic", "pt"])
    check("(g') Pt 変更なし + --scale-ic pt → 停止", code == 2 and "Pt を変えていない" in err)
    done(ref, new)

    # ------------------------------------------------------------------ (h)
    ref, new = make_ref()
    code, _, err = run_main([ref, new, "--forge", os.path.join(TMP, "no_such_forge")])
    check("(h) restart_field を失敗させる → 非ゼロ・NEW_RUN が残らない",
          code == 1 and not os.path.exists(new) and "restart_field.py が失敗" in err, err[-200:])
    done(ref, new)

    # ------------------------------------------------------------------ (i')
    ref, new = make_ref()
    with h5py.File(os.path.join(ref, RES), "r+") as f:
        del f["VALUE/roOmega"]
    code, _, err = run_main([ref, new])
    check("(i') SRC に roOmega が無い → 拒否", code == 2 and "roOmega" in err and not os.path.exists(new), err[-200:])
    done(ref, new)
    ref, new = make_ref()
    with h5py.File(os.path.join(ref, RES), "r+") as f:
        a = np.asarray(f["VALUE/roUx"]); a[10] = np.nan; f["VALUE/roUx"][...] = a
    code, _, err = run_main([ref, new])
    check("(i') SRC に非有限 → 拒否", code == 2 and "非有限" in err and not os.path.exists(new))
    done(ref, new)
    ref, new = make_ref()
    with h5py.File(os.path.join(ref, "nozzle.h5"), "r+") as f:
        f["VALUE"].create_dataset("P", data=np.ones(f["VALUE/ro"].shape, np.float32))
    code, _, err = run_main([ref, new])
    check("(i') DST /VALUE に集合外の量 → 拒否", code == 2 and "以外がある" in err and not os.path.exists(new))
    done(ref, new)

    # ------------------------------------------------------------------ (j)
    ref, new = make_ref(edit_cfg=cfg_three_species_tracer, edit_bc=bc_three_species)
    add_three_species_fields(ref)
    meta_three_species(ref)
    cfg = yaml.safe_load(open(os.path.join(ref, "solverConfig.yaml")))
    req = rc.required_conserved_from_cfg(cfg)
    check("(j) 3 種 + トレーサの必要保存量に roY2・roXi", "roY2" in req and "roXi" in req, req)
    p = plan_of([ref, new, "--Pt", "4.4e6", "--Ps", "1789.6", "--scale-ic", "pt"])
    check("(j) 3 種 + トレーサで --scale-ic pt の計画が通る (f = 0.8)", abs(p["f"] - 0.8) < 1e-15)
    dst = os.path.join(TMP, "scale_dst.h5")
    shutil.copy2(os.path.join(ref, "nozzle.h5"), dst)
    rel = rc.scale_fields(dst, os.path.join(ref, RES), req, 0.8)
    with h5py.File(os.path.join(ref, RES), "r") as s, h5py.File(dst, "r") as d:
        ok_all = all(np.allclose(np.asarray(d["VALUE"][k], np.float64), 0.8 * np.asarray(s["VALUE"][k], np.float64),
                                 rtol=1e-6, atol=0) for k in req)
        z2 = np.asarray(s["VALUE/roY2"]) == 0
        zx = np.asarray(s["VALUE/roXi"]) == 0
        zero_ok = bool(np.all(np.asarray(d["VALUE/roY2"])[z2] == 0) and np.all(np.asarray(d["VALUE/roXi"])[zx] == 0))
    with h5py.File(os.path.join(ref, "nozzle.h5"), "r") as o, h5py.File(dst, "r") as d:
        wd = np.array_equal(np.asarray(d["VALUE/wall_dist"]), np.asarray(o["VALUE/wall_dist"]))
    check("(j) scale_fields: 必要保存量 (roY2・roXi 含む) が f 倍 (rtol 1e-6, atol 0)", ok_all, rel)
    check("(j) ゼロ成分はゼロのまま (roY2・roXi)", zero_ok and z2.any() and zx.any())
    check("(j) wall_dist は触らない", wd)
    os.remove(dst)
    if forge is None:
        skip("(j) 作成経路 (restart_field → scale → RERUN_CONDITIONS.json の検査)", forge_why)
    else:
        # SRC/DST に 3 種 config の解決記録と属性を付ける (宛先と同じ config なので互換ハッシュが一致する)
        r3 = fsp.resolve_species(ref, forge, inplace=True)
        att = {"species_hash": r3["hash"], "species_record_sha256": r3["record"]["integrity"],
               "species_record_file": r3["record_file"], "species_input_unverified": 0}
        for fn in (RES, "nozzle.h5"):
            fsp.write_species_attrs(os.path.join(ref, fn), att)
        code, out, err = run_main([ref, new, "--Pt", "4.4e6", "--Ps", "1789.6", "--scale-ic", "pt"])
        check("(j) 作成: 終了コード 0", code == 0, err[-400:])
        if code == 0:
            rec = json.load(open(os.path.join(new, "RERUN_CONDITIONS.json")))
            check("(j) 記録にスケール検査 (roY2・roXi)", "roY2" in rec["scale_check_max_rel"] and "roXi" in rec["scale_check_max_rel"])
            with h5py.File(os.path.join(ref, RES), "r") as s, h5py.File(os.path.join(new, "nozzle.h5"), "r") as d:
                check("(j) 作成後の roY2・roXi が f 倍・ゼロはゼロ",
                      all(np.allclose(np.asarray(d["VALUE"][k], np.float64), 0.8 * np.asarray(s["VALUE"][k], np.float64),
                                      rtol=1e-6, atol=0) for k in ("roY2", "roXi")))
    done(ref, new)

    # ------------------------------------------------------------------ (k)
    cases = [
        ("inletProfile: 1", None, lambda t: t.replace("outputHDFflg: 0, ints: ,", "outputHDFflg: 0, ints: {inletProfile: 1},", 1),
         "inletProfile"),
        ("X{s} 形式", None, lambda t: t.replace("Y0: 0.91420000, Y1: 0.08580000", "X0: 0.86, X1: 0.14", 1), "X{s}"),
        ("inlet 2 本", None, lambda t: t.replace(
            "outlet:", "inlet2: {physID: 5, kind: inlet_Pressure,   outputHDFflg: 0, ints: , floats: {Y0: 0.9142, Y1: 0.0858, "
            "Pt: 5500000.0, Tt: 1600.0, k: 1.0, omega: 18000.0}}\noutlet:", 1), "inlet_Pressure が 2 本"),
        ("valueFileName が別名", lambda t: t.replace('valueFileName: "nozzle.h5"', 'valueFileName: "init.h5"', 1), None,
         "valueFileName"),
        ("外部参照 (speciesDBFile が run 外)", lambda t: t.replace("thermoHrefTemp: 298.15}",
                                                               'thermoHrefTemp: 298.15, speciesDBFile: "../db.yaml"}', 1), None,
         "外部参照"),
    ]
    for label, ec, eb, key in cases:
        ref, new = make_ref(edit_cfg=ec, edit_bc=eb)
        code, _, err = run_main([ref, new])
        check(f"(k) {label} → 作成前に拒否", code == 2 and key in err and not os.path.exists(new), err.strip()[-160:])
        done(ref, new)

    # ------------------------------------------------------------------ (l)
    ref, new = make_ref()
    code, _, err = run_main([ref, new, "--Pt", "4.4e6"])
    check("(l) --Pt のみ → 停止", code == 2 and not os.path.exists(new))
    pe = plan_of([ref, new, "--Pt", "4.4e6", "--keep-Ps"])["P_exit_ref"]
    check("(l) Ps/(f·P_exit_ref) を表示 (--keep-Ps なら 2237/(0.8·P_exit_ref))",
          f"Ps/(f·P_exit_ref) = {2237.0 / (0.8 * pe):.6g}" in err and f"P_exit_ref = {pe}" in err, err.strip()[-300:])
    code, _, _ = run_main([ref, new, "--Pt", "4.4e6", "--keep-Ps", "--dry-run"])
    check("(l) --keep-Ps なら通る (dry-run)", code == 0 and not os.path.exists(new))
    done(ref, new)

    # ------------------------------------------------------------------ (m)
    ref, new = make_ref()
    check("(m) 無変更 → recommended_stages none", plan_of([ref, new])["recommended_stages"]["stages"] == "none")
    check("(m) --steps/--cfl だけ → none (条件ではない)",
          plan_of([ref, new, "--steps", "12000", "--cfl", "3"])["recommended_stages"]["stages"] == "none")
    check("(m) --Tt → full", plan_of([ref, new, "--Tt", "1500"])["recommended_stages"]["stages"] == "full")
    check("(m) 参照と同じ値の --Pt は変更に数えない", plan_of([ref, new, "--Pt", "5.5e6"])["changes"] == {})
    p = plan_of([ref, new, "--cfl", "3"])
    check("(m) --cfl は `cfl: X, cfl_pseudo: X` を 1 回だけ書き換える",
          p["new_cfg_text"] == p["cfg_text"].replace("cfl: 5.0, cfl_pseudo: 5.0", "cfl: 3.0, cfl_pseudo: 3.0"))
    done(ref, new)

    # ------------------------------------------------------------------ 作成経路の周辺ロジック (restart_field を模擬)
    test_execute_with_mock_restart()

    # ------------------------------------------------------------------ run_staged_ns の段終了ゲート (§4.9)
    test_stage_gate()

    shutil.rmtree(TMP, ignore_errors=True)
    print(f"\nSKIP 件数: {len(skips)}")
    for n, w in skips:
        print(f"  - {n}: {w}")
    print(f"FAIL 件数: {fails}")
    return 1 if fails else 0


def test_execute_with_mock_restart():
    """**restart_field を模擬した**作成経路の試験 (複製の許可リスト・書き換えの書き出し・スケール・記録・prepare_info)。
    restart_field 自体の照合 (化学種) は通していないので、(a)(j) の作成経路の代わりにはならない — 周辺ロジックの回帰用。"""
    real_run = rc.subprocess.run

    def mock_run(cmd, *a, **k):
        if len(cmd) > 1 and str(cmd[1]).endswith("restart_field.py"):
            src, dst = cmd[2], cmd[3]
            with h5py.File(src, "r") as s, h5py.File(dst, "r+") as d:
                for n in d["VALUE"]:
                    if n != "wall_dist" and n in s["VALUE"]:
                        d["VALUE"][n][...] = np.asarray(s["VALUE"][n])
                for key in ("species_hash", "species_record_sha256", "species_record_file", "species_input_unverified"):
                    if key in s.attrs:
                        d.attrs[key] = s.attrs[key]
            return subprocess.CompletedProcess(cmd, 0, "VERDICT: OK (mock)\n", "")
        return real_run(cmd, *a, **k)

    rc.subprocess.run = mock_run
    try:
        ref, new = make_ref()
        code, out, err = run_main([ref, new, "--Pt", "4.4e6", "--Ps", "1789.6", "--scale-ic", "pt", "--steps", "12000"])
        check("[mock restart] 作成: 終了コード 0", code == 0, err[-300:])
        if code == 0:
            names = set(os.listdir(new))
            check("[mock restart] 許可リスト外 (VERDICT・PROVENANCE・res) を持ち込まない",
                  not ({"CONVERGENCE_VERDICT.txt", "RUN_PROVENANCE.txt", RES, "res_outlet_2_6000.h5"} & names), sorted(names))
            check("[mock restart] resolved_species_*・wall_*.csv・prepare_info を複製",
                  any(n.startswith("resolved_species_") for n in names) and any(n.startswith("wall_") for n in names)
                  and "prepare_info.json" in names)
            y = yaml.safe_load(open(os.path.join(new, "bcondConfig.yaml")))
            check("[mock restart] bcond を書き出した", y["inlet"]["floats"]["Pt"] == 4.4e6 and y["outlet"]["floats"]["Ps"] == 1789.6)
            check("[mock restart] solverConfig を書き出した", "nStepOuter: 12000" in open(os.path.join(new, "solverConfig.yaml")).read())
            rec = json.load(open(os.path.join(new, "RERUN_CONDITIONS.json")))
            with h5py.File(os.path.join(ref, RES), "r") as s, h5py.File(os.path.join(new, "nozzle.h5"), "r") as d:
                ok = all(np.allclose(np.asarray(d["VALUE"][k], np.float64), 0.8 * np.asarray(s["VALUE"][k], np.float64),
                                     rtol=1e-6, atol=0) for k in rec["required_conserved"])
            check("[mock restart] --scale-ic pt: 必要保存量が f = 0.8 倍", ok and abs(rec["f"] - 0.8) < 1e-15)
            check("[mock restart] 記録: 変更前後・P_exit_ref・Ps/(f·P_exit_ref)・recommended_stages・commit",
                  rec["changes"]["Pt"] == [5500000.0, 4400000.0] and abs(rec["P_exit_ref"] - 2242.0) < 4.0
                  and abs(rec["Ps_over_f_P_exit_ref"] - 1789.6 / (0.8 * rec["P_exit_ref"])) < 1e-12 and rec["recommended_stages"]["stages"] == "full"
                  and rec["tool_commit"]["head"] != "unknown" and rec["restart_field_verdict"].startswith("VERDICT: OK"))
            pi = json.load(open(os.path.join(new, "prepare_info.json")))
            pr = json.load(open(os.path.join(ref, "prepare_info.json")))
            check("[mock restart] prepare_info: 幾何据え置き・ic_from・rerun_of",
                  pi["x_E"] == pr["x_E"] and pi["ic_from"].endswith(RES) and pi.get("rerun_of"))
        done(ref, new)
    finally:
        rc.subprocess.run = real_run


def test_stage_gate():
    """forge を起動せず run_forge を差し替えて、soft 段の最終 res に NaN があれば次段へ進まないことを見る。"""
    from forge_design.evaluate import runner_axismach as ram
    rd = os.path.join(TMP, "staged")
    os.makedirs(rd)
    for fn in ("solverConfig.yaml", "bcondConfig.yaml"):
        shutil.copy2(os.path.join(BASE, fn), os.path.join(rd, fn))
    cfg_text = open(os.path.join(rd, "solverConfig.yaml")).read()

    nan_res = os.path.join(TMP, "nan_res.h5")
    shutil.copy2(os.path.join(BASE, RES), nan_res)
    with h5py.File(nan_res, "r+") as f:
        a = np.asarray(f["VALUE/roUx"]); a[123] = np.nan; f["VALUE/roUx"][...] = a
    check("段ゲート: NaN を入れた res は不合格", any("roUx" in s and "非有限" in s for s in ram.stage_gate(nan_res, cfg_text)))
    check("段ゲート: 参照 res は合格", ram.stage_gate(os.path.join(BASE, RES), cfg_text) == [])
    neg = os.path.join(TMP, "neg_res.h5")
    shutil.copy2(os.path.join(BASE, RES), neg)
    with h5py.File(neg, "r+") as f:
        a = np.asarray(f["VALUE/ro"]); a[5] = -1.0; f["VALUE/ro"][...] = a
    check("段ゲート: ρ ≤ 0 は不合格", any("ρ ≤ 0" in s for s in ram.stage_gate(neg, cfg_text)))
    check("convMethod: 2 も 1 次化 (convMethod: 12 は触らない)",
          ram._first_order("space: {convMethod: 2, limiter: 2}") == "space: {convMethod: 0, limiter: 2}"
          and ram._first_order("convMethod: 1,") == "convMethod: 0," and ram._first_order("convMethod: 12") == "convMethod: 12")

    calls = SimpleNamespace(forge=[], restart=[])

    def fake_forge_factory(src):
        def fake_forge(run_dir):
            run_dir = str(run_dir)
            t = open(os.path.join(run_dir, "solverConfig.yaml")).read()
            n = int(re.search(r"nStepOuter: (\d+)", t).group(1))
            calls.forge.append(t)
            shutil.copy2(src, os.path.join(run_dir, f"res_{n}.h5"))
            with open(os.path.join(run_dir, "residual_history.csv"), "w") as f:
                f.write("step,rms_ro\n0,1.0\n%d,0.5\n" % n)
            return 0
        return fake_forge

    orig = (ram.run_forge, ram._restart_same_mesh)
    try:
        ram._restart_same_mesh = lambda res, mesh: calls.restart.append(str(res))
        ram.run_forge = fake_forge_factory(nan_res)
        raised = None
        try:
            ram.run_staged_ns(rd, stages="full")
        except RuntimeError as e:
            raised = str(e)
        check("run_staged_ns: soft 段の res に NaN → RuntimeError で停止", raised is not None and "段終了ゲート" in (raised or ""),
              (raised or "")[:120])
        check("run_staged_ns: 次段 (mid) の forge を起動しない", len(calls.forge) == 1, len(calls.forge))
        check("run_staged_ns: NaN の場を restart_field に渡さない", calls.restart == [])
        man = json.load(open(os.path.join(rd, "stage_manifest.json")))
        check("run_staged_ns: 失敗した段も manifest と residual_history_S1_soft.csv に残る",
              [s["tag"] for s in man["stages"]] == ["S1_soft"] and os.path.exists(os.path.join(rd, "residual_history_S1_soft.csv")))
        check("run_staged_ns: soft 段は convMethod: 0", "convMethod: 0" in calls.forge[0] and "cfl: 0.5, cfl_pseudo: 0.5" in calls.forge[0])

        # 健全な場なら 3 段とも進み、manifest に 3 段・段ごとの履歴が残る (段の CFL・step 数は従来どおり)
        shutil.rmtree(rd)
        os.makedirs(rd)
        for fn in ("solverConfig.yaml", "bcondConfig.yaml"):
            shutil.copy2(os.path.join(BASE, fn), os.path.join(rd, fn))
        calls.forge.clear(); calls.restart.clear()
        ram.run_forge = fake_forge_factory(os.path.join(BASE, RES))
        r = ram.run_staged_ns(rd, stages="full")
        man = json.load(open(os.path.join(rd, "stage_manifest.json")))
        check("run_staged_ns (健全): 3 段・restart 2 回・rc 0", r == 0 and len(calls.forge) == 3 and len(calls.restart) == 2)
        check("run_staged_ns (健全): manifest の段 = S1_soft, S2_mid, main",
              [s["tag"] for s in man["stages"]] == ["S1_soft", "S2_mid", "main"])
        check("run_staged_ns (健全): 段ごとの残差履歴",
              all(os.path.exists(os.path.join(rd, f"residual_history_{t}.csv")) for t in ("S1_soft", "S2_mid", "main")))
        check("run_staged_ns (健全): 段の CFL・step 数は従来どおり (0.5/3000, 1.0/3000, 本段 5.0/6000)",
              "cfl: 0.5, cfl_pseudo: 0.5" in calls.forge[0] and "nStepOuter: 3000" in calls.forge[0]
              and "cfl: 1.0, cfl_pseudo: 1.0" in calls.forge[1] and "nStepOuter: 3000" in calls.forge[1]
              and calls.forge[2] == cfg_text)
        sys.path.insert(0, HERE)
        import stage_manifest as smod
        segs = smod.segments(man)
        check("run_staged_ns (健全): 判定区間 (最後の区間) は本段だけ", [s["tag"] for s in segs[-1]] == ["main"],
              [[s["tag"] for s in g] for g in segs])
    finally:
        ram.run_forge, ram._restart_same_mesh = orig


if __name__ == "__main__":
    try:
        sys.exit(main())
    finally:
        shutil.rmtree(TMP, ignore_errors=True)
