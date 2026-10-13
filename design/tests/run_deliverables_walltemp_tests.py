#!/usr/bin/env python3
"""成果物の壁温分布 (plan tooling-nozzle-wall-temperature-profile §5.1 #3 の成果物、§4.4・§6 V3・V1d の成果物側) の試験。

  V3  分布の run (case/70 c2b の問題の壁温を s の表にした prepare-only の run) で成果物を書き出し、
      `fluent/wall_temperature.prof` を**試験側の小さなパーサ**で読み直す。Fluent メッシュ (`read_fluent_ascii` で読み直したもの) の
      壁ゾーンの面中心に**最近傍で割り当てた値**が、run の JSON の実寸の表を面中心で `np.interp` した値と相対 1e-9 (点数 = 壁面数)。
      prepare_info の元の表 (s) を問題 YAML の長さ・r_t・メッシュの壁の終点で試験側が写した表とも相対 1e-9。
      境界条件表の壁の行が profile を指し、conditions.csv・xlsx (シート `wall_temperature`)・README・MANIFEST に表と記録が載る。
      壁ノードで書く誤実装なら最近傍の割り当てで面の値がずれる (試験の感度の確認)。表の範囲外の面中心は例外 (端値で延長しない)。
  V1d (成果物側の負例) CSV を 1 値書き換えた run・別の run の JSON と CSV を写した run・スロートの温度が同じで分布だけ違う表で作った
      run の JSON と差し替えた run で、`deliverables.export` が出力を作る前に例外で止まる (出力ディレクトリができない)。
      part を単体で呼んでも止まる。外の問題 YAML (一様 350 K) を渡しても壁温は run の JSON から書く。
  後方互換: 断熱 (M6 fixture) と一様 350 K (case/70 run_0010) の run で、変更前のコード (REF_COMMIT の design/forge_design と
      solver_density_cuda/tools を git archive) と今のコードの `part_fluent`・`part_excel` の出力 (境界条件表・conditions.csv・
      Fluent メッシュ・CFD 壁面量表ほか、xlsx 以外の全ファイル) がバイト同一、xlsx に渡す conditions・notes とシートの並びが同一。

分布の run は FORGE_WALLTEMP_RUNS (既定: このセッションの scratchpad の wt_e/test_out; `run_wall_temperature_tests.py` を
WALL_TEMP_TEST_OUT つきで回すと残る `v1d_run_prof_s`・`v1d_run_prof_s2` と問題 YAML) を読むだけで使う。無ければ一時ディレクトリに
`deltastar_loop --init-integral --prepare-only` で作る (約 1 分)。run には書かない (負例は小さいファイルだけを写した複製で行う)。
実行 (forge-deliv の直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_deliverables_walltemp_tests.py
"""
import copy
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import yaml
from scipy.spatial import cKDTree

DESIGN = Path(__file__).resolve().parents[1]
REPO = DESIGN.parent
BUILD = REPO / "solver_density_cuda" / ".build-native" / "relwithdebinfo"
os.environ.setdefault("FORGE_CONVERTER", str(BUILD / "convertGmshToForge"))
os.environ.setdefault("FORGE_BIN", str(BUILD / "forge"))
os.environ.setdefault("FORGE_NO_DESIGN_REPORT", "1")
sys.path.insert(0, str(DESIGN))

from forge_design import wall_temperature as WTP  # noqa: E402
from forge_design.export import deliverables as dl  # noqa: E402
from forge_design.export import fluent_msh as fm  # noqa: E402
from forge_design.export import mesh_structure as ms  # noqa: E402

SCRATCH = "/tmp/claude-1000/-home-sano-work-forge/bd18adac-85f4-4627-811b-8f3186867995/scratchpad"
RUNS = Path(os.environ.get("FORGE_WALLTEMP_RUNS", f"{SCRATCH}/wt_e/test_out"))
FIX = Path(os.environ.get("FORGE_DELIV_FIXTURE", f"{SCRATCH}/fixtures/m6_prod"))       # 断熱 (case/45 の生産の壁)
C70 = REPO / "case" / "70.m5_air_wt_demo" / "problem_m5_air_c2b.yaml"                   # 一様 350 K の問題
C70_RUN = REPO / "case" / "70.m5_air_wt_demo" / "run_0010_ns_c2b_final"                  # 一様 350 K の既存 run (読むだけ)
CASE45_PROBLEM = REPO / "case" / "45.isobutane_m6_d155" / "problem_d155_ns_prod.yaml"
REF_COMMIT = "5a6dd76f62d3d3103af1054d3cd70fd6b7974844"   # 成果物に壁温分布を入れる前 (配管だけの commit)
# run_wall_temperature_tests.py と同じ表 (スロート s = 0 の温度はどちらも 500 K)
PROF_S = {"s": [-2.0, -1.0, -0.5, -0.2, 0.0, 0.1, 0.4, 1.0], "Tw_K": [320.0, 340.0, 400.0, 470.0, 500.0, 480.0, 380.0, 330.0]}
PROF_S2 = {"s": [-2.0, -1.0, -0.3, 0.0, 0.3, 1.0], "Tw_K": [300.0, 310.0, 450.0, 500.0, 420.0, 350.0]}
SMALL = {"prepare_info.json", "bcondConfig.yaml", WTP.RUN_JSON, "wall_profile_3.csv"}   # 負例で書き換えるので実体を写す
fails = 0
SUMMARY = {}


def check(name, ok, info=""):
    global fails
    print(("ok   " if ok else "FAIL ") + name + (f" ({info})" if info != "" else ""), flush=True)
    fails += (not ok)


def skip(name, why):
    print(f"skip {name} ({why})", flush=True)


def raises(fn, *substr):
    """fn() が ValueError を出し、メッセージが substr を全部含むか。(結果, メッセージ)。"""
    try:
        fn()
    except ValueError as ex:
        msg = str(ex)
        return all(s in msg for s in substr), msg
    except Exception as ex:          # noqa: BLE001 (別の例外は不合格として内容を出す)
        return False, f"別の例外 {type(ex).__name__}: {str(ex)[:300]}"
    return False, "例外が出ない"


def sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ---------------------------------------------------------------------------------------------------------------------
# 分布の run の用意 (読むだけ / 無ければ作る) と複製
# ---------------------------------------------------------------------------------------------------------------------

def _prepare(problem: Path, run_dir: Path):
    env = dict(os.environ, FORGE_ROOT=str(REPO), FORGE_NO_DESIGN_REPORT="1")
    r = subprocess.run([sys.executable, "-m", "forge_design.feedback.deltastar_loop", "--problem", str(problem), "--euler-ref", "none",
                        "--run-dir", str(run_dir), "--init-integral", "--prepare-only"],
                       cwd=DESIGN, env=env, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def ensure_run(tag: str, prof: dict, tmp: Path):
    """(run, 問題 YAML, 由来)。RUNS の run が同じ表で照合を通れば読むだけで使い、無ければ tmp に prepare-only で作る。"""
    run, prob = RUNS / f"v1d_run_prof_{tag}", RUNS / f"c70_prof_{tag}.yaml"
    if run.is_dir() and prob.is_file():
        try:
            rec = WTP.load_run_profile(run)
        except ValueError:
            rec = None
        if rec is not None and rec["coord"] == "s" and rec["coords"] == prof["s"] and rec["Tw_K"] == prof["Tw_K"]:
            return run, prob, "既存 (読むだけ)"
    raw = yaml.safe_load(C70.read_text())
    raw["spec"]["wall_thermal"] = {"mode": "isothermal", "Tw_profile": copy.deepcopy(prof)}
    prob = tmp / f"c70_prof_{tag}.yaml"
    prob.write_text(yaml.safe_dump(raw, sort_keys=False, allow_unicode=True))
    run = tmp / f"run_prof_{tag}"
    rc, log = _prepare(prob, run)
    if rc != 0:
        raise RuntimeError(f"prepare-only が失敗 ({tag}): {log[-1500:]}")
    return run, prob, "作成 (prepare-only)"


def clone_run(src: Path, dst: Path) -> Path:
    """小さいファイル (正本・記録・bcond) は実体を写し、ほかはリンクする (成果物は run に書かない; 書き換えは実体だけ)。"""
    dst.mkdir(parents=True)
    for f in sorted(src.iterdir()):
        if f.is_file():
            if f.name in SMALL:
                shutil.copy2(f, dst / f.name)
            else:
                os.symlink(f.resolve(), dst / f.name)
    return dst


# ---------------------------------------------------------------------------------------------------------------------
# 試験側の読み込み (書き出し器のコードを使わない)
# ---------------------------------------------------------------------------------------------------------------------

def parse_point_profile(path: Path):
    """Fluent の point profile を S 式として読む: 戻り (1 行目, [名前, "point", N], {場の名前: ndarray})。"""
    txt = Path(path).read_text(encoding="ascii")
    toks = txt.replace("(", " ( ").replace(")", " ) ").split()

    def node(i):
        if toks[i] != "(":
            raise ValueError(f"字句 {i} が '(' でない: {toks[i]!r}")
        out, i = [], i + 1
        while toks[i] != ")":
            if toks[i] == "(":
                sub, i = node(i)
                out.append(sub)
            else:
                out.append(toks[i])
                i += 1
        return out, i + 1

    tree, end = node(0)
    if end != len(toks):
        raise ValueError(f"profile の後ろに余分な字句 ({len(toks) - end} 個)")
    head = tree[0]
    fields = {}
    for f in tree[1:]:
        if not isinstance(f, list) or not f or isinstance(f[0], list) or any(isinstance(v, list) for v in f[1:]):
            raise ValueError(f"場の書式が (名前 v1 ... vN) でない: {str(f)[:80]}")
        fields[f[0]] = np.array([float(v) for v in f[1:]])
    return txt.splitlines()[0], head, fields


def parse_md_table(path, header0):
    """`| header0 | ...` で始まる Markdown の表を [{列名: 値}] で返す。"""
    rows, hdr = [], None
    for line in Path(path).read_text().splitlines():
        if not line.startswith("|"):
            if hdr is not None and rows:
                break
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if hdr is None:
            if cells[0] == header0:
                hdr = cells
            continue
        if set(cells[0]) <= {"-", " "}:
            continue
        rows.append(dict(zip(hdr, cells)))
    return rows


def fluent_wall_face_centers(msh: Path, zone_name: str = "wall"):
    """Fluent メッシュを `read_fluent_ascii` で読み直し、壁ゾーン (名前 `zone_name`・型 wall) の面中心 (2 節点の中点) と節点番号。"""
    rd = fm.read_fluent_ascii(msh)
    ids = [zid for zid, (t, n) in rd["zone_names"].items() if n == zone_name and t == "wall"]
    if len(ids) != 1:
        raise ValueError(f"壁ゾーン {zone_name} が 1 つでない: {rd['zone_names']}")
    fz = [z for z in rd["face_zones"] if z["id"] == ids[0]]
    if len(fz) != 1 or fz[0]["bctype"] != 3:
        raise ValueError(f"壁ゾーン {zone_name} の 13 節が 1 つでない・型 3 でない")
    f = fz[0]["faces"]
    xy = rd["xy"][:, :2]
    return 0.5 * (xy[f[:, 0]] + xy[f[:, 1]]), f[:, :2], xy


def ref_x_from_s(coords_s, L_pipe_m, L_U_m, x_e_m):
    """試験側で書いた s → x [m] の写像 (直管 [−2,−1]・縮流部 [−1,0]・拡大部 [0,1])。"""
    xi, xc = -(L_pipe_m + L_U_m), -L_U_m
    out = []
    for s in coords_s:
        if s < -1.0:
            out.append(xi + (s + 2.0) * (xc - xi))
        elif s < 0.0:
            out.append(xc + (s + 1.0) * (0.0 - xc))
        else:
            out.append(s * x_e_m)
    return np.asarray(out)


# ---------------------------------------------------------------------------------------------------------------------
# 後方互換: 変更前と今のコードで part_fluent・part_excel を回す
# ---------------------------------------------------------------------------------------------------------------------

BC_SCRIPT = r'''
import json, sys
from pathlib import Path
design, run, out, problem = sys.argv[1:5]
sys.path.insert(0, design)
from forge_design.export import deliverables as dl
from forge_design.export import wall_tables as wt
cap = {}
real = wt.write_xlsx
def spy(*a, **k):
    cap["conditions"] = a[3] if len(a) > 3 else k.get("conditions")
    cap["notes"] = a[4] if len(a) > 4 else k.get("notes")
    cap["thermal_use"] = k.get("thermal_use")
    cap["wall_temperature"] = k.get("wall_temperature")
    r = real(*a, **k)
    cap["sheets"] = r.get("sheets")
    return r
wt.write_xlsx = spy
out = Path(out)
rf = dl.part_fluent(Path(run), out, None if problem == "-" else problem)
rx = dl.part_excel(Path(run), out)
cap["module"] = dl.__file__
cap["fluent_keys"] = sorted(rf)
cap["excel_keys"] = sorted(rx)
cap["cgns_stats"] = {k: rf["cgns"][k] for k in ("msh", "cgns", "match", "tags_preserved")}
(out / "xlsx_args.json").write_text(json.dumps(cap, ensure_ascii=False, default=str, indent=1))
'''
# thermal_use の判定ツールは一時ディレクトリ (乱数の名前) で回るので、VERDICT 行にその道が入る (コードの変更と無関係の揺れ)
_TMP_RE = re.compile(r"/tmp/thermal_use_(conv|qs)_[A-Za-z0-9_]+")


def _norm(b: bytes) -> bytes:
    return _TMP_RE.sub(lambda m: f"/tmp/thermal_use_{m.group(1)}_<TMP>", b.decode("utf-8", errors="surrogateescape")).encode(
        "utf-8", errors="surrogateescape")


def run_parts(design_dir: Path, run: Path, out: Path, problem):
    q = subprocess.run([sys.executable, "-I", "-c", BC_SCRIPT, str(design_dir), str(run), str(out), str(problem) if problem else "-"],
                       capture_output=True, text=True, cwd=str(out.parent))
    if q.returncode != 0:
        raise RuntimeError(f"part の書き出しが失敗 ({design_dir}): {q.stderr[-2000:]}")
    cap = json.loads((out / "xlsx_args.json").read_text())
    cap["notes"] = [_norm(str(t).encode()).decode() for t in (cap.get("notes") or [])]
    return cap


def compare_outputs(a: Path, b: Path):
    """2 つの出力ディレクトリの全ファイルを比べる。xlsx (作成時刻が入る)・CGNS (gmsh の HDF5 に時刻が入り、同じコードでも
    バイト同一にならない; 中身は gmsh の読み直しの統計で別に比べる)・記録 json は除く。thermal_use の一時ディレクトリの道は
    `<TMP>` に置き換えて比べる。戻り (比べた数, 違うファイル, 片方だけ, 置き換えたファイル)。"""
    def files(d):
        return {str(p.relative_to(d)) for p in d.rglob("*")
                if p.is_file() and p.suffix not in (".xlsx", ".cgns") and p.name != "xlsx_args.json"}
    fa, fb = files(a), files(b)
    both = sorted(fa & fb)
    diff, normed = [], []
    for k in both:
        ba, bb = (a / k).read_bytes(), (b / k).read_bytes()
        if ba != bb:
            if _norm(ba) == _norm(bb):
                normed.append(k)
            else:
                diff.append(k)
    return len(both), diff, sorted(fa ^ fb), normed


# =====================================================================================================================

def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="deliv_walltemp_"))
    keep = False
    try:
        return _main(tmp)
    except Exception:
        keep = True
        raise
    finally:
        if keep or fails:
            print(f"作業ディレクトリを残した: {tmp}")
        else:
            shutil.rmtree(tmp, ignore_errors=True)


def _main(tmp: Path) -> int:
    if not Path(os.environ["FORGE_CONVERTER"]).exists():
        skip("全体", f"変換器が無い ({os.environ['FORGE_CONVERTER']})")
        return 0
    runS, probS, howS = ensure_run("s", PROF_S, tmp)
    runS2, probS2, howS2 = ensure_run("s2", PROF_S2, tmp)
    print(f"分布の run: S = {runS} ({howS})、S2 = {runS2} ({howS2})", flush=True)
    recS = json.loads((runS / WTP.RUN_JSON).read_text())
    recS2 = json.loads((runS2 / WTP.RUN_JSON).read_text())
    piS = json.loads((runS / "prepare_info.json").read_text())

    # ================================================================ V3
    print("== V3 Fluent の壁温 profile・境界条件表・Excel・README (分布の run)", flush=True)
    run = clone_run(runS, tmp / "runS")
    out = tmp / "deliv_S"
    man = dl.export(run, out, skip=("cad", "mesh", "dimensions"), problem=probS)
    pf, px = man["parts"].get("fluent", {}), man["parts"].get("excel", {})
    check("書き出しが通る (fluent・excel が ok・errors なし)", pf.get("ok") is True and pf.get("pass") is True and px.get("ok") is True
          and not man["errors"], f"errors {man['errors']}・fluent pass {pf.get('pass')}")
    prof_path = out / "fluent" / "wall_temperature.prof"
    check("fluent/wall_temperature.prof がある", prof_path.is_file())
    if not prof_path.is_file():
        return 1
    first, head, fields = parse_point_profile(prof_path)
    print("     profile の先頭 10 行:")
    for ln in prof_path.read_text().splitlines()[:10]:
        print(f"       {ln}")
    cen, fnodes, xyF = fluent_wall_face_centers(out / "fluent" / "nozzle_2d_axi.msh")
    n_faces = len(cen)
    check(f"profile の 1 行目が ((wall-temperature point N)・場は x・y・temperature の 3 つ・各 N 点 (N = 壁面数 {n_faces})",
          first == f"((wall-temperature point {n_faces})" and head == ["wall-temperature", "point", str(n_faces)]
          and list(fields) == ["x", "y", "temperature"] and all(len(v) == n_faces for v in fields.values()),
          f"{first!r}・場 {list(fields)}・点 {[len(v) for v in fields.values()]}")
    pts = np.c_[fields["x"], fields["y"]]
    dist, idx = cKDTree(pts).query(cen)
    got = fields["temperature"][idx]
    want = np.interp(cen[:, 0], np.asarray(recS["x_m"], float), np.asarray(recS["Tw_x_K"], float))
    rel = np.abs(got - want) / np.abs(want)
    SUMMARY["V3"] = {"n_points": int(len(pts)), "n_wall_faces": int(n_faces), "max_rel": float(rel.max()),
                     "max_abs_K": float(np.abs(got - want).max()), "max_nn_dist_m": float(dist.max()),
                     "unique_assigned": int(len(np.unique(idx)))}
    check(f"V3: Fluent メッシュの壁面中心に最近傍で割り当てた値が run の JSON の実寸の表 (np.interp) と相対 1e-9 "
          f"(最大 {rel.max():.2e}・最大差 {np.abs(got - want).max():.2e} K・最近傍の距離の最大 {dist.max():.1e} m)",
          bool(rel.max() <= 1e-9) and len(np.unique(idx)) == n_faces, SUMMARY["V3"])
    # prepare_info の元の表 (s) を試験側で写す (区切り: 問題 YAML の L_pipe・L_U × r_t、壁の終点 = Fluent メッシュの壁の最大の x)
    g = yaml.safe_load(probS.read_text())["geometry"]
    scale = float(piS["scale_m"])
    wt0 = piS["wall_thermal"]["Tw_profile"]
    x_e_mesh = float(np.max(xyF[np.unique(fnodes.ravel()), 0]))
    x_ind = ref_x_from_s(wt0["coords"], float(g.get("L_pipe", 0.5)) * scale, float(g.get("L_U", 3.5)) * scale, x_e_mesh)
    want2 = np.interp(cen[:, 0], x_ind, np.asarray(wt0["Tw_K"], float))
    rel2 = np.abs(got - want2) / np.abs(want2)
    SUMMARY["V3"]["max_rel_vs_original_s_table"] = float(rel2.max())
    check(f"V3: prepare_info の元の表 (s) を試験側で写した表とも相対 1e-9 (最大 {rel2.max():.2e}; 壁の終点 メッシュ {x_e_mesh!r} m・"
          f"JSON {recS['anchors_m']['x_exit_m']!r} m)",
          bool(rel2.max() <= 1e-9) and wt0["coord"] == "s" and abs(x_e_mesh - recS["anchors_m"]["x_exit_m"]) <= 1e-9)
    # 感度: 壁ノードで書く誤実装なら最近傍の割り当てで面の値がずれる
    nodes = np.unique(fnodes.ravel())
    xn = xyF[nodes]
    tn = np.interp(xn[:, 0], np.asarray(recS["x_m"], float), np.asarray(recS["Tw_x_K"], float))
    _, idn = cKDTree(xn).query(cen)
    bad_abs = np.abs(tn[idn] - want)
    SUMMARY["V3"]["node_profile_max_abs_K"] = float(bad_abs.max())
    SUMMARY["V3"]["node_profile_faces_over_0p1K"] = int(np.sum(bad_abs > 0.1))
    check(f"V3 (感度): 壁ノードで書いた profile なら最大 {bad_abs.max():.3g} K ずれる ({int(np.sum(bad_abs > 0.1))}/{n_faces} 面が 0.1 K 超) "
          "— この試験は誤実装を落とす", bool(np.max(bad_abs / np.abs(want)) > 1e-6))
    rec_m = pf.get("wall_temperature_profile") or {}
    top = man.get("wall_temperature_profile") or {}
    check("MANIFEST の part 結果に profile の記録 (ヘッダ・点数・x 範囲・T 範囲・sha256) があり実ファイルと一致",
          rec_m.get("header") == first and rec_m.get("n_points") == n_faces and rec_m.get("sha256") == sha(prof_path)
          and man["files_sha256"].get("fluent/wall_temperature.prof") == sha(prof_path)
          and rec_m.get("x_range_m") == [float(fields["x"].min()), float(fields["x"].max())]
          and rec_m.get("T_range_K") == [float(fields["temperature"].min()), float(fields["temperature"].max())],
          {k: rec_m.get(k) for k in ("header", "n_points", "x_range_m", "T_range_K")})
    fi = piS["wall_thermal"]["files"]
    check("MANIFEST の最上位に正本 (JSON・CSV) の sha256 と profile の記録",
          top.get("run_json") == fi["json"] and top.get("run_json_sha256") == fi["json_sha256"]
          and top.get("run_csv_sha256") == fi["csv_sha256"] and (top.get("fluent_profile") or {}).get("sha256") == sha(prof_path), top)
    # 境界条件表
    rows = parse_md_table(out / "fluent" / "boundary_conditions.md", "ゾーン名")
    wall = [r for r in rows if r["ゾーン名"] == "wall"]
    check("境界条件表: 壁の行が profile `wall-temperature` `temperature` を指す (値の欄に点数と温度の範囲、出典に .prof と run の JSON)",
          len(wall) == 1 and "Thermal: Temperature = profile `wall-temperature` `temperature`" in wall[0]["項目 (Fluent の入力欄)"]
          and wall[0]["絶対値"].startswith(f"profile ({n_faces} 点, T_w ") and wall[0]["単位"] == "K"
          and "fluent/wall_temperature.prof" in wall[0]["出典"] and WTP.RUN_JSON in wall[0]["出典"]
          and "floats.Ts" not in wall[0]["出典"], wall)
    bc_txt = (out / "fluent" / "boundary_conditions.md").read_text()
    check("境界条件表: bcond の Ts は代表値として注記に回す", "代表値" in bc_txt and "wall_profile_3.csv" in bc_txt)
    # conditions.csv
    crow = {r[0]: r[1:] for r in csv.reader(open(out / "excel" / "conditions.csv"))}
    n_tab = len(recS["coords"])
    ok_pts = all(float(crow[f"wall_Tw_profile_{i}_orig_s"][0]) == recS["coords"][i]
                 and float(crow[f"wall_Tw_profile_{i}_orig_Tw_K"][0]) == recS["Tw_K"][i]
                 and float(crow[f"wall_Tw_profile_{i}_x_m"][0]) == recS["x_m"][i]
                 and float(crow[f"wall_Tw_profile_{i}_Tw_x_K"][0]) == recS["Tw_x_K"][i] for i in range(n_tab)) \
        if all(f"wall_Tw_profile_{i}_x_m" in crow for i in range(n_tab)) else False
    ok_anc = all(float(crow[f"wall_Tw_profile_anchor_{k}"][0]) == v for k, v in recS["anchors_m"].items()) \
        if all(f"wall_Tw_profile_anchor_{k}" in crow for k in recS["anchors_m"]) else False
    check(f"conditions.csv: 元の座標の種類 (s)・元の表・実寸の表 ({n_tab} 点)・区切りの点が JSON と一致 (3 列のまま)",
          crow.get("wall_Tw_profile_coord", [None])[0] == "s" and int(float(crow["wall_Tw_profile_n_points"][0])) == n_tab
          and ok_pts and ok_anc and crow.get("wall_Tw_profile_json_sha256", [None])[0] == fi["json_sha256"]
          and all(len(v) == 2 for v in crow.values()))
    xl = px.get("xlsx") or {}
    if xl.get("written"):
        with zipfile.ZipFile(out / "excel" / "nozzle_tables.xlsx") as z:
            names = re.findall(r'sheet name="([^"]+)"', z.read("xl/workbook.xml").decode())
            sx = z.read(f"xl/worksheets/sheet{names.index('wall_temperature') + 1}.xml").decode() if "wall_temperature" in names else ""
        # xlsxwriter は数値を 16 有効桁 (%.16g) で書くので、17 桁目が要る値は相対 1e-15 で照合する (表の正本は JSON と conditions.csv)
        vals = np.array([float(v) for v in re.findall(r"<v>([^<]+)</v>", sx) if re.fullmatch(r"[-+0-9.eE]+", v)])
        want_v = recS["x_m"] + recS["Tw_x_K"] + recS["coords"] + list(recS["anchors_m"].values())
        miss = [v for v in want_v if not np.any(np.abs(vals - v) <= 1e-15 * max(1.0, abs(v)))] if vals.size else want_v
        check(f"xlsx: シート {names} の最後に wall_temperature・元の表・実寸の表・区切りの点の値が載る (相対 1e-15)",
              names[-1:] == ["wall_temperature"] and not miss
              and (xl.get("info") or {}).get("wall_temperature", {}).get("n_rows") == n_tab, (miss, xl.get("info")))
        notes_ok = any("wall_temperature" in str(t) for t in (px.get("wall_temperature_profile") or {}).get("notes", []))
        check("xlsx の notes に壁温分布の行 (シートと正本の所在)", notes_ok, px.get("wall_temperature_profile"))
    else:
        skip("xlsx", f"書いていない: {xl}")
    rd_txt = (out / "README.md").read_text()
    sec = rd_txt.split("## 壁温分布")[1].split("\n## ")[0] if "## 壁温分布" in rd_txt else ""
    check("README: 壁温分布の節に Fluent での設定手順 (profile の読み込み・Thermal で wall-temperature temperature) と「読み込みは未確認」",
          "wall_temperature.prof" in sec and "`wall-temperature` `temperature`" in sec and "Fluent での読み込みは未確認" in sec
          and "手元に Fluent が無い" in sec)
    # 外の問題 YAML (一様 350 K) を渡しても壁温は run の JSON から書く
    out_u = tmp / "deliv_S_problem_uniform"
    man_u = dl.export(run, out_u, skip=("cad", "mesh", "excel", "dimensions"), problem=C70)
    rows_u = parse_md_table(out_u / "fluent" / "boundary_conditions.md", "ゾーン名")
    wu = [r for r in rows_u if r["ゾーン名"] == "wall"]
    check("外の問題 YAML (一様 350 K) を渡しても profile と壁の行は run の JSON のまま (sha256 同一)",
          sha(out_u / "fluent" / "wall_temperature.prof") == sha(prof_path) and len(wu) == 1
          and "profile `wall-temperature`" in wu[0]["項目 (Fluent の入力欄)"] and "350" not in wu[0]["絶対値"],
          (man_u["errors"], wu))
    # 表の範囲外の面中心は例外 (端値で延長しない)
    coords, quads, bedges, tags = ms.read_msh41_mesh2d(run / "nozzle.msh")
    ok, msg = raises(lambda: fm.wall_profile_points_2d(coords, bedges["wall"], recS["x_m"][:-1], recS["Tw_x_K"][:-1]), "表の範囲")
    check("表の範囲外の面中心があれば例外 (出口側の点を 1 つ落とした表)", ok, msg[:160])
    ok, msg = raises(lambda: fm.wall_profile_points_2d(coords, bedges["wall"], recS["x_m"][1:], recS["Tw_x_K"][1:]), "表の範囲")
    check("表の範囲外の面中心があれば例外 (入口側の点を 1 つ落とした表)", ok, msg[:160])

    # ================================================================ V1d
    print("== V1d 成果物側の正本の照合 (負例)", flush=True)
    check("(前提) S と S2 はスロートの温度が同じ (500 K) で分布が違う",
          recS["Tw_representative_K"] == recS2["Tw_representative_K"] == 500.0 and recS["Tw_x_K"] != recS2["Tw_x_K"])
    negs = []
    t1 = clone_run(runS, tmp / "neg_csv")
    lines = (t1 / "wall_profile_3.csv").read_text().splitlines()
    xv, tv = lines[3].split()
    lines[3] = f"{xv} {float(tv) + 1.0!r}"
    (t1 / "wall_profile_3.csv").write_text("\n".join(lines) + "\n")
    negs.append(("CSV を 1 値書き換えた (+1 K)", t1, "CSV が書き換えられた"))
    t2 = clone_run(runS, tmp / "neg_other_run")
    shutil.copy2(runS2 / WTP.RUN_JSON, t2 / WTP.RUN_JSON)
    shutil.copy2(runS2 / "wall_profile_3.csv", t2 / "wall_profile_3.csv")
    negs.append(("別の run (S2) の JSON と CSV を写した", t2, "prepare_info と合わない"))
    t3 = clone_run(runS, tmp / "neg_same_throat_json")
    shutil.copy2(runS2 / WTP.RUN_JSON, t3 / WTP.RUN_JSON)
    negs.append(("スロートの温度が同じで分布だけ違う表で作った run の JSON と差し替えた", t3, "書き換えられた"))
    msgs = {}
    for name, rdir, key in negs:
        o = tmp / f"out_{rdir.name}"
        ok, msg = raises(lambda: dl.export(rdir, o, skip=("cad", "mesh", "dimensions"), problem=probS), "壁温分布の正本の照合で不合格", key)
        msgs[name] = msg
        check(f"V1d: {name} → export が例外で止まり、出力ディレクトリを作らない", ok and not o.exists(), msg[:200])
        okp, msgp = raises(lambda: dl.part_fluent(rdir, tmp / f"pf_{rdir.name}", probS), key)
        okx, msgx = raises(lambda: dl.part_excel(rdir, tmp / f"px_{rdir.name}"), key)
        check(f"V1d: {name} → part_fluent・part_excel を単体で呼んでも止まる (Fluent の profile を書かない)",
              okp and okx and not (tmp / f"pf_{rdir.name}" / "fluent" / "wall_temperature.prof").exists(), (msgp[:120], msgx[:120]))
    SUMMARY["V1d"] = msgs

    # ================================================================ 後方互換
    print("== 後方互換 (断熱・一様は変更前とバイト同一)", flush=True)
    ref_root = tmp / "ref_code"
    ref_root.mkdir()
    q = subprocess.run(f"git -C {REPO} archive {REF_COMMIT} design/forge_design solver_density_cuda/tools | tar -x -C {ref_root}",
                       shell=True, capture_output=True, text=True)
    if q.returncode != 0:
        skip("後方互換", f"変更前のコードを取り出せない ({REF_COMMIT[:8]}): {q.stderr[-300:]}")
    else:
        cases = []
        if (FIX / "nozzle.msh").exists():
            cases.append(("断熱 (M6 fixture)", FIX, CASE45_PROBLEM if CASE45_PROBLEM.is_file() else None))
        else:
            skip("後方互換 断熱", f"fixture が無い ({FIX})")
        if (C70_RUN / "nozzle.msh").exists():
            cases.append(("一様 350 K (case/70 run_0010)", C70_RUN, C70))
        else:
            skip("後方互換 一様", f"run が無い ({C70_RUN})")
        bc = {}
        for i, (label, rdir, prob) in enumerate(cases):
            check(f"(前提) {label}: load_run_profile は None (分布でない)", WTP.load_run_profile(rdir) is None)
            # 同じ出力先の名前で順に回す (status.json は出力先の絶対パスを持つので、名前を揃えないと中身が違って見える)
            work, oo, on = tmp / f"bc_work_{i}", tmp / f"bc_old_{i}", tmp / f"bc_new_{i}"
            work.mkdir()
            co = run_parts(ref_root / "design", rdir, work, prob)
            work.rename(oo)
            work.mkdir()
            cn = run_parts(DESIGN, rdir, work, prob)
            work.rename(on)
            n, diff, only, normed = compare_outputs(oo, on)
            check(f"{label}: 変更前のコードが変更前の場所から読まれた", str(ref_root) in co["module"] and str(DESIGN) in cn["module"],
                  (co["module"], cn["module"]))
            check(f"{label}: part_fluent・part_excel の出力 ({n} ファイル; 境界条件表・conditions.csv・Fluent メッシュほか) が変更前とバイト同一"
                  + (f" (一時ディレクトリの道だけ置き換えて比べた: {normed})" if normed else ""),
                  n >= 4 and not diff and not only and (on / "fluent" / "boundary_conditions.md").exists()
                  and (on / "excel" / "conditions.csv").exists()
                  and (oo / "fluent" / "boundary_conditions.md").read_bytes() == (on / "fluent" / "boundary_conditions.md").read_bytes()
                  and (oo / "excel" / "conditions.csv").read_bytes() == (on / "excel" / "conditions.csv").read_bytes(),
                  f"違う {diff}・片方だけ {only}")
            check(f"{label}: CGNS の中身 (gmsh で読み直した節点数・要素数・物理名) が変更前と同一", co["cgns_stats"] == cn["cgns_stats"])
            check(f"{label}: xlsx に渡す conditions・notes・シートの並びが変更前と同一 (壁温のシートなし)",
                  co["conditions"] == cn["conditions"] and co["notes"] == cn["notes"] and co["sheets"] == cn["sheets"]
                  and cn["wall_temperature"] is None and "wall_temperature" not in (cn["sheets"] or []),
                  f"シート {cn['sheets']}")
            check(f"{label}: part の結果のキーが変更前と同じ (壁温の記録を足さない)",
                  co["fluent_keys"] == cn["fluent_keys"] and co["excel_keys"] == cn["excel_keys"],
                  (sorted(set(cn["fluent_keys"]) ^ set(co["fluent_keys"])), sorted(set(cn["excel_keys"]) ^ set(co["excel_keys"]))))
            bc[label] = {"files_compared": n, "diff": diff, "normalized_tmp_paths": normed, "excluded": "xlsx (作成時刻)・cgns (HDF5 の時刻)"}
        SUMMARY["backcompat"] = bc
    print("SUMMARY " + json.dumps(SUMMARY, ensure_ascii=False, default=str))
    print(f"FAIL 件数: {fails}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
