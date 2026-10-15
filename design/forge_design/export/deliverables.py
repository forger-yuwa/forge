r"""風洞ノズルの成果物一式を最終 NS の run から書き出す (plans/accepted/tooling-nozzle-wt-deliverables.md §4.2)。

`<run>/deliverables/` に次を作り、各成果物に独立の読み直し検査の結果を付ける:

- `cad/`    2D の壁曲線の STEP (`wall_step.export_run`)・3D の内面の STEP と STL (`wall_step.write_surface`)・gmsh での読み直し
- `mesh/`   解析メッシュ (msh4.1 のまま・msh2.2)・構造ファイル・STEP から同じ構造で作り直す gmsh スクリプトとその照合
- `fluent/` 2D 軸対称の Fluent メッシュ (`.msh`)・CGNS・境界条件表 (壁温分布の run は壁温の point profile `wall_temperature.prof` も)
- `excel/`  形状表・CFD 壁面量表 (状態つき)・条件・`.xlsx` (壁温分布の run は壁温の表も)
- `DIMENSIONS.json` 寸法の出荷ゲート (§4.2a)、`MANIFEST.json` (sha256・検査の数値・飛ばしたものと理由)、`README.md`

入力が欠けた成果物は飛ばして理由を MANIFEST に書く (黙って欠けさせない)。寸法ゲートが FAIL なら MANIFEST の `verdict` を FAIL にし、
README の冒頭に書く (書き出しは止めない — 何が外れたかを見せるため)。壁温分布の run は、run の正本 (`wall_temperature_profile.json`・
`wall_profile_<id>.csv`) が `prepare_info.json` の sha256 と合わなければ何も書かずに止める (plan tooling-nozzle-wall-temperature-profile §4.2)。

usage:
  design/.venv-opt/bin/python -m forge_design.export.deliverables RUN [--card CARD] [--problem YAML] [--out DIR] [--step N] [--skip cad,fluent,...]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FORGE_ROOT = HERE.parents[2]
GMSH_PYTHON = "/home/sano/work/forge/.venv-mesh/bin/python"
PPTX_PYTHON = "/home/sano/work/forge/.venv-pptx/bin/python"
PARTS = ("cad", "mesh", "fluent", "excel", "dimensions")


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _jsonable(o):
    return json.loads(json.dumps(o, default=lambda v: v.tolist() if isinstance(v, np.ndarray) else float(v)))


class _Log:
    def __init__(self):
        self.parts, self.skipped, self.errors = {}, {}, {}

    def run(self, name, fn):
        t0 = time.time()
        try:
            r = fn()
            self.parts[name] = {"ok": True, "seconds": round(time.time() - t0, 2), **(r or {})}
        except _Skip as e:
            self.skipped[name] = str(e)
        except Exception as e:  # noqa: BLE001 — 失敗は MANIFEST に残して他の成果物は続ける
            import traceback
            self.errors[name] = f"{type(e).__name__}: {e}"
            self.parts[name] = {"ok": False, "seconds": round(time.time() - t0, 2), "traceback": traceback.format_exc()[-3000:]}


class _Skip(Exception):
    pass


def _need(run: Path, names):
    miss = [n for n in names if not (run / n).exists()]
    if miss:
        raise _Skip(f"入力が無い: {miss}")


def _gmsh(script: str, *args, timeout=1800) -> dict:
    """gmsh の python (.venv-mesh) で短いスクリプトを走らせ、最後の行の JSON を返す。"""
    q = subprocess.run([GMSH_PYTHON, "-c", script, *map(str, args)], capture_output=True, text=True, timeout=timeout)
    if q.returncode != 0:
        raise RuntimeError(f"gmsh の作業が失敗 (rc={q.returncode}):\n{q.stdout[-1500:]}\n{q.stderr[-1500:]}")
    return json.loads(q.stdout.strip().splitlines()[-1])


_MSH22 = r"""
import sys, json, gmsh
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
gmsh.open(sys.argv[1]); gmsh.option.setNumber("Mesh.MshFileVersion", 2.2); gmsh.option.setNumber("Mesh.Binary", 0)
gmsh.write(sys.argv[2]); gmsh.clear(); gmsh.open(sys.argv[2])
nt, _, _ = gmsh.model.mesh.getNodes()
names = {gmsh.model.getPhysicalName(d, t): d for d, t in gmsh.model.getPhysicalGroups()}
print(json.dumps({"n_nodes": len(nt), "physical": names})); gmsh.finalize()
"""


def part_cad(run: Path, out: Path, freecadcmd=None) -> dict:
    """cad/: 2D の STEP (転送誤差の検査つき)、3D の内面の STEP と STL、独立の読み直し (§6 V2)。
    合格は**全部の検査が実行されて合格したときだけ** (失敗・未実行・結果の欠落を合格にしない; codex result レビュー 2026-10-13 M2):
    2D の転送 PASS・回転面が作れる、3D の書き出しの FreeCAD の読み直しが有効で面 1 枚・ノット区間の和の面積が ∫2πr ds と相対 ≤ 1e-6、
    gmsh の読み直しが成功・求積が収束・面積が相対 ≤ 1e-6・bbox の x 範囲が曲線の定義域と ≤ 1e-6 mm (+ OCC の許容差 1e-7 mm)。"""
    from . import wall_step
    _need(run, ["wall_repr.json"])
    rec = json.loads((run / "wall_repr.json").read_text())
    if rec.get("physical_wall_repr") != "single_bspline":
        raise _Skip(f"物理壁の表現が single_bspline でない ({rec.get('physical_wall_repr')!r}) — STEP は最終 NS 用の雛形 (axismach_final_v1) の run から")
    d = out / "cad"
    sc = wall_step.export_run(run, d, freecadcmd)
    tr = sc["readback"]["transfer"]
    res = {"step2d": "cad/wall_physical.step", "transfer_pass": bool(tr["pass"]), "transfer_pos_max_mm": tr.get("pos_max_mm"),
           "revolve_ok": bool(sc["readback"]["revolve_ok"]), "chord_minus_curve": sc["cad_vs_cfd"].get("chord_minus_curve")}
    checks = {"transfer_2d": res["transfer_pass"], "revolve_2d": res["revolve_ok"]}
    # 解析の面積 ∫2πr ds (保存した曲線を 40 万区間の台形; 打ち切り誤差 ≲1e-10 相対) と曲線の定義域 [mm]
    from ..geometry.wall_axismach import load_wall_file
    W = load_wall_file(run)
    s = float(W["scale_m"]) * 1000.0
    x = np.linspace(*W["domain"], 400001)
    r = np.asarray(W["physical"].r(x), dtype=float)
    A = float(np.sum(np.pi * (r[1:] + r[:-1]) * np.hypot(np.diff(x), np.diff(r)))) * s * s
    dom = [float(W["domain"][0]) * s, float(W["domain"][1]) * s]
    res["surface_area_analytic_mm2"], res["domain_mm"] = A, dom
    sf = wall_step.write_surface(d / "wall_physical.step", d / "wall_surface.step", d / "wall_surface.stl")
    res["surface"] = sf
    fc_area = sf.get("readback_area_knot_split_mm2")
    fc_bb = sf.get("readback_bbox_mm") or [np.nan] * 6
    checks["freecad_readback_valid"] = sf.get("readback_valid") is True
    checks["freecad_one_face"] = sf.get("readback_faces") == 1
    checks["freecad_area"] = fc_area is not None and abs(fc_area - A) / A <= 1e-6
    checks["freecad_bbox_x"] = bool(abs(fc_bb[0] - dom[0]) <= 1e-6 and abs(fc_bb[3] - dom[1]) <= 1e-6)
    chk = HERE / "cad_check_gmsh.py"
    q = subprocess.run([GMSH_PYTHON, str(chk), str(d / "wall_surface.step")], capture_output=True, text=True, timeout=600)
    g = None
    if q.returncode == 0 and q.stdout.strip():
        try:
            g = json.loads(q.stdout.strip().splitlines()[-1])
        except ValueError:
            g = None
    res["surface_gmsh"] = g if g is not None else {"error": (q.stdout + q.stderr)[-800:], "exit": q.returncode}
    g = g or {}
    g_bb = g.get("bbox_mm") or [np.nan] * 6
    checks["gmsh_readback"] = q.returncode == 0 and "area_mm2" in g
    checks["gmsh_quadrature_converged"] = g.get("area_quadrature_converged") is True
    checks["gmsh_area"] = "area_mm2" in g and abs(float(g["area_mm2"]) - A) / A <= 1e-6
    checks["gmsh_bbox_x"] = bool(abs(g_bb[0] - dom[0]) <= 1.1e-6 and abs(g_bb[3] - dom[1]) <= 1.1e-6)
    res["surface_area_rel_err"] = {"freecad": (abs(fc_area - A) / A if fc_area is not None else None),
                                   "gmsh": (abs(float(g["area_mm2"]) - A) / A if "area_mm2" in g else None)}
    res["checks"] = checks
    res["pass"] = bool(all(checks.values()))
    return res

def _converter() -> Path | None:
    """convertGmshToForge の場所 (runner と同じ規則: FORGE_CONVERTER > FORGE_BIN の隣 > 既定のビルド)。この作業ツリーの native ビルドも探す。"""
    from ..evaluate.runner import converter_path
    for c in (converter_path(), FORGE_ROOT / "solver_density_cuda" / ".build-native" / "relwithdebinfo" / "convertGmshToForge"):
        if Path(c).exists():
            return Path(c)
    return None


def _mesh_quality(msh: Path, run: Path, ar_max: float) -> dict:
    """primal の品質 (cell 変換のコピーで測る; runner_axismach の prepare と同じ作り方) を check_mesh_quality.py で判定する。"""
    import tempfile
    conv = _converter()
    if conv is None:
        return {"verdict": "未判定", "reason": "convertGmshToForge が見つからない (FORGE_CONVERTER で指定)"}
    if not (run / "solverConfig.yaml").exists() or not (run / "bcondConfig.yaml").exists():
        return {"verdict": "未判定", "reason": "run に solverConfig.yaml / bcondConfig.yaml が無い (変換の設定を作れない)"}
    cfg = (run / "solverConfig.yaml").read_text()
    cfg = cfg.replace('discretization: "node"', 'discretization: "cell"').replace(", axisCentroidShift: 1, nodeWallDirichlet: 1", "")
    import re
    cfg = re.sub(r'meshFileName: "[^"]*"', 'meshFileName: "nozzle_qc.h5"', cfg)
    cfg = re.sub(r'valueFileName: "[^"]*"', 'valueFileName: "nozzle_qc.h5"', cfg)
    with tempfile.TemporaryDirectory(prefix="deliv_qc_") as td:
        t = Path(td)
        (t / "solverConfig.yaml").write_text(cfg)
        shutil.copy2(run / "bcondConfig.yaml", t / "bcondConfig.yaml")
        shutil.copy2(msh, t / "nozzle.msh")
        c = subprocess.run([str(conv), "nozzle.msh", "nozzle_qc.h5"], cwd=t, capture_output=True, text=True, timeout=3600)
        if c.returncode != 0:
            return {"verdict": "FAIL", "reason": "変換器が失敗", "log": (c.stdout + c.stderr)[-800:]}
        q = subprocess.run([sys.executable, str(FORGE_ROOT / "solver_density_cuda" / "tools" / "check_mesh_quality.py"),
                            "nozzle_qc.h5", "--ar-max", str(int(ar_max))], cwd=t, capture_output=True, text=True, timeout=1800)
    lines = [ln for ln in q.stdout.splitlines() if ln.strip()]
    v = next((ln for ln in reversed(lines) if ln.startswith("VERDICT")), "VERDICT: ?")
    return {"verdict": "PASS" if v.startswith("VERDICT: PASS") else "FAIL", "verdict_line": v, "ar_max": ar_max,
            "summary": [ln for ln in lines if ("AR" in ln or "skew" in ln)][:6], "converter": str(conv)}


def _wall_slope(run: Path, coords, bedges) -> tuple:
    """壁法線の近壁層の格子 (構造ファイル v2) に渡す station ごとの壁の傾き r′ = dr/dx (plan tooling-nozzle-card-local-mesh §5.1 #5)。
    run の壁ファイル (`load_wall_file(run)["physical"]` の r(x, 1)) を壁の節点の x [m] / scale_m (r_t 単位) で評価する。
    壁の節点の x が定義域の外に丸めの分 (定義域の幅の 1e-9 以内) だけ出ていれば端で評価する。戻り値は (傾き, 出典 dict, None)、
    作れないときは (None, None, 理由)。縦線の格子では extract_structure が傾きを使わない (v1 の出力は渡しても同じ)。"""
    from ..geometry.wall_axismach import load_wall_file
    try:
        W = load_wall_file(run)
        ni, nj = len(bedges["wall"]) + 1, len(bedges["inlet"]) + 1
        if coords.shape[0] != ni * nj:
            raise ValueError(f"節点数 {coords.shape[0]} が ni·nj = {ni * nj} でない")
        scale = float(W["scale_m"])
        x = np.asarray(coords[:, 0], dtype=float).reshape(ni, nj)[:, -1] / scale
        lo, hi = (float(v) for v in W["domain"])
        tol = 1e-9 * (hi - lo)
        if x.min() < lo - tol or x.max() > hi + tol:
            raise ValueError(f"壁の節点の x [{x.min():.12g}, {x.max():.12g}] r_t が壁ファイルの定義域 [{lo:.12g}, {hi:.12g}] の外")
        xc = np.clip(x, lo, hi)
        rp = np.asarray(W["physical"].r(xc, 1), dtype=float)
        if rp.shape != x.shape or not np.all(np.isfinite(rp)):
            raise ValueError("壁の傾きが非有限か形が違う")
    except Exception as e:  # noqa: BLE001 — 作れなければ理由を返す (縦線の格子は傾きなしで v1 を作る)
        return None, None, f"{type(e).__name__}: {e}"
    src = {"file": Path(W["path"]).name, "sha256": _sha(Path(W["path"])), "kind": W["kind"], "scale_m": scale,
           "how": "load_wall_file(run)['physical'].r(x, 1)、x = 壁の節点の x [m] / scale_m (r_t 単位)",
           "n_clipped_to_domain": int(np.sum(xc != x))}
    return rp, src, None


def _rule_params(run: Path, msh_sha256: str) -> tuple:
    """構造ファイル v2 の規則の照合に渡す mesh2d の規則のパラメータ (plan tooling-nozzle-card-local-mesh §5.1 #5 の判断)。
    run の prepare_info.json の `mesh.params` (生成に使った Mesh2DParams そのもの; `mesh.block` は問題 YAML の mesh 節で、
    wall_normal_end_stations を書かないことがある) の `wall_normal_layer`・`wall_normal_end_stations` (キーが無ければ null =
    Mesh2DParams の既定) を取る。出典に prepare_info.json の sha256 と、そこに記録された msh の sha256 が run の nozzle.msh と
    一致するかを書く。戻り値は (rule_params, None)、取れないときは (None, 理由)。"""
    p = run / "prepare_info.json"
    try:
        info = json.loads(p.read_text())
        mesh = info["mesh"]
        prm = mesh["params"]
        wnl = prm.get("wall_normal_layer")
        if wnl is None:
            raise ValueError("prepare_info.json の mesh.params に wall_normal_layer が無い (null)")
    except Exception as e:  # noqa: BLE001 — 取れなければ理由を返す (規則の照合は未実施になる)
        return None, f"{type(e).__name__}: {e}"
    rec = (mesh.get("hashes") or {}).get("msh_sha256")
    src = {"file": p.name, "sha256": _sha(p), "key": "mesh.params",
           "msh_sha256_recorded": rec, "msh_sha256_matches": None if rec is None else rec == msh_sha256}
    return {"wall_normal_layer": list(wnl), "wall_normal_end_stations": prm.get("wall_normal_end_stations"), "source": src}, None


def part_mesh(run: Path, out: Path) -> dict:
    from . import mesh_structure as ms
    _need(run, ["nozzle.msh"])
    d = out / "mesh"
    d.mkdir(parents=True, exist_ok=True)
    shutil.copy2(run / "nozzle.msh", d / "nozzle.msh")
    res = {"msh41": "mesh/nozzle.msh", "msh41_sha256": _sha(d / "nozzle.msh")}
    res["msh22"] = _gmsh(_MSH22, d / "nozzle.msh", d / "nozzle_v22.msh")
    coords, quads, bedges, tags = ms.read_msh41_mesh2d(d / "nozzle.msh")
    # 壁の傾きは壁法線の近壁層の格子 (縦線でない) の構造ファイル v2 に要る。縦線の格子は渡しても v1 のまま (出力は今までと同じ)
    slope, slope_src, slope_why = _wall_slope(run, coords, bedges)
    # 規則のパラメータ (v2 の規則の照合; 取れなければ照合は未実施と記録される)
    rule, rule_why = _rule_params(run, res["msh41_sha256"])
    try:
        st = ms.extract_structure(coords, quads, bedges, wall_slope=slope, wall_slope_source=slope_src, rule_params=rule)
    except Exception as e:  # noqa: BLE001 — 規則に従わない格子・壁の傾きを作れない非縦線の格子は構造ファイルを作らず理由を書く
        res["structure"] = f"作らない: {type(e).__name__}: {e}" + (f" (壁の傾きを作れない: {slope_why})" if slope_why else "")
        return res
    ms.write_structure_json(st, d / "mesh_structure.json")
    shutil.copy2(HERE / "rebuild_from_step.py", d / "rebuild_from_step.py")
    res["structure"] = "mesh/mesh_structure.json"
    res["structure_format"] = st["format"]
    if st["format"] == ms.FORMAT_V2:                 # 壁法線の近壁層の格子: 節点の規則・自己検査・規則の照合・出典
        res["structure_node_rule"] = st["node_rule"]
        res["structure_self_check"] = st["checks"]["self_check"]
        res["structure_rule_check"] = st["checks"]["rule_check"]
        res["structure_wall_slope_source"] = st["wall_slope_source"]
        res["structure_rule_params"] = st["rule_params"]
        if rule is None:
            res["structure_rule_params_missing"] = rule_why
    step = out / "cad" / "wall_physical.step"
    if step.exists():
        q = subprocess.run([GMSH_PYTHON, str(d / "rebuild_from_step.py"), str(step), str(d / "mesh_structure.json"),
                            str(d / "rebuilt_from_step"), "--compare", str(d / "nozzle.msh")],
                           capture_output=True, text=True, timeout=3600)
        cj = d / "rebuilt_from_step_compare.json"
        res["rebuild"] = json.loads(cj.read_text()) if cj.exists() else {"error": (q.stdout + q.stderr)[-1500:]}
        # V3(d): 元と作り直し (forge 形の msh) の両方が明示した閾値で品質 PASS
        ar_max = 1000.0
        try:
            ar_max = float(json.loads((run / "prepare_info.json").read_text())["mesh"]["block"].get("ar_max", 1000))
        except Exception:  # noqa: BLE001
            pass
        res["quality"] = {"original": _mesh_quality(d / "nozzle.msh", run, ar_max)}
        rb = d / "rebuilt_from_step.msh"
        if rb.exists():
            res["quality"]["rebuilt"] = _mesh_quality(rb, run, ar_max)
        qv = [v.get("verdict") for v in res["quality"].values()]
        res["pass"] = bool(isinstance(res["rebuild"], dict) and res["rebuild"].get("verdict") == "PASS"
                           and all(v == "PASS" for v in qv) and len(qv) == 2)
    else:
        res["rebuild"] = "STEP が無いので作り直しの照合は飛ばした"
    return res


def _wall_profile(run: Path):
    """run の壁温分布の正本 (`wall_temperature_profile.json`・`wall_profile_<id>.csv`) を `prepare_info.json` の sha256 と照合して読む
    (plan tooling-nozzle-wall-temperature-profile §4.2・V1d)。壁温は run の JSON から読み、外の問題 YAML や bcond の代表 `Ts` で
    代替しない。一様・断熱の run は None。不一致 (書き換え・別の run の写し・bcond の wallProfile なし) は ValueError。
    戻り: {"rec": JSON の中身, "files": prepare_info の wall_thermal.files (名前と sha256)}。"""
    from ..wall_temperature import load_run_profile
    try:
        rec = load_run_profile(run)
    except ValueError as ex:
        raise ValueError(f"成果物: 壁温分布の正本の照合で不合格 — 書き出さない: {ex}") from ex
    if rec is None:
        return None
    files = json.loads((run / "prepare_info.json").read_text())["wall_thermal"]["files"]
    return {"rec": rec, "files": files}


def part_fluent(run: Path, out: Path, problem=None) -> dict:
    """fluent/: 2D 軸対称の Fluent メッシュ・CGNS・境界条件表。壁温分布の run は、壁ゾーンの面中心で run の正本の表を評価した
    point profile `wall_temperature.prof` も書き、境界条件表の壁の行をその profile にする (plan tooling-nozzle-wall-temperature-profile §4.4)。
    正本の照合と、面中心が表の範囲内かの検査は何かを書く前に行う (不合格なら fluent/ に何も書かない)。"""
    from . import fluent_msh as fm
    from . import mesh_structure as ms
    _need(run, ["nozzle.msh", "bcondConfig.yaml", "prepare_info.json"])
    wtp = _wall_profile(run)
    coords, quads, bedges, tags = ms.read_msh41_mesh2d(run / "nozzle.msh")
    wall_zone = prof_pts = None
    if wtp is not None:
        wid = int(wtp["rec"]["wall_id"])
        names = [n for n, t in tags.items() if int(t) == wid]
        if len(names) != 1 or names[0] not in bedges:
            raise ValueError(f"壁温分布の壁 (physID {wid}) が msh の境界に 1 つでない (物理名 {names}、境界 {sorted(bedges)})")
        wall_zone = names[0]
        prof_pts = fm.wall_profile_points_2d(coords, bedges[wall_zone], wtp["rec"]["x_m"], wtp["rec"]["Tw_x_K"])
    d = out / "fluent"
    d.mkdir(parents=True, exist_ok=True)
    fp = d / "nozzle_2d_axi.msh"
    zs = fm.write_fluent_msh_2d(fp, coords, quads, bedges)
    chk = fm.check_fluent_msh_2d(fp, coords, quads, bedges)
    cg = fm.write_cgns_via_gmsh(run / "nozzle.msh", d / "nozzle_2d_axi.cgns", GMSH_PYTHON)
    prob = _problem_dict(problem)
    wp = None
    if wtp is not None:
        pr = fm.write_point_profile_2d(d / "wall_temperature.prof", *prof_pts)
        files = wtp["files"]
        n_faces = int(zs["zones"][wall_zone]["faces"])
        wp = {"file": "fluent/wall_temperature.prof", **pr, "wall_zone": wall_zone, "wall_id": int(wtp["rec"]["wall_id"]),
              "n_wall_faces": n_faces, "points_match_wall_faces": bool(pr["n_points"] == n_faces),
              "points": "Fluent メッシュの壁ゾーンの面中心 (2 節点の中点; x = 軸方向・y = 半径 [m]、メッシュと同じ座標)",
              "values": "run の正本の実寸の表 (x_m・Tw_x_K) を面中心の x で 1D 線形補間 (表の範囲外は例外)",
              "run_json": files["json"], "run_json_sha256": files["json_sha256"],
              "run_csv": files["csv"], "run_csv_sha256": files["csv_sha256"],
              "table_coord": wtp["rec"]["coord"], "table_T_range_K": [float(min(wtp["rec"]["Tw_x_K"])), float(max(wtp["rec"]["Tw_x_K"]))],
              "fluent_readback": "未確認 (手元に Fluent が無い)"}
    bt = fm.write_fluent_bc_table(d / "boundary_conditions.md", prob, run / "bcondConfig.yaml",
                                  json.loads((run / "prepare_info.json").read_text()),
                                  species_names=_species_names(run), msh_physical=tags, wall_profile=wp)
    res = {"msh": "fluent/nozzle_2d_axi.msh", "check": chk, "cgns": cg, "bc_table": "fluent/boundary_conditions.md",
           "bc_rows": len(bt.get("rows", [])), "problem_yaml": str(problem) if problem else "無し (組成・物性は「不明」と書いた)",
           "pass": bool(chk.get("ok", False) and cg.get("match", False)),
           "fluent_readback": "未確認 (手元に Fluent が無い。plan §5.1 #10 でユーザに確認を依頼)"}
    if wp is not None:
        res["wall_temperature_profile"] = wp
        res["pass"] = bool(res["pass"] and wp["points_match_wall_faces"])
    return res


def _species_names(run: Path):
    """solverConfig の化学種の並び (bcond の Y0, Y1, ... の順)。`species_meta.yaml` の species を使う。無ければ None。"""
    import yaml
    p = run / "species_meta.yaml"
    if not p.exists():
        return None
    sp = (yaml.safe_load(p.read_text()) or {}).get("species")
    return [str(v) for v in sp] if isinstance(sp, list) else None


def _problem_dict(problem) -> dict | None:
    """問題 YAML を読む (prepare_info.json は問題の所在を記録しないので、呼び出し側が --problem で渡す)。無ければ None。"""
    import yaml
    if problem is None:
        return None
    return yaml.safe_load(Path(problem).read_text())


def part_excel(run: Path, out: Path, step=None) -> dict:
    """excel/: 形状表・CFD 壁面量表 (状態つき)・条件・`.xlsx`。壁温分布の run は conditions.csv に壁温の表 (元の座標の種類・元の表・
    実寸の表・区切りの点; `wall_Tw_profile_*`) を足し、xlsx に専用シート `wall_temperature` を足す (run の正本の JSON から;
    plan tooling-nozzle-wall-temperature-profile §4.4)。CFD 壁面量表の `T_w_K` は壁ダンプの `Ts` のまま (分布が入ったことが見える)。"""
    from . import wall_tables as wt
    wtp = _wall_profile(run)                                  # 正本の照合 (不合格なら excel/ に何も書かない)
    d = out / "excel"
    d.mkdir(parents=True, exist_ok=True)
    res = {}
    if (run / "wall_repr.json").exists():
        res["shape"] = wt.shape_table(run, d / "wall_shape_mm.csv")
    else:
        res["shape"] = "飛ばした: wall_repr.json が無い"
    walls = sorted(run.glob("res_wall_*.h5"))
    if walls:
        res["cfd"] = wt.cfd_wall_table(run, d / "wall_cfd_nodes.csv", step=step, allow_diagnostic=True)
        try:
            res["thermal_use"] = wt.thermal_use(run, step)
        except Exception as e:  # noqa: BLE001
            res["thermal_use"] = {"verdict": "未評価", "reason": f"{type(e).__name__}: {e}"}
    else:
        res["cfd"] = "飛ばした: res_wall_*.h5 が無い (CFD 前の run)"
        res["thermal_use"] = {"verdict": "未評価", "reason": "CFD の結果が無い"}
    # 幾何の面積との積分の整合 (§6 V5 の再定義): q × A_node は Q に厳密に戻る。幾何のリング面積 (多角形の台形) との差は面積の定義の差
    # (ソルバの半割面の 1/4 点の r の重み) で、|Σq·ring − ΣQ|/ΣQ ≤ 1e-4 を合格とする (codex result レビュー 2026-10-13 M5)
    cf = d / "wall_cfd_nodes.csv"
    if cf.exists():
        import csv as _csv
        rows = list(_csv.DictReader(open(cf)))
        qk = "q_into_wall_W_m2" if rows and "q_into_wall_W_m2" in rows[0] else None
        if qk and all(r_[qk] not in ("", None) for r_ in rows):
            Qs = np.array([float(r_["Q_node_W"]) for r_ in rows]); qs = np.array([float(r_[qk]) for r_ in rows])
            ring = np.array([float(r_["ring_area_m2"]) for r_ in rows]); An = np.array([float(r_["A_node_m2"]) for r_ in rows])
            # 正規化は Σ|Q| (総熱量が 0 に近い・符号が相殺する壁でも意味を保つ; codex result レビュー 2 回目 M2)
            den = float(np.sum(np.abs(Qs)))
            rel = float(abs(np.sum(qs * ring) - Qs.sum()) / den) if den > 0 else float("nan")
            res["cfd_integral"] = {"sum_Q_W": float(Qs.sum()), "sum_abs_Q_W": den, "sum_q_ring_W": float(np.sum(qs * ring)),
                                   "rel_to_sum_abs_Q": rel, "max_abs_ring_over_A_minus_1": float(np.max(np.abs(ring / An - 1.0))),
                                   "tol": 1e-4, "pass": bool(np.isfinite(rel) and rel <= 1e-4)}
            res["pass"] = res["cfd_integral"]["pass"]           # 総合判定へ伝える (MANIFEST の failed_parts に入る)
    if isinstance(res.get("thermal_use"), dict) and isinstance(res.get("cfd_integral"), dict) and not res["cfd_integral"]["pass"]:
        res["thermal_use"] = dict(res["thermal_use"], verdict="FAIL",
                                  reason=f"幾何の面積との積分の整合が許容外 (相対 {res['cfd_integral']['rel_to_sum_abs_Q']:.3g} > 1e-4)")
    # thermal_use を表の付属 JSON と xlsx に専用の欄で残す (表だけを渡されても使えるかが分かるように; codex result レビュー 2026-10-13 M3)
    tu = res["thermal_use"] if isinstance(res.get("thermal_use"), dict) else {"verdict": "未評価"}
    tu_short = {"verdict": tu.get("verdict"), "step": tu.get("step"), "rule": tu.get("rule"), "reason": tu.get("reason")}
    for k in ("convergence", "quasisteady", "wall_resolution"):
        v = tu.get(k) or {}
        tu_short[k] = {"result": v.get("result"), "verdict_line": v.get("verdict_line"), "reason": v.get("reason")}
    stj = d / "wall_cfd_nodes.status.json"
    if stj.exists():
        sj = json.loads(stj.read_text())
        sj["thermal_use"] = tu_short
        sj["use_for_thermal_design"] = (tu_short["verdict"] == "PASS")
        stj.write_text(json.dumps(_jsonable(sj), indent=1, ensure_ascii=False))
    conds = _conditions(run)
    if wtp is not None:
        wconds = _wall_profile_conditions(wtp)
        import yaml
        bcd = yaml.safe_load((run / "bcondConfig.yaml").read_text()) or {}
        for name, v in bcd.items():                            # 壁の bcond の Ts は代表値 (分布は下の行) と明記する
            if isinstance(v, dict) and v.get("physID") == wtp["rec"]["wall_id"] and f"bc_{name}_Ts" in conds:
                conds[f"bc_{name}_Ts"] = (conds[f"bc_{name}_Ts"][0], "K (代表値 = 設計スロートの値; 壁温は分布 wall_Tw_profile_*)")
        conds.update(wconds)
    conds["thermal_use_verdict"] = (tu_short["verdict"], "PASS のときだけ壁面量表の熱量を熱設計に使ってよい")
    (d / "conditions.csv").write_text("key,value,unit\n" + "".join(f"{k},{v[0]},{v[1]}\n" for k, v in conds.items()))
    wnotes = []
    if wtp is not None:
        f = wtp["files"]
        wnotes = [f"壁温分布: run の正本 {f['json']} (sha256 {f['json_sha256'][:16]}…) の表 → シート wall_temperature・"
                  "conditions の wall_Tw_profile_*。Fluent は fluent/wall_temperature.prof。cfd_wall の T_w_K は壁ダンプの Ts (分布が入ったことが見える)"]
        res["wall_temperature_profile"] = {"conditions_rows": len(wconds), "xlsx_sheet": "wall_temperature", "notes": wnotes,
                                           "run_json": f["json"], "run_json_sha256": f["json_sha256"]}
    if Path(PPTX_PYTHON).exists() and hasattr(wt, "write_xlsx"):
        notes = [f"thermal_use (熱量を熱設計に使えるか): {tu_short['verdict']}  — PASS 以外なら cfd_wall の熱量を使わない",
                 f"  収束: {tu_short['convergence']['result']}  {tu_short['convergence']['verdict_line'] or tu_short['convergence']['reason'] or ''}",
                 f"  熱量の準定常 (総熱量と壁の全節点の熱量の最悪; 16 区間の和は補助): {tu_short['quasisteady']['result']}  {tu_short['quasisteady']['verdict_line'] or tu_short['quasisteady']['reason'] or ''}",
                 f"  壁解像: {tu_short['wall_resolution']['result']}  {tu_short['wall_resolution']['verdict_line'] or tu_short['wall_resolution']['reason'] or ''}",
                 f"  対象の step: {tu_short['step']}"] + wnotes
        res["xlsx"] = wt.write_xlsx(d / "nozzle_tables.xlsx", (d / "wall_shape_mm.csv") if (d / "wall_shape_mm.csv").exists() else None,
                                    (d / "wall_cfd_nodes.csv") if (d / "wall_cfd_nodes.csv").exists() else None,
                                    {k: f"{v[0]} {v[1]}" for k, v in conds.items()}, notes, PPTX_PYTHON,
                                    thermal_use=str(tu_short["verdict"]),
                                    wall_temperature=_wall_profile_sheet(wtp) if wtp is not None else None)
    else:
        res["xlsx"] = "飛ばした: xlsxwriter (.venv-pptx) が無い"
    return res


def _wall_profile_conditions(wtp: dict) -> dict:
    """壁温分布の表を conditions.csv の行 {key: (value, unit)} にする: 正本の所在と sha256、元の座標の種類、区切りの点、各点の
    元の座標・元の温度・実寸の x・実寸の表の温度 (`wall_Tw_profile_<i>_*`)。値・単位の欄にカンマを入れない (引用符なしの 3 列)。"""
    rec, f = wtp["rec"], wtp["files"]
    c = rec["coord"]
    cu = "m" if c == "x_m" else "-"
    out = {"wall_Tw_profile_json": (f["json"], "run の壁温の正本 (実行条件)"),
           "wall_Tw_profile_json_sha256": (f["json_sha256"], ""),
           "wall_Tw_profile_csv": (f["csv"], "forge が読む壁温 (列 x [m] と Ts [K])"),
           "wall_Tw_profile_csv_sha256": (f["csv_sha256"], ""),
           "wall_Tw_profile_coord": (c, "元の座標の種類 (s = 区間で正規化: 入口直管 -2→-1・縮流部 -1→0・拡大部 0→1 / "
                                        "x_m = 実寸 [m]・原点は設計スロート)"),
           "wall_Tw_profile_n_points": (len(rec["coords"]), "点"),
           "wall_Tw_profile_Tw_representative": (rec["Tw_representative_K"], "K (代表値 = 設計スロートの値; bcond の Ts)"),
           "wall_Tw_profile_scale": (rec["scale_m"], "m (写したときの r_t)")}
    labels = {"x_in_m": "m (区切り: 入口端)", "x_contraction_m": "m (区切り: 縮流部の始まり)",
              "x_throat_m": "m (区切り: 設計スロート)", "x_exit_m": "m (区切り: 壁の終点)"}
    for k, v in rec["anchors_m"].items():
        out[f"wall_Tw_profile_anchor_{k}"] = (float(v), labels.get(k, "m"))
    for i, (s, t, x, tx) in enumerate(zip(rec["coords"], rec["Tw_K"], rec["x_m"], rec["Tw_x_K"])):
        out[f"wall_Tw_profile_{i}_orig_{c}"] = (float(s), f"{cu} (元の表)")
        out[f"wall_Tw_profile_{i}_orig_Tw_K"] = (float(t), "K (元の表)")
        out[f"wall_Tw_profile_{i}_x_m"] = (float(x), "m (実寸の表)")
        out[f"wall_Tw_profile_{i}_Tw_x_K"] = (float(tx), "K (実寸の表)")
    bad = [k for k, v in out.items() if "," in f"{v[0]}{v[1]}"]
    if bad:
        raise ValueError(f"conditions.csv の行にカンマが入る: {bad}")
    return out


def _wall_profile_sheet(wtp: dict) -> dict:
    """xlsx のシート `wall_temperature` の中身 (`wall_tables.write_xlsx` の `wall_temperature`)。"""
    rec, f = wtp["rec"], wtp["files"]
    c = rec["coord"]
    return {"lines": [f"壁温分布 (run の正本 {f['json']}、sha256 {f['json_sha256']}; plan tooling-nozzle-wall-temperature-profile §4.2)",
                      f"元の座標の種類: {c} ({'実寸 [m]、原点 = 設計スロート' if c == 'x_m' else '区間で正規化: 入口直管 -2→-1・縮流部 -1→0・拡大部 0→1'})",
                      f"forge は {f['csv']} を壁ノードの x で 1D 線形補間する。Fluent は fluent/wall_temperature.prof (壁面中心で同じ表を評価した point profile)",
                      f"写したときの r_t = {rec['scale_m']!r} m、代表値 (設計スロート) = {rec['Tw_representative_K']!r} K"],
            "anchors": [[k, float(v)] for k, v in rec["anchors_m"].items()],
            "columns": [f"orig_{c}", "orig_Tw_K", "x_m", "Tw_x_K"],
            "rows": [[float(s), float(t), float(x), float(tx)]
                     for s, t, x, tx in zip(rec["coords"], rec["Tw_K"], rec["x_m"], rec["Tw_x_K"])]}


def _conditions(run: Path) -> dict:
    import yaml
    out = {}
    info = json.loads((run / "prepare_info.json").read_text()) if (run / "prepare_info.json").exists() else {}
    if "scale_m" in info:
        out["r_t_design"] = (info["scale_m"], "m")
    bc = yaml.safe_load((run / "bcondConfig.yaml").read_text()) if (run / "bcondConfig.yaml").exists() else {}
    for name, v in (bc or {}).items():
        if isinstance(v, dict):
            out[f"bc_{name}_kind"] = (v.get("kind"), "")
            for k, val in (v.get("floats") or {}).items():
                unit = {"Pt": "Pa(abs)", "Ps": "Pa(abs)", "Tt": "K", "Ts": "K"}.get(k, "")
                out[f"bc_{name}_{k}"] = (val, unit)
    return out


def part_dimensions(run: Path, out: Path, card_path) -> dict:
    import yaml
    from . import dimension_gate as dg
    step = out / "cad" / "wall_physical.step"
    if not step.exists():
        raise _Skip("STEP が無い (cad を先に)")
    sd = dg.step_dimensions(step)
    md = dg.mesh_dimensions(run / "nozzle.msh") if (run / "nozzle.msh").exists() else None
    res = {"step_dimensions": sd, "mesh_dimensions": md}
    if card_path is None:
        res["verdict"] = "未判定"
        res["reason"] = "要求カード (--card) が無いので公差と比べていない (測った寸法だけ)"
    else:
        card = yaml.safe_load(Path(card_path).read_text())
        dt = dg.design_throat_diameter_mm(run) if "diameter_design_m" in (card.get("throat") or {}) else None
        # 弦の形状誤差は配布するメッシュ (run の nozzle.msh = mesh/nozzle.msh) と配布する STEP から測る (codex result レビュー 2 回目 M1)
        chord = dg.chord_from_mesh(step, run / "nozzle.msh") if (run / "nozzle.msh").exists() else None
        res["chord"] = chord
        # 直管と縮流部の接合: STEP の曲線から決めた接合を越える区間の壁節点の半径を曲線と照合 (plan tooling-nozzle-card-local-mesh §5.1 #3j)
        junction = dg.junction_from_mesh(step, run / "nozzle.msh", sd["straight"]["x_end_mm"]) if (run / "nozzle.msh").exists() else None
        res["junction"] = junction
        g = dg.evaluate_gate(card, sd, md, dt, chord, junction)
        res.update(g)
        res["card"] = str(card_path)
        res["card_sha256"] = _sha(Path(card_path))
    (out / "DIMENSIONS.json").write_text(json.dumps(_jsonable(res), indent=1, ensure_ascii=False))
    return {"verdict": res["verdict"], "file": "DIMENSIONS.json"}


def _stl_setting(man: dict) -> str:
    sf = ((man.get("parts") or {}).get("cad") or {}).get("surface") or {}
    if not isinstance(sf, dict) or "stl_linear_deflection_mm" not in sf:
        return "記録なし"
    return (f"線の偏差 {sf['stl_linear_deflection_mm']} mm・角度の偏差 {sf.get('stl_angular_deflection_rad')} rad "
            f"(三角形 {sf.get('stl_triangles')} 枚)")


def _readme_wall_profile(man: dict, wtp: dict) -> list:
    """README の壁温分布の節 (分布の run だけ): 正本・Fluent での設定手順・未確認であること・Excel の所在。"""
    rec, f = wtp["rec"], wtp["files"]
    wp = (man["parts"].get("fluent") or {}).get("wall_temperature_profile")
    anc = "・".join(f"{k} {v!r} m" for k, v in rec["anchors_m"].items())
    L = ["## 壁温分布 (等温壁、分布)", "",
         f"- 正本: run の `{f['json']}` (sha256 `{f['json_sha256'][:16]}…`、元の座標 `{rec['coord']}`・{len(rec['coords'])} 点、"
         f"表の T_w {min(rec['Tw_x_K']):.6g}〜{max(rec['Tw_x_K']):.6g} K; 区切りの点 {anc})。forge は `{f['csv']}` を壁ノードの x で 1D 線形補間する。"
         "成果物の壁温はこの JSON から作った (問題 YAML や bcond の代表 `Ts` では作らない)。書き出しの前に `prepare_info.json` の sha256 と照合した。"]
    if isinstance(wp, dict):
        L += [f"- Fluent: `{wp['file']}` (point profile `{wp['name']}`、場 `{wp['field']}`、{wp['n_points']} 点 = 壁ゾーン `{wp['wall_zone']}` の"
              f"面中心、x = 軸方向・y = 半径 [m] でメッシュと同じ座標、T_w {wp['T_range_K'][0]:.6g}〜{wp['T_range_K'][1]:.6g} K)。設定手順:",
              f"  1. `File > Read > Profile...` で `{Path(wp['file']).name}` を読む (TUI: `/file/read-profile {Path(wp['file']).name}`)。",
              f"  2. `Boundary Conditions > {wp['wall_zone']} > Thermal` で `Temperature` を選び、値の欄で `{wp['name']}` `{wp['field']}` を選ぶ。",
              "  3. point profile は各面に最も近い点の値が入る (既定)。点は壁面中心なので、面ごとに表の値 (面中心で線形補間) がそのまま入る。",
              "- **Fluent での読み込みは未確認** (手元に Fluent が無い)。読み込めない・値の欄に出ない・面の値が表と違うときは教えてほしい。"]
    else:
        L += [f"- Fluent の profile は作っていない (fluent を飛ばした・失敗した: {(man['skipped'].get('fluent') or man['errors'].get('fluent') or '理由不明')})。"
              "**Fluent での読み込みは未確認** (手元に Fluent が無い)。"]
    L += ["- Excel: `excel/conditions.csv` の `wall_Tw_profile_*` (元の表・実寸の表・区切りの点) と `nozzle_tables.xlsx` のシート `wall_temperature`。"
          "CFD 壁面量表 (`wall_cfd_nodes.csv`) の `T_w_K` は壁ダンプの `Ts` (分布が入ったことが見える)。", ""]
    return L


def _readme(run: Path, out: Path, man: dict, wtp: dict | None = None) -> str:
    dim = man["parts"].get("dimensions", {}).get("verdict", "未判定")
    th = (man["parts"].get("excel", {}).get("thermal_use") or {})
    thv = th.get("verdict", "未評価") if isinstance(th, dict) else str(th)
    L = [f"# ノズルの成果物一式 — `{run.name}`", "",
         f"- **寸法の出荷ゲート**: **{dim}** (`DIMENSIONS.json`: STEP とメッシュから測り直した寸法と要求カードの公差)",
         f"- **壁面の熱量を熱設計に使えるか (`thermal_use`)**: **{thv}**",
         "- **Fluent での読み込み**: 未確認 (手元に Fluent が無い。`fluent/nozzle_2d_axi.msh` を読んで `Mesh > Check` の結果を教えてほしい)",
         f"- 総合: **{man['verdict']}**", "",
         "## 単位・座標", "",
         "- 長さ: CAD (STEP・形状表) は mm、メッシュ (msh・Fluent・CGNS) と壁面量表は m。",
         "- x = 流れ方向 (回転軸)、y (= r) = 半径。原点 = 設計スロート (非粘性の設計壁の x = 0)。物理スロート (境界層込みの最小径) の位置は `DIMENSIONS.json`。",
         "- 2D 軸対称の上半分 (r ≥ 0)。Fluent は 2D Axisymmetric、x 軸が対称軸。", "",
         "## フォルダ", "",
         "| フォルダ | 中身 |", "| --- | --- |",
         "| `cad/` | `wall_physical.step` (2D の壁曲線、5 次 B-spline 1 本)・`wall_surface.step` (3D の内面)・`wall_surface.stl`・検査の JSON |",
         "| `mesh/` | `nozzle.msh` (解析に使った msh4.1)・`nozzle_v22.msh` (msh2.2、番号は gmsh が振り直す)・`mesh_structure.json` (構造ファイル: 半径方向の格子線の格子は v1、壁法線の近壁層の格子は v2 = 節点の規則 `wall_normal_blend_v1` と重み β・壁の傾き)・`rebuild_from_step.py` (STEP から同じ構造で作り直す: `python rebuild_from_step.py STEP mesh_structure.json OUT [--compare nozzle.msh]`、gmsh が要る)・その実行例 `rebuilt_from_step*` |",
         "| | 作り直しの出力: `OUT.msh` = forge の変換器に入れる形 (元と同じ番号・書式、座標 17 桁)、`OUT_cad.msh` = 壁の節点が STEP の曲線に分類された gmsh の msh4.1 (**forge の変換器には入れない**: ファイル内の節点の順が番号順でない)、`OUT_v22.msh`・`OUT.cgns` (番号・物理タグは gmsh が振り直す。物理名で対応付ける) |",
         "| `fluent/` | `nozzle_2d_axi.msh` (Fluent 形式、境界の名前と型つき)・`nozzle_2d_axi.cgns`・`boundary_conditions.md` (動作圧 0 Pa で絶対圧を入れる) |",
         "| `excel/` | `wall_shape_mm.csv` (形状表)・`wall_cfd_nodes.csv` (壁節点ごとの圧力・温度・せん断・熱量)・`conditions.csv`・`nozzle_tables.xlsx` |", "",
         "## 注意", "",
         "- forge は節点中心、Fluent はセル中心。同じメッシュでも壁の第 1 未知数までの距離は Fluent が約半分になる。",
         "- 熱量の列は `.status.json` の状態 (`valid` / `not_applicable` / `missing` / `diagnostic` / `nonfinite`) と `thermal_use` を見てから使う。"
         "`diagnostic` (再構成勾配の `qwall`) は正本でない。",
         "- 軸対称の壁への熱量は `iface_Qf_eff` (1 rad あたりの積分値) × 2π から作った (`iface_q_eff` は軸対称では r·q で熱流束ではない)。",
         "- 壁の両端の節点 (入口・出口の角) の熱量は信頼できない (別の壁との境目で再構成勾配との比が 2〜4 倍になった例がある。原因は未確認)。"
         "`excel/wall_cfd_nodes.status.json` の `end_nodes` を見る。",
         "- STEP を CAD で編集したら `rebuild_from_step.py` で作り直し、`check_mesh_quality.py` で格子品質を確かめる。",
         "- 壁法線の近壁層の格子 (`mesh_structure.json` が v2) の作り直しは、STEP の曲線の傾きで近壁の節点を壁の法線上に置く。"
         "入口・出口の面は直線 (x 一定) のままなので、端の列の節点が直線から外れることがある (外れる量は `mesh/rebuilt_from_step_compare.json` の `end_columns`)。"
         "構造ファイルを作るときに、格子の生成の記録 (run の `prepare_info.json` の `mesh.params`) の規則のパラメータで節点を置き直して元と照合した"
         " (`MANIFEST.json` の `structure_rule_check`。記録が無ければ「未実施」)。",
         f"- STL (`cad/wall_surface.stl`) は表示・形の確認用の近似: 弦誤差の設定 {_stl_setting(man)}。寸法の正本は STEP。",
         "- CGNS (`fluent/nozzle_2d_axi.cgns`) を Fluent で読むと**境界の型は付かない** (入口・出口・壁・軸を Fluent 側で付け直す)。"
         "gmsh で読み直すと物理タグの番号も振り直されるので、境界は**名前**で対応付ける。型つきの Fluent 形式は `fluent/nozzle_2d_axi.msh`。", ""]
    if wtp is not None:                                       # 壁温分布の run だけ (一様・断熱の README は変えない)
        L.insert(L.index(f"- 総合: **{man['verdict']}**") + 1,
                 "- **壁温**: 分布 (`fluent/wall_temperature.prof` の point profile; 下の「壁温分布」。Fluent での読み込みは未確認)")
        i = next(k for k, s in enumerate(L) if s.startswith("| `fluent/` |"))
        L[i] = L[i][:-2] + "・`wall_temperature.prof` (壁温の point profile) |"
        L[L.index("## 注意"):L.index("## 注意")] = _readme_wall_profile(man, wtp)
    if man["skipped"]:
        L += ["## 飛ばしたもの", ""] + [f"- `{k}`: {v}" for k, v in man["skipped"].items()] + [""]
    if man["errors"]:
        L += ["## 失敗したもの", ""] + [f"- `{k}`: {v}" for k, v in man["errors"].items()] + [""]
    L += ["## 出典", "", f"- run: `{run}`", f"- 作成: {man['created']}、forge_design の git commit `{man.get('git_commit')}`", ""]
    return "\n".join(L)


def export(run_dir, out_dir=None, card=None, step=None, skip=(), freecadcmd=None, problem=None) -> dict:
    """成果物一式を書き出す。壁温分布の run は、何かを書く前 (出力ディレクトリを作る前) に run の正本 (JSON・CSV) を
    `prepare_info.json` の sha256 と照合し、合わなければ ValueError で止める (部分的な出力を残さない;
    plan tooling-nozzle-wall-temperature-profile §4.2・V1d)。"""
    run = Path(run_dir).resolve()
    wtp = _wall_profile(run)
    out = Path(out_dir).resolve() if out_dir else run / "deliverables"
    out.mkdir(parents=True, exist_ok=True)
    log = _Log()
    if "cad" not in skip:
        log.run("cad", lambda: part_cad(run, out, freecadcmd))
    if "mesh" not in skip:
        log.run("mesh", lambda: part_mesh(run, out))
    if "fluent" not in skip:
        log.run("fluent", lambda: part_fluent(run, out, problem))
    if "excel" not in skip:
        log.run("excel", lambda: part_excel(run, out, step))
    if "dimensions" not in skip:
        log.run("dimensions", lambda: part_dimensions(run, out, card))
    for s in skip:
        log.skipped.setdefault(s, "--skip で指定")
    files = {str(p.relative_to(out)): _sha(p) for p in sorted(out.rglob("*")) if p.is_file() and p.name not in ("MANIFEST.json", "README.md")}
    dimv = log.parts.get("dimensions", {}).get("verdict")
    part_fail = [k for k, v in log.parts.items() if v.get("ok") is False or v.get("pass") is False]
    verdict = "FAIL" if (dimv == "FAIL" or part_fail) else ("PASS" if dimv == "PASS" else "未判定")
    try:
        commit = subprocess.run(["git", "-C", str(HERE), "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    except Exception:  # noqa: BLE001
        commit = None
    man = {"format": "forge.deliverables/v1", "run": str(run), "created": time.strftime("%Y-%m-%d %H:%M:%S"),
           "git_commit": commit, "verdict": verdict, "failed_parts": part_fail, "parts": log.parts,
           "skipped": log.skipped, "errors": log.errors, "files_sha256": files}
    if wtp is not None:                                       # 壁温分布の run だけ (一様・断熱の MANIFEST のキーは変えない)
        f = wtp["files"]
        man["wall_temperature_profile"] = {
            "run_json": f["json"], "run_json_sha256": f["json_sha256"], "run_csv": f["csv"], "run_csv_sha256": f["csv_sha256"],
            "verified": "書き出しの前に prepare_info.json の sha256 と照合した (wall_temperature.load_run_profile)",
            "coord": wtp["rec"]["coord"],
            "fluent_profile": (log.parts.get("fluent") or {}).get("wall_temperature_profile")
            or f"作っていない ({log.skipped.get('fluent') or log.errors.get('fluent') or '理由不明'})"}
    man = _jsonable(man)
    (out / "MANIFEST.json").write_text(json.dumps(man, indent=1, ensure_ascii=False))
    (out / "README.md").write_text(_readme(run, out, man, wtp))
    return man


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("--card", default=None, help="要求カード (forge.wt_card/v1)。無ければ寸法ゲートは「未判定」")
    ap.add_argument("--problem", default=None, help="run を作った問題 YAML (Fluent の境界条件表のガス組成・物性の記載に使う)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--step", type=int, default=None, help="壁面量に使う res_wall の step (既定: 最後)")
    ap.add_argument("--skip", default="", help=f"飛ばす成果物 (カンマ区切り: {','.join(PARTS)})")
    ap.add_argument("--freecadcmd", default=None)
    a = ap.parse_args(argv)
    man = export(a.run_dir, a.out, a.card, a.step, tuple(s for s in a.skip.split(",") if s), a.freecadcmd, a.problem)
    print(json.dumps({"verdict": man["verdict"], "failed_parts": man["failed_parts"], "skipped": man["skipped"],
                      "errors": man["errors"]}, indent=1, ensure_ascii=False))
    print(f"VERDICT: {man['verdict']}")
    return 0 if man["verdict"] in ("PASS", "未判定") and not man["errors"] else 1


if __name__ == "__main__":
    sys.exit(main())
