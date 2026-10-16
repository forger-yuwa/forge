r"""Excel 向けの壁の表 (plans/accepted/tooling-nozzle-wt-deliverables.md §4.5、§5.1 #3・#6)。

- `shape_table`: 物理壁 (`wall_repr.json` の B-spline、STEP と同じ曲線) の形状表 `wall_shape_mm.csv` [mm]。
- `cfd_wall_table`: CFD の壁節点ごとの壁面量表 `wall_cfd_nodes.csv` [SI] と状態 `wall_cfd_nodes.status.json`。
- `wall_heat_series`: 全 `res_<壁名>_<physID>_<step>.h5` の総熱量の時系列 (`check_quasisteady.py --series-csv` の入力)。
- `thermal_use`: 壁面量を熱設計に使えるかの判定 (収束・熱量の準定常 [総熱量と壁の全節点の熱量、16 区間の和は補助]・壁解像・熱量の状態)。
  収束の道具は run_case.sh・受け渡しの関門と同じ規則 (`settle.select_judge`; plan tooling-settle-convergence-and-monitor §5.1 #7b):
  `settle_criteria.yaml` のある run の本段は `check_settled.py --bind` (停止規則 v2 の再評価と照合)、それ以外は今どおり `check_convergence.py`。
- `write_xlsx`: 上の CSV と条件・注記を `.xlsx` に (`.venv-pptx` の xlsxwriter を subprocess で呼ぶ)。
- `table_x_mm`: 形状表の x の点 (粗 1 mm・スロート付近 0.1 mm)。設計壁と物理壁の比較の表 (`wall_compare.py`) も同じ点を使う。

**熱量の正本** (plan §3・§4.5、codex plan レビュー M1・M3): 壁へ入る熱は積分済みの節点熱量 `iface_Qf_eff` (固体向き正、
`methods/boundary.md`「界面熱量の 4 つの定義」) から作る。軸対称 (`mesh.isAxisymmetric: 1`・`axisymMethod: 0`) では流体の面ベクトルが
デバイス側で $r$ 倍 (`solver_density_cuda/variables.cpp` の `ss[ip] *= r_face`) なので `iface_Qf_eff` は 1 rad あたり [W/rad]。
全周の節点熱量 $Q_{node}=2\pi\,$`iface_Qf_eff`、熱流束 $q=Q_{node}/A_{node}$、$A_{node}=2\pi\,\max(r_f, r_{floor})\,|S_f|$
($r_f$ = 壁の境界面の重心の $y$、$r_{floor}$ はソルバと同じ `mesh.axisRFloor` か 1e-20)。
`iface_q_eff` は使わない (main のソルバでは軸対称で $r\,q$ になる; ソルバの版で意味が変わる量に依存しない)。
平面 (非軸対称) では $Q_{node}=$`iface_Qf_eff` [W/m]、$A_{node}=|S_f|$ [m²/m]。
`iface_*` は等温壁・定常・`output.interfaceDiag: 1` のときだけ値が入り、それ以外は**データセットはあって中身が NaN**。
状態は `valid` / `not_applicable` (全 NaN) / `missing` (データセットなし) / `nonfinite` (一部 NaN) / `diagnostic` (`-qwall` を使った)。
正本が無効なとき**黙って診断値に置き換えない**: `allow_diagnostic=True` のときだけ `-qwall` を列名 `q_into_wall_diag_W_m2` で出す。

$y_1^+$ は `solver_density_cuda/tools/check_wall_resolution.py` の関数 (`dof_neighbors`・`wall_first_distance`・`mu_of`) を import して
同じ定義 (第 1 内部節点までの壁法線距離・接線応力・ソルバの `vis_lam` を優先) で出す。ソルバの `ypls` は載せない。

usage (design/.venv-opt):
  python -m forge_design.export.wall_tables shape WALL_DIR_OR_JSON OUT.csv
  python -m forge_design.export.wall_tables cfd RUN_DIR OUT.csv [--step N] [--allow-diagnostic] [--walls 3,4]
  python -m forge_design.export.wall_tables series RUN_DIR OUT.csv [--walls 3]
  python -m forge_design.export.wall_tables thermal RUN_DIR [--step N] [--walls 3]
"""
from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
TOOLS = REPO / "solver_density_cuda" / "tools"
def _pptx_python() -> str:
    """xlsxwriter を入れた Python: 環境変数 FORGE_PPTX_PY → リポジトリ直下の .venv-pptx → 主の作業ツリーの .venv-pptx (以前はこのマシンの絶対パス)。"""
    root = Path(__file__).resolve().parents[3]
    if os.environ.get("FORGE_PPTX_PY"):
        return os.environ["FORGE_PPTX_PY"]
    for r in (root, root.parent / "forge"):
        c = r / ".venv-pptx" / "bin" / "python"
        if c.exists():
            return str(c)
    return str(root / ".venv-pptx" / "bin" / "python")


PPTX_PY = _pptx_python()
WALL_KINDS = ("wall", "wall_isothermal")
SHAPE_COLUMNS = ("x_mm", "r_mm", "wall_angle_deg", "curvature_radius_mm", "s_mm", "perimeter_mm", "lateral_area_cum_mm2")
# 曲率がこれ以下なら直線とみなして曲率半径を inf にする [1/mm] (STEP の転送誤差の曲率の許容 wall_step.TRANSFER_TOL と同じ値)
KAPPA_ZERO_PER_MM = 1e-9
HEAT_STATUSES = ("valid", "not_applicable", "missing", "nonfinite", "diagnostic")
_DUMP_RE = re.compile(r"res_(.+)_(\d+)_(\d+)\.h5$")


def _fmt(v) -> str:
    """CSV の数値: 最短で往復する表記 (repr)。NaN は空欄、±inf は `inf`/`-inf`。"""
    if v is None:
        return ""
    if isinstance(v, (int, np.integer)) and not isinstance(v, bool):
        return str(int(v))
    f = float(v)
    if math.isnan(f):
        return ""
    if math.isinf(f):
        return "inf" if f > 0 else "-inf"
    return repr(f)


def _write_csv(path: Path, columns, cols_data) -> None:
    """列名だけのヘッダ 1 行 + 数値行 (コメント行を入れない: Excel が崩れる)。列名は ASCII に限る。"""
    for c in columns:
        if not c.isascii():
            raise ValueError(f"CSV の列名は ASCII に限る: {c!r}")
    n = len(cols_data[0]) if cols_data else 0
    with open(path, "w", newline="") as f:
        f.write(",".join(columns) + "\n")
        for i in range(n):
            f.write(",".join(_fmt(c[i]) for c in cols_data) + "\n")


# ---------------------------------------------------------------------------------------------------------------------
# 形状表
# ---------------------------------------------------------------------------------------------------------------------

def table_x_mm(X0: float, X1: float, s_mm: float, coarse_mm: float = 1.0, fine_per_mm: int = 10,
               throat_half_width_mm: float | None = None) -> tuple:
    """表の x の点 [mm] (形状表と壁の比較の表で共通): `coarse_mm` ごと + スロート付近 ($|x| \\le$ hw) は 1/`fine_per_mm` mm ごと、
    両端 X0・X1 を含む。hw の既定は $\\max(50\\,\\mathrm{mm}, 2 r_t)$ (`s_mm` = r_t [mm])。戻り: (x の点, hw)。"""
    hw = max(50.0, 2.0 * s_mm) if throat_half_width_mm is None else float(throat_half_width_mm)
    k0, k1 = math.ceil(X0 / coarse_mm), math.floor(X1 / coarse_mm)
    xc = np.arange(k0, k1 + 1, dtype=float) * coarse_mm
    f0, f1 = math.ceil(max(-hw, X0) * fine_per_mm), math.floor(min(hw, X1) * fine_per_mm)
    xf = np.arange(f0, f1 + 1, dtype=float) / float(fine_per_mm)  # k/10 は 10 進の 0.1 刻みに最も近い浮動小数 (repr で短く出る)
    xm = np.unique(np.r_[xc, xf, X0, X1])
    return xm[(xm >= X0) & (xm <= X1)], hw


def shape_table(wall_file_dir_or_json, out_csv, coarse_mm: float = 1.0, fine_per_mm: int = 10,
                throat_half_width_mm: float | None = None, n_sub: int = 16) -> dict:
    """物理壁の形状表 (plan §4.5): 1 mm ごと + スロート付近 ($|x| \\le \\max(50\\,\\mathrm{mm}, 2 r_t)$) は 0.1 mm ごと、両端を含む。

    列: `x_mm, r_mm, wall_angle_deg, curvature_radius_mm, s_mm, perimeter_mm, lateral_area_cum_mm2`。
    - $r$・$r'$・$r''$ は `load_wall_file(...)["physical"].r` (保存した B-spline、r_t 単位) を mm に換算して評価する。
    - 曲率半径は符号つき $1/\\kappa$、$\\kappa = r''/(1+r'^2)^{3/2}$ [1/mm]。正 = 流路側から見て凸 (スロートの周り)。
      $|\\kappa| \\le$ `KAPPA_ZERO_PER_MM` は `inf` (直管)。
    - 弧長 $s$ と累積側面積 $\\int 2\\pi r\\,ds$ は、**表の点の間を `n_sub` 分割した細かい点の台形積分の累積** (入口端から)。
      精度の目安は同じ点の 1 つおき (`n_sub/2` 分割) との差 / 3 (Richardson) を戻り値に入れる。
    戻り: 説明 (列の意味・単位・原点・刻み・積分の精度) の dict。呼び出し側が README に書く。"""
    from ..geometry.wall_axismach import load_wall_file
    if n_sub < 2 or n_sub % 2:
        raise ValueError(f"n_sub は 2 以上の偶数 (受け取った {n_sub})")
    W = load_wall_file(wall_file_dir_or_json)
    phys = W["physical"]
    s_mm = 1000.0 * float(W["scale_m"])                       # r_t → mm
    x0, x1 = (float(v) for v in W["domain"])                  # r_t
    X0, X1 = x0 * s_mm, x1 * s_mm
    xm, hw = table_x_mm(X0, X1, s_mm, coarse_mm, fine_per_mm, throat_half_width_mm)

    def ev(x_mm, d=0):
        return phys.r(np.clip(np.asarray(x_mm, dtype=float) / s_mm, x0, x1), d)

    r = ev(xm) * s_mm
    r1 = ev(xm, 1)
    r2 = ev(xm, 2) / s_mm
    ang = np.degrees(np.arctan(r1))
    kap = r2 / (1.0 + r1 ** 2) ** 1.5
    big = np.abs(kap) > KAPPA_ZERO_PER_MM
    rc = np.where(big, 1.0 / np.where(big, kap, 1.0), np.inf)
    # 細かい点: 各区間 [x_i, x_{i+1}] を n_sub 等分
    dx = np.diff(xm)
    j = np.arange(n_sub) / n_sub
    xs = np.r_[(xm[:-1, None] + dx[:, None] * j[None, :]).ravel(), xm[-1]]
    rs = ev(xs) * s_mm
    gs = np.sqrt(1.0 + ev(xs, 1) ** 2)

    def cum(f, every):
        """各区間の台形積分 (every 点おき) の累積。"""
        m = len(xm) - 1
        fi = f[:-1].reshape(m, n_sub)[:, ::every]
        fe = np.c_[fi[:, 1:], f[n_sub::n_sub][:, None]]
        h = (dx / (n_sub // every))[:, None]
        return np.r_[0.0, np.cumsum((0.5 * (fi + fe) * h).sum(axis=1))]

    s16, a16 = cum(gs, 1), cum(2.0 * np.pi * rs * gs, 1)
    s8, a8 = cum(gs, 2), cum(2.0 * np.pi * rs * gs, 2)
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(out_csv, SHAPE_COLUMNS, [xm, r, ang, rc, s16, 2.0 * np.pi * r, a16])
    th = W["throat"]
    return {
        "csv": str(out_csv), "n_rows": int(len(xm)), "columns": list(SHAPE_COLUMNS),
        "source": str(W["path"]), "physical_wall_repr": W["kind"], "scale_mm_per_rt": s_mm,
        "units": "mm (角度は deg)", "origin": "設計スロート (設計壁の x = 0)。x = 流れ方向 (回転軸)、r = 半径",
        "domain_mm": [X0, X1], "spacing": {"coarse_mm": coarse_mm, "fine_mm": 1.0 / fine_per_mm,
                                           "fine_zone_mm": [max(-hw, X0), min(hw, X1)],
                                           "fine_zone_rule": "|x| ≤ max(50 mm, 2 r_t)"},
        "physical_throat_mm": {"x": float(th["x"]) * s_mm, "r": float(th["r"]) * s_mm},
        "column_notes": {
            "x_mm": "軸方向の位置 [mm]。原点 = 設計スロート",
            "r_mm": "物理壁の半径 [mm] (wall_repr.json の物理壁 = STEP と同じ B-spline を評価)",
            "wall_angle_deg": "壁の傾き atan(dr/dx) [deg]。正 = 下流で広がる",
            "curvature_radius_mm": f"符号つき曲率半径 1/κ [mm]、κ = r''/(1+r'^2)^1.5。正 = 流路側から見て凸。|κ| ≤ {KAPPA_ZERO_PER_MM:g} /mm は inf (直線)",
            "s_mm": "入口端からの壁の弧長 [mm]",
            "perimeter_mm": "周長 2πr [mm]",
            "lateral_area_cum_mm2": "入口端からの累積側面積 ∫2πr ds [mm²]"},
        "integration": {"method": f"表の点の間を {n_sub} 分割した台形積分の累積",
                        "s_total_mm": float(s16[-1]), "lateral_area_total_mm2": float(a16[-1]),
                        "s_err_est_mm": float(abs(s16[-1] - s8[-1]) / 3.0),
                        "area_err_est_mm2": float(abs(a16[-1] - a8[-1]) / 3.0),
                        "area_err_est_rel": float(abs(a16[-1] - a8[-1]) / 3.0 / a16[-1]),
                        "err_est_note": f"{n_sub} 分割と {n_sub // 2} 分割の差 / 3 (Richardson)"},
    }


# ---------------------------------------------------------------------------------------------------------------------
# CFD の壁面量表
# ---------------------------------------------------------------------------------------------------------------------

def _load_yaml(p: Path):
    import yaml
    return yaml.safe_load(p.read_text()) if p.exists() else None


def _run_config(run_dir: Path) -> dict:
    """run の solverConfig と壁の bcond、メッシュの場所、軸対称の設定 (ソルバと同じ読み方: mesh ブロックが正本、physProp は旧)。"""
    cfg = _load_yaml(run_dir / "solverConfig.yaml") or {}
    bc = _load_yaml(run_dir / "bcondConfig.yaml") or {}
    walls = {}
    for name, v in bc.items():
        if isinstance(v, dict) and str(v.get("kind", "")) in WALL_KINDS:
            walls[int(v["physID"])] = {"name": str(name), "kind": str(v["kind"])}
    m = cfg.get("mesh") or {}
    pp = cfg.get("physProp") or {}
    axi = int(m.get("isAxisymmetric", pp.get("isAxisymmetric", 0)) or 0)
    meth = int(m.get("axisymMethod", pp.get("axisymMethod", 0)) or 0)
    rfl = float(m.get("axisRFloor", 0.0) or 0.0)
    # メッシュ名は run の config から取る (check_wall_resolution.py と同じ: mesh.meshFileName、無ければ mesh.h5)
    mesh = None
    mf = m.get("meshFileName")
    for cand in ([run_dir / mf] if mf else []) + [run_dir / "mesh.h5"]:
        if cand.exists():
            mesh = cand
            break
    return {"cfg": cfg, "walls": walls, "isAxisymmetric": axi, "axisymMethod": meth,
            "r_floor": rfl if rfl > 0.0 else 1.0e-20, "mesh": mesh}


def _wall_dumps(run_dir: Path, walls: dict) -> dict:
    """{physID: {step: path}}。ファイル名 `res_<壁名>_<physID>_<step>.h5` の壁名が bcond の名前と一致するものだけ
    (同じ physID の別物 `res_solid_<physID>_*` などを拾わない)。"""
    out = {pid: {} for pid in walls}
    for p in sorted(run_dir.glob("res_*_*_*.h5")):
        m = _DUMP_RE.search(p.name)
        if not m:
            continue
        name, pid, st = m.group(1), int(m.group(2)), int(m.group(3))
        if pid in walls and walls[pid]["name"] == name:
            if st in out[pid]:
                raise ValueError(f"同じ壁・同じ step の壁ダンプが 2 つある: {out[pid][st]} と {p}")
            out[pid][st] = p
    return out


def _select_step(dumps: dict, step) -> int:
    """使う step を 1 つに決める。step 未指定なら各壁の最終 step が揃っていることを要求する (異なる step の混入は例外)。"""
    pids = sorted(dumps)
    if not pids:
        raise ValueError("対象の壁が無い (bcondConfig.yaml の kind が wall / wall_isothermal の physID)")
    for pid in pids:
        if not dumps[pid]:
            raise FileNotFoundError(f"壁 physID {pid} の壁ダンプ res_<壁名>_{pid}_<step>.h5 が無い")
    if step is None:
        last = {pid: max(dumps[pid]) for pid in pids}
        if len(set(last.values())) != 1:
            raise ValueError(f"壁ごとの最終 step が違う (異なる step の混入): {last}。step を指定すること")
        return int(next(iter(last.values())))
    step = int(step)
    miss = [pid for pid in pids if step not in dumps[pid]]
    if miss:
        raise FileNotFoundError(f"step {step} の壁ダンプが無い壁: {miss} (異なる step を混ぜない)")
    return step


def _read_dump(path: Path, step: int) -> tuple[dict, np.ndarray, bool]:
    """壁ダンプの値・座標 (倍精度に直す) と、座標が float32 で書かれていたか (forge の float ビルド)。"""
    import h5py
    with h5py.File(path, "r") as f:
        st = f.attrs.get("step")
        if st is not None and int(st) != int(step):
            raise ValueError(f"{path.name}: 属性 step {int(st)} ≠ ファイル名の step {step} (異なる step の混入)")
        V = {k: np.array(f["VALUE"][k]) for k in f["VALUE"]}
        xyz32 = f["MESH/COORD"].dtype == np.float32
        xyz = np.array(f["MESH/COORD"], dtype=float).reshape(-1, 3)
    return V, xyz, bool(xyz32)


def _chain(edges: np.ndarray, xyz: np.ndarray) -> list:
    """壁の辺 (節点の組) を 1 本の開いた折れ線にたどる。始点は x の小さい端。分岐・閉路・複数本は例外。"""
    adj: dict = {}
    for a, b in edges:
        a, b = int(a), int(b)
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    ends = [n for n, v in adj.items() if len(v) == 1]
    if any(len(v) > 2 for v in adj.values()) or len(ends) != 2:
        raise ValueError(f"壁の辺が 1 本の開いた折れ線でない (端 {len(ends)} 個、分岐 {sum(len(v) > 2 for v in adj.values())} 個)")
    start = min(ends, key=lambda n: (xyz[n, 0], xyz[n, 1]))
    out, prev, cur = [start], None, start
    while True:
        nxt = [n for n in adj[cur] if n != prev]
        if not nxt:
            break
        prev, cur = cur, nxt[0]
        out.append(cur)
    if len(out) != len(adj):
        raise ValueError(f"壁の辺が 1 本につながらない ({len(out)} / {len(adj)} 節点)")
    return out


def _wall_block(run_dir: Path, rc: dict, pid: int, path: Path, step: int, mesh_cache: dict, with_y1plus: bool) -> dict:
    """1 つの壁の節点ごとの量 (折れ線の順)。"""
    import h5py
    from scipy.spatial import cKDTree
    mesh = rc["mesh"]
    if "nodes" not in mesh_cache:
        with h5py.File(mesh, "r") as f:
            mesh_cache["nodes"] = np.array(f["MESH/COORD"], dtype=float).reshape(-1, 3)
            mesh_cache["sa"] = np.array(f["PLANES/surfArea"], dtype=float)
            mesh_cache["pc"] = np.array(f["PLANES/centCoords"], dtype=float).reshape(-1, 3)
            mesh_cache["bconds"] = sorted(f["BCONDS"].keys()) if "BCONDS" in f else []
    nodes = mesh_cache["nodes"]
    if str(pid) not in mesh_cache["bconds"]:
        raise ValueError(f"メッシュ {mesh} に BCONDS/{pid} が無い")
    with h5py.File(mesh, "r") as f:
        g = f[f"BCONDS/{pid}"]
        icells = np.array(g["iCells"]).astype(np.int64)
        iplanes = np.array(g["iPlanes"]).astype(np.int64)
        vn = np.array(g["vizBfaceNodes"]).astype(np.int64)
        vs = np.array(g["vizBfaceSizes"]).astype(np.int64)
    if len(np.unique(icells)) != len(icells):
        raise ValueError(f"BCONDS/{pid}: 同じ節点に境界面が 2 つ以上ある (node 方式の 2D 壁だけ対応)")
    if icells.max() >= len(nodes):
        raise ValueError(f"BCONDS/{pid}: iCells が節点数を超える (node 方式のメッシュだけ対応)")
    if not np.all(vs == 2):
        raise ValueError(f"BCONDS/{pid}: 境界面が線分でない (2D の壁だけ対応)")
    V, wxyz, wxyz32 = _read_dump(path, step)
    n = len(V["ro"]) if "ro" in V else len(wxyz)
    if n != len(icells) or len(wxyz) != len(icells):
        raise ValueError(f"{path.name}: 壁ダンプの点数 {n} (座標 {len(wxyz)}) ≠ BCONDS/{pid} の境界面数 {len(icells)} (node 方式の壁ダンプだけ対応)")
    # 壁ダンプの点を座標で境界面 (ib) に対応づける (順序を仮定しない; check_wall_resolution.py と同じ許容)
    d, jmap = cKDTree(nodes[icells]).query(wxyz)
    scale = float(np.median(np.linalg.norm(np.diff(wxyz[:min(len(wxyz), 200)], axis=0), axis=1))) if len(wxyz) > 1 else 1.0
    tol = 1e-6 * max(scale, 1e-12) + 1e-12
    # forge の float ビルドの壁ダンプは座標が float32 で、倍精度のメッシュの節点から丸め幅 (≤ 0.5 ulp ≈ 6e-8 |x|) だけずれる
    # (2026-10-16 試走 case/90 run_0005_ns_p1: 最大 2.4e-7 m で止まった)。check_wall_resolution.py と同じく、丸め幅を許容に足すのは
    # 壁の節点の最小間隔の 1/4 未満のとき (対応づけが一意である範囲) だけ
    if wxyz32 and len(icells) > 1:
        tol32 = 2.0 * float(np.finfo(np.float32).eps) * float(np.abs(wxyz).max())
        dd, _ = cKDTree(nodes[icells]).query(nodes[icells], k=2)
        if tol32 < 0.25 * float(dd[:, 1].min()):
            tol = max(tol, tol32)
    if float(d.max()) > tol or len(np.unique(jmap)) != len(jmap):
        raise ValueError(f"{path.name}: 壁ダンプの点を節点に 1 対 1 で対応づけられない (最大 {float(d.max()):.3e} m)")
    k_of_ib = np.empty(len(icells), dtype=np.int64)
    k_of_ib[jmap] = np.arange(len(jmap))
    chain = _chain(vn.reshape(-1, 2), nodes)
    if set(chain) != set(int(c) for c in icells):
        raise ValueError(f"BCONDS/{pid}: 折れ線の節点と iCells が一致しない")
    ib_of_node = {int(c): i for i, c in enumerate(icells)}
    ib = np.array([ib_of_node[c] for c in chain], dtype=np.int64)
    k = k_of_ib[ib]
    P = nodes[chain, :2]
    x, r = P[:, 0], P[:, 1]
    seg = np.diff(P, axis=0)
    L = np.hypot(seg[:, 0], seg[:, 1])
    s = np.r_[0.0, np.cumsum(L)]
    ring_e = np.pi * (r[:-1] + r[1:]) * L
    ring = np.zeros(len(chain))
    ring[:-1] += 0.5 * ring_e
    ring[1:] += 0.5 * ring_e
    # 壁の接線 (折れ線の中心差分、端は片側)。向きは折れ線の向き (x の小さい端から) = 流れ方向 +x
    tg = np.empty_like(P)
    tg[1:-1] = P[2:] - P[:-2]
    tg[0], tg[-1] = P[1] - P[0], P[-1] - P[-2]
    tg /= np.linalg.norm(tg, axis=1)[:, None]
    if "twall_x" in V and "twall_y" in V:
        tau = V["twall_x"].astype(float)[k] * tg[:, 0] + V["twall_y"].astype(float)[k] * tg[:, 1]
    else:
        tau = np.full(len(chain), np.nan)                   # 壁ダンプに twall が無い → 空欄
    sa, pc = mesh_cache["sa"][iplanes[ib]], mesh_cache["pc"][iplanes[ib]]
    axi = rc["isAxisymmetric"] == 1
    A_node = 2.0 * np.pi * np.maximum(pc[:, 1], rc["r_floor"]) * sa if axi else sa.copy()
    # 熱量の正本 iface_Qf_eff
    heat = {"dataset": "iface_Qf_eff", "file": path.name}
    if axi and rc["axisymMethod"] != 0:
        heat.update(status="not_applicable",
                    reason=f"axisymMethod {rc['axisymMethod']} (平面の幾何 + 1/y ソース) の iface_Qf_eff は単位が未検証")
        Q = np.full(len(chain), np.nan)
    elif "iface_Qf_eff" not in V:
        heat.update(status="missing", reason="壁ダンプに iface_Qf_eff が無い (output.interfaceDiag: 1 でない run)")
        Q = np.full(len(chain), np.nan)
    else:
        qf = V["iface_Qf_eff"].astype(float)[k]
        nfin = int(np.isfinite(qf).sum())
        heat.update(n=int(len(qf)), n_finite=nfin)
        if nfin == len(qf):
            heat["status"] = "valid"
        elif nfin == 0:
            heat.update(status="not_applicable", reason="iface_Qf_eff が全 NaN (断熱壁・非定常・等温でない壁など)")
        else:
            heat.update(status="nonfinite", reason=f"iface_Qf_eff の {len(qf) - nfin} / {len(qf)} 点が非有限")
        Q = (2.0 * np.pi * qf) if axi else qf
    qdiag = (-V["qwall"].astype(float)[k]) if "qwall" in V else None
    blk = {"pid": pid, "name": rc["walls"][pid]["name"], "kind": rc["walls"][pid]["kind"], "file": path.name,
           "x": x, "r": r, "s": s, "ring": ring, "A_node": A_node, "p": V["Ps"].astype(float)[k] if "Ps" in V else np.full(len(k), np.nan),
           "T": V["Ts"].astype(float)[k] if "Ts" in V else np.full(len(k), np.nan), "tau": tau,
           "Q": Q, "q": Q / A_node, "qdiag": qdiag, "heat": heat, "polygon_lateral_area": float(ring_e.sum()),
           "x_monotone": bool(np.all(np.diff(x) > 0.0)), "y1plus": np.full(len(chain), np.nan), "y1plus_note": None}
    if with_y1plus:
        blk["y1plus"], blk["y1plus_note"] = _y1plus(run_dir, rc, pid, step, V, wxyz, icells, jmap, k, mesh_cache)
    return blk


def _cwr():
    """check_wall_resolution.py を import (式を書き写さない)。"""
    if str(TOOLS) not in sys.path:
        sys.path.insert(0, str(TOOLS))
    import check_wall_resolution as cwr  # noqa: PLC0415
    return cwr


def _y1plus(run_dir, rc, pid, step, V, wxyz, icells, jmap, k, mesh_cache, align_min: float = 0.5):
    """check_wall_resolution.py の main と同じ定義の y1+ (壁ダンプの点の順) を折れ線の順 k で返す。"""
    import h5py
    need = [c for c in ("twall_x", "twall_y", "twall_z", "ro", "Ts") if c not in V]
    if need:
        return np.full(len(k), np.nan), f"壁ダンプに {need} が無い → 空欄 (check_wall_resolution.py も判定不能にする)"
    cwr = _cwr()
    mesh = str(rc["mesh"])
    if "nb" not in mesh_cache:
        mesh_cache["nb"] = cwr.dof_neighbors(mesh)
    uniq = icells[jmap]                                      # 壁ダンプの各点に対応する DOF
    wd = cwr.wall_first_distance(mesh, pid, mesh_cache["nb"], align_min, coords=mesh_cache["nodes"])
    y1 = np.array([wd.get(int(c), (np.nan, None, False))[0] for c in uniq])
    okc = np.array([wd.get(int(c), (np.nan, None, False))[2] for c in uniq])
    nvec = np.array([wd.get(int(c), (np.nan, np.zeros(3), False))[1]
                     if wd.get(int(c), (np.nan, None, False))[1] is not None else np.zeros(3) for c in uniq])
    tw = np.stack([V["twall_x"], V["twall_y"], V["twall_z"]], axis=1)
    tn = np.sum(tw * nvec, axis=1)[:, None] * nvec
    tt = np.linalg.norm(tw - tn, axis=1)                    # 接線成分
    mu, src = None, None
    vpath = Path(run_dir) / f"res_{int(step)}.h5"
    if vpath.exists():
        with h5py.File(vpath, "r") as fv:
            if "VALUE/vis_lam" in fv:
                vl = np.array(fv["VALUE/vis_lam"])
                if len(vl) > int(np.max(uniq)):
                    mu, src = vl[uniq].astype(float), f"{vpath.name} の vis_lam"
    if mu is None:
        mu = cwr.mu_of(V["Ts"].astype(float), rc["cfg"])
        src = "mu_of (viscMethod の式)" if mu is not None else None
    if mu is None:
        return np.full(len(k), np.nan), "粘性を解決できない (vis_lam も viscMethod の式も無い) → 空欄"
    ro = V["ro"].astype(float)
    yp = np.where(okc, y1 * np.sqrt(np.maximum(ro, 0) * tt) / np.maximum(mu, 1e-30), np.nan)
    return yp[k], f"check_wall_resolution.py と同じ定義 (第 1 内部節点の壁法線距離・接線応力、粘性 = {src})"


def _build_table(run_dir, step=None, allow_diagnostic: bool = False, wall_ids=None, with_y1plus: bool = True) -> dict:
    run_dir = Path(run_dir)
    rc = _run_config(run_dir)
    if rc["mesh"] is None:
        raise FileNotFoundError(f"{run_dir}: メッシュ h5 が無い (solverConfig の mesh.meshFileName と mesh.h5 を見た)")
    walls = rc["walls"]
    if wall_ids is not None:
        bad = [w for w in wall_ids if int(w) not in walls]
        if bad:
            raise ValueError(f"physID {bad} は bcondConfig.yaml の壁 (kind {WALL_KINDS}) でない (壁: {sorted(walls)})")
        walls = {int(w): walls[int(w)] for w in wall_ids}
        rc = dict(rc, walls=walls)
    dumps = _wall_dumps(run_dir, walls)
    st = _select_step(dumps, step)
    cache: dict = {}
    blocks = [_wall_block(run_dir, rc, pid, dumps[pid][st], st, cache, with_y1plus) for pid in sorted(walls)]
    blocks.sort(key=lambda b: (float(b["x"][0]), b["pid"]))
    per = {b["pid"]: b["heat"] for b in blocks}
    stats = [h["status"] for h in per.values()]
    if all(s_ == "valid" for s_ in stats):
        overall = "valid"
    else:
        overall = next(s_ for s_ in ("nonfinite", "missing", "not_applicable") if s_ in stats)
    use_diag = False
    if overall != "valid" and allow_diagnostic:
        if all(b["qdiag"] is not None for b in blocks):
            use_diag, overall = True, "diagnostic"
    return {"run_dir": run_dir, "rc": rc, "step": st, "blocks": blocks, "status": overall, "per_wall": per,
            "use_diag": use_diag, "dumps": dumps}


def cfd_wall_table(run_dir, out_csv, step=None, allow_diagnostic: bool = False, wall_ids=None,
                   with_y1plus: bool = True) -> dict:
    """CFD の壁面量表 (plan §4.5): 壁節点ごとに
    `x_m, r_m, s_m, ring_area_m2, A_node_m2, p_w_Pa, T_w_K, tau_w_Pa, q_into_wall_W_m2, Q_node_W, y1plus, wall_id`。

    - 壁 = `bcondConfig.yaml` の kind が `wall`/`wall_isothermal` の physID (`wall_ids` で絞れる)。壁ごとに x の小さい端からの折れ線の順
      (= x 昇順。壁が複数なら始点の x の順に並べる)。
    - `ring_area_m2`: 隣り合う壁節点を結ぶ辺の側面積 π(r_a+r_b)|ab| を両端に半分ずつ配った持ち分 (総和 = 折れ線の側面積)。
    - `A_node_m2`: ソルバの面 (壁の境界面の r 重み面積 2π max(r_f, r_floor)|S_f|、平面は |S_f|)。q の分母。
    - `tau_w_Pa`: `twall_{x,y}` の折れ線の接線への射影 (流れ方向 +x を正)。
    - 熱量: `iface_Qf_eff` から (モジュールの説明)。状態が `valid` でなければ q・Q は空欄。`allow_diagnostic=True` のときだけ
      `-qwall` を `q_into_wall_diag_W_m2`・`Q_node_diag_W` (= q_diag × A_node) で出し状態 `diagnostic`。
    - 使う `res_wall` の step は 1 つ (`step` か各壁の最終。異なる step の混入は例外)。
    CSV と同じ名前の `.status.json` に状態・出典 step・データセット名を書く。戻り: その内容。"""
    T = _build_table(run_dir, step, allow_diagnostic, wall_ids, with_y1plus)
    blocks, use_diag = T["blocks"], T["use_diag"]
    qn, Qn = ("q_into_wall_diag_W_m2", "Q_node_diag_W") if use_diag else ("q_into_wall_W_m2", "Q_node_W")
    cols = ("x_m", "r_m", "s_m", "ring_area_m2", "A_node_m2", "p_w_Pa", "T_w_K", "tau_w_Pa", qn, Qn, "y1plus", "wall_id")
    data = {c: [] for c in cols}
    for b in blocks:
        if use_diag:
            q = b["qdiag"]
            Q = q * b["A_node"]
        elif b["heat"]["status"] == "valid":
            q, Q = b["q"], b["Q"]
        else:
            q = Q = np.full(len(b["x"]), np.nan)
        for c, v in (("x_m", b["x"]), ("r_m", b["r"]), ("s_m", b["s"]), ("ring_area_m2", b["ring"]),
                     ("A_node_m2", b["A_node"]), ("p_w_Pa", b["p"]), ("T_w_K", b["T"]), ("tau_w_Pa", b["tau"]),
                     (qn, q), (Qn, Q), ("y1plus", b["y1plus"]), ("wall_id", np.full(len(b["x"]), float(b["pid"])))):
            data[c].extend(np.asarray(v, dtype=float).tolist())
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    # wall_id は整数で書く
    _write_csv(out_csv, cols, [data[c] if c != "wall_id" else [int(v) for v in data[c]] for c in cols])
    rc = T["rc"]
    axi = rc["isAxisymmetric"] == 1
    walls_info = {}
    for b in blocks:
        h = dict(b["heat"])
        valid = h["status"] == "valid"
        rel_ring = np.abs(b["A_node"] - b["ring"]) / np.maximum(b["ring"], 1e-300)
        walls_info[str(b["pid"])] = {
            "name": b["name"], "kind": b["kind"], "file": b["file"], "n_nodes": int(len(b["x"])), "heat": h,
            "x_range_m": [float(b["x"].min()), float(b["x"].max())], "x_monotone": b["x_monotone"],
            "polygon_lateral_area_m2": b["polygon_lateral_area"], "ring_area_sum_m2": float(np.sum(b["ring"])),
            "A_node_sum_m2": float(np.sum(b["A_node"])),
            "A_node_vs_ring_rel_max": float(rel_ring.max()), "A_node_vs_ring_rel_median": float(np.median(rel_ring)),
            "Q_total_W": (float(np.sum(b["Q"])) if valid else None),
            "Q_total_from_q_ring_W": (float(np.sum(b["q"] * b["ring"])) if valid else None),
            "y1plus_note": b["y1plus_note"],
            "y1plus_mean": (float(np.nanmean(b["y1plus"])) if np.isfinite(b["y1plus"]).any() else None),
            "y1plus_max": (float(np.nanmax(b["y1plus"])) if np.isfinite(b["y1plus"]).any() else None)}
        if valid and len(b["x"]) >= 3:
            # 端節点 (入口・出口・隣の境界と共有する節点) の q は隣の節点と大きく違うことがある (case/63 run_0013 の加熱壁で
            # −qwall の 3.8 倍・1.9 倍、2026-10-13 測定)。原因は未確定なので値は変えず、隣との比を記録する。
            walls_info[str(b["pid"])]["end_nodes"] = [
                {"x_m": float(b["x"][i]), "q_W_m2": float(b["q"][i]), "q_neighbor_W_m2": float(b["q"][j]),
                 "ratio_to_neighbor": float(b["q"][i] / b["q"][j]) if b["q"][j] != 0.0 else None,
                 "minus_qwall_W_m2": (float(b["qdiag"][i]) if b["qdiag"] is not None else None)}
                for i, j in ((0, 1), (len(b["x"]) - 1, len(b["x"]) - 2))]
    st = {"status": T["status"], "step": T["step"], "run_dir": str(T["run_dir"]), "mesh": str(rc["mesh"]),
          "csv": str(out_csv), "n_rows": len(data["x_m"]), "columns": list(cols),
          "heat_dataset": "iface_Qf_eff" if not use_diag else "qwall (診断値、符号を反転)",
          "q_column": qn, "Q_column": Qn, "allow_diagnostic": bool(allow_diagnostic),
          "geometry": ("軸対称 (isAxisymmetric 1, axisymMethod %d)" % rc["axisymMethod"]) if axi else "平面 2D",
          "units": {"x_m": "m", "r_m": "m", "s_m": "m", "ring_area_m2": "m²", "A_node_m2": "m²" if axi else "m²/m (奥行き 1 m あたり)",
                    "p_w_Pa": "Pa", "T_w_K": "K", "tau_w_Pa": "Pa", qn: "W/m²", Qn: "W" if axi else "W/m (奥行き 1 m あたり)",
                    "y1plus": "-"},
          "definitions": {
              "Q_node": ("2π · iface_Qf_eff (全周の節点熱量、固体向き正)" if axi else "iface_Qf_eff (固体向き正)") if not use_diag
              else "q_diag × A_node (診断値)",
              "q": "Q_node / A_node" if not use_diag else "−qwall (qwall は流体向き正の診断値)",
              "A_node": "2π max(r_f, r_floor) |S_f| (壁の境界面の r 重み面積 = ソルバの面)" if axi else "|S_f|",
              "ring_area": "π (r_a + r_b) |ab| を両端に半分ずつ (折れ線の側面積の持ち分)",
              "tau_w": "twall の折れ線の接線への射影 (流れ方向 +x を正)",
              "sign": "q・Q は流体 → 壁が正",
              "not_used": "iface_q_eff (main のソルバでは軸対称で r·q。ソルバの版で意味が変わる)"},
          "values_note": ("状態が valid でない壁の q・Q は空欄 (正本が無効なとき診断値に置き換えない)" if not use_diag
                          else "診断値 (−qwall)。熱設計の正本ではない"),
          "walls": walls_info}
    sp = out_csv.with_suffix(".status.json")
    sp.write_text(json.dumps(st, indent=1, ensure_ascii=False))
    st["status_json"] = str(sp)
    return st


HEAT_BANDS = 16   # 局所の熱量の定常性を見る区間の数 (壁ごと、x の順に節点数で等分; codex 諮問 2026-10-13 result の採用)


def wall_heat_series(run_dir, out_csv, wall_ids=None, n_bands: int = HEAT_BANDS) -> dict:
    """全 `res_<壁名>_<physID>_*.h5` から step ごとの熱量を時系列 CSV に (`check_quasisteady.py --series-csv` の入力): 列 `step`・
    `Q_total_W`・壁ごとの `Q_wall<physID>_W` と、**壁を x の順に節点数で `n_bands` 等分した区間ごとの熱量 `Qband<physID>_<k>_W`**
    (総熱量が一定のまま局所の熱流束が入れ替わる変化を見落とさないため; codex 諮問 2026-10-13 result)。正本 (`iface_Qf_eff`) が
    有限でない step・壁のダンプが欠けた step は空欄 (非有限として残す。落とさない)。"""
    import h5py
    run_dir = Path(run_dir)
    rc = _run_config(run_dir)
    walls = rc["walls"] if wall_ids is None else {int(w): rc["walls"][int(w)] for w in wall_ids}
    dumps = _wall_dumps(run_dir, walls)
    axi = rc["isAxisymmetric"] == 1
    ok_geom = not (axi and rc["axisymMethod"] != 0)
    steps = sorted(set().union(*[set(d) for d in dumps.values()])) if dumps else []
    pids = sorted(walls)
    rows, stat = [], {}
    nb, nn = {}, {}
    for st in steps:
        tot, per, bands = 0.0, [], []
        for pid in pids:
            p = dumps[pid].get(st)
            Q = np.nan
            bq = None
            if p is not None and ok_geom:
                with h5py.File(p, "r") as f:
                    a = f.attrs.get("step")
                    if a is not None and int(a) != st:
                        raise ValueError(f"{p.name}: 属性 step {int(a)} ≠ ファイル名の step {st}")
                    if "VALUE/iface_Qf_eff" in f:
                        qf = np.array(f["VALUE/iface_Qf_eff"], dtype=float)
                        xw = np.asarray(f["MESH/COORD"], dtype=float).reshape(-1, 3)[:, 0]
                        if np.all(np.isfinite(qf)):
                            fac = 2.0 * np.pi if axi else 1.0
                            Q = float(fac * qf.sum())
                            k = nb.setdefault(pid, max(1, min(int(n_bands), len(qf))))
                            o = np.argsort(xw, kind="stable")
                            bq = [float(fac * qf[ix].sum()) for ix in np.array_split(o, k)]
                            nn.setdefault(pid, len(qf))
                            bq = (bq, [float(fac * v) for v in qf[o]])          # 区間の和 (補助) と節点ごと (判定の本体)
            per.append(Q)
            tot += Q
            bands.append((pid, bq))
        rows.append((st, tot, per, bands))
        stat[st] = "valid" if np.isfinite(tot) else "invalid"
    band_cols = [f"Qband{pid}_{k:02d}_W" for pid in pids for k in range(nb.get(pid, 0))]
    node_cols = [f"Qnode{pid}_{i:05d}_W" for pid in pids for i in range(nn.get(pid, 0))]
    cols = ["step", "Q_total_W"] + [f"Q_wall{pid}_W" for pid in pids] + band_cols + node_cols
    flat = []
    for st, tot, per, bands in rows:
        bvals, nvals = [], []
        for pid, bq in bands:
            n, m = nb.get(pid, 0), nn.get(pid, 0)
            ok_ = bq is not None and len(bq[0]) == n and len(bq[1]) == m
            bvals += (bq[0] if ok_ else [np.nan] * n)
            nvals += (bq[1] if ok_ else [np.nan] * m)
        flat.append([st, tot] + per + bvals + nvals)
    rows = flat
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        f.write(",".join(cols) + "\n")
        for r_ in rows:
            f.write(",".join([str(int(r_[0]))] + [_fmt(v) for v in r_[1:]]) + "\n")
    n_valid = sum(1 for v in stat.values() if v == "valid")
    return {"csv": str(out_csv), "columns": cols, "band_columns": band_cols, "node_columns": node_cols, "n_bands": nb,
            "n_nodes": nn, "n_steps": len(steps),
            "n_valid": n_valid, "steps": steps, "geometry_ok": ok_geom, "unit": "W" if axi else "W/m", "per_step": stat}


# ---------------------------------------------------------------------------------------------------------------------
# 熱設計に使えるかの判定
# ---------------------------------------------------------------------------------------------------------------------

def _run_tool(args, python=None, timeout=1800):
    exe = python or sys.executable
    q = subprocess.run([exe, *[str(a) for a in args]], capture_output=True, text=True, timeout=timeout)
    return q.returncode, q.stdout + (("\n" + q.stderr) if q.stderr.strip() else "")


def _mirror_for_segment(run_dir: Path, td: Path) -> list:
    """`check_convergence.py --segment` は区間の連結 CSV を run ディレクトリに書くので、判定に要る入力だけを一時ディレクトリへ
    写して、そこで回す (他の run ディレクトリに何も書かない)。写すもの: stage_manifest.json・forge_launches.jsonl・
    residual_history.csv・manifest の各段の履歴と、`solverConfig.yaml` (`check_convergence.transition_active` が残差 CSV の隣から
    読む。写さないと遷移モデルの run で必須列 rms_roGamma/rms_roReth の検査が黙って外れる)・`bcondConfig.yaml`。"""
    man = json.loads((run_dir / "stage_manifest.json").read_text())
    names = {"stage_manifest.json", "forge_launches.jsonl", "residual_history.csv", "solverConfig.yaml", "bcondConfig.yaml"}
    stages = man.get("stages", man) if isinstance(man, dict) else man
    if isinstance(stages, list):
        for s_ in stages:
            if isinstance(s_, dict) and s_.get("history"):
                names.add(str(s_["history"]))
    copied = []
    for nm in sorted(names):
        src = run_dir / nm
        if src.is_file():
            dst = td / nm
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            copied.append(nm)
    return copied


def _judge_selection(run_dir: Path) -> dict:
    """収束の判定の道具の選択 (run_case.sh・stop_rule.postprocess・受け渡しの関門と同じ `settle.select_judge`)。"""
    if str(TOOLS) not in sys.path:
        sys.path.insert(0, str(TOOLS))
    import settle
    return settle.select_judge(run_dir)


def _settled_convergence(run_dir: Path, st: int, last_step: int, sel: dict, python=None) -> dict:
    """新方式 (settle_criteria.yaml のある run の本段・手で組んだ run) の収束: `check_settled.py <run> --bind --end-step <st> --no-write
    --json` (本段の起動・基準との照合の後に判定; run ディレクトリに書かない)。終了コード 0 → PASS、1 (NOT_SETTLED)・2 (DIVERGED) →
    FAIL、3 (INDETERMINATE: 照合の不一致・判定不能) ほか → 未評価。停止規則 v2 の再評価 (STOP_RULE_STATE.json の reeval) が同じ終端に
    あれば、量の判定 (digest)・verdict が一致することも確かめ、食い違えば未評価 (記録どうしが合わない場を熱設計に使わない)。
    verdict_line には道具の名前を入れる (成果物の xlsx の注記・status.json にどちらの道具かが残る)。"""
    conv = {"tool": "check_settled.py --bind", "selection": {k: sel.get(k) for k in ("tool", "mode", "why", "problems")}}
    if st != last_step:
        conv.update(result="未評価", reason=f"スナップショット step {st} が最終の壁ダンプ {last_step} でない (判定の終端を切り出さない)")
        return conv
    args = [TOOLS / "check_settled.py", run_dir, "--bind", "--no-write", "--json", "--end-step", st]
    rc_, txt = _run_tool(args, python)
    conv["cmd"] = " ".join(str(a) for a in args)
    conv["exit"] = rc_
    lines = txt.splitlines()
    vl = next((ln.strip() for ln in reversed(lines) if ln.startswith("VERDICT:")), None)
    res = None
    try:
        a = lines.index("{")
        b = max(i for i, ln in enumerate(lines) if ln == "}")
        res = json.loads("\n".join(lines[a:b + 1]))
    except ValueError:
        res = None
    if res is not None:
        conv["judged"] = {"end_step": res.get("end_step"), "terminal_step": res.get("terminal_step"), "stop_step": res.get("stop_step"),
                          "verdict": res.get("verdict"), "criteria_sha256": (res.get("criteria") or {}).get("sha256"),
                          "quantities": {str(q.get("name")): q.get("verdict") for q in res.get("quantities") or []}}
    where = (f"終端 {res.get('terminal_step')}" + (f"・停止の観測 {res.get('stop_step')}" if res.get("stop_step") is not None else "")
             if res is not None else f"end_step {st}")
    conv["verdict_line"] = f"{vl or 'VERDICT: (行なし)'} [check_settled.py --bind、{where}]"
    result = {0: "PASS", 1: "FAIL", 2: "FAIL"}.get(rc_, "未評価")
    if result == "未評価":
        conv["reason"] = vl or f"check_settled.py の終了コード {rc_} (VERDICT の行が無い)"
    # 停止規則 v2 の再評価との一致 (同じ終端のとき)
    try:
        from ..evaluate import stop_rule as SR
        stt = SR.read_state(run_dir)
    except Exception:  # noqa: BLE001 — 読めなければ照合しない (記録に残す)
        stt = None
    ve = (stt or {}).get("reeval") if isinstance(stt, dict) and stt.get("method") == "v2" else None
    if isinstance(ve, dict) and ve.get("end_step") == st:
        same = res is not None and SR.quantities_digest(res) == ve.get("digest") and res.get("verdict") == ve.get("verdict")
        conv["stop_rule_reeval"] = {"end_step": ve.get("end_step"), "verdict": ve.get("verdict"), "matches": bool(same),
                                    "state_status": stt.get("status")}
        if not same and result == "PASS":
            result = "未評価"
            conv["reason"] = (f"停止規則 v2 の再評価 (終端 {ve.get('end_step')}、{ve.get('verdict')}) と今の判定が食い違う "
                              "(記録どうしが合わない場を熱設計に使わない)")
    else:
        conv["stop_rule_reeval"] = {"compared": False, "why": ("停止規則 v2 の記録が無い" if not isinstance(ve, dict)
                                                               else f"再評価の終端 {ve.get('end_step')} が対象の step {st} と違う")}
    conv["result"] = result
    return conv


def thermal_use(run_dir, step=None, wall_ids=None, python=None) -> dict:
    """壁面量を熱設計に使えるか (plan §4.5): 対象 snapshot と同じ区間の `check_convergence.py` (`stage_manifest.json` があれば
    `--segment`)・総熱量の時系列に `check_quasisteady.py --series-csv`・`check_wall_resolution.py` の 3 つが合格し、熱量の状態が
    `valid` のときだけ `PASS`。不合格があれば `FAIL`、判定材料が無ければ `未評価`。各ツールの VERDICT 行を返す。
    ツールは run ディレクトリに書かせない (収束の --segment は入力を一時ディレクトリに写して回す)。
    収束の道具は `settle.select_judge` で選ぶ (plan tooling-settle-convergence-and-monitor §5.1 #7b): settle_criteria.yaml のある run の
    本段の記録のある run・手で組んだ run は `_settled_convergence` (check_settled.py --bind; 照合の不一致・起動の記録の欠け・前の段と確かめ
    られない run は判定不能)、それ以外 (基準の無い run・前の段と確かめられた run) は今どおり (plan §5.1 #7c)。
    基準のある run は判定不能でも check_convergence.py へ落とさない。"""
    run_dir = Path(run_dir)
    T = _build_table(run_dir, step, False, wall_ids, with_y1plus=False)
    st, heat = T["step"], T["status"]
    names = [b["name"] for b in T["blocks"]]
    last_step = max(max(d) for d in T["dumps"].values())
    out = {"run_dir": str(run_dir), "step": st, "walls": {str(b["pid"]): b["name"] for b in T["blocks"]},
           "heat_status": heat, "heat_per_wall": {str(k): v["status"] for k, v in T["per_wall"].items()}}
    # 1. 収束
    sel = _judge_selection(run_dir) if os.path.lexists(run_dir / "settle_criteria.yaml") else None
    if sel is not None and sel["tool"] == "check_settled.py":
        conv = _settled_convergence(run_dir, st, last_step, sel, python)
    else:
        conv = _legacy_convergence(run_dir, st, last_step, python)
        if sel is not None:
            conv["selection"] = {k: sel.get(k) for k in ("tool", "mode", "why")}
    out["convergence"] = conv
    _thermal_use_rest(out, run_dir, T, st, heat, names, python)
    return out


def _legacy_convergence(run_dir: Path, st, last_step, python=None) -> dict:
    """今までの収束 (`check_convergence.py`; stage_manifest.json があれば --segment を一時ディレクトリで)。"""
    conv = {"tool": "check_convergence.py"}
    if st != last_step:
        conv.update(result="未評価", reason=f"スナップショット step {st} が最終の壁ダンプ {last_step} でない (区間を切り出さない)")
    else:
        with tempfile.TemporaryDirectory(prefix="thermal_use_conv_") as td:
            if (run_dir / "stage_manifest.json").exists():
                copied = _mirror_for_segment(run_dir, Path(td))
                args = [TOOLS / "check_convergence.py", td, "--segment"]
                conv["mirrored"] = copied
            else:
                args = [TOOLS / "check_convergence.py", run_dir]
            rc_, txt = _run_tool(args, python)
        conv["cmd"] = " ".join(str(a) for a in args)
        conv["exit"] = rc_
        conv["verdict_line"] = next((ln.strip() for ln in txt.splitlines() if ln.startswith("===") and "->" in ln), None)
        conv["overall_line"] = next((ln.strip() for ln in txt.splitlines() if ln.startswith("OVERALL")), None)
        if rc_ == 0:
            conv["result"] = "PASS"
        elif any(t in txt for t in ("NO residual_history.csv", "判定不能", "empty residual file")):
            conv["result"] = "未評価"
            conv["reason"] = next((ln.strip() for ln in txt.splitlines()
                                   if "NO residual_history" in ln or "判定不能" in ln or "empty residual" in ln), None)
        else:
            conv["result"] = "FAIL"
    return conv


def _thermal_use_rest(out: dict, run_dir: Path, T: dict, st, heat, names, python=None) -> None:
    """thermal_use の 2. 総熱量の準定常・3. 壁解像・総合の判定 (out に書く; 収束は out["convergence"] に入っている)。"""
    conv = out["convergence"]
    # 2. 総熱量の準定常
    qs = {"tool": "check_quasisteady.py --series-csv"}
    if heat != "valid":
        qs.update(result="未評価", reason=f"熱量の状態が {heat} (正本の総熱量の時系列が作れない)")
    else:
        with tempfile.TemporaryDirectory(prefix="thermal_use_qs_") as td:
            ser = wall_heat_series(run_dir, Path(td) / "wall_heat_series.csv", [b["pid"] for b in T["blocks"]])
            # 総熱量と**節点ごとの熱量** (表と 1 次元の壁伝導が使う量; codex result レビュー 2 回目 M3: 区間の和は区間内の入れ替わりを
            # 相殺するので補助の診断に留め、判定は節点ごとの最悪で行う)
            qcols = ["Q_total_W"] + ser["node_columns"]
            args = [TOOLS / "check_quasisteady.py", "--series-csv", ser["csv"], "--series-cols", ",".join(qcols)]
            rc_, txt = _run_tool(args, python)
        qs.update(cmd=f"check_quasisteady.py --series-csv <tmp>/wall_heat_series.csv --series-cols Q_total_W,<{len(qcols) - 1} 節点の熱量>",
                  exit=rc_, n_bands=ser["n_bands"],
                  not_steady=[ln.strip() for ln in txt.splitlines() if ln.strip().startswith(("Q_total_W", "Qnode"))
                              and not ln.rstrip().endswith(" STEADY")][:8],
                  n_steps=ser["n_steps"], n_valid=ser["n_valid"],
                  verdict_line=next((ln.strip() for ln in txt.splitlines() if ln.startswith("===") and "->" in ln), None),
                  detail_line=next((ln.strip() for ln in txt.splitlines() if ln.strip().startswith("Q_total_W")), None))
        if rc_ == 0:
            qs["result"] = "PASS"
        elif "snapshot(s) (<" in txt or "ERROR" in txt:
            qs.update(result="未評価", reason="スナップショットが足りない (判定材料なし)")
        else:
            qs["result"] = "FAIL"
    out["quasisteady"] = qs
    # 3. 壁解像
    wr = {"tool": "check_wall_resolution.py"}
    args = [TOOLS / "check_wall_resolution.py", run_dir, "--step", st, "--groups", ",".join(names)]
    rc_, txt = _run_tool(args, python)
    wr.update(cmd=" ".join(str(a) for a in args), exit=rc_,
              verdict_line=next((ln.strip() for ln in txt.splitlines() if ln.startswith("VERDICT")), None))
    wr["result"] = {0: "PASS", 1: "FAIL"}.get(rc_, "未評価")
    out["wall_resolution"] = wr
    res = [conv["result"], qs["result"], wr["result"]]
    if "FAIL" in res or heat == "nonfinite":
        out["verdict"] = "FAIL"
    elif all(r_ == "PASS" for r_ in res) and heat == "valid":
        out["verdict"] = "PASS"
    else:
        out["verdict"] = "未評価"
    out["rule"] = ("収束・熱量の準定常 (総熱量と壁の全節点の熱量; 最悪で判定)・壁解像の 3 つが PASS かつ熱量の状態 valid で PASS、"
                   "不合格 (熱量 nonfinite を含む) があれば FAIL、他は未評価")


# ---------------------------------------------------------------------------------------------------------------------
# xlsx
# ---------------------------------------------------------------------------------------------------------------------

_XLSX_JOB = r'''
import csv, json, math, sys
job = json.load(open(sys.argv[1]))
def done(d):
    json.dump(d, open(job["result"], "w"), ensure_ascii=False)
try:
    import xlsxwriter
except ImportError as e:
    done({"written": False, "skipped": True, "reason": "xlsxwriter が無い: %s" % e}); sys.exit(0)
wb = xlsxwriter.Workbook(job["out"])
bold = wb.add_format({"bold": True})
inp = wb.add_format({"bg_color": "#FFF2CC", "border": 1})
sheets, info = [], {}
def put_csv(name, path):
    ws = wb.add_worksheet(name)
    rows = list(csv.reader(open(path)))
    ws.write_row(0, 0, rows[0], bold)
    for i, r in enumerate(rows[1:], 1):
        for j, v in enumerate(r):
            if v == "":
                continue
            try:
                f = float(v)
            except ValueError:
                ws.write_string(i, j, v); continue
            if math.isfinite(f):
                ws.write_number(i, j, f)
            else:
                ws.write_string(i, j, v)
    ws.freeze_panes(1, 0)
    sheets.append(name)
    info[name] = {"columns": rows[0], "n_rows": len(rows) - 1}
    return rows[0], len(rows) - 1
if job.get("shape_csv"):
    put_csv("shape", job["shape_csv"])
cfd_hdr = None
if job.get("cfd_csv"):
    cfd_hdr, n_cfd = put_csv("cfd_wall", job["cfd_csv"])
ws = wb.add_worksheet("conditions"); sheets.append("conditions")
ws.write_row(0, 0, ["key", "value"], bold)
for i, (k, v) in enumerate(job["conditions"], 1):
    ws.write_string(i, 0, k)
    if isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v):
        ws.write_number(i, 1, v)
    else:
        ws.write_string(i, 1, str(v))
ws = wb.add_worksheet("notes"); sheets.append("notes")
for i, t in enumerate(job["notes"]):
    ws.write_string(i, 0, t)
if cfd_hdr is not None:
    from xlsxwriter.utility import xl_col_to_name
    qcol = next((c for c in cfd_hdr if c.startswith("q_into_wall")), None)
    need = ["x_m", "r_m", "T_w_K"] + ([qcol] if qcol else [])
    miss = [c for c in need if c not in cfd_hdr]
    if miss or qcol is None:
        info["wall_1d"] = {"skipped": True, "reason": "cfd_wall の列が足りない: %s" % (miss or ["q_into_wall*"])}
    else:
        ws = wb.add_worksheet("wall_1d"); sheets.append("wall_1d")
        ws.write_string(0, 0, "t_m", bold); ws.write_string(0, 2, "肉厚 t [m] (入力)")
        ws.write_string(1, 0, "k_W_mK", bold); ws.write_string(1, 2, "壁の熱伝導率 k [W/(m K)] (入力)")
        for rr, key in ((0, "t_m"), (1, "k_W_mK")):
            v = job["wall_1d"].get(key)
            if v is None:
                ws.write_blank(rr, 1, None, inp)
            else:
                ws.write_number(rr, 1, float(v), inp)
        ws.write_string(2, 0, "q_column", bold); ws.write_string(2, 1, qcol)
        ws.write_string(3, 0, "等温壁の単純な抵抗の式 T_back = T_w - q t / k (q は cfd_wall の %s、流体→壁が正)。連成は cht_loop + solid_shell" % qcol)
        if qcol.endswith("_diag_W_m2"):
            ws.write_string(4, 0, "注意: q は診断値 (-qwall) で熱設計の正本ではない")
        tuv = job.get("thermal_use")
        if tuv is not None:
            ws.write_string(5, 0, "thermal_use (この熱量を熱設計に使えるか): %s%s" % (
                tuv, "" if tuv == "PASS" else " — PASS でないので、この表の熱量で壁の温度を決めない (収束・準定常・壁解像を conditions/notes で確認)"), bold)
        hdr = ["x_m", "r_m", "T_w_K", "q_W_m2", "T_back_K"]
        h0 = 6
        ws.write_row(h0, 0, hdr, bold)
        L = {c: xl_col_to_name(cfd_hdr.index(c)) for c in need}
        for i in range(n_cfd):
            src = i + 2                  # cfd_wall の行 (1 始まり、ヘッダの次)
            row = h0 + 1 + i
            ex = row + 1
            ws.write_formula(row, 0, "=cfd_wall!%s%d" % (L["x_m"], src))
            ws.write_formula(row, 1, "=cfd_wall!%s%d" % (L["r_m"], src))
            ws.write_formula(row, 2, "=cfd_wall!%s%d" % (L["T_w_K"], src))
            ws.write_formula(row, 3, '=IF(ISNUMBER(cfd_wall!%s%d),cfd_wall!%s%d,"")' % (L[qcol], src, L[qcol], src))
            ws.write_formula(row, 4, '=IF(AND(ISNUMBER($B$1),ISNUMBER($B$2),ISNUMBER(D%d),ISNUMBER(C%d)),IF($B$2>0,C%d-D%d*$B$1/$B$2,""),"")'
                             % (ex, ex, ex, ex))
        ws.freeze_panes(h0 + 1, 0)
        info["wall_1d"] = {"n_rows": n_cfd, "q_column": qcol}
wtj = job.get("wall_temperature")
if wtj:
    # 壁温分布の表 (run の正本の JSON から; 説明の行・区切りの点・元の表と実寸の表)
    ws = wb.add_worksheet("wall_temperature"); sheets.append("wall_temperature")
    r0 = 0
    for t in wtj["lines"]:
        ws.write_string(r0, 0, t); r0 += 1
    r0 += 1
    ws.write_row(r0, 0, ["anchor", "x_m"], bold); r0 += 1
    for k, v in wtj["anchors"]:
        ws.write_string(r0, 0, k); ws.write_number(r0, 1, float(v)); r0 += 1
    r0 += 1
    ws.write_row(r0, 0, wtj["columns"], bold); h1 = r0; r0 += 1
    for row_ in wtj["rows"]:
        for j, v in enumerate(row_):
            ws.write_number(r0, j, float(v))
        r0 += 1
    ws.freeze_panes(h1 + 1, 0)
    info["wall_temperature"] = {"n_rows": len(wtj["rows"]), "columns": wtj["columns"]}
wcj = job.get("wall_compare")
if wcj:
    # 設計壁と物理壁の比較の表 (CSV のまま) と、表の右に注記 (区間・定義域・空欄の意味)
    hdr_c, _ = put_csv("wall_compare", wcj["csv"])
    ws = wb.get_worksheet_by_name("wall_compare")
    nc = len(hdr_c) + 1
    ws.write_string(0, nc, "注記", bold)
    for i, t in enumerate(wcj["notes"], 1):
        ws.write_string(i, nc, t)
    info["wall_compare"]["notes_column"] = nc
    info["wall_compare"]["n_notes"] = len(wcj["notes"])
wb.close()
done({"written": True, "skipped": False, "sheets": sheets, "info": info, "xlsxwriter": xlsxwriter.__version__})
'''


def _flatten(d, pre=""):
    out = []
    for k, v in (d or {}).items():
        key = f"{pre}{k}"
        if isinstance(v, dict):
            out.extend(_flatten(v, key + "."))
        elif isinstance(v, (list, tuple)):
            out.append((key, json.dumps(v, ensure_ascii=False, default=float)))
        elif isinstance(v, (bool, np.bool_)):
            out.append((key, str(bool(v))))
        elif isinstance(v, (int, float, np.integer, np.floating)):
            out.append((key, float(v)))
        else:
            out.append((key, "" if v is None else str(v)))
    return out


def write_xlsx(out_xlsx, shape_csv=None, cfd_csv=None, conditions: dict | None = None, notes: list | None = None,
               python: str | None = None, wall_1d: dict | None = None, timeout: int = 600, thermal_use: str | None = None,
               wall_temperature: dict | None = None, wall_compare: dict | None = None) -> dict:
    """形状表・CFD 壁面量表・条件・注記を `.xlsx` に書く (plan §4.5)。xlsxwriter は design の venv に無いので
    `.venv-pptx` の python (`python`、既定 `PPTX_PY`) を subprocess で呼ぶ。python か xlsxwriter が無ければ書かずに
    `{"written": False, "skipped": True, "reason": ...}` を返す (呼び出し側が飛ばしたことを MANIFEST・README に書く)。
    シート: shape・cfd_wall・conditions・notes と、cfd 表があれば 1 次元の壁伝導の記入用 wall_1d (入力セル t_m・k_W_mK、
    各壁節点の T_back_K = T_w − q t / k)。`wall_1d` = {"t_m": .., "k_W_mK": ..} で入力セルの初期値 (既定は空欄)。
    `wall_temperature` (壁温分布の run だけ; plan tooling-nozzle-wall-temperature-profile §4.4) = {"lines": [説明の行],
    "anchors": [[名前, x_m]], "columns": [列名], "rows": [[数値]]} を渡すと最後に専用シート `wall_temperature` を足す
    (None なら今までどおり; ほかのシートの並びは変えない)。
    `wall_compare` (設計壁と物理壁の比較; plan tooling-nozzle-card-local-mesh §5.1 #3u) = {"csv": 比較の表の CSV, "notes": [注記の行]}
    を渡すと最後にシート `wall_compare` (表と、その右の列に注記) を足す (None なら今までどおり)。"""
    exe = python or PPTX_PY
    if not Path(exe).exists():
        return {"written": False, "skipped": True, "reason": f"python が無い: {exe}"}
    nts = list(notes or []) + ["wall_1d: 等温壁の単純な抵抗の式 T_back = T_w - q t / k。連成が要るなら cht_loop + solid_shell を使う"]
    out_xlsx = Path(out_xlsx)
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="wall_xlsx_") as td:
        jp, rp, sp = Path(td) / "job.json", Path(td) / "result.json", Path(td) / "xlsx_job.py"
        sp.write_text(_XLSX_JOB)
        job = {"out": str(out_xlsx.resolve()), "result": str(rp),
               "shape_csv": str(Path(shape_csv).resolve()) if shape_csv else None,
               "cfd_csv": str(Path(cfd_csv).resolve()) if cfd_csv else None,
               "conditions": _flatten(conditions), "notes": nts,
               "wall_1d": dict(wall_1d or {}), "thermal_use": thermal_use}
        if wall_temperature is not None:
            job["wall_temperature"] = wall_temperature
        if wall_compare is not None:
            job["wall_compare"] = {"csv": str(Path(wall_compare["csv"]).resolve()), "notes": [str(t) for t in wall_compare["notes"]]}
        jp.write_text(json.dumps(job, ensure_ascii=False))
        q = subprocess.run([exe, "-I", str(sp), str(jp)], capture_output=True, text=True, timeout=timeout, cwd=td)
        if not rp.exists():
            raise RuntimeError(f"xlsx の作業が結果を返さない (rc={q.returncode}):\n{q.stdout[-2000:]}\n{q.stderr[-2000:]}")
        res = json.loads(rp.read_text())
    res["python"] = exe
    if res.get("written"):
        res["xlsx"] = str(out_xlsx)
    return res


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a1 = sub.add_parser("shape")
    a1.add_argument("wall")
    a1.add_argument("out")
    for nm in ("cfd", "series", "thermal"):
        a = sub.add_parser(nm)
        a.add_argument("run_dir")
        if nm != "thermal":
            a.add_argument("out")
        if nm != "series":
            a.add_argument("--step", type=int, default=None)
        a.add_argument("--walls", default=None, help="physID をカンマ区切り (既定は全部の壁)")
        if nm == "cfd":
            a.add_argument("--allow-diagnostic", action="store_true")
    a = ap.parse_args(argv)
    walls = [int(v) for v in a.walls.split(",")] if getattr(a, "walls", None) else None
    if a.cmd == "shape":
        r = shape_table(a.wall, a.out)
    elif a.cmd == "cfd":
        r = cfd_wall_table(a.run_dir, a.out, a.step, a.allow_diagnostic, walls)
    elif a.cmd == "series":
        r = wall_heat_series(a.run_dir, a.out, walls)
    else:
        r = thermal_use(a.run_dir, a.step, walls)
    print(json.dumps(r, indent=1, ensure_ascii=False, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
