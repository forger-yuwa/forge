r"""axis-Mach チェーン A4: Hall + 5次 Hermite 軸 Mach → 逆 MOC → node Euler 評価。

plan: plans/accepted/tooling-nozzle-axismach-chain.md §6 A4。既存経路 (runner_wt の
モード F / runner_walldriven) には触れない。

チェーン: HallThroat (初期値線 + 軸アンカー) → QuinticHermiteAxisLaw
(自由度 L_c のみ) → inverse_design (E→F 閉包込み) → wall_qa →
AxisMachCFDWall (直管 + U→T Hermite + 設計壁 spline。旧・縦線構成のみ
[T, x0] に骨接放物線が入る) →
TFI メッシュ → node Euler (段階起動)。

CFD-in-the-loop アンカー更新 (A5) は problem YAML の geometry キーで受ける:
  x_reach_cfd:     壁始点発 C⁻ の軸着地点 (node 基準 run から)
  x_reach_anchor:  [M, M', M''] (同 run の軸平滑化フィットから)
  axis_segment_run: [x0, x_reach) の実測軸 M を target_moc に渡す元 run
これらが無ければ初回 (Hall アンカー、x_A = x0)。
初期値線は geometry.start_line で選ぶ ('throat_char' = スロート特性線 [A8] /
'vertical' = M_start の縦線 [旧構成])。壁の決め方は geometry.wall_mode
('flux' = 断面の質量流束閉包 [A9] / 'streamline' = 流線積分 [旧構成])。

格子: Euler (`prepare`) は問題 YAML の `mesh_euler`、NS (`prepare_ns`) は `mesh` を読む (2026-10-07 から。混ぜて補完しない。
`mesh_euler` の無い問題は Euler の prepare で移行先を示して止まる; `mesh_euler_block`)。

使い方:
  design/.venv-opt/bin/python -m forge_design.evaluate.runner_axismach \
      case/45.isobutane_m6_d155/problem_d155_euler_pin_G1_recal_mono.yaml run_dir [--prepare-only]
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

from ..geometry.axis_law import QuinticHermiteAxisLaw, KnotQuinticAxisLaw
from ..geometry.moc_inverse import inverse_design
from ..geometry.transonic import HallThroat
from ..geometry.wall_axismach import (PHYSICAL_WALL_REPRS, PW_UPSTREAM_DEFAULT, PW_UPSTREAMS, AxisMachCFDWall,
                                      area_ratio_isentropic, wall_qa)
from ..meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, node_count, radial_node_count, write_msh41_2d
from ..probdef import Problem, dv_value, load_problem
from ..wall_temperature import check_covers, integral_thermal_bc, load_run_profile, to_x_m, write_run_files
from ..wall_temperature import anchors_from_design as _wt_anchors
from .ic import paste_isentropic_ic, stamp_isentropic_ic_species
from .runner import run_converter, FORGE_TOOLS, PROBE_STUB, _ENV, converter_path, run_forge
from .runner_wt import (_bcond, _config_euler, _config_euler_node,
                        _config_sst_node)


def axis_curve_node(run_dir, scale: float, lam: float = 1e-5, mode: str = "evenfit"):
    r"""node run の軸列から平滑化スプライン M(x) を作る (§15 の実装)。

    mode="evenfit" (既定): 軸ノードを除外し r>0 の 4 点で M=a₀+a₂r² を最小二乗した a₀ を軸値に
    使う (`axis_extract.axis_curve_evenfit`)。旧生産設定 (`nodeAxisDirichlet: 1`, 撤去済) の軸ノードは第一内点の
    コピーで ½M_rr r₁² のバイアス (case/41 で ~0.08% M_d) を持っていたため。現行 (軸 DOF) でも外挿の方が精度が良い。
    mode="node": 従来の軸ノード直読。末尾 3 スナップ平均 → `make_smoothing_spline`。
    M, M', M'' は**同一スプライン**から評価する (生の 2 階差分禁止)。
    戻り値: (spline, x_min, x_max) — x は r_t 単位。"""
    if mode == "evenfit":
        from .axis_extract import axis_curve_evenfit
        return axis_curve_evenfit(run_dir, scale, lam=lam)
    import h5py
    from scipy.interpolate import make_smoothing_spline
    rd = Path(run_dir)
    res = sorted(rd.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    res = [f for f in res if int("".join(c for c in f.stem if c.isdigit())) > 0]
    if not res:
        raise FileNotFoundError(f"{rd} に res_*.h5 が無い")
    with h5py.File(rd / "nozzle.h5") as nz:
        nc = nz["/MESH/COORD"][:].reshape(-1, 3)
    Ms = []
    for rf in res[-3:]:
        with h5py.File(rf) as f:
            Ux, Uy, son = f["/VALUE/Ux"][:], f["/VALUE/Uy"][:], f["/VALUE/sonic"][:]
        if len(Ux) != len(nc):
            raise ValueError("node run でない (VALUE 長 != 節点数)")
        Ms.append(np.hypot(Ux, Uy) / np.maximum(son, 1e-9))
    M = np.mean(Ms, axis=0)
    ax = nc[:, 1] < 1e-12 * max(float(nc[:, 1].max()), 1.0)
    x = nc[ax, 0] / scale
    o = np.argsort(x)
    x, Ma = x[o], M[ax][o]
    spl = make_smoothing_spline(x, Ma, lam=lam)
    return spl, float(x.min()), float(x.max())


def _gam_or_gas(p: Problem):
    """MOC/幾何関数に渡す γ: cpg なら float、semiperfect ならガスモデル
    (cfd_gas: cpg のときは設計も CPG に落とす — 熱力学の一致)。"""
    if str(p.evaluate.get("cfd_gas", "same")) == "cpg":
        return p.gamma
    return p.gas_model if p.is_semiperfect else p.gamma


def _apply_gas_to_config(cfg: str, p: Problem, run_dir, viscous: bool = False) -> str:
    """semi-perfect のとき forge config を TP (thermalMethod 2) に書き換える。cpg なら無変更。
    `physProp.species` には種名と lump の構成 (構成種と全桁の lump 内モル分率, basis: mole) だけを書き、lump の係数は
    ソルバが起動時に合成する (plan thermophysics-solver-owned-species-db §4.7 #9; 合成済み擬似種の species_db.yaml は作らない)。
    ソルバ内蔵で解決できない実種 (外部 DB `gas.species_db` 由来など) があるときだけ、その生エントリを
    species_db_external.yaml に置いて speciesDBFile で渡す。設計 (MOC) と CFD は同じ共通データの係数を使う。
    `viscous=True` (NS/SST の config) では種ごとの輸送物性を使う: `viscMethod: 1` (空気の Sutherland) を `viscMethod: 2` に、
    `gas.transport` を実種ごとの `physProp.transport` にする (plan thermophysics-solver-owned-species-db #9b; `gas.transport` 必須)。
    `thermCondMethod` は viscMethod 2 では読まれないので落とす。`visc` (dt と陰解法対角の剛性見積り) と必須キーの `thermCond`、
    `prandtlLam` (SST 壁関数の回復係数) は残す。Euler (`viscous=False`) と CPG の config は従来と同じ。"""
    # 切り分け用: evaluate.axisym_method で CPG でも SU2 流軸対称に切替可
    if int(p.evaluate.get("axisym_method", 0)) == 1:
        cfg = cfg.replace("isAxisymmetric: 1", "isAxisymmetric: 1, axisymMethod: 1", 1)
    if not p.is_semiperfect or str(p.evaluate.get("cfd_gas", "same")) == "cpg":
        # cfd_gas: cpg = 設計は semi-perfect のまま CFD だけ CPG(γ*, cp 参照値) で回す
        # (TP × node 軸対称の forge 側発散 [case/42 run_0001] の回避。相対比較には十分)
        return cfg
    from ..gas.composition import (physprop_species_flow, solver_species_config, species_db_raw_yaml,
                                   write_species_meta)
    # 統一スキーマ (plan thermophysics-cea-mole-fraction-species §4.5): evaluate.tp_species {mode: full|lumped, lumps, keep}
    # (旧 pseudo / split_h2o は別名) を解決済み DB で輸送種配置に解決し、config の species (lump 記法) と species_meta.yaml を書く
    layout = p.species_layout()
    species_list = list(layout.species)
    items, external = solver_species_config(layout)
    transport = p.transport_for_ns(layout) if viscous else None
    write_species_meta(layout, run_dir, transport)
    db_key = ""
    if external:
        (Path(run_dir) / "species_db_external.yaml").write_text(species_db_raw_yaml(external))
        db_key = ', speciesDBFile: "species_db_external.yaml"'
    # thermalMethod 0 → 2、species (/speciesDBFile) を physProp に追加 (cp/gamma は参照値のまま
    # 残すが TP では NASA-9 が優先される)
    cfg = cfg.replace("thermalMethod: 0", "thermalMethod: 2", 1)
    # [2026-08-16] nodeAxisDirichlet は撤去済み (node は軸ノードを DOF として解く整合セットが常時 ON、
    # TP スカラーの軸ピン問題は消滅)。axisRFloor は evaluate.axis_r_floor 指定時のみ従来どおり付ける。
    floor = float(p.evaluate.get("axis_r_floor", 0.0))
    if floor > 0:
        cfg = cfg.replace("isAxisymmetric: 1", f"isAxisymmetric: 1, axisRFloor: {floor}", 1)
    # thermoHrefTemp (sensible-enthalpy datum): 絶対基準 (生成エンタルピー込み) のままだと
    # 陰解法 Jacobian の χ_eos = c² − κh が桁違いになり block-DPLUR が軸近傍で発散する
    # (case/42 run_0020–0025 で切り分け: 一定 cp 種/陽解法は完走、実 NASA-9 + 陰解法だけ発散、
    #  thermoHrefTemp 298.15 で完走)。IC の roe も同じ datum で作る (paste_isentropic_ic の h_ref)。
    href = float(p.evaluate.get("thermo_href_temp", 298.15))
    sp_txt = physprop_species_flow(items)   # 引用符付き: NO/N/Y は無引用だと YAML 1.1 で真偽値になる (codex result M1)
    cfg = cfg.replace("cp: %s, gamma: %s}" % (p.cp, p.gamma),
                      "cp: %s, gamma: %s,\n           species: %s%s, thermoHrefTemp: %s}"
                      % (p.cp, p.gamma, sp_txt, db_key, href), 1)
    if f"species: {sp_txt}" not in cfg:
        raise RuntimeError("_apply_gas_to_config: physProp の書き換えに失敗 (テンプレート変更?)")
    if transport is not None:
        from ..gas.composition import physprop_transport_flow
        n_vm = cfg.count("viscMethod: 1,")
        n_tcm = cfg.count(", thermCondMethod: 1")
        if n_vm != 1 or n_tcm > 1:
            raise RuntimeError(f"_apply_gas_to_config: NS の physProp が想定外 (viscMethod: 1 が {n_vm} 個; テンプレート変更?)")
        cfg = cfg.replace("viscMethod: 1,", "viscMethod: 2,", 1).replace(", thermCondMethod: 1", "", 1)
        tr_txt = physprop_transport_flow(transport)
        cfg = cfg.replace(f"species: {sp_txt}{db_key}, ", f"species: {sp_txt}{db_key},\n           transport: {tr_txt}, ", 1)
        if f"transport: {tr_txt}" not in cfg:
            raise RuntimeError("_apply_gas_to_config: physProp.transport の書き込みに失敗 (テンプレート変更?)")
    # 凝縮 (evaluate.condensation: dict) — forge の condensation ブロックをそのまま通す
    cond = p.evaluate.get("condensation")
    if cond:
        cond = dict(cond)
        # 凝縮種は名前が正本 (gas.condensing_species): index は生成値。手書き condGasSpecies があれば一致検査
        idx = layout.cond_index
        if idx is None:
            raise ValueError(f"凝縮 ON だが凝縮種 {layout.condensing_species} が輸送種 {species_list} に無い (tp_species.keep に入れる)")
        if "condGasSpecies" in cond and int(cond["condGasSpecies"]) != idx:
            raise ValueError(f"evaluate.condensation.condGasSpecies {cond['condGasSpecies']} が凝縮種 {layout.condensing_species} の index {idx} と不一致")
        cond["condGasSpecies"] = idx
        cond["condensationSpecies"] = layout.condensing_species
        cfg = cfg.rstrip("\n") + "\ncondensation: {" + ", ".join(f"{k}: {v}" for k, v in cond.items()) + "}\n"
    return cfg


def _tp_species_list(p: Problem) -> list:
    """TP の species リスト (統一スキーマ `evaluate.tp_species` を解決した輸送種順序)。"""
    return list(p.species_layout().species)


def _tp_species_Y(p: Problem):
    """IC/BC 用の組成 [Y_s] (species リスト順)。単一輸送種 (旧 pseudo = MIX) や CPG なら None (roY 不要)。"""
    if not p.is_semiperfect or str(p.evaluate.get("cfd_gas", "same")) == "cpg":
        return None
    layout = p.species_layout()
    if layout.n == 1:
        return None
    return layout.Y_transport("inflow")


def _stamp_ic_species(p: Problem, run_dir) -> str | None:
    """新規初期場 (paste_isentropic_ic) に化学種の属性を付ける (TP のときだけ; plan thermophysics-solver-owned-species-db
    §4.3 #3b)。IC と同じ gas・datum・輸送種の順序と MW を宛先の `forge --resolve-species` の記録と照合し、違えば例外で止める。"""
    if not p.is_semiperfect or str(p.evaluate.get("cfd_gas", "same")) == "cpg":
        return None
    layout = p.species_layout()
    species = list(layout.species)
    return stamp_isentropic_ic_species(Path(run_dir) / "nozzle.h5", run_dir, p.gas_model,
                                       float(p.evaluate.get('thermo_href_temp', 298.15)), species,
                                       [float(layout.entries[s].MW) for s in species], _tp_species_Y(p))


def _restart_same_mesh(res_h5, mesh_h5) -> None:
    """同一メッシュの段間引き継ぎ: `restart_field.py` (保存量の index コピー、SRC とビット一致を検査; 化学種の属性を継承)。
    `interp_field.py` (cross-mesh 用) は原始量から保存量を組み直すので同一メッシュには使わない (AGENTS.md「メッシュ変更後の restart」)。"""
    r = subprocess.run([sys.executable, str(FORGE_TOOLS / "restart_field.py"), str(res_h5), str(mesh_h5)],
                       env=_ENV, capture_output=True, text=True)
    with (Path(mesh_h5).parent / "restart_field.log").open("a") as f:
        f.write(r.stdout + r.stderr)
    if r.returncode != 0:
        raise RuntimeError(f"restart_field.py が失敗 ({res_h5} -> {mesh_h5}):\n{r.stdout[-2000:]}{r.stderr[-2000:]}")


def _species_info(p: Problem) -> dict | None:
    """prepare_info.json 用: 組成 (X / Y / 入力総和 / DB の出典) と輸送種配置の要約 (plan §4.1)。"""
    if not p.is_semiperfect:
        return None
    Y, X, tot = p.gas_composition
    out = {"composition_basis": p.composition_basis, "Y": Y, "X": X, "input_sum": tot,
           "species_db": p.raw.get("gas", {}).get("species_db"), "condensing_species": p.condensing_species}
    if str(p.evaluate.get("cfd_gas", "same")) != "cpg":
        L = p.species_layout()
        out.update({"mode": L.mode, "transported": list(L.species), "keep": list(L.keep), "condensing_index": L.cond_index,
                    "Y_transport": L.Y_transport("inflow"), "db_source": {s: L.entries[s].source for s in L.species}})
    return out


def _bcond_with_species(txt: str, Ys) -> str:
    """bcond の inlet floats に Y0.. を追記 (多成分 TP の入口組成)。"""
    if Ys is None:
        return txt
    ins = ", ".join(f"Y{i}: {y:.8f}" for i, y in enumerate(Ys))
    return txt.replace("floats: {Pt:", "floats: {" + ins + ", Pt:", 1)



# 問題 YAML のキーの選択肢。先頭が既定: 2026-10-07 ユーザ決定 (生産採用) で analytic・converge を既定にした (旧方式は legacy・fixed2 を
# 明示して再現する。plan discretization-moc-axis-limit-and-corrector §9)。MOC の関数 (moc_kernel・moc_inverse) の引数の既定は legacy・fixed2 のまま
_MOC_KEYS = {"moc_axis_limit": ("analytic", "legacy"), "moc_corrector": ("converge", "fixed2")}


def _moc_keys(geometry: dict) -> tuple:
    """geometry.moc_axis_limit / geometry.moc_corrector を読む (plans/accepted/discretization-moc-axis-limit-and-corrector.md
    §4.3)。キーが無ければ既定 (analytic, converge — 2026-10-07 から。それ以前は legacy, fixed2)。値は文字列で選択肢のどれかに完全一致すること — null・大文字違い・
    前後の空白・数値・真偽値は既定に読み替えず例外にする (黙って既定で設計しない)。"""
    out = []
    for key, choices in _MOC_KEYS.items():
        if key not in geometry:
            out.append(choices[0])
            continue
        v = geometry[key]
        if not isinstance(v, str) or v not in choices:
            raise ValueError(f"geometry.{key} は {' / '.join(repr(c) for c in choices)} のどれか "
                             f"(受け取った値: {v!r}。既定にするならキーを書かない)")
        out.append(v)
    return tuple(out)


def _physical_wall_repr(geometry: dict):
    """geometry.physical_wall_repr を読む (plans/accepted/tooling-nozzle-wall-single-bspline.md §4.4)。明示の値だけを返し、キーが無ければ None。
    キー無しの実効の表現は `build_physical_wall` が決める: poly の壁は 1 本で表せれば `single_bspline` (2026-10-07 ユーザ決定)、
    それ以外は今の区分表現 (壁ファイルを書かない)。値は 'legacy' / 'single_bspline' に完全一致
    すること — null・大文字違い・前後の空白・数値・真偽値は既定に読み替えず例外にする。"""
    if "physical_wall_repr" not in geometry:
        return None
    v = geometry["physical_wall_repr"]
    if not isinstance(v, str) or v not in PHYSICAL_WALL_REPRS:
        raise ValueError(f"geometry.physical_wall_repr は {' / '.join(repr(c) for c in PHYSICAL_WALL_REPRS)} のどれか "
                         f"(受け取った値: {v!r}。既定にするならキーを書かない)")
    return v


def _pw_upstream(geometry: dict) -> dict:
    """geometry.pw_upstream を読んで検査する (plans/accepted/tooling-nozzle-upstream-poly-and-throat-sizing.md §4.1)。重い処理・
    run dir を作る前に呼ぶ。戻り: {"value": 物理壁に渡す値 (None = 解析経路でない), "source": "explicit" | "default" | None,
    "requested": 書かれた値 (キー無しは None)}。

    - joint 壁 (`wall_repr: joint`) の物理壁は解析経路で、キー無しは既定 `poly` (ユーザ決定 2026-10-07)。`ramp` は明示したときだけ。
    - 値は 'ramp' / 'poly' に完全一致すること — null・大文字違い・前後の空白・数値・真偽値は例外 (既定に読み替えない)。
    - `poly` (明示・既定とも) と `pw_ramp` の併記は例外 (`pw_ramp` のキーがあれば値が null でも併記とみなす)。
    - 解析経路でない (joint でない) 壁に `poly` を明示したら例外。キー無しは今の振る舞いのまま (上流は従来経路)。
    - `physical_wall_repr: single_bspline` は `poly` だけ (`ramp` との組み合わせは例外; §4.1b)。"""
    joint = str(geometry.get("wall_repr", "interp")) == "joint"
    if "pw_upstream" in geometry:
        v = geometry["pw_upstream"]
        if not isinstance(v, str) or v not in PW_UPSTREAMS:
            raise ValueError(f"geometry.pw_upstream は {' / '.join(repr(c) for c in PW_UPSTREAMS)} のどれか "
                             f"(受け取った値: {v!r}。既定 '{PW_UPSTREAM_DEFAULT}' にするならキーを書かない)")
        requested, source = v, "explicit"
    else:
        requested, source = None, ("default" if joint else None)
    if not joint:
        if requested == "poly":
            raise ValueError("geometry.pw_upstream: poly は joint 壁 (wall_repr: joint) の物理壁の解析経路専用 "
                             f"(wall_repr = {geometry.get('wall_repr', 'interp')!r})")
        return {"value": None, "source": (None if requested is None else "explicit (解析経路でないので無効)"),
                "requested": requested}
    value = requested if requested is not None else PW_UPSTREAM_DEFAULT
    if value == "poly" and "pw_ramp" in geometry:
        raise ValueError(f"geometry.pw_upstream: poly ({'明示' if requested else '既定'}) と geometry.pw_ramp "
                         f"{geometry['pw_ramp']!r} の併記は不可 — 旧来のランプを使うなら pw_upstream: ramp を明示する、"
                         "poly にするなら pw_ramp を消す (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1)")
    if value == "ramp" and _physical_wall_repr(geometry) == "single_bspline":
        raise ValueError("geometry.physical_wall_repr: single_bspline は pw_upstream: poly の壁だけ (ramp のランプ区間は 5 次の "
                         "B-spline で厳密に表せない; plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1b)")
    return {"value": value, "source": source, "requested": requested}


# 物理壁だけの上流の形 `geometry.physical_upstream` のキー (geometry の r_inlet・L_U・L_pipe と同じ名前・同じ単位 r_t)
PHYSICAL_UPSTREAM_KEYS = ("r_inlet", "L_U", "L_pipe")


def _physical_upstream(geometry: dict) -> dict | None:
    """geometry.physical_upstream を読んで検査する (plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #34)。キー無しは None
    (今と同じ: 物理壁の直管・多項式 Q の始点は geometry の r_inlet・L_U・L_pipe)。重い処理・run dir を作る前に呼ぶ。

    `{r_inlet, L_U, L_pipe}` (r_t 単位、3 つとも必須) は**物理壁だけ**に効く上流の形: 直管の半径 r_inlet ([x_in, −L_U))・
    多項式 Q の始点 (−L_U, r_inlet, 0, 0) (Q の終点は今と同じ x = 0 の S + δ_r の (r, r′, r″))。geometry の r_inlet・L_U・L_pipe は
    設計壁 (MOC・CFD ピンの初期線の凍結元との照合、`design_chain`) の値のまま残す — 照合を弱めず、初期線を取り直さない。
    条件 (違えば ValueError): mapping でキーがちょうど 3 つ・有限の数値 (真偽値は不可)・r_inlet > 1・L_U > 0・L_pipe > 0、
    joint 壁 (wall_repr: joint) の pw_upstream poly (既定・明示とも; ramp は不可)、入口の位置 −(L_U + L_pipe) が設計の
    −(geometry.L_U + geometry.L_pipe) と同じ (格子の区間・壁ファイルの定義域を変えない)、名目の U→T Hermite の
    μ = L_U²/(R Δr) ≤ 20 (設計壁と同じ単調収縮の条件、R は geometry.R)。
    戻り: {"r_inlet", "L_U", "L_pipe", "x_in", "mu", "frozen_design": {r_inlet, L_U, L_pipe}, "source": "geometry.physical_upstream"}。"""
    if "physical_upstream" not in geometry:
        return None
    v = geometry["physical_upstream"]
    if not isinstance(v, dict):
        raise ValueError(f"geometry.physical_upstream は {{r_inlet, L_U, L_pipe}} の mapping (受け取った値: {v!r})")
    unknown = sorted((str(k) for k in v if k not in PHYSICAL_UPSTREAM_KEYS))
    missing = [k for k in PHYSICAL_UPSTREAM_KEYS if k not in v]
    if unknown or missing:
        raise ValueError(f"geometry.physical_upstream のキーが不正 (未知 {unknown}、欠け {missing}; {list(PHYSICAL_UPSTREAM_KEYS)} の 3 つとも書く)")
    vals = {}
    for k in PHYSICAL_UPSTREAM_KEYS:
        x = v[k]
        if isinstance(x, bool) or not isinstance(x, (int, float)) or not np.isfinite(float(x)):
            raise ValueError(f"geometry.physical_upstream.{k} は有限の数値 (受け取った値: {x!r})")
        vals[k] = float(x)
    if not (vals["r_inlet"] > 1.0 and vals["L_U"] > 0.0 and vals["L_pipe"] > 0.0):
        raise ValueError(f"geometry.physical_upstream: r_inlet > 1・L_U > 0・L_pipe > 0 であること ({vals})")
    if str(geometry.get("wall_repr", "interp")) != "joint":
        raise ValueError("geometry.physical_upstream は joint 壁 (wall_repr: joint) の物理壁の解析経路専用 "
                         f"(wall_repr = {geometry.get('wall_repr', 'interp')!r})")
    if _pw_upstream(geometry)["value"] != "poly":
        raise ValueError("geometry.physical_upstream は pw_upstream: poly の物理壁だけ (ramp は設計壁の縮流部に δ_r を足す作り方なので、"
                         "別の形の上流には当てられない)")
    fd = {"r_inlet": float(geometry.get("r_inlet", 2.5)), "L_U": float(geometry.get("L_U", 3.5)),
          "L_pipe": float(geometry.get("L_pipe", 0.5))}
    x_in_d = -fd["L_U"] - fd["L_pipe"]
    x_in = -vals["L_U"] - vals["L_pipe"]
    if abs(x_in - x_in_d) > 1e-12 * max(1.0, abs(x_in_d)):
        raise ValueError(f"geometry.physical_upstream: 入口の位置 −(L_U + L_pipe) = {x_in!r} が設計の −(geometry.L_U + geometry.L_pipe) "
                         f"= {x_in_d!r} と違う (入口の位置を動かすのは未対応: 格子の区間・壁ファイルの定義域を設計壁と揃える)")
    from ..geometry.wall_walldriven import UpstreamThroatPoly
    up = UpstreamThroatPoly(r_U=vals["r_inlet"], R_t=float(geometry.get("R", 2.0)), L_U=vals["L_U"])
    bad = up.validate()
    if bad:
        raise ValueError("geometry.physical_upstream の名目の U→T Hermite が設計の条件を満たさない: " + "; ".join(bad))
    return {**vals, "x_in": x_in_d, "mu": float(up.mu), "frozen_design": fd, "source": "geometry.physical_upstream"}


def _sizing_spec(spec: dict) -> dict | None:
    """spec.sizing (任意) を読む: 寸法 (spec.r_throat) をどちらで決めたかの記録 (plan tooling-nozzle-upstream-poly-and-throat-sizing
    §4.2)。{method: exit | throat, target_m: 目標の出口半径 / 物理スロート半径 [m], note: 任意}。prepare_ns が prepare_info.json の
    `sizing` に、実際の壁の値と目標との差を並べて書く。無ければ None (寸法の決め方は未記録)。不正な値は例外。"""
    if "sizing" not in spec:
        return None
    s = spec["sizing"]
    if not isinstance(s, dict) or s.get("method") not in ("exit", "throat"):
        raise ValueError(f"spec.sizing は {{method: exit | throat, target_m: 数値}} (受け取った値: {s!r})")
    t = s.get("target_m")
    if isinstance(t, bool) or not isinstance(t, (int, float)) or not np.isfinite(float(t)) or float(t) <= 0.0:
        raise ValueError(f"spec.sizing.target_m は正の有限の数値 [m] (受け取った値: {t!r})")
    return {"method": s["method"], "target_m": float(t), "note": s.get("note")}


def environment_record() -> dict:
    """設計チェーンを回した Python の環境 (prepare_info.json の `environment`)。積分法の δ_r (scipy の RK45) は numpy・scipy の版で
    1e-6 r_t の桁で動く (2026-10-07 実測: 同じ問題・同じコードで scipy 1.11.4 と 1.18.0 の差が出口付近で 7e-6 r_t)。壁・メッシュの
    ビット同一を要する照合は同じ環境どうしで行う (procedures/nozzle-design-workflow.md「計算環境」)。"""
    import platform
    import sys as _sys
    import scipy
    from ..geometry.moc_inverse import moc_kernel_record
    return {"python": platform.python_version(), "executable": _sys.executable, "numpy": np.__version__, "scipy": scipy.__version__,
            "host": platform.node(), "moc_kernel": moc_kernel_record()}


def require_moc_gate(p: Problem, d: dict) -> None:
    """計算準備 (`prepare`・`prepare_ns`) の入口で MOC の単位過程のゲートを必須にする (plan discretization-moc-axis-limit-and-corrector
    §4.2「反復失敗が 1 対でもあれば検証・生産は不合格」、2026-10-07 result 段レビュー M1: 以前は design_chain の診断に記録するだけで、
    不合格の MOC でもメッシュ生成へ進んでいた)。設計の経路 (design_chain) は止めない (診断を取れるように残す)。
    - `moc_corrector: converge` のとき: 診断 (`d["moc"]["gate"]`) が無い・`applicable` でない・`pass` が True でなければ ValueError。
    - `fixed2` (収束を判定しない) のときは合否を出さないので通す (ゲートの `pass` は None)。"""
    _, corrector = _moc_keys(p.geometry)
    moc = d.get("moc")
    gate = (moc or {}).get("gate")
    if corrector != "converge":
        return
    if not isinstance(gate, dict) or not gate.get("applicable"):
        raise ValueError(f"MOC のゲートの診断が無い (moc_corrector: converge; d['moc']['gate'] = {gate!r}) — 計算準備を止める")
    if gate.get("pass") is not True:
        raise ValueError(f"MOC のゲートが不合格 ({'; '.join(gate.get('reasons') or []) or gate.get('pass')}) — 計算準備を止める "
                         "(plan discretization-moc-axis-limit-and-corrector §4.2)")


def design_chain(p: Problem) -> dict:
    """Hall (+CFD アンカー) → Hermite law → 逆 MOC → 壁 QA → CFD 壁。決定的。

    ガス: `p.gas_model` (cpg なら float γ と等価 / semiperfect なら NASA-9 テーブル)。
    Hall 遷音速級数は定数 γ 前提なので**スロートの局所 γ* を渡す。MOC の ν↔M・
    質量流束・面積比はガスモデル経由 (`moc_kernel._is_gas` 規約)。"""
    if "wall_fit_mono_r2" in p.geometry and str(p.geometry.get("wall_repr", "interp")) != "joint":
        # 単調拘束 (plan tooling-nozzle-throat-monotone-r2 §4.2) は joint 壁の当てはめのオプション。他の壁表現で黙って無視しない
        raise ValueError("geometry.wall_fit_mono_r2 は wall_repr: joint 専用 "
                         f"(wall_repr = {p.geometry.get('wall_repr', 'interp')!r})")
    moc_axis_limit, moc_corrector = _moc_keys(p.geometry)     # 不正値は重い処理 (CFD ピンの読込) の前に例外
    gas = p.gas_model
    # cfd_gas: cpg のときは**設計も** CPG(γ 参照値) で作る。設計だけ semi-perfect にすると
    # 壁 (A/A*=13.1) と CFD (CPG γ=1.309 なら A/A*=15.3) の熱力学が食い違い、出口 M が
    # 3.86 で止まる (case/42 run_0003–0011 で実測、破棄)。設計と CFD の熱力学は必ず一致。
    if str(p.evaluate.get("cfd_gas", "same")) == "cpg":
        from ..gas import GasCPG
        gas = GasCPG(p.gamma, p.cp)
    use_gas = getattr(gas, "kind", "cpg") == "semiperfect"
    g = gas if use_gas else p.gamma                   # MOC/面積比に渡す「γ or ガス」
    g_hall = float(gas.gamma_throat(float(p.spec["Tt"])))
    # 報告・評価の M_d は spec.M_design のまま。MOC・軸 law・壁 QA には出口較正 Md_moc_offset を足した値を渡す
    # (plan tooling-nozzle-cfd-pinned-initial-line §5.1 #6; 既定 0 で従来と同一)
    Md_report = float(p.spec["M_design"])
    Md_moc_offset = float(p.geometry.get("Md_moc_offset", 0.0))
    Md = Md_report + Md_moc_offset if Md_moc_offset != 0.0 else Md_report
    R = float(p.geometry.get("R", 2.0))
    M_start = float(p.geometry.get("M_start", 1.05))
    n_start = int(p.geometry.get("n_start", 41))
    # 初期線の出所 (methods/design/overview.md「初期線の出所: Hall / CFD ピン」):
    # 'hall' (既定) = Hall 級数 / 'cfd' = 凍結源の node Euler 場 (CFDPinnedThroat, 線・軸アンカー・場の M)
    initial_line = str(p.geometry.get("initial_line", "hall"))
    il_src = None
    if initial_line == "hall":
        ht = HallThroat(R=R, gamma=g_hall)
    elif initial_line == "cfd":
        from ..feedback.cfd_initial_line import pinned_factory
        il_run = p.geometry.get("initial_line_run")
        if not il_run:
            raise ValueError("geometry.initial_line: cfd には geometry.initial_line_run が必須")
        il_run = Path(str(il_run))
        if not il_run.is_absolute() and p.path:
            il_run = (Path(p.path).resolve().parent / il_run).resolve()
        il_res = p.geometry.get("initial_line_res")
        if not il_res:
            raise ValueError("geometry.initial_line: cfd には geometry.initial_line_res (凍結源の snapshot) が必須")
        # 凍結源の形・ガスを使う側と照合 (codex result M2): 縮流部 (r_U・L_U・L_pipe・R の U→T Hermite)・組成・Tt・γ_Hall
        expect = {"gas": (gas.summary() if hasattr(gas, "summary") else {"kind": "cpg"}),
                  "r_U": float(p.geometry.get("r_inlet", 2.5)), "L_U": float(p.geometry.get("L_U", 3.5)),
                  "L_pipe": float(p.geometry.get("L_pipe", 0.5))}
        # 実効入口 (BC の Pt・Tt・組成) と熱力学条件 (thermalMethod・species・thermoHrefTemp) も照合 (codex diagnose 2026-10-06)
        from ..feedback.cfd_initial_line import expected_inlet_thermo
        expect["inlet_thermo"] = expected_inlet_thermo(p)
        ht = pinned_factory(il_run, il_res, expect=expect)(R, g_hall)
        il_src = dict(ht.source)
    else:
        raise ValueError("geometry.initial_line は 'hall' か 'cfd'")
    # 初期値線: 'throat_char' = スロート壁点発 C⁻ (CONTUR 流、壁が T から MOC 出力に
    # なる) / 'vertical' = M_start の縦線 (旧構成、回帰対照)。x0 = その軸着地点。
    start_line = str(p.geometry.get("start_line", "vertical"))
    if initial_line == "cfd" and start_line != "throat_char":
        raise ValueError("geometry.initial_line: cfd は start_line: throat_char 専用 (CFD 線の軸アンカーは線の軸着地でだけ定義)")
    if start_line == "throat_char":
        x0 = float(ht.throat_characteristic(n=n_start)[0][0])
    elif start_line == "vertical":
        x0 = ht.x_axis_of_mach(M_start)
    else:
        raise ValueError("geometry.start_line は 'throat_char' か 'vertical'")

    # --- アンカー: 初回 = Hall 解析値 (x_A = x0)。反復 = CFD 実測 (x_A = x_reach) ---
    x_reach_cfd = p.geometry.get("x_reach_cfd")
    seg = None
    if x_reach_cfd is not None:
        x_A = float(x_reach_cfd)
        anc = p.geometry.get("x_reach_anchor")
        if anc is None or len(anc) != 3:
            raise ValueError("x_reach_cfd 指定時は x_reach_anchor: [M,M',M''] が必須")
        M_A, Mp_A, Mpp_A = float(anc[0]), float(anc[1]), float(anc[2])
        seg_run = p.geometry.get("axis_segment_run")
        if seg_run:
            seg, _, _ = axis_curve_node(seg_run, float(p.spec["r_throat"]))
    else:
        x_A = x0
        M_A, Mp_A, Mpp_A = ht.axis_anchor(x0)
    # アンカーの出所を成分別に (codex result m6): CFD ピンは M・M′ が CFD 場、M″ は Hall の式を CFD の x0 で評価
    if x_reach_cfd is not None:
        anchor_src = {"M": "cfd_reach", "Mp": "cfd_reach", "Mpp": "cfd_reach"}
    elif initial_line == "cfd":
        anchor_src = {"M": "cfd", "Mp": "cfd", "Mpp": "hall@x0_cfd"}
    else:
        anchor_src = {"M": "hall", "Mp": "hall", "Mpp": "hall"}

    # --- 軸 M 則: 'quintic' = 単一 5 次 Hermite (DOF L_c) / 'knot' = 内部 knot 1 個の
    # 区分 C² (A6, DOF L_c + M_knot)。高マッハ (M6) では単一 quintic の L_c 上限
    # (単調性) が短すぎて壁角 θ_w > μ_w の fold を起こすので knot を使う。
    # 'bspline_M' (B) / 'bspline_dnu' (C): axislaw-smoothness 比較 (plans/accepted/
    # tooling-nozzle-axislaw-smoothness.md)。$L_c$ に上限はなく常に Lc_mode: explicit。
    axis_law = str(p.geometry.get("axis_law", "quintic"))
    M_K = None
    if axis_law == "quintic":
        lo, hi = QuinticHermiteAxisLaw.admissible_Lc_range(x_A, M_A, Mp_A, Mpp_A, Md)
    elif axis_law == "knot":
        # 既定 knot Mach = 急膨張の終わり。M_A + 0.25ΔM (M4.2 で 1.9、M6 で 2.4)
        M_K = float(p.geometry.get("M_knot", M_A + 0.25 * (Md - M_A)))
        lo, hi = KnotQuinticAxisLaw.admissible_Lc_range(x_A, M_A, Mp_A, Mpp_A, Md, M_K)
    elif axis_law in ("bspline_M", "bspline_dnu", "onepoint"):
        lo, hi = 1e-2, float("inf")
    else:
        raise ValueError("geometry.axis_law は 'quintic'/'knot'/'bspline_M'/'bspline_dnu'/'onepoint'")
    # --- L_c の決め方 (Lc_mode): 許容窓は常に計算する ---
    #   'explicit'    = dv.L_c を直接与える (窓内検査)
    #   'max'         = 窓上限×0.98
    #   'from_length' = dv.L_total (スロート x=0 → 物理出口 x_F の長さ, r_t 単位) を
    #                   設計変数とし、軸 M 極大点 x_E (= x_A + L_c) は終端特性線込みの
    #                   設計パスから逆算する (E→F 一様化区間は物理で決まるため)。
    #                   plans/accepted/tooling-nozzle-axismach-length-dv.md

    def _make_law(L_c: float):
        if axis_law == "knot":
            return KnotQuinticAxisLaw(x_A, L_c, M_A, Mp_A, Mpp_A, Md, M_K)
        if axis_law == "bspline_M":
            from ..geometry.axis_law_bspline import MonotoneBSplineAxisLaw
            return MonotoneBSplineAxisLaw(
                x_A, x_A + L_c, M_A, Mp_A, Mpp_A, Md,
                n_interior=int(p.geometry.get("bspline_n_interior", 15)),
                exit_curvature=str(p.geometry.get("bspline_exit_curvature", "hard")),
                spread_x_frac=float(p.geometry.get("bspline_spread_x_frac", 0.5)),
                spread_M_frac=p.geometry.get("bspline_spread_M_frac", 0.75))
        if axis_law == "bspline_dnu":
            from ..geometry.axis_law_bspline import NonnegDnuBSplineAxisLaw
            return NonnegDnuBSplineAxisLaw(
                x_A, x_A + L_c, M_A, Mp_A, Mpp_A, Md, gas,
                n_interior=int(p.geometry.get("bspline_n_interior", 20)),
                exit_slope=str(p.geometry.get("bspline_exit_curvature", "hard")),
                spread_x_frac=float(p.geometry.get("bspline_spread_x_frac", 0.5)),
                spread_M_frac=p.geometry.get("bspline_spread_M_frac", 0.75))
        if axis_law == "onepoint":
            # D 案 (plans/accepted/tooling-nozzle-axislaw-onepoint.md): 端点アンカー固定 +
            # 内部補間点 1 点 (ξ_P, η_P) の C⁴ 区分 5 次。実験用選択肢。
            from ..geometry.axis_law_onepoint import OnePointC4AxisLaw
            return OnePointC4AxisLaw(x_A, L_c, M_A, Mp_A, Mpp_A, Md,
                                     float(p.geometry["onepoint_xi_P"]),
                                     float(p.geometry["onepoint_eta_P"]))
        return QuinticHermiteAxisLaw(x_A, L_c, M_A, Mp_A, Mpp_A, Md)

    def _checked_law(L_c: float):
        law = _make_law(L_c)
        v = law.gates()["violations"]
        if v:
            raise ValueError("軸 Mach law ゲート不合格: " + "; ".join(v))
        return law

    rF_pred = float(np.sqrt(area_ratio_isentropic(Md, g)))
    # 逆 MOC の単位過程 (plans/accepted/discretization-moc-axis-limit-and-corrector.md §4.3):
    # geometry.moc_axis_limit = analytic (既定、軸則から解析極限 θ_r) | legacy (軸端点の sinθ/r は相手の値で代用)
    # geometry.moc_corrector = fixed2 (既定、予測 1 + 修正 2 回) | converge (更新量 ≤ 1e-12 まで、上限 50 回)
    # キーが無ければ従来とビット同一。キーの検査は design_chain の冒頭 (`_moc_keys`)

    def _run_inverse(law):
        # target: [x0, x_A) は実測 (反復時) / [x_A, x_E] law / 以降 M_d
        def target_moc(x: float) -> float:
            x = float(x)
            if x < x_A:
                if seg is not None:
                    return float(seg(np.float64(x)))
                return float(ht.mach(x, 0.0)) if x_reach_cfd is None else float(law(x_A))
            return float(law(x))

        def target_dM(x: float) -> float:
            """target_moc の dM/dx (analytic の θ_r 用。軸則の解析微分)。区間の分け方は target_moc と同じ。"""
            x = float(x)
            if x < x_A:
                if seg is not None:
                    return float(seg(np.float64(x), 1))
                if x_reach_cfd is None:
                    # x_A = x0 なので軸節点 (x > x0) では通らない。Hall の級数なら解析微分、他 (CFD ピン) は未定義
                    return float(ht.axis_anchor(x)[1]) if type(ht) is HallThroat else float("nan")
                return 0.0                                  # target = law(x_A) の定数
            return float(law.deriv(x, 1))
        x_end = law.x_E + float(p.geometry.get("x_end_margin", 2.3)) \
            * rF_pred * float(np.sqrt(Md * Md - 1.0))
        kw_moc = {}
        if (moc_axis_limit, moc_corrector) != ("legacy", "fixed2"):
            kw_moc = dict(axis_limit=moc_axis_limit, corrector=moc_corrector, target_dM=target_dM,
                          axis_anchor=(x_A, M_A, Mp_A))
        return inverse_design(ht, target_moc, x_axis_end=float(x_end),
                              n_axis=int(p.geometry.get("n_axis_inv", 500)),
                              n_start=n_start, gamma=g,
                              dx_wall=float(p.geometry.get("dx_wall", 0.02)),
                              th_wall0=float(np.arctan(x0 / R)), M_start=M_start,
                              exit_mode=str(p.geometry.get("exit_mode", "characteristic")),
                              x_E=law.x_E, M_d=Md, start_line=start_line,
                              wall_mode=str(p.geometry.get("wall_mode", "streamline")),
                              blend_width=float(p.geometry.get("wall_blend_width", 1.0)),
                              axis_dx0=p.geometry.get("axis_dx0"), **kw_moc)

    mode = str(p.geometry.get("Lc_mode", "explicit"))
    solve_diag = None
    if mode == "max":
        if axis_law in ("knot", "bspline_M", "bspline_dnu", "onepoint") and hi >= 400.0:
            raise ValueError(f"{axis_law} 則の L_c に上限はない (窓が開放) — "
                             "Lc_mode: explicit で L_c を与えること")
        L_c = hi * 0.98
        law = _checked_law(L_c)
        res = _run_inverse(law)
    elif mode == "explicit":
        L_c = float(dv_value(p, "L_c"))
        if not (lo <= L_c <= hi):
            raise ValueError(f"L_c = {L_c:.4g} が許容窓 ({lo:.4g}, {hi:.4g}) 外 "
                             "(M'≥0 単調ゲート)")
        law = _checked_law(L_c)
        res = _run_inverse(law)
    elif mode == "from_length":
        # 逆問題: x_F(L_c) = L_total を L_c について解く。x_F − x_E ≈ r_F√(Md²−1)
        # (終端 Mach line 近似, 実測 0.04% 一致) なので写像はほぼ**傾き 1** の単調。
        # ただし離散 MOC の x_F には解像度依存のノイズ床がある (n_axis 500 で
        # ~0.05 r_t — 隣接 L_c 間で階段状。secant の局所勾配推定は壊れる) ため、
        # ステップ = −残差 の固定点反復 (縮小率 |1 − dx_F/dL_c| ≈ 0.05) を使い、
        # 最良反復点を採用する。既定 tol 0.05 はこのノイズ床相当 — より詰めるなら
        # n_axis_inv を上げて L_total_tol を下げる。
        L_target = float(dv_value(p, "L_total"))
        tol = float(p.geometry.get("L_total_tol", 0.05))
        n_iter_max = int(p.geometry.get("L_total_maxiter", 8))
        # 窓端は数値マージンを取ってクランプ (単調ゲート境界での law 破綻を避ける)
        lo_c = lo * 1.001
        hi_c = hi * 0.999 if np.isfinite(hi) else hi

        def _clamp(v: float) -> float:
            return float(min(max(v, lo_c), hi_c))

        n_eval = 0

        def _eval(L_c: float):
            nonlocal n_eval
            law = _checked_law(L_c)
            res = _run_inverse(law)
            xF = float(res["exit"].get("x_F", float("nan")))
            if not np.isfinite(xF):
                raise ValueError("Lc_mode: from_length — 終端特性線の x_F 追跡に失敗 "
                                 "(場が短い: geometry.x_end_margin を増やす)")
            n_eval += 1
            return xF - L_target, law, res

        L_c = _clamp(L_target - x_A - rF_pred * float(np.sqrt(Md * Md - 1.0)))
        f0, law, res = _eval(L_c)
        best = (abs(f0), L_c, f0, law, res)
        for _ in range(n_iter_max):
            if abs(f0) <= tol:
                break
            L_new = _clamp(L_c - f0)
            if L_new == L_c:               # 窓端に張り付いた — 窓内で到達不能
                raise ValueError(
                    f"dv.L_total = {L_target:.4g} は L_c 許容窓 ({lo:.4g}, {hi:.4g}) "
                    f"内で実現不能 (窓端 L_c = {L_c:.4g} で x_F = {L_target + f0:.4g})")
            f0, law, res = _eval(L_new)
            L_c = L_new
            if abs(f0) < best[0]:
                best = (abs(f0), L_c, f0, law, res)
            else:                          # 改善停止 = ノイズ床に到達
                break
        _, L_c, f0, law, res = best
        if abs(f0) > tol:
            raise ValueError(
                f"Lc_mode: from_length 未収束 (|x_F − L_total| = {abs(f0):.3g} > "
                f"tol {tol:.3g}) — 逆 MOC の離散化ノイズ床の可能性: "
                "geometry.n_axis_inv を上げるか L_total_tol を緩める")
        solve_diag = {"L_total_target": float(L_target),
                      "xF_residual": float(f0), "n_design_evals": int(n_eval),
                      "tol": float(tol)}
    else:
        raise ValueError("geometry.Lc_mode は 'explicit' / 'max' / 'from_length'")
    gates = law.gates()
    qa = wall_qa(res["wall"], Md, law.x_E, g, R=R)
    if qa["violations"]:
        raise ValueError("壁 QA 不合格: " + "; ".join(qa["violations"]))
    # 壁表現 (A14): 'interp' = 補間 5 次 B-spline (現行・比較基準) /
    # 'lsq' = 制約付き最小二乗 B-spline (n_cp は誤差ゲートで自動、または wall_ncp) /
    # 'joint' = 位置 + 壁角の同時当てはめ (V0 型壁, plan tooling-nozzle-cfd-pinned-initial-line §4.5)
    wall_repr = str(p.geometry.get("wall_repr", "interp"))
    wkw = dict(R=R, r_U=float(p.geometry.get("r_inlet", 2.5)),
               L_U=float(p.geometry.get("L_U", 3.5)),
               L_pipe=float(p.geometry.get("L_pipe", 0.5)))
    if wall_repr == "lsq":
        from ..geometry.wall_axismach import LSQBsplineCFDWall
        wall = LSQBsplineCFDWall(res["wall"], n_cp=p.geometry.get("wall_ncp"),
                                 tol_dr=float(p.geometry.get("wall_lsq_tol_dr", 5e-4)),
                                 tol_dtheta_deg=float(p.geometry.get("wall_lsq_tol_dtheta", 0.05)),
                                 **wkw)
    elif wall_repr == "interp":
        wall = AxisMachCFDWall(res["wall"], **wkw)
    elif wall_repr == "joint":
        from ..geometry.wall_axismach import JointFitCFDWall
        # geometry.wall_fit_mono_r2 = [a, b]: [a, b] r_t で r″ を単調非増加に拘束 (plan tooling-nozzle-throat-monotone-r2 §4.2)。
        # 無い (または null) なら拘束なし = 従来とビット同一
        wall = JointFitCFDWall(res["wall"], mono_r2=p.geometry.get("wall_fit_mono_r2"), **wkw)
    else:
        raise ValueError("geometry.wall_repr は 'interp' / 'lsq' / 'joint'")
    msgs = wall.validate()
    if msgs:
        raise ValueError("axis-Mach 壁フィルタ不合格: " + "; ".join(msgs))
    return {"wall": wall, "wall_inv": res["wall"], "law": law, "qa": qa,
            "exit": {k: v for k, v in res["exit"].items() if k != "term_path"},
            "x0": float(x0), "x_A": float(x_A), "x_E": float(law.x_E),
            "L_c": float(L_c), "Lc_window": (float(lo), float(hi)),
            "Lc_mode": mode, "Lc_solve": solve_diag,
            "axis_law": axis_law, "M_knot": M_K,
            "x_K": (float(law.x_K) if axis_law == "knot" else None),
            "anchor": (float(M_A), float(Mp_A), float(Mpp_A)),
            "anchor_source": anchor_src,
            "start_line": start_line, "wall_mode": res["wall_mode"],
            "wall_repr": wall_repr,
            # 逆 MOC の単位過程の診断 (plan discretization-moc-axis-limit-and-corrector §4.2): キーの値・対の 5 分類・反復回数・
            # 最終残差・源項の分岐 (AXIS_LIMIT_FRAC の発火の数と位置)・θ_r の出所と軸端の接続検査・ゲート (converge のとき合否)
            "moc": res.get("moc"),
            "gas": (gas.summary() if hasattr(gas, "summary")
                    else {"kind": "cpg", "gamma": p.gamma, "cp": p.cp}),
            "gamma_hall": g_hall,
            "wall_fit": getattr(wall, "fit_diag", None),
            "Md": Md_report, "Md_moc_offset": Md_moc_offset, "Md_moc": Md, "R": R, "gates": gates,
            "initial_line": {"source": initial_line,
                             "run": (il_src or {}).get("run"), "res": (il_src or {}).get("res"),
                             "x0": float(x0), "anchor": [float(M_A), float(Mp_A), float(Mpp_A)],
                             "mstar": float(res["mdot_start"]),
                             "sha256_16": (il_src or {}).get("sha256_16"),
                             "config": (il_src or {}).get("config"),
                             "match": (il_src or {}).get("match")},
            # cplus 閉包では構成的に 1 になる循環指標なので出さない (A9 の教訓)
            "mdot_ratio_moc": (None if not np.isfinite(res["mdot_exit"])
                               else float(res["mdot_exit"] / res["mdot_start"])),
            "cd_series": float(ht.cd_series())}


def _mesh_params_from(m: dict, scale, ni, nj, wall_first_frac) -> Mesh2DParams:
    """格子のブロック (dict) → Mesh2DParams。既定値は呼び出し側の ni/nj/wall_first_frac だけが違う
    (2026-10-06 codex diagnose: 以前の Euler 経路は ni/nj/wall_first_frac/throat_refine しか渡さず、throat_width・
    wall_first_frac_throat・前後ブレンド等を黙って無視していた)。"""
    opt = lambda k: None if m.get(k) is None else float(m[k])  # noqa: E731
    return Mesh2DParams(ni=int(m.get("ni", ni)), nj=int(m.get("nj", nj)),
                        wall_first_frac=float(m.get("wall_first_frac", wall_first_frac)),
                        throat_refine=float(m.get("throat_refine", 3.0)),
                        throat_width=float(m.get("throat_width", 1.5)),
                        wall_first_frac_throat=opt("wall_first_frac_throat"),
                        wall_first_blend_x0=float(m.get("wall_first_blend_x0", 0.5)),
                        wall_first_blend_x1=float(m.get("wall_first_blend_x1", 6.0)),
                        wall_first_up_x0=opt("wall_first_up_x0"), wall_first_up_x1=opt("wall_first_up_x1"),
                        axis_gap_frac=opt("axis_gap_frac"), axis_cap_frac=opt("axis_cap_frac"),
                        wall_first_frac_table=m.get("wall_first_frac_table"), x_density_table=m.get("x_density_table"),
                        wall_normal_layer=m.get("wall_normal_layer"),
                        # 区間の局所の並べ直し (NS の mesh だけ。MESH_EULER_KEYS に無いので mesh_euler には書けない;
                        # plan tooling-nozzle-isothermal-wall-chain §4.8-1・§5.1 #29)
                        x_respace=m.get("x_respace"), scale=scale)


def mesh_params(p, scale, ni, nj, wall_first_frac):
    """problem の mesh ブロック → Mesh2DParams。**NS (`prepare_ns`) の格子**。
    2026-10-07 (plan verification-case45-euler-total-enthalpy §4「E3 以降の対処」) から Euler (`prepare`) は `mesh_euler` を読む
    (`mesh_params_euler`)。それ以前は Euler と NS が同じ mesh を読んでいた。"""
    return _mesh_params_from(p.mesh, scale, ni, nj, wall_first_frac)


# Euler 専用の格子の設定 `mesh_euler` (2026-10-07, plan verification-case45-euler-total-enthalpy §4「E3 以降の対処」、
# 諮問 notes/reviews/2026-10-07-euler-grid-switch-plan-diagnose.md):
# - prepare (Euler) は mesh_euler、prepare_ns (NS) は mesh を読む。2 つのブロックを混ぜて補完しない (Euler への mesh の暗黙の継承を禁止)。
# - 既定は全断面で wall_first_frac 0.005 の等比の配点、スロートの別指定なし、軸側の cap なし (ni・nj の既定は従来の Euler の 321 × 65)。
#   NS 向けに壁へ寄せた配点 (case/45 の G1: 1.3e-5・スロート 4.5e-6) の Euler は、スロート付近の全温が Tt を数百 K 超えたまま整定しなかった
#   (同 plan §9 E2。全域 0.005 の配点は同じ窓で全領域 |T₀ − Tt| ≤ 0.103 K)。
# - mesh だけがあって mesh_euler が無い問題は移行先を示して止める (黙って NS の配点を使わない・指定を無視して既定に落とさない)。
MESH_EULER_DEFAULTS = {"ni": 321, "nj": 65, "wall_first_frac": 5.0e-3}
# mesh_euler に書けるキー (_mesh_params_from が読むもの + 品質検査の ar_max + node 固定の discretization + 節点数の上限の許可)。
# それ以外は書き誤りとして止める
MESH_EULER_KEYS = ("ni", "nj", "wall_first_frac", "throat_refine", "throat_width", "wall_first_frac_throat",
                   "wall_first_blend_x0", "wall_first_blend_x1", "wall_first_up_x0", "wall_first_up_x1",
                   "axis_gap_frac", "axis_cap_frac", "ar_max", "discretization", "max_nodes", "max_nodes_reason")

# 2D 軸対称の格子の節点数の上限 (2026-10-14, plan tooling-nozzle-mesh-size-guard §4.2): Euler・NS とも 30 万点 (Euler の ni は明示が無ければ
# NS と同じなので同じ値; 雛形 1600〜2000 × 97 = 15.5〜19.4 万、M6 生産の NS 19.4 万の 1.5 倍)。超える格子は問題 YAML のブロック
# (mesh / mesh_euler) に max_nodes (正の整数) と max_nodes_reason (空でない文字列) を組で書けば、その値まで許す (ユーザの了承を取ってから書く)。
# 判定は run dir を作る前・設計チェーンの前に、実効の節点数 (mesh2d.node_count; axis_gap_frac の導出 nj を含む) で行う。
MAX_NODES_2D = 300_000
NS_STEPS_DEFAULT = 60_000      # NS の本段の既定の step 数 (2026-10-14 ユーザ決定、手順書 §2 ⑥ の最終 NS の標準。evaluate.ns_steps で上書き)


def node_budget(m: dict, block: str) -> dict:
    """格子のブロックの節点数の上限 {max_nodes, max_nodes_reason, max_nodes_source ("default" / "explicit")}。
    max_nodes と max_nodes_reason は組で書く (片方だけ・正の整数でない値・空の理由は ValueError)。"""
    m = m or {}
    v, why = m.get("max_nodes"), m.get("max_nodes_reason")
    if v is None and why is None:
        return {"max_nodes": MAX_NODES_2D, "max_nodes_reason": None, "max_nodes_source": "default"}
    errs = []
    if v is None:
        errs.append(f"{block}.max_nodes_reason だけがあり {block}.max_nodes (許す節点数) が無い")
    elif isinstance(v, bool) or not isinstance(v, (int, np.integer)) or int(v) <= 0:
        errs.append(f"{block}.max_nodes は正の整数 (受け取った値: {v!r})")
    if why is None:
        errs.append(f"{block}.max_nodes だけがあり理由 {block}.max_nodes_reason が無い")
    elif not isinstance(why, str) or not why.strip():
        errs.append(f"{block}.max_nodes_reason は空でない文字列 (受け取った値: {why!r})")
    if errs:
        raise ValueError("節点数の上限の許可は max_nodes と max_nodes_reason を組で書く: " + "; ".join(errs))
    return {"max_nodes": int(v), "max_nodes_reason": why, "max_nodes_source": "explicit"}


def check_node_budget(mp: Mesh2DParams, m: dict, block: str, allow_hint: str | None = None) -> dict:
    """格子を作る前の節点数の判定 (plan tooling-nozzle-mesh-size-guard §4.2)。実効の節点数 (mesh2d.node_count、上限を超えたら数えるのを
    打ち切る) が上限 (`node_budget`) を超えたら ValueError (実効の節点数・ni・nj・上限・許可の書き方)。allow_hint は許可の書き方の文言
    (既定は問題 YAML のブロックに書く形; 文言の {n_nodes} は数えた節点数に置き換える)。
    戻り値 {n_nodes, ni, nj, max_nodes, max_nodes_reason, max_nodes_source}。"""
    b = node_budget(m, block)
    lim = b["max_nodes"]
    ni = int(mp.ni)
    n = node_count(mp, limit=lim)
    nj, exact = radial_node_count(mp, max_nj=lim // max(ni, 1))
    if n > lim:
        how = (allow_hint or f"問題 YAML の {block} に max_nodes: {{n_nodes}} と max_nodes_reason: \"<理由>\" を組で書く"
               ).replace("{n_nodes}", str(n) if exact else f"<{n} 以上>")
        src = "既定 MAX_NODES_2D" if b["max_nodes_source"] == "default" else f"明示: {b['max_nodes_reason']}"
        more = "" if exact else " 以上"
        derived = "、nj は axis_gap_frac からの導出値" if mp.axis_gap_frac is not None else ""
        raise ValueError(
            f"{block} の格子の実効の節点数 {n}{more} (ni {ni} × nj {int(nj)}{more}{derived}) が上限 {lim} ({src}) を超える — "
            f"格子を作らずに止めた。上限を超える格子はユーザの了承を取ってから {how} (了承なしに書かない)。"
            "超えないようにするなら ni・nj・axis_gap_frac (・要求カードなら第 1 層と ar_max) を見直す")
    return {"n_nodes": int(n), "ni": ni, "nj": int(nj), **b}


def _check_generated_nodes(coords, nodes: dict, block: str) -> None:
    """格子を作った後の照合: 実際の節点数が事前の数え (mesh2d.node_count) と違えば止める (事前の判定が当てにならないため)。"""
    if int(coords.shape[0]) != int(nodes["n_nodes"]):
        raise RuntimeError(f"{block}: 生成した格子の節点数 {coords.shape[0]} が事前の数え (mesh2d.node_count) {nodes['n_nodes']} と"
                           " 一致しない — 事前の節点数の判定が当てにならないので止める")


def mesh_euler_block(p) -> dict:
    """Euler の格子のブロック `mesh_euler` (dict の写し)。無い・dict でない・未知のキー・node でない discretization は ValueError。
    mesh だけの問題には、旧格子を再現する写し方と新しい既定の書き方を示して止める。"""
    m = (p.raw or {}).get("mesh_euler")
    if m is None:
        old = p.raw.get("mesh") or {}
        flow = lambda d: "{" + ", ".join(f"{k}: {v}" for k, v in d.items()) + "}"  # noqa: E731  (YAML の flow 形式で示す)
        axial = {k: old[k] for k in ("ni", "nj", "throat_refine", "throat_width") if k in old}
        raise ValueError(
            "prepare (Euler) は問題 YAML の mesh_euler を読む (2026-10-07 から。plan verification-case45-euler-total-enthalpy §4)。"
            "この問題には mesh だけがあり mesh_euler が無い — NS の mesh を Euler に流用しない。移行先: "
            "(1) 新しい既定 (全断面 wall_first_frac 0.005・スロートの別指定なし・軸側の cap なし) なら "
            f"mesh_euler: {flow({**axial, 'wall_first_frac': 0.005})}、"
            f"(2) 旧格子の記録を再現するなら mesh の全キーを mesh_euler に写す (mesh_euler: {flow(old)})")
    if not isinstance(m, dict):
        raise ValueError(f"mesh_euler は辞書 (受け取った値: {m!r})")
    bad = sorted(str(k) for k in m if k not in MESH_EULER_KEYS)
    if bad:
        raise ValueError(f"mesh_euler の未知のキー {bad} (書けるのは {list(MESH_EULER_KEYS)}; 黙って無視しない)")
    if "discretization" in m and str(m["discretization"]) != "node":
        raise ValueError(f"mesh_euler.discretization は node だけ (axis-Mach の Euler は node; 受け取った値: {m['discretization']!r})")
    return dict(m)


def mesh_params_euler(p, scale) -> Mesh2DParams:
    """problem の mesh_euler ブロック → Mesh2DParams。**Euler (`prepare`) の格子** (既定は MESH_EULER_DEFAULTS と Mesh2DParams の
    スロートの別指定なし・軸側の cap なし)。mesh は読まない。"""
    return _mesh_params_from(mesh_euler_block(p), scale, **MESH_EULER_DEFAULTS)


# 格子の座標と接続のハッシュ (prepare_info.json の mesh 欄)。式は case/45.isobutane_m6_d155/euler_t0_e2.py の E2_MESH.json
# (`_sha_bytes(/MESH/COORD)`・`_topology_digest`) と同じ (データセット名・dtype・形・バイト列) なので、E2 の記録と直接比べられる
_MESH_TOPOLOGY = ("MESH/CONNE", "VIZMESH/CONNE", "CELLS/STRUCT", "PLANES/STRUCT", "CELLS/regionId")
_MESH_BCOND_TOPO = ("iBPlanes", "iCells", "iPlanes", "vizBfaceNodes", "vizBfaceSizes")


def h5_mesh_hashes(h5_path) -> dict:
    """格子の h5 の座標 (/MESH/COORD) と接続 (_MESH_TOPOLOGY + 境界の位相) の sha256 (準備の記録 `prepare_info.json` の
    `mesh.hashes.coord_sha256`・`topology_sha256` と同じ式; `wt_qualify` が計算の格子と準備の記録の照合に使う)。"""
    import hashlib

    import h5py

    with h5py.File(h5_path, "r") as f:
        c = f["/MESH/COORD"][:]
        hc = hashlib.sha256()
        hc.update(str(c.dtype).encode() + str(c.shape).encode() + np.ascontiguousarray(c).tobytes())
        ht = hashlib.sha256()
        names = [t for t in _MESH_TOPOLOGY if t in f]
        names += sorted(f"BCONDS/{b}/{k}" for b in f["BCONDS"] for k in _MESH_BCOND_TOPO if k in f["BCONDS"][b])
        for n in names:
            a = np.asarray(f[n])
            ht.update(n.encode() + str(a.dtype).encode() + str(a.shape).encode() + a.tobytes())
    return {"coord_sha256": hc.hexdigest(), "coord_dtype": str(c.dtype), "topology_sha256": ht.hexdigest(), "topology_datasets": names}


def _mesh_hashes(run_dir: Path, coords, quads, bedges) -> dict:
    """変換後の nozzle.h5 の座標 (/MESH/COORD) と接続 (_MESH_TOPOLOGY + 境界の位相) の sha256、nozzle.msh の sha256、
    生成器の出力 (coords float64・quads・境界辺; 変換器に依らない) の sha256。"""
    import hashlib

    def _upd(h, name, a):
        a = np.ascontiguousarray(a)
        h.update(name.encode() + str(a.dtype).encode() + str(a.shape).encode() + a.tobytes())

    hh = h5_mesh_hashes(run_dir / "nozzle.h5")
    hg = hashlib.sha256()
    _upd(hg, "coords", np.asarray(coords, dtype=np.float64))
    _upd(hg, "quads", np.asarray(quads, dtype=np.int64))
    for g in ("inlet", "outlet", "wall", "axis"):
        _upd(hg, g, np.asarray(bedges[g], dtype=np.int64))
    hm = hashlib.sha256((run_dir / "nozzle.msh").read_bytes())
    return {"coord_sha256": hh["coord_sha256"], "coord_dtype": hh["coord_dtype"], "topology_sha256": hh["topology_sha256"],
            "topology_datasets": hh["topology_datasets"], "msh_sha256": hm.hexdigest(), "generator_sha256": hg.hexdigest()}


def _mesh_record(p, mp: Mesh2DParams, source: str, nj: int, run_dir: Path, coords, quads, bedges, nodes: dict) -> dict:
    """prepare_info.json の mesh 欄: 従来の ni・nj・wall_first_frac (読む側: nozzle_report・metrics.deltastar・cfd_initial_line) に、
    全 Mesh2DParams・採用元のブロック名・問題に書かれたブロックの有無と中身・座標と接続のハッシュを足す
    (plan verification-case45-euler-total-enthalpy §4、2026-10-07)。nodes (`check_node_budget` の戻り値) から節点数と上限・理由・
    上限の出所 (default / explicit) も書く (plan tooling-nozzle-mesh-size-guard §4.2、2026-10-14)。"""
    import dataclasses
    return {"ni": mp.ni, "nj": int(nj), "wall_first_frac": mp.wall_first_frac,
            "source": source, "blocks_present": {b: (b in (p.raw or {})) for b in ("mesh", "mesh_euler")},
            "block": dict((p.raw or {}).get(source) or {}), "params": dataclasses.asdict(mp),
            "n_nodes": int(coords.shape[0]),
            **{k: nodes[k] for k in ("max_nodes", "max_nodes_reason", "max_nodes_source")},
            "hashes": _mesh_hashes(run_dir, coords, quads, bedges)}


def prepare(problem_path, run_dir, nsteps=None, ic_from=None, cfl_main=None, implicit_relax=None) -> dict:
    """**Euler (slip 壁) の run** を準備する。格子は `mesh_euler`、壁は `euler_reference_wall`:

    - `geometry.physical_upstream` が無い (今と同じ): 設計壁 (geometry の r_inlet・L_U・L_pipe の直管と縮流部 + MOC の設計壁)。
      `physical_wall_repr: single_bspline`・`pw_upstream: poly` の明示は物理壁の指定なので例外 (黙って無視しない)。
    - ある (plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #36b、#36e で案 A に): **H の上流を付けた設計壁**
      (`UpstreamDesignWall`: x ≥ 0 は凍結した設計壁 S、x < 0 は physical_upstream の直管と S(0) への 5 次多項式。下流は設計壁のゲート
      `monotone_gate`、上流は単調な収縮・正の半径・一意な大域最小・C² の継ぎ目で判定)。δ* ループの固定の Euler の参照 (抽出器の δ_in =
      NS の壁 − この壁) に使う。`pw_upstream` は poly だけ (_physical_upstream が ramp を拒む)。物理壁ではないので `physical_wall_repr` の
      明示は例外。壁ファイル (`wall_repr.json`、種類 `design_upstream`) を書いて読み直し、復元した壁が元と一致することを確かめる。
      prepare_info.json の `physical_upstream` と `euler_wall` に記録する。設計壁と CFD ピンの照合は geometry の値のまま (design_chain)。"""
    p = load_problem(problem_path)
    if p.type != "wind_tunnel_axisym_axismach":
        raise ValueError("runner_axismach は wind_tunnel_axisym_axismach 専用")
    # 物理壁だけの上流の形 (§5.1 #34) の検査 (不正値・joint でない壁・ramp は例外)。ある問題は H の上流を付けた設計壁で格子を作る
    # (§5.1 #36b、#36e の案 A)
    pu = _physical_upstream(p.geometry)
    pwu = None
    if pu is None:
        if _physical_wall_repr(p.geometry) == "single_bspline":
            # 物理壁の表現 (plan tooling-nozzle-wall-single-bspline §4.2) は物理壁だけ。Euler の設計壁で黙って無視しない
            raise ValueError("geometry.physical_wall_repr: single_bspline は物理壁の指定 — Euler の prepare は physical_upstream の無い問題では"
                             "設計壁で格子を作る")
        if "pw_upstream" in p.geometry:
            # 上流の作り方 (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1) は joint 壁の物理壁の解析経路だけ。
            # 値は検査し、poly の明示は例外 (Euler の設計壁で黙って無視しない)。ramp (移行で明示した旧来の作り方) は通す
            v = p.geometry["pw_upstream"]
            if not isinstance(v, str) or v not in PW_UPSTREAMS:
                raise ValueError(f"geometry.pw_upstream は {' / '.join(repr(c) for c in PW_UPSTREAMS)} のどれか (受け取った値: {v!r})")
            if v == "poly":
                raise ValueError("geometry.pw_upstream: poly は物理壁の指定 — Euler の prepare は physical_upstream の無い問題では設計壁で格子を作る")
    else:
        pwu = _pw_upstream(p.geometry)            # poly (既定・明示) だけ (_physical_upstream が ramp を拒む)
        if _physical_wall_repr(p.geometry) is not None:
            # 案 A (§5.1 #36e): Euler の参照の壁は H の上流を付けた設計壁 (区分表現) で物理壁ではない。物理壁の表現の指定を黙って無視しない
            raise ValueError(f"geometry.physical_wall_repr: {_physical_wall_repr(p.geometry)} は物理壁の指定 — Euler の prepare は "
                             "physical_upstream の問題を H の上流を付けた設計壁 (壁ファイル design_upstream) で準備する")
    # Euler の格子は mesh_euler (mesh は読まない)。無い・不正なら run dir を作る前・設計チェーンの前に止める (2026-10-07)
    m_eu = mesh_euler_block(p)
    # 節点数の上限 (plan tooling-nozzle-mesh-size-guard §4.2): 実効の節点数で run dir を作る前・設計チェーンの前に判定する
    # (節点数は scale によらないので scale 1 の Mesh2DParams で数える)
    nodes = check_node_budget(mesh_params_euler(p, 1.0), m_eu, "mesh_euler")
    run_dir = Path(run_dir)
    if run_dir.exists():
        raise FileExistsError(f"{run_dir} が既にある")
    d = design_chain(p)
    require_moc_gate(p, d)                        # MOC の不合格は run dir を作る前に止める (2026-10-07)
    scale = float(p.spec["r_throat"])
    # 壁: physical_upstream が無ければ設計壁そのもの (今と同じ)、あれば H の上流を付けた設計壁 (ゲートの不合格は run dir を作る前に例外)
    wall, euler_wall = euler_reference_wall(p, d, scale, pwu=pwu)
    run_dir.mkdir(parents=True, exist_ok=False)
    mp = mesh_params_euler(p, scale)
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
    _check_generated_nodes(coords, nodes, "mesh_euler")
    write_msh41_2d(run_dir / "nozzle.msh", coords, quads, bedges)
    # 記録: 目標軸分布 (x0 → x_E) と設計壁
    law = d["law"]
    xs = np.linspace(d["x0"], d["x_E"], 400)
    tgt = [float(law(x)) if x >= d["x_A"] else float("nan") for x in xs]
    np.savetxt(run_dir / "target_axis_M.csv", np.c_[xs * scale, tgt],
               delimiter=",", header="x_m,M_target", comments="")
    np.savetxt(run_dir / "wall_design.csv",
               np.c_[d["wall_inv"] * [scale, scale, 1.0, 1.0]], delimiter=",",
               header="x_m,r_m,theta_rad,M_wall", comments="")
    n = nsteps or int(p.evaluate.get("nStepOuter", 12000))
    out_int = int(p.evaluate.get("outStepInterval", max(n // 3, 1)))
    (run_dir / "bcondConfig.yaml").write_text(_bcond_with_species(_bcond(p, euler=True), _tp_species_Y(p)))
    (run_dir / "probe.yaml").write_text(PROBE_STUB)
    # 品質検査は cell 変換の一時コピー (品質ツールは node CONNE 非対応)
    (run_dir / "solverConfig.yaml").write_text(
        _apply_gas_to_config(_config_euler(p, n, out_int, 4.0, 1), p, run_dir))
    run_converter([str(converter_path()), "nozzle.msh", "nozzle_qc.h5"], cwd=run_dir, env=_ENV)
    q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"),
                            "nozzle_qc.h5", "--ar-max", str(int(m_eu.get("ar_max", 1000)))], cwd=run_dir, env=_ENV,
                       capture_output=True, text=True)
    (run_dir / "MESH_QUALITY.txt").write_text(
        "# cell 変換コピーで検査 (品質は primal の性質)\n" + q.stdout + q.stderr)
    (run_dir / "nozzle_qc.h5").unlink()
    if q.returncode != 0:
        raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
    # node 変換 (node config を書いてから変換 — 必須) + 等エントロピー IC
    cfg_e = _apply_gas_to_config(_config_euler_node(p, n, out_int, 4.0, 1), p, run_dir)
    # Euler の本段の cfl: 引数 (--cfl) > evaluate.cfl_euler (要求カードの雛形; NS の本段が B0 で cfl_main を書かないので Euler の cfl を
    # 別に持つ、plan tooling-settle-convergence-and-monitor §4.5・codex m55-three-coef M4)。どちらも無ければ今までどおり設定は 4.0 のまま
    # (cfl_main は run_staged の呼び出し側が渡す)
    if cfl_main is None and p.evaluate.get("cfl_euler") is not None:
        cfl_main = euler_cfl(p)
    if cfl_main is not None:
        cfg_e = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", f"cfl: {float(cfl_main)}, cfl_pseudo: {float(cfl_main)}", cfg_e)
    if implicit_relax is not None:
        cfg_e = cfg_e.replace("blockDPLUR: 1,", f"blockDPLUR: 1, implicitRelax: {float(implicit_relax)},", 1)
    (run_dir / "solverConfig.yaml").write_text(cfg_e)
    run_converter([str(converter_path()), "nozzle.msh", "nozzle.h5"], cwd=run_dir, env=_ENV)
    # IC を入れる前 (座標・接続は以後変わらない)。nj は実際の値 (axis_gap_frac では導出値になる)
    mesh_rec = _mesh_record(p, mp, "mesh_euler", coords.shape[0] // mp.ni, run_dir, coords, quads, bedges, nodes)
    paste_isentropic_ic(run_dir / "nozzle.h5", wall, scale,
                        float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp,
                        gas=(None if str(p.evaluate.get('cfd_gas', 'same')) == 'cpg' else p.gas_model),
                        h_ref_T=float(p.evaluate.get('thermo_href_temp', 298.15)),
                        species_Y=_tp_species_Y(p))
    _stamp_ic_species(p, run_dir)          # 新規初期場の化学種属性 (ic_from ならこの後 interp_field が継承/消去を決める)
    if ic_from is not None:
        src = sorted(Path(ic_from).glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
        subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                        str(src), str(run_dir / "nozzle.h5")],
                       env=_ENV, check=True, capture_output=True, text=True)
    info = {"chain": "axismach", "discretization": "node",
            "x0": d["x0"], "x_A": d["x_A"], "x_E": d["x_E"], "L_c": d["L_c"],
            "Lc_window": list(d["Lc_window"]), "anchor": list(d["anchor"]),
            "Lc_mode": d["Lc_mode"], "Lc_solve": d["Lc_solve"],
            "anchor_source": d["anchor_source"], "start_line": d["start_line"],
            "wall_mode": d["wall_mode"], "wall_repr": d["wall_repr"],
            "axis_law": d["axis_law"], "M_knot": d["M_knot"], "x_K": d["x_K"],
            "L_U": float(p.geometry.get("L_U", 3.5)),
            "gas": d["gas"], "species": _species_info(p), "gamma_hall": d["gamma_hall"],
            "wall_fit": d["wall_fit"],
            "Md": d["Md"], "Md_moc_offset": d["Md_moc_offset"], "R": d["R"],
            "initial_line": d["initial_line"],
            "moc": d["moc"],
            "qa": {k: v for k, v in d["qa"].items() if k != "violations"},
            "exit": d["exit"],
            "mdot_ratio_moc": d["mdot_ratio_moc"], "cd_series": d["cd_series"],
            "nStepOuter": n, "scale_m": scale, "ic_from": str(ic_from) if ic_from else None,
            "mesh": mesh_rec, "environment": environment_record()}
    if euler_wall is not None:
        # H の上流を付けた設計壁 (§5.1 #36e の案 A): 壁ファイル (書いて読み直し、復元した壁が今の壁と一致することを確かめる) と記録。
        # キー無しは書かない
        from ..geometry.wall_axismach import save_wall_file
        wpath, wrec, wsha = save_wall_file(run_dir, wall, scale, euler_wall["repr"])
        euler_wall["wall_file"] = {"repr": euler_wall["repr"], "file": wpath.name, "sha256": wsha,
                                   **{k: wrec[k] for k in ("format", "version", "units", "origin", "domain", "domain_m", "physical_wall")}}
        info["physical_upstream"] = euler_wall["physical_upstream"]
        info["euler_wall"] = euler_wall
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


# 段階起動の段 config の変更 (plan tooling-rerun-conditions §4.9、codex result 段 #2)。
# 旧実装は `cfl: [\d.]+, cfl_pseudo: [\d.]+` などの正規表現置換で、`convMethod:  2` (空白 2)・`cfl: 5.0e+0`
# (指数表記)・block 形式の deltaT では黙って当たらず、soft 段が 2 次・CFL 5 のまま回りえた。
# ここでは YAML 上の位置 (_CFG_PATHS) で値を読み、その値トークンだけを書き換え、読み直して実効値と
# 「他の値が変わっていないこと」を検査する (コメント・書式は保つ)。
_CFG_PATHS = {"cfl": ("time", "deltaT", "cfl"), "cfl_pseudo": ("time", "deltaT", "cfl_pseudo"),
              "convMethod": ("space", "convMethod"), "nStepOuter": ("time", "last", "nStepOuter"),
              "outStepInterval": ("time", "outStepInterval"), "nStepInner": ("time", "nStepInner")}
_CFG_INT_KEYS = ("convMethod", "nStepOuter", "outStepInterval", "nStepInner")


def _cfg_num(v):
    """YAML の値を数値として読む (`5e+0` は PyYAML では文字列だが forge は数値として読む)。数値でなければ None。"""
    if v is None or isinstance(v, bool):
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _yaml_strict():
    """重複キーを拒否する YAML ローダー (`solver_density_cuda/tools/yaml_strict.py`)。solver の yaml-cpp は先勝ち・
    PyYAML は後勝ちなので、重複キーのある段 config を PyYAML で検査すると solver と違う値を照合してしまう
    (plan tooling-rerun-conditions、codex result 段 3 回目 #1)。"""
    if str(FORGE_TOOLS) not in sys.path:
        sys.path.insert(0, str(FORGE_TOOLS))
    import yaml_strict
    return yaml_strict


def _cfg_load(cfg_text: str) -> dict:
    """solverConfig を読む (重複キーは yaml_strict.DuplicateKeyError)。辞書でなければ ValueError。"""
    doc = _yaml_strict().load(cfg_text)
    if not isinstance(doc, dict):
        raise ValueError("solverConfig が YAML の辞書でない")
    return doc


def _cfg_get(doc: dict, key: str):
    """YAML 上の位置 _CFG_PATHS[key] の値 (無ければ None)。"""
    v = doc
    for k in _CFG_PATHS[key]:
        if not isinstance(v, dict) or k not in v:
            return None
        v = v[k]
    return v


def _cfg_value(cfg_text: str, key: str):
    """config 本文の key の実効値 (数値)。無い・数値でなければ ValueError。"""
    v = _cfg_num(_cfg_get(_cfg_load(cfg_text), key))
    if v is None:
        raise ValueError(f"solverConfig の {'.'.join(_CFG_PATHS[key])} が無い・数値でない")
    return int(v) if key in _CFG_INT_KEYS else v


def _cfg_set(cfg_text: str, updates: dict) -> str:
    """updates {key: 値} を YAML 上の位置 _CFG_PATHS[key] に書く (yaml_strict.replace_scalars: その位置の値トークンだけを
    置換)。書き換え後に読み直し、(1) 各 key の実効値が要求値、(2) それ以外の値が変わっていない ことを検査する。
    違反は ValueError (黙って当たらない置換を起こさない)。"""
    ys = _yaml_strict()
    before = _cfg_load(cfg_text)
    for key in updates:
        if _cfg_get(before, key) is None:
            raise ValueError(f"solverConfig に {'.'.join(_CFG_PATHS[key])} が無い (段の {key} を書けない)")
    # 値トークンの位置は YAML の構造から決める (正規表現だとコメント・引用符付きキー・`key :` に当たり外れが出る;
    # codex result 段 4 回目 #3)。replace_scalars が読み直して要求値・他の値の不変を検査する
    toks = {_CFG_PATHS[k]: (str(int(v)) if k in _CFG_INT_KEYS else repr(float(v))) for k, v in updates.items()}
    text = ys.replace_scalars(cfg_text, toks)
    after = _cfg_load(text)
    for key, val in updates.items():
        got = _cfg_num(_cfg_get(after, key))
        if got is None or got != float(val):
            raise ValueError(f"段の config を読み直したら {'.'.join(_CFG_PATHS[key])} = {_cfg_get(after, key)!r} (要求 {val})")
    return text


def _cfg_check(cfg_text: str, want: dict, label: str) -> None:
    """起動前の検査: config の実効値 (YAML 上の位置で読む) が want {key: 要求値} と一致しなければ RuntimeError。"""
    doc = _cfg_load(cfg_text)
    bad = {k: _cfg_get(doc, k) for k, v in want.items() if _cfg_num(_cfg_get(doc, k)) != float(v)}
    if bad:
        raise RuntimeError(f"段 {label} の config の実効値が要求と違う — forge を起動しない: "
                           + ", ".join(f"{'.'.join(_CFG_PATHS[k])} = {bad[k]!r} (要求 {want[k]})" for k in bad))


def run_staged(run_dir, cfl_main: float | None = None, mid_stage: bool = False, stages: str = "full",
               mesh_h5: str = "nozzle.h5") -> int:
    """soft 段 (1次+cfl0.5, 3000 step) → [mid 段 (2次+cfl1, 3000 step)] → 本段。
    walldriven W3 と同方式・段階起動必須。`cfl_main` で本段 CFL を上書き
    (semi-perfect TP は cfl4 で本段 step ~60 に爆発 — case/42 実測。TP は cfl≤2 が実績
    [[cutler-cpg-vs-tp-dplur-sst]])。`mid_stage=True` で 2 次化と CFL 上げを分離する。
    `run_staged_ns` と同じく、段ごとの実効設定を `stage_manifest.json` に記録し、段の残差履歴を
    `residual_history_<tag>.csv` に残す (tag = soft / mid / main)。段の最終 res の必要保存量が非有限・ρ≤0 なら
    次段へ進まず RuntimeError (段終了ゲート `stage_gate`; plan tooling-rerun-conditions §4.9、codex result 段 3 回目 #3)。"""
    import shutil
    run_dir = Path(run_dir)
    cfg_main = (run_dir / "solverConfig.yaml").read_text()
    bc_path = run_dir / "bcondConfig.yaml"
    bc_text = bc_path.read_text() if bc_path.exists() else ""
    if cfl_main is not None:
        cfg_main = _cfg_set(cfg_main, {"cfl": cfl_main, "cfl_pseudo": cfl_main})
    main_want = {k: _cfg_value(cfg_main, k) for k in ("convMethod", "cfl", "cfl_pseudo", "nStepOuter", "outStepInterval")}
    sm = _stage_manifest_cls()(run_dir)

    def _record(cfg, tag):
        """段の残差履歴を段名つきで残し、manifest に段を足して書く (失敗した段も診断用に残す; run_staged_ns と同じ)。"""
        hist = run_dir / "residual_history.csv"
        if hist.exists():
            shutil.copy(hist, run_dir / f"residual_history_{tag}.csv")
        sm.add(tag, cfg, bc_text, history=f"residual_history_{tag}.csv")
        sm.write()

    def _stage(cfg, nsteps, label, want):
        cfg = _cfg_set(cfg, {"nStepOuter": nsteps, "outStepInterval": nsteps})
        _cfg_check(cfg, {**want, "nStepOuter": nsteps, "outStepInterval": nsteps}, label)
        (run_dir / "solverConfig.yaml").write_text(cfg)
        (run_dir / "residual_history.csv").unlink(missing_ok=True)   # 前段の履歴を別段の名前で写さない
        rc = run_forge(run_dir)
        _record(cfg, label)
        res = sorted(run_dir.glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        if rc != 0 or not res or int("".join(c for c in res[-1].stem if c.isdigit())) < nsteps:
            raise RuntimeError(f"{label} 段が失敗 (発散切り分けは res_nan_*.h5 を見る)")
        probs = stage_gate(res[-1], cfg)                    # restart_field のビット一致検査は Inf・負密度を排除しない
        if probs:
            raise RuntimeError(f"段 {label} の最終場 {res[-1].name} が段終了ゲートで不合格 — 次段へ進まない:\n    "
                               + "\n    ".join(probs))
        _restart_same_mesh(res[-1], run_dir / mesh_h5)      # 同一メッシュ: index コピー (旧: interp_field.py)
        for f in run_dir.glob("res_*"):
            f.unlink()

    if stages != "none":
        soft = _first_order(_cfg_set(cfg_main, {"cfl": 0.5, "cfl_pseudo": 0.5}))   # 旧: `convMethod: 1` だけを 0 に置換
        _stage(soft, 3000, "soft", {"convMethod": 0, "cfl": 0.5, "cfl_pseudo": 0.5})
        if mid_stage and stages == "full":
            mid = _cfg_set(cfg_main, {"cfl": 1.0, "cfl_pseudo": 1.0})
            _stage(mid, 3000, "mid", {"convMethod": main_want["convMethod"], "cfl": 1.0, "cfl_pseudo": 1.0})
    _cfg_check(cfg_main, main_want, "main")
    (run_dir / "solverConfig.yaml").write_text(cfg_main)
    (run_dir / "residual_history.csv").unlink(missing_ok=True)
    rc = run_forge(run_dir)
    _record(cfg_main, "main")
    return rc


def collect(problem_path, run_dir) -> dict:
    """軸 M vs 目標 (x_A ≤ x ≤ x_E−0.3 マスク) + 出口一様性 + overshoot。"""
    from ..metrics.extract import axis_mach, exit_uniformity

    p = load_problem(problem_path)
    run_dir = Path(run_dir)
    res = sorted(run_dir.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
    info = json.loads((run_dir / "prepare_info.json").read_text())
    scale = info["scale_m"]
    Md = float(info["Md"])
    x_ach, M_ach = axis_mach(run_dir / "nozzle.h5", res[-1],
                             axis_band=0.9 * scale * 0.1)
    tgt = np.loadtxt(run_dir / "target_axis_M.csv", delimiter=",", skiprows=1)
    ok = np.isfinite(tgt[:, 1])
    xs = np.linspace(max(info["x_A"], tgt[ok, 0].min() / scale) * scale + 1e-9,
                     (info["x_E"] - 0.3) * scale, 200)
    Mt = np.interp(xs, tgt[ok, 0], tgt[ok, 1])
    Ma = np.interp(xs, x_ach, M_ach)
    dM = Ma - Mt
    # overshoot: x_E 以降〜出口手前の軸 M の最大値
    m_dn = x_ach > info["x_E"] * scale
    M_max_dn = float(np.max(M_ach[m_dn])) if m_dn.any() else float("nan")
    try:
        uni = exit_uniformity(run_dir / "nozzle.h5", res[-1], Md,
                              x_d=info["x_E"] * scale, gamma=p.gamma)
    except Exception as e:  # noqa: BLE001 — 一様性は診断 (核心は軸 M)
        uni = {"error": str(e)}
    mf = None
    if info.get("euler_ref"):
        try:
            from ..metrics.deltastar import massflow_ratio
            # 一定 x の断面で積分 (plan tooling-nozzle-card-local-mesh §5.1 #3i): 列のままは壁法線の近壁層の格子で断面の流量でない
            # (case/71 run_0005_ns: 列 120.19 kg/s・散らばり 2.1 % → 一定 x 121.350 kg/s・6.4e−5)。列が一定 x の格子では値は同じ
            mf = massflow_ratio(run_dir, info["euler_ref"], remap_constant_x=True)
        except Exception as e:  # noqa: BLE001 — 帳簿は診断
            mf = {"error": str(e)}
    out = {"res_file": res[-1].name,
           "massflow": mf,
           "dM_max": float(np.max(np.abs(dM))),
           "dM_max_rel_Md": float(np.max(np.abs(dM)) / Md),
           "dM_rms": float(np.sqrt(np.mean(dM ** 2))),
           "M_axis_exit": float(np.interp(info["x_E"] * scale, x_ach, M_ach)),
           "M_axis_max_downstream": M_max_dn,
           "overshoot_rel": float((M_max_dn - Md) / Md) if np.isfinite(M_max_dn) else None,
           "exit_uniformity": uni,
           "mdot_ratio_moc": info["mdot_ratio_moc"]}
    np.savetxt(run_dir / "achieved_vs_target.csv", np.c_[xs, Mt, Ma],
               delimiter=",", header="x_m,M_target,M_achieved", comments="")
    (run_dir / "metrics.json").write_text(json.dumps(out, indent=1))
    return out


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="axis-Mach チェーン評価 (node Euler)")
    ap.add_argument("problem")
    ap.add_argument("run_dir")
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--cfl", type=float, default=None, help="本段 cfl (YAML evaluate.cfl_main を上書き)")
    ap.add_argument("--implicit-relax", type=float, default=None, help="implicitRelax を deltaT に挿入")
    ap.add_argument("--stages", default="full", choices=("full", "soft", "none"), help="full=soft/mid/本段, soft=soft+本段, none=本段のみ")
    ap.add_argument("--ic-from", default=None)
    a = ap.parse_args(argv)
    info = prepare(a.problem, a.run_dir, nsteps=a.steps, ic_from=a.ic_from, cfl_main=a.cfl, implicit_relax=a.implicit_relax)
    info["stages"] = a.stages
    (Path(a.run_dir) / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    print(json.dumps(info, indent=1, default=str))
    if a.prepare_only:
        return 0
    pp = load_problem(a.problem)
    import time as _t
    t0 = _t.time()
    rc = run_staged(a.run_dir, cfl_main=(a.cfl if a.cfl is not None else euler_cfl(pp)),     # cfl_euler > cfl_main
                    mid_stage=bool(pp.evaluate.get("mid_stage", pp.is_semiperfect)), stages=a.stages)
    print(f"forge exit={rc} (wall {_t.time() - t0:.0f} s, stages={a.stages})")
    print(json.dumps(collect(a.problem, a.run_dir), indent=1))
    return rc



# --- A12: 粘性 δ* 補正 (RANS 経路) ------------------------------------------------
def delta_r_from_table(x, d):
    """δ_r の表 (P-spline 平滑化済みの `delta_r_next.csv`) を壁に渡す関数にする: 5 次補間スプライン、範囲外は端値。

    2026-10-05 (plan verification-m6-axis-wave-mesh-su2 §5.1 #8a): 旧実装は np.interp (直線補間) で、表の点ごと
    (0.09 r_t) の傾きの折れ目を補間 5 次 B-spline の壁が通るため r″ が点間隔の周期で波打っていた (高周波 3〜6e-4 [1/r_t])。
    平滑化前の生値を渡すと 5 次補間はリンギングするので、平滑化済みの表に限る。

    `f(xq, deriv)` で導関数 (deriv=1..3) も返す (plan tooling-nozzle-cfd-pinned-initial-line §5.1 #6b: joint 壁の
    物理壁は r′・r″・r‴ を解析的に足す)。表の範囲外は値が端値クリップなので導関数は 0。"""
    from scipy.interpolate import make_interp_spline
    x = np.asarray(x, dtype=float); d = np.asarray(d, dtype=float)
    spl = make_interp_spline(x, d, k=5)
    lo, hi = float(x[0]), float(x[-1])

    def f(xq, deriv: int = 0, _s=spl, _lo=lo, _hi=hi):
        xq = np.asarray(xq, dtype=float)
        if deriv == 0:
            return _s(np.clip(xq, _lo, _hi))
        if deriv not in (1, 2, 3):
            raise ValueError("deriv は 0..3")
        return np.where((xq >= _lo) & (xq <= _hi), _s(np.clip(xq, _lo, _hi), deriv), 0.0)
    f.supports_deriv = True
    f.spline = spl          # ノット (区分の境界) を形状の厳密評価に渡すため (plan tooling-nozzle-throat-monotone-r2 §6 S6)
    f.x_range = (lo, hi)    # 表の範囲 (1 本の B-spline の物理壁は、ランプ開始〜出口を覆うことを要求する; plan tooling-nozzle-wall-single-bspline §4.1)
    return f


def _closure_kw(p: Problem, init_cfg: dict) -> dict:
    """deltastar_initializer.closure_version (既定 contur_v1 = 従来どおり、何も渡さない)。contur_v2 は版と NS と同じ μ(T) を渡す。"""
    v = str(init_cfg.get("closure_version", "contur_v1"))
    if v == "contur_v1":
        return {}
    if v != "contur_v2":
        raise ValueError(f"deltastar_initializer.closure_version {v!r} は contur_v1 | contur_v2")
    return {"closure_version": v, "mu_fn": contur_mu_fn(p)}


def contur_mu_fn(p: Problem):
    """積分法 `contur_v2` の粘性 μ(T): NS と同じ輸送物性 (plan tooling-nozzle-isothermal-wall-chain §4.7-3)。
    多成分 TP の CFD (viscMethod 2) なら gas.transport の混合気 (無ければ例外)、それ以外は NS の Sutherland (1.716e-5, 273, 111)。"""
    if p.uses_tp_cfd:
        tr = p.gas_transport
        if not tr:
            raise ValueError("contur_v2: TP の問題は gas.transport (NS と同じ輸送物性) が要る — Sutherland に落とさない")
        from ..gas.transport import MixtureViscosity
        Y = p.gas_composition[0]
        return MixtureViscosity({str(k).upper(): float(v) for k, v in Y.items()}, tr)
    return lambda T: 1.716e-5 * (np.asarray(T, dtype=float) / 273.0) ** 1.5 * (273.0 + 111.0) / (np.asarray(T, dtype=float) + 111.0)


def integral_delta_r(p: Problem, d: dict, init_cfg: dict, scale: float | None = None, rtol: float | None = None):
    """積分法初期壁の δ_r (`prepare_ns` の initializer 経路): `integral_bl` → 5 次 P-spline 平滑化 → 壁に渡す δ_r 関数。
    戻り: (res_init, delta_r_x, init_info)。`prepare_ns` から切り出したもの (振る舞いは同一; plan
    tooling-nozzle-throat-monotone-r2 §6 S6 の形状ゲートが同じ経路で物理壁を作るために共有する)。
    scale: スロート半径 r_t [m] (None = spec.r_throat)。寸法の逆算 (`deltastar_loop.solve_rt`・`solve_rt_throat`) が、反復のたびに
    同じ経路 (k_f の cf_scale・熱条件・平滑化) で δ_r を作り直すために渡す (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.2)。
    rtol: **診断用** (plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U2c)。`integral_bl` の RK45 の相対許容差を明示の引数で
    注入する。None = 渡さない (`integral_bl` の既定 1e-6、今の振る舞いのまま)。YAML のキーにはしない。実効値は
    `res_init["solve_ivp"]["rtol"]` (solve_ivp に渡した値) に残る。
    `geometry.physical_upstream` がある問題は ValueError (積分法は設計壁の縮流部を入口から前進積分するので、物理壁の別の形の上流の
    境界層の履歴にならない。未対応のまま黙って設計壁の縮流部で数えない; plan tooling-nozzle-isothermal-wall-chain §5.1 #34)。"""
    from ..feedback.deltastar_integral import integral_bl
    if _physical_upstream(p.geometry) is not None:
        raise ValueError("積分法の初期壁 (deltastar_initializer) と geometry.physical_upstream の併用は未対応 — 積分法は設計壁の縮流部を "
                         "前進積分するので物理壁の上流と別の形になる。δ_r は delta_r_csv で与える")
    scale = float(p.spec["r_throat"]) if scale is None else float(scale)
    wall_inv = d["wall_inv"]
    model = str(init_cfg.get("model", "contur"))
    if model not in ("contur", "contur_momentum_integral"):
        raise ValueError(f"deltastar_initializer.model={model!r} は未対応 (contur のみ)")
    # 熱境界条件は spec.wall_thermal が単一ソース (plan tooling-nozzle-isothermal-wall-chain §4.1)。
    # initializer/YAML の thermal_bc 指定は無視し、食い違えば警告する (NS と積分法が別の壁温を読む状態を作らない)。
    # 壁温分布の表は**ここで解決した同じ scale** で換算する (integral_bl に渡す scale と同じ; 寸法の逆算の試しの r_t を含む。
    # plan tooling-nozzle-wall-temperature-profile §4.2 の契約)。一様・断熱は wall_thermal_bc_integral と同じ dict
    tbc = integral_thermal_bc(p, d, scale)
    if init_cfg.get("thermal_bc") and dict(init_cfg["thermal_bc"]) != tbc:
        print(f"[prepare_ns] warning: initializer.thermal_bc={init_cfg['thermal_bc']} は無視 (spec.wall_thermal={p.wall_thermal} を使用)")
    res_init = integral_bl(d["wall"], wall_inv, _gam_or_gas(p), p.cp, float(p.spec["Pt"]), float(p.spec["Tt"]),
                           scale, thermal_bc=tbc,
                           theta0_m=init_cfg.get("theta0_m"), x_virtual_m=init_cfg.get("x_virtual_m"),
                           a_crocco=float(init_cfg.get("a_crocco", 1.0)), closure=str(init_cfg.get("closure", "contur")),
                           cf_scale=float(init_cfg.get("cf_scale", 1.0)), n_scale=float(init_cfg.get("n_scale", 1.0)),
                           **({} if rtol is None else {"rtol": float(rtol)}), **_closure_kw(p, init_cfg))
    # 積分法の出力も同じ 5 次 P-spline で平滑化 (N(Re) テーブルの折れ目などを壁曲率に持ち込まない)
    from ..metrics.deltastar import smooth_delta_quintic
    f_s, sm_diag = smooth_delta_quintic(res_init["x"], res_init["delta_r"], knot_spacing=2.0, lam=1.0,
                                        positive=(p.wall_thermal["mode"] == "adiabatic"))  # 等温は符号付き
    res_init["delta_r_raw_integral"] = res_init["delta_r"].copy()
    res_init["delta_r"] = f_s(res_init["x"])
    # 2 係数の較正 (2026-10-15, plan tooling-nozzle-card-local-mesh §4.3、ユーザ決定): δ = c·δ_C + δ₀/scale (r_t 単位; c = delta_scale、
    # δ₀ = delta_offset_m [m] を全長一定)。NS の δ_E のスロートと出口の 2 点に合わせた CONTUR の倍率とずらし (負の δ も表せる)。
    # 実効の scale (寸法の逆算の試しの r_t を含む) で換算し、値・導関数・表・記録を同じ変換の結果から作る。既定 (1, 0) は何もしない
    # 3 係数 (同 §4.3 の末尾、ユーザ決定「bがよさそう」): δ = (c + c₁·x_m)·δ_C + δ₀/scale (x_m = x·scale [m]、c₁ = delta_scale_slope_per_m
    # [1/m]、既定 0 = 2 係数とビット同一の式)
    two = _two_coef(init_cfg)
    if two is not None:
        c2, d0_m = two
        c1 = _slope_coef(init_cfg)
        f_C = f_s
        res_init["delta_r_contur"] = res_init["delta_r"].copy()          # 係数を当てる前の CONTUR (再較正の基底)
        if c1 == 0.0:
            def f_s(xx, _f=f_C, _c=c2, _d=d0_m / scale):
                return _c * _f(xx) + _d
        else:
            def f_s(xx, _f=f_C, _c=c2, _c1=c1 * scale, _d=d0_m / scale):
                xx = np.asarray(xx, dtype=float)
                return (_c + _c1 * xx) * _f(xx) + _d
        res_init["delta_r"] = f_s(res_init["x"])
    # 壁には平滑化関数そのものを渡す (2026-10-05, plan verification-m6-axis-wave-mesh-su2 §5.1 #8a):
    # 1500 点の表を np.interp で渡すと点ごとの傾きの折れ目を補間 5 次スプラインが通り、r″ が点間隔で波打つ
    delta_r_x = f_s
    if getattr(d["wall"], "wall_repr", None) == "joint" or type(d["wall"]).__name__ == "JointFitCFDWall":
        # joint 壁の物理壁 (解析経路) は δ_r の導関数を要る (plan tooling-nozzle-cfd-pinned-initial-line §5.1 #6b)。
        # 平滑化済みの値 (5 次 P-spline) を表にし、導関数を返せる 5 次補間 (delta_r_from_table) で渡す。
        delta_r_x = delta_r_from_table(res_init["x"], f_s(res_init["x"]))
    init_info = dict(res_init["settings"])
    init_info["smooth"] = {"kind": "quintic_pspline", **sm_diag}
    init_info["delta_r_throat"] = float(np.interp(0.0, res_init["x"], res_init["delta_r"]))
    init_info["delta_r_exit"] = float(res_init["delta_r"][-1])
    if two is not None:
        init_info["delta_scale"], init_info["delta_offset_m"] = two
        if c1 != 0.0:
            init_info["delta_scale_slope_per_m"] = c1
        init_info["two_coef"] = {"delta_scale": two[0], "delta_offset_m": two[1], "scale_m": scale,
                                 **({"delta_scale_slope_per_m": c1} if c1 != 0.0 else {}),
                                 "delta_r_throat_contur": float(np.interp(0.0, res_init["x"], res_init["delta_r_contur"])),
                                 "delta_r_exit_contur": float(res_init["delta_r_contur"][-1]),
                                 "plan": "plans/active/tooling-nozzle-card-local-mesh.md §4.3"}
    return res_init, delta_r_x, init_info


def _slope_coef(init_cfg: dict) -> float:
    """`deltastar_initializer.delta_scale_slope_per_m` (3 係数の c₁ [1/m]、既定 0) を読む。数でない・非有限は ValueError。"""
    v = init_cfg.get("delta_scale_slope_per_m", 0.0)
    try:
        v = float(v)
    except (TypeError, ValueError):
        raise ValueError(f"deltastar_initializer.delta_scale_slope_per_m は数 ({init_cfg.get('delta_scale_slope_per_m')!r})") from None
    if not np.isfinite(v):
        raise ValueError(f"deltastar_initializer.delta_scale_slope_per_m は有限 ({v!r})")
    return v


def _two_coef(init_cfg: dict):
    """`deltastar_initializer` の 2 係数 (delta_scale, delta_offset_m) を読む。どちらも既定 (1, 0) で c₁ (delta_scale_slope_per_m) も 0 なら
    None (何もしない)。数でない・非有限・c ≤ 0 は ValueError (黙って既定に戻さない)。"""
    c = init_cfg.get("delta_scale", 1.0); d = init_cfg.get("delta_offset_m", 0.0)
    try:
        c, d = float(c), float(d)
    except (TypeError, ValueError):
        raise ValueError(f"deltastar_initializer.delta_scale・delta_offset_m は数 ({init_cfg.get('delta_scale')!r}, "
                         f"{init_cfg.get('delta_offset_m')!r})") from None
    if not (np.isfinite(c) and np.isfinite(d)) or c <= 0.0:
        raise ValueError(f"deltastar_initializer.delta_scale は正の有限・delta_offset_m は有限 ({c!r}, {d!r})")
    return None if (c == 1.0 and d == 0.0 and _slope_coef(init_cfg) == 0.0) else (c, d)


def build_physical_wall(p: Problem, d: dict, scale: float, delta_r_x=None, dstar_x=None, offset: str = "normal",
                        pwu: dict | None = None):
    """物理壁の構築 (`prepare_ns` と寸法の逆算 `deltastar_loop.solve_rt`・`solve_rt_throat` が共有する; plan
    tooling-nozzle-upstream-poly-and-throat-sizing §4.2「同じ壁の構築」)。`PhysicalNozzleWall` に `geometry.pw_ramp`・
    `geometry.pw_upstream` (`_pw_upstream` の解決済みの値) を渡し、`physical_wall_repr: single_bspline` なら 1 本の B-spline に
    作り直す。scale: r_t [m] (相関 δ* の診断に使う。壁の形は r_t 単位)。
    `geometry.physical_upstream` (`_physical_upstream`、plan tooling-nozzle-isothermal-wall-chain §5.1 #34) があれば、物理壁の直管と
    多項式 Q の始点だけをその値にする (設計壁 `d["wall"]` は凍結した形のまま)。キー無しは今と同じ (引数も渡さない)。"""
    from ..geometry.wall_axismach import PhysicalNozzleWall, SingleBSplinePhysicalWall
    pwu = _pw_upstream(p.geometry) if pwu is None else pwu
    # joint 壁の物理壁 (解析経路) の δ_r ランプ区間: geometry.pw_ramp (None = 既定、default_pw_ramp)。ゲート不合格は例外
    pw_ramp = p.geometry.get("pw_ramp")
    pu = _physical_upstream(p.geometry)
    pu_kw = {} if pu is None else {"physical_upstream": {"r_U": pu["r_inlet"], "L_U": pu["L_U"], "L_pipe": pu["L_pipe"]}}
    wall = PhysicalNozzleWall(d["wall"], d["wall_inv"], float(scale), float(p.spec["Pt"]),
                              float(p.spec["Tt"]), _gam_or_gas(p), p.cp, dstar_x=dstar_x,
                              offset=offset, delta_r_x=delta_r_x,
                              ramp=(None if pw_ramp is None else tuple(float(v) for v in pw_ramp)),
                              upstream=pwu["value"], **pu_kw)
    req = _physical_wall_repr(p.geometry)
    if req == "single_bspline":
        # 今の物理壁 (pw_upstream poly) を全域 1 本の 5 次 B-spline にノット挿入で作り直す (許容誤差・ゲート不合格は例外)
        wall = SingleBSplinePhysicalWall(wall)
    elif req is None and pwu["value"] == "poly":
        # キー無しの poly の壁は既定で 1 本の B-spline (2026-10-07 ユーザ決定、plan tooling-nozzle-wall-single-bspline §9)。
        # 1 本で表せない構成なら区分表現のまま作り、理由を壁に残す (prepare_ns が prepare_info に記録する。黙って落とさない)
        why = SingleBSplinePhysicalWall.applicability(wall)
        if why is None:
            wall = SingleBSplinePhysicalWall(wall)
        else:
            wall.single_bspline_default_skipped = why
    return wall


def delta_r_zero(x_lo: float, x_hi: float):
    """補正 δ_r ≡ 0 の関数 (`delta_r_from_table` と同じ形: 値と導関数・`spline`・`x_range`)。#36b では Euler の参照の物理壁が NS の物理壁と
    同じ構築 (`build_physical_wall`) を δ_r = 0 で通るために使った (§5.1 #36e の案 A からは `euler_reference_wall` は使わない。試験・診断で
    「δ_r = 0 の物理壁は今も物理壁のゲートで落ちる」ことを確かめるのに残す)。係数がすべて 0 の 5 次の
    補間スプラインなので値も導関数も厳密に 0 (S + 0 は S とビット同一)。表の範囲 [x_lo, x_hi] は壁の定義域を渡す
    (1 本の B-spline の壁は表が [0, x_e] を覆うことを要求する)。"""
    xt = np.linspace(float(x_lo), float(x_hi), 7)
    f = delta_r_from_table(xt, np.zeros_like(xt))
    f.is_zero = True
    return f


def euler_reference_wall(p: Problem, d: dict, scale: float, pwu: dict | None = None):
    """Euler (`prepare`) の格子の壁 (plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #36b、§5.1 #36e で案 A に)。戻り (wall, rec)。

    - `geometry.physical_upstream` が無い: 設計壁 `d["wall"]` そのもの (今と同じ)。rec は None。
    - ある: **H の上流を付けた設計壁** (`UpstreamDesignWall`; codex diagnose 3 `notes/reviews/2026-10-14-h370-registration-diagnose.md`
      の案 A)。x ≥ 0 は凍結した設計壁 S そのもの、x < 0 は直管 r_inlet と (−L_U, r_inlet, 0, 0) → x = 0 の S の (r, r′, r″) の 5 次
      多項式 (物理壁の `poly` の上流と同じ式)。δ* ループの抽出器は δ_in = NS の壁 − この壁 と定義するので、境界層の補正を足さない固定の
      Euler の参照に使う (NS の補正済みの物理壁を slip にした壁ではない — δ_in の基準が変わる; codex diagnose 2026-10-14)。
      ゲート: 下流 (x ≥ 0) は**設計壁のゲート** `monotone_gate`、上流は単調な収縮・正の半径、全域の一意な大域最小、x = −L_U・0 の C²
      (`UpstreamDesignWall`)。**物理壁のクラス (`PhysicalNozzleWall`) とその厳密なゲート (`THROAT_R1_TOL`) は使わず、緩めもしない**
      (#36b の「S の出口端の微小な減少が物理壁のゲートで落ちる」判断待ちは案 A で決着: この壁は物理壁ではない)。設計壁と CFD ピンの
      照合は geometry の r_inlet・L_U・L_pipe のまま (`design_chain` が照合する。照合を外さない)。ゲート・validate の不合格は ValueError。
    pwu: `_pw_upstream` の解決済みの値 (poly だけ。記録に使う)。
    rec: prepare_info.json の `euler_wall` (種類・δ_r = 0・上流の形とゲート・スロート・壁ファイルの種類 `design_upstream`)。
    壁ファイルは呼び出し側が書く。"""
    from ..geometry.wall_axismach import DESIGN_UPSTREAM_KIND, UpstreamDesignWall
    pu = _physical_upstream(p.geometry)
    if pu is None:
        return d["wall"], None
    pwu = _pw_upstream(p.geometry) if pwu is None else pwu
    if pwu["value"] != "poly":
        raise ValueError(f"Euler の参照の壁 (physical_upstream) は pw_upstream poly だけ (受け取った {pwu['value']!r})")
    wall = UpstreamDesignWall(d["wall"], {"r_U": pu["r_inlet"], "L_U": pu["L_U"], "L_pipe": pu["L_pipe"]})
    msgs = wall.validate()
    if msgs:
        raise ValueError("Euler の参照の壁 (H の上流を付けた設計壁) が不合格: " + "; ".join(msgs))
    rec = {"kind": "design_with_physical_upstream", "delta_r": 0.0,
           "note": ("H の上流を付けた設計壁 (案 A、§5.1 #36e): x ≥ 0 は設計壁 S、x < 0 は physical_upstream の直管と S(0) への 5 次多項式。"
                    "物理壁ではない (境界層の補正なし)"),
           "physical_upstream": {"requested": pu, "wall": wall.physical_upstream},
           "pw_upstream": {"value": wall.pw_upstream, "source": pwu["source"], "requested": pwu["requested"]},
           "gate": wall.gate,
           "throat": {"x": float(wall.x_throat), "r": float(wall.r_throat), "kappa": float(wall.kappa_throat)},
           "repr": DESIGN_UPSTREAM_KIND}
    return wall, rec


# --- 本段だけの陰解法 B0 (`evaluate.main_line`; plan tooling-nozzle-isothermal-wall-chain §4.8-2・§5.1 #30) ---------------
# case/45 の生産の陰解法 B0 (壁法線のライン陰解法・方向別の擬似 dt とその伸びの上限・implicitThermalJacobian・cfl・implicitRelax・
# リミッタの基準値の固定・出力の量) を、段階起動の前の段 (S1・S2) には付けずに本段だけに書く。キーと書き方は case/45 の
# `cold_cfl.py prep --line dir --itj … --cap … --cfl … --limiter-ref-from … --extra …` と同じ (implicitRelax だけは prepare_ns の
# `implicit_relax` と同じく blockDPLUR の直後に挿入する。B0 の系譜 run_0182 → run_0183 → run_0483 の本段の書き方)。
# 前の段の元の設定 (B0 を当てる前) は run の PRESTAGE_CONFIG に置き、run_staged_ns が S1・S2 をそこから作る (本段の設定から
# CFL・次数・nStepInner だけを変えて作ると S1 にもライン陰解法・上限・ITJ が残る; codex diagnose 2026-10-13)。
PRESTAGE_CONFIG = "solverConfig_prestage.yaml"
MAIN_STAGE_STARTED = "MAIN_STAGE_STARTED"
_MAIN_LINE_KEYS = ("line", "cap", "itj", "cfl", "relax", "limiter_ref", "extra_fields")
# limiter_ref の省略 = 初期場の元が NS の run ならその基準値を継承、それ以外は自動 (forge が本段の開始の場で決めて run の limiter_ref.json に
# 書く; plan tooling-settle-convergence-and-monitor §4.6・§5.1 #6・#6b、2026-10-16)。`auto` の明示は継承しない自動。extra_fields (出力の量) も任意
_MAIN_LINE_REQUIRED = ("line", "cap", "itj", "cfl", "relax")
LIMITER_REF_AUTO = "auto"                                     # main_line.limiter_ref: auto (初期場の元を見ずに自動)
_LIMITER_REF_KEYS = (("ro", "limiterRoRef"), ("p", "limiterPRef"), ("a", "limiterARef"))
_LIMITER_REF_OPT_KEYS = (("L", "limiterRefLength"),)          # 明示の数値で任意 (省略 = 格子の境界箱; from_run は 4 つとも写す)
# B0 だけが書くパス (前の段の設定に残っていてはいけない)。cfl・cfl_pseudo・implicitRelax は前の段にもありうるので含めない
_MAIN_LINE_ONLY_PATHS = ({("time", "deltaT", k) for k in ("lineImplicit", "lineDtDirectional", "lineDtDirectionalCap",
                                                         "implicitThermalJacobian")}
                         | {("space", k) for _, k in _LIMITER_REF_KEYS + _LIMITER_REF_OPT_KEYS} | {("output",)})
# 前の段の元の設定と本段の設定の差として許すパス (B0 のパス + cfl・cfl_pseudo・implicitRelax)
_MAIN_LINE_PATHS = (_MAIN_LINE_ONLY_PATHS
                    | {("time", "deltaT", k) for k in ("cfl", "cfl_pseudo", "implicitRelax")})
_CFG_MISSING = object()


def _cfg_diff_paths(a, b, path=()) -> dict:
    """読み込んだ 2 つの設定の差 {パス: (a の値, b の値)}。片側にしか無いキーは _CFG_MISSING
    (case/45 の make_ns_n012_problems.diff_paths と同じ定義)。"""
    out = {}
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            out.update(_cfg_diff_paths(a.get(k, _CFG_MISSING), b.get(k, _CFG_MISSING), path + (k,)))
        return out
    if a is _CFG_MISSING or b is _CFG_MISSING or type(a) is not type(b) or a != b:
        out[path] = (a, b)
    return out


def _cfg_path_get(doc, path):
    """読み込んだ設定の path (キーの tuple) の値。無ければ _CFG_MISSING。"""
    v = doc
    for k in path:
        if not isinstance(v, dict) or k not in v:
            return _CFG_MISSING
        v = v[k]
    return v


def _main_line_mode(p: Problem) -> dict:
    """本段の陰解法の選び方の記録 (prepare_info の `main_line_mode`)。`evaluate.main_line` が mapping → B0、`point` (文字列) → point の
    明示 (本段の cfl は `cfl_main`; 要求カードの雛形の既定 B0 から外すとき)、キー無し・null → point (今までどおり)。
    それ以外の値は ValueError (`_main_line_spec` と同じ)。"""
    ml = p.evaluate.get("main_line")
    if isinstance(ml, dict):
        return {"value": "B0", "source": "evaluate.main_line"}
    if ml == "point":
        return {"value": "point", "source": "evaluate.main_line: point (明示)"}
    if ml is None:
        return {"value": "point", "source": "evaluate.main_line が無い (今までどおり)"}
    raise ValueError(f"evaluate.main_line は mapping (B0) か 'point' であること ({ml!r})")


def _limiter_tools():
    """solver_density_cuda/tools/stage_manifest.py (limiter_ref.json の読み・検査と fnv1a64)。"""
    if str(FORGE_TOOLS) not in sys.path:
        sys.path.insert(0, str(FORGE_TOOLS))
    import stage_manifest
    return stage_manifest


def limiter_ref_from_run(run, base_dir=None) -> tuple:
    """別 run のリミッタの基準値の記録を継承する (`main_line.limiter_ref: {from_run: <run>}`; plan tooling-settle-convergence-and-monitor
    §4.6)。戻り ({limiterRoRef, limiterPRef, limiterARef, limiterRefLength: float}, 記録)。
    記録の正本は forge が起動ごとに書く `<run>/limiter_ref.json` (%.17g)。その cfg_fnv が `<run>/solverConfig.yaml` (最後に回した設定 =
    run_staged_ns の run なら本段) と一致すること、基準値を使った起動であること (computed)、`<run>/stage_manifest.json` があれば最後の段の
    記録と値が同じことを確かめる (違えば ValueError; 黙って別の値を継承しない)。相対パスは base_dir (問題 YAML の場所) から。"""
    SM = _limiter_tools()
    src = Path(str(run)).expanduser()
    if not src.is_absolute() and base_dir is not None:
        src = Path(base_dir) / src
    cfgp = src / "solverConfig.yaml"
    if not cfgp.is_file():
        raise ValueError(f"main_line.limiter_ref.from_run: {cfgp} が無い")
    try:
        rec = SM.read_limiter_ref(src, cfgp.read_text())
    except SM.LimiterRecordError as e:
        raise ValueError(f"main_line.limiter_ref.from_run: 継承する記録を読めない — {e}") from e
    if not rec["computed"]:
        raise ValueError(f"main_line.limiter_ref.from_run: {src} は基準値を使わない起動の記録 (limiterScaled 0) — 継承できない")
    vals = {k: float(rec[k]) for k in SM.LIMITER_REF_NAMES}
    man = src / "stage_manifest.json"
    mrec = None
    if man.is_file():
        st = (json.loads(man.read_text()).get("stages") or [])
        mrec = (st[-1].get("limiter_ref") if st else None)
        if isinstance(mrec, dict) and mrec.get("values") is not None and mrec.get("values") != vals:
            raise ValueError(f"main_line.limiter_ref.from_run: {src} の limiter_ref.json {vals} と stage_manifest の最後の段 "
                             f"({st[-1].get('tag')}) の記録 {mrec.get('values')} が違う — どちらを継承するか決めない")
    return vals, {"mode": "from_run", "run": str(src), "file": SM.LIMITER_REF_FILE, "sha256": rec["sha256"],
                  "cfg_fnv": rec["cfg_fnv"], "auto_in_source": rec.get("auto"), "values": vals,
                  "stage_manifest_checked": bool(isinstance(mrec, dict) and mrec.get("values") is not None)}


# --- 系列のリミッタの基準値の自動の継承 (plan tooling-settle-convergence-and-monitor §4.6・§5.1 #6b、2026-10-16) -------------------------
# NS の初期場を別の NS run から作るとき (`prepare_ns(ic_from=…)`; deltastar_loop の run_pass・run_from_table・run_pass0_integral もここに
# 集まる) は、`main_line.limiter_ref` の省略をその run の基準値の継承と読む。B0 の本段 (main_line が mapping) だけ — point の本段には効かせない。
#   - 元の run の `limiter_ref.json` が本段の設定の起動のもの (cfg_fnv が一致) で computed → 4 つ (ρ・P・a・L_ref) を固定値で継承
#   - 記録が無い・cfg_fnv が本段の設定と違う → 元の run の本段の設定の明示 ρ・P・a を継承 (L_ref は明示があれば写し、無ければ新しい run でも
#     forge が格子の境界箱から自動で決める = 元の run と同じ扱い; 既存の固定値の run (case/45 の B0 の系譜・case/71 の p2/p3/p4) は L が自動)。
#     stage_manifest の「値が分かる」(4 つ) より緩いが、新しい run の limiter_ref.json に L の実際の値が残るので区間の鍵は決まる
#   - 記録が壊れている・computed false (元の本段が基準値を使わない)・ρ・P・a の明示が欠ける → ValueError (黙って自動に戻さない;
#     `limiter_ref: auto` か明示で通る)
#   - 元の run の solverConfig.yaml が本段の設定であることは MAIN_STAGE_STARTED (無ければ stage_manifest の最後の段が main) で確かめる。
#     S1・S2 で落ちた run (どちらも本段を示さない) は ValueError、段の記録の無い旧い run は「本段の確認なし」と記録して solverConfig.yaml を使う
_IC_CHAIN_KIND = {"axismach_ns": "ns", "axismach": "euler"}   # prepare_info.json の chain → 種類 (prepare_ns / prepare)


def _ic_source_kind(src: Path) -> dict:
    """初期場の元の run が NS か Euler か (§5.1 #6b)。`<src>/solverConfig.yaml` (最後に回した設定) の turbulence.model・physProp.viscMethod
    (model が none 以外か viscMethod が 0 以外 → NS、model none かつ viscMethod 0 → Euler) と `<src>/prepare_info.json` の chain
    (axismach_ns → NS、axismach → Euler、ほかは判定に使わない) で決める。食い違う・どちらでも決められない・読めない → ValueError。
    戻り {"kind": "ns" | "euler", "by": {判定に使ったファイル: 値と種類}}。"""
    by = {}
    cfgp = src / "solverConfig.yaml"
    if cfgp.is_file():
        try:
            doc = _cfg_load(cfgp.read_text())
        except Exception as e:  # noqa: BLE001 (重複キー・YAML の誤り)
            raise ValueError(f"初期場の元の run の {cfgp} を読めない ({type(e).__name__}: {e})") from e
        tb, pp = doc.get("turbulence"), doc.get("physProp")
        model = tb.get("model") if isinstance(tb, dict) else None
        vm = pp.get("viscMethod") if isinstance(pp, dict) else None
        vm_i = vm if (isinstance(vm, int) and not isinstance(vm, bool)) else None
        model_s = "none" if model is None else str(model).strip().lower()
        kind = "ns" if (model_s != "none" or (vm_i is not None and vm_i != 0)) else ("euler" if vm_i == 0 else None)
        by["solverConfig.yaml"] = {"turbulence.model": model, "viscMethod": vm, "kind": kind}
    pip = src / "prepare_info.json"
    if pip.is_file():
        try:
            chain = json.loads(pip.read_text()).get("chain")
        except Exception as e:  # noqa: BLE001
            raise ValueError(f"初期場の元の run の {pip} を読めない ({type(e).__name__}: {e})") from e
        by["prepare_info.json"] = {"chain": chain, "kind": _IC_CHAIN_KIND.get(chain)}
    kinds = {v["kind"] for v in by.values() if v["kind"] is not None}
    if len(kinds) > 1:
        raise ValueError(f"初期場の元の run {src} が NS か Euler かが食い違う ({by}) — 基準値を継承するか決めない "
                         f"(evaluate.main_line.limiter_ref に {LIMITER_REF_AUTO} か {{ro, p, a}} を明示する)")
    if not kinds:
        raise ValueError(f"初期場の元の run {src} が NS か Euler かを決められない (solverConfig.yaml の turbulence.model・viscMethod、"
                         f"prepare_info.json の chain から; {by or 'どちらのファイルも無い'}) — evaluate.main_line.limiter_ref に "
                         f"{LIMITER_REF_AUTO} か {{ro, p, a}} を明示する")
    return {"kind": kinds.pop(), "by": by}


IC_MAIN_UNCONFIRMED = "本段の確認なし"                      # 元の run の本段の設定を確かめられなかった旧い run の記録 (§5.1 #6b)


def _ic_main_config(src: Path) -> tuple:
    """初期場の元の run の本段の設定 (§5.1 #6b): `<src>/solverConfig.yaml`。本段であることを次の順で確かめ、記録の `main_confirmed` に書く。
    1. `<src>/MAIN_STAGE_STARTED` (run_staged_ns が本段の forge の直前に書く) があれば、その main_config の sha256 (と #6 以後の run なら fnv) が
       solverConfig.yaml と一致すること (違えば ValueError — 最後に回した設定が本段でない・書き換えられた)。
    2. 無ければ `<src>/stage_manifest.json` の最後の段が tag `main` で、その cfg_fnv (あれば) が solverConfig.yaml と一致すること。最後の段が
       main でない (S1・S2 で落ちた run: solverConfig.yaml は前の段の設定のまま) なら ValueError (`limiter_ref: auto` か明示で通る)。
    3. どちらも無い (段の記録の無い旧い run) なら solverConfig.yaml を本段とみなし、`main_confirmed` に IC_MAIN_UNCONFIRMED と書く (from_run と同じ扱い)。
    戻り (設定の文字列, 記録)。"""
    SM = _limiter_tools()
    cfgp = src / "solverConfig.yaml"
    hint = f"(evaluate.main_line.limiter_ref に {LIMITER_REF_AUTO} か {{ro, p, a}} を明示する)"
    if not cfgp.is_file():
        raise ValueError(f"初期場の元の run {src} は NS だが solverConfig.yaml が無い — 本段のリミッタの基準値が分からない {hint}")
    raw = cfgp.read_bytes()
    text = raw.decode()
    sha, fnv = _sha256_bytes(raw), SM.fnv1a64(text)
    rec = {"file": "solverConfig.yaml", "sha256": sha, "fnv": fnv, "main_stage_started": None, "main_confirmed": None}
    ms = src / MAIN_STAGE_STARTED
    if ms.is_file():
        try:
            mc = json.loads(ms.read_text()).get("main_config") or {}
        except Exception as e:  # noqa: BLE001
            raise ValueError(f"初期場の元の run の {ms} を読めない ({type(e).__name__}: {e})") from e
        if mc.get("sha256") != sha:
            raise ValueError(f"初期場の元の run {src} の solverConfig.yaml (sha256 {sha[:16]}…) が {MAIN_STAGE_STARTED} の本段の設定 "
                             f"(sha256 {str(mc.get('sha256'))[:16]}…) と違う — 本段の設定が分からないので基準値を継承しない")
        if mc.get("fnv") is not None and mc["fnv"] != fnv:
            raise ValueError(f"初期場の元の run {src} の solverConfig.yaml (fnv {fnv}) が {MAIN_STAGE_STARTED} の本段の設定 (fnv {mc['fnv']}) と違う")
        rec["main_stage_started"] = {"sha256_checked": True, "fnv_checked": mc.get("fnv") is not None}
        rec["main_confirmed"] = f"{MAIN_STAGE_STARTED} の main_config"
        return text, rec
    man = src / "stage_manifest.json"
    stages = []
    if man.is_file():
        try:
            stages = json.loads(man.read_text()).get("stages") or []
        except Exception as e:  # noqa: BLE001
            raise ValueError(f"初期場の元の run の {man} を読めない ({type(e).__name__}: {e})") from e
    if stages:
        last = stages[-1] if isinstance(stages[-1], dict) else {}
        if last.get("tag") != "main":
            raise ValueError(f"初期場の元の run {src} に {MAIN_STAGE_STARTED} が無く、stage_manifest の最後の段が {last.get('tag')!r} (main でない) — "
                             f"本段の前の段で終わった run で、solverConfig.yaml は本段の設定でない。基準値を継承しない {hint}")
        if last.get("cfg_fnv") is not None and last["cfg_fnv"] != fnv:
            raise ValueError(f"初期場の元の run {src} の stage_manifest の最後の段 main (cfg_fnv {last['cfg_fnv']}) が solverConfig.yaml (fnv {fnv}) "
                             f"と違う — 本段の設定が分からないので基準値を継承しない {hint}")
        rec["main_confirmed"] = "stage_manifest の最後の段 (tag main)"
        return text, rec
    rec["main_confirmed"] = (f"{IC_MAIN_UNCONFIRMED} ({MAIN_STAGE_STARTED} も stage_manifest の段も無い旧い run; solverConfig.yaml を本段と"
                             "みなした)")
    return text, rec


def _ic_manifest_main_check(src: Path, vals: dict, fnv: str) -> dict:
    """継承する値を元の run の stage_manifest.json の tag `main` の段 (最後のもの) の記録と照らす (§5.1 #6b; 最後の段でなく本段)。
    main の段の cfg_fnv が本段の設定と違う・継承する値と記録の値が違う → ValueError。照らせなければ (ファイル・main の段・値が無い) 理由を返す。"""
    man = src / "stage_manifest.json"
    if not man.is_file():
        return {"checked": False, "reason": "stage_manifest.json が無い"}
    try:
        st = json.loads(man.read_text()).get("stages") or []
    except Exception as e:  # noqa: BLE001
        raise ValueError(f"初期場の元の run の {man} を読めない ({type(e).__name__}: {e})") from e
    mains = [s for s in st if isinstance(s, dict) and s.get("tag") == "main"]
    if not mains:
        return {"checked": False, "reason": "tag main の段が無い (本段の終わる前・記録の前の run)"}
    m = mains[-1]
    if m.get("cfg_fnv") is not None and m["cfg_fnv"] != fnv:
        raise ValueError(f"初期場の元の run {src} の stage_manifest の main の段 (cfg_fnv {m['cfg_fnv']}) が本段の設定 (fnv {fnv}) と違う "
                         "— どちらが本段か決めない")
    lr = m.get("limiter_ref")
    mv = lr.get("values") if isinstance(lr, dict) else None
    if not isinstance(mv, dict):
        return {"checked": False, "reason": f"main の段に基準値の値が無い (status {lr.get('status') if isinstance(lr, dict) else None})"}
    bad = {k: (mv.get(k), v) for k, v in vals.items() if mv.get(k) != v}
    if bad:
        raise ValueError(f"初期場の元の run {src} から継承する基準値と stage_manifest の main の段の記録が違う ({bad}; (記録, 継承)) "
                         "— どちらを継承するか決めない")
    return {"checked": True, "tag": "main", "order": m.get("order"), "status": lr.get("status")}


def limiter_ref_inherit_from_ic(ic_from) -> tuple:
    """初期場の元の run (`prepare_ns` の ic_from) のリミッタの基準値を継承する (plan tooling-settle-convergence-and-monitor §4.6・§5.1 #6b;
    `main_line.limiter_ref` の省略のとき `_main_line_spec` が呼ぶ)。戻り ({space のキー: float}, 出所の記録)。空の dict = 自動。
    - 元の run が Euler (`_ic_source_kind`) → 自動 (新しい物理問題)。
    - 元の run が NS → 本段の設定 (`_ic_main_config`) の起動の `limiter_ref.json` (cfg_fnv 一致・computed) の 4 つ、無ければ本段の設定の明示
      ρ・P・a (+ 明示の L)。規則は上の区切りのコメント。基準値が分からないときは ValueError。"""
    SM = _limiter_tools()
    src = Path(str(ic_from)).expanduser()
    kind = _ic_source_kind(src)
    if kind["kind"] == "euler":
        return {}, {"mode": "auto", "reason": "初期場の元が Euler の run (新しい物理問題)", "ic_from": str(src), "ic_kind": kind,
                    "note": "forge が本段の開始の場で決めて run の limiter_ref.json に書く (stage_manifest の本段の段に写る)"}
    text, main_rec = _ic_main_config(src)
    base = {"mode": "inherit", "run": str(src), "ic_kind": kind, "main_confirmed": main_rec["main_confirmed"], "main_config": main_rec}
    try:
        doc = _cfg_load(text)
    except Exception as e:  # noqa: BLE001
        raise ValueError(f"初期場の元の run の solverConfig.yaml を読めない ({type(e).__name__}: {e})") from e
    sp = doc.get("space") if isinstance(doc.get("space"), dict) else {}

    def _posf(v):
        return float(v) if (not isinstance(v, bool) and isinstance(v, (int, float)) and np.isfinite(float(v)) and float(v) > 0.0) else None
    expl = {k: _posf(sp[k]) for k in SM.LIMITER_REF_NAMES if k in sp}
    bad_expl = sorted(k for k, v in expl.items() if v is None)
    if bad_expl:
        raise ValueError(f"初期場の元の run {src} の本段の設定の基準値 {bad_expl} が正の有限の数でない — 継承しない")
    recp = src / SM.LIMITER_REF_FILE
    skipped = f"{SM.LIMITER_REF_FILE} が無い"
    if recp.is_file():
        try:
            rec = SM.read_limiter_ref(src)              # 形式・値の検査 (起動の照合は下で; 不一致は固定値へ回す)
        except SM.LimiterRecordError as e:
            raise ValueError(f"初期場の元の run の記録が壊れている — 黙って固定値に落とさない: {e}") from e
        if rec["cfg_fnv"] != main_rec["fnv"]:
            skipped = (f"{SM.LIMITER_REF_FILE} は本段の設定の起動のものでない (cfg_fnv {rec['cfg_fnv']}、本段の設定 {main_rec['fnv']}; "
                       f"sha256 {rec['sha256']})")
        elif not rec["computed"]:
            raise ValueError(f"初期場の元の run {src} の本段は基準値を使わない起動 (limiter_ref.json computed false、limiterScaled "
                             f"{rec.get('limiterScaled')}) — 継承できない (evaluate.main_line.limiter_ref に {LIMITER_REF_AUTO} か "
                             "{ro, p, a} を明示する)")
        else:
            vals = {k: float(rec[k]) for k in SM.LIMITER_REF_NAMES}
            clash = sorted(k for k, v in expl.items() if v != vals[k] or rec["auto"][k])
            if clash:
                raise ValueError(f"初期場の元の run {src} の本段の設定の明示 {clash} と limiter_ref.json の値・auto が食い違う — 継承しない")
            man = _ic_manifest_main_check(src, vals, main_rec["fnv"])
            return vals, {**base, "basis": SM.LIMITER_REF_FILE, "file": SM.LIMITER_REF_FILE, "sha256": rec["sha256"],
                          "cfg_fnv": rec["cfg_fnv"], "auto_in_source": dict(rec["auto"]), "values": vals,
                          "L_ref": "fixed (記録の値)", "stage_manifest": man}
    if sp.get("limiterScaled") == 0 and not sp.get("limiterDiag"):
        raise ValueError(f"初期場の元の run {src} の本段は limiterScaled 0 (基準値を使わない) — 継承できない "
                         f"(evaluate.main_line.limiter_ref に {LIMITER_REF_AUTO} か {{ro, p, a}} を明示する)")
    missing = [sk for _, sk in _LIMITER_REF_KEYS if sk not in expl]
    if missing:
        raise ValueError(f"初期場の元の run {src} は NS だが基準値が分からない ({skipped}、本段の設定に {missing} の明示が無い) — 黙って自動に"
                         f"戻さない (evaluate.main_line.limiter_ref に {LIMITER_REF_AUTO} か {{ro, p, a}} を明示する)")
    vals = {sk: expl[sk] for _, sk in _LIMITER_REF_KEYS + _LIMITER_REF_OPT_KEYS if sk in expl}
    man = _ic_manifest_main_check(src, vals, main_rec["fnv"])
    return vals, {**base, "basis": "solverConfig.yaml", "file": "solverConfig.yaml", "sha256": main_rec["sha256"], "values": vals,
                  "record_skipped": skipped,
                  "L_ref": ("fixed (本段の設定の明示)" if "limiterRefLength" in vals
                            else "auto (元の run と同じく forge が新しい run の格子の境界箱から決める; 値は新しい run の limiter_ref.json)"),
                  "stage_manifest": man}


def _main_line_spec(p: Problem, cfl_main=None, implicit_relax=None, source_out: dict | None = None, ic_from=None) -> dict | None:
    """`evaluate.main_line` を読んで検査し、解決した値を返す (キー無し・null・`point` は None = point の本段、今と同じ設定・同じ段)。
    {line: dir, cap, itj, cfl, relax, limiter_ref, extra_fields: [...]}。本段の cfl・implicitRelax は main_line が
    決めるので、引数・YAML の `cfl_main` と同時に指定したら ValueError (どちらが効くかを曖昧にしない)。`implicit_relax` は
    前の段 (S1・S2) にも効く値なので、main_line.relax と同じ値のときだけ併記を許す (2026-10-13 親の決定: 生産の系譜 run_0182 は
    `prepare_ns(implicit_relax=0.7)` で S1・S2 も 0.7。違う値は ValueError)。
    line は dir (lineImplicit 1 + lineDtDirectional 1) だけ。未知のキー・欠けたキー・範囲外の値も ValueError。
    **limiter_ref** (plan tooling-settle-convergence-and-monitor §4.6・§5.1 #6b、2026-10-16):
      - 省略 (null): `ic_from` (prepare_ns の初期場の元の run) が NS の run なら**その基準値を継承** (`limiter_ref_inherit_from_ic`;
        記録の 4 つ、無ければ本段の設定の明示 ρ・P・a (+ L)、分からなければ ValueError)。`ic_from` が無い・Euler の run なら **自動**
        (本段の設定に基準値を書かず、forge が本段の開始の場で決めて run の `limiter_ref.json` に書く; stage_manifest の本段の段にも写る;
        戻りの limiter_ref は {})。ic_from を渡さない呼び出し (case/45 の台本の登録の照合など) は今までどおり自動。
      - `auto` (文字列) = 初期場の元を見ずに自動 (継承しない明示)。
      - `{ro, p, a}` (+ 任意の `L` = limiterRefLength) = 明示の固定値 (今までどおり; 継承より優先)。
      - `{from_run: <run>}` = 別 run の記録を継承 (`limiter_ref_from_run`; ρ・P・a・L_ref の 4 つを固定値で書く; 初期場の元より優先)。
    基準値の由来 (mode auto / inherit / explicit / from_run と記録の sha256 など) は `source_out` (dict を渡したとき) の `limiter_ref` に書く
    (戻りの dict のキーは今までと同じ — case/45 の台本が登録の dict と比べる)。"""
    ml = p.evaluate.get("main_line")
    if ml is None or _main_line_mode(p)["value"] == "point":
        return None
    if not isinstance(ml, dict):
        raise ValueError(f"evaluate.main_line は mapping であること ({ml!r})")
    clash = [w for w, on in (("引数 cfl_main", cfl_main is not None), ("evaluate.cfl_main", "cfl_main" in p.evaluate)) if on]
    if clash:
        raise ValueError(f"evaluate.main_line と {', '.join(clash)} は同時に指定できない "
                         "(本段の cfl は main_line の cfl で決める。前の段は S1 0.5・S2 1.0 のまま)")
    unknown = sorted(set(ml) - set(_MAIN_LINE_KEYS), key=str)
    missing = [k for k in _MAIN_LINE_REQUIRED if k not in ml]
    if unknown or missing:
        raise ValueError(f"evaluate.main_line のキーが不正 (未知 {unknown}、欠け {missing}; 使えるのは {list(_MAIN_LINE_KEYS)})")
    if ml["line"] != "dir":
        raise ValueError(f"evaluate.main_line.line は dir だけ (lineImplicit 1 + lineDtDirectional 1; 指定 {ml['line']!r})")

    def _pos(name, v, hi=None):
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not np.isfinite(float(v)) or float(v) <= 0.0 \
                or (hi is not None and float(v) > hi):
            raise ValueError(f"evaluate.main_line.{name} は有限の正の数{'' if hi is None else f' (≤ {hi})'}であること ({v!r})")
        return float(v)
    itj = ml["itj"]
    if isinstance(itj, bool) or not isinstance(itj, int) or not 0 <= itj <= 7 or ((itj & 4) and not (itj & 1)):
        # ソルバの検査 (solverConfig.cpp: 0..7、ビット 4 はビット 1 が要る) と同じ
        raise ValueError(f"evaluate.main_line.itj は 0..7 の整数 (ビット 4 はビット 1 と一緒に) であること ({itj!r})")
    lr = ml.get("limiter_ref")
    need = {k for k, _ in _LIMITER_REF_KEYS}
    opt = {k for k, _ in _LIMITER_REF_OPT_KEYS}
    auto_note = "forge が本段の開始の場で決めて run の limiter_ref.json に書く (stage_manifest の本段の段に写る)"
    if lr is None and ic_from is not None:
        refs, lr_src = limiter_ref_inherit_from_ic(ic_from)       # NS の元なら継承、Euler の元なら自動 (分からなければ ValueError)
    elif lr is None:
        refs = {}
        lr_src = {"mode": "auto", "reason": "初期場の元の run が無い (新しい問題)", "note": auto_note}
    elif lr == LIMITER_REF_AUTO:
        refs = {}
        lr_src = {"mode": "auto", "reason": f"evaluate.main_line.limiter_ref: {LIMITER_REF_AUTO} の明示 (初期場の元を見ない)",
                  "ic_from": (None if ic_from is None else str(ic_from)), "note": auto_note}
    elif isinstance(lr, dict) and set(lr) == {"from_run"}:
        refs, lr_src = limiter_ref_from_run(lr["from_run"], base_dir=(Path(p.path).resolve().parent if getattr(p, "path", None) else None))
    elif isinstance(lr, dict) and need <= set(lr) <= need | opt:
        refs = {sk: _pos(f"limiter_ref.{k}", lr[k]) for k, sk in _LIMITER_REF_KEYS + _LIMITER_REF_OPT_KEYS if k in lr}
        lr_src = {"mode": "explicit", "values": dict(refs)}
    else:
        raise ValueError(f"evaluate.main_line.limiter_ref は 省略 (初期場の元の NS run から継承、無ければ自動)・{LIMITER_REF_AUTO}・"
                         f"{{ro, p, a}} (+ 任意の L)・{{from_run: <run>}} のどれかであること ({lr!r})")
    ex = ml.get("extra_fields")
    if ex is not None and (not isinstance(ex, list) or not ex or len(set(map(str, ex))) != len(ex)
                           or any(not isinstance(s, str) or not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", s) for s in ex)):
        raise ValueError(f"evaluate.main_line.extra_fields は重複の無い量の名前の空でないリストであること ({ex!r})")
    relax = _pos("relax", ml["relax"], hi=1.0)
    if source_out is not None:
        source_out["limiter_ref"] = lr_src
    ir = implicit_relax if implicit_relax is not None else p.evaluate.get("implicit_relax")
    if ir is not None and float(ir) != relax:
        raise ValueError(f"implicit_relax ({ir!r}、前の段にも効く) は evaluate.main_line.relax ({relax!r}) と同じ値に限る "
                         "(違う値を前の段と本段に分けて書く形は決めていない)")
    return {"line": "dir", "cap": _pos("cap", ml["cap"]), "itj": int(itj), "cfl": _pos("cfl", ml["cfl"]),
            "relax": relax, "prestage_relax": (None if ir is None else float(ir)), "limiter_ref": refs,
            "extra_fields": (list(ex) if ex else None)}


def _apply_main_line(cfg: str, ml: dict) -> str:
    """前の段の元の設定 cfg に B0 (`_main_line_spec` の解決済みの値) を当てた本段の設定を返す。キー・書き方・前提の検査は
    case/45 の `cold_cfl.py prep` と同じ (cfl・cfl_pseudo の値の置換 → output → space の先頭にリミッタの基準値 →
    deltaT の先頭に lineImplicit・lineDtDirectional → implicitThermalJacobian → lineDtDirectionalCap の順)。implicitRelax は
    prepare_ns の `implicit_relax` と同じく `blockDPLUR: 1,` の直後に挿入する。読み直して、元との差が B0 のパスだけで値が
    指定どおりであることを確かめる (違えば ValueError)。"""
    doc = _cfg_load(cfg)
    tm = doc.get("time")
    dt = tm.get("deltaT") if isinstance(tm, dict) else None
    sp = doc.get("space")
    if not isinstance(dt, dict) or not isinstance(sp, dict):
        raise ValueError("solverConfig に time.deltaT・space の mapping が無い — main_line を書けない")
    if int(dt.get("blockDPLUR", 0)) != 1 or int(tm.get("timeIntegration", 0)) != 11 or int(dt.get("lowMachPrecond", 0)) >= 2:
        raise ValueError("lineImplicit は timeIntegration 11 + blockDPLUR 1 + lowMachPrecond < 2 専用 — main_line を書けない")
    pre_relax = _cfg_path_get(doc, ("time", "deltaT", "implicitRelax"))
    if pre_relax is not _CFG_MISSING and float(pre_relax) != float(ml["relax"]):
        raise ValueError(f"前の段の元の設定の implicitRelax ({pre_relax!r}) が main_line.relax ({ml['relax']!r}) と違う")
    have = [".".join(pth) for pth in sorted(_MAIN_LINE_ONLY_PATHS)
            if _cfg_path_get(doc, pth) is not _CFG_MISSING and (pth != ("output",) or ml["extra_fields"])]
    if have:
        raise ValueError(f"前の段の元の設定に既に {have} がある — main_line を重ねて書けない")
    if cfg.count("deltaT: {") != 1 or cfg.count("space: {") != 1 or cfg.count("blockDPLUR: 1,") != 1:
        raise ValueError("deltaT・space が 1 つのフロー形式でないか `blockDPLUR: 1,` が 1 か所でない — main_line の足し方を決めていない")
    text = _cfg_set(cfg, {"cfl": ml["cfl"], "cfl_pseudo": ml["cfl"]})
    if ml["extra_fields"]:
        text = text.rstrip("\n") + "\noutput: {level: 1, extraFields: [" + ", ".join(ml["extra_fields"]) + "]}\n"
    if ml["limiter_ref"]:                  # 空 = 自動 (forge が本段の開始の場で決めて limiter_ref.json に書く) なら書かない
        text = text.replace("space: {", "space: {" + ", ".join(f"{k}: {float(v)!r}" for k, v in ml["limiter_ref"].items()) + ", ", 1)
    text = text.replace("deltaT: {", "deltaT: {lineImplicit: 1, lineDtDirectional: 1, ", 1)
    text = text.replace("deltaT: {", f"deltaT: {{implicitThermalJacobian: {int(ml['itj'])}, ", 1)
    text = text.replace("deltaT: {", f"deltaT: {{lineDtDirectionalCap: {float(ml['cap'])!r}, ", 1)
    if pre_relax is _CFG_MISSING:          # 前の段に implicitRelax が無いときだけ本段に挿入 (あれば同じ値を引き継ぐ)
        text = text.replace("blockDPLUR: 1,", f"blockDPLUR: 1, implicitRelax: {float(ml['relax'])},", 1)
    want = {("time", "deltaT", "lineImplicit"): 1, ("time", "deltaT", "lineDtDirectional"): 1,
            ("time", "deltaT", "lineDtDirectionalCap"): float(ml["cap"]),
            ("time", "deltaT", "implicitThermalJacobian"): int(ml["itj"]),
            ("time", "deltaT", "cfl"): float(ml["cfl"]), ("time", "deltaT", "cfl_pseudo"): float(ml["cfl"]),
            ("time", "deltaT", "implicitRelax"): float(ml["relax"]),
            **{("space", k): float(v) for k, v in ml["limiter_ref"].items()}}
    if ml["extra_fields"]:
        want[("output",)] = {"level": 1, "extraFields": list(ml["extra_fields"])}
    after = _cfg_load(text)
    stray = sorted(set(_cfg_diff_paths(doc, after)) - set(want), key=str)
    bad = {".".join(pth): _cfg_path_get(after, pth) for pth, v in want.items()
           if _cfg_path_get(after, pth) is _CFG_MISSING or _cfg_path_get(after, pth) != v}
    if stray or bad:
        raise ValueError(f"main_line を当てた本段の設定を読み直したら想定と違う (B0 以外の差 {stray}、値の違い {bad})")
    return text


def _ns_configs(p: Problem, run_dir, n: int, out_int: int, cfl_main=None, implicit_relax=None, main_line: dict | None = None):
    """NS の (前の段の元の設定, 本段の設定) を返す。`main_line` (`_main_line_spec` の戻り値) が None なら 2 つは同じ文字列で、
    従来の prepare_ns の組み立てと同じ (`_config_sst_node` → `_apply_gas_to_config` → implicit_relax の挿入)。
    `_apply_gas_to_config` は run_dir に species_meta.yaml などを書く (従来どおり 1 回だけ呼ぶ)。"""
    cfl_main = float(cfl_main if cfl_main is not None else p.evaluate.get("cfl_main", 1.0))
    implicit_relax = implicit_relax if implicit_relax is not None else p.evaluate.get("implicit_relax")
    cfg_ns = _apply_gas_to_config(_config_sst_node(p, n, out_int, cfl_main), p, run_dir, viscous=True)
    if implicit_relax is not None:
        # 陰解法の緩和 (cfl 6 + implicitRelax 0.7 が生産推奨: case/45 run_0018)。deltaT ブロックに挿入
        cfg_ns = cfg_ns.replace("blockDPLUR: 1,", f"blockDPLUR: 1, implicitRelax: {float(implicit_relax)},", 1)
    if main_line is None:
        return cfg_ns, cfg_ns
    return cfg_ns, _apply_main_line(cfg_ns, main_line)


# --- 本段の監視用の出力と判定の基準 (`evaluate.ns_output`・`evaluate.settle`; plan tooling-settle-convergence-and-monitor §4.5・§5.1 #6) -----
# 要求カードの雛形 (axismach_final_v1・axismach_explore_v1) が既定で書く。キー無しの問題は今と同じ設定・同じファイル。
# 出力のキーは**本段の設定だけ**に書き (B0 と同じ)、前の段 (S1・S2) は PRESTAGE_CONFIG (出力のキーの無い元の設定) から作る。
NS_OUTPUT_KEYS = ("outStepInterval", "fieldOutStepInterval", "boundary_h5", "boundaryMassflux", "interfaceDiag")
NS_BOUNDARY_NAMES = ("inlet", "outlet", "wall", "axis")       # runner_wt._bcond の境界の名前
SETTLE_NODES_FILE = "settle_nodes.txt"                        # forge の output.nodeProbeFile (settle_presets.probe_node_ids)
# 本段だけに書く出力のパス (前の段の設定に残っていてはいけない; 本段と前の段の元の設定の差として許す)
_NS_OUTPUT_PATHS = ({("time", "fieldOutStepInterval"), ("output",)}
                    | {("output", k) for k in ("boundaryMassflux", "nodeProbeFile", "interfaceDiag")})


def _settle_presets():
    from . import settle_presets
    return settle_presets


def _ns_output_spec(p: Problem) -> dict | None:
    """`evaluate.ns_output` (NS だけが読む; Euler の `prepare` は読まない) を検査して返す (キー無しは None = 今と同じ)。
    {outStepInterval: 2000 (境界の h5・節点・流量・残差の判定の格子; `_ns_steps` が evaluate.outStepInterval より優先),
     fieldOutStepInterval: 20000 (全場 res_<step>.h5 だけ; forge の time.fieldOutStepInterval、本段だけ),
     boundary_h5: [inlet, outlet, wall, axis] (outputHDFflg を 1 にする境界; 書いていない境界は今のまま),
     boundaryMassflux: 1 (本段の output.boundaryMassflux)}。知らないキー・不正な値は ValueError。"""
    no = p.evaluate.get("ns_output")
    if no is None:
        return None
    if not isinstance(no, dict):
        raise ValueError(f"evaluate.ns_output は mapping であること ({no!r})")
    bad = sorted(set(no) - set(NS_OUTPUT_KEYS), key=str)
    if bad:
        raise ValueError(f"evaluate.ns_output の知らないキー {bad} (使えるのは {list(NS_OUTPUT_KEYS)})")

    def _pint(k):
        v = no.get(k)
        if v is None:
            return None
        if isinstance(v, bool) or not isinstance(v, int) or v < 1:
            raise ValueError(f"evaluate.ns_output.{k} は 1 以上の整数であること ({v!r})")
        return int(v)
    oi, fi = _pint("outStepInterval"), _pint("fieldOutStepInterval")
    if fi is not None and oi is not None and fi % oi:
        raise ValueError(f"evaluate.ns_output.fieldOutStepInterval {fi} が outStepInterval {oi} の倍数でない (全場が出力の格子に乗らない)")
    bh = no.get("boundary_h5", [])
    if not isinstance(bh, list) or len(set(map(str, bh))) != len(bh) or any(b not in NS_BOUNDARY_NAMES for b in bh):
        raise ValueError(f"evaluate.ns_output.boundary_h5 は {list(NS_BOUNDARY_NAMES)} の重複の無いリストであること ({bh!r})")
    bm = no.get("boundaryMassflux", 0)
    if bm not in (0, 1) or isinstance(bm, float):
        raise ValueError(f"evaluate.ns_output.boundaryMassflux は 0 か 1 ({bm!r})")
    idg = no.get("interfaceDiag", 0)
    if idg not in (0, 1) or isinstance(idg, (bool, float)):
        raise ValueError(f"evaluate.ns_output.interfaceDiag は 0 か 1 ({idg!r})")
    # interfaceDiag: 壁の熱流の判定に使わなくても成果物の壁の熱流 (iface_Qf_eff) のために出す (2026-10-16、wt_nozzle_v3 は壁の熱流を判定から外した)
    return {"outStepInterval": oi, "fieldOutStepInterval": fi, "boundary_h5": list(bh), "boundaryMassflux": int(bm), "interfaceDiag": int(idg)}


def _settle_spec(p: Problem) -> tuple:
    """`evaluate.settle` を読む: (問題 YAML の evaluate.settle の dict か None, prepare_info の `settle` の記録)。
    - キー無し: (None, 未指定)。今と同じ (settle_criteria.yaml を書かない)。
    - `false`: (None, 外した・理由)。**明示のときだけ外す** (前提を満たせない問題は settle_presets.resolve の例外で止まる)。
    - mapping: キーとプリセットの名前をここで検査 (run dir を作る前)。中身の解決は prepare_ns が run の設定で行う。
    それ以外 (true など) は ValueError。"""
    es = p.evaluate.get("settle", None)
    if "settle" not in p.evaluate:
        return None, {"status": "unset", "reason": "evaluate.settle が無い (settle_criteria.yaml を書かない; 今までどおり)"}
    if es is False:
        return None, {"status": "off", "reason": "evaluate.settle: false (問題 YAML で明示して外した)"}
    if not isinstance(es, dict):
        raise ValueError(f"evaluate.settle は mapping (preset ほか) か false であること ({es!r})")
    SP = _settle_presets()
    unknown = [k for k in es if k not in SP.SETTLE_KEYS]
    if unknown or es.get("preset") not in SP.PRESETS:
        raise ValueError(f"evaluate.settle が不正 (知らないキー {unknown}、preset {es.get('preset')!r}; キー {list(SP.SETTLE_KEYS)}、"
                         f"プリセット {list(SP.PRESETS)})")
    return dict(es), {"status": "on", "preset": es["preset"]}


def _bcond_set_output(txt: str, names) -> str:
    """bcondConfig (runner_wt._bcond の 1 行 1 境界の flow 形式) の names の境界の outputHDFflg を 1 にする。読み直して、差が
    その境界の outputHDFflg だけであることを確かめる (違えば ValueError)。"""
    import yaml
    before = yaml.safe_load(txt)
    out = txt
    for nm in names:
        pat = re.compile(rf"^({re.escape(nm)}:\s*\{{.*?\boutputHDFflg:\s*)(\d+)", re.M)
        if len(pat.findall(out)) != 1:
            raise ValueError(f"bcondConfig に境界 {nm!r} の outputHDFflg の行がちょうど 1 つでない")
        out = pat.sub(lambda m: m.group(1) + "1", out, count=1)
    after = yaml.safe_load(out)
    want = {nm: {**before[nm], "outputHDFflg": 1} for nm in names}
    if after != {**before, **want}:
        raise ValueError("bcondConfig の outputHDFflg を書き換えたら想定と違う")
    return out


def _settle_needs(crit: dict) -> dict:
    """解決した settle_criteria の量から、本段の出力に要るもの {nodeProbeFile, boundaryMassflux, interfaceDiag}。"""
    kinds = {q["kind"] for q in crit["quantities"] if "series_csv" not in q}
    return {"nodeProbeFile": bool(kinds & {"exit_core_mach", "line_mach", "exit_flow_angle"}),
            "boundaryMassflux": "mass_balance" in kinds, "interfaceDiag": "wall_heat_nodes" in kinds}


def _apply_ns_outputs(cfg: str, no: dict | None, needs: dict | None) -> str:
    """本段の設定に監視用の出力を足す (plan tooling-settle-convergence-and-monitor §4.1・§4.5)。no = `_ns_output_spec`、needs =
    `_settle_needs` (settle が無ければ None)。time.fieldOutStepInterval を outStepInterval の行の直後に、output のキー
    (boundaryMassflux・nodeProbeFile・interfaceDiag) を output の mapping に (無ければ最後に 1 行) 書く。outStepInterval は
    `_ns_steps` が既に設定に入れている。読み直して、元との差が _NS_OUTPUT_PATHS だけで値が要求どおりであることを確かめる。
    足すものが無ければ cfg をそのまま返す (キー無しの問題は今と同じ文字列)。"""
    no = no or {}
    needs = needs or {}
    outk = {}
    if no.get("boundaryMassflux") or needs.get("boundaryMassflux"):
        outk["boundaryMassflux"] = 1
    if needs.get("nodeProbeFile"):
        outk["nodeProbeFile"] = SETTLE_NODES_FILE
    if needs.get("interfaceDiag") or no.get("interfaceDiag"):
        outk["interfaceDiag"] = 1
    fi = no.get("fieldOutStepInterval")
    if not outk and fi is None:
        return cfg
    doc = _cfg_load(cfg)
    have = [".".join(pth) for pth in sorted(_NS_OUTPUT_PATHS - {("output",)}) if _cfg_path_get(doc, pth) is not _CFG_MISSING]
    if have:
        raise ValueError(f"前の段の元の設定に既に {have} がある — 本段の出力を重ねて書けない")
    text = cfg
    if fi is not None:
        pat = re.compile(r"^([ \t]*)outStepInterval:[ \t]*\S+[ \t]*$", re.M)
        if len(pat.findall(text)) != 1:
            raise ValueError("time.outStepInterval の行がちょうど 1 つの block 形式でない — fieldOutStepInterval の足し方を決めていない")
        text = pat.sub(lambda m: m.group(0) + "\n" + m.group(1) + f"fieldOutStepInterval: {int(fi)}", text, count=1)
    if outk:
        flow = ", ".join(f"{k}: {v}" for k, v in outk.items())
        if _cfg_path_get(doc, ("output",)) is _CFG_MISSING:
            text = text.rstrip("\n") + "\noutput: {" + flow + "}\n"
        elif len(re.findall(r"^output: \{", text, re.M)) == 1:
            text = re.sub(r"^output: \{", "output: {" + flow + ", ", text, count=1, flags=re.M)
        else:
            raise ValueError("output が 1 行の flow 形式でない — 本段の出力の足し方を決めていない")
    want = {("output", k): v for k, v in outk.items()}
    if fi is not None:
        want[("time", "fieldOutStepInterval")] = int(fi)
    after = _cfg_load(text)
    stray = sorted(set(_cfg_diff_paths(doc, after)) - _NS_OUTPUT_PATHS, key=str)
    bad = {".".join(pth): _cfg_path_get(after, pth) for pth, v in want.items() if _cfg_path_get(after, pth) != v}
    if stray or bad:
        raise ValueError(f"本段の出力を足した設定を読み直したら想定と違う (出力以外の差 {stray}、値の違い {bad})")
    return text


def write_settle_files(run_dir, settle_es: dict, info: dict, cfg_main: str, bc_text: str) -> dict:
    """run に判定の基準 `settle_criteria.yaml` と節点の表 `settle_nodes.txt` を書く (plan tooling-settle-convergence-and-monitor §4.2・
    §5.1 #6; prepare_ns が本段の設定と bcondConfig を書いた後に呼ぶ)。info は prepare_info の中身 (Md・scale_m・mesh の ni・nj・x_E・
    physical_wall)、cfg_main は本段の設定 (出力のキー入り)。節点の座標は本段の設定の mesh.meshFileName (変換した格子) から。
    `grid.segment` は本段の開始の通知 (`settle_presets.main_stage_segment`)。基準の前提 (本段の出力の設定) が 1 つでも満たされて
    いなければ ValueError (書いた本段の設定の誤り)。戻り: prepare_info の `settle` に足す記録。"""
    import yaml
    SP = _settle_presets()
    run_dir = Path(run_dir)
    main_doc = _cfg_load(cfg_main)
    bc_doc = yaml.safe_load(bc_text)
    seg = SP.main_stage_segment()
    crit = SP.resolve(settle_es, info, main_doc, bc_doc, segment=seg)
    nodes_rec = None
    if _settle_needs(crit)["nodeProbeFile"]:
        nf = str((main_doc.get("output") or {}).get("nodeProbeFile", ""))
        if nf != SETTLE_NODES_FILE:
            raise ValueError(f"本段の設定の output.nodeProbeFile ({nf!r}) が {SETTLE_NODES_FILE} でない")
        ni_m, nj_m = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
        X, R = SP.mesh_xr(run_dir / str(main_doc["mesh"]["meshFileName"]), ni_m, nj_m, float(info["scale_m"]))
        ids = SP.probe_node_ids(ni_m, nj_m, X, R, crit)
        if len(ids) == 0:
            raise ValueError("判定の基準の節点の量の節点が 0 個 (forge は空の節点の表を拒む)")
        SP.write_node_probe_file(run_dir / SETTLE_NODES_FILE, ids,
                                 header=(f"settle_presets.probe_node_ids (preset {settle_es['preset']}; 構造格子 i·nj + j、ni {ni_m}・nj {nj_m})\n"
                                         "runner_axismach.prepare_ns が書いた — 手で直さない"))
        nodes_rec = {"file": SETTLE_NODES_FILE, "n": int(len(ids)), "sha256": _sha256_file(run_dir / SETTLE_NODES_FILE)}
    w = SP.write_settle_criteria(run_dir, settle_es, info, main_doc, bc_doc, segment=seg)
    unmet = [c for c in w["criteria"]["preconditions"] if not c["satisfied"]]
    if unmet:
        raise ValueError(f"判定の基準の前提を本段の設定が満たさない (書いた出力の誤り): {unmet}")
    return {"file": Path(w["path"]).name, "sha256": w["sha256"], "nodes": nodes_rec,
            "quantities": [q["name"] for q in w["criteria"]["quantities"]], "excluded": w["criteria"]["excluded"],
            "evaluate_settle": settle_es}


def euler_cfl(p: Problem):
    """Euler の本段の cfl (`evaluate.cfl_euler` > `evaluate.cfl_main` > None = 設定の既定 4.0)。要求カードの雛形は NS の本段を B0
    (`main_line`、cfl_main と併記できない) にしたので、Euler の cfl を別に持つ (codex diagnose 2026-10-15 m55-three-coef M4: cfl_main を
    消すだけでは Euler が cfl 4 に変わる)。Euler の台本 (run_staged の呼び出し側) はこれを使う。"""
    v = p.evaluate.get("cfl_euler", p.evaluate.get("cfl_main"))
    if v is None:
        return None
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not np.isfinite(float(v)) or float(v) <= 0.0:
        raise ValueError(f"evaluate.cfl_euler (または cfl_main) は有限の正の数であること ({v!r})")
    return float(v)


def _ns_steps(p: Problem, nsteps=None) -> tuple:
    """NS の本段の step 数・出力の間隔とその出所 (n, out_int, n の出所)。prepare_ns と停止規則の設定 (`_stop_rule_settings`) が共有する
    (M6 の版が prepare_ns の中の式を関数に出した; plan tooling-nozzle-isothermal-wall-chain §5.1 #36a)。
    step 数は 2026-10-14 のユーザ決定「1 ok」(forge-f5 fb07daf8、plan tooling-nozzle-cooled-wall-ns-settling §4.5) のとおり、引数 nsteps
    (`--steps`) > `evaluate.ns_steps` > `NS_STEPS_DEFAULT` (60,000)。`evaluate.nStepOuter` は Euler (`prepare`) だけが読む
    (雛形の 12,000 が NS にも効いて本段が短すぎた)。停止規則を使うときの本段の長さは上限 (`stop_rule.resolve_settings` の cap) で、
    この n は使わない・使えないときの長さ (`fallback_steps`)。
    出力の間隔は `evaluate.ns_output.outStepInterval` (NS だけ; plan tooling-settle-convergence-and-monitor §4.5) > `evaluate.outStepInterval`
    (Euler と共有) > n // 6。ns_output の間隔は判定の格子なので、n を割り切らなければ直さず ValueError。"""
    if nsteps:
        n, src = int(nsteps), "引数 (--steps)"
    elif p.evaluate.get("ns_steps") is not None:
        n, src = int(p.evaluate["ns_steps"]), "evaluate.ns_steps"
    else:
        n, src = NS_STEPS_DEFAULT, f"既定 {NS_STEPS_DEFAULT}"
    no = _ns_output_spec(p)
    if no is not None and no["outStepInterval"] is not None:
        # NS だけの出力の間隔 (evaluate.ns_output; plan tooling-settle-convergence-and-monitor §4.5)。判定の格子なので黙って直さない
        out_int = no["outStepInterval"]
        if n % out_int:
            raise ValueError(f"NS の本段 {n} step ({src}) が evaluate.ns_output.outStepInterval {out_int} の倍数でない")
        return n, out_int, src
    out_int = int(p.evaluate.get("outStepInterval", max(n // 6, 1)))
    if n % out_int:
        # forge は outStepInterval の倍数でしか res を書かない → 最終 step の res が残るよう n を割り切る間隔に直す
        out_int = next(n // k for k in (3, 2, 4, 6, 1) if n % k == 0)
    return n, out_int, src


def _stop_rule_settings(p: Problem, nsteps, euler_ref, stop_rule: bool | None) -> dict:
    """本段の停止規則の設定 (`stop_rule.resolve_settings`; plan tooling-nozzle-cooled-wall-ns-settling §4.6・
    plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #36a)。**既定で使う** (問題 YAML に `evaluate.stop_rule` が無く引数 stop_rule も
    None)。使えない (euler_ref が無い・見つからない・固定のスナップショットの sha が合わない、出力の間隔が窓に 9 点を入れられない、
    上限が出力の間隔の倍数でない) とき:
      - **明示** (`evaluate.stop_rule: true` か引数 stop_rule=True; 記録 `explicit`): ValueError (prepare_ns は run dir を作る前に呼ぶ)。
      - **既定**: `use` = False と理由 (`reason`) を記録して返す (本段は fallback_steps の長さ; run_staged_ns が UNAVAILABLE を書く)。
    `evaluate.stop_rule: false`・引数 stop_rule=False は使わない (requested = False)。"""
    from . import stop_rule as SR
    n, out_int, n_src = _ns_steps(p, nsteps)
    sr = SR.resolve_settings(p.evaluate, nsteps_arg=nsteps, n_main=n, n_main_source=n_src, out_int=out_int,
                             euler_ref=euler_ref, requested=stop_rule)
    if sr["requested"]:
        why = [sr["reason"]] if not sr["use"] else []
        if euler_ref and SR.resolve_path(euler_ref) is None:
            why.append(f"euler_ref が見つからない ({euler_ref!r})")
        elif isinstance(euler_ref, dict):
            try:
                SR.check_euler_ref(euler_ref)            # 固定のスナップショットの sha256 (§5.1 #36e (6))
            except ValueError as e:
                why.append(str(e))
        if why and sr["explicit"]:
            raise ValueError(f"本段の停止規則を要求した ({sr['requested_source']}) が使えない: {'; '.join(why)}")
        if why:
            sr.update(use=False, reason="; ".join(why))
    return sr


def prepare_ns(problem_path, run_dir, nsteps=None, ic_from=None,
               dstar_csv=None, dstar_blend=(6.0, 9.0),
               delta_r_csv=None, offset: str = "normal", euler_ref=None,
               omega: float | None = None, prev_run=None, initializer=None,
               cfl_main: float | None = None, implicit_relax: float | None = None,
               stop_rule: bool | None = None) -> dict:
    r"""**物理壁 (inviscid + δ*) の RANS run** を準備する (A12)。

    plan: plans/active/tooling-nozzle-axismach-viscous-deltastar.md。

    - `dstar_csv` なし → v1: `feedback.deltastar.deltastar_offset` (相関) で法線オフセット
    - `dstar_csv` あり → v2: CSV (x_rt, dstar_rt) を補間して法線オフセット
      (CFD 抽出 δ* の固定点反復)。`dstar_blend=(x_lo, x_hi)` [r_t] で相関→CSV の
      smoothstep 区間を指定 (既定 (6, 9) = case/41 の規約。CSV が全域を覆うなら
      (-1, -0.5) 等で CSV を全域採用 [case/44 v3])
    - メッシュ/段階起動レシピは B8 系 NS v1 (run_0028-0030) で確立したものを流用。
      coarse 中継 (y+~50) は YAML の mesh/evaluate 設定だけの違いで同じ関数で作る
    - **`delta_r_csv` (生産経路, plans/active/tooling-nozzle-deltastar-core-matched-euler.md)**:
      半径方向補正 δ_r(x) [r_t] の CSV (列 x_rt, delta_r; `feedback.deltastar_loop` が作る
      `delta_r_next.csv`) を全域そのまま使い、`offset="radial"` で壁を作る。`dstar_csv`/`dstar_blend`
      (旧 v3 継ぎはぎ) とは排他。`euler_ref` (固定 Euler 参照 run) は帳簿用に prepare_info へ記録。
    - **`geometry.physical_wall_repr`** (plans/accepted/tooling-nozzle-wall-single-bspline.md §4.4): 物理壁の表現。キー無し = poly の壁は
      1 本で表せれば `single_bspline` (2026-10-07 ユーザ決定の既定)、それ以外 (ramp・1 本で表せない構成) は今の区分表現で壁ファイルを書かない /
      `legacy` = 今の区分表現 + 壁ファイル (復元に要る全要素) / `single_bspline` = 入口から出口まで 1 本の 5 次 B-spline に作り直した壁
      (`SingleBSplinePhysicalWall`) + 壁ファイル。1 本の B-spline か明示のときは壁ファイル
      `wall_repr.json` (形式の版・表現の種類・有効域・単位・設計壁と物理壁の係数) を run に置き、prepare_info.json の
      `physical_wall` にも写す。joint 壁 + 物理壁の解析経路 (offset radial) 専用で、それ以外にキーを書いたら例外。
    - **`geometry.pw_upstream`** (plans/accepted/tooling-nozzle-upstream-poly-and-throat-sizing.md §4.1): joint 壁の物理壁の上流の作り方。
      キー無し = `poly` (配管〜設計スロートの 5 次多項式、上流に δ_r を足さない) / `ramp` = 旧来の δ_r ランプ (`pw_ramp`)。
      `poly` と `pw_ramp` の併記・不正値・joint でない壁への `poly`・`ramp` + `single_bspline` は run dir を作る前に例外。
      解決済みの値は prepare_info.json の `pw_upstream` に、`poly` のゲートは `pw_upstream_gate` に書く。
    - **`spec.sizing`** (任意、§4.2): 寸法の決め方 {method: exit | throat, target_m}。prepare_info.json の `sizing` に、実際の物理壁の
      物理スロート半径・出口半径 [m] と目標との差を書く (キー無しは method = null = 未記録)。
    - **`evaluate.main_line`** (任意、plan tooling-nozzle-isothermal-wall-chain §4.8-2・§5.1 #30): 本段だけに B0 (ライン陰解法・方向別の
      擬似 dt の上限・implicitThermalJacobian・cfl・implicitRelax・リミッタの基準値・出力の量) を当てる (`_main_line_spec`・
      `_apply_main_line`)。B0 を当てる前の設定を `PRESTAGE_CONFIG` に置き (run_staged_ns が S1・S2 をそこから作る)、変換・初期場の
      処理はその設定のまま行い、最後に solverConfig.yaml を本段の設定に替える。`cfl_main`・`implicit_relax` との併記は run dir を
      作る前に例外。キー無しは今と同じ設定・同じファイル。prepare_info.json の `main_line` と `ns_stage_configs` (stages=full の
      S1・S2・本段の設定の sha256) に記録する。
    - **`geometry.physical_upstream`** (任意、plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #34): 物理壁だけの上流の形
      {r_inlet, L_U, L_pipe} (`_physical_upstream`)。設計壁・CFD ピンの照合は geometry の r_inlet・L_U・L_pipe のまま。
      δ_r は `delta_r_csv` で与えること (積分法の初期壁・相関 δ* の経路は設計壁の縮流部で数えるので、run dir を作る前に例外)。
      壁温分布 (Tw_profile) の区切りの点も設計壁の値なので、併用は run dir を作る前に例外。キー無しは今と同じ (prepare_info.json にも
      書かない)。指定したときは prepare_info.json の `physical_upstream` に解決済みの値と凍結した設計の値を書く。第一層の厚さの表
      (`mesh.wall_first_frac_table`) は局所半径との比なので、上流の半径を変えたら表も別に直す (このキーは格子の表を書き換えない)。
    - **本段の停止規則** (plan tooling-nozzle-cooled-wall-ns-settling §4.6・plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #36a、
      `stop_rule.py`; 旧 m9_watch と同等とは言わない新しい規則、`_stop_rule_settings`): **既定で使う** (問題 YAML に `evaluate.stop_rule` が
      無ければ使う、forge-f5 の版・ユーザ決定)。`evaluate.stop_rule: false` か引数 `stop_rule=False` (CLI `--no-stop-rule`) で使わない。
      `evaluate.stop_rule_cap` (上限 = 到達 + 保持を含む本段の step 数、既定 200,000; 引数 nsteps があればそれ)・`evaluate.stop_rule_thresholds`
      (閾値・窓・断面の比 x_fracs の上書き)・`evaluate.stop_rule_hold_steps` (到達の後の継続、既定 `stop_rule.HOLD_STEPS_DEFAULT` = 40,000)。
      使えない (euler_ref が無い・見つからない・固定の sha が合わない、出力の間隔、上限が出力の間隔の倍数でない、設計壁の出口
      `physical_wall.domain` の無い物理壁) とき: `evaluate.stop_rule: true` か引数 `stop_rule=True` の**明示**なら例外 (壁の出口の検査だけは
      物理壁を作った後なので run dir を作った後、他は run dir を作る前; M6 の版)、**既定**なら `use` = False と理由を記録して本段を
      fallback_steps (`--steps` > `evaluate.ns_steps` > 60,000) で回す (run_staged_ns が UNAVAILABLE を記録; forge-f5 の版)。決めた設定は
      prepare_info.json の `stop_rule` に書き、使うときは `run_staged_ns` が本段の nStepOuter を上限に書き換えて見張る (ここで書く solverConfig の
      nStepOuter は fallback_steps のまま。上限と違うと `ns_stage_configs.main` の sha256 は実際に回した本段と違うので、実際の本段は
      MAIN_STAGE_STARTED を見る)。θ_r の断面は設計壁の出口 (`physical_wall.domain[1]`) の比。
    - **`euler_ref`** (固定の Euler の参照): run dir (str; 今までどおり、抽出は最新の出力を読む) か、スナップショットの固定 {run, res, sha256}
      (plan tooling-nozzle-isothermal-wall-chain §5.1 #36e (6); `metrics.deltastar.resolve_euler_ref`)。固定は run dir を作る前に sha256 を照合し
      (違えば ValueError)、prepare_info.json の `euler_ref` (と停止規則の `stop_rule.euler_ref`) にそのまま残す。停止規則の帯・見張りの中の
      出力ごとの δ の記録・受け渡しの関門の抽出はその固定のファイルを読む。
    - **本段の監視用の出力と判定の基準** (plan tooling-settle-convergence-and-monitor §4.5・§5.1 #6; 要求カードの雛形が既定で書く。キー無しは今と同じ):
      `evaluate.ns_output` {outStepInterval, fieldOutStepInterval, boundary_h5, boundaryMassflux} (`_ns_output_spec`) — outStepInterval は NS の
      出力の間隔 (`_ns_steps`; evaluate.outStepInterval は Euler だけ)、boundary_h5 の境界は bcondConfig の outputHDFflg 1、fieldOutStepInterval・
      boundaryMassflux は**本段の設定だけ** (`_apply_ns_outputs`)。`evaluate.settle` {preset, ...} (`_settle_spec`) — run の設定で解決して
      `settle_criteria.yaml` (`settle_presets.write_settle_criteria`、`grid.segment` = 本段の開始の通知) と節点の表 `settle_nodes.txt`
      (`settle_presets.probe_node_ids`) を書き、要る出力 (output.nodeProbeFile・boundaryMassflux、等温壁なら interfaceDiag) を本段の設定に足す。
      解決できない設定 (凝縮のモーメントの列など) は例外 (`evaluate.settle: false` の明示だけで外す; 外した理由は prepare_info の `settle`)。
      本段の設定が前の段の元の設定と違えば (B0 か出力) PRESTAGE_CONFIG を書き、S1・S2 はそこから作る。判定の基準が無い run で
      fieldOutStepInterval を出力の間隔より粗くし、停止規則 (旧方式 = 出力ごとに全場を読む) を使う組み合わせは run dir を作る前に例外。
    - **`main_line.limiter_ref`** (同 plan §4.6・§5.1 #6b): 省略 = `ic_from` が NS の run なら**その基準値を継承** (`limiter_ref_inherit_from_ic`:
      元の run の limiter_ref.json の 4 つ、無ければ本段の設定の明示 ρ・P・a (+ L)、分からなければ run dir を作る前に ValueError)、`ic_from` が
      無い・Euler の run なら自動 (forge が本段の開始の場で決めて `limiter_ref.json` に書き、stage_manifest の本段の段に写る)、`auto` = 継承しない自動、
      `{ro, p, a[, L]}` = 明示、`{from_run: <run>}` = 別 run の記録を継承 (`limiter_ref_from_run`)。明示・from_run は初期場の元より優先。B0 の本段
      だけ (point の本段は今までどおり自動)。出所は prepare_info の `main_line.limiter_ref_source` と MAIN_STAGE_STARTED の `limiter_ref.source`。
      `evaluate.main_line: point` は B0 を使わない明示 (prepare_info の `main_line_mode`)。
    """
    from ..feedback.deltastar import _sutherland
    from . import stop_rule as SR
    from ..geometry.wall_axismach import PhysicalNozzleWall, SingleBSplinePhysicalWall, check_required_attrs, save_wall_file
    p = load_problem(problem_path)
    if p.type != "wind_tunnel_axisym_axismach":
        raise ValueError("runner_axismach は wind_tunnel_axisym_axismach 専用")
    # 物理壁の表現 (§4.4): 不正値と未対応の設計壁は run dir を作る前・設計チェーンの前に例外
    pw_repr = _physical_wall_repr(p.geometry)
    if pw_repr is not None and str(p.geometry.get("wall_repr", "interp")) != "joint":
        raise ValueError(f"geometry.physical_wall_repr: {pw_repr} は joint 壁 (wall_repr: joint) + 物理壁の解析経路専用 "
                         f"(wall_repr = {p.geometry.get('wall_repr', 'interp')!r})")
    # 上流の作り方 (pw_upstream) と寸法の記録 (spec.sizing): 不正値・併記は run dir を作る前に例外
    pwu = _pw_upstream(p.geometry)
    sizing = _sizing_spec(p.spec)
    # 種ごとの輸送物性 (#9b): TP の NS は gas.transport 必須。run dir を作る前・設計チェーンの前に検査する
    transport = p.transport_for_ns()
    # 節点数の上限 (plan tooling-nozzle-mesh-size-guard §4.2): 実効の節点数で run dir を作る前・設計チェーンの前に判定する
    # (節点数は scale によらないので scale 1 で数える。既定の ni・nj・第 1 層は下の格子と同じ ns_mesh_defaults)
    ns_mesh_defaults = {"ni": 561, "nj": 97, "wall_first_frac": 4.5e-5}
    nodes = check_node_budget(mesh_params(p, 1.0, **ns_mesh_defaults), p.mesh, "mesh")
    # 本段だけの B0 (§5.1 #30): 不正値・cfl_main / implicit_relax との併記は run dir を作る前に例外
    # (limiter_ref の省略 = 初期場の元の NS run から継承・無ければ自動、from_run の記録の読み・検査もここ; plan
    # tooling-settle-convergence-and-monitor §4.6・§5.1 #6b)
    ml_src = {}
    main_line = _main_line_spec(p, cfl_main, implicit_relax, source_out=ml_src, ic_from=ic_from)
    main_line_mode = _main_line_mode(p)
    # 本段の監視用の出力と判定の基準 (plan tooling-settle-convergence-and-monitor §4.5・§5.1 #6): キーとプリセットの名前、
    # 流量の収支の入口・出口の判別 (bcond の kind) は run dir を作る前に検査する。中身の解決は格子と設定ができてから
    ns_out = _ns_output_spec(p)
    settle_es, settle_rec = _settle_spec(p)
    bc_ns = _bcond_with_species(_bcond(p, euler=False, profile_ok=True), _tp_species_Y(p))
    if ns_out is not None and ns_out["boundary_h5"]:
        bc_ns = _bcond_set_output(bc_ns, ns_out["boundary_h5"])
    if settle_es is not None:
        import yaml as _yaml
        _settle_presets().classify_bconds(_yaml.safe_load(bc_ns))
    # 物理壁だけの上流の形 (§5.1 #34): 不正値・δ_r の経路 (delta_r_csv だけ)・壁温分布との併用は run dir を作る前に例外
    pu = _physical_upstream(p.geometry)
    if pu is not None:
        if delta_r_csv is None:
            raise ValueError("geometry.physical_upstream は delta_r_csv の経路だけ (積分法の初期壁・相関 δ* は設計壁の縮流部で数えるので、"
                             "物理壁の別の形の上流の境界層の履歴にならない)")
        if "Tw_profile" in p.wall_thermal:
            raise ValueError("geometry.physical_upstream と spec.wall_thermal.Tw_profile の併用は未対応 (壁温分布の区切りの点 −L_U は "
                             "設計壁の値なので、物理壁の縮流部の始まりと食い違う)")
    # 固定の Euler の参照 (§5.1 #36e (6)): スナップショットの指定 {run, res, sha256} は run dir を作る前に sha を照合する (ディレクトリだけは今と同じ)
    from ..metrics.deltastar import resolve_euler_ref
    euler_ref_rec = (resolve_euler_ref(euler_ref, verify=isinstance(euler_ref, dict))["record"] if euler_ref else None)
    # 本段の停止規則 (§5.1 #36a): 既定で使う。不正値と、明示 (evaluate.stop_rule: true) なのに使えない設定は run dir を作る前に例外
    # (既定なら使えない理由を記録して fallback_steps で回す)
    sr_settings = _stop_rule_settings(p, nsteps, euler_ref, stop_rule)
    # 判定の基準の無い run (evaluate.settle が無い・false) の停止規則は旧方式 (v1; 出力ごとに全場を読む)。全場の間隔を出力の間隔より
    # 粗くすると窓に全場が揃わないので、その組み合わせは run dir を作る前に止める (plan tooling-settle-convergence-and-monitor §4.5)
    if (settle_es is None and sr_settings.get("requested") and stop_rule is not False and ns_out is not None
            and ns_out["fieldOutStepInterval"] not in (None, _ns_steps(p, nsteps)[1])):
        raise ValueError("evaluate.ns_output.fieldOutStepInterval が出力の間隔と違うが、判定の基準 (evaluate.settle) が無いので本段の停止規則は"
                         "出力ごとに全場を読む旧方式になる — evaluate.settle を書くか、evaluate.stop_rule: false にするか、"
                         "fieldOutStepInterval を外す")
    run_dir = Path(run_dir)
    if run_dir.exists():
        raise FileExistsError(f"{run_dir} が既にある")
    d = design_chain(p)
    require_moc_gate(p, d)                        # MOC の不合格は run dir を作る前に止める (2026-10-07)
    run_dir.mkdir(parents=True, exist_ok=False)
    scale = float(p.spec["r_throat"])
    wall_inv = d["wall_inv"]                      # (n,4) [x,r,th,M] r_t 単位
    # 物理壁 (A13): 上流履歴込み δ* + 真のスロート探索 + 上流 Hermite 再生成。
    # v2 (dstar_csv) は「下流は CFD 抽出、x<6 は履歴相関、[6,9] smoothstep」の合成。
    dstar_x = None
    if dstar_csv is not None:
        tbl = np.loadtxt(dstar_csv, delimiter=",", skiprows=1)
        base = PhysicalNozzleWall(d["wall"], wall_inv, scale, float(p.spec["Pt"]),
                                  float(p.spec["Tt"]), _gam_or_gas(p), p.cp)

        b_lo, b_hi = float(dstar_blend[0]), float(dstar_blend[1])

        def dstar_x(x, _tbl=tbl, _corr=base._dstar_hist):
            x = np.asarray(x, dtype=float)
            w = np.clip((x - b_lo) / max(b_hi - b_lo, 1e-9), 0.0, 1.0)
            w = w * w * (3.0 - 2.0 * w)
            csv = np.interp(np.clip(x, _tbl[0, 0], _tbl[-1, 0]), _tbl[:, 0], _tbl[:, 1])
            return (1.0 - w) * _corr(x) + w * csv
    delta_r_x = None
    delta_r_table_ends = None          # 表の端の照合 (delta_r_csv のときだけ; §5.1 #3h)
    init_info = None
    init_cfg = initializer if initializer is not None else p.raw.get("deltastar_initializer")
    if init_cfg and delta_r_csv is None and dstar_csv is None:
        # 積分法初期壁 (plan §4.1): 初回 NS 専用。断熱 / 指定壁温は thermal_bc で。
        res_init, delta_r_x, init_info = integral_delta_r(p, d, init_cfg)
        offset = "radial"
        (run_dir / "delta_r_initial.csv").write_text("")   # 後で上書き (run_dir は下で作る)
    if delta_r_csv is not None:
        if dstar_csv is not None:
            raise ValueError("delta_r_csv と dstar_csv は排他")
        tbl_r = np.loadtxt(delta_r_csv, delimiter=",", skiprows=1)
        if not np.all(np.isfinite(tbl_r[:, 1])):
            raise ValueError("delta_r_csv に非有限値がある (deltastar_loop.extract_and_merge で前回値保持済みの CSV を渡す)")
        # 表の端が設計の入口端・出口から ε_f32 以内なら揃える (格子の x [m] / r_t の丸めで 1 ulp 短い表; plan
        # tooling-nozzle-card-local-mesh §5.1 #3h)。超えた端はそのまま (今までどおり壁の構築が判定する)
        from ..feedback.deltastar_loop import design_anchors, snap_table_ends_same_design
        xs_r, delta_r_table_ends = snap_table_ends_same_design(tbl_r[:, 0], design_anchors(d["wall"]), scale)
        delta_r_x = delta_r_from_table(xs_r, tbl_r[:, 1])
        offset = "radial"
    # 物理壁 (pw_ramp・pw_upstream・physical_wall_repr に従う。寸法の逆算と同じ構築; build_physical_wall)。ゲート不合格は例外
    wall = build_physical_wall(p, d, scale, delta_r_x=delta_r_x, dstar_x=dstar_x, offset=offset, pwu=pwu)
    if wall.pw_upstream is not None:
        check_required_attrs(wall)          # 方式別の必須属性 (ramp: ramp_gate / poly: upstream_gate・upstream_poly)
    if init_info is not None:
        np.savetxt(run_dir / "delta_r_initial.csv",
                   np.c_[res_init["x"], res_init["delta_r"], res_init["dstar_n"], res_init["theta"] * scale,
                         res_init["H"], res_init["N"], res_init["Cf"], res_init["M"], res_init["Tw"], res_init["Taw"]],
                   delimiter=",", comments="",
                   header="x_rt,delta_r_rt,dstar_n_rt,theta_m,H,N,Cf,M_e,Tw,Taw")
        (run_dir / "delta_r_initial.json").write_text(json.dumps(init_info, indent=1, default=str))
    if init_info is not None:
        dstar_src = f"integral_bl:{init_info['model']} thermal_bc={init_info['thermal_bc']} (radial)"
    elif delta_r_csv is not None:
        dstar_src = f"delta_r_csv:{delta_r_csv} (radial, core-matched Euler ref {euler_ref_rec}, omega {omega})"
    else:
        dstar_src = ("correlation_hist_v1" if dstar_csv is None else f"{dstar_csv} (blend {dstar_blend})") + f" ({offset})"
    msgs = wall.validate()
    if msgs:
        raise ValueError("物理壁フィルタ不合格: " + "; ".join(msgs))
    # 壁温分布 (plan tooling-nozzle-wall-temperature-profile §4.2): 表が最終の物理壁 (入口端〜出口端) を覆うか検査し、run に実行条件の
    # 正本 2 つ (wall_temperature_profile.json・wall_profile_3.csv; physID 3 = _bcond の壁) を書いて sha256 を prepare_info に記録する。
    # 一様・断熱は今までどおり p.wall_thermal をそのまま記録する
    wt_rec = p.wall_thermal
    if "Tw_profile" in wt_rec:
        x_tab, _ = to_x_m(wt_rec["Tw_profile"], _wt_anchors(p, d, scale))
        check_covers(x_tab, float(wall.x_in) * scale, float(wall.x_e) * scale, f"最終の物理壁 (入口端〜出口端、scale = {scale!r} m)")
        wt_rec = {**wt_rec, "files": write_run_files(run_dir, p, d, scale, wall_id=3)}
    pw_info = None
    # 実効の表現: 明示 (legacy / single_bspline) か、キー無しの poly の既定 (1 本で表せれば single_bspline)
    pw_source = "explicit" if pw_repr is not None else None
    if pw_repr is None and isinstance(wall, SingleBSplinePhysicalWall):
        pw_repr, pw_source = "single_bspline", "default"
    if pw_repr is not None:
        # 壁ファイル (保存した壁の復元の取り決め §4.2): 書いて読み直し、復元した物理壁が今の壁と一致することを確かめる
        wpath, wrec, wsha = save_wall_file(run_dir, wall, scale, pw_repr)
        pw_info = {"repr": pw_repr, "source": pw_source, "file": wpath.name, "sha256": wsha,
                   **{k: wrec[k] for k in ("format", "version", "units", "origin", "domain", "domain_m", "physical_wall")}}
        if pw_repr == "single_bspline":
            pw_info["fit"] = wall.fit_diag
    elif getattr(wall, "single_bspline_default_skipped", None):
        pw_info = {"repr": None, "source": "default (1 本で表せないので区分表現のまま)", "reason": wall.single_bspline_default_skipped}
    mp = mesh_params(p, scale, **ns_mesh_defaults)
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
    _check_generated_nodes(coords, nodes, "mesh")
    # msh の座標の桁数 (mesh.msh_digits、既定 10 = 従来どおり; 冷却壁の薄い第一層は 17、plan tooling-nozzle-isothermal-wall-chain §5.1 #16)
    write_msh41_2d(run_dir / "nozzle.msh", coords, quads, bedges, digits=int(p.mesh.get("msh_digits", 10)))
    law = d["law"]
    xs = np.linspace(d["x0"], d["x_E"], 400)
    tgt = [float(law(x)) if x >= d["x_A"] else float("nan") for x in xs]
    np.savetxt(run_dir / "target_axis_M.csv", np.c_[xs * scale, tgt],
               delimiter=",", header="x_m,M_target", comments="")
    np.savetxt(run_dir / "wall_design.csv",
               np.c_[d["wall_inv"] * [scale, scale, 1.0, 1.0]], delimiter=",",
               header="x_m,r_m,theta_rad,M_wall", comments="")
    xs_p = np.linspace(wall.x_throat, wall.x_e, 1500)
    np.savetxt(run_dir / "wall_physical.csv", np.c_[xs_p * scale, wall.r(xs_p) * scale],
               delimiter=",", header="x_m,r_m", comments="")
    # NS の本段の step 数 (2026-10-14 ユーザ決定「1 ok」: 何も指定しないときは手順書の標準の 60,000 step)。`evaluate.nStepOuter` は Euler (`prepare`) だけが読む
    # (雛形の 12,000 が NS にも効いて本段が短すぎた)。優先: 引数 nsteps (`--steps`) > `evaluate.ns_steps` > NS_STEPS_DEFAULT (`_ns_steps`)。
    # 停止規則を使うときは run_staged_ns が本段を上限 (stop_rule.cap) に書き換えるので、ここの n は使わない・使えないときの長さ
    n, out_int, n_src = _ns_steps(p, nsteps)     # forge は outStepInterval の倍数でしか res を書かないので n を割り切る間隔に直す
    if not nsteps and p.evaluate.get("nStepOuter") is not None and p.evaluate.get("ns_steps") is None:
        print(f"[prepare_ns] 本段 {n} step ({n_src})。evaluate.nStepOuter {p.evaluate['nStepOuter']} は Euler だけが読む (2026-10-14 から)", flush=True)
    cfl_main = float(cfl_main if cfl_main is not None else p.evaluate.get("cfl_main", 1.0))
    implicit_relax = implicit_relax if implicit_relax is not None else p.evaluate.get("implicit_relax")
    # 壁温分布なら壁の行に wallProfile: 1 (CSV は上で書いた。この runner だけが profile_ok=True)。evaluate.ns_output.boundary_h5 の
    # 境界は outputHDFflg 1 (plan tooling-settle-convergence-and-monitor §4.5; 全段に効く出力だけ — 前の段の出力は段の終わりに消える)
    (run_dir / "bcondConfig.yaml").write_text(bc_ns)
    (run_dir / "probe.yaml").write_text(PROBE_STUB)
    # 品質は cell 変換コピーで検査 (品質ツールは node CONNE 非対応)
    (run_dir / "solverConfig.yaml").write_text(
        _apply_gas_to_config(_config_euler(p, n, out_int, 4.0, 1), p, run_dir))
    run_converter([str(converter_path()), "nozzle.msh", "nozzle_qc.h5"], cwd=run_dir, env=_ENV)
    q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"),
                            "nozzle_qc.h5", "--ar-max", str(int(p.mesh.get("ar_max", 1000)))], cwd=run_dir, env=_ENV,
                       capture_output=True, text=True)
    (run_dir / "MESH_QUALITY.txt").write_text(
        "# cell 変換コピーで検査 (品質は primal の性質)\n" + q.stdout + q.stderr)
    (run_dir / "nozzle_qc.h5").unlink()
    if q.returncode != 0:
        raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
    # node/SST 変換 (config を先に書く — wall_dist は no-slip 壁で作られる)。main_line があっても変換・初期場の処理は B0 を
    # 当てる前の設定で行い、本段の設定 (cfg_main) への切り替えは最後に行う (§5.1 #30)
    cfg_ns, cfg_main = _ns_configs(p, run_dir, n, out_int, cfl_main, implicit_relax, main_line)
    # 本段の監視用の出力 (plan tooling-settle-convergence-and-monitor §4.1・§4.5): 判定の基準 (evaluate.settle) を run の設定で一度解決して
    # 要る出力 (節点・流量・界面の診断) を決め、evaluate.ns_output と合わせて本段の設定だけに書く。前提を満たせない設定 (凝縮のモーメントの
    # 列・非定常・入口か出口の判別ができない境界など) は settle_presets.resolve の例外で止まる (黙って外さない; 外すのは evaluate.settle: false)
    settle_pi = {"Md": d["Md"], "scale_m": scale, "mesh": {"ni": mp.ni, "nj": coords.shape[0] // mp.ni}, "x_E": d["x_E"],
                 **({"physical_wall": pw_info} if pw_info is not None else {})}
    settle_needs = None
    if settle_es is not None:
        import yaml as _yaml
        settle_needs = _settle_needs(_settle_presets().resolve(settle_es, settle_pi, _cfg_load(cfg_main), _yaml.safe_load(bc_ns)))
    cfg_main = _apply_ns_outputs(cfg_main, ns_out, settle_needs)
    (run_dir / "solverConfig.yaml").write_text(cfg_ns)
    run_converter([str(converter_path()), "nozzle.msh", "nozzle.h5"], cwd=run_dir, env=_ENV)
    # 記録だけ (NS の格子は mesh のまま、座標・接続は変えない; 2026-10-07 plan verification-case45-euler-total-enthalpy §4)
    mesh_rec = _mesh_record(p, mp, "mesh", coords.shape[0] // mp.ni, run_dir, coords, quads, bedges, nodes)
    paste_isentropic_ic(run_dir / "nozzle.h5", wall, scale,
                        float(p.spec["Pt"]), float(p.spec["Tt"]), p.gamma, p.cp,
                        gas=(None if str(p.evaluate.get('cfd_gas', 'same')) == 'cpg' else p.gas_model),
                        h_ref_T=float(p.evaluate.get('thermo_href_temp', 298.15)),
                        species_Y=_tp_species_Y(p))
    _stamp_ic_species(p, run_dir)          # 新規初期場の化学種属性 (ic_from ならこの後 interp_field が継承/消去を決める)
    if ic_from is not None:
        src = sorted(Path(ic_from).glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
        subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                        str(src), str(run_dir / "nozzle.h5")],
                       env=_ENV, check=True, capture_output=True, text=True)
        import h5py as _h5
        # ソースが Euler (k/omega~0) のときだけ入口相当を貼る (runner_wt と同じ規約)
        with _h5.File(src) as fs:
            om_src_max = float(np.max(fs["/VALUE/omega"][:])) if "/VALUE/omega" in fs else 0.0
        if om_src_max < 1.0:
            with _h5.File(run_dir / "nozzle.h5", "r+") as f:
                ro = f["/VALUE/ro"][:]
                f["/VALUE/roK"][:] = ro * 1.0
                f["/VALUE/roOmega"][:] = ro * 18000.0
        # 近壁 omega の粘性底層フロア omega = 6 nu / (beta1 y^2) (run_0030 の教訓)
        g = p.gamma
        R_gas = p.cp * (g - 1.0) / g
        with _h5.File(run_dir / "nozzle.h5", "r+") as f:
            ro = f["/VALUE/ro"][:]
            roUx, roUy = f["/VALUE/roUx"][:], f["/VALUE/roUy"][:]
            roe = f["/VALUE/roe"][:]
            wd = f["/VALUE/wall_dist"][:]
            T = np.maximum((roe - 0.5 * (roUx ** 2 + roUy ** 2) / np.maximum(ro, 1e-12))
                           * (g - 1.0) / (np.maximum(ro, 1e-12) * R_gas), 50.0)
            nu = _sutherland(T) / np.maximum(ro, 1e-12)
            om_floor = 6.0 * nu / (0.075 * np.maximum(wd, 1e-9) ** 2)
            f["/VALUE/roOmega"][:] = np.maximum(f["/VALUE/roOmega"][:], ro * om_floor)
    info = {"chain": "axismach_ns", "discretization": "node", "viscous": True,
            "dstar_source": dstar_src, "mdot_ratio_moc": None,
            "throat_physical": {"x": wall.x_throat, "r": wall.r_throat,
                                "kappa": wall.kappa_throat,
                                "dstar_throat_correlation": float(wall._dstar_hist(0.0)),
                                "delta_r_throat_applied": float(wall.r_throat - 1.0)},
            "offset": wall.offset_mode, "euler_ref": euler_ref_rec,
            "pw_ramp_gate": (wall.ramp_gate if wall.pw_upstream == "ramp" else None),
            "initializer": init_info,
            "omega": omega, "prev_run": (str(prev_run) if prev_run else None),
            "delta_r_csv": (str(delta_r_csv) if delta_r_csv else None),
            "x0": d["x0"], "x_A": d["x_A"], "x_E": d["x_E"], "L_c": d["L_c"],
            "Lc_mode": d["Lc_mode"], "Lc_solve": d["Lc_solve"],
            "anchor": list(d["anchor"]), "anchor_source": d["anchor_source"],
            "start_line": d["start_line"], "wall_mode": d["wall_mode"],
            "wall_repr": d["wall_repr"], "initial_line": d["initial_line"], "moc": d["moc"],
            "wall_fit": d["wall_fit"],      # 設計壁の当てはめ (joint: spline・mono_r2 ほか; plan tooling-nozzle-throat-monotone-r2 §5.1 #5 M4)
            "Md": d["Md"], "Md_moc_offset": d["Md_moc_offset"], "R": d["R"],
            "qa": {k: v for k, v in d["qa"].items() if k != "violations"},
            "nStepOuter": n, "nStepOuter_source": n_src, "cfl_main": cfl_main, "implicit_relax": implicit_relax, "scale_m": scale,
            "stop_rule": sr_settings,       # 停止規則 (plan tooling-nozzle-cooled-wall-ns-settling §4.6): 使うか・上限・閾値・使わない理由
            "wall_thermal": wt_rec,         # 分布なら Tw_profile と正本 2 つの sha256 (files) も
            "ic_from": str(ic_from) if ic_from else None,
            "mesh": {**mesh_rec, "axis_gap_frac": mp.axis_gap_frac}}
    if transport is not None:
        # 来歴: 種ごとの輸送物性の指定 (solverConfig の physProp.transport と同じ; 解決結果はソルバの resolved_species 記録)
        info["transport"] = {"source": "gas.transport", "viscMethod": 2, "models": transport}
    if pw_info is not None:
        # 物理壁の表現と係数 (明示、またはキー無しの poly の既定で 1 本にしたとき。1 本で表せず区分表現のままのときは理由)
        info["physical_wall"] = pw_info
    # 上流の作り方の解決済みの値とゲート、寸法の記録 (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1・§4.2)
    info["pw_upstream"] = {"value": wall.pw_upstream, "source": pwu["source"], "requested": pwu["requested"]}
    info["pw_upstream_gate"] = (wall.upstream_gate if wall.pw_upstream == "poly" else None)
    if pu is not None:
        # 物理壁だけの上流の形 (§5.1 #34): 問題の値と、壁が使った値 (名目の U→T Hermite・凍結した設計の値)。キー無しは書かない
        info["physical_upstream"] = {"requested": pu, "wall": wall.physical_upstream}
    info["sizing"] = _sizing_record(sizing, wall, scale)
    if delta_r_table_ends is not None:
        info["delta_r_table_ends"] = delta_r_table_ends
    info["environment"] = environment_record()     # 照合は同じ環境どうしで (2026-10-07)
    # 本段だけの B0 (§5.1 #30) と本段の監視用の出力 (plan tooling-settle-convergence-and-monitor §4.5): 前の段の元の設定を PRESTAGE_CONFIG に
    # 置き、solverConfig.yaml を本段の設定に替える (cfl_main・implicit_relax の記録は本段の値にする)。どちらも無ければファイルを書かない
    info["main_line"] = None
    info["main_line_mode"] = main_line_mode
    if cfg_main != cfg_ns:
        (run_dir / PRESTAGE_CONFIG).write_text(cfg_ns)
        (run_dir / "solverConfig.yaml").write_text(cfg_main)
    if main_line is not None:
        info["cfl_main"], info["implicit_relax"] = main_line["cfl"], main_line["relax"]
        info["main_line"] = {**main_line, "source": "evaluate.main_line", "prestage_config": PRESTAGE_CONFIG,
                             "limiter_ref_source": ml_src.get("limiter_ref")}
    info["ns_output"] = (None if ns_out is None else
                         {**ns_out, "outStepInterval_effective": out_int, "source": "evaluate.ns_output", "main_stage_only": True})
    info["ns_stage_configs"] = ns_stage_sha256(cfg_main, cfg_ns if cfg_main != cfg_ns else None)
    # 判定の基準 settle_criteria.yaml と節点の表 settle_nodes.txt (plan tooling-settle-convergence-and-monitor §4.2・§5.1 #6)
    if settle_es is not None:
        settle_rec = {**settle_rec, **write_settle_files(run_dir, settle_es, info, cfg_main, bc_ns)}
    info["settle"] = settle_rec
    # 本段の停止規則 (§5.1 #36a; run_staged_ns が読む)。θ_r の断面は設計壁の出口 (physical_wall.domain[1]) の比なので、壁ファイルを
    # 書かない区分表現の物理壁では使えない: 明示なら例外 (M6 の版)、既定なら使えない理由を記録して fallback_steps で回す (forge-f5 の版)
    if sr_settings["use"] and not (isinstance(pw_info, dict) and isinstance(pw_info.get("domain"), (list, tuple))):
        why = ("θ_r の断面を設計壁の出口 (prepare_info.physical_wall.domain[1]) の比で決めるが、この物理壁は壁ファイルを書かない"
               "区分表現 (physical_wall.domain が無い)")
        if sr_settings["explicit"]:
            raise ValueError(f"本段の停止規則は {why} — physical_wall_repr を指定する")
        sr_settings.update(use=False, reason=why)
    info["stop_rule"] = sr_settings
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1))
    return info


def _sizing_record(sizing: dict | None, wall, scale: float) -> dict:
    """prepare_info.json の `sizing`: 寸法の決め方 (spec.sizing、無ければ未記録) と、実際の物理壁の物理スロート半径・出口半径 [m]、
    目標との差 (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.2)。"""
    x_e = float(wall.x_e)
    got = {"throat": float(wall.r_throat) * scale, "exit": float(wall.r(np.array([x_e]))[0]) * scale}
    rec = {"method": None if sizing is None else sizing["method"], "target_m": None if sizing is None else sizing["target_m"],
           "note": None if sizing is None else sizing.get("note"), "r_throat_design_m": float(scale),
           "physical_throat_radius_m": got["throat"], "physical_throat_x_m": float(wall.x_throat) * scale,
           "exit_radius_m": got["exit"], "residual_m": None,
           "status": "未記録 (spec.sizing が無い)" if sizing is None else "記録あり"}
    if sizing is not None:
        rec["residual_m"] = got[sizing["method"]] - sizing["target_m"]
    return rec


def _first_order(cfg: str) -> str:
    """段階起動の前段用に空間 1 次化する (space.convMethod 1 / 2 → 0、0 はそのまま)。旧実装は `convMethod: 1` しか置換せず、
    2 次 (`convMethod: 2`) の config では前段が 2 次のまま回っていた (plan tooling-rerun-conditions §4.9)。
    YAML 上の位置で読み書きする (空白数・書式に依らない; codex result 段 #2)。convMethod が無い・0/1/2 以外は ValueError。"""
    conv = _cfg_num(_cfg_get(_cfg_load(cfg), "convMethod"))
    if conv == 0.0:
        return cfg
    if conv not in (1.0, 2.0):
        raise ValueError(f"space.convMethod = {_cfg_get(_cfg_load(cfg), 'convMethod')!r} — 前段を 1 次化できない")
    return _cfg_set(cfg, {"convMethod": 0})


def stage_gate(res_h5, cfg_text: str) -> list:
    """段終了ゲート: 段の最終 res で必要保存量 (config から決める; `rerun_conditions.required_conserved_from_cfg`)
    が揃い、有限で ρ>0 か。問題のリストを返す (空なら次段へ進んでよい)。
    forge の rc と最終 step だけを見ていた旧実装は、非有限の場を `restart_field` で次段の初期場へ写しえた (§4.9)。"""
    if str(FORGE_TOOLS) not in sys.path:
        sys.path.insert(0, str(FORGE_TOOLS))
    from rerun_conditions import field_problems, required_conserved_from_cfg
    req = required_conserved_from_cfg(_cfg_load(cfg_text))
    probs, _ = field_problems(res_h5, req, species_bounds=False)
    return probs


def _stage_manifest_cls():
    if str(FORGE_TOOLS) not in sys.path:
        sys.path.insert(0, str(FORGE_TOOLS))
    from stage_manifest import StageManifest
    return StageManifest


def _ns_prestage(cfg: str, cfl: float) -> str:
    """前段: CFL を cfl に、空間 1 次、nStepInner 5 → 10 (従来どおり; 5 以外は据え置き)。"""
    cfg = _first_order(_cfg_set(cfg, {"cfl": cfl, "cfl_pseudo": cfl}))
    if _cfg_num(_cfg_get(_cfg_load(cfg), "nStepInner")) == 5.0:
        cfg = _cfg_set(cfg, {"nStepInner": 10})
    return cfg


def ns_stage_configs(cfg_main: str, stages: str = "full", ramp=None, ramp_steps: int = 1000, cfg_pre: str | None = None):
    """`run_staged_ns` の段の設定を、forge を回さずに導く。戻り: (前の段 [(tag, 実効設定, step 数, 起動前の要求値)],
    本段の実効設定, 本段の起動前の要求値)。前の段の実効設定は step 数・出力の間隔まで書いた起動直前の文字列。
    `cfg_pre` (main_line・本段の監視用の出力の run の PRESTAGE_CONFIG = B0・出力を当てる前の設定) があれば S1・S2 はそこから作り
    (§5.1 #30)、本段との差が B0 のパス・本段の出力のパス (_NS_OUTPUT_PATHS; plan tooling-settle-convergence-and-monitor §4.5) と
    段の長さだけであること、S1・S2 に B0・本段の出力だけのキーが残っていないことを検査する
    (違えば ValueError)。`cfg_pre` があるときの stages="ramp" は定義していないので ValueError。"""
    n_main = _cfg_value(cfg_main, "nStepOuter")
    main_want = {k: _cfg_value(cfg_main, k) for k in ("convMethod", "cfl", "cfl_pseudo", "outStepInterval")}
    if cfg_pre is not None:
        len_paths = {("time", "last", "nStepOuter"), ("time", "outStepInterval")}
        stray = sorted(set(_cfg_diff_paths(_cfg_load(cfg_pre), _cfg_load(cfg_main))) - _MAIN_LINE_PATHS - _NS_OUTPUT_PATHS - len_paths,
                       key=str)
        if stray:
            raise ValueError(f"{PRESTAGE_CONFIG} と本段の設定の差が B0 (main_line)・本段の出力のパスと段の長さ以外にもある: "
                             + ", ".join(".".join(map(str, pth)) for pth in stray) + " — 段の設定を作らない")
        if stages == "ramp":
            raise ValueError("main_line (前の段の元の設定がある run) の stages='ramp' は定義していない — 'full' か 'none' を使う")
    pre = []

    def _add(cfg, nsteps, tag, want):
        cfg = _cfg_set(cfg, {"nStepOuter": nsteps, "outStepInterval": nsteps})
        pre.append((tag, cfg, nsteps, {**want, "nStepOuter": nsteps, "outStepInterval": nsteps}))

    if stages == "full":
        src = cfg_main if cfg_pre is None else cfg_pre
        _add(_ns_prestage(src, 0.5), 3000, "S1_soft", {"convMethod": 0, "cfl": 0.5, "cfl_pseudo": 0.5})
        _add(_ns_prestage(src, 1.0), 3000, "S2_mid", {"convMethod": 0, "cfl": 1.0, "cfl_pseudo": 1.0})
        if cfg_pre is not None:
            # S1・S2 にライン陰解法・上限・ITJ・リミッタの基準値・全場の間隔が無く、出力の指定は前の段の元の設定のまま
            # (本段の extraFields・節点・流量・界面の診断が無い)
            out_pre = _cfg_path_get(_cfg_load(cfg_pre), ("output",))
            for tag, cfg, _, _ in pre:
                doc = _cfg_load(cfg)
                left = [".".join(pth) for pth in sorted(pth for pth in _MAIN_LINE_ONLY_PATHS | _NS_OUTPUT_PATHS if pth[0] != "output")
                        if _cfg_path_get(doc, pth) is not _CFG_MISSING]       # output は下で前の段の元の設定と丸ごと比べる
                if _cfg_path_get(doc, ("output",)) != out_pre:
                    left.append("output")
                if left:
                    raise ValueError(f"前の段 {tag} の設定に B0 のキー {left} が残っている — 段の設定を作らない")
    elif stages == "ramp":
        ramp = tuple(ramp or (1.0, 2.0, 3.5))
        for i, c in enumerate(ramp):
            _add(_cfg_set(cfg_main, {"cfl": c, "cfl_pseudo": c}), int(ramp_steps), f"R{i + 1}_cfl{c:g}",
                 {"convMethod": main_want["convMethod"], "cfl": c, "cfl_pseudo": c})
        n_main = max(n_main - int(ramp_steps) * len(ramp), int(ramp_steps))
        # 最終 res が書かれるよう outStepInterval の倍数に丸める (forge は outStepInterval の倍数でしか res を書かない)
        out_int = main_want["outStepInterval"]
        n_main = max((n_main // out_int) * out_int, out_int)
        cfg_main = _cfg_set(cfg_main, {"nStepOuter": n_main})
    elif stages != "none":
        raise ValueError("stages は 'full' / 'none' / 'ramp'")
    return pre, cfg_main, {**main_want, "nStepOuter": n_main}


def _sha256_bytes(b: bytes) -> str:
    import hashlib
    return hashlib.sha256(b).hexdigest()


def _sha256_file(path) -> str | None:
    """ファイルの sha256 (無ければ None)。"""
    import hashlib
    path = Path(path)
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ns_stage_sha256(cfg_main: str, cfg_pre: str | None = None) -> dict:
    """prepare_info.json の `ns_stage_configs`: prepare の時点の本段の設定から stages="full" で導いた S1・S2・本段の
    設定の sha256 (`ns_stage_configs` と同じ導出、UTF-8 の文字列の sha256 = 書いたファイルの sha256)。呼び出し側が後で
    solverConfig.yaml を書き換えたら実際に回した設定は MAIN_STAGE_STARTED (`previous_stages`・`main_config`) を見る。"""
    pre, main, _ = ns_stage_configs(cfg_main, "full", cfg_pre=cfg_pre)
    return {"stages": "full", "at": "prepare",
            **{tag: _sha256_bytes(cfg.encode()) for tag, cfg, _, _ in pre},
            "main": _sha256_bytes(main.encode()),
            "prestage_config": (None if cfg_pre is None else _sha256_bytes(cfg_pre.encode()))}


def _write_main_stage_started(run_dir: Path, cfg_main: str, stages: str, done: list, cfg_pre: str | None) -> dict:
    """本段の通知 MAIN_STAGE_STARTED (JSON) を書く。本段の設定への切り替え (solverConfig.yaml) と前の段の出力の処理 (最終場の
    restart と res_* の削除) の後、forge を起動する直前に呼ぶ。外の見張り (case/45 の m9_watch.py など) はこのファイルが
    出てから起動する (前の段の res_3000.h5 を本段の出力と取り違えない; codex diagnose 2026-10-13)。一時ファイルに書いて
    rename する (読みかけの JSON を見せない)。中身: 時刻、本段の設定と本段の初期の場 (mesh.valueFileName) の sha256、
    前の段の tag・step 数・設定の sha256、前の段の元の設定 (main_line の run だけ)、この時点で残っている res_*。
    判定区間の識別 (plan tooling-settle-convergence-and-monitor §4.2・§5.1 #6): 本段の設定の FNV-1a 64 (`main_config.fnv`; forge_launches.jsonl・
    limiter_ref.json・stage_manifest の本段の段の cfg_fnv と同じ定義) と、判定の基準 settle_criteria.yaml の sha256 (`settle_criteria`; ある
    run だけ) を書く。基準の `grid.segment` (準備のときに書く) はこの通知を指す。前の段のリミッタの基準値の記録 limiter_ref.json
    (前の段の stage_manifest の段に写し済み) は消す (本段の forge が書き直すまで前の段の値を本段のものと読ませない)。
    リミッタの基準値の出所 (§5.1 #6b): prepare_info.json の `main_line.limiter_ref_source` (継承元の run・ファイル・sha256・値、自動・明示の別) を
    `limiter_ref.source` に写し、本段の設定の明示の値 (`config_values`) とその値が同じか (`source_matches_config`; 自動なら明示が無いこと) を書く
    (prepare_info・main_line の無い run は source null・一致 null)。"""
    import datetime
    import os
    import time
    doc_main = _cfg_load(cfg_main)
    vf = str((doc_main.get("mesh") or {}).get("valueFileName", "nozzle.h5"))
    SMt = _limiter_tools()
    pi = run_dir / "prepare_info.json"
    ml_pi = json.loads(pi.read_text()).get("main_line") if pi.is_file() else None
    lim_src = ml_pi.get("limiter_ref_source") if isinstance(ml_pi, dict) else None
    sp_main = doc_main.get("space") if isinstance(doc_main.get("space"), dict) else {}
    cfg_vals = {k: sp_main[k] for k in SMt.LIMITER_REF_NAMES if k in sp_main}
    lim_match = (None if not isinstance(lim_src, dict)
                 else cfg_vals == ({} if lim_src.get("mode") == "auto" else (lim_src.get("values") or {})))
    lim_old = run_dir / SMt.LIMITER_REF_FILE
    lim_removed = lim_old.is_file()
    lim_old.unlink(missing_ok=True)
    crit = run_dir / _settle_presets().criteria_file_name()
    rec = {"event": "main_stage_started", "time": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
           "time_unix": time.time(), "stages": stages,
           "main_config": {"file": "solverConfig.yaml", "sha256": _sha256_file(run_dir / "solverConfig.yaml"),
                           "nStepOuter": _cfg_value(cfg_main, "nStepOuter"), "fnv": SMt.fnv1a64(cfg_main)},
           "initial_field": {"file": vf, "sha256": _sha256_file(run_dir / vf)},
           "previous_stages": list(done),
           "prestage_config": (None if cfg_pre is None
                               else {"file": PRESTAGE_CONFIG, "sha256": _sha256_file(run_dir / PRESTAGE_CONFIG)}),
           "res_files_present": sorted(f.name for f in run_dir.glob("res_*")),
           "start_step": 0,
           "settle_criteria": ({"file": crit.name, "sha256": _sha256_file(crit)} if crit.is_file() else None),
           "limiter_ref": {"file": SMt.LIMITER_REF_FILE, "previous_removed": lim_removed,
                           "note": "本段の forge が起動時に書く (cfg_fnv = main_config.fnv)。stage_manifest の本段の段の limiter_ref に写る",
                           "source": lim_src, "config_values": cfg_vals, "source_matches_config": lim_match}}
    tmp = run_dir / (MAIN_STAGE_STARTED + ".tmp")
    tmp.write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    os.replace(tmp, run_dir / MAIN_STAGE_STARTED)
    return rec


class StagedRC(int):
    """`run_staged_ns` の戻り (plan tooling-nozzle-isothermal-wall-chain §5.1 #36a): 今までどおりの終了コード (int として比べる・書ける) に
    本段の停止の状態 `stop_rule` (STOP_RULE_STATE.json と同じ中身から version・plan を除いた dict; 停止規則を要求していない run は None) を
    付けたもの。**終了コードで停止の状態を決めない** (停止規則の本段は CENSORED・一部の UNAVAILABLE でも 0)。"""

    def __new__(cls, rc, stop_rule: dict | None = None):
        obj = super().__new__(cls, int(rc))
        obj.stop_rule = stop_rule
        return obj

    @property
    def status(self) -> str | None:
        """停止の状態 (REACHED / CENSORED / DIVERGED / EXEC_ERROR / DATA_ERROR / UNAVAILABLE / NOT_USED)。要求していない run は None。"""
        return None if self.stop_rule is None else self.stop_rule.get("status")


def _stop_rule_for_main(run_dir: Path, stop_rule: bool | None) -> tuple:
    """本段で停止規則を使うか: (prepare_info の stop_rule, 使う, 使えない理由)。
    - prepare_info・その stop_rule が無い (2026-10-14 より前の準備)・要求していない (`evaluate.stop_rule: false`・`--no-stop-rule` で準備)・
      引数 stop_rule=False は (sr, False, None) で今までどおり。
    - 要求したのに使えない (prepare の記録が use でない・euler_ref が見つからない・固定の sha が合わない・θ_r の断面の基準
      physical_wall.domain や構造格子が読めない) とき: **明示** (`stop_rule.is_explicit`: `evaluate.stop_rule: true`・引数 True で準備) なら
      RuntimeError — S1 の前・forge を起動する前 (M6 の版)。**既定** (キー無しで準備) なら (sr, False, 理由) — 本段は prepare の長さ
      (fallback_steps) で回し、UNAVAILABLE を記録する (forge-f5 の版)。
    - stop_rule=True は prepare が要求していなければ ValueError (上限・閾値の記録が無い)。"""
    from . import stop_rule as SR
    pi = run_dir / "prepare_info.json"
    sr = json.loads(pi.read_text()).get("stop_rule") if pi.exists() else None
    requested = isinstance(sr, dict) and bool(sr.get("requested"))
    if stop_rule is True and not requested:
        raise ValueError(f"{run_dir}: run_staged_ns(stop_rule=True) だが prepare_info.json に停止規則の要求が無い "
                         "(問題 YAML の evaluate.stop_rule を false にしないで準備する)")
    if not requested or stop_rule is False:
        return sr, False, None
    why = None
    if not sr.get("use"):
        why = f"prepare で使えないと記録 ({sr.get('reason')})"
    elif SR.resolve_path(sr.get("euler_ref")) is None:
        why = f"euler_ref が見つからない ({sr.get('euler_ref')!r})"
    elif isinstance(sr.get("euler_ref"), dict):
        try:
            SR.check_euler_ref(sr["euler_ref"])          # 固定のスナップショットの sha256 (§5.1 #36e (6))
        except ValueError as ex:
            why = f"Euler の参照の固定が合わない ({ex})"
    if why is None:
        try:
            SR.RunGeometry(run_dir)
        except (ValueError, KeyError, OSError) as ex:
            why = f"格子・断面の基準を読めない ({type(ex).__name__}: {ex})"
    if why is None:
        return sr, True, None
    if SR.is_explicit(sr) or stop_rule is True:
        raise RuntimeError(f"{run_dir}: 本段の停止規則を明示で要求したが使えない: {why} — forge を起動しない")
    return sr, False, why


def _stop_rule_v2_for_main(run_dir: Path, stop_rule: bool | None, cfg_main: str) -> dict | None:
    """本段で停止規則 v2 (plan tooling-settle-convergence-and-monitor §4.4・§4.8 の 1・§5.1 #5、`stop_rule.watch_main_stage_v2`) を使うか。
    run に settle_criteria.yaml があり、prepare_info の stop_rule が要求していて、引数 stop_rule が False でなければ v2 — 戻り: prepare_info の
    stop_rule (本段の上限 cap は旧方式と同じ `--steps` > `evaluate.stop_rule_cap` > 200,000)。それ以外は None (旧方式の `_stop_rule_for_main`
    が決める: settle_criteria.yaml の無い run、要求なし・stop_rule=False は今までどおり)。
    v2 で使えない (`stop_rule.preflight_v2`: 制御のファイル STOP_REQUEST・STOPPED_AT・SNAPSHOT_REQUEST の残り (消さない)・基準の誤用や
    未設定・出力の設定・節点の表) とき、上限 cap の記録が無いときは、明示・既定を問わず S1 の前・forge を起動する前に RuntimeError
    (settle_criteria.yaml を置いた run は v2 を明示したとみなし、旧方式の UNAVAILABLE + fallback_steps には落とさない。v2 は Euler の参照を
    判定に使わないので euler_ref が無くても止めない — δ_E の抽出だけを飛ばして記録する)。stop_rule=True で prepare が要求していなければ
    ValueError (旧方式と同じ)。"""
    from . import stop_rule as SR
    if not SR.criteria_present(run_dir):
        return None
    pi = run_dir / "prepare_info.json"
    sr = json.loads(pi.read_text()).get("stop_rule") if pi.exists() else None
    requested = isinstance(sr, dict) and bool(sr.get("requested"))
    if stop_rule is True and not requested:
        raise ValueError(f"{run_dir}: run_staged_ns(stop_rule=True) だが prepare_info.json に停止規則の要求が無い "
                         "(問題 YAML の evaluate.stop_rule を false にしないで準備する)")
    if not requested or stop_rule is False:
        return None
    if sr.get("cap") is None:
        raise RuntimeError(f"{run_dir}: {SR.CRITERIA_FILE} の run (停止規則 v2) だが prepare_info.json の stop_rule に上限 cap が無い "
                           "— forge を起動しない")
    pf = SR.preflight_v2(run_dir, _cfg_load(cfg_main), int(sr["cap"]))
    if pf["problems"]:
        raise RuntimeError(f"{run_dir}: {SR.CRITERIA_FILE} の run (停止規則 v2) だが使えない — forge を起動しない: " + "; ".join(pf["problems"]))
    return sr


def _stop_state_record(run_dir: Path, w: dict, sr: dict) -> dict:
    """見張り (`stop_rule.watch_main_stage`) の戻りと詳しい記録 stop_rule.json から、STOP_RULE_STATE.json に書く終了の状態を作る。"""
    from . import stop_rule as SR
    rec = json.loads((run_dir / SR.RECORD).read_text())
    return {"method": "v1", "status": w["status"], "stop_step": w["stop_step"], "last_output": rec.get("last_output"),
            "forge_rc": w["forge_rc"], "run_case_rc": rec.get("run_case_rc"), "rc": int(w["rc"]),
            "cap": int(sr["cap"]), "fallback_steps": int(sr["fallback_steps"]), "hold_steps": int(sr.get("hold_steps") or 0),
            "hold": rec.get("hold"), "stop_sha256": rec.get("stop_sha256"), "note": rec.get("note"), "unavailable": rec.get("unavailable"),
            "record": SR.RECORD, "record_sha256": _sha256_file(run_dir / SR.RECORD), "finished": rec.get("finished")}


def run_staged_ns(run_dir, stages: str = "full", ramp=None, ramp_steps: int = 1000, stop_rule: bool | None = None,
                  stop_rule_watch: dict | None = None) -> int:
    """NS の起動。stages:
    - "full" (既定・run_0030 レシピ): soft (1次 cfl0.5 ni10, 3000 step) → mid (1次 cfl1, 3000) → 本段 (2次 cfl_main)。
    - "none": 本段だけ (収束済み NS 場からの warm start 用)。
    - "ramp": 2 次のまま cfl を `ramp` (例 (1, 2, 3.5)) の順に各 ramp_steps だけ回して本段 cfl_main へ
      (forge に CFL ランプ機能は無いので restart で段階化する。本段の step 数はランプ分を差し引く)。
    各段の最終場を IC に引き継ぐ。段の設定は `ns_stage_configs` で導く。run に PRESTAGE_CONFIG (prepare_ns が
    `evaluate.main_line` のときに書く、B0 を当てる前の設定) があれば S1・S2 はそこから作り、本段だけが B0 になる (§5.1 #30)。
    段ごとの実効設定を `stage_manifest.json` に記録し (`check_convergence.py --segment` 用)、段の残差履歴を
    `residual_history_<tag>.csv` に残す (段の res_* は従来どおり消す)。段の最終 res の必要保存量が非有限・ρ≤0 なら
    次段へ進まず RuntimeError (段終了ゲート `stage_gate`; plan tooling-rerun-conditions §4.9)。
    本段の設定への切り替えと前の段の出力の処理の後、forge を起動する直前に通知 MAIN_STAGE_STARTED を書く
    (`_write_main_stage_started`; 前の呼び出しの通知は最初に消す)。
    壁温分布の run は**各段の起動前** (前段の restart の後を含む) に run の正本 (wall_temperature_profile.json・wall_profile_<id>.csv) を
    prepare_info.json の sha256 と照合し、合わなければ config を書く前・forge を起動する前に ValueError で止める
    (`wall_temperature.load_run_profile`; plan tooling-nozzle-wall-temperature-profile §4.2)。一様・断熱の run は何もしない。

    **本段の停止規則** (plan tooling-nozzle-cooled-wall-ns-settling §4.6・plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #36a、`stop_rule.py`;
    旧 m9_watch と同等とは言わない新しい規則): prepare_info.json の `stop_rule` が要求していて使える (`use`) なら、本段の nStepOuter を上限
    (`cap`: `--steps` > `evaluate.stop_rule_cap` > 200,000、到達 + 保持を含む) に書き換えて solverConfig・stage_manifest.json・MAIN_STAGE_STARTED
    に記録し、forge を子プロセスで持ったまま見張る (`stop_rule.watch_main_stage`; REACHED・DIVERGED で止めて後の出力を消す、詳しい記録は
    stop_rule.json)。**前の段 (S1・S2・ramp) には当てない**。本段の後に終了の状態を STOP_RULE_STATE.json に 1 度だけ書き、戻り値 (`StagedRC`) の
    `stop_rule` にも付ける。要求したのに使えない (`_stop_rule_for_main`) とき: **明示** (`evaluate.stop_rule: true` で準備) なら S1 の前・forge を
    起動する前に RuntimeError (M6 の版)、**既定** (キー無しで準備) なら stop_rule.json に UNAVAILABLE と理由を書き、本段は prepare の長さ
    (fallback_steps = `--steps` > `evaluate.ns_steps` > 60,000) で今までどおり回して、状態 UNAVAILABLE を STOP_RULE_STATE.json に書く (forge-f5 の版)。
    `stop_rule=False` (CLI `--no-stop-rule` の後に呼ぶとき) は要求を打ち消して今までの経路で回し、状態 NOT_USED を記録する。prepare_info.json に
    `stop_rule` の無い run (2026-10-14 より前の準備) と要求していない run (`evaluate.stop_rule: false`) は今までと同じ (STOP_RULE_STATE.json も
    書かない)。`stop_rule_watch` は見張りの時間の設定 (試験用; poll_s・settle_s ほか)。
    **停止規則 v2** (plan tooling-settle-convergence-and-monitor §4.4・§4.8 の 1・§5.1 #5): run に settle_criteria.yaml があり停止規則を要求して
    いれば、上の旧方式の代わりに `stop_rule.watch_main_stage_v2` で見張る (`_stop_rule_v2_for_main`; 使えなければ S1 の前に RuntimeError)。
    本段の nStepOuter を上限 cap に書き換えるのは同じ。3 条件 (残差・流量の収支・目的の量) の保持 → STOP_REQUEST → forge の協調停止
    (STOPPED_AT) → 止めた step で再評価 (SETTLE_VERDICT.json) → δ_E を 1 回。状態は REACHED / CENSORED / FAILED、STOP_RULE_STATE.json の
    `method` が v2 (旧方式は v1)。`stop_rule_watch` は v2 の見張りの時間の設定 (poll_s・stop_wait_s ほか)。
    戻り: forge の終了コード (`StagedRC`、int のまま比べられる)。停止規則のときは REACHED・CENSORED・一部の UNAVAILABLE = 0、
    それ以外は非 0 — **合否は終了コードでなく状態 (`rc.stop_rule["status"]`・STOP_RULE_STATE.json) で決める**。"""
    import shutil
    from . import stop_rule as SR
    run_dir = Path(run_dir)
    cfg_main = (run_dir / "solverConfig.yaml").read_text()
    pre_path = run_dir / PRESTAGE_CONFIG
    cfg_pre = pre_path.read_text() if pre_path.is_file() else None
    bc_path = run_dir / "bcondConfig.yaml"
    bc_text = bc_path.read_text() if bc_path.exists() else ""
    pre, cfg_main, main_want = ns_stage_configs(cfg_main, stages, ramp, ramp_steps, cfg_pre)
    # 停止規則 v2 (settle_criteria.yaml の run): 使えなければここで例外 (S1 の前・forge を起動しない)。無ければ旧方式
    sr_v2 = _stop_rule_v2_for_main(run_dir, stop_rule, cfg_main)
    if sr_v2 is not None:
        sr, sr_use, sr_why = sr_v2, False, None
    else:
        # 停止規則: 明示の要求なのに使えなければここで例外 (S1 の前・forge を起動しない)。既定の要求で使えなければ sr_why に理由 (本段で UNAVAILABLE)
        sr, sr_use, sr_why = _stop_rule_for_main(run_dir, stop_rule)
    (run_dir / MAIN_STAGE_STARTED).unlink(missing_ok=True)   # 前の呼び出しの通知で見張りを早く起動させない
    (run_dir / SR.STATE).unlink(missing_ok=True)             # 前の呼び出しの停止の状態を今回の結果と取り違えない
    sm = _stage_manifest_cls()(run_dir)

    def _check_wall_profile(tag):
        """壁温分布の正本の照合 (一様・断熱の run では load_run_profile が None を返すだけ)。"""
        try:
            load_run_profile(run_dir)
        except ValueError as ex:
            raise ValueError(f"段 {tag}: 壁温分布の照合で不合格 — forge を起動しない: {ex}") from ex

    def _record(cfg, tag):
        """段の残差履歴を段名つきで残し、manifest に段を足して書く (失敗した段も診断用に残す)。"""
        hist = run_dir / "residual_history.csv"
        if hist.exists():
            shutil.copy(hist, run_dir / f"residual_history_{tag}.csv")
        sm.add(tag, cfg, bc_text, history=f"residual_history_{tag}.csv")
        sm.write()

    def _stage(cfg, nsteps, tag, want):
        _check_wall_profile(tag)                                     # 壁温分布の正本 (前段の restart の後も) — config を書く前 (main の §4.2)
        _cfg_check(cfg, want, tag)                                   # 起動前に実効値を検査
        (run_dir / "solverConfig.yaml").write_text(cfg)
        sha = _sha256_file(run_dir / "solverConfig.yaml")
        (run_dir / "residual_history.csv").unlink(missing_ok=True)   # 前段の履歴を別段の名前で写さない
        rc = run_forge(run_dir)
        _record(cfg, tag)
        res = sorted(run_dir.glob("res_[0-9]*.h5"),
                     key=lambda f: int("".join(c for c in f.stem if c.isdigit())))
        if rc != 0 or not res or int("".join(c for c in res[-1].stem if c.isdigit())) < nsteps:
            raise RuntimeError(f"段階起動が失敗 (rc={rc}, res={res[-1].name if res else None})")
        probs = stage_gate(res[-1], cfg)
        if probs:
            raise RuntimeError(f"段 {tag} の最終場 {res[-1].name} が段終了ゲートで不合格 — 次段へ進まない:\n    "
                               + "\n    ".join(probs))
        _restart_same_mesh(res[-1], run_dir / "nozzle.h5")    # 同一メッシュ: index コピー (旧: interp_field.py)
        for f in run_dir.glob("res_*"):
            f.unlink()
        return {"tag": tag, "nsteps": nsteps, "config_sha256": sha}

    done = [_stage(cfg, nsteps, tag, want) for tag, cfg, nsteps, want in pre]
    if sr_use or sr_v2 is not None:
        # 停止規則の上限を本段の nStepOuter に (manifest・通知の本段の記録も同じ値になる)。ramp の差し引きはしない (上限は本段の step 数)
        cfg_main = _cfg_set(cfg_main, {"nStepOuter": int(sr["cap"])})
        main_want = {**main_want, "nStepOuter": int(sr["cap"])}
    _check_wall_profile("main")
    _cfg_check(cfg_main, main_want, "main")
    (run_dir / "solverConfig.yaml").write_text(cfg_main)
    (run_dir / "residual_history.csv").unlink(missing_ok=True)
    _write_main_stage_started(run_dir, cfg_main, stages, done, cfg_pre)
    state = None
    if sr_v2 is not None:
        # 停止規則 v2: 協調停止で止め、止めた step で再評価する (forge を kill しない; 非常時を除く)
        w = SR.watch_main_stage_v2(run_dir, sr, _cfg_load(cfg_main), **(stop_rule_watch or {}))
        rc, state = w["rc"], w["state"]
        print(f"[stop_rule v2] 本段: {w['status']} (停止の step {w['stop_step']}、forge rc {w['forge_rc']}、上限 {sr['cap']}) → "
              f"{run_dir / SR.STATE}", flush=True)
    elif sr_use:
        w = SR.watch_main_stage(run_dir, sr, **(stop_rule_watch or {}))
        rc = w["rc"]
        state = _stop_state_record(run_dir, w, sr)
        print(f"[stop_rule] 本段: {w['status']} (止めた/最後の出力 {w['stop_step']}、forge rc {w['forge_rc']}、上限 {sr['cap']}) → "
              f"{run_dir / SR.STATE}", flush=True)
    elif sr_why is not None:
        # 既定の要求 (キー無しで準備) で使えない: 理由を記録し、本段は prepare の長さ (fallback_steps) で今までどおり回す (forge-f5 の版)
        n_main = int(main_want["nStepOuter"])
        urec = {"version": 1, "plan": SR.PLAN, "settings": sr, "status": "UNAVAILABLE", "reason": sr_why,
                "action": f"停止規則を使わず本段を nStepOuter {n_main} で回す", "started": SR._now()}
        SR.save_record(run_dir, urec)
        print(f"[stop_rule] 使えない ({sr_why}) — 本段は {n_main} step", flush=True)
        rc = run_forge(run_dir)
        forge_rc = SR._forge_exit(run_dir)
        urec.update(forge_rc=forge_rc, run_case_rc=int(rc), finished=SR._now())
        SR.save_record(run_dir, urec)
        state = {"method": "v1", "status": "UNAVAILABLE", "reason": sr_why, "rc": int(rc), "forge_rc": forge_rc, "stop_step": None, "cap": None,
                 "fallback_steps": int(sr["fallback_steps"]), "hold_steps": int(sr.get("hold_steps") or 0),
                 "record": SR.RECORD, "record_sha256": _sha256_file(run_dir / SR.RECORD), "finished": urec["finished"]}
    else:
        rc = run_forge(run_dir)
        if sr is not None and sr.get("requested"):
            # prepare は要求したが引数 stop_rule=False で使わなかった: 打ち消したことを状態として残す (受け渡しの関門は NOT_USED と見る)
            state = {"method": "v2" if SR.criteria_present(run_dir) else "v1",
                     "status": "NOT_USED", "reason": "引数 stop_rule=False (prepare は要求した)", "rc": int(rc), "stop_step": None,
                     "cap": None, "fallback_steps": int(sr["fallback_steps"])}
    _record(cfg_main, "main")
    if state is not None:
        SR.write_state(run_dir, state)
    return StagedRC(rc, state)


if __name__ == "__main__":
    sys.exit(main())
