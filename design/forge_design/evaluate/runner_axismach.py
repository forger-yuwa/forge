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
from ..meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d
from ..probdef import Problem, dv_value, load_problem
from .ic import paste_isentropic_ic, stamp_isentropic_ic_species
from .runner import FORGE_TOOLS, PROBE_STUB, _ENV, converter_path, run_forge
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
    return {"python": platform.python_version(), "executable": _sys.executable, "numpy": np.__version__, "scipy": scipy.__version__,
            "host": platform.node()}


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
                        wall_normal_layer=m.get("wall_normal_layer"), scale=scale)


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
# mesh_euler に書けるキー (_mesh_params_from が読むもの + 品質検査の ar_max + node 固定の discretization)。それ以外は書き誤りとして止める
MESH_EULER_KEYS = ("ni", "nj", "wall_first_frac", "throat_refine", "throat_width", "wall_first_frac_throat",
                   "wall_first_blend_x0", "wall_first_blend_x1", "wall_first_up_x0", "wall_first_up_x1",
                   "axis_gap_frac", "axis_cap_frac", "ar_max", "discretization")


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


def _mesh_hashes(run_dir: Path, coords, quads, bedges) -> dict:
    """変換後の nozzle.h5 の座標 (/MESH/COORD) と接続 (_MESH_TOPOLOGY + 境界の位相) の sha256、nozzle.msh の sha256、
    生成器の出力 (coords float64・quads・境界辺; 変換器に依らない) の sha256。"""
    import hashlib

    import h5py

    def _upd(h, name, a):
        a = np.ascontiguousarray(a)
        h.update(name.encode() + str(a.dtype).encode() + str(a.shape).encode() + a.tobytes())

    with h5py.File(run_dir / "nozzle.h5", "r") as f:
        c = f["/MESH/COORD"][:]
        hc = hashlib.sha256()
        hc.update(str(c.dtype).encode() + str(c.shape).encode() + np.ascontiguousarray(c).tobytes())
        ht = hashlib.sha256()
        names = [t for t in _MESH_TOPOLOGY if t in f]
        names += sorted(f"BCONDS/{b}/{k}" for b in f["BCONDS"] for k in _MESH_BCOND_TOPO if k in f["BCONDS"][b])
        for n in names:
            a = np.asarray(f[n])
            ht.update(n.encode() + str(a.dtype).encode() + str(a.shape).encode() + a.tobytes())
    hg = hashlib.sha256()
    _upd(hg, "coords", np.asarray(coords, dtype=np.float64))
    _upd(hg, "quads", np.asarray(quads, dtype=np.int64))
    for g in ("inlet", "outlet", "wall", "axis"):
        _upd(hg, g, np.asarray(bedges[g], dtype=np.int64))
    hm = hashlib.sha256((run_dir / "nozzle.msh").read_bytes())
    return {"coord_sha256": hc.hexdigest(), "coord_dtype": str(c.dtype), "topology_sha256": ht.hexdigest(),
            "topology_datasets": names, "msh_sha256": hm.hexdigest(), "generator_sha256": hg.hexdigest()}


def _mesh_record(p, mp: Mesh2DParams, source: str, nj: int, run_dir: Path, coords, quads, bedges) -> dict:
    """prepare_info.json の mesh 欄: 従来の ni・nj・wall_first_frac (読む側: nozzle_report・metrics.deltastar・cfd_initial_line) に、
    全 Mesh2DParams・採用元のブロック名・問題に書かれたブロックの有無と中身・座標と接続のハッシュを足す
    (plan verification-case45-euler-total-enthalpy §4、2026-10-07)。"""
    import dataclasses
    return {"ni": mp.ni, "nj": int(nj), "wall_first_frac": mp.wall_first_frac,
            "source": source, "blocks_present": {b: (b in (p.raw or {})) for b in ("mesh", "mesh_euler")},
            "block": dict((p.raw or {}).get(source) or {}), "params": dataclasses.asdict(mp),
            "hashes": _mesh_hashes(run_dir, coords, quads, bedges)}


def prepare(problem_path, run_dir, nsteps=None, ic_from=None, cfl_main=None, implicit_relax=None) -> dict:
    p = load_problem(problem_path)
    if p.type != "wind_tunnel_axisym_axismach":
        raise ValueError("runner_axismach は wind_tunnel_axisym_axismach 専用")
    if _physical_wall_repr(p.geometry) == "single_bspline":
        # 物理壁の表現 (plan tooling-nozzle-wall-single-bspline §4.2) は NS の物理壁 (prepare_ns) だけ。Euler の設計壁で黙って無視しない
        raise ValueError("geometry.physical_wall_repr: single_bspline は prepare_ns (物理壁) 専用 — Euler の prepare には物理壁が無い")
    if "pw_upstream" in p.geometry:
        # 上流の作り方 (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1) は joint 壁の物理壁の解析経路だけ。
        # 値は検査し、poly の明示は例外 (Euler の設計壁で黙って無視しない)。ramp (移行で明示した旧来の作り方) は通す
        v = p.geometry["pw_upstream"]
        if not isinstance(v, str) or v not in PW_UPSTREAMS:
            raise ValueError(f"geometry.pw_upstream は {' / '.join(repr(c) for c in PW_UPSTREAMS)} のどれか (受け取った値: {v!r})")
        if v == "poly":
            raise ValueError("geometry.pw_upstream: poly は prepare_ns (joint 壁の物理壁) 専用 — Euler の prepare には物理壁が無い")
    # Euler の格子は mesh_euler (mesh は読まない)。無い・不正なら run dir を作る前・設計チェーンの前に止める (2026-10-07)
    m_eu = mesh_euler_block(p)
    run_dir = Path(run_dir)
    if run_dir.exists():
        raise FileExistsError(f"{run_dir} が既にある")
    d = design_chain(p)
    require_moc_gate(p, d)                        # MOC の不合格は run dir を作る前に止める (2026-10-07)
    run_dir.mkdir(parents=True, exist_ok=False)
    wall = d["wall"]
    scale = float(p.spec["r_throat"])
    mp = mesh_params_euler(p, scale)
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
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
    subprocess.run([str(converter_path()), "nozzle.msh", "nozzle_qc.h5"],
                   cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
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
    if cfl_main is not None:
        cfg_e = re.sub(r"cfl: [\d.]+, cfl_pseudo: [\d.]+", f"cfl: {float(cfl_main)}, cfl_pseudo: {float(cfl_main)}", cfg_e)
    if implicit_relax is not None:
        cfg_e = cfg_e.replace("blockDPLUR: 1,", f"blockDPLUR: 1, implicitRelax: {float(implicit_relax)},", 1)
    (run_dir / "solverConfig.yaml").write_text(cfg_e)
    subprocess.run([str(converter_path()), "nozzle.msh", "nozzle.h5"],
                   cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    # IC を入れる前 (座標・接続は以後変わらない)。nj は実際の値 (axis_gap_frac では導出値になる)
    mesh_rec = _mesh_record(p, mp, "mesh_euler", coords.shape[0] // mp.ni, run_dir, coords, quads, bedges)
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
            mf = massflow_ratio(run_dir, info["euler_ref"])
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
    rc = run_staged(a.run_dir, cfl_main=(a.cfl if a.cfl is not None else pp.evaluate.get("cfl_main")),
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


def integral_delta_r(p: Problem, d: dict, init_cfg: dict, scale: float | None = None, rtol: float | None = None):
    """積分法初期壁の δ_r (`prepare_ns` の initializer 経路): `integral_bl` → 5 次 P-spline 平滑化 → 壁に渡す δ_r 関数。
    戻り: (res_init, delta_r_x, init_info)。`prepare_ns` から切り出したもの (振る舞いは同一; plan
    tooling-nozzle-throat-monotone-r2 §6 S6 の形状ゲートが同じ経路で物理壁を作るために共有する)。
    scale: スロート半径 r_t [m] (None = spec.r_throat)。寸法の逆算 (`deltastar_loop.solve_rt`・`solve_rt_throat`) が、反復のたびに
    同じ経路 (k_f の cf_scale・熱条件・平滑化) で δ_r を作り直すために渡す (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.2)。
    rtol: **診断用** (plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U2c)。`integral_bl` の RK45 の相対許容差を明示の引数で
    注入する。None = 渡さない (`integral_bl` の既定 1e-6、今の振る舞いのまま)。YAML のキーにはしない。実効値は
    `res_init["solve_ivp"]["rtol"]` (solve_ivp に渡した値) に残る。"""
    from ..feedback.deltastar_integral import integral_bl
    scale = float(p.spec["r_throat"]) if scale is None else float(scale)
    wall_inv = d["wall_inv"]
    model = str(init_cfg.get("model", "contur"))
    if model not in ("contur", "contur_momentum_integral"):
        raise ValueError(f"deltastar_initializer.model={model!r} は未対応 (contur のみ)")
    # 熱境界条件は spec.wall_thermal が単一ソース (plan tooling-nozzle-isothermal-wall-chain §4.1)。
    # initializer/YAML の thermal_bc 指定は無視し、食い違えば警告する (NS と積分法が別の壁温を読む状態を作らない)。
    tbc = p.wall_thermal_bc_integral
    if init_cfg.get("thermal_bc") and dict(init_cfg["thermal_bc"]) != tbc:
        print(f"[prepare_ns] warning: initializer.thermal_bc={init_cfg['thermal_bc']} は無視 (spec.wall_thermal={p.wall_thermal} を使用)")
    res_init = integral_bl(d["wall"], wall_inv, _gam_or_gas(p), p.cp, float(p.spec["Pt"]), float(p.spec["Tt"]),
                           scale, thermal_bc=tbc,
                           theta0_m=init_cfg.get("theta0_m"), x_virtual_m=init_cfg.get("x_virtual_m"),
                           a_crocco=float(init_cfg.get("a_crocco", 1.0)), closure=str(init_cfg.get("closure", "contur")),
                           cf_scale=float(init_cfg.get("cf_scale", 1.0)), n_scale=float(init_cfg.get("n_scale", 1.0)),
                           **({} if rtol is None else {"rtol": float(rtol)}))
    # 積分法の出力も同じ 5 次 P-spline で平滑化 (N(Re) テーブルの折れ目などを壁曲率に持ち込まない)
    from ..metrics.deltastar import smooth_delta_quintic
    f_s, sm_diag = smooth_delta_quintic(res_init["x"], res_init["delta_r"], knot_spacing=2.0, lam=1.0,
                                        positive=(p.wall_thermal["mode"] == "adiabatic"))  # 等温は符号付き
    res_init["delta_r_raw_integral"] = res_init["delta_r"].copy()
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
    return res_init, delta_r_x, init_info


def build_physical_wall(p: Problem, d: dict, scale: float, delta_r_x=None, dstar_x=None, offset: str = "normal",
                        pwu: dict | None = None):
    """物理壁の構築 (`prepare_ns` と寸法の逆算 `deltastar_loop.solve_rt`・`solve_rt_throat` が共有する; plan
    tooling-nozzle-upstream-poly-and-throat-sizing §4.2「同じ壁の構築」)。`PhysicalNozzleWall` に `geometry.pw_ramp`・
    `geometry.pw_upstream` (`_pw_upstream` の解決済みの値) を渡し、`physical_wall_repr: single_bspline` なら 1 本の B-spline に
    作り直す。scale: r_t [m] (相関 δ* の診断に使う。壁の形は r_t 単位)。"""
    from ..geometry.wall_axismach import PhysicalNozzleWall, SingleBSplinePhysicalWall
    pwu = _pw_upstream(p.geometry) if pwu is None else pwu
    # joint 壁の物理壁 (解析経路) の δ_r ランプ区間: geometry.pw_ramp (None = 既定、default_pw_ramp)。ゲート不合格は例外
    pw_ramp = p.geometry.get("pw_ramp")
    wall = PhysicalNozzleWall(d["wall"], d["wall_inv"], float(scale), float(p.spec["Pt"]),
                              float(p.spec["Tt"]), _gam_or_gas(p), p.cp, dstar_x=dstar_x,
                              offset=offset, delta_r_x=delta_r_x,
                              ramp=(None if pw_ramp is None else tuple(float(v) for v in pw_ramp)),
                              upstream=pwu["value"])
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


def prepare_ns(problem_path, run_dir, nsteps=None, ic_from=None,
               dstar_csv=None, dstar_blend=(6.0, 9.0),
               delta_r_csv=None, offset: str = "normal", euler_ref=None,
               omega: float | None = None, prev_run=None, initializer=None,
               cfl_main: float | None = None, implicit_relax: float | None = None) -> dict:
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
    """
    from ..feedback.deltastar import _sutherland
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
        delta_r_x = delta_r_from_table(tbl_r[:, 0], tbl_r[:, 1])
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
        dstar_src = f"delta_r_csv:{delta_r_csv} (radial, core-matched Euler ref {euler_ref}, omega {omega})"
    else:
        dstar_src = ("correlation_hist_v1" if dstar_csv is None else f"{dstar_csv} (blend {dstar_blend})") + f" ({offset})"
    msgs = wall.validate()
    if msgs:
        raise ValueError("物理壁フィルタ不合格: " + "; ".join(msgs))
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
    mp = mesh_params(p, scale, ni=561, nj=97, wall_first_frac=4.5e-5)
    coords, quads, bedges = generate_axisym_mesh(wall, mp)
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
    n = nsteps or int(p.evaluate.get("nStepOuter", 48000))
    out_int = int(p.evaluate.get("outStepInterval", max(n // 6, 1)))
    if n % out_int:
        # forge は outStepInterval の倍数でしか res を書かない → 最終 step の res が残るよう n を割り切る間隔に直す
        out_int = next(n // k for k in (3, 2, 4, 6, 1) if n % k == 0)
    cfl_main = float(cfl_main if cfl_main is not None else p.evaluate.get("cfl_main", 1.0))
    implicit_relax = implicit_relax if implicit_relax is not None else p.evaluate.get("implicit_relax")
    (run_dir / "bcondConfig.yaml").write_text(_bcond_with_species(_bcond(p, euler=False), _tp_species_Y(p)))
    (run_dir / "probe.yaml").write_text(PROBE_STUB)
    # 品質は cell 変換コピーで検査 (品質ツールは node CONNE 非対応)
    (run_dir / "solverConfig.yaml").write_text(
        _apply_gas_to_config(_config_euler(p, n, out_int, 4.0, 1), p, run_dir))
    subprocess.run([str(converter_path()), "nozzle.msh", "nozzle_qc.h5"],
                   cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    q = subprocess.run([sys.executable, str(FORGE_TOOLS / "check_mesh_quality.py"),
                            "nozzle_qc.h5", "--ar-max", str(int(p.mesh.get("ar_max", 1000)))], cwd=run_dir, env=_ENV,
                       capture_output=True, text=True)
    (run_dir / "MESH_QUALITY.txt").write_text(
        "# cell 変換コピーで検査 (品質は primal の性質)\n" + q.stdout + q.stderr)
    (run_dir / "nozzle_qc.h5").unlink()
    if q.returncode != 0:
        raise RuntimeError(f"メッシュ品質 FAIL:\n{q.stdout}")
    # node/SST 変換 (config を先に書く — wall_dist は no-slip 壁で作られる)
    cfg_ns = _apply_gas_to_config(_config_sst_node(p, n, out_int, cfl_main), p, run_dir, viscous=True)
    if implicit_relax is not None:
        # 陰解法の緩和 (cfl 6 + implicitRelax 0.7 が生産推奨: case/45 run_0018)。deltaT ブロックに挿入
        cfg_ns = cfg_ns.replace("blockDPLUR: 1,", f"blockDPLUR: 1, implicitRelax: {float(implicit_relax)},", 1)
    (run_dir / "solverConfig.yaml").write_text(cfg_ns)
    subprocess.run([str(converter_path()), "nozzle.msh", "nozzle.h5"],
                   cwd=run_dir, env=_ENV, check=True, capture_output=True, text=True)
    # 記録だけ (NS の格子は mesh のまま、座標・接続は変えない; 2026-10-07 plan verification-case45-euler-total-enthalpy §4)
    mesh_rec = _mesh_record(p, mp, "mesh", coords.shape[0] // mp.ni, run_dir, coords, quads, bedges)
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
            "offset": wall.offset_mode, "euler_ref": (str(euler_ref) if euler_ref else None),
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
            "nStepOuter": n, "cfl_main": cfl_main, "implicit_relax": implicit_relax, "scale_m": scale,
            "wall_thermal": p.wall_thermal,
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
    info["sizing"] = _sizing_record(sizing, wall, scale)
    info["environment"] = environment_record()     # 照合は同じ環境どうしで (2026-10-07)
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


def run_staged_ns(run_dir, stages: str = "full", ramp=None, ramp_steps: int = 1000) -> int:
    """NS の起動。stages:
    - "full" (既定・run_0030 レシピ): soft (1次 cfl0.5 ni10, 3000 step) → mid (1次 cfl1, 3000) → 本段 (2次 cfl_main)。
    - "none": 本段だけ (収束済み NS 場からの warm start 用)。
    - "ramp": 2 次のまま cfl を `ramp` (例 (1, 2, 3.5)) の順に各 ramp_steps だけ回して本段 cfl_main へ
      (forge に CFL ランプ機能は無いので restart で段階化する。本段の step 数はランプ分を差し引く)。
    各段の最終場を IC に引き継ぐ。
    段ごとの実効設定を `stage_manifest.json` に記録し (`check_convergence.py --segment` 用)、段の残差履歴を
    `residual_history_<tag>.csv` に残す (段の res_* は従来どおり消す)。段の最終 res の必要保存量が非有限・ρ≤0 なら
    次段へ進まず RuntimeError (段終了ゲート `stage_gate`; plan tooling-rerun-conditions §4.9)。"""
    import shutil
    run_dir = Path(run_dir)
    cfg_main = (run_dir / "solverConfig.yaml").read_text()
    bc_path = run_dir / "bcondConfig.yaml"
    bc_text = bc_path.read_text() if bc_path.exists() else ""
    n_main = _cfg_value(cfg_main, "nStepOuter")
    main_want = {k: _cfg_value(cfg_main, k) for k in ("convMethod", "cfl", "cfl_pseudo", "outStepInterval")}
    sm = _stage_manifest_cls()(run_dir)

    def _record(cfg, tag):
        """段の残差履歴を段名つきで残し、manifest に段を足して書く (失敗した段も診断用に残す)。"""
        hist = run_dir / "residual_history.csv"
        if hist.exists():
            shutil.copy(hist, run_dir / f"residual_history_{tag}.csv")
        sm.add(tag, cfg, bc_text, history=f"residual_history_{tag}.csv")
        sm.write()

    def _stage(cfg, nsteps, tag, want):
        cfg = _cfg_set(cfg, {"nStepOuter": nsteps, "outStepInterval": nsteps})
        _cfg_check(cfg, {**want, "nStepOuter": nsteps, "outStepInterval": nsteps}, tag)   # 起動前に実効値を検査
        (run_dir / "solverConfig.yaml").write_text(cfg)
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

    def _pre(cfg, cfl):
        """前段: CFL を cfl に、空間 1 次、nStepInner 5 → 10 (従来どおり; 5 以外は据え置き)。"""
        cfg = _first_order(_cfg_set(cfg, {"cfl": cfl, "cfl_pseudo": cfl}))
        if _cfg_num(_cfg_get(_cfg_load(cfg), "nStepInner")) == 5.0:
            cfg = _cfg_set(cfg, {"nStepInner": 10})
        return cfg

    if stages == "full":
        _stage(_pre(cfg_main, 0.5), 3000, "S1_soft", {"convMethod": 0, "cfl": 0.5, "cfl_pseudo": 0.5})
        _stage(_pre(cfg_main, 1.0), 3000, "S2_mid", {"convMethod": 0, "cfl": 1.0, "cfl_pseudo": 1.0})
    elif stages == "ramp":
        ramp = tuple(ramp or (1.0, 2.0, 3.5))
        for i, c in enumerate(ramp):
            st = _cfg_set(cfg_main, {"cfl": c, "cfl_pseudo": c})
            _stage(st, int(ramp_steps), f"R{i + 1}_cfl{c:g}",
                   {"convMethod": main_want["convMethod"], "cfl": c, "cfl_pseudo": c})
        n_main = max(n_main - int(ramp_steps) * len(ramp), int(ramp_steps))
        # 最終 res が書かれるよう outStepInterval の倍数に丸める (forge は outStepInterval の倍数でしか res を書かない)
        out_int = main_want["outStepInterval"]
        n_main = max((n_main // out_int) * out_int, out_int)
        cfg_main = _cfg_set(cfg_main, {"nStepOuter": n_main})
    elif stages != "none":
        raise ValueError("stages は 'full' / 'none' / 'ramp'")
    _cfg_check(cfg_main, {**main_want, "nStepOuter": n_main}, "main")
    (run_dir / "solverConfig.yaml").write_text(cfg_main)
    (run_dir / "residual_history.csv").unlink(missing_ok=True)
    rc = run_forge(run_dir)
    _record(cfg_main, "main")
    return rc


if __name__ == "__main__":
    sys.exit(main())
