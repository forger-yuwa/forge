"""v2: Euler 帰還ループ (親計画 §4.7 の 2 段目)。

1 パス = forge Euler 収束計算 → 軸 ΔM 評価 → **v1 特性線網の凍結 C⁻ マップ**
(壁点 → その波が着地する軸点) で残差を壁へ帰還 (PM 換算 Δθ_w = ω·Δν) →
壁再構築 (θ_w に平滑化增分を加え r=∫tanθ で再積分, 接合点アンカー) →
同一トポロジ再メッシュ → 前パス場から warm restart。

マスク (§4.7(c)): 壁足が円弧接合近傍に落ちる station (帰還で直せない —
接合波スパイクを含む) と目標終端近傍はランプで重みゼロ。
収束: masked ‖ΔM‖∞ ≤ tol (既定 0.5% Md) or パス上限。

使い方:
  design/.venv-opt/bin/python -m forge_design.feedback.euler_loop \
      case/41.wind_tunnel_design/problem_m4.yaml <work_dir> [--passes 8]
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import LinearNDInterpolator

from ..geometry.moc_kernel import pm_nu
from ..geometry.cminus_map import build_cminus_map
from ..geometry.wall_modef import ModeFWall
from ..meshing.mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d
from ..probdef import load_problem
from .. import evaluate
from ..evaluate.ic import paste_isentropic_ic
from ..evaluate import judge as J
from ..evaluate.runner import FORGE_BUILD, FORGE_TOOLS, PROBE_STUB, _ENV, run_forge
from ..evaluate import runner_wt


def build_cminus_map_cfd(run_dir, wall_pts, scale: float, gamma: float,
                         m_floor: float = 1.05) -> np.ndarray:
    r"""**v3**: そのパスの CFD 最終場の中で C⁻ を壁→軸へ積分してマップを作る。

    注: 「収束場」とは呼ばない — 収束判定は `evaluate.health` が
    `check_convergence.py` の VERDICT として別途記録し、読み替えはしない (A3)。

    v2 の凍結マップ (v1 の MOC 設計場= 完全気体・等エントロピー・単一組成) と違い、
    局所の流れ角 θ と マッハ角 μ を**実場から**読むため、粘性による有効形状の変化・
    組成分布・物性の温度依存が自動的に入る (ユーザ指摘 2026-08-14)。
    ガード: トレースが $M<m\_floor$ 域に入った壁点は棄却 (§4.7(c))。
    戻り: (n,2) [x_w, x_axis] (無次元 r*)。"""
    run_dir = Path(run_dir)
    res = sorted(run_dir.glob("res_[0-9]*.h5"),
                 key=lambda f: int("".join(c for c in f.stem if c.isdigit())))[-1]
    with h5py.File(run_dir / "nozzle.h5") as nz:
        cc = nz["CELLS/centCoords"][:].reshape(-1, 3)
        node_c = (nz["MESH/COORD"][:].reshape(-1, 3) if "MESH/COORD" in nz else None)
    with h5py.File(res) as f:
        Ux, Uy, son = f["VALUE/Ux"][:], f["VALUE/Uy"][:], f["VALUE/sonic"][:]
    # node は軸上に DOF を持つので**真の節点座標**を使い r=0 まで積分できる
    # (`CELLS/centCoords` は node では双対 CV 重心で最小 r が 0 でない)
    is_node = node_c is not None and len(node_c) == len(Ux)
    if is_node:
        cc = node_c
    x, r = cc[:, 0] / scale, cc[:, 1] / scale          # → r* 単位
    M = np.hypot(Ux, Uy) / np.maximum(son, 1e-9)
    th = np.arctan2(Uy, Ux)
    itp_M = LinearNDInterpolator(np.c_[x, r], M)
    itp_th = LinearNDInterpolator(np.c_[x, r], th)
    # **cell モード限定の回避策**: cell 中心は軸上に無いので補間の凸包が r_min で
    # 切れ、軸到達前に必ず NaN になり全点棄却される (2026-08-14 実測)。そこまで
    # 積分して最後の C⁻ 勾配で r=0 へ外挿する (θ→0 で勾配は緩やかなので誤差は小)。
    # **node では外挿不要** — 真の節点座標を使うので凸包が r=0 まで届く (2026-08-17)。
    r_axis = 0.0 if is_node else float(r.min()) * 1.5
    rows = []
    for xw, rw in wall_pts[::4, :2]:
        xx, rr = float(xw), float(rw) * 0.985          # 壁のわずか内側から出発
        ok, slope = True, None
        for _ in range(4000):
            if rr <= r_axis:
                break
            Mv, tv = float(itp_M(xx, rr)), float(itp_th(xx, rr))
            if not (np.isfinite(Mv) and np.isfinite(tv)) or Mv < m_floor:
                ok = False
                break
            mu = np.arcsin(1.0 / max(Mv, 1.0 + 1e-9))
            dr = -min(0.02, rr - r_axis)
            slope = np.tan(tv - mu)
            if abs(slope) < 1e-9:
                ok = False
                break
            xx += dr / slope
            rr += dr
        if ok and slope is not None and rr <= r_axis * 1.01:
            rows.append((float(xw), xx - rr / slope))   # r=0 へ外挿
    return np.asarray(rows)


def junction_taper(x: np.ndarray, x_j: float, width: float,
                   order: int = 5) -> np.ndarray:
    r"""接合点 $x_j$ で $\Delta\theta$ を滑らかに立ち上げるテーパ。

    壁補正を流れ角 $\Delta\theta$ で与えるとき、接合点で
    $\Delta\theta(x_j)=0$ かつ $\Delta\theta'(x_j)=0$ を課せば、$r'=\tan\theta$ より
    **$r'$ と $r''$ の両方が保たれる** (円弧との C2 接続を壊さない)。再構築は
    $r_{new}'=\tan(\theta+\Delta\theta)$、
    $r_{new}''=\sec^2(\theta+\Delta\theta)\,(\theta'+\Delta\theta')$
    なので、この 2 条件から保存は**構成的に厳密**。

    $x\le x_j$ では 0 (円弧区間は完全凍結)、$x_j+width$ 以降で 1。
    `order=3` は smoothstep $t^2(3-2t)$ ($f=f'=0$ @0)、`order=5` は quintic
    $t^3(10-15t+6t^2)$ ($f=f'=f''=0$ @0) で曲率の立ち上がりも滑らか (既定)。

    幅マスクで広範囲を殺すのではなく**接合点でのみ拘束**するのが要点
    (2026-08-14 の指摘): 広幅マスクは接合波の帯ごと帰還を無効化してしまう。
    """
    t = np.clip((np.asarray(x, dtype=float) - x_j) / max(width, 1e-12), 0.0, 1.0)
    if order == 3:
        return t * t * (3.0 - 2.0 * t)
    return t ** 3 * (10.0 - 15.0 * t + 6.0 * t * t)


def _rebuild_wall(wall_pts: np.ndarray, wall_th: np.ndarray, dtheta_fn,
                  x_anchor: float, taper_width: float = 0.0,
                  taper_order: int = 5) -> tuple:
    r"""壁角 $\theta$ を**主変数として持ち回り**、$r=\int\tan\theta\,dx$ で再構築する。

    **重要 (2026-08-14)**: 以前は $r$ を `np.gradient` で数値微分して $\theta$ を作り
    再積分していたが、この往復は $r$ に $O(h^2)$ の形状誤差を入れ、2 階微分では
    $O(1)$ の誤差になる — 実測で**接合点の $r''$ が 0.50→0.25 と半減し、しかも
    解像度に依らなかった**。$\theta$ を配列として持ち回れば微分が不要になり、
    $r'(x_j)=\tan\theta[0]$、$r''(x_j)=\sec^2\theta\,\theta'(x_j)$ が入力どおり保たれる。

    `taper_width>0` のとき接合点で $\Delta\theta=\Delta\theta'=0$ (テーパ) となり、
    円弧との $r,r',r''$ 連続が**構成的に**保たれる。戻り: (wall_pts_new, th_new)。
    """
    x, r = wall_pts[:, 0], wall_pts[:, 1]
    th = np.asarray(wall_th, dtype=float)
    dth = np.asarray(dtheta_fn(x), dtype=float)
    if taper_width > 0.0:
        dth = dth * junction_taper(x, x_anchor, taper_width, order=taper_order)
    th_new = th + dth
    r_new = np.empty_like(r)
    i0 = int(np.argmin(np.abs(x - x_anchor)))
    r_new[i0] = r[i0]
    tn = np.tan(th_new)
    for i in range(i0 + 1, len(x)):
        r_new[i] = r_new[i - 1] + 0.5 * (tn[i - 1] + tn[i]) * (x[i] - x[i - 1])
    for i in range(i0 - 1, -1, -1):
        r_new[i] = r_new[i + 1] - 0.5 * (tn[i] + tn[i + 1]) * (x[i + 1] - x[i])
    return np.c_[x, r_new], th_new


CHUNK = 6000          # warm 継続チャンク長 [step] (= axis_M_series_steady の n_ref)
MAX_WARM_CHUNKS = 6   # warm パスの継続上限 (実測: 壁更新後の遅い波の定常化に
                      # 24000 step で trend 0.008 vs tol 0.005 の僅差不足が
                      # 2 run で再現 [run_0023 pass6 / run_0024 pass8] → 36000 まで許可)


def _steps_of(f: Path) -> int:
    return int("".join(c for c in f.stem if c.isdigit()))


def _chunk_dirs(rd: Path) -> list:
    """パスのチャンク列 [rd, rd/cont_1, ...] (存在するものだけ・番号順)。"""
    return [rd] + sorted(rd.glob("cont_[0-9]*"),
                         key=lambda d: int(d.name.split("_")[1]))


def _snaps_of(rd: Path) -> list:
    """全チャンクを通しの step に直したスナップショット列 [(step, res, mesh)]。

    各チャンクの res_0 は直前チャンク最終場の貼り直し (同一メッシュの index コピー = ビット一致; 2026-10-15 まで
    最近傍補間) なので通し step は「オフセット + チャンク内 step」。判定側は dstep>0 の対のみ
    使うため res_0 の重複は無害。"""
    out, off = [], 0
    for d in _chunk_dirs(rd):
        res = sorted(d.glob("res_[0-9]*.h5"), key=_steps_of)
        out += [(off + _steps_of(f), f, d / "nozzle.h5") for f in res]
        if res:
            off += _steps_of(res[-1])
    return out


def start_cont_chunk(rd: Path, j: int, last_res: Path, cfg_chunk: str) -> Path:
    """継続チャンク `rd/cont_<j>` を準備する (forge は起動しない; 呼び出し側が宣言を書いて起動する)。パスの格子 (`rd/nozzle.h5`;
    メッシュは全チャンク同一) と bcond・probe を写し、config を書き、直前のチャンクの最終場 `last_res` を段終了ゲート
    (必要保存量・有限・ρ>0; res を作ったチャンクの config で) に通してから `restart_field.py` (保存量の index コピー、ビット一致の検査、
    `wall_dist` は宛先のまま) で初期場に写す (`runner_wt._staged_handoff`)。ゲート・`restart_field.py` に落ちたら RuntimeError
    (チャンクを起動しない)。2026-10-15 まで同じ格子へ `interp_field.py` (cross-mesh 用、原始量から保存量を組み直す) を使っていた
    (AGENTS.md「メッシュ変更後の restart」の規則違反; plan tooling-generator-defaults-and-launch-gate §4.5・§5.1 #19)。戻り: チャンクの dir。"""
    import shutil
    cd = rd / f"cont_{j}"
    cd.mkdir(exist_ok=True)
    shutil.copy(rd / "nozzle.h5", cd / "nozzle.h5")  # メッシュは全チャンク同一
    (cd / "solverConfig.yaml").write_text(cfg_chunk)
    for cf in ("bcondConfig.yaml", "probe.yaml"):
        shutil.copy(rd / cf, cd / cf)
    last_res = Path(last_res)
    runner_wt._staged_handoff(cd, last_res, (last_res.parent / "solverConfig.yaml").read_text(), last_res.parent.name, "nozzle.h5")
    return cd


def pass_config(p, n: int, use_node: bool | None = None) -> tuple:
    """パス (と継続チャンク) の solverConfig: node は `runner_wt._config_euler_node`、cell は `_config_euler` (n step、出力 n/4)。
    機種のプロファイル (wt_feedback = node) に合えば明示すべきキーを旧い実効値で足す (plan tooling-generator-defaults-and-launch-gate
    §5.1 #11、`config_builder.finish_config`; cell は書き写したプロファイルが無いので変えない)。戻り (config, 仕上げの記録)。
    書き写しの台本 (recipes/regen_profiles.py) も使う。"""
    from ..evaluate import config_builder as CB
    if use_node is None:
        use_node = str(p.mesh.get("discretization", "node")).lower() == "node"
    oi = max(n // 4, 1)
    cfg = runner_wt._config_euler_node(p, n, oi, 4.0, 1) if use_node else runner_wt._config_euler(p, n, oi, 4.0, 1)
    return CB.finish_config(cfg, runner_wt._bcond(p, euler=True), "wt_feedback", tag="main")


def run_loop(problem_path, work_dir, n_pass: int = 8, omega: float = 0.4,
             tol_rel: float = 0.005, clip_deg: float = 0.5, nsteps=None,
             map_mode: str = "frozen", gate_mode: str = "formal",
             mask_mode: str = "wide", taper_width: float = 0.3,
             drop: float = 2.0) -> dict:
    """map_mode: "frozen" = v2 (v1 MOC 設計場の凍結マップ) /
    "cfd" = **v3** (毎パスその回の CFD 収束場から C⁻ を引き直す)。

    **実行方針 (2026-08-15, レビュー反映)**: pass 1 は準 1D IC から cold で
    `nsteps` (既定 24000) 回し、正式ゲート (PASS + 準定常 + 軸 M 定常) の
    アンカーとする。pass 2 以降は**前パス最終場を cross-mesh restart**
    (`interp_field.py`, AGENTS のメッシュ変更後 restart 必須ルール) して
    CHUNK step ずつ回し、**軸 M 定常 (レート正規化) になり次第停止**する
    (適応停止)。壁は 1 パスに数十 µm しか動かないので過渡は小さく、
    cold 24000 step/パスの大半は無駄だった (ユーザ指摘)。"""
    p = load_problem(problem_path)
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    # 同一 work dir を formal / exploratory で使い回さない (結果の由来が混ざるため)
    man_f = work / "run_manifest.json"
    if man_f.exists():
        man = json.loads(man_f.read_text())
        if man.get("gate_mode") != gate_mode:
            raise RuntimeError(
                f"work dir {work} は gate_mode='{man.get('gate_mode')}' で作成済み。"
                f"'{gate_mode}' との混用は不可 — 別の run ディレクトリを使うこと")
    else:
        man_f.write_text(json.dumps({"gate_mode": gate_mode,
                                     "problem": str(problem_path),
                                     "map_mode": map_mode}, indent=1))
    d = runner_wt.design_chain(p, viscous=False)
    Md, xd, x0 = d["Md"], d["xd"], d["x0"]
    scale = float(p.spec["r_throat"])
    cmap = build_cminus_map(d["pts"], d["wall_inv"], p.gamma)
    x_w0 = float(d["wall_inv"][0, 0])
    # **x_reach はループ中 固定** (B10): design_chain が config の
    # `geometry.x_reach_cfd` (node 基準 run 由来) を返すのでパス毎に再導出しない。
    # 旧 B8 経路では逆 MOC マップ由来の値が入る。
    x_reach = float(d.get("x_reach", x_w0))
    x_reach_mode = d.get("x_reach_mode", "moc")
    # node/cell は problem YAML の `mesh.discretization` で選ぶ (既定 node — B10)
    use_node = str(p.mesh.get("discretization", "node")).lower() == "node"
    from ..evaluate import config_builder as CB

    def _declare(d, fin):
        """パス・継続チャンクの段の宣言 RECIPE.json (forge の起動の直前; plan tooling-generator-defaults-and-launch-gate §5.1 #11)。"""
        CB.declare_stage(d, "main", "wt_feedback", finished=fin, problem=problem_path, generator="feedback.euler_loop.run_loop")

    def make_ax2w(cm):
        r"""軸→壁の逆引き。`np.interp` は範囲外を端値へ**黙ってクランプ**するため、
        「壁が届かない (物理)」と「マップの被覆が単に足りない (測定限界)」を
        区別できない (2026-08-15 レビュー指摘)。呼び出し側が `.lo`/`.hi` で
        被覆範囲を診断できるよう関数に持たせる — 挙動 (クランプ) 自体は変えず、
        `mask_out_of_coverage_frac` として記録するに留める。"""
        o = np.argsort(cm[:, 1])
        f = lambda xi: np.interp(xi, cm[o, 1], cm[o, 0])
        f.lo, f.hi = float(cm[o, 1].min()), float(cm[o, 1].max())
        return f

    ax2w = make_ax2w(cmap)   # v3 では各パスの CFD 場で毎回置き換える

    # 評価目標: 到達不能域 [x0,x_reach) は実測こぶ曲線・以降は設計 Bézier
    # (B7 最終形 — 設計用 d["target"] と分離。無ければ後方互換で design を使う)
    tgt_fn = d.get("target_eval", d["target"])
    # === C3 (2026-08-16): 評価窓をリップ−マージンまで延長 ===
    # 旧窓は x_d−0.2 で切れており、その C⁻ 足より下流の壁 (整流区間〜リップ側)
    # は一度も帰還を受けなかった。x_d 以降の目標は M_d 一定 (target が返す)。
    # メッシュはリップまでしか無いので軸点もそこまで — C⁻ 足がメッシュ内に
    # 着地する壁 (x_w ≲ 13) までが帰還可能で、以遠は最終 Δθ の定数外挿。
    x_lip = float(d["wall_inv"][-1, 0])
    x_eval_end = x_lip - 0.5
    # **評価窓の始点** (B10): target が定義されるのは x >= x_reach のみなので、
    # そこから始める。上流は target 誤差ではなく物理診断として別途記録する。
    x_eval_lo = (x_reach if x_reach_mode == "cfd" else x0)
    tgt = np.array([[x, tgt_fn(x)] for x in np.linspace(x_eval_lo, x_eval_end, 800)])
    wall_pts = d["wall_inv"][:, :2].copy()
    wall_th = d["wall_inv"][:, 2].copy()   # θ を主変数として持ち回る (数値微分を避ける)
    # **初期壁を仮壁で置き換える** (B10 手順 6: 「現在の仮壁を初期値として帰還
    # ループで J→リップを更新する」)。新 target から逆設計し直した壁で始めると、
    # 既存 run が積んだ改良を捨てることになる (実測: 逆設計壁は仮壁と最大 170 µm
    # 違い、うねり 0.0214→0.0279・出口 |θ|max 0.0150→0.0474° と物理が後退した)。
    # θ は主変数なので壁点列から復元する (r の数値微分を避けるため、設計 θ を
    # 仮壁の勾配で置き換えるのではなく、仮壁の x 上で θ=arctan(dr/dx) を
    # スプライン微分で作る — ModeFWall と同じクランプ付き 5 次表現を使う)。
    init_wall = p.geometry.get("initial_wall_csv")
    if init_wall:
        w0 = np.loadtxt(init_wall, delimiter=",", skiprows=1)
        mw0 = ModeFWall(w0, R=d["R"],
                        r_inlet=float(p.geometry.get("r_inlet", 2.5)),
                        L_pipe=float(p.geometry.get("L_pipe", 0.5)),
                        L_contract=float(p.geometry.get("L_contract", 3.0)),
                        phi_u=np.deg2rad(float(p.geometry.get("phi_u_deg", 25.0))),
                        blend_len=0.0)
        msgs0 = mw0.validate()
        if msgs0:
            raise ValueError("initial_wall_csv の壁が不合格: " + "; ".join(msgs0))
        wall_pts = w0[:, :2].copy()
        wall_th = np.arctan(mw0.drdx(w0[:, 0]))
        print(f"[init] 仮壁を初期値に採用: {init_wall} ({len(wall_pts)} 点)", flush=True)
    prev_res = None
    hist = []
    # trust-region 型の決定的下り探索: 悪化したら最良壁に巻き戻し ω 半減
    best = {"wall": wall_pts.copy(), "th": wall_th.copy(), "dM_inf": np.inf,
            "resid": None, "ax2w": ax2w}   # C1: マップも最良状態と一緒に保持し reject 時に戻す
    for k in range(1, n_pass + 1):
        rd = work / f"pass_{k:02d}"
        if not (rd / "metrics_pass.json").exists():
            rd.mkdir(parents=True, exist_ok=True)
            wall = ModeFWall(wall_pts, R=d["R"],
                             r_inlet=float(p.geometry.get("r_inlet", 2.5)),
                             L_pipe=float(p.geometry.get("L_pipe", 0.5)),
                             L_contract=float(p.geometry.get("L_contract", 3.0)),
                             phi_u=np.deg2rad(float(p.geometry.get("phi_u_deg", 25.0))),
                     blend_len=float(p.geometry.get("blend_len", 0.0)))
            # design_chain の初回壁と同じ検査をパス毎の再構築壁にも課す
            # (以前は呼ばれておらず、帰還が壁を壊しても気づけなかった — レビュー指摘)
            wmsgs = wall.validate()
            if wmsgs:
                raise ValueError(f"pass {k} 壁フィルタ不合格: " + "; ".join(wmsgs))
            mp = Mesh2DParams(ni=int(p.mesh.get("ni", 321)), nj=int(p.mesh.get("nj", 65)),
                              wall_first_frac=float(p.mesh.get("wall_first_frac", 5.0e-3)),
                              throat_refine=float(p.mesh.get("throat_refine", 3.0)),
                              scale=scale)
            coords, quads, bedges = generate_axisym_mesh(wall, mp)
            write_msh41_2d(rd / "nozzle.msh", coords, quads, bedges)
            warm = (k > 1 and prev_res is not None)
            # pass 1 (cold アンカー): 24000 step 固定 — 残差は step ~371 で全列
            # 2 桁到達後プラトーだが、軸 M の定常化に 24000 要る (実測 run_0008:
            # 射影ドリフト 12k=0.0079 → 24k=0.0036)。
            # pass ≥2 (warm): CHUNK ずつの適応継続なのでここは 1 チャンク分。
            n = CHUNK if warm else max(nsteps or 24000, 12000)
            cfg_pass, fin_pass = pass_config(p, n, use_node)
            (rd / "solverConfig.yaml").write_text(cfg_pass)
            (rd / "bcondConfig.yaml").write_text(runner_wt._bcond(p, euler=True))
            (rd / "probe.yaml").write_text(PROBE_STUB)
            subprocess.run([str(FORGE_BUILD / "convertGmshToForge"), "nozzle.msh", "nozzle.h5"],
                           cwd=rd, env=_ENV, check=True, capture_output=True, text=True)
            if warm:
                # メッシュ変更後の cross-mesh restart (AGENTS 必須ルール):
                # 前パス収束場を最近傍で移植 — uniform/準1D IC からやり直さない
                subprocess.run([sys.executable, str(FORGE_TOOLS / "interp_field.py"),
                                str(prev_res), str(rd / "nozzle.h5")],
                               env=_ENV, check=True, capture_output=True, text=True)
            else:
                paste_isentropic_ic(rd / "nozzle.h5", wall, scale,
                                    float(p.spec["Pt"]), float(p.spec["Tt"]),
                                    p.gamma, p.cp)
            _declare(rd, fin_pass)
            # 計算の後の判定 (plan tooling-convergence-and-wall-resolution-gates §4.6・§5.1 #12): 仕様は起動の前に、総合は起動の後に 1 度だけ
            J.write_spec(rd, J.spec_snapshot("feedback.euler_loop.run_loop (パス)"))
            with J.Staged(rd, "feedback.euler_loop.run_loop (パス)"):
                run_forge(rd)
            # === 適応停止: 軸 M 判定が合格 (STEADY / OSCILLATING+平均精度) に
            # なるまで CHUNK 継続 === cold パスにも適用する (正式 PASS は第 1
            # チャンクの過渡で判定するので継続しても壊れない — conv/eff 分離)。
            # DRIFTING (系統移動) と OSCILLATING で平均がまだ粗い場合に継続 —
            # 後者は追いサンプルで平均の標準誤差が締まる。
            from ..evaluate.health import axis_M_series_steady, axis_ok
            for j in range(1, MAX_WARM_CHUNKS):
                ax = axis_M_series_steady(
                    rd, scale, (x0 + 0.3) * scale, (x_eval_end - 0.2) * scale,
                    snaps=_snaps_of(rd))
                if axis_ok(ax):
                    break
                last = _snaps_of(rd)[-1]   # 直近チャンクの最終場 (継続 IC)
                # 継続チャンクは常に CHUNK step (cold 第 1 チャンクが長くても)
                cfg_chunk, fin_chunk = pass_config(p, CHUNK, use_node)
                # 同じ格子の継続: 段終了ゲート + restart_field.py (旧: interp_field.py; plan §5.1 #19)
                cd = start_cont_chunk(rd, j, last[1], cfg_chunk)
                _declare(cd, fin_chunk)
                # 継続チャンクは直前のチャンクと同じ数値設定なら収束をその床から判定する (`--from-floor`)
                prev_d = Path(last[1]).parent
                same = J.same_numerics((prev_d / "solverConfig.yaml").read_text(), (cd / "solverConfig.yaml").read_text(),
                                       (prev_d / "bcondConfig.yaml").read_text(), (cd / "bcondConfig.yaml").read_text())
                J.write_spec(cd, J.spec_snapshot("feedback.euler_loop.run_loop (継続チャンク)",
                                                 from_floor=(os.path.relpath(prev_d, cd) if same else None)))
                with J.Staged(cd, "feedback.euler_loop.run_loop (継続チャンク)"):
                    run_forge(cd)
        eff = _chunk_dirs(rd)[-1]   # 最終チャンク = 評価に使う場
        res = sorted(eff.glob("res_[0-9]*.h5"), key=_steps_of)
        prev_res = res[-1]
        # === A3 ゲートは最終場を得た**直後**に判定する ===
        # 不合格ならマップ生成・軸残差評価・正式 metrics の**いずれにも進まない**
        # (未収束/非定常の場から特性線マップや壁を作らないため — レビュー指摘)。
        # 軸 M 定常 (レート正規化) を axis_window で**実接続** (2026-08-15 —
        # 以前は渡し忘れで一度も実行されていなかった: レビュー指摘)。
        # warm パスは低下桁数要件なし (gate docstring 参照, cold アンカーが前提)。
        from ..evaluate.health import gate
        hl = gate(rd, M_design=Md, x_d=xd * scale, mode=gate_mode, drop=drop,
                  axis_window=((x0 + 0.3) * scale, (x_eval_end - 0.2) * scale), scale=scale,
                  warm=(k > 1), snaps=_snaps_of(rd), eff_dir=eff)
        if not hl["ok"]:
            hist.append({"pass": k, "gate_ok": False,
                         "gate_reasons": hl["gate_reasons"],
                         "convergence_verdict": hl["convergence_verdict"],
                     "convergence_drop_threshold": hl.get("convergence_drop_threshold"),
                         "quasisteady_verdict": hl["quasisteady_verdict"],
                         "eps_series_verdict": hl.get("eps_steady", {}).get("verdict"),
                         "note": "ゲート不合格のため残差評価・マップ生成を行わず中止"})
            (rd / "metrics_pass.json").write_text(json.dumps(hist[-1]))
            print(f"[pass {k}] **中止**: ゲート不合格 ({'; '.join(hl['gate_reasons'])}) "
                  f"— 未収束/非定常の場から壁もマップも作らない", flush=True)
            break
        if map_mode == "cfd":
            cm = build_cminus_map_cfd(eff, wall_pts, scale, p.gamma)
            if len(cm) >= 10:
                ax2w = make_ax2w(cm)
                np.savetxt(rd / "cminus_map.csv", cm, delimiter=",",
                           header="x_wall,x_axis_landing", comments="")
            else:
                print(f"[pass {k}] CFD マップ生成失敗 ({len(cm)} 点) — 前回の表を継続",
                      flush=True)
        # 軸 ΔM 評価 — **末尾スナップショットの平均**を使う (2026-08-15)。
        # 中域 x/rt≈9 に片振幅 ~0.004 の弱い音響リミットサイクルがあり (run_0017
        # pass_01 実測)、単一スナップだと位相次第で ±0.004 の測定ノイズが乗る。
        # STEADY 場では平均は恒等なので常時平均で害はない (AGENTS:
        # OSCILLATING は平均±振幅で扱う、の実体化)。窓はゲートと同じ末尾 5 個 —
        # ゲートの mean_se がこの平均の精度をそのまま保証する。
        from ..metrics.extract import axis_mach
        sn_tail = [s for s in _snaps_of(rd) if s[0] > 0][-5:]
        x_ach = None
        M_acc = []
        for _st, _rf, _mh in sn_tail:
            _x, _M = axis_mach(_mh, _rf, axis_band=0.09 * scale)
            x_ach = _x
            M_acc.append(_M)
        M_ach = np.mean(M_acc, axis=0)
        xi = np.linspace(max(x_eval_lo, x0 + 0.3), x_eval_end, 440)  # C3+B10
        Mt = np.interp(xi, tgt[:, 0], tgt[:, 1])
        Ma = np.interp(xi * scale, x_ach, M_ach)
        dM = Mt - Ma
        # **マスク非依存の共通指標** (2026-08-15 レビュー指摘): dM_inf_masked は
        # wide/narrow で評価集合が違い A/B 比較に使えない。固定窓の最大 |ΔM| を
        # 全構成共通で記録する (接合帯 = 接合波スパイク含む / 中域 = 帰還の主対象)。
        eA = np.abs(dM)
        common = {"dM_inf_common": float(eA.max()),
                  "dM_inf_junction_band": float(eA[xi <= 2.5].max()),
                  "dM_inf_mid_band": float(eA[(xi >= 4.0) & (xi <= 12.0)].max()),
                  # C3: x_d 以降の下流帯 (目標 = M_d 一定、旧窓では無評価だった)
                  "dM_inf_downstream_band": float(eA[xi >= xd].max())
                                            if np.any(xi >= xd) else None}
        # マスク/重みランプ: 壁足が接合近傍 or 目標終端近傍
        xw = ax2w(xi)
        # **被覆外クランプの明示記録** (2026-08-15 レビュー指摘): マスクの挙動
        # (w=0) は変えないが、「物理的に到達不能」と「マップの被覆がそもそも
        # 無い (測定限界)」を区別できるよう、被覆外だった軸ステーション数を
        # 別途記録する。実測 (run_0018): x/rt<1.62 はマップ被覆外で、そこは
        # 独立に実 CFD 場を後退トレースして円弧着地 (到達不能) と確認済み。
        out_of_cov = (xi < ax2w.lo) | (xi > ax2w.hi)
        if mask_mode == "narrow":
            # **凍結された円弧に着地する station のみ**マスクする (本来の §4.7(c))。
            # 設計壁は可動なので殺さない。接合の C2 保存は Δθ のテーパで担保する。
            w = (xw > x_w0).astype(float) * np.clip((x_eval_end - 0.5 - xi) / 1.0, 0.0, 1.0)
        else:   # "wide" (従来): 接合から 0.3 + ランプ 1.0 r* まで設計壁も殺す
            w = np.clip((xw - (x_w0 + 0.3)) / 1.0, 0.0, 1.0) * np.clip((x_eval_end - 0.5 - xi) / 1.0, 0.0, 1.0)
        w *= (Mt > 1.15)
        # **被覆外は w=0 を強制する** (2026-08-15): クランプされた xw は最初の
        # マップ点 (> x_w0) を返すため、到達不能な境界ステーションが w=1 に化け、
        # 直せない境界スパイク (+0.058) が masked 指標を頭打ちにし、かつテーパ帯の
        # 壁へ誤帰還されて壁を歪めていた (run_0023 で実測)。被覆外 = 壁足が不明 =
        # 帰還不能、なので w=0 が正しい。
        w *= (~out_of_cov).astype(float)
        # **x_reach 未満も w=0** (2026-08-15, run_0024 で発見): x_reach は
        # ブートストラップ設計のマップ由来、ax2w は最終設計のマップ由来で被覆が
        # 僅かに違い、x∈[map.lo, x_reach) の隙間ステーションが「到達可能」判定に
        # なる。だがそこの評価目標は凍結アンカー実測 (設計意図でない) で、責任壁は
        # テーパ帯 (係数 ~0) — 7 パス押しても 0.058 のまま動かず masked を頭打ちに
        # した。評価目標の分割点とマスクを一致させるのが整合的。
        w *= (xi >= x_reach).astype(float)
        mask_diag = {"mask_out_of_coverage_frac": float(np.mean(out_of_cov)),
                     "mask_out_of_coverage_and_zero_frac": float(np.mean(out_of_cov & (w < 0.5)))}
        dM_inf = float(np.max(np.abs(dM * (w > 0.5))))
        rejected = bool(dM_inf > best["dM_inf"] * 1.02)
        # 出口一様性は**診断として**毎パス記録する (安価)。役割の区別 (2026-08-14):
        #   本ループの目的 = ΔM→0 = 「宣言した目標 M_c(x) を実際に達成する壁を作る」。
        #     これは出口品質のためでなく**評価の忠実性**のため — dv が「達成された
        #     目標分布」を意味して初めて dv→メトリクス写像が決定的になり、
        #     サロゲートが劣化しない (親計画 §4.7(c))。
        #   ε_M/ε_θ = **最適化ループ (§5 ステップ 7) の目的** — Bézier CP・L・R を
        #     振って最小化する量であり、1 評価内で追うものではない。したがって
        #     ΔM と ε が単調に一致しなくても (実測: pass9 ε_M=0.020% < pass13 0.039%)
        #     本ループの欠陥ではない。
        from ..metrics.extract import exit_uniformity
        try:
            eu = exit_uniformity(eff / "nozzle.h5", res[-1], Md,
                                 x_d=xd * scale, core_mode="traced")
        except Exception as e:  # 測定面が取れない等
            eu = {"error": f"{type(e).__name__}: {e}"}
        hist.append({"pass": k, "dM_inf_masked": dM_inf, "omega": omega,
                     "rejected": rejected, **common, **mask_diag,
                     "n_chunks": len(_chunk_dirs(rd)),
                     "total_steps": _snaps_of(rd)[-1][0],
                     "axis_M_verdict": hl.get("axis_M_steady", {}).get("verdict"),
                     "axis_M_trend": hl.get("axis_M_steady", {}).get("trend_dM_per_nref"),
                     "axis_M_amplitude": hl.get("axis_M_steady", {}).get("amplitude"),
                     "dM_rms": float(np.sqrt(np.mean((dM * w) ** 2))),
                     "eps_M_rms": eu.get("eps_M_rms"),
                     "eps_theta_max_deg": eu.get("eps_theta_max_deg"),
                     "convergence_verdict": hl["convergence_verdict"],
                     "convergence_drop_threshold": hl.get("convergence_drop_threshold"),
                     "quasisteady_verdict": hl["quasisteady_verdict"],
                     "eps_series_verdict": hl.get("eps_steady", {}).get("verdict"),
                     "gate_ok": hl["ok"], "gate_reasons": hl["gate_reasons"]})
        print(f"[pass {k}] masked ‖ΔM‖∞ = {dM_inf:.4f} ({dM_inf/Md*100:.2f}% Md)"
              f"  共通: 全域 {common['dM_inf_common']:.4f} / 接合帯 "
              f"{common['dM_inf_junction_band']:.4f} / 中域 {common['dM_inf_mid_band']:.4f}"
              f"  [{hl['convergence_verdict']}/axis:{hl.get('axis_M_steady',{}).get('verdict','n/a')}"
              f"/{_snaps_of(rd)[-1][0]}step]"
              f"{' [reject→ω/2]' if rejected else ''}", flush=True)
        np.savetxt(rd / "achieved_vs_target.csv", np.c_[xi, Mt, Ma, w],
                   delimiter=",", header="x_rt,M_target,M_achieved,weight", comments="")
        (rd / "metrics_pass.json").write_text(json.dumps(hist[-1]))
        if rejected:
            omega *= 0.5           # 過ゲイン → 最良壁から小さい ω でやり直し
            xi_c, Mt_c, Ma_c, w_c = best["resid"]
            ax2w = best["ax2w"]    # C1: 棄却場のマップで最良場の残差を戻さない
        else:
            best = {"wall": wall_pts.copy(), "th": wall_th.copy(), "dM_inf": dM_inf,
                    "resid": (xi, Mt, Ma, w), "ax2w": ax2w}
            if dM_inf <= tol_rel * Md:
                break
            xi_c, Mt_c, Ma_c, w_c = xi, Mt, Ma, w
        # 帰還: Δν → Δθ_w (壁座標へ写像し平滑化・クリップ)
        dnu = np.array([pm_nu(max(float(m1), 1.001)) - pm_nu(max(float(m2), 1.001))
                        for m1, m2 in zip(Mt_c, Ma_c)])
        dth_i = np.clip(omega * dnu * w_c, -np.deg2rad(clip_deg), np.deg2rad(clip_deg))
        xw_s, o2 = np.unique(ax2w(xi_c), return_index=True)
        dth_w = dth_i[o2]
        for _ in range(3):
            dth_w = np.convolve(dth_w, np.ones(9) / 9.0, mode="same")

        def dtheta_fn(x):
            return np.interp(x, xw_s, dth_w, left=0.0, right=dth_w[-1])

        wall_pts, wall_th = _rebuild_wall(
            best["wall"], best["th"], dtheta_fn, x_anchor=x_w0,
            taper_width=(taper_width if mask_mode == "narrow" else 0.0))
        np.savetxt(rd / "wall_next.csv", wall_pts, delimiter=",",
                   header="x_rt,r_rt", comments="")
    last = hist[-1] if hist else {}
    # 収束の主張には ΔM 条件だけでなく**ゲート合格**を必須とする (warm パスの
    # ゲートは「DIVERGED でない + 軸 M 定常」— gate() docstring 参照。cold の
    # PASS アンカーは pass 1 が担う)
    converged = bool(
        last.get("gate_ok") is True
        and last.get("dM_inf_masked") is not None
        and last["dM_inf_masked"] <= tol_rel * Md)
    out = {"history": hist, "passes": len(hist), "converged": converged,
           "converged_criteria": "gate_ok (cold: PASS+準定常+軸M定常 / warm: "
                                 "非DIVERGED+軸M定常) AND dM_inf_masked<=tol",
           "gate_mode": gate_mode, "tol": tol_rel * Md}
    (work / "loop_summary.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1), flush=True)
    return out


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="v2 Euler 帰還ループ")
    ap.add_argument("problem")
    ap.add_argument("work_dir")
    ap.add_argument("--passes", type=int, default=8)
    ap.add_argument("--omega", type=float, default=0.4)
    ap.add_argument("--drop", type=float, default=2.0,
                    help="check_convergence の要求低下桁数 (既定 2 + 準定常)")
    ap.add_argument("--mask-mode", choices=("wide", "narrow"), default="wide",
                    help="wide = 従来 (接合から 1.3 r* を無効化) / narrow = 円弧のみ + Δθ テーパ")
    ap.add_argument("--taper-width", type=float, default=0.3,
                    help="narrow 時に Δθ=O((x-xj)^2) を立ち上げる幅 [r*]")
    ap.add_argument("--gate-mode", choices=("formal", "exploratory"), default="formal",
                    help="formal = PASS 以外は停止 (正式) / exploratory = 発散以外は継続")
    ap.add_argument("--map-mode", choices=("frozen", "cfd"), default="frozen",
                    help="frozen = v2 (v1 MOC 場の凍結マップ) / cfd = v3 (毎パス実場)")
    ap.add_argument("--steps", type=int, default=None)
    a = ap.parse_args(argv)
    run_loop(a.problem, a.work_dir, n_pass=a.passes, omega=a.omega, nsteps=a.steps,
             map_mode=a.map_mode, gate_mode=a.gate_mode,
             mask_mode=a.mask_mode, taper_width=a.taper_width, drop=a.drop)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
