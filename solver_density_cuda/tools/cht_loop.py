#!/usr/bin/env python3
r"""CHT 外部弱連成ループ (Phase 1): forge の壁熱流束 → 固体シェル → `wallProfile` → forge、を反復する。

仕様は methods/boundary.md「共役熱伝達 (CHT)」、設計判断は
plans/active/boundary-conjugate-heat-transfer.md。固体側は `solid_shell.py`。

1 反復 = **forge 1 回** (CFD 評価は高価なので line search はしない)。
更新は固定点を保存する形 $(A_s+D_f)T^{k+1}=b_s+Q_f(T^k)+D_fT^k$ + Anderson 加速で、
受理は固定重みのメリット関数 $\Phi$、棄却したら最後に良かった状態へ退避して $D_f$ を倍にする。

ディレクトリ構成 (1 ループ = 1 run ディレクトリ):
    <run_dir>/it_000/ … 各反復の完全な forge run (mesh.h5, config, res_*, 壁ダンプ)
    <run_dir>/cht_history.csv … 反復ごとの T_w 統計・残差・Q 合計 (収束判定の一次情報)

usage:
  python3 solver_density_cuda/tools/cht_loop.py <run_dir> \
      --template <template_dir> --forge <forge binary> --solid solid.json \
      --phys-id 4 [--phys-name wall] [--steps 2000] [--max-iter 20] [--flux q_eff]

<template_dir> には mesh.h5 / solverConfig.yaml / bcondConfig.yaml (+ probe.yaml) を置く。
対象壁の bcond には `ints: {wallProfile: 1}` が要る (無ければ本スクリプトが落とす)。
solverConfig の `output` には `interfaceDiag: 1` が要る。
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from solid_shell import ShellOperator, SolidModel   # noqa: E402
from solid_fem2d import Fem2DOperator                # noqa: E402

XDMF_NNODE = {2: 2, 4: 3, 5: 4}     # Polyline / Triangle / Quadrilateral


# ------------------------------------------------------------------ wall dump
def read_wall_dump(path: Path):
    """壁ダンプから (coords, faces, values) を読む。node 可視化 (Center='Node') 前提。"""
    with h5py.File(path, "r") as f:
        coords = np.array(f["MESH/COORD"]).reshape(-1, 3)
        conne = np.array(f["MESH/CONNE"]).astype(np.int64)
        vals = {k: np.array(f["VALUE/" + k]) for k in f["VALUE"].keys()}
    faces, i = [], 0
    while i < conne.size:
        et = int(conne[i]); i += 1
        if et == 2:                      # Polyline: 次に節点数が来る
            nn = int(conne[i]); i += 1
        elif et in XDMF_NNODE:
            nn = XDMF_NNODE[et]
        elif et == 1:                    # Polyvertex
            nn = int(conne[i]); i += 1
        else:
            raise ValueError(f"unsupported XDMF element type {et} in {path}")
        faces.append([int(x) for x in conne[i:i + nn]]); i += nn
    return coords, faces, vals


def latest_wall_dump(run: Path, phys_name: str, phys_id: int) -> Path:
    pat = re.compile(rf"^res_{re.escape(phys_name)}_{phys_id}_(\d+)\.h5$")
    best, best_step = None, -1
    for p in run.glob(f"res_{phys_name}_{phys_id}_*.h5"):
        m = pat.match(p.name)
        if m and int(m.group(1)) > best_step:
            best, best_step = p, int(m.group(1))
    if best is None:
        raise FileNotFoundError(f"no wall dump res_{phys_name}_{phys_id}_*.h5 in {run}")
    return best


def last_wall_dumps(run: Path, phys_name: str, phys_id: int, n: int):
    """壁ダンプを step 順に並べ、**最後の n 枚**を返す (step 0 は初期状態なので除く)。"""
    pat = re.compile(rf"^res_{re.escape(phys_name)}_{phys_id}_(\d+)\.h5$")
    fs = []
    for p in run.glob(f"res_{phys_name}_{phys_id}_*.h5"):
        m = pat.match(p.name)
        if m and int(m.group(1)) > 0:
            fs.append((int(m.group(1)), p))
    fs.sort()
    return [p for _, p in fs[-n:]]


def latest_field(run: Path) -> Path:
    best, best_step = None, -1
    for p in run.glob("res_*.h5"):
        m = re.match(r"^res_(\d+)\.h5$", p.name)
        if m and int(m.group(1)) > best_step:
            best, best_step = p, int(m.group(1))
    if best is None:
        raise FileNotFoundError(f"no res_<step>.h5 in {run}")
    return best


# ------------------------------------------------------------------ loop
def write_wall_profile(path: Path, coords: np.ndarray, Tw: np.ndarray):
    with open(path, "w") as f:
        f.write("x y z Ts\n")
        for (x, y, z), T in zip(coords, Tw):
            f.write(f"{x:.10e} {y:.10e} {z:.10e} {T:.10e}\n")



def build_fem2d(spec: dict, wall_coords: np.ndarray):
    """固体 npz + 孔ごとの Robin から `Fem2DOperator` を作り、**壁ダンプ順 → 界面節点順**の
    並べ替え index を返す。

    界面は**流体の壁節点と 1 対 1** であることを要求する (座標一致 < 1e-7 m。壁ダンプの座標は
    float32 なので 1e-8 m 級の丸めが乗る。節点間隔 0.15 mm に対しては十分に厳しい)。
    `gen_solid_mesh.py --outer-from <compare_h が書く wall_nodes_ordered.csv>` で作った
    固体メッシュを渡すこと。一致しない節点があれば落とす (黙って内挿しない)。
    """
    d = np.load(spec["mesh_npz"])
    nodes, tris = d["nodes"], d["tris"]
    outer_edges = d["outer_edges"]
    iface = np.array(sorted(set(outer_edges.ravel().tolist())), int)
    holes = [k for k in d.files if k.startswith("hole")]
    robin = []
    hs = spec["holes"]
    if len(hs) not in (1, len(holes)):
        sys.exit(f"[cht_loop] solid.json の holes が {len(hs)} 個、メッシュの孔は {len(holes)} 個")
    for i, hk in enumerate(sorted(holes, key=lambda s_: int(s_[4:]))):
        hp = hs[i] if len(hs) > 1 else hs[0]
        for (n0, n1) in d[hk]:
            robin.append((int(n0), int(n1), float(hp["h"]), float(hp["T_c"])))
    k_solid = tuple(spec["k_table"]) if "k_table" in spec else float(spec["k_solid"])
    op = Fem2DOperator(nodes, tris, iface, [tuple(e) for e in outer_edges], robin, k_solid)

    # 壁ダンプ順 -> 界面節点順 の対応 (座標一致を要求する)
    W = np.asarray(wall_coords, float)[:, :2]
    S = op.coords[:, :2]
    if len(W) != len(S):
        sys.exit(f"[cht_loop] 壁節点 {len(W)} と固体界面節点 {len(S)} の数が違う "
                 "(gen_solid_mesh.py --outer-from で合わせること)")
    perm = np.empty(len(S), int)
    for i, p in enumerate(S):
        k = int(np.argmin(np.hypot(*(W - p).T)))
        dmin = float(np.hypot(*(W[k] - p)))
        if dmin > 1e-7:      # 壁ダンプの座標は float32 なので ~1e-8 m の丸めが乗る
            sys.exit(f"[cht_loop] 固体界面節点 {i} {p} に一致する壁節点が無い (最近傍 {dmin:.3e} m)")
        perm[i] = k
    if len(set(perm.tolist())) != len(perm):
        sys.exit("[cht_loop] 壁節点と固体界面節点の対応が 1 対 1 でない")
    return op, perm


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--template", required=True, help="mesh.h5 と config を置いたディレクトリ")
    ap.add_argument("--forge", required=True)
    ap.add_argument("--solid", required=True,
                    help="固体モデル JSON。`shell2d` は solid_shell.SolidModel、"
                         "`fem2d` は {mesh_npz, k_solid|k_table, holes:[{h,T_c}|...], T_init} を読む")
    ap.add_argument("--solid-mode", default="shell2d", choices=["shell2d", "fem2d"],
                    help="固体バックエンド。fem2d は一般 2D 断面 (冷却孔を Robin 辺で持つ)")
    ap.add_argument("--phys-id", type=int, required=True)
    ap.add_argument("--phys-name", default="wall")
    ap.add_argument("--max-iter", type=int, default=20)
    ap.add_argument("--flux", default="q_eff", choices=["q_compact", "q_recon", "q_2nd", "q_eff"],
                    help="界面に渡す熱流束の定義 (既定 q_eff = 壁半 CV に実際に入った保存形。"
                         "plan boundary-conjugate-heat-transfer §4.3 の正本。"
                         "**どれを使ったか履歴に残す**)")
    ap.add_argument("--tol-K", type=float, default=1.0e-3,
                    help="max|dTw| の収束許容 [K]。--flux-avg >= 2 のときは**壁温平均の傾きの許容 [K/反復]** として使う")
    ap.add_argument("--tol-rel", type=float, default=1.0e-3, help="max|r|/スケール の収束許容")
    ap.add_argument("--n-consec", type=int, default=2, help="収束と見なす連続回数")
    # ---- 界面ゲート G-if (plan §6、codex result M5) ----
    # **事前登録**して渡す。dT と res_rel だけでは、$D_f$ を上げて更新が小さくなっただけの
    # 状態を収束と認めてしまう (反例は `test_solid_shell.py` T7)。
    ap.add_argument("--tol-abs-W", type=float, default=None,
                    help="界面残差の**絶対**許容 [W] (単位奥行きなら W/m)。ケースごとに事前登録する")
    ap.add_argument("--tol-solid", type=float, default=None,
                    help="固体**内部**残差の許容 [W] (同上)。`fem2d` のみ評価される")
    ap.add_argument("--flux-avg", type=int, default=1,
                    help="界面熱流束を 1 反復の**最後の N 枚の壁ダンプで平均**する (既定 1 = 最後の 1 枚、従来どおり)。"
                         "N>=2 なら節点ごとの平均の標準誤差も出してドライバに渡し、受理・収束判定が"
                         "ノイズを知った形になる (plan boundary-conjugate-heat-transfer §5.1 #63)。"
                         "テンプレートの outStepInterval を、後半に N 枚以上入るように設定すること")
    ap.add_argument("--early-steps", type=int, default=None,
                    help="最初の --early-iters 反復だけ forge の nStepOuter をこの値にする (序盤は壁温が大きく動くので "
                         "流れを緩和させきる意味が薄い。plan boundary-conjugate-heat-transfer §5.1 #64)。"
                         "outStepInterval も N // (2*flux_avg) に合わせる。既定は無効")
    ap.add_argument("--early-iters", type=int, default=0)
    ap.add_argument("--anderson", type=int, default=5)
    ap.add_argument("--Tw-init", type=float, default=None, help="初期壁温 [K] (既定 = 固体の背面温度)")
    ap.add_argument("--Tg", type=float, default=None,
                    help="界面感度 D_f の初期値を**熱伝達係数** h=q/(Tg-Tw) から作るときの駆動温度 [K] "
                         "(既定は第一セルの伝導 k_eff A/d1)。")
    ap.add_argument("--Df-safety", type=float, default=2.0,
                    help="D_f 初期値の安全率 (過大なら遅いだけ、過小だと発散しうる)")
    ap.add_argument("--axisym", action="store_true")
    a = ap.parse_args()

    run = Path(a.run_dir).resolve()
    tpl = Path(a.template).resolve()
    run.mkdir(parents=True, exist_ok=True)
    if a.solid_mode == "shell2d":
        model = SolidModel.from_json(a.solid)
    else:
        model = json.loads(Path(a.solid).read_text())
    shutil.copy(a.solid, run / "solid.json")

    env = dict(os.environ)
    env.setdefault("LD_LIBRARY_PATH", "/usr/lib/x86_64-linux-gnu/hdf5/serial")

    # --- 前提の検査 (黙って効かない設定で回さない) ---
    bcond = (tpl / "bcondConfig.yaml").read_text()
    if "wallProfile" not in bcond:
        sys.exit(f"[cht_loop] {tpl}/bcondConfig.yaml: 対象壁に ints: {{wallProfile: 1}} が無い")
    solver_cfg = (tpl / "solverConfig.yaml").read_text()
    if "interfaceDiag" not in solver_cfg:
        sys.exit(f"[cht_loop] {tpl}/solverConfig.yaml: output に interfaceDiag: 1 が無い")

    hist_path = run / "cht_history.csv"
    hist = open(hist_path, "w", newline="")
    wr = csv.writer(hist)
    wr.writerow(["iter", "flux", "Tw_min", "Tw_max", "Tw_mean", "dTw_max", "res_rel",
                 "res_abs_W", "res_solid_W", "Q_total_W", "Df_mean", "used", "rejected",
                 "converged", "phi", "phi_noise", "at_floor", "sigma_rms_W", "sigma_floor_W", "drift_K_per_it"])
    hist.flush()

    op = drv = None
    Tw = None
    prev = None
    for it in range(a.max_iter):
        itd = run / f"it_{it:03d}"
        itd.mkdir(exist_ok=True)
        for f in ("mesh.h5", "solverConfig.yaml", "bcondConfig.yaml", "probe.yaml"):
            if (tpl / f).exists():
                shutil.copy(tpl / f, itd / f)
        early = (a.early_steps is not None) and (it < a.early_iters)
        if early:
            cfg_p = itd / "solverConfig.yaml"; txt = cfg_p.read_text()
            txt, n1 = re.subn(r"nStepOuter:\s*\d+", f"nStepOuter: {a.early_steps}", txt, count=1)
            txt, n2 = re.subn(r"outStepInterval:\s*\d+", f"outStepInterval: {max(a.early_steps // (2 * max(a.flux_avg, 1)), 1)}", txt, count=1)
            if n1 != 1 or n2 != 1:
                sys.exit(f"[cht_loop] --early-steps: {cfg_p} に nStepOuter / outStepInterval が見つからない")
            cfg_p.write_text(txt)
        if prev is not None:      # warm start (同一メッシュなので index コピー)
            subprocess.run([sys.executable, str(HERE / "interp_field.py"),
                            str(latest_field(prev)), str(itd / "mesh.h5")],
                           check=True, env=env, stdout=subprocess.DEVNULL)
        prof = itd / f"wall_profile_{a.phys_id}.csv"
        if Tw is not None:
            write_wall_profile(prof, op.coords, Tw)   # coords は固体界面節点の座標 (順不同で可)
        elif (tpl / f"wall_profile_{a.phys_id}.csv").exists():
            shutil.copy(tpl / f"wall_profile_{a.phys_id}.csv", prof)
        else:
            # 初回は壁の節点座標をまだ知らない (壁ダンプを読んで初めて分かる) ので、
            # **1 行だけの CSV** で一様な初期壁温を与える (3D 最近傍なので全面が同じ値になる)。
            T_init = a.Tw_init if a.Tw_init is not None else model.T_b
            write_wall_profile(prof, np.zeros((1, 3)), np.array([T_init]))

        print(f"[cht_loop] iter {it}: forge in {itd.relative_to(run.parent)}")
        with open(itd / "forge_run.log", "w") as log:
            rc = subprocess.run([a.forge], cwd=itd, stdout=log, stderr=subprocess.STDOUT, env=env)
        if rc.returncode != 0:
            sys.exit(f"[cht_loop] forge failed in {itd} (exit {rc.returncode}); see forge_run.log")

        coords, faces, vals = read_wall_dump(latest_wall_dump(itd, a.phys_name, a.phys_id))
        key = "iface_" + a.flux
        if key not in vals:
            sys.exit(f"[cht_loop] wall dump has no {key} (output.interfaceDiag: 1 が要る)")
        if "iface_ok" in vals and np.any(vals["iface_ok"] < 0.5):
            n_bad = int(np.sum(vals["iface_ok"] < 0.5))
            sys.exit(f"[cht_loop] {n_bad} wall nodes have no first interior point (iface_ok=0). "
                     "角・斜交で評価不能。幾何を見直すか --align-min を検討すること。")

        if op is None:
            if a.solid_mode == "shell2d":
                op = ShellOperator(coords, faces, model, axisym=a.axisym)
                T_back = model.T_b
            else:
                op, perm = build_fem2d(model, coords)
                T_back = float(model.get("T_init", np.mean([h["T_c"] for h in model["holes"]])))
            # **初期壁温は「実際に課した分布」から取る** (codex result M2, 2026-09-20)。
            # 旧実装はテンプレートの `wall_profile_*.csv` (実測分布) を forge に課しながら、
            # ドライバには一様値を持たせていた (実 run: 課 512-612 K / 仮定 566 K)。
            # 壁ダンプの `Ts` はまさに forge がその反復で使った壁温なので、これを T0 にする。
            if "Ts" in vals:
                T0 = np.asarray(vals["Ts"], float)
                if a.solid_mode == "fem2d":
                    T0 = T0[perm]
                print(f"[cht_loop]   T0 <- 壁ダンプの Ts ({T0.min():.2f}..{T0.max():.2f} K)")
            else:
                T0 = np.full(op.n, a.Tw_init if a.Tw_init is not None else T_back)
            # D_f の初期推定 = k_eff A / d1 (**上界ではない**。受理判定と退避で守る)
            keff, d1 = np.asarray(vals["iface_keff"], float), np.asarray(vals["iface_d1"], float)
            if a.solid_mode == "fem2d":
                keff, d1 = keff[perm], d1[perm]       # 壁ダンプ順 -> 界面節点順
            if a.Tg is not None:
                # **D_f は $\partial Q_f/\partial T_w$ = 熱伝達係数 × 面積**である。
                # 第一セルの伝導 $k_{\rm eff}A/d_1$ を使うと、境界層の厚み分だけ過大になる
                # ($k/d_1$ vs $k/\delta_T$)。実測 (case/53, $d_1$=2 µm): 13.7 vs 0.44 W/K で **31 倍**。
                # 過大な $D_f$ は安定だが**更新が止まる**ので、棄却のたびに倍加すると収束しない。
                qq0 = np.asarray(vals["iface_" + a.flux], float)
                if a.solid_mode == "fem2d":
                    qq0 = qq0[perm]
                dT = np.maximum(a.Tg - np.asarray(T0, float), 10.0)
                Df0 = a.Df_safety * np.maximum(np.abs(qq0) * op.area / dT, 1e-12)
            else:
                Df0 = a.Df_safety * np.maximum(keff * op.area / np.maximum(d1, 1e-12), 1e-12)
            drv = op.driver(T0, Df0=Df0, anderson=a.anderson)
            Tw = drv.T.copy()
            # 初回は壁温が config の一様値なので、そのまま 1 回目の Q_f を使う
        q = np.asarray(vals[key], float)                 # [W/m2] 固体向き正
        q_sem = None
        # **積分済み荷重** (ソルバが出す `iface_Qf_eff`)。あればこちらを正本にする。
        # **`--flux q_eff` のときだけ**正本にする (3 巡目 M3: 無条件に使うと `--flux q_compact` を
        # 指定した A/B でも q_eff の荷重が渡り、履歴には指定した名前が残って食い違う)。
        Qdirect = (np.asarray(vals["iface_Qf_eff"], float)
                   if (a.flux == "q_eff" and "iface_Qf_eff" in vals) else None)
        Qd_sem = None
        if a.flux_avg >= 2:
            dumps = last_wall_dumps(itd, a.phys_name, a.phys_id, a.flux_avg)
            if len(dumps) < a.flux_avg:
                sys.exit(f"[cht_loop] --flux-avg {a.flux_avg} だが壁ダンプが {len(dumps)} 枚しか無い "
                         f"({itd})。テンプレートの outStepInterval を細かくすること")
            allv = [read_wall_dump(d_)[2] for d_ in dumps]
            qs = np.array([np.asarray(v[key], float) for v in allv])
            q = qs.mean(axis=0)
            q_sem = qs.std(axis=0, ddof=1) / np.sqrt(len(dumps))
            if Qdirect is not None and all("iface_Qf_eff" in v for v in allv):
                Qs = np.array([np.asarray(v["iface_Qf_eff"], float) for v in allv])
                Qdirect = Qs.mean(axis=0)
                Qd_sem = Qs.std(axis=0, ddof=1) / np.sqrt(len(dumps))
        if a.solid_mode == "fem2d":
            q = q[perm]                                  # 壁ダンプ順 -> 固体界面節点順
            if q_sem is not None:
                q_sem = q_sem[perm]
        # **積分済み荷重があればそれを使う** (codex result 2 巡目 M4)。面積で割って集中辺長を
        # 掛け直すと、両者が違う角で荷重が歪む (Mark II 後縁で +29.9 %)。ソルバ内連成は
        # `iface_Qf_eff` を直接渡しており、外部ループも同じ契約に揃える。
        if Qdirect is not None:
            Qf = Qdirect[perm] if a.solid_mode == "fem2d" else Qdirect
            Qf_sigma = None if Qd_sem is None else (Qd_sem[perm] if a.solid_mode == "fem2d" else Qd_sem)
        else:
            Qf = q * op.area                             # 節点荷重 [W] (平面 2D は W/m)
            Qf_sigma = None if q_sem is None else q_sem * op.area
        Tw_new, info = drv.advance(Qf, tol_K=a.tol_K, tol_rel=a.tol_rel, n_consec=a.n_consec,
                                   tol_abs_W=a.tol_abs_W, tol_solid=a.tol_solid, Qf_sigma=Qf_sigma)
        wr.writerow([it, a.flux, f"{Tw.min():.6f}", f"{Tw.max():.6f}", f"{Tw.mean():.6f}",
                     f"{info['dT']:.6e}", f"{info['res_rel']:.6e}",
                     f"{info['res_abs']:.6e}", f"{info['res_solid']:.6e}",
                     f"{np.sum(Qf):.6e}",
                     f"{info['Df_mean']:.6e}", info["used"], int(info["rejected"]), int(info["converged"]),
                     f"{info['phi']:.6e}", f"{info.get('phi_noise', float('nan')):.6e}", int(info.get("at_floor", False)),
                     ("nan" if Qf_sigma is None else f"{float(np.sqrt(np.mean(Qf_sigma**2))):.6e}"),
                     f"{info.get('sigma_floor', 0.0):.6e}", f"{info.get('drift_K_per_it', float('nan')):.4e}"])
        hist.flush()
        print(f"[cht_loop]   Tw {Tw.min():.3f}..{Tw.max():.3f} K | dTw {info['dT']:.3e} K | "
              f"res {info['res_abs']:.3e} W ({info['res_rel']:.3e}) | "
              f"solid {info['res_solid']:.2e} | Q {np.sum(Qf):.4g} | {info['used']}"
              + (" REJECTED" if info["rejected"] else "") + (" [at noise floor]" if info.get("at_floor") else ""))
        prev, Tw = itd, Tw_new
        # 序盤 (短い流体評価) が収束判定の窓に入っている間は収束を宣言しない
        if info["converged"] and a.early_steps is not None and it < a.early_iters + max(a.n_consec, 4):
            info["converged"] = False
        if info["converged"]:
            if a.flux_avg >= 2:
                print(f"[cht_loop] CONVERGED at iter {it} (直近 {max(a.n_consec, 4)} 反復で壁温平均の傾き < {a.tol_K:g} K/反復、"
                      f"かつメリット関数が下げ止まり。壁温はその間の平均 -> Tw_final.csv)")
            else:
                print(f"[cht_loop] CONVERGED at iter {it} (dTw < {a.tol_K} K, res_rel < {a.tol_rel}, "
                      f"{a.n_consec} 回連続)")
            np.savetxt(run / "Tw_final.csv",
                       np.column_stack([op.coords, Tw]), delimiter=",",
                       header="x,y,z,Tw", comments="")
            hist.close()
            return 0
    hist.close()
    print(f"[cht_loop] NOT CONVERGED in {a.max_iter} iterations (see {hist_path})")
    return 1


if __name__ == "__main__":
    sys.exit(main())
