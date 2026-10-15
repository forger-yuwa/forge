r"""RANS 場からの排除厚 $\delta^*(x)$ 抽出 (A12 = axis-Mach 粘性補正)。

計画: plans/active/tooling-nozzle-axismach-viscous-deltastar.md §4。

定義 (圧縮性・**軸対称の質量収支** — 2026-09-01 に平面近似から置換):

$$\delta^*\left(1-\frac{\delta^*}{2r_w}\right) = \int_0^{y_e}
  \left(1 - \frac{\rho u}{(\rho u)_e}\right)\left(1-\frac{y}{r_w}\right) dy
  \;\Rightarrow\; \delta^* = r_w\left(1-\sqrt{1-2I/r_w}\right)$$

壁から $y$ 入った欠損は半径 $r_w-y$ の円環を流れるため $(1-y/r_w)$ の重みが付く。
旧・平面近似 ($\delta^*=\int(1-\rho u/(\rho u)_e)dy$) は欠損重心が壁近くに集中する
ため過大評価は小さい (M6 出口 $\delta/r_w\approx0.14$ で **+1.4 %**、出口 M 換算
+0.05 % — case/45 run_0007 実測)。帳簿の正確さのため軸対称形を正とする。

**教訓 (run_0030 で実測、繰り返し禁止)**:
- 探索窓はフリーストリームに届くまで取る (窓 0.35 r* は下流で BL 縁の手前で切れ、
  「δ* が下流ほど減る」という非物理な結論を出した — 撤回済み)。既定 1.5 r*。
- スロート近傍 (x ≲ 8 r_t) はコア流が半径方向に未一様で「ρu 極大 = BL 縁」の
  判定が破綻する。呼び出し側で相関へフォールバックすること。
"""
from __future__ import annotations

import h5py
import numpy as np
from scipy.interpolate import LinearNDInterpolator


def deltastar_from_run(mesh_h5, res_h5, wall_xy, scale: float,
                       x_stations, window: float = 1.5, n_prof: int = 400,
                       edge_frac: float = 0.995) -> dict:
    """node RANS run から各 x ステーションの δ* を測る。

    wall_xy: (n,2) 物理壁 [m 単位 or r_t 単位 — scale と整合していれば可]。
    x_stations: r_t 単位。戻り: dict(x, dstar, ue_rel, y_edge, ok) 全て r_t 単位。
    """
    with h5py.File(mesh_h5) as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    with h5py.File(res_h5) as f:
        ro = f["/VALUE/ro"][:]
        Ux = f["/VALUE/Ux"][:]
    if len(ro) != len(nc):
        raise ValueError("node run でない (VALUE 長 != 節点数)")
    x_n, r_n = nc[:, 0] / scale, nc[:, 1] / scale
    itp = LinearNDInterpolator(np.c_[x_n, r_n], ro * Ux)
    wx, wr = np.asarray(wall_xy)[:, 0], np.asarray(wall_xy)[:, 1]
    out = {"x": [], "dstar": [], "y_edge": [], "ok": []}
    for xs in np.atleast_1d(x_stations):
        rw = float(np.interp(xs, wx, wr))
        y = np.linspace(0.0, min(window, 0.95 * rw), n_prof)   # 壁 → 内側
        q = itp(np.full(n_prof, xs), rw - y)
        good = np.isfinite(q)
        if good.sum() < n_prof // 2:
            out["x"].append(float(xs)); out["dstar"].append(np.nan)
            out["y_edge"].append(np.nan); out["ok"].append(False)
            continue
        q = np.where(good, q, 0.0)
        qe = float(np.max(q))
        # 縁 = 壁から見て最初に edge_frac·(ρu)_max へ達する y
        idx = np.argmax(q >= edge_frac * qe)
        y_e = float(y[idx])
        m = y <= y_e
        # 軸対称の質量収支: 欠損に円環重み (1 - y/r_w) を掛け、
        # δ*(1 - δ*/(2 r_w)) = I を δ* について解く (docstring 参照)
        d = 1.0 - q[m] / max(qe, 1e-30)
        I = float(np.trapezoid(d * (1.0 - y[m] / rw), y[m]))
        integ = float(rw * (1.0 - np.sqrt(max(1.0 - 2.0 * I / rw, 0.0))))
        out["x"].append(float(xs)); out["dstar"].append(integ)
        out["y_edge"].append(y_e); out["ok"].append(bool(idx > 3))
    return {k: np.asarray(v) for k, v in out.items()}


# --- 固定 Euler 基準のコア整合質量欠損抽出 (plans/active/tooling-nozzle-deltastar-core-matched-euler.md) ----
REMAP_EDGES = ("clamp", "missing")
_EULER_REF_KEYS = ("run", "res", "sha256")
EULER_RANGE_PAD = 1e-9          # 既定 (remap_edge="clamp") の Euler の範囲の余裕 [r_t] (2026-10-14 までの全経路の値)
EULER_EDGE_TOL_FRAC = 1.0e-4    # remap_edge="missing" の端の外側の幾何の許容 = Euler の端の列の間隔 × これ (§5.1 #36g)


def euler_cover_range(xE, edge_tol: bool = True) -> dict:
    """Euler の参照の「補間できる範囲」[lo, hi] (r_t) を決める共通の関数 (plan tooling-nozzle-isothermal-wall-chain §5.1 #36g)。
    抽出器 (`deltastar_from_core_matched_euler` の列の被覆と `not_extracted` の理由)・δ_in の照合 (`delta_in_check`)・停止規則の帯の生成と
    δ の記録・受け渡しの関門の抽出と被覆 (いずれも抽出器を通る)・V-h370 の判定の A/B の比べる列が、すべてこれを呼ぶ。

    edge_tol=True (remap_edge="missing" の経路): lo = xE[0] − tol_lo、hi = xE[−1] + tol_hi、tol = EULER_EDGE_TOL_FRAC (1e−4) × その端の
    列の間隔 (tol_lo = 1e−4 (xE[1] − xE[0])、tol_hi = 1e−4 (xE[−1] − xE[−2]))。許容は**端の外側だけ**に足し、内側の列の被覆と評価
    (列 k と重み w の補間) は変えない。許容の内側で端を越える列は、抽出器の補間 (`euler_q_at`: 重み w を [0, 1] に丸める) により
    **端の列の値**で評価される。数学的には端の外側の微小区間 (幅 tol) への**定数の延長**であり、これを「外挿ではない」とは呼ばない
    (plan §5.1 #36i (v)、codex diagnose `notes/reviews/2026-10-14-h370-prelaunch-review-diagnose.md`)。
    **この許容は、同一の端点の丸め (本来同じ座標の端の列が、格子のファイルの桁数の違いで 1e−9 r_t 程度ずれる) を吸収するためのもので、
    範囲の不足を隠すためではない**。2026-10-14 の B の dry-run で、変換後の座標の NS の最後の列 x = 95.29740274407894 が Euler の参照の端
    95.29740274239427 を 1.68e−9 r_t 越えて被覆が落ちた (Euler の格子は msh の座標を既定の 10 桁で書くので、出口の x × r_t =
    7.306010214129155 m が 7.306010214 m に丸まり、÷ r_t = 95.29740274239427 と一致する。NS の格子は 17 桁なので生成の値のまま)。この種の
    差 (≲ 1e−8 r_t) は許容 (端の間隔の 1e−4 倍) より十分小さく、範囲の不足 (旧い run_0174 の終端 95.2276 と NS の 95.2974 の 0.07 r_t) は
    許容より桁違いに大きいので、後者は今どおり欠測になる。
    edge_tol=False (既定の remap_edge="clamp" の経路): 今までどおり ±EULER_RANGE_PAD (1e−9 r_t) — 既定の呼び出し元の結果は変えない。
    戻り: {"lo", "hi", "x0", "x1", "tol_lo", "tol_hi", "edge_tol", "frac" (edge_tol のとき)}。"""
    xE = np.asarray(xE, dtype=float)
    if xE.ndim != 1 or len(xE) < 2:
        raise ValueError("euler_cover_range: Euler の列の x は 2 列以上の 1 次元の配列")
    if edge_tol:
        tl = EULER_EDGE_TOL_FRAC * float(xE[1] - xE[0])
        th = EULER_EDGE_TOL_FRAC * float(xE[-1] - xE[-2])
        if not (tl > 0.0 and th > 0.0):
            raise ValueError("euler_cover_range: Euler の端の列の間隔が正でない (x が単調増加でない)")
    else:
        tl = th = EULER_RANGE_PAD
    out = {"lo": float(xE[0]) - tl, "hi": float(xE[-1]) + th, "x0": float(xE[0]), "x1": float(xE[-1]), "tol_lo": tl, "tol_hi": th,
           "edge_tol": bool(edge_tol)}
    if edge_tol:
        out["frac"] = EULER_EDGE_TOL_FRAC
    return out


def euler_covers(x, rng: dict):
    """x (スカラーか配列) が `euler_cover_range` の範囲 [lo, hi] の内か (閉区間)。スカラーなら bool、配列なら bool の配列。"""
    c = (np.asarray(x, dtype=float) >= rng["lo"]) & (np.asarray(x, dtype=float) <= rng["hi"])
    return bool(c) if np.ndim(c) == 0 else c


def _sha256_file(path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def resolve_euler_ref(spec, verify: bool = True) -> dict:
    """Euler の参照の指定を解決する (plan tooling-nozzle-isothermal-wall-chain §5.1 #36e (6)、codex diagnose 3: 抽出器は参照の
    ディレクトリの最新の出力を読むので、計算中の Euler を参照すると停止の帯と後段の抽出の参照が変わりうる)。

    - str / Path (**ディレクトリだけ、今の挙動**): {"run": Path(spec) (そのまま), "res": None (最新の res_<N>.h5), "pinned": False}。
    - mapping {"run": ディレクトリ, "res": "res_<N>.h5" (ファイル名だけ), "sha256": 64 桁の 16 進} (**固定のスナップショット**): run は
      そのまま (相対なら cwd、無ければリポジトリ直下から)。verify なら <run>/<res> の sha256 を計算して指定と一致することを確かめ、
      無い・違えば ValueError (sha が合わなければ使わない)。
    戻り: {"run", "res", "sha256", "pinned", "record" (JSON に書ける形: ディレクトリだけなら str、固定なら mapping)}。"""
    import re
    from pathlib import Path
    if isinstance(spec, (str, Path)):
        return {"run": Path(spec), "res": None, "sha256": None, "pinned": False, "record": str(spec)}
    if not isinstance(spec, dict):
        raise ValueError(f"Euler の参照の指定は run dir (str) か {{run, res, sha256}} の mapping (受け取った {spec!r})")
    if set(spec) != set(_EULER_REF_KEYS):
        raise ValueError(f"Euler の参照の固定の指定はキー {list(_EULER_REF_KEYS)} ちょうど (受け取った {sorted(spec)})")
    run, res, sha = spec["run"], spec["res"], spec["sha256"]
    if not isinstance(run, (str, Path)) or not str(run):
        raise ValueError(f"Euler の参照の run が不正 ({run!r})")
    if not isinstance(res, str) or not re.fullmatch(r"res_\d+\.h5", res):
        raise ValueError(f"Euler の参照の res は res_<N>.h5 のファイル名 (ディレクトリを含まない) ({res!r})")
    if not isinstance(sha, str) or not re.fullmatch(r"[0-9a-f]{64}", sha):
        raise ValueError(f"Euler の参照の sha256 は 64 桁の小文字の 16 進 ({sha!r})")
    cand = [Path(run)] + ([] if Path(run).is_absolute() else [Path(__file__).resolve().parents[3] / run])
    rd = next((c for c in cand if c.is_dir()), Path(run))
    if verify:
        f = rd / res
        if not f.is_file():
            raise ValueError(f"Euler の参照の固定のスナップショット {f} が無い")
        got = _sha256_file(f)
        if got != sha:
            raise ValueError(f"Euler の参照の固定のスナップショット {f} の sha256 {got[:16]}… が指定 {sha[:16]}… と違う — 使わない")
    return {"run": rd, "res": res, "sha256": sha, "pinned": True, "record": {"run": str(run), "res": res, "sha256": sha}}


def _load_structured(run_dir, remap_constant_x: bool = False, remap_edge: str = "clamp", res: str | None = None,
                     res_name: str | None = None):
    """node 構造格子 run (index = i*nj + j) を (info, x[ni,nj], r[ni,nj], q=ρUx[ni,nj]) で返す (r_t 単位)。
    res: 読む出力のファイル名 (None = 最新の res_<N>.h5; 今の挙動)。Euler の参照の固定 (`resolve_euler_ref`) が使う。

    remap_constant_x (2026-10-08, plan tooling-nozzle-isothermal-wall-chain §5.1 #21、codex diagnose): 近壁を壁法線に沿わせた格子
    (`mesh2d.wall_normal_layer`) では列 i の中で x が変わる (x ≈ 40 で列内の幅 0.15 r_t)。True なら各行 j に沿って x 方向に線形補間し、
    列 i を壁の station の x (= 軸の x、x[i, 0]) の一定 x の断面にそろえる。行の x の範囲の外 (端の補間範囲の不足) の扱いは remap_edge:
    "clamp" (既定、今の挙動) は端値 (両端の数列だけ)、"missing" (§5.1 #36e (4)、codex diagnose 3「端の補間範囲の不足を定数の外挿で
    隠して全域合格にしない」) はその列を欠測にする (列の r・q を NaN にし、`remap_missing_cols` に列の番号と x を記録。抽出器はその列を
    抽出せず not_extracted に記録する)。False (既定) は従来どおり列をそのまま使う (列内で x が一定の格子では両者は一致する)。
    res_name (2026-10-12, plan tooling-nozzle-bell-two-mode §5.1 #7 (7b-1) M3): `res` の別名 (③ベルの δ_r の時系列が使う)。指定したファイルが無ければ FileNotFoundError。"""
    import json
    from pathlib import Path
    if remap_edge not in REMAP_EDGES:
        raise ValueError(f"remap_edge={remap_edge!r} は {REMAP_EDGES} のどれか")
    run_dir = Path(run_dir)
    info = json.loads((run_dir / "prepare_info.json").read_text())
    ni, nj = int(info["mesh"]["ni"]), int(info["mesh"]["nj"])
    S = float(info["scale_m"])
    with h5py.File(run_dir / "nozzle.h5") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    if res is None and res_name is not None:
        res = res_name
    if res is None:
        res = sorted(run_dir.glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))[-1]
    else:
        res = run_dir / str(res)
        if not res.is_file():
            raise FileNotFoundError(f"{run_dir}: スナップショット {res.name} が無い")
    with h5py.File(res) as f:
        q = f["/VALUE/ro"][:] * f["/VALUE/Ux"][:]
    if len(q) != ni * nj or len(nc) != ni * nj:
        raise ValueError(f"{run_dir}: node 構造格子 run でない (VALUE {len(q)}, COORD {len(nc)}, ni*nj {ni*nj})")
    x = nc[:, 0].reshape(ni, nj) / S; r = nc[:, 1].reshape(ni, nj) / S; q = q.reshape(ni, nj)
    out = dict(info=info, S=S, res=res.name, x=x, r=r, q=q)
    if remap_constant_x:
        xt = x[:, 0].copy()                                  # 列の目標の x (軸の節点 = 壁の station の x)
        rr = np.empty_like(r); qq = np.empty_like(q)
        clamped = 0
        outside = np.zeros(ni, dtype=bool)
        for j in range(nj):
            xj = x[:, j]
            if not np.all(np.diff(xj) > 0):
                raise ValueError(f"{run_dir}: 行 {j} の x が単調でない (remap_constant_x)")
            oj = (xt < xj[0]) | (xt > xj[-1])
            clamped += int(np.sum(oj))
            outside |= oj
            rr[:, j] = np.interp(xt, xj, r[:, j]); qq[:, j] = np.interp(xt, xj, q[:, j])
        out.update(x=np.repeat(xt[:, None], nj, axis=1), r=rr, q=qq, remap_constant_x=True, remap_clamped_points=clamped,
                   remap_max_column_dx=float(np.max(np.ptp(x, axis=1))), remap_edge=remap_edge)
        if remap_edge == "missing":
            miss = np.where(outside)[0]
            rr[miss, :] = np.nan; qq[miss, :] = np.nan
            out["remap_missing_cols"] = {"i": [int(i) for i in miss], "x": [float(v) for v in xt[miss]]}
    return out


def _exit_x(info: dict):
    """流量の帳簿の区間の下流端 x_E (r_t 単位) と出所。①風洞の run は prepare_info の `x_E`。③ベルの run (prepare_info に
    `bell_mode` があり `x_E` が無い — 2026-10-12, plan tooling-nozzle-bell-two-mode §5.1 #7 (7b-1)) は壁の設計の L
    (`wall_design.L` = 出口面、領域の下流端) を使う。どちらも無ければ (None, 理由)。"""
    if "x_E" in info:
        return float(info["x_E"]), "prepare_info.x_E"
    wd = info.get("wall_design") or {}
    if info.get("bell_mode") and wd.get("L") is not None:
        return float(wd["L"]), "prepare_info.wall_design.L (bell の出口面; x_E 無し)"
    return None, "prepare_info に x_E も (bell の) wall_design.L も無い"


def massflow_ratio(ns_run, euler_run, ns_res: str | None = None, remap_constant_x: bool = False) -> dict:
    """NS/Euler 質量流量比 → 有効スロート半径・質量流量由来の等価スロート補正量 (r_t 単位)。
    列ごとに 2π∫ρUx r dr (台形則)、x∈(−2.5, x_E) の中央値を採用。euler_run は run dir か固定の指定 (`resolve_euler_ref`)。
    x_E は `_exit_x` (①は prepare_info.x_E、③ベルは wall_design.L)。どちらも無い run は止めずに帳簿を省き、
    `{"skipped": True, "reason": ...}` を返す (黙って省かない)。ns_res: NS のスナップショット名 (None = 最新; `_load_structured`)。
    remap_constant_x: 列を一定 x の断面に直してから積分する (`_load_structured` の同名; False = 従来どおり列のまま)。列の中で x が
    変わる格子では 2π∫ρU_x r dr は断面の流量でない (−ρU_r dx の項が落ちる): ③ L100 の NS (列内 0.28 r_t) で列のままは
    260.17 kg/s・断面間の散らばり 1.8 %、一定 x は 268.34 kg/s・0.03 % (run の出口の流量 268.36; 2026-10-12)。"""
    out = {}
    eref = resolve_euler_ref(euler_run)
    for tag, rd, res in (("ns", ns_run, ns_res), ("euler", eref["run"], eref["res"])):
        d = _load_structured(rd, remap_constant_x=remap_constant_x, res=res)
        S = d["S"]
        # 実寸の流量 2π∫ρU_x r dr [kg/s]: 被積分 q·(r S) と積分変数 (r S) で S は 2 乗 (2026-10-08 修正: 以前は被積分にも S² を掛けて
        # 3 乗になり、絶対値が S 倍に小さく、r_t の違う run どうしの比も S の比だけずれていた。診断の帳簿だけで使う値)
        m = 2 * np.pi * np.trapezoid(d["q"] * d["r"] * S, d["r"] * S, axis=1)
        xs = d["x"][:, 0]
        x_E, x_E_src = _exit_x(d["info"])
        if x_E is None:
            skip = {"skipped": True, "reason": f"{tag} run {rd}: {x_E_src} — 流量の帳簿 (x∈(−2.5, x_E) の中央値) を省いた"}
            skip.update({k: v for k, v in out.items() if k.startswith("res_")})
            skip[f"res_{tag}"] = d["res"]
            return skip
        if x_E_src != "prepare_info.x_E":                   # ①の run (x_E あり) では鍵を増やさない (戻り値を従来とビット同一に保つ)
            out[f"x_E_{tag}"] = x_E
            out[f"x_E_{tag}_source"] = x_E_src
        if remap_constant_x:
            out[f"remap_constant_x_{tag}"] = True
        sel = (xs > -2.5) & (xs < x_E)
        out[f"mdot_{tag}"] = float(np.median(m[sel]))
        out[f"mdot_{tag}_spread"] = float(m[sel].std() / np.median(m[sel]))
        out[f"res_{tag}"] = d["res"]
        if tag == "ns":
            r_tW = float(d["info"].get("throat_physical", {}).get("r", np.nan)) if d["info"].get("throat_physical") else float(np.min(d["r"][:, -1]))
    ratio = out["mdot_ns"] / out["mdot_euler"]
    out["mdot_ratio"] = float(ratio)
    out["r_throat_eff"] = float(np.sqrt(ratio))            # 設計 r_t 単位
    out["r_throat_phys"] = r_tW
    out["cd_rel"] = float(ratio / r_tW ** 2)
    out["delta_r_throat_eff"] = float(r_tW - np.sqrt(ratio))
    out["delta_r_throat_given"] = float(r_tW - 1.0)
    return out



def _negative_delta_r(D: float, rwN: float, q_wall_ref: float) -> float:
    r"""符号付き δ_r (2026-09-12, plan tooling-nozzle-isothermal-wall-chain §4.1 / codex M1):
    質量欠損 $D<0$ (冷却壁で壁近傍の ρu が参照より大きい = 過剰) のときは、参照 ρu を壁値で壁の外へ延長し
    $2\pi\int_{r_w}^{r_{eff}} q_{ref,w}\, r\,dr = -D$ を解く: $r_{eff}^2 = r_w^2 - D/(\pi q_{ref,w})$、
    $\delta_r = r_w - r_{eff} < 0$ (壁が実効的に外へ動く)。旧実装の線形化 $D/(2\pi r_w \bar q)$ を置換 (診断値だった)。"""
    qw = max(float(q_wall_ref), 1e-30)
    return float(rwN - np.sqrt(rwN * rwN - D / (np.pi * qw)))


def core_matched_deficit(r, q_ns, q_e, rw_e: float, core_frac: float = 0.30,
                         n_axis_skip: int = 1, outer_frac: float = 0.25, refine: int = 4):
    r"""1 断面の**コア整合質量欠損 → 半径方向等価排除厚** (純関数、単体試験対象)。

    r: 軸→壁の半径 (NS 格子, 最後が NS 壁)、q_ns: NS の ρUx、q_e: 同じ r に写した Euler の ρUx
    (Euler 壁 rw_e より外は壁値) — 配列 (r と同長) か callable(r)。
    refine>1 で NS 区間を細分 (q_ns は区分線形、q_e は callable なら細分点で評価、配列なら線形補間)。
    壁集中メッシュはコアに節点が ~5 個しか無いので、積分と α を安定させるために使う。
    戻り: dict または None (コア点数不足)。"""
    r0 = np.asarray(r, dtype=float); q0 = np.asarray(q_ns, dtype=float)
    if refine > 1:
        t = np.linspace(0.0, 1.0, refine, endpoint=False)
        r = np.concatenate([(r0[:-1, None] + (r0[1:] - r0[:-1])[:, None] * t[None, :]).ravel(), [r0[-1]]])
        qN = np.interp(r, r0, q0)
        qE = np.asarray(q_e(r), dtype=float) if callable(q_e) else np.interp(r, r0, np.asarray(q_e, dtype=float))
    else:
        r, qN = r0, q0
        qE = np.asarray(q_e(r), dtype=float) if callable(q_e) else np.asarray(q_e, dtype=float)
    rwN = float(r[-1])
    mc = r <= core_frac * rw_e
    if mc.sum() < n_axis_skip + 3:
        return None
    num = np.trapezoid(qN[mc] * r[mc], r[mc]); den = np.trapezoid(qE[mc] * r[mc], r[mc])
    alpha = float(num / den) if den > 0 else np.nan
    qref = alpha * qE
    rel = qN[mc] / np.maximum(qref[mc], 1e-30) - 1.0
    core_rms = float(np.sqrt(np.mean(rel ** 2)))
    core_rms_noaxis = float(np.sqrt(np.mean(rel[n_axis_skip:] ** 2)))
    core_maxdev = float(np.max(np.abs(rel[n_axis_skip:])))
    integrand = (qref - qN) * r
    D = float(2 * np.pi * np.trapezoid(integrand, r))
    # 壁からの累積 F(y) = 2π∫_{rw-y}^{rw} q_ref r' dr' (y = 壁からの距離, 単調増加)
    seg = 0.5 * (qref[1:] * r[1:] + qref[:-1] * r[:-1]) * np.diff(r)
    F_from_wall = 2 * np.pi * np.concatenate([[0.0], np.cumsum(seg[::-1])])
    r_desc = r[::-1]
    reason = []
    if D < 0:
        delta_r = _negative_delta_r(D, rwN, qref[-1])                      # 符号付き (負値 = 過剰質量)
        reason.append("negative_deficit")
    elif D > F_from_wall[-1]:
        delta_r = float("nan"); reason.append("no_root")
    else:
        r_eff = float(np.interp(D, F_from_wall, r_desc))
        delta_r = rwN - r_eff
    segD = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r)
    G_from_wall = 2 * np.pi * np.concatenate([[0.0], np.cumsum(segD[::-1])])
    y_desc = rwN - r_desc
    G_outer = float(np.interp(outer_frac * rwN, y_desc, G_from_wall))
    outer_share = G_outer / D if D > 0 else float("nan")
    return dict(r_wall_euler=float(rw_e), r_wall_ns=rwN, delta_in=rwN - float(rw_e), alpha=alpha,
                core_rms=core_rms, core_rms_noaxis=core_rms_noaxis, core_maxdev=core_maxdev,
                mass_deficit=D, delta_r=float(delta_r), outer_share=outer_share, reason=reason)



def band_local_deficit(r, q_ns, q_e, rw_e: float, delta_in: float | None = None,
                       k_band: float = 1.5, band_min_frac: float = 0.02, band_width: float = 0.5,
                       n_axis_skip: int = 1, core_frac: float = 0.30, refine: int = 4, _retry: int = 0,
                       slope_max: float = 0.01, band_growth: float = 1.25, max_retry: int = 14,
                       fit_from: float = 1.0, y_b_fixed: float | None = None):
    r"""1 断面の**帯局所参照・質量欠損 → 半径方向等価排除厚** (生産抽出、純関数)。

    plan §4.3 の修正 (2026-09-04 V1 pass1 の知見): 内側 30 % コアの単一倍率 α では、上流の壁 δ 誤差が
    特性線で下流の軸へ運ぶ波 (コア形状差) を欠損に取り込む。そこで参照 $q_{ref}$ は
    **境界層のすぐ外側の帯** $y\in[y_b, (1+w)y_b]$ で比 $q_{NS}/q_E$ を $y$ の 1 次で LSQ フィットし、
    境界層域 $y<y_b$ へ延長して作る。$y_b$ は $\max(k\,\delta_{in},\ f_{min} r_w)$ (既定 k=1.5 = 境界層の中) から始め、
    帯内の比の変化が slope_max (1 %) を下回るまで band_growth (1.25) 倍ずつ広げる (適応帯: 4 つの異なる壁の NS 場で
    δ_r のばらつき 3 %、固定 k=4/w=1 では 6.5 % — case/45 run_0004/0019/0020/0021)。
    欠損 $D = 2\pi\int_{y<y_b}(q_{ref}-q_{NS})r\,dr$、$\delta_r = r_w - r_{eff}$。
    $r_{eff}$ は壁〜$y_b$ の区間の累積流量 $F$ だけで引く (2026-10-13, plan tooling-nozzle-deltastar-band-root-zone §4.1)。
    探索域で $q_{ref}$ が非正・非有限なら `nonpos_ref`、区間の中で $F$ が減れば `nonmono_zone` (どちらも hard、δ_r = NaN)。
    内側 30 % の α・コア RMS は診断として併記する。q_e は配列または callable(r)。"""
    r0 = np.asarray(r, dtype=float); q0 = np.asarray(q_ns, dtype=float)
    if refine > 1:
        t = np.linspace(0.0, 1.0, refine, endpoint=False)
        rr = np.concatenate([(r0[:-1, None] + (r0[1:] - r0[:-1])[:, None] * t[None, :]).ravel(), [r0[-1]]])
        qN = np.interp(rr, r0, q0)
        qE = np.asarray(q_e(rr), dtype=float) if callable(q_e) else np.interp(rr, r0, np.asarray(q_e, dtype=float))
    else:
        rr, qN = r0, q0
        qE = np.asarray(q_e(rr), dtype=float) if callable(q_e) else np.asarray(q_e, dtype=float)
    rwN = float(rr[-1]); y = rwN - rr
    d_in = float(rwN - rw_e) if delta_in is None else float(delta_in)
    if y_b_fixed is not None:                      # 2 パス目: 隣接列で平滑化した帯位置を固定して使う (帯のトグルを防ぐ)
        y_b0 = float(y_b_fixed); max_retry = 0
    else:
        y_b0 = max(k_band * max(d_in, 0.0), band_min_frac * rwN) * (band_growth ** _retry)   # 適応探索する「縁」候補
    y_b = fit_from * y_b0                                                                 # 参照フィット帯の下端 = 積分域の上端
    y_t = min((1.0 + band_width) * y_b, 0.85 * rwN)
    if y_b >= 0.8 * rwN:
        return None
    mb = (y >= y_b) & (y <= y_t)
    if mb.sum() < 4:
        return None
    ratio = qN[mb] / np.maximum(qE[mb], 1e-30)
    A = np.c_[np.ones(mb.sum()), y[mb] - y_b]
    coef, *_ = np.linalg.lstsq(A, ratio, rcond=None)
    a, b = float(coef[0]), float(coef[1])
    band_rms = float(np.sqrt(np.mean((ratio - (a + b * (y[mb] - y_b))) ** 2)))
    # 縁候補帯 [y_b0, (1+w) y_b0] での比の変化 (適応判定用; fit_from=1 なら上と同じ帯)
    mb0 = (y >= y_b0) & (y <= min((1.0 + band_width) * y_b0, 0.85 * rwN))
    if mb0.sum() >= 4:
        c0, *_ = np.linalg.lstsq(np.c_[np.ones(mb0.sum()), y[mb0] - y_b0], qN[mb0] / np.maximum(qE[mb0], 1e-30), rcond=None)
        slope0 = float(c0[1]) * (min((1.0 + band_width) * y_b0, 0.85 * rwN) - y_b0)
        a0 = float(c0[0])
    else:
        slope0, a0 = 0.0, a
    # 内側コアの単一倍率 (診断 + 「帯が境界層の深部にある」検出: 縁候補の比がコア比より 5 % 以上低い)
    mc = rr <= core_frac * rw_e
    if mc.sum() >= n_axis_skip + 3:
        num = np.trapezoid(qN[mc] * rr[mc], rr[mc]); den = np.trapezoid(qE[mc] * rr[mc], rr[mc])
        alpha = float(num / den) if den > 0 else float("nan")
        rel = qN[mc] / np.maximum(alpha * qE[mc], 1e-30) - 1.0
        core_rms = float(np.sqrt(np.mean(rel[n_axis_skip:] ** 2)))
    else:
        alpha, core_rms = float("nan"), float("nan")
    if np.isfinite(alpha) and a0 < 0.95 * alpha and _retry < max_retry:
        return band_local_deficit(r, q_ns, q_e, rw_e, delta_in, k_band, band_min_frac, band_width,
                                  n_axis_skip, core_frac, refine, _retry + 1, slope_max=slope_max,
                                  band_growth=band_growth, max_retry=max_retry, fit_from=fit_from)
    # 帯が境界層の中にあると比 q_NS/q_E が帯の幅で数 % 以上変わる (コアの波は < 1 %)。
    # 線形フィットは BL の形も吸収してしまい残留欠損では検出できないので、帯内の比の変化量で判定し y_b を倍にする。
    if abs(slope0) > slope_max and _retry < max_retry:
        return band_local_deficit(r, q_ns, q_e, rw_e, delta_in, k_band, band_min_frac, band_width,
                                  n_axis_skip, core_frac, refine, _retry + 1, slope_max=slope_max,
                                  band_growth=band_growth, max_retry=max_retry, fit_from=fit_from)
    ratio_fn = a + b * (y - y_b)
    qref = ratio_fn * qE
    mz = y <= y_b                                     # 境界層域 (壁側)
    integrand = (qref - qN) * rr
    D = float(2 * np.pi * np.trapezoid(integrand[mz], rr[mz]))
    # 帯の中の残留欠損 (参照の整合診断): 本来 ~0。大きければ帯が境界層の中にある → y_b を 2 倍して再試行
    # (δ99/δ* は 2〜10 と場所で変わるので k_band 固定では足りない。閾値でなく自己整合で縁を外す)
    D_band = float(2 * np.pi * np.trapezoid(integrand[mb], rr[mb]))
    if D > 0 and abs(D_band) > 0.1 * D and _retry < max_retry:
        return band_local_deficit(r, q_ns, q_e, rw_e, delta_in, k_band, band_min_frac, band_width,
                                  n_axis_skip, core_frac, refine, _retry + 1, slope_max=slope_max,
                                  band_growth=band_growth, max_retry=max_retry, fit_from=fit_from)
    seg = 0.5 * (qref[1:] * rr[1:] + qref[:-1] * rr[:-1]) * np.diff(rr)
    F_from_wall = 2 * np.pi * np.concatenate([[0.0], np.cumsum(seg[::-1])])
    r_desc = rr[::-1]
    reason = []
    # 探索域 (壁〜帯の下端 y_b) の参照の検査 (plan tooling-nozzle-deltastar-band-root-zone §4.1、codex plan レビュー M2): 帯のフィットを
    # 壁側へ外挿した q_ref が探索域で非正・非有限なら、D の符号によらず等価半径は定義できない (D < 0 の枝は q_ref の壁値を 1e-30 に
    # 置き換えて、codex の再現例で δ_r = −4.98e13 を soft で返していた) → hard (nonpos_ref) で NaN。正の参照に対する負の δ_r (冷却壁) は
    # そのまま
    qz = qref[mz]
    # 探索域の上端 y_b そのもの (格子点の間にあることが多い) の参照値も見る: 帯のフィットは y_b で a、q_E は線形補間 (codex result レビュー 2026-10-13 M1;
    # 格子点だけでは y_b の手前で q_ref が負に横切る断面を見逃し、D < 0 の枝で soft の negative_deficit を返していた)
    qz = np.append(qz, a * float(np.interp(rwN - y_b, rr, qE)))
    if not (np.all(np.isfinite(qz)) and np.all(qz > 0.0)):
        delta_r = float("nan"); reason.append("nonpos_ref")
    elif D < 0:
        delta_r = _negative_delta_r(D, rwN, qref[-1])                      # 符号付き (負値 = 過剰質量)
        reason.append("negative_deficit")
    else:
        F_zone = float(np.interp(rwN - y_b, r_desc[::-1], F_from_wall[::-1]))  # 帯下端までの容量
        if D > F_zone:
            if _retry < max_retry:
                return band_local_deficit(r, q_ns, q_e, rw_e, delta_in, k_band, band_min_frac, band_width,
                                          n_axis_skip, core_frac, refine, _retry + 1, slope_max=slope_max,
                                          band_growth=band_growth, max_retry=max_retry, fit_from=fit_from)
            delta_r = float("nan"); reason.append("no_root")
        else:
            # 等価半径は壁〜y_b の区間の F だけで引く (§4.1)。全半径の F で引くと、帯の傾きが大きくコアで q_ref < 0 の断面では F が
            # 途中から減って末尾が D を下回り、np.interp が軸の半径を返して δ_r = r_w になっていた (case/70 の x ≈ 1〜3)。
            # 区間の端は y_b で線形補間した点 (F_zone) で、それは y_b を挟む 2 格子点を結ぶ線分の上にある。根はその線分 (全半径の探索と
            # 同じ 2 点) で引くので、全半径で F が単調な断面では結果はビット同一
            y_desc = rwN - r_desc
            kz = int(np.searchsorted(y_desc, y_b, side="right"))           # 壁から y ≤ y_b の格子点の数
            hi = kz + 1 if (kz < len(y_desc) and y_desc[kz - 1] < y_b) else kz   # y_b が格子点の間なら、それを挟む外側の点まで
            if np.any(np.diff(F_from_wall[:hi]) < 0.0):
                delta_r = float("nan"); reason.append("nonmono_zone")      # 区間の中で F が減る (hard)
            else:
                r_eff = float(np.interp(D, F_from_wall[:hi], r_desc[:hi]))
                delta_r = rwN - r_eff
    return dict(r_wall_euler=float(rw_e), r_wall_ns=rwN, delta_in=rwN - float(rw_e), alpha=alpha,
                core_rms=core_rms, core_rms_noaxis=core_rms, core_maxdev=float("nan"),
                mass_deficit=D, delta_r=float(delta_r), outer_share=float("nan"),
                band_y_b=y_b0, band_fit_from=fit_from, band_ratio_a=a, band_slope=slope0, band_rms=band_rms,
                band_deficit_share=(D_band / D if D > 0 else float("nan")), band_retry=_retry,
                reason=reason)


def pspline_uniform(x, v, knot: float, lam: float = 1.0, weights=None, k: int = 3, penalty_order: int = 2):
    """一様ノット (間隔 knot) の B-spline + 係数差分ペナルティ。weights=0 / 非有限の点は無視。戻り: x での値。"""
    from scipy.interpolate import BSpline
    x = np.asarray(x, dtype=float); v = np.asarray(v, dtype=float)
    w = np.ones_like(v) if weights is None else np.asarray(weights, dtype=float)
    ok = np.isfinite(v) & (w > 0)
    lo, hi = float(x.min()), float(x.max())
    n_int = max(int(np.ceil((hi - lo) / knot)) - 1, 1)
    inner = np.linspace(lo, hi, n_int + 2)[1:-1]
    t = np.concatenate([[lo] * (k + 1), inner, [hi] * (k + 1)])
    n_cp = len(t) - k - 1
    B = BSpline.design_matrix(x[ok], t, k).toarray()
    D = np.diff(np.eye(n_cp), n=penalty_order, axis=0)
    W = w[ok][:, None]
    c = np.linalg.solve(B.T @ (W * B) + lam * (D.T @ D), B.T @ (w[ok] * v[ok]))
    return BSpline(t, c, k, extrapolate=True)(x)


def edge_band_positions(x, y_ladders, d_ladders, d_in, r_w, eps: float = 0.003, c: float = 1.25,
                        window: float = 3.0, knot: float = 6.0, edge_ratio: float = 1.25,
                        ref_min_frac: float = 0.0133, lam: float = 1.0, crossing: str = "quadratic",
                        hold_lo: float | None = 3.0, hold_hi_from_end: float | None = 3.0, persist: int = 1) -> dict:
    r"""帯の位置を「境界層の縁の連続推定 → x 方向平滑化 → その c 倍」で決める (band_select="edge"、純関数)。

    plan verification-m6-axis-wave-mesh-su2 §5.1 #7a (diagnostician 2026-10-04 推薦 E)。各断面 i で、帯を固定した
    δ_r の梯子 (y_ladders[i], d_ladders[i]) から $R(y)=\delta_r(\text{edge\_ratio}\,y)/\delta_r(y)-1$ を作り
    (log y で線形補間)、**そこから persist 段続けて** $|R|<\varepsilon$ を満たす最初の段を交点の段とし、その段と 1 つ内側の段の間を
    **格子点間の線形補間で**解いて縁 $y^*$ とする (無ければ NaN)。persist = 1 は「$|R|<\varepsilon$ を初めて満たす段」(2026-10-13 までの式)。
    既定は persist = 1。persist > 1 (opt-in) は一番内側の段で R が偶然 ±ε の中に入るだけの交点を拾わない (冷却壁のスロート付近で y* の有無が
    場のわずかな違いで切り替わっていた) が、K = 4 は帯の 1 段の揺れを消せず不採用 (plan tooling-nozzle-cooled-wall-ns-settling §4.2b・§4.4)。
    $\kappa=\ln(y^*/y_{ref})$, $y_{ref}=\max(\delta_{in}, f\,r_w)$ を x 方向に移動中央値 (半幅 window) → 3 次 P-spline
    (ノット knot, λ) で平滑化し、$y_b=c\,y_{ref}e^{\tilde\kappa}$。NaN の断面は y_b 側をスプラインで埋める (δ_r 側の穴にしない)。"""
    persist = int(persist)
    if persist < 1:
        raise ValueError(f"persist={persist} は 1 以上")
    x = np.asarray(x, dtype=float); n = len(x)
    d_in = np.asarray(d_in, dtype=float); r_w = np.asarray(r_w, dtype=float)
    y_ref = np.maximum(np.maximum(d_in, 0.0), ref_min_frac * r_w)
    y_star = np.full(n, np.nan)
    for i in range(n):
        y = np.asarray(y_ladders[i], dtype=float); d = np.asarray(d_ladders[i], dtype=float)
        fin = np.isfinite(d) & (d > 0)
        if fin.sum() < 3:
            continue
        ly, dd = np.log(y[fin]), d[fin]
        lq = ly + np.log(edge_ratio)
        m = lq <= ly[-1] + 1e-12
        if m.sum() < 2:
            continue
        Rabs = np.abs(np.interp(lq[m], ly, dd) / dd[m] - 1.0) - eps
        hit = Rabs < 0.0
        if persist > 1:
            # 持続: 段 j から persist 段 (j … j+persist−1) がすべて |R| < ε のときだけ j を候補にする (梯子の外へははみ出さない)。
            # 選んだ j の 1 つ内側の段は必ず |R| ≥ ε (そうでなければ j−1 が先の候補) なので、下の交点の補間は今と同じ 2 段の間で解く
            run = np.zeros_like(hit)
            if len(hit) >= persist:
                run[:len(hit) - persist + 1] = np.lib.stride_tricks.sliding_window_view(hit, persist).all(axis=1)
            hit = run
        j = int(np.argmax(hit)) if np.any(hit) else -1
        if j < 0:
            continue
        L = ly[m]
        if j == 0:
            y_star[i] = float(np.exp(L[0]))
            continue
        lin = None
        a, b = Rabs[j - 1], Rabs[j]
        t = a / (a - b) if a != b else 0.0
        lin = L[j - 1] + t * (L[j] - L[j - 1])
        y_star[i] = float(np.exp(lin))
        if crossing == "quadratic" and len(L) >= 3:
            # 隣接 3 点 (交点を挟む 2 点 + もう 1 点) の 2 次フィットの交点 (位相に鈍感にする; diagnostician 2 回目)
            k0 = j - 1 if j + 1 < len(L) else j - 2
            k0 = max(k0, 0)
            pts = slice(k0, k0 + 3)
            cf = np.polyfit(L[pts] - L[j - 1], Rabs[pts], 2)
            rts = np.roots(cf)
            span = L[j] - L[j - 1]
            ok = [float(np.real(r_)) for r_ in rts if abs(np.imag(r_)) < 1e-12 and -1e-9 <= np.real(r_) <= span + 1e-9]
            if ok:
                y_star[i] = float(np.exp(L[j - 1] + min(ok, key=lambda v: abs(v - (lin - L[j - 1])))))
    kappa = np.log(y_star / y_ref)
    med = np.full(n, np.nan)
    fin = np.isfinite(kappa)
    for i in range(n):
        w = fin & (np.abs(x - x[i]) <= window)
        if w.sum() >= 3:
            med[i] = float(np.median(kappa[w]))
    ks = pspline_uniform(x, med, knot=knot, lam=lam, weights=np.isfinite(med).astype(float))
    # 端の値保持 (P-spline の端を無拘束にしない; diagnostician 2 回目)
    if hold_lo is not None and np.any(x >= hold_lo):
        ks = np.where(x < hold_lo, float(np.interp(hold_lo, x, ks)), ks)
    if hold_hi_from_end is not None:
        xh = float(x.max()) - hold_hi_from_end
        ks = np.where(x > xh, float(np.interp(xh, x, ks)), ks)
    return dict(y_b=c * y_ref * np.exp(ks), y_star=y_star, kappa=kappa, kappa_med=med, kappa_s=ks, y_ref=y_ref)


def deltastar_from_core_matched_euler(ns_run, euler_run, core_frac: float = 0.30,
                                      core_frac_sens=(0.25, 0.35), n_axis_skip: int = 1,
                                      outer_frac: float = 0.25, smooth_lam: float = 1e-3,
                                      gate_core_rms: float = 0.01, gate_sens: float = 0.05,
                                      gate_sens_floor: float = 1e-3, refine: int = 4, out_dir=None,
                                      method: str = "band", k_band: float = 1.5, k_band_sens=(1.125, 2.25),
                                      band_min_frac: float = 0.02, band_width: float = 0.5,
                                      gate_band_rms: float = 0.005, slope_max: float = 0.01,
                                      band_growth: float = 1.25, max_retry: int = 14, fit_from: float = 1.0,
                                      band_smooth_x: float = 0.0, band_select: str = "adaptive",
                                      plateau_tol: float = 0.003, plateau_kmax: int = 9,
                                      edge_eps: float = 0.003, edge_c: float = 1.25, edge_window: float = 3.0,
                                      edge_knot: float = 6.0, edge_step: float = 1.25 ** 0.125, edge_y0_fac: float = 1.2,
                                      edge_min_frac: float = 0.016, edge_ymax_frac: float = 0.5, edge_ref_min_frac: float = 0.0133,
                                      edge_y0_shift: float = 1.0, edge_persist: int = 1, return_ladders: bool = False,
                                      remap_constant_x: bool = False, remap_edge: str = "clamp", ns_res: str | None = None,
                                      yb_sens_factors=None) -> dict:
    r"""**固定 Euler 基準・コア整合**の半径方向等価排除厚 $\delta_r(x)$ を NS 全列で抽出する。

    定義 (plan §4.2–4.4):

    - 対応は**同じ物理 $r$** (設計 $r_t$ 単位)。NS の有効壁 $r_{ref}=r_{w,NS}-\delta_{in}$ は、
      $\delta_{in}$ を「NS 壁 − Euler 壁」で定義するので $r_{ref}=r_{w,E}$ となり、η 写像は恒等になる。
      Euler 壁の外側 ($r>r_{w,E}$) は Euler 壁値で定数外挿。
    - コア $r \le \eta_c r_{w,E}$ の面積重み流量を一致させる倍率 $\alpha$、$q_{ref}=\alpha q_E$。
    - 符号付き欠損 $D = 2\pi\int_0^{r_{w,NS}}(q_{ref}-q_{NS})\,r\,dr$ (クリップなし)。
      $2\pi\int_{r_{eff}}^{r_{w,NS}} q_{ref}\,r\,dr = D$ の根から $\delta_r = r_{w,NS}-r_{eff}$。
    - ゲート (診断; 値は捨てない): コア相対 RMS、コア範囲感度、$D<0$、根なし、欠損のコア分布。

    帯の選び方 `band_select` (method="band" のみ):

    - "adaptive" (既定): `band_local_deficit` の自動選択 (帯内比の変化・残留欠損で 1.25 倍ずつ広げる)。
    - "plateau" (2026-10-04, plan verification-m6-axis-wave-mesh-su2 §4.4): 帯を
      $y_{b,k}=\max(k_{band}\delta_{in}, f_{min} r_w)\,g^k$ ($k=0..$plateau_kmax) と外へ並べて各 k で δ_r を測り、
      **連続 2 段の相対変化がともに plateau_tol 未満になった最初の k** の中央 δ_r,k+1 を採る (感度は δ_r,k / δ_r,k+2)。
      見つからなければ `no_plateau` (hard、前回値保持)。adaptive は帯の比の形で止めるので、境界層の縁にかかった帯を
      通してしまい x 方向に 1.25 倍刻みで跳ぶ (case/45 x≈20 で δ_r −1.8 %)。plateau は「答えが帯の位置に依存しなくなったか」で止める。
      **不採用** (離散 k のまま、no_plateau の穴が平滑化段を壊す; plan §9 2026-10-04)。比較用に残置。
    - "edge" (2026-10-04, diagnostician 推薦 E, plan §5.1 #7a): 帯を細かい梯子 (比 edge_step) で固定して δ_r(y) を測り、
      `edge_band_positions` で縁 y* を連続推定 → x 平滑 → y_b = edge_c × 縁 に帯を固定して測る。感度欄は c×0.8 / c×1.25。
      edge_y0_shift は梯子の位相試験用 (始点の倍率)。return_ladders=True で梯子と縁を戻り値に含める (試験用)。
      edge_persist は縁の交点に要求する持続の段数 (`edge_band_positions` の persist; 既定 1 = 2026-10-13 までの抽出とビット同一。
      > 1 は opt-in、plan tooling-nozzle-cooled-wall-ns-settling §4.2b・§4.4)。

    euler_run: run dir (最新の出力を読む; 今の挙動) か固定のスナップショットの指定 {run, res, sha256} (`resolve_euler_ref`、sha が合わなければ
    ValueError; plan tooling-nozzle-isothermal-wall-chain §5.1 #36e (6))。remap_edge: `_load_structured` の remap_constant_x の端の扱い
    ("clamp" 既定 = 今の挙動 / "missing" = 補間範囲の足りない列を欠測にして抽出しない; §5.1 #36e (4))。Euler の範囲の外の列の判定は
    `euler_cover_range` (§5.1 #36g): "missing" の経路は Euler の端の外側だけに幾何の許容 (端の列の間隔 × 1e−4) を足し、許容の内側で
    端を越える列は端の列の値で評価する (同一の端点の丸めを吸収する許容 = 微小区間の定数の延長; 範囲の不足は今どおり欠測)。
    "clamp" は ±1e−9 のまま。

    帯位置の感度 `yb_sens_factors` (2026-10-12, plan tooling-nozzle-bell-two-mode §5.1 #7 (7b-1)、codex diagnose M2; method="band" のみ):
    同じ壁 (delta_in = 0) では帯の始点が $\max(k_{band}\delta_{in}, f_{min}r_w)=f_{min}r_w$ になり、k_band の感度系列は主系列と
    ビット一致して感度にならない。そこで各断面で**実際に選ばれた帯位置** $y_b$ (`band_y_b`) を倍率 f 倍に固定して測り直し
    (`y_b_fixed`)、`delta_r_sens_yb` (n, 倍率の数) に返す。主系列との相対差の最大が gate_sens (下限 gate_sens_floor) を超えた
    断面は soft の理由 `yb_sens` を付ける (全倍率が NaN なら `yb_sens_nan`)。None (既定) は計算しない (戻り値は従来どおり)。
    ns_res: NS のスナップショット名 (None = 最新; `_load_structured`)。δ_r の時系列用 (流量の帳簿にも同じ名前を渡す)。

    戻り値: dict of arrays (x, r_wall_euler, r_wall_ns, delta_in, alpha, core_rms, core_rms_noaxis,
    core_maxdev, mass_deficit, delta_r_raw, delta_r_smooth, delta_r_sens (n,2), ok, reason) と mdot 帳簿、抽出しなかった NS の列
    `not_extracted` ({n, i, x, reason}; reason は outside_euler (Euler の範囲の外)・remap_missing (端の補間範囲の不足)・no_result
    (帯の計算が結果を返さない)。"missing" の経路では範囲 `euler_cover` と許容の内側で端を越えた列 `edge_tol_cols` も付く)。
    `out_dir` を与えると `delta_r_equiv.csv` / `delta_r_equiv_diag.json` を書く。"""
    import json
    from pathlib import Path
    from scipy.interpolate import make_smoothing_spline
    if yb_sens_factors is not None and method != "band":
        raise ValueError("yb_sens_factors は method='band' でのみ使える (帯位置 y_b の感度)")
    N = _load_structured(ns_run, remap_constant_x=remap_constant_x, remap_edge=remap_edge, res=ns_res)
    eref = resolve_euler_ref(euler_run)
    E = _load_structured(eref["run"], res=eref["res"])
    miss_cols = set(N.get("remap_missing_cols", {}).get("i", []))      # 端の補間範囲の不足で欠測にした列 (remap_edge="missing")
    xE = E["x"][:, 0]
    rwE_col = E["r"][:, -1]
    # Euler の補間できる範囲 (§5.1 #36g、`euler_cover_range`): "missing" の経路は端の外側だけに幾何の許容 (端の列の間隔 × 1e−4) を足し、
    # 許容の内側で端を越える列は euler_q_at (w を [0, 1] に丸める) で端の列の値になる (同一の端点の丸めを吸収する許容 = 微小区間の
    # 定数の延長、§5.1 #36i (v))。"clamp" (既定) は今までどおり ±1e−9
    cover = euler_cover_range(xE, edge_tol=(remap_edge == "missing"))

    def covered(x) -> bool:
        return cover["lo"] <= x <= cover["hi"]

    def euler_q_at(x, r):
        """Euler 場を (x, r) へ: 列間線形 × 列内線形 (r>壁は壁値)。"""
        k = int(np.clip(np.searchsorted(xE, x) - 1, 0, len(xE) - 2))
        w = (x - xE[k]) / max(xE[k + 1] - xE[k], 1e-30)
        w = float(np.clip(w, 0.0, 1.0))
        q0 = np.interp(r, E["r"][k], E["q"][k])          # np.interp は範囲外を端値で埋める
        q1 = np.interp(r, E["r"][k + 1], E["q"][k + 1])
        rw = (1 - w) * rwE_col[k] + w * rwE_col[k + 1]
        return (1 - w) * q0 + w * q1, rw

    def extract_one(i, fc, y_b_fixed=None):
        x = float(N["x"][i, 0]); r0 = N["r"][i]; q0 = N["q"][i]
        if not covered(x) or i in miss_cols:
            return None
        _, rwE = euler_q_at(x, r0[:1])
        if method == "band":
            res = band_local_deficit(r0, q0, lambda rr: euler_q_at(x, rr)[0], float(rwE), k_band=fc,
                                     band_min_frac=band_min_frac, band_width=band_width,
                                     n_axis_skip=n_axis_skip, core_frac=core_frac, refine=refine,
                                     slope_max=slope_max, band_growth=band_growth, max_retry=max_retry, fit_from=fit_from,
                                     y_b_fixed=y_b_fixed)
        else:
            res = core_matched_deficit(r0, q0, lambda rr: euler_q_at(x, rr)[0], float(rwE), core_frac=fc,
                                       n_axis_skip=n_axis_skip, outer_frac=outer_frac, refine=refine)
        if res is None:
            return None
        res["x"] = x
        res["i"] = int(i)
        res["i_col"] = i                                # 帯位置の感度 (yb_sens_factors) で同じ列を測り直すため
        return res

    main_par = k_band if method == "band" else core_frac
    sens_par = list(k_band_sens) if method == "band" else list(core_frac_sens)
    yb_fixed = None
    if method == "band" and band_smooth_x > 0:
        # 1 パス目: 各列の適応 y_b。隣接列で帯が 1.25 倍刻みにトグルすると δ_r が ±1 % ジグザグになる
        # (run_0022 x=90–94 で実測) ので、log y_b を x 方向に移動中央値 (窓 band_smooth_x [r_t]) で平滑化して固定する。
        xs_all = N["x"][:, 0]; yb_all = np.full(len(xs_all), np.nan)
        for i in range(len(xs_all)):
            b0 = extract_one(i, main_par)
            if b0 is not None:
                yb_all[i] = b0.get("band_y_b", np.nan)
        fin = np.isfinite(yb_all)
        lyb = np.log(np.where(fin, yb_all, np.nan))
        sm = np.full_like(lyb, np.nan)
        for i in range(len(xs_all)):
            w = fin & (np.abs(xs_all - xs_all[i]) <= 0.5 * band_smooth_x)
            if w.sum() >= 3:
                sm[i] = np.median(lyb[w])
        yb_fixed = np.exp(sm)
    if band_select not in ("adaptive", "plateau", "edge"):
        raise ValueError(f"band_select={band_select!r} は 'adaptive' / 'plateau' / 'edge'")
    if band_select in ("plateau", "edge") and (method != "band" or band_smooth_x > 0):
        raise ValueError(f"band_select={band_select!r} は method='band' かつ band_smooth_x=0 でのみ使える")
    edge = None
    if band_select == "edge":
        ii, xs_e, din_e, rw_e_l, yl, dl = [], [], [], [], [], []
        for i in range(N["x"].shape[0]):
            x = float(N["x"][i, 0])
            if not covered(x) or i in miss_cols:
                continue
            _, rwE = euler_q_at(x, N["r"][i][:1])
            rwN = float(N["r"][i][-1]); d_in = rwN - float(rwE)
            y0 = max(edge_y0_fac * max(d_in, 0.0), edge_min_frac * rwN) * edge_y0_shift
            ys = []
            y = y0
            while y <= edge_ymax_frac * rwN:
                ys.append(y); y *= edge_step
            ds = []
            for yb in ys:
                rk = extract_one(i, main_par, y_b_fixed=yb)
                ds.append(rk["delta_r"] if (rk is not None and np.isfinite(rk["delta_r"])) else np.nan)
            ii.append(i); xs_e.append(x); din_e.append(d_in); rw_e_l.append(rwN); yl.append(np.array(ys)); dl.append(np.array(ds))
        pos = edge_band_positions(np.array(xs_e), yl, dl, np.array(din_e), np.array(rw_e_l), eps=edge_eps, c=edge_c,
                                  window=edge_window, knot=edge_knot, edge_ratio=1.25, ref_min_frac=edge_ref_min_frac,
                                  persist=edge_persist)
        edge = dict(i=ii, x=np.array(xs_e), y_b=pos["y_b"], y_star=pos["y_star"], kappa_s=pos["kappa_s"],
                    kappa_med=pos["kappa_med"], y_ref=pos["y_ref"], ladders_y=yl, ladders_d=dl)
        edge_map = {i: k for k, i in enumerate(ii)}

    def extract_plateau(i):
        x = float(N["x"][i, 0])
        if not covered(x) or i in miss_cols:
            return None
        _, rwE = euler_q_at(x, N["r"][i][:1])
        rwN = float(N["r"][i][-1]); d_in = rwN - float(rwE)
        y0 = max(k_band * max(d_in, 0.0), band_min_frac * rwN)
        lad, res_k = [], []
        for k in range(plateau_kmax + 1):
            rk = extract_one(i, main_par, y_b_fixed=y0 * band_growth ** k)
            res_k.append(rk)
            lad.append(rk["delta_r"] if (rk is not None and np.isfinite(rk["delta_r"])) else np.nan)
        lad = np.array(lad)
        for k in range(len(lad) - 2):
            v = lad[k:k + 3]
            if np.all(np.isfinite(v)) and v[0] > 0 and abs(v[1] / v[0] - 1) < plateau_tol and abs(v[2] / v[1] - 1) < plateau_tol:
                res = dict(res_k[k + 1]); res["band_retry"] = k + 1
                return res, [float(v[0]), float(v[2])]
        # プラトー無し: 中央段の値を入れ、hard 不合格として前回値を保持させる
        fin = [kk for kk in range(len(lad)) if np.isfinite(lad[kk]) and res_k[kk] is not None]
        if not fin:
            return None
        k = min(fin, key=lambda kk: abs(kk - len(lad) // 2))
        base = dict(res_k[k])
        if base is None:
            return None
        base["reason"] = list(base["reason"]) + ["no_plateau"]
        return base, [np.nan, np.nan]

    rows = []
    for i in range(N["x"].shape[0]):
        if band_select == "edge":
            if i not in edge_map:
                continue
            yb = float(edge["y_b"][edge_map[i]])
            base = extract_one(i, main_par, y_b_fixed=yb)
            if base is None:
                continue
            sens = [extract_one(i, main_par, y_b_fixed=yb * f) for f in (0.8, 1.25)]
            base["delta_r_sens"] = [s_["delta_r"] if s_ else np.nan for s_ in sens]
            base["edge_y_star"] = float(edge["y_star"][edge_map[i]])
            rows.append(base)
            continue
        if band_select == "plateau":
            got = extract_plateau(i)
            if got is None:
                continue
            base, sens_v = got
            base["delta_r_sens"] = sens_v
            rows.append(base)
            continue
        if yb_fixed is not None and np.isfinite(yb_fixed[i]):
            base = extract_one(i, main_par, y_b_fixed=float(yb_fixed[i]))
            sens = [extract_one(i, main_par, y_b_fixed=float(yb_fixed[i]) * f) for f in (0.8, 1.25)]
        else:
            base = extract_one(i, main_par)
            sens = [extract_one(i, fc) for fc in sens_par]
        if base is None:
            continue
        base["delta_r_sens"] = [s["delta_r"] if s else np.nan for s in sens]
        rows.append(base)
    if yb_sens_factors is not None:
        # 帯位置の感度 (codex M2): 実際に選ばれた y_b を倍率で動かして固定し測り直す
        for r_ in rows:
            yb_sel = float(r_.get("band_y_b", np.nan))
            ss = [extract_one(int(r_["i_col"]), main_par, y_b_fixed=yb_sel * float(f)) if np.isfinite(yb_sel) else None
                  for f in yb_sens_factors]
            r_["delta_r_sens_yb"] = [s_["delta_r"] if s_ else np.nan for s_ in ss]
            r_["_yb_sens_rows"] = ss                    # 感度帯の診断 (停止条件の量) を返すため
    x = np.array([r_["x"] for r_ in rows])
    draw = np.array([r_["delta_r"] for r_ in rows])
    sens = np.array([r_["delta_r_sens"] for r_ in rows])
    sens_yb = (np.array([r_["delta_r_sens_yb"] for r_ in rows], dtype=float).reshape(len(rows), len(yb_sens_factors))
               if yb_sens_factors is not None else None)
    # negative_deficit は 2026-09-12 から soft (符号付き δ_r を採用: 冷却壁では物理)。hard は根なし/NaN と、帯局所の探索域の参照が
    # 非正 (nonpos_ref)・探索区間の F が非単調 (nonmono_zone) (2026-10-13, plan tooling-nozzle-deltastar-band-root-zone §4.1)
    HARD = ("no_root", "nan", "no_plateau", "nonpos_ref", "nonmono_zone")
    ok = np.ones(len(rows), bool); hard_ok = np.ones(len(rows), bool); reasons = []
    for k, r_ in enumerate(rows):
        rs = list(r_["reason"])
        if method == "band":
            if r_.get("band_rms", 0.0) > gate_band_rms: rs.append("band_fit")
            if np.isfinite(r_.get("band_deficit_share", np.nan)) and abs(r_["band_deficit_share"]) > 0.1: rs.append("band_not_clean")
        elif r_["core_rms_noaxis"] > gate_core_rms: rs.append("core_shape")
        if np.isfinite(draw[k]) and draw[k] > 0 and np.all(np.isfinite(sens[k])):
            if np.max(np.abs(sens[k] - draw[k])) > max(gate_sens * draw[k], gate_sens_floor): rs.append("core_frac_sens")
        if sens_yb is not None and np.isfinite(draw[k]) and draw[k] > 0:
            fin_yb = np.isfinite(sens_yb[k])
            if not fin_yb.any():
                rs.append("yb_sens_nan")
            elif np.max(np.abs(sens_yb[k][fin_yb] - draw[k])) > max(gate_sens * draw[k], gate_sens_floor):
                rs.append("yb_sens")
        if np.isfinite(r_["outer_share"]) and r_["outer_share"] < 0.8: rs.append("deficit_not_near_wall")
        if not np.isfinite(draw[k]): rs.append("nan")
        ok[k] = len(rs) == 0
        hard_ok[k] = not any(h in rs for h in HARD)
        reasons.append(",".join(rs))
    # 平滑化: hard-ok の生値だけを使い、範囲外は端値 (外挿しない)。soft 不合格の値は使う (plan §4.4)。
    dsm = np.full_like(draw, np.nan)
    if hard_ok.sum() > 10:
        xs_, ds_ = x[hard_ok], draw[hard_ok]
        spl = make_smoothing_spline(xs_, ds_, lam=smooth_lam)
        dsm = spl(np.clip(x, xs_[0], xs_[-1]))
    # 壁更新に渡す値: hard-ok なら平滑化値、hard 不合格なら NaN (呼び出し側が前回値を保持する)
    duse = np.where(hard_ok, dsm, np.nan)
    out = dict(x=x, r_wall_euler=np.array([r_["r_wall_euler"] for r_ in rows]),
               r_wall_ns=np.array([r_["r_wall_ns"] for r_ in rows]),
               delta_in=np.array([r_["delta_in"] for r_ in rows]),
               alpha=np.array([r_["alpha"] for r_ in rows]),
               core_rms=np.array([r_["core_rms"] for r_ in rows]),
               core_rms_noaxis=np.array([r_["core_rms_noaxis"] for r_ in rows]),
               core_maxdev=np.array([r_["core_maxdev"] for r_ in rows]),
               mass_deficit=np.array([r_["mass_deficit"] for r_ in rows]),
               outer_share=np.array([r_["outer_share"] for r_ in rows]),
               band_y_b=np.array([r_.get("band_y_b", np.nan) for r_ in rows]),
               band_slope=np.array([r_.get("band_slope", np.nan) for r_ in rows]),
               band_rms=np.array([r_.get("band_rms", np.nan) for r_ in rows]),
               band_deficit_share=np.array([r_.get("band_deficit_share", np.nan) for r_ in rows]),
               band_retry=np.array([r_.get("band_retry", 0) for r_ in rows]),
               delta_r_raw=draw, delta_r_smooth=dsm, delta_r_use=duse, delta_r_sens=sens,
               ok=ok, hard_ok=hard_ok, reason=np.array(reasons),
               settings=dict(method=method, k_band=k_band, k_band_sens=list(k_band_sens), band_min_frac=band_min_frac,
                             band_width=band_width, slope_max=slope_max, band_growth=band_growth, max_retry=max_retry, fit_from=fit_from,
                             band_smooth_x=band_smooth_x, band_select=band_select,
                             plateau_tol=plateau_tol, plateau_kmax=plateau_kmax,
                             edge=(dict(eps=edge_eps, c=edge_c, window=edge_window, knot=edge_knot, step=edge_step,
                                        y0_fac=edge_y0_fac, min_frac=edge_min_frac, ymax_frac=edge_ymax_frac,
                                        ref_min_frac=edge_ref_min_frac, y0_shift=edge_y0_shift, persist=int(edge_persist))
                                   if band_select == "edge" else None),
                             gate_band_rms=gate_band_rms,
                             core_frac=core_frac, core_frac_sens=list(core_frac_sens), n_axis_skip=n_axis_skip,
                             outer_frac=outer_frac, smooth_lam=smooth_lam, gate_core_rms=gate_core_rms,
                             gate_sens=gate_sens, gate_sens_floor=gate_sens_floor, refine=refine,
                             ns_res=N["res"], euler_res=E["res"],
                             **({"remap_constant_x": True, "remap_max_column_dx": N["remap_max_column_dx"],
                                 "remap_clamped_points": N["remap_clamped_points"]} if remap_constant_x else {}),
                             **({"remap_edge": remap_edge, "remap_missing_cols": N.get("remap_missing_cols")}
                                if remap_constant_x and remap_edge != "clamp" else {}),
                             **({"euler_ref_pinned": eref["record"]} if eref["pinned"] else {}),
                             ns_run=str(ns_run), euler_run=(str(eref["run"]) if eref["pinned"] else str(euler_run))))
    if sens_yb is not None:                         # 既定 (None) では鍵を増やさない (戻り値・diag JSON を従来どおりに保つ)
        out["delta_r_sens_yb"] = sens_yb
        out["settings"]["yb_sens_factors"] = [float(f) for f in yb_sens_factors]
        # 帯の診断 (2026-10-12, plan tooling-nozzle-bell-two-mode §5.1 #7 の edge の A/B): 主系列の帯の比の切片 a (適応の停止条件
        # a ≥ 0.95 α に使う量) と、感度帯 (n, 倍率の数) ごとの停止条件の量。固定帯では再試行しないので、条件の違反は理由に残らない —
        # ここで量を返し、判定は呼び出し側 (`bell_deltastar`) が行う
        out["band_ratio_a"] = np.array([r_.get("band_ratio_a", np.nan) for r_ in rows])
        for key in ("band_y_b", "band_slope", "band_rms", "band_deficit_share", "band_ratio_a", "alpha", "mass_deficit"):
            out["yb_sens_" + key] = np.array([[(s_[key] if s_ is not None else np.nan) for s_ in r_["_yb_sens_rows"]] for r_ in rows],
                                             dtype=float).reshape(len(rows), len(yb_sens_factors))
        out["yb_sens_reason"] = np.array([[(",".join(s_["reason"]) if s_ is not None else "none") for s_ in r_["_yb_sens_rows"]]
                                          for r_ in rows]).reshape(len(rows), len(yb_sens_factors))
    # 抽出しなかった NS の列 (§5.1 #36e: 受け渡しの関門が全域の被覆・使用域の欠測を判定する)
    done = {r_["i"] for r_ in rows}
    xs_all = N["x"][:, 0]
    ne = [i for i in range(len(xs_all)) if i not in done]
    why = ["remap_missing" if i in miss_cols else ("outside_euler" if not covered(float(xs_all[i])) else "no_result") for i in ne]
    out["not_extracted"] = {"n": len(ne), "n_cols_ns": int(len(xs_all)), "i": [int(i) for i in ne],
                            "x": [float(xs_all[i]) for i in ne], "reason": why}
    if cover["edge_tol"]:
        # 幾何の許容 (§5.1 #36g) の内側で Euler の端を越え、端の列の値で評価した NS の列 (記録; 範囲の不足を隠していないことを後で見る)
        et = [i for i in range(len(xs_all)) if covered(float(xs_all[i])) and not (cover["x0"] <= float(xs_all[i]) <= cover["x1"])]
        rec_et = {"n": len(et), "i": [int(i) for i in et], "x": [float(xs_all[i]) for i in et],
                  "beyond": [float(max(cover["x0"] - xs_all[i], xs_all[i] - cover["x1"])) for i in et]}
        out["not_extracted"].update(euler_cover=cover, edge_tol_cols=rec_et)
        out["settings"].update(euler_cover=cover, edge_tol_cols=rec_et)
    # 流量の帳簿も抽出と同じ断面で (2026-10-15, plan tooling-nozzle-card-local-mesh §5.1 #3i): 壁法線の近壁層の格子では列が一定 x でなく、
    # 列のままの 2π∫ρU_x r dr は断面の流量でない (case/71 run_0005_ns で −0.96 %、散らばり 2.1 %)
    out["massflow"] = massflow_ratio(ns_run, euler_run, ns_res=ns_res, remap_constant_x=remap_constant_x)
    if edge is not None:
        out["edge_y_star"] = np.array([r_.get("edge_y_star", np.nan) for r_ in rows])
        if return_ladders:
            out["edge_detail"] = edge
    if out_dir is not None:
        od = Path(out_dir); od.mkdir(parents=True, exist_ok=True)
        np.savetxt(od / "delta_r_equiv.csv",
                   np.c_[x, draw, dsm, duse, out["delta_in"], out["alpha"], out["core_rms_noaxis"], sens,
                         out["band_y_b"], out["band_slope"], out["band_rms"], out["band_deficit_share"],
                         ok.astype(int), hard_ok.astype(int)],
                   delimiter=",", comments="",
                   header="x_rt,delta_r_raw,delta_r_smooth,delta_r_use,delta_in,alpha,core_rms_noaxis,delta_r_sens1,delta_r_sens2,band_y_b,band_slope,band_rms,band_deficit_share,ok,hard_ok")
        if edge is not None:
            np.savetxt(od / "delta_r_edge.csv", np.c_[edge["x"], edge["y_ref"], edge["y_star"], edge["kappa_med"], edge["kappa_s"], edge["y_b"]],
                       delimiter=",", comments="", header="x_rt,y_ref,y_star,kappa_med,kappa_s,y_b")
        diag = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in out.items() if k != "edge_detail"}
        (od / "delta_r_equiv_diag.json").write_text(json.dumps(diag))
    return out


# --- δ_r(x) の 5 次 P-spline 平滑化 (壁曲率を滑らかにする; 2026-09-04 ユーザ指摘) ----------------------------
def smooth_delta_quintic(x, d, weights=None, knot_spacing: float = 2.0, lam: float = 1.0,
                         x_stretch: float = 4.0, penalty_order: int = 3, positive: bool = True):
    r"""抽出した $\delta_r(x)$ を **5 次 B-spline + 3 階差分ペナルティ (P-spline)** で平滑化する。

    壁は $r_{inv}(x)+\delta_r(x)$ を補間 5 次 B-spline で通すので、δ_r の 2 階微分のノイズがそのまま壁曲率に乗る。
    3 次平滑化スプラインは曲率まで滑らかにしないため、$\int(\delta_r''')^2$ に相当する係数 3 階差分ペナルティで
    **曲率が滑らかな** δ_r を作る (ユーザ提案)。ノットは $\xi=\sqrt{x-x_0+x_s}$ で等間隔 (スロート側で密:
    δ_r は x=0→8 で 20 倍に育つ)。weights=0 の点 (hard 不合格) は無視。戻り: (callable δ_s(x), 診断 dict)。"""
    from scipy.interpolate import BSpline
    x = np.asarray(x, dtype=float); d = np.asarray(d, dtype=float)
    w = np.ones_like(d) if weights is None else np.asarray(weights, dtype=float)
    ok = np.isfinite(d) & (w > 0)
    x0 = float(x[ok].min())
    xi = lambda xx: np.sqrt(np.maximum(xx - x0, 0.0) + x_stretch)
    k = 5
    xi_all = xi(x[ok]); xi_lo, xi_hi = float(xi_all.min()), float(xi_all.max())
    # ノット間隔: 物理 x で knot_spacing 相当 (ξ 空間では dξ = dx/(2ξ) → 中央値で換算)
    dxi = knot_spacing / (2.0 * np.median(xi_all))
    n_int = max(int(np.ceil((xi_hi - xi_lo) / dxi)) - 1, 4)
    inner = np.linspace(xi_lo, xi_hi, n_int + 2)[1:-1]
    t = np.concatenate([[xi_lo] * (k + 1), inner, [xi_hi] * (k + 1)])
    n_cp = len(t) - k - 1
    B = BSpline.design_matrix(xi_all, t, k).toarray()
    D = np.diff(np.eye(n_cp), n=penalty_order, axis=0)
    W = w[ok][:, None]
    A = B.T @ (W * B) + lam * (D.T @ D)
    c = np.linalg.solve(A, B.T @ (w[ok] * d[ok]))
    spl = BSpline(t, c, k, extrapolate=False)

    def f(xx):
        xx = np.asarray(xx, dtype=float)
        q = np.clip(xi(xx), xi_lo, xi_hi)
        v = spl(q)
        return np.maximum(v, 0.0) if positive else v

    fit = f(x[ok])
    rel = (fit - d[ok]) / np.maximum(np.abs(d[ok]), 1e-6)
    diag = dict(n_cp=int(n_cp), knot_spacing=knot_spacing, lam=lam, penalty_order=penalty_order,
                resid_rel_rms=float(np.sqrt(np.mean(rel ** 2))), resid_abs_max=float(np.max(np.abs(fit - d[ok]))),
                resid_rel_p95=float(np.percentile(np.abs(rel), 95)))
    return f, diag
