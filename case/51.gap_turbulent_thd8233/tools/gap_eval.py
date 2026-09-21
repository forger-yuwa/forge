#!/usr/bin/env python3
"""case/51 (NASA TN D-8233) T2' 評価 — すきま壁の $h/h_{fp}$ を深さ別に出し、3 分類を暫定判定する。

plan `plans/active/case-hypersonic-gap-heating-validation.md` §5.1 #49 の抽出部。
**合否基準は自分で決めない** — `case/51/acceptance.json` (事前登録) と `case/51/tolerances.json`
から読む。ここで実装しているのは「実測との突き合わせ」と「区間の一部を使った暫定分類」だけで、
分類の確定 (case/49 への受け渡し) には下の「未実装の区間」が埋まっている必要がある。

幾何 (`geometry.json`):
  すきまは x 方向に幅 W = 2.29 mm、y 方向に深さ D = 45.72 mm。開口は y = 0。
  **後壁 (downstream / せん断層が当たる側) が x = 0、前壁 (upstream) が x = -W**。
  床は y = -D。physID 6 (`gap`) はこの 3 面 (gen_mesh.py の Curve 17/18/19)。
  **深さの無次元化は z/W = (開口からの深さ)/W** で、ここでの「z」は深さ方向の記号であって
  座標軸 z ではない (y = -z_depth)。判定深さ z/W = 1.4 / 2.2 / 3.0 は y = -3.206 / -5.038 / -6.870 mm。

量 (`conditions.json` の `reference`):
  h = q_w / (T_aw - T_w)、T_aw = 995.486 K (`derived.json` の Re' 1.47e6 系列)、T_w = 壁の Ts。
  分母 h_fp = C Re'^0.69、C = 7.95e-5。**Re' は `--re-eff` で受ける** (既定は台帳の 1.47e6/m)。
  forge の `viscMethod: 1` が Sutherland 定数直書きで実効 Re' がずれる件 (plan §5.1 #48 の未決) を
  どちらでも評価できるようにするため。使った Re' と h_fp は必ず出力に書く。

**符号規約 (重要)**:
  forge の `/VALUE/qwall` は **壁 → 流体が正** (`viscousFlux_d.cu`:654 の診断出力)。
  本ツールはこれを反転し、**流体 → 壁に入る向きを正** として q_w を扱う。
  冷却壁 (T_w 300 K < T_aw 995 K) なので q_w > 0 が正常。出力・CSV の列名にも明記する。

usage:
  python3 tools/gap_eval.py run_0001_t2p_re147_tw300
  python3 tools/gap_eval.py run_0001_t2p_re147_tw300 --series-csv
  python3 tools/gap_eval.py <任意の run> --probe-schema        # 読むデータセットの存在確認だけ
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

# numpy 2.0 で np.trapz が削除された (AWS は numpy 2.x)。名前だけの差なので薄く吸収する。
_trapz = getattr(np, "trapezoid", None) or np.trapz
import h5py

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = CASE.parents[1]
QS_TOOL = ROOT / "solver_density_cuda" / "tools" / "check_quasisteady.py"

# 読むデータセット。無ければ「何が無いか」を報告して止まる (別の量に差し替えない)。
REQUIRED = ["/MESH/COORD", "/VALUE/qwall", "/VALUE/Ts"]
OPTIONAL = ["/VALUE/Tt", "/VALUE/utau", "/VALUE/ypls", "/VALUE/twall_x"]

# acceptance.json `gates[T2p-C].classes` の文言から転記した閾値。
# 文言そのものは実行時に印字して突き合わせられるようにしてある (勝手に緩めないこと)。
ALPHA_MAX = 1.0 / 10.0          # alpha: 上端でも forge/実測 < 1/10
BETA_BAND = (1.0 / 3.0, 3.0)    # beta : 比が [1/3, 3]
ALPHA_DEPTHS = (2.2, 3.0)
BETA_DEPTHS = (1.4, 2.2)
QUAL_ZW_MAX = 4.0               # T2p-Q: z/W < 4 で 後壁 h > 前壁 h

# 時系列 (--series-csv) に出す深さ。判定深さ + 参考 0.7。
SERIES_DEPTHS = [0.7, 1.4, 2.2, 3.0]


# ---------------------------------------------------------------- 入出力まわり
def resolve_run(arg):
    p = Path(arg)
    if p.is_dir():
        return p.resolve()
    q = CASE / arg
    if q.is_dir():
        return q.resolve()
    raise SystemExit(f"REFUSED: run ディレクトリが無い: {arg}")


def wall_dumps(rd, pid):
    """`res_<bcond 名>_<physID>_<step>.h5` を step でひいた dict。"""
    out = {}
    for f in rd.glob(f"res_*_{pid}_*.h5"):
        try:
            out[int(f.stem.rsplit("_", 1)[1])] = f
        except ValueError:
            continue
    return dict(sorted(out.items()))


def probe_schema(rd, pid):
    """読むデータセットが実在するかだけを確認する (値は評価しない・何も書かない)。"""
    print(f"--probe-schema: {rd}")
    nan = sorted(rd.glob("res_nan_*.h5"))
    if nan:
        print(f"  注意: res_nan_*.h5 が {len(nan)} 個ある (発散履歴)")
    dumps = wall_dumps(rd, pid)
    if not dumps:
        print(f"  MISSING: res_*_{pid}_*.h5 が 1 つも無い")
        print("  → bcondConfig の physID %d の bcond に `outputHDFflg: 1` が要る" % pid)
        print("     (case/51 は tools/make_case.py の BC テンプレートで gap: outputHDFflg: 1 済み)")
        return 1
    step, f = max(dumps.items())
    print(f"  壁ダンプ {len(dumps)} 枚 (step {min(dumps)}..{max(dumps)})、点検対象 = {f.name}")
    missing = []
    with h5py.File(f, "r") as h:
        for name in REQUIRED:
            if name in h:
                d = h[name]
                print(f"  OK       {name:18s} shape={d.shape} dtype={d.dtype}")
            else:
                missing.append(name)
                print(f"  MISSING  {name}")
        for name in OPTIONAL:
            print(f"  {'(任意) あり' if name in h else '(任意) 無し'}  {name}")
        if "/VALUE/qwall" in h:
            q = h["/VALUE/qwall"][:].astype(float)
            nz = int(np.count_nonzero(q))
            print(f"  qwall: 非ゼロ {nz}/{q.size} 点、範囲 [{q.min():.4g}, {q.max():.4g}] "
                  "(壁→流体が正なので冷却壁では負が正常)")
            if nz == 0:
                print("  → 全点ゼロ。低 Re 壁 (wallTreatmentSST: 0) の qwall 診断出力が"
                      " 無いビルドか、壁が断熱。`viscousFlux_d.cu`:654 を含む版で回すこと")
                missing.append("/VALUE/qwall (全点ゼロ)")
    if missing:
        print("\n判定: NG — " + ", ".join(missing))
        print("case/51 の run で要るもの: bcondConfig の `gap` (physID 6) と `plate` (physID 4) を"
              " `outputHDFflg: 1` にする。**`output.extraFields` への追加は不要**"
              " (qwall/Ts は境界ダンプ固有の量で、`output.level`/`extraFields` は体積場の話)。")
        return 1
    print("\n判定: OK — gap_eval.py が読む 3 データセットは揃っている")
    return 0


# ---------------------------------------------------------------- 台帳の読み込み
def load_registry(re_eff_cli, t_aw_cli, run_dir):
    acc = json.loads((CASE / "acceptance.json").read_text(encoding="utf-8"))
    tol = json.loads((CASE / "tolerances.json").read_text(encoding="utf-8"))
    geom = json.loads((CASE / "geometry.json").read_text(encoding="utf-8"))
    cond = json.loads((CASE / "conditions.json").read_text(encoding="utf-8"))
    der = json.loads((CASE / "derived.json").read_text(encoding="utf-8"))
    setup = {}
    p = run_dir / "case_setup.json"
    if p.exists():
        setup = json.loads(p.read_text(encoding="utf-8"))

    re_nom = float(setup.get("re_m", der["series"][0]["Re_m"]))
    re_eff = float(re_eff_cli) if re_eff_cli else re_nom
    C = float(cond["reference"]["C_digitized"])
    hfp = C * re_eff ** 0.69
    t_aw = float(t_aw_cli) if t_aw_cli else float(
        setup.get("T_aw", der["series"][0]["T_aw"]))
    return dict(acc=acc, tol=tol, geom=geom, cond=cond, der=der, setup=setup,
                re_nom=re_nom, re_eff=re_eff, C=C, hfp=hfp, t_aw=t_aw)


def tol_item(tol, tid):
    for it in tol["items"]:
        if it["id"] == tid:
            return it
    raise SystemExit(f"REFUSED: tolerances.json に id={tid} が無い")


def ref_measured(csv_path, t_ratio=1.00):
    """ref CSV から {wall: {z/W: h/h_fp}} を作る。solid=1 (黒塗り) は除外。"""
    out = {"down": {}, "up": {}}
    for line in csv_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("T_ratio"):
            continue
        tr, wall, zw, hr, solid = line.split(",")
        if abs(float(tr) - t_ratio) > 1e-9 or int(solid) != 0:
            continue
        out[wall][round(float(zw), 3)] = float(hr)
    return out


# ---------------------------------------------------------------- 壁データの抽出
def _dedupe(d, v):
    """同一深さ (角の共有ノード等) を平均して単調増加にする。"""
    o = np.argsort(d, kind="stable")
    d, v = d[o], v[o]
    du, idx = np.unique(np.round(d, 12), return_inverse=True)
    vu = np.zeros_like(du)
    np.add.at(vu, idx, v)
    cnt = np.bincount(idx, minlength=du.size)
    return du, vu / cnt


def read_wall(path, W, D, t_aw, hfp, tw_override=None, tol=1e-6):
    """壁ダンプ 1 枚を後壁 / 前壁 / 床に分けて返す。

    **符号**: `/VALUE/qwall` は壁→流体が正なので反転し、**流体→壁が正**で持つ。
    """
    with h5py.File(path, "r") as h:
        missing = [n for n in REQUIRED if n not in h]
        if missing:
            raise SystemExit(
                f"REFUSED: {path.name} に {', '.join(missing)} が無い。"
                " 推測で別の量に差し替えない。bcond の outputHDFflg と forge の版を確認すること")
        c = h["/MESH/COORD"][:].reshape(-1, 3).astype(float)
        qw = -h["/VALUE/qwall"][:].astype(float)     # 壁→流体正 → 流体→壁正
        ts = h["/VALUE/Ts"][:].astype(float)
    if not np.count_nonzero(qw):
        raise SystemExit(
            f"REFUSED: {path.name} の qwall が全点ゼロ。低 Re 壁の qwall 診断出力が無い版か"
            " 断熱壁。別の量で代用しない")
    tw = np.full_like(ts, float(tw_override)) if tw_override else ts
    dT = t_aw - tw
    if np.any(np.abs(dT) < 1e-9):
        raise SystemExit("REFUSED: T_aw - T_w が 0 の点がある (h が定義できない)")
    h_conv = qw / dT

    x, y = c[:, 0], c[:, 1]
    m_down = (np.abs(x) < tol) & (y < tol)                 # 後壁 = x = 0
    m_up = (np.abs(x + W) < tol) & (y < tol)               # 前壁 = x = -W
    m_floor = np.abs(y + D) < tol                          # 床  = y = -D
    walls = {}
    for key, m, in (("down", m_down), ("up", m_up)):
        if not m.any():
            raise SystemExit(f"REFUSED: {key} 壁の節点が拾えない (W/D/座標規約を確認)")
        d, q = _dedupe(-y[m], qw[m])
        _, hh = _dedupe(-y[m], h_conv[m])
        _, tt = _dedupe(-y[m], tw[m])
        walls[key] = dict(depth=d, q=q, h=hh, tw=tt, ratio=hh / hfp, n=int(m.sum()),
                          nzero=int(np.count_nonzero(q == 0.0)))
    if m_floor.any():
        o = np.argsort(x[m_floor])
        walls["floor"] = dict(x=x[m_floor][o], q=qw[m_floor][o], h=h_conv[m_floor][o],
                              tw=tw[m_floor][o], ratio=h_conv[m_floor][o] / hfp,
                              n=int(m_floor.sum()),
                              nzero=int(np.count_nonzero(qw[m_floor] == 0.0)))
    return walls


def at_zw(w, zw, W):
    """z/W の位置の h/h_fp (線形補間)。範囲外は NaN。"""
    d = zw * W
    dep = w["depth"]
    if d < dep[0] - 1e-12 or d > dep[-1] + 1e-12:
        return float("nan")
    return float(np.interp(d, dep, w["ratio"]))


def band_mean(w, W, zw_max=QUAL_ZW_MAX):
    """z/W <= zw_max の h/h_fp の深さ平均 (台形則)。代表値として時系列に出す。"""
    m = w["depth"] <= zw_max * W + 1e-12
    if m.sum() < 2:
        return float("nan")
    return float(_trapz(w["ratio"][m], w["depth"][m]) /
                 (w["depth"][m][-1] - w["depth"][m][0]))


# ---------------------------------------------------------------- 区間 (伝播)
def ratio_band(r_forge, r_meas, denom_hi, dig_abs):
    """比 forge/実測 の取りうる範囲。

    tolerances.json の 2 項目だけを入れる:
      - `denominator` [1.0, 1.8]: forge も実測も同じ相関分母で割るので**一次では相殺**する
        (= 下端 1.0)。相関自体が 3 次元化したトンネル壁 BL 由来なので**上端側だけ残す**
        (= 上端 1.8)。tolerances.json の `direction` の文言どおりの入れ方。
      - `digitize` 絶対 ±0.03: 実測値の読み取り。比を最大にする側は実測 -0.03。
    acceptance.json の `band_rule` は**上端で判定する** (= 欠損が最も縮む側に倒す)。
    """
    nom = r_forge / r_meas
    lo = r_forge / (r_meas + dig_abs)
    den = r_meas - dig_abs
    hi = float("inf") if den <= 0 else r_forge * denom_hi / den
    return lo, nom, hi


# ---------------------------------------------------------------- 表示
def print_wall_table(name, w, W, every):
    print(f"\n--- {name} (点 {w['n']}, q_w == 0.0 の点 {w['nzero']}) ---")
    if "depth" in w:
        print(f"{'深さ[mm]':>9} {'z/W':>7} {'T_w[K]':>8} {'q_w[W/m2]':>12} "
              f"{'h[W/m2K]':>10} {'h/h_fp':>9}")
        idx = list(range(0, len(w["depth"]), max(1, every)))
        if idx[-1] != len(w["depth"]) - 1:
            idx.append(len(w["depth"]) - 1)
        for i in idx:
            print(f"{w['depth'][i]*1e3:9.3f} {w['depth'][i]/W:7.3f} {w['tw'][i]:8.2f} "
                  f"{w['q'][i]:12.4g} {w['h'][i]:10.4g} {w['ratio'][i]:9.4f}")
    else:
        print(f"{'x[mm]':>9} {'深さ[mm]':>9} {'T_w[K]':>8} {'q_w[W/m2]':>12} "
              f"{'h[W/m2K]':>10} {'h/h_fp':>9}   (床は深さ一定)")
        idx = list(range(0, len(w["x"]), max(1, every)))
        if idx[-1] != len(w["x"]) - 1:
            idx.append(len(w["x"]) - 1)
        for i in idx:
            print(f"{w['x'][i]*1e3:9.4f} {'-':>9} {w['tw'][i]:8.2f} "
                  f"{w['q'][i]:12.4g} {w['h'][i]:10.4g} {w['ratio'][i]:9.4f}")


def fmt(v, w=9, p=4):
    return f"{'nan':>{w}}" if not np.isfinite(v) else (
        f"{'inf':>{w}}" if np.isinf(v) else f"{v:{w}.{p}f}")


# ---------------------------------------------------------------- 時系列
def write_series(dumps, W, D, reg, args, out_path):
    cols = ["step"] + [f"down_zW{z:g}".replace(".", "p") for z in SERIES_DEPTHS] \
        + ["down_mean_zW0_4", "up_mean_zW0_4"]
    rows = []
    for step, f in dumps.items():
        w = read_wall(f, W, D, reg["t_aw"], reg["hfp"], args.tw)
        rows.append([step] + [at_zw(w["down"], z, W) for z in SERIES_DEPTHS]
                    + [band_mean(w["down"], W), band_mean(w["up"], W)])
    a = np.array(rows, float)
    np.savetxt(out_path, a, delimiter=",", header=",".join(cols), comments="")
    print(f"\n時系列 CSV: {out_path}  ({len(rows)} 枚, step {a[0,0]:.0f}..{a[-1,0]:.0f})")
    print("  列 = " + ", ".join(cols[1:]) + "   (すべて h/h_fp。q_w は流体→壁が正)")
    print(f"  判定はこれを check_quasisteady.py に渡す:\n"
          f"    python3 {QS_TOOL} --series-csv {out_path} \\\n"
          f"        --series-cols {','.join(cols[1:])}")
    if len(rows) < 4:
        print("  注意: スナップショットが 4 枚未満なので check_quasisteady は判定できない")


# ---------------------------------------------------------------- 本体
def main():
    ap = argparse.ArgumentParser(
        description="case/51 すきま壁の h/h_fp 抽出と T2' の暫定分類 "
                    "(基準は acceptance.json / tolerances.json から読む)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="符号規約: forge の /VALUE/qwall は壁→流体が正。本ツールは反転して"
               "**流体→壁が正** (冷却壁で q_w > 0) で扱う。")
    ap.add_argument("run", help="run ディレクトリ (case/51 からの相対名でも絶対パスでも可)")
    ap.add_argument("--step", type=int, default=None, help="評価する step (既定は最終)")
    ap.add_argument("--phys-id", type=int, default=6, help="すきま壁の physID (既定 6 = gap)")
    ap.add_argument("--re-eff", type=float, default=None,
                    help="分母 h_fp = C Re'^0.69 に使う実効 Re' [1/m]。"
                         "既定は case_setup.json / derived.json の 1.47e6")
    ap.add_argument("--t-aw", type=float, default=None,
                    help="T_aw [K] (既定は case_setup.json / derived.json の 995.486)")
    ap.add_argument("--tw", type=float, default=None,
                    help="壁温 T_w [K] を上書き (既定は壁ダンプの /VALUE/Ts)")
    ap.add_argument("--t-ratio", type=float, default=1.00,
                    help="突き合わせる実測系列の T_surf/T_gap (既定 1.00)")
    ap.add_argument("--print-every", type=int, default=8, help="壁の表を何点おきに印字するか")
    ap.add_argument("--series-csv", nargs="?", const="AUTO", default=None,
                    help="全スナップショットの時系列 CSV を書く (既定の出力先は <run>/gap_series.csv)")
    ap.add_argument("--probe-schema", action="store_true",
                    help="読むデータセットの存在確認だけして終わる (値は評価しない・何も書かない)")
    a = ap.parse_args()

    rd = resolve_run(a.run)
    if a.probe_schema:
        return probe_schema(rd, a.phys_id)

    nan_files = sorted(rd.glob("res_nan_*.h5"))
    if nan_files:
        raise SystemExit(f"REFUSED: {rd.name} に res_nan_*.h5 が {len(nan_files)} 個ある (発散)")
    dumps = wall_dumps(rd, a.phys_id)
    if not dumps:
        raise SystemExit(
            f"REFUSED: {rd}/res_*_{a.phys_id}_*.h5 が無い。"
            " bcondConfig の physID %d を outputHDFflg: 1 にして回すこと" % a.phys_id)
    step = a.step if a.step is not None else max(dumps)
    if step not in dumps:
        raise SystemExit(f"REFUSED: step {step} の壁ダンプが無い (あるのは {list(dumps)})")

    reg = load_registry(a.re_eff, a.t_aw, rd)
    W = reg["geom"]["gap"]["width"] * 1e-3
    D = reg["geom"]["gap"]["depth"] * 1e-3
    walls = read_wall(dumps[step], W, D, reg["t_aw"], reg["hfp"], a.tw)

    print(f"run   : {rd}")
    print(f"壁    : {dumps[step].name}  (physID {a.phys_id})")
    print("読んだデータセット: /MESH/COORD, /VALUE/qwall, /VALUE/Ts")
    print("符号規約: /VALUE/qwall は**壁→流体が正**。本ツールは反転して"
          "**流体→壁が正** (冷却壁で q_w > 0) で扱う")
    print(f"幾何  : W = {W*1e3:.3f} mm, D = {D*1e3:.3f} mm, D/W = {D/W:.2f}、"
          f"後壁 x = 0 / 前壁 x = -W / 床 y = -D")
    print(f"分母  : h_fp = {reg['C']:.4g} * Re'^0.69 = {reg['hfp']:.4f} W/(m^2 K)"
          f"   [Re' = {reg['re_eff']:.4g} /m"
          + ("" if abs(reg["re_eff"] - reg["re_nom"]) < 1e-6
             else f"  ← --re-eff 指定。台帳値は {reg['re_nom']:.4g}") + "]")
    tws = walls["down"]["tw"]
    print(f"駆動差: T_aw = {reg['t_aw']:.3f} K, T_w = {tws.min():.2f}..{tws.max():.2f} K"
          f"  → T_aw - T_w = {reg['t_aw']-tws.mean():.2f} K")

    for key, label in (("down", "後壁 (downstream, x = 0, せん断層が当たる側)"),
                       ("up", "前壁 (upstream, x = -W)"),
                       ("floor", "床 (y = -D)")):
        if key in walls:
            print_wall_table(label, walls[key], W, a.print_every)

    # ---- 判定深さ ------------------------------------------------------
    meas = ref_measured(CASE / "ref" / "th76_fig5a_re147.csv", a.t_ratio)
    accC = [g for g in reg["acc"]["gates"] if g["id"] == "T2p-C"][0]
    reg_depths = [float(z) for z in accC["depths_zW"]]
    reg_ref = {float(k): float(v) for k, v in accC["reference_values"].items()
               if not k.startswith("_")}
    for z, v in reg_ref.items():                     # 事前登録値と ref CSV の突き合わせ
        if abs(meas["down"].get(round(z, 3), float("nan")) - v) > 1e-9:
            print(f"\n警告: acceptance.json の実測 z/W={z} ({v}) と ref CSV が一致しない")
    denom = tol_item(reg["tol"], "denominator")["interval"]
    dig = float(tol_item(reg["tol"], "digitize")["interval_abs"])

    print("\n=== 判定深さの h/h_fp (後壁, T_surf/T_gap = %.2f) ===" % a.t_ratio)
    print(f"伝播に入れた区間: denominator {denom} (上端のみ残す) / digitize 絶対 ±{dig}")
    print(f"{'z/W':>5} {'深さ[mm]':>9} {'forge':>9} {'実測':>7} {'比 下端':>9} "
          f"{'比 公称':>9} {'比 上端':>9}  備考")
    hi_at = {}
    show = sorted(set([0.7] + reg_depths))
    for z in show:
        fv = at_zw(walls["down"], z, W)
        mv = meas["down"].get(round(z, 3))
        note = "参考 (スキン滲みの拡散長内・acceptance の excluded)" if z == 0.7 else ""
        if mv is None:
            print(f"{z:5.1f} {z*W*1e3:9.3f} {fmt(fv)} {'-':>7} "
                  f"{'-':>9} {'-':>9} {'-':>9}  実測なし")
            continue
        lo, nom, hi = ratio_band(fv, mv, denom[1], dig)
        if z in reg_depths:
            hi_at[z] = hi
        print(f"{z:5.1f} {z*W*1e3:9.3f} {fmt(fv)} {mv:7.2f} {fmt(lo)} {fmt(nom)} "
              f"{fmt(hi)}  {note}")

    # ---- T2p-Q (定性) --------------------------------------------------
    accQ = [g for g in reg["acc"]["gates"] if g["id"] == "T2p-Q"][0]
    zs = np.linspace(0.1, QUAL_ZW_MAX, 40)
    dn = np.array([at_zw(walls["down"], z, W) for z in zs])
    up = np.array([at_zw(walls["up"], z, W) for z in zs])
    ok = np.isfinite(dn) & np.isfinite(up)
    print("\n=== T2p-Q (定性ゲート) ===")
    print("  基準 (acceptance.json): " + accQ["criterion"])
    if not ok.any():
        print("  VERDICT: 判定不能 (z/W 0.1..4 に壁データが無い)")
    else:
        bad = zs[ok][dn[ok] <= up[ok]]
        frac = 100.0 * (dn[ok] > up[ok]).sum() / ok.sum()
        print(f"  z/W 0.1..{QUAL_ZW_MAX:g} の {ok.sum()} 点中 {frac:.1f} % で 後壁 > 前壁")
        if bad.size:
            print(f"  不成立の z/W: {np.round(bad,2).tolist()}")
        print(f"  T2p-Q: {'PASS' if bad.size == 0 else 'FAIL'}")

    # ---- T2p-C (3 分類) ------------------------------------------------
    print("\n=== T2p-C (3 分類, 暫定) ===")
    for k, v in accC["classes"].items():
        print(f"  {k:6s}: {v}")
    print(f"  band_rule: {accC['band_rule']}")
    need = set(ALPHA_DEPTHS) | set(BETA_DEPTHS)
    if any(not np.isfinite(hi_at.get(z, float("nan"))) for z in need):
        verdict = "判定不能"
        why = "判定深さのどれかで上端が有限でない (壁データ範囲外 / 実測-0.03 <= 0)"
    elif all(hi_at[z] < ALPHA_MAX for z in ALPHA_DEPTHS):
        verdict, why = "alpha", f"z/W {ALPHA_DEPTHS} の上端がいずれも < {ALPHA_MAX:.3f}"
    elif all(BETA_BAND[0] <= hi_at[z] <= BETA_BAND[1] for z in BETA_DEPTHS):
        verdict, why = "beta", f"z/W {BETA_DEPTHS} の上端がいずれも [1/3, 3]"
    else:
        verdict, why = "gamma", "alpha / beta のいずれも成立しない"
    print(f"\n  暫定分類: {verdict}   ({why})")
    print("  上端: " + ", ".join(f"z/W {z} → {fmt(hi_at.get(z, float('nan')), 0, 4)}"
                                 for z in sorted(need)))

    # ---- 未実装の区間 (必ず出す) ---------------------------------------
    print("\n  **この分類は暫定である** — tolerances.json の区間のうち、"
          "実装したのは次の 2 つだけ:")
    print("    - denominator  [%.1f, %.1f]  (上端側のみ残す)" % (denom[0], denom[1]))
    print("    - digitize     絶対 ±%.2f" % dig)
    print("  **未実装 (この判定に入っていない) 区間**:")
    print("    - skin_smear : すきま薄肉 (304SS 0.406 mm) の横方向伝導による実測の嵩上げ。"
          "case/50/tools/skin_smear.py 相当の順方向モデルの移植が要る (別作業)。"
          "向きは**実測を下げる側 = 欠損が縮む側**なので、入れると上端はさらに上がる")
    print("    - inflow_BL  : 入口 BL の再構成 (dstar [0.7,1.4] / u_tau [0.7,1.4])。両側")
    print("    - grid       : 開口幅方向 2 段格子の差 (run が 2 本要る)")
    print("  → **未実装の区間を『無い』ものとして分類を確定させないこと**。"
          "case/49 へ渡す分類は、上の 3 つを入れてから決める (plan §5.1 #49)")
    print("  参考: T2p-M (数値要件) はこのツールの対象外。"
          "check_mesh_quality / check_wall_resolution / check_quasisteady / "
          "check_convergence の VERDICT を別に揃えること")

    # ---- 時系列 --------------------------------------------------------
    if a.series_csv:
        out = rd / "gap_series.csv" if a.series_csv == "AUTO" else Path(a.series_csv)
        write_series(dumps, W, D, reg, a, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
