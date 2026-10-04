#!/usr/bin/env python3
"""二相拡散 A/B (plan condensation-two-phase-transport §5.1 #1b) の報告量を res_*.h5 ごとに 1 行ずつ CSV 化する。

    python3 twophase_ab_series.py RUN_DIR [--mask-from H5] [--gmin-mask 1e-6] [--out CSV]
      -> RUN_DIR/twophase_series.csv (check_quasisteady.py --series-csv 用; 列 `step` + 下の報告量)

報告量の定義 (A・B で同じ関数を通す。twophase_ab_compare.py も本モジュールの metrics() を使う):
  場所の取り方: case/16 のメッシュは x 一定の構造列 (列ごとに 120 節点, y∈[−12.7, 12.7] mm) を前提とする。
  - 中心線値: 各 x 列で y ≥ 0 の最下点と y < 0 の最上点を y について線形補間した y = 0 の値
    (このメッシュには y = 0 の節点が無い; compare_condfix.py の「|y| 最小の節点」は列ごとに上下が入れ替わるので使わない)。
  - 上壁: wall_dist ≤ 0 かつ y > 0 の節点のうち各 x 列で y 最大のもの (compare_condfix.py と同じ定義)。
  onset_c_g1e4_mm : 中心線 g (= g_0, 液質量分率) が x 上流から見て初めて 1e-4 を超える x [mm] (隣接列間の線形補間)。
  onset_c_g1e3_mm : 同じく 1e-3 (cond_series.csv の onset_mm と同じ閾値; 中心線の定義は上のとおり違う)。
  g_exit_mw       : 出口列 (x = x_max の全節点) の質量流束重み平均 Σ ρU_x g w / Σ ρU_x w。w は y 方向の双対長
                    (隣接節点との中点間の長さ; 端点は半分)。平面 2D (isAxisymmetric 0) なので周長重みは掛けない。
  g_exit_c        : 出口の中心線 g (上の補間)。
  pw_mean_x10     : 上壁 p/p0 (p0 = 入口全圧 59070 Pa) を Wyslouzil Fig.3 の実験点の x (x ≥ 10 mm の 12 点, wyslouzil_fig3_pp0.csv)
                    で補間した値の算術平均。
  pw21 / pw42 / pw52 : 上壁 p/p0 の x = 21 / 42 / 52 mm での補間値。
  dev_pct         : 上記 12 点で (p/p0 − 実験 1.00 kPa 凝縮値)/実験値 の平均 [%] (cond_series.csv の dev_pct と同じ定義・同じ壁節点)。
  Tw_mean_x10_K   : 上壁の静温を同じ 12 点で補間した平均 [K]。**壁熱流束の代わり**: 本ケースの壁は断熱
                    (bcondConfig `kind: wall`, 壁温指定なし) で、res_*.h5 に壁熱流束の出力は無い (qwall なし)。
                    断熱壁なので壁熱流束は境界条件で 0 であり、拡散作用素の差は壁温 (回復温度) に出る。
  T_cond0_vmean_K : 固定マスク M0 上の体積重み平均静温 [K]。M0 = 共通初期場 (run の valueFileName = restart で写した
                    run_0482 res_48000) で g = ρg/ρ > 1e-6 の節点。A・B の初期場は同一なので M0 は両者共通で時間不変
                    (各時刻の「両者とも g > 1e-6」の共通域は twophase_ab_compare.py が最終時刻で出す)。体積は VALUE/volume
                    (output.level 2 か extraFields: [volume] が要る)。
  g_cond0_vmean   : M0 上の体積重み平均 g。
  gmax            : 全域の g 最大。
  nonfinite       : 読んだ量 (P, T, ro, Ux, Uy, g_0) の非有限値の数 (判定列に入れない; 監視用)。
  検査列 (#1b-pre (5); 定常性の判定ではなく各スナップショットで 0 / 非負であることを直接見る):
  nonfinite_cons  : 全保存量 (VALUE の ro, roUx, roUy, roUz, roe, roK, roOmega, roY*, rog_*, roQ*_* のうち存在するもの) の非有限値の数。
  min_rv / min_rg / min_rQ2 / min_rQ1 / min_rQ0 : 蒸気 ρY_w − ρg (w = condGasSpecies)・液 ρg・各 ρQ の全節点最小 (float の格納値を double で差)。
  neg_rv / neg_rg / neg_rQ : 同じ量が負の節点数 (Q は 3 本の合計)。受入には 0 を要求する (plan §5.1 #1b)。
壁熱流束の列は作らない (上記のとおり出力が無く、断熱壁で 0)。

使い方の例 (判定は check_quasisteady.py。許容は事前登録で決める):
  python3 solver_density_cuda/tools/check_quasisteady.py --series-csv RUN/twophase_series.csv \\
      --series-cols onset_c_g1e4_mm,g_exit_mw,pw_mean_x10,Tw_mean_x10_K,T_cond0_vmean_K --drift 0.001 --osc 0.001
"""
import argparse, glob, os, sys
import h5py, numpy as np, yaml

HERE = os.path.dirname(os.path.abspath(__file__))
P0 = 59070.0                     # 入口全圧 [Pa] (run_0482 bcondConfig inlet Pt)
COLS = ["step", "onset_c_g1e4_mm", "onset_c_g1e3_mm", "g_exit_mw", "g_exit_c", "pw_mean_x10", "pw21", "pw42", "pw52",
        "dev_pct", "Tw_mean_x10_K", "T_cond0_vmean_K", "g_cond0_vmean", "gmax", "nonfinite",
        "nonfinite_cons", "min_rv", "min_rg", "min_rQ2", "min_rQ1", "min_rQ0", "neg_rv", "neg_rg", "neg_rQ"]
CONS_FIXED = ("ro", "roUx", "roUy", "roUz", "roe", "roK", "roOmega")

_exp = np.genfromtxt(os.path.join(HERE, "wyslouzil_fig3_pp0.csv"), delimiter=",", skip_header=1)[:, :3]
XE = _exp[:, 0]*10.0             # [mm]
COND_EXP = _exp[:, 2]
SEL10 = XE >= 10.0


def res_files(run):
    fs = glob.glob(os.path.join(run, "res_[0-9]*.h5"))
    out = []
    for f in fs:
        b = os.path.basename(f)[4:-3]
        if b.isdigit():
            out.append((int(b), f))
    return sorted(out)


def value_file(run):
    cfg = yaml.safe_load(open(os.path.join(run, "solverConfig.yaml")))
    m = cfg["mesh"]
    return os.path.join(run, m.get("valueFileName", m["meshFileName"]))


class Geometry:
    """列構造・中心線補間・上壁・出口の索引 (1 run で 1 回作る)。"""

    def __init__(self, coord, wall_dist):
        c = coord.reshape(-1, 3)
        self.xy = c[:, :2].astype(np.float64)
        xs = np.round(c[:, 0], 6)
        self.ux = np.unique(xs)
        self.x_mm = self.ux*1e3
        y = self.xy[:, 1]
        lo, hi, wlo, whi, wall = [], [], [], [], []
        is_wall = (wall_dist <= 0) & (y > 0)
        for xv in self.ux:
            cand = np.where(xs == xv)[0]
            up = cand[y[cand] >= 0]; dn = cand[y[cand] < 0]
            if not len(up) or not len(dn):
                sys.exit(f"x = {xv} の列に y≥0 と y<0 の両方の節点が無い (中心線を補間できない)")
            iu = up[np.argmin(y[up])]; idn = dn[np.argmax(y[dn])]
            yu, yd = y[iu], y[idn]
            t = (0.0 - yd)/(yu - yd)
            lo.append(idn); hi.append(iu); wlo.append(1.0 - t); whi.append(t)
            cw = cand[is_wall[cand]]
            wall.append(cw[np.argmax(y[cw])] if len(cw) else -1)
        self.c_lo, self.c_hi = np.array(lo), np.array(hi)
        self.c_wlo, self.c_whi = np.array(wlo), np.array(whi)
        wall = np.array(wall)
        self.wall = wall[wall >= 0]
        self.xw_mm = self.xy[self.wall, 0]*1e3
        # 出口列: x = x_max の全節点を y で並べ、双対長 w を付ける
        out = np.where(xs == self.ux[-1])[0]
        out = out[np.argsort(y[out])]
        yo = y[out]
        w = np.empty(len(out))
        w[1:-1] = 0.5*(yo[2:] - yo[:-2]); w[0] = 0.5*(yo[1] - yo[0]); w[-1] = 0.5*(yo[-1] - yo[-2])
        self.out, self.out_w = out, w

    def centerline(self, f):
        return self.c_wlo*f[self.c_lo] + self.c_whi*f[self.c_hi]


def onset(x_mm, g, thr):
    """上流から初めて g > thr となる x を隣接点の線形補間で返す (compare_condfix.py と同じ)。無ければ NaN。"""
    i = np.where(g > thr)[0]
    if not len(i):
        return np.nan
    k = i[0]
    if k == 0 or g[k] == g[k-1]:
        return x_mm[k]
    return x_mm[k-1] + (thr - g[k-1])/(g[k] - g[k-1])*(x_mm[k] - x_mm[k-1])


def load_mask(path, gmin):
    with h5py.File(path, "r") as f:
        V = f["VALUE"]
        if "rog_0" not in V or "ro" not in V:
            sys.exit(f"{path}: VALUE/rog_0 か VALUE/ro が無い (共通初期場の凝縮域 M0 を作れない)。"
                     " restart で液を写していない値ファイルの可能性 (prepare_twophase_ab.sh を確認)。--mask-from で指定もできる")
        g0 = np.asarray(V["rog_0"], np.float64)/np.maximum(np.asarray(V["ro"], np.float64), 1e-30)
    return g0 > gmin


def is_conserved(k):
    return k in CONS_FIXED or (k.startswith("roY") and k[3:].isdigit()) or k.startswith("rog_") or \
        (k.startswith("roQ") and "_" in k and k[3:k.index("_")].isdigit())


def read(fn, iw=1):
    with h5py.File(fn, "r") as f:
        V = f["VALUE"]
        d = {k: np.asarray(V[k], np.float64) for k in ("P", "T", "ro", "Ux", "Uy", "g_0", "wall_dist", "volume") if k in V}
        d["coord"] = np.asarray(f["MESH/COORD"])
        cons = [k for k in V if is_conserved(k)]
        d["nonfinite_cons"] = sum(int(np.count_nonzero(~np.isfinite(np.asarray(V[k])))) for k in cons)
        d["n_cons"] = len(cons)
        need = [f"roY{iw}", "rog_0", "roQ2_0", "roQ1_0", "roQ0_0"]
        miss = [k for k in need if k not in V]
        if miss:
            sys.exit(f"{fn}: VALUE/{miss} が無い (非負の検査列を作れない; output.level 2 で出る)")
        # float の格納値をそのまま double に上げて差を取る (蒸気 = 総水分 − 液)
        rw, rg = np.asarray(V[f"roY{iw}"], np.float64), np.asarray(V["rog_0"], np.float64)
        d["rv"] = rw - rg; d["rg"] = rg
        d["rQ"] = [np.asarray(V[k], np.float64) for k in ("roQ2_0", "roQ1_0", "roQ0_0")]
    for k in ("P", "T", "ro", "Ux", "Uy", "g_0", "volume"):
        if k not in d:
            sys.exit(f"{fn}: VALUE/{k} が無い (output.level 2 か extraFields に {k} を入れる)")
    return d


def metrics(d, geo, mask0):
    """1 スナップショットの報告量 (dict)。定義はモジュール docstring。"""
    g, P, T = d["g_0"], d["P"], d["T"]
    gc = geo.centerline(g)
    pw = P[geo.wall]/P0
    Tw = T[geo.wall]
    pi = np.interp(XE, geo.xw_mm, pw)
    Twi = np.interp(XE, geo.xw_mm, Tw)
    o, w = geo.out, geo.out_w
    mf = d["ro"][o]*d["Ux"][o]*w
    vol = d["volume"]
    m0v = vol[mask0]
    nonfin = sum(int(np.count_nonzero(~np.isfinite(d[k]))) for k in ("P", "T", "ro", "Ux", "Uy", "g_0"))
    return dict(
        onset_c_g1e4_mm=onset(geo.x_mm, gc, 1e-4),
        onset_c_g1e3_mm=onset(geo.x_mm, gc, 1e-3),
        g_exit_mw=float(np.sum(mf*g[o])/np.sum(mf)) if np.sum(mf) != 0 else np.nan,
        g_exit_c=float(gc[-1]),
        pw_mean_x10=float(pi[SEL10].mean()),
        pw21=float(np.interp(21.0, geo.xw_mm, pw)), pw42=float(np.interp(42.0, geo.xw_mm, pw)),
        pw52=float(np.interp(52.0, geo.xw_mm, pw)),
        dev_pct=float(((pi - COND_EXP)/COND_EXP*100.0)[SEL10].mean()),
        Tw_mean_x10_K=float(Twi[SEL10].mean()),
        T_cond0_vmean_K=float(np.sum(T[mask0]*m0v)/np.sum(m0v)) if m0v.sum() > 0 else np.nan,
        g_cond0_vmean=float(np.sum(g[mask0]*m0v)/np.sum(m0v)) if m0v.sum() > 0 else np.nan,
        gmax=float(np.max(g)), nonfinite=nonfin,
        nonfinite_cons=d.get("nonfinite_cons", np.nan),
        min_rv=float(np.min(d["rv"])), min_rg=float(np.min(d["rg"])),
        min_rQ2=float(np.min(d["rQ"][0])), min_rQ1=float(np.min(d["rQ"][1])), min_rQ0=float(np.min(d["rQ"][2])),
        neg_rv=int(np.count_nonzero(d["rv"] < 0)), neg_rg=int(np.count_nonzero(d["rg"] < 0)),
        neg_rQ=int(sum(np.count_nonzero(q < 0) for q in d["rQ"])))


def cond_gas_species(run):
    """run の solverConfig の condensation.condGasSpecies (総水分の化学種 index)。"""
    cfg = yaml.safe_load(open(os.path.join(run, "solverConfig.yaml")))
    iw = (cfg.get("condensation") or {}).get("condGasSpecies", 1)
    if not isinstance(iw, int):
        sys.exit(f"{run}: condensation.condGasSpecies が整数でない ({iw!r}; 名前指定は未対応)")
    return iw


def setup(run, mask_from=None, gmin_mask=1e-6, first=None):
    """run の Geometry と M0 を作る。first は最初に読むスナップショット (wall_dist/COORD の取得用)。"""
    fs = res_files(run)
    if not fs:
        sys.exit(f"{run}: res_*.h5 が無い")
    d0 = read(first or fs[0][1], cond_gas_species(run))
    if "wall_dist" not in d0:
        with h5py.File(value_file(run), "r") as f:
            d0["wall_dist"] = np.asarray(f["VALUE/wall_dist"], np.float64)
    geo = Geometry(d0["coord"], d0["wall_dist"])
    mask0 = load_mask(mask_from or value_file(run), gmin_mask)
    if len(mask0) != len(d0["T"]):
        sys.exit(f"M0 の節点数 {len(mask0)} != 場の節点数 {len(d0['T'])}")
    return fs, geo, mask0, d0["coord"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    ap.add_argument("--mask-from", default=None, help="M0 を作る h5 (既定: run の valueFileName = 共通初期場)")
    ap.add_argument("--gmin-mask", type=float, default=1e-6, help="M0 の g 閾値 (既定 1e-6)")
    ap.add_argument("--out", default=None, help="出力 CSV (既定 RUN/twophase_series.csv)")
    a = ap.parse_args()
    fs, geo, mask0, coord0 = setup(a.run, a.mask_from, a.gmin_mask)
    rows = []
    iw = cond_gas_species(a.run)
    for st, fn in fs:
        d = read(fn, iw)
        if not np.array_equal(d["coord"], coord0):
            sys.exit(f"{fn}: MESH/COORD が他のスナップショットと違う")
        m = metrics(d, geo, mask0)
        rows.append([st] + [m[c] for c in COLS[1:]])
    out = a.out or os.path.join(a.run, "twophase_series.csv")
    np.savetxt(out, np.array(rows, dtype=float), delimiter=",", header=",".join(COLS), comments="", fmt="%.9g")
    print(f"{out}  ({len(rows)} rows, steps {fs[0][0]}..{fs[-1][0]}; M0 = {int(mask0.sum())} nodes)")
    arr = np.array(rows, dtype=float)
    ci = {c: i for i, c in enumerate(COLS)}
    bad = int(np.sum(arr[:, ci["nonfinite_cons"]] != 0) + np.sum(arr[:, ci["neg_rv"]] != 0)
              + np.sum(arr[:, ci["neg_rg"]] != 0) + np.sum(arr[:, ci["neg_rQ"]] != 0))
    print(f"検査列 (全スナップショット): 保存量の非有限 最大 {int(arr[:, ci['nonfinite_cons']].max())}・"
          f"負の節点 蒸気 {int(arr[:, ci['neg_rv']].max())} 液 {int(arr[:, ci['neg_rg']].max())} Q {int(arr[:, ci['neg_rQ']].max())} (最大値)、"
          f"min ρv {arr[:, ci['min_rv']].min():.3e} ρg {arr[:, ci['min_rg']].min():.3e} → {'OK' if bad == 0 else 'NG (' + str(bad) + ' スナップショット·項目)'}")


if __name__ == "__main__":
    main()
