"""ノズル δ* 総点検 (notes/investigations/nozzle-deltastar-throat-review.md) の再現スクリプト。

1. 同設計の Euler run (設計壁) と NS run (物理壁) の質量流量比 → 有効スロート半径・等価 δ*_t
2. 信頼域にある既存の ρUx 最大縁抽出と、与えた δ* を比較
3. 図 notes/investigations/figs/nozzle-deltastar-throat-review.png (PNG は .gitignore 対象、本スクリプトで再生成)

使い方 (リポジトリルートで):  design/.venv-opt/bin/python notes/investigations/scripts/nozzle_deltastar_throat_review.py
追加計算は不要 (既存 run の nozzle.h5 / res_*.h5 / wall_*.csv / prepare_info.json を読むだけ)。
"""
from __future__ import annotations
import json, sys, warnings
from pathlib import Path
import h5py, numpy as np

ROOT = Path(__file__).resolve().parents[3]          # forge リポジトリルート
CASE = ROOT / "case"
sys.path.insert(0, str(ROOT / "design"))
warnings.filterwarnings("ignore")

# (NS run, 同設計の Euler run) — 設計キー (scale_m, x_E, Md) が一致する組
PAIRS = [
    ("42.isobutane_wt/run_0068_ib_R2_LU6_Lc8_h2o5_split_ns_dry", "42.isobutane_wt/run_0064_ib_R2_LU6_Lc8_h2o5_split"),
    ("42.isobutane_wt/run_0100_ib_m5_R3_LU9_Lc14_ns", "42.isobutane_wt/run_0096_ib_m5_R3_LU9_Lc14_tp"),
    ("42.isobutane_wt/run_0102_ib_m5_R3_LU9_Lc14_ns_v3", "42.isobutane_wt/run_0096_ib_m5_R3_LU9_Lc14_tp"),
    ("42.isobutane_wt/run_0104_ib_m6_R2_Lc50_ns", "42.isobutane_wt/run_0055_ib_m6_R2_Lc50_MK2.5_tp"),
    ("42.isobutane_wt/run_0105_ib_m6_R2_Lc50_ns_v3", "42.isobutane_wt/run_0055_ib_m6_R2_Lc50_MK2.5_tp"),
    ("42.isobutane_wt/run_0107_ib_m6_R2_Lc50_ns_yp1", "42.isobutane_wt/run_0055_ib_m6_R2_Lc50_MK2.5_tp"),
    ("44.vitiated_air_wt/run_0017_va_R2_LU6_Lc8_ns_v1", "44.vitiated_air_wt/run_0005_va_R2_LU6_Lc8"),
    ("44.vitiated_air_wt/run_0021_va_R2_LU6_Lc8_ns_v3", "44.vitiated_air_wt/run_0005_va_R2_LU6_Lc8"),
    ("45.isobutane_m6_d155/run_0003_ns_v1", "45.isobutane_m6_d155/run_0001_euler_shortest_dry"),
    ("45.isobutane_m6_d155/run_0004_ns_v3", "45.isobutane_m6_d155/run_0001_euler_shortest_dry"),
    ("45.isobutane_m6_d155/run_0006_ns_trim_v1", "45.isobutane_m6_d155/run_0005_euler_trim"),
    ("45.isobutane_m6_d155/run_0007_ns_trim_v3", "45.isobutane_m6_d155/run_0005_euler_trim"),
    ("45.isobutane_m6_d155/run_0011_cpg_ns", "45.isobutane_m6_d155/run_0009_cpg_euler"),
    ("45.isobutane_m6_d155/run_0013_cpg_ns_plainsst", "45.isobutane_m6_d155/run_0009_cpg_euler"),
]


def load(rd: Path) -> dict:
    """構造格子 (node, index = i*nj + j) として場を読む。"""
    info = json.loads((rd / "prepare_info.json").read_text())
    ni, nj = info["mesh"]["ni"], info["mesh"]["nj"]
    S = float(info["scale_m"])
    with h5py.File(rd / "nozzle.h5") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))[-1]
    with h5py.File(res) as f:
        ro = f["/VALUE/ro"][:]; Ux = f["/VALUE/Ux"][:]
    if len(ro) != ni * nj or len(nc) != ni * nj:
        raise ValueError(f"{rd}: 構造格子 node run でない")
    d = dict(info=info, S=S, res=res.name)
    d["x"] = nc[:, 0].reshape(ni, nj) / S
    d["r"] = nc[:, 1].reshape(ni, nj) / S
    d["q"] = (ro * Ux).reshape(ni, nj)
    # 列ごとの質量流量 2π∫ρUx r dr [kg/s]
    d["mdot"] = 2 * np.pi * np.trapezoid(d["q"] * d["r"] * S, d["r"] * S, axis=1)
    return d


def mdot_median(d: dict) -> tuple[float, float]:
    xs = d["x"][:, 0]
    sel = (xs > -2.5) & (xs < d["info"]["x_E"])
    m = np.median(d["mdot"][sel])
    return float(m), float(d["mdot"][sel].std() / m)


def given_dstar(rd: Path, S: float, x):
    """この run が実際に与えた δ*(x) = 物理壁 − 設計壁 (半径差, r_t 単位)。"""
    wp = np.loadtxt(rd / "wall_physical.csv", delimiter=",", skiprows=1) / S
    wd = np.loadtxt(rd / "wall_design.csv", delimiter=",", skiprows=1); wd[:, :2] /= S
    x = np.asarray(x, dtype=float)
    g = np.interp(x, wp[:, 0], wp[:, 1]) - np.interp(x, wd[:, 0], wd[:, 1])
    return np.where(x >= wd[0, 0], g, np.nan), wp, wd


def sensitivity_n(Md: float, gam: float = 1.3) -> float:
    """n = d ln(A/A*) / d ln M (1D 等エントロピー)。"""
    return -1.0 + (gam + 1) * Md ** 2 / (2 + (gam - 1) * Md ** 2)


def table_massflow() -> None:
    print("### 1. 質量流量比 → 有効スロート")
    hdr = f"{'NS run':52s} {'Euler ref':38s} {'Md':>5s} {'mNS/mE':>8s} {'r_t,W':>7s} {'Cd_rel':>7s} {'δ*t given':>9s} {'δ*t eq':>8s} {'ΔM/M est':>8s}"
    print(hdr)
    for ns, eu in PAIRS:
        N = load(CASE / ns); E = load(CASE / eu)
        iN, iE = N["info"], E["info"]
        if not (abs(iN["scale_m"] - iE["scale_m"]) < 1e-9 and abs(iN["x_E"] - iE["x_E"]) < 1e-3 and iN["Md"] == iE["Md"]):
            print(f"{ns:52s} !! 設計キー不一致 ({eu})"); continue
        mN, _ = mdot_median(N); mE, _ = mdot_median(E)
        ratio = mN / mE
        rp = iN["throat_physical"]["r"]
        n = sensitivity_n(iN["Md"])
        print(f"{ns:52s} {eu.split('/')[1]:38s} {iN['Md']:5.2f} {ratio:8.5f} {rp:7.4f} {ratio / rp**2:7.4f} "
              f"{rp - 1:9.5f} {rp - np.sqrt(ratio):8.5f} {-(ratio - 1) / n * 100:+7.3f}%")


def figure(out: Path) -> None:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    try: fm.fontManager.addfont("/home/sano/.fonts/NotoSansCJKjp-Regular.otf")
    except Exception: pass
    matplotlib.rcParams["font.family"] = ["Noto Sans CJK JP", "DejaVu Sans"]
    from forge_design.metrics.deltastar import deltastar_from_run
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    runs = [("45.isobutane_m6_d155/run_0004_ns_v3", "case/45 M6 d1.55 (run_0004 v3)"),
            ("42.isobutane_wt/run_0105_ib_m6_R2_Lc50_ns_v3", "case/42 M6 R2/Lc50 (run_0105 v3)")]
    for ax, (run, title) in zip(axes, runs):
        rd = CASE / run; d = load(rd); S = d["S"]
        xg = np.linspace(-0.8, 30, 400)
        g, wp, _ = given_dstar(rd, S, xg)
        xr = np.arange(14.0, 30.01, 2.0)
        res = sorted(rd.glob("res_[0-9]*.h5"), key=lambda p: int(p.stem.split("_")[1]))[-1]
        dr = deltastar_from_run(rd / "nozzle.h5", res, wp, S, xr)
        ax.plot(xg, g, "-", lw=2, label="与えた δ* (物理壁 − 設計壁, v3)")
        ax.plot(xr, dr["dstar"], "s", ms=4, mfc="none", label="生産抽出 deltastar_from_run (ρu 最大縁, x≥x_lo=14)")
        ax.axhline(0.0015, color="gray", lw=0.8, ls="--")
        ax.text(12, 0.00165, "質量流量から得た等価スロート補正 ≈ 0.0015", fontsize=8, color="gray")
        ax.axvline(14, color="gray", lw=0.6, ls=":"); ax.text(14.3, 0.32, "x_lo", fontsize=8, color="gray")
        ax.set_yscale("log"); ax.set_ylim(5e-4, 0.5); ax.set_xlim(-1, 30)
        ax.set_xlabel("x / r_t"); ax.set_ylabel("δ* / r_t"); ax.set_title(title)
        ax.grid(alpha=.3, which="both"); ax.legend(fontsize=8, loc="upper left")
    fig.suptitle("与えた排除厚さ vs 信頼域の CFD 抽出 (x<x_lo の局所抽出値は未確定)", fontsize=11)
    fig.tight_layout(); out.parent.mkdir(parents=True, exist_ok=True); fig.savefig(out, dpi=120)
    print(f"\nsaved {out}")


if __name__ == "__main__":
    table_massflow()
    figure(ROOT / "notes/investigations/figs/nozzle-deltastar-throat-review.png")
