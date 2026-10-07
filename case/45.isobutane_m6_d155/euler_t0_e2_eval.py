"""plan verification-case45-euler-total-enthalpy §6 E2 (2026-10-07 登録、改訂 a2259768・f158bac2・再改訂 a59e61b5。結果を見る前に実装) の評価器。
半径方向の配点プロファイルへの感度: A = run_0161_euler_t0cluster_g1 (G1 のまま) / B = run_0162_euler_t0cluster_u5em3 (両 wall_first_frac 0.005)。
起動は run_euler_t0_e2.sh (準備・実行は euler_t0_e2.py)。この評価器は保存場を読むだけで、forge は起動しない。

全温の復元: E1 と同じ (euler_t0_stage_ab.recon = V5c の経路 B の Euler 版)。保存量 ro・roU*・roe・roY* と、その run の解決済み熱物性から
  float64 で T を Newton で解き、h0 = (roe + P)/ro、h_mix(T0) = h0 で全温 T0 (保存された P・T・h0 の値は使わない)。偏差は T0 − 1600 K。

長さ (再改訂 a59e61b5): soft 3000 + 本段 54000 step (出力 1000 ごと)。本段の step 数・出力間隔・窓は下の MAIN_NSTEPS・OUT_INTERVAL から
  作り、実行側 (euler_t0_e2.py の prep・run) も同じ定数を読む。評価では各 run の本段の solverConfig.yaml・stage_manifest の実効値と、本段の残差の
  最終 step (= MAIN_NSTEPS − 1) を照合し、食い違えば固定の条件の不成立 (判別不能) にする。
見る時点: soft 段の後 (`_soft_stage/res_3000.h5`、euler_t0_e2.py run が runner の段の引き継ぎの前に退避したもの) と、本段の res_0〜res_54000
  (1000 ごと、55 枚)。判定に使うのは本段 42000〜54000 の 13 枚 (評価窓) と 50000〜54000 の 5 枚。

領域 (登録): 軸方向 x/r_t < −4 (up)・−4 ≤ x/r_t ≤ 1 (thr)・x/r_t > 1 (dn) と、半径方向 η = r/r_w < 0.1 (axis)・0.1 ≤ η ≤ 0.9 (core)・
  η > 0.9 (wall) の直積 9 領域 (鍵は "thr_wall" のように書く)。

主指標 (登録: 両側で共通の固定の評価点に同じ補間で移した T0 − 1600 の、領域ごとの最大・最小・|·| の 99 % 点・±1 K と ±100 K を超える割合):
  評価点 (実装時の読み、結果を見る前): 列 i (A と B で x が同一であることを確かめる) ごとに、A の節点の η と B の節点の η の和集合
    (np.unique)。補間は列内の η の区分線形 (np.interp)。両側とも同じ評価点・同じ補間。自分の節点では補間は節点の値そのものなので、
    評価点の最大・最小は各格子の節点の最大・最小に一致し、99 % 点と割合は両側で同じ標本密度で比べられる。
  99 % 点は numpy.quantile の既定 (linear)。割合は評価点の数に対する割合で、+1 K を超える (> +1)・−1 K を下回る (< −1)・|·| > 1 を別に出す
    (100 K も同じ)。
補助: 生の節点の最大・最小とその位置 (i, j、壁からの層、x/r_t、η。正負を分ける)・±1 K/±100 K を超える節点の数、
  h0 − h_mix(1600 K, Y) (J/kg と、cp_mix(1600 K, Y) で割った K 換算)。

判定の前提 (登録 f158bac2、P2 は再改訂 a59e61b5):
  P1 本段だけの区間の check_convergence (投入スクリプトが `--segment` で書いた CONVERGENCE_VERDICT_segment.txt。区間が `main` だけであること) が
     PASS、または全残差列の内訳を見て不合格の理由が停滞 (STALLED/plateau) だけ。RISING・DIVERGED・判定不能 (入力不備)・still converging を
     含む列が 1 つでもあれば不成立。総合表示 (`check_convergence.py:413–414` は停滞の列が 1 つあれば総合を plateau にする) は使わず、列ごとに読む。
     停滞だけの場合は「停滞のみ (NOT CONVERGED)」と記録し、PASS とは書かない。
  P2 主指標は偏差のまま (9 領域 × 最大・最小・|偏差| の 99 % 点 = 27 列)。両腕・全 27 列について、本段 42000〜54000 の 13 枚と
     50000〜54000 の 5 枚のそれぞれで、窓の時点がそろい、全値が有限、かつ時間方向の最大 − 最小 ≤ 0.1 K (窓の CSV に書いた値で測る)。
     原系列の `check_quasisteady --series-csv --tail 1` の VERDICT は併記するが、STEADY は必須にしない (平均で正規化するので、近零の偏差では
     絶対の幅の条件より過度に厳しくなる)。これは指定の窓の中での絶対の変動幅による感度の診断で、定常化・収束の認定ではない。
     T0 そのものを渡して STEADY を得る方法や、近零の列だけ判定を切り替える方法は採らない (登録)。比較は float の ≤ そのままで許容を足さない
     (実装時の読み: 300.1 − 300 = 0.1 + 2.3e-14 は不成立)。
  固定の条件 (登録「固定」の確認; 不成立なら比較が登録のものでないので判別不能): 両 run の forge の sha256 (本段と soft 段の RUN_PROVENANCE) が
     同じで run_0143 の RUN_PROVENANCE と同じ / stage_manifest の段が soft → main / 起動前の IC の検査 (E2_PREP.json) が両側 OK /
     入口の Tt = 1600 K / A と B の列の x と壁の r が同一 / 本段の step 数・出力間隔の実効値 (solverConfig.yaml・stage_manifest) と本段の残差の
     最終 step が MAIN_NSTEPS・OUT_INTERVAL と一致 / 評価窓の 13 枚がそろい、復元の Newton が全節点で収束し、値が有限。
判定 (登録の 3 区分。評価窓の全保存時点 = 13 枚):
  支持: 前提を満たし、A のスロートの区間 (thr の 3 領域を合わせた主指標) で |T0 − 1600| の最大 > 100 K が 13 枚すべてで再現し、
        B が全 9 領域で最大 ≤ +1 K かつ最小 ≥ −1 K を 13 枚すべてで保つ
        → 配点の介入への依存を支持し、「総点数だけで異常が決まる」説を退ける (この観測区間について)。
  退ける: 前提を満たし、両側でスロートの区間の > 100 K の異常が 13 枚すべてで残り、全 9 領域の主指標 (最大・最小・99 % 点) の A − B が
        13 枚すべてで |·| ≤ 1 K → 「今回の配点の変更がこの観測区間の異常を支配的に変える」を退ける (収束解に対する配点の依存までは否定しない)。
  判別不能: それ以外 (前提の未達・中間的な改善)。スロートだけが改善しても下流を含めた単一の原因とは認定しない。窓の変更や自動の延長で救済しない。
  「異常」は登録の「> 100 K の異常」を、スロートの 3 領域を合わせた主指標の max(|最大|, |最小|) > 100 K と読む (各領域すべてで超える要求では
  ない; a59e61b5 の読みの補足)。割合は評価点の数の割合。割合は K の量でないので「差 1 K」の対象に入れない
  (登録の「領域ごとの最大・最小・99 % 点」)。plateau を PASS と読み替えない。差が小さいことを「収束解が一致」と書かない。
出力: `_band_ab/euler_t0_e2_eval.json` (評価器の sha256・登録の commit・前提の定義と実効の判定・各時点の主指標と補助)、各 run の
  `euler_t0_e2_series.csv` (全時点)・`euler_t0_e2_series_{win13,tail5}.csv` (窓)・`euler_t0_e2_qs_{win13,tail5}.txt` (check_quasisteady の出力)。
usage: python3 euler_t0_e2_eval.py [case_dir] [--out PATH]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOOLS = ROOT / "solver_density_cuda/tools"

PLAN = "plans/active/verification-case45-euler-total-enthalpy.md §6 E2"
PLAN_REG_COMMIT = "a59e61b5"                     # 再改訂 (P2 の絶対の幅・本段 54000・窓 42000〜54000)。初稿 3073451a、改訂 a2259768、前提 f158bac2
PLAN_REG_HISTORY = ("3073451a", "a2259768", "f158bac2", "a59e61b5")
RUNS = {"A": "run_0161_euler_t0cluster_g1", "B": "run_0162_euler_t0cluster_u5em3"}
REF_RUN_BIN = "run_0143_euler_wallfit_monoG1_r1"  # forge の sha256 の基準 (登録「solver のバイナリ」を固定)
TT_REG = 1600.0                                  # [K] 登録の全温の基準
SOFT_DIR = "_soft_stage"                         # soft 段の出力の退避先 (euler_t0_e2.py run が書く; runner の glob("res_*") に当たらない)
SOFT_STEPS = 3000
MAIN_NSTEPS = 54000                              # 本段の step 数 (再改訂 a59e61b5)。実行側 euler_t0_e2.py もこの値を使う
OUT_INTERVAL = 1000                              # 本段の出力間隔
MAIN_STEPS = tuple(range(0, MAIN_NSTEPS + 1, OUT_INTERVAL))                        # 本段の出力 (55 枚)
WIN13 = tuple(range(MAIN_NSTEPS - 12 * OUT_INTERVAL, MAIN_NSTEPS + 1, OUT_INTERVAL))  # 評価窓 42000〜54000 (13 枚)
TAIL5 = tuple(range(MAIN_NSTEPS - 4 * OUT_INTERVAL, MAIN_NSTEPS + 1, OUT_INTERVAL))   # 末尾 50000〜54000 (5 枚)
WINDOWS = {"win13": WIN13, "tail5": TAIL5}
QS_TAIL = 1.0                                    # check_quasisteady --tail 1 (窓全体を評価)
RANGE_TOL_K = 0.1                                # [K] 主指標の時間方向の最大 − 最小の上限 (登録 f158bac2; a59e61b5 で P2 の唯一の幅の条件)
ANOMALY_K = 100.0                                # [K] 「> 100 K の異常」
WITHIN_K = 1.0                                   # [K] 「±1 K 以内」
DIFF_K = 1.0                                     # [K] 「差が 1 K 以内」
EXCEED_K = (1.0, 100.0)                          # 割合・節点数を数えるしきい値
X_EDGES_RT = (-4.0, 1.0)                         # up: x < −4 / thr: −4 ≤ x ≤ 1 / dn: x > 1
ETA_EDGES = (0.1, 0.9)                           # axis: η < 0.1 / core: 0.1 ≤ η ≤ 0.9 / wall: η > 0.9
X_BANDS = ("up", "thr", "dn")
ETA_BANDS = ("axis", "core", "wall")
REGION_KEYS = tuple(f"{xb}_{eb}" for xb in X_BANDS for eb in ETA_BANDS)
THROAT_REGIONS = tuple(f"thr_{eb}" for eb in ETA_BANDS)
QS_STATS = ("max", "min", "q99")                 # 前提 P2 と「差 1 K」の対象 (登録「領域ごとの最大・最小・99 % 点」)
FRAC_STATS = ("fpos1", "fneg1", "fabs1", "fpos100", "fneg100", "fabs100")
SEGMENT_VERDICT_FILE = "CONVERGENCE_VERDICT_segment.txt"   # 投入スクリプトが check_convergence --segment の出力を書く (V5 と同じ名前)
SEGMENT_MAIN = "main"                            # 本段だけの区間 (stage_manifest の最後の区間)
REQUIRED_RESID_COLS = ("rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe")   # check_convergence.REQUIRED_COLS と同じ
CONV_DEFINITION = ("本段だけの区間 (stage_manifest の最後の区間 = main) の check_convergence が PASS、または全残差列の内訳で不合格の理由が"
                   "停滞 (STALLED/plateau) だけ (RISING・DIVERGED・判定不能・still converging の列があれば不成立。総合表示は使わない)")
QS_DEFINITION = (f"主指標 (偏差のまま、9 領域 × {'・'.join(QS_STATS)}) の時系列について、両腕・全 27 列で、本段 {WIN13[0]}〜{WIN13[-1]} の "
                 f"13 枚と {TAIL5[0]}〜{TAIL5[-1]} の 5 枚のそれぞれの時点がそろい、全値が有限、かつ時間方向の最大 − 最小 ≤ {RANGE_TOL_K} K。"
                 f"check_quasisteady --series-csv --tail {QS_TAIL:g} の VERDICT は併記するが STEADY は必須にしない")
P2_NOTE = "指定の窓の中での絶対の変動幅による感度の診断で、定常化・収束の認定ではない"
SERIES_CSV = "euler_t0_e2_series.csv"
OUT_JSON = "_band_ab/euler_t0_e2_eval.json"
STAGES_EXPECTED = ["soft", "main"]

LBL_SUPPORT = "支持: 配点の介入への依存を支持し、「総点数だけで異常が決まる」説を退ける (この観測区間について)"
LBL_REJECT = "退ける: 「今回の配点の変更がこの観測区間の異常を支配的に変える」を退ける (収束解に対する配点の依存までは否定しない)"
LBL_UNDET = "判別不能"
SCOPE_NOTE = (f"結論は登録した計算手順・固定の予算 (soft {SOFT_STEPS} + 本段 {MAIN_NSTEPS})・観測区間 (本段 {WIN13[0]}〜{WIN13[-1]}) についてのもの。"
              "P2 は指定の窓の中での絶対の変動幅による感度の診断で、定常化・収束の認定ではない。停滞のみの "
              "NOT CONVERGED は PASS ではなく、収束の証明でもない。差が小さいことを「収束解が一致」と読まない。改善を float32 の幾何の誤差の証明としない")


def _sha(p) -> str | None:
    p = Path(p)
    if not p.is_file():
        return None
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


# --- 領域と評価点 ------------------------------------------------------------------------------------------------------------
def x_band(x_rt) -> np.ndarray:
    """0 = up (x < −4)、1 = thr (−4 ≤ x ≤ 1)、2 = dn (x > 1)。"""
    x = np.asarray(x_rt, dtype=np.float64)
    return np.where(x < X_EDGES_RT[0], 0, np.where(x <= X_EDGES_RT[1], 1, 2))


def eta_band(eta) -> np.ndarray:
    """0 = axis (η < 0.1)、1 = core (0.1 ≤ η ≤ 0.9)、2 = wall (η > 0.9)。"""
    e = np.asarray(eta, dtype=np.float64)
    return np.where(e < ETA_EDGES[0], 0, np.where(e <= ETA_EDGES[1], 1, 2))


def region_id(x_rt, eta) -> np.ndarray:
    """REGION_KEYS の番号 (3·x の帯 + η の帯)。"""
    return 3 * x_band(x_rt) + eta_band(eta)


def column_eta(R) -> np.ndarray:
    """(ni, nj) の r から η = r/r_w (j = nj−1 が壁)。各列で 0 から 1 へ狭義単調増加であることを確かめる。"""
    R = np.asarray(R, dtype=np.float64)
    if np.any(R[:, -1] <= 0):
        raise ValueError("壁の r が正でない列がある")
    eta = R / R[:, -1:]
    if not (np.all(eta[:, 0] == 0.0) and np.all(np.diff(eta, axis=1) > 0)):
        raise ValueError("列の η が 0 から始まる狭義単調増加でない (j = 0 が軸、j = nj−1 が壁の構造格子でない)")
    return eta


def eval_points(XA, RA, XB, RB) -> dict:
    """共通の固定の評価点。X・R は (ni, nj) [r_t 単位]。列の x と壁の r が A と B で同一であることを要求し (違えば ValueError)、
    列ごとに A と B の節点の η の和集合を評価点にする。戻り: {col, eta, x_rt, region, start (列 i の評価点は start[i]:start[i+1])}。"""
    XA, XB = np.asarray(XA, dtype=np.float64), np.asarray(XB, dtype=np.float64)
    if XA.shape != XB.shape or np.asarray(RA).shape != XA.shape or np.asarray(RB).shape != XB.shape:
        raise ValueError(f"A と B の格子の形が違う ({XA.shape} / {XB.shape})")
    if not (np.array_equal(XA, XB) and np.all(XA == XA[:, :1])):
        raise ValueError("列の x が A と B で同一でない (または列の中で x が一定でない)")
    if not np.array_equal(np.asarray(RA)[:, -1], np.asarray(RB)[:, -1]):
        raise ValueError("壁の r が A と B で同一でない")
    eA, eB = column_eta(RA), column_eta(RB)
    etas, cols, start = [], [], [0]
    for i in range(XA.shape[0]):
        u = np.unique(np.concatenate([eA[i], eB[i]]))
        etas.append(u)
        cols.append(np.full(u.size, i, dtype=np.int64))
        start.append(start[-1] + u.size)
    eta = np.concatenate(etas)
    col = np.concatenate(cols)
    x_rt = XA[col, 0]
    return {"col": col, "eta": eta, "x_rt": x_rt, "region": region_id(x_rt, eta), "start": np.asarray(start, dtype=np.int64)}


def eval_points_digest(ep: dict) -> str:
    h = hashlib.sha256()
    for k in ("col", "eta"):
        a = np.ascontiguousarray(ep[k])
        h.update(k.encode() + str(a.dtype).encode() + a.tobytes())
    return h.hexdigest()


def to_points(V, eta_mesh, ep: dict) -> np.ndarray:
    """節点の値 V (ni, nj) を評価点へ (列内の η の区分線形; 両側で同じ関数)。"""
    V = np.asarray(V, dtype=np.float64)
    out = np.empty(ep["eta"].size)
    st = ep["start"]
    for i in range(V.shape[0]):
        s, e = st[i], st[i + 1]
        out[s:e] = np.interp(ep["eta"][s:e], eta_mesh[i], V[i])
    return out


def indicators(vals, rid) -> dict:
    """主指標: 領域ごとの最大・最小・|·| の 99 % 点・割合 (評価点の数に対して)。空の領域は n = 0 と None。"""
    vals = np.asarray(vals, dtype=np.float64)
    out = {}
    for k, key in enumerate(REGION_KEYS):
        v = vals[rid == k]
        if v.size == 0:
            out[key] = {"n": 0, **{s: None for s in QS_STATS + FRAC_STATS}}
            continue
        a = np.abs(v)
        rec = {"n": int(v.size), "max": float(v.max()), "min": float(v.min()), "q99": float(np.quantile(a, 0.99))}
        for t in EXCEED_K:
            tg = f"{t:g}"
            rec[f"fpos{tg}"] = float(np.mean(v > t))
            rec[f"fneg{tg}"] = float(np.mean(v < -t))
            rec[f"fabs{tg}"] = float(np.mean(a > t))
        out[key] = rec
    return out


def _loc(i, j, X, R, eta) -> dict:
    nj = X.shape[1]
    return {"i": int(i), "j": int(j), "j_from_wall": int(nj - 1 - j), "x_rt": float(X[i, j]), "r_rt": float(R[i, j]),
            "eta": float(eta[i, j])}


def raw_stats(dev, X, R, eta) -> dict:
    """補助: 生の節点の統計 (領域は各節点の x/r_t と自分の格子の η)。最大・最小は位置つき (正負を分ける)。"""
    dev = np.asarray(dev, dtype=np.float64)
    rid = region_id(X, eta)
    out = {}
    for k, key in list(enumerate(REGION_KEYS)) + [(-1, "all")]:
        m = np.ones(dev.shape, dtype=bool) if k < 0 else (rid == k)
        n = int(m.sum())
        if n == 0:
            out[key] = {"n_nodes": 0}
            continue
        d = np.where(m, dev, np.nan)
        imax = np.unravel_index(int(np.nanargmax(d)), d.shape)
        imin = np.unravel_index(int(np.nanargmin(d)), d.shape)
        v = dev[m]
        rec = {"n_nodes": n, "max": float(v.max()), "argmax": _loc(*imax, X, R, eta),
               "min": float(v.min()), "argmin": _loc(*imin, X, R, eta)}
        for t in EXCEED_K:
            tg = f"{t:g}"
            rec[f"n_above_plus{tg}K"] = int(np.sum(v > t))
            rec[f"n_below_minus{tg}K"] = int(np.sum(v < -t))
        out[key] = rec
    return out


def h0_stats(dh, cp, X, eta) -> dict:
    """補助: h0 − h_mix(1600 K, Y) の領域ごとの最大・最小・|·| の 99 % 点 (J/kg と K 換算)。"""
    dh = np.asarray(dh, dtype=np.float64)
    dk = dh / np.asarray(cp, dtype=np.float64)
    rid = region_id(X, eta)
    out = {}
    for k, key in list(enumerate(REGION_KEYS)) + [(-1, "all")]:
        m = np.ones(dh.shape, dtype=bool) if k < 0 else (rid == k)
        if not m.any():
            out[key] = None
            continue
        a, b = dh[m], dk[m]
        out[key] = {"max_Jkg": float(a.max()), "min_Jkg": float(a.min()), "q99_abs_Jkg": float(np.quantile(np.abs(a), 0.99)),
                    "max_K": float(b.max()), "min_K": float(b.min()), "q99_abs_K": float(np.quantile(np.abs(b), 0.99))}
    return out


# --- 本段の step 数の同期 (実行側 euler_t0_e2.py と共有) --------------------------------------------------------------------
def window_problems() -> list:
    """評価器の窓の定数が MAIN_NSTEPS・OUT_INTERVAL から作られた形か (片方だけ書き換えた食い違いを拾う)。"""
    bad = []
    if MAIN_STEPS != tuple(range(0, MAIN_NSTEPS + 1, OUT_INTERVAL)):
        bad.append(f"MAIN_STEPS が 0〜{MAIN_NSTEPS} の {OUT_INTERVAL} ごとでない")
    if len(WIN13) != 13 or WIN13[-1] != MAIN_NSTEPS or any(b - a != OUT_INTERVAL for a, b in zip(WIN13, WIN13[1:])):
        bad.append(f"評価窓 {WIN13[0]}〜{WIN13[-1]} ({len(WIN13)} 枚) が本段の最後の 13 枚でない (本段 {MAIN_NSTEPS})")
    if len(TAIL5) != 5 or TAIL5 != WIN13[-5:]:
        bad.append(f"末尾の窓 {TAIL5} が評価窓の最後の 5 枚でない")
    return bad


def config_steps(cfg_text: str) -> dict:
    """solverConfig の本段の step 数と出力間隔の実効値 (yaml_strict: 重複キーは例外; solver の先勝ちと PyYAML の後勝ちの食い違いを避ける)。"""
    if str(TOOLS) not in sys.path:
        sys.path.insert(0, str(TOOLS))
    import yaml_strict
    doc = yaml_strict.load(cfg_text)
    t = (doc or {}).get("time") or {}
    n = (t.get("last") or {}).get("nStepOuter")
    o = t.get("outStepInterval")
    return {"nStepOuter": (int(n) if isinstance(n, (int, float)) and not isinstance(n, bool) else n),
            "outStepInterval": (int(o) if isinstance(o, (int, float)) and not isinstance(o, bool) else o)}


def step_problems(eff: dict) -> list:
    """実効値 {nStepOuter, outStepInterval[, manifest_main_nStepOuter, residual_last_step]} と評価器の定数の食い違い。"""
    bad = list(window_problems())
    if eff.get("nStepOuter") != MAIN_NSTEPS:
        bad.append(f"本段の nStepOuter {eff.get('nStepOuter')!r} が評価器の MAIN_NSTEPS {MAIN_NSTEPS} と違う")
    if eff.get("outStepInterval") != OUT_INTERVAL:
        bad.append(f"本段の outStepInterval {eff.get('outStepInterval')!r} が評価器の OUT_INTERVAL {OUT_INTERVAL} と違う")
    if "manifest_main_nStepOuter" in eff and eff["manifest_main_nStepOuter"] != MAIN_NSTEPS:
        bad.append(f"stage_manifest の main の nStepOuter {eff['manifest_main_nStepOuter']!r} が {MAIN_NSTEPS} と違う")
    if "residual_last_step" in eff and eff["residual_last_step"] != MAIN_NSTEPS - 1:
        bad.append(f"本段の残差の最終 step {eff['residual_last_step']!r} が {MAIN_NSTEPS - 1} でない (本段が途中で終わった・別の長さ)")
    return bad


def run_steps(rd: Path, conv: dict | None = None) -> dict:
    """run の本段の step 数の実効値: solverConfig.yaml (run_staged が最後に書いた本段の config)、stage_manifest の main、本段の残差の最終 step。"""
    eff = {}
    try:
        eff.update(config_steps((rd / "solverConfig.yaml").read_text()))
    except Exception as e:  # noqa: BLE001
        eff.update(nStepOuter=f"読めない: {type(e).__name__}: {e}", outStepInterval=None)
    try:
        man = json.loads((rd / "stage_manifest.json").read_text())
        mains = [st for st in man.get("stages", []) if st.get("tag") == SEGMENT_MAIN]
        v = str((mains[-1].get("soft") or {}).get("nStepOuter", "")) if mains else ""
        m = re.match(r"\d+", v)
        eff["manifest_main_nStepOuter"] = int(m.group(0)) if m else f"読めない ({v!r})"
    except (OSError, ValueError, AttributeError, TypeError) as e:
        eff["manifest_main_nStepOuter"] = f"読めない: {type(e).__name__}"
    if conv is not None:
        eff["residual_last_step"] = conv.get("last_step")
    return eff


# --- 前提 P1: check_convergence の列ごとの内訳 ---------------------------------------------------------------------------------
def parse_convergence(text) -> dict:
    """check_convergence.py <run> --segment の出力 → {status, head, segment, columns {列: {msg, class}}, reasons}。
    status: pass (全列合格) / stall_only (不合格の列がすべて STALLED/plateau、総合は stalled/plateau) / fail (それ以外の理由の列がある) /
            undeterminable (区間・全体行・必須列が読めない) / missing。前提 P1 は pass か stall_only。"""
    if text is None:
        return {"status": "missing", "reasons": [f"{SEGMENT_VERDICT_FILE} が無い"], "columns": {}}
    lines = text.splitlines()
    out = {"status": None, "head": None, "segment": None, "columns": {}, "reasons": []}
    seg = [re.search(r"\[segment\] 判定区間 = (.+?)\s+\(", l) for l in lines]
    seg = [m.group(1).strip() for m in seg if m]
    if len(seg) != 1:
        out.update(status="undeterminable", reasons=[f"--segment の判定区間の行が {len(seg)} 行 (1 行であること)"])
        return out
    out["segment"] = seg[0]
    heads = [i for i, l in enumerate(lines) if l.startswith("=== ") and "-> " in l]
    if len(heads) != 1:
        out.update(status="undeterminable", reasons=[f"判定の全体行が {len(heads)} 行 (1 run 1 行であること)"])
        return out
    h = heads[0]
    out["head"] = lines[h].strip()
    m = re.search(r"\[last step (-?\d+)\]", lines[h])
    out["last_step"] = int(m.group(1)) if m else None
    overall = lines[h].split("-> ", 1)[1].rstrip(" =")
    for l in lines[h + 1:]:
        if l.startswith("OVERALL") or l.startswith("=== "):
            break
        m = re.match(r"^  (\S+)\s*: (.*)$", l)
        if not m:
            continue
        col, msg = m.group(1), m.group(2).strip()
        if "<--" in msg:
            mark = msg.split("<--", 1)[1].strip()
            cls = "stalled" if mark.startswith("STALLED (plateau)") else "other"
        else:
            cls = "ok"
        out["columns"][col] = {"msg": msg, "class": cls}
    cols = out["columns"]
    R = out["reasons"]
    if seg[0] != SEGMENT_MAIN:
        R.append(f"判定区間が {seg[0]!r} (本段だけの区間 {SEGMENT_MAIN!r} でない)")
    miss = [c for c in REQUIRED_RESID_COLS if c not in cols]
    if miss:
        R.append(f"必須の残差列の行が無い: {miss}")
    if R:
        out["status"] = "undeterminable"
        return out
    other = {c: v["msg"] for c, v in cols.items() if v["class"] == "other"}
    stalled = [c for c, v in cols.items() if v["class"] == "stalled"]
    if other:
        out["status"] = "fail"
        R.extend(f"{c}: {msg}" for c, msg in other.items())
    elif overall.startswith("PASS") and not stalled:
        out["status"] = "pass"
    elif overall.startswith("NOT CONVERGED (stalled/plateau") and stalled:
        out["status"] = "stall_only"
        out["stalled_columns"] = stalled
    else:
        out["status"] = "undeterminable"
        R.append(f"全体行 {overall!r} と列の内訳 (停滞 {len(stalled)} 列・他の不合格 0 列) が食い違う")
    return out


def convergence_ok(parsed: dict) -> bool:
    return parsed.get("status") in ("pass", "stall_only")


def convergence_label(parsed: dict) -> str:
    """記録用の表示。停滞のみを PASS と書かない。"""
    st = parsed.get("status")
    return {"pass": "PASS", "stall_only": "停滞のみ (NOT CONVERGED、全不合格列が STALLED/plateau; PASS ではない)",
            "fail": "不成立 (停滞以外の理由の列がある)", "undeterminable": "判定不能", "missing": "判定ファイルなし"}.get(st, str(st))


# --- 前提 P2: 主指標の時系列 -------------------------------------------------------------------------------------------------
def qs_columns() -> list:
    return [f"{s}_{r}" for r in REGION_KEYS for s in QS_STATS]


def all_columns() -> list:
    return [f"{s}_{r}" for r in REGION_KEYS for s in QS_STATS + FRAC_STATS + ("n",)]


def series_row(stage: str, step: int, ind: dict) -> dict:
    row = {"stage": stage, "step": int(step)}
    for r in REGION_KEYS:
        for s in QS_STATS + FRAC_STATS + ("n",):
            row[f"{s}_{r}"] = ind[r][s]
    return row


def write_series_csv(path: Path, rows: list, cols: list) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["stage", "step"] + cols)
        for r in rows:
            w.writerow([r["stage"], r["step"]] + [("" if r.get(c) is None else f"{r[c]:.10g}") for c in cols])


def read_csv_cols(path: Path, cols: list) -> tuple:
    """窓の CSV を読み直す (check_quasisteady が読む値そのもので幅を測る)。戻り (steps, {列: 値の配列})。空欄は NaN。"""
    with open(path) as f:
        rows = list(csv.DictReader(f))
    steps = [int(float(r["step"])) for r in rows]
    vals = {c: np.array([float(r[c]) if r.get(c) not in (None, "", "None") else math.nan for r in rows]) for c in cols}
    return steps, vals


def run_quasisteady(csv_path: Path, cols: list, out_txt: Path) -> str:
    cmd = [sys.executable, str(TOOLS / "check_quasisteady.py"), "--series-csv", str(csv_path), "--series-cols", ",".join(cols),
           "--tail", f"{QS_TAIL:g}"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    text = f"$ {' '.join(cmd[1:])}\n(exit {r.returncode})\n" + r.stdout + r.stderr
    out_txt.write_text(text)
    return text


def parse_quasisteady(text: str, cols) -> dict:
    """check_quasisteady --series-csv の出力 → {列: VERDICT} (行が無ければ "UNKNOWN")。moc_v5_euler_eval.parse_quasisteady と同じ読み方。"""
    out = {}
    for c in cols:
        m = re.search(rf"^\s+{re.escape(c)}\s*:.*\s(\S+)\s*$", text, re.M)
        out[c] = m.group(1) if m else "UNKNOWN"
    return out


def window_check(steps, vals: dict, qs_text: str, cols, want_steps) -> dict:
    """1 窓の前提 P2 (再改訂 a59e61b5): 時点がそろい、各列の全値が有限、かつ時間方向の最大 − 最小 ≤ RANGE_TOL_K。
    check_quasisteady の VERDICT は列ごとに併記するだけ (STEADY を必須にしない)。"""
    ver = parse_quasisteady(qs_text, cols)
    rec = {"steps": list(steps), "columns": {}, "problems": []}
    if list(steps) != list(want_steps):
        rec["problems"].append(f"窓の時点 {list(steps)} が登録 {list(want_steps)} と違う")
    m = re.search(r"\[(\d+) rows, steps", qs_text or "")
    if not m or int(m.group(1)) != len(want_steps):
        rec["problems"].append(f"check_quasisteady の行数 {m.group(1) if m else '読めない'} が {len(want_steps)} でない")
    for c in cols:
        v = np.asarray(vals[c], dtype=np.float64)
        fin = bool(np.all(np.isfinite(v))) and v.size > 0
        rng = float(v.max() - v.min()) if fin else None
        ok = fin and rng <= RANGE_TOL_K
        rec["columns"][c] = {"verdict_quasisteady_recorded": ver[c], "finite": fin, "range_K": rng, "ok": bool(ok)}
    bad = [c for c, v in rec["columns"].items() if not v["ok"]]
    rec["n_bad"] = len(bad)
    rec["bad"] = bad
    rec["ok"] = (not rec["problems"]) and not bad
    rec["verdict_counts_recorded"] = {v: sum(1 for c in cols if ver[c] == v) for v in sorted(set(ver.values()))}
    rec["note"] = P2_NOTE
    return rec


# --- 判定 ----------------------------------------------------------------------------------------------------------------------
def anomaly_K(ind: dict) -> float:
    """スロートの区間の主指標の max(|最大|, |最小|) (3 領域の最大)。"""
    return max(max(abs(ind[r]["max"]), abs(ind[r]["min"])) for r in THROAT_REGIONS)


def within_K(ind: dict) -> float:
    """全 9 領域の主指標の max(最大, −最小) (±1 K 以内の判定に使う)。"""
    return max(max(ind[r]["max"], -ind[r]["min"]) for r in REGION_KEYS)


def max_diff_K(indA: dict, indB: dict) -> float:
    return max(abs(indA[r][s] - indB[r][s]) for r in REGION_KEYS for s in QS_STATS)


def _complete(ind: dict) -> bool:
    try:
        return all(ind[r]["n"] > 0 and all(ind[r][s] is not None and math.isfinite(ind[r][s]) for s in QS_STATS) for r in REGION_KEYS)
    except (KeyError, TypeError):
        return False


def judge(ind: dict, pre_ok: bool, pre_reasons=()) -> dict:
    """登録の 3 区分。ind = {"A": {step: 主指標}, "B": {step: 主指標}} (評価窓 WIN13 の 13 枚)。pre_ok = 前提と固定の条件の成立。"""
    facts = {"per_step": {}}
    missing = [f"{a}:{s}" for a in ("A", "B") for s in WIN13 if s not in ind.get(a, {}) or not _complete(ind[a][s])]
    if missing:
        return {"verdict": LBL_UNDET, "reasons": list(pre_reasons) + [f"評価窓の主指標が欠ける・非有限: {missing}"], "facts": facts}
    for s in WIN13:
        a, b = ind["A"][s], ind["B"][s]
        facts["per_step"][s] = {"A_throat_anomaly_K": anomaly_K(a), "B_throat_anomaly_K": anomaly_K(b),
                                "A_within_K": within_K(a), "B_within_K": within_K(b), "max_abs_A_minus_B_K": max_diff_K(a, b)}
    ps = facts["per_step"].values()
    A_anom = all(v["A_throat_anomaly_K"] > ANOMALY_K for v in ps)
    B_anom = all(v["B_throat_anomaly_K"] > ANOMALY_K for v in ps)
    B_within = all(v["B_within_K"] <= WITHIN_K for v in ps)
    diff_ok = all(v["max_abs_A_minus_B_K"] <= DIFF_K for v in ps)
    facts.update(A_throat_anomaly_all_steps=A_anom, B_throat_anomaly_all_steps=B_anom, B_within_1K_all_steps=B_within,
                 A_minus_B_within_1K_all_steps=diff_ok,
                 A_throat_anomaly_min_over_steps_K=min(v["A_throat_anomaly_K"] for v in ps),
                 B_within_max_over_steps_K=max(v["B_within_K"] for v in ps),
                 max_abs_A_minus_B_over_steps_K=max(v["max_abs_A_minus_B_K"] for v in ps))
    if not pre_ok:
        return {"verdict": LBL_UNDET, "reasons": ["前提の未達"] + list(pre_reasons), "facts": facts}
    if A_anom and B_within:
        return {"verdict": LBL_SUPPORT, "reasons": [], "facts": facts}
    if A_anom and B_anom and diff_ok:
        return {"verdict": LBL_REJECT, "reasons": [], "facts": facts}
    why = []
    if not A_anom:
        why.append(f"A のスロートの区間で > {ANOMALY_K:g} K の異常が評価窓の全時点で再現しない")
    else:
        if not B_within:
            why.append(f"B が評価窓の全時点で全領域 ±{WITHIN_K:g} K 以内を保たない")
        if not (B_anom and diff_ok):
            why.append(f"両側の異常の残存かつ全領域の主指標の差 ≤ {DIFF_K:g} K が評価窓の全時点で成り立たない (中間的)")
    return {"verdict": LBL_UNDET, "reasons": why, "facts": facts}


# --- 保存場の読み (試験では使わない) ------------------------------------------------------------------------------------------------
def load_geometry(run_dir: Path) -> tuple:
    """X, R [r_t 単位] (ni, nj) と S [m]。節点の並び (i·nj + j、j = nj−1 が壁) は VIZMESH/CONNE で確かめる (ic_index_map.structured_shape)。"""
    import h5py
    sys.path.insert(0, str(HERE))
    from ic_index_map import structured_shape
    info = json.loads((run_dir / "prepare_info.json").read_text())
    S = float(info["scale_m"])
    with h5py.File(run_dir / "nozzle.h5", "r") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3).astype(np.float64)
        ni, nj = structured_shape(f["VIZMESH/CONNE"], nc.shape[0])
    if ni != int(info["mesh"]["ni"]):
        raise ValueError(f"{run_dir.name}: ni {ni} が prepare_info の {info['mesh']['ni']} と違う")
    return (nc[:, 0] / S).reshape(ni, nj), (nc[:, 1] / S).reshape(ni, nj), S


def eval_snapshot(run_dir: Path, h5: Path, X, R, eta, ep: dict) -> dict:
    """1 時点: 全温の復元 (E1 と同じ) → 主指標・補助。問題は problems に入れる (既定値で通さない)。"""
    rec = {"file": str(h5.relative_to(run_dir.parent)), "problems": []}
    if not h5.is_file():
        rec["problems"].append("ファイルが無い")
        return rec
    rec["sha256"] = _sha(h5)
    try:
        sys.path.insert(0, str(HERE))
        import euler_t0_stage_ab as E1
        r = E1.recon(run_dir, h5)
        T0 = np.asarray(r["T0"], dtype=np.float64)
        if T0.size != X.size:
            raise ValueError(f"節点数 {T0.size} が格子 {X.size} と違う")
        dev = (T0 - TT_REG).reshape(X.shape)
        gas, Y = r["gas"], r["Y"]
        full = np.full(T0.shape, TT_REG)
        dh = (np.asarray(r["h0"], dtype=np.float64) - gas.h(Y, full)).reshape(X.shape)
        cp = np.asarray(gas.cp(Y, full), dtype=np.float64).reshape(X.shape)
    except Exception as e:  # noqa: BLE001 — 読めない・復元できないことも記録して判別不能へ
        rec["problems"].append(f"{type(e).__name__}: {e}")
        return rec
    n_bad = int(np.sum(~np.isfinite(dev)))
    if n_bad:
        rec["problems"].append(f"全温が非有限の節点 {n_bad}")
    if not r["newton_ok"]:
        rec["problems"].append("復元の Newton が収束しない節点がある")
    rec["thermo_source"] = r["src"]
    rec["newton_ok"] = bool(r["newton_ok"])
    if n_bad:
        return rec
    rec["indicators"] = indicators(to_points(dev, eta, ep), ep["region"])
    rec["raw"] = raw_stats(dev, X, R, eta)
    rec["h0_minus_hmix_Tt"] = h0_stats(dh, cp, X, eta)
    return rec


def _psha(p: Path):
    try:
        for l in p.read_text().splitlines():
            if l.startswith("forge_sha256"):
                return l.split(":", 1)[1].strip()
    except OSError:
        return None
    return None


def fixed_conditions(case: Path, geom_ok: bool, geom_why: str | None) -> dict:
    """固定の条件の確認 (登録「固定」)。"""
    sys.path.insert(0, str(HERE))
    probs = []
    ref = _psha(case / REF_RUN_BIN / "RUN_PROVENANCE.txt")
    shas = {}
    for arm, rn in RUNS.items():
        rd = case / rn
        shas[f"{arm}_main"] = _psha(rd / "RUN_PROVENANCE.txt")
        shas[f"{arm}_soft"] = _psha(rd / SOFT_DIR / "RUN_PROVENANCE.txt")
    if not ref:
        probs.append(f"{REF_RUN_BIN}/RUN_PROVENANCE.txt の forge_sha256 を読めない")
    bad = {k: v for k, v in shas.items() if not v or v != ref}
    if bad:
        probs.append(f"forge の sha256 が {REF_RUN_BIN} ({ref}) と違う・読めない: {bad}")
    stages, ic, tts = {}, {}, {}
    for arm, rn in RUNS.items():
        rd = case / rn
        try:
            man = json.loads((rd / "stage_manifest.json").read_text())
            st = man.get("stages") if isinstance(man, dict) else man
            stages[arm] = [s.get("tag") for s in st]
        except (OSError, ValueError, AttributeError, TypeError) as e:
            stages[arm] = f"読めない: {type(e).__name__}"
        if stages[arm] != STAGES_EXPECTED:
            probs.append(f"{arm}: stage_manifest の段が {stages[arm]} ({STAGES_EXPECTED} でない)")
        try:
            pr = json.loads((rd / "E2_PREP.json").read_text())
            ic[arm] = (pr.get("ic_check") or {}).get("VERDICT")
        except (OSError, ValueError):
            ic[arm] = None
        if ic[arm] != "OK":
            probs.append(f"{arm}: 起動前の IC の検査 (E2_PREP.json) が OK でない ({ic[arm]!r})")
        try:
            import moc_v5c_thermo_ab as V
            tts[arm] = V.read_Tt(rd)
        except Exception as e:  # noqa: BLE001
            tts[arm] = f"読めない: {type(e).__name__}"
        if tts[arm] != TT_REG:
            probs.append(f"{arm}: 入口の Tt が {tts[arm]!r} ({TT_REG} でない)")
    if not geom_ok:
        probs.append(f"A と B の格子の照合が不成立: {geom_why}")
    return {"ok": not probs, "problems": probs, "forge_sha256": shas, "forge_sha256_ref": ref, "stages": stages,
            "ic_check": ic, "Tt": tts}


def evaluate(case: Path) -> dict:
    case = Path(case).resolve()
    out = {"plan": PLAN, "plan_reg_commit": PLAN_REG_COMMIT, "plan_reg_history": list(PLAN_REG_HISTORY),
           "evaluator": Path(__file__).name, "evaluator_sha256": _sha(Path(__file__)),
           "deps_sha256": {**{n: _sha(HERE / n) for n in ("euler_t0_stage_ab.py", "moc_v5c_thermo_ab.py")},
                           **{n: _sha(TOOLS / n) for n in ("check_quasisteady.py", "check_convergence.py", "forge_species.py")}},
           "runs": RUNS, "Tt_reg": TT_REG,
           "precondition_definition": {"convergence": CONV_DEFINITION, "quasisteady": QS_DEFINITION, "p2_note": P2_NOTE,
                                       "windows": {k: list(v) for k, v in WINDOWS.items()}, "qs_tail": QS_TAIL,
                                       "range_tol_K": RANGE_TOL_K, "stats": list(QS_STATS)},
           "thresholds": {"anomaly_K": ANOMALY_K, "within_K": WITHIN_K, "diff_K": DIFF_K},
           "regions": {"x_edges_rt": list(X_EDGES_RT), "eta_edges": list(ETA_EDGES), "keys": list(REGION_KEYS)},
           "scope_note": SCOPE_NOTE}
    # 格子と評価点
    geom, geom_ok, geom_why = {}, True, None
    try:
        for arm, rn in RUNS.items():
            geom[arm] = load_geometry(case / rn)
        (XA, RA_, SA), (XB, RB_, SB) = geom["A"], geom["B"]
        if SA != SB:
            raise ValueError(f"scale_m が違う ({SA} / {SB})")
        ep = eval_points(XA, RA_, XB, RB_)
        out["eval_points"] = {"n": int(ep["eta"].size), "digest": eval_points_digest(ep),
                              "n_by_region": {k: int(np.sum(ep["region"] == i)) for i, k in enumerate(REGION_KEYS)},
                              "definition": "列ごとに A と B の節点の η の和集合、列内の η の区分線形補間"}
    except Exception as e:  # noqa: BLE001
        geom_ok, geom_why = False, f"{type(e).__name__}: {e}"
        ep = None
    out["fixed"] = fixed_conditions(case, geom_ok, geom_why)
    pre_reasons = [f"[固定] {p}" for p in out["fixed"]["problems"]]
    out["steps"] = {"evaluator": {"MAIN_NSTEPS": MAIN_NSTEPS, "OUT_INTERVAL": OUT_INTERVAL, "windows": {k: list(v) for k, v in WINDOWS.items()}},
                    "runs": {}}
    for p in window_problems():
        pre_reasons.append(f"[固定] 評価器の窓: {p}")
    ind = {"A": {}, "B": {}}
    out["per_run"] = {}
    for arm, rn in RUNS.items():
        rd = case / rn
        pr = {"run": rn, "snapshots": {}}
        out["per_run"][arm] = pr
        # P1
        segf = rd / SEGMENT_VERDICT_FILE
        conv = parse_convergence(segf.read_text() if segf.is_file() else None)
        conv["label"] = convergence_label(conv)
        pr["convergence"] = conv
        eff = run_steps(rd, conv)
        sp = [p for p in step_problems(eff) if p not in window_problems()]
        out["steps"]["runs"][arm] = {"effective": eff, "problems": sp}
        pre_reasons += [f"[固定] {arm}: {p}" for p in sp]
        if not convergence_ok(conv):
            pre_reasons.append(f"[P1] {arm}: {conv['label']} — {conv.get('reasons')}")
        if ep is None:
            pre_reasons.append(f"[P2] {arm}: 評価点が作れない ({geom_why})")
            continue
        X, R, _ = geom[arm]
        eta = column_eta(R)
        snaps = [("soft_end", SOFT_STEPS, rd / SOFT_DIR / f"res_{SOFT_STEPS}.h5")] + [("main", s, rd / f"res_{s}.h5") for s in MAIN_STEPS]
        rows = []
        for stage, step, h5 in snaps:
            rec = eval_snapshot(rd, h5, X, R, eta, ep)
            pr["snapshots"][f"{stage}:{step}"] = rec
            if rec.get("indicators") is not None and not rec["problems"]:
                rows.append(series_row(stage, step, rec["indicators"]))
                if stage == "main":
                    ind[arm][step] = rec["indicators"]
            elif stage == "main" and step in WIN13:
                pre_reasons.append(f"[P2] {arm}: 評価窓の res_{step}.h5 を評価できない ({rec['problems']})")
            else:
                pr.setdefault("record_gaps", []).append(f"{stage}:{step} {rec['problems']}")
        write_series_csv(rd / SERIES_CSV, rows, all_columns())
        pr["series_csv"] = {"file": f"{rn}/{SERIES_CSV}", "sha256": _sha(rd / SERIES_CSV)}
        # P2 (窓ごとに CSV を書いて check_quasisteady に渡し、同じ CSV の値で幅を測る)
        cols = qs_columns()
        pr["quasisteady"] = {}
        for wn, wsteps in WINDOWS.items():
            wrows = [r for r in rows if r["stage"] == "main" and r["step"] in wsteps]
            wcsv = rd / f"euler_t0_e2_series_{wn}.csv"
            write_series_csv(wcsv, wrows, cols)
            txt = run_quasisteady(wcsv, cols, rd / f"euler_t0_e2_qs_{wn}.txt")
            steps, vals = read_csv_cols(wcsv, cols)
            wc = window_check(steps, vals, txt, cols, wsteps)
            wc.update(csv=f"{rn}/{wcsv.name}", csv_sha256=_sha(wcsv), qs_output=f"{rn}/euler_t0_e2_qs_{wn}.txt")
            pr["quasisteady"][wn] = wc
            if not wc["ok"]:
                pre_reasons.append(f"[P2] {arm} {wn}: {wc['problems']} 不成立の列 {wc['n_bad']}/{len(cols)} "
                                   f"(例 {[(c, wc['columns'][c]['range_K']) for c in wc['bad'][:3]]})")
    pre_ok = not pre_reasons
    out["precondition"] = {"ok": pre_ok, "reasons": pre_reasons,
                           "convergence_labels": {a: out["per_run"][a]["convergence"]["label"] for a in out["per_run"]}}
    out["judgment"] = judge(ind, pre_ok, pre_reasons)
    out["VERDICT"] = out["judgment"]["verdict"]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("case", nargs="?", default=str(HERE))
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    case = Path(a.case).resolve()
    out = evaluate(case)
    op = Path(a.out) if a.out else case / OUT_JSON
    op.parent.mkdir(parents=True, exist_ok=True)
    op.write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float))
    for arm, pr in out["per_run"].items():
        qs = "  ".join(f"{w} " + ("OK" if v["ok"] else "不成立 %d 列" % v["n_bad"]) for w, v in pr.get("quasisteady", {}).items())
        print(f"{arm} {pr['run']}: 収束 {pr['convergence']['label']}  準定常 {qs}")
    f = out["judgment"]["facts"]
    if "A_throat_anomaly_min_over_steps_K" in f:
        print(f"窓 {WIN13[0]}〜{WIN13[-1]}: A スロート異常の最小 {f['A_throat_anomaly_min_over_steps_K']:.3f} K / "
              f"B 全領域の最大偏差 {f['B_within_max_over_steps_K']:.3f} K / |A − B| の最大 {f['max_abs_A_minus_B_over_steps_K']:.3f} K")
    print("前提:", "成立" if out["precondition"]["ok"] else "未達", f"(P2 は {P2_NOTE})")
    for r in out["precondition"]["reasons"][:12]:
        print("  -", r)
    print("VERDICT:", out["VERDICT"])
    for r in out["judgment"]["reasons"]:
        print("  -", r)
    print("注:", SCOPE_NOTE)
    print("→", op)
    return 0


if __name__ == "__main__":
    sys.exit(main())
