"""plan discretization-moc-axis-limit-and-corrector §6 V5c (別登録 2026-10-07、commit 4a484270、V5b の保留を受けて結果を見る前) の
0 step A/B。根拠: notes/reviews/2026-10-07-moc-v5b-interpretation-diagnose.md。CFD は回さない (保存場を読むだけ)。
V5・V5b の出力 (_band_ab/moc_v5_euler_eval.json・moc_v5b_ext_eval.json、各 run の wallfit_series_*) は書き換えない。

目的: Euler の最終場でスロート付近の壁際の全温が Tt を数百 K 超えることが、保存状態そのものにあるのか、出力・後処理 (原始量・
  全エンタルピー・熱物性の復元) の不整合なのかを判別する。比較に使う Euler 場の健全性の確認で、P 傾きの差の真因を決めるものではない。
変えるのは熱力学状態の取得経路だけ:
  経路 A: 保存された VALUE/P・T・Y*・h0。全温は solver_density_cuda/tools/total_quantities.total_state (既存の後処理そのもの)。
    P 傾きは保存された P (nozzle_report.load_field) から。
  経路 B: 同じ節点の VALUE/ro・roUx・roUy・roUz・roe・roY* (h0_includes_k = 1 のときだけ roK も) と、その run の解決済み熱物性
    (forge_species.run_thermo: total_state と同じ読み元・同じ優先順) とエンタルピーの基準 (thermoHrefTemp) から float64 で独立に復元する。
      Y = roY/ro (負は 0 にして和で規格化: ソルバの dependentVariables_d.cu と同じ)、u = roU/ro、e = roe/ro − |u|²/2
      (roe は sstEnergyIncludesK でも平均流エネルギーのまま: ransTransport_d.cu の sst_energy_k_correction)、
      T: e_mix(T) = h_mix(T) − R_mix T = e を Newton (初期値 1000 K。保存された T は使わない)、P = ρ R_mix T、
      h0 = (roe + P)/ro (+ k = max(roK/ro, 0) は h0_includes_k = 1 のときだけ。output.cpp の h0 = Ht [+ k] と同じ約束)。
      全温は h_mix(T0) = h0 を A と同じ熱物性関数 (forge_species.thermo_gas = total_quantities._TPGas) で逆算する。
    経路 B は保存された P・T・Y*・h0 の値を読まない (読むのは VALUE/h0 の属性 h0_includes_k だけ。read_conserved の読み取り一覧)。
    ρ・P の床 (physProp.roMin・pMin) は当てない (当たる節点の数を記録する)。B の場はファイルに書かない。
対象 (登録、計 56 標本): run_0157 (腕 M r1 の延長) と run_0160 (ISEN の延長) の窓 A (通算 24000〜36000 = 子の local 6000〜18000) と
  窓 B (通算 42000〜54000 = local 24000〜36000) の各 13 枚、親 run_0150・run_0153 の res_0.h5 (初期場) と res_18000.h5 (最終場)。
P 傾き: η = 0.1、V5b と同じ定義・座標・標本。eval_wallfit_euler.quantities の P_slope_eta0.1 と同じ式 (grid(0.05) の標本、
  窓 WIN_T の 1 次回帰の傾き × 窓幅、P_REF = 2237 Pa)。評価座標は _band_ab/wallfit_series_v5b.json の geometry。定義は
  moc_v5_euler_eval._wallfit_defs で eval_wallfit_euler.py から取り出す (あのモジュールは import すると旧 A/B が走る)。
  定義の再現の検査 (前提): eval_wallfit_euler.py の sha256 が V5b の時系列の記録と同じ / 経路 A の値が quantities の値と同じ
  (差 ≤ 1e-12) / 時系列の CSV (sha256 が記録と同じもの) に同じ step があれば、その値と同じ (差 ≤ 1e-9 %pt。CSV は 10 桁)。
記録 (標本ごと): P 傾きの A − B、全温の差 (最大・分位・位置)、h0 の整合 (保存 h0 − B の h0、保存 h0 − (roe + 保存 P)/ro)、
  保存された P・T・Y と B の差、壁からの層 (j_w = nj − 1 − j、0 = 壁の節点) ごとの全温の最大・分位値 (50・99 %)・Tt + 1 K を超える
  節点の数と最大の x (A と B それぞれ)、層ごとの |全温の A − B| の最大、全温の最大の位置 (A・B・|A − B|)、B の Newton の収束と床。
判定 (事前登録の閾値: P 傾き 0.003 %pt、全温 1 K。「閾値内」は |A − B| ≤ 閾値、全温は全節点の最大):
  退ける: 全標本で両経路が閾値内、かつ数百 K の超過が両経路で再現する (超過のある標本では A・B とも超過があり、そのような標本が
    1 つ以上ある) → 「出力・復元の経路だけが原因」を退け、保存状態側の診断へ進む (次の手は改めて諮る)。
  棄却: 閾値を超えて食い違う標本がある → 今の復元・出力の整合を棄却し、その修正を先に行う。両者の差を物理的な IC 依存と解釈しない。
  どちらでもない: 上のどちらでもない (例: 全標本で閾値内だが超過が数百 K に達しない・超過が片方の経路にだけある)。
  判定不能: 標本の欠損・読めない・非有限値・B の Newton の未収束・定義の再現の不成立・熱物性の読み元が A と B で違う。
  どの判定にも、超過が両経路・A だけ・B だけ・どちらにも無い標本の数と、超過の最大を事実として添える。
  「数百 K」は登録に数値が無いので、下限 100 K (max T0 − Tt ≥ 100 K) と読む (実装時の読み。全温の超過の最大はそのまま記録する)。
  --only で一部の標本だけ回したときは「部分標本」と明記し、登録の判定には使わない。
やらないこと (登録): 再延長、窓の移動、2SE による救済、出口較正の変更、全温のクリップ。
出力: _band_ab/moc_v5c_thermo_ab.json (--only のときの既定は _band_ab/moc_v5c_thermo_ab_partial.json)。評価器自身と依存の sha256、
  登録の commit (4a484270) を入れる。
usage:
  python3 moc_v5c_thermo_ab.py [case_dir] [--only RUN:STEP[,RUN:STEP...]] [--out PATH] [--Tt K] [--species-source auto|record|speciesDBFile]
    STEP は res_<STEP>.h5 の番号 (その run の local step)。--Tt を省くと各 run の bcondConfig.yaml の入口の Tt (全 run で同じこと)。
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TOOLS = ROOT / "solver_density_cuda/tools"
for _p in (str(HERE), str(TOOLS), str(ROOT / "design")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import moc_v5_euler_eval as EV  # noqa: E402  (eval_wallfit_euler の定義の取り出しと評価座標の設定だけを使う)
import total_quantities as TQ  # noqa: E402  (経路 A の全温)
import forge_species as FSP  # noqa: E402  (熱物性の読み元。A と B で同じ関数)
from forge_design.report.nozzle_report import eta_line, load_field  # noqa: E402

PLAN = "plans/active/discretization-moc-axis-limit-and-corrector.md §6 V5c"
PLAN_REG_COMMIT = "4a484270"                      # §6 V5c の別登録 (2026-10-07、結果を見る前) の commit
TOL_PSLOPE = 0.003                                # [%pt] P 傾きの A − B の閾値 (事前登録)
TOL_T0 = 1.0                                      # [K] 全温の A − B の閾値 (事前登録)
EXCEED_MARGIN = 1.0                               # [K] 「Tt + 1 K を超える節点」の数え方
EXCESS_HUNDREDS = 100.0                           # [K] 「数百 K の超過」の下限 (登録に数値が無いので 100 K と読む; 実装時)
ETA = 0.1
SERIES_H = 0.05                                   # [r_t] eval_wallfit_euler の時系列の評価刻み (quantities(F, 0.05))
PSLOPE_KEY = "P_slope_eta0.1"
IDENTITY_TOL = 1e-12                              # [%pt] 経路 A の P 傾きと quantities の値の差 (同じ式なので 0 のはず)
SERIES_TOL = 1e-9                                 # [%pt] 経路 A の P 傾きと時系列の CSV (10 桁) の差
NEWTON_T_INIT = 1000.0                            # [K] 経路 B の温度の Newton の初期値 (保存された T を使わない)
NEWTON_TOL = 1e-7                                 # [K] 経路 B の Newton の収束 (|ΔT|)
NEWTON_MAX_IT = 200
T_SOLVER_RANGE = (50.0, 6000.0)                   # ソルバの温度反転のクランプ (dependentVariables_d.cu DEPVAR_TMIN/TMAX)。B はクランプしない (外れた数を記録)
PARENT_END = 18000                                # 親の終了 step (子の通算 = local + 18000)
CHILDREN = {"run_0157_euler_wallfit_mocG1_r1_ext36k": "run_0150_euler_wallfit_mocG1_r1",
            "run_0160_euler_icdep_mocG1_isen_ext36k": "run_0153_euler_icdep_mocG1_isen"}
WIN_A_LOCAL = list(range(6000, 18001, 1000))      # 通算 24000〜36000 (13 枚)
WIN_B_LOCAL = list(range(24000, 36001, 1000))     # 通算 42000〜54000 (13 枚)
PARENT_STEPS = (0, 18000)                         # 親の初期場と最終場
GEOM_RECORD = "_band_ab/wallfit_series_v5b.json"  # V5b の時系列の記録 (評価座標と評価器の sha256)
SERIES_RECORDS = {"v5b": "_band_ab/wallfit_series_v5b.json", "v5": "_band_ab/wallfit_series_v5.json"}
OUT_JSON, OUT_JSON_PARTIAL = "_band_ab/moc_v5c_thermo_ab.json", "_band_ab/moc_v5c_thermo_ab_partial.json"
B_READS = ("ro", "roUx", "roUy", "roUz", "roe")   # 経路 B が値を読む保存量 (+ roY{i}、h0_includes_k = 1 なら roK)
LBL_REJECT_PATH = "退ける: 出力・復元の経路だけが原因"
LBL_REJECT_CONS = "棄却: 今の復元・出力の整合"
LBL_NEITHER = "どちらでもない"
LBL_UNDET = "判定不能"
PARTIAL_PREFIX = "部分標本 (登録の全標本でない; 登録の判定に使わない): "


def _sha(p: Path):
    p = Path(p)
    if not p.is_file():
        return None
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


# --- 標本 -------------------------------------------------------------------------------------------------------------
def registered_samples() -> list:
    """登録の 56 標本: [{run, role, local_step, total_step}]。子は窓 A・B の各 13 枚、親は初期場と最終場。"""
    out = []
    for child, parent in CHILDREN.items():
        for st in PARENT_STEPS:
            out.append({"run": parent, "role": "parent_initial" if st == 0 else "parent_final", "local_step": st, "total_step": st})
        for role, steps in (("child_winA", WIN_A_LOCAL), ("child_winB", WIN_B_LOCAL)):
            for st in steps:
                out.append({"run": child, "role": role, "local_step": st, "total_step": st + PARENT_END})
    return out


def sample_of(run: str, step: int) -> dict:
    """--only の 1 標本。登録の run なら役割と通算 step を付ける (登録外の run は role other、通算 = local)。"""
    for s in registered_samples():
        if s["run"] == run and s["local_step"] == step:
            return dict(s)
    if run in CHILDREN:
        return {"run": run, "role": "child_other", "local_step": step, "total_step": step + PARENT_END}
    return {"run": run, "role": "parent_other" if run in CHILDREN.values() else "other", "local_step": step, "total_step": step}


def parse_only(text: str) -> list:
    out = []
    for tok in [t for t in text.split(",") if t]:
        run, _, st = tok.rpartition(":")
        if not run or not st.isdigit():
            raise SystemExit(f"--only の {tok!r} は RUN:STEP の形でない")
        out.append(sample_of(run, int(st)))
    return out


# --- 経路 B: 保存量からの独立な復元 -----------------------------------------------------------------------------------------
def read_conserved(res_path: Path, n_species: int) -> dict:
    """経路 B が読むもの: VALUE の B_READS・roY{i} (・roK) の値と、VALUE/h0 の属性 h0_includes_k だけ (P・T・Y*・h0 の値は読まない)。
    欠けた保存量は KeyError (既定値で埋めない)。単一種で roY0 が無ければ Y = 1。"""
    with h5py.File(res_path, "r") as f:
        if "VALUE/h0" not in f:
            raise KeyError(f"{res_path.name}: VALUE/h0 が無い (属性 h0_includes_k を読めない)")
        inc_k = int(f["VALUE/h0"].attrs.get("h0_includes_k", 0))
        need = list(B_READS) + (["roK"] if inc_k else [])
        miss = [k for k in need if f"VALUE/{k}" not in f]
        if miss:
            raise KeyError(f"{res_path.name}: 保存量 {miss} が無い")
        d = {k: f[f"VALUE/{k}"][:].astype(np.float64) for k in need}
        if n_species == 1 and "VALUE/roY0" not in f:
            d["roY"] = None
        else:
            miss = [f"roY{i}" for i in range(n_species) if f"VALUE/roY{i}" not in f]
            if miss:
                raise KeyError(f"{res_path.name}: 保存量 {miss} が無い")
            d["roY"] = [f[f"VALUE/roY{i}"][:].astype(np.float64) for i in range(n_species)]
    d["includes_k"] = inc_k
    return d


def newton_T(fun, dfun, target, T_init, tol=NEWTON_TOL, max_it=NEWTON_MAX_IT):
    """fun(T) = target を節点ごとに Newton で解く (1 反復の刻みは ±0.5 T まで、T ≥ 1 K)。
    戻り値 (T, 最後の |ΔT|, 反復回数)。未収束は呼び出し側が最後の |ΔT| で数える。"""
    T = np.array(T_init, dtype=np.float64, copy=True)
    step = np.full(T.shape, np.inf)
    it = 0
    for it in range(1, max_it + 1):
        dT = (fun(T) - target) / dfun(T)
        dT = np.clip(dT, -0.5 * T, 0.5 * T)
        T = np.maximum(T - dT, 1.0)
        step = np.abs(dT)
        if np.all(step < tol):
            break
    return T, step, it


def path_b(run_dir: Path, res_path: Path, species_source: str = "auto") -> dict:
    """経路 B: 保存量と解決済み熱物性から float64 で原始量・全エンタルピー・全温を復元する (保存された P・T・Y*・h0 は使わない)。"""
    th = FSP.run_thermo(str(run_dir), res_path=str(res_path), source=species_source)
    if th is None:
        raise ValueError("V5c は TP (thermalMethod 2) の run が対象 (CPG は経路 B を実装していない)")
    gas = FSP.thermo_gas(th)
    names = list(th["names"])
    c = read_conserved(Path(res_path), len(names))
    ro = c["ro"]
    if c["roY"] is None:
        Yraw = [np.ones_like(ro)]
    else:
        Yraw = [r / ro for r in c["roY"]]
    ysum_raw = np.sum(Yraw, axis=0)
    Yc = [np.where(y < 0.0, 0.0, y) for y in Yraw]
    ysum = np.sum(Yc, axis=0)
    Y = [y / np.maximum(ysum, 1e-30) for y in Yc]
    ux, uy, uz = c["roUx"] / ro, c["roUy"] / ro, c["roUz"] / ro
    ek = 0.5 * (ux * ux + uy * uy + uz * uz)
    e = c["roe"] / ro - ek
    Rm = gas.Rmix(Y)
    T, stepT, itT = newton_T(lambda t: gas.h(Y, t) - Rm * t, lambda t: np.maximum(gas.cp(Y, t) - Rm, 1e-2 * Rm), e,
                             np.full(ro.shape, NEWTON_T_INIT))
    P = ro * Rm * T
    k = np.maximum(c["roK"] / ro, 0.0) if c["includes_k"] else np.zeros_like(ro)
    h0 = (c["roe"] + P) / ro + k
    T0, stepT0, itT0 = newton_T(lambda t: gas.h(Y, t), lambda t: np.maximum(gas.cp(Y, t), 1.0), h0, T)
    resid_e = np.abs(gas.h(Y, T) - Rm * T - e) / (gas.cp(Y, T) - Rm)
    resid_h = np.abs(gas.h(Y, T0) - h0) / gas.cp(Y, T0)
    cfg = yaml.safe_load((Path(run_dir) / "solverConfig.yaml").read_text())
    pp = cfg.get("physProp") or {}
    roMin, pMin = float(pp.get("roMin", 1.0e-4)), float(pp.get("pMin", 1.0))
    src = f"{th['source']} ({th['how']}{': ' + Path(th['path']).name if th['path'] else ''})"
    diag = {"newton_T_iterations": int(itT), "newton_T_last_step_max_K": float(np.max(stepT)),
            "newton_T_n_nonconverged": int(np.sum(~(stepT < NEWTON_TOL))), "newton_T_resid_max_K": float(np.max(resid_e)),
            "newton_T0_iterations": int(itT0), "newton_T0_last_step_max_K": float(np.max(stepT0)),
            "newton_T0_n_nonconverged": int(np.sum(~(stepT0 < NEWTON_TOL))), "newton_T0_resid_max_K": float(np.max(resid_h)),
            "sumY_raw_minus_1_max_abs": float(np.max(np.abs(ysum_raw - 1.0))), "n_negative_Y": int(sum(np.sum(y < 0.0) for y in Yraw)),
            "n_ro_below_roMin": int(np.sum(ro < roMin)), "n_P_below_pMin": int(np.sum(P < pMin)),
            "n_T_outside_solver_range": int(np.sum((T < T_SOLVER_RANGE[0]) | (T > T_SOLVER_RANGE[1]))),
            "roMin": roMin, "pMin": pMin,
            "n_nonfinite_conserved": int(sum(int(np.sum(~np.isfinite(v))) for v in
                                             [c[k_] for k_ in c if isinstance(c[k_], np.ndarray)] + list(c["roY"] or [])))}
    return {"T": T, "P": P, "h0": h0, "T0": T0, "Y": Y, "ek": ek, "k": k, "e": e, "includes_k": c["includes_k"],
            "names": names, "gas": gas, "thermo_source": src, "diag": diag, "reads": list(B_READS) + ["roY*"]
            + (["roK"] if c["includes_k"] else []) + ["VALUE/h0 の属性 h0_includes_k"]}


def load_geometry(run_dir: Path) -> tuple:
    """経路 B の座標: prepare_info.json の scale_m・mesh.ni と nozzle.h5 の /MESH/COORD (nozzle_report.load_field と同じ取り方)。"""
    info = json.loads((Path(run_dir) / "prepare_info.json").read_text())
    S = float(info["scale_m"]); ni = int(info["mesh"]["ni"])
    with h5py.File(Path(run_dir) / "nozzle.h5", "r") as f:
        nc = f["/MESH/COORD"][:].reshape(-1, 3)
    nj = nc.shape[0] // ni
    return (nc[:, 0] / S).reshape(ni, nj), (nc[:, 1] / S).reshape(ni, nj), ni, nj


# --- P 傾き (V5b と同じ定義) ------------------------------------------------------------------------------------------------
class _NpCompat:
    """numpy 1 には np.trapezoid が無い (quantities の出口コア M が使う)。定義の一致の検査で quantities を呼ぶときだけ np.trapz で補う。"""

    def __getattr__(self, k):
        if k == "trapezoid":
            return getattr(np, "trapezoid", None) or np.trapz
        return getattr(np, k)


def wallfit_ns(geom: dict) -> dict:
    """eval_wallfit_euler.py の grid・quantities と P_REF を取り出し、評価座標を V5b の geometry に固定する。"""
    ns = EV._wallfit_defs()
    ns["eta_line"] = eta_line
    if not hasattr(np, "trapezoid"):
        ns["np"] = _NpCompat()
    EV._set_geometry(ns, geom)
    return ns


def p_slope(ns: dict, F: dict, eta: float = ETA) -> float:
    """eval_wallfit_euler.quantities の P_slope_eta<η> と同じ式 (grid(0.05) の標本、窓 WIN_T の 1 次回帰の傾き × 窓幅)。F は X・R・V["P"] だけ使う。"""
    xq = ns["grid"](SERIES_H)
    WIN_T = ns["WIN_T"]
    wt = (xq >= WIN_T[0]) & (xq <= WIN_T[1])
    dp = 100 * (ns["eta_line"](F, "P", eta, xq) / ns["P_REF"] - 1)
    cf = np.polyfit(xq[wt], dp[wt], 1)
    return float(cf[0] * (WIN_T[1] - WIN_T[0]))


def series_value(case: Path, smp: dict) -> dict:
    """記録用の照合: V5b (子) / V5 (親) の時系列の CSV の同じ local step の P_slope_eta0.1。CSV の sha256 が時系列の記録と同じときだけ使う。"""
    tag = "v5b" if smp["run"] in CHILDREN else "v5"
    rec_p = case / SERIES_RECORDS[tag]
    try:
        want = ((json.loads(rec_p.read_text()).get("runs") or {}).get(smp["run"]) or {}).get("csv_sha256")
    except (OSError, ValueError):
        want = None
    for p in (case / smp["run"] / f"wallfit_series_{tag}.csv", case / "_band_ab" / f"{tag}_series" / smp["run"] / f"wallfit_series_{tag}.csv"):
        if not p.is_file():
            continue
        sha = _sha(p)
        out = {"file": str(p.relative_to(case)), "sha256": sha, "sha256_matches_record": bool(want) and sha == want}
        if not out["sha256_matches_record"]:
            out["status"] = "記録と違う CSV (照合に使わない)"
            return out
        rows = [r for r in csv.DictReader(open(p)) if r.get("step") not in (None, "") and int(float(r["step"])) == smp["local_step"]]
        if not rows:
            out["status"] = f"step {smp['local_step']} の行が無い (照合なし)"
            return out
        out.update(status="ok", value=float(rows[0][PSLOPE_KEY]))
        return out
    return {"status": "CSV が無い (照合なし)"}


# --- 集計 -------------------------------------------------------------------------------------------------------------------
def _stats(a) -> dict:
    """符号付きの差の統計 (|a| の分位を含む)。非有限があれば数だけ。"""
    a = np.asarray(a, dtype=np.float64).ravel()
    bad = int(np.sum(~np.isfinite(a)))
    if bad:
        return {"n_nonfinite": bad}
    ab = np.abs(a)
    q = np.quantile(ab, [0.5, 0.9, 0.99, 0.999])
    return {"max_abs": float(ab.max()), "min": float(a.min()), "max": float(a.max()),
            "q50_abs": float(q[0]), "q90_abs": float(q[1]), "q99_abs": float(q[2]), "q999_abs": float(q[3]), "n_nonfinite": 0}


def _loc(idx: int, X, R) -> dict:
    ni, nj = X.shape
    i, j = divmod(int(idx), nj)
    return {"i": i, "j": j, "j_from_wall": nj - 1 - j, "x_rt": float(X[i, j]), "r_rt": float(R[i, j]),
            "eta": float(R[i, j] / R[i, -1]) if R[i, -1] != 0 else None}


def layer_stats(T0: np.ndarray, X, Tt: float) -> dict:
    """壁からの層 (j_w = 0 が壁の節点) ごとの全温の最大・50 %・99 % 分位・Tt + 1 K を超える節点の数・最大の x。リストは j_w の順。"""
    G = T0.reshape(X.shape)
    imax = np.argmax(G, axis=0)
    q = np.quantile(G, [0.5, 0.99], axis=0)
    rev = slice(None, None, -1)
    return {"max": G.max(axis=0)[rev].tolist(), "q50": q[0][rev].tolist(), "q99": q[1][rev].tolist(),
            "n_exceed": (G > Tt + EXCEED_MARGIN).sum(axis=0)[rev].astype(int).tolist(),
            "x_at_max": X[imax, np.arange(X.shape[1])][rev].tolist()}


def read_Tt(run_dir: Path):
    """bcondConfig.yaml の入口 (kind が inlet で始まる最初の境界) の Tt (total_quantities の main と同じ探し方)。無ければ None。"""
    p = Path(run_dir) / "bcondConfig.yaml"
    if not p.is_file():
        return None
    bc = yaml.safe_load(p.read_text()) or {}
    for v in bc.values():
        if isinstance(v, dict) and str(v.get("kind", "")).startswith("inlet") and "Tt" in (v.get("floats") or {}):
            return float(v["floats"]["Tt"])
    return None


def evaluate_sample(case: Path, smp: dict, ns: dict, Tt: float, species_source: str) -> dict:
    """1 標本の A/B。問題があれば problems に入れて status を ok 以外にする (既定値で通さない)。"""
    run, step = smp["run"], smp["local_step"]
    rd = case / run
    res = rd / f"res_{step}.h5"
    out = {**smp, "res": f"{run}/res_{step}.h5", "problems": []}
    if not res.is_file():
        out.update(status="missing", problems=[f"{out['res']} が無い"])
        return out
    out["res_sha256"] = _sha(res)
    try:
        # 経路 A: 保存された P・T・Y*・h0 と既存の後処理
        stA = TQ.total_state(str(rd), str(res), species_source=species_source)
        FA = load_field(rd, res.name)
        T0A = np.asarray(stA["T0"], dtype=np.float64)
        PA, TA = np.asarray(stA["P"], dtype=np.float64), np.asarray(stA["T"], dtype=np.float64)
        sA = p_slope(ns, FA)
        sA_q = float(ns["quantities"](FA, SERIES_H)[0][PSLOPE_KEY])
        # 経路 B: 保存量からの独立な復元 (上の A の値は渡さない)
        B = path_b(rd, res, species_source)
        XB, RB, ni, nj = load_geometry(rd)
        if XB.size != B["T"].size:
            raise ValueError(f"節点数 {B['T'].size} が格子 ni·nj = {XB.size} と違う")
        sB = p_slope(ns, {"X": XB, "R": RB, "V": {"P": B["P"].reshape(ni, nj)}})
        # 比較 (ここで初めて保存された h0・Y・roe を読む: 整合の記録用)
        with h5py.File(res, "r") as f:
            h0S = f["VALUE/h0"][:].astype(np.float64)
            roe, ro = f["VALUE/roe"][:].astype(np.float64), f["VALUE/ro"][:].astype(np.float64)
            YS = [f[f"VALUE/Y{i}"][:].astype(np.float64) for i in range(len(B["names"]))] if "VALUE/Y0" in f else None
    except (Exception, SystemExit) as e:  # noqa: BLE001 — total_state は h0 の欠損を SystemExit で返す
        out.update(status="error", problems=[f"{type(e).__name__}: {e}"])
        return out
    gas = B["gas"]
    cpB0 = gas.cp(B["Y"], B["T0"])
    dT0 = T0A - B["T0"]
    h0_cons_S = (roe + PA) / ro + B["k"]
    YSl = YS if YS is not None else [np.ones_like(ro)]
    residA = np.abs(gas.h(YSl, T0A) - h0S) / gas.cp(YSl, T0A)
    X, R = FA["X"], FA["R"]
    excA, excB = float(T0A.max() - Tt), float(B["T0"].max() - Tt)
    ser = series_value(case, smp)
    out.update(
        thermo_source_A=stA.get("species_source"), thermo_source_B=B["thermo_source"], h0_includes_k=int(stA["includes_k"]),
        b_reads=B["reads"], geometry_B_equals_A=bool(np.array_equal(XB, X) and np.array_equal(RB, R)),
        P_slope={"eta": ETA, "A": sA, "B": sB, "A_minus_B": sA - sB, "absA_minus_absB": abs(sA) - abs(sB),
                 "A_via_quantities": sA_q, "A_identity_diff": sA - sA_q, "series": ser,
                 "A_minus_series": (sA - ser["value"]) if ser.get("status") == "ok" else None},
        T0={"Tt": Tt, "A_max": float(T0A.max()), "B_max": float(B["T0"].max()), "A_min": float(T0A.min()), "B_min": float(B["T0"].min()),
            "A_max_minus_Tt": excA, "B_max_minus_Tt": excB,
            "A_n_exceed": int(np.sum(T0A > Tt + EXCEED_MARGIN)), "B_n_exceed": int(np.sum(B["T0"] > Tt + EXCEED_MARGIN)),
            "A_minus_B": _stats(dT0), "argmax_A": _loc(np.argmax(T0A), X, R), "argmax_B": _loc(np.argmax(B["T0"]), X, R),
            "argmax_abs_A_minus_B": _loc(np.argmax(np.abs(dT0)), X, R) if np.all(np.isfinite(dT0)) else None,
            "A_newton_resid_max_K": float(np.max(residA))},
        h0={"stored_minus_B_Jkg": _stats(h0S - B["h0"]), "stored_minus_B_K": _stats((h0S - B["h0"]) / cpB0),
            "stored_minus_roeP_stored_Jkg": _stats(h0S - h0_cons_S), "stored_minus_roeP_stored_K": _stats((h0S - h0_cons_S) / cpB0)},
        P={"stored_minus_B_rel": _stats((PA - B["P"]) / B["P"]), "stored_minus_B_Pa": _stats(PA - B["P"])},
        T={"stored_minus_B_K": _stats(TA - B["T"])},
        Y={"stored_minus_B_max_abs": (float(max(np.max(np.abs(a - b)) for a, b in zip(YS, B["Y"]))) if YS is not None else None)},
        B_diag=B["diag"],
        layers={"j_from_wall": list(range(nj)), "A": layer_stats(T0A, X, Tt), "B": layer_stats(B["T0"], X, Tt),
                "abs_A_minus_B_max": np.abs(dT0).reshape(X.shape).max(axis=0)[::-1].tolist()})
    # 前提 (標本ごと)
    pr = out["problems"]
    for nm, v in (("経路 A の全温", T0A), ("経路 B の全温", B["T0"]), ("経路 B の P", B["P"])):
        n = int(np.sum(~np.isfinite(v)))
        if n:
            pr.append(f"{nm}に非有限値が {n} 節点")
    if B["diag"]["n_nonfinite_conserved"]:
        pr.append(f"保存量に非有限値が {B['diag']['n_nonfinite_conserved']} 個")
    if not (math.isfinite(sA) and math.isfinite(sB)):
        pr.append(f"P 傾きが非有限 (A {sA}, B {sB})")
    for k_ in ("newton_T_n_nonconverged", "newton_T0_n_nonconverged"):
        if B["diag"][k_]:
            pr.append(f"経路 B の Newton が {B['diag'][k_]} 節点で未収束 ({k_}、|ΔT| ≥ {NEWTON_TOL:g} K)")
    if not abs(sA - sA_q) <= IDENTITY_TOL:
        pr.append(f"経路 A の P 傾き {sA!r} が eval_wallfit_euler.quantities の値 {sA_q!r} と違う (定義の再現が不成立)")
    if ser.get("status") == "ok" and not abs(sA - ser["value"]) <= SERIES_TOL:
        pr.append(f"経路 A の P 傾き {sA:.10g} が時系列の記録 {ser['file']} の値 {ser['value']:.10g} と違う (定義・座標の再現が不成立)")
    if out["thermo_source_A"] != out["thermo_source_B"]:
        pr.append(f"熱物性の読み元が A ({out['thermo_source_A']}) と B ({out['thermo_source_B']}) で違う")
    if not out["geometry_B_equals_A"]:
        pr.append("経路 B の座標が load_field の座標と一致しない")
    out["status"] = "ok" if not pr else "problem"
    if out["status"] == "ok":
        out["agree_P_slope"] = bool(abs(sA - sB) <= TOL_PSLOPE)
        out["agree_T0"] = bool(out["T0"]["A_minus_B"]["max_abs"] <= TOL_T0)
        out["agree"] = out["agree_P_slope"] and out["agree_T0"]
        out["excess_A"] = bool(excA >= EXCESS_HUNDREDS)
        out["excess_B"] = bool(excB >= EXCESS_HUNDREDS)
    return out


def _tag(s: dict) -> str:
    return f"{s['run']}@{s['local_step']} (通算 {s['total_step']})"


def judge(samples: list, global_problems: list, partial: bool, expected: list) -> dict:
    """事前登録の判定。標本の問題・欠損・全体の前提の不成立は判定不能。事実 (超過の内訳・最大の差) は常に添える。"""
    ok = [s for s in samples if s.get("status") == "ok"]
    bad = [s for s in samples if s.get("status") != "ok"]
    got = {(s["run"], s["local_step"]) for s in samples}
    missing = [] if partial else [e for e in expected if (e["run"], e["local_step"]) not in got]
    exA = {_tag(s) for s in ok if s["excess_A"]}
    exB = {_tag(s) for s in ok if s["excess_B"]}
    both, onlyA, onlyB = sorted(exA & exB), sorted(exA - exB), sorted(exB - exA)
    dis_p = [s for s in ok if not s["agree_P_slope"]]
    dis_t = [s for s in ok if not s["agree_T0"]]
    facts = {
        "n_samples": len(samples), "n_ok": len(ok), "n_expected": len(expected), "partial": partial,
        "n_excess_both": len(both), "n_excess_A_only": len(onlyA), "n_excess_B_only": len(onlyB),
        "n_excess_neither": len(ok) - len(both) - len(onlyA) - len(onlyB),
        "excess_A_only": onlyA, "excess_B_only": onlyB,
        "max_excess_A_K": max((s["T0"]["A_max_minus_Tt"] for s in ok), default=None),
        "max_excess_B_K": max((s["T0"]["B_max_minus_Tt"] for s in ok), default=None),
        "max_abs_P_slope_A_minus_B": max((abs(s["P_slope"]["A_minus_B"]) for s in ok), default=None),
        "max_abs_T0_A_minus_B_K": max((s["T0"]["A_minus_B"]["max_abs"] for s in ok), default=None),
        "n_disagree_P_slope": len(dis_p), "n_disagree_T0": len(dis_t),
        "disagree_P_slope": [_tag(s) for s in dis_p], "disagree_T0": [_tag(s) for s in dis_t],
        "problem_samples": {_tag(s): s.get("problems") for s in bad}, "missing_samples": [_tag(e) for e in missing],
        "global_problems": global_problems,
    }
    fact_txt = (f"超過 (max T0 − Tt ≥ {EXCESS_HUNDREDS:g} K) の標本: 両経路 {len(both)}・A だけ {len(onlyA)}・B だけ {len(onlyB)}・"
                f"どちらにも無い {facts['n_excess_neither']} (評価できた {len(ok)} 標本中)。超過の最大 A {facts['max_excess_A_K']}・"
                f"B {facts['max_excess_B_K']} K。|A − B| の最大: P 傾き {facts['max_abs_P_slope_A_minus_B']} %pt、全温 "
                f"{facts['max_abs_T0_A_minus_B_K']} K。")
    if global_problems or bad or missing:
        why = list(global_problems) + [f"{_tag(s)}: {'; '.join(s.get('problems') or [])}" for s in bad] + \
            ([f"標本の欠損 {len(missing)} 件 ({', '.join(facts['missing_samples'][:4])}{' …' if len(missing) > 4 else ''})"] if missing else [])
        label, text = LBL_UNDET, "判定不能 (前提不成立): " + " / ".join(why) + "。参考の事実: " + fact_txt
    elif dis_p or dis_t:
        side = ("超過は A にだけある標本がある。" if onlyA else "") + ("超過は B にだけある標本がある。" if onlyB else "")
        label = LBL_REJECT_CONS
        text = (f"閾値を超えて食い違う標本がある (P 傾き > {TOL_PSLOPE} %pt: {len(dis_p)} 標本、全温 > {TOL_T0:g} K: {len(dis_t)} 標本)。"
                f"今の復元・出力の整合を棄却し、その修正を先に行う。両者の差を物理的な IC 依存と解釈しない。{side}" + fact_txt)
    elif both and not onlyA and not onlyB:
        label = LBL_REJECT_PATH
        text = (f"両経路が全 {len(ok)} 標本で閾値内 (P 傾き ≤ {TOL_PSLOPE} %pt、全温 ≤ {TOL_T0:g} K) で、数百 K の超過が両経路で再現した "
                f"({len(both)} 標本)。「出力・復元の経路だけが原因」を退ける。保存状態側の診断へ進む (次の手は改めて諮る)。" + fact_txt)
    else:
        if onlyA or onlyB:
            what = (f"超過が A にだけある標本 {len(onlyA)}" if onlyA else "") + ("・" if onlyA and onlyB else "") + \
                (f"超過が B にだけある標本 {len(onlyB)}" if onlyB else "")
        else:
            what = f"どの標本にも数百 K (≥ {EXCESS_HUNDREDS:g} K) の超過が無い"
        label = LBL_NEITHER
        text = (f"どちらでもない: 両経路は全 {len(ok)} 標本で閾値内だが、{what}。「出力・復元の経路だけが原因」を退ける条件"
                f" (数百 K の超過が両経路で再現) を満たさない。" + fact_txt)
    if partial:
        label, text = PARTIAL_PREFIX + label, PARTIAL_PREFIX + text
    return {"label": label, "text": text, "facts": facts}


def build_context(case: Path) -> tuple:
    """評価座標 (V5b の記録) と eval_wallfit_euler の定義。全体の前提の問題を返す。"""
    problems, info = [], {}
    rp = case / GEOM_RECORD
    rec = None
    try:
        rec = json.loads(rp.read_text())
    except (OSError, ValueError) as e:
        problems.append(f"V5b の時系列の記録 {GEOM_RECORD} を読めない ({type(e).__name__}: {e})")
    geom = (rec or {}).get("geometry")
    info.update(geometry_record=GEOM_RECORD, geometry_record_sha256=_sha(rp), geometry=geom)
    gfail = EV.geometry_failures(geom) if rec is not None else []
    problems += [f"V5b の評価座標: {g}" for g in gfail]
    ns = None
    if rec is not None and not gfail:
        ns = wallfit_ns(geom)
        info["eval_wallfit_euler_sha256"] = ns["_source_sha256"]
        info["series_evaluator_sha256_v5b"] = rec.get("evaluator_sha256")
        if ns["_source_sha256"] != rec.get("evaluator_sha256"):
            problems.append(f"eval_wallfit_euler.py の sha256 {ns['_source_sha256'][:12]} が V5b の時系列の評価器 "
                            f"{str(rec.get('evaluator_sha256'))[:12]} と違う (P 傾きの定義が同じと言えない)")
        try:
            g5 = json.loads((case / SERIES_RECORDS["v5"]).read_text()).get("geometry")
            info["v5_geometry_equal_v5b"] = (g5 == geom)
        except (OSError, ValueError):
            info["v5_geometry_equal_v5b"] = None
    return ns, problems, info


def evaluate(case: Path, samples: list, partial: bool, Tt_override=None, species_source="auto") -> dict:
    case = Path(case).resolve()
    ns, gprob, ginfo = build_context(case)
    Tts = {}
    for run in sorted({s["run"] for s in samples}):
        Tts[run] = read_Tt(case / run)
    if Tt_override is not None:
        Tt = float(Tt_override)
    else:
        vals = {v for v in Tts.values() if v is not None}
        if len(vals) != 1 or any(v is None for v in Tts.values()):
            gprob.append(f"入口の Tt が run で一つに決まらない ({Tts}); --Tt で与える")
            Tt = float("nan")
        else:
            Tt = vals.pop()
    res = []
    for k, smp in enumerate(samples):
        if ns is None:
            r = {**smp, "res": f"{smp['run']}/res_{smp['local_step']}.h5", "status": "skipped", "problems": ["評価座標が無い (全体の前提)"]}
        else:
            r = evaluate_sample(case, smp, ns, Tt, species_source)
        res.append(r)
        if r["status"] == "ok":
            print(f"[{k + 1}/{len(samples)}] {_tag(r)}: P 傾き A {r['P_slope']['A']:.6f} B {r['P_slope']['B']:.6f} "
                  f"(A−B {r['P_slope']['A_minus_B']:+.2e}) | T0 最大 A {r['T0']['A_max']:.2f} B {r['T0']['B_max']:.2f} "
                  f"|A−B|max {r['T0']['A_minus_B']['max_abs']:.2e} K | 超過節点 A {r['T0']['A_n_exceed']} B {r['T0']['B_n_exceed']}", flush=True)
        else:
            print(f"[{k + 1}/{len(samples)}] {_tag(r)}: {r['status']} — {'; '.join(r['problems'])}", flush=True)
    verdict = judge(res, gprob, partial, registered_samples())
    deps = {n: _sha(p) for n, p in (("total_quantities.py", TOOLS / "total_quantities.py"), ("forge_species.py", TOOLS / "forge_species.py"),
                                     ("nozzle_report.py", ROOT / "design/forge_design/report/nozzle_report.py"),
                                     ("moc_v5_euler_eval.py", HERE / "moc_v5_euler_eval.py"),
                                     ("eval_wallfit_euler.py", HERE / "eval_wallfit_euler.py"))}
    return {"plan": PLAN, "plan_reg_commit": PLAN_REG_COMMIT, "evaluator": Path(__file__).name,
            "evaluator_sha256": _sha(Path(__file__).resolve()), "dependencies_sha256": deps, "case": str(case),
            "thresholds": {"P_slope_pt": TOL_PSLOPE, "T0_K": TOL_T0, "exceed_margin_K": EXCEED_MARGIN,
                           "excess_hundreds_K": EXCESS_HUNDREDS, "excess_hundreds_note": "「数百 K」は登録に数値が無いので下限 100 K と読む (実装時)"},
            "newton": {"T_init_K": NEWTON_T_INIT, "tol_K": NEWTON_TOL, "max_it": NEWTON_MAX_IT}, "species_source": species_source,
            "Tt": Tt, "Tt_by_run": Tts, "Tt_override": Tt_override, "partial": partial, **ginfo,
            "verdict": verdict, "samples": res}


def main(argv=None):
    ap = argparse.ArgumentParser(description="plan §6 V5c 熱力学の整合の 0 step A/B")
    ap.add_argument("case_dir", nargs="?", default=str(HERE))
    ap.add_argument("--only", default=None, help="RUN:STEP[,RUN:STEP...] (部分標本。登録の判定に使わない)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--Tt", type=float, default=None)
    ap.add_argument("--species-source", default="auto", choices=("auto", "record", "speciesDBFile"))
    a = ap.parse_args(argv)
    case = Path(a.case_dir).resolve()
    partial = a.only is not None
    samples = parse_only(a.only) if partial else registered_samples()
    out = evaluate(case, samples, partial, a.Tt, a.species_source)
    dst = Path(a.out) if a.out else case / (OUT_JSON_PARTIAL if partial else OUT_JSON)
    if not dst.is_absolute():
        dst = case / dst
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(out, indent=1, ensure_ascii=False, allow_nan=True))
    print(f"VERDICT: {out['verdict']['label']}\n  {out['verdict']['text']}\n  -> {dst}")
    return out


if __name__ == "__main__":
    main()
