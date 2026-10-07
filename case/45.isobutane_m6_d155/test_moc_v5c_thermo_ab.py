"""moc_v5c_thermo_ab.py の試験 (plan discretization-moc-axis-limit-and-corrector §6 V5c)。forge・AWS・V5/V5b の run の結果は使わない。
一時ディレクトリの模擬 case (TP 2 種 MIX + H2O、speciesDBFile、40 × 7 節点の構造格子) に、合成の res_*.h5 を書いて回す。
保存量 (ro・roU*・roe・roY*) と原始量 (P・T・Y*・h0) は同じ真の状態から float32 で作る (一致の場合)。
  標本: 登録の 56 標本 (窓 A・B の各 13 枚と親の初期場・最終場、通算 step)。
  一致: 両経路が閾値内で、壁際の超過が両方にある → 退ける / 超過が無い → どちらでもない / B が真の T・P・h0 を復元する。
  一致しない: 保存 T だけ +2 K (登録の 2 量には出ない、記録の T の差に出る) / 原始量を T + 2 K の状態から作る / 保存 h0 + 2 cp → 棄却。
  超過が片方だけ: 保存 h0 だけ壁際 +300 cp (A だけ) / roe だけ壁際 +300 cv (B だけ) → 棄却、事実に「A にだけ」「B にだけ」。
  独立性: 保存 P・T・Y*・h0 の値を NaN にしても経路 B の結果がビット同一。h0_includes_k = 1 なら B の h0 は k を含む (属性に従う)。
  判定不能: 標本の欠損・保存量の非有限値・時系列の CSV との不一致・評価座標の記録の欠損・評価器の sha256 の違い・Tt が決まらない。
  P 傾きの定義: 経路 A の値が eval_wallfit_euler.quantities と同じ (差 0)、CSV が記録と違えば照合に使わない。出力 JSON の記録。
usage: python3 test_moc_v5c_thermo_ab.py
"""
import hashlib
import json
import shutil
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

C = Path(__file__).resolve().parent
sys.path.insert(0, str(C))
import moc_v5c_thermo_ab as V  # noqa: E402

FSP = V.FSP
fails = 0


def check(name, cond):
    global fails
    print(("ok   " if cond else "FAIL ") + name)
    fails += 0 if cond else 1


def sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


NI, NJ, SCALE, TT = 40, 7, 0.05, 1600.0
GEOM = {"X_E": 3.0, "X_F": 10.0, "WIN_T": [5.0, 9.0], "WIN_O": [-12.0, 10.0]}
CHILD, PARENT = "run_0157_euler_wallfit_mocG1_r1_ext36k", "run_0150_euler_wallfit_mocG1_r1"
SPECIES_DB = {   # MIX は case/44 の species_db.yaml、H2O は解決済み記録 (CEA) の係数
    "MIX": {"MW": 0.029080607436051467, "Tlo": 200.0, "Tmid": 1000.0, "Thi": 6000.0,
            "nasa9_low": [8300.191110009942, -162.22037972364785, 4.747084537742655, -0.004636579666164405, 9.186775751715496e-06,
                          -6.950235465013877e-09, 1.927164588238703e-12, -3436.5657689852787, -3.079713666222943],
            "nasa9_high": [245924.4954317055, -1252.909301603011, 5.166006679945892, -6.12911885105234e-05, 2.704127627261617e-08,
                           -4.853694170836424e-12, 3.7129709260291397e-16, 3330.9623564588733, -8.774252243539319]},
    "H2O": {"MW": 0.01801528, "Tlo": 200.0, "Tmid": 1000.0, "Thi": 6000.0,
            "nasa9_low": [-39479.6083, 575.573102, 0.931782653, 0.00722271286, -7.34255737e-06, 4.95504349e-09, -1.336933246e-12,
                          -33039.7431, 17.24205775],
            "nasa9_high": [1034972.096, -2412.698562, 4.64611078, 0.002291998307, -6.83683048e-07, 9.42646893e-11, -4.82238053e-15,
                           -13842.86509, -7.97814851]},
}


def mesh_xr():
    x = np.linspace(0.0, 10.0, NI); eta = np.linspace(0.0, 1.0, NJ)
    X = np.repeat(x[:, None], NJ, axis=1)
    R = eta[None, :] * (1.0 + 0.1 * x[:, None])
    return X, R


def wall_mask(X):
    """壁際 (j_w ≤ 1) かつ 2 ≤ x ≤ 3 の節点 (平坦化した bool)。"""
    m = np.zeros(X.shape, bool)
    m[(X[:, 0] >= 2.0) & (X[:, 0] <= 3.0), NJ - 2:] = True
    return m.ravel()


def mock_run(case: Path, run: str, Tt=TT):
    d = case / run
    d.mkdir(parents=True, exist_ok=True)
    (d / "solverConfig.yaml").write_text(
        "physProp: {thermalMethod: 2, cp: 1360.0, gamma: 1.27, species: [MIX, H2O], speciesDBFile: species_db.yaml, thermoHrefTemp: 298.15}\n")
    import yaml
    (d / "species_db.yaml").write_text(yaml.safe_dump(SPECIES_DB))
    (d / "bcondConfig.yaml").write_text(
        f"inlet: {{physID: 1, kind: inlet_Pressure, floats: {{Y0: 0.9142, Y1: 0.0858, Pt: 5500000.0, Tt: {Tt}}}}}\n"
        "outlet: {physID: 2, kind: outlet_statPress, floats: {Ps: 2237.0, Pt: 2237.0, Tt: 300.0}}\n")
    (d / "prepare_info.json").write_text(json.dumps({"scale_m": SCALE, "mesh": {"ni": NI}}))
    X, R = mesh_xr()
    coord = np.stack([X.ravel() * SCALE, R.ravel() * SCALE, np.zeros(X.size)], axis=1).astype(np.float32)
    with h5py.File(d / "nozzle.h5", "w") as f:
        f.create_dataset("/MESH/COORD", data=coord.ravel())
    return d


def gas_of(d: Path):
    return FSP.thermo_gas(FSP.run_thermo(str(d)))


def true_state(excess=False, T_fn=None):
    X, R = mesh_xr()
    x, eta = X.ravel(), (R / R[:, -1:].clip(min=1e-30)).ravel()
    T = 1500.0 - 20.0 * x + 30.0 * eta if T_fn is None else T_fn(x, eta)
    if excess:
        T = np.where(wall_mask(X), T + 400.0, T)
    P = 2237.0 * (1.0 + 0.02 * x + 0.003 * eta)
    yw = 0.0858 + 0.002 * eta
    return {"T": T, "P": P, "Y": [1.0 - yw, yw], "Ux": 400.0 + 30.0 * x, "Uy": 20.0 * eta, "Uz": np.zeros_like(x)}


def write_res(path: Path, gas, st: dict, inc_k=0, k=0.0, prim_dT=0.0, roe_add=None, mutate=None, h0_with_k=None):
    """真の状態 st から保存量と原始量を作って float32 で書く。prim_dT: 原始量 (P・T・h0・sonic) だけを T + prim_dT の状態から作る。
    roe_add: 保存量 roe への加算。mutate {名前: 関数}: 書く直前の配列の書き換え。h0_with_k: 保存 h0 に k を入れるか (既定 inc_k)。"""
    Y, T = st["Y"], st["T"]
    Rm = gas.Rmix(Y)
    ro = st["P"] / (Rm * T)
    ek = 0.5 * (st["Ux"] ** 2 + st["Uy"] ** 2 + st["Uz"] ** 2)
    kk = np.full_like(T, k)
    roe = ro * (gas.h(Y, T) - Rm * T + ek)
    Tp = T + prim_dT
    cp = gas.cp(Y, Tp)
    vals = {"ro": ro, "roUx": ro * st["Ux"], "roUy": ro * st["Uy"], "roUz": ro * st["Uz"], "roe": roe if roe_add is None else roe + roe_add,
            "roK": ro * kk, "k": kk, "P": ro * Rm * Tp, "T": Tp, "Ux": st["Ux"], "Uy": st["Uy"], "Uz": st["Uz"],
            "h0": gas.h(Y, Tp) + ek + (kk if (inc_k if h0_with_k is None else h0_with_k) else 0.0),
            "sonic": np.sqrt(cp / (cp - Rm) * Rm * Tp)}
    for i, y in enumerate(Y):
        vals[f"roY{i}"] = ro * y
        vals[f"Y{i}"] = y
    for nm, fn in (mutate or {}).items():
        vals[nm] = fn(vals[nm])
    with h5py.File(path, "w") as f:
        for nm, v in vals.items():
            ds = f.create_dataset(f"/VALUE/{nm}", data=np.asarray(v, dtype=np.float32))
            if nm == "h0":
                ds.attrs["h0_includes_k"] = int(inc_k)
    return vals


def mock_case(root: Path, ev_sha=None, geom=GEOM, write_geom=True, Tt_parent=TT) -> Path:
    case = root / "case"
    (case / "_band_ab").mkdir(parents=True)
    if write_geom:
        rec = {"geometry": geom, "evaluator_sha256": ev_sha or sha(C / "eval_wallfit_euler.py"), "runs": {}}
        (case / "_band_ab/wallfit_series_v5b.json").write_text(json.dumps(rec))
        (case / "_band_ab/wallfit_series_v5.json").write_text(json.dumps(rec))
    mock_run(case, CHILD)
    mock_run(case, PARENT, Tt=Tt_parent)
    return case


SAMPLES = [V.sample_of(CHILD, 36000), V.sample_of(PARENT, 18000)]


def run_case(st_child=None, st_parent=None, **kw):
    """模擬 case を作り、子 res_36000 と親 res_18000 を書いて部分標本で評価する。kw は write_res へ (両方に同じ変更)。"""
    tmp = Path(tempfile.mkdtemp(prefix="v5c_"))
    try:
        case = mock_case(tmp)
        g = gas_of(case / CHILD)
        write_res(case / CHILD / "res_36000.h5", g, st_child or true_state(excess=True), **kw)
        write_res(case / PARENT / "res_18000.h5", g, st_parent or true_state(excess=True), **kw)
        return V.evaluate(case, SAMPLES, partial=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def label(out):
    return out["verdict"]["label"].replace(V.PARTIAL_PREFIX, "")


# --- 標本 -------------------------------------------------------------------------------------------------------------
reg = V.registered_samples()
check("登録の標本は 56", len(reg) == 56)
for child, parent in V.CHILDREN.items():
    a = [s["total_step"] for s in reg if s["run"] == child and s["role"] == "child_winA"]
    b = [s["total_step"] for s in reg if s["run"] == child and s["role"] == "child_winB"]
    check(f"{child[:8]} 窓 A は通算 24000〜36000 の 13 枚", a == list(range(24000, 36001, 1000)))
    check(f"{child[:8]} 窓 B は通算 42000〜54000 の 13 枚", b == list(range(42000, 54001, 1000)))
    check(f"{child[:8]} 窓 A の local は 6000〜18000",
          [s["local_step"] for s in reg if s["run"] == child and s["role"] == "child_winA"] == list(range(6000, 18001, 1000)))
    check(f"{parent[:8]} は初期場 0 と最終場 18000", sorted(s["local_step"] for s in reg if s["run"] == parent) == [0, 18000])
check("登録の run は run_0157・run_0160・run_0150・run_0153",
      {s["run"] for s in reg} == {"run_0157_euler_wallfit_mocG1_r1_ext36k", "run_0160_euler_icdep_mocG1_isen_ext36k",
                                  "run_0150_euler_wallfit_mocG1_r1", "run_0153_euler_icdep_mocG1_isen"})
check("--only の解釈 (子の local 36000 は窓 B・通算 54000)", V.parse_only(f"{CHILD}:36000")[0]["total_step"] == 54000
      and V.parse_only(f"{CHILD}:36000")[0]["role"] == "child_winB")

# --- 一致 (超過あり) → 退ける -----------------------------------------------------------------------------------------
out = run_case()
s0 = out["samples"][0]
check("一致・超過あり: 退ける", label(out) == V.LBL_REJECT_PATH)
check("一致・超過あり: 部分標本と明記", out["verdict"]["label"].startswith(V.PARTIAL_PREFIX) and out["partial"])
check("一致・超過あり: 全温の |A − B| ≤ 0.01 K", out["verdict"]["facts"]["max_abs_T0_A_minus_B_K"] <= 0.01)
check("一致・超過あり: P 傾きの |A − B| ≤ 1e-5 %pt", out["verdict"]["facts"]["max_abs_P_slope_A_minus_B"] <= 1e-5)
check("一致・超過あり: 2 標本とも両経路で超過", out["verdict"]["facts"]["n_excess_both"] == 2)
check("P 傾きの定義: A は quantities と同じ (差 0)", s0["P_slope"]["A_identity_diff"] == 0.0)
check("P 傾き: 合成の場で有限・非零", np.isfinite(s0["P_slope"]["A"]) and abs(s0["P_slope"]["A"]) > 0.1)
check("熱物性の読み元が A と B で同じ", s0["thermo_source_A"] == s0["thermo_source_B"] and "speciesDBFile" in s0["thermo_source_A"])
check("B の座標が load_field と同じ", s0["geometry_B_equals_A"])
check("B の Newton が全節点で収束", s0["B_diag"]["newton_T_n_nonconverged"] == 0 and s0["B_diag"]["newton_T0_n_nonconverged"] == 0)
check("層: 壁 (j_w 0・1) に超過、j_w ≥ 2 には無い (A・B)",
      all(s0["layers"][p]["n_exceed"][0] > 0 and s0["layers"][p]["n_exceed"][1] > 0 and sum(s0["layers"][p]["n_exceed"][2:]) == 0
          for p in ("A", "B")))
check("層: 層の数は nj", len(s0["layers"]["A"]["max"]) == NJ and s0["layers"]["j_from_wall"] == list(range(NJ)))
check("最大の位置: 壁際 (j_w ≤ 1) で 2 ≤ x ≤ 3", s0["T0"]["argmax_A"]["j_from_wall"] <= 1 and 2.0 <= s0["T0"]["argmax_A"]["x_rt"] <= 3.0
      and s0["T0"]["argmax_B"]["j_from_wall"] <= 1)
check("超過節点の数: A と B が同じ (壁際 2 層 × 該当列)", s0["T0"]["A_n_exceed"] == s0["T0"]["B_n_exceed"] == int(wall_mask(mesh_xr()[0]).sum()))
check("h0 の整合: 保存 h0 − (roe + 保存 P)/ro は float32 の丸め程度 (≤ 0.01 K)", s0["h0"]["stored_minus_roeP_stored_K"]["max_abs"] <= 0.01)
check("判定文に事実 (超過の内訳) がある", "両経路 2・A だけ 0・B だけ 0" in out["verdict"]["text"])

# B が真の状態を復元する (低温・高温の範囲を含む)
tmp = Path(tempfile.mkdtemp(prefix="v5c_"))
case = mock_case(tmp)
g = gas_of(case / CHILD)
stw = true_state(T_fn=lambda x, eta: 220.0 + 280.0 * x + 50.0 * eta)      # 220〜3070 K (NASA-9 の区切り 1000 K をまたぐ)
write_res(case / CHILD / "res_36000.h5", g, stw)
B = V.path_b(case / CHILD, case / CHILD / "res_36000.h5")
Rm = g.Rmix(stw["Y"]); ek = 0.5 * (stw["Ux"] ** 2 + stw["Uy"] ** 2)
check("B: T を真値から 0.01 K 以内に復元 (220〜3070 K)", float(np.max(np.abs(B["T"] - stw["T"]))) < 0.01)
check("B: P を相対 1e-5 以内に復元", float(np.max(np.abs(B["P"] / stw["P"] - 1))) < 1e-5)
check("B: h0 を真値から 0.01 cp 以内に復元", float(np.max(np.abs(B["h0"] - (g.h(stw["Y"], stw["T"]) + ek)) / g.cp(stw["Y"], stw["T"]))) < 0.01)
check("B: Newton が収束 (初期値 1000 K から)", B["diag"]["newton_T_n_nonconverged"] == 0 and B["diag"]["newton_T0_n_nonconverged"] == 0)
check("B: 読むものの一覧に P・T・h0 の値が無い", not any(r in ("P", "T", "h0") or r.startswith("Y") for r in B["reads"]))

# 独立性: 保存 P・T・Y*・h0 の値を NaN にしても B はビット同一
p2 = case / CHILD / "res_0.h5"
write_res(p2, g, stw, mutate={"P": lambda a: a * np.nan, "T": lambda a: a * np.nan, "h0": lambda a: a * np.nan,
                              "Y0": lambda a: a * np.nan, "Y1": lambda a: a * np.nan})
B2 = V.path_b(case / CHILD, p2)
check("独立性: 保存 P・T・Y*・h0 を NaN にしても B の T・P・h0・T0 がビット同一",
      all(np.array_equal(B[k], B2[k]) for k in ("T", "P", "h0", "T0")))
with h5py.File(p2, "r+") as f:
    del f["VALUE/roe"]
try:
    V.path_b(case / CHILD, p2); ok = False
except KeyError:
    ok = True
check("B: 保存量 roe が無ければ例外 (既定値で埋めない)", ok)
shutil.rmtree(tmp, ignore_errors=True)

# --- 一致 (超過なし) → どちらでもない --------------------------------------------------------------------------------------
out = run_case(st_child=true_state(), st_parent=true_state())
check("一致・超過なし: どちらでもない", label(out) == V.LBL_NEITHER)
check("一致・超過なし: 判定文に「超過が無い」", "超過が無い" in out["verdict"]["text"])

# --- 一致しない ----------------------------------------------------------------------------------------------------------
out = run_case(mutate={"T": lambda a: a + 2.0})
s0 = out["samples"][0]
check("保存 T だけ +2 K: 登録の 2 量は閾値内 (全温は h0 から、P 傾きは P から) → 退けるのまま", label(out) == V.LBL_REJECT_PATH)
check("保存 T だけ +2 K: 記録の T の差は +2 K", abs(s0["T"]["stored_minus_B_K"]["min"] - 2.0) < 1e-3 and abs(s0["T"]["stored_minus_B_K"]["max"] - 2.0) < 1e-3)

out = run_case(prim_dT=2.0)
s0 = out["samples"][0]
check("原始量を T + 2 K の状態から作る: 棄却", label(out) == V.LBL_REJECT_CONS)
check("原始量を T + 2 K: 全温の A − B は約 +2 K", 1.5 < s0["T0"]["A_minus_B"]["max"] < 2.5 and s0["T0"]["A_minus_B"]["min"] > 1.5)
check("原始量を T + 2 K: P の相対差が記録される (> 1e-3)", s0["P"]["stored_minus_B_rel"]["max_abs"] > 1e-3)
check("原始量を T + 2 K: h0 の整合 (保存 h0 − B の h0) も約 2 K", 1.5 < s0["h0"]["stored_minus_B_K"]["max"] < 2.5)

out = run_case(mutate={"h0": lambda a: a + 2.0 * 1360.0})
s0 = out["samples"][0]
check("保存 h0 + 2 cp: 棄却", label(out) == V.LBL_REJECT_CONS)
check("保存 h0 + 2 cp: 全温だけが食い違い、P 傾きは閾値内", out["verdict"]["facts"]["n_disagree_T0"] == 2 and out["verdict"]["facts"]["n_disagree_P_slope"] == 0)
check("保存 h0 + 2 cp: 判定文に IC 依存と解釈しない旨", "IC 依存と解釈しない" in out["verdict"]["text"])

# --- 超過が片方だけ ------------------------------------------------------------------------------------------------------
Xm = mesh_xr()[0]
wm = wall_mask(Xm)
out = run_case(st_child=true_state(), st_parent=true_state(), mutate={"h0": lambda a: np.where(wm, a + 300.0 * 1360.0, a)})
f_ = out["verdict"]["facts"]
check("超過が A だけ: 棄却", label(out) == V.LBL_REJECT_CONS)
check("超過が A だけ: 事実 (A だけ 2・B だけ 0)", f_["n_excess_A_only"] == 2 and f_["n_excess_B_only"] == 0)
check("超過が A だけ: 判定文に「A にだけ」", "A にだけ" in out["verdict"]["text"])

st = true_state()
tmp = Path(tempfile.mkdtemp(prefix="v5c_"))
case = mock_case(tmp)
g = gas_of(case / CHILD)
Rm = g.Rmix(st["Y"]); ro = st["P"] / (Rm * st["T"]); cv = g.cp(st["Y"], st["T"]) - Rm
for run, step in ((CHILD, 36000), (PARENT, 18000)):
    write_res(case / run / f"res_{step}.h5", g, st, roe_add=np.where(wm, ro * cv * 300.0, 0.0))
out = V.evaluate(case, SAMPLES, partial=True)
f_ = out["verdict"]["facts"]
check("超過が B だけ: 棄却", label(out) == V.LBL_REJECT_CONS)
check("超過が B だけ: 事実 (A だけ 0・B だけ 2)", f_["n_excess_A_only"] == 0 and f_["n_excess_B_only"] == 2)
check("超過が B だけ: 判定文に「B にだけ」", "B にだけ" in out["verdict"]["text"])
check("超過が B だけ: B の超過は約 300 K", 250.0 < f_["max_excess_B_K"] < 450.0 and f_["max_excess_A_K"] < 100.0)
shutil.rmtree(tmp, ignore_errors=True)

# --- k の扱い (h0_includes_k に従う) --------------------------------------------------------------------------------------
out = run_case(inc_k=1, k=5000.0)
s0 = out["samples"][0]
check("h0_includes_k = 1 (保存 h0 も k 込み): 一致して退ける", label(out) == V.LBL_REJECT_PATH and s0["h0_includes_k"] == 1)
check("h0_includes_k = 1: B は roK を読む", "roK" in s0["b_reads"])
out = run_case(inc_k=1, k=5000.0, h0_with_k=False)
check("属性 = 1 なのに保存 h0 が k を含まない: 全温が食い違って棄却", label(out) == V.LBL_REJECT_CONS)
out = run_case(inc_k=0, k=5000.0)
check("h0_includes_k = 0 (roK があっても): B は k を足さず一致", label(out) == V.LBL_REJECT_PATH and "roK" not in out["samples"][0]["b_reads"])

# --- 判定不能 ------------------------------------------------------------------------------------------------------------
tmp = Path(tempfile.mkdtemp(prefix="v5c_"))
case = mock_case(tmp)
g = gas_of(case / CHILD)
write_res(case / CHILD / "res_36000.h5", g, true_state(excess=True))
write_res(case / PARENT / "res_18000.h5", g, true_state(excess=True))
out = V.evaluate(case, V.registered_samples(), partial=False)
f_ = out["verdict"]["facts"]
check("登録の全標本で 2 枚だけ: 判定不能", out["verdict"]["label"] == V.LBL_UNDET)
check("登録の全標本で 2 枚だけ: 欠損 54 は missing として数える", f_["n_ok"] == 2 and len(f_["problem_samples"]) == 54)
check("登録の全標本で 2 枚だけ: 部分標本の印は付かない", not out["partial"] and not out["verdict"]["label"].startswith(V.PARTIAL_PREFIX))
# 子の res を消すと --only でも判定不能
(case / CHILD / "res_36000.h5").unlink()
out = V.evaluate(case, SAMPLES, partial=True)
check("--only で res が無い: 判定不能 (部分標本)", label(out) == V.LBL_UNDET and out["samples"][0]["status"] == "missing")
# 保存量の非有限値
write_res(case / CHILD / "res_36000.h5", g, true_state(excess=True), mutate={"roe": lambda a: np.where(np.arange(a.size) == 5, np.nan, a)})
out = V.evaluate(case, SAMPLES, partial=True)
check("保存量 roe に NaN: 判定不能", label(out) == V.LBL_UNDET and out["samples"][0]["status"] == "problem")
check("保存量 roe に NaN: 理由に非有限値", any("非有限" in p for p in out["samples"][0]["problems"]))

# 時系列の CSV との照合
write_res(case / CHILD / "res_36000.h5", g, true_state(excess=True))
out = V.evaluate(case, SAMPLES, partial=True)
sA = out["samples"][0]["P_slope"]["A"]
check("CSV が無い: 照合なしで通る", out["samples"][0]["P_slope"]["series"]["status"].startswith("CSV が無い") and label(out) == V.LBL_REJECT_PATH)
cols = "step,P_slope_eta0.1,exit_core_M\n"
csvp = case / CHILD / "wallfit_series_v5b.csv"


def put_csv(val, register=True):
    csvp.write_text(cols + f"35000,0.1,6.0\n36000,{val:.10g},6.0\n")
    rp = case / "_band_ab/wallfit_series_v5b.json"
    rec = json.loads(rp.read_text())
    rec["runs"] = {CHILD: {"csv_sha256": sha(csvp) if register else "0" * 64}}
    rp.write_text(json.dumps(rec))


put_csv(sA)
out = V.evaluate(case, SAMPLES, partial=True)
check("CSV (記録と同じ) の値と一致: 照合 ok", out["samples"][0]["P_slope"]["series"]["status"] == "ok"
      and abs(out["samples"][0]["P_slope"]["A_minus_series"]) <= V.SERIES_TOL and label(out) == V.LBL_REJECT_PATH)
put_csv(sA + 0.01)
out = V.evaluate(case, SAMPLES, partial=True)
check("CSV (記録と同じ) の値と違う: 判定不能 (定義・座標の再現が不成立)", label(out) == V.LBL_UNDET
      and any("時系列の記録" in p for p in out["samples"][0]["problems"]))
put_csv(sA + 0.01, register=False)
out = V.evaluate(case, SAMPLES, partial=True)
check("CSV の sha256 が記録と違う: 照合に使わない (判定は通常どおり)", out["samples"][0]["P_slope"]["series"]["status"].startswith("記録と違う")
      and label(out) == V.LBL_REJECT_PATH)
csvp.unlink()
shutil.rmtree(tmp, ignore_errors=True)

# 評価座標の記録・評価器の sha256・Tt
for name, kw in (("評価座標の記録が無い", {"write_geom": False}), ("評価器の sha256 が V5b の記録と違う", {"ev_sha": "f" * 64}),
                 ("評価座標が不正 (X_E ≥ X_F)", {"geom": {**GEOM, "X_E": 11.0}}), ("Tt が run で違う", {"Tt_parent": 1700.0})):
    tmp = Path(tempfile.mkdtemp(prefix="v5c_"))
    case = mock_case(tmp, **kw)
    g = gas_of(case / CHILD)
    write_res(case / CHILD / "res_36000.h5", g, true_state(excess=True))
    write_res(case / PARENT / "res_18000.h5", g, true_state(excess=True))
    out = V.evaluate(case, SAMPLES, partial=True)
    check(f"{name}: 判定不能", label(out) == V.LBL_UNDET and out["verdict"]["facts"]["global_problems"])
    if name == "Tt が run で違う":
        out = V.evaluate(case, SAMPLES, partial=True, Tt_override=1600.0)
        check("Tt が run で違う: --Tt を与えれば判定できる", label(out) == V.LBL_REJECT_PATH and out["Tt"] == 1600.0)
    shutil.rmtree(tmp, ignore_errors=True)

# --- 出力 JSON (main) ----------------------------------------------------------------------------------------------------
tmp = Path(tempfile.mkdtemp(prefix="v5c_"))
case = mock_case(tmp)
g = gas_of(case / CHILD)
write_res(case / CHILD / "res_36000.h5", g, true_state(excess=True))
V.main([str(case), "--only", f"{CHILD}:36000"])
jp = case / V.OUT_JSON_PARTIAL
d = json.loads(jp.read_text()) if jp.is_file() else {}
check("main --only: 部分標本の既定の出力先に書く (登録の出力先には書かない)", jp.is_file() and not (case / V.OUT_JSON).exists())
check("出力 JSON: 登録の commit 4a484270", d.get("plan_reg_commit") == "4a484270")
check("出力 JSON: 評価器自身の sha256", d.get("evaluator_sha256") == sha(C / "moc_v5c_thermo_ab.py"))
check("出力 JSON: 依存の sha256 (total_quantities・forge_species・eval_wallfit_euler)",
      all(d.get("dependencies_sha256", {}).get(k) for k in ("total_quantities.py", "forge_species.py", "eval_wallfit_euler.py")))
check("出力 JSON: 閾値 (0.003 %pt・1 K) と「数百 K」の読み", d.get("thresholds", {}).get("P_slope_pt") == 0.003
      and d.get("thresholds", {}).get("T0_K") == 1.0 and d.get("thresholds", {}).get("excess_hundreds_K") == 100.0)
check("出力 JSON: 評価座標は V5b の記録のもの", d.get("geometry") == GEOM)
shutil.rmtree(tmp, ignore_errors=True)

print(f"FAIL 件数: {fails}")
sys.exit(1 if fails else 0)
