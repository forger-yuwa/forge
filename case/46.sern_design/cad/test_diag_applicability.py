#!/usr/bin/env python3
"""diag_applicability.py の否定・肯定試験 (plan convection-slau-wall-normal-chi-usage-rule codex result M1/M2)。

**AWS で再照合する前に** commit する (ツールを直してから測る順を痕跡に残す)。合成配列で純関数を叩く。
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from diag_applicability import check_config, check_faces
from diag_wall_cv_budget import slau_mdot

fails = 0
def check(name, cond):
    global fails
    print(("PASS " if cond else "FAIL ") + name); fails += (not cond)

# ---- 合成メッシュ: CV 0 (壁) - CV 1, CV 0 - CV 2 の内部面 2 枚 + CV 0 の境界半割面 1 枚 ----
own = np.array([0, 0, 0]); nei = np.array([1, 2, -1])
S = np.array([[1e-4, 0, 0], [0, 1e-4, 0], [0, 0, 1e-4]])
st = {"ro": np.array([1e-3, 1e-2, 1e-2]), "Ux": np.array([0.0, 1100.0, 300.0]), "Uy": np.array([0.0, 50.0, 0.0]),
      "Uz": np.zeros(3), "P": np.array([60.0, 3500.0, 3000.0]), "sonic": np.array([700.0, 700.0, 700.0])}
def kernel_mf(st):
    mf = np.zeros(3)
    for k in (0, 1):
        A = np.linalg.norm(S[k]); n = S[k] / A
        L = {q: st[q][own[k]] for q in st}; R = {q: st[q][nei[k]] for q in st}
        mf[k] = slau_mdot(A, *n, L, R)[0]
    return mf
wall = {0}; bk = {2: "wall:wall"}
mf = kernel_mf(st)
ok, bad, per, sm = check_faces(own, nei, S, mf, st, [0], wall, bk)
check("肯定: 一致する流束 + 境界面 0 → 診断可能", ok and sm["n_int"] == 2 and sm["n_bnd"] == 1)
check("肯定: 接続面は内部 2 + 境界 1 = 3 面すべて数える", per[0]["faces"] == 3)

st_nan = {k: v.copy() for k, v in st.items()}; st_nan["P"][1] = np.nan
ok, bad, _, sm = check_faces(own, nei, S, mf, st_nan, [0], wall, bk)
check("否定: 状態に NaN → 診断不能", (not ok) and sm["n_nonfinite"] >= 1)
mf_nan = mf.copy(); mf_nan[0] = np.nan
ok, _, _, _ = check_faces(own, nei, S, mf_nan, st, [0], wall, bk)
check("否定: カーネル流束に NaN → 診断不能", not ok)
ok, bad, _, _ = check_faces(own, nei, S, mf, st, [5], wall, bk)
check("否定: CV 範囲外 → 診断不能", not ok)
ok, bad, _, _ = check_faces(own, nei, S, mf, st, [], wall, bk)
check("否定: 対象 CV 空集合 (照合 0 面) → 診断不能", (not ok) and any("0 面" in b for b in bad))
ok, bad, _, _ = check_faces(own, nei, S, mf[:2], st, [0], wall, bk)
check("否定: massflux の面数不一致 (欠落面) → 診断不能", not ok)
st_neg = {k: v.copy() for k, v in st.items()}; st_neg["ro"][1] = -1e-3
ok, _, _, _ = check_faces(own, nei, S, kernel_mf(st), st_neg, [0], wall, bk)
check("否定: 負の密度 → 診断不能", not ok)
mf_b = mf.copy(); mf_b[2] = 1e-3
ok, _, _, _ = check_faces(own, nei, S, mf_b, st, [0], wall, bk)
check("否定: 境界半割面に有意な流束 → 診断不能", not ok)
mf_off = mf.copy(); mf_off[0] *= 1.5
ok, _, _, _ = check_faces(own, nei, S, mf_off, st, [0], wall, bk)
check("否定: 内部面がツールと食い違う → 診断不能", not ok)

# ---- キー検査 ----
cfg = {"mesh": {"discretization": "node", "nodeWallDirichlet": 1}, "solver": "SLAU",
       "space": {"convMethod": 0}, "physProp": {"thermalMethod": 2}}
echo = {"solver": "SLAU", "space.convMethod": "0", "physProp.thermalMethod": "2", "space.slauWallNormalChi": "0"}
_, bad, unconf = check_config(cfg, echo, {"wall": {"physID": 3, "kind": "wall"}}, True)
check("肯定: 既知構成の設定 → キー OK", bad == [])
check("肯定: エコーの無いキーは『未確認』に列挙", "space.slauContactFloor" in unconf)
_, bad, _ = check_config(cfg, echo, {}, False)
check("否定: ログ欠落 → 診断不能", bool(bad))
_, bad, _ = check_config(cfg, dict(echo, **{"physProp.thermalMethod": "0"}), {}, True)
check("否定: thermalMethod の yaml/エコー不一致 → 診断不能", bool(bad))
c2 = {**cfg, "physProp": {"thermalMethod": 2, "isAxisymmetric": 1}}
_, bad, _ = check_config(c2, echo, {}, True)
check("否定: 旧配置 physProp.isAxisymmetric: 1 → 診断不能", bool(bad))
_, bad, _ = check_config(cfg, echo, {"side": {"physID": 7, "kind": "periodic"}}, True)
check("否定: 周期 bcond → 診断不能", bool(bad))
c3 = {**cfg, "space": {"convMethod": 0, "slauContactFloor": 0.01}}
_, bad, _ = check_config(c3, echo, {}, True)
check("否定: slauContactFloor 0.01 → 診断不能", bool(bad))
_, bad, _ = check_config(cfg, dict(echo, **{"space.convMethod": "1"}), {}, True)
check("否定: convMethod の yaml 0 / エコー 1 → 診断不能", bool(bad))

from diag_applicability import log_echo
e = log_echo("'slauWallNormalChi' effective: 1 (auto: node+nodeWallDirichlet+SLAU)  (wall-adjacent ...)\n")
check("新エコー形式を読める (実効 1)", e.get("space.slauWallNormalChi") == "1")
_, bad, _ = check_config(cfg, dict(echo, **{"space.slauWallNormalChi": "1"}), {}, True)
check("否定: 実効 1 (auto で有効) の run → 診断不能 (診断は flag 0 の run で行う)", bool(bad))
_, bad, _ = check_config(cfg, {k: v for k, v in echo.items() if k != "space.slauWallNormalChi"}, {}, True)
check("否定: 省略かつエコー無し → 実効値未確定で診断不能", bool(bad))
_, bad, _ = check_config(cfg, dict(echo, **{"space.slauWallNormalChi": "0"}), {}, True)
check("肯定: エコーで実効 0 → キー OK", bad == [])
print("VERDICT:", "PASS" if fails == 0 else f"FAIL ({fails})")
sys.exit(1 if fails else 0)
