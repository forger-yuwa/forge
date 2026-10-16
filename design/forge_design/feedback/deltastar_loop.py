r"""排除厚さの固定点反復ドライバ (固定 Euler 基準・コア整合抽出)。

計画: plans/active/tooling-nozzle-deltastar-core-matched-euler.md §4.5–4.6。

1 pass = 前 pass の NS run から $\delta_r(x)$ を抽出 (`metrics.deltastar.deltastar_from_core_matched_euler`)
→ 緩和して次の入力 $\delta_{in}^{k+1} = (1-\omega)\delta_{in}^k + \omega\,\delta_{use}^k$ (hard 不合格の断面は前回値保持)
→ 半径方向オフセットの物理壁で `prepare_ns` → 段階起動 NS → `collect` (質量流量帳簿込み)
→ 新 run からも抽出して固定点差を記録 (`fixed_point.json`)。

使い方 (リポジトリルートで):
  design/.venv-opt/bin/python -m forge_design.feedback.deltastar_loop \
      --problem case/45.isobutane_m6_d155/problem_d155_ns.yaml \
      --euler-ref case/45.isobutane_m6_d155/run_0001_euler_shortest_dry \
      --prev case/45.isobutane_m6_d155/run_0004_ns_v3 \
      --run-dir case/45.isobutane_m6_d155/run_0019_ns_cm_pass1 [--omega 0.5] [--prepare-only]

補正の表で NS を立てる (2026-10-13, plan tooling-nozzle-deltastar-band-root-zone §4.3 3; 表は `wt_card --ns-correction` が問題 YAML の横に
書く `<stem>_delta_r.csv`。YAML の `_card.ns_correction` があれば表の sha256 を照合する):
  design/.venv-opt/bin/python -m forge_design.feedback.deltastar_loop --problem P.yaml --euler-ref EULER_RUN \
      --delta-r-csv P_delta_r.csv --prev PREV_RUN --run-dir RUN [--stages none --cfl 1 --steps 60000] [--prepare-only]

受け渡しの関門 (2026-10-14, plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #36a、codex diagnose 2026-10-14): `run_pass`・`run_from_table`・
`run_pass0_integral` は NS の本段の後 (`run_pass` は前の pass の run から次の壁を作る前にも) `handoff_gate` を通し、
(a) 本段の停止の状態 (`STOP_RULE_STATE.json`; run_staged_ns の終了コードでは決めない) が REACHED、(b) `check_convergence.py --segment`
(本段の区間) が PASS のときだけ抽出 (`band_select="edge"` を明示) と次の壁の入力 (`delta_r_next.csv`・`fixed_point.json`) へ進む。
判定は `HANDOFF_GATE.json`・`HANDOFF_CONVERGENCE.txt` に残す。明示の上書き (理由の文字列が必須、HANDOFF_GATE.json の override に残す;
§5.1 #36h): `--allow-not-reached "<理由>"` (停止規則を使わない run・CENSORED・UNAVAILABLE から進む)・`--allow-unconverged "<理由>"` (収束の
判定が NOT CONVERGED でも進む)・`--allow-unsettled "<理由>"` (保存則・δ の 2 窓の実測が閾値を超えても進む; 冷却壁の case/70 のように
落ち着かないのが普通の流れ用)。DIVERGED・EXEC_ERROR・DATA_ERROR・状態の記録なし・収束の判定の DIVERGED・判定不能、保存則・2 窓の判定
不能 (記録の欠け・窓を覆わない・点数の不足・非有限) は上書きでも進まない。落ち着きの条件 (停止の状態の REACHED・保存則・2 窓) を上書きで
通した受け渡しは「落ち着いていない場からの設計」として `settled: false`・`unsettled_overrides` を HANDOFF_GATE.json と抽出の要約
(`fixed_point.json`・`delta_r_extract_summary.json`) に残す。止めたら CLI は終了コード 3。**壁の更新を次の pass へ自動でつながない** (1 pass ずつ
呼び出し側が確かめる)。2026-10-14 (§5.1 #36e、codex diagnose 3・§6 V-h370) から関門は §4.9 の「更新に進む条件」の残りも見る
(保存則と 2 窓は実測の超過だけ #36h の `--allow-unsettled` で通せる、ほかは上書きでは通さない): 保存則 (stop_rule.json の r の窓)・
局所 y₁⁺ (`check_wall_resolution.py`)・δ_E と δ_target の隣り合う 2 窓
(見張りの中の出力ごとの記録)・Euler の参照と ω が本段の記録と同じ、それを満たしたら固定の条件 (`handoff_extract_settings`: edge・一定 x の
断面への補間・端は欠測・ω 0.5・λ 1.0・昇格の上限 ×1000) で抽出して全域の被覆・使用域の held と欠測・次の壁の幾何の健全性。全部満たしたときだけ
抽出の表を run の直下に写す。Euler の参照は `--euler-ref-res res_<N>.h5 --euler-ref-sha256 HEX` でスナップショットを固定できる (sha が
合わなければ例外; 無ければ今までどおり最新の出力)。

**関門の新方式** (2026-10-16、plan tooling-settle-convergence-and-monitor §4.9 の 1・2、codex diagnose settle-migration-7 を全件採用):
run に `settle_criteria.yaml` があれば (分岐は基準のファイルの有無だけ; 判定不能・照合の不一致で旧方式へ落とさない) `_handoff_gate_v2`:
停止規則 v2 の `STOP_RULE_STATE.json` が REACHED、停止の step の再評価 (`SETTLE_VERDICT.json`) の全体が SETTLED (残差・流量の収支・目的の量)、
状態・再評価・STOPPED_AT・停止の step・本段の起動 (MAIN_STAGE_STARTED の main_config の sha256・fnv と最後の forge の起動)・今の基準の sha256・
見張りの記録の照合、δ_E は停止規則が止めた step で抽出した `stop_rule_v2/delta_E_<step>/` を契約 (状態 ok・step・場と表の sha256・
抽出の設定・Euler の参照・ω) を確かめて使う (抽出し直さない; 帯の感度は記録だけ)。(d) y₁⁺・(f) 参照と ω・(g) 被覆・(h) 次の壁は今のまま。
合格したら表を旧方式と同じ名前で run の直下に写す (`delta_r_next.csv` ほか)。HANDOFF_GATE.json の `method` (v1 / v2) にどちらを使ったかを書く。

推奨 (2026-09-04, case/45 run_0026 で検証): 収束済み NS 場からの warm start では
  `--stages none --cfl 5 --implicit-relax 0.7 --steps 12000`
(soft/mid 6000 step は不要、cfl 5 + relax 0.7 で 5000〜9000 step で cfl1/24000 step と同じ残差水準・出口 M は 8000 step で凍結。
 NS 1 本 ≈ 95 s)。cold start (Euler → 中継 → 初回 NS) では `--stages full` (既定) か `--stages ramp --ramp 1,2,3.5`。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from ..metrics.deltastar import deltastar_from_core_matched_euler, massflow_ratio


# delta_r_next.csv の列 (2026-10-05, plan tooling-nozzle-cfd-pinned-initial-line §5.1 #8 — codex result M1):
# 第 2 列は**緩和・再平滑化後の次 pass の壁入力** δ_target、第 4 列が**未緩和の平滑化抽出値** δ_E。
# 2026-10-05 以前のファイルは同じ並びで名前だけ旧称 (delta_r / delta_r_use_smoothed)。列の位置は変えていないので、
# 第 2 列を壁入力として読む既存の読み手 (prepare_ns の delta_r_csv 等) の挙動は変わらない。
DELTA_R_NEXT_COLUMNS = ("x_rt", "delta_target", "delta_in_prev", "delta_E", "held")
_DELTA_R_NEXT_OLD = ("x_rt", "delta_r", "delta_in_prev", "delta_r_use_smoothed", "held")

# 受け渡しの関門 (plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #36a、#36e で §4.9 の更新に進む条件の残りを足した)
HANDOFF_GATE = "HANDOFF_GATE.json"            # 関門の判定の記録 (run ごと)
HANDOFF_CONVERGENCE = "HANDOFF_CONVERGENCE.txt"   # 関門で回した check_convergence.py --segment の出力
HANDOFF_WALL_RESOLUTION = "HANDOFF_WALL_RESOLUTION.txt"   # 関門で回した check_wall_resolution.py の出力 (§5.1 #36e)
HANDOFF_EXTRACT_DIR = "handoff_extract"       # 関門の抽出の置き場 (合格したときだけ表を run の直下へ写す; §5.1 #36e)
HANDOFF_BAND_SELECT = "edge"                  # 関門を通った後の抽出の帯の選び方 (2026-10-14 まで未指定 = 既定の adaptive)
_STOP_OVERRIDABLE = ("NOT_USED", "CENSORED", "UNAVAILABLE")   # allow_not_reached で通せる停止の状態 (それ以外の非 REACHED は通さない)
# allow_unsettled で通せる条件 (§5.1 #36h、forge-f5 と合意): 実測が閾値を超えたときだけ (各 check の measured が真)。判定不能 (記録の欠け・
# 窓を覆わない・点数の不足・非有限・δ の記録の held と欠測) は通さない
_UNSETTLED_OVERRIDABLE = ("conservation", "two_windows")

# 診断専用の run の印 (plan tooling-nozzle-isothermal-wall-chain §6 V-h370-diag、codex diagnose 2026-10-14 h370-diag-ns M1)。
# この印のある run は受け渡しの関門に回さない (資格の無い Euler 参照で回した診断の場から壁の入力を作らない)。上書きでは通さない。
DIAG_ONLY_MARKS = ("H370_DIAG_ONLY.json", "DIAGNOSTIC_ONLY.json")
# 受け渡しの抽出の固定の条件 (§5.1 #36e (4)(5)、codex diagnose 3・§6 V-h370)。停止規則の帯の生成・見張りの中の出力ごとの δ の記録・
# 関門の抽出で共通に使う (prepare_info の stop_rule.extract に写す):
#   - 一定 x の断面への補間 remap_constant_x=True、端の補間範囲の不足は欠測 remap_edge="missing" (定数で外挿しない)
#   - ω 0.5 (§6 V-h370 の登録)、平滑化は 2026-10-14 まで deltastar_loop が使っていた値 (extract_and_merge の既定) を明示して固定:
#     λ 1.0・ノット 2.0 r_t・単調性ガードの λ の昇格は ×1, ×10, ×100, ×1000 の段 (extract_and_merge の段そのもの) で、その最後 ×1000 を
#     上限として明示する (2026-10-14 までは上限の検査なし = ×1000 でも単調にならなければ黙って進んだ。今は抽出の例外 → 関門で止める)
HANDOFF_EXTRACT = {"band_select": HANDOFF_BAND_SELECT, "remap_constant_x": True, "remap_edge": "missing"}
HANDOFF_OMEGA = 0.5
HANDOFF_SMOOTH = {"smooth_lam": 1.0, "knot_spacing": 2.0, "max_lam_factor": 1000.0}
HANDOFF_WINDOW_TOL_MM = 0.01                  # 2 窓の δ_E・δ_target の最大−最小・窓の平均の差の上限 [mm] (§6 V-h370)
HANDOFF_CONV_DROP = 3.0                       # check_convergence.py --segment --drop (§6 V-h370; ツールの既定と同じ値を明示)


def handoff_extract_settings() -> dict:
    """受け渡しの抽出の固定の条件 (帯の選び方・一定 x の断面への補間と端の扱い・ω・平滑化の λ とノット・昇格の上限)。
    `extract_and_merge(run, euler_ref, out_dir=..., **handoff_extract_settings())` でそのまま渡せる形。"""
    return {**HANDOFF_EXTRACT, "omega": HANDOFF_OMEGA, **HANDOFF_SMOOTH}


class HandoffBlocked(RuntimeError):
    """受け渡しの関門が不合格 (前の pass の run から次の壁を作らない)。"""


def read_delta_r_next(path) -> dict:
    """`delta_r_next.csv` を列名で読む。戻り: x_rt / delta_target (緩和後の壁入力) / delta_in_prev /
    delta_E (未緩和の平滑化抽出) / held と header_kind ("named" か "legacy")。並びが想定と違えば拒否する。"""
    path = Path(path)
    with open(path) as f:
        head = f.readline().strip()
    names = tuple(c.strip().split(" ")[0] for c in head.split(","))
    if names[:5] == DELTA_R_NEXT_COLUMNS:
        kind = "named"
    elif names[:5] == _DELTA_R_NEXT_OLD:
        kind = "legacy"
    else:
        raise ValueError(f"{path}: delta_r_next.csv の列が想定外 ({names})")
    a = np.loadtxt(path, delimiter=",", skiprows=1)
    out = {k: a[:, i] for i, k in enumerate(DELTA_R_NEXT_COLUMNS)}
    out["header_kind"] = kind
    return out


def extract_and_merge(prev_run, euler_run, omega: float = 0.5, smooth_lam: float = 1.0,
                      knot_spacing: float = 2.0, out_dir=None, max_lam_factor: float | None = None, **kw) -> dict:
    """前 pass の NS run から抽出し、次 pass の入力 δ_r(x) を作る。
    出力 (out_dir、既定 prev_run): delta_r_equiv.csv / delta_r_equiv_diag.json / delta_r_next.csv / delta_r_extract_summary.json。
    out_dir は同じ場から別設定 (例 band_select) で抽出するとき互いに上書きしないために使う。"""
    prev_run = Path(prev_run)
    od = Path(out_dir) if out_dir is not None else prev_run
    od.mkdir(parents=True, exist_ok=True)
    d = deltastar_from_core_matched_euler(prev_run, euler_run, out_dir=od, **kw)
    d_in = d["delta_in"]
    d_use = d["delta_r_use"]
    held = ~np.isfinite(d_use)
    # 平滑化 (2026-09-04 ユーザ指摘: 3 次平滑化では壁曲率が凸凹): 抽出 raw を 5 次 P-spline (3 階差分ペナルティ,
    # ノット 2 r_t, λ=1) で「曲率が滑らか」な δ_ext(x) にしてから緩和する。hard 不合格の点は重み 0。
    from ..metrics.deltastar import smooth_delta_quintic
    # 符号付き δ_r (2026-09-12): 等温壁 (prepare_info.json の wall_thermal) では負値を許す。断熱は従来どおり ≥0 に丸める
    try:
        _wt = json.loads((prev_run / "prepare_info.json").read_text()).get("wall_thermal", {"mode": "adiabatic"})
    except Exception:
        _wt = {"mode": "adiabatic"}
    positive = (str(_wt.get("mode", "adiabatic")) == "adiabatic")
    x = d["x"]; r_inv = d["r_wall_euler"]
    raw = d["delta_r_raw"]; wts = d["hard_ok"].astype(float) & np.isfinite(raw) if False else (d["hard_ok"] & np.isfinite(raw)).astype(float)
    lam_used = None; diag = None
    for lam in (smooth_lam, smooth_lam * 10, smooth_lam * 100, smooth_lam * 1000):
        f_s, diag = smooth_delta_quintic(x, raw, weights=wts, knot_spacing=knot_spacing, lam=lam, positive=positive)
        d_ext = f_s(x)
        d_next = (1.0 - omega) * d_in + omega * d_ext
        if omega < 1.0:
            # 前回入力 δ_in (旧方式の壁など) の凸凹を引き継がないよう、緩和後も同じ P-spline を通す
            f_b, _ = smooth_delta_quintic(x, d_next, knot_spacing=knot_spacing, lam=lam, positive=positive)
            d_next = f_b(x)
        # 単調性ガード: 新しい物理壁 r_inv + δ_next がスロート下流で非単調なら λ を 10 倍して再平滑化
        r_new = r_inv + d_next
        m = x >= 0.5
        if np.all(np.diff(r_new[m]) > -1e-9):
            lam_used = lam; break
    else:
        lam_used = f"{lam} (still non-monotone)"
    # 単調性ガードの λ 昇格に上限 (2026-10-04, diagnostician): 昇格しすぎた平滑化は抽出値を追わない壁になるので pass を失敗にする
    if max_lam_factor is not None and (not isinstance(lam_used, float) and not isinstance(lam_used, int)
                                       or lam_used > smooth_lam * max_lam_factor):
        raise RuntimeError(f"extract_and_merge: 単調性ガードで λ が {lam_used} まで上がった (上限 smooth_lam×{max_lam_factor})")
    d_use = d_ext                                   # 以降の帳簿は平滑化後の抽出値で取る
    held = ~np.isfinite(d["delta_r_use"])
    np.savetxt(od / "delta_r_next.csv", np.c_[d["x"], d_next, d_in, d_use, held.astype(int)],
               delimiter=",", comments="",
               header=",".join(DELTA_R_NEXT_COLUMNS) + f" (delta_target = relaxed wall input; delta_E = unrelaxed smoothed extraction; omega={omega}; quintic P-spline knot={knot_spacing} lam={lam_used}; resid_rel_rms={diag['resid_rel_rms']:.4f})")
    fin = np.isfinite(d_use) & (d_in > 1e-4)
    ratio = d_use[fin] / d_in[fin]
    summary = {"omega": omega, "smooth": {"kind": "quintic_pspline", "knot_spacing": knot_spacing, "lam": str(lam_used), **diag},
               "max_lam_factor": max_lam_factor, "not_extracted": d.get("not_extracted"),
               "n_stations": int(len(d["x"])), "n_ok": int(d["ok"].sum()),
               "n_hard_ok": int(d["hard_ok"].sum()), "n_held": int(held.sum()),
               "use_over_in_median": float(np.median(ratio)) if fin.any() else None,
               "use_over_in_p10_p90": [float(np.percentile(ratio, 10)), float(np.percentile(ratio, 90))] if fin.any() else None,
               "max_abs_change": float(np.nanmax(np.abs(d_use - d_in))) if fin.any() else None,
               "massflow": d["massflow"],
               "delta_r_throat_use": float(np.interp(0.0, d["x"], np.where(np.isfinite(d_use), d_use, d_in))),
               "delta_r_throat_in": float(np.interp(0.0, d["x"], d_in))}
    (od / "delta_r_extract_summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def _sizing_delta_supplier(p, init_cfg: dict, rtol: float | None = None):
    """寸法の逆算の CFD 前の δ_r の供給: `prepare_ns` と同じ経路 (`integral_delta_r`: k_f の cf_scale・熱条件・5 次 P-spline の平滑化・
    `delta_r_from_table`) で、反復のたびに r_t を変えて作り直す (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.2)。
    rtol: 診断用 (同 plan §6 U2c)。`integral_delta_r(rtol=...)` へそのまま渡す (None = `integral_bl` の既定)。
    戻り: supplier(d, rt) → (delta_r_x, 記録)。記録の `integral_rtol` は solve_ivp に実際に渡った値。"""
    from ..evaluate.runner_axismach import integral_delta_r

    def supplier(d, rt):
        res, drx, info = integral_delta_r(p, d, init_cfg, scale=rt, rtol=rtol)
        return drx, {"kind": "integral_delta_r (prepare_ns と同じ経路)", "settings": {**{k: info.get(k) for k in
                     ("model", "thermal_bc", "cf_scale", "n_scale", "a_crocco", "closure", "theta0_m", "x_virtual_m")},
                     **{k: info[k] for k in ("delta_scale", "delta_offset_m", "delta_scale_slope_per_m") if k in info},   # 2・3 係数 (card-local-mesh §4.3) は使うときだけ
                     **({"closure_version": info["closure_version"]} if "closure_version" in info else {})},
                     "smooth": info.get("smooth", {}).get("kind"),
                     "integral_rtol": float(res["solve_ivp"]["rtol"]), "integral_rtol_injected": rtol is not None}
    return supplier


# --- NS 後の補正の表: 表の端の照合と前の設計 → 今の設計への写し (plan tooling-nozzle-deltastar-band-root-zone §4.2) --------------
# 抽出の表 (delta_r_next.csv) の x は前の設計 (prev run) の r_t 単位。実寸の長さを保って作り直す設計では r_t 単位の x_e も変わるので、
# 区間で正規化した座標 s (壁温分布 `wall_temperature` と同じ区切り: 入口直管 [−2,−1]・縮流部 [−1,0]・拡大部 [0,1]) で今の設計の x に写す。
# 端は端に写るので覆い漏れが起きない。表の端は抽出の float32 の丸め分だけ設計の端からずれる (#18) ので、照合したうえで揃える
ANCHOR_KEYS = ("x_in", "x_contraction", "x_throat", "x_exit")    # 区切りの点 [r_t]: 入口端・縮流部の始まり (−L_U)・設計スロート (0)・出口


def _sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def eps_f32(x: float, S: float) -> float:
    r"""表の端の照合の許容 [r_t]: $\varepsilon_{f32} = \mathrm{ulp}_{f32}(|x| S)/S + \mathrm{ulp}_{f32}(|x|)$ (格子の座標 [m] の float32 と、
    抽出が x を r_t 単位の float32 で持つ分の 2 回の丸め; plan §3 #18 の実測: 短い分 1.54e-6 r_t ≤ ε 8.6e-6 r_t)。"""
    ax, S = abs(float(x)), float(S)
    return float(np.spacing(np.float32(ax * S))) / S + float(np.spacing(np.float32(ax)))


def prev_design_anchors(prev_run) -> dict:
    """前の設計 (prev run の壁ファイル `wall_repr.json`) の区切りの点 [r_t、前の設計] と S_prev [m]。入口端・出口は `domain`、縮流部の始まり
    −L_U は `design_wall.upstream_hermite.L_U` (`design_wall.pipe.domain` の終点・`design_wall` の x_in/x_e と照合)、スロートは 0。
    S_prev は `units.scale_m` (prev run の `prepare_info.json` の scale_m と照合)。無い・欠けている・食い違う・並びが
    x_in < −L_U < 0 < x_e でなければ ValueError。戻り: 区切りの点・scale_m・壁ファイルと prepare_info.json のパスと sha256。"""
    from ..geometry.wall_axismach import WALL_FILE
    prev_run = Path(prev_run)
    path = prev_run / WALL_FILE
    if not path.exists():
        raise ValueError(f"{path} が無い — NS 後の補正の表は前の設計の区切りの点 (入口端・縮流部の始まり・スロート・出口) で今の設計へ写す "
                         "(plan tooling-nozzle-deltastar-band-root-zone §4.2)。壁ファイルのある run (physical_wall_repr: single_bspline / legacy) を使う")
    rec = json.loads(path.read_text())
    try:
        dom = [float(v) for v in rec["domain"]]
        S, unit = float(rec["units"]["scale_m"]), rec["units"]["length"]
        dw = rec["design_wall"]
        L_U = float(dw["upstream_hermite"]["L_U"])
        pipe_dom = [float(v) for v in dw["pipe"]["domain"]]
        dw_in, dw_e = float(dw["x_in"]), float(dw["x_e"])
    except (KeyError, TypeError, ValueError, IndexError) as ex:
        raise ValueError(f"{path}: 前の設計の区切りの点を読めない ({type(ex).__name__}: {ex}) — domain・units.scale_m・design_wall の x_in・x_e・"
                         "pipe.domain・upstream_hermite.L_U が要る") from ex
    if unit != "r_t":
        raise ValueError(f"{path}: 長さの単位 {unit!r} は未対応 (r_t)")
    if dom != [dw_in, dw_e] or pipe_dom != [dw_in, -L_U]:
        raise ValueError(f"{path}: domain {dom}・設計壁の [x_in, x_e] {[dw_in, dw_e]}・直管の区間 {pipe_dom} (終点は −L_U = {-L_U!r}) が食い違う")
    x_in, xc, xe = dom[0], -L_U, dom[1]
    if not (x_in < xc < 0.0 < xe):
        raise ValueError(f"{path}: 区切りの点の並びが x_in {x_in!r} < −L_U {xc!r} < 0 < x_e {xe!r} でない (入口直管の長さ 0 などは写せない)")
    pi_path = prev_run / "prepare_info.json"
    if not pi_path.exists():
        raise ValueError(f"{pi_path} が無い (S_prev の照合に要る)")
    S_pi = float(json.loads(pi_path.read_text())["scale_m"])
    if S_pi != S:
        raise ValueError(f"{path} の units.scale_m {S!r} m と {pi_path} の scale_m {S_pi!r} m が違う (別の run の壁ファイル)")
    return {"x_in": x_in, "x_contraction": xc, "x_throat": 0.0, "x_exit": xe, "scale_m": S,
            "wall_repr": str(path), "wall_repr_sha256": _sha256(path), "prepare_info": str(pi_path), "prepare_info_sha256": _sha256(pi_path)}


def design_anchors(wall) -> dict:
    """今の設計 (`design_chain` の `d["wall"]`、joint 壁) の区切りの点 [r_t]: 入口端 x_in・縮流部の始まり −L_U・スロート 0・出口 x_e
    (壁ファイルの `domain`・`design_wall.upstream_hermite.L_U` と同じ量)。"""
    try:
        x_in, L_U, xe = float(wall.x_in), float(wall.up.L_U), float(wall.x_e)
    except AttributeError as ex:
        raise ValueError(f"design_anchors: 設計壁 {type(wall).__name__} に x_in・up.L_U・x_e が無い (joint 壁だけ)") from ex
    if not (x_in < -L_U < 0.0 < xe):
        raise ValueError(f"design_anchors: 区切りの点の並びが x_in {x_in!r} < −L_U {-L_U!r} < 0 < x_e {xe!r} でない")
    return {"x_in": x_in, "x_contraction": -L_U, "x_throat": 0.0, "x_exit": xe}


def x_to_s(x, a: dict):
    """x [r_t] → 区間で正規化した座標 s (`wall_temperature.to_x_m` の逆: [x_in, −L_U] → [−2, −1]・[−L_U, 0] → [−1, 0]・[0, x_e] → [0, 1]、
    外側は同じ区間の写像を延長)。"""
    x = np.asarray(x, dtype=float)
    xi, xc, xt, xe = (float(a[k]) for k in ANCHOR_KEYS)
    return np.where(x < xc, -2.0 + (x - xi) / (xc - xi), np.where(x < xt, -1.0 + (x - xc) / (xt - xc), (x - xt) / (xe - xt)))


def s_to_x(s, a: dict):
    """s → x [r_t] (`wall_temperature.to_x_m` と同じ式・同じ区間の分け方)。"""
    s = np.asarray(s, dtype=float)
    xi, xc, xt, xe = (float(a[k]) for k in ANCHOR_KEYS)
    return np.where(s < -1.0, xi + (s + 2.0) * (xc - xi), np.where(s < 0.0, xc + (s + 1.0) * (xt - xc), xt + s * (xe - xt)))


def map_x_between_designs(x, a_prev: dict, a_new: dict) -> tuple:
    """前の設計の x → s → 今の設計の x。区切りの点が全部同じ (同じ設計) なら恒等 (コピーを返す; ビット同一)。戻り: (x, 恒等か)。"""
    if all(float(a_prev[k]) == float(a_new[k]) for k in ANCHOR_KEYS):
        return np.array(x, dtype=float, copy=True), True
    return s_to_x(x_to_s(x, a_prev), a_new), False


def snap_table_ends(x, a_prev: dict, S_prev: float) -> tuple:
    """表の始点・終点を前の設計の入口端・出口に揃える (#18; plan §4.2、codex plan レビュー M3)。**両側** |Δ| ≤ ε_f32 (`eps_f32`) なら端の
    x をその値に置き換える (値はそのまま、点は足さない)。それより大きい差 (内側・外側とも) は ValueError — 範囲の足りない表・別の設計の表は
    揃えない。戻り: (揃えた x のコピー, 記録 {start, end: 表の x・設計の x・差・ε・置き換えたか})。"""
    x = np.array(x, dtype=float, copy=True)
    S = float(S_prev)
    rec = {}
    for name, i, key, what in (("start", 0, "x_in", "入口端"), ("end", -1, "x_exit", "出口")):
        xd, xt = float(a_prev[key]), float(x[i])
        diff, eps = xt - xd, eps_f32(xd, S)
        if not abs(diff) <= eps:
            raise ValueError(f"表の{'始点' if name == 'start' else '終点'} {xt!r} r_t と前の設計の{what} {xd!r} r_t の差 {diff:+.3e} r_t "
                             f"({diff * S:+.3e} m) が ε_f32 {eps:.3e} r_t (float32 の 2 回の丸め) を超える — 表が前の設計を覆わない・"
                             "別の設計の表 (plan tooling-nozzle-deltastar-band-root-zone §4.2)")
        rec[name] = {"x_table_rt": xt, "x_design_rt": xd, "diff_rt": diff, "diff_m": diff * S, "eps_f32_rt": eps, "snapped": diff != 0.0}
        x[i] = xd
    if not np.all(np.diff(x) > 0.0):
        raise ValueError("表の端を揃えると x が狭義単調増加でなくなる")
    return x, rec


def snap_table_ends_same_design(x, a: dict, S: float) -> tuple:
    """同じ設計の壁に表を当てるとき (`prepare_ns(delta_r_csv=...)`・関門の `next_wall_check`) の端の照合: 端ごとに |Δ| ≤ ε_f32 (`eps_f32`、
    `snap_table_ends` と同じ許容) なら端の x を設計の入口端・出口に置き換える。超えた端は**そのまま** (例外にしない; 覆わない表の扱いは
    今までどおり壁の構築に任せる)。2026-10-15 case/71 run_0003_ns_cont: 表の終点が x_e より 1 ulp 短く (格子の x [m] / r_t の丸め)、
    1 本の B-spline の物理壁が「表が出口まで覆わない」で止まった (plan tooling-nozzle-card-local-mesh §5.1 #3h)。
    戻り: (x のコピー, 記録 {start, end: 表の x・設計の x・差・ε・置き換えたか})。"""
    x = np.array(x, dtype=float, copy=True)
    S = float(S)
    rec = {}
    for name, i, key in (("start", 0, "x_in"), ("end", -1, "x_exit")):
        xd, xt = float(a[key]), float(x[i])
        diff, eps = xt - xd, eps_f32(xd, S)
        snap = diff != 0.0 and abs(diff) <= eps
        rec[name] = {"x_table_rt": xt, "x_design_rt": xd, "diff_rt": diff, "eps_f32_rt": eps, "snapped": bool(snap)}
        if snap:
            x[i] = xd
    if not np.all(np.diff(x) > 0.0):
        raise ValueError("表の端を揃えると x が狭義単調増加でなくなる")
    return x, rec


def load_ns_correction(prev_run, delta_next=None) -> dict:
    """NS 後の補正の入力を読んで検査する (plan tooling-nozzle-deltastar-band-root-zone §4.2・§4.3 3、codex plan レビュー M1・M3):
    prev run の `delta_r_next.csv` (delta_next で別の抽出を指定可) の δ_E と held、前の設計の区切りの点と S_prev (`prev_design_anchors`)。

    - δ_E に非有限値があれば ValueError。表の始点・終点を前の設計の入口端・出口に揃える (`snap_table_ends`、合わなければ ValueError)。
    - **使用域 [0, x_e] (前の設計の出口) に held (hard 不合格で前回値を保持した点) が 1 点でもあれば ValueError** (平滑化は hard の点を重み 0 で
      飛ばして全点で値を作るので、δ_E の列だけを見ると不合格の点が素通りする)。使用域の外 (x < 0) の held は記録だけ。
      held の列が無い読み (旧形式) は検査しない。
    戻り: x (揃えた後、前の設計の r_t 単位)・delta_E・S_prev_m・anchors_prev・ends・held・出典 (パスと sha256)。"""
    prev_run = Path(prev_run)
    dn = Path(delta_next) if delta_next is not None else prev_run / "delta_r_next.csv"
    if not dn.exists():
        raise ValueError(f"solve_rt_throat (NS 後): {dn} が無い (extract_and_merge で作る — δ_E の全分布が要る)")
    nx = read_delta_r_next(dn)
    dE = np.asarray(nx["delta_E"], dtype=float)
    if not np.all(np.isfinite(dE)):
        raise ValueError(f"{dn} の delta_E に非有限値がある")
    a_prev = prev_design_anchors(prev_run)
    S_prev = a_prev["scale_m"]
    x, ends = snap_table_ends(np.asarray(nx["x_rt"], dtype=float), a_prev, S_prev)
    held_rec = None
    if nx.get("held") is not None:
        held = np.asarray(nx["held"], dtype=float) != 0.0
        use = held & (x >= 0.0) & (x <= a_prev["x_exit"])
        if use.any():
            xs = x[use]
            raise ValueError(f"{dn}: 使用域 [0, x_e = {a_prev['x_exit']:.9g}] r_t に held (hard 不合格で前回値を保持した点) が {int(use.sum())} 点 "
                             f"(x = {', '.join(f'{v:.4f}' for v in xs[:10])}{' …' if len(xs) > 10 else ''} r_t) — 不合格の点の δ_E "
                             "(平滑化の値) で寸法と壁を作らない (plan tooling-nozzle-deltastar-band-root-zone §4.1、codex plan レビュー M1)")
        held_rec = {"n_held": int(held.sum()), "n_held_in_use": 0, "x_held_outside_use_rt": [float(v) for v in x[held][:20]]}
    return {"prev_run": str(prev_run), "delta_next": str(dn), "delta_next_sha256": _sha256(dn), "header_kind": nx["header_kind"],
            "x": x, "delta_E": dE, "S_prev_m": S_prev, "anchors_prev": {k: a_prev[k] for k in ANCHOR_KEYS}, "ends": ends,
            "held": held_rec, "wall_repr": a_prev["wall_repr"], "wall_repr_sha256": a_prev["wall_repr_sha256"],
            "prepare_info": a_prev["prepare_info"], "prepare_info_sha256": a_prev["prepare_info_sha256"]}


def ns_correction_table(corr: dict, anchors_new: dict, rt: float) -> dict:
    r"""補正の表: x は前の設計 → s → 今の設計 (`map_x_between_designs`)、値は今の `solve_rt_throat --prev` と同じ大きさ
    $\delta = \delta_E\,(r_t/S_{prev})^{-0.2}$。戻り: x・delta_r・factor・identity (写しが恒等か)。"""
    fac = (rt / corr["S_prev_m"]) ** -0.2
    x_new, ident = map_x_between_designs(corr["x"], corr["anchors_prev"], anchors_new)
    if not np.all(np.diff(x_new) > 0.0):
        raise ValueError("ns_correction_table: 今の設計へ写した x が狭義単調増加でない")
    return {"x": x_new, "delta_r": corr["delta_E"] * fac, "factor": fac, "identity": ident}


def write_ns_correction_table(path, tbl: dict, corr: dict, rt: float, anchors_new: dict) -> None:
    """補正の表を書く (列 x_rt, delta_r; `prepare_ns(delta_r_csv=...)` が読む形)。寸法の逆算 (`solve_rt_throat`) と
    要求カードの NS 補正 (`wt_card --ns-correction`) の共通の書き手。"""
    e = corr["ends"]
    ap_, an_ = (", ".join(repr(float(a[k])) for k in ANCHOR_KEYS) for a in (corr["anchors_prev"], anchors_new))
    np.savetxt(path, np.c_[tbl["x"], tbl["delta_r"]], delimiter=",", comments="",
               header=f"x_rt,delta_r (= delta_E x (r_t/S_prev)^-0.2; r_t={rt!r} m, S_prev={corr['S_prev_m']!r} m, factor={tbl['factor']!r}; "
                      f"from {corr['delta_next']} col delta_E; x mapped from the previous design (x_in, -L_U, 0, x_e) = ({ap_}) r_t to this design "
                      f"({an_}) r_t via s{' [identity]' if tbl['identity'] else ''}; table ends snapped to the previous design: "
                      f"start {e['start']['diff_rt']:+.3e} r_t, end {e['end']['diff_rt']:+.3e} r_t (eps_f32 {e['start']['eps_f32_rt']:.3e} / "
                      f"{e['end']['eps_f32_rt']:.3e}); solve_rt_throat)")


def ns_correction_record(corr: dict, anchors_new: dict) -> dict:
    """NS 後の補正の出典と写しの記録 (JSON にそのまま書ける形)。"""
    return {"prev_run": corr["prev_run"], "delta_next": corr["delta_next"], "delta_next_sha256": corr["delta_next_sha256"],
            "delta_column": f"delta_E ({corr['header_kind']} header)",
            "wall_repr": corr["wall_repr"], "wall_repr_sha256": corr["wall_repr_sha256"],
            "prepare_info": corr["prepare_info"], "prepare_info_sha256": corr["prepare_info_sha256"],
            "S_prev_m": corr["S_prev_m"], "anchors_prev_rt": dict(corr["anchors_prev"]), "anchors_new_rt": dict(anchors_new),
            "x_map": ("identity (同じ区切りの点)" if all(float(corr["anchors_prev"][k]) == float(anchors_new[k]) for k in ANCHOR_KEYS)
                      else "前の設計の x → s (区間で正規化; wall_temperature と同じ区切り) → 今の設計の x"),
            "table_ends": corr["ends"], "held": corr["held"]}


SIZING_TOL_M = 1e-9     # 寸法の逆算の数値解法の許容差 [m] (plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U2・U2b の合格条件と同じ値)。
                         # 既定で使うのは NS 後の `solve_rt_throat` (δ_E の全分布) と試験用の `_delta_supplier`。CFD 前の既定は下の 2 つ
SIZING_MAX_ITER = 30
# CFD 前の寸法 (積分法の δ_r) は**初期見積もり**: その既定の許容差 [m] (plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U2d、
# 諮問 notes/reviews/2026-10-07-upoly-sizing-after-u2c-diagnose.md の P1)。CFD 前の δ_r は r_t に対して滑らかでなく、局所の散らばり
# (rtol 1e-6: スロート半径 7.41e-9 m・出口半径 2.64e-6 m、U2c) が SIZING_TOL_M より大きくて届かないため。真値の誤差の保証ではなく、
# 精密な寸法・感度の評価には未検証 (同 plan §6 既知の制約)。SIZING_TOL_M (共通) は変えない
SIZING_TOL_PRE_CFD_THROAT_M = 1e-7
SIZING_TOL_PRE_CFD_EXIT_M = 1e-5


def _resolve_sizing_tol(tol_R_m, target: str, cfd_before: bool) -> tuple[float, str]:
    """寸法の逆算の実効の許容差 [m] とその出典。明示 (tol_R_m) が優先。既定は CFD 前なら初期見積もりの定数 (target 別)、
    それ以外 (NS 後・試験用の供給) は SIZING_TOL_M。"""
    if tol_R_m is not None:
        if isinstance(tol_R_m, bool) or not np.isfinite(float(tol_R_m)) or float(tol_R_m) <= 0.0:
            raise ValueError(f"tol_R_m は正の有限の数値 [m] (受け取った値: {tol_R_m!r})")
        return float(tol_R_m), "explicit"
    if cfd_before:
        if target == "exit":
            return SIZING_TOL_PRE_CFD_EXIT_M, "default: SIZING_TOL_PRE_CFD_EXIT_M (CFD 前の初期見積もり)"
        return SIZING_TOL_PRE_CFD_THROAT_M, "default: SIZING_TOL_PRE_CFD_THROAT_M (CFD 前の初期見積もり)"
    return SIZING_TOL_M, "default: SIZING_TOL_M"


class SizingNotConverged(RuntimeError):
    """寸法の逆算が反復の上限で許容差に収まらなかった (不合格)。`history` に各反復の r_t・残差、`tol_R_m` に実効の許容差を持つ。"""

    def __init__(self, msg: str, history: list, tol_R_m: float | None = None) -> None:
        super().__init__(msg)
        self.history = history
        self.tol_R_m = tol_R_m


def _fixed_point_rt(p, d, R_target_m: float, target: str, supplier, rt0: float, max_iter: int, tol_R_m: float) -> dict:
    """r_t を固定点反復で解く: 反復ごとに δ_r を作り直して物理壁を組み (`build_physical_wall`、prepare_ns と同じ)、
    target = 'exit' なら r_t·r_W(x_e)、'throat' なら r_t·min r_W (物理壁の大域最小) を目標に合わせる。
    各反復で今の r_t の残差 rt·r(rt) − R を評価し、|残差| ≤ tol_R_m で止める (返す r_t の残差は最後に評価した値そのもの)。
    max_iter 回で収まらなければ例外 (plan §4.2「反復の上限に達したら不合格」)。"""
    from ..evaluate.runner_axismach import build_physical_wall
    rt = float(rt0)
    hist = []
    for k in range(int(max_iter)):
        drx, src = supplier(d, rt)
        W = build_physical_wall(p, d, rt, delta_r_x=drx, offset="radial")
        x_e = float(W.x_e)
        r_e = float(W.r(np.array([x_e]))[0])
        r_val = r_e if target == "exit" else float(W.r_throat)
        res = rt * r_val - R_target_m
        hist.append({"iter": k, "r_t_m": rt, "r_rt": r_val, "residual_m": res, "x_throat_rt": float(W.x_throat),
                     "r_throat_rt": float(W.r_throat), "exit_r_rt": r_e, "delta_r_exit_rt": float(drx(np.array([x_e]))[0]),
                     "delta_r_x0_rt": float(drx(np.array([0.0]))[0])})
        if abs(res) <= tol_R_m:
            return {"r_t_m": rt, "residual_m": res, "wall": W, "delta_r_x": drx, "delta_r_source": src, "iters": hist,
                    "n_iter": k + 1, "tol_R_m": float(tol_R_m)}
        rt = R_target_m / r_val
    res = [h["residual_m"] for h in hist]
    raise SizingNotConverged(f"寸法の逆算 ({target}): {max_iter} 回で |r_t·r − R| ≤ {tol_R_m:g} m に収まらない "
                             f"(最後の残差 {res[-1]:.3e} m、後半の |残差| の範囲 {min(map(abs, res[len(res) // 2:])):.2e}〜"
                             f"{max(map(abs, res[len(res) // 2:])):.2e} m) — 不合格 (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.2)",
                             hist, tol_R_m=float(tol_R_m))


def _sizing_result(fp: dict, p, R_target_m: float, target: str, tol_source: str) -> dict:
    """`solve_rt`・`solve_rt_throat` の戻りの共通部分: 解いた r_t、最終の寸法での残差、実効の許容差とその出典、物理スロートの位置と半径、
    出口半径、δ_r の出典。"""
    W, rt = fp["wall"], fp["r_t_m"]
    x_e = float(W.x_e)
    return {"r_t_m": rt, "sizing": target, "target_m": R_target_m, "residual_m": fp["residual_m"], "n_iter": fp["n_iter"],
            "tol_R_m": fp["tol_R_m"], "tol_R_m_source": tol_source,
            "pw_upstream": W.pw_upstream, "physical_wall_repr": type(W).__name__,
            "physical_throat": {"x_rt": float(W.x_throat), "r_rt": float(W.r_throat), "kappa_rt": float(W.kappa_throat),
                                "x_m": float(W.x_throat) * rt, "r_m": float(W.r_throat) * rt},
            "exit_radius_m": float(W.r(np.array([x_e]))[0]) * rt, "x_e_rt": x_e, "length_m": x_e * rt,
            "delta_r_source": fp["delta_r_source"], "iter_history": fp["iters"]}


def solve_rt(problem, R_exit_m: float, prev_run=None, euler_run=None, n_iter: int = 6,
             max_iter: int = SIZING_MAX_ITER, tol_R_m: float | None = None, integral_rtol: float | None = None) -> dict:
    r"""出口の物理半径 $R$ を仕様に合わせる **スロート半径 $r_t$ の 1 変数解**。

    設計は $r_t$ 無次元で不変なので $R = r_t\,r_W(x_e; r_t)$ の $r_t$ だけを解く。
    - prev_run なし (CFD 前): `prepare_ns` と同じ δ_r の経路 (`integral_delta_r`: problem の `deltastar_initializer` の k_f・熱条件・
      平滑化) と同じ壁の構築 (`build_physical_wall`、pw_upstream に従う) で、反復のたびに δ_r と物理壁を作り直し、
      $r_t\,r_W(x_e) = R$ を |残差| ≤ tol_R_m まで解く (2026-10-07 変更、plan tooling-nozzle-upstream-poly-and-throat-sizing §4.2:
      以前は未較正・未平滑化の `integral_bl` を直接呼んでいて、生産の壁と δ_r が違った)。max_iter 回で収まらなければ例外
      (`SizingNotConverged`)。tol_R_m の既定 (None) は `SIZING_TOL_PRE_CFD_EXIT_M` (1e-5 m)、明示で渡せる。
      **CFD 前の寸法は初期見積もり。精密な寸法・感度の評価には未検証** (同 plan §6 U2d・既知の制約): `integral_bl` (RK45 rtol 1e-6) の
      δ_r は r_t に対して滑らかでなく (r_t を 1e-12 m 変えただけで δ_r(x_e) が ~1e-5 r_t 揺れる)、r_t の ±1e-10 m の 11 点での出口半径の
      局所の散らばりは 2.64e-6 m (1 次の傾向を除いた max − min、U2c `_band_ab/upoly/U2c.json`。rtol 1e-10 でも 5.09e-8 m)。
      `SIZING_TOL_M` (1e-9 m) には届かない (U2b)。原因は未確定、真値の誤差は未評価。
    - prev_run あり (NS 後の最終補正、変えていない): その run の抽出 δ_r(x_F) を使い、$r_t$ 依存は $Re^{-0.2}$ で補正。n_iter 回の反復。
      **この経路は固定回数 (n_iter) の反復で、tol_R_m も max_iter も判定しない** (渡しても使わない。残差も返さない)。NS 後の 1e-9 m の
      保証は `solve_rt_throat` の δ_E の全分布の経路だけ (同 plan §6「NS 後の経路の保証の範囲の訂正」)。
    - integral_rtol: **診断用** (plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U2c)。CFD 前の経路の `integral_bl` の RK45 の
      相対許容差を明示の引数で注入する (None = 既定 1e-6、今の振る舞いのまま)。積分を使わない NS 後の経路で指定したら例外
      (渡っていない口を作らない)。実効値は戻りの `delta_r_source.integral_rtol` (solve_ivp に渡った値)。
    戻り: dict(r_t_m, delta_exit_rt, source, delta_column [prev_run ありで読んだ列], iters, …)。CFD 前は残差・実効の許容差
    (`tol_R_m`・`tol_R_m_source`)・物理スロート・出口半径・δ_r の出典も (`_sizing_result`)。"""
    from ..probdef import load_problem
    from ..evaluate.runner_axismach import design_chain
    if integral_rtol is not None and prev_run is not None:
        raise ValueError("solve_rt: integral_rtol は CFD 前の経路 (prev_run なし) 専用 — NS 後の経路は積分法を使わない")
    p = load_problem(problem); d = design_chain(p)
    S0 = float(p.spec["r_throat"])
    rF = float(d["wall_inv"][-1, 1]); xF = float(d["wall_inv"][-1, 0])
    hist = []
    if prev_run is None:
        init = resolve_integral_initializer(problem)
        tol, tol_src = _resolve_sizing_tol(tol_R_m, "exit", cfd_before=True)
        fp = _fixed_point_rt(p, d, float(R_exit_m), "exit", _sizing_delta_supplier(p, init, rtol=integral_rtol), S0, max_iter, tol)
        out = _sizing_result(fp, p, float(R_exit_m), "exit", tol_src)
        last = fp["iters"][-1]
        out.update(r_t_prev_m=S0, delta_exit_rt=last["delta_r_exit_rt"], r_F_rt=rF, x_F_rt=xF, R_exit_m=R_exit_m,
                   source="integral_delta_r (prepare_ns と同じ経路: deltastar_initializer の k_f・熱条件・5 次 P-spline 平滑化) + "
                          "build_physical_wall", delta_column=None, delta_measured_xF_rt=None,
                   iters=[(h["r_t_m"], h["delta_r_exit_rt"], h["residual_m"]) for h in fp["iters"]],
                   initializer=init)
        return out
    else:
        prev_run = Path(prev_run)
        if not (prev_run / "delta_r_equiv.csv").exists():
            if euler_run is None: raise ValueError("prev_run に抽出結果が無く euler_run も未指定")
            deltastar_from_core_matched_euler(prev_run, euler_run, out_dir=prev_run)
        # 未緩和の平滑化抽出値 δ_E (delta_r_next.csv の delta_E 列) の x_F 端を使う。2026-10-05 訂正 (codex result M1):
        # 以前は第 2 列 (ω 緩和後の次 pass 壁入力 δ_target) を読んでいて、出口の δ を約半分しか動かしていなかった。
        # delta_r_next.csv が無ければ生抽出の x_F−0.3
        if (prev_run / "delta_r_next.csv").exists():
            nx = read_delta_r_next(prev_run / "delta_r_next.csv")
            d_meas = float(np.interp(xF, nx["x_rt"], nx["delta_E"]))
            d_col = f"delta_E (unrelaxed smoothed extraction; delta_r_next.csv {nx['header_kind']} header, col 4)"
        else:
            e = np.genfromtxt(prev_run / "delta_r_equiv.csv", delimiter=",", names=True)
            d_meas = float(np.interp(xF - 0.3, e["x_rt"], e["delta_r_raw"]))
            d_col = "delta_r_raw at x_F-0.3 (delta_r_equiv.csv)"
        S_prev = float(json.loads((prev_run / "prepare_info.json").read_text())["scale_m"])
        rt = S_prev
        for k in range(n_iter):
            de = d_meas * (rt / S_prev) ** -0.2
            rt_new = R_exit_m / (rF + de); hist.append((rt, de, rt_new)); rt = rt_new
        src = f"measured delta_r from {prev_run.name} (Re^-0.2 scaling)"
    return dict(r_t_m=float(rt), r_t_prev_m=S0, delta_exit_rt=float(hist[-1][1]), r_F_rt=rF, x_F_rt=xF,
                R_exit_m=R_exit_m, source=src, delta_column=d_col, delta_measured_xF_rt=d_meas,
                iters=hist, length_m=float(xF * rt))


def solve_rt_throat(problem, R_throat_m: float, prev_run=None, delta_r_out=None, delta_next=None,
                    max_iter: int = SIZING_MAX_ITER, tol_R_m: float | None = None, _delta_supplier=None) -> dict:
    r"""**物理スロート径から寸法を決める**: $r_t\cdot\min_x r_W(x; r_t) = R_{throat}$ を $r_t$ について解く
    (plan tooling-nozzle-upstream-poly-and-throat-sizing §4.2、ユーザ決定 2026-10-07)。$r_W$ は物理壁 (r_t 単位、pw_upstream に従う)、
    最小は物理壁の大域最小 (`wall_global_min`)。スロート径を固定するときは出口径を同時に固定条件にしない (出口半径は戻り値に記録するだけ)。

    - prev_run なし (CFD 前): 反復のたびに `prepare_ns` と同じ経路 (`integral_delta_r`: k_f の cf_scale・熱条件・5 次 P-spline の
      平滑化・`delta_r_from_table`) で δ_r を作り、同じ壁の構築 (`build_physical_wall`) で最小半径を求める。k_f を較正し直したら
      寸法も解き直すこと (k_f はスロートの δ_r も変える)。
    - prev_run あり (NS 後): その run の `delta_r_next.csv` (delta_next で別の抽出 [例 `_extract_edge/delta_r_next.csv`] を指定可。
      S_prev は prev_run の prepare_info.json の scale_m) の **δ_E の全分布** (未緩和の平滑化抽出) から補正関数
      $\delta(x; r_t) = \delta_E(x)\,(r_t/S_{prev})^{-0.2}$ を作り、同じ関数で寸法と壁を作る ($Re^{-0.2}$ の換算は予測の近似)。
      解いた r_t での補正関数の表を delta_r_out (必須) に書く — 次の壁はこの表を `prepare_ns(delta_r_csv=...)` に渡して作る
      (寸法と壁で同じ関数)。delta_r_next.csv が無ければ例外。2026-10-13 (plan tooling-nozzle-deltastar-band-root-zone §4.2、
      codex plan レビュー M1・M3; `load_ns_correction`): 使用域 [0, x_e] に held (hard 不合格) の点があれば例外。表の始点・終点を
      前の設計 (prev run の `wall_repr.json`、無ければ例外) の入口端・出口と照合し、両側 ε_f32 以内なら揃える (超えれば例外)。表の x は
      前の設計 → 区間で正規化した s → この問題の設計へ写す (同じ設計なら恒等)。戻りの `ns_correction` に出典・揃えた差・写しを記録する。
    - 反復: 各反復で今の r_t の残差 r_t·min r_W − R を評価し、|残差| ≤ tol_R_m で止める (返す残差 = 最終の寸法で評価し直した値)。
      max_iter 回で収まらなければ例外 (`SizingNotConverged`、不合格)。tol_R_m の既定 (None) は経路で違う: CFD 前は
      `SIZING_TOL_PRE_CFD_THROAT_M` (1e-7 m)、NS 後と `_delta_supplier` は `SIZING_TOL_M` (1e-9 m、今までどおり)。明示で渡せる。
      **CFD 前の寸法は初期見積もり。精密な寸法・感度の評価には未検証** (plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U2d・既知の
      制約): CFD 前の δ_r (`integral_bl`、RK45 rtol 1e-6) は r_t を 1e-12 m 変えただけで出口で ~1e-5 r_t・スロートで ~1e-7 r_t 揺れ、
      r_t の ±1e-10 m の 11 点での物理スロート半径の局所の散らばりは 7.41e-9 m (U2c)。1e-9 m には届かない。原因は未確定、真値の誤差は未評価。
    - `_delta_supplier` (試験用): supplier(d, rt) → (delta_r_x, 記録) で δ_r の供給経路だけを差し替える (判別 A/B の腕 A)。
    戻り: r_t_m・残差・実効の許容差 (`tol_R_m`・`tol_R_m_source`)・反復の履歴・δ_r の出典と設定・物理スロートの位置と半径・
    そのとき決まる出口半径 ほか。"""
    from ..probdef import load_problem
    from ..evaluate.runner_axismach import delta_r_from_table, design_chain
    corr = None
    if _delta_supplier is None and prev_run is not None:
        if delta_r_out is None:
            raise ValueError("solve_rt_throat (NS 後): delta_r_out (次の壁に渡す補正関数の表の出力先) が必要 — 寸法と壁を同じ関数で作るため")
        # 入力の検査 (使用域の held・δ_E の非有限・表の端と前の設計の照合) と前の設計の区切りの点 (plan tooling-nozzle-deltastar-band-root-zone
        # §4.2、codex plan レビュー M1・M3)。設計チェーンの前に検査する
        corr = load_ns_correction(prev_run, delta_next)
    p = load_problem(problem); d = design_chain(p)
    S0 = float(p.spec["r_throat"])
    R = float(R_throat_m)
    extra = {}
    if _delta_supplier is not None:
        supplier, rt0 = _delta_supplier, S0
        extra["delta_r_route"] = "test supplier (_delta_supplier)"
        tol, tol_src = _resolve_sizing_tol(tol_R_m, "throat", cfd_before=False)
    elif prev_run is None:
        init = resolve_integral_initializer(problem)
        supplier, rt0 = _sizing_delta_supplier(p, init), S0
        extra.update(delta_r_route="cfd_before: integral_delta_r (prepare_ns と同じ経路)", initializer=init)
        tol, tol_src = _resolve_sizing_tol(tol_R_m, "throat", cfd_before=True)
    else:
        tol, tol_src = _resolve_sizing_tol(tol_R_m, "throat", cfd_before=False)
        # 表の x は前の設計 → s → この問題の設計 (同じ設計なら恒等) へ写す
        a_new = design_anchors(d["wall"])
        out_csv = Path(delta_r_out)

        def supplier(d_, rt, _c=corr, _a=a_new, _out=out_csv):
            # 表を書いて読み直し、prepare_ns(delta_r_csv=...) と同じ読み方・同じ関数で壁を作る (寸法と壁で同じ関数)
            tb = ns_correction_table(_c, _a, rt)
            write_ns_correction_table(_out, tb, _c, rt, _a)
            tbl = np.loadtxt(_out, delimiter=",", skiprows=1)
            return delta_r_from_table(tbl[:, 0], tbl[:, 1]), {
                "kind": "delta_E 全分布 × (r_t/S_prev)^-0.2 (Re^-0.2 換算は予測の近似)", "prev_run": _c["prev_run"],
                "delta_next": _c["delta_next"], "delta_column": f"delta_E ({_c['header_kind']} header)", "S_prev_m": _c["S_prev_m"],
                "factor": tb["factor"], "table": str(_out), "x_map_identity": tb["identity"]}
        rt0 = corr["S_prev_m"]
        extra.update(delta_r_route="ns_after: delta_E 全分布", delta_r_out=str(out_csv), ns_correction=ns_correction_record(corr, a_new))
    fp = _fixed_point_rt(p, d, R, "throat", supplier, rt0, max_iter, tol)
    out = _sizing_result(fp, p, R, "throat", tol_src)
    out.update(r_t_problem_m=S0, R_throat_m=R, **extra)
    return out


def stop_state(run_dir) -> dict:
    """run の本段の停止の状態 (plan §5.1 #36a)。`STOP_RULE_STATE.json` (run_staged_ns が本段の後に書く) があればその status。
    無ければ prepare_info.json の stop_rule を見て、要求していれば MISSING (本段が終わっていないか記録が失われた)、要求していない
    (`evaluate.stop_rule: false`・`--no-stop-rule` で準備・2026-10-14 より前の準備; キー無しの問題は統合後は要求する) なら NOT_USED。
    run_staged_ns の終了コードは見ない (CENSORED・一部の UNAVAILABLE でも 0 なので)。"""
    from ..evaluate import stop_rule as SR
    run_dir = Path(run_dir)
    st = SR.read_state(run_dir)
    if st is not None:
        return {"status": st.get("status"), "source": SR.STATE, "stop_step": st.get("stop_step"), "forge_rc": st.get("forge_rc")}
    try:
        sr = json.loads((run_dir / "prepare_info.json").read_text()).get("stop_rule")
    except (OSError, ValueError):
        sr = None
    if isinstance(sr, dict) and sr.get("requested"):
        return {"status": "MISSING", "source": f"prepare_info.json は停止規則を要求しているのに {SR.STATE} が無い "
                                                "(本段が終わっていないか記録が無い)"}
    return {"status": "NOT_USED", "source": "prepare_info.json (停止規則を使わない run)"}


def convergence_verdict(run_dir) -> dict:
    """`check_convergence.py --segment --drop 3` (stage_manifest.json の最後の区間 = 方程式・BC・離散化が同じ本段、全列) の判定 (plan §5.1 #36a;
    `--drop 3` はツールの既定と同じ値を §6 V-h370 の登録どおり明示 #36e)。
    kind: PASS (終了コード 0) / DIVERGED (NaN/Inf) / UNDETERMINED (判定不能: stage_manifest・残差の履歴が無い ほか) / NOT_CONVERGED。
    2026-10-14 までの `run_pass` ほかは区間を指定せずに回して結果を捨てていた。"""
    from ..evaluate.runner import FORGE_TOOLS
    cmd = [sys.executable, str(FORGE_TOOLS / "check_convergence.py"), "--segment", "--drop", f"{HANDOFF_CONV_DROP:g}", str(run_dir)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    out = r.stdout + r.stderr
    lines = [ln.strip() for ln in out.splitlines()
             if ln.startswith("===") or ln.startswith("OVERALL") or "判定不能" in ln or "NO residual_history" in ln or "empty residual" in ln]
    if r.returncode == 0:
        kind = "PASS"
    elif any("-> DIVERGED" in ln for ln in lines):
        kind = "DIVERGED"
    elif r.returncode != 1 or any(("判定不能" in ln or "NO residual_history" in ln or "empty residual" in ln) for ln in lines):
        kind = "UNDETERMINED"
    else:
        kind = "NOT_CONVERGED"
    return {"cmd": cmd, "rc": r.returncode, "kind": kind, "verdict": lines, "output": out}


# --- 受け渡しの関門の残りの条件 (§4.9 の「更新に進む条件」、§5.1 #36e (5)・§6 V-h370) ------------------------------------------
def wall_resolution_verdict(run_dir, step=None) -> dict:
    """局所 y₁⁺ (`check_wall_resolution.py`、目標は局所 y₁⁺ ≤ 1・超過の許容はツールの既定のまま)。step があれば `--step step` (停止の出力の
    壁のダンプ)。kind: PASS (終了コード 0) / FAIL (1) / UNDETERMINED (それ以外 = 判定不能; 合格ではない)。"""
    from ..evaluate.runner import FORGE_TOOLS
    cmd = [sys.executable, str(FORGE_TOOLS / "check_wall_resolution.py"), str(run_dir)] + (["--step", str(int(step))] if step is not None else [])
    r = subprocess.run(cmd, capture_output=True, text=True)
    out = r.stdout + r.stderr
    kind = {0: "PASS", 1: "FAIL"}.get(r.returncode, "UNDETERMINED")
    return {"cmd": cmd, "rc": r.returncode, "kind": kind, "verdict": [ln.strip() for ln in out.splitlines() if ln.startswith("VERDICT")],
            "output": out}


def _read_stop_record(run_dir) -> dict | None:
    from ..evaluate import stop_rule as SR
    p = Path(run_dir) / SR.RECORD
    try:
        return json.loads(p.read_text()) if p.exists() else None
    except (OSError, ValueError):
        return None


def conservation_check(rec: dict | None, step) -> dict:
    """保存則 (質量の収支): stop_rule.json の行の r = (ṁ_in − ṁ_out)/ṁ_in の、停止の出力 s で終わる窓 [s − W, s] の平均と振幅
    ((最大 − 最小)/2) が停止規則の閾値 (r_mean_max・r_amp_max) 以内。窓が本段の履歴に覆われ (最初の出力 ≤ s − W)、min_points 点以上で、
    r が全部有限であること。記録が無い・s が無い・足りないときは不合格 (判定不能は合格ではない)。
    measured: 平均と振幅を測れた (不合格なら閾値の超過; 関門の allow_unsettled はこのときだけ通す、§5.1 #36h)。"""
    if rec is None or step is None:
        return {"pass": False, "why": "stop_rule.json か停止の出力 s が無い (停止規則を使っていない run は保存則を判定できない)"}
    st = rec.get("settings") or {}
    W, mp = int(st.get("window", 0)), int(st.get("min_points", 0))
    rows = [r for r in rec.get("rows") or [] if isinstance(r.get("step"), int)]
    out = {"step": int(step), "window": W, "min_points": mp, "r_mean_max": st.get("r_mean_max"), "r_amp_max": st.get("r_amp_max"),
           "pass": False}
    if not rows or not W or rows[0]["step"] > step - W:
        out["why"] = f"本段の履歴が窓 [s − {W}, s] を覆わない"
        return out
    win = [r for r in rows if step - W <= r["step"] <= step]
    rv = np.array([r.get("r", np.nan) for r in win], dtype=float)
    out["n"] = len(win)
    if len(win) < mp or not np.all(np.isfinite(rv)):
        out["why"] = f"窓の点数 {len(win)} (≥ {mp} が要る) か r に非有限値"
        return out
    out["r_mean"], out["r_amp"] = float(rv.mean()), float((rv.max() - rv.min()) / 2.0)
    out["measured"] = True
    out["pass"] = bool(abs(out["r_mean"]) <= float(st["r_mean_max"]) and out["r_amp"] <= float(st["r_amp_max"]))
    if not out["pass"]:
        out["why"] = f"|r の平均| {abs(out['r_mean']):.3e} (≤ {st['r_mean_max']}) か r の振幅 {out['r_amp']:.3e} (≤ {st['r_amp_max']})"
    return out


def two_window_check(rec: dict | None, step, tol_mm: float = HANDOFF_WINDOW_TOL_MM) -> dict:
    """δ_E・δ_target の隣り合う 2 窓の条件 (§6 V-h370): 停止の出力 s に対して W1 = [s − 2W, s − W]・W2 = [s − W, s] (閉区間、W は停止規則の
    窓)。stop_rule.json の行の `delta` (見張りの中の出力ごとの抽出、`stop_rule.delta_snapshot`) を読み、計算の後に抽出し直さない。
    各窓の全使用域 (delta_series.x の全列) で δ_E・δ_target の最大 − 最小 ≤ tol_mm、窓の平均どうしの差も全列で ≤ tol_mm [mm]。
    各窓は min_points 点以上 (停止規則と同じ)、窓の全出力に欠測・held・抽出の失敗・飛ばしの無い記録があること。満たさなければ不合格。
    窓に入るべき出力は記録にある行からではなく、出力の間隔 (settings.out_step_interval) から `stop_rule.required_window_steps` で列挙し、
    行そのものが欠けた出力も不合格にする (§5.1 #36g、codex diagnose 2026-10-14 の Major「判定窓の完全性」; V-h370 では 17 出力)。
    measured: 両方の量の全部の振れを測れた (不合格なら閾値の超過; 関門の allow_unsettled はこのときだけ通す、§5.1 #36h)。"""
    from ..evaluate.stop_rule import required_window_steps
    if rec is None or step is None:
        return {"pass": False, "why": "stop_rule.json か停止の出力 s が無い (出力ごとの δ の記録が要る)"}
    st, ds = rec.get("settings") or {}, rec.get("delta_series") or {}
    W, mp = int(st.get("window", 0)), int(st.get("min_points", 0))
    xs = np.asarray(ds.get("x") or [], dtype=float)
    S = ds.get("scale_m")
    out = {"step": int(step), "windows": [[step - 2 * W, step - W], [step - W, step]], "tol_mm": float(tol_mm), "n_cols": int(len(xs)),
           "pass": False}
    rows = [r for r in rec.get("rows") or [] if isinstance(r.get("step"), int)]
    if not W or not len(xs) or S is None:
        out["why"] = "stop_rule.json に窓・delta_series (列の x・scale_m) が無い"
        return out
    if not rows or rows[0]["step"] > step - 2 * W:
        out["why"] = f"本段の履歴が 2 窓 [s − {2 * W}, s] を覆わない (最初の出力 {rows[0]['step'] if rows else None})"
        return out
    oi = st.get("out_step_interval")
    try:
        need = required_window_steps(int(step), W, oi)
    except (TypeError, ValueError) as e:
        out["why"] = f"2 窓に入るべき出力を列挙できない (settings.out_step_interval {oi!r}: {e})"
        return out
    have = {r["step"] for r in rows}
    out["required"] = {"n": len(need), "first": need[0], "last": need[-1], "out_step_interval": int(oi)}
    lack = [k for k in need if k not in have]
    if lack:
        out["why"] = f"2 窓に入るべき出力 {len(need)} のうち記録の行が無い出力 {len(lack)}: {lack[:8]}" + (" …" if len(lack) > 8 else "")
        out["missing_rows"] = lack
        return out
    by = {r["step"]: r for r in rows}
    wins = [[by[k] for k in need if lo <= k <= hi] for lo, hi in out["windows"]]
    out["n"] = [len(w) for w in wins]
    if min(out["n"]) < mp:
        out["why"] = f"窓の点数 {out['n']} (各 ≥ {mp} が要る)"
        return out
    bad = []
    for w in wins:
        for r in w:
            d = r.get("delta")
            if not isinstance(d, dict) or "delta_E" not in d:
                bad.append(f"{r['step']}: " + ("記録なし" if not isinstance(d, dict) else str(d.get("error") or d.get("skipped") or "δ なし")))
            elif len(d["delta_E"]) != len(xs) or d.get("n_missing") or d.get("n_held"):
                bad.append(f"{r['step']}: 列の数 {len(d['delta_E'])} / 欠測 {d.get('n_missing')} / held {d.get('n_held')}")
    if bad:
        out["why"] = "窓の出力の δ の記録が使えない: " + "; ".join(bad[:6]) + (" …" if len(bad) > 6 else "")
        return out
    mm = float(S) * 1e3
    worst, ok = {}, True
    for q in ("delta_E", "delta_target"):
        A = [np.array([r["delta"][q] for r in w], dtype=float) for w in wins]
        if not all(np.all(np.isfinite(a)) for a in A):
            out["why"] = f"{q} に非有限値"
            return out
        for name, v in (("range_W1", np.ptp(A[0], axis=0)), ("range_W2", np.ptp(A[1], axis=0)),
                        ("mean_diff", np.abs(A[1].mean(axis=0) - A[0].mean(axis=0)))):
            i = int(np.argmax(v))
            worst[f"{q}.{name}"] = {"max_mm": float(v[i] * mm), "x_at_max": float(xs[i])}
            ok &= bool(v[i] * mm <= tol_mm)
    out["worst"] = worst
    out["measured"] = True
    out["pass"] = ok
    if not ok:
        out["why"] = "; ".join(f"{k} {v['max_mm']:.4g} mm (x {v['x_at_max']:.4g})" for k, v in worst.items() if v["max_mm"] > tol_mm)
    return out


def coverage_check(summary: dict, delta_next, x_exit: float) -> dict:
    """抽出の全域の被覆と使用域の held・欠測 (§4.9): (i) NS のすべての列を抽出した (`not_extracted` が空: Euler の範囲の外・端の補間範囲の
    不足 (欠測)・帯の計算の失敗が 1 列も無い)、(ii) 使用域 [0, x_e] の delta_r_next.csv に held が無く δ_E・δ_target が有限。
    Euler の範囲は抽出器の `deltastar.euler_cover_range` (§5.1 #36g: 端の外側だけの幾何の許容、許容の内側で端を越えた列は端の列の値;
    その列は `edge_tol_cols` として記録に残す)。"""
    ne = (summary or {}).get("not_extracted")
    out = {"x_exit": float(x_exit), "pass": False}
    if not isinstance(ne, dict):
        out["why"] = "抽出の記録に not_extracted が無い (被覆を判定できない)"
        return out
    out["not_extracted"] = {"n": ne["n"], "n_cols_ns": ne.get("n_cols_ns"), "x": ne["x"][:20], "reason": ne["reason"][:20],
                            "euler_cover": ne.get("euler_cover"), "edge_tol_cols": ne.get("edge_tol_cols")}
    nx = read_delta_r_next(delta_next)
    x = np.asarray(nx["x_rt"], dtype=float)
    use = (x >= 0.0) & (x <= float(x_exit) + 1e-9)
    held = (np.asarray(nx["held"]) != 0) & use
    nonfin = use & ~(np.isfinite(nx["delta_E"]) & np.isfinite(nx["delta_target"]))
    out.update(n_use=int(use.sum()), n_held_use=int(held.sum()), n_nonfinite_use=int(nonfin.sum()),
               x_held_use=[float(v) for v in x[held][:20]])
    why = []
    if ne["n"]:
        why.append(f"抽出しなかった NS の列 {ne['n']} (理由 {sorted(set(ne['reason']))}、x {', '.join(f'{v:.4f}' for v in ne['x'][:6])})")
    if out["n_held_use"]:
        why.append(f"使用域の held {out['n_held_use']} 点")
    if out["n_nonfinite_use"]:
        why.append(f"使用域の δ_E・δ_target の非有限 {out['n_nonfinite_use']} 点")
    out["pass"] = not why
    if why:
        out["why"] = "; ".join(why)
    return out


def next_wall_check(problem, delta_next) -> dict:
    """次の壁の幾何の健全性: delta_r_next.csv の次の壁の入力 (第 2 列 δ_target) から `prepare_ns(delta_r_csv=...)` と同じ構築
    (`design_chain` → `delta_r_from_table` → `build_physical_wall`、offset radial) で物理壁を作り、既存の壁のゲート (上流のゲート: 正の半径・
    一意な最小・前後の単調性・継ぎ目、1 本の B-spline の作り直しの許容、必須属性、validate) を通るか。例外・validate の不合格は不合格。"""
    from ..evaluate.runner_axismach import build_physical_wall, delta_r_from_table, design_chain, load_problem, require_moc_gate
    from ..geometry.wall_axismach import check_required_attrs
    if problem is None:
        return {"pass": False, "why": "問題 YAML が渡されていない (次の壁を作れない)"}
    try:
        p = load_problem(problem)
        d = design_chain(p)
        require_moc_gate(p, d)
        tbl = np.loadtxt(delta_next, delimiter=",", skiprows=1)
        # prepare_ns と同じ: 表の端が設計の端から ε_f32 以内なら揃える (plan tooling-nozzle-card-local-mesh §5.1 #3h)
        xs, ends = snap_table_ends_same_design(tbl[:, 0], design_anchors(d["wall"]), float(p.spec["r_throat"]))
        W = build_physical_wall(p, d, float(p.spec["r_throat"]), delta_r_x=delta_r_from_table(xs, tbl[:, 1]), offset="radial")
        check_required_attrs(W)
        msgs = W.validate()
    except Exception as e:  # noqa: BLE001 — 壁のゲートの例外は不合格として記録する
        return {"pass": False, "why": f"次の壁の構築が例外 ({type(e).__name__}: {e})"[:1200]}
    g = getattr(W, "upstream_gate", None) or {}
    out = {"pass": not msgs, "wall": type(W).__name__, "pw_upstream": getattr(W, "pw_upstream", None),
           "throat": {"x": float(W.x_throat), "r": float(W.r_throat), "kappa": float(W.kappa_throat)},
           "upstream_gate_pass": g.get("pass"), "validate": msgs, "table_ends": ends}
    if msgs:
        out["why"] = "validate: " + "; ".join(msgs)
    return out


def _same_euler_ref(a, b) -> bool:
    """2 つの Euler の参照の指定が同じ (run のディレクトリ・固定の res と sha256)。"""
    from ..evaluate import stop_rule as SR
    pa, pb = SR.resolve_path(a), SR.resolve_path(b)
    if pa is None or pb is None or pa.resolve() != pb.resolve():
        return False
    fa = (a.get("res"), a.get("sha256")) if isinstance(a, dict) else (None, None)
    fb = (b.get("res"), b.get("sha256")) if isinstance(b, dict) else (None, None)
    return fa == fb


def override_reason(name: str, value) -> str | None:
    """受け渡しの関門の明示の上書き (`allow_not_reached`・`allow_unconverged`・`allow_unsettled`) の値を理由の文字列にする
    (2026-10-14 統合 integrate/design-flow-1 と同じ形、forge-f5 の要望: 上書きは理由の文字列が必須; §5.1 #36h)。None・False は上書きなし
    (None を返す; 2026-10-14 までの呼び出し元が渡す False をそのまま受ける)。空でない文字列はその理由 (前後の空白を除く)。True や空の文字列・
    文字列でない値は ValueError (理由の無い上書きを通さない)。"""
    if value is None or value is False:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} は上書きの理由の文字列で書く (受け取った値: {value!r}; 理由の無い上書きは通さない)")
    return value.strip()


def _caller_delta_e_source(extractor, delta_e_source) -> dict | None:
    """呼び出し側の δ_E の作り方 (§5.1 #36h): 関門の抽出を呼び出し側の `extractor` (`extract_and_merge` と同じ呼び方・戻り) に替えるときの
    作り方の記録。extractor と delta_e_source は一緒に渡す (片方だけは ValueError: 関門の抽出は最新の出力 1 枚なので、作り方の記録だけを
    書き替えない)。delta_e_source は mapping で kind (空でない文字列、例 "window_mean")・n_snapshots (1 以上の整数) が要る。
    5 枚の窓の平均の作り方 (forge-f5、plan tooling-nozzle-cooled-wall-ns-settling §4.4・§5.1 #3b) はこの口から渡す (ここには移植しない)。"""
    if extractor is None and delta_e_source is None:
        return None
    if extractor is None or delta_e_source is None:
        raise ValueError("extractor と delta_e_source は一緒に渡す (関門の抽出は最新の出力 1 枚; 作り方の記録と抽出を食い違わせない)")
    if not callable(extractor):
        raise ValueError(f"extractor は extract_and_merge と同じ呼び方の関数 (受け取った値: {extractor!r})")
    n = delta_e_source.get("n_snapshots") if isinstance(delta_e_source, dict) else None
    if (not isinstance(delta_e_source, dict) or not isinstance(delta_e_source.get("kind"), str) or not delta_e_source["kind"].strip()
            or not isinstance(n, int) or isinstance(n, bool) or n < 1):
        raise ValueError(f"delta_e_source は kind (空でない文字列) と n_snapshots (1 以上の整数) を持つ mapping (受け取った値: {delta_e_source!r})")
    return {**delta_e_source, "provider": "caller (handoff_gate の extractor)"}


def _single_snapshot_source(od: Path, stop_step) -> dict:
    """関門の既定の抽出 (`extract_and_merge`、NS の run の最新の出力 1 枚) の δ_E の作り方の記録 (§5.1 #36h)。読んだ出力は抽出器の記録
    (delta_r_equiv_diag.json の settings.ns_res) から取る (無ければ None)。"""
    try:
        res = json.loads((od / "delta_r_equiv_diag.json").read_text())["settings"]["ns_res"]
    except (OSError, ValueError, KeyError, TypeError):
        res = None
    try:
        step = int(Path(res).stem.split("_")[1]) if res else None
    except (IndexError, ValueError):
        step = None
    return {"kind": "single_snapshot", "n_snapshots": 1, "res": res, "steps": [step] if step is not None else None,
            "stop_step": stop_step, "is_stop_output": (step is not None and stop_step is not None and step == int(stop_step)),
            "provider": "handoff_gate (extract_and_merge、NS の run の最新の出力 1 枚)"}


def _unsettled_entry(name: str, c: dict, reason: str) -> dict:
    """allow_unsettled で通した条件の記録 (実測の値; §5.1 #36h)。保存則: r の窓の平均と振幅・閾値・窓。2 窓: 変動の最大とその量・窓・x、
    閾値を超えた量の全部。"""
    if name == "conservation":
        return {"override": "allow_unsettled", "condition": "conservation", "reason": reason, "r_mean": c.get("r_mean"), "r_amp": c.get("r_amp"),
                "r_mean_max": c.get("r_mean_max"), "r_amp_max": c.get("r_amp_max"),
                "window": [c["step"] - c["window"], c["step"]], "n": c.get("n")}
    w = c.get("worst") or {}
    k = max(w, key=lambda k_: w[k_]["max_mm"]) if w else None
    return {"override": "allow_unsettled", "condition": "two_windows", "reason": reason, "windows": c.get("windows"), "tol_mm": c.get("tol_mm"),
            "max_mm": (w[k]["max_mm"] if k else None), "quantity": k, "x_at_max": (w[k]["x_at_max"] if k else None),
            "exceeded": {k_: v for k_, v in w.items() if v["max_mm"] > float(c.get("tol_mm", HANDOFF_WINDOW_TOL_MM))}, "n": c.get("n")}


def handoff_gate(run_dir, *, problem=None, euler_ref=None, omega: float | None = None, allow_not_reached: str | None = None,
                 allow_unconverged: str | None = None, allow_unsettled: str | None = None, write: bool = True, extract: bool = True,
                 extractor=None, delta_e_source: dict | None = None) -> dict:
    """受け渡しの関門 (plan tooling-nozzle-isothermal-wall-chain §4.9・§5.1 #36a・#36e・#36h、codex diagnose 2026-10-14・diagnose 3、§6 V-h370):
    その run の場から抽出して次の壁の入力 (`delta_r_next.csv`・`fixed_point.json`) を作ってよいか (§4.9 の「更新に進む条件」)。
    (a) 本段の停止の状態 (`stop_state`) が REACHED。`allow_not_reached` で NOT_USED・CENSORED・UNAVAILABLE も通す (明示の上書き)。
        DIVERGED・EXEC_ERROR・DATA_ERROR・MISSING は通さない。
    (b) `convergence_verdict` (`check_convergence.py --segment --drop 3`) が PASS。`allow_unconverged` で NOT_CONVERGED も通す (明示の上書き)。
        DIVERGED・UNDETERMINED は通さない (判定不能は合格ではない)。
    (c) 保存則 (`conservation_check`、stop_rule.json の r の窓)。(e) δ_E・δ_target の隣り合う 2 窓 (`two_window_check`、見張りの中の出力ごとの
        記録)。`allow_unsettled` (§5.1 #36h) で実測が閾値を超えた不合格 (check の measured が真) だけ通す。判定不能 (記録の欠け・窓を覆わない・
        点数の不足・非有限・δ の記録の held と欠測) は通さない。
    以下は **上書きでは通さない** (#36e):
    (d) 局所 y₁⁺ (`wall_resolution_verdict` の PASS)。(f) Euler の参照が本段の記録 (停止規則の euler_ref) と同じ、ω が固定の値
        (`HANDOFF_OMEGA`) と同じ。
    上書き 3 つは**理由の文字列** (`override_reason`; True・空の文字列・文字列でない値は ValueError) で渡し、HANDOFF_GATE.json の override に
    理由と実際に使った上書き (applied) を残す。**落ち着きの条件** ((a)・(c)・(e)) を全部満たせば `settled: true`、上書きで通したら
    `settled: false` と `unsettled_overrides` (理由・実測の値: 停止の状態、r の平均と振幅と窓、2 窓の変動の最大とその量・窓・x) を残し、
    抽出の要約 (`delta_r_extract_summary.json`・`_after_ns` の `fixed_point.json`) にも同じ印を写す (「落ち着いていない場からの設計」を後で
    探せるように; `grep '"settled": false'`)。
    (a)〜(f) を満たした (上書きで通したものを含む) ときだけ抽出する (`extract`; run の HANDOFF_EXTRACT_DIR に、固定の条件
    `handoff_extract_settings` = edge・一定 x の断面への補間・端の欠測・ω・平滑化の λ・昇格の上限、Euler の参照は (f) の固定のスナップショット)
    → (g) 抽出の全域の被覆と使用域の held・欠測なし (`coverage_check`)、(h) 次の壁の幾何の健全性 (`next_wall_check`、problem が要る)。
    (g)(h) は上書きでは通さない。抽出の例外 (昇格の上限を超えた単調性ガード・sha の不一致ほか) も不合格。全部満たしたら抽出の表を run の
    直下に写す (`delta_r_next.csv` ほか; 不合格なら写さない)。
    δ_E の作り方 (`delta_e_source`): 既定は `extract_and_merge` で NS の run の最新の出力 1 枚 (`kind: single_snapshot`、読んだ res を記録)。
    呼び出し側が別の作り方 (例: 5 枚の窓の平均) で抽出するときは `extractor` (`extract_and_merge` と同じ呼び方・戻り、out_dir に
    delta_r_next.csv を書く) と `delta_e_source` (kind・n_snapshots を持つ mapping) を一緒に渡す (`_caller_delta_e_source`)。(g)(h) は
    その抽出にも同じに掛かる。
    戻り: {proceed, blocked_by, settled, unsettled_overrides, delta_e_source, stop_rule, convergence, checks, override, extract (settings・made・
    summary)}。write なら run に HANDOFF_GATE.json (原子的に)・HANDOFF_CONVERGENCE.txt・HANDOFF_WALL_RESOLUTION.txt を書く (抽出は write で
    なくても HANDOFF_EXTRACT_DIR に書く)。"""
    import datetime
    import shutil
    from ..evaluate import stop_rule as SR
    run_dir = Path(run_dir)
    marks = [m for m in DIAG_ONLY_MARKS if (run_dir / m).exists()]
    if marks:
        raise HandoffBlocked(f"{run_dir} は診断専用の run ({marks[0]}) なので受け渡しの関門に回さない (§6 V-h370-diag; 上書きでは通さない)")
    why_nr = override_reason("allow_not_reached", allow_not_reached)
    why_uc = override_reason("allow_unconverged", allow_unconverged)
    why_us = override_reason("allow_unsettled", allow_unsettled)
    src = _caller_delta_e_source(extractor, delta_e_source)
    if SR.criteria_present(run_dir):
        # 新方式 (settle_criteria.yaml のある run; plan tooling-settle-convergence-and-monitor §4.9 の 1): 停止規則 v2 の記録と照合し、
        # 停止の step で抽出済みの δ_E を使う (抽出し直さない)。判定不能・照合の不一致で旧方式へ落とさない
        if src is not None:
            raise ValueError("新方式 (settle_criteria.yaml のある run) の δ_E は停止規則 v2 が止めた step で抽出したものだけを使う — "
                             "extractor・delta_e_source は旧方式の run だけ")
        return _handoff_gate_v2(run_dir, problem=problem, euler_ref=euler_ref, omega=omega, why_nr=why_nr, why_uc=why_uc, why_us=why_us,
                                write=write, extract=extract)
    st = stop_state(run_dir)
    cv = convergence_verdict(run_dir)
    blocked, unsettled, applied = [], [], []
    s = st["status"]
    if s != "REACHED":
        if why_nr is not None and s in _STOP_OVERRIDABLE:
            unsettled.append({"override": "allow_not_reached", "condition": "stop_state", "reason": why_nr, "status": s,
                              "stop_step": st.get("stop_step")})
            applied.append("allow_not_reached")
        else:
            blocked.append(f"停止の状態 {s} (REACHED が要る"
                           + (f"; allow_not_reached で通せるのは {'・'.join(_STOP_OVERRIDABLE)} だけ)" if s not in _STOP_OVERRIDABLE
                              else "; 明示の上書き allow_not_reached が無い)"))
    k = cv["kind"]
    if k != "PASS":
        if why_uc is not None and k == "NOT_CONVERGED":
            applied.append("allow_unconverged")
        else:
            blocked.append(f"収束の判定 {k} (check_convergence.py --segment の終了コード {cv['rc']}"
                           + ("; 明示の上書き allow_unconverged が無い)" if k == "NOT_CONVERGED" else "; 上書きでも通さない)"))
    # --- (c)〜(f): 記録から ((c)(e) の実測の超過だけ allow_unsettled で通す、ほかは上書きでは通さない) ---
    try:
        info = json.loads((run_dir / "prepare_info.json").read_text())
    except (OSError, ValueError):
        info = {}
    sr = info.get("stop_rule") if isinstance(info.get("stop_rule"), dict) else {}
    settings = dict(sr.get("extract") or handoff_extract_settings())
    settings.pop("source", None)
    rec_sr = _read_stop_record(run_dir)
    step = st.get("stop_step")
    if step is None and rec_sr is not None:
        step = rec_sr.get("stop_step")
    checks = {"stop_step": step}
    checks["conservation"] = conservation_check(rec_sr, step)
    wr = wall_resolution_verdict(run_dir, step)
    checks["wall_resolution"] = {k_: v for k_, v in wr.items() if k_ != "output"}
    checks["wall_resolution"]["pass"] = wr["kind"] == "PASS"
    checks["two_windows"] = two_window_check(rec_sr, step)
    run_eref = sr.get("euler_ref") or info.get("euler_ref")
    use_eref = run_eref if run_eref else euler_ref
    eref_ok = use_eref is not None and (euler_ref is None or run_eref is None or _same_euler_ref(euler_ref, run_eref))
    checks["euler_ref"] = {"run": run_eref, "caller": (euler_ref if isinstance(euler_ref, dict) or euler_ref is None else str(euler_ref)),
                           "pinned": isinstance(use_eref, dict), "pass": bool(eref_ok)}
    om = HANDOFF_OMEGA if omega is None else float(omega)
    checks["omega"] = {"caller": om, "fixed": settings.get("omega"), "pass": bool(om == settings.get("omega"))}
    for name, label in (("conservation", "保存則 (質量の収支)"), ("wall_resolution", "局所 y₁⁺ (check_wall_resolution)"),
                        ("two_windows", "δ_E・δ_target の 2 窓"), ("euler_ref", "Euler の参照 (本段の記録と同じ)"), ("omega", "ω (固定の値)")):
        c = checks[name]
        if c["pass"]:
            continue
        if name in _UNSETTLED_OVERRIDABLE and why_us is not None and c.get("measured"):
            unsettled.append(_unsettled_entry(name, c, why_us))
            if "allow_unsettled" not in applied:
                applied.append("allow_unsettled")
            continue
        blocked.append(f"{label}: " + str(c.get("why") or c.get("kind") or c)
                       + (" (判定不能・記録の欠けは allow_unsettled でも通さない)" if name in _UNSETTLED_OVERRIDABLE and why_us is not None else ""))
    # 落ち着きの条件 ((a)・(c)・(e)) を上書きなしで満たしたか (§5.1 #36h)
    settled = bool(s == "REACHED" and checks["conservation"]["pass"] and checks["two_windows"]["pass"])
    mark = {"settled": settled, "unsettled_overrides": unsettled, "delta_e_source": src}
    # --- (g)(h): 抽出 (a)〜(f) を満たしたときだけ (上書き不可) ---
    ext = {"band_select": HANDOFF_BAND_SELECT, "settings": settings, "dir": HANDOFF_EXTRACT_DIR, "made": None, "summary": None}
    if blocked or not extract:
        why = "前の条件で止めた" if blocked else "extract=False"
        checks["coverage"] = {"pass": False, "evaluated": False, "why": why}
        checks["next_wall"] = {"pass": False, "evaluated": False, "why": why}
        if not blocked:
            blocked.append("抽出をしていない (extract=False) — 被覆・held・次の壁を判定していない")
    else:
        od = run_dir / HANDOFF_EXTRACT_DIR
        if od.exists():
            shutil.rmtree(od)
        od.mkdir(parents=True)
        try:
            summ = (extractor or extract_and_merge)(run_dir, SR.euler_ref_spec(use_eref) or use_eref, out_dir=od, **settings)
            if not isinstance(summ, dict):
                raise TypeError(f"抽出の戻りが要約の mapping でない ({type(summ).__name__})")
            ext["summary"] = summ
        except Exception as e:  # noqa: BLE001 — 抽出の例外 (昇格の上限・sha の不一致ほか) は関門の不合格
            summ = None
            blocked.append(f"抽出が例外 ({type(e).__name__}: {e})"[:600])
            checks["coverage"] = {"pass": False, "evaluated": False, "why": "抽出が例外"}
            checks["next_wall"] = {"pass": False, "evaluated": False, "why": "抽出が例外"}
        if summ is not None:
            if src is None:
                mark["delta_e_source"] = _single_snapshot_source(od, step)
            # 抽出の要約に落ち着きの印を写す (表の隣の記録; run の直下へ写すときもこの要約を写す)
            summ.update(mark)
            (od / "delta_r_extract_summary.json").write_text(json.dumps(summ, indent=1, default=str))
            dn = od / "delta_r_next.csv"
            x_e = ((info.get("physical_wall") or {}).get("domain") or [None, None])[1]
            if x_e is None:
                x_e = float(np.max(np.loadtxt(dn, delimiter=",", skiprows=1)[:, 0]))
            checks["coverage"] = coverage_check(summ, dn, float(x_e))
            checks["next_wall"] = next_wall_check(problem, dn)
            for name, label in (("coverage", "抽出の全域の被覆・使用域の held と欠測"), ("next_wall", "次の壁の幾何の健全性")):
                if not checks[name]["pass"]:
                    blocked.append(f"{label}: {checks[name].get('why')}")
            if not blocked:
                made = []
                for f in ("delta_r_next.csv", "delta_r_equiv.csv", "delta_r_equiv_diag.json", "delta_r_extract_summary.json", "delta_r_edge.csv"):
                    if (od / f).exists():
                        shutil.copy2(od / f, run_dir / f)
                        made.append(f)
                ext["made"] = made
    rec = {"version": 3, "method": "v1",
           "method_note": "旧方式 (settle_criteria.yaml の無い run): 停止規則 v1・check_convergence.py --segment・stop_rule.json の窓・関門の抽出",
           "plan": "plans/active/tooling-nozzle-isothermal-wall-chain.md §4.9・§5.1 #36a・#36e・#36h (§6 V-h370)",
           "run": str(run_dir), "time": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
           "proceed": not blocked, "blocked_by": blocked,
           "settled": mark["settled"], "unsettled_overrides": mark["unsettled_overrides"], "delta_e_source": mark["delta_e_source"],
           "stop_rule": st, "convergence": {k_: v for k_, v in cv.items() if k_ != "output"}, "checks": checks,
           "override": {"allow_not_reached": why_nr is not None, "allow_unconverged": why_uc is not None, "allow_unsettled": why_us is not None,
                        "allow_not_reached_reason": why_nr, "allow_unconverged_reason": why_uc, "allow_unsettled_reason": why_us,
                        "applied": applied,
                        "scope": "上書きは (a) 停止の状態 (allow_not_reached)・(b) 収束の判定 (allow_unconverged)・保存則と 2 窓の実測の超過 "
                                 "(allow_unsettled) だけ。保存則・2 窓の判定不能 (記録の欠け・覆わない・点数・非有限・held・欠測)・y₁⁺・参照・ω・"
                                 "被覆・held・次の壁は上書きでは通さない"},
           "extract": ext,
           "note": "proceed は §4.9 の「更新に進む条件」の全部を満たしたか (壁の更新は次の pass へ自動でつながない)。0.05 mm の固定点の条件は設計を終える用で、ここでは見ない。"
                   "settled: false は落ち着きの条件を上書きで通した「落ち着いていない場からの設計」 (後で作り直す対象)"}
    if write:
        (run_dir / HANDOFF_CONVERGENCE).write_text(cv["output"])
        (run_dir / HANDOFF_WALL_RESOLUTION).write_text(wr["output"])
        _write_json_atomic(run_dir / HANDOFF_GATE, rec)
    return rec


# --- 受け渡しの関門の新方式 (settle_criteria.yaml のある run; plan tooling-settle-convergence-and-monitor §4.9 の 1・2、
#     codex diagnose 2026-10-16 settle-migration-7 の M1・M2・M5) ---------------------------------------------------------------------
HANDOFF_V2_PLAN = ("plans/active/tooling-settle-convergence-and-monitor.md §4.9 の 1・2 (§4.8 の 2、codex diagnose 2026-10-16 "
                   "settle-migration-7: 照合・δ_E の契約・分岐の規則)")
_V2_TABLE_FILES = ("delta_r_next.csv", "delta_r_equiv.csv", "delta_r_equiv_diag.json", "delta_r_extract_summary.json", "delta_r_edge.csv")
_V2_SENS_FILES = ("band_sensitivity.json", "band_sensitivity.csv")
_V2_NEED_FILES = ("delta_r_next.csv", "delta_r_equiv_diag.json", "delta_r_extract_summary.json")


def _settle_group(qs: list, kinds, label: str) -> dict:
    """再評価 (SETTLE_VERDICT.json) の量のうち kind が kinds のものの判定。全部 SETTLED → kind PASS。
    不合格の分類: DIVERGED (どれかが DIVERGED) / UNDETERMINED (どれかが INDETERMINATE か、理由に閾値の超過以外がある・量が無い) /
    NOT_SETTLED (閾値の超過だけ = measured; 上書きの対象)。"""
    sel = [q for q in qs if isinstance(q, dict) and q.get("kind") in kinds]
    out = {"label": label, "quantities": {str(q.get("name")): q.get("verdict") for q in sel}, "pass": False}
    if not sel:
        out.update(kind="UNDETERMINED", why=f"再評価に {label} の量 ({'・'.join(kinds)}) が無い")
        return out
    bad = [q for q in sel if q.get("verdict") != "SETTLED"]
    if not bad:
        out.update(kind="PASS", **{"pass": True})
        return out
    vs = {q.get("verdict") for q in bad}
    if "DIVERGED" in vs:
        out["kind"] = "DIVERGED"
    elif "INDETERMINATE" in vs or vs - {"NOT_SETTLED"} or any(set(q.get("categories") or []) - {"threshold"} for q in bad):
        out["kind"] = "UNDETERMINED"
    else:
        out["kind"] = "NOT_SETTLED"
        out["measured"] = True
    out["failed"] = {str(q.get("name")): {"verdict": q.get("verdict"), "categories": q.get("categories"),
                                          "checks": [c.get("check") for c in q.get("checks") or [] if not c.get("ok")][:6],
                                          "reasons": [r.get("text") for r in q.get("reasons") or []][:3]} for q in bad}
    out["why"] = "; ".join(f"{q.get('name')}: {q.get('verdict')}" + (f" ({(q.get('reasons') or [{}])[0].get('text')})" if q.get("reasons") else "")
                           for q in bad)[:800]
    return out


def _v2_binding(run_dir: Path, st: dict | None, S) -> dict:
    """照合 (§4.9 の 1、codex M1): STOP_RULE_STATE・SETTLE_VERDICT (停止の step の再評価)・STOPPED_AT・停止の step・本段の起動
    (MAIN_STAGE_STARTED の main_config の sha256・fnv と最後の forge の起動)・今の settle_criteria.yaml の sha256・stop_rule.json が
    互いに一致すること。どれかが照合できない・食い違えば problems (受け渡しの関門は上書きでも通さない)。
    戻り: {pass, problems, items, verdict (SETTLE_VERDICT.json の中身か None), stop_step, criteria_sha256}。"""
    from ..evaluate import stop_rule as SR
    probs, items = [], {}
    sel = S.select_judge(run_dir)
    items["select"] = {k: sel.get(k) for k in ("tool", "mode", "why", "problems")}
    b = sel.get("binding") or {}
    items["main_stage"] = {k: b.get(k) for k in ("record", "main_config_fnv", "main_config_sha256", "config_fnv", "config_sha256",
                                                 "launch_cfg_fnv", "record_criteria_sha256", "criteria_sha256")}
    if sel["mode"] != "main":
        probs.append(f"本段の起動として照合できない ({sel['why']} [mode {sel['mode']}])")
    probs += [f"本段の起動: {p}" for p in sel["problems"]]
    if sel["mode"] == "main" and b.get("launch_cfg_fnv") is None:
        probs.append(f"{S.LAUNCHES_FILE} の起動の記録が無い (最後の forge の起動を本段と照合できない)")
    crit_sha = b.get("criteria_sha256")
    if crit_sha is None:
        try:
            lc = S.load_criteria(run_dir)
            crit_sha = lc.get("sha256")
        except S.CriteriaError as e:
            probs.append(f"{S.CRITERIA_FILE} の誤り: {e}")
    out = {"pass": False, "problems": probs, "items": items, "verdict": None, "stop_step": None, "criteria_sha256": crit_sha}
    if not isinstance(st, dict):
        probs.append(f"{SR.STATE} が無い (停止規則 v2 の記録が要る)")
        return out
    if st.get("method") != "v2":
        probs.append(f"{SR.STATE} の method が v2 でない ({st.get('method')!r})")
    sc = (st.get("criteria") or {}).get("sha256")
    if sc != crit_sha:
        probs.append(f"{SR.STATE} の基準の sha256 {str(sc)[:16]}… が今の {S.CRITERIA_FILE} ({str(crit_sha)[:16]}…) と違う")
    try:
        note = json.loads((run_dir / S.MAIN_STAGE_RECORD).read_text())
        n_main = ((note.get("main_config") or {}).get("nStepOuter")) if isinstance(note, dict) else None
    except (OSError, ValueError):
        n_main = None
    if st.get("cap") is None or n_main != st.get("cap"):
        probs.append(f"本段の記録の nStepOuter {n_main!r} が停止規則の上限 {st.get('cap')!r} と違う")
    k = st.get("stop_step")
    items["stop_step"] = k
    if not isinstance(k, int) or isinstance(k, bool):
        probs.append(f"{SR.STATE} に停止の step が無い ({k!r}; 状態 {st.get('status')})")
        return out
    out["stop_step"] = k
    # STOPPED_AT (REACHED は協調停止の通知、CENSORED は上限まで回って通知なし)
    stopped = S.read_stopped_at(run_dir)
    sa_rec = st.get("stopped_at") or {}
    sa_path = run_dir / SR.STOPPED_AT
    items["stopped_at"] = {"present": stopped is not None, "step": (stopped or {}).get("step"), "state": sa_rec}
    if st.get("status") == "CENSORED":
        if stopped is not None:
            probs.append(f"CENSORED なのに {SR.STOPPED_AT} がある")
        if k != st.get("cap"):
            probs.append(f"CENSORED の停止の step {k} が上限 {st.get('cap')} と違う")
    else:
        if stopped is None or "error" in stopped:
            probs.append(f"{SR.STOPPED_AT} が無い・読めない ({(stopped or {}).get('error', '無い')})")
        else:
            if stopped["step"] != k:
                probs.append(f"{SR.STOPPED_AT} の step {stopped['step']} が {SR.STATE} の停止の step {k} と違う")
            if sa_rec.get("step") != stopped["step"] or sa_rec.get("sha256") != _sha256(sa_path):
                probs.append(f"{SR.STOPPED_AT} (step {stopped['step']}、sha256 {_sha256(sa_path)[:16]}…) が {SR.STATE} の記録 "
                             f"(step {sa_rec.get('step')}、sha256 {str(sa_rec.get('sha256'))[:16]}…) と違う")
    # 停止の場
    res = run_dir / f"res_{k}.h5"
    res_sha = _sha256(res) if res.is_file() else None
    items["stop_field"] = {"res": res.name, "sha256": res_sha, "state": st.get("stop_sha256")}
    if res_sha is None:
        probs.append(f"停止の場 {res.name} が無い")
    elif res_sha != st.get("stop_sha256"):
        probs.append(f"停止の場 {res.name} の sha256 {res_sha[:16]}… が {SR.STATE} の記録 {str(st.get('stop_sha256'))[:16]}… と違う")
    # 見張りの記録 stop_rule.json (δ_E の抽出の設定の出所)
    rp = run_dir / SR.RECORD
    rec_sha = _sha256(rp) if rp.is_file() else None
    if rec_sha is None or rec_sha != st.get("record_sha256"):
        probs.append(f"{SR.RECORD} (sha256 {str(rec_sha)[:16]}…) が {SR.STATE} の記録 ({str(st.get('record_sha256'))[:16]}…) と違う・無い")
    # 停止の step の再評価 SETTLE_VERDICT.json
    ve = st.get("reeval") or {}
    vp = run_dir / S.VERDICT_FILE
    try:
        v = json.loads(vp.read_text())
        if not isinstance(v, dict):
            raise ValueError("JSON の対応表でない")
    except (OSError, ValueError) as e:
        probs.append(f"{S.VERDICT_FILE} を読めない ({e})")
        return {**out, "pass": False}
    out["verdict"] = v
    dig = SR.quantities_digest(v)
    items["reeval"] = {"file": S.VERDICT_FILE, "end_step": v.get("end_step"), "verdict": v.get("verdict"), "digest": dig,
                       "state": {kk: ve.get(kk) for kk in ("end_step", "verdict", "digest", "sha256")},
                       "file_sha256": _sha256(vp), "file_sha256_matches_state": _sha256(vp) == ve.get("sha256")}
    if ve.get("error") or not ve:
        probs.append(f"{SR.STATE} に再評価の記録が無い ({ve.get('error') if ve else '無い'})")
    if v.get("end_step") != k or ve.get("end_step") != k:
        probs.append(f"再評価の終端 ({S.VERDICT_FILE} {v.get('end_step')!r}・{SR.STATE} {ve.get('end_step')!r}) が停止の step {k} と違う")
    if (v.get("criteria") or {}).get("sha256") != crit_sha:
        probs.append(f"{S.VERDICT_FILE} の基準の sha256 が今の {S.CRITERIA_FILE} と違う")
    if dig != ve.get("digest") or v.get("verdict") != ve.get("verdict"):
        probs.append(f"{S.VERDICT_FILE} の判定 ({v.get('verdict')}、量の sha256 {dig[:16]}…) が {SR.STATE} の再評価 "
                     f"({ve.get('verdict')}、{str(ve.get('digest'))[:16]}…) と違う")
    vsa = v.get("stopped_at")
    if st.get("status") == "CENSORED":
        if vsa is not None:
            probs.append(f"CENSORED の再評価に STOPPED_AT がある")
    elif not isinstance(vsa, dict) or vsa.get("step") != k:
        probs.append(f"再評価が読んだ STOPPED_AT の step ({(vsa or {}).get('step') if isinstance(vsa, dict) else vsa!r}) が停止の step {k} と違う")
    ll = S.last_launch(run_dir)
    vl = (v.get("launch") or {}).get("forge_launches_last")
    if ll is None or "error" in (ll or {}) or vl != ll:
        probs.append(f"再評価の起動の記録 (forge_launches.jsonl の最後の行) が今の最後の起動と違う・無い")
    out["pass"] = not probs
    return out


def _v2_delta_e(run_dir: Path, st: dict | None, rec_sr: dict | None, k, run_eref, om: float) -> dict:
    """δ_E の受け渡しの契約 (§4.9 の 1、codex M2): REACHED は δ_E が使える保証ではないので、停止規則 v2 が止めた step で抽出した
    `stop_rule_v2/delta_E_<k>/` を照合する: 状態の delta_E.status ok・抽出の step = 停止の step・停止の場の sha256・delta_r_next.csv の
    sha256 (状態・見張りの記録・ファイル)・抽出の設定 (handoff_extract_settings と同じ: 見張りの記録の repr、抽出器の diag の帯の選び方と
    一定 x の断面・端、要約の ω・昇格の上限・ノット)・Euler の参照 (抽出器の diag と本段の記録)・ω (呼び出し側 = 固定の値 = 抽出)。
    帯の感度は合否に使わず、記録の欠け・非有限は warnings に書く (「問題なし」とは扱わない)。
    戻り: {pass, why, problems, warnings, dir, summary (delta_r_extract_summary.json), delta_r_next_sha256, band_sensitivity}。"""
    from ..evaluate import stop_rule as SR
    probs, warns = [], []
    out = {"pass": False, "problems": probs, "warnings": warns, "dir": None, "summary": None}
    de = (st or {}).get("delta_E")
    if not isinstance(de, dict):
        probs.append(f"{SR.STATE} に δ_E の抽出の記録が無い ({de!r}; 停止規則は REACHED・CENSORED のときだけ抽出する)")
        out["why"] = "; ".join(probs)
        return out
    out["state"] = {kk: de.get(kk) for kk in ("status", "step", "dir", "res_sha256", "delta_r_next_sha256", "seconds")}
    if de.get("status") != "ok":
        probs.append(f"δ_E の抽出の状態が ok でない ({de.get('status')!r}: {de.get('error') or de.get('reason')})")
    if not isinstance(k, int) or de.get("step") != k:
        probs.append(f"δ_E の抽出の step {de.get('step')!r} が停止の step {k!r} と違う")
    if de.get("res_sha256") is None or de.get("res_sha256") != (st or {}).get("stop_sha256"):
        probs.append(f"δ_E の抽出元の場の sha256 {str(de.get('res_sha256'))[:16]}… が停止の場の記録 {str((st or {}).get('stop_sha256'))[:16]}… と違う")
    if probs:
        out["why"] = "; ".join(probs)
        return out
    dname = f"delta_E_{k}"
    dpath = Path(str(de.get("dir") or ""))
    if dpath.name != dname or dpath.parent.name != SR.V2_DIR:
        probs.append(f"δ_E の抽出の置き場 {de.get('dir')!r} が {SR.V2_DIR}/{dname} でない")
    od = run_dir / SR.V2_DIR / dname
    out["dir"] = str(od)
    miss = [f for f in _V2_NEED_FILES if not (od / f).is_file()]
    if miss:
        probs.append(f"δ_E の抽出の置き場 {od} に {miss} が無い")
        out["why"] = "; ".join(probs)
        return out
    dn_sha = _sha256(od / "delta_r_next.csv")
    out["delta_r_next_sha256"] = dn_sha
    rde = (rec_sr or {}).get("delta_E") if isinstance(rec_sr, dict) else None
    if dn_sha != de.get("delta_r_next_sha256") or not isinstance(rde, dict) or dn_sha != rde.get("delta_r_next_sha256"):
        probs.append(f"{dname}/delta_r_next.csv の sha256 {dn_sha[:16]}… が {SR.STATE}・{SR.RECORD} の記録と違う")
    if not isinstance(rde, dict) or rde.get("status") != "ok" or rde.get("step") != k or rde.get("res_sha256") != de.get("res_sha256"):
        probs.append(f"{SR.RECORD} の δ_E の記録が {SR.STATE} と合わない (状態・step・場の sha256)")
    # 抽出の設定 (見張りの記録は repr で残す)
    want = handoff_extract_settings()
    ex = (rde or {}).get("extract") or {}
    bad_ex = [kk for kk, vv in want.items() if ex.get(kk) != repr(vv)]
    if ex.get("function") != "feedback.deltastar_loop.extract_and_merge" or bad_ex:
        probs.append(f"δ_E の抽出の設定が受け渡しの固定の条件と違う (関数 {ex.get('function')!r}、違うキー {bad_ex})")
    try:
        summ = json.loads((od / "delta_r_extract_summary.json").read_text())
        diag = json.loads((od / "delta_r_equiv_diag.json").read_text())
        dset = diag.get("settings") or {}
    except (OSError, ValueError) as e:
        probs.append(f"δ_E の抽出の記録を読めない ({e})")
        out["why"] = "; ".join(probs)
        return out
    out["summary"] = summ
    if not (dset.get("band_select") == want["band_select"] and dset.get("remap_constant_x") is True
            and dset.get("remap_edge") == want["remap_edge"]):
        probs.append(f"抽出器の記録の帯・断面の条件が固定の条件と違う (band_select {dset.get('band_select')!r}・remap_constant_x "
                     f"{dset.get('remap_constant_x')!r}・remap_edge {dset.get('remap_edge')!r})")
    sm = summ.get("smooth") if isinstance(summ.get("smooth"), dict) else {}
    if summ.get("omega") != want["omega"] or summ.get("max_lam_factor") != want["max_lam_factor"] or sm.get("knot_spacing") != want["knot_spacing"]:
        probs.append(f"抽出の要約の ω {summ.get('omega')!r}・昇格の上限 {summ.get('max_lam_factor')!r}・ノット {sm.get('knot_spacing')!r} が"
                     f"固定の条件 ({want['omega']}・{want['max_lam_factor']}・{want['knot_spacing']}) と違う")
    if om != want["omega"] or om != summ.get("omega"):
        probs.append(f"ω: 呼び出し側 {om} と抽出 {summ.get('omega')!r} と固定 {want['omega']} が揃わない")
    # Euler の参照 (抽出器が読んだもの = 本段の記録の参照)
    pinned = dset.get("euler_ref_pinned")
    if isinstance(run_eref, dict):
        ok_e = (isinstance(pinned, dict) and pinned.get("res") == run_eref.get("res") and pinned.get("sha256") == run_eref.get("sha256")
                and _same_euler_ref(pinned, run_eref))
    else:
        rp_ = SR.resolve_path(run_eref) if run_eref else None
        ok_e = (pinned is None and rp_ is not None and dset.get("euler_run") is not None
                and Path(str(dset.get("euler_run"))).resolve() == rp_.resolve())
    out["euler_ref"] = {"extraction": pinned if pinned is not None else dset.get("euler_run"), "run": run_eref, "pass": bool(ok_e)}
    if not ok_e:
        probs.append(f"δ_E の抽出の Euler の参照 ({out['euler_ref']['extraction']!r}) が本段の記録 ({run_eref!r}) と違う")
    # 帯の感度 (記録だけ; 欠け・非有限は警告)
    bs = de.get("band_sensitivity")
    bsf = None
    if (od / "band_sensitivity.json").is_file():
        try:
            bsf = json.loads((od / "band_sensitivity.json").read_text())
        except ValueError:
            bsf = None
    out["band_sensitivity"] = bs
    if not isinstance(bs, dict) or bsf is None:
        warns.append("帯の位置の感度の記録が欠けている (状態の band_sensitivity・band_sensitivity.json; 「感度に問題なし」とは扱わない)")
    else:
        for kk in ("max_abs_B_minus_A_mm", "A_reproduces_delta_E_max_abs_rt"):
            vv = bs.get(kk)
            if not isinstance(vv, (int, float)) or not np.isfinite(vv):
                warns.append(f"帯の位置の感度の {kk} が無い・非有限 ({vv!r}; 「感度に問題なし」とは扱わない)")
        if bs.get("n_nonfinite_B_use"):
            warns.append(f"帯の位置の感度: 使用域で B (y_b × {bs.get('factor')}) が非有限の列 {bs.get('n_nonfinite_B_use')}")
    out["pass"] = not probs
    if probs:
        out["why"] = "; ".join(probs)
    return out


def _v2_convergence_text(k, v: dict | None, groups: dict, crit_sha) -> str:
    """HANDOFF_CONVERGENCE.txt (新方式): 停止の step の再評価 (SETTLE_VERDICT.json) の量ごとの結果。"""
    lines = [f"=== 受け渡しの関門 (新方式): 停止の step {k} の再評価 SETTLE_VERDICT.json (settle.judge、基準 sha256 {str(crit_sha)[:16]}…) ==="]
    if v is None:
        lines.append("再評価を読めない")
        return "\n".join(lines) + "\n"
    lines.append(f"OVERALL: {v.get('verdict')} (終端 {v.get('terminal_step')}、停止の観測 {v.get('stop_step')}、end_step {v.get('end_step')})")
    for q in v.get("quantities") or []:
        lines.append(f"  [{q.get('verdict')}] {q.get('name')} ({q.get('kind')})")
        for c in q.get("checks") or []:
            if not c.get("ok"):
                lines.append(f"      NG {c.get('check')}: {c.get('value')} (許容 {c.get('limit')})")
    for name, g in groups.items():
        lines.append(f"{g['label']}: {g['kind']}" + (f" — {g.get('why')}" if g.get("why") else ""))
    return "\n".join(lines) + "\n"


def _handoff_gate_v2(run_dir: Path, *, problem, euler_ref, omega, why_nr, why_uc, why_us, write: bool, extract: bool) -> dict:
    """受け渡しの関門の新方式 (`handoff_gate` から、settle_criteria.yaml のある run; plan tooling-settle-convergence-and-monitor §4.9 の 1)。
    (a) STOP_RULE_STATE.json (method v2) が REACHED。CENSORED・NOT_USED は allow_not_reached で通す (今どおり理由つき)、FAILED・MISSING は通さない。
    照合 (`_v2_binding`): 状態・再評価・STOPPED_AT・停止の step・本段の起動・今の基準・見張りの記録の一致。**上書きでも通さない**。
    (b)(c) 停止の step の再評価 (SETTLE_VERDICT.json) の**全体が SETTLED**: (b) 残差 (residual)・(c) 流量の収支 (mass_balance)・目的の量
    (出口コア M・線・流れ角・壁の熱流)。閾値の超過だけの NOT_SETTLED は (b) allow_unconverged・(c) と目的の量 allow_unsettled で通す
    (今の上書きの対応: 収束の判定 / 落ち着きの条件)。DIVERGED・判定不能は通さない。基準に残差が無く非有限の見張り (residual_finite) が
    ある run (wt_nozzle_v3; plan §4.4・§5.1 #8) は (b) を対象外 (止めない) にし、見張り (b′) の DIVERGED・判定不能は上書きでも通さない。
    (d) 局所 y₁⁺ (今のまま、停止の step の壁のダンプ)。(e) δ_E の契約 (`_v2_delta_e`; 抽出し直さない)。(f) Euler の参照・ω (今のまま + 抽出との一致)。
    (a)〜(f) を満たした (上書きを含む) ときだけ、停止規則の抽出の表を HANDOFF_EXTRACT_DIR に写して (g) 被覆・(h) 次の壁を今と同じ関数で
    判定し、全部満たしたら run の直下に旧方式と同じ名前で写す (`delta_r_next.csv` ほか; 次の pass がそのまま読む)。
    落ち着きの印 settled は (a) REACHED・(c)・目的の量が上書きなしで SETTLED (旧方式の (a)・(c)・(e) に対応)。"""
    import datetime
    import shutil
    from ..evaluate import stop_rule as SR
    S = SR._settle()
    st_full = SR.read_state(run_dir)
    st = stop_state(run_dir)
    blocked, unsettled, applied, warnings = [], [], [], []
    s = st["status"]
    # (a) 停止の状態
    if s != "REACHED":
        if why_nr is not None and s in ("NOT_USED", "CENSORED"):
            unsettled.append({"override": "allow_not_reached", "condition": "stop_state", "reason": why_nr, "status": s,
                              "stop_step": st.get("stop_step")})
            applied.append("allow_not_reached")
        else:
            blocked.append(f"停止の状態 {s} (REACHED が要る"
                           + ("; allow_not_reached で通せるのは NOT_USED・CENSORED だけ)" if s not in ("NOT_USED", "CENSORED")
                              else "; 明示の上書き allow_not_reached が無い)"))
    # 照合 (上書きでは通さない)
    bind = _v2_binding(run_dir, st_full, S)
    if not bind["pass"]:
        blocked.append("照合 (上書きでは通さない): " + "; ".join(bind["problems"])[:1500])
    k = bind["stop_step"] if bind["stop_step"] is not None else st.get("stop_step")
    checks = {"stop_step": k, "binding": {kk: vv for kk, vv in bind.items() if kk != "verdict"}}
    # (b)(c) 再評価の全体
    v = bind["verdict"]
    qs = (v or {}).get("quantities") or []
    kinds_v = {q.get("kind") for q in qs if isinstance(q, dict)}
    watch_kind = getattr(S, "KIND_RESIDUAL_FINITE", "residual_finite")
    if S.KIND_RESIDUAL not in kinds_v and watch_kind in kinds_v:
        # 残差を止める判定から外した基準 (wt_nozzle_v3; plan §4.4・§5.1 #8): (b) は対象外で止めない。発散は非有限の見張りで見る (下)
        g_res = {"label": "(b) 残差", "kind": "NOT_APPLICABLE", "pass": True, "quantities": {},
                 "why": "基準に残差の量が無い (残差を止める判定から外したプリセット、例 wt_nozzle_v3) — 対象外"}
    else:
        g_res = _settle_group(qs, (S.KIND_RESIDUAL,), "(b) 残差")
    g_watch = _settle_group(qs, (watch_kind,), "(b′) 残差の非有限の見張り") if watch_kind in kinds_v else None
    g_mass = _settle_group(qs, (S.KIND_MASS,), "(c) 流量の収支")
    g_tgt = _settle_group(qs, S.TARGET_KINDS, "目的の量")
    other = [q.get("name") for q in qs if q.get("kind") not in (S.KIND_RESIDUAL, S.KIND_MASS, watch_kind) + tuple(S.TARGET_KINDS)]
    overall = (v or {}).get("verdict")
    checks["reeval"] = {"verdict": overall, "by_quantity": {str(q.get("name")): q.get("verdict") for q in qs},
                        "pass": bool(overall == "SETTLED" and g_res["pass"] and (g_watch is None or g_watch["pass"])
                                     and g_mass["pass"] and g_tgt["pass"] and not other)}
    if other:
        blocked.append(f"再評価に知らない種類の量 {other} (上書きでは通さない)")
    if g_watch is not None:
        checks["residual_watch"] = g_watch
        if not g_watch["pass"]:
            blocked.append(f"(b′) 残差の非有限の見張り {g_watch['kind']}: {g_watch.get('why')} (上書きでも通さない)")
    cv_kind = {"PASS": "PASS", "NOT_SETTLED": "NOT_CONVERGED"}.get(g_res["kind"], g_res["kind"])
    cv = {"tool": f"{S.VERDICT_FILE} (停止の step の再評価、settle.judge; 新方式)", "kind": cv_kind, "rc": None,
          "verdict": [f"OVERALL: {overall}", f"残差: {g_res['kind']}"] + ([f"残差の非有限の見張り: {g_watch['kind']}"] if g_watch else []),
          "end_step": (v or {}).get("end_step"),
          "criteria_sha256": bind["criteria_sha256"]}
    checks["residual"] = g_res
    if not g_res["pass"]:          # 対象外 (NOT_APPLICABLE、v3) は止めない
        if why_uc is not None and g_res["kind"] == "NOT_SETTLED":
            applied.append("allow_unconverged")
        else:
            blocked.append(f"(b) 残差の再評価 {g_res['kind']}: {g_res.get('why')}"
                           + ("; 明示の上書き allow_unconverged が無い" if g_res["kind"] == "NOT_SETTLED" else "; 上書きでも通さない"))
    for name, g in (("mass_balance", g_mass), ("targets", g_tgt)):
        checks[name] = g
        if g["pass"]:
            continue
        if why_us is not None and g.get("measured"):
            unsettled.append({"override": "allow_unsettled", "condition": name, "reason": why_us, "end_step": k,
                              "quantities": g["quantities"], "failed": g.get("failed")})
            if "allow_unsettled" not in applied:
                applied.append("allow_unsettled")
            continue
        blocked.append(f"{g['label']}の再評価 {g['kind']}: {g.get('why')}"
                       + (" (判定不能・非有限は allow_unsettled でも通さない)" if why_us is not None else ""))
    # (d) 局所 y₁⁺
    wr = wall_resolution_verdict(run_dir, k)
    checks["wall_resolution"] = {k_: v_ for k_, v_ in wr.items() if k_ != "output"}
    checks["wall_resolution"]["pass"] = wr["kind"] == "PASS"
    if not checks["wall_resolution"]["pass"]:
        blocked.append("局所 y₁⁺ (check_wall_resolution): " + str(wr.get("kind")))
    # (e)(f) δ_E の契約・Euler の参照・ω
    try:
        info = json.loads((run_dir / "prepare_info.json").read_text())
    except (OSError, ValueError):
        info = {}
    sr = info.get("stop_rule") if isinstance(info.get("stop_rule"), dict) else {}
    run_eref = sr.get("euler_ref") or info.get("euler_ref")
    use_eref = run_eref if run_eref else euler_ref
    eref_ok = use_eref is not None and (euler_ref is None or run_eref is None or _same_euler_ref(euler_ref, run_eref))
    checks["euler_ref"] = {"run": run_eref, "caller": (euler_ref if isinstance(euler_ref, dict) or euler_ref is None else str(euler_ref)),
                           "pinned": isinstance(use_eref, dict), "pass": bool(eref_ok)}
    om = HANDOFF_OMEGA if omega is None else float(omega)
    rec_sr = _read_stop_record(run_dir)
    de = _v2_delta_e(run_dir, st_full, rec_sr, k, run_eref, om)
    warnings += de["warnings"]
    checks["delta_E"] = {kk: vv for kk, vv in de.items() if kk != "summary"}
    # ω: 呼び出し側 = 固定の値、抽出の要約があればその ω も同じ (要約が無いときは δ_E の契約の側で止まる)
    om_ex = (de.get("summary") or {}).get("omega") if de.get("summary") is not None else None
    checks["omega"] = {"caller": om, "fixed": HANDOFF_OMEGA, "extraction": om_ex,
                       "pass": bool(om == HANDOFF_OMEGA and (de.get("summary") is None or om_ex == HANDOFF_OMEGA))}
    for name, label in (("euler_ref", "Euler の参照 (本段の記録と同じ)"), ("omega", "ω (固定の値・抽出と同じ)")):
        if not checks[name]["pass"]:
            blocked.append(f"{label}: {checks[name]}"[:600])
    if not de["pass"]:
        blocked.append("(e) δ_E の受け渡しの契約 (上書きでは通さない): " + str(de.get("why"))[:1200])
    settled = bool(s == "REACHED" and bind["pass"] and g_mass["pass"] and g_tgt["pass"])
    src = {"kind": "stop_rule_v2_stop_snapshot", "n_snapshots": 1, "res": (f"res_{k}.h5" if k is not None else None),
           "steps": ([k] if k is not None else None), "stop_step": k, "is_stop_output": k is not None,
           "res_sha256": (st_full or {}).get("stop_sha256"), "dir": f"{SR.V2_DIR}/delta_E_{k}",
           "delta_r_next_sha256": de.get("delta_r_next_sha256"),
           "provider": "stop_rule.watch_main_stage_v2 (extract_delta_e、止めた step の全場から 1 回; 関門は抽出し直さない)"}
    mark = {"settled": settled, "unsettled_overrides": unsettled, "delta_e_source": src}
    # (g)(h): 停止規則の抽出の表を写して判定 ((a)〜(f) を満たしたときだけ; 上書き不可)
    ext = {"band_select": HANDOFF_BAND_SELECT, "settings": handoff_extract_settings(), "dir": HANDOFF_EXTRACT_DIR,
           "source": src["dir"], "made": None, "summary": None, "delta_r_next_sha256": de.get("delta_r_next_sha256")}
    if blocked or not extract:
        why = "前の条件で止めた" if blocked else "extract=False"
        checks["coverage"] = {"pass": False, "evaluated": False, "why": why}
        checks["next_wall"] = {"pass": False, "evaluated": False, "why": why}
        if not blocked:
            blocked.append("抽出の表を判定していない (extract=False) — 被覆・held・次の壁を判定していない")
    else:
        od = run_dir / HANDOFF_EXTRACT_DIR
        if od.exists():
            shutil.rmtree(od)
        od.mkdir(parents=True)
        srcd = Path(de["dir"])
        for f in _V2_TABLE_FILES + _V2_SENS_FILES:
            if (srcd / f).is_file():
                shutil.copy2(srcd / f, od / f)
        dn = od / "delta_r_next.csv"
        summ = dict(de["summary"])
        summ.update(mark)
        (od / "delta_r_extract_summary.json").write_text(json.dumps(summ, indent=1, default=str))
        ext["summary"] = summ
        if _sha256(dn) != de.get("delta_r_next_sha256"):
            blocked.append(f"写した {dn} の sha256 が停止規則の抽出と違う (写しの失敗)")
            checks["coverage"] = {"pass": False, "evaluated": False, "why": "写しの失敗"}
            checks["next_wall"] = {"pass": False, "evaluated": False, "why": "写しの失敗"}
        else:
            x_e = ((info.get("physical_wall") or {}).get("domain") or [None, None])[1]
            if x_e is None:
                x_e = float(np.max(np.loadtxt(dn, delimiter=",", skiprows=1)[:, 0]))
            checks["coverage"] = coverage_check(summ, dn, float(x_e))
            checks["next_wall"] = next_wall_check(problem, dn)
            for name, label in (("coverage", "抽出の全域の被覆・使用域の held と欠測"), ("next_wall", "次の壁の幾何の健全性")):
                if not checks[name]["pass"]:
                    blocked.append(f"{label}: {checks[name].get('why')}")
        if not blocked:
            made = []
            for f in _V2_TABLE_FILES:
                if (od / f).exists():
                    shutil.copy2(od / f, run_dir / f)
                    made.append(f)
            ext["made"] = made
    rec = {"version": 4, "method": "v2",
           "method_note": ("新方式 (settle_criteria.yaml のある run): 停止規則 v2 の REACHED・停止の step の再評価 (SETTLE_VERDICT.json) の全体・"
                           "本段の起動との照合・停止規則の δ_E の抽出 (抽出し直さない)"),
           "plan": HANDOFF_V2_PLAN, "run": str(run_dir), "time": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
           "proceed": not blocked, "blocked_by": blocked, "warnings": warnings,
           "settled": mark["settled"], "unsettled_overrides": mark["unsettled_overrides"], "delta_e_source": mark["delta_e_source"],
           "stop_rule": st, "convergence": cv, "checks": checks,
           "override": {"allow_not_reached": why_nr is not None, "allow_unconverged": why_uc is not None, "allow_unsettled": why_us is not None,
                        "allow_not_reached_reason": why_nr, "allow_unconverged_reason": why_uc, "allow_unsettled_reason": why_us,
                        "applied": applied,
                        "scope": "上書きは (a) 停止の状態 NOT_USED・CENSORED (allow_not_reached)・(b) 残差の閾値の超過 (allow_unconverged)・"
                                 "(c) 流量の収支と目的の量の閾値の超過 (allow_unsettled) だけ。照合 (状態・再評価・STOPPED_AT・停止の step・"
                                 "本段の起動・基準・見張りの記録)・判定不能・非有限・FAILED・y₁⁺・δ_E の契約・参照・ω・被覆・held・次の壁は"
                                 "上書きでは通さない"},
           "extract": ext,
           "note": "proceed は §4.9 の 1 の条件の全部を満たしたか (壁の更新は次の pass へ自動でつながない)。判定は「設定された量・区間について "
                   "SETTLED」で、旧の低下の桁数の判定との一致は求めない。帯の位置の感度は抽出の方法への感度で、時間の整定や誤差の上限の証明では"
                   "ない (合否に使わない)。settled: false は落ち着きの条件を上書きで通した「落ち着いていない場からの設計」"}
    if write:
        (run_dir / HANDOFF_CONVERGENCE).write_text(_v2_convergence_text(k, v, {"residual": g_res, **({"residual_watch": g_watch} if g_watch else {}),
                                                                                "mass": g_mass, "targets": g_tgt},
                                                                        bind["criteria_sha256"]))
        (run_dir / HANDOFF_WALL_RESOLUTION).write_text(wr["output"])
        _write_json_atomic(run_dir / HANDOFF_GATE, rec)
    return rec


def _write_json_atomic(path: Path, rec: dict) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=str))
    os.replace(tmp, path)


def _print_gate(label: str, gate: dict) -> None:
    ov = gate["override"]
    used = "; ".join("{} ({})".format(k, ov.get(k + "_reason")) for k in ov.get("applied") or [])
    print(f"handoff gate ({label}, 方式 {gate.get('method', 'v1')}): {'進む' if gate['proceed'] else '止める'} — 停止の状態 "
          f"{gate['stop_rule']['status']}、収束の判定 {gate['convergence']['kind']}"
          + (f"、理由: {'; '.join(gate['blocked_by'])}" if gate["blocked_by"] else "")
          + (f"、上書き {used}" if used else "")
          + ("" if gate.get("settled", True) else "、settled: false (落ち着いていない場からの設計)"), flush=True)


def _check_overrides(allow_not_reached, allow_unconverged, allow_unsettled) -> None:
    """上書きの理由 (`override_reason`) を NS の前に検査する (理由の無い上書きで本段を回してから止めない)。"""
    override_reason("allow_not_reached", allow_not_reached)
    override_reason("allow_unconverged", allow_unconverged)
    override_reason("allow_unsettled", allow_unsettled)


def _after_ns(run_dir, euler_ref, omega, m: dict, label: str, allow_not_reached: str | None, allow_unconverged: str | None,
              problem=None, allow_unsettled: str | None = None) -> dict:
    """NS の後段 (`run_pass`・`run_pass0_integral`・`run_from_table` 共通): 受け渡しの関門 (`handoff_gate`; 条件を満たせば関門の中で抽出し、
    表を run の直下に写す) → 通れば fixed_point.json (抽出の要約 = 落ち着きの印 settled・unsettled_overrides・delta_e_source を含む) →
    設計の標準出力。関門で止めたら fixed_point.json・delta_r_next.csv を run の直下に作らない (標準出力は今までどおり作る)。m (metrics) に
    関門の記録 `handoff_gate` を足して返す。新方式 (settle_criteria.yaml のある run、plan tooling-settle-convergence-and-monitor §4.9 の 1) では
    関門は抽出し直さず、停止規則 v2 が止めた step で抽出した表 (`stop_rule_v2/delta_E_<step>/`) を照合して同じ名前で run の直下に写すので、
    ここ (fixed_point.json = その要約 + 落ち着きの印) は旧方式と同じ流れ。"""
    run_dir = Path(run_dir)
    gate = handoff_gate(run_dir, problem=problem, euler_ref=euler_ref, omega=omega, allow_not_reached=allow_not_reached,
                        allow_unconverged=allow_unconverged, allow_unsettled=allow_unsettled)
    print((run_dir / HANDOFF_CONVERGENCE).read_text(), flush=True)      # 収束の判定の出力 (今までは check_convergence がそのまま出した)
    _print_gate("new", gate)
    if gate["proceed"]:
        fp = gate["extract"]["summary"]
        (run_dir / "fixed_point.json").write_text(json.dumps(fp, indent=1, default=str))
        print(f"{label}:", json.dumps({k: v for k, v in fp.items() if k != "massflow"}, default=str), flush=True)
        print("massflow(new):", json.dumps(fp["massflow"], default=str), flush=True)
        gate["extract"]["made"] = ["fixed_point.json"] + list(gate["extract"]["made"] or [])
        _write_json_atomic(run_dir / HANDOFF_GATE, gate)
    else:
        print("handoff gate: 抽出の表 (fixed_point.json・delta_r_next.csv) を run の直下に作らない — 判定は " + str(run_dir / HANDOFF_GATE),
              flush=True)
    m["handoff_gate"] = gate
    _design_report(run_dir, euler_ref)
    return m


def run_pass(problem, euler_ref, prev_run, run_dir, omega: float = 0.5, ic_from=None,
             prepare_only: bool = False, nsteps=None, cfl_main=None, implicit_relax=None,
             stages: str = "full", ramp=None, ramp_steps: int = 1000, stop_rule: bool | None = None,
             allow_not_reached: str | None = None, allow_unconverged: str | None = None, allow_unsettled: str | None = None) -> dict:
    """前の pass の run (prev_run) から抽出して次の壁を作り、NS を準備・実行する 1 pass (自動で次の pass へはつながない)。
    prev_run の受け渡しの関門 (`handoff_gate`、HANDOFF_GATE.json は prev_run に書かない; 条件を満たせば関門の中で固定の条件で抽出し、表を
    prev_run の直下に写す) が不合格なら準備の前に `HandoffBlocked` (新しい run dir を作らない)。関門の結果は新しい run の prepare_info.json
    の `handoff_from_prev` に残す (前の run の関門が落ち着きの条件を上書きで通したら、その `settled: false`・`unsettled_overrides` も
    ここに入る = 落ち着いていない場から作った壁の run; §5.1 #36h)。NS の後は `_after_ns`。stop_rule: prepare_ns へ (None = 問題 YAML の
    evaluate.stop_rule、キー無しは使う; False = 使わない)。allow_not_reached・allow_unconverged・allow_unsettled は上書きの理由の文字列 (`override_reason`;
    前の run と新しい run の関門の両方に効く)。"""
    from ..evaluate.runner_axismach import prepare_ns, run_staged_ns, collect
    _check_overrides(allow_not_reached, allow_unconverged, allow_unsettled)
    prev_run = Path(prev_run); run_dir = Path(run_dir)
    gate_prev = handoff_gate(prev_run, problem=problem, euler_ref=euler_ref, omega=omega, allow_not_reached=allow_not_reached,
                             allow_unconverged=allow_unconverged, allow_unsettled=allow_unsettled, write=False)
    _print_gate("prev", gate_prev)
    if not gate_prev["proceed"]:
        raise HandoffBlocked(f"前の pass の run {prev_run} が受け渡しの関門で不合格 — 抽出・次の壁の作成へ進まない: "
                             + "; ".join(gate_prev["blocked_by"]))
    summ = gate_prev["extract"]["summary"]
    print("extract(prev):", json.dumps({k: v for k, v in summ.items() if k != "massflow"}, default=str), flush=True)
    print("massflow(prev):", json.dumps(summ["massflow"], default=str), flush=True)
    info = prepare_ns(problem, run_dir, nsteps=nsteps, ic_from=ic_from or prev_run,
                      delta_r_csv=prev_run / "delta_r_next.csv", offset="radial",
                      euler_ref=euler_ref, omega=omega, prev_run=prev_run,
                      cfl_main=cfl_main, implicit_relax=implicit_relax, stop_rule=stop_rule)
    info["stages"] = {"stages": stages, "ramp": (list(ramp) if ramp else None), "ramp_steps": ramp_steps}
    info["handoff_from_prev"] = gate_prev
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    print(json.dumps({k: info[k] for k in ("throat_physical", "dstar_source", "mesh", "nStepOuter", "cfl_main")}, indent=1), flush=True)
    print(Path(run_dir, "MESH_QUALITY.txt").read_text().splitlines()[-1], flush=True)
    if prepare_only:
        return info
    import time as _t
    t0 = _t.time()
    rc = run_staged_ns(run_dir, stages=stages, ramp=ramp, ramp_steps=ramp_steps)
    print(f"forge rc {rc}  (NS wall time {_t.time() - t0:.0f} s, stages={stages}; 停止の状態 {stop_state(run_dir)['status']})", flush=True)
    m = collect(problem, run_dir)
    m["ns_wall_time_s"] = _t.time() - t0
    (run_dir / "metrics.json").write_text(json.dumps(m, indent=1))
    print(json.dumps(m, indent=1), flush=True)
    # 受け渡しの関門 → 固定点差: 新 run から抽出し、入力 (= 新 run の壁 − Euler 壁) と比較
    return _after_ns(run_dir, euler_ref, omega, m, "fixed_point(new)", allow_not_reached, allow_unconverged, problem=problem,
                     allow_unsettled=allow_unsettled)


def resolve_integral_initializer(problem) -> dict:
    """`--init-integral` の積分法初期化: problem YAML の `deltastar_initializer` があればそれ、無ければ断熱 contur (cf_scale 1)。
    2026-10-07 までは YAML を読まずに常に {"model": "contur"} を渡していた (help・docstring と食い違い)。
    C2 で較正した k_f (cf_scale) を YAML に書いた生産問題で、prep_c2pin.py (k_f を明示) と同じ物理壁にするため
    (plan tooling-nozzle-throat-monotone-r2 §9 2026-10-07: 読まないと出口半径が 2.4 mm ずれた)。"""
    from ..evaluate.runner_axismach import load_problem
    cfg = load_problem(problem).raw.get("deltastar_initializer")
    if cfg is None:
        return {"model": "contur"}
    if not isinstance(cfg, dict) or "model" not in cfg:
        raise ValueError(f"deltastar_initializer は model を持つ mapping であること ({cfg!r})")
    return dict(cfg)


def run_pass0_integral(problem, euler_ref, run_dir, ic_from, initializer=None, prepare_only: bool = False,
                       nsteps=None, omega: float = 0.5, cfl_main=None, implicit_relax=None,
                       stages: str = "full", ramp=None, ramp_steps: int = 1000, stop_rule: bool | None = None,
                       allow_not_reached: str | None = None, allow_unconverged: str | None = None,
                       allow_unsettled: str | None = None) -> dict:
    """pass 0: 積分法 (CONTUR) 初期壁で NS を立て、終了後に受け渡しの関門 (`_after_ns`) を通れば抽出して次 pass 用 delta_r_next.csv まで作る。
    initializer=None なら problem YAML の `deltastar_initializer` (無ければ断熱 contur)。上書きの理由 (`override_reason`) は NS の前に検査する。"""
    from ..evaluate.runner_axismach import prepare_ns, run_staged_ns, collect
    _check_overrides(allow_not_reached, allow_unconverged, allow_unsettled)
    run_dir = Path(run_dir)
    # thermal_bc は prepare_ns 側で spec.wall_thermal から決める (ここでは model のみ)
    init = initializer if initializer is not None else resolve_integral_initializer(problem)
    info = prepare_ns(problem, run_dir, nsteps=nsteps, ic_from=ic_from, initializer=init,
                      euler_ref=euler_ref, omega=omega, cfl_main=cfl_main, implicit_relax=implicit_relax, stop_rule=stop_rule)
    info["stages"] = {"stages": stages, "ramp": (list(ramp) if ramp else None), "ramp_steps": ramp_steps}
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    print(json.dumps({k: info[k] for k in ("throat_physical", "dstar_source", "initializer", "mesh", "nStepOuter", "cfl_main")},
                     indent=1, default=str), flush=True)
    print(Path(run_dir, "MESH_QUALITY.txt").read_text().splitlines()[-1], flush=True)
    if prepare_only:
        return info
    import time as _t
    t0 = _t.time()
    rc = run_staged_ns(run_dir, stages=stages, ramp=ramp, ramp_steps=ramp_steps)
    print(f"forge rc {rc}  (NS wall time {_t.time() - t0:.0f} s, stages={stages}; 停止の状態 {stop_state(run_dir)['status']})", flush=True)
    m = collect(problem, run_dir)
    m["ns_wall_time_s"] = _t.time() - t0
    (run_dir / "metrics.json").write_text(json.dumps(m, indent=1))
    print(json.dumps(m, indent=1), flush=True)
    return _after_ns(run_dir, euler_ref, omega, m, "extract(new)", allow_not_reached, allow_unconverged, problem=problem,
                     allow_unsettled=allow_unsettled)


def check_table_against_card(problem, table) -> dict:
    """`--delta-r-csv TABLE` の照合 (plan tooling-nozzle-deltastar-band-root-zone §4.3 3): 問題 YAML の `_card.ns_correction`
    (`wt_card --ns-correction` が書く NS 後の補正の記録) があれば、TABLE の sha256 がその `table_sha256` と一致すること (違えば ValueError —
    寸法を解いた表と別の表で壁を作らない)。記録が無ければ照合しない。戻り: {table, sha256, card_table, card_table_sha256, checked}。"""
    from ..probdef import load_problem
    table = Path(table)
    if not table.exists():
        raise ValueError(f"--delta-r-csv: {table} が無い")
    sha = _sha256(table)
    rec = (load_problem(problem).raw.get("_card") or {}).get("ns_correction")
    out = {"table": str(table), "sha256": sha, "card_table": None, "card_table_sha256": None, "checked": False}
    if rec is None:
        return out
    want = rec.get("table_sha256") if isinstance(rec, dict) else None
    if not want:
        raise ValueError(f"{problem}: _card.ns_correction に table_sha256 が無い (壊れた記録; wt_card --ns-correction で作り直す)")
    out.update(card_table=rec.get("table"), card_table_sha256=want, checked=True)
    if sha != want:
        raise ValueError(f"--delta-r-csv {table} の sha256 {sha[:16]}… が問題 YAML の記録 (_card.ns_correction.table {rec.get('table')!r}、"
                         f"sha256 {want[:16]}…) と違う — 寸法を解いた補正の表と別の表で壁を作らない (wt_card --ns-correction で作り直す)")
    return out


def run_from_table(problem, euler_ref, table, run_dir, ic_from=None, prev_run=None, prepare_only: bool = False, nsteps=None,
                   omega: float = 0.5, cfl_main=None, implicit_relax=None, stages: str = "full", ramp=None, ramp_steps: int = 1000,
                   stop_rule: bool | None = None, allow_not_reached: str | None = None, allow_unconverged: str | None = None,
                   allow_unsettled: str | None = None) -> dict:
    """補正の表 TABLE (`wt_card --ns-correction` が問題 YAML の横に書く `<stem>_delta_r.csv`、または `--solve-rt-throat --prev` の
    `--delta-r-out`) で物理壁を作り、NS を準備・実行する (plan tooling-nozzle-deltastar-band-root-zone §4.3 3):
    `prepare_ns(delta_r_csv=TABLE, offset="radial")` → `run_staged_ns` → collect → 受け渡しの関門 (`_after_ns`、plan
    tooling-nozzle-isothermal-wall-chain §5.1 #36a) → 通れば抽出 (`fixed_point.json`) → 設計の標準出力 (`run_pass` と同じ後段)。
    表は緩和しない (omega は後段の抽出の delta_r_next.csv にだけ使う)。問題 YAML の `_card.ns_correction` が
    あれば TABLE の sha256 を照合する (`check_table_against_card`、違えば run dir を作る前に例外)。照合の結果は prepare_info.json の
    `delta_r_table` に書く。表は呼び出し側が渡すので、前の run の関門はここでは見ない (prev_run は記録と初期場だけ)。
    上書きの理由 (`override_reason`) は NS の前に検査する。"""
    from ..evaluate.runner_axismach import prepare_ns, run_staged_ns, collect
    _check_overrides(allow_not_reached, allow_unconverged, allow_unsettled)
    run_dir = Path(run_dir)
    chk = check_table_against_card(problem, table)
    info = prepare_ns(problem, run_dir, nsteps=nsteps, ic_from=ic_from, delta_r_csv=table, offset="radial",
                      euler_ref=euler_ref, omega=None, prev_run=prev_run, cfl_main=cfl_main, implicit_relax=implicit_relax,
                      stop_rule=stop_rule)
    info["stages"] = {"stages": stages, "ramp": (list(ramp) if ramp else None), "ramp_steps": ramp_steps}
    info["delta_r_table"] = chk
    (run_dir / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    print(json.dumps({k: info[k] for k in ("throat_physical", "dstar_source", "delta_r_table", "mesh", "nStepOuter", "cfl_main")},
                     indent=1, default=str), flush=True)
    print(Path(run_dir, "MESH_QUALITY.txt").read_text().splitlines()[-1], flush=True)
    if prepare_only:
        return info
    import time as _t
    t0 = _t.time()
    rc = run_staged_ns(run_dir, stages=stages, ramp=ramp, ramp_steps=ramp_steps)
    print(f"forge rc {rc}  (NS wall time {_t.time() - t0:.0f} s, stages={stages}; 停止の状態 {stop_state(run_dir)['status']})", flush=True)
    m = collect(problem, run_dir)
    m["ns_wall_time_s"] = _t.time() - t0
    (run_dir / "metrics.json").write_text(json.dumps(m, indent=1))
    print(json.dumps(m, indent=1), flush=True)
    return _after_ns(run_dir, euler_ref, omega, m, "extract(new)", allow_not_reached, allow_unconverged, problem=problem,
                     allow_unsettled=allow_unsettled)


def _gate_rc(res) -> int:
    """CLI の終了コード: 関門で止めた (run の後段の handoff_gate が proceed でない) なら 3、それ以外 0 (prepare-only は 0)。"""
    gate = res.get("handoff_gate") if isinstance(res, dict) else None
    return 3 if (gate is not None and not gate["proceed"]) else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="δ_r 固定点反復 (固定 Euler 基準)")
    ap.add_argument("--problem", required=True)
    ap.add_argument("--euler-ref", required=True, help="固定の Euler の参照の run dir")
    ap.add_argument("--euler-ref-res", default=None, metavar="RES",
                    help="Euler の参照のスナップショットを固定する: 読む出力のファイル名 (res_<N>.h5)。--euler-ref-sha256 と一緒に書く "
                         "(plan tooling-nozzle-isothermal-wall-chain §5.1 #36e (6)。無ければ今までどおり最新の出力)")
    ap.add_argument("--euler-ref-sha256", default=None, metavar="HEX", help="--euler-ref-res の sha256 (合わなければ例外)")
    ap.add_argument("--prev", default=None, help="前 pass の NS run (抽出元・IC 既定)。--init-integral では不要")
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--omega", type=float, default=0.5)
    ap.add_argument("--ic-from", default=None)
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--prepare-only", action="store_true")
    ap.add_argument("--extract-only", action="store_true",
                    help="抽出と delta_r_next.csv だけ作る (受け渡しの関門を通さない手作業の抽出; 帯は extract_and_merge の既定 adaptive)")
    ap.add_argument("--init-integral", action="store_true",
                    help="pass 0: 積分法 (CONTUR) 初期壁で NS を立てる (YAML の deltastar_initializer を使う)")
    ap.add_argument("--init-thermal", default=None, help="廃止 (2026-09-12): 壁温は problem YAML の spec.wall_thermal が単一ソース")
    ap.add_argument("--cfl", type=float, default=None, help="本段 cfl (YAML evaluate.cfl_main を上書き)")
    ap.add_argument("--implicit-relax", type=float, default=None, help="implicitRelax (YAML evaluate.implicit_relax を上書き)")
    ap.add_argument("--stages", default="full", choices=("full", "none", "ramp"), help="起動: full=soft/mid/本段, none=本段のみ, ramp=cfl を段階的に上げて本段")
    ap.add_argument("--ramp", default=None, help="--stages ramp の cfl 列 (例 '1,2,3.5')")
    ap.add_argument("--ramp-steps", type=int, default=1000)
    ap.add_argument("--solve-rt", type=float, default=None, metavar="R_EXIT_M",
                    help="出口物理半径 [m] を与えて r_t を解くだけ (--prev があれば実測 δ_r で、無ければ prepare_ns と同じ積分法の経路で)")
    ap.add_argument("--solve-rt-throat", type=float, default=None, metavar="R_THROAT_M",
                    help="物理スロート半径 [m] を与えて r_t を解くだけ (--prev があれば δ_E の全分布から作る補正関数で。"
                         "そのとき --delta-r-out に次の壁に渡す表を書く)")
    ap.add_argument("--delta-r-out", default=None, help="--solve-rt-throat --prev のとき、解いた r_t での補正関数の表の出力先")
    ap.add_argument("--delta-next", default=None, help="--solve-rt-throat --prev のとき、δ_E を読む delta_r_next.csv (既定 PREV/delta_r_next.csv)")
    ap.add_argument("--delta-r-csv", default=None, metavar="TABLE",
                    help="補正の表 (wt_card --ns-correction が YAML の横に書く <stem>_delta_r.csv など) で物理壁を作って NS を準備・実行する "
                         "(prepare_ns(delta_r_csv=TABLE, offset=radial) → run_staged_ns)。YAML の _card.ns_correction があれば TABLE の "
                         "sha256 を照合する。IC は --ic-from (無ければ --prev)。--init-integral とは排他")
    ap.add_argument("--no-stop-rule", action="store_true",
                    help="本段の停止規則を使わない (問題 YAML の evaluate.stop_rule を打ち消す; plan tooling-nozzle-cooled-wall-ns-settling §4.6・"
                         "plan tooling-nozzle-isothermal-wall-chain §5.1 #36a)。既定は使う = 上限 --steps か evaluate.stop_rule_cap か 200,000 "
                         "(到達 + 保持を含む) まで回し、r と θ_r の条件を満たして保持 (evaluate.stop_rule_hold_steps、既定 40,000) を終えた出力で止める。"
                         "使わないときの本段は --steps > evaluate.ns_steps > 60,000 step")
    ap.add_argument("--allow-not-reached", default=None, metavar="理由",
                    help="受け渡しの関門の明示の上書き (理由の文字列が必須、HANDOFF_GATE.json に残す): 本段の停止の状態が REACHED でない run "
                         "(停止規則を使わない NOT_USED・CENSORED・UNAVAILABLE) からも抽出・次の壁の入力を作る。DIVERGED・EXEC_ERROR・"
                         "DATA_ERROR・記録なしは通さない。使ったら settled: false (落ち着いていない場からの設計) を記録する")
    ap.add_argument("--allow-unconverged", default=None, metavar="理由",
                    help="受け渡しの関門の明示の上書き (理由の文字列が必須、HANDOFF_GATE.json に残す): check_convergence.py --segment が "
                         "NOT CONVERGED でも抽出・次の壁の入力を作る (DIVERGED・判定不能は通さない)")
    ap.add_argument("--allow-unsettled", default=None, metavar="理由",
                    help="受け渡しの関門の明示の上書き (理由の文字列が必須、HANDOFF_GATE.json に残す; plan tooling-nozzle-isothermal-wall-chain "
                         "§5.1 #36h): 保存則 (停止規則の窓の質量の収支) と δ_E・δ_target の直近 2 窓の変動が閾値を超えても抽出・次の壁の入力を作る。"
                         "実測の値を記録し settled: false (落ち着いていない場からの設計) にする。判定不能 (記録の欠け・窓を覆わない・点数・"
                         "非有限・held・欠測)・y₁⁺・参照・ω・被覆・held・次の壁は通さない")
    a = ap.parse_args(argv)
    for flag, v in (("--allow-not-reached", a.allow_not_reached), ("--allow-unconverged", a.allow_unconverged),
                    ("--allow-unsettled", a.allow_unsettled)):
        if v is not None and not v.strip():
            ap.error(f"{flag} には上書きの理由の文字列を書く (空は不可)")
    if (a.euler_ref_res is None) != (a.euler_ref_sha256 is None):
        ap.error("--euler-ref-res と --euler-ref-sha256 は一緒に書く (Euler の参照のスナップショットの固定)")
    if a.euler_ref_res is not None:
        # 固定のスナップショット {run, res, sha256} (§5.1 #36e (6)): ここで形と sha を確かめる (違えば例外)
        from ..metrics.deltastar import resolve_euler_ref
        a.euler_ref = {"run": a.euler_ref, "res": a.euler_ref_res, "sha256": a.euler_ref_sha256}
        resolve_euler_ref(a.euler_ref, verify=True)
    sr_arg = False if a.no_stop_rule else None      # None = 問題 YAML の evaluate.stop_rule (無ければ使う)
    gate_kw = dict(allow_not_reached=a.allow_not_reached, allow_unconverged=a.allow_unconverged, allow_unsettled=a.allow_unsettled)
    if a.delta_r_csv is not None:
        clash = [n for n, v in (("--init-integral", a.init_integral), ("--extract-only", a.extract_only),
                                ("--solve-rt", a.solve_rt is not None), ("--solve-rt-throat", a.solve_rt_throat is not None)) if v]
        if clash:
            ap.error(f"--delta-r-csv と {', '.join(clash)} は排他")
        ramp = tuple(float(v) for v in a.ramp.split(",")) if a.ramp else None
        res = run_from_table(a.problem, a.euler_ref, a.delta_r_csv, a.run_dir, ic_from=(a.ic_from or a.prev), prev_run=a.prev,
                             prepare_only=a.prepare_only, nsteps=a.steps, omega=a.omega, cfl_main=a.cfl, implicit_relax=a.implicit_relax,
                             stages=a.stages, ramp=ramp, ramp_steps=a.ramp_steps, stop_rule=sr_arg, **gate_kw)
        return _gate_rc(res)
    if a.solve_rt is not None:
        r = solve_rt(a.problem, a.solve_rt, prev_run=a.prev, euler_run=a.euler_ref)
        print(json.dumps(r, indent=1, default=str))
        return 0
    if a.solve_rt_throat is not None:
        r = solve_rt_throat(a.problem, a.solve_rt_throat, prev_run=a.prev, delta_r_out=a.delta_r_out, delta_next=a.delta_next)
        print(json.dumps(r, indent=1, default=str))
        return 0
    if a.extract_only:
        s = extract_and_merge(a.prev, a.euler_ref, omega=a.omega)
        print(json.dumps(s, indent=1))
        return 0
    if a.init_integral:
        init = None
        if a.init_thermal:
            ap.error("--init-thermal は廃止。spec.wall_thermal (problem YAML) で指定する")
        ramp = tuple(float(v) for v in a.ramp.split(",")) if a.ramp else None
        res = run_pass0_integral(a.problem, a.euler_ref, a.run_dir, a.ic_from, initializer=init,
                                 prepare_only=a.prepare_only, nsteps=a.steps, omega=a.omega, cfl_main=a.cfl,
                                 implicit_relax=a.implicit_relax, stages=a.stages, ramp=ramp, ramp_steps=a.ramp_steps,
                                 stop_rule=sr_arg, **gate_kw)
        return _gate_rc(res)
    if not a.prev:
        ap.error("--prev が必要 (--init-integral でなければ)")
    ramp = tuple(float(v) for v in a.ramp.split(",")) if a.ramp else None
    try:
        res = run_pass(a.problem, a.euler_ref, a.prev, a.run_dir, omega=a.omega, ic_from=a.ic_from,
                       prepare_only=a.prepare_only, nsteps=a.steps, cfl_main=a.cfl, implicit_relax=a.implicit_relax,
                       stages=a.stages, ramp=ramp, ramp_steps=a.ramp_steps, stop_rule=sr_arg, **gate_kw)
    except HandoffBlocked as e:
        print(f"handoff gate: {e}", flush=True)
        return 3
    return _gate_rc(res)



def _design_report(run_dir, euler_ref):
    """ノズル設計の標準出力 (procedures/nozzle-design-outputs.md) を NS pass の後に自動で作る。失敗しても chain は止めない。"""
    if os.environ.get("FORGE_NO_DESIGN_REPORT"):
        return
    try:
        from ..report.nozzle_report import make_report
        # 固定のスナップショットの指定 {run, res, sha256} は run dir を渡す (標準出力の図・評価量は診断で、Euler は最新の出力を読む)
        print("design report:", make_report(run_dir, euler_ref["run"] if isinstance(euler_ref, dict) else euler_ref), flush=True)
    except BaseException as e:  # noqa: BLE001
        print(f"design report: 失敗 ({type(e).__name__}: {e}) — chain は続行", flush=True)


if __name__ == "__main__":
    sys.exit(main())
