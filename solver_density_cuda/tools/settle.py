#!/usr/bin/env python3
r"""収束の判定の共通の部品 — 「残差が準定常・流量の収支が取れる・目的の量が準定常」の AND
(plan `plans/accepted/tooling-settle-convergence-and-monitor.md` §4.2・§4.3・§4.8 の 1 と 4、§5.1 #2〜#4)。

道具 `check_settled.py` (CLI) と停止規則 v2 (§5.1 #5) が同じこの部品を使う。設計チェーンには依存しない
(プリセットの解決は `design/forge_design/evaluate/settle_presets.py`)。

## 入力は run の `settle_criteria.yaml` だけ (§4.2)

- 無い・目的の量 (`exit_core_mach`・`line_mach`・`exit_flow_angle`・`wall_heat_nodes`) が 1 つも無い・残差 (`residual`) も
  残差の非有限の見張り (`residual_finite`) も無い・流量の収支 (`mass_balance`) が無い → `INDETERMINATE (評価対象が未設定)`
  (例外にしない; 置き場所と量の種類の一覧を文言に付ける)。
- 知らない種類・キー、閾値の欠け、YAML の重複キー → `CriteriaError` (黙って既定値で埋めない)。CLI は終了コード 3。
- 書式 (ひな形は `settle_criteria.template.yaml`):

      schema: forge.settle_criteria/v1
      source: {...}                       # 出典 (プリセット・版・sha256)。記録だけ
      structured: {ni, nj, scale_m, mesh_file}   # 節点の量に要る (構造格子の節点番号 = i·nj + j、j = 0 が軸、j = nj−1 が壁;
                                          # mesh_file は任意 = 座標を読む格子の h5)
      grid: {out_step_interval, start_step, segment}
      window: {W, min_points}
      hold: 40000                         # 停止規則用 (記録だけ)
      quantities: [ {kind: ..., thresholds: {..., provisional: true}}, ... ]
      excluded: [ {kind, reason} ]        # 外した量と理由 (記録だけ)
      preconditions: [ {quantity, requires, satisfied} ]   # 記録だけ
      notes: [...]

## 判定の格子 (§4.8 の 1)

- 観測は `grid.out_step_interval` (= forge の `time.outStepInterval`) の**倍数の step だけ**。終端 T = `end_step` 以下の最後の格子点。
  `end_step` は CLI の `--end-step` > run の `STOPPED_AT` (協調停止) の step > 入力のデータの最後の step。
- 窓 W は端点込みで点を数える (W 20,000・間隔 2,000 → 11 点): C = {T−W, …, T}。直前の窓 P = {T−2W, …, T−W} (同じく 11 点、T−W を共有)。
  P は残差の上昇の検査だけが使う。`W` は間隔の倍数、`W/間隔 + 1 ≥ min_points` でなければ `CriteriaError`。
- 格子点 s の観測に要るデータ: 残差は行 s−1 (forge は場を `iStep+1` で書き、残差を `iStep` で記録する)、流量は行 s、節点は
  `res_nodes_s.h5`、壁は `res_<name>_<physID>_s.h5`。窓の格子点のどれかが欠ける・`grid.start_step` より前 → 判定不能 (データ不足)。
  **本段の開始の場 (s = start_step) には残差の行 s−1 も節点・流量の出力も無い** (forge は step 0 に全場 res_0 だけを書く) ので、
  開始から 2W の終端 (例 40,000) では P が欠けて判定不能になり、最初に判定できる終端は start + 2W + 間隔 (例 42,000)。
- `end_step` が格子に乗らない (協調停止の 80,001 など) とき、その観測は窓に入れず、別の検査だけに使う: 非有限なし、直前の格子の観測
  (T) からの差が各量の許容の内、窓に加えても各量の許容の内。
- step の単調・重複を検査する (CSV の step の重複・逆行、h5 の属性 `step` とファイル名の不一致 → 判定不能)。

## 量の種類 (§4.2・§4.3)

各量の結果に verdict (SETTLED / NOT_SETTLED / DIVERGED / INDETERMINATE)・使った値・理由の分類 (`threshold` 閾値の超過 /
`data` データ不足 / `unset` 未設定 / `nonfinite` 非有限) を書く。全体は AND で、優先は DIVERGED > INDETERMINATE > NOT_SETTLED > SETTLED。
判定不能は合格にしない。

- `residual` (`file` 既定 residual_history.csv、`.gz` 可; `outer_end` の行): `required` の列ごとに y = log10(rms)。
  格子に間引いた系列 (行 s−1) の C で 傾き |s|·W ≤ `slope_dec` (最小二乗)・振幅 p95 − p5 ≤ `amp_dec`、mean(C) − mean(P) ≤ `rise_dec`。
  **生の系列** (行 T−W−1〜T−1 の全行) で max(y) − median(y) ≤ `excursion_dec` (末尾の 1 点の跳ね・格子の間の短いバーストを落とす)。
  区間の行に非有限 (どの数値の列でも) → DIVERGED。`required` の列が無い・恒等ゼロ・一部ゼロ (負も) → 判定不能。`inactive_ok` の列は
  恒等ゼロなら外し、そうでなければ required と同じに判定する。**どちらにも無い列が値を持つ → 判定不能 (未設定)** (黙って見ない扱いにしない)。
  `rms_dq_*` は非有限の検査だけ。
  **`series: raw`** (プリセット `wt_nozzle_v2`; 既定 `grid` は上の通りで今と同じ): 傾き・振幅・上昇を格子点の行でなく**毎 step の
  outer_end の全行**で測る。区間 (T−W, T] と直前の区間 (T−2W, T−W] の場の step を行 s−1 の対応で行 T−W〜T−1・T−2W〜T−W−1
  (各 W 行) にし、傾き |s|·W (全行の最小二乗)・振幅 p95 − p5・上昇 mean(区間) − mean(直前の区間) を同じ閾値で判定する。
  区間ごとに整数 step の行が**重複なく全部**あることを検査し、欠け・重複は判定不能 (データ不足)。直前の区間が本段の開始の
  直後から始まる終端 (T − 2W = 開始) も判定不能 (開始直後の行を直前の区間に入れない; 最初に判定できる終端は grid と同じ
  開始 + 2W + 間隔)。跳ね・停止の観測の検査は grid と同じ。
- `residual_finite` (`file` 既定 residual_history.csv; プリセット `wt_nozzle_v3`、plan §4.4「残差を止める判定から外す」(2026-10-16
  ユーザ決定)・§5.1 #8): **残差の非有限の見張りで、判定の条件ではない**。区間 (outer_end の行 start〜終端 − 1) のどの数値の列にも
  非有限が無ければ SETTLED、非有限 → DIVERGED。履歴が無い・区間の最後の行が無い → 判定不能。`residual` を外した基準はこれが要る。
  **見張りの被覆** (§5.1 #10、result 段レビュー 2 回目 M2; 閾値は戻さない): 区間の整数 step の outer_end の行が**重複なく全部**
  (forge は毎 step 書く; `residual` の `series: raw` と同じ定義)・step の逆行なし、**必要な列** (`required`; 書かなければ
  `RESIDUAL_FINITE_DEFAULT_REQUIRED` = forge が常に書く保存量の 4 列 `residual` の既定と同じ名前) が全部ある、を必須にする。
  欠け・重複・逆行・列の欠け → 判定不能 (データ不足)。以前は在る数値だけを見て最終行があれば合格にしたので、最終 1 行だけ・
  rms_ro の列だけ・残差の列の無い `step,inner,phase` の 1 行でも SETTLED になった。乱流 (rms_roK・rms_roOmega)・化学種などの列は
  基準の `required` に書いたときだけ要る (書くのは生成器 `settle_presets.py` の役目)。
- `mass_balance` (`file` 既定 boundary_massflux.csv; forge の `output.boundaryMassflux: 1`): 列 `mdot_<name>_<physID>` [kg/s] から
  ṁ_in = −Σ_inlet mdot、ṁ_out = Σ_outlet mdot、r = (ṁ_in − ṁ_out)/ṁ_in。C で |mean r| ≤ `mean_abs`、max r − min r ≤ `swing`、
  (max − min)/mean の ṁ_in・ṁ_out ≤ `mdot_change_rel`。入口・出口以外の境界 (壁・軸) の mdot は記録だけ。ṁ_in ≤ 0 は判定不能 (r の
  分母が無い)。ṁ_in > 0 で ṁ_out ≤ 0 (流出が無い・逆流・符号の反転) は理由つきの NOT_SETTLED (r ≥ 1 は測れている; ṁ_out の相対の
  変化は計算しない; §5.1 #7c、result 段レビュー M3)。停止の観測も同じ規則。
- 節点の量 (`nodes_file` 既定 `res_nodes_{step}.h5`; forge の `output.nodeProbeFile`): 構造格子の節点番号 id = i·nj + j から
  (i, j) = divmod(id, nj)。座標と M は `design/forge_design/report/nozzle_report.py` の `load_field` と**演算の順まで同じ**:
  X = COORD[:, 0] / scale_m・R = COORD[:, 1] / scale_m (COORD の型のまま割る)、M = np.hypot(Ux, Uy) / sonic (VALUE を float64 に
  してから)。COORD は `structured.mesh_file` (格子の h5、load_field が読む `nozzle.h5`) があればその /MESH/COORD を NODE_ID で引いた値
  (res_nodes の COORD がそれを res_nodes の型に丸めた値と一致することを確かめる)、無ければ res_nodes の COORD。変換器は COORD を
  double で書き、forge は節点の座標を geom_float で書くので、FP64 のビルドは両者が同じ、**float のビルドは mesh_file が無いと
  load_field と一致しない** (res_nodes の COORD は float32)。節点の値に非有限 → DIVERGED。
  - `exit_core_mach`: 最後の列 i = ni−1 の**全節点**が要る。η = R/R[j=nj−1] ∈ `eta_range` の節点の M の単純平均
    (`design/forge_design/evaluate/stop_rule.py` の `recorded_quantities` の `exit_core_M`・`nozzle_report.metrics` の `M_mean` と同じ)。
    C の max − min ≤ `tol` (絶対値)。
  - `line_mach`: `nozzle_report.eta_line` と同じ計算を、列ごとに (軸の節点 j=0 の x、η を挟む 2 節点、壁の節点 j=nj−1) の部分集合で行う。
    xq = linspace(x_range_m[0]/scale_m, x_range_m[1]/scale_m, n)、d = 100 (M/Md − 1) [%]。量: maxabs = max|d|、
    halfamp = 100 (max M − min M)/(2 Md) (`case/71.m55_combustion_wt/_runs/run_0007_ns/v5c_series.py` と同じ式)、位置ごとの d。
    判定: C の maxabs の幅 ≤ `tol_maxabs`、halfamp の幅 ≤ `tol_halfamp`、**位置ごとの d の幅の最大 ≤ `tol_pos`**
    (位置が動くが maxabs・halfamp が一定の線を落とす)。部分集合で全節点と同じ値になる条件 (列の対が全格子で隣り合う・η を挟む 2 節点が
    隣り合う・軸と壁の節点がある) を毎回確かめ、満たさなければ判定不能 (データ不足)。
  - `exit_flow_angle`: `exit_core_mach` と同じ節点の θ = degrees(atan2(Uy, Ux))、節点ごとの C の幅の最大 ≤ `tol` [deg]。
- `wall_heat_nodes`: 壁の境界の h5 `res_<name>_<physID>_<step>.h5` (forge `outputBconds_H5_XDMF`、bcond の `outputHDFflg: 1`) の
  `VALUE/iface_Qf_eff` (積分済みの節点の荷重、固体向き正; `output.interfaceDiag: 1`・定常・等温壁のときだけ値を持ち、それ以外は NaN;
  `solver_density_cuda/conjugateWall.cpp` の `fillInterfaceDiagnostics`)。node の壁の出力は値が**局所の節点の順** (`inodes_l2g`) に
  並び、全 step で同じ順なので、節点は配列の位置で対応させる (COORD が step の間で同一であることを確かめる; 値の数 = COORD の節点の数で
  なければ node の壁の出力でないとして判定不能)。許容 (§4.8 の 4): 窓の全格子点で |Q_i − Q_ref,i| ≤ q_tol·A_i + eps_rel·|Q_ref,i|、
  Q_ref = 終端 T の値。**A_i の取り方 (2026-10-15 に決めた)**: 同じ h5 の `MESH/COORD`・`MESH/CONNE` (XDMF Mixed、node の壁は
  Polyline [2, 2, a, b]) から、節点に接する辺の半分の面積の和。軸対称 (`axisymmetric: true`、r = y) は 1 rad あたり
  A_i = Σ_e (|e|/2)·(3 r_i + r_j)/4 (半分の辺の上の ∫ r dl を厳密に; `iface_Qf_eff` は 1 rad あたりの [W/rad] なので単位が揃う;
  `methods/design/overview.md` の A_node = 2π Σ r_f |S_f| の 1 rad 分)、平面は A_i = Σ_e |e|/2 [m²/m]。三角形・四角形 (3D) は面積/節点数の
  和 (軸対称では判定不能)。forge の `msh.planes[ip].surfArea` (平面の面積) とは別の量で、許容の尺度にだけ使う。
  `iface_Qf_eff` に非有限 → 判定不能 (適用外の NaN と発散を区別できないので; 流れの発散は残差・節点の量が DIVERGED にする)。
  **`measure: range`** (プリセット `wt_nozzle_v2`; 既定 `terminal_diff` は上の「終端の値からの差」で今と同じ): 節点ごとの窓の幅
  max_s Q_i − min_s Q_i ≤ q_tol·A_i + eps_rel·|Q_i(T)| (相対の参照は終端 T の値; 終端からの差では窓の中の往復を見落とす —
  codex diagnose 2026-10-16 settle-v2-span10k の例 [0.25, 1.75, 1, 1, 1, 1] は差の最大 0.75・幅 1.5)。停止の観測は終端からの差と
  窓に加えた幅を同じ許容で検査する。

## 系列の入力の代替 (§6 V3(a) の再生用)

`exit_core_mach` は `series_csv: {path, step_col, value_col}`、`line_mach` は `{path, step_col, maxabs_col, halfamp_col}`、
`mass_balance` は `{path, step_col, value_col (= r), [mdot_in_col, mdot_out_col]}` でも与えられる。そのとき位置ごと (line_mach の
`tol_pos`) と、ṁ の列が無い mass_balance の `mdot_change_rel` は判定不能 (データ不足)。path は run からの相対か絶対。

## 結果 `SETTLE_VERDICT.json` (§4.3)

基準の中身と sha256・起動の識別 (forge_launches.jsonl の最後の行・STOPPED_AT)・格子・判定区間・終端の step・窓の格子点・使ったファイルの
一覧 (大きさ・sha256)・各量の結果。閾値に `provisional: true` があれば VERDICT の行に「閾値は暫定 (plan §6 V3 で決める)」と出す。
`judge` は起動の識別を**記録するだけで照合しない**。本段の起動との照合は `select_judge` (下) と受け渡しの関門が行う。

## 道具の選択と本段の起動との照合 (`select_judge`、§4.9 の 1・2)

基準のファイルの無い run と、**前の段であることを記録と独立に確かめられた** run (本段の開始の記録 `MAIN_STAGE_STARTED` が無く、
`stage_manifest.json` に本段の段が無い) は `check_convergence.py`。本段の記録がある run は必ず本段として照合し (記録と今の
solverConfig.yaml・最後の起動・基準の sha256)、食い違い・起動の記録の欠けは `check_settled.py` の INDETERMINATE (旧の道具へ切り替えない;
§5.1 #7c、result 段レビュー M1)。記録が無く前の段とも確かめられない run も INDETERMINATE。手で組んだ run (基準が記録を指さず、記録も
無い) は `check_settled.py`。基準があって判定不能になったときに旧の道具へ落とさない。run_case.sh・`stop_rule.postprocess`・受け渡しの
関門が同じ規則を使う (`check_settled.py --select`; `check_settled.py` は既定で照合する)。
"""
from __future__ import annotations

import csv
import datetime as _dt
import gzip
import hashlib
import json
import math
import os
import re
from pathlib import Path

import numpy as np

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

SCHEMA = "forge.settle_criteria/v1"
VERDICT_SCHEMA = "forge.settle_verdict/v1"
CRITERIA_FILE = "settle_criteria.yaml"
VERDICT_FILE = "SETTLE_VERDICT.json"
TEMPLATE = "solver_density_cuda/tools/settle_criteria.template.yaml"
PROVISIONAL_TEXT = "閾値は暫定 (plan §6 V3 で決める)"

VERDICTS = ("SETTLED", "NOT_SETTLED", "INDETERMINATE", "DIVERGED")
_PRIORITY = {"SETTLED": 0, "NOT_SETTLED": 1, "INDETERMINATE": 2, "DIVERGED": 3}
EXIT_CODE = {"SETTLED": 0, "NOT_SETTLED": 1, "DIVERGED": 2, "INDETERMINATE": 3}

# 理由の分類
CAT_THRESHOLD = "threshold"
CAT_DATA = "data"
CAT_UNSET = "unset"
CAT_NONFINITE = "nonfinite"
CAT_GRID = "grid"
CAT_LABEL = {CAT_THRESHOLD: "閾値の超過", CAT_DATA: "データ不足", CAT_UNSET: "未設定", CAT_NONFINITE: "非有限",
             CAT_GRID: "格子と線の定義の不整合"}
# 線の補間で軸の節点 (j = 0) が持ってよい重みの上限 (codex diagnose 2026-10-16 settle-v1-and-line-height Major 3)。
# η が第 1 層の η₁ と丸めの差で一致する列では重みは ~1e−8、η₁ > η の格子では軸の節点の値が線に混ざる (軸の節点を避けて
# 近軸の線を取る意図が崩れる) ので判定不能にする。
AXIS_WEIGHT_MAX = 1.0e-6

# 残差の傾き・振幅・上昇を測る系列 (`residual.series`; plan §4.4 の「判定の区間を 1 万 step にする」、§5.1 #5d)。
# grid (既定・今の判定): 格子点の行 s−1 だけ。raw: 区間の毎 step の outer_end の全行 (欠け・重複は判定不能)。
SERIES_GRID = "grid"
SERIES_RAW = "raw"
# 壁の熱流の窓の量 (`wall_heat_nodes.measure`; §5.1 #5d、codex diagnose 2026-10-16 settle-v2-span10k)。
# terminal_diff (既定・今の判定): 終端 T の値からの差の最大。range: 節点ごとの窓の幅 max − min。どちらも許容の相対の参照は終端 T の値。
WALL_TERMINAL = "terminal_diff"
WALL_RANGE = "range"

KIND_RESIDUAL = "residual"
KIND_MASS = "mass_balance"
# 残差の非有限の見張り (判定の条件ではない; plan §4.4「残差を止める判定から外す」(2026-10-16 ユーザ決定)、§5.1 #8)。
# 残差 (residual) を外した基準 (プリセット wt_nozzle_v3) は代わりにこれを要る (発散を見張らずに止めない)。
KIND_RESIDUAL_FINITE = "residual_finite"
TARGET_KINDS = ("exit_core_mach", "line_mach", "exit_flow_angle", "wall_heat_nodes")
ALL_KINDS = (KIND_RESIDUAL, KIND_MASS) + TARGET_KINDS + (KIND_RESIDUAL_FINITE,)
NODE_KINDS = ("exit_core_mach", "line_mach", "exit_flow_angle")

# 種類ごとの閾値のキー (すべて必須; `provisional` は任意の記録)
THRESHOLD_KEYS = {
    "residual": ("slope_dec", "amp_dec", "rise_dec", "excursion_dec"),
    "mass_balance": ("mean_abs", "swing", "mdot_change_rel"),
    "exit_core_mach": ("tol",),
    "line_mach": ("tol_maxabs", "tol_halfamp", "tol_pos"),
    "exit_flow_angle": ("tol",),
    "wall_heat_nodes": ("q_tol", "eps_rel"),
    "residual_finite": (),          # 見張りだけ (閾値なし; thresholds は {} か provisional だけ)
}
# 種類ごとの許されるキー (`kind`・`name`・`thresholds`・`notes` は共通)
_COMMON_KEYS = ("kind", "name", "thresholds", "notes")
KIND_KEYS = {
    "residual": ("file", "required", "inactive_ok", "series"),
    "mass_balance": ("file", "series_csv", "inlet", "outlet"),
    "exit_core_mach": ("nodes_file", "series_csv", "eta_range"),
    "line_mach": ("nodes_file", "series_csv", "eta", "x_range_m", "n", "Md"),
    "exit_flow_angle": ("nodes_file", "eta_range"),
    "wall_heat_nodes": ("walls", "file_pattern", "dataset", "axisymmetric", "measure"),
    "residual_finite": ("file", "required"),
}
# residual_finite の必要な列の既定 (基準に `required` が無いとき; §5.1 #10)。forge の残差の履歴は保存量の列を常に書く
# (`main.cpp` の kResidualEquationNames: ro・roUx・roUy・roUz・roe)。名前と範囲は `residual` の既定の required
# (`settle_presets.residual_columns` の最初の 4 つ; rms_roUz は inactive_ok なので要らない) に揃える。
RESIDUAL_FINITE_DEFAULT_REQUIRED = ("rms_ro", "rms_roUx", "rms_roUy", "rms_roe")
SERIES_KEYS = {
    "exit_core_mach": (("path", "step_col", "value_col"), ()),
    "line_mach": (("path", "step_col", "maxabs_col", "halfamp_col"), ()),
    "mass_balance": (("path", "step_col", "value_col"), ("mdot_in_col", "mdot_out_col")),
}
TOP_KEYS = ("schema", "source", "structured", "grid", "window", "hold", "quantities", "excluded", "preconditions", "notes")
DEFAULT_NODES_FILE = "res_nodes_{step}.h5"
DEFAULT_WALL_PATTERN = "res_{name}_{physID}_{step}.h5"
DEFAULT_WALL_DATASET = "iface_Qf_eff"


class CriteriaError(ValueError):
    """settle_criteria.yaml の誤用 (知らない種類・キー、閾値の欠け、型の誤り、重複キー)。黙って既定値で埋めない。"""


class _DataMissing(Exception):
    """判定に要るデータが無い・読めない・整合しない (判定不能、分類 data)。"""


class _NonFinite(Exception):
    """入力に非有限がある (DIVERGED、分類 nonfinite)。"""


# ============================================================ 基準のファイル
def unset_message(run_dir, why: str) -> str:
    p = Path(run_dir) / CRITERIA_FILE
    return (f"INDETERMINATE (評価対象が未設定): {why}。置き場所: {p} (手で組む run はひな形 {TEMPLATE} を写す; 設計チェーンの run は"
            f"問題 YAML の evaluate.settle を準備のときに解決して写す)。量の種類: {', '.join(ALL_KINDS)} "
            f"(残差 residual か残差の非有限の見張り residual_finite、流量の収支 mass_balance と、目的の量 {', '.join(TARGET_KINDS)} の "
            "1 つ以上が要る)")


def _yaml_load_strict(text: str):
    """重複キーを例外にする safe_load (PyYAML は黙って後勝ちにする)。"""
    if yaml is None:
        raise CriteriaError("PyYAML が無い (settle_criteria.yaml を読めない)")

    class _Loader(yaml.SafeLoader):
        pass

    def _construct_mapping(loader, node, deep=False):
        keys = set()
        for k_node, _ in node.value:
            k = loader.construct_object(k_node, deep=deep)
            if k in keys:
                raise CriteriaError(f"重複キー '{k}' (行 {k_node.start_mark.line + 1})")
            keys.add(k)
        return yaml.SafeLoader.construct_mapping(loader, node, deep=deep)

    _Loader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_mapping)
    try:
        return yaml.load(text, Loader=_Loader)
    except CriteriaError:
        raise
    except yaml.YAMLError as e:
        raise CriteriaError(f"YAML として読めない: {e}") from e


def _num(v, where, *, integer=False, positive=False, nonneg=False):
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise CriteriaError(f"{where}: 数でない ({v!r})")
    if integer and not (isinstance(v, int) or float(v).is_integer()):
        raise CriteriaError(f"{where}: 整数でない ({v!r})")
    x = int(v) if integer else float(v)
    if not math.isfinite(float(x)):
        raise CriteriaError(f"{where}: 有限でない ({v!r})")
    if positive and not x > 0:
        raise CriteriaError(f"{where}: 正でない ({v!r})")
    if nonneg and not x >= 0:
        raise CriteriaError(f"{where}: 負 ({v!r})")
    return x


def _keys(d, allowed, where, required=()):
    if not isinstance(d, dict):
        raise CriteriaError(f"{where}: 対応表 (mapping) でない ({type(d).__name__})")
    unknown = [k for k in d if k not in allowed]
    if unknown:
        raise CriteriaError(f"{where}: 知らないキー {unknown} (許されるキー: {list(allowed)})")
    missing = [k for k in required if k not in d]
    if missing:
        raise CriteriaError(f"{where}: キーが欠けている {missing}")


def _str_list(v, where, *, allow_empty=True):
    if not isinstance(v, list) or not all(isinstance(x, str) and x for x in v):
        raise CriteriaError(f"{where}: 文字列のリストでない ({v!r})")
    if not allow_empty and not v:
        raise CriteriaError(f"{where}: 空")
    if len(set(v)) != len(v):
        raise CriteriaError(f"{where}: 重複がある ({v!r})")
    return list(v)


def _int_list(v, where):
    if not isinstance(v, list) or not v:
        raise CriteriaError(f"{where}: 整数の空でないリストでない ({v!r})")
    out = [_num(x, where, integer=True, nonneg=True) for x in v]
    if len(set(out)) != len(out):
        raise CriteriaError(f"{where}: 重複がある ({v!r})")
    return out


def _thresholds(kind, th, where):
    keys = THRESHOLD_KEYS[kind]
    _keys(th, keys + ("provisional",), where + ".thresholds", required=keys)
    out = {k: _num(th[k], f"{where}.thresholds.{k}", nonneg=True) for k in keys}
    if "provisional" in th:
        if not isinstance(th["provisional"], bool):
            raise CriteriaError(f"{where}.thresholds.provisional: true/false でない ({th['provisional']!r})")
        out["provisional"] = th["provisional"]
    return out


def _eta_range(v, where):
    if not (isinstance(v, list) and len(v) == 2):
        raise CriteriaError(f"{where}: [lo, hi] でない ({v!r})")
    lo, hi = (_num(x, where) for x in v)
    if not (0.0 <= lo < hi <= 1.0):
        raise CriteriaError(f"{where}: 0 ≤ lo < hi ≤ 1 でない ({v!r})")
    return [lo, hi]


def _series(kind, s, where):
    req, opt = SERIES_KEYS[kind]
    _keys(s, req + opt, where, required=req)
    for k in req + opt:
        if k in s and not (isinstance(s[k], str) and s[k]):
            raise CriteriaError(f"{where}.{k}: 空でない文字列でない ({s[k]!r})")
    if kind == "mass_balance" and (("mdot_in_col" in s) != ("mdot_out_col" in s)):
        raise CriteriaError(f"{where}: mdot_in_col と mdot_out_col は組で書く")
    return dict(s)


def validate_criteria(obj) -> dict:
    """基準の中身を検査して正規化した dict を返す。誤用は CriteriaError。未設定 (量が足りない) はここでは検査しない
    (`load_criteria` が INDETERMINATE にする)。"""
    _keys(obj, TOP_KEYS, "settle_criteria", required=("schema", "grid", "window"))
    if obj["schema"] != SCHEMA:
        raise CriteriaError(f"schema が '{obj['schema']}' (期待 '{SCHEMA}')")
    c = {"schema": SCHEMA}
    g = obj["grid"]
    _keys(g, ("out_step_interval", "start_step", "segment"), "grid", required=("out_step_interval", "start_step"))
    dt = _num(g["out_step_interval"], "grid.out_step_interval", integer=True, positive=True)
    st = _num(g["start_step"], "grid.start_step", integer=True, nonneg=True)
    c["grid"] = {"out_step_interval": dt, "start_step": st, "segment": g.get("segment")}
    w = obj["window"]
    _keys(w, ("W", "min_points"), "window", required=("W", "min_points"))
    W = _num(w["W"], "window.W", integer=True, positive=True)
    mp = _num(w["min_points"], "window.min_points", integer=True, positive=True)
    if W % dt != 0:
        raise CriteriaError(f"window.W {W} が grid.out_step_interval {dt} の倍数でない (窓の端が格子に乗らない)")
    if W // dt + 1 < mp:
        raise CriteriaError(f"窓 W {W} (間隔 {dt}) の端点込みの点 {W // dt + 1} が min_points {mp} に足りない")
    if mp < 3:
        raise CriteriaError(f"window.min_points {mp} < 3 (傾き・振幅を測れない)")
    c["window"] = {"W": W, "min_points": mp}
    if "hold" in obj:
        c["hold"] = _num(obj["hold"], "hold", integer=True, nonneg=True)
    for k in ("source", "notes"):
        if k in obj:
            c[k] = obj[k]
    if "structured" in obj:
        s = obj["structured"]
        _keys(s, ("ni", "nj", "scale_m", "mesh_file"), "structured", required=("ni", "nj", "scale_m"))
        c["structured"] = {"ni": _num(s["ni"], "structured.ni", integer=True, positive=True),
                           "nj": _num(s["nj"], "structured.nj", integer=True, positive=True),
                           "scale_m": _num(s["scale_m"], "structured.scale_m", positive=True)}
        if "mesh_file" in s:
            if not (isinstance(s["mesh_file"], str) and s["mesh_file"]):
                raise CriteriaError("structured.mesh_file: 空でない文字列でない")
            c["structured"]["mesh_file"] = s["mesh_file"]
        if c["structured"]["ni"] < 2 or c["structured"]["nj"] < 2:
            raise CriteriaError("structured.ni・nj は 2 以上")
    ex = obj.get("excluded", [])
    if not isinstance(ex, list):
        raise CriteriaError("excluded: リストでない")
    for i, e in enumerate(ex):
        _keys(e, ("kind", "name", "reason"), f"excluded[{i}]", required=("kind", "reason"))
        if e["kind"] not in ALL_KINDS:
            raise CriteriaError(f"excluded[{i}].kind: 知らない種類 '{e['kind']}' (種類: {list(ALL_KINDS)})")
        if not (isinstance(e["reason"], str) and e["reason"].strip()):
            raise CriteriaError(f"excluded[{i}].reason: 理由が空 (黙って外さない)")
    c["excluded"] = ex
    pc = obj.get("preconditions", [])
    if not isinstance(pc, list):
        raise CriteriaError("preconditions: リストでない")
    for i, p in enumerate(pc):
        _keys(p, ("quantity", "requires", "satisfied"), f"preconditions[{i}]", required=("quantity", "requires", "satisfied"))
    c["preconditions"] = pc
    qs = obj.get("quantities", [])
    if qs is None:
        qs = []
    if not isinstance(qs, list):
        raise CriteriaError("quantities: リストでない")
    names = set()
    out = []
    for i, q in enumerate(qs):
        where = f"quantities[{i}]"
        if not isinstance(q, dict) or "kind" not in q:
            raise CriteriaError(f"{where}: kind が無い")
        kind = q["kind"]
        if kind not in ALL_KINDS:
            raise CriteriaError(f"{where}: 知らない種類 '{kind}' (種類: {list(ALL_KINDS)})")
        _keys(q, _COMMON_KEYS + KIND_KEYS[kind], where, required=("thresholds",))
        name = q.get("name", kind)
        if not (isinstance(name, str) and name):
            raise CriteriaError(f"{where}.name: 空でない文字列でない")
        if name in names:
            raise CriteriaError(f"{where}: 名前 '{name}' が重複 (同じ種類を 2 つ書くときは name を付ける)")
        names.add(name)
        o = {"kind": kind, "name": name, "thresholds": _thresholds(kind, q["thresholds"], where)}
        if "notes" in q:
            o["notes"] = q["notes"]
        if kind == "residual":
            if "required" not in q:
                raise CriteriaError(f"{where}: required (必要な残差の列) が無い")
            o["file"] = q.get("file", "residual_history.csv")
            o["required"] = _str_list(q["required"], f"{where}.required", allow_empty=False)
            o["inactive_ok"] = _str_list(q.get("inactive_ok", []), f"{where}.inactive_ok")
            both = set(o["required"]) & set(o["inactive_ok"])
            if both:
                raise CriteriaError(f"{where}: required と inactive_ok の両方にある列 {sorted(both)}")
            for col in o["required"] + o["inactive_ok"]:
                if not col.startswith("rms_") or col.startswith("rms_dq_"):
                    raise CriteriaError(f"{where}: 残差の列 '{col}' は rms_<方程式> の形でない (rms_dq_* は判定に使わない)")
            if "series" in q:      # 書いたときだけ残す (書かない基準の中身は今と同じ)
                if q["series"] not in (SERIES_GRID, SERIES_RAW):
                    raise CriteriaError(f"{where}.series: '{q['series']}' は {SERIES_GRID} か {SERIES_RAW}")
                o["series"] = q["series"]
        elif kind == KIND_RESIDUAL_FINITE:
            o["file"] = q.get("file", "residual_history.csv")
            if "required" in q:    # 書いたときだけ残す (書かない基準の中身は今と同じ; 見張りは既定の 4 列を要る)
                o["required"] = _str_list(q["required"], f"{where}.required", allow_empty=False)
                for col in o["required"]:
                    if not col.startswith("rms_") or col.startswith("rms_dq_"):
                        raise CriteriaError(f"{where}: 残差の列 '{col}' は rms_<方程式> の形でない (rms_dq_* は必要な列にしない)")
        elif kind == "mass_balance":
            if "series_csv" in q:
                if "file" in q:
                    raise CriteriaError(f"{where}: series_csv と file は同時に書けない")
                o["series_csv"] = _series(kind, q["series_csv"], f"{where}.series_csv")
            else:
                o["file"] = q.get("file", "boundary_massflux.csv")
            for k in ("inlet", "outlet"):
                if k not in q:
                    raise CriteriaError(f"{where}: {k} (physID のリスト) が無い")
                o[k] = _int_list(q[k], f"{where}.{k}")
            if set(o["inlet"]) & set(o["outlet"]):
                raise CriteriaError(f"{where}: inlet と outlet の両方にある physID")
        elif kind in ("exit_core_mach", "exit_flow_angle"):
            if "eta_range" not in q:
                raise CriteriaError(f"{where}: eta_range が無い")
            o["eta_range"] = _eta_range(q["eta_range"], f"{where}.eta_range")
            if "series_csv" in q:
                if "nodes_file" in q:
                    raise CriteriaError(f"{where}: series_csv と nodes_file は同時に書けない")
                o["series_csv"] = _series(kind, q["series_csv"], f"{where}.series_csv")
            else:
                o["nodes_file"] = q.get("nodes_file", DEFAULT_NODES_FILE)
        elif kind == "line_mach":
            for k in ("eta", "x_range_m", "n", "Md"):
                if k not in q:
                    raise CriteriaError(f"{where}: {k} が無い")
            o["eta"] = _num(q["eta"], f"{where}.eta")
            if not 0.0 <= o["eta"] <= 1.0:
                raise CriteriaError(f"{where}.eta: [0, 1] の外 ({q['eta']!r})")
            xr = q["x_range_m"]
            if not (isinstance(xr, list) and len(xr) == 2):
                raise CriteriaError(f"{where}.x_range_m: [x0, x1] でない")
            o["x_range_m"] = [_num(xr[0], f"{where}.x_range_m"), _num(xr[1], f"{where}.x_range_m")]
            if not o["x_range_m"][0] < o["x_range_m"][1]:
                raise CriteriaError(f"{where}.x_range_m: x0 < x1 でない")
            o["n"] = _num(q["n"], f"{where}.n", integer=True, positive=True)
            if o["n"] < 2:
                raise CriteriaError(f"{where}.n: 2 以上")
            o["Md"] = _num(q["Md"], f"{where}.Md", positive=True)
            if "series_csv" in q:
                if "nodes_file" in q:
                    raise CriteriaError(f"{where}: series_csv と nodes_file は同時に書けない")
                o["series_csv"] = _series(kind, q["series_csv"], f"{where}.series_csv")
            else:
                o["nodes_file"] = q.get("nodes_file", DEFAULT_NODES_FILE)
        elif kind == "wall_heat_nodes":
            if "walls" not in q or "axisymmetric" not in q:
                raise CriteriaError(f"{where}: walls と axisymmetric が要る")
            walls = q["walls"]
            if not isinstance(walls, list) or not walls:
                raise CriteriaError(f"{where}.walls: 空でないリストでない")
            ws = []
            for k, wl in enumerate(walls):
                _keys(wl, ("name", "physID"), f"{where}.walls[{k}]", required=("name", "physID"))
                if not (isinstance(wl["name"], str) and wl["name"]):
                    raise CriteriaError(f"{where}.walls[{k}].name: 空でない文字列でない")
                ws.append({"name": wl["name"], "physID": _num(wl["physID"], f"{where}.walls[{k}].physID", integer=True, nonneg=True)})
            o["walls"] = ws
            if not isinstance(q["axisymmetric"], bool):
                raise CriteriaError(f"{where}.axisymmetric: true/false でない")
            o["axisymmetric"] = q["axisymmetric"]
            o["file_pattern"] = q.get("file_pattern", DEFAULT_WALL_PATTERN)
            o["dataset"] = q.get("dataset", DEFAULT_WALL_DATASET)
            if "measure" in q:     # 書いたときだけ残す (書かない基準の中身は今と同じ)
                if q["measure"] not in (WALL_TERMINAL, WALL_RANGE):
                    raise CriteriaError(f"{where}.measure: '{q['measure']}' は {WALL_TERMINAL} か {WALL_RANGE}")
                o["measure"] = q["measure"]
        if kind in NODE_KINDS and "series_csv" not in o and "structured" not in c:
            raise CriteriaError(f"{where}: 節点の量 ({kind}) には structured: {{ni, nj, scale_m}} が要る")
        out.append(o)
    c["quantities"] = out
    return c


def _sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def load_criteria(run_dir, path=None) -> dict:
    """run の settle_criteria.yaml を読む。path を渡すと run の settle_criteria.yaml の代わりにそのファイルを読む (再生で別の基準を
    当てるとき; `settle_virtual_stop.py --criteria`。基準の中の相対のパスは今どおり run から解く)。

    戻り: {'status': 'ok', 'criteria': <検証済み>, 'path', 'sha256'} か
          {'status': 'unset', 'verdict': 'INDETERMINATE', 'reason': ..., 'message': ...} (ファイルが無い・量が足りない)。
    誤用 (知らない種類・キー、閾値の欠けなど) は CriteriaError。"""
    p = Path(path) if path is not None else Path(run_dir) / CRITERIA_FILE
    if not p.is_file():
        why = f"{p} が無い"
        return {"status": "unset", "verdict": "INDETERMINATE", "reason": why, "message": unset_message(run_dir, why), "path": str(p)}
    raw = p.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as e:
        raise CriteriaError(f"{p}: UTF-8 で読めない ({e})") from e
    obj = _yaml_load_strict(text)
    if obj is None:
        why = f"{p} が空"
        return {"status": "unset", "verdict": "INDETERMINATE", "reason": why, "message": unset_message(run_dir, why), "path": str(p),
                "sha256": _sha256_bytes(raw)}
    crit = validate_criteria(obj)
    kinds = [q["kind"] for q in crit["quantities"]]
    lack = []
    if not any(k in TARGET_KINDS for k in kinds):
        lack.append("目的の量 (" + ", ".join(TARGET_KINDS) + ") が 1 つも無い")
    if KIND_RESIDUAL not in kinds and KIND_RESIDUAL_FINITE not in kinds:
        # 残差を止める判定から外した基準 (wt_nozzle_v3) でも、発散の見張り (residual_finite) は要る (plan §4.4、§5.1 #8)
        lack.append("残差 (residual) も残差の非有限の見張り (residual_finite) も無い (発散を見張れない)")
    if KIND_MASS not in kinds:
        lack.append("流量の収支 (mass_balance) が無い")
    if kinds.count(KIND_RESIDUAL) > 1:
        raise CriteriaError("residual が 2 つ以上ある")
    if kinds.count(KIND_RESIDUAL_FINITE) > 1:
        raise CriteriaError("residual_finite が 2 つ以上ある")
    if lack:
        why = f"{p}: " + "、".join(lack)
        return {"status": "unset", "verdict": "INDETERMINATE", "reason": why, "message": unset_message(run_dir, why), "path": str(p),
                "sha256": _sha256_bytes(raw), "criteria": crit}
    return {"status": "ok", "criteria": crit, "path": str(p), "sha256": _sha256_bytes(raw)}


def any_provisional(crit: dict) -> bool:
    return any(q["thresholds"].get("provisional") for q in crit.get("quantities", []))


# ============================================================ 入力のデータ (読み込みは 1 回だけ、--scan で使い回す)
def _open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r", newline="")


def _resolve(run: Path, p: str) -> Path:
    q = Path(p)
    return q if q.is_absolute() else run / q


class RunData:
    """run の入力のデータ。読んだファイルを `used` に残す (結果のファイルの一覧)。"""

    def __init__(self, run_dir, crit: dict):
        self.run = Path(run_dir)
        self.crit = crit
        self.used: dict[str, Path] = {}
        self._resid = {}
        self._mf = {}
        self._series = {}
        self._nodes = {}
        self._wall = {}
        self._mesh_coord = None

    def _use(self, p: Path):
        self.used[str(p)] = p

    # ---------- 残差
    def residual(self, fname: str) -> dict:
        if fname in self._resid:
            r = self._resid[fname]
            if isinstance(r, Exception):
                raise r
            return r
        p = _resolve(self.run, fname)
        try:
            if not p.is_file():
                raise _DataMissing(f"残差の履歴 {p} が無い")
            self._use(p)
            with _open_text(p) as f:
                rd = csv.reader(f)
                try:
                    h = next(rd)
                except StopIteration:
                    raise _DataMissing(f"{p} が空") from None
                if len(h) < 3 or h[0] != "step" or h[2] != "phase":
                    raise _DataMissing(f"{p}: 見出しが step,inner,phase,... でない ({h[:4]})")
                ic = [i for i, k in enumerate(h) if i >= 3]
                names = [h[i] for i in ic]
                steps, vals, bad = [], [], []
                for row in rd:
                    if len(row) < 3 or row[2] != "outer_end":
                        continue
                    try:
                        s = int(row[0])
                    except ValueError:
                        bad.append(None)
                        continue
                    if len(row) != len(h):
                        bad.append(s)
                        continue
                    try:
                        vals.append([float(row[i]) for i in ic])
                    except ValueError:
                        bad.append(s)
                        continue
                    steps.append(s)
            st = np.asarray(steps, dtype=np.int64)
            v = np.asarray(vals, dtype=float).reshape(len(steps), len(ic))
            r = {"path": str(p), "steps": st, "cols": {n: v[:, k] for k, n in enumerate(names)}, "bad": bad}
        except Exception as e:  # noqa: BLE001
            self._resid[fname] = e if isinstance(e, _DataMissing) else _DataMissing(f"{p}: 読めない ({e})")
            raise self._resid[fname] from None
        self._resid[fname] = r
        return r

    # ---------- 境界の質量流束
    def massflux(self, fname: str) -> dict:
        if fname in self._mf:
            r = self._mf[fname]
            if isinstance(r, Exception):
                raise r
            return r
        p = _resolve(self.run, fname)
        try:
            if not p.is_file():
                raise _DataMissing(f"境界の質量流束 {p} が無い (forge の output.boundaryMassflux: 1 が要る)")
            self._use(p)
            with _open_text(p) as f:
                lines = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
            if not lines:
                raise _DataMissing(f"{p}: 見出しが無い")
            rd = csv.reader(lines)
            h = next(rd)
            if h[0] != "step":
                raise _DataMissing(f"{p}: 見出しの最初が step でない ({h[0]})")
            mcols = {}
            for i, k in enumerate(h):
                m = re.fullmatch(r"mdot_(.+)_(\d+)", k)
                if m:
                    mcols[int(m.group(2))] = (m.group(1), i)
            steps, rows, bad = [], [], []
            for row in rd:
                try:
                    s = int(row[0])
                except (ValueError, IndexError):
                    bad.append(None)
                    continue
                if len(row) != len(h):
                    bad.append(s)
                    continue
                try:
                    rows.append([float(x) for x in row[1:]])
                except ValueError:
                    bad.append(s)
                    continue
                steps.append(s)
            v = np.asarray(rows, dtype=float).reshape(len(steps), len(h) - 1)
            r = {"path": str(p), "steps": np.asarray(steps, dtype=np.int64),
                 "mdot": {pid: (nm, v[:, i - 1]) for pid, (nm, i) in mcols.items()}, "all": v, "bad": bad}
        except Exception as e:  # noqa: BLE001
            self._mf[fname] = e if isinstance(e, _DataMissing) else _DataMissing(f"{p}: 読めない ({e})")
            raise self._mf[fname] from None
        self._mf[fname] = r
        return r

    # ---------- 系列の CSV (再生用)
    def series(self, spec: dict) -> dict:
        key = json.dumps(spec, sort_keys=True)
        if key in self._series:
            r = self._series[key]
            if isinstance(r, Exception):
                raise r
            return r
        p = _resolve(self.run, spec["path"])
        try:
            if not p.is_file():
                raise _DataMissing(f"系列の CSV {p} が無い")
            self._use(p)
            with _open_text(p) as f:
                lines = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
            rd = csv.DictReader(lines)
            cols = [spec[k] for k in spec if k.endswith("_col") and k != "step_col"]
            for c in [spec["step_col"]] + cols:
                if rd.fieldnames is None or c not in rd.fieldnames:
                    raise _DataMissing(f"{p}: 列 '{c}' が無い (列: {rd.fieldnames})")
            steps, vals, bad = [], {c: [] for c in cols}, []
            for row in rd:
                try:
                    s = int(float(row[spec["step_col"]]))
                    vv = {c: float(row[c]) for c in cols}
                except (TypeError, ValueError):
                    bad.append(row.get(spec["step_col"]))
                    continue
                steps.append(s)
                for c in cols:
                    vals[c].append(vv[c])
            r = {"path": str(p), "steps": np.asarray(steps, dtype=np.int64), "cols": {c: np.asarray(v, dtype=float) for c, v in vals.items()},
                 "bad": bad}
        except Exception as e:  # noqa: BLE001
            self._series[key] = e if isinstance(e, _DataMissing) else _DataMissing(f"{p}: 読めない ({e})")
            raise self._series[key] from None
        self._series[key] = r
        return r

    # ---------- 節点の状態
    def nodes(self, pattern: str, step: int) -> dict:
        key = (pattern, step)
        if key in self._nodes:
            r = self._nodes[key]
            if isinstance(r, Exception):
                raise r
            return r
        p = self.run / pattern.format(step=step)
        try:
            r = _read_nodes(p, step, self.crit["structured"], self.mesh_coord())
            self._use(p)
        except (_DataMissing, _NonFinite) as e:
            self._nodes[key] = e
            raise
        except Exception as e:  # noqa: BLE001
            self._nodes[key] = _DataMissing(f"{p}: 読めない ({e})")
            raise self._nodes[key] from None
        self._nodes[key] = r
        return r

    # ---------- 格子の座標 (structured.mesh_file; load_field と同じ COORD)
    def mesh_coord(self):
        mf = self.crit.get("structured", {}).get("mesh_file")
        if not mf:
            return None
        if self._mesh_coord is None:
            import h5py
            p = _resolve(self.run, mf)
            try:
                if not p.is_file():
                    raise _DataMissing(f"structured.mesh_file の格子 {p} が無い")
                with h5py.File(p, "r") as f:
                    c = f["/MESH/COORD"][:].reshape(-1, 3)
                st = self.crit["structured"]
                if c.shape[0] != st["ni"] * st["nj"]:
                    raise _DataMissing(f"{p}: 節点数 {c.shape[0]} が ni × nj = {st['ni'] * st['nj']} でない")
                self._use(p)
                self._mesh_coord = c
            except _DataMissing as e:
                self._mesh_coord = e
            except Exception as e:  # noqa: BLE001
                self._mesh_coord = _DataMissing(f"{p}: 読めない ({e})")
        if isinstance(self._mesh_coord, Exception):
            raise self._mesh_coord
        return self._mesh_coord

    # ---------- 壁
    def wall(self, pattern: str, dataset: str, wall: dict, step: int) -> dict:
        key = (pattern, dataset, wall["name"], wall["physID"], step)
        if key in self._wall:
            r = self._wall[key]
            if isinstance(r, Exception):
                raise r
            return r
        p = self.run / pattern.format(name=wall["name"], physID=wall["physID"], step=step)
        try:
            r = _read_wall(p, step, dataset)
            self._use(p)
        except _DataMissing as e:
            self._wall[key] = e
            raise
        except Exception as e:  # noqa: BLE001
            self._wall[key] = _DataMissing(f"{p}: 読めない ({e})")
            raise self._wall[key] from None
        self._wall[key] = r
        return r

    # ---------- データの最後の step (end_step の省略時)
    def last_step(self) -> tuple[int | None, list]:
        cands = []
        for q in self.crit["quantities"]:
            try:
                if q["kind"] == "residual":
                    r = self.residual(q["file"])
                    if len(r["steps"]):
                        cands.append((int(r["steps"].max()) + 1, f"{q['name']}: 残差の最後の行 + 1"))
                elif "series_csv" in q:
                    r = self.series(q["series_csv"])
                    if len(r["steps"]):
                        cands.append((int(r["steps"].max()), f"{q['name']}: 系列の CSV の最後"))
                elif q["kind"] == "mass_balance":
                    r = self.massflux(q["file"])
                    if len(r["steps"]):
                        cands.append((int(r["steps"].max()), f"{q['name']}: 流量の CSV の最後"))
                elif q["kind"] in NODE_KINDS:
                    s = _last_file_step(self.run, q["nodes_file"])
                    if s is not None:
                        cands.append((s, f"{q['name']}: 節点の h5 の最後"))
                elif q["kind"] == "wall_heat_nodes":
                    for wl in q["walls"]:
                        s = _last_file_step(self.run, q["file_pattern"].replace("{name}", wl["name"]).replace("{physID}", str(wl["physID"])))
                        if s is not None:
                            cands.append((s, f"{q['name']}: 壁の h5 の最後"))
            except _DataMissing:
                continue
        if not cands:
            return None, []
        return max(c[0] for c in cands), cands


def _last_file_step(run: Path, pattern: str) -> int | None:
    pat = re.escape(pattern).replace(re.escape("{step}"), r"(\d+)")
    best = None
    for f in run.iterdir() if run.is_dir() else []:
        m = re.fullmatch(pat, f.name)
        if m:
            s = int(m.group(1))
            best = s if best is None else max(best, s)
    return best


def _read_nodes(p: Path, step: int, st: dict, mesh_coord=None) -> dict:
    """res_nodes_<step>.h5 を読む (節点の番号・座標・値)。座標・M の作り方は load_field と同じ順 (docstring)。
    mesh_coord (structured.mesh_file の /MESH/COORD) があれば座標はそれを使い、res_nodes の COORD がそれを res_nodes の型に
    丸めた値と一致することを確かめる (float のビルドは節点の座標を float32 で書くが、load_field は格子の h5 の double を読む)。"""
    import h5py
    if not p.is_file():
        raise _DataMissing(f"節点の出力 {p.name} が無い (forge の output.nodeProbeFile が要る)")
    ni, nj, S = st["ni"], st["nj"], st["scale_m"]
    with h5py.File(p, "r") as f:
        if "step" in f.attrs and int(f.attrs["step"]) != step:
            raise _DataMissing(f"{p.name}: 属性 step {int(f.attrs['step'])} がファイル名の step {step} と違う")
        ids = np.asarray(f["/MESH/NODE_ID"][:], dtype=np.int64)
        coord = f["/MESH/COORD"][:].reshape(-1, 3)
        need = ("Ux", "Uy", "sonic")
        for k in need:
            if "/VALUE/" + k not in f:
                raise _DataMissing(f"{p.name}: VALUE/{k} が無い")
        V = {k: f["/VALUE/" + k][:] for k in f["/VALUE"].keys()}
        probe_file = f.attrs.get("nodeProbeFile")
    if coord.shape[0] != ids.shape[0] or any(V[k].shape[0] != ids.shape[0] for k in V):
        raise _DataMissing(f"{p.name}: NODE_ID・COORD・VALUE の長さが揃わない")
    if ids.size == 0 or ids.min() < 0 or ids.max() >= ni * nj or len(np.unique(ids)) != ids.size:
        raise _DataMissing(f"{p.name}: 節点の番号が範囲外 (0〜{ni * nj - 1}) か重複")
    nonfin = {k: int(np.count_nonzero(~np.isfinite(a))) for k, a in V.items() if a.dtype.kind == "f"}
    nonfin = {k: n for k, n in nonfin.items() if n}
    coord_src = "res_nodes"
    if mesh_coord is not None:
        mc = mesh_coord[ids]
        if not np.array_equal(coord, mc.astype(coord.dtype)):
            raise _DataMissing(f"{p.name}: 節点の座標が structured.mesh_file の COORD と合わない (別の格子?)")
        coord = mc
        coord_src = "mesh_file"
    # load_field と同じ: X = COORD[:, 0] / S (COORD の型のまま)、M = hypot(Ux, Uy) / sonic (float64 にしてから)
    X = coord[:, 0] / S
    R = coord[:, 1] / S
    ux = V["Ux"].astype(float)
    uy = V["Uy"].astype(float)
    a = V["sonic"].astype(float)
    M = np.hypot(ux, uy) / a
    ci, cj = np.divmod(ids, nj)
    cols = {}
    order = np.lexsort((cj, ci))
    for k in order:
        cols.setdefault(int(ci[k]), []).append(int(k))
    cols = {i: np.asarray(ix, dtype=np.int64) for i, ix in cols.items()}
    return {"path": str(p), "step": step, "ids": ids, "j": cj, "cols": cols, "X": X, "R": R, "M": M, "ux": ux, "uy": uy,
            "nonfinite": nonfin, "coord_dtype": str(coord.dtype), "coord_source": coord_src, "nodeProbeFile": probe_file}


def _read_wall(p: Path, step: int, dataset: str) -> dict:
    import h5py
    if not p.is_file():
        raise _DataMissing(f"壁の出力 {p.name} が無い (bcond の outputHDFflg: 1 が要る)")
    with h5py.File(p, "r") as f:
        if "step" in f.attrs and int(f.attrs["step"]) != step:
            raise _DataMissing(f"{p.name}: 属性 step {int(f.attrs['step'])} がファイル名の step {step} と違う")
        if "/VALUE/" + dataset not in f:
            raise _DataMissing(f"{p.name}: VALUE/{dataset} が無い (forge の output.interfaceDiag: 1 が要る; 等温壁・定常のときだけ値を持つ)")
        Q = np.asarray(f["/VALUE/" + dataset][:], dtype=float)
        coord = np.asarray(f["/MESH/COORD"][:]).reshape(-1, 3)
        conne = np.asarray(f["/MESH/CONNE"][:], dtype=np.int64)
    if Q.shape[0] != coord.shape[0]:
        raise _DataMissing(f"{p.name}: 値の数 {Q.shape[0]} が節点の数 {coord.shape[0]} と違う (node の壁の出力でない)")
    return {"path": str(p), "Q": Q, "coord": coord, "conne": conne}


def wall_node_area(coord: np.ndarray, conne: np.ndarray, axisymmetric: bool) -> tuple[np.ndarray, str]:
    """壁の節点の面積 A_i (docstring の「A_i の取り方」)。戻り (A, 説明)。"""
    c = np.asarray(coord, dtype=float)
    n = c.shape[0]
    A = np.zeros(n)
    k = 0
    nseg = nface = 0
    L = len(conne)
    while k < L:
        t = int(conne[k])
        if t in (1, 2):                       # Polyvertex / Polyline: [type, nn, nodes...]
            nn = int(conne[k + 1])
            nodes = conne[k + 2:k + 2 + nn]
            k += 2 + nn
            if t == 2:
                for a, b in zip(nodes[:-1], nodes[1:]):
                    a, b = int(a), int(b)
                    le = float(np.linalg.norm(c[b] - c[a]))
                    if axisymmetric:
                        ra, rb = c[a, 1], c[b, 1]
                        A[a] += 0.5 * le * (3.0 * ra + rb) / 4.0
                        A[b] += 0.5 * le * (3.0 * rb + ra) / 4.0
                    else:
                        A[a] += 0.5 * le
                        A[b] += 0.5 * le
                    nseg += 1
        elif t in (4, 5):                     # Triangle / Quadrilateral: [type, nodes...]
            nn = 3 if t == 4 else 4
            nodes = [int(x) for x in conne[k + 1:k + 1 + nn]]
            k += 1 + nn
            if axisymmetric:
                raise _DataMissing("軸対称の壁の出力に三角形・四角形がある (2D の Polyline の壁だけに対応)")
            p0 = c[nodes[0]]
            area = 0.0
            for q in range(1, nn - 1):
                area += 0.5 * float(np.linalg.norm(np.cross(c[nodes[q]] - p0, c[nodes[q + 1]] - p0)))
            for a in nodes:
                A[a] += area / nn
            nface += 1
        else:
            raise _DataMissing(f"壁の出力の CONNE に知らない要素の型 {t}")
    if axisymmetric:
        how = (f"壁の h5 の COORD・CONNE (Polyline {nseg} 本) から A_i = Σ_e (|e|/2)(3 r_i + r_j)/4 [m²/rad] (1 rad あたり、r = y; "
               "iface_Qf_eff [W/rad] と単位を揃える)")
    else:
        how = f"壁の h5 の COORD・CONNE (Polyline {nseg} 本、面 {nface} 枚) から A_i = Σ (辺の長さ/2 または 面積/節点数)"
    return A, how


# ============================================================ 1 つの量の結果
class _Q:
    def __init__(self, q: dict, input_kind: str = "native"):
        self.r = {"name": q["name"], "kind": q["kind"], "input": input_kind, "thresholds": q["thresholds"],
                  "provisional": bool(q["thresholds"].get("provisional")), "checks": [], "reasons": [], "values": {}}

    def reason(self, cat: str, text: str):
        self.r["reasons"].append({"category": cat, "label": CAT_LABEL[cat], "text": text})

    def check(self, name: str, value, limit, ok: bool, *, detail: str | None = None):
        c = {"check": name, "value": value, "limit": limit, "ok": bool(ok)}
        if detail:
            c["detail"] = detail
        self.r["checks"].append(c)
        if not ok:
            self.reason(CAT_THRESHOLD, f"{name}: {_fmt(value)} > {_fmt(limit)}" + (f" ({detail})" if detail else ""))

    def finish(self) -> dict:
        cats = {x["category"] for x in self.r["reasons"]}
        if CAT_NONFINITE in cats:
            v = "DIVERGED"
        elif CAT_DATA in cats or CAT_UNSET in cats or CAT_GRID in cats:
            v = "INDETERMINATE"
        elif CAT_THRESHOLD in cats:
            v = "NOT_SETTLED"
        else:
            v = "SETTLED"
        self.r["verdict"] = v
        self.r["categories"] = sorted(cats)
        return self.r


def _fmt(v):
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def _span(a) -> float:
    a = np.asarray(a, dtype=float)
    return float(a.max() - a.min())


def _check_points(Q: _Q, pts, start: int, what: str):
    before = [s for s in pts if s < start]
    if before:
        Q.reason(CAT_DATA, f"{what}: 格子点 {before[0]}〜{before[-1]} が判定区間の始まり (step {start}) より前 — 履歴が窓を覆わない")
        return False
    return True


# ============================================================ 残差
def _raw_intervals(Q: _Q, st, sel, T: int, W: int, start: int, dt: int, ok_grid: bool, P):
    """残差の毎 step の系列 (series raw) の区間の行 (`_eval_residual` のコメント)。戻り (rawC, rawP, missC, missP):
    rawC・rawP は区間の行の位置 (判定できなければ None)、missC・missP は欠けた行の step (先頭 6 つ)。
    ok_grid が偽 (窓の格子点が判定区間の始まりより前; 理由は _check_points が書いた) なら何もしない。"""
    if not ok_grid:
        return None, None, [], []
    if P[0] <= start:
        Q.reason(CAT_DATA, f"直前の区間 (行 {T - 2 * W}〜{T - W - 1}) が本段の開始 (step {start}) の直後から始まる — 開始直後の行は直前の区間に"
                           f"入れない (§5.1 #2; grid と同じく最初に判定できる終端は開始 + 2W + 間隔 = {start + 2 * W + dt})")
        return None, None, [], []
    out = []
    ok = True
    for lab, lo, hi in (("区間", T - W, T - 1), ("直前の区間", T - 2 * W, T - W - 1)):
        rows = sel[(st[sel] >= lo) & (st[sel] <= hi)]
        s = st[rows]
        u = np.unique(s)
        ndup = int(len(s) - len(u))
        nmiss = int((hi - lo + 1) - len(u))
        miss = [int(x) for x in np.setdiff1d(np.arange(lo, hi + 1, dtype=np.int64), u, assume_unique=True)[:6]] if nmiss else []
        if nmiss or ndup:
            ok = False
            Q.reason(CAT_DATA, f"毎 step の系列の{lab} (行 {lo}〜{hi}、{hi - lo + 1} 行が要る) の outer_end の行の欠け {nmiss} 行"
                               + (f" (例 step {miss}{'…' if nmiss > len(miss) else ''})" if nmiss else "")
                               + (f"・重複 {ndup} 行" if ndup else "") + " — 全行がそろわないと傾き・振幅・上昇を測れない")
        out.append((rows, miss))
    if not ok:
        return None, None, out[0][1], out[1][1]
    return out[0][0], out[1][0], [], []


def _eval_residual(q: dict, ctx: dict) -> dict:
    Q = _Q(q)
    th = q["thresholds"]
    W, T, C, P, start, stop = ctx["W"], ctx["T"], ctx["C"], ctx["P"], ctx["start"], ctx["stop"]
    try:
        R = ctx["data"].residual(q["file"])
    except _DataMissing as e:
        Q.reason(CAT_DATA, str(e))
        return Q.finish()
    st = R["steps"]
    end_row = (stop if stop is not None else T) - 1
    seg = (st >= start) & (st <= end_row)
    s_seg = st[seg]
    Q.r["values"]["file"] = R["path"]
    Q.r["values"]["rows_in_segment"] = int(seg.sum())
    bad_in = [b for b in R["bad"] if b is None or start <= b <= end_row]
    if bad_in:
        Q.reason(CAT_DATA, f"読めない・欠けた outer_end の行が {len(bad_in)} 行 (step {bad_in[:5]})")
    if len(s_seg) > 1:
        d = np.diff(s_seg)
        if np.any(d == 0):
            Q.reason(CAT_DATA, f"step の重複 (例 {int(s_seg[1:][d == 0][0])})")
        if np.any(d < 0):
            Q.reason(CAT_DATA, f"step が単調でない (例 {int(s_seg[1:][d < 0][0])} の前が {int(s_seg[:-1][d < 0][0])})")
    # 非有限: 区間の全数値列
    nonfin = {k: int(np.count_nonzero(~np.isfinite(v[seg]))) for k, v in R["cols"].items()}
    nonfin = {k: n for k, n in nonfin.items() if n}
    if nonfin:
        first = {}
        for k in nonfin:
            v = R["cols"][k][seg]
            first[k] = int(s_seg[np.argmax(~np.isfinite(v))])
        Q.reason(CAT_NONFINITE, f"非有限の残差: {nonfin} (最初の行の step {first})")
        Q.r["values"]["nonfinite"] = {"count": nonfin, "first_step": first}
        return Q.finish()
    rms = [k for k in R["cols"] if k.startswith("rms_") and not k.startswith("rms_dq_")]
    listed = set(q["required"]) | set(q["inactive_ok"])
    missing = [k for k in q["required"] if k not in R["cols"]]
    if missing:
        Q.reason(CAT_DATA, f"必要な残差の列が無い: {missing} (ある列: {rms})")
    judged = []
    col_res = {}
    for k in rms:
        v = R["cols"][k][seg]
        nz = np.count_nonzero(v == 0.0)
        neg = np.count_nonzero(v < 0.0)
        if k not in listed:
            if len(v) and nz == len(v):
                col_res[k] = {"status": "inactive (基準に無い・恒等ゼロ)"}
                continue
            Q.reason(CAT_UNSET, f"残差の列 {k} は required にも inactive_ok にも無いが値を持つ (基準に書く)")
            col_res[k] = {"status": "unlisted"}
            continue
        if len(v) == 0:
            continue
        if nz == len(v):
            if k in q["inactive_ok"]:
                col_res[k] = {"status": "inactive (恒等ゼロ、inactive_ok)"}
                continue
            Q.reason(CAT_DATA, f"必要な列 {k} が恒等ゼロ (活動していない; 基準の required と run の設定が合わない)")
            col_res[k] = {"status": "zero"}
            continue
        if nz or neg:
            Q.reason(CAT_DATA, f"列 {k} の一部がゼロ・負 ({nz + neg} 行; 例 step {int(s_seg[np.argmax(v <= 0.0)])})")
            col_res[k] = {"status": "partly zero"}
            continue
        judged.append(k)
    if not len(s_seg):
        Q.reason(CAT_DATA, f"判定区間 (step {start}〜{end_row}) の outer_end の行が無い")
    if not judged and not Q.r["reasons"]:
        Q.reason(CAT_DATA, "判定できる残差の列が無い")
    okC = _check_points(Q, C, start, "窓 W")
    okP = _check_points(Q, P, start, "直前の窓 (上昇の検査)")
    idx = {int(s): i for i, s in enumerate(st)}
    sel = np.flatnonzero(seg)
    lo_raw, hi_raw = T - W - 1, T - 1
    raw_w = sel[(st[sel] >= lo_raw) & (st[sel] <= hi_raw)]
    rowsC = [idx.get(s - 1) for s in C]
    rowsP = [idx.get(s - 1) for s in P]
    raw = q.get("series", SERIES_GRID) == SERIES_RAW
    rawC = rawP = None
    if raw:
        # 毎 step の系列 (§4.4・§5.1 #5d): 区間 (T−W, T] と直前の区間 (T−2W, T−W] の場の step を、行 s−1 の対応 (§4.8 の 1) で
        # 残差の行 T−W〜T−1 と T−2W〜T−W−1 (各 W 行) にし、その全行を使う。被覆 (整数 step が重複なく全部) を区間ごとに検査し、
        # 欠け・重複は判定不能 (codex diagnose 2026-10-16 settle-v2-span10k: 格子点の行だけの検査では途中の 1 行の欠け・2,000 ごとの
        # 行だけの履歴も SETTLED になる)。直前の区間が本段の開始の直後から始まる終端 (T − 2W = 開始) も判定不能にする: 開始直後の
        # 行を直前の区間に入れない (§5.1 #2 の決定) ので、grid と同じく最初に判定できる終端は開始 + 2W + 間隔。
        rawC, rawP, missC, missP = _raw_intervals(Q, st, sel, T, W, start, ctx["dt"], okC and okP, P)
    else:
        missC = [s for s, r in zip(C, rowsC) if r is None]
        missP = [s for s, r in zip(P, rowsP) if r is None]
        if okC and missC:
            Q.reason(CAT_DATA, f"窓の格子点 {missC[:6]}{'…' if len(missC) > 6 else ''} の残差の行 (s−1) が無い ({len(missC)}/{len(C)} 点の欠測)")
        if okP and missP:
            Q.reason(CAT_DATA, f"直前の窓の格子点 {missP[:6]}{'…' if len(missP) > 6 else ''} の残差の行 (s−1) が無い ({len(missP)}/{len(P)} 点)"
                               + (" — 本段の開始の場には残差の行が無い (最初に判定できる終端は開始 + 2W + 間隔)" if missP[0] == start else ""))
    stop_rows = None
    if stop is not None:
        stop_rows = sel[(st[sel] >= T) & (st[sel] <= stop - 1)]
        if idx.get(stop - 1) is None:
            Q.reason(CAT_DATA, f"停止の step {stop} の残差の行 {stop - 1} が無い")
            stop_rows = None
    Q.r["values"]["columns"] = col_res
    if raw:
        Q.r["values"]["series"] = SERIES_RAW
        Q.r["values"]["raw_rows"] = {"C": [T - W, T - 1], "P": [T - 2 * W, T - W - 1],
                                     "n_C": None if rawC is None else int(len(rawC)), "n_P": None if rawP is None else int(len(rawP))}
        have = rawC is not None and rawP is not None
        xs = st[rawC].astype(float) if have else None
        rdet = f"毎 step の系列 (行 {T - W}〜{T - 1} の {W} 行、直前の区間は行 {T - 2 * W}〜{T - W - 1})"
    else:
        have = okC and okP and not missC and not missP
        xs = np.asarray(C, dtype=float)
        rdet = None
    for k in judged:
        v = R["cols"][k]
        cr = {"status": "judged", "last": float(v[sel[-1]]) if len(sel) else None}
        if have:
            yC = np.log10(v[rawC if raw else rowsC])
            yP = np.log10(v[rawP if raw else rowsP])
            xc = xs - xs.mean()
            slope = float(np.dot(xc, yC - yC.mean()) / np.dot(xc, xc))
            amp = float(np.percentile(yC, 95) - np.percentile(yC, 5))
            rise = float(yC.mean() - yP.mean())
            if raw:
                cr.update(slope_dec_per_W=slope * W, amp_dec=amp, rise_dec=rise, n_C=int(len(yC)), n_P=int(len(yP)), series=SERIES_RAW)
            else:
                cr.update(slope_dec_per_W=slope * W, amp_dec=amp, rise_dec=rise, n_C=len(C), n_P=len(P))
            Q.check(f"{k} 傾き |s|·W [桁]", abs(slope) * W, th["slope_dec"], abs(slope) * W <= th["slope_dec"], detail=rdet)
            Q.check(f"{k} 振幅 p95−p5 [桁]", amp, th["amp_dec"], amp <= th["amp_dec"], detail=rdet)
            Q.check(f"{k} 直前の窓からの上昇 [桁]", rise, th["rise_dec"], rise <= th["rise_dec"], detail=rdet)
        if len(raw_w):
            yr = np.log10(v[raw_w])
            med = float(np.median(yr))
            exc = float(yr.max() - med)
            cr.update(excursion_dec=exc, n_raw=int(len(raw_w)), raw_max_step=int(st[raw_w][np.argmax(yr)]))
            Q.check(f"{k} 生の系列の跳ね max−median [桁]", exc, th["excursion_dec"], exc <= th["excursion_dec"],
                    detail=f"step {int(st[raw_w][np.argmax(yr)])}")
            if stop_rows is not None and len(stop_rows) and have:
                ys = np.log10(v[stop_rows])
                y_stop = float(np.log10(v[idx[stop - 1]]))
                jump = abs(y_stop - float(np.log10(v[rowsC[-1]])))
                exc_s = float(ys.max() - med)
                cr.update(stop_jump_dec=jump, stop_excursion_dec=exc_s)
                Q.check(f"{k} 停止の観測 (step {stop}) の格子の終端からの差 [桁]", jump, th["excursion_dec"], jump <= th["excursion_dec"])
                Q.check(f"{k} 停止の観測までの跳ね max−median [桁]", exc_s, th["excursion_dec"], exc_s <= th["excursion_dec"])
        elif okC:
            Q.reason(CAT_DATA, f"生の系列の窓 (行 {lo_raw}〜{hi_raw}) に行が無い")
        col_res[k] = cr
    return Q.finish()


def _eval_residual_finite(q: dict, ctx: dict) -> dict:
    """残差の非有限の見張り (`residual_finite`; plan §4.4「残差を止める判定から外す」(2026-10-16 ユーザ決定)、§5.1 #8)。
    **判定の条件ではない** (傾き・振幅・上昇・跳ねは見ない): 残差の履歴の本段の区間 (outer_end の行 start〜終端 − 1、停止の観測が
    あれば停止の step − 1 まで) のどの数値の列にも非有限が無ければ SETTLED (見張りの合格)、非有限 → DIVERGED (`residual` と同じ検査)。
    履歴が無い・読めない行・区間の行が無い・最後の行 (終端 − 1) が無い (見張りが終端まで届かない) → 判定不能。
    **被覆と列** (§5.1 #10、result 段レビュー 2 回目 M2): 区間の整数 step の行が重複なく全部・逆行なし、必要な列 (基準の `required`、
    無ければ `RESIDUAL_FINITE_DEFAULT_REQUIRED`) が全部ある、を必須にする (欠け・重複・逆行・列の欠け → 判定不能)。閾値は見ない。"""
    Q = _Q(q)
    T, start, stop = ctx["T"], ctx["start"], ctx["stop"]
    try:
        R = ctx["data"].residual(q["file"])
    except _DataMissing as e:
        Q.reason(CAT_DATA, str(e))
        return Q.finish()
    st = R["steps"]
    end_row = (stop if stop is not None else T) - 1
    seg = (st >= start) & (st <= end_row)
    s_seg = st[seg]
    req = list(q.get("required") or RESIDUAL_FINITE_DEFAULT_REQUIRED)
    Q.r["values"].update(file=R["path"], rows_in_segment=int(seg.sum()), rows=[start, end_row],
                         role="非有限の見張り (判定の条件ではない)", required=req,
                         required_source=("基準の required" if "required" in q
                                          else "既定 (forge が常に書く保存量の列; 乱流・化学種などの列は基準の required に書く)"))
    bad_in = [b for b in R["bad"] if b is None or start <= b <= end_row]
    if bad_in:
        Q.reason(CAT_DATA, f"読めない・欠けた outer_end の行が {len(bad_in)} 行 (step {bad_in[:5]})")
    if len(s_seg) > 1:         # `residual` と同じ検査 (重複・逆行は判定不能)
        d = np.diff(s_seg)
        if np.any(d == 0):
            Q.reason(CAT_DATA, f"step の重複 (例 {int(s_seg[1:][d == 0][0])})")
        if np.any(d < 0):
            Q.reason(CAT_DATA, f"step が単調でない (例 {int(s_seg[1:][d < 0][0])} の前が {int(s_seg[:-1][d < 0][0])})")
    nonfin = {k: int(np.count_nonzero(~np.isfinite(v[seg]))) for k, v in R["cols"].items()}
    nonfin = {k: n for k, n in nonfin.items() if n}
    if nonfin:
        first = {k: int(s_seg[np.argmax(~np.isfinite(R["cols"][k][seg]))]) for k in nonfin}
        Q.reason(CAT_NONFINITE, f"非有限の残差: {nonfin} (最初の行の step {first})")
        Q.r["values"]["nonfinite"] = {"count": nonfin, "first_step": first}
        return Q.finish()
    if not len(s_seg):
        Q.reason(CAT_DATA, f"区間 (行 {start}〜{end_row}) の outer_end の行が無い (見張れない)")
    else:
        if int(s_seg.max()) < end_row:
            Q.reason(CAT_DATA, f"区間の最後の行 (step {end_row}) が無い (最後の行 {int(s_seg.max())}; 見張りが終端まで届かない)")
        # 区間の毎 step の行 (forge は outer_end の行を毎 step 書く; `residual` の series raw の被覆と同じ定義)
        miss = np.setdiff1d(np.arange(start, end_row + 1, dtype=np.int64), s_seg)
        if len(miss):
            Q.reason(CAT_DATA, f"区間 (行 {start}〜{end_row}) の outer_end の行が {len(miss)} 行欠ける (step {miss[:6].tolist()}"
                               f"{' …' if len(miss) > 6 else ''}; 見張りは区間の毎 step の行を要る)")
    rms = [k for k in R["cols"] if k.startswith("rms_") and not k.startswith("rms_dq_")]
    miss_cols = [k for k in req if k not in R["cols"]]
    if miss_cols:
        Q.reason(CAT_DATA, f"必要な残差の列が無い: {miss_cols} (ある列: {rms})")
    return Q.finish()


# ============================================================ 流量の収支
def _eval_mass(q: dict, ctx: dict) -> dict:
    series = "series_csv" in q
    Q = _Q(q, "series_csv" if series else "native")
    th = q["thresholds"]
    C, T, start, stop = ctx["C"], ctx["T"], ctx["start"], ctx["stop"]
    if not _check_points(Q, C, start, "窓 W"):
        return Q.finish()
    pts = C + ([stop] if stop is not None else [])
    try:
        if series:
            S = ctx["data"].series(q["series_csv"])
            sp = q["series_csv"]
            steps = S["steps"]
            r_all = S["cols"][sp["value_col"]]
            mi_all = S["cols"].get(sp.get("mdot_in_col")) if "mdot_in_col" in sp else None
            mo_all = S["cols"].get(sp.get("mdot_out_col")) if "mdot_out_col" in sp else None
            others = {}
            bad = S["bad"]
        else:
            M = ctx["data"].massflux(q["file"])
            steps = M["steps"]
            for pid in q["inlet"] + q["outlet"]:
                if pid not in M["mdot"]:
                    raise _DataMissing(f"{M['path']}: physID {pid} の列 mdot_<name>_{pid} が無い (列: {sorted(M['mdot'])})")
            mi_all = -sum(M["mdot"][pid][1] for pid in q["inlet"])
            mo_all = sum(M["mdot"][pid][1] for pid in q["outlet"])
            with np.errstate(divide="ignore", invalid="ignore"):
                r_all = (mi_all - mo_all) / mi_all
            others = {f"{nm}_{pid}": arr for pid, (nm, arr) in M["mdot"].items() if pid not in q["inlet"] + q["outlet"]}
            bad = M["bad"]
    except _DataMissing as e:
        Q.reason(CAT_DATA, str(e))
        return Q.finish()
    endp = stop if stop is not None else T
    seg = (steps >= start) & (steps <= endp)
    s_seg = steps[seg]
    if [b for b in bad if b is None or (isinstance(b, int) and start <= b <= endp)]:
        Q.reason(CAT_DATA, "読めない・欠けた行がある")
    if len(s_seg) > 1 and np.any(np.diff(s_seg) <= 0):
        Q.reason(CAT_DATA, "step が単調でない・重複")
    if not series:
        allv = ctx["data"].massflux(q["file"])["all"][seg]
        if np.any(~np.isfinite(allv)):
            Q.reason(CAT_NONFINITE, f"非有限の流量 (step {int(s_seg[np.argmax(np.any(~np.isfinite(allv), axis=1))])})")
            return Q.finish()
    elif np.any(~np.isfinite(r_all[seg])) or any(a is not None and np.any(~np.isfinite(a[seg])) for a in (mi_all, mo_all)):
        Q.reason(CAT_NONFINITE, "系列の CSV の r・ṁ に非有限")
        return Q.finish()
    idx = {int(s): i for i, s in enumerate(steps)}
    miss = [s for s in pts if s not in idx]
    if miss:
        Q.reason(CAT_DATA, f"格子点 {miss[:6]}{'…' if len(miss) > 6 else ''} の流量の行が無い ({len(miss)}/{len(pts)})")
        return Q.finish()
    rows = [idx[s] for s in C]
    r = r_all[rows]
    vals = {"r": [float(x) for x in r], "steps": C}
    # 分母の成立 (result 段レビュー M3: 出口の流量 0 で ZeroDivisionError、オンラインでは見張りの例外で forge を止める経路だった)。
    #   ṁ_in ≤ 0: r = (ṁ_in − ṁ_out)/ṁ_in の分母が成り立たない (流入が無い・入口の逆流・physID の取り違え) → 判定不能 (データ不足)。
    #   ṁ_in > 0 で ṁ_out ≤ 0 (流出が無い・逆流・符号の反転): r = 1 − ṁ_out/ṁ_in ≥ 1 が測れて収支が取れていないことは確か
    #     → 理由つきの不合格 (NOT_SETTLED)。ṁ_out の相対の変化 (max − min)/mean は定義できないので計算しない。
    out_ok = True
    if mi_all is not None:
        mi = mi_all[rows]
        mo = mo_all[rows]
        if np.any(mi <= 0):
            bad_s = [C[i] for i in np.flatnonzero(mi <= 0)]
            Q.reason(CAT_DATA, f"入口の流量 ṁ_in が正でない格子点がある (step {bad_s[:6]}{'…' if len(bad_s) > 6 else ''}; 流入が無い — "
                               "r の分母が成り立たない。inlet・outlet の physID か流れを確かめる)")
            return Q.finish()
        vals.update(mdot_in=[float(x) for x in mi], mdot_out=[float(x) for x in mo])
        out_ok = _mdot_out_positive(Q, mo, C, "窓の")
    Q.r["values"] = vals
    if not series and others:
        Q.r["values"]["others_at_T"] = {k: float(a[idx[T]]) for k, a in others.items()}
    mr = float(r.mean())
    Q.check("|mean r|", abs(mr), th["mean_abs"], abs(mr) <= th["mean_abs"])
    Q.check("max r − min r", _span(r), th["swing"], _span(r) <= th["swing"])
    if mi_all is not None:
        for lab, a in (("ṁ_in", mi), ("ṁ_out", mo)):
            if lab == "ṁ_out" and not out_ok:
                continue
            rel = _span(a) / float(a.mean())
            Q.check(f"(max − min)/mean {lab}", rel, th["mdot_change_rel"], rel <= th["mdot_change_rel"])
    else:
        Q.reason(CAT_DATA, "ṁ_in・ṁ_out の系列が無い (series_csv の r だけ) — mdot_change_rel を判定できない")
    if stop is not None:
        k = idx[stop]
        if mi_all is not None and not float(mi_all[k]) > 0:
            Q.r["values"]["stop"] = {"step": stop, "mdot_in": float(mi_all[k]), "mdot_out": float(mo_all[k])}
            Q.reason(CAT_DATA, f"停止の観測 (step {stop}) の入口の流量 ṁ_in {float(mi_all[k]):.6g} が正でない (r の分母が成り立たない)")
            return Q.finish()
        rs = r_all[k]
        Q.r["values"]["stop"] = {"step": stop, "r": float(rs)}
        Q.check(f"停止の観測 (step {stop}) の r の格子の終端からの差", abs(float(rs - r[-1])), th["swing"], abs(float(rs - r[-1])) <= th["swing"])
        r2 = np.append(r, rs)
        Q.check("停止の観測を加えた |mean r|", abs(float(r2.mean())), th["mean_abs"], abs(float(r2.mean())) <= th["mean_abs"])
        Q.check("停止の観測を加えた max r − min r", _span(r2), th["swing"], _span(r2) <= th["swing"])
        if mi_all is not None:
            # (正常な経路の values・checks は変えない: 量の digest を停止規則の再評価と照合するため)
            stop_out_ok = _mdot_out_positive(Q, np.asarray([mo_all[k]]), [stop], "停止の観測の")
            for lab, a, aa in (("ṁ_in", mi, mi_all), ("ṁ_out", mo, mo_all)):
                if lab == "ṁ_out" and not (out_ok and stop_out_ok):
                    continue
                a2 = np.append(a, aa[k])
                rel = _span(a2) / float(a2.mean())
                Q.check(f"停止の観測を加えた (max − min)/mean {lab}", rel, th["mdot_change_rel"], rel <= th["mdot_change_rel"])
    return Q.finish()


def _mdot_out_positive(Q: _Q, mo, steps, where: str) -> bool:
    """出口の流量 ṁ_out が全点で正か (流出がある・符号が反転しない)。正でない点があれば理由つきの不合格 (閾値の超過の分類) を足して
    False。判定不能にしない理由: 呼ぶのは ṁ_in > 0 のときだけで、r = 1 − ṁ_out/ṁ_in ≥ 1 は測れており収支が取れていないことは確か
    (データの欠けではない)。ṁ_out の相対の変化は定義できないので、呼び出し側は計算しない。"""
    mo = np.asarray(mo, dtype=float)
    bad = np.flatnonzero(~(mo > 0))
    if bad.size == 0:
        return True
    bad_s = [steps[i] for i in bad]
    mn = float(mo.min())
    Q.r["checks"].append({"check": f"{where}出口の流量 ṁ_out が正 (流出がある)", "value": mn, "limit": "> 0", "ok": False,
                          "detail": f"正でない点 {len(bad_s)}/{len(mo)} (step {bad_s[:6]}{'…' if len(bad_s) > 6 else ''})"})
    Q.reason(CAT_THRESHOLD, f"{where}出口の流量 ṁ_out が正でない点がある (最小 {mn:.6g}、step {bad_s[:6]}{'…' if len(bad_s) > 6 else ''}; "
                            "流出が無い・逆流・符号の反転 — r = 1 − ṁ_out/ṁ_in ≥ 1 で流量の収支が取れていない。ṁ_out の相対の変化は定義できない)")
    return False


# ============================================================ 節点の量
def _exit_column(obs: dict, ni: int, nj: int, eta_range):
    ix = obs["cols"].get(ni - 1)
    if ix is None or len(ix) != nj or not np.array_equal(obs["j"][ix], np.arange(nj)):
        raise _DataMissing(f"{Path(obs['path']).name}: 出口の列 i = {ni - 1} の全節点 (nj = {nj}) が無い")
    R = obs["R"][ix]
    eta = R / R[-1]
    core = (eta >= eta_range[0]) & (eta <= eta_range[1])
    if not core.any():
        raise _DataMissing(f"η ∈ {eta_range} の節点が出口の列に無い")
    return ix, core


def exit_core_mach_of(obs: dict, ni: int, nj: int, eta_range) -> float:
    """出口コア M (stop_rule.recorded_quantities の exit_core_M と同じ: M[-1][core].mean())。"""
    ix, core = _exit_column(obs, ni, nj, eta_range)
    return float(obs["M"][ix][core].mean())


def exit_flow_angle_of(obs: dict, ni: int, nj: int, eta_range) -> np.ndarray:
    ix, core = _exit_column(obs, ni, nj, eta_range)
    return np.degrees(np.arctan2(obs["uy"][ix][core], obs["ux"][ix][core]))


def eta_line_nodes(obs: dict, key: str, eta: float, xq, ni: int, nj: int, diag: dict | None = None) -> np.ndarray:
    """`nozzle_report.eta_line` と同じ計算を節点の部分集合で行う。全節点と同じ値になる条件を確かめ、満たさなければ _DataMissing。

    eta_line: xa = X[:, 0]; i = clip(searchsorted(xa, xq) − 1, 0, len − 2); w = (xq − xa[i])/(xa[i+1] − xa[i]);
    列 ii の値 = np.interp(eta, R[ii]/R[ii, −1], Q[ii]); out = (1 − w) v(i) + w v(i+1)。
    diag (dict) を渡すと、列ごとの補間の対 (j の組) の数と、軸の節点 (j = 0) の重みの最大 `axis_weight_max` を入れる。"""
    xq = np.asarray(xq)
    cols = obs["cols"]
    jj = obs["j"]
    avail = np.asarray(sorted(i for i, ix in cols.items() if jj[ix[0]] == 0 and jj[ix[-1]] == nj - 1), dtype=np.int64)
    if len(avail) < 2:
        raise _DataMissing("軸と壁の節点を持つ列が 2 つ以上無い")
    X = obs["X"]
    xa = np.asarray([X[cols[int(i)][0]] for i in avail], dtype=X.dtype)
    if not np.all(np.diff(xa) > 0):
        raise _DataMissing("列の軸の x が単調に増えない")
    s0 = np.searchsorted(xa, xq)
    i = np.clip(s0 - 1, 0, len(xa) - 2)
    c0, c1 = avail[i], avail[i + 1]
    if np.any(c1 - c0 != 1):
        k = int(np.argmax(c1 - c0 != 1))
        raise _DataMissing(f"xq {float(xq[k]):.6g} を挟む列 {int(c0[k])}・{int(c1[k])} が全格子で隣り合わない (節点の表に列が足りない)")
    if np.any(s0 == 0) and avail[0] != 0:
        raise _DataMissing("xq が節点の表の最初の列より上流にある (全格子の最初の列が無い)")
    if np.any(s0 == len(xa)) and avail[-1] != ni - 1:
        raise _DataMissing("xq が節点の表の最後の列より下流にある (全格子の最後の列が無い)")
    w = (xq - xa[i]) / (xa[i + 1] - xa[i])
    Qv = obs[key]
    R = obs["R"]
    vcol = {}
    for col in np.unique(np.concatenate([c0, c1])):
        ix = cols[int(col)]
        Rc = R[ix]
        xp = Rc / Rc[-1]
        js = jj[ix]
        k = int(np.searchsorted(xp, eta, side="right")) - 1
        if k < 0:
            if js[0] != 0:
                raise _DataMissing(f"列 {int(col)}: 軸の節点が無い")
        elif k < len(xp) - 1 and js[k + 1] != js[k] + 1:
            raise _DataMissing(f"列 {int(col)}: η {eta} を挟む節点 j = {int(js[k])}・{int(js[k + 1])} が隣り合わない (節点の表に足りない)")
        if diag is not None and 0 <= k < len(xp) - 1:
            pair = (int(js[k]), int(js[k + 1]))
            diag.setdefault("pairs", {})
            diag["pairs"][pair] = diag["pairs"].get(pair, 0) + 1
            if js[k] == 0:   # 軸の節点と第 1 層の間で補間する列: 軸の節点の重み 1 − (η − η_0)/(η_1 − η_0)
                wax = float(1.0 - (eta - xp[k]) / (xp[k + 1] - xp[k]))
                diag["axis_weight_max"] = max(diag.get("axis_weight_max", 0.0), wax)
        vcol[int(col)] = np.interp(eta, xp, Qv[ix])
    v0 = np.asarray([vcol[int(c)] for c in c0])
    v1 = np.asarray([vcol[int(c)] for c in c1])
    return (1 - w) * v0 + w * v1


def line_quantities(M_line, Md: float) -> dict:
    """v5c_series.py と同じ式: maxabs = max|100 (M/Md − 1)|、halfamp = 100 (max M − min M)/(2 Md)、d = 100 (M/Md − 1)。"""
    M = np.asarray(M_line)
    d = 100.0 * (M / Md - 1.0)
    return {"maxabs": float(np.max(np.abs(d))), "halfamp": float(100.0 * (M.max() - M.min()) / (2.0 * Md)), "d": d}


def _node_obs(ctx, q, s, Q: _Q):
    obs = ctx["data"].nodes(q["nodes_file"], s)
    if obs["nonfinite"]:
        raise _NonFinite(f"{Path(obs['path']).name}: 非有限の節点の値 {obs['nonfinite']}")
    return obs


def _eval_node_scalar_series(q, ctx, Q, cols_needed):
    """series_csv の系列から格子点の値を返す。戻り ({col: array over C}, {col: stop の値} or None)。"""
    S = ctx["data"].series(q["series_csv"])
    steps = S["steps"]
    if len(steps) > 1 and np.any(np.diff(steps) <= 0):
        raise _DataMissing(f"{S['path']}: step が単調でない・重複")
    idx = {int(s): i for i, s in enumerate(steps)}
    pts = ctx["C"] + ([ctx["stop"]] if ctx["stop"] is not None else [])
    miss = [s for s in pts if s not in idx]
    if miss:
        raise _DataMissing(f"{S['path']}: 格子点 {miss[:6]} の行が無い ({len(miss)}/{len(pts)})")
    out = {}
    stop = {}
    for c in cols_needed:
        a = S["cols"][q["series_csv"][c]]
        seg = (steps >= ctx["start"]) & (steps <= pts[-1])
        if np.any(~np.isfinite(a[seg])):
            raise _NonFinite(f"{S['path']}: 列 {q['series_csv'][c]} に非有限")
        out[c] = np.asarray([a[idx[s]] for s in ctx["C"]])
        if ctx["stop"] is not None:
            stop[c] = float(a[idx[ctx["stop"]]])
    return out, (stop if ctx["stop"] is not None else None)


def _eval_exit_core(q: dict, ctx: dict) -> dict:
    series = "series_csv" in q
    Q = _Q(q, "series_csv" if series else "native")
    th = q["thresholds"]
    C, stop = ctx["C"], ctx["stop"]
    if not _check_points(Q, C, ctx["start"], "窓 W"):
        return Q.finish()
    try:
        if series:
            v, sv = _eval_node_scalar_series(q, ctx, Q, ["value_col"])
            vals = v["value_col"]
            vstop = sv["value_col"] if sv else None
        else:
            st = ctx["crit"]["structured"]
            vals = np.asarray([exit_core_mach_of(_node_obs(ctx, q, s, Q), st["ni"], st["nj"], q["eta_range"]) for s in C])
            vstop = exit_core_mach_of(_node_obs(ctx, q, stop, Q), st["ni"], st["nj"], q["eta_range"]) if stop is not None else None
    except _NonFinite as e:
        Q.reason(CAT_NONFINITE, str(e))
        return Q.finish()
    except _DataMissing as e:
        Q.reason(CAT_DATA, str(e))
        return Q.finish()
    Q.r["values"] = {"steps": C, "exit_core_M": [float(x) for x in vals], "eta_range": q["eta_range"]}
    Q.check("出口コア M の幅 max − min", _span(vals), th["tol"], _span(vals) <= th["tol"])
    if vstop is not None:
        Q.r["values"]["stop"] = {"step": stop, "exit_core_M": vstop}
        jump = abs(vstop - float(vals[-1]))
        Q.check(f"停止の観測 (step {stop}) の格子の終端からの差", jump, th["tol"], jump <= th["tol"])
        v2 = np.append(vals, vstop)
        Q.check("停止の観測を加えた幅", _span(v2), th["tol"], _span(v2) <= th["tol"])
    return Q.finish()


def _eval_line(q: dict, ctx: dict) -> dict:
    series = "series_csv" in q
    Q = _Q(q, "series_csv" if series else "native")
    th = q["thresholds"]
    C, stop = ctx["C"], ctx["stop"]
    if not _check_points(Q, C, ctx["start"], "窓 W"):
        return Q.finish()
    D = None
    dstop = None
    try:
        if series:
            v, sv = _eval_node_scalar_series(q, ctx, Q, ["maxabs_col", "halfamp_col"])
            maxabs, halfamp = v["maxabs_col"], v["halfamp_col"]
            stop_vals = (sv["maxabs_col"], sv["halfamp_col"]) if sv else None
        else:
            st = ctx["crit"]["structured"]
            S = st["scale_m"]
            xq = np.linspace(q["x_range_m"][0] / S, q["x_range_m"][1] / S, q["n"])
            res = []
            diags = []
            for s in C + ([stop] if stop is not None else []):
                obs = _node_obs(ctx, q, s, Q)
                dg = {}
                res.append(line_quantities(eta_line_nodes(obs, "M", q["eta"], xq, st["ni"], st["nj"], diag=dg), q["Md"]))
                diags.append(dg)
            rc = res[:len(C)]
            maxabs = np.asarray([r["maxabs"] for r in rc])
            halfamp = np.asarray([r["halfamp"] for r in rc])
            D = np.vstack([r["d"] for r in rc])
            if stop is not None:
                stop_vals = (res[-1]["maxabs"], res[-1]["halfamp"])
                dstop = res[-1]["d"]
            else:
                stop_vals = None
            Q.r["values"]["xq_r_t"] = [float(xq[0]), float(xq[-1]), int(len(xq))]
            wax = max(d.get("axis_weight_max", 0.0) for d in diags)
            Q.r["values"]["interp_pairs"] = {f"{a}-{b}": n for (a, b), n in sorted(diags[-1].get("pairs", {}).items())}
            Q.r["values"]["axis_weight_max"] = wax
            if wax > AXIS_WEIGHT_MAX:
                Q.reason(CAT_GRID, f"線 η {q['eta']} の補間に軸の節点が重み {wax:.3g} で入る (> {AXIS_WEIGHT_MAX:g}; 第 1 層の η が線より上)"
                         " — 近軸の線の定義が格子に合わない。η を第 1 層以上にするか格子を見直す")
                return Q.finish()
    except _NonFinite as e:
        Q.reason(CAT_NONFINITE, str(e))
        return Q.finish()
    except _DataMissing as e:
        Q.reason(CAT_DATA, str(e))
        return Q.finish()
    Q.r["values"].update(steps=C, maxabs=[float(x) for x in maxabs], halfamp=[float(x) for x in halfamp], eta=q["eta"], Md=q["Md"])
    Q.check("最大偏差 maxabs の幅 [%pt]", _span(maxabs), th["tol_maxabs"], _span(maxabs) <= th["tol_maxabs"])
    Q.check("片振幅 halfamp の幅 [%pt]", _span(halfamp), th["tol_halfamp"], _span(halfamp) <= th["tol_halfamp"])
    if D is not None:
        pos = D.max(axis=0) - D.min(axis=0)
        kmax = int(np.argmax(pos))
        Q.r["values"]["pos_span_max"] = float(pos[kmax])
        Q.check("位置ごとの d の幅の最大 [%pt]", float(pos[kmax]), th["tol_pos"], float(pos[kmax]) <= th["tol_pos"],
                detail=f"xq の位置 {kmax}")
    else:
        Q.reason(CAT_DATA, "位置ごとの d が無い (series_csv の maxabs・halfamp だけ) — tol_pos を判定できない")
    if stop_vals is not None:
        Q.r["values"]["stop"] = {"step": stop, "maxabs": stop_vals[0], "halfamp": stop_vals[1]}
        for lab, a, sv_, tol in (("maxabs", maxabs, stop_vals[0], th["tol_maxabs"]), ("halfamp", halfamp, stop_vals[1], th["tol_halfamp"])):
            jump = abs(sv_ - float(a[-1]))
            Q.check(f"停止の観測 (step {stop}) の {lab} の格子の終端からの差", jump, tol, jump <= tol)
            Q.check(f"停止の観測を加えた {lab} の幅", _span(np.append(a, sv_)), tol, _span(np.append(a, sv_)) <= tol)
        if dstop is not None:
            D2 = np.vstack([D, dstop])
            p2 = float((D2.max(axis=0) - D2.min(axis=0)).max())
            Q.check("停止の観測を加えた位置ごとの d の幅の最大", p2, th["tol_pos"], p2 <= th["tol_pos"])
    return Q.finish()


def _eval_angle(q: dict, ctx: dict) -> dict:
    Q = _Q(q)
    th = q["thresholds"]
    C, stop = ctx["C"], ctx["stop"]
    if not _check_points(Q, C, ctx["start"], "窓 W"):
        return Q.finish()
    st = ctx["crit"]["structured"]
    try:
        A = np.vstack([exit_flow_angle_of(_node_obs(ctx, q, s, Q), st["ni"], st["nj"], q["eta_range"]) for s in C])
        Astop = exit_flow_angle_of(_node_obs(ctx, q, stop, Q), st["ni"], st["nj"], q["eta_range"]) if stop is not None else None
    except _NonFinite as e:
        Q.reason(CAT_NONFINITE, str(e))
        return Q.finish()
    except _DataMissing as e:
        Q.reason(CAT_DATA, str(e))
        return Q.finish()
    span = A.max(axis=0) - A.min(axis=0)
    k = int(np.argmax(span))
    Q.r["values"] = {"steps": C, "n_nodes": int(A.shape[1]), "span_max_deg": float(span[k]), "theta_at_T_deg": [float(x) for x in A[-1]]}
    Q.check("節点ごとの流れ角の幅の最大 [deg]", float(span[k]), th["tol"], float(span[k]) <= th["tol"], detail=f"コアの {k} 番目の節点")
    if Astop is not None:
        jump = float(np.abs(Astop - A[-1]).max())
        Q.check(f"停止の観測 (step {stop}) の格子の終端からの差の最大 [deg]", jump, th["tol"], jump <= th["tol"])
        A2 = np.vstack([A, Astop])
        s2 = float((A2.max(axis=0) - A2.min(axis=0)).max())
        Q.check("停止の観測を加えた幅の最大 [deg]", s2, th["tol"], s2 <= th["tol"])
    return Q.finish()


# ============================================================ 壁の熱流
def _eval_wall(q: dict, ctx: dict) -> dict:
    Q = _Q(q)
    th = q["thresholds"]
    C, T, stop = ctx["C"], ctx["T"], ctx["stop"]
    if not _check_points(Q, C, ctx["start"], "窓 W"):
        return Q.finish()
    worstC = None       # 窓の格子点の最悪 (|ΔQ_i| / 許容_i の最大)
    worstS = None       # 停止の観測の最悪
    vals = {"walls": [], "steps": C}
    rng = q.get("measure", WALL_TERMINAL) == WALL_RANGE
    worstR = None       # range: 節点ごとの窓の幅 / 許容 の最大
    worstR2 = None      # range: 停止の観測を加えた幅 / 許容 の最大
    try:
        for wl in q["walls"]:
            ref = ctx["data"].wall(q["file_pattern"], q["dataset"], wl, T)
            if not np.all(np.isfinite(ref["Q"])):
                nbad = int(np.count_nonzero(~np.isfinite(ref["Q"])))
                raise _DataMissing(f"{Path(ref['path']).name}: {q['dataset']} に非有限 {nbad}/{len(ref['Q'])} 点 "
                                   "(適用外の NaN — 断熱・非定常・interfaceDiag の条件外 — か発散かを区別できない)")
            A, how = wall_node_area(ref["coord"], ref["conne"], q["axisymmetric"])
            Qref = ref["Q"]
            tol = th["q_tol"] * A + th["eps_rel"] * np.abs(Qref)
            vals["walls"].append({"name": wl["name"], "physID": wl["physID"], "n_nodes": int(len(Qref)), "area_method": how,
                                  "area_sum": float(A.sum()), "Q_sum_at_T": float(np.sum(Qref))})
            QC = []             # range: 窓の格子点の値
            Qstop = None
            for s in C + ([stop] if stop is not None else []):
                ob = ctx["data"].wall(q["file_pattern"], q["dataset"], wl, s)
                if not np.array_equal(ob["coord"], ref["coord"]):
                    raise _DataMissing(f"{Path(ob['path']).name}: 壁の節点の COORD が終端 {T} と違う (節点の順を対応させられない)")
                if not np.all(np.isfinite(ob["Q"])):
                    nbad = int(np.count_nonzero(~np.isfinite(ob["Q"])))
                    raise _DataMissing(f"{Path(ob['path']).name}: {q['dataset']} に非有限 {nbad}/{len(ob['Q'])} 点 "
                                       "(適用外の NaN — 断熱・非定常・interfaceDiag の条件外 — か発散かを区別できない)")
                dQ = np.abs(ob["Q"] - Qref)
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio = np.where(tol > 0, dQ / tol, np.where(dQ > 0, np.inf, 0.0))
                k = int(np.argmax(ratio))
                cand = {"wall": wl["name"], "step": s, "node": k, "ratio": float(ratio[k]), "dQ": float(dQ[k]), "tol": float(tol[k]),
                        "Q_ref": float(Qref[k]), "x": float(ref["coord"][k, 0]), "r": float(ref["coord"][k, 1])}
                if s == stop:
                    Qstop = ob["Q"]
                    if worstS is None or cand["ratio"] > worstS["ratio"]:
                        worstS = cand
                else:
                    if rng:
                        QC.append(ob["Q"])
                    if worstC is None or cand["ratio"] > worstC["ratio"]:
                        worstC = cand
            if rng:
                # 節点ごとの窓の幅 max − min (v2; 終端からの差では窓の中の往復を見落とす — codex の例 [0.25, 1.75, 1, 1, 1, 1] は
                # 終端からの差の最大 0.75 だが幅 1.5)。許容の相対の参照は終端 T の値 |Q_i(T)|
                for lab, M_ in (("C", np.vstack(QC)), ("C+stop", np.vstack(QC + [Qstop]) if Qstop is not None else None)):
                    if M_ is None:
                        continue
                    span = M_.max(axis=0) - M_.min(axis=0)
                    with np.errstate(divide="ignore", invalid="ignore"):
                        ratio = np.where(tol > 0, span / tol, np.where(span > 0, np.inf, 0.0))
                    k = int(np.argmax(ratio))
                    cand = {"wall": wl["name"], "node": k, "ratio": float(ratio[k]), "span": float(span[k]), "tol": float(tol[k]),
                            "Q_ref": float(Qref[k]), "Q_max": float(M_[:, k].max()), "Q_min": float(M_[:, k].min()),
                            "x": float(ref["coord"][k, 0]), "r": float(ref["coord"][k, 1])}
                    if lab == "C":
                        if worstR is None or cand["ratio"] > worstR["ratio"]:
                            worstR = cand
                    elif worstR2 is None or cand["ratio"] > worstR2["ratio"]:
                        worstR2 = cand
    except _DataMissing as e:
        Q.reason(CAT_DATA, str(e))
        return Q.finish()
    if rng:
        vals["measure"] = WALL_RANGE
        vals["worst_window"] = worstR
        vals["worst_terminal_diff"] = worstC          # 記録だけ (v1 の量)
        Q.r["values"] = vals
        w = worstR
        Q.check("窓の節点の熱量の幅 / 許容 (max_s Q_i − min_s Q_i) / (q_tol·A_i + eps_rel·|Q_i(T)|) の最大", w["ratio"], 1.0, w["ratio"] <= 1.0,
                detail=f"壁 {w['wall']} 節点 {w['node']} (x {w['x']:.6g}, r {w['r']:.6g}): 幅 {w['span']:.4g} "
                       f"({w['Q_min']:.4g}〜{w['Q_max']:.4g}), 許容 {w['tol']:.4g} (相対の参照は終端 {T} の値)")
        if worstS is not None:
            vals["worst_stop"] = worstS
            vals["worst_window_with_stop"] = worstR2
            Q.check(f"停止の観測 (step {stop}) の節点の熱量の格子の終端からの差 / 許容 の最大", worstS["ratio"], 1.0, worstS["ratio"] <= 1.0,
                    detail=f"壁 {worstS['wall']} 節点 {worstS['node']}")
            Q.check("停止の観測を加えた節点の熱量の幅 / 許容 の最大", worstR2["ratio"], 1.0, worstR2["ratio"] <= 1.0,
                    detail=f"壁 {worstR2['wall']} 節点 {worstR2['node']}")
        return Q.finish()
    vals["worst_window"] = worstC
    Q.r["values"] = vals
    w = worstC
    Q.check("窓の節点の熱量の差 / 許容 |ΔQ_i| / (q_tol·A_i + eps_rel·|Q_ref,i|) の最大", w["ratio"], 1.0, w["ratio"] <= 1.0,
            detail=f"壁 {w['wall']} 節点 {w['node']} (x {w['x']:.6g}, r {w['r']:.6g}) step {w['step']}: ΔQ {w['dQ']:.4g}, 許容 {w['tol']:.4g}")
    if worstS is not None:
        vals["worst_stop"] = worstS
        Q.check(f"停止の観測 (step {stop}) の節点の熱量の差 / 許容 の最大", worstS["ratio"], 1.0, worstS["ratio"] <= 1.0,
                detail=f"壁 {worstS['wall']} 節点 {worstS['node']}")
    return Q.finish()


EVALUATORS = {
    "residual": _eval_residual,
    "mass_balance": _eval_mass,
    "exit_core_mach": _eval_exit_core,
    "line_mach": _eval_line,
    "exit_flow_angle": _eval_angle,
    "wall_heat_nodes": _eval_wall,
    KIND_RESIDUAL_FINITE: _eval_residual_finite,
}


# ============================================================ 判定
def read_stopped_at(run_dir) -> dict | None:
    p = Path(run_dir) / "STOPPED_AT"
    if not p.exists():
        return None
    try:
        d = json.loads(p.read_text())
        if not isinstance(d, dict) or d.get("status") != "ok" or not isinstance(d.get("step"), int):
            return {"error": f"{p}: status ok・整数の step が無い"}
        return d
    except Exception as e:  # noqa: BLE001
        return {"error": f"{p}: 読めない ({e})"}


def _launch_info(run: Path) -> dict:
    out = {}
    p = run / "forge_launches.jsonl"
    if p.is_file():
        lines = [ln for ln in p.read_text().splitlines() if ln.strip()]
        if lines:
            try:
                out["forge_launches_last"] = json.loads(lines[-1])
            except json.JSONDecodeError:
                out["forge_launches_last"] = lines[-1]
            out["forge_launches_n"] = len(lines)
    p = run / "RUN_PROVENANCE.txt"
    if p.is_file():
        out["RUN_PROVENANCE_sha256"] = _sha256_bytes(p.read_bytes())
    return out


def _overall(results) -> str:
    v = "SETTLED"
    for r in results:
        if _PRIORITY[r["verdict"]] > _PRIORITY[v]:
            v = r["verdict"]
    return v


def judge(run_dir, end_step: int | None = None, *, data: RunData | None = None, loaded: dict | None = None) -> dict:
    """run を判定する。戻りは結果の dict (`write_verdict` で SETTLE_VERDICT.json に書ける)。
    基準の誤用は CriteriaError (呼び出し側で INDETERMINATE として扱う)。"""
    run = Path(run_dir)
    lc = loaded or load_criteria(run)
    base = {"schema": VERDICT_SCHEMA, "run_dir": str(run.resolve()), "created": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
            "tool": "solver_density_cuda/tools/settle.py", "tool_sha256": _sha256_bytes(Path(__file__).read_bytes())}
    if lc["status"] != "ok":
        base.update(verdict="INDETERMINATE", reason_category=CAT_UNSET, reason=lc["message"],
                    criteria={"path": lc["path"], "sha256": lc.get("sha256")}, quantities=[], provisional=False)
        return base
    crit = lc["criteria"]
    data = data or RunData(run, crit)
    dt = crit["grid"]["out_step_interval"]
    start = crit["grid"]["start_step"]
    W = crit["window"]["W"]
    stopped = read_stopped_at(run)
    if end_step is not None:
        end, src = int(end_step), "--end-step"
    elif stopped is not None and "error" not in stopped:
        end, src = int(stopped["step"]), "STOPPED_AT"
    else:
        end, cands = data.last_step()
        src = "データの最後 (" + "; ".join(f"{c[1]} {c[0]}" for c in cands) + ")" if cands else "データが無い"
    base.update(criteria={"path": lc["path"], "sha256": lc["sha256"], "content": crit},
                launch=_launch_info(run), stopped_at=stopped, grid=crit["grid"], window=crit["window"],
                provisional=any_provisional(crit), end_step=end, end_source=src)
    if stopped is not None and "error" in stopped and end_step is None:
        base.update(verdict="INDETERMINATE", reason_category=CAT_DATA, reason=stopped["error"], quantities=[])
        return base
    if end is None:
        base.update(verdict="INDETERMINATE", reason_category=CAT_DATA, reason="入力のデータが 1 つも無い (残差・流量・節点・壁)",
                    quantities=[_missing_all(q, "入力のデータが無い") for q in crit["quantities"]])
        return base
    T = (end // dt) * dt
    stop = end if end != T else None
    C = list(range(T - W, T + 1, dt))
    P = list(range(T - 2 * W, T - W + 1, dt))
    ctx = {"crit": crit, "data": data, "T": T, "C": C, "P": P, "W": W, "dt": dt, "start": start, "stop": stop, "end": end}
    results = [EVALUATORS[q["kind"]](q, ctx) for q in crit["quantities"]]
    v = _overall(results)
    base.update(terminal_step=T, stop_step=stop, window_steps={"C": C, "P": P}, quantities=results, verdict=v,
                excluded=crit.get("excluded", []), preconditions=crit.get("preconditions", []))
    return base


def _missing_all(q, why):
    Q = _Q(q, "series_csv" if "series_csv" in q else "native")
    Q.reason(CAT_DATA, why)
    return Q.finish()


def scan(run_dir, every: int, *, end_step: int | None = None) -> dict:
    """格子点ごとに「その時点まで」で判定して最初に SETTLED になる step を返す。every は格子の間隔の倍数。"""
    lc = load_criteria(run_dir)
    if lc["status"] != "ok":
        return {"first_settled": None, "reason": lc["message"], "rows": []}
    crit = lc["criteria"]
    dt = crit["grid"]["out_step_interval"]
    if every <= 0 or every % dt != 0:
        raise CriteriaError(f"--scan {every} は格子の間隔 {dt} の正の倍数でない")
    data = RunData(run_dir, crit)
    last, _ = data.last_step()
    if end_step is not None:
        last = end_step
    rows = []
    first = None
    if last is None:
        return {"first_settled": None, "rows": rows}
    start = crit["grid"]["start_step"]
    g0 = ((start + every - 1) // every) * every
    for g in range(max(g0, every), last + 1, every):
        r = judge(run_dir, g, data=data, loaded=lc)
        rows.append({"step": g, "verdict": r["verdict"],
                     "by_quantity": {q["name"]: q["verdict"] for q in r["quantities"]}})
        if r["verdict"] == "SETTLED" and first is None:
            first = g
    return {"first_settled": first, "every": every, "rows": rows}


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return [_jsonable(v) for v in o.tolist()]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        f = float(o)
        return f if math.isfinite(f) else str(f)
    if isinstance(o, float) and not math.isfinite(o):
        return str(o)
    if isinstance(o, bytes):
        return o.decode("utf-8", "replace")
    return o


def files_record(data: RunData | None, extra=()) -> list:
    out = []
    paths = list(data.used.values()) if data else []
    paths += [Path(p) for p in extra]
    seen = set()
    for p in sorted(paths, key=str):
        if str(p) in seen or not p.exists():
            continue
        seen.add(str(p))
        h = hashlib.sha256()
        with open(p, "rb") as f:
            for b in iter(lambda: f.read(1 << 20), b""):
                h.update(b)
        out.append({"path": str(p), "bytes": p.stat().st_size, "sha256": h.hexdigest()})
    return out


def write_verdict(run_dir, result: dict, data: RunData | None = None, out: str | os.PathLike | None = None) -> Path:
    """結果を SETTLE_VERDICT.json に書く (一時ファイル → rename)。使ったファイルの一覧 (大きさ・sha256) を足す。"""
    run = Path(run_dir)
    p = Path(out) if out else run / VERDICT_FILE
    res = dict(result)
    extra = [run / CRITERIA_FILE] if (run / CRITERIA_FILE).exists() else []
    if (run / "STOPPED_AT").exists():
        extra.append(run / "STOPPED_AT")
    res["files"] = files_record(data, extra)
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text(json.dumps(_jsonable(res), ensure_ascii=False, indent=1) + "\n")
    os.replace(tmp, p)
    return p


# ============================================================ 判定の道具の選択と本段の起動との照合
# plan tooling-settle-convergence-and-monitor §4.9 の 1・2 (codex diagnose 2026-10-16 settle-migration-7 の M1・M5)。
# run_case.sh の後処理 (`check_settled.py --select`)・停止規則の後処理 (`stop_rule.postprocess`)・受け渡しの関門
# (`deltastar_loop.handoff_gate`) が同じこの規則を使う。
MAIN_STAGE_RECORD = "MAIN_STAGE_STARTED"       # 設計チェーンの本段の開始の記録 (runner_axismach.run_staged_ns が本段の forge の直前に書く)
LAUNCHES_FILE = "forge_launches.jsonl"         # forge が起動ごとに 1 行足す (cfg_fnv = solverConfig.yaml の FNV-1a 64)
TOOL_NEW, TOOL_OLD = "check_settled.py", "check_convergence.py"


def last_launch(run_dir) -> dict | None:
    """forge_launches.jsonl の最後の非空行 (= 最後の forge の起動)。ファイルが無い → None、読めない・JSON でない → {"error": ...}。"""
    p = Path(run_dir) / LAUNCHES_FILE
    if not p.is_file():
        return None
    try:
        lines = [ln for ln in p.read_text().splitlines() if ln.strip()]
        if not lines:
            return {"error": f"{p} に行が無い"}
        d = json.loads(lines[-1])
        return d if isinstance(d, dict) else {"error": f"{p} の最後の行が JSON の対応表でない"}
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as e:
        return {"error": f"{p} を読めない ({e})"}


def _dotted(d, key: str):
    for k in key.split("."):
        if not isinstance(d, dict) or k not in d:
            return None
        d = d[k]
    return d


def _previous_stage_evidence(run: Path, seg: dict) -> tuple[bool, str]:
    """本段の開始の記録が無い run が前の段 (S1・S2 など) であることを、記録と独立に確かめる (§5.1 #7c、result 段レビュー M1)。
    確かめられるのは、`stage_manifest.json` (run_staged_ns が段ごとに forge の後で書く; 基準の grid.segment.stage_manifest) が読めて
    段を 1 つ以上持ち、どの段も本段の tag (grid.segment.stage_tag、既定 main) でないときだけ。戻り (確かめられた, 理由)。
    run_staged_ns の最初の段 (S1) の後処理の時点ではまだ段の記録が無いので確かめられない (判定不能になる; 段の終了の関門は場で見る)。"""
    man_name = str(seg.get("stage_manifest") or "stage_manifest.json")
    tag = str(seg.get("stage_tag") or "main")
    p = run / man_name
    if not p.is_file():
        return False, f"{man_name} が無い (前の段であることを確かめられない)"
    try:
        man = json.loads(p.read_text())
        stages = man.get("stages") if isinstance(man, dict) else None
    except (OSError, UnicodeDecodeError, ValueError) as e:
        return False, f"{man_name} を読めない ({e})"
    if not isinstance(stages, list) or not stages:
        return False, f"{man_name} に段が無い (前の段であることを確かめられない)"
    tags = [s.get("tag") if isinstance(s, dict) else None for s in stages]
    if tag in tags:
        return False, (f"{man_name} に本段の段 (tag {tag}) がある (段 {tags}) — 本段の開始の記録が欠けているか、前の呼び出しの残りで、"
                       "前の段と確かめられない")
    return True, f"{man_name} の段 {tags} に本段 (tag {tag}) が無い"


def select_judge(run_dir) -> dict:
    """run の今の起動を判定する道具を選ぶ。戻り {tool, mode, why, problems, binding}。
    mode:
      legacy      基準のファイル (settle_criteria.yaml) が無い → check_convergence.py (今どおり)。
      invalid     基準のファイルが読めない・誤り・評価対象が未設定・grid.segment の stage が main でない → check_settled.py
                  (判定不能を返す; **旧の道具へ落とさない**)。
      manual      基準が本段の開始の記録 (grid.segment.start_record) を指さず、run に本段の記録 (MAIN_STAGE_STARTED) も無い
                  = 手で組んだ run → check_settled.py。
      pre_main    基準は本段の開始の記録を指すが記録が無く、**前の段であることを独立に確かめられた** (`_previous_stage_evidence`:
                  stage_manifest.json に本段の段が無い) → check_convergence.py。
      unconfirmed 記録が無く、前の段とも確かめられない (stage_manifest.json が無い・読めない・本段の段がある) → check_settled.py の
                  判定不能 (problems に記録の欠け; §5.1 #7c、result 段レビュー M1)。
      main        本段の記録がある (基準が記録を指さなくても、run に MAIN_STAGE_STARTED があれば既定のキーで) → check_settled.py。
                  照合 (binding、どれかが食い違う・記録が欠けると problems): 記録の event、main_config の sha256・fnv = 今の
                  solverConfig.yaml (main_config.file)、**起動の記録 (forge_launches.jsonl) があり**最後の起動の cfg_fnv = main_config.fnv、
                  記録の settle_criteria.sha256 = 今の基準のファイルの sha256。**最後の起動が本段の記録と違う・起動の記録が無いときも
                  旧の道具へ切り替えず problems** (以前の not_main は result 段レビュー M1 で廃止: 照合の回避の経路だった)。
    記録のキーは基準の grid.segment の config_fnv_key・criteria_sha256_key (既定 main_config.fnv・settle_criteria.sha256)。
    problems があるときの判定は INDETERMINATE (`check_settled.py` は既定で照合する; 受け渡しの関門は上書きでも通さない)。"""
    run = Path(run_dir)
    out = {"tool": TOOL_NEW, "mode": None, "why": None, "problems": [], "binding": None}
    if not os.path.lexists(run / CRITERIA_FILE):
        return {**out, "tool": TOOL_OLD, "mode": "legacy", "why": f"{CRITERIA_FILE} が無い"}
    try:
        lc = load_criteria(run)
    except CriteriaError as e:
        return {**out, "mode": "invalid", "why": f"{CRITERIA_FILE} の誤り", "problems": [f"{CRITERIA_FILE} の誤り: {e}"]}
    if lc["status"] != "ok":
        return {**out, "mode": "invalid", "why": f"{CRITERIA_FILE} の評価対象が未設定", "problems": [lc["reason"]]}
    seg = lc["criteria"]["grid"].get("segment")
    rec_name = seg.get("start_record") if isinstance(seg, dict) else None
    via_default = False
    if not rec_name:
        if not os.path.lexists(run / MAIN_STAGE_RECORD):
            return {**out, "mode": "manual", "why": f"手で組んだ run ({CRITERIA_FILE} の grid.segment が本段の開始の記録を指さない)"}
        # 基準が記録を指さなくても、本段の記録がある run は本段として照合する (基準の segment を外して照合を避ける経路を塞ぐ)
        seg, rec_name, via_default = {}, MAIN_STAGE_RECORD, True
    if seg.get("stage", "main") != "main":
        return {**out, "mode": "invalid", "why": "grid.segment の stage が main でない",
                "problems": [f"grid.segment.stage {seg.get('stage')!r} (本段 main の基準だけを扱う)"]}
    rp = run / str(rec_name)
    if not os.path.lexists(rp):
        ok, ev = _previous_stage_evidence(run, seg)
        if ok:
            return {**out, "tool": TOOL_OLD, "mode": "pre_main",
                    "why": f"本段の開始の記録 {rec_name} が無く、前の段と確かめた ({ev})"}
        return {**out, "mode": "unconfirmed", "why": f"本段の開始の記録 {rec_name} が無く、前の段とも確かめられない",
                "problems": [f"本段の開始の記録 {rec_name} が無い・前の段と確かめられない: {ev} (記録の欠け; 旧の道具へ切り替えない)"]}
    fnv_key = str(seg.get("config_fnv_key") or "main_config.fnv")
    sha_key = str(seg.get("criteria_sha256_key") or "settle_criteria.sha256")
    probs = []
    why_main = (f"本段の起動 ({rec_name} と照合" + (f"; {CRITERIA_FILE} の grid.segment は記録を指さないが run に {MAIN_STAGE_RECORD} がある"
                                                   if via_default else "") + ")")
    try:
        note = json.loads(rp.read_text())
        if not isinstance(note, dict):
            raise ValueError("JSON の対応表でない")
    except (OSError, UnicodeDecodeError, ValueError) as e:
        return {**out, "mode": "main", "why": why_main, "problems": [f"{rec_name} を読めない ({e})"]}
    if note.get("event") != "main_stage_started":
        probs.append(f"{rec_name} の event が main_stage_started でない ({note.get('event')!r})")
    cfg_name = str(_dotted(note, "main_config.file") or "solverConfig.yaml")
    main_fnv = _dotted(note, fnv_key)
    main_sha = _dotted(note, "main_config.sha256")
    crit_sha = _dotted(note, sha_key)
    try:
        cfg_raw = (run / cfg_name).read_bytes()
    except OSError as e:
        cfg_raw = None
        probs.append(f"本段の設定 {cfg_name} を読めない ({e})")
    import stage_manifest as _SM     # 同じ tools の FNV-1a 64 (forge の appendLaunchRecord と同じ定義)
    cfg_fnv = _SM.fnv1a64(cfg_raw) if cfg_raw is not None else None
    cfg_sha = _sha256_bytes(cfg_raw) if cfg_raw is not None else None
    ll = last_launch(run)
    launch_fnv = ll.get("cfg_fnv") if isinstance(ll, dict) and "error" not in ll else None
    binding = {"record": str(rec_name), "config_file": cfg_name, "main_config_fnv": main_fnv, "main_config_sha256": main_sha,
               "config_fnv": cfg_fnv, "config_sha256": cfg_sha, "launch": ll, "launch_cfg_fnv": launch_fnv,
               "record_criteria_sha256": crit_sha, "criteria_sha256": lc["sha256"], "keys": {"config_fnv": fnv_key, "criteria_sha256": sha_key}}
    # 起動の照合 (result 段レビュー M1): 起動の記録が無い・最後の起動が本段の記録と違うときも旧の道具へ切り替えず problems にする
    # (本段の記録がある run は、前の段であることを独立に確かめられない)
    if ll is None:
        probs.append(f"起動の記録 {LAUNCHES_FILE} が無い (最後の forge の起動を本段と照合できない)")
    elif "error" in ll:
        probs.append(f"起動の記録: {ll['error']}")
    elif launch_fnv is None:
        probs.append(f"{LAUNCHES_FILE} の最後の起動に cfg_fnv が無い")
    elif main_fnv is not None and launch_fnv != main_fnv:
        probs.append(f"最後の forge の起動 (cfg_fnv {launch_fnv}) が {rec_name} の {fnv_key} {main_fnv} と違う "
                     "(本段でない起動か記録の不一致; 前の段と確かめられないので旧の道具へ切り替えない)")
    if main_fnv is None:
        probs.append(f"{rec_name} に {fnv_key} が無い (本段の起動を照合できない)")
    elif cfg_fnv is not None and cfg_fnv != main_fnv:
        probs.append(f"今の {cfg_name} の fnv {cfg_fnv} が {rec_name} の {fnv_key} {main_fnv} と違う")
    if main_sha is None:
        probs.append(f"{rec_name} に main_config.sha256 が無い")
    elif cfg_sha is not None and cfg_sha != main_sha:
        probs.append(f"今の {cfg_name} の sha256 {cfg_sha[:16]}… が {rec_name} の main_config.sha256 {str(main_sha)[:16]}… と違う")
    if crit_sha is None:
        probs.append(f"{rec_name} に {sha_key} が無い (基準を本段に結び付けられない)")
    elif crit_sha != lc["sha256"]:
        probs.append(f"今の {CRITERIA_FILE} の sha256 {lc['sha256'][:16]}… が {rec_name} の {sha_key} {str(crit_sha)[:16]}… と違う")
    return {**out, "mode": "main", "why": why_main, "problems": probs, "binding": binding}


def verdict_line(result: dict) -> str:
    v = result["verdict"]
    s = f"VERDICT: {v}"
    if v == "INDETERMINATE" and result.get("reason_category") == CAT_UNSET:
        s += " (評価対象が未設定)"
    if result.get("provisional"):
        s += f" ({PROVISIONAL_TEXT})"
    return s
