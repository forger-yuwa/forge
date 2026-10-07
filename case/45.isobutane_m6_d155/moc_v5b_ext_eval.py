"""plan discretization-moc-axis-limit-and-corrector §6 V5b (別登録 2026-10-07、commit 4dd94af2) の判定と、延長 run の準備・検査。
根拠: notes/reviews/2026-10-07-moc-v5-hold-next-step-diagnose.md。V5 の結果 (保留、_band_ab/moc_v5_euler_eval.json) は書き換えない。
出力は _band_ab/moc_v5b_ext_eval.json だけ (V5 の評価器 moc_v5_euler_eval.py は import して関数を使うが、evaluate は呼ばない)。

run (親 → 子。親の res_18000.h5 から同一メッシュの restart で追加 36000 step、通算 54000。soft 段なし、再開は全 run 1 回):
  腕 B  run_0143〜0145 → run_0154〜0156 (…_monoG1_r{1,2,3}_ext36k)
  腕 M  run_0150〜0152 → run_0157〜0159 (…_mocG1_r{1,2,3}_ext36k)
  ISEN  run_0153       → run_0160       (…_icdep_mocG1_isen_ext36k)
時系列の接続: 親の wallfit_series_v5.csv (本段 step 1000〜18000) と子の wallfit_series_v5b.csv (local step 1000〜36000) を、親の終了 step
  PARENT_END = 18000 を明示して「通算 = local + 18000」で繋ぐ (子の restart 元が <親>/res_18000.h5 であることは延長の記録で照合する)。
  重複 (子の local 0 = 親の 18000 と同じ場を含む)・欠落・余分・逆順・非有限は接続を拒否する (その run は保留)。
  既存の延長の経路 (eval_wallfit_euler.series_files の、_ext6k を探して 12000 を固定加算する経路) は使わない。
  親と子の時系列は同じ評価器 (2 つの時系列の記録の evaluator_sha256 が同じ) と同じ評価座標 (記録の geometry が完全一致) で作ったこと。
窓: 再開後 6000 step (通算 18000〜24000) は使わない。窓 A = 通算 24000〜36000、窓 B = 通算 42000〜54000、各 13 枚 (1000 ごと)。
量: V5 の 6 量 (Δq = M 波 0.001・P 波 0.010・オーバーシュート 0.003・出口規格化オーバーシュート 0.003・|P 傾き| 0.03 [%pt]、
  |出口コア M − 6| 1.8e-4) と符号付きの出口コア M − 6 (列 exit_core_M。準定常の判定は V5 と同じく exit_core_M そのものに当てる)。
  許容 tol_q = Δq/10 (6 量)、符号付き出口 M は 1e-5 (登録の値)。
記録: 各窓の平均・幅 (max − min)・回帰の傾き (step あたりと窓 12000 step あたり)、窓の末尾 5 枚と全 13 枚の check_quasisteady の VERDICT
  (drift 0.05・osc 0.10。moc_v5_euler_eval.quasisteady_window と同じ呼び方)、通算 1000〜54000 の連続系列の極大・極小 (下の転回点、
  h = tol_q) と最大・最小。
判定 (量ごと・run ごと):
  減衰側: 窓 B が V5 の追加登録の準定常の条件を満たす (窓 B の末尾 5 枚・全 13 枚と、子の記録の準定常ファイル QUASISTEADY_wallfit_v5b.txt
    [子の系列の末尾 5 枚 = 窓 B の末尾 5 枚] が STEADY。例外は exit_M_dev だけで、moc_v5_euler_eval.exit_exception の「絶対許容内」を
    窓 B の符号付き出口 M − 6 に当てる)、かつ窓 A → B の平均の移動 |m_B − m_A| ≤ tol_q → その精度で持続振動説を退ける。
  持続振動側: 下の周期の定義で、通算 24000〜54000 の 31 枚に完全な周期が 3 以上あり解像の条件を満たし、隣り合う周期の平均の差と
    半振幅の差がすべて ≤ tol_q、かつ振幅が減衰しない (最後の周期の半振幅 ≥ 最初の周期の半振幅 − tol_q) → 減衰過渡説を退ける。
    この場合も V5 は保留のまま (平均 ± 振幅による採否は別に登録する)。
  判別不能: どちらにも当たらない。両方に当たる場合も判別不能とする (登録はこの場合を決めていないので、どちらの説も退けない側に倒す)。
  run の分類: 前提が不成立・早期停止・欠損なら「保留」。7 量すべてが減衰側なら「減衰側」。持続振動側の量があれば「持続振動側を含む」。
    それ以外は「判別不能」。
周期の数え方 (実装時の定義。登録の「少なくとも 3 周期を解像」を具体化したもの。結果を見る前に決めた):
  (1) 転回点: 振れのしきい値 h = tol_q のヒステリシスで極大・極小を交互に確定する (zigzag)。極大の候補から h を超えて下がったら
      その候補を極大として確定し、極小の候補から h を超えて上がったら極小として確定する。区間の最初の点は転回点にしない
      (区間の外の値が分からないので、端が極値かどうか決まらない)。最後の未確定の候補も転回点にしない。
  (2) 半周期 = 隣り合う転回点の間、周期 = 同じ種類の転回点の間。転回点 T0, T1, T2, … に対し周期 j = [T_2j, T_2j+2] (重ならない)。
      完全な周期の数 = ⌊(転回点の数 − 1)/2⌋。
  (3) 解像の条件: 完全な周期が 3 以上、かつすべての半周期が出力間隔の 2 倍 (2000 step) 以上 (1 間隔ごとの上下は、出力間隔より
      短い振動の折り返しと区別できない)。振幅の分解能は h (転回点の間の振れが h を超えるものだけを数える)。
  (4) 周期 j の平均 = T_2j 以上 T_2j+2 未満の枚の平均。半振幅 = T_2j〜T_2j+2 (両端を含む) の (最大 − 最小)/2。
窓 B での比較 (V5b の採否): 7 本すべてが減衰側で前提が成立したときだけ、窓 B の 13 枚で V5 と同じ比較を行い V5b の結果とする
  (moc_v5_euler_eval の arm_stats・judge_one_sided・judge_two_sided・judge_exit_calibration・_ic_overall: 6 量・Δq・平均と SE の
  定義・D ± 2SE・IC 依存の両側判定と出口コア M そのものの差 ≤ 1e-4・出口較正の 3 区分)。1 本でも減衰側でなければ、窓 B の差は
  参考値として記録するだけ (判定の欄は「参考値」、前提成立時の判定は別欄)。
前提 (欠損・読めない証拠は既定値で通さず不成立):
  V5 の静的な前提: 親の run について moc_v5_euler_eval.preconditions をそのまま呼ぶ (時系列の記録 wallfit_series_v5.json と親の CSV、
    壁の証拠、MOC のゲート、IC の写像と IC の検査 C1〜C5 の 5 キー [V5 の記録 IC_INSPECTION.json を読む]、ISEN の壁と格子、親どうしの
    実効設定)。ただし V5 の判定区間 (親の本段) の残差判定の理由は V5b の前提に使わない (V5b の判定区間は子の本段。記録には残す)。
    V5 の評価器が登録版 (sha256 V5_EVAL_SHA256) でなければ不成立 (理由の文の読み分けが登録版の文言に依存するため)。
  子 (settings): 延長の記録 V5B_EXTENSION.json (verify-child) が OK・乾式でない・親と restart 元 (<親>/res_18000.h5) が一致・
    restart_field が SRC とビット一致・restart 元の保存量 (ro・roU*・roe・roY* など) がすべて移り値もビット一致・nozzle.h5 の /VALUE
    以外が親とビット一致。solverConfig.yaml が親と YAML として time.last.nStepOuter (18000 → 36000) だけ違い出力間隔 1000、
    bcondConfig.yaml が同じ、forge の sha256 が run_0143 の RUN_PROVENANCE と同じ、起動時の実効値の行と環境変数の行 (FORGE_BIN を除く)
    が親と同じ、stage_manifest.json が main の 1 段で hard キーが親の main と同じ。時系列の記録 (wallfit_series_v5b.json) にあり
    status ok・CSV と準定常ファイルの sha256 が記録と同じ・自身の X_E・X_F と基準の差が 1e-6 r_t 以下で親と同じ。
  子 (data): RUN_RC が 0、EARLY_STOP.txt が無い (起動スクリプトが残差の NaN・Inf で止めた印)、子の判定区間の残差判定
    (CONVERGENCE_VERDICT_segment.txt) が DIVERGED・欠損・判定不能でない、接続した系列が 1000〜54000 でそろい有限。
  比較の前提には、V5 と同じく子の判定区間の残差判定が pass か既知の plateau であることを加える。
SE は自己相関を補正しない実務の指標で、振動の平均の不確かさには流用しない (諮問 Major)。
出力 JSON に評価器自身の sha256、plan の登録 commit (4dd94af2)、V5 の評価器の sha256 を記録する。
usage:
  python3 moc_v5b_ext_eval.py [case_dir]                       判定 (→ _band_ab/moc_v5b_ext_eval.json と標準出力の表)
  python3 moc_v5b_ext_eval.py pairs                            子・親・腕の一覧 (起動スクリプト用。1 行 1 run)
  python3 moc_v5b_ext_eval.py prep-child <親> <子> [--dry]     子の入力 (親の入力の複製と nStepOuter 36000、prepare_info) を作る
  python3 moc_v5b_ext_eval.py verify-child <親> <子> [--dry]   restart_field の後に検査し <子>/V5B_EXTENSION.json を書く
  python3 moc_v5b_ext_eval.py run <子>                         forge を本段だけで回す (runner_axismach.run_staged stages="none")
"""
import csv
import hashlib
import json
import math
import re
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
TOOLS = HERE.parents[1] / "solver_density_cuda/tools"
import moc_v5_euler_eval as EV  # noqa: E402  (V5 の判定関数・定数を共有する。evaluate は呼ばない)

PLAN = "plans/active/discretization-moc-axis-limit-and-corrector.md §6 V5b"
PLAN_REG_COMMIT = "4dd94af2"                      # §6 V5b の別登録 (2026-10-07、延長 run の結果を誰も見る前) の commit
V5_EVAL_SHA256 = "f7efa2ca7b13e0a7e55a8b536beac7295f36837903f6cbe53830c833abbedb4c"   # V5 の評価器の登録版 (plan §9)
TAG = "v5b"
PARENT_END, EXT_STEPS, OUT_INT = 18000, 36000, 1000   # 親の終了 step (本段)・追加 step・出力間隔
TOTAL_END = PARENT_END + EXT_STEPS                # 54000
SKIP_AFTER_RESTART = 6000                         # 再開後 6000 step は比較に使わない
WIN_A, WIN_B = (24000, 36000), (42000, 54000)     # 通算 step
OSC_RANGE = (PARENT_END + SKIP_AFTER_RESTART, TOTAL_END)   # 周期を数える区間 (窓 A の始め〜窓 B の終わり)


def _steps(lo: int, hi: int) -> list:
    return list(range(lo, hi + 1, OUT_INT))


WIN_A_STEPS, WIN_B_STEPS, OSC_STEPS = _steps(*WIN_A), _steps(*WIN_B), _steps(*OSC_RANGE)   # 13・13・31 枚
QTYS = list(EV.QS_COLS)                           # V5 の 6 量 + exit_core_M (符号付き出口 M − 6 として扱う)
SIGNED_QTY = "exit_core_M"
SIGNED_LABEL = "符号付き出口コア M − 6"
TOL = {**{q: EV.DQ[q] / 10.0 for q in EV.DQ}, SIGNED_QTY: 1e-5}   # 窓 A → B の平均の移動・周期の差・転回点のしきい値
MIN_PERIODS = 3
MIN_HALF_INTERVALS = 2                            # 半周期 ≥ 2 出力間隔 (2000 step)
CHILDREN = {                                      # 子 → 親 (親が決めた run 名。plan §6 V5b)
    "run_0154_euler_wallfit_monoG1_r1_ext36k": "run_0143_euler_wallfit_monoG1_r1",
    "run_0155_euler_wallfit_monoG1_r2_ext36k": "run_0144_euler_wallfit_monoG1_r2",
    "run_0156_euler_wallfit_monoG1_r3_ext36k": "run_0145_euler_wallfit_monoG1_r3",
    "run_0157_euler_wallfit_mocG1_r1_ext36k": "run_0150_euler_wallfit_mocG1_r1",
    "run_0158_euler_wallfit_mocG1_r2_ext36k": "run_0151_euler_wallfit_mocG1_r2",
    "run_0159_euler_wallfit_mocG1_r3_ext36k": "run_0152_euler_wallfit_mocG1_r3",
    "run_0160_euler_icdep_mocG1_isen_ext36k": "run_0153_euler_icdep_mocG1_isen",
}
_C = list(CHILDREN)
ARMS = {"B": _C[0:3], "M": _C[3:6]}
ISEN = _C[6]
ROLE = {**{c: "B" for c in ARMS["B"]}, **{c: "M" for c in ARMS["M"]}, ISEN: "ISEN"}
if list(CHILDREN.values()) != EV.ARMS["B"] + EV.ARMS["M"] + [EV.ISEN]:
    raise RuntimeError("CHILDREN の親が V5 の腕 (moc_v5_euler_eval.ARMS・ISEN) と一致しない")
GEOM_REF_PARENT = EV.GEOM_REF                     # 評価座標の基準 = 腕 B の r1 (run_0143)
EXT_RECORD = "V5B_EXTENSION.json"
EARLY_STOP_FILE = "EARLY_STOP.txt"                # 起動スクリプトが残差の NaN・Inf で止めたときに書く
RUN_RC_FILE = "RUN_RC"
SERIES_CSV, QS_FILE = f"wallfit_series_{TAG}.csv", f"QUASISTEADY_wallfit_{TAG}.txt"
SERIES_REC_V5, SERIES_REC = f"_band_ab/wallfit_series_{EV.TAG}.json", f"_band_ab/wallfit_series_{TAG}.json"
OUT_JSON = "_band_ab/moc_v5b_ext_eval.json"
V5_SEGMENT_PREFIX = "判定区間の残差判定が "         # moc_v5_euler_eval.run_preconditions の親の判定区間の理由の書き出し (登録版の文言)
NSTEP_PATH = ("time", "last", "nStepOuter")
# 子へ複製する親の入力 (出力・判定ファイル・親の IC の記録は持ち込まない)
COPY_REQUIRED = ("solverConfig.yaml", "bcondConfig.yaml", "nozzle.h5", "prepare_info.json", "wall_design.csv")
COPY_OPTIONAL = ("probe.yaml", "species_meta.yaml", "species_db_external.yaml", "species_db.yaml", "MESH_QUALITY.txt",
                 "target_axis_M.csv")
COPY_GLOBS = ("resolved_species_*.yaml",)
OUTPUT_LIKE = re.compile(r"^(res_|residual_history|forge_run|RUN_|CONVERGENCE_VERDICT|stage_manifest|wallfit_series|QUASISTEADY|"
                         r"IC_MAP|IC_INSPECTION|ic_index_map|run_stdout|run_case_stdout|EARLY_STOP|V5B_|metrics|achieved_vs_target)")
CONS_RE = re.compile(r"^(ro|roU[xyz]|roe|roK|roOmega|roGamma|roReth|roXi|roY\d+|rog_.+|roQ[012]_.+)$")   # ic_index_map と同じ
DECAY, OSC, UNDET, HOLD = "減衰側", "持続振動側", "判別不能", "保留"
RUN_OSC = "持続振動側を含む"
REF_LABEL = "参考値 (7 本すべてが減衰側でない; V5b の判定に使わない)"


def _sha(p: Path):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest() if Path(p).is_file() else None


def _json(p: Path) -> dict:
    return json.loads(Path(p).read_text())


def _yaml_strict():
    if str(TOOLS) not in sys.path:
        sys.path.insert(0, str(TOOLS))
    import yaml_strict
    return yaml_strict


def _num(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


# --- 時系列の接続 ---------------------------------------------------------------------------------------------------
class JoinError(ValueError):
    """親と子の時系列を接続できない (重複・欠落・余分・逆順・非有限・列の欠損)。"""


def read_rows(path: Path) -> list:
    with open(path) as f:
        return list(csv.DictReader(f))


def _int_steps(rows: list, who: str) -> list:
    out = []
    for r in rows:
        s = r.get("step")
        try:
            x = float(s)
        except (TypeError, ValueError):
            raise JoinError(f"{who}の時系列の step {s!r} が数でない") from None
        if not math.isfinite(x) or x != int(x):
            raise JoinError(f"{who}の時系列の step {s!r} が有限の整数でない")
        out.append(int(x))
    return out


def join_series(parent_rows: list, child_rows: list, cols, parent_end: int = PARENT_END, ext_steps: int = EXT_STEPS) -> tuple:
    """親 (step 1000〜parent_end) と子 (local step 1000〜ext_steps) を「通算 = local + parent_end」で接続する。
    戻り値 (通算 step の ndarray, {列: ndarray})。重複・逆順・欠落・余分・非有限・列の欠損は JoinError (既定値で埋めない)。"""
    ps, cs = _int_steps(parent_rows, "親"), _int_steps(child_rows, "子")
    for who, s in (("親", ps), ("子", cs)):
        dup = sorted({x for x in s if s.count(x) > 1})
        if dup:
            raise JoinError(f"{who}の時系列に重複した step がある ({dup[:5]})")
        back = [i for i in range(1, len(s)) if s[i] <= s[i - 1]]
        if back:
            raise JoinError(f"{who}の時系列が逆順 (step {s[back[0] - 1]} の後に {s[back[0]]})")
    if 0 in cs:
        raise JoinError(f"子の local step 0 は親の終了 step {parent_end} と同じ場 (重複)")
    for who, s, want in (("親", ps, _steps(OUT_INT, parent_end)), ("子", cs, _steps(OUT_INT, ext_steps))):
        if s != want:
            lack = sorted(set(want) - set(s))
            extra = sorted(set(s) - set(want))
            raise JoinError(f"{who}の時系列が {want[0]}〜{want[-1]} の {OUT_INT} ごとでない"
                            + (f" (欠落 {lack[:5]}{'…' if len(lack) > 5 else ''}" if lack else " (")
                            + (f"; 余分 {extra[:5]}" if extra else "") + ")")
    g = ps + [s + parent_end for s in cs]
    if g != _steps(OUT_INT, parent_end + ext_steps):               # 上の検査で成り立つはず (念のため)
        raise JoinError("接続した系列が 1000 ごとに連続しない")
    vals = {}
    for c in cols:
        v = []
        for st, r in zip(g, parent_rows + child_rows):
            if c not in r:
                raise JoinError(f"列 {c} が無い")
            x = r[c]
            try:
                f = float(x) if x not in ("", None) else math.nan
            except ValueError:
                f = math.nan
            if not math.isfinite(f):
                raise JoinError(f"{c} に非有限値・空欄 (通算 step {st})")
            v.append(f)
        vals[c] = np.array(v)
    return np.array(g), vals


def window(steps, vals: dict, win_steps: list) -> dict:
    """接続した系列から窓の枚を取り出す ({列: ndarray})。窓の step が 1 つでも無ければ JoinError。"""
    idx = {int(s): i for i, s in enumerate(steps)}
    miss = [s for s in win_steps if s not in idx]
    if miss:
        raise JoinError(f"窓の step {miss[:5]} が無い")
    ii = [idx[s] for s in win_steps]
    return {c: np.asarray(v)[ii] for c, v in vals.items()}


def window_stats(win_steps: list, v) -> dict:
    v = np.asarray(v, dtype=float)
    s = np.asarray(win_steps, dtype=float)
    slope = float(np.polyfit(s, v, 1)[0])
    return {"n": int(v.size), "steps": [int(win_steps[0]), int(win_steps[-1])], "mean": float(v.mean()), "width": float(np.ptp(v)),
            "min": float(v.min()), "max": float(v.max()), "slope_per_step": slope, "change_over_window": slope * float(s[-1] - s[0])}


# --- 周期 (持続振動側) と極大・極小 ---------------------------------------------------------------------------------
def turning_points(y, h: float) -> list:
    """しきい値 h のヒステリシスで確定した極大・極小 [(index, "max" | "min")] (交互)。区間の最初の点と最後の未確定の候補は含めない。"""
    y = np.asarray(y, dtype=float)
    if not (_num(h) and h > 0):
        raise ValueError(f"h = {h!r} は正の有限値であること")
    if not np.all(np.isfinite(y)):
        raise ValueError("非有限の値がある")
    n = len(y)
    tps = []
    if n < 2:
        return tps
    imax = imin = 0
    direction, cand = 0, None
    for k in range(1, n):
        if direction == 0:
            if y[k] > y[imax]:
                imax = k
            if y[k] < y[imin]:
                imin = k
            if y[imax] - y[imin] > h:
                if imin < imax:                     # 極小 → 極大へ上がっている
                    tps.append((imin, "min"))
                    direction, cand = 1, imax
                else:
                    tps.append((imax, "max"))
                    direction, cand = -1, imin
        elif direction == 1:
            if y[k] > y[cand]:
                cand = k
            elif y[cand] - y[k] > h:
                tps.append((cand, "max"))
                direction, cand = -1, k
        else:
            if y[k] < y[cand]:
                cand = k
            elif y[k] - y[cand] > h:
                tps.append((cand, "min"))
                direction, cand = 1, k
    return [t for t in tps if t[0] != 0]          # 区間の最初の点は極値かどうか決まらない


def period_analysis(steps, y, h: float) -> dict:
    """持続振動側の判定 (docstring の周期の定義)。戻り値に転回点・周期・解像の条件・隣接差・判定 persistent と不成立の理由。"""
    steps = [int(s) for s in steps]
    y = np.asarray(y, dtype=float)
    tps = turning_points(y, h)
    gaps = [b[0] - a[0] for a, b in zip(tps, tps[1:])]
    n_per = (len(tps) - 1) // 2 if len(tps) >= 3 else 0
    periods = []
    for j in range(n_per):
        a, b = tps[2 * j][0], tps[2 * j + 2][0]
        periods.append({"start_step": steps[a], "end_step": steps[b], "kind": tps[2 * j][1], "n": b - a,
                        "mean": float(y[a:b].mean()), "half_amp": float((y[a:b + 1].max() - y[a:b + 1].min()) / 2.0)})
    dmean = [abs(q["mean"] - p["mean"]) for p, q in zip(periods, periods[1:])]
    damp = [abs(q["half_amp"] - p["half_amp"]) for p, q in zip(periods, periods[1:])]
    why = []
    if n_per < MIN_PERIODS:
        why.append(f"完全な周期 {n_per} < {MIN_PERIODS} (転回点 {len(tps)})")
    short = [g for g in gaps if g < MIN_HALF_INTERVALS]
    if short:
        why.append(f"半周期が {MIN_HALF_INTERVALS} 出力間隔 ({MIN_HALF_INTERVALS * OUT_INT} step) 未満のものが {len(short)} 個 (解像していない)")
    if n_per >= 2:
        if not all(d <= h for d in dmean):
            why.append(f"隣り合う周期の平均の差 最大 {max(dmean):.3g} > {h:.3g}")
        if not all(d <= h for d in damp):
            why.append(f"隣り合う周期の半振幅の差 最大 {max(damp):.3g} > {h:.3g}")
        if not periods[-1]["half_amp"] >= periods[0]["half_amp"] - h:
            why.append(f"振幅が減衰している (最初 {periods[0]['half_amp']:.3g} → 最後 {periods[-1]['half_amp']:.3g}、差 > {h:.3g})")
    return {"range": [steps[0], steps[-1]], "n": len(steps), "h": float(h),
            "turning_points": [{"step": steps[i], "value": float(y[i]), "kind": k} for i, k in tps],
            "half_period_intervals": gaps, "n_periods": n_per, "periods": periods,
            "adjacent_mean_diff": dmean, "adjacent_half_amp_diff": damp,
            "resolved": n_per >= MIN_PERIODS and not short, "persistent": not why, "reasons": why}


def extrema_record(steps, y, h: float) -> dict:
    """連続系列の極大・極小 (しきい値 h の転回点) と最大・最小 (記録用)。"""
    steps = [int(s) for s in steps]
    y = np.asarray(y, dtype=float)
    tps = turning_points(y, h)
    return {"range": [steps[0], steps[-1]], "h": float(h),
            "turning_points": [{"step": steps[i], "value": float(y[i]), "kind": k} for i, k in tps],
            "max": {"step": steps[int(np.argmax(y))], "value": float(y.max())},
            "min": {"step": steps[int(np.argmin(y))], "value": float(y.min())}}


# --- 量ごと・run ごとの判定 -------------------------------------------------------------------------------------------
def classify_run(steps, vals: dict, qfile_text) -> dict:
    """1 run の分類 (接続済みの系列)。qfile_text = 子の QUASISTEADY_wallfit_v5b.txt の本文 (無ければ None → 全量 UNKNOWN)。
    戻り値 {"quantities": {量: {label・窓 A/B の記録・平均の移動・周期・連続系列の極値}}, "label" (run の分類), "by_label",
    "quasisteady_B_reasons", "exit_exception_B", "window_B_values" (窓 B の 13 枚。比較に使う)}。"""
    from throat_mono_judge import ABS_TOL_WITHIN
    A, B = window(steps, vals, WIN_A_STEPS), window(steps, vals, WIN_B_STEPS)
    qsw = {"A": EV.quasisteady_window(A), "B": EV.quasisteady_window(B)}   # 13 枚の並びだけを使う (step の平行移動に依らない)
    qfile = EV.parse_quasisteady(qfile_text or "", QTYS)
    exc = EV.exit_exception(np.asarray(B[SIGNED_QTY], dtype=float) - EV.M_TARGET)
    why_b, status_b = EV.judge_window_quasisteady(qsw["B"], qfile, exc)
    osc_w = window(steps, vals, OSC_STEPS)
    out_q = {}
    for q in QTYS:
        sign = EV.M_TARGET if q == SIGNED_QTY else 0.0           # 符号付き出口 M は M − 6 で記録 (平均の移動・周期の差は平行移動に依らない)
        sa, sb = window_stats(WIN_A_STEPS, A[q] - sign), window_stats(WIN_B_STEPS, B[q] - sign)
        shift = sb["mean"] - sa["mean"]
        tol = TOL[q]
        qs_ok = status_b.get(q) == "STEADY" or (q == EV.EXIT_EXC_QTY and status_b.get(q) == ABS_TOL_WITHIN)
        decay_why = []
        if not qs_ok:
            decay_why.append(f"窓 B の準定常が {status_b.get(q)} (末尾 5 枚 {qsw['B'][q]['tail5']}・全 13 枚 {qsw['B'][q]['full13']}・"
                             f"記録のファイル {qfile.get(q)})")
        if not (math.isfinite(shift) and abs(shift) <= tol):
            decay_why.append(f"窓 A → B の平均の移動 {shift:+.3g} の絶対値 > {tol:.3g}")
        pa = period_analysis(OSC_STEPS, osc_w[q] - sign, tol)
        decaying, persistent = not decay_why, pa["persistent"]
        if decaying and not persistent:
            label = DECAY
        elif persistent and not decaying:
            label = OSC
        else:
            label = UNDET
        out_q[q] = {"label": label, "tol": tol, "name": (SIGNED_LABEL if q == SIGNED_QTY else q),
                    "decaying": decaying, "decaying_reasons": decay_why, "persistent": persistent,
                    "both": decaying and persistent,
                    "window_A": {**sa, "tail5": qsw["A"][q]["tail5"], "full13": qsw["A"][q]["full13"]},
                    "window_B": {**sb, "tail5": qsw["B"][q]["tail5"], "full13": qsw["B"][q]["full13"],
                                 "file_tail5": qfile.get(q), "quasisteady_status": status_b.get(q)},
                    "mean_shift_A_to_B": shift, "oscillation": pa,
                    "extrema_full": extrema_record(steps, np.asarray(vals[q]) - sign, tol)}
    labels = {q: r["label"] for q, r in out_q.items()}
    if all(v == DECAY for v in labels.values()):
        run_label = DECAY
    elif any(v == OSC for v in labels.values()):
        run_label = RUN_OSC
    else:
        run_label = UNDET
    return {"quantities": out_q, "label": run_label, "by_label": labels, "quasisteady_B_reasons": why_b,
            "exit_exception_B": exc, "window_B_values": {q: [float(x) for x in B[q]] for q in QTYS}}


# --- 実効設定と延長の記録 ---------------------------------------------------------------------------------------------
def config_diff(a_text: str, b_text: str) -> list:
    """2 つの solverConfig を YAML として (重複キー・merge key は拒否) 比べ、違う葉 [(パス, a の値, b の値)]。"""
    ys = _yaml_strict()
    A, B = ys.load(a_text), ys.load(b_text)
    out = []
    absent = "<無い>"

    def walk(x, y, path):
        if isinstance(x, dict) and isinstance(y, dict):
            for k in sorted(set(x) | set(y), key=str):
                walk(x.get(k, absent), y.get(k, absent), path + (k,))
        elif x != y or type(x) is not type(y):
            out.append((path, x, y))
    walk(A, B, ())
    return out


def config_failures(parent_text: str, child_text: str) -> tuple:
    """子の solverConfig が親と time.last.nStepOuter (PARENT_END → EXT_STEPS) だけ違い、出力間隔が OUT_INT か。戻り値 (理由, 差分)。"""
    ys = _yaml_strict()
    why = []
    try:
        diff = config_diff(parent_text, child_text)
        P, Cc = ys.load(parent_text), ys.load(child_text)
    except Exception as e:  # noqa: BLE001 — 読めないことは不成立
        return [f"solverConfig を YAML として読めない: {type(e).__name__}: {e}"], None
    want = [(NSTEP_PATH, PARENT_END, EXT_STEPS)]
    if diff != want:
        why.append(f"solverConfig の親との差が time.last.nStepOuter {PARENT_END} → {EXT_STEPS} だけでない "
                   f"({[('.'.join(map(str, p)), a, b) for p, a, b in diff]})")
    oi = (((Cc or {}).get("time") or {}).get("outStepInterval"))
    if not (isinstance(oi, int) and not isinstance(oi, bool) and oi == OUT_INT):
        why.append(f"子の outStepInterval {oi!r} が {OUT_INT} でない")
    return why, [['.'.join(map(str, p)), a, b] for p, a, b in diff]


def provenance(rd: Path) -> dict:
    """RUN_PROVENANCE.txt → forge_sha256・起動時の実効値の行・環境変数の行 (FORGE_BIN を除く)。sha256 が読めなければ ValueError。"""
    lines = (Path(rd) / "RUN_PROVENANCE.txt").read_text().splitlines()
    sha = [l.split(":", 1)[1].strip() for l in lines if l.startswith("forge_sha256")]
    if len(sha) != 1 or not re.fullmatch(r"[0-9a-f]{64}", sha[0]):
        raise ValueError(f"RUN_PROVENANCE.txt の forge_sha256 を読めない ({sha})")
    env = sorted(l.strip() for l in lines if l.startswith("env") and "FORGE_BIN=" not in l)
    return {"forge_sha256": sha[0], "effective": sorted(l.strip() for l in lines if l.startswith("'") and " effective" in l),
            "env": env}


def _manifest_main(rd: Path, want_tags) -> dict:
    st = _json(Path(rd) / "stage_manifest.json")["stages"]
    tags = [s["tag"] for s in st]
    if tags != list(want_tags):
        raise ValueError(f"段の並び {tags} が {list(want_tags)} でない")
    return st[-1]["key"]


def nonvalue_digest(h5path: Path) -> str:
    """h5 の /VALUE 以外の全データセット (名前・型・形・値) の sha256 (メッシュが同じことの照合用。属性は含めない)。"""
    import h5py
    h = hashlib.sha256()
    with h5py.File(h5path, "r") as f:
        names = []
        f.visit(lambda n: names.append(n) if isinstance(f[n], h5py.Dataset) and not (n == "VALUE" or n.startswith("VALUE/")) else None)
        for n in sorted(names):
            d = f[n][()]
            a = np.asarray(d)
            h.update(n.encode() + b"\0" + str(a.dtype).encode() + b"\0" + str(a.shape).encode() + b"\0")
            h.update(a.tobytes() if a.dtype.kind not in "OSU" else repr(a.tolist()).encode())
    return h.hexdigest()


def _parse_list_line(text: str, head: str):
    m = re.search(rf"^{re.escape(head)}\s*:\s*(\[.*\])\s*$", text, re.M)
    if not m:
        return None
    try:
        import ast
        v = ast.literal_eval(m.group(1))
        return [str(x) for x in v] if isinstance(v, list) else None
    except (ValueError, SyntaxError):
        return None


def verify_child(parent: Path, child: Path, dry: bool = False) -> dict:
    """restart_field の後の子の検査 (書き換えはしない。記録 <子>/V5B_EXTENSION.json だけを書く)。戻り値 = 記録 (VERDICT OK / REFUSED)。"""
    import h5py
    parent, child = Path(parent).resolve(), Path(child).resolve()
    src = parent / f"res_{PARENT_END}.h5"
    why = []
    rec = {"plan": PLAN, "tool": "moc_v5b_ext_eval.verify-child", "parent": parent.name, "child": child.name,
           "src_res": f"{parent.name}/{src.name}", "dry": bool(dry), "parent_end_step": PARENT_END, "ext_steps": EXT_STEPS}
    if CHILDREN.get(child.name) != parent.name:
        why.append(f"子 {child.name} の親は {CHILDREN.get(child.name)!r} (登録の対応と違う)")
    # 親: 本段 18000 の最終場が最後の res、本段の段の並び
    res = sorted(parent.glob("res_[0-9]*.h5"), key=lambda f: int(re.findall(r"\d+", f.name)[0]))
    rec["parent_last_res"] = res[-1].name if res else None
    if not src.is_file():
        why.append(f"親の {src.name} が無い")
    elif res[-1].name != src.name:
        why.append(f"親の最後の res が {res[-1].name} ({src.name} であること)")
    try:
        rec["parent_main_key"] = _manifest_main(parent, ("soft", "main"))
    except Exception as e:  # noqa: BLE001
        why.append(f"親の stage_manifest.json: {type(e).__name__}: {e}")
    # 設定 (nStepOuter だけ)
    try:
        cw, cdiff = config_failures((parent / "solverConfig.yaml").read_text(), (child / "solverConfig.yaml").read_text())
        why += cw
        rec["config_diff"] = cdiff
    except OSError as e:
        why.append(f"solverConfig を読めない: {e}")
    try:
        ys = _yaml_strict()
        if ys.load((parent / "bcondConfig.yaml").read_text()) != ys.load((child / "bcondConfig.yaml").read_text()):
            why.append("bcondConfig.yaml が親と違う")
    except Exception as e:  # noqa: BLE001
        why.append(f"bcondConfig を読めない: {type(e).__name__}: {e}")
    # restart_field の記録
    log_p = child / "restart_field.log"
    log = log_p.read_text() if log_p.is_file() else ""
    vline = [l for l in log.splitlines() if l.startswith("VERDICT:")]
    moved, kept, missing = (_parse_list_line(log, "移した保存量"), _parse_list_line(log, "DST を残す"),
                            _parse_list_line(log, "SRC に無く据置"))
    sp_line = [l for l in log.splitlines() if l.startswith("species 属性")]
    rec["restart_field"] = {"log": log_p.name, "verdict_line": vline[-1] if vline else None, "moved": moved, "kept_from_dst": kept,
                            "src_missing_kept": missing, "species_line": sp_line[-1] if sp_line else None,
                            "n_moved": (len(moved) if moved is not None else None)}
    if not vline or not vline[-1].startswith("VERDICT: OK") or "SRC とビット一致" not in vline[-1] or "丸め" in vline[-1]:
        why.append(f"restart_field の VERDICT が「OK・SRC とビット一致」でない ({vline[-1] if vline else '行が無い'})")
    if moved is None:
        why.append("restart_field の「移した保存量」の行を読めない (移した量の数が分からない)")
    if not dry and not (sp_line and "継承" in sp_line[-1]):
        why.append(f"化学種の属性が継承されていない ({sp_line[-1] if sp_line else '行が無い'})")
    # 独立の検査: restart 元の保存量がすべて子の nozzle.h5 にあり、ビット一致 (化学種 roY* の取りこぼしを拾う)
    try:
        with h5py.File(src, "r") as s, h5py.File(child / "nozzle.h5", "r") as d:
            sv, dv = s["VALUE"], d["VALUE"]
            cons = sorted(n for n in sv if CONS_RE.match(n))
            lack = [n for n in cons if n not in dv]
            diff = [n for n in cons if n in dv and not (sv[n].dtype == dv[n].dtype and np.array_equal(sv[n][()], dv[n][()]))]
            common = sorted(n for n in dv if n in sv and n != "wall_dist")
        rec.update(src_conserved=cons, n_src_conserved=len(cons), conserved_missing_in_child=lack, conserved_not_identical=diff,
                   n_common_value_datasets=len(common))
        if not cons:
            why.append("restart 元に保存量が無い")
        if lack:
            why.append(f"restart 元の保存量 {lack} が子の nozzle.h5 に無い (移っていない; 化学種なら skill forge-aws-run §3)")
        if diff:
            why.append(f"保存量 {diff} が restart 元とビット一致しない")
        if moved is not None and sorted(moved) != common:
            why.append(f"restart_field の移した量 {len(moved)} が共通のデータセット {len(common)} と違う")
        if moved is not None and not set(cons) <= set(moved):
            why.append("restart_field の移した量に restart 元の保存量が全部は入っていない")
    except Exception as e:  # noqa: BLE001
        why.append(f"restart 元と子の nozzle.h5 を照合できない: {type(e).__name__}: {e}")
    try:
        rec["nonvalue_digest_child"] = nonvalue_digest(child / "nozzle.h5")
        rec["nonvalue_digest_parent"] = nonvalue_digest(parent / "nozzle.h5")
        if rec["nonvalue_digest_child"] != rec["nonvalue_digest_parent"]:
            why.append("子の nozzle.h5 の /VALUE 以外 (メッシュ) が親とビット一致しない")
    except Exception as e:  # noqa: BLE001
        why.append(f"nozzle.h5 のメッシュを照合できない: {type(e).__name__}: {e}")
    try:
        info = _json(child / "prepare_info.json")
        if bool(info.get("DRY")) != bool(dry):
            why.append(f"prepare_info の DRY {info.get('DRY')!r} が今回の指定 {dry} と違う")
        if info.get("extends") != parent.name:
            why.append(f"prepare_info の extends {info.get('extends')!r} が親 {parent.name} でない")
    except Exception as e:  # noqa: BLE001
        why.append(f"子の prepare_info.json を読めない: {type(e).__name__}: {e}")
    rec.update(VERDICT=("OK" if not why else "REFUSED"), failures=why, nozzle_sha256_after=_sha(child / "nozzle.h5"))
    (child / EXT_RECORD).write_text(json.dumps(rec, indent=1, ensure_ascii=False, default=str))
    return rec


def _referenced_inputs(*texts) -> list:
    """config が参照するファイル名 (yaml・yml・csv・h5・dat・txt) のうち出力らしくないもの。"""
    names = set()
    for t in texts:
        for m in re.finditer(r"[\"']?([A-Za-z0-9_.\-]+\.(?:yaml|yml|csv|h5|dat|txt))[\"']?", t or ""):
            n = m.group(1)
            if "/" not in n and not OUTPUT_LIKE.match(n):
                names.add(n)
    return sorted(names)


def prep_child(parent: Path, child: Path, dry: bool = False) -> dict:
    """子の入力を作る: 親の入力の複製 (出力・判定ファイル・親の IC の記録は持ち込まない)、solverConfig の nStepOuter を 36000 に
    (yaml_strict.replace_scalars で値トークンだけ置換し、親との差がそれだけであることを検査)、prepare_info に延長の記録。
    子が既にあれば止める (既存 run は消さない)。restart_field は起動スクリプトが呼ぶ。"""
    parent, child = Path(parent).resolve(), Path(child).resolve()
    if CHILDREN.get(child.name) != parent.name:
        raise SystemExit(f"子 {child.name} の親は {CHILDREN.get(child.name)!r} (登録の対応と違う) — 止める")
    if child.exists():
        raise SystemExit(f"{child} が既にある — 番号の衝突。止める (既存 run は消さない)")
    for f in COPY_REQUIRED:
        if not (parent / f).is_file():
            raise SystemExit(f"親の {f} が無い — 止める")
    if not (parent / f"res_{PARENT_END}.h5").is_file():
        raise SystemExit(f"親の res_{PARENT_END}.h5 が無い — 止める")
    ptext = (parent / "solverConfig.yaml").read_text()
    ys = _yaml_strict()
    P = ys.load(ptext)
    if ((P.get("time") or {}).get("last") or {}).get("nStepOuter") != PARENT_END:
        raise SystemExit(f"親の nStepOuter が {PARENT_END} でない — 止める")
    if (P.get("time") or {}).get("outStepInterval") != OUT_INT:
        raise SystemExit(f"親の outStepInterval が {OUT_INT} でない — 止める (登録は 1000 ごとの出力)")
    ctext = ys.replace_scalars(ptext, {NSTEP_PATH: str(EXT_STEPS)})
    why, _ = config_failures(ptext, ctext)
    if why:
        raise SystemExit("子の solverConfig の検査が不成立 — 止める: " + "; ".join(why))
    child.mkdir(parents=True, exist_ok=False)
    copied = []
    btext = (parent / "bcondConfig.yaml").read_text()
    extra = [n for n in _referenced_inputs(ptext, btext) if n not in COPY_REQUIRED + COPY_OPTIONAL]
    for f in list(COPY_REQUIRED) + list(COPY_OPTIONAL) + extra:
        if (parent / f).is_file():
            shutil.copy2(parent / f, child / f)
            copied.append(f)
    for g in COPY_GLOBS:
        for p in sorted(parent.glob(g)):
            shutil.copy2(p, child / p.name)
            copied.append(p.name)
    miss = [n for n in _referenced_inputs(ptext, btext) if not (child / n).is_file()]
    (child / "solverConfig.yaml").write_text(ctext)
    info = _json(parent / "prepare_info.json")
    info.update(stages="none", extends=parent.name, restart_from=f"{parent.name}/res_{PARENT_END}.h5 (restart_field)",
                plan=PLAN, nStepOuter=EXT_STEPS, parent_end_step=PARENT_END,
                v5b={"copied": copied, "referenced_missing": miss, "ic_of_parent": info.get("ic")})
    if dry:
        info["DRY"] = True
    else:
        info.pop("DRY", None)
    (child / "prepare_info.json").write_text(json.dumps(info, indent=1, default=str))
    return {"child": child.name, "copied": copied, "referenced_missing": miss}


def run_child(rd: Path) -> int:
    """forge を本段だけで回す (runner_axismach.run_staged stages="none": 子の solverConfig のまま、段の記録は main の 1 段)。"""
    rd = Path(rd).resolve()
    info = _json(rd / "prepare_info.json")
    if info.get("DRY"):
        raise SystemExit(f"{rd} は乾式確認の準備 — 回さない")
    rec = _json(rd / EXT_RECORD)
    if rec.get("VERDICT") != "OK" or rec.get("dry") is not False:
        raise SystemExit(f"{rd}: 延長の記録 {EXT_RECORD} が OK でない ({rec.get('VERDICT')}, dry {rec.get('dry')}) — 回さない")
    if any((rd / f).exists() for f in ("stage_manifest.json", "residual_history.csv")) or list(rd.glob("res_[0-9]*.h5")):
        raise SystemExit(f"{rd}: 既に出力がある — 回さない (既存 run は消さない)")
    if not rec.get("nozzle_sha256_after") or _sha(rd / "nozzle.h5") != rec.get("nozzle_sha256_after"):
        raise SystemExit(f"{rd}: nozzle.h5 が延長の記録 (verify-child) の後に変わった・記録が無い — 回さない")
    sys.path.insert(0, str(HERE.parents[1] / "design"))
    from forge_design.evaluate import runner_axismach as RA
    rc = RA.run_staged(rd, cfl_main=None, mid_stage=False, stages="none")
    last = sorted(rd.glob("res_[0-9]*.h5"), key=lambda f: int(re.findall(r"\d+", f.name)[0]))
    print(f"forge exit={rc} last_res={last[-1].name if last else None}")
    return rc


# --- 子の前提 -----------------------------------------------------------------------------------------------------------
def child_checks(case: Path, child: str, parent: str, srec5: dict, srec5b: dict, ref_prov) -> tuple:
    """1 本の子の前提。戻り値 (settings の理由, data の理由, 記録)。"""
    from throat_mono_judge import SEGMENT_VERDICT_FILE, parse_segment_verdict
    cd, pd = case / child, case / parent
    s_why, d_why, meta = [], [], {"parent": parent, "role": ROLE[child]}
    # 延長の記録
    try:
        ext = _json(cd / EXT_RECORD)
        meta["extension"] = {k: ext.get(k) for k in ("VERDICT", "dry", "src_res", "n_src_conserved", "restart_field", "config_diff")}
        if ext.get("VERDICT") != "OK" or ext.get("dry") is not False or ext.get("parent") != parent \
                or ext.get("src_res") != f"{parent}/res_{PARENT_END}.h5":
            s_why.append(f"延長の記録 {EXT_RECORD} が OK・乾式でない・親 {parent} の res_{PARENT_END}.h5 からでない "
                         f"({ext.get('VERDICT')}, dry {ext.get('dry')!r}, {ext.get('src_res')}): {ext.get('failures')}")
        # 起動直前の照合は run_child が行う (違えば回さない)。ここでは記録だけ (判定の後に nozzle.h5 を消すことがある; skill §5)
        meta["nozzle_sha256_matches_record_now"] = (_sha(cd / "nozzle.h5") == ext.get("nozzle_sha256_after"))
    except Exception as e:  # noqa: BLE001
        s_why.append(f"延長の記録 {EXT_RECORD} を読めない: {type(e).__name__}: {e}")
    try:
        info = _json(cd / "prepare_info.json")
        if info.get("DRY") or info.get("DRY_NO_IC"):
            s_why.append("乾式確認 (DRY) の準備から作った run")
    except Exception as e:  # noqa: BLE001
        s_why.append(f"prepare_info.json を読めない: {type(e).__name__}: {e}")
    # 実効設定 (親と nStepOuter だけ違う・同じバイナリ・同じ実効値と環境変数・main の 1 段で hard キーが同じ)
    try:
        cw, cdiff = config_failures((pd / "solverConfig.yaml").read_text(), (cd / "solverConfig.yaml").read_text())
        s_why += cw
        meta["config_diff"] = cdiff
        ys = _yaml_strict()
        if ys.load((pd / "bcondConfig.yaml").read_text()) != ys.load((cd / "bcondConfig.yaml").read_text()):
            s_why.append("bcondConfig.yaml が親と違う")
    except Exception as e:  # noqa: BLE001
        s_why.append(f"設定を読めない: {type(e).__name__}: {e}")
    try:
        pc, pp = provenance(cd), provenance(pd)
        meta["forge_sha256"] = pc["forge_sha256"]
        if ref_prov is None or pc["forge_sha256"] != ref_prov["forge_sha256"]:
            s_why.append(f"forge の sha256 {pc['forge_sha256'][:12]} が {GEOM_REF_PARENT} の RUN_PROVENANCE "
                         f"({(ref_prov or {}).get('forge_sha256', '読めない')[:12]}) と違う")
        if pc["effective"] != pp["effective"]:
            s_why.append(f"起動時の実効値の行が親と違う ({pc['effective']} / {pp['effective']})")
        if pc["env"] != pp["env"]:
            s_why.append(f"環境変数の行が親と違う ({pc['env']} / {pp['env']})")
    except Exception as e:  # noqa: BLE001
        s_why.append(f"RUN_PROVENANCE を読めない: {type(e).__name__}: {e}")
    try:
        if _manifest_main(cd, ("main",)) != _manifest_main(pd, ("soft", "main")):
            s_why.append("子の main 段の hard キーが親の main と違う")
    except Exception as e:  # noqa: BLE001
        s_why.append(f"段の記録: {type(e).__name__}: {e}")
    # 時系列の記録
    r5, r5b = ((srec5.get("runs") or {}).get(parent)), ((srec5b.get("runs") or {}).get(child))
    if r5b is None:
        s_why.append(f"時系列の記録 ({SERIES_REC}) に無い")
    else:
        if r5b.get("status") != "ok":
            d_why.append(f"子の時系列の評価が {r5b.get('status')!r} ({r5b.get('reason')})")
        for key, name in (("csv_sha256", SERIES_CSV), ("quasisteady_sha256", QS_FILE)):
            got = _sha(cd / name)
            if not r5b.get(key) or got != r5b.get(key):
                s_why.append(f"{name} が時系列の記録と違う (sha256 記録 {str(r5b.get(key))[:12]} / 実際 {str(got)[:12]})")
        for k in ("dX_E", "dX_F"):
            v = r5b.get(k)
            if not (_num(v) and abs(v) <= EV.GEOM_TOL):
                s_why.append(f"自身の {k[1:]} と評価の基準の差 {v!r} が{EV.GEOM_TOL_LABEL} {EV.GEOM_TOL:g} r_t を超える・欠損")
            elif r5 is None or r5.get(k) != v:
                s_why.append(f"自身の {k[1:]} と評価の基準の差 {v!r} が親 ({(r5 or {}).get(k)!r}) と違う")
        meta["series"] = {k: r5b.get(k) for k in ("dX_E", "dX_F", "n_snaps", "last_step", "status")}
    if r5 is None:
        s_why.append(f"親 {parent} が V5 の時系列の記録 ({SERIES_REC_V5}) に無い")
    else:
        got = _sha(pd / f"wallfit_series_{EV.TAG}.csv")
        if not r5.get("csv_sha256") or got != r5.get("csv_sha256"):
            s_why.append(f"親の wallfit_series_{EV.TAG}.csv が V5 の記録と違う (古い成果物か別の呼び出し)")
    # 走り終わり (data)
    rc_p = cd / RUN_RC_FILE
    rc = rc_p.read_text().strip() if rc_p.is_file() else None
    meta["run_rc"] = rc
    if rc != "0":
        d_why.append(f"RUN_RC が 0 でない ({rc!r})")
    if (cd / EARLY_STOP_FILE).exists():
        meta["early_stop"] = (cd / EARLY_STOP_FILE).read_text()[:2000]
        d_why.append(f"早期停止 ({EARLY_STOP_FILE}: 残差の NaN・Inf)")
    p = cd / SEGMENT_VERDICT_FILE
    seg = parse_segment_verdict(p.read_text() if p.is_file() else None)
    meta["segment"] = seg
    if seg["status"] in ("diverged", "missing", "undeterminable"):
        d_why.append(f"子の判定区間の残差判定が {seg['status']} ({seg.get('reason') or seg.get('line')})")
    return s_why, d_why, meta


def v5_static_preconditions(case: Path, srec5: dict, wall_check=None) -> dict:
    """V5 の静的な前提 (親の run)。moc_v5_euler_eval.preconditions をそのまま呼び、親の判定区間の残差判定の理由だけを外す (記録には残す)。"""
    actual = _sha(Path(EV.__file__).resolve())
    out = {"v5_evaluator_sha256": actual, "v5_evaluator_registered": V5_EVAL_SHA256, "failed": {}, "exceptions": [],
           "segment_reasons_not_used": {}, "meta": {}}
    if actual != V5_EVAL_SHA256:
        out["failed"]["共通"] = [f"V5 の評価器の sha256 {str(actual)[:12]} が登録版 {V5_EVAL_SHA256[:12]} でない (理由の読み分けが保証されない)"]
        return out
    bad, exc, metas = EV.preconditions(case, srec5, wall_check or EV.default_wall_check)
    for run, reasons in bad.items():
        keep = [r for r in reasons if not r.startswith(V5_SEGMENT_PREFIX)]
        drop = [r for r in reasons if r.startswith(V5_SEGMENT_PREFIX)]
        if keep:
            out["failed"][run] = keep
        if drop:
            out["segment_reasons_not_used"][run] = drop
    out["exceptions"], out["meta"] = exc, metas
    return out


# --- 窓 B での比較 (V5 と同じ) ---------------------------------------------------------------------------------------
def compare_window(S: dict) -> dict:
    """S = {子: {量: 窓 B の 13 枚}}。V5 の evaluate と同じ統計と判定 (腕 B 対 腕 M、IC 依存 = ISEN 対 腕 M、出口較正)。"""
    rows, verdicts, ic_rows, ic_verdicts = [], {}, [], {}
    for k in QTYS:
        st = {arm: EV.arm_stats([np.asarray(S[r][k], dtype=float) for r in runs]) for arm, runs in ARMS.items()}
        si = EV.arm_stats([np.asarray(S[ISEN][k], dtype=float)])
        D = st["M"]["mean"] - st["B"]["mean"]
        SE = float(np.hypot(st["B"]["se"], st["M"]["se"]))
        Dic = si["mean"] - st["M"]["mean"]
        SEic = float(np.hypot(st["M"]["se"], si["se"]))
        row = {"qty": k, "B": st["B"], "M": st["M"], "D_M_minus_B": float(D), "SE_D": SE, "detected": bool(abs(D) > 2 * SE)}
        irow = {"qty": k, "ISEN": si, "D_ISEN_minus_M": float(Dic), "SE_ic": SEic}
        if k in EV.DQ:
            dq = EV.DQ[k]
            row.update(dq=dq, verdict=EV.judge_one_sided(D, SE, dq), D_plus_2SE_over_dq=float((D + 2 * SE) / dq))
            irow.update(dq=dq, verdict=EV.judge_two_sided(Dic, SEic, dq), absD_plus_2SE_over_dq=float((abs(Dic) + 2 * SEic) / dq))
            verdicts[k] = row["verdict"]
        else:
            row["verdict"] = "記録のみ"
            irow.update(dq=EV.IC_EXIT_M_TOL, signed=True, verdict=EV.judge_two_sided(Dic, SEic, EV.IC_EXIT_M_TOL),
                        absD_plus_2SE_over_dq=float((abs(Dic) + 2 * SEic) / EV.IC_EXIT_M_TOL))
        ic_verdicts[k] = irow["verdict"]
        rows.append(row)
        ic_rows.append(irow)
    ex = next(r for r in rows if r["qty"] == SIGNED_QTY)
    exit_cal = {"qty": SIGNED_QTY, "D_M_minus_B": ex["D_M_minus_B"], "SE_D": ex["SE_D"], "tol": EV.EXIT_CAL_TOL,
                "absD_plus_2SE": abs(ex["D_M_minus_B"]) + 2 * ex["SE_D"], "absD_minus_2SE": abs(ex["D_M_minus_B"]) - 2 * ex["SE_D"],
                "verdict": EV.judge_exit_calibration(ex["D_M_minus_B"], ex["SE_D"]), "note": EV.EXIT_CAL_NOTE}
    return {"rows": rows, "verdicts": verdicts, "ic_rows": ic_rows, "ic_verdicts": ic_verdicts, "exit_calibration": exit_cal}


def _comparison_overall(cmp_: dict, pre_bad: dict) -> tuple:
    """V5 の evaluate と同じ順の総合 (前提 → IC 依存 → 許容幅内 / 悪化 / 保留)。戻り値 (総合, IC 依存の総合)。"""
    ic_overall = EV._ic_overall(cmp_["ic_verdicts"])
    if pre_bad:
        for irow in cmp_["ic_rows"]:
            irow["verdict_if_preconditions_held"] = irow["verdict"]
            irow["verdict"] = EV.IC_REF_LABEL
        ic_overall = "保留 (前提未達: IC の差は参考値で、IC 依存あり/なしを確定しない)"
    v = cmp_["verdicts"]
    if pre_bad:
        overall = "保留 (前提不成立): " + "; ".join(f"{r}: {', '.join(w)}" for r, w in pre_bad.items())
    elif ic_overall != EV.IC_WIN_LABEL:
        overall = f"保留 (IC 依存の前提が成り立たない: {ic_overall})"
    elif all(x == EV.WIN_LABEL for x in v.values()):
        det = [r["qty"] for r in cmp_["rows"] if r.get("detected") and r["qty"] in EV.DQ]
        overall = "許容幅内 (Euler の窓 B で全量の D + 2·SE ≤ Δq)。差を検出した量: " + (", ".join(det) if det else "なし")
    elif any(x == EV.WORSE_LABEL for x in v.values()):
        overall = "悪化 (D − 2·SE ≥ Δq の量: " + ", ".join(k for k, x in v.items() if x == EV.WORSE_LABEL) + ")"
    else:
        overall = "保留 (ユーザ判断: " + ", ".join(k for k, x in v.items() if x != EV.WIN_LABEL) + ")"
    ec = cmp_["exit_calibration"]
    if pre_bad or ic_overall != EV.IC_WIN_LABEL:
        ec["verdict_if_preconditions_held"] = ec["verdict"]
        ec["verdict"] = "保留 (前提不成立・IC 依存の前提; 較正は変えない。較正判断は未完了)"
    ec["decision_complete"] = not ec["verdict"].startswith("保留")
    return overall, ic_overall


def _mark_reference(cmp_: dict) -> None:
    """参考値の表示: 判定の欄を参考値にし、前提成立時の判定は別欄に残す。"""
    for r in cmp_["rows"] + cmp_["ic_rows"]:
        r["verdict_if_registered"] = r["verdict"]
        r["verdict"] = REF_LABEL
    ec = cmp_["exit_calibration"]
    ec["verdict_if_registered"] = ec["verdict"]
    ec["verdict"] = REF_LABEL
    ec["decision_complete"] = False


# --- 判定 ---------------------------------------------------------------------------------------------------------------
def evaluate(case: Path, wall_check=None, v5_static=None) -> dict:
    """§6 V5b の判定。wall_check は nozzle.h5 を読む部分 (V5 の静的な前提)、v5_static は V5 の静的な前提の差し替え (試験用)。"""
    case = Path(case).resolve()
    common = []
    recs = {}
    for name, rel in (("v5", SERIES_REC_V5), ("v5b", SERIES_REC)):
        try:
            recs[name] = _json(case / rel)
        except (OSError, ValueError) as e:
            recs[name] = {}
            common.append(f"時系列の記録 {rel} を読めない ({type(e).__name__})")
    srec5, srec5b = recs["v5"], recs["v5b"]
    if srec5b and srec5b.get("tag") != TAG:
        common.append(f"時系列の記録の tag {srec5b.get('tag')!r} が {TAG!r} でない")
    if srec5b and srec5b.get("geom_ref") != GEOM_REF_PARENT:
        common.append(f"時系列の記録の X_E・X_F の基準 {srec5b.get('geom_ref')!r} が {GEOM_REF_PARENT} (腕 B の r1) でない")
    if srec5 and srec5b:
        if EV.geometry_failures(srec5b.get("geometry")) or srec5b.get("geometry") != srec5.get("geometry"):
            common.append("子の時系列の評価座標 (geometry) が V5 の記録と完全一致しない (同じ評価座標で接続できない)")
        if not srec5b.get("evaluator_sha256") or srec5b.get("evaluator_sha256") != srec5.get("evaluator_sha256"):
            common.append(f"時系列の評価器 eval_wallfit_euler.py の sha256 が V5 と違う ({str(srec5b.get('evaluator_sha256'))[:12]} / "
                          f"{str(srec5.get('evaluator_sha256'))[:12]}; 同じ量の定義で接続できない)")
        if srec5b.get("fixed_coef") != srec5.get("fixed_coef"):
            common.append(f"時系列の fixed_coef が V5 と違う ({srec5b.get('fixed_coef')!r} / {srec5.get('fixed_coef')!r})")
    try:
        ref_prov = provenance(case / GEOM_REF_PARENT)
    except Exception as e:  # noqa: BLE001
        ref_prov = None
        common.append(f"{GEOM_REF_PARENT} の RUN_PROVENANCE を読めない: {type(e).__name__}: {e}")
    runs, settings_bad, data_bad, S = {}, {}, {}, {}
    for child, parent in CHILDREN.items():
        r = {"parent": parent, "role": ROLE[child]}
        runs[child] = r
        if not (case / child).is_dir() or not (case / parent).is_dir():
            data_bad[child] = [f"run dir が無い ({child if not (case / child).is_dir() else parent})"]
            r["label"] = HOLD
            continue
        try:
            s_why, d_why, meta = child_checks(case, child, parent, srec5, srec5b, ref_prov)
        except Exception as e:  # noqa: BLE001
            s_why, d_why, meta = [f"前提を読めない: {type(e).__name__}: {e}"], [], {}
        r["meta"] = meta
        # 接続と分類 (読めない・接続できない run は保留。0 で埋めない)
        try:
            steps, vals = join_series(read_rows(case / parent / f"wallfit_series_{EV.TAG}.csv"),
                                      read_rows(case / child / SERIES_CSV), QTYS)
            r["join"] = {"ok": True, "parent_end_step": PARENT_END, "steps": [int(steps[0]), int(steps[-1])], "n": int(len(steps))}
            qp = case / child / QS_FILE
            cl = classify_run(steps, vals, qp.read_text() if qp.is_file() else None)
            r["classification"] = cl
            S[child] = {q: np.asarray(cl["window_B_values"][q], dtype=float) for q in QTYS}
        except Exception as e:  # noqa: BLE001
            r["join"] = {"ok": False, "reason": f"{type(e).__name__}: {e}"}
            d_why.append(f"時系列を接続・分類できない: {type(e).__name__}: {e}")
        if s_why:
            settings_bad[child] = s_why
        if d_why:
            data_bad[child] = d_why
        if s_why or d_why:
            r["label"] = HOLD
        else:
            r["label"] = r["classification"]["label"]
    labels = {c: runs[c]["label"] for c in CHILDREN}
    all_decay = all(v == DECAY for v in labels.values())
    use_b = all_decay and not common                 # 共通の前提 (記録・評価座標・評価器・基準のバイナリ) が不成立なら判定に使わない
    # V5 の静的な前提 (親) と比較の前提
    try:
        st = v5_static(case, srec5) if v5_static else v5_static_preconditions(case, srec5, wall_check)
    except Exception as e:  # noqa: BLE001
        st = {"failed": {"共通": [f"V5 の静的な前提を判定できない: {type(e).__name__}: {e}"]}, "exceptions": [],
              "segment_reasons_not_used": {}, "meta": {}}
    pre_bad = {}
    for run, w in (st.get("failed") or {}).items():
        pre_bad[f"親 {run}"] = list(w)
    if common:
        pre_bad["共通"] = list(common)
    for c in CHILDREN:
        w = settings_bad.get(c, []) + data_bad.get(c, [])
        seg = ((runs[c].get("meta") or {}).get("segment") or {})
        if seg and seg.get("status") not in ("pass", "plateau") and not any("判定区間の残差判定" in x for x in w):
            w = w + [f"子の判定区間の残差判定が {seg.get('status')} ({seg.get('reason') or seg.get('line')}; 比較の前提は pass か既知の plateau)"]
        if w:
            pre_bad[c] = w
    complete = all(c in S for c in CHILDREN)
    comparison = None
    summary = ", ".join(f"{c}: {labels[c]}" for c in CHILDREN)
    if complete:
        comparison = compare_window(S)
        if use_b:
            overall, ic_overall = _comparison_overall(comparison, pre_bad)
            comparison.update(mode="V5b の判定 (窓 B)", overall=overall, ic_overall=ic_overall)
        else:
            # 参考値: 前提成立時の判定 (V5 の手順) も別欄に残す
            ov_if, ic_if = _comparison_overall(comparison, pre_bad)
            _mark_reference(comparison)
            comparison.update(mode=("参考値 (共通の前提が不成立)" if all_decay else "参考値 (7 本すべてが減衰側でない)"),
                              overall=REF_LABEL, overall_if_registered=ov_if, ic_overall=REF_LABEL, ic_overall_if_registered=ic_if)
    if common:
        v5b_overall = ("保留 (前提不成立: 共通 — " + "; ".join(common) + ") / 分類 (参考) " + summary)
    elif not all_decay:
        v5b_overall = ("保留 (V5b: 7 本すべてが減衰側ではない — 窓 B の差は参考値。V5 は保留のまま) / " + summary
                       + ("" if complete else " / 窓 B の値がそろわない run があり参考値も無い"))
    elif not complete:                               # 全 run が減衰側なら窓 B はそろっている (念のため)
        v5b_overall = "保留 (欠損)"
    else:
        v5b_overall = comparison["overall"]
    out = {"plan": PLAN, "plan_registration_commit": PLAN_REG_COMMIT,
           "evaluator": Path(__file__).name, "evaluator_sha256": _sha(Path(__file__).resolve()),
           "v5_evaluator_sha256": st.get("v5_evaluator_sha256", _sha(Path(EV.__file__).resolve())),
           "v5_evaluator_registered": V5_EVAL_SHA256,
           "series_records": {"v5": SERIES_REC_V5, "v5b": SERIES_REC, "geometry": srec5b.get("geometry"),
                              "evaluator_sha256_v5": srec5.get("evaluator_sha256"), "evaluator_sha256_v5b": srec5b.get("evaluator_sha256")},
           "constants": {"parent_end_step": PARENT_END, "ext_steps": EXT_STEPS, "total_end": TOTAL_END, "out_interval": OUT_INT,
                         "skip_after_restart": SKIP_AFTER_RESTART, "window_A": list(WIN_A), "window_B": list(WIN_B),
                         "oscillation_range": list(OSC_RANGE), "n_window": len(WIN_A_STEPS), "n_oscillation": len(OSC_STEPS),
                         "tol": TOL, "dq": EV.DQ, "min_periods": MIN_PERIODS, "min_half_intervals": MIN_HALF_INTERVALS,
                         "quasisteady": {"drift": EV.QS_DRIFT, "osc": EV.QS_OSC, "tail": EV.QS_TAIL, "full": len(WIN_B_STEPS),
                                         "exception_qty": EV.EXIT_EXC_QTY, "exception_full_tol": EV.EXIT_EXC_FULL_TOL}},
           "children": CHILDREN, "arms": ARMS, "isen": ISEN, "labels": labels, "all_decaying": all_decay, "window_B_comparison_is_v5b_result": bool(use_b and complete),
           "runs": runs, "preconditions_failed": pre_bad, "v5_static": {k: st.get(k) for k in ("failed", "exceptions", "segment_reasons_not_used", "meta")},
           "comparison": comparison, "overall": v5b_overall,
           "limits": ["SE は自己相関を補正しない実務の指標で、振動の平均の不確かさには使わない (諮問 Major)",
                      "減衰側と持続振動側の両方に当たる量は判別不能とした (登録はこの場合を決めていない)",
                      "周期の数え方 (転回点のしきい値 h = tol_q、半周期 ≥ 2000 step、周期は同じ種類の転回点の間) は実装時の定義",
                      "窓 B の準定常は moc_v5_euler_eval の関数を 13 枚の並びに当てたもの (step の平行移動に依らない)",
                      "V5 の判定区間 (親の本段) の残差判定は V5b の前提に使わない (V5b の判定区間は子の本段)",
                      "評価座標の感度 (V5 の evaluation_coordinates) は記録していない (子は親と同じ評価座標で、差の上限だけを検査)",
                      "36000 step は診断の予算であって収束時間の予測ではない。再延長はしない"]}
    (case / "_band_ab").mkdir(exist_ok=True)
    (case / OUT_JSON).write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str))
    return out


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "pairs":
        for c, p in CHILDREN.items():
            print(c, p, ROLE[c])
        return 0
    if argv and argv[0] in ("prep-child", "verify-child") and len(argv) in (3, 4):
        dry = len(argv) == 4 and argv[3] == "--dry"
        if len(argv) == 4 and not dry:
            print(__doc__)
            return 2
        if argv[0] == "prep-child":
            r = prep_child(Path(argv[1]), Path(argv[2]), dry=dry)
            print(f"prep {r['child']}: 複製 {len(r['copied'])} 件 {r['copied']}" + (f" / config が参照するが無い {r['referenced_missing']}"
                                                                                     if r["referenced_missing"] else ""))
            return 0
        r = verify_child(Path(argv[1]), Path(argv[2]), dry=dry)
        rf = r.get("restart_field") or {}
        print(f"verify {r['child']}: {r['VERDICT']}  移した量 {rf.get('n_moved')} (保存量 {r.get('n_src_conserved')}: {r.get('src_conserved')})"
              + ("" if r["VERDICT"] == "OK" else "\n  " + "\n  ".join(r["failures"])))
        return 0 if r["VERDICT"] == "OK" else 1
    if len(argv) == 2 and argv[0] == "run":
        return run_child(Path(argv[1]))
    if len(argv) > 1 or (argv and argv[0].startswith("-")):
        print(__doc__)
        return 2
    case = Path(argv[0]).resolve() if argv else HERE
    out = evaluate(case)
    print(f"評価器 {out['evaluator']} sha256 {out['evaluator_sha256']}  (plan の登録 {out['plan_registration_commit']}; "
          f"V5 の評価器 {str(out['v5_evaluator_sha256'])[:12]})")
    short = {q: q.replace("_eta0.1", "") for q in QTYS}
    print(f"{'run':42s} {'分類':10s} " + " ".join(f"{short[q][:14]:>14s}" for q in QTYS))
    for c in CHILDREN:
        r = out["runs"][c]
        by = ((r.get("classification") or {}).get("by_label")) or {}
        print(f"{c:42s} {r['label']:10s} " + " ".join(f"{by.get(q, '-'):>14s}" for q in QTYS))
    for c in CHILDREN:
        cl = out["runs"][c].get("classification")
        if not cl:
            continue
        for q in QTYS:
            x = cl["quantities"][q]
            print(f"  {c[:8]} {short[q][:22]:22s} 窓A {x['window_A']['mean']:+.6g} 窓B {x['window_B']['mean']:+.6g} "
                  f"移動 {x['mean_shift_A_to_B']:+.2e} (tol {x['tol']:.1e}) 準定常B {x['window_B']['quasisteady_status']} "
                  f"周期 {x['oscillation']['n_periods']} → {x['label']}")
    for r, w in out["preconditions_failed"].items():
        print(f"前提不成立 {r}: " + " / ".join(w))
    cmp_ = out["comparison"]
    if cmp_:
        print(f"窓 B の比較 [{cmp_['mode']}]")
        ic = {x["qty"]: x for x in cmp_["ic_rows"]}
        for row in cmp_["rows"]:
            i = ic[row["qty"]]
            print(f"  {row['qty']:28s} D=M−B {row['D_M_minus_B']:+10.2e} 2SE {2 * row['SE_D']:9.2e} {row['verdict']}"
                  f" ({row.get('verdict_if_registered', '')})  | ISEN−M {i['D_ISEN_minus_M']:+10.2e} 2SE {2 * i['SE_ic']:9.2e} {i['verdict']}")
        print(f"  出口較正: {cmp_['exit_calibration']['verdict']}")
    print("総合:", out["overall"])
    print("->", case / OUT_JSON)
    return 0


if __name__ == "__main__":
    sys.exit(main())
