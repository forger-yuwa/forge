#!/usr/bin/env python3
"""V-h370 の投入の台本 (h370_prep.sh・h370_run.sh・h370_run.py・h370_judge.py; plans/active/tooling-nozzle-isothermal-wall-chain.md §5.1 #36f・
§6 V-h370) の手元の試験。forge・変換器・AWS は使わない (場は合成、設定は本物の生成器で作る)。

  A. 登録の値の写し: h370_run.py・h370_judge.py の定数が h370_problems.py・問題 YAML (停止規則・B0・Euler の本段・格子・H)・
     V-tw370 (tw370_run.py のバイナリ・環境変数・numpy) と同じ。台本の構文 (bash -n)、h370_run.sh の環境変数が tw370_run.sh と同じ。
  B. 壁の出力の検査 (wall_output_check) と NS の段の設定の検査 (ns_config_checks) を、問題 YAML から本物の生成器 (_bcond・_ns_configs) で
     作った bcond・設定に当てる: 合格 (壁 wall_isothermal・outputHDFflg 1・Ts 370、本段の output.level 1、wallTreatmentSST 0、間隔 2,500)。
     負例 (outputHDFflg 0・output.level 0・wallTreatmentSST 1・prepare_info の sha の不一致) で不合格。
  C. Euler の段の設定の導出 (euler_stage_configs) が runner_axismach.run_staged が実際に起動する設定 (soft・mid・本段) と sha で一致。
     Observer が導いた段の外の設定を拒む (forge を起動しない)。
  D. Euler の予備の資格 (h370_judge.py --euler-prequalify、§5.1 #36g (1)・#36i (M1)) を合成の Euler の run で: 合格 → 固定 (res_60000・sha)・
     σ_ref の代わりの参照 (末尾の窓の全 20 出力 40,000〜59,000 と sha)・旧い代表 (5,000 ごと) と tail_field_dev は記録だけ・スクリーニングの
     出力、2 回目は同じ固定を使う。負例: 出口 M・流量のドリフト (DRIFTING)、場の NaN、残差の停滞 (check_convergence NOT CONVERGED)、既にある
     違う固定、run の欠け、代わりの参照の出力が消された (旧い代表の 45,000・代表の外の 41,000 のどちらも判定不能)。
  E. 抽出の参照の A/B の分類 (classify_ab)・参照の最終の資格 (reference_final)・1 出力の比較 (ab_pair) の単体の試験 (結果 A・B・判定不能の
     境界、§5.1 #36g (3)〜(5))。codex の判定の基準の諮問の合成の反例 3 つ (中央 1 列の欠測・停止の 1 行だけ・途中の 1 出力だけの参照の感度
     0.030 mm) が「結果 A」にならず理由つきで拒まれる。最終の資格 (§5.1 #36i): 判定不能 (U_time を含む) は「あり」にしない (M2)、
     σ_ref,use > 0.01 mm は A/B の分類に関係なく資格なし (M3)、σ_ref,use を測れなければ資格なし。
  F. NS の判定 (h370_judge.py) を合成の NS・Euler の run (本物の抽出 stop_rule.delta_snapshot・extract_and_merge、関門と次の壁は偽物) で:
     結果 A (腕 A = 範囲の短い同じ参照: 比べる列は座標と範囲から先に固定、範囲の外の列は除く; σ_ref は判定窓の全出力 × 末尾の窓の全 20 出力、
     σ_ref,use とその位置も記録; 並列 3 と直列 1 で記録がビット同一)、
     結果 B (腕 A の参照の形を変える)、判定不能 (σ_ref が大きい → 参照の最終の資格なし → 更新に進めない; 代わりの参照の sha が違う →
     σ_ref を測れない)、行が欠けても res があれば抽出し直して補う、反例 3 つの組み込みの版、比べる列の固定に出口の幾何の許容が効く、
     関門で止まる、本段の直前の照合の記録が無い → 判定不能。結果 A 以外の判定の多くは試験の時間のため Euler の出力の間隔を 5,000 に替えて
     代わりの参照を 4 つにする (登録の 20 出力は A・D・結果 A・反例・並列で確かめる)。代わりの参照の sha が違う → 抽出の前に止める。
  J. 投入前のレビューの諮問 (2026-10-14-h370-prelaunch-review-diagnose.md) の反例 4 つと判別 A/B を、判定器の全経路 (ab_compare →
     classify_ab → reference_final → update_conditions_met) に偽物の抽出 (速い) で通す: NS 82,500 × Euler 42,000 だけ 0.030 mm・D が窓の中で
     0.008/0.012 で σ_ref 0.030 (と σ 0 の U_time だけの版)・比べる列の外の使用域だけ σ_ref,use 0.030・全出力で D 0.100 と σ_ref 0.030 —
     すべて更新不可。判別 A/B: 同じ入力で代わりの参照を 5,000 ごとの代表 4 つにすると σ_ref は 0 (取りこぼし; 旧い判定器なら結果 A・資格あり)
     だが、今の判定器は代わりの参照が全 20 出力でないので σ_ref を測れず更新不可、全 20 出力では σ_ref 0.030 で更新不可。
  K. プロセス並列の抽出 (run_snaps) が直列とビット同一 (δ_E・δ_target・held・missing・例外の扱い)、並列数の決め方 (snap_jobs)。
  G. h370_run.py の段取り (--dry-run・本番) を偽の準備 (本物の生成器の設定・偽の格子) で: 投入前の検査の合格・不成立・壁の出力が無い
     (終了コード 4)・既存の run で止まる・dry-run と不成立で forge を回さない・資格未達で NS に進まない・環境変数と固定の欠けで止まる・
     h370_problems.py check を作業ツリーの版で回す。
  L. 診断の経路 (§5.1 #36l・§6 V-h370-diag; h370_run.py --phase ns-diag・h370_judge.py --diag-sens、forge なし):
     (a) 診断の経路・感度の判定・診断の run への正式の判定のどれでも受け渡しの関門 (HOOKS と deltastar_loop.handoff_gate を例外を投げる偽物に
     差し替え) が呼ばれない、正式の判定 (h370_judge.py の既定) も呼ばない; (b) 正式の記録 (予備の資格・固定) が書き換わらない (バイトで照合)・
     診断の参照の固定は別のファイル; (c) 正式の --phase ns は変わらない (G・H のほか、診断の記録があっても read_pin・phase_ns は資格 FAIL で
     止まる)、診断の経路は資格ありの記録・正式の固定・参照の sha の違いで止まる; (d) 保存の容量の不足で forge を回さない (見積もりの式も照合);
     (e) 感度の分類 (結果 A・B・判定不能、17 出力の欠け、合否の変化、D が 0.01 mm をまたぐ時刻があっても最大で決まる) を単体・偽物の抽出の
     全経路・本物の抽出 (diag_snapshot は stop_rule.delta_snapshot と同じ値) で; (f) 途中の出力を消さない (診断の経路の後に本段の全出力と
     境界の出力が残り sha の一覧に入る、stop_rule.delete_after は停止の出力より後だけを消す、診断の関数・run_staged_ns の本段の部分に消す処理が無い)。
  M. 帯を共通にした A/B (§5.1 #36n; h370_judge.py --diag-sens-fixed-band、forge なし): 分類 (classify_fixed: 結果 A・B・判定不能、0.01 mm の境目、
     F の中の held は有効な列の差が大きくても判定不能、腕 A の再現の不成立は判定不能)・区間ごとの最大 (section_max) を単体で; 偽物の抽出の全経路で
     結果 A・B・判定不能 (held・抽出の不成立)、腕 A の再現の不成立 (1e−9 r_t・非有限) で判定不能かつ腕 B を回さない、丸め程度 (5e−13 r_t) は再現、
     帯の固定が両腕に同じ値で渡る (腕 A に #36l の npz の A_y_b、腕 B に腕 A が使った全列の帯)、#36l の記録との照合・#36l と同じ照合の不成立で
     抽出しない; 本物の抽出 (diag_snapshot_fixed) で、自分の帯を固定すると diag_snapshot とビット同一・使用域の外は自分の帯・腕 B に腕 A の帯が
     そのまま入る・差し替えが元に戻る、本物の全経路 (L の #36l の記録) で腕 A の再現がビット同一・34 抽出・両腕の帯が同じ。
  N. 帯の時間方向の扱いの A/B (§5.1 #36p; h370_judge.py --diag-time-band、forge なし): 分類 (classify_time: 結果 A = D_time > 0.01 mm・結果 B・
     判定不能 (腕 A の再現・先頭の時刻の再現・held・例外・出力の欠け)) を単体で; 梯子の交差の区間の写し (ladder_crossings) を本物の
     deltastar.edge_band_positions の y* と照合 (途中の交差・最初の点・交差なし・点の不足・区間 2 つ); 偽物の抽出の全経路で結果 B (両腕の時間変化が
     同じ)・結果 A (腕 B だけ x 3〜5 で 0.02 mm)・判定不能 (held・腕 A の非有限と 1e−9 r_t の再現の不成立 (腕 B を回さない)・先頭の時刻の再現の不成立・
     帯の固定の不成立・例外・入力の照合)、帯の固定が全時刻に同じ値 (腕 B の 17 回すべてに k₀ の腕 A の全列の帯)・2 窓・代表の列の check_quasisteady・
     梯子の変化の数え方; 本物の抽出で record_ladder が値を変えず梯子の写しが全列で合う、本物の全経路で腕 A の再現・先頭の時刻の再現がビット同一。
  O. 受け渡しの関門の局所 y₁⁺ (§5.1 #36q; deltastar_loop.wall_resolution_verdict と本物の check_wall_resolution.py、forge なし) を 2 枚の
     no-slip 壁と出口の出力を持つ合成の run (solver_density_cuda/tools/test_gate_bad_input.make_wall_run) で: 事前登録の受入 — 対照 (両壁とも
     停止の step に有効な出力) は PASS (--step s・--groups は bcond の wall* の群だけ・--weight area)、片方の壁の停止の step の出力だけを欠くと
     UNDETERMINED (ツールを回さない)、停止の step が無い (古い step だけ) と UNDETERMINED、面積と点数で合否が逆になる入力は面積で判定 (両向き)、
     出口の出力があっても壁だけで判定 (群なしのツールは出口で判定不能になる入力)。ほか: 壁の出力が 1 つも無い壁・要る量 (twall) の欠け・
     属性 step の時刻違い・停止の step なし・no-slip 壁なし・bcondConfig なし → UNDETERMINED、ツールの壁ごとの行の欠け・step 違い → 終了コード 0 でも
     UNDETERMINED、受け渡しの関門 (handoff_gate、他の条件は偽物で合格) で UNDETERMINED は 3 つの上書きでも止まり対照は止めない。run_0951 の
     停止の出力の形の合成 (壁 4,719 点・超過 205 点・面積の割合 2.932 %・出口の出力・古い step 97,500) は FAIL・面積の割合を写す (点数なら
     4.344 %)。環境変数 H370_RUN0951_DIR に run_0951 の停止の出力の写し (bcondConfig・solverConfig・nozzle.h5・res_152500.h5・
     res_wall_3_152500.h5・res_outlet_2_152500.h5) があれば、それも新しい関門に通して FAIL を確かめる。
usage: /home/sano/work/forge/design/.venv-opt/bin/python case/45.isobutane_m6_d155/test_h370_launch.py   (所要 約 8〜12 分; 終了コード 0 = ALL PASS)
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import types
from pathlib import Path

import h5py
import numpy as np

CASE = Path(__file__).resolve().parent
ROOT = CASE.parents[1]
sys.path.insert(0, str(ROOT / "design"))
sys.path.insert(0, str(CASE))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
os.environ.setdefault("FORGE_NO_DESIGN_REPORT", "1")

import h370_judge as HJ  # noqa: E402
import h370_run as HR  # noqa: E402
from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.evaluate import stop_rule as SR  # noqa: E402
from forge_design.feedback import deltastar_loop as DL  # noqa: E402
from forge_design.metrics import deltastar as DS  # noqa: E402

# プロセス並列の抽出 (h370_judge.run_snaps、spawn) の子はこの台本を __mp_main__ として読み直すので、子では一時ディレクトリを作らない
TMP = Path(tempfile.mkdtemp(prefix="h370_launch_")) if __name__ != "__mp_main__" else None
FAIL = 0
TOOLS = ROOT / "solver_density_cuda" / "tools"


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


def raises(fn, exc=Exception):
    try:
        fn()
    except exc as e:  # noqa: BLE001
        return f"{type(e).__name__}: {e}"
    return None


class Patch:
    """属性の差し替え (with で戻す)。"""

    def __init__(self, *items):
        self.items, self.old = items, []

    def __enter__(self):
        for obj, name, val in self.items:
            self.old.append((obj, name, getattr(obj, name)))
            setattr(obj, name, val)
        return self

    def __exit__(self, *a):
        for obj, name, val in reversed(self.old):
            setattr(obj, name, val)


# ====================================================================================== A. 登録の値の写し
def test_constants():
    print("== A. 登録の値の写し・台本の構文", flush=True)
    import h370_problems as HP
    import tw370_run as TR
    check("δ_r の CSV の sha256・停止規則の閾値と断面・Euler の上限と throat_refine・壁温が h370_problems.py と同じ",
          HR.CSV_SHA == HP.CSV_SHA and HR.THRESHOLDS == HP.STOP_RULE_THRESHOLDS and HR.X_FRACS == HP.STOP_RULE_X_FRACS
          and HR.EU_STEPS == HP.EU_NSTEP and HR.EU_MESH["throat_refine"] == HP.EU_THROAT_REFINE and HR.TW == HP.TW_NEW
          and HR.NS_STAGE_WANT["main"]["nStepOuter"] == 200000)
    pn, pe = RA.load_problem(CASE / HR.NS_PROBLEM), RA.load_problem(CASE / HR.EU_PROBLEM)
    ml = RA._main_line_spec(pn)
    ev = pn.evaluate
    check("NS の問題 YAML: 停止規則 (上限・保持・閾値)・出力の間隔・B0・格子 (Gref 4719 × 160)・H の上流・壁温が台本の定数と同じ",
          ev["stop_rule"] is True and ev["stop_rule_cap"] == HR.NS_CAP and ev["stop_rule_hold_steps"] == HR.NS_HOLD
          and ev["outStepInterval"] == HR.NS_OUT and ev["nStepOuter"] == HR.NS_CAP and ev["stop_rule_thresholds"] == HR.THRESHOLDS
          and {k: ml[k] for k in HR.MAIN_LINE} == HR.MAIN_LINE and (pn.mesh["ni"], pn.mesh["nj"]) == (HR.NS_MESH["ni"], HR.NS_MESH["nj"])
          and HR.phys_up_ok({"physical_upstream": {"requested": RA._physical_upstream(pn.geometry)}})
          and pn.wall_thermal == {"mode": "isothermal", "Tw": HR.TW},
          str({k: ml[k] for k in HR.MAIN_LINE}))
    me = RA.mesh_euler_block(pe)
    check("Euler の問題 YAML: 本段 60,000・出力 1,000・cfl_main 2.0・mid_stage・格子 2000 × 97・throat_refine 5・H の上流",
          (pe.evaluate["nStepOuter"], pe.evaluate["outStepInterval"], pe.evaluate["cfl_main"], pe.evaluate["mid_stage"])
          == (HR.EU_STEPS, HR.EU_OUT, HR.EU_CFL, True) and (me["ni"], me["nj"], me["throat_refine"]) == tuple(HR.EU_MESH.values())
          and HR.phys_up_ok({"physical_upstream": {"requested": RA._physical_upstream(pe.geometry)}})
          and pe.geometry["initial_line_run"] == HR.CFD_PIN_RUN)
    check("physical_upstream の照合: 解決済みの値 (x_in・mu・frozen_design 付き) の 3 値で比べる; 違う半径は不合格",
          not HR.phys_up_ok({"physical_upstream": {"requested": {**RA._physical_upstream(pe.geometry), "r_inlet": 6.485}}})
          and not HR.phys_up_ok({}))
    check("バイナリ (forge 7e95e5bbd8f50c9e・変換器 711f4438)・FORGE_* の環境変数・numpy 2.5.2 が V-tw370 (tw370_run.py) と同じ",
          (HR.BIN, HR.BIN_SHA, HR.CONV_SHA, HR.ENV_WANT, HR.NUMPY_WANT) == (TR.BIN, TR.BIN_SHA, TR.CONV_SHA, TR.ENV_WANT, TR.NUMPY_WANT))
    check("h370_judge.py の run の名前・固定と資格の記録の場所・バイナリ・A/B の許容 (関門の 2 窓の予算 0.01 mm) が h370_run.py と同じ",
          (HJ.EU_RUN, HJ.NS_RUN, HJ.NS_PROBLEM, HJ.OLD_EULER, HJ.PIN, HJ.QUAL, HJ.EU_STEPS, HJ.BIN_SHA)
          == (HR.EU_RUN, HR.NS_RUN, HR.NS_PROBLEM, HR.OLD_EULER, HR.PIN, HR.QUAL, HR.EU_STEPS, HR.BIN_SHA)
          and HJ.TOL_MM == DL.HANDOFF_WINDOW_TOL_MM)
    check("σ_ref の代わりの参照の列挙 (§5.1 #36i M1): h370_judge.sigma_ref_required() = h370_run.SIGMA_REF_STEPS = 40,000〜59,000 の 20 出力 "
          "(Euler の出力の間隔・末尾の窓が h370_run.py と同じ)",
          HJ.sigma_ref_required() == HR.SIGMA_REF_STEPS == list(range(40000, 60000, 1000)) and (HJ.EU_OUT, HJ.EU_TAIL) == (HR.EU_OUT, HR.EU_TAIL))
    sh = {n: subprocess.run(["bash", "-n", str(CASE / n)], capture_output=True, text=True) for n in ("h370_prep.sh", "h370_run.sh")}
    check("h370_prep.sh・h370_run.sh の構文 (bash -n)", all(r.returncode == 0 for r in sh.values()), str({k: v.stderr for k, v in sh.items()}))
    env = lambda p: [ln.strip() for ln in (CASE / p).read_text().splitlines() if ln.strip().startswith(("export FORGE", "for v in $(compgen"))]  # noqa: E731
    check("h370_run.sh の FORGE_* の設定 (export・unset) が tw370_run.sh と同じ行", env("h370_run.sh") == env("tw370_run.sh") and len(env("h370_run.sh")) == 4,
          str(env("h370_run.sh")))


# ====================================================================================== B. 壁の出力・NS の段の設定
def gen_ns_configs(rd: Path):
    """問題 YAML から prepare_ns と同じ生成器で bcond・前の段の元の設定・本段の設定を作る (格子は偽物)。"""
    p = RA.load_problem(CASE / HR.NS_PROBLEM)
    ml = RA._main_line_spec(p)
    n, out_int, _ = RA._ns_steps(p, None)
    rd.mkdir(parents=True)
    (rd / "bcondConfig.yaml").write_text(RA._bcond_with_species(RA._bcond(p, euler=False, profile_ok=True), RA._tp_species_Y(p)))
    cfg_ns, cfg_main = RA._ns_configs(p, rd, n, out_int, None, None, ml)
    (rd / "solverConfig.yaml").write_text(cfg_main)
    (rd / RA.PRESTAGE_CONFIG).write_text(cfg_ns)
    with h5py.File(rd / "nozzle.h5", "w") as f:              # 偽の格子 (mesh_identity が読む最小の中身)
        f["MESH/COORD"] = np.arange(12.0)
        f["VALUE/wall_dist"] = np.ones(4)
        f["VALUE/ro"] = np.ones(4)
    (rd / "nozzle.msh").write_text("$MeshFormat 偽物\n")
    (rd / "probe.yaml").write_text("")
    info = {"ns_stage_configs": RA.ns_stage_sha256(cfg_main, cfg_ns), "main_line": {**ml, "source": "evaluate.main_line"},
            "delta_r_csv": str(HR.HERE / HR.CSV), "physical_upstream": {"requested": RA._physical_upstream(p.geometry), "wall": {}},
            "wall_thermal": dict(p.wall_thermal), "mesh": {"ni": 4719, "nj": 160, "hashes": {"generator_sha256": "x" * 64}},
            "initial_line": {"match": {"checked": ["R", "L_U"], "unchecked": []}}}
    (rd / "prepare_info.json").write_text(json.dumps(info))
    return p, info


def test_wall_output_and_configs():
    print("== B. 壁の出力の検査・NS の段の設定の検査 (本物の生成器の設定)", flush=True)
    rd = TMP / "cfg_ok"
    gen_ns_configs(rd)
    w = HR.wall_output_check(rd)
    check(f"壁の出力: 合格 ({w.get('wall_dump_name')}・output.level {w.get('output_level')}・間隔 {w.get('outStepInterval')}・wallTreatmentSST "
          f"{w.get('wallTreatmentSST')})", w["ok"] and w.get("wall_dump_name") == "res_wall_3_<step>.h5" and w.get("output_level") == 1
          and w.get("outStepInterval") == 2500 and w.get("wallTreatmentSST") == 0, str(w.get("failures")))
    c = HR.ns_config_checks(RA, rd)
    check("NS の段の設定: 合格 (S1 0.5・S2 1.0 は 1 次・3,000 step・B0 のキーなし、本段 cfl 4・200,000・2,500・B0、3 段の implicitRelax 0.7、"
          "prepare_info の sha と一致)", c["ok"] and c["stages"]["main"]["lineDtDirectionalCap"] == 50.0 and c["stages"]["S1_soft"]["lineImplicit"] is None
          and all(c["prepare_info_match"].values()), str(c.get("failures")))
    for tag, fn, want in (("outputHDFflg 0", lambda d: (d / "bcondConfig.yaml").write_text((d / "bcondConfig.yaml").read_text().replace(
                              "kind: wall_isothermal,  outputHDFflg: 1", "kind: wall_isothermal,  outputHDFflg: 0")), "outputHDFflg"),
                          ("output.level 0", lambda d: (d / "solverConfig.yaml").write_text((d / "solverConfig.yaml").read_text().replace(
                              "output: {level: 1,", "output: {level: 0,")), "output.level"),
                          ("wallTreatmentSST 1", lambda d: (d / "solverConfig.yaml").write_text((d / "solverConfig.yaml").read_text().replace(
                              "wallTreatmentSST: 0", "wallTreatmentSST: 1")), "wallTreatmentSST")):
        d = TMP / f"cfg_bad_{tag.split()[0]}"
        shutil.copytree(rd, d)
        before = (d / "bcondConfig.yaml").read_text() + (d / "solverConfig.yaml").read_text()
        fn(d)
        changed = before != (d / "bcondConfig.yaml").read_text() + (d / "solverConfig.yaml").read_text()
        r = HR.wall_output_check(d)
        check(f"壁の出力: 負例 {tag} → 不合格", changed and not r["ok"] and any(want in f for f in r["failures"]), str(r.get("failures")))
    d = TMP / "cfg_bad_sha"
    shutil.copytree(rd, d)
    info = json.loads((d / "prepare_info.json").read_text())
    info["ns_stage_configs"]["main"] = "0" * 64
    (d / "prepare_info.json").write_text(json.dumps(info))
    r = HR.ns_config_checks(RA, d)
    check("NS の段の設定: prepare_info の本段の sha が違う → 不合格", not r["ok"] and not r["prepare_info_match"]["main"])
    r = HR.ns_prepare_info_checks(json.loads((rd / "prepare_info.json").read_text()), None, prelaunch=True)
    check("NS の prepare_info (投入前の準備): 合格 (δ_r の CSV・H・壁温・B0・格子・CFD ピンの照合)", r["ok"], str(r["failures"]))
    bad = {**json.loads((rd / "prepare_info.json").read_text()), "wall_thermal": {"mode": "isothermal", "Tw": 300.0}}
    r = HR.ns_prepare_info_checks(bad, None, prelaunch=True)
    check("NS の prepare_info: 壁温 300 K → 不合格", not r["ok"] and any("370" in f for f in r["failures"]))


# ====================================================================================== C. Euler の段の設定
def test_euler_stage_configs():
    print("== C. Euler の段の設定の導出 = run_staged が起動する設定", flush=True)
    pe = RA.load_problem(CASE / HR.EU_PROBLEM)
    rd = TMP / "eu_cfg"
    rd.mkdir()
    cfg = RA._apply_gas_to_config(RA._config_euler_node(pe, HR.EU_STEPS, HR.EU_OUT, 4.0, 1), pe, rd)
    (rd / "solverConfig.yaml").write_text(cfg)
    (rd / "bcondConfig.yaml").write_text(RA._bcond_with_species(RA._bcond(pe, euler=True), RA._tp_species_Y(pe)))
    want = HR.euler_stage_configs(RA, cfg, HR.EU_CFL, True)
    got = []

    def fake_forge(run_dir):
        t = (Path(run_dir) / "solverConfig.yaml").read_text()
        got.append(HR.sha256_text(t))
        (Path(run_dir) / f"res_{RA._cfg_value(t, 'nStepOuter')}.h5").write_bytes(b"")
        return 0
    stub_res = lambda *a, **k: {"ok": True, "failures": []}  # noqa: E731
    with Patch((RA, "run_forge", fake_forge), (RA, "stage_gate", lambda *a: []), (RA, "_restart_same_mesh", lambda *a: None), (HR, "resources", stub_res)):
        tags = {HR.sha256_text(c): t for t, c in want.items()}
        obs = HR.Observer(RA, rd, tags, "main")
        obs.install()
        try:
            rc = RA.run_staged(rd, cfl_main=HR.EU_CFL, mid_stage=True, stages="full")
        finally:
            obs.uninstall()
        check("run_staged が起動した soft・mid・本段の設定の sha = euler_stage_configs の導出 (Observer は 3 段とも通す)",
              rc == 0 and got == [HR.sha256_text(want[t]) for t in ("soft", "mid", "main")], f"{len(got)} 段")
        ev = [json.loads(ln) for ln in (rd / "h370_stages.jsonl").read_text().splitlines()]
        check("Observer の記録: 段の start・end と本段の前の資源", [e.get("tag") for e in ev if e["event"] == "start"] == ["soft", "mid", "main"]
              and any(e["event"] == "resources" for e in ev))
        rd2 = TMP / "eu_cfg2"
        shutil.copytree(rd, rd2, ignore=shutil.ignore_patterns("res_*", "h370_stages.jsonl", "stage_manifest.json", "residual_history_*"))
        (rd2 / "solverConfig.yaml").write_text(cfg)
        got.clear()
        obs = HR.Observer(RA, rd2, {k: v for k, v in tags.items() if v != "mid"}, "main")
        obs.install()
        try:
            e = raises(lambda: RA.run_staged(rd2, cfl_main=HR.EU_CFL, mid_stage=True, stages="full"), RuntimeError)
        finally:
            obs.uninstall()
        check(f"Observer: 導いた段の外の設定 (mid を外した表) → forge を起動しない ({(e or '')[:40]}…)", e is not None and len(got) == 1)
    st = {t: {k: RA._cfg_value(c, k) for k in ("convMethod", "cfl", "nStepOuter", "outStepInterval")} for t, c in want.items()}
    check("Euler の段: soft は 1 次・cfl 0.5・3,000、mid は cfl 1.0・3,000、本段は cfl 2.0・60,000・出力 1,000",
          st["soft"]["convMethod"] == 0 and st["soft"]["cfl"] == 0.5 and st["soft"]["nStepOuter"] == 3000 and st["mid"]["cfl"] == 1.0
          and st["mid"]["nStepOuter"] == 3000 and (st["main"]["cfl"], st["main"]["nStepOuter"], st["main"]["outStepInterval"]) == (2.0, 60000, 1000),
          str(st))


# ====================================================================================== D. Euler の資格
EU_NI, EU_NJ, EU_S = 80, 9, 0.0767


def make_eu_run(case: Path, drift: float = 0.0, nan_step=None, flat_resid: bool = False, steps=None):
    """合成の Euler の run (構造格子 80 × 9、本段の出力 1,000〜60,000、段の残差の履歴と stage_manifest、RUN_RC 0、投入前の検査の合格)。"""
    from stage_manifest import StageManifest
    run = case / HJ.EU_RUN
    run.mkdir(parents=True)
    x = np.linspace(-12.5, 95.3, EU_NI)
    rw = np.where(x < 0, 1.0 + 0.017 * x * x, 1.0 + 0.09 * x)
    eta = np.linspace(0.0, 1.0, EU_NJ)
    X, R = np.repeat(x[:, None], EU_NJ, 1), rw[:, None] * eta[None, :]
    with h5py.File(run / "nozzle.h5", "w") as f:
        f["/MESH/COORD"] = np.c_[X.ravel() * EU_S, R.ravel() * EU_S, np.zeros(X.size)].ravel()
    M0 = np.where(x[:, None] < 0, 0.2 + 0.8 * np.exp(x[:, None] / 3.0), 1.0 + 5.0 * (1.0 - np.exp(-np.maximum(x[:, None], 0.0) / 20.0))) \
        * (1.0 - 0.002 * eta[None, :] ** 2)
    for st in (steps or range(1000, 60001, 1000)):
        k = 1.0 + 1e-3 * np.exp(-st / 5000.0) + drift * st / 60000.0
        M = M0 * k
        T = 1600.0 / (1.0 + 0.135 * M ** 2)
        a = np.sqrt(1.27 * 287.0 * T)
        P = 5.5e6 * (T / 1600.0) ** (1.27 / 0.27)
        ro = P / (287.0 * T)
        with h5py.File(run / f"res_{st}.h5", "w") as f:
            for n_, v in (("ro", ro), ("Ux", M * a), ("Uy", 0.0 * M), ("Uz", 0.0 * M), ("P", P), ("T", T), ("sonic", a), ("roe", P / 0.27)):
                v = np.asarray(v, dtype=float).ravel().copy()
                if nan_step == st and n_ == "Ux":
                    v[5] = np.nan
                f[f"/VALUE/{n_}"] = v
    (run / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": EU_NI, "nj": EU_NJ}, "scale_m": EU_S, "x_E": 40.0, "Md": 6.0}))
    cfgs = {"soft": "space: {convMethod: 0, limiter: 2}\ntime:\n  deltaT: {cfl: 0.5, cfl_pseudo: 0.5}\n",
            "mid": "space: {convMethod: 1, limiter: 2}\ntime:\n  deltaT: {cfl: 1.0, cfl_pseudo: 1.0}\n",
            "main": "space: {convMethod: 1, limiter: 2}\ntime:\n  deltaT: {cfl: 2.0, cfl_pseudo: 2.0}\n"}
    sm = StageManifest(run)
    for tag, n in (("soft", 3000), ("mid", 3000), ("main", 60000)):
        s = np.arange(1, n + 1, 10)
        r = np.full(len(s), 1e-3) if (flat_resid and tag == "main") else 1e-2 * 10.0 ** (-5.0 * s / n)
        with open(run / f"residual_history_{tag}.csv", "w") as fh:
            fh.write("step,rms_ro,rms_roUx,rms_roUy,rms_roUz,rms_roe\n")
            for si, ri in zip(s, r):
                fh.write(f"{si},{ri:.6e},{2 * ri:.6e},{0.5 * ri:.6e},0,{3 * ri:.6e}\n")
        sm.add(tag, cfgs[tag], "wall: {kind: slip}\n", history=f"residual_history_{tag}.csv")
    sm.write()
    (run / "RUN_RC").write_text("0\n")
    (run / "H370_PRECHECK.json").write_text(json.dumps({"ok": True}))
    return run


def qualify(case: Path) -> int:
    return HJ.main(["--euler-prequalify", "--case-dir", str(case), "--design", str(ROOT / "design"), "--tools", str(TOOLS)])


def test_euler_qualify():
    print("== D. Euler の予備の資格 (合成の run)", flush=True)
    case = TMP / "q_ok"
    run = make_eu_run(case)
    rc = qualify(case)
    q, pin = json.loads((case / HJ.QUAL).read_text()), HJ.jload(case / HJ.PIN)
    v = {c: d["verdict"] for c, d in q["quasisteady"]["verdicts"].items()}
    check(f"合格: 予備の資格あり (rc {rc})、資格の量 {len(v)} 列すべて STEADY、check_convergence --segment PASS、記録の名前は予備の資格 (prequalified)",
          rc == 0 and q["prequalified"] and "qualified" not in q and set(v) == set(HJ.QUAL_COLS) and set(v.values()) == {"STEADY"}
          and q["convergence"]["rc"] == 0 and HJ.QUAL.endswith("h370_euler_prequalify.json"), str(q["failures"] + q["undecidable"]))
    check("固定: 本段の最後の出力 res_60000.h5 とその sha256 (h370_euler_ref_pin.json)",
          pin is not None and pin["res"] == "res_60000.h5" and pin["sha256"] == HJ.sha256_file(run / "res_60000.h5") and pin["run"] == HJ.EU_RUN)
    sr = q["sigma_ref_snapshots"]
    check("σ_ref の代わりの参照 (§5.1 #36i M1): 末尾の窓 [40,000, 60,000] の全出力 40,000〜59,000 の 20 出力 (固定の 60,000 を除く)、sha つき; "
          "末尾の窓の全 21 出力の sha も記録",
          [o["step"] for o in sr] == list(range(40000, 60000, 1000)) and all(o["sha256"] == HJ.sha256_file(run / o["res"]) for o in sr)
          and [o["step"] for o in q["tail_outputs"]] == list(range(40000, 60001, 1000)) and q["sigma_ref_selection"]["required"] == list(range(40000, 60000, 1000))
          and q["sigma_ref_selection"]["tail_window"] == [40000, 60000] and "340 抽出" in q["sigma_ref_selection"]["why"], str([o["step"] for o in sr]))
    td = q.get("tail_field_dev") or {}
    check("旧い代表 (5,000 ごと 40,000・45,000・50,000・55,000) と末尾の窓の 20 出力の場の差 (tail_field_dev) は記録だけ",
          [o["step"] for o in q["sigma_ref_representatives_record"]] == [40000, 45000, 50000, 55000] and len(td.get("per_output", [])) == 20
          and td.get("max_rel_representatives") is not None and td.get("max_rel_others") is not None, str(td)[:200])
    check("スクリーニングの出力: check_quasisteady の末尾の窓 (既定 0.4 → 24 出力) の始まり 37,000 と中ほど 49,000、sha つき",
          [s_["step"] for s_ in q["screening_snapshots"]] == [37000, 49000] and q["quasisteady_tail_window"] == [37000, 60000]
          and all(len(s_["sha256"]) == 64 for s_ in q["screening_snapshots"]))
    check("記録だけの量 (nozzle_report.metrics の η 0.1 の線の wave・overshoot・slope・range) の VERDICT も残る",
          set(q["quasisteady_record"]["record_cols (既定の閾値)"]["verdicts"]) == set(HJ.REC_COLS))
    rc2 = qualify(case)
    check("2 回目: 同じ固定をそのまま使う (rc 0、固定のファイルは変えない)", rc2 == 0 and HJ.jload(case / HJ.PIN) == pin
          and json.loads((case / HJ.QUAL).read_text()).get("pin_note") == "既にある同じ固定をそのまま使う")
    for name, kw, want_rc, want in (("drift", {"drift": 0.2}, 1, "STEADY"), ("nan", {"nan_step": 45000}, 1, "非有限"),
                                    ("flat", {"flat_resid": True}, 1, "check_convergence")):
        c = TMP / f"q_{name}"
        make_eu_run(c, **kw)
        rc = qualify(c)
        q = json.loads((c / HJ.QUAL).read_text())
        check(f"負例 {name}: 未達 (rc {rc})・固定を書かない", rc == want_rc and not q["prequalified"] and not (c / HJ.PIN).exists()
              and any(want in f for f in q["failures"]), str(q["failures"])[:300])
    c = TMP / "q_conflict"
    shutil.copytree(case, c)
    p = json.loads((c / HJ.PIN).read_text())
    p["sha256"] = "0" * 64
    (c / HJ.PIN).write_text(json.dumps(p))
    rc = qualify(c)
    check("既にある違う固定 → 判定不能 (rc 2)・固定を消さない", rc == 2 and json.loads((c / HJ.PIN).read_text())["sha256"] == "0" * 64)
    c = TMP / "q_missing"
    c.mkdir()
    check("run が無い → 判定不能 (rc 2)", qualify(c) == 2)
    c = TMP / "q_short"
    make_eu_run(c, steps=range(1000, 50001, 1000))
    rc = qualify(c)
    check("本段が上限 60,000 まで回っていない → 未達 (rc 1)", rc == 1 and any("上限" in f for f in json.loads((c / HJ.QUAL).read_text())["failures"]))
    c = TMP / "q_gap"
    make_eu_run(c, steps=[s_ for s_ in range(1000, 60001, 1000) if s_ != 45000])
    rc = qualify(c)
    q = json.loads((c / HJ.QUAL).read_text())
    check("σ_ref の代わりの参照の出力 (45,000) が消されている → 判定不能 (rc 2)・固定を書かない", rc == 2 and not (c / HJ.PIN).exists()
          and any("45000" in u for u in q["undecidable"]), str(q["undecidable"])[:200])
    c = TMP / "q_gap41"
    make_eu_run(c, steps=[s_ for s_ in range(1000, 60001, 1000) if s_ != 41000])
    rc = qualify(c)
    q = json.loads((c / HJ.QUAL).read_text())
    check("§5.1 #36i M1: 旧い代表の外の出力 (41,000) が消されている → 判定不能 (rc 2)・固定を書かない (旧い判定器は代表だけを見て通した)",
          rc == 2 and not (c / HJ.PIN).exists() and any("41000" in u for u in q["undecidable"]), str(q["undecidable"])[:200])


# ====================================================================================== E. A/B の分類
def test_classify_ab():
    print("== E. 抽出の参照の A/B の分類・参照の最終の資格・1 出力の比較 (単体)", flush=True)
    req2 = [0, 1]
    per = lambda Ds, pf=None: [{"step": i, "D_mm": d, "pf_A": True if pf is None else pf[i], "pf_B": True} for i, d in enumerate(Ds)]  # noqa: E731

    def sig(mm, alts_pf=None, status="ok"):
        """σ_ref の記録の形 (per_step の alts の pf_F は代わりの参照の時刻の腕 B の合否)。"""
        ps = None if alts_pf is None else {k: {"alts": [{"res": f"res_{e}.h5", "pf_F": v} for e, v in enumerate(alts_pf)]} for k in req2}
        return {"status": status, "mm": mm, "per_step": ps, "why": "試験"}
    cases = [("全出力 0.003 mm・σ 0.001 → 結果 A", per([0.003, 0.002]), 0.001, "結果 A", "A"),
             ("全出力 0.02 mm・σ 0.001 → 結果 B", per([0.02, 0.03]), 0.001, "結果 B", "B_D"),
             ("窓の中で 0.01 mm をまたぐ (0.008・0.012) → 判定不能 (時間の変動)", per([0.008, 0.012]), 0.0, "判定不能", "U_time"),
             ("D 0.008・σ 0.005 → 判定不能 (参照の未整定)", per([0.008, 0.008]), 0.005, "判定不能", "U_sigma_straddle"),
             ("D 0.012・σ 0.005 → 判定不能 (参照の未整定)", per([0.012, 0.013]), 0.005, "判定不能", "U_sigma_straddle"),
             ("合否が変わる出力・代わりの参照の全部で違いが保たれる → 結果 B (D が小さくても)",
              [{"step": 0, "D_mm": None, "pf_A": False, "pf_B": True}, {"step": 1, "D_mm": 0.001, "pf_A": True, "pf_B": True}],
              sig(0.0, [True, True]), "結果 B", "B_pf"),
             ("合否が変わる出力・代わりの参照で確かめられない (per_step が無い) → 判定不能",
              [{"step": 0, "D_mm": None, "pf_A": False, "pf_B": True}, {"step": 1, "D_mm": 0.001, "pf_A": True, "pf_B": True}],
              0.0, "判定不能", "U_pf_unconfirmed"),
             ("合否が変わる出力・代わりの参照の 1 つで違いが消える → 判定不能",
              [{"step": 0, "D_mm": None, "pf_A": False, "pf_B": True}, {"step": 1, "D_mm": 0.001, "pf_A": True, "pf_B": True}],
              sig(0.0, [True, False]), "判定不能", "U_pf_unconfirmed"),
             ("両腕とも比べる列の中に欠測 (合否は同じ) → 判定不能 (比較を合格にしない)",
              [{"step": 0, "D_mm": None, "pf_A": False, "pf_B": False}, {"step": 1, "D_mm": 0.001, "pf_A": True, "pf_B": True}],
              0.0, "判定不能", "U_fixed_cols"),
             ("ちょうど 0.01 mm (max + σ) → 結果 A (≤)", per([0.009, 0.009]), 0.001, "結果 A", "A"),
             ("抽出の失敗 → 判定不能", per([0.001]) + [{"step": 1, "error_A": "x"}], 0.0, "判定不能", "U_error"),
             ("σ を測れない → 判定不能", per([0.001, 0.001]), None, "判定不能", "U_sigma_unmeasurable"),
             ("σ が判定不能 (代わりの参照で合否が変わる) → 判定不能", per([0.001, 0.001]), sig(0.001, status="undecidable"), "判定不能",
              "U_sigma_undecidable"),
             ("出力が無い → 判定不能", [], 0.0, "判定不能", "U_empty")]
    for name, p_, sg, want, code in cases:
        c = HJ.classify_ab(p_, sg, req2)
        check(name, c["result"] == want and c["code"] == code, f"{c['result']} {c['code']}: {c['why'][:120]}")
    c = HJ.classify_ab(per([0.001, 0.001]), 0.0, None)
    check("必要な出力の集合が無い → 判定不能", c["result"] == "判定不能" and c["code"] == "U_required")
    # --- codex の判定の基準の諮問の合成の反例 (2026-10-14-h370-judge-criteria-diagnose.md) ---
    req17 = list(range(22500, 62501, 2500))
    check("必要な出力の列挙: 停止 62,500・窓 20,000・間隔 2,500 の 2 窓 → 17 出力 (22,500〜62,500)",
          SR.required_window_steps(62500, 20000, 2500) == req17 and len(req17) == 17)
    c = HJ.classify_ab([{"step": 62500, "D_mm": 0.0, "pf_A": True, "pf_B": True}], 0.0, req17)
    check(f"反例 2 (停止時の 1 行だけ): 必要な 17 出力のうち 16 を比べていない → 判定不能 ({c['why'][:60]}…)", c["result"] == "判定不能"
          and c["code"] == "U_incomplete")
    per17 = [{"step": k, "D_mm": 0.005, "pf_A": True, "pf_B": True} for k in req17]
    ps = {k: {"alts": [{"res": "res_45000.h5", "pf_F": True, "max_mm_F": (0.030 if k == 42500 else 0.0)}]} for k in req17}
    sg3 = {"status": "ok", "mm": max(a["max_mm_F"] for v in ps.values() for a in v["alts"]), "per_step": ps}
    c = HJ.classify_ab(per17, sg3, req17)
    rf = HJ.reference_final({"classify": c, "sigma_ref": sg3}, True)
    check(f"反例 3 (17 出力で D 0.005 mm、停止の出力の感度 0・途中の 1 出力だけ代わりの参照で 0.030 mm): σ_ref = 全出力の最大 0.030 mm → "
          f"判定不能 (参照の未整定)、参照の最終の資格なし ({c['code']})",
          sg3["mm"] == 0.030 and c["result"] == "判定不能" and c["code"] == "U_sigma_straddle" and rf["qualified"] is False)
    # 反例 1 (中央 1 列の欠測): 比べる列 F を先に固定 → 腕 A の F の中の欠測は合否の違い (結果 A にしない)
    x3 = np.array([10.0, 11.0, 12.0])
    F3 = np.ones(3, dtype=bool)
    dB = {"delta_E": [1e-3] * 3, "delta_target": [2e-3] * 3, "held": [0] * 3, "missing": []}
    dA = {"delta_E": [1e-3, float("nan"), 1e-3], "delta_target": [2e-3, float("nan"), 2e-3], "held": [0] * 3, "missing": [1]}
    r1 = HJ.ab_pair(62500, dA, dB, 76.6, "記録", x3, F3)
    check(f"反例 1 (中央 1 列だけ腕 A が欠測): ab_pair は比べる列 3 を保ち、腕 A は合否 否 (欠測 1・x {r1.get('x_bad_A')})・D は出さない",
          r1["n_fixed"] == 3 and r1["pf_A"] is False and r1["pf_B"] is True and r1["D_mm"] is None and r1["check_A"]["n_missing"] == 1
          and r1.get("x_bad_A") == [11.0], str(r1))
    c1 = HJ.classify_ab([r1], {"status": "ok", "mm": 0.0, "per_step": {62500: {"alts": [{"res": "r", "pf_F": True}]}}}, [62500])
    c1b = HJ.classify_ab([r1], 0.0, [62500])
    check(f"反例 1: 代わりの参照でも腕 B が合 → 結果 B (合否の違い) / 代わりの参照で確かめられなければ判定不能 — どちらも結果 A ではない "
          f"({c1['code']}・{c1b['code']})", c1["result"] == "結果 B" and c1["code"] == "B_pf" and c1b["result"] == "判定不能")
    r_both = HJ.ab_pair(62500, dA, dA, 76.6, "記録", x3, F3)
    check("両腕とも同じ列が欠測 → 合否は同じ (否)・D を出さない → 判定不能 (F の中の欠測)", r_both["D_mm"] is None
          and HJ.classify_ab([r_both], 0.0, [62500])["code"] == "U_fixed_cols")
    r_out = HJ.ab_pair(62500, dA, dB, 76.6, "記録", x3, np.array([True, False, True]))
    check("欠測の列が座標で先に決めた比べる列の外 (範囲の外) なら比べない: 合否は同じ・D 0", r_out["pf_A"] and r_out["D_mm"] == 0.0)
    # 参照の最終の資格 (reference_final; §5.1 #36i M2・M3)
    mk = lambda code, use=0.001, st="ok": {"classify": {"code": code, "result": "x", "why": "試験"},  # noqa: E731
                                           "sigma_ref": {"mm": 0.001, "mm_use": use, "status": st}}
    check("参照の最終の資格: 結果 A・B (σ_ref,use ≤ 0.01 mm) → 資格あり / σ を測れない・判定不能・参照の未整定 → 資格なし / 時間の変動の判定不能 "
          "(U_time、§5.1 #36i M2 で「あり」から外した)・出力の欠け・抽出の失敗・F の中の欠測 → 確かめられない (None) / --no-ab → 確かめられない",
          all(HJ.reference_final(mk(k), True)["qualified"] is True for k in ("A", "B_D", "B_pf"))
          and all(HJ.reference_final(mk(k), True)["qualified"] is False for k in ("U_sigma_unmeasurable", "U_sigma_undecidable", "U_sigma_straddle"))
          and all(HJ.reference_final(mk(k), True)["qualified"] is None for k in ("U_time", "U_incomplete", "U_error", "U_fixed_cols", "U_pf_unconfirmed"))
          and HJ.reference_final(None, False)["qualified"] is None and "U_time" not in HJ.REF_OK)
    check("§5.1 #36i M3: σ_ref,use > 0.01 mm は A/B の分類 (結果 A・B) に関係なく資格なし、ちょうど 0.01 mm は資格あり、σ_ref,use を測れない "
          "(None・status が ok でない) → 資格なし",
          all(HJ.reference_final(mk(k, use=0.0101), True)["qualified"] is False for k in ("A", "B_D", "B_pf", "U_time"))
          and HJ.reference_final(mk("A", use=0.01), True)["qualified"] is True
          and HJ.reference_final(mk("A", use=None), True)["qualified"] is False
          and HJ.reference_final(mk("A", st="undecidable"), True)["qualified"] is False)
    # 投入前のレビューの諮問の反例 (2026-10-14-h370-prelaunch-review-diagnose.md) の分類・資格の単体の版 (全経路の版は J)
    for name, Ds, sg_, use, want_code in (
            ("反例 1: 17 出力で D 0.005、全代替 20 出力で σ_ref 0.030 (NS 82,500 × Euler 42,000)", [0.005] * 17, 0.030, 0.030, "U_sigma_straddle"),
            ("反例 2: D 0.008/0.012 が窓の中で混在・σ_ref 0.030", [0.008, 0.012] * 8 + [0.008], 0.030, 0.030, "U_time"),
            ("反例 2′ (M2 だけ): D 0.008/0.012・σ_ref 0", [0.008, 0.012] * 8 + [0.008], 0.0, 0.0, "U_time"),
            ("反例 3: 比べる列 F の σ_ref 0、F の外の使用域だけ σ_ref,use 0.030", [0.005] * 17, 0.0, 0.030, "A"),
            ("反例 4: 全出力で D 0.100・σ_ref 0.030", [0.100] * 17, 0.030, 0.030, "B_D")):
        c = HJ.classify_ab(per17[:0] + [{"step": k, "D_mm": d, "pf_A": True, "pf_B": True} for k, d in zip(req17, Ds)],
                           {"status": "ok", "mm": sg_, "per_step": None}, req17)
        rf = HJ.reference_final({"classify": c, "sigma_ref": {"status": "ok", "mm": sg_, "mm_use": use}}, True)
        check(f"{name} → 分類 {c['code']}・資格 {rf['qualified']} (更新不可)", c["code"] == want_code and rf["qualified"] is not True, rf["why"][:120])


# ====================================================================================== F. NS の判定
NS_S, NI, NJ, X0, X1 = 0.0766, 56, 48, -2.0, 20.0
S_STOP, W_TEST = 62500, 5000


def _rwE(x):
    return 1.0 + 0.002 * x * x


def _dl(x):
    return 0.01 + 0.002 * (x + 2.0)


def _core(x, r, rw):
    return 100.0 * (1.0 + 0.02 * np.tanh(x / 10.0)) * (1.0 - 0.01 * (r / rw) ** 2)


def _radial(n):
    return 1.0 - (1.0 - np.linspace(0.0, 1.0, n)) ** 2


def write_field_run(rd: Path, X, R, Q, info: dict, steps, perturb=None):
    rd.mkdir(parents=True, exist_ok=True)
    ni, nj = X.shape
    with h5py.File(rd / "nozzle.h5", "w") as f:
        f["/MESH/COORD"] = np.c_[X.ravel() * NS_S, R.ravel() * NS_S, np.zeros(ni * nj)].ravel()
    for st in steps:
        k = 1.0 if perturb is None else perturb(st)
        with h5py.File(rd / f"res_{st}.h5", "w") as f:
            f["/VALUE/ro"] = np.ones(ni * nj)
            f["/VALUE/Ux"] = (Q * k).ravel()
            f["/VALUE/Uy"] = np.zeros(ni * nj)
            f["/VALUE/sonic"] = np.full(ni * nj, 300.0)
    (rd / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": ni, "nj": nj}, "scale_m": NS_S, "x_E": 15.0, **info}))


def synth_ns():
    x = np.linspace(X0, X1, NI)
    rw = _rwE(x) + _dl(x)
    R = rw[:, None] * _radial(NJ)[None, :]
    d = rw[:, None] - R
    slope = (0.004 * x + 0.002) / np.sqrt(1.0 + (0.004 * x + 0.002) ** 2)
    beta = np.clip(1.0 - d / (0.1 * rw[:, None]), 0.0, 1.0)
    X = x[:, None] + beta * d * slope[:, None]
    y = rw[:, None] - R
    return X, R, _core(X, R, rw[:, None]) * np.minimum(1.0, np.maximum(y / (5.0 * _dl(X)), 0.0) ** (1.0 / 7.0))


def synth_eu(x1: float, shape: float = 0.0):
    n = int(round(69 * (x1 - X0 + 0.5) / (X1 - X0 + 1.0))) + 1
    x = np.linspace(X0 - 0.5, x1, n)
    rw = _rwE(x)
    R = rw[:, None] * _radial(NJ)[None, :]
    X = np.repeat(x[:, None], NJ, axis=1)
    return X, R, _core(X, R, rw[:, None]) * (1.0 + shape * (R / rw[:, None]) ** 8)


def write_eu_shape(rd: Path, step: int, shape: float):
    X, R, Q = synth_eu(X1 + 0.5, shape)
    with h5py.File(rd / f"res_{step}.h5", "w") as f:
        f["/VALUE/ro"] = np.ones(Q.size)
        f["/VALUE/Ux"] = Q.ravel()
        f["/VALUE/Uy"] = np.zeros(Q.size)
        f["/VALUE/sonic"] = np.full(Q.size, 300.0)


ALT_SHAPES = {k: 1e-4 * (60000 - k) / 20000 for k in range(40000, 60000, 1000)}   # σ_ref の代わりの参照 (末尾の窓の全 20 出力) の形の小さな違い
SCREEN_SHAPES = {37000: 1e-4}                                           # スクリーニングの出力 (37,000; もう 1 つの 49,000 は代わりの参照と同じ出力)
REPS4 = [40000, 45000, 50000, 55000]                                    # 旧い代表 (試験の時間のため Euler の出力の間隔を 5,000 にした判定で使う)


def set_alts(case: Path, steps):
    """予備の資格の記録の σ_ref の代わりの参照を steps の出力 (今の sha) に書き換える。"""
    eb = case / HJ.EU_RUN
    q = json.loads((case / HJ.QUAL).read_text())
    q["sigma_ref_snapshots"] = [{"res": f"res_{k}.h5", "step": k, "sha256": HJ.sha256_file(eb / f"res_{k}.h5")} for k in steps]
    (case / HJ.QUAL).write_text(json.dumps(q))


def alts4():
    """試験の時間のため Euler の出力の間隔を 5,000 にする (代わりの参照 = 旧い代表 4 つ; 登録の 20 出力は A・D・F の結果 A・J・K で確かめる)。"""
    return Patch((HJ, "EU_OUT", 5000))


def make_ns_case(case: Path):
    """合成の NS の判定の case dir: Euler B (run_0950、固定 res_60000、σ_ref の代わりの参照 res_40000〜59000 の 20 出力、スクリーニング
    res_37000・res_49000)・腕 A (run_0174、範囲 x ≤ 19 の同じ形)・NS (run_0951、出力 52,500〜62,500、停止の出力 62,500、窓 5,000 の 2 窓)・
    記録一式。"""
    case.mkdir(parents=True)
    eb = case / HJ.EU_RUN
    X, R, Q = synth_eu(X1 + 0.5)
    write_field_run(eb, X, R, Q, {}, (60000,))
    for st, sh in {**ALT_SHAPES, **SCREEN_SHAPES}.items():
        write_eu_shape(eb, st, sh)
    pin = {"run": HJ.EU_RUN, "res": "res_60000.h5", "sha256": HJ.sha256_file(eb / "res_60000.h5"), "step": 60000}
    (case / "_band_ab/half_inlet").mkdir(parents=True)
    (case / HJ.PIN).write_text(json.dumps(pin))
    snap = lambda steps: [{"res": f"res_{s_}.h5", "step": s_, "sha256": HJ.sha256_file(eb / f"res_{s_}.h5")} for s_ in steps]  # noqa: E731
    qual = {"prequalified": True, "pin": pin, "sigma_ref_snapshots": snap(sorted(ALT_SHAPES)), "screening_snapshots": snap([37000, 49000]),
            "sigma_ref_selection": HJ.sigma_ref_selection(sorted([*ALT_SHAPES, 60000]))}
    (case / HJ.QUAL).write_text(json.dumps(qual))
    ea = case / HJ.OLD_EULER
    X, R, Q = synth_eu(19.0)
    write_field_run(ea, X, R, Q, {}, (12000,))
    spec_b = {"run": str(eb), "res": "res_60000.h5", "sha256": pin["sha256"]}
    ns = case / HJ.NS_RUN
    sr = SR.resolve_settings({"stop_rule": True, "stop_rule_thresholds": {"window": W_TEST, "min_points": 2}}, nsteps_arg=None, n_main=200000,
                             n_main_source="試験", out_int=2500, euler_ref=spec_b)
    steps = list(range(S_STOP - 2 * W_TEST, S_STOP + 1, 2500))
    X, R, Q = synth_ns()
    write_field_run(ns, X, R, Q, {"wall_thermal": {"mode": "isothermal", "Tw": 370.0}, "physical_wall": {"domain": [X0, X1]},
                                  "euler_ref": spec_b, "stop_rule": sr}, steps, perturb=lambda st: 1.0 + 1e-6 * (st - S_STOP) / 2500.0)
    uc = SR.use_columns(ns)
    ext = {k: v for k, v in sr["extract"].items() if k != "source"}
    rows = []
    for st in steps:
        row = {"step": st, "nonfinite_all": 0, "r": 1e-5 * np.sin(st), "exit_core_M": 5.99 + 1e-6 * np.cos(st), "mdot_rt": 10.0,
               SR.theta_key(0.5): 0.01}
        if st != steps[1]:                                      # 1 出力は記録なし (判定が抽出し直す)
            row["delta"] = SR.delta_snapshot(ns, st, spec_b, ext, uc["x"])
        rows.append(row)
    rec = {"version": 1, "settings": {**sr, "x_fracs": [0.5]}, "status": "REACHED", "rows": rows, "stop_step": S_STOP,
           "sections": [{"key": SR.theta_key(0.5), "x_frac": 0.5, "x_requested": 10.0}],
           "delta_series": {"x": uc["x"].tolist(), "use_domain": [0.0, uc["x_w"]], "scale_m": uc["scale_m"], "extract": sr["extract"],
                            "euler_ref": sr["euler_ref"], "span_steps": 2 * W_TEST}}
    SR.save_record(ns, rec)
    SR.write_state(ns, {"status": "REACHED", "stop_step": S_STOP, "forge_rc": 143, "stop_sha256": HJ.sha256_file(ns / f"res_{S_STOP}.h5")})
    (ns / f"res_wall_3_{S_STOP}.h5").write_bytes(b"")
    (ns / "H370_PRECHECK.json").write_text(json.dumps({"ok": True}))
    (ns / "H370_MAIN_START.json").write_text(json.dumps({"ok": True}))
    (ns / "RUN_RC").write_text("0\n")
    (ns / "RUN_PROVENANCE.txt").write_text("forge_sha256: 7e95e5bbd8f50c9e0000\nenv         : FORGE_CUDA_BLOCKSIZE=128\n")
    for t in ("S1_soft", "S2_mid", "main"):
        (ns / f"residual_history_{t}.csv").write_text("step,rms_ro\n1,1e-3\n2,1e-4\n")
    (case / HJ.NS_PROBLEM).write_text("# 試験の偽物 (関門・次の壁は偽物で切り替える)\n")
    return case


def stub_gate(ok=True):
    def g(run, **kw):
        return {"proceed": ok, "blocked_by": [] if ok else ["収束の判定 NOT_CONVERGED (試験の偽物)"],
                "stop_rule": {"status": "REACHED"}, "convergence": {"kind": "PASS" if ok else "NOT_CONVERGED"},
                "checks": {"two_windows": {"pass": ok}, "wall_resolution": {"pass": ok, "kind": "PASS"}}, "extract": {"made": []},
                "override": {"allow_not_reached": False}, "_kw": {k: (str(v) if isinstance(v, Path) else v) for k, v in kw.items()}}
    return g


def stub_next_wall(problem, table):
    a = np.loadtxt(table, delimiter=",", skiprows=1, ndmin=2)
    return {"pass": True, "throat": {"x": 0.0, "r": 1.0 + float(np.interp(0.0, a[:, 0], a[:, 1])), "kappa": 1.0}}


def judge(case: Path, extra=()) -> tuple:
    rc = HJ.main(["--case-dir", str(case), "--design", str(ROOT / "design"), "--tools", str(TOOLS), *extra])
    return rc, HJ.jload(case / HJ.OUT)


def snap_hook(fn):
    """本物の stop_rule.delta_snapshot の結果を fn(run, step, spec, d) で書き換える試験の差し替え (反例を組み込みで作る)。"""
    def g(run, step, spec, ext, x_use):
        return fn(run, step, spec, SR.delta_snapshot(run, step, spec, ext, x_use))
    return g


def mid_col(R):
    """比べる列 F の中ほどの列の番号 (delta_series.x の番号)。"""
    x = np.asarray(json.loads((R["_case"] / HJ.NS_RUN / SR.RECORD).read_text())["delta_series"]["x"])
    lo, hi = R["ab"]["columns"]["x_range"]
    return int(np.argmin(np.abs(x - 0.5 * (lo + hi))))


TIME_KEYS = ("seconds", "sigma_ref_seconds", "parallel", "created", "children_maxrss_mib")


def norm_record(d, case: Path):
    """判定の記録を比べる形に: case dir のパスを置き換え、時間・並列の記録 (TIME_KEYS) を除く (値はそのまま; 浮動小数を丸めない)。"""
    def walk(o):
        if isinstance(o, dict):
            return {k: walk(v) for k, v in o.items() if k not in TIME_KEYS}
        if isinstance(o, list):
            return [walk(v) for v in o]
        if isinstance(o, str):
            return o.replace(str(case.resolve()), "<case>").replace(str(case), "<case>")
        return o
    return json.loads(json.dumps(walk(d), default=str))


def first_diff(a, b, path=""):
    """2 つの記録の最初の違い (試験の表示用)。"""
    if type(a) is not type(b):
        return f"{path}: 型 {type(a).__name__} ≠ {type(b).__name__}"
    if isinstance(a, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                return f"{path}/{k}: 片方だけ"
            r = first_diff(a[k], b[k], f"{path}/{k}")
            if r:
                return r
        return ""
    if isinstance(a, list):
        if len(a) != len(b):
            return f"{path}: 長さ {len(a)} ≠ {len(b)}"
        for i, (x, y) in enumerate(zip(a, b)):
            r = first_diff(x, y, f"{path}[{i}]")
            if r:
                return r
        return ""
    return "" if a == b else f"{path}: {a!r} ≠ {b!r}"


def test_ns_judge():
    print("== F. NS の判定 (合成の NS・Euler の run、本物の抽出、関門と次の壁は偽物)", flush=True)
    base = make_ns_case(TMP / "ns_base")
    hooks = {"handoff_gate": stub_gate(True), "next_wall_check": stub_next_wall}
    req = SR.required_window_steps(S_STOP, W_TEST, 2500)
    with Patch((HJ, "HOOKS", dict(hooks))):
        case = TMP / "ns_A"
        shutil.copytree(base, case)
        rc, R = judge(case)
        R["_case"] = case
        ab = R["ab"]
        srcs = {r.get("B_source", "")[:2] for r in ab["per_output"]}
        sg = ab["sigma_ref"]
        check(f"結果 A: 腕 A = 範囲の短い同じ参照 (D max {ab.get('D_mm', {}).get('max', float('nan')):.4f} mm、σ_ref {sg['mm']:.2e} mm "
              f"(使用域 {sg['mm_use']:.2e})、{ab['classify']['code']})", rc == 0 and ab["result"] == "結果 A" and ab["D_mm"]["max"] < HJ.TOL_MM
              and sg["status"] == "ok", ab["why"])
        cols = ab["columns"]
        check(f"比べる列は抽出の前に座標と範囲から固定: 使用域 {cols['n_use']} 列のうち F {cols['n']} 列 (腕 A の範囲の外 {cols['excluded']['A']['n']} 列 "
              f"x {cols['excluded']['A']['x_range']} は比べない)", cols["excluded"]["A"]["n"] == 3 and cols["excluded"]["A"]["x_range"][0] > 19.0
              and cols["excluded"]["B"]["n"] == 0 and cols["n"] == cols["n_use"] - 3 and cols["x_range"][1] <= 19.0 and cols["edge_tol"])
        check(f"判定窓の完全性: 入るべき出力 {ab['completeness']['n_required']} (出力の間隔 2,500 と 2 窓から列挙) をすべて比べた; 腕 B の値は記録 "
              f"(sha が同じ) と、記録の無い出力は抽出し直し", ab["completeness"]["required"] == req and [r["step"] for r in ab["per_output"]] == req
              and srcs == {"記録", "抽出"} and ab["completeness"]["reextracted_B"] == [req[1]], str(srcs))
        check(f"σ_ref: 判定窓の全 {len(req)} 出力 × 末尾の窓の全 20 出力 (40,000〜59,000) で抽出し直した (合否・列の印の変化なし; "
              f"抽出 {ab['extractions']['sigma_ref']} 回、{ab['extractions']['parallel']['mode']} {ab['extractions']['parallel']['jobs']})",
              sorted(int(k) for k in sg["per_step"]) == req and all(len(v["alts"]) == 20 and all(a_["pf_F"] and a_["pf_use"] and not a_["change"] for a_ in v["alts"])
                                                    for v in sg["per_step"].values()) and sg["alts"] == [f"res_{k}.h5" for k in sorted(ALT_SHAPES)]
              and ab["extractions"]["sigma_ref"] == 20 * len(req) and sg["alternatives_check"]["ok"])
        check(f"σ_ref,use (使用域の全体 {len(json.loads((case / HJ.NS_RUN / SR.RECORD).read_text())['delta_series']['x'])} 列) と σ_ref (F) の最大の位置 "
              f"(σ_ref {sg['mm']:.2e} at {sg['at_F']}、σ_ref,use {sg['mm_use']:.2e} at {sg['at_use']}) を記録・最終の資格の記録にも写す",
              sg["mm_use"] >= sg["mm"] and sg["at_use"]["step"] in req and sg["at_use"]["euler_res"] in sg["alts"]
              and R["reference_final"]["sigma_ref_use_mm"] == sg["mm_use"] and R["reference_final"]["sigma_ref_use_at"] == sg["at_use"]
              and R["reference_final"]["sigma_ref_at"] == sg["at_F"])
        check("スクリーニング (停止の出力 × 37,000・49,000) は記録だけ", ab["screening"]["label"].startswith("スクリーニング")
              and len(ab["screening"]["per_snapshot"]) == 2)
        rf = R["reference_final"]
        check(f"参照の最終の資格: あり ({rf['why'][:60]}…)", rf["qualified"] is True and rf["ab_code"] == "A" and rf["sigma_ref_use_mm"] <= HJ.TOL_MM)
        check(f"判定の全体: 「更新に進める条件をすべて満たす / 抽出の参照の A/B: 結果 A」({R['verdict'][:60]}…)",
              R["verdict"].startswith("更新に進める条件をすべて満たす") and "結果 A" in R["verdict"] and R["update_conditions_met"])
        st = ab["at_stop"]
        check("停止の出力の両腕の全体の抽出: 共通の列・帯の位置 (y_b・y*)・感度・held・δ_target・未緩和の固定点の残差・次の壁の真のスロート",
              st["compare"]["n_common_use"] > 0 and "band_y_b_max_diff_mm" in st["compare"] and "band_sensitivity_rel" in st["compare"]
              and "B_fixed_point_residual_mm" in st["compare"] and "throat_diff" in st["next_wall"]
              and all((case / HJ.NS_RUN / HJ.AB_DIR / d / "delta_r_next.csv").is_file() for d in ("A", "B", "A_on_common")), str(st.get("compare"))[:300])
        q = R.get("quasisteady") or {}
        check("check_quasisteady の記録 (2 窓の r・θ_r・出口 M・流量・δ_E・δ_target)", q.get("n") == 5 and "dE_40" in q["default"]["verdicts"]
              and "exit_core_M" in q["default"]["verdicts"], str(q.get("default", {}).get("missing")))
        check("関門は固定の Euler の参照・ω 0.5・上書きなしで呼ばれる",
              R["gate"]["proceed"] and R["stop"]["stop_sha256_match"] and R["stop"]["wall_dump_at_stop"] == [f"res_wall_3_{S_STOP}.h5"])
        j_mid = mid_col(R)
        mmpr = NS_S * 1e3
        # 直列 (--jobs 1) と並列 (既定) で判定の記録 (A/B・σ_ref・σ_ref,use・位置・分類・資格) がビット同一 (§5.1 #36i 計算量)
        R_par, case_par = R, case
        case = TMP / "ns_A_serial"
        shutil.copytree(base, case)
        rc, R1 = judge(case, ["--jobs", "1"])
        pp, ps = R_par["ab"]["extractions"]["parallel"], R1["ab"]["extractions"]["parallel"]
        a_par, a_ser = norm_record(R_par["ab"], case_par), norm_record(R1["ab"], case)
        check(f"並列 ({pp['mode']} {pp['jobs']}、{pp['seconds']} s) と直列 ({ps['mode']} {ps['jobs']}、{ps['seconds']} s) で A/B・σ_ref・σ_ref,use・"
              f"最大の位置・分類の記録がビット同一、最終の資格・更新の条件も同じ",
              pp["jobs"] > 1 and pp["mode"].startswith("process") and ps["mode"].startswith("process") and ps["jobs"] == 1 and a_par == a_ser
              and norm_record(R_par["reference_final"], case_par) == norm_record(R1["reference_final"], case)
              and R_par["update_conditions_met"] == R1["update_conditions_met"],
              "" if a_par == a_ser else first_diff(a_par, a_ser))
        check(f"並列数の記録: min(上限 3, 空きメモリ {R_par['parallel_plan']['mem_available_mib']} MiB / {HJ.SNAP_PEAK_MIB} MiB) = "
              f"{R_par['parallel_plan']['jobs']}・子の最大メモリ {pp['children_maxrss_mib']} MiB", R_par["parallel_plan"]["cap"] == HJ.JOBS_CAP
              and R_par["parallel_plan"]["jobs"] == pp["jobs"] and R1["parallel_plan"]["cap"] == 1)
        # 結果 B: 腕 A の参照の形を変える (以下、試験の時間のため代わりの参照を旧い代表 4 つにする: alts4・set_alts)
        case = TMP / "ns_B"
        shutil.copytree(base, case)
        set_alts(case, REPS4)
        shutil.rmtree(case / HJ.OLD_EULER)
        X, Rr, Q = synth_eu(19.0, 0.01)
        write_field_run(case / HJ.OLD_EULER, X, Rr, Q, {}, (12000,))
        with alts4():
            rc, R = judge(case)
        check(f"結果 B: 腕 A の参照の形を変える (D min {R['ab'].get('D_mm', {}).get('min', float('nan')):.4f} mm > 0.01 mm)",
              rc == 0 and R["ab"]["result"] == "結果 B" and R["ab"]["classify"]["code"] == "B_D" and R["reference_final"]["qualified"] is True,
              R["ab"]["why"])
        # 判定不能: 代わりの参照 (45,000) の形が大きく違う → σ_ref が大きい → 参照の最終の資格なし → 更新に進めない
        case = TMP / "ns_sigma"
        shutil.copytree(base, case)
        eb = case / HJ.EU_RUN
        write_eu_shape(eb, 45000, 0.003)
        set_alts(case, REPS4)
        with alts4():
            rc, R = judge(case)
        check(f"判定不能: σ_ref {R['ab']['sigma_ref']['mm']:.4f} mm で 0.01 mm をまたぐ (参照の未整定) → 参照の最終の資格なし → "
              "関門が進めても更新に進める条件を満たさない",
              R["ab"]["result"] == "判定不能" and R["ab"]["classify"]["code"] == "U_sigma_straddle" and R["reference_final"]["qualified"] is False
              and R["gate"]["proceed"] and not R["update_conditions_met"] and R["verdict"].startswith("更新に進める条件を満たさない")
              and "参照の最終の資格" in R["verdict"], R["ab"]["why"])
        case = TMP / "ns_sigma_sha"
        shutil.copytree(base, case)
        q = json.loads((case / HJ.QUAL).read_text())
        q["sigma_ref_snapshots"][1]["sha256"] = "0" * 64
        (case / HJ.QUAL).write_text(json.dumps(q))
        rc, R = judge(case, ["--no-ab"])
        check("--no-ab: A/B は未実施・参照の最終の資格を確かめていない → 更新に進める条件を満たさない",
              R["ab"]["result"].startswith("未実施") and R["reference_final"]["qualified"] is None and not R["update_conditions_met"])
        rc, R = judge(case)
        ac = R["ab"]["sigma_ref"]["alternatives_check"]
        check(f"代わりの参照 (41,000) の sha が違う → 抽出の前の照合で止める (σ_ref の抽出 {R['ab']['extractions']['sigma_ref']} 回) → σ_ref を測れない → "
              "判定不能・参照の最終の資格なし", R["ab"]["result"] == "判定不能"
              and R["ab"]["sigma_ref"]["status"] == "unmeasurable" and R["ab"]["sigma_ref"]["mm"] is None and ac["sha_mismatch"] == [41000]
              and not ac["ok"] and R["ab"]["extractions"]["sigma_ref"] == 0
              and R["reference_final"]["qualified"] is False and not R["update_conditions_met"], str(ac.get("why")))
        # 行が欠けても res があれば抽出し直して補う (結果 A のまま)
        case = TMP / "ns_rows"
        shutil.copytree(base, case)
        rp = case / HJ.NS_RUN / SR.RECORD
        rec = json.loads(rp.read_text())
        rec["rows"] = [r for r in rec["rows"] if r["step"] == S_STOP]
        rp.write_text(json.dumps(rec))
        set_alts(case, REPS4)
        with alts4():
            rc, R = judge(case)
        check("stop_rule.json に停止の行だけ (res は全部ある) → 入るべき出力を列挙して抽出し直して補う → 結果 A",
              R["ab"]["result"] == "結果 A" and R["ab"]["completeness"]["rows_missing"] == req[:-1]
              and R["ab"]["completeness"]["reextracted_B"] == req[:-1], R["ab"]["why"])
        # 反例 2 (組み込み): 停止の行だけ・他の出力の res も無い → 判定不能
        for k in req[:-1]:
            (case / HJ.NS_RUN / f"res_{k}.h5").unlink()
        with alts4():
            rc, R = judge(case)
        check(f"反例 2 (組み込み): 停止の 1 行だけで他の res も無い (再抽出で補えない) → 判定不能 ({R['ab']['classify']['code']})",
              R["ab"]["result"] == "判定不能" and R["ab"]["classify"]["code"] == "U_incomplete" and R["ab"]["completeness"]["res_missing"] == req[:-1]
              and not R["update_conditions_met"], R["ab"]["why"])
    # 反例 1 (組み込み): 腕 A の中ほどの 1 列だけ欠測 (F の中) → 結果 A にしない
    def hole(run, step, spec, d):
        if str(Path(spec["run"]).resolve()).endswith(HJ.OLD_EULER):
            d = {**d, "delta_E": list(d["delta_E"]), "missing": list(d["missing"]) + [j_mid]}
            d["delta_E"][j_mid] = float("nan")
        return d
    with Patch((HJ, "HOOKS", {**hooks, "delta_snapshot": snap_hook(hole)})), alts4():
        case = TMP / "ns_cx1"
        shutil.copytree(base, case)
        set_alts(case, REPS4)
        rc, R = judge(case)
        per = R["ab"]["per_output"]
        check(f"反例 1 (組み込み): 腕 A の比べる列の中ほど (x {per[0].get('x_bad_A')}) だけ欠測 → 比べる列は減らさず腕 A は合否 否 → "
              f"代わりの参照の全部でも腕 B は合 → 結果 B ({R['ab']['classify']['code']}、結果 A ではない)",
              R["ab"]["result"] == "結果 B" and R["ab"]["classify"]["code"] == "B_pf" and all(r["pf_A"] is False and r["n_fixed"] == R["ab"]["columns"]["n"]
                                                                                             for r in per), R["ab"]["why"][:200])
    # 反例 3 (組み込み): 途中の 1 出力だけ代わりの参照 (50,000) で δ_E が 0.030 mm 動く → 判定不能・参照の最終の資格なし
    k_mid = req[2]

    def bump(run, step, spec, d):
        if spec.get("res") == "res_50000.h5" and step == k_mid:
            d = {**d, "delta_E": list(d["delta_E"])}
            d["delta_E"][j_mid] += 0.030 / mmpr
        return d
    with Patch((HJ, "HOOKS", {**hooks, "delta_snapshot": snap_hook(bump)})), alts4():
        case = TMP / "ns_cx3"
        shutil.copytree(base, case)
        set_alts(case, REPS4)
        rc, R = judge(case)
        sg = R["ab"]["sigma_ref"]
        at = [a_["max_mm_F"] for a_ in sg["per_step"][str(k_mid)]["alts"] if a_["res"] == "res_50000.h5"]
        check(f"反例 3 (組み込み): 停止の出力 {S_STOP} の感度は小さく、途中の出力 {k_mid} × res_50000 だけ 0.030 mm → σ_ref {sg['mm']:.4f} mm → "
              f"判定不能 (参照の未整定)・参照の最終の資格なし・更新に進めない",
              0.029 < sg["mm"] < 0.031 and at and abs(at[0] - sg["mm"]) < 1e-12 and R["ab"]["classify"]["code"] == "U_sigma_straddle"
              and R["reference_final"]["qualified"] is False and not R["update_conditions_met"], R["ab"]["why"][:200])
    # 比べる列の固定に出口の幾何の許容 (Euler の端の外側だけ、端の列の間隔 × 1e−4) が効く
    xu = np.asarray(json.loads((base / HJ.NS_RUN / SR.RECORD).read_text())["delta_series"]["x"])
    for name, short, want in (("1.68e−9 r_t 短い (変換の丸め)", 1.68e-9, True), ("0.07 r_t 短い (範囲の不足)", 0.07, False)):
        rd = TMP / f"fx_{name[:4]}"
        X, Rr, Q = synth_eu(X1)
        X = X.copy()
        X[-1, :] = xu[-1] - short
        write_field_run(rd, X, Rr, Q, {}, (60000,))
        fx = HJ.fixed_columns(xu, {"B": rd}, edge_tol=True)
        fx0 = HJ.fixed_columns(xu, {"B": rd}, edge_tol=False)
        check(f"比べる列の固定: Euler の端が NS の最後の列より {name} → 最後の列は {'F に入る' if want else '入らない (欠測のまま)'}"
              f" (許容 {fx['cover']['B']['tol_hi']:.2e} r_t; 許容なし (±1e−9) では入らない)",
              bool(fx["mask"][-1]) == want and not fx0["mask"][-1] and fx["mask"][:-2].all())
    with Patch((HJ, "HOOKS", {"handoff_gate": stub_gate(False), "next_wall_check": stub_next_wall})):
        case = TMP / "ns_blocked"
        shutil.copytree(base, case)
        rc, R = judge(case, ["--no-ab"])
        check("関門で止まる → 「更新に進める条件を満たさない: 理由」", rc == 0 and R["verdict"].startswith("更新に進める条件を満たさない")
              and "NOT_CONVERGED" in R["verdict"] and not R["update_conditions_met"])
        case = TMP / "ns_nomain"
        shutil.copytree(base, case)
        (case / HJ.NS_RUN / "H370_MAIN_START.json").unlink()
        (case / HJ.NS_RUN / "DISK_STOP").write_text("試験")
        rc, R = judge(case, ["--no-ab"])
        check("本段の直前の照合の記録が無い・DISK_STOP → 判定不能 (rc 2)", rc == 2 and R["verdict"].startswith("判定不能")
              and any("H370_MAIN_START" in u for u in R["undecidable"]) and any("DISK_STOP" in u for u in R["undecidable"]))
        case = TMP / "ns_nowall"
        shutil.copytree(base, case)
        (case / HJ.NS_RUN / f"res_wall_3_{S_STOP}.h5").unlink()
        rc, R = judge(case, ["--no-ab"])
        check("停止の出力の壁の出力 (res_wall_3_<s>.h5) が無い → 判定不能 (局所 y₁⁺ を判定できない)", rc == 2
              and any("res_wall_3_" in u for u in R["undecidable"]))
        case = TMP / "ns_prov"
        shutil.copytree(base, case)
        (case / HJ.NS_RUN / "RUN_PROVENANCE.txt").write_text("forge_sha256: deadbeef\n")
        rc, R = judge(case, ["--no-ab"])
        check("RUN_PROVENANCE.txt のバイナリ・環境変数が違う → 判定不能", rc == 2 and any("RUN_PROVENANCE" in u for u in R["undecidable"]))


# ====================================================================================== J. 投入前のレビューの諮問の反例 (全経路、偽物の抽出)
MMPR = NS_S * 1e3                                                       # 1 r_t あたりの mm (合成の NS の scale_m)
S_FAST, W_FAST = 102500, 20000                                          # 諮問の反例の NS: 停止 102,500・窓 20,000 → 17 出力 (62,500〜102,500)


def make_fast_case(base: Path, case: Path):
    """反例の case (抽出は偽物 fake_snap なので速い): base を写し、判定窓を 2 窓 × 20,000 (17 出力、停止 102,500; 諮問の NS 82,500 を含む)
    にする。NS の res は base の res_62500 の写し (偽物の抽出は中身を読まない; sha は判定が照合)。行は δ を持たない (腕 B も偽物で抽出)。"""
    shutil.copytree(base, case)
    ns = case / HJ.NS_RUN
    req = SR.required_window_steps(S_FAST, W_FAST, 2500)
    for k in req:
        if not (ns / f"res_{k}.h5").exists():
            shutil.copy(ns / f"res_{S_STOP}.h5", ns / f"res_{k}.h5")
    rp = ns / SR.RECORD
    rec = json.loads(rp.read_text())
    rec["settings"]["window"] = W_FAST
    rec["stop_step"] = S_FAST
    rec["rows"] = [{"step": k, "nonfinite_all": 0, "r": 1e-5, "exit_core_M": 5.99, "mdot_rt": 10.0, SR.theta_key(0.5): 0.01} for k in req]
    rp.write_text(json.dumps(rec))
    SR.write_state(ns, {"status": "REACHED", "stop_step": S_FAST, "forge_rc": 143, "stop_sha256": HJ.sha256_file(ns / f"res_{S_FAST}.h5")})
    (ns / f"res_wall_3_{S_FAST}.h5").write_bytes(b"")
    return req


def fake_snap(effect):
    """偽物の抽出: δ_E = δ_target = 1e-3 r_t + effect(腕, NS の出力, Euler の出力, x) [mm] / MMPR。腕 A (旧い参照、x ≤ 19 の範囲) は範囲の外を欠測。"""
    def g(run, step, spec, ext, x_use):
        x = np.asarray(x_use, dtype=float)
        arm = "A" if str(spec["run"]).rstrip("/").endswith(HJ.OLD_EULER) else "B"
        e = 1e-3 + np.asarray(effect(arm, int(step), spec.get("res"), x), dtype=float) / MMPR
        if arm == "A":
            e = np.where(x > 19.0, np.nan, e)
        miss = np.flatnonzero(~np.isfinite(e))
        return {"delta_E": e.tolist(), "delta_target": e.tolist(), "held": [0] * len(x), "missing": miss.tolist(), "n_missing": int(len(miss)),
                "n_held": 0, "lam": 1.0, "res_sha256": None, "seconds": 0.0}
    return g


def bump_at(x, j, mm):
    v = np.zeros(len(x))
    v[j] = mm
    return v


def test_prelaunch_counterexamples():
    print("== J. 投入前のレビューの諮問の反例 4 つと判別 A/B (判定器の全経路、偽物の抽出)", flush=True)
    base = TMP / "ns_base"
    if not base.exists():
        make_ns_case(base)
    xu = np.asarray(json.loads((base / HJ.NS_RUN / SR.RECORD).read_text())["delta_series"]["x"])
    j_in = int(np.argmin(np.abs(xu - 9.0)))                              # 比べる列 F の中 (両腕の範囲の中)
    j_out = int(np.flatnonzero(xu > 19.0)[0])                            # 使用域の中・F の外 (腕 A の範囲の外)
    hooks = {"handoff_gate": stub_gate(True), "next_wall_check": stub_next_wall}

    def run_case(name, effect, alts=None):
        case = TMP / f"cx_{name}"
        req = make_fast_case(base, case)
        if alts is not None:
            set_alts(case, alts)
        with Patch((HJ, "HOOKS", {**hooks, "delta_snapshot": fake_snap(effect)})):
            rc, R = judge(case)
        return req, R

    def D(dmm):                                                            # 腕 A を腕 B から dmm [mm] ずらす (D_k = dmm)
        return lambda arm, k, res, x: -dmm * np.ones(len(x)) if arm == "A" else np.zeros(len(x))

    def plus(f, g):
        return lambda arm, k, res, x: f(arm, k, res, x) + g(arm, k, res, x)

    def alt_bump(k0, res0, j, mm=0.030):
        return lambda arm, k, res, x: bump_at(x, j, mm) if (arm == "B" and k == k0 and res == res0) else np.zeros(len(x))

    def show(R):
        ab, rf = R["ab"], R["reference_final"]
        sg = ab.get("sigma_ref") or {}
        return (f"A/B {ab.get('classify', {}).get('code')}・σ_ref {sg.get('mm')} at {sg.get('at_F')}・σ_ref,use {sg.get('mm_use')} at "
                f"{sg.get('at_use')}・資格 {rf.get('qualified')}・更新 {R.get('update_conditions_met')}")
    # 反例 1: 17 出力で D 0.005 mm、NS 82,500 × Euler 42,000 だけ 0.030 mm
    eff1 = plus(D(0.005), alt_bump(82500, "res_42000.h5", j_in))
    req, R = run_case("1_full", eff1)
    sg = R["ab"]["sigma_ref"]
    check(f"反例 1 (全 20 出力): NS 82,500 × Euler 42,000 だけ 0.030 mm → σ_ref 0.030 mm (位置も記録) → 判定不能 (参照の未整定)・資格なし・更新不可 "
          f"({show(R)})", len(req) == 17 and 82500 in req and abs(sg["mm"] - 0.030) < 1e-9 and sg["at_F"]["step"] == 82500
          and sg["at_F"]["euler_res"] == "res_42000.h5" and abs(sg["at_F"]["x"] - xu[j_in]) == 0.0 and R["ab"]["classify"]["code"] == "U_sigma_straddle"
          and R["reference_final"]["qualified"] is False and not R["update_conditions_met"] and R["ab"]["extractions"]["sigma_ref"] == 340)
    # 判別 A/B: 同じ入力で代わりの参照の選び方だけを 5,000 → 1,000 に
    req, Ra = run_case("1_reps", eff1, alts=REPS4)
    ac = Ra["ab"]["sigma_ref"]["alternatives_check"]
    # 旧い判定器の帰結 (代表 4 つで σ_ref を測ると取りこぼす): 同じ偽物の抽出で σ_ref を代表 4 つだけで測る
    case = TMP / "cx_1_reps"
    rec = json.loads((case / HJ.NS_RUN / SR.RECORD).read_text())
    F = np.ones(len(xu), dtype=bool)
    F[j_out:] = False
    snap = fake_snap(eff1)
    spec_b = {"run": str(case / HJ.EU_RUN), "res": "res_60000.h5", "sha256": "x"}
    dB = {k: snap(case / HJ.NS_RUN, k, spec_b, {}, xu) for k in req}
    alts4_ = [{"res": f"res_{k}.h5", "step": k, "sha256": "x"} for k in REPS4]
    s4 = HJ.sigma_ref_measure(case / HJ.NS_RUN, req, spec_b, alts4_, {}, xu, F, dB, MMPR, snap)
    per = [{"step": k, "D_mm": 0.005, "pf_A": True, "pf_B": True} for k in req]
    c4 = HJ.classify_ab(per, s4, req)
    old_ok = c4["code"] in ("A", "B_D", "B_pf", "U_time")                 # 旧い REF_OK (σ_ref,use の条件なし)
    check(f"判別 A/B (A = 間引き 5,000): 代表 4 つの σ_ref は {s4['mm']} mm で 42,000 の変動を取りこぼし、旧い判定器なら {c4['code']}・資格あり "
          f"({old_ok}) — 今の判定器は代わりの参照が全 20 出力でない (無い {len(ac['missing_in_record'])}) ので σ_ref の抽出を回さず測れない → "
          f"資格なし・更新不可 ({show(Ra)})", s4["mm"] == 0.0 and c4["code"] == "A" and old_ok and not ac["ok"]
          and len(ac["missing_in_record"]) == 16 and Ra["ab"]["sigma_ref"]["status"] == "unmeasurable" and Ra["ab"]["extractions"]["sigma_ref"] == 0
          and Ra["reference_final"]["qualified"] is False and not Ra["update_conditions_met"])
    check("判別 A/B (B = 全末尾 1,000): 上の反例 1 のとおり σ_ref 0.030 mm で資格なし — 両腕とも今の判定器では更新不可 (取りこぼし仮説は修正版で反証)",
          R["reference_final"]["qualified"] is False and not R["update_conditions_met"])
    # 反例 2: D 0.008/0.012 が窓の中で混在、σ_ref 0.030 mm
    dmix = lambda arm, k, res, x: (-(0.008 if (k // 2500) % 2 else 0.012) * np.ones(len(x))) if arm == "A" else np.zeros(len(x))  # noqa: E731
    req, R = run_case("2", plus(dmix, alt_bump(72500, "res_51000.h5", j_in)))
    Ds = [r["D_mm"] for r in R["ab"]["per_output"]]
    check(f"反例 2: D {min(Ds):.3f}/{max(Ds):.3f} mm が窓の中で混在・σ_ref 0.030 mm → 判定不能 (U_time) → 資格なし (σ_ref,use 0.030 > 0.01)・更新不可 "
          f"({show(R)})", abs(min(Ds) - 0.008) < 1e-9 and abs(max(Ds) - 0.012) < 1e-9 and R["ab"]["classify"]["code"] == "U_time"
          and R["reference_final"]["qualified"] is False and not R["update_conditions_met"])
    req, R = run_case("2b", dmix)
    check(f"反例 2′ (M2 だけ): 同じ D の混在で σ_ref 0 → U_time は資格を確かめられない (旧い判定器は資格あり) → 更新不可 ({show(R)})",
          R["ab"]["classify"]["code"] == "U_time" and R["ab"]["sigma_ref"]["mm_use"] == 0.0 and R["reference_final"]["qualified"] is None
          and not R["update_conditions_met"] and "判定不能" in R["reference_final"]["why"])
    # 反例 3: 比べる列 F の σ_ref は 0、F の外の使用域の列だけ σ_ref,use 0.030 mm
    req, R = run_case("3", plus(D(0.005), alt_bump(92500, "res_57000.h5", j_out)))
    sg = R["ab"]["sigma_ref"]
    check(f"反例 3: F の σ_ref 0・F の外 (x {xu[j_out]:.2f}) だけ σ_ref,use 0.030 mm → 結果 A でも資格なし・更新不可 ({show(R)})",
          R["ab"]["classify"]["code"] == "A" and sg["mm"] == 0.0 and abs(sg["mm_use"] - 0.030) < 1e-9 and sg["at_use"]["x"] == float(xu[j_out])
          and sg["at_use"]["step"] == 92500 and sg["at_use"]["euler_res"] == "res_57000.h5" and R["reference_final"]["qualified"] is False
          and "σ_ref,use" in R["reference_final"]["why"] and not R["update_conditions_met"])
    # 反例 4: 全出力で D 0.100 mm・σ_ref 0.030 mm
    req, R = run_case("4", plus(D(0.100), alt_bump(S_FAST, "res_59000.h5", j_in)))
    check(f"反例 4: 全出力で D 0.100 mm・σ_ref 0.030 mm → 結果 B (旧い参照を退ける) でも新しい参照は σ_ref,use 0.030 > 0.01 で資格なし・更新不可 "
          f"({show(R)})", R["ab"]["classify"]["code"] == "B_D" and abs(R["ab"]["sigma_ref"]["mm_use"] - 0.030) < 1e-9
          and R["reference_final"]["qualified"] is False and not R["update_conditions_met"])
    # 対照: 変動が無ければ (D 0.005・σ 0) 全経路で資格あり (反例の拒否が一律の拒否でないこと)
    req, R = run_case("ok", D(0.005))
    check(f"対照: D 0.005 mm・σ_ref 0 → 結果 A・資格あり・更新の条件を満たす ({show(R)})", R["ab"]["classify"]["code"] == "A"
          and R["reference_final"]["qualified"] is True and R["update_conditions_met"] and R["ab"]["extractions"]["sigma_ref"] == 340)


# ====================================================================================== K. プロセス並列の抽出
def test_parallel_extraction():
    print("== K. プロセス並列の抽出 (run_snaps) と直列のビット同一・並列数", flush=True)
    base = TMP / "ns_base"
    if not base.exists():
        make_ns_case(base)
    rec = json.loads((base / HJ.NS_RUN / SR.RECORD).read_text())
    ds = rec["delta_series"]
    ext = {k: v for k, v in ds["extract"].items() if k != "source"}
    eb = base / HJ.EU_RUN
    ns = base / HJ.NS_RUN
    tasks = []
    for k, e in ((S_STOP, 40000), (S_STOP - 2500, 47000), (S_STOP - 5000, 59000), (S_STOP, 60000)):
        tasks.append((ns, k, {"run": str(eb), "res": f"res_{e}.h5", "sha256": HJ.sha256_file(eb / f"res_{e}.h5")}, ext, ds["x"]))
    tasks.append((ns, S_STOP, {"run": str(eb), "res": "res_41000.h5", "sha256": "0" * 64}, ext, ds["x"]))     # sha の不一致 → 例外
    r1, i1 = HJ.run_snaps(tasks, SR.delta_snapshot, 1)
    r3, i3 = HJ.run_snaps(tasks, SR.delta_snapshot, 3)
    strip = lambda r: {k: v for k, v in r["ok"].items() if k != "seconds"} if "ok" in r else r  # noqa: E731
    same = [strip(a) == strip(b) for a, b in zip(r1, r3)]
    check(f"run_snaps: 直列 ({i1['mode']} {i1['jobs']}、{i1['seconds']} s) と並列 ({i3['mode']} {i3['jobs']}、{i3['seconds']} s) で δ_E・δ_target・held・"
          f"missing・λ・出力の sha がビット同一、例外 (sha の不一致) も同じ文で返る (どちらも子プロセス・BLAS 1 スレッド)",
          all(same) and i1["mode"].startswith("process") and i1["jobs"] == 1 and i3["mode"].startswith("process") and i3["jobs"] == 3
          and i1.get("blas_threads") == i3.get("blas_threads") == 1 and "exc" in r1[-1]
          and r1[-1] == r3[-1] and "sha256" in r1[-1]["exc"] and all("ok" in r for r in r1[:-1]), str(same))
    rin = [HJ._call_snap(SR.delta_snapshot, t) for t in tasks[:2]]
    check("記録: この試験の小さい場では、このプロセス (BLAS 既定のスレッド数) の抽出も子の抽出とビット同一 (実寸の場では BLAS のスレッド数で "
          "δ_E が最大 6.9e−18 r_t 違うので、判定は常に子で回す)", [strip(a) == strip(b) for a, b in zip(rin, r1[:2])] == [True, True])
    check("run_snaps: 抽出が差し替え (試験の偽物) なら並列を指定してもこのプロセスで直列 (子に渡せない)",
          HJ.run_snaps(tasks[:2], fake_snap(lambda *a: 0.0), 3)[1]["mode"] == "in-process (serial)")
    with Patch((HJ, "mem_available_mib", lambda: 600.0)):
        j_lo = HJ.snap_jobs(None)
    with Patch((HJ, "mem_available_mib", lambda: 64000.0)):
        j_hi, j_4 = HJ.snap_jobs(None), HJ.snap_jobs(4)
    with Patch((HJ, "mem_available_mib", lambda: None)):
        j_none = HJ.snap_jobs(None)
    check(f"並列数 = min(上限 (既定 3、--jobs 1〜4), 空きメモリ / {HJ.SNAP_PEAK_MIB} MiB): 空き 600 MiB → {j_lo['jobs']}、64 GB → {j_hi['jobs']}・"
          f"--jobs 4 → {j_4['jobs']}、空きが読めない → {j_none['jobs']}; --jobs 5 は拒む",
          (j_lo["jobs"], j_hi["jobs"], j_4["jobs"], j_none["jobs"]) == (2, 3, 4, 1) and raises(lambda: HJ.snap_jobs(5), ValueError) is not None)


# ====================================================================================== G. h370_run.py の段取り
def fake_prepare(problem, run):
    run = Path(run)
    if run.exists():
        raise FileExistsError(run)
    run.mkdir(parents=True)
    pe = RA.load_problem(problem)
    cfg = RA._apply_gas_to_config(RA._config_euler_node(pe, HR.EU_STEPS, HR.EU_OUT, 4.0, 1), pe, run)
    (run / "solverConfig.yaml").write_text(cfg)
    (run / "bcondConfig.yaml").write_text("wall: {kind: slip}\n")
    for f in ("nozzle.h5", "nozzle.msh"):
        (run / f).write_bytes(b"")
    (run / "MESH_QUALITY.txt").write_text("VERDICT: PASS (AR<= 1000, skew<= 0.90)\n")
    info = {"mesh": {"ni": 2000, "nj": 97, "params": {"throat_refine": 5.0}, "hashes": {"generator_sha256": "g" * 64}},
            "physical_upstream": {"requested": RA._physical_upstream(pe.geometry), "wall": {}}, "scale_m": 0.0766, "nStepOuter": HR.EU_STEPS,
            "euler_wall": {"kind": "design_with_physical_upstream", "delta_r": 0.0, "gate": {"pass": True}, "throat": {"x": 0.0}},
            "initial_line": {"match": {"checked": ["R"], "unchecked": []}}}
    (run / "prepare_info.json").write_text(json.dumps(info))
    return info


def make_fake_prepare_ns(mutate=None):
    def f(problem, run, delta_r_csv=None, offset=None, euler_ref=None, stop_rule=None, **kw):
        run = Path(run)
        if run.exists():
            raise FileExistsError(run)
        p, info = gen_ns_configs(run)
        info.update(delta_r_csv=str(delta_r_csv), stop_rule={"requested": False, "use": False})
        (run / "prepare_info.json").write_text(json.dumps(info))
        (run / "MESH_QUALITY.txt").write_text("VERDICT: PASS (AR<= 5000, skew<= 0.90)\n")
        if mutate:
            mutate(run)
        return info
    return f


def ok_item(*a, **k):
    return {"ok": True, "failures": [], "judged": {"max_abs_mm": 7.96e-4}, "check_mesh_quality_nozzle_h5": "VERDICT: PASS (試験の偽物)"}


def bad_item(*a, **k):
    return {"ok": False, "failures": ["δ_in の照合が不合格 (試験の偽物)"], "judged": {"max_abs_mm": 2e-3}}


def fake_case() -> Path:
    """偽の B の case dir (TMP/repo/case/45x): 問題 YAML・δ_r の CSV・CFD ピンの偽物・h370_judge.py。repo/design・solver_density_cuda は本物へ
    symlink (h370_run.py の HERE.parents[1] の代わり)。"""
    repo = TMP / "repo"
    if not repo.exists():
        (repo / "case").mkdir(parents=True)
        (repo / "design").symlink_to(ROOT / "design")
        (repo / "solver_density_cuda").symlink_to(ROOT / "solver_density_cuda")
    c = repo / "case" / f"c{len(list((repo / 'case').iterdir()))}"
    (c / "_band_ab/half_inlet").mkdir(parents=True)
    for n in (HR.EU_PROBLEM, HR.NS_PROBLEM, HR.CSV, "h370_judge.py"):
        shutil.copy(CASE / n, c / n)
    (c / HR.CFD_PIN_RUN).mkdir()
    for n in HR.CFD_PIN_FILES:
        (c / HR.CFD_PIN_RUN / n).write_text("偽物")
    return c


def args(**kw):
    return types.SimpleNamespace(**{"commit": "abc1234", "dry_run": None, "no_judge": True, "phase": "euler", **kw})


def test_run_orchestration():
    print("== G. h370_run.py の段取り (偽の準備・本物の生成器の設定)", flush=True)
    calls = []
    no_forge = lambda *a, **k: calls.append("forge") or (_ for _ in ()).throw(AssertionError("forge を回した"))  # noqa: E731
    base = [(HR, "resources", lambda *a, **k: {"ok": True, "failures": []}), (HR, "WT", TMP / "no_worktree"),
            (HR, "quality_check", ok_item), (HR, "euler_wall_restore_check", ok_item), (HR, "delta_in_check", ok_item),
            (RA, "prepare", fake_prepare), (RA, "prepare_ns", make_fake_prepare_ns()), (RA, "run_staged", no_forge), (RA, "run_staged_ns", no_forge)]
    c = fake_case()
    with Patch(*base, (HR, "HERE", c)):
        rc = HR.phase_euler(args(dry_run="t1"), RA, {})
        pc = json.loads((c / "_h370_dry_t1_euler" / "H370_PRECHECK.json").read_text())
        check(f"--dry-run: Euler と NS の投入前の検査が合格 (rc {rc})・forge を回さない・run_0950 を作らない",
              rc == 0 and pc["ok"] and not calls and not (c / HR.EU_RUN).exists() and (c / "_h370_dry_t1_ns").is_dir()
              and set(pc["items"]) >= {"E1_wall_restore", "E2_quality", "E3_prepare", "N1_quality", "N2_delta_in", "N3_wall_output", "N4_config", "N5_prepare"},
              str(pc["failures"]))
        check("--dry-run の記録: Euler の段の sha (soft・mid・本段)・NS の段の sha・壁の出力の名前",
              set(pc["items"]["E3_prepare"]["sha256"]) == {"soft", "mid", "main"} and set(pc["items"]["N4_config"]["sha256"]) >= {"S1_soft", "S2_mid", "main"}
              and pc["items"]["N3_wall_output"]["wall_dump_name"] == "res_wall_3_<step>.h5")
        e = raises(lambda: HR.phase_euler(args(dry_run="t1"), RA, {}), SystemExit)
        check("同じ名前の dry-run の準備が既にある → 止める (消さない)", e is not None)
    c = fake_case()
    with Patch(*base, (HR, "HERE", c), (HR, "delta_in_check", bad_item)):
        rc = HR.phase_euler(args(), RA, {})
        check("本番: δ_in の照合が不合格 → 前提の検査の不成立 (rc 2)・PRECHECK_FAILED・forge を回さない",
              rc == 2 and (c / HR.EU_RUN / "PRECHECK_FAILED").is_file() and not calls)
    c = fake_case()
    mut = lambda run: (run / "bcondConfig.yaml").write_text((run / "bcondConfig.yaml").read_text().replace("outputHDFflg: 1, ints: , floats: {Ux", "outputHDFflg: 0, ints: , floats: {Ux"))  # noqa: E731
    with Patch(*base, (HR, "HERE", c), (RA, "prepare_ns", make_fake_prepare_ns(mut))):
        rc = HR.phase_euler(args(), RA, {})
        check("本番: NS が停止の出力で壁の出力を書かない設定 → 終了コード 4 (親に返す)・forge を回さない", rc == HR.EXIT_WALL_OUTPUT and not calls)
    # 合格 → Euler を回す (偽物) → 資格の判定 (本物の h370_judge.py を子プロセスで) が判定不能 → NS に進まない
    c = fake_case()

    def fake_staged(run, **kw):
        calls.append(("run_staged", kw))
        return 0
    with Patch(*base, (HR, "HERE", c), (RA, "run_staged", fake_staged), (RA, "collect", lambda *a: {"偽物": True})):
        rc = HR.phase_euler(args(phase="all"), RA, {})
        rr = json.loads((c / HR.EU_RUN / "RUN_END.json").read_text())
        check(f"本番: 検査の合格 → run_staged を cfl_main 2.0・mid_stage・full で呼ぶ → 資格が無い (偽の Euler) → rc 2 で NS に進まない",
              rc == 2 and calls and calls[-1][1] == {"cfl_main": 2.0, "mid_stage": True, "stages": "full"} and rr["run_rc"] == "0"
              and (c / HR.QUAL).is_file() and not json.loads((c / HR.QUAL).read_text())["prequalified"] and not (c / HR.PIN).exists(), str(calls[-1:]))
        calls.clear()
        e = raises(lambda: HR.phase_ns(args(phase="ns"), RA, {}), SystemExit)
        check("固定が無い → NS を準備しない (止める)", e is not None and not (c / HR.NS_RUN).exists())
    # 環境・prep.done の検査
    with Patch((HR, "HERE", c)):
        old = {k: v for k, v in os.environ.items() if k.startswith("FORGE_")}
        try:
            for k in old:
                os.environ.pop(k)
            os.environ.update(HR.ENV_WANT)
            os.environ["FORGE_EXTRA"] = "1"
            e1 = raises(lambda: HR.code_checks(args()), SystemExit)
            os.environ.pop("FORGE_EXTRA")
            e2 = raises(lambda: HR.code_checks(args()), SystemExit)
        finally:
            for k in [k for k in os.environ if k.startswith("FORGE_")]:
                os.environ.pop(k)
            os.environ.update(old)
        check("FORGE_* が登録と違う → 止める; h370_prep.done が無い → 止める", e1 is not None and e2 is not None)
    # h370_problems.py check を作業ツリーの版で
    c = fake_case()
    (c / "_band_ab").mkdir(exist_ok=True)
    with Patch((HR, "HERE", c), (HR, "WT", ROOT), (HR, "load_forge_design", lambda commit: (RA, {}))):
        e = raises(lambda: HR.check_problems(args(phase="check-problems")), SystemExit)
        check("check-problems: V5d の起動の記録 (moc_v5d_run.log) が B に無い → 止める", e is not None)
        shutil.copy(CASE / "_band_ab" / "moc_v5d_run.log", c / "_band_ab" / "moc_v5d_run.log")
        rc = HR.check_problems(args(phase="check-problems"))
        out = HJ.jload(c / "_band_ab/half_inlet/h370_problems_B.json")
        check("check-problems: 作業ツリーの h370_problems.py check が PASS、出力は h370_problems_B.json (commit した json を上書きしない)",
              rc == 0 and out and out["VERDICT"] == "PASS")


# ====================================================================================== H. NS の段取り・部品
def test_ns_phase_and_parts():
    print("== H. NS の段取り (偽の準備・偽の run_staged_ns で本段の直前の照合)・検査の部品", flush=True)
    from forge_design.feedback.deltastar_loop import handoff_extract_settings
    c = fake_case()
    eb = c / HR.EU_RUN
    eb.mkdir()
    (eb / "res_60000.h5").write_bytes(b"euler")
    sha = HR.sha256_file(eb / "res_60000.h5")
    for k in HR.SIGMA_REF_STEPS:
        (eb / f"res_{k}.h5").write_bytes(f"euler {k}".encode())
    (c / HR.PIN).write_text(json.dumps({"run": HR.EU_RUN, "res": "res_60000.h5", "sha256": sha}))
    (c / HR.QUAL).write_text(json.dumps({"prequalified": True, "pin": {"sha256": sha},
                                         "sigma_ref_snapshots": [{"res": f"res_{k}.h5", "step": k, "sha256": HR.sha256_file(eb / f"res_{k}.h5")}
                                                                 for k in HR.SIGMA_REF_STEPS]}))
    fake_ns = make_fake_prepare_ns()

    def fake_prepare_ns_real(problem, run, delta_r_csv=None, offset=None, euler_ref=None, stop_rule=None, **kw):
        info = fake_ns(problem, run, delta_r_csv=delta_r_csv)
        sr = {"requested": True, "use": True, "cap": HR.NS_CAP, "hold_steps": HR.NS_HOLD, "out_step_interval": HR.NS_OUT, **HR.THRESHOLDS,
              "extract": {**handoff_extract_settings(), "source": "試験"}, "euler_ref": euler_ref}
        info.update(stop_rule=sr, euler_ref=euler_ref)
        (Path(run) / "prepare_info.json").write_text(json.dumps(info))
        return info
    seen = {}

    def fake_staged_ns(run, stages="full", **kw):
        """run_staged_ns の段の導出 → 本段の通知 (Observer が照合) だけをまねる (forge は回さない)。"""
        run = Path(run)
        cfg_main = (run / "solverConfig.yaml").read_text()
        cfg_pre = (run / RA.PRESTAGE_CONFIG).read_text()
        pre, main, _ = RA.ns_stage_configs(cfg_main, stages, cfg_pre=cfg_pre)
        done = [{"tag": t, "nsteps": n, "config_sha256": HR.sha256_text(cfg)} for t, cfg, n, _ in pre]
        main = RA._cfg_set(main, {"nStepOuter": HR.NS_CAP})
        (run / "solverConfig.yaml").write_text(main)
        seen["rec"] = RA._write_main_stage_started(run, main, stages, seen.get("done_override", done), cfg_pre)
        return type("RC", (int,), {"stop_rule": {"status": "REACHED"}})(0)
    judge_calls = []
    real_run = subprocess.run

    def spy_run(cmd, *a, **k):
        if any("h370_judge.py" in str(x) for x in cmd):
            judge_calls.append(cmd)
            return types.SimpleNamespace(returncode=0, stdout="", stderr="")
        return real_run(cmd, *a, **k)
    pin_rec = {"run": str(c / HR.EU_RUN), "res": "res_60000.h5", "sha256": sha}
    with Patch((HR, "resources", lambda *a, **k: {"ok": True, "failures": []}), (HR, "WT", TMP / "no_worktree"), (HR, "HERE", c),
               (HR, "quality_check", ok_item), (HR, "delta_in_check", ok_item), (RA, "prepare_ns", fake_prepare_ns_real),
               (RA, "run_staged_ns", fake_staged_ns), (subprocess, "run", spy_run)):
        rc = HR.phase_ns(args(phase="ns", no_judge=False), RA, {})
    run = c / HR.NS_RUN
    pc = json.loads((run / "H370_PRECHECK.json").read_text())
    ms = HJ.jload(run / "H370_MAIN_START.json") or {}
    n7 = pc["items"]["N7_same_as_prelaunch"]
    other = [f for f in pc["failures"] if not f.startswith("(i) の投入前の準備")]
    check(f"NS: (i) の投入前の準備 {HR.NS_PRELAUNCH} が無い → 照合できない → 前提の検査の不成立 (rc {rc})・forge を回さない",
          rc == 2 and not n7["ok"] and not other and (run / "PRECHECK_FAILED").is_file() and "rec" not in seen and not judge_calls, str(other)[:300])
    # (i) の準備があるとき: 合格 → 本段の直前の照合 → 判定を呼ぶ
    c2 = fake_case()
    shutil.copytree(eb, c2 / HR.EU_RUN)
    for f in (HR.PIN, HR.QUAL):
        shutil.copy(c / f, c2 / f)
    pin_rec = {"run": str(c2 / HR.EU_RUN), "res": "res_60000.h5", "sha256": sha}
    with Patch((HR, "HERE", c2)):
        fake_ns(c2 / HR.NS_PROBLEM, c2 / HR.NS_PRELAUNCH, delta_r_csv=str(c2 / HR.CSV))
    base = [(HR, "resources", lambda *a, **k: {"ok": True, "failures": []}), (HR, "WT", TMP / "no_worktree"), (HR, "HERE", c2),
            (HR, "quality_check", ok_item), (HR, "delta_in_check", ok_item), (RA, "prepare_ns", fake_prepare_ns_real),
            (RA, "run_staged_ns", fake_staged_ns), (subprocess, "run", spy_run)]
    with Patch(*base):
        rc = HR.phase_ns(args(phase="ns", no_judge=False), RA, {})
    run = c2 / HR.NS_RUN
    pc = json.loads((run / "H370_PRECHECK.json").read_text())
    ms = HJ.jload(run / "H370_MAIN_START.json") or {}
    check(f"NS: 投入前の検査 (停止規則・Euler の参照の固定・設定・壁の出力・(i) の準備とビット同一) → 本段の直前の照合が合格 (rc {rc})・判定を呼ぶ",
          rc == 0 and pc["ok"] and ms.get("ok") and pc["items"]["N5_prepare"]["euler_ref"] == pin_rec and pc["items"]["N7_same_as_prelaunch"]["ok"]
          and judge_calls, str(pc["failures"] + ms.get("failures", []))[:300])
    check("本段の直前の照合の記録: 本段の設定の sha = 準備の時点で導いた本段 (上限 200,000)・前の段 S1・S2 の sha",
          ms.get("main_stage_started", {}).get("main_config", {}).get("nStepOuter") == HR.NS_CAP
          and [s["tag"] for s in ms["main_stage_started"]["previous_stages"]] == ["S1_soft", "S2_mid"])
    # 負例: 前の段の記録が導いた S1・S2 と違う → WATCH_ABORT・本段を起動させない
    shutil.rmtree(run)
    seen["done_override"] = [{"tag": "S1_soft", "nsteps": 3000, "config_sha256": "0" * 64}]
    with Patch(*base):
        rc = HR.phase_ns(args(phase="ns"), RA, {})
    ms = HJ.jload(run / "H370_MAIN_START.json") or {}
    rr = json.loads((run / "RUN_END.json").read_text())
    check("本段の直前の照合: 前の段の記録が導いた S1・S2 と違う → WATCH_ABORT・例外 (本段の forge を起動させない)",
          not ms.get("ok") and (run / "WATCH_ABORT").is_file() and rr["run_rc"] == "EXC" and "WATCH_ABORT" in rr.get("exception", ""))
    seen.pop("done_override")
    # 部品
    with Patch((HR, "HERE", c2)):
        r = HR.read_pin()
        check("read_pin: 固定 {run, res, sha256} を run の絶対パスで返す (σ_ref の代わりの参照 = 末尾の窓の全 20 出力も在って sha が合う)",
              r["spec"] == pin_rec)
        for k in (45000, 41000):
            (c2 / HR.EU_RUN / f"res_{k}.h5").unlink()
            check(f"read_pin: σ_ref の代わりの参照の出力 ({k}{'、旧い代表の外' if k % 5000 else ''}) が消されている → NS を回す前に止める",
                  raises(HR.read_pin, SystemExit) is not None)
            (c2 / HR.EU_RUN / f"res_{k}.h5").write_bytes(f"euler {k}".encode())
        q_ = json.loads((c2 / HR.QUAL).read_text())
        (c2 / HR.QUAL).write_text(json.dumps({**q_, "sigma_ref_snapshots": []}))
        check("read_pin: 予備の資格の記録に σ_ref の代わりの参照が無い → 止める", raises(HR.read_pin, SystemExit) is not None)
        (c2 / HR.QUAL).write_text(json.dumps({**q_, "sigma_ref_snapshots": [a_ for a_ in q_["sigma_ref_snapshots"] if a_["step"] % 5000 == 0]}))
        e = raises(HR.read_pin, SystemExit)
        check("read_pin: 予備の資格の記録の代わりの参照が旧い代表 4 つだけ (§5.1 #36g の記録) → 全 20 出力でないので止める", e is not None)
        (c2 / HR.QUAL).write_text(json.dumps(q_))
        check("read_pin: 元に戻せば通る", HR.read_pin()["spec"] == pin_rec)
        (c2 / HR.EU_RUN / "res_60000.h5").write_bytes(b"changed")
        check("read_pin: 固定した出力の sha が違う → 止める", raises(HR.read_pin, SystemExit) is not None)
    a, b = TMP / "mi_a.h5", TMP / "mi_b.h5"
    for p_, ro in ((a, 1.0), (b, 2.0)):
        with h5py.File(p_, "w") as f:
            f["MESH/COORD"] = np.arange(6.0)
            f["VALUE/wall_dist"] = np.ones(2)
            f["VALUE/ro"] = np.full(2, ro)
    r = HR.mesh_identity(a, b)
    check("mesh_identity: 格子・wall_dist が同じで初期の場だけ違う → 合格 (初期の場の違いは記録)", r["ok"] and r["initial_field_differ (記録だけ)"] == ["VALUE/ro"])
    with h5py.File(b, "r+") as f:
        f["MESH/COORD"][0] = 1e-12
    check("mesh_identity: 座標が 1 か所でも違う → 不成立", not HR.mesh_identity(a, b)["ok"])
    rd = TMP / "timing_run"
    rd.mkdir()
    (rd / "forge_run.log").write_text("step 1900\nHighFive::DataSetException: Write failed\n")
    HR.Timing(rd).run_once()
    check("Timing: forge_run.log の Write failed → WRITE_FAILED を残す・timing.tsv を書く", (rd / "WRITE_FAILED").is_file() and (rd / "timing.tsv").is_file())
    fake_tools = TMP / "fake_tools"
    fake_tools.mkdir()
    rq = TMP / "q_run"
    rq.mkdir()
    (rq / "MESH_QUALITY.txt").write_text("# cell\nVERDICT: PASS (AR<= 5000, skew<= 0.90)\n")
    for verdict, want in (("PASS (AR<= 5000, skew<= 0.90)", True), ("SOFT-PASS (<0.1% outliers) (AR<= 5000, skew<= 0.90)", False)):
        (fake_tools / "check_mesh_quality.py").write_text(f"print('VERDICT: {verdict}')\n")
        r = HR.quality_check(types.SimpleNamespace(FORGE_TOOLS=fake_tools), rq, 5000)
        check(f"quality_check: check_mesh_quality の VERDICT {verdict.split()[0]} → {'合格' if want else '不成立 (SOFT-PASS は通さない)'}", r["ok"] == want)
    (rq / "MESH_QUALITY.txt").write_text("VERDICT: PASS (AR<= 1000, skew<= 0.90)\n")
    (fake_tools / "check_mesh_quality.py").write_text("print('VERDICT: PASS (AR<= 5000, skew<= 0.90)')\n")
    check("quality_check: prepare の MESH_QUALITY.txt の --ar-max が違う (1000) → 不成立",
          not HR.quality_check(types.SimpleNamespace(FORGE_TOOLS=fake_tools), rq, 5000)["ok"])
    # δ_in の照合 (本物の check_runs、合成の構造格子)
    tbl = np.loadtxt(CASE / HR.CSV, delimiter=",", skiprows=1)
    drx = RA.delta_r_from_table(tbl[:, 0], tbl[:, 1])
    x = np.linspace(-12.5, 95.29, 300)
    rE = 1.0 + 0.0009 * x * x
    for name, shift_mm, want in (("ok", 0.0, True), ("bad", 0.002, False)):
        for kind, rw in (("eu", rE), ("ns", rE + drx(x) + shift_mm * 1e-3 / 0.0767)):
            d = TMP / f"din_{name}_{kind}"
            d.mkdir()
            R2 = rw[:, None] * np.linspace(0.0, 1.0, 5)[None, :]
            with h5py.File(d / "nozzle.h5", "w") as f:
                f["/MESH/COORD"] = np.c_[np.repeat(x, 5) * 0.0767, R2.ravel() * 0.0767, np.zeros(R2.size)].ravel()
            (d / "prepare_info.json").write_text(json.dumps({"mesh": {"ni": 300, "nj": 5}, "scale_m": 0.0767}))
        with Patch((HR, "HERE", CASE)):
            r = HR.delta_in_check(TMP / f"din_{name}_ns", TMP / f"din_{name}_eu")
        check(f"delta_in_check (本物の check_runs): NS の壁 = Euler の壁 + δ_r{'' if want else ' + 0.002 mm'} → {'合格' if want else '不合格'} "
              f"(max {r.get('judged', {}).get('max_abs_mm')} mm)", r["ok"] == want and r.get("csv_sha256") == HR.CSV_SHA, str(r.get("failures")))


def find_pin():
    for c in (CASE / HR.CFD_PIN_RUN, Path("/home/sano/work/forge/case/45.isobutane_m6_d155") / HR.CFD_PIN_RUN):
        if (c / "prepare_info.json").is_file():
            return c.resolve()
    return None


def test_euler_wall_restore():
    print("== I. 正式の Euler の壁の保存・復元の照合 (本物の設計・CFD ピン)", flush=True)
    pin = find_pin()
    if pin is None:
        check("CFD ピンの凍結元 run_0062 が手元に無い — I を回せない", False)
        return
    import h370_problems as HP
    from forge_design.geometry.wall_axismach import save_wall_file
    from forge_design.meshing.mesh2d import generate_axisym_mesh, write_msh41_2d
    td = TMP / "pinned"
    td.mkdir()
    pe = HP.load_pinned(HP.EU, pin, td)
    de = RA.design_chain(pe)
    scale = float(pe.spec["r_throat"])
    wE, er = RA.euler_reference_wall(pe, de, scale)
    c = fake_case()
    run = c / HR.EU_RUN
    run.mkdir()
    mp = RA.mesh_params_euler(pe, scale)
    co, q, b = generate_axisym_mesh(wE, mp)
    write_msh41_2d(run / "nozzle.msh", co, q, b)
    wpath, wrec, wsha = save_wall_file(run, wE, scale, er["repr"])
    info = {"scale_m": scale, "euler_wall": {**er, "wall_file": {"repr": er["repr"], "file": wpath.name, "sha256": wsha}},
            "mesh": {"hashes": {"generator_sha256": HR.generator_sha256(co, q, b)}}}
    (run / "prepare_info.json").write_text(json.dumps(info, default=str))
    shutil.copy(CASE / HR.PREP_CHECK_REC, c / HR.PREP_CHECK_REC)
    with Patch((HR, "HERE", c)):
        r = HR.euler_wall_restore_check(RA, run)
        loc = r.get("local_prep_check (記録だけ)", {})
        check(f"壁ファイル (design_upstream) → 復元した壁で生成した格子の msh が nozzle.msh とバイト一致・生成器の sha が prepare_info と同じ "
              f"(手元の prep-check の記録と同じ: 壁 {loc.get('wall_file_sha256_same')}・格子 {loc.get('generator_sha256_same')})",
              r["ok"] and r["msh_bytes_identical"] and loc.get("wall_file_sha256_same") and loc.get("generator_sha256_same"), str(r.get("failures")))
        m = (run / "nozzle.msh").read_text().splitlines()
        m[len(m) // 2] = m[len(m) // 2] + " "
        (run / "nozzle.msh").write_text("\n".join(m) + "\n")
        r = HR.euler_wall_restore_check(RA, run)
        check("負例: nozzle.msh が 1 文字違う → 不成立", not r["ok"] and not r["msh_bytes_identical"])
        info["euler_wall"]["wall_file"]["sha256"] = "0" * 64
        (run / "prepare_info.json").write_text(json.dumps(info, default=str))
        r = HR.euler_wall_restore_check(RA, run)
        check("負例: 壁ファイルの sha が prepare_info の記録と違う → 不成立", not r["ok"])


# ====================================================================================== L. 診断の経路 (§5.1 #36l・§6 V-h370-diag)
S_DIAG = 62500                                                          # 合成の診断の NS の停止の時刻 (窓 [22,500, 62,500] の 17 出力)
ARM_RUNS = {"A": "run_0952_h370_euler_cfl2", "B": "run_0953_h370_euler_cfl05"}


class GateCalled(AssertionError):
    pass


def gate_boom(calls):
    """受け渡しの関門の偽物: 呼ばれたら記録して例外 (診断の経路では呼ばれてはいけない)。"""
    def g(*a, **k):
        calls.append(("handoff_gate", a, k))
        raise GateCalled("受け渡しの関門が呼ばれた")
    return g


def diag_fake_case(qual_prequalified=False, pin=False) -> Path:
    """診断の経路の偽の B の case dir: Euler run_0950 (res_60000)・腕 A・B (res_40000)・予備の資格の記録 (prequalified 偽、末尾の出力の sha)・
    (i) の投入前の準備 _h370_prelaunch_ns (本物の生成器の設定・偽の格子)。"""
    c = fake_case()
    for rn, res in ((HR.EU_RUN, "res_60000.h5"), *HR.DIAG_SENS_REFS.values()):
        (c / rn).mkdir(exist_ok=True)
        (c / rn / res).write_bytes(f"euler {rn} {res}".encode())
    sha60 = HR.sha256_file(c / HR.EU_RUN / "res_60000.h5")
    (c / HR.QUAL).write_text(json.dumps({"prequalified": qual_prequalified, "rc": 1, "failures": ["check_convergence.py --segment --drop 3 が PASS でない"],
                                         "undecidable": [], "tail_outputs": [{"res": "res_60000.h5", "step": 60000, "sha256": sha60}]}))
    if pin:
        (c / HR.PIN).write_text(json.dumps({"run": HR.EU_RUN, "res": "res_60000.h5", "sha256": sha60}))
    with Patch((HR, "HERE", c)):
        make_fake_prepare_ns()(c / HR.NS_PROBLEM, c / HR.NS_PRELAUNCH, delta_r_csv=str(c / HR.CSV))
    return c


def fake_prepare_ns_sr(problem, run, delta_r_csv=None, offset=None, euler_ref=None, stop_rule=None, **kw):
    """prepare_ns の偽物 (本物の生成器の設定・偽の格子) に停止規則と Euler の参照の記録を足す (test H と同じ形)。"""
    from forge_design.feedback.deltastar_loop import handoff_extract_settings
    info = make_fake_prepare_ns()(problem, run, delta_r_csv=delta_r_csv)
    sr = {"requested": True, "use": True, "cap": HR.NS_CAP, "hold_steps": HR.NS_HOLD, "out_step_interval": HR.NS_OUT, **HR.THRESHOLDS,
          "extract": {**handoff_extract_settings(), "source": "試験"}, "euler_ref": euler_ref}
    info.update(stop_rule=sr, euler_ref=euler_ref)
    (Path(run) / "prepare_info.json").write_text(json.dumps(info))
    return info


def make_fake_staged_ns_diag(status="REACHED", last=S_DIAG, raise_exc=False, seen=None):
    """run_staged_ns の偽物: 段の導出 → 本段の通知 (Observer が照合) → 本段の出力 (場・壁・出口の境界の出力、2,500 ごと) を last まで書き、
    停止規則の記録と状態を書く (forge は回さない)。何も消さない。"""
    def f(run, stages="full", **kw):
        run = Path(run)
        cfg_main = (run / "solverConfig.yaml").read_text()
        cfg_pre = (run / RA.PRESTAGE_CONFIG).read_text()
        pre, main, _ = RA.ns_stage_configs(cfg_main, stages, cfg_pre=cfg_pre)
        done = [{"tag": t, "nsteps": n, "config_sha256": HR.sha256_text(cfg)} for t, cfg, n, _ in pre]
        main = RA._cfg_set(main, {"nStepOuter": HR.NS_CAP})
        (run / "solverConfig.yaml").write_text(main)
        RA._write_main_stage_started(run, main, stages, done, cfg_pre)
        if raise_exc:
            raise RuntimeError("偽の見張りの例外")
        steps = list(range(HR.NS_OUT, last + 1, HR.NS_OUT))
        for k in steps:
            (run / f"res_{k}.h5").write_bytes(f"field {k}".encode())
            (run / f"res_{k}.xmf").write_text(f"xmf {k}")
            (run / f"res_wall_3_{k}.h5").write_bytes(f"wall {k}".encode())
            (run / f"res_outlet_2_{k}.h5").write_bytes(f"outlet {k}".encode())
        (run / "stop_rule_band").mkdir(exist_ok=True)
        (run / "stop_rule_band" / "delta_r_edge.csv").write_text("x_rt,y_ref,y_star,kappa_med,kappa_s,y_b\n")
        hold = {"hold_steps": HR.NS_HOLD, "start": last - HR.NS_HOLD if status == "REACHED" else None, "reached": [10000, last - HR.NS_HOLD],
                "breaks": [{"start": 10000, "broken_at": 12500, "why": "試験"}]}
        SR.save_record(run, {"version": 1, "status": status, "rows": [{"step": k} for k in steps], "stop_step": last, "hold": hold,
                             "band": {"step": 2500, "res_sha256": "x", "euler_ref": "試験", "out_dir": str(run / "stop_rule_band")}})
        state = {"status": status, "stop_step": last, "last_output": last, "forge_rc": 143 if status == "REACHED" else 0, "hold": hold}
        SR.write_state(run, state)
        if seen is not None:
            seen["steps"] = steps
        return type("RC", (int,), {"stop_rule": state})(0)
    return f


def write_diag_case(case: Path, shape_b: float = 0.0, status: str = "REACHED", s_stop: int = S_DIAG) -> Path:
    """診断の感度の判定の合成の case dir: Euler の参照 aux (run_0950・res_60000)・腕 A (run_0952・res_40000)・腕 B (run_0953・res_40000、形
    shape_b)、診断の NS (run_0951、出力 17,500〜s_stop、窓 2 × 20,000 → 17 出力、停止規則の記録・状態・印・投入前の検査 (ns-diag・版)・
    本段の直前の照合・sha の一覧・終了の記録)、診断の参照の固定 (DIAG_REF)。"""
    case.mkdir(parents=True)
    refs = {}
    for key, rn, step, shape in (("aux", HJ.EU_RUN, 60000, 0.0), ("A", ARM_RUNS["A"], 40000, 0.0), ("B", ARM_RUNS["B"], 40000, shape_b)):
        X, R, Q = synth_eu(X1 + 0.5, shape)
        write_field_run(case / rn, X, R, Q, {}, (step,))
        refs[key] = {"run": rn, "res": f"res_{step}.h5", "sha256": HJ.sha256_file(case / rn / f"res_{step}.h5")}
    (case / "_band_ab/half_inlet").mkdir(parents=True)
    (case / HJ.DIAG_REF).write_text(json.dumps({"euler_ref": {**refs["aux"], "step": 60000}, "sens_refs": {"A": refs["A"], "B": refs["B"]},
                                                "qualified": False}))
    spec_aux = {"run": str(case / HJ.EU_RUN), "res": "res_60000.h5", "sha256": refs["aux"]["sha256"]}
    ns = case / HJ.NS_RUN
    sr = SR.resolve_settings({"stop_rule": True, "stop_rule_thresholds": {"window": 20000}}, nsteps_arg=None, n_main=200000,
                             n_main_source="試験", out_int=2500, euler_ref=spec_aux)
    steps = list(range(s_stop - 45000, s_stop + 1, 2500))
    X, R, Q = synth_ns()
    write_field_run(ns, X, R, Q, {"wall_thermal": {"mode": "isothermal", "Tw": 370.0}, "physical_wall": {"domain": [X0, X1]},
                                  "euler_ref": spec_aux, "stop_rule": sr}, steps, perturb=lambda st: 1.0 + 1e-6 * (st - s_stop) / 2500.0)
    uc = SR.use_columns(ns)
    rows = [{"step": st, "r": 1e-5, "exit_core_M": 5.99, "mdot_rt": 10.0} for st in steps]
    rec = {"version": 1, "settings": sr, "status": status, "rows": rows, "stop_step": s_stop,
           "delta_series": {"x": uc["x"].tolist(), "use_domain": [0.0, uc["x_w"]], "scale_m": uc["scale_m"], "extract": sr["extract"],
                            "euler_ref": sr["euler_ref"], "span_steps": 40000}}
    SR.save_record(ns, rec)
    SR.write_state(ns, {"status": status, "stop_step": s_stop, "forge_rc": 143, "stop_sha256": HJ.sha256_file(ns / f"res_{s_stop}.h5")})
    (ns / HJ.DIAG_MARK).write_text(json.dumps({"diagnostic_only": True, "update_conditions_met": False}))
    (ns / "H370_PRECHECK.json").write_text(json.dumps({"ok": True, "phase": "ns-diag", "versions": HJ.module_versions()}))
    (ns / "H370_MAIN_START.json").write_text(json.dumps({"ok": True}))
    (ns / "RUN_RC").write_text("0\n")
    (ns / "RUN_PROVENANCE.txt").write_text("forge_sha256: 7e95e5bbd8f50c9e0000\nenv         : FORGE_CUDA_BLOCKSIZE=128\n")
    (ns / HJ.DIAG_END).write_text(json.dumps({"end_state": status, "stop_step": s_stop}))
    (ns / HJ.DIAG_SHA).write_text(json.dumps(HR.files_sha256(ns, skip=(HJ.DIAG_SHA, HJ.DIAG_END))))
    return case


def fake_diag_snap(effect=None, mutate=None, calls=None):
    """偽物の診断の抽出 (速い): δ_E = δ_target = 1e-3 r_t + effect(腕, 出力, x) [mm] / MMPR、帯の位置は一定。mutate(腕, 出力, d) で held・欠測・
    例外を作る。腕は Euler の run の名前で見分ける (aux = run_0950)。"""
    def g(run, step, spec, ext, x_use):
        x = np.asarray(x_use, dtype=float)
        name = Path(spec["run"]).name
        arm = "aux" if name == HJ.EU_RUN else ("A" if name == ARM_RUNS["A"] else "B")
        if calls is not None:
            calls.append((arm, int(step)))
        e = 1e-3 + (np.asarray(effect(arm, int(step), x), dtype=float) / MMPR if effect else 0.0) + 0.0 * x
        d = {"delta_E": e.tolist(), "delta_target": e.tolist(), "held": [0] * len(x), "missing": [], "n_missing": 0, "n_held": 0, "lam": 1.0,
             "res_sha256": None, "seconds": 0.0, "y_b": (0.1 + 0.0 * x).tolist(), "y_star": (0.08 + 0.0 * x).tolist()}
        if mutate:
            d = mutate(arm, int(step), d) or d
        return d
    return g


def diag_judge(case: Path, extra=()) -> tuple:
    rc = HJ.main(["--diag-sens", "--case-dir", str(case), "--design", str(ROOT / "design"), "--tools", str(TOOLS), *extra])
    return rc, HJ.jload(case / HJ.DIAG_OUT)


def test_diag_path():
    print("== L. 診断の経路 (--phase ns-diag・--diag-sens; §5.1 #36l・§6 V-h370-diag)", flush=True)
    import inspect
    gate_calls, judge_calls = [], []
    real_run = subprocess.run

    def spy_run(cmd, *a, **k):
        if any("h370_judge.py" in str(x) for x in cmd):
            judge_calls.append(cmd)
            return types.SimpleNamespace(returncode=0, stdout="", stderr="")
        return real_run(cmd, *a, **k)
    # --- 登録の値の写し ---
    p_args = None
    try:
        HR.main(["--phase", "ns-diag", "--commit", "zz"])
    except SystemExit as e:
        p_args = e.code
    check("L 定数: h370_judge.py の DIAG_* (参照の固定・印・終了の記録・sha の一覧・版のファイル) が h370_run.py と同じ・判定の窓 2 × 20,000・"
          "出力 2,500・17 出力・腕 A/B の run が CFL の A/B (h370_euler_cfl_ab.sh) と同じ・--phase ns-diag を受け付ける",
          (HJ.DIAG_REF, HJ.DIAG_MARK, HJ.DIAG_END, HJ.DIAG_SHA, HJ.VERSION_FILES) == (HR.DIAG_REF, HR.DIAG_MARK, HR.DIAG_END, HR.DIAG_SHA, HR.VERSION_FILES)
          and HJ.DIAG_WINDOW == HR.THRESHOLDS["window"] and HJ.DIAG_OUT_INT == HR.NS_OUT
          and len(SR.required_window_steps(S_DIAG, HJ.DIAG_WINDOW, HJ.DIAG_OUT_INT)) == HJ.DIAG_N_REQ == 17
          and {a_: r_ for a_, (r_, _) in HR.DIAG_SENS_REFS.items()} == ARM_RUNS
          and all(r_ in (CASE / "h370_euler_cfl_ab.sh").read_text() for r_ in ARM_RUNS.values())
          and HR.DIAG_REF != HR.QUAL and HR.DIAG_REF != HR.PIN and p_args == 1, f"--phase ns-diag の終了 {p_args}")
    # --- 診断の経路の本番 (偽の準備・偽の本段): (a) 関門を呼ばない・(b) 正式の記録を書き換えない・(f) 途中の出力を消さない ---
    c = diag_fake_case()
    qual_before = (c / HR.QUAL).read_bytes()
    seen = {}
    base = [(HR, "resources", lambda *a, **k: {"ok": True, "failures": []}), (HR, "WT", TMP / "no_worktree"), (HR, "HERE", c),
            (HR, "quality_check", ok_item), (HR, "delta_in_check", ok_item), (HR, "disk_free", lambda: 64 * 2 ** 30),
            (RA, "prepare_ns", fake_prepare_ns_sr), (subprocess, "run", spy_run), (DL, "handoff_gate", gate_boom(gate_calls)),
            (HJ, "HOOKS", {"handoff_gate": gate_boom(gate_calls)})]
    with Patch(*base, (RA, "run_staged_ns", make_fake_staged_ns_diag(seen=seen))):
        rc = HR.phase_ns_diag(args(phase="ns-diag", no_judge=False), RA, {})
    run = c / HR.NS_RUN
    pc = json.loads((run / "H370_PRECHECK.json").read_text())
    end = json.loads((run / HR.DIAG_END).read_text())
    shal = json.loads((run / HR.DIAG_SHA).read_text())
    dref = json.loads((c / HR.DIAG_REF).read_text())
    check(f"(a) 診断の経路の本番 (rc {rc}): 投入前の検査が合格 (phase ns-diag・版・保存の容量)・本段の直前の照合が合格・受け渡しの関門を呼ばない "
          f"({len(gate_calls)} 回)・正式の判定 h370_judge.py を呼ばない ({len(judge_calls)} 回)・update_conditions_met は偽・診断の印",
          rc == 0 and pc["ok"] and pc["phase"] == "ns-diag" and pc["update_conditions_met"] is False and set(pc["versions"]) == set(HR.VERSION_FILES)
          and all(pc["versions"].values()) and pc["items"]["D1_capacity"]["ok"] and json.loads((run / "H370_MAIN_START.json").read_text())["ok"]
          and not gate_calls and not judge_calls and end["update_conditions_met"] is False and (run / HR.DIAG_MARK).is_file()
          and json.loads((run / HR.DIAG_MARK).read_text())["update_conditions_met"] is False, str(pc["failures"])[:300])
    check("(a) Euler の参照は run_0950 の res_60000 (sha) を未適格のまま診断専用に: prepare_ns・停止規則・投入前の検査の euler_ref が DIAG_REF の固定、"
          "記録の状態は「未適格 (予備の資格 FAIL)」",
          pc["euler_ref"]["res"] == "res_60000.h5" and pc["euler_ref"]["sha256"] == dref["euler_ref"]["sha256"]
          == HR.sha256_file(c / HR.EU_RUN / "res_60000.h5") and pc["items"]["N5_prepare"]["euler_ref"] == pc["euler_ref"]
          and dref["qualified"] is False and "未適格" in dref["kind"] and dref["prequalify_record"]["cross_check"] == "一致"
          and set(dref["sens_refs"]) == {"A", "B"} and all(len(v["sha256"]) == 64 for v in dref["sens_refs"].values()))
    check("(b) 正式の記録 (予備の資格 h370_euler_prequalify.json) はバイトで同じ・正式の固定 h370_euler_ref_pin.json は作らない・終了の記録の照合も "
          "「書き換えなし」・診断の参照の固定は別のファイル (h370_diag_euler_ref.json)",
          (c / HR.QUAL).read_bytes() == qual_before and not (c / HR.PIN).exists() and end["formal_records"]["unchanged"] is True
          and (c / HR.DIAG_REF).is_file() and HR.DIAG_REF.endswith("h370_diag_euler_ref.json"))
    main_files = [f"res_{k}.h5" for k in seen["steps"]] + [f"res_wall_3_{k}.h5" for k in seen["steps"]] + [f"res_outlet_2_{k}.h5" for k in seen["steps"]]
    check(f"(f) 本段の全出力 {len(seen['steps'])} 枚 (場・壁・出口の境界の出力・xmf) が診断の経路の後も残り、sha の一覧 ({shal['n']} ファイル) に入る・"
          "停止規則の帯・記録・stage の記録も一覧に入る",
          all((run / n).is_file() for n in main_files) and all(n in shal["files"] for n in main_files)
          and "stop_rule_band/delta_r_edge.csv" in shal["files"] and "stop_rule.json" in shal["files"]
          and shal["files"][f"res_{S_DIAG}.h5"] == HR.sha256_file(run / f"res_{S_DIAG}.h5"))
    w = end.get("window") or {}
    check(f"終了の記録: REACHED・停止の時刻 {end.get('stop_step')}・判定の窓 17 出力 (欠けなし)・本段の出力の完全性 (欠けなし)・到達と中断の履歴・"
          f"check_convergence --segment の VERDICT (記録)",
          end["end_state"] == "REACHED" and end["stop_step"] == S_DIAG and w.get("n") == 17 and w["steps"][0] == S_DIAG - 40000
          and not w["missing_res"] and not w["missing_wall"] and w["kind"].startswith("停止の時刻の窓") and end["main_outputs"]["complete"]
          and end["main_outputs"]["gaps"] == [] and end["stop_rule"]["hold"]["breaks"] and "rc" in end["convergence (記録)"]
          and (run / HR.DIAG_CONV).is_file(), str(end.get("main_outputs")))
    # (a) 正式の判定を診断の run に当てても関門を呼ばない
    with Patch((HJ, "HOOKS", {"handoff_gate": gate_boom(gate_calls)}), (DL, "handoff_gate", gate_boom(gate_calls))):
        rcj, Rj = judge(c, ["--no-ab"])
    check("(a) 正式の判定 (h370_judge.py の既定) を診断の run に当てる → 判定不能 (rc 2)・関門を呼ばない・update_conditions_met 偽",
          rcj == 2 and not gate_calls and Rj["verdict"].startswith("判定不能") and "診断専用" in Rj["verdict"] and Rj["update_conditions_met"] is False)
    # (c) 正式の --phase ns は変わらない: 診断の記録があっても資格 FAIL で止まる (read_pin の資格の検査は迂回されない)
    with Patch((HR, "HERE", c)):
        e1 = raises(HR.read_pin, SystemExit)
        e2 = raises(lambda: HR.phase_ns(args(phase="ns"), RA, {}), SystemExit)
    check("(c) 正式の --phase ns: 診断の参照の固定があっても read_pin・phase_ns は予備の資格 FAIL で止まる (迂回は診断の経路の中だけ)",
          e1 is not None and e2 is not None)
    # (c) 同じ診断の経路をもう一度 → 既存の run で止める (消さない)
    with Patch(*base, (RA, "run_staged_ns", make_fake_staged_ns_diag())):
        e3 = raises(lambda: HR.phase_ns_diag(args(phase="ns-diag"), RA, {}), SystemExit)
    check("診断の経路: run_0951 が既にある → 止める (消さない)", e3 is not None and (run / f"res_{S_DIAG}.h5").is_file())
    # CENSORED・実行の失敗
    for name, kw, want in (("CENSORED", {"status": "CENSORED", "last": HR.NS_CAP}, "CENSORED"), ("失敗", {"raise_exc": True}, "実行の失敗")):
        c2 = diag_fake_case()
        with Patch(*base[:2], (HR, "HERE", c2), *base[3:], (RA, "run_staged_ns", make_fake_staged_ns_diag(**kw))):
            rc2 = HR.phase_ns_diag(args(phase="ns-diag"), RA, {})
        e2_ = json.loads((c2 / HR.NS_RUN / HR.DIAG_END).read_text())
        rr = (c2 / HR.NS_RUN / "RUN_RC").read_text().strip()
        ok_ = (e2_["end_state"] == want and (rr == "EXC" if want == "実行の失敗" else rr == "0")
               and (e2_["stop_step"] is None if want == "実行の失敗" else (e2_["stop_step"] == HR.NS_CAP and e2_["window"]["kind"].startswith("打切りの時点")
                                                                        and e2_["window"]["steps"][0] == HR.NS_CAP - 40000 and e2_["main_outputs"]["n"] == 80)))
        check(f"終了の状態の区別: {name} → {e2_['end_state']} (RUN_RC {rr}{'、打切りの時点の診断の窓 = 末尾の 17 出力' if want == 'CENSORED' else ''})・"
              "関門を呼ばない", rc2 == 0 and ok_ and not gate_calls, str({k: e2_.get(k) for k in ("end_state", "stop_step", "exception")})[:300])
    # (c) 診断の経路の入口の検査: 資格ありの記録・正式の固定・記録の無い予備の資格・参照の sha の違い → 止める
    for name, mk in (("予備の資格の記録が資格あり (prequalified 真)", lambda: diag_fake_case(qual_prequalified=True)),
                     ("正式の固定 PIN がある", lambda: diag_fake_case(pin=True))):
        c3 = mk()
        with Patch((HR, "HERE", c3)):
            e = raises(HR.diag_euler_ref, SystemExit)
        check(f"(c) 診断の経路の入口: {name} → 止める (正式の --phase ns を使う・記録の食い違い)・診断の固定を書かない",
              e is not None and not (c3 / HR.DIAG_REF).exists())
    c3 = diag_fake_case()
    (c3 / HR.QUAL).unlink()
    with Patch((HR, "HERE", c3)):
        e = raises(HR.diag_euler_ref, SystemExit)
    check("(c) 予備の資格の記録が無い → 止める", e is not None)
    c3 = diag_fake_case()
    with Patch((HR, "HERE", c3)):
        r1 = HR.diag_euler_ref()
        (c3 / ARM_RUNS["B"] / "res_40000.h5").write_bytes(b"changed")
        e = raises(HR.diag_euler_ref, SystemExit)
    check("(c) 診断の参照の固定の後に腕 B の res_40000 の中身が変わる → 止める (固定を書き換えない)",
          r1["spec"]["res"] == "res_60000.h5" and e is not None
          and json.loads((c3 / HR.DIAG_REF).read_text())["sens_refs"]["B"]["sha256"] == r1["record"]["sens_refs"]["B"]["sha256"])
    # --- (d) 保存の容量 ---
    c4 = diag_fake_case()
    calls4 = []
    staged4 = lambda *a, **k: calls4.append("forge") or 0  # noqa: E731
    with Patch(*base[:2], (HR, "HERE", c4), *base[3:5], (HR, "disk_free", lambda: 1 * 2 ** 30), *base[6:], (RA, "run_staged_ns", staged4)):
        rc4 = HR.phase_ns_diag(args(phase="ns-diag"), RA, {})
    pc4 = json.loads((c4 / HR.NS_RUN / "H370_PRECHECK.json").read_text())
    d1 = pc4["items"]["D1_capacity"]
    per = max(d1["nozzle_h5_bytes"].values())
    want_need = int(np.ceil((HR.DIAG_N_OUT + HR.DIAG_EXTRA_OUT) * (per + per * HR.DIAG_WALL_FRAC) * HR.DIAG_MARGIN)) + HR.DISK_MIN
    check(f"(d) 保存の容量の不足 (空き 1 GiB < 要る {d1['need_gib']} GiB) → 前提の不成立 (rc {rc4})・PRECHECK_FAILED・forge を回さない; 見積もり = "
          f"(80 + 2) × (nozzle.h5 + 1 %) × 1.1 + 2 GiB",
          rc4 == 2 and not calls4 and (c4 / HR.NS_RUN / "PRECHECK_FAILED").is_file() and not d1["ok"] and d1["need_bytes"] == want_need
          and d1["n_outputs"] == 82 and HR.NS_PRELAUNCH in d1["nozzle_h5_bytes"] and any("保存の容量" in f for f in pc4["failures"]))
    with Patch((HR, "disk_free", lambda: want_need)):
        ok_eq = HR.capacity_check(c4 / HR.NS_RUN, c4 / HR.NS_PRELAUNCH)["ok"]
    with Patch((HR, "disk_free", lambda: want_need - 1)):
        ng_eq = HR.capacity_check(c4 / HR.NS_RUN, c4 / HR.NS_PRELAUNCH)["ok"]
    check("(d) 容量の境目: 空き = 要る → 合格、1 バイト足りない → 不成立; 投入前の準備が無い → 不成立",
          ok_eq and not ng_eq and not HR.capacity_check(c4 / HR.NS_RUN, c4 / "no_prelaunch")["ok"])
    c5 = diag_fake_case()
    q5 = (c5 / HR.QUAL).read_bytes()
    calls5 = []
    with Patch(*base[:2], (HR, "HERE", c5), *base[3:], (RA, "run_staged_ns", lambda *a, **k: calls5.append("forge") or 0)):
        rc5 = HR.phase_ns_diag(args(phase="ns-diag", dry_run="d1"), RA, {})
    pc5 = json.loads((c5 / "_h370_dry_d1_nsdiag" / "H370_PRECHECK.json").read_text())
    check("--dry-run (診断の経路): 投入前の検査の合格 (rc 0)・forge を回さない・run_0951 を作らない・診断の印は dry の準備に・正式の記録は同じ",
          rc5 == 0 and pc5["ok"] and not calls5 and not (c5 / HR.NS_RUN).exists() and (c5 / "_h370_dry_d1_nsdiag" / HR.DIAG_MARK).is_file()
          and (c5 / HR.QUAL).read_bytes() == q5 and not (c5 / HR.PIN).exists(), str(pc5.get("failures"))[:200])
    # --- (f) 消す処理 ---
    dd = TMP / "delete_after"
    dd.mkdir()
    names = [f"res_{k}.h5" for k in (2500, 5000, 7500, 10000)] + [f"res_{k}.xmf" for k in (5000, 7500)] + \
            [f"res_wall_3_{k}.h5" for k in (5000, 7500)] + ["res_nan_7500.h5", "stop_rule.json", "residual_history.csv"]
    for n in names:
        (dd / n).write_text("x")
    gone = SR.delete_after(dd, 5000)
    left = sorted(p.name for p in dd.iterdir())
    check(f"(f) stop_rule.delete_after(s = 5,000) は s より後の出力 (場・xmf・境界) だけを消す (消した {gone})、s 以前・res_nan_*・記録は残す",
          sorted(gone) == sorted(["res_7500.h5", "res_10000.h5", "res_7500.xmf", "res_wall_3_7500.h5"])
          and set(left) == set(names) - set(gone))
    srcs = {f"h370_run.{f.__name__}": inspect.getsource(f) for f in (HR.phase_ns_diag, HR.diag_end_record, HR.diag_euler_ref, HR.capacity_check,
                                                                   HR.files_sha256, HR.write_diag_mark)}
    srcs.update({f"h370_judge.{f.__name__}": inspect.getsource(f) for f in (HJ.diag_sens, HJ.diag_snapshot, HJ.diag_pair, HJ.classify_diag)})
    bad = [n for n, s_ in srcs.items() if any(w_ in s_ for w_ in ("unlink(", "rmtree(", "os.remove(", "delete_after("))]
    rs = inspect.getsource(RA.run_staged_ns)
    tail_main = rs[rs.index("done = [_stage("):]
    un = [ln.strip() for ln in tail_main.splitlines() if ".unlink(" in ln]
    check(f"(f) 診断の関数に消す処理が無い ({bad or 'なし'})・run_staged_ns の本段の部分で消すのは前の段の残差の履歴の写し (residual_history.csv) だけ "
          f"({un}; 前の段の res_* を消すのは本段の前の _stage の中)",
          not bad and all("residual_history.csv" in ln for ln in un) and "for f in run_dir.glob(\"res_*\"):" in rs[:rs.index("done = [_stage(")])
    # --- (e) 感度の分類 (単体) ---
    req = SR.required_window_steps(S_DIAG, 20000, 2500)
    P = lambda Ds, pfA=None, pfB=None, part=None: [{"step": k, "D_mm": (d if (pfA is None or pfA[i]) and (pfB is None or pfB[i]) else None),  # noqa: E731
                                                     "pf_A": True if pfA is None else pfA[i], "pf_B": True if pfB is None else pfB[i],
                                                     **({"D_partial_mm": part[i]} if part is not None else {})}
                                                    for i, (k, d) in enumerate(zip(req, Ds))]
    mixed = [0.008, 0.012] * 8 + [0.008]
    one_pf = [True] * 17
    one_pf[5] = False
    both = [True] * 17
    both[3] = False
    for name, per_, want, code in (
            ("全 17 出力で D 0.003 mm → 結果 A", P([0.003] * 17), "結果 A", "A"),
            ("ちょうど 0.01 mm → 結果 A (≤)", P([0.01] * 17), "結果 A", "A"),
            ("D が窓の中で 0.008/0.012 をまたぐ → 最大 0.012 で結果 B (時刻ごとに判定不能にしない; classify_ab なら U_time)", P(mixed), "結果 B", "B_D"),
            ("1 出力だけ 0.0101 mm → 結果 B (最大で決まる)", P([0.001] * 16 + [0.0101]), "結果 B", "B_D"),
            ("1 出力で入れ替えにより合否が変わる (A 否・B 合、D は小さい) → 結果 B", P([0.001] * 17, pfA=one_pf), "結果 B", "B_pf"),
            ("両腕とも合否 否の出力 (有効な列の D 0.002) → 判定不能", P([0.001] * 17, pfA=both, pfB=both, part=[0.002] * 17), "判定不能", "U_both_invalid"),
            ("両腕とも合否 否の出力でも有効な列の D 0.02 > 0.01 → 結果 B (D の下限が超える)", P([0.001] * 17, pfA=both, pfB=both, part=[0.02] * 17),
             "結果 B", "B_D"),
            ("17 出力のうち 1 つが無い → 判定不能 (縮めない)", P([0.001] * 17)[:-1] + [{"step": req[-1], "error_missing": "無い"}], "判定不能", "U_incomplete"),
            ("16 出力しか渡されない → 判定不能", P([0.001] * 17)[1:], "判定不能", "U_incomplete"),
            ("抽出の例外 (腕 B の 1 出力) → 判定不能", P([0.001] * 16) + [{"step": req[-1], "error_B": "RuntimeError"}], "判定不能", "U_error")):
        cl = HJ.classify_diag(per_, req)
        check(f"(e) 分類: {name}", cl["result"] == want and cl["code"] == code, f"{cl['result']} {cl['code']}: {cl['why'][:120]}")
    check("(e) 分類: 必要な出力の集合が無い → 判定不能", HJ.classify_diag(P([0.001] * 17), [])["code"] == "U_required")
    x3, F3 = np.array([10.0, 11.0, 12.0]), np.ones(3, dtype=bool)
    dB = {"delta_E": [1e-3] * 3, "delta_target": [2e-3] * 3, "held": [0] * 3, "missing": [], "y_b": [0.1] * 3, "y_star": [0.08] * 3}
    dA = {"delta_E": [1e-3, float("nan"), 1e-3 + 0.02 / 76.6], "delta_target": [2e-3] * 3, "held": [0] * 3, "missing": [1], "y_b": [0.1, 0.1, 0.11],
          "y_star": [0.08] * 3}
    r_ = HJ.diag_pair(S_DIAG, dA, dB, 76.6, x3, F3)
    r2 = HJ.diag_pair(S_DIAG, dA, dA, 76.6, x3, F3)
    check(f"(e) diag_pair: 腕 A の中ほどの欠測 → 合否 否・D なし・有効な列の差 {r_.get('D_partial_mm', 0):.4f} mm と帯の位置 y_b の差 "
          f"{r_.get('band_y_b_max_diff_mm', 0):.3f} mm を記録 / 両腕とも同じ欠測 → D なし・有効な列の差 0",
          r_["pf_A"] is False and r_["pf_B"] is True and r_["D_mm"] is None and abs(r_["D_partial_mm"] - 0.02) < 1e-9
          and abs(r_["band_y_b_max_diff_mm"] - 0.766) < 1e-9 and r2["D_mm"] is None and r2["D_partial_mm"] == 0.0
          and HJ.classify_diag([r2], [S_DIAG])["code"] == "U_both_invalid" and HJ.classify_diag([r_], [S_DIAG])["code"] == "B_pf")
    # --- (e) 感度の判定の全経路 (偽物の抽出) ---
    dbase = write_diag_case(TMP / "diag_base")
    hooks = lambda snap: {"handoff_gate": gate_boom(gate_calls), "diag_snapshot": snap}  # noqa: E731

    def run_diag(name, snap, prep=None, extra=()):
        case = TMP / f"diag_{name}"
        shutil.copytree(dbase, case, symlinks=True)
        if prep:
            prep(case)
        with Patch((HJ, "HOOKS", hooks(snap)), (DL, "handoff_gate", gate_boom(gate_calls))):
            return diag_judge(case, extra)
    calls = []
    rc, R = run_diag("A", fake_diag_snap(lambda arm, k, x: 0.003 * np.ones(len(x)) if arm == "B" else 0.0 * x, calls=calls))
    check(f"(e) 全経路 (偽物の抽出): 両腕の差 0.003 mm → 結果 A (rc {rc})・17 出力 × (A・B・補助) = 51 抽出・比べる列は抽出の前に固定・update_conditions_met 偽・"
          "関門を呼ばない・配列の npz",
          rc == 0 and R["result"] == "結果 A" and R["code"] == "A" and len(calls) == 51 and sorted({k for _, k in calls}) == req
          and abs(R["D_mm"]["max"] - 0.003) < 1e-9 and R["columns"]["n"] == R["columns"]["n_use"] and R["update_conditions_met"] is False
          and not gate_calls and (TMP / "diag_A" / HJ.DIAG_NPZ).is_file() and R["window_kind"].startswith("停止の時刻の窓")
          and all("max_abs_vs_A_mm" in a_ for a_ in R["aux_run_0950"]["per_output"]), str(R.get("why"))[:300])
    npz = np.load(TMP / "diag_A" / HJ.DIAG_NPZ)
    check("(e) 配列: 出力 17 × 使用域の列の δ_E・δ_target・held・y_b・y* (腕 A・B・補助)", npz["A_delta_E"].shape == (17, len(npz["x_use"]))
          and {"B_y_b", "aux_delta_target", "A_held", "F", "steps"} <= set(npz.files))
    rc, R = run_diag("mixed", fake_diag_snap(lambda arm, k, x: (0.012 if (k // 2500) % 2 else 0.008) * np.ones(len(x)) if arm == "B" else 0.0 * x))
    check(f"(e) 全経路: D が窓の中で 0.008/0.012 をまたぐ → 最大 {R['D_mm']['max']:.3f} mm で結果 B ({R['code']}; 時刻ごとの判定不能にしない)",
          R["result"] == "結果 B" and R["code"] == "B_D" and abs(R["D_mm"]["max"] - 0.012) < 1e-9 and abs(R["D_mm"]["min"] - 0.008) < 1e-9
          and 0 < R["D_mm"]["n_over_tol"] < 17, str(R.get("why"))[:200])

    def held_b(arm, k, d):
        if arm == "B" and k == req[4]:
            d["held"] = list(d["held"])
            d["held"][len(d["held"]) // 2] = 1
            d["n_held"] = 1
        return d
    rc, R = run_diag("pf", fake_diag_snap(mutate=held_b))
    check(f"(e) 全経路: 1 出力で腕 B だけ比べる列の中に held → 合否が変わる → 結果 B ({R['code']})", R["result"] == "結果 B" and R["code"] == "B_pf",
          str(R.get("why"))[:200])

    def exc_a(arm, k, d):
        if arm == "A" and k == req[7]:
            raise RuntimeError("extract_and_merge: 単調性ガードで λ が上限を超えた (試験)")
        return d
    rc, R = run_diag("exc", fake_diag_snap(mutate=exc_a))
    check(f"(e) 全経路: 抽出の例外 (腕 A の 1 出力) → 判定不能 ({R['code']}、rc {rc})", rc == 2 and R["result"] == "判定不能" and R["code"] == "U_error")
    for name, prep, want in (
            ("窓の 1 出力の res が無い", lambda cs: (cs / HJ.NS_RUN / f"res_{req[8]}.h5").unlink(), "が無い"),
            ("窓の 1 出力の中身が sha の一覧と違う", lambda cs: shutil.copy(cs / HJ.NS_RUN / f"res_{req[0]}.h5", cs / HJ.NS_RUN / f"res_{req[1]}.h5"), "sha"),
            ("腕 B の参照の sha が固定と違う", lambda cs: (cs / ARM_RUNS["B"] / "res_40000.h5").write_bytes(b"changed"), "参照 B"),
            ("抽出器の版が投入の時点と違う", lambda cs: (cs / HJ.NS_RUN / "H370_PRECHECK.json").write_text(json.dumps(
                {"ok": True, "phase": "ns-diag", "versions": {**HJ.module_versions(), "metrics/deltastar.py": "0" * 64}})), "版"),
            ("診断の印が無い (正式の run)", lambda cs: (cs / HJ.NS_RUN / HJ.DIAG_MARK).unlink(), "診断の印"),
            ("停止の状態 DIVERGED", lambda cs: SR.write_state(cs / HJ.NS_RUN, {"status": "DIVERGED", "stop_step": None}), "DIVERGED"),
            ("抽出の条件が #36e の固定と違う (ω 0.7)", lambda cs: (lambda p_, r_: p_.write_text(json.dumps({**r_, "delta_series": {
                **r_["delta_series"], "extract": {**r_["delta_series"]["extract"], "omega": 0.7}}})))(
                cs / HJ.NS_RUN / SR.RECORD, json.loads((cs / HJ.NS_RUN / SR.RECORD).read_text())), "抽出の条件")):
        calls = []
        n_u = len(list(TMP.glob("diag_u*")))
        rc, R = run_diag(f"u{n_u}", fake_diag_snap(calls=calls), prep=prep)
        check(f"(e) 全経路: {name} → 判定不能 (rc {rc}、{R.get('code')})・抽出しない ({len(calls)} 回)", rc == 2 and R["result"] == "判定不能"
              and R["code"] == "U_input" and not calls and any(want in u for u in R["undecidable"]) and R["update_conditions_met"] is False,
              str(R.get("undecidable"))[:300])
    ccase = write_diag_case(TMP / "diag_cens_base", status="CENSORED", s_stop=200000)
    with Patch((HJ, "HOOKS", hooks(fake_diag_snap())), (DL, "handoff_gate", gate_boom(gate_calls))):
        rc, R = diag_judge(ccase)
    check(f"(e) 全経路: CENSORED (最後の出力 200,000) → 打切りの時点の診断の窓 (末尾の 17 出力 160,000〜200,000) で判定 ({R['result']})",
          R["result"] == "結果 A" and R["required"] == list(range(160000, 200001, 2500)) and R["window_kind"].startswith("打切りの時点"))
    # --- (e) 本物の抽出 (diag_snapshot、プロセス並列): stop_rule.delta_snapshot と同じ値・結果 A・B ---
    ns = dbase / HJ.NS_RUN
    rec = json.loads((ns / SR.RECORD).read_text())
    ext = {k: v for k, v in rec["delta_series"]["extract"].items() if k != "source"}
    specB = {"run": str(dbase / ARM_RUNS["B"]), "res": "res_40000.h5", "sha256": HJ.sha256_file(dbase / ARM_RUNS["B"] / "res_40000.h5")}
    d_sr = SR.delta_snapshot(ns, S_DIAG, specB, ext, rec["delta_series"]["x"])
    d_dg = HJ.diag_snapshot(ns, S_DIAG, specB, ext, rec["delta_series"]["x"])
    same = all(np.array_equal(np.asarray(d_sr[k], dtype=float), np.asarray(d_dg[k], dtype=float), equal_nan=True) for k in ("delta_E", "delta_target", "held")) \
        and all(d_sr[k] == d_dg[k] for k in ("missing", "n_missing", "n_held", "lam", "res_sha256"))
    check(f"(e) diag_snapshot は stop_rule.delta_snapshot と δ_E・δ_target・held・missing・λ・出力の sha が同じ値で、帯の位置 y_b・y* を足す "
          f"(y_b の有限 {int(np.isfinite(d_dg['y_b']).sum())}/{len(d_dg['y_b'])})", same and np.isfinite(d_dg["y_b"]).all() and np.isfinite(d_dg["y_star"]).all())
    with Patch((HJ, "HOOKS", {"handoff_gate": gate_boom(gate_calls)}), (DL, "handoff_gate", gate_boom(gate_calls))):
        rc, R = diag_judge(dbase)
        rB_case = write_diag_case(TMP / "diag_real_B", shape_b=0.01)
        rcB, RB = diag_judge(rB_case)
    px = R.get("extractions", {}).get("parallel", {})
    check(f"(e) 本物の抽出 (diag_snapshot、{px.get('mode')} {px.get('jobs')}、{px.get('seconds')} s): 腕 A・B の Euler が同じ場 → D 0 mm → 結果 A; "
          f"腕 B の形を変える → D {RB.get('D_mm', {}).get('max', float('nan')):.4f} mm > 0.01 → 結果 B; 補助 (run_0950) の差も記録・関門を呼ばない",
          rc == 0 and R["result"] == "結果 A" and R["D_mm"]["max"] == 0.0 and px.get("mode", "").startswith("process") and rcB == 0
          and RB["result"] == "結果 B" and RB["code"] == "B_D" and RB["D_mm"]["max"] > HJ.TOL_MM and not gate_calls
          and all(a_.get("max_abs_vs_A_mm") == 0.0 for a_ in R["aux_run_0950"]["per_output"])
          and all(r_.get("band_y_b_max_diff_mm") is not None for r_ in RB["per_output"]), str((R.get("why"), RB.get("why")))[:300])


# ====================================================================================== M. 帯を共通にした A/B (§5.1 #36n)
X_OUT = (-1.5, -0.5)                                                    # 偽物の抽出の使用域の外の列 (縮流部; #36l の記録に無い列)


def fake_fixed_snap(effect=None, mutate=None, calls=None, a_offset=None):
    """偽物の帯を共通にした抽出 (速い; diag_snapshot_fixed の形): δ_E = δ_target = 1e-3 r_t + effect(腕, 出力, x) [mm] / MMPR (腕 A は
    fake_diag_snap と同じ値 = #36l の偽物の記録を再現) + a_offset(出力) [r_t] (腕 A の再現の不成立を作る)。帯を決める列 = 使用域の外 X_OUT + 使用域。
    帯は band の表を x の完全一致で引き、表に無い列は fill_own なら自分の帯 (選び直した帯 = 0.1、腕 B は 0.11)、そうでなければ例外 (本物と同じ)。
    calls に (腕, 出力, band) を足す。mutate(腕, 出力, d) で held・例外を作る。"""
    def g(run, step, spec, ext, x_use, band):
        x = np.asarray(x_use, dtype=float)
        arm = "A" if Path(spec["run"]).name == ARM_RUNS["A"] else "B"
        if calls is not None:
            calls.append((arm, int(step), band))
        xa = np.concatenate([np.asarray(X_OUT, dtype=float), x])
        lut = dict(zip(np.asarray(band["x"], dtype=float).tolist(), np.asarray(band["y_b"], dtype=float).tolist()))
        nat = (0.11 if arm == "B" else 0.1) + 0.0 * xa
        hit = np.array([v in lut for v in xa.tolist()])
        if not hit.all() and not band.get("fill_own"):
            raise ValueError("固定する帯の表に無い列 (偽物)")
        yb = np.array([lut[v] if h else nat[i] for i, (v, h) in enumerate(zip(xa.tolist(), hit))])
        e = 1e-3 + (np.asarray(effect(arm, int(step), x), dtype=float) / MMPR if effect else 0.0) + 0.0 * x
        if arm == "A" and a_offset:
            e = e + a_offset(int(step))
        n = len(X_OUT)
        d = {"delta_E": e.tolist(), "delta_target": e.tolist(), "held": [0] * len(x), "missing": [], "n_missing": 0, "n_held": 0, "lam": 1.0,
             "res_sha256": None, "seconds": 0.0, "y_b": yb[n:].tolist(), "y_star": (0.08 + 0.0 * x).tolist(), "y_b_natural": nat[n:].tolist(),
             "delta_r_raw": (0.9 * e).tolist(), "band_all": {"x": xa.tolist(), "y_b": yb.tolist(), "y_b_natural": nat.tolist()},
             "band_fill": {"n_table": int(hit.sum()), "n_own": int((~hit).sum()), "x_own_range": None, "edge_table_matches_fixed": True}}
        if mutate:
            d = mutate(arm, int(step), d) or d
        return d
    return g


def fixed_judge(case: Path, extra=()) -> tuple:
    rc = HJ.main(["--diag-sens-fixed-band", "--case-dir", str(case), "--design", str(ROOT / "design"), "--tools", str(TOOLS), *extra])
    return rc, HJ.jload(case / HJ.FIXED_OUT)


def test_fixed_band():
    print("== M. 帯を共通にした A/B (--diag-sens-fixed-band; §5.1 #36n)", flush=True)
    import inspect
    gate_calls = []
    req = SR.required_window_steps(S_DIAG, 20000, 2500)
    # --- 分類 (単体) ---
    P = lambda Ds, err=None, bad=None: [{"step": k, "D_mm": (None if bad and i in bad else d), "pf_A": True, "pf_B": not (bad and i in bad),  # noqa: E731
                                         **({"error_B": "RuntimeError"} if err and i in err else {}),
                                         **({"D_partial_mm": 0.02, "x_at_D_valid": 1.0} if bad and i in bad else {}), "x_at_D": 2.0}
                                        for i, (k, d) in enumerate(zip(req, Ds))]
    for name, per_, repro, want, code in (
            ("全 17 出力で D_fixed 0.003 mm → 結果 A", P([0.003] * 17), True, "結果 A", "A"),
            ("ちょうど 0.01 mm → 結果 A (≤)", P([0.01] * 17), True, "結果 A", "A"),
            ("1 出力だけ 0.0101 mm → 結果 B", P([0.001] * 16 + [0.0101]), True, "結果 B", "B_D"),
            ("腕 A の再現の不成立 (D は小さい) → 判定不能", P([0.001] * 17), False, "判定不能", "U_repro"),
            ("1 出力で腕 B の F の中に held (有効な列の差 0.02 > 0.01 でも) → 判定不能 (登録: 固定の帯で欠測・非有限・抽出の不成立)",
             P([0.001] * 17, bad={4}), True, "判定不能", "U_invalid"),
            ("抽出の例外 (腕 B の 1 出力) → 判定不能", P([0.001] * 17, err={7}), True, "判定不能", "U_error"),
            ("17 出力のうち 1 つが無い → 判定不能 (縮めない)", P([0.001] * 17)[1:], True, "判定不能", "U_incomplete")):
        cl = HJ.classify_fixed(per_, req, repro)
        check(f"M 分類: {name}", cl["result"] == want and cl["code"] == code, f"{cl['result']} {cl['code']}: {cl['why'][:120]}")
    check("M 分類: 必要な出力の集合が無い → 判定不能", HJ.classify_fixed(P([0.001] * 17), [], True)["code"] == "U_required")
    xs = np.array([0.5, 2.0, 4.0, 7.0, 20.0, 60.0, 94.0, 95.3])
    Dk = np.array([[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], [0.0, 0.9, 0.0, 0.0, np.nan, 0.0, 0.0, 0.85]])
    sm = HJ.section_max(xs, Dk, np.array([True] * 7 + [False]), [100, 200])
    check("M 区間ごとの最大: 区間 [lo, hi)・出口までの最後の区間・比べる列 F の外と非有限を除く・最大の出力と x",
          [r_["max_mm"] for r_ in sm] == [0.1, 0.9, 0.3, 0.4, 0.5, 0.6, 0.7] and sm[1]["step"] == 200 and sm[1]["x"] == 2.0
          and sm[6]["x_range"] == [94.0, 94.0] and len(sm) == len(HJ.FIXED_SECTIONS) == 7, str(sm))
    # --- 全経路 (偽物の抽出): #36l の偽物の記録 (腕 B − A = 0.02 mm) → 帯を共通にした A/B ---
    mbase = write_diag_case(TMP / "fixed_base")
    with Patch((HJ, "HOOKS", {"handoff_gate": gate_boom(gate_calls), "diag_snapshot": fake_diag_snap(
            lambda arm, k, x: 0.02 * np.ones(len(x)) if arm == "B" else 0.0 * x)}), (DL, "handoff_gate", gate_boom(gate_calls))):
        rc36, R36 = diag_judge(mbase)
    check(f"M 前提: #36l の偽物の記録 (結果 B、D 0.02 mm、rc {rc36})", rc36 == 0 and R36["result"] == "結果 B" and (mbase / HJ.DIAG_NPZ).is_file())
    z36 = np.load(mbase / HJ.DIAG_NPZ)

    def run_fixed(name, snap, prep=None):
        case = TMP / f"fixed_{name}"
        shutil.copytree(mbase, case, symlinks=True)
        if prep:
            prep(case)
        with Patch((HJ, "HOOKS", {"handoff_gate": gate_boom(gate_calls), "diag_snapshot_fixed": snap}), (DL, "handoff_gate", gate_boom(gate_calls))):
            return fixed_judge(case)
    calls = []
    eff_a = lambda arm, k, x: np.where((x >= 3.0) & (x < 5.0), 0.003, 0.0) if arm == "B" else 0.0 * x  # noqa: E731
    rc, R = run_fixed("A", fake_fixed_snap(eff_a, calls=calls))
    arms = [a_ for a_, _, _ in calls]
    sec = {tuple(s_["x_range"][:1]): s_ for s_ in R.get("sections", [])}
    check(f"M 全経路 (偽物の抽出): 帯を共通にすると差 0.003 mm (x 3〜5 だけ) → 結果 A (rc {rc})・34 抽出 (腕 A 17 → 腕 B 17 の順)・"
          "update_conditions_met 偽・関門を呼ばない・#36l の記録と照合 (npz の sha・比べる列)・腕 A の再現 (ビット同一、最大差 0)",
          rc == 0 and R["result"] == "結果 A" and R["code"] == "A" and arms == ["A"] * 17 + ["B"] * 17 and sorted({k for _, k, _ in calls}) == req
          and abs(R["D_fixed_mm"]["max"] - 0.003) < 1e-12 and R["update_conditions_met"] is False and not gate_calls
          and R["columns"]["same_as_36l"] is True and R["source_36l"]["result"] == "結果 B" and R["repro_A"]["ok"] is True
          and R["repro_A"]["bitwise_all"] is True and R["repro_A"]["max_abs_delta_E_rt"] == 0.0 and R["extractions"]["n"] == 34,
          str(R.get("why"))[:300])
    check("M 区間ごと: D_fixed は 3〜5 だけ 0.003 mm、ほかの区間は 0・#36l の D (帯を選び直した) は全区間 0.02 mm を並べる・帯の移動 (0.01 r_t) を記録",
          len(R["sections"]) == 7 and abs(sec[(3.0,)]["D_fixed"]["max_mm"] - 0.003) < 1e-12 and sec[(1.0,)]["D_fixed"]["max_mm"] == 0.0
          and all(abs(s_["D_36l"]["max_mm"] - 0.02) < 1e-12 for s_ in R["sections"] if s_["n_cols"])
          and all(abs(s_["band_move_B"]["max_mm"] - 0.01 * MMPR) < 1e-9 for s_ in R["sections"] if s_["n_cols"])
          and all(abs(v - 0.02) < 1e-12 for v in R["D_fixed_mm"]["D_36l_per_output"].values()), str(R["sections"])[:300])
    bA = [b_ for a_, _, b_ in calls if a_ == "A"]
    bB = [b_ for a_, _, b_ in calls if a_ == "B"]
    stepsA = [k for a_, k, _ in calls if a_ == "A"]
    okA = all(b_["fill_own"] is True and np.array_equal(np.asarray(b_["x"]), z36["x_use"])
              and np.array_equal(np.asarray(b_["y_b"]), z36["A_y_b"][req.index(k)]) for b_, k in zip(bA, stepsA))
    xa_all = np.concatenate([np.asarray(X_OUT), z36["x_use"]])
    okB = all(b_["fill_own"] is False and np.array_equal(np.asarray(b_["x"]), xa_all)
              and np.array_equal(np.asarray(b_["y_b"]), np.concatenate([[0.1] * len(X_OUT), z36["A_y_b"][0]])) for b_ in bB)
    check("M 帯の固定が両腕に同じ値で渡る: 腕 A には #36l の npz の A_y_b(k) (使用域の列・fill_own; 使用域の外は自分の帯)、腕 B には腕 A が同じ出力で"
          "使った全列の帯 (使用域の外を含む・fill_own なし)、両腕の使った帯 (使用域の列) が全出力で同じ (band_common)",
          okA and okB and len(bA) == len(bB) == 17 and all(r_["band_common"] is True for r_ in R["per_output"]))
    npz = np.load(TMP / "fixed_A" / HJ.FIXED_NPZ)
    check("M 配列: 腕 A・B の δ_E・δ_target・held・使った帯・選び直した帯・y*・生の δ_r (出力 17 × 使用域の列)・D_fixed と #36l の D と帯の移動 "
          "(出力 × 列・列ごとの最大)・帯を決めた全列の帯",
          npz["A_delta_r_raw"].shape == (17, len(npz["x_use"])) and {"B_y_b_natural", "D_fixed_kx_mm", "D_36l_x_mm", "band_move_B_x_mm",
                                                                      "band_all_x", "A_band_all_y_b", "B_band_all_y_b", "A_delta_E"} <= set(npz.files)
          and np.array_equal(npz["A_delta_E"], z36["A_delta_E"]))
    eff_b = lambda arm, k, x: np.where((x >= 1.0) & (x < 3.0), 0.02, 0.0) if arm == "B" else 0.0 * x  # noqa: E731
    rc, R = run_fixed("B", fake_fixed_snap(eff_b))
    check(f"M 全経路: 帯を共通にしても差 0.02 mm (x 1〜3) → 結果 B ({R['code']}、D_fixed {R['D_fixed_mm']['max']:.4f} mm、区間 1〜3 の最大)",
          rc == 0 and R["result"] == "結果 B" and R["code"] == "B_D" and abs(R["D_fixed_mm"]["max"] - 0.02) < 1e-12
          and 1.0 <= R["D_max"]["x"] < 3.0 and R["over_tol_columns"]["D_fixed"]["n"] > 0, str(R.get("why"))[:200])

    def held_b(arm, k, d):
        if arm == "B" and k == req[4]:
            d["held"] = list(d["held"])
            d["held"][len(d["held"]) // 2] = 1
        return d
    rc, R = run_fixed("held", fake_fixed_snap(eff_b, mutate=held_b))
    check(f"M 全経路: 固定の帯で腕 B の 1 出力の比べる列の中に held (差 0.02 mm でも) → 判定不能 ({R['code']}、rc {rc})",
          rc == 2 and R["result"] == "判定不能" and R["code"] == "U_invalid")

    def exc_b(arm, k, d):
        if arm == "B" and k == req[9]:
            raise ValueError("固定する帯の表に無い列 (試験)")
        return d
    rc, R = run_fixed("exc", fake_fixed_snap(mutate=exc_b))
    check(f"M 全経路: 腕 B の 1 出力の抽出の不成立 (帯の表に無い列) → 判定不能 ({R['code']}、rc {rc})", rc == 2 and R["code"] == "U_error")
    calls = []
    rc, R = run_fixed("repro", fake_fixed_snap(calls=calls, a_offset=lambda k: 1e-9 if k == req[3] else 0.0))
    check(f"M 腕 A の再現の不成立 (1 出力で δ_E が #36l と 1e−9 r_t 違う) → 判定不能 ({R['code']}、rc {rc})・最大差を記録・腕 B は回さない "
          f"({sum(1 for a_, _, _ in calls if a_ == 'B')} 回)",
          rc == 2 and R["result"] == "判定不能" and R["code"] == "U_repro" and abs(R["repro_A"]["max_abs_delta_E_rt"] - 1e-9) < 1e-15
          and R["repro_A"]["ok"] is False and not any(a_ == "B" for a_, _, _ in calls) and R["update_conditions_met"] is False)
    rc, R = run_fixed("repro_tol", fake_fixed_snap(a_offset=lambda k: 5e-13 if k == req[3] else 0.0))
    check(f"M 腕 A の再現: 丸め程度の差 (5e−13 r_t ≤ 許容 {HJ.REPRO_TOL_RT}) は再現とする (ビット同一ではないと記録) → {R['result']}",
          rc == 0 and R["repro_A"]["ok"] is True and R["repro_A"]["bitwise_all"] is False and R["result"] == "結果 A")

    def a_nan(arm, k, d):
        if arm == "A" and k == req[2]:
            d["delta_E"] = list(d["delta_E"])
            d["delta_E"][3] = float("nan")
        return d
    rc, R = run_fixed("repro_nan", fake_fixed_snap(mutate=a_nan))
    check(f"M 腕 A の再現: 固定の帯で腕 A の 1 列が非有限 (#36l は有限) → 再現の不成立で判定不能 ({R['code']})",
          rc == 2 and R["code"] == "U_repro" and R["repro_A"]["per_output"][2]["finite_pattern_same"] is False)
    # --- 入力の照合 (#36l の記録・#36l と同じ照合) → 判定不能・抽出しない ---

    def set_json(p_, f_):
        p_.write_text(json.dumps(f_(json.loads(p_.read_text()))))

    def npz_change(cs, **kw):
        with np.load(cs / HJ.DIAG_NPZ) as f_:
            a_ = {k: f_[k] for k in f_.files}
        a_.update(kw)
        np.savez_compressed(cs / HJ.DIAG_NPZ, **a_)
        set_json(cs / HJ.DIAG_OUT, lambda r_: {**r_, "arrays": {**r_["arrays"], "sha256": HJ.sha256_file(cs / HJ.DIAG_NPZ)}})
    F_bad = z36["F"].copy()
    F_bad[3] = False
    yb_nan = z36["A_y_b"].copy()
    yb_nan[5, 10] = np.nan
    for name, prep, want in (
            ("#36l の記録が無い", lambda cs: (cs / HJ.DIAG_OUT).unlink(), "#36l の記録"),
            ("#36l の配列の中身が記録の sha と違う", lambda cs: set_json(cs / HJ.DIAG_OUT, lambda r_: {**r_, "arrays": {**r_["arrays"], "sha256": "0" * 64}}),
             "sha256"),
            ("#36l の比べる列 F が座標から決め直した F と違う", lambda cs: npz_change(cs, F=F_bad), "比べる列"),
            ("共通の帯 (#36l の腕 A の y_b) に非有限", lambda cs: npz_change(cs, A_y_b=yb_nan), "非有限"),
            ("#36l の 17 出力が今の窓と違う", lambda cs: set_json(cs / HJ.DIAG_OUT, lambda r_: {**r_, "required": r_["required"][1:]}), "17 出力"),
            ("#36l と同じ照合: 窓の 1 出力の res が無い", lambda cs: (cs / HJ.NS_RUN / f"res_{req[8]}.h5").unlink(), "が無い"),
            ("#36l と同じ照合: 腕 B の参照の sha が固定と違う", lambda cs: (cs / ARM_RUNS["B"] / "res_40000.h5").write_bytes(b"changed"), "参照 B")):
        calls = []
        n_u = len(list(TMP.glob("fixed_u*")))
        rc, R = run_fixed(f"u{n_u}", fake_fixed_snap(calls=calls), prep=prep)
        check(f"M 入力の照合: {name} → 判定不能 (rc {rc}、{R.get('code')})・抽出しない ({len(calls)} 回)", rc == 2 and R["result"] == "判定不能"
              and R["code"] == "U_input" and not calls and any(want in u for u in R["undecidable"]) and R["update_conditions_met"] is False,
              str(R.get("undecidable"))[:300])
    srcs = {f.__name__: inspect.getsource(f) for f in (HJ.diag_inputs, HJ.diag_sens_fixed_band, HJ.diag_snapshot_fixed, HJ.fixed_band,
                                                        HJ.classify_fixed, HJ.section_max)}
    bad = [n for n, s_ in srcs.items() if any(w_ in s_ for w_ in ("unlink(", "rmtree(", "os.remove(", "delete_after(", "handoff_gate", "update_wall"))]
    check(f"M 帯を共通にした A/B の関数に消す処理・関門・壁の更新の呼び出しが無い ({bad or 'なし'})", not bad and not gate_calls)
    # --- 本物の抽出 (diag_snapshot_fixed、プロセス並列) ---
    rb = TMP / "diag_real_B"                                    # L の本物の #36l の記録 (腕 B の形 0.01) があれば使う
    if not (rb / HJ.DIAG_OUT).is_file():
        rb = write_diag_case(TMP / "fixed_real_src", shape_b=0.01)
        with Patch((HJ, "HOOKS", {"handoff_gate": gate_boom(gate_calls)}), (DL, "handoff_gate", gate_boom(gate_calls))):
            diag_judge(rb)
    ns = rb / HJ.NS_RUN
    rec = json.loads((ns / SR.RECORD).read_text())
    ext = {k: v for k, v in rec["delta_series"]["extract"].items() if k != "source"}
    xu = np.asarray(rec["delta_series"]["x"], dtype=float)
    spec = {a_: {"run": str(rb / r_), "res": "res_40000.h5", "sha256": HJ.sha256_file(rb / r_ / "res_40000.h5")} for a_, r_ in ARM_RUNS.items()}
    dA = HJ.diag_snapshot(ns, S_DIAG, spec["A"], ext, xu)
    fA = HJ.diag_snapshot_fixed(ns, S_DIAG, spec["A"], ext, xu, {"x": xu, "y_b": np.asarray(dA["y_b"]), "fill_own": True})
    same = all(np.array_equal(np.asarray(dA[k], dtype=float), np.asarray(fA[k], dtype=float), equal_nan=True) for k in ("delta_E", "delta_target", "held", "y_b"))
    bfA = fA["band_fill"]
    check(f"M 本物の抽出: 腕 A に自分の帯 (使用域の列) を固定 → diag_snapshot と δ_E・δ_target・held・帯がビット同一・使用域の外の {bfA['n_own']} 列 "
          f"(x {bfA['x_own_range']}) は自分の帯・差し替えが抽出の表まで届く",
          same and bfA["n_own"] > 0 and bfA["x_own_range"][1] < 0.0 and bfA["n_table"] == len(xu) and bfA["edge_table_matches_fixed"] is True
          and np.array_equal(np.asarray(fA["y_b_natural"]), np.asarray(fA["y_b"])) and np.isfinite(fA["delta_r_raw"]).all())
    dB = HJ.diag_snapshot(ns, S_DIAG, spec["B"], ext, xu)
    fB = HJ.diag_snapshot_fixed(ns, S_DIAG, spec["B"], ext, xu, {"x": fA["band_all"]["x"], "y_b": fA["band_all"]["y_b"], "fill_own": False})
    mv = float(np.nanmax(np.abs(np.asarray(fB["y_b_natural"]) - np.asarray(fA["y_b"])))) * MMPR
    check(f"M 本物の抽出: 腕 B に腕 A の全列の帯 → 使った帯 (使用域) が腕 A と同じ・選び直した帯は固定しない抽出の帯と同じ (帯の移動 {mv:.3f} mm)・"
          "表に無い列があれば例外 (fill_own なし)",
          np.array_equal(np.asarray(fB["y_b"]), np.asarray(fA["y_b"])) and fB["band_fill"]["n_own"] == 0
          and np.array_equal(np.asarray(fB["y_b_natural"]), np.asarray(dB["y_b"]), equal_nan=True) and mv > 0.0
          and raises(lambda: HJ.diag_snapshot_fixed(ns, S_DIAG, spec["B"], ext, xu, {"x": xu, "y_b": np.asarray(dA["y_b"]), "fill_own": False}),
                     ValueError) is not None)
    import forge_design.metrics.deltastar as DSm
    check("M 本物の抽出: 帯の固定は抽出の後に戻る (deltastar.edge_band_positions・deltastar_loop の抽出器が元の関数)",
          DSm.edge_band_positions is DS.edge_band_positions and DSm.edge_band_positions.__qualname__ == "edge_band_positions"
          and DL.deltastar_from_core_matched_euler.__qualname__ == "deltastar_from_core_matched_euler")
    rcase = TMP / "fixed_real_B"
    shutil.copytree(rb, rcase, symlinks=True)
    with Patch((HJ, "HOOKS", {"handoff_gate": gate_boom(gate_calls)}), (DL, "handoff_gate", gate_boom(gate_calls))):
        rc, R = fixed_judge(rcase)
    px = (R.get("extractions") or {}).get("B") or {}
    dfx = (R.get("D_fixed_mm") or {}).get("max")
    check(f"M 本物の全経路 (腕 B の形 0.01、{px.get('mode')} {px.get('jobs')}): 腕 A の再現 (最大差 {R.get('repro_A', {}).get('max_abs_delta_E_rt')} r_t、"
          f"ビット同一 {R.get('repro_A', {}).get('bitwise_all')})・34 抽出・両腕の帯が同じ・D_fixed {dfx} mm で分類 {R.get('result')} "
          f"(#36l の D {R.get('source_36l', {}).get('D_max', {}).get('mm')})・関門を呼ばない",
          R.get("repro_A", {}).get("ok") is True and R["repro_A"]["bitwise_all"] is True and R["repro_A"]["max_abs_delta_E_rt"] == 0.0
          and R["extractions"]["n"] == 34 and px.get("mode", "").startswith("process") and all(r_["band_common"] for r_ in R["per_output"])
          and R["result"] == ("結果 B" if dfx > HJ.TOL_MM else "結果 A") and rc == 0 and not gate_calls
          and R["repro_A"]["per_output"][0]["band_fill"]["n_own"] > 0 and R["per_output"][0]["band_fill_B"]["n_own"] == 0,
          str(R.get("why"))[:300])


# ====================================================================================== N. 帯の時間方向の扱いの A/B (§5.1 #36p)
TB_K0 = S_DIAG - 40000                                                  # 合成の判定の窓の先頭 k₀ (22,500)


def tb_aux_mm(k, x):
    """#36l の偽物の記録の aux (run_0950) の δ_E の時間変化 [mm] (腕 A がこれを再現する)。"""
    return 0.004 * (k - TB_K0) / 40000.0 + 0.0 * np.asarray(x, dtype=float)


def tb_yb(k):
    """#36l の偽物の記録の aux の帯 (各出力で選び直した帯、出力で動く) [r_t]。"""
    return 0.1 + 0.001 * (k - TB_K0) / 2500.0


def fake_time_snap(effB=None, mutate=None, calls=None, a_offset=None, jump_col=None):
    """偽物の #36p の抽出 (速い; diag_snapshot_fixed の record_ladder の形)。腕は band の fill_own で見分ける (腕 A = 空の表・fill_own、
    腕 B = 表・fill_own なし)。選び直した帯 = tb_yb(k) (全列)、使った帯 = 表から x の完全一致で引く (表に無い列は fill_own なら選び直した帯、
    そうでなければ例外 = 本物と同じ)。δ_E = δ_target = 1e-3 + [tb_aux_mm(k) + (腕 B なら effB(k, x))] / MMPR + (腕 A なら a_offset(k)) [r_t]。
    band_detail: 交差の区間 k_hi = 10 (jump_col の列は k ≥ k₀ + 20,000 で 13)。calls に (腕, 出力, band) を足す。mutate(腕, 出力, d) で held・
    例外・帯の崩れを作る。"""
    def g(run, step, spec, ext, x_use, band):
        k = int(step)
        x = np.asarray(x_use, dtype=float)
        arm = "A" if band.get("fill_own") else "B"
        if calls is not None:
            calls.append((arm, k, band))
        xa = np.concatenate([np.asarray(X_OUT, dtype=float), x])
        lut = dict(zip(np.asarray(band["x"], dtype=float).tolist(), np.asarray(band["y_b"], dtype=float).tolist()))
        nat = tb_yb(k) + 0.0 * xa
        hit = np.array([v in lut for v in xa.tolist()])
        if not hit.all() and not band.get("fill_own"):
            raise ValueError("固定する帯の表に無い列 (偽物)")
        yb = np.array([lut[v] if h else nat[i] for i, (v, h) in enumerate(zip(xa.tolist(), hit))])
        e = 1e-3 + (tb_aux_mm(k, x) + (np.asarray(effB(k, x), dtype=float) if (effB and arm == "B") else 0.0)) / MMPR
        if arm == "A" and a_offset:
            e = e + a_offset(k)
        n, na = len(X_OUT), len(xa)
        khi = np.full(na, 10)
        if jump_col is not None and k >= TB_K0 + 20000:
            khi[n + jump_col] = 13
        det = {"y_star": [0.08] * na, "kappa": [0.0] * na, "kappa_med": [0.0] * na, "kappa_s": [0.0] * na, "y_ref": [0.05] * na,
               "y_lo": [0.079] * na, "y_hi": [0.081] * na, "k_lo": (khi - 1).tolist(), "k_hi": khi.tolist(), "n_cross": [1] * na,
               "k_second": [-1] * na, "n_fin": [40] * na, "fin_lo": [0] * na, "fin_hi": [39] * na, "j": khi.tolist(), "copy_consistent": [True] * na}
        d = {"delta_E": e.tolist(), "delta_target": e.tolist(), "held": [0] * len(x), "missing": [], "n_missing": 0, "n_held": 0, "lam": "1.0",
             "res_sha256": None, "seconds": 0.0, "y_b": yb[n:].tolist(), "y_star": (0.08 + 0.0 * x).tolist(), "y_b_natural": nat[n:].tolist(),
             "delta_r_raw": (0.9 * e).tolist(), "band_all": {"x": xa.tolist(), "y_b": yb.tolist(), "y_b_natural": nat.tolist()},
             "band_fill": {"n_table": int(hit.sum()), "n_own": int((~hit).sum()), "x_own_range": None, "edge_table_matches_fixed": True},
             "band_detail": det}
        if mutate:
            d = mutate(arm, k, d) or d
        return d
    return g


def time_judge(case: Path, extra=()) -> tuple:
    rc = HJ.main(["--diag-time-band", "--case-dir", str(case), "--design", str(ROOT / "design"), "--tools", str(TOOLS), *extra])
    return rc, HJ.jload(case / HJ.TIME_OUT)


def test_time_band():
    print("== N. 帯の時間方向の扱いの A/B (--diag-time-band; §5.1 #36p)", flush=True)
    import inspect
    gate_calls = []
    req = SR.required_window_steps(S_DIAG, 20000, 2500)
    # --- 分類 (単体) ---
    P = lambda Ds, err=None, bad=None: [{"step": k, "Dt_mm": (None if bad and i in bad else d), "pf_A": True, "pf_B": not (bad and i in bad),  # noqa: E731
                                         **({"error_B": "RuntimeError"} if err and i in err else {}), "x_at_Dt": 4.0}
                                        for i, (k, d) in enumerate(zip(req, Ds))]
    for name, per_, repro, head, want, code in (
            ("全 17 出力で D_time 0.003 mm → 結果 B (選び直しだけが 0.01 mm 超の時間変化を生む、を退ける)", P([0.003] * 17), True, True, "結果 B", "B"),
            ("ちょうど 0.01 mm → 結果 B (≤)", P([0.01] * 17), True, True, "結果 B", "B"),
            ("1 出力だけ 0.0101 mm → 結果 A (選び直しが予算を超える時間変化を生む)", P([0.0] + [0.001] * 15 + [0.0101]), True, True, "結果 A", "A_D"),
            ("腕 A の再現 (#36l の aux) の不成立 → 判定不能", P([0.001] * 17), False, True, "判定不能", "U_repro"),
            ("先頭の時刻 k₀ の再現の不成立 → 判定不能 (D が小さくても)", P([0.001] * 17), True, False, "判定不能", "U_head"),
            ("1 出力で held・非有限 (D_time なし) → 判定不能", P([0.001] * 17, bad={5}), True, True, "判定不能", "U_invalid"),
            ("抽出の例外・帯の固定の不成立 → 判定不能", P([0.001] * 17, err={8}), True, True, "判定不能", "U_error"),
            ("17 出力のうち 1 つが無い → 判定不能 (縮めない)", P([0.001] * 17)[:-1], True, True, "判定不能", "U_incomplete")):
        cl = HJ.classify_time(per_, req, repro, head)
        check(f"N 分類: {name}", cl["result"] == want and cl["code"] == code, f"{cl['result']} {cl['code']}: {cl['why'][:120]}")
    check("N 分類: 必要な出力の集合が無い → 判定不能", HJ.classify_time(P([0.001] * 17), [], True, True)["code"] == "U_required")
    # --- 梯子の交差の区間の写し (本物の deltastar.edge_band_positions と照合) ---
    y = 0.02 * (1.25 ** (1 / 8)) ** np.arange(120)

    def col(t):
        if t == 0:                                                       # 途中で交差 (梯子の 4 点目が非有限 = 元の添字の対応を試す)
            d = 1.0 - 0.5 * np.exp(-y / 0.03)
            d[3] = np.nan
            return d
        if t == 1:                                                       # 最初の点で満たす
            return np.ones_like(y)
        if t == 2:                                                       # 交差なし
            return y.copy()
        if t == 3:                                                       # 正・有限の点が 2 つだけ
            d = np.full_like(y, np.nan)
            d[:2] = 1.0
            return d
        return np.where(y < 0.05, 1.0, np.where(y < 0.1, 1 + (y - 0.05) / 0.05, 2.0))   # |R| < ε の区間が 2 つ
    kinds = [i % 5 for i in range(15)]
    yl, dl = [y.copy() for _ in kinds], [col(t) for t in kinds]
    pos = DS.edge_band_positions(np.arange(15.0), yl, dl, np.zeros(15), np.ones(15))
    lc = HJ.ladder_crossings(yl, dl)
    cc = HJ.ladder_copy_consistent({**lc, "y_star": pos["y_star"]})
    check(f"N 梯子の交差の区間の写し: 本物の edge_band_positions の y* が全列で区間に入る ({int(cc.sum())}/15)・途中の交差 (k_lo {lc['k_lo'][0]}, "
          f"k_hi {lc['k_hi'][0]}、4 点目の非有限を飛ばした元の添字)・最初の点 (k_lo −1, k_hi 0)・交差なし (n_cross 0, y* NaN)・点の不足 (n_cross −1)・"
          f"区間 2 つ (n_cross 2, k_second {lc['k_second'][4]})",
          cc.all() and lc["k_hi"][0] == lc["k_lo"][0] + 1 and lc["k_lo"][0] > 3 and lc["n_fin"][0] == 119
          and (lc["k_lo"][1], lc["k_hi"][1]) == (-1, 0) and lc["n_cross"][2] == 0 and not np.isfinite(pos["y_star"][2])
          and lc["n_cross"][3] == -1 and lc["n_cross"][4] == 2 and lc["k_second"][4] > lc["k_hi"][4])
    bad = dict(lc)
    bad["y_star"] = np.where(np.arange(15) == 0, pos["y_star"] * 1.1, pos["y_star"])
    check("N 梯子の写しの検査: y* が区間の外なら合わないと出る", not HJ.ladder_copy_consistent(bad)[0] and HJ.ladder_copy_consistent(bad)[1:].all())
    # --- 全経路 (偽物の抽出): #36l の偽物の記録 (aux の δ_E・帯が出力で動く) ---
    tbase = write_diag_case(TMP / "time_base")

    def mut36(arm, k, d):
        if arm == "aux":
            d["y_b"] = [tb_yb(k)] * len(d["y_b"])
        return d
    with Patch((HJ, "HOOKS", {"handoff_gate": gate_boom(gate_calls), "diag_snapshot": fake_diag_snap(
            lambda arm, k, x: tb_aux_mm(k, x) if arm == "aux" else 0.0 * x, mutate=mut36)}), (DL, "handoff_gate", gate_boom(gate_calls))):
        rc36, R36 = diag_judge(tbase)
    check(f"N 前提: #36l の偽物の記録 (aux の δ_E が窓で 0.004 mm 動き、帯も動く; rc {rc36})", rc36 == 0 and (tbase / HJ.DIAG_NPZ).is_file())

    def run_time(name, snap, prep=None):
        case = TMP / f"time_{name}"
        shutil.copytree(tbase, case, symlinks=True)
        if prep:
            prep(case)
        with Patch((HJ, "HOOKS", {"handoff_gate": gate_boom(gate_calls), "diag_snapshot_time": snap}), (DL, "handoff_gate", gate_boom(gate_calls))):
            return time_judge(case)
    calls = []
    rc, R = run_time("B", fake_time_snap(calls=calls, jump_col=5))
    arms = [a_ for a_, _, _ in calls]
    bB = [b_ for a_, _, b_ in calls if a_ == "B"]
    xa_all = np.concatenate([np.asarray(X_OUT), np.load(tbase / HJ.DIAG_NPZ)["x_use"]])
    okB = all(b_["fill_own"] is False and np.array_equal(np.asarray(b_["x"], dtype=float), xa_all)
              and np.array_equal(np.asarray(b_["y_b"], dtype=float), np.full(len(xa_all), tb_yb(TB_K0))) for b_ in bB)
    tw = R.get("two_windows", {})
    check(f"N 全経路 (偽物の抽出): 両腕の時間変化が同じ (腕 A は帯が動き、腕 B は k₀ の帯) → D_time {R.get('D_time_mm', {}).get('max')} mm → 結果 B "
          f"(rc {rc})・34 抽出 (腕 A 17 → 腕 B 17)・腕 A の再現 (ビット同一)・先頭の時刻の再現・update_conditions_met 偽・関門を呼ばない",
          rc == 0 and R["result"] == "結果 B" and R["code"] == "B" and arms == ["A"] * 17 + ["B"] * 17 and R["D_time_mm"]["max"] == 0.0
          and R["repro_A"]["ok"] is True and R["repro_A"]["bitwise_all"] is True and R["head_repro"]["ok"] is True and R["k0"] == TB_K0
          and R["update_conditions_met"] is False and not gate_calls and R["extractions"]["n"] == 34, str(R.get("why"))[:300])
    check("N 帯の固定が全時刻に同じ値: 腕 B の 17 回の抽出すべてに k₀ の腕 A の帯 (上流の列を含む全列、fill_own なし) が渡り、使った帯も全出力で同じ "
          f"(band_B_fixed_ok {sum(1 for r_ in R['per_output'] if r_.get('band_B_fixed_ok'))}/17)・腕 A には固定の表を渡さない",
          okB and len(bB) == 17 and all(r_.get("band_B_fixed_ok") is True for r_ in R["per_output"]) and R["band_B"]["n_cols_x_lt_0"] == len(X_OUT)
          and all(len(b_["x"]) == 0 and b_["fill_own"] is True for a_, _, b_ in calls if a_ == "A"))
    npz = np.load(TMP / "time_B" / HJ.TIME_NPZ)
    check(f"N 併記: 2 窓 (腕 A の δ_E の窓の最大 − 最小 {tw.get('A', {}).get('delta_E', {}).get('W1', {}).get('range_max_mm')} mm・平均の差 "
          f"{tw.get('A', {}).get('delta_E', {}).get('mean_diff', {}).get('max_mm')} mm = 窓 20,000 で 0.002 mm)・代表の列の check_quasisteady (両腕)・"
          f"梯子の変化 (飛び {R.get('ladder_A', {}).get('flags_total', {}).get('jump_gt1')} 列)・配列",
          abs(tw["A"]["delta_E"]["W1"]["range_max_mm"] - 0.002) < 1e-9 and abs(tw["A"]["delta_E"]["mean_diff"]["max_mm"] - 0.002) < 1e-9
          and abs(tw["B"]["delta_E"]["W2"]["range_max_mm"] - 0.002) < 1e-9 and tw["A"]["delta_E"]["W1"]["n_outputs"] == 9
          and all(set(R["quasisteady"][a_]["verdicts"]) == set(R["quasisteady"]["cols"]) for a_ in ("A", "B"))
          and {"dE_40", "dE_70", "dE_94", "dE_x3to5_Dtmax"} == set(R["quasisteady"]["cols"])
          and R["ladder_A"]["flags_total"]["jump_gt1"] == 1 and R["ladder_A"]["flags_total"]["interval_changed"] == 1
          and R["ladder_B_equals_A"]["k_hi"] is True and {"A_lad_k_hi", "B_band_kappa_s", "D_time_kx_mm", "A_delta_r_raw", "B_y_b_natural",
                                                          "A_delta_E_range_W1_x_mm"} <= set(npz.files)
          and np.array_equal(npz["A_delta_E"], np.load(tbase / HJ.DIAG_NPZ)["aux_delta_E"]), str(R.get("quasisteady"))[:300])
    effA = lambda k, x: np.where((x >= 3.0) & (x < 5.0), 0.02 * (k - TB_K0) / 40000.0, 0.0)  # noqa: E731
    rc, R = run_time("A", fake_time_snap(effB=effA))
    sec = {tuple(s_["x_range"][:1]): s_ for s_ in R.get("sections", [])}
    check(f"N 全経路: 腕 B だけ x 3〜5 で窓の終わりに 0.02 mm 動く → D_time {R.get('D_time_mm', {}).get('max')} mm → 結果 A ({R.get('code')})・"
          f"区間 3〜5 が最大・代表の列 dE_x3to5_Dtmax が x 3〜5",
          rc == 0 and R["result"] == "結果 A" and R["code"] == "A_D" and abs(R["D_time_mm"]["max"] - 0.02) < 1e-12 and 3.0 <= R["D_max"]["x"] < 5.0
          and abs(sec[(3.0,)]["D_time"]["max_mm"] - 0.02) < 1e-12 and sec[(1.0,)]["D_time"]["max_mm"] == 0.0
          and 3.0 <= R["quasisteady"]["cols"]["dE_x3to5_Dtmax"]["x"] < 5.0, str(R.get("why"))[:300])

    def held_b(arm, k, d):
        if arm == "B" and k == req[4]:
            d["held"] = list(d["held"])
            d["held"][len(d["held"]) // 2] = 1
        return d
    rc, R = run_time("held", fake_time_snap(mutate=held_b))
    check(f"N 全経路: 腕 B の 1 出力の比べる列の中に held → 判定不能 ({R['code']}、rc {rc})", rc == 2 and R["result"] == "判定不能" and R["code"] == "U_invalid")

    def nan_a(arm, k, d):
        if arm == "A" and k == req[10]:
            d["delta_E"] = list(d["delta_E"])
            d["delta_E"][7] = float("nan")
        return d
    rc, R = run_time("nan", fake_time_snap(mutate=nan_a))
    check(f"N 全経路: 腕 A の 1 列が非有限 (#36l は有限) → 腕 A の再現の不成立で判定不能 ({R['code']})", rc == 2 and R["code"] == "U_repro")
    calls = []
    rc, R = run_time("repro", fake_time_snap(calls=calls, a_offset=lambda k: 1e-9 if k == req[3] else 0.0))
    check(f"N 腕 A の再現の不成立 (1 出力で δ_E が #36l の aux と 1e−9 r_t 違う) → 判定不能 ({R['code']}、rc {rc})・腕 B は回さない "
          f"({sum(1 for a_, _, _ in calls if a_ == 'B')} 回)",
          rc == 2 and R["code"] == "U_repro" and abs(R["repro_A"]["max_abs_delta_E_rt"] - 1e-9) < 1e-15 and not any(a_ == "B" for a_, _, _ in calls))
    rc, R = run_time("head", fake_time_snap(effB=lambda k, x: 0.001 + 0.0 * x))
    check(f"N 先頭の時刻の再現の不成立 (腕 B が全出力で 0.001 mm ずれる = k₀ でも違う; D_time は 0) → 判定不能 ({R['code']}、rc {rc})",
          rc == 2 and R["code"] == "U_head" and R["head_repro"]["ok"] is False)

    def unfix(arm, k, d):
        if arm == "B" and k == req[6]:
            d["band_all"] = {**d["band_all"], "y_b": [v + 1e-3 for v in d["band_all"]["y_b"]]}
        return d
    rc, R = run_time("unfix", fake_time_snap(mutate=unfix))
    check(f"N 帯の固定の不成立 (腕 B の 1 出力の使った帯が k₀ の帯と違う) → 判定不能 ({R['code']}、rc {rc})",
          rc == 2 and R["code"] == "U_error" and any("error_band_B" in r_ for r_ in R["per_output"]))

    def exc_b(arm, k, d):
        if arm == "B" and k == req[9]:
            raise ValueError("固定する帯の表に無い列 (試験)")
        return d
    rc, R = run_time("exc", fake_time_snap(mutate=exc_b))
    check(f"N 腕 B の 1 出力の抽出の不成立 → 判定不能 ({R['code']}、rc {rc})", rc == 2 and R["code"] == "U_error")
    for name, prep, want in (
            ("#36l の記録が無い", lambda cs: (cs / HJ.DIAG_OUT).unlink(), "#36l の記録"),
            ("参照 aux (run_0950) の sha が固定と違う", lambda cs: (cs / HJ.EU_RUN / "res_60000.h5").write_bytes(b"changed"), "参照 aux"),
            ("#36l と同じ照合: 窓の 1 出力の res が無い", lambda cs: (cs / HJ.NS_RUN / f"res_{req[8]}.h5").unlink(), "が無い")):
        calls = []
        n_u = len(list(TMP.glob("time_u*")))
        rc, R = run_time(f"u{n_u}", fake_time_snap(calls=calls), prep=prep)
        check(f"N 入力の照合: {name} → 判定不能 (rc {rc}、{R.get('code')})・抽出しない ({len(calls)} 回)", rc == 2 and R["code"] == "U_input"
              and not calls and any(want in u for u in R["undecidable"]), str(R.get("undecidable"))[:300])
    srcs = {f.__name__: inspect.getsource(f) for f in (HJ.diag_time_band, HJ.classify_time, HJ.ladder_crossings, HJ.ladder_copy_consistent,
                                                        HJ.window_stats, HJ.ladder_change_summary, HJ._cmp_arr)}
    bad_ = [n for n, s_ in srcs.items() if any(w_ in s_ for w_ in ("unlink(", "rmtree(", "os.remove(", "delete_after(", "handoff_gate", "update_wall"))]
    check(f"N 帯の時間方向の A/B の関数に消す処理・関門・壁の更新の呼び出しが無い ({bad_ or 'なし'})", not bad_ and not gate_calls)
    # --- 本物の抽出 (diag_snapshot_fixed の record_ladder、プロセス並列) ---
    rb = TMP / "diag_real_B"                                    # L の本物の #36l の記録があれば使う
    if not (rb / HJ.DIAG_OUT).is_file():
        rb = write_diag_case(TMP / "time_real_src", shape_b=0.01)
        with Patch((HJ, "HOOKS", {"handoff_gate": gate_boom(gate_calls)}), (DL, "handoff_gate", gate_boom(gate_calls))):
            diag_judge(rb)
    ns = rb / HJ.NS_RUN
    rec = json.loads((ns / SR.RECORD).read_text())
    ext = {k: v for k, v in rec["delta_series"]["extract"].items() if k != "source"}
    xu = np.asarray(rec["delta_series"]["x"], dtype=float)
    spec = {"run": str(rb / HJ.EU_RUN), "res": "res_60000.h5", "sha256": HJ.sha256_file(rb / HJ.EU_RUN / "res_60000.h5")}
    f0 = HJ.diag_snapshot_fixed(ns, S_DIAG, spec, ext, xu, {"x": [], "y_b": [], "fill_own": True})
    f1 = HJ.diag_snapshot_fixed(ns, S_DIAG, spec, ext, xu, {"x": [], "y_b": [], "fill_own": True, "record_ladder": True})
    same = all(np.array_equal(np.asarray(f0[q], dtype=float), np.asarray(f1[q], dtype=float), equal_nan=True)
               for q in ("delta_E", "delta_target", "held", "y_b", "y_star", "delta_r_raw"))
    bd = f1.get("band_detail") or {}
    check(f"N 本物の抽出: record_ladder は値を変えない (δ_E・δ_target・held・帯・y*・生の δ_r がビット同一)・band_detail は帯を決めた全列 "
          f"({len(bd.get('k_hi', []))} 列) で、交差の区間の写しが抽出器の y* と全列で合う ({int(np.sum(bd.get('copy_consistent', [])))})・"
          "record_ladder なしでは band_detail を付けない",
          same and "band_detail" not in f0 and len(bd.get("k_hi", [])) == len(f1["band_all"]["x"]) and all(bd["copy_consistent"])
          and np.array_equal(np.asarray(bd["y_star"], dtype=float)[np.isin(f1["band_all"]["x"], xu)], np.asarray(f1["y_star"], dtype=float), equal_nan=True))
    rcase = TMP / "time_real"
    shutil.copytree(rb, rcase, symlinks=True)
    with Patch((HJ, "HOOKS", {"handoff_gate": gate_boom(gate_calls)}), (DL, "handoff_gate", gate_boom(gate_calls))):
        rc, R = time_judge(rcase)
    px = (R.get("extractions") or {}).get("B") or {}
    dt = (R.get("D_time_mm") or {}).get("max")
    check(f"N 本物の全経路 ({px.get('mode')} {px.get('jobs')}): 腕 A が #36l の aux をビット同一で再現 ({R.get('repro_A', {}).get('bitwise_all')})・"
          f"先頭の時刻の再現 ({R.get('head_repro', {}).get('ok')})・34 抽出・腕 B の帯が全出力で k₀ の帯・梯子の写しが全出力で合う "
          f"(合わない列 {R.get('ladder_A', {}).get('flags_total', {}).get('copy_inconsistent')})・腕 B の梯子 = 腕 A・D_time {dt} mm で分類 {R.get('result')}",
          rc == 0 and R["repro_A"]["ok"] is True and R["repro_A"]["bitwise_all"] is True and R["head_repro"]["ok"] is True
          and R["head_repro"]["delta_E"]["bitwise"] is True and R["extractions"]["n"] == 34 and px.get("mode", "").startswith("process")
          and all(r_.get("band_B_fixed_ok") is True for r_ in R["per_output"]) and R["ladder_A"]["flags_total"]["copy_inconsistent"] == 0
          and all(v for q, v in R["ladder_B_equals_A"].items() if q != "note")
          and R["result"] == ("結果 A" if dt > HJ.TOL_MM else "結果 B") and not gate_calls and R["band_B"]["n_cols_x_lt_0"] > 0,
          str(R.get("why"))[:300])


# ====================================================================================== O. 受け渡しの関門の局所 y₁⁺ (§5.1 #36q)
RUN0951_OVER_AREA_M2, RUN0951_AREA_M2 = 0.850252596933461, 28.996474362873958   # h370_wallres_split.json (§5.1 #36o) の wall の値


def run0951_like(d: Path, S: int = 152500, old: int = 97500) -> Path:
    """run_0951 の停止の出力の形の合成: 壁 1 枚 (wall・physID 3・wall_isothermal) 4,719 点のうち先頭 205 点が y₁⁺ 1.5 (残り 0.5)、
    超過の面積の割合が記録 (h370_wallres_split.json の 0.850252596933461 / 28.996474362873958 m² = 2.932 %) と同じになるよう超過の区間と
    残りの間隔を変える (半径一定の直線の壁: 超過の面積 = L1 + h2/2、全長 1)。出口の出力あり・古い step (97,500) の出力もある。
    実際の壁の形・y₁⁺ の分布ではない (関門の経路と重みの確認用)。"""
    import test_gate_bad_input as TGB
    n, n_over = 4719, 205
    m = n - n_over
    f = RUN0951_OVER_AREA_M2 / RUN0951_AREA_M2
    L1 = (f - 0.5 / m) / (1.0 - 0.5 / m)
    xs = np.concatenate([np.linspace(0.0, L1, n_over), L1 + (1.0 - L1) * np.arange(1, m + 1) / m])
    y1p = np.where(np.arange(n) < n_over, 1.5, 0.5)
    return Path(TGB.make_wall_run(str(d), {3: ("wall", "wall_isothermal", 0.2486, xs, y1p)}, {"wall": [old, S]}))


def test_wall_resolution_gate():
    print("== O. 受け渡しの関門の局所 y₁⁺ (§5.1 #36q)", flush=True)
    import h5py
    import test_gate_bad_input as TGB
    from forge_design.evaluate import runner as RN
    C = TGB.wall_cases()
    S, OLD = 2000, 1000
    both = {"wall_a": [OLD, S], "wall_b": [OLD, S]}
    root = TMP / "wallres_gate"

    def mk(name, walls, steps, **kw) -> Path:
        return Path(TGB.make_wall_run(str(root / name), walls, steps, **kw))

    def arg(g, key):
        c = g.get("cmd") or []
        return c[c.index(key) + 1] if key in c else None

    # --- 事前登録の受入 ---
    d = mk("control", C["control"], both)
    g = DL.wall_resolution_verdict(d, S)
    check("O 対照 (2 枚の no-slip 壁とも停止の step に有効な出力・出口の出力あり) → PASS; --step s・--groups は bcond の wall* だけ "
          "(wall_a,wall_b)・--weight area・ツールの行が両壁とも step s",
          g["kind"] == "PASS" and g["walls"] == {"wall_a": 3, "wall_b": 5} and arg(g, "--groups") == "wall_a,wall_b"
          and arg(g, "--weight") == "area" and arg(g, "--step") == str(S) and g["rows"] == {"wall_a": S, "wall_b": S}
          and all(e["ok"] for e in g["dumps"].values()), f"{g['kind']} {g['verdict']} rows={g.get('rows')}")
    rc_all, v_all, _ = TGB.run_cwr(str(d), "--step", S, "--weight", "area")
    check("O 出口の出力があっても壁だけで判定 (同じ入力をツールに群なしで渡すと出口で判定不能 = 2026-10-14 までの関門の呼び方)",
          g["kind"] == "PASS" and "outlet" not in (arg(g, "--groups") or "") and rc_all == 2 and "outlet" in v_all, f"群なし rc {rc_all} {v_all}")
    d = mk("b_old_only", C["control"], {"wall_a": [OLD, S], "wall_b": [OLD]})
    g = DL.wall_resolution_verdict(d, S)
    check("O 片方の壁 (wall_b) の停止の step の出力だけを欠く (古い step はある) → UNDETERMINED、ツールは回さない (古い step に戻らない)",
          g["kind"] == "UNDETERMINED" and g["cmd"] is None and g["dumps"]["wall_a"]["ok"] and not g["dumps"]["wall_b"]["ok"]
          and g["dumps"]["wall_b"]["steps_present"] == [OLD] and "wall_b" in g["why"], g.get("why"))
    d = mk("old_only", C["control"], {"wall_a": [OLD], "wall_b": [OLD]})
    g = DL.wall_resolution_verdict(d, S)
    check("O 停止の step が無い (両壁とも古い step だけ) → UNDETERMINED", g["kind"] == "UNDETERMINED" and g["cmd"] is None
          and not any(e["ok"] for e in g["dumps"].values()), g.get("why"))
    for name, want, want_pt in (("inv_point_fail", "PASS", 1), ("inv_area_fail", "FAIL", 0)):
        d = mk(name, C[name], both)
        g = DL.wall_resolution_verdict(d, S)
        rp, vp, _ = TGB.run_cwr(str(d), "--step", S, "--groups", "wall_a,wall_b", "--weight", "point")
        check(f"O 面積の割合と点数の割合で合否が逆になる入力 ({name}) → 面積で判定 {want} (超過の面積 {g.get('over_area_pct')} %; "
              f"同じ入力を点数で数えると rc {rp})",
              g["kind"] == want and rp == want_pt and arg(g, "--weight") == "area" and g.get("over_area_allow_pct") == 2.0, f"{g['verdict']} / point {vp}")
    # --- ほかの判定不能 ---
    d = mk("b_absent", C["control"], {"wall_a": [OLD, S]})
    g = DL.wall_resolution_verdict(d, S)
    check("O 壁 (wall_b) の出力が 1 つも無い → UNDETERMINED", g["kind"] == "UNDETERMINED" and g["dumps"]["wall_b"]["n_steps_present"] == 0, g.get("why"))
    d = mk("no_twall", C["control"], both)
    with h5py.File(d / f"res_wall_b_5_{S}.h5", "a") as f:
        del f["VALUE/twall_x"]
    g = DL.wall_resolution_verdict(d, S)
    check("O 停止の step の壁の出力に twall_x が無い → UNDETERMINED", g["kind"] == "UNDETERMINED" and "twall_x" in g["dumps"]["wall_b"]["why"], g.get("why"))
    d = mk("attr_step", C["control"], both)
    with h5py.File(d / f"res_wall_a_3_{S}.h5", "a") as f:
        f.attrs["step"] = np.int64(OLD)
    g = DL.wall_resolution_verdict(d, S)
    with h5py.File(d / f"res_wall_a_3_{S}.h5", "a") as f:
        f.attrs["step"] = np.int64(S)
    g2 = DL.wall_resolution_verdict(d, S)
    check("O 壁の出力の属性 step がファイル名の step と違う (時刻違い) → UNDETERMINED、同じなら PASS",
          g["kind"] == "UNDETERMINED" and "時刻違い" in g["dumps"]["wall_a"]["why"] and g2["kind"] == "PASS", g.get("why"))
    d = mk("control2", C["control"], both)
    g = DL.wall_resolution_verdict(d, None)
    check("O 停止の step が無い (step=None) → UNDETERMINED (最新の出力で判定しない)", g["kind"] == "UNDETERMINED" and g["cmd"] is None, g.get("why"))
    (d / "bcondConfig.yaml").write_text((d / "bcondConfig.yaml").read_text().replace("kind: wall_isothermal", "kind: slip").replace("kind: wall,", "kind: slip,"))
    g = DL.wall_resolution_verdict(d, S)
    check("O bcondConfig に no-slip 壁 (kind が wall*) が無い → UNDETERMINED", g["kind"] == "UNDETERMINED" and g["walls"] == {}, g.get("why"))
    (d / "bcondConfig.yaml").unlink()
    g = DL.wall_resolution_verdict(d, S)
    check("O bcondConfig.yaml が無い → UNDETERMINED", g["kind"] == "UNDETERMINED" and g["walls"] is None, g.get("why"))
    # --- ツールの壁ごとの行の照合 (偽のツール: 終了コード 0・VERDICT PASS) ---
    fake = root / "fake_tools"
    fake.mkdir(parents=True, exist_ok=True)
    d = root / "control"
    for label, rows in (("wall_b の行が無い", [("wall_a", S)]), ("wall_b の行が古い step", [("wall_a", S), ("wall_b", OLD)])):
        body = "".join(f'print("  {n:<12s} step {st:6d}  y1 = 1.000e-05 m  y1+ 平均   0.500 / p99   0.500 / 最大   0.500")\n' for n, st in rows)
        (fake / "check_wall_resolution.py").write_text(body + 'print("VERDICT: PASS (偽物)")\n')
        with Patch((RN, "FORGE_TOOLS", fake)):
            g = DL.wall_resolution_verdict(d, S)
        check(f"O ツールの壁ごとの行が揃わない ({label}) → 終了コード 0 でも UNDETERMINED", g["kind"] == "UNDETERMINED" and g["rc"] == 0
              and "揃っていない" in g["why"], g.get("why"))
    # --- 受け渡しの関門 (他の条件は偽物で合格) ---
    def gate(run):
        with Patch((DL, "stop_state", lambda r: {"status": "REACHED", "stop_step": S, "source": "試験の偽物"}),
                   (DL, "convergence_verdict", lambda r: {"cmd": ["試験の偽物"], "rc": 0, "kind": "PASS", "verdict": [], "output": "偽物\n"}),
                   (DL, "conservation_check", lambda rec, s: {"pass": True}), (DL, "two_window_check", lambda rec, s: {"pass": True})):
            return DL.handoff_gate(run, euler_ref={"run": "試験", "res": "res_1.h5", "sha256": "0"}, allow_not_reached="試験",
                                   allow_unconverged="試験", allow_unsettled="試験", extract=False)
    g_ok = gate(root / "control")
    g_ng = gate(root / "b_old_only")
    txt = (root / "b_old_only" / DL.HANDOFF_WALL_RESOLUTION).read_text()
    check("O 受け渡しの関門: 停止の step の壁の出力の欠け (UNDETERMINED) は 3 つの上書きでも止まる・対照は局所 y₁⁺ で止めない・"
          "HANDOFF_WALL_RESOLUTION.txt に理由",
          any("局所 y₁⁺" in b for b in g_ng["blocked_by"]) and not g_ng["proceed"] and not any("局所 y₁⁺" in b for b in g_ok["blocked_by"])
          and g_ok["checks"]["wall_resolution"]["pass"] and "UNDETERMINED" in txt and "wall_b" in txt,
          str([b for b in g_ng["blocked_by"] if "局所" in b])[:300])
    # --- run_0951 の停止の出力の形の合成 ---
    d = run0951_like(root / "run0951_like")
    g = DL.wall_resolution_verdict(d, 152500)
    prof = root / "run0951_like_prof.csv"
    TGB.run_cwr(str(d), "--step", 152500, "--groups", "wall", "--weight", "area", "--profile-csv", prof)
    a = np.loadtxt(root / "run0951_like_prof_wall.csv", delimiter=",", skiprows=1)
    over_area = 100.0 * float(a[a[:, 4] > 1.0, 9].sum() / a[:, 9].sum())
    over_pt = 100.0 * float(np.mean(a[:, 4] > 1.0))
    rc_old, v_old, _ = TGB.run_cwr(str(d), "--step", 152500)
    check(f"O run_0951 の停止の出力の形の合成 (壁 4,719 点・超過 205 点・出口・古い step) → 新しい関門は FAIL (面積の割合 {over_area:.3f} %、"
          f"ツールの表示 {g.get('over_area_pct')} %; 点数なら {over_pt:.3f} %)、2026-10-14 までの呼び方 (群・重みなし) は rc {rc_old}",
          g["kind"] == "FAIL" and abs(over_area - 2.932) < 5e-4 and abs(over_pt - 4.344) < 5e-4 and g.get("over_area_pct") == 2.9
          and rc_old == 2 and arg(g, "--groups") == "wall", f"{g['verdict']} / 旧 {v_old}")
    real = os.environ.get("H370_RUN0951_DIR")
    if real and Path(real, "res_wall_3_152500.h5").is_file():
        g = DL.wall_resolution_verdict(Path(real), 152500)
        check(f"O run_0951 の停止の出力の写し ({real}) → 新しい関門は FAIL (超過の面積 {g.get('over_area_pct')} %)",
              g["kind"] == "FAIL" and g.get("over_area_pct") == 2.9 and arg(g, "--groups") == "wall", str(g["verdict"]))
    else:
        print("     (H370_RUN0951_DIR が無いので run_0951 の実データは見ない; 合成の形で代えた)", flush=True)


def main() -> int:
    for t in (test_constants, test_wall_output_and_configs, test_euler_stage_configs, test_euler_qualify, test_classify_ab, test_ns_judge,
              test_prelaunch_counterexamples, test_parallel_extraction, test_run_orchestration, test_ns_phase_and_parts, test_euler_wall_restore,
              test_diag_path, test_fixed_band, test_time_band, test_wall_resolution_gate):
        try:
            t()
        except Exception:  # noqa: BLE001
            import traceback
            traceback.print_exc()
            check(f"{t.__name__} が例外", False)
    print(f"\n{'ALL PASS' if FAIL == 0 else f'FAIL {FAIL}'}  (一時ファイル {TMP})")
    if FAIL == 0:
        shutil.rmtree(TMP, ignore_errors=True)
    return 0 if FAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
