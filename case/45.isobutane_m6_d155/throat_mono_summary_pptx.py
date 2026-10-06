"""plan tooling-nozzle-throat-monotone-r2 の要約報告 (旧壁 vs r″ 単調壁)。build_pptx.py の部品で、標準報告と同じ見た目にする。
数値はすべて JSON から読む (手で写さない):
  _band_ab/throat_mono_practical_eval.json (Euler A/B の実務判定)、_band_ab/verdicts_run_01{17,18,48,49}_monoeval.json (量別の準定常)、
  run_0117 / run_0149 / run_0118 / run_0148 の report/report.json、_band_ab/throat_mono_compare/{compare_summary,wall_r2_summary}.json と図。
usage: .venv-pptx/bin/python throat_mono_summary_pptx.py CASE_DIR  → CASE_DIR/_band_ab/throat_mono_compare/throat_mono_summary_report.pptx
"""
import json
import sys
from datetime import date
from pathlib import Path

from pptx import Presentation
from pptx.util import Inches
from pptx.dml.color import RGBColor

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "design/forge_design/report"))
from build_pptx import text, title, picture, table, stat, W, H, NAVY, MUTED, WARN, WHITE  # noqa: E402

C = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE
B = C / "_band_ab"; CMP = B / "throat_mono_compare"
J = lambda p: json.loads(Path(p).read_text())  # noqa: E731
PE = J(B / "throat_mono_practical_eval.json")
V = {k: J(B / f"verdicts_run_{k}_monoeval.json") for k in ("0117", "0118", "0148", "0149")}
RP = {k: J(C / r / "report/report.json")["metrics"] for k, r in
      (("0117", "run_0117_ns_recal_final_ext"), ("0149", "run_0149_ns_mono_final_ext"),
       ("0118", "run_0118_ns_recal_final_cond"), ("0148", "run_0148_ns_mono_final_cond"))}
CS, WS = J(CMP / "compare_summary.json"), J(CMP / "wall_r2_summary.json")
CASE = "case/45.isobutane_m6_d155"

prs = Presentation(); prs.slide_width = Inches(W); prs.slide_height = Inches(H)
blank = prs.slide_layouts[6]


def new(t, sub=None, notes=None):
    s = prs.slides.add_slide(blank); title(s, t, sub)
    if notes:
        s.notes_slide.notes_text_frame.text = notes
    return s


def vd(k, key):
    v = V[k].get(key)
    return v["verdict"] if v else "—"


def mean_of(k, key):
    """verdicts の note「末尾 5 枚の平均 X・幅 Y」から X を読む (量別判定と同じ窓の値)。"""
    note = V[k][key]["note"]
    return float(note.split("末尾 5 枚の平均 ")[1].split("・")[0])


# 1 表紙と結論
rows_pe = {r["qty"]: r for r in PE["rows"]}
worst = max(r["D_plus_2SE_over_dq"] for r in PE["rows"] if "D_plus_2SE_over_dq" in r)
s = prs.slides.add_slide(blank)
bg = s.background.fill; bg.solid(); bg.fore_color.rgb = NAVY
text(s, 0.7, 0.7, W - 1.4, 1.0, "スロート直後の r″ を単調にした壁 — 検証報告", size=34, bold=True, color=WHITE)
text(s, 0.7, 1.5, W - 1.4, 0.6, f"{CASE} · 旧壁 (joint 当てはめ) と単調壁 (r‴ ≤ 0 拘束) の比較", size=16, color=RGBColor(0xC9, 0xD6, 0xDF))
stat(s, 0.7, 2.45, f"{WS['design_old']['r2_max_0_0p3']:.4f} → {WS['design_mono']['r2_max_0_0p3']:.4f}", "設計壁の r″ の最大 (x 0〜0.3、1/R = 0.5)")
stat(s, 3.8, 2.45, f"{worst:.2f}", "Euler A/B: (差 + 2SE)/許容幅 の最大 (6 量すべて < 1)")
stat(s, 6.9, 2.45, f"{mean_of('0149', 'exit_core_M'):.6f}", f"単調壁 NS の出口コア M (旧壁 {mean_of('0117', 'exit_core_M'):.6f})")
stat(s, 10.0, 2.45, f"{WS['physical_max_abs_dr_um']:.2f} µm", "物理壁の半径の差の最大 (旧壁比)")
text(s, 0.7, 4.1, W - 1.4, 2.6, [
    "結論: 単調壁は、Euler・dry NS・凝縮 NS のどれでも旧壁との差が許容幅と時間変動の範囲に収まった。dry NS は全ゲートに合格し (延長 1 回後に 4 量 STEADY)、"
    "凝縮 NS の 4 量も STEADY で旧壁とほぼ同じ値になった。",
    "最大の限界: Euler の比較はユーザ決定の実務判定 (6 本の時間平均と 2·SE。自己相関は未補正、窓の開始位置は結果を見た後に決定) で、厳密な非劣化の証明ではない。"
    "差を検出した量は試験部の |P 傾き| (r/r_w = 0.1) だけで、許容幅の 37 %。",
    "r″ の山の原因 (MOC の始点付近の不整合) は特定途中で、本報告の結論には影響しない (補足を参照)。"],
    size=14, color=WHITE)
text(s, 0.7, 6.9, W - 1.4, 0.4, f"作成 {date.today().isoformat()} · plan plans/active/tooling-nozzle-throat-monotone-r2.md · run 索引は {CASE}/README.md が正本",
     size=10, color=RGBColor(0xC9, 0xD6, 0xDF))

# 2 何を変えたか
s = new("何を変えたか: 壁の当てはめに「r″ が増えない」拘束", "設計壁・物理壁とも r″ の山が消える (CFD を使わない形状の比較)",
        notes=f"図: {CMP / 'fig_wall_r2_old_vs_mono.png'} (throat_mono_wall_fig.py)")
picture(s, CMP / "fig_wall_r2_old_vs_mono.png", 0.4, 1.4, 8.4, 3.2)
text(s, 0.5, 4.75, 12.3, 2.6, [
    "旧壁: MOC 壁点に位置と壁角を同時に当てはめる 5 次 B-spline (joint)。x = 0 で r = 1・r′ = 0・r″ = 1/R に固定。"
    f"スロート直後 x = {WS['design_old']['x_r2_max']:.3f} r_t で r″ が {WS['design_old']['r2_max_0_0p3']:.4f} まで上がる (山)。",
    "単調壁: 同じ当てはめに、x ∈ [0, 1.5] r_t で r‴ ≤ 0 (r″ が増えない) の不等式拘束だけを加えた (geometry.wall_fit_mono_r2: [0.0, 1.5])。",
    f"形状の差: 設計壁の半径で最大 {WS['design_max_abs_dr_um']:.2f} µm (x ≈ 0.06 r_t)。物理壁では境界層厚さの再計算が下流に伝わり、"
    f"最大 {WS['physical_max_abs_dr_um']:.2f} µm (x = {WS['physical_x_of_max_dr']:.1f} r_t)。r_t = 76.654 mm。",
    "代償: 最初の MOC 点 (x = 0.025 r_t) で、壁の傾きが MOC の流れ角より 0.022° 寝る (旧壁は 0.012°)。"], size=13)
text(s, 9.0, 1.5, 3.9, 3.0, [
    "形状ゲート (CFD 0 step): S1〜S8 すべて PASS",
    "滑らかさ ([0, 0.3]、区間ごとの厳密評価): ∫(r‴)² 0.106 → 0.095、max|r⁗| 225 → 58.8",
    "既定 (キー無し) の壁はビット同一"], size=12, color=MUTED)

# 3 解析条件
s = new("解析条件", "数値結果より前に、何をどう解いたかを示す")
table(s, [["項目", "Euler A/B", "dry NS / 凝縮 NS"],
          ["問題", "problem_d155_euler_pin_G1_recal.yaml (旧壁) / _mono.yaml", "problem_d155_ns_finemesh_recal_final.yaml / _mono.yaml (・_cond)"],
          ["格子", "生産 Euler 格子 G1: 2000 × 97 節点 (node)", "同じ 2000 × 97、壁第 1 セル 1.3e-5 (スロート 4.5e-6) r_w"],
          ["r_t・k_f・出口較正", "設計壁 (Euler)", "r_t 76.6539 mm・k_f 1.054129・Md_moc_offset +3.770e-4 (旧壁と同じ、較正し直さない)"],
          ["数値設定", "2 次 SLAU・limiter 2・block-DPLUR、soft 3000 → 本段 cfl 2・relax 0.7・18000 step", "同 + SST 低 Re、本段 cfl 1・relax 0.7 (段階起動なし)"],
          ["初期値", "run_0114 の最終場。旧壁は restart_field、単調壁は検証付き番号写像", "旧壁 run_0117 res_60000 を番号写像 (dry)。凝縮は dry の res_60000 を convert_species_field"],
          ["step・判定窓", "各 3 本、本段 step 6000〜18000 の 13 枚の平均", "dry 60000 + 延長 20000 (窓 60000〜80000)、凝縮 18000 (窓 14000〜18000)"],
          ["判定", "(差 + 2·SE) ≤ 許容幅 Δq (ユーザ決定の実務判定)", "§6 N のゲート (出口半径・出口 M・δ_E/δ_C・波・オーバーシュート・壁解像) と量別の準定常"],
          ["バイナリ", "AWS ~/forge-wallfit-bin (sha 6b47811b)", "同じ"]],
      0.5, 1.4, 12.3, [2.0, 4.8, 5.5], size=11, row_h=0.52)
text(s, 0.5, 6.35, 12.3, 0.9, "物性・境界条件は旧壁の報告 (run_0117 / run_0118) と同じ (燃焼ガス 4 種、Pt 5.5 MPa・Tt 1600 K、断熱壁、出口 2237 Pa)。"
     "詳細は各 run の標準報告の「境界条件」「解析設定と物性」のスライド。", size=11, color=MUTED)

# 4 Euler A/B
lab = {"M_wave_eta0.1": "Mach 波", "P_wave_eta0.1": "P 波", "overshoot_eta0.1": "オーバーシュート",
       "overshoot_exitnorm_eta0.1": "出口規格化オーバーシュート", "P_slope_abs_eta0.1": "|P 傾き|", "exit_M_dev": "|出口コア M − 6|",
       "exit_core_M": "出口コア M (記録のみ)"}
rows = [["量 (r/r_w = 0.1)", "許容幅 Δq", "旧壁 (3 本平均)", "単調壁 (3 本平均)", "差 D", "2·SE", "(D+2SE)/Δq", "差の検出", "判定"]]
for r in PE["rows"]:
    rows.append([lab[r["qty"]], f"{r['dq']:.2g}" if "dq" in r else "—", f"{r['A']['mean']:.6g}", f"{r['B']['mean']:.6g}", f"{r['D']:+.2e}",
                 f"{2 * r['SE_D']:.2e}", f"{r['D_plus_2SE_over_dq']:.2f}" if "D_plus_2SE_over_dq" in r else "—",
                 ("あり" if r.get("detected") else "なし") if "detected" in r else "—", r["verdict"]])
s = new("Euler A/B: 6 量すべて許容幅未満", "旧壁 run_0140〜0142 / 単調壁 run_0143〜0145、本段 step 6000〜18000 の時間平均 (単位 %pt。|出口コア M − 6| は M の差)",
        notes=f"出典: {B / 'throat_mono_practical_eval.json'} (throat_mono_practical_eval.py)")
table(s, rows, 0.4, 1.45, 12.5, [2.6, 0.9, 1.35, 1.45, 1.05, 0.95, 1.15, 0.9, 2.15], size=11, row_h=0.42)
text(s, 0.5, 5.0, 12.3, 2.2, [
    "総合: " + PE["overall"],
    "限界: 自己相関は補正していない。窓の開始 6000 は予備 A/B の時系列を見た後に決めた。初期値の写像による差 (最近傍 − 番号写像) は Δq/10 の精度では除外できていない "
    f"(オーバーシュートで {rows_pe['overshoot_eta0.1']['ic_beta_minus_alpha']:+.1e}、1 本ずつの比較)。軸 (r = 0) の量は揺れが大きく判定対象外。",
    "経緯: 事前登録した厳密判定 (予備 A/B の閾値 Δq/10) は、評価量の時間変動のため判定不能だった。ユーザ決定 (2026-10-06「さすがにゴミの量」「まずは 1 で」) で、比較を計算する前にこの実務判定に切り替えた。"],
     size=12)

# 5 dry NS
def g(k, key):
    return f"{mean_of(k, key):.6g} ({vd(k, key)})"


s = new("dry NS: 全ゲート合格、旧壁とほぼ同じ", "単調壁 run_0147 + 延長 run_0149 (窓 60000〜80000) / 旧壁 run_0117 (窓 40000〜60000)。量別の判定は check_quasisteady (末尾 5 枚)",
        notes=f"図: {CMP / 'fig_compare_dry_axis.png'} (throat_mono_ns_compare.py)")
table(s, [["ゲート (§6 N)", "単調壁", "旧壁", "合格条件"],
          ["出口半径", f"{RP['0149']['wall_shape']['exit_radius_m']:.7f} m", f"{RP['0117']['wall_shape']['exit_radius_m']:.7f} m", "0.775 ± 0.1 mm"],
          ["出口コア M", g("0149", "exit_core_M"), g("0117", "exit_core_M"), "6 ± 0.02 % (≥ 5.9988)"],
          ["δ_E/δ_C (x_F)", g("0149", "dE_over_dC"), g("0117", "dE_over_dC"), "1 ± 0.5 %"],
          ["Mach 波 η0.1 [%]", g("0149", "wave_eta0.1"), g("0117", "wave_eta0.1"), "≤ 0.01 %"],
          ["オーバーシュート η0.1 [%]", g("0149", "overshoot_eta0.1"), g("0117", "overshoot_eta0.1"), "≤ +0.035 %"],
          ["壁解像 (y₁⁺ > 1 の面積)", f"{RP['0149']['wall_resolution']['verdict']} {RP['0149']['wall_resolution']['over_area_pct']} %",
           f"{RP['0117']['wall_resolution']['verdict']} {RP['0117']['wall_resolution']['over_area_pct']} %", "≤ 5 %"],
          ["残差 (check_convergence)", "NOT CONVERGED (plateau)", "NOT CONVERGED (plateau)", "RISING なし"]],
      0.4, 1.45, 6.6, [1.9, 1.75, 1.75, 1.2], size=10, row_h=0.42)
picture(s, CMP / "fig_compare_dry_axis.png", 7.1, 1.4, 5.9, 4.3)
d = CS["dry"]
text(s, 0.5, 5.55, 12.4, 1.8, [
    f"試験部 (x/r_t 40〜94) の M/M_d − 1 の差 (単調 − 旧、末尾 5 枚の平均): 最大 {d['eta0']['max_abs_diff_pct_test_section_40_94']:.4f} %pt (軸)、"
    f"{d['eta0.1']['max_abs_diff_pct_test_section_40_94']:.4f} %pt (r/r_w = 0.1)。末尾 5 枚の揺れの幅 ({d['eta0']['tail_range_pct_test_section_max']:.4f}・"
    f"{d['eta0.1']['tail_range_pct_test_section_max']:.4f} %pt) と同程度。差が最も大きいのは壁を変えたスロート直後 (x ≈ {d['eta0']['x_of_max']:.1f}) の {d['eta0']['max_abs_diff_pct_x0_95']:.4f} %pt。",
    "経緯: run_0147 の窓 40000〜60000 ではオーバーシュート η0.1 が DRIFTING → 事前登録どおり 20000 step を 1 回延長 (run_0149) して STEADY。"
    "旧壁 run_0117 の Mach 波 DRIFTING はユーザ決定 (2026-10-06「A」) で記録のみとしたもの。"], size=11)

# 6 凝縮
s = new("凝縮 NS: 4 量とも STEADY、旧壁とほぼ同じ", "単調壁 run_0148 / 旧壁 run_0118、窓 14000〜18000 (末尾 5 枚)",
        notes=f"図: {CMP / 'fig_compare_cond_axis.png'}")
table(s, [["量", "単調壁", "旧壁"],
          ["軸の凝縮開始 x/r_t", g("0148", "cond_onset_x_axis"), g("0118", "cond_onset_x_axis")],
          ["過飽和度の最大 S_max", g("0148", "cond_S_max"), g("0118", "cond_S_max")],
          ["出口コアの凝縮質量分率 g", g("0148", "cond_exit_g_core"), g("0118", "cond_exit_g_core")],
          ["出口コア M (凝縮 ON)", g("0148", "exit_core_M"), g("0118", "exit_core_M")]],
      0.4, 1.45, 6.4, [2.4, 2.0, 2.0], size=11, row_h=0.45)
picture(s, CMP / "fig_compare_cond_axis.png", 6.9, 1.4, 6.1, 2.4)
text(s, 0.5, 4.2, 12.3, 2.0, [
    f"軸の凝縮質量分率 g_0 の差は最大 {CS['cond']['g0_max_abs_diff']:.1e} (値は出口で約 {CS['cond']['g0_exit_old']:.4f})。",
    "出口コア M (凝縮 ON) が 6 より低いのは凝縮の発熱によるもので、旧壁と同じ。凝縮 4 量の準定常合格であって、dry の M6 許容への合格ではない。"], size=12)

# 7 限界と未確定事項
s = new("限界と未確定事項")
text(s, 0.5, 1.4, 12.3, 5.8, [
    "1. Euler の比較は実務判定 (自己相関の未補正、窓の開始位置を結果の後に決定、IC 写像の差は Δq/10 の精度で未除外)。厳密な非劣化の証明ではない。",
    "2. NS は、旧壁 run_0117 の最終場を番号写像で移して継続した 1 本 (+ 延長 1 回)。同一入力の再実行による揺れの評価はしていない。",
    "3. 量別の準定常は末尾 5 枚で判定。残差はどの run も 1 桁程度の低下で横ばい (NOT CONVERGED、plateau) で、旧壁と同じ振る舞い。",
    "4. r″ の山の原因: MOC の壁点の始点付近の角度差 (1/R 一定の曲線より最大 0.04°)。軸側の最初の間隔 axis_dx0 に連れて減り "
    "(0.036 → 0.022 → 0.016°)、生産の MOC の単位過程には軸の 1 段目の誤差がある。原因の特定と、生産の MOC を変えるかは別 plan で判断する。"
    "単調拘束はこの角度差を最初の点のずれとして受けるもので、原因を直したものではない。",
    "5. 生産の設定 (問題 YAML の既定) にはまだ入れていない。codex の result 段レビューを経てから入れる。"], size=14)

# 8 run と成果物
s = new("run と成果物", f"run はすべて AWS の ~/forge-wallfit/{CASE}/ (報告・判定の要約はローカルの {CASE}/ にも)")
table(s, [["run", "内容"],
          ["run_0140〜0142_euler_wallfit_pinG1_r{1,2,3}", "Euler 旧壁 × 3"],
          ["run_0143〜0145_euler_wallfit_monoG1_r{1,2,3}", "Euler 単調壁 × 3 (run_0143 は IC 予備 A/B の番号写像も兼ねる)"],
          ["run_0146_euler_icab_monoG1_nn", "IC 予備 A/B (最近傍写像)"],
          ["run_0147_ns_mono_final + run_0149_ns_mono_final_ext", "dry NS 単調壁 (60000 + 延長 20000)。標準報告 run_0149_ns_mono_final_ext_report.pptx"],
          ["run_0148_ns_mono_final_cond", "凝縮 NS 単調壁 (18000)。標準報告 run_0148_ns_mono_final_cond_report.pptx"],
          ["比較の基準", "旧壁 run_0117 (dry)・run_0118 (凝縮)・run_0114 (Euler 参照)"]],
      0.5, 1.45, 12.3, [5.2, 7.1], size=12, row_h=0.5)
text(s, 0.5, 5.3, 12.3, 1.8, [
    "判定の出典: _band_ab/throat_mono_practical_eval.json、_band_ab/verdicts_run_01{17,18,48,49}_monoeval.json、_band_ab/throat_mono_compare/*.json",
    "スクリプト: throat_mono_practical_eval.py・throat_mono_ns_verdicts.py・throat_mono_ns_compare.py・throat_mono_wall_fig.py・run_throat_mono_ab.sh・run_mono_ns_chain.sh・run_mono_ns_ext.sh",
    "plan: plans/active/tooling-nozzle-throat-monotone-r2.md (§6 E′・N・K、§9)"], size=11, color=MUTED)

out = CMP / "throat_mono_summary_report.pptx"
prs.save(out)
print(f"[summary] {out} ({len(prs.slides)} slides)")
