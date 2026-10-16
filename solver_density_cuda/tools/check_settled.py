#!/usr/bin/env python3
r"""収束の判定の道具 — 「残差が準定常・流量の収支が取れる・目的の量が準定常」の AND
(plan `plans/active/tooling-settle-convergence-and-monitor.md` §4.2・§4.3・§4.8、§5.1 #2〜#4。判定の中身は共通の部品 `settle.py`)。

**契約** (低下の桁数は見ない):
- 入力は run の `settle_criteria.yaml` だけ。無い・目的の量が無い → `INDETERMINATE (評価対象が未設定)` (置き場所と量の種類の一覧を出す)。
  知らない種類・キー・閾値の欠け → 基準のファイルの誤りとして INDETERMINATE (黙って既定値で埋めない)。ひな形は
  `settle_criteria.template.yaml`、設計チェーンの run は `design/forge_design/evaluate/settle_presets.py` が問題 YAML の
  `evaluate.settle` から作る。
- 判定は通常の出力の格子 (`time.outStepInterval` の倍数の step) の点だけで行う。窓 W は端点込み (W 20,000・間隔 2,000 → 11 点)、
  残差の上昇は直前の窓 (T−2W〜T−W) と比べる。格子点 s の観測 = 残差の行 s−1・流量の行 s・`res_nodes_s.h5`・壁の h5。
  窓の格子点のどれかが欠ければ判定不能 (データ不足)。本段の開始の場には残差の行も節点・流量の出力も無いので、最初に判定できる
  終端は開始 + 2W + 間隔。協調停止などで終端が格子に乗らなければ、その観測は別の検査 (非有限・直前の格子からの差・窓に加えた許容)
  だけに使う。
- 量の種類: residual・mass_balance・exit_core_mach・line_mach・exit_flow_angle・wall_heat_nodes と、残差を止める判定から外した基準
  (wt_nozzle_v3) の残差の非有限の見張り residual_finite (判定の条件ではない; 非有限 → DIVERGED) (詳細は settle.py の docstring)。
  各量に verdict・使った値・不合格の理由の分類 (閾値の超過 / データ不足 / 未設定 / 非有限) を書く。全体は AND
  (DIVERGED > INDETERMINATE > NOT_SETTLED > SETTLED)。判定不能は合格にしない。
- 結果は `<run>/SETTLE_VERDICT.json` (基準の中身と sha256・起動の識別・格子・判定区間・終端・使ったファイルの一覧と sha256・各量)。
  閾値に `provisional: true` があれば VERDICT の行に「閾値は暫定 (plan §6 V3 で決める)」と出す。

使い方:
  python3 check_settled.py <run_dir> [--end-step N] [--scan N] [--json] [--no-write | --out PATH] [--no-bind]
  --end-step N: N までのデータで判定する (省略時は STOPPED_AT の step、無ければ入力のデータの最後)。
  --scan N: 格子点のうち N の倍数ごとに「その時点まで」で判定し、最初に SETTLED になった step を報告する (N は格子の間隔の倍数)。
  --json: 結果の JSON を標準出力に出す。
  **本段の起動との照合は既定** (plan §4.9 の 2、§5.1 #7c で既定にした — result 段レビュー M1: 引数なしの通常のコマンドが照合を省いていた)。
          `settle.select_judge` が道具を選び、本段の記録 (MAIN_STAGE_STARTED) のある run は記録・今の設定・最後の起動・基準の sha256 と
          照合する。食い違い・記録の欠け (起動の記録が無い、記録が無く前の段とも確かめられない) があれば判定せず INDETERMINATE
          (SETTLE_VERDICT.json も書かない)。**基準が無効 (`mode invalid`: 誤り・評価対象が未設定・grid.segment の stage が main でない)
          のときも mode によらず判定せず INDETERMINATE** (§5.1 #9、result 段レビュー 2 回目 M1: 以前は stage を書き換えた基準で判定して
          SETTLED を返した)。前の段と確かめられた run (check_convergence.py の対象) も INDETERMINATE。手で組んだ run
          (基準が本段の開始の記録を指さず、記録も無い) と基準の無い run はそのまま判定する (後者は「評価対象が未設定」)。
  --bind: 互換のために残す (既定と同じ)。run_case.sh・停止規則の後処理・成果物が付けて呼ぶ。
  --no-bind: 照合を外して判定する (調査用)。外したことを VERDICT の行と結果の JSON (`binding.skipped`) に書く。
  python3 check_settled.py --select <run_dir>
          判定の道具の選択だけ (1 行目 check_settled.py / check_convergence.py、2 行目 理由、続いて照合の問題; 終了コード 0)。
          基準のファイルの無い run と前の段と確かめられた run は check_convergence.py、それ以外は check_settled.py (照合の不一致・
          記録の欠け・判定不能でも旧の道具へ落とさない)。

終了コード: 0 SETTLED / 1 NOT_SETTLED / 2 DIVERGED / 3 INDETERMINATE (基準のファイルの誤り・引数の誤り・照合の不一致も 3; argparse の既定の
2 は DIVERGED と紛れるので使わない)。
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import settle  # noqa: E402


class _Parser(argparse.ArgumentParser):
    def error(self, message):  # 引数の誤りは INDETERMINATE (3) にする (既定の 2 は DIVERGED と紛れる)
        self.print_usage(sys.stderr)
        print(f"{self.prog}: 引数の誤り: {message}", file=sys.stderr)
        print("VERDICT: INDETERMINATE (引数の誤り)")
        sys.exit(settle.EXIT_CODE["INDETERMINATE"])


def _print(res: dict) -> None:
    run = res.get("run_dir", "?")
    T = res.get("terminal_step")
    hdr = f"=== {run}  終端 {T if T is not None else '-'}"
    if res.get("stop_step") is not None:
        hdr += f" (停止の観測 {res['stop_step']} は別の検査だけ)"
    hdr += f"  end_step {res.get('end_step')} [{res.get('end_source', '-')}]"
    print(hdr + " ===")
    if res.get("reason"):
        print(f"  {res['reason']}")
    for q in res.get("quantities", []):
        print(f"  [{q['verdict']:13s}] {q['name']} ({q['kind']}, 入力 {q['input']})")
        for c in q.get("checks", []):
            print(f"      {'ok ' if c['ok'] else 'NG '} {c['check']}: {settle._fmt(c['value'])} (許容 {settle._fmt(c['limit'])})"
                  + (f"  {c['detail']}" if c.get("detail") else ""))
        for r in q.get("reasons", []):
            if r["category"] != settle.CAT_THRESHOLD:
                print(f"      [{r['label']}] {r['text']}")
        if q["kind"] == "residual":
            for k, c in (q.get("values", {}).get("columns") or {}).items():
                if c.get("status") != "judged":
                    print(f"      列 {k}: {c.get('status')}")
    for e in res.get("excluded", []) or []:
        print(f"  外した量: {e['kind']} — {e['reason']}")
    for p in res.get("preconditions", []) or []:
        if not p.get("satisfied"):
            print(f"  前提が満たされていない: {p['quantity']} — {p['requires']}")


def _rejection(run: Path, sel: dict) -> tuple[str, str | None]:
    """判定に進まない選択 (mode invalid か照合の問題) の VERDICT の括弧の文言と、添える案内 (無ければ None)。
    invalid の中身は基準を読み直して分ける: 誤り → 「基準のファイルの誤り」(判定の経路の文言と同じ)、評価対象が未設定 → 置き場所と
    量の種類の一覧 (`settle.load_criteria` の message; §4.2 の契約)、それ以外 (grid.segment の stage が main でない等) → 「基準が本段の
    判定に使えない」。"""
    if sel["mode"] == "main":
        return "本段の起動との照合の不一致", None
    if sel["mode"] == "unconfirmed":
        return "本段の開始の記録の欠け", None
    if sel["mode"] != "invalid":
        return f"照合の問題 (mode {sel['mode']})", None
    try:
        lc = settle.load_criteria(run)
    except settle.CriteriaError as e:
        return "基準のファイルの誤り", f"ERROR: 基準のファイルの誤り ({run / settle.CRITERIA_FILE}): {e}"
    if lc["status"] != "ok":
        return "評価対象が未設定", lc["message"]
    return "基準が本段の判定に使えない", None


def main(argv=None) -> int:
    ap = _Parser(description=__doc__.split("\n\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("--end-step", type=int, default=None)
    ap.add_argument("--scan", type=int, default=0)
    ap.add_argument("--json", action="store_true")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--no-write", action="store_true", help="SETTLE_VERDICT.json を書かない")
    g.add_argument("--out", default=None, help="結果の JSON の書き先 (既定 <run>/SETTLE_VERDICT.json)")
    m = ap.add_mutually_exclusive_group()
    m.add_argument("--select", action="store_true", help="判定の道具の選択だけを出す (1 行目 道具、2 行目 理由; 終了コード 0)")
    m.add_argument("--bind", action="store_true", help="本段の起動との照合をしてから判定する (既定と同じ; 互換のため残す)")
    m.add_argument("--no-bind", action="store_true", help="照合を外して判定する (調査用; 外したことを VERDICT と結果に書く)")
    a = ap.parse_args(argv)
    run = Path(a.run_dir)
    if a.select:
        sel = settle.select_judge(run)
        print(sel["tool"])
        print(f"{sel['why']} [mode {sel['mode']}]")
        for p in sel["problems"]:
            print(f"照合の問題: {p}")
        return 0
    if not run.is_dir():
        print(f"ERROR: {run} はディレクトリでない")
        print("VERDICT: INDETERMINATE (データ不足)")
        return settle.EXIT_CODE["INDETERMINATE"]
    binding = {"skipped": True, "by": "--no-bind"} if a.no_bind else None
    if a.no_bind:
        print("注意: --no-bind — 本段の起動との照合を外して判定する (この VERDICT は照合の裏付けを持たない)")
    else:
        sel = settle.select_judge(run)
        binding = {"skipped": False, **{k: sel.get(k) for k in ("tool", "mode", "why", "problems", "binding")}}
        if sel["mode"] == "legacy":
            # 基準のファイルが無い run: 判定は「評価対象が未設定」(置き場所と量の種類の一覧を出す; 旧の道具の対象であることも書く)
            print(f"照合: {sel['why']} (この run の判定の道具は {sel['tool']}) [mode legacy]")
        elif sel["tool"] != settle.TOOL_NEW:
            print(f"この起動は {settle.TOOL_NEW} の対象でない: {sel['why']} (判定は {sel['tool']})")
            print("VERDICT: INDETERMINATE (本段の起動でない)")
            return settle.EXIT_CODE["INDETERMINATE"]
        elif sel["mode"] == "invalid" or sel["problems"]:
            # 基準が無効 (mode invalid) か照合に問題があれば、mode によらず判定に進まない (result 段レビュー 2 回目 M1、plan §5.1 #9:
            # 以前は main・unconfirmed だけを拒否し、基準の grid.segment.stage を main 以外に書き換えた run (invalid; 基準の sha256 と
            # 本段の記録の照合も通らない経路) を判定して SETTLED・終了コード 0 を返していた)
            what, msg = _rejection(run, sel)
            print(f"{what} ({sel['why']} [mode {sel['mode']}]; 判定しない):")
            for p in sel["problems"]:
                print(f"  {p}")
            if msg:
                print(msg)
            print(f"VERDICT: INDETERMINATE ({what})")
            return settle.EXIT_CODE["INDETERMINATE"]
        elif sel["mode"] == "main":
            b = sel["binding"]
            print(f"照合: {b['record']} の本段の設定 (fnv {b['main_config_fnv']}) = 今の {b['config_file']} = 最後の起動の cfg_fnv"
                  f"、基準の sha256 {b['criteria_sha256'][:16]}… = 記録")
        else:
            print(f"照合: {sel['why']} [mode {sel['mode']}]")
    try:
        loaded = settle.load_criteria(run)
        data = settle.RunData(run, loaded["criteria"]) if loaded["status"] == "ok" else None
        res = settle.judge(run, a.end_step, data=data, loaded=loaded)
        sc = settle.scan(run, a.scan, end_step=a.end_step) if a.scan else None
    except settle.CriteriaError as e:
        print(f"ERROR: 基準のファイルの誤り ({run / settle.CRITERIA_FILE}): {e}")
        print("VERDICT: INDETERMINATE (基準のファイルの誤り)")
        return settle.EXIT_CODE["INDETERMINATE"]
    if sc is not None:
        res["scan"] = sc
    res["binding"] = binding
    _print(res)
    if sc is not None:
        print(f"最初に SETTLED になった step (--scan {a.scan}): {sc['first_settled']}")
    if not a.no_write:
        p = settle.write_verdict(run, res, data, a.out)
        print(f"結果: {p}")
    if a.json:
        print(json.dumps(settle._jsonable(res), ensure_ascii=False, indent=1))
    print(settle.verdict_line(res) + (" (本段の起動との照合を外した: --no-bind)" if a.no_bind else ""))
    return settle.EXIT_CODE[res["verdict"]]


if __name__ == "__main__":
    sys.exit(main())
