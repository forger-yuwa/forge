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
- 量の種類: residual・mass_balance・exit_core_mach・line_mach・exit_flow_angle・wall_heat_nodes (詳細は settle.py の docstring)。
  各量に verdict・使った値・不合格の理由の分類 (閾値の超過 / データ不足 / 未設定 / 非有限) を書く。全体は AND
  (DIVERGED > INDETERMINATE > NOT_SETTLED > SETTLED)。判定不能は合格にしない。
- 結果は `<run>/SETTLE_VERDICT.json` (基準の中身と sha256・起動の識別・格子・判定区間・終端・使ったファイルの一覧と sha256・各量)。
  閾値に `provisional: true` があれば VERDICT の行に「閾値は暫定 (plan §6 V3 で決める)」と出す。

使い方:
  python3 check_settled.py <run_dir> [--end-step N] [--scan N] [--json] [--no-write | --out PATH]
  --end-step N: N までのデータで判定する (省略時は STOPPED_AT の step、無ければ入力のデータの最後)。
  --scan N: 格子点のうち N の倍数ごとに「その時点まで」で判定し、最初に SETTLED になった step を報告する (N は格子の間隔の倍数)。
  --json: 結果の JSON を標準出力に出す。
  --bind: 本段の起動との照合をしてから判定する (run_case.sh の後処理が使う; plan §4.9 の 2)。`settle.select_judge` が
          本段 (main) と照合し、食い違い・記録の欠けがあれば判定せず INDETERMINATE (SETTLE_VERDICT.json も書かない)。前の段・本段でない起動
          (check_convergence.py の対象) も INDETERMINATE。手で組んだ run (基準が本段の開始の記録を指さない) はそのまま判定する。
  python3 check_settled.py --select <run_dir>
          判定の道具の選択だけ (1 行目 check_settled.py / check_convergence.py、2 行目 理由、続いて照合の問題; 終了コード 0)。
          基準のファイルの無い run・前の段・本段でない起動は check_convergence.py、それ以外は check_settled.py (判定不能でも旧の道具へ落とさない)。

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
    m.add_argument("--bind", action="store_true", help="本段の起動との照合をしてから判定する (run_case.sh の後処理)")
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
    if a.bind:
        sel = settle.select_judge(run)
        if sel["tool"] != settle.TOOL_NEW:
            print(f"この起動は {settle.TOOL_NEW} の対象でない: {sel['why']} (判定は {sel['tool']})")
            print("VERDICT: INDETERMINATE (本段の起動でない)")
            return settle.EXIT_CODE["INDETERMINATE"]
        if sel["mode"] == "main":
            if sel["problems"]:
                print(f"本段の起動との照合の不一致 ({sel['why']}; 判定しない):")
                for p in sel["problems"]:
                    print(f"  {p}")
                print("VERDICT: INDETERMINATE (本段の起動との照合の不一致)")
                return settle.EXIT_CODE["INDETERMINATE"]
            b = sel["binding"]
            print(f"照合: {b['record']} の本段の設定 (fnv {b['main_config_fnv']}) = 今の {b['config_file']}"
                  + (f" = 最後の起動の cfg_fnv" if b["launch_cfg_fnv"] is not None else " (起動の記録なし)")
                  + f"、基準の sha256 {b['criteria_sha256'][:16]}… = 記録")
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
    _print(res)
    if sc is not None:
        print(f"最初に SETTLED になった step (--scan {a.scan}): {sc['first_settled']}")
    if not a.no_write:
        p = settle.write_verdict(run, res, data, a.out)
        print(f"結果: {p}")
    if a.json:
        print(json.dumps(settle._jsonable(res), ensure_ascii=False, indent=1))
    print(settle.verdict_line(res))
    return settle.EXIT_CODE[res["verdict"]]


if __name__ == "__main__":
    sys.exit(main())
