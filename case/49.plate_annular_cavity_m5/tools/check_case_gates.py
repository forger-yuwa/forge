#!/usr/bin/env python3
"""1 つの run について**報告してよいか**を全ゲートで判定する (AGENTS.md の各ルールの束ね)。

これが通らない run の数値は報告しない。個別ツールを呼び出して VERDICT を集約する:

| # | 何を | ツール | 合格条件 |
| --- | --- | --- | --- |
| 1 | NaN/発散 | 残差 csv + 最終 res | 非有限ゼロ |
| 2 | 残差の収束 | `check_convergence.py --segment` | 判定区間で PASS |
| 3 | 結論量の準定常 | `check_cavity_steady.py` | 全量 STEADY |
| 4 | 壁解像 | `check_wall_resolution.py` | y1+ 超過面積 <= --yplus-frac (**既定は非ブロッキング**) |
| 5 | 保存性 | `cavity_eval.py` | 正味/片道 <= --mass-tol, CV 収支 <= --budget-tol |

usage:
  python3 tools/check_case_gates.py <run_dir> [--yplus-frac 2] [--mass-tol 0.015] [--budget-tol 0.05]

**ブロッキングと非ブロッキングを分ける** (2026-09-19): 計算を延長して直るのは 1/2/3/5 だけで、
**4 (壁解像) と メッシュ品質は step を増やしても変わらない**。それらを不合格扱いにすると
自動延長ループが無駄に回るので、既定では「**報告時に必ず添える制約**」として出力し、
終了コードには入れない (`--yplus-blocking` で厳格化できる)。

終了コード 0 = ブロッキングゲート PASS (数値を報告してよい)。1 = 不合格。2 = 判定不能。
"""
import argparse
import glob
import json
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASE = HERE.parent
ROOT = CASE.parents[1]
STOOLS = ROOT / "solver_density_cuda" / "tools"


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def parent_run(rd):
    """`CONTINUED_FROM` が指す引き継ぎ元の run ディレクトリ。"""
    f = rd / "CONTINUED_FROM"
    if not f.exists():
        return None
    src = Path(f.read_text().strip().splitlines()[0])
    par = src.parent
    if not par.is_absolute():
        par = (rd.parent / par).resolve()
    return par if (par / "residual_history.csv").exists() else None


def tail_levels(csv, tail=0.2):
    """残差列ごとの**末尾窓の中央値** (水準の比較用)。"""
    import csv as _csv
    try:
        rows = list(_csv.reader(open(csv)))
    except OSError:
        return {}
    if len(rows) < 3:
        return {}
    hdr = [c.strip() for c in rows[0]]
    idx = [i for i, c in enumerate(hdr) if c.startswith("rms_")]
    body = [r for r in rows[1:] if r and len(r) > max(idx, default=0)]
    if not body:
        return {}
    n = max(2, int(len(body) * tail))
    out = {}
    for i in idx:
        v = []
        for r in body[-n:]:
            try:
                x = abs(float(r[i]))
            except (ValueError, IndexError):
                continue
            if x == x and x != float("inf"):
                v.append(x)
        if v:
            v.sort()
            out[hdr[i]] = v[len(v) // 2]
    return out


BUDGET_KEYS = ("budget_residual_all_alldepth", "budget_residual_alldepth", "budget_residual")


def budget_key(fld):
    """収支残差として読むキーを選ぶ。**全流束 (対流+伝導+粘性仕事) を含むものが正本**。

    2026-09-21 codex result M2: 従来は `budget_residual_alldepth` (= 対流のみ) を優先し、
    伝導・粘性仕事を落とした残差で合否を決めていた。対流だけ 0・全項 20 % の入力でも
    PASS になる。伝導込みのキーが無い JSON は**古い評価版**なので合格にしない。
    """
    return next((k for k in BUDGET_KEYS if k in (fld or {})), None)


def residual_inconclusive(out, rc):
    """継続 run の残差判定が**判定不能**か。

    2026-09-21 codex result M1: 継続 run は収束場から始まるので `rc == 0` は要求できないが
    (`NOT CONVERGED` が正常)、残差列の欠落・末尾窓 0 のような**入力不備**は別物で、
    `RISING`/`DIVERGED` を含まないため合格として通っていた。
    """
    return ("判定不能" in (out or "")) or (rc is not None and rc >= 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--yplus-frac", type=float, default=2.0, help="y1+>1 を許す面積 [%]")
    ap.add_argument("--mass-tol", type=float, default=0.015,
                    help="開口の |正味/片道| (末尾平均) 許容。**これは離散スキームの保存性ではなく"
                         "後処理の面積分精度の指標**なので、本当の保存ゲートは CV エネルギー収支の方。"
                         "**plan の 0.1 %% はこの求積では到達不能** (2026-09-20 実測: 181x30 -2.46 %% / "
                         "361x60 -1.44 %% / 721x120 -1.00 %% / 1441x240 -0.88 %% と約 -0.85 %% へ漸近)。"
                         "既定 1.5 %% は漸近値に余裕を見た値。厳密検算は残作業 #4 (離散流束)")
    # **plan §6.4 の 5 % に合わせる** (2026-09-21 codex result M8)。既定が 0.07 で、
    # plan (5 %) と docstring (0.05) の両方とずれていた。
    ap.add_argument("--budget-tol", type=float, default=0.05, help="CV エネルギー収支 残差 許容")
    ap.add_argument("--yplus-blocking", action="store_true",
                    help="壁解像の不合格でも報告を止める (既定は制約として併記するのみ。"
                         "step を増やしても y1+ は変わらないため)")
    a = ap.parse_args()
    rd = Path(a.run)
    fails, notes, caveats = [], [], []
    py = sys.executable

    print("======== ゲート判定: %s ========" % a.run)

    # --- 1. NaN / 発散 ---
    bad = 0
    for f in list(rd.glob("residual_history*.csv")):
        for ln in f.read_text().splitlines()[1:]:
            if re.search(r"(?i)\b(nan|inf)\b", ln):
                bad += 1
    if list(rd.glob("res_nan_*.h5")):
        bad += 1
    print("[1] NaN/発散       : %s" % ("OK" if bad == 0 else "**FAIL** (%d 箇所 / NaN ダンプ有り)" % bad))
    if bad:
        fails.append("NaN/発散")

    # --- 2. 残差の収束 (判定区間) ---
    # **継続 run (`--main-only`) は「低下桁数」で判定しない**。収束場から始まるので 3 桁落ちることは
    # 原理的に無く、正しい問いは「床に留まっているか (上昇していないか)」である
    # (2026-09-19: これを見落として自動延長 run が必ず不合格になる設計ミスをした)。
    cont = (rd / "CONTINUED_FROM").exists()
    seg = ["--segment"] if (rd / "stage_manifest.json").exists() else []
    rc, out = run([py, str(STOOLS / "check_convergence.py"), str(rd)] + seg)
    head = next((l for l in out.splitlines() if l.startswith("===")), "(出力なし)")
    # **残差列がそもそも出ていない run を合格にしない** (2026-09-19)。
    # インスタンス再起動で forge が 1 step も回らなかった run
    # (`residual_history.csv` 無し / 空) を「継続 run は上昇の有無で判定」の経路が
    # 「上昇が無い」=OK と判定して通していた。判定不能は不合格。
    rh = rd / "residual_history.csv"
    nline = len(rh.read_text().strip().splitlines()) if rh.exists() else 0
    if nline < 2:
        print("[2] 残差の収束     : **判定不能** (residual_history.csv が %s)"
              % ("無い" if nline == 0 else "ヘッダのみ"))
        return 2
    if cont:
        if head == "(出力なし)":                 # 判定行が出ていないのも判定不能
            print("[2] 残差の収束     : **判定不能** (check_convergence が判定行を出さなかった)")
            return 2
        rising = [l for l in out.splitlines() if "RISING" in l or "DIVERGED" in l]
        # **文字列検査だけで通さない** (2026-09-20 codex result M2)。従来は `RISING` /
        # `DIVERGED` が出ていなければ合格で、**終了コードも親の残差水準も見ていなかった**。
        # 継続 run は収束場から始まるので低下桁数では測れないが、
        # **親の末尾水準より悪化していないこと**は測れる。
        worse = []
        par = parent_run(rd)
        if par is not None:
            lv_c = tail_levels(rd / "residual_history.csv")
            lv_p = tail_levels(par / "residual_history.csv")
            for k in sorted(set(lv_c) & set(lv_p)):
                if lv_p[k] > 0 and lv_c[k] > 1.5 * lv_p[k]:
                    worse.append("%s %.2e -> %.2e (x%.2f)" % (k, lv_p[k], lv_c[k], lv_c[k] / lv_p[k]))
        nan_bad = "NaN/Inf present" in out
        # **判定不能・入力不備・異常終了は合格にしない** (2026-09-21 codex result M1)。
        # 継続 run は低下桁数で測れないので `rc == 0` は要求できないが (収束場から
        # 始まるので `check_convergence.py` は NOT CONVERGED を返すのが正常)、
        # **「判定不能」= 残差列の欠落・末尾窓 0 などの入力不備**は別物で、
        # 従来はこれが RISING/DIVERGED を含まないため合格として通っていた。
        incon = residual_inconclusive(out, rc)
        ok2 = (not rising) and (not nan_bad) and (not worse) and (not incon)
        print("[2] 残差の収束     : %s  (継続 run: 低下桁数でなく**上昇の有無 + 親の水準比**で判定)"
              % ("OK" if ok2 else "**FAIL**"))
        print("      %s" % head.strip())
        if par is None:
            print("      [注意] 親 run を特定できず水準比を取れていない (CONTINUED_FROM を確認)")
        if rising:
            print("      上昇している列: %s" % " / ".join(l.split(":")[0].strip() for l in rising))
        if worse:
            print("      親より悪化した列: %s" % " / ".join(worse))
        if incon:
            print("      **判定不能** (入力不備 / check_convergence が異常終了 rc=%d)" % rc)
    else:
        ok2 = (rc == 0)
        print("[2] 残差の収束     : %s" % ("OK" if ok2 else "**FAIL**"))
        print("      %s" % head.strip())
    if not ok2:
        fails.append("残差")
        notes.append(out)

    # --- 3. 結論量の準定常 ---
    rc3, out3 = run([py, str(HERE / "check_cavity_steady.py"), str(rd)])
    v3 = next((l for l in out3.splitlines() if l.startswith("VERDICT")), "")
    ok3 = "STEADY (全量)" in v3
    print("[3] 結論量の準定常 : %s   %s" % ("OK" if ok3 else "**FAIL**", v3.strip()))
    if not ok3:
        fails.append("準定常")
        notes.append(out3)

    # --- 4. 壁解像 ---
    rc4, out4 = run([py, str(STOOLS / "check_wall_resolution.py"), str(rd),
                     "--over-frac", str(a.yplus_frac)])
    v4 = next((l for l in out4.splitlines() if l.startswith("VERDICT")), "")
    tag4 = "OK" if rc4 == 0 else ("**FAIL**" if a.yplus_blocking else "**制約あり**")
    print("[4] 壁解像 y1+     : %s   %s" % (tag4, v4.strip()))
    if rc4 != 0:
        # **step を増やしても変わらない**ので既定ではブロックしない。報告時に添える制約にする。
        worst_lines = [l.strip() for l in out4.splitlines() if "y1+ >" in l]
        for l in worst_lines:
            print("      %s" % l)
        (fails if a.yplus_blocking else caveats).append("壁解像 (y1+ 超過)")

    # --- 5. 保存性 (開口の質量収支と CV エネルギー収支) ---
    j = rd / "cavity_eval.json"
    if not j.exists():
        run([py, str(HERE / "cavity_eval.py"), str(rd)])
    if j.exists():
        d = json.loads(j.read_text())
        # **時間平均で見る**。振動する (キャビティが呼吸する) 流れでは瞬時の正味流束は 0 にならない。
        imb = abs(float(d["field"].get("mdot_imbalance", float("nan"))))
        sc = rd / "cavity_series.csv"
        if sc.exists():
            import csv as _csv
            rows = list(_csv.DictReader(open(sc)))
            vals = [float(r["mdot_imbalance"]) for r in rows if r.get("mdot_imbalance")]
            if len(vals) >= 3:
                n = max(2, int(len(vals) * 0.4))
                imb = abs(sum(vals[-n:]) / n)
        # **欠損はスキップでなく判定不能** (2026-09-20 codex result M4)。
        # 従来は `budget_residual` が JSON に無く、この検査が常に飛んでいた。
        # **`field` の下を見る** (2026-09-20 修正)。トップレベルを見ていたため、
        # `cavity_eval.py` が `field.budget_residual` に書いていても常に「無い」と判定し、
        # 準定常 STEADY の run を 2 回ずつ無駄に延長していた。
        fld = d.get("field", {})
        # **全深さ Σq で正規化した残差**を使う (2026-09-20)。評価面より下の壁入熱で割ると、
        # 非一様壁温では熱い壁の放熱と冷たい壁の吸熱が打ち消して分母が小さくなり、
        # 同じ絶対差 (実測 0.9-3.0 W) でも相対値が跳ねる (mixA: 面下基準 -14.7 % / 全深さ基準 -5.8 %)。
        # **全流束 (対流+伝導+粘性仕事) の残差で判定する** (2026-09-21 codex result M2)。
        # 従来は `budget_residual_alldepth` を優先していたが、これは**対流のみ**の残差で、
        # 伝導・粘性仕事を落としていた (対流だけ 0・全項 20 % の入力でも PASS を再現)。
        # 分母は壁温分布に依らない全深さ Σq のままにする。
        key = budget_key(fld)
        bud = abs(float(fld[key])) if key else None
        if bud is None:
            print("[5] 保存性         : **判定不能** (cavity_eval.json に budget_residual が無い"
                  " — cavity_eval.py を新しい版で回し直すこと)")
            return 2
        if key != "budget_residual_all_alldepth":
            # 伝導・粘性を含む残差が無い JSON は**古い評価版**。合格にしない。
            print("[5] 保存性         : **判定不能** (%s しか無い = 対流のみの残差。"
                  "伝導・粘性込みの budget_residual_all_alldepth が要る"
                  " — cavity_eval.py を新しい版で回し直すこと)" % key)
            return 2
        ok5 = imb <= a.mass_tol
        print("[5] 保存性         : %s   開口の正味/片道 %.3e (許容 %.3g)"
              % ("OK" if ok5 else "**FAIL**", imb, a.mass_tol))
        if not ok5:
            fails.append("質量収支")
        if bud is not None:
            ok5b = bud <= a.budget_tol
            print("                      CV エネルギー収支 残差 %.3f (許容 %.3g) %s"
                  % (bud, a.budget_tol, "OK" if ok5b else "**FAIL**"))
            if not ok5b:
                fails.append("エネルギー収支")
    else:
        print("[5] 保存性         : **判定不能** (cavity_eval.json が作れない)")
        fails.append("保存性(判定不能)")

    print("\nGATES: %s" % ("PASS (報告してよい)" if not fails
                           else "FAIL — " + " / ".join(fails) + " (数値を報告しないこと)"))
    if caveats:
        print("CAVEATS (報告時に必ず添えること): " + " / ".join(caveats))
    if fails and notes:
        print("--- 不合格の詳細 ---")
        for n in notes:
            print("\n".join(n.strip().splitlines()[-12:]))
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
