"""case/45 の生産の問題 YAML を NS の N2 の問題から作る (2026-10-07 ユーザ決定「生産採用: 上流を多項式に・MOC の新方式・全域 1 本の
B スプライン壁」、plans/accepted/tooling-nozzle-upstream-poly-and-throat-sizing.md §9)。

  problem_d155_ns_n012_N2.yaml      → problem_d155_ns_prod.yaml
  problem_d155_ns_n012_N2_cond.yaml → problem_d155_ns_prod_cond.yaml

N2 との差は name・冒頭のコメント・`geometry.physical_wall_repr: single_bspline` の 1 行だけ (全域 1 本の plan の W3・W4 で
ソルバ入力がビット同一 = CFD のやり直しは要らない)。較正値・r_t・k_f・MOC のキー (analytic・converge を明示)・格子・設定は N2 のまま。
生産の NS は run_0167 + 延長 run_0179、凝縮は run_0170 (どちらも N2)。

usage: python3 make_prod_problems.py [--check]   (--check は書かずに、今のファイルが生成結果と一致するかだけを見る)
"""
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAIRS = (("problem_d155_ns_n012_N2.yaml", "problem_d155_ns_prod.yaml"),
         ("problem_d155_ns_n012_N2_cond.yaml", "problem_d155_ns_prod_cond.yaml"))
UP = "  pw_upstream: poly"


def build(src: str, dst: str) -> str:
    txt = (HERE / src).read_text()
    sha = hashlib.sha256(txt.encode()).hexdigest()
    lines = txt.splitlines(keepends=True)
    names = [i for i, l in enumerate(lines) if l.startswith("name: ")]
    ups = [i for i, l in enumerate(lines) if l.startswith(UP)]
    if len(names) != 1 or len(ups) != 1 or "physical_wall_repr" in txt:
        raise SystemExit(f"{src}: name・pw_upstream の行が 1 つずつでない、または physical_wall_repr が既にある — 止める")
    lines[names[0]] = f"name: {Path(dst).stem}\n"
    lines.insert(ups[0] + 1, "  physical_wall_repr: single_bspline   # 全域 1 本の 5 次 B-spline (ソルバ入力は区分表現とビット同一; plan tooling-nozzle-wall-single-bspline)\n")
    head = (f"# case/45 の生産の問題 (2026-10-07 ユーザ決定「生産採用: 上流を多項式に・MOC の新方式・全域 1 本の B スプライン壁」)。\n"
            f"#   生成: make_prod_problems.py (手で編集しない)。元: {src} (sha256 {sha[:16]}…)。元との差は name とこの冒頭と physical_wall_repr の行だけ。\n"
            f"#   生産の NS は N2 の run_0167 + 延長 run_0179、凝縮は run_0170 (plans/accepted/tooling-nozzle-upstream-poly-and-throat-sizing.md §9)。\n"
            f"# ---- 以下は元の YAML ----\n")
    return head + "".join(lines)


def main():
    check = "--check" in sys.argv[1:]
    bad = 0
    for src, dst in PAIRS:
        out = build(src, dst)
        p = HERE / dst
        if check:
            ok = p.is_file() and p.read_text() == out
            print(f"{dst}: {'一致' if ok else '不一致'}")
            bad += 0 if ok else 1
        else:
            if p.exists() and p.read_text() != out:
                raise SystemExit(f"{dst} が既にあり内容が違う (上書きしない)")
            p.write_text(out)
            print(f"wrote {dst}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
