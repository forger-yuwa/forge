"""plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U4「NS の 3 条件」・plan discretization-moc-axis-limit-and-corrector §6 V5d の
順序 (2026-10-07 諮問 notes/reviews/2026-10-07-euler-grid-switch-plan-diagnose.md) の問題 YAML を作る。

生産の問題 (単調壁) を写し、出口較正 (Euler の E4 の結果) の Md_moc_offset を引数で書き込む。3 条件で較正値・k_f・r_t・NS の格子・
実効の設定は共通 (元の YAML のまま。r_t だけは --r-throat を渡すと 6 本すべてに同じ値を書く)。変えるのは次だけ:
  N0 (対照)   : Md_moc_offset だけ。pw_upstream: ramp (pw_ramp [−11, −6] のまま)、MOC は legacy・fixed2 (既定と同じ値を明示)。
  N1 (U4)     : N0 に pw_upstream: poly (pw_ramp を外す)。
  N2 (V5′)   : N1 に geometry.moc_axis_limit: analytic・moc_corrector: converge。
dry と凝縮 (_cond) の両方を作る。元: problem_d155_ns_finemesh_recal_final_mono.yaml / …_mono_cond.yaml。
--r-throat <m の repr> (任意; plan §6 U4「r_t の解き直し」): 新しい較正値の N0 の問題で CFD 前の solve_rt (prepare_ns と共通の δ_r の経路、
  許容差 1e-5 m、k_f は変えない) により出口半径 0.775 m に解き直した r_t (主セッションが計算して渡す)。spec.r_throat を 3 条件・dry と
  _cond の 6 本すべてに同じ値で書く。省略時は元の YAML の r_t のまま (spec.r_throat の差は許さない)。

検査 (1 つでも破れば何も書かずに止める):
  - 元の YAML が想定どおり (geometry に pw_ramp・pw_upstream: ramp・Md_moc_offset が 1 行ずつ、moc の 2 キーが無い)。
  - 生成した YAML を読み直し (重複キーは拒否)、元との差が name・geometry の Md_moc_offset / pw_upstream / pw_ramp /
    moc_axis_limit / moc_corrector (--r-throat のときは spec.r_throat も) だけで、値が条件どおり。Md_moc_offset と r_throat は
    引数の値とビット一致する float (YAML のトークンを読み直して往復一致を検査)。
  - dry と _cond の差が、元の dry と _cond の差と同じ (name を除く)。
  - runner_axismach の読み取り (_pw_upstream・_moc_keys・load_problem) での実効値が条件どおり (import するだけ、design/ は変えない)。
既存のファイルと中身が同じなら何もしない。違えば止める (--overwrite で上書き。準備済みの run があるときは使わない)。
記録: _band_ab/ns_n012_problems.json (較正値・元と生成物の sha256・実効値・検査)。投入スクリプトはこれで照合する。

usage: python3 make_ns_n012_problems.py --md-offset <値> [--r-throat <m>] [--overwrite] [--out-dir DIR (既定: このディレクトリ)] [--check]
  --check: 書かずに、既存の 6 本と記録が --md-offset・--r-throat の値・今の元の YAML と一致するかだけを調べる (終了コード 0 / 2)。
    --r-throat を省いた --check は「元の r_t のまま」の 6 本を期待する (r_t を書いた 6 本とは一致しない)。
"""
import argparse
import copy
import hashlib
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
TOOLS = HERE.parents[1] / "solver_density_cuda/tools"
DESIGN = HERE.parents[1] / "design"

PLAN_U4 = "plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md §6 U4 (NS の 3 条件)"
PLAN_V5 = "plans/active/discretization-moc-axis-limit-and-corrector.md §6 V5d (順序)・V5′"
BASE = {"dry": "problem_d155_ns_finemesh_recal_final_mono.yaml", "cond": "problem_d155_ns_finemesh_recal_final_mono_cond.yaml"}
PREFIX = "problem_d155_ns_n012"
RECORD = "_band_ab/ns_n012_problems.json"
CONDS = ("N0", "N1", "N2")
# 条件ごとの実効値 (runner_axismach の読み取り結果と照合する)
SPEC = {
    "N0": {"pw_upstream": "ramp", "moc_axis_limit": "legacy", "moc_corrector": "fixed2", "pw_ramp": True,
           "role": "対照: 今の生産 (単調壁・ramp・legacy MOC) で Md_moc_offset だけを新しい較正値に"},
    "N1": {"pw_upstream": "poly", "moc_axis_limit": "legacy", "moc_corrector": "fixed2", "pw_ramp": False,
           "role": "U4: N0 に pw_upstream poly (pw_ramp を外す)。N1 − N0 で上流の多項式化を評価"},
    "N2": {"pw_upstream": "poly", "moc_axis_limit": "analytic", "moc_corrector": "converge", "pw_ramp": False,
           "role": "V5′: N1 に MOC の analytic・converge。N2 − N1 で MOC の軸処理を評価"},
}
# 元との差として許すパス (これ以外の差は止める)
ALLOWED = {("name",), ("geometry", "Md_moc_offset"), ("geometry", "pw_upstream"), ("geometry", "pw_ramp"),
           ("geometry", "moc_axis_limit"), ("geometry", "moc_corrector")}
RT_PATH = ("spec", "r_throat")          # --r-throat のときだけ許す差
R_THROAT_RANGE_M = (0.01, 1.0)          # --r-throat の受け付け範囲 [m] (mm で書いた・桁違いを止めるための粗い枠。物理の判定ではない)
BASE_PW_RAMP = [-11.0, -6.0]
_MISSING = object()


def sha256_file(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def out_name(cond: str, kind: str) -> str:
    return f"{PREFIX}_{cond}{'_cond' if kind == 'cond' else ''}.yaml"


def yaml_load(text: str):
    """重複キーを拒否して読む (solver_density_cuda/tools/yaml_strict.load; 値は yaml.safe_load と同じ)。"""
    if str(TOOLS) not in sys.path:
        sys.path.insert(0, str(TOOLS))
    import yaml_strict
    return yaml_strict.load(text)


def offset_token(value: float) -> str:
    """Md_moc_offset の YAML トークン (float_token)。"""
    return float_token(value, "Md_moc_offset")


def float_token(value: float, what: str) -> str:
    """float の YAML トークン。repr で丸めずに書き、PyYAML (YAML 1.1) が float と読む形にする
    (`1e-05` は小数点が無いと文字列になるので `1.0e-05` にする)。読み直して型と値が往復一致しなければ ValueError。"""
    if not (isinstance(value, float) and math.isfinite(value)):
        raise ValueError(f"{what} は有限の float (受け取った値: {value!r})")
    tok = repr(value)
    if "e" in tok and "." not in tok.split("e")[0]:
        m, e = tok.split("e")
        tok = f"{m}.0e{e}"
    if "e" in tok and not re.search(r"e[-+]", tok):
        tok = tok.replace("e", "e+")
    back = yaml_load(f"v: {tok}\n")["v"]
    if type(back) is not float or back != value:
        raise ValueError(f"トークン {tok!r} を読み直すと {back!r} ({type(back).__name__}) で、{value!r} と一致しない")
    return tok


def parse_offset(text: str) -> float:
    try:
        v = float(text)
    except (TypeError, ValueError):
        raise ValueError(f"--md-offset {text!r} は数でない") from None
    if not math.isfinite(v):
        raise ValueError(f"--md-offset {text!r} は有限でない")
    if abs(v) > 0.01:
        raise ValueError(f"--md-offset {v!r} は |値| ≤ 0.01 の範囲の外 (出口較正の補正量としては大きすぎる。単位・桁を確かめる)")
    return v


def parse_r_throat(text: str) -> float:
    """--r-throat [m]。有限で R_THROAT_RANGE_M の中 (mm で書いた・桁違いを止める)。"""
    try:
        v = float(text)
    except (TypeError, ValueError):
        raise ValueError(f"--r-throat {text!r} は数でない") from None
    lo, hi = R_THROAT_RANGE_M
    if not (math.isfinite(v) and lo < v < hi):
        raise ValueError(f"--r-throat {text!r} は {lo} < r_t < {hi} [m] の外 (単位は m。mm で書いていないか)")
    return v


def diff_paths(a, b, path=()) -> dict:
    """2 つの読み込み結果の差 {パス: (a の値, b の値)}。辺の片側にしか無いキーは _MISSING。"""
    out = {}
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            out.update(diff_paths(a.get(k, _MISSING), b.get(k, _MISSING), path + (k,)))
        return out
    if a is _MISSING or b is _MISSING or type(a) is not type(b) or a != b:
        out[path] = (a, b)
    return out


def _block_lines(lines, block: str):
    """トップレベルのブロック (例 geometry:) の行番号の範囲 [i0, i1)。"""
    i0 = [i for i, l in enumerate(lines) if l.rstrip("\n") == f"{block}:"]
    if len(i0) != 1:
        raise ValueError(f"トップレベルの {block}: が {len(i0)} 行 (1 行であること)")
    i = i0[0] + 1
    while i < len(lines) and (lines[i].startswith((" ", "#")) or not lines[i].strip()):
        i += 1
    return i0[0] + 1, i


def _key_line(lines, lo, hi, key: str, block: str = "geometry") -> int:
    hits = [i for i in range(lo, hi) if re.match(rf"^  {re.escape(key)}:(\s|$)", lines[i])]
    if len(hits) != 1:
        raise ValueError(f"{block} の {key}: の行が {len(hits)} 行 (1 行であること)")
    return hits[0]


def check_base(text: str, kind: str) -> dict:
    """元の YAML が想定どおりか (行と値)。戻り値 = 読み込み結果。"""
    d = yaml_load(text)
    g = d.get("geometry") or {}
    why = []
    if g.get("pw_upstream") != "ramp":
        why.append(f"geometry.pw_upstream が ramp でない ({g.get('pw_upstream')!r})")
    if g.get("pw_ramp") != BASE_PW_RAMP:
        why.append(f"geometry.pw_ramp が {BASE_PW_RAMP} でない ({g.get('pw_ramp')!r})")
    for k in ("moc_axis_limit", "moc_corrector"):
        if k in g:
            why.append(f"geometry.{k} が既にある ({g[k]!r}; 元は既定の legacy・fixed2 のはず)")
    if not isinstance(g.get("Md_moc_offset"), float):
        why.append(f"geometry.Md_moc_offset が float でない ({g.get('Md_moc_offset')!r})")
    if str(g.get("wall_repr")) != "joint":
        why.append(f"geometry.wall_repr が joint でない ({g.get('wall_repr')!r})")
    if not isinstance((d.get("spec") or {}).get("r_throat"), float):
        why.append(f"spec.r_throat が float でない ({(d.get('spec') or {}).get('r_throat')!r})")
    if why:
        raise ValueError(f"元の YAML ({BASE[kind]}) が想定と違う: " + "; ".join(why))
    lines = text.splitlines(keepends=True)
    lo, hi = _block_lines(lines, "geometry")
    for k in ("pw_ramp", "pw_upstream", "Md_moc_offset"):
        _key_line(lines, lo, hi, k)
    slo, shi = _block_lines(lines, "spec")
    _key_line(lines, slo, shi, "r_throat", "spec")
    return d


def render(base_text: str, base_name: str, base_sha: str, cond: str, kind: str, value: float, tok: str,
           rt: float | None = None, rt_tok: str | None = None) -> str:
    """元の YAML のテキストを行単位で書き換える (コメントと書式は残す)。rt (--r-throat) があれば spec.r_throat の行も書き換える。"""
    sp = SPEC[cond]
    lines = base_text.splitlines(keepends=True)
    if rt is not None:
        slo, shi = _block_lines(lines, "spec")
        i_rt = _key_line(lines, slo, shi, "r_throat", "spec")
        lines[i_rt] = (f"  r_throat: {rt_tok}      # [m] N0 の問題で CFD 前の solve_rt (prepare_ns と共通の δ_r の経路、許容差 1e-5 m、k_f は不変) で"
                       f"出口半径 0.775 m に解き直した値 (3 条件で共通; plan tooling-nozzle-upstream-poly-and-throat-sizing §6 U4「r_t の解き直し」。"
                       f"元の YAML の値 {yaml_load(lines[i_rt])['r_throat']!r})\n")
    lo, hi = _block_lines(lines, "geometry")
    i_md = _key_line(lines, lo, hi, "Md_moc_offset")
    i_up = _key_line(lines, lo, hi, "pw_upstream")
    i_ramp = _key_line(lines, lo, hi, "pw_ramp")
    new = list(lines)
    new[i_md] = (f"  Md_moc_offset: {tok}       # 出口較正 (Euler の E4 の結果; make_ns_n012_problems.py --md-offset で書き込み。"
                 f"元の YAML の値 {yaml_load(lines[i_md])['Md_moc_offset']!r} は G1 の Euler の較正)\n")
    moc = [f"  moc_axis_limit: {sp['moc_axis_limit']}"
           + ("       # 既定と同じ値を明示 (コードの既定が変わっても N0・N1 の MOC を固定する)\n" if sp["moc_axis_limit"] == "legacy"
              else "     # 軸上の端点の sinθ/r を解析極限 θ_r に (plan discretization-moc-axis-limit-and-corrector §4.1)\n"),
           f"  moc_corrector: {sp['moc_corrector']}"
           + ("        # 既定と同じ値を明示\n" if sp["moc_corrector"] == "fixed2"
              else "      # 予測修正を収束まで (同 plan §4.2)\n")]
    new[i_md] = new[i_md] + "".join(moc)
    if sp["pw_upstream"] == "poly":
        new[i_up] = ("  pw_upstream: poly           # 配管〜設計スロートの 5 次多項式、上流に δ_r を足さない "
                     "(plan tooling-nozzle-upstream-poly-and-throat-sizing §4.1。pw_ramp は併記不可なので外した)\n")
        new[i_ramp] = ""
    i_name = [i for i, l in enumerate(new) if re.match(r"^name:\s", l)]
    if len(i_name) != 1:
        raise ValueError(f"トップレベルの name: が {len(i_name)} 行")
    new[i_name[0]] = f"name: {out_name(cond, kind)[:-5]}\n"
    head = [f"# {PLAN_U4}・{PLAN_V5}: 条件 {cond} ({'凝縮 ON' if kind == 'cond' else 'dry'})。\n",
            f"#   {sp['role']}。\n",
            f"#   生成: make_ns_n012_problems.py --md-offset {tok}{f' --r-throat {rt_tok}' if rt is not None else ''} (手で編集しない)。"
            f"元: {base_name} (sha256 {base_sha[:16]}…)。\n",
            "#   元との差: name・Md_moc_offset・moc_axis_limit / moc_corrector の明示"
            + ("・pw_upstream poly (pw_ramp を外す)" if sp["pw_upstream"] == "poly" else "")
            + (f"・spec.r_throat {rt_tok} (--r-throat、3 条件で共通)。k_f・格子・設定は元のまま。\n" if rt is not None
               else "。r_t・k_f・格子・設定は元のまま。\n"),
            "# ---- 以下は元の YAML のコメント ----\n"]
    return "".join(head + new)


def check_generated(text: str, base_doc: dict, cond: str, kind: str, value: float, rt: float | None = None) -> dict:
    """生成物の検査。戻り値 = 読み込み結果。破れれば ValueError (全項目をまとめて)。rt (--r-throat) が None なら spec.r_throat は
    元の値のまま (差を許さない)。"""
    d = yaml_load(text)
    sp = SPEC[cond]
    why = []
    diffs = diff_paths(base_doc, d)
    allowed = ALLOWED | ({RT_PATH} if rt is not None else set())
    extra = sorted(p for p in diffs if p not in allowed)
    if extra:
        why.append(f"元との差に許していないパスがある: {[('.'.join(map(str, p))) for p in extra]}")
    g = d.get("geometry") or {}
    want = {"name": out_name(cond, kind)[:-5]}
    if d.get("name") != want["name"]:
        why.append(f"name が {d.get('name')!r} ({want['name']!r} であること)")
    md = g.get("Md_moc_offset")
    if type(md) is not float or md != value:
        why.append(f"Md_moc_offset が {md!r} ({type(md).__name__}; {value!r} の float であること)")
    for k in ("pw_upstream", "moc_axis_limit", "moc_corrector"):
        if g.get(k) != sp[k]:
            why.append(f"geometry.{k} が {g.get(k)!r} ({sp[k]!r} であること)")
    if sp["pw_ramp"] and g.get("pw_ramp") != BASE_PW_RAMP:
        why.append(f"geometry.pw_ramp が {g.get('pw_ramp')!r} ({BASE_PW_RAMP} であること)")
    if not sp["pw_ramp"] and "pw_ramp" in g:
        why.append("geometry.pw_ramp が残っている (poly とは併記不可)")
    got_rt = (d.get("spec") or {}).get("r_throat")
    want_rt = rt if rt is not None else base_doc["spec"]["r_throat"]
    if type(got_rt) is not float or got_rt != want_rt:
        why.append(f"spec.r_throat が {got_rt!r} ({type(got_rt).__name__}; {want_rt!r} の float であること"
                   f"{'' if rt is not None else ' = 元の値'})")
    if why:
        raise ValueError(f"{out_name(cond, kind)}: " + "; ".join(why))
    return d


def effective(path: Path) -> dict:
    """runner_axismach の読み取りでの実効値 (import するだけ)。"""
    if str(DESIGN) not in sys.path:
        sys.path.insert(0, str(DESIGN))
    from forge_design.evaluate import runner_axismach as RA
    p = RA.load_problem(path)
    pwu = RA._pw_upstream(p.geometry)
    moc = RA._moc_keys(p.geometry)
    sp = p.spec
    return {"pw_upstream": pwu["value"], "pw_upstream_source": pwu["source"], "moc_axis_limit": moc[0], "moc_corrector": moc[1],
            "Md_moc_offset": float(p.geometry["Md_moc_offset"]), "r_throat": float(sp["r_throat"]),
            "cf_scale": float((p.raw.get("deltastar_initializer") or {}).get("cf_scale", float("nan"))),
            "mesh": copy.deepcopy(p.mesh)}


def build(value: float, out_dir: Path, rt: float | None = None) -> dict:
    """6 本のテキストを作って検査する (書かない)。戻り値 {ファイル名: テキスト}, 記録。rt = --r-throat [m] (None = 元の r_t)。"""
    tok = offset_token(value)
    rt_tok = float_token(rt, "r_throat") if rt is not None else None
    base_text = {k: (HERE / BASE[k]).read_text() for k in BASE}
    base_sha = {k: sha256_file(HERE / BASE[k]) for k in BASE}
    base_doc = {k: check_base(base_text[k], k) for k in BASE}
    base_pair = {p: v for p, v in diff_paths(base_doc["dry"], base_doc["cond"]).items() if p != ("name",)}
    texts, docs = {}, {}
    for c in CONDS:
        for k in BASE:
            t = render(base_text[k], BASE[k], base_sha[k], c, k, value, tok, rt, rt_tok)
            docs[(c, k)] = check_generated(t, base_doc[k], c, k, value, rt)
            texts[out_name(c, k)] = t
        pair = {p: v for p, v in diff_paths(docs[(c, "dry")], docs[(c, "cond")]).items() if p != ("name",)}
        if pair != base_pair:
            raise ValueError(f"{c}: dry と _cond の差が元の dry と _cond の差と違う ({sorted(pair)} / 元 {sorted(base_pair)})")
    # 3 条件の間の差は条件のキーだけ
    for a, b in (("N0", "N1"), ("N1", "N2")):
        for k in BASE:
            dd = {p for p in diff_paths(docs[(a, k)], docs[(b, k)]) if p != ("name",)}
            want = ({("geometry", "pw_upstream"), ("geometry", "pw_ramp")} if (a, b) == ("N0", "N1")
                    else {("geometry", "moc_axis_limit"), ("geometry", "moc_corrector")})
            if dd != want:
                raise ValueError(f"{a} と {b} ({k}) の差が {sorted(dd)} ({sorted(want)} だけであること)")
    rec = {"plan": [PLAN_U4, PLAN_V5], "tool": "make_ns_n012_problems.py", "tool_sha256": sha256_file(Path(__file__)),
           "md_offset": value, "md_offset_token": tok,
           "base": {k: {"file": BASE[k], "sha256": base_sha[k]} for k in BASE},
           "base_md_offset": base_doc["dry"]["geometry"]["Md_moc_offset"],
           "r_throat": rt, "r_throat_token": rt_tok, "base_r_throat": base_doc["dry"]["spec"]["r_throat"],
           "r_throat_source": (None if rt is None else "--r-throat (N0 の問題で CFD 前の solve_rt、許容差 1e-5 m・k_f 不変; plan §6 U4「r_t の解き直し」)"),
           "conditions": {c: {"role": SPEC[c]["role"], "files": {k: out_name(c, k) for k in BASE}} for c in CONDS}}
    return texts, rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--md-offset", required=True, help="出口較正の Md_moc_offset (Euler の E4 の結果)")
    ap.add_argument("--r-throat", default=None, help="spec.r_throat [m] (任意; 3 条件・6 本に同じ値。省略時は元の r_t)")
    ap.add_argument("--overwrite", action="store_true", help="中身の違う既存の生成物を上書きする")
    ap.add_argument("--out-dir", default=str(HERE))
    ap.add_argument("--check", action="store_true", help="書かずに既存の生成物と記録を照合する")
    a = ap.parse_args(argv)
    out_dir = Path(a.out_dir).resolve()
    try:
        value = parse_offset(a.md_offset)
        rt = parse_r_throat(a.r_throat) if a.r_throat is not None else None
        texts, rec = build(value, out_dir, rt)
    except ValueError as e:
        print(f"[make_ns_n012_problems] 止める (何も書いていない): {e}")
        return 2
    rec_path = out_dir / RECORD
    if a.check:
        bad = []
        for n, t in texts.items():
            p = out_dir / n
            if not p.is_file() or p.read_text() != t:
                bad.append(f"{n} が無い・今の --md-offset・--r-throat と元の YAML から作る中身と違う")
        try:
            old = json.loads(rec_path.read_text())
            if old.get("md_offset") != value or old.get("base") != rec["base"]:
                bad.append(f"{RECORD} の較正値・元の sha256 が違う ({old.get('md_offset')!r})")
            if old.get("r_throat") != rt:
                bad.append(f"{RECORD} の r_throat {old.get('r_throat')!r} が --r-throat {rt!r} と違う (None = 元の r_t)")
            for n in texts:
                if old.get("files", {}).get(n) != hashlib.sha256(texts[n].encode()).hexdigest():
                    bad.append(f"{RECORD} の {n} の sha256 が違う")
        except (OSError, ValueError) as e:
            bad.append(f"{RECORD} を読めない: {e}")
        print("CHECK: " + ("OK" if not bad else "NG\n  " + "\n  ".join(bad)))
        return 0 if not bad else 2
    clash = [n for n, t in texts.items() if (out_dir / n).is_file() and (out_dir / n).read_text() != t]
    if clash and not a.overwrite:
        print(f"[make_ns_n012_problems] 止める: 中身の違う既存のファイル {clash} (較正値を変えたなら --overwrite。"
              "準備済みの run があるときは上書きしない)")
        return 2
    for n, t in texts.items():
        (out_dir / n).write_text(t)
    eff = {}
    for c in CONDS:
        for k in BASE:
            e = effective(out_dir / out_name(c, k))
            sp = SPEC[c]
            bad = [f"{x} = {e[x]!r} ({sp[x]!r} であること)" for x in ("pw_upstream", "moc_axis_limit", "moc_corrector") if e[x] != sp[x]]
            if e["Md_moc_offset"] != value:
                bad.append(f"Md_moc_offset = {e['Md_moc_offset']!r}")
            want_rt = rt if rt is not None else rec["base_r_throat"]
            if e["r_throat"] != want_rt:
                bad.append(f"r_throat = {e['r_throat']!r} ({want_rt!r} であること)")
            if bad:
                for n in texts:
                    (out_dir / n).unlink(missing_ok=True)
                print(f"[make_ns_n012_problems] 止める (書いた 6 本を消した): {out_name(c, k)} の実効値が条件と違う: " + "; ".join(bad))
                return 2
            eff[out_name(c, k)] = {x: e[x] for x in ("pw_upstream", "pw_upstream_source", "moc_axis_limit", "moc_corrector",
                                                      "Md_moc_offset", "r_throat", "cf_scale")}
    # 3 条件で r_t・k_f・格子が同じ (元のまま)
    keys = {n: (v["r_throat"], v["cf_scale"]) for n, v in eff.items()}
    if len(set(keys.values())) != 1:
        print(f"[make_ns_n012_problems] 止める: r_t・k_f が条件で違う {keys}")
        return 2
    rec.update(generated=datetime.now(timezone.utc).isoformat(timespec="seconds"),
               files={n: hashlib.sha256(t.encode()).hexdigest() for n, t in texts.items()}, effective=eff)
    rec_path.parent.mkdir(parents=True, exist_ok=True)
    rec_path.write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    print(f"[make_ns_n012_problems] Md_moc_offset {rec['md_offset_token']} (元 {rec['base_md_offset']!r})、"
          f"r_throat {rec['r_throat_token'] or '元のまま'} (元 {rec['base_r_throat']!r}) で 6 本を書いた: "
          + ", ".join(texts) + f" → 記録 {RECORD}")
    for n, v in eff.items():
        print(f"  {n}: pw_upstream {v['pw_upstream']} ({v['pw_upstream_source']}), MOC {v['moc_axis_limit']}+{v['moc_corrector']}, "
              f"r_t {v['r_throat']}, k_f {v['cf_scale']}")
    if value == rec["base_md_offset"]:
        print("  注意: 較正値が元の YAML と同じ (N0 は今の生産と同じ問題になる)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
