#!/usr/bin/env python3
"""格子の節点数の上限 (plan tooling-nozzle-mesh-size-guard §5.1 #3、§6 V2・V3, 2026-10-14) の試験。

(a) `mesh2d.node_count` が実際に `generate_axisym_mesh` が作る座標の数と一致する (4 経路: 一様の nj・wall_first_frac_throat・
    axis_cap_frac・axis_gap_frac)。`_radial_fracs_capped` (内訳を `_capped_layout` に切り出した) が変更前の写しとビット一致する。
(b) 極小の axis_gap_frac (1e-7) と、等比の部分が長い指定 (nj 10 万) でも limit を渡せば 1 秒以内に limit + 1 以上を返す。
(c) 許可の値と理由の組み合わせ 4 通り (無し / 値と理由 / 値だけ / 理由だけ) と型違い。mesh_euler に 2 つのキーが書ける。
(d) V3 の実効の節点数の A/B (この試験だけ既定の上限を 25 万にする): Euler の ni 2000・nj 97・wall_first_frac 0.005 で
    A: axis_gap_frac 未指定 → 194,000 で prepare の判定を通る (設計チェーンに入る)、B: axis_gap_frac 0.005 → 40 万以上と数えて
    格子の生成の前・run dir を作る前に止まる (生成の関数と設計チェーンは呼んだら失敗する差し替え)。
(e) prepare (Euler) と prepare_ns (NS) が上限を超える問題で run dir を作らずに止まる・max_nodes と理由を書けば判定を通る・理由が
    無ければ止まる・上限以下は判定を通る (case/70 の problem_m5_air_c2b_euler.yaml: 1800 × 97 = 174,600)。事前の数えと実際の
    節点数が違えば格子を作った後に止まる。prepare_info.json の mesh 欄 (`_mesh_record`) に節点数・上限・理由・出所。
(f) V2 要求カードの生成器 (case/70 の card_m5_air_demo_c2b.yaml、r_t は --r-throat で固定して速くする):
    上限以下のカードは、判定を外した生成と節点数のコメント行を除いて同じ生成物で、YAML のコメント・標準出力・card_log に
    NS・Euler の節点数と上限が出る。NS の ni が 6000 になるカード (第 1 層を薄く) は `_mesh_ar_full` (全格子) を呼ぶ前に止まる
    (呼んだら失敗する差し替え)。mesh.max_nodes・mesh_euler.max_nodes と理由を書けば通る (NS だけ許すと Euler で止まる、理由が無ければ
    止まる)。明示の mesh_euler.ni は保持される (上限以下・超過; 超過は許可を書けば通る)。
    **変更前の生成器との比較は 1 回だけ実施** (2026-10-14、plan §9): c2b のカードを r_t の固定点込みで生成し、YAML・Euler の写し・
    card_log が節点数のコメント行・時間の値・新しいキー nodes を除いて一致した。

実行 (forge-deliv の直下で): /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_mesh_size_guard_tests.py
所要: 約 2 分 (生成器 6 回 × 約 13〜20 s、設計チェーン 1 回)。
"""
import contextlib
import copy
import io
import json
import re
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import yaml

DESIGN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(DESIGN))

from forge_design.evaluate import runner_axismach as RA  # noqa: E402
from forge_design.intake import wt_card as W  # noqa: E402
from forge_design.meshing import mesh2d as M  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

REPO = DESIGN.parent
C70 = REPO / "case" / "70.m5_air_wt_demo"
CARD_C2B = C70 / "card_m5_air_demo_c2b.yaml"
PROB_EU = C70 / "problem_m5_air_c2b_euler.yaml"
PROB_NS = C70 / "problem_m5_air_c2b.yaml"
RT_C2B = 0.024973576891821576          # c2b のカードの r_t の固定点 (2026-10-14 の生成; --r-throat で与えて解かない)
FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail != "" else ""), flush=True)
    if not cond:
        FAIL += 1


def raises(fn, exc=ValueError):
    """fn() が exc を出せばそのメッセージ、出さなければ None (他の例外はそのまま上げる)。"""
    try:
        fn()
    except exc as ex:
        return str(ex)
    return None


class Wall:
    """解析的な仮の壁 (r_t 単位): 入口 x_in = −12、出口 x_e = 95、スロート x = 0 で r = 1 (run_mesh_euler_tests と同じ)。"""
    x_in, x_e = -12.0, 95.0

    def r(self, x):
        x = np.asarray(x, dtype=float)
        return np.where(x < 0, 1.0 + 0.04 * x * x, 1.0 + 0.09 * x - 0.0003 * x * x)


class Reached(Exception):
    """差し替えた関数が呼ばれた (判定を通った/通ってはいけない所に来た)。"""


@contextlib.contextmanager
def patched(obj, name, value):
    old = getattr(obj, name)
    setattr(obj, name, value)
    try:
        yield
    finally:
        setattr(obj, name, old)


def stub(label):
    def f(*a, **k):
        raise Reached(label)
    return f


TMP = Path(tempfile.mkdtemp(prefix="mesh_size_guard_tests_"))
t_all = time.time()

# ---------------------------------------------------------------------------------------------- (a) node_count と実際の座標数
print("== (a) node_count と実際に作った座標の数 (4 経路)")
WL = Wall()
paths = [
    ("一様の nj", dict(ni=300, nj=97, wall_first_frac=0.005)),
    ("wall_first_frac_throat", dict(ni=300, nj=97, wall_first_frac=1.15e-5, wall_first_frac_throat=2.75e-6, wall_first_up_x0=-9.0,
                                    wall_first_up_x1=-4.0, wall_first_blend_x0=1.0, wall_first_blend_x1=17.0)),
    ("axis_cap_frac", dict(ni=300, nj=97, wall_first_frac=1.3e-5, wall_first_frac_throat=4.5e-6, axis_cap_frac=0.02)),
    ("axis_gap_frac (nj は導出値)", dict(ni=300, nj=97, wall_first_frac=0.005, axis_gap_frac=0.005)),
    ("axis_gap_frac (壁に寄せた第 1 層)", dict(ni=300, nj=97, wall_first_frac=4.5e-5, axis_gap_frac=0.025)),
]
table = []
for name, kw in paths:
    prm = M.Mesh2DParams(**kw)
    c, q, b = M.generate_axisym_mesh(WL, prm)
    n = M.node_count(prm)
    nj_eff, exact = M.radial_node_count(prm)
    table.append((name, kw["ni"], kw["nj"], nj_eff, n, c.shape[0]))
    check(f"node_count = 実際の座標数: {name}", n == c.shape[0] and exact and n == kw["ni"] * nj_eff, f"{n} / {c.shape[0]} (nj {nj_eff})")
check("axis_gap_frac 0.005・第 1 層 0.005 の導出 nj は 201 (指定 97 の 2 倍超)", table[3][3] == 201, table[3])
for cmb, why in ((dict(ni=20, nj=97, axis_cap_frac=0.02, axis_gap_frac=0.01), "axis_cap_frac と axis_gap_frac"),
                 (dict(ni=20, nj=97, axis_gap_frac=0.01, wall_first_frac_throat=1e-5), "axis_gap_frac と wall_first_frac_throat")):
    e1 = raises(lambda cmb=cmb: M.node_count(M.Mesh2DParams(**cmb)))
    e2 = raises(lambda cmb=cmb: M.generate_axisym_mesh(WL, M.Mesh2DParams(**cmb)))
    check(f"併用できない指定 ({why}) は node_count も生成も同じ ValueError", e1 is not None and e1 == e2, e1)


def _old_radial_fracs_capped(nj_base, first_frac, cap):
    """変更前 (commit 88675245) の mesh2d._radial_fracs_capped の写し。"""
    s0 = M._radial_fracs(nj_base, first_frac)
    g0 = np.diff(s0)[::-1]
    q = g0[1] / g0[0]
    gaps = []
    g = g0[0]
    while g < cap and sum(gaps) + g < 1.0:
        gaps.append(g)
        g *= q
    rest = 1.0 - sum(gaps)
    m = int(np.ceil(rest / cap - 1e-12))
    gaps = gaps + [rest / m] * m
    s = np.concatenate([[0.0], np.cumsum(gaps[::-1])])
    return s / s[-1]


with np.errstate(divide="ignore"):       # 等間隔で残りが丸めの 0 になる組 (m = 0) は新旧とも rest/0 を作って捨てる (同じ振る舞い)
    bad = [(nj, fr, cap) for nj in (10, 65, 97) for fr in (0.2, 5e-3, 4.5e-5, 1e-6) for cap in (0.5, 0.05, 0.025, 0.005, 1e-3)
           if not np.array_equal(M._radial_fracs_capped(nj, fr, cap), _old_radial_fracs_capped(nj, fr, cap))]
check("_radial_fracs_capped が変更前の写しとビット一致 (60 通り)", not bad, bad[:3])

# ---------------------------------------------------------------------------------------------- (b) 打ち切り
print("== (b) 極小の axis_gap_frac・長い等比の部分で limit を渡すと速く返る")
t = time.perf_counter()
n_tiny = M.node_count(M.Mesh2DParams(ni=2000, nj=97, wall_first_frac=0.005, axis_gap_frac=1e-7), limit=300_000)
dt = time.perf_counter() - t
check("axis_gap_frac 1e-7: limit 30 万で 1 秒以内に limit + 1 以上", n_tiny > 300_000 and dt < 1.0, f"{n_tiny}、{dt * 1e3:.2f} ms")
t = time.perf_counter()
prm_long = M.Mesh2DParams(ni=2000, nj=100_000, wall_first_frac=2e-5, axis_gap_frac=0.5)   # first·n ≥ 1 で等間隔 → 等比の部分 ≈ 10 万
n_long = M.node_count(prm_long, limit=300_000)
nj_long, exact_long = M.radial_node_count(prm_long, max_nj=300_000 // 2000)
dt = time.perf_counter() - t
check("nj 10 万・axis_gap_frac 0.5 (等比の部分が長い): 等比の数え上げを打ち切って 1 秒以内に limit + 1 以上",
      n_long > 300_000 and not exact_long and nj_long > 150 and dt < 1.0, f"{n_long}、nj ≥ {nj_long}、{dt * 1e3:.1f} ms")
# result 段レビュー M1 (2026-10-14): 指定 nj が極端に大きいと `_capped_layout` が数える前に nj 個の配列を作っていた → 指定 nj − 1 の下限で配列の前に打ち切る
with patched(M, "_radial_fracs", lambda *a, **k: (_ for _ in ()).throw(AssertionError("_radial_fracs を呼んだ (配列を確保)"))):
    t = time.perf_counter()
    n_huge = M.node_count(M.Mesh2DParams(ni=2000, nj=1_000_000_000, wall_first_frac=0.005, axis_gap_frac=0.005), limit=300_000)
    dt = time.perf_counter() - t
check("指定 nj 10 億・axis_gap_frac 0.005: `_radial_fracs` (配列) を呼ばずに limit + 1 以上", n_huge > 300_000 and dt < 1.0, f"{n_huge}、{dt * 1e3:.2f} ms")
viol = []
for nj in (3, 10, 65, 97, 201, 1000):
    for fr in (0.4, 0.05, 5e-3, 4.5e-5, 1e-6):
        for cap in (0.9, 0.5, 0.05, 0.005, 1e-3):
            g0 = np.diff(M._radial_fracs(nj, fr))[::-1]
            if not (np.all(np.isfinite(g0)) and g0[0] > 0 and np.isfinite(g0[1] / g0[0])):
                continue    # 等比の配点が壊れる組 (nj が大きく第 1 層が薄いと Newton が退化し第 1 層 0; 既存の欠陥、plan §5.1 #8) は下限の前提の外
            try:
                nj_eff, ok = M.radial_node_count(M.Mesh2DParams(ni=10, nj=nj, wall_first_frac=fr, axis_gap_frac=cap))
            except Exception:  # noqa: BLE001 — 実現不能な組は数えない
                continue
            if nj_eff < nj - 1:
                viol.append((nj, fr, cap, nj_eff))
check("打ち切りの下限 (実効の nj ≥ 指定 nj − 1) が等比の配点の正常な組で成り立つ", not viol, viol[:5])
t = time.perf_counter()
e = raises(lambda: RA.check_node_budget(M.Mesh2DParams(ni=2000, nj=97, wall_first_frac=0.005, axis_gap_frac=1e-7), {}, "mesh_euler"))
dt = time.perf_counter() - t
check("check_node_budget も極小の axis_gap_frac で速く止まる", e is not None and "300000" in e and dt < 1.0, f"{dt * 1e3:.2f} ms")

# ---------------------------------------------------------------------------------------------- (c) 許可の組み合わせ
print("== (c) 許可の値と理由の組み合わせ")
b0 = RA.node_budget({}, "mesh")
check("無し → 既定 30 万 (出所 default)", b0 == {"max_nodes": 300_000, "max_nodes_reason": None, "max_nodes_source": "default"}
      and RA.MAX_NODES_2D == 300_000, b0)
b1 = RA.node_budget({"max_nodes": 600_000, "max_nodes_reason": "格子収束の細格子 (ユーザ了承 2026-10-14)"}, "mesh")
check("値と理由 → その値 (出所 explicit)", b1["max_nodes"] == 600_000 and b1["max_nodes_source"] == "explicit"
      and b1["max_nodes_reason"].startswith("格子収束"), b1)
e = raises(lambda: RA.node_budget({"max_nodes": 600_000}, "mesh"))
check("値だけ → 止まる (理由が無い)", e is not None and "mesh.max_nodes_reason" in e, e)
e = raises(lambda: RA.node_budget({"max_nodes_reason": "細格子"}, "mesh_euler"))
check("理由だけ → 止まる (値が無い)", e is not None and "mesh_euler.max_nodes" in e, e)
for v, why in (("600000", "文字列"), (True, "真偽値"), (0, "0"), (-5, "負"), (6.0e5, "浮動小数")):
    e = raises(lambda v=v: RA.node_budget({"max_nodes": v, "max_nodes_reason": "x"}, "mesh"))
    check(f"max_nodes が{why}なら止まる", e is not None and "正の整数" in e, e)
for r, why in (("", "空"), ("   ", "空白だけ"), (123, "数値")):
    e = raises(lambda r=r: RA.node_budget({"max_nodes": 600_000, "max_nodes_reason": r}, "mesh"))
    check(f"理由が{why}なら止まる", e is not None and "空でない文字列" in e, e)
check("mesh_euler に max_nodes・max_nodes_reason が書ける (MESH_EULER_KEYS)",
      {"max_nodes", "max_nodes_reason"} <= set(RA.MESH_EULER_KEYS))
big = M.Mesh2DParams(ni=4000, nj=97, wall_first_frac=0.005)
e = raises(lambda: RA.check_node_budget(big, {}, "mesh"))
check("上限超過のメッセージ: 実効の節点数・ni・nj・上限・許可の書き方", e is not None and all(s in e for s in (
    "388000", "ni 4000", "nj 97", "上限 300000", "max_nodes: 388000", "max_nodes_reason", "ユーザの了承")), e)
r_ok = RA.check_node_budget(big, {"max_nodes": 400_000, "max_nodes_reason": "試験"}, "mesh")
check("許可を書けば通り、記録に節点数と上限・理由・出所", r_ok == {"n_nodes": 388_000, "ni": 4000, "nj": 97, "max_nodes": 400_000,
                                                       "max_nodes_reason": "試験", "max_nodes_source": "explicit"}, r_ok)


# ---------------------------------------------------------------------------------------------- 問題 YAML の写し
def write_problem(src: Path, name: str, edit) -> Path:
    doc = yaml.safe_load(src.read_text())
    edit(doc)
    out = TMP / name
    out.write_text(yaml.safe_dump(doc, allow_unicode=True, sort_keys=False))
    return out


def set_block(blk, **kw):
    def edit(d):
        for k, v in kw.items():
            if v is None:
                d[blk].pop(k, None)
            else:
                d[blk][k] = v
    return edit


# ---------------------------------------------------------------------------------------------- (d) V3 実効の節点数の A/B
print("== (d) V3: 実効の節点数の A/B (Euler の ni 2000・nj 97・wall_first_frac 0.005、既定の上限を 25 万にする)")
base_eu = {"ni": 2000, "nj": 97, "wall_first_frac": 0.005, "throat_refine": 4.0, "throat_width": 3.0}
yA = write_problem(PROB_EU, "v3_A.yaml", lambda d: d.__setitem__("mesh_euler", dict(base_eu)))
yB = write_problem(PROB_EU, "v3_B.yaml", lambda d: d.__setitem__("mesh_euler", dict(base_eu, axis_gap_frac=0.005)))
pA, pB = load_problem(yA), load_problem(yB)
nA, nB = M.node_count(RA.mesh_params_euler(pA, 1.0)), M.node_count(RA.mesh_params_euler(pB, 1.0))
check("A の実効の節点数 194,000・B は 40 万以上 (axis_gap_frac で nj が導出値になる)", nA == 194_000 and nB >= 400_000, f"A {nA}、B {nB}")
with patched(RA, "MAX_NODES_2D", 250_000), patched(RA, "design_chain", stub("design_chain")), \
        patched(RA, "generate_axisym_mesh", stub("generate_axisym_mesh")):
    rdA, rdB = TMP / "run_v3_A", TMP / "run_v3_B"
    try:
        RA.prepare(yA, rdA)
        rA = "例外なし"
    except Reached as ex:
        rA = f"Reached {ex}"
    check("A: 194,000 ≤ 25 万で判定を通り設計チェーンに入る", rA == "Reached design_chain" and not rdA.exists(), rA)
    try:
        RA.prepare(yB, rdB)
        rB = "例外なし"
    except Reached as ex:
        rB = f"Reached {ex} (判定より先に進んだ)"
    except ValueError as ex:
        rB = str(ex)
    check("B: 40 万以上と数えて、生成の関数・設計チェーンを呼ぶ前・run dir を作る前に止まる",
          rB.startswith("mesh_euler の格子の実効の節点数") and str(nB) in rB and "nj 201" in rB and "上限 250000" in rB
          and "axis_gap_frac からの導出値" in rB and not rdB.exists(), rB[:200])

# ---------------------------------------------------------------------------------------------- (e) prepare / prepare_ns
print("== (e) prepare (Euler) と prepare_ns (NS) の判定")


def run_until_chain(fn, problem, rd):
    """設計チェーンを呼んだら止まる差し替えで prepare* を呼ぶ。'chain' (判定を通った) か ValueError のメッセージ。"""
    with patched(RA, "design_chain", stub("design_chain")):
        try:
            fn(problem, rd)
        except Reached:
            return "chain"
        except ValueError as ex:
            return str(ex)
    return "例外なし"


r = run_until_chain(RA.prepare, PROB_EU, TMP / "run_eu_c70")
check("prepare: case/70 の problem_m5_air_c2b_euler.yaml (1800 × 97 = 174,600) は判定を通る", r == "chain", r[:120])
r = run_until_chain(RA.prepare_ns, PROB_NS, TMP / "run_ns_c70")
check("prepare_ns: case/70 の problem_m5_air_c2b.yaml (1800 × 97) は判定を通る", r == "chain", r[:120])
y = write_problem(PROB_EU, "eu_over.yaml", set_block("mesh_euler", ni=4000))
rd = TMP / "run_eu_over"
r = run_until_chain(RA.prepare, y, rd)
check("prepare: mesh_euler.ni 4000 (388,000) は run dir を作らずに止まる", "388000" in r and "mesh_euler" in r and not rd.exists(), r[:120])
y = write_problem(PROB_EU, "eu_allow.yaml", set_block("mesh_euler", ni=4000, max_nodes=400_000, max_nodes_reason="試験"))
r = run_until_chain(RA.prepare, y, TMP / "run_eu_allow")
check("prepare: mesh_euler.max_nodes 400000 と理由を書けば判定を通る", r == "chain", r[:120])
y = write_problem(PROB_EU, "eu_noreason.yaml", set_block("mesh_euler", ni=4000, max_nodes=400_000))
rd = TMP / "run_eu_noreason"
r = run_until_chain(RA.prepare, y, rd)
check("prepare: 理由が無ければ止まる", "max_nodes_reason" in r and not rd.exists(), r[:120])
y = write_problem(PROB_NS, "ns_over.yaml", set_block("mesh", ni=4000))
rd = TMP / "run_ns_over"
r = run_until_chain(RA.prepare_ns, y, rd)
check("prepare_ns: mesh.ni 4000 (388,000) は run dir を作らずに止まる", r.startswith("mesh の格子") and "388000" in r and not rd.exists(),
      r[:120])
y = write_problem(PROB_NS, "ns_allow.yaml", set_block("mesh", ni=4000, max_nodes=400_000, max_nodes_reason="試験"))
r = run_until_chain(RA.prepare_ns, y, TMP / "run_ns_allow")
check("prepare_ns: mesh.max_nodes と理由を書けば判定を通る", r == "chain", r[:120])
y = write_problem(PROB_NS, "ns_noreason.yaml", set_block("mesh", ni=4000, max_nodes_reason="理由だけ"))
rd = TMP / "run_ns_noreason"
r = run_until_chain(RA.prepare_ns, y, rd)
check("prepare_ns: 理由だけ (値が無い) は止まる", "mesh.max_nodes" in r and not rd.exists(), r[:120])
# 事前の数えと実際の節点数の照合 (数えを 1 つずらすと、格子を作った後・msh を書く前に止まる)
real_count = RA.node_count
with patched(RA, "node_count", lambda prm, limit=None: real_count(prm, limit) - 1), \
        patched(RA, "write_msh41_2d", stub("write_msh41_2d")):
    try:
        RA.prepare(PROB_EU, TMP / "run_eu_mismatch")
        r = "例外なし"
    except RuntimeError as ex:
        r = str(ex)
    except Reached as ex:
        r = f"Reached {ex}"
check("事前の数えと実際の節点数が違えば格子を作った後・msh を書く前に止まる", "一致しない" in r and "174600" in r, r[:120])
# prepare_info.json の mesh 欄 (ハッシュは変換後の h5 が要るので差し替え)
pe = load_problem(PROB_EU)
mp = RA.mesh_params_euler(pe, 1.0)
cq = M.generate_axisym_mesh(WL, mp)
nodes = RA.check_node_budget(mp, RA.mesh_euler_block(pe), "mesh_euler")
with patched(RA, "_mesh_hashes", lambda *a: {"stub": True}):
    rec = RA._mesh_record(pe, mp, "mesh_euler", cq[0].shape[0] // mp.ni, TMP, *cq, nodes)
check("prepare_info.json の mesh 欄に n_nodes・max_nodes・max_nodes_reason・max_nodes_source",
      (rec["n_nodes"], rec["max_nodes"], rec["max_nodes_reason"], rec["max_nodes_source"], rec["ni"], rec["nj"])
      == (174_600, 300_000, None, "default", 1800, 97), {k: rec[k] for k in ("n_nodes", "max_nodes", "max_nodes_reason", "max_nodes_source")})

# ---------------------------------------------------------------------------------------------- (f) V2 要求カードの生成器
print("== (f) V2: 要求カードの生成器 (c2b、r_t 固定)")
CARD = yaml.safe_load(CARD_C2B.read_text())


def card_with(name, **ov):
    c = copy.deepcopy(CARD)
    c["overrides"] = dict(c.get("overrides") or {}, **ov)
    p = TMP / f"{name}.yaml"
    p.write_text(yaml.safe_dump(c, allow_unicode=True, sort_keys=False))
    return p


def norm_yaml(text: str) -> str:
    """比べる前の正規化: 生成日 (1 行目・_card.generated)、時間の値 (time_s)、節点数のコメント行、出力先のパスを除く。"""
    out = []
    for ln in text.splitlines():
        if ln.startswith("# 節点数"):
            continue
        ln = re.sub(r"\(\d{4}-\d{2}-\d{2}\)", "(DATE)", ln)
        ln = re.sub(r"generated: '\d{4}-\d{2}-\d{2}'", "generated: DATE", ln)
        ln = re.sub(r"time_s: [0-9.e+-]+", "time_s: T", ln)
        out.append(ln)
    return "\n".join(out)


NODE_LINE = ("節点数 (格子を作る前の判定): NS (mesh) ni 1800 × nj 97 = 174600 (上限 300000、既定)・"
             "Euler (mesh_euler) ni 1800 × nj 97 = 174600 (上限 300000、既定)")
# 上限以下 (CLI 経由: 標準出力も見る)
outA = TMP / "gen_on" / "problem_c2b.yaml"
outA.parent.mkdir()
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    rc = W.main([str(CARD_C2B), "--out", str(outA), "--r-throat", repr(RT_C2B)])
stdout = buf.getvalue()
textA = outA.read_text()
logA = json.loads(Path(str(outA) + ".card_log.json").read_text())
rawA = yaml.safe_load(textA)
check("上限以下: 生成できる (NS・Euler とも 1800 × 97)", rc == 0 and rawA["mesh"]["ni"] == 1800 and rawA["mesh_euler"]["ni"] == 1800, rc)
check("YAML のコメントに NS・Euler の節点数と上限", ("# " + NODE_LINE) in textA, [ln for ln in textA.splitlines() if "節点数" in ln][:1])
check("標準出力に NS・Euler の節点数と上限", NODE_LINE in stdout, [ln for ln in stdout.splitlines() if "節点数" in ln][:1])
check("card_log の nodes (NS・Euler) と mesh_estimate.nodes", logA["nodes"]["mesh"]["n_nodes"] == 174_600
      and logA["nodes"]["mesh_euler"]["n_nodes"] == 174_600 and logA["nodes"]["mesh"]["max_nodes_source"] == "default"
      and logA["mesh_estimate"]["nodes"] == logA["nodes"]["mesh"], logA["nodes"])
check("YAML のデータ (_card を含む) に節点数のキーを足していない", "nodes" not in rawA["_card"].get("mesh_estimate", {})
      and "max_nodes" not in rawA["mesh"] and "max_nodes" not in rawA["mesh_euler"])
# 判定を外した生成 (check_node_budget を記録だけ返す差し替え) と、節点数のコメント行を除いて同じ
outN = TMP / "gen_off" / "problem_c2b.yaml"


def no_guard(mp, m, block, allow_hint=None):
    return {"n_nodes": -1, "ni": int(mp.ni), "nj": int(mp.nj), "max_nodes": -1, "max_nodes_reason": None, "max_nodes_source": "default"}


with patched(RA, "check_node_budget", no_guard):
    W.build_problem(CARD_C2B, outN, r_throat_m=RT_C2B, log=lambda *a, **k: None)
same_p = norm_yaml(textA) == norm_yaml(outN.read_text())
eA, eN = outA.with_name("problem_c2b_euler.yaml"), outN.with_name("problem_c2b_euler.yaml")
same_e = eA.is_file() and eN.is_file() and norm_yaml(eA.read_text()) == norm_yaml(eN.read_text())
check("上限以下: 判定を外した生成と同じ生成物 (問題 YAML・Euler の写し; 生成日・時間の値・節点数のコメント行を除く)", same_p and same_e,
      f"問題 {same_p}、Euler の写し {same_e}")

# NS の ni が 6000 になるカード (第 1 層を薄く)
thin = {"mesh.wall_first_frac": 2.0e-6, "mesh.wall_first_frac_throat": 5.0e-7}
c6k = card_with("card_ns6000", **thin)
with patched(W, "_mesh_ar_full", stub("_mesh_ar_full")):
    try:
        W.build_problem(c6k, TMP / "p6k.yaml", r_throat_m=RT_C2B, log=lambda *a, **k: None)
        r = "例外なし"
    except Reached as ex:
        r = f"Reached {ex} (全格子を作った)"
    except ValueError as ex:
        r = str(ex)
check("NS の ni 6000 (第 1 層を薄く): 全格子 (_mesh_ar_full) を呼ぶ前に止まる", "estimate_mesh" in r and "582000" in r
      and "ni 6000" in r and "NI_CAP" in r and "overrides: {mesh.max_nodes: 582000" in r and not (TMP / "p6k.yaml").exists(), r[:240])
allow = {"mesh.max_nodes": 600_000, "mesh.max_nodes_reason": "試験 (第 1 層を薄くしたカード)",
         "mesh_euler.max_nodes": 600_000, "mesh_euler.max_nodes_reason": "試験 (NS の ni の写し)"}
out6 = TMP / "p6k_allow.yaml"
log6 = W.build_problem(card_with("card_ns6000_allow", **thin, **allow), out6, r_throat_m=RT_C2B, log=lambda *a, **k: None)
raw6 = yaml.safe_load(out6.read_text())
check("mesh.max_nodes・mesh_euler.max_nodes と理由を書けば通る (NS・Euler とも 6000 × 97)",
      raw6["mesh"]["ni"] == 6000 and raw6["mesh_euler"]["ni"] == 6000 and log6["nodes"]["mesh"]["n_nodes"] == 582_000
      and log6["nodes"]["mesh_euler"]["max_nodes_source"] == "explicit" and raw6["mesh_euler"]["max_nodes"] == 600_000
      and "明示: 試験 (第 1 層を薄くしたカード)" in out6.read_text(), log6["nodes"])
check("許可を書いた生成物は prepare の判定も通る (mesh_euler に max_nodes のキー)",
      run_until_chain(RA.prepare, out6.with_name("p6k_allow_euler.yaml"), TMP / "run_p6k_eu") == "chain")
ns_only = {k: v for k, v in allow.items() if k.startswith("mesh.")}
try:
    W.build_problem(card_with("card_ns6000_nsonly", **thin, **ns_only), TMP / "p6k_ns.yaml", r_throat_m=RT_C2B, log=lambda *a, **k: None)
    r = "例外なし"
except ValueError as ex:
    r = str(ex)
check("NS だけ許すと Euler (mesh.ni の写し) で止まる", r.startswith("要求カード → 問題 (Euler の格子 ni 6000 は mesh.ni 6000 の写し)")
      and "mesh_euler.max_nodes" in r and not (TMP / "p6k_ns.yaml").exists(), r[:200])
t = time.perf_counter()
r = raises(lambda: W.build_problem(card_with("card_noreason", **{"mesh.max_nodes": 600_000}), TMP / "p_nr.yaml", r_throat_m=RT_C2B,
                                   log=lambda *a, **k: None))
check("理由が無ければ止まる (重い処理の前: card_to_problem)", r is not None and "mesh.max_nodes_reason" in r
      and time.perf_counter() - t < 5.0, f"{(r or '')[:120]} ({time.perf_counter() - t:.2f} s)")
# 明示の mesh_euler.ni は保持される
outU = TMP / "p_eu1000.yaml"
W.build_problem(card_with("card_eu1000", **{"mesh_euler.ni": 1000}), outU, r_throat_m=RT_C2B, log=lambda *a, **k: None)
rawU = yaml.safe_load(outU.read_text())
check("明示の mesh_euler.ni 1000 (上限以下) は保持される (mesh.ni は 1800)", rawU["mesh_euler"]["ni"] == 1000 and rawU["mesh"]["ni"] == 1800,
      (rawU["mesh_euler"]["ni"], rawU["mesh"]["ni"]))
t = time.perf_counter()
r = raises(lambda: W.build_problem(card_with("card_eu4000", **{"mesh_euler.ni": 4000}), TMP / "p_eu4000.yaml", r_throat_m=RT_C2B,
                                   log=lambda *a, **k: None))
check("明示の mesh_euler.ni 4000 (388,000、超過) は重い処理の前に止まる", r is not None and "明示の mesh_euler.ni 4000" in r and "388000" in r
      and "mesh_euler.max_nodes" in r and time.perf_counter() - t < 5.0, f"{(r or '')[:160]} ({time.perf_counter() - t:.2f} s)")
outO = TMP / "p_eu4000_allow.yaml"
logO = W.build_problem(card_with("card_eu4000_allow", **{"mesh_euler.ni": 4000, "mesh_euler.max_nodes": 400_000,
                                                         "mesh_euler.max_nodes_reason": "試験"}), outO, r_throat_m=RT_C2B, log=lambda *a, **k: None)
rawO = yaml.safe_load(outO.read_text())
check("明示の mesh_euler.ni 4000 は許可を書けば通り、上書きされない", rawO["mesh_euler"]["ni"] == 4000 and rawO["mesh"]["ni"] == 1800
      and logO["nodes"]["mesh_euler"]["n_nodes"] == 388_000, logO["nodes"]["mesh_euler"])

print(f"\n所要 {time.time() - t_all:.0f} s")
print("ALL OK" if FAIL == 0 else f"{FAIL} FAIL")
sys.exit(1 if FAIL else 0)
