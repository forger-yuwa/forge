"""格子の段の契約 (plan tooling-design-orchestrator §4.6、§5.1 #8c; 2026-10-18 codex 諮問 mesh-limits)。

契約: 「問題 YAML の mesh 節 + 格子の不変部分のハッシュ + 品質の VERDICT + 依頼ごとの承認済みの節点数の上限と承認元」。
実装は 2 つ。手組み (case/71) は `case/71.m55_combustion_wt/_mesh/hand_mesh.py build --approval`・`check` (recipe v1 の格子の段)。
雛形 (`wt_card.estimate_mesh` が mesh 節を決める問題、例 case/70) はこの `check_template`:

1. 承認の照合 (`approval.check_card`: カードが承認の原本か承認の記録を指す試行の入力か、カードだけで決まる保護の条件) と、問題 YAML の照合
   (`approval.check`: 問題 YAML がカードから作られたか・保護の条件・節点数の上限 ≤ 承認の mesh_limits)、Euler の写しの上限。どれかが
   通らなければ何も作らずに `ApprovalError`。
2. 実効の節点数を承認の上限で数える (`runner_axismach.check_node_budget`; 格子の配列は作らない)。超えれば何も作らずに `ValueError`。
3. 実際の準備の経路で格子を作る (forge は起動しない): Euler は `runner_axismach.prepare`、NS は `prepare_ns` (積分法の初期壁
   `deltastar_loop.resolve_integral_initializer`、停止規則なし — 停止規則は Euler の参照の結果が要るので、この段では使わない)。
   生成の手順を写さない (codex 諮問 mesh-limits (5))。
4. 記録 (`MESH_STAGE.json`): 問題 YAML の sha256、節点数と承認の上限・承認元、準備の記録の格子 (`prepare_info.json` の mesh: 節点数・
   ハッシュ)、品質の VERDICT (`MESH_QUALITY.txt` の最後の行)。

部品の検証であり、campaign の通しの成功とは言わない (recipe v1 には入れない)。
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

FORMAT = "forge.mesh_stage/v1"
RECORD = "MESH_STAGE.json"


def _sha256(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _quality_line(run_dir: Path) -> str | None:
    p = Path(run_dir) / "MESH_QUALITY.txt"
    if not p.is_file():
        return None
    return next((ln for ln in reversed(p.read_text().splitlines()) if ln.startswith("VERDICT")), None)


def companion_of(problem) -> Path:
    """Euler の写し (`<stem>_euler.yaml`、`wt_card.write_euler_companion`) があればそれ、無ければ問題 YAML そのもの。"""
    p = Path(problem)
    eo = p.with_name(p.stem + "_euler" + p.suffix)
    return eo if eo.is_file() else p


def precheck(problem, *, approval, card, problem_euler=None) -> dict:
    """1・2 (何も作らない)。通らなければ ApprovalError / ValueError。戻り {limits, nodes_ns, nodes_euler, problem_euler}。"""
    from ..evaluate import runner_axismach as RA
    from ..intake import approval as AP
    from ..probdef import load_problem
    import yaml
    ck = AP.check_card(approval, card)
    if ck["status"] != "OK":
        raise AP.ApprovalError("カードが承認と結び付かない — 格子を作らない: " + "; ".join(ck["problems"])[:600])
    lim = AP.mesh_limits_for_card(approval, card)
    cp = AP.check(approval, card, problem)
    if cp["status"] != "OK":
        raise AP.ApprovalError("問題 YAML が承認と違う — 格子を作らない: " + "; ".join(cp["problems"] + cp["diffs"])[:600])
    pe = Path(problem_euler) if problem_euler is not None else companion_of(problem)
    le = AP.mesh_limit_problems({"max_nodes_ns": lim["max_nodes_ns"], "max_nodes_euler": lim["max_nodes_euler"]},
                                yaml.safe_load(pe.read_text()))
    if le:
        raise AP.ApprovalError(f"Euler の問題 YAML ({pe.name}) の節点数の上限が承認と合わない — 格子を作らない: " + "; ".join(le)[:600])
    src = f"承認の記録 {Path(lim['approval_file']).name} (sha256 {lim['approval_sha256'][:16]}…) の mesh_limits"
    p = load_problem(problem)
    pe_p = load_problem(pe)
    nodes_ns = RA.check_node_budget(RA.mesh_params(p, 1.0, ni=561, nj=97, wall_first_frac=4.5e-5),
                                    {**dict(p.mesh), "max_nodes": lim["max_nodes_ns"], "max_nodes_reason": src + ".max_nodes_ns"}, "mesh")
    nodes_eu = RA.check_node_budget(RA.mesh_params_euler(pe_p, 1.0),
                                    {**dict(RA.mesh_euler_block(pe_p)), "max_nodes": lim["max_nodes_euler"],
                                     "max_nodes_reason": src + ".max_nodes_euler"}, "mesh_euler")
    return {"limits": lim, "nodes_ns": nodes_ns, "nodes_euler": nodes_eu, "problem_euler": str(pe)}


def check_template(problem, out_dir, *, approval, card, problem_euler=None) -> dict:
    """雛形の格子の段 (docstring の 1〜4)。out_dir/euler・out_dir/ns に準備の run (forge なし) を作り、out_dir/MESH_STAGE.json を書く。
    品質 FAIL は準備の関数の例外 (RuntimeError) で止まる。戻り: 記録の dict。"""
    from ..evaluate import runner_axismach as RA
    from ..feedback.deltastar_loop import resolve_integral_initializer
    out_dir = Path(out_dir)
    pre = precheck(problem, approval=approval, card=card, problem_euler=problem_euler)
    out_dir.mkdir(parents=True, exist_ok=True)
    rec = {"format": FORMAT, "problem": str(Path(problem).resolve()), "problem_sha256": _sha256(problem),
           "problem_euler": str(Path(pre["problem_euler"]).resolve()), "problem_euler_sha256": _sha256(pre["problem_euler"]),
           "card": str(Path(card).resolve()), "card_sha256": _sha256(card),
           "limits": {k: pre["limits"][k] for k in ("max_nodes_ns", "max_nodes_euler", "approval_file", "approval_sha256", "user_ref", "kind")},
           "nodes": {"ns": pre["nodes_ns"]["n_nodes"], "euler": pre["nodes_euler"]["n_nodes"]}}
    rd_e = out_dir / "euler"
    RA.prepare(pre["problem_euler"], rd_e)
    rd_n = out_dir / "ns"
    RA.prepare_ns(problem, rd_n, initializer=resolve_integral_initializer(problem), stop_rule=False)
    for tag, rd in (("euler", rd_e), ("ns", rd_n)):
        info = json.loads((rd / "prepare_info.json").read_text()) if (rd / "prepare_info.json").is_file() else {}
        m = info.get("mesh") or {}
        rec[tag] = {"run_dir": str(rd), "quality": _quality_line(rd), "n_nodes": m.get("n_nodes"), "max_nodes": m.get("max_nodes"),
                    "hashes": m.get("hashes"), "block": m.get("block")}
    (out_dir / RECORD).write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str) + "\n")
    return rec
