"""初期場の固定・継承する値の解決・適合性 (plan tooling-design-orchestrator §4.1.1 の 1、§5.1 #4a)。

実行役は準備の**前**に初期場を決め、計画 (`execution_plan.json`) に書く。準備 (`runner_axismach.prepare`・`prepare_ns`) には
固定したファイル (`ic_file`) と sha256 (`ic_sha256`) を渡し、元の run を準備の中で読み直さない。

- `resolve_ic(ic_file, problem_path, for_ns)`: 初期場のファイル (パス・step・sha256)、元の run の種類 (NS / Euler)、
  NS の本段が継承するリミッタの基準値 (`limiter_ref_inherit_from_ic` を計画の時点で呼び、明示の値 `{ro, p, a, L}` に直す)、
  元の run の記録 (準備の記録・本段の設定の sha256) を返す。
- 適合性 (`compat_problems`): 元の run と今の問題で、熱力学のモデル (thermalMethod)・エネルギー基準 (thermoHrefTemp)・
  輸送種の配置 (physProp.species の列、まとめた種の中身まで)・設計マッハ数が同じであること。違えば止める (格子は cross-mesh の補間が前提なので一致を求めず記録だけ)。

試験: design/tests/run_campaign_ic_tests.py。
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import yaml


def _sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _limiter_explicit(refs: dict) -> dict | str:
    """`limiter_ref_inherit_from_ic` の戻り ({space のキー: 値}) を問題 YAML の明示の形 {ro, p, a, L} に直す。空 = 自動 ("auto")。"""
    from ..evaluate import runner_axismach as RA
    if not refs:
        return RA.LIMITER_REF_AUTO
    inv = {sk: k for k, sk in RA._LIMITER_REF_KEYS + RA._LIMITER_REF_OPT_KEYS}
    out = {}
    for sk, v in refs.items():
        if sk not in inv:
            raise ValueError(f"継承したリミッタの基準値に知らないキー {sk!r}")
        out[inv[sk]] = float(v)
    return out


def compat_problems(src_run: Path, p) -> list:
    """元の run と今の問題の適合性 (熱力学のモデル・エネルギー基準・輸送種の配置・設計マッハ数)。違いの列。"""
    from ..evaluate import runner_axismach as RA
    probs = []
    try:
        cfg = yaml.safe_load((src_run / "solverConfig.yaml").read_text()) or {}
        info = json.loads((src_run / "prepare_info.json").read_text())
    except (OSError, ValueError) as e:
        return [f"元の run {src_run} の solverConfig.yaml・prepare_info.json を読めない ({e})"]
    pp = cfg.get("physProp") or {}
    tp_now = bool(p.is_semiperfect) and str(p.evaluate.get("cfd_gas", "same")) != "cpg"
    tm_src = pp.get("thermalMethod")
    if (tm_src == 2) != tp_now:
        probs.append(f"熱力学のモデル: 元の run の thermalMethod {tm_src!r}、今の問題は {'TP (2)' if tp_now else 'CPG (0)'}")
    if tp_now:
        h_now = float(p.evaluate.get("thermo_href_temp", 298.15))
        h_src = pp.get("thermoHrefTemp")
        if h_src is None or abs(float(h_src) - h_now) > 1e-9:
            probs.append(f"エネルギー基準: 元の run の thermoHrefTemp {h_src!r}、今の問題は {h_now!r}")
        # 輸送種の配置は元の run の solverConfig の physProp.species (NS の準備の記録は species を持たない) と、今の問題から書き手と同じ規則
        # (`solver_species_config`) で作った列を、まとめた種の中身 (モル分率) まで含めて比べる
        from ..gas.composition import solver_species_config
        items_now, _ = solver_species_config(p.species_layout())
        if json.dumps(pp.get("species"), sort_keys=True) != json.dumps(items_now, sort_keys=True):
            probs.append(f"輸送種の配置: 元の run {pp.get('species')!r}、今の問題 {items_now!r}")
    md_src = info.get("Md")
    if md_src is None or abs(float(md_src) - float(p.spec["M_design"])) > 1e-12:
        probs.append(f"設計マッハ数: 元の run の Md {md_src!r}、今の問題の M_design {p.spec.get('M_design')!r}")
    return probs


def resolve_ic(ic_file, problem_path, for_ns: bool) -> dict:
    """計画の時点で初期場を解決する。適合しなければ ValueError。戻り:
    {snapshot: {path, step, sha256, how}, source_run, kind, limiter_ref (NS の本段の明示の値 or "auto" or None), limiter_ref_source,
     source_record: {prepare_info_sha256, solverConfig_sha256, mesh_hashes}, compat: "OK"}。"""
    from ..evaluate import runner_axismach as RA
    from ..probdef import load_problem
    f, snap = RA.select_ic(None, ic_file, None)
    src = f.resolve().parent
    p = load_problem(problem_path)
    probs = compat_problems(src, p)
    if probs:
        raise ValueError("初期場が今の問題に適合しない: " + "; ".join(probs))
    kind = RA._ic_source_kind(src)
    out = {"snapshot": snap, "source_run": str(src), "kind": kind["kind"], "kind_by": kind.get("by"), "compat": "OK",
           "source_record": {"prepare_info_sha256": _sha256_file(src / "prepare_info.json"),
                             "solverConfig_sha256": _sha256_file(src / "solverConfig.yaml"),
                             "mesh_hashes": (json.loads((src / "prepare_info.json").read_text()).get("mesh") or {}).get("hashes")}}
    out["limiter_ref"], out["limiter_ref_source"] = None, None
    if for_ns:
        ml = (p.evaluate.get("main_line") or None) if isinstance(p.evaluate, dict) else None
        if isinstance(ml, dict) and ml.get("limiter_ref") is None:
            refs, rec = RA.limiter_ref_inherit_from_ic(src)
            out["limiter_ref"], out["limiter_ref_source"] = _limiter_explicit(refs), rec
        elif isinstance(ml, dict):
            out["limiter_ref_source"] = {"mode": "problem", "value": ml.get("limiter_ref")}
    return out
