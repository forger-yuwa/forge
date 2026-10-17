"""承認の記録: 保護する条件の解決・preview・write・check (plan tooling-design-orchestrator §4.7.1 (1)〜(4)、§5.1 #2a)。

AI が変えてはいけない「ユーザが決めた欄」と「判定の条件」を、生成した後の実効の条件 (保護する条件) で決め、ユーザが承認した値を
承認の記録 (`card_approval.json`) に固定する。試行の入力 (較正の係数を足した `card_p1.yaml` など) は承認の記録を参照し、
「試行の入力から解決した保護する条件 = 承認の値」で検査する。格子・較正の係数だけの変更は再承認なしで通り、保護する条件の変更は止まる。

保護する条件 (`protected_conditions`):

| 節 | 中身 | 出どころ |
| --- | --- | --- |
| requirements | 組成・入口条件・M_exit・スロート・入口の形・拡大部の長さ・壁の熱条件・背圧・寸法の公差 | カード (検査した後の値) |
| qualification | V5b・V5c の公差、評価の線 (η・x の窓の規則・点数)、無くてもよい成果物 | カード → 既定 (出所つき) |
| run_policy | 反復の上限・較正し直しの上限 | カード → 既定 |
| mesh_limits | 格子の節点数の上限 (照合は §5.1 #8) | カード → 既定 |
| budget | GPU 時間・ディスク | カード (既定なし; 無ければ preview が問題にする) |
| requirements_without_gate | 判定器の無い要求 | カード → 既定 [] |
| physics | 問題の型・ガスのモデル・輸送・組成・壁の熱条件・物理に効く evaluate のキー | 問題 YAML |
| evaluation | 停止規則のプリセット (中身の sha256)・evaluate.settle・区間・閾値・評価の線 | 問題 YAML + プリセット |

カードの出所 (`sources`) の検査 (`check_sources`): カードに書いたユーザ欄の葉は、最も細かく一致する `sources` の項目が `user` (ref つき) か
`recipe_default` (既定の提供元の値と一致) でなければならない。`proposal`・`derived`・出所なし・既定と違う値・未知のパスは問題。
カードに書かず既定で埋めた欄は `recipe_default` として、提供元と一緒に記録する。

承認の記録は真正性を証明しない (sha256 と user_ref は「指定した承認の記録との整合」だけを示す; 強制は plan §5.1 #14)。

使い方:
  python -m forge_design.intake.approval preview --card CARD --problem PROBLEM --out PREVIEW.json
  python -m forge_design.intake.approval write --preview PREVIEW.json --user-ref "ユーザの発言の引用 (日付)" --out card_approval.json
  python -m forge_design.intake.approval check --approval card_approval.json --card CARD --problem PROBLEM [--run RUN_DIR]
終了コード: 0 OK、1 問題あり (preview の問題・check の不一致)、2 入力の誤り。
"""
from __future__ import annotations

import argparse
import copy
import datetime as _dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import yaml

from . import wt_card as W

FORMAT = "forge.card_approval/v1"
PREVIEW_FORMAT = "forge.card_approval_preview/v1"
PLAN = "plans/active/tooling-design-orchestrator.md §4.7.1"

# ユーザ欄 (カードの最上位の節; plan §4.7.1 (1))。recipe・overrides は AI が選ぶ欄 (保護する条件を変えない範囲で)
REQUIREMENT_SECTIONS = ("gas", "conditions", "M_exit", "throat", "inlet", "divergent_length_m", "wall", "p_ambient", "tolerances")
POLICY_SECTIONS = ("qualification", "run_policy", "mesh_limits", "budget", "requirements_without_gate")
USER_SECTIONS = REQUIREMENT_SECTIONS + POLICY_SECTIONS
# 問題 YAML の evaluate のうち物理に効くキー (runner_axismach が読むもの; 数値の設定 cfl・段・step 数は含めない)
PHYS_EVALUATE_KEYS = ("cfd_gas", "condensation", "thermo_href_temp", "tp_species", "axisym_method", "model")
LIST_LEAF_SECTIONS = ("requirements_without_gate",)


class ApprovalError(ValueError):
    """承認の記録・照合の誤り (wt_qualify の正式な入口はこれで止める)。"""


# ---------------------------------------------------------------------------------------------------- 小道具

def _sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_file(path) -> str:
    return _sha256_bytes(Path(path).read_bytes())


def _canon(o):
    """JSON に通して正規化する (タプル → 列、キーの順を固定)。値の比較と sha256 に使う。"""
    return json.loads(json.dumps(o, sort_keys=True, ensure_ascii=False, default=str))


def sha256_json(o) -> str:
    return _sha256_bytes(json.dumps(_canon(o), sort_keys=True, ensure_ascii=False).encode())


def _leaves(d, prefix: str = "") -> dict:
    """mapping を葉 (パス → 値) にする。列は葉。"""
    out = {}
    if isinstance(d, dict) and d:
        for k, v in d.items():
            p = f"{prefix}.{k}" if prefix else str(k)
            out.update(_leaves(v, p))
    else:
        out[prefix] = d
    return out


def diff_paths(a, b, prefix: str = "") -> list:
    """2 つの (正規化した) 値の違うパスの列。"""
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b), key=str):
            p = f"{prefix}.{k}" if prefix else str(k)
            if k not in a:
                out.append(f"{p} (承認に無い)")
            elif k not in b:
                out.append(f"{p} (今の値に無い)")
            else:
                out += diff_paths(a[k], b[k], p)
        return out
    return [] if a == b else [prefix or "(全体)"]


def _git_state(repo_dir) -> dict:
    """判定器の版の記録 (完全な HEAD・追跡ファイルの変更の有無)。照合は判定専用の環境 (plan §5.1 #2b)。"""
    def run(*args):
        r = subprocess.run(["git", "-C", str(repo_dir), *args], capture_output=True, text=True)
        return r.stdout.strip() if r.returncode == 0 else None
    head = run("rev-parse", "HEAD")
    dirty = run("status", "--porcelain", "--untracked-files=no")
    return {"commit": head, "tracked_clean": (dirty == "") if dirty is not None else None}


# ---------------------------------------------------------------------------------------------------- 既定の提供元

def default_values(prob: dict | None) -> dict:
    """カードに書かなかったユーザ欄の葉の既定: {パス: (値, 提供元)}。提供元はコードの定数・雛形のプリセット。"""
    from ..evaluate import runner_axismach as RA
    from ..evaluate import settle_presets as SP
    from ..feedback import wt_qualify as Q
    out = {"p_ambient": ("design", "wt_card.validate_card の既定 (p_ambient: design = 設計出口静圧)")}
    for k, v in W.DEFAULT_TOLERANCES.items():
        out[f"tolerances.{k}"] = (v, "wt_card.DEFAULT_TOLERANCES")
    for k, (v, src) in Q.DEFAULT_TOLS.items():
        out[f"qualification.{k}"] = (v, f"wt_qualify.DEFAULT_TOLS ({src})")
    es = settle_spec(prob)
    if es is not None:
        P = SP.PRESETS[es["preset"]]
        if "line_eta" in es:
            out["qualification.line_eta"] = (float(es["line_eta"]), "問題 YAML の evaluate.settle.line_eta")
        else:
            out["qualification.line_eta"] = (float(P["line"]["eta"]), f"停止規則のプリセット {es['preset']} の線 η")
        out["qualification.line_n"] = (int(P["line"]["n"]), f"停止規則のプリセット {es['preset']} の線の点数")
        if "x_range_m" in es:
            out["qualification.line_x_range_m"] = ([float(x) for x in es["x_range_m"]], "問題 YAML の evaluate.settle.x_range_m")
        else:
            out["qualification.line_x_range_m"] = ({"rule": "preset", "x_offsets_rt": list(P["line"]["x_offsets_rt"])},
                                                   f"停止規則のプリセット {es['preset']} の規則 (x_E + o0, x_F + o1) r_t × scale_m")
    out["qualification.optional_parts"] = ([], "wt_qualify の既定 (成果物の部品は全部必要)")
    out["run_policy.max_iterations"] = (Q.MAX_ITERATIONS, "wt_qualify.MAX_ITERATIONS (手順書 §8 の --max-iterations 3)")
    out["run_policy.max_recalibrations"] = (Q.MAX_RECALIBRATIONS, "wt_qualify.MAX_RECALIBRATIONS (plan tooling-nozzle-card-local-mesh §4.3)")
    out["mesh_limits.max_nodes_ns"] = (RA.MAX_NODES_2D, "runner_axismach.MAX_NODES_2D")
    out["mesh_limits.max_nodes_euler"] = (RA.MAX_NODES_2D, "runner_axismach.MAX_NODES_2D")
    out["requirements_without_gate"] = ([], "既定 (判定器の無い要求は無い)")
    return out


def settle_spec(prob: dict | None) -> dict | None:
    if not isinstance(prob, dict):
        return None
    es = ((prob.get("evaluate") or {}).get("settle"))
    return copy.deepcopy(es) if isinstance(es, dict) and es.get("preset") else None


# ---------------------------------------------------------------------------------------------------- 出所の検査

def _match_source(path: str, sources: dict):
    """path に最も細かく一致する sources の項目 (パス, 項目)。無ければ (None, None)。"""
    best = None
    for sp in sources:
        if path == sp or path.startswith(sp + "."):
            if best is None or len(sp) > len(best):
                best = sp
    return (best, sources[best]) if best is not None else (None, None)


def raw_user_leaves(raw: dict) -> dict:
    out = {}
    for sec in USER_SECTIONS:
        if sec not in raw:
            continue
        if sec in LIST_LEAF_SECTIONS:
            out[sec] = raw[sec]
        else:
            out.update(_leaves(raw[sec], sec) if isinstance(raw[sec], dict) else {sec: raw[sec]})
    return out


def check_sources(raw: dict, prob: dict | None) -> tuple[list, dict]:
    """戻り: (問題の列, {葉: 出所の記録})。"""
    probs, table = [], {}
    sources = raw.get("sources") or {}
    defaults = default_values(prob)
    leaves = raw_user_leaves(raw)
    for leaf, val in leaves.items():
        sp, ent = _match_source(leaf, sources)
        if ent is None:
            probs.append(f"{leaf}: 出所が無い (sources に user か recipe_default を書く)")
            table[leaf] = {"kind": None}
            continue
        kind = ent.get("kind")
        if kind == "user":
            table[leaf] = {"kind": "user", "ref": ent.get("ref"), "via": sp}
        elif kind == "recipe_default":
            dv = defaults.get(leaf)
            if dv is None:
                probs.append(f"{leaf}: recipe_default だが既定の提供元が無い (ユーザの指定なら kind user)")
            elif _canon(dv[0]) != _canon(val):
                probs.append(f"{leaf}: recipe_default だが値 {val!r} が既定 {dv[0]!r} ({dv[1]}) と違う")
            table[leaf] = {"kind": "recipe_default", "provider": None if dv is None else dv[1], "via": sp}
        else:
            probs.append(f"{leaf}: 出所 {kind!r} は承認に使えない (user か recipe_default; 提案・導出の値はユーザに確かめる)")
            table[leaf] = {"kind": kind, "via": sp}
    known = set(leaves) | set(defaults)
    for sp in sources:
        if sp in ("sources", "approval", "recipe", "overrides", "format", "name") or sp.split(".")[0] not in USER_SECTIONS:
            probs.append(f"sources.{sp}: ユーザ欄でない (書けるのは {list(USER_SECTIONS)} の下)")
            continue
        if not any(k == sp or k.startswith(sp + ".") for k in known):
            probs.append(f"sources.{sp}: 未知のパス (カードにも既定にも無い)")
    for leaf, (v, prov) in defaults.items():
        if leaf not in leaves and not any(k.startswith(leaf + ".") for k in leaves) and not any(leaf.startswith(k + ".") for k in leaves):
            table[leaf] = {"kind": "recipe_default", "provider": prov, "filled": True}
    return probs, table


# ---------------------------------------------------------------------------------------------------- 保護する条件

def _resolved_policy(card: dict, prob: dict | None, section: str) -> dict:
    defaults = default_values(prob)
    got = copy.deepcopy(card.get(section) or {})
    for leaf, (v, _src) in defaults.items():
        sec, _, key = leaf.partition(".")
        if sec == section and key and key not in got:
            got[key] = copy.deepcopy(v)
    return got


def resolved_settle(prob: dict | None) -> dict | None:
    """問題 YAML の evaluate.settle をプリセットで解決した、run に依らない部分 (区間・保持・閾値・評価の線・出口の η の範囲)。"""
    from ..evaluate import settle_presets as SP
    es = settle_spec(prob)
    if es is None:
        return None
    name = es["preset"]
    if name not in SP.PRESETS:
        raise ApprovalError(f"evaluate.settle.preset {name!r} を知らない (プリセット: {list(SP.PRESETS)})")
    P = SP.PRESETS[name]
    th = copy.deepcopy(P["thresholds"])
    for nm, kv in (es.get("thresholds") or {}).items():
        th[nm] = {**th.get(nm, {}), **kv}
    line = {"eta": float(es.get("line_eta", P["line"]["eta"])), "n": int(P["line"]["n"])}
    if "x_range_m" in es:
        line["x_range_m"] = [float(x) for x in es["x_range_m"]]
    else:
        line["x_offsets_rt"] = list(P["line"]["x_offsets_rt"])
    return {"preset": name, "preset_version": P.get("version"), "preset_sha256": SP._sha256_json(P),
            "evaluate_settle": es, "evaluate_settle_sha256": SP._sha256_json(es),
            "window": dict(P["window"]), "hold": int(P["hold"]), "residual": P.get("residual"), "wall_measure": P.get("wall_measure"),
            "exit_eta_range": list(P["exit_eta_range"]), "line": line, "thresholds": th,
            "quantities_remove": list(es.get("quantities_remove") or []), "quantities_add": copy.deepcopy(es.get("quantities_add") or [])}


MESH_AR_MAX_DEFAULT = 1000.0      # AGENTS.md「メッシュ品質チェック」の既定 (壁法線の構造の境界層は 5000 まで; 問題 YAML の mesh.ar_max)
MESH_SKEW_MAX = 0.9


def mesh_quality_spec(prob: dict | None) -> dict | None:
    """計算の格子の品質の閾値 (判定器の閾値なので承認で固定する; plan §4.7.1 (6))。"""
    if not isinstance(prob, dict):
        return None
    ar = (prob.get("mesh") or {}).get("ar_max")
    return {"ar_max": float(ar) if ar is not None else MESH_AR_MAX_DEFAULT, "skew_max": MESH_SKEW_MAX,
            "source": "問題 YAML の mesh.ar_max" if ar is not None else "AGENTS.md の既定 (AR ≤ 1000)"}


def protected_conditions(card: dict, prob: dict | None) -> dict:
    """検査した後のカード (validate_card の戻り) と問題 YAML (dict) から保護する条件を作る。"""
    req = {k: copy.deepcopy(card.get(k)) for k in REQUIREMENT_SECTIONS}
    csv = (card.get("wall") or {}).get("Tw_profile_csv")
    if isinstance(csv, str) and card.get("_card_dir"):
        p = Path(card["_card_dir"]) / csv
        req["wall_Tw_profile_csv_sha256"] = sha256_file(p) if p.is_file() else None
    phys = None
    if isinstance(prob, dict):
        g = prob.get("gas") or {}
        ev = prob.get("evaluate") or {}
        phys = {"type": prob.get("type"), "gas_model": g.get("model"), "gas_transport": g.get("transport"),
                "gas_species": g.get("species"), "gas_composition_basis": g.get("composition_basis"),
                "wall_thermal": (prob.get("spec") or {}).get("wall_thermal"),
                "evaluate": {k: ev[k] for k in PHYS_EVALUATE_KEYS if k in ev}}
    out = {"requirements": req,
           "qualification": _resolved_policy(card, prob, "qualification"),
           "run_policy": _resolved_policy(card, prob, "run_policy"),
           "mesh_limits": _resolved_policy(card, prob, "mesh_limits"),
           "budget": copy.deepcopy(card.get("budget")),
           "requirements_without_gate": copy.deepcopy(card.get("requirements_without_gate") or []),
           "physics": phys,
           "evaluation": {"settle": resolved_settle(prob), "mesh_quality": mesh_quality_spec(prob)}}
    return _canon(out)


def problem_matches_card(card: dict, prob: dict) -> list:
    """問題 YAML がカードから作られたか (入口条件・M・壁の熱条件・組成)。違いの列。"""
    probs = []
    sp = prob.get("spec") or {}
    for k_card, k_prob in (("Pt", "Pt"), ("Tt", "Tt")):
        a, b = (card.get("conditions") or {}).get(k_card), sp.get(k_prob)
        if a is None or b is None or float(a) != float(b):
            probs.append(f"spec.{k_prob} {b!r} がカードの conditions.{k_card} {a!r} と違う")
    if sp.get("M_design") is None or float(sp["M_design"]) != float(card["M_exit"]):
        probs.append(f"spec.M_design {sp.get('M_design')!r} がカードの M_exit {card['M_exit']!r} と違う")
    sp_card = W.card_species(card)
    sp_prob = (prob.get("gas") or {}).get("species")
    if not isinstance(sp_prob, dict) or {k.upper(): float(v) for k, v in sp_prob.items()} != {k.upper(): float(v) for k, v in sp_card.items()}:
        probs.append(f"gas.species {sp_prob!r} がカードの組成 {sp_card!r} と違う")
    w = card.get("wall") or {}
    wt = sp.get("wall_thermal") or {}
    if w.get("thermal") == "isothermal" and "Tw" in w:
        if wt.get("mode") != "isothermal" or wt.get("Tw") is None or float(wt["Tw"]) != float(w["Tw"]):
            probs.append(f"spec.wall_thermal {wt!r} がカードの壁 {w!r} と違う")
    elif w.get("thermal") == "adiabatic" and wt and wt.get("mode") not in (None, "adiabatic"):
        probs.append(f"spec.wall_thermal {wt!r} がカードの断熱壁と違う")
    return probs


def _load(card_path, problem_path):
    raw = yaml.safe_load(Path(card_path).read_text())
    card = W.load_card(card_path)
    card["_card_dir"] = str(Path(card_path).resolve().parent)
    prob = yaml.safe_load(Path(problem_path).read_text()) if problem_path else None
    return raw, card, prob


# ---------------------------------------------------------------------------------------------------- preview・write・check

def preview(card_path, problem_path, repo_dir=None, judge_env: bool = True, judge_python=None) -> dict:
    """承認の前に見せる表。problems が空でなければ write できない。judge_env: 判定専用の環境 (判定の Python・外部の道具) を記録する
    (plan §4.7.1 (5); 既定は今の Python)。作業ツリーの追跡ファイルに変更があれば問題にする (承認の commit を固定するため)。"""
    raw, card, prob = _load(card_path, problem_path)
    probs = []
    if card["format"] != W.CARD_FORMAT_V2:
        probs.append(f"カードの format が {W.CARD_FORMAT_V2!r} でない ({card['format']!r}) — 承認には版 2 が要る")
    if "approval" in raw:
        probs.append("承認の原本のカードに approval (承認の記録への参照) がある — 参照は試行の入力にだけ書く")
    sprobs, table = check_sources(raw, prob)
    probs += sprobs
    if prob is None:
        probs.append("問題 YAML が無い (物理モデル・評価の条件を決められない)")
    else:
        probs += problem_matches_card(card, prob)
        if settle_spec(prob) is None:
            probs.append("問題 YAML に evaluate.settle (停止規則のプリセット) が無い")
        q = prob.get("qualification")
        if q:
            probs.append(f"問題 YAML に qualification がある ({q!r}) — 判定の条件はカードの qualification に書く")
    if not card.get("budget"):
        probs.append("budget (GPU 時間・ディスクの承認値) が無い — 正式な入口で必須 (design-intake §8 の見積もり)")
    prot = protected_conditions(card, prob)
    repo = Path(repo_dir) if repo_dir else Path(__file__).resolve().parents[3]
    gate = _git_state(repo)
    if gate.get("tracked_clean") is not True:
        probs.append("判定器の版: 作業ツリーの追跡ファイルに変更がある — commit してから preview する (承認の commit を固定するため)")
    if judge_env:
        from ..feedback import gate_env as GE
        try:
            gate["env"] = GE.snapshot(judge_python)
        except GE.GateEnvError as e:
            probs.append(f"判定専用の環境を記録できない ({e})")
            gate["env"] = None
        else:
            if gate["env"]["python"].get("editable_dists"):
                probs.append(f"判定の Python に editable の配布物がある {gate['env']['python']['editable_dists']}")
    return {"format": PREVIEW_FORMAT, "plan": PLAN, "created": _dt.datetime.now().isoformat(timespec="seconds"),
            "card": {"path": str(card_path), "name": card["name"], "format": card["format"], "sha256": sha256_file(card_path),
                     "text": Path(card_path).read_text()},
            "problem": None if problem_path is None else {"path": str(problem_path), "sha256": sha256_file(problem_path)},
            "protected": prot, "protected_sha256": sha256_json(prot), "sources": table, "problems": probs,
            "gate": gate, "gate_sha256": sha256_json(gate)}


def preview_table(pv: dict) -> list[str]:
    """ユーザに見せる表 (葉・値・出所)。"""
    L = [f"承認の対象: カード {pv['card']['name']} ({pv['card']['path']}, sha256 {pv['card']['sha256'][:16]}…)",
         f"判定器の版: {pv['gate'].get('commit')} (追跡ファイルの変更なし: {pv['gate'].get('tracked_clean')})", "",
         "| 欄 | 値 | 出所 |", "| --- | --- | --- |"]
    prot = pv["protected"]
    rows = {}
    rows.update(_leaves(prot["requirements"]))
    for sec in POLICY_SECTIONS:
        v = prot.get(sec)
        if sec in LIST_LEAF_SECTIONS or not isinstance(v, dict):
            rows[sec] = v
        else:
            rows.update(_leaves(v, sec))
    for k, v in rows.items():
        s = pv["sources"].get(k) or next((pv["sources"][p] for p in pv["sources"] if k.startswith(p + ".")), None) or {}
        src = s.get("kind") or "—"
        if s.get("ref"):
            src += f": {s['ref']}"
        if s.get("provider"):
            src += f" ({s['provider']})"
        L.append(f"| {k} | {json.dumps(v, ensure_ascii=False)} | {src} |")
    ph = prot.get("physics") or {}
    L += ["", f"物理モデル (問題 YAML): {json.dumps(ph, ensure_ascii=False)}"]
    st = (prot.get("evaluation") or {}).get("settle") or {}
    L.append(f"停止規則: プリセット {st.get('preset')} (sha256 {str(st.get('preset_sha256'))[:16]}…)、区間 {st.get('window')}、"
             f"評価の線 {st.get('line')}")
    if pv["problems"]:
        L += ["", "**承認できない (直してから preview し直す)**:"] + [f"- {p}" for p in pv["problems"]]
    return L


def write(preview_path, user_ref: str, out_path) -> dict:
    """preview のファイルの内容をそのまま承認の記録にする (表示の後に解決し直さない)。既にあれば止める。"""
    out_path = Path(out_path)
    if out_path.exists():
        raise ApprovalError(f"{out_path} が既にある — 承認の記録は上書きしない (承認し直すなら別の名前で)")
    if not (isinstance(user_ref, str) and user_ref.strip()):
        raise ApprovalError("user_ref (ユーザの承認の発言の引用と日付) が空")
    pv = json.loads(Path(preview_path).read_text())
    if pv.get("format") != PREVIEW_FORMAT:
        raise ApprovalError(f"{preview_path} は preview でない (format {pv.get('format')!r})")
    if pv["problems"]:
        raise ApprovalError("preview に問題がある — 承認できない: " + "; ".join(pv["problems"]))
    if sha256_json(pv["protected"]) != pv["protected_sha256"] or sha256_json(pv["gate"]) != pv.get("gate_sha256"):
        raise ApprovalError("preview の protected・gate と sha256 が合わない (preview を書き換えた?)")
    rec = {"format": FORMAT, "plan": PLAN, "approved": _dt.datetime.now().isoformat(timespec="seconds"), "user_ref": user_ref,
           "preview_sha256": sha256_file(preview_path), "card": pv["card"], "problem": pv["problem"], "protected": pv["protected"],
           "protected_sha256": pv["protected_sha256"], "sources": pv["sources"], "gate": pv["gate"], "gate_sha256": pv["gate_sha256"],
           "note": "承認の記録は真正性を証明しない (指定した承認の記録との整合だけを示す; plan §4.7.1 (4))"}
    out_path.write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    return rec


def load_approval(path) -> dict:
    ap = json.loads(Path(path).read_text())
    if ap.get("format") != FORMAT:
        raise ApprovalError(f"{path} は承認の記録でない (format {ap.get('format')!r})")
    if sha256_json(ap["protected"]) != ap.get("protected_sha256") or sha256_json(ap.get("gate")) != ap.get("gate_sha256"):
        raise ApprovalError(f"{path} の protected・gate と sha256 が合わない (承認の記録を書き換えた?)")
    return ap


def check(approval_path, card_path, problem_path) -> dict:
    """カード (承認の原本か、承認の記録を参照する試行の入力) と問題 YAML が承認の値と一致するか。
    戻り: {status: OK / MISMATCH, kind: original / trial, diffs, problems}。"""
    ap = load_approval(approval_path)
    raw, card, prob = _load(card_path, problem_path)
    probs = []
    sha = sha256_file(card_path)
    if sha == ap["card"]["sha256"]:
        kind = "original"
    else:
        kind = "trial"
        ref = raw.get("approval")
        if not isinstance(ref, dict):
            probs.append("承認の原本と違うカードに approval (承認の記録への参照) が無い — 試行の入力は {file, sha256} で承認の記録を指す")
        elif ref.get("sha256") != sha256_file(approval_path):
            probs.append(f"カードの approval.sha256 {str(ref.get('sha256'))[:16]}… が承認の記録 {Path(approval_path).name} の sha256 と違う")
    if prob is not None:
        probs += problem_matches_card(card, prob)
        if prob.get("qualification"):
            probs.append(f"問題 YAML に qualification がある ({prob['qualification']!r}) — 判定の条件は承認の記録だけ")
    prot = protected_conditions(card, prob)
    diffs = diff_paths(ap["protected"], prot)
    return {"status": "OK" if not diffs and not probs else "MISMATCH", "kind": kind, "card_sha256": sha,
            "approval_sha256": sha256_file(approval_path), "diffs": diffs, "problems": probs}


def _norm_quantity(q: dict) -> dict:
    keep = ("name", "kind", "thresholds", "eta", "n", "x_range_m", "eta_range", "Md", "required", "inlet", "outlet", "walls", "measure")
    return {k: q[k] for k in keep if k in q}


def check_run(approval, run_dir) -> dict:
    """run の停止規則の基準 (settle_criteria.yaml) が承認した評価の条件から作られたか (plan §4.7.1 (1)、§5.1 #2a の判別 A/B (ii))。
    承認した evaluate.settle を今のコードで run の設定に解決し直した基準と、run の基準を比べる。プリセットの中身が承認の後に変わったら
    preset_sha256 の不一致で止まる。戻り: {status: OK / MISMATCH, diffs, problems}。"""
    from ..evaluate import settle_presets as SP
    ap = load_approval(approval) if not isinstance(approval, dict) else approval
    run_dir = Path(run_dir)
    st = (ap["protected"].get("evaluation") or {}).get("settle")
    probs, diffs = [], []
    if not st:
        return {"status": "MISMATCH", "diffs": [], "problems": ["承認の記録に停止規則の条件が無い"]}
    cf = run_dir / SP.criteria_file_name()
    if not cf.is_file():
        return {"status": "MISMATCH", "diffs": [], "problems": [f"{cf.name} が無い"]}
    crit = yaml.safe_load(cf.read_text())
    src = crit.get("source") or {}
    for k in ("preset", "preset_sha256", "evaluate_settle_sha256"):
        if src.get(k) != st.get(k):
            probs.append(f"基準の source.{k} {str(src.get(k))[:24]!r} が承認の値 {str(st.get(k))[:24]!r} と違う")
    if SP._sha256_json(SP.PRESETS.get(st["preset"], {})) != st["preset_sha256"]:
        probs.append(f"今のコードのプリセット {st['preset']} の中身が承認の後に変わった (sha256 が違う)")
    info = json.loads((run_dir / "prepare_info.json").read_text())
    md_ap = (ap["protected"].get("requirements") or {}).get("M_exit")
    if md_ap is None or float(info.get("Md")) != float(md_ap):
        probs.append(f"prepare_info の Md {info.get('Md')!r} が承認の M_exit {md_ap!r} と違う")
    try:
        scfg = yaml.safe_load((run_dir / "solverConfig.yaml").read_text())
        bcfg = yaml.safe_load((run_dir / "bcondConfig.yaml").read_text())
        exp = SP.resolve(st["evaluate_settle"], info, scfg, bcfg, segment=(crit.get("grid") or {}).get("segment"))
    except (OSError, ValueError, KeyError) as e:
        return {"status": "MISMATCH", "diffs": [], "problems": probs + [f"承認した evaluate.settle を run の設定で解決できない ({e})"]}
    for k in ("window", "hold"):
        if _canon(crit.get(k)) != _canon(exp.get(k)):
            diffs.append(f"{k}: run {crit.get(k)!r} / 承認から {exp.get(k)!r}")
    qa = {q["name"]: _norm_quantity(q) for q in crit.get("quantities") or []}
    qe = {q["name"]: _norm_quantity(q) for q in exp.get("quantities") or []}
    for nm in sorted(set(qa) | set(qe)):
        if nm not in qa:
            diffs.append(f"quantities.{nm}: run の基準に無い")
        elif nm not in qe:
            diffs.append(f"quantities.{nm}: 承認の条件に無い")
        else:
            diffs += [f"quantities.{nm}.{d}" for d in diff_paths(_canon(qe[nm]), _canon(qa[nm]))]
    ex_a = sorted(x.get("name") for x in crit.get("excluded") or [])
    ex_e = sorted(x.get("name") for x in exp.get("excluded") or [])
    if ex_a != ex_e:
        diffs.append(f"excluded: run {ex_a} / 承認から {ex_e}")
    ql = st["line"]
    lq = qa.get("line_mach") or {}
    if lq and (int(lq.get("n", -1)) != int(ql["n"]) or float(lq.get("eta", -1)) != float(ql["eta"])):
        diffs.append(f"評価の線: run (η {lq.get('eta')}, 点数 {lq.get('n')}) / 承認 (η {ql['eta']}, 点数 {ql['n']})")
    return {"status": "OK" if not diffs and not probs else "MISMATCH", "diffs": diffs, "problems": probs}


# ---------------------------------------------------------------------------------------------------- CLI

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m forge_design.intake.approval",
                                 description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog="plan tooling-design-orchestrator §4.7.1。AI はユーザに preview の表を見せて承認を得てから write する。")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("preview", help="保護する条件と出所の表を作る (承認の前に見せる)")
    p.add_argument("--card", required=True)
    p.add_argument("--problem", required=True, help="カードから生成した問題 YAML")
    p.add_argument("--out", required=True, help="preview の JSON (既にあれば止める)")
    w = sub.add_parser("write", help="preview の内容をそのまま承認の記録にする")
    w.add_argument("--preview", required=True)
    w.add_argument("--user-ref", required=True, help="ユーザの承認の発言の引用と日付")
    w.add_argument("--out", required=True, help="承認の記録 (card_approval.json; 既にあれば止める)")
    c = sub.add_parser("check", help="カード・問題 YAML (・run の停止規則の基準) が承認の値と一致するか")
    c.add_argument("--approval", required=True)
    c.add_argument("--card", required=True)
    c.add_argument("--problem", required=True)
    c.add_argument("--run", default=None, help="run のディレクトリ (settle_criteria.yaml を承認の評価の条件と照合)")
    a = ap.parse_args(argv)
    try:
        if a.cmd == "preview":
            out = Path(a.out)
            if out.exists():
                raise ApprovalError(f"{out} が既にある")
            pv = preview(a.card, a.problem)
            out.write_text(json.dumps(pv, indent=1, ensure_ascii=False))
            print("\n".join(preview_table(pv)))
            print(f"\n書いた: {out}")
            return 0 if not pv["problems"] else 1
        if a.cmd == "write":
            rec = write(a.preview, a.user_ref, a.out)
            print(f"書いた: {a.out} (protected sha256 {rec['protected_sha256'][:16]}…)")
            return 0
        r = check(a.approval, a.card, a.problem)
        if a.run:
            rr = check_run(a.approval, a.run)
            r["run"] = rr
            if rr["status"] != "OK":
                r["status"] = "MISMATCH"
        print(json.dumps(r, indent=1, ensure_ascii=False))
        print(f"VERDICT: {r['status']}")
        return 0 if r["status"] == "OK" else 1
    except (ApprovalError, ValueError, OSError) as e:
        print(f"approval: 止めた — {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
