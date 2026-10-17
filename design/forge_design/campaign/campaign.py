"""search の無い campaign `campaign.yaml` (`forge.campaign/v1`; plan tooling-design-orchestrator §4.1.1、§5.1 #4c、親 plan
tooling-design-problem-campaign-recipe の限定の実装 #6a)。

承認した依頼 1 つを、決まった recipe で回すための最小の記述。親 plan の search の無い campaign (`evaluations` を持つ形) とは
「evaluation 1 つ = 承認した依頼 1 つ」で対応する。探索 (`search`)・複数の評価は親 plan の #6b (探索の需要が出てから) なので拒否する。

    format: forge.campaign/v1
    name: m55_trial                         # 依頼の名前 (表示だけ; 識別は台帳の campaign の ID)
    role: fixed                             # fixed (決まった流れだけ・直さない) / fixer (止まった後に §4.4 の範囲で直してよい)
    approval: {file: card_approval.json, sha256: <64 桁>}   # 承認の記録 (case のディレクトリからの相対パス可)
    recipe: {name: axismach.wt_card_chain, version: v1}
    case_dir: .                             # case のディレクトリ (campaign.yaml の所在からの相対パス可)
    ledger: _ledger                         # 台帳の置き場所 (case のディレクトリからの相対パス)
    description: 任意の説明                 # 任意
    estimates:                              # 段の種類ごとの見積もり (plan §4.5.1「段の見積もりの契約」): 依頼と実行のホストに合わせて書く
      euler: {gpu_hours: 0.1, disk_gb: 1.0, basis: "根拠の run・格子の規模・計算と出力の条件・ホスト"}

知らないキー・`search`・`evaluations`・知らない recipe とその版・承認の記録の sha256 の不一致・GPU の段 (`REQUIRED_ESTIMATES`) の見積もりの欠けは ValueError。
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import yaml

FORMAT = "forge.campaign/v1"
REQUIRED = ("format", "name", "role", "approval", "recipe", "case_dir", "ledger")
ROLES = ("fixed", "fixer")      # fixed = 知識の少ないモデル (決まった流れだけ)、fixer = 賢い AI (止まった後に直してよい); 台帳を作った後は変えられない
OPTIONAL = ("description", "estimates")
# 段の種類 (recipe axismach.wt_card_chain/v1 の段の名前の `_p` の前)。見積もりはこの種類ごとに書く
STEP_KINDS = {"axismach.wt_card_chain": ("problem", "mesh", "euler", "ns", "calibrate", "recalibrate", "qualify", "deliverables")}
ESTIMATE_KEYS = ("gpu_hours", "disk_gb", "basis")
# 見積もりが必須の段 (GPU の段; plan §4.5.1「段の見積もりの契約」、2026-10-18 codex 諮問 v3 と判別の A/B
# notes/investigations/2026-10-18-estimate-ab): 無いと起動の前の予約の検査が何も止めない)
REQUIRED_ESTIMATES = {"axismach.wt_card_chain": ("euler", "ns")}
KNOWN_RECIPES = {"axismach.wt_card_chain": ("v1",)}
DEFERRED = {"search": "探索は未対応 (親 plan tooling-design-problem-campaign-recipe §5.1 #6b、探索の需要が出てから)",
            "evaluations": "複数の評価の列は未対応 (1 つの campaign = 承認した依頼 1 つ; 親 plan #6b)",
            "qualification": "判定の条件は承認の記録に置く (campaign には書かない; plan tooling-design-orchestrator §4.7.1)"}


def _sha256(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def validate(doc, base_dir=None) -> dict:
    """campaign の dict を検査して、パスを解決した写しを返す。不正なら理由をまとめて ValueError。"""
    errs = []
    if not isinstance(doc, dict):
        raise ValueError(f"campaign: 最上位は mapping (受け取った値: {type(doc).__name__})")
    for k, why in DEFERRED.items():
        if k in doc:
            errs.append(f"{k}: {why}")
    unknown = [k for k in doc if k not in REQUIRED + OPTIONAL + tuple(DEFERRED)]
    if unknown:
        errs.append(f"知らないキー {unknown} (書けるのは {list(REQUIRED + OPTIONAL)})")
    missing = [k for k in REQUIRED if k not in doc]
    if missing:
        errs.append(f"キー {missing} が無い")
    if doc.get("format") not in (None, FORMAT) and "format" in doc:
        errs.append(f"format は {FORMAT!r} (受け取った値: {doc.get('format')!r})")
    if "name" in doc and not (isinstance(doc["name"], str) and doc["name"].strip()):
        errs.append("name は空でない文字列")
    if "role" in doc and doc["role"] not in ROLES:
        errs.append(f"role は {list(ROLES)} のどれか (受け取った値: {doc['role']!r})")
    rc = doc.get("recipe")
    if "recipe" in doc:
        if not (isinstance(rc, dict) and set(rc) == {"name", "version"}):
            errs.append(f"recipe は {{name, version}} (受け取った値: {rc!r})")
        elif rc["name"] not in KNOWN_RECIPES or rc["version"] not in KNOWN_RECIPES[rc["name"]]:
            errs.append(f"知らない recipe {rc['name']}/{rc['version']} (使えるのは {KNOWN_RECIPES})")
    est = doc.get("estimates")
    need = REQUIRED_ESTIMATES.get((rc or {}).get("name") if isinstance(rc, dict) else None, ())
    lack = [k for k in need if not isinstance(est, dict) or k not in est]
    if lack:
        errs.append(f"estimates に GPU の段 {lack} の見積もり {{gpu_hours, disk_gb, basis}} が無い (依頼と実行のホストに合わせて書く; "
                    "無いと予算の起動の前の検査が止めない、plan §4.5.1)")
    if est is not None:
        kinds = STEP_KINDS.get((rc or {}).get("name") if isinstance(rc, dict) else None, ())
        if not isinstance(est, dict):
            errs.append(f"estimates は段の種類 → {{gpu_hours, disk_gb, basis}} の mapping (受け取った値: {est!r})")
        else:
            for kind, e in est.items():
                if kind not in kinds:
                    errs.append(f"estimates.{kind}: 知らない段の種類 (使えるのは {list(kinds)})")
                elif not isinstance(e, dict) or set(e) != set(ESTIMATE_KEYS):
                    errs.append(f"estimates.{kind} は {{{', '.join(ESTIMATE_KEYS)}}} (受け取った値: {e!r})")
                else:
                    for k in ("gpu_hours", "disk_gb"):
                        v = e[k]
                        if isinstance(v, bool) or not isinstance(v, (int, float)) or not v >= 0:
                            errs.append(f"estimates.{kind}.{k} は 0 以上の数 (受け取った値: {v!r})")
                    if not (isinstance(e["basis"], str) and e["basis"].strip()):
                        errs.append(f"estimates.{kind}.basis (根拠) は空でない文字列")
    base = Path(base_dir) if base_dir else Path.cwd()
    out = dict(doc)
    case_dir = None
    if isinstance(doc.get("case_dir"), str):
        case_dir = (base / doc["case_dir"]).resolve()
        if not case_dir.is_dir():
            errs.append(f"case_dir {case_dir} が無い")
        out["case_dir"] = str(case_dir)
    elif "case_dir" in doc:
        errs.append("case_dir は文字列")
    ap = doc.get("approval")
    if "approval" in doc:
        if not (isinstance(ap, dict) and set(ap) == {"file", "sha256"} and isinstance(ap.get("file"), str)):
            errs.append(f"approval は {{file, sha256}} (受け取った値: {ap!r})")
        elif case_dir is not None:
            f = (case_dir / ap["file"]).resolve()
            if not f.is_file():
                errs.append(f"承認の記録 {f} が無い")
            elif _sha256(f) != ap["sha256"]:
                errs.append(f"承認の記録 {f.name} の sha256 が campaign の値と違う (承認し直したなら campaign も作り直す)")
            out["approval"] = {"file": str(f), "sha256": ap["sha256"]}
    if isinstance(doc.get("ledger"), str) and case_dir is not None:
        out["ledger"] = str((case_dir / doc["ledger"]).resolve())
    elif "ledger" in doc:
        errs.append("ledger は文字列")
    if errs:
        raise ValueError("campaign: " + "; ".join(errs))
    return out


def load(path) -> dict:
    p = Path(path)
    return validate(yaml.safe_load(p.read_text()), base_dir=p.resolve().parent)
