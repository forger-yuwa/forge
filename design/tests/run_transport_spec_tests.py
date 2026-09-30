#!/usr/bin/env python3
"""problem YAML `gas.transport` (種ごとの輸送物性の指定) の検査と NS config への翻訳 (plan thermophysics-solver-owned-species-db #9b)。

- semiperfect TP の NS/SST config: `viscMethod: 2` + `physProp.transport` (lump 構成種を含む全実種)、`thermCondMethod` なし
- `gas.transport` が無い TP の NS は prepare でエラー (必要な実種と書き方の例を示す)
- 指定漏れ・余分な種名・lump 名・未知モデル・custom の誤用・重複・cpg での指定は load_problem で拒否
- Euler の config と species_meta は `gas.transport` の有無で変わらない、CPG / `cfd_gas: cpg` の NS は Sutherland のまま
design_chain・メッシュ・ソルバは使わない (config 生成関数を直接呼ぶ)。
"""
import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "design"))
from forge_design.evaluate import runner_axismach as A  # noqa: E402
from forge_design.evaluate import runner_wt as W  # noqa: E402
from forge_design.probdef import load_problem  # noqa: E402

FAIL = 0
NS42 = ROOT / "case/42.isobutane_wt/problem_ib_R2_LU6_Lc8_h2o5_split_ns.yaml"      # split_h2o: [MIXDRY(N2,CO2,O2), H2O]
VA3 = ROOT / "case/44.vitiated_air_wt/problem_va3_M4.19_Lc8_dry_lumpX.yaml"          # lump に AR を含む
CPG42 = ROOT / "case/42.isobutane_wt/problem_ib_R1.5_LU4_Lc6.yaml"                   # semiperfect + cfd_gas: cpg
CPG41 = ROOT / "case/41.wind_tunnel_design/problem_m4_axismach.yaml"                 # gas.model cpg
TMP = Path(tempfile.mkdtemp(prefix="transport_spec_"))


def _chk(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""))
    if not cond:
        FAIL += 1


def with_transport(src, tr, tag, model=None):
    raw = yaml.safe_load(Path(src).read_text())
    if tr is None:
        raw["gas"].pop("transport", None)
    else:
        raw["gas"]["transport"] = tr
    if model is not None:
        raw["gas"]["model"] = model
    p = TMP / f"{tag}.yaml"
    p.write_text(yaml.safe_dump(raw, sort_keys=False, allow_unicode=True))
    return p


def rejects(name, fn, must=()):
    try:
        fn()
    except ValueError as ex:
        msg = str(ex)
        _chk(name, all(m in msg for m in must), msg.splitlines()[0][:160])
        return
    _chk(name, False, "拒否されなかった")


def ns_cfg(p, tag):
    d = TMP / tag
    d.mkdir()
    return A._apply_gas_to_config(W._config_sst_node(p, 100, 50, 1.0), p, d, viscous=True), d


def euler_cfg(p, tag):
    d = TMP / tag
    d.mkdir()
    return A._apply_gas_to_config(W._config_euler_node(p, 100, 50, 4.0, 1), p, d), d


GOOD = {"N2": "cea", "CO2": "cea", "O2": "kinetic", "H2O": "custom:h2o_iapws_cea_v1"}

# (a) 正例: NS config に全実種の指定、viscMethod 2、thermCondMethod なし、species_meta に来歴
p = load_problem(with_transport(NS42, GOOD, "good"))
cfg, d = ns_cfg(p, "good_ns")
pp = yaml.safe_load(cfg)["physProp"]
_chk("NS: viscMethod 2", pp.get("viscMethod") == 2)
_chk("NS: physProp.transport = lump 構成種を含む全実種 (輸送種の実種順)",
     pp.get("transport") == {"N2": "cea", "CO2": "cea", "O2": "kinetic", "H2O": "custom:h2o_iapws_cea_v1"}
     and list(pp["transport"]) == ["N2", "CO2", "O2", "H2O"], str(pp.get("transport")))
_chk("NS: thermCondMethod を残さない (viscMethod 2 では読まれない)", "thermCondMethod" not in pp)
_chk("NS: 必須キー visc / thermCond と prandtlLam は残る", all(k in pp for k in ("visc", "thermCond", "prandtlLam")))
lump_members = [k for it in pp["species"] if isinstance(it, dict) for k in it["lump"]]
_chk("NS: transport の綴り = physProp.species の lump 構成種の綴り", set(lump_members) <= set(pp["transport"]), str(lump_members))
meta = yaml.safe_load((d / "species_meta.yaml").read_text())
_chk("species_meta.yaml に輸送指定の来歴", meta.get("transport", {}).get("models") == pp["transport"])

# (b) 大小文字違いの綴りは設計側の規則 (大文字化) で同じ種になる
p = load_problem(with_transport(VA3, {"N2": "cea", "O2": "cea", "Ar": "cea", "CO2": "cea", "H2O": "cea"}, "alias"))
cfg, _ = ns_cfg(p, "alias_ns")
_chk("`Ar` は AR として受理し config には AR で書く", '"AR": "cea"' in cfg)

# (c) Euler・CPG は gas.transport の有無で変わらない
p0 = load_problem(with_transport(NS42, None, "none"))
e0, d0 = euler_cfg(p0, "none_eu")
e1, d1 = euler_cfg(load_problem(with_transport(NS42, GOOD, "good2")), "good_eu")
_chk("Euler config は gas.transport の有無で同一 (transport を書かない)", e0 == e1 and "transport" not in e1)
_chk("Euler の species_meta は gas.transport の有無で同一",
     (d0 / "species_meta.yaml").read_bytes() == (d1 / "species_meta.yaml").read_bytes())
for src, tag in ((CPG42, "cfdcpg"), (CPG41, "cpg")):
    q = load_problem(src)
    c = A._apply_gas_to_config(W._config_sst_node(q, 100, 50, 1.0), q, TMP, viscous=True)
    _chk(f"{tag}: NS は Sutherland のまま (生テンプレートと同一)", c == A._apply_gas_to_config(W._config_sst_node(q, 100, 50, 1.0), q, TMP)
         and "viscMethod: 1," in c and "transport" not in c)

# (d) 負例
rejects("TP の NS で gas.transport が無い → 実種と例を示して拒否",
        lambda: ns_cfg(p0, "none_ns"), must=("gas.transport が必要", "N2", "H2O", "transport: {"))
rejects("prepare_ns も run dir を作る前に拒否",
        lambda: A.prepare_ns(with_transport(NS42, None, "none3"), TMP / "never"), must=("gas.transport が必要",))
_chk("prepare_ns の拒否で run dir を作らない", not (TMP / "never").exists())
rejects("指定漏れ (H2O)", lambda: load_problem(with_transport(NS42, {"N2": "cea", "CO2": "cea", "O2": "cea"}, "miss")),
        must=("指定の無い実種 ['H2O']",))
rejects("余分な種名 (AR は組成に無い)", lambda: load_problem(with_transport(NS42, dict(GOOD, AR="cea"), "extra")),
        must=("輸送種に無い種 ['AR']",))
rejects("lump 名 (MIXDRY) を書く", lambda: load_problem(with_transport(NS42, dict(GOOD, MIXDRY="cea"), "lumpname")),
        must=("MIXDRY",))
rejects("未知モデル", lambda: load_problem(with_transport(NS42, dict(GOOD, N2="sutherland"), "unk")), must=("未知",))
rejects("モデル名の大小文字違い (CEA)", lambda: load_problem(with_transport(NS42, dict(GOOD, N2="CEA"), "upper")), must=("未知",))
rejects("custom:h2o_iapws_cea_v1 を H2O 以外に", lambda: load_problem(with_transport(NS42, dict(GOOD, N2="custom:h2o_iapws_cea_v1"), "cust")),
        must=("H2O 専用",))
rejects("未知の custom", lambda: load_problem(with_transport(NS42, dict(GOOD, H2O="custom:h2o_iapws_cea_v2"), "cust2")), must=("未知",))
rejects("同じ種を大小文字違いで 2 回", lambda: load_problem(with_transport(NS42, dict(GOOD, o2="cea"), "dup")), must=("2 回",))
rejects("空の mapping", lambda: load_problem(with_transport(NS42, {}, "empty")), must=("空でない",))
rejects("gas.model cpg で gas.transport", lambda: load_problem(with_transport(CPG41, {"N2": "cea"}, "cpgtr")),
        must=("semiperfect | frozen_tp でのみ有効",))
# frozen_tp (SERN、2026-09-30 に lump 記法へ切り替え = R8) は gas.transport を受け付け、NS の config が viscMethod 2 + transport になる
try:
    from forge_design.evaluate import runner_sern as _R2
    _ps = load_problem(ROOT / "case/46.sern_design/problem_moo_frozen_tp_cycle3op_transport.yaml")
    _R2.design_snapshot(_ps); _R2.select_operating_point(_ps, "m6_on")
    _cfg = _R2._solver_config(_ps, 10, 10, 1.0, 2851.0)
    _ok = ("viscMethod: 2" in _cfg and 'transport: {"N2": "cea", "H2O": "custom:h2o_iapws_cea_v1"' in _cfg
           and "thermCondMethod" not in _cfg)
    _p0 = load_problem(ROOT / "case/46.sern_design/problem_moo_frozen_tp_cycle3op.yaml")
    _R2.design_snapshot(_p0); _R2.select_operating_point(_p0, "m6_on")
    _ok = _ok and "viscMethod: 1" in _R2._solver_config(_p0, 10, 10, 1.0, 2851.0)
    _chk("frozen_tp (SERN): gas.transport があれば viscMethod 2 + 全実種の transport、無ければ viscMethod 1 のまま", _ok)
except Exception as e:
    _chk(f"frozen_tp (SERN) の gas.transport ({e})", False)

import shutil  # noqa: E402
shutil.rmtree(TMP, ignore_errors=True)
print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
