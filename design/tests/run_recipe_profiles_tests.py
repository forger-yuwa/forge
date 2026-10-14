#!/usr/bin/env python3
"""機種ごとのプロファイルの書き写しの試験 (plan tooling-generator-defaults-and-launch-gate §5.1 #6)。

(a) 書き写しの台本 (`solver_density_cuda/tools/recipes/regen_profiles.py`) の出力 = `recommended_settings.yaml` の生成区間
    (生成器の出力が変わったら、ここで落ちる。直し方は `regen_profiles.py --write` で書き直し、差を plan に書く)。
(b) 正例: 各生成器を代表の問題で回した**全段の起動列**の config が、それぞれのプロファイルの段に PASS (FAIL も WARN も無し)。
(c) 負例: 生成器の本段の config を前段 (soft) として宣言すると FAIL (段の不一致を名指し)、別の機種のプロファイルで宣言すると FAIL。
forge・変換器は使わない (段の config は偽の forge で取る)。recipe_settings.py の単体の負例は
solver_density_cuda/tests/unit/test_recipe_settings.py。
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools" / "recipes"))
sys.path.insert(0, str(ROOT / "solver_density_cuda" / "tools"))
import recipe_settings as RS  # noqa: E402
import regen_profiles as RG  # noqa: E402

FAIL = 0


def check(name, cond, detail=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{detail}]" if detail else ""), flush=True)
    if not cond:
        FAIL += 1


text = RG.YAML_PATH.read_text(encoding="utf-8")
recipes = RS.load_recipes(text=text)
profiles, configs = RG.build_profiles(recipes)
cur, _, _ = RG.current_block(text)
check("(a) 書き写しの台本の出力 = recommended_settings.yaml の生成区間", cur == RG.render_block(profiles),
      "違えば regen_profiles.py --write で書き直し、生成器の出力が変わった理由を plan に書く")
check("(a) 全機種を書き写した", set(profiles) == {sp["name"] for sp in RG.SPECS}, sorted(profiles))

for sp in RG.SPECS:
    name = sp["name"]
    stages = configs[name]
    for st, cfg in stages:
        r = RS.check_config(recipes, cfg, name, st, bcond=None)
        check(f"(b) {name}/{st}: 生成器の config がプロファイルに PASS", not r["fails"] and not r["warns"],
              f"fails {r['fails'][:3]} warns {r['warns'][:3]}")
    main = [c for st, c in stages if st == "main"][-1]
    if any(st == "soft" for st, _ in stages):
        r = RS.check_config(recipes, main, name, "soft")
        check(f"(c) {name}: 本段の config を soft として宣言すると FAIL", bool(r["fails"]),
              sorted({k for k, _ in r["fails"]})[:5])
    other = "bell_tp_ns_sst" if name != "bell_tp_ns_sst" else "axismach_ns_sst_tp"
    r = RS.check_config(recipes, main, other, "main")
    check(f"(c) {name}: 別の機種のプロファイル ({other}) で宣言すると FAIL", bool(r["fails"]), sorted({k for k, _ in r["fails"]})[:5])

print(f"\n{'ALL PASS' if FAIL == 0 else f'FAIL {FAIL}'}")
sys.exit(1 if FAIL else 0)
