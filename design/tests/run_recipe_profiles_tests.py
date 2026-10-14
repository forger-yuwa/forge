#!/usr/bin/env python3
"""機種ごとのプロファイルの書き写しの試験 (plan tooling-generator-defaults-and-launch-gate §5.1 #6・#16)。

(a) 書き写しの台本 (`solver_density_cuda/tools/recipes/regen_profiles.py`) の出力 = `recommended_settings.yaml` の生成区間
    (生成器の出力が変わったら、ここで落ちる。直し方は `regen_profiles.py --write` で書き直し、差を plan に書く)。
(b) 各生成器を代表の問題で回した**全段の起動列**の config を、それぞれのプロファイルの段に照合する。#16 (2026-10-15) から
    `fixed`・`explicit` の省略は FAIL なので、今の生成器の config は**その段の明示すべきキー (`explicit`) の省略だけ**で FAIL する
    (これは #11 で生成器に書かせる例外表のキー。黙って通さず、FAIL のキーの集合 = 例外表であることを確かめる)。
    例外表のキーを旧い実効値で足した config (#11 で生成器が書く形) は PASS (FAIL も WARN も無し)。
(c) 負例: 生成器の本段の config を前段 (soft) として宣言すると FAIL (段の不一致を名指し)、別の機種のプロファイルで宣言すると FAIL。
(d) #16 の判別 A/B (2026-10-15 codex diagnose、forge 0 step): SERN 2D 本段の生成器の config (例外表のキーを足した形) で
    A = `space.limiterScaled: 1` を明記 → PASS、B = そのキーだけ削除 → 同じキーを名指しで FAIL。
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
EXCEPTIONS = 0


def with_keys(cfg, add=None, drop=()):
    """config の文字列に完全修飾パスのキーを足す・消す (#11 の例外表の形を作る)。"""
    import yaml
    doc = yaml.safe_load(cfg)
    for k, v in (add or {}).items():
        cur, parts = doc, k.split(".")
        for p in parts[:-1]:
            cur = cur.setdefault(p, {})
        cur[parts[-1]] = v
    for k in drop:
        cur, parts = doc, k.split(".")
        for p in parts[:-1]:
            cur = cur[p]
        del cur[parts[-1]]
    return yaml.safe_dump(doc, sort_keys=False)


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
    res = RS.resolve(recipes, name)
    for st, cfg in stages:
        ex = res["stages"][st]["explicit"]
        r = RS.check_config(recipes, cfg, name, st, bcond=None)
        fk = {k for k, _ in r["fails"]}
        EXCEPTIONS += len(ex)
        check(f"(b) {name}/{st}: 生成器の config の FAIL は #11 の例外表のキー (明示すべきキー {len(ex)} 個) の省略だけ",
              fk == set(ex) and all("明示すべきキー" in m for _, m in r["fails"]) and not r["warns"],
              f"FAIL のうち例外表に無い {sorted(fk - set(ex))[:3]}、例外表で FAIL しない {sorted(set(ex) - fk)[:3]}、warns {r['warns'][:3]}")
        r = RS.check_config(recipes, with_keys(cfg, add=ex), name, st, bcond=None)
        check(f"(b) {name}/{st}: 例外表のキーを旧い実効値で足した config (#11 の形) がプロファイルに PASS",
              not r["fails"] and not r["warns"], f"fails {r['fails'][:3]} warns {r['warns'][:3]}")
    main = [c for st, c in stages if st == "main"][-1]
    if any(st == "soft" for st, _ in stages):
        r = RS.check_config(recipes, main, name, "soft")
        check(f"(c) {name}: 本段の config を soft として宣言すると FAIL", bool(r["fails"]),
              sorted({k for k, _ in r["fails"]})[:5])
    other = "bell_tp_ns_sst" if name != "bell_tp_ns_sst" else "axismach_ns_sst_tp"
    r = RS.check_config(recipes, main, other, "main")
    check(f"(c) {name}: 別の機種のプロファイル ({other}) で宣言すると FAIL", bool(r["fails"]), sorted({k for k, _ in r["fails"]})[:5])

# (d) #16 の判別 A/B: SERN 2D 本段の生成器の config (例外表のキーを足した形)
SERN = "sern2d_finite_base_frozen_tp_sst"
main = [c for st, c in configs[SERN] if st == "main"][-1]
ex = RS.resolve(recipes, SERN)["stages"]["main"]["explicit"]
A = with_keys(main, add=ex)
B = with_keys(main, add=ex, drop=["space.limiterScaled"])
rA = RS.check_config(recipes, A, SERN, "main")
rB = RS.check_config(recipes, B, SERN, "main")
check("(d) A/B の前提: 生成器の SERN 2D 本段は space.limiterScaled: 1 を書く (fixed)",
      "limiterScaled: 1" in main and RS.resolve(recipes, SERN)["stages"]["main"]["fixed"].get("space.limiterScaled") == 1)
check("(d) A (space.limiterScaled: 1 を明記): PASS", not rA["fails"] and not rA["warns"], f"fails {rA['fails'][:3]}")
check("(d) B (space.limiterScaled だけ削除): 同じキーを名指しで FAIL",
      {k for k, _ in rB["fails"]} == {"space.limiterScaled"} and not rB["warns"], f"fails {[k for k, _ in rB['fails']]}")

# (e) #17 代表の問題 YAML の固定の写し: 書き写しは写しを読み、元 (case) との違いは別の検査で WARN、相対パスは元の所在で解決
import contextlib  # noqa: E402
import io  # noqa: E402
import shutil  # noqa: E402
import tempfile  # noqa: E402

import yaml  # noqa: E402

man = RG.verify_fixtures()
check("(e) 全機種の写しが MANIFEST の記録どおり (元のパス = SPECS の problem、写しの sha256)",
      set(man) == {sp["name"] for sp in RG.SPECS} and all(man[sp["name"]]["source"] == sp["problem"] for sp in RG.SPECS))
st = RG.fixture_status()
print("     元と写しの比較 (情報、WARN は失敗にしない): " + ", ".join(f"{n}={s}" for n, s, _ in st))
BELL = next(sp for sp in RG.SPECS if sp["name"] == "bell_tp_ns_sst")
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    fx = td / "fixtures"
    shutil.copytree(RG.FIXTURE_DIR, fx)
    # 相対パスは元の問題 YAML の所在で解決する (写しの所在・作業ディレクトリではない)
    p = RG.materialize(BELL, td / "w1", fixture_dir=fx)
    raw = yaml.safe_load(p.read_text())
    base = (RG.ROOT / BELL["problem"]).resolve().parent
    got = {k: raw[s][k] for s, k in (("gas", "exhaust"), ("geometry", "rao_wall_file"), ("geometry", "rao_offset_file"))}
    want = {k: str((base / v).resolve()) for k, v in
            {"exhaust": "exhaust_lox_lh2_pc100_of6.yaml",
             "rao_wall_file": "rao_m/rao_corrprod_parts/L100_corr_vext_n321_wall.npz",
             "rao_offset_file": "rao_m/rao_corrprod_parts/L100_corr_vext_n321_offset.npz"}.items()}
    check("(e) 写しの相対パス (gas.exhaust・rao_wall_file・rao_offset_file) が元の問題 YAML の所在の絶対パスになる",
          got == want and all(Path(v).is_file() for v in got.values()), f"{got}")
    from forge_design.probdef import load_problem  # noqa: E402
    pp = load_problem(p)
    check("(e) 作業ディレクトリに置いた写しを load_problem が読める (gas.exhaust の解決を含む)",
          pp.name == yaml.safe_load((fx / "bell_tp_ns_sst.yaml").read_text())["name"])
    # 書き写しは写しの中身を使う (元ではない): 写しの name だけ変えると置いた問題の name が変わる
    fake = yaml.safe_load((fx / "bell_tp_ns_sst.yaml").read_text())
    fake["name"] = "regen_fixture_probe"
    (fx / "bell_tp_ns_sst.yaml").write_text(yaml.safe_dump(fake, allow_unicode=True, sort_keys=False))
    raw2 = yaml.safe_load(RG.materialize(BELL, td / "w2", fixture_dir=fx).read_text())
    check("(e) materialize は写しの中身を使う (元の case の YAML を読まない)", raw2["name"] == "regen_fixture_probe")
    # 写しを手で書き換えると (MANIFEST の sha256 と違う) 書き写しの前に止まる
    try:
        RG.verify_fixtures(fixture_dir=fx)
        check("(e) 写しの手での書き換えは止まる", False, "例外にならない")
    except RuntimeError as e:
        check("(e) 写しの手での書き換えは機種を名指しで止まる", "bell_tp_ns_sst" in str(e) and "sha256" in str(e), str(e)[:200])
    # 既知のキーの外に、元の所在に実在するファイルの相対パスがあれば止まる (黙って作業ディレクトリ基準で解決しない)
    fake = yaml.safe_load((RG.FIXTURE_DIR / "bell_tp_ns_sst.yaml").read_text())
    fake["evaluate"]["some_table_file"] = "rao_m/rao_corrprod_parts/L100_corr_vext_n321_wall.npz"
    (fx / "bell_tp_ns_sst.yaml").write_text(yaml.safe_dump(fake, allow_unicode=True, sort_keys=False))
    try:
        RG.materialize(BELL, td / "w3", fixture_dir=fx)
        check("(e) 未知のキーの相対パスは止まる", False, "例外にならない")
    except RuntimeError as e:
        check("(e) 未知のキーの相対パス (元の所在に実在) は名指しで止まる", "evaluate.some_table_file" in str(e), str(e)[:200])
    # 写しの更新: 元のバイト列を写し、MANIFEST を書き直す (更新の後は記録どおり)
    done = RG.refresh_fixtures(["bell_tp_ns_sst"], fixture_dir=fx, today="2026-10-15")
    check("(e) --refresh-fixture: 元のバイト列の写しと MANIFEST の sha256 を書き直す",
          (fx / "bell_tp_ns_sst.yaml").read_bytes() == (RG.ROOT / BELL["problem"]).read_bytes()
          and done[0][2] == man["bell_tp_ns_sst"]["sha256"] and RG.verify_fixtures(fixture_dir=fx) is not None)
    # 元と違う写しは別の検査 (--check-fixtures) で WARN、終了 0。元が無いのも WARN
    root = td / "root"
    src = root / BELL["problem"]
    src.parent.mkdir(parents=True)
    src.write_bytes((RG.ROOT / BELL["problem"]).read_bytes() + "\n# case 側の編集\n".encode())
    saved = (RG.ROOT, RG.FIXTURE_DIR)
    RG.ROOT, RG.FIXTURE_DIR = root, fx
    try:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = RG.main(["--check-fixtures"])
        out = buf.getvalue()
        st2 = dict((n, s) for n, s, _ in RG.fixture_status())
    finally:
        RG.ROOT, RG.FIXTURE_DIR = saved
    check("(e) 元と違う写しは --check-fixtures で WARN・終了 0 (元の無いものも WARN)",
          rc == 0 and st2["bell_tp_ns_sst"] == "differs" and "WARN bell_tp_ns_sst" in out
          and all(s == "missing" for n, s in st2.items() if n != "bell_tp_ns_sst"), f"rc {rc} {out[:300]}")

print(f"\n#11 の例外表 (明示すべきキー) の延べ数: {EXCEPTIONS} (機種 {len(RG.SPECS)}・全段)")
print(f"\n{'ALL PASS' if FAIL == 0 else f'FAIL {FAIL}'}")
sys.exit(1 if FAIL else 0)
