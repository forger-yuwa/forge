#!/usr/bin/env python3
"""共通 species データ化 (plans/active/thermophysics-solver-owned-species-db.md §5.1 #4) の値のビット一致試験。

移行前に固定した基準 (`tests/unit/data/species_builtin_baseline_v0.json`) と、現在のソースの値を 16 進浮動小数で比べる。

- C++: `tests/unit/dump_species_builtin.cpp` を `input/speciesDB.cpp` とビルドして実行し、
  `speciesDB_builtin()` の全キー (別名込み)・名前解決 (config の綴りと大小文字違い)・内蔵種だけの互換性ハッシュを比べる。
  共通データ (`data/species/forge_species_v1.yaml`) の埋め込みヘッダは `cmake/embed_species_data.cmake` で生成する。
- Python: `design/forge_design/gas/semiperfect.py` の `SPECIES_NASA9`・`LJ_PARAMS`・`T_MID`、
  `composition.py` の `BUILTIN_ATOMS`・`ResolvedSpeciesDB.builtin()` の全エントリ (キーの順序・値)。
- #13-2 (内蔵拡大, 2026-10-01) の拡張:
  (E1) C++ 内蔵のキーは基準の全キー (値はビット一致) + 共通データの気相 id のうち基準に無いもの**だけ**。増えたキーの値
       (MW・区切り・全係数・LJ; LJ: null は 0) が共通データ (YAML を Python で読んだ double) とビット一致。
  (E2) 名前解決: 基準で見つかった綴りは不変、基準で "not found" の綴りが見つかるようになったら、それは大小文字無視で
       新しい内蔵種 1 つに当たり、その値が共通データと一致すること (例 H2, CO, Co → CO; 大小文字無視の解決は plan #8 まで現行どおり)。
  (E3) 共通データの名前が大小文字無視で一意でないと起動を拒否する: 共通データに CO の写しを id `Co` で足したヘッダで
       dump を作り、非ゼロ終了とメッセージ ("differs only in letter case") を確認する (負例 1 件)。
  (E4) LJ: null の内蔵種は kinetic 輸送と LJ の混合平均拡散に使うと拒否、LJ のある種・拡散を使わない設定は通る (dump の lj_use)。

使い方:
  python3 solver_density_cuda/tests/unit/test_species_data_bitexact.py            # 比較 (ALL PASS / FAIL)
  python3 solver_density_cuda/tests/unit/test_species_data_bitexact.py --write-baseline   # 基準の作成 (移行前に 1 回だけ; 既存は上書きしない)
規約: [PASS]/[FAIL] を出し、失敗があれば非ゼロ終了。
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SOLVER = os.path.normpath(os.path.join(HERE, "..", ".."))
REPO = os.path.dirname(SOLVER)
BASELINE = os.path.join(HERE, "data", "species_builtin_baseline_v0.json")
EMBED_SCRIPT = os.path.join(SOLVER, "cmake", "embed_species_data.cmake")
DATA_FILE = os.path.join(SOLVER, "data", "species", "forge_species_v1.yaml")

FAIL = 0


def check(ok, what):
    global FAIL
    print(("[PASS] " if ok else "[FAIL] ") + what)
    if not ok:
        FAIL += 1


def hx(v):
    return float(v).hex()


def build_dump(workdir, data_file=DATA_FILE, tag="gen"):
    """dump_species_builtin を現在のソースと data_file の埋め込みヘッダでビルドして実行ファイルのパスを返す。"""
    inc = []
    if os.path.exists(EMBED_SCRIPT):
        gen = os.path.join(workdir, tag)
        subprocess.run(["cmake", f"-DIN={data_file}", f"-DTRANS={os.path.join(os.path.dirname(DATA_FILE), 'forge_transport_v1.yaml')}",
                        f"-DOUT={os.path.join(gen, 'forge_species_data.hpp')}", "-P", EMBED_SCRIPT], check=True)
        inc = ["-I", gen]
    exe = os.path.join(workdir, "dump_species_builtin_" + tag)
    subprocess.run(["g++", "-O1", "-std=c++17", "-I", SOLVER, *inc,
                    os.path.join(HERE, "dump_species_builtin.cpp"), os.path.join(SOLVER, "input", "speciesDB.cpp"),
                    os.path.join(SOLVER, "input", "speciesTransportDB.cpp"),   # speciesDB_resolve(cfg) が輸送の解決を呼ぶ (#5t2)
                    os.path.join(SOLVER, "input", "solverConfig.cpp"),         # (E4) の solverConfig (#13-2)
                    "-lyaml-cpp", "-o", exe], check=True)
    return exe


def dump_cpp(workdir):
    """dump_species_builtin を現在のソースでビルドして JSON を返す (C の %a を float.hex に揃える)。"""
    exe = build_dump(workdir)
    out = subprocess.run([exe], check=True, capture_output=True, text=True).stdout
    d = json.loads(out)

    def norm(e):
        if isinstance(e, str):
            return e
        return {k: ([hx(float.fromhex(x)) for x in v] if isinstance(v, list) else hx(float.fromhex(v))) for k, v in e.items()}

    d["builtin"] = {k: norm(v) for k, v in d["builtin"].items()}
    d["resolve"] = {k: (v if isinstance(v, str) else {"source": v["source"], "entry": norm(v["entry"])})
                    for k, v in d["resolve"].items()}
    return d


def dump_py():
    sys.path.insert(0, os.path.join(REPO, "design"))
    import forge_design.gas.semiperfect as sp
    import forge_design.gas.composition as comp
    out = {
        "SPECIES_NASA9": [[k, {kk: ([hx(x) for x in vv] if isinstance(vv, list) else hx(vv)) for kk, vv in e.items()}]
                          for k, e in sp.SPECIES_NASA9.items()],
        "LJ_PARAMS": [[k, [hx(x) for x in v]] for k, v in sp.LJ_PARAMS.items()],
        "T_MID": hx(sp.T_MID),
        "BUILTIN_ATOMS": [[k, [[a, type(n).__name__, repr(n)] for a, n in v.items()]] for k, v in comp.BUILTIN_ATOMS.items()],
        "builtin_db": [],
    }
    for k, e in comp.ResolvedSpeciesDB.builtin().entries.items():
        out["builtin_db"].append([k, {
            "name": e.name, "MW": hx(e.MW), "low": [hx(x) for x in e.low], "high": [hx(x) for x in e.high],
            "T": [hx(e.Tlo), hx(e.Tmid), hx(e.Thi)], "LJ": [hx(e.LJ_sigma), hx(e.LJ_eps_kB)],
            "atoms": [[a, type(n).__name__, repr(n)] for a, n in e.atoms.items()], "source": e.source}])
    return out


def data_gas_entries():
    """共通データの気相エントリ (Python double) を、dump の JSON と同じ正規化 (float.hex) で {id: エントリ} と {別名: id} にする。"""
    import yaml
    raw = yaml.safe_load(open(DATA_FILE, encoding="utf-8"))
    ents, alias = {}, {}
    for e in raw["species"]:
        if e.get("phase") != "gas":
            continue
        iv = e["intervals"]
        lj = e.get("LJ") or {}
        sig, eps = float(lj.get("sigma", 0.0)), float(lj.get("eps_kB", 0.0))
        bounds = [float(iv[0]["Tlo"])] + [float(x["Thi"]) for x in iv]
        if len(iv) == 2:
            d = {"MW": hx(e["MW"]), "Tlo": hx(bounds[0]), "Tmid": hx(bounds[1]), "Thi": hx(bounds[2]),
                 "LJ_sigma": hx(sig), "LJ_eps_kB": hx(eps),
                 "low": [hx(x) for x in iv[0]["coeffs"]], "high": [hx(x) for x in iv[1]["coeffs"]]}
        else:
            d = {"MW": hx(e["MW"]), "Tbounds": [hx(x) for x in bounds], "LJ_sigma": hx(sig), "LJ_eps_kB": hx(eps)}
            for k, x in enumerate(iv):
                d[f"coef{k}"] = [hx(v) for v in x["coeffs"]]
        ents[str(e["id"])] = d
        for a in e.get("aliases") or []:
            alias[str(a)] = str(e["id"])
    return ents, alias


def compare_extended(base, cur):
    """(E1)(E2): 基準のキーは値ごと不変、増えたキーは共通データの新しい気相種 (とその別名) だけで値が一致。"""
    ents, alias = data_gas_entries()
    bb, cb = base["cpp"]["builtin"], cur["cpp"]["builtin"]
    missing = [k for k in bb if k not in cb]
    check(not missing, f"cpp.builtin: 基準の {len(bb)} キーが残っている" + (f" — 欠落 {missing}" if missing else ""))
    changed = [k for k in bb if k in cb and bb[k] != cb[k]]
    check(not changed, f"cpp.builtin: 基準の {len(bb)} キーの値がビット一致" + (f" — 変化 {changed}" if changed else ""))
    base_ids = {alias.get(k, k) for k in bb}
    expect_new = {k for k in ents if k not in base_ids} | {a for a, i in alias.items() if i not in base_ids}
    added = set(cb) - set(bb)
    check(added == expect_new, f"cpp.builtin: 増えたキー {len(added)} 個 = 共通データの新しい気相種と別名 {len(expect_new)} 個"
          + ("" if added == expect_new else f" — 余分 {sorted(added - expect_new)} / 不足 {sorted(expect_new - added)}"))
    bad = [k for k in added & expect_new if cb[k] != ents[alias.get(k, k)]]
    n3 = sum(1 for k in added if "Tbounds" in cb[k])
    nlj0 = sum(1 for k in added if float.fromhex(cb[k]["LJ_sigma"]) == 0.0)
    check(not bad, f"cpp.builtin: 増えた {len(added)} 種 (3 区間 {n3}, LJ なし {nlj0}) の MW・区切り・全係数・LJ が共通データとビット一致"
          + (f" — 不一致 {sorted(bad)}" if bad else ""))
    br, cr = base["cpp"]["resolve"], cur["cpp"]["resolve"]
    check(list(br) == list(cr), "cpp.resolve: 綴りの並びが同じ")
    ids_ci = {}
    for i in ents:
        ids_ci.setdefault(i.upper(), []).append(i)
    for a, i in alias.items():
        ids_ci.setdefault(a.upper(), []).append(i)
    for k in br:
        if br[k] != "not found" or cr.get(k) == "not found":
            check(br[k] == cr.get(k), f"cpp.resolve.{k}" + ("" if br[k] == cr.get(k) else f": base {br[k]} vs current {cr.get(k)}"))
            continue
        hit = sorted(set(ids_ci.get(k.upper(), [])))
        ok = (len(hit) == 1 and hit[0] not in base_ids and isinstance(cr.get(k), dict)
              and cr[k]["source"] == "builtin" and cr[k]["entry"] == ents[hit[0]])
        check(ok, f"cpp.resolve.{k}: 基準 not found → 新しい内蔵種 {hit} (#13-2)" + ("" if ok else f": current {cr.get(k)}"))


def check_lj_use(cur):
    """(E4) LJ: null の内蔵種の使用可否。"""
    lu = cur["cpp"].get("lj_use", {})
    expect = {
        "kinetic Kr (LJ null)": "rejected: .*no Lennard-Jones data",
        "kinetic N (LJ from the Cantera table)": "ok",
        "kinetic O2": "ok",
        "diffusion N2+Kr visc1 diff1": "rejected: .*Kr have no Lennard-Jones data",
        "diffusion N2+Kr visc1 diff0": "ok",
        "diffusion N2+Kr visc0 diff1": "ok",
        "diffusion MIXKR(N2,Kr)+O2 visc1 diff1": "rejected: .*MIXKR have no Lennard-Jones data",
        "diffusion N2+N visc1 diff1": "ok",
    }
    import re
    for k, pat in expect.items():
        v = lu.get(k, "(missing)")
        ok = re.match(pat, v) is not None
        check(ok, f"lj_use.{k}: {v[:140]}")


def check_case_collision(workdir):
    """(E3) 共通データに CO の写しを id `Co` で足すと、名前の大小文字無視の衝突で起動を拒否する (負例 1 件)。"""
    import yaml
    text = open(DATA_FILE, encoding="utf-8").read()
    raw = yaml.safe_load(text)
    ids = [str(e["id"]) for e in raw["species"]]
    assert "CO" in ids and "Co" not in ids
    lines = text.splitlines(keepends=True)
    s = next(k for k, L in enumerate(lines) if L.strip() == '- id: "CO"')
    e = s + 1
    while e < len(lines) and not lines[e].lstrip().startswith("- id:") and not lines[e].lstrip().startswith("# ----"):
        e += 1
    block = "".join(lines[s:e]).replace('- id: "CO"', '- id: "Co"', 1)
    bad = os.path.join(workdir, "forge_species_v1_co_collision.yaml")
    with open(bad, "w", encoding="utf-8") as f:
        f.write(text.rstrip("\n") + "\n" + block)
    exe = build_dump(workdir, bad, "gen_co")
    r = subprocess.run([exe], capture_output=True, text=True)
    msg = (r.stderr or "") + (r.stdout or "")
    ok = r.returncode != 0 and "differs only in letter case" in msg and "'Co'" in msg and "'CO'" in msg
    tail = [ln for ln in msg.splitlines() if "letter case" in ln]
    check(ok, f"(E3) CO と Co が並ぶ共通データは起動を拒否 (rc {r.returncode}): {tail[0].strip()[:200] if tail else msg[-200:]}")


def compare(tag, base, cur):
    if isinstance(base, dict) and isinstance(cur, dict):
        ok = list(base) == list(cur)
        check(ok, f"{tag}: keys/order {'same' if ok else f'differ: base {list(base)} vs current {list(cur)}'}")
        for k in base:   # 値 (1 種分の辞書・ハッシュ文字列) は丸ごと一致を見る
            if k in cur:
                check(base[k] == cur[k], f"{tag}.{k}" + ("" if base[k] == cur[k] else f": base {base[k]} vs current {cur[k]}"))
        return
    if isinstance(base, list) and isinstance(cur, list) and base and isinstance(base[0], list) and len(base[0]) == 2 \
            and isinstance(base[0][0], str):
        bk, ck = [x[0] for x in base], [x[0] for x in cur]
        check(bk == ck, f"{tag}: keys/order " + ("same" if bk == ck else f"differ: base {bk} vs current {ck}"))
        cd = dict((x[0], x[1]) for x in cur)
        for k, v in base:
            check(cd.get(k) == v, f"{tag}.{k}" + ("" if cd.get(k) == v else f": base {v} vs current {cd.get(k)}"))
        return
    check(base == cur, tag + ("" if base == cur else f": base {base} vs current {cur}"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write-baseline", action="store_true", help="基準 JSON を書く (既存なら拒否)")
    a = ap.parse_args()
    with tempfile.TemporaryDirectory() as td:
        cur = {"cpp": dump_cpp(td), "py": dump_py()}
        if not a.write_baseline:
            check_case_collision(td)
    if a.write_baseline:
        if os.path.exists(BASELINE):
            raise SystemExit(f"{BASELINE} exists; refusing to overwrite the pre-migration baseline")
        os.makedirs(os.path.dirname(BASELINE), exist_ok=True)
        with open(BASELINE, "w") as f:
            json.dump(cur, f, indent=1)
            f.write("\n")
        print(f"wrote {BASELINE}")
        return
    with open(BASELINE) as f:
        base = json.load(f)
    # #13-2: 内蔵が共通データの全気相種に広がったので、キー集合の完全一致でなく (E1)(E2) で見る
    compare_extended(base, cur)
    check_lj_use(cur)
    compare("cpp.compat_hash", base["cpp"]["compat_hash"], cur["cpp"]["compat_hash"])
    for k in ("SPECIES_NASA9", "LJ_PARAMS", "T_MID", "BUILTIN_ATOMS", "builtin_db"):
        compare(f"py.{k}", base["py"][k], cur["py"][k])
    print("ALL PASS" if FAIL == 0 else f"FAIL ({FAIL})")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
