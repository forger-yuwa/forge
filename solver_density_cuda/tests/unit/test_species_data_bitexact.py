#!/usr/bin/env python3
"""共通 species データ化 (plans/active/thermophysics-solver-owned-species-db.md §5.1 #4) の値のビット一致試験。

移行前に固定した基準 (`tests/unit/data/species_builtin_baseline_v0.json`) と、現在のソースの値を 16 進浮動小数で比べる。

- C++: `tests/unit/dump_species_builtin.cpp` を `input/speciesDB.cpp` とビルドして実行し、
  `speciesDB_builtin()` の全キー (別名込み)・名前解決 (config の綴りと大小文字違い)・内蔵種だけの互換性ハッシュを比べる。
  共通データ (`data/species/forge_species_v1.yaml`) の埋め込みヘッダは `cmake/embed_species_data.cmake` で生成する。
- Python: `design/forge_design/gas/semiperfect.py` の `SPECIES_NASA9`・`LJ_PARAMS`・`T_MID`、
  `composition.py` の `BUILTIN_ATOMS`・`ResolvedSpeciesDB.builtin()` の全エントリ (キーの順序・値)。

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


def dump_cpp(workdir):
    """dump_species_builtin を現在のソースでビルドして JSON を返す (C の %a を float.hex に揃える)。"""
    inc = []
    if os.path.exists(EMBED_SCRIPT):
        gen = os.path.join(workdir, "gen")
        subprocess.run(["cmake", f"-DIN={DATA_FILE}", f"-DOUT={os.path.join(gen, 'forge_species_data.hpp')}",
                        "-P", EMBED_SCRIPT], check=True)
        inc = ["-I", gen]
    exe = os.path.join(workdir, "dump_species_builtin")
    subprocess.run(["g++", "-O1", "-std=c++17", "-I", SOLVER, *inc,
                    os.path.join(HERE, "dump_species_builtin.cpp"), os.path.join(SOLVER, "input", "speciesDB.cpp"),
                    "-lyaml-cpp", "-o", exe], check=True)
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
    compare("cpp.builtin", base["cpp"]["builtin"], cur["cpp"]["builtin"])
    compare("cpp.resolve", base["cpp"]["resolve"], cur["cpp"]["resolve"])
    compare("cpp.compat_hash", base["cpp"]["compat_hash"], cur["cpp"]["compat_hash"])
    for k in ("SPECIES_NASA9", "LJ_PARAMS", "T_MID", "BUILTIN_ATOMS", "builtin_db"):
        compare(f"py.{k}", base["py"][k], cur["py"][k])
    print("ALL PASS" if FAIL == 0 else f"FAIL ({FAIL})")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
