#!/usr/bin/env python3
"""共通 species データ化 (plans/active/thermophysics-solver-owned-species-db.md §5.1 #4) の値のビット一致試験。

基準 (`tests/unit/data/species_builtin_baseline_v1.json`) と、現在のソースの値を 16 進浮動小数で比べる。
基準は段 3 (#13-3, 2026-10-01) で v0 (移行前に固定した値) から v1 (既存種を CEA thermo.inp そのものにした値) へ明示的に張り替えた。
v0 は削除せず、(T) で「v0 → 現在」の差が段 3 の想定 (Δ 表) どおりであることを毎回確かめる。

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

- #13-3 (既存種の CEA 化, 2026-10-01) の張り替え検査 (T): v0 → 現在の差分は次だけ。
  C++ 内蔵 (v0 の 13 キー): N2・O2・CO2 は MW・先頭 2 区間の区切りと全係数がビット一致 (第 3 区間 6000–20000 K が増えるだけ)、
  H2O (h2o・WATER) は MW だけ、He (HE) は MW と第 3 区間、Ar (AR) は第 2 区間の係数と第 3 区間、AIR (Air・air) は不変。
  変わった種の 200 ≤ T < 6000 K の max |Δcp|/cp・|Δh|・|Δs°| が段 3 の Δ 表 (plan #13-3 の事前確認) と相対 1e-3 で一致。
  Python (SPECIES_NASA9・builtin_db) は H2O の MW と AR の high だけが変わり、LJ_PARAMS・T_MID・BUILTIN_ATOMS は不変。
  互換性ハッシュは AIR だけの組 (air|*) が不変、他は全部変わる (旧→新を表示)。
- #13-3 の設計側 (plan #13-3 の追加合格条件 (1)–(3)):
  (D1) 設計側 5 モジュール (gas.semiperfect・gas.composition・gas.frozen・probdef・evaluate.ic) が段 3 のデータで import でき、
       `_load_design_species` (import 時に走る読み込み) は 1 区間の種・非標準の区切りの種を ValueError で拒否する (負例 2 件)。
  (D2) 内蔵種の cp/h/s° は 6000.0001 K で例外、5999.9999 K で値を返す (lump も同じ)。T < 200 K は従来どおり
       (cp は 200 K の値、h は 200 K から線形、s° は対数)。
  (D3) 設計側 11 種と `--eval` のソルバ値 (thermo_cp_mass / thermo_h_mass / thermo_s0_mass) が V2 格子
       (200–6000 K の 1000 点 + 999.99994/1000/1000.00006/5999.9999 K; 1000 点の端点 6000 K は除く) で相対 ≤4e-16
       (V2 と同じく cp・s° は点ごと、h は格子上の max|h| に対する相対)。6000 K ちょうどの差は info
       (3 区間の種はソルバが第 3 区間を選ぶので #13-0 (6) の境界段差になる)。

使い方:
  python3 solver_density_cuda/tests/unit/test_species_data_bitexact.py            # 比較 (ALL PASS / FAIL)
  python3 solver_density_cuda/tests/unit/test_species_data_bitexact.py --write-baseline   # 基準 v1 の作成 (既存は上書きしない)
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
BASELINE = os.path.join(HERE, "data", "species_builtin_baseline_v1.json")      # 段 3 (#13-3) 以降の基準
BASELINE_V0 = os.path.join(HERE, "data", "species_builtin_baseline_v0.json")   # 移行前の値 (張り替え検査 (T) だけが読む)
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


DUMP_EXE = None


def dump_cpp(workdir):
    """dump_species_builtin を現在のソースでビルドして JSON を返す (C の %a を float.hex に揃える)。"""
    global DUMP_EXE
    exe = build_dump(workdir)
    DUMP_EXE = exe
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


# 段 3 の Δ 表 (plan #13-3 の事前確認; 200 ≤ T < 6000 K, 1 K 刻み + 区切り直前直後, 質量あたり):
#   種: (max |Δcp|/cp, max |Δh| [J/kg], max |Δs°| [J/(kg K)])
DELTA_TABLE = {"N2": (0.0, 0.0, 0.0), "O2": (0.0, 0.0, 0.0), "CO2": (0.0, 0.0, 0.0),
               "H2O": (1.110e-06, 1.510e+01, 2.024e-02), "Ar": (4.870e-07, 7.314e-02, 1.266e-05),
               "He": (4.997e-07, 1.480e+01, 2.354e-02)}
V0_ALIAS = {"AR": "Ar", "HE": "He", "WATER": "H2O", "h2o": "H2O", "Air": "AIR", "air": "AIR"}


def _ivs(e):
    """dump の 1 エントリ (hex) → (MW, [(Tlo, Thi, [a0..a8]), ...]) (2 区間の書式と区間可変の書式の両方)。"""
    f = float.fromhex
    if "Tbounds" in e:
        b = [f(x) for x in e["Tbounds"]]
        return f(e["MW"]), [(b[k], b[k + 1], [f(x) for x in e[f"coef{k}"]]) for k in range(len(b) - 1)]
    return f(e["MW"]), [(f(e["Tlo"]), f(e["Tmid"]), [f(x) for x in e["low"]]), (f(e["Tmid"]), f(e["Thi"]), [f(x) for x in e["high"]])]


def check_v0_transition(v0, cur):
    """(T) 基準 v0 → 現在の差分が段 3 の想定どおり (docstring)。"""
    sys.path.insert(0, os.path.join(SOLVER, "tools"))
    from cea_thermo_to_forge_species import delta_props
    b0, bc = v0["cpp"]["builtin"], cur["cpp"]["builtin"]
    expect = {"N2": set(), "O2": set(), "CO2": set(), "H2O": {"MW"}, "He": {"MW"}, "Ar": {"iv1"}, "AIR": set()}
    rows = []
    for k in b0:
        sid = V0_ALIAS.get(k, k)
        if k not in bc:
            check(False, f"(T) cpp.builtin.{k}: 現在の内蔵に無い")
            continue
        (m0, iv0), (mc, ivc) = _ivs(b0[k]), _ivs(bc[k])
        diff = set()
        if m0.hex() != mc.hex():
            diff.add("MW")
        for n in range(2):
            if [iv0[n][0], iv0[n][1]] != [ivc[n][0], ivc[n][1]] or [x.hex() for x in iv0[n][2]] != [x.hex() for x in ivc[n][2]]:
                diff.add(f"iv{n}")
        extra = [(lo, hi) for lo, hi, _ in ivc[2:]]
        lj = (b0[k]["LJ_sigma"], b0[k]["LJ_eps_kB"]) == (bc[k]["LJ_sigma"], bc[k]["LJ_eps_kB"])
        want_extra = [] if sid == "AIR" or sid == "H2O" else [(6000.0, 20000.0)]
        ok = diff == expect[sid] and extra == want_extra and lj and len(iv0) == 2
        check(ok, f"(T) cpp.builtin.{k} ({sid}): 先頭 2 区間・MW の変化 {sorted(diff) or 'なし'} (想定 {sorted(expect[sid]) or 'なし'})、"
                  f"追加区間 {extra} (想定 {want_extra})、LJ {'不変' if lj else '変化'}")
        if k == sid and sid in DELTA_TABLE:
            got = delta_props(iv0, m0, ivc, mc)
            g = (got["cp"][0], got["h"][0], got["s"][0])
            exp = DELTA_TABLE[sid]
            ok = all((e == 0.0 and x == 0.0) or (e != 0.0 and abs(x - e) <= 1e-3 * e) for x, e in zip(g, exp))
            rows.append(sid)
            check(ok, f"(T) {sid}: 200 ≤ T < 6000 K の max |Δcp|/cp {g[0]:.3e}・|Δh| {g[1]:.3e} J/kg・|Δs°| {g[2]:.3e} J/(kg K)"
                      f" = Δ 表 ({exp[0]:.3e}, {exp[1]:.3e}, {exp[2]:.3e}; 相対 1e-3)")
    check(sorted(rows) == sorted(DELTA_TABLE), f"(T) Δ 表の全種 {sorted(DELTA_TABLE)} を照合した")
    # Python: H2O の MW と AR の high だけが変わる
    allowed = {("H2O", "MW"), ("AR", "high")}
    for tag in ("SPECIES_NASA9", "builtin_db"):
        p0, pc = dict((x[0], x[1]) for x in v0["py"][tag]), dict((x[0], x[1]) for x in cur["py"][tag])
        check(list(p0) == list(pc), f"(T) py.{tag}: キーと順序が v0 と同じ")
        d = {(k, f) for k in p0 for f in p0[k] if p0[k][f] != pc.get(k, {}).get(f)}
        check(d == allowed, f"(T) py.{tag}: v0 からの変化 {sorted(d)} = 想定 {sorted(allowed)}")
    for tag in ("LJ_PARAMS", "T_MID", "BUILTIN_ATOMS"):
        check(v0["py"][tag] == cur["py"][tag], f"(T) py.{tag}: v0 と同じ")
    h0, hc = v0["cpp"]["compat_hash"], cur["cpp"]["compat_hash"]
    for k in h0:
        same = h0[k] == hc.get(k)
        want_same = k.startswith("air|")
        check(same == want_same, f"(T) compat_hash[{k}]: {h0[k][:16]} → {str(hc.get(k))[:16]} ({'不変' if same else '変化'}; "
                                 f"想定 {'不変' if want_same else '変化'})")


def check_design_side(workdir):
    """(D1)–(D3) (docstring)。"""
    import importlib
    import math
    import yaml
    sys.path.insert(0, os.path.join(REPO, "design"))
    mods = ["forge_design.gas.semiperfect", "forge_design.gas.composition", "forge_design.gas.frozen",
            "forge_design.probdef", "forge_design.evaluate.ic"]
    bad = []
    for m in mods:
        try:
            importlib.import_module(m)
        except Exception as ex:   # noqa: BLE001 — 失敗の中身を出す
            bad.append(f"{m}: {ex}")
    check(not bad, f"(D1) 設計側 {len(mods)} モジュールが段 3 のデータで import できる" + (f" — {bad}" if bad else ""))
    import forge_design.gas.semiperfect as sp
    import forge_design.gas.composition as comp
    raw = yaml.safe_load(open(DATA_FILE, encoding="utf-8"))
    n2 = next(e for e in raw["species"] if str(e["id"]) == "N2")
    check(len(n2["intervals"]) == 3 and sp.SPECIES_NASA9["N2"]["low"] == [float(x) for x in n2["intervals"][0]["coeffs"]]
          and sp.SPECIES_NASA9["N2"]["high"] == [float(x) for x in n2["intervals"][1]["coeffs"]],
          "(D1) N2 (共通データ 3 区間) の設計側 low/high = 先頭 2 区間")
    for tag, edit in (("1 区間の種", lambda e: e.__setitem__("intervals", e["intervals"][:1])),
                      ("非標準の区切り (1000→1500 K)", lambda e: (e["intervals"][0].__setitem__("Thi", 1500.0),
                                                                e["intervals"][1].__setitem__("Tlo", 1500.0)))):
        r2 = yaml.safe_load(open(DATA_FILE, encoding="utf-8"))
        edit(next(e for e in r2["species"] if str(e["id"]) == "N2"))
        path = os.path.join(workdir, "design_neg.yaml")
        with open(path, "w", encoding="utf-8") as f:
            yaml.safe_dump(r2, f, allow_unicode=True)
        try:
            sp._load_design_species(path)
            msg, ok = "通った", False
        except ValueError as ex:
            msg, ok = str(ex).split(": ", 1)[-1][:120], True
        check(ok, f"(D1) 負例 {tag}: 設計側の読み込みが ValueError ({msg})")
    # (D2) 6000 K 超のガード・200 K 未満の扱い
    db = comp.ResolvedSpeciesDB.builtin()
    e = db["N2"]
    raised = []
    for fn in (db.species_cp_R, db.species_h_RT, db.species_s0_R):
        try:
            fn(e, 6000.0001)
        except ValueError:
            raised.append(fn.__name__)
    vals = [float(fn(e, 5999.9999)[0]) for fn in (db.species_cp_R, db.species_h_RT, db.species_s0_R)]
    check(len(raised) == 3 and all(math.isfinite(v) for v in vals),
          f"(D2) 内蔵種 N2 の cp/h/s° は 6000.0001 K で例外 {raised}、5999.9999 K で値 {[f'{v:.6g}' for v in vals]}")
    L = comp.lump_entry("MIXT", {"N2": 0.7, "O2": 0.3}, db)
    try:
        db.species_cp_R(L, 6000.0001)
        ok = False
    except ValueError:
        ok = True
    check(ok and L.T_eval_max == sp.DESIGN_T_MAX, f"(D2) 内蔵種の lump も 6000.0001 K で例外 (T_eval_max {L.T_eval_max})")
    try:
        sp.check_design_T([300.0, 6000.0001], "IC")
        ok = False
    except ValueError:
        ok = True
    sp.check_design_T([200.0, 6000.0], "IC")
    check(ok, "(D2) evaluate/ic.py が使う check_design_T は 6000.0001 K で例外、6000 K ちょうどは通す")
    lo = [float(x) for x in n2["intervals"][0]["coeffs"]]
    T = 150.0
    cp200, h200 = float(sp._cp_R_raw(lo, 200.0)), float(sp._h_RT_raw(lo, 200.0)) * 200.0
    s200 = None
    from forge_design.gas.frozen import _s0_R_raw
    s200 = float(_s0_R_raw(__import__("numpy").asarray(lo), 200.0))
    got = (float(db.species_cp_R(e, T)[0]), float(db.species_h_RT(e, T)[0]) * T, float(db.species_s0_R(e, T)[0]))
    exp = (cp200, h200 + cp200 * (T - 200.0), s200 + cp200 * math.log(T / 200.0))
    rel = max(abs(g - x) / abs(x) for g, x in zip(got, exp))
    check(rel <= 4e-16, f"(D2) T < 200 K (150 K) は従来どおり cp 固定・h 線形・s° 対数 (相対 {rel:.1e})")
    # (D3) 設計側 vs ソルバ
    names = list(sp.SPECIES_NASA9)
    # V2 格子の 1000 点 (200–6000 K の等分) の端点 6000 K ちょうどは区間選択の規約で決まる点なので厳密比較から外し、下の info に回す
    grid = [200.0 + (6000.0 - 200.0) * k / 999 for k in range(999)] + [999.99994, 1000.0, 1000.00006, 5999.9999]
    out = subprocess.run([DUMP_EXE, "--eval", ",".join(names)] + [repr(t) for t in grid + [6000.0]],
                         check=True, capture_output=True, text=True).stdout
    sol = {k: [[float.fromhex(x) for x in row] for row in v] for k, v in json.loads(out).items()}
    RUv = sp.RU
    worst, info6000 = (0.0, None), []
    for k in names:
        ek = db[k]
        Ta = __import__("numpy").asarray(grid + [6000.0])
        cp = db.species_cp_R(ek, Ta[:-1]) * RUv / ek.MW
        h = db.species_h_RT(ek, Ta[:-1]) * RUv * Ta[:-1] / ek.MW
        s0 = db.species_s0_R(ek, Ta[:-1]) * RUv / ek.MW
        # 相対の取り方は V2 (test_species_lump_solver.py) と同じ: cp・s° は点ごと、h は格子上の max|h| で割る
        # (絶対エンタルピーは 0 を横切るので点ごとの相対は丸めの桁を表さない)
        hscale = max(max(abs(x) for x in h), max(abs(r_[1]) for r_ in sol[k][:len(grid)]))
        for i in range(len(grid)):
            for q, (d, c) in enumerate(zip((cp[i], h[i], s0[i]), sol[k][i])):
                den = hscale if q == 1 else max(abs(d), abs(c))
                r = 0.0 if d == c else abs(d - c) / den
                if r > worst[0]:
                    worst = (r, (k, ("cp", "h", "s0")[q], grid[i]))
        # 6000 K ちょうど (設計側は第 2 区間の上端で評価できる; ソルバは 3 区間の種で第 3 区間)
        dh = float(db.species_h_RT(ek, 6000.0)[0]) * RUv * 6000.0 / ek.MW - sol[k][-1][1]
        info6000.append(f"{k} {dh:+.3e}")
    check(worst[0] <= 4e-16, f"(D3) 設計側 {len(names)} 種 vs ソルバ (V2 格子 {len(grid)} 点) の cp/h/s° 相対差の最大 {worst[0]:.2e} "
                             f"@ {worst[1]} (許容 4e-16)")
    print("[INFO] (D3) 6000 K ちょうどの h 差 (設計 − ソルバ) [J/kg]: " + ", ".join(info6000))
    # 差の出所の切り分け: 同じ係数を**ソルバと同じ演算順** (RU*(a0*Ti2 + a1*Ti + ... + a6*T*T*T*T))/MW で Python 評価すると
    # ソルバの cp と一致するか (一致すれば差は設計側 `_cp_R_raw` の演算順 (a0/T**2, a5*T**3 …) による丸めだけ)
    np_ = __import__("numpy")
    Ta = np_.asarray(grid)
    ro = 0.0
    for k in names:
        ek = db[k]
        a_ = np_.where((Ta < ek.Tmid)[:, None], np_.asarray(ek.low), np_.asarray(ek.high))
        Ti = 1.0 / Ta
        Ti2 = Ti * Ti
        cps = RUv * (a_[:, 0] * Ti2 + a_[:, 1] * Ti + a_[:, 2] + a_[:, 3] * Ta + a_[:, 4] * Ta * Ta + a_[:, 5] * Ta * Ta * Ta
                     + a_[:, 6] * Ta * Ta * Ta * Ta) / ek.MW
        ref = np_.asarray([r_[0] for r_ in sol[k][:len(grid)]])
        ro = max(ro, float(np_.max(np_.abs(cps - ref) / np_.abs(ref))))
    print(f"[INFO] (D3) 同じ係数をソルバと同じ演算順で Python 評価した cp とソルバ cp の相対差の最大 {ro:.2e} "
          "(0 なら上の差は設計側の演算順による丸め; 段 3 のデータとは無関係)")


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
            check_design_side(td)
    if a.write_baseline:
        if os.path.exists(BASELINE):
            raise SystemExit(f"{BASELINE} exists; refusing to overwrite the baseline")
        os.makedirs(os.path.dirname(BASELINE), exist_ok=True)
        with open(BASELINE, "w") as f:
            json.dump(cur, f, indent=1)
            f.write("\n")
        print(f"wrote {BASELINE}")
        return
    with open(BASELINE) as f:
        base = json.load(f)
    with open(BASELINE_V0) as f:
        check_v0_transition(json.load(f), cur)
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
