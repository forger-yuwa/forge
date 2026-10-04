#!/usr/bin/env python3
"""lump (擬似種) の起動時合成と config 記法の試験 (plans/active/thermophysics-solver-owned-species-db.md §5.1 #6a, §6 V2)。GPU 不要。

  python3 solver_density_cuda/tests/unit/test_species_lump_solver.py --forge BIN [--base-forge BASE_BIN] [--keep]

`forge --resolve-species` だけを使う (計算はしない)。case/44 の run は読むだけで書かない (config は一時ディレクトリへ複製して書き換える)。

  (V2)  case/44 va3 の lump (run_0510 prepare_info.json の species.X から H2O を除いて正規化した**全桁の**モル分率) を
        config の lump 記法で書き、ソルバが合成した MIXDRY (解決済み記録) が run_0510 の生成 species_db.yaml の MIXDRY と
        係数・MW・LJ で相対 1e-12 以内。cp/h/s° を 200–6000 K の 1000 点で相対 1e-12 かつ絶対 cp 1e-9 J/(kg K)・h 1e-6 J/kg・
        s° 1e-9 J/(kg K) 以内 (評価式は同じ numpy の NASA-9; h の相対は |h| の範囲最大で規格化 — 下記)。
  (V2i) 独立検算: 合成 MIXDRY の cp/h/s° (質量基準) = 構成種 (記録の構成種係数) の質量分率加重和 (同じ許容差)。
  (M)   basis: mass (run_0510 の Y から H2O を除いて正規化) で書いても、合成係数・MW・LJ が mole 版と相対 1e-12 以内。
  (R)   記録: lump の構成 (構成種名・basis・入力分率・正規化モル分率・構成種の係数) が出ている、Python load_record の再計算が
        記録の互換性ハッシュと一致、記録の x を 1 つ書き換えると compare_signatures が該当構成種を示す。
  (H)   文字列リストの config は互換性ハッシュ・記録とも不変: --base-forge (本変更前のバイナリ) と新バイナリで
        case/44 run_0509 (外部 DB; 4378b7d78339ba27) と内蔵のみ 3 構成のハッシュと記録ファイルがバイト一致
        (#14-L1 から記録の provenance に LJ の出所 lj_source・lj_resolved が増えたので、その 2 欄を除いた本文で比べる)。
  (L)   run_0509 の config を lump 記法 (speciesDBFile なし) にすると --resolve-species が通り、ハッシュは外部 DB 版と違う (構成情報が入るため)。
  (B6)  区切りの違う構成種 (外部 DB で Tmid=1500 の試験種) は区切りの和集合で合成される (#6b → #13-1; 以前は拒否): 区切り 200/1000/1500/6000、
        合成規約は和集合の文字列、記録の schema は forge_resolved_species_v1_nint (値の検査は test_thermo_intervals_host.cpp の V3)。
  (N)   拒否: 凝縮種を lump に入れる / lump 名の衝突 (内蔵・外部 DB・大小文字違い) /
        分率 0・負・非有限 / basis なし・未知 / 未知の構成種 / 構成種の重複 (別名 AR と Ar) / lump に未知キー。
        分率の総和が 1 から外れると警告して通る。
規約: [PASS]/[FAIL]、失敗があれば非ゼロ終了。

h の相対誤差について: 絶対エンタルピー (生成エンタルピー込み) は MIXDRY で 200–6000 K の間に符号が変わる (0 を通る) ので、
点ごとの |Δh|/|h| は 0 付近で意味を持たない。相対は範囲内の max|h| で規格化し、絶対 1e-6 J/kg も併せて見る。
"""
import argparse, json, math, os, re, shutil, subprocess, sys, tempfile

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.normpath(os.path.join(HERE, "..", "..", "tools"))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, TOOLS)
import forge_species as fs  # noqa: E402

RU = 8.314462618
CASE44 = os.path.join(REPO, "case", "44.vitiated_air_wt")
RUN0510 = os.path.join(CASE44, "run_0510_va3_M4.19_Lc8_noneq_lumpX")
RUN0509 = os.path.join(CASE44, "run_0509_va3_M4.19_Lc8_dry_lumpX")
FAIL = 0


def check(ok, what):
    global FAIL
    print(("[PASS] " if ok else "[FAIL] ") + what, flush=True)
    if not ok:
        FAIL += 1


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


# ---- NASA-9 (thermo_d.cuh と同じ式; 区間 k = Tb[k] <= T < Tb[k+1] (区切りちょうどは上)。評価範囲 200–6000 K は外挿域を含まない) ----
def ivs(e):
    """記録 / 外部 DB のエントリ → (境界 [Tlo, 区切り..., Thi], 係数の列)。2 区間の書式と区間可変の書式 (#13-1) の両方。"""
    if "Tbounds" in e:
        return [float(x) for x in e["Tbounds"]], [[float(c) for c in r] for r in e["nasa9_intervals"]]
    return [float(e["Tlo"]), float(e["Tmid"]), float(e["Thi"])], [[float(c) for c in e["nasa9_low"]], [float(c) for c in e["nasa9_high"]]]


def _coef(e, T):
    Tb, co = ivs(e)
    k = np.zeros(len(T), dtype=int)
    for j in range(1, len(co)):
        k = np.where(T >= Tb[j], j, k)
    return np.asarray(co)[k]


def cp_mass(e, T):
    a = _coef(e, T)
    return RU / e["MW"] * (a[:, 0] / T**2 + a[:, 1] / T + a[:, 2] + a[:, 3] * T + a[:, 4] * T**2 + a[:, 5] * T**3 + a[:, 6] * T**4)


def h_mass(e, T):
    a = _coef(e, T)
    hRT = (-a[:, 0] / T**2 + a[:, 1] * np.log(T) / T + a[:, 2] + a[:, 3] * T / 2 + a[:, 4] * T**2 / 3 + a[:, 5] * T**3 / 4
           + a[:, 6] * T**4 / 5 + a[:, 7] / T)
    return RU / e["MW"] * T * hRT


def s0_mass(e, T):
    a = _coef(e, T)
    sR = (-a[:, 0] / (2 * T**2) - a[:, 1] / T + a[:, 2] * np.log(T) + a[:, 3] * T + a[:, 4] * T**2 / 2 + a[:, 5] * T**3 / 3
          + a[:, 6] * T**4 / 4 + a[:, 8])
    return RU / e["MW"] * sR


def compare_props(A, B, tag, T):
    """A, B: callable(T) -> (cp, h, s) 質量基準。V2 の許容差で判定。"""
    cpA, hA, sA = A(T)
    cpB, hB, sB = B(T)
    rcp = np.max(np.abs(cpA - cpB) / np.maximum(np.abs(cpA), np.abs(cpB)))
    rs = np.max(np.abs(sA - sB) / np.maximum(np.abs(sA), np.abs(sB)))
    hscale = max(np.max(np.abs(hA)), np.max(np.abs(hB)))
    rh = np.max(np.abs(hA - hB)) / hscale
    acp, ah, as_ = np.max(np.abs(cpA - cpB)), np.max(np.abs(hA - hB)), np.max(np.abs(sA - sB))
    ok = rcp <= 1e-12 and rh <= 1e-12 and rs <= 1e-12 and acp <= 1e-9 and ah <= 1e-6 and as_ <= 1e-9
    check(ok, f"{tag}: 200-6000 K x{len(T)}  cp rel {rcp:.2e} abs {acp:.2e} J/kgK | h rel(max|h|) {rh:.2e} abs {ah:.2e} J/kg "
              f"| s0 rel {rs:.2e} abs {as_:.2e} J/kgK")


class Ctx:
    def __init__(self, a, root):
        self.a, self.root = a, root
        self.n = 0

    def make(self, seed_run, species_yaml, drop_db=True, db=None, extra_sub=None):
        """seed_run の solverConfig.yaml を複製し physProp.species を species_yaml (flow 形式の文字列) に置き換える。"""
        self.n += 1
        d = os.path.join(self.root, f"c{self.n:02d}")
        os.makedirs(d)
        t = open(os.path.join(seed_run, "solverConfig.yaml")).read()
        t, k = re.subn(r"species: *\[[^\]]*\]", "species: " + species_yaml, t, count=1)
        assert k == 1, "physProp.species not found"
        if drop_db:
            t = re.sub(r",? *speciesDBFile: *\"?[^,}\"]*\"?", "", t, count=1)
        if db is not None:
            t = t.replace("thermoHrefTemp:", "speciesDBFile: \"test_db.yaml\", thermoHrefTemp:", 1)
            with open(os.path.join(d, "test_db.yaml"), "w") as f:
                yaml.safe_dump(db, f, sort_keys=False)
        if extra_sub:
            for pat, rep in extra_sub:
                t = re.sub(pat, rep, t)
        with open(os.path.join(d, "solverConfig.yaml"), "w") as f:
            f.write(t)
        return d

    def copy_run_config(self, seed_run, with_db=True):
        self.n += 1
        d = os.path.join(self.root, f"c{self.n:02d}")
        os.makedirs(d)
        shutil.copy(os.path.join(seed_run, "solverConfig.yaml"), d)
        if with_db and os.path.exists(os.path.join(seed_run, "species_db.yaml")):
            shutil.copy(os.path.join(seed_run, "species_db.yaml"), d)
        return d

    def resolve(self, d, exe=None):
        p = subprocess.run([exe or self.a.forge, "--resolve-species"], cwd=d, capture_output=True, text=True)
        lines = p.stdout.strip().splitlines()
        h = lines[-1].strip() if lines else ""
        m = re.search(r"\[species\] record (\S+) \(sha256 ([0-9a-f]{64})\)", p.stderr)
        rec = fs.load_record(os.path.join(d, m.group(1))) if (p.returncode == 0 and m) else None
        return p.returncode, h, rec, p.stderr


def flow_lump(name, fr, basis):
    return "{name: %s, lump: {%s}, basis: %s}" % (name, ", ".join(f"{k}: {v!r}" for k, v in fr.items()), basis)


def strip_lj_provenance(text):
    """記録本文から provenance の lj_source / lj_resolved (#14-L1 で増えた来歴) を除く (旧バイナリの記録と比べるため)。"""
    out, skip = [], False
    for L in text.splitlines(keepends=True):
        if L.startswith("  lj_source:"):
            continue
        if L.startswith("  lj_resolved:"):
            skip = True
            continue
        if skip and L.startswith("    - "):
            continue
        skip = False
        out.append(L)
    return "".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--forge", default=os.environ.get("FORGE_BIN"))
    ap.add_argument("--base-forge", default=None, help="本変更前のバイナリ (文字列リスト config のハッシュ不変の確認; 省略時は (H) を飛ばす)")
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    if not a.forge or not os.path.exists(a.forge):
        raise SystemExit("--forge BIN (or FORGE_BIN) is required")
    a.forge = os.path.abspath(a.forge)
    root = tempfile.mkdtemp(prefix="forge_species_lump_")
    C = Ctx(a, root)
    try:
        info = json.load(open(os.path.join(RUN0510, "prepare_info.json")))["species"]
        X = {k: float(v) for k, v in info["X"].items() if k != "H2O"}
        tot = sum(X.values())
        Xn = {k: v / tot for k, v in X.items()}
        Y = {k: float(v) for k, v in info["Y"].items() if k != "H2O"}
        ty = sum(Y.values())
        Yn = {k: v / ty for k, v in Y.items()}
        gen = yaml.load(open(os.path.join(RUN0510, "species_db.yaml")).read(), Loader=fs._StrSafeLoader)["MIXDRY"]
        ref = {"MW": float(gen["MW"]), "Tlo": float(gen["Tlo"]), "Tmid": float(gen["Tmid"]), "Thi": float(gen["Thi"]),
               "LJ_sigma": float(gen["LJ_sigma"]), "LJ_eps_kB": float(gen["LJ_eps_kB"]),
               "nasa9_low": [float(x) for x in gen["nasa9_low"]], "nasa9_high": [float(x) for x in gen["nasa9_high"]]}
        # 段 3 (#13-3): run_0510 の生成 DB は段 3 前のデータ (Ar の第 2 区間 = 単原子理想の低温係数の流用、第 3 区間なし) で作った。
        # 内蔵が CEA そのものになったので、参照の第 2 区間に Ar の差 x_Ar (a_Ar,CEA − a_Ar,旧) を足し (NASA-9 はモル分率加重で厳密に混合)、
        # 第 3 区間は構成種の CEA 第 3 区間のモル分率加重和を独立に作る。N2/O2/CO2 の先頭 2 区間は段 3 で不変。
        spd = {str(e["id"]): e for e in yaml.safe_load(open(os.path.join(REPO, "solver_density_cuda", "data", "species",
                                                                         "forge_species_v1.yaml")))["species"]}
        canon = {"AR": "Ar"}
        ar = spd["Ar"]["intervals"]
        dAr = [float(n_) - float(o_) for n_, o_ in zip(ar[1]["coeffs"], ar[0]["coeffs"])]   # 旧第 2 区間 = 第 1 区間の係数
        ref["nasa9_high"] = [h_ + Xn["AR"] * d_ for h_, d_ in zip(ref["nasa9_high"], dAr)]
        ref3 = [sum(Xn[k] * float(spd[canon.get(k, k)]["intervals"][2]["coeffs"][i]) for k in Xn) for i in range(9)]
        T = np.linspace(200.0, 6000.0, 1000)
        T6 = T[:-1]   # 6000 K ちょうどは区間選択の規約で決まる点 (段 3 で第 3 区間になる) なので 2 区間の参照との比較から外す

        # ---- V2 ----
        d = C.make(RUN0510, "[" + flow_lump("MIXDRY", Xn, "mole") + ", H2O]")
        rc, h, rec, err = C.resolve(d)
        check(rc == 0 and rec is not None and rec["consistent"] and rec["compat_recomputed"] == h,
              f"V2 lump config (run_0510 noneq, H2O condensing, mole basis) resolves; Python load_record recomputes the hash ({h[:16]})")
        if rec is None:
            print(err[-2000:])
            raise SystemExit(1)
        mix = rec["species"][0]
        mTb, mco = ivs(mix)
        rTb, rco = ivs(ref)
        worst = max([rel(mix[k], ref[k]) for k in ("MW", "LJ_sigma", "LJ_eps_kB")]
                    + [rel(x, y) for j in range(2) for x, y in zip(mco[j], rco[j])])
        check(worst <= 1e-12, f"V2 MIXDRY coefficients (first 2 intervals)/MW/LJ vs run_0510 species_db.yaml + x_Ar·Δa_Ar (#13-3): "
                              f"max rel {worst:.2e} (<= 1e-12)")
        w3 = max(rel(x, y) for x, y in zip(mco[2], ref3)) if len(mco) == 3 else float("inf")
        check(w3 <= 1e-12, f"V2 MIXDRY 3rd interval 6000–20000 K = Σ x_k a_k,CEA (#13-3): max rel {w3:.2e} (<= 1e-12)")
        for k in ("MW", "LJ_sigma", "LJ_eps_kB"):
            print(f"       {k}: solver {mix[k]!r} generated {ref[k]!r} rel {rel(mix[k], ref[k]):.2e}")
        check(mTb == [200.0, 1000.0, 6000.0, 20000.0], f"V2 MIXDRY breakpoints 200/1000/6000/20000 (#13-3; constituents are CEA 3-interval): {mTb}")
        compare_props(lambda t: (cp_mass(mix, t), h_mass(mix, t), s0_mass(mix, t)),
                      lambda t: (cp_mass(ref, t), h_mass(ref, t), s0_mass(ref, t)), "V2 cp/h/s0 solver lump vs generated (+x_Ar·Δa_Ar)", T6)
        # V2i: 構成種の質量分率加重和 (独立検算)
        mem = mix["lump"]["members"]
        Yk = [m["x"] * m["MW"] / mix["MW"] for m in mem]

        def members_sum(t):
            return (sum(y * cp_mass(m, t) for y, m in zip(Yk, mem)), sum(y * h_mass(m, t) for y, m in zip(Yk, mem)),
                    sum(y * s0_mass(m, t) for y, m in zip(Yk, mem)))
        compare_props(lambda t: (cp_mass(mix, t), h_mass(mix, t), s0_mass(mix, t)), members_sum,
                      "V2i cp/h/s0 solver lump vs mass-weighted constituents (s0 without mixing term)", T)

        # ---- R: 記録の中身 ----
        lp = mix["lump"]
        check(lp["basis"] == "mole" and [m["name"] for m in lp["members"]] == list(Xn)
              and all(m["fraction_input"] == Xn[m["name"]] for m in lp["members"])
              and all(rel(m["x"], Xn[m["name"]]) <= 1e-15 for m in lp["members"])
              and all(m["source"] == "builtin" for m in lp["members"]) and mix["source"] == "lump",
              f"R record carries the lump composition (basis mole, members {list(Xn)}, input = x (normalized), source builtin)")
        check("PROVISIONAL" in err and "lump MIXDRY" in err and "T=298.15 K: cp=" in err,
              "R startup log shows lump contents, provisional LJ and cp/h at reference temperatures")
        sig = fs.signature_from_record(rec)
        rec2 = json.loads(json.dumps(rec))
        rec2["species"][0]["lump"]["members"][1]["x"] *= (1 + 1e-9)
        diff = fs.compare_signatures(sig, fs.signature_from_record(rec2))
        check(len(diff) == 1 and "MIXDRY.lump.O2.x" in diff[0], f"R compare_signatures detects a changed lump mole fraction: {diff}")
        check(fs.compare_signatures(sig, fs.signature_from_record(rec)) == [], "R compare_signatures: same record -> no difference")

        # ---- M: basis mass ----
        dm = C.make(RUN0510, "[" + flow_lump("MIXDRY", Yn, "mass") + ", H2O]")
        rcm, hm, recm, errm = C.resolve(dm)
        if recm is None:
            check(False, f"M mass-basis lump resolves: {errm[-800:]}")
        else:
            mm = recm["species"][0]
            w = max([rel(mm[k], mix[k]) for k in ("MW", "LJ_sigma", "LJ_eps_kB")]
                    + [rel(x, y) for ra, rb in zip(ivs(mm)[1], ivs(mix)[1]) for x, y in zip(ra, rb)]
                    + ([0.0] if ivs(mm)[0] == ivs(mix)[0] else [float("inf")]))
            check(rcm == 0 and w <= 1e-12 and recm["species"][0]["lump"]["basis"] == "mass",
                  f"M basis: mass (Y from prepare_info) -> same synthesized MIXDRY as mole basis: max rel {w:.2e}")

        # ---- H: 文字列リストのハッシュ不変 ----
        cfgs = [("case44 run_0509 (external DB MIXDRY/H2O)", C.copy_run_config(RUN0509))]
        for sp in (["N2"], ["N2", "H2O"], ["H2O", "N2", "O2", "AR", "CO2"]):
            cfgs.append((f"builtin {sp}", C.make(RUN0509, "[" + ", ".join(f'"{s}"' for s in sp) + "]")))
        h0509 = None
        for tag, dd in cfgs:
            rcn, hn, recn, _ = C.resolve(dd)
            if tag.startswith("case44"):
                h0509 = hn
                check(rcn == 0 and hn.startswith("4378b7d78339ba27"), f"H {tag}: hash {hn[:16]} == 4378b7d78339ba27 (#3a/#4)")
            if a.base_forge:
                db_ = dd + "_base"
                shutil.copytree(dd, db_, ignore=shutil.ignore_patterns("resolved_species_*"))
                rcb, hb, recb, _ = C.resolve(db_, exe=os.path.abspath(a.base_forge))
                same = (rcn == 0 and rcb == 0 and hn == hb and recn is not None and recb is not None
                        and strip_lj_provenance(open(recn["path"]).read()) == strip_lj_provenance(open(recb["path"]).read()))
                if tag.startswith("case44"):
                    # 外部 DB だけで種が決まる config は段 3 (#13-3) でも不変
                    check(same, f"H {tag}: hash {hn[:16]} and record bytes (without the #14-L1 LJ provenance) identical to base binary ({hb[:16]})")
                else:
                    # 内蔵種の config は段 3 で内蔵が CEA そのものになったので、段 3 前のバイナリとはハッシュが違うのが正しい
                    # (--base-forge が段 3 後のバイナリなら従来どおりバイト一致を求める)
                    base_pre13_3 = rcb == 0 and recb is not None and any(
                        e.get("name") == "N2" and "Tbounds" not in e for e in recb["species"])
                    ok = (rcn == 0 and rcb == 0 and hn != hb) if base_pre13_3 else same
                    check(ok, f"H {tag}: hash {hb[:16]} -> {hn[:16]} ("
                              + ("changed by #13-3, base binary is pre-#13-3" if base_pre13_3
                                 else "identical to base binary, record without the #14-L1 LJ provenance") + ")")

        # ---- L: run_0509 を lump 記法に ----
        X9 = json.load(open(os.path.join(RUN0509, "prepare_info.json")))["species"]["X"]
        X9 = {k: float(v) for k, v in X9.items() if k != "H2O"}
        t9 = sum(X9.values())
        X9 = {k: v / t9 for k, v in X9.items()}
        dl = C.make(RUN0509, "[" + flow_lump("MIXDRY", X9, "mole") + ", H2O]")
        rcl, hl, recl, errl = C.resolve(dl)
        check(rcl == 0 and recl is not None and recl["species"][0].get("lump") and hl != h0509,
              f"L run_0509 config in lump notation (no speciesDBFile) resolves; hash {hl[:16]} != external-DB {str(h0509)[:16]} (expected)")

        # ---- N: 拒否 ----
        base = "[" + flow_lump("MIXDRY", Xn, "mole") + ", H2O]"

        def expect_fail(tag, species_yaml, needle, seed=RUN0510, db=None, drop_db=True):
            dd = C.make(seed, species_yaml, db=db, drop_db=drop_db)
            rcx, _, _, errx = C.resolve(dd)
            check(rcx != 0 and needle in errx, f"N {tag}: rejected ('{needle}')" + ("" if (rcx != 0 and needle in errx) else f" rc={rcx} err={errx[-600:]}"))

        n2 = yaml.safe_load(open(os.path.join(REPO, "solver_density_cuda", "data", "species", "forge_species_v1.yaml")))["species"][0]
        assert n2["id"] == "N2"
        tdb = {"N2B2": {"MW": n2["MW"], "LJ_sigma": 3.621, "LJ_eps_kB": 97.53, "Tlo": 200.0, "Tmid": 1000.0, "Thi": 6000.0,
                        "nasa9_low": n2["intervals"][0]["coeffs"], "nasa9_high": n2["intervals"][1]["coeffs"]},
               "TMID1500": {"MW": n2["MW"], "LJ_sigma": 3.621, "LJ_eps_kB": 97.53, "Tlo": 200.0, "Tmid": 1500.0, "Thi": 6000.0,
                            "nasa9_low": n2["intervals"][0]["coeffs"], "nasa9_high": n2["intervals"][1]["coeffs"]},
               "TESTGAS": {"MW": n2["MW"], "LJ_sigma": 3.621, "LJ_eps_kB": 97.53,
                           "nasa9_low": n2["intervals"][0]["coeffs"], "nasa9_high": n2["intervals"][1]["coeffs"]}}
        # B6: 区切りの違う構成種は和集合で合成 (#13-1)。内蔵 N2 は段 3 (#13-3) から 3 区間 (…6000/20000) で TMID1500 との和集合が
        # 4 区間になり拒否されるので、2 区間の N2 (N2B2, 内蔵の先頭 2 区間) を外部 DB で与える
        db6 = C.make(RUN0510, "[" + flow_lump("MIXDRY", {"N2B2": 0.8, "TMID1500": 0.2}, "mole") + ", H2O]", db=tdb)
        rc6, _, rec6, err6 = C.resolve(db6)
        m6 = rec6["species"][0] if rec6 is not None else {}
        check(rc6 == 0 and m6.get("Tbounds") == [200.0, 1000.0, 1500.0, 6000.0] and "union" in m6.get("lump", {}).get("synthesis", "")
              and rec6["schema"] == fs.SPECIES_RECORD_SCHEMA_NINT and rec6["consistent"],
              "B6 different breakpoints (external TMID1500, Tmid 1500 K) -> synthesized over the union 200/1000/1500/6000"
              + ("" if rc6 == 0 else f" rc={rc6} err={err6[-600:]}"))
        expect_fail("condensing species H2O inside the lump (condensation: 1)",
                    "[" + flow_lump("MIXDRY", {"N2": 0.9, "H2O": 0.1}, "mole") + ", H2O]", "condensing species")
        expect_fail("condensing species via alias WATER inside the lump",
                    "[" + flow_lump("MIXDRY", {"N2": 0.9, "WATER": 0.1}, "mole") + ", H2O]", "condensing species")
        expect_fail("lump name collides with built-in N2", "[" + flow_lump("N2", Xn, "mole") + ", H2O]", "collides with built-in")
        expect_fail("lump name collides with built-in AIR (case-insensitive 'air')", "[" + flow_lump("air", Xn, "mole") + ", H2O]", "collides with built-in")
        expect_fail("lump name collides with external DB species", "[" + flow_lump("TESTGAS", Xn, "mole") + ", H2O]",
                    "collides with speciesDBFile", db=tdb)
        expect_fail("zero fraction", "[" + flow_lump("MIXDRY", {**Xn, "CO2": 0.0}, "mole") + ", H2O]", "must be positive")
        expect_fail("negative fraction", "[" + flow_lump("MIXDRY", {**Xn, "O2": -0.1}, "mole") + ", H2O]", "must be positive")
        expect_fail("non-finite fraction", "[{name: MIXDRY, lump: {N2: .nan, O2: 0.2}, basis: mole}, H2O]", "must be positive")
        expect_fail("missing basis", "[{name: MIXDRY, lump: {N2: 0.8, O2: 0.2}}, H2O]", "needs 'basis")
        expect_fail("unknown basis", "[{name: MIXDRY, lump: {N2: 0.8, O2: 0.2}, basis: volume}, H2O]", "basis must be")
        expect_fail("unknown constituent", "[{name: MIXDRY, lump: {N2: 0.8, XENON: 0.2}, basis: mole}, H2O]", "constituent 'XENON' not found")
        expect_fail("duplicate constituent via alias (AR and Ar)", "[{name: MIXDRY, lump: {N2: 0.8, AR: 0.1, Ar: 0.1}, basis: mole}, H2O]",
                    "listed twice")
        expect_fail("unknown key in lump entry", "[{name: MIXDRY, lump: {N2: 0.8, O2: 0.2}, basis: mole, frac: 1}, H2O]", "unknown key")
        # 総和が 1 から外れる → 警告して通る
        dw = C.make(RUN0510, "[" + flow_lump("MIXDRY", X, "mole") + ", H2O]")
        rcw, hw, recw, errw = C.resolve(dw)
        check(rcw == 0 and "differs from 1 by >= 1e-3" in errw and recw is not None
              and max(rel(m["x"], Xn[m["name"]]) for m in recw["species"][0]["lump"]["members"]) <= 1e-15,
              f"N un-normalized X (sum {tot:.6f}) -> warning, normalized within the lump")
    finally:
        if a.keep:
            print(f"kept: {root}")
        else:
            shutil.rmtree(root, ignore_errors=True)
    print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
