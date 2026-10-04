#!/usr/bin/env python3
"""CEA `trans.inp` → forge 共通輸送データ (`solver_density_cuda/data/species/forge_transport_v1.yaml`) の生成器。

plans/active/thermophysics-solver-owned-species-db.md §4.3b・§4.3c・§5.1 #5t2 (段 1)。仕様 methods/thermophysics.md。

  python3 solver_density_cuda/tools/cea_trans_to_forge_transport.py [--trans .venv-cea/nasa_cea/trans.inp] [--out PATH]
  python3 solver_density_cuda/tools/cea_trans_to_forge_transport.py --check     # 配布ファイルが trans.inp から再生成したものと一字一句同じか

入力の読み方は CEA 本体 (cea2.f UTRAN) と同じ固定桁:
  見出し行 FORMAT (2A16,2X,A1,I1,A1,I1)  = 種名 1 (1–16 桁)・種名 2 (17–32 桁; 空なら単成分)・'V' nV・'C' nC
  データ行 FORMAT (1X,A1,2F9.2,4E15.8)   = 'V'|'C'・Tlo・Thi・A・B・C・D
  (E15.8 の指数部の空白は Fortran の既定どおり 0 とみなす: "0.61205763E 00" = 0.61205763E+00)
フィットの形: ln f = A ln T + B/T + C/T² + D (T [K])。単位: V (粘性) は f = η [μP]、C (熱伝導) は f = λ [μW/(cm·K)]。
区間の選び方 (cea2.f TRANIN の kt): 上端 Thi を下から順に見て T ≤ Thi となる最初の区間、どれにも入らなければ最後の区間
(= 下端より低い T は最初の区間、上端より高い T は最後の区間で外挿; 区間の境界ちょうどは下側の区間)。

出力の数値は Python の float (IEEE double) を repr で書く (C++ yaml-cpp / Python yaml が同じ double に戻る)。
"""
import argparse
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SOLVER = os.path.normpath(os.path.join(HERE, ".."))
REPO = os.path.normpath(os.path.join(SOLVER, ".."))
DEFAULT_TRANS = os.path.join(REPO, ".venv-cea", "nasa_cea", "trans.inp")
DEFAULT_OUT = os.path.join(SOLVER, "data", "species", "forge_transport_v1.yaml")
SCHEMA = "forge_transport_data_v1"


def _fnum(s):
    """Fortran の E/F 欄 (空白は 0 とみなす) を float に。"""
    t = s.replace(" ", "")
    if t == "":
        return 0.0
    return float(t)


def parse_trans(text):
    """trans.inp を読み、[{names: (a, b), nV, nC, reference, V: [[Tlo,Thi,A,B,C,D],...], C: [...]}] を返す (ファイル順)。"""
    lines = text.splitlines()
    if not lines or not lines[0].lower().startswith("transport property"):
        raise ValueError("trans.inp: first line must be 'transport property coefficients'")
    out = []
    i = 1
    while i < len(lines):
        h = lines[i]
        n1, n2 = h[0:16].strip(), h[16:32].strip()
        if n1 in ("end", "LAST"):
            break
        tag = h[34:38]
        if len(tag) != 4 or tag[0] != "V" or tag[2] != "C" or not tag[1].isdigit() or not tag[3].isdigit():
            raise ValueError(f"trans.inp line {i + 1}: header '{h}' does not have V<n>C<n> at columns 35-38")
        nv, nc = int(tag[1]), int(tag[3])
        if nv > 3 or nc > 3:
            raise ValueError(f"trans.inp line {i + 1}: more than 3 intervals (CEA limit)")
        rec = {"names": (n1, n2), "nV": nv, "nC": nc, "reference": h[40:].strip(), "V": [], "C": []}
        for k in range(nv + nc):
            L = lines[i + 1 + k]
            kind = L[1:2]
            if kind not in ("V", "C"):
                raise ValueError(f"trans.inp line {i + 2 + k}: expected 'V' or 'C' in column 2: '{L}'")
            row = [_fnum(L[2:11]), _fnum(L[11:20])] + [_fnum(L[20 + 15 * m:35 + 15 * m]) for m in range(4)]
            rec[kind].append(row)
        if len(rec["V"]) != nv or len(rec["C"]) != nc:
            raise ValueError(f"trans.inp line {i + 1}: V/C line counts do not match the header")
        for kind in ("V", "C"):
            iv = rec[kind]
            for m, r in enumerate(iv):
                if not r[0] < r[1]:
                    raise ValueError(f"trans.inp {n1}/{n2} {kind}: interval {m} has Tlo >= Thi")
                if m > 0 and not iv[m - 1][1] <= r[0]:
                    raise ValueError(f"trans.inp {n1}/{n2} {kind}: intervals overlap or are not increasing")
        out.append(rec)
        i += 1 + nv + nc
    return out


def _r(x):
    # PyYAML (YAML 1.1) は "1e-05" のように仮数に小数点の無い指数表記を float と解釈しないので "1.0e-05" にする
    s = repr(float(x))
    if "e" in s and "." not in s.split("e")[0]:
        m, e = s.split("e")
        s = m + ".0e" + e
    return s


def _rows(rows, indent):
    return "".join(f"{indent}- [{', '.join(_r(v) for v in r)}]\n" for r in rows)


def _q(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def render(recs, sha, trans_name):
    # 大小文字だけが違う種名の衝突を拒否 (canonical ID は大小文字を区別する; plan §4.1)
    seen = {}
    for r in recs:
        if r["names"][1]:
            continue
        k = r["names"][0].upper()
        if k in seen:
            raise ValueError(f"trans.inp: species '{r['names'][0]}' and '{seen[k]}' differ only in case")
        seen[k] = r["names"][0]
    o = []
    o.append("# =============================================================================\n")
    o.append("# forge 共通輸送データ (CEA trans.inp の種別フィットと相互作用フィット)。生成ファイル: 手で編集しない。\n")
    o.append("#   生成器 solver_density_cuda/tools/cea_trans_to_forge_transport.py (--check で再生成と照合)。\n")
    o.append("#   plans/active/thermophysics-solver-owned-species-db.md §4.3b・§4.3c・§5.1 #5t2、仕様 methods/thermophysics.md。\n")
    o.append("#   C++ はビルド時にバイナリへ埋め込む (cmake/embed_species_data.cmake; forge_species_v1.yaml と同じヘッダ)。\n")
    o.append("#\n")
    o.append("# スキーマ (forge_transport_data_v1):\n")
    o.append("#   species[].id           trans.inp の種名 (大小文字を区別)\n")
    o.append("#   species[].V / C        区間ごとに [Tlo, Thi, A, B, C, D]:  ln f = A ln T + B/T + C/T^2 + D  (T [K])\n")
    o.append("#                          V: f = 粘性 [microPoise] (= 1e-7 Pa s)、C: f = 熱伝導率 [microW/(cm K)] (= 1e-4 W/(m K))\n")
    o.append("#   interactions[].pair    異種の組 (trans.inp の並び)。V は相互作用粘性 eta_ij (単位・形は上と同じ)。C は CEA が使わない\n")
    o.append("#   区間の選び方 (cea2.f TRANIN, 行 5466 付近の kt): T <= Thi となる最初の区間、無ければ最後の区間\n")
    o.append("#                          (範囲外は最寄りの区間の式で外挿。区間の境界ちょうどは下側の区間)\n")
    o.append("# =============================================================================\n")
    o.append(f"schema: {SCHEMA}\n")
    o.append("provenance:\n")
    o.append("  cea_trans_inp:\n")
    o.append(f"    file: {_q(trans_name)}\n")
    o.append(f"    sha256: {_q(sha)}\n")
    o.append("    description: \"NASA Glenn CEA trans.inp (McBride & Gordon, NASA RP-1311; per-species references on each entry)\"\n")
    o.append("units:\n")
    o.append("  V: \"ln(viscosity [microPoise])\"\n")
    o.append("  C: \"ln(thermal conductivity [microW/(cm K)])\"\n")
    o.append("form: \"ln f = A ln T + B/T + C/T^2 + D\"\n")
    o.append("interval_selection: \"first interval with T <= Thi, else the last (cea2.f TRANIN kt); outside the data range the nearest interval is extrapolated\"\n")
    o.append("species:\n")
    for r in recs:
        if r["names"][1]:
            continue
        o.append(f"  - id: {_q(r['names'][0])}\n")
        o.append(f"    reference: {_q(r['reference'])}\n")
        o.append("    V:\n" + _rows(r["V"], "      ") if r["V"] else "    V: []\n")
        o.append("    C:\n" + _rows(r["C"], "      ") if r["C"] else "    C: []\n")
    o.append("interactions:\n")
    for r in recs:
        if not r["names"][1]:
            continue
        o.append(f"  - pair: [{_q(r['names'][0])}, {_q(r['names'][1])}]\n")
        o.append(f"    reference: {_q(r['reference'])}\n")
        o.append("    V:\n" + _rows(r["V"], "      ") if r["V"] else "    V: []\n")
        o.append("    C:\n" + _rows(r["C"], "      ") if r["C"] else "    C: []\n")
    return "".join(o)


def generate(trans_path):
    raw = open(trans_path, "rb").read()
    sha = hashlib.sha256(raw).hexdigest()
    recs = parse_trans(raw.decode("latin-1"))
    return render(recs, sha, os.path.basename(trans_path))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trans", default=DEFAULT_TRANS)
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--check", action="store_true", help="--out が trans.inp からの再生成と一致するか確認するだけ (書かない)")
    a = ap.parse_args()
    text = generate(a.trans)
    if a.check:
        cur = open(a.out, encoding="utf-8").read() if os.path.exists(a.out) else None
        if cur != text:
            print(f"[FAIL] {a.out} differs from regeneration from {a.trans}")
            return 1
        print(f"[PASS] {a.out} == regeneration from {a.trans}")
        return 0
    with open(a.out, "w", encoding="utf-8") as f:
        f.write(text)
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
