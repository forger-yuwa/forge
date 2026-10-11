"""plan architecture-float-state-double-geometry §6.34 の前提 1 (SASS の照合)。
使い方: qa_sass_cmp.py <基準 forge の sass> <新しい forge の sass> <出力 json> [--variants <dir>]
- 関数ごとに SASS を切り、無名名前空間のパス由来ハッシュ (_GLOBAL__N__<hex>_、_cu_<hex>) を正規化して比べる。
- enforceAxisSymmetry_d はシグネチャが変わるので名前が変わる (許す)。それ以外で一致しない関数は、--variants の dir にある
  <関数を含む .cu の名前>_base_<n>.sass (基準のソースを同じフラグで再コンパイルした SASS) のどれかと一致すれば
  「nvcc のコード生成の揺れ」として一致とみなす (§6.34 の前提 1 の追記、2026-10-11)。それでも外れたら不一致として残す。"""
import sys, re, json, collections, hashlib
from pathlib import Path
def norm(t):
    t = re.sub(r"_GLOBAL__N__[0-9a-f]{8}_", "_GLOBAL__N__X_", t)
    return re.sub(r"(_cu_)[0-9a-f]{8}", r"\1X", t)
def split(path):
    funcs = collections.defaultdict(list); cur = None; buf = []
    for line in norm(Path(path).read_text(errors="replace")).splitlines(keepends=True):
        m = re.match(r"\s*Function : (\S+)", line)
        if m:
            if cur is not None: funcs[cur].append("".join(buf))
            cur = m.group(1); buf = []; continue
        if line.startswith("Fatbin") or re.match(r"\s*code for sm_", line) or line.startswith("\t\t....."):
            if cur is not None: funcs[cur].append("".join(buf)); cur = None; buf = []
            continue
        if cur is not None: buf.append(line)
    if cur is not None: funcs[cur].append("".join(buf))
    return {k: sorted(hashlib.sha256(x.encode()).hexdigest() for x in v) for k, v in funcs.items()}
a = split(sys.argv[1]); b = split(sys.argv[2]); out = Path(sys.argv[3])
vdir = Path(sys.argv[sys.argv.index("--variants") + 1]) if "--variants" in sys.argv else None
axis = lambda k: "enforceAxisSymmetry_d" in k
onlyA = sorted(k for k in set(a) - set(b)); onlyB = sorted(k for k in set(b) - set(a))
diff = sorted(k for k in set(a) & set(b) if a[k] != b[k])
resolved, unresolved = {}, []
variants = {}
if vdir and vdir.exists():
    for f in sorted(vdir.glob("*_base_*.sass")):
        for k, h in split(f).items(): variants.setdefault(k, set()).add(tuple(h))
for k in diff:
    if tuple(b[k]) in variants.get(k, set()): resolved[k] = f"基準のソースの再コンパイル {len(variants[k])} 通りのうちの 1 つと一致"
    else: unresolved.append(k)
bad_only = [k for k in onlyA + onlyB if not axis(k)]
res = dict(n_base=len(a), n_new=len(b), n_identical=sum(1 for k in set(a) & set(b) if a[k] == b[k]),
           only_base=onlyA, only_new=onlyB, differ=diff, differ_resolved_as_nvcc_variant=resolved,
           sass_mismatch_except_axis=sorted(unresolved + bad_only), variants_dir=str(vdir) if vdir else None)
print(json.dumps(res, ensure_ascii=False, indent=1))
old = json.loads(out.read_text()) if out.exists() else {}
old.update(res); out.write_text(json.dumps(old, ensure_ascii=False, indent=1))
