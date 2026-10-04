#!/usr/bin/env bash
# 二相拡散 A/B (plan condensation-two-phase-transport §5.1 #1b; 設計メモ
# notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md §14.4) の 2 run を準備する。**forge は起動しない**。
#
#   prepare_twophase_ab.sh SRC_RUN RUN_A RUN_B [NSTEP=4000] [OUTINT=200]
#     環境変数 SRC_STEP (既定 48000): 共通初期場にする SRC_RUN/res_<SRC_STEP>.h5
#
# やること (A・B で同じ扱い。違いは B の condensation.condTwoPhaseDiffusion: 1 だけ):
#   1. RUN_A / RUN_B を新規作成 (既にあれば失敗して止まる)。SRC_RUN の solverConfig.yaml・bcondConfig.yaml・probe.yaml・
#      species_db.yaml と、solverConfig の meshFileName / valueFileName の h5 を複製する。
#   2. species_db.yaml から H2O を外す (このブランチでは凝縮 ON の H2O は内蔵項目を使う。run_0482 の DB の H2O は MW が
#      内蔵と違い拒否される; メモ §14.4 の 2)。外した後に他の種の内容が変わっていないことを検査する。
#   3. solverConfig.yaml: time.last.nStepOuter = NSTEP、time.outStepInterval = OUTINT (A・B 共通)。B だけ
#      condensation.condTwoPhaseDiffusion: 1 (condTwoPhaseRelax は書かない = 既定 1)。書換え後に読み直し、
#      意図した差分以外が無いことを検査する (anchor・重複キー・節の取り違えの防止)。
#   4. 共通初期場: SRC_RUN/res_<SRC_STEP>.h5 を restart_field.py で各 run の値ファイルへ index コピー
#      (--force-species: SRC は種の属性を持たない旧 run、宛先は H2O を外した DB で種ハッシュが変わるため。属性なし=未検証のまま)。
#      **値ファイルに凝縮モーメントの保存量 (rog_0, roQ2_0, roQ1_0, roQ0_0) のデータセットが無いと restart_field.py は
#      写さず (DST にある名前だけ写す)、forge は液 0 の dry restart になる** (variables.cpp の凝縮モーメント読込)。
#      run_0482 の値ファイルはこれらを持たないので、SRC にあって DST に無い ro*/rog_*/roQ*_* を先に作ってから写す。
#   5. 検査: A と B の値ファイルの VALUE が全データセットで一致すること、凝縮モーメントが SRC と一致すること。
#
# 起動 (親が AWS で行う。両 run とも同じ許可が要る):
#   FORGE_ALLOW_UNVERIFIED_SPECIES=1 FORGE_BIN=<このブランチのクリーンビルド> solver_density_cuda/tools/run_case.sh RUN_A
#   (RUN_B も同じ)。B の起動ログに `[twophase] condTwoPhaseDiffusion 1: ...` と、終了時に `[twophase-audit] VERDICT:` が出る。
set -euo pipefail
if [ $# -lt 3 ]; then
  sed -n '2,30p' "$0"; exit 2
fi
SRC_RUN="$(cd "$1" && pwd)"; RUN_A="$2"; RUN_B="$3"; NSTEP="${4:-4000}"; OUTINT="${5:-200}"
SRC_STEP="${SRC_STEP:-48000}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
TOOLS="$ROOT/solver_density_cuda/tools"
SRC_RES="$SRC_RUN/res_${SRC_STEP}.h5"
[ -f "$SRC_RES" ] || { echo "SRC の res が無い: $SRC_RES"; exit 1; }
for d in "$RUN_A" "$RUN_B"; do
  [ -e "$d" ] && { echo "既に存在する (上書きしない): $d"; exit 1; }
done
# OUTINT が NSTEP を割り切らないと最終状態の res が出ない
if [ $((NSTEP % OUTINT)) -ne 0 ]; then echo "NSTEP ($NSTEP) は OUTINT ($OUTINT) で割り切れること"; exit 1; fi

for d in "$RUN_A" "$RUN_B"; do mkdir "$d"; done
RUN_A="$(cd "$RUN_A" && pwd)"; RUN_B="$(cd "$RUN_B" && pwd)"

python3 - "$SRC_RUN" "$RUN_A" "$RUN_B" "$NSTEP" "$OUTINT" <<'PY'
import copy, os, shutil, sys, yaml
src, ra, rb, nstep, outint = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
cfg0 = yaml.safe_load(open(os.path.join(src, "solverConfig.yaml")))
m = cfg0["mesh"]
h5s = sorted({m["meshFileName"], m.get("valueFileName", m["meshFileName"])})
db0 = yaml.safe_load(open(os.path.join(src, cfg0["physProp"]["speciesDBFile"])))
if "H2O" not in db0:
    sys.exit("species_db.yaml に H2O が無い (想定と違う; 手で確認すること)")
db = {k: v for k, v in db0.items() if k != "H2O"}

def walk_diff(a, b, path=""):
    """a と b の差分パスを列挙 (dict は再帰、その他は ==)。"""
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b), key=str):
            p = f"{path}.{k}" if path else str(k)
            if k not in a: out.append(p + " (追加)")
            elif k not in b: out.append(p + " (削除)")
            else: out += walk_diff(a[k], b[k], p)
        return out
    return [] if a == b else [path]

for run, twophase in ((ra, False), (rb, True)):
    for f in ["bcondConfig.yaml", "probe.yaml"] + h5s:
        shutil.copy2(os.path.join(src, f), os.path.join(run, f))
    dbf = os.path.join(run, cfg0["physProp"]["speciesDBFile"])
    yaml.safe_dump(db, open(dbf, "w"), sort_keys=False)
    if walk_diff(db0, yaml.safe_load(open(dbf))) != ["H2O (削除)"]:
        sys.exit(f"{dbf}: H2O 以外の差がある: {walk_diff(db0, yaml.safe_load(open(dbf)))}")
    cfg = copy.deepcopy(cfg0)
    cfg["time"]["last"]["nStepOuter"] = nstep
    cfg["time"]["outStepInterval"] = outint
    want = set()
    if cfg0["time"]["last"]["nStepOuter"] != nstep:
        want.add("time.last.nStepOuter")
    if "outStepInterval" not in cfg0["time"]:
        want.add("time.outStepInterval (追加)")
    elif cfg0["time"]["outStepInterval"] != outint:
        want.add("time.outStepInterval")
    if twophase:
        if "condTwoPhaseDiffusion" in cfg["condensation"]:
            sys.exit("SRC の solverConfig に既に condTwoPhaseDiffusion がある (A が現行にならない)")
        cfg["condensation"]["condTwoPhaseDiffusion"] = 1
        want |= {"condensation.condTwoPhaseDiffusion (追加)"}
    elif "condTwoPhaseDiffusion" in cfg.get("condensation", {}):
        sys.exit("SRC の solverConfig に condTwoPhaseDiffusion がある (A が現行にならない)")
    if "bndFirstOrder" in cfg.get("mesh", {}):
        sys.exit("mesh.bndFirstOrder が入っている (使用禁止; AGENTS.md)")
    if cfg.get("turbulence", {}).get("wallTreatmentSST", 0) != 0:
        sys.exit("wallTreatmentSST != 0 (壁関数は使用禁止)")
    cf = os.path.join(run, "solverConfig.yaml")
    yaml.safe_dump(cfg, open(cf, "w"), sort_keys=False, default_flow_style=None)
    got = set(walk_diff(cfg0, yaml.safe_load(open(cf))))
    if got != want:
        sys.exit(f"{cf}: 想定外の差分 {sorted(got ^ want)}")
    print(f"[prepare] {os.path.basename(run)}: config 差分 {sorted(got)}; species_db から H2O を削除")
PY

# 4. 共通初期場 (凝縮モーメントのデータセットを先に作ってから index コピー)
VALF="$(python3 -c "import yaml,sys; m=yaml.safe_load(open('$SRC_RUN/solverConfig.yaml'))['mesh']; print(m.get('valueFileName', m['meshFileName']))")"
for d in "$RUN_A" "$RUN_B"; do
  python3 - "$SRC_RES" "$d/$VALF" <<'PY'
import re, sys, h5py, numpy as np
src, dst = sys.argv[1], sys.argv[2]
with h5py.File(src, "r") as s, h5py.File(dst, "r+") as d:
    sv, dv = s["VALUE"], d["VALUE"]
    n = dv["ro"].shape[0]
    add = [k for k in sv if re.fullmatch(r"ro(g|Q[0-9])_[0-9]+", k) and k not in dv]
    for k in add:
        dv.create_dataset(k, data=np.zeros(n, dtype=dv["ro"].dtype))   # 値は直後の restart_field.py が SRC から写す
    print(f"[prepare] {dst}: 凝縮モーメントのデータセットを作成 {add}")
PY
  python3 "$TOOLS/restart_field.py" "$SRC_RES" "$d/$VALF" --dst-run "$d" --force-species
done

# 5. 検査: A と B の初期場が一致し、凝縮モーメントが SRC と一致する
python3 - "$SRC_RES" "$RUN_A/$VALF" "$RUN_B/$VALF" <<'PY'
import sys, h5py, numpy as np
src, fa, fb = sys.argv[1:]
with h5py.File(src, "r") as s, h5py.File(fa, "r") as a, h5py.File(fb, "r") as b:
    va, vb, sv = a["VALUE"], b["VALUE"], s["VALUE"]
    if set(va) != set(vb): sys.exit(f"A と B の VALUE の集合が違う {sorted(set(va) ^ set(vb))}")
    bad = [k for k in va if not np.array_equal(np.asarray(va[k]), np.asarray(vb[k]))]
    if bad: sys.exit(f"A と B の初期場が違う: {bad}")
    need = ["ro", "roUx", "roUy", "roe", "roK", "roOmega", "roY0", "roY1", "rog_0", "roQ2_0", "roQ1_0", "roQ0_0"]
    miss = [k for k in need if k not in va]
    if miss: sys.exit(f"値ファイルに無い保存量: {miss}")
    nd = [k for k in need if not np.array_equal(np.asarray(va[k]), np.asarray(sv[k]))]
    if nd: sys.exit(f"SRC と一致しない: {nd}")
    g = np.asarray(va["rog_0"], np.float64)/np.asarray(va["ro"], np.float64)
    print(f"[prepare] 検査 OK: A・B の VALUE {len(va)} 量が一致、保存量 {len(need)} 量が SRC とビット一致; g max {g.max():.6g}, g>1e-6 の節点 {int((g>1e-6).sum())}")
PY

for d in "$RUN_A" "$RUN_B"; do
  {
    echo "plan condensation-two-phase-transport §5.1 #1b (二相拡散 A/B; prepare_twophase_ab.sh)"
    echo "IC: $SRC_RES を restart_field.py --force-species で index コピー (凝縮モーメント rog_0/roQ*_0 を含む; A・B 同一)"
    echo "edits: species_db.yaml から H2O を削除 (内蔵項目を使う); time.last.nStepOuter=$NSTEP; time.outStepInterval=$OUTINT"
    [ "$d" = "$RUN_B" ] && echo "       condensation.condTwoPhaseDiffusion=1 (B)" || echo "       (A: 現行; condTwoPhaseDiffusion なし)"
    echo "launch: FORGE_ALLOW_UNVERIFIED_SPECIES=1 (A・B 同じ)"
  } > "$d/IC_FROM.txt"
done
echo "[prepare] 完了: $RUN_A / $RUN_B (forge は起動していない)"
