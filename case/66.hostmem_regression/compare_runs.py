#!/usr/bin/env python3
"""case/66 回帰の比較 (plan architecture-solver-host-memory §6 の (a)(b)(c) と出力互換性)。AWS 上で回す (h5py・numpy)。

    # 構成ごと (registry.tsv から run を選ぶ。base だけなら同ビルド内のばらつきだけを出す)
    python3 compare_runs.py --cfg c36node [--base-build base] [--new-build new] [--out report.txt]
    python3 compare_runs.py --all [--new-build new] --out-dir compare/       # 全構成 (構成ごとの報告 + summary.txt)
    # run を直接指定 (予定 N は --steps で与える。--reps 既定 3、--require で追加の必要ファイル)
    python3 compare_runs.py --base RUN RUN RUN [--new RUN RUN RUN] --steps N [--out-interval K] [--require F ...]
    # 2 ファイルの単純比較 (dual-time の分割 res_100 と連続 res_200 など)
    python3 compare_runs.py --diff2 A.h5 B.h5

完全性検査 (plan §5.1 #6、2026-10-07 result レビュー M1。A・B とも比較の前に構成ごとに行い、満たさなければ**その構成を FAIL**、
値の比較はしない): 予定は matrix_spec.planned() (直接指定は --steps 等) から決め、run に「存在するもの」からは決めない。
  - 予定の反復数 (base・new とも、既定 3) の run がそろう (registry の rep 1..reps が 1 本ずつ、.exclude を除く。ディレクトリが実在)
  - 各 run が実行成功 (.state が done・rc 0 [変換器は既知の終了時 rc 1 も可、converted.h5 に /VALUE・/MESH、written=0 は不可]・NANCHECK: PASS)
  - 必要ファイル (初期出力 res_0・res_{k·out}・最終出力 res_N・残差 CSV・構成ごとの境界出力/probe/診断 CSV) がある
  - 最終出力の step = 予定 N (res_N がある・N より後の res が無い・solverConfig の nStepOuter = N・残差 CSV の step が 0..N−1 で
    最後の行が outer_end・境界出力の属性 step = N)。物理時刻 (CHECKPOINT の totalTime/dt 等の属性) が全 run で一致
  - 残差 CSV の必須列 (step・inner・phase・rms_ro/roUx/roUy/roUz/roe、SST なら rms_roK/roOmega、遷移なら rms_roGamma/roReth)
  - 比較した量の数が 0 なら FAIL

判定の定義 (§6。結果を見てから変えない):
  (a) step 0 の残差行 (residual_history.csv の最初の step の全行) を全列ビット一致で比べる。base 反復の間で値が割れる列
      (既知の 1 ulp の非決定性、atomicAdd) は「new の値が base で観測した値のどれか」で可。
      初期出力 (最初の res_*.h5 = writeInitialOutputs、残差の組立より前) では**保存量・原始量・幾何量だけ**をビット一致で比べる
      (残差・補正量・勾配・リミッタ・診断は初期出力では比べない。level 2 の初期出力には未初期化の残差が入る: §6 M4)。
  (b) 最後の出力 (最終 res_*.h5・境界出力・CSV 出力・残差履歴) の各データセット (列) について
      m(A,B) = max|A−B| / max|A| (A は組の前側、max|A| = 0 なら差が 0 のとき 0、そうでなければ inf)。
      S = 同ビルド内ペア差 (base 3 回の 3 対 + new 3 回の 3 対) の最大、D = base×new の 9 対の最大。合格は D ≤ 2·S
      (S = 0 のときは D = 0)。new が無いときは S_base だけを出す (= 基準のばらつき)。
  (c) 各 run の NaN/Inf (run_matrix.py が残した NANCHECK.txt と、ここでの再検査)。
  出力互換: 全 h5 のデータセット集合・shape・dtype・属性 (値も) が一致すること。属性は全 run で比べる。
  ログ: 起動時の名前登録・checkpoint 復元・遷移初期化・ψ 退避件数・警告の行が全 run で一致すること。
"""
import argparse
import csv
import fnmatch
import glob
import itertools
import os
import time
import re
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import matrix_spec as ms  # noqa: E402

# 初期出力で比べる量 (保存量・原始量・幾何量)。/VALUE 以外 (/MESH・/CHECKPOINT・/BCONDS 等) は幾何・状態として全部比べる。
T0_VALUE = re.compile(
    r"^(ro|roUx|roUy|roUz|roe|roK|roOmega|roY\d+|roXi|roGamma|roReth|ro[gQ]\w*_\d+|"
    r"P|T|Ux|Uy|Uz|k|omega|Y\d+|Xi|h0|[gQ][0-2]?_\d+|wall_dist|volume|ccx|ccy|ccz)$")
# 比較対象のファイル (run が書いたもの)。入力 (INPUT_FILES) と xmf・ログは除く。
OUT_GLOBS = ["*.h5", "*.csv", "*.out"]
# mem_samples.csv はハーネスの採取 (run_matrix.py --memwatch) で forge の出力ではない
SKIP_FILES = {"residual_history.png", "mem_samples.csv"}
# ログから拾う行 (時間・速度を含む行は拾わない)
LOG_PAT = re.compile(r"(registered|restored|history|\[variables\]|snapshot:|\[FORGE_OUT_RESIDUALS\]|\[FORGE_RESID_SNAP\]|"
                     r"psi-dualeval\] ON|WARNING|warning|警告|ignored|無視|lineImplicit|line-implicit|"
                     r"\[transition\]|transition .* (read|initiali)|output: |'output')", re.I)
# 時間・速度の行と、FORGE_MEMLOG の計測行 ([memlog] は RSS/HWM の実測値が入り、"lineImplicit" 等で LOG_PAT に掛かる) は比べない
LOG_DROP = re.compile(r"(ms/step|elapsed|eta |wall|Time = |sec|秒|\[memlog\])", re.I)
# 変更後ビルドが設計どおり新しく出す情報行 (ログ一致の対象から外す。理由は README「判定」)
LOG_NEW_INFO = re.compile(r"\[variables\] host cell arrays \(gpu: \d\): \d+ of \d+ registered")


# ------------------------------------------------------------------ run の選択
def read_registry(root=HERE):
    p = os.path.join(root, "registry.tsv")
    if not os.path.exists(p):
        return []
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def runs_for(cfg, build):
    out = []
    for r in read_registry():
        if r["cfg"] != cfg or r["build"] != build:
            continue
        d = os.path.join(HERE, r["run"])
        st = open(os.path.join(d, ".state")).read() if os.path.exists(os.path.join(d, ".state")) else ""
        if os.path.exists(os.path.join(d, ".exclude")):   # 破棄予定 (起動失敗など) は比較に使わない
            continue
        if st.startswith("done"):
            out.append(d)
    return out


def input_files(d):
    p = os.path.join(d, "INPUT_FILES")
    return set(open(p).read().split()) if os.path.exists(p) else set()


def out_files(d):
    inp = input_files(d)
    fs = set()
    for g in OUT_GLOBS:
        for p in glob.glob(os.path.join(d, g)):
            b = os.path.basename(p)
            if b not in inp and b not in SKIP_FILES:
                fs.add(b)
    return fs


def res_steps(d):
    st = []
    for p in glob.glob(os.path.join(d, "res_*.h5")):
        m = re.match(r"res_(\d+)\.h5$", os.path.basename(p))
        if m:
            st.append(int(m.group(1)))
    return sorted(st)


# ------------------------------------------------------------------ 完全性検査 (plan §5.1 #6、result レビュー M1)
# 予定 (plan) は matrix_spec.plan_from() の dict: cfg・kind・reps・N・out・required。
RES_REQUIRED_COLS = ("step", "inner", "phase", "rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe")
RES_SST_COLS = ("rms_roK", "rms_roOmega")             # main.cpp scalarResidualEnabled (SST) のとき
RES_TRANSITION_COLS = ("rms_roGamma", "rms_roReth")   # 遷移モデルが有効なとき (check_convergence.py と同じ)
# 最終出力から拾う「物理時刻・step」の属性 (全 run で一致すること。step / step_abs は境界出力の属性で、step は N と一致すること)
TIME_ROOT_ATTRS = ("step", "step_abs", "time", "totalTime")
TIME_CKPT_ATTRS = ("totalTime", "dt")


def select_runs(cfg, build, plan, root=HERE):
    """registry から構成 cfg・ビルド build の予定の反復 (rep 1..reps) を 1 本ずつ選ぶ。
    返り値 (run ディレクトリの list, FAIL の理由の list, 情報の list)。.exclude の run は使わない。
    ある rep が無い・除外されていない run が同じ rep に複数あるのは FAIL (runs_for と違い、終わっていない run も黙って落とさない)。"""
    dirs, why, notes = [], [], []
    by_rep = {}
    for r in read_registry(root):
        if r["cfg"] != cfg or r["build"] != build:
            continue
        d = os.path.join(root, r["run"])
        if os.path.exists(os.path.join(d, ".exclude")):
            notes.append(f"{r['run']}: 除外 (.exclude) — 比較に使わない")
            continue
        try:
            k = int(r["rep"])
        except (TypeError, ValueError):
            why.append(f"{r['run']}: registry の rep が整数でない ({r.get('rep')!r})")
            continue
        by_rep.setdefault(k, []).append(d)
    for k in range(1, plan["reps"] + 1):
        c = by_rep.get(k, [])
        if not c:
            why.append(f"{build} r{k} の run が registry に無い (除外分を除く; 予定 {plan['reps']} 本)")
        elif len(c) > 1:
            why.append(f"{build} r{k} が複数ある ({', '.join(os.path.basename(x) for x in c)}) — どれを使うか決まらない")
        else:
            dirs.append(c[0])
    extra = sorted(k for k in by_rep if not 1 <= k <= plan["reps"])
    if extra:
        notes.append(f"{build}: 予定外の反復 r{extra} は比較に使わない")
    return dirs, why, notes


def parse_state(d):
    """run_matrix.py のワーカーが書く .state ("done rc=0 nan=PASS [written=1] wall=…s") を読む。無ければ None。"""
    p = os.path.join(d, ".state")
    if not os.path.exists(p):
        return None
    raw = open(p).read().strip()
    m_rc = re.search(r"\brc=(-?\d+)", raw)
    m_w = re.search(r"\bwritten=(\d+)", raw)
    return dict(raw=raw, done=raw.startswith("done"), rc=int(m_rc.group(1)) if m_rc else None,
                written=int(m_w.group(1)) if m_w else None)


def _find_key(obj, key):
    """入れ子の dict から最初に見つかった key の値 (無ければ None)。"""
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]
        for v in obj.values():
            r = _find_key(v, key)
            if r is not None:
                return r
    return None


def solver_config_info(d):
    """run の solverConfig.yaml から nStepOuter・SST・遷移の有無。読めなければ None。PyYAML が無ければ正規表現で読む。"""
    p = os.path.join(d, "solverConfig.yaml")
    if not os.path.exists(p):
        return None
    txt = open(p, encoding="utf-8", errors="replace").read()
    try:
        import yaml
        c = yaml.safe_load(txt) or {}
        n = _find_key(c, "nStepOuter")
        tu = c.get("turbulence") or {}
        model = str(tu.get("model", "")).strip().lower()
        sst = model == "sst" or (tu.get("LESorRANS") == 2 and tu.get("RANSmodel") == 1)
        tr = str(tu.get("transition", "none")).strip().lower()
        return dict(N=int(n) if n is not None else None, sst=sst, transition=tr not in ("none", "0", "false", ""))
    except ImportError:
        t = re.sub(r"#.*", "", txt)
        m = re.search(r"\bnStepOuter\s*:\s*(\d+)", t)
        sst = bool(re.search(r"\bmodel\s*:\s*[\"']?sst\b", t, re.I)) or bool(
            re.search(r"\bLESorRANS\s*:\s*2\b", t) and re.search(r"\bRANSmodel\s*:\s*1\b", t))
        mt = re.search(r"\btransition\s*:\s*[\"']?([A-Za-z0-9_]+)", t)
        return dict(N=int(m.group(1)) if m else None, sst=sst,
                    transition=bool(mt) and mt.group(1).lower() not in ("none", "0", "false"))


def _required_present(d, name):
    if any(ch in name for ch in "*?["):
        return bool(fnmatch.filter(os.listdir(d), name))
    return os.path.exists(os.path.join(d, name))


def check_run(d, plan):
    """1 本の run の完全性。返り値 (FAIL の理由の list [空なら合格], 最終出力の物理時刻・step の属性 dict)。"""
    if not os.path.isdir(d):
        return [f"ディレクトリが無い ({d})"], {}
    why, tinfo = [], {}
    kind, n = plan["kind"], plan.get("N")
    # ---- 実行成功 (.state の rc と NANCHECK)
    st = parse_state(d)
    if st is None:
        why.append(".state が無い (ワーカーが終了を記録していない)")
    elif not st["done"]:
        why.append(f"終了していない (.state: {st['raw']})")
    elif kind == "convert":
        # 変換器は終了時の cudaFree で exit 1 になる既知の罠 (出力は完全)。完了は rc でなく出力 h5 で判定する (run_matrix.py の扱いのまま):
        # written=0 は FAIL。written の無い .state (記録を足す前のワーカー) は下の converted.h5 の検査 (/VALUE・/MESH) だけで判定する
        if st["written"] == 0:
            why.append(f"変換結果が書かれていない (.state: {st['raw']})")
        if st["rc"] not in (0, 1):
            why.append(f"終了コード {st['rc']} (変換器で許すのは 0 と既知の終了時 1 だけ)")
    elif st["rc"] != 0:
        why.append(f"終了コード {st['rc']} (.state: {st['raw']})")
    p = os.path.join(d, "NANCHECK.txt")
    if not os.path.exists(p):
        why.append("NANCHECK.txt が無い")
    else:
        lines = open(p).read().strip().splitlines()
        last = lines[-1] if lines else ""
        if last.strip() != "NANCHECK: PASS":
            why.append(f"NaN 検査が合格でない ({last or '空'})")
    # ---- 必要ファイル
    if kind == "forge" and n is None:
        why.append("予定 N が無い (構成表 matrix_spec か --steps で与える)")
    miss = [f for f in plan["required"] if not _required_present(d, f)]
    if miss:
        why.append(f"必要ファイルが無い: {', '.join(miss)}")
    if kind == "convert":
        try:
            with h5py.File(os.path.join(d, "converted.h5"), "r") as f:
                if "VALUE" not in f or "MESH" not in f:
                    why.append("converted.h5 に /VALUE か /MESH が無い")
        except Exception as e:  # 開けない h5 は合格にしない
            if "converted.h5" not in miss:
                why.append(f"converted.h5 が開けない ({e})")
        return why, tinfo
    if n is None:
        return why, tinfo
    # ---- 最終出力の step = 予定 N (存在する最大 step を最終と見なさない)。入力として複製した h5 (INPUT_FILES) は数えない
    inp = input_files(d)
    steps = [s for s in res_steps(d) if f"res_{s}.h5" not in inp]
    late = [s for s in steps if s > n]
    if late:
        why.append(f"予定 N {n} より後の出力がある (res_{late}.h5) — 予定と違う run")
    if n not in steps:
        why.append(f"最終出力 res_{n}.h5 が無い (存在する res の step: {steps or 'なし'})")
    sc = solver_config_info(d)
    if sc is None:
        why.append("solverConfig.yaml が無い")
    elif sc["N"] != n:
        why.append(f"solverConfig.yaml の nStepOuter {sc['N']} が予定 N {n} と違う")
    # ---- 残差 CSV: 必須列と step の範囲
    rp = os.path.join(d, "residual_history.csv")
    if os.path.exists(rp):
        with open(rp) as f:
            rd = csv.reader(f)
            hdr = [h.strip() for h in next(rd, [])]
            rows = [r for r in rd if r]
        need = list(RES_REQUIRED_COLS)
        if sc and sc["sst"]:
            need += RES_SST_COLS
        if sc and sc["transition"]:
            need += RES_TRANSITION_COLS
        lack = [c for c in need if c not in hdr]
        if lack:
            why.append(f"residual_history.csv に必須の列が無い: {', '.join(lack)}")
        if not rows:
            why.append("residual_history.csv にデータ行が無い")
        elif "step" in hdr:
            i = hdr.index("step")
            try:
                ss = {int(r[i]) for r in rows}
                if ss != set(range(n)):
                    why.append(f"residual_history.csv の step が 0..{n - 1} でない (最小 {min(ss)}、最大 {max(ss)}、{len(ss)} 種類)")
            except (ValueError, IndexError):
                why.append("residual_history.csv の step 列が整数でない")
            if "phase" in hdr:
                j = hdr.index("phase")
                if len(rows[-1]) <= j or rows[-1][j].strip() != "outer_end":
                    why.append(f"residual_history.csv の最後の行が outer_end でない ({rows[-1][:3]})")
    # ---- 最終出力の物理時刻・step の属性
    finals = sorted(b for b in (os.path.basename(x) for x in glob.glob(os.path.join(d, f"*_{n}.h5"))) if b not in inp)
    for fn in finals:
        try:
            with h5py.File(os.path.join(d, fn), "r") as f:
                for k in TIME_ROOT_ATTRS:
                    if k in f.attrs:
                        tinfo[f"{fn}:/@{k}"] = _attr(f.attrs[k])
                if "CHECKPOINT" in f:
                    for k in TIME_CKPT_ATTRS:
                        if k in f["CHECKPOINT"].attrs:
                            tinfo[f"{fn}:/CHECKPOINT@{k}"] = _attr(f["CHECKPOINT"].attrs[k])
        except Exception as e:
            why.append(f"{fn} が開けない ({e})")
    for k, v in tinfo.items():
        if k.endswith("/@step") and v != n:
            why.append(f"{k} = {v} が予定 N {n} と違う")
    return why, tinfo


def completeness(plan, base, new, pre_why=(), need_new=True):
    """構成の完全性。base / new は run ディレクトリの list (選択済み)。pre_why は選択の段階の FAIL 理由。
    返り値 dict(verdict "PASS"/"FAIL", reasons, lines)。"""
    why = list(pre_why)
    lines = [f"[完全性] 予定: kind {plan['kind']}、反復 {plan['reps']} 本 (base"
             + (" と new" if need_new else "") + f")、N {plan.get('N')}、出力間隔 {plan.get('out')}、"
             f"必要ファイル {', '.join(plan['required'])}"]
    for label, dirs in (("base", base), ("new", new)):
        if label == "new" and not need_new:
            if dirs:
                why.append(f"new の run が {len(dirs)} 本あるが new を比較しない指定")
            continue
        if len(dirs) != plan["reps"]:
            why.append(f"{label} {len(dirs)} 本 (予定 {plan['reps']} 本)")
    real = [os.path.realpath(d) for d in base + new]
    if len(set(real)) != len(real):
        why.append("同じ run が 2 回以上指定されている")
    tinfo = {}
    for d in base + new:
        w, ti = check_run(d, plan)
        lines.append(f"  {os.path.basename(os.path.normpath(d))}: {'OK' if not w else 'FAIL — ' + '; '.join(w)}")
        why += [f"{os.path.basename(os.path.normpath(d))}: {x}" for x in w]
        if not w:
            tinfo[d] = ti
    if tinfo:
        d0 = next(iter(tinfo))
        for d, ti in tinfo.items():
            if ti != tinfo[d0]:
                ks = sorted(k for k in set(ti) | set(tinfo[d0]) if ti.get(k) != tinfo[d0].get(k))
                why.append(f"{os.path.basename(d)}: 物理時刻・step の属性が {os.path.basename(d0)} と違う "
                           f"({', '.join(f'{k} {tinfo[d0].get(k)} vs {ti.get(k)}' for k in ks[:4])})")
        lines.append(f"  最終出力の物理時刻・step の属性 ({os.path.basename(d0)}): "
                     + (", ".join(f"{k}={v}" for k, v in sorted(tinfo[d0].items())) or "(属性なし)"))
    verdict = "PASS" if not why else "FAIL"
    lines += [f"  FAIL の理由: {x}" for x in why]
    lines.append(f"  >> COMPLETE: {verdict}" + (f" ({len(why)} 件)" if why else ""))
    return dict(verdict=verdict, reasons=why, lines=lines)


# ------------------------------------------------------------------ h5
def h5_items(path):
    """{データセット名: (shape, dtype, attrs)} と {グループ名: attrs}。"""
    ds, grp = {}, {}
    with h5py.File(path, "r") as f:
        grp["/"] = {k: _attr(v) for k, v in f.attrs.items()}

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                ds[name] = (obj.shape, str(obj.dtype), {k: _attr(v) for k, v in obj.attrs.items()})
            else:
                grp[name] = {k: _attr(v) for k, v in obj.attrs.items()}
        f.visititems(visit)
    return ds, grp


def _attr(v):
    if isinstance(v, bytes):
        return v.decode(errors="replace")
    if isinstance(v, np.ndarray):
        return tuple(v.ravel().tolist())
    if isinstance(v, np.generic):
        return v.item()
    return v


def h5_data(path, names):
    with h5py.File(path, "r") as f:
        return {n: f[n][()] for n in names}


def metric(a, b):
    """max|a−b| / max|a| (整数は不一致の個数を返さず、一致なら 0、不一致なら inf)。"""
    if a.shape != b.shape:
        return float("inf")
    if a.dtype.kind in "iub":
        return 0.0 if np.array_equal(a, b) else float("inf")
    a64 = np.asarray(a, dtype=np.float64)
    b64 = np.asarray(b, dtype=np.float64)
    if a64.size == 0:
        return 0.0
    d = np.nanmax(np.abs(a64 - b64)) if np.all(np.isfinite(a64)) and np.all(np.isfinite(b64)) else float("inf")
    s = np.max(np.abs(a64))
    if s == 0.0:
        return 0.0 if d == 0.0 else float("inf")
    return float(d / s)


def bit_equal(a, b):
    return a.shape == b.shape and a.dtype == b.dtype and a.tobytes() == b.tobytes()


# ------------------------------------------------------------------ CSV
def read_table(path):
    """数値表を {列名: ndarray} にする。区切りは ',' (probe は ' , ')。行キー (step,inner,phase) を別に返す。"""
    with open(path) as f:
        lines = [x.rstrip("\n") for x in f if x.strip()]
    if not lines:
        return {}, []
    # 区切りは ',' (残差・probe の ' , ' を含む)。',' の無い表 (conjugate_Tw_*.csv は空白区切り) は空白で切る
    if "," in lines[0]:
        hdr = [h.strip() for h in lines[0].split(",")]
        rows = [[c.strip() for c in ln.split(",")] for ln in lines[1:]]
    else:
        hdr = lines[0].split()
        rows = [ln.split() for ln in lines[1:]]
    keys = []
    cols = {h: [] for h in hdr}
    keycols = [i for i, h in enumerate(hdr) if h in ("step", "inner", "phase", "Step", "var", "physID")]
    for r in rows:
        keys.append(tuple(r[i] if i < len(r) else "" for i in keycols))
        for i, h in enumerate(hdr):
            cols[h].append(r[i] if i < len(r) else "")
    out = {}
    for h, v in cols.items():
        try:
            out[h] = np.array([float(x) for x in v])
        except ValueError:
            out[h] = np.array(v, dtype=object)
    return out, keys


def ulp_key(x):
    """float32 の値を単調な整数に写す (ulp 距離 = 整数の差)。CSV の値は float32 を double で書いたもの。"""
    i = int(np.array([float(x)], dtype=np.float32).view(np.int32)[0])
    return i if i >= 0 else -(i & 0x7FFFFFFF)


def ulp_dist(a, b):
    try:
        return abs(ulp_key(a) - ulp_key(b))
    except ValueError:
        return 0 if a == b else -1


def step0_rows(path):
    with open(path) as f:
        rd = csv.reader(f)
        hdr = next(rd)
        rows = list(rd)
    if not rows:
        return hdr, []
    s0 = rows[0][0]
    return hdr, [r for r in rows if r[0] == s0]


# ------------------------------------------------------------------ 比較本体
class Report:
    def __init__(self):
        self.lines = []
        self.verdicts = {}

    def p(self, s=""):
        self.lines.append(s)

    def v(self, key, val):
        self.verdicts[key] = val
        self.p(f"  >> {key}: {val}")


def compare(base, new, rep, cfgname, plan):
    """登録判定 A。plan (matrix_spec.plan_from の dict) の N を最終出力の step に使う (完全性検査を通った run だけを渡す)。"""
    allruns = base + new
    tag = {d: ("B" if d in base else "N") + str((base if d in base else new).index(d) + 1) for d in allruns}
    rep.p(f"=== 構成 {cfgname}")
    for d in allruns:
        rep.p(f"  {tag[d]}: {os.path.relpath(d, HERE)}")

    # ---- (c) NaN/Inf
    rep.p("\n[c] NaN/Inf (各 run)")
    nan_bad = 0
    for d in allruns:
        p = os.path.join(d, "NANCHECK.txt")
        v = open(p).read().strip().splitlines()[-1] if os.path.exists(p) else "NANCHECK: (無し)"
        if "PASS" not in v:
            nan_bad += 1
        rep.p(f"  {tag[d]}: {v}")
    rep.v("NAN", "PASS" if nan_bad == 0 else f"FAIL ({nan_bad} run)")

    # ---- 出力ファイル集合
    rep.p("\n[構造] 出力ファイルの集合")
    fsets = {d: out_files(d) for d in allruns}
    ref = fsets[allruns[0]]
    fs_bad = [d for d in allruns if fsets[d] != ref]
    for d in fs_bad:
        rep.p(f"  {tag[d]}: 余分 {sorted(fsets[d] - ref)} / 欠け {sorted(ref - fsets[d])}")
    rep.p(f"  {len(ref)} ファイル: {', '.join(sorted(ref))}")
    common = set.intersection(*fsets.values())
    h5s = sorted(f for f in common if f.endswith(".h5"))

    # ---- h5 の構造 (集合・shape・dtype・属性)
    rep.p("\n[構造] h5 のデータセット集合・shape・dtype・属性")
    st_bad = len(fs_bad)
    items = {}
    for fn in h5s:
        items[fn] = {d: h5_items(os.path.join(d, fn)) for d in allruns}
        r_ds, r_grp = items[fn][allruns[0]]
        for d in allruns[1:]:
            ds, grp = items[fn][d]
            if set(ds) != set(r_ds):
                st_bad += 1
                rep.p(f"  {fn} {tag[d]}: データセット集合が違う 余分 {sorted(set(ds) - set(r_ds))[:8]} 欠け {sorted(set(r_ds) - set(ds))[:8]}")
            for n in set(ds) & set(r_ds):
                if ds[n][:2] != r_ds[n][:2]:
                    st_bad += 1
                    rep.p(f"  {fn}:{n} {tag[d]}: shape/dtype {ds[n][:2]} vs {r_ds[n][:2]}")
                if ds[n][2] != r_ds[n][2]:
                    st_bad += 1
                    rep.p(f"  {fn}:{n} {tag[d]}: 属性 {ds[n][2]} vs {r_ds[n][2]}")
            for g in set(grp) | set(r_grp):
                if grp.get(g) != r_grp.get(g):
                    st_bad += 1
                    rep.p(f"  {fn} グループ {g} {tag[d]}: 属性 {grp.get(g)} vs {r_grp.get(g)}")
        rep.p(f"  {fn}: {len(r_ds)} データセット")
    rep.v("STRUCT", "PASS" if st_bad == 0 else f"FAIL ({st_bad} 件)")

    # ---- ログの行
    rep.p("\n[ログ] 名前登録・復元・初期化・警告の行")
    logs = {}
    for d in allruns:
        p = os.path.join(d, "forge_run.log") if os.path.exists(os.path.join(d, "forge_run.log")) else os.path.join(d, "convert.log")
        seen = []
        if os.path.exists(p):
            with open(p, errors="replace") as f:
                for ln in f:
                    if LOG_PAT.search(ln) and not LOG_DROP.search(ln) and not LOG_NEW_INFO.search(ln):
                        x = ln.rstrip()
                        if x not in seen:
                            seen.append(x)
        logs[d] = seen
    lg_bad = 0
    # 行の**集合**で比べる (順序は見ない): 変更後は環境変数の出力登録を確保の前へ移すので、
    # [FORGE_OUT_RESIDUALS] 等の行の位置が変わるのは設計どおり (監査 §2 の (i))
    for d in allruns[1:]:
        if set(logs[d]) != set(logs[allruns[0]]):
            lg_bad += 1
            a, b = set(logs[allruns[0]]), set(logs[d])
            rep.p(f"  {tag[d]}: B1 に無い行 {sorted(b - a)[:5]} / B1 にだけある行 {sorted(a - b)[:5]}")
    for x in logs[allruns[0]][:40]:
        rep.p(f"    | {x[:200]}")
    rep.v("LOG", "PASS" if lg_bad == 0 else f"FAIL ({lg_bad} run)")

    # ---- (a) step 0 の残差行
    rep.p("\n[a] step 0 の残差行 (全列ビット一致; base で割れた列は new が base の観測値のどれか)")
    rc = [d for d in allruns if os.path.exists(os.path.join(d, "residual_history.csv"))]
    if len(rc) == len(allruns) and rc:
        hdr0, _ = step0_rows(os.path.join(rc[0], "residual_history.csv"))
        rows = {d: step0_rows(os.path.join(d, "residual_history.csv"))[1] for d in allruns}
        nrows = {len(v) for v in rows.values()}
        a_bad = 0
        split_cols = []
        ulp_base = 0     # 情報: base 反復の間の最大 ulp 幅
        ulp_new = 0      # 情報: new の値と最近傍の base 値の最大 ulp 距離
        if len(nrows) != 1:
            a_bad += 1
            rep.p(f"  step 0 の行数が違う: {[len(rows[d]) for d in allruns]}")
        else:
            for i in range(nrows.pop()):
                for j, h in enumerate(hdr0):
                    bvals = {rows[d][i][j] for d in base}
                    if len(bvals) > 1:
                        split_cols.append((i, h, len(bvals)))
                        bl = sorted(bvals)
                        ulp_base = max(ulp_base, max(ulp_dist(x, y) for x in bl for y in bl))
                    for d in new:
                        if rows[d][i][j] not in bvals:
                            a_bad += 1
                            u = min(ulp_dist(rows[d][i][j], x) for x in bvals)
                            ulp_new = max(ulp_new, u)
                            if a_bad <= 30:
                                rep.p(f"  行 {i} 列 {h}: {tag[d]} {rows[d][i][j]} ∉ base {sorted(bvals)} (最近傍まで {u} ulp)")
        rep.p(f"  base で値が割れた (行, 列, 値の数): {split_cols[:12] if split_cols else 'なし (base 3 回でビット一致)'}"
              f"{' …' if len(split_cols) > 12 else ''}")
        rep.p(f"  (情報) step 0 の行数 {len(rows[allruns[0]])}、base 内の最大 ulp 幅 {ulp_base}"
              + (f"、new の値と最近傍 base 値の最大 ulp 距離 {ulp_new} (不一致 {a_bad} 値)" if new else ""))
        rep.ulp = (ulp_base, ulp_new)
        if new:
            rep.verdicts["STEP0_ULP"] = f"base 幅 {ulp_base} / new 最近傍 {ulp_new} ulp"
        rep.v("STEP0", ("PASS" if a_bad == 0 else f"FAIL ({a_bad} 値)") if new else
              f"BASE-ONLY (割れた列 {len(split_cols)})")
    else:
        rep.p("  residual_history.csv が無い (変換器など)")

    # ---- 初期出力 (保存量・原始量・幾何量のビット一致)。初期出力は res_0、最終出力は予定 N (存在する最大 step ではない)
    steps = [0, plan["N"]] if plan["kind"] == "forge" else []
    if steps:
        fn0 = f"res_{steps[0]}.h5"
        rep.p(f"\n[a] 初期出力 {fn0}: 保存量・原始量・幾何量のビット一致 (他の量は情報として不一致の数だけ)")
        ds0 = items.get(fn0, {}).get(allruns[0], ({}, {}))[0]
        cmp_names = [n for n in ds0 if (not n.startswith("VALUE/")) or T0_VALUE.match(n.split("/", 1)[1])]
        info_names = [n for n in ds0 if n not in cmp_names]
        ref = h5_data(os.path.join(allruns[0], fn0), cmp_names + info_names)
        t0_bad, info_diff = [], {}
        for d in allruns[1:]:
            x = h5_data(os.path.join(d, fn0), cmp_names + info_names)
            for n in cmp_names:
                if not bit_equal(ref[n], x[n]):
                    t0_bad.append((tag[d], n, metric(ref[n], x[n])))
            for n in info_names:
                if not bit_equal(ref[n], x[n]):
                    info_diff.setdefault(n, []).append(tag[d])
        rep.p(f"  比べた量 {len(cmp_names)}、比べない量 {len(info_names)}")
        for t, n, m in t0_bad[:20]:
            rep.p(f"  不一致 {t} {n}: m={m:.3e}")
        if info_diff:
            rep.p(f"  (情報) 比べない量で B1 と違うもの: {sorted(info_diff)[:20]}")
        rep.v("INIT_OUT", ("PASS" if not t0_bad else f"FAIL ({len(t0_bad)} 件)") if cmp_names else "FAIL (比べた量 0)")

    # ---- (b) 最後の出力・境界出力・CSV・残差履歴
    rep.p(f"\n[b] N step 後 (予定 N = {plan.get('N')}): m = max|A−B|/max|A|。S = 同ビルド内ペアの最大、D = base×new の最大、合格 D ≤ 2·S")
    finals = []
    if steps:
        last = steps[-1]
        finals = [f for f in h5s if re.search(rf"_{last}\.h5$", f) and f != f"res_{steps[0]}.h5"]
    if not steps:
        finals = h5s   # 変換器
    tables = sorted(f for f in common if f.endswith((".csv", ".out")))
    pairs_b = [(base[i], base[j]) for i in range(len(base)) for j in range(i + 1, len(base))]
    pairs_n = [(new[i], new[j]) for i in range(len(new)) for j in range(i + 1, len(new))]
    pairs_x = [(b, n) for b in base for n in new]
    rows_out = []   # (ファイル, 量, S_base, S_new, D, verdict)
    absinfo = {}    # FAIL の読み解き用 (判定には使わない): (ファイル, 量) -> (max|B1|, max|N1|, D の組の max|A−B|)

    def eval_set(getter, names, fname):
        cache = {d: getter(d) for d in allruns}
        for n in names:
            def m(p):
                a, b = cache[p[0]].get(n), cache[p[1]].get(n)
                if a is None or b is None:
                    return float("inf")
                return metric(a, b)
            sb = max((m(p) for p in pairs_b), default=float("nan"))
            sn = max((m(p) for p in pairs_n), default=float("nan"))
            dd = max((m(p) for p in pairs_x), default=float("nan"))
            if new:
                s = max(sb, sn)
                ok = (dd == 0.0) if s == 0.0 else (dd <= 2.0 * s)
                vd = "PASS" if ok else "FAIL"
                if not ok:
                    def amax(x):
                        return float(np.max(np.abs(np.asarray(x, dtype=np.float64)))) if x is not None and np.size(x) else float("nan")
                    worst = max(pairs_x, key=m)
                    a, b = cache[worst[0]].get(n), cache[worst[1]].get(n)
                    dabs = float(np.max(np.abs(np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64)))) \
                        if a is not None and b is not None and np.shape(a) == np.shape(b) else float("nan")
                    absinfo[(fname, n)] = (amax(cache[base[0]].get(n)), amax(cache[new[0]].get(n)), dabs,
                                           f"{tag[worst[0]]}-{tag[worst[1]]}")
            else:
                vd = "-"
            rows_out.append((fname, n, sb, sn, dd, vd))

    for fn in finals:
        names = sorted(items[fn][allruns[0]][0])
        eval_set(lambda d, fn=fn, names=names: h5_data(os.path.join(d, fn), names), names, fn)
    for fn in tables:
        def get(d, fn=fn):
            t, keys = read_table(os.path.join(d, fn))
            return {k: v for k, v in t.items() if v.dtype != object and k not in ("step", "inner", "Step")}
        t0, _ = read_table(os.path.join(allruns[0], fn))
        names = [k for k, v in t0.items() if v.dtype != object and k not in ("step", "inner", "Step")]
        eval_set(get, names, fn)

    b_fail = [r for r in rows_out if r[5] == "FAIL"]
    rep.p(f"  {'ファイル':28s} {'量':28s} {'S_base':>10s} {'S_new':>10s} {'D':>10s}  判定")
    # 主要なものから: 最終場の保存量、残差履歴、その他。S_base の大きい順に上位を出し、FAIL は全部出す
    def fmt(x):
        return "      -   " if x != x else f"{x:10.3e}"
    shown = 0
    for r in sorted(rows_out, key=lambda r: (r[5] != "FAIL", -(r[2] if r[2] == r[2] else 0))):
        if r[5] == "FAIL" or shown < 25:
            rep.p(f"  {r[0][:28]:28s} {r[1][:28]:28s} {fmt(r[2])} {fmt(r[3])} {fmt(r[4])}  {r[5]}")
            shown += 1
    for (fn_, n_), (ab, an, dabs, pr) in absinfo.items():
        rep.p(f"  (情報) FAIL {fn_}:{n_}: max|B1| {ab:.3e}, max|N1| {an:.3e}, D の組 {pr} の max|A−B| {dabs:.3e}")
    nz = sum(1 for r in rows_out if r[2] == r[2] and r[2] > 0)
    rep.p(f"  ({len(rows_out)} 量。base 内で差が 0 でない量 {nz})")
    key = {}
    for r in rows_out:
        if re.match(r"res_\d+\.h5$", r[0]) and r[1] in ("VALUE/ro", "VALUE/roUx", "VALUE/roe", "VALUE/P", "VALUE/T"):
            key[r[1]] = r[2]
        if r[0] == "residual_history.csv" and r[1] in ("rms_ro", "rms_roe"):
            key["csv:" + r[1]] = r[2]
    rep.p("  主要量の S_base: " + ", ".join(f"{k}={v:.2e}" for k, v in key.items()))
    if not rows_out:
        rep.v("NSTEP", "FAIL (比較した量 0)")     # 量 0・FAIL 0 を合格にしない
    elif new:
        rep.v("NSTEP", "PASS" if not b_fail else f"FAIL ({len(b_fail)} 量)")
        # 情報: D/(2S) が最大のデータセット (1 を超えると FAIL)。S = 0 で D = 0 の量は除く
        cand = [(r[4] / (2 * max(r[2], r[3])), r) for r in rows_out
                if r[4] == r[4] and max(r[2], r[3]) > 0]
        if cand:
            q, r = max(cand, key=lambda x: x[0])
            rep.verdicts["WORST"] = f"{r[0]}:{r[1]} D/2S={q:.2f} (D {r[4]:.2e}, S {max(r[2], r[3]):.2e})"
    else:
        smax = max((r[2] for r in rows_out if r[2] == r[2]), default=float("nan"))
        rep.v("NSTEP", f"BASE-ONLY (S_base 最大 {smax:.3e})")
    rep.rows = rows_out
    return rep


# ================================================================== 追加診断 B (plan §6.2、2026-10-07 事前登録)
# 尺度だけを替える: d(A,B) = max|A − B| (float64、絶対 L∞)。ペアに対称で、反復の並べ替えに不変。
# S_abs = 同ビルド内 6 対 (base 3 対 + new 3 対) の最大、D_abs = ビルド間 9 対の最大。合格は D_abs ≤ 2·S_abs (S_abs = 0 なら D_abs = 0)。
# 全入力が有限で shape・列・行キーが対応すること (比較不能は FAIL)。整数・文字列は全 run で厳密一致。
# 登録判定 A (上の metric / compare) は書き換えない。
TABLE_KEYS = ("step", "inner", "phase", "Step", "var", "physID")


def dabs(a, b):
    """d(A,B) = max|A − B| を float64 で。shape が違う・非有限を含むときは None (比較不能)。"""
    a64 = np.asarray(a, dtype=np.float64)
    b64 = np.asarray(b, dtype=np.float64)
    if a64.shape != b64.shape:
        return None
    if not (np.all(np.isfinite(a64)) and np.all(np.isfinite(b64))):
        return None
    if a64.size == 0:
        return 0.0
    return float(np.max(np.abs(a64 - b64)))


def judge_abs(base_arrs, new_arrs, dfun=None):
    """B の判定。base_arrs / new_arrs は run ごとの配列 (欠落は None)。
    返り値 dict(S, D, verdict, reason, worst_S, worst_D) — worst_* は (base/new の列での位置の組)。"""
    nan = float("nan")
    dfun = dfun or dabs
    allv = list(base_arrs) + list(new_arrs)
    nb = len(base_arrs)
    out = dict(S=nan, D=nan, verdict="FAIL", reason="", worst_S=None, worst_D=None)
    if nb < 2 or len(new_arrs) < 1:
        out["reason"] = "比較不能: 反復が足りない"
        return out
    if any(x is None for x in allv):
        out["reason"] = "比較不能: 欠落した run がある"
        return out
    arrs = [np.asarray(x) for x in allv]
    if len({a.shape for a in arrs}) != 1:
        out["reason"] = "比較不能: shape が違う"
        return out
    if any(a.dtype.kind != "f" for a in arrs):
        ok = all(a.dtype == arrs[0].dtype and np.array_equal(a, arrs[0]) for a in arrs[1:])
        out.update(S=0.0 if ok else nan, D=0.0 if ok else nan, verdict="PASS" if ok else "FAIL",
                   reason="整数・文字列は厳密一致" + ("" if ok else ": 不一致"))
        return out
    if any(not np.all(np.isfinite(a)) for a in arrs):
        out["reason"] = "比較不能: 非有限を含む"
        return out
    idx = list(range(len(arrs)))
    within = [(i, j) for g in (idx[:nb], idx[nb:]) for k, i in enumerate(g) for j in g[k + 1:]]
    cross = [(i, j) for i in idx[:nb] for j in idx[nb:]]
    dv = {}
    for p in within + cross:
        v = dfun(arrs[p[0]], arrs[p[1]])
        if v is None:
            out["reason"] = "比較不能: d が定義できない"
            return out
        dv[p] = v
    ws = max(within, key=lambda p: dv[p]) if within else None
    wd = max(cross, key=lambda p: dv[p])
    S = dv[ws] if ws else nan
    D = dv[wd]
    ok = (D == 0.0) if S == 0.0 else (D <= 2.0 * S)
    out.update(S=S, D=D, verdict="PASS" if ok else "FAIL", reason="" if ok else ("S_abs = 0 で D_abs > 0" if S == 0.0 else "D_abs > 2·S_abs"),
               worst_S=ws, worst_D=wd)
    return out


def self_check_abs(base_arrs, new_arrs, ref):
    """自己検査: (i) base・new それぞれの全順列 (3!×3! = 36 通り) と (ii) base/new の交換で S・D・判定が ref と同じか。
    d は (配列, 配列) の**順序つき**の組でキャッシュする (逆順の組は別に計算するので、対称でなければ検出される)。"""
    cache = {}

    def dfun(a, b):
        k = (id(a), id(b))
        if k not in cache:
            cache[k] = dabs(a, b)
        return cache[k]

    def same(r):
        def eq(x, y):
            return (x != x and y != y) or x == y
        return eq(r["S"], ref["S"]) and eq(r["D"], ref["D"]) and r["verdict"] == ref["verdict"]

    arrs_b = [None if x is None else np.asarray(x) for x in base_arrs]
    arrs_n = [None if x is None else np.asarray(x) for x in new_arrs]
    n_perm = n_bad = 0
    for pb in itertools.permutations(range(len(arrs_b))):
        for pn in itertools.permutations(range(len(arrs_n))):
            r = judge_abs([arrs_b[i] for i in pb], [arrs_n[i] for i in pn], dfun)
            n_perm += 1
            n_bad += 0 if same(r) else 1
    sw = judge_abs(arrs_n, arrs_b, dfun)        # (ii) base と new を入れ替える (S は同じ集合、D は逆順の組)
    swap_ok = same(sw)
    return n_perm, n_bad, swap_ok


def final_files(d, n=None):
    """B の比較対象 (A と同じ保存時点): 最終ステップの h5 (初期出力を除く) と CSV/probe 出力。変換器は全 h5。
    n = 予定 N (完全性検査と同じ値)。n を渡さないと「存在する最大 step」を最終と見なす旧い挙動 (fixedwidth の凍結の再現用だけ)。"""
    fs = out_files(d)
    steps = [0, n] if n is not None else res_steps(d)
    h5s = sorted(f for f in fs if f.endswith(".h5"))
    if steps:
        last = steps[-1]
        finals = [f for f in h5s if re.search(rf"_{last}\.h5$", f)]
        if len(steps) > 1:
            finals = [f for f in finals if f != f"res_{steps[0]}.h5"]
    else:
        finals = h5s
    tables = sorted(f for f in fs if f.endswith((".csv", ".out")))
    return finals, tables


def compare_abs(base, new, cfgname, plan):
    """追加診断 B。返り値 (報告の行, 量ごとの行の list[dict], 自己検査の集計 dict)。
    plan の N を最終出力の step に使う (完全性検査を通った run だけを渡す)。比較した量が 0 なら FAIL の行を 1 つ足す。"""
    allruns = base + new
    tag = {d: ("B" if d in base else "N") + str((base if d in base else new).index(d) + 1) for d in allruns}
    lines = [f"=== 構成 {cfgname} (追加診断 B: d = max|A−B|、S_abs = 同ビルド内 6 対の最大、D_abs = ビルド間 9 対の最大、合格 D_abs ≤ 2·S_abs)"]
    for d in allruns:
        lines.append(f"  {tag[d]}: {d}")
    rows = []
    sc = dict(quantities=0, orderings=0, perm_bad=0, swap_bad=0)

    def add(fname, name, barrs, narrs, struct_bad=""):
        r = judge_abs(barrs, narrs)
        if struct_bad:
            r.update(verdict="FAIL", reason=("構造: " + struct_bad + ("; " + r["reason"] if r["reason"] else "")))
        n_perm, n_bad, swap_ok = self_check_abs(barrs, narrs, judge_abs(barrs, narrs))
        sc["quantities"] += 1
        sc["orderings"] += n_perm
        sc["perm_bad"] += n_bad
        sc["swap_bad"] += 0 if swap_ok else 1
        S, D = r["S"], r["D"]
        ratio = (D / (2 * S)) if (S == S and S > 0) else (0.0 if (D == 0.0) else float("inf"))
        det = ""
        if r["verdict"] == "FAIL":
            mx = []
            for d, x in zip(allruns, list(barrs) + list(narrs)):
                if x is None:
                    mx.append(f"{tag[d]} 欠落")
                elif np.asarray(x).dtype.kind == "f":
                    a = np.asarray(x, dtype=np.float64)
                    mx.append(f"{tag[d]} max|x| {np.max(np.abs(a)) if a.size else 0:.3e}"
                              + ("" if np.all(np.isfinite(a)) else " (非有限あり)"))
                else:
                    mx.append(f"{tag[d]} (整数・文字列)")
            if r["worst_D"]:
                i, j = r["worst_D"]
                det += f"最悪のビルド間ペア {tag[allruns[i]]}-{tag[allruns[j]]} d={D:.3e}; "
            if r["worst_S"]:
                i, j = r["worst_S"]
                det += f"最悪の同ビルド内ペア {tag[allruns[i]]}-{tag[allruns[j]]} d={S:.3e}; "
            det += "; ".join(mx)
        rows.append(dict(cfg=cfgname, file=fname, name=name, S=S, D=D, ratio=ratio, verdict=r["verdict"],
                         reason=r["reason"], detail=det))

    # ---- ファイル集合
    sets = {d: final_files(d, plan["N"]) for d in allruns}
    f_ref = sets[allruns[0]]
    for d in allruns[1:]:
        if sets[d] != f_ref:
            rows.append(dict(cfg=cfgname, file="(ファイル集合)", name=tag[d], S=float("nan"), D=float("nan"),
                             ratio=float("inf"), verdict="FAIL", reason="比較不能: 出力ファイルの集合が違う",
                             detail=f"{sets[d]} vs {f_ref}"))
    h5f = sorted(set.intersection(*[set(sets[d][0]) for d in allruns]))
    tbf = sorted(set.intersection(*[set(sets[d][1]) for d in allruns]))
    # ---- h5
    for fn in h5f:
        its = {d: h5_items(os.path.join(d, fn)) for d in allruns}
        names = sorted(set().union(*[set(its[d][0]) for d in allruns]))
        grp_bad = "" if all(its[d][1] == its[allruns[0]][1] for d in allruns) else "グループ属性が違う"
        if grp_bad:
            rows.append(dict(cfg=cfgname, file=fn, name="(グループ属性)", S=float("nan"), D=float("nan"),
                             ratio=float("inf"), verdict="FAIL", reason="構造: " + grp_bad, detail=""))
        for n in names:
            meta = [its[d][0].get(n) for d in allruns]
            sb = ""
            if any(m is None for m in meta):
                sb = "データセットが無い run がある"
            elif any(m[:2] != meta[0][:2] for m in meta):
                sb = "shape/dtype が違う"
            elif any(m[2] != meta[0][2] for m in meta):
                sb = "属性が違う"
            data = []
            for d in allruns:
                data.append(h5_data(os.path.join(d, fn), [n])[n] if its[d][0].get(n) is not None else None)
            add(fn, n, data[:len(base)], data[len(base):], sb)
    # ---- 表 (CSV・probe)
    for fn in tbf:
        tabs = {d: read_table(os.path.join(d, fn)) for d in allruns}
        cols0 = list(tabs[allruns[0]][0].keys())
        keys0 = tabs[allruns[0]][1]
        bad = ""
        if any(list(tabs[d][0].keys()) != cols0 for d in allruns):
            bad = "列が違う"
        elif any(tabs[d][1] != keys0 for d in allruns):
            bad = "行キー (step・inner・phase 等) が対応しない"
        names = sorted(set().union(*[set(tabs[d][0].keys()) for d in allruns]) - set(TABLE_KEYS))
        for n in names:
            data = [tabs[d][0].get(n) for d in allruns]
            add(fn, n, data[:len(base)], data[len(base):], ("比較不能: " + bad) if bad else "")

    if sc["quantities"] == 0:   # 量 0・FAIL 0 を合格にしない
        rows.append(dict(cfg=cfgname, file="(比較量)", name="-", S=float("nan"), D=float("nan"), ratio=float("inf"),
                         verdict="FAIL", reason="比較した量が 0", detail=""))
    nfail = sum(1 for r in rows if r["verdict"] == "FAIL")
    lines.append(f"\n  量 {len(rows)}、FAIL {nfail}")
    lines.append(f"  自己検査: {sc['quantities']} 量 × 並べ替え {sc['orderings'] // max(sc['quantities'], 1)} 通り = {sc['orderings']} 評価で"
                 f"判定・S_abs・D_abs が変わったもの {sc['perm_bad']}、base/new 交換で変わったもの {sc['swap_bad']}")
    lines.append(f"\n  {'ファイル':28s} {'量':30s} {'S_abs':>11s} {'D_abs':>11s} {'D/2S':>8s}  判定  理由")
    for r in sorted(rows, key=lambda r: (r["verdict"] != "FAIL", -(r["ratio"] if r["ratio"] == r["ratio"] else 1e300))):
        def f(x):
            return "        -  " if x != x else f"{x:11.3e}"
        lines.append(f"  {r['file'][:28]:28s} {r['name'][:30]:30s} {f(r['S'])} {f(r['D'])} {r['ratio']:8.3f}  {r['verdict']:4s}  {r['reason']}")
        if r["detail"]:
            lines.append(f"      {r['detail']}")
    return lines, rows, sc


def diff2(a, b):
    da, _ = h5_items(a)
    db, _ = h5_items(b)
    names = sorted(set(da) & set(db))
    xa, xb = h5_data(a, names), h5_data(b, names)
    print(f"A={a}\nB={b}\n集合の差: A のみ {sorted(set(da) - set(db))[:10]} / B のみ {sorted(set(db) - set(da))[:10]}")
    for n in names:
        m = metric(xa[n], xb[n])
        if m != 0.0:
            print(f"  {n:40s} m={m:.3e}")
    print(f"  ({len(names)} 共通データセット、差 0 でないものだけ表示)")


def build_jobs(a, need_new):
    """比較の単位 (構成) ごとに予定と run を決める。返り値は dict(cfg, plan, base, new, why, notes) の list。
    --all / --cfg は構成表 (matrix_spec.planned) の予定と registry の rep 1..reps。--base/--new の直接指定は --steps 等の予定。"""
    jobs = []
    direct_opts = [o for o, v in (("--steps", a.steps), ("--out-interval", a.out_interval), ("--require", a.require),
                                  ("--kind", a.kind), ("--reps", a.reps)) if v]
    root = os.path.abspath(os.path.expanduser(a.root)) if a.root else HERE
    if a.all or a.cfg:
        if direct_opts:
            sys.exit(f"{', '.join(direct_opts)} は --base/--new の直接指定のときだけ使う (構成表の予定を上書きしない)")
        for c in (list(ms.CONFIGS) if a.all else [a.cfg]):
            why, notes = [], []
            if c not in ms.CONFIGS:
                plan = ms.plan_from(c)
                why.append(f"構成 {c} が構成表 (matrix_spec.CONFIGS) に無い")
            else:
                plan = ms.planned(c)
            base, wb, nb = select_runs(c, a.base_build, plan, root)
            why += wb
            notes += nb
            new = []
            if need_new:
                if not a.new_build:
                    why.append("new のビルド名 (--new-build) が無い")
                else:
                    new, wn, nn = select_runs(c, a.new_build, plan, root)
                    why += wn
                    notes += nn
            jobs.append(dict(cfg=c, plan=plan, base=base, new=new, why=why, notes=notes))
    else:
        plan = ms.plan_from(a.tag, a.kind or "forge", a.reps or ms.REPS, a.steps, a.out_interval, a.require or [])
        why = [] if a.base else ["--base の run が指定されていない"]
        jobs.append(dict(cfg=a.tag, plan=plan, base=[os.path.abspath(os.path.expanduser(x)) for x in a.base],
                         new=[os.path.abspath(os.path.expanduser(x)) for x in a.new], why=why, notes=[]))
    return jobs


def write_completeness_tsv(path, results):
    with open(path, "w") as f:
        f.write("cfg\tverdict\tbase\tnew\treps\tN\treasons\n")
        for c, comp, nb, nn, plan in results:
            f.write(f"{c}\t{comp['verdict']}\t{nb}\t{nn}\t{plan['reps']}\t{plan.get('N')}\t{' / '.join(comp['reasons'])}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cfg")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--base", nargs="*", default=[])
    ap.add_argument("--new", nargs="*", default=[])
    ap.add_argument("--base-build", default="base")
    ap.add_argument("--new-build")
    ap.add_argument("--out")
    ap.add_argument("--out-dir")
    ap.add_argument("--diff2", nargs=2)
    ap.add_argument("--metric", choices=["m", "abs"], default="m",
                    help="m = 登録判定 A (max|A−B|/max|A|)、abs = 追加診断 B (max|A−B|、plan §6.2)")
    ap.add_argument("--tag", default="(指定)", help="--base/--new 直接指定のときの構成名")
    ap.add_argument("--root", help="--all/--cfg: registry.tsv と run_* の置き場 (既定はこのスクリプトのディレクトリ)")
    # 直接指定のときの予定 (完全性検査)。--all/--cfg では構成表 (matrix_spec.planned) を使い、これらは受け付けない
    ap.add_argument("--steps", type=int, help="直接指定: 予定 N (最終出力 res_N.h5)。無ければ完全性 FAIL")
    ap.add_argument("--out-interval", type=int, help="直接指定: 出力間隔 (既定 N)")
    ap.add_argument("--reps", type=int, help=f"直接指定: 予定の反復数 (base・new とも、既定 {ms.REPS})")
    ap.add_argument("--kind", choices=["forge", "convert"], help="直接指定: forge (既定) / convert")
    ap.add_argument("--require", nargs="*", default=[], help="直接指定: 追加の必要ファイル ({N} を展開、*?[ は glob)")
    a = ap.parse_args()
    if a.diff2:
        diff2(*a.diff2)
        return
    if a.metric == "abs":
        main_abs(a)
        return
    need_new = bool(a.new_build) if (a.all or a.cfg) else bool(a.new)
    jobs = build_jobs(a, need_new)
    summary, comps = [], []
    for j in jobs:
        c, plan, base, new = j["cfg"], j["plan"], j["base"], j["new"]
        comp = completeness(plan, base, new, j["why"], need_new)
        comps.append((c, comp, len(base), len(new), plan))
        head = [f"=== 構成 {c} — 完全性検査 (比較の前。FAIL なら値の比較はしない)"] + [f"  (情報) {x}" for x in j["notes"]] + comp["lines"]
        if comp["verdict"] != "PASS":
            txt = "\n".join(head + ["  (完全性 FAIL のため値の比較はしない)"]) + "\n"
            r0 = comp["reasons"]
            summary.append(f"{c:22s} base {len(base)} new {len(new)} | COMPLETE FAIL ({len(r0)} 件: {r0[0]}{' …' if len(r0) > 1 else ''})"
                           " | 比較しない")
        else:
            rep = Report()
            rep.lines += head + [""]
            rep.verdicts["COMPLETE"] = "PASS"
            rep = compare(base, new, rep, c, plan)
            txt = "\n".join(rep.lines) + "\n"
            summary.append(f"{c:22s} base {len(base)} new {len(new)} | " +
                           " | ".join(f"{k} {v}" for k, v in rep.verdicts.items()))
        if a.out_dir:
            os.makedirs(a.out_dir, exist_ok=True)
            open(os.path.join(a.out_dir, f"{c}.txt"), "w").write(txt)
        elif a.out:
            open(a.out, "w").write(txt)
        else:
            print(txt)
    nfc = sum(1 for x in comps if x[1]["verdict"] != "PASS")
    summary.append(f"\n完全性検査 (plan §5.1 #6): PASS {len(comps) - nfc} 構成、FAIL {nfc} 構成"
                   + (" — FAIL の構成は値を比較していない (構成ごとの報告に理由)" if nfc else ""))
    s = "\n".join(summary) + "\n"
    if a.out_dir:
        open(os.path.join(a.out_dir, "summary.txt"), "w").write(s)
        write_completeness_tsv(os.path.join(a.out_dir, "completeness.tsv"), comps)
    print(s)


def main_abs(a):
    """追加診断 B を回し、構成ごとの報告・全量の表 (all_quantities.tsv)・完全性の表 (completeness.tsv)・summary.txt・
    selfcheck.txt を out-dir に書く。完全性 FAIL の構成は値を比較せず、all_quantities.tsv に「(完全性)」の FAIL 行を 1 つ置く。"""
    t0 = time.time()
    jobs = build_jobs(a, True)
    od = a.out_dir or "."
    os.makedirs(od, exist_ok=True)
    allrows, summ, scl, comps = [], [], [], []
    nq_all = 0
    for j in jobs:
        c, plan, base, new = j["cfg"], j["plan"], j["base"], j["new"]
        t1 = time.time()
        comp = completeness(plan, base, new, j["why"], True)
        comps.append((c, comp, len(base), len(new), plan))
        head = [f"  (情報) {x}" for x in j["notes"]] + comp["lines"]
        if comp["verdict"] != "PASS":
            lines = [f"=== 構成 {c} (追加診断 B)"] + head + ["  (完全性 FAIL のため値の比較はしない)"]
            rows = [dict(cfg=c, file="(完全性)", name="-", S=float("nan"), D=float("nan"), ratio=float("inf"), verdict="FAIL",
                         reason="完全性: " + " / ".join(comp["reasons"]), detail="")]
            sc = dict(quantities=0, orderings=0, perm_bad=0, swap_bad=0)
        else:
            lines, rows, sc = compare_abs(base, new, c, plan)
            lines = lines[:1] + head + lines[1:]
        open(os.path.join(od, f"{c}.txt"), "w").write("\n".join(lines) + "\n")
        allrows += rows
        nq_all += sc["quantities"]
        nf = sum(1 for r in rows if r["verdict"] == "FAIL")
        fin = [r for r in rows if r["ratio"] == r["ratio"] and r["file"] not in ("(完全性)", "(比較量)")]
        worst = max(fin, key=lambda r: r["ratio"]) if fin else None
        ws = f"{worst['file']}:{worst['name']} D/2S={worst['ratio']:.3f}" if worst else "-"
        cv = comp["verdict"] + ("" if comp["verdict"] == "PASS" else f" ({len(comp['reasons'])} 件: {comp['reasons'][0]})")
        summ.append(f"{c:30s} base {len(base)} new {len(new)} | 完全性 {cv} | 量 {sc['quantities']:4d} | FAIL {nf:3d} | 最大 {ws}")
        scl.append(f"{c:30s} 量 {sc['quantities']:4d} | 並べ替え {sc['orderings']:6d} 評価で変化 {sc['perm_bad']} | "
                   f"base/new 交換で変化 {sc['swap_bad']} | {time.time() - t1:.1f} s")
        print(summ[-1], flush=True)
    with open(os.path.join(od, "all_quantities.tsv"), "w") as f:
        f.write("cfg\tfile\tname\tS_abs\tD_abs\tD_over_2S\tverdict\treason\tdetail\n")
        for r in allrows:
            f.write(f"{r['cfg']}\t{r['file']}\t{r['name']}\t{r['S']:.6e}\t{r['D']:.6e}\t{r['ratio']:.6g}\t{r['verdict']}\t{r['reason']}\t{r['detail']}\n")
    write_completeness_tsv(os.path.join(od, "completeness.tsv"), comps)
    fails = [r for r in allrows if r["verdict"] == "FAIL"]
    nfc = sum(1 for x in comps if x[1]["verdict"] != "PASS")
    top = sorted([r for r in allrows if r["ratio"] == r["ratio"] and r["verdict"] == "PASS"], key=lambda r: -r["ratio"])[:10]
    tail = [f"\n完全性検査 (plan §5.1 #6): PASS {len(comps) - nfc} 構成、FAIL {nfc} 構成 (completeness.tsv)",
            f"全体: {nq_all} 量、FAIL {len(fails)} (完全性 FAIL の構成は値を比較せず 1 行の FAIL として数える)、所要 {time.time() - t0:.0f} s"]
    for r in fails:
        tail.append(f"  FAIL {r['cfg']} {r['file']}:{r['name']} S_abs {r['S']:.3e} D_abs {r['D']:.3e} D/2S {r['ratio']:.3f} ({r['reason']}) {r['detail']}")
    tail.append("D/2S の上位 10 (PASS の量):")
    for r in top:
        tail.append(f"  {r['cfg']} {r['file']}:{r['name']} D/2S {r['ratio']:.3f} (S_abs {r['S']:.3e}, D_abs {r['D']:.3e})")
    open(os.path.join(od, "summary.txt"), "w").write("\n".join(summ + tail) + "\n")
    open(os.path.join(od, "selfcheck.txt"), "w").write(
        "比較器の自己検査 (追加診断 B): 各量で base・new の全順列 (3!×3!) と base/new の交換で S_abs・D_abs・判定が変わらないこと\n"
        + "\n".join(scl) + "\n")
    print("\n".join(tail))


if __name__ == "__main__":
    main()
