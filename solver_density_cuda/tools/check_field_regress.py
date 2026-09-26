#!/usr/bin/env python3
"""場の非退行判定 (plan config-key-pruning §6.2' の実体化)。

**同一バイナリ・同一設定の反復 run からノイズ床を作り、候補 run の差がその何倍かで判定する**。CFD の残差は
float の `atomicAdd` で集積するので、無変更でもビット一致しない。したがって「絶対差がいくつ」ではなく
「同じものを 2 回回した差に対して何倍か」で見る。基準は run を回す前に固定すること。

  - 比較量は **全量を判定**する (1 量でも外れたら不合格)。既定は保存量 + 原始量、あれば乱流量も。
  - `--boundary` で**境界出力ファイル** (`res_<名前>_<physID>_<step>.h5`, `outputHDFflg: 1` の bcond) も比較する。
    壁せん断応力 `twall_*`・壁熱流束 `qwall`・`utau`・`ypls` は保存量に出ない**出力専用量**なので、壁経路の
    回帰ではこちらが本体。量の名前は `<接頭辞>/<量>` で表示する。
  - ノルムは相対 L2 と相対 L∞ の**両方**。正規化は基準 run の L2 ノルム / 最大絶対値。
  - ノイズ床は反復 run の**全ペアの最大**。候補は**全反復に対する最大差**で測る。
  - `--candidate` は**複数指定できる**。指定するとノイズ床は「基準側の全ペア」と「候補側の全ペア」の**大きい方**、
    比較は「基準×候補の全ペアの最大」になる。**カオス的な run (LES/DES) では片側 3 本の床が桁で足りない**ことが
    あるので、両側 3 本以上にすること (実例: 周期丘 DDES 400 step の壁せん断 `twall_x` は、新側 3 本の床
    1.81e-3 に対し 5 本にすると 2.70e-2 と 15 倍になり、3 本での判定は偽の不合格を出した)。
  - ノイズ床が厳密に 0 (ビット一致) の量は、候補にもビット一致を要求する。
  - **数値判定の前にデータ不備を検査する**: 欠落・形状不一致・非有限値 (NaN/Inf) はどの run のどの量でも
    その場で終了コード 2 とし、数値判定に進まない (`perf_regress.py` と同じ作法; `max(0.0, NaN)` が `0.0` になって
    異常が消えるため、比較の中で検出することはできない)。**検査は比較する全量に掛ける** — 任意量 (`roK` など) も、
    **どれか 1 本の run にあれば全 run に要る**ものとして扱う (片側だけ欠けるのは異常であって「比較しない」理由にならない)。
  - **数値的にゼロの量** (最大絶対値 < `--zero-scale`, 既定 1e-20) は相対ノルムが意味を持たないので `zero` と表示して
    判定から外す (例: 平面 2D の `roUz` は 1e-36 の非正規化数で、相対差は幾らでも大きくなる)。判定から外したことは
    表に残す。候補側も同じ閾値を下回ることを確認する。

使い方:
  check_field_regress.py --repeat RUN1 RUN2 [RUN3 ...] --candidate RUN [--step N] [--factor 2.0]
                         [--quantities ro,roUx,...]
VERDICT PASS / FAIL, exit 0/1。
"""
import argparse
import glob
import os
import sys

DEFAULT_Q = ['ro', 'roUx', 'roUy', 'roUz', 'roe', 'P', 'T']
OPTIONAL_Q = ['roK', 'roOmega', 'h0']
# 境界出力 (res_<名前>_<physID>_<step>.h5) で壁経路の回帰に要る量。存在するものだけ比較する。
BOUNDARY_Q = ['twall_x', 'twall_y', 'twall_z', 'qwall', 'utau', 'ypls', 'Ps', 'Ts', 'ro', 'roUx', 'roUy', 'roUz', 'roe']


def last_step(run):
    steps = [int(os.path.basename(f)[4:-3]) for f in glob.glob(os.path.join(run, 'res_*.h5'))
             if os.path.basename(f)[4:-3].isdigit()]
    return max(steps) if steps else None


def boundary_files(run, step):
    """{接頭辞: パス} — 境界出力 res_<名前>_<physID>_<step>.h5 を拾う (res_<step>.h5 は除く)。"""
    out = {}
    for f in sorted(glob.glob(os.path.join(run, f'res_*_{step}.h5'))):
        stem = os.path.basename(f)[4:-(len(str(step)) + 4)].rstrip('_')
        if stem:
            out[stem] = f
    return out


def load_boundary(run, step, exclude=frozenset()):
    """境界出力を読む。`exclude` は 'ypls' (全境界) か 'wall_4/ypls' (その境界だけ) の集合。

    境界量は BOUNDARY_Q 固定で `--quantities` からは選べない (2026-09-26, codex plan M4)。
    **意図的に定義を変えた量を非退行判定から外す**ための口。外した量は別ゲートで検査すること。
    """
    import h5py
    out = {}
    for stem, path in boundary_files(run, step).items():
        with h5py.File(path, 'r') as f:
            if 'VALUE' not in f:
                continue
            for q in BOUNDARY_Q:
                if q in exclude or f'{stem}/{q}' in exclude:
                    continue
                if q in f['VALUE']:
                    out[f'{stem}/{q}'] = f['VALUE'][q][()].astype('float64').ravel()
    return out


def load(run, step, quantities):
    import h5py
    p = os.path.join(run, f'res_{step}.h5')
    out = {}
    with h5py.File(p, 'r') as f:
        g = f['VALUE']
        for q in quantities:
            if q in g:
                out[q] = g[q][()].astype('float64').ravel()
    return out


def data_check(data, runs, required):
    """数値判定の**前**の不備検査。戻り: 問題の一覧 (空なら健全)。

    NaN は比較の中では検出できない (`max(0.0, NaN)` が `0.0` になり異常が消える) ので、ここで落とす。
    検査対象は `required` に加えて、**どれか 1 本の run に出ている量すべて**。片側だけ欠けるのは異常であって
    「比較から外す」理由にならない (codex result M1)。"""
    import numpy as np
    problems = []
    ref = runs[0]
    shapes = {}
    seen = set(required)
    for r in runs:
        seen |= set(data[r])
    for r in runs:
        for q in sorted(seen):
            if q not in data[r]:
                problems.append(f'{r}: 比較量 {q} が出力に無い'
                                + ('' if q in required else ' (他の run には出ている)'))
                continue
            v = data[r][q]
            shapes.setdefault(q, (ref, v.shape))
            if v.shape != shapes[q][1]:
                problems.append(f'{r}: {q} の形状 {v.shape} が {shapes[q][0]} の {shapes[q][1]} と違う')
            nn = int((~np.isfinite(v)).sum())
            if nn:
                problems.append(f'{r}: {q} に非有限値 (NaN/Inf) が {nn} 個')
    return problems


def norms(a, b):
    """(相対 L2, 相対 L∞)。正規化は a 側 (基準) の大きさで行う。"""
    import numpy as np
    d = a - b
    l2 = float(np.sqrt((d * d).sum()))
    li = float(np.abs(d).max()) if d.size else 0.0
    n2 = float(np.sqrt((a * a).sum()))
    ni = float(np.abs(a).max()) if a.size else 0.0
    return (l2 / n2 if n2 > 0 else l2), (li / ni if ni > 0 else li)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repeat', nargs='+', required=True, help='同一バイナリ・同一設定の反復 run (2 本以上)')
    ap.add_argument('--candidate', nargs='+', required=True,
                    help='比較する run (旧バイナリ / 変更後 など)。複数指定すると候補側の反復ペアもノイズ床に入る')
    ap.add_argument('--step', type=int, default=None, help='比較する step (既定: 反復 run 共通の最終 step)')
    ap.add_argument('--factor', type=float, default=2.0, help='許容 = ノイズ床 × factor (既定 2)')
    ap.add_argument('--boundary', action='store_true',
                    help='境界出力ファイル (壁せん断応力・壁熱流束・utau・y+) も比較する。壁経路の回帰では必須')
    ap.add_argument('--zero-scale', type=float, default=1e-20,
                    help='最大絶対値がこれ未満の量は「数値的にゼロ」として判定から外す (既定 1e-20)')
    ap.add_argument('--quantities', default=None, help='カンマ区切りで上書き (**体積量のみ**)')
    ap.add_argument('--exclude-boundary', default='',
                    help="境界量を判定から外す。カンマ区切りで 'ypls' か 'wall_4/ypls' の形。"
                         '境界量は BOUNDARY_Q 固定で --quantities では選べないため。'
                         '**外した量は別ゲートで検査すること** (外した旨は出力に出る)')
    a = ap.parse_args()

    cands = list(a.candidate)
    runs = list(a.repeat) + cands
    step = a.step
    if step is None:
        steps = [last_step(r) for r in runs]
        if None in steps or len(set(steps)) != 1:
            print(f'FAIL: run ごとに最終 step が違う {dict(zip(runs, steps))}。--step で指定すること')
            sys.exit(1)
        step = steps[0]
    explicit = a.quantities.split(',') if a.quantities else None
    quantities = explicit if explicit else DEFAULT_Q + OPTIONAL_Q
    data = {r: load(r, step, quantities) for r in runs}
    required = list(explicit or DEFAULT_Q)
    if a.boundary:
        excl_b = {t.strip() for t in a.exclude_boundary.split(',') if t.strip()}
        if excl_b:
            print('境界量のうち判定から外したもの: %s  (**別ゲートで検査すること**)'
                  % ', '.join(sorted(excl_b)))
        bnd = {r: load_boundary(r, step, exclude=excl_b) for r in runs}
        names = sorted(set().union(*[set(v) for v in bnd.values()])) if bnd else []
        if not names:
            print(f'FAIL: --boundary を指定したが、step {step} の境界出力 (res_<名前>_<physID>_{step}.h5) が無い。'
                  ' bcondConfig の outputHDFflg を 1 にすること')
            sys.exit(2)
        for r in runs:
            data[r].update(bnd[r])
        quantities = list(quantities) + names
        required += names
    problems = data_check(data, runs, required)
    if problems:
        print(f'--- データ不備 {len(problems)} 件; 数値判定は行わない')
        for q in problems:
            print('  ' + q)
        sys.exit(2)
    common = sorted(set().union(*[set(data[r]) for r in runs]),
                    key=lambda q: (quantities.index(q) if q in quantities else len(quantities), q))

    print(f'step {step} / 基準側 {len(a.repeat)} 本・候補側 {len(cands)} 本 / 許容 = ノイズ床 × {a.factor:g}')
    if len(cands) == 1:
        print('  注意: 候補側が 1 本なのでノイズ床は基準側だけから測っている。カオス的な run では過小評価になる')
    print(f"{'量':<10} {'ノイズ床 L2':>12} {'候補 L2':>12} {'比':>7} "
          f"{'ノイズ床 L∞':>12} {'候補 L∞':>12} {'比':>7}  判定")
    bad = 0
    for q in common:
        import numpy as np
        scale = max(float(np.abs(data[r][q]).max()) for r in runs)
        if scale < a.zero_scale:
            print(f'{q:<10} {"":>12} {"":>12} {"":>7} {"":>12} {"":>12} {"":>7}  zero (max|.|={scale:.2e})')
            continue
        floor2 = floori = 0.0
        for group in (a.repeat, cands):
            for i in range(len(group)):
                for j in range(i + 1, len(group)):
                    n2, ni = norms(data[group[i]][q], data[group[j]][q])
                    floor2, floori = max(floor2, n2), max(floori, ni)
        cand2 = candi = 0.0
        for r in a.repeat:
            for c in cands:
                n2, ni = norms(data[r][q], data[c][q])
                cand2, candi = max(cand2, n2), max(candi, ni)
        ok2 = (cand2 == 0.0) if floor2 == 0.0 else (cand2 <= a.factor * floor2)
        oki = (candi == 0.0) if floori == 0.0 else (candi <= a.factor * floori)
        r2 = (cand2 / floor2) if floor2 > 0 else float('inf') if cand2 > 0 else 0.0
        ri = (candi / floori) if floori > 0 else float('inf') if candi > 0 else 0.0
        if not (ok2 and oki):
            bad += 1
        print(f'{q:<10} {floor2:12.4e} {cand2:12.4e} {r2:7.2f} {floori:12.4e} {candi:12.4e} {ri:7.2f}  '
              f'{"ok" if (ok2 and oki) else "**NG**"}')
    print(f'\nVERDICT: {"PASS" if bad == 0 else f"FAIL ({bad} 量が許容外)"}')
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
