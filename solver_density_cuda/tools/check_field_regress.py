#!/usr/bin/env python3
"""場の非退行判定 (plan config-key-pruning §6.2' の実体化)。

**同一バイナリ・同一設定の反復 run からノイズ床を作り、候補 run の差がその何倍かで判定する**。CFD の残差は
float の `atomicAdd` で集積するので、無変更でもビット一致しない。したがって「絶対差がいくつ」ではなく
「同じものを 2 回回した差に対して何倍か」で見る。基準は run を回す前に固定すること。

  - 比較量は **全量を判定**する (1 量でも外れたら不合格)。既定は保存量 + 原始量、あれば乱流量と壁量も。
  - ノルムは相対 L2 と相対 L∞ の**両方**。正規化は基準 run の L2 ノルム / 最大絶対値。
  - ノイズ床は反復 run の**全ペアの最大**。候補は**全反復に対する最大差**で測る。
  - ノイズ床が厳密に 0 (ビット一致) の量は、候補にもビット一致を要求する。
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
OPTIONAL_Q = ['roK', 'roOmega', 'Tau_Wall', 'Qw_Wall', 'h0']


def last_step(run):
    steps = [int(os.path.basename(f)[4:-3]) for f in glob.glob(os.path.join(run, 'res_*.h5'))
             if os.path.basename(f)[4:-3].isdigit()]
    return max(steps) if steps else None


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
    ap.add_argument('--candidate', required=True, help='比較する run (旧バイナリ / 変更後 など)')
    ap.add_argument('--step', type=int, default=None, help='比較する step (既定: 反復 run 共通の最終 step)')
    ap.add_argument('--factor', type=float, default=2.0, help='許容 = ノイズ床 × factor (既定 2)')
    ap.add_argument('--zero-scale', type=float, default=1e-20,
                    help='最大絶対値がこれ未満の量は「数値的にゼロ」として判定から外す (既定 1e-20)')
    ap.add_argument('--quantities', default=None, help='カンマ区切りで上書き')
    a = ap.parse_args()

    runs = list(a.repeat) + [a.candidate]
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
    common = [q for q in quantities if all(q in data[r] for r in runs)]
    missing = [q for q in (explicit or DEFAULT_Q) if q not in common]
    if missing:
        print(f'FAIL: 必須の比較量が出力に無い: {missing} (output.level を上げること)')
        sys.exit(1)

    print(f'step {step} / 反復 {len(a.repeat)} 本 / 許容 = ノイズ床 × {a.factor:g}')
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
        for i in range(len(a.repeat)):
            for j in range(i + 1, len(a.repeat)):
                n2, ni = norms(data[a.repeat[i]][q], data[a.repeat[j]][q])
                floor2, floori = max(floor2, n2), max(floori, ni)
        cand2 = candi = 0.0
        for r in a.repeat:
            n2, ni = norms(data[r][q], data[a.candidate][q])
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
