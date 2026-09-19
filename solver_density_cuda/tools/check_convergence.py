#!/usr/bin/env python3
"""
収束判定ツール (AGENTS.md「収束確認 (必須)」の実体化)。

forge の run ディレクトリの residual_history.csv を読み、**全保存量の残差列**
(rms_ro, rms_roUx, rms_roUy, rms_roUz, rms_roe, RANS時 rms_roK/rms_roOmega, 化学種 rms_roY*,
凝縮 rms_rog_*/rms_roQ*_*) について
初期値・最終値・低下桁数・末尾トレンド (falling/flat/rising) を出し、明確な VERDICT を返す。

目的: 「rms_ro と NaN だけ見て収束と判断する」ことを防ぐ (AGENTS.md 違反の常習を防止)。
結果を「収束した」「一致した」と報告する前に必ず本ツールを通すこと。

使い方:
  python3 tools/check_convergence.py <run_dir> [<run_dir2> ...]
  python3 tools/check_convergence.py --drop 4 --tail 0.2 run_a run_b   # 判定閾値を調整
終了コード: 全 run が PASS なら 0、1つでも未収束/NaN があれば 1 (CI/スクリプトで使える)。

判定基準 (既定):
  - NaN/Inf が無い。
  - 各非ゼロ保存量残差が**系列のピークから** >= --drop 桁 (既定 3) 低下している
    (step 0 基準ではない: IC の作り方で step 0 が過渡ピークより桁違いに小さい
     成分があり [node の準1D IC は Uy≈0]、初期比だと誤って停滞判定になる)。
  - 末尾 (--tail, 既定 20%) が rising でない (flat/falling)。
"all-zero" 列 (例: 2D の rms_roUz) は判定から除外する。

収束場からの restart (交差 restart など) 用: --from-floor REF_RUN (別名 --reference-floor)
  ピークからの低下桁数は収束場から再開した run では原理的に 0 なので、代わりに **参照 run の末尾床**
  (REF_RUN の residual_history.csv の末尾 --tail 平均) を基準に、各残差列が run の全期間で床の
  --floor-factor 倍 (既定 1.5) 以内に留まり (ピーク ≤ factor×床)、末尾平均も factor×床 以内なら PASS。
  列ごとの床比 (末尾平均/床, ピーク/床) を表示する。NaN/Inf は DIVERGED、超過は NOT CONVERGED。
  例: python3 tools/check_convergence.py --from-floor run_0476 run_0478 run_0479
"""
import csv, math, sys, argparse, os

CONSERVED = ['rms_ro', 'rms_roUx', 'rms_roUy', 'rms_roUz', 'rms_roe',
             'rms_roK', 'rms_roOmega']


def load_series(path):
    rows = [r for r in csv.DictReader(open(path)) if r.get('phase') == 'outer_end']
    if not rows:
        # outer_end が無い構成 (純 explicit) では最終 inner / 全行を使う
        rows = list(csv.DictReader(open(path)))
    cols = {}
    for c in CONSERVED:
        if rows and c in rows[0]:
            cols[c] = [float(r[c]) for r in rows]
    # 化学種 (rms_roY*)・受動トレーサ (rms_roXi) と凝縮 (rms_rog_*, rms_roQ{0,1,2}_*) の保存量残差も検査する (存在時)。
    # 凝縮 run で NS 5 本 + SST 2 本だけ見て「収束」と判定していた穴 (codex 指摘 2026-09-10) を塞ぐ。
    if rows:
        for c in rows[0].keys():
            if c in cols or not c.startswith('rms_') or c.startswith('rms_dq_'):
                continue
            if c.startswith('rms_roY') or c.startswith('rms_rog_') or c.startswith('rms_roQ') or c == 'rms_roXi':
                cols[c] = [float(r[c]) for r in rows]
    return rows, cols


def _tail_rise(ser, tail_frac):
    """末尾 2*tail_frac 窓を 2 通りに評価し `(spike, slow)` を返す。

    - `spike`: **末端の急増**。窓末尾の数点の最大が窓中央値の 10 倍を超えるか。
      緩やかなトレンドの検定では捕まらない「最後の数 step での爆発」を独立に見る
      (実例: `case/37.pintle_nozzle/run_0009` の `rms_roUx` は末尾 5 点で 2.77e-4 → 7.58e22 と
      26 桁跳ねるのに、窓 444 点の回帰では上昇 0.13 桁・散らばり 1.58 桁に埋もれる)。
    - `slow`: **持続的な緩い上昇**。log10|値| の線形回帰の傾きが有意 (t > 3) で、窓全体の
      上昇が 10 % を超えるか。散らばりが大きいことは上昇が無い証明にならないので、
      「上昇 > 散らばり」ではなく回帰の有意性で見る (codex plan レビュー Major 1)。

    点数が足りない・分散が無い場合は **None** を返す (判定は既存の 2 窓平均比だけに委ねる。
    数 step で落ちた run の系列がこれに当たる)。"""
    w = [abs(x) for x in ser[int(len(ser) * (1 - 2 * tail_frac)):] if x > 0.0]
    m = len(w)
    if m < 8:
        return None
    y = [math.log10(v) for v in w]
    xb = (m - 1) / 2.0
    yb = sum(y) / m
    sxx = sum((i - xb) ** 2 for i in range(m))
    if sxx <= 0.0:
        return None
    slope = sum((i - xb) * (yi - yb) for i, yi in enumerate(y)) / sxx
    resid = [yi - (yb + slope * (i - xb)) for i, yi in enumerate(y)]
    var = sum(r * r for r in resid) / (m - 2) if m > 2 else 0.0
    se = math.sqrt(var / sxx) if var > 0.0 and sxx > 0.0 else 0.0
    rise = slope * (m - 1)
    tval = (rise / (m - 1)) / se if se > 0.0 else (math.inf if slope > 0.0 else 0.0)
    slow = tval > 3.0 and rise > math.log10(1.10)

    ntip = max(3, int(math.ceil(0.02 * m)))
    srt = sorted(w)
    med = srt[m // 2] if m % 2 else 0.5 * (srt[m // 2 - 1] + srt[m // 2])
    tip = max(w[-ntip:])
    spike = med > 0.0 and tip > 10.0 * med
    return spike, slow


# **必須の保存量残差列**。これが 1 つも無い / 欠けている CSV を合格にしてはいけない
# (2026-09-19 codex: `step,phase` だけの CSV や `rms_ro` だけの CSV が ok=True になっていた)。
REQUIRED_COLS = ('rms_ro', 'rms_roUx', 'rms_roUy', 'rms_roUz', 'rms_roe')


def analyze(path, min_drop, tail_frac):
    rows, cols = load_series(path)
    if not rows:
        return None
    laststep = rows[-1].get('step', '?')
    report = {}
    ok = True
    # --- 入力の完全性検査 (合格の前提。欠けていたら「判定不能」= ok=False) ---
    if not cols:
        report['(入力)'] = ('残差列が 1 つも無い  <-- 判定不能', False)
        return laststep, report, False, False, False, False
    missing = [c for c in REQUIRED_COLS if c not in cols]
    if missing:
        report['(入力)'] = ('必須の保存量残差列が無い: %s  <-- 判定不能' % ', '.join(missing), False)
        ok = False
    any_nan = False
    any_stalled = False
    any_converging = False
    for c, ser in cols.items():
        nz = [v for v in ser if v != 0.0]
        if not nz:
            report[c] = ('all-zero (inactive, skip)', True)
            continue
        if any(math.isnan(v) or math.isinf(v) for v in ser):
            report[c] = ('NaN/Inf present  <-- DIVERGED', False)
            ok = False; any_nan = True
            continue
        init, fin = ser[0], ser[-1]
        n = len(ser)
        a = ser[int(n * (1 - tail_frac)):]
        b = ser[int(n * (1 - 2 * tail_frac)):int(n * (1 - tail_frac))] or a
        ma = sum(abs(x) for x in a) / len(a)
        mb = sum(abs(x) for x in b) / len(b)
        # rising = 窓平均が 5% 増 **かつ** プラトーを実際に離脱している (系列最小の
        # 2 倍超)。後者のガードが無いと、深く収束したプラトー (例: 3.4 桁低下後の
        # 1e-7 台) のリミットサイクル呼吸 ±5% を「発散傾向」と誤判定する
        # (run 毎の非決定性で同一設定が PASS/NOT CONVERGED に割れる実害が出た
        # 2026-08-15)。本物のリバウンド発散は最小値の 2 倍を速やかに超えるので
        # 検出力は保たれる。
        smin = min(abs(x) for x in ser if x != 0.0) if any(x != 0.0 for x in ser) else 0.0
        # rising の 3 条件目: **末端の急増か、有意な緩い上昇があること**。これが無いと、
        # プラトー自体が数倍の幅で揺れている列では 2 窓平均の大小がジッタの位相だけで決まり、
        # 判定が再現しない。実証 (2026-09-19): case/46 の run_0193 と run_0195 は設定重複で
        # 同一形状・同一レシピになっており rms_roY1 の分布も同一 (後半中央値 1.37e-6 / 1.33e-6、
        # p5-p95 一致、帯 8.5 倍) だったのに、ma/mb が 0.961 と 1.053 に割れて flat / rising に
        # 分かれ、ゲートが PASS / FAIL に反転した。判定は `_tail_rise` に分離してある。
        rise = _tail_rise(ser, tail_frac)
        over_jitter = True if rise is None else (rise[0] or rise[1])
        trend = ('rising' if (ma > mb * 1.05 and ma > 2.0 * smin and over_jitter)
                 else ('flat' if ma > mb * 0.9 else 'falling'))

        # init==0 (例: アライン格子で Uy が初期厳密 0) も下のピーク基準で判定する。旧特例 (rising でなければ
        # 合格) は [0,1,1,...] のように立ち上がって落ちない列を合格にしてしまった (codex 指摘 2026-09-10)。
        # 全期間ゼロの列だけを上の all-zero で除外する。

        # 低下桁数は **step 0 ではなく系列のピーク**から測る。IC の作り方によっては
        # ある成分の step 0 残差が過渡ピークより桁違いに小さいことがあり (node の
        # 準 1D IC は半径方向速度がほぼ厳密 0 なので rms_roUy の init が 4.5e-6、
        # 一方 step 2 のピークは 4.8e-2 = 4 桁上)、初期比だと実際に 5.0 桁落ちて
        # いる収束列を「0.9 桁で停滞」と誤判定する (2026-08-17 実測)。ピークは
        # 「その成分が実際にどこから落ちたか」なので物理的にも正しい尺度。
        # step 0 がピークの通常ケースでは値は変わらない (後方互換)。
        peak = max(abs(x) for x in ser)
        # 低下桁数は **末尾窓の代表値**で測る。最終 1 点だけを使うと、[1]*19+[0] のように
        # 末尾が 1 点だけ 0 の系列で drop=inf になり合格してしまう (2026-09-19 codex)。
        # 代表値は末尾窓の |値| の中央値。窓に 0 が混じっても中央値は 0 になりにくく、
        # 本当に全体が 0 まで落ちた列は下の all-zero か「窓中央値 0」で判定不能にする。
        tail_abs = sorted(abs(x) for x in a)
        fin_rep = tail_abs[len(tail_abs) // 2]
        if fin_rep == 0.0:
            report[c] = ('末尾窓の代表値が 0 (残差が数値的にゼロ)  <-- 判定不能', False)
            ok = False
            continue
        drop = math.log10(peak / fin_rep)
        col_ok = drop >= min_drop and trend != 'rising'
        ok = ok and col_ok
        # 未達の理由を区別: falling=収束途中(あと steps)、flat=停滞、rising=発散傾向
        status = ''
        if not col_ok:
            if trend == 'rising':   status = '  <-- RISING (divergent)'
            elif trend == 'flat':   status = '  <-- STALLED (plateau)'; any_stalled = True
            else:                   status = '  <-- still converging'; any_converging = True
        pk = "" if abs(peak - init) <= 1e-30 * max(abs(init), 1.0) else f" peak={peak:.2e}"
        report[c] = (f"init={init:.2e}{pk} fin={fin:.2e} drop={drop:4.1f}dec {trend:7s}{status}", col_ok)
    return laststep, report, ok, any_nan, any_stalled, any_converging


def build_segment_csv(run_dir):
    """`stage_manifest.json` の**最後の区間**を step オフセット付きで連結した CSV を書いて返す。

    段名でなく**実効設定 (hard キー)** が同じ連続区間だけを繋ぐので、別の方程式・BC の過渡を
    本段の低下桁数の基準にしてしまう事故が起きない (2026-09-19 codex Major 1)。
    列が段で違う場合は**共通列に落とさず**、区間内で列集合が一致することを要求する
    (共通列への縮退は検出したい誤合格を再導入する: 同 Major 2)。
    """
    import csv as _csv
    import json as _json
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from stage_manifest import segments
    p = os.path.join(run_dir, 'stage_manifest.json')
    if not os.path.exists(p):
        return None
    segs = segments(_json.load(open(p)))
    if not segs:
        return None
    seg = segs[-1]
    rows, hdr, off = [], None, 0
    used = []
    for st in seg:
        f = os.path.join(run_dir, st['history'])
        if not os.path.exists(f):
            continue
        r = list(_csv.reader(open(f)))
        if len(r) < 2:
            continue
        h = [c.strip() for c in r[0]]
        if hdr is None:
            hdr = h
        elif h != hdr:
            print(f"  [segment] 段 {st['tag']} の列が区間内で一致しない -> 連結しない "
                  f"(共通列に落とすと誤合格を再導入する)")
            return None
        body = [x for x in r[1:] if x and x[0].strip().lstrip('-').isdigit()]
        if not body:
            continue
        n1 = int(body[-1][0])
        for x in body:
            x = list(x); x[0] = str(int(x[0]) + off); rows.append(x)
        used.append(st['tag'])
        off += n1 + 1
    if not rows:
        return None
    out = os.path.join(run_dir, 'residual_history_segment.csv')
    with open(out, 'w', newline='') as f:
        w = _csv.writer(f); w.writerow(hdr); w.writerows(rows)
    print(f"  [segment] 判定区間 = {' -> '.join(used)}  ({len(rows)} 行) -> {out}")
    return out


def reference_floor(path, tail_frac):
    """参照 run の各残差列の末尾床 (末尾 tail_frac の |値| 平均)。
    戻り値 (floor, zero_cols): floor は非ゼロ列だけ、zero_cols は all-zero (非活性) 列の集合。"""
    rows, cols = load_series(path)
    floor = {}
    zero_cols = set()
    for c, ser in cols.items():
        if not any(v != 0.0 for v in ser):
            zero_cols.add(c)
            continue
        a = ser[int(len(ser) * (1 - tail_frac)):]
        floor[c] = sum(abs(x) for x in a) / len(a)
    return floor, zero_cols


def analyze_from_floor(path, floor, zero_cols, factor, tail_frac):
    """収束場からの restart 判定: 全期間ピークと末尾平均が参照床の factor 倍以内なら列 PASS。
    列対応 (codex result-2 m1): 参照に無い列・参照で all-zero だった列が対象で非ゼロなら FAIL (前提不成立を PASS にしない)。
    参照の非ゼロ列が対象に無い場合も FAIL。"""
    rows, cols = load_series(path)
    if not rows:
        return None
    laststep = rows[-1].get('step', '?')
    report = {}
    ok = True
    any_nan = False
    for c in floor:
        if c not in cols:
            report[c] = ('missing in target (reference has a nonzero floor)  <-- COLUMN MISMATCH', False)
            ok = False
    for c, ser in cols.items():
        if not any(v != 0.0 for v in ser):
            report[c] = ('all-zero (inactive, skip)', True)
            continue
        if any(math.isnan(v) or math.isinf(v) for v in ser):
            report[c] = ('NaN/Inf present  <-- DIVERGED', False)
            ok = False; any_nan = True
            continue
        if c not in floor:
            if c in zero_cols:
                report[c] = ('reference column was all-zero but target is nonzero  <-- COLUMN MISMATCH', False)
            else:
                report[c] = ('no such column in reference  <-- COLUMN MISMATCH', False)
            ok = False
            continue
        peak = max(abs(x) for x in ser)
        a = ser[int(len(ser) * (1 - tail_frac)):]
        tail_mean = sum(abs(x) for x in a) / len(a)
        r_tail = tail_mean / floor[c] if floor[c] > 0 else float('inf')
        r_peak = peak / floor[c] if floor[c] > 0 else float('inf')
        col_ok = (r_tail <= factor) and (r_peak <= factor)
        ok = ok and col_ok
        status = '' if col_ok else ('  <-- ABOVE FLOOR (peak)' if r_peak > factor else '  <-- ABOVE FLOOR (tail)')
        report[c] = (f"ref_floor={floor[c]:.2e} tail={tail_mean:.2e} (x{r_tail:.2f}) peak={peak:.2e} (x{r_peak:.2f}) "
                     f"fin={ser[-1]:.2e}{status}", col_ok)
    return laststep, report, ok, any_nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_dirs', nargs='+')
    ap.add_argument('--drop', type=float, default=3.0, help='required orders-of-magnitude drop')
    ap.add_argument('--tail', type=float, default=0.2, help='tail fraction for trend check')
    ap.add_argument('--from-floor', '--reference-floor', dest='from_floor', default=None, metavar='REF_RUN',
                    help='収束場からの restart 判定: REF_RUN の末尾床 (tail 平均) の --floor-factor 倍以内に全期間留まれば PASS')
    ap.add_argument('--floor-factor', type=float, default=1.5, help='--from-floor の許容倍率 (既定 1.5)')
    ap.add_argument('--segment', action='store_true',
                    help='段階起動の run で **stage_manifest.json の最後の区間** '
                         '(方程式・BC・空間離散化が同一の連続区間) を連結して判定する。'
                         '段名でなく実効設定で区間を決めるので、別 BC の過渡を基準にしない')
    args = ap.parse_args()

    floor = None
    if args.from_floor:
        ref_path = args.from_floor if args.from_floor.endswith('.csv') else os.path.join(args.from_floor, 'residual_history.csv')
        if not os.path.exists(ref_path):
            print(f"[{args.from_floor}] NO residual_history.csv (reference)"); sys.exit(1)
        # 参照 run 自身が通常判定 (--drop/--tail) を PASS していることを要求する (codex result-2 m1: 未収束 run を参照にすると
        # 「床に留まった」だけで PASS になる)。不成立なら判定を拒否 (exit 2)。
        ref_res = analyze(ref_path, args.drop, args.tail)
        if ref_res is None:
            print(f"[{args.from_floor}] empty reference residual file"); sys.exit(2)
        ref_laststep, ref_report, ref_ok, ref_nan, _, _ = ref_res
        if not ref_ok:
            print(f"=== reference {args.from_floor}  [last step {ref_laststep}]  -> "
                  f"{'DIVERGED (NaN/Inf)' if ref_nan else 'NOT CONVERGED'} ===")
            for c, (msg, _) in ref_report.items():
                print(f"  {c:12s}: {msg}")
            print("\nREFUSED: --from-floor requires a reference run that PASSES the normal criterion "
                  f"(--drop {args.drop}, --tail {args.tail}); its residual floor is not a converged floor.")
            sys.exit(2)
        floor, zero_cols = reference_floor(ref_path, args.tail)
        print(f"reference floor from {args.from_floor} (normal criterion PASS at step {ref_laststep}; tail {args.tail:.0%} mean), "
              f"factor {args.floor_factor}; columns {sorted(floor)}; inactive {sorted(zero_cols)}")

    all_pass = True
    for rd in args.run_dirs:
        if args.segment and not rd.endswith('.csv'):
            path = build_segment_csv(rd)
            if path is None:
                print(f"[{rd}] stage_manifest.json が無い -> --segment は使えない "
                      f"(判定区間を人が明示すること)")
                worst = max(worst, 2)
                continue
        else:
            path = rd if rd.endswith('.csv') else os.path.join(rd, 'residual_history.csv')
        if not os.path.exists(path):
            print(f"[{rd}] NO residual_history.csv"); all_pass = False; continue
        if floor is not None:
            res = analyze_from_floor(path, floor, zero_cols, args.floor_factor, args.tail)
            if res is None:
                print(f"[{rd}] empty residual file"); all_pass = False; continue
            laststep, report, ok, any_nan = res
            verdict = ('PASS (within %.1fx of reference floor %s)' % (args.floor_factor, args.from_floor) if ok
                       else 'DIVERGED (NaN/Inf)' if any_nan
                       else 'NOT CONVERGED (residual left the reference floor, or column mismatch)')
            print(f"\n=== {rd}  [last step {laststep}]  -> {verdict} ===")
            for c, (msg, _) in report.items():
                print(f"  {c:12s}: {msg}")
            all_pass = all_pass and ok
            continue
        res = analyze(path, args.drop, args.tail)
        if res is None:
            print(f"[{rd}] empty residual file"); all_pass = False; continue
        laststep, report, ok, any_nan, any_stalled, any_converging = res
        if ok:
            verdict = 'PASS (converged)'
        elif any_nan:
            verdict = 'DIVERGED (NaN/Inf)'
        elif any_stalled:
            verdict = 'NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)'
        elif any_converging:
            verdict = 'NOT CONVERGED (still converging — run more steps)'
        else:
            verdict = 'NOT CONVERGED'
        print(f"\n=== {rd}  [last step {laststep}]  -> {verdict} ===")
        for c, (msg, _) in report.items():
            print(f"  {c:12s}: {msg}")
        all_pass = all_pass and ok
    print(f"\nOVERALL: {'ALL PASS' if all_pass else 'CHECK FAILURES ABOVE'}")
    sys.exit(0 if all_pass else 1)


if __name__ == '__main__':
    main()
