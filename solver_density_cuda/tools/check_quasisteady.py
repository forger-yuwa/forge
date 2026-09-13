#!/usr/bin/env python3
"""
準定常 (quasi-steady) 判定ツール (AGENTS.md「準定常確認 (必須)」の実体化)。

`check_convergence.py` が **残差** の収束を見るのに対し、本ツールは **報告する派生量そのもの**
(衝撃位置 shock / 上下非対称 asym / 最大マッハ machmax / 最大静圧 pmax 等) を、run ディレクトリ内の
全 `res_*.h5` スナップショット時系列で評価し、その量が **頭打ち (定常化) したか** を判定する。

目的: 「残差プラトー = 報告していい」と誤認して、**過渡ピーク / ドリフト中の量を定常値として報告する**ことを
防ぐ (本ツールが無かったため、過渡 0.25 の非対称を定常偏りと誤報告した事例があった)。衝撃位置・非対称・
CL/CD・massflux・推力・peak μt 等を「○○だ」と報告する応答は、必ず本ツールの VERDICT を引用すること。

使い方:
  python3 tools/check_quasisteady.py <run_dir> [--quantity shock,asym] [--mesh mesh.h5]
  python3 tools/check_quasisteady.py <run_dir> --tail 0.4 --drift 0.05 --osc 0.10
  python3 tools/check_quasisteady.py --series-csv <run_dir>/force_history.csv --series-cols C_T_with_shear,C_L,C_M
    (CSV 系列モード: `step` 列を持つ CSV の指定列を **同じ classify** で判定する。設計チェーンの力係数履歴
     [`forge_design.metrics.sern_forces` が書く `force_history.csv`] など、res_*.h5 から直接抽出できない
     派生量を正式ツールの VERDICT で報告するための入口。非有限値を含む列は NONFINITE (最重症) にする)
終了コード: 全量が STEADY なら 0、1つでも DRIFTING/UNSETTLED があれば 1 (CI/スクリプトで使える)。
OSCILLATING (リミットサイクル) は 0 扱いだが「平均±振幅」で報告すること (瞬時値で報告しない)。

判定 (各量の末尾 tail について):
  DRIFTING            : 末尾が単調トレンドで |傾き×tail幅|/平均 > --drift (まだ動いている。run を伸ばす)
  OSCILLATING         : 末尾の振れ (max-min)/平均 > --osc だがトレンド小 (リミットサイクル → 平均±振幅で報告)
  TRANSIENT-UNSETTLED : スナップショットが少なすぎる / 全系列の極値が末尾にある (過渡が減衰しきっていない)
  STEADY              : 末尾が許容内で平坦
"""
import sys, os, glob, argparse, math
import numpy as np
import h5py

GAMMA = 1.4


def find_mesh(run_dir, explicit):
    if explicit:
        return explicit if os.path.isabs(explicit) else os.path.join(run_dir, explicit)
    # 入力メッシュ = res_ でない .h5 で /CELLS/centCoords を持つもの
    for f in sorted(glob.glob(os.path.join(run_dir, '*.h5'))):
        if os.path.basename(f).startswith('res_'):
            continue
        try:
            with h5py.File(f, 'r') as h:
                if '/CELLS/centCoords' in h:
                    return f
        except Exception:
            pass
    return None


def res_files(run_dir):
    # 主スナップショット res_<step>.h5 のみ対象 (res_wall_*/res_outlet_*/res_nan_* 等の境界・
    # 診断ファイルは除外)。旧実装は全数字連結で step を作っており、拡張子の「5」や bcond id まで
    # step に混入していた (res_0.h5→5, res_outlet_2_12000.h5→2120005; 2026-08-11 レビュー指摘)。
    import re
    pat = re.compile(r'^res_(\d+)\.h5$')
    pairs = []
    for f in glob.glob(os.path.join(run_dir, 'res_*.h5')):
        m = pat.match(os.path.basename(f))
        if m:
            pairs.append((int(m.group(1)), f))
    pairs.sort()
    return [f for _, f in pairs], [s for s, _ in pairs]


def centroids(mesh):
    return h5py.File(mesh, 'r')['/CELLS/centCoords'][:].reshape(-1, 3)


def mirror_index(cc):
    try:
        from scipy.spatial import cKDTree
    except Exception:
        return None, None
    x, y = cc[:, 0], cc[:, 1]
    t = cKDTree(np.column_stack([x, y]))
    d, idx = t.query(np.column_stack([x, -y]), k=1)
    return idx, d.max()


# ---- 量の抽出子: (mesh cc, res VALUE) -> scalar ----
def q_shock(cc, V, yband=1e-3):
    x, y = cc[:, 0], cc[:, 1]
    P = V['P'][:]
    m = np.abs(y) < yband
    o = np.argsort(x[m]); xx = x[m][o]; PP = P[m][o]
    if len(xx) < 10:
        return float('nan')
    p0 = np.median(PP[xx < 0.05]) if np.any(xx < 0.05) else PP[0]
    thr = p0 + 0.15 * (PP.max() - p0)
    idx = np.where(PP > thr)[0]
    return xx[idx[0]] * 1e3 if len(idx) else float('nan')  # mm


def q_machmax(cc, V):
    return float((np.sqrt(V['Ux'][:]**2 + V['Uy'][:]**2 + V['Uz'][:]**2) / V['sonic'][:]).max())


def q_pmax(cc, V):
    return float(V['P'][:].max())


# ---- 平板専用: 運動量厚さ theta と Cf/Cf_KS (Karman-Schoenherr 外部相関) ----
# 前提: 壁は y=0 で x>0 が平板 (case/26 系)。x<0 は slip の助走区間なので除外する。
# Cf は「壁関数/壁解像を問わず同じ定義」にするため Reichardt 則の逆解きで u_tau を求める
# (壁解像 y+<1 でも Reichardt は u+=y+ に縮退するので同一定義で扱える)。
# theta は積分核が十分ゼロになる外部流まで積分する (delta99 で打ち切らない)。
KAPPA_WL = 0.41

# 平板の抽出・壁法則・運動量積分は **正式後処理と共通のモジュール**を使う
# (かつて式を複製していて fit 窓/端条件/音速が食い違い、同じ run に別の値を報告していた)。
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flatplate_bl as _fbl


def _reichardt_uplus(yp):
    return np.log(1.0 + KAPPA_WL * yp) / KAPPA_WL + 7.8 * (
        1.0 - np.exp(-yp / 11.0) - (yp / 11.0) * np.exp(-yp / 3.0))


B_WL = 5.0


def _spalding_uplus(yp):
    """Spalding は陰形式 y+=H(u+) なので u+ を Newton で求める (kappa=0.41, B=5.0)。"""
    up = np.log(max(yp, 1e-12)) / KAPPA_WL + B_WL
    for _ in range(200):
        e = KAPPA_WL * up
        f = up + np.exp(-KAPPA_WL * B_WL) * (np.exp(e) - 1 - e - e**2 / 2 - e**3 / 6) - yp
        df = 1 + np.exp(-KAPPA_WL * B_WL) * KAPPA_WL * (np.exp(e) - 1 - e - e**2 / 2)
        st = f / df
        up -= st
        if abs(st) < 1e-13:
            break
    return up


def _uplus(yp, law):
    if law == 'spalding':
        return _spalding_uplus(yp)
    return _reichardt_uplus(yp)


def _solve_utau(Ut, y, nu, law='reichardt'):
    """壁法則の逆解き (共通実装)。**ソルバが使った法則と揃えること**。"""
    return _fbl.solve_utau(Ut, y, nu, law)


cf_karman_schoenherr = _fbl.cf_karman_schoenherr


def _plate_column(cc, V, xs):
    """平板上 x=xs 最近傍の壁法線カラムを壁から昇順で返す。"""
    x, y = cc[:, 0], cc[:, 1]
    on = x > 1e-6
    if not np.any(on):
        return None
    xa = np.unique(np.round(x[on], 6))
    xc = xa[np.argmin(np.abs(xa - xs))]
    idx = np.sort(np.where(np.abs(x - xc) < 1e-4)[0])
    o = np.argsort(y[idx])
    return idx[o]


def _theta_and_utau(cc, V, xs, ytop=None, law='reichardt'):
    col = _plate_column(cc, V, xs)
    if col is None or len(col) < 5:
        return float('nan'), float('nan'), float('nan')
    yy = cc[col, 1]
    ro = V['ro'][:][col]
    ux = V['roUx'][:][col] / ro
    mu = V['vis_lam'][:][col]
    # 壁点 (y=0, u=0) を必ず含める: node は既にある / cell は補う
    if yy[0] > 1e-12:
        yy = np.concatenate(([0.0], yy))
        ux = np.concatenate(([0.0], ux))
        ro = np.concatenate((ro[:1], ro))
        mu = np.concatenate((mu[:1], mu))
        first = 1
    else:
        first = 1                      # 壁点の次が第一 DOF
    m = np.ones(len(yy), dtype=bool) if ytop is None else (yy <= ytop)
    ue = ux[m].max()
    roe = ro[m][int(np.argmax(ux[m]))]
    core = (ro[m] / roe) * (ux[m] / ue) * (1.0 - ux[m] / ue)
    theta = float(np.trapz(core, yy[m]))
    nu1 = mu[first] / ro[first]
    utau = _solve_utau(ux[first], yy[first], nu1, law)
    return theta, utau, ue


def _dstar(cc, V, xs, ytop=None):
    col = _plate_column(cc, V, xs)
    if col is None or len(col) < 5:
        return float('nan')
    yy = cc[col, 1]
    ro = V['ro'][:][col]
    ux = V['roUx'][:][col] / ro
    if yy[0] > 1e-12:
        yy = np.concatenate(([0.0], yy)); ux = np.concatenate(([0.0], ux))
        ro = np.concatenate((ro[:1], ro))
    m = np.ones(len(yy), bool) if ytop is None else (yy <= ytop)
    ue = ux[m].max(); roe = ro[m][int(np.argmax(ux[m]))]
    return float(np.trapz(1.0 - (ro[m] * ux[m]) / (roe * ue), yy[m]))


def make_q_cf_momentum(xs, ytop, window=0.08, order=2, xmin=0.1, xmax=0.95):
    """運動量積分 Cf (壁出力非依存) の時系列判定。**正式後処理と同一実装** (flatplate_bl)。"""
    def f(cc, V):
        Vv = dict(u=V['roUx'][:] / V['ro'][:], ro=V['ro'][:], mu=V['vis_lam'][:],
                  P=V['P'][:] if 'P' in V else None)
        cf, ret = _fbl.cf_momentum(cc, Vv, xs, xmin, xmax, ytop,
                                   fit_window=window, fit_order=order)
        ks = _fbl.cf_karman_schoenherr(ret)
        return float(cf / ks) if np.isfinite(ks) and ks > 0 and np.isfinite(cf) else float('nan')
    return f


def make_q_flatplate(xs, ytop, law='reichardt'):
    def q_theta(cc, V):
        return _theta_and_utau(cc, V, xs, ytop, law)[0]

    def q_cf_retheta(cc, V):
        theta, utau, ue = _theta_and_utau(cc, V, xs, ytop, law)
        col = _plate_column(cc, V, xs)
        if col is None:
            return float('nan')
        ro = V['ro'][:][col]
        mu = V['vis_lam'][:][col]
        ux = V['roUx'][:][col] / ro
        j = int(np.argmax(ux))
        roe, mue = ro[j], mu[j]
        re_theta = roe * ue * theta / mue
        cf = 2.0 * (utau / ue) ** 2
        ks = cf_karman_schoenherr(re_theta)
        return float(cf / ks) if np.isfinite(ks) and ks > 0 else float('nan')

    return q_theta, q_cf_retheta


def make_q_asym(cc):
    idx, dmax = mirror_index(cc)
    if idx is None or dmax > 1e-5:   # scipy 無し or 非対称メッシュ
        return None
    def f(cc_, V):
        M = np.sqrt(V['Ux'][:]**2 + V['Uy'][:]**2) / V['sonic'][:]
        return float(np.linalg.norm(M - M[idx]) / max(np.linalg.norm(M), 1e-30))
    return f


def classify(steps, vals, tail_frac, drift_tol, osc_tol, min_snaps):
    s = np.array(steps, float); v = np.array(vals, float)
    good = np.isfinite(v)
    s, v = s[good], v[good]
    n = len(v)
    if n < min_snaps:
        return 'TRANSIENT-UNSETTLED', f"only {n} snapshot(s) (<{min_snaps})", None
    k = max(3, int(math.ceil(tail_frac * n)))
    st, vt = s[-k:], v[-k:]
    mean = float(np.mean(vt)); scale = max(abs(mean), 1e-30)
    span = float(vt.max() - vt.min())
    fluct = span / scale
    # 末尾の線形トレンド (傾き×幅 / 平均)
    if st.max() > st.min():
        slope = np.polyfit(st, vt, 1)[0]
        drift = abs(slope * (st.max() - st.min())) / scale
    else:
        drift = 0.0
    # 全系列の極値が末尾末端にある = まだ成長/減衰中
    extremum_at_end = (np.argmax(v) >= n - 2) or (np.argmin(v) >= n - 2)
    detail = f"tail mean={mean:.4g}  drift={drift*100:.1f}%/tail  fluct={fluct*100:.1f}%"
    amp = span / 2.0
    if drift > drift_tol:
        return 'DRIFTING', detail + ("  (extremum at tail-end)" if extremum_at_end else ""), (mean, amp)
    if fluct > osc_tol:
        return 'OSCILLATING', detail + f"  -> report {mean:.4g} +/- {amp:.2g}", (mean, amp)
    if extremum_at_end and drift > drift_tol * 0.5:
        return 'TRANSIENT-UNSETTLED', detail + "  (still trending at tail-end)", (mean, amp)
    return 'STEADY', detail, (mean, amp)


SEV = {'STEADY': 0, 'OSCILLATING': 1, 'TRANSIENT-UNSETTLED': 2, 'DRIFTING': 3, 'NONFINITE': 4}


def classify_series(steps, vals, tail_frac, drift_tol, osc_tol, min_snaps):
    """classify の非有限値を **黙って落とさない** 版 (CSV 系列モード / 外部呼び出し用)。
    classify は NaN を除いて判定するので、[1,1,1,NaN] のような発散末尾を STEADY にしてしまう
    (case/46 の `steadiness` で実害: codex 指摘 2026-09-09)。ここでは非有限値が 1 つでもあれば
    NONFINITE を返し、詳細に個数を残す。戻り値は classify と同形 (verdict, detail, (mean, amp))。"""
    v = np.asarray(vals, float)
    bad = int(np.count_nonzero(~np.isfinite(v)))
    if bad:
        return 'NONFINITE', f"{bad}/{len(v)} non-finite value(s) in series", None
    return classify(steps, vals, tail_frac, drift_tol, osc_tol, min_snaps)


def analyze_series_csv(path, cols, tail_frac, drift_tol, osc_tol, min_snaps):
    """`step` 列を持つ CSV の指定列を classify_series で判定する。戻り値は worst の SEV。"""
    import csv
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        print(f"\n=== {path}  -> ERROR: empty CSV ==="); return 4
    missing = [c for c in cols if c not in rows[0]]
    if missing or 'step' not in rows[0]:
        print(f"\n=== {path}  -> ERROR: missing column(s) {missing + ([] if 'step' in rows[0] else ['step'])} "
              f"(available: {sorted(rows[0].keys())}) ==="); return 4
    steps = [float(r['step']) for r in rows]
    worst = 0; lines = []
    for c in cols:
        vals = [float(r[c]) if r[c] not in ('', 'None') else float('nan') for r in rows]
        verdict, detail, _ = classify_series(steps, vals, tail_frac, drift_tol, osc_tol, min_snaps)
        worst = max(worst, SEV[verdict])
        lines.append(f"  {c:16s}: {detail:55s} {verdict}")
    overall = [k for k, vv in SEV.items() if vv == worst][0]
    print(f"\n=== {path}  [{len(rows)} rows, steps {steps[0]:g}..{steps[-1]:g}]  -> {overall} ===")
    for l in lines:
        print(l)
    return worst
BUILTIN = ['shock', 'asym', 'machmax', 'pmax']
# 平板専用 (明示指定のときだけ有効): --quantity theta,cf_retheta
FLATPLATE = ['theta', 'cf_retheta', 'cf_momentum']
# 壁応力 (res_wall_<physID>_<step>.h5 の **twall 接線成分**を弧長 (×周長) 重みで平均)。
# **正式評価関数 case/26 tools/wall_tau_eval.py と同一定義**。twall が無い古い出力だけ
# rho*utau^2 へフォールバックする。
# **壁関数 run 専用** — wallTreatmentSST=0 (壁解像) では utau が出力されず全ゼロになるので
# NOT APPLICABLE を返す。壁解像の壁応力が要るなら接線 traction を別途算出する必要がある。
# 範囲は --wall-xmin/--wall-xmax、対象壁は --wall-phys-id で指定する (ケース固有値をハードコードしない)。
WALL = ['wall_model_tau']


def analyze(run_dir, want, tail_frac, drift_tol, osc_tol, min_snaps, mesh_arg, cf_x=0.6, cf_ytop=None,
            wall_phys_id=None, wall_xmin=None, wall_xmax=None, wall_law='reichardt',
            wall_axisym=False):
    mesh = find_mesh(run_dir, mesh_arg)
    if mesh is None:
        print(f"\n=== {run_dir}  -> NO mesh (.h5 with /CELLS/centCoords) ==="); return 3
    fs, steps = res_files(run_dir)
    if len(fs) < 2:
        print(f"\n=== {run_dir}  -> TRANSIENT-UNSETTLED (only {len(fs)} res_*.h5) ==="); return 2
    cc = centroids(mesh)
    extractors = {'shock': q_shock, 'machmax': q_machmax, 'pmax': q_pmax}
    if any(q in want for q in FLATPLATE):
        qt, qc = make_q_flatplate(cf_x, cf_ytop, wall_law)
        extractors['theta'] = qt
        extractors['cf_retheta'] = qc
        if 'cf_momentum' in want:
            extractors['cf_momentum'] = make_q_cf_momentum(cf_x, cf_ytop)
    qa = make_q_asym(cc)
    if qa:
        extractors['asym'] = qa
    # 未知の量を黙って無視すると「評価していないのに ALL STEADY」になるので必ずエラーにする。
    # ただし判定は **固定集合** に対して行う (run 依存で適用不能な asym 等を unknown 扱いしない)。
    known = set(BUILTIN) | set(FLATPLATE) | set(WALL)
    unknown = [q for q in want if q not in known]
    if unknown:
        print(f"\n=== {run_dir}  -> ERROR: unknown quantity {unknown} "
              f"(available: {sorted(known)}) ===")
        return 3
    # 適用不能な量 (非対称メッシュの asym 等) は従来どおり skip する
    quantities = [q for q in want if q in extractors]
    series = {q: [] for q in quantities}
    for f in fs:
        V = h5py.File(f, 'r')['VALUE']
        for q in quantities:
            try:
                series[q].append(extractors[q](cc, V))
            except Exception:
                series[q].append(float('nan'))
    worst = 0
    lines = []
    if 'wall_model_tau' in want:
        import re as _re
        pat_w = _re.compile(r'^res_wall_(\d+)_(\d+)\.h5$')
        by_step = {}
        for f in glob.glob(os.path.join(run_dir, 'res_wall_*.h5')):
            m = pat_w.match(os.path.basename(f))
            if not m:
                continue
            pid, st = int(m.group(1)), int(m.group(2))
            if wall_phys_id is not None and pid != wall_phys_id:
                continue
            by_step.setdefault(st, []).append(f)
        wsteps, wvals = [], []
        allzero = True
        for st in sorted(by_step):
            num = den = 0.0                    # 同一 step の複数壁面は面積(点数)重みで集約
            for f in by_step[st]:
                with h5py.File(f, 'r') as hw:
                    if 'VALUE/utau' not in hw or 'VALUE/ro' not in hw:
                        continue
                    Cw = hw['MESH/COORD'][:].reshape(-1, 3)
                    o = np.argsort(Cw[:, 0])
                    xw, yw = Cw[o, 0], Cw[o, 1]
                    # 正式評価関数と同一定義: twall の接線成分を弧長 (×周長) 重みで平均
                    # (case/26 tools/wall_tau_eval.py と揃える。旧実装の rho*utau^2 節点算術
                    #  平均は別の汎関数で、A/B の値がツール間で食い違う原因になっていた)。
                    if 'VALUE/twall_x' in hw:
                        tx = hw['VALUE/twall_x'][:][o]
                        ty = hw['VALUE/twall_y'][:][o]
                        tg = np.empty((len(xw), 2))
                        if len(xw) >= 3:
                            tg[1:-1] = np.column_stack([xw[2:] - xw[:-2], yw[2:] - yw[:-2]])
                            tg[0] = [xw[1] - xw[0], yw[1] - yw[0]]
                            tg[-1] = [xw[-1] - xw[-2], yw[-1] - yw[-2]]
                        else:
                            tg[:] = [1.0, 0.0]
                        tg /= np.maximum(np.linalg.norm(tg, axis=1, keepdims=True), 1e-30)
                        nrm = np.column_stack([-tg[:, 1], tg[:, 0]])
                        t2 = np.column_stack([tx, ty])
                        tau = np.linalg.norm(t2 - np.sum(t2 * nrm, axis=1)[:, None] * nrm, axis=1)
                    else:
                        tau = hw['VALUE/ro'][:][o] * hw['VALUE/utau'][:][o] ** 2
                    ds = np.empty(len(xw))
                    if len(xw) >= 3:
                        ds[1:-1] = 0.5 * np.hypot(xw[2:] - xw[:-2], yw[2:] - yw[:-2])
                        ds[0] = np.hypot(xw[1] - xw[0], yw[1] - yw[0])
                        ds[-1] = np.hypot(xw[-1] - xw[-2], yw[-1] - yw[-2])
                    else:
                        ds[:] = 1.0
                    ww = ds * (2.0 * np.pi * np.maximum(yw, 1e-30) if wall_axisym else 1.0)
                    mm = np.ones_like(xw, bool)
                    if wall_xmin is not None:
                        mm &= (xw >= wall_xmin)
                    if wall_xmax is not None:
                        mm &= (xw <= wall_xmax)
                    if mm.sum() == 0:
                        continue
                    if np.any(hw['VALUE/utau'][:][o][mm] > 0.0):
                        allzero = False
                    num += float((tau[mm] * ww[mm]).sum()); den += float(ww[mm].sum())
            if den > 0:
                wsteps.append(st); wvals.append(num / den)
        if allzero and wsteps:
            lines.append(f"  {'wall_model_tau':14s}: utau is all zero -> NOT APPLICABLE "
                         f"(wall-function run only; wallTreatmentSST=0 does not output utau)")
            worst = max(worst, 3)
        elif len(wsteps) < 2:
            lines.append(f"  {'wall_model_tau':14s}: only {len(wsteps)} usable res_wall_*.h5 "
                         f"-> TRANSIENT-UNSETTLED")
            worst = max(worst, 2)
        else:
            v, d, _ = classify(wsteps, wvals, tail_frac, drift_tol, osc_tol, min_snaps)
            worst = max(worst, SEV[v])
            lines.append(f"  {'wall_model_tau':14s}: {d:50s} {v}")
    for q in quantities:
        verdict, detail, _ = classify(steps, series[q], tail_frac, drift_tol, osc_tol, min_snaps)
        worst = max(worst, SEV[verdict])
        lines.append(f"  {q:9s}: {detail:55s} {verdict}")
    overall = [k for k, vv in SEV.items() if vv == worst][0]
    print(f"\n=== {run_dir}  [{len(fs)} snapshots, steps {steps[0]}..{steps[-1]}]  -> {overall} ===")
    for l in lines:
        print(l)
    if 'asym' not in quantities and 'asym' in want:
        print("  (asym skipped: mesh not up-down symmetric or scipy missing)")
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_dirs', nargs='*')
    ap.add_argument('--series-csv', default=None,
                    help='CSV 系列モード: `step` 列を持つ CSV の --series-cols を同じ classify で判定する')
    ap.add_argument('--series-cols', default=None, help='--series-csv で判定する列名 (カンマ区切り)')
    ap.add_argument('--quantity', default=','.join(BUILTIN),
                    help='comma list of: ' + ','.join(BUILTIN + FLATPLATE + WALL) +
                         ' (default all applicable; theta/cf_retheta are flat-plate specific '
                         'and assume the wall is y=0 with the plate at x>0)')
    ap.add_argument('--cf-x', type=float, default=0.6,
                    help='streamwise station [m] for theta / cf_retheta (default 0.6)')
    ap.add_argument('--wall-phys-id', type=int, default=None,
                    help='wall_model_tau: 対象壁の physID (res_wall_<physID>_*.h5)。既定は全壁を集約')
    ap.add_argument('--wall-xmin', type=float, default=None,
                    help='wall_model_tau: 集約する x の下限 [m] (前縁・淀み域を除くのに使う)')
    ap.add_argument('--wall-xmax', type=float, default=None, help='wall_model_tau: x の上限 [m]')
    ap.add_argument('--wall-axisym', action='store_true',
                    help='wall_model_tau: 面積重みに周長 2*pi*r を掛ける (軸対称。正式評価関数と揃える)')
    ap.add_argument('--cf-ytop', type=float, default=None,
                    help='upper limit [m] of the theta integral (default: full column). '
                         'Use to check sensitivity of theta to the truncation height.')
    ap.add_argument('--wall-law', default='reichardt', choices=['reichardt', 'spalding'],
                    help='cf_retheta の逆解きに使う壁法則。**ソルバが使った法則と揃えること** '
                         '(揃えないと A/B の符号すら逆に見える)。cf_momentum は壁出力非依存なので無関係')
    ap.add_argument('--mesh', default=None, help='input mesh h5 (auto-detected if omitted)')
    ap.add_argument('--tail', type=float, default=0.4, help='tail fraction of snapshots for the steadiness check')
    ap.add_argument('--drift', type=float, default=0.05, help='max fractional trend across tail for STEADY')
    ap.add_argument('--osc', type=float, default=0.10, help='max fractional fluctuation across tail for STEADY')
    ap.add_argument('--min-snaps', type=int, default=4, help='min snapshots required to judge')
    args = ap.parse_args()
    want = [q.strip() for q in args.quantity.split(',') if q.strip()]
    worst = 0
    if args.series_csv:
        if not args.series_cols:
            ap.error('--series-csv には --series-cols が必要')
        cols = [c.strip() for c in args.series_cols.split(',') if c.strip()]
        worst = max(worst, analyze_series_csv(args.series_csv, cols, args.tail, args.drift, args.osc, args.min_snaps))
    elif not args.run_dirs:
        ap.error('run_dir か --series-csv を指定すること')
    for rd in args.run_dirs:
        worst = max(worst, analyze(rd, want, args.tail, args.drift, args.osc, args.min_snaps,
                                   args.mesh, args.cf_x, args.cf_ytop,
                                   args.wall_phys_id, args.wall_xmin, args.wall_xmax,
                                   args.wall_law, args.wall_axisym))
    print(f"\nOVERALL: {'ALL STEADY' if worst == 0 else 'NOT ALL STEADY (see above)'}")
    sys.exit(0 if worst == 0 else 1)


if __name__ == '__main__':
    main()
