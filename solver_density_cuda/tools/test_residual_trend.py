"""`check_convergence.py` の rising 判定 (末端スパイク + 有意な緩い上昇) の回帰試験。

codex plan レビュー (2026-09-19, tooling-nozzle-sern-startup Major 1) が出した反例を含む。
判定を触るときは必ずこれを通すこと。リポジトリルートから `python3 solver_density_cuda/tools/test_residual_trend.py`。
"""
import sys, math, numpy as np
sys.path.insert(0,'solver_density_cuda/tools')
import check_convergence as cc
TF=0.2
def trend(ser):
    n=len(ser); a=ser[int(n*(1-TF)):]; b=ser[int(n*(1-2*TF)):int(n*(1-TF))] or a
    ma=sum(abs(x) for x in a)/len(a); mb=sum(abs(x) for x in b)/len(b)
    nz=[abs(x) for x in ser if x!=0.0]; smin=min(nz) if nz else 0.0
    r=cc._tail_rise(ser,TF); oj = True if r is None else (r[0] or r[1])
    return 'rising' if (ma>mb*1.05 and ma>2.0*smin and oj) else ('flat' if ma>mb*0.9 else 'falling'), r
ok=True
def chk(name, ser, want):
    global ok
    t,r=trend(list(ser))
    good = t==want
    ok = ok and good
    print(f'  [{"OK " if good else "NG "}] {name:52s} -> {t:8s} (期待 {want:8s}) spike/slow={r}')

# (a) 実 run: pintle の非有限になる直前まで (codex 反例)
_,cols=cc.load_series('case/37.pintle_nozzle/run_0009_slau_outflow_full/residual_history.csv')
s=np.array(cols['rms_roUx']); k=int(np.argmax(~np.isfinite(s))) if (~np.isfinite(s)).any() else len(s)
chk('pintle run_0009 rms_roUx (Inf 直前で打ち切り)', s[:k], 'rising')
# (b) 上昇 + 振動 (codex の合成反例): 末尾 400 点で 0.3 桁上昇、±0.4 桁の交互振動
n=1000; base=np.linspace(0,0.0,n); base[-400:]=np.linspace(0,0.306,400)
osc=0.4*((-1.0)**np.arange(n)); chk('合成: 0.3 桁上昇 + ±0.4 桁交互振動', 10**(base+osc)*1e-6, 'rising')
# (c) 定常振動のみ (上昇なし)
chk('合成: ±0.4 桁交互振動のみ (上昇なし)', 10**osc*1e-6, 'flat')
# (d) 単調低下
chk('合成: 単調低下 3 桁', 10**np.linspace(0,-3,1000)*1e-3, 'falling')
# (e) 実 run: 同一形状の重複ペア (R-h の元の症状)
for run,want in [('run_0193_rb_ref_A','flat'),('run_0195_rb_ref_B','flat')]:
    _,c=cc.load_series(f'case/46.sern_design/{run}/residual_history.csv')
    chk(f'{run} rms_roY1 (同一設定の重複ペア)', c['rms_roY1'], want)
# (f) 実 run: 数 step で爆発
_,c=cc.load_series('case/46.sern_design/run_0008_smoke_sst_node/residual_history.csv')
chk('run_0008 rms_roe (4 step で発散)', c['rms_roe'], 'rising')
print('\n全体:', 'PASS' if ok else 'FAIL')
raise SystemExit(0 if ok else 1)
