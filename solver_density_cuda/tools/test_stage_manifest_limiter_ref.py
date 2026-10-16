#!/usr/bin/env python3
"""stage_manifest の区間の鍵にリミッタの基準値 ρ・P・a・L_ref を入れる (plan tooling-settle-convergence-and-monitor §4.8 の 3、
§6 V5、codex plan 段 2026-10-15 M4: limiterRoRef 1 → 10 の 2 段が 1 区間に連結されていた)。

合格 (測る前に固定):
  (1) 同じ値の継承は同じ区間 — 自動で決めた段 (forge の記録 limiter_ref.json) と、その値を固定値で継承した段は 1 区間。
  (2) どれかの基準値の変更は別の区間 — ρ・P・a・L_ref のどれか 1 つを変えると 2 区間 (codex の再現 1 → 10 を含む)。
  (3) 自動の値が分からなければ連結しない — 記録が無い・別の起動の記録 (cfg_fnv が違う)・記録と明示の食い違いは unknown、
      unknown の段は同じ設定の段とも、前後の値の分かる段とも連結しない。
  (4) 基準値を使わない段 (forge の記録 computed false・limiterScaled 0 の明示) どうしは連結する。
  (5) キーの無い旧い manifest は従来どおり (連結の規則を変えない)。
  (6) 4 つとも明示の段は記録が無くても値が分かる (status config)。
  (7) 記録の読み: %.17g の値が Python の float で同じ値に戻り、repr を YAML に書いて読み直しても同じ値。
forge・GPU は使わない (記録は forge の書式どおりに試験が書く; forge 自身の書き出しは plan §5.1 #6 の短い run で確かめた)。
"""
import json
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stage_manifest as SM  # noqa: E402

BC = "wall: {physID: 3, kind: wall}\n"
BASE = ('solver: "SLAU"\nmesh: {discretization: "node", nodeWallDirichlet: 1}\n'
        'space: {%sconvMethod: 1, limiter: 2}\n')
X = {"limiterRoRef": 0.59298852861870699, "limiterPRef": 192416.41762690683, "limiterARef": 327.13657637494532,
     "limiterRefLength": 9.0296406421583377}

fails = 0


def check(name, cond, detail=""):
    global fails
    print(("PASS " if cond else "FAIL ") + name + (("  [%s]" % detail) if detail and not cond else ""))
    fails += (not cond)


def cfg(explicit=None, extra=""):
    """space の先頭に明示の基準値を書いた設定 (explicit: {名前: 値})。"""
    s = "".join("%s: %r, " % (k, float(v)) for k, v in (explicit or {}).items()) + extra
    return BASE % s


def write_rec(run, cfg_text, vals=None, auto=None, computed=True, fnv=None, scaled=1):
    """forge の writeLimiterRefRecord と同じ形の limiter_ref.json を書く (値は %.17g)。"""
    lines = ['{', '  "schema": "forge.limiter_ref/v1",', '  "time": 1,',
             '  "cfg_fnv": "%s",' % (fnv or SM.fnv1a64(cfg_text)), '  "limiterScaled": %d,' % scaled, '  "limiterDiag": 0,',
             '  "computed": %s,' % ("true" if computed else "false")]
    if computed:
        for k in SM.LIMITER_REF_NAMES:
            lines.append('  "%s": %s,' % (k, "%.17g" % vals[k]))
        a = auto if auto is not None else {k: True for k in SM.LIMITER_REF_NAMES}
        lines.append('  "auto": {%s},' % ", ".join('"%s": %s' % (k, "true" if a[k] else "false") for k in SM.LIMITER_REF_NAMES))
    else:
        lines += ['  "%s": null,' % k for k in SM.LIMITER_REF_NAMES] + ['  "auto": null,']
    lines += ['  "venkatK": 0.050000000000000003,', '  "h_i": "sqrt(A_planar)"', '}']
    with open(os.path.join(run, SM.LIMITER_REF_FILE), "w") as f:
        f.write("\n".join(lines) + "\n")


def nseg(run, stages):
    """stages: [(tag, cfg_text, 記録を書く関数 or None (記録を消す) or 'keep')] を順に add (各段の forge の後に呼ぶ形)。"""
    sm = SM.StageManifest(run)
    for tag, c, rec in stages:
        if rec is None:
            p = os.path.join(run, SM.LIMITER_REF_FILE)
            if os.path.exists(p):
                os.remove(p)
        elif rec != "keep":
            rec()
        sm.add(tag, c, BC)
    sm.write()
    man = json.load(open(os.path.join(run, "stage_manifest.json")))
    return len(SM.segments(man)), man


TMP = tempfile.mkdtemp(prefix="sm_limref_")
try:
    def run(name):
        d = os.path.join(TMP, name)
        os.makedirs(d)
        return d

    # (1) 自動 → 同じ値を固定で継承
    r = run("inherit")
    cA, cB = cfg(), cfg(X)
    n, man = nseg(r, [("main", cA, lambda: write_rec(r, cA, X)),
                      ("cont", cB, lambda: write_rec(r, cB, X, auto={k: False for k in SM.LIMITER_REF_NAMES}))])
    check("(1) 自動の段 → その値を固定で継承した段: 1 区間", n == 1, n)
    check("(1) 段の記録 status forge・値は記録のまま (自動 / 固定の印も)",
          [s["limiter_ref"]["status"] for s in man["stages"]] == ["forge", "forge"]
          and man["stages"][0]["limiter_ref"]["values"] == X and man["stages"][0]["limiter_ref"]["auto"]["limiterRoRef"] is True
          and man["stages"][1]["limiter_ref"]["auto"]["limiterRoRef"] is False and man["manifest_version"] == 3)
    r = run("same_auto")
    n, _ = nseg(r, [("a", cA, lambda: write_rec(r, cA, X)), ("b", cA, lambda: write_rec(r, cA, X))])
    check("(1) 同じ設定・同じ自動の値の 2 段: 1 区間", n == 1, n)

    # (2) 基準値の変更
    for k in SM.LIMITER_REF_NAMES:
        r = run("change_" + k)
        Y = dict(X)
        Y[k] = X[k] * 1.0000001
        n, _ = nseg(r, [("a", cA, lambda: write_rec(r, cA, X)), ("b", cA, lambda Y=Y: write_rec(r, cA, Y))])
        check("(2) %s だけ相対 1e−7 変えた自動の値: 2 区間" % k, n == 2, n)
    # codex の再現: 他を同じにして limiterRoRef を 1 → 10 (記録あり・記録なしの両方)
    c1, c10 = cfg({"limiterRoRef": 1.0}), cfg({"limiterRoRef": 10.0})
    r = run("codex_rec")
    X1, X10 = dict(X, limiterRoRef=1.0), dict(X, limiterRoRef=10.0)
    fixed_ro = {k: k != "limiterRoRef" for k in SM.LIMITER_REF_NAMES}      # ρ だけ固定 (明示)、他は自動
    n, _ = nseg(r, [("a", c1, lambda: write_rec(r, c1, X1, auto=fixed_ro)), ("b", c10, lambda: write_rec(r, c10, X10, auto=fixed_ro))])
    check("(2) codex の再現 limiterRoRef 1 → 10 (forge の記録あり): 2 区間", n == 2, n)
    r = run("codex_norec")
    n, man = nseg(r, [("a", c1, None), ("b", c10, None)])
    check("(2) codex の再現 limiterRoRef 1 → 10 (記録なし、他は自動): 2 区間 (unknown)", n == 2
          and all(s["key"][SM.LIMITER_KEY] == SM.LIMITER_UNKNOWN for s in man["stages"]), n)
    c4a, c4b = cfg(X), cfg(dict(X, limiterRoRef=10.0))
    r = run("codex_all4")
    n, _ = nseg(r, [("a", c4a, None), ("b", c4b, None)])
    check("(2) 4 つとも明示で limiterRoRef だけ違う 2 段 (記録なし): 2 区間", n == 2, n)

    # (3) 値が分からない
    r = run("unknown_same")
    n, man = nseg(r, [("a", cA, None), ("b", cA, None)])
    check("(3) 自動・記録なしの同じ設定の 2 段: 2 区間 (連結しない)", n == 2, n)
    check("(3) その段の status unknown と理由", all(s["limiter_ref"]["status"] == "unknown" and "記録が無い" in s["limiter_ref"]["reason"]
                                                for s in man["stages"]))
    r = run("unknown_between")
    n, _ = nseg(r, [("a", cB, None), ("b", cA, None), ("c", cB, None)])
    check("(3) 値の分かる段 (4 つ明示) の間の unknown: 3 区間 (前後とも連結しない)", n == 3, n)
    r = run("stale")
    n, man = nseg(r, [("a", cA, lambda: write_rec(r, cA, X)), ("b", cfg(extra="venkatK: 0.05, "), "keep")])
    check("(3) 前の段の記録が残っている (cfg_fnv が違う): その段は unknown・2 区間",
          n == 2 and man["stages"][1]["limiter_ref"]["status"] == "unknown" and "別の起動" in man["stages"][1]["limiter_ref"]["reason"], n)
    r = run("clash")
    cR = cfg({"limiterRoRef": 0.5})
    n, man = nseg(r, [("a", cR, lambda: write_rec(r, cR, X))])
    check("(3) 明示 (0.5) と記録 (自動 0.5929…) の食い違い: unknown",
          man["stages"][0]["limiter_ref"]["status"] == "unknown" and "食い違う" in man["stages"][0]["limiter_ref"]["reason"])
    r = run("broken")
    with open(os.path.join(r, SM.LIMITER_REF_FILE), "w") as f:
        f.write('{"schema": "forge.limiter_ref/v1", "cfg_fnv": "%s", "computed": true, "limiterRoRef": 1.0}' % SM.fnv1a64(cA))
    n, man = nseg(r, [("a", cA, "keep")])
    check("(3) 値の欠けた記録: unknown", man["stages"][0]["key"][SM.LIMITER_KEY] == SM.LIMITER_UNKNOWN)

    # (4) 基準値を使わない
    c0 = cfg(extra="limiterScaled: 0, ")
    r = run("unused_rec")
    n, man = nseg(r, [("a", c0, lambda: write_rec(r, c0, computed=False, scaled=0)), ("b", c0, lambda: write_rec(r, c0, computed=False, scaled=0))])
    check("(4) forge の記録 computed false の 2 段: 1 区間 (unused)", n == 1 and man["stages"][0]["key"][SM.LIMITER_KEY] == SM.LIMITER_UNUSED, n)
    r = run("unused_cfg")
    n, man = nseg(r, [("a", c0, None), ("b", c0, None)])
    check("(4) limiterScaled 0 の明示・記録なしの 2 段: 1 区間 (unused)", n == 1 and man["stages"][1]["limiter_ref"]["status"] == "config", n)

    # (5) 旧い manifest
    old = {"manifest_version": 2, "stages": [{"tag": t, "key": SM.stage_key(cA, BC)} for t in ("a", "b")]}
    check("(5) キーの無い旧い manifest (同じ設定 2 段): 従来どおり 1 区間", len(SM.segments(old)) == 1)

    # (6) 4 つとも明示 → 記録なしでも値が分かる
    r = run("config_known")
    n, man = nseg(r, [("a", cB, None), ("b", cB, None)])
    check("(6) 4 つとも明示・記録なしの同じ設定の 2 段: 1 区間 (status config)",
          n == 1 and man["stages"][0]["limiter_ref"]["status"] == "config" and man["stages"][0]["limiter_ref"]["values"] == X, n)

    # (7) 往復の精度
    r = run("roundtrip")
    write_rec(r, cA, X)
    rec = SM.read_limiter_ref(r, cA)
    back = {k: float(rec[k]) for k in SM.LIMITER_REF_NAMES}
    import yaml
    y = yaml.safe_load("space: {" + ", ".join("%s: %r" % (k, v) for k, v in back.items()) + "}")["space"]
    check("(7) %.17g → float → repr → YAML → float で同じ値 (ビット)", back == X and y == X, (back, y))
    try:
        SM.read_limiter_ref(r, cfg(extra="venkatK: 0.05, "))
        check("(7) cfg_fnv の違う設定で読むと例外", False)
    except SM.LimiterRecordError as e:
        check("(7) cfg_fnv の違う設定で読むと例外", "別の起動" in str(e))
finally:
    shutil.rmtree(TMP, ignore_errors=True)

print("VERDICT:", "PASS" if fails == 0 else "FAIL (%d)" % fails)
sys.exit(1 if fails else 0)
