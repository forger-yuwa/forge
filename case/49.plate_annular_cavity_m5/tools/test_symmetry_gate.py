#!/usr/bin/env python3
"""旋回ゲートの回帰試験 (codex result-1 M8 の反例)。

旧実装は **終点 1 点**で `|末尾-基準| / |擾乱-基準|` を見ていたので、
  (a) 減衰したあと再成長する
  (b) 振動して終点だけ偶然基準に近い
のどちらも合格していた。判定を触るときは必ずこれを通すこと。

usage: python3 tools/test_symmetry_gate.py
"""
import sys
import numpy as np

sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
import check_symmetry  # noqa: E402

ok = True


def verdict(dev, decay=0.2):
    """**本番判定を直接呼ぶ** (2026-09-20 codex result m10)。
    以前は試験側に同じロジックを書き写していたので、本番の退行を検出できなかった。"""
    return check_symmetry.swirl_verdict(dev, decay)[0]


def endpoint_only(dev, decay=0.2):
    """旧実装 (終点 1 点)。反例が本当に旧実装を通ることを示すため残す。"""
    dev = np.asarray(dev, float)
    if dev[0] < 1e-6:
        return "NO-PERTURBATION"
    frac = dev[-1] / max(dev[0], 1e-30)
    if frac <= decay:
        return "SYMMETRIC"
    if frac >= 1.0:
        return "ASYMMETRIC"
    return "PARTIAL-DECAY"


def chk(name, dev, want, want_old=None):
    global ok
    got = verdict(dev)
    old = endpoint_only(dev)
    good = got == want
    ok = ok and good
    note = ""
    if want_old is not None:
        note = "   [旧実装は %s%s]" % (old, " ← 見逃していた" if old != want else "")
    print("  [%s] %-42s -> %-14s (期待 %s)%s"
          % ("OK " if good else "NG ", name, got, want, note))


def main():
    n = 30
    # (a) 減衰したあと再成長する (codex M8 の反例)
    dev = np.concatenate([np.geomspace(1.0, 0.05, 15), np.geomspace(0.05, 0.4, 15)])
    chk("減衰 -> 再成長", dev, "ASYMMETRIC", want_old=True)
    # (b) 振動していて終点だけ偶然小さい
    base = np.full(n, 0.5)
    osc = base * (1 + 0.9 * np.cos(np.arange(n) * 1.1))
    osc[-1] = 0.05
    chk("振動、終点だけ小さい", osc, "PARTIAL-DECAY", want_old=True)
    # (c) 素直な単調減衰 (合格すべき)
    chk("単調に 2 桁減衰", np.geomspace(1.0, 0.01, n), "SYMMETRIC")
    # (d) 減衰しない
    chk("横ばい", np.full(n, 1.0), "ASYMMETRIC")
    # (e) 成長
    chk("単調に成長", np.geomspace(1.0, 5.0, n), "ASYMMETRIC")
    # (f) 半端な減衰
    chk("5 割しか減衰しない", np.geomspace(1.0, 0.5, n), "PARTIAL-DECAY")
    # (g) 擾乱が入っていない
    chk("擾乱ゼロ", np.zeros(n), "NO-PERTURBATION")
    print("\nVERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
