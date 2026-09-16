#!/usr/bin/env python3
"""確定場 (最終 res_<int>.h5) の受動種ゲート: 0 ≤ roXi/ro ≤ 1、モーメント ≥0、実現可能性 (solver と同じ条件; passive_gate_common.check_field)。
使い方: check_passive_field.py RUN_DIR   (VERDICT PASS/FAIL, exit 0/1)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from passive_gate_common import load_config, check_field

if __name__ == '__main__':
    run = sys.argv[1]
    cfg = load_config(run)
    ok, probs = check_field(run, cfg)
    for p in probs: print('  -', p)
    print('VERDICT:', 'PASS' if ok else 'FAIL'); sys.exit(0 if ok else 1)
