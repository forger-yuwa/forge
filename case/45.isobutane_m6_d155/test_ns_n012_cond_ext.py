"""ns_n012_cond_ext.residual_gate の試験 (2026-10-07 result 段レビュー M1: 判定不能・列の欠落・RISING を合格にしない)。
usage: python3 test_ns_n012_cond_ext.py (FAIL 件数 0 で終了コード 0)"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ns_n012_cond_ext as X  # noqa: E402

FAIL = 0


def check(label, cond):
    global FAIL
    print(("ok   " if cond else "FAIL ") + label)
    FAIL += 0 if cond else 1


HDR = ["step", "inner", "phase", "rms_ro", "rms_roUx", "rms_roUz", "rms_rog_0", "rms_dq_ro"]
OK = """=== /tmp/x  [last step 37999]  -> NOT CONVERGED (stalled/plateau — needs scheme change, not more steps) ===
  rms_ro      : init=1 fin=1 drop= 0.1dec flat     <-- STALLED (plateau)
  rms_roUx    : init=1 fin=1 drop= 0.1dec flat     <-- STALLED (plateau)
  rms_roUz    : all-zero (inactive, skip)
  rms_rog_0   : init=1 fin=1 drop= 0.1dec flat     <-- STALLED (plateau)
"""
r = X.residual_gate(OK, HDR)
check("plateau・全列あり → 合格 (rms_roUz の inactive は許す、rms_dq_* は要らない)", r["ok"] and r["status"] == "plateau" and not r["missing_cols"])
r = X.residual_gate(OK + "  必須列の欠損: 判定不能\n", HDR)
check("判定不能の行がある → 不合格", not r["ok"] and r["status"] == "undeterminable")
r = X.residual_gate("\n".join(l for l in OK.splitlines() if "rog_0" not in l), HDR)
check("凝縮のモーメントの列の行が無い → 不合格", not r["ok"] and r["missing_cols"] == ["rms_rog_0"])
r = X.residual_gate(OK.replace("rms_rog_0   : init=1 fin=1 drop= 0.1dec flat     <-- STALLED (plateau)",
                               "rms_rog_0   : init=1 fin=1 drop= 0.1dec rising   <-- RISING (divergent)"), HDR)
check("RISING の列 → 不合格", not r["ok"] and r["status"] == "rising")
r = X.residual_gate("", HDR)
check("出力が空 → 不合格", not r["ok"])
r = X.residual_gate(OK.replace("rms_rog_0   : init=1 fin=1 drop= 0.1dec flat     <-- STALLED (plateau)", "rms_rog_0   : all-zero (inactive, skip)"), HDR)
check("rms_roUz 以外の inactive → 不合格", not r["ok"] and "rms_rog_0" in r["missing_cols"])
print(f"FAIL 件数: {FAIL}")
sys.exit(1 if FAIL else 0)
