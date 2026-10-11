"""受理検査の版と残差の入力の検査済み要約の試験用の値 (plan tooling-sern-te-wake-grid §5.1 R2k)。

学習の採否 (driver_sern._learnable) は 2026-10-12 から、各作動点に現行の受理検査の版 (sern_gates.ACCEPTANCE_VERSION) と
残差の入力の検査済み要約 (sern_gates.residual_evidence_summary) の合格を要求する。学習の他の条件 (処置・評価方式・床事象・判定区間・
scalarGradient・flag_policy) を試す既存の試験は、手で組んだ作動点要約・偽の collect の gates に「現行のゲートで判定した」記録を足して、
試験の意図 (他の条件による選別) を保つ。値は residual_health が完全な残差 (outer_end が 0..N−1、最終場 res_N) に返す形と同じ。
"""
from forge_design.metrics import sern_gates as G

N_STEPS = 20000
FLUID = ("rms_ro", "rms_roUx", "rms_roUy", "rms_roUz", "rms_roe")
SST = ("rms_roK", "rms_roOmega")


def accepted_residual(verdict: str = "NOT CONVERGED (stalled/plateau)", n_steps: int = N_STEPS, columns=FLUID + SST) -> dict:
    """完全な残差に residual_health が返す結果 (受理・延長と台帳の要約が読む欄)。"""
    cols = list(columns)
    return {"ok": True, "converged": verdict.startswith("PASS"), "verdict": verdict, "nan": False, "rising": [],
            "stalled": cols if "stalled" in verdict else [], "converging": [],
            "columns": {c: {"msg": "", "ok": False, "trend": "flat", "drop": 1.0} for c in cols},
            "last_step": str(n_steps - 1), "unverifiable": False, "input": {"ok": True, "reasons": [], "notes": []},
            "required_columns": cols, "missing_columns": [],
            "series": {"n_steps": n_steps, "n_rows": n_steps, "n_outer_end": n_steps, "first": 0, "last": n_steps - 1,
                       "final_field": f"res_{n_steps}.h5"},
            "reasons": []}


def accepted_gate_fields(**kw) -> dict:
    """偽の collect の gates に足す欄: 受理検査の版と完全な残差の判定。"""
    return {G.ACCEPTANCE_VERSION_KEY: G.ACCEPTANCE_VERSION, "residual": accepted_residual(**kw)}


def accepted_op_fields(**kw) -> dict:
    """手で組んだ作動点要約に足す欄: 受理検査の版と残差の入力の検査済み要約。"""
    return {G.ACCEPTANCE_VERSION_KEY: G.ACCEPTANCE_VERSION, "residual_evidence": G.residual_evidence_summary(accepted_residual(**kw))}
