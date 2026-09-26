# plan gradient-scalar-lsq-unification #6 の確認: 既定 lsq (省略) と 明示 lsq の 1 step 比較、起動エコー
import os, shutil, subprocess, sys, yaml
sys.path.insert(0, os.path.expanduser("~/forge-pgrad-new/case/09.Taylor-Green/_g0_lsq_seam"))
import gharness as G
B = os.path.expanduser("~/sglsq/forge_24365d72"); S = os.path.expanduser("~/sglsq/defsw")
os.makedirs(S, exist_ok=True)
def mk(case, tag, edit, dump):
    d = os.path.join(S, f"{case}_{tag}")
    if os.path.exists(d): shutil.rmtree(d)
    shutil.copytree(os.path.expanduser(f"~/sglsq/s1h/s1_prep/{case}"), d)
    p = os.path.join(d, "solverConfig.yaml"); c = yaml.safe_load(open(p)); edit(c)
    yaml.safe_dump(c, open(p, "w"), sort_keys=False, default_flow_style=None)
    if dump: os.environ["FORGE_DUMP_PREGATHER"] = os.path.join(d, "pregather")
    try: G.run_forge(d, expect_merge=(case == "tgv"), forge_bin=B)
    except SystemExit as e: print("run failed", d, e)
    os.environ.pop("FORGE_DUMP_PREGATHER", None)
    log = open(os.path.join(d, "forge_run.log"), errors="replace").read()
    print(d, [l for l in log.splitlines() if "scalarGradient" in l])
def exp_lsq(c): c["mesh"]["scalarGradient"] = "lsq"
def omit(c): c["mesh"].pop("scalarGradient", None)
def exp_gg(c): c["mesh"]["scalarGradient"] = "gg"
for case in ("tgv", "case48"):
    for r in "abc":
        mk(case, f"explsq_{r}", exp_lsq, False); mk(case, f"default_{r}", omit, True)
    mk(case, "expgg", exp_gg, False)
