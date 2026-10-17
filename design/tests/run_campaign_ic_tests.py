#!/usr/bin/env python3
"""初期場の固定・継承する値の解決・適合性 (`forge_design.campaign.ic`、`runner_axismach.select_ic`) の試験
(plan tooling-design-orchestrator §4.1.1 の 1、§5.1 #4a)。手元の試走 A (`case/90.qs_trial_a`、git 追跡外) を読む (書かない)。

  (P1) NS p0 (run_0003) の res_42060.h5 を NS p1 の初期場として解決: 元は NS・適合・ファイルの sha256、継承するリミッタの基準値は
       明示の値 {ro, p, a, L} になり、実際に run_0005_ns_p1 が使った値 (solverConfig の space) と同じ
  (P2) select_ic(ic_from) は今までどおり最大番号の場を選び、選んだファイルを記録する
  (N1) 計画を確定した後に元の run に新しい場 (番号の大きい res) を足す → 固定したファイル (ic_file・ic_sha256) は変わらない。
       ic_from (ディレクトリ) だけだと新しい場を拾う (直す前の振る舞い)
  (N2) 固定したファイルの中身が計画の後に変わる (sha256 が違う) → 止める
  (N3) ic_file と ic_from の run が違う → 止める
  (N4) 適合しない: エネルギー基準 (thermo_href_temp 300)・熱力学のモデル (cfd_gas cpg)・設計マッハ数 (M_design 5.0) → 止める

実行: /home/sano/work/forge/design/.venv-opt/bin/python design/tests/run_campaign_ic_tests.py
所要: 数秒。試走 A が無ければ飛ばす (rc 0 と「飛ばした」)。
"""
import os
import shutil
import sys
import tempfile
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from forge_design.campaign import ic as IC  # noqa: E402
from forge_design.evaluate import runner_axismach as RA  # noqa: E402

QA = Path(os.environ.get("FORGE_QS_TRIAL_A", HERE.parents[1] / "case" / "90.qs_trial_a"))
R3, R5, PROB1 = QA / "run_0003_ns_p0", QA / "run_0005_ns_p1", QA / "problem_p1.yaml"
FAIL = 0


def check(name, cond, info=""):
    global FAIL
    print(("ok  " if cond else "FAIL") + " " + name + (f"  [{info}]" if info else ""))
    if not cond:
        FAIL += 1


def raises(fn, needle):
    try:
        fn()
    except ValueError as e:
        return needle in str(e)
    return False


need = [R3 / "res_42060.h5", R3 / "prepare_info.json", R3 / "solverConfig.yaml", R5 / "solverConfig.yaml", PROB1]
if any(not p.exists() for p in need):
    print(f"飛ばした: 試走 A が無い ({[str(p) for p in need if not p.exists()]})")
    sys.exit(0)

print("== 正例")
r = IC.resolve_ic(R3 / "res_42060.h5", PROB1, for_ns=True)
used = (yaml.safe_load((R5 / "solverConfig.yaml").read_text()).get("space") or {})
want = {"ro": used.get("limiterRoRef"), "p": used.get("limiterPRef"), "a": used.get("limiterARef"), "L": used.get("limiterRefLength")}
check("(P1) 元は NS・適合・ファイルの step と sha256", r["kind"] == "ns" and r["compat"] == "OK" and r["snapshot"]["step"] == 42060
      and len(r["snapshot"]["sha256"]) == 64, str({k: r[k] for k in ("kind", "compat")}))
check("(P1) 継承するリミッタの基準値が明示の値 {ro, p, a, L} になり、run_0005 が実際に使った値と同じ",
      isinstance(r["limiter_ref"], dict) and r["limiter_ref"] == {k: float(v) for k, v in want.items()},
      f"{r['limiter_ref']} / {want}")
f, rec = RA.select_ic(R3)
check("(P2) select_ic(ic_from) は最大番号の場を選び、選んだファイルを記録する", f.name == "res_42060.h5" and rec["how"].startswith("ic_from")
      and rec["sha256"] == r["snapshot"]["sha256"])

print("== 負例")
with tempfile.TemporaryDirectory() as tdn:
    td = Path(tdn)
    src = td / "run_src"
    src.mkdir()
    for fn in ("prepare_info.json", "solverConfig.yaml", "bcondConfig.yaml", "limiter_ref.json", "stage_manifest.json", "MAIN_STAGE_STARTED"):
        if (R3 / fn).exists():
            shutil.copy2(R3 / fn, src / fn)
    (src / "res_42060.h5").symlink_to(R3 / "res_42060.h5")
    pinned = IC.resolve_ic(src / "res_42060.h5", PROB1, for_ns=True)["snapshot"]
    (src / "res_50000.h5").write_bytes(b"newer field")
    f1, _ = RA.select_ic(None, src / "res_42060.h5", pinned["sha256"])
    f2, _ = RA.select_ic(src)
    check("(N1) 新しい場を足しても固定したファイルは変わらない、ic_from だけだと新しい場を拾う (直す前の振る舞い)",
          f1.name == "res_42060.h5" and f2.name == "res_50000.h5", f"{f1.name} / {f2.name}")
    (src / "res_42060.h5").unlink()
    (src / "res_42060.h5").write_bytes(b"changed")
    check("(N2) 固定したファイルの中身が計画の後に変わる → 止める",
          raises(lambda: RA.select_ic(None, src / "res_42060.h5", pinned["sha256"]), "sha256"))
    check("(N3) ic_file と ic_from の run が違う → 止める", raises(lambda: RA.select_ic(R3, src / "res_42060.h5", None), "run が違う"))
    base = yaml.safe_load(PROB1.read_text())
    for nm, edit, needle in (("エネルギー基準 (thermo_href_temp 300)", lambda d: d["evaluate"].update(thermo_href_temp=300.0), "エネルギー基準"),
                             ("熱力学のモデル (cfd_gas cpg)", lambda d: d["evaluate"].update(cfd_gas="cpg"), "熱力学のモデル"),
                             ("設計マッハ数 (M_design 5.0)", lambda d: d["spec"].update(M_design=5.0), "設計マッハ数")):
        d = yaml.safe_load(PROB1.read_text())
        edit(d)
        pp = td / "problem_edit.yaml"
        pp.write_text(yaml.safe_dump(d, sort_keys=False, allow_unicode=True))
        check(f"(N4) 適合しない: {nm} → 止める", raises(lambda: IC.resolve_ic(R3 / "res_42060.h5", pp, for_ns=True), needle))
    del base

print(f"\n{'ALL PASS' if FAIL == 0 else f'{FAIL} FAIL'}")
sys.exit(1 if FAIL else 0)
