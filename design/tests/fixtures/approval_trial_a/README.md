# 承認の記録の試験の固定の入力 (試走 A)

`design/tests/run_approval_tests.py` (plan `tooling-design-orchestrator` §4.7.1・§5.1 #2a) が使う。

- `card_p0.yaml`・`card_p1.yaml`: 手順書の試走 A のカード (`case/90.qs_trial_a/`、sha256 は `case/71.m55_combustion_wt/_runs/qs_trial_a_evidence/` と同じ)。
  p1 は p0 に較正の 3 係数を足しただけ。試験はこれを版 2 に書き換えて使う。
- `problem_p0.yaml`・`problem_p1.yaml`: それぞれのカードから生成した問題 YAML (AWS B の `~/qs_trial/forge/case/90.qs_trial_a/` から写した)。
  違いは較正の係数と r_t から決まる値だけ。
- `run_0003_ns_p0/`: NS p0 の run の停止規則の基準 (`settle_criteria.yaml`)・`solverConfig.yaml`・`bcondConfig.yaml` と、
  `prepare_info.json` のうち基準の解決に使う欄 (Md・scale_m・x_E・physical_wall.domain・mesh の ni と nj)。場のファイルは含めない。
