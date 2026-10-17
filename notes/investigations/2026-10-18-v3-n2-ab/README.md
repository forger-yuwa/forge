# V3 13a の N2 の FAIL の判別の A/B (2026-10-18)

plan [`tooling-design-orchestrator`](../../../plans/active/tooling-design-orchestrator.md) §6「N2 の FAIL の扱いと判別の A/B」(事前登録 2fe8a5ad、
codex 諮問 [v3-n2-fail](../../reviews/2026-10-18-orchestrator-v3-n2-fail-diagnose-v3-n2-fail.md))。台本 [`n2_ab.py`](n2_ab.py)、結果 [`result.json`](result.json)。
手元、forge・全格子の変換・品質なし。

変えた因子はコードの版だけ (A = a36a87ff = 試走 A の p0 の版、B = b1a02511)。カード `case/90.qs_trial_a/card_p0.yaml` (v1)、`design/.venv-opt`、
`--axis-cap-frac 0.03`・`--radial-ratio-max 1.2` は同じ。NS・Euler の節点数は今の数え方 (`check_node_budget`) で数え直した。

| 量 | A (a36a87ff) | B (b1a02511) |
| --- | --- | --- |
| `r_throat` [m] | 0.10490286547916376 | 0.10490283359806168 (相対 −3.0e−7) |
| `M_knot` | 2.65 | 2.65 |
| `n_axis_inv` (生産の網) | 4800 (自動で倍にした) | 2400 |
| `ni_min` → `ni` (10 の倍数に切り上げ) | 4470 → 4470 | 4471 → 4480 |
| `nj` | 92 | 92 |
| NS の節点数 | **411,240** (旧の記録と同じ) | **412,160** (今の値と同じ) |
| Euler の節点数 | **433,590** | **434,560** |
| 第一層の表・密度の表 | (B と違う) | |

**判定 (事前登録の文言)**: A が旧の値、B が今の値を再現 → コードの版の変更で差を説明できる。差の出所は MOC の生産の網 (`n_axis_inv` 4800 vs 2400) から
r_t・壁・表が少し動き、`ni_min` が 1 増えて 10 の倍数の切り上げをまたいだこと。eabbab97 (#6e) で網の自動の倍化をやめた経路と合うが、
#6e 単独への帰属は保留 (事前登録どおり)。N2 の FAIL と元の上限 {412000, 434000} は変えない。
