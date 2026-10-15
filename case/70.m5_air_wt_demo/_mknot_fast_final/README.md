# M_knot 高速化の最終の版の証拠 (2026-10-15、AWS B)

plan [`tooling-nozzle-mknot-fast`](../../../plans/active/tooling-nozzle-mknot-fast.md) §5.1 #10・§6.6 (result 段 NO-GO の M1 への対応)。
最終の版 = feature/mknot-fast **cba2914f** (逆 MOC の単位過程の numba の版、Newton の修正の後)。比較の相手 = 引き継ぎ時点 **9468e0d0** (numpy の版)。

- 機械: AWS B (g5.xlarge、AMD EPYC 7R32、4 vCPU) — `env.txt`。Python の環境は自分の venv (`~/mknot_venv`: numpy 2.5.2・scipy 1.18.0・numba 0.68.0・llvmlite 0.50.0)。
- 台本: `run_final.sh` (全部を順に)、`time_eval.py` (N 2400 の 1 評価の計時)、`mem_one.py` (1 評価のメモリのピーク)、`check_card_v2.py` (カードの出力の照合)。
  B の作業 `~/mknot_time/` に、`base/` (9468e0d0 の `git archive`)・`final/` (cba2914f の `git archive`) を置いて回した。BLAS は 1 スレッド。

| ファイル | 中身 |
| --- | --- |
| `time_base.log` / `time_final.log` | N 2400 の 1 評価 (a) case/70 c2b (M_K 2.05)・(b) case/71 (r_t 0.1 m、M_K 2.80) を各 3 回。`numpy`/`numba` の行は `FORGE_MOC_KERNEL` (9468e0d0 には numba の経路が無いので両方 numpy) |
| `mem.txt` | 新しいプロセスで 1 評価したときの最大 RSS (N 600 は L_total_tol 0.1 r_t) |
| `card_a/` | M5 c2b のカード (`../card_m5_air_demo_c2b.yaml`) の生成: `gen.log` (標準出力と `/usr/bin/time -v`)・`problem.yaml`・card_log・Euler の写し・`check.txt` (照合) |
| `card_b/` | M5.5 の最終のカード (`card_m55_final.yaml` = M6ノズルの設計の feature/card-local-mesh 94790833 の `case/71.m55_combustion_wt/card_m55.yaml` の写し。`overrides: {mesh.ni: 3000}` は格子の見積もり (節点の上限) を避けるためだけで、上限は緩めていない。寸法・M_knot は ni に依らない) の生成: 同上 |
| `card_times.txt` | カードの生成の壁時計 |

結果の要約は plan §6.6。
