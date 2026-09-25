# 引き継ぎ: `convection-slau-wall-normal-chi` (2026-09-24 〜 25)

**残作業の正本は [`plans/active/convection-slau-wall-normal-chi.md`](../../plans/active/convection-slau-wall-normal-chi.md) §5.1 の表。**
本文書は写しとポインタだけ (`notes/sessions/` は使い捨て。恒久的な残作業をここに置かない)。

- ブランチ `feature/sern-design`、本文書のコミット `bb2b075d` 時点で push 済み・未コミット無し。
- CHT セッションが `feature/cht-phase2-fem2d` をこのブランチへマージ済み (`60bd88f`)。

## いま何が残っているか

**受入セット 6 項目のうち 5 つ完了。残るは #10 (case/16) だけ。**
**#8 の codex `--stage result` 3 回目は #10 が揃ってから**回す (揃う前に出すと同じ M1 で NO-GO になり 1 周捨てる)。

| # | 状態 |
| --- | --- |
| ~~7~~ 周期・軸対称 | 済 (全項目 PASS) |
| ~~9~~ 格子感度 3 水準 | 済 |
| ~~10b~~ SERN 2D 受入 | 済 (**PASS**) |
| ~~10c~~ 固定点 | 済 (**CFL 依存だが許容の 1/5 以下**。「同一固定点」とは書かない) |
| ~~10d~~ 実カーネル照合 | 済 (P0–P5 全て PASS) |
| **10** case/16 の V2・V3 | **条件確定・実行が残っている** |
| 11 既定化の判断 (F) | #10 完了 + $C_L$ 0.433 % / $C_M$ 0.316 % / 衝撃足 0.842 % の正否を独立に判定してから |

## #10 を再開する手順 (plan §6 V3 の case/16 条件に全部書いてある)

1. `mesh/nozzle_user_2d_planar.geo` を gmsh → **no-slip 壁付き bcond で `convertGmshToForge`**
   (2 引数: `<in.msh> <out.h5>`)。41,640 節点・`check_mesh_quality` PASS を確認済み。
2. **既製レシピ `case/16.nozzle_wys/run_user_profile.py`** を `--disc node --phys sst --laminar`、
   出口を `outflow` に `--bc-sub` して **flag 0 を 1 本**起動する。**自前の段構成を組まない**。
   - 既知の躓き: `run_user_profile.py:174` が存在しない `run_0050_*/probe.yaml` を copy する → 1 行修正が要る。
3. **前提ゲート**: flag 0 の壁 $p/p_0$ が既知値 (x=16.4: 0.367 / 45.6: 0.255 / 85: 0.187) と ±2 % 以内、
   出口中心線 $M \gtrsim 1.7$。**これを通すまで V3 を判定しない**。
4. flag 0 の収束場から **2 本を分岐** (別々に起動しない)。`slauWallNormalChi: 1` は `space:` 配下、
   起動ログの `'slauWallNormalChi' in 'space': 1` で入ったことを確認する。
5. V3 判定 → V2 (ノイズ床、両側 3 本以上)。

**期待値は「差 ≈ 0」**。壁クラスタでは壁隣接内点の $M \ll \sqrt2$ なので $\chi_n \approx \chi$。
**V3 は「害が無い」の回帰で、効果の検証ではない** (効果側は V1/V5)。差 0 を「実装が効いていない」と読まないこと。

## この 2 日で確立した道具 (次回そのまま使える)

- **AWS を自分で起こせる**: `solver_density_cuda/tools/aws_instance.sh {status|ip|start|stop|ssh}` (26 秒)。
  鍵は `/home/sano/.aws-wsl/` (700/600)。**`/mnt/c` 配下は chmod が効かず 777 になる**。
  `config` だけ `[profile forge]`、`credentials` は `[forge]`。skill `/aws-gpu` に罠をまとめた。
  **インスタンスが止まるのは正常** (`idle_autostop.sh`: GPU 0 % + forge 無し + ログイン無しが 30 分)。回避しない。
  **インスタンスは 1 台を複数セッションで共有している**。「前回自動停止したから止まっている」と決めつけない
  (別セッションが起こして run を回していることがある)。触る前に `status`、running なら `busy`
  (`forge=… logins=… gpu=…`) で使用状況を見る。`stop` は使用中なら拒否する (2026-09-25 追加)。
  既に running でも自分が起こしたことにせず、他の run と GPU を取り合う投入をしない。
- **面流束で測る**: `FORGE_DUMP_MASSFLUX=<path>` で第 1 評価の `massflux` と
  **カーネルが読んだ状態** (`<path>.state`) を書く。場は 1 step でもビット再現しない (`atomicAdd`)。
  **`res_*.h5` の `P` は流束が読む `P` と 1 ulp 違う** (EOS の往復が流束前に 2 回走る) ので `.state` を使う。
- 分析スクリプト: `case/46.sern_design/{v3sern_series,v3sern_single,v3sern_ct}.py`、
  `case/46.sern_design/cad/{v7_dist,v7_closure,v7_cos}.py`。
- **力係数は自前で組まない**: `forge_design.metrics.sern_forces.force_history` を使う (`twall` の cell/node 規約込み)。

## 次のセッションが踏みやすい穴 (この 2 日で全部踏んだ)

1. **`initial:` キーは run 時に効かない**。`main.cpp:1193` が h5 の VALUE を読むだけで、
   `setInitial()` は**変換時** (`convertGmshToForge.cpp:59`) にしか呼ばれない。
2. **一様亜音速 IC + `outflow` はノズルを起動できない** (「$u=0$・$P=P_t$ 全域」が両 BC を満たす定常解)。
3. **`check_quasisteady.py` の `--drift`/`--osc` は割合**、`fluct` は**系列自身の平均に対する比**。
   **差の系列に「許容の 1/5」を渡すと数百倍厳しくなる** (#10b で踏んだ)。量は run の量に当てる。
4. **既定閾値の判定と §6 閾値の判定を混ぜない** (「ALL STEADY」が既定 5 %/10 % の話だった)。
5. **dump の一覧は `sort -n`** (文字列ソートは 9000 > 27000)。
6. `pgrep -f "<文字列>"` は**自分のコマンドラインに一致する**。until ループが永久に待つ。
7. **収束判定は系譜連結で**。単体 run は既に下がった場から再開するので桁が小さく出る。

## 関連文書

- plan §6 V5 (実カーネル照合) / V6 (周期・軸対称) / V7 (固定点) — **判定基準はすべて測る前に固定してある**。
- `procedures/cloud-aws-gpu.md` (起動・鍵・idle 停止) / `procedures/development-environment.md` (レジスタ上限)。
- case README: `case/46.sern_design/README.md` (run_0445–0449, `_v3sern/*_ext`)、
  `case/16.nozzle_wys/README.md` (run_0471 の失敗 2 回)。
