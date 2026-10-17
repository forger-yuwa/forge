# 諮問ブリーフ: §6 V1a (事前登録した比較) の FAIL の扱い — 実行役の差かコードの版の差か (2026-10-17)

- 依頼者: 主セッション `M6ノズルの設計`
- エスカレーション条件 3 (事前に書いた比較が FAIL)・4 (承認済みの手順に無い再実行の前)・原因を書く前。
- plan: `plans/active/tooling-design-orchestrator.md` §6 の「V1a の事前登録 (2026-10-17、§5.1 #6d、結果を見る前)」、§5.1 #6d
- commit: 28945473 (事前登録と、旧の台本との突き合わせで見つけた 2 つの穴の修正) + 比較の台本 (未 commit)

## 観測事実

- 実行役の側: recipe の NS の段と同じ CLI を手元で `--prepare-only` (`deltastar_loop --init-integral --ic-file run_0008_ns_b0/res_90000.h5 --ic-sha256 c64a73f1… --euler-ref run_0009_euler
  --euler-ref-res res_12000.h5 --euler-ref-sha256 2b2d02ce…`、問題 `problem_m55_p5_settle_v2.yaml`、変換器は手元の card4a のビルド)。一時の run に準備 (28 s) → `campaign.plan.confirm`。
- 比べる相手: V6 `case/71.m55_combustion_wt/run_0014_chain_v2/` (`_runs/run_0014_chain_v2/p6_chain.sh` が `run_pass0_integral(prob, eref, run, ic_from=run_0008_ns_b0, stop_rule=True)`、
  AWS B、git_head f20cded9、FP64 のビルド)。
- 比較の台本: `case/71.m55_combustion_wt/_runs/v1a_executor/compare.py`、結果 `result.json`。
- 結果 (事前登録の 9 項目):

| 項目 | 結果 | 差 |
| --- | --- | --- |
| (1) 全段の起動直前の設定の FNV (3 段) | 一致 | — |
| (2) bcondConfig.yaml・(3) solverConfig_prestage.yaml | バイト一致 | — |
| (4) settle_criteria.yaml | FAIL | `source.generator_sha256` (許した) と `source.written` (書いた時刻; **事前登録の正規化に書き漏らした**) |
| (5) settle_nodes.txt | FAIL | 節点の番号の列が違う |
| (6) 格子 (nozzle.h5 の座標とトポロジ・nozzle.msh の sha256) | FAIL | 座標とトポロジのハッシュ・msh が違う (ni 4810・nj 92 は同じ) |
| (7) 準備の記録の主な欄 | FAIL | initializer の数値 (delta_r_exit 0.41354049 vs 0.41356168、delta_r_throat −0.0024278873 vs −0.0024278710、平滑の残差、2 係数の CONTUR の値) |
| (8) 初期場の固定 (ic_snapshot の sha256・step 90000) | 一致 | — |
| (9) バイナリ・変換器の sha256 の記録 | 記録あり | 値の一致は求めない (事前登録) |

- 形状の量: x_E 40.633257 (実行役) vs 40.630797 (V6)、L_c 40.302422 vs 40.299962 (差 0.00246 r_t)、physical_wall.domain の終わり 79.805673 vs 79.809687、x0・x_A・scale_m は同じ。
- `git log f20cded9..HEAD -- design/forge_design/geometry/` は eabbab97 (#6e: 終端 C⁺ を x_E を挟む 2 節点で補間、長さ合わせの固定点の位置が変わる; ユーザ了承「全部の壁がビット単位で変わります」に「ぜんぜんいいよ」) の 1 つだけ。

## 仮説

- (5)〜(7) の差は実行役でなく、V6 の後の形状のコードの変更 (#6e) による壁の違い (x_E・L_c が 0.0025 r_t 動く → 壁 → 格子 → 節点の表・初期壁の数値)。(1)〜(3) の設定の文字列は形状に依らないので一致する。
- (4) の `source.written` は時刻で、事前登録の書き漏らし。

## 案 (まだ回していない)

- 判別の A/B: 同じ今のコード・同じ手元のビルドで、A = 旧の台本の呼び方 (`run_pass0_integral(prob, eref, run, ic_from=<run_0008 の写しのディレクトリ>, stop_rule=True, prepare_only=True)`)、
  B = 実行役の呼び方 (済み) を準備だけで比べる。A と B で (1)〜(7) が一致すれば (時刻の欄を除く)、差はコードの版によるもので実行役は同じ入力を作る。A と V6 の差が (5)〜(7) と同じなら #6e の支持。
- V1a の判定の書き方: 事前登録の文言どおり V1a は「FAIL」と記録し、原因の切り分け (A/B) と、V6 との比較は今のコードで作り直した基準 (A) に対して行う、を追記する (事前登録の正規化は広げない)。

## 問い

1. この A/B で切り分けてよいか。足りないもの (例: #6e を戻したコードで A を回す、AWS B の V6 と同じビルドで回す) はあるか。
2. V1a の合否の扱い: 事前登録の文言どおり FAIL にし、A/B の結果で「実行役の入力の決定性」を別に判定するのでよいか。`source.written` の書き漏らしの扱い。
3. 旧の台本との対応表 (plan §5.1 #6d) で見つけた 2 つの穴 (準備の後のファイルの変更・初期壁の係数) の直し方 (28945473) に問題はあるか。

## 読んでよいファイル

- `plans/active/tooling-design-orchestrator.md` §6 (V1a の事前登録)・§5.1 #6d
- `case/71.m55_combustion_wt/_runs/v1a_executor/compare.py`・`result.json`、`case/71.m55_combustion_wt/_runs/run_0014_chain_v2/p6_chain.sh`
- `case/71.m55_combustion_wt/_runs/run_0011_ns_v2/p5_ns.py`・`p5_all.sh`
- `design/forge_design/feedback/deltastar_loop.py` (`run_pass0_integral`)、`design/forge_design/campaign/plan.py`、`design/forge_design/campaign/recipes/wt_card_chain.py`
- `git show eabbab97`
