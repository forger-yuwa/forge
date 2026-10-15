# 諮問ブリーフ: 収束の判定の新方式の呼び出し側の移行 (#7) と通しの検証 (V6) の設計 (2026-10-16)

- 依頼者: 主セッション `M6ノズルの設計`
- エスカレーション条件: 1 (plan §4・§6 を書いた)。§5.1 #7 は担当 F。
- plan: `plans/active/tooling-settle-convergence-and-monitor.md`
  - §4.9 (新規): 呼び出し側の移行の 1〜5 と V6
  - 関連: §4.7・§4.8 の 2、§6 (1″)・(2′) の結果
- ユーザの指示: 「ガンガン進めて」 (2026-10-16)
- 実装は並行して implementer に出している。指摘は取り込む。

## 観測事実 (ここまで)

- `wt_nozzle_v2` (区間 1 万・保持なし) の再生 (§6 (1″))
  - run_0012・0013 とも A。
  - t_s 46,000。span で全量が許容の 0.1 以下、全条件 SETTLED。
- V4 (§6 (2′)): `run_0011_ns_v2` で、次のとおり条件 (i)〜(vi) すべて合格。
  - 46,000 で要求し、46,202 で協調停止、REACHED。
  - オンラインと再生の判定が一致した。
  - 再評価は SETTLED。
  - δ_E を 1 回抽出した (71.6 s)。
- 今の受け渡しの関門 (`deltastar_loop.handoff_gate`) は、v2 の run を必ず止める。
  - 旧方式の (b) `check_convergence --segment`・(c) 折れ線の流量・(e) δ の 2 窓を要求するため。

## 設計 (§4.9)

- 新方式の関門で置き換えるもの:
  - (a) STOP_RULE_STATE の v2 REACHED
  - (b)(c) 停止の step の再評価の残差・流量 SETTLED
  - (e) 停止規則の δ_E をそのまま使う。帯の感度は記録だけ。
- 残すもの: (d)・(f)・(g)・(h)。
- 呼び出し側の変更:
  - `run_case.sh`: 本段で基準があれば `check_settled.py` を回す。
  - Stop フック: 両方の出力を受け付ける。
  - AGENTS.md: 判定の正本を 2 つにする。
- 雛形の既定を `wt_nozzle_v2` にする。
- V6: case/71 を `deltastar_loop.run_pass0_integral` で 1 pass 回す。新方式の関門 → `delta_r_next.csv` → 成果物。

## 問い

1. §4.9 の移行に穴はあるか。特に:
   - (a) 停止の step の再評価の残差・流量を、関門の (b)(c) の代わりにしてよいか。区間 1 万の判定は、旧の「低下の桁数」とは別の量を見ている。
   - (b) δ_E を停止の step の 1 回で受け渡すときの注意。帯の感度 0.164 mm は、前に「設計の細かい誤差の根拠にしない」と決めた大きさ。
   - (c) 前の段 (S1・S2) の `run_case.sh` の後処理を今の道具のままにする分岐の正しさ。
   - (d) Stop フックのセッション別の claim との関係。
   - (e) AGENTS.md の正本を 2 つにする書き方。
2. 雛形の既定を `wt_nozzle_v2` にしてよいか。閾値は暫定で、確かめたのは case/71 の 1 条件、同じ条件の 2 軌道だけ。雛形から作る他の風洞ノズルにも効く。
3. V6 の合格条件と、run の設計に穴はあるか。

## 読んでよいファイル

- `plans/active/tooling-settle-convergence-and-monitor.md`
- `design/forge_design/feedback/deltastar_loop.py`
- `design/forge_design/evaluate/stop_rule.py`
- `solver_density_cuda/tools/run_case.sh`
- `solver_density_cuda/tools/hook_convergence_gate.py`
- `solver_density_cuda/tools/settle.py`
- `design/forge_design/evaluate/settle_presets.py`
- `AGENTS.md`
