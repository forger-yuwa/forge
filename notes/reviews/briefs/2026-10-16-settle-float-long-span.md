# 諮問ブリーフ: float の v3 の長い区間の検証 (#14) を回すか、回すならその事前登録 (2026-10-16)

- 依頼者: 主セッション `M6ノズルの設計`
- エスカレーション条件:
  - 4: plan §5.1 #14 の新しい run (AWS B、約 30 分) に進む前。
  - 1: §6 に事前登録を書く前。
- plan: `plans/active/tooling-settle-convergence-and-monitor.md`
  - §5.1 #8 (v3 の再生の結果)・#14・#15、§6 の「V3(b)・V4 の実施」(FP64 の事前登録と結果)、§4.4 (2026-10-16 のユーザ決定: float、残差と節点ごとの壁の熱流を監視から外す)

## 観測事実

- FP64 の自由継続 (停止規則なし) の再生 (v3、`case/71.m55_combustion_wt/_runs/settle_evidence/v3_replay/`):
  - run_0012・0013 とも **A**。t_s 32,000、span 32,000〜112,000 で t_s からの変化 / 許容の最大は、出口コア M 0.027・線 0.018・流れ角 0.0016。
- float の run_0015・0016 (停止規則 v2 で回したが、#8 の時点で 60,000 まで; 判定の記録は `case/71…/run_0015_ns_float_a/`・`run_0016_ns_float_b/`):
  - v3b (壁の熱流も外した) の再生で t_s 32,000 (FP64 と同じ)。
  - 32,000〜60,000 の変化 / 許容の最大は、出口コア M 0.069/0.050・線 0.14/0.12・流れ角 0.008/0.007。
  - span の終わり 112,000 までデータが無いので **未成立**。
- float と FP64 の比較 (card plan V8): 60,000 の全場で目的の量の差は 6e−5 % 以下、δ_E の差 ≤ 0.04 mm。
- 手順書の試走 A (別の入口条件・寸法、float、v3) は p0・p1 とも REACHED (本段 42,060・30,088 step)。停止の後の区間は回していない。

## 期待値と出典

- V3 の結果 A の条件 (FP64 の事前登録と同じ形): t_s から t_s + 80,000 まで、目的の量の t_s からの差が固定の許容の内、各格子点で全条件が SETTLED。
- ユーザは float に決めた (2026-10-16「ひとまずは、設計を早く進めるのが第一だから、floatで良い」)。

## 案 (主セッション)

- run: `case/71.m55_combustion_wt/run_0017_ns_float_free130`。run_0015 と同じ入力 (`problem_m55_p5_settle_v2.yaml`、float のビルド `~/forge-cardmesh-float`、初期場・Euler の参照も同じ) で、停止規則なし・本段 130,000 step (約 33 分)。
- 判定: `settle_virtual_stop.py --criteria <v3b の基準>` で再生 → t_s と span (t_s〜t_s + 80,000) の A/B/未成立。細則は V3(b) の (1) と同じ。t_s > 50,000 → 未成立。
- 閾値は v3 の暫定値に固定する (較正しない)。

## 問い

1. #14 は回す価値があるか。それとも FP64 の A と、float と FP64 の 60,000 までの一致で足りるとして、#14 を「回さない (理由)」で閉じてよいか。
   - 判断の材料: float の 32,000〜60,000 の変化 / 許容は FP64 の 2〜8 倍 (出口コア M 0.069 対 0.027、線 0.14 対 0.018)。ただし許容の内。
2. 回すなら、上の事前登録 (本段 130,000、t_s の上限 50,000、閾値の固定) に穴はないか。
3. run_0015 (停止規則 v2 で回った run) を restart で延ばす案は、再開の直後の跳ねがあるので採らない。この判断でよいか。

## 読んでよいファイル

- `plans/active/tooling-settle-convergence-and-monitor.md`
- `case/71.m55_combustion_wt/_runs/settle_evidence/` (README・v3_replay・v3b_replay)
- `case/71.m55_combustion_wt/run_0015_ns_float_a/`・`run_0016_ns_float_b/` (判定の記録)
- `solver_density_cuda/tools/settle_virtual_stop.py`
- `case/71.m55_combustion_wt/_runs/run_0011_ns_v2/p5_all.sh`・`p5_ns.py`
