# 諮問ブリーフ: 収束の判定の新方式の run_0012 の結果の解釈と次の手、#1c の不一致の扱い (2026-10-16)

- 依頼者: 主セッション `M6ノズルの設計`
- エスカレーション条件:
  - 3: 事前登録の比較が「未成立」。#1c は事前の条件「差 0」に不一致。
  - 7: result 段の解釈を確定する前。
  - 4: plan に無い再実行・修正へ進む前。
- plan: `plans/active/tooling-settle-convergence-and-monitor.md`
  - §6「V3(b)・V4 の実施」: 事前登録、(1) の判定の道具と細則、(1) の結果
  - §5.1 #1c (結果)・#5b・#5c
- 前の諮問:
  - `notes/reviews/2026-10-16-settle-v3b-v4-preregistration-diagnose-settle-v3b-v4.md` (事前登録の NO-GO と書き直し)
  - `notes/reviews/2026-10-16-settle-v1-and-line-height-diagnose-settle-v1-and-line-height.md` (#1c の発案)

## A. run_0012 の結果

### 観測事実

全文は plan §6 の「(1) の結果」を読むこと。要点:

#### 判定

- VERDICT: **未成立**
  - t_s = 136,000 で、事前に決めた上限 120,000 を超えた。
  - span は 136,000〜216,000 で、データ (〜200,000) に収まらない。
- 格子点ごとの経過:

| 格子点 | 全体の判定 | 理由 |
| --- | --- | --- |
| 42,000〜54,000 | NOT_SETTLED | rms_roK |
| 56,000〜92,000 | SETTLED | 保持は 96,000 で満たすはずだった |
| 94,000 の 1 点 | NOT_SETTLED | rms_roOmega の傾き \|s\|·W = 0.02283 > 許容 0.02119 (8 % 超) |
| 96,000〜 | SETTLED | 保持をやり直し、t_s 136,000 |

#### 読めた span (136,000〜200,000) の中身

- (a) t_s からの差 / 許容の最大:

| 量 | 比の最大 |
| --- | --- |
| 出口コア M | 0.00043 |
| 線 | 0.0043 |
| 流れ角 | 9.0e−5 |
| 壁の熱流 | 0.055 |

- (b) 全体の判定はすべて SETTLED。

#### 参考の記録 (判定に使わない)

- 目的の量だけの停止候補 t_s′ = 82,000 (保持は 42,000 から)。
- t_s′ から 162,000 まで超過なし (比の最大 0.071)。
- 目的の量だけなら SETTLED なのに残差・流量で落ちた格子点は 8。すべて t_s の前。

#### ほかのゲート

- 非有限 0。
- `check_convergence --segment` NOT CONVERGED (plateau、2.3〜3.0 桁)。
- 壁解像 PASS。
- 全圧の検査は既知の型。
- V5b FAIL、V5c の片振幅 FAIL。いずれも run_0010 と同じ値。

### 期待値と出典

- 事前登録 (§6、codex の採用):
  - 結果 A: span の全期間で (a)・(b) とも許容の内。
  - 結果 B: どれかが超過。
  - 未成立: t_s > 120,000 など。
  - 閾値は暫定値に固定する。実測に合わせて較正しない。
- 次の手は事前登録に書いていない。V4 (`run_0011_ns_v2`) は結果 A の後と書いた。

### 問い A

1. 未成立の後の次の手は何がよいか。候補:
   - (i) 同じ設定の自由継続をもっと長く回す (新しい run、上限 260,000 など、約 1.8 時間) → 同じ規則で再判定する。
     - 非決定性で t_s がまた変わりうる。
     - 120,000 の上限は「span がデータに収まる」ための値で、上限を延ばせば 136,000 も検証できる。
   - (ii) run_0012 を同じ格子で restart_field で継いで 216,000 以上まで延ばす。
     - 再開直後に残差が跳ねる既知の現象があり、span の中の (b) を乱す。
   - (iii) 停止規則の設計を見直す。
     - 1 つの格子点の僅かな超過 (傾き 8 %) で 38,000 step の保持がやり直しになるのは脆い。
     - 例: 保持は目的の量・流量で数え、残差は要求の時点と再評価で見る。または、残差の各検査に連続 k 点の超過を要する。
     - これは設計の変更で、事前登録の後に結果を見て決めることになる。
   - (iv) どれも採らず、未成立として記録し、小さいノズルなど別の条件に進む。
2. この run は「残差の条件が独立の情報を足すか」に対して否の材料か。
   - 目的の量だけで止めても超過が無かった。
   - 残差の条件は停止を 54,000 step 遅らせた。
3. V4 (run_0011) は今回の結果のもとで回してよいか。
   - 事前登録では結果 A の後としていた。

## B. #1c (状態不変 A/B) の不一致

### 観測事実

plan §5.1 #1c の「結果」を読むこと。要点:

- device: 3 ケースとも全腕でビット同一。
- host: 出力の D2H が device の値 (S0 とビット同一) を host の写しへ写すので、43〜58 個の配列が変わる。
- 数値の経路が host の写しを読むのは、点検した範囲で 2 箇所。
  - farfield の流束が host の `bc.bvar` の自由流を読む。device 側で書かれないので同じ値の写しになる。未試験。
  - CHT は使う前に D2H し直す。
- 変わった host の写しを NaN で埋めて続けても NaN は出ない。

### 問い B

4. 事前の合格条件「device と host の両方で差 0」に当たらないこの結果を、V1 (f′)・(g′) の合格と見なしてよいか。
   - 足りなければ、何を足すか。例: farfield のケースでの A/B、host の写しの読み手の網羅的な点検。

## 読んでよいファイル

- `plans/active/tooling-settle-convergence-and-monitor.md`
- `solver_density_cuda/tools/settle_virtual_stop.py`
- `solver_density_cuda/tools/settle.py`
- `design/forge_design/evaluate/stop_rule.py` (v2)
- `case/71.m55_combustion_wt/run_0012_ns_free_judge/settle_virtual_stop/`
  - `VERDICT.txt`
  - `grid_verdicts.csv`
  - `grid_checks.csv`
  - `span_deviation.csv`
  - `span_deviation_ts_prime.csv`
  - `virtual_stop.json`
- `solver_density_cuda/main.cpp` (`FORGE_DIAG_OUTPUT_INVARIANCE` の診断)
- `solver_density_cuda/cuda_forge/convectiveFlux_d.cu:110-130`
