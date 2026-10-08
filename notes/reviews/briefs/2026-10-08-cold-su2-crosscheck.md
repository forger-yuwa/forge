# 諮問: 冷却ノズルの SU2 との照合 (CPG) の事前登録の点検

日付 2026-10-08。諮問先 codex (diagnose)。エスカレーション条件 1 (plan §6 の検証計画を新規に書く) と 4 (承認済みの手順に無い再実行の前)。
plan: `plans/active/tooling-nozzle-isothermal-wall-chain.md` §5.1 #21〜#26、§6 V-c45 (全文を読むこと)。作業ツリー `/home/sano/work/forge-integ-1005` (commit ba14bb2a)。

## 背景 (観測事実)
- 冷却壁の対 (run_0181 断熱 / run_0182・延長 run_0183 300 K、燃焼ガス TP、FP64、冷却壁用の格子 57 万節点、dilatationCorrection 2) で、δ_E の比 R が 0.726 (x = 40) 〜 0.793 (x = 94)。
- CONTUR はどの形でも 0.98〜1.05 で、冷却の効果を再現しない (V-c45 は判定不能; 記録は §6)。
- 抽出の座標の A/B では R の差が 0.06 % (#21、棄却)。
- C_f の比は NS 1.49 vs CONTUR 1.44〜1.51 で合う (#23)。
- θ の比は NS 2.85 → 2.53 vs CONTUR 約 1.9。δ/θ の比は NS 0.23〜0.32 vs CONTUR 0.52 (#24)。
- 係数を 1 つ動かしても、CONTUR は NS の x 方向の分布を再現しない (#25)。
- ユーザ: 「forge が間違ってる説はある？」「SU2 の計算も同時に進めて。SU2 と比較するためにね。dilatation 補正は本来は入れるのが正しいと思う」。
- SU2 v8.5 は燃焼ガスの多成分を forge と同じ形で解けない。CPG の設計 (`cfd_gas: cpg`) は、CFD ピンの凍結源が TP なので拒否された。そこで **壁・格子は TP の対と同じまま、CFD だけを CPG** にする。
- 前回の断熱 CPG の照合 (case/45 run_0011/0012/0013、別格子・FP32): forge の素の SST と SU2 で δ99 ≤ 3 %・θ ≤ 0.4 %・δ* ≤ 1.3 %。生産の SST (dilatationCorrection 2) は SU2 比で δ99 −17〜21 %・θ −16 %・δ* −5 %。

## 事前登録の案 (plan §5.1 #26)
- run:
  - forge (FP64、AWS): CPG × {断熱, 300 K} × {素の SST (dilat 0・KL 0), 生産の SST (dilat 2・KL 1)} の 4 本。準備済みの TP の run の設定を CPG に書き換える (thermalMethod 0・viscMethod 1・cp 1360・γ 1.27354)。初期値は収束した TP の解の ρ・U・P・k・ω を CPG の保存量に組み直し、段階起動 full。
  - SU2 (手元の PC): CPG × {断熱, 300 K} の 2 本。SST V2003m・ROE + MUSCL・Sutherland 273/111。格子は同じ msh を 17 桁のまま su2 形式に変換する。
- 比べる量: SU2 の解を同じ節点で forge の形式に写して比べる (Euler の参照が要らない量にした)。
  - (a) 共通の帯の外縁 y_b(x) (TP の断熱 run_0181 の抽出から) での δ_loc・θ_r
  - (b) 壁の C_f・q_w
  - (c) x = 40・70・94 の断面の ρu・T
- 判定: forge の素の SST と SU2 で、試験部 [40, 94] で次の 3 つを満たせば「一致」。超えたら forge 側の原因を調べる。
  - R_loc の差 ≤ 2 %
  - 各腕の δ_loc・θ_r の差 ≤ 3 %
  - C_f の差 ≤ 3 %
- 前提: NaN なし・RISING/DIVERGED なし・比べる量が 5 枚で STEADY。満たさなければ両方を同じ回数だけ延長する。

## 問い
1. 「壁・格子は TP のまま、CFD だけ CPG」に替えて SU2 と比べることで、「forge の冷却の効果は解き方の誤りか」という問いに答えられるか。TP と CPG の違いで失うものは何か。
2. 比べる量を δ_E (Euler の参照つき) から、共通の y_b の δ_loc・θ_r に替えたのは妥当か。y_b を断熱の TP run から取ることの問題 (冷却で境界層の厚さが変わる、CPG で外縁が変わる) はないか。
3. 許容差 (R 2 %・δ/θ 3 %・C_f 3 %) は、前回の照合の実績と冷却の効果 (約 25 %) から見て適切か。
4. 初期値 (TP の収束解から CPG に組み直し) と段階起動で問題ないか。SU2 は自由流の一様な初期値から始めるので、比べる前の収束・定常の確認で足りるか。
5. 生産の SST (dilat 2) と素の SST (dilat 0) の両方を回す必要はあるか (ユーザは dilat 補正は入れるのが正しいと考えている)。

## 読んでよいもの
- 上記 plan、`case/45.isobutane_m6_d155/cold_pair.py` (`theta_diag`)、`run_0012_su2_sst/sst.cfg`、`compare_bl_su2.py`
- `procedures/su2-cross-check.md`
- `notes/reviews/2026-10-08-cold-pair-result-diagnose.md`
