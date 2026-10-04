# 諮問: farfield V3 (SERN の側方幅) の結果の解釈 (2026-09-29、条件 7)

plan: `plans/active/boundary-node-farfield-characteristic.md` §6 V3 (事前規則)、§5.1 #4 (結果記録)、#3g (V2 の残りをユーザ判断で限定付き記録して V3 へ)。
関連: `plans/active/tooling-nozzle-sern-3d.md` §5.1 R4d (slip のまま広げた幅系列、必要幅 3.42 H)。
run: AWS `case/46.sern_design/`、集計 `case/46.sern_design/V3_EVAL.txt` (`v3_farfield_eval.py`)。作成 `v3_farfield_setup.py`。

## 観測事実
- 全 run: g3、生産設定 (TP EXH/AIR、SST 低 Re、等温壁、block-DPLUR cfl 0.25) と同一。新バイナリ `build-ff` (ソース `a6ceee0b` + 出力専用の変更)、リミッタ基準値を run_0986 の自動値に固定 (L 2.136214053、ρ 0.04544068206、P 4411.497933、a 368.6412254)、20000 step・500 step 出力。side_far 以外の BC は R4d の同じ幅の run と同一 (diff 確認)。farfield の自由流 = inlet_ext と同じ外気。全 run GATES PASS、farfield の置換・HLL 退避 0、`check_convergence` はプラトー (R4d と同じ)。
- 規則 (R4d と同じ): 末尾 10000 step 平均、a = max|値 − 平均|、D = |Δ平均| + a₁ + a₂、窓条件 = 前 10000 step の平均との差 ≤ 0.1ε、未達なら +20000 延長。ε: C_T・C_T_with_shear・C_L 5e-4、C_M 5e-3。

| run | 幅・side_far | 初期場 | C_L | C_M | 窓 |
|---|---|---|---|---|---|
| run_0986 (R4d、旧バイナリ・基準値自動) | 2.50 H slip | run_0971 | 0.0625905 | −1.6441445 | OK |
| run_0992 | 2.50 H slip | run_0986 最終 | 0.0625928 | −1.6441951 | OK |
| run_0993 | 2.50 H farfield | run_0986 最終 | 0.0618214 | −1.6211579 | C_M NG (8.4e-4) |
| run_0996 | 2.50 H farfield (延長) | run_0993 最終 | 0.0618221 | −1.6211744 | OK |
| run_0994 | 3.42 H farfield | run_0988 最終 (3.42 H slip) | 0.0618801 | −1.6228598 | OK |
| run_0995 | 4.35 H farfield | run_0989 の入力場 (3.42 H slip 最終の共通部 index コピー + 追加層最近傍) | 0.0618780 | −1.6228136 | OK |
| run_0988 (R4d) | 3.42 H slip | | 0.0618810 | −1.6228715 | |
| run_0990 (R4d) | 6.19 H slip | | 0.0618870 | −1.6229828 | |

- D: slip 2.50 vs 旧 run_0986 = C_L 2.6e-5・C_M 5.5e-4 (バイナリ・基準値固定の影響小)。V3a (2.50 H slip vs farfield) = C_L 7.9e-4・C_M 0.0234 (> ε、記録のみ)。
  2.50 ff vs 3.42 ff = C_L 8.2e-5・C_M 2.2e-3、2.50 ff vs 4.35 ff = 8.3e-5・2.2e-3、3.42 ff vs 4.35 ff = 3.4e-5・7.4e-4、全量 ≤ ε。C_T・C_T_with_shear は全組 D ≤ 1.2e-5 (V3a の slip vs ff は 3.0e-5)。
  参考: 2.50 ff vs 6.19 H slip = 8.4e-5・2.2e-3、3.42 ff vs 6.19 H slip = 3.0e-5・6.2e-4。

## 期待値と出典
- §6 V3b の事前規則: 候補幅をそれより広いすべての試験済み幅 (少なくとも 1 つ) と比べ、全 4 量 D ≤ ε の最小幅を「farfield での必要幅」。結論は「試験した幅系列で許容内」に限定。

## 仮説 / 解釈案 (確定前)
- I1: 事前規則により farfield の必要側方幅 = 2.50 H (現行の生産幅)。slip では 3.42 H が要った (R4d)。
- I2: 2.50 H の slip → farfield の C_L・C_M の変化は、slip 側面の反射 (または境界が閉じていること) の効果で、farfield はそれをほぼ取り除き、広い領域の答え (3.42–6.19 H) に 2.50 H で近づけた。
- 注意点: (a) C_M の 2.50 H ff は広い側と約 1.7e-3 の系統差 (ε の 45 %、a を含まない平均差でも 1.6e-3)。(b) V2 は限定付き (V2c 判定不能、V2d-2 時間精度未確認、独立参照未収束) のまま、ユーザ判断で V3 に進んだ。(c) 3.42/4.35 H の farfield は slip の最終場から、2.50 H は slip 最終場からの継続で、初期場の履歴は幅ごとに違う (V2c で角の非一意性を見た)。(d) 反射の低減を直接測ってはいない (力の係数の幅依存だけ)。

## 問い
1. I1 を事前規則どおりの結論として確定してよいか。付けるべき限定は何か (C_M の系統差 45 %、V2 の限定、初期場の履歴、g3 のみ・M6 作動点のみ)。
2. I2 の「反射の除去」という解釈はこのデータで言えるか。言えないなら何を書くべきか。
3. 次の手: (a) 生産の side_far を farfield に切り替える提案 (ユーザ判断事項) の前に必要な確認、(b) plan を result 段の codex レビューに出して accepted にする前に足りないもの (V2 の限定付き項目をどう扱うか)、(c) R4d の格子差 (g1/g3/g4 を 2.50 H で取った) との関係。
