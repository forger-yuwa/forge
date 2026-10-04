# codex 諮問 (diagnose): farfield-4f-r9-result

- **brief**: [`notes/reviews/briefs/2026-10-02-farfield-4f-r9-result.md`](../../notes/reviews/briefs/2026-10-02-farfield-4f-r9-result.md)
- **plan**: [`plans/active/boundary-node-farfield-characteristic.md`](../../plans/active/boundary-node-farfield-characteristic.md)
- **date**: 2026-10-02
- **commit**: `c5c5c36e` (feature/sern-design)
- **codex**: effort `high`, 4.5 min, rc=0
- **結論**: **限定運用を維持し、同じ #14 後バイナリ・同じ g4 初期場から `physProp.ljSource` だけを変える A/B を1組行い、持ち越しの未検証部分を確かめる。**
- **extra**: `plans/active/tooling-nozzle-sern-chain.md`, `case/46.sern_design/README.md`, `case/46.sern_design/v3_farfield_eval.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

対象の `run_1008`–`run_1027` は手元にありません。以下の数値は **plan・run 台帳の報告値**として扱い、設定仕様・評価器・コード差分を照合しました。数値の PASS は独立再認定していません。

| 判断 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| #14 前の新輸送で、必要幅 2.50 H・格子差込みで §8 内 | **限定付き採用** | `run_1017/1021/1022/1023` の報告値では幅比較は全量 D≤ε、最大の C_M の \|G\|＋D は 0.0188＜0.05。「当該形状・m6_on・固定基準値・試験した幅と g3/g4 間の感度指標が許容内」とする。**2 水準の差は真値に対する格子誤差の上限ではない**。[farfield plan:160](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:160) |
| `run_1027` 1 本で #14 後にも幅・G の結論をそのまま持ち越す | **要再検証・Major** | C_M の D＝5.8e−4≤τ は **g3・2.50 H の更新前後比較**。g4・広幅の応答は未測定。各条件の更新差を δ とすると、G の変化は δ_g4−δ_g3、幅間差の変化は δ_広幅−δ_2.50。δ_g3 だけでは双方を拘束できない。**事前規則を満たした事実は保持するが、規則の証明範囲を訂正する**。まず下記の g4 対照を取る。[R9:433](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-chain.md:433) |
| R9・R9b は力係数を「ノイズ床以下でしか動かさない」、B 群のばらつき増大は放置可 | **一般化を却下・Major** | N は A・B の群内最大差から作るため、B のばらつきが増えると判定も緩くなる。C_M の群内差は 3.6e−4→1.5e−3、約4.2倍。ただし R9 の平均差 4.3e−5 は旧 A だけから作る N_A＝1.08e−3 よりも小さく、**平均差の PASS 自体を撤回する根拠はない**。「登録 N 以下の群平均差」と「再現性幅の増大・原因未確認」を併記する。3D の τ はノイズ床ではない。各反復の平均・振幅・前後窓差を回収し、過渡混入と反復間差を分ける。[R9:433](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-chain.md:433) |
| G がほぼ不変なのは輸送切替が g3・g4 に同程度に効いたため | **記述的解釈は採用、因果断定は保留・Minor** | 報告平均から、旧→新の変化は g3/g4 で ΔC_L＝＋1.400e−4/＋1.096e−4、ΔC_M＝−3.474e−3/−2.692e−3。**同符号・同程度だが、g4 は g3 の約78%**。その差が ΔG＝−3.04e−5/＋7.82e−4 になる。「報告値は同方向の応答と整合」と書けるが、g4 の輸送単独 A/B がないため機序は未同定とする。[farfield plan:158](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:158)、[run 台帳:361](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:361) |
| farfield を accepted にし、生産設定の限定運用も解除できる | **却下・Major** | V2c は必要条件未達、V2d-2 は時間精度差が1.13%→2.18%で未達、独立参照も精度不足。#4f はこれらを相殺しない。farfield の受入れには残る検証・帳簿実装・docs・result レビューが必要。生産全体では R4d の他方向の領域感度・同時拡大確認も残る。[farfield plan:149](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:149)、[未完了実装:162](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:162)、[R4d:1942](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1942) |

結論: **限定運用を維持し、同じ #14 後バイナリ・同じ g4 初期場から `physProp.ljSource` だけを変える A/B を1組行い、持ち越しの未検証部分を確かめる。**

第 1 仮説: **LJ 更新による4係数の変化は、g4・2.50 H でも τ 以下に収まる。** 確度: **中**

  根拠: 報告上、EXH の σ・ε/k の変化は約0.04%、AMB は不変で、g3 の比較は全量 D≤τ。コードでも lump の LJ は構成種の質量分率平均、拡散係数はその解決済み LJ を読む。このため小さい応答は妥当な予測だが、g4 での上限保証にはならない。[speciesDB.cpp:509](/home/sano/work/forge-sern-design/solver_density_cuda/input/speciesDB.cpp:509)、[thermo_d.cuh:571](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/thermo_d.cuh:571)

  反証条件: 下記 A/B が必要条件を満たしたうえで、いずれかの係数で **|平均_B−平均_A|−a_A−a_B＞τ** になること。

第 2 仮説: **LJ 更新の応答が格子解像度に依存し、g4 では τ を超える。** 確度: **低・未確認**。現在の g3 1 点では除外できない。

判別 A/B:

- 起点は `case/46.sern_design/run_1021_ff4f_g4_2p50_cont20k/` の最終場。同一格子の保存量を `restart_field.py` で別々の新規 run にコピーする。
- バイナリは両方 `f28ca2fa` の同一実体。**変更点は生成後の `solverConfig.yaml` の `physProp.ljSource` だけ**。A＝`[legacy_v1]`、B＝`[gri30, svehla1962]`。BC・輸送指定・CFL・リミッタ基準値・初期保存量を照合する。旧 LJ の再現指定は仕様に存在する。[solver-settings.md:260](/home/sano/work/forge-sern-design/procedures/solver-settings.md:260)
- 各20000 step、500 step間隔。前後10000 stepの既存窓条件を満たさなければ各＋20000、なお未達なら判定不能。4係数の平均・a・窓差、残差と準定常性の VERDICT、NaN・床・farfield 置換回数を残す。
- **小差側**：全量 D(A,B)≤τ → g4 で τ 超の応答という第2仮説を、この測定精度で棄却する。A が旧 `run_1021` も D≤τ で再現することを確認する。
- **大差側**：いずれかで |Δ平均|−a_A−a_B＞τ → 第1仮説を棄却し、#14 後の系列で G・幅判定を取り直す。両者の中間は判定不能。
- **小差側でも広幅への持ち越しは未証明**。この1組は g4 の疑義を解消する試験であり、幅系列全体の代用にはしない。

やらない方がよいこと: 事前登録済みという理由で持ち越し規則の妥当性まで保証されたと扱うこと、N を使って再現性悪化を無視すること、`GATES PASS` を残差収束と言い換えること。ゲートは `require_residual_pass` が無効なら残差 PASS を必須にしない。[sern_gates.py:292](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:292)

呼び出し側の前提への異議: **観測された小差と、未測定条件でも小差であるという推論が混ざっている。** 前回の道筋は基本的に維持する。#4e の輸送対照と #14 前の #4f は進展として認め、今回の持ち越し検証を挟んで、V2 未達・独立参照・残る実装 → result レビューへ進む。限定運用解除の条件はまだ満たしていない。

不足情報: 対象 run の実効 config、バイナリ完全ハッシュ、初期場照合記録、`force_history.csv`、全残差列、判定区間付き `CONVERGENCE_VERDICT.txt`・準定常性の原本、メッシュ品質・置換回数の記録。R9/R9b は反復ごとの値が必要で、群平均と最大差だけではばらつき増大を診断できない。

**plan 未反映。** 指示どおりファイル変更・forge 起動は行っていない。反映先は farfield plan §5.1 #4f・§6.1、SERN chain §5.1 R9。
