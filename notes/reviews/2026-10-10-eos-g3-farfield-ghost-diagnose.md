# codex 諮問 (diagnose): eos-g3-farfield-ghost

- **brief**: [`notes/reviews/briefs/2026-10-11-eos-g3-farfield-ghost.md`](../../notes/reviews/briefs/2026-10-11-eos-g3-farfield-ghost.md)
- **plan**: [`plans/active/tooling-sern-te-wake-grid.md`](../../plans/active/tooling-sern-te-wake-grid.md)
- **date**: 2026-10-10
- **commit**: `4935db8b` (feature/sern-design)
- **codex**: effort `high`, 2.5 min, rc=0
- **結論**: **v2 を緩めず、farfield ghost の EOS 入力を有限で整合した状態にした共通再生入力で、旧版／新版の EOS 1 回比較を行う。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 重大度 | 判断対象 | 採否・根拠・対案 |
|---|---|---|
| **Major** | 現在の g3 比較を合格とする | **却下。`INVALID` を維持する。** 5 組とも差分 0 バイトだが、`roe/T` の pre と `roe/T/P/Ht/sonic/gamma/cp` の post に各 51,143 点の非有限がある（[比較記録:12](/home/sano/work/forge-sern-design/notes/investigations/2026-10-08-eos-dump/g3/w/cmp_old_r1_vs_new_r1.txt:12)、[同:42](/home/sano/work/forge-sern-design/notes/investigations/2026-10-08-eos-dump/g3/w/cmp_old_r1_vs_new_r1.txt:42)）。EOS は ghost を含め `roe/T` を実際に読む（[dependentVariables_d.cu:76](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/dependentVariables_d.cu:76)、同:93・126）。前回の書き込み専用配列の除外とは異なる。**有限で整合した ghost を持つ共通入力で再生比較する。** |
| **Major** | farfield ghost の問題を #8 に含める | **採用。ただし同一原因の確定は保留。** 初期化 EOS は境界条件より先に走り（[main.cpp:1846](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1846)、同:1857）、farfield の状態境界処理は何もしない（[boundaryCond.cpp:498](/home/sano/work/forge-sern-design/solver_density_cuda/boundaryCond.cpp:498)）。#8 を「ghost の熱力学状態・物性の初期化と参照順序」に広げ、farfield の `roe/T` 残留を子項目にする。 |
| **Major** | `GATES PASS` を理由に #4 へ進む | **現時点では却下。** g3 比較の合格が時間積分の前提として明記されている（[plan:93](/home/sano/work/forge-sern-design/plans/active/tooling-sern-te-wake-grid.md:93)）。また farfield 流束が ghost を読まなくても、EOS・`gasProperties` は読む（[gasProperties_d.cu:84](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/gasProperties_d.cu:84)、同:111・121）。**まず下記 A/B を完了する。合格しても #8 の無害性を証明したことにはしない。** |

結論: **v2 を緩めず、farfield ghost の EOS 入力を有限で整合した状態にした共通再生入力で、旧版／新版の EOS 1 回比較を行う。**

第 1 仮説: **初期化 EOS が不完全な ghost 状態から非有限を作り、farfield では主要状態を再充填しないため残留している。カウンタ変更による新規発生ではない。** 確度: **高。ただし最初の発生時点は未確認。**

  根拠:

- 保存量と組成の入力読み込みは実節点まで（[variables.cpp:801](/home/sano/work/forge-sern-design/solver_density_cuda/variables.cpp:801)、同:863）。EOS は ghost まで処理する。
- 組成が全成分 0 なら、正規化後も 0 のまま（[dependentVariables_d.cu:105](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/dependentVariables_d.cu:105)）。温度反転の研磨では `cp = R = 0` に対して分母の保護値も 0 になる（[thermo_d.cuh:742](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/thermo_d.cuh:742)）。これは具体的な NaN 発生候補である。
- 保存記録では、新版の初期化入力の非有限 ghost は **0**、密度床 ghost は **165,659**。step 1 の EOS 入力では非有限 ghost が **51,143** に増えている（[SUMMARY.json:46](/home/sano/work/forge-sern-design/notes/investigations/2026-10-08-eos-dump/g3/w/SUMMARY.json:46)、同:80）。旧版にも同数の非有限がある。
- farfield の `roY` は後から Neumann コピーされる一方、主要状態と RANS ghost は充填されない（[speciesTransport_d.cu:886](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/speciesTransport_d.cu:886)、[ransBoundary_d.cu:254](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/ransBoundary_d.cu:254)）。**「farfield ghost は何も埋まらない」ではなく、配列ごとに処理が異なる。**

  反証条件: 初期化 EOS 直後には対象 ghost が有限で、別処理が初めて非有限を生成している、または非有限 index 集合が farfield ghost 集合と一致しない。

第 2・第 3 仮説: **非有限演算が有限入力時の旧新差を隠している可能性。確度: 低・未確認。** 現在の差分 0 は支持材料ではないが、有限入力での非干渉まで保証しない。第 3 仮説は置かない。

判別 A/B: **変えるのは EOS の実装だけ。A = 旧版、B = 新版カウンタ有効。各腕 2 回、時間更新 0 回、EOS 1 回。**

- 共通入力の準備として、`BCONDS/10` の接続から ghost と owner を対応付け、**非有限 index 集合との完全一致**を確認する。件数一致だけでは足りない。
- 対象 ghost の読む入力一式 `ro, roU*, roe, roK, roOmega, roY*, T` を、対応する有限な owner の状態で一括充填する。実節点・他境界・物性 DB・設定は変更しない。`roe/T` だけのゼロ置換はしない。
- その同一ファイルを両版へ再生し、再ピンや別 EOS を挟まず評価する。再生位置は既存実装にある（[eosDump.cpp:119](/home/sano/work/forge-sern-design/solver_density_cuda/eosDump.cpp:119)）。
- 判定は既存 v2 のまま、**読む入力の有限性・同一性、全書き込み先 post の有限性、差分 0 バイト、未変更入力の不変、腕内再現性**を要求する。

→ **結果 A：全条件成立**なら、この有限化した g3 入力における EOS の数値干渉仮説を退ける。元の `INVALID` は保存する。  
→ **結果 B：有限・同一入力から再現する旧新差**なら、第 2 仮説を支持し、非干渉の受入れを止める。非有限が残る場合は試験不成立であり、旧新差とは解釈しない。

やらない方がよいこと: **ghost を有限性判定から外す、NaN のビット一致を合格に読み替える、先に通常 BC や EOS を修正して比較対象まで変えること。** 再生入力の変更と製品コードの修正を混ぜない。

呼び出し側の前提への異議:

- **「床事象 0」は実節点に限定する必要がある。** 初期化の ghost 密度床は 165,659 件ある。
- **`GATES PASS` は収束・ghost 無害性の証拠ではない。** `case/46.sern_design/run_1079_tewake_B10_m10/` の保存判定は **`NOT CONVERGED (stalled/plateau)`**（[判定の写し:2](/home/sano/work/forge-sern-design/notes/investigations/2026-10-08-te-wake-ab/run_1079_tewake_B10_m10/CONVERGENCE_VERDICT.txt:2)）。力量は別判定で **`OVERALL: ALL STEADY`**（[準定常判定](/home/sano/work/forge-sern-design/notes/investigations/2026-10-08-te-wake-ab/run_1079_tewake_B10_m10/QUASISTEADY_forces.txt)）。この区別を維持する。

不足情報: ローカルの写しには元 EOS ダンプ・対象メッシュがなく、**全非有限点が farfield ghost であるという index 単位の照合は独立再検証できなかった**。初期化直後の状態も未記録。今回確認したのは保存比較記録とコードである。**ファイル変更なし・plan 未反映**。呼び出し側で §5.1 #2・#8、§6.1 に採否と追加再生入力の適用範囲を記録する。
