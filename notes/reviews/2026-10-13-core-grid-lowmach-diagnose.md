# codex 諮問 (diagnose): core-grid-lowmach

- **brief**: [`notes/reviews/briefs/2026-10-13-core-grid-lowmach.md`](../../notes/reviews/briefs/2026-10-13-core-grid-lowmach.md)
- **plan**: [`plans/active/tooling-nozzle-core-grid.md`](../../plans/active/tooling-nozzle-core-grid.md)
- **date**: 2026-10-13
- **commit**: `ba069ed6` (feature/nozzle-core-grid)
- **codex**: effort `high`, 2.5 min, rc=0
- **結論**: **(a) 同一の停滞場から、B0 と `lowMachPrecond=1` だけを変えた各 20,000 step の A/B を行う。CFL の腕追加と Γ の実装は、この試験では行わない。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **(a) 同一の停滞場から、B0 と `lowMachPrecond=1` だけを変えた各 20,000 step の A/B を行う。CFL の腕追加と Γ の実装は、この試験では行わない。**

ユーザへの説明: 「今回の結果だけでは、低マッハ補正が欠けたピースという見立ては否定できません。温度差が主に育っていても、補正が変える質量流束は温度を運ぶエネルギー流束にも効くためです。ただし、圧力の市松が停滞を起こしているという証拠もまだなく、補正だけで改善するかを直接比較する段階です。」

第 1 仮説: **B0 の停滞には圧力差に対する質量流束の結合不足が関与し、既存の `=1` 単独でも実用上有意に緩和できる。** 確度: **低**

根拠: 停滞した Gref の流れ方向の P の a/q は 1.9e−4、上限 1・Grefc は約 7.8e−6。ただし、これは因果関係を示さない観測上の関連である（[plan:1147](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1147)）。コードでは `c_diss` が `mdot` を変え、その `mdot` が連続・運動量・エネルギー流束に共通して使われる。したがって「差がエントロピー主体」は、この仮説の反証にならない（[SLAU:587](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:587)、[SLAU:617](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:617)）。

反証条件: 補正の作動を確認した適格な比較で、下記の改善基準を明確に満たさないこと。発散した場合も「現在の B0 に `=1` を足すだけで改善できる」は棄却するが、**低マッハ機構そのものは棄却しない**。

第 2・第 3 仮説:
- 第 2: 圧力結合を変える効果はあるが、既存の左辺・擬似時間更新では安定して利用できない。確度: 中。ただし B0 では未確認。backstep の `=1` は CFL 1 で発散、0.3 で圧力振幅を減らしても残差が停滞した（[調査記録:106](/home/sano/work/forge-coregrid/notes/investigations/backstep-lowmach-checkerboard.md:106)）。
- 第 3: 圧力の短波長成分は主に停滞の結果であり、支配要因は別の反復上の結合にある。確度: 低・未確認。

判別 A/B:

- **起点:** `case/45.isobutane_m6_d155/run_0608_cg_gref_cont_k1/res_60000.h5`。同一メッシュの保存量を `restart_field.py` で引き継ぎ、新規 run 2 本の入力を一致させる。
- **唯一の変更:** `time.deltaT.lowMachPrecond: 0 / 1`。CFL 4、上限 50、`implicitRelax 0.7`、`nStepInner 5`、ITJ 5、`precondEps 0.15` は共通。正しいキー階層は [solverConfig.cpp:473](/home/sano/work/forge-coregrid/solver_density_cuda/input/solverConfig.cpp:473)。
- **長さ:** 各 20,000 step。序盤の非有限・床の発動を確認し、発散時は打ち切る。判定窓は **10k〜15k と 15k〜20k の両方**。
- **見る量:** 同じ島の体積重み付き残差ノルム（ρ・ρE）、全保存量の残差、両方向の P の a/q・η・f。島外への残差の移動も確認する。双子間の差は今回の判定量にしない。
- **資格:** 対照が停滞を再現し、補正側の圧力差係数が実際に変わっていること。両腕の `check_convergence.py` の VERDICT と、報告する派生量の `check_quasisteady.py` の VERDICT を添える。

**結果 A:** 両判定窓で補正側の島の残差が対照の **1/3 以下（ρ・ρE とも）**となり、他の保存量や島外に悪化を移していない → 第 1 仮説を支持。「`=1` 単独では有意に改善しない」を退ける。ただし、変更した離散方程式の反復改善であり、元の解の正確さや Δx 依存の機構を証明したとはしない。

**結果 B:** 両判定窓で両残差の比が **0.8 以上**、または補正側が発散 → 第 1 仮説の「`=1` 単独で緩和できる」を退ける。特に P の a/q が半減しても残差が改善しなければ、「測った圧力成分を抑えるだけで停滞が解ける」を退けられる。

閾値は今回提案する事前登録値である。中間・指標間の不一致・資格不足は判定不能とし、A/B の二択に押し込まない。

やらない方がよいこと: 補正側だけ CFL を下げること、失敗後に同じ試験の閾値を緩めること、`=1` の失敗を理由に低マッハ仮説全体を除外すること。Γ のライン実装や調査終了も、今回の観測だけからは選ばない。

呼び出し側の前提への異議:

- **Major — エントロピー主体から補正の無効性は導けない。** 成分比は「差の内訳」であり「差を駆動する結合」の同定ではない。根拠は上記 `mdot` の共用。対案は、成分解析を保持したまま直接 A/B を行うこと。
- **Major — backstep の失敗は重要だが、本ケースでの失敗も成功も確定しない。** 同記録の `=1` は圧力振幅を減らしながら停滞している。対案は今回の問いを「現在の B0 に `=1` だけを足して改善するか」に限定し、失敗時に Γ 必須とも低マッハ無関係とも断定しないこと。
- **Major — 「時間項は解を変えないから振動を直せない」は、反復中の振動には適用できない。** [theory.md:237](/home/sano/work/forge-coregrid/methods/time_integration/theory.md:237) の説明では、離散定常解と反復軌道を区別する必要がある。固定点を保存しても、その安定性は変わり得る。対案は「定常解に残る空間的市松」と「擬似時間反復の増幅」を分けること。
- **Minor — a は一般には市松の振幅そのものではない。** [cg_lowmach_saw.py:33](/home/sano/work/forge-coregrid/case/45.isobutane_m6_d155/cg_lowmach_saw.py:33) は添字の 4 階差分で、非交互の波や滑らかな曲率変化にも応答する。対案は §4.79 の「判定不能」を維持し、約 24 倍を「同じ物理波長の不安定モードが 24 倍」と解釈しないこと。

不足情報: 生データ・実効 config・バイナリ・VERDICT 原本は、今回の読み取り制限に従い独立検証していない。数値は指定 plan の記録に依拠する。ファイル変更・forge 起動は行っていない。**plan 未反映**。採用する場合は、呼び出し側が `plans/active/tooling-nozzle-core-grid.md` の §4.79 に続く事前登録と §5.1 に反映する。
