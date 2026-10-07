# codex 諮問 (diagnose): m6-finemesh-divergence

- **brief**: [`notes/reviews/briefs/2026-10-05-m6-finemesh-divergence.md`](../../notes/reviews/briefs/2026-10-05-m6-finemesh-divergence.md)
- **plan**: [`plans/active/tooling-nozzle-cfd-pinned-initial-line.md`](../../plans/active/tooling-nozzle-cfd-pinned-initial-line.md)
- **date**: 2026-10-05
- **commit**: `754f02d8` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.3 min, rc=0
- **結論**: **メッシュは変更せず、`run_0098` の本段開始場を共通起点に、2 次のまま実効 CFL 5 対 1 の A/B を行う。**
- **extra**: `case/45.isobutane_m6_d155/README.md`, `procedures/divergence-and-startup.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| H1「薄い壁セル→ω→T 暴走」を主因とする | **要再検証・Major** | `run_0095` では `rms_roe` が先に増加し、`run_0098` の ω 増大も因果順序を示していない。同ケースには乱流更新凍結・層流・Euler でも発散した記録がある（[README:112](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:112)）。今回は旧試験と格子・コード世代が異なるため H1 を除外もできない。まず同一起点で CFL だけを分離する。 |
| A、C を次手として採る | **却下・Major** | メッシュ分布と起動を同時に動かすと帰属できない。また y₁⁺≤1 は出口 δ_E の格子依存が小さい証拠ではない。今回の診断では細分格子を固定する。 |
| B の CFL 低減 | **採用・Major：固定 CFL 1 に限定** | `prep_c2pin.py` は `cfl_main=5.0` を明示指定しており、problem YAML の値だけ変えても効かない（[prep_c2pin.py:11](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/prep_c2pin.py:11)）。生成済み `solverConfig.yaml` の `time.deltaT.cfl_pseudo` を変更し、構造として差分確認する。 |
| 完走後、そのまま②の δ_E 較正へ進む | **却下・Major** | 現スクリプトは①の完走後、収束・準定常のゲートを挟まず `c2pin_solve.py` を呼ぶ（[run_finemesh_pin.sh:19](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/run_finemesh_pin.sh:19)、[同:28](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/run_finemesh_pin.sh:28)）。診断試行は①だけで停止し、判定後に後工程へ進める。 |

結論: **メッシュは変更せず、`run_0098` の本段開始場を共通起点に、2 次のまま実効 CFL 5 対 1 の A/B を行う。**

第 1 仮説: **H2 のうち、本段への CFL 引上げが細分格子上の反復安定限界を超えている。** 確度: **中**。  
「1 次段が短すぎる」まで支持する証拠はない。

- 根拠: ブリーフによれば、`case/45.isobutane_m6_d155/run_0098_ns_finemesh_pass_staged/` は 1 次・CFL 0.5/1 を各 3000 step 完走し、2 次・CFL 5 への変更後 35 step で発散した。ただし実装では次数・CFL に加え、設定文字列が該当すれば `nStepInner` も 10→5 に戻る（[runner_axismach.py:1022](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1022)）。現在の比較は単因子ではない。
- 同ケースの過去の記録も、壁際の ω 爆発を流れ側の反復不安定の症状として説明している。ただし、今回の機構が同じだと確定する根拠にはしない。
- 反証条件: 同一の健全な開始場から **2 次・CFL 1 でも同じ前駆症状を伴って発散する**なら、「5→1 の低減だけで今回の即発散を回避できる」という仮説は棄却する。

第 2 仮説: **細分による格子分布・近壁勾配・SST/粘性項の変化が、低 CFL でも残る不安定を作っている。** 確度: **低〜中、未確認**。ω が大きいだけでは支持できない。

第 3 仮説: **cross-mesh 初期化後の近壁状態の不整合が残っている。** 確度: **低、未確認**。`prepare_ns` は補間後、一定 γ のエネルギー式と Sutherland 粘性で ω の床を再設定する（[runner_axismach.py:954](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:954)）。TP の実温度・輸送物性と一致する保証はなく、単なる「run_0092 の補間場」ではない。ただし、これが今回の発散原因という証拠はない。

判別 A/B:

- **共通起点**: `run_0098` の mid 終了場。実装は各段末場を `nozzle.h5` に保存してから `res_*` を削除する（[runner_axismach.py:1018](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1018)）。したがって AWS 側の `nozzle.h5` が候補になる。後から上書きされていないことと、保存量の有限性を確認して凍結する。
- **変更は実効 CFL だけ**: A=5、B=1。両腕とも `convMethod: 1`、`limiter: 2`、`implicitRelax: 0.7`、本段の同じ `nStepInner`、同じバイナリ・BC・種物性・メッシュを使う。表示用 `cfl` も対応値に揃える。`prepare_ns` の再実行や `stages="ramp"` は使わない。
- **長さ**: 各 1000 step 上限。最初の 100 step は毎 step 保存して発散前の順序を追う。これは即発散の診断であり、生産合格の判定期間ではない。
- **見る量**: 全 `rms_*`、ρ/P/T の極値と床・上限到達数、k/ω/μt の極値。異常節点の座標・壁から何列目かを記録し、壁ノードと第一内部ノードを分ける。
- **事前の解釈**:
  - A が発散を再現し、B が健全に 1000 step 維持 → CFL 依存を支持。「即発散を避けるためにメッシュ修正が必須」は退ける。ただし SST と流れ側の機構分離まではできない。
  - A、B とも発散 → CFL 1 への低減だけで足りるという第 1 仮説を棄却。第 2・第 3 仮説を残す。
  - A が再現しない → 起点・設定・バイナリの同一性を再確認し、原因判定は保留する。

**3 回目の停止・合格条件**は次のように登録する。

1. NaN/Inf、ρ/P/T の非正値、新たな EOS 床・温度上限への到達、または任意の有効残差列が初期 10 step の中央値の 100 倍を 5 step 続けて超えたら、その腕を停止する。これは診断用の打切り条件。
2. B の短期合格は上記停止条件なしで 1000 step 到達し、全有効残差列の末尾 200 step の中央値が直前 200 step 比で +0.1 桁以内。「収束」とは呼ばない。
3. 生産へ進む前には、2 次区間だけの `check_convergence.py`、δ_E を含む報告量の `check_quasisteady.py`、正式な壁解像・品質 VERDICT を揃える。`NOT CONVERGED` はそのまま記録し、1 次段と連結して低下桁数を稼がない。CFL 1 で 12000 step 回しただけでは、従来 CFL 5 の到達状態を保証しない。

やらない方がよいこと:

- ω の最大値や `detectNaN` が最初に表示した変数から原因を決めること。終端の `res_nan_35.h5` だけでは、P/ρ/T と ω の時間的な先後は復元できない。
- `stages="ramp"` を「full の後にランプ」と解釈すること。実装上は別分岐で、最初から 2 次のランプになる（[runner_axismach.py:1033](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1033)）。
- 診断と同時に下流格子、緩和率、乱流モデルなども変更すること。

呼び出し側の前提への異議:

- **Major — Q2: A で確認できる格子感度の範囲を限定すべき。** 下流第一層を旧値へ戻して壁解像 PASS を得ても、出口 δ_E の壁法線解像への独立性は示せない。#11 の「選んだ細分前後の差 ≤1 %」は評価できるが、全体の格子独立性とは区別する。比較は **同じ run_0092 の物理壁・r_t・k_f 上で、①と旧格子の準定常場の未緩和 δ_E** を使う。②③で形状を較正し直した後との差には形状変更も混ざる（[plan:92](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:92)）。
- **Major — 第一セルだけの変更ではない。** `_radial_fracs` は第一セル厚から断面全体の等比率を解き直す（[mesh2d.py:67](/home/sano/work/forge-integ-1005/design/forge_design/meshing/mesh2d.py:67)）。同関数の読み取り実行では、`nj=97`、第一セル比 4.5e−5→1.3e−5 により軸側間隔は **0.07518→0.08876 r_w** に増える。壁側の細分が断面全体の細分にはなっていない。対案は実際の局所間隔分布を比較対象として記録すること。
- **Minor — 旧格子の点数に不整合。** ブリーフは `nj=65`、現在の [problem YAML:62](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/problem_d155_ns_c2pin_final.yaml:62) は `nj=97`。旧 run の実メッシュと `prepare_info.json` で確定する必要がある。

不足情報: AWS 上の実効 config、mid 終了場、段別残差・VERDICT、発散前の場、実行バイナリの SHA-256。対象 run はローカルに存在せず、今回の run 数値はブリーフの報告値として扱った。ファイルサイズと commit だけではバイナリ同一性を確認できない。

ファイル変更・forge 実行は行っていない。**plan 未反映**。呼び出し側で本提案を `plans/active/tooling-nozzle-cfd-pinned-initial-line.md` §5.1 #11・§6 に事前登録してから試行する。
