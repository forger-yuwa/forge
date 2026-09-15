# codex レビュー: thermophysics-cea-mole-fraction-species (plan)

- **plan**: [`plans/active/thermophysics-cea-mole-fraction-species.md`](../../plans/active/thermophysics-cea-mole-fraction-species.md)
- **stage**: `plan`
- **date**: 2026-09-15
- **commit**: `f1e40f41` (feature/sern-design)
- **codex**: effort `high`, 6.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M6/m1
- **focus**: §4.5 の SERN 統一 (full|lumped を run ごとに選択; lumps スキーマ; full での排気率を元素混合分率で定義; SERN runner の SPECIES_ORDER 撤去) の妥当性と見落とし。前回レビュー (M1–M8) は反映済みなので再指摘は不要、§4.5 と §6 の SERN 項に集中
- **extra**: `design/forge_design/evaluate/runner_sern.py`, `design/forge_design/gas/frozen.py`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

## 判定: **GO-with-changes**

`full | lumped` の統一は妥当です。既存の `split_h2o` と SERN の固定 `[EXH,AIR]` を一般化する目的にも重複はありません。  
ただし、**元素混合分率を排気由来の割合として無条件に使う設計は成立しません**。段階 restart にも、種保存量を落とす既存経路があります。  
以下を実装前に計画へ反映してください。前回 M1–M8 の再指摘は除いています。

## 指摘一覧

### 1. **Major — 元素混合分率では、粘性 `full` の流入元と組成を一意に復元できない**

**根拠:** [speciesTransport_d.cu:248](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:248) は種ごとに異なる拡散係数を使い、[263行](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:263) で質量保存補正を加えます。`speciesDiffusionMethod` の既定は [solverConfig.hpp:387](/home/sano/work/forge/solver_density_cuda/input/solverConfig.hpp:387) の `1` です。

元素の保存式は成立しても、その拡散流束は `Σ_s w_{e,s} J_s` です。**frozen は差動拡散を止めません**。したがって、

- 元素ごとに求めた `ξ` が異なり得る。
- `Y = ξ Y_exh + (1−ξ) Y_air` から外れた組成を、単一の `ξ` では復元できない。
- `full` の場をこの式で再構成する warm restart は、差動拡散による組成変化を消す。

これは前回 M6 の「両モードの解が異なる」という話に加え、**今回導入するアクセサの用途自体の問題**です。

**対案:** 流入元の識別には、輸送則を明示した独立トレーサを使ってください。元素混合分率は二流体混合からのずれを調べる診断量に限定します。同一条件 restart は全 `roY{s}` を保持し、作動点変更時の組成再構成は、情報を落とす初期化操作として別契約にしてください。

### 2. **Major — power-off の退化処理が物理的にも float32 にも不適切**

**根拠:** 計画 §4.5 は組成が同じなら「排気なし」とします。しかし、[実際の `m4_off` 入力:44](/home/sano/work/forge/case/46.sern_design/problem_moo_frozen_tp_cycle3op.yaml:44) は燃料遮断後もノズル入口流を持ちます。燃料がないことと、入口由来の流れがないことは別です。

さらに、その入力の空気組成と [既定外気 `AIR_MOLE`:19](/home/sano/work/forge/design/forge_design/gas/frozen.py:19) は微量成分が異なります。現行 DB と換算関数で再計算した結果は以下です。

| 確認量 | 値 |
|---|---:|
| 最大差として選ばれる元素 | N |
| `Z_N,exh − Z_N,air` | `−4.049095×10⁻⁴` |
| 逆算の増幅率 `1/abs(ΔZ)` | `2469.69` |
| 線形混合組成を float32 保存しただけの最大 `ξ` 誤差 | **`7.36×10⁻⁵`** |

最後の値は `ξ=0…1` の10,001点で確認しました。CFD 誤差を加える前から `1e-6` を超えます。

**対案:** `powered` と流入元ラベルを分離し、power-off でも入口トレーサを維持してください。元素診断には完全退化・近退化の判定と誤差増幅の上限を設け、成立しない場合は「定義不能」を返します。`0` への置換はしません。

### 3. **Major — `stream` lump と `keep` を併用したときの質量配分が未定義**

**根拠:** [plan:129](/home/sano/work/forge/plans/active/thermophysics-cea-mole-fraction-species.md:129) は排気全体を `EXH` に畳み、同時に `keep: [H2O]` を許します。一方、[139行](/home/sano/work/forge/plans/active/thermophysics-cea-mole-fraction-species.md:139) は常に `ξ=Y_EXH` としています。

`keep` を lump から除けば、純排気でも `Y_EXH≠1` です。実入力 `m6_on` の換算値は、

- `Y_H2O = 0.2411091186`
- H₂Oを除いた lump の割合は `0.7588908814`

となります。除かなければ H₂O を二重計上します。また、排気と外気が同じ実種を含むこと自体は正常なので、単純な「種集合の重複禁止」でも解決しません。

**対案:** **流れごとの質量配分**を正本にしてください。

- `keep` を各流れから先に取り出す。
- 残りをその流れに対応する lump に配分する。
- 各流れについて未配分・二重配分・空 lump を検査する。
- lump→実種の展開行列と、流れ→輸送種の入口ベクトルを保存する。

排気率を擬似種名 `EXH` に依存させず、指摘1の流入元トレーサで解決します。

### 4. **Major — 通常の段階 restart が `roY{s}` をコピーしていない**

**根拠:** [runner_sern.py:416](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:416) の `restart_by_index()` は7変数だけをコピーし、種保存量を除外しています。warm 適応段の直後にも [531行](/home/sano/work/forge/design/forge_design/evaluate/runner_sern.py:531) から呼ばれます。3D の [warm_from_same_mesh():150](/home/sano/work/forge/design/forge_design/evaluate/runner_sern3d.py:150) にも同じ欠落があります。

実関数をメモリ上の入出力に接続して再現しました。

```text
コピー元: ρ=2, roY=[0.8, 1.2]
コピー先の旧IC: ρ=1, roY=[1, 0]
実行後: ρ=2, roY=[1, 0] → ΣY=0.5
```

後から再正規化しても `[1,0]` となり、元の `[0.4,0.6]` は戻りません。前回 M3 の `interp_field.py` 対策だけでは、この独自経路を修正できません。

**対案:** §5に両関数を明記し、全種・トレーサの保存量とメタデータを引き継いでください。検証には、**段階切替直前・直後の全 `roY{s}`、`ΣρY/ρ`、温度の保持**を追加します。後方互換比較からは、この既存バグの修正による差を分離してください。

### 5. **Major — 元素計算に必要なメタデータを取得・保存する設計がない**

**根拠:** [plan:89](/home/sano/work/forge/plans/active/thermophysics-cea-mole-fraction-species.md:89) の `ResolvedSpeciesDB` は名前・MW・係数・出典だけです。[CEA parser:55](/home/sano/work/forge/solver_density_cuda/tools/cea_thermo_to_species_db.py:55) も元素組成を読み取らず、[出力:90](/home/sano/work/forge/solver_density_cuda/tools/cea_thermo_to_species_db.py:90) に残しません。

これでは `w_{e,s}` の正本がありません。任意の CEA 種名や擬似種名から原子組成を推測する実装は不適切です。コメントだけでは通常の YAML 読込で情報が消えます。

**対案:** 実種の原子組成、lump の構成比、各流れの正規化済み組成、輸送種順序を**機械可読メタデータ**として保存してください。CEA の元素欄から取得し、別名解決後も維持します。後処理・restart は元の問題 YAML を再解釈せず、run に保存した情報を使います。

### 6. **Major — SERN の検証が主要な変更経路を通らず、比較精度と定常性閾値も釣り合わない**

**根拠:** [plan:191](/home/sano/work/forge/plans/active/thermophysics-cea-mole-fraction-species.md:191) は1作動点の node Euler 比較です。これでは power-off、作動点変更、`keep` 併用、粘性の差動拡散、3D restart を検証できません。

また、既存 SERN ゲートには次の違いがあります。

- [sern_gates.py:124](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:124): `require_residual_pass=False` が既定。
- [sern_forces.py:30](/home/sano/work/forge/design/forge_design/metrics/sern_forces.py:30): 定常性閾値は drift **2%**、fluctuation **5%**。計画の力比較 **0.1%** より大幅に緩い。
- [sern_gates.py:30](/home/sano/work/forge/design/forge_design/metrics/sern_gates.py:30): 場の健全性検査に `roY{s}` がない。

実際に `case/46.sern_design/run_0032_cycle_euler_m4off/` を再確認すると、

```text
check_convergence: NOT CONVERGED (stalled/plateau)
check_quasisteady: C_T / C_L / C_M = STEADY
```

でした。これは旧 CPG run の判定例であり、`full/lumped` の検証結果ではありません。README 記載の frozen 基準 `run_0094–0099` はこの checkout になく、実測を再確認できませんでした。

**対案:** §6に以下を追加してください。

- powered／power-off、微小組成差、`keep`、種順序変更の単体試験。
- 同一条件 restart と `m6_on→m10_on` の組成変更試験。
- 小型粘性二流体ケースと3D経路の確認。
- 比較する両 run の全種残差を含む `PASS`、全種の有限性・非負性。
- ノズル力・機体力・出口運動量それぞれの時系列を保存し、**比較許容差より十分小さい変動幅**で `STEADY` を要求。
- ゼロ近傍の力には `F_ideal` 基準の絶対許容差を定義。

### 7. **Minor — §4.5追加後も旧スキーマと完了条件が残っている**

**根拠:** [plan:84](/home/sano/work/forge/plans/active/thermophysics-cea-mole-fraction-species.md:84) は文字列 `tp_species`＋`tp_lump`、[125行](/home/sano/work/forge/plans/active/thermophysics-cea-mole-fraction-species.md:125) は mapping です。[212行](/home/sano/work/forge/plans/active/thermophysics-cea-mole-fraction-species.md:212) の完了条件も「検証1–5」のままで、SERN が明示されていません。

**対案:** 内部表現を mapping に一本化し、旧形式は入力時の別名変換として列挙してください。競合指定は拒否し、SERN の検証を完了条件へ追加します。

## 推奨

**統一モードを採用し、組成表現と流入元ラベルを分離する設計に修正してから実装する**ことを推奨します。

実装前の優先順は、①アクセサの用途・退化処理、②lump配分とメタデータ、③全 restart 経路、④SERN検証ゲート、⑤スキーマ・完了条件の同期です。独立トレーサには既存の [汎用スカラ輸送コア](/home/sano/work/forge/solver_density_cuda/cuda_forge/scalarTransport_d.cuh:11) を利用できますが、登録・BC・陰解法更新・残差・restart まで計画に含める必要があります。

**plan 未反映**：依頼どおりファイルは変更していません。

指摘数: Critical 0 / Major 6 / Minor 1
