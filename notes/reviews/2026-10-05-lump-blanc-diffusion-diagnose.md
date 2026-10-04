# codex 諮問 (diagnose): lump-blanc-diffusion

- **brief**: [`notes/reviews/briefs/2026-10-05-lump-blanc-diffusion.md`](../../notes/reviews/briefs/2026-10-05-lump-blanc-diffusion.md)
- **plan**: [`plans/active/thermophysics-solver-owned-species-db.md`](../../plans/active/thermophysics-solver-owned-species-db.md)
- **date**: 2026-10-05
- **commit**: `2d5dac2c` (feature/species-transport)
- **codex**: effort `high`, 6.1 min, rc=0
- **結論**: **GPU 実装に入る前に、質量加重縮約と従来の対称 Blanc 縮約だけを切り替える 0-step 流束 A/B を行い、質量流束の一致とエネルギー・重複の近似誤差を分離して §4.4・V4 を確定してください。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

**採否表**

主案は、**非重複 lump の質量流束を縮約する方法として採用できます。ただし、重なり・エネルギー流束・検証基準が未確定なので、そのまま実装へ進むのは不可**です。

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 実種へのモル分率展開と `D_L = Σ y_{r｜L} D_r` | **採用** | `J = −ρD∇Y` に対する集約なので質量重みが正しい。`y_{r｜L} = E[L,r] M_r/M_L` と明記する。現行の勾配・補正は [speciesTransport_d.cu:323](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:323)。 |
| 「補正後の外部種流束は `full` 比 +0.50 %」の継承 | **却下／Major M1** | この数値は旧縮約の検算で、新案の期待値ではない。非重複・固定内部組成・同じ面状態なら、新案は補正後の外部種流束まで一致する。対案は V4d を一致試験に変更し、エネルギーの縮約誤差を別判定にする。[plan:168](/home/sano/work/forge-species/plans/active/thermophysics-solver-owned-species-db.md:168) |
| lump 同士・構成種の重なり | **要再検証／Major M2** | 実種組成を合算しても、輸送ラベルへの流束の配分は決まらない。下記のゼロ組成勾配反例がある。「自然に決まる」は撤回し、重複時の縮約規約と許容誤差を明記する。既存 §4.4 の二重和とは別のモデルなので、置き換えを設計判断として記録する。[plan:171](/home/sano/work/forge-species/plans/active/thermophysics-solver-owned-species-db.md:171) |
| lump なしは現行関数を呼ぶ | **採用。ただしビット一致の保証範囲を限定** | 固定入力の係数評価には妥当。一方、残差には `atomicAdd` があるため、CFD 全場のビット一致は分岐だけでは保証できない。係数のビット試験と、反復ノイズを含む場の回帰を分ける。[speciesTransport_d.cu:345](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:345) |
| 性能対策 | **表引きを採用** | 面状態で評価する現在の離散化を保持できる。節点で評価して平均すると別の離散化になる。式の確定後、`D_rq P` を表引きし、各非対角対を面ごとに一度だけ評価する。既存表の境界処理を踏襲する。[transportTables_d.cuh:11](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/transportTables_d.cuh:11) |
| V4a〜d | **現状では不足／Major M3** | エネルギー・重複・純成分極限・lump を含む二相経路・保存収支が不足している。既存 V4 自体も収支式と許容差の事前登録を要求している。[plan:339](/home/sano/work/forge-species/plans/active/thermophysics-solver-owned-species-db.md:339) |
| 旧平均 LJ への切替キー | **新設しない。ただし記録だけで済ませない／Major M4** | 平均 LJ 廃止は既決方針。旧結果は固定した旧バイナリ・入力で再現する。新しい拡散モデルの版・縮約規約を解決済み記録と互換性照合へ反映する。現在のハッシュ対象には lump の合成規約が入る。[speciesDB.cpp:1245](/home/sano/work/forge-species/solver_density_cuda/input/speciesDB.cpp:1245) |

結論: **GPU 実装に入る前に、質量加重縮約と従来の対称 Blanc 縮約だけを切り替える 0-step 流束 A/B を行い、質量流束の一致とエネルギー・重複の近似誤差を分離して §4.4・V4 を確定してください。**

第 1 仮説: **主な設計上の問題は質量重みではなく、異なる縮約モデルの期待値を混ぜていることです。** 確度: 高  
  根拠: 非重複 lump では、内部質量分率を a_r とすると、

  `Y_r = a_r Y_L`  
  `Σ_{r∈L} J_r = −ρ(Σ a_r D_r)∇Y_L = J_L`

  したがって全種の補正前流束総和も同じになり、

  `Σ_{r∈L} J_r* = J_L*`

  が成立します。外部独立種の補正後流束も同じです。これは**同じ状態・勾配における流束評価の恒等式**であり、時間発展した `full` と lump の場が一致するという意味ではありません。

  一方、現行エネルギー結合は輸送種の `h_s J_s*` です。[speciesTransport_d.cu:340](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:340)  
  同じ datum で `h_L = Σ a_r h_r` としても、一般には

  `q_full − q_lump = −ρ Σ_L [Σ_{r∈L} a_r h_r D_r − h_L D_L] ∇Y_L`

  が残ります。**内部組成が固定でも、この共分散項はゼロになりません。**

  Python の独立計算で確認した合成例は次のとおりです。実在気体や CFD run の測定値ではありません。

  - 分子量比 `(2, 4, 8)`、実種モル分率 `(0.4, 0.4, 0.2)`。
  - `D₁₂/D₀ = 1`、`D₁₃/D₀ = 2`、`D₂₃/D₀ = 4`。
  - lump は種 1・2 を等モルで含むため、内部質量重みは `(1/3, 2/3)`。
  - `Y_L = 0.6`、`Y₃ = 0.4`、`∇Y_L = 1 m⁻¹`。
  - エンタルピーは共通単位 h₀ で `(1, 5, 0)`。

  | 評価 | 外部種 `J₃*` の規格化値 | エンタルピー流束の規格化値 |
  |---|---:|---:|
  | 実種 `full` | 2.115555556 | −7.875555556 |
  | 提案の質量加重縮約 | 2.115555556 | −7.757037037 |
  | 対称 Blanc 係数を両輸送種に使用 | 2.666666667 | −9.777777778 |

  対案: 質量流束は上の恒等式を検証し、エネルギーは `Σ h_s J_s*` という縮約モデル自身の独立参照と照合する。`full` との差は別に定量記録する。エネルギーだけを無条件に `full` の値へ置き換えることは勧めません。

  反証条件: 非重複・固定内部組成・同一面状態・同一勾配・同一補正で、質量加重縮約の外部種流束が独立 double 参照と規格化誤差 `1e−12` を超えて異なること。その場合は基底変換か流束実装を再点検します。

第 2 仮説: **重複構成では、実種の合算後にスカラー `D_s` を各ラベルへ戻す操作が、実種組成に存在しない拡散を作ります。** 確度: 高  
  上と同じ合成例を、輸送種 `[L₁₂, 種1, 種2, 種3]`、質量分率 `(0.3, 0.1, 0.2, 0.4)` で表します。輸送種勾配を `(1, −1/3, −2/3, 0) m⁻¹` とすると、**全実種の組成勾配はゼロ**です。しかし提案の縮約流束を実種へ戻すと、規格化値は `(-0.02962963, +0.02962963, 0)` になりました。総質量保存は満たすため、`ΣJ*=0` の試験だけでは検出できません。

  対案: 重複時は「固定組成ラベルの近似モデル」と明記し、実種へ戻した流束誤差を検証対象にする。SERN の重複組成で許容差を事前登録するまでは、既定投入を合格にしない。また、`n_real=1`・輸送ラベル複数の場合の自己拡散を明示すること。現行関数をそのまま使うと `n==1` でゼロを返します。[thermo_d.cuh:685](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/thermo_d.cuh:685)

第 3 仮説: 無し。

判別 A/B: **変えるのは lump の係数規約だけ**です。上記の合成入力を固定し、時間積分なし・面流束を各一回評価します。

- **A**: 外部種との Blanc 係数を lump 側にも対称に使う。
- **B**: 実種展開後の `D_L = Σ a_r D_r` を使う。
- 見る量: 外部種と lump 合計の補正後流束、`ΣJ*`、`ΣhJ*`。独立 `full` 参照を別計算する。
- **A のみ質量流束が一致するなら**、第 1 仮説の縮約導出・対応付けを棄却する。
- **B のみ一致するなら**、旧 `+0.50 %` を新案の期待値とする前提を棄却する。
- 両方不一致なら実装へ進まず、基底・面状態・補正の対応を調べる。

本レビューの独立代数計算では **B が一致し、A は不一致**でした。次はこの判別を実装予定の評価経路で固定する段階です。

やらない方がよいこと:

- `D_i` だけの一致や、組成勾配のない CFD で V4 を合格にすること。
- 速度対策として節点評価へ変更し、モデル変更と離散化変更を同時に入れること。
- lump なしの回帰だけで二相経路まで保証したとすること。係数の利用先は通常拡散・二相面入力・監査に分かれています。[speciesTransport_d.cu:2228](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2228)、[同:2466](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:2466)
- 旧方式の再現用として残っている `ljSource: [legacy_v1]` を、本件の廃止対象へ広げること。

呼び出し側の前提への異議:

1. **「Blanc と同じ」の成立範囲が広すぎます。** 非重複の外部種から見た係数については成立します。一般の lump 同士については、提案は対称二元係数の二重和ではなく、組成依存の縮約です。二つの lump だけなら、補正後は `J_A* = −ρ(Y_B D_A + Y_A D_B)∇Y_A` と書けます。この式を §4.4 の新しい契約として明記すべきです。
2. **現行 `full` は物理的な厳密解ではありません。** 現行は `(1−X_i)/Σ X_j/D_ij` を `∇Y_i` に掛けています。標準的なモル勾配用・質量勾配用の混合平均係数は区別されます。[Cantera 公式仕様](https://www.cantera.org/3.2/cxx/d8/d58/classCantera_1_1GasTransport.html) 本件では既存式を保存し、検証を「forge の現行 `full` 演算への一致」と限定してください。
3. **性能表の合格基準も先に固定する必要があります。** 推奨は、同じ float 入力の独立 double 参照に対し、二元係数の表引き誤差 `≤2e−6`、縮約後係数 `≤1e−5`。`T*=0.3/100` の両側・微量種・純成分・上限 32 種・表範囲外を含め、範囲外は元の式へ戻す。これは事前登録する提案値であり、達成確認済みではありません。
4. plan #7 の「`viscMethod ≠ 2` の μ・λ も平均 LJ」という記述は古いです。`viscMethod 1` は組成によらない Sutherland です。[solver-settings.md:233](/home/sano/work/forge-species/procedures/solver-settings.md:233)

不足情報: SERN の**混合域を含む**代表組成・勾配と流束誤差の許容値、lump を含む二相試験、V4 の種・エネルギー収支式と正規化・許容差が未提示です。入口の係数変化表だけでは、この不足を埋められません。

ファイル変更・`forge` 起動は行っていません。**plan 未反映**です。呼び出し側で `plans/active/thermophysics-solver-owned-species-db.md` の §4.4・§5.1 #7・§6 V4 に反映してください。
