# codex 諮問 (diagnose): transport-stage2-gpu

- **brief**: [`notes/reviews/briefs/2026-09-27-transport-stage2-gpu.md`](../../notes/reviews/briefs/2026-09-27-transport-stage2-gpu.md)
- **plan**: [`plans/active/thermophysics-solver-owned-species-db.md`](../../plans/active/thermophysics-solver-owned-species-db.md)
- **date**: 2026-09-27
- **commit**: `4eb97d15` (feature/species-transport)
- **codex**: effort `xhigh`, 4.3 min, rc=0
- **結論**: **展開をモル基底に統一した共通評価経路を作り、独立参照付きの full/lump・CFD 0 step 試験を最初に通す。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（修正の最優先は、展開行列の基底の明示）：

| 項目 | 採否 | 根拠・対案 |
|---|---|---|
| 素案の「Y → 展開行列 → 実種Y」 | **要再検証・Major** | 既存行列はモル基底。[speciesTransportDB.hpp:53](/home/sano/work/forge-species/solver_density_cuda/input/speciesTransportDB.hpp:53)。**輸送種Y → 輸送種X → 既存行列で実種X**に統一する。 |
| 新指定がある場合だけ切替、セル・壁で評価関数を共有 | **採用** | 旧経路を維持するユーザ決定と整合。[plan:152](/home/sano/work/forge-species/plans/active/thermophysics-solver-owned-species-db.md:152)。起動拒否の解除は接続試験後。 |
| host の同じ `transport_mix` との一致だけで合格 | **却下・Major〔検証条件〕** | CPU/GPU の共通誤りを検出できない。plan §5.1 #5t2 の独立参照条件を維持し、組成変換から独立評価する。 |
| ηᵢⱼ を各組1回評価 | **採用・Minor〔実装との差〕** | 現コードは i→j と j→i で再評価する。[transportMix_d.cuh:227](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/transportMix_d.cuh:227)。各組を1回評価し、両方向の分母へ加算する。 |

結論: **展開をモル基底に統一した共通評価経路を作り、独立参照付きの full/lump・CFD 0 step 試験を最初に通す。**

第 1 仮説: **既存の `expand` を質量分率に直接掛けると、lump のGPU物性が別組成の値になる。** 確度: **高〔素案を文字どおり実装した場合〕**

  根拠: 行列は lump 内モル分率であり、host も `Xreal[r] += Xs[s] * expand[s,r]` と使っている。[speciesTransportDB.cpp:494](/home/sano/work/forge-species/solver_density_cuda/input/speciesTransportDB.cpp:494)。現在のセル経路も Y→X を混合評価前に行う。[gasProperties_d.cu:77](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/gasProperties_d.cu:77)。

  独立した代数検算では、`L={N2:0.5, He:0.5}`、輸送種モル分率 `X_L=0.8, X_He=0.2` の正解は実種 `X_He=0.6`。既存行列をYへ直接掛けてから実種Xへ変換すると **0.88730848** になる。これはCFD実測ではなく、MWを N2=0.0280134、He=0.0040026 kg/mol とした検算。

  対案: 正規化済み輸送種Yから `X_s=(Y_s/M_s)/Σ(Y/M)` を作り、`Xreal_r=Σ_s X_s·expand[s,r]` とする。重複実種は加算し、そのまま `transport_mix` に渡す。セル・壁ともこの処理を共有する。

  反証条件: 実装が既にモル基底で展開している、または別途正しく質量基底へ変換した行列を使い、下記試験で展開後組成と物性が独立参照に一致すること。

第 2・第 3 仮説: 追加しない。GPU接続前なので、発散原因を推定する段階ではない。

判別 A/B: **変えるのは組成の表現だけ**。T=400 K、ρ=1、両実種 `cea`、凝縮OFF、時間更新なしで、実際のセル経路・ghost・壁評価を各1回呼ぶ。

- **A:** full の `X_N2=0.4, X_He=0.6`。
- **B:** `L={N2:0.5, He:0.5}` と独立 `He`、`X_L=0.8, X_He=0.2`。Heの重複合算も同時に検査する。
- **A合格・Bだけ不合格なら**、展開・重複合算経路を原因候補として残す。Bの `X_He≈0.88730848` は基底取り違えの具体的な署名。
- **A/Bとも合格なら**、この入力で第1仮説を棄却する。両方不合格なら lump 固有原因へ帰属せず、共通の係数転送・評価経路を調べる。

段2の合格条件は、この試験を含めて次に固定する。

| 対象 | 事前に固定する条件 |
|---|---|
| double評価 | **実際にGPUへ渡した入力**を独立参照へ渡し、μ・λの最大相対誤差 ≤1e−12。展開後Xの最大絶対誤差 ≤1e−12。 |
| float格納 | 独立double参照比 ≤1e−5。さらに参照をfloatへ丸めた値との差 ≤2 ULPを確認する。 |
| full/lump・列挙順 | double入力では ≤1e−12。別々にfloat化した `roY` からの物性比較は ≤1e−5とし、各入力に対する独立参照との照合を必須にする。 |
| セル・ghost・壁 | 同一の実入力で上記基準を全経路が満たす。未更新要素・NaN/Inf・非正のμ/λは0件。ghost末尾まで確認する。 |
| 試験範囲 | 単成分、重複lump、実種数12と上限32、ゼロ分率を含める。T=200/253.15/400/500/600/700/1000/2000 Kと各フィット境界を確認。GPUの境界両側は隣接float値を使う。実種上限は既に32。[transportMix_d.cuh:37](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/transportMix_d.cuh:37) |
| 旧経路 | 輸送指定なしでは、**同一の固定入力**に対する `vis_lam`・`thermCond` のビット差0。時間発展後の異なる入力同士を比較しない。 |

やらない方がよいこと: **case/44 のEuler起動成功や凝縮200 stepを、GPU物性接続の合格根拠にしない。** 素案自身が前者ではμ・λ不使用、後者では合否判定しないと定義している。NS・凝縮の妥当性は段3に残し、今回 `cfg.visc` の意味や気相組成処理まで同時変更しない。

呼び出し側の前提への異議:

- **Major〔比較精度〕:** 「同じ物理組成」でも、full/lumpで別々にfloat化した `roY` は同一入力ではない。上の検算例では、正しい変換でも `X_He` は full **0.60000000561**、lump **0.59999999826**、差 **7.35e−9**。したがって実配列経由の比較へ一律1e−12を要求するのは誤り。上表のように入力丸めと評価誤差を分ける。
- **Minor〔性能見積り〕:** 現在の `transport_mix` は全分率正なら12種でη評価**132回**。66組すべてがCEAフィットとも限らない。[transportMix_d.cuh:199](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/transportMix_d.cuh:199)。性能予算としては、対象GPU・入力を固定し、旧step時間に対する物性更新時間の増分を**10%以内**とする案を推す。これは提案値であり、現時点の実測に裏付けられた上限ではない。

  μ・λの利用箇所については、次を接続確認表へ加える。**以下は確認済みの依存関係であり、すべてを今回改修せよという指摘ではない。**

| 経路 | コードから確認した内容・扱い |
|---|---|
| 粘性流束・熱流束 | [viscousFlux_d.cu:164](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/viscousFlux_d.cu:164)、[同:265](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/viscousFlux_d.cu:265) は局所配列を読む。 |
| SST・遷移・スカラー拡散 | [ransSource_d.cu:97](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/ransSource_d.cu:97)、[transition_d.cu:149](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/transition_d.cu:149)、[scalarTransport_d.cu:120](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/scalarTransport_d.cu:120) に新μが伝わる。 |
| 化学種・受動種の拡散 | LJの混合平均係数は据え置けるが、**定数Sc経路は新μにより変わる**。[speciesTransport_d.cu:259](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/speciesTransport_d.cu:259)、[passiveKernels_d.cuh:171](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/passiveKernels_d.cuh:171)。 |
| 軸対称・共役壁熱伝達 | [axisymmetricSource_d.cu:239](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/axisymmetricSource_d.cu:239)、[conjugateWall.cpp:242](/home/sano/work/forge-species/solver_density_cuda/conjugateWall.cpp:242) は局所μ・λを読む。 |
| SST壁関数 | 局所μ・λによるPrは更新されるが、回復係数は **`cfg.prandtlLam` のまま**。[ransWallFunction_d.cu:347](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/ransWallFunction_d.cu:347)、[同:583](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/ransWallFunction_d.cu:583)。自動的に全閉包が新物性へ統一されるわけではない。 |
| CFL・陰解法 | [setDT_d.cu:82](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/setDT_d.cu:82)、[timeIntegration_d.cu:813](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/timeIntegration_d.cu:813)、[同:1496](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1496) は定数粘性を使う。FP64 accumulator側にも渡される。既存の近似として明記し、段2の物性一致から安定性まで保証しない。 |

不足情報: 確認時のHEADは **`4eb97d15`** で、ブリーフの `89f0a5c7` と異なる。試験対象commit・バイナリハッシュを固定する必要がある。段1の実測結果、GPU性能基準、境界更新後の物性再評価順序は今回の確認では保証していない。runデータは未読、forge未起動、ファイル変更なし。**plan未反映**；呼び出し側が §4.3c・§5.1 #5t2・§6へ反映する。
