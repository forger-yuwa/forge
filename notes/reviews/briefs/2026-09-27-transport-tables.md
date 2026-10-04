# 諮問ブリーフ: 輸送物性の表引き化 (#5t2-3) — 表の刻み・範囲・精度基準の事前固定 (2026-09-27)

作業ツリー `/home/sano/work/forge-species` (ブランチ `feature/species-transport`, HEAD cf10ecf4 以降)。関連 plan: `plans/active/thermophysics-solver-owned-species-db.md` §5.1 #5t2-2・#5t2-2r・§10。
AGENTS.md エスカレーション 6 (数値カーネルの変更前)。

## 問い (1 つ)

段 2 で GPU に接続した種ごとの輸送物性 (double で CEA フィット・IAPWS・Chapman–Enskog を毎 step 評価) が遅いので、
**種別 μᵢ・λᵢ と組ごとの ηᵢⱼ を ln T の表 (区分 3 次) にして float で引く**ことにした (ユーザ決定)。
下の素案の刻み・範囲・補間・精度基準・性能基準で進めてよいか。**穴があれば 1 つに絞って直し方を示し、事前に固定すべき合格条件を示してほしい。**

## 観測事実

- 段 2 の実装 (`cuda_forge/transportMix_d.cuh`: `TransportTemp`, `transport_species`, `transport_eta_pair`, `transport_mix`, `transport_mix_Y`; 組成はモル基底で展開) は正しさの合格条件をすべて満たした
  (`tests/unit/test_transport_gpu.py`: double 3.5e-15・float 格納 5.5e-8 かつ 0 ULP、既定経路ビット一致)。
- **性能** (RTX 3060, FP64 は FP32 の 1/64, 23725 CV, case/44 メッシュ, 300 step × 2): gas_properties が旧 0.062 ms/step → 実種 5: 0.409 (+14 % of step)・12: 2.567 (+94 %)・32: 19.37 (+800 %)。
  提案値は「物性更新の増分が旧 step 時間の 10 % 以内」。ln T と 1/T の使い回しは入れた後の値。
- 前例: 凝縮経路の float 化 (`cuda_forge/condensationTables_d.cuh`, plan `plans/accepted/condensation-float-speedup.md`): double で作った区分 3 次表を float で引く。
  N2 20–125.6 K h=0.1 K、H2O 120.15–1200.15 K h=0.25 K、**物性の接続点を格子点に置く**。式を直接 float 化すると H2O 潜熱が 9e-3 ずれたので表方式にした。範囲外は旧 double 関数へ委譲。
- 輸送の式の性質:
  - CEA フィット `ln f = A lnT + B/T + C/T² + D` は区間ごと (多くは 200–1000–5000(or 6000)–15000 K、H2O は 373.2–1073.2–5000–15000 K)。**区間の境界で値・傾きが連続とは限らない** (CEA 本体も区間を切り替えるだけ)。
  - H2O custom: 253.15 K 未満は IAPWS 勾配の冪 (C¹)、253.15–500 K は IAPWS、500–700 K は ln 空間 smoothstep (C¹)、700 K 以上は CEA (区間境界 1073.2・5000 K あり)。
  - kinetic: Chapman–Enskog (Neufeld Ω22、T* を [0.3, 100] にクランプ; クランプ点で傾き不連続)、Brokaw 補正、λ は修正 Eucken で cp (NASA-9, Tmid=1000 K で区間切替) に依存。
  - ηᵢⱼ: 両 kinetic なら二元 CE (T* クランプ)、CEA 相互作用 (区間あり、例 H2O–N2 300–1000–5000 K)、剛体球近似 (μᵢ・μⱼ に依存 → 種別表から作れる)。
- 混合 (φᵢⱼ・ψᵢⱼ の和と割り算) は O(n²) の四則演算。現状 double。

## 素案 (当方)

1. 表: 各実種の ln μᵢ・ln λᵢ と、各組 (CE または CEA 相互作用のもの) の ln ηᵢⱼ を、**ln T の等間隔格子**で double から作り、各節点で値と d/d(ln T) を持つ **3 次 Hermite** (区間境界は連続でないので、**すべての式の区間境界・接続点・T* クランプ点・NASA Tmid を節点に置き、節点で左右の片側値を別に持てるようにする** か、境界ごとに表を分割)。剛体球の ηᵢⱼ は種別表の μ から実行時に計算。
2. 範囲: 150–15000 K (表外は端の式で double 評価、または端の値で延長 — 凝縮と同じく「範囲外は double 関数へ委譲」を推す)。刻み: Δ(ln T) = 1/256 程度 (≈4.6 k 点 × 実種 32 × 2 物性 + 組 496 × … でメモリ数 MB)。
3. 評価: float で表引き (区間探索は等間隔なので直接 index)、混合は float (FMA)。μ・λ は float 格納。
4. 合格条件 (案): (i) 表引き単体が独立 double 参照と **相対 ≤2e-6** (全範囲・全節点の中点・各境界の両側)、(ii) GPU の格納値が独立 double 参照と **≤1e-5** (段 2 と同じ基準; ULP 条件は外す)、
   (iii) 既定経路ビット一致 (段 2 と同じ)、(iv) 性能: 物性更新の増分が旧 step 時間の 10 % 以内 (実種 5・12 で; 32 は記録のみ)、(v) 段 2 の A/B・lump/full・列挙順・セル/ghost/壁の試験をそのまま再実行。

## 当方の懸念 (棄却してよい)

- 区間境界で不連続な式を 1 本の Hermite 表にすると境界付近で誤差が大きくなる。境界を節点に置く方法と、表を区間ごとに分ける方法のどちらが簡単で安全か。
- float の混合で φ・ψ の和が桁落ちしないか (X が 1e-6 級の微量種)。
- 性能 10 % が n=12 で達成できるか (O(n²) の float 四則 + 表引き 66 組)。達成できないときに何を諦めるか。
- 表引きに変えたことで「全 CEA 指定なら FCEA2 と ≤0.1 %」は自動で保たれるはず (表誤差 ≤2e-6) だが、試験に残すべきか。

## 禁止事項 (厳守)

- ファイルを変更しない。**`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md`, `trans.inp`・`thermo.inp` 全体を読まない**。下の範囲以外のファイル読みをしない (grep は可)。
- 推奨は 1 つに絞る。根拠は `ファイル:行` か本ブリーフで示す。

## 読んでよいもの (パスは `/home/sano/work/forge-species/` 基準)

- `sed -n '1,360p' solver_density_cuda/cuda_forge/transportMix_d.cuh`
- `sed -n '1,120p' solver_density_cuda/cuda_forge/condensationTables_d.cuh`
- `sed -n '1,80p' plans/accepted/condensation-float-speedup.md`
- `sed -n '60,260p' plans/active/thermophysics-solver-owned-species-db.md`
- `sed -n '1,80p' solver_density_cuda/tests/unit/test_transport_gpu.py`
