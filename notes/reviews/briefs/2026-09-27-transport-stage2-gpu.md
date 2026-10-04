# 諮問ブリーフ: 種ごとの輸送物性を GPU に接続する (#5t2 段 2) — 実装前の設計確認 (2026-09-27)

作業ツリー `/home/sano/work/forge-species` (ブランチ `feature/species-transport`, HEAD 89f0a5c7)。関連 plan: `plans/active/thermophysics-solver-owned-species-db.md` §4.3c・§5.1 #5t2 (段 1 完了行 `5t2-1`)。
AGENTS.md エスカレーション 6 (`cuda_forge/` の数値の振る舞いを変える編集の前)。

## 問い (1 つ)

段 1 で host に実装した種ごとの輸送物性 (μᵢ・λᵢ・ηᵢⱼ と CEA frozen 混合則) を GPU の物性計算と壁モデルに接続する、下の素案で進めてよいか。
**穴があれば 1 つに絞って直し方を示し、段 2 の合格条件 (事前に固定する数値) と、見落としている μ・λ の利用箇所があれば挙げてほしい。**

## 段 1 までの状態 (commit 済み)

- 式: `solver_density_cuda/cuda_forge/transportMix_d.cuh` (244 行, `THERMO_HD`, double)。`SpeciesTransportD` (種別モデル・フィット・LJ・双極子)、`TransportPairD` (ηᵢⱼ の出所: CE / CEA 相互作用 / 剛体球)、
  `transport_species` (cea / fit / kinetic (Brokaw 補正, λ 修正 Eucken) / custom:h2o_iapws_cea_v1 (≤500 K IAPWS、≥700 K CEA、間は ln 空間 smoothstep、253.15 K 未満は IAPWS 勾配の冪))、
  `transport_eta_pair`、`transport_mix` (CEA φ・ψ)。**どのカーネルからも呼んでいない。**
- 解決: `input/speciesTransportDB.{hpp,cpp}`。`physProp.transport` があるときだけ解決し、lump を実種へ展開 (同じ実種は合算) した展開行列、実種表、組の表を作る。記録 schema v2 と互換性ハッシュ。
- `main.cpp:1189-1196`: `physProp.transport` がある config は計算起動を拒否 (段 2 で外す)。
- 試験: `tests/unit/test_species_transport.py` ALL PASS (全 CEA 16 状態 vs FCEA2 μ 0.0012 %・λ 0.013 %、選択モデル vs 独立参照 `transport_reference.py` 最大 4.2e-15、H2O 接続点 500/700/253.15 K で値・勾配連続、負例 18 件、既存 config の記録バイト不変)。
- 現行の GPU 経路: `cuda_forge/gasProperties_d.cu:68-85` (`viscMethod == 2`: `roY` → Y 正規化 → X → `thermo_mu_mix`/`thermo_lambda_mix` (内部 float) → `vis_lam`・`thermCond` 配列に float で格納; nCells_all = ghost 込み)、
  `cuda_forge/wmlesWallModel_d.cu:53-54` (壁で同じ関数)。拡散係数 `speciesTransport_d.cu:259` (`thermo_Dmix_species_f`, LJ) は今回変えない。

## 決定済み事項 (plan §4.3c)

- 混合則は CEA 形 (両モード共通)。ηᵢⱼ: 両種 kinetic → 二元 CE、他は CEA 相互作用、無ければ剛体球。
- `viscMethod: 2` は案 C: 実装は 2 を新方式に置き換える形だが、**種別指定なしを起動時エラーにする既定の切り替えは後** (このブランチでは可、共有ブランチへは並行セッションの区切りでマージ)。
- 段 3 (NS 統合) は凝縮 plan の気相組成処理を固定してから。段 2 では組成の作り方 (`roY` そのまま) を変えない。
- 精度: double で評価し float で格納 (現行は内部 float → 精度変更として扱う)。

## 素案 (当方)

1. `physProp.transport` があり `viscMethod == 2` のとき、`gasProperties_d` と `wmlesWallModel_d` は新経路: 輸送種の Y → 展開行列で実種の質量分率 → 実種 X → `transport_mix` (double) → float 格納。
   `physProp.transport` が無いときは現行 (`thermo_mu_mix` 等) のまま (既定の切り替えは別コミット)。
2. デバイスへ: 実種表 (`SpeciesTransportD[nReal]`)、組の表 (`TransportPairD[nReal(nReal−1)/2]`)、展開行列 (`nTransported × nReal`)。nReal の上限を定数で持ち起動時に検査 (現行 `THERMO_MAX_SPECIES 16`)。
   ηᵢⱼ は各非対角組を 1 回ずつ評価、φᵢᵢ = ψᵢᵢ = 1。セルと壁が同じデバイス関数を使う。
3. `main.cpp` の起動拒否を外す。
4. 試験: (i) GPU のセル値 (ghost 含む) と壁値が、同じ T・組成で host の `transport_mix` と double 評価 ≤1e-12 相対・float 格納値 ≤1e-5 (float の丸め幅) で一致 — 1 step 実行 or 物性だけを計算する試験ハーネスで、
   (ii) 重複実種を含む lump と full 展開が一致、(iii) 種の列挙順を変えても一致 (≤1e-12)、(iv) `physProp.transport` の無い config で `vis_lam`/`thermCond` が旧バイナリとビット一致 (既定経路不変)、
   (v) case/44 lump 記法 config に `physProp.transport` を足して短い run (Euler なので μ・λ は使われない; 起動と NaN なしの確認のみ)、(vi) N2+H2O の粘性 run を 1 本 (例 case/16 の node SST 凝縮を収束済み場から 200 step) で起動・NaN なし・μ の変化量を記録 (合否にしない)。

## 当方の懸念 (棄却してよい)

- μ・λ の利用箇所の見落とし: 陰解法の粘性対角 (`physProp.visc` を読む箇所がある; メモ「visc is the stiffness estimate」)、dt の粘性制限、SST の `vis_lam` 利用、境界 ghost の物性、FP64 accumulator 経路、`thermoFloat` の温度反転との順序。
- n² のコスト: nReal = 12 (SERN の EXH 構成種) で 1 セルあたり 66 組 × CEA フィット評価 (exp・log) を毎 step。性能劣化の上限をどう決めるか。
- 凝縮 run で `roY` の H2O が総水分 (液を含む) である問題は段 3 で扱う (段 2 では現行と同じ組成の作り方)。

## 禁止事項 (厳守)

- ファイルを変更しない。**`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md`, `trans.inp`・`thermo.inp` 全体を読まない**。下の範囲以外のファイル読みをしない (grep は可)。
- 推奨は 1 つに絞る。根拠は `ファイル:行` か本ブリーフで示す。

## 読んでよいもの (パスは `/home/sano/work/forge-species/` 基準)

- `sed -n '1,244p' solver_density_cuda/cuda_forge/transportMix_d.cuh`
- `sed -n '1,83p' solver_density_cuda/input/speciesTransportDB.hpp`
- `sed -n '1,120p' solver_density_cuda/input/speciesTransportDB.cpp`
- `sed -n '1,120p' solver_density_cuda/cuda_forge/gasProperties_d.cu`
- `sed -n '1,90p' solver_density_cuda/cuda_forge/wmlesWallModel_d.cu`
- `sed -n '1180,1200p' solver_density_cuda/main.cpp`
- `sed -n '60,230p' plans/active/thermophysics-solver-owned-species-db.md`
- grep で `vis_lam`・`thermCond`・`cfg.visc` の利用箇所を調べてよい (`solver_density_cuda/cuda_forge/`, `solver_density_cuda/*.cpp`)
