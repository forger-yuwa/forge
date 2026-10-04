# 諮問ブリーフ: lump を含む化学種拡散係数 (plan #7) の実装設計

エスカレーション条件 1 (§4.4 の未確定部分を確定させる)、6 (`cuda_forge/` の拡散係数の数値の振る舞いを変える)。
plan: `plans/active/thermophysics-solver-owned-species-db.md` §4.4、§5.1 #7、§6 V4。ユーザ指示 (2026-10-05): 「ガンガンやって」。

## 観測事実 (現状のコード)
- μ・λ: `viscMethod 2` は `physProp.transport` 必須で、lump は実種に展開済み (#5t2-2、`transportMix_d.cuh` `transport_expand_X`、モル基底)。
  旧 kinetic μ・λ 経路は廃止済み (`procedures/solver-settings.md` 「physProp.viscMethod」)。**lump の平均 LJ が残っているのは化学種の分子拡散係数だけ**。
- 拡散係数: `thermo_d.cuh` `thermo_Dmix_species_f` (float、補数形 D_i = Σ_{j≠i}X_j / Σ_{j≠i} X_j/D_ij、D_ij = Chapman–Enskog + Neufeld Ω(1,1))。
  呼び出しは面ごと: `speciesTransport_d.cu:328` (species_diffusion_d)、`:2228` (二相拡散の面入力)、`:2466` (独立監査)、`:2848`、probe `:997`。
  輸送種 s の LJ は `SpeciesThermoF.sigma_LJ/eps_kB`。lump の LJ は構成種の**質量分率平均** (`speciesDB.cpp:1047` に PROVISIONAL 表示)。
- 補正 `J_i* = J_i − Y_i ΣJ` と `Σh_i J_i*` は輸送種単位 (`speciesTransport_d.cu:271` 付近)。
- 構成種の情報は host の `ResolvedLump` (`input/speciesDB.hpp:36`: members・x (lump 内モル分率)・memberSpecies (LJ・MW 入り))。
  device には lump の構成は渡っていない (`physProp.transport` のときだけ輸送表側に展開行列がある)。
- 拡散は `viscMethod ≠ 0` かつ `speciesDiffusionMethod 1` (既定) で、`viscMethod 1` (Sutherland) + 多成分でも使われる → 輸送表の有無と独立に展開が要る。
- 使用例: SERN (`case/46`、lump `EXH` (CEA 凍結組成) と `AMB` (空気) の 2 lump)、case/44 va3 (`MIXDRY` lump + H2O)、case/16 (N2 + H2O、lump なし)。

## 設計案 (主セッション)
**実種展開の混合平均**: 輸送種のモル分率 X_s → 実種のモル分率 X_r = Σ_s X_s·E[s,r] (E = lump 内モル分率、非 lump は単位行; 重複実種は加算)。
各実種 r の混合平均 D_r = Σ_{q≠r} X_q / Σ_{q≠r} X_q/D_rq (補数形を維持)。
- 非 lump の輸送種 i (実種 r(i)): **D_i = D_{r(i)}**。lump の外の種 j から見た lump 内寄与は Σ_{q∈L} X_L x_q/D_iq = X_L/D_{i,L}^{Blanc} なので §4.4 の Blanc と同じ。
  r(i) が lump の構成にも含まれる場合 (例: H2O を輸送種と EXH の構成の両方に持つ) も、同じ分子として全量 X_r で評価するので重なりの定義が自然に決まる。
- lump の輸送種 L: lump 内の質量分率が固定なら J_L = Σ_{q∈L} J_q = −ρ (Σ_q y_{q|L} D_q) ∇Y_L (構成種 q が L にしか無いとき厳密) → **D_L = Σ_{q∈L} y_{q|L} D_q** (lump 内質量分率重み)。
  §4.4 の「lump 同士の二重和」を別に作らず、この定義で lump 同士の相互拡散も表す。構成種が他の輸送種と重なるときは近似 (記録する)。
- lump の無い config では E = 単位行列で、現行の `thermo_Dmix_species_f` と同じ演算順にしてビット不変を保つ (分岐で現行関数を呼ぶ)。
- device: 実種の (MW, σ, ε) と E を `thermo_init_db` で上げる (n_real ≤ 32)。
- コスト: 面ごとに n_real² の D_rq (powf/expf)。SERN の実種数は ~6〜10 で現行 (2 輸送種で 2 回) の数十倍。対案: (i) D_rq·P を ln T の区分 3 次表で引く (輸送表 #5t2-3 と同じ方式)、(ii) 節点で D を作って面は平均。

## 検証案 (§6 V4 の具体化、事前登録用)
- V4a 単体: lump を含む config と、同じ実種組成を `full` (全実種を輸送種) で書いた config で、外部種の D_i が double 参照で ≤1e-12、float で ≤1e-6 一致。lump の D_L は定義どおり。
- V4b 不変: lump の無い config (case/16 湿り凝縮 ON/OFF、乾き) は旧バイナリ ×3 vs 新 ×2 の `check_field_regress` PASS (ビット一致を期待)。
- V4c 変化量の記録 (不変を合格にしない): case/44 va3 (MIXDRY + H2O) の H2O の D の変化 (検算 −0.7〜−1.4 %) と報告量の変化。SERN 入口組成の D_i,mix 変化表 (SERN へ連絡、SERN 側の回帰で確認)。
- V4d 補正後流束の lump vs full の差 (codex 検算 +0.50 %) を記録。

## 問い
1. 実種展開 (D_i = D_{r(i)}) と lump の D_L = Σ y_{q|L} D_q の定義は §4.4 (Blanc・縮約拡散・lump 同士・重なり) の要件を満たすか。モル重みでなく質量重みでよいか (J は質量流束、勾配は ∇Y)。重なりの扱いに穴は。
2. ビット不変の保ち方 (lump 無しは現行関数) は妥当か。
3. コストへの対処: 表引き・節点評価・そのまま、のどれを今やるべきか。精度基準は。
4. V4a〜d は十分か。SERN への影響の扱い (既定変更として §9 に旧挙動の明示キーを残すか、残さず記録だけか)。

読んでよいファイル: plan 上記、`solver_density_cuda/cuda_forge/thermo_d.cuh` (455–700 付近)、`speciesTransport_d.cu` (200–340・2200–2480)、
`input/speciesDB.{hpp,cpp}`、`cuda_forge/transportMix_d.cuh`・`transportTables_d.cuh`、`procedures/solver-settings.md`。
