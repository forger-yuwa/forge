# 諮問ブリーフ: 種ごとの輸送物性 (μ・λ) と CEA 形混合則の実装設計 (2026-09-27)

関連 plan: `plans/active/thermophysics-solver-owned-species-db.md` §4.3b・§4.3c・§5.1 #5t (② ③)・§10。
AGENTS.md エスカレーション 1 (§4 の設計方針) と 6 (`cuda_forge/` の数値の振る舞いを変える編集の前) に当たるので、実装前に諮る。

## 問い (1 つ)

下の「決定済み事項」を満たす実装設計 (下の素案) で進めてよいか。**穴があれば 1 つに絞って直し方を示し、実装を最小の段階に分けた順序と、各段階の合格条件 (事前に固定する数値) を 1 つ示してほしい。**
特に次の論点に答えてほしい (縛られなくてよい):
- (a) 既存の `viscMethod: 2` (LJ + Wilke/Mason–Saxena 共用 φ) の run をどう扱うか: 新しい `viscMethod: 3` を足して 2 を旧経路として残すか、2 を置き換えて種ごとの指定を必須にするか (既存 config は起動時エラーで案内)。
- (b) 種ごとの指定をどこに書くか: config の `physProp.transport: {種: モデル}` か、species データのエントリ側か。lump の構成種の指定の持ち方。
- (c) GPU 実装の精度と構造: 種ごとの評価 (CEA フィット・IAPWS・Chapman–Enskog)、相互作用 ηᵢⱼ の n² 評価、float/double の使い分け。
- (d) 解決済み記録・互換性ハッシュへの入れ方 (既存 TP の記録を変えないこと)。

## 決定済み事項 (ユーザ決定 2026-09-27; plan §4.3c)

1. **混合則は必ず CEA 形 (frozen)** にする: φᵢⱼ = 2Mⱼ μᵢ/[(Mᵢ+Mⱼ) ηᵢⱼ]、λ には ψᵢⱼ = φᵢⱼ {1 + 2.41 (Mᵢ−Mⱼ)(Mᵢ−0.142 Mⱼ)/(Mᵢ+Mⱼ)²}、
   μ = Σ Xᵢ μᵢ / Σⱼ φᵢⱼ Xⱼ、λ = Σ Xᵢ λᵢ / Σⱼ ψᵢⱼ Xⱼ (`.venv-cea/nasa_cea/cea2.f:5599-5634`)。ηᵢⱼ は `trans.inp` の相互作用データ、無い組は CEA の剛体球近似 (`cea2.f:5565-5570`)。
   **A/B 済み**: 種別 μ・λ を CEA にそろえて N2–H2O 16 状態 (T 400/600/1000/2000 K × X_H2O 0/0.1/0.5/1) を FCEA2 と比較し、CEA 形は最大 μ 0.001 %・λ 0.013 %、forge 現行 (Wilke φ を λ にも共用) は最大 μ 13.1 %・λ 10.7 %
   (`notes/investigations/2026-09-27-cea-vs-forge-properties/mixing_ab.py`・`mixing_ab_result.txt`・`fcea_mix/`)。
2. **μ・λ の出所は種ごとに必ず指定する** (全体モード・既定値は廃止)。候補: `cea` (trans.inp) / `kinetic` (LJ + Chapman–Enskog; 極性分子は双極子補正) / `fit` (ユーザ指定の CEA 形フィット) / H2O の合成。
   指定の無い種は起動時エラー。config は runner が problem から生成する。
3. **kinetic モードでも混合則は CEA 形**、ηᵢⱼ は LJ の二元 Chapman–Enskog (組み合わせ則 σᵢⱼ = (σᵢ+σⱼ)/2, εᵢⱼ = √(εᵢεⱼ))。kinetic の λ は Warnatz 式 (Chemkin/Cantera 系: 内部自由度分離 + Parker の回転緩和数) を候補
   (当方検算: N2・O2・CO2 で CEA と ±6 %; 修正 Eucken は CO2 で −13 %; `h2o_blend_and_lambda.py`)。
4. **H2O の μ・λ**: 600 K より下は IAPWS 希薄気体の式 (粘性 2008 μ₀、熱伝導 2011 λ₀; 公式適用域 253–1173 K)、上は CEA。**500–700 K で log μ・log λ を smoothstep (3s²−2s³) で移す** (この範囲の CEA/IAPWS 差は 0.4 % 以内)。
   253 K 未満は IAPWS の外挿 (μ は ~134 K、λ は ~100 K まで単調) として記録。H2O はどの簡易 kinetic λ でも合わない (Warnatz +9〜+43 %) のでデータを使う。
   H2O–N2 などの相互作用 ηᵢⱼ は trans.inp のもの (H2O–N2 は 300–1000 K 区間、300 K 未満は外挿)。
5. 拡散係数は LJ (二元 Chapman–Enskog、lump は Blanc) を継続し今回は変えない。

## 観測事実 (現行コード)

- `solver_density_cuda/cuda_forge/gasProperties_d.cu:52-85`: `viscMethod` 0 定数 / 1 Sutherland / 2 kinetic。2 は `roY` から Y→X を作り `thermo_mu_mix`・`thermo_lambda_mix` (double 入力を float で評価) を呼ぶ。
- `solver_density_cuda/cuda_forge/thermo_d.cuh:379-452`: 単成分 `thermo_mu_species` (Chapman–Enskog, Neufeld Ω22)、`thermo_lambda_species` (修正 Eucken)、Wilke `thermo_wilke_phi` を μ と λ で共用。`THERMO_MAX_SPECIES 16` (`:35`)。
  壁モデル `wmlesWallModel_d.cu:53-54` も同じ関数を呼ぶ。
- 拡散 `thermo_Dbinary`/`thermo_Dmix_species(_f)` (`thermo_d.cuh:456-575`) は `speciesTransport_d.cu:259` から。
- 種データ: 共通データ `solver_density_cuda/data/species/forge_species_v1.yaml` (ビルド時埋め込み; LJ 欄あり) と外部 `speciesDBFile`。lump はソルバが起動時に構成種から合成 (#6a 完了) し、LJ は暫定で質量分率平均 (#7 で構成実種へ展開予定)。
- 解決済み記録 `resolved_species_*.yaml` と互換性ハッシュ (#3a) があり、lump が無い config のハッシュ本文はバイト不変という約束で運用している。
- 凝縮域では液を除いた気相組成で μ・λ を評価する別 plan がある (`plans/active/condensation-two-phase-transport.md` §4.1)。組成を作る場所 (`gasProperties_d.cu:71-83`) を両 plan が触る。

## 素案 (当方)

- 共通データに CEA `trans.inp` の種別フィット (V/C, 区間) と相互作用 ηᵢⱼ を取り込み、LJ に加えて双極子モーメント・幾何 (原子/線形/非線形)・Z_rot(298) を持たせる。外部 `speciesDBFile` でも同じ欄を書ける (ユーザ定義 `fit`)。
- config: `physProp.transport: {N2: cea, H2O: iapws_cea, CO: kinetic, ...}` (必須、lump は構成種名で書く)。新しい `viscMethod: 3` で有効化し、`viscMethod: 2` は旧経路として残す (既存 run の再現)。
- ソルバ起動時に、輸送用の「実種リスト」(lump を構成種へ展開、輸送種 s → 実種 k の分率行列) と、種ごとの評価子の表・ηᵢⱼ の組の表を作り GPU へ置く。`gasProperties_d` は Y (輸送種) → 実種 X → 種別 μᵢ・λᵢ・ηᵢⱼ → CEA 混合、を double で評価して float で格納 (現行と同じ精度方針)。
- 記録・ハッシュ: `viscMethod: 3` のときだけ輸送の指定・係数を記録と互換性ハッシュに追記 (それ以外はバイト不変)。
- 検証: ① Python 参照 (`mixing_ab.py` の mixB) との double 一致 ≤1e-12・float ≤1e-5 (16 状態 + lump + H2O つなぎ区間の両側)、② FCEA2 との一致 ≤0.1 % (既存 16 状態)、
  ③ H2O つなぎの C¹ 連続 (500/700 K の両側で値・傾き)、④ lump 展開 = full 実種の結果、⑤ NS run 1 本 (case/16 Wysłouzil の SST 凝縮, 収束済み場から) で変化量を記録 (不変を合格にしない; 収束・準定常 VERDICT)。

## 禁止事項 (厳守)

- ファイルを変更しない。**`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md`, `thermo.inp` 全体・`trans.inp` 全体を読まない**。下の範囲以外のファイル読みをしない (grep は可)。
- 推奨は 1 つに絞る。根拠は `ファイル:行` か本ブリーフで示す。

## 読んでよいもの

- `sed -n '370,460p' solver_density_cuda/cuda_forge/thermo_d.cuh`
- `sed -n '40,90p' solver_density_cuda/cuda_forge/gasProperties_d.cu`
- `sed -n '40,65p' solver_density_cuda/cuda_forge/wmlesWallModel_d.cu`
- `sed -n '5595,5640p' .venv-cea/nasa_cea/cea2.f`
- `sed -n '5530,5572p' .venv-cea/nasa_cea/cea2.f`
- `sed -n '1,120p' notes/investigations/2026-09-27-cea-vs-forge-properties/mixing_ab.py`
- `sed -n '1,80p' notes/investigations/2026-09-27-cea-vs-forge-properties/h2o_blend_and_lambda.py`
- `sed -n '60,200p' plans/active/thermophysics-solver-owned-species-db.md`
- `sed -n '1,60p' solver_density_cuda/data/species/forge_species_v1.yaml`
- `sed -n '40,90p' plans/active/condensation-two-phase-transport.md`
