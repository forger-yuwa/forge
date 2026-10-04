# 諮問ブリーフ: 輸送物性 (μ・λ・D) の正本と解決順の設計 (2026-09-27)

関連 plan: `plans/active/thermophysics-solver-owned-species-db.md` §4.3b・§5.1 #5・#7。

## 問い (1 つ)

ユーザ方針「CEA を正解とする」と、ユーザ承認済みの案「`species_db.yaml` にも μ・λ のフィット係数を直接書ける欄を設け、解決順を
**ユーザ指定フィット → CEA `trans.inp` → LJ (最後の手段)** とする。拡散係数は LJ (または Schmidt 一定) のまま」について、
**この設計で進めてよいか、穴があれば 1 つに絞って直し方を示してほしい**。特に下の論点 (a)〜(e) に答えてほしい (縛られなくてよい)。
あわせて、実装と検証を最小の段階に分けた順序を 1 つ示してほしい。

- (a) 混合則: CEA は相互作用データ (41 組, 例 `H2O–N2` の相互作用粘性) を持つ。forge の Wilke / Mason–Saxena にどう入れるか (CEA 本体の混合則に合わせるべきか)。
- (b) CEA の範囲外: **H2O は 373.2 K 以上しかない** (`trans.inp:290-296`)。凝縮域 200–300 K の水蒸気をどう扱うか (延長 / LJ に落とす / 別文献)。区間端での連続性。
- (c) 拡散係数: CEA に無い。μ・λ を CEA にし D を LJ のままにすると、Schmidt・Lewis 数の整合が崩れないか。相互作用粘性から D を作る道はあるか。
- (d) lump: 構成実種に展開して混合する方針 (plan §4.4, #7) と、CEA の種別フィットの組み合わせ。
- (e) 検証: 何を正解として合否を決めるか (CEA 本体の混合物 μ・λ 出力との一致? `FCEA2` がローカルにある)。

## 観測事実

- 現行の実装 (`solver_density_cuda/cuda_forge/gasProperties_d.cu:52-85`): `viscMethod` 0 = 一定、1 = Sutherland (空気の係数固定)、
  2 = 種ごとに Chapman–Enskog (`thermo_mu_species`, `thermo_d.cuh:379`; LJ σ・ε、Neufeld Ω22) → Wilke (`thermo_mu_mix`, `:418`)、
  λ は修正 Eucken λ=μ(cp+1.25R) (`thermo_lambda_species`, `:391`) → Mason–Saxena (`thermo_lambda_mix`, `:436`)。
  拡散は二元 Chapman–Enskog (`thermo_Dbinary`, `:456`) → 混合平均 (`thermo_Dmix_species`, `:472`) か Schmidt 一定。float 版 (`_f`) もある。
- 種ごとの物性で輸送に効くのは LJ の 2 つだけ。ユーザが `species_db.yaml` で種を定義しても μ(T)・λ(T) を直接与える手段は無い。
- CEA `.venv-cea/nasa_cea/trans.inp`: 66 種・相互作用 41 組。形式 `ln η[μP] = A lnT + B/T + C/T² + D`、λ[μW/(cm·K)] 同形、種ごと 2〜3 区間
  (例 N2 200–1000–5000–15000 K、H2O 373.2–1073.2–5000–15000 K)。相互作用は V (粘性) のみの組が多い (`H2O–N2` は V3C0)。拡散係数は無い。
- 現行 forge (単成分, LJ + Chapman–Enskog + 修正 Eucken) と CEA の比 (当方の Python 検算, 2026-09-27):

  | 種 | μ 200/300/600/1000/2000 K | λ 300/600/1000/2000 K |
  |---|---|---|
  | N2 | +1.5 / +0.9 / +0.3 / −0.4 / −0.9 % | −2.5 / −2.5 / −1.5 / −4.6 % |
  | O2 | +0.5 / −0.5 / −1.6 / −2.4 / −3.0 % | −3.5 / −3.9 / −3.8 / −6.8 % |
  | Ar | +2.1 / +1.8 / +1.1 / +0.1 / −1.3 % | +1.5 / +0.7 / −0.2 / −1.7 % |
  | CO2 | −0.7 / +0.4 / −0.5 / −0.7 / −1.2 % | −1.8 / −9.5 / −11.2 / −13.0 % |
  | H2O | (CEA 範囲外) / (同) / **+22.9 / +13.6** / +5.0 % | (範囲外) / **+47 / +26** / +3.4 % |

  H2O の LJ (2.605 Å, 572.4 K) は極性分子の値を双極子補正なしで使っているためと推定 (未検証)。
- 影響する run: `viscMethod: 2` の NS run (case/16 Wysłouzil の N2+H2O SST 凝縮、case/42・45 の NS 凝縮、SERN の燃焼生成物など)。Euler (visc 0) は影響なし。
  凝縮域では液を除いた気相組成で評価する方針 (`plans/active/condensation-two-phase-transport.md` §4.1) と組み合わさる。
- 共通データ `solver_density_cuda/data/species/forge_species_v1.yaml` (#4 で導入, ビルド時埋め込み) に LJ 欄がある。ここに CEA の輸送フィットも入れる想定。

## 期待値と出典

- ユーザ方針: 「CEA が正解と考えて実装してくれたらいい」(2026-09-27)、「species_db.yaml で直接物性を指定することもできる」を維持。
- CEA の輸送: Gordon & McBride NASA RP-1311 (1994/1996) の transport 部、`trans.inp` の各エントリの出典 (Svehla 1994, Bich 1990 など)。

## 実施済みの操作

- 上表の検算のみ。コード変更なし。

## 仮説 (棄却してよい)

- H1: μ・λ の正本を CEA フィットにし、ユーザ定義フィット → CEA → LJ の順で解決、混合は Wilke/Mason–Saxena のまま (CEA の相互作用粘性は初回は使わない) で十分。
- H2: H2O の 373 K 未満は、CEA の最低区間の式を下へ延長するより、別の信頼できる相関 (IAPWS の蒸気粘性など) を入れるべき。

## 禁止事項 (厳守)

- ファイルを変更しない。**`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md`, `thermo.inp` 全体を読まない**。下の範囲以外のファイル読みをしない (grep は可)。
- 推奨は 1 つに絞る。根拠は `ファイル:行` か本ブリーフで示す。

## 読んでよいもの

- `sed -n '340,480p' solver_density_cuda/cuda_forge/thermo_d.cuh`
- `sed -n '40,90p' solver_density_cuda/cuda_forge/gasProperties_d.cu`
- `sed -n '1,30p' .venv-cea/nasa_cea/trans.inp`
- `sed -n '285,310p' .venv-cea/nasa_cea/trans.inp`
- `sed -n '390,410p' .venv-cea/nasa_cea/trans.inp`
- `sed -n '60,140p' plans/active/thermophysics-solver-owned-species-db.md`
- `sed -n '1,60p' solver_density_cuda/data/species/forge_species_v1.yaml`
- `sed -n '40,80p' plans/active/condensation-two-phase-transport.md`
