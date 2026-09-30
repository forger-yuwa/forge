# 化学種の熱物性をソルバが持ち、run には組成 (モル分率) だけを書く

## メタ

- **area**: `thermophysics`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/thermophysics.md`](../../methods/thermophysics.md) (種 DB・擬似種・モル分率入力の現在仕様)
  - [`methods/condensation.md`](../../methods/condensation.md) (§潜熱 L(T)、§L(T) が液相の熱力学を決めている)
  - [`procedures/solver-settings.md`](../../procedures/solver-settings.md) (bcond `X{s}`、`speciesDBFile`)
- **related_plans**:
  - 前身: [`thermophysics-cea-mole-fraction-species.md`](../accepted/thermophysics-cea-mole-fraction-species.md) (モル分率入力・`full|lumped`。「`thermo_d.cu` 内蔵 DB の変更」をスコープ外にした結果、設計側が DB を合成して run に書く構造になった — 本 plan でその判断を改める)
  - 関連: [`tooling-design-problem-campaign-recipe.md`](tooling-design-problem-campaign-recipe.md) (problem には組成・lump の指定だけを置き、生成された数値を run に持ち込まない方針と同じ向き)
- **created**: `2026-09-27`
- **owner**: `Claude (主セッション) / ユーザ`

## 1. 目的

ユーザ要望 (2026-09-27): 「problem から `species_db.yaml` が生成されて、わけのわからない値 (擬似種の NASA-9 係数・`atoms` の平均原子数) が入るのはやめてほしい。
CEA の情報はソルバが持ち、モル分率で指定できるようにしたい。外部 DB で上書きできるのは残してよい」。
完了時には、run の入力は**種名と lump の中身のモル分率** (と境界組成のモル分率) だけになり、熱物性の正本はソルバ配布物の CEA 由来データ 1 か所になる。
ソルバは使った物性を**解決済み記録として出力**し (入力ではない)、restart はその内容で照合する。

## 2. スコープ

- **やる**:
  - CEA (`thermo.inp`, McBride–Gordon 2002) 由来の**版とハッシュを固定した共通データ**をソルバ配布物に置き、C++ (ソルバ・変換器) と Python (設計・後処理) が同じものを読む。外部 DB による上書きは残す。
  - config で lump を「名前 + 構成種のモル分率」で書けるようにし、係数は**ソルバが起動時に合成**する (温度区間は構成種の全区切りの和集合で分割)。
  - 使用した物性の解決済み記録 (ソルバ出力) と、**内容による** restart 照合。
  - lump の輸送物性 (粘性・熱伝導) は構成種へ展開して実種の Wilke / Mason–Saxena 混合で評価する。
  - `atoms` は入力から外す (共通データと解決済み記録には残す)。
  - 設計 runner は `species_db.yaml` を生成しない (lump の指定を config に書く)。
- **やらない** (初回):
  - `species_db.yaml` を先に消すこと (照合の仕組みができるまで残す; codex 2026-09-27 諮問)。
  - 内蔵化と同時に CEA の版・MW・外挿規約を変えること (差は別項目で判断する; §5.1 #4)。
  - `full` の強制 (`full|lumped` の選択は前身 plan の既存判断, `thermophysics-cea-mole-fraction-species.md:155`)。
  - 液相 H2O の物性モデル (密度・表面張力・飽和圧・273.15 K 未満/373.15 K 超の延長規約) の変更。
  - 移設と同時に N2 等の 6000–20000 K 区間を有効化すること (外挿規約の変更なので別判断; §5.1 #5)。

## 3. 関連 docs と前提 (観測事実, 2026-09-27)

根拠の詳細はブリーフ [`notes/reviews/briefs/2026-09-27-solver-owned-species-db.md`](../../notes/reviews/briefs/2026-09-27-solver-owned-species-db.md)。

- **CEA 係数が 4 か所**: ソルバ内蔵 7 種 (`solver_density_cuda/input/speciesDB.cpp:83-157`)、設計側 `SPECIES_NASA9` 11 種 (`design/forge_design/gas/semiperfect.py`)、
  `cea_thermo_to_species_db.py` の thermo.inp 直読み、凝縮潜熱 `h2o_latent` の H2O 気相再ハードコード (`condensationProperties_d.cuh:239-262`)。
  既知の不一致: H2O MW 0.0180153 vs 0.01801528、AR 高温 a0 0 vs 20.105。
- ソルバの DB 読み込みは `MW, LJ, Tlo/Tmid/Thi, nasa9_low/high` のみで `atoms` は読まない (`speciesDB.cpp:173-190`)。温度区間は種ごとに持てるが **2 区間まで** (`thermo_d.cuh:46-50`)。
- `atoms` は設計側の元素質量分率診断 (`composition.py:248-255`) 用で、呼び出しは単体テストだけ。lump の値は構成種のモル加重平均 = 平均原子数 (誤りではない)。
- lump 合成 (`composition.py:259-280`) は NASA-9 のモル加重線形混合 (固定組成なら cp/h 厳密)。区切りが 200/1000/6000 K でない種は拒否。
  LJ は質量分率の単純平均で、ソルバの Wilke 混合 (`thermo_d.cuh:346-435`) とは一致しない。codex の物性式検算で平均 LJ の粘性は実種 Wilke に対し **200 K で −1.30 %、300 K で −0.98 %** (CFD 実測ではない; Euler では不使用)。
- `s°` は等エントロピー変換で使われる (`thermo_d.cuh:801`)。固定組成では混合エントロピー項が差分で消えるので lump で省略してよいが、区間・外挿規約は保持する (codex Minor)。
- **restart の種署名は内蔵種の係数変更を検出しない**: `solver_density_cuda/tools/forge_species.py:189` が内蔵種の係数・区切りを `None` にし、`:221` が両側 builtin の比較を省く。
  codex が比較関数を単独実行して確認 (内蔵 `N2.nasa9_low[2]` 0→1 で不一致リスト `[]`、file 由来なら検出)。`interp_field.py:67` の拒否判定もこれを使う。
- 液相 H2O はソルバにハードコード: 飽和圧 Murphy–Koop、密度 `1000−0.12(277−T)` (下限 920)、表面張力 IAPWS 外挿、潜熱 `L=h_v−h_l` (h_l は CEA H2O(L)、273.15 K 未満は cp_l 4228 一定、373.15 K 以上はクランプ)、
  `e_l = e_v + R_v T − L`。H2O の液比熱は h_l の傾きで決まる (潜熱からの逆算ではない) が、`h2o_latent` の h_v が種 DB と別ソース・別外挿なので、EOS が暗黙に使う液エンタルピーに差が入る (codex Major; 120 K で 2.39 kJ/kg, `methods/condensation.md:587-590`)。
  N2 (CPG 経路) は液物性を持たず `c_l = c_p,v − dL/dT` が L のフィットから暗黙に決まる。
- 既定の化学種拡散は `speciesDiffusionMethod = 1` (kinetic 混合平均, `solverConfig.hpp:579`) で、二元拡散係数も LJ を使う (`thermo_d.cuh:452-462`)。種流束の補正と `Σh_s J_s` のエネルギー項は `speciesTransport_d.cu:254` 付近。
- `species_db.yaml` を読む後処理・restart reader がソルバ外にもある: `solver_density_cuda/tools/total_quantities.py:119` (無ければ `species_db.yaml` を開く, `_TPGas` は 2 区間前提 `:36-48`)、
  種変換器 (`convert_species_field.py` が `_TPGas` を import)、SERN runner の `_species_signature` (`design/forge_design/evaluate/runner_sern.py:571`)。
- 現状の `restart_field.py` は種署名を見ずに配列をコピーし (`restart_field.py:36`)、ソルバは `valueFileName` を直接読む (`main.cpp:1193`)。
- CEA `thermo.inp` は 2,012 エントリあり、`CO` (MW 28.0101) と `Co` (58.9332) のように大小文字だけが違う別種がある。Python は種名を大文字化 (`composition.py:35`)、C++ も大小文字を同一視 (`speciesDB.cpp:75`)。
  `H2O(L)` など凝縮相も含む。LJ 表は 23 種分しかなく、現行生成器は LJ 不明種に N2 相当値を仮置きする (`cea_thermo_to_species_db.py:105`)。
- float の熱物性は `SpeciesThermoF` と専用評価関数 (`thermo_d.cuh:64`)、datum は `low[7]`/`high[7]` の 2 か所に焼き込み (`thermo_d.cu:69`)。
- 方向の判断: codex (diagnose) に諮った — [`notes/reviews/2026-09-27-solver-owned-species-db-diagnose.md`](../../notes/reviews/2026-09-27-solver-owned-species-db-diagnose.md)
  — 結論「移行する。(b) を温度区間の和集合に拡張した起動時係数合成。最初の変更は内蔵種を含む解決済み物性の保存・内容照合」。採否は同記録末尾。

## 4. 設計方針

### 4.1 熱物性の正本 (共通データ)

- `thermo.inp` から生成した**全種表** (区間数可変・MW・Hf・LJ・元素組成) を、CEA の版と `thermo.inp` の SHA-256 付きでソルバ配布物に置く (例 `solver_density_cuda/data/species/cea2002.yaml`、生成器は `cea_thermo_to_species_db.py` を拡張)。
  C++ の内蔵 DB と Python の `SPECIES_NASA9` はこれを読む (ハードコードを撤去)。LJ は CEA に無いので、出典付きで同じファイルに持つ。
- 外部 DB (`speciesDBFile`) による上書き・追加は残す。
- **種の ID は大小文字を区別する canonical ID** (CEA の表記) とし、互換の別名 (`AR`→`Ar`、`WATER`→`H2O` 等) は明示的な alias 表で持つ。大文字化による同一視はやめる (`CO`/`Co` の衝突; codex plan レビュー M3)。
- 各エントリに**相** (気相/凝縮相)・**熱力学の利用可否**・**輸送データ (LJ) の有無**を持たせ、凝縮相を気相 EOS に使う・LJ 無しの種を粘性/拡散に使う、は使用時に拒否する (N2 相当値の仮置きは廃止)。
- 現行の `AIR` (cp/R 3.5 一定の擬似空気) は CEA の `Air` と別の互換擬似種として残す。
- 移行の第一段では**現行内蔵値をそのまま共通データに移す** (値を変えない)。CEA 直読みとの差 (H2O MW・AR a0) をどちらに寄せるかは §5.1 #4 で別に判断する。

### 4.2 config での lump 指定と起動時合成

```yaml
physProp:
  species:
    - {name: MIXDRY, lump: {N2: 0.708873, O2: 0.230376, AR: 0.00850387, CO2: 0.0522474}, basis: mole}
    - H2O
```

- ソルバは起動時に lump の係数を合成する: 構成種の**全温度区切りの和集合**で区間を分け、区間ごとにモル加重で NASA-9 係数を足す (NASA-9 は係数に線形なので厳密)。
  端の外挿 (定 cp) も区間として表す。`SpeciesThermo` を 2 区間固定から**区間数可変**に拡張する (上限は定数; GPU 側は区間表のオフセットで持つ)。
- datum (`thermoHrefTemp`) は**全区間**に適用する (現行は `low[7]`/`high[7]` の 2 か所)。float 表 (`SpeciesThermoF`) も同じ区間表から作る。
- 輸送種数・lump 展開後の実種数・区間数の上限を分けて定数で持ち、超過は起動時に拒否する。
- 起動ログに合成結果 (MW・区間・参照温度での cp/h) と lump の中身を出す。
- 凝縮種は lump に入れられない (独立種として残す)。
- 境界の組成は従来どおり種 (lump を含む) のモル分率 `X{s}` で書く。将来は実種のモル分率を書いて lump へ写す入力も検討 (初回はやらない)。

### 4.3 解決済み記録と内容照合 (最初の一歩)

- ソルバは使用した全種 (内蔵種を含む) の**解決済み物性** (係数・区間・MW・datum・lump の中身・データの版とハッシュ) を run に**出力**する (例 `resolved_species.yaml`)。これは入力ではなく記録。
- **記録を保存場に結び付ける** (codex plan レビュー M1): 使用物性の内容ハッシュを**各 `res_*.h5` の属性**に書き、解決済み記録 (ハッシュ名で不変に保存、上書きしない) と対応させる。
  lump の構成実種の係数・LJ も記録対象。
- 宛先 (これから回す run) の物性は**起動前に同じ resolver で解決**し、次の全入口で保存場のハッシュと照合する: ソルバの `valueFileName` 直接読込、
  同一メッシュ restart (`restart_field.py`)、補間 (`interp_field.py`)、種変換 (`convert_species_field.py`)、設計 runner の段間引き継ぎ・warm start。
- **ハッシュは誰が付けるか** (2026-09-27 codex diagnose [`notes/reviews/2026-09-27-species-hash-attachment-diagnose.md`](../../notes/reviews/2026-09-27-species-hash-attachment-diagnose.md) — 補強した案 A を採用):
  属性は**その保存量を実際に生成した処理**が付ける。新規初期場 = IC 生成処理 (解決済み物性・datum で保存量を作り、成功時に付与; 既存場への後付け認証はしない)、
  ソルバ保存場 = ソルバ、コピー・restart・補間 = SRC の記録を検証し宛先と互換を確認してから継承 (SRC が未検証なら DST の既存属性を消す; 宛先のハッシュで埋めない)、
  種変換 = 入力の検証と変換の成功後に変換先のハッシュを付与。変換器 (`convertGmshToForge`) の属性で完成した TP 場を証明しない (変換器の解決と IC のエネルギー生成は別経路, `ic.py:68-85`)。
- **GPU 不要の resolve-only モード** (`forge --resolve-species` 相当) を #3 に含める。IC 生成・runner・種変換はこれで宛先の物性とハッシュを得る (V1 (d))。resolve-only の結果を既存場へ貼るだけで検証済みにしない。
- **属性・対応記録の無い場は既定で拒否**。許可は**その呼び出しだけ**の明示指定 (環境変数など) とし、生成 config に恒常的な許可を埋め込まない。未検証入力だった履歴は記録に残す。
- **ハッシュを 2 本に分ける**: 互換性ハッシュ (種の順序・canonical ID・相・MW・絶対基準の全係数と温度区間・外挿規約・LJ・lump の正規化済み構成と構成実種の物性・`thermoHrefTemp` と datum 規約・スキーマ版) と、
  記録全文の完全性ハッシュ。`source`・ファイルパス・配布 DB 全体の版は来歴として記録するが互換性ハッシュには入れない。CPG は照合対象外。液相は #10 まで含めない (それまで V1 (e) は未達)。
- 比較は `source` (builtin/file) でなく**内容**で行う (`forge_species.py:189,221` の builtin 省略を撤廃)。過去 run の署名を現在の内蔵表から再生成しない。
  記録が無い旧 run は「照合不能」として扱い、明示フラグでのみ許可する。

### 4.3b 輸送物性の正本は CEA `trans.inp` (2026-09-27 ユーザ決定)

- ユーザ方針: 「CEA が正解と考えて実装する」。CEA の熱物性 (`thermo.inp`) に LJ は無く、輸送物性は `.venv-cea/nasa_cea/trans.inp` に**粘性係数・熱伝導率の温度フィット**
  (`ln η[μP] = A lnT + B/T + C/T² + D`, `ln λ[μW/(cm·K)]` 同形; 種ごとに 2〜3 区間) として 66 種・相互作用 41 組がある (出典 Svehla 1994, Bich et al. 1990 など)。拡散係数は無い。
- 現行 forge (LJ + Chapman–Enskog + 修正 Eucken) との照合 (当方 2026-09-27、単成分):
  μ は N2 +1.5〜−0.9 %、O2 +0.5〜−3.0 %、Ar +2.1〜−1.3 %、CO2 ±1.2 % (200–2000 K)、**H2O +22.9 % (600 K)・+13.6 % (1000 K)・+5.0 % (2000 K)**;
  λ は N2 −1.5〜−4.6 %、O2 −3.5〜−6.8 %、Ar ±1.7 %、**CO2 −9.5〜−13 % (600–2000 K)**、**H2O +47 % (600 K)・+26 % (1000 K)**。
  H2O は極性分子の LJ 値を双極子補正なしで使っているためと推定 (未検証)。
- **CEA の H2O は 373.2 K 以上しかない** (凝縮域 200–300 K の水蒸気は範囲外)。CEA に無い種・範囲は LJ に落とす/延長規約を決める必要がある (§5.1 #5)。
- 方針: μ・λ は CEA `trans.inp` の種別フィットを正本とし、混合は Wilke / Mason–Saxena (CEA に相互作用データがある組はそれを使う)。拡散係数は CEA に無いので LJ (Blanc) のまま。
  NS の結果を変える数値変更なので、実装前に上位へ諮る (#5・#7)。
- **ユーザ定義の輸送フィット** (2026-09-27 ユーザ承認、codex 諮問中): `species_db.yaml` (外部 DB) と共通データに μ・λ のフィット係数を直接書ける欄を設ける (CEA と同じ形 `ln μ = A lnT + B/T + C/T² + D`、温度区間つき)。解決順は **ユーザ指定フィット → CEA `trans.inp` → LJ (最後の手段)**。拡散係数は LJ (または Schmidt 一定) のまま。
- **codex diagnose (2026-09-27, [`notes/reviews/2026-09-27-transport-source-design-diagnose.md`](../../notes/reviews/2026-09-27-transport-source-design-diagnose.md))**: 解決順は維持。
  **混合則も CEA の frozen 混合則全体を採用する** — forge は λ にも粘性用の φ を使う (`thermo_d.cuh:448`) が、CEA は相互作用粘性 ηᵢⱼ から φᵢⱼ = 2Mⱼμᵢ/[(Mᵢ+Mⱼ)ηᵢⱼ] を作り、
  λ には ψᵢⱼ = φᵢⱼ{1 + 2.41(Mᵢ−Mⱼ)(Mᵢ−0.142Mⱼ)/(Mᵢ+Mⱼ)²} を使う (`cea2.f:5613`; 相互作用データの無い組は CEA と同じ推定。`V3C0` でも ηᵢⱼ は ψ を介して λ に効く)。
  単成分値を両方 CEA にしても混合則だけで 600 K・X_H2O 0.5 で μ −11.4 %・λ −10.0 % (codex の double 検算; CFD・FCEA2 実行ではない)。
  比較の正解は FCEA2 の**同一 T・同一気相組成の μ と frozen λ** (平衡反応寄与込みの λ は別物, `cea2.f:5635,5738`)。
  **H2O の 373.2 K 未満を黙って LJ に切り替えない** (373.2 K で LJ は CEA より +31.6 %; CEA 最低区間の下方外挿は ~218.6 K 以下で dμ/dT<0; CEA 本体も最低区間を外挿する `cea2.f:5466`)。
  低温は μ・λ 両方の参照データ・許容誤差・接続規約が揃うまで「対応済み」にしない (§10)。拡散係数は LJ を継続し今回は変えない (Sc・Le が変わるのは不整合ではない)。
  lump は全実種へ展開してから CEA 混合を一度だけ行い、同じ実種が複数 lump に出れば分率を合算。μ・λ ごとの係数・単位・区間・範囲外規約を記録と互換性ハッシュに含める。

### 4.3c 輸送モデルの選択と混合則 (2026-09-27 ユーザ決定)

- **混合則は必ず変える**: forge の現行 (Wilke の φ を μ と λ で共用) をやめ、CEA 形の frozen 混合則 (相互作用粘性 ηᵢⱼ から φᵢⱼ、λ には ψᵢⱼ) にする。
  **両モード共通の混合則とし、kinetic モードでは ηᵢⱼ を LJ の二元 Chapman–Enskog (組み合わせ則で σᵢⱼ, εᵢⱼ) から作る** (2026-09-27 ユーザ承認)。
- ~~**μ・λ は 2 モードから選ぶ** (全体で cea / kinetic、既定は検証まで kinetic)~~ → **改訂 (2026-09-27 ユーザ決定): 種ごとに輸送物性の出所を必ず指定する** (全体モードと既定は廃止)。
  出所の候補: `cea` (trans.inp) / `kinetic` (LJ + Chapman–Enskog; 極性分子は双極子補正つき) / `fit` (ユーザ指定フィット) / H2O 用の合成 (例 `cea+iapws` = 案 B)。lump は構成種ごとの指定を使う。指定の無い種は起動時にエラー。config は runner が problem から生成する。混合則はどの組み合わせでも CEA 形 (FCEA2 で確認済み, #5t ①)。キー名・書式は #5t ② の設計で決める。
- **H2O の μ・λ (2026-09-27 ユーザ決定・確定)**: **600 K 付近でつなぐ** — 600 K より下は IAPWS 希薄気体の式 (粘性 IAPWS 2008 μ₀・熱伝導 IAPWS 2011 λ₀; 公式適用域 253–1173 K)、
  上は CEA `trans.inp`。**500–700 K で log μ・log λ を smoothstep (3s²−2s³) で IAPWS → CEA へ移す** (C¹ 連続; この範囲の CEA/IAPWS 差は μ・λ とも 0.4 % 以内)。
  高温は CEA の範囲内。**253.15 K (IAPWS の公式適用域の下端) 未満は、253.15 K での IAPWS の対数勾配に合わせた冪乗則で延長** (2026-09-27 ユーザ決定; 指数は式から計算: μ 0.7486・λ 0.9754、C¹)。IAPWS の式そのままの延長は μ₀ が 202.17 K で最小・134.12 K に極で非物理、kinetic+双極子は 253 K の傾き (1.08) が IAPWS (0.75) と合わず、CEA 下限 (373.2 K) からの延長は保証範囲内で IAPWS と μ −4.8 %・λ +7.3 % ずれるので採らない。200 K での 3 案の幅は μ で約 ±8 % だが、混合物 (水蒸気数 %) の μ への影響は 0.1 % 以下。
  (経緯: 同日「全温度 IAPWS」案を検討 — IAPWS の外挿は 4000 K まで CEA と μ ±1 %・λ ±3.6 % だったが、保証範囲を優先して 600 K つなぎに戻した。)
- **実装設計 (codex diagnose 2026-09-27, [`notes/reviews/2026-09-27-transport-implementation-design-diagnose.md`](../../notes/reviews/2026-09-27-transport-implementation-design-diagnose.md))** — 採用:
  - 選択は config (`physProp.transport`、実種ごとに必須)、係数は species データ (共通データ・外部 DB)。lump は展開後に同じ実種の分率を合算。
  - ηᵢⱼ の規約 (対称): 両種が `kinetic` なら二元 Chapman–Enskog、それ以外は CEA 相互作用データ、無ければ CEA 剛体球近似 (混在モデルの精度は未検証)。
  - GPU: μ・λ を同時評価、ηᵢⱼ は非対角組を 1 回ずつ、φᵢᵢ = ψᵢᵢ = 1。double で評価し float で格納 (**現行は内部 float なので精度変更として扱う**)。セルと壁 (`wmlesWallModel_d.cu:53`) が同じ評価関数を使う。
  - 記録・ハッシュ: **新しい輸送指定を使うときだけ**新スキーマで輸送ブロックを追記 (既存 TP の本文はバイト不変)。種別係数に加え、解決後の ηᵢⱼ の出所・係数、混合則の版、H2O の接続温度、外挿規約、展開行列を含める。
  - **合格条件を 2 系統に分ける** (Major): 「全種 CEA 指定」は FCEA2 と ≤0.1 %、「種ごとに選んだモデル」(H2O の IAPWS 接続など) は同じ MW・係数・接続規約の独立評価と double ≤1e-12。
    H2O を IAPWS にすると 400 K 純 H2O で CEA 比 μ +0.57 %・λ −2.26 % が正しい値なので FCEA2 基準では判定しない。参照側の H2O MW (18.01528) と DB (18.0153) の差も揃える。
  - ~~`viscMethod: 2` の扱い~~ → **決着 (2026-09-27 ユーザ, 案 C)**: 実装は `viscMethod: 2` を新方式 (種ごとの指定必須・CEA 形混合則) に置き換える形で進めるが、**既定の切り替え (種別指定なしを起動時エラーにする時期) は並行セッションの区切りを待つ** (照合の厳密化 #3c と同じ扱い)。それまでは旧経路が動く状態を保つ。2 の利用: case/27・28・50・55・56 (50/55/56 は並行セッションが作業中) と SERN runner。旧結果の再現は旧バイナリで。
  - **特別なモデルは `custom:` 名前空間** (2026-09-27 ユーザ決定): 標準の出所 `cea` / `kinetic` / `fit` と区別し、特定の種のために作った合成モデルは `custom:<名前>_v<版>` と書く。H2O は `custom:h2o_iapws_cea_v1` (600 K より下 IAPWS 希薄気体、上 CEA、500–700 K smoothstep)。定義 (接続温度・式) を変えたら版を上げる。版は記録と互換性ハッシュに入る。対象外の種に使えば起動時エラー。
- **kinetic モードの λ** (当方の検算, 同スクリプト): 修正 Eucken は CO2 で最大 −13 %。Warnatz 式 (Chemkin/Cantera 系: 内部自由度分離 + Parker の回転緩和数) は N2・O2・CO2 で ±6 % 以内。
  **H2O はどちらの簡易式でも合わない** (Warnatz +9〜+43 %) — 極性分子の熱伝導は簡易 kinetic では表せないので H2O の λ はデータ (案 B) を使う。kinetic は CEA にデータの無い種の予備で、λ は Warnatz 式を候補とする (#5t ② で決定)。
- **H2O の低温は特別なケアが要る可能性** (ユーザ): CEA モードで 373 K 未満をどう扱うか (CEA 外挿 / IAPWS つなぎ = 案 B) は §10 の未決事項のまま。比較は `notes/investigations/2026-09-27-cea-vs-forge-properties/` と Artifact ページ。

### 4.4 lump の輸送物性

- lump は構成種の組成を保持し、粘性・熱伝導は**全実種に展開して** Wilke / Mason–Saxena 混合で評価する (lump 内で粘性を作ってさらに混ぜる方式は採らない)。
- **化学種拡散** (2026-09-27 ユーザ指摘で改訂): これまで lump の LJ は構成種の**質量分率平均** (`composition.py:276`, 根拠の記載なし) で、
  H2O の混合平均拡散は「H2O と平均 LJ の擬似分子」の二元係数で評価していた。lump の外の種 i と lump の二元拡散係数は、lump の内部組成が固定なら
  **Blanc の法則 $1/D_{i,\mathrm{lump}}=\sum_{j\in\mathrm{lump}} x_j/D_{ij}$ (x_j は lump 内モル分率) で構成実種から厳密に作れる** (混合平均拡散の式と整合)。
  平均 LJ との差は H2O–MIXDRY (case/44 va3 組成) で 200 K −1.4 %、300 K −1.2 %、1000 K −0.7 % (Chapman–Enskog + Neufeld Ω の比のみの検算, 当方 2026-09-27)。
  → lump と kinetic 混合平均拡散 (`speciesDiffusionMethod: 1`, 既定) は**併用可**とし、二元係数を Blanc で作る。失うのは lump **内部**の構成種どうしの差動拡散だけで、
  これは「lump 内の組成は固定」という lump の定義そのもの (必要なら `full`)。
  codex plan レビュー M2 の「lump + kinetic を拒否」は、lump 外の種との拡散まで失うとみなした過剰な制約として撤回 (§6.1)。
- **縮約拡散モデルとして明記する** (codex plan レビュー 2 回目 M1): Blanc で集約できるのは混合平均係数の分母までで、補正 `J_i* = J_i − Y_i ΣJ`
  (`speciesTransport_d.cu:271`) は lump 内各種の流束に依存するので、**補正後の外部種流束は `full` と厳密には一致しない**
  (codex 検算: va3 乾き組成・300 K・Y_H2O 0.04 で H2O の補正後実効係数が `full` 比 +0.50 %)。これは「lump 内組成固定」の近似誤差として扱い V4 で記録する。
  **lump 同士** (SERN の `EXH`/`AIR` のような複数 lump) の二元係数は、両 lump の構成実種どうしの二重和 $1/D_{AB}=\sum_{j\in A}\sum_{k\in B}\ldots$ の形を
  実装前に式として確定し、構成実種が lump 間で重なる場合 (自己拡散を含む) の扱いも定義する (§5.1 #7 の前提)。エネルギー流束 `Σh_i J_i*` も同じ縮約で評価する。
- 平均 LJ の擬似分子は廃止する。NS の結果は変わる (codex 検算で粘性 −1 % 級の是正) ので、変更量を記録する (不変を合格条件にしない)。

### 4.5 `atoms`

- run の入力には書かない。共通データ (CEA の元素欄) と解決済み記録には残し、元素診断はそこから作る。

### 4.6 Python 側の共通 API と reader の移行

- Python の物性解決・区間評価・署名比較を 1 つの共通 API (共通データを読む) に集約し、`total_quantities.py` (`_TPGas`)、`convert_species_field.py`、`forge_species.py`、
  SERN runner の `_species_signature`、設計側 `gas/` をこれに載せ替える (codex plan レビュー M4)。2 区間前提を撤廃する。

### 4.7 設計 runner

- `species_db.yaml` を生成しない。problem の組成・lump 指定を §4.2 の config に翻訳するだけにする。設計側 (MOC・IC) の熱物性も共通データを読むので、設計と CFD の熱力学が同じ正本から来る。

### 4.8 潜熱 (凝縮種の液相エンタルピー) — 2026-09-27 ユーザ方針

- 液相は**気相と同じ datum 定数でエンタルピーをシフト**し、気液差 (CEA の H2O と H2O(L) の絶対エンタルピー差) を保つ。潜熱は別に持たず、
  $L(T)=h_v(T)-h_l(T)$ をその差として作る (気相 $h_v$ は種 DB と**同じ評価・同じ外挿規約**)。これで `h2o_latent` の気相係数の二重ソース (120 K で 2.39 kJ/kg の差) が消える。
- **気液ペアの基準契約** (2 回目 M2): datum 不変性は「気液データの正しい組合せ」を保証しない (外部 DB の気相 H2O に定数 +100 kJ/kg が入っていても datum 不変試験は通り、L だけ +100 kJ/kg ずれる)。
  凝縮種の気相エントリと液相エントリは**同じ絶対基準 (CEA) のペア**であることを共通データで宣言し、整合を確認できない外部 DB による凝縮種 (気相) の上書きは凝縮 ON で拒否する。
  解決済み記録には液相係数・MW・延長規約・気相と共有するシフトを含める。
- **潜熱を使う全経路を移す** (2 回目 M3): `h2o_latent` だけでなく、`condFloat: 1` の潜熱表生成 (`condensationTables_d.cuh:106`)、表と double 退避の切替 (`condensationSourceF_d.cuh:283`)、
  `dL/dT` に依存する二相熱容量・音速 (`condensationEOS_d.cuh:184`)、Python の潜熱再ハードコード (`convert_species_field.py:72`, `:119` で湿り場の保存エネルギー再構築に使用)。
  新 L は現行比 120 K で −2387 J/kg、150 K で −466 J/kg (codex 検算) なので湿り場の変換に効く。
- codex plan レビュー M6 の反例 (sensible の気相から絶対基準の液相を引くと L(300 K) が 2.44 → 15.9 MJ/kg) は、液相に同じシフトを掛けない誤った置換の例であり、本方針では起きない。
  これを単体試験で固定する (§6 V7: datum を変えても L 不変)。
- 潜熱を独立の関数として持つ方式 (N2 の CPG 経路の現状: 液比熱が $c_l=c_{p,v}-dL/dT$ で暗黙に決まる) は H2O では採らない。
- **N2 も同じ方式に統一するかは未決** (2026-09-27 ユーザ: 「統一したい」→ 同日「変えなくてもいいのかも、今後判断」; §10 未確定事項)。統一する場合は、液相エンタルピーを気相と同じ datum で持ち、L をその差で作る。
  CEA の `N2(L)` は多項式を持たず **77.352 K の 1 点 (生成エンタルピー −12107 J/mol) だけ** (`thermo.inp:15618`) なので、液相は
  「CEA の 77.352 K 点を基準に、液の比熱モデル (現 `condN2LiquidCp` 2000 J/kg/K、63–77 K 実測 ≈2.0 kJ/kg/K) で積分」する形になる。
  現行の Lin 2014 の L フィット (`n2_latent`) と低温線形外挿 (`condN2LatentLowT`)、飽和圧の Clausius–Clapeyron 再構成 (`condN2PsatLowT`) との整合、
  CPG carrier 経路 (`thermalMethod 0`, 空気凝縮) での datum の扱いを決める必要がある。
- **湿り場の変換 (codex diagnose 2026-09-27, `notes/reviews/2026-09-27-h2o-latent-datum-diagnose.md`)**: `convert_species_field.py` の補正 `roe += ρ(e_gas,dst − e_gas,src)` (`:617-619`) は液相項を含まない。
  L が変わると `g(R_wT − L)` (`:124-127`) が変わるので、補正を `Δ(roe) = ρ{Δe_gas + g_dst(R_w,dst T − L_dst) − g_src(R_w,src T − L_src)}` にし、**液相モデルだけが違う**場合も再構成を発火させる。
  SRC/DST の二相 EOS はそれぞれの解決済み記録から作る。旧記録に液相情報が無ければ現在の resolver で補完せず、旧モデル (現行 `h2o_latent`) を明示指定する移行手順か拒否にする。
  潜熱の差 (150 K で ΔL −465.9 J/kg) は「統一分として記録するだけ」にしない — 場の変換の正しさは別に試験する。
- 本 plan の範囲では N2 の現行 L・飽和圧・CPG 経路を**保持**し、#10 の共通化 (ディスパッチ変更) で N2 の結果が変わらないことだけ確認する (2 回目 m6)。

### 4.9 熱物性を CEA thermo.inp から生成しソルバ内蔵にする (2026-09-30 ユーザ決定・diagnostician 諮問済み)

- **目的**: 熱物性 (NASA-9) も CEA thermo.inp から生成してソルバ内蔵にし、輸送物性 (trans.inp 66 種) と揃える。`legacy_builtin` (過渡欄) による内蔵の絞り込みは撤去する。
  実体は **61 種** (trans 66 種のうち thermo.inp に同名があるもの; `CCL2F2, CCL3F, CCLF3, CHCL2F, CHCLF2` の 5 種は thermo.inp に無い → 輸送のみ・熱力学不可として扱う)。区間数は 2 区間 39 種・3 区間 22 種、境界が標準 (200/1000/6000/20000) でないのは `e-` (Tlo 298.15 K) のみ → `e-` は内蔵から除外。
- **区間可変 (#6b) は schema 変更として扱う** (diagnostician 最重の指摘): 「2 区間」は device 構造体だけでなく、互換性テキスト (`speciesDB.cpp:872-875`, Python `forge_species.py:282-293`)、ハッシュ本文の schema・外挿規約文字列 (`speciesDB.cpp:1020-1022,1228,1240`; `forge_species.py:252`)、記録 reader (`speciesDB.cpp:976-995`, `forge_species.py:326-328,440-446`)、外部 DB (`speciesDB.cpp:322-325`)、後処理 `_TPGas` (`total_quantities.py:38-53`)、lump 合成の区切り検査 (`speciesDB.cpp:431-434`) に焼き込まれている。
  規約: **`nInt == 2` の種は直列化を一字一句現行と同じにし**、3 区間以上の種だけ拡張行を出す。schema・外挿規約の文字列は記録に `nInt ≠ 2` の種が 1 つでもあるときだけ別値 (#5t2-1 の `transport_compat` と同じ条件付き発行)。C++ と Python の鏡像 (`compat_text`・reader・`_TPGas`・`interp_field` 署名・外部 DB reader) は同じ段で直す。
- **既存 run との互換 (決定)**: 両立はしない。**全種を CEA そのもの (3 区間・CEA の MW・Ar 高温区間) に揃え、互換性ハッシュは最後の段 (段 3) で一斉に 1 回だけ変える**。7 種を 2 区間のまま据え置く案は採らない (Ar 高温区間の既知の誤りを恒久化し、1 つの表に 2 つの規約が混ざる)。過去 run の記録は run 側に残るので読める・照合できる。変わるのは新バイナリへの継続だけで、移行は #6a と同じ「明示許可 1 回」。
- 20000 K 区間はデータとして持つ (評価は各種の Thi でクランプ)。既存ケースは 6000 K を超えないので結果は不変。
- 名前: `parseBuiltinData` の `claim()` (`speciesDB.cpp:110-117`) を大小文字無視でも一意に強化 (thermo.inp には `CO`/`Co`・`CS2`/`Cs2` が実在; 66 種の範囲内では衝突なし)。CEA の `Air` は取り込まない (空気は lump で表す; 内蔵 `AIR` との衝突回避)。
- GPU: `SpeciesThermo`/`SpeciesThermoF` は device global にポインタ渡しなので構造体の拡大は律速でない。リスクは区間選択によるレジスタ増 → 段 1 G1-e で測る。

## 5. 実装ステップ

1. 内容照合 (§4.3) — 最初の一歩。
2. 共通データと生成器、C++/Python の読み込み (§4.1)。
3. 区間数可変の `SpeciesThermo` と起動時合成、config の lump 指定 (§4.2)。
4. lump の輸送物性展開 (§4.4)。
5. 設計 runner の切り替え (§4.6)、docs。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~plan 段 codex レビュー~~ | 完了 2026-09-27 (GO-with-changes, C0/M7/m1, 全件採用 → §2・§3・§4・#2–#11・§6 に反映; §6.1) | O |
| 2 | ~~仕様文書の先行更新~~ | 完了 2026-09-27: `methods/thermophysics.md` §1b (共通データ・canonical ID・lump の起動時合成・輸送展開と Blanc・解決済み記録と内容照合) | O |
| 3a | ~~記録・互換性ハッシュ・resolve-only・ソルバ入口~~ | 完了 2026-09-27: C++ (`input/speciesDB.{cpp,hpp}` の SHA-256・互換テキスト %.17g・記録・`diffRecord`・入口照合、`main.cpp` の `--resolve-species` と `valueFileName` 照合、`output/output.cpp` の res 属性 4 つ)、Python `forge_species.py` の内容比較 (builtin 省略・source 拒否を撤廃)。試験 `tests/unit/test_species_record_host.cpp` (29 項目)・`test_forge_species_record.py`・`test_species_record_solver.py` (17 項目) ALL PASS: V1 の A・B・(b)・(c)・旧場・source 差・resolve 一致。数値不変: `case/44.vitiated_air_wt/run_0516`–`0519` (基準/新バイナリ各 2 回, 200 step) `check_field_regress` PASS (比 ≤1.27)。**過渡期の既定**: 属性なしの場は警告 + `species_input_unverified=1` で通す (並行セッションの TP run の継続を止めないため)。`FORGE_REQUIRE_VERIFIED_SPECIES=1` で最終方針 (停止) を先取り。**#3b 完了時に既定を厳密へ切り替える** (#3c)。残る過渡の影響: `interp_field.py` は記録のない内蔵種を含む TP run 同士で unverifiable REFUSED (`--force-species` で通る) | O |
| 3c | ~~未検証の既定を厳密へ切り替え~~ | 完了 2026-09-27 (`feature/species-transport`): 属性なしの場は既定で停止、許可はその実行だけの `FORGE_ALLOW_UNVERIFIED_SPECIES=1`、`FORGE_REQUIRE_VERIFIED_SPECIES` をソルバから撤去、警告文の残課題 (#6a) も解消。試験 `test_species_record_solver.py` (19)・`test_species_record_host.cpp` PASS。Python ツールも 2026-09-27 に厳密化 (`forge_species.refuse_unverified`・`plan_inherit`・`plan_convert`・`stamp_new_field`・`interp_field` の署名照合・`runner_sern.warm_from_run`; 未検証 SRC・宛先を解決できない (旧バイナリ・config なし) は停止、許可は `FORGE_ALLOW_UNVERIFIED_SPECIES=1` / `--force-species` で属性なし)。`FORGE_REQUIRE_VERIFIED_SPECIES` を tools/tests/design から撤去。試験 entry (31)・interp_field_tracer・convert_fail・tpgas_lowT・design `run_species_attrs_ic_tests --no-cfd` PASS (`test_species_record_solver.py` は CFD を含むので未実行)。**残 (#3d)**: 印付きの場 (`species_input_unverified=1` かつハッシュ一致) をソルバは通して印を継承するがツールは停止する不一致 — 揃えるかは未決 (methods/thermophysics.md §1b.4 に現状記載) | O |
| 3b | ~~引き継ぎ・種変換・IC 生成の入口~~ | 完了 2026-09-27: 共通 API (`tools/forge_species.py` の `find_forge`/`resolve_species`/`plan_inherit`+`commit_inherit`/`check_ic_against_record`/`stamp_new_field`/`plan_convert`)、`restart_field.py`・`interp_field.py`・`convert_species_field.py`・`runner_sern.py` (`restart_by_index`・領域 IC・`warm_from_run`)・`design/forge_design/evaluate/ic.py` (`stamp_isentropic_ic_species`)・`runner_axismach.py` (IC 付与、段間を restart_field へ)。試験 `tests/unit/test_species_attrs_entry.py` (19 項目)・`design/tests/run_species_attrs_ic_tests.py` (7 項目) ALL PASS: V1 (a)(d)(f)、未検証 SRC で DST 属性を消す、`--force-species` は属性なし、種変換。既存試験も PASS。配管 run `case/44.vitiated_air_wt/run_0520_species_attrs_runner` (IC→段間 2 回→本段 60 step で unverified=0、restart_field はビット一致; 同じ場を interp_field で写すと roUx 22950/23725 点で最大 1.4e-5 相対の差)。解釈として決めた点: 宛先を解決できない (旧バイナリ等) ときは過渡期は警告で属性なし・strict で停止 / 属性はあるが記録なし・完全性不一致の SRC は過渡期でも停止 (`--force-species` で属性なし) / SRC 記録を DST の隣へ複製 / interp_field の内蔵種の照合不能は過渡期は警告 (#3a の残課題を解消)。未了: runner 経路 (段の手順が run_0509 と違う) の本段 24000 step の数値基準、`warm_from_run` の実 run、`runner_wt`/`runner.py`/case 内 `gen_*_ic.py` (CPG か範囲外) | O |
| 4 | ~~共通データ化 (値は変えない)~~ | 完了 2026-09-27: `solver_density_cuda/data/species/forge_species_v1.yaml` (13 種、canonical ID・別名・相・2 区間係数・LJ・元素組成・出典・`deviations`・過渡欄 `legacy_builtin`)、C++ はビルド時埋め込み (`cmake/embed_species_data.cmake`, `speciesDB_builtin()` のハードコード撤去)、Python `semiperfect.py`/`composition.py` は同ファイルを読む (名前は案 (a) のまま)。値のビット一致 `tests/unit/test_species_data_bitexact.py` ALL PASS (98 項目, 変異試験で FAIL を確認)。`--resolve-species` のハッシュは新旧で一致 (case/44 4378b7d78339ba27、内蔵のみ 3 構成) し記録ファイルもバイト一致。既存 species 試験は新バイナリで全 PASS。残: `tools/forge_species.py:69` の `BUILTIN_MW` の写し (#8)、He の atoms を新設 (読み手なし)、既存の `test_solver_config_species.cpp` は fixture の廃止キー `mesh.meshFormat` で FAIL 16 (本件と無関係、未修正) | O |
| 5 | CEA 直読みとの差・20000 K 区間・LJ の寄せ先 | H2O MW・AR 高温 a0・**He MW (0.0040026 vs 4.002602 g/mol, 相対 5e-7)** の寄せ先、**LJ の出典と寄せ先 (Python `LJ_PARAMS` と `cea_thermo_to_species_db.py` の Cantera 由来表が H2 2.827/59.7 vs 2.920/38.0、H 2.708/37.0 vs 2.050/145.0、O 3.050/106.7 vs 2.750/80.0、OH 3.147/79.8 vs 2.750/80.0、NO 3.492/116.7 vs 3.621/97.53、CO 3.690/91.7 vs 3.650/98.10 で食い違う)**、6000–20000 K 区間を有効にするか、CEA 全種を入れたとき現行 `AIR` の alias `Air` と CEA `Air` の衝突をどう解くか。値を変えるなら case/44 と #6 の小型ケースの報告量変化を記録 。**加えて (2026-09-27 ユーザ方針 §4.3b)**: 輸送物性の正本を CEA `trans.inp` にする設計 — 種別フィットの取り込み、相互作用データの混合則への入れ方、CEA 範囲外 (H2O < 373 K など) と CEA に無い種の扱い、拡散係数 (LJ 継続) との整合、NS 結果の変化量の記録 | F |
| 5t | 輸送物性の CEA 化 (段階) | codex diagnose の順: ① ~~**混合則の A/B (CFD 0 step)**~~ **完了 2026-09-27: B (CEA frozen 混合則) 採用** — 種別 μ・λ を CEA にそろえ混合則だけ変えて FCEA2 と 16 状態比較: A (forge 現行) 最大 μ 13.1 %・λ 10.7 % で不合格、B 最大 μ 0.001 %・λ 0.013 % で全点合格 (`notes/investigations/2026-09-27-cea-vs-forge-properties/mixing_ab.py`・`mixing_ab_result.txt`・`fcea_mix/`)。当初の設計: — 単成分フィット・MW・T・組成を固定し A = 現行 (φ 共用) / B = CEA の ηᵢⱼ・φ・ψ、T = 400/600/1000/2000 K × X_H2O = 0/0.1/0.5/1 の 16 状態を FCEA2 の μ・frozen λ と比較 (A 全点 0.1 % 以内なら「現行で不足」を棄却、A 失敗・B 全点合格なら B 採用、両方失敗なら GPU 実装へ進まない) → ② 共通 resolver と単成分評価 (ユーザ定義フィット → CEA → LJ、単位・区間端・記録) → ③ CEA 混合・lump 展開・CUDA 評価 (独立 double 基準に double ≤1e-12・float ≤1e-5) → ④ 低温 H2O モデル確定後の NS 検証 (収束・準定常 VERDICT 必須) | F |
| 5t2 | 輸送の実装 (codex diagnose 2026-09-27 の 3 段) | **段 1 CPU resolver・単成分・独立参照・記録**: 全 CEA の 16 状態を FCEA2 と μ・frozen λ ≤0.1 %、選択モデルは独立評価と ≤1e-12、H2O 500/700 K の左右極限で値の相対差 ≤1e-12・無次元勾配 T·d ln f/dT の差 ≤1e-10、未指定種の拒否漏れ 0 件、既存 TP 記録本文 0 byte 変化。前段の判別 A/B: 400 K 純 H2O で出所だけ `cea`/`custom:h2o_iapws_cea_v1` に変え、それぞれ CEA 値・IAPWS 値と ≤1e-12 (同値なら分岐が効いていない)。**段 2 混合・lump・GPU・壁**: 同じ入力で double ≤1e-12・float 格納値 ≤1e-5、重複実種を含む lump = full、種の列挙順を変えても同値、セルと壁で同一 T・組成なら同値。CPU 参照は実装側の混合関数を共有しない。**段 3 NS 統合**: 凝縮 plan の気相組成処理を固定してから輸送モデルの変更だけを評価。NaN/Inf 0、同一数値設定区間の `check_convergence` PASS、事前指定の報告量が `check_quasisteady --tail 0.4 --drift 0.01 --osc 0.01` STEADY。旧結果との不変は合格条件にしない。未決: case/16 の基準 run・判定区間・報告量、Warnatz・極性補正の独立参照 | F |
| 5t2-1 | ~~段 1 (CPU resolver・単成分・独立参照・記録)~~ | 完了 2026-09-27 (`feature/species-transport`): `data/species/forge_transport_v1.yaml` (trans.inp 66 種・相互作用 41 組、生成器 `tools/cea_trans_to_forge_transport.py`、ビルド時埋め込み)、`input/speciesTransportDB.{hpp,cpp}` (解決・検査・lump 展開・記録)、`cuda_forge/transportMix_d.cuh` (評価式; どのカーネルからも未使用)、`physProp.transport` の読み込み、記録 schema v2 (`transport_compat`)、H2O に双極子 1.844 D。試験 `tests/unit/test_species_transport.py` ALL PASS: 全 CEA 16 状態 vs FCEA2 μ 0.0012 %・λ 0.013 %、選択モデル 10 構成 × 21 温度で独立参照 (`transport_reference.py`) と最大 4.2e-15、判別 A/B (400 K 純 H2O: cea/custom が CEA 値/IAPWS 値と一致、差 μ +0.57 %・λ −2.26 %)、接続点 500/700/150 K で値 ≤2e-15・勾配 ≤2.5e-11、負例 18 件拒否、`physProp.transport` の無い config の記録・ハッシュはバイト不変 (run_0509 4378b7d78339ba27)。**`physProp.transport` を書いた config は計算の起動を拒否** (GPU は旧経路のままなので; 段 2 で外す)。H2O の冪外挿は 253.15 K からに変更済み (試験 (J) 253.15 K: 値 2.2e-16・μ 勾配 8.3e-11・λ 5.6e-12)。残: 両 kinetic の組に極性補正なし, `cea` で V・C の無い種は拒否 (UF6 のみ) | O |
| 5t2-2 | ~~段 2 GPU 接続~~ **(正しさは完了 2026-09-27、性能は未達)** (codex diagnose 2026-09-27, `notes/reviews/2026-09-27-transport-stage2-gpu-diagnose.md`) | **設計**: `physProp.transport` があるときだけ `gasProperties_d`・`wmlesWallModel_d` を新経路に (同じデバイス関数を共有)。**組成はモル基底で展開**: 正規化した輸送種 Y → 輸送種 X = (Y/M)/Σ(Y/M) → `Xreal_r = Σ_s X_s·expand[s,r]` (展開行列は lump 内モル分率, `speciesTransportDB.hpp:53`; Y に直接掛けると N2/He lump で X_He 0.6 → 0.887 になる — codex 検算)。重複実種は加算。ηᵢⱼ は各組を 1 回評価して両方向の分母へ (現 `transportMix_d.cuh:227` は i→j と j→i で再評価)。起動拒否 (`main.cpp:1189`) は接続試験の後で外す。**合格条件 (事前固定)**: 実際に GPU へ渡した入力を独立参照 (`transport_reference.py`) に渡し μ・λ ≤1e-12 相対・展開後 X ≤1e-12 絶対 (double 評価); float 格納値は独立 double 参照比 ≤1e-5 かつ参照を float へ丸めた値と ≤2 ULP; full/lump・列挙順は double 入力で ≤1e-12、別々に float 化した `roY` 経由は ≤1e-5 (各入力の独立参照との照合必須); セル・ghost (末尾まで)・壁が同じ基準を満たし未更新・NaN/Inf・非正 0 件; 範囲は単成分・重複 lump・実種 12 と上限 32・ゼロ分率、T = 200/253.15/400/500/600/700/1000/2000 K と各フィット境界 (境界両側は隣接 float); `physProp.transport` なしでは同一の固定入力に対し `vis_lam`・`thermCond` がビット一致。判別 A/B (T 400 K, ρ 1, 両実種 cea, 凝縮 OFF, 時間更新なし, セル・ghost・壁を各 1 回): A = full X_N2 0.4・X_He 0.6 / B = lump {N2 0.5, He 0.5} + 独立 He, X_L 0.8・X_He 0.2。性能 (提案値): 固定 GPU・入力で物性更新の増分が旧 step 時間の 10 % 以内。case/44 Euler 起動や凝縮 200 step は合格根拠にしない (段 3)。**μ・λ の利用箇所 (接続確認表; 今回すべて改修はしない)**: 粘性・熱流束 `viscousFlux_d.cu:164,265`、SST・遷移・スカラー拡散 `ransSource_d.cu:97`・`transition_d.cu:149`・`scalarTransport_d.cu:120`、定数 Sc の化学種・受動種拡散は新 μ で変わる `speciesTransport_d.cu:259`・`passiveKernels_d.cuh:171`、軸対称・CHT `axisymmetricSource_d.cu:239`・`conjugateWall.cpp:242`、SST 壁関数は局所 Pr が更新されるが回復係数は `cfg.prandtlLam` のまま `ransWallFunction_d.cu:347,583`、CFL・陰解法は定数 `physProp.visc` `setDT_d.cu:82`・`timeIntegration_d.cu:813,1496` (既存の近似として明記) | O |
| 5t2-2r | 段 2 の結果と残り | 実装: `transportMix_d.cuh` に `TransportTemp`・`TransportTableD`・`transport_expand_X`・`transport_mix_Y` (組成はモル基底で展開)、ηᵢⱼ は各組 1 回、`thermo_init_db` が device 表を上げる、`gasProperties_d`・`wmlesWallModel_d` が同じ関数を使う、物性だけ評価する 0 step ハーネス (`FORGE_TRANSPORT_PROBE`) と `tests/unit/test_transport_gpu.py` (43 項目 ALL PASS: double 3.5e-15・float 5.5e-8 かつ 0 ULP・A/B で X_He 0.600000000000000・full/lump/列挙順・セル/ghost/壁・既定経路ビット一致 4 構成; 基底取り違えの変異版で 14 FAIL)。**性能は提案値 10 % を超過** (RTX 3060, 23725 CV): 物性更新が旧 0.06 ms → 実種 5 で 0.41 ms (+14 %)・12 で 2.57 ms (+94 %)・32 で 19.4 ms (+800 %)。FP64 の exp/log を O(n²) で毎 step 評価するため。対策は未決 (§10)。`physProp.transport` と `viscMethod ≠ 2` の併用は起動拒否 (記録と計算の食い違いを防ぐ; implementer の判断、要確認)。μ・λ の読み手に `turbulent_viscosity_d.cu:201,359`・`ransBoundary_d.cu:43` も追加。CFL (`setDT_d.cu:82`)・陰解法対角 (`timeIntegration_d.cu:813,1475,1496,1567`) は定数 `physProp.visc` のまま (既存の近似) | O |
| 5t2-3 | ~~輸送の表引き化~~ **(完了 2026-09-27)** (codex diagnose 2026-09-27, `notes/reviews/2026-09-27-transport-tables-diagnose.md`) | **設計**: 種別 ln μᵢ・ln λᵢ と組 ln ηᵢⱼ (CE・CEA 相互作用; 剛体球は種別表の μ から実行時) を、**式の区間境界・H2O の接続点・T* クランプ点・NASA Tmid で表を分割**し、各分割区間の中で ln T 等間隔 (Δln T ≤ 1/256、約 1179 区間/物性) の 3 次 Hermite (両端の片側値・片側微分から係数; 右区間の左端は右側の式を明示して評価)。**区間選択は元の T と元の境界値で行い (現行 `T <= Thi` と同じ所属; `transportMix_d.cuh:118`)、float の ln T では選ばない** (999.99994/1000/1000.00006 K の ln T は float で同値になり、1 % 段差の fit で 9.9e-3 の誤りになる — codex 算術反例)。範囲 150–15000 K、範囲外は現行の double 評価へ委譲 (端値クランプしない)。表引き・混合とも float。刻みは初期値で、精度不合格なら該当区間を細分。**合格条件 (事前固定)**: 実際の float 入力を独立 double 参照へ渡し、単体 μᵢ・λᵢ・ηᵢⱼ ≤2e-6 相対、混合後の GPU 格納値 ≤1e-5; 検査点は各小区間の 17 等分点・全境界の直前直後 8 個の float・表範囲外; NaN/Inf・非正・未更新 0 件; 既定経路ビット一致、full/lump・列挙順・セル/ghost/壁を維持 (新 float 経路には旧 double の ≤1e-12 を課さない); 全 CEA 指定の FCEA2 比 ≤0.1 % を残す; 判別 A/B (区間選択だけ変える: A = float ln T で選択 / B = 元の T で選択; 1000 K 以下 μ 1e-5・超過 1.01e-5 Pa·s の単成分 fit を 1000 K と両隣の float で 0 step 評価; A 不合格・B 合格が期待)。**性能**: n = 5 (va3 lump 記法 MIXDRY + H2O を実種 5 に展開) と n = 12 (段 2 試験と同じ構成) で `(新物性時間 − 旧物性時間)/旧 step 時間 ≤ 0.10`、同じ入力・GPU・ビルド条件で暖機後 300 回 × 5 反復、反復ごとの比の最大値で判定、32 は記録のみ。不合格でも精度・微量種を削らず #5t2-3 を未達とする。境界付近の誤差を刻みだけで直す・元の段差を連続化する・微量種を切るはしない | O |
| 5t2-3r | 表引き化の結果 | 実装: `cuda_forge/transportTables_d.cuh` (分割区間ごとの ln T 3 次 Hermite、区間は元の T の閾値で選択、範囲外と NaN は段 2 の double へ委譲)、`thermo_d.cu` で構築・アップロード (失敗なら起動停止)、セル・壁が `transport_mix_Y_tab` を共有、`FORGE_TRANSPORT_TABLE=0` で double 評価に戻せる。試験 `tests/unit/test_transport_gpu.py` ALL PASS: 単体 ≤4.7e-7 (約 1230 万点、境界 ±8 float・範囲外含む)、混合 ≤4.1e-7、FCEA2 μ 0.0012 %・λ 0.013 %、判別 A/B は A 9.9e-3 不合格・B 4e-8 合格、既定経路ビット一致 4 構成、full/lump・列挙順・セル/ghost/壁。性能 (RTX 3060, 23725 CV, 暖機後 300 回 × 5): 物性時間 n5 0.049 ms (旧 Wilke 0.062)・n12 0.20 (旧 0.31)・n32 1.34 (旧 0.44) → (新−旧)/旧 step の最大 n5 −0.006・n12 −0.028 で合格、n32 +0.21 (記録のみ)。表のメモリ n12 0.76 MB。注記: n12/n32 の計測場は run_0509 の res_24000 に固定組成の ρY を載せたもの (3 経路に同じ入力); 150 K 未満・15000 K 超のセルは double 評価に戻る (その分遅い)。T2 格子点検査は n12・n32 で間引き (単体 T1 は全小区間) | O |
| 9b | ~~runner の輸送指定~~ **(axismach 完了 2026-09-27、SERN 未)** (2026-09-27 ユーザ指示で段 3 より先) | 観測: 設計 runner の NS/SST は `viscMethod: 1` (空気の Sutherland; `runner_wt._config_sst_node`・`runner_sern.py:246`・`runner.py:103`) で、TP 混合物でも μ は空気の値。**やること**: problem YAML `gas.transport: {実種: モデル}` (lump 構成種を含め実種ごとに必須、書いたときだけ) を検査し、NS/SST 経路の config を `viscMethod: 2` + `physProp.transport` にする。Euler 経路には書かない (viscMethod≠2 との併用は拒否)。**改訂 (2026-09-27 ユーザ決定): semiperfect TP の NS/SST は `viscMethod: 2` + `physProp.transport` を既定とし、`gas.transport` が無ければ prepare でエラー** (CPG と Euler 経路は従来どおり・バイト不変; 既存の semiperfect NS 問題 YAML は `gas.transport` の追記が要る)。合格: `--prepare-only` で config に全実種の指定が入り `forge --resolve-species` が通る、指定漏れ・誤用は prepare で拒否、CPG・Euler の生成 config がバイト不変、NS の短い run (50 step 程度) が NaN なく起動。**未決 (§10)**: NS の既定を Sutherland から種ごとの物性へ切り替えるか (δ* 生産レシピ case/42・44・45 の結果が変わる) | O |
| 9b-r | runner の輸送指定の結果 | 実装: `gas.transport` の検査 (`design/forge_design/probdef.py`・`gas/composition.py` の `parse_gas_transport`/`resolve_transport`; 実種ごと必須、lump 名・余分な種・未知モデル・custom の誤用・大小文字違いの重複を拒否)、`runner_axismach.prepare_ns` は TP のとき `viscMethod: 2` + `physProp.transport` (`thermCondMethod` は落とす; `visc`・`thermCond`・`prandtlLam` は残す) を書き、`gas.transport` が無ければ run dir を作る前にエラー (必要な実種と書き方の例を表示)。`species_meta.yaml`・`prepare_info.json` に記録。試験 `design/tests/run_transport_spec_tests.py` (正例 11・負例 13) ALL PASS。CPG・`cfd_gas: cpg`・Euler 経路と runner.py の生成物は変更前とバイト一致、case 配下 343 件の problem YAML で load 失敗 0。case/42 split_ns のコピーで prepare_ns → `--resolve-species` 通過 → 50 step (soft 段相当, 表引き経路) NaN なし (scratch、削除済み)。**残**: SERN (`runner_sern`, frozen_tp) は lump 記法への切り替え (§5.2) 待ちで未対応 (frozen_tp で `gas.transport` を書くと拒否、NS は Sutherland のまま); `runner.py` (thruster_bell)・`runner_wt.prepare(euler=False)` は常に CPG なので対象外; 既存の semiperfect NS 問題 YAML は 2026-09-27 に追記済み (case/42 6 件・44 7 件・45 9 件、N2/O2/AR/CO2 = cea、H2O = custom:h2o_iapws_cea_v1; 全 343 件の検査で追記漏れ 0); `prepare_ns` の近壁 ω 床の IC は Sutherland と CPG の T で ν を見積もる (設計ツール plan `tooling-design-problem-campaign-recipe.md` #3 の熱力学修正と同じ箇所) | O |
| 6a | ~~起動時 lump 合成 + config 指定 (区切りが揃う種のみ)~~ | 完了 2026-09-27: `physProp.species` の mapping 記法 (`input/solverConfig.cpp`, `input/speciesLump.hpp`)、`speciesDB_resolve` での合成・検査 7 種・起動ログ・記録と互換性ハッシュ (lump なしはバイト不変)、Python `forge_species.py` の lump 対応。試験 `tests/unit/test_species_lump_solver.py` ALL PASS (31): **V2 = 生成 DB と係数・MW・LJ 相対 ≤4e-16、cp/h/s° (200–6000 K 1000 点) 相対 ≤8e-16**、既存 config のハッシュ・記録バイト一致、負例すべて拒否。起動確認 `case/44.vitiated_air_wt/run_0521_species_lump_startup` (lump 記法で 200 step、NaN 0)。basis は必須。**lump 記法と外部 DB 版はハッシュが違うので既存場からの restart は照合で止まる — 正しい挙動として採用** (1 ulp の差があり、外部 DB 擬似種の出自は場から確かめられない; 移行は明示許可 1 回)。残: lump の mapping を読めない Python reader (`total_quantities.py`・`convert_species_field.py`・`passive_gate_common.py`・`gen_inlet_profile.py`・`runner_sern._species_signature`・design `probdef`) は #8/#9。過渡期既定で通すときの警告文が環境変数を立てたように読める (#3c で直す) | O |
| 6b | ~~区間可変 (区切りの違う種を畳む)~~ → **#13 段 1 に統合** (2026-09-30) | §4.9。schema 変更として扱う | — |
| 13-0 | CEA 熱物性の内蔵化 段 0: データ監査 (コード変更なし) | §4.9。生成器の dry run で表を出す: 61 種の区間数・境界、無い 5 種、`e-` の扱い、既存 7 種の現行値と thermo.inp パース値の差分 (係数ビット一致・MW・Ar 高温)、6 種 (CO/H2/OH/H/NO/O) と SERN 外部 DB の係数差分。合格: 表が揃い plan に貼る。判断: 2026-09-30 diagnostician — 段分けと各段の合格条件を採用 | O |
| 13-1 | 段 1: 区間可変の構造のみ (データ不変) | §4.9。**G1-a** 既存記録ハッシュ不変 (`run_0509` の compat_hash 4378b7d78339ba27、`test_species_data_bitexact.py` 98 項目、`test_species_record_*`・`test_forge_species_record.py`)。**G1-b** 7 種 × T = {200, 999.99994, 1000, 1000.00006, 6000, 7000} K で cp/h/s°/T_from_e (double・float) が新旧 0 ulp。**G1-c** V3 (相対 ≤1e-12、区切り ±1e-9 K、段差増分 ≤ V2 許容; 試験種は 1000/6000/20000 の 3 区間型と Tlo 298.15 型)。**G1-d** V3f (`tools/test_thermo_float.cpp` の `errHyb/T < 3e-8` を 20000 K まで)、上限超過の起動拒否 1 件。**G1-e** SLAU の REG が 136 を超えない・起動失敗 0、case/44 run_0509 config 200 step で `check_field_regress` PASS (比 ≤1.27)。**G1-f** 3 区間種を含む記録の書き→読み→再ハッシュ一致、2 区間のみの記録は schema 文字列不変。`cuda_forge/thermo_d` の変更 (エスカレーション 6 は本諮問で済み) | O |
| 13-2 | 段 2: 生成器 + 内蔵拡大 (既存 7 種 + H2O(L) の値は不変) | §4.9。生成器 (`tools/cea_thermo_to_forge_species.py` 等) 往復で 61 種の全区間・全係数・MW が thermo.inp パース値とビット一致。既存 8 エントリの YAML ブロック diff 0 行 → G1-a 再実行で不変。6 種は段 0 の差分表どおり (一致なら SERN run の `--resolve-species` ハッシュ不変を 1 run で確認)。`claim()` を大小文字無視で一意に (負例 `CO`/`Co` で起動拒否)。Python `semiperfect.py`/`composition.py` の値ビット不変 (bitexact 試験拡張)。LJ の出典を 1 つに決めて `source` に (7 種の LJ 変更は段 3) | O |
| 13-3 | 段 3: 既存 7 種を CEA そのものへ (**ハッシュが一斉に変わる唯一の段**) | §4.9。先に Δ 表 (種ごとの 200–6000 K max \|Δcp\|/cp・\|Δh\| J/kg・\|Δs°\|) を出し予測と照合 (N2/O2/CO2 は往復一致なら 0、H2O は MW 由来 ≤1.2e-6 相対、He ≤5e-7、Ar は高温区間差で h 数 J/kg を実測) — 予測超過なら止まる。case/44 V5 (i) 再実行で V0 許容 (ṁ 1e-4 相対・M 1e-5・T 0.01 K・軸 M 目標差 1e-5)、報告量 STEADY。旧記録からの継続は係数不一致キーを示して停止し明示許可でのみ通る (V1(c) 型、runner の warm start も同じ)。影響範囲 (全 TP run) と旧→新ハッシュ対応を変更ログに | O (結果の解釈は F) |
| 7 | lump の輸送物性展開と Blanc 拡散 | §4.4 (粘性・熱伝導は実種展開、lump を含む二元拡散係数は Blanc)。合格は §6 V4 | O |
| 8a | ~~reader の移行 (記録から読む)~~ | 完了 2026-09-27: 共通読み出し `forge_species.run_thermo` (res 属性の記録 → run の記録 → `speciesDBFile` → `--resolve-species`)、`total_quantities.py`・`convert_species_field.py`・`gen_inlet_profile.py`・`runner_sern._species_signature` を移行。V6 PASS (下の #9) | O |
| 8 | Python 共通 API の残り、canonical ID への移行 | 残り: `forge_species.species_info` の lump MW (`BUILTIN_MW` の写し)・`species_signature` (lump を照合不能扱い)・`interp_field` の署名。 §4.6。加えて (2026-09-27): Python の種名の大文字化をやめ canonical ID + alias 表へ、C++ `ResolvedSpeciesDB::index()` (`speciesDB.cpp:79-85`, 大小文字無視; 重複検査 `:210` も使う) の完全一致化と、tracer・凝縮種など名前で引く箇所の影響調査。**互換性ハッシュには config の名前がそのまま入る (`speciesDB.cpp:461`) ので、canonical 化で既存記録と不一致にならない規約 (ハッシュには canonical ID を入れ、既存記録は移行ツールで読み替え等) を設計してから**。合格は §6 V6 | O |
| 9 | ~~設計 runner の切り替え~~ (axismach 完了、SERN 未) | 完了 2026-09-27 (axismach): `_apply_gas_to_config` は `physProp.species` を lump 記法 (全桁の正規化モル分率) で書き **`species_db.yaml` を作らない**。内蔵に無い種・外部 DB (`gas.species_db`) が内蔵値を上書きする種だけ生エントリを `species_db_external.yaml` に置く (合成物は置かない)。変換器は `FORGE_BIN` と同じビルドのもの (`runner.converter_path()`)。**V5 (i)** `case/44.vitiated_air_wt/run_0522_species_nodb_lumpX_v5` (ref): V0 (run_0509/0513–0515) との最大差 ṁ_in 6.1e-7 相対・ṁ_out 2.6e-7・出口 M 9.5e-7・出口 T 1.2e-4 K・軸 M 出口 **8.6e-6 (許容 1e-5, 反復差の約 4.5 倍; 原因未切り分け)**・軸 M 目標差 3.6e-7 → 許容内。本段区間 NOT CONVERGED (plateau, run_0509 と同型)・series ALL STEADY、NaN 0、メッシュ PASS。**V6**: DB ファイルなしで prepare → 段間 restart_field (記録継承) → 本段 → `total_quantities` (旧経路比 T0 1.2e-15・P0 2.6e-14) → lump→full5 変換 `run_0523_species_nodb_full5_convert` (ρY 保存差 0、変換後 200 step で照合一致・NaN 0) → `gen_inlet_profile` (CSV バイト一致)。**SERN は未切替**: SERN の lump 名 `AIR` がソルバ内蔵の擬似種 `AIR` と衝突し起動時に拒否される (`speciesDB.cpp:329-334`)。lump 名の変更か内蔵 `AIR` の扱い (#5 の Air 衝突と同根) を決めてから | O |
| 10 | 潜熱: 液相を気相と同じ datum でシフトし差で L を作る | §4.8。`h2o_latent` の H2O 気相再ハードコードを撤去し、液相 H2O(L) (共通データの凝縮相エントリ) に**気相 H2O と同じ datum 定数**を適用して `L = h_v − h_l` を作る (気液差は CEA のまま保たれる)。273.15 K 未満/373.15 K 超の液の延長規約は現行のまま。**移行対象**: 潜熱表生成・範囲外退避・二相熱容量/音速 (`dL/dT`)・面流束・二相反転・Python `convert_species_field.py` の潜熱 (§4.8)。外部 DB の気液ペア契約と拒否。合格は §6 V7・V1 (液相だけ変えた restart の拒否)。判断: 2026-09-27 codex diagnose (`notes/reviews/2026-09-27-h2o-latent-datum-diagnose.md`) — 同 datum・同 MW の設計は採用、Python 変換式に液相項の差を入れる (§4.8, V7(e′))、表分割は判定で決める (V7(c′))。実装 2026-09-28: 共通データに `H2O(L)` (pair_of H2O)、`h2o_latent` の気相再ハードコード撤去 (種 DB の気相と液相の差、液相 datum は定数加算、kernel へはポインタ)、外部 DB の気相 H2O は内蔵とビット一致のときだけ許可、記録・互換性ハッシュに `condensed[0]` (凝縮 OFF はバイト不変)、`convert_species_field.py` の液相項補正と `--src-latent legacy-v0` 移行。V7 (a) 2.1e-15・(b) 7.6e-16・(c′) L 7.9e-8/ΔL′ 比 4.5e-3・(d) T 2.9e-10 K・(e′) A −0.0058022 K / B +4.6592006 J/kg (誤差 3.2e-10) PASS、N2 ビット一致。ΔL (新−旧) 150 K −465.92・120 K −2387.19 J/kg、200 K 以上 3e-15。**残**: V7(f) 湿潤回帰 (基準 run 未固定, AWS)、ソルバ起動経路 (`[cond] H2O latent heat` ログ・device 検査) の実行確認、`test_species_record_solver.py` (1 step)。既存 H2O 凝縮 TP 場の restart はハッシュ変更で拒否 → `convert_species_field.py --src-latent legacy-v0` で移行 | O |
| 12 | N2 の潜熱を同方式に統一するかの判断 (未決) | §4.8 末尾・§10。H2O (#10) の後に、統一するか現行 (Lin フィット + 低温外挿) のままにするかをユーザと決める。統一する場合: 液相 = CEA `N2(L)` 77.352 K 点 + 液比熱モデル、L は差。飽和圧の低温再構成・CPG carrier 経路との整合を決めてから。合格: V7 と同じ datum 不変試験 + 現行 Lin フィットとの差を 45–120 K で記録 + 空気凝縮 run (**case/34 Arthur**; plan condensation-air の検証先。旧記載の case/28 は He/空気同軸ジェットで誤り) の onset 変化を記録。**H2O (#10) 完了の阻害条件にしない**。凝縮カーネル変更なので編集前に上位へ諮る | F |
| 11 | docs 同期 (完了時) | `procedures/solver-settings.md` (`physProp.species` の lump 形、`speciesDBFile` の位置づけ、lump と拡散の制約)、`recommended-settings.md` §3、`design/CAPABILITIES.md` | O |

### 5.2 SERN への連絡事項 (2026-09-27 起票, 2026-09-30 改訂: 案 (a) で決定)

SERN セッション (`feature/sern-design`) へ渡す内容。変更は `feature/species-transport` (HEAD 8077e509 時点)。**2026-09-30 ユーザ決定: 名前の衝突は (a) SERN の lump 名を変える** (内蔵 `AIR` は既存 config が使うので残す)。

1. **lump 名の変更 (決定)**: SERN の外気 lump `AIR` を **`AMB`** (ambient; 外気の物理的な意味。内蔵種・別名と衝突しないことを確認済み 2026-09-30) に改名する。内蔵擬似種 `AIR` (cp/R 3.5 一定) と衝突し、lump 記法では起動時に拒否されるため (`speciesDB.cpp:329-334` 付近)。
   対象: `runner_sern.py` の `SPECIES_ORDER`・`tp_species` 既定 (`:88,110`)・`FrozenGas.from_mole(..., "AIR", ...)` (`:108`)・段間継承のコメント/式 (`:729-732`)、meta、後処理 (`plot_sern*.py` 等の種名参照)、problem YAML に `AIR` と書いた箇所。順序 (`EXH`, `AMB`) は変えない。
2. **何が変わったか**: ソルバが config の lump 記法 `physProp.species: [{name: X, lump: {構成種: モル分率}, basis: mole}, ...]` から NASA-9 を起動時に合成する。合成済み擬似種の `species_db.yaml` は不要。~~EXH の構成種 CO/H2/OH/H/NO/O も内蔵データにある~~ → **訂正 (2026-09-30, SERN セッション指摘)**: 6 種は共通データにあるが `legacy_builtin: [design]` のみで、ソルバの内蔵 (`speciesDB.cpp:118`, `legacy_builtin: solver` だけ読む) には入らない。生の係数を `species_db_external.yaml` で渡す (axis-Mach ランナーと同じ形; SERN は `5eee3ff2` でこの形に実装済み)。当初の誤りは輸送データ (`forge_transport_v1.yaml`) だけを確認して熱物性側を見なかったため。
   ソルバは `resolved_species_<hash>.yaml` を書き、res に `species_hash` を付け、起動時に入力場と照合する。**未検証 (属性なし) の場は既定で停止** (ソルバ・Python ツールとも)。許可はその実行だけ `FORGE_ALLOW_UNVERIFIED_SPECIES=1` か `--force-species` (許可して書いた場は属性なし)。`FORGE_REQUIRE_VERIFIED_SPECIES` は撤去済み。
   `forge --resolve-species` で GPU なしに宛先ハッシュ。後処理は `forge_species.run_thermo` で記録から読む。
3. **輸送物性 (必須)**: `viscMethod: 2` は `physProp.transport` (実種ごとの μ・λ の出所) が必須になり、無ければ起動時エラー。混合則は CEA frozen。設計 runner は semi-perfect TP の NS/SST で `viscMethod: 2` を既定にし、problem YAML に **`gas.transport`** が要る。
   SERN の実種 (N2, O2, AR, CO2, H2O, CO, H2, OH, H, NO, O) は、H2O を `custom:h2o_iapws_cea_v1`、他を `cea` にする (いずれも CEA trans.inp のデータあり)。書き方は case/44 の problem YAML を参照。現行 `viscMethod: 1` (空気 Sutherland) からは **NS/SST の結果が変わる**。
4. **切り替え作業**: `runner_sern` の config 生成を lump 記法へ (`runner_axismach._apply_gas_to_config`、`composition.solver_species_config`/`physprop_species_flow` が参考)。`write_species_db` (`:528`) の擬似種 NASA-9 書き出しはやめる。
5. **バイナリ**: lump 記法・`physProp.transport` の config は旧いソルバ・変換器で読めない。`FORGE_BIN` で新ビルドを指定 (TP の準備で必須; runner は同じビルドの変換器を使う)。AWS の clone も取り込み・再ビルドが要る。
6. **既存 run からの継続**: 改名と lump 記法で互換性ハッシュが変わるので、旧場からの restart は照合で止まる。種の順序が同じ (EXH, 外気) なので、移行は `restart_field.py --force-species` (またはその実行だけ `FORGE_ALLOW_UNVERIFIED_SPECIES=1`) で**1 回だけ**行い、以後は新しい記録で継承する。
7. **回帰**: 切り替え前後の比較は 2 段に分ける。(i) 改名 + lump 記法のみ (輸送は旧のまま比較できる構成で) で、代表作動点の推力・モーメントが事前に決めたノイズ床以内。(ii) 輸送物性の切り替え (`viscMethod: 2` + transport) は結果が変わる前提で、変化量を記録する (不変を合格にしない)。
8. **取り込み方**: 種 DB・輸送関係の commit だけを `feature/sern-design` に取り込むか main 経由で合わせるかは SERN 側の都合で決める。凝縮関係 (潜熱 #10・気相組成の輸送) は SERN が凝縮 OFF なら影響しない。

### 5.3 gap-heating セッション (`feature/gap-heating-precision`) との合流準備 (2026-09-27)

- **git 上の衝突なし**: merge-base `d8c11f52` 以降、gap 側 3 commit (case/60・notes・plan) と本ブランチ 20 commit に共通ファイルなし (`git merge-tree` で無衝突を確認)。
- **実行時の非互換 (合流前に gap 側で対応が要る)**: case/50・55・56 の run config は `viscMethod: 2` + `species: [N2, O2, CO2, H2O, AR]` + `speciesDBFile` で **`physProp.transport` が無い** → 本ブランチのバイナリでは起動時エラー (§4.3c 案 C)。
  追記する指定: `physProp.transport: {N2: cea, O2: cea, CO2: cea, AR: cea, H2O: "custom:h2o_iapws_cea_v1"}`。
- **結果が変わる**: 混合則が Wilke 共用 φ → CEA frozen (λ は ψ) になり、種別 μ・λ も CEA/IAPWS に変わる。gap 側の熱流束・Eckert 比の既報値は旧物性に基づくので、合流後の run は「物性変更による差」を分けて記録する (旧 run の再現が要るなら合流前のバイナリを残す)。変化量はその mixture・温度域で gap 側が合流時に見積もる。
- **合流の時期**: gap 側の区切り (ユーザ判断) を待つ。本ブランチから gap 側の config は書き換えない (並行セッションの run を触らない)。

## 6. 検証

事前に決める合格条件 (結果を見てから変えない)。**数値を変える作業 (#6 以降) の前に**、比較の基準 (V0) を固定する:
[`tooling-design-problem-campaign-recipe.md`](tooling-design-problem-campaign-recipe.md) §5.1 #2・#4 の抽出関数・許容差と、**旧バイナリ** (本 plan 着手前の commit で build し sha256 を記録) を固定してから V5 を回す。
この依存が満たされるまで #6 以降の CFD 比較をしない。

**V0 取得済み (2026-09-27)**: 基準バイナリ `~/forge-ref-bin/species-db-baseline/forge` (sha256 45cf9822…, ソース b04e7de9 と数値的に同一、2026-09-24 01:45 build)、抽出は `case/44.vitiated_air_wt/lumpX_series_csv.py` (版 = commit b04e7de9 時点) の最終スナップショット。`run_0513`–`0515` (run_0509 と同一入力・同一 prepare、md5 一致メッシュ) の 3 反復: ṁ_in 幅 1.1e-6 相対、出口 M 4.8e-7、出口 T 1.2e-4 K、軸 M 出口 1.9e-6、軸 M 目標差 max 6.0e-7 (run_0509 との最大差も同程度)。**許容差 = 反復差 × 3 と下限の大きい方 = ṁ 1e-4 相対・M 1e-5・T 0.01 K・軸 M 目標差 1e-5** (すべて下限が効く)。3 本とも本段区間 NOT CONVERGED (run_0509 と同じ plateau)。

- **V1 (保存場との照合, CFD 0 step)**: codex 諮問の判別 A/B (作成時記録を固定、現在の内蔵 `N2.nasa9_low[2]` を A = 同値 / B = +0.001) に加え、
  (a) 保存場だけを別 run にコピー、(b) 記録ファイルの取り違え、(c) 外部 DB の係数変更、(d) 起動前の宛先解決、を試す。
  合格: A と (d) の一致ケースは許可、B・(a) の不一致・(b)・(c) は**該当係数を示して拒否**、記録の無い旧場は「照合不能」で止まり明示フラグでのみ通る。これを
  ソルバ直接読込・`restart_field.py`・`interp_field.py`・`convert_species_field.py`・runner の段間引き継ぎの**全入口**で確認。
  (e) 凝縮種の**液相エントリだけ**を変えた場合も拒否されること (2 回目 M2; #10 完了まで未達)。
  期待エラーは「照合不能 (属性・記録なし)」と「係数不一致 (差のある種・キーを示す)」を分ける (記録が失われた場は差を特定できない)。
  (f) IC 生成の datum の結合: 宛先 `thermoHrefTemp=298.15` 固定で、IC 生成の `h_ref_T` だけ A=298.15 / B=0 → A 許可・B 拒否 (codex 判別 A/B)。
- **V2 (起動時合成 = 現行生成 DB, double)**: case/44 va3 の lump をソルバで合成した係数が `run_0510` の生成 `species_db.yaml` と相対 1e-12 以内、
  200–6000 K の 1000 点で cp/h/s° が相対 1e-12 かつ絶対 cp 1e-9 J/(kg·K)・h 1e-6 J/kg・s° 1e-9 J/(kg·K) 以内。
- **V3 (区切りの違う種を畳む, double)**: Tmid ≠ 1000 K の試験種を含む lump を合成し、構成種ごとの重み付き和と 100–20000 K の全点 (**各区切り温度そのものとその両側 ±1e-9 K** を含む、外挿域を含む) で V2 と同じ許容差。
  区切り温度での h の段差は、**現行係数が持つ段差を基準として保持**し (1000 K で N2 1.38e-4、H2O 1.90e-2、CO2 2.72e-3 J/kg; codex 検算)、datum 適用・合成による**段差の増分**が V2 の許容差以内であることを確認する (連続化のための係数修正は本 plan に混ぜない; 2 回目 M5)。
- **V3f (float 経路)**: float 表 (`SpeciesThermoF`) の cp/h と e↔T 往復が、既存 `tools/test_thermo_float.cpp` の基準 (`errHyb/T < 3e-8` 等) を区間可変後も満たす。上限超過が起動時に拒否されること。
- **V4 (輸送物性と拡散)**: lump を含む混合の粘性・熱伝導が、全実種で直接評価した Wilke / Mason–Saxena と機械精度で一致。旧方式 (平均 LJ) との差を 200/300/1000 K で記録。
  lump を含む二元拡散係数が Blanc の式 (実種の $D_{ij}$ から) と機械精度で一致し、旧方式 (平均 LJ) との差を記録。補正後の `J_i*` と `Σh_i J_i*` を独立な参照計算 (`full` 展開の double 実装) で検査し、`full` との差を縮約近似の誤差として記録 (合否は「lump 内組成が一様な場で差 ≤1e-12」、組成勾配のある場は記録のみ)。
  lump + kinetic 拡散で**組成勾配を持つ試験** (2 流入の混合層など小型ケース) の化学種・エネルギー収支は、**収支式・正規化・許容差・実行コマンドを実装前にここへ書いてから**回す (node の開境界のピンによる交換を含める; `check_passive_budget.py` は受動種用で代用しない; 2 回目 M4)。
- **V5 (CFD 回帰)**:
  - (i) case/44 va3 M4.19 L_c8 dry を新 config (lump 指定、`species_db.yaml` なし) で回し、固定済み V0 の抽出関数・許容差で旧バイナリの同条件 run と比較。
    これは**準定常回帰** (残差は plateau で NOT CONVERGED のまま; 生の VERDICT を記録し、報告量が `check_quasisteady` STEADY であることを条件にする)。
  - (ii) 全残差 `check_convergence` **PASS** が得られる小型 node TP ケース (前身 plan で使った case/16 の 5 種 run 系列など、着手時に 1 つ選んで §6 に追記) で、新旧の場と報告量を比較。
  - (iii) NS: 対象ケース・輸送設定 (lump + Schmidt 拡散)・保存収支・判定量を**着手前にここへ書いてから**回す (#7 の輸送変更による変化量を記録し、不変は合格条件にしない)。
  - いずれもメッシュ品質 VERDICT、IC と datum の整合 (step 0 で T が跳ばない)、段階起動の段と判定区間 (`stage_manifest.json`) を明記。
- **V6 (DB ファイルなしの一貫経路)**: 新 config で prepare → 段階 restart → `total_quantities.py` → `convert_species_field.py` までを、生成 `species_db.yaml` なしで通す。
  全温・全圧が旧経路と V2 の許容差で一致。SERN runner の段間署名も同じ API で通る。

- **V7 (潜熱)** — double 試験の許容差 1e-12 は double に限る:
  - (a) datum 不変性: 下記。(b) **既知の気液差**: 298.15 K の L が CEA の H2O と H2O(L) の絶対エンタルピー差と一致 (datum 不変性とは別に組合せの正しさを見る)。
  - (c) float: `condFloat=0/1` の L・`dL/dT` が区切り両側・表範囲外で既存 `tests/unit/test_cond_float.cpp` の基準 (L 相対 2e-6 ほか, `:64`) を満たし、湿潤反転試験 (`:195`) が通る。
  - (c′) 200 K (気相の区間境界) とその両隣の float、境界を含む表小区間で L と `|ΔL′| ≤ 2e-4|L′| + 0.1 J/(kg·K)` を直接判定 (`test_cond_float.cpp:50-58` の接続点除外に頼らない)。表の分割は必須としない (気相の外挿は h・h′ 連続) — 判定で決める。
  - (d) `g>0` の場で保存エネルギーを datum 変換した後、T・P・二相音速が保たれる。(e) 湿り場の種変換・restart 前後で T とエネルギーが整合。
  - (e′) **気相同一・L だけ異なる変換の 0 step A/B** (diagnose): 150 K・N2/H2O 0.95/0.05・g 0.01 の固定状態で、A = 気相差のみの現行式 (補正 0・T −0.005802 K を再現)、B = 気液を含む全差で補正 +4.6592 J/kg。**B の合格: 補正誤差 ≤1e-6 J/kg、T 相対誤差 ≤1e-8**。旧記録 (液相情報なし) からの変換は移行手順か拒否。
  - (f) 湿潤回帰の基準 run (正確な run パス・バイナリ・判定区間)・onset の抽出定義・g の報告量・記録する変化量を**実装前に固定**し (不変を合格にしない) (候補: case/44 va3 入口 Tt 分布 noneq の `run_0127` 系; `run_0510` は g≡0 なので不可)、系列を `check_quasisteady --series-csv` で判定する。
- (旧 V7 の datum 不変性) `thermoHrefTemp` を 0 / 298.15 / 任意値に変えても $L(T)$ が 120–400 K の全点で相対 1e-12 以内で不変。
  新 $L$ と現行 `h2o_latent` の差を 150–373 K で記録 (200 K 以上は係数同一なので丸め程度、200 K 未満は外挿規約の統一分)。凝縮 run (case/44 va3 入口 Tt 分布の noneq) で onset・g の変化量を記録。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan (2 回目, 改訂 3 点に集中) | 2026-09-27 | [2026-09-27-thermophysics-solver-owned-species-db-plan-2.md](../../notes/reviews/2026-09-27-thermophysics-solver-owned-species-db-plan-2.md) | GO-with-changes, C0/M5/m1 | Blanc による併用可・同一 datum の潜熱共通化・N2 保留はいずれも支持。全件採用: M1 縮約拡散の明記・補正後流束の差 (+0.50 %) を近似誤差として記録・複数 lump の式を実装前に確定 → §4.4・V4。M2 気液ペアの基準契約・外部 DB 拒否 → §4.8・#10・V1(e)・V7(b)。M3 潜熱の全経路 (表・float・dL/dT・Python) → §4.8・#10・V7(c–e)。M4 V4 収支と V7 湿潤回帰の合否を実装前に固定 → V4・V7(f)。M5 区切りの段差は現行段差を基準に増分で判定 → V3。m6 N2 検証先は case/34 (case/28 は誤り)・#12 は H2O の阻害条件にしない → §4.8・#12 |
| plan | 2026-09-27 | [2026-09-27-thermophysics-solver-owned-species-db-plan.md](../../notes/reviews/2026-09-27-thermophysics-solver-owned-species-db-plan.md) | GO-with-changes, C0/M7/m1 | 全件採用。M1 (記録を保存場に結び付け全入口で照合) → §4.3・#3・V1。M2 (lump と kinetic 拡散) → 当初採用 (拒否) したが 2026-09-27 ユーザ指摘で**改訂**: lump 外の種との二元係数は Blanc の法則で厳密に作れるので併用可、失うのは lump 内部の差動拡散だけ (§4.4・#7・V4)。M3 (canonical ID・相・LJ 有無・AIR 互換) → §4.1・#4。M4 (後処理・restart reader の移行) → §4.6・#8・V6 (`total_quantities.py:119` を当方で確認)。M5 (float 経路・全区間 datum・上限) → §4.2・#6・V3/V3f。M6 (潜熱共通化は datum 設計が要る) → 当初は後続 plan へ切り出したが、2026-09-27 ユーザ方針 (液相を気相と同じ datum でシフトし差を保つ) で datum 整合が決まったので本 plan の #10・§4.8・V7 に戻す (反例は同じシフトを掛けない置換で、本方針では起きない)。M7 (V0 依存・PASS ケース・NS の事前指定) → §6 冒頭・V5。m8 (仕様文書を実装前に) → #2。判断役 (codex) 自身の指摘で却下が無いため、採否の別途諮問は省略 |

## 7. 影響範囲

- `solver_density_cuda/input/speciesDB.{cpp,hpp}`, `cuda_forge/thermo_d.{cu,cuh}` (区間可変・輸送展開), `cuda_forge/condensationProperties_d.cuh` (#8), 新規 `solver_density_cuda/data/species/`
- `solver_density_cuda/tools/{forge_species.py, interp_field.py, convert_species_field.py, cea_thermo_to_species_db.py, restart_field.py}`
- `design/forge_design/gas/{semiperfect.py, composition.py}`, `design/forge_design/evaluate/runner_axismach.py` ほか runner
- 既存 run: `species_db.yaml` 付きの旧 run は引き続き読める (外部 DB 上書き経路)。restart 照合は旧 run で「照合不能」になり得る (§4.3)

## 8. 完了条件

- [ ] 関連 `methods/` の現在仕様を更新済み
- [ ] 実装・検証完了 (本計画の §6 を満たす)
- [ ] codex レビュー 2 回 (`plan` / `result`) を §6.1 に記録し、Critical / Major の採否を残作業表に反映済み
- [ ] 本計画の `status` を `done` に変更し、§9 に変更ログを記載
- [ ] ファイルを `plans/active/` → `plans/accepted/` へ移動
- [ ] [`plans/README.md`](../README.md) の一覧を同期

## 10. 未確定事項

- ~~設計 runner の NS/SST の既定を Sutherland から種ごとの輸送物性へ切り替えるか~~ → **決着 (2026-09-27 ユーザ)**: semiperfect TP の NS/SST は `viscMethod: 2` を既定 (`gas.transport` 必須)。CPG は Sutherland のまま。δ* 生産レシピ (case/42・44・45) の NS 結果は変わるので、再実行時に変化量を記録する。
- ~~輸送の GPU 性能と viscMethod≠2 併用~~ → **決着 (2026-09-27 ユーザ)**: 性能は (a) 種別 μᵢ・λᵢ と組 ηᵢⱼ を ln T の表 (区分 3 次) にして float で引く (凝縮 `condFloat` と同じ方式)。表の刻み・範囲・精度基準は codex 諮問のうえ #5t2-3 で固定。`physProp.transport` と `viscMethod ≠ 2` の併用は起動拒否で確定。
- H2O の 373 K 未満の比較 (2026-09-27, `notes/investigations/2026-09-27-cea-vs-forge-properties/h2o_lowT_options.py`・`_result.txt`): 案 B (CEA ≥373.2 K + IAPWS 253–373 K + 冪外挿) に対し、単成分で CEA 外挿は μ +1〜+23 %・λ +12〜+88 % (300→200 K)、LJ (現行) μ +26〜+33 %・λ +59〜+74 %、**LJ + 双極子補正 (Brokaw, δ* 1.22) は μ −4〜+7 % だが λ は +25〜+36 % (修正 Eucken が極性分子に合わない)**。混合物 N2 + 蒸気 (CEA 混合則) では μ の差は全案 ≤0.1 % (X_H2O ≤ 0.06)、λ は X_H2O 0.017 で ≤1 %・0.06 で ≤3.5 %・0.2 で最大 12 %。当方の推奨は案 B (水は IAPWS が標準)。ユーザ判断待ち。
- ~~輸送モデルの既定値と kinetic モードの混合則~~ → **決着 (2026-09-27 ユーザ, §4.3c)**: 既定は検証まで `kinetic`・検証後 `cea`、kinetic モードも CEA 形の混合則 (ηᵢⱼ は LJ から)。残る論点: キー名、CEA モードで CEA に無い種/範囲の扱い。
- 凝縮域 (200–373.2 K) の希薄水蒸気の μ・λ: CEA `trans.inp` に無く、LJ は 373.2 K で CEA より +31.6 %、CEA 最低区間の外挿は ~219 K 以下で dμ/dT<0。参照データ (IAPWS 等は 200 K での妥当性を確認要)・許容誤差・接続規約を決める (codex diagnose 2026-09-27)。
  当方の検算 (2026-09-27): CEA 最低区間 (373.2–1073.2 K) を下へ外挿した値を IAPWS の希薄気体項 (粘性 IAPWS 2008 μ₀、熱伝導 IAPWS 2011 λ₀; 公式の適用域は 253.15 K 以上) と比べると、μ は 300 K +1.1 %・273 K +2.6 %・250 K +4.6 %・200 K +12.6 %、**λ は 300 K +12 %・273 K +19 %・250 K +29 %・200 K +83 %** (λ の外挿は ~250 K 以下で T を下げると増える非物理な形)。LJ (現行) は μ で +30 % 前後。373.2 K では CEA と IAPWS が μ −0.4 %・λ +3.6 % で接続できる。混合物への影響は X_H2O に比例して薄まる (凝縮域では蒸気のモル分率は数 % 以下)。
- N2 の潜熱を H2O と同じ方式 (液相を気相と同じ datum で持ち差で L) に統一するか、現行の L フィット方式のままにするか (2026-09-27 ユーザ「今後判断」)。
  判断材料: H2O (#10) の実装と V7 の結果、現行 N2 方式の既知の不整合 (L フィットが液比熱を暗黙に決める)、空気凝縮 run への影響の見込み。

- ~~凝縮域の輸送物性・拡散で液相をどう扱うか~~ → **決着 (2026-09-27)**: 別 plan [`condensation-two-phase-transport.md`](condensation-two-phase-transport.md) へ移管 (懸濁は無視、μ・λ は気相組成、拡散は蒸気の勾配・液は分子拡散なし・微小な液 Sc は任意・乱流は同じ Sc_t)。

## 9. 変更ログ

- `2026-09-30` — ユーザ決定: 熱物性も CEA thermo.inp から生成してソルバ内蔵にする。diagnostician 諮問 (ブリーフ `notes/reviews/briefs/2026-09-30-cea-thermo-builtin.md`): 区間可変は schema 変更として扱い 2 区間の直列化はバイト不変、ハッシュは段 3 で一斉に 1 回だけ変える、61 種 (5 CFC は thermo.inp に無い)・`e-` 除外・CEA `Air` は取り込まない → §4.9・§5.1 #13-0〜13-3 (#6b は統合)。
- `2026-09-30` — §5.2 項 2 を訂正: EXH の構成種 6 種はソルバ内蔵でない (legacy_builtin: design のみ)、外部 DB で渡す。
- `2026-09-30` — ユーザ決定: SERN の lump 名衝突は案 (a) (SERN 側を `AMB` に改名; 当初案 `EXT` から同日ユーザ指定で変更)。§5.2 を現状 (厳密化・輸送物性必須) に合わせて改訂。
- `2026-09-28` — #10 実装 (§5.1 #10 行)。振る舞いの変化: H2O の L が 200 K 未満で変わる (150 K −466 J/kg)、旧 H2O 凝縮 TP 場の restart は拒否 (移行 `--src-latent legacy-v0`)、CPG の H2O 凝縮は記録を持たないので 200 K 未満の L が黙って変わる、凝縮 ON で外部 DB の気相 H2O が内蔵と違えば起動拒否。
- `2026-09-27` — #10 を codex diagnose に諮問 (`notes/reviews/2026-09-27-h2o-latent-datum-diagnose.md`): 湿り場変換の液相項補正を §4.8 に追加、V7(c′)(e′)(f) を具体化、#10 を実装可 (O) に。
- `2026-09-27` — #3c 残を実装: Python ツールも既定で厳密 (未検証 SRC・宛先解決不能・記録破損で停止、許可はその実行だけ、許可時は属性なし)。設計 runner の TP 準備は `FORGE_BIN` (`--resolve-species` 対応バイナリ) が必須になった。残: 印付きの場の扱いのソルバ/ツール不一致 (#3c 行)。
- `2026-09-27` — #3c と案 C の既定切り替えを実装 (`feature/species-transport`): 属性なしの場は既定で停止; `viscMethod: 2` は `physProp.transport` 必須 (無ければ起動時エラー、旧 kinetic 経路をセル・壁から削除)。`viscMethod: 0/1` と transport ありは 0 step で旧バイナリとビット一致。影響: transport なしの `viscMethod: 2` の config (case/05・27・28・50・55・56) はこのブランチのバイナリでは起動しない (旧結果は e2daaba8 までのビルドで再現)。
- `2026-09-27` — 既存の semiperfect NS 問題 YAML 22 件に `gas.transport` を追記 (N2/O2/AR/CO2 は cea、H2O は custom:h2o_iapws_cea_v1)。
- `2026-09-27` — #9b 実装 (上表)。ユーザの PC 停止で git オブジェクト 5 個が空になり、GitHub から取り直して復旧 (fsck クリーン)。
- `2026-09-27` — ユーザ決定: 設計 runner の semiperfect TP の NS/SST は viscMethod 2 + physProp.transport を既定 (gas.transport 必須、CPG は Sutherland のまま)。
- `2026-09-27` — ユーザ指示で runner の輸送指定 (#9b) を段 3 より先に。runner の NS が Sutherland (空気) であることを確認し、`gas.transport` がある問題だけ新経路 (既定の切り替えは §10)。
- `2026-09-27` — #5t2-3 表引き化を実装・検証 (上表)。精度・A/B・既定経路・性能の合格条件をすべて満たし、表引きは旧 Wilke 経路より速い (n12 で 0.20 vs 0.31 ms)。次は段 3 (NS 統合; 凝縮 plan の気相組成を先に固定)。
- `2026-09-27` — 表引き化の設計を codex diagnose で確認し全件採用: 区間境界で表を分割、区間選択は元の T、合格条件と性能の測り方を事前固定 (#5t2-3)。
- `2026-09-27` — ユーザ決定: 輸送の性能対策は表引き (a)、viscMethod≠2 との併用は拒否で確定。
- `2026-09-27` — #5t2-2 段 2 を実装: 正しさの合格条件はすべて PASS、性能は提案値 10 % を超過 (実種 5/12/32 で +14/+94/+800 %)。対策と viscMethod≠2 併用の扱いは §10。
- `2026-09-27` — 段 2 (GPU 接続) の設計を codex diagnose で確認し全件採用: モル基底の展開、独立参照による合格条件、入力丸めと評価誤差の分離、μ・λ 利用箇所の確認表 (#5t2-2)。
- `2026-09-27` — ユーザ決定: H2O の冪外挿は 253.15 K (IAPWS 適用域下端) から、傾きは IAPWS の式から計算。実装・試験 ALL PASS。
- `2026-09-27` — #5t2 段 1 を実装 (上表)。設計の穴: IAPWS μ₀ は 202.17 K で最小・134.12 K に極 (§4.3c の「~134 K まで単調」は誤りで訂正)。冪外挿の開始温度は 253.15 K 案をユーザに提案中。
- `2026-09-27` — **作業ブランチを移動**: 並行セッション (D-7275 すきま加熱) が `feature/gap-heating-precision` の同じ作業ツリーで作業中のため、以後の種 DB・輸送の作業は `feature/species-transport` (作業ツリー `/home/sano/work/forge-species`, d8c11f52 から分岐) で行う。既定の切り替え (#3c、`viscMethod: 2` の置き換え) もこのブランチで進め、並行セッションは区切りでマージする。
- `2026-09-27` — ユーザ決定: `viscMethod: 2` は案 C (置き換えで実装、既定の切り替えは並行セッションの区切りを待つ)。H2O は特別なモデルとして `custom:h2o_iapws_cea_v1` (`custom:` 名前空間 + 版)。
- `2026-09-27` — 輸送の実装設計を codex diagnose で確認: 合格条件を全 CEA と選択モデルで分ける、ηᵢⱼ の規約、記録の追記条件、3 段の実装順 (#5t2)。`viscMethod: 2` の扱いはユーザ判断待ち。
- `2026-09-27` — ユーザ決定 (確定): H2O は 600 K 付近でつなぐ (IAPWS ↔ CEA、500–700 K smoothstep)。全温度 IAPWS 案は保証範囲を理由に取り下げ。
- `2026-09-27` — ユーザ決定: H2O の μ・λ は全温度で IAPWS 希薄気体の式 (CEA とのつなぎは廃止; 4000 K まで CEA と μ ±1 %・λ ±3.6 %)。
- `2026-09-27` — ユーザ決定: H2O は案 B。つなぎを 500–700 K の smoothstep (log 空間) に、kinetic の λ は Warnatz 式を候補に (§4.3c)。
- `2026-09-27` — ユーザ決定: 輸送物性の出所は全体モードでなく種ごとに必ず指定 (§4.3c 改訂)。H2O 低温の案比較を §10 に記録。
- `2026-09-27` — #5t ① 混合則の A/B 完了: CEA frozen 混合則 (B) が FCEA2 と全 16 状態で ≤0.013 %、現行 (A) は最大 13 %。B を採用。次は ②③ (輸送カーネルの変更、実装前に上位へ諮る)。
- `2026-09-27` — ユーザ承認: 既定は検証まで kinetic → 検証後 cea、kinetic モードも CEA 形の混合則。#5t ① に着手。
- `2026-09-27` — ユーザ決定: 混合則は必ず CEA 形へ、μ・λ は `cea` / `kinetic` の 2 モード選択、H2O 低温は別途ケア (§4.3c)。
- `2026-09-27` — 輸送の正本設計を codex diagnose で確認: 解決順維持、CEA の frozen 混合則全体を採用、H2O 低温は未決 (§10)、実装順 #5t。
- `2026-09-27` — #9 (axismach) と #8a 完了 (上表)。**ユーザの目的 (problem のモル分率・lump → config に lump のモル分率、生成 species_db.yaml なし) が case/44 の設計チェーンで実現** (V5・V6 PASS)。SERN は lump 名 AIR の衝突で未切替。
- `2026-09-27` — #6a 完了 (上表)。
- `2026-09-27` — #6 を #6a (区切りが揃う種の起動時合成、cuda_forge 不変) と #6b (区間可変) に分割。ユーザの目的 (config に lump の中身をモル分率で書く) を #6a + #8 の一部 + #9 で先に実現する。
- `2026-09-27` — #4 完了 (上表)。He の MW と thermo.inp の差を #5 に追加。
- `2026-09-27` — ユーザ方針「CEA を正解とする」を §4.3b に記録。CEA に LJ は無く `trans.inp` が μ・λ のフィットを持つこと、現行 LJ 由来の値が H2O で μ +23 %・λ +47 % (600 K) ずれることを確認。#5 を拡張。
- `2026-09-27` — #4 の委譲で implementer が停止 (編集なし): C++ 内蔵と Python の共通 5 種は一致、LJ は Python `LJ_PARAMS` と CEA ツールの表が 6 種で食い違い、canonical 化は Python の大文字化・ハッシュ内の名前と衝突。#4 の範囲を「現行内蔵経路の値の移設、名前は案 (a)」に絞り、LJ の寄せ先は #5、canonical 化は #8 へ。
- `2026-09-27` — #3b 実装・検証 (上表)。runner_axismach の段間引き継ぎを interp_field → restart_field に変更 (tooling-design-problem-campaign-recipe #5 と同じ変更)。#3c (既定を厳密へ) は並行セッションの TP run と case 内 IC スクリプトの移行を見てから。
- `2026-09-27` — #3a 実装・検証 (上表)。並行セッション (case/56 TP) の継続を止めないよう、属性なしの場は #3b 完了まで警告で通す過渡期の既定にした (係数不一致は常に停止)。
- `2026-09-27` — #3 の委譲で implementer が「属性の無い場を拒否すると新規 TP 初期場も全部止まる」穴を発見 (編集なし)。codex diagnose で補強した案 A を採用し §4.3・#3a/#3b・V1 (f) を改訂。V0 取得済み (§6 冒頭)。
- `2026-09-27` — codex plan 段レビュー 2 回目 (改訂 3 点, GO-with-changes C0/M5/m1) を全件採用。凝縮域の輸送物性で液を蒸気扱いしている件を §10 に追加 (ユーザ判断待ち)。
- `2026-09-27` — N2 の統一はユーザが保留に戻した (「変えなくてもいいのかな、今後判断」)。§4.8・#12 を「判断」項目に、§10 未確定事項に追加。
- `2026-09-27` — ユーザ決定 (のち保留に変更): N2 の潜熱も同じ方式 (液を気相と同じ datum で持ち差で L) に統一する、ただし H2O の後で優先度低 (§4.8 末尾, §5.1 #12)。
- `2026-09-27` — ユーザ指摘で 2 点改訂: (1) lump の拡散は平均 LJ でなく Blanc の法則で二元係数を作り kinetic 拡散と併用可 (M2 の拒否を撤回)、(2) 潜熱は液相を気相と同じ datum でシフトして差で作る (§4.8, 本 plan 内に戻す)。
- `2026-09-27` — codex plan 段レビュー (GO-with-changes, C0/M7/m1) を全件採用し §2–§6 を改訂。最初の実装は「保存場に結び付いた解決済み記録と全 restart 経路の照合」に限定。潜熱の共通化は後続 plan。
- `2026-09-27` — 初稿。ユーザ要望と codex diagnose 諮問 (`notes/reviews/2026-09-27-solver-owned-species-db-diagnose.md`) の推奨から起票。
