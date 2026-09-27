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
- 本 plan の範囲では N2 の現行 L・飽和圧・CPG 経路を**保持**し、#10 の共通化 (ディスパッチ変更) で N2 の結果が変わらないことだけ確認する (2 回目 m6)。

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
| 3c | 未検証の既定を厳密へ切り替え | #3b 完了後、`main.cpp` の既定を「属性なしは停止」に (環境変数 `FORGE_REQUIRE_VERIFIED_SPECIES` を撤去し `FORGE_ALLOW_UNVERIFIED_SPECIES` だけ残す)。`procedures/solver-settings.md` と試験 (0t) を更新 | O |
| 3b | ~~引き継ぎ・種変換・IC 生成の入口~~ | 完了 2026-09-27: 共通 API (`tools/forge_species.py` の `find_forge`/`resolve_species`/`plan_inherit`+`commit_inherit`/`check_ic_against_record`/`stamp_new_field`/`plan_convert`)、`restart_field.py`・`interp_field.py`・`convert_species_field.py`・`runner_sern.py` (`restart_by_index`・領域 IC・`warm_from_run`)・`design/forge_design/evaluate/ic.py` (`stamp_isentropic_ic_species`)・`runner_axismach.py` (IC 付与、段間を restart_field へ)。試験 `tests/unit/test_species_attrs_entry.py` (19 項目)・`design/tests/run_species_attrs_ic_tests.py` (7 項目) ALL PASS: V1 (a)(d)(f)、未検証 SRC で DST 属性を消す、`--force-species` は属性なし、種変換。既存試験も PASS。配管 run `case/44.vitiated_air_wt/run_0520_species_attrs_runner` (IC→段間 2 回→本段 60 step で unverified=0、restart_field はビット一致; 同じ場を interp_field で写すと roUx 22950/23725 点で最大 1.4e-5 相対の差)。解釈として決めた点: 宛先を解決できない (旧バイナリ等) ときは過渡期は警告で属性なし・strict で停止 / 属性はあるが記録なし・完全性不一致の SRC は過渡期でも停止 (`--force-species` で属性なし) / SRC 記録を DST の隣へ複製 / interp_field の内蔵種の照合不能は過渡期は警告 (#3a の残課題を解消)。未了: runner 経路 (段の手順が run_0509 と違う) の本段 24000 step の数値基準、`warm_from_run` の実 run、`runner_wt`/`runner.py`/case 内 `gen_*_ic.py` (CPG か範囲外) | O |
| 4 | ~~共通データ化 (値は変えない)~~ | 完了 2026-09-27: `solver_density_cuda/data/species/forge_species_v1.yaml` (13 種、canonical ID・別名・相・2 区間係数・LJ・元素組成・出典・`deviations`・過渡欄 `legacy_builtin`)、C++ はビルド時埋め込み (`cmake/embed_species_data.cmake`, `speciesDB_builtin()` のハードコード撤去)、Python `semiperfect.py`/`composition.py` は同ファイルを読む (名前は案 (a) のまま)。値のビット一致 `tests/unit/test_species_data_bitexact.py` ALL PASS (98 項目, 変異試験で FAIL を確認)。`--resolve-species` のハッシュは新旧で一致 (case/44 4378b7d78339ba27、内蔵のみ 3 構成) し記録ファイルもバイト一致。既存 species 試験は新バイナリで全 PASS。残: `tools/forge_species.py:69` の `BUILTIN_MW` の写し (#8)、He の atoms を新設 (読み手なし)、既存の `test_solver_config_species.cpp` は fixture の廃止キー `mesh.meshFormat` で FAIL 16 (本件と無関係、未修正) | O |
| 5 | CEA 直読みとの差・20000 K 区間・LJ の寄せ先 | H2O MW・AR 高温 a0・**He MW (0.0040026 vs 4.002602 g/mol, 相対 5e-7)** の寄せ先、**LJ の出典と寄せ先 (Python `LJ_PARAMS` と `cea_thermo_to_species_db.py` の Cantera 由来表が H2 2.827/59.7 vs 2.920/38.0、H 2.708/37.0 vs 2.050/145.0、O 3.050/106.7 vs 2.750/80.0、OH 3.147/79.8 vs 2.750/80.0、NO 3.492/116.7 vs 3.621/97.53、CO 3.690/91.7 vs 3.650/98.10 で食い違う)**、6000–20000 K 区間を有効にするか、CEA 全種を入れたとき現行 `AIR` の alias `Air` と CEA `Air` の衝突をどう解くか。値を変えるなら case/44 と #6 の小型ケースの報告量変化を記録 。**加えて (2026-09-27 ユーザ方針 §4.3b)**: 輸送物性の正本を CEA `trans.inp` にする設計 — 種別フィットの取り込み、相互作用データの混合則への入れ方、CEA 範囲外 (H2O < 373 K など) と CEA に無い種の扱い、拡散係数 (LJ 継続) との整合、NS 結果の変化量の記録 | F |
| 5t | 輸送物性の CEA 化 (段階) | codex diagnose の順: ① **混合則の A/B (CFD 0 step)** — 単成分フィット・MW・T・組成を固定し A = 現行 (φ 共用) / B = CEA の ηᵢⱼ・φ・ψ、T = 400/600/1000/2000 K × X_H2O = 0/0.1/0.5/1 の 16 状態を FCEA2 の μ・frozen λ と比較 (A 全点 0.1 % 以内なら「現行で不足」を棄却、A 失敗・B 全点合格なら B 採用、両方失敗なら GPU 実装へ進まない) → ② 共通 resolver と単成分評価 (ユーザ定義フィット → CEA → LJ、単位・区間端・記録) → ③ CEA 混合・lump 展開・CUDA 評価 (独立 double 基準に double ≤1e-12・float ≤1e-5) → ④ 低温 H2O モデル確定後の NS 検証 (収束・準定常 VERDICT 必須) | F |
| 6a | ~~起動時 lump 合成 + config 指定 (区切りが揃う種のみ)~~ | 完了 2026-09-27: `physProp.species` の mapping 記法 (`input/solverConfig.cpp`, `input/speciesLump.hpp`)、`speciesDB_resolve` での合成・検査 7 種・起動ログ・記録と互換性ハッシュ (lump なしはバイト不変)、Python `forge_species.py` の lump 対応。試験 `tests/unit/test_species_lump_solver.py` ALL PASS (31): **V2 = 生成 DB と係数・MW・LJ 相対 ≤4e-16、cp/h/s° (200–6000 K 1000 点) 相対 ≤8e-16**、既存 config のハッシュ・記録バイト一致、負例すべて拒否。起動確認 `case/44.vitiated_air_wt/run_0521_species_lump_startup` (lump 記法で 200 step、NaN 0)。basis は必須。**lump 記法と外部 DB 版はハッシュが違うので既存場からの restart は照合で止まる — 正しい挙動として採用** (1 ulp の差があり、外部 DB 擬似種の出自は場から確かめられない; 移行は明示許可 1 回)。残: lump の mapping を読めない Python reader (`total_quantities.py`・`convert_species_field.py`・`passive_gate_common.py`・`gen_inlet_profile.py`・`runner_sern._species_signature`・design `probdef`) は #8/#9。過渡期既定で通すときの警告文が環境変数を立てたように読める (#3c で直す) | O |
| 6b | 区間可変 (区切りの違う種を畳む) | §4.2 の和集合区間。`SpeciesThermo`・float 表・datum の区間可変化。`cuda_forge/thermo_d` の変更なので実装前に上位へ諮る。合格は §6 V3・V3f | F |
| 7 | lump の輸送物性展開と Blanc 拡散 | §4.4 (粘性・熱伝導は実種展開、lump を含む二元拡散係数は Blanc)。合格は §6 V4 | O |
| 8a | ~~reader の移行 (記録から読む)~~ | 完了 2026-09-27: 共通読み出し `forge_species.run_thermo` (res 属性の記録 → run の記録 → `speciesDBFile` → `--resolve-species`)、`total_quantities.py`・`convert_species_field.py`・`gen_inlet_profile.py`・`runner_sern._species_signature` を移行。V6 PASS (下の #9) | O |
| 8 | Python 共通 API の残り、canonical ID への移行 | 残り: `forge_species.species_info` の lump MW (`BUILTIN_MW` の写し)・`species_signature` (lump を照合不能扱い)・`interp_field` の署名。 §4.6。加えて (2026-09-27): Python の種名の大文字化をやめ canonical ID + alias 表へ、C++ `ResolvedSpeciesDB::index()` (`speciesDB.cpp:79-85`, 大小文字無視; 重複検査 `:210` も使う) の完全一致化と、tracer・凝縮種など名前で引く箇所の影響調査。**互換性ハッシュには config の名前がそのまま入る (`speciesDB.cpp:461`) ので、canonical 化で既存記録と不一致にならない規約 (ハッシュには canonical ID を入れ、既存記録は移行ツールで読み替え等) を設計してから**。合格は §6 V6 | O |
| 9 | ~~設計 runner の切り替え~~ (axismach 完了、SERN 未) | 完了 2026-09-27 (axismach): `_apply_gas_to_config` は `physProp.species` を lump 記法 (全桁の正規化モル分率) で書き **`species_db.yaml` を作らない**。内蔵に無い種・外部 DB (`gas.species_db`) が内蔵値を上書きする種だけ生エントリを `species_db_external.yaml` に置く (合成物は置かない)。変換器は `FORGE_BIN` と同じビルドのもの (`runner.converter_path()`)。**V5 (i)** `case/44.vitiated_air_wt/run_0522_species_nodb_lumpX_v5` (ref): V0 (run_0509/0513–0515) との最大差 ṁ_in 6.1e-7 相対・ṁ_out 2.6e-7・出口 M 9.5e-7・出口 T 1.2e-4 K・軸 M 出口 **8.6e-6 (許容 1e-5, 反復差の約 4.5 倍; 原因未切り分け)**・軸 M 目標差 3.6e-7 → 許容内。本段区間 NOT CONVERGED (plateau, run_0509 と同型)・series ALL STEADY、NaN 0、メッシュ PASS。**V6**: DB ファイルなしで prepare → 段間 restart_field (記録継承) → 本段 → `total_quantities` (旧経路比 T0 1.2e-15・P0 2.6e-14) → lump→full5 変換 `run_0523_species_nodb_full5_convert` (ρY 保存差 0、変換後 200 step で照合一致・NaN 0) → `gen_inlet_profile` (CSV バイト一致)。**SERN は未切替**: SERN の lump 名 `AIR` がソルバ内蔵の擬似種 `AIR` と衝突し起動時に拒否される (`speciesDB.cpp:329-334`)。lump 名の変更か内蔵 `AIR` の扱い (#5 の Air 衝突と同根) を決めてから | O |
| 10 | 潜熱: 液相を気相と同じ datum でシフトし差で L を作る | §4.8。`h2o_latent` の H2O 気相再ハードコードを撤去し、液相 H2O(L) (共通データの凝縮相エントリ) に**気相 H2O と同じ datum 定数**を適用して `L = h_v − h_l` を作る (気液差は CEA のまま保たれる)。273.15 K 未満/373.15 K 超の液の延長規約は現行のまま。**移行対象**: 潜熱表生成・範囲外退避・二相熱容量/音速 (`dL/dT`)・面流束・二相反転・Python `convert_species_field.py` の潜熱 (§4.8)。外部 DB の気液ペア契約と拒否。合格は §6 V7・V1 (液相だけ変えた restart の拒否)。`cuda_forge` の凝縮カーネル変更なので編集前に上位へ諮る (AGENTS.md エスカレーション 6) | F |
| 12 | N2 の潜熱を同方式に統一するかの判断 (未決) | §4.8 末尾・§10。H2O (#10) の後に、統一するか現行 (Lin フィット + 低温外挿) のままにするかをユーザと決める。統一する場合: 液相 = CEA `N2(L)` 77.352 K 点 + 液比熱モデル、L は差。飽和圧の低温再構成・CPG carrier 経路との整合を決めてから。合格: V7 と同じ datum 不変試験 + 現行 Lin フィットとの差を 45–120 K で記録 + 空気凝縮 run (**case/34 Arthur**; plan condensation-air の検証先。旧記載の case/28 は He/空気同軸ジェットで誤り) の onset 変化を記録。**H2O (#10) 完了の阻害条件にしない**。凝縮カーネル変更なので編集前に上位へ諮る | F |
| 11 | docs 同期 (完了時) | `procedures/solver-settings.md` (`physProp.species` の lump 形、`speciesDBFile` の位置づけ、lump と拡散の制約)、`recommended-settings.md` §3、`design/CAPABILITIES.md` | O |

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
  - (d) `g>0` の場で保存エネルギーを datum 変換した後、T・P・二相音速が保たれる。(e) 湿り場の種変換・restart 前後で T とエネルギーが整合。
  - (f) 湿潤回帰の基準 run・onset の定義・g の報告量を**実装前に固定**し (候補: case/44 va3 入口 Tt 分布 noneq の `run_0127` 系; `run_0510` は g≡0 なので不可)、系列を `check_quasisteady --series-csv` で判定する。
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

- 凝縮域 (200–373.2 K) の希薄水蒸気の μ・λ: CEA `trans.inp` に無く、LJ は 373.2 K で CEA より +31.6 %、CEA 最低区間の外挿は ~219 K 以下で dμ/dT<0。参照データ (IAPWS 等は 200 K での妥当性を確認要)・許容誤差・接続規約を決める (codex diagnose 2026-09-27)。
- N2 の潜熱を H2O と同じ方式 (液相を気相と同じ datum で持ち差で L) に統一するか、現行の L フィット方式のままにするか (2026-09-27 ユーザ「今後判断」)。
  判断材料: H2O (#10) の実装と V7 の結果、現行 N2 方式の既知の不整合 (L フィットが液比熱を暗黙に決める)、空気凝縮 run への影響の見込み。

- ~~凝縮域の輸送物性・拡散で液相をどう扱うか~~ → **決着 (2026-09-27)**: 別 plan [`condensation-two-phase-transport.md`](condensation-two-phase-transport.md) へ移管 (懸濁は無視、μ・λ は気相組成、拡散は蒸気の勾配・液は分子拡散なし・微小な液 Sc は任意・乱流は同じ Sc_t)。

## 9. 変更ログ

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
