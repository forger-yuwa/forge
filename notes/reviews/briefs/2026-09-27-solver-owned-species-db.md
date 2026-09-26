# 諮問ブリーフ: 化学種の熱物性 (CEA) をソルバが持ち、run には組成 (モル分率) だけを書く構造にしたい (2026-09-27)

## ユーザの気持ち (原文の要旨)

- 「`species_db.yaml` が problem から生成されて、わけのわからない値 (擬似種 MIXDRY の NASA-9 係数・MW・LJ・`atoms: {N: 1.41775, O: 0.565246, …}`) が入ること、正直やめてほしい」
- 「入れることができるようにしている (外部 DB で上書きできる) のは良い。でも CEA の情報はソルバーで持っておいて、モル分率を指定できる仕組みにしたい」
- 「Tlo/Tmid/Thi もガス種によって変わることもあるでしょう」
- 「`atoms` は意味が分からない。どういう役割?」「H2O 液はどんな情報を持っている? cp は潜熱から逆算している?」
- 背景: 同日、設計ツールの入力を problem / campaign / recipe に分ける plan を起票済み
  (`plans/active/tooling-design-problem-campaign-recipe.md`)。その議論で「生成された config にわけのわからない数値が入る」ことにユーザが違和感を示した。

## 問い (1 つ)

**化学種の熱物性の正本をソルバに置き、run の入力は「種名 (または lump の名前と中身のモル分率) + 境界の組成 (モル分率)」だけにする構造へ移るべきか。
移るなら、lump (擬似種) をどこで・どう作るか (ソルバが実行時に構成種の和として評価 / ソルバ起動時に係数を合成 / 設計側で合成して渡す現状維持) を 1 つに絞って推奨し、
最初の一歩として何を 1 つだけ変えるべきかを示してほしい。**

候補 (縛られなくてよい):
- (a) 現状維持 + 見せ方だけ改善 (生成 DB を隠す・コメントを充実)
- (b) ソルバ内蔵 DB を CEA `thermo.inp` 直読み (または生成済み全種表) に一本化し、config に `species: [{name: MIXDRY, lump: {N2: .., O2: .., AR: .., CO2: ..}, basis: mole}, H2O]` のように書く。lump の係数合成はソルバ起動時 (温度区間が揃う種だけ)
- (c) (b) と同じ入力で、lump の熱物性を**実行時に構成種の重み付き和として評価**する (温度区間が違う種も畳める。輸送方程式は lump 数のまま、熱物性評価のコストだけ構成種数に比例)
- (d) lump をやめて常に `full` (全種輸送) にする (case/44 で 2 種 → 5 種で step +23 %, 回帰は full ≡ lumped が 5 桁一致)

## 観測事実 (2026-09-27 セッションで確認, `ファイル:行`)

1. **CEA 係数が 4 か所にある (完全一致しない)**
   - ソルバ内蔵 DB `solver_density_cuda/input/speciesDB.cpp:83-157` (`speciesDB_builtin()`): N2, O2, AR, CO2, HE, H2O, AIR (擬似空気 cp/R 3.5 一定) の 7 種。
   - 設計側 `design/forge_design/gas/semiperfect.py` の `SPECIES_NASA9`: N2 O2 CO2 H2O AR H2 OH H NO O CO の 11 種 (CEA 転記)。
   - `solver_density_cuda/tools/cea_thermo_to_species_db.py`: CEA `thermo.inp` 直読み (ローカル `.venv-cea/nasa_cea/` 同梱)。
   - 凝縮の潜熱 `solver_density_cuda/cuda_forge/condensationProperties_d.cuh:239-262` (`h2o_latent`) が H2O 気相 200–1000 K 係数を**再度ハードコード**。
   - 既知の不一致 (accepted plan `plans/accepted/thermophysics-cea-mole-fraction-species.md:66-67`): H2O MW 0.0180153 vs 0.01801528、AR 高温域 a0 0 vs 20.105。
   - 同 plan は「`thermo_d.cu` 内蔵 DB の変更」をスコープ外にした (`同:49`)。その結果、設計側が DB を合成して run に `species_db.yaml` を書き、ソルバはそれで内蔵を上書きする構造になった。
2. **ソルバの DB 読み込み** (`speciesDB.cpp:158-235`): `speciesDBFile` があれば同名種を上書き/追加。読むキーは `MW, LJ_sigma, LJ_eps_kB, Tlo, Tmid, Thi, nasa9_low[9], nasa9_high[9]` のみ。**`atoms` は読まない**。
   温度区間は種ごとに持てる (`thermo_d.cuh:46-50`) が **2 区間 (low/high) まで**。範囲外は端でクランプ + 定 cp 外挿 (`thermo_d.cuh:101-107`)。
3. **`atoms` の役割**: 設計側の元素質量分率診断 `element_mass_fractions` (`design/forge_design/gas/composition.py:248-255`) と `species_meta.yaml` の記録 (`composition.py:566`) のみ。
   呼び出しは単体テスト (`design/tests/run_gas_tests.py:176`) だけ。lump の `atoms` は構成種のモル分率加重平均 (`composition.py:269-274`) なので非整数になる。
   accepted plan §2 (`:52-53`) は SERN の `full` で「元素質量分率から混合分率を作る」用途を想定していた。ソルバの計算には一切効かない。
4. **lump の合成** (`composition.py:259-280`, `lump_entry`): NASA-9 係数を lump 内のモル分率 (= y·MW_mix/MW_k) で線形混合。組成が空間一様なら cp/h は厳密。
   温度区切りが 200/1000/6000 K でない種は拒否 (`composition.py:265-267`)。**LJ パラメータは質量分率の単純平均** (`composition.py:276`)。
   ソルバの TP 粘性は Chapman–Enskog 単成分 + Wilke 混合 (`thermo_d.cuh:346-435`) なので、NS で lump を使うと lump の粘性は「平均 LJ の仮想分子」になり、構成種の Wilke 混合とは一致しない (Euler では不使用)。
   s° の混合エントロピー項は lump 係数に入らない (Euler/NS の流れには不使用のはず; 要確認)。
5. **入口組成のモル分率入力は既にソルバが受ける**: bcond `floats: {X0: .., X1: ..}` を MW で Y に換算 (`procedures/solver-settings.md:226`, `boundaryCond.cpp:90`)。
   case/44 `run_0509`–`0511` で lump の X (MIXDRY 0.939014829007 / H2O 0.0609851709927) を与え、起動ログの Y_H2O 0.03769539669 が質量分率入力 (0.0376953964) と 1e-8 で一致。
   ただし「lump の中身のモル分率」を config に書く手段は無く、lump の中身は生成された `species_db.yaml` の係数に埋め込まれている。
6. **液相 H2O の物性はすべてソルバにハードコード** (`condensationProperties_d.cuh`; species_db には入らない):
   - 飽和圧: Murphy & Koop 2005 の過冷却液式 (`:205-212`, 123–332 K)。
   - 液密度: `1000 − 0.12(277 − T)`, 下限 920 kg/m³ の「ゆるい近似」(`:215-221`)。
   - 表面張力: IAPWS R1-76 を過冷却へ外挿 (`:268-275`)。
   - 潜熱: `L = h_v − h_l`。h_v は CEA H2O 気相 (200–1000 K) を再ハードコード、h_l は CEA H2O(L) (273.15–373.15 K)。**273.15 K 未満は h_l を cp_l(273.15)=4228 J/kgK 一定で線形外挿、373.15 K 以上は h_l をクランプ (実質 cp_l=0)**、L は 1.5–3.5 MJ/kg にクランプ (`:239-266`, `methods/condensation.md:479-490, 583-592`)。
   - 二相 EOS の液内部エネルギーは `e_l = e_v + R_v T − L` (= h_l, pv 項を落とす) (`condensationEOS_d.cuh:50-65`)。
   - **cp は潜熱から逆算していない (H2O)**: L のほうを液相 CEA エンタルピーから作っている。液の比熱は h_l の傾き (273 K 以上は CEA、未満は 4228 一定)。
     一方 **N2 (CPG 経路) は液物性を独立に持たず、L(T) のフィットから c_l = c_p,v − dL/dT が暗黙に決まる** (`methods/condensation.md:564-574`)。
   - 既知の課題: h_v の二重ソース (種 DB は定 cp 外挿、`h2o_latent` は生多項式; 120 K で 2.39 kJ/kg の差, `methods/condensation.md:587-590`)。
7. 性能: case/44 で 2 種 (lump) 2.45 ms/step、5 種 (full) 3.01 ms/step (+23 %) (`case/44.vitiated_air_wt/README.md` の species 節)。`THERMO_MAX_SPECIES` の上限あり (`speciesDB.cpp:196-199`)。
8. 再現性の契約: restart は species 署名 (順序・MW・係数・区切り・datum・tracer) 不一致で拒否 (accepted plan の完了記録)。DB の出所を変えると既存 run の restart 互換に影響する。

## 期待値と出典

- 熱物性の正本は CEA (McBride–Gordon 2002 `thermo.inp`)。forge TP は CEA 凍結流と ≤0.05 % (case/44 `run_0091`)。
- 移行しても case/44 の報告量 (ṁ, 出口 M/T, 軸 M) が現状と許容差内で一致することは検証条件にする (許容差は反復ノイズ床から事前に決める)。

## 実施済みの操作と結果

- 上記 1〜8 はコード・文書の読み取り。新しい run は回していない。
- 同日 case/44 `run_0509`–`0511` (lump + bcond X 入力) は全 NaN 0・報告量 ALL STEADY・残差 plateau (NOT CONVERGED)。

## 仮説 (当方の見立て。棄却してよい)

- H1: (b) か (c) が本筋。熱物性の正本をソルバ 1 か所 (CEA thermo.inp 由来の生成表をソースに埋め込む or 同梱ファイルを直読み) に寄せ、設計側 Python も同じ表を読む。
  lump は config に「名前 + 中身のモル分率」を書き、係数はソルバが作る (起動ログに合成結果を出す)。外部 DB による上書きは残す。
- H2: lump の実行時評価 (c) は温度区間の違いと LJ 平均の問題を同時に消すが、熱物性評価が構成種数倍になる (TP の温度反転 Newton の中で呼ばれる) のでコスト次第。
- H3: `atoms` はソルバに不要。元素診断が必要なら CEA 表の元素欄から後処理で作ればよく、run の DB ファイルに書く必要はない。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md`, `thermo.inp` 全体を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か本ブリーフの番号で示す。

## 読んでよいもの

- `sed -n '1,66p' solver_density_cuda/input/speciesDB.hpp`
- `sed -n '60,240p' solver_density_cuda/input/speciesDB.cpp`
- `sed -n '30,120p' solver_density_cuda/cuda_forge/thermo_d.cuh`
- `sed -n '340,440p' solver_density_cuda/cuda_forge/thermo_d.cuh`
- `sed -n '195,310p' solver_density_cuda/cuda_forge/condensationProperties_d.cuh`
- `sed -n '40,70p' solver_density_cuda/cuda_forge/condensationEOS_d.cuh`
- `sed -n '1,130p' design/forge_design/gas/composition.py`
- `sed -n '240,300p' design/forge_design/gas/composition.py`
- `sed -n '520,580p' design/forge_design/gas/composition.py`
- `sed -n '1,160p' plans/accepted/thermophysics-cea-mole-fraction-species.md`
- `sed -n '479,600p' methods/condensation.md`
- `sed -n '215,232p' procedures/solver-settings.md`
- `sed -n '1,40p' case/44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX/species_db.yaml`
- `sed -n '1,120p' plans/active/tooling-design-problem-campaign-recipe.md`
