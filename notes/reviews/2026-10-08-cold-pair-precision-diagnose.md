# codex 諮問 (diagnose): cold-pair-precision

- **brief**: [`notes/reviews/briefs/2026-10-08-cold-pair-precision.md`](../../notes/reviews/briefs/2026-10-08-cold-pair-precision.md)
- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **date**: 2026-10-08
- **commit**: `fb3144b6` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 6.6 min, rc=0
- **結論**: **まず座標読み込みだけを `stof`／`stod` で変える FP64 変換器の CFD 0 step A/B を行い、入力精度の保持を確認してください。**
- **extra**: `case/45.isobutane_m6_d155/cold_pair_mesh.py`, `case/45.isobutane_m6_d155/cold_pair.py`, `solver_density_cuda/flowFormat.hpp`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：**FP64 を使う方向は採用しますが、「typedef 4 行だけ変更した (A)」の投入は却下します。座標読み込みに FP32 への丸めが残っています。**

| 重大度 | 採否 | 根拠と対案 |
|---|---|---|
| **Major** | **現状の (A)：却下。入力経路を直した FP64 版：採用** | [`gmshReader.hpp:504`](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/gmshReader.hpp:504) は座標を `stof` で読みます。**生産元の commit `e2696d8f0` でも同じ**でした。`geom_float = double` にしても、float に丸めた値を double に格納するだけです。対案は座標を `stod` で読み、**変換器と forge の両方**を同じ修正ソースからビルドすることです。 |
| **Minor** | **msh 出力桁数の変更：採用** | [`mesh2d.py:269`](/home/sano/work/forge-integ-1005/design/forge_design/meshing/mesh2d.py:269) は `.10g` 出力です。現行冷却格子の再生成で、このテキスト丸めだけでも第一層厚が最大 **0.2245%** 変わりました。`.17g` にして入力精度を保持してください。 |
| **Major** | **品質ゲート：修正が必要** | [`ns_n012.py:217`](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/ns_n012.py:217) の `"PASS" in line` は **`SOFT-PASS` も通します**。品質ツールは外れ値が0.1%以下なら exit 0 です（[`check_mesh_quality.py:280`](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_mesh_quality.py:280)）。今回の補正版格子を FP32 に丸めると、AR > 5000 は **327セル、0.0578%**。現ゲートでは通り得ます。V-c45 の条件どおり、厳密な `VERDICT: PASS` を要求してください。 |
| **Major** | **縮流部の AR 補正：採用。ただし AR 5000 例外の適用は要再検証** | [`cold_pair_mesh.py:57`](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_pair_mesh.py:57) には既に勾配補正が入っています。しかし生成格子は半径方向配置です（[`mesh2d.py:220`](/home/sano/work/forge-integ-1005/design/forge_design/meshing/mesh2d.py:220)）。再生成した壁第一層の内角は **73.46〜129.60°**。全域を「壁法線・無スキュー」と呼ぶのは誤りで、[`AGENTS.md:137`](/home/sano/work/forge-integ-1005/AGENTS.md:137) の例外を自動適用できません。対案は非直交部を既定 AR 制限に収める設計とし、5000 は条件を満たす層に限定することです。 |
| **Major** | **壁解像ゲート：修正が必要。設計目標 0.8：採用** | [`check_wall_resolution.py:378`](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_wall_resolution.py:378) は**節点数割合を「面積割合」と表示**しています。対案は軸対称の壁面積 `2πr ds` で重み付けすること。0.8 は設計余裕として妥当ですが、CONTUR の摩擦比による予測であり保証ではありません。入口の `0.05 r_t` 除外（[`cold_pair_mesh.py:47`](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cold_pair_mesh.py:47)）は設計用の外挿に限定し、結果では入口帯も別集計して残してください。 |
| **Major** | **「同じ FP64 なら比への影響は消える」：却下** | 比の数値誤差は概ね冷却側と断熱側の相対誤差の差です。同じバイナリでも、異なる熱境界条件で相殺する保証はありません。さらに現行判定の `u` は時間変動だけです（[`plan:218`](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-isothermal-wall-chain.md:218)）。対案は、この対をまず**当該格子・FP64 条件での比較**として扱い、時間変動を全数値誤差と呼ばないことです。 |
| **Major** | **「生産格子は約40 ULPなので影響なし」：却下** | 生産格子を再生成し、FP32 化した座標が保存済み `_band_ab/prod_confirm/prep/nozzle.h5` の **194000節点と完全一致**することを確認しました。丸め前に対する第一層厚の誤差は **−2.094〜+1.984%**。これは δ_E/δ_C の誤差そのものではありませんが、無影響の証明にもなりません。対案は「幾何量子化あり、流れの評価量への影響は未評価」と記録することです。 |

V-c45 の事前登録は、次のように書き換えるのが妥当です。

- 生産バイナリとの sha256 一致を、**固定した FP64 用ソース commit・全差分・ビルド設定・forge と実変換器それぞれの sha256** の照合に置き換える。「typedef 4 行だけ」とは記録しない。両 run の格子・物性DB・実効設定を照合し、熱境界条件以外を固定する。
- 幾何検査は HDF5 の dtype だけで済ませず、生成時の double 座標に対する第一層厚の保持を確認する。同じ変換器で再変換して一致する検査は、**同じ丸めを再現しても通る**ため補助扱いにする。
- FP64 冷却格子の断熱結果と生産 FP32 の差は、**「格子と精度を合わせて変更した感度」**と記録する。「格子感度」だけとは呼ばない。同じ生産格子の FP32/FP64 比較は精度変更の影響を調べる資料にはなりますが、冷却／断熱比の誤差相殺は証明しません。
- 静定は R_NS だけでなく、**分子・分母それぞれの δ_E** も確認する。同じ割合でドリフトすれば比だけは静定します。判定する全 x 点と判定窓を揃え、`check_quasisteady.py` は5枚全体を使う設定を明記する。`NOT CONVERGED (stalled/plateau)` を許容する場合も、その VERDICT を残して「収束」と表現しない。

結論: **まず座標読み込みだけを `stof`／`stod` で変える FP64 変換器の CFD 0 step A/B を行い、入力精度の保持を確認してください。**

第 1 仮説: **今回の座標精度不足は、typedef の変更だけでは解消せず、`stof` による入力時の FP32 丸めが残る。** 確度: **高**

  根拠: 上記コードに加え、保存された生産壁と現行の格子表から、**4719 × 121節点**をメモリ上で再生成しました。測定演算を double に固定した結果です。

  | 座標の経路 | AR 最大 | AR > 5000 | 第一層厚の最大相対誤差 |
  |---|---:|---:|---:|
  | 生成時の double 座標 | 4572.89 | 0 | 基準 |
  | 現行 `.10g` → double 読み込み | 4573.62 | 0 | 0.2245% |
  | 現行 `.10g` → float 読み込み → double | 5836.57 | 327 | 24.609% |

  実例では、入力半径 `0.1265365707` と `0.12653663` の差は double 読み込みで **5.9300e−8 m**、float 読み込みで **4.47035e−8 m** になります。**潰れた第一層は0本**でした。「第一セルが表せない」より「第一層厚の相対精度が不足する」が正確です。これは再生成座標の検査であり、未取得の FP64 変換後 HDF5 の品質 VERDICT ではありません。

  反証条件: 対象実変換器が `stof` 経路を通らず、typedef 変更だけで生成時の第一層厚を相対 `1e−6` 以内に保持すること。

第 2 仮説: **座標入力を直した後も、非直交部の格子誤差が熱閉包の数%差の判別に残る。** 確度: 中。非直交性は確認済み、δ_E への影響は未確認。

第 3 仮説: 追加しません。

判別 A/B: **両腕とも同じ FP64 typedef・同じ `.17g` の msh・同じ変換設定とし、座標読み込みだけを A=`stof`、B=`stod` に変える。forge は0 step。** 生成時座標を基準に、第一層厚、非正の層厚の有無、AR、HDF5 の実データ精度を比較します。

- **A が不合格、B が第一層厚の相対誤差 ≤ `1e−6`**なら、第1仮説を支持し、座標入力の修正を採用。
- **B でも超過する**なら、`stof` 修正だけで解消するという仮説を棄却し、残る切り詰め箇所を追う。冷却 NS は投入しない。

やらない方がよいこと: **FP32 で作った HDF5 を double にキャストして使うこと、sha256 制約だけ外して現状の (A) を投入すること。** (B) は必要な壁解像を犠牲にするので選びません。(C) も `stof` 問題を共有し、混合精度経路が未検証なので今回の推奨にはしません。入力経路を直した FP64 版を選びます。

呼び出し側の前提への異議: **AR 超過をすべて FP32 のせいにはできません。** 初版には壁勾配の見落としがあり、現在はその補正が入っています。また、ブリーフの4496節点列に対し、手元の現行表は4719列です。補正版で得た数値を、削除済みの初版 run の実測値と混同してはいけません。

不足情報: 対象の FP64 ビルド、削除済み冷却格子の HDF5、生産 `run_0167`／`run_0179` の元場・全残差・判定ファイルは手元にありません。したがって、ビルドの成立、元 run の収束、冷却比の妥当性は今回認定していません。

ファイル変更・forge 起動はしていません。**plan 未反映**です。反映先は `tooling-nozzle-isothermal-wall-chain.md` §5.1 #13・#14、§6 V-c45 です。
