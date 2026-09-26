# codex 諮問 (diagnose): solver-owned-species-db

- **brief**: [`notes/reviews/briefs/2026-09-27-solver-owned-species-db.md`](../../notes/reviews/briefs/2026-09-27-solver-owned-species-db.md)
- **date**: 2026-09-27
- **commit**: `7c1ca0ce` (feature/gap-heating-precision)
- **codex**: effort `xhigh`, 4.6 min, rc=0
- **結論**: **移行する――推奨は (b) を温度区間の和集合に拡張した「ソルバ起動時の係数合成」とし、最初の変更は内蔵種を含む解決済み物性の保存・内容照合を導入すること一つに絞る。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（コード確認と、ファイルを書かない関数単体の検算に基づく）

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| H1：熱物性の正本をソルバ側へ集約 | **採用 / Major** | ソルバと設計側は別々に内蔵表を構築している（[speciesDB.cpp:83](/home/sano/work/forge/solver_density_cuda/input/speciesDB.cpp:83)、[composition.py:76](/home/sano/work/forge/design/forge_design/gas/composition.py:76)）。**版を固定した CEA 由来の共通データをソルバ配布物に置き、C++・Python が共有する**。外部 DB 上書きは残す。ユーザ入力と、再現用に保存する解決済み物性を分ける。 |
| H2：温度区間が違うなら実行時の構成種評価が必要 | **却下 / Major** | NASA-9 は係数に線形（[thermo_d.cuh:39](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:39)）。**構成種の全区切りの和集合で区間を分ければ、起動時に区間ごとの係数を合成できる**。端の定 cp 外挿も定数 cp の係数と積分定数で表現できる。ただし現行の `low/high` 固定構造は拡張が必要。「区間不一致＝合成不能」は現実装の制約である。 |
| H2：実行時の重み付き和で LJ 平均問題も解消 | **却下 / Major** | `cp/h` の線形混合と Wilke 混合は別物（[composition.py:275](/home/sano/work/forge/design/forge_design/gas/composition.py:275)、[thermo_d.cuh:406](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:406)）。対案は**lump の構成情報を保持し、輸送物性では実種へ展開して全実種の混合則を評価する**こと。lump 内で求めた粘性をさらに混ぜる方法も同等とは限らない。 |
| 観測8：既存の署名で内蔵 DB 更新も保護される | **却下 / Major** | [forge_species.py:189](/home/sano/work/forge/solver_density_cuda/tools/forge_species.py:189) は内蔵種の係数・区切りを `None` にし、[同:221](/home/sano/work/forge/solver_density_cuda/tools/forge_species.py:221) は両側 `builtin` の比較を省略する。対案は**使用時点の解決済み内容を保存し、内容で照合する**こと。現在の内蔵表から過去 run の署名を再生成してはいけない。 |
| H3：`atoms` を run の手入力から外す | **採用 / Minor** | ソルバの DB 読込は `atoms` を使用しない（[speciesDB.cpp:173](/home/sano/work/forge/solver_density_cuda/input/speciesDB.cpp:173)）。ただし元素診断には使う（[composition.py:248](/home/sano/work/forge/design/forge_design/gas/composition.py:248)）。**入力から外し、共通 DB と解決済みメタデータには残す**。非整数は平均原子数であり、それ自体は誤りではない。 |

結論: **移行する――推奨は (b) を温度区間の和集合に拡張した「ソルバ起動時の係数合成」とし、最初の変更は内蔵種を含む解決済み物性の保存・内容照合を導入すること一つに絞る。**

第 1 仮説: 内蔵 DB への移行で先に解消すべき障害は、使用した物性を既存の署名が固定できないことである。 確度: **高**

  根拠: 上記 `forge_species.py:189,221`。比較関数そのものを抽出して単独実行し、名前・MW・datum を同一にした模擬署名の `nasa9_low[2]` だけを `0→1` に変更した結果は以下だった。

  - `source="builtin"`：不一致リスト `[]`
  - `source="file"`：`['N2.nasa9_low[2] 0.0 vs 1.0']`

  この関数は [interp_field.py:67](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:67) の実際の拒否判定にも使われる。**「署名検査がある」は確認できたが、「内蔵係数の変更を検出する」は成立しない。**

  反証条件: 過去 run 作成時の内蔵係数・区間を保持する別の必須照合が存在し、今回の変更を転送前に拒否することが示されれば、移行を妨げるという評価を撤回する。

第 2 仮説: 固定内部組成の気相 lump は、区間を細分して起動時に合成すれば、構成種ごとの `cp/h` 評価と同じ熱力学を表せる。 確度: **高〔代数上〕、実装・GPU 性能は未確認**。根拠は `composition.py:268–272` のモル重み合成と NASA-9 の線形性。凝縮種は独立種として残す。

第 3 仮説: 現行 lump と実種展開の粘性差は、DB の置き場所を変更しても残る。 確度: **高**。`case/44.vitiated_air_wt/run_0510_va3_M4.19_Lc8_noneq_lumpX/species_db.yaml` の丸め済み構成比コメントと閲覧した式による float64 検算では、平均 LJ の粘性は実種 Wilke に対して **200 K で −1.30%、300 K で −0.98%**。これは物性式の検算であり、CFD の実測値ではない。

判別 A/B: **最初の変更に対する host の照合試験を一つ行う。** 作成時の物性記録を固定し、現在の内蔵 DB の `N2.nasa9_low[2]` だけを、A は同値、B は `+0.001` とする。名前・MW・DB 名・datum は固定。長さは各1回、CFD は **0 step**。見る量は照合の許可／拒否と差分理由。  
→ **A のみ許可・B は該当係数を示して拒否**なら、内蔵更新を検出できないという障害は解消。**B も許可**なら障害が残り、生成 `species_db.yaml` を廃止する段階へ進めない。

やらない方がよいこと: `species_db.yaml` を先に消すこと、内蔵化と同時に CEA の版・MW・外挿規約まで変えること、`full` を強制すること。`full/lumped` の選択維持は [accepted plan:155](/home/sano/work/forge/plans/accepted/thermophysics-cea-mole-fraction-species.md:155) の既存判断であり、今回これを覆す証拠はない。

呼び出し側の前提への異議:

- **Major：H2O の「液比熱は液相 CEA の傾き」は無条件には成立しない。** `h2o_latent` は独自の気相係数と外挿を使う（[condensationProperties_d.cuh:239](/home/sano/work/forge/solver_density_cuda/cuda_forge/condensationProperties_d.cuh:239)）。種 DB の気相エンタルピーと異なれば、EOS が暗黙に使う液相エンタルピーにも差が入る。対案は共通の気相評価と整合した潜熱評価。ただしこれは外挿モデルの変更も伴い得るため、DB の移設とは分けて検証する。
- **Minor：「s° は流れに不使用」は誤り。** [thermo_d.cuh:801](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:801) は等エントロピー変換に使う。固定組成では混合エントロピー項が差分で消えることが省略の根拠であり、移行時も `s°` の区間・外挿規約を保持する。
- **Major：`full` の +23% は候補 (c) の費用を示さない。** 輸送方程式数も変わる比較だからである。また `run_0509–0511` の `NOT CONVERGED / ALL STEADY` は今回はブリーフの申告のみで、収束解の一致を判断する証拠には採用していない。対案は物性・署名の決定的な照合を先行し、CFD の移行判定は別途 VERDICT と判定区間を揃えて行う。

不足情報: 採用する CEA データの版・ハッシュ、移行回帰の具体的許容差、対象 run の VERDICT と判定区間。指定された閲覧制限に従い、残差・場・ログは読んでいない。**ファイル変更なし。plan 未反映であり、上記は呼び出し側が検証・記録するための設計提案である。**
