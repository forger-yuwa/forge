# codex レビュー: tooling-nozzle-sern-3d (plan)

- **plan**: [`plans/active/tooling-nozzle-sern-3d.md`](../../plans/active/tooling-nozzle-sern-3d.md)
- **stage**: `plan`
- **date**: 2026-09-21
- **commit**: `b7ad1dc4` (feature/sern-design)
- **codex**: effort `high`, 7.1 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M7/m1
- **focus**: §4.44 R5q の設計だけを見てほしい: 側壁に物理厚み + z 分布の station 依存化 + カウル側端テーパ廃止、の 3 点同時で厚さ0の交線を消す案。トポロジ不変・dup2 の条件も k_sw 索引も変えない前提が成立するか、検証項目 5 件で足りるか、見落とした接続 (入口面・後縁 i_sw・ext_top ブロック・面分類の zm) がないか。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

§4.44 の目的は妥当ですが、「座標だけ変更し、`dup2` とトポロジを維持する」設計では接続を保証できません。  
特に `L_sw < L_cowl` の露出カウル側端、側壁下端、機体側面との接合が未定義です。検証5件も補強が必要です。

現行の `python3 -B design/tests/run_sern_mesh3d_tests.py` は **ALL PASS**。以下はコード読解と、提案の変更をメモリ内だけで再現した幾何検査によります。ファイル変更・CFD実行はしていません。

1. **Major — 側壁後縁より下流の、有限厚カウル側端を閉じる面がない。**

   生産設定は `L_sw=0.8`、`L_cowl=1.2` です。[設定ファイル](/home/sano/work/forge/case/46.sern_design/problem_3d_prod_m6on_wallres.yaml:106)

   `dup1` は `i<i_te` まで存在する一方、側壁境界は `i+1<=i_sw` までです。`sz` をステップにすると、その間ではカウル上下端が別座標になるのに、側端面が追加されません。[ノード生成](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:235)、[境界生成](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:334)

   縮小検査例（`L_cowl=1`、`L_sw=0.8`、カウル厚 `0.005H`）では、未分類の単独所有面が **12面**ありました。このID上の不整合は旧版にもありますが、旧版でゼロだった側端の上下間隔が、提案版では `x/H=0.8324325` で **0.00418919H** に開きます。既存欠陥をそのまま引き継げません。

   **対案:** `L_sw` から `L_cowl` までの露出カウル側端を物理壁として定義し、外部流セルの接続、境界タグ、壁距離、力への算入まで設計する。「トポロジ不変」は撤回する。

2. **Major — `dup2` の `j>jm` 条件を維持すると、側壁下端に指定厚みを与えられない。**

   `node(..., j=jm, side="up_out")` は `dup2` を使わず、`base` を返します。その座標は提案では内側の `z=W/2−t_sw/2`、高さは `cowl_out` の `ym−t_c/2` です。外側の単一中間線 `(ym, W/2+t_sw/2)` にはなりません。[根拠](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:263)、[ノード選択](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:275)

   `W=2`、`t_sw=0.02`、`t_c=0.005` の入口断面で、提案操作の実座標は次でした（単位 `H`）。

   | 節点 | y | z |
   |---|---:|---:|
   | 側壁内面下端 | 0.0025 | 0.99 |
   | 側壁外面下端 | −0.0025 | **0.99** |
   | 外面の次のj節点 | 0.004 | 1.01 |

   外面は最初のj区間で斜めに立ち上がり、接合形状が第一層厚に依存します。

   **対案:** カウル上下端・側壁内外面・カウル下外部流の接合断面を先に固定する。必要な独立節点と接続面を追加し、`j=jm` を単純な座標移動で済ませない。

3. **Major — `ext_top` と機体側面バンドが「変更4箇所」から漏れている。**

   `_add_ext_top3d` は新節点を `zs[k]` で生成し、`_add_vehicle_side3d` も内部節点を `zs[k]` に置きます。一方、後者の下端は移動済みの `dup2` を共有します。[上部ブロック](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:433)、[側面バンド](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:501)

   したがって4箇所だけ変更すると、側壁外面上端の `W/2+t_sw/2` から、側面バンド内部の旧 `W/2` へ、格子のj区間で戻る形になります。IDによる閉性が通っても、機体接合形状は格子依存です。また、共通の `zsx[i][k_sw]` を渡すだけでは、内面用と外面用の異なる座標を表現できません。

   **対案:** ランプ端、側壁外面上端、`vehicle_side`、`vehicle_top` の接合輪郭を定義し、ブロック別の座標写像と共有IDを設計する。

4. **Major — `i_sw` は物理的な `L_sw` のstationではない。**

   現行コードは `i_sw=argmin(abs(xs−L_sw))` で、`L_sw` をstation列へ挿入していません。[根拠](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:137)

   上記縮小例では、`L_sw=0.8` に対して共有開始位置が **0.7734683014H**。指定した厚み則はそこで **0.0033164623H** を要求するのに、`dup2` は既に存在しません。共有IDと非ゼロ厚みの要求が衝突します。

   **対案:** `0.8L_sw` と `L_sw` を物理stationとして挿入してから索引を作る。厚みゼロとノード共有開始を同じstationで成立させ、`ni_noz` を変えても後縁位置が動かないことを検証する。

5. **Major — 「既定0でビット一致」と「`sz` 廃止」が矛盾する。**

   現行の `sz` はカウル厚が正なら側端でゼロになります。[根拠](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:244)

   側壁厚をゼロにしても、`sz` をステップへ替えるとカウル側端のy座標が変わります。`L_sw=L_cowl=1`、カウル厚 `0.005H` の縮小例では、**1,108節点**の座標が変化しました。§4.44 の検証1は現設計では通りません。[要件](/home/sano/work/forge/plans/active/tooling-nozzle-sern-3d.md:1812)

   **対案:** `sidewall_thickness==0` は旧座標生成を保持する互換分岐にする。正値の場合だけ、接続を完成させた新形状を使う。

6. **Major — 入口幅が変わるのに、入力経路と推力の帳簿が変更対象に入っていない。**

   `runner_sern3d.prepare` はパラメータを明示列挙して構築するので、新キーをメッシャへ追加するだけではYAMLから届きません。[根拠](/home/sano/work/forge/design/forge_design/evaluate/runner_sern3d.py:116)

   また、`half_W_m=W H/2` を保存し、`forces3d` はその幅で入口運動量・圧力項と理想推力を計算します。[幅の保存](/home/sano/work/forge/design/forge_design/evaluate/runner_sern3d.py:166)、[入口項](/home/sano/work/forge/design/forge_design/evaluate/runner_sern3d.py:287)

   新形状では入口流路幅が縮むため、旧幅による入口項は実際の入口面と整合しません。格子列を取り直しても、この系統誤差は消えません。

   **対案:** YAML→生成→`prepare_info` を変更範囲に追加する。入口項は実際の `inlet_nozzle` 面積・流束に基づかせ、係数の基準幅と理想推力の定義を明記する。基準幅を旧値に固定する場合でも、実入口項とは分離する。

7. **Major — 既存幾何試験を「そのまま掛ける」だけでは、新接続を保証できない。**

   現行試験の該当ケースはカウル厚ゼロ、`L_sw` 未指定です。[試験条件](/home/sano/work/forge/design/tests/run_sern_mesh3d_tests.py:73)  
   Jacobian検査は各hexの一頂点からの三辺だけで、座標一致数の検査も「全 `dup1+dup2` が一致する」前提です。[検査実装](/home/sano/work/forge/design/tests/run_sern_mesh3d_tests.py:125)

   **対案:** 少なくとも以下を新モードの必須試験にする。

   - 有限カウル厚×有限側壁厚で、`L_sw<L_cowl`／同位置／station間指定を検査する。
   - 入口・側壁下端・側壁後縁・機体接合部で、設計断面、面所有、境界タグ、入口BCの節点重複を照合する。
   - x・y・zの解像度を個別に変え、物理輪郭と面積が不変であることを確認する。
   - float32化後もhexの全頂点・内部評価点のJacobian、変換後の双対体積・閉性を確認する。
   - 厚みの有効範囲を検査し、内外スパンが非正になる入力を拒否する。

   CFD受理には既存の準定常性・床・壁解像ゲートを引き継ぐ。幾何修正だけで密度床の原因が解消すると確定してはいけません。

8. **Minor — `zm` の分類条件についての説明がコードと違う。**

   `ramp` 判定は `k<k_sw` です。`zm` が使われるのは、`vs_on` が偽の経路での **`W_vehicle/2`** との比較です。[根拠](/home/sano/work/forge/design/forge_design/meshing/mesh_sern3d.py:321)  
   「`W/2` 比較だから分類不変」は根拠になりません。またstation依存化後の面重心は、一方のstationだけから計算できません。

   **対案:** 索引による分類はそのまま明記し、物理座標による分類が必要な面は、実際の4頂点から重心を計算する。`ext_top=False`／`vehicle_side=False` も回帰対象にする。

**推奨は、有限厚化を維持しつつ「トポロジ不変」の制約を外し、接合部を設計し直すことです。** 優先順は、露出カウル側端と側壁下端の接続 → 機体ブロックとの接合 → 物理station固定 → 互換分岐・入口帳簿 → 幾何検証です。接続の小規模検証が通ってから生産格子列へ進めてください。

本レビューの提案は **plan未反映**です。依頼の読み取り専用条件に従い、採否・更新は行っていません。

指摘数: Critical 0 / Major 7 / Minor 1
