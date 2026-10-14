# codex レビュー: tooling-nozzle-wall-temperature-profile (plan)

- **plan**: [`plans/active/tooling-nozzle-wall-temperature-profile.md`](../../plans/active/tooling-nozzle-wall-temperature-profile.md)
- **stage**: `plan`
- **date**: 2026-10-13
- **commit**: `138e83b2` (feature/nozzle-deliverables)
- **codex**: effort `high`, 4.1 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m2
- **focus**: §4 の配管 (単一ソース・実効 r_t・CSV・カード・成果物) と §6 の合格条件が自己一致になっていないか。methods は読まなくてよい

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

課題の同定と、実効 `r_t` を使う設計は妥当です。既存機能の重複実装ではなく、設計チェーンへの接続が未完です。  
ただし、Fluent 出力には補間方式の不整合があり、範囲外補間・成果物の出典・格子見積もりを現行の合格条件では取り逃がします。

1. **Major — Fluent の壁ノード profile は、指定した線形分布を境界面に再現しない**

   **根拠:** [plan:83](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wall-temperature-profile.md:83) は壁ノードの `point profile`、[plan:115](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wall-temperature-profile.md:115) はその読み戻しを検査します。しかし Fluent の既定は、境界面に最近傍の profile 値を与える方式です。ファイルの往復一致では、この適用誤差を検出できません。[Ansys 公式仕様](https://ansyshelp.ansys.com/public/views/secured/corp/v252/en/flu_ug/flu_ug_bcs_sec_prof.html)

   実測として、`case/70.m5_air_wt_demo/run_0006_ns_pass0_iface_60k/nozzle.msh` 上で、入口320 K→設計スロート500 K→出口330 Kを与え、既定補間を再現しました。**1,599壁面中362面で線形分布との差が0.1 K超、最大0.236215 K**です。これは既存メッシュ上の補間評価であり、Fluent 実行結果ではありません。

   **対案:** 同梱する Fluent メッシュの**壁面中心**で原表を評価し、そこで `point profile` を生成する。V3 はファイル往復に加え、境界面への適用値を独立に計算した原表の値と照合する。対象メッシュを固定し、Fluent 実機での読み込み未確認は残す。

2. **Major — 範囲チェックが遅く、V1b も範囲外クランプを合格にできる**

   **根拠:** [plan:59](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wall-temperature-profile.md:59) の検査場所は `prepare_ns` ですが、カード生成はそれ以前に [wt_card.py:1235](/home/sano/work/forge-deliv/design/forge_design/intake/wt_card.py:1235) で `solve_rt_throat` を実行します。その反復中も積分法が呼ばれます。

   [deltastar_integral.py:138](/home/sano/work/forge-deliv/design/forge_design/feedback/deltastar_integral.py:138) の `np.interp` は範囲外をクランプします。既存関数で確認すると、表の終端1.2 mに対して1.44 mを評価しても330 Kを返し、原表側も `np.interp` で評価するV1b型の比較は**差0 K**になります。「同じ温度」と「定義域内」を混同しています。

   **対案:** 各積分呼び出しの前に、実効 `scale` で換算した積分区間の包含を検査する。本計画では試行半径も覆う表を要求し、不足は必要範囲付きで拒否する。最終物理壁の検査も別途維持する。V1b に「Bだけ範囲外になる場合は例外」の負例を追加する。

3. **Major — SHA記録だけでは、NSと成果物の単一ソースが保証されない**

   **根拠:** [plan:70](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wall-temperature-profile.md:70) はCSVのSHAを記録しますが、起動・再開・成果物生成時に照合する契約がありません。既存の成果物コードには複数の出典があります。

   - [fluent_msh.py:693](/home/sano/work/forge-deliv/design/forge_design/export/fluent_msh.py:693)：外から渡した問題YAMLを `prepare_info` より優先。
   - [deliverables.py:339](/home/sano/work/forge-deliv/design/forge_design/export/deliverables.py:339)：条件表は `bcondConfig` のスカラー値から作成。
   - [stage_manifest.py:274](/home/sano/work/forge-deliv/solver_density_cuda/tools/stage_manifest.py:274)：BC識別はYAMLのハッシュで、壁温CSVの内容を含まない。

   同じスロート温度のまま分布だけ変わると、代表 `Ts` とBC識別は変わりません。V3の書出し・読戻し一致だけでは、別の表を忠実に輸出した誤りを検出できません。

   **対案:** run内に固定した正規化表・CSV・SHAを実行条件の正本とし、起動・再開・exportで照合する。不一致は停止し、外部YAMLや代表 `Tw` への代替を禁止する。CSV変更、別問題YAML、同じ代表温度で異なる分布を負例として登録する。

4. **Major — 局所壁温を使う格子見積もりに、対応する合格条件がない**

   **根拠:** [plan:77](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wall-temperature-profile.md:77) は局所温度を要求しますが、V1〜V4には `estimate_mesh` の分布対応を数値検査する条件がありません。既存実装は [wt_card.py:665](/home/sano/work/forge-deliv/design/forge_design/intake/wt_card.py:665) で一様温度を作り、`:670`以降で参照温度・壁密度・粘性・第一層厚へ伝えます。

   分布を代表温度に置き換える実装でも、NS・積分法・Fluentだけ正しく配管すれば登録試験を通せます。

   **対案:** 同じスロート温度で上流・下流温度が異なる二つの表を使い、選定断面の `Tw` と丸め前の第一層厚を、独立に計算した期待値と照合する。合格条件は例えば相対誤差≤1e−8。代表温度固定・`x` の無次元値をそのまま使う誤実装が落ちることも確認する。最終的な格子品質・壁解像の判定は別に残す。

5. **Minor — 入力契約に有限値・配列長・CSV解決規則が不足している**

   **根拠:** [plan:58](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wall-temperature-profile.md:58) の「狭義単調増加・正」だけでは、例えば正の無限大を排除できません。[plan:75](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wall-temperature-profile.md:75) も、表とCSVの併記、相対パスの基準を定めていません。

   **対案:** `x_m`・`Tw_K`は同長の一次元有限数値列、`scale`は正の有限値とする。`Tw`・表・CSVは厳密に一つ。CSV相対パスはカード所在基準とし、NaN/Inf・長さ不一致・併記・異なる作業ディレクトリからの読込みを試験する。

6. **Minor — V2は配管検査として妥当だが、起動条件と根拠の扱いを明確にすべき**

   **根拠:** [nodeWallDirichlet_d.cu:84](/home/sano/work/forge-deliv/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:84) は `T[ic] = Tw` を直接実行します。実際、`case/70.m5_air_wt_demo/run_0006_ns_pass0_iface_60k/` の最終場では**壁1,600ノードすべて350 K、指定値との差0 K**でした。一方、今回再実行した収束判定は：

   `NOT CONVERGED (still converging — run more steps)`

   したがって壁温一致は配管検査には使えますが、流れの妥当性を保証しません。計画が収束を対象外としている点は正しいです。また、0.068 Kの根拠である `case/48.flat_plate_cooled_m4/run_0015_wallprofile/` はこの作業コピーになく、数値を再検証できませんでした。

   **対案:** V2を「温度拘束の配管検査」と明記して維持する。warm start元、同一／別メッシュの引継ぎ方法、壁温変更後のsoft段、メッシュ品質、序盤と最終の有限性検査を登録する。[既存等温壁計画:116](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-isothermal-wall-chain.md:116) に対応手順があります。収束VERDICTは合否と分離して必ず記録する。

**推奨は、既存ソルバを再利用する現在の方針を維持し、上記1→2→3→4→5→6の順で計画と合格条件を修正してから実装することです。** 実効 `r_t` のV1bは有効な検査として残してください。今回の指摘は新しい数値スキームの必要性ではなく、入力から適用先までの契約と検証の不足です。ファイル変更禁止に従い、**plan未反映**です。

指摘数: Critical 0 / Major 4 / Minor 2
