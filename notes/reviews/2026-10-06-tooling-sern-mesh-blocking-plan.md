# codex レビュー: tooling-sern-mesh-blocking (plan)

- **plan**: [`plans/active/tooling-sern-mesh-blocking.md`](../../plans/active/tooling-sern-mesh-blocking.md)
- **stage**: `plan`
- **date**: 2026-10-06
- **commit**: `b29c6ecc` (feature/sern-design)
- **codex**: effort `high`, 3.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M4/m1

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
全ヘキサ方式と「B4a–B4c はメッシュ検証、B4d＋B5 後に CFD」という順序を支持します。  
ただし、壁層検査には誤合格を再現できる欠落があり、形状・細分・restart の契約も実装前に補う必要があります。

依頼ファイル全文と現行コードを照合しました。依頼文の「上流板厚 0.02 H＋テーパ」は、現行 plan では「生産板厚 0.005 H 一定」に訂正されています。この訂正済み事項は未解決指摘に数えていません。

1. **Major — 壁層検査が、壁上の節点を「第一内部点」として合格させる。**

   現行検査は対向節点を選びますが、その節点が別の壁上にあるかを除外していません。層数も生成設定を記録するだけで、実際の接続から検証していません。根拠: [対向節点の選択](/home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py:493)、[距離・合否判定](/home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py:507)。

   読み取り実行で、**全６面が壁・内部節点ゼロの単一立方体**を既存 `check()` に渡すと、６壁すべてが `ratio_min=ratio_max=1.0, ok=True` になりました。これは B0c が要求する「選択点が内部点」と矛盾します。

   **対案:** B0c の残作業を「間隔比だけ」から改め、内部点の資格、期待する壁タグ集合、壁節点の全数被覆、実層数も検査対象に戻してください。稜では隣接する複数壁を考慮して内部点を選び、取得不能なら不合格にします。この反例を拒否できることを B4a の前提にしてください。

2. **Major — 第一層と成長率の受入仕様が、採用済み設定・実装と一致しない。**

   [§6](/home/sano/work/forge-sern-design/plans/active/tooling-sern-mesh-blocking.md:837) は端面も共通の `h₁` に対して判定しますが、§4.12 は一般壁 16 µm・端面 250 µm を採用し、コードも別の `h1e` を使います。両者の比は **15.625** です。B4a の「各端面も同じ指定」には、その指定値がありません。

   また [細分処理](/home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py:547) は `G=1.2^(1/scale)` とするため、`scale=0.5/0.71` では **1.44/1.2928**。さらに [生成側の合格条件](/home/sano/work/forge-sern-design/case/46.sern_design/cad/hex_junction_model.py:399) は `G×1.02` まで許し、§6 の実測比 ≤1.2 と一致しません。

   **対案:** B4a の共通入力を数値表で固定してください。壁タグ別の第一内部点距離、層数の決定規則、後流区間、間隔上限、成長率を明記し、端面も含めた３水準の細分規則を定義します。全水準で実測比 ≤1.2 を満たす系列へ変更し、旧模型の PASS は新仕様へ持ち越さないでください。

3. **Major — B4a の合格条件では、目的の形状を作ったことを保証できない。**

   [§6.2](/home/sano/work/forge-sern-design/plans/active/tooling-sern-mesh-blocking.md:853) は §6 の 1・3・4・5 群を要求しますが、MOC 輪郭への誤差、法線板厚、端面位置の定量検査がありません。誤った輪郭でも品質・閉性・壁距離は合格できます。解像度間で同じ形状という条件も、参照形状への正確さを保証しません。

   特に [法線オフセットと `x=L_cowl` の平端面](/home/sano/work/forge-sern-design/plans/active/tooling-sern-mesh-blocking.md:762) の接続規則が未定義です。生産カウル角 5°・板厚 0.5 mm では、内壁端点をそのまま法線オフセットすると外壁端点の x が **約43.6 µm** 動きます。

   **対案:** B4a(1) の接続表に、オフセット曲線の延長・切断方法、端面との交点、上流端の閉じ方を追加してください。`ramp_fillet` 適用後の参照輪郭、法線板厚、物理端点に対する誤差許容を事前登録し、B4a の独立した形状ゲートにします。

4. **Major — 新メッシュの双子節点検査だけでは、旧形状からの restart を保証できない。**

   [B4d](/home/sano/work/forge-sern-design/plans/active/tooling-sern-mesh-blocking.md:822) は `interp_field` と双子節点不在の検査を予定しています。しかし旧メッシャは [側壁内外の双子節点](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:268) を生成し、[同一座標を与えます](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:310)。移植元に曖昧さが残ります。

   [`interp_field.py`](/home/sano/work/forge-sern-design/solver_density_cuda/tools/interp_field.py:212) は領域を区別せず最近傍を１点選びます。同一座標に排気側・外気側の２状態を置いた小例では、両側の移植先が同じ状態を選択しました。化学種署名の一致では、この空間的な取り違えを検出できません。

   **対案:** SRC/DST 両方の領域・壁側情報で候補を制限する restart 契約を追加してください。新たに流体となる領域の初期化も別途定義し、同一座標に異なる組成を持つ試験で取り違えを拒否することを B4d の受入条件にします。

5. **Minor — 現行方針と履歴が混在し、作業範囲を誤読させる。**

   表題は Salome 先行のまま、[§7](/home/sano/work/forge-sern-design/plans/active/tooling-sern-mesh-blocking.md:871) は runner を変更しないとしていますが、B4d/B5 は変更を要求します。B1d も旧診断の状態を残し、後続の accepted plan と R1 に訂正が分散しています。

   **対案:** 冒頭に現行の方式・対象フェーズ・優先する仕様を集約し、§2・§7・§8 と残作業表を同期してください。履歴は残し、失効・訂正先を明示します。

**推奨は、全ヘキサ方式を維持し、上記 Major を plan に反映してから B4a に進むことです。** 優先順は、壁層ゲートの修復 → 数値仕様の固定 → 形状ゲートの追加 → restart 契約の修正です。B1b の修正や accepted 済みの SLAU 対策を再開する必要はありません。B6 では準定常性と残差収束を別々に判定する現在の方針を維持してください。

ファイル変更・CFD・実メッシュ生成は行っていません。検証はコード照合と書込みを伴わない再現例までです。**plan 未反映**です。

指摘数: Critical 0 / Major 4 / Minor 1
