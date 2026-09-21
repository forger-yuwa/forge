# codex レビュー: tooling-sern-mesh-blocking (plan)

- **plan**: [`plans/active/tooling-sern-mesh-blocking.md`](../../plans/active/tooling-sern-mesh-blocking.md)
- **stage**: `plan`
- **date**: 2026-09-21
- **commit**: `186ec51d` (feature/sern-design)
- **codex**: effort `high`, 6.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M7/m1
- **focus**: §4.1–§4.4 の形状モデルとブロッキングが全ヘキサで閉じるか (特に: バタフライ C リングと H ブロックの共有面整合、フィレット r_f→0 の絞りで位相が保てるか、スロット後流/カウル後流ブロックの接続、L_sw<x<L_cowl の露出カウル側端)。gmsh geo カーネルの transfinite 体積で x 変化断面が実現できるか。§6 の受入条件で方式確定として十分か。

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
有限厚形状と gmsh 全ヘキサ化の方向は妥当ですが、現案には CAD エンティティの退化、後流の共有面分割、受入判定の欠陥があります。  
まず接続設計と検査仕様を修正すべきです。方式そのものを撤回する必要はありません。

`plans/README.md` と `plans/accepted/` を確認し、このメッシャを実現済みの計画は見つかりませんでした。親計画から独立した実証フェーズにするスコープも妥当です。

以下はコード読解と、gmsh 4.15.2・既存品質関数を使ったメモリ上の試験に基づきます。ファイル変更・CFD 実行はありません。根拠とされる `run_0421`・`run_0422` はローカルに存在せず、その実測値と収束判定は再検証できていません。

1. **Major — 「円弧＋直線を1辺とする」説明では、`r_f→0` を含む gmsh の固定位相は成立しません。**

   **根拠:** [plan:71](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:71) は論理的な辺長が正であることを根拠にしていますが、[plan:119](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:119) の実装では直線と円弧は別エンティティです。`geo.addSurfaceFilling` は実際の境界曲線数が3または4である必要があります。[Gmsh公式仕様](https://gmsh.info/doc/texinfo/#index-gmsh_002fmodel_002fgeo_002faddSurfaceFilling)

   メモリ上で次を再現しました。

   - 5曲線のループ: `5 borders instead of 3 or 4`
   - 半径ゼロの円弧: `Zero radius in circle or ellipse`

   断面を `addPlaneSurface` に替えても、円弧を別エンティティとして残す側面パッチの退化は解消しません。case/49 は多辺の平面に角を指定して**押し出す**実装であり、今回の接続の実証ではありません（[build_hex_mesh.py:255](/home/sano/work/forge/case/49.plate_annular_cavity_m5/cad/build_hex_mesh.py:255)）。

   **対案:** B1 より前に、各ブロックの頂点・曲線・6面の対応を列挙してください。フィレット終端は、消滅する円弧をそのままブロック辺にしない専用の3D遷移として設計する必要があります。近似曲線・曲面を使う場合は、U3 の形状誤差と接線誤差の許容値を先に確定してください。

2. **Major — カウル後流で N／スロット後流／U の面分割が一致しません。**

   **根拠:** [plan:76](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:76) の U の z 分割には `W/2` がありません。一方、[plan:94](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:94) のスロット後流は `[W/2, z_out]` を占めます。

   `x>L_cowl` ではカウル後流の上面が、N 側の `[0,W/2]` とスロット側の `[W/2,z_out]` に分かれます。下面の U は `[0,z_out]` の一面です。このままでは、6面の transfinite ブロック同士の一対一接続になりません。「同じ分割数」はこの面分割の不一致を解決しません。

   **対案:** `z=W/2` の分割をカウル後流と U に伝播させ、`L_sw`・`L_cowl` の前後で全共有面を一覧化してください。共有条件は、分割数に加えて**同一 surface ID、節点ID、分布、向き**まで規定すべきです。露出した `cowl_side` とその下流への接続も、この一覧で閉じることを確認してください。

3. **Major — 断面リングの第一層指定だけでは、全壁・端面の ±5% を満たせません。**

   **根拠:** [plan:73](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:73) はリング横断方向を壁法線としていますが、x 方向に傾く壁の3D法線とは異なります。例えば傾斜22°のランプで y 方向間隔を `h1` にすると、法線距離は `h1 cos(22°)=0.92718 h1`。既に −7.28% で受入範囲外です。

   また、`sidewall_end`・`cowl_base` の法線は x 方向ですが、[plan:86](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:86) にはこれらの直後の第一層配置がありません。`r_f=0` の45°接続でも、辺長と壁法線距離は一致しません。

   **対案:** 壁面ごとに3D法線に対する目標距離を定義し、端面下流の x 分布を追加してください。共有辺上では複数壁への距離を個別に検査する設計が必要です。全壁一律の `h1` が不要なら、壁別の目標値を明示して受入条件を修正してください。

4. **Major — 入口面積不変の主張と、格子間 `1e-9` 一致条件が誤っています。**

   **根拠:** [plan:57](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:57) の外向き板厚は幅を維持しますが、入口から存在する上下の内隅フィレットは流路面積を減らします。通常の接円形状なら、半スパン断面の減少量は

   \[
   \Delta A=2(1-\pi/4)r_f^2.
   \]

   `r_f=0.10H` では `0.00429204 H²` です。現行の入口項は矩形面積に基づく `half_W` 倍のままです（[runner_sern3d.py:287](/home/sano/work/forge/design/forge_design/evaluate/runner_sern3d.py:287)）。

   さらに、同じ円弧上に節点を置いても、線形要素の境界は弦になります。半径 `0.10H` の上下四分円を各8・16・32分割した場合、解析面積との差の絶対値はそれぞれ `1.00738e-4`・`2.52208e-5`・`6.30748e-6 H²`。CAD が同一でも離散面積は変わります。[plan:147](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:147) は両者を混同しています。

   **対案:** 「CAD形状・端点の不変」と「離散境界の面積誤差」を別ゲートにしてください。入口解析面積にはフィレットを含め、離散面積は解析値への収束と実用上の許容誤差で判定します。入口帳簿の改修は B5 まで延期できますが、不要になったという記述は撤回すべきです。

5. **Major — `check_wall_resolution.py --geometry-only` は、記載された壁層受入試験を実施しません。**

   **根拠:** [check_wall_resolution.py:186](/home/sano/work/forge/solver_density_cuda/tools/check_wall_resolution.py:186) は評価できた点の統計を出し、指定厚さとの比較をせず `VERDICT: MEASURED` を返します。一部の点が評価不能でも、その壁で1点以上評価できれば失敗になりません。

   また、[同:134](/home/sano/work/forge/solver_density_cuda/tools/check_wall_resolution.py:134) は法線への整列度が最大の隣接点を選びますが、その点が別の壁上にあることを除外していません。今回重点とする隅・端部で、第一**内部**点の保証になりません。

   **対案:** B2 に壁別・点別の自動判定を追加してください。指定値との誤差、評価不能数、選択点の内部所属、必要層数を検査し、未評価点を合格扱いしないことが必要です。`MEASURED` や終了コード0を受入成功にしてはいけません。

6. **Major — float32 の健全性ゲートは、既存ツールの成功だけでは保証できません。**

   **根拠:** [check_mesh_quality.py:79](/home/sano/work/forge/solver_density_cuda/tools/check_mesh_quality.py:79) が調べるのは四面体分割の**体積和**です。局所 Jacobian の正値性ではありません。既存関数に float32 の変形ヘキサを与え、次の反例を確認しました。

   | 指標 | 値 |
   |---|---:|
   | 体積和 | `0.19999999` |
   | AR | `1.14891` |
   | skew | `0.490353` |
   | 最小頂点 Jacobian | **`−0.175`** |

   加えて、双対体積生成は小四面体の絶対体積を加算するため、`dualVolume>0` だけでは反転を排除できません（[gmshReader.hpp:1913](/home/sano/work/forge/solver_density_cuda/mesh/gmshReader.hpp:1913)）。変換器の閉性ログも、各CVの局所面積ではなく**全域最大面積**で正規化しています（[同:2082](/home/sano/work/forge/solver_density_cuda/mesh/gmshReader.hpp:2082)）。

   **対案:** 本番の物理座標を float32 化した後の局所 Jacobian 検査を独立に実装し、双対閉性は各CVについて `|ΣS|/Σ|S|` で判定してください。品質判定の `SOFT-PASS` は拒否し、AR 5000 の適用は近直交の壁法線構造層に限定すべきです。45°の隅セルまで自動的に緩和対象にはできません。

7. **Major — 実証する3解像度は、宣言した `scale` 細分列ではありません。**

   **根拠:** [plan:106](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:106) は全分割数を `scale` 倍、`h1` をその逆数としています。しかし次行の試験は、第一層が16倍変わる間に層数・接線分割は2倍しか変わりません。壁解像の感度試験にはなりますが、宣言した系統細分機能の検証にはなりません。

   また、[plan:154](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:154) の「短/長」には数値がなく、最も退化しやすい `fillet_taper_len` と `L_cowl-L_sw` の下限近傍も試験対象が定まっていません。

   **対案:** 方式受入には実際の `scale` による3段階を使い、第一層4倍刻みの系列は別の感度試験にしてください。dv の具体値・許容範囲・拒否境界を表にし、短いテーパと短い露出カウル区間を含める必要があります。

8. **Minor — 「有限厚化で密度床を根治」は、このフェーズでは検証できません。**

   **根拠:** [plan:4](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:4) は根治を主張しますが、[plan:29](/home/sano/work/forge/plans/active/tooling-sern-mesh-blocking.md:29) 以降のスコープはメッシュのみです。共有点での境界処理や離散演算が残因である可能性は、メッシュ検査では排除できません。

   **対案:** 目的を「厚さゼロ交線を除去し、床張り付き解消を検証できる形状を作る」に限定してください。根治の判定は B6 の時系列・床ノード検査・正式な収束／準定常 VERDICT に委ねます。

**推奨は、gmsh 方式を維持し、B1 の前に接続設計を完成させることです。** 優先順は、①非退化な CAD と共有面一覧、②全壁の層配置、③実行可能な受入ゲート、④系統細分と dv 試験表です。

x 方向に曲がる単純ブロック自体は、メモリ上の試験で `h1=4e-5H`・40層、全1,280要素がヘキサ、最小 Jacobian `3.04e-8>0` となりました。したがって問題は「gmsh では不可能」ではなく、**今回の接続模型が閉じる設計になっていないこと**です。メッシュのみの方式実証に CFD・周期・軸対称回帰を追加する必要はありませんが、block-DPLUR の安定性や生産条件での受理は B6 まで未検証です。

本レビューは **plan 未反映**です。指示どおりファイルは変更していません。

指摘数: Critical 0 / Major 7 / Minor 1
