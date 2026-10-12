# codex レビュー: tooling-nozzle-wt-deliverables (plan)

- **plan**: [`plans/active/tooling-nozzle-wt-deliverables.md`](../../plans/active/tooling-nozzle-wt-deliverables.md)
- **stage**: `plan`
- **date**: 2026-10-13
- **commit**: `90a09685` (feature/nozzle-deliverables)
- **codex**: effort `high`, 5.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m1
- **focus**: §4.2a 寸法ゲート・§4.3 構造ファイル・§4.4 Fluent 書き出し (2D の向きは実物で測定済み)・§4.5 熱流束・§6 の合格条件が自己一致になっていないか。methods/design/overview.md は「実寸の要求カードと成果物一式」節だけ読む (grep -n で位置を探す)
- **extra**: `methods/design/overview.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
目的と既存機能の再利用方針は妥当です。ただし、軸対称熱流束の単位、寸法ゲート、CFD データの採用条件に実装前の修正が必要です。  
特に V5 は、熱量が誤っていても面積と列名の検査だけで合格します。

ファイルは変更していません。指定された既存 run の出力・CAD 検証環境はこのワークスペースに無いため、コード照合と書き込みを伴わない数値反例で確認しました。CFD の収束・準定常性や Fluent での読み込みを確認したという主張はしていません。

1. **Major — §4.5: 軸対称の `iface_q_eff` をそのまま W/m² として使えない**

   根拠: [`variables.cpp:589`](/home/sano/work/forge-deliv/solver_density_cuda/variables.cpp:589) は `axisymMethod: 0` の流束面積に半径を掛けます。一方、[`conjugateWall.cpp:229`](/home/sano/work/forge-deliv/solver_density_cuda/conjugateWall.cpp:229) は、その流束・残差から作った熱量を、半径を掛けていない `msh.planes[ip].surfArea` で割ります。したがって、この経路の `iface_q_eff` は物理的な熱流束への換算が済んでいません。

   コードの式を円筒壁で評価すると、半径 0.025 m、区間長 0.001 m、熱流束 1000 W/m² に対してダンプ相当値は 25。計画どおりリング面積を掛けると 0.003927 W となり、正しい 0.15708 W の **1/40** です。これは実 run の測定値ではなく、実装式の数値反例です。

   **対案:** [`conjugateWall.cpp:292`](/home/sano/work/forge-deliv/solver_density_cuda/conjugateWall.cpp:292) の積分済み `iface_Qf_eff` を起点に、離散化・軸対称方式ごとの単位換算を明文化してください。半径重み方式では全周の節点熱量を `2π·iface_Qf_eff` とし、CSV 用の熱流束はこれを明示した節点持ち分面積で割ります。V5 に、異なる半径の円筒と傾斜壁での**局所熱流束・積分熱量・符号**の検証を必須追加すべきです。平面の case/48 だけでは検出できません。

2. **Major — §4.2a: 直管長の測定方法が正しい形状を不合格にする**

   根拠: [`plan:144`](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wt-deliverables.md:144) は、入口半径との差が 10 µm を超える位置を縮流開始点としています。しかし既存の縮流曲線は接合点で勾配・二階微分がともにゼロです（[`wall_walldriven.py:79`](/home/sano/work/forge-deliv/design/forge_design/geometry/wall_walldriven.py:79)）。半径差が閾値に達する位置は、接合点より下流になります。

   このクラスで計算した反例は次のとおりです。

   | 条件 | 真の直管長 | 計画の測定値 | 誤差 |
   |---|---:|---:|---:|
   | M5 カードの寸法、設計半径 25 mm、`R=3` | 200 mm | 214.924 mm | +14.924 mm |
   | case/45 の YAML にある設計寸法 | 38.404 mm | 56.018 mm | +17.615 mm |

   縮流部長も同じ量だけ短く判定され、長さ公差 1 mm を超えます。最終物理壁の上流多項式も同じ平坦な接合条件を持ちます（[`wall_axismach.py:451`](/home/sano/work/forge-deliv/design/forge_design/geometry/wall_axismach.py:451)）。

   **対案:** STEP の区間係数から定半径区間の終端を同定し、形状表現の数値誤差と寸法公差を分離してください。メッシュ側は接合点を挟む区間と測定不確かさを返す方式にします。また、`diameter_design_m` 指定時は物理最小径との比較を必須にせず、設計径の検証と物理径の参考測定を分ける必要があります。

3. **Major — §4.2・§4.5・V5/V6: データセットの存在と、熱設計に使える結果であることを区別していない**

   根拠: [`plan:191`](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wt-deliverables.md:191) の分岐は `iface_q_eff` の有無だけです。しかし [`conjugateWall.cpp:192`](/home/sano/work/forge-deliv/solver_density_cuda/conjugateWall.cpp:192) は配列を NaN で初期化し、適用条件を満たす場合だけ値を設定して、適用外でもデータセットを出します。例えば断熱壁では「存在するが無効」が起きます。

   また、計画には壁解像ゲートはありますが、選択した `--res` と結び付いた収束・熱流束の準定常判定がありません。[`plan:249`](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wt-deliverables.md:249) の V5 は、面積と列名が正しければ通ります。

   **対案:** 欠損・非有限・適用外・診断値・検証済みを区別し、無効な正本から診断値への黙った置換を禁止してください。熱設計用の合格には、対象 snapshot と同一区間の `check_convergence.py`、対象熱量時系列の `check_quasisteady.py`、壁解像の判定を要求します。CFD 前や未収束の書き出しは許しても、MANIFEST は「未評価／未合格」とし、熱伝導シートにも状態を伝えるべきです。NaN・未収束・異なる step の混入を負例に追加してください。

4. **Major — §4.3: CAD 編集後も `u=x` が成立するという前提がない**

   根拠: [`plan:161`](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wt-deliverables.md:161) は助変数範囲だけを調べ、`u=x_i` で評価します。`x(u)=u` は既存 exporter がグレビル点で作る曲線の性質です（[`wall_step.py:46`](/home/sano/work/forge-deliv/design/forge_design/export/wall_step.py:46)）。同節が約束する CAD 編集後の曲線には保証されません。

   例えば同じ端点・助変数範囲 `[0,100]` を持つ単調曲線でも、`x(u)=u+0.1u(100−u)/100`、`r(u)=10+0.1u` なら、station 50 mm に `u=50` を使うと半径が **0.24938 mm** ずれます。範囲検査は通ります。

   **対案:** STEP の実座標から `x(u)=x_i` を解き、得た助変数を壁節点にも設定してください。平面性・正半径・x の単調性・全 station の包含を検査し、多価曲線を拒否します。再パラメータ化した同一形状を正例、折り返す曲線を負例に追加すべきです。

5. **Major — §4.4: Fluent の圧力入力に必要な基準圧が抜けている**

   根拠: [`plan:179`](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wt-deliverables.md:179) は全圧・背圧・初期ゲージ静圧を列挙しますが、`Operating Pressure` と絶対圧からの換算を定義していません。Fluent の入口全圧は operating pressure に対するゲージ値です。[Ansys 公式境界条件リファレンス](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/flu_ug/flu_ug_bcs_sec_bound_cond.html)

   **対案:** 配布する設定表では `Operating Pressure = 0 Pa` を明示し、絶対圧と Fluent 入力欄の値を併記してください。V4 に圧力換算の検査を加えます。面の向きについては、実物の 2D メッシュを基準にした改訂方針を維持してよく、今回の指摘はその規則の否定ではありません。

6. **Minor — §4.2a と V3 で「弦誤差」と「再構築誤差」が混同されている**

   根拠: [`plan:150`](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wt-deliverables.md:150) はメッシュと STEP の寸法差を「§6 V3 の弦誤差」で判定するとしています。しかし [`plan:240`](/home/sano/work/forge-deliv/plans/active/tooling-nozzle-wt-deliverables.md:240) の V3 が測るのは、元メッシュと再構築メッシュの**節点位置差**です。節点間の曲線と弦の差や、最小半径の取り逃しは測れません。また「品質 VERDICT が元と同じ」では両方 FAIL でも条件を満たします。

   **対案:** CAD 転送誤差、再構築節点誤差、曲線と弦の形状誤差を別の検査量にしてください。品質判定は、元と再構築後の双方が明示した閾値で PASS することを要求します。

推奨は、**既存 exporter を再利用する構成を維持し、上記を計画に反映してから実装すること**です。実装前の優先順は、①熱量の単位・保存性、②寸法測定、③結果の採用条件、④STEP の逆評価、⑤Fluent の圧力基準、⑥検証量の整理です。`plans/README.md` と関連 accepted plan を確認した範囲では、成果物統合は未解決であり、既存 STEP 機能の重複実装を避ける現在の分担は妥当です。

レビュー専用の指示に従い、これらの提案は **plan 未反映** です。

指摘数: Critical 0 / Major 5 / Minor 1
