# codex レビュー: boundary-cht-axisymmetric-fem2d (plan)

- **plan**: [`plans/active/boundary-cht-axisymmetric-fem2d.md`](../../plans/active/boundary-cht-axisymmetric-fem2d.md)
- **stage**: `plan`
- **date**: 2026-09-27
- **commit**: `02b77df4` (feature/cht-axisym-fem2d)
- **codex**: effort `high`, 5.1 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m1

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**GO-with-changes**

`r` 重みの弱形式・Robin 行列・積分済み荷重を直接渡す方針は妥当です。親計画からの切り出しにも重複はありません。  
ただし、診断修正の検出力、丸め床判定、再開・後処理への軸対称情報の伝達が不足しており、以下を実装前に計画へ反映すべきです。

1. **Major — FP64・軸対称の残差を、単精度・平面用の丸め床で判定しようとしている。**

   根拠: [plan:165](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:165) が参照する `solver_density_cuda/tools/rounding_floor.py` は存在しません。見つかった [case/57 の実装:19](/home/sano/work/forge-cht/case/57.transition_flat_plate/tools/rounding_floor.py:19) は未加重の面ベクトルを読み、[同:31](/home/sano/work/forge-cht/case/57.transition_flat_plate/tools/rounding_floor.py:31) で `eps=1.1920929e-7` を固定しています。対象も質量・エネルギーだけで、運動量残差を評価しません。

   同じ流束尺度でも、この `eps` は FP64 の **536,870,912 倍**です。「全残差が丸め床の10倍以内」という合格条件を、このツールでは検証できません。

   **対案:** 現段階では丸め床による合格例外を外し、`NOT CONVERGED` は不合格にする。例外を必要とするなら、実際の精度・`r` 重み・各保存量の残差尺度に対応した評価器の実装と検証を、先行する必須作業として登録してください。

2. **Major — `iface_q_eff` の修正漏れを、現在の試験では検出できない。**

   根拠: [conjugateWall.cpp:715](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:715) の `fem2d` 経路は `iface_Qf_eff` を直接使い、`iface_q_eff` を経由しません。一方、[V-ax2:158](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:158) は壁温・収支・界面残差を評価し、`q` には準定常性しか要求していません。

   登録した円環の解析値は **91.85535 W/m²**です。既存の平面面積除算を残すと、半径 `0.01 m` の壁では数値が **0.9185535**になります。しかし積分済み荷重と壁温は変わらず、誤った `q` も定常なら準定常試験を通ります。

   **対案:** V-ax2 に全節点の `iface_q_eff` 対解析値の誤差判定を追加する。別途、半径可変面・`axisRFloor` 適用面で `q_eff A_fluid^r = Qf_eff` を直接照合し、`q_eff_raw` も検査する。**分母修正を意図的に外した実装が FAIL すること**まで要求してください。

3. **Major — 再開ファイルが平面と軸対称を識別できず、単位の違う荷重履歴を復元する。**

   根拠: [solid_mesh_to_h5.py:109](/home/sano/work/forge-cht/solver_density_cuda/tools/solid_mesh_to_h5.py:109) の `content_sha1` は座標・接続・物性・Robin 条件から作られ、軸対称設定を含みません。[conjugateWall.cpp:493](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:493) はそのハッシュを照合した後、[同:512](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:512) で `SOLID/QBUF` を復元します。

   実際に `case/58.conjugate_slot/run_0012_v6p_df5_i50_nlog/` の固体・再開ファイルの属性を確認しましたが、幾何方式・荷重単位の識別子はありません。同一固体メッシュで平面→軸対称を切り替えると、現契約では **[W/m] の履歴を [W/rad] として受理**します。

   **対案:** 再開状態に幾何方式・荷重単位・契約バージョンを保存し、不一致は拒否する。軸対称では識別属性のない旧状態も拒否する。同一方式での連続実行対再開、方式切替の拒否を検証へ追加してください。

4. **Major — 既存の収支評価器と外部連成経路に、平面の計算が残る。**

   根拠: [check_cht_balance.py:163](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:163) は `Fem2DOperator` を幾何方式なしで生成し、さらに [同:272](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:272) で Robin 放熱を平面式から直接計算します。オラクルへ `axisym=False` を追加するだけでは直りません。

   V-ax2 の解析温度をこの式へ代入すると、正しい **0.001837107 W/rad** に対し、旧式は数値として **0.09185535**、つまり50倍を返します。また、[cht_loop.py:281](/home/sano/work/forge-cht/solver_density_cuda/tools/cht_loop.py:281) の `--axisym` は `shell2d` にしか渡されず、`fem2d` は平面のままです。

   **対案:** `check_cht_balance.py` を§5・§7の必須変更対象に追加し、保存状態と幾何方式を照合して重み付き積分を行う。外部 `cht_loop.py` の軸対称 `fem2d` は、本計画で対応しないなら明示的に拒否してください。C++ 出力だけの修正では不十分です。

5. **Major — 最重要の局所整合性試験 V-ax2b に、格子条件と細分化判定がない。**

   根拠: [V-ax2b:172](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:172) は寸法を指定していますが、流体・固体の分割数、界面節点配置、三角形の分割方向を指定していません。今回の不整合は、まさにその離散幾何に依存します。

   読み取り専用の独立計算で、同節の固体円板を厚さ方向16分割とし、一様熱流束を流体半辺積分の荷重として与えました。界面温度の最大誤差は、半径方向 **8／16／32分割で 0.9421／0.2983／0.08796 %**でした。consistent FEM 荷重なら相対誤差は約 `1e-14` です。これは **forge の連成 run ではなく、荷重移送だけを分離した試験**ですが、同じ設計が格子選択だけで0.5%判定を跨ぐことを示します。

   **対案:** 分割数・節点配置・三角形方向を事前固定し、半径方向の細分化系列と非一様格子を追加する。端点を含む最大誤差とその減少を判定してください。既存方針どおり、結果を見て荷重を面積比で掛け直す対応は採用しないこと。

6. **Minor — FP64で検証する範囲と、機能を提供する範囲を明記すべき。**

   根拠: [plan:127](/home/sano/work/forge-cht/plans/active/boundary-cht-axisymmetric-fem2d.md:127) は全判定を FP64 に限定しますが、受理条件には精度の限定がありません。現行の [flowFormat.hpp:6](/home/sano/work/forge-cht/solver_density_cuda/flowFormat.hpp:6) は流体・幾何とも `float` です。

   **対案:** 今回の保証範囲を「定常・node・`axisymMethod: 0`・FP64で検証済み」と明記し、FP32での精度保証は未検証とする。cell・dual-time 等の既存拒否も維持すると§4.4に明記してください。node のみを検証する方針自体は、現行の検証手順と整合します。

**推奨は、現行の重み付き FEM と積分済み荷重の直接受け渡しを維持し、上記1→2→3→4→5の順で検証・データ契約を補ってから実装へ進むことです。** 数学的な中核を取り替える必要はありません。Robin 行列は独立した数値積分でも相対差 `5.4e-16` 以下、円環の登録解析値も再計算で確認しました。

ファイルは変更していません。提案は **plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
