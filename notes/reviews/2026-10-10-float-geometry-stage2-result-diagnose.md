# codex 諮問 (diagnose): float-geometry-stage2-result

- **brief**: [`notes/reviews/briefs/2026-10-10-float-geometry-stage2-result.md`](../../notes/reviews/briefs/2026-10-10-float-geometry-stage2-result.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-10
- **commit**: `3de780a4` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 3.3 min, rc=0
- **結論**: **次の計算は、同じFP64のB0保存場から`mesh.hoopAreaFromClosure`だけを0／1に変えた上限2000 stepの短い比較とし、長い到達時間比較はまだ行わない。**
- **extra**: `plans/active/axisymmetric-freestream-hoop-gauge.md`, `plans/active/time_integration-line-implicit-speed.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表です。**段②の無条件合格は保留し、旧基準のFAILを残して改訂基準で再評価する**ことを推奨します。対象runと判定JSONはローカルに存在しないため、以下の実測値は依頼文の報告値として扱い、コードから確認できた事実と区別します。

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | H1「case/48は物差しの問題」：**採用。ただし実装全体の無欠陥までは証明しない** | 成分相対差の分母はほぼゼロの係数も拾う定義です（[fg2_an.py:21](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/fg2_an.py:21)）。報告されたcase/48のベクトル相対差最大8.63e−8は、このFAILが微小成分に由来する説明を支持します。**元の0.0722／FAILを保存し、改訂版の判定を別記**してください。 |
| **Major** | H2「ベクトル相対差に変更」：**採用。ただし追加条件が必要** | ベクトル精度だけでは、伸縮格子で小さい成分が担う勾配方向の精度まで保証できません。改訂基準は非ゼロincidenceの **最大 ‖Δc‖/‖c64‖ ≤ 1e−6**、全係数の有限性、ゼロベクトル・境界incidenceの別検査とします。さらに各節点の線形再現行列 Σ c⊗d を評価し、解像される方向のFP64との差を最大1e−6以下で確認してください。結果を見た後の再評価であることを明記し、以後の段の基準として固定します。 |
| **Major** | §6.3の3・3′・4・5：**限定的に採用** | 3′はfloatの有効なclosure経路の非劣化を補っています。4・5も登録された20 step比較の範囲では報告値が基準内です。ただし **FP64のclosure有効時の新旧比較は未実施**。また、無効時にはclosure配列を計算しないため、ダンプの「float対FP64差0」は精度の証拠になりません（[variables.cpp:665](/home/sano/work/forge-integ-1005/solver_density_cuda/variables.cpp:665)、[geomStage2Dump.cpp:231](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/geomStage2Dump.cpp:231)）。有効経路で実面ベクトル和と配列を照合してください。 |
| **Major** | 判定器のPASSの信頼性：**要再検証** | `neq`は欠損を`None`にし、呼び出し側の`or 0`が不一致0として扱います。本番配列との一致属性も表示だけです（[fg2_an.py:18](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/fg2_an.py:18)、同:45）。3′の非有限件数も合否に入りません（[fg2d_an.py:21](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/fg2d_an.py:21)）。欠損・形状不一致・非有限・本番との不一致を判定不能またはFAILにして、保存済み成果物を再集計するのが先です。 |
| **Major** | H3「float壁際はatomicAddが原因で段③④では減らない」：**要再検証** | 圧力流束の積と原子加算はともに`flow_float`です（[convectiveFlux_slau_d.inc.cuh:622](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:622)、同:668）。しかし積の丸め、面圧力のばらつき、加算の丸めは未分離です。**圧力流束の評価・集計の精度が有力**とは言えますが、atomicAdd単独への帰属はできません。段③④がこの静止場誤差を直接修復するとは期待せず、V3で面流束と集計後残差を分けて測ってください。 |
| **Major** | H4「静止場の偽の力がB0の床を作る」：**要再検証** | 静止場の対象は半径運動量、B0の床の位置は密度残差です。最大位置もx/r_t≈−3.02と−1.79で異なります。さらにclosure切替は **周方向粘性応力に掛ける面積と軸方向補正も変更**します（[axisymmetricSource_d.cu:51](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/axisymmetricSource_d.cu:51)、同:111）。同じB0保存場で変更項の大きさ・符号と、その後の密度残差の応答を測る必要があります。 |

段③へは、上記の係数再評価と判定器の欠落を解消したうえで**条件付きで進めてよい**と判断します。§6.4の`implicitSolvePrecision: 1`による同一性検査は妥当ですが、20 step後の状態比較だけでV2全体を代替しないでください。同一入力での本番面流束・commit前残差・dqの照合を残し、既定のISP 0についても有限性だけでなく係数・ガード・対象面の変更を検査する必要があります。未実行経路の実行検証をV6へ残す扱いは妥当です。

結論: **次の計算は、同じFP64のB0保存場から`mesh.hoopAreaFromClosure`だけを0／1に変えた上限2000 stepの短い比較とし、長い到達時間比較はまだ行わない。**

第1仮説: **FP64でも残る静止場の偽の力は、折れた双対面を集約してから半径を掛ける幾何近似に由来する。** 確度: 中。

根拠: 変換器は各区間の面ベクトルを足し、重心を区間長で平均します（[gmshReader.hpp:1456](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/gmshReader.hpp:1456)）。ソルバはその集約重心の半径を合成面ベクトルに掛けます（[variables.cpp:630](/home/sano/work/forge-integ-1005/solver_density_cuda/variables.cpp:630)）。一般には、

`r̄·Σ Sₖ ≠ Σ rₖSₖ`

です。多角形の厳密な閉性は、各直線区間で半径を掛けて積分した場合に成立します。**「多角形だから現在の集約表現でも厳密」という前提は受け入れません。** 実際、二つの折れた区間を用いた算術例では、doubleでも集約積と区間積分に約0.9%の差が出ました。これはcase/45の実測再現ではありません。

`axisCentroidShift`はCV中心を選ぶ処理で、同じ箇所では面ベクトル・面重心・`dualVolume`を変更していません（[gmshReader.hpp:2228](/home/sano/work/forge-integ-1005/solver_density_cuda/mesh/gmshReader.hpp:2228)）。静止場の圧力不釣合いの直接原因としては優先度を下げます。

反証条件: 実メッシュの双対区間から再構成した`Σ rₖSₖ`と、現行の`r̄ΣSₖ`の差が、j=2〜4で観測された閉性欠損の符号・大きさを説明しないこと。比較は同じCVについて両正規化、`A_planar`と`Σ|rS|`、で行います。説明できなければ、この仮説を下げて面積・向き・集約接続を点検します。

第2仮説: case/48の係数FAILは、ほぼゼロの成分で割る指標の問題である。確度: 高。ただし線形再現性は未確認。

第3仮説: float壁際の残差は、大きな圧力流束の評価・集計の丸めが支配する。確度: 中。積・面圧力・加算の内訳は未確認。

判別A/B:

- **出発点**：`case/45.isobutane_m6_d155/run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2/res_60000.h5`。同一メッシュの`restart_field.py`で二つの新規runへ複製し、同じ固定済み段②FP64バイナリを使用します。A=`hoopAreaFromClosure: 0`、B=`1`。`pRef: 0`を含む他設定は固定します。
- **長さ**：最初のcommit前残差を保存し、その後2000 stepまで。場は1・100・500・1000・1500・2000 stepで保存。これは**床付近の局所応答の診断**であり、元の初期場から水準に届くstep数の比較ではありません。
- **同一状態での照合**：非軸CVの差が、`ΔR_y = (P + pk − pRef − τθθ)(A_closure_y − A_planar)`、`ΔR_x = (P + pk − pRef)A_closure_x`で説明できることを確認します。圧力寄与と粘性寄与は別記します。同じ状態の密度残差は、このソース切替だけでは直接変わりません。
- **事前判定**：領域`x/r_t∈[−5,0)`を固定し、j=2〜4とj=20〜59を別集計。末尾500 stepの全`rms_*`、領域内Σ|res_ro|、最大位置・値、θ_r(40/70/94)・符号付きQ_wの変化を比較します。領域内Σ|res_ro|が対照の半分以下になり、再実行差を超える持続的な低下なら、closure補正が床へ強く影響する仮説を支持します。比が0.9以上なら「短時間で床を大幅に下げる」という期待を棄却し、延長しません。中間は判別不能です。
- θ_r差0.05%・Q_w差0.1%超は**解への影響のフラグ**とし、改善・悪化を断定しません。各腕の`check_convergence`と、当該区間の`check_quasisteady`のVERDICTを記録します。短い比較から定常値や総到達時間は主張しません。

やらない方がよいこと: case/48の旧FAILをPASSへ上書きすること、`hoopAreaFromClosure`と`pRef`を同時に変更すること、静止場の改善だけでB0の床の原因を確定すること、直ちに12万〜20万stepを回すこと、`axisCentroidShift`を先に切り替えること。

呼び出し側の前提への異議: closureの無効配列の差0は精度の証拠ではありません。また、静止場の`res_roUy`が消えても、移流・粘性・乱流があるB0の`res_ro`の床が消えるとは限りません。段③④による幾何差分の改善と、圧力流束の丸め・状態commitの丸めは別々に評価すべきです。

不足情報: AWS上の実効設定・判定JSON・有限性検査結果、微小係数の再現可能な内訳集計、case/45の区間別幾何モーメントが不足しています。現在の作業ツリーには段③の変更もあるため、測定済み段②バイナリに対応するソースの固定も必要です。ファイル変更・forge起動は行っていません。**本レビューの提案はplan未反映**であり、呼び出し側でfloat-geometry plan §6.3〜§6.5と、hoop-gauge planへ採否を記録してください。
