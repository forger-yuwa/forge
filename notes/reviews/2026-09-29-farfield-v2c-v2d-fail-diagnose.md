# codex 諮問 (diagnose): farfield-v2c-v2d-fail

- **brief**: [`notes/reviews/briefs/2026-09-29-farfield-v2c-v2d-fail.md`](../../notes/reviews/briefs/2026-09-29-farfield-v2c-v2d-fail.md)
- **plan**: [`plans/active/boundary-node-farfield-characteristic.md`](../../plans/active/boundary-node-farfield-characteristic.md)
- **date**: 2026-09-29
- **commit**: `71442b34` (feature/sern-design)
- **codex**: effort `high`, 4.4 min, rc=0
- **結論**: **V2c・V2d-2 の旧判定は FAIL のまま保持し、次は V2c の低・高領域へ同一の局所状態を与え、角近傍の初回残差だけを比較する A/B を行う。**
- **extra**: `case/58.farfield_verification/README.md`, `case/58.farfield_verification/setup_v2c.py`, `case/58.farfield_verification/setup_v2d.py`, `solver_density_cuda/cuda_forge/convection/farfieldFlux_d.inc.cuh`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **V2c・V2d-2 の旧判定は FAIL のまま保持し、次は V2c の低・高領域へ同一の局所状態を与え、角近傍の初回残差だけを比較する A/B を行う。**

採否表:

| 重大度 | 諮問事項 | 採否・根拠・対案 |
|---|---|---|
| **Major** | V2c を「境界の実力 ≈6 %反射」と確定し V3 へ進む | **却下**。評価器が測るのは壁圧の領域間差であり、入射・反射波を分離した反射係数ではない（[eval_v2c.py:43](/home/sano/work/forge-sern-design/case/58.farfield_verification/eval_v2c.py:43)）。また V3 は V1–V2 合格後と明記される（[plan:185](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:185)）。対案は「基準超過、原因未分離」と記録して角の交絡を先に調べること。 |
| **Major** | 角を評価から除外する | **却下**。事前条件は「壁面、全線」（[plan:179](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:179)）。A≈B は上面 BC への感度が小さい証拠にはなるが、C との比較が正しい証拠にはならない。原因を同定・修正して全線判定をやり直す。除外版は補助診断として併記できるが、旧 FAIL を置き換えない。 |
| **Major** | `Z=ρU/√(M²−1)` へ変更する | **却下（現段階）**。既存評価器の衝撃波関係式から再計算すると、`Δp/P₁=0.864`、衝撃後 `M₂=2.086`、上面に対する `Mₙ≈0.362`。小振幅・接線平均流の仮定から外れる。6.52 %という見積もりの算術は再現するが、本試験への適用は未検証。現行 TRRS は維持し、誤差の所在を先に確かめる。多次元の波の分離が必要という理論的背景は [Giles の線形解析](https://people.maths.ox.ac.uk/~gilesm/files/aiaa90.pdf)とも整合する。 |
| **Major** | V2d-2 を「右端擾乱ゼロ」として合格に変更する | **要再検証**。`--left-hot` は左端の接触面自体を消す操作（[setup_v2d.py:105](/home/sano/work/forge-sern-design/case/58.farfield_verification/setup_v2d.py:105)）。右端流出を隔離する対照としては採用できるが、元の左右冷気条件を合格にはできない。旧 FAIL と「左端高温対照では表示精度内で擾乱なし」を別記する。左端流入の原因帰属も未確定とする。 |

第 1 仮説: **V2c の角の差には、領域拡張に伴う番号付け・変換後の局所離散作用素の差が混入している。** 確度: **中**
  
  根拠: 座標生成式は共通（[setup_v2c.py:38](/home/sano/work/forge-sern-design/case/58.farfield_verification/setup_v2c.py:38)）だが、節点番号は `ny` に依存する（[make_box_msh.py:13](/home/sano/work/forge-sern-design/case/58.farfield_verification/make_box_msh.py:13)）。双対面は節点番号をキーとして構築・整列される（[gmshReader.hpp:1806](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/gmshReader.hpp:1806)）。これは不具合の証明ではないが、「同一座標だから同一作用素」とは確認できていない。ブリーフの A≈B、C のみ角で約 7618 Pa 違うという観測とも矛盾しない。

  反証条件: 同じ保存量を座標対応で与えた初回評価において、角の局所幾何・再構成状態・面流束・残差が下記許容内なら、**その状態で局所空間作用素が異なるという仮説**を棄却する。以後は反復過程・解の履歴を調べる。

第 2 仮説: **V2d-2 は左端で生成される接触面の TP 保存形混合が音波源になっている。** 確度: **中**。境界を含まない既存ホスト関数で `0.9×高温保存量＋0.1×外気保存量` を再計算し、単成分 AIR で **+14.634 Pa（+0.513 %）**、多成分で **+18.641 Pa（+0.654 %）**を確認した（[farfield_proto1d_tp.py:230](/home/sano/work/forge-sern-design/solver_density_cuda/tools/farfield_proto1d_tp.py:230)）。ただし、これは run の **2.49 Pa**を再現した結果ではなく、境界実装固有の寄与は未除外。

第 3 仮説: **V2c 下流の差には現行境界の斜入射に対する反射が含まれる。** 確度: **中**。現仕様もその限界を明記する（[methods/boundary.md:150](/home/sano/work/forge-sern-design/methods/boundary.md:150)）。ただし、壁圧差から反射率 **6 %**と定量化する部分は未確認。

判別 A/B: **変更するのは領域高さだけ**とする。A＝低領域・上面 slip、B＝高領域・上面 slip。既存 C の有限なスナップショットを共通の入力場とし、低領域には座標が一致する節点の保存量を正確に移す。収束解とは仮定せず、作用素検査の入力として使う。

- 各新規 run は **1 step、判定は更新前の初回残差評価だけ**。同一バイナリ・数値設定・CUDA block sizeを用いる。
- `x=0.195–0.230 m` の下面と再構成に必要な近傍を対象に、体積・接続・面積ベクトルを照合し、帳簿の `after_eos_bc`、面ごとの再構成状態・流束、`res_after_conv` を比較する。既存診断で取得できる（[main.cpp:1469](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1469)、[convectiveFlux_d.cu:615](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:615)）。
- 判別閾値は事前に、各保存量の残差差を接続面の絶対流束和で規格化して **`1e-5`**。ゼロ規模の成分には非零の自由流基準を使う。

→ **差が閾値を超えるなら**局所比較の交絡を採用し、最初に異なる幾何・再構成・面流束へ絞る。**全項目が許容内なら**第 1 仮説の直接的な作用素差を棄却し、反復過程へ切り分けを移す。どちらでも、この 1 評価から収束や境界性能は判定しない。

やらない方がよいこと: **角をマスクして合格化する、2 %を観測後に緩める、6 %に合うよう境界式を調整する、V3 で V2 の未受入れを代替すること。** また、亜音速の既存入口へ交換する試験は閉包・ピン処理まで変わるため、TP 混合だけの判別にはならない。

呼び出し側の前提への異議:

- 同じ `.msh` の再変換は**再現性**の確認であり、決定的な変換誤差・番号依存性の除外ではない。CFL 半減も全ての経路依存性を除外しない。
- 角の圧力が一定でも、全場の収束・壁圧差の準定常性は未確認。`eval_v2c.py` は最終スナップショットだけを使う（同ファイル:24）。
- V2d の「反射 PASS」も再評価が必要。評価器は全時間窓で `max|P_long−P∞|` を入射振幅とする（[eval_v2d.py:65](/home/sano/work/forge-sern-design/case/58.farfield_verification/eval_v2d.py:65)）。設定から、左端発の音波到達は `0.6154/cᵢ`、右端からの反射到達は `0.7473/cᵢ`。**反射評価前に左端擾乱が入る**ため、短長差の相殺だけでは隔離試験にならない。
- plan §5.1 #3 には旧 FAIL、角の原因未分離、TP 混合仮説、VERDICT 原本未確認を残すべき。§6 の旧基準は保持する。**plan 未反映（依頼どおりファイル変更なし）。**

不足情報: AWS の対象 run 本体、実効設定・バイナリ識別情報、メッシュ品質判定、`CONVERGENCE_VERDICT.txt`、判定区間、壁圧時系列に対する `check_quasisteady.py` の VERDICT。手元は入力スクリプトのみ（[case README:5](/home/sano/work/forge-sern-design/case/58.farfield_verification/README.md:5)）であり、ブリーフの run 数値は独立検証できていない。
