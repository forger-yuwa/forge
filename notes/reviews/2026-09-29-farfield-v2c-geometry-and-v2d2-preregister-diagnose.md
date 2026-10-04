# codex 諮問 (diagnose): farfield-v2c-geometry-and-v2d2-preregister

- **brief**: [`notes/reviews/briefs/2026-09-29-farfield-v2c-geometry-and-v2d2-preregister.md`](../../notes/reviews/briefs/2026-09-29-farfield-v2c-geometry-and-v2d2-preregister.md)
- **plan**: [`plans/active/boundary-node-farfield-characteristic.md`](../../plans/active/boundary-node-farfield-characteristic.md)
- **date**: 2026-09-29
- **commit**: `55de2fb9` (feature/sern-design)
- **codex**: effort `high`, 4.8 min, rc=0
- **結論**: **新形状の投入より先に、既存の `run_0097`／`run_0098` の帳簿を「slip 法線残差の除去なし」で再集計する。現在の残差評価器は、ソルバが更新する運動量成分を捨てている。**
- **extra**: `case/58.farfield_verification/setup_v2c.py`, `case/58.farfield_verification/v2c_residual_ab.py`, `case/58.farfield_verification/setup_v2d.py`, `case/58.farfield_verification/eval_v2d.py`, `notes/reviews/2026-09-29-farfield-v2c-corner-and-v2d2-redesign-diagnose.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **新形状の投入より先に、既存の `run_0097`／`run_0098` の帳簿を「slip 法線残差の除去なし」で再集計する。現在の残差評価器は、ソルバが更新する運動量成分を捨てている。**

採否表（Critical 0、Major 4、Minor 1）。run の実測値は原本未確認であり、以下ではコードで確認した事実と区別する。

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | 絶対残差の接線射影：**却下** | [v2c_residual_ab.py:93](/home/sano/work/forge-sern-design/case/58.farfield_verification/v2c_residual_ab.py:93) は「更新時に捨てられる」として `slip` の法線運動量残差を除く。しかし、[mesh.cpp:979](/home/sano/work/forge-sern-design/solver_density_cuda/mesh/mesh.cpp:979) が強拘束対象にするのは `wall`／`wall_isothermal` だけで、[残差射影:56](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:56) と [陰解法の行拘束:1033](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1033) もそのフラグに限定される。`slip_d` は外側・面状態の速度を作る処理である（[boundaryCond_d.cu:68](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/boundaryCond_d.cu:68)）。**帳簿の `res_final` を追加射影せず全5成分で評価する。** |
| **Major** | 「絶対残差 ≤ `1e-5` **または** 収束 PASS」：**却下** | [ブリーフ:26](/home/sano/work/forge-sern-design/notes/reviews/briefs/2026-09-29-farfield-v2c-geometry-and-v2d2-preregister.md:26) は、前回の診断閾値を収束認定の代替にしている。評価器自身も「収束基準の代用ではない」と明記する（[同スクリプト:7](/home/sano/work/forge-sern-design/case/58.farfield_verification/v2c_residual_ab.py:7)）。**修正した局所残差検査、`check_convergence.py`、壁圧の準定常判定を別々の必要条件として残す。** 継続計算の床判定も、参照自身の通常 PASS が必要（[check_convergence.py:359](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_convergence.py:359)）。 |
| **Major** | 新形状・2本立て：**採用**。凸角を「試験形状側の問題」と確定：**却下** | ランプ延長と共通領域の座標式は提案に対応している（[setup_v2c.py:21](/home/sano/work/forge-sern-design/case/58.farfield_verification/setup_v2c.py:21)、同:52）。ただし報告された最大位置は折れ点 `x=0.700` そのものではなく `x=0.715`。これは**凸角近傍に残差が残った観測**であり、幾何・再構成・反復・境界のどれが原因かは未分離。新形状は追加の受入れ配置として使い、旧 FAIL と「凸角近傍の未解決問題」を §5.1 #3c に保持する。node slip の確定欠陥ともまだ呼ばない。 |
| **Major** | V2d-2 の窓付き定義：**採用**。現状の実行長：**要修正** | 正しい中心時刻は `t_inc=0.4/(c+u)`、`t_ref=0.6/(c+u)+0.2/(c−u)`。しかし [setup_v2d.py:75](/home/sano/work/forge-sern-design/case/58.farfield_verification/setup_v2d.py:75) の終了時刻は提案窓を覆わない。コードの定数から再計算すると CPG は **終了 1.64552 ms、窓末尾 1.67017 ms、28 step 不足**。現行の [acoustic_win:83](/home/sano/work/forge-sern-design/case/58.farfield_verification/eval_v2d.py:83) もこれを拒否する。短・長・パルスなし対照すべてを窓末尾まで生成する。 |
| **Minor** | 出口高さによる「閉塞回避の保証」：**要再検証** | 衝撃波で全圧が低下するため、入口の等エントロピー臨界高さ `0.1896 m` だけでは保証にならない。[NASA の斜め衝撃波関係式](https://www.grc.nasa.gov/WWW/BGH/oblique.html) による幾何見積りでは、上面到達 `x=1.00483`、反射衝撃の下面到達 `x=1.51569`、C の上面到達 `x=1.97063` で、配置意図は妥当。**幾何案は採用し、実際の出口で `outflow` に必要な法線超音速流出を確認する**（[methods/boundary.md:161](/home/sano/work/forge-sern-design/methods/boundary.md:161)）。 |

事前登録には、次も具体化する必要がある。

- **V2c**：壁圧の「変動幅」は、固定した末尾区間における  
  `max_x(max_t p(x,t) − min_t p(x,t)) < 0.001 Δp`  
  と定義する。各評価点の系列を判定し、空間最大圧力だけの系列で代用しない。`check_quasisteady.py` は平均値で規格化するため、`STEADY` 単独ではこの絶対許容を保証しない（[同ツール:280](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_quasisteady.py:280)）。500 step 間隔に加えて、末尾窓・必要標本数・延長刻み・打切り条件を数値で固定する。
- **V2c の2系列**：一様 IC 系列と C 共通初期場系列は、それぞれ全線 `0.02 Δp` と対照 A の `0.1 Δp` を判定する。片方の成功で他方を置換しない。C の定常性確認前に「収束済み共通初期場」と扱わない。
- **V2d-2**：短長で実際のプローブ節点・共通領域の離散化・物性・時間刻みを一致させ、両窓の完全なデータ被覆、有限性、非零の入射振幅を検査する。パルスなし対照の許容は **0.02851 Pa**。V2d 条件での時間精度確認も必要で、V2a の CPG 結果だけでは TP を保証しない。
- **共通**：メッシュ品質、NaN/Inf、真空置換・HLL 退避・非有限流束置換ゼロの既存条件を維持する。独立参照比較を境界の合否に入れない整理は採用するが、§6 で要求された独立参照の精度確認と比較記録は残作業として残す。

第 1 仮説: **圧縮角を残差上の問題から除外した判断には、後処理による法線運動量残差の隠蔽が含まれている可能性がある。** 確度: **中**。  
根拠: 上表のソルバと評価器の不一致。射影が不適切なことはコードで確認できたが、隠された残差の大きさは未確認。  
反証条件: 射影前の圧縮角 `x=0.195–0.230` でも、両状態の全5成分が同じ規格化で `1e-5` 以下なら、この局所診断への影響仮説を棄却する。

第 2 仮説: **凸角を除去しても、上面 `farfield` の斜入射への応答による壁圧差は残り得る。** 確度: **中**。仕様も斜入射・強い衝撃波の限界を認めている（[methods/boundary.md:150](/home/sano/work/forge-sern-design/methods/boundary.md:150)）。今回の誤差への寄与率は未確認。

判別 A/B: **ソルバ再実行なし。既存帳簿の追加射影だけを ON／OFF にする。** 対象は `case/58.farfield_verification/run_0097_v2c_resab_fromU/` と `run_0098_v2c_resab_fromC/` の初回評価。分母・対象節点・閾値を固定し、全域と圧縮角・凸角近傍・出口交線を別集計する。  
→ **射影前だけ圧縮角で閾値超過なら**、「圧縮角は残差上問題なし」を撤回する。  
→ **射影前も圧縮角で全成分が許容内なら**、第1仮説をその精度で棄却する。どちらでも、密度残差による旧系列の FAIL は維持する。

やらない方がよいこと: **診断閾値を収束認定へ昇格する、凸角除去を原因同定と呼ぶ、評価線や時間窓を結果に合わせて切る、旧 FAIL を新試験で上書きすること。**

呼び出し側の前提への異議: 「角の全5成分が小さい」は現状では**追加射影後の値**についての主張である。「最大残差が凸角近傍にある」から「原因は凸角だけ」への推論も受け入れない。上記修正後に §6 へ事前登録して進める方針は妥当。**plan 未反映（依頼どおりファイル変更なし）。**

不足情報: 手元の対象 `run_*` は **0件**で、実体は AWS にある（[case README:5](/home/sano/work/forge-sern-design/case/58.farfield_verification/README.md:5)）。帳簿原本、実効設定とバイナリ識別、全残差履歴、場・壁圧時系列、各 VERDICT と判定区間が不足している。**今回、run の数値再現・収束認定は行っていない。**
