# codex 諮問 (diagnose): farfield-v2c-corner-and-v2d2-redesign

- **brief**: [`notes/reviews/briefs/2026-09-29-farfield-v2c-corner-and-v2d2-redesign.md`](../../notes/reviews/briefs/2026-09-29-farfield-v2c-corner-and-v2d2-redesign.md)
- **plan**: [`plans/active/boundary-node-farfield-characteristic.md`](../../plans/active/boundary-node-farfield-characteristic.md)
- **date**: 2026-09-29
- **commit**: `71442b34` (feature/sern-design)
- **codex**: effort `high`, 3.1 min, rc=0
- **結論**: **次は同じ低領域 B の格子・BC に、`run_0064_v2c_B` と `run_0094_v2c_B_fromC` の最終保存量をそれぞれ与え、初回の絶対残差を比較する A/B で、「残差を残した停止」と「異なる離散定常解の候補」を分ける。**
- **extra**: `case/58.farfield_verification/v2c_operator_ab.py`, `case/58.farfield_verification/eval_v2c.py`, `case/58.farfield_verification/setup_v2d.py`, `notes/reviews/2026-09-29-farfield-v2c-v2d-fail-diagnose.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（**Major 4 件、Minor 1 件**。run 数値はブリーフ記載値であり、原本による独立検証はできていない）。

| 重大度 | 諮問事項 | 採否・根拠・対案 |
|---|---|---|
| **Major** | V2c を C 共通初期場で判定し直す | **追加の隔離試験として採用。旧試験の訂正・置換は要再検証**。`fromC()` は保存量を座標対応で移しており、初期場を揃える操作として妥当（[v2c_operator_ab.py:145](/home/sano/work/forge-sern-design/case/58.farfield_verification/v2c_operator_ab.py:145)）。ただし「異なる離散定常解がある」は未証明。作用素 A/B が測ったのは**残差の差**であり、各状態の残差がゼロに近いかではない（同ファイル:121）。対案は旧 FAIL を保持し、「履歴依存を制御する追加試験」として §6 に事前登録すること。全線・0.02 の基準は維持する。 |
| **Major** | 全線 0.0210 を「反射約 2 %」と記録する／出口端を避ける／境界式を変更する | **現段階ではいずれも却下**。評価器は最終スナップショットの壁圧差を測るだけで、反射波を分離していない（[eval_v2c.py:23](/home/sano/work/forge-sern-design/case/58.farfield_verification/eval_v2c.py:23)、同:43）。記録できるのは「共通初期場試験でも、最終場の全線誤差が基準超過」。超過幅は `0.001 Δp ≈ 87.5 Pa` であり、この幅より十分小さい時間変動の確認が必要。対案は下記の残差 A/B を先行させること。出口延長は後続の原因切り分けにはなり得るが、反射到達前で評価線を切る案は検出対象そのものを外す。 |
| **Major** | `--left-hot` を正式な V2d-2 にする | **右端流出の隔離試験として採用。合格判定は要再検証**。コードは左端だけ高温にし、右端の冷たい自由流を保つため、内外物性差を持つ右端の検証目的は残る（[setup_v2d.py:102](/home/sano/work/forge-sern-design/case/58.farfield_verification/setup_v2d.py:102)）。旧配置の FAIL は別記する。対案は短・長の双方でパルスなし対照を通し、CPG／単成分 TP／多成分 TP について、入射窓で振幅を測り、反射到達窓で短長差を測ること。現評価器の全時間最大値による分母・分子は修正が必要（[eval_v2d.py:61](/home/sano/work/forge-sern-design/case/58.farfield_verification/eval_v2d.py:61)）。独立参照の精度確認も §6 のまま残す。 |
| **Major** | V2 未受入れのまま V3 に進む | **却下**。V3 は V1–V2 合格後という明示条件がある（[plan:187](/home/sano/work/forge-sern-design/plans/active/boundary-node-farfield-characteristic.md:187)）。対案は V2c の定常性・誤差帰属と、再設計した V2d-2 の受入れを先に終えること。幅系列 4 本ではこの未解決事項を代替できない。 |
| **Minor** | 作用素比較器の `VERDICT` を単独で根拠にする | **要修正**。幾何差・状態差・未対応面数は表示されるが、`ok` に入るのは残差差だけ。ゼロ規模も宣言された自由流基準ではなく `1.0` になる（[v2c_operator_ab.py:86](/home/sano/work/forge-sern-design/case/58.farfield_verification/v2c_operator_ab.py:86)、同:103、同:127）。対案は全項目・対応の完全性・有限性を判定に含めること。今回報告された「全項目差 0、未対応なし」が原本で確認できれば、その観測まで否定する理由はない。 |

結論: **次は同じ低領域 B の格子・BC に、`run_0064_v2c_B` と `run_0094_v2c_B_fromC` の最終保存量をそれぞれ与え、初回の絶対残差を比較する A/B で、「残差を残した停止」と「異なる離散定常解の候補」を分ける。**

第 1 仮説: **角の履歴差は、少なくとも一方が非零の離散残差を残す反復停滞であり、二つの定常解を示していない。** 確度: **低**
  
  根拠: ブリーフでは旧系列の残差低下が 2.4–2.5 桁でプラトー、新系列は未確認。`case/58.farfield_verification/run_0092_v2c_opab_low`／`run_0093_v2c_opab_high` の報告値 `|ΔR|/Σ|F| ≤ 7.7e-8` は、両者の `|R|` を制限しない。`eval_v2c.py:25` も最終場しか読まない。**停滞の直接証拠は未確認**であり、優先する理由は二解説を安く反証できるため。

  反証条件: 同一 B 作用素で、異なる両状態とも全保存量の局所残差が下記の診断閾値内なら、「閾値を超える残差を残した停止」という仮説を棄却する。ただし、それだけで厳密な二解・収束を認定しない。

第 2 仮説: **共通初期場でも残る下流壁圧差には、上面 `farfield` の斜入射への応答が含まれる。** 確度: **中**。現仕様も斜入射の限界を明記する（[methods/boundary.md:150](/home/sano/work/forge-sern-design/methods/boundary.md:150)）。ただし **0.0210 の出口処理への帰属、反射率への読み替えは未確認**。

第 3 仮説: **旧 V2d-2 の擾乱源は、左端から入る接触面の TP 保存形混合である。** 確度: **中**。`--left-hot` が接触面を消す操作であることはコードで確認したが、境界固有の寄与は未除外。周期内の段差から「同じ桁の 2.5 Pa」が出るだけでは原因帰属できない。周期配置に伴うもう一つの接触面との干渉、混合率・局所 EOS 復元圧力・波の発生位置を区別する必要がある。

判別 A/B: **変えるのは初期保存量だけ。各 1 step、判定は更新前の初回評価。**

- 格子・BC・数値設定・バイナリ・CUDA block size を `run_0064_v2c_B` と同一にする。
- A 入力＝`run_0064_v2c_B/res_6000.h5`、B 入力＝`run_0094_v2c_B_fromC/res_6000.h5`。いずれも `case/58.farfield_verification/` 配下。同一格子であることを確認して index コピーする。
- 全節点・全 5 保存量の**更新に渡る残差**を比較する。帳簿の `res_final` は壁射影等の後に取得される（[main.cpp:1555](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1555)）。角・出口交線の値も別集計する。
- 診断閾値は事前に `|R_q| / Σ面|F_q| ≤ 1e-5`。分母には境界面も含め、ゼロ規模には単位の合う非零の自由流基準を使う。**これは今回の診断閾値であり、収束基準の代用ではない。**

→ **少なくとも片方が超えるなら**「二つとも離散定常解」という説明を棄却し、非零残差が更新されない箇所を追う。  
→ **両方とも全項目で許容内なら**第 1 仮説をその精度で棄却し、異なる離散平衡の候補として履歴依存を扱う。既存系列の `check_convergence.py` と壁圧時系列の `check_quasisteady.py --series-csv` の判定は、別途必要。

やらない方がよいこと: **出口端の除外で合格化する、0.0210 に合わせて許容値を緩める、壁圧差を反射係数と呼ぶ、未確認の C を収束床に指定すること。** `--from-floor` は参照 run 自身の通常判定 PASS を要求する（[check_convergence.py:359](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_convergence.py:359)）。

呼び出し側の前提への異議: 観測は「初期場を変えると、6000 step 後の角圧力が別の値を保持した」。解釈は「少なくとも二つの定常解」「領域高さが解を選ぶ」「出口端差は反射と出口処理の重なり」であり、後者はまだ確定できない。局所作用素差の除外は、報告された入力状態・対象領域・精度に限定して受け入れ、同じ試験を繰り返す提案はしない。**plan 未反映（依頼どおり変更なし）。**

不足情報: 手元の対象 `run_*` は **0 件**（[case README:5](/home/sano/work/forge-sern-design/case/58.farfield_verification/README.md:5)）。AWS の実効 YAML、バイナリ識別、帳簿原本、全残差履歴、保存場時系列、メッシュ品質・収束・準定常性の VERDICT と判定区間が必要。したがって今回、**run 数値の再現確認も収束認定も行っていない**。
