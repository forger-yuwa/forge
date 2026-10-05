# codex 諮問 (diagnose): sern-r5h-m10-ab

- **brief**: [`notes/reviews/briefs/2026-10-05-sern-r5h-m10-ab.md`](../../notes/reviews/briefs/2026-10-05-sern-r5h-m10-ab.md)
- **plan**: [`plans/active/tooling-nozzle-sern-chain.md`](../../plans/active/tooling-nozzle-sern-chain.md)
- **date**: 2026-10-05
- **commit**: `2347ebc6` (feature/sern-design)
- **codex**: effort `high`, 5.0 min, rc=0
- **結論**: **2面介入や20000 step継続へ進む前に、既存の1面A/Bについて「残差→実保存量更新→次回EOS」の帳簿を時間整合させ、517199が冷える過程を100 step以内で切り分ける。**
- **extra**: `plans/active/tooling-nozzle-sern-3d.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **2面介入や20000 step継続へ進む前に、既存の1面A/Bについて「残差→実保存量更新→次回EOS」の帳簿を時間整合させ、517199が冷える過程を100 step以内で切り分ける。**

第 1 仮説: **517160側の速度再構成は同点の低温維持に寄与し、介入後の517199の低温化は、同じ外部流側の近壁層どうしの連成応答である。** 確度: **中**。

根拠: 報告された`run_1061/1062`では、内部エネルギー射影が−0.801→+10.745、対流の`res_roe`が−6.17→+5.02に変わる。介入は実装上も指定面・片側の速度だけで、流束は両CVへ逆符号で加算される。[介入実装:279](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:279)、[残差加算:661](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:661)、[run索引:379](/home/sano/work/forge-sern-design/case/46.sern_design/README.md:379)

ただし、**残差への寄与は支持されても、実更新とEOSまで含めた機序は未確認**。517199を冷やした面・処理も未特定である。

反証条件: 時間を揃えた保存量から、介入による517160の内部エネルギー増加が確認できないこと。または517199の冷却が、対流を含む更新では説明できず、EOS・境界等の状態上書きで生じていること。

第 2 仮説: **陰解法の増分・組成変化・EOS射影の連成が、517199の床到達を生む。** 確度: **中・未確認**。今回のEOS経路は温度反転後に`roe`を書き換えるため、旧m6_on試験での除外を持ち越せない。[EOS:208](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/dependentVariables_d.cu:208)

第 3 仮説: 追加しない。

判別 A/B: **変更する設定は、従来どおり`FORGE_DIAG_FACE_VEL_CELL=1535833:517160`の有無だけ。**

- 同じ`run_1055`最終保存量を起点に、A＝介入なし、B＝既存1面介入。バイナリ・CFL・`nStepInner`等は揃える。既存帳簿で足りる部分は再解析し、不足分だけ再取得する。
- **まず各2 step**で、初回残差、その残差による保存量のcommit、次回EOS前後を対応させる。続けても**各100 stepまで**とし、517199が初めて床へ到達する更新を捕捉する。
- 対象は517160・517199・517200、介入面の相手522113と接続先。面流束和、段別残差、実保存量増分、`roY*`・`roK`、EOS・壁／境界処理の変更を分ける。内部エネルギー射影と有限の実増分は別々に評価する。
- 閉合許容は既登録の**介入差の1%または丸め誤差上限の大きい方**を維持する。

**→ 結果A:** 更新収支が閉じ、517160の加熱と517199の冷却がEOS前の更新で説明されるなら、近壁層間の連成を支持し、EOS射影だけで冷点が生じる説明を棄却する。  
**→ 結果B:** 冷却が更新では説明されず、EOS等の上書きで発生するなら、対流再構成だけで説明する仮説を棄却する。閉合不足・寄与を分離できない場合は判定不能とする。

やらない方がよいこと: **現時点での20000 step継続、2面・後縁列全体への介入拡大、有限厚後縁への直行。** 初回更新の検証が不足しており、介入範囲を増やしても原因の切り分けが進まない。床ゲート除外や診断用ID指定の生産採用も行わない。

呼び出し側の前提への異議:

1. **Major — 事前登録の総合判定は「判定不能」。局所残差への効果だけ採用する。**  
   100 stepでは床離脱条件を満たさず、準定常後の床残存も確認していない。「機序支持として合格」「1面では永久に解消しない」のどちらも未成立。前回登録は、過渡未終了を明示的に判定不能としている。対案は上記の更新収支補完。[事前登録](/home/sano/work/forge-sern-design/notes/reviews/2026-10-05-sern-r7b-m10-te-floor-diagnose.md)

2. **Major — 「1 step後もTが同じ」を実更新の無応答と読むのは却下する。**  
   コードの順序は、EOS・温度計算→残差構築→陰解法commit→出力である。内部節点の温度はcommit後に再計算されないため、初回保存場のTは更新前の値になり得る。また同じcallの`entry→after_eos_bc`は、そのcallでこれから計算する残差への応答ではない。対案は、**call 1の残差をcall 2の`entry`保存量と結び、次回EOS前後を別に比較すること**。[前処理:1878](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1878)、[commit:2363](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:2363)、[出力:2514](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:2514)

3. **Major — 「板の反対側の双子」という解釈は却下する。**  
   公開設定では`NJ=127`、`nz=39`、カウル線は`j=54`。番号式`(i×NJ+j)×nz+k`を逆算すると、517160＝`(104,52,20)`、517199＝`(104,53,20)`、後縁壁点は517238＝`(104,54,20)`となる。**両点ともカウル下側で、第2層から第1層への移動に対応する。** 後縁では上下コピーも共有される。実メッシュの再番号付け有無・境界所属による最終確認は必要だが、双子を前提とした「対応面への追加介入」は根拠がない。[設定:124](/home/sano/work/forge-sern-design/case/46.sern_design/problem_3d_prod_3op_wallres_lswx08.yaml:124)、[番号・共有規則:265](/home/sano/work/forge-sern-design/design/forge_design/meshing/mesh_sern3d.py:265)

4. **Major — 壁解像の「約4倍の悪化」は、g3とg4の列の取り違えである。**  
   §4.42のramp **0.549**・cowl_in **1.234**は**g4**。同表の**g3**は **2.167・4.919**で、今回のg3 **2.05・4.69**に近い。対案は比較先をg3へ訂正し、局所解像不足は別途保持すること。[§4.42:1704](/home/sano/work/forge-sern-design/plans/active/tooling-nozzle-sern-3d.md:1704)  
   **R7b④の順序を、この誤比較を理由に変更する必要はない。** m10_on診断を先に完了し、4-3以降の本投入は保留する。後続の壁解像評価・緩和判断には、現在条件の格子感度が必要である。

5. **Major — 壁解像ツールの「面積割合」は、その実装では個数割合である。**  
   `good.sum()/len(yp)`と`count_nonzero(yp>target)/good.sum()`で計算しており、面積重みがない。したがって57%・97%等を超過**面積**として使えない。対案は、対応する壁面積で重み付けして再集計し、現行値は標本個数割合と明記すること。平均値だけで壁解像合格へ読み替えない。[集計:363](/home/sano/work/forge-sern-design/solver_density_cuda/tools/check_wall_resolution.py:363)

不足情報: **対象run原本はこのcheckoutにない。** 今回の数値はブリーフ・run索引の報告値で、独立再計算した値ではない。必要なのは`run_1061–1064`のcall番号付き帳簿・保存量時系列・実効設定とバイナリ来歴、実メッシュの境界所属、収束／準定常VERDICT原本である。今回、収束・定常性を新たに認定していない。

ファイル変更・forge起動なし。**plan未反映**。呼び出し側の反映先は`tooling-nozzle-sern-chain.md` §5.1 R7b 4-2と`tooling-nozzle-sern-3d.md` §5.1 R5h。
