# codex 諮問 (diagnose): line-viscous-heatK

- **brief**: [`notes/reviews/briefs/2026-10-09-line-viscous-heatK.md`](../../notes/reviews/briefs/2026-10-09-line-viscous-heatK.md)
- **plan**: [`plans/active/time_integration-line-viscous-jacobian.md`](../../plans/active/time_integration-line-viscous-jacobian.md)
- **date**: 2026-10-09
- **commit**: `3b73a0d2` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.0 min, rc=0
- **結論**: **既存E1の検査不足を埋めた後、同一条件で「熱伝導Kの密度列だけを残す／除く」A/Bを一組行う。**
- **extra**: `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表:

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 「熱伝導Kを除けば破綻を回避できる」と確定する | **要再検証／Major** | [plan:261](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:261)は事前確認が外れれば無効と規定し、[277行](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:277)では保存場の非有限・非正チェックも未完了。現状は「その解釈を支持する探索的観測」とする。登録上の合格へ読み替えない。 |
| floatのLHSに相対許容1e−12を置いたのは設定ミス | **採用。ただし説明の訂正が必要／Major** | [組立精度の分岐:1743](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1743)からISP 0はfloat。ただし、記載値−4050.89付近のfloat32の1 ULPは **2.44140625e−4**。差4.9e−4は約**2 ULP**であり、「1 ULP」は不正確。原ダンプの成分別差を確認し、演算順序に基づく丸め誤差の許容を次回用に事前登録する。今回の無効判定は残す。 |
| 現行の事前検査を起動ゲートとして使う | **却下／Major** | [e1_compare.py:38](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/e1_compare.py:38)は無効でも終了コードを失敗にせず、[e1.sh:22](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/e1.sh:22)も判定を確認せず計算へ進む。検査失敗・準備失敗で起動を止める。また、現在の検査は拘束フラグの一致や、行4の差が期待する熱伝導式そのものかを検査していないため追加する。 |
| H1の「TPの負のe・基準のずれが原因」という説明 | **現段階では却下／Major** | [熱伝導K:140](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:140)の全列を合わせた作用はκδT。密度列の係数だけの大きさ・符号では原因を判定できない。対案は、密度列の**必要性**とエネルギー基準の問題を分けて検査すること。 |
| マスク5を加速候補として評価する | **条件付き採用／Major** | [plan:276](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-viscous-jacobian.md:276)では2000 step有限だが、最大残差は開始比ρ650倍・ρv1060倍。短期安定・収束・速度改善の証明ではない。下記の監査後に候補として評価し、既定変更・本線採用は保留する。 |

結論: **既存E1の検査不足を埋めた後、同一条件で「熱伝導Kの密度列だけを残す／除く」A/Bを一組行う。**

第 1 仮説: **今回の早期破綻には、熱伝導Kの密度列を通じた流れブロックとの連成が必要である。** 確度: **低、未確認**。H1をこの狭い意味に限定して先に判別する。  
　根拠: [コード:144](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:144)で密度補正が隣のエネルギー行へ入る。plan記録では、`case/45.isobutane_m6_d155/run_0305_e1_m7/`は150 stepで非有限、`run_0306_e1_m5/`は2000 step有限。ただし、この比較は熱伝導Kの全列を除いており、密度列を特定する証拠ではない。先に調べる理由は、変更を一項に限定できるため。  
　反証条件: 密度列の熱伝導寄与を除いても、同じ場所・同じ前駆過程の早期破綻が再現すること。「この密度列が必要」を棄却する。

第 2 仮説: **H2――熱伝導の近傍結合によって温度補正への減衰が弱まり、凍結物性・省略した勾配項・対流を含む実残差との不整合が顕在化する。** 確度: 低。実残差には[LHSが省略する熱伝導の勾配項](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/viscousFlux_d.cu:282)があるが、支配項は未確認。

第 3 仮説: **H3――熱伝導Kと壁拘束との相互作用。** 確度: 低。壁拘束単独の比較では、この相互作用は除外できない。ただし、[等温壁隣のK除去](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:982)と[ΔT_w＝0の拘束](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1147)は代数的には整合しており、先に壁からN層という任意の範囲を導入する根拠はない。

判別 A/B:

- **A:** 値3・マスク7の全熱伝導K。
- **B:** Aから、**熱伝導寄与のK[4][0]だけ**を除く。運動量列・エネルギー列、粘性仕事のK[4][0]、D、壁拘束は残す。候補(a)の密度列と運動量列を同時に除く案より、H1の判別範囲を狭める。
- 同一新バイナリ、同一メッシュ、`run_0183_ns_coldmesh_tw300_ext/res_100000.h5`から保存量を一致させて開始する。両側ともキー5、方向別、上限0、`cfl_pseudo: 4`、緩和0.7、5 sweep、**ISP 0のまま最大2000 step**。新しいrunを使い、逐次実行する。
- 初回の同一状態で、D・拘束行・rhsが一致し、Kの差が指定した熱伝導の密度項だけであることを確認する。rhsの再実行差と、決定的な行列組立差は区別する。
- 全残差を毎step記録する。序盤200 stepは既存帳簿に加え、全域の非有限・ρ/P/T非正・EOS床到達を監視する。密度項、運動量項、エネルギー項それぞれの `K_heat ΔQ` と、その和κδTを、壁拘束適用後の補正で比較する。これは過渡の診断量として扱う。

**→ Aが既知の破綻を再現し、Bが2000 step有限・非物理値なし:** 密度列の除去がこの条件・期間の破綻回避に十分。H1の限定版を支持する。ただし、式の誤り・TP基準の問題・恒久修正の妥当性は証明しない。  
**→ 両側で同じ前駆過程の破綻:** 密度列が必要という仮説を棄却する。H2/H3を確定したとは言わない。  
**→ Aが再現しない、またはBだけ別の破綻を起こす:** 判別不能。旧runを対照に置き換えない。

やらない方がよいこと:

- **列を除いて安定化したことを、その列の微分が誤りという証拠にしない。** 全熱伝導Kは、速度・温度一定の補正δQ＝δρ(1,u,v,w,E)を消す。密度列だけ、またはエネルギー列以外を除くと、この性質を壊す。今回の変更は診断専用とする。
- `FORGE_FREEZE_TURB`を熱伝導率の凍結とみなさない。[SST更新の停止](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:2376)とは別に、[物性](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:1906)と[渦粘性](/home/sano/work/forge-integ-1005/solver_density_cuda/main.cpp:1955)は再計算される。このスイッチだけではH2を除外できない。
- 検査を通すためだけにISP 1へ変えない。元の比較条件を変えてしまう。

呼び出し側の前提への異議:

**［Major］観測と解釈を分ける必要がある。** 記録された観測は「7が非有限、5は残差が有限で2000 step到達」。非物理値の検査と事前ゲートが未完了なので、「分岐1が成立」「運動量・仕事だけでは壊れない」はまだ確定できない。対案は、既存成果物を監査し、登録判定と探索的解釈を別欄に残すこと。

マスク5の長期評価に進む前には、少なくとも全残差列、保存された全時点の`VALUE/*`、ρ/P/Tの正値性、EOS床到達履歴、実効設定・マスク・バイナリSHA・restart一致を確認する。記録上の判定は、`run_0306_e1_m5/`の区間0〜1999で **`NOT CONVERGED (still converging)`**。低下3.0桁を収束と扱わない。

監査後はマスク5を候補として残してよい。ただし性能の結論は、上限なし・上限50・既存方式について、**同じ終了条件までのpoint仕上げを含む総壁時計**で出す。全残差の判定と、対象量の`check_quasisteady`、point切り戻し後の変化を要求する。2000 stepの残差低下だけで長期評価の優先候補にはしない。

不足情報: 対象`run_0305〜0310`の原本は、この作業ツリーおよび確認した`/home/sano/work/forge`にない。残差CSV・HDF5・行列ダンプ・`CONVERGENCE_VERDICT.txt`を独立確認できず、上記run数値はplanの記録を引用した。マスク0の差が丸めだけという説明も、全成分については未検証である。

ファイル変更・forge起動なし。**plan未反映**。呼び出し側で§6.9に判定の留保と採否を、次節にこのA/Bの事前分岐を反映すること。
