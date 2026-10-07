# codex 諮問 (diagnose): u4-v5prime-result-interpretation

- **brief**: [`notes/reviews/briefs/2026-10-07-u4-v5prime-result-interpretation.md`](../../notes/reviews/briefs/2026-10-07-u4-v5prime-result-interpretation.md)
- **plan**: [`plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md`](../../plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md)
- **date**: 2026-10-07
- **commit**: `ad29c8fd` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.7 min, rc=0
- **結論**: **N0 の凝縮を同一条件で一度だけ延長し、対照の未達を判定してから、主張を限定した版を result 段レビューへ進める。**
- **extra**: `plans/active/discretization-moc-axis-limit-and-corrector.md`, `case/45.isobutane_m6_d155/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| H1「上流の多項式化は測れるほど量を動かさない」 | **要再検証／Major** | dry の登録 4 量の差が窓幅以下であることは確認できる。しかし凝縮の g は差 −2.256e−7、窓幅最大 1.957e−7 で、N0 自身も残差ゲート未達。[評価 JSON:1985](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/ns_n012_eval.json:1985)。**「N1 は自身の dry・凝縮ゲートに合格、N0 との差は参考値」**に限定する。幾何の利点は採用理由にできるが、流れの同等性を証明したとは言えない。 |
| H2「V5′ は改善し、悪化側もゲート内」 | **数値差の記録は採用、性能改善の確定は却下／Major** | g −3.153e−6、M +1.545e−4、時間幅との比 17・22 倍は再計算できる。ただし時間幅は信頼区間ではない。[評価器:457](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/ns_n012_eval.py:457)。また凝縮 K は **4 量 STEADY・残差 RISING なし**であり、開始位置や S_max の悪化幅の上限はない。[monotone plan:268](/home/sano/work/forge-integ-1005/plans/accepted/tooling-nozzle-throat-monotone-r2.md:268)。対案は「N2 は登録ゲート合格、差の方向と大きさは参考値」。 |
| H3「N0 は対照なので延長不要、差の解釈にも効かない」 | **却下／Major** | N0 の未達は N1・N2 自身の合否や N2−N1 を直接無効にはしないが、**N1−N0 による U4 の解釈には効く**。さらに §6 K は「未達は延長 1 回」と明記し、[評価 JSON:762](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/ns_n012_eval.json:762) も同じ指示を保持している。対案は N0 凝縮だけを一度延長すること。 |
| V5′ 合格をもって MOC の生産採用条件が揃ったとする | **要再検証／Major** | V5d と出口較正判断は依然「保留（前提不成立）」。腕 M の全温の時間幅 0.113〜0.126 K は条件 0.1 K を超える。[MOC plan:251](/home/sano/work/forge-integ-1005/plans/active/discretization-moc-axis-limit-and-corrector.md:251)。NS の合格で Euler の未達を上書きしない。採用を諮る際は、この保留を残した限定採用・要件変更であることを明示する。 |
| 凝縮開始位置の差 0／−0.0539 r_t の読み方 | **訂正を採用／Minor** | 抽出は軸方向 2401 点のうち g > 1e−4 となる最初の点で、閾値交差の補間ではない。[cond_series.py:11](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/cond_series.py:11)。N2 の窓幅 0.0448865 r_t は抽出刻みと整合する。対案は「同じ抽出点」「約 1 刻み上流」と記録し、物理的な差ゼロや高精度の開始位置変化とは扱わない。 |

結論: **N0 の凝縮を同一条件で一度だけ延長し、対照の未達を判定してから、主張を限定した版を result 段レビューへ進める。**

第 1 仮説: N0 の `RISING` は、凝縮モーメントの遅い緩和を現在の判定窓が捉えている可能性がある。確度: **中**。

根拠: 対象は `case/45.isobutane_m6_d155/run_0168_ns_n012_N0_cond/`。保存された判定は次のとおり。

- 凝縮 14000〜18000 step：4 量とも **`STEADY`**。
- 残差の本段区間：**`NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)`**。ただし列別判定には **`RISING`** があり、評価器は正しく未達としている。[評価 JSON:752](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/_band_ab/ns_n012_eval.json:752)。
- plan に記録された `rms_rog_0` のブロック平均は 1.26→1.39→1.45→1.51e−10。これは急増とは異なるが、無害な丸め誤差の床と認定する証拠でもない。[upstream plan:288](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-upstream-poly-and-throat-sizing.md:288)。

判定器は単なる前後平均の大小に加え、緩い持続上昇も調べる。「絶対値が小さい」「N1・N2 と同じ帯」だけでは判定を退けられない。[check_convergence.py:60](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_convergence.py:60)。

反証条件: 下記の固定予算の延長でも `RISING` が残る、または凝縮 4 量が非定常になる場合、「今回の延長で解消する緩和過渡」という限定仮説を棄却する。真因の同定とは区別する。

第 2 仮説: 残差の床付近の揺れが、判定窓によって上昇として現れている。確度: **低・未確認**。同じ残差水準という情報だけでは、遅い過渡と区別できない。

判別 A/B: **変える要因は N0 凝縮の計算長だけ**。A＝既存 18000 step、B＝同条件で追加 20000 step、通算 38000 step。

- 追加 20000 step は今回の提案値。§6 K は延長回数を定めるが、長さは明記していないため、投入前に登録する。
- `run_0168` の `res_18000.h5` から、新しい `run_NNNN_ns_n012_N0_cond_ext/` へ `restart_field.py` で保存量を継ぐ。凝縮を再初期化せず、バイナリ・メッシュ・CFL・物理設定を固定する。
- 1000 step ごとに保存し、凝縮 4 量は通算 **34000〜38000 の 5 枚**で判定。残差は親子の同一設定区間を明示して連結し、全列を判定する。
- **B で RISING が消え、4 量も STEADY なら**、今回の延長では未達が解消しないという仮説を棄却し、K の合格として記録する。残差が plateau なら「収束」とは書かない。
- **B でも未達なら**、追加 20000 step で解消できるという仮説を棄却し、再延長せず未達のままユーザ判断へ戻す。

これは N0 の未達を診断する試験である。延長後の N0 と既存 18000 step の N1 の差を、同じ判定窓による比較として扱わない。

やらない方がよいこと: 残差の絶対値を理由に N0 を合格へ変更すること、時間幅との比を有意差・非劣化の証明にすること、V5d を参考値の良さで合格へ読み替えること、今回の診断と同時に較正値・k_f・CFL を変更すること。

呼び出し側の前提への異議:

- **ユーザ決定の 2 件は有効な変更として受け入れる。** dry は通算 80000〜100000 の窓で新しいゲートに合格している。ただしオーバーシュートの元の VERDICT は 3 条件とも **`DRIFTING`**、残差は **`NOT CONVERGED`** のままである。
- U4 は、N1 自身の合格と幾何の利点を根拠にした限定採用の材料がある。しかし「凝縮への影響なし」を根拠にした採用諮問には不足がある。
- V5′ は、N1・N2 の凝縮 4 量がともに **`STEADY`**、残差はともに **`NOT CONVERGED (stalled/plateau)`・RISING なし**という登録上の合格を支持する。**MOC 全体の生産採用には V5d の保留の扱いと result レビューが残る。**「全検証合格」としては諮れないが、保留を明示した要件変更の判断を諮ることはできる。

不足情報: 対象 run の HDF5・残差 CSV・全時系列 CSV・元の列別 VERDICT はこの作業ツリーに無く、独立再判定はできなかった。確認できたのは評価 JSON、評価コード、plan の記録であり、`ns_n012_eval.py` と `moc_v5d_eval.py` の SHA256 は各保存記録と一致した。result 段には、元の列別判定・実効設定と区間の証拠、および採用予定の生産入力が検証した壁を再現する記録を揃えること。

run の恒久索引は [case README](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:113)。ファイル変更・forge 起動は行っていない。**本提案は plan 未反映**で、採否と延長条件は呼び出し側が upstream plan §6 U4・§9、MOC plan §6・§9 に記録する。
