# codex レビュー: tooling-nozzle-dual-bell (plan)

- **plan**: [`plans/active/tooling-nozzle-dual-bell.md`](../../plans/active/tooling-nozzle-dual-bell.md)
- **stage**: `plan`
- **date**: 2026-10-11
- **commit**: `07447940` (feature/thruster-nozzle-design)
- **codex**: effort `high`, 4.6 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M8/m0
- **focus**: §4.2 モード判定の契約 (方向別閾値・適格データ・TP 写像)、§4.3 の制約と事前判定、§4.4 の形状 (全長 CP の延長長の従属決定、衝撃波のないベースのみ受理)、§6 V1〜V5 の較正と独立検証の分離。前回 diagnose (notes/reviews/2026-10-11-dual-bell-plan-diagnose.md) の指摘が反映されたかも確認。*.pdf・*.log・plans/README.md は読まない

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
目的と兄弟 plan との分担は妥当ですが、現稿のまま実装に進むことは勧めません。  
前回指摘は概ね反映されています。ただし、保証状態の契約、履歴を含む制約、全長 CP の受理条件には誤受理を防げない穴が残っています。

ファイルは変更していません。指定の禁止ファイルも未読です。関連する `accepted/` と実装を確認しましたが、デュアルベル全体を解決済みの機能は見つかりませんでした。`case/69.dual_bell/` はまだ存在せず、以下の数値は**既存関数・合成入力の実行結果であり、CFD run の実測ではありません**。原典 PDF の内容は再検証していません。

1. **Major — 「較正済み・未検証」が保証判定へどう接続するか未定義**

   **根拠:** [plan:71](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-dual-bell.md:71) は TP の外挿幅未定を「較正済み・未検証」、[V1:121](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-dual-bell.md:121) は同じ条件を「保証未判定」とします。一方、[MOO:92](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-dual-bell.md:92) が明示的に除外するのは「未判定」だけです。「未検証」を別状態として可行解へ通す余地があります。

   **対案:** `validation_status` と `guarantee_verdict` を分離し、必要閾値の欠測・外挿幅未定・検証未受理は必ず `guarantee_verdict=UNKNOWN`、可行解への採用不可と定義してください。全候補が `UNKNOWN` でも `INFEASIBLE` とは判定しません。TP の保証付き MOO は、この証拠ゲートの通過後に置くべきです。

2. **Major — §4.3 はスロットルを含む任意履歴を誤受理する**

   **根拠:** [plan:80](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-dual-bell.md:80) は通常圧で上げ側完了閾値を超え、最低圧で下げ側開始閾値を超えれば高空モードを認めます。合成例 `N_up,end=60`、`N_dn,start=45`、通常 NPR=80、最低 NPR=50 はこの条件を通ります。しかし最低圧のまま剥離状態から上昇した履歴では、全付着へ到達していません。

   また、[plan:76](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-dual-bell.md:76) の「壁圧と剥離圧の比が正」では、`p_w/p_sep=0.5` も通ります。必要な最低燃焼室圧も [YAML:62](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-dual-bell.md:62) の低高度側入力にはありません。

   **対案:** `arbitrary` では、高高度域の最小 NPR 自体に `N_up,end^hi` 到達条件を課してください。下げ側閾値だけでよいのは、先行する全付着到達を保証した履歴に限定します。ベース余裕は `min(p_w/p_sep,upper − 1)>0` と定義し、必要な圧力上下限を入力契約へ追加してください。事前判定も「NPR 比 ≤ 係数比なら不成立、超えても成立未確定」と明記すべきです。

3. **Major — V3/V4 の分離だけでは較正の識別性と独立性を確保できない**

   **根拠:** [plan:110](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-dual-bell.md:110) の台帳には、係数 `c=N_event/N₀` の分母を決める **延長壁面 Mach または壁圧比の出所・不確かさ**がありません。[調査ノート:109](/home/sano/work/forge-thruster/notes/investigations/nozzle-flow-separation-validation-survey.md:109) では P6.2 のベース設計 Mach が不明で、[V2:122](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-dual-bell.md:122) 自身も幾何一致は壁面 Mach の証明ではないと認めています。再構成誤差を係数のばらつきへ混入させる危険があります。

   **対案:** 台帳に `M_e`／`p_ext/p_c`、模型 ID、試験系列、重複掲載の対応、圧力掃引速度を追加し、入力側の不確かさも予測幅へ伝播させてください。較正・検証は点の無作為分割ではなく、少なくとも模型・試験系列単位で分離します。適格性は事象ごとに判定し、「全付着しなかった」試験でも有効な開始点・閾値の下限情報まで一律に捨てないでください。

4. **Major — θ_end の感度記録では「全長 CP」を検証できない**

   **根拠:** [plan:141](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-dual-bell.md:141) は 0.5°で止め、0.25°・1°で遷移判定への影響を記録します。しかし理想 CP で `p_ext` と壁面 Mach が同じなら、[遷移式:48](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-dual-bell.md:48) は長さを含まず、**切断位置を変えても予測 NPR は変わりません**。これでは切り詰めによる相関の適用外を検出できません。

   **対案:** 全長終端の決定と、その数値近似誤差の受理を独立したゲートにしてください。θ→0 への継続・外挿で残存長を評価し、長さと終端状態の許容差を事前登録します。終端を確認できない候補は `UNKNOWN` とします。0.5°は探索停止用の仮値であり、保証適格性の根拠にはできません。

5. **Major — B-spline・δ*・メッシュ変換後に変曲点を保持する契約がない**

   **根拠:** [plan:87](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-dual-bell.md:87) は当てはめと δ* 補正を行いますが、最終物理壁の左右勾配・角度差・終端条件を規定していません。既存の [物理壁処理:354](/home/sano/work/forge-thruster/design/forge_design/geometry/wall_axismach.py:354) は法線オフセット後に全体を滑らかな spline にします。[メッシャ:145](/home/sano/work/forge-thruster/design/forge_design/meshing/mesh2d.py:145) にも変曲点を必須節点とする入力はありません。合成配置では変曲点 `x=6.2` に対して最寄り節点は `6.257426` でした。

   **対案:** ベースと延長を区分表現し、変曲点では位置だけを共有、左右の勾配を別々に拘束してください。メッシュにも変曲点を必須節点として渡します。δ* 補正・当てはめ・変換後に、角度差、終端、寸法、CP 設計からの逸脱を再検査し、その最終形状に保証判定を紐付ける必要があります。

6. **Major — 衝撃波なしの受理ゲートと軸対称 MOC の検証が不足**

   **根拠:** [plan:85](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-dual-bell.md:85) の受理方針は適切ですが、既存の [MOC:29](/home/sano/work/forge-thruster/design/forge_design/geometry/moc_kernel.py:29) は亜音速入力をクリップします。実行すると `M=0.8 → ν=0 → M≈1.0000000000002` となりました。[前進処理:382](/home/sano/work/forge-thruster/design/forge_design/geometry/moc_kernel.py:382) も一部の失敗点を `continue` で捨てます。さらに [同:151](/home/sano/work/forge-thruster/design/forge_design/geometry/moc_kernel.py:151) は、軸上解析極限の改善がスカラー前進経路には入っていないと明記しています。

   **対案:** クリップ前の状態検査、内部欠損、反復未収束、同族特性線の交差・網の反転を明示的な不受理理由にしてください。V2 に軸対称放射源流の厳密解試験と、衝撃波を必要とする壁・亜音速入力の負例試験を追加します。指定値である壁圧の一定性と平面単純波だけでは、軸対称適合式の正しさを検証できません。

7. **Major — 面積比が変わる設計で既存 η を目的にすると推力順位が逆転する**

   **根拠:** [plan:30](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-dual-bell.md:30) は η を目的としますが、既存の [定義:43](/home/sano/work/forge-thruster/design/forge_design/metrics/extract.py:43) は同面積比の理想推力で割っています。全長 CP の出口面積比は形状生成の結果として変わります。既存 `cf_ideal` を使った合成例では、γ=1.4・真空で次の順位逆転がありました。

   | 面積比 | 仮定した η | 得られる C_F |
   |---|---:|---:|
   | 20 | 0.99 | 1.67303 |
   | 50 | 0.97 | 1.67728 |

   **対案:** 推力性能を最大化する目的は、各高度の `C_F=F/(p_c A_t)` に統一し、η は損失の診断量にしてください。併せて「③と共通の寸法 spec」が固定するのは ε_b か ε_e かを明記し、全長 CP の従属決定と拘束が衝突しないよう整理すべきです。

8. **Major — V5 には閾値を数値として受理する手順と基準ケースの検証がない**

   **根拠:** [plan:125](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-dual-bell.md:125) は合否未登録で、開始・完了を判別する位置許容差、NPR の刻み依存、各方向の初期状態確認がありません。兄弟 plan の定常性ゲートだけでは、この部分は埋まりません。また [plan:143](/home/sano/work/forge-thruster/plans/active/tooling-nozzle-dual-bell.md:143) の「外気圧が有限なのでそのまま扱える」は、NPR≈1.25×10⁵ の TP 外部流計算を検証した根拠にはなりません。

   **対案:** 各閾値を「最後に旧状態だった NPR／最初に新状態だった NPR」の区間で返し、位置判定・NPR 刻み・格子・領域の不確かさを登録してください。各端点に同一実効設定区間の `check_convergence.py` の **PASS** と、剥離位置・推力の `check_quasisteady.py` の **STEADY** を必須とします。加えて、MOO 前に基準 TP の h₁・h₂ を固定形状で検証するゲートを §5.1 に置いてください。

**推奨は、保証契約と最終形状の受理を先行させ、TP の保証付き MOO を後段へ分離する段階実装です。** 実装前の修正優先順は、①保証状態 → ②制約式 → ③データ適格性 → ④全長終端 → ⑤変曲点の保持 → ⑥MOC 検証 → ⑦目的関数 → ⑧CFD 受理手順です。ソルバ無変更の範囲では、周期・cell の回帰を追加するより、この node 軸対称経路を直接検証する方が適切です。

読み取り専用の依頼に従い、**plan 未反映**です。

指摘数: Critical 0 / Major 8 / Minor 0
