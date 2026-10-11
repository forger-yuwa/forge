# codex 諮問 (diagnose): m64-and-float

- **brief**: [`notes/reviews/briefs/2026-10-10-m64-and-float.md`](../../notes/reviews/briefs/2026-10-10-m64-and-float.md)
- **plan**: [`plans/active/time_integration-line-implicit-speed.md`](../../plans/active/time_integration-line-implicit-speed.md)
- **date**: 2026-10-10
- **commit**: `71fc9e92` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.8 min, rc=0
- **結論**: **float 化より先に、既存 M64 の終盤区間を対照とする全長ライン B0 の再生を追加 1 本だけ行い、R3 が反復方式による持続的な差かを判別する。本番は当面 B0 を維持する。**
- **extra**: `notes/reviews/2026-10-10-line-implicit-speed-results-diagnose.md`, `methods/time_integration/implementation.md`, `case/45.isobutane_m6_d155/tt_judge.py`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **float 化より先に、既存 M64 の終盤区間を対照とする全長ライン B0 の再生を追加 1 本だけ行い、R3 が反復方式による持続的な差かを判別する。本番は当面 B0 を維持する。**

第 1 仮説: M64 はライン外に戻した領域の平均流の反復を変え、B0 と異なる残差床・到達状態を作っている。**確度: 中**

- 根拠: 提示された `case/45.isobutane_m6_d155/run_0378_tt_M64/` の集計では、Σ|res_ro|·2π は 0.493、B0 の `run_0252_ns_coldmesh_tw300_linedir_tj5_cap50_ext2/res_60000.h5` は 2.417 kg/s。j=20〜56 の寄与も 16 % 対 80 % と異なる。ただし時点も異なるため、これだけで床の因果関係は確定しない。
- コード上、ライン面は近傍補正の lag 参照を外し、Thomas で解く。ライン外は point solve に戻る。この変更は実際に効く。根拠: [timeIntegration_d.cu:889](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:889)、[同:1171](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1171)。
- **「point に近づいたから改善した」は支持しない。** R は未収束の参照状態であり、正解ではない。記録上の M64 の `check_convergence` は **NOT CONVERGED**。M64 の到達窓について `check_quasisteady` の VERDICT は提示されていない。
- 反証条件: 共通の終盤保存場から B0 を再生しても、従来の B0 側へ残差分布・ω残差が戻らず、M64 と同程度の状態を維持するなら、「B0 固有の持続的な床が R3 の主因」という説明を棄却する。

第 2 仮説: ω の差は、平均流の変化と共用する Δτ を通じた SST の応答である。**確度: 中、寄与の内訳は未確認。** SST は B0 でも M64 でも point-implicit 更新であり、「ω も Thomas から point に戻った」という説明は成立しない。[implementation.md:382](/home/sano/work/forge-integ-1005/methods/time_integration/implementation.md:382)

第 3 仮説: 状態 float では、幾何を修復しても残差評価の丸め・commit の消失が到達限界になる。**確度: 中、case/45 での支配性は未確認。** |dq| < ½ ULP(Q) の更新が消える機構は既知で、現在の `qAccumulatorFP64` は軸対称に対応しない。[implementation.md:219](/home/sano/work/forge-integ-1005/methods/time_integration/implementation.md:219)、[同:245](/home/sano/work/forge-integ-1005/methods/time_integration/implementation.md:245)

判別 A/B: **既存 M64 の通算 100000〜140000 step を A とし、その起点 `res_100000.h5` から全長ライン B0 を 40000 step 再生したものを B とする。追加計算は B の 1 本だけ。**

- 変更点は `FORGE_LINE_MAXLEN=64` の解除だけ。同一バイナリ・同一メッシュ・同一保存量から開始し、sweep 5、緩和 0.7、キー 5、方向別 Δτ・上限 50 を固定する。同一メッシュの引き継ぎは `restart_field.py` を使う。
- 出力間隔は 2500 step、判定窓は両者の通算 120000〜140000。B の費用は提示単価で約 **21 分＋起動・出力**。これは終盤の状態差を調べる試験であり、既存の総時間比較を置き換えない。
- 見る量は全 `rms_*`、同じ j 区分での Σ|res_ro|、θ_r の 3 断面、Q_w。ω の解釈には、ω残差の局在と k・ω・μt/μ の局所変化も必要。
- **結果 A**: B0 再生で j=20〜56 の残差寄与が M64 の 2 倍超へ増え、ω残差中央値が M64 の ½ 未満へ下がる傾向が末尾窓で持続する  
  → 第 1 仮説を支持し、「従来の異なる履歴だけで差が出た」という説明を退ける。
- **結果 B**: 全残差中央値・帯別残差が M64 の ½〜2 倍以内、θ_r 差 ≤0.05 %、Q_w 差 ≤0.1 % を維持し、対象量の準定常判定も満たす  
  → 「B0 固有の床」という強い説明を棄却し、履歴・到達時点への依存を候補に戻す。
- 中間結果や `DRIFTING` は判別不能とする。この診断用閾値は事前登録し、既存 R3 の採否を遡って変更しない。`check_convergence` と `check_quasisteady` は比較窓を明示して実行する。
- **起点ファイルの AWS 上での残存確認が必要。** 残っていなければ、この「追加 1 本で共通起点」の試験は成立しない。最終場からの切替試験を同等の A/B と扱わない。

やらない方がよいこと:

- M64 を「point に近い」「SU2 許容差 3 % より小さい」だけで本番採用すること。**診断候補には残すが、R3 解除は保留**する。
- M48 や M64+S3 を今すぐ長時間回すこと。M48 の単価改善は M64 比では約 **6.6 %**で、到達 step が約 **7 %**増えると消える。S2 の発散は S3 の発散を証明しないが、M64+S3 の長期安定性・総時間は未評価。まず現在の品質フラグを処理する。
- float の 1000 step 完走や θ_r の小さなドリフトを、FP64 と同等の到達性能の証拠にすること。更新が丸めで消えた状態も、小さなドリフトを作る。
- `qAccumulatorFP64` の軸対称制限を外すだけで試すこと。

呼び出し側の前提への異議:

| 重大度 | 採否・根拠 | 対案 |
|---|---|---|
| **Major** | **却下:** 「R に近づいたので R3 は悪化ではない」。R を正解扱いしないというユーザ決定と矛盾する。[plan:419](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:419) | R3 は品質の悪化を証明するものでもない。診断フラグとして維持し、上記の共通起点比較で扱う。 |
| **Major** | **要再検証:** H2 の「point に戻したため Δτ が小さくなって遅れた」。コードはライン面を CFL の最大値から除外するが、その面が律速でなければ Δτ は変わらない。plan 自身も縮流部の対象領域では point と同じと記載する。[setDT_d.cu:158](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/setDT_d.cu:158)、[plan:489](/home/sano/work/forge-integ-1005/plans/active/time_integration-line-implicit-speed.md:489) | 同じ状態から算出した Δτ 比と律速面を位置別に確認し、変化した領域だけに帰属させる。 |
| **Major** | **却下:** 「混成 typedef だけで幾何 double の試験になる」。幾何と状態は同じ `var.c_d` / `var.p_d` 経由で渡され、呼出先では `geom_float*` と `flow_float*` を要求する。さらに座標差は **各座標を ST に落としてから**計算し、ISP 0 は ST=float。[setDT_d.cu:299](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/setDT_d.cu:299)、[timeIntegration_d.cu:952](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:952)、[同:1746](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1746) | float 化に着手するなら **案 (b) を先行**する。絶対座標の差・必要な幾何係数を double で作り、差を取った後に float へ渡す。入力・格納・利用箇所まで追跡し、第一層距離と幾何閉包を検査する。 |
| **Minor** | **要再検証:** H3 の「Thomas は 12.6 ms のまま」。内部演算は double でも、D・K・rhs の読込みと dq の保存は `flow_float` で、型変更により転送量・変換が変わる。[timeIntegration_d.cu:2312](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2312)、[同:2372](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/timeIntegration_d.cu:2372) | 19.4→7.9 ms という差引きは内訳の実測として扱わず、混成版で測る。 |

float の θ_r・Q_w が用途上の許容内に入る可能性はあるが、**FP64 と同じ残差床になる見込みは現状の資料から判断できない**。案 (b) でも、float にした面ベクトル・体積の整合性と状態更新の精度は別に残る。

将来の float 試験で最低限必要なのは、幾何係数の FP64 比、同じ保存場での全残差、各保存量の |dq|/ULP(Q) と実際に反映された更新量、そして末尾 2 万 step の θ_r・Q_w・全残差の判定である。短い試験は不適合の早期検出に使えるが、同等性の合格には使えない。commit が支配的なら、代案は **保存量・commit を double に保ち、費用の大きい演算だけを float 化する方式**になる。

不足情報:

- AWS の元 run、実効設定、判定 JSON は手元になく、実測値は提示記録に基づく。M64 の `check_quasisteady`、ω残差の空間分布、再生起点の残存が未確認。
- 設計チェーンが θ_r・δ_r に要求する許容差。差 0.1 % を保証された誤差限界と扱う根拠もない。
- R1 の「速い」は登録規則では成立するが、旧 B0 の到達 step 再利用と専有単価による**推定総時間**である。分解能の和を差し引いた余裕は約 133 秒で、再現ばらつきは含まれていない。

ファイル変更・forge 起動は行っていない。**plan 未反映**。呼び出し側で `plans/active/time_integration-line-implicit-speed.md` §5.1 #20・#23、§6.20 に採否と試験条件を反映すること。
