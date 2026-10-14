# codex 諮問 (diagnose): core-grid-line-reach

- **brief**: [`notes/reviews/briefs/2026-10-11-core-grid-line-reach.md`](../../notes/reviews/briefs/2026-10-11-core-grid-line-reach.md)
- **plan**: [`plans/active/tooling-nozzle-core-grid.md`](../../plans/active/tooling-nozzle-core-grid.md)
- **date**: 2026-10-11
- **commit**: `12869fc7` (feature/nozzle-core-grid)
- **codex**: effort `high`, 2.7 min, rc=0
- **結論**: Gref の B0 の場から、ラインを有効に保ち、`time.deltaT.lineDtDirectionalCap` だけを **50 対 1** に変える各 2 万 step の A/B を行い、まず停滞から抜ける設定かを判別する。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: Gref の B0 の場から、ラインを有効に保ち、`time.deltaT.lineDtDirectionalCap` だけを **50 対 1** に変える各 2 万 step の A/B を行い、まず停滞から抜ける設定かを判別する。

第 1 仮説: 方向別 dt の拡大が、Gref では反復の減衰を悪化させている。ライン結合そのものを外さなくても、dt の拡大を抑えれば残差の停滞を崩せる。確度: 中。

- 根拠: [`setDT_d.cu:167`](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/setDT_d.cu:167) はライン内部面を CFL の最大値から除外し、[同:194](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/setDT_d.cu:194) で `cfl ≥ cfl_all / R` を課す。`R=1` はこの面スペクトルの集計を全面の最大値へ戻す。上限は平均流と SST の両方に効く（[`solver-settings.md:76`](/home/sano/work/forge-coregrid/procedures/solver-settings.md:76)）。
- 観測との関係: B0 と point の比較では dt・ライン結合・熱 Jacobian を同時に変えているため、ライン単独への帰属はできない（[`tooling-nozzle-core-grid.md:591`](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:591)）。現時点では、変更を最も明確に限定できる dt を先に調べる。
- 反証条件: dt が実際に縮んだことを確認しても、下記の残差低下が出なければ、「上限 1 への変更だけで 2 万 step 内に停滞を崩せる」という限定仮説を退ける。dt が原因ではない、とまでは断定しない。

第 2 仮説: ライン結合または `implicitThermalJacobian 5` による反復行列が遅いモードを残している。確度: 低、個別の寄与は未確認。

第 3 仮説: dt に依存する残差処理・拘束・クリッピングなどが方式依存を作っている。確度: 低、今回確認したコードには立証材料がない。

判別 A/B: **A＝上限 50、B＝上限 1**。両腕とも次の条件で固定する。

- 起点は `case/45.isobutane_m6_d155/run_0608_cg_gref_cont_k1/res_60000.h5`。同一メッシュの `restart_field.py` で保存量をビット一致させ、新しい run を二つ作る。
- `lineImplicit 1`・`lineDtDirectional 1`・`implicitThermalJacobian 5`・実効 CFL 4・緩和 0.7・格子・BC・空間離散化・バイナリを固定する。各 20,000 step、1,000 step ごとの場を保持する。速度は測らない。
- 同じ初期場で、実効設定と局所 dt の変更を確認する。dt が変わらなければ、この A/B は仮説を検査できていない。
- 全保存量の `rms_*` と、登録済みの 8 量を追う。`rms_dq_*` の低下を残差低下の代用にしない。

**事前の読み方**:

- **B のみで明瞭に残差が下がる場合**: 10k〜15k と 15k〜20k の両窓で、平均流の非ゼロ残差列それぞれの中央値が A の 1/3 以下なら、dt 拡大の抑制を有望と判定する。この倍率は今回提案する診断閾値であり、収束基準ではない。SST 残差も必ず併記し、増大していれば「平均流の停滞は緩和、連成系は未解決」とする。
- **実効 dt は縮むが低下しない場合**: 上限 1 を「軽い計算で停滞を崩す設定」とする提案は支持されない。第 2 仮説は残るが、確定しない。
- A も同様に下がる、列ごとに結果が割れる、非有限が出る場合は判定不能または発散として区別する。自動延長しない。

これはまず**残差への応答を測る試験**である。量が point 側へ動いたことだけでは到達としない。「同じ結果」の判定を行うなら、8 量の `check_quasisteady` の `STEADY`、窓内変動、参照との差を併せて確認する。今回の識別用には、相対量 0.1 %、軸・η 0.1 の M は 0.005 %pt を暫定の差の上限として事前登録することを提案する。これは生産の許容差の変更ではない。残差については別途 `check_convergence` の VERDICT をそのまま報告する。

やらない方がよいこと:

- **[Major] 上限キーを外して保守側の試験とすること。** 既定 0 は上限なしで、dt の拡大を許す方向である。対案は明示的な `1`（[`solver-settings.md:76`](/home/sano/work/forge-coregrid/procedures/solver-settings.md:76)）。
- **[Major] 短い時系列で point の値を横切ったことを到達とすること。** 既存 point 継続でも θ_r は目標側を越えてから戻っている。対案は全時系列と準定常判定を使う（[`tooling-nozzle-core-grid.md:603`](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:603)）。
- 熱 Jacobian・CFL・緩和を同時に変更すること。今回の A/B では固定する。

呼び出し側の前提への異議:

- **[Major] T1 の「戻る → 右辺が違う」という判別は成立しない。** 同じ R(U) の根でも反復写像が不安定なら離れるし、今回の point 場は厳密な根ではない。plan 記載の VERDICT は、point 通算 0〜200k が **`NOT CONVERGED (プラトー)`**、量の 160k〜200k が **8 量 `STEADY`** である（[`tooling-nozzle-core-grid.md:663`](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:663)）。対案は今回の一設定 A/B。右辺変更の立証には、別途、同一保存量で更新前の残差を直接比較する必要がある。
- **[Major] 格子比較の許容差 1.5 %を、そのまま「line で同じ結果」の判定に使えない。** 既存の B0–point 差は θ_r 0.572〜0.629 %・Q_w 0.514 %で、改善前から通ってしまう。対案は上記の識別用閾値を別に登録する（[`tooling-nozzle-core-grid.md:669`](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:669)）。

不足情報: 実行バイナリと確認したソースの対応、局所 dt の実測、全残差列の窓統計、参照量の実値・漸近値。指定された閲覧制限に従い、生の run データと VERDICT ファイルは再検証していない。数値・VERDICT は plan の記録による。別セッションの「全部入りは止まる」は判断材料に使わない。

plan 未反映。ファイル変更・計算実行は行っていない。呼び出し側で `plans/active/tooling-nozzle-core-grid.md` §4.32 に試験条件を登録する。
