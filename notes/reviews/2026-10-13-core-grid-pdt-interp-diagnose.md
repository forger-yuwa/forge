# codex 諮問 (diagnose): core-grid-pdt-interp

- **brief**: [`notes/reviews/briefs/2026-10-13-core-grid-pdt-interp.md`](../../notes/reviews/briefs/2026-10-13-core-grid-pdt-interp.md)
- **plan**: [`plans/active/tooling-nozzle-core-grid.md`](../../plans/active/tooling-nozzle-core-grid.md)
- **date**: 2026-10-13
- **commit**: `ef43866d` (feature/nozzle-core-grid)
- **codex**: effort `high`, 2.9 min, rc=0
- **結論**: **候補 (a) を採用し、同一バイナリ・同一起点・`precondEps: 1.0` で `lowMachPrecond: 0` 対 `3` を比較して、Γ が恒等行列でも前処理カーネルが反復を悪化させるかを先に判別する。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **候補 (a) を採用し、同一バイナリ・同一起点・`precondEps: 1.0` で `lowMachPrecond: 0` 対 `3` を比較して、Γ が恒等行列でも前処理カーネルが反復を悪化させるかを先に判別する。**

第 1 仮説: **Γ 以外のカーネル差が、刻みの拡大を止めた S0 の残差低減を妨げている。** 確度: 中  
  根拠: 実装差は確認できた。ただし、このケースへの寄与は未確認である。

- **`Ht` の扱いが違う。** 標準経路は `Ht[ic]` をそのまま使うが、前処理経路は `max(Ht[ic], 1e−8)` に置き換える。この値が対流 Jacobian と軸対称ソース、Γ に渡るため、β = 1 でも差が残り得る。根拠: [timeIntegration_d.cu:843](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/timeIntegration_d.cu:843)、[同:1338](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1338)。負の `roe` の記録だけから、今回の `Ht` も負だとは断定しない。
- **軸対称ソース Jacobian が違う。** `isAxisymmetric == 1` の密度列には、前処理経路だけ `chi_hoop = c² − (γ−1)(Ht−|u|²/2)` が入る。β に依存しない差である。根拠: [同:1113](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1113)、[同:1427](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1427)。どちらの近似が今回の反復に有利かは別途判別が必要。
- **演算精度も同一とは限らない。** 標準経路は `implicitSolvePrecision` で `float`／`double` を選ぶ一方、前処理経路の行列・解法は `flow_float`。FP64 バイナリという名称だけでは両経路の精度一致を保証しない。根拠: [同:1878](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1878)、[同:1352](/home/sano/work/forge-coregrid/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1352)。

  反証条件: β = 1・拡大なしの前処理経路が、標準経路と下記の許容差内で同程度に残差を減らすなら、「Γ が恒等でも残るカーネル差が大きな悪化を生む」という仮説を、この設定・予算内で退ける。β < 1 との相互作用まで除外したことにはしない。

第 2・第 3 仮説:  
  第 2: β < 1 の時間項と、有限回の block 反復・境界拘束・SST の分離更新との組み合わせが減衰を悪くしている。確度: 低、未確認。  
  第 3: 追加しない。既存の境界凍結を、根拠なく再び真因に据えない。

判別 A/B: **両腕共通で `time.deltaT.precondEps: 1.0`、変更する設定は `time.deltaT.lowMachPrecond` の 0／3 だけ。**

- 同じ診断バイナリ、`run_0608` の `res_60000`、point、`cfl_pseudo: 1`、`implicitRelax: 0.7`、`nStepInner: 5`。両腕とも `FORGE_DIAG_PRECOND_NO_DT_SCALE=1` とする。精度設定や境界設定は同時に変えない。既存 C は背景資料とし、対照も新しく作る。
- 各 **20,000 step、1,000 step ごとに記録**。主指標は縮流部の ρ・ρE の q。窓は従来どおり W1 = 10k〜14k、W2 = 16k〜20k。`throatwall` と他の保存量は悪化の移動を監視する。
- 資格は、初期保存量のビット一致、更新前の全残差の相対差 ≤ 1e−6、初期の局所刻みの一致、前処理経路の β = 1、非有限なし。`precondEps` は実効値を確認する。初回更新後の差は測定対象なので、一致を資格にしない。
- **結果 A:** 両窓・両主成分で q̄₃/q̄₀ ≥ 2 → Γ が恒等でも悪化するため、「残る悪化は β < 1 の時間項だけが原因」を却下し、カーネル差を優先する。どの差が原因かまでは決めない。
- **結果 B:** 両窓・両主成分の q̄₃/q̄₀ と減衰率比 d₃/d₀ がともに 0.9〜1.1、かつ他領域へ明確な悪化が移らない → 第 1 仮説を上記の範囲で退け、β < 1 を必要とする第 2 仮説を優先する。
- 中間・成分間の不一致・資格不足は判定不能。これらは**今回提案する診断閾値**であり、収束基準ではない。両腕の `check_convergence.py` の VERDICT と判定区間を記録する。

やらない方がよいこと: `precondEps` の複数値や CFL を先に振ること、`Ht`・ソース項・精度をまとめて修正すること、今回の結果だけで低マッハの仮説全体を閉じること。まず上の比較を一つだけ行う。

呼び出し側の前提への異議:

- **[Major・採用] §4.87 の「判定不能」は維持する。** S0 の d は `throatwall` で 0.95／1.02、縮流部で 1.05／1.03。登録した減衰条件を満たさない。対案は「この窓では拡大を止めると増加率が小さくなった」と観測を限定し、成功判定へ読み替えないこと。根拠: [plan §4.87](/home/sano/work/forge-coregrid/plans/active/tooling-nozzle-core-grid.md:1239)。
- **[Major・却下] S0 と C の差を Γ 単独の効果とする解釈。** 同じ離散方程式を解くことは、同じ反復作用素であることを意味しない。`case/45.isobutane_m6_d155/run_0693_pdt_off/` の q = 817 と `case/45.isobutane_m6_d155/run_0688_pp_c/` の 49 は、step 20000 の未収束残差の比較である。対案は上記の β = 1 比較。
- **[Major・要再検証] 「β = 1 なら標準カーネルに戻る」という前提。** `lowMachUr` の式から Γ = I になることは確認できるが、上述の実装差は消えない。対案は、β = 1 を同値性の保証ではなく、Γ の寄与を消す診断条件として使うこと。
- **[Minor・表現修正] 「増加は消えた」「横ばい」を全区間へ広げない。** S0 の縮流部の ρ の q は 477 → 679 → 817。対案は「末尾二窓の比は 1.05」と書くこと。§4.87 記録の VERDICT は、S+ が **NOT CONVERGED（`rms_ro`・`rms_roe` が RISING）**、S0 も **NOT CONVERGED（全列 falling、`rms_ro` 0.9 桁）**である。全域残差の下降は局所残差の下降を保証しない。

元の問いとの関係では、今回調べているのは固定した Gref 上の **point 反復の悪化**であり、流れ方向の格子幅もライン結合も変えていない。したがって「流れ方向に細かいとライン陰解法が遅い機構」はまだ特定できない。今回の比較は、その説明に低マッハ前処理の失敗を使う前に、専用カーネルの差を切り分ける位置づけである。

不足情報: 実行バイナリでの `flow_float` と線形 solve の実効精度、`thermalMethod`・`axisymMethod`・壁条件の実効値、初期場で `Ht` の下限が発動する節点数・位置が未確認。禁止された生データは読まず、数値と VERDICT は提示資料の記録を根拠とした。**ファイル変更・forge 起動は行っていない。plan 未反映であり、事前登録への反映は呼び出し側が行う。**
