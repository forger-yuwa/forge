# codex 諮問 (diagnose): limiter-inlet-next-ab

- **brief**: [`notes/reviews/briefs/2026-10-03-limiter-inlet-next-ab.md`](../../notes/reviews/briefs/2026-10-03-limiter-inlet-next-ab.md)
- **plan**: [`plans/active/limiter-inlet-column-oscillation.md`](../../plans/active/limiter-inlet-column-oscillation.md)
- **date**: 2026-10-03
- **commit**: `f45c1b54` (feature/species-transport)
- **codex**: effort `high`, 4.6 min, rc=0
- **結論**: **新しい計算を投入する前に、既存の毎 step 出力で「ρ を同じ step から読む／一つ前の step から読む」だけを変えた再計算 A/B を行い、ρ の不一致を解消してください。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 対象・採否 | 根拠と対案 |
|---|---|---|
| **Major** | 「ρ の不一致は保留して先へ進む」：**却下** | `face_probe.py` は同じ保存ファイルの ρ・勾配・ψ を組み合わせています（[face_probe.py:28](/home/sano/work/forge-species/notes/investigations/2026-10-03-limiter-inlet/face_probe.py:28)）。しかし定常陰解法は残差評価後に保存量を更新して出力します（[main.cpp:2083](/home/sano/work/forge-species/solver_density_cuda/main.cpp:2083)、同:2130、2296）。**まず ρ の参照時点を揃えて再計算してください。** |
| **Major** | 「typedef だけ変更して全域 FP64」：**却下** | `prim` 経路には型固定の `float4` が残ります（[limiter_d.cu:329](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/limiter_d.cu:329)）。生成側も `float4` で書きます（[calcGradient_d.cu:750](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/calcGradient_d.cu:750)）。`flow_float=double` では配列の歩幅と書込み形式が食い違い、例えば limiter の `prim[ic1*8+4]` は生成側の P の格納位置と一致しません。**この経路が有効なら比較自体が壊れます。** 精度比較には型固定箇所の監査が先です。 |
| **Major** | 「K=0.2 で H-R を判別」：**却下** | K 変更は ε̂ を 8 倍、ε̂² を **64 倍**にし、同じ入力でも ψ を変えます（[limiterFunctions_d.cuh:43](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/limiterFunctions_d.cuh:43)）。改善しても、丸め起因と有限精度によらない非線形振動を分離できません。**K 感度試験という解釈に限定**すべきです。 |
| **Major** | 「リミッタ関数だけ double」：**因果試験として却下** | `face_probe.py` は座標・ε̂²を float64 にし、面増分と有理式を既に主として double で再計算しています（同:13、29、38–42）。報告どおり Ux/Uy/P が 2e−7 以内なら、**その入力に対する関数評価の丸めだけで ψ の大振幅を説明することは困難**です。入力の float32 量子化は double 評価では戻りません。指定コードに関数だけ精度を切り替える設定もありません。 |
| **Major** | 「1/5 以下＝主因、2 倍以内＝否定」を他の介入へ転用：**却下** | 元の基準は CFL 半減への応答です（[plan:104](/home/sano/work/forge-species/plans/active/limiter-inlet-column-oscillation.md:104)）。効果量の基準としては使えても、**H-R の採否基準にはなりません**。また「2 倍以内」は `0.5 ≤ 比 ≤ 2` と明記しないと、1/5 以下とも重なります。 |

結論: **新しい計算を投入する前に、既存の毎 step 出力で「ρ を同じ step から読む／一つ前の step から読む」だけを変えた再計算 A/B を行い、ρ の不一致を解消してください。**

第 1 仮説: **ρ の最大差 1 は、更新後の密度と更新前の勾配・ψ を組み合わせた後処理の時点不一致による。** 確度: **中**。

  根拠:
  
- `prim` の先頭は **ρ そのもの**です。生成側は `[ro,Ux,Uy,Uz,P,T,0,0]` を格納します。逆密度や別の変数ではありません（[calcGradient_d.cu:744](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/calcGradient_d.cu:744)）。
- limiter は owner の `Q0=ro` と、近傍の `prim` 先頭を比較します（[limiter_d.cu:319](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/limiter_d.cu:319)、同:329、480）。ただし `prim` 使用条件は `primPack!=0 && gradLSQ==2` で、対象 run の実効設定は未確認です。
- 原始量評価→BC→勾配→limiter→保存量更新→出力という順序です（[main.cpp:1828](/home/sano/work/forge-species/solver_density_cuda/main.cpp:1828)、同:1847、1866、2130、2296）。ρ は更新対象の保存量でもあるため、**Ux/Uy/P は再現でき、ρ だけ再現できないという報告と整合します**。時点混在の注意は既に [limiter_d.cu:258](/home/sano/work/forge-species/solver_density_cuda/cuda_forge/limiter_d.cu:258) にもあります。

  反証条件: ρ の時点を一つ前へ揃えても、同じ節点・step で O(1) の ψ 不一致が残れば、**単純な 1 step のずれだけが原因という仮説**は棄却します。BC/EOS による途中の密度書換えまで除外したことにはなりません。

第 2 仮説: **H-R：入力場・勾配の float32 変動が、狭い遷移幅のリミッタで増幅されている。** 確度: **中〜低、run 未確認**。報告された数 ulp の極値差とは整合しますが、床への因果寄与は未証明です。

第 3 仮説: **有限精度とは別に、入口 BC・再構成・外反復の非線形結合が振動を維持している。** 確度: **低、未確認**。CFL 半減への弱い応答だけでは除外できません。

判別 A/B:

- **変更点は後処理の ρ の参照 step のみ**。対照は現在の `ρₙ, ∇ρₙ, ψₙ`、比較側は `ρₙ₋₁, ∇ρₙ, ψₙ`。owner と全近傍の ρ をともにずらし、幾何・参照値・K・演算方法を固定します。Ux/Uy/P はずらしません。
- **forge の追加計算は 0 step**。開始・終了×CFL 2・1 の既存 4 窓について、連続保存された全隣接ペアを使います。`res_0` を使わなくても各 100 枚から 99 ペアを検査できます。
- 入口集合 S を固定し、保存された `limiter_ro` に対する誤差の最大値・p99・`2e−7` 超過率を、入口境界列／2 列目に分けて比較します。これは再現性の検査であり、収束判定ではありません。
- **A：時点を揃えた側だけ全窓で既報の再現精度 2e−7 以内に入る** → ρ の不一致は後処理で説明でき、これを実装不整合の証拠から外せます。
- **B：同じ箇所に O(1) の不一致が残る** → 単純な時点ずれ説を棄却します。中間的な改善は「一部説明」に留め、合格にしません。

やらない方がよいこと: **現段階で 48000 step の K／FP64 比較を先行させること、ρ のみ除外して因果を確定すること、ψ 凍結を導入すること。** 関数だけ double にして変化がなくても H-R は否定できません。全域 FP64 で改善しても、流束・勾配・EOS・更新のどこが効いたかは特定できません。

呼び出し側の前提への異議:

- **「H-R と整合」は条件付きで採用**します。ただし観測は「δ± が小さい／ゼロ、ψ が変動する」であり、「丸めによって極値判定が切り替わる」は解釈です。δ±=0 は正確な局所極値でも生じます。
- **Major：制限面の入替率を原因の証拠にしないでください。** `face_probe.py:42` はクランプ前の `argmin` を記録します。最小値が同値・近接している場合や、最終 ψ が 1 にクランプされる場合にも面番号は変わります。対案は、候補面の ψ の差と、実際の再構成補正 **ψδm** の変動を併記することです。
- **Major：因果の向きは現在の指標では決まりません。** ψ と残差は同じ場から作られ、残差による更新が次の ψ を変える閉ループです（`main.cpp:1866,2083,2130`）。必要な証拠は、状態・BC・勾配の変化から分離した ψ の変更に対する残差応答です。現時点では「相互に結合している可能性」に留めてください。
- **Minor：末尾窓の検査が記載より狭いです。** `body_compare.py:9` は CSV に対して末尾 10%／20% を計算しますが、`R_入口` は同:13、20 の **末尾20%だけ**です。対案は `R_入口` も両区間で集計することです。

不足情報: 対象の `case/16.nozzle_wys/run_0524_floor_dry_L1/`、`run_0527_liminlet_ctrl_cfl2/`、`run_0528_liminlet_cfl1/` と毎 step 窓はローカルにありません。実効 config、HDF5、残差履歴、各 `CONVERGENCE_VERDICT.txt`・準定常 VERDICT、終了窓の run パスが必要です。**0.81 倍・0.83 倍・ALL STEADY はブリーフの報告値であり、今回は独立検証していません。** 関連する主要コードは `8ce8da3a` と現 HEAD の差分がないことを確認しました。ファイル変更・forge 起動は行っておらず、提案は **plan 未反映**です。
