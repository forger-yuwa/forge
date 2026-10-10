# codex レビュー: time_integration-line-viscous-jacobian-faceh (plan)

- **plan**: [`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `dec9e6af` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 3.1 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2
- **focus**: §6.3 (2026-10-10 に追加した U-J の列ごと・向き・壁拘束の照合の事前登録) と、それを実装した solver_density_cuda/tools/test_line_visc_jacobian.cpp の差分 (作業ツリーの未 commit の変更。git diff で見る) だけを点検する。§6・§6.2 は実施済みで対象外。試験はまだ実行していない。点検してほしいこと: (1) 列ごとの判定 (s_c で割る、h と h/2 の再現、零列の無次元の絶対誤差、解像しない列 1 % まで) が諮問の判別 A/B を満たすか、抜けや甘さはないか。(2) 面の向きの整合 (1″) と速度だけ・温度だけ固定の隣 (3′) の組み方が正しいか。(3) 実装が §6.3 の文言どおりか (乱数の 200 組と端の場合、FAIL と判別不能の区別、旧 (1) を記録だけにした扱い)。(4) この試験が PASS/FAIL のときに言えること・言えないことの書き方。ソルバは変えない。
- **extra**: `solver_density_cuda/tools/test_line_visc_jacobian.cpp`, `solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh`, `notes/reviews/2026-10-10-lvc-faceh-result-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
列別照合へ進む目的と、面交換・壁拘束の微分式は妥当です。  
ただし現状の判定器は、未確認や閾値超過を PASS に含めます。以下を正式実行前に修正してください。

1. **Major — 壁拘束の照合を全部スキップしても PASS にできる。**

   [test_line_visc_jacobian.cpp:208](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:208) は解像しない方向を除外しますが、[最終判定:261](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:261) は `n_fix_unres` を参照しません。400 方向すべてが解像しなくても、誤差の初期値 0 が残り、他項目が通れば PASS です。

   また、[壁試験:190](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:190) は各状態で乱数方向を一本照合するだけです。[§5.1 #4:61](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:61) の「自由度を消去した系」の全列照合にはなっておらず、小さい自由度の誤りを大きい作用で隠す余地が残ります。

   **対案:** 速度固定は自由変数 `(ρ, ρE)` の2列、温度固定は `(ρ, ρu, ρv, ρw)` の4列について、拘束の写像 B を使った `K B` と独立差分を列別照合する。必須の拘束列が未解像なら `INDETERMINATE` とし、乱数方向は補助確認にする。

2. **Major — NaN/Inf が不合格になる保証がない。**

   [列集計:70](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:70) は `std::max` と大小比較だけで、有限性を検査していません。`std::max(有限値, NaN)` は有限値を残すため、非有限の差分が誤差 0 として消えます。[最終判定:262](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:262) の `誤差 > 閾値` も NaN を捕捉しません。

   判定式を Python で再現した人工入力では、`J = I`、両差分行列がすべて NaN の場合に、**解像5列・不一致0列・最大誤差0**となりました。これは実標本に NaN が出たという報告ではなく、判定器の反例です。

   **対案:** 状態・流束・D/K・差分・正規化量・線形解を集計前に `isfinite` で検査し、一件でも非有限なら PASS を禁止する。差分評価不能と製品 Jacobian の非有限を理由付きで区別する。

3. **Major — 閾値をまたぐ列が、FAIL にも判別不能にもならない。**

   [check_cols:83](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:83) は再現性を通過後、**両方の誤差が** `1e-6` を超えた場合だけ失敗にします。

   人工入力で `s_c = 1`、誤差が h 側 `1.04e-6`、h/2 側 `0.96e-6` の場合、差 `8e-8` は再現性条件を満たし、**最大誤差が許容外なのに不一致0列**になります。[§6.3:133](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:133) の合格条件との間に空白があります。

   **対案:** 列の分類を明示する。両幅で許容内なら合格、再現性を満たして両幅で許容外なら FAIL、片側だけ許容外なら `INDETERMINATE`。壁方向にも同じ分類を適用する。

4. **Major — 零列の処理が登録と異なり、再現性も検査していない。**

   [零列分岐:74](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:74) は F1 だけを検査し、F2 を無視します。人工入力では、零列について **F1 = 0、F2 = 1 でも合格**しました。また、登録の `|q_c|` を実装では `max(|q_c|, 1e-3)` に変更しています。

   [登録式:134](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:134) の分母は異なる残差行の最大値です。運動量行の分子をエネルギー行の尺度で割る場合、物理的に無次元ではありません。さらに通常の標本は β・κ が正で、運動量対角と熱伝導のエネルギー列が非零になるため、零列分岐を確実に通す標本がありません。

   **対案:** 保存量・残差の基準尺度を事前固定して無次元化し、零列でも両差分幅の誤差と再現性を検査する。β = 0 の静止状態など、既知の零列を持つ標本と、判定器への人工入力試験を追加する。

5. **Major — 全体の1％許容では、重要な未解像列を覆い隠せる。PASS/FAIL の解釈を狭める必要がある。**

   [標本生成:146](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:146) は200組＋端の30組なので、[判定:261](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:261) は **2300列中23列**を未解像のまま許容します。例えば特定の端条件6組の D/K のエネルギー列12本が全滅しても、総数条件だけなら通ります。出力も集計値だけで、未解像の標本・側・列を特定できません。

   したがって、[§6.3:147](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:147) の「第2仮説を退ける」は、未解像部分まで含めて読める点で強すぎます。また実装の FAIL には短いラインや float 比較の失敗も含まれ、すべてを double の微分不整合と解釈できません。

   **対案:** 端条件・零列・拘束列は全件解像を要求し、乱数標本の許容分は標本ID・D/K・列・両誤差・再現性を保存する。PASS は「解像した対象で不一致を検出しなかった。未解像部分は保留」と記す。FAIL は失敗項目を分け、独立差分との再現する不一致だけを微分実装への反証とする。

6. **Minor — 標本分布の説明が実装と違う。**

   [生成:95](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:95) は速度の**各成分**を ±1700 に制限するので、速度ノルムの上限は約2944 m/sです。また、立方体内の一様乱数を正規化する [法線生成:102](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:102) は球面一様分布ではありません。[登録:131](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:131) と一致しません。

   **対案:** 既存分布を維持する方針に合わせ、文言を「各速度成分」「立方体一様乱数の正規化」に訂正する。200組＋端30組の追加自体は実装されています。

7. **Minor — float/double の列別記録が未実装。**

   [§6.3:144](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:144) は列別値を記録するとしていますが、[比較:139](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:139) と出力は全行列正規化の最大値だけです。

   **対案:** float/double の列別値を記録に追加する。旧 (1) の `e_fd` を記録専用にした変更は、登録どおりで問題ありません。

数学的な組み方については、**(1″) の D(j,i) = K(i,j)、K(j,i) = D(i,j) は正しい**です。現在のモデルは P = I + nnᵀ/3 なので、法線を反転しなくても同じ P になります。ただし、確認できるのは共通関数の交換対称性までで、製品側の面重み・`Kprev/Knext` 配置ではありません。

**(3′) の速度固定方向と温度固定方向も正しい**です。温度固定の式を ∂T/∂Q に代入するとゼロになります。`Q ± hδQ` は有限幅で厳密な等温状態ではありませんが、接方向の中心差分としては妥当です。問題は上記の検査範囲と判定処理です。

**推奨は、上記1〜5を優先順に事前登録へ反映し、単体試験だけを修正して実行することです。** 修正後の PASS でも、実残差・CUDA経路・float32の列別精度・周期／軸対称・case/45 の安定性は保証しません。FAIL でも、その不一致が case/45 の破綻原因であるとは言えません。

確認したのは未 commit 差分、関連仕様・親 plan、構文検査、判定式の人工入力検算です。C++単体試験本体と `forge` は未実行、§6・§6.2 の再評価はしていません。ファイル変更なし、**plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 2
