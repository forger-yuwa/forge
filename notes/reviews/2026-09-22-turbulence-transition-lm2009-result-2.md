# codex レビュー: turbulence-transition-lm2009 (result)

- **plan**: [`plans/active/turbulence-transition-lm2009.md`](../../plans/active/turbulence-transition-lm2009.md)
- **stage**: `result` (diff base `ded72e17`)
- **date**: 2026-09-22
- **commit**: `23a4ea02` (feature/sern-design)
- **codex**: effort `high`, 7.2 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M6/m2
- **focus**: result 2 巡目。1 巡目 (notes/reviews/2026-09-22-turbulence-transition-lm2009-result.md, NO-GO M5/m4) の指摘が実際に解消されたかの確認と、その後に足した結果の審査。見てほしい点: (1) M1 初期出力/restart、M2 --from-floor、M3 単体検査の強化が diff と plan §5.1/§6.1 の記述どおりか (main.cpp の初期経路、tools/check_convergence.py analyze_from_floor、tools/check_lm_kernel.py)。(2) plan §6.2 の書き直し (未収束・準定常・基準内、丸め床、前縁説の撤回、y1+ 全壁) が過不足ないか。(3) plan §6.3/§6.4 の翼の結論『負圧面層流域の過大は消えるが全体の一致は入口粘性比次第、遷移が遅く急なのはモデル固有 (SU2 と -0.4〜+1.4 % で一致、遷移位置 0.337 で同一)』は、提示した数値から言ってよい範囲か。言い過ぎがあれば具体的に。(4) status を done にして accepted へ移してよいか、それとも CHT 未実施・T3B 実験値未入手・SU2 T3A の量の時系列 (#4b) 未取得を理由に active のままにすべきか、どちらか一方を推奨してほしい。case/53.c3x_vane_cht/README.md と case/54.markii_vane_cht/README.md は巨大なので『遷移モデル』で grep した周辺だけ読むこと。*.log, residual_history.csv, res_*.h5, plans/README.md は読まないこと
- **extra**: `case/57.transition_flat_plate/README.md`, `case/53.c3x_vane_cht/su2_smooth_lm/COMPARE_FORGE.txt`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
前回 M1・M2 の修正は確認できましたが、M3 の判定器には偽合格が残っています。  
T3B の検証記録は報告値と整合せず、翼の「モデル固有と確定」も証拠を超えています。`active/in_progress` の維持を推奨します。

指定の diff を取得し、実装・保存済み VERDICT・比較コードを確認しました。前回 M1 は [main.cpp:1236](/home/sano/work/forge/solver_density_cuda/main.cpp:1236) の初期化経路と [interp_field.py:133](/home/sano/work/forge/solver_density_cuda/tools/interp_field.py:133) の保存量優先で修正されています。M2 はメモリ上の回帰確認で、遷移残差列欠落・遷移有無の異なる参照を拒否し、正常入力を通しました。`Pk_base`、両経路の `F₁` 補正、`2(μ+μt)`、周期の部分体積・gather・mirror に新たな誤りは見つかりませんでした。

1. **Major — 単体検査が NaN・必須診断の欠落を `PASS` にする。前回 M3 は未解消。**

   **根拠:** [check_lm_kernel.py:57](/home/sano/work/forge/solver_density_cuda/tools/check_lm_kernel.py:57) は `ratio > 1` の個数だけで判定するため、NaN が違反に数えられません。対角と反復数も50・64行で任意扱いです。

   実ファイルを使わない合成入力で、以下を再現しました。

   - `lmPgamma` の1節点を NaN → `最大 |差|/許容 = nan`、超過節点0、**`VERDICT: PASS`、終了コード0**。
   - `src_jac_gamma`・`src_jac_reth`・`lmCorrIter` を削除 → **`VERDICT: PASS`、終了コード0**。

   また、44行の評価対象判定が検査される側の `lmFlength > 0` に依存し、誤ってゼロになった節点を除外できます。

   **対案:** 必須場・形状・全節点の有限性を先に検査し、欠落・非有限・評価対象なしを FAIL にする。評価対象は壁・入口・速度など独立した入力から定義する。この反例を回帰試験に追加してから `run_0014`・`run_0154` を再判定してください。

2. **Major — T3B の `STEADY` は、報告した遷移位置を判定していない。`Cf` の座標ラベルにも欠陥がある。**

   **根拠:** [run_0013 の cf_series.csv:2](/home/sano/work/forge/case/57.transition_flat_plate/run_0013_t3b_lm/cf_series.csv:2) は全10時刻で `x_onset=x_end=1.5 m`。再実行した `check_quasisteady.py` も、その定数列に対して **`STEADY`** を返しました。plan の開始 **0.090 m**・終了 **0.189 m** の定常性は未確認です。

   さらに [cf_plate.py:80](/home/sano/work/forge/case/57.transition_flat_plate/tools/cf_plate.py:80) の `%.1f` は、測定位置 **0.05 と0.1を両方 `cf_x0.1`** にします。保存列を名前どおり読むと、最終 `Cf(0.1)=0.00520010` に対し、SU2 の `restart_flow.csv` から同じ抽出関数で再計算した値は **0.00478174**、差は **+8.75%**。記載された **+1.7%** と整合しません。座標誤記の可能性があるため、これだけで流れ場の不一致とは断定できません。

   **対案:** 座標を失わない一意な列名に直し、現行の極値抽出で時系列を再生成する。5測定位置・遷移開始終了を両コードで同一定義に揃え、数値表と VERDICT を更新するまで §5.1 #5 を完了扱いにしないでください。

3. **Major — C3X の比較は「モデル固有と確定」を支えない。**

   **根拠:** [COMPARE_FORGE.txt:1](/home/sano/work/forge/case/53.c3x_vane_cht/su2_smooth_lm/COMPARE_FORGE.txt:1) の領域平均差 **−0.4〜+1.4%** は確認できました。しかし比較コードは最終場の領域平均を取り、遷移位置は **41点移動平均後、閾値を初めて超える節点**を返します。[compare_su2_lm.py:30](/home/sano/work/forge/case/53.c3x_vane_cht/tools/compare_su2_lm.py:30)  
   **0.337 が同じ**とは、この抽出・格子解像度で同じ位置という意味です。

   forge の [準定常判定:2](/home/sano/work/forge/case/53.c3x_vane_cht/run_0153_lm_cf0_su2ctrl/QUASISTEADY_VERDICT.txt:2) は **`OSCILLATING`**。`SS_post` は **1.125 ± 0.073ポイント**で、遷移位置・遷移幅は判定対象にありません。SU2 の対応する局所量の VERDICT もありません。加えて、この比較は壁解像 FAIL の2 µm格子・一様壁温・入口粘性比10に限定されています。

   **対案:** 結論を「この条件の最終場では、両実装に同様の傾向と小さい領域平均差が見られる」に限定する。両者の局所 `h`、遷移位置・幅、報告する実験偏差 rms の時系列を判定し、共通の格子・入口条件の影響を切り分けるまで「モデル固有と確定」を撤回してください。Mark II の「そのまま合う」も、対照 `run_0036` が壁解像 FAIL である範囲の機構診断に留めるべきです。

4. **Major — 「残差床は単精度の丸め」とする新しい根拠は、次元が合っていない。**

   **根拠:** [plan:181](/home/sano/work/forge/plans/active/turbulence-transition-lm2009.md:181) は `rms_ro/⟨ρ⟩` を機械イプシロンと比較しています。しかし [SLAU:564](/home/sano/work/forge/solver_density_cuda/cuda_forge/convection/convectiveFlux_slau_d.inc.cuh:564) は面積込み質量流束を残差に加算し、[residualMonitor_d.cu:151](/home/sano/work/forge/solver_density_cuda/cuda_forge/residualMonitor_d.cu:151) はその RMS を出しています。`rms_ro` は密度の更新量ではなく、密度で割っても無次元になりません。

   [RESIDUAL_MAP.txt:3](/home/sano/work/forge/case/57.transition_flat_plate/run_0015_t3a_lm_resmap/RESIDUAL_MAP.txt:3) は前縁寄与が小さいことを示しますが、原因が丸めであることや「全域に一様」までは示しません。

   **対案:** 前縁説の撤回は維持し、丸めは仮説に戻す。全保存量について更新量を状態量・ULPと比較するか、同一条件の倍精度対照で床と報告量への影響を測る。`NOT CONVERGED` の扱いは維持してください。

5. **Major — T3A の数値誤差ゲートを「基準内」とする総括は、壁解像条件と矛盾する。**

   **根拠:** [plan:149](/home/sano/work/forge/plans/active/turbulence-transition-lm2009.md:149) の事前条件は全格子で局所 `y₁⁺≤1`。一方、[結果:170](/home/sano/work/forge/plans/active/turbulence-transition-lm2009.md:170) は最大 **1.67／1.31／1.02**、超過壁長 **0.34／0.12／0.08%** と記載しながら、ゲート全体を「基準内」としています。

   3格子の幾何品質は保存記録で **`VERDICT: PASS`** ですが、壁解像条件とは別です。

   **対案:** 「選定した報告量の格子差は基準内、全壁の解像条件は未達」と分ける。前縁を除外・許容するなら、その範囲と局所摩擦への影響を検証して例外を記録する。最大値の表見出しにも `x>0.02 m` を明記してください。

6. **Major — 当初スコープの未完了項目が残り、一部は完了行に埋め込まれている。**

   **根拠:** [plan §2:24](/home/sano/work/forge/plans/active/turbulence-transition-lm2009.md:24) は T3B の実験照合も対象です。しかし [残作業表:124](/home/sano/work/forge/plans/active/turbulence-transition-lm2009.md:124) では、SU2 T3A の時系列と CHT は未完、T3B は「実験値未入手」のまま取り消し線付きです。§6 の翼の入口粘性比 **1／10／100** に対し、Mark II は10／100の記録です。

   **対案:** T3B 実験照合、両コードの報告量時系列、Mark II の未実施感度、CHT を独立した未完行として残す。入口長さスケールの不確かさは CHT の結果解釈を制約しますが、未実施を完了に変える理由にはなりません。

7. **Minor — SU2 の「積分熱流束は後半5桁不変」は実測と違う。**

   **根拠:** [plan:211](/home/sano/work/forge/plans/active/turbulence-transition-lm2009.md:211) に対し、`case/53.c3x_vane_cht/su2_smooth_lm/history.csv` の後半20000反復の `HF` は **37072.48268〜37105.69628**、相対幅 **0.0896%**。5桁不変ではありません。`check_quasisteady.classify_series` の既定条件では **`STEADY`** でした。

   **対案:** 不変桁数を削除し、評価区間・変動幅・VERDICT を記載する。積分量の判定を局所量へ拡張しないでください。

8. **Minor — 現在仕様と plan に旧記述が残る。**

   **根拠:** [plan:134](/home/sano/work/forge/plans/active/turbulence-transition-lm2009.md:134) の単体許容「相対 `1e−5`」は現行の `1e−3·|参照|+1e−5·最大値` と異なります。§4.2 の陰的対角も変更ログの「項ごとの負部」と未同期。§5.1 #5b・§6.3 は SU2 を「実行中」とし、§7 は旧ファイル名 `transitionSource_d.cu` のままです。

   **対案:** 現行仕様を本文へ反映し、旧方針は履歴として明示する。`methods/index.md` の既存リンクは有効で、遷移レシピ・SU2との上下限差の追記は確認しました。`plans/README.md` は指定に従い未確認です。

**推奨は、`active/in_progress` のまま維持することです。** 優先順は、判定器の偽合格修正 → T3B の抽出・時系列再検証 → 翼の結論と丸め床説明の訂正 → 壁解像・未完検証の完了です。負圧面前半の平均過大が縮小する傾向、全体誤差が入口粘性比に依存する傾向は残せますが、現段階で計画全体を `done` にする根拠はありません。

ファイル変更・禁止対象の読み込みは行っていません。以上はレビュー提案であり、**plan 未反映**です。

指摘数: Critical 0 / Major 6 / Minor 2
