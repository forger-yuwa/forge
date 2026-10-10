# codex レビュー: time_integration-line-viscous-jacobian-faceh (plan)

- **plan**: [`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `fe878c4b` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 3.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m1
- **focus**: §6.9 (2026-10-10 に追加した、LHS の座標の差だけを変える A/B の事前登録) と、台本 lvcgeom.sh・判定 lvcgeom_judge.py (作業ツリーの未 commit の新規ファイル) だけを点検する。諮問 (lvcaudit-result-diagnose) は同じバイナリに診断の切替を足す形を勧めたが、元のセッションの float 化の段 ② (LHS は ST(ccx[o]) − ST(ccx[ic])、952 行) と段 ③ (ST(±ge_x)、965 行) の FP64 のバイナリの組で行うことにした (段 ③ が本番に入る変更そのもの、元のセッションの記録で FP64 の残差の面の流束が段 ② と段 ③ でビット一致)。この置き換えが判別 A/B として妥当か、抜けている確認は何か。介入の成立の確認 (1 step の書き出しで状態のビット一致と壁際の K の差 > 1e-3) が十分か。ほかは §6.2 の判定器と同じ分類・ゲート。元のセッションのソースは /home/sano/work/forge-integ-1005 (ブランチ feature/nozzle-wall-fit-and-pipeline) にあり、float 化の plan は plans/active/architecture-float-state-double-geometry.md の §6.5・§6.7 (段 ②・段 ③ の結果)。
- **extra**: `case/45.isobutane_m6_d155/lvcgeom.sh`, `case/45.isobutane_m6_d155/lvcgeom_judge.py`, `notes/reviews/2026-10-10-lvcaudit-result-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

段②・段③のバイナリを使う方針は条件付きで妥当です。ただし、残差の不変性を裏付ける記録の解釈と、介入ゲートに不足があります。以下を実装・登録してから本試験へ進めてください。

1. **Major — `run_0419` は、段②・段③の別バイナリ間のビット一致を証明していない。**

   [plan:302](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:302) の記述は、根拠の測定より強い主張です。[fg4.sh:26](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/fg4.sh:26) は段③のバイナリだけを起動し、その内部で旧式・新式の粘性および k/ω 拡散を評価しています。[geomAbDiag.cpp:235](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/geomAbDiag.cpp:235) もこの動作を明記しています。入力も今回の `run_0183/res_100000` ではなく、B0 の保存場です。

   ソース差分では block DPLUR の変更は意図した座標差の置換ですが、段③には本条件でも動く粘性・スカラー拡散の変更が含まれます。既存記録だけで、別バイナリのコード生成差まで除外できません。

   **対案:** 記述を「段③内の旧式・新式の比較」に訂正する。両バイナリのソース・ビルド条件を固定し、今回と同じ凍結入力で、本番経路の面流束と最初の更新前の残差・RHSを比較する。原子加算による差は旧腕の再実行で測り、許容を事前登録する。この確認を独立したゲートにする。

2. **Major — 「K が変わった」は「幾何の誤差が修復された」の代わりにならない。**

   [判定器:323](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcgeom_judge.py:323) は、合成された `Kprev/Knext` のどこか一列で相対差が `1e-3` を超えれば合格にします。誤差が減った方向も、変化が幾何に由来することも検査しません。人工入力で **`K_new = −K_old` でも全介入ゲートを通る**ことを確認しました。

   また、幾何の変更は K に加えて薄層の D、値3のスカラー対角、ライン外のキー5の熱伝導対角にも効きます。K だけでは介入全体を確認できません。

   **対案:** [元の諮問:45](/home/sano/work/forge-faceh/notes/reviews/2026-10-10-lvcaudit-result-diagnose.md:45) の「対象係数の double 参照に対する差 ≤ `1e-5`」を成立条件へ戻す。生の座標・面積・物性から独立に作った参照に対し、対象面の `δ/dcc`、β・κ、関連する対角係数を照合する。K の差は補助記録とし、D/K の変化がこの係数変更から説明できることを確認する。

3. **Major — 書き出し2本の入力・形・有限性が十分に検査されていない。**

   [入力ゲート:270](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcgeom_judge.py:270) の対象は本試験4本だけです。書き出し2本は [intervention:311](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcgeom_judge.py:311) でバイナリのハッシュを確認するだけで、出発場・設定・格子を検査していません。要求した5節点・5本のライン、配列間の行数、状態の有限性も必須条件ではありません。

   人工入力では、**要求外の1節点だけの記録も、両腕の状態が同じ NaN の記録も全ゲートを通りました**。整数表現での一致は有限性を保証しません。

   **対案:** 書き出しにも共通の入力ゲートを適用し、本試験との差を step 数・出力間隔・診断設定に限定する。要求5節点と5ライン、重複なし、配列形 `N×2`・`N×7`・`N×25`、有限性・必要な正値を検査する。壁側の起点は既存の壁フラグと接続で確認する。同一入力の確認には、保存量に加えて乱流・組成・物性、`dt_vol`、拘束フラグも含める。

4. **Minor — 介入ゲートが本試験の前に実行されない。**

   [台本:72](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcgeom.sh:72) は書き出し2本の直後に本試験4本を起動します。`dump1` が確認するのは終了コードと書き出しログだけで、状態一致や K の条件は未評価です。「介入を先に確認する」という登録になっていません。

   **対案:** 判定器に事前確認専用モードを設け、書き出し2本の後に呼び、合格した場合だけ本試験へ進む。証拠不備の `INVALID` と、正常な証拠で修復成立条件を満たさない「判別不能」も分ける。

**推奨:** 段②・段③のバイナリ組を維持し、上記1〜3を事前ゲートとして実装、4で実行順序を保証する方針に絞ります。各腕2本・2000 step、Aの破綻再現を要求する分岐、有限を収束と呼ばない制限は妥当です。今回の目的に cell・周期・速度評価を追加する必要はありません。

シェル構文確認、Python構文確認、メモリ上の人工入力試験を実施しました。AWS の実バイナリ・run データは再測定していません。ファイル変更・forge 起動なし。**plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 1
