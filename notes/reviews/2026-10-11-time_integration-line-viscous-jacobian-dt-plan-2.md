# codex レビュー: time_integration-line-viscous-jacobian-dt (plan)

- **plan**: [`plans/active/time_integration-line-viscous-jacobian-dt.md`](../../plans/active/time_integration-line-viscous-jacobian-dt.md)
- **stage**: `plan`
- **date**: 2026-10-11
- **commit**: `0e2228b1` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 4.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m3
- **focus**: §6.3 (2026-10-11 に追加した、上限 50 での全部入り L7・生産 L0・粘性入り L5 の 3 本を、速度 plan §6.15 の水準に届く step 数で比べる事前登録) と、台本 lvcmp.sh・事前のゲート lvcmp_pregate.py・本判定 lvcmp_judge.py (作業ツリーの未 commit の新規ファイル) だけを点検する。§6 (dt の A/B) の lvcdt_*.py と faceh の lvc75_*.py を元にし、見張りは速度 plan の m9_watch.py (--phase line) をそのまま使う。V(b) の運動量の薄層の項の独立計算、比較の閾値 15000 step・各腕 1 本という設計、分岐と言える範囲が妥当か、抜けている確認は何か。設計の諮問は diagnostician (notes/reviews/2026-10-11-lvcdt-result-diagnose.md の Q2)。
- **extra**: `case/45.isobutane_m6_d155/lvcmp.sh`, `case/45.isobutane_m6_d155/lvcmp_pregate.py`, `case/45.isobutane_m6_d155/lvcmp_judge.py`, `case/45.isobutane_m6_d155/m9_watch.py`, `notes/reviews/2026-10-11-lvcdt-result-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

§6.3 の同一バイナリ・上限 50 による比較と、V(b) の運動量 Jacobian の式は妥当です。
ただし、本判定器が不完全な証拠を採用する穴があり、各腕 1 本で主張できる範囲も修正が必要です。
以下はコード照合とメモリ内の反例試験による指摘です。実 run の HDF5 は手元になく、監査値の再計算はしていません。

1. **Major — 不完全な系列・終端記録から、有効な比較結果を作れてしまう。**

   根拠: [lvcmp_judge.py:159](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcmp_judge.py:159)。`assess()` は `run`・`phase`・`budget`・`consec`、系列の欠落・重複、到達 HDF5 の存在と `reach_sha256` を照合しません。`CENSORED` では予算までの完走、`DIVERGED` では発散の証拠を確認しません。

   実際に、空の `rows` が `CENSORED`／`DIVERGED` として採用されました。また、step 2500・22500 の **2 点だけ**でも、欠損と θ_r が一定なら `reach_level()` は 22500 到達を返します。前者は L0 が到達していれば誤った NO-GO に直結します。

   **対案:** 本判定前に記録のスキーマと実行条件を照合し、2500 step ごとの系列を終端まで要求する。`CENSORED` は `RUN_RC=0`・150000 到達・完全な系列、`DIVERGED` は残差または非有限出力の実証を必須にする。`REACHED` は保存された到達 HDF5 のハッシュと量も再確認する。

2. **Major — 「判定不能」の残差を、許容される `NOT CONVERGED` と取り違える。**

   根拠: [lvcmp_judge.py:104](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcmp_judge.py:104) と同ファイル178行。本文を捨てて VERDICT だけを取り出し、`DIVERGED` でなければ通します。一方、[check_convergence.py:126](/home/sano/work/forge-faceh/solver_density_cuda/tools/check_convergence.py:126) は残差列欠落を「判定不能」と記録しますが、最終見出しは `NOT CONVERGED` になります。この組合せを反例で確認しました。

   **対案:** 「入力が有効で未収束」と「入力不足で判定不能」を分離する。必須の平均流・SST 残差列、有限性、判定区間の終端までの記録を検査し、本文に判定不能があれば INVALID にする。準定常判定も、非定常という有効な結果と、実行失敗・VERDICT 欠落を区別する。

3. **Major — 15000 step は実用上の採否幅には使えるが、各腕 1 本で再実行の揺れを除いた差とは判断できない。**

   根拠: [plan §6.3:181](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-dt.md:181)、同213行、および[諮問 Q2:38](/home/sano/work/forge-faceh/notes/reviews/2026-10-11-lvcdt-result-diagnose.md:38)。出力6回分という説明は観測間隔の説明であり、**到達 step の再実行変動**の上限ではありません。1 step の再実行による行列・残差の照合も、その代わりにはなりません。

   **対案:** 各腕1本・15000 step は維持し、今回を探索的な比較と明記する。分岐の文言は「今回の実行では登録幅以上早かった／差を検出しなかった／15万 step 内に届かなかった」に限定する。「手法の利点を支持」「熱伝導 K が到達を変える」という再現性を含む結論は、別途確認後にする。追加投資を止める判断としての NO-GO は可能です。

4. **Minor — V(b) は運動量 K の成立確認であり、L5−L0 全体の介入確認ではない。**

   根拠: [lvcmp_pregate.py:323](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcmp_pregate.py:323) の補間、βP∂u_j/∂Q_j、拘束処理は[共通関数:120](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh:120) と整合しています。ただし L0 の D とエネルギー K は照合対象外です。さらに[timeIntegration_d.cu:1263](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1263) により、キー5でも L5/L7 では等温壁の D[4][0] が −ρE/ρ、L0 では0になります。

   **対案:** この拘束差を §6.3 に明記し、少なくとも壁の D 行を期待値と照合する。V(b) の合格範囲は「運動量 K の係数変更が成立」と限定する。L7−L0 は構成全体の比較、熱伝導 K の分離は L7−L5 として扱う。

5. **Minor — 登録したスナップショット保持方針と実装が異なる。**

   根拠: [plan:191](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-dt.md:191) は「5万 step ごと・最新5つ・到達」を保持するとしていますが、[m9_watch.py:123](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/m9_watch.py:123) は最新1つと指定された到達出力だけを残します。

   **対案:** 既存 watcher をそのまま使う方針なら、plan を実装に合わせる。その場合、削除前に各出力のハッシュ、全保存量の有限性、ρ・P・T の正値確認を系列に保存する。現在の非有限数だけでは正値確認を代替できません。

6. **Minor — 見張りの異常終了が台本の失敗として伝播しない。**

   根拠: [lvcmp.sh:100](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcmp.sh:100)。`wrc` は記録するだけで、`arm()` の最後は `echo` なので、見張りが失敗しても次の腕へ進みます。80行でも既存の `EXEC_ERROR`・`DATA_ERROR` を正常な完了と同様に扱います。

   **対案:** 正常に分類された `REACHED`・`CENSORED`・`DIVERGED` だけを継続可能とし、見張りの失敗・実行エラー・データエラーでは台本を非ゼロ終了させる。

推奨は、**1→2→3の順で本判定と主張範囲を直し、4〜6を同期してから、既定の3腕を探索的比較として実行すること**です。終了水準の変更や point 仕上げの追加は不要です。ファイルは変更しておらず、指摘は plan 未反映です。

指摘数: Critical 0 / Major 3 / Minor 3
