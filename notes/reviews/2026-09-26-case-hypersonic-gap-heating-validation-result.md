# codex レビュー: case-hypersonic-gap-heating-validation (result)

- **plan**: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-26
- **commit**: `813403db` (feature/gap-heating-precision)
- **codex**: effort `high`, 8.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M9/m1

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

段階起動・比較点・観測モデルに不整合があり、D4 の原因判定は再現できません。  
「実験がモデルの適用範囲外」「CFD が下界を与える」という結論を維持して `accepted` へ移すことはできません。

1. **Major — 3D の段階起動で、前段の解を引き継いでいない。**  
   [make_case3d.py:223](/home/sano/work/forge/case/56.gap_tp1187/tools/make_case3d.py:223) は各段の終了後に `seed_turb()` を呼ぶだけで、出力を次段の `mesh.h5` に移していません。設定は全段 `valueFileName: mesh.h5` です。実際、`case/56.gap_tp1187/run_0063_ab3d_L30/` の本段 `res_0.h5` は、全 1,622,178 点で `ro/roUx/roUy/roUz` が初期 `mesh.h5` と一致します。`roUy/roUz` は全点ゼロですが、前段 `res_2500.h5` では非ゼロです。  
   **対案:** 各段の指定最終出力を `restart_field.py` で index コピーし、次段開始時の保存量一致を確認する。現行結果を「梯子を経た結果」とする記述・累積反復数も訂正する。

2. **Major — 熱電対の局所測定を、別位置・帯平均と比較している。**  
   [skin_forward3d.py:32](/home/sano/work/forge/case/56.gap_tp1187/tools/skin_forward3d.py:32) は TC93/94 を上面の 0.25/0.51 cm に配置しています。しかし原報 Fig.5 と登録内容では TC93/94 は半径部、TC95 が上面側です。TC95 自体も出力から欠けています。さらに同ファイルの 76 行では `z≤6 mm` の節点平均を採用していますが、原報本文 pp.4–5 は測定中心線を縦すきま中心線へ合わせ、位置のばらつきを ±0.03 cm としています。  
   `run_0052_ab3d_open_long` の TC92 相当の生熱流束は、帯平均 −32.46 kW/m² に対し中心線 −48.96 kW/m²です。比較演算だけで約51%違います。  
   **対案:** 熱電対ごとの座標台帳を確定し、その位置へ補間する。帯平均は診断として分離し、TC95 を含めて D1/D2 を再評価する。

3. **Major — 0.25秒の「平均化窓」を「加熱経過時間」と取り違えている。**  
   [skin_forward3d.py:115](/home/sano/work/forge/case/56.gap_tp1187/tools/skin_forward3d.py:115) は `t=0.25 s` まで加熱して瞬時の温度変化率を返します。一方、原報本文 p.6 は、模型挿入が約2秒強、解析時刻はトンネル中心到達時、温度を0.25秒区間で平均して区間間の上昇率を求めた、と記しています。現在の演算はこの測定手順を再現していません。  
   同じ保存場でも、モデルの経過時間を0.25→2秒に変えた感度試算では、TC90 の観測後 L30/L15 比が **0.0197→約0.537** に変わります。2秒が正解という意味ではなく、時間の不確かさが判定を左右する証拠です。  
   **対案:** 加熱履歴・解析時刻・平均化窓を別々に定義し、原報と同じ平均・差分演算を実装する。未確定な履歴は感度幅として残す。

4. **Major — D4 だけ観測モデルを省略し、「読み C」を確定している。**  
   [acceptance.json:318](/home/sano/work/forge/case/56.gap_tp1187/acceptance.json:318) の比は生の `qwall` です。観測モデルは空間的に熱を移すため、比を取っても相殺されません。次の比較対を現行 `skin_forward3d.py` で再計算しました。

   - L15: `case/56.gap_tp1187/run_0052_ab3d_open_long/`
   - L30: `case/56.gap_tp1187/run_0064_ab3d_L30_long/`

   | 比較量 | TC92 | TC91 | TC90 |
   |---|---:|---:|---:|
   | 報告された生熱流束比 | 0.320 | 0.0236 | 0.0122 |
   | 現行観測モデル適用後の比 | **0.6342** | **0.1483** | **0.01969** |

   後者の6枚の系列は `check_quasisteady.classify_series` で各列 **`STEADY`** でした。実測との不一致は残りますが、「3点とも R≤0.35」という原因分類は崩れます。  
   **対案:** D4 も測定と同じ観測量で比較し、「読み C で確定」を撤回する。上記再計算値も、指摘2・3を直すまで暫定値とする。

5. **Major — D4 の比較対が壁解像ゲートを満たしていない。**  
   両 run に `check_wall_resolution.py --target 1 --over-frac 0` を実行した結果は、いずれも **`VERDICT: FAIL`** でした。すきまの最大 y₁⁺ は L15 が **175.6**、L30 が **220.1**。平板では評価点の **97.8%／98.1%** が1を超えています。  
   この割合は現行ツール内部では点数割合であり、面積割合ではありません。それでも超過許容ゼロの契約に対する不合格は変わりません。  
   **対案:** 比較領域の局所壁解像を満たし、実面積による超過割合を報告する。未達の格子対から輸送機構やモデル限界を確定しない。

6. **Major — 2格子の外挿値を、格子誤差の上限として扱っている。**  
   [plan:919](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:919) は、仮定した次数 p=1／2 の外挿値から「どの単調収束でも実測に届かない」としています。しかし2点から次数も漸近域も決まりません。同じ TC90 の値 0.060／0.098 と格子比約1.4146でも、p=0.04なら外挿値は約2.82です。これは次数の推奨ではなく、単調性だけでは上限を証明できない反例です。  
   また [acceptance.json:218](/home/sano/work/forge/case/56.gap_tp1187/acceptance.json:218) は、D4 不合格時に格子系列を追加してから結論することを明記しています。  
   **対案:** 0.19／0.14を「仮定次数による参考外挿」と表記し、登録した追加格子確認を実施する。それを行わない場合は格子依存を未解決として残す。

7. **Major — 給気形状の簡略が未検証なのに、モデルクラスを原因と断定している。**  
   [gen_mesh_gap3d.py:41](/home/sano/work/forge/case/56.gap_tp1187/gen_mesh_gap3d.py:41) は、縦すきま入口のエッジ半径を省略しています。原報の r=2.5 mm は W=1.8 mm より大きく、微小な省略ではありません。これを検証せず、取り込むガスが冷たいという結果から「定常 RANS の枠外」へ進むことはできません。  
   また [tt_audit.py:70](/home/sano/work/forge/case/56.gap_tp1187/tools/tt_audit.py:70) が出す壁手前の帯内**中央値**は、供給ガスの最大全エンタルピーや熱流束上限の証明ではありません。  
   **対案:** 結論を「現在の簡略形状・設定では過小予測」に限定する。半径を含む給気形状の感度を確認するまで、モデルクラスの否定と case/49 への「下界」の受け渡しを取り下げる。

8. **Major — 準定常判定を、未達の収束条件の代わりにしている。**  
   上記 D4 比較対を再判定すると、両方とも残差は **`NOT CONVERGED (stalled/plateau)`**、浅部熱流束・壁圧系列は **`OVERALL: ALL STEADY`** でした。保存された中間・最終体積場には NaN/Inf はありません。これらは別々の事実です。  
   [acceptance.json:252](/home/sano/work/forge/case/56.gap_tp1187/acceptance.json:252) の「事前登録どおり派生量で判定」という説明は、同ファイルの収束に関する失格条件と plan G9 に整合しません。G8 の物理収支診断も [plan:1004](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:1004) で未実施です。  
   **対案:** 残差・準定常・収支の判定を別欄に置き、未達を明記する。判定契約を変更する場合は事後変更として記録し、従来契約の達成とは扱わない。

9. **Major — `--dry` を含む再実行で既存 run の成果物を削除する。**  
   [make_case3d.py:149](/home/sano/work/forge/case/56.gap_tp1187/tools/make_case3d.py:149) は既存ディレクトリを受け入れ、`res_*`・ログ・残差履歴を削除してから、後段で `--dry` を判定します。run 名の再指定で検証根拠を失います。  
   **対案:** 計算済みディレクトリが存在したら削除前に失敗させる。継続は新規 run と `extend_run.py` に限定する。

10. **Minor — 残作業表と現在状態の同期が完了していない。**  
    [plan:1009](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:1009) は、既に存在する前壁置換実装や `acceptance.json` を未作成としています。3D の行には終了宣言と、既に実施した投入・事前登録を「未了」とする記述が同居しています。[README:463](/home/sano/work/forge/case/56.gap_tp1187/README.md:463) では `run_0060` が実行中・完走の2行に重複しています。  
    **対案:** 現在の判定・残作業を一つの表に整理し、過去の判断は日付付き履歴へ移す。`plans/README.md` の到達点も同期する。

**推奨は、plan を `active` に残し、比較手順を修正して再レビューすることです。** 不合格という結果自体は成果になり得ますが、現在の実装と比較では、その原因や下界まで確定できません。

指定差分と関連実装、原報PDF、ローカル保存場を確認しました。`run_0061_prod3d_cfl4` はローカルに存在せず、AWS 生産場の直接再検証はできていません。ファイルは変更していません。

指摘数: Critical 0 / Major 9 / Minor 1
