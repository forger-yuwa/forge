# codex レビュー: tooling-rerun-conditions (plan)

- **plan**: [`plans/active/tooling-rerun-conditions.md`](../../plans/active/tooling-rerun-conditions.md)
- **stage**: `plan`
- **date**: 2026-10-06
- **commit**: `e2210e0e` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.2 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M7/m1

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
形状を固定した条件変更ツールは必要で、campaign 層との分離も妥当です。  
ただし、保存量の扱い、Pt スケーリングの適用条件、既存 runner への引き渡し、検証ゲートは実装前に修正が必要です。

確認範囲: コード・仕様・既存計画を照合し、後処理関数の合成入力試験を実行しました。参照 `run_0094`・`run_0104` の実データはこの作業ツリーに存在せず、`check_convergence.py` は両方に `NO residual_history.csv`、全体に `OVERALL: CHECK FAILURES ABOVE` を返しました。以下の既存 run 数値は台帳の記録であり、今回の再測定ではありません。

1. **Major — 「9 保存量」と `restart_field` の OK だけでは、正しい初期場を保証できない。**

   根拠: [plan:55](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:55) は固定の「9 保存量」を前提とします。しかし [restart_field.py:81](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/restart_field.py:81) は宛先に存在する項目だけを走査し、SRC にない項目は据え置きます。最後の OK は「コピーした項目」の一致であり、必要量の完備を保証しません。3 種以上・トレーサ付きでは必要量が増え、追加の `roY*` や `roXi` をスケールしなければ Y・ξ が変わります。

   また `d_new/d_ref=f` はゼロ成分で破綻します。float32 の `[0,1,−2]` を 0.8 倍する実行確認では、場は全要素有限でも比は `[NaN,0.8,0.8]` でした。

   **対案:** config から必要保存量集合を決定し、SRC・DST 双方の存在、shape、有限性、ρ>0、組成範囲をコピー前に検査する。全対象量をスケールし、`d_new ≈ f·d_ref` をゼロ対応の誤差基準で比較する。未知の輸送量は拒否する。欠落保存量、3 種、トレーサ、ゼロ運動量の試験を追加してください。

2. **Major — Pt だけの変更を「Euler では厳密な相似」とするのは、境界条件込みでは誤り。`stages=none` の一般化も早い。**

   根拠: [plan:56](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:56) は全保存量を f 倍しますが、背圧は `--Ps` 指定時しか変えません。f=0.8 なら Ps/Pt は **1.25 倍**になります。[boundaryCond_d.cu:610](/home/sano/work/forge-integ-1005/solver_density_cuda/cuda_forge/boundaryCond_d.cu:610) は超音速流出では背圧を無視する一方、亜音速・逆流では指定 Ps を使います。したがって相似なのは状態変換であり、任意のノズル境界値問題の解ではありません。

   [plan:63](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:63) の「±25% 以内なら段階起動なし」は、粗格子・減圧側 1 点では裏付けられません。圧力・密度フロアの作動も相似を壊します。

   **対案:** 「T・U・Y を保つ初期場変換」と説明を訂正し、背圧はユーザー指定どおり保持する。Pt スケーリングは明示指定の opt-in、条件変更時の段階起動を既定にする。`none` は同条件継続または検証済み条件に限定してください。

3. **Major — 複製許可リストに、実効入力の依存ファイルを確認する契約がない。**

   根拠: [plan:48](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:48) は `inlet_profile_<physID>.csv` をコピーしません。しかし [boundaryCond.cpp:284](/home/sano/work/forge-integ-1005/solver_density_cuda/boundaryCond.cpp:284) は `inletProfile: 1` なら同ファイルを必須とします。単に追加コピーしても、CSV が設定する Pt/Tt/Y に対して YAML の floats だけを変更する意味が未定義です。

   同様に、既存の `valueFileName` が別の res を指す場合、`NEW/nozzle.h5` の編集が実行に使われる保証がありません。X 組成入力へ Y を追記すると、[speciesDB.cpp:972](/home/sano/work/forge-integ-1005/solver_density_cuda/input/speciesDB.cpp:972) の X/Y 混在拒否に当たります。

   **対案:** v1 の対応入力を明文化し、未対応の profile、複数入口、外部参照、X 入力形式は作成前に拒否する。`meshFileName`・`valueFileName` が新 run 内の意図したファイルへ解決されることを検査する。バイト一致はこの契約を満たす入力に限定してください。

4. **Major — 既存の `run_staged_ns` に委ねるだけでは、段階起動の検証記録が成立しない。**

   根拠: [runner_axismach.py:1012](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1012) の `_stage` は終了コードと最終出力番号を確認し、restart 後に `res_*` を削除します。`run_staged_ns` 内には `StageManifest` の記録も NaN/準定常ゲートもありません。さらに [同:1025](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:1025) は `convMethod: 1` だけを文字列置換するため、2 なら soft 段でも高次のままです。

   計画の CFD 試験はすべて `stages="none"` なので、Tt・組成変更で推奨する `full` 経路を検証しません。

   **対案:** runner の対応設定と段間記録を整える作業を前提依存として §5 に追加する。各段の実効設定・残差・診断結果を保存し、失敗時には次段へ進めないことを確認する。Tt/Y/Ps の変更についても、準備から実効入力・段階起動までの結合試験を追加してください。

5. **Major — 旧 Pt の Euler 参照を使う δ_E 評価には、圧力変更と抽出アルゴリズムの交絡がある。**

   根拠: [deltastar.py:220](/home/sano/work/forge-integ-1005/design/forge_design/metrics/deltastar.py:220) は q_NS/q_E をフィットし、[同:248](/home/sano/work/forge-integ-1005/design/forge_design/metrics/deltastar.py:248) でその傾きの**絶対値**を固定閾値と比較して抽出帯を変更します。この判定は、NS の密度だけを定数倍したとき不変ではありません。

   合成入力 `r∈[0,10]`、q_E=1、q_NS=1−exp[−(10−r)]、Euler 壁半径 9.8 で既存関数を実行すると、流速分布形状を変えず q_NS だけ 0.8 倍しても、δ は **0.920740 → 0.883294（−4.07%）**になりました。両方を 0.8 倍すると 0.920740 を保ちます。これは実 run の結果ではなく、抽出関数の反例です。

   **対案:** δ_E 用 Euler 参照は Pt 変更時も条件を合わせる。`mdot_ratio_vs_euler≈0.799` は旧 Pt 基準の診断量として分離する。推測による δ∝Re⁻¹/⁵ の帯や M の増減方向は、ツールの合否条件から外してください。

6. **Major — `STEADY` の閾値が、比較許容差より桁違いに緩い。**

   根拠: [check_quasisteady.py:542](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_quasisteady.py:542) の既定値は drift **5%**、fluctuation **10%**です。計画は具体的なオプションを指定していません。

   13 点の合成系列 Mᵢ=6+10⁻⁴i を既定設定で判定すると **`STEADY`**でした。末尾窓の変化は **5×10⁻⁴**で、計画の A−B 許容差 **5×10⁻⁵の10倍**です。両腕が同方向に動いていれば、平均差だけでも通り得ます。

   **対案:** δ_E・M・流量それぞれについて、評価列、窓、`--drift`・`--osc`、末尾幅、直前窓との平均差を事前登録する。時間変動の幅を比較許容差の内側に収める。プラトー run は `NOT CONVERGED` として保持し、再実行の再現性試験と物理解の収束検証を分けてください。`--from-floor` も未収束参照には使えません（[check_convergence.py:359](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/check_convergence.py:359)）。

7. **Major — 粗格子を固定する試験に、根拠なく壁解像 PASS を必須としている。**

   根拠: [case README:77](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:77) の `run_0094_ns_c2pin_pass2_ext6k` は、**y₁⁺>1 が29.4%、壁解像 FAIL**と記録されています。元 plan の [変更ログ:164](/home/sano/work/forge-integ-1005/plans/active/tooling-nozzle-cfd-pinned-initial-line.md:164) には最大 **13.4**とも記録されています。それなのに本計画は同じ格子で Pt を 20% 下げる試験へ壁解像 PASS を課しています。

   PASS に変わる根拠がなく、実装が正しくても既知の格子不足で失格になる検証です。

   **対案:** 粗格子は初期場・入力変更・再実行の再現性を測る用途に限定し、既知の FAIL を明示して保持する。物理結果の受入れには、壁解像を満たす参照 run を別に選び、同じく形状・格子を固定して試験してください。

8. **Minor — 文書と run 台帳を維持する作業が完了条件から漏れている。**

   根拠: [plan:7](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:7) の `related_docs` に `methods/` がなく、[実装ステップ:68](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:68) に仕様更新と case README の同期責任がありません。[AGENTS.md:166](/home/sano/work/forge-integ-1005/AGENTS.md:166) は実装前の仕様更新、[同:59](/home/sano/work/forge-integ-1005/AGENTS.md:59) は run 作成・破棄時の台帳同期を要求します。また §1 の「乾き成分の組成変更」は、この参照構成では v1 が拒否する lump 変更に当たります。

   **対案:** `methods/design/overview.md` に対応入力・保存量変換・拒否条件を追記する工程と、台帳更新の担当を §5.1 に追加する。§1 は「輸送種の定義を維持した入口分率変更」と限定してください。

**推奨は、v1 を対応入力を厳密に検査する run 準備ツールとして実装し、Pt スケーリングと段階起動省略の自動判断を既定から外すことです。** `plans/README.md`・関連 accepted plan の確認範囲では、この準備機能が既に完成済みという重複は見当たりません。ソルバを変更しない構成は妥当で、現行の node 主体の検証方針にも合います。cell・周期への一般保証は付けないでください。

実装前の修正優先順は、①必要保存量と入力依存関係、②スケーリング条件と runner 契約、③Euler 参照・準定常・壁解像の検証設計、④仕様文書と台帳責任です。ファイルは変更しておらず、**plan 未反映**です。

指摘数: Critical 0 / Major 7 / Minor 1
