# codex レビュー: tooling-nozzle-mknot-fast (plan)

- **plan**: [`plans/active/tooling-nozzle-mknot-fast.md`](../../plans/active/tooling-nozzle-mknot-fast.md)
- **stage**: `plan`
- **date**: 2026-10-15
- **commit**: `528c1c4d` (feature/mknot-fast)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 別 (依頼側 `claude-opus-5-5`、実効 `gpt-6-astra`)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 3.9 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-15-tooling-nozzle-mknot-fast-plan-jit.codex.log` (codex: ok)
- **判定**: **GO-with-changes**, 指摘 C0/M2/m2
- **focus**: §4.4 (網の充填の JIT 化・L_c の使い回し) と §6 V5・V6 の受入の妥当性。他の節は前回までのレビューで見た

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

判定: GO-with-changes
JIT 化の方向は妥当です。ただし、L_c の独立性は一般には成立せず、現行 V6 では評価順依存を検出できません。
実装前に、再利用の適用範囲・決定性と、JIT 側の異常系・診断の受入を明文化してください。

1. **Major — L_c の再利用は、長さ残差だけでは元の評価と同等にならない。**

   根拠: [runner_axismach.py:649](/home/sano/work/forge-mknot/design/forge_design/evaluate/runner_axismach.py:649) は初期推定から反復し、残差が許容内に入った点を採用します。初期値を変えると、別の L_c でも受理できます。また、許容窓自体が M_K に依存します（[axis_law.py:287](/home/sano/work/forge-mknot/design/forge_design/geometry/axis_law.py:287)）。

   **実測反例**: `case/70.m5_air_wt_demo/problem_m5_air_ns2_euler.yaml` を基に、`N=600`、`wall_mode=streamline`、`L_total_tol=0.02` として確認しました。未実装の seed 入力を試すため、既存の `explicit` 経路で再利用時の初回評価を模擬しました。ファイルは変更していません。

   | M_K=2.10 の評価方法 | L_c [r_t] | x_F−L_total [r_t] |
   | --- | ---: | ---: |
   | 通常の `from_length` | 21.497084419019913 | −1.0867460×10⁻⁵ |
   | M_K=2.05 の L_c を再利用 | 21.497099185457813 | +3.8871581×10⁻⁶ |

   両方とも長さ許容差内ですが、**L_c の相対差は 6.869×10⁻⁷**。V6 の 10⁻⁹ を超えます。これは生産の `cplus` が失敗した証拠ではありませんが、§4.4 の無条件な独立性の主張への反例です。

   さらに、プロセス内キャッシュでは各 worker が受け取った先行評価が seed になります。現在の [mknot_select.py:426](/home/sano/work/forge-mknot/design/forge_design/evaluate/mknot_select.py:426) には、worker と親プロセスの両方で評価する経路があり、5 点を一順序で比較する V6 だけでは決定性を保証できません。

   **対案**: 再利用を検証済みの `cplus`／`characteristic` 条件に限定し、固定した代表評価から seed を作って共有する仕様にしてください。適用外・seed 不正・初回検査失敗は元の初期推定から再計算します。V6 に昇順・降順・ランダム順、worker 数 1/4、設定変更による無効化、最終選択値と壁形状の比較を追加してください。数値許容差は量別の絶対＋相対とし、時間・反復回数は比較対象から明示的に除外します。

2. **Major — V5 の「分類件数一致＋既存試験通過」では、JIT に移す診断の正しさを十分に検証できない。**

   根拠: [moc_inverse.py:98](/home/sano/work/forge-mknot/design/forge_design/geometry/moc_inverse.py:98) は最終残差と棄却位置を保持し、[moc_inverse.py:802](/home/sano/work/forge-mknot/design/forge_design/geometry/moc_inverse.py:802) のゲートは残差上限と、**棄却された対が壁の内側か**を判定します。分類件数だけ同じでも、この判定は変わります。

   既存の異常系試験には、Python の `pm_mach_vec` を差し替えて `interior_vec` を直接呼ぶものがあります（[run_moc_axis_limit_tests.py:176](/home/sano/work/forge-mknot/design/tests/run_moc_axis_limit_tests.py:176)）。`fill_levels` 全体を別の JIT 関数にすると、その試験が通っても JIT 内部の異常処理を通した証拠にはなりません。

   **対案**: 同じ異常入力を**実際の両 backend**に通す受入を追加してください。入力欠損、反復中の非有限、平行特性線、軸下への棄却、反復上限、軸極限の各分岐を含め、対ごとの状態・源項分岐・棄却位置・残差・`moc_gate` を比較します。`_geom_all` は全件保持が必要です。JIT の数値契約には float64、`fastmath=False`、CPG の Newton と semiperfect の補間の区別も明記してください。既存の空間誤差許容差を緩める必要はありません。

3. **Minor — 速度・メモリの受入で、JIT と L_c 再利用の効果を分離する条件が不足している。**

   根拠: [plan:160](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:160) の「1 評価が 3 倍以上」は、JIT コンパイル前後と L_c キャッシュ有無を指定していません。また「呼出回数と同じ比率でピークメモリが減る」という見込みは、既に捕捉を 2 件に制限した [mknot_select.py:241](/home/sano/work/forge-mknot/design/forge_design/evaluate/mknot_select.py:241) からは導けません。呼出回数と同時生存する配列数は別です。

   **対案**: V5 は L_c 再利用を無効にして測定し、初回コンパイル込みとコンパイル済みを分けて報告してください。最終構成では worker 起動・コンパイルを含めて V4 を再実行します。メモリ削減率は予測で確定せず、JIT 化後のピーク実測で並列数の見積もりを更新してください。

4. **Minor — §2・§4.3 が追加スコープと矛盾している。**

   根拠: [plan:25](/home/sano/work/forge-mknot/plans/active/tooling-nozzle-mknot-fast.md:25) は JIT を対象外とし、§4.3 も後回しですが、§4.4 と残作業 #4 は実装対象です。

   **対案**: JIT と L_c 再利用を現行スコープへ移し、実装順を「仕様・受入更新 → JIT 単独の V5 → 再利用単独の V6 → 統合後の V3・V4」に統一してください。

**推奨は、上記を計画へ反映してから、JIT と限定的な L_c 再利用を分離実装することです。** 実装前の優先順位は 1 → 2 → 3 → 4 です。

目的は妥当で、`plans/README.md` と accepted の `mknot-auto`・長さ逆算・軸極限の計画を確認した範囲では、既存機能の重複実装には当たりません。今回の対象は CPU 側の float64 の設計処理であり、CUDA の float32、node/cell、block-DPLUR、周期境界を直接変更しません。壁と診断の同等性を守る変更として、CFD 0 step の受入は妥当です。V3(b) の寸法ゲート未判定は、引き続き計画全体の完了を阻みます。

証拠: 足りる

確認対象は `HEAD=528c1c4d`。既存コード・関連仕様・計画を照合し、上記の MOC 反例を実行しました。`run_mknot_v2_tests.py` は通常実行で `ALL OK`（`--real` は未実行）。今回は計画レビューであり、CFD の収束・準定常性や AWS の既報性能を新たに認定していません。作業ツリーの変更はありません。

指摘数: Critical 0 / Major 2 / Minor 2
