# codex 諮問 (diagnose): sern-judgment-path

- **brief**: [`notes/reviews/briefs/2026-10-12-sern-judgment-path.md`](../../notes/reviews/briefs/2026-10-12-sern-judgment-path.md)
- **plan**: [`plans/accepted/tooling-sern-te-wake-grid.md`](../../plans/accepted/tooling-sern-te-wake-grid.md)
- **date**: 2026-10-12
- **commit**: `fa816325` (feature/sern-design)
- **codex**: effort `high`, 4.7 min, rc=0
- **結論**: **次は R2 の受入れを止め、必須判定量・標本・判定結果を問題設定から固定して、`collect` から `_learnable` まで欠測・未達を拒否する経路を完成させる。**
- **extra**: `notes/reviews/2026-10-11-tooling-sern-te-wake-grid-result.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（**現時点では受入れ・`status: done` は不可**）

| 対象 | 採否 | 理由・対案 |
|---|---|---|
| 1. 同じ run ディレクトリ・別 session | **採用** | 床記録は起動ごとの session を追記し、最後の session を判定する仕様。別ディレクトリ化は必須でない。ただし既存の検証 run を上書きせず、新規 run 内の段として使う。 |
| 2. index コピー＋一致検査 | **採用、検査は修正** | 同一格子の index コピーは妥当。ただし現在の検査は SST 保存量の欠落を検出しない。下記 Minor。 |
| 3. 共通窓 W・N ≥ 2W・20標本以上 | **採用、実装は要修正** | 設定の制約は妥当。しかし「必要な列・予定した標本が揃ったこと」を保証していない。下記 Major 1。 |
| 4. 最後の `metrics.json` が正本 | **採用** | 延長前の記録を退避する実装と整合。ただし学習側にも同じ判定要件が必要。下記 Major 2。 |
| 5. 未定常だけなら一度延長 | **採用** | 低温・床・欠測・発散を延長で救済しない方針は妥当。`judgment_extendable` は先頭の `fail_class` だけでなく個別ゲートも検査している。 |
| 6. 延長後の未達は判定不能 | **採用** | 生産評価に採らない判断は妥当。ただし旧行の取り込みまで含む除外は未完成。 |
| 7. 固定領域に合わない形状は判定不能 | **採用** | 今回の固定形状の受入れに限定する。一般の設計探索が使えるという主張は不可。 |
| 8. 2D は判定不能 | **採用** | 2D の独立検証は移管済み。本 plan を閉じるために実装範囲を戻す必要はない。 |
| R1・R3 | **修正方針を採用** | 床要件をキャンペーン側に置くこと、`L_sw` の値を方式 ID から外すことは妥当。ただし R1 だけでは判定区間全体の学習要件にならない。 |
| Q2：偽 forge の延長 PASS だけで受入れ | **却下** | 実 GPU の成功側経路が未検証。下記 Major 3。 |
| Q4：result 再レビュー | **修正・追加確認後に採用** | 現状のまま再提出しても NO-GO。以下を §5.1 R2・§6 に反映してから再提出する。 |

**Major 1 — 必須の判定量が欠けると、検査対象を減らして PASS する。**

根拠：[sern_gates.py:339](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:339) は、最後の履歴に `C_T_with_shear` がなければ3量へ切り替えます。

保存された `case/46.sern_design/run_1096_tewake_lsw10_m10_cont20k/` の履歴をメモリ上で再判定し、次を再現しました。

- 元の4量：`judgment_gate = PASS`。
- **`C_T_with_shear` 列だけ削除：同じく PASS**。検査対象は `C_T / C_L / C_M` の3量。
- 元履歴への `check_quasisteady.py --tail 0.5` 再実行は **`OVERALL: ALL STEADY`**。

さらに [同ファイル:333](/home/sano/work/forge-sern-design/design/forge_design/metrics/sern_gates.py:333) は標本数しか検査しません。合成入力で、各10000 step窓の**後半5000 stepだけに20標本**を置いても PASS しました。

**対案：** 必須量は `evaluate.model` 等の設定から決め、SST は4量を必須にする。予定 step 列を N・W・出力間隔から作り、欠測・重複・終了 step 不一致を拒否する。Euler の3量への切替は、列の欠落ではなくモデル指定で行う。

**Major 2 — 追加した温度・窓の判定が、旧行の学習採否には接続されていない。**

根拠：[driver_sern.py:203](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:203) の `_learnable` は床判定までしか照合しません。`judgment` は [同ファイル:405](/home/sano/work/forge-sern-design/design/forge_design/opt/driver_sern.py:405) で保存しますが、採否には使いません。

`case/46.sern_design/run_1097_tewake_lsw10_m6_cont20k/` の保存記録では、

- 旧総合判定：`GATES PASS`。
- 床：**PASS、(10000, 20000]、10000更新**。
- 再計算した C_M の d：**1.9722766e−3 > 5e−4**。
- 温度 CSV への `check_quasisteady.py` 再実行：**`TRANSIENT-UNSETTLED`、`OVERALL: NOT ALL STEADY`**。

この記録から作った作動点要約に `judgment: UNSETTLED` を明示しても、方式 ID を照合条件に合わせると **`_learnable == True`** でした。

現在は旧 ID の不一致が取り込みを止めています。したがって「既に現行系列へ混入した」とは断定しません。しかし、**(d) の移行・持ち越しで露出する穴**です。

**対案：** 判定区間を要求するキャンペーンでは、各作動点に #1〜#3 の PASS・判定窓・必須量の充足を要求する。旧 ID の承認だけでは解除せず、保存系列の再判定を必須にする。

**Major 3 — GPU 接続試験の PASS は、延長後の正常終了を保証していない。**

根拠：[run_sern_judgment_connection_test.py:181](/home/sano/work/forge-sern-design/design/tests/run_sern_judgment_connection_test.py:181) は、2回目が完走した場合だけ完了 session・監視行数を検査します。発散時は「FAIL として処理できた」ことで試験を通します。

これは失敗側の接続確認として有効です。しかし、**「発散は流れだけの問題で経路とは無関係」は未確認**です。

**対案：** 受入れ条件に、実 GPU・本来の N=20000／W=10000で、強制分岐なしの「未定常のみ → 自動延長 → PASS」を1件追加する。保存した本段最終場が使えれば、段階起動からの再実行は不要です。延長後の床 session、#2・#3、残差 VERDICT、退避記録と最終 metrics の対応を残す。これは今後の受入れ条件として事前登録してください。

**Minor — 「保存量のビット一致検査」という説明は実装より強い。**

根拠：[runner_sern.py:1754](/home/sano/work/forge-sern-design/design/forge_design/evaluate/runner_sern.py:1754) は両ファイルの共通キーだけを比較し、必須検査は流体5保存量だけです。メモリ内 HDF5 で、転送先から `roK / roOmega` を欠落させても `_restart_mismatch = []` を再現しました。また `np.array_equal` は +0 と −0 を区別しません。

**対案：** 設定から必須保存量集合を決め、存在・shape・dtype・バイト一致を検査する。通常の生産 run で転送誤りが起きた証拠ではありませんが、現在の検査では除外できません。

残る穴の扱いは次の一案に絞ります。

| 穴 | 受入れ前の扱い・移管先 |
|---|---|
| (a) 形状からの監視領域生成 | 固定形状限定を明記して後回し。一般化は2D plan と共通化し、親 R7b の探索再開条件に置く。 |
| (b) 2D driver | 2D plan。今回の3D受入れを止める理由にはしない。 |
| (c) `L_sw` が x にない | 親 R7b。異なる `L_sw` を同じ学習器へ入れる前に必須。ID統合だけで入力変数不足は解決しない。 |
| (d) 旧方式 ID | 移行は親 R7b。ただし**再判定なしの取り込み拒否は今回塞ぐ**。元記録を保存し、旧レシピからの変換と #1〜#3 の再判定を別記録にする。現在のコードで prepare し直すだけでは過去の来歴を証明できない。 |
| (e) 環境依存 | 起動契約・実際の環境値の記録は受入れ前に必要。`FORGE_ALLOW_UNVERIFIED_SPECIES=1` の無条件既定化はしない。[solver-settings.md:324](/home/sano/work/forge-sern-design/procedures/solver-settings.md:324) の実行単位の許可に従う。 |

結論: **次は R2 の受入れを止め、必須判定量・標本・判定結果を問題設定から固定して、`collect` から `_learnable` まで欠測・未達を拒否する経路を完成させる。**

第 1 仮説: 主な残存リスクは格子や物理ではなく、必要な証拠が欠けても判定対象を減らす／旧 PASS を採る採否経路にある。 確度: **高**  
  根拠: 上記の列削除による PASS と、`run_1097` 相当行の `_learnable == True` を再現。  
  反証条件: 設定上必須の証拠を欠落させた入力が、判定・台帳・学習の全経路で拒否されること。

第 2 仮説: GPU 延長時の発散には restart 経路の問題が含まれる可能性が残る。 確度: **低・未確認**。保存量コピーの検査だけでは、ソルバ起動後の状態再構成まで除外できない。

判別 A/B: **0 step の保存データ試験**。A＝`run_1096` の元履歴、B＝同じ履歴から `C_T_with_shear` 列だけを削除。SST設定・温度系列・窓は固定する。→ **Aだけ PASS、Bが判定不能なら列欠測の穴は解消／両方 PASS なら第1仮説を支持**。現行実装では後者を再現済みなので、修正後の受入れ試験にする。

やらない方がよいこと: 窓の短縮、格子の再調整、旧 ID の一括付け替え、発散した接続試験を成功側経路の証明として扱うこと。

呼び出し側の前提への異議: 「134 PASS／16 PASS」は登録した試験の通過であり、必須量の充足や延長後 PASS の独立実証ではない。また、方式 ID の統一と、同じ学習入力で扱えることは別問題である。

不足情報: 実 GPU の自然な延長後 PASS の記録と、その起点の保存場・バイナリ・環境の対応。今回確認した12本は保存 CSV の再判定であり、GPU実行の再現ではない。ファイル変更・forge起動はしておらず、**本判断は plan 未反映**。
