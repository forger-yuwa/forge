# codex レビュー: case-plate-annular-cavity-m5 (result)

- **plan**: [`plans/active/case-plate-annular-cavity-m5.md`](../../plans/active/case-plate-annular-cavity-m5.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-21
- **commit**: `843d6734` (feature/sern-design)
- **codex**: effort `high`, 7.7 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M8/m0
- **focus**: 前回 result レビュー (2026-09-20, NO-GO, M1-M9/m1) の 10 指摘それぞれについて、対応が実際に入ったか・主張が裏付けられているかを最優先で確認してほしい。特に (a) TP 前駆への切り替えと熱回路 6 係数の TP 再同定 (§4.3.1 / §4.7.7 / §4.7.11)、(b) ゲート類の厳格化 (check_case_gates.py の終了コード・入力不備・収支残差 JSON 保存・stage_manifest の turbulence.model)、(c) §4.14 の表がまだ CPG 入口系列であることの但し書きが十分か、(d) 偏心係数を run_0427 の収支残差 -7.52% を理由に暫定扱いにした判断の妥当性。
- **extra**: `case/49.plate_annular_cavity_m5/README.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**  
前回指摘の一部は修正されていますが、残差ゲートの誤合格と、伝導・粘性仕事を除外した収支判定が残っています。  
偏心係数を暫定扱いにする判断は妥当ですが、同心係数を確定扱いにする根拠も不足しています。

指定の diff、未コミット差分、前回レビュー、コード、ローカル成果物を確認しました。TP 前駆 `run_0004_precursor_m9_tp` と生産 `run_0419`〜`run_0428` の実データはローカルになく、それらの最終値・VERDICTは独立再検証できていません。反例試験はメモリ内で実施し、ファイルは変更していません。

1. **Major — 継続 run の残差ゲートは、現在も不合格を通す**

   [check_case_gates.py:138](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/check_case_gates.py:138) の継続経路は、取得した終了コード `rc` を使いません。親との比較も共通列だけで、親自身の収束確認を要求していません。

   必須の運動量・エネルギー残差列を欠く入力について、本物の `check_convergence.py` が **終了コード1・判定不能**を返した出力をゲートに渡すと、他の検査を合格に固定した反例で **`GATES: PASS`・終了コード0**になりました。

   実データでも、`case/49.plate_annular_cavity_m5/_local_s2/run_0401_s2_off000_tw20_ext1/` は再判定で次の状態です。

   ```text
   check_convergence.py: NOT CONVERGED (stalled/plateau)
   check_cavity_steady.py: VERDICT: STEADY (全量)
   check_wall_resolution.py: VERDICT: FAIL
   ```

   **対案:** 入力不備・非ゼロ終了を拒否し、通常判定を通過した参照に対する `--from-floor`、または同一設定区間の正規判定へ置き換える。前回 M2 は未解消です。

2. **Major — エネルギー収支ゲートが伝導・粘性仕事を落としている**

   [cavity_eval.py:910](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/cavity_eval.py:910) の `budget_residual_alldepth` は、分子が **壁入熱−対流流入だけ**です。伝導・粘性込みの値は別キー `budget_residual_all` に保存されますが、[check_case_gates.py:221](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/check_case_gates.py:221) は前者を優先します。

   対流のみの残差を0、全項の残差を20%とした入力でも、ゲートが **`PASS`**を返すことを再現しました。欠損・NaNの拒否とJSON保存自体は修正されていますが、検査対象が誤っています。

   [plan §4.8.6](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:1216) の同心 mixA は、壁20.09 W・対流23.04 W・伝導2.85 Wです。掲載の−5.80%は対流のみの差に対応し、同じ正規化で伝導まで加えると**約−11.4%**になります〔掲載丸め値からの計算、粘性仕事未加算〕。

   **対案:** 全流束を含む残差を統一した尺度で保存・判定し、全TP系列を再評価する。**偏心だけを暫定、同心を確定とする区別はいったん撤回する。**

3. **Major — TP 前駆の初期エネルギーが CPG のまま**

   TP物性ブロックとDBコピーは追加されています。しかし [precursor/gen_runs.py:99](/home/sano/work/forge/case/49.plate_annular_cavity_m5/precursor/gen_runs.py:99) は、TPでも
   `roe = P/(gamma−1) + ρ|U|²/2`
   を書き、[同:179](/home/sano/work/forge/case/49.plate_annular_cavity_m5/precursor/gen_runs.py:179) で無条件に呼びます。

   `case_m9.json` の採用物性と `thermoHrefTemp=298.15` で読み戻すと、初期静温は **624.65 K**、静圧は **15.785 kPa**になります。指定値は **216.65 K・5.475 kPa**です。

   **対案:** 生産側と同様に、TPのエンタルピー基準から内部エネルギーを組み直す。既存TP前駆の最終解が誤りとまでは断定しませんが、入口分布の妥当性は起動完走と切り離して検証する必要があります。

4. **Major — `turbulence.model` の修正が特定のYAML表記にしか効かない**

   [stage_manifest.py:59](/home/sano/work/forge/solver_density_cuda/tools/stage_manifest.py:59) は正規表現による抽出です。生成器の二重引用符付きflow形式では層流とSSTを区別しますが、次の両形式では `none` と `sst` が**同一キーになる**ことを再現しました。

   ```yaml
   turbulence:
     model: "sst"
   ```
   ```yaml
   turbulence: {model: 'sst'}
   ```

   **対案:** YAMLを構造として解析し、実効設定を正規化して比較する。前回 M3 は生成器の現行書式についてのみ修正済みで、共有ツールとしては未解消です。

5. **Major — 「熱回路係数の移動は1〜3%」が掲載値と一致しない**

   [同心係数の説明](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:791) と[偏心係数の説明](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:968) にある旧値・新値から再計算すると、次の差です。

   | 係数 | 同心 | 偏心 |
   |---|---:|---:|
   | 外筒−円柱 | +3.1% | +2.7% |
   | 外筒−底面 | **−5.2%** | **−5.5%** |
   | 円柱−底面 | **−9.7%** | **−9.1%** |

   底面の外部流係数だけを例外としても、「1〜3%」には収まりません。§4.14の但し書きにも同じ過小評価が入っています。

   **対案:** 同定入力・係数・再現誤差・旧系列との差を同じ計算から再生成する。「壁間結合が重要」という定性的結論と、係数変化が小さいという定量的主張を分ける。

6. **Major — TP再同定の正本と、取得できるFEM成果物が同期していない**

   plan #38 は同心FEM束をTPで再生成済みとしますが、取得できる [同心 `network.json`:50](/home/sano/work/forge/case/49.plate_annular_cavity_m5/_fem_bc/fem_bc_conc/network.json:50) の出典は依然 `run_0401`〜`run_0408` です。底面係数も **0.00108691 W/K**で、planの **0.00075 W/K**と異なります。

   出力器には `smeared: false` と未加工の説明が入りました。一方、既存の[出力README:16](/home/sano/work/forge/case/49.plate_annular_cavity_m5/_fem_bc/fem_bc_conc/README.md:16) は、現在も「均してある」と説明しています。

   **対案:** 正本のTP入力・判定記録・同定結果をレビュー可能な場所へ同期し、旧束は旧版と明示する。ユーザ決定どおり均し実装は保留で構いませんが、**既存成果物の虚偽の加工説明は残せません**。レビュー中に追加されたcase READMEのTP run一覧は確認しました。

7. **Major — 必須の物理検証が未完で、非定常性の断定も残っている**

   格子不確かさの転用撤回と半割検証の限定は改善です。しかし、対象のM9・形状2に対する格子検証、全周URANSの四条件、塞ぎ3DのBL移送、離散流束収支は未完です。[完了条件](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:1867) を満たしていません。

   また [§4.12](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:1720) は、定常反復の振動を物理的非定常性と断定したままです。`case/49.plate_annular_cavity_m5/_local_off150/run_0205_off150_s14_ext/` の再判定は、

   ```text
   check_convergence.py: NOT CONVERGED (stalled/plateau)
   check_cavity_steady.py: VERDICT: DRIFTING
   ```

   でした。擬似時間の振動から物理的リミットサイクルは確定できません。

   **対案:** 必須検証を対象条件で完了するまで、同心・偏心とも暫定結果とする。失われた感度時系列は復元または再計算し、未解決作業を実在する残作業行へ集約する。現在も #21〜#30 の行がなく、本文の #28〜#30 参照が宙に浮いています。

8. **Major — 準定常の絶対温度許容が、planと実装で異なる**

   [plan §6.4](/home/sano/work/forge/plans/active/case-plate-annular-cavity-m5.md:1835) は温度drift **2 K以下**を要求します。一方、[check_cavity_steady.py:134](/home/sano/work/forge/case/49.plate_annular_cavity_m5/tools/check_cavity_steady.py:134) は `max(2 K, 平均値の2%)` を使います。

   開口温度差が単調増加し、末尾窓で **6 K**動く反例でも、許容が **9.1 K**となり **`VERDICT: STEADY (全量)`**を返しました。また、エネルギー収支の許容も §6.4は5%、実装は7%です。

   **対案:** 合否条件を単一の設定から読み、planと一致させる。許容を変更するなら根拠と対象量を明記し、旧VERDICTをそのまま新条件の証拠に使わないこと。

前回10指摘の対応確認は以下です。

| 前回 | 今回の確認 |
|---|---|
| M1：入口EOS | TP対応は追加。初期エネルギー不整合が残り、TP実データは未取得 |
| M2：継続残差 | **未解消。誤合格を再現** |
| M3：段キー | 現行生成書式では修正。別の有効なYAML表記で再現 |
| M4：収支ゲート | JSON保存・欠損拒否は修正。**全流束を判定しない問題が残存** |
| M5：局所比熱 | **修正確認**。1273.15 Kで1184.809、直接評価1184.811 J/(kg·K) |
| M6：格子不確かさ | 転用撤回を確認。対象条件での検証は未完 |
| M7：半割・非定常 | 限定追記を確認。必須比較と本文の断定修正は未完 |
| M8：リップ加工表示 | 出力器は修正。既存成果物は未更新 |
| M9：残作業・正本同期 | 一部改善。未完行・係数・成果物の不整合が残存 |
| m1：旋回回帰試験 | **修正確認**。本番関数を直接呼ぶ7試験が `VERDICT: PASS` |

**推奨:** `active` に留め、まず残差・収支・準定常ゲートとTP初期化を修正し、TP正本系列を再判定してください。その結果から係数・表・成果物を同期し、対象条件の必須検証を完了した後にresultレビューを再実施することを推奨します。§4.14のCPG入口という但し書きは必要ですが、見出しの「確定値」とも矛盾しており、移管を認める根拠にはなりません。

診断用 `qwall` 追加の確認した分岐には、直接の残差変更を認めませんでした。ただし、ビット一致の実測やGPU回帰試験は実施していません。本レビューは依頼どおり**ファイル未変更・plan未反映**です。

指摘数: Critical 0 / Major 8 / Minor 0
