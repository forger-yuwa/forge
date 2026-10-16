# codex 諮問 (diagnose): settle-v2-span10k

- **brief**: [`notes/reviews/briefs/2026-10-16-settle-v2-span10k.md`](../../notes/reviews/briefs/2026-10-16-settle-v2-span10k.md)
- **plan**: [`plans/active/tooling-settle-convergence-and-monitor.md`](../../plans/active/tooling-settle-convergence-and-monitor.md)
- **date**: 2026-10-16
- **commit**: `cf41a318` (feature/card-local-mesh)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 3.8 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-16-settle-v2-span10k-diagnose-settle-v2-span10k.codex.log` (codex: ok)
- **結論**: #5d に raw 残差の全行被覆検査と壁熱流の幅の定義を補い、固定した `wt_nozzle_v2` で run_0012 の同じ軌道を再生する。
- **extra**: `notes/investigations/2026-10-16-output-host-readers.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（重大度は、記載した問題を残した場合の評価）

| 対象 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 1 万 step・保持なし | **採用** | ユーザ決定を維持する。許容の 1/8 化や保持の復活は提案しない。ただし現行の `settle_presets.py` はまだ `wt_nozzle_v1` のみであり、以下は V2 の実装・検証要件である。 |
| 毎 step の残差判定 | **要再検証・Major** | [settle.py:1024](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/settle.py:1024) は出力格子の行の存在を検査するが、その間の全行の被覆を保証しない。既存関数をメモリ上の一定系列で実行すると、W=10,000 で「途中の 1 行を削除」「2,000 step ごとの行だけ残す」の両方が `SETTLED`。**raw モードでは、現在窓・直前窓とも必要な整数 step が重複なく全て存在することを検査し、欠測は `INDETERMINATE` にする。** |
| 壁熱流の「区間の幅」 | **要再検証・Major** | [settle.py:1434](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/settle.py:1434) は `max−min` ではなく、終端値からの最大差を使う。6 点を `[0.25, 1.75, 1, 1, 1, 1]`、許容を 1.001 とした既存関数の再現では、**幅 1.5 なのに `SETTLED`**、判定比は 0.74925。V2 は節点ごとの `max(Q)−min(Q)` を判定し、相対許容の参照値を明記する。V1 の既存結果は書き換えない。 |
| 停止候補から 8 万 step の検証 | **採用。ただし保証範囲の拡大は却下・Major** | [settle_virtual_stop.py:268](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/settle_virtual_stop.py:268) の固定した t_s からの差は、窓ごとの許容内に隠れる累積ドリフトを検出できる。ただし目的量は 2,000 step ごとの観測であり、8 万 step でも端点込み 41 点。例えば周期 2,000 step の変動は全観測点で同じ値になり得る。**「この軌道・観測間隔・8 万 step で反例なし」に限定し、プリセット全体の `provisional: false` は保留する。** |
| B の解釈 | **対応方針は採用、原因の断定は却下・Major** | [settle_virtual_stop.py:425](/home/sano/work/forge-cardmesh/solver_density_cuda/tools/settle_virtual_stop.py:425) は、目的量の累積超過がなくても残差条件の再逸脱だけで B を返す。したがって **B＝目的量を早く止めすぎた証明、とは限らない**。「(a) 累積許容超過」「(b) 全条件の持続性違反」を分けて報告する。B を合格へ読み替えず、閾値を緩めずユーザに判断を戻す方針は妥当。 |
| V1 (f′)・(g′) | **範囲付きで採用。無限定の合格は却下・Major** | 診断原本を照合し、float の host 差 43・58・57、FP64 の 43・58・59、全ケースの device 差 0 を確認した。farfield の 96 配列と CHT の `Ts` も B/C 差・`s0twin` が 0。CHT の読み直しは [conjugateWall.cpp:930](/home/sano/work/forge-cardmesh/solver_density_cuda/conjugateWall.cpp:930)、farfield の host 読みは [convectiveFlux_d.cu:119](/home/sano/work/forge-cardmesh/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:119) と整合する。**現行コード、`gpu: 1`・node・定常陰解法、試験済みの Euler／B0 NS／CHT `q_eff`／farfield、float・FP64 の範囲で合格としてよい。** dual-time・RK・`q_compact` の動作確認へは拡張しない。 |
| 将来の読み手の監視・原本保存 | **採用・Minor** | [調査記録:141](/home/sano/work/forge-cardmesh/notes/investigations/2026-10-16-output-host-readers.md:141) の網羅性は現行コードに限られる。今回の合格に新規 lint の完成は要求しないが、host 読み手・H2D・出力経路の変更時に分類表と状態不変試験を更新することを維持条件にする。原本は現在 `/tmp` にあるため、合格記録とともに恒久保存する。 |

残差の位相については、定常陰解法は 1 outer step ごとに非線形更新し、出力後に `outer_end` を記録する構造である（[main.cpp:3111](/home/sano/work/forge-cardmesh/solver_density_cuda/main.cpp:3111)、[main.cpp:3156](/home/sano/work/forge-cardmesh/solver_density_cuda/main.cpp:3156)）。`outer_end` のみに限定する方針は妥当。ただし、これは更新後の場から残差を再計算した行ではなく、記録済みの残差を参照する行である。場 s と残差 s−1 の対応、完成行の確認を維持し、内反復行の混入・出力境界での欠落を試験すること。

結論: **#5d に raw 残差の全行被覆検査と壁熱流の幅の定義を補い、固定した `wt_nozzle_v2` で run_0012 の同じ軌道を再生する。**

第 1 仮説: V2 の実装を既存処理の設定変更だけで済ませると、ユーザが決めた「毎 step」「区間の幅」より弱い判定になる。確度: **高**。  
根拠: 上記の欠測系列と壁熱流系列が、既存判定関数でいずれも `SETTLED` になった。  
反証条件: V2 の経路で欠測例が `INDETERMINATE`、幅超過例が `NOT_SETTLED`、完全な一定系列が `SETTLED` となり、オンラインと再生の判定が一致すること。

第 2 仮説: V2 が、この軌道で許容を超える累積ドリフトの前に停止する。確度: **未確認**。V1 の集計から V2 の停止候補や合否は決められない。

判別 A/B: **変更する設定は判定プリセットだけ**とし、run_0012 の元データを `wt_nozzle_v2` で再生する。最初の合格点 t_s から 80,000 step、目的量の t_s からの差／元の許容と、全格子点の全条件を測る。**A なら**この軌道・観測間隔・検証長で早すぎる停止の反例を退ける。**B なら**登録した検証契約を反証し、累積超過か条件の再逸脱かを区別してユーザへ報告する。候補なし・期間不足・反例のない欠測は未成立。GPU の追加実行はこの一手には不要。

やらない方がよいこと: 94,000 step の結果に合わせた閾値調整、6 点の観測から観測間の不変性を主張すること、毒入れ後の場の比較を読み手の網羅の代用にすること、#1c の原本 `OUTPUT_INVARIANCE_VERDICT: FAIL` を PASS に書き換えること。契約上の範囲付き合格は別に記録する。

呼び出し側の前提への異議: 「94,000 の格子点の 1 点」は**空間の節点ではなく判定時刻**である。`run_0012_ns_free_judge/settle_virtual_stop/grid_checks.csv` で確認した逸脱は、`rms_roOmega` の窓の傾き **0.022828284759136806 > 0.021189299069938092**。生残差の単発異常とは断定できない。また、B だけで目的量の早すぎる停止を断定する解釈は採らない。

不足情報: V2 の実データ再生結果は未取得。これは次の検証で得る情報であり、今回の設計判断と V1 の範囲付き採否を妨げない。

証拠: 足りる
