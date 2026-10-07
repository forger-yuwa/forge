# codex 諮問 (diagnose): hostmem-result-interpretation

- **brief**: [`notes/reviews/briefs/2026-10-07-hostmem-result-interpretation.md`](../../notes/reviews/briefs/2026-10-07-hostmem-result-interpretation.md)
- **plan**: [`plans/active/architecture-solver-host-memory.md`](../../plans/active/architecture-solver-host-memory.md)
- **date**: 2026-10-07
- **commit**: `f8c989d4` (feature/sern-design)
- **codex**: effort `xhigh`, 5.4 min, rc=0
- **結論**: ?
- **extra**: `case/66.hostmem_regression/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 諮問事項 | 判断 | 根拠・対案 |
|---|---|---|
| B の比較器修正と集計結果 | **採用** | 保存された TSV は計3402量すべて PASS、厳密一致185量。単体試験8件も再実行して PASS。ただし、並べ替え不変性は比較器の性質であり、回帰の検出力を保証しない。 |
| B 全件 PASS を根拠に R1–R3 を `done` にする | **要再検証〔Major〕** | `S` に new 内の差を含めるため、変更後だけの異常も許容幅を広げる。[`judge_abs`:514](/home/sano/work/forge-sern-design/case/66.hostmem_regression/compare_runs.py:514) を実行すると、人工入力 base=`[0,0,0]`、new=`[0,0,1e23]` は **S=D≈1e23、PASS**。対案は、評価する new から独立した固定幅での確認。 |
| 7量を「ゼロ近傍・診断量なので無害」と扱う | **却下〔Major〕** | `condClampCorrQ_0` は `maxₖ |ΔQₖ|/max(|Qₖ,before|,1e−30)` の最大値で、絶対補正量ではない。[実装:313](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:313)。一方、`roUz` と履歴は保存量。前者は補正前後のモーメント、後者は代表運動量に対する絶対許容幅で評価し、除外しない。 |
| `twall_z` を「float32 の刻みの問題」として決着 | **却下〔Minor〕** | **B の境界上 PASS 自体は正しい**。しかし S=2⁻⁷、D=2⁻⁶という数値だけでは原因を特定できない。[比較報告:13](/home/sano/work/forge-sern-design/case/66.hostmem_regression/results/2026-10-07_base9c9f623c_vs_new93e55957_abs/sern_g3/sern_g3.txt:13)。最悪差の位置・6本の値・その位置の ULP 距離を記録し、「丸めだけ」との断定を外す。 |
| ホストメモリ削減の実測 | **採用** | g3 の約5157→2630 MiBとGPU使用量不変は[工程別記録](/home/sano/work/forge-sern-design/case/66.hostmem_regression/results/2026-10-07_base9c9f623c_vs_new93e55957/sern_g3_memlog.txt:34)で確認できる。1436 B/節点は切片込みの単点値、1371 B/節点はローカル2格子の傾きとして区別する。 |

**結論:** R1–R3 の `done` は保留し、まず `c44dual_ckpt100` について、new を許容幅の算定に使わない固定幅の独立 A/B を一括で行う。

**第1仮説:** 凝縮の外れ値には、既存の非決定性と微小分母による診断量の増幅が寄与している。**確度: 中**。

根拠: `case/66.hostmem_regression/run_0106_c44dual_ckpt100_new_r1/` の `condClampCorrQ_0` が、SとDの両方を約1.462e23にしている。一方、別構成の `run_0077_c44dual_pindiag_base_r3/` にも約1.45e23が記録されている（[run別の整理](/home/sano/work/forge-sern-design/case/66.hostmem_regression/README.md:207)）。前述の分母処理は増幅機構を説明する。ただし、pin診断の有無が異なるため、後者を同一構成の対照反復として混ぜてはいけない。

反証条件: 同一入力・同一診断条件で、独立した base 対照が固定幅を満たす一方、new だけが補正前後の保存量や固定幅を継続的に逸脱する場合、「既存ノイズだけ」という説明を棄却する。

**第2仮説:** ホスト確保・初期化順序の変更に関連する間欠的な差を、B の許容幅が吸収している。**確度: 低、未確認**。カーネルソース不変更と初期出力のビット一致は支持材料だが、この可能性の除外にはならない。

**判別 A/B:** 変更点は **実行バイナリだけ**。`c44dual_ckpt100` の同一入力から、base=`9c9f623c`、new=`93e55957` を各6本、各100 step。順序は事前固定で交互にし、新しい run に保存する。

- 実行前に、既存の同構成 base 3本だけから各量の `S₀=max‖Bᵢ−Bⱼ‖∞` を計算し、診断幅 `T=2S₀` を凍結する。`S₀=0` は差0を要求する。
- 新規各runと既存base各runの差を同じ T で評価する。対象は既存の全196量。初期出力・構造は厳密一致、非有限は失敗。残差の行キーと物理時刻を固定する。
- **A：新規base・newがともに全量で幅内**なら、「追加標本でもこの固定幅を超える差は検出されなかった」。限定した支持材料として result レビューへ渡す。
- **B：新規baseは幅内、newだけ超過**なら、変更起因の差を優先して追う。
- **新規baseも超過**なら、基準3本では再現性を捉えられていない。判定不能とし、その場でTや反復数を増やして合格にしない。

これは**受入れ幅の妥当性を確認する診断**である。物理的に許せる誤差の上限が未定義なので、全量が幅内でも、それだけで包括的な同等性を宣言しない。

**やらない方がよいこと:** Sが増えるまで反復を追加する、7量を除外する、`condClampCorrQ_0` の分母床を変えて今回の差を消す、未収束の100–200 step比較を定常解の一致と呼ぶこと。

**呼び出し側の前提への異議:** 「Bでは区別できなかった」は、限定事項として書くだけで受入れ根拠に転化できるものではない。記録は **「AはFAILを保存／Bは定義した条件で全量PASS／7量の変更起因性は未解決／`twall_z` は境界上PASS・原因未確定」** とする。

段階2では、固定した許容幅・入力・バイナリと対照runを先に確定する。ホストピークが `readMesh` と `setStructuralVariables` 付近で並ぶ実測から、進める場合は **R4＋R5を先に検討し、R6は別の変更・検証として扱う**。Bの方式をそのまま恒久的な受入れゲートにはしない。

**不足情報:** 対象runの生HDF5・残差CSV・`CONVERGENCE_VERDICT.txt` はローカルに無く、独立再集計はできていない。台帳では `NOT CONVERGED`。凝縮補正が最大になった時点の補正前後値、`twall_z` の最悪位置の値、量別の許容誤差上限も不足している。

ファイル変更・forge実行はしていない。**plan 未反映**。呼び出し側で `architecture-solver-host-memory.md` §6.2・§5.1・§8へ反映するための診断として返す。
