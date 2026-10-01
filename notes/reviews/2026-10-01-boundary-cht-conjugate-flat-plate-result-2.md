# codex レビュー: boundary-cht-conjugate-flat-plate (result)

- **plan**: [`plans/active/boundary-cht-conjugate-flat-plate.md`](../../plans/active/boundary-cht-conjugate-flat-plate.md)
- **stage**: `result` (diff base `feature/cht-axisym-graetz`)
- **date**: 2026-10-01
- **commit**: `deee16e8` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 6.2 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M0/m3
- **extra**: `notes/reviews/2026-10-01-boundary-cht-conjugate-flat-plate-result.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**

n32/n64 の限定結果は支持できます。前回の Major 2 件は実質的に解消されています。
ただし、上境界延長の端点と、閉鎖時の文書・残作業表に Minor 3 件が残っています。

1. **Minor — 「上境界 1.5 倍」の実装が、登録した端点に届いていない**

   **根拠:** [eval_conj.py:99](/home/sano/work/forge-cht/case/64.conjugate_pipe_wall/eval_conj.py:99) の `np.arange` は端点を切り捨てます。実際の延長高さは n32 で **1.4663708 × H_TOP**、n64 で **1.499683 × H_TOP** でした。

   書き込みなしで端点を正確に `1.5 * H_TOP` に補い、C1/C2 × n32/n64 を再計算すると、感度は最大で **θ 約 7.4e-13、q/q_mean 約 6.4e-11**。結論への影響はありません。

   また、[plan:31](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-flat-plate.md:31) と変更ログの「θ・q とも 1e-11 以下」は保存結果と不一致です。`case/65.conjugate_flat_plate/run_0023_c2_n64_lim0/EVAL_CONJ.txt` の q 感度は **6.8568e-11** です。

   **対案:** 最終節点を登録端点へ一致させ、実際の延長高さを出力する。延長参照の評価記録・台帳を更新し、要約は「θ 約 1e-12、q/q_mean 約 1e-10 以下」とする。流体計算の再走は不要です。

2. **Minor — 限定閉鎖の判断と、現在の完了条件・索引が揃っていない**

   **根拠:** [plan:19](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-flat-plate.md:19) は「完了条件は変えない」、[plan:258](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-flat-plate.md:258) は依然「6 本すべて合格」を現在の完了条件として掲げています。一方、§4.6.1 と結論は C2 n16 FAIL を保持した限定閉鎖です。

   [plans/README.md:27](/home/sano/work/forge-cht/plans/README.md:27) も `draft`・旧結果のまま。[case README:3](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/README.md:3) の「登録文そのまま」は事後改訂と整合しません。

   **対案:** 元の条件は「旧登録条件・未達」として保存し、現在の閉鎖根拠を §4.6.1 の限定判断へ明示的に接続する。移動時に索引、case README、`methods/boundary.md` のリンクを同期してください。C2 n16 FAIL と旧 `limiter: 2` の判定不能は維持します。

3. **Minor — `limiter: 2` の未解決原因が残作業表に独立して残っていない**

   **根拠:** [plan:226](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-flat-plate.md:226) は丸め分離と `limiter: 2` の停滞原因を保留していますが、[§5.1:252](/home/sano/work/forge-cht/plans/active/boundary-cht-conjugate-flat-plate.md:252) の保留行は④の丸め分離だけです。停滞原因の調査は、完了済み診断の記録に埋もれています。

   **対案:** §5.1 に「生産用リミッタ経路の停滞原因：保留・限定閉鎖の対象外」を追加し、既に決めた再開トリガ、担当、根拠 `run_0005〜0010`・B3 診断への参照を記載してください。

照合できた結果は次のとおりです。

- [run 一覧](/home/sano/work/forge-cht/case/65.conjugate_flat_plate/README.md:9) の `case/65.conjugate_flat_plate/run_0022〜0027` を対象に、圧縮残差から公式 `check_convergence.py` の判定を再実行し、全件 **`PASS (converged)`** を再現しました。判定区間は各再走単独の step 0〜149999。能動列は下降中で、床への到達を示す結果ではありません。
- 保存された全節点系列を `check_quasisteady.py` の公式関数で再判定し、全件 **`STEADY`**。G-if は登録閾値で **`PASS`**、旧閾値では④のみ **`NOT CONVERGED`** を再現しました。保存 G-cons は全件 **`VERDICT: PASS`**。最終場の `VALUE/*` に NaN/Inf はありません。
- C2 の追加評価をファイル保存なしで再実行しました。

  | run（`case/65.conjugate_flat_plate/` 配下） | Δθ_ax | U_ax | D_f − U_θ | 主判定 |
  |---|---:|---:|---:|---|
  | `run_0027_c2_n32_lim0` | 0.013032 | 0.0002680 | 0.012372 | PASS |
  | `run_0023_c2_n64_lim0` | 0.013457 | 0.0001274 | 0.013322 | PASS |

  登録した固体効果を検出できています。C2 n16 の保存された4水準判定も **`FAIL`（2.3536 %＋0.97918 %＞3 %）** で、限定閉鎖の記述と整合します。
- 3格子の品質を再判定し **`VERDICT: PASS`**、参照解の自己検査も **`VERDICT: PASS`**。台帳の手元回収分198件はハッシュ一致し、圧縮残差7件も展開後の元ハッシュと一致しました。
- 指定 diff に CUDA/C++ 本体、`design`、既定設定の変更はありません。共有準定常ツールの既定動作は旧版との200系列比較で一致しました。今回の証拠が支える範囲は **FP64・node・登録した平面条件・凍結流れ場参照**です。FP32、cell、周期 seam の精度検証には拡張できません。

**推奨は、上記 1 → 2 → 3 の順で補正し、既存 run を維持したまま n32/n64 の限定結果として `accepted` に移すことです。** ソルバ本体の変更や全面再走を要求する根拠はありません。

ファイル変更なし。指摘・提案は plan 未反映です。

指摘数: Critical 0 / Major 0 / Minor 3
