# codex 諮問 (diagnose): wall-orthogonality-judgments

- **brief**: [`notes/reviews/briefs/2026-10-14-wall-orthogonality-judgments.md`](../../notes/reviews/briefs/2026-10-14-wall-orthogonality-judgments.md)
- **plan**: [`plans/active/tooling-generator-defaults-and-launch-gate.md`](../../plans/active/tooling-generator-defaults-and-launch-gate.md)
- **date**: 2026-10-14
- **commit**: `570790e8` (HEAD)
- **分類**: `ok`
- **backend**: `codex` — 選んだ理由: 既定 (auto: codex を試し、unavailable のときだけ claude)
- **要求したモデル**: `既定 (~/.codex/config.toml)`
- **実効のモデル**: `gpt-6-astra` (codex のログの見出し)
- **CLI の版**: `codex-cli 0.162.1`
- **依頼側との関係**: 不明 (`--requester-model` 未指定)
- **permission_denials**: 対象外 (codex は read-only サンドボックス)
- **実行**: effort `high`, 3.8 min, rc=0
- **試行のログ**: `notes/reviews/2026-10-14-wall-orthogonality-judgments-diagnose.codex.log` (codex: ok)
- **結論**: 関門への接続前に、§4.4 の判定単位を `physID` ごとの帯から、境界の分割に依存しない列ごとの全層へ改訂する。

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（検証対象 `HEAD 570790e8`。ファイル変更・`forge` 起動なし。再現はメモリ上で既存の判定関数と `main()` を実行。）

| 判断 | 採否 | 根拠と対案 |
|---|---|---|
| 1. 格納座標で判定 | **採用** | 丸めを差し引いて合格にしない実装は妥当（[check_mesh_quality.py:665](/home/sano/work/forge-review-h370/solver_density_cuda/tools/check_mesh_quality.py:665)）。既存の合成格子を再測定すると、float32 は φ 最大 14.80°・AR 最大 5458.74、float64 は φ 最大 0.00000121°・AR 最大 4666.67。格納精度で形が変わっている。ただし、この検査だけでソルバ内部の幾何精度まで保証したことにはならない。 |
| 2. `physID` ごとの帯全体で緩和 | **却下 — Major** | [check_mesh_quality.py:674](/home/sano/work/forge-review-h370/solver_density_cuda/tools/check_mesh_quality.py:674) は壁ごとの最大深さまで一括判定する。**同じ座標・接続・壁種でも、壁の分割だけで `FAIL → WARN`、終了コード `1 → 0` を再現した**。対案は、壁面から対象セルまでの**列ごとの連続した全層**を検査する方式。`physID` は報告単位に使い、幾何の合否を決める単位にしない。 |
| 3. 複数 `physID` の接点を角部として除外 | **却下 — Major** | [check_mesh_quality.py:543](/home/sano/work/forge-review-h370/solver_density_cuda/tools/check_mesh_quality.py:543) は同一直線上の境界分割も角部とし、[同:728](/home/sano/work/forge-review-h370/solver_density_cuda/tools/check_mesh_quality.py:728) で警告から除く。AR 100・φ 40° の同じ格子で、1 個の `physID` なら警告 20 セル、20 個に分けると警告 0 セル、**`WARN → PASS`**。対案は、境界ラベルの継ぎ目と幾何学的な角を分離し、真の角部も最大角度・位置を報告すること。45° は暫定の分類値として保持できるが、警告免除の根拠にはしない。 |
| 4. `--bcond` なしの旧判定を残す | **採用（互換用途に限定）** | [check_mesh_quality.py:925](/home/sano/work/forge-review-h370/solver_density_cuda/tools/check_mesh_quality.py:925) の注意付き旧動作は残してよい。ただし、その `PASS` を新関門の合格証拠として使わない。関門側で `--bcond` の欠落を拒否する負例試験を必須にする。 |
| φ・θ_tol とも暫定 5° | **採用（構成の暫定基準）** | [plan:171](/home/sano/work/forge-review-h370/plans/active/tooling-generator-defaults-and-launch-gate.md:171) の精度予算との分離を維持する。今回見つかった判定の変化は閾値変更では直らない。「0.07° の約70倍」は安全率の根拠にはならず、5°以下でも生成方式や物理解の精度は保証しない。 |
| float32 旧メッシュの移行 | **採用（元入力からの再変換）** | [plan:204](/home/sano/work/forge-review-h370/plans/active/tooling-generator-defaults-and-launch-gate.md:204) に従い、元座標と変換条件を確認し、幾何を再計算して再判定する。**旧 HDF5 の座標を float64 に型変換するだけでは不可**。元入力の精度が失われていれば再生成が必要。新しいハッシュ・来歴で判定し、旧キャッシュを引き継がない。 |

**追加指摘 — Major:** `--bcond` を渡しても AR 5000 の上限をすり抜ける。[check_mesh_quality.py:651](/home/sano/work/forge-review-h370/solver_density_cuda/tools/check_mesh_quality.py:651) は `AR ≤ ar_max` だけを緩和検査に入れ、[同:950](/home/sano/work/forge-review-h370/solver_density_cuda/tools/check_mesh_quality.py:950) は上限超過を `SOFT-PASS` にできる。2000 セル中 1 セルだけ AR 6000.3 の入力で、**`VERDICT: SOFT-PASS`、終了コード 0** を実測した。対案は、関門では AR > 5000、および AR > 1000 で緩和条件を確認できないセルを、件数によらず `FAIL` にすること。

結論: 関門への接続前に、§4.4 の判定単位を `physID` ごとの帯から、境界の分割に依存しない列ごとの全層へ改訂する。

第 1 仮説: 判定の主要な欠陥は角度閾値ではなく、境界ラベルを幾何学的な帯・角部の定義に使っていること。 確度: 高  
　根拠: 上記コードと再現結果。直線壁の5セル（AR 約4000、4000、400、400、401）では、下流の低ARセルの φ 40° により、単一 `physID` で直交した上流2セルまで緩和不可になる。壁を分割するだけでその2セルが緩和可となった。  
　反証条件: 座標・接続・壁の面集合・壁種を固定し、`physID` の分割だけを変えても、セルごとの緩和可否と警告対象が不変であること。

第 2・第 3 仮説: 追加しない。今回の次手は上記の再現だけで決められる。

判別 A/B: **変えるのは壁の `physID` 分割だけ**。上記5セルを A＝単一ID、B＝上流2面と下流3面で分割し、キャッシュなし・CFD 0 step で判定する。見る量は緩和確認セル数・警告セル数・VERDICT・終了コード。**差が出れば第1仮説を支持、一致すればこの反例では反証**。現実装では A＝確認0・`FAIL`、B＝確認2・`WARN` を確認済み。修正後の受入値は両方とも確認2・警告1・`WARN` と事前登録する。

やらない方がよいこと: φ の閾値を緩める、丸め補正で合格にする、`physID` を細分化して通す、既存試験の `ALL PASS` だけで関門へ接続すること。

呼び出し側の前提への異議: 「帯全体」は `physID` 全体と同義ではない。「2つの `physID` が会う」は幾何学的な角の証拠ではない。また、角部は件数だけという説明に対し、実装は角度・位置も出すが、警告と通常部の最大値から除外している。case/45 の生産格子の数値は本レビューでは独立確認しておらず、今回の結論の根拠に使っていない。

不足情報: 生産移行を承認する段階では、case/45 の元メッシュ・変換条件・再変換後の品質VERDICTが必要。今回の設計判断は妨げない。

証拠: 足りる
