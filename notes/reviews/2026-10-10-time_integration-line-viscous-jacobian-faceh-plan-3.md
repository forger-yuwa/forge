# codex レビュー: time_integration-line-viscous-jacobian-faceh (plan)

- **plan**: [`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `b196885b` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 2.9 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m3
- **focus**: §6.5 (2026-10-10 に追加した、参照の差分の演算精度だけを変える A/B の事前登録) と、それを実装した test_line_visc_jacobian.cpp の差分 (作業ツリーの未 commit の変更、git diff で見る。引数 --ref/--eta/--dump、mp_flux・fd_jac_mp・dir_fd_mp・ulp の計算・apply_eta・全列の書き出し) だけを点検する。§6.3・§6.4 は実施・記録済みで対象外。多倍長のモードはまだ回していない。引数なしの出力は登録の試験の出力 (notes/reviews/briefs/2026-10-10-uj-colwise-output.txt) と完全一致を確かめた。点検してほしいこと: (1) A と B で変わるのが参照の差分の精度だけになっているか (標本・J・差分幅・拘束の方向)。(2) η の定義と閾値 (1e-8、零列 1e-12) と、分類への組み込み方。(3) mp100 の実装の正しさ (h の取り方、double への変換の位置、ulp)。(4) 判定の分岐と、PASS/FAIL で言えること。
- **extra**: `solver_density_cuda/tools/test_line_visc_jacobian.cpp`, `notes/reviews/2026-10-10-uj-colwise-result-diagnose.md`, `notes/reviews/briefs/2026-10-10-uj-colwise-output.txt`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
多倍長参照による切り分けは妥当ですが、入力・差分幅の固定、非有限の扱い、PASS/FAIL の解釈を実行前に修正してください。  
対象は §6.5 と未 commit 差分のみです。Critical はありません。

1. **Major — A/B で保存変数と差分幅の数値が固定されていない。**

   **根拠:** A は `to_q` で double の積を作り、double で h を計算します。一方 B は `mp_q` で積から多倍長計算し、h も再計算しています。[`test_line_visc_jacobian.cpp:41`](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:41)、[同:93](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:93)、[同:110](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:110)。

   人工例 `ρ=0.1, u=1700, hs=1e-6` では、A の運動量は `170`、B は `170.000000000000009436…`。h の相対差も約 `9.79e-17` になります。小さい差であり、旧不一致13列の原因だとは主張しませんが、§6.5 の厳密な固定条件には合いません。

   **対案:** double の基準 Q と各幅 h を共通で生成し、B はその値を MP に持ち上げて使ってください。温度拘束の `dq` は現在どおり共通値を使用します。Q±h の演算精度を変えることは、参照計算の精度変更に含むと明記します。

2. **Major — η の非有限が消え、非有限を含む試験が PASS になり得る。**

   **根拠:** `apply_eta` は有限性を確認せず `std::max` を使います。[同:180](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:180)。同じ比較演算を人工入力で確認すると、`U1=U2=全NaN` でも η は `0` のままです。Inf は通常の「解像しない」に落ちます。

   また集計は `rand` の `C_NF_FD` を通常の未解像枠に含めています。[同:415](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:415)。他が合格なら、非有限1列でも `1 ≤ 2000/100` により PASS 可能です。これは [§6.5:188](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:188) の「非有限があれば PASS にしない」に反します。

   **対案:** η の入力・正規化後の値を明示的に有限検査し、新登録では非有限を1%枠から除外してください。J の非有限は FAIL、参照・η の非有限は PASS 禁止とします。NaN/Inf、閾値直前・同値・直後、零列を人工入力で確認してください。

3. **Major — B の総合判定だけでは、旧13列の原因を判定できない。**

   **根拠:** [§6.5:189–192 の分岐](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:189) は、A を判定条件にせず B の PASS から丸め仮説を支持します。しかし B は rand の未解像を20列まで許します。**旧不一致の rand 5列が全部未解像でも、他が合格なら総合 PASS は可能**です。

   逆に総合 FAIL には、精度切替と無関係な短いラインや float 比較も含まれます。[試験コード:420](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:420)。それらの失敗から「旧13列は丸めだけでは説明できない」とは言えません。

   **対案:** 総合 VERDICT と原因判定を分離してください。丸め仮説の支持条件を、**A のη適用前の旧不一致再現＋同じ13列すべての B での解像・合格＋B の総合 PASS**にします。旧13列に B の解像済み不一致が残れば丸めだけの説明を棄却し、別項目の FAIL はその項目の問題として報告します。

4. **Minor — `ulp_mp` は実際の ulp ではない。**

   **根拠:** [同:104](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:104) の `|x|·epsilon` は、二進浮動小数点の隣接間隔と異なります。例えば `x=1.5` では実際の間隔の1.5倍です。通常の非零値では保守側ですが、A の `nextafter` による定義と揃っていません。

   **対案:** MP の隣接値との差で ulp を定義し、零・非有限を明示処理してください。η の比まで MP で計算して最後に double 化すると、途中変換による過小評価も避けられます。なお、**FD 本体を差・除算の後で double 化している位置は適切**です。

5. **Minor — 引数の誤りが検出されず、指定した試験と違うものを実行できる。**

   **根拠:** [同:192](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:192) は `--ref mp100 --ref double` でも MP のままです。`--eta` は `atof` のため、不正文字列が0になりゲートを無効化します。NaN・Inf・負値も拒否しません。

   **対案:** 重複指定を拒否し、`--eta` は文字列全体の変換成功と有限・正値を検査してください。今回の登録コマンド自体は正しいため、これは実行時の取り違え防止です。

6. **Minor — 全列 CSV は追加されたが、閾値近傍の再判定と保存成功を保証しない。**

   **根拠:** [同:443](/home/sano/work/forge-faceh/solver_density_cuda/tools/test_line_visc_jacobian.cpp:443) は数値を `%.6e` で保存します。例えば `1.00000004e-6` は `1.000000e-6` となり、閾値超過が記録から消えます。`fprintf`・`fclose` の失敗も確認していません。

   **対案:** double は `%.17g` 等で保存し、η適用前後の分類を併記してください。書き込み・close の失敗は非成功終了にします。全3560列を列挙するループ自体は適切です。

**推奨は、上記を修正したうえで、この host A/B を実施することです。** 優先順は①入力・h の共通化、②非有限ゲート、③旧13列に結び付いた原因判定、④補助実装の修正です。η の `1e-8`／零列 `1e-12` は各誤差許容の1/100であり、量子化を理由に判定を保留する閾値として維持して構いません。内部演算誤差の保証には使えません。

確認できた範囲では、乱数生成・標本生成・拘束方向は旧コードと同一で、`mp_flux` の式も元の流束と整合しています。構文検査 `g++ -fsyntax-only` は成功しました。B 本試験、引数なし出力の再実行、forge・収束／準定常判定は実施していません。ファイル変更なし、**plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 3
