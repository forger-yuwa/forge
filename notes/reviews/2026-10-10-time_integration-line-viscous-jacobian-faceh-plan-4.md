# codex レビュー: time_integration-line-viscous-jacobian-faceh (plan)

- **plan**: [`plans/active/time_integration-line-viscous-jacobian-faceh.md`](../../plans/active/time_integration-line-viscous-jacobian-faceh.md)
- **stage**: `plan`
- **date**: 2026-10-10
- **commit**: `e0ab3edb` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 6.0 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2
- **focus**: §6.7 (2026-10-10 に追加した、製品の経路の照合の事前登録) と、その実装だけを点検する: solver_density_cuda/cuda_forge/timeIntegration_d.cu の #if defined(FORGE_LINE_AUDIT) の塊 (作業ツリーの未 commit の変更、git diff で見る。通常のビルドでは塊を取り除いたファイルが変更前と完全一致することを確かめた)、solver_density_cuda/tools/line_audit_helper.cpp、case/45.isobutane_m6_d155/lvcaudit_judge.py。まだビルド・採取していない。点検してほしいこと: (1) 記録の位置と中身で、係数の生成・float の演算・K の置き場所・拘束の行・D の組立を照合できるか (抜けている記録、間違った番地、重なり)。(2) 監査の書き込みが通常の数値を変えないか (監査用のビルドの中で)。(3) 判定の段 G・A・B・C・D の式と許容 (1e-6、完全一致、A の 1e-5・1e-3) が妥当か、甘さ・抜けはないか。(4) A (元の入力から double で計算した係数との差) の分類と、PASS・FAIL・A の結論で言えることの範囲。特に、ISP 0 で dcc を float の座標の差で作る点 (timeIntegration_d.cu の粘性の幾何) を、A で捉えられるか。
- **extra**: `solver_density_cuda/tools/line_audit_helper.cpp`, `case/45.isobutane_m6_d155/lvcaudit_judge.py`, `notes/reviews/2026-10-10-uj-mp-result-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
製品の係数生成から組立までを照合する目的と順序は妥当です。  
ただし現状の判定器には、未検証の入力を PASS にする穴と、正常な float 演算を FAIL にする穴があります。採取前に修正してください。

1. **Major — G に採取条件・網羅性・有限性のゲートがなく、空の監査が PASS します。**

   根拠: [`lvcaudit_judge.py:59`](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:59) は記録された節点数をそのまま採用し、[同:250](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:250) は検査件数 0 を拒否しません。I/O をメモリ上に置き換えて実際の `main()` を実行したところ、次の人工入力がともに終了コード 0 でした。

   - 節点 0・面 0 → `PASS`、A も「係数は 1e-3 以内」。
   - 節点 1・面 0、`Kprev/Knext`・状態・物性・`dt_vol` が NaN → 同じく `PASS`。

   また、実効 ISP・マスク・キー・factor 回数・FP64 の検査がありません。現在の [`flowFormat.hpp:6`](/home/sano/work/forge-faceh/solver_density_cuda/flowFormat.hpp:6) は `float` です。誤って通常の float ビルドで採取すると、A の「生座標」も既に丸められ、今回検出したい差を失います。

   **対案:** 最初に独立した採取ゲートを設け、入力・メッシュ・バイナリのハッシュ、型の幅、実効設定、指定した 5 ラインの全節点・両向きの面、sweep 0〜4、必須配列の形状・有限性・正値条件を検査してください。欠落・重複・未知の枝・読めない記録は理由付き `INVALID` とし、G 不合格なら helper を呼ばず停止します。

2. **Major — 接続・拘束・状態の「完全一致」が、独立した入力との一致になっていません。**

   根拠: [`lvcaudit_judge.py:85`](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:85) はライン面の個数と隣の集合への所属だけを確認します。`lp` 側を二重に記録して `ln` 側を欠落させる場合、端点の未使用 K、逆向き面との対応を保証しません。

   [`line_audit_helper.cpp:59`](/home/sano/work/forge-faceh/solver_density_cuda/tools/line_audit_helper.cpp:59) は、記録された隣の状態・`cp`・`jVel/jTemp`・マスクをそのまま採用します。これらを隣節点の原入力・境界フラグと照合していません。誤った拘束フラグを CUDA と helper が共有すれば、その誤りを再現して合格します。`dt_vol` も読み込むだけで、時間項の入力との対応を検査しません。

   **対案:** `lp/ln` を順序付きで確定し、各向きに面がちょうど一つあること、面番号・両端・逆向きの対応を検査してください。さらに、共通関数へ渡した全引数を原入力の cast・floor と照合し、特に `jVel/jTemp` を**相手節点の境界フラグから独立に**求めます。存在しない隣の K も検査対象に含めてください。

3. **Major — B はライン面の薄層係数しか検査せず、D を構成する他の係数生成が未監査です。**

   根拠: [`lvcaudit_judge.py:131`](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:131) と 134 行の `continue` により、A・B は `isLine && br == 1` だけです。ライン外面の `k_face/cfac`、スカラー粘性対角、軸対称の `hoop` は記録値を信頼して組み直しています。B は面積・法線も照合しません。

   したがって、例えば間違った法線を CUDA と helper が共有しても、B・C はその誤りを検出できません。A は法線差を記録しますが、分類には使いません。

   **対案:** 選択節点に接続する全対象面へ係数検査を広げ、必要なライン外隣接節点の物性も採取してください。`k_face/cfac`、`ν_eff/viscous_diag`、軸対称の分岐・`A_pl/r_eff/hoop` を原入力から再計算します。スカラー側の層流粘性は現実装では `cfg.visc` を使うため、節点の `vis_lam` と取り違えず、その実引数も記録する必要があります。根拠は [`timeIntegration_d.cu:852`](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:852) と [同:1162](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1162) です。

4. **Major — C・D は製品の float 加算を再現しておらず、固定の 1e-6 では誤判定を防げません。**

   根拠: [`timeIntegration_d.cu:1141`](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1141) は `D_after − D_before` を **ST で引いてから** double に変換します。[`lvcaudit_judge.py:186`](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:186) はそれを double で足し戻します。

   float32 の人工反例を実行すると、`D_before = 100000000`、正しい逐次加算後の `D_after = 1` に対し、保存増分は `−100000000`、判定器の復元値は `0` になりました。**正しい組立でも相対誤差 1 で FAIL** します。

   また製品は「対流 → 値 3 のスカラー → 薄層」の順ですが、判定器は「対流 → 薄層 → スカラー」です。[製品:1029](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1029)、[判定器:118](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:118)。監査の薄層 D は実際の累積 D の寄与を採取したものではなく、ゼロ行列へ別途再計算したものです。[製品:1070](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/timeIntegration_d.cu:1070)。

   さらにローカルの `nvcc --help` で FMA contraction の既定は有効と確認しました。NumPy の個別演算や指定の host ビルドと、演算の丸め方が同じとは限りません。

   **対案:** 実際の D の各段階を保存し、helper で同じ初期値・面順・加算順を再現してください。増分を残すなら最低限 `double(after) − double(before)` とします。許容は項別の尺度と丸め誤差の上限を登録し、相殺や FMA 差で説明できる未解像の比較は「判別不能」にします。大きな総 D を分母にした比較だけで、小さい寄与まで検証済みにしないでください。

5. **Major — 監査ビルドが通常経路を変えていないことが、採用条件になっていません。**

   根拠: [plan:241](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:241) は別ソース版の `lineM_fp64` との比較を「記録だけ」とし、[`lvcaudit_judge.py:229`](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:229) は比較配列が欠けても黙って省略します。不一致でも総合 PASS は変わりません。

   監査追加による製品配列への直接の書き込みは見つかりませんでした。ただし追加の一時配列・再計算・store はコード生成を変え得るため、ソース上の無変更だけでは監査ビルド内の数値不変を証明できません。

   **対案:** 最初から**同一ソース・同一コンパイル条件で `FORGE_LINE_AUDIT` の有無だけが違う二本**を対照にしてください。必須比較配列と節点順を固定し、不一致・欠落は通常経路への結論を保留する条件にします。ビット一致は `np.array_equal` ではなく整数 view 等で確認してください。実際、現在の比較は `+0.0` と `−0.0` を一致とします。

6. **Minor — A は座標差の劣化を捉えられますが、分類名と結論が測定範囲を超えています。**

   A の基本式は、floor が作動しない条件では妥当です。δ/dcc = S²/|Δx·S| なので、**float 化した絶対座標から差を取ることで生じる、投影距離の誤差**を β・κ の差として検出できます。[`lvcaudit_judge.py:143`](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:143)。

   人工例 `y₀=0.01`、`y₁=0.010000001` では、double の距離は約 `1.0e−9`、float の距離は `9.3132257e−10`、同じ物性で β の相対差は **7.374%** でした。これは現在の A で「許容外」になります。ただし実 run の測定値ではありません。

   問題は、[同:243](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:243) が β・κ だけで分類しながら、[同:257](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcaudit_judge.py:257) が係数全般について「1e-3 以内」と述べる点です。また `≤1e-5` は丸めが原因である証明ではなく、`>1e-3` でも原因が座標の丸めという場合があります。

   **対案:** 分類名を「差 ≤1e-5／≤1e-3／>1e-3」とし、結論を対象面の β・κ に限定してください。dcc・法線等は別に判定します。**再現 PASS と A の許容外は両立する**ため、その場合は「製品演算は再現できたが、double 基準との差は残る」と報告し、発散原因の断定や係数生成の問題の除外には使わないでください。

7. **Minor — §2・§7 と現在の監査スコープが矛盾しています。**

   根拠: [plan:28](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:28) と [同:261](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-faceh.md:261) は「コードは変えない」のままです。§6.7 は CUDA の監査コードと host helper を追加しています。

   **対案:** 監査ビルド限定の変更と未検証範囲を明記し、§5.1 を「記録・判定器の検証 → 同一版の監査有無比較 → 本採取・判定」の順へ更新してください。今回の case/45・1 step は局所組立監査として適切ですが、cell・周期・収束性への保証は含めません。

**推奨は、製品の数値設定を維持したまま、監査を修正してから同一版の監査有無一組を採取することです。** 採取前の修正優先順は **1 → 2 → 3 → 4 → 5**、続いて 6・7 の文言を同期してください。

確認済みなのは、監査 16 ブロックを除去したファイルと `HEAD` のバイト一致、記録枠内の番地に重なりがないこと、helper の C++ 構文検査、判定器の構文と上記人工入力です。CUDA ビルド・実採取は未実施で、数値不変や実 run の PASS は未確認です。ファイル変更なし、**plan 未反映**です。

指摘数: Critical 0 / Major 5 / Minor 2
