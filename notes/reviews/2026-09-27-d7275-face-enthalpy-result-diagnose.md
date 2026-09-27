# codex 諮問 (diagnose): d7275-face-enthalpy-result

- **brief**: [`notes/reviews/briefs/2026-09-27-d7275-face-enthalpy-result.md`](../../notes/reviews/briefs/2026-09-27-d7275-face-enthalpy-result.md)
- **plan**: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md)
- **date**: 2026-09-27
- **commit**: `3a17e638` (feature/gap-heating-precision)
- **codex**: effort `high`, 3.8 min, rc=0
- **結論**: **(c) を推奨する。T4-0b-HF は「面エンタルピー精度への大きな感度を観測、事前登録上は判別保留」と記録し、G15 未達の探索結果として他項目へ進む提案をユーザに出す。**
- **extra**: `case/60.flatplate_d7275_m7/README.md`, `case/60.flatplate_d7275_m7/acceptance.json`, `case/60.flatplate_d7275_m7/tools/hface_double.patch`, `plans/accepted/performance-3d-node-sst-speedup.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

結論: **(c) を推奨する。T4-0b-HF は「面エンタルピー精度への大きな感度を観測、事前登録上は判別保留」と記録し、G15 未達の探索結果として他項目へ進む提案をユーザに出す。**

採否表:

| 判断 | 採否 | 根拠・対案 |
|---|---|---|
| 自由流域の床の主因を「支持」と確定する | **要再検証** | **Major:** 延長後の窓間変化 −13.5〜−30% は保留条件に該当する。下降方向への例外はない。加えて、両腕ともプローブの分解能以下なので「局所振幅も減る」を確認できていない。**大幅な残差低下という観測と、登録上の判定を分ける。** [acceptance.json:238](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:238) |
| (a) opt-in 化と一様初期値からの再計算へ進む | **現時点では却下** | 保留時は追加の手当てを連鎖させない規則がある。再初期化で低下桁数が増えても、残差床の原因を解決した証拠にはならない。診断パッチは保存し、製品化は別の判断にする。[acceptance.json:245](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:245) |
| (b) ω の床を続けて調べる | **現時点では却下** | 未解決課題として残すが、同じ打ち切り規則に従う。 |
| (c) G15 未達のまま他項目へ進む提案 | **採用** | 登録済みの `stop_rule` に直接対応する。G15 の免除・合格への変更は伴わない。 |
| ω を既知の数値床として G15 から免除する | **却下** | **Major:** 大きな項の相殺と残差の停滞だけでは、許容できる丸め床や熱流束の反復誤差上限を証明できない。全残差の判定を維持する。[G14・G15](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:1030) |

第 1 仮説: **SLAU 面エンタルピーの float 評価が、今回の自由流域の残差水準を支配する主要因である。** 確度: **中（有力だが、実測原本は未検証）**。

根拠:

- 現行コードと指定コミット `4687c3c` の双方で、double 入力の `Y`・`T` を float に落として `thermo_h_mix_f` を評価することを確認した。[thermo_d.cuh:959](/home/sano/work/forge/solver_density_cuda/cuda_forge/thermo_d.cuh:959)
- 提示パッチは、同じ面状態に対するエンタルピー評価を double 係数・double 演算へ変更している。切替箇所は仮説に対応している。[hface_double.patch:1](/home/sano/work/forge/case/60.flatplate_d7275_m7/tools/hface_double.patch:1)
- 記録上、A の `case/60.flatplate_d7275_m7/run_0019_t26_hf_f/` → `run_0021_t26_hf_f_ext/` と、B の `case/60.flatplate_d7275_m7/run_0020_t26_hf_d/` → `run_0022_t26_hf_d_ext/` では、自由流域の流れ・化学種残差 B/A が約 **3.4e−4〜1.3e−3**。5k・10k の両時点で同じ桁の差がある。[acceptance.json:248](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:248)

これは強い感度の証拠であり、保留だから無効になるわけではない。ただし、**定常的な床の帰属・局所振動の減衰・G15 達成は別の主張**である。

反証条件: 実バイナリの切替と入力同一性を確認した比較で、過渡収束後の対象残差・分解可能な局所振幅が A の 0.8〜1.2 倍へ戻るなら、主因説を棄却する。

第 2 仮説: **再開操作が過渡を再励起し、延長区間の下降傾向に寄与している。** 確度: **中**。記録された B の `rms_ro` は再開前 6.8e−12 → 最大 5.5e−9 → 末尾 6.1e−12。ただし、失われた状態が乱流粘性なのか、別の初期化処理なのかは未確認。保存量のビット一致だけでは完全な再開忠実性を保証しない。[acceptance.json:251](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:251)

第 3 仮説: 追加しない。第一内部列の ω が残る原因は未特定として保持する。

判別 A/B: **今回指定された `stop_rule` を優先し、新規 A/B は提案しない。** 判別対象は実施済みの面エンタルピー float／double 一組、各 5k＋延長 5k とする。

- B/A ≤ 0.1、局所振幅減少、窓間変化 ≤5%が揃う → 支持。
- 過渡収束後も B/A = 0.8〜1.2 → 棄却。
- 今回 → **残差比の条件は満たすが、窓間変化が未達、振幅減少は測定不能なので保留**。

再開の跳ねを理由に登録済み条件を変更するなら、それは事後の規則変更であり、今回の「支持確定」には使えない。

やらない方がよいこと:

- **低下桁数を稼ぐ目的の一様初期化、異なる方程式・離散化の履歴連結、未収束 run の参照床化。** `--from-floor` も参照 run 自身の通常判定 PASS を要求する。[check_convergence.py:353](/home/sano/work/forge/solver_density_cuda/tools/check_convergence.py:353)
- 熱流束の小さな A/B 差を、反復誤差上限や ω の免除根拠にすること。
- 今回の結果だけで既定 float を変更すること。

呼び出し側の前提への異議:

- **Major — 「第一内部列の ω は面エンタルピーでは変わらない」は強すぎる。** 記録された B/A = 0.58 は **42%低下**である。対案は「寄与はあるが、1/10 以下の改善には届かず、残差停滞も解消していない」。[acceptance.json:253](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:253)
- **相殺残り ≈80 を、そのまま丸め床と認定できない。** `omg_trans` と各ソースは体積込みの収支量であり、実装上の意味は確認できた。しかし 80/(6e8) ≈ 1.3e−7 は局所的な釣り合いの尺度であって、G15 の許容や熱流束誤差への換算ではない。[ransSource_d.cu:238](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransSource_d.cu:238)
- **Minor — 「各窓11枚」は集計コードと異なる。** 標準指定では前窓 `[3000,4000)` は10枚、後窓 `[4000,5000]` は11枚になる。記録を実際の標本数に訂正する。[hf_ab.py:58](/home/sano/work/forge/case/60.flatplate_d7275_m7/tools/hf_ab.py:58)

ω の G15 での扱いは、現在の判定をそのまま残す。記録された延長区間の判定は以下である。

```text
run_0021_t26_hf_f_ext: NOT CONVERGED（全列横ばい）
run_0022_t26_hf_d_ext: NOT CONVERGED（流れ・化学種 falling、roOmega stalled/plateau）
```

出典は [acceptance.json:254](/home/sano/work/forge/case/60.flatplate_d7275_m7/acceptance.json:254)。**ω の絶対残差値だけから合否や将来の到達不能を断定せず、G15 未達を維持する。**

不足情報: ローカルの `case/60.flatplate_d7275_m7/` には **run ディレクトリ・残差履歴・VERDICT 原本がいずれも存在しない**。したがって数値は記録に基づく評価であり、NaN/Inf、実効設定、バイナリ同一性、領域別時系列、準定常 VERDICT は独立再検証できていない。`hf_ab.py` の連続プローブ集計も T・P のみで、要求された高精度 ρ・P・k・ω の振幅確認を満たさない。[hf_ab.py:81](/home/sano/work/forge/case/60.flatplate_d7275_m7/tools/hf_ab.py:81)

ファイル変更・`forge` 起動なし。**plan 未反映**。呼び出し側の反映先は、関連 plan §4.12・§5.1 #61・§6 G15 と [case README の run 一覧](/home/sano/work/forge/case/60.flatplate_d7275_m7/README.md:24)。
