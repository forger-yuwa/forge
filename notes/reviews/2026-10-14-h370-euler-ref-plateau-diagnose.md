# codex 諮問 (diagnose): h370-euler-ref-plateau

- **brief**: [`notes/reviews/briefs/2026-10-14-h370-euler-ref-plateau.md`](../../notes/reviews/briefs/2026-10-14-h370-euler-ref-plateau.md)
- **plan**: [`plans/active/tooling-nozzle-isothermal-wall-chain.md`](../../plans/active/tooling-nozzle-isothermal-wall-chain.md)
- **date**: 2026-10-14
- **commit**: `fdf3b728` (HEAD)
- **codex**: effort `high`, 3.9 min, rc=0
- **結論**: **案 B を採用し、`run_0951` を保留したまま、同じ `res_60000` から実効 CFL だけを変える Euler の診断 A/B を別登録してください。**
- **extra**: `plans/accepted/discretization-moc-axis-limit-and-corrector.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表（Critical 0 / Major 3 / Minor 1）

| 論点 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| 案 A：残差条件を「有界な横ばい＋STEADY」に置換 | **却下・Major 1** | `σ_ref` は保存した参照間の変動を測り、時間的に一定な偏りは検出しません。[h370_judge.py:559](/home/sano/work/forge-review-h370/case/45.isobutane_m6_d155/h370_judge.py:559)。残差の振幅上限・傾きの条件だけでは、この不足を補えません。**現行 FAIL を維持して原因を判別する案 B を採用**します。 |
| `res_60000` を固定して NS を先行投入 | **却下・Major 2** | 登録は Euler の残差 PASS を要求し、実装も不合格なら固定を書かずに終了します。[計画:448](/home/sano/work/forge-review-h370/plans/active/tooling-nozzle-isothermal-wall-chain.md:448)、[h370_judge.py:362](/home/sano/work/forge-review-h370/case/45.isobutane_m6_d155/h370_judge.py:362)。**`run_0951` は保留**。後段の `σ_ref,use ≤ 0.01 mm` は追加の必要条件であり、未達の残差条件の代替ではありません。 |
| 旧 `run_0174` も同じ状態なので比較上は問題ない | **却下・Major 3** | V5d は残差停滞だけでなく、腕 M の全温の時間幅 **0.113〜0.126 K > 0.1 K**、一部量の準定常未達で保留でした。[旧計画:252](/home/sano/work/forge-review-h370/plans/accepted/discretization-moc-axis-limit-and-corrector.md:252)。旧参照は**保存された旧入力との感度比較**には使えますが、適格な正解として扱えません。小差から両者の妥当性を、大差から入口形状だけの効果を結論しないでください。 |
| `drift 0.0 %・fluct 0.0 %` の扱い | **要再検証・Minor 1** | 表示は小数点以下1桁です。[check_quasisteady.py:291](/home/sano/work/forge-review-h370/solver_density_cuda/tools/check_quasisteady.py:291)。読み取りだけの合成試験で、実際の変動幅 **0.03994 %** が `fluct=0.0%` と表示されました。報告値の「表示の桁で」は正しく、これを変動ゼロと解釈せず、生の数値を残してください。 |

結論: **案 B を採用し、`run_0951` を保留したまま、同じ `res_60000` から実効 CFL だけを変える Euler の診断 A/B を別登録してください。**

第 1 仮説: 現設定の擬似時間反復が残差の振動を維持しており、実効 CFL の低減で残差床を下げられる。確度: **低・未確認**。リミッタの切り替えが真因とはまだ特定できません。

根拠: ブリーフが報告する `case/45.isobutane_m6_d155/run_0950_h370_euler_ref/` は、本段5,000 step以降の残差が横ばいで、判定区間 **mid → main** の VERDICT は **`NOT CONVERGED (stalled/plateau)`**、8量は **`STEADY`** です。これは「現設定で停滞する」という観測であり、「この二次離散化では3桁低下が構造的に不可能」という証拠ではありません。判定器はピーク／末尾中央値を使い、3桁低下後の横ばいは許可します。[check_convergence.py:182](/home/sano/work/forge-review-h370/solver_density_cuda/tools/check_convergence.py:182)。

反証条件: 下記の低 CFL 側でも残差条件に届かず、両腕が同程度の床で停滞すること。その場合、「CFL 低減だけで今回の資格未達を解消できる」という仮説を退けます。

第 2 仮説: CFL に鈍感な演算誤差または局所的な離散化・境界処理の問題が床を作っている。確度: **低・未確認**。FP64 という申告だけでは全演算経路を確認できませんが、既に FP64 とされる run に対して、単に「倍精度で再実行」と提案する根拠もありません。

判別 A/B: **変更する数値設定は `time.deltaT.cfl_pseudo` のみ、2 → 0.5。**

- 両腕とも `run_0950/res_60000.h5` から、別々の新規 run に同一メッシュの `restart_field.py` 経由で開始し、保存量のビット一致を確認する。格子、BC、バイナリ、リミッタと基準値、熱物性、内部反復数は固定する。
- **対照腕 A：2、低 CFL 腕 B：0.5、各40,000 step、保存1,000 step**。比較窓は追加区間の前半・後半各20,000 stepに固定し、自動延長しない。表示用 `cfl` は対応する値に揃える。定常計算で実際に効くキーは `cfl_pseudo` です。[solver-settings.md:36](/home/sano/work/forge-review-h370/procedures/solver-settings.md:36)。
- 全残差列の中央値・振幅・傾きと、同じ8量の `check_quasisteady.py` VERDICTを見る。収束判定は、実効設定の同一性を確認した **元の mid → main → 各継続区間**で `--segment --drop 3` を適用する。未合格の親を `--from-floor` の参照にはしない。
- **A が停滞、B が全残差 PASS** → 「二次 Euler では不可避の床」を退け、CFL 依存の反復問題を支持する。ただし、リミッタ原因の確定でも、元の `res_60000` の適格化でもありません。
- **A・B とも同程度の床で停滞し不合格** → 今回の CFL 低減だけで解消する仮説を退ける。丸め誤差が確定したとは扱わない。
- B がまだ減衰中、A が元の停滞を再現しない、非有限値が出る場合は判別不能。短い過渡同士の場の差を「一致」と報告しない。

やらない方がよいこと: **観測した約2倍の振幅や ±0.003 dec/10,000 step に合わせて合格閾値を作ること**。傾きが有意でないことも、許容内でゼロに近いことの証明にはなりません。また、一次化・リミッタ凍結は離散作用素を変えるため、今回の最初の対照には選びません。`soft` を判定区間に混ぜて低下桁数を稼ぐこと、固定ファイルを手作業で作って NS の関門を迂回することも避けてください。

呼び出し側の前提への異議: **「残差が構造的に3桁下がらない」「設計入力として必要なのは `σ_ref,use` だけ」は未検証の解釈です。** 一方、3桁という相対低下条件自体も設計精度の保証ではありません。将来条件を改訂するなら、残差水準と抽出誤差の対応を独立に裏付け、今回の FAIL を残した別登録で検証する必要があります。

不足情報: この作業ツリーには `run_0950`・`run_0174` の実データがありません。上記 run の数値と VERDICT は**ブリーフの報告であり、独立再判定はしていません**。全段の残差 CSV、`stage_manifest.json`、実効 config、予備資格 JSON、丸め前の量の時系列が必要です。また、#36i の修正版はこのツリーでは確認できず、現版には代表4出力の選択と `U_time` の資格許可が残っています。

ファイル変更・forge 起動なし。**plan 未反映**。呼び出し側の反映先は `plans/active/tooling-nozzle-isothermal-wall-chain.md` §5.1 #36f・#36i、§6 の診断追加登録です。
