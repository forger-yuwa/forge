# codex レビュー: time_integration-line-viscous-jacobian-dt (plan)

- **plan**: [`plans/active/time_integration-line-viscous-jacobian-dt.md`](../../plans/active/time_integration-line-viscous-jacobian-dt.md)
- **stage**: `plan`
- **date**: 2026-10-11
- **commit**: `8a8ee7f8` (feature/faceh-audit-viscjac-close)
- **codex**: effort `high`, 4.4 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M5/m2
- **focus**: §6 (方向別 dt の上限の A/B の事前登録) と、台本 lvcdt.sh・事前のゲート lvcdt_pregate.py・本判定 lvcdt_judge.py (作業ツリーの未 commit の新規ファイル) を点検する。faceh plan §6.11〜§6.13 の lvc75.sh・lvc75_pregate.py・lvc75_judge.py (plan-6 レビューと 2 回の諮問を採用済み、実行して PASS・支持) を元にした。腕の違いは time.deltaT.lineDtDirectionalCap (なし/50) と、条件付きの C の lineDtDirectional (なし)。事前のゲートの V で、K の不変・Δτ の比 r・D の対角の変化が時間項 V/Δτ の差で説明できることを確かめる設計が十分か、分岐 (帳簿の第一層の Δρ/ρ の条件を含む) と言える範囲が妥当か、条件付きの C の起動の仕方、抜けている確認は何か。設計の諮問は diagnostician (notes/reviews/2026-10-11-lvcdt-design-diagnose.md)。
- **extra**: `case/45.isobutane_m6_d155/lvcdt.sh`, `case/45.isobutane_m6_d155/lvcdt_pregate.py`, `case/45.isobutane_m6_d155/lvcdt_judge.py`, `notes/reviews/2026-10-11-lvcdt-design-diagnose.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**
上限 50 を試す目的と単因子 A/B は妥当です。ただし、介入ゲートと帳簿判定に誤通過を再現しました。
以下を実装前に修正してください。現状の分岐では、dt 短縮の効果から H1〜H3 の機構まで判別できません。

目的・スコープは既存 plan と重複していません。親 plan の上限 50 は値 2、今回は値 3・マスク 7です。固定バイナリによる case/45 の診断として適切で、cell・周期・FP32 全般への保証にはなりません。

既存証拠を読み直し、`check_convergence.py` も変更なしで実行しました。元の配置はいずれも `case/45.isobutane_m6_d155/`、手元の証拠は `/home/sano/forge-evidence/2026-10-10-faceh/` にあります。

| run | 確認結果 |
|---|---|
| `run_0570_lvc75_m7_a1` | detectNaN: step 348、`DIVERGED (NaN/Inf)` |
| `run_0572_lvc75_m7_a2` | detectNaN: step 280、`DIVERGED (NaN/Inf)` |
| `run_0571_lvc75_m5_b1`・`run_0573_lvc75_m5_b2` | step 1999、両方 `NOT CONVERGED` |

全４本の保存済みメッシュ判定は `PASS (AR<=5000, skew<=0.90)`。索引は [case README:194](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/README.md:194) です。

1. **Major — 帳簿が欠損すると「支持」になる。**

   根拠: [lvcdt_judge.py:151](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcdt_judge.py:151) の `ledger_check()` は call と node の集合だけを確認します。[同:256](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcdt_judge.py:256) の `ledger_drho()` は `after_eos_bc` がなくても初期値 `mx=0.0` を返します。

   **再現:** 既存 B の帳簿から `after_eos_bc` 行だけをメモリ内で除去すると、`ledger_check()` は問題なし、`ledger_drho()` は **0.0**。既存 A の実測値 **0.544228／0.566672** と組み合わせると `decide()` は **「支持」**を返しました。ファイルは変更していません。

   対案: 判定対象の全 `(call, node, tag, field)` の存在・一意性・有限性を検査する。基準密度は正値を必須とし、欠損・重複・非有限は `INVALID`。A が200 call前に壊れた場合も、比較窓を黙って短縮せず事前登録してください。`inf` を比較の分母側に使って支持させてはいけません。

2. **Major — V は dt の変更を検出するが、「上限 50」の実効性を保証しない。**

   根拠: [lvcdt_pregate.py:289](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcdt_pregate.py:289) は縮小方向・第一層の `r≥10`・C/B の大小を確認するだけです。D の検査も、その場で読んだ dt に追従していれば通ります。

   **再現:** `run_0577` の実配列を基に、メモリ内で **B=C=A/50** とし、非拘束対角だけを時間項差で変更しました。現在の **V 全15項目が合格**しました。第一層の解像条件も許容の約 **4.9万〜8.8万倍**で通ります。これは指定した B/C の区別が失われても通る反例です。

   対案: 今回の設定では、同一状態の C を基準に全605節点で **Δτ_B ≈ min(Δτ_A, 50Δτ_C)** を、実装の演算精度に応じた許容で検査する。`axisTimestepBeta`・粘性 CFL 割引など、この関係の前提となる設定も固定してください。軸項を有効にする場合は [setDT_d.cu:194](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/setDT_d.cu:194) 以降の加算順序を含めて再構成する必要があります。

   K 不変、拘束行不変、D 差の時間項照合という方針自体は正しいです。**64 ulp＋解像条件は保持し、dt の用量検査を追加する**のが適切です。

3. **Major — 第一層の「1/10以下」は、同じ過渡を遅く進んだ可能性を排除しない。H1 の支持・棄却が強すぎる。**

   根拠: [plan:111](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-dt.md:111)、[lvcdt_judge.py:289](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcdt_judge.py:289)。

   初期変化が概ね Δρ ∝ ΣΔτ なら、dt を1/10以下にするだけでこの条件を満たせます。増幅機構が残っていても成立します。既存7/5の再計算では最大変化が **0.544228／0.566672 対 0.076076／0.081993**で、現判定式の比は **0.150659**。この例を排除できることは、別の遅延過渡を排除できる証明にはなりません。

   逆側も、B の非有限だけでは「熱拡散数≫1」を棄却できません。ゲートの `r≥10` では、引用値23を使っても2.3までしか下がらず、**実際の熱拡散数≤0.5を検査していません**。物性・状態も実行中に変わります。さらに dt は SST・化学種にも作用します。

   対案: 今回の主判定を **「上限50による2000 step内の非有限化回避／回避せず」**に限定する。Δρ/ρ は追加観測に下げ、H1 の機構的支持・棄却とは分離してください。C 有限→H3、C 発散→H2も、**次の調査順序**としてのみ扱う。H2とH3はともに時間項の強化に応答し得るため、Cだけでは分離できません。将来の行列調査では、単位依存の生の対角優位性だけでなく、スケーリングと実際の反復写像を評価すべきです。

4. **Major — 登録した累積擬似時間を、現在の台本では測れない。**

   根拠: [plan:125](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-dt.md:125) は累積擬似時間の併記を要求しますが、[lvcdt.sh:70](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcdt.sh:70) の腕の追加出力は `res_ro,volume` のみ。帳簿の状態変数一覧にも `dt_local` はありません（[convectiveFlux_d.cu:701](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/convection/convectiveFlux_d.cu:701)）。既存112節点帳簿でも不在を確認しました。

   dt は状態に応じて再計算されるため、**初回の dt×step は累積値ではありません**。また、介入ダンプには支持条件の節点4959・6169が含まれていません。

   対案: 判定対象節点をダンプと帳簿で揃え、更新に実際に使った `dt_local` を各反復で記録して ΣΔτ を作る。固定バイナリの制約で取得できないなら、事前登録から「実測累積擬似時間」を外し、未測定と明記してください。初回値による換算を代用してはいけません。

5. **Major — 条件付き C の起動が、現在の A/B 証拠に結び付いていない。**

   根拠: [lvcdt.sh:107](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcdt.sh:107) は保存済み JSON の `main_v3.verdict == "棄却"` だけを読みます。直前の `pregate --verify` が検査するのはダンプ側で、A/B の CSV・帳簿・場・判定器の変更は検知しません。

   古い判定後に A/B の証拠が欠けたり変わったりしても、C を開始できます。「本判定の現在の証拠が分岐2」という条件を満たしていません。

   対案: C 起動直前に **A/B限定の再判定**を行い、入力・証拠ゲート合格、４本とも `DIVERGED` を要求する。あるいは本判定に証拠と判定器のハッシュを保存し、起動時に検証する。途中まで作られた C の存在が A/B 再判定を妨げない入口も必要です。

6. **Minor — 登録した10倍到達 step を判定器が出していない。**

   根拠: [plan:121](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-dt.md:121) に対し、[lvcdt_judge.py:249](/home/sano/work/forge-faceh/case/45.isobutane_m6_d155/lvcdt_judge.py:249) は `first_step_over_3x` を計算します。

   対案: `first_step_over_10x` を追加し、既存 `run_0570/0572` の **roUy=4、roOmega=14、ro=16、roe=34**を再現することを確認する。第一・第二層の指定 call の量、非有限節点集合についても、出力担当と成果物名を§5.1に明記してください。

7. **Minor — ライン被覆100%から「全節点のdtが変わる」は導けない。**

   根拠: [plan:37](/home/sano/work/forge-faceh/plans/active/time_integration-line-viscous-jacobian-dt.md:37)。実装は内部ライン面を max から除くだけです（[setDT_d.cu:167](/home/sano/work/forge-faceh/solver_density_cuda/cuda_forge/setDT_d.cu:167)）。残る面が元々最大なら dt は変わりません。壁境界半割面も残ります。

   対案: 「全節点が方向別処理の対象。実際にdtが変わる節点は実測する」と訂正する。605節点の検査結果を全570999節点の実測として扱わないでください。

**推奨は、A/B/Cの構成を維持したまま、判定を短期の非有限化回避に限定して進めることです。** 実装前の優先順は **①帳簿の欠損拒否 → ②上限50の定量照合 → ③機構に関する判定の限定と記録可能量の整合 → ④C起動時の再検証 → ⑤記録項目の修正**。新しい数値スキームや広範な回帰試験を先に増やす必要はありません。

ファイルは変更していません。提案は **plan未反映**です。

指摘数: Critical 0 / Major 5 / Minor 2
