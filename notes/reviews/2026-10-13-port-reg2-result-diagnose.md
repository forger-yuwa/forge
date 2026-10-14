# codex 諮問 (diagnose): port-reg2-result

- **brief**: [`notes/reviews/briefs/2026-10-13-port-reg2-result.md`](../../notes/reviews/briefs/2026-10-13-port-reg2-result.md)
- **plan**: [`plans/active/architecture-float-state-double-geometry.md`](../../plans/active/architecture-float-state-double-geometry.md)
- **date**: 2026-10-13
- **commit**: `c967a19f` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 4.0 min, rc=0
- **結論**: **main への統合を保留し、float の移植前後を固定加算順で比較する短い診断A/Bを、別試験として事前登録する。**
- **extra**: `case/45.isobutane_m6_d155/port_reg2_an.py`, `case/45.isobutane_m6_d155/_band_ab/cold_pair/port_reg2_judge.json`, `notes/investigations/2026-10-13-main-integration/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 重大度 | 判断 | 根拠と対案 |
|---|---|---|
| Major | **main への統合は現時点では却下** | 統合条件は R1・R3 合格だが、§6.60 は総合 `判別不能`。[統合記録:26](/home/sano/work/forge-integ-1005/notes/investigations/2026-10-13-main-integration/README.md:26)。結果Bが無いことを合格に読み替えず、float の未判定を解消する別試験を事前登録する。 |
| — | **FP64 の限定的な結果Aは採用** | 保存済みJSONを再集計し、最大 `(D+E_A+E_B)/τ` は point **0.07241**、B0 **0.21977**、平板 **0.0001207**。書けるのは「指定入力・10 step・登録指標で結果A」まで。収束解の一致、float の非回帰、統合全体の合格には拡張しない。[plan:1544](/home/sano/work/forge-integ-1005/plans/active/architecture-float-state-double-geometry.md:1544) |
| Major | **同じ判定式のまま本数を増やす案は却下** | `E` は平均の標準誤差ではなく標本の最大幅。[判定コード:48](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/port_reg2_an.py:48)。B0 float の `vis_turb` は **E_A/τ＝10.68、E_B/τ＝12.33**。τを固定すれば、Dがゼロでも合格しない。既存標本への追加で幅は縮まらないため、先に再現性を制御する診断へ変える。 |
| Major | **確認4は未検証として残す扱いを採用。ただし合格の代用にはしない** | 旧版・新版とも入力読込みで停止しており、その先の準備処理を比較できていない。[統合記録:32](/home/sano/work/forge-integ-1005/notes/investigations/2026-10-13-main-integration/README.md:32)。参照入力を復旧して比較する。代替入力を使う場合は別試験とし、生産入力の再現確認が未了であることを維持する。 |
| Minor | **ブリーフの「B0 float の `roUy` は許容内」は却下** | `case/45.isobutane_m6_d155/run_0792〜0797_q1_b0_f32_*` の集計は、D＝0.07309、E_A＝0.14325、E_B＝0.13040、τ＝0.16453。**登録式では約2.108τ**。D単独がτ未満という記述に訂正する。 |
| Minor | **判定器の入力検査は要修正** | 全runで空の量集合を返す模擬入力でも `A_ok=True`、全runが同じ1行の履歴でも合格した。[判定コード:39](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/port_reg2_an.py:39)。必須量・空集合・予定step列を検査する。今回のJSONには量があるため、これを今回の判別不能の原因とはしない。 |

結論: **main への統合を保留し、float の移植前後を固定加算順で比較する短い診断A/Bを、別試験として事前登録する。**

第 1 仮説: float の判別不能は、既存の実行間ばらつきが許容を消費していることが主因である。移植による差の不存在は未証明。 確度: **中**

　根拠: B0 float の `vis_turb` は D＝1.1723e−5、E_A＋E_B＝5.0355e−5、τ＝2.1886e−6。群内ばらつきだけで約23τある。ただし、Dがその幅より小さいことは同等性の証明ではない。

　反証条件: 同一入力・固定加算順で各腕が再現するにもかかわらず、移植前後に差が残ること。

第 2 仮説: 移植による実効入力・ホスト側処理の差がばらつきに隠れている。**未確認**。SASS の同一性だけでは、渡す配列や呼出し条件まで除外できない。

判別 A/B: **移植前／後だけを変える**。R1 B0 float の共通入力で、両版に同じ固定順集計の診断処理を適用し、まず **1 step、各2回**。初期化後の状態、残差、実効dt、更新後の全保存量を比較する。各腕内のビット再現と有限性を前提とする。

- **A：腕内・腕間ともビット同一** → この1 stepについて、第2仮説の決定的な差を棄却する。
- **B：腕内は再現するが腕間で差がある** → 「実行間ばらつきだけ」という説明を棄却し、最初に差が出る処理を調べる。
- 腕内が再現しなければ診断は不成立。Aでも§6.60を合格に書き換えない。統合ゲートへの採用範囲は実行前に明記する。

やらない方がよいこと: 観測された最大比28倍に合わせたτの拡大、合格するまでの追加反復、`D < E_A+E_B` による合格扱い、確認4の双方停止を再現成功とすること。

呼び出し側の前提への異議: 「float のばらつきが大きい」は観測として支持されるが、「移植の影響が無い」は未確定。R2の結果Aも、β＝1で確認した境界修復の範囲に限る。`speciesPrecondDt` のβ＜1・多成分検証は別に未了と記載されている。[設定仕様:168](/home/sano/work/forge-integ-1005/procedures/solver-settings.md:168)

不足情報: 実runのHDF5・CSV・設定・実行証跡はこの作業ツリーに無く、今回はJSONの再集計と判定コードの確認まで。固定加算順の診断機能の有無も未確認で、無ければ追加実装が必要。統合確認1〜3・5の実行結果も提示記録にはない。**ファイル変更なし・plan未反映**。採用する判断は呼び出し側が§5.1 #26と§6.60後続節へ反映すること。
