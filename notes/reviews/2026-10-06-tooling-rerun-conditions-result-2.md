# codex レビュー: tooling-rerun-conditions (result)

- **plan**: [`plans/active/tooling-rerun-conditions.md`](../../plans/active/tooling-rerun-conditions.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-10-06
- **commit**: `3c613995` (feature/nozzle-wall-fit-and-pipeline)
- **codex**: effort `high`, 5.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M3/m2

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

主要 run の監査可能化と訂正後の準定常判定は確認できました。ただし、既存 Euler CLI の実行不能と、起動推奨の誤分類・適用範囲の誤りが残っています。現状での `accepted` 移動には反対します。

1. **Major — 既存 Euler CLI が `NameError` で停止する**

   根拠: [runner_axismach.py:748](/home/sano/work/forge-integ-1005/design/forge_design/evaluate/runner_axismach.py:748) は新設の `_first_order` を呼びますが、CLI は同ファイル837行で `main()` に入り、`_first_order` の定義は1092行です。既定の `--stages full` と `soft` では定義に到達する前に呼び出します。

   CLI と同じ定義順序をメモリ上で再現し、ファイルを書かずに次を確認しました。

   ```text
   line 748, in run_staged
   NameError: name '_first_order' is not defined
   ```

   **対案:** `if __name__ == "__main__"` を全関数定義の後へ移動し、モジュールの import 経由だけでなく CLI 経由の `full`・`soft` を回帰試験に含めること。

2. **Major — 層流 NS を Euler と誤分類し、推奨設定の強制チェックを外す**

   根拠: [rerun_conditions.py:694](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:694) の判定は、`prepare_info.viscous` または SST 名だけです。任意ファイルである `prepare_info.json` がない層流 run は、粘性が有効でも Euler 扱いになります。

   `run_0094` の読み取り内容をメモリ上で「`turbulence.model: none`、`prepare_info` なし」に変更して確認しました。`physProp.viscMethod: 2` のまま、Pt 変更を次の状態で受理します。

   ```text
   runner: run_staged
   stages: none
   config_effective: cfl=5, nStepOuter=6000, override=false
   ```

   **対案:** 粘性・熱伝導・乱流の実効 config から分類し、補助メタデータは照合用にすること。判別できない入力は作成前に拒否する。定数粘性の `viscMethod: 0` もあるため、単純な `viscMethod != 0` 判定では不十分です。

3. **Major — 検証した条件を超えて起動レシピを推奨し、一部は実行不可能**

   根拠: [rerun_conditions.py:674](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:674) と [同:699](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:699)。

   - Euler は変更キーが `{Pt, Ps}` 内なら、スケール有無・背圧比を確認せず `stages: none` を推奨します。実際に `run_0086 --Pt 4.4e6 --keep-Ps`、既定の `scale_ic: none` でもこの分岐に入ります。腕 E が検証したのは、**保存量と背圧をともに0.8倍した条件**です。
   - NS で Pt と Tt を同時変更すると `--scale-ic pt` を推奨しますが、その指定は [同:568](/home/sano/work/forge-integ-1005/solver_density_cuda/tools/rerun_conditions.py:568) の併用禁止条件で拒否されます。組成変更・凝縮 block との組合せも同様です。

   **対案:** 推奨の適用条件をコードで明示すること。Euler の `none` は検証済みのスケール条件に限定し、それ以外は `full`・未検証とする。スケール禁止条件ではスケール推奨を出さず、複合変更の整定実績がないことを記録する。

4. **Minor — 「6000 step ごとの増分」と窓差を取り違えている**

   根拠: [plan:181](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:181) の流量増分 `+2.9 → +1.4 → +0.95 → +0.82` は、[rerun_audit.py:25](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/rerun_audit.py:25) が計算する**各 run 内の隣接5枚窓の平均差**です。

   `run_0135`〜`run_0138` の CSV から再計算した、6000 step 離れた末尾平均同士の差は次です。

   | 量 | blk1→2 | blk2→3 | blk3→4 |
   |---|---:|---:|---:|
   | `mdot` | +4.26328 | +2.42617 | +2.04570 |
   | `delta_E` | −0.00374432 | −0.00181572 | −0.00113889 |

   **対案:** 窓幅と差の定義を訂正する。「残り約0.1 %」は未確定の漸近値との差ではなく、簡易予想値17150との差として表記すること。`DRIFTING` の記録から真の整定余量は確定できません。

5. **Minor — 最新の訂正が現在仕様と残作業表に行き渡っていない**

   根拠: [methods/design/overview.md:977](/home/sano/work/forge-integ-1005/methods/design/overview.md:977) は Euler の実績を一般化した記述のままで、Tt・組成変更時の未確立条件を説明していません。[plan §5.1:86](/home/sano/work/forge-integ-1005/plans/active/tooling-rerun-conditions.md:86) は原データ監査可能化を未了としたままです。また #11 には、今回判明した「Euler 参照を整定させてから δ_E を再抽出する」という作業が明記されていません。

   **対案:** 現在仕様・手順と推奨分岐を同期し、#7 を実績に合わせて更新する。#11 に Euler 参照の準定常確認と δ_E 再評価を追加すること。生産利用への持ち越し自体は、記録されたユーザ決定どおりで構いません。

数値の照合結果は以下です。原データは `/home/sano/work/forge/case/45.isobutane_m6_d155/`、run 索引は [case README](/home/sano/work/forge-integ-1005/case/45.isobutane_m6_d155/README.md:96) です。`check_convergence.py` と指定閾値の `check_quasisteady.py` を再実行しました。

| `case/45.isobutane_m6_d155/` 配下の run | 残差 VERDICT | 対象量の VERDICT |
|---|---|---|
| `run_0119_rerun_ctrl` | `NOT CONVERGED` | δ_E・出口 M・流量 `STEADY` |
| `run_0120_rerun_euler_pt08` | `NOT CONVERGED` | 出口 M・流量 `STEADY` |
| `run_0132_rerun_pt08_scale_cfl1_ext` | `NOT CONVERGED` | δ_E・出口 M・流量 `STEADY` |
| `run_0134_rerun_euler_tt1500` | `NOT CONVERGED` | 出口 M・流量 `DRIFTING` |
| `run_0138_rerun_fullpath_blk4` | `NOT CONVERGED` | 出口 M `STEADY`、δ_E・流量 `DRIFTING` |
| `run_0139_rerun_pt08_noscale_cfl1_ext` | `NOT CONVERGED` | 出口 M・流量 `DRIFTING` |

取得済み主要5場の `VALUE/*` に非有限値はなく、出口 M の再抽出値は CSV と一致しました。B3 最終場の逆流779節点も確認しました。したがって、**A3 の対象量の準定常化と、B3 の規定時間内の準定常未達という限定した結論は支持します**。

**推奨は、移動を保留し、1→2→3の順で修正・回帰確認した後、4・5を同期して再レビューすることです。** 数値カーネルの変更は今回の diff にありませんが、既存 CLI の回帰があるため受入れ不可です。ファイルは変更していません。書き込みを伴う単体試験一式は再実行していません。

指摘数: Critical 0 / Major 3 / Minor 2
