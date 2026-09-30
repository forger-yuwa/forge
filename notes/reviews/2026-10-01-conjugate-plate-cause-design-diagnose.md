# codex 諮問 (diagnose): conjugate-plate-cause-design

- **brief**: [`notes/reviews/briefs/2026-10-01-conjugate-plate-cause-design.md`](../../notes/reviews/briefs/2026-10-01-conjugate-plate-cause-design.md)
- **plan**: [`plans/active/boundary-cht-conjugate-flat-plate.md`](../../plans/active/boundary-cht-conjugate-flat-plate.md)
- **date**: 2026-10-01
- **commit**: `be84a446` (feature/cht-conjugate-benchmarks)
- **codex**: effort `high`, 5.5 min, rc=0
- **結論**: **n64 C1 の同一チェックポイントから「連成継続／壁温分布固定」の B1 を対照付きで行い、壁温更新が流体停滞の維持に必要かを判別する。**
- **extra**: `case/65.conjugate_flat_plate/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表：

| 論点 | 採否・重大度 | 根拠と対案 |
|---|---|---|
| B1：壁温分布を固定 | **採用。ただし再開条件の補完が必要〔Major〕** | 手元の `run_0007_c1_n64/mesh.h5` の保存量は `float32`。通常のコピーでは FP64 場を丸めるため、[`restart_field.py:31`](/home/sano/work/forge-cht/solver_density_cuda/tools/restart_field.py:31) の `--keep-src-dtype` が必要。さらに壁温分布と固体状態・更新位相を同時刻から復元する。固体状態の復元は [`conjugateWall.cpp:578`](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:578)。 |
| B2：`Df_scale` を20へ | **要再検証。初手には選ばない〔Major〕** | `Df` は更新の前処理に入る（[`conjugateWall.cpp:957`](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:957)）。壁温の振幅が小さくなっても、流体側の揺れを減衰させただけかもしれない。先に B1 で、壁温更新なしでも流体停滞が持続するかを判別する。 |
| B3：slip を no-slip に変更、領域短縮 | **今回の判別案として却下〔Major〕** | 上流も no-slip にすると板前縁へ流入する境界層が変わり、領域短縮も入口・出口の影響を変える。登録問題は [`boundary-cht-conjugate-benchmarks.md:71`](/home/sano/work/forge-cht/plans/accepted/boundary-cht-conjugate-benchmarks.md:71)。改善しても slip 実装の欠陥を特定できない。境界条件を保持した B1 を使う。 |
| H1：④の丸め床と閾値相対化 | **丸めの寄与は採用。「達成不能」の断定と即時相対化は却下〔Major〕** | ④は非界面節点の `Au−b` の最大値で、大きな Robin 項を差し引く（[`conjugateWall.cpp:953`](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:953)、[`solidFem2d.cpp:205`](/home/sano/work/forge-cht/solver_density_cuda/conjugate/solidFem2d.cpp:205)）。まず残差評価の丸めと保存温度の量子化を分離する。閾値を改訂する場合も旧判定を残し、**事後改訂した別基準**として扱う。 |

結論: **n64 C1 の同一チェックポイントから「連成継続／壁温分布固定」の B1 を対照付きで行い、壁温更新が流体停滞の維持に必要かを判別する。**

第 1 仮説: **固定した壁温の下でも流体側の停滞が持続する。slip 欠陥はその候補だが、まだ特定できない。** 確度: **中**

  根拠: `case/65.conjugate_flat_plate/run_0006_c1_n32/` では、G-if 判定区間 `596000..599950` の `dTw_max ≤ 7.00e-5 K` に対し、流体は `rms_ro` 1.6桁、`rms_roUy` 1.8桁、`rms_roe` 1.4桁で停滞し、保存判定は **`NOT CONVERGED`**。大きな壁温振動がなくても同型の流体停滞がある。ただし、これだけでは連成から独立とは証明できない。

  反証条件: 正しく再開した連成対照では停滞が再現する一方、壁温固定側だけで停滞していた全残差が持続的に減衰し、局所的な場の変動も消えること。

第 2 仮説: **連成フィードバックが n64 の振れを維持・増幅する**。確度: 中。`run_0007_c1_n64/` の末尾400更新を再集計すると `res_abs_Wm2` は **13.7～741.8、平均182.8**。発生源が連成か流体かは未確認。

第 3 仮説: **④には Robin 項の丸めが支配的に寄与する**。確度: 高。ただし①～③・流体停滞を説明する仮説ではない。実メッシュと C++ の組立て・加算順を Python で再現した数値検査では、厳密平衡の一様310 K・界面荷重ゼロでも、非界面残差は n16 C1 **`3.725e-9 W/m`**、n64 C1 **`9.313e-10 W/m`**となった。これは実 run の固体温度を再評価した結果ではない。

判別 A/B:

- **比較する変更は壁温更新の有無だけ。** 同一の `res_600000.h5`、固体チェックポイント、壁温分布、同一 FP64 バイナリから開始する。両側で `wallProfile: 1` により同じ壁温分布を読み込む。単に `conjugate: 1` を外すと、現在の設定では一様 `Ts: 300` に戻り、別実験になる。
- **A＝連成継続、B＝壁温固定。** B は `conjugate.warmup` を再開後の累積終了stepより大きくする方法で実現できる。更新前にこの条件で戻るため、変更を1キーに限定できる（[`conjugateWall.cpp:1080`](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:1080)）。
- **長さはまず各60000 step**。20000 stepずつの3区間で評価する。提案先は `case/65.conjugate_flat_plate/run_0011_ab_c1_n64_coupled/` と `run_0012_ab_c1_n64_fixed_tw/`（未作成）。全残差、界面全節点の熱流束、前後縁付近の `P/T/U_n` の時系列を採る。既存の12010 step間隔だけでは短い診断には疎なので、両側共通で出力頻度を上げる。
- **事前登録する判別条件**：A が元の停滞を再現することを前提とし、B の末尾2区間で、停滞していた3列 `rms_ro/rms_roUy/rms_roe` がすべて A の **1/10以下**となり、局所変動も減衰すれば、**「流体側だけで同程度の停滞を維持する」を棄却**する。反対に、B でも3列が A の **1/2以上**に残り、最後の3区間の代表値が10%以内の横ばいなら、**「動的連成だけが停滞を維持する」を棄却**する。後者でも slip 原因とは断定しない。
- この比率は**診断用で、収束合格条件ではない**。両側を `check_convergence.py` で判定し、準定常も対象量を `check_quasisteady.py` に渡す。元 run は未収束なので `--from-floor run_0007...` は使えない（[`check_convergence.py:353`](/home/sano/work/forge-cht/solver_density_cuda/tools/check_convergence.py:353)）。中間的な結果・減衰途中なら120000 stepまで同じ比較を延長し、それでも分かれなければ判定不能とする。

やらない方がよいこと: **B1 と同時に `Df_scale`、slip、CFL、閾値を変更しない。** ④を `hT_hΔx` で割って合格化しない。この尺度は大きな相殺項であり、検証したい伝熱荷重そのものではない。残差評価を温度差ベースにして高精度評価と照合し、それでも基準改訂が必要なら、旧 FAIL と新判定を併記する。

呼び出し側の前提への異議:

- **〔Major〕slip 原因を示す観測と、既知欠陥の症状がまだ一致していない。** 既知メモは「市松状の**接線速度**」（[`調査メモ:9`](/home/sano/work/forge-cht/notes/investigations/node-slip-tangential-density-spurious-flow.md:9)）。今回挙げられたのは法線速度である。また `slip_d` は境界状態の法線速度を除去するが、内側 DOF を直接ゼロにする処理ではない（[`boundaryCond_d.cu:68`](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/boundaryCond_d.cu:68)）。節点 `Uy` だけを壁面漏れの証拠にせず、座標・法線・境界流束で確認する。
- **〔Minor〕「`2^-31` の整数倍」は不正確。** `6.985e-10` はその1.5倍。実履歴の離散値は `2^-32` の整数倍と整合する。n64 C1 は後半更新の **96.7%**が④の閾値内なので、「単発でも達成不能」ではない。問題は80回連続判定の頑健性である。
- 評価窓の保存系列は再実行した `check_quasisteady.py` でも **`OVERALL: ALL STEADY`**。これは窓外の前後縁を除外しない。主参照も forge の流れを固定した伝熱検証なので、流れの正しさを保証しない（[`発注元plan:50`](/home/sano/work/forge-cht/plans/accepted/boundary-cht-conjugate-benchmarks.md:50)）。
- 「亜音速に超音速専用入口を使用」は今回は採らない。実装には亜音速用の特性分岐がある（[`boundaryCond_d.cu:808`](/home/sano/work/forge-cht/solver_density_cuda/cuda_forge/boundaryCond_d.cu:808)）。手元メッシュも `z=0` の平面で、押し出し2列ではなかった。

不足情報: **同時刻の流体・固体チェックポイント、元の `residual_history.csv`、界面節点ログ、前後縁の場の時系列。** 手元の `check_convergence.py` 再実行は `NO residual_history.csv` となり、流体収束は保存済み `CONVERGENCE_CHECK.txt` の確認まで。局所変動位置と NaN/Inf の全場検査は独立に再確認できていない。AWS は起動していない。

**plan 未反映**。ファイル変更禁止のため、呼び出し側で対象 plan の §4・§6・§5.1 に事前登録する提案である。
