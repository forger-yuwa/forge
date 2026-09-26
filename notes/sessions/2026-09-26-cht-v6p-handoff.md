# 引き継ぎ: CHT の V6′ を閉じる (2026-09-26)

**残作業の正本は [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md)
の §5.1 (特に #100 / #101) と §6 V6′。本ノートは写しとポインタだけ**なので、着手前に plan を読むこと。

## 0. いまどこにいるか

- **ブランチ `feature/cht-phase2-fem2d`** (ワークツリー `/home/sano/work/forge-cht`)。`main` 未マージ。
  `feature/sern-design` へのマージは 2026-09-24 に済み、`main` はユーザ判断で保留。
- **CHT plan を閉じる条件は 3 つに固定されていて、2 つは埋まった**:
  1. **V2 の `case/48` 段の判定** → **済**。(i) 0.5009 % PASS / (ii) 1001 点中 1 点 FAIL / (iii) 0.573 % PASS。
     **判定は FAIL のまま残し、(ii) の 1 点をユーザ決定 (2026-09-26) で例外として受け入れた** (§6「閉じる条件」1)。
  2. **V6′** → **判定不能**。**これが唯一の残り**。
  3. **放射** → ユーザ決定「まだやらない・将来も否定しない」で決着 (§8-3)。
- **V6′ 以外の §5.1 は「記録」が大半**。95→101 行あるが、閉じる条件は上の 3 つだけ。
  **表が長いことを「残作業が 100 件ある」と読み違えないこと** (整理そのものも未着手の項目)。

## 1. V6′ の現状 — 「判定不能」であって FAIL ではない

ケースは **`case/58.conjugate_slot`** (今回作った自己完結の検証)。深いスロット
($W$=1 mm × $D$=20 mm、$D/W$=20) の**前壁だけ共役**、後壁は固定等温 500 K、底と板は断熱、M=2 層流・定数物性。
漸近解 $q_*$=**4709.0 W/m²**、$T_{w1,*}$=**394.180 K**、固体の温度上昇 **94.18 K**。

**FP64 run `run_0010_v6p_fp64` と感度 2×2 は `eval_v6p.py` で 6 条件 PASS**。深部は解析解と
$T_{w1}$ **+0.009 K**・$q_{w1}$ **−0.025 W/m²**。それでも総合を PASS にしていないのは 2 つの理由:

1. **(e) 「連成の保存」の定義が誤っていた** — `conjugateWall.cpp`:818 で **`q_iface` は固体が受け持つ熱ではなく
   流体荷重 $Q_f$ のコピー**。現行 (e) は**渡した荷重を渡した荷重と比べていた**。
   正しい量 $Q_{\rm sol}=(K_su-b_s)_{\rm iface}$ で組み直すと `run_0010` 0.0025 % だが
   **`df20_i200` は 0.4607 %** で許容 0.1 % を超える。
2. **登録した「帯内 G-if の 80 更新連続」の証拠が無い** — 節点ごとの界面残差は
   `conjugateWall.cpp`:741 で計算されるが**保存されるのは全域 max だけ**。

## 2. 次の一手 (正本は §5.1 #101。codex diagnose の全件採用)

1. **`eval_v6p.py` の (e) を直す** (`case/58.conjugate_slot/eval_v6p.py`:144)。固体の物理作用素から
   $Q_{\rm sol}=(K_su-b_s)_{\rm iface}$ を**符号を保って**組み、$Q_f$ と比べる。
   参考実装は `solver_density_cuda/conjugate/` と、外部ループ側の `Fem2DOperator.assemble_full`。
2. **帯内 G-if を毎更新で測る診断出力**を足す。各更新の $r_i, Q_{f,i}, A_i, \Delta T_i$ を
   **節点 ID + 更新番号つき**で記録し、固定した帯の**最後の 80 更新**で判定する。
   **5000 step 間隔の固体ダンプに残差を足すだけでは不足** (更新は 50 step ごと)。
   固体内部残差の検査は全域で維持。
3. **`check_quasisteady.py --series-csv`** (`:519`) で帯平均**と局所量**を判定する。
   **拡張は不要** (本体を触る必要はない)。温度は `Tw−300` を渡し、394 K を分母にして許容を緩めない。
   **帯平均の `ALL STEADY` は局所の静定を保証しない** — 実測で `20/50` と `5/200` は各 13 点、
   `20/200` は温度 7 点・熱流束 6 点が `DRIFTING`/`TRANSIENT-UNSETTLED`。
4. **感度を共通帯・節点ごとに取り直す** (基準は `Df_scale` 5 / `interval` 50)。
   共通帯 211 行での現状: 温度差 `20/50` 0.0231 / `5/200` 0.0211 / **`20/200` 0.2186 %**。
5. そのあと **codex の result レビュー 8 巡目** → `status: done` → **`plans/accepted/` へ移動** +
   `plans/README.md` 同期。

## 3. 踏み抜きやすい罠 (今回すべて実際に踏んだ)

| | 内容 |
| --- | --- |
| **IC は変換器が書く** | `initial` を変えたら **`convertGmshToForge` をやり直さないと効かない**。残差が 1 桁も変わらないのが手がかり |
| **IC と入口の不整合** | `uniform_p101325_u10` は P=101325 / u=10。M=2 の入口に使うと $p$ が 2.18e7 Pa に。case ごとに `setInitial.hpp` へ追加する (`slot_m2` が例) |
| **安全停止の向きの前提** | 固体温度は **[min($T_c$)−20, 流体の最大全温+20]** に制限される (`conjugateWall.cpp`:774-800)。**背面加熱でガスより熱い固体は拒否される** |
| **restart 元と本段の壁温** | 違うと最初の結合更新が過大 (固体が 712 K に跳ねた)。**spinup の壁温を本段と揃え、`warmup` を入れる** |
| **共有角ガード** | 角 CV を複数の等温壁が矛盾する $T_s$ で共有すると起動時に拒否 (正しい挙動)。片方を断熱にする |
| **`Df_scale` 1.0 は発散** | lip の再入角で $\partial Q_f/\partial T_w$ が $D_f$ の 15 倍。**`Df_scale` 5 以上**。**`interval` を短くするのは逆効果** (ソルバのエラーメッセージのヒントはこの機構では誤り) |
| **`check_solver_config.py` の偽陽性** | `conjugate.mode/flux/gate` を「無視される」と警告するが**すべて読まれている** (`solverConfig.cpp`:509-511)。別項目として未起票 |
| **`codex_review.py` の版** | `forge-cht` は cherry-pick で `--stage diagnose` に対応済み (`ba5a42aa`/`6d7e1af5`/`e19f8acb`)。**本体ツリーのスクリプトを直接呼ばない** (codex が本体ツリーを読む) |

## 4. 諮問の経路 (2026-09-26 時点)

**既定は Fable の `diagnostician` サブエージェント、代替が codex (`--stage diagnose`)。**
`~/.config/forge/diagnose-backend` が `"codex"` のとき、または Fable が usage 上限で失敗したときだけ codex を使う。
**現在のスイッチは `codex`** なので `codex_review.py --stage diagnose --brief <brief.md> [plan]` で回す。
記録は `notes/reviews/<日付>-<slug>-diagnose.md`、応答には
「codex (diagnose) に諮った: `<記録>` — 結論 1 行」を書く。

## 5. V6′ の後に控えているもの (正本は §5.1 #95 / #96、§8-0)

**ユーザ決定 (2026-09-26): 軸対称 (#95) → 3D (#96) の順。本計画は V6′ で閉じ、#95/#96 は別 plan に切り出す。**

- **#95 軸対称の `fem2d`** — これが入ると**ノズル CHT** ができる (壁厚 + 外表面に HTC を貼る = `hole` の Robin)。
  `case/40` が壁温 1400 ± 15 K を**仮定**している所を解にできる。**円筒殻の 1 次元伝導に解析解がある**ので
  V1 と同型の対照が取れる。**最初にやる確認**: 現状の軸対称が「起動時に拒否」なのか「黙って間違う」のかが
  **未確定** (`conjugateWall.cpp`/`solidFem2d.cpp` は `isAxisymmetric` を一度も見ていない)。
  「黙って間違う」なら**即座に拒否を入れる**。
- **#96 3D の `fem2d`** — `case/49` 環状キャビティ・SERN・3D ノズルを閉ざしている最大の欠落。
  **複数壁の連成** (`solid:` が単一文字列なので今は 1 枚だけ) と**放射**も揃わないと `case/49` には届かない。

## 6. 別系統の未完 (この plan の外)

- **`ypls`** の #12f (カーネルとツールの幾何・未評価仕様の統一。直角二壁の交点で反例あり) と
  #12g (番兵 −1 の利用側対応。`cavity_eval.py` は case/49 = 別セッション所有なので申し送り) —
  正本は [`tooling-convergence-and-wall-resolution-gates.md`](../../plans/active/tooling-convergence-and-wall-resolution-gates.md) §5.1
- **#94** `qAccumulatorFP64` との比較 (`feature/gap-heating-precision` のマージ後)
- **受け渡し** — [`2026-09-25-cht-gap-application-handoff.md`](2026-09-25-cht-gap-application-handoff.md) を
  case/49 側に渡す話は**書いてあるだけで、まだ声をかけていない**
