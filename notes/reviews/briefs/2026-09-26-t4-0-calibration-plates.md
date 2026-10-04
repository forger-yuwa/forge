# 諮問ブリーフ: T4-0 校正平板の再現の設計 (plan §4.12・§6 G14–G16) を確定してよいか

- 日付: 2026-09-26 / ブランチ `feature/gap-heating-precision`
- plan: `plans/active/case-hypersonic-gap-heating-validation.md` の **§4.12** (新設)、§6 の **G14–G16** (新設)、§5.1 **#61**
- 精読メモ: `notes/investigations/tp1187-test-system-reading.md`
- 一次資料 (PDF、読んでよい。必要な頁だけ): `papers/gap_heating/Deveikis_Hunt_1973_NASA-TN-D-7275_large_flat_plate_M7_8ftHTST.pdf`
  (本文 p.4–21、Table II p.26 相当、Fig 19 p.58 相当、Fig 22 p.61 相当)、
  `papers/gap_heating/Cary_1970_NASA-TN-D-5863_turbulent_heat_transfer_cooling_M6.pdf` (本文・結論・Table II)、
  `papers/gap_heating/Avery_1978_NASA-TP-1187_gap_heating_laminar_turbulent.pdf` (p.5–9、Fig 9/10)。
- 読まないこと: plan の §4.12・§6・§5.1 #60/#61 以外 (plan は 1000 行超、#60 の行は数万文字)。

## 1. 観測事実

- case/56 の 3D すきま照合 (TP-1187 90° 配列) は観測モデルを直しても D1 が全 6 点で実測の 1/2.2〜1/25。原因未特定。
- 精読で判明: TP-1187 の**乱流 q_FP は TP-1187 の測定ではない** (p.8「obtained from the experimental data of reference 5」、ref.5 = TN D-7275)。
  run 8 と 14 は Tt が違うのに q_FP が同値 63.44、α 7.5° の 6 run も同値 62.99 → D-7275 Fig 22(b) のパネル平均 carpet plot の内挿とみられる (状況証拠)。
- D-7275: 乱流の中心線 St は Eckert 参照温度法より 10–30 % 低い。Cary TN D-5863 (M6、冷却平板、Tw/Tt 0.19–0.70): 冷壁で Van Driest・T' 法が大きく過大、Spalding–Chi が合う (仮想原点を加熱ピークに置いた場合)。
- forge は理論側: case/48 (M4.19, Tw/Taw 0.26) で C_f/VD-II 0.97–0.99、case/56 の 2D 平板 (M7) で q_w/Eckert 1–2 %、実測 q_FP に +19/+22 %。
- D-7275 試験 26 (TP-1187 run 8 とほぼ同条件) は M∞ 6.64。case/56 Gate A は M 7 固定。

## 2. 期待値と出典

上記の各報告。plan §4.12 に表として整理した。

## 3. 再現条件 / 4. 実施済み

- まだ計算していない。§4.12 は設計のみ。
- ユーザの指示と承認: 「まず calibration panel の結果や q_FP を再現できるかから」「おｋ頼む」(D-7275 + Cary の 2 本立てを提案して承認)。

## 5. 仮説

- H-model: SST (低 Re、壁関数なし) は冷壁の乱流熱伝達を過大に出す。
- H-facility: 当時の測定・自由流還元・トリップに共通の系統誤差で実測が低い。
- 両方が同時に効いている可能性もある。

## 6. 聞きたいこと

1. **§4.12 の設計と §6 G14–G16 の判定規則で H-model と H-facility を分けられるか**。特に G14 の「R の Tw/Tt 依存 ≤10 % / >10 %」の閾値と、
   G15 の「|B/R(0.2) − 1| ≤ 0.10」の閾値が妥当か。事前に決めるべき量・前提ゲートの抜けがあれば指摘してほしい。
2. **Cary を判別に使うことの落とし穴** (低温 T∞≈65 K の空気、自然遷移の扱い、仮想原点の置き方、液体窒素冷却の壁温分布など)。
3. **T4-0b で 2D・トリップ無し (起点 2 通り) で挟む近似**が許されるか。
4. 順番 (Cary → D-7275 → TP-1187 層流) に反対があれば。
Critical/Major は採否を plan に書く。1 つに絞る必要はないが、重い順に並べてほしい。
