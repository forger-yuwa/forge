# 諮問: 出口較正 (`Md_moc_offset`) を細分格子の Euler でやり直す計画 (§5.1 #11f 草案) の妥当性

日付 2026-10-06。codex (diagnose)。エスカレーション条件 1 (検証計画の新設) と 4。ユーザ決定「出口較正を細分格子に合う Euler でやり直す」(案 1)。
plan: `plans/active/tooling-nozzle-cfd-pinned-initial-line.md` §4.8 (出口較正)、§5.1 #11〜#11f、§9 末尾 (#11e まで)。前回諮問 `notes/reviews/2026-10-05-m6-finemesh-exitM-diagnose.md`。

## 観測事実 (要約)

- 生産 NS (細分格子 ni 2000 × nj 97) の最終 ③ `run_0107`/延長 `run_0109`: 出口コア M 5.998329 (腕 A、STEADY) で登録 6.000 ± 0.02 % を FAIL。粗格子 NS (run_0094) は 5.99926。他ゲート (出口半径・δ_E/δ_C 1.0008・波 η0.1 0.0066 %・オーバーシュート η0.1 −0.002 %・壁解像) は合格。
- 標本位置 A/B (#11e): 粗格子の η 節点へ補間すると +0.000194 (5.99852)。不足の約 1 割。
- 現行の出口較正: Euler (slip、ni 1100 × nj 65、第 1 セル 5e-3、2 次 SLAU・cfl 2→6?・soft 段) で CFD ピン壁の出口コア M 6.000416 を 1 係数で 6 に戻した (`Md_moc_offset` −4.16e-4、run_0083〜0088 で 6.000000)。凍結線は run_0062 (同じ粗い Euler 格子) から。
- Euler の出口コア M は格子で変わるか未測定。NS は粗→細で −0.00093。

## 草案 (§5.1 #11f)

E1: 現行壁の Euler を G0 (較正格子、既存)、G1 (生産 NS 細分格子と同じ格子パラメータ — 出口断面の節点配置を一致させる)、G2 (第三水準 3000 × 145) で各 1 本。予測: H1 (Euler 較正の格子持越し) なら M_E(G1) − M_E(G0) ≈ −0.00093 ± 50 %。
判定 (a) 予測帯 → 新 offset = −4.16e-4 − (M_E(G1) − 6)、E2 で新壁 Euler(G1) 6.0000 ± 1e-4 → NS 連鎖 ①②③ (#11 と同じゲート) ; (b) |G1 − G0| < 2e-4 → H1 棄却・諮問 ; (c) その他 → 諮問。
凍結線 (run_0062、粗い Euler) は変えない (§4.4 の定義「同じ形・ガスで Hall 初期線の V0 型壁を Euler で解いた場」)。

## 問い

Q1. この計画で H1 を判別し、出口 M を合わせられるか。抜けている交絡 (例: 凍結線の格子依存、Euler の第 1 セルが極薄で slip 壁の数値誤差が変わる、cfl 設定の違い、バイナリ変更 [旧 93960c0b → 新 6b47811b、同じソース])。
Q2. G1 に「NS と同じ格子パラメータ (極薄の第 1 セル)」を使うのは妥当か、Euler 用の格子 (2000 × 97・第 1 セル 5e-3) の方がよいか。
Q3. 予測帯・判定閾値・E2 の許容差 (1e-4) と NS 連鎖の予測 (r_t・k_f は #11 ② と ±0.3 % 内) の妥当性。G0 を新バイナリで回し直す必要はあるか。
Q4. 新 offset は 1 係数の線形較正 (V3′ と同じ) でよいか (ΔMd と出口 M の関係が 1:1 である根拠は V3′ の実測)。

## 読んでよいもの

上記 plan、`case/45.isobutane_m6_d155/{README.md,prep_wallfit_euler.py,problem_d155_euler_c2final_n2400_pincal.yaml,problem_d155_ns_finemesh_pin.yaml,run_finemesh_final.sh,exitM_sampling_ab.py}`、`design/forge_design/evaluate/runner_axismach.py`、`design/forge_design/report/nozzle_report.py`。run は AWS のみ。編集禁止。
