# 諮問: B4 — 全ヘキサ接続模型を SERN 全体 (有限厚カウル後縁込み) へ広げる設計と段階分け (2026-10-06)

関連 plan: `plans/active/tooling-sern-mesh-blocking.md` (§4.10 接続設計 19 ブロック/断面 × 3 区間、§4.11–4.13、§5.1 B1b 済・R1・B4–B6)、
`plans/active/tooling-nozzle-sern-chain.md` §5.1 R7b 4-2 (ユーザ決定「カウル後縁を有限厚」、経路 T2)、`plans/active/tooling-nozzle-sern-3d.md` (現行メッシャ・領域・帳簿)。
前回諮問: `notes/reviews/2026-10-05-sern-cowl-blunt-te-route-diagnose.md` (T2 採用・t_te 0.005 H 暫定・受入条件)、`2026-10-05-sern-b1b-regression-diagnose.md`。

## 事実
- 接続模型 `case/46.sern_design/cad/hex_junction_model.py` (550 行、gmsh transfinite): ランプは 2 次式 y_r/H = 1 + 0.30 x/H − 0.02 (x/H)²、x ∈ [0, XEND 2.4 H]、側壁 (厚み TSW 0.05 H) とカウル (TC 0.02 H) は有限厚、端面 (sidewall_end・cowl_base) に壁層、断面の隅フィレット (PROFILE)、z は対称面 0 〜 ZFAR 2 H、y は YBOT −1 H 〜。ブロック: N (ダクト内バタフライ)・SW (側壁の跡)・CW1/CW2 (カウルの跡)・U (カウル下)・S (側壁の外)。B1b 修正後の変換器で閉性 PASS。scale 0.5/0.71 で品質・閉性・Jacobian PASS、規模は第一層 16 µm・端面 0.25 mm で 245 万節点 (B1c)。上流の外部流・機体上面/側面/ベース・プルーム全長・MOC 輪郭は持たない。
- 現行テンソル積メッシャ `design/forge_design/meshing/mesh_sern3d.py` の領域とタグ: inlet_nozzle・inlet_ext (上流面)・outlet・ramp・cowl_in・cowl_out・bottom・top_out・sym・side_far (farfield)・sidewall_in/out・vehicle・vehicle_top・underside_far・vehicle_side・vehicle_base。上流 L_up 0.5 H、プルーム x_out = L_ramp + x_out_extra 2 H、下 bot_depth 3 H、上 ext_top (機体上面 + 自由流バンド)、側方 Z_ext 1.5 H (+ farfield、必要幅 2.50 H)、機体ベース t_base 0.02 H。生産形状は MOC 設計 (L_ramp ≈ 10 H、L_cowl 1.2 H、L_sw 0.8 H exact)。
- 帳簿 `runner_sern3d.forces3d` は新端面タグを集計せず入口項は矩形面積 (B5 で改修)。
- 受入条件 (前回諮問で確定): float32 後の全頂点 Jacobian 正・双対体積正・全 CV 閉性 ≤ 1e-5・タグ漏れ/重複 0・skew ≤ 0.90・AR 既定、後流はベース直後 t_te の範囲で Δx ≤ t_te/5・間隔比 ≤ 1.2、端面の第一内部点距離と y₁⁺ を別検査、細分列でも端面解像を確実に変える。
- AWS g5 (RAM 16 GB、3D 変換は 1 本ずつ、available ≥ 9 GB)。

## 問い
1. B4 の最小の段階分け: 例えば (B4a) 接続模型のランプ/カウルを MOC 輪郭 (設計の ramp_xy/cowl_xy、L_ramp ~10 H) に置き換え、有限厚カウル後縁 t_te を入れる → (B4b) 上流 (inlet_nozzle/inlet_ext) と下流プルーム全長・出口 → (B4c) 機体上面/側面/ベース (ext_top 相当) と側方 farfield 幅 2.50 H → (B4d) runner への `mesh3d.backend: blocks` 接続。各段の受入ゲートと、段ごとに CFD を挟むか (メッシュだけで進めるか)。
2. 生産規模の見積もり方法 (節点数・変換メモリ) と、16 GB で載らないときの選択肢 (AWS の大きいインスタンス、変換器の省メモリ化、ブロック粗化)。
3. 現行メッシャとのタグ・帳簿の互換 (同じ physID を使うか、新タグ cowl_base/cowl_side/sidewall_end の扱い)、形状パラメータの持ち方 (問題 YAML の mesh3d 節)。
4. 断面隅フィレットを初版で入れるか (前回: 初版でフィレットは同時に変えない)。
