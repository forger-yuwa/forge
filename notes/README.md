# forge 調査・作業ノート (`notes/`)

`notes/` は恒久仕様 (`methods/`) でも変更単位の設計判断 (`plans/`) でもない、**調査メモ・作業ログ・一時的知見**の逃がし先。
ここに置いたものは「現在地」を示す `plans/` の一覧を汚さない。

- `investigations/` — 技術調査・サーベイ・外部ソルバ棚卸し (コード変更を伴わない研究)。結論は該当 `plans/` 計画や `methods/` に反映する。
- `sessions/` — 使い捨ての作業セッション用プロンプト・引き継ぎメモ。恒久参照を意図しない。

## investigations

| ノート | area | 概要 |
| --- | --- | --- |
| [chemistry-finite-rate-h2-survey.md](investigations/chemistry-finite-rate-h2-survey.md) | `thermophysics / chemistry` | 有限速度化学 (H₂ 燃焼・ノズル上流の化学非平衡) 導入の文献調査と方針: 機構選定 (Jachimowski→Burke)、剛性処理 (decoupled point-implicit 種ブロック)、sensible datum + 反応熱陽注入、Phase 0–4 計画と検証 (Cantera/CEA/Burrows–Kurkov) (2026-09-04) |
| [condensation-droplet-hypersonic-survey.md](investigations/condensation-droplet-hypersonic-survey.md) | `condensation / multiphase` | 極超音速燃焼風洞・飛翔体における凝縮液滴/雨滴の挙動 (風洞凝縮の蒸発vs成長・燃焼器到達 end-to-end・迎角依存) |
| [node-slip-tangential-density-spurious-flow.md](investigations/node-slip-tangential-density-spurious-flow.md) | `boundary / discretization` | node slip 境界 + 接線密度勾配で市松状スプリアス接線流が定在する未修正バグの切り分け記録 (case/24 純伝導検証の副産物) |
| [condensation-nitrogen-nozzle-flows-summary.md](investigations/condensation-nitrogen-nozzle-flows-summary.md) | `condensation` | "On Nitrogen Condensation in Hypersonic Nozzle Flows" 論文精読まとめ + Arthur 1952 形状の forge dry 再現 (PDF: [papers/condensation/](../papers/condensation/)) |
| [nozzle-pintle-literature-review.md](investigations/nozzle-pintle-literature-review.md) | `nozzle` | ピントル/ニードル弁式 超音速ノズル推力調整機構 — 文献・特許まとめ (PDF/成果物: [papers/pintle_nozzle/](../papers/pintle_nozzle/)) |
| [convection-central-scheme-oscillation-control.md](investigations/convection-central-scheme-oscillation-control.md) | `convection` | LES/DES における中心差分 (KEEP) のスプリアス振動抑制 技術調査 (振動源4分類の統括) |
| [convection-pep-scheme-survey.md](investigations/convection-pep-scheme-survey.md) | `convection` | PEP (Pressure-Equilibrium-Preserving) 系スキーム技術調査と forge 実装方針 |
| [nozzle-optimization-tool-survey.md](investigations/nozzle-optimization-tool-survey.md) | `—` | 超音速・極超音速ノズル最適化ツール — 技術動向調査と開発フロー提案 |
| [nozzle-top-internal-shock-diagnosis.md](investigations/nozzle-top-internal-shock-diagnosis.md) | `verification / nozzle` | TOP ベルの内部衝撃波と軸 M 過大 (+22%) は物理と確証 — SU2 同一メッシュ 2×2 (RANS/Euler) で forge と定量一致、Rd/Ru スイープで円弧曲率起因も確認 (2026-08-13) |
| [nozzle-deltastar-throat-review.md](investigations/nozzle-deltastar-throat-review.md) | `design / nozzle` | 排除厚さ補正の総点検: 与えたスロート δ\* が NS 実効値 (0.0015 r_t) の 3〜12 倍過大で NS 質量流量が Euler 設計比 +0.8〜3.7 % → 試験部 M −0.2〜−0.7 % の主因。Md トリム/law 側帰還の順序問題と打ち手 (質量流量ゲート・全域 CFD 抽出・積分法) (2026-09-04) |
| [sern-design-method-survey.md](investigations/sern-design-method-survey.md) | `design / nozzle` | SERN (single ramp) 設計手法の系統調査と ⑤ 方針再検討: 主流は最大推力理論 (Rao 2D) の key point 逆設計で壁圧規定は剥離位置制御の道具、NS 帰還不要、入口は燃焼器出口 (遷音速不要)、2D 設計+3D 確認 (2026-09-04) |
| [nozzle-throat-curvature-shape-representation-survey.md](investigations/nozzle-throat-curvature-shape-representation-survey.md) | `design / nozzle` | スロート曲率を制御する形状表現の文献調査 (case/41 接合こぶ対策) — CONTUR は壁接合を作らず軸分布 C2 で曲率連続を保証、R=2 は遷音速理論の下限、推奨は R=5 判別 A/B → κ(s) 表現 (2026-08-17) |
| [su2-nemo-contact-thermo-investigation.md](investigations/su2-nemo-contact-thermo-investigation.md) | `convection / thermophysics` | SU2-NEMO contact/interface thermo 取り扱い調査 (forge mixed-order face-state 比較) |
| [turbulence-des-flux-survey.md](investigations/turbulence-des-flux-survey.md) | `turbulence` | DES/DDES/IDDES 用 低散逸対流 flux 設計 技術調査（圧縮性 LES/DES） |
| [turbulence-des-wmles-survey.md](investigations/turbulence-des-wmles-survey.md) | `turbulence` | DES / WMLES 最新動向サーベイ（超音速ノズル・ピントルバルブ適用向け） |

## sessions

| ノート | 概要 |
| --- | --- |
| [2026-09-05-handover-sern-design.md](sessions/2026-09-05-handover-sern-design.md) | ⑤ SERN 設計チェーンの引き継ぎ (branch feature/sern-design, case/46): 状態表・使い方・罠 8 件・未解決 9 件 |
| [condensation-nonequilibrium-session-prompt.md](sessions/condensation-nonequilibrium-session-prompt.md) | 次セッション用プロンプト: 非平衡凝縮 (4 モーメント方程式) の forge 実装 |
| [cutler-tp-multispecies-cfl-analysis-prompt.md](sessions/cutler-tp-multispecies-cfl-analysis-prompt.md) | 分析依頼プロンプト: 多成分 TP の陰解法 (block-DPLUR) で安定 CFL が ~1 に頭打ちする原因 |
| [implicit-acceleration-session-prompt.md](sessions/implicit-acceleration-session-prompt.md) | 別セッション用プロンプト: 陰解法の安定 CFL 引き上げ |
| [precision-mixed-axisym-session-prompt.md](sessions/precision-mixed-axisym-session-prompt.md) | 新セッション引き継ぎプロンプト — 軸対称 近軸の陰解法を混合精度で root-fix |
| [species-in-dplur-session-prompt.md](sessions/species-in-dplur-session-prompt.md) | 引継ぎプロンプト: 化学種 `roY_s` を block-DPLUR に結合する (roe↔roY coupling) |
