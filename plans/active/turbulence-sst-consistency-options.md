# SST 整合オプション (codex レビュー副次指摘の実装と検証)

## メタ

- **area**: `turbulence`
- **status**: `in_progress` (実装済 2026-09-08、A/B 検証中、codex レビュー待ち)
- **related_docs**: `methods/turbulence/implementation.md`, `procedures/solver-settings.md` (turbulence.sst*)
- **related_plans**: `accepted/turbulence-sst-node-corner-heating.md` (§3.5 の副次指摘が出所)
- **created**: `2026-09-08`
- **owner**: `CFD Dev`

## 1. 目的

codex レビュー (2026-09-08) が挙げた SST 実装の標準形からのずれ 5 件 + 単体試験 1 件を、**既定挙動を変えないオプション**として
実装し (例外: node 壁 k ピンは正しさの修正なので既定 ON)、A/B で影響を定量化して既定値を決める。

## 2. 実装 (commit 0fcde108)

| # | キー | 既定 | 内容 | 変更箇所 |
|---|---|---|---|---|
| 3 | `sstNodeWallKPin` | **1** | node 低 Re 壁で k/roK=0 を壁ノードにピン、壁ノードの k/ω 残差・対角を 0 化 | `ransBoundary_d.cu`, `ransSource_d.cu` |
| 1 | `sstOmegaProdFromPk` | 0 | P_ω = α P_k/ν_t (リミッタ後 P_k と整合) | `ransSource_d.cu` |
| 2 | `sstSigmaBlend` | 0 | σ_k/σ_ω を F1 ブレンド (F1 は `sstF1` に保存、1 step 遅れ) | `scalarTransport_d.cu(h)`, `ransTransport_d.cu`, `variables.hpp` |
| 4 | `sstIsotropicStress` | 0 | 応力に −(2/3)ρk δᵢⱼ (内部面のみ。境界面は未対応 = 出口面で ρk 分の力の不整合が残る) | `viscousFlux_d.cu` |
| 5 | `sstEnergyKSource` | 0 | エネルギー式に −(P_k−D_k)V | `ransSource_d.cu` |
| 6 | `tools/test_scale_invariance.py` | — | 相似メッシュ (座標 α 倍, μ/λ/dt/ω を次元どおりスケール) で場の一致を見る | tools |

## 3. 検証結果

### 3.1 case/16 2D node SST (run_0213 の収束場から 3000 step, run_0231_ab2d_*)

壁 p/p0 の base 比: omegaPk +0.03〜+0.11 %、sigmaBlend ≤0.01 %、isoStress +0.09〜+0.20 %、energyK −0.05〜−0.10 %、
all +0.06〜+0.23 %。壁 k ピン (base vs kpin0) は差なし (壁 k は既に 2e-5 以下)。壁 T は全て Tt−0.03 K。
→ 2D 平面壁では全て 0.25 % 以下の影響 (F1=1 域でリミッタ非発動、ρk ≪ P)。

### 3.2 case/26 平板 (run_node_sst_muscl_cont の res_90000 から 5000 step, run_0025_sstopt_*)

Cf/Schlichting @x=0.3/0.6/0.9 (node, 第一内点勾配 `tools/cf_node.py`): kpin0 0.8899/0.9292/0.9553 = **base と完全一致** (壁 k は
既に ~0)、all 0.8897/0.9289/0.9552 (差 ≤0.03 %)。参照 res_90000: 0.8927/0.9260/0.9570。→ 平板 (M0.2, y1+≈0.35) では全て無影響。
(config は旧 `LESorRANS` 体系だったため `model: "sst"` 体系へ書き換え、dilatation 0 は維持。)

### 3.3 case/16 3D 半幅 SST (run_0228 の場から 6000 step, AWS run_0232_ab3d_base / _all)

| | 壁 p/p0 @x=16 / 46 / 85 | T>Tt+1 ノード | 角対角線 (wd 0.6 / 2.5 / 4.8 / 19 µm) k, ω, μt/μ |
|---|---|---|---|
| 0228 (開始場) | 0.3761 / 0.2727 / 0.2063 | 0 | k 6.8e-9 / 2.5e-5 / 1.7e-3 / 6.6; ω 1.8e10 / 1.4e9 / 3.5e8 / 2.1e7 |
| base (既定, +6000) | 0.3761 / 0.2723 / 0.2049 | 0 | 同上 (±3 %) |
| all (+6000) | 0.3770 / 0.2737 / 0.2067 | 0 | 同上 (±3 %) |

all は base 比 +0.2〜+0.9 % (等方応力 + P_ω 整合の分)。角の ω/k 分布・Tt 超えは不変 (角部問題は相対ガードで解決済みで、
本オプション群は寄与しない)。rms_ro は両者 2.2〜2.3e-11 まで低下。

### 3.4 相似メッシュ試験 (`test_scale_invariance.py`, case/16 2D node SST, 1 次風上, 50 step)

α = 0.1 / 0.01: ρ, P, T, u, k, μt は 1e-5 以下で一致 (スケール不変)。**ω だけ喉付近の壁ノードで 2〜3 % 差** (原因未特定: 壁 ω ピン
の y_eff か float32 幾何)。α = 1e-3 は step 1 で NaN (float32 幾何の限界と思われる、未追跡)。旧バイナリ (絶対ガード 1e-12) は
α=0.01 で d·S ≈ 2e-14 < 1e-12 となり必ず不一致になる (フックの都合で旧バイナリ直接実行は未実施)。

### 3.5 codex 検証レビュー (gpt-6-astra, reasoning high, 2026-09-08; ログ scratchpad/codex_review_sstopts.log) と対応

- **P1 α=1e-3 の NaN は幾何でなく Kato–Launder の `sqrt(S_sq*Om_sq)` の float32 中間積 overflow** (res_nan_1 で S_sq·Om_sq>FLT_MAX の 5768 ノード = roOmega 非有限ノードと一致) → `sqrt(S_sq)*sqrt(Om_sq)` に修正。修正後 α=1e-3 **PASS** (全量 1e-5 以下、ω 1.0e-3)。
- **P1 相似試験ツールの合格条件不足** (終了コード・要求 step 到達・非有限・保存量) → 修正 (res_nan があれば FAIL、res_{steps} 必須、roe/roK/roOmega も比較、rel が NaN なら FAIL)。既定 tol 2e-3。
- **P1 `sstIsotropicStress` は内部面のみで境界閉包 (入口・対称面・slip・壁関数壁) と軸対称 (τθθ, method 1 エネルギー幾何項) が未対応** → 既定 OFF 維持、残作業。
- **P1 `sstSigmaBlend` の cell 境界 ghost F1 未設定** → 拡散カーネルで ghost 側は内部側の F1 を使うよう修正。ラグは「前回の残差評価」(RK stage / dual-time subiter 単位) — 定常では Picard ラグとして許容、非定常は要確認 (残作業)。軸対称 method 1 の 1/y ソースは現行 F1 で常時ブレンド (評価時点の整合は残作業)。
- **P1 `sstEnergyKSource` の符号は正しいが、分離解法では有限刻みで ΣV(ρE_mean+ρk) が厳密には保存しない** (k 更新の対角減衰と E 側の源の不一致) → 既定 OFF 維持、時間離散収支の試験が要る (残作業)。
- **P2 `sstOmegaProdFromPk` の ν_t は `mu_t_eff` (1e-12 で a1 付きの代替値に切替) でなく closure の正本 `vis_turb` を使うべき** → 修正 (μt<1e-6μ では αρS_prod にフォールバック)。**「OFF=SST-2003 / ON=SST-1994」の説明は誤り**: NASA TMR は 2003 論文の αρS² を誤植とし、訂正式をリミッタ後の αP_k/ν_t としている → 訂正 (config コメント・solver-settings)。等方項 (dilatation 2) をリミッタ**後**に足す現行順序では最終 P_k が 10β*ρkω 以下とは限らない (残作業)。
- **P2 `sstNodeWallKPin` は定常 point-implicit では整合。RK (N/M 状態の再利用) と dual-time (残差 0 化後の BDF 項) では旧 restart の壁 k≠0 が時間履歴から再注入され得る** → 残作業 (非ゼロ壁 k restart のピン試験)。既定 ON には賛成。
- **P2 壁 ω の 2〜3 % スケール依存の主因は壁距離の評価点が双対 CV 重心 (axisCentroidShift) で、値位置 (ノード) と不一致** → `calcWallDistance_kdtree.cpp` を node ではノード座標で計算するよう修正。case/16 2D では第一層 wall_dist が最大 30 % 変わり (双対重心間 1.34 µm → 節点間 1.74 µm)、壁 p/p0 −0.1〜−0.3 %、ω ピン 2.0e11→1.2e11 (x=46)。平板 (case/26) は不変 (メッシュ均一で重心=節点)。
- **P1 A/B は回帰確認にはなるが既定値決定には不足** (全 run NOT CONVERGED、pmax/cf は 2 スナップショットで TRANSIENT-UNSETTLED、平板メッシュは AR 4828 で FAIL)。「差が小さいから OFF」は根拠にならない。追加検証案: KL overflow 試験 / 一定 ρk 境界付き領域で等方応力の偽残差 0 / F1=0,1,中間の拡散流束と ghost / 生産制限発動場での Pk/Pw 照合 / 一様乱流減衰の ΣV(ρE+ρk) 保存 / 非ゼロ壁 k restart のピン / 品質合格メッシュでの平板再検証 / node と cell 両方。
- **既定値の推奨**: `sstNodeWallKPin` ON 維持; `sstOmegaProdFromPk`・`sstSigmaBlend` は上記整理後 **ON を目指す**; `sstIsotropicStress`・`sstEnergyKSource` は境界/軸対称/時間離散収支の完備後に再判断。

## 4. 未決定事項

- 既定値: 3 は ON 確定。1/2/4/5 は A/B (§3.2, §3.3) と codex レビューを見て決める。4 は境界面未対応のまま既定 ON にしない。
- ω 壁ピンの 2〜3 % スケール依存 (§3.4) の原因。

## 5. 残作業表

| 優先 | 項目 | 状態 |
|---|---|---|
| 1 | codex (gpt-6-astra, reasoning high) に検証結果のレビューを受ける | 進行中 (2026-09-08) |
| 2 | 既定値の決定 (codex 推奨: OmegaProdFromPk / SigmaBlend を整理後 ON) と docs 反映 | 未 (ユーザ判断待ち) |
| 3 | `sstIsotropicStress` の境界面 (入口/対称面/slip/壁関数壁) と軸対称 (τθθ, method 1) 対応 | 未 |
| 4 | ~~§3.4 の ω 2〜3 % / α=1e-3 NaN の原因追跡~~ → 決着 (KL overflow / 壁距離の評価点, §3.5) | 済 |
| 5 | `sstEnergyKSource` の時間離散収支試験 (一様乱流減衰で ΣV(ρE+ρk)) | 未 |
| 6 | `sstNodeWallKPin` の RK / dual-time 経路での再注入試験 (非ゼロ壁 k restart) | 未 |
| 7 | codex 提案の追加検証 (等方応力の偽残差 0 試験, F1 中間値の拡散/ghost, 生産制限発動場の Pk/Pw, 品質合格メッシュでの平板, node/cell 両方) | 未 |
| 8 | 等方項 (dilatation 2) をリミッタ後に足す順序の見直し | 未 |

## 6. 変更ログ

- `2026-09-08` — 実装 (0fcde108)、2D/平板/3D A/B (§3.1–3.3)、相似試験ツール (§3.4)。**結論候補**: 3 (壁 k ピン) 既定 ON 維持、1/2/4/5 は影響 ≤1 % で既定 OFF のまま (SU2 との整合が要るケースで個別 ON)。codex レビューへ。
- `2026-09-08` — codex 検証レビュー (§3.5) → KL overflow 修正・ν_t を vis_turb に・cell ghost F1・node 壁距離をノード座標に・相似試験ツール強化。α=1e-3 PASS。既定値は ON=壁 k ピンのみのまま (残りはユーザ判断)。
