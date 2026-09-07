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

## 4. 未決定事項

- 既定値: 3 は ON 確定。1/2/4/5 は A/B (§3.2, §3.3) と codex レビューを見て決める。4 は境界面未対応のまま既定 ON にしない。
- ω 壁ピンの 2〜3 % スケール依存 (§3.4) の原因。

## 5. 残作業表

| 優先 | 項目 | 状態 |
|---|---|---|
| 1 | codex (gpt-6-astra, reasoning high) に検証結果のレビューを受ける | 進行中 (2026-09-08) |
| 2 | 既定値の決定と docs (methods/turbulence) 反映 | 未 |
| 3 | `sstIsotropicStress` の境界面対応 | 未 |
| 4 | §3.4 の ω 2〜3 % / α=1e-3 NaN の原因追跡 | 未 |

## 6. 変更ログ

- `2026-09-08` — 実装 (0fcde108)、2D/平板/3D A/B (§3.1–3.3)、相似試験ツール (§3.4)。**結論候補**: 3 (壁 k ピン) 既定 ON 維持、1/2/4/5 は影響 ≤1 % で既定 OFF のまま (SU2 との整合が要るケースで個別 ON)。codex レビューへ。
