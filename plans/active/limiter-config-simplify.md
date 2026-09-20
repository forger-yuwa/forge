# リミッタまわりの設定を整理する (凍結を入れる前の判断)

## メタ

- **area**: `limiter`
- **status**: `in_progress`
- **related_docs**:
  - `methods/limiter.md` (現在仕様)
  - `methods/boundary.md`
- **related_plans**: `plans/active/convection-node-wall-reconstruction.md` (本件の発端。§4.20–§4.37)
- **created**: `2026-09-20`
- **owner**: `CFD Dev`

## 1. 目的

リミッタの調査 ([convection-node-wall-reconstruction](convection-node-wall-reconstruction.md)) で
**config キーが 8 個増え、plan が 134 KB になった**。ユーザから「**あまり複雑化したくない**」という
明示の要望が出ている。ここで一度止まり、**何を既定に畳み、何を消し、何を入れないか**を決める。

完了時の状態: リミッタの config 表面が最小になり、「② を既定にするか」の判断に必要な作業だけが残る。

## 2. スコープ

- **やる**: 既存キーの取捨 (既定に畳む / opt-in で残す / 消す) の判断と、その実施。
- **やる**: 「定常 2 次のプラトー対策 (ψ 凍結) を**そもそも入れるべきか**」の判断。
- **やらない**: ψ 凍結の実装そのもの (入れると決まったら別 plan)。
- **やる**: G4 (SU2 クロスチェック) を既定化ゲートから外す判断 → §4.5 で決着 (ユーザ承認 2026-09-20)。

## 3. 関連 docs と前提

- `methods/limiter.md` に現在仕様を同期済み (次元不整合・無次元化形・評価点・診断・フォールバック)。
- 実測の根拠はすべて [convection-node-wall-reconstruction](convection-node-wall-reconstruction.md) にある。

## 4. 設計方針

### 4.1 いま増えたキー (8 個)

| キー | 既定 | 実測で分かっていること |
| --- | --- | --- |
| `space.limiterMatchRecon` | 0 | 1 で `Ux`/`P` の近傍逸脱が 0 になる。**`limiterScaled>0` は 1 を要求する** |
| `space.limiterScaled` | 0 | 1 (無次元化) で `ro`/`Uy` の逸脱が 0。2 (比の形) は**棄却済み** |
| `space.venkatK` | 1.0 | **0.05 が正**。1.0 は Sod でも SERN でも悪い |
| `space.limiterRefLength` | 0 (自動) | 効果は小さい (基準値 ×2 で解は 2.1e-4) |
| `space.limiterRoRef`/`PRef`/`ARef` | 0 (自動) | 同上。restart ドリフトは restart 非忠実性に埋もれる |
| `space.limiterDiag` | 0 | 診断専用。生産では走らない |
| `space.badReconDiag` | 0 | 診断専用 |
| `space.badReconFallback` | 0 | **負の結果** (SERN ベース発散に効かない)。opt-in 残置 |

### 4.2 畳む案 (codex plan レビューを反映した確定版)

| キー | 処置 | 理由 |
| --- | --- | --- |
| `limiterScaled: 2` (比の形) | **削除** | §4.22 で棄却済み (残差床が 2〜5 倍)。残す理由が無い |
| `limiterScaled` | **公開値を 0 / 1 に畳む** | `0` = 現行の旧経路 (**G4 完了まで既定**)、`1` = 評価点・増分の整合を含む修正版 (Venkat は無次元化) |
| `limiterMatchRecon` | **廃止する。ただし「同義キーの削除」ではない** | `scaled=0, matchRecon=1` は**独立した有効経路** (`limiter_d.cu:373`: 評価点・増分だけ直して旧 Venkat 式を使う。Barth にも効く)。**明示的な機能打ち切り**として扱う |
| `venkatK` | **`limiterScaled: 1` のときの既定を 0.05 に** | 1.0 は Sod でも SERN でも悪い。**旧経路の K は `limiterFunctions_d.cuh:13` で `1.f` 固定**なので、既定を変えても旧経路は直らない (旧式・化学種側の固定 K は**今回の範囲外**) |
| `limiterRefLength` / `limiterRoRef` / `PRef` / `ARef` | **4 つとも残す** | 親 plan §4.27 は「既定化時の保存・継承は必須」に戻している。参照値の変化は $q_{\rm ref}^2(Kh_i/L_{\rm ref})^3$ = **作用素の変化**である。新規計算は自動でよいが、**G4・A/B・段階起動・継続では最初の値を既存キーへ固定して継承する** (新しいキーも checkpoint 機構も要らない) |
| `limiterDiag` / `badReconDiag` | **残す** | これが無いと今回の調査ができなかった |
| `badReconFallback` | **opt-in で残置** | [[optin-retained-is-not-rejected]] |

**⚠ 訂正**: §4.1 で `limiterRefLength` の感度根拠に使った「基準値 ×2 で 2.1e-4」は、
実際には `RoRef`/`PRef`/`ARef` を 2 倍にした試験 (`run_0503_w1f_ref_x2`) であり、
**`limiterRefLength` は指定していない**。長さ尺度の感度は未測定。

**削減は 8 → 7 キー** (比の形は値であってキーではない。`limiterMatchRecon` の廃止だけがキー削減)。

### 4.3 ψ 凍結は**入れない** (codex plan レビューの推奨を採用)

**判断**: **今回は ψ 凍結を実装しない。** プラトー run は「**未収束の準定常評価**」として扱い、
**目的量・局所場・保存収支・CFL/内部反復感度が事前に決めた誤差予算内にある場合だけ受理する**。

**根拠と、一般化しすぎない範囲**:

- リミッタを外すと収束する (`run_0056_pl_nolim` は `check_convergence` **PASS**、4.6〜4.7 桁) のは
  **この条件でリミッタが停滞に関与する強い証拠**だが、**作用素自体を変えている**ので
  「動的な ψ がある限り標準ケースは収束不能」「凍結が唯一の対策」までは示していない。
- 親 plan §4.37 の「CFL に比例」は**不正確**。`cfl_pseudo` 0.5 → 1 では残差は**下がる** (1.393e-3 → 1.184e-3)。
- **⑤ SERN の `run_0308` は ② の運用実績ではない** — その `solverConfig_main.yaml` は
  **`matchRecon=0, scaled=0`** である。力係数の安定 (`C_T_with_shear` 末尾平均 0.915413807、半幅 7.09e-8) は実在するが、
  これは旧経路での話。
- `sern_gates.py:281` 自身が「プラトーだけでは更新消失・クランプによる停滞と区別できず、
  無次元残差上限・保存収支が必要」と書いている。**現在の `GATES PASS` はそれらの検証完了を意味しない**。
- `case/08` の `check_quasisteady` は保存場 3 枚しかなく **TRANSIENT-UNSETTLED** (判定不足)。
  プラトー運用を名乗るなら**スナップショットを増やす**必要がある。

### 4.4 G4 はプラトーが残っていても実施できる (codex)

ただし条件つき: **両ソルバの比較対象量が `STEADY` で、反復誤差・格子誤差を含めて差を識別できること**。
親 plan が要求する **~0.025 %** の識別に対し、`check_quasisteady.py:543` の既定閾値は
**drift 5 % / 変動幅 10 %** である。**`ALL STEADY` の文字列だけで既定化は決められない。**
→ G4 では**誤差予算を先に固定し、必要なら閾値を絞る**。

### 4.5 G4 (SU2 クロスチェック) を既定化のゲートから外す (2026-09-20, ユーザ承認)

**判断**: **② を既定にする判定は「厳密解との比較」で行い、SU2 クロスチェックはゲートにしない。**

**理由 — SU2 比較の分解能が足りない**。リポジトリ内の実績:

| ケース | 比較量 | 一致度 |
| --- | --- | --- |
| `26.flat_plate_sst` 壁解像 | 壁 $C_f$/KS | forge node **0.943** vs SU2 **0.938–0.947** = **0.5 %** ← 実績最良 |
| `29.bell_vs_conical` 層流軸対称 | 中心線 Mach | 2〜5 % |
| `36.passive_pseudoshock_control` SST | peak $\mu_t/\mu$ | forge が SU2 より **20〜26 % 低い** |
| `48.flat_plate_cooled_m4` | SST 派生形どうし | ±1〜2 % |

**0.5 % はモデル差を殺した条件 (壁解像・層流近く) での最良値**で、乱流モデルが効くと 20 % 台になる。
`case/36` の 20〜26 % は **forge のバグではなく `dilatationCorrection` のモデル形式差**
([[sst-dilatation-bl-thickness]])。また `case/26` では **SU2 の壁関数版が運動量収支 21 % 破綻で棄却**されており、
**SU2 は「真値」ではなく独立した第二意見**である。

**判別したい差は 0.54 %** = 実績上の限界値そのもので、しかも SERN は SST + 多成分でモデル差が大きい側。
**このまま比較しても答えは出ない。**

**代わりに使う証拠** (すべて取得済み):

1. **厳密 Riemann 解との比較** (`case/05`, 3 解像度・CFL 固定, 親 plan §4.30):
   L1(ro) 収束次数 ① **0.76/1.00** vs ② **0.82/1.04**、衝撃後 ρ は全解像度で ② が厳密値 2.2240 に近い。
   **細かくするほど ② の優位が開く。**
2. **有界性** — ② は近傍逸脱 0、① は数千 (親 plan §4.32)。
3. **スケール不変性** — ① は ×1024 で実質 OFF、② はビット同一 (§4.21)。
4. **格子細分で悪化しない** (§4.30)。

**SERN の +0.54 % は「① が間違っていたぶん動いた」と解釈する。**
なお 0.54 % はノイズではない — 設計チェーン自身の再現性許容は `C_T_with_shear` 相対差 **2.4e-6** (R-b 検証) なので
**その 2000 倍**であり、「答えが変わった」と言える大きさである。だが**どちらが正しいかは厳密解で決める**。

**SU2 比較は降格**: 「やるなら後から、Euler か層流に落としてスキームの切り分けをする目的で」。
ゲートからは外す。

### 4.6 S1–S4 の実施結果 (2026-09-20)

**S1 ソルバ** (`input/solverConfig.cpp`, `cuda_forge/limiter_d.cu`, `limiterPeriodic_d.cuh`, `limiterFunctions_d.cuh`):

- `limiterMatchRecon` を**移行先を示して拒否**。メッセージに「**中間の `matchRecon: 1, scaled: 0` は機能打ち切りであって改名ではない**」と明記。
- `limiterScaled: 2` (比の形) を**理由つきで拒否**し、**カーネル分岐と `venkata_limiter_ratio` を削除**。
- `limiterScaled: 1` が評価点の一致を**内包**する (内部フラグを config 側で立てる)。
- **`venkatK` の既定を経路で分ける**: `scaled: 1` → **0.05**、`scaled: 0` → 1.0
  (旧経路の K は `limiterFunctions_d.cuh` で `1.f` 固定なので既定を変えても効かない)。

**S1 生成器** (`design/forge_design/evaluate/runner_sern.py`): `limiter_match_recon` を拒否、
`limiter_scaled` を 0/1 に限定、`venkat_k` の既定を **0.05 (scaled=1) / 1.0 (scaled=0)** に、
`limiterMatchRecon` の出力を停止。

**S3 段の同一性** (`tools/stage_manifest.py`): `HARD_PATTERNS` に **`limiterScaled` と `venkatK`** を追加。
リミッタの式が変わった段を同一区間として連結しなくなった。

**S2/S4 検証**:

| 確認 | 結果 |
| --- | --- |
| 廃止キーの拒否 | `limiterMatchRecon: 1` → 移行先つきで拒否 ✅ / `limiterScaled: 2` → 理由つきで拒否 ✅ |
| 新経路の既定 | `limiterScaled: 1` + `venkatK` 省略 → ログに **`K=0.05`** ✅ |
| **既定パスの退行 (対照)** | `case/05` **1.06e-4** (ノイズ床 1.51e-4) / `case/44` **9.28e-6** (床 3.38e-4) / `case/09` **0.000e+00 (ビット同一)** → **退行なし** |
| 設計チェーン | 生成 config が `space: {..., limiterScaled: 0, venkatK: 1.0}` になり forge が起動する ✅ |

**過去の計測 run は書き換えていない** ([`config-key-pruning.md`](../accepted/config-key-pruning.md) の方針)。
旧キーを持つ run は**起動時に拒否され、移行先が示される**。再実行したい場合だけコピーを移行すること。

**キー数: 8 → 7** (`limiterMatchRecon` の廃止。比の形は値なのでキー数に入らない)。

## 5. 実装ステップ

### 5.1 残作業

| 項目 | 内容 |
| --- | --- |
| S1 | ~~§4.2 の畳み込み~~ **完了 (2026-09-20) → §4.6**。**ソルバだけでなく設定生成器も同時に直す** — `design/forge_design/evaluate/runner_sern.py:210` は `limiter_match_recon`/`limiter_scaled`/`venkat_k` を独自に読み **`venkat_k` の既定を 1.0** にしており、:257 は廃止予定の `limiterMatchRecon` を**常に生成**する。ソルバだけ直すと起動不能か、新既定 0.05 が効かない (codex M2) |
| S2 | ~~廃止キー・廃止値の拒否~~ **完了 → §4.6**。 ([`config-key-pruning.md`](../accepted/config-key-pruning.md) の方針)。**過去の計測 run を黙って書き換えない** — 原設定を保存し、再実行用コピーだけ移行する。旧 `scaled=1` で `venkatK` 省略の run は再現用に `1.0` を明示する |
| S3 | **段階の同一性判定にリミッタ設定を含める** (codex M5)。`stage_manifest.py` が今回のキーを見ていないので、リミッタ設定が変わった段を同一区間として連結してしまう |
| S4 | ~~標準ケース回帰~~ **完了 → §4.6** (既定パスは 05/09/44 でノイズ床以下、09 はビット同一)。**ノイズ床比較だけでは意図した変更と退行を区別できない** (codex M6) ので、**既定を変えない対照**と**変えた対照**を両方走らせ、差が意図した経路に由来することを示す |
| S5 | ~~ψ 凍結の要否判断~~ **決着: 入れない** (§4.3)。プラトー運用の受理条件を `procedures/` に書く |

## 6. 検証

- S1 は既定値変更を伴うので、標準ケースで**変更前後の場をノイズ床基準で比較**する。
- S1/S2 は既定値変更を伴うので、**変更前後の場を「既定を変えない対照」と併せて比較**する (ノイズ床だけでは足りない)。
- プラトー運用を名乗る run は **`check_quasisteady` が判定できるだけのスナップショット数**を残す (3 枚では TRANSIENT-UNSETTLED)。

### 6.1 レビュー記録 (codex)

| stage | 日付 | 記録 | 判定 | 採否 |
| --- | --- | --- | --- | --- |
| plan | 2026-09-20 | [2026-09-20-limiter-config-simplify-plan.md](../../notes/reviews/2026-09-20-limiter-config-simplify-plan.md) | GO-with-changes, M6/m2 | **全件採用 → §4.2–§4.4, §5.1**。**M1** → `matchRecon=1, scaled=0` は独立した有効経路 (`limiter_d.cu:373`、Barth にも効く) なので廃止は**機能打ち切り**と明記。旧経路の K は `limiterFunctions_d.cuh:13` で `1.f` 固定なので既定変更では直らない。**M2** → 設定生成器 (`runner_sern.py:210/257`) を S1 に追加。**M3** → **ψ 凍結は入れない**。ただし「凍結が唯一の対策」も「プラトーなら運用可」も一般化しないよう限定。`run_0308` は `matchRecon=0, scaled=0` で ② の運用実績ではないと明記。**M4** → 参照 4 キーは残し G4/A-B/段階起動では固定継承。§4.1 の `limiterRefLength` 感度の誤帰属を訂正。**M5** → 段の同一性判定にリミッタ設定を含める (S3)。**M6** → 対照つき比較に変更 (S4)。**m1/m2** → 削減数を 8→7 に訂正、責務と現在仕様の同期 |

## 7. 影響範囲

`solver_density_cuda/input/solverConfig.{hpp,cpp}`, `cuda_forge/limiter_d.cu`,
`cuda_forge/limiterPeriodic_d.cuh`, `methods/limiter.md`, 既存 run の `solverConfig.yaml`、
**`design/forge_design/evaluate/runner_sern.py` / `runner_sern3d.py` (設定生成器)**、
**問題 YAML**、**`solver_density_cuda/tools/stage_manifest.py`** (段の同一性判定)。

## 8. 完了条件

config 表面が最小化され、ψ 凍結の要否が決まり、標準ケースで退行が無いこと。

## 9. 変更ログ

- 2026-09-20: 起票。
