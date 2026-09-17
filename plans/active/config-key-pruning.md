---
title: solverConfig キーの整理 (採用があり得ないパラメータの削除)
related_docs:
  - procedures/recommended-settings.md
  - procedures/solver-settings.md
status: in_progress
updated: 2026-09-17
---

- **status**: `in_progress`
- **起票**: 2026-09-17 (ユーザ指示: 「パラメータがとても多いように思っていまして、採用があり得ないパラメータは消していってほしい」)

## 1. 目的

`solverConfig.yaml` のキーが 155 個まで増え、**どれが生産設定でどれが一時的な A/B スイッチか区別できない**状態になっている。
新規 config を組むときの誤用 (旧既定・非推奨の組合せ) の温床でもある。**採用があり得ないキーを削除**し、残すキーは
「生産で選ぶ」「診断・A/B で要る」のどちらかに分類して文書に位置づけを書く。

## 2. スコープ

- 対象: `solver_density_cuda/input/solverConfig.{hpp,cpp}` が読むキーと、それを参照する solver コード・文書・case config。
- 対象外: `bcondConfig.yaml` の境界キー、設計チェーン (`problem_*.yaml`)、メッシュ変換器のキー。別途。

## 3. 現状 (棚卸し, 2026-09-17)

`tools/config_key_inventory.py` による機械集計 (case 配下の run config 2503 本と `procedures/` `methods/` を突き合わせ)。

| 指標 | 数 |
| --- | --- |
| 読み込むキー | 155 |
| どの run でも設定されていない | 25 |
| 設定されていても既定値のみ (実質 OFF のまま) | 3 |
| `recommended-settings.md` に出てこない | 78 |

## 4. 設計方針 (削除の基準)

**削除する**のは次のいずれかに当てはまり、かつ**後継または既定で置き換えられる**もの。

- **A. 実測で効果が無い / 有害と判定済み**: 計測結果が plan・notes にあり、採用しない結論が出ているスイッチ。
- **B. 後継に置換された旧経路で、A/B 期間が終わったもの**: 回帰が新既定で揃っており、旧挙動の再現要求が無いもの。
- **C. 決定で不採用になったモード**: ユーザ決定・codex レビューで「使わない」と決めたモード値やフラグ。
- **D. 一度も使われず、既定から動かす根拠も無い内部チューニング定数**: 定数としてコードに埋める。

**残す**のは次のいずれか。位置づけを `solver-settings.md` に 1 行で書く。

- **P. 生産で選ぶ**: 解析種別ごとに値が変わるもの (`convMethod`, `cfl_pseudo`, `turbulence.model` …)。
- **V. 診断・A/B で要る**: 切り分けに使う実績があり、今後も使うもの (`detectNaN`, `extraFields` …)。

**削除の手順** (キーごと): コード分岐を削除 → 既定側の挙動だけを残す → 旧キーが config にあったら**起動時エラー**で気づかせる
(黙って無視しない) → `solver-settings.md` の該当節を削除し `recommended-settings.md` §9 (旧設定) に 1 行残す。

## 5. 実装ステップ

1. 棚卸しツールを `solver_density_cuda/tools/config_key_inventory.py` として恒久化 (再実行できる形に)。
2. 全キーを A–D / P / V に分類した表を §5.1 に作る (根拠 1 行つき)。
3. codex plan レビュー。分類の妥当性と「消してはいけないもの」の指摘を受ける。
4. 分類 A–D を削除。1 キー 1 コミットではなく、まとまり (line 系, SST 分離型, 旧凝縮経路 …) 単位。
5. 回帰: §6。
6. `solver-settings.md` / `recommended-settings.md` / skill `forge-config` を同期。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | 棚卸しツールの恒久化 | `config_key_inventory.py` (キー定義・既定値・run での使用実績・文書での言及を突き合わせ、分類表の素を出す) |
| 2 | ~~全 155 キーの分類表~~ | 済 (2026-09-17): §5.2。削除候補 31 (A 14 / B 12 相当 / C 5 / D 10)、保留 7、残す 117。付随して文書バグ 5 件 (§5.3) |
| 3 | codex plan レビュー | 分類表ができた時点で `--stage plan` |
| 4 | 削除の実装 | まとまり単位。旧キーは起動時エラーにする |
| 5 | 回帰検証 | §6 |
| 6 | 文書・skill の同期 | `solver-settings.md` の節削除、`recommended-settings.md` §9 への 1 行、`check_solver_config.py` に旧キー検出を追加 |

### 5.2 分類 (2026-09-17 調査, 3 領域を分担して根拠づけ)

分類基準は §4。**根拠は全て `ファイル:行` か run の実測**。全 155 キーを見て、削除候補は 31、保留 7、残す 117。

#### 削除候補 A: 実測で効果なし / 有害と判定済み

| キー | 既定 | 非既定 run | 根拠 |
| --- | --- | --- | --- |
| `implicitRelaxSST` | −1 | 3 | 独立 3 件の A/B が全て無効 (line-implicit v2 plan「ω 収縮はむしろ悪化」、軸対称 plan「プラトー不変」、node ノズル plan「発散時期を変えない」) |
| `updateGuardAlpha` | 0.0 | 0 | `time_integration-update-positivity-guard.md`「guard 単独 α=0.5: cfl 4/8/16/32 全て発散」「負の結果で done」 |
| `lineViscCoupling` | 0 | 30 (診断のみ) | line-implicit v2 plan「粘性結合・dt 割引はこのケースでは僅差 (roe +0.02 桁)」「壁法線律速は音響 (λ_visc/λ_ac≈0.02)」 |
| `lineViscousDtRelief` | 0.0 | 29 (診断のみ) | 同上 |
| `lineDtWallRelief` | 0 | 1 (診断) | 同 plan「結果は分岐③: 発散 (step ~80–100 で ro 非有限)」。ヘッダも「(診断)」 |
| `blockDPLURDiagCache` | 0 | 0 | `performance-3d-node-sst-speedup.md`「実装したが不採用: 44.0→46.5 ms/step と遅化」 |
| `blockDPLURDqPack` | 0 | 0 | 同 plan「不採用で確定」、`architecture/performance.md`「+0.7〜+2.7 ms/step の逆効果」 |
| `lowMachThornber` | 0 | 2 | `convection/implementation.md`「検証結果 (負)。無効〜僅かに悪化…根治用途では使わない」 |
| `multispeciesRhoYCommonLimiter` | 0 | 3 | `convection-multispecies-contact-pressure.md`「依然 S2 より明確に悪い」「cfl4 で C も発散 step 550」。S3 の生産化で救済目的自体が消滅 |
| `primPack` | 0 | 0 | `performance-3d-node-sst-speedup.md`「不採用で確定…差なし〜微増」 |
| `condGyarmathyC` | 3.18 | 4 | case/16「この系では Gyarmathy 標準係数 3.18 がほぼ最適」(C=0/1.59/6.36/12.72 で悪化) → 定数化 |
| `condTwoTemp` | 0 | 2 | `condensation-nonequilibrium.md`「希薄水/N2 では影響 <1 % で一温度が既定。追加開発は不要」。float 経路も無効化する副作用つき |
| `condEvapKelvin` | 0 | 0 | `methods/condensation.md`「正帰還で離散的に Q1→0 崩壊」「質量収支には効かない」 |
| `condN2LiquidCp` | 2000 | 2 | 掃引完了 (1500/2000/2500 で onset ±0.6〜0.7 K の同符号) → 定数化 |

#### 削除候補 B: 後継に置換された旧経路で A/B 期間が終わったもの

| キー | 既定 | 非既定 run | 後継 / 根拠 |
| --- | --- | --- | --- |
| `sstIsotropicStress` | 0 | 3 | `sstEnergyIncludesK`。`turbulence-sst-energy-includes-k.md`「既定 1 に切替、`sstEnergyKSource`/`sstIsotropicStress` を削除」 |
| `sstEnergyKSource` | 0 | 3 | 同上。`solverConfig.cpp` が `sstEnergyIncludesK=1` 時に強制 0 にしており後継関係がコードにも出ている |
| `sstNodeWallKPin` | 1 | 2 | 0 は「node 境界半割面拡散 skip で壁 k=0 が効いていなかった」旧バグ挙動。A/B は 2 run で終了 |
| `nodeWallDirichlet` | 1 | 4 | 0 は壁速度が立たない旧挙動。切り分け 4 run のみ |
| `nodeWallStressEdgeKernel` | 1 | 2 | 0 は旧 twall。出力専用で場は不変、A/B 2 run で終了 |
| `turbulentSchmidt` | (別名) | 0 | `physProp.Sc_t` の同義キー。どの run も使っていない |
| `physProp.isAxisymmetric` / `physProp.axisymMethod` | 0 | — | 正本は `mesh` ブロック。`solverConfig.cpp` が既に deprecated 警告を出している |
| `gpu` | 1 | 15 (旧 case) | `methods/boundary.md`「CPU 経路 (`gpu == 0`) は未対応 (`applyBconds` 冒頭で exit)」= 事実上 1 値 |
| `condN2LatentLowT` | 1 | 3 | `condensation-air.md` §9 表で A/B 決着 (旧 onset 2.367 in → 新 2.112 in、理論線 +1.0→+2.1 K で合格) |
| `condN2PsatLowT` | 1 | 2 | 同上 (run_0016 行で決着、onset 不変・g_exit 0.081→0.059) |
| `condEqRelax` / `condEqDTmax` / `condEqDgMax` | — | 0 | `condEquilibrium==1` 分岐でしか読まれない。mode 1 は mode 2 (EOS 拘束形) に置換済み (固定点同一・θ 遅れ無し) |

#### 削除候補 C: 決定で不採用になったモード値 (キーは残し、値の受理をやめる)

| 対象 | 根拠 |
| --- | --- |
| `bndFirstOrder` (キーごと) | 使用禁止。削除計画 [architecture-bndfirstorder-removal](architecture-bndfirstorder-removal.md) が既存 → 本 plan からはリンクのみ |
| `condEquilibrium: 1` (緩和形) | `condensation-equilibrium-eos.md`「mode 2 は onset 以降 全軸で S=1.0000 (緩和形は S≤1.18 が ~2 r_t 続く)、下流は |ΔM|≤1e-3 で一致」 |
| `sstThermalWallFunction: 2` | `turbulence/implementation.md`「動的試験で壁ノード単調冷却が止まらず現形は未採用」 |
| `timeIntegration: 1` | `3` と完全に同一分岐 (`timeIntegration_d.cu`, `scalarTransport_d.cu` とも `== 1 or == 3`)。14 run のみ |
| `mesh.axisymMethod: 1` | 「実装済みだが implicit 深収束未達」「出口角で発散 — 既知」 |

#### 削除候補 D: 一度も使われず既定から動かす根拠も無い内部定数 → コードに埋める

`C_DES_kw` (0.78) / `C_DES_ke` (0.61) / `wmlesNewtonTol` (1e-6) / `wmlesNewtonMaxIt` (20) / `wmlesPrt` (0.9) /
`gradLSQDegenThresh` (1e-2) / `condEvapRmin` (1e-9) / `passiveFctTolAbs` (1e-30) /
`physProp.chemistry.tMaxReaction` (6000) / `freezeBelowT` (0) の 10 キー。いずれも全ワークツリーの run config で使用ゼロ、掃引の実測も無い。

#### 保留 (判断が要る)

| キー | 保留の理由 |
| --- | --- |
| `ducrosLimiter` | 単純削除できない。`enable==0` で `ducros` 場を 0 に潰すため、KEEP の blend `max(ducros, 1−ψ)` から項が永久に消える。[turbulence-iddes-sst](turbulence-iddes-sst.md)「中期で改良 Ducros を別途復活」と要突き合わせ |
| `condLimiterMode` | 旧経路 0 に RK 陽解法・`passiveScalarScheme 0` の dual-time が**自動降格で依存**している。followups F-cf9 の前提 (RK の更新クランプ試験) が未了 |
| `nodeOmegaWfDirichlet` | `turbulence-sst-su2-taw-coupling.md` が「剛性対策オプションとして残置」と**明示的に決定済み**。削除は決定の見直しが要る |
| `sstOmegaProdFromPk` / `sstSigmaBlend` | `recommended-settings.md` §9 では旧設定扱いだが、case/44 の active な SU2 クロスチェック run が 0 を使用中 = A/B 期間が終わっていない |
| `keepDissCbCoeff` / `keepDissCbEps` | 「本番格子で off と同一 or NaN」「動機だった丘頂鋸歯は抽出アーチファクトで撤回」だが、plan が `draft` のまま |
| `isCompressible: 0` (SMAC 非圧縮) | 使用は旧 4 run のみだが `methods/poisson.md`「コードベースに残っており利用できる」。経路ごと廃止するかはユーザ判断 |
| `passiveFctSweeps` / `passiveFctTol` | 2026-09-16 新設で「使用 0」は年齢の反映。凝縮 dual-time の生産期に `[passiveFct] WARNING` が出ないことを確認してから定数化 |

#### 残すキー (P: 生産で選ぶ / V: 診断・A/B で要る)

上記以外の 117 キー。`solver-settings.md` に位置づけ 1 行が無いものは追記する (§5.1 #6)。

### 5.3 付随して見つかった不具合 (削除とは別に対処)

| # | 内容 | 状態 |
| --- | --- | --- |
| a | `recommended-settings.md` の SST レシピが**存在しないキー** `kInf` / `omegaInf` を推奨していた。実キーは `kInit` / `omegaInit` で、書いても黙って無視され k=ω=0 の初期値になる (node SST が step 0 で壊れる既知の原因) | **修正済 (2026-09-17)** |
| b | 同 §3 の `speciesPrecondDt: 1 (既定)` は化学ブランチ限定で main / sern に無い | **注記済 (2026-09-17)** |
| c | `thermCondMethod` (297 run が 1 を使用) と `prandtlLam` が `procedures/` `methods/` に一度も出てこない | §5.1 #6 で `solver-settings.md` に追記 |
| d | `initial` は solver 実行時には読むだけで使われない (`setInitial` の呼び出しは変換器のみ) のに全 run config で必須キー | 変換時キーへ移すか、必須をやめる。§5.1 #4 に含める |
| e | 存在しないキーを書いても黙って無視される (a の原因)。`check_solver_config.py` に「未知キーの検出」を足すべき | §5.1 #6 に追加 |

## 6. 検証

1. **ビット不変**: 削除したキーを**既定値で**書いていた run (case/16 `run_0476`, case/44 `run_0312`, case/09 `run_0156` 等) を
   削除前後のバイナリで再実行し、保存量がビット一致すること (既定の挙動を変えていないことの証明)。
2. **起動時エラー**: 削除したキーを書いた config が起動時に落ちること (黙って無視されない) を単体試験で。
3. **config lint**: `check_solver_config.py` に旧キー検出を足し、case 配下の run config を走査して残存を洗い出す。
4. **ビルドと単体試験**: 既存の単体試験が全て PASS。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| (未実施) | | | | |

## 7. 影響範囲

solver の config 読込と分岐、`procedures/` の設定文書、skill `forge-config`、既存 run config (削除キーを書いているものは
再実行時に起動エラーになる → case README に注記)。

## 8. 完了条件

- [ ] 分類表 (§5.1 #2) が全キーを覆い、codex plan レビューを通っている
- [ ] A–D の削除と旧キーの起動時エラー化
- [ ] §6 の 1–4 を満たす
- [ ] `solver-settings.md` / `recommended-settings.md` / skill の同期

## 9. 変更ログ

- `2026-09-17` — 起票 (ユーザ指示)。棚卸し (155 キー / 未使用 22 / 既定のみ 3 / recommended 未記載 78; 全 9 worktree の run config 5163 本を対象)。
- `2026-09-17` — 全キーの分類を根拠つきで作成 (§5.2)。付随して `recommended-settings.md` の存在しないキー 2 件を修正・注記 (§5.3 a/b)。

## 10. 未確定事項

- 「診断・A/B で要る」の線引き: 過去 1 回しか使っていないスイッチを残すか。**案**: 切り分けの再現手順が plan/notes に
  書かれているものだけ残し、それ以外は削除して必要になったら復活させる (git にある)。
- 既定値を変えるか: `lineKFreeze` のように「dual-time では常用推奨」なのに既定 0 のものは、削除ではなく既定変更が適切。
