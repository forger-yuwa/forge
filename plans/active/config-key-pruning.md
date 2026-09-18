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
- **対象に含む (codex plan M4)**: 同じ `solverConfig::read()` を使う**メッシュ変換器** (`mesh/convertGmshToForge.cpp`; `cfg.gpu` を確保処理に渡す)
  と、`solverConfig.yaml` を**生成する側** (`design/forge_design/evaluate/runner*.py` の 3 本が `gpu` / `nodeWallDirichlet` を書く)。
  生成器を直さずに起動時拒否だけ入れると、新規の設計ケースが起動できなくなる。
- 対象外: `bcondConfig.yaml` の境界キー、設計チェーンの問題 YAML (`problem_*.yaml`) の仕様。

## 3. 現状 (棚卸し, 2026-09-17)

`tools/config_key_inventory.py` による機械集計 (case 配下の run config 2503 本と `procedures/` `methods/` を突き合わせ)。

**改訂 (2026-09-18, codex plan M3)**: 初版の抽出はヘルパ経由の呼び出ししか見ておらず `config["mesh"][...]` 形を取りこぼし、
同名の末端キー (例 `control` が `time.deltaT` と `time.last` の両方にある) を 1 つに統合していた。**完全修飾パスで再抽出**した。

| 指標 | 初版 (名称ベース) | 改訂 (完全修飾パス) |
| --- | --- | --- |
| 読み込むキー | 155 | **188** |
| どの run でも設定されていない | 25 | 26 |
| 設定されていても既定値のみ | 3 | 3 |
| `recommended-settings.md` に出てこない | 78 | 89 |

集計対象は全 9 worktree の run config 5163 本 (本ツリーだけだと `forge-cond` / `forge-chem` の実使用を見落とす)。

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
| 1 | 棚卸しツールの恒久化 | `config_key_inventory.py`。**再オープン (2026-09-18)**: 使用数をキー名ベースでなく**完全修飾パスで PyYAML 集計**に直す (誤判定の原因; §5.3 j) |
| 2 | 全キーの分類表 (live な値キー 167 パス) | 初版 (155 名称) は素案。**codex plan M3 / plan-2 M3 で抽出の取りこぼしと誤抽出が出たので再オープン**。修正後の棚卸しは **値キー 167 / 節 10 / 拒否専用 12 / 起動時拒否 8**、run config 4030 本。確定した削除は §5.2 の 2 パス + 定数化 5 件 + 撤去 1 件のみ。残りは**この一覧で**分類をやり直す |
| 3 | codex plan レビュー | 分類表ができた時点で `--stage plan` |
| 4 | 削除の実装 (第 1 陣) | ~~確定 10 パス + 定数化 7 件~~ **済 (2026-09-18)**: config 読みとメンバを削除し、呼び出し側は既定値を直接渡す。旧キーは**起動時エラー** (どこへ移ったかを言う)。`check_solver_config.py` に未知キー検出を追加。第 2 陣は分類やり直し後 |
| 5 | 回帰検証 | **済 (2026-09-18 やり直し)**: §6.2' の事前確定基準で 4 経路 (SST 定常陰解法 / DDES / line-implicit + DES 診断 / WMLES) + 拒否 8 件 + 受理 9 件 + opt-in 3 件。**全 PASS** (§6.4')。初回の 3 run は全量が反復幅以内でなく (`P` 2.125 vs 2.0) 収束も準定常も未達だったので破棄 |
| 7 | 棚卸し・検査ツールの修正 | **済 (2026-09-18, codex plan-2 M2/M3/M4/M5)**: 入れ子 config の探索、コメント除去、節・拒否専用参照の分離、必須キーの既定値なし、別名解決 (代入先メンバ名で引く)、絶対許容差の撤去、記載/非既定の分離、読み取り失敗の報告。`check_solver_config.py` の未知キー検出を**完全修飾パス**に置換 (誤配置を WARN、起動時拒否を FAIL) |
| 8 | 文書の同期 | **済 (2026-09-18, codex plan-2 M6)**: `recommended-settings.md` §9.1 を最終対象に合わせ、残置した 9 件を明記 |
| 6 | 文書・skill の同期 | ~~第 1 陣ぶん~~ **済 (2026-09-18)**: `recommended-settings.md` §9.1 に削除キーの表 (理由と移行先)、`solver-settings.md` の該当記述を削除。第 2 陣は分類後 |

### 5.2 分類 (2026-09-18 改訂, codex plan レビュー反映)

初版 (2026-09-17) は 3 領域を分担して根拠づけたが、codex plan レビューで **4 件の誤分類と抽出の取りこぼし**が出たので改訂した。
改訂の方針は「**根拠が確定したものだけ削除し、限定的な実測から機能廃止へ飛躍しない**」。

**分類のやり直し (§5.1 #2 を再オープン)**: 完全修飾パス 188 件を単位に、各項目へ
`実効既定値 / 読む場所 / 分類 / 根拠 / 移行先 / 検証` を持たせて作り直す。初版の 155 名称ベースの表は下の「確定分」以外は素案として扱う。

#### 確定: 今回削除するもの

根拠が「全用途で不要」まで言えるものだけ。**いずれも run 使用ゼロか診断 run のみ**で、accepted plan に不採用が明記されている。

| パス | 既定 | 根拠 |
| --- | --- | --- |
| `time.deltaT.lineDtWallRelief` | 0 | line-implicit v2 plan「結果は分岐③: 発散 (step ~80–100 で ro 非有限)」。ヘッダも「(診断)」。診断 1 run のみ |
| `time.deltaT.implicitRelaxSST` | −1 | 独立 3 件の A/B が全て無効 (ω 収縮は悪化 / プラトー不変 / 発散時期不変)。**注意 (M9)**: −1 は `implicitRelax` の継承指定なので、削除は「SST も `implicitRelax` に従う」への統合であって定数置換ではない |

**削除から戻したもの (codex plan-2 M1, 2026-09-18)**: `time.deltaT.blockDPLURDiagCache` /
`time.deltaT.blockDPLURDqPack` / `mesh.primPack` の 3 件。いずれも元 plan
([performance-3d-node-sst-speedup](../accepted/performance-3d-node-sst-speedup.md):160, `solverConfig.hpp`) の処置が
「**不採用で確定、opt-in 残置**」であり、「既定に採用しない」と「設定機能を廃止する」を取り違えていた
(`updateGuardAlpha` 等で 2 度自己訂正したのと同じ誤り。**これで 3 度目**)。さらに棚卸しの「使用 0」も誤りで、
入れ子探索を入れると `forge-perf` の nested run が実際に 1 を書いている (それぞれ 1 / 4 / 4 run が非既定)。
**廃止するなら、再現用途を終了する理由と影響する run を明記して残置決定を置き換える判断が要る**。

**未使用の内部定数 (D)** は、**「定数として固定する」判断を明記したうえで**埋め込む: `turbulence.C_DES_kw` / `C_DES_ke` /
`wmlesNewtonTol` / `wmlesNewtonMaxIt` / `mesh.gradLSQDegenThresh` の 5 件
(全 worktree で使用 0、掃引の実測も無し)。`turbulence.wmlesPrt` は**定数化ではなく撤去** (codex plan-2 m8):
Kader 原式への修正 (`wallLaw_d.cuh:130`) で壁法則が Pr_t を使わなくなっており、カーネル引数 `Prt` を誰も読んでいない。
メンバ・引数ごと削除し、`turbulentPrandtl` への移送案内も撤回した (壁法則とは別物)。`condGyarmathyC` と `condN2LiquidCp` は掃引済みだが、**感度試験機能を捨てる**判断になるので
§5.4 の保留へ移した。

#### 保留に移したもの (初版で削除候補にしていたが、根拠が足りない)

| 対象 | 初版の誤り | 保留解除の条件 |
| --- | --- | --- |
| `time.timeIntegration: 1` | 「3 と同一分岐」は誤り。**1 段 Euler と 3 段 TVD RK で段数・係数が別** (`solverConfig.cpp:907/919`) | Euler 陽解法を機能として廃止するかの判断 (次数・安定性・計算量の変更を伴う) |
| `turbulence.sstIsotropicStress` / `sstEnergyKSource` | 引用先の**最終決定と逆**。accepted plan §4.0 は `sstEnergyIncludesK` 既定 0 + 分離型の個別利用可、残作業表も「撤去せず現状維持」 | 分離型を廃止する新しい判断 (境界未完備・離散保存せずを理由に) を別途行う |
| `condensation.condEquilibrium: 1` と `condEqRelax`/`condEqDTmax`/`condEqDgMax` | 「固定点同一」は一般には誤り。新しい accepted plan が**緩和形の固定点は Δτ 依存**と明記、mode 2 は単一凝縮種限定 | 旧モデルの受付範囲を明示したうえでの廃止判断 |
| `mesh.axisymMethod: 1` | 「不採用決定済み」ではない。active plan が **opt-in 保持・線形化改善後に再評価**と書いている | 当該 plan の再評価が終わること |
| `time.deltaT.lineViscCoupling` / `lineViscousDtRelief` | 実測は「**このケースでは**僅差」で、粘性律速になる条件も同 plan に書かれている | 粘性律速側 (低 Re・微小 Δn) での確認、または適用範囲を限定した廃止判断 |
| `condensation.condTwoTemp` | 実測は**希薄水/N2 限定**。全用途で不要の証拠ではない | 液滴負荷の高い条件での確認 |
| `physProp.chemistry.tMaxReaction` / `freezeBelowT` | 内部反復定数ではなく、**反応源の温度評価と停止条件を変える** (`chemistry_d.cu:86`)。使用実績ゼロは廃止の証明にならない | 化学ブランチ側の判断と同時に |
| `condensation.condGyarmathyC` / `condN2LiquidCp` | 掃引で既定が最良と分かったことは、**感度試験機能を捨てる**理由にはならない | 標準モデルとして固定する判断を明示できたら |

| `time.deltaT.updateGuardAlpha` | **元 plan の「処置」が「opt-in のまま残置 (無害・発散遅延の診断的価値はある)」**。負の結果は「CFL 上限を上げない」という限定 | LHS 整合化で cfl 上限の真因が解消すること |
| `time.deltaT.lowMachThornber` | 元 plan は「**この症状 (ノズル limit cycle) に無効**。機能は opt-in で残置」。低マッハ LES/解像用途は未検証 | 低マッハ解像用途での確認、または機能として持たない決定 |
| `time.deltaT.multispeciesRhoYCommonLimiter` | 元 plan が「**診断オプションとして残置**」と決定済み。実測 (case/28 TP, cfl 1/2/4) は確かだが、削除は残置決定を覆す新しい判断になる | その判断を明示的に行うこと |
| `time.deltaT.passiveFctTolAbs` | `passiveFctSweeps` / `passiveFctTol` と**同一 plan・同一日・同一の受入判定**。3 つとも同じ扱いにする | 3 キー一括で判断 (凝縮 dual-time 生産期の確認後) |

| `turbulence.turbulentSchmidt` | **161 run が使用中** (chem ブランチ: case/48 が 0.7 ×124、case/47 が **0.5** ×37)。`physProp.Sc_t` の既定 0.7 と違う値があるので、単なる別名撤去ではなく**設定値の移送**が要る | chem ブランチ側の run config を `physProp.Sc_t` に書き換えること |
| `physProp.isAxisymmetric` / `physProp.axisymMethod` | **1064 run が使用中** (値 1 が 593)。deprecated 読み (警告つき) は生きており、拒否すると既存 run の再実行が全部落ちる | 既存 run の移送、または警告のまま残す判断 |

#### 保留 (初版から継続)

| キー | 保留の理由 | 保留解除の条件 |
| --- | --- | --- |
| `ducrosLimiter` | `enable==0` で `ducros` 場を 0 に潰すため、KEEP の blend `max(ducros, 1−ψ)` から項が永久に消える | [turbulence-iddes-sst](turbulence-iddes-sst.md)「中期で改良 Ducros を別途復活」の方針が決着すること |
| `condensation.condLimiterMode` | 旧経路 0 に RK 陽解法・`passiveScalarScheme 0` の dual-time が**自動降格で依存**している | followups F-cf9 の前提 (RK の更新クランプ試験) が済むこと |
| `turbulence.nodeOmegaWfDirichlet` | accepted plan が「剛性対策オプションとして残置」と**明示的に決定済み** | その決定を見直すこと |
| `turbulence.sstOmegaProdFromPk` / `sstSigmaBlend` | case/44 の active な SU2 クロスチェック run が 0 を使用中 = A/B 期間が終わっていない | 当該 run の役割が終わること |
| `keepDissCbCoeff` / `keepDissCbEps` | 「本番格子で off と同一 or NaN」「動機だった鋸歯は抽出アーチファクトで撤回」だが、plan が `draft` | 当該 plan の決着 |
| `physProp.isCompressible: 0` (SMAC 非圧縮) | `methods/poisson.md`「コードベースに残っており利用できる」 | 経路ごと廃止するかのユーザ判断 |
| `time.deltaT.passiveFctSweeps` / `passiveFctTol` | 2026-09-16 新設で「使用 0」は年齢の反映 | 凝縮 dual-time の生産期に `[passiveFct] WARNING` が出ないことの確認 |
| `mesh.bndFirstOrder` | 使用禁止キー。削除は既存の [architecture-bndfirstorder-removal](architecture-bndfirstorder-removal.md) の責務 | (本 plan の対象外) |

#### 残すキー (P: 生産で選ぶ / V: 診断・A/B で要る)

**V の判定は「文書があるか」ではなく (codex plan M7)、いま必要な切り分け目的・実行手順・解除条件が言えるか**で行う。
言えないものは用途を確認してから分類し、「過去の手順が残っているだけ」を残置理由にしない。
`solver-settings.md` に位置づけ 1 行が無いキーは追記する (§5.1 #6)。

### 5.3 付随して見つかった不具合 (削除とは別に対処)

| # | 内容 | 状態 |
| --- | --- | --- |
| a | `recommended-settings.md` の SST レシピが**存在しないキー** `kInf` / `omegaInf` を推奨していた。実キーは `kInit` / `omegaInit` で、書いても黙って無視され k=ω=0 の初期値になる (node SST が step 0 で壊れる既知の原因) | **修正済 (2026-09-17)** |
| b | 同 §3 の `speciesPrecondDt: 1 (既定)` は化学ブランチ限定で main / sern に無い | **注記済 (2026-09-17)** |
| c | `thermCondMethod` (297 run が 1 を使用) と `prandtlLam` が `procedures/` `methods/` に一度も出てこない | §5.1 #6 で `solver-settings.md` に追記 |
| d | `initial` は solver 実行時には読むだけで使われない (`setInitial` の呼び出しは変換器のみ) のに全 run config で必須キー | 変換時キーへ移すか、必須をやめる。§5.1 #4 に含める |
| e | 存在しないキーを書いても黙って無視される (a の原因)。`check_solver_config.py` に「未知キーの検出」を足すべき | **済 (2026-09-18)**: solverConfig のソースにキー名が現れなければ WARN (偽陽性なし)。`kInf`/`omegaInf` を実際に検出 |
| f | `time.last.control` は**必須キーなのに消費者がゼロ**、`time.last.time` は**読む場所が無いのに 1988 run が書いている** (終了条件を時刻で指定する旧機能の残骸) | 第 2 陣で削除。全 run が書いているので段階移行 (受理 → 警告 → 削除) |
| g | `precondEps` (184 run) / `monitorInterval` (390 run) / `axisTimestepBeta` (39 run) が `procedures/` `methods/` に無い | §5.1 #6 で追記 |
| h | `speciesFaceReconstruction` を `solverConfig.cpp` が**二重に読んでいる** (値は同じで無害だが、片方を消すと齟齬) | 第 2 陣で整理 |
| i | 棚卸しツールが `config["time"]` 形を取りこぼし | **済 (2026-09-18)** |
| j | **棚卸しツールの使用数がキー名ベースで、完全修飾パス別になっていない** (節をまたいで合算・1 ファイル内の複数出現を重複計上)。これが「使用 0」の誤判定を生み、実際には 161 run / 1064 run が使うキーを削除しかけた | 第 2 陣の前に **PyYAML でパスごとに数える**実装へ置き換える (§5.1 #1 を再オープン) |
| k | 誤った場所に書かれて黙って無視されているキーが多数 | **検出は済 (2026-09-18)**: `check_solver_config.py` を完全修飾パス化し、本ツリーの 3455 config を掃いた結果が下表。修正は case README に注記して順次 |
| l | `physProp.isCompressible` と `physProp.ro`、`time.last.control` は**必須キーだが下流に消費者が無い** (パーサは読む; `solverConfig.cpp:312/615/646`)、`mesh.meshFormat` は合法値が `hdf5` 1 つだけ | 第 2 陣で**任意化 + 無効である旨の警告**から始める (codex plan-2 の助言)。`meshFormat` は省略時 `hdf5`・不正値は拒否のまま (`main.cpp:1094`) |
| m | `detectNaN` / `detectNaNInterval` が**トップレベルと `time.deltaT` の 2 か所**から読める (`solverConfig.cpp:399-403`, トップレベルが後勝ち)。実績はトップレベル 182 run / `time.deltaT` **0 run** | 第 2 陣で `time.deltaT` 側の綴りを落とす (使用ゼロなので移送不要) |
| n | `physProp.prandtlLam` は 797 run が書いているが**全て既定値 0.72**、しかも `procedures/` `methods/` に一度も出てこない | §5.1 #6 で `solver-settings.md` に追記 (削除でなく文書化。既定を変えたい用途は実在しうる) |

**§5.3 k の実測 (本ツリー `case/**/solverConfig.yaml` 3455 本, 2026-09-18)**

| 起動時に拒否される (FAIL) | run 数 | | 節の位置が違う (WARN) | run 数 | | どこにも無い (WARN) | run 数 |
| --- | ---: | --- | --- | ---: | --- | --- | ---: |
| `turbulence.LESorRANS` / `LESmodel` | 1259 | | トップレベル `lowMachPrecond` | 22 | | `time.last.time` | 1292 |
| `turbulence.RANSmodel` | 904 | | `physProp.lowMachPrecond` | 10 | | `turbulence.kInf` / `omegaInf` | 889 |
| `mesh.nodeAxisDirichlet` | 387 | | `condensation.condRealizProject` | 5 | | `space.keepDissipation` | 163 |
| `turbulence.DESmode` | 63 | | `space.keepDissType` ほか 3 件 | 各 4 | | `mesh.nodeWallViscGradFlux` | 104 |
| `mesh.nodeValueAtNode` ほか 3 件 | 9–48 | | `physProp.speciesImplicitCoupling` | 2 | | `time.implicit` | 48 |

拒否される側は**再実行すると落ちる**ので実害がすぐ出る。WARN 側は**黙って無視される**ので、
`lowMachPrecond` を書いたつもりの 32 run は前処理が効いていない。これが本 plan を起票した動機そのもの。

## 6. 検証 (2026-09-18 改訂, codex plan M8/M9)

初版の 3 run は全て node / SLAU / `turbulence.model: none` / `viscMethod 0` で、**削除する経路 (SST・粘性壁・line-implicit・cell・陽解法・KEEP) を一つも通らない**。
削除項目ごとに、その分岐が実際に作動する run を選ぶ。

### 6.1' 経路別の最小回帰表 (2026-09-18 改訂, codex plan-2 M7)

最終的な削除対象は 2 パス + 定数化 5 件 + 撤去 1 件だけなので、回帰も**その分岐を通る run に限る**。
戻した 3 件 (`blockDPLURDiagCache` / `blockDPLURDqPack` / `primPack`) は削除対象ではないので、
回帰ではなく**受理試験** (書いた config が起動して opt-in 経路に入ること) の対象にする。

| 削除項目 | 作動させる経路 | 回帰に使う run (複製して新 `run_NNNN_*` で) |
| --- | --- | --- |
| `implicitRelaxSST` (継承への統合) | SST + 陰解法 (block-DPLUR) | case/16 node SST ノズル |
| `lineDtWallRelief` | line-implicit (壁境界半割面) | case/39 の line-implicit 生産相当 (DDES) |
| `C_DES_kw` / `C_DES_ke` の定数化 | DES (`model: sst-ddes`) | 同上 |
| `wmlesNewtonTol` / `wmlesNewtonMaxIt` の定数化、`wmlesPrt` の撤去 | WMLES 壁モデル (`wallModelLES`) | WMLES の周期チャネル |
| `gradLSQDegenThresh` の定数化 | LSQ 勾配 (`mesh.gradLSQ: 2`) | node LSQ の run |

**等価性の主張 (回帰の前に済ませる審査)**: 上の 5 件はいずれも、キーを書かない既定 config に対して
**呼び出し側が渡す値が削除前と文字どおり同じ**である (`lineDtWallRelief` は既定 0 を 0 に、`implicitRelaxSST` は
`-1 → implicitRelax` の解決結果を `implicitRelax` に、定数 3 件は同じ数値に、`wmlesPrt` は誰も読まない引数の削除)。
したがって回帰は「許容差の判定」ではなく**この主張の反証試験**であり、差が出たらそれは丸めでなく欠陥である。
削除パスごとに「呼び出し側の旧式 / 新式 / 既定値」を 1 行で表に残す。

### 6.2' 判定基準 (M9, 2026-09-18 に事前確定 — codex plan-2 M7)

**比較量・ノルム・許容差・反復数を run を回す前に固定する**。後から基準を選ばない。

- **比較量 (全量を判定する)**: `ro`, `roUx`, `roUy`, `roUz`, `roe`, `P`, `T`。SST 系は `roK`, `roOmega` を追加。
  壁量 `Tau_Wall`, `Qw_Wall` (出力専用経路は保存量に出ないため必須)。**1 量でも外れたら不合格**。
- **ノルム**: 各量について相対 L2 差と相対 L∞ 差の両方。正規化は当該量の基準 run の L2 ノルム / 最大絶対値。
  絶対差は使わない (`P` の 2.125 のような値は桁が分からない)。
- **反復数とノイズ床**: **同一バイナリ・同一設定を 3 本**回し、**全ペアの最大**をノイズ床とする
  ([[air-condensation-cpg-carrier]] と同じ規約)。旧新差が**ノイズ床の 2 倍以内**なら合格。
- **決定的な局所試験はビット一致**: 単体試験と、1 step の決定的経路 (`atomicAdd` を通らない量) は旧新でビット一致を要求する。
- **定常量を報告する run は VERDICT 2 つ**: `check_convergence` PASS と `check_quasisteady` STEADY の両方。
  どちらかが欠ける run は「同時刻の非定常回帰」としてのみ使い、定常解の一致の根拠にしない。
- **キーの受理／拒否試験は完全修飾パスごと**に行う。削除した各パスが**起動時に落ちる**こと、および
  **残置した 9 件** (`blockDPLURDiagCache`, `blockDPLURDqPack`, `mesh.primPack`, `updateGuardAlpha`,
  `lowMachThornber`, `multispeciesRhoYCommonLimiter`, `passiveFctTolAbs`, `turbulence.turbulentSchmidt`,
  `physProp.isAxisymmetric`) が**受理されて意図した分岐に入る**ことを、それぞれ確認する。
  残す同名キー・残すモード値も対象 (`time.deltaT.control` を消さずに `time.last.control` を残す、など)。
- run は**新しい `run_NNNN_<slug>` に複製**して実行し、メッシュ品質・IC・出力先・case README の run 一覧・残差図まで
  通常の運用手順に乗せる。

**初回 (2026-09-18) の検証は不合格**: `run_0467` → `run_0482` / `run_0483` の比較は絶対差だけを見ており、
`P` (2.125 vs 2.0) と `roUy` (3.84e-3 vs 3.38e-3) が反復幅に収まっていなかった。3 run とも
`check_convergence: NOT CONVERGED` / `check_quasisteady --quantity pmax: TRANSIENT-UNSETTLED` で、
いずれも `model: none` / `viscMethod: 0` のため SST・DES・WMLES を通していない。
`run_0481_prune_reject_test` は**いま残置する** `lowMachThornber` を拒否した旧実装の記録なので、拒否試験としても使えない。
上の基準でやり直す。

### 6.4' 実施結果 (2026-09-18 やり直し)

基準バイナリは削除前の `f8bdec3d` を worktree `/home/sano/work/forge-prune-base` でビルドしたもの。
候補は HEAD (`build-native`)。判定は `solver_density_cuda/tools/check_field_regress.py` (本 §6.2' の実体化;
反復 3 本の全ペア最大をノイズ床、許容はその 2 倍、**数値的にゼロの量** (平面 2D の `roUz` 等) は `zero` として判定外)。

| 経路 (通す削除項目) | run (3 反復 + 基準) | 比較量 | 最大比 (ノイズ床比) | VERDICT |
| --- | --- | ---: | ---: | --- |
| SST 定常陰解法 (`implicitRelaxSST`, `gradLSQDegenThresh`) | `case/26.flat_plate_sst/run_0081`–`0084` | 15 | **1.17** | **PASS** |
| DDES (`C_DES_kw`, `C_DES_ke`) | `case/39.periodic_hills/run_0024`–`0027` | 10 | **1.50** | **PASS** |
| line-implicit + DDES 診断 (`lineDtWallRelief`, `C_DES_*`) | `case/39.periodic_hills/run_0028`–`0031` | 15 | **1.23** | **PASS** |
| WMLES 壁モデル (`wmlesNewtonTol`, `wmlesNewtonMaxIt`, `wmlesPrt` 撤去) | `case/38.channel_wmles/run_0024`–`0027` | 8 | **1.08** | **PASS** |

- **出力専用経路も判定に入れた**: SST は `wf_pk` (壁関数生産) と `vis_turb`、DES は `delta_les` / `l_des` / `rd_des` /
  `fd_shield`。`delta_les` と `wall_dist` は旧新でビット一致。
- **拒否試験 8/8 PASS**: `case/26.flat_plate_sst/run_0085_prune_keytest` (`logs/reject_*.log`)。削除した各パスを書いた
  config が `no longer supported` を出して非ゼロ終了する。
- **受理試験 9/9 PASS**: 同 run (`logs/accept_*.log`)。残置した 9 件をそれぞれ書いた config が起動して完走する。
- **opt-in 受理 3/3 PASS**: `case/26.flat_plate_sst/run_0086_optin_*`。復元した性能スイッチ 3 件が既定経路と同じ場を出す
  (全量 0.96〜1.02 倍)。**この run は反復間でもビット一致しない** ので、plan が言う「ビット同一」は
  `atomicAdd` 経由の残差集積がある run では検証できず、ノイズ比で判定した。
- **VERDICT の但し書き**: SST 回帰は収束済みの場からの継続なので `check_convergence` は旧新とも
  `NOT CONVERGED (plateau)` (残差が初手から床)。`check_quasisteady` は旧新とも `ALL STEADY`。
  したがって**同一 IC・同一 step の場の一致**として読むこと。定常解の収束を主張する試験ではない。
- **意図的な逸脱**: DDES の回帰 config は `implicitRelax: 0.5` のまま回した (`check_solver_config` は
  dual-time なので FAIL を出す)。1.0 にすると旧実装の `implicitRelaxSST: -1 → implicitRelax` の継承が
  恒等になり、継承統合の誤りを検出できないため。生産設定としては §6 の dual-time レシピに従うこと。

### 6.3' 生成器・変換器 (M4)

- `design/forge_design/evaluate/runner*.py` の 3 本が生成する `solverConfig.yaml` から削除キーを外し、**生成器を直してから**起動時拒否を入れる。
- `mesh/convertGmshToForge.cpp` は同じ `solverConfig::read()` を使い `cfg.gpu` を確保に渡すので、**solver の GPU 固定と変換器のホスト実行を分けて**扱う。
- `initial` は solver 実行時には使われない (`setInitial` の呼び出しは変換器のみ) ので、**呼出側ごとに必須性を定義**する。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-17` | [2026-09-17-config-key-pruning-plan.md](../../notes/reviews/2026-09-17-config-key-pruning-plan.md) | **GO-with-changes**, C0/M9/m0 | **全採用 (2026-09-18, §5.2 を改訂)**: M1 `timeIntegration: 1` は 1 段 Euler で 3 (3 段 TVD RK) の別名ではない (`solverConfig.cpp:907/919` で段数・係数が別) → 削除対象から除外; M2 `sstIsotropicStress`/`sstEnergyKSource` の引用先は最終決定と逆 (accepted plan §4.0 は既定 0 + 個別利用可、残作業表も「撤去せず現状維持」) → 保留; M3 「全 155 キー分類済み」が不成立 (棚卸しが `config["mesh"][...]` 形を取りこぼし、同名末端キーを統合していた; 分類表の内訳も 31 と合わない) → **完全修飾パスで再抽出 (188 パス)** し分類をやり直す; M4 共有パーサ・生成器が影響範囲から漏れ (`convertGmshToForge` が同じ `solverConfig::read()` を使い `cfg.gpu` を確保に渡す、design の runner 3 本が `gpu`/`nodeWallDirichlet` を生成) → スコープに追加; M5 `condEquilibrium 1→2` の「固定点同一」は一般には誤り (新しい accepted plan が緩和形の固定点は Δτ 依存と明記、mode 2 は単一凝縮種限定) → 保留; M6 `mesh.axisymMethod: 1` は不採用決定ではなく opt-in 保持で再評価待ち → 保留; M7 A/D/V が限定的な実測から機能廃止へ飛躍 (`lineVisc*` は「このケースでは僅差」、`condTwoTemp` は希薄水/N2 限定、化学 2 キーは反応源の温度評価・停止条件を変える) → 保留に移し、A は「検証した適用範囲」と「廃止する対応範囲」を対応づける; M8 §6 の回帰 3 run が全て node/SLAU/`model: none`/`viscMethod 0` で、SST・粘性壁・line-implicit・cell・陽解法・KEEP の経路を通せない → 経路別の最小回帰表に作り替え; M9 「長時間 run のビット一致」は非退行判定として実行不能 (残差が float atomicAdd で集積、出力専用経路は保存量に出ない、`implicitRelaxSST=-1` は継承指定で定数置換と別物) → 決定的な局所試験はビット一致、CFD は同一バイナリ反復幅 + 事前に定めた物理量許容。保留方針 (`ducrosLimiter`/`condLimiterMode`/`nodeOmegaWfDirichlet`) は同意を得た |
| plan | `2026-09-18` | [2026-09-18-config-key-pruning-plan.md](../../notes/reviews/2026-09-18-config-key-pruning-plan.md) | **GO-with-changes**, C0/M7/m1 | **全採用 (2026-09-18)**: M1 性能スイッチ 3 件 (`blockDPLURDiagCache`/`blockDPLURDqPack`/`primPack`) は元 plan が「不採用で確定、**opt-in 残置**」と決めており、残置決定の読み違いが 3 度目 → **コードごと復元** (§5.2)。M2 棚卸しが入れ子の run config を探索せず「使用 0」を再び誤判定 (826 config が対象外、`forge-perf` が実際に 3 件とも 1 を書いている) → `case/**` を再帰探索し、`<case>/<相対パス>` を識別子に、内容違いの同名は別設定、読み取り失敗は報告 (2796 → **4030 本**)。M3 「live 181 パス」が受理キー一覧になっていない (コメントを走査、節・拒否専用参照を live に混入、必須キーの既定値に節名が入る) → コメント除去 + 節/拒否専用/必須の分離で **値キー 167 / 節 10 / 拒否専用 12 / 起動時拒否 8**。M4 既定値比較の絶対許容差 1e-12 が 10 桁の変更を「既定と同じ」にする + 別名 (`mesh.renumber`→`meshRenumber` 等) を解決できない → 正規化した厳密比較へ、既定値は**代入先メンバ名**で引き、解決できないものは「不明」とし、記載 run 数と非既定 run 数を分離。M5 未知キー検出が末端名照合で誤配置を拾えない → **完全修飾パス**照合に置換 (誤配置は WARN、起動時拒否は FAIL)。実例 (トップレベル `lowMachPrecond` 22 run、`physProp.lowMachPrecond` 10、`condensation.condRealizProject` 5、`space.keepDiss*` 各 4、`physProp.speciesImplicitCoupling` 2) を検出、回帰試験にも追加。M6 正本文書が復元済みキーを「起動時エラー」と案内 → `recommended-settings.md` §9.1 を最終対象に同期し残置 9 件を明記、plan の確定表・残作業表・検証表も同期。M7 非退行の主張が全量では成立せず受入試験も不足 → §6.1'/§6.2' を**事前確定の基準**に書き換え検証を再オープン (§5.1 #5)。m8 `wmlesPrt` は定数化でなく**消費者のない引数の撤去**で、`turbulentPrandtl` への移送案内も誤り → メンバ・カーネル引数ごと削除し案内を撤回。C の段階移行方針は妥当との評価。`isCompressible`/`ro`/`time.last.control` は「読まれない」でなく「パーサは読むが下流の消費者が無い」、`mesh.meshFormat` は省略時 `hdf5` + 不正値拒否を維持、という助言も採用 (§5.3 l) |

## 7. 影響範囲

solver の config 読込と分岐、`procedures/` の設定文書、skill `forge-config`、既存 run config (削除キーを書いているものは
再実行時に起動エラーになる → case README に注記)。

## 8. 完了条件

- [ ] 分類表 (§5.1 #2) が完全修飾パス 188 件を覆い、codex plan レビュー 2 回目を通っている
- [ ] A–D の削除と旧キーの起動時エラー化
- [ ] §6 の 1–4 を満たす
- [ ] `solver-settings.md` / `recommended-settings.md` / skill の同期

## 9. 変更ログ

- `2026-09-18` — **検証をやり直して全 PASS** (§6.4', codex plan-2 M7)。削除前バイナリ `f8bdec3d` を worktree `forge-prune-base` でビルドし、4 経路 (SST 定常陰解法 / DDES / line-implicit + DES 診断 / WMLES 壁モデル) で「同一バイナリ 3 反復のノイズ床 × 2」を事前に基準として判定。全量が 0.70〜1.50 倍で **4 経路とも PASS**、出力専用の `wf_pk` / `l_des` / `rd_des` / `fd_shield` も含む。削除キーの拒否 8/8、残置キーの受理 9/9、復元した opt-in スイッチ 3/3 も PASS。判定は恒久ツール `solver_density_cuda/tools/check_field_regress.py` に実体化した。
- `2026-09-18` — codex plan レビュー **2 回目 GO-with-changes (C0/M7/m1)** を全採用 (§6.1)。**M1 で 3 度目の残置決定の読み違い**が出た: `blockDPLURDiagCache` / `blockDPLURDqPack` / `mesh.primPack` は元 plan が「不採用で確定、**opt-in 残置**」と決めたキーで、「既定に採用しない」を「機能を廃止する」と取り違えていた → コード (config 読み・メンバ・`calcGradient_d.cu` / `limiter_d.cu` / `timeIntegration_d.cu` の分岐) ごと復元し、削除前のファイルと diff ゼロを確認。**棚卸しも 2 度目の「使用 0」誤判定**で、`case/*/run_*/solverConfig.yaml` に限定していた探索を `case/**` の再帰に直すと 2796 → 4030 本になり、`forge-perf` の nested run が 3 件とも実際に 1 を書いていた。ツールはさらに、コメント除去・節/拒否専用参照/必須キーの分離 (181 → **値キー 167**)・別名を代入先メンバ名で解決・絶対許容差 1e-12 の撤去 (1e-20 と 1e-30 を同一視していた)・記載/非既定の分離・読み取り失敗の報告を入れた。`check_solver_config.py` の未知キー検出は**完全修飾パス**照合に置き換え、誤配置 (`lowMachPrecond` を トップレベルに 22 run 等) を WARN、起動時拒否を FAIL で拾うようにして実例を回帰試験に追加。`wmlesPrt` は定数化でなく**消費者のない引数の撤去**に分類し直し (Kader 原式で Pr_t が不要)、メンバ・カーネル引数ごと削除。`recommended-settings.md` §9.1 を最終対象に同期し、**残置した 9 件**を明記。**最終的な削除は 2 パス (`lineDtWallRelief`, `implicitRelaxSST`) + 定数化 5 件 + 撤去 1 件**。検証は §6.2' の事前確定基準でやり直す (再オープン)。

- `2026-09-18` — **自己訂正 2**: さらに 2 件を戻した。`turbulence.turbulentSchmidt` は「全 worktree で使用 0」が誤りで**実際は 161 run が使用中** (chem ブランチ case/48 が 0.7 ×124、case/47 が **0.5** ×37)。`physProp.Sc_t` の既定と違う値があるので別名撤去には設定値の移送が要る。`physProp.isAxisymmetric` は **1064 run (値 1 が 593)** が使っており、拒否にすると既存 run の再実行が全部落ちる (deprecated 読みは生きている)。**原因**: 棚卸しツールの使用数集計が完全修飾パスでなくキー名ベースで、しかも節をまたいで合算・重複計上していたため「使用 0」を誤って出していた。**第 1 陣で実際に削除できたのは 5 パス + 定数化 6 件**に縮小。
- `2026-09-18` — **自己訂正**: 第 1 陣で削除した 17 件のうち **4 件を戻した** (`updateGuardAlpha`, `lowMachThornber`, `multispeciesRhoYCommonLimiter`, `passiveFctTolAbs`)。前 3 者は**元 plan の「処置」が明示的に「opt-in のまま残置」**で、codex plan-1 M2 が指摘したのと同じ「古い実装ステップや限定的な実測を最終決定と取り違える」誤りを繰り返していた。`passiveFctTolAbs` は同一 plan・同一日・同一機構の `passiveFctSweeps`/`passiveFctTol` を保留にしながら 1 つだけ削っており §5.2 内で不整合だった。4 件とも保留表へ移し、解除条件を書いた。併せて棚卸しツールが `getOptionalValidatedValue<int>(config["time"], ...)` 形を取りこぼしていた (`nStepInner` / `nSubIterDualTime` / `bdfOrder` が抜けていた) のを修正。
- `2026-09-18` — **新発見 (分類やり直し)**: `time.last.control` は**必須キーなのに solver 全域で消費者がゼロ**、`time.last.time` は**読む場所が存在しないのに 1988 run が書いている**。§5.3 f として記録し、削除は「全 run が書いているので段階移行」とする。
- `2026-09-18` — **第 1 陣の削除を実装**: 確定 5 パス (`blockDPLURDiagCache`, `blockDPLURDqPack`, `primPack`, `lineDtWallRelief`, `implicitRelaxSST`) と定数化 6 件 (`C_DES_kw`, `C_DES_ke`, `wmlesNewtonTol`, `wmlesNewtonMaxIt`, `wmlesPrt`, `gradLSQDegenThresh`)。旧キーは起動時エラー (移行先つき) にした。**検証**: (a) `lowMachThornber: 1` を書いた config が起動時に落ちることを実 run で確認 (case/44 `run_0481`)、(b) 既定挙動の非退行 — 削除前 `run_0467` と削除後 `run_0482` の差が同一設定の反復ノイズ `run_0483` 以下 (ρ 1.19e-5 vs 1.29e-5、roe 6.1 vs 9.7、roQ0 6.6e10 vs 1.28e11)、(c) 単体試験 6 本 ALL PASS。`check_solver_config.py` に**未知キー検出**を追加 (solverConfig のソースにキー名が現れなければ WARN; `kInf`/`omegaInf` を実際に検出)。
- `2026-09-18` — codex plan レビュー 1 回目 **GO-with-changes (M9)** を全採用 (§6.1): 誤分類 4 件 (`timeIntegration: 1` は 1 段 Euler で 3 の別名ではない、SST 分離型 2 キーは引用先の最終決定と逆、`condEquilibrium 1→2` の固定点は Δτ 依存、`mesh.axisymMethod: 1` は再評価待ち) を保留へ移し、限定的な実測から飛躍していた 4 群 (`lineVisc*`, `condTwoTemp`, 化学 2 キー, 感度係数 2 件) も保留に。棚卸しを**完全修飾パス**で再抽出 (155 名称 → 188 パス) し分類を再オープン。スコープに変換器と設計チェーンの生成器を追加。§6 を経路別の最小回帰表と反復幅基準に作り替え。**今回確定した削除は 10 パス + 定数化 7 件**。
- `2026-09-17` — 起票 (ユーザ指示)。棚卸し (155 キー / 未使用 22 / 既定のみ 3 / recommended 未記載 78; 全 9 worktree の run config 5163 本を対象)。
- `2026-09-17` — 全キーの分類を根拠つきで作成 (§5.2)。付随して `recommended-settings.md` の存在しないキー 2 件を修正・注記 (§5.3 a/b)。

## 10. 未確定事項

- 「診断・A/B で要る」の線引き: 過去 1 回しか使っていないスイッチを残すか。**案**: 切り分けの再現手順が plan/notes に
  書かれているものだけ残し、それ以外は削除して必要になったら復活させる (git にある)。
- 既定値を変えるか: `lineKFreeze` のように「dual-time では常用推奨」なのに既定 0 のものは、削除ではなく既定変更が適切。
