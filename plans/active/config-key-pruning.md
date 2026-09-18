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

### 4.1 必須キーの段階移行 (第 2 陣, 2026-09-18)

第 1 陣は「使用ゼロの opt-in スイッチ」だったので、いきなり起動時エラーにしてよかった。**第 2 陣は全 run が書いている
必須キー**なので同じ手は使えない — 一足飛びに拒否すると 4000 本超の既存 run が再実行できなくなる。段階を分ける。

| 段階 | 挙動 | いつ次へ進むか |
| --- | --- | --- |
| **S1 任意化 + 警告** (今回やる) | キーが無くても起動する。**書いてあると 1 行警告**を出す (「このキーは読まれていない / 合法値が 1 つしかない」と理由つき)。値は無視するか、これまでと同じ既定を使う | 生産チェーン (`design/forge_design/evaluate/runner*.py`) と `procedures/` の雛型から外し、新規 run に混入しなくなってから |
| **S2 起動時エラー** (次回以降) | 第 1 陣と同じ `removed[]` 行に移す | S1 の警告が実運用のログから消えてから |

**S1 で守ること** (codex plan-4 M2/M3/M4 反映):

- **既定の挙動を変えない**。
  - `mesh.meshFormat`: 省略時 `hdf5`。**不正値の拒否は維持する** (`main.cpp:1094` の分岐を残す)。
    「合法値は 1 つ」は**現在の solver 側の制約**であって変換器の制約ではない (変換器は Gmsh を直接読み、この分岐を通らない)。
  - 消費者のない 3 件 (`physProp.isCompressible` / `physProp.ro` / `time.last.control`): 読むのをやめて**メンバごと撤去**する。
    どこからも読まれていないことは grep で確認済み (§5.2 第 2 陣の表)。**架空の既定値を置く必要はない**。
  - `time.deltaT.detectNaNInterval`: **読みを維持する**。この別名は無効ではなく、下位だけに書けば検査間隔が実際に変わる
    (`solverConfig.cpp:402`, `main.cpp:1900`)。解決順序は **トップレベル → 旧別名 → 既定 1**、最後に**下限 1 への補正**
    (従来どおり)。S1 は「トップレベルへ移してください」の警告だけ。拒否は S2。`detectNaN` の両綴りは触らない。
- **警告は「キーごと・1 プロセスの設定読込みにつき 1 回」**。毎 step 出さない。変換器と solver が同じ run でそれぞれ
  読むときは、それぞれ 1 回ずつ出る (定義をこう置くと曖昧さが無い)。
- `check_solver_config.py` が同じキーを **WARN** で拾えるようにする。**既知の値キーは無条件に通る実装**なので、
  パーサ側に警告を足すだけでは足りない — S1 対象の明示リストを検査器に持たせる。
- `recommended-settings.md` §9.1 に「S1 段階のキー」の表を作り、削除済み (起動時エラー) と**区別して**書く。

**S1 の作業順序** (codex plan-4 M3。第 1 陣とは逆向き):

1. パーサを任意化する (solver 側)。
2. **実際に使う `forge` と `convertGmshToForge` を両方リビルドして確認**する。生産チェーンの runner は
   `solver_density_cuda/build/` を見る (`design/forge_design/evaluate/runner.py:31`) ので、`build-native` だけ更新しても
   チェーンの確認にはならない。
3. そのあとで生成器 (`runner.py` / `runner_wt.py` / `runner_sern3d.py`) と `procedures/` の雛型から旧キーを外す。

**先に生成器を直すと起動できなくなる** (旧パーサはまだ必須にしている)。§6.3' の「生成器を直してから起動時拒否」は
**S2 向けの順序**であって S1 には当てはまらない。今回対象外の `gpu`・`initial` の必須性は**混ぜない**。

**S2 へ進む条件** (m5): 「新規生成から警告が消えた」だけでは足りない。**保持している再実行対象 run の移送確認**か、
**互換性を打ち切る対象・時期・移行手順の明示**のどちらかを済ませること。S1 完了から自動で S2 へ進めない。**別マイルストーンとする**。

### 4.2 S2 (起動時エラー化) の設計 (2026-09-18, ユーザ指示で着手)

**条件の満たし方は「移送」と「打ち切りの明示」の両方**を行う (片方でよいが、既存 run が 7000 本あるので両方やる)。

1. **移送ツールを作る**: `solver_density_cuda/tools/migrate_solver_config.py`。`case/**/solverConfig.yaml` から
   S1 の 5 キーを外す。**節を特定してから編集する**ので `time.last.control` を消すつもりで `time.deltaT.control`
   を消す事故が起きない。書き戻す前に**編集後の YAML が「元 − 消したキー」と厳密一致**することを確かめ、
   1 か所でも違えば書かない。コメントと整形は保つ (行単位の最小編集; `--reformat` で PyYAML 整形も可)。
   - `time.deltaT.detectNaNInterval` は**値を捨てない**。トップレベルに無ければ**移送**する
     (この別名は無効キーではなく、書いた値が実際に効く)。トップレベルに既にあれば、そちらが優先なので捨ててよい。
2. **このツリーの `case/` を移送する**。他の worktree は各ブランチの作業物なので触らない — ツールを渡して、
   そのブランチで必要になったときに回してもらう (並行セッションの作業を壊さないため)。
3. **打ち切りを明示する**: `recommended-settings.md` §9.2 を「S1 (警告)」から「削除済み (起動時エラー)」へ移し、
   **対象 5 キー・日付・移行コマンド 1 行**を書く。
4. **実装**: 5 キーを `removed[]` の起動時拒否テーブルへ移す。`check_solver_config.py` の `STAGED_KEYS` は
   **FAIL 側へ移す** (WARN では投入前に止められない)。`main.cpp` の `meshFormat` 分岐は**残す**
   (省略時 `hdf5`、不正値は拒否、という契約は S1 のまま変わらない)。

**S2 の合否は §6.7'**。

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
| 2 | 全キーの分類表 (live な値キー 167 パス) | **済 (2026-09-18)**: 修正後の棚卸し (値キー 167 / 節 10 / 拒否専用 12 / 起動時拒否 8、run config 4063 本) で分類を完了。素データは [notes/investigations/config-key-inventory-2026-09-18.md](../../notes/investigations/config-key-inventory-2026-09-18.md)、結論は §5.2「分類の結論」。**消してよいと言えるのは第 2 陣の候補 6 件だけ** |
| 9 | 第 2 陣 S1 (任意化 + 警告) | **済 (2026-09-18)**。対象 5 件を任意化し 1 回だけ警告。消費者のない 3 件はメンバごと撤去、`meshFormat` は省略時 `hdf5`、`detectNaNInterval` の別名は読みを維持。受入は §6.6' で全 PASS |
| 10 | 第 2 陣 S2 (起動時エラー) | **別マイルストーン**。S1 完了では進めない — 保持する再実行対象の移送確認、または互換打ち切りの対象・時期・手順の明示が条件 (m5) |
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

#### 分類の結論 (2026-09-18, 修正後の棚卸しで完了 — §5.1 #2)

素データは [notes/investigations/config-key-inventory-2026-09-18.md](../../notes/investigations/config-key-inventory-2026-09-18.md)
(全 worktree の `case/**/solverConfig.yaml` **4063 本**)。内訳:

| 区分 | 件数 | 処置 |
| --- | ---: | --- |
| live な値キー | **167** | 下の 5 行に分解される |
| ├ 非既定の使用がある | 119 | **残す** (P/V)。生産で選ばれているか、A/B の実績がある |
| ├ 必須キー (既定なし) | 30 | うち 4 件が第 2 陣の候補 (下表)、他は残す |
| ├ 既定値が解決できない (リスト・マップ値) | 5 | `bodyForce` `mesh.wallDistExtraPhysIDs` `output.extraFields` `physProp.species` `space.uRef`。いずれも実使用があり**残す** (「既定と同じか」を機械判定できないだけ) |
| ├ どの run にも書かれていない | 10 | 9 件は保留表に載っているもの。残り 1 件が第 2 陣の候補 |
| └ 書かれているが既定値のみ | 3 | `keepDissCbEps` (保留)、`physProp.prandtlLam` (**残す**: 文書化した。壁関数の回復係数でも読む)、`time.deltaT.speciesImplicitRelax` (第 2 陣で要判断) |
| 拒否専用の参照 | 12 | 起動時に落ちる旧キー。**現状維持** (`nodeAxisDirichlet` 系 5・旧乱流 4 キー・`nInnerLoop` 系 3) |
| 起動時拒否テーブル `removed[]` | 8 | 第 1 陣で削除したもの |
| 節そのもの | 10 | キーではない |

**分類は完了**。「消してよい」と言えるのは、上の 167 のうち**第 2 陣の候補 6 件だけ**である。

#### 第 2 陣の候補 (5 件, 2026-09-18 に codex plan-4 で確定)

| パス | 実績 | 根拠 | 処置案 |
| --- | ---: | --- | --- |
| `physProp.isCompressible` | 必須・4052 run | メンバ `isCompressible` を**パーサ外で読む場所がゼロ** (`solverConfig.cpp:615` で格納するだけ。対照の `dtControl` は `main.cpp:510` ほかで読む) | **任意化 + 無効である旨の警告**。SMAC 非圧縮経路を廃止するかの判断は別 |
| `physProp.ro` | 必須・4052 run | 同上 (`solverConfig.cpp:646`)。圧縮性では密度は EOS で決まる | 同上 |
| `time.last.control` | 必須・4052 run | メンバ `endTimeControl` を読む場所がゼロ (`solverConfig.cpp:313`)。終了条件を時刻で指定する旧機能の残骸 | 同上 |
| `mesh.meshFormat` | 必須・4052 run | 合法値が `hdf5` 1 つだけ (`main.cpp:1094`) | **省略時 `hdf5`** とし、不正値の拒否は維持する |
| `time.deltaT.detectNaNInterval` | **0 run** | トップレベル `detectNaNInterval` (182 run) と同じものを読む 2 つ目の綴り (`solverConfig.cpp:402-403`)。使用ゼロなので移送は要らない。**`detectNaN` の方は 2781 run が使っているので触らない** (§5.3 m) | **S1 では読みを維持**し「トップレベルへ移せ」と警告するだけ (codex plan-4 M2: この別名は無効ではなく、下位だけに書けば実際に間隔が変わる)。拒否は S2 |

**対象から外したもの (codex plan-4 M1)**: `time.deltaT.speciesImplicitRelax` は「10 run すべて既定 1.0」だが、
**無効キーではない** — `speciesTransport_d.cu:835` が読み、`scalarTransport_d.cu:290` で化学種更新の増分に掛かる
(`relax * (res·dt/V)/fac`)。非既定値を無視すれば `speciesImplicitCoupling: 0` の反復写像が変わる。
accepted plan [species-passive-scalar-unification](../accepted/species-passive-scalar-unification.md) が独立キーとして設計し、
効果確認が残件として残っている。**現状維持**とし、今回の段階移行から外す (掃引を新たにやる必要も無い)。
「使用実績が既定値だけ」は**不要性の実測ではない** — これは §4 の D 基準 (一度も使われず、既定から動かす根拠も無い
内部チューニング定数) を、消費者の有無を見ずに適用しかけた誤り。

**いずれも段階移行**にする (受理 → 無効である旨の警告 → 削除)。全 run が書いている必須キーを一足飛びに拒否すると、
既存 run の再実行が全部落ちる。`time.last.time` は**そもそも live キーでない** (1292 run が書いているが solver は読まない)
ので削除対象が無い。`check_solver_config.py` が WARN で拾うので、順次 config から外す。

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
| m | `detectNaN` / `detectNaNInterval` が**トップレベルと `time.deltaT` の 2 か所**から読める (`solverConfig.cpp:399-403`, トップレベルが後勝ち) | **訂正 (2026-09-18, codex plan-3 M1)**: 当初「`time.deltaT` 側は 0 run」と書いたが、**2 パスを取り違えていた**。`time.deltaT.detectNaN` は**記載 2879 / 非既定 2781 run** (例: `case/44.vitiated_air_wt/run_0117_va_ns_ar5k_iso300_samewall/solverConfig.yaml:12`)、0 run なのは `time.deltaT.detectNaNInterval` の方。**削除決定は撤回**し、`detectNaN` は互換読み (トップレベル優先) を維持する。統合するなら入力の移送が要る |
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

| 経路 (通す削除項目) | run | 比較量 | 最大比 (ノイズ床比) | VERDICT |
| --- | --- | ---: | ---: | --- |
| SST 定常陰解法 (`implicitRelaxSST`, `gradLSQDegenThresh`) | `case/26.flat_plate_sst/run_0081`–`0084` (新 3 + 旧 1) | 24 | **1.27** | **PASS** |
| DDES (`C_DES_kw`, `C_DES_ke`) | `case/39.periodic_hills/run_0024`–`0027`, `0032`–`0035` (新 5 + 旧 3) | 28 | **1.27** | **PASS** |
| line-implicit + DDES 診断 (`lineDtWallRelief`, `C_DES_*`) | `case/39.periodic_hills/run_0028`–`0031` (新 3 + 旧 1) | 33 | **1.25** | **PASS** |
| WMLES 壁モデル (`wmlesNewtonTol`, `wmlesNewtonMaxIt`, `wmlesPrt` 撤去) | `case/38.channel_wmles/run_0024`–`0027` (新 3 + 旧 1) | 26 | **1.42** | **PASS** |

**壁量は境界出力から採る (codex plan-3 M3)**: `Tau_Wall` / `Qw_Wall` は `output_cellValNames` に無いので `res_<step>.h5`
には出ない。実体は `outputHDFflg: 1` の bcond が書く `res_<名前>_<physID>_<step>.h5` の `twall_x/y/z`・`qwall`・`utau`・`ypls`
で、`check_field_regress.py --boundary` がこれを比較する。上表の比較量にはこれらが入っている
(4 経路とも壁せん断応力・摩擦速度・y⁺ が判定対象、WMLES は `qwall` も)。

**ノイズ床は両側から測る (2026-09-18 の修正)**: 当初は新バイナリ 3 反復だけで床を作っていたが、**カオス的な run では
桁で足りない**。周期丘 DDES 400 step の壁せん断 `yhi_4/twall_x` は、新 3 本の床 1.81e-3 に対し 5 本にすると 2.70e-2 と
**15 倍**になり、3 本での判定は比 2.00 の**偽の不合格**を出した。新 5 本 + 旧 3 本で測り直すと、旧内 3.36e-2・新内 2.70e-2 に対し
交差 4.23e-2 = 1.26 倍で、全量が 1.27 倍以内に収まる。ツールは `--candidate` を複数取れるようにし、
**床 = max(基準側の全ペア, 候補側の全ペア)、比較 = 基準×候補の全ペア**とした。**LES/DES は両側 3 本以上**を要求する。

- **出力専用経路も判定に入れた**: SST は `wf_pk` (壁関数生産) と `vis_turb`、DES は `delta_les` / `l_des` / `rd_des` /
  `fd_shield`、全経路で境界出力の壁量。`delta_les` と `wall_dist` は旧新でビット一致。
- **判定器の不備検査 (codex plan-3 M2)**: 比較の中では NaN を検出できない (`max(0.0, NaN)` が `0.0` になる) ため、
  必須量の欠落・形状不一致・非有限値を**数値判定の前に**検査して終了コード 2 で落とす。反例は単体試験
  `solver_density_cuda/tests/unit/test_check_field_regress.py` に入れた。
- **拒否試験 8/8 PASS**: `case/26.flat_plate_sst/run_0085_prune_keytest` (`logs/reject_*.log`)。削除した各パスを書いた
  config が `no longer supported` を出して非ゼロ終了する。
- **受理試験 9/9 PASS (パーサ受理)**: 同 run (`logs/accept_*.log`)。残置した 9 件をそれぞれ書いた config が起動して完走する。
  **ただしこの run は単成分・`viscMethod 0`・定常なので、3 件は分岐に届いていない** (codex plan-3 M4)。下の到達確認で補った。
- **分岐到達の確認 3/3 PASS** (受理と分けて記録する):

  | キー | 到達の観測 | run |
  | --- | --- | --- |
  | `multispeciesRhoYCommonLimiter` | 診断印字 `RHOYLIM call=` がログに出る (2 件) | `case/44.vitiated_air_wt/run_0487_reach_rhoylim` (2 成分・`speciesFaceReconstruction 2`) |
  | `passiveFctTolAbs` | FCT 低次解の相対残差が **4.32e-8 → 9.73e-22** (14 桁) 変わる。この値は分母の加算項なので、読まれなければ動かない | `run_0488_reach_fcttolabs_zero` (0) と `run_0489_reach_fcttolabs_big` (1e30)、いずれも dual-time + SLAU + SFR 2 |
  | `turbulence.turbulentSchmidt` | `0.2` にするとトレーサ `roXi` が既定 (0.7) から**相対 L2 で 1.41e-1** 動く。`physProp.Sc_t: 0.2` との差は **5.6e-6** = 別名として等価 | `case/16.nozzle_wys/run_0502`–`0504` (SST + `viscMethod 1` + トレーサ) |
- **opt-in 受理 3/3 PASS**: `case/26.flat_plate_sst/run_0086_optin_*`。復元した性能スイッチ 3 件が既定経路と同じ場を出す
  (全量 0.96〜1.02 倍)。**この run は反復間でもビット一致しない** ので、plan が言う「ビット同一」は
  `atomicAdd` 経由の残差集積がある run では検証できず、ノイズ比で判定した。
- **VERDICT の但し書き**: SST 回帰は収束済みの場からの継続なので `check_convergence` は旧新とも
  `NOT CONVERGED (plateau)` (残差が初手から床)。`check_quasisteady` は旧新とも `ALL STEADY`。
  したがって**同一 IC・同一 step の場の一致**として読むこと。定常解の収束を主張する試験ではない。
- **cell の回帰は組まない**: SST 緩和 (`update_d.cu:356`) と WMLES カーネルは cell/node 共有だが、
  **2026-09-16 のユーザ決定「cell はもう使わない」**により検証・回帰は node のみとする (codex plan-3 M5 は棄却)。
  `procedures/verification/README.md` の「node/cell 両方」規則をこの決定に合わせて書き換えた。
  **残る risk**: cell 側の共有変更は未検証のまま。cell を再び使うなら、その時点でやり直す (§10 に記載)。
- **意図的な逸脱**: DDES の回帰 config は `implicitRelax: 0.5` のまま回した (`check_solver_config` は
  dual-time なので FAIL を出す)。1.0 にすると旧実装の `implicitRelaxSST: -1 → implicitRelax` の継承が
  恒等になり、継承統合の誤りを検出できないため。生産設定としては §6 の dual-time レシピに従うこと。

### 6.5' 第 2 陣 (S1) の受入表 (2026-09-18, codex plan-4 M4)

第 1 陣の合否基準 (§6.2') は「削除キーが**起動時に落ちる**こと」を要求するが、S1 で要るのは逆で
「**旧入力を受理し、省略した入力も同じように動く**」こと。第 2 陣専用に置く。

| 対象 | 必須の判定 |
| --- | --- |
| 消費者のない 3 件 (`physProp.isCompressible` / `physProp.ro` / `time.last.control`) | 個別省略・一括省略・旧値記載のいずれも**受理**。残る設定値が旧実装と同じ。**旧キーを書いたときだけ**警告が出る |
| `mesh.meshFormat` | 省略と `hdf5` で実効値が同じ。**不正文字列は既定値へ置換せず拒否**。共有パーサで検証し、変換器でも同じ契約 |
| `detectNaNInterval` の別名 | 両方省略 → 1 / 旧別名のみ → 旧値 / トップのみ → トップ値 / 併記 → **トップ優先** / 0 以下 → 1 |
| 残置キー (巻き込み事故の検出) | `time.deltaT.control`、**両方の** `detectNaN`、`speciesImplicitRelax: 0.7` の実効値が変わっていないこと |
| 検査ツール | S1 対象**だけ** WARN。省略時はその WARN が出ない。S1 キーを削除済み (FAIL) に混入させない |
| 生成器・実行経路 | runner 3 本のテンプレートとそれを再利用する経路から旧キーが消えていること。**生成 → 変換 → node solver 起動**を通す |

**数値回帰** (M4): 実効設定の厳密比較に加え、**標準 node ケース 1 件の短い固定 step 比較**で足りる。
同一 IC で「変更前 + 旧入力」「変更後 + 旧入力」「変更後 + 省略入力」を比べ、§6.2' の非有限検査と反復ノイズ基準を当てる。
**第 1 陣の SST / DDES / WMLES 一式を再実行する理由は無い** (今回の変更は分岐を通らない)。
短時間回帰を収束の証明にはしない。定常量を主張するときだけ §6.2' の収束・準定常 VERDICT を別途要求する。

### 6.6' 第 2 陣 S1 の実施結果 (2026-09-18)

基準バイナリは S1 前の `d1501dc7` (worktree `forge-s1-base`)、候補は S1 後。**生産チェーンが使う
`solver_density_cuda/build/` の `forge` と `convertGmshToForge` を両方リビルドして確認**した (plan-4 M3)。

| 受入項目 | 結果 | 根拠 |
| --- | --- | --- |
| 旧入力 (5 件記載) を受理、記載したキーだけ警告、警告は 1 回ずつ | **PASS** (11 項目) | `case/26.flat_plate_sst/run_0087_s1_keytest/logs/old_input.log` |
| 個別省略 5 件・一括省略とも受理、省略したキーの警告は出ない | **PASS** (12 項目) | 同 `logs/drop_*.log`, `logs/drop_all.log` |
| `mesh.meshFormat`: 省略と `hdf5` が同じ、不正値 `cgns` は**既定へ置換せず拒否** | **PASS** | 同 `logs/mf_*.log` (`unknown mesh format` で非ゼロ終了) |
| `detectNaNInterval` の解決順序: 両省略 1 / 旧別名のみ 7 / トップのみ 9 / 併記 9 / 0 以下 1 | **PASS** (5 水準とも実効値一致) | 同 `logs/dni_*.log`。起動ログに `[config] detectNaNInterval effective: N` を 1 行足して観測点にした |
| 残置キーの巻き込み無し (`time.deltaT.control` / 両方の `detectNaN` / `speciesImplicitRelax: 0.7`) | **PASS** | 同 `logs/retained.log` (3 件とも読み取りログに出る、S1 警告は出ない) |
| 検査ツールが S1 対象**だけ** WARN、省略時は出ない、削除済み FAIL に混ぜない | **PASS** | `check_solver_config.py` の `STAGED_KEYS`、単体試験 `test_check_solver_config.py` |
| 生成器 3 本から旧キーが消えている | **PASS** | `runner.py` / `runner_wt.py` / `runner_sern.py`。生成した config を PyYAML で読み S1 キー 4 種が無いことを確認 |
| 生成 → 変換 → node solver 起動を通す | **PASS** | `case/26.flat_plate_sst/run_0091_s1_chain_e2e` (S1 キーを外した config で `convertGmshToForge` が exit 0・警告 0、その h5 で forge が 20 step 完走・NaN なし) |

**数値回帰** (§6.5'、同一 IC・300 step・`--boundary` 込み 24 量):

| 比較 | ノイズ床の測り方 | 最大比 | VERDICT |
| --- | --- | ---: | --- |
| 旧バイナリ + 旧入力 → 新バイナリ + 旧入力 | 両側 3 反復 | **1.14** | **PASS** |
| 新バイナリ + 旧入力 → 新バイナリ + **省略入力** | 基準側 3 反復 | **1.14** | **PASS** |

run は `case/26.flat_plate_sst/run_0088`–`0090` (と各 `_rep1`/`_rep2`)。
第 1 陣の SST / DDES / WMLES 一式は再実行していない (今回の変更はそれらの分岐を通らない — plan-4 M4)。
短時間回帰なので収束の証明としては使わない。

**変換器の注意 (この作業中に踏んだ)**: `convertGmshToForge` は引数 2 つ (`入力.msh 出力.h5`) を要る。
1 つしか渡さないと `argv[2]` が null で `basic_string: construction from null is not valid` と abort する。
S1 前のバイナリでも同じなので S1 とは無関係。

### 6.7' 第 2 陣 S2 の受入表 (2026-09-18, codex plan-5 M5)

| 対象 | 必須の判定 |
| --- | --- |
| 移送ツール (`migrate_solver_config.py`) | **`--dry-run` と通常実行が同じ編集・照合を行い、書き込みだけが違う** (移送できない入力を dry-run で見つけられること)。**対象外の変更がゼロ** (編集後 YAML が「元 − 消したキー (+ 移送)」と厳密一致)。**再実行で変更ゼロ** (冪等)。自動移送できない入力を**明示して報告**する |
| 移送ツールの拒否 | anchor / alias / merge key、重複キー、`meshFormat` が `hdf5` 以外 — いずれも**書き換えず手動対応として報告**する。反例を単体試験に入れる |
| 移送ツールの対象範囲 | 既定は**このツリーだけ**。他 worktree は明示したパスでしか触らない。**未知オプションを拒否**する (`--dry-run` の誤記が書き込みにならないこと) |
| 5 キーの拒否 | それぞれを書いた config が **solver / 変換器**の起動時に落ちる。`check_solver_config.py` が **FAIL** で拾う (WARN では投入前に止められない) |
| 移送後の受理 | 移送した config が solver・変換器とも通る。残置キー (`time.deltaT.control`・両方の `detectNaN`・`speciesImplicitRelax`) の実効値が変わらない。`detectNaNInterval` はトップレベルのみになり、**優先順位と下限 1 補正**が従来どおり |
| 雛型と生成器 | `recommended-settings.md` §1 の雛型と runner 3 本の出力に対象 5 キーが**残っていない** (これが S2 へ切り替える条件) |
| 通し確認 | 生産用 `build/` の `forge` と `convertGmshToForge` で「生成 → 変換 → node 起動」を通す |
| 数値回帰 | node 1 ケースの固定 step。比較量は §6.2' の既定 + `--boundary`、**基準側 3 反復**、許容は**ノイズ床 × 2** (§6.2' と同じ数字に揃える)。移送前後の config を同一バイナリ・同一 IC で比較する。短時間回帰から収束・定常性は主張しない |

**S1 と同じく、SST / DDES / WMLES 一式の再実行は不要** (S2 はそれらの分岐を通らない)。

### 6.8' 第 2 陣 S2 の実施結果 (2026-09-18)

| 受入項目 | 結果 | 根拠 |
| --- | --- | --- |
| 5 キーそれぞれが solver の起動時に落ちる | **PASS** (5/5) | `case/26.flat_plate_sst/run_0092_s2_keytest/logs/reject_*.log` |
| 変換器も同じキーを拒否する | **PASS** | 同 `logs/_conv` (`convertGmshToForge` が rc=1 + `no longer supported`) |
| 移送後の config が通る (警告もエラーも無し) | **PASS** | 同 `logs/migrated.log` |
| 残置キーの実効値が変わらない (`time.deltaT.control` / 両方の `detectNaN` / `speciesImplicitRelax`) | **PASS** | 同 `logs/retained.log` |
| `detectNaNInterval` はトップレベルのみ、下限 1 補正は従来どおり | **PASS** (省略 1 / 9 → 9 / 0 → 1) | 同 `logs/dni_*.log` |
| `check_solver_config.py` が **FAIL** で拾う | **PASS** | 5 キーとも FAIL、移送後は出ない |
| 移送ツールが残すキーを壊さない | **PASS** | 単体試験 `tests/unit/test_migrate_solver_config.py` (anchor / 重複キー / 非 hdf5 の拒否、ブロック形の末尾コメント、フロー形、冪等、トップレベル null) |
| 移送の実施と冪等性 | **PASS** | このツリー 3497 本中 **3491 本を移送**、既に清潔 3 本、**手動対応 2 本** (`case/01.cavity/run_heat`, `case/02.airfoil` は `meshFormat: gmsh` の旧 smac 入力)、失敗 0 本。2 回目は変更ゼロ |
| 雛型と生成器に対象キーが残っていない | **PASS** | `recommended-settings.md` §1 の雛型から除去、runner 3 本は S1 で除去済み |
| 生成 → 変換 → node 起動の通し | **PASS** | `case/26.flat_plate_sst/run_0095_s2_chain_e2e` (変換 rc=0、solver 20 step 完走、NaN なし) |
| 数値回帰 (基準側 3 反復、許容 = ノイズ床 × 2、`--boundary` 込み 24 量) | **PASS** (最大比 **1.15**) | 移送後 config (S2 バイナリ) × 3 vs 移送前 config (S1 バイナリ `fbf305cc`)。`run_0093_s2_migrated`(+`_rep1`/`_rep2`) と `run_0094_s2_premigration` |

**副産物の確認**: 非圧縮 SMAC 経路は `main.cpp` に呼び出しが無く (`solvePoisson*` の include はコメントアウト)、
`"smac"` の文字列も solver のソースに存在しない。§5.2 保留表の `physProp.isCompressible: 0 (SMAC 非圧縮)` の
保留理由 (`methods/poisson.md`「コードベースに残っており利用できる」) は**現状と合っていない**。
経路そのものの扱いは本 plan の対象外なので、記録だけ残す。

### 6.3' 生成器・変換器 (M4)

- `design/forge_design/evaluate/runner*.py` の 3 本が生成する `solverConfig.yaml` から削除キーを外し、**生成器を直してから**起動時拒否を入れる。
- `mesh/convertGmshToForge.cpp` は同じ `solverConfig::read()` を使い `cfg.gpu` を確保に渡すので、**solver の GPU 固定と変換器のホスト実行を分けて**扱う。
- `initial` は solver 実行時には使われない (`setInitial` の呼び出しは変換器のみ) ので、**呼出側ごとに必須性を定義**する。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan | `2026-09-17` | [2026-09-17-config-key-pruning-plan.md](../../notes/reviews/2026-09-17-config-key-pruning-plan.md) | **GO-with-changes**, C0/M9/m0 | **全採用 (2026-09-18, §5.2 を改訂)**: M1 `timeIntegration: 1` は 1 段 Euler で 3 (3 段 TVD RK) の別名ではない (`solverConfig.cpp:907/919` で段数・係数が別) → 削除対象から除外; M2 `sstIsotropicStress`/`sstEnergyKSource` の引用先は最終決定と逆 (accepted plan §4.0 は既定 0 + 個別利用可、残作業表も「撤去せず現状維持」) → 保留; M3 「全 155 キー分類済み」が不成立 (棚卸しが `config["mesh"][...]` 形を取りこぼし、同名末端キーを統合していた; 分類表の内訳も 31 と合わない) → **完全修飾パスで再抽出 (188 パス)** し分類をやり直す; M4 共有パーサ・生成器が影響範囲から漏れ (`convertGmshToForge` が同じ `solverConfig::read()` を使い `cfg.gpu` を確保に渡す、design の runner 3 本が `gpu`/`nodeWallDirichlet` を生成) → スコープに追加; M5 `condEquilibrium 1→2` の「固定点同一」は一般には誤り (新しい accepted plan が緩和形の固定点は Δτ 依存と明記、mode 2 は単一凝縮種限定) → 保留; M6 `mesh.axisymMethod: 1` は不採用決定ではなく opt-in 保持で再評価待ち → 保留; M7 A/D/V が限定的な実測から機能廃止へ飛躍 (`lineVisc*` は「このケースでは僅差」、`condTwoTemp` は希薄水/N2 限定、化学 2 キーは反応源の温度評価・停止条件を変える) → 保留に移し、A は「検証した適用範囲」と「廃止する対応範囲」を対応づける; M8 §6 の回帰 3 run が全て node/SLAU/`model: none`/`viscMethod 0` で、SST・粘性壁・line-implicit・cell・陽解法・KEEP の経路を通せない → 経路別の最小回帰表に作り替え; M9 「長時間 run のビット一致」は非退行判定として実行不能 (残差が float atomicAdd で集積、出力専用経路は保存量に出ない、`implicitRelaxSST=-1` は継承指定で定数置換と別物) → 決定的な局所試験はビット一致、CFD は同一バイナリ反復幅 + 事前に定めた物理量許容。保留方針 (`ducrosLimiter`/`condLimiterMode`/`nodeOmegaWfDirichlet`) は同意を得た |
| plan | `2026-09-18` | [2026-09-18-config-key-pruning-plan.md](../../notes/reviews/2026-09-18-config-key-pruning-plan.md) | **GO-with-changes**, C0/M7/m1 | **全採用 (2026-09-18)**: M1 性能スイッチ 3 件 (`blockDPLURDiagCache`/`blockDPLURDqPack`/`primPack`) は元 plan が「不採用で確定、**opt-in 残置**」と決めており、残置決定の読み違いが 3 度目 → **コードごと復元** (§5.2)。M2 棚卸しが入れ子の run config を探索せず「使用 0」を再び誤判定 (826 config が対象外、`forge-perf` が実際に 3 件とも 1 を書いている) → `case/**` を再帰探索し、`<case>/<相対パス>` を識別子に、内容違いの同名は別設定、読み取り失敗は報告 (2796 → **4030 本**)。M3 「live 181 パス」が受理キー一覧になっていない (コメントを走査、節・拒否専用参照を live に混入、必須キーの既定値に節名が入る) → コメント除去 + 節/拒否専用/必須の分離で **値キー 167 / 節 10 / 拒否専用 12 / 起動時拒否 8**。M4 既定値比較の絶対許容差 1e-12 が 10 桁の変更を「既定と同じ」にする + 別名 (`mesh.renumber`→`meshRenumber` 等) を解決できない → 正規化した厳密比較へ、既定値は**代入先メンバ名**で引き、解決できないものは「不明」とし、記載 run 数と非既定 run 数を分離。M5 未知キー検出が末端名照合で誤配置を拾えない → **完全修飾パス**照合に置換 (誤配置は WARN、起動時拒否は FAIL)。実例 (トップレベル `lowMachPrecond` 22 run、`physProp.lowMachPrecond` 10、`condensation.condRealizProject` 5、`space.keepDiss*` 各 4、`physProp.speciesImplicitCoupling` 2) を検出、回帰試験にも追加。M6 正本文書が復元済みキーを「起動時エラー」と案内 → `recommended-settings.md` §9.1 を最終対象に同期し残置 9 件を明記、plan の確定表・残作業表・検証表も同期。M7 非退行の主張が全量では成立せず受入試験も不足 → §6.1'/§6.2' を**事前確定の基準**に書き換え検証を再オープン (§5.1 #5)。m8 `wmlesPrt` は定数化でなく**消費者のない引数の撤去**で、`turbulentPrandtl` への移送案内も誤り → メンバ・カーネル引数ごと削除し案内を撤回。C の段階移行方針は妥当との評価。`isCompressible`/`ro`/`time.last.control` は「読まれない」でなく「パーサは読むが下流の消費者が無い」、`mesh.meshFormat` は省略時 `hdf5` + 不正値拒否を維持、という助言も採用 (§5.3 l) |
| plan | `2026-09-18` | [2026-09-18-config-key-pruning-plan-2.md](../../notes/reviews/2026-09-18-config-key-pruning-plan-2.md) | **GO-with-changes**, C0/M5/m1 | **M1–M4・m6 を採用、M5 は棄却 (2026-09-18)**: M1 §5.3 m が `time.deltaT.detectNaN` を「0 run」としたのは**2 パスの取り違え**で、実際は記載 2879 / 非既定 2781 run (0 run なのは `detectNaNInterval` の方) → 削除決定を撤回し互換読みを維持 (§5.3 m)。M2 `check_field_regress.py` が NaN 入りの候補を PASS にしていた (`max(0.0, NaN)` が `0.0`) → 必須量の欠落・形状不一致・非有限値を**数値判定の前**に検査して終了コード 2、反例を単体試験に追加。M3 §6.2' が必須とした壁量 `Tau_Wall`/`Qw_Wall` が §6.4' の比較から脱落していた (これらは `res_<step>.h5` に出ず、境界出力ファイルにある) → `--boundary` を足して `twall_*`/`qwall`/`utau`/`ypls` を全経路で判定し、4 経路とも再判定して PASS。M4 受理 9 件は分岐到達を証明していない (単成分・`viscMethod 0`・定常) → 「パーサ受理」と明記し、`multispeciesRhoYCommonLimiter` (診断印字)・`passiveFctTolAbs` (残差 4.3e-8→9.7e-22)・`turbulentSchmidt` (トレーサ rel L2 1.4e-1, `Sc_t` との差 5.6e-6) の**到達確認**を追加。m6 棚卸しの「既定のみ」に既定不明・必須が混入 → 必須 30 / 既定不明 5 / 未記載 10 / 既定のみ 3 に分離し、既定不明の非既定数は `null` に。**M5 (cell 回帰が無い) は棄却**: 2026-09-16 のユーザ決定「cell はもう使わない」により検証・回帰は node のみとする。指摘自体は正しい (SST 緩和と WMLES は共有コード) ので、`procedures/verification/README.md` の「node/cell 両方」規則を決定に合わせて書き換え、未検証のまま残る risk を §10 に記載した。**なお M3 の対応中に、ノイズ床を新バイナリ 3 反復だけで測るとカオス的 DDES で 15 倍の過小評価になり偽の不合格を出すことが分かり**、床を両側 (基準側・候補側) の全ペアから測る形に直した |
| plan | `2026-09-18` | [2026-09-18-config-key-pruning-plan-3.md](../../notes/reviews/2026-09-18-config-key-pruning-plan-3.md) | **GO-with-changes**, C0/M4/m1 | **全採用 (2026-09-18)**: M1 `speciesImplicitRelax` は無効キーではない (`speciesTransport_d.cu:835` が読み `scalarTransport_d.cu:290` で化学種更新の増分に掛かる。`relax * (res·dt/V)/fac`) → **第 2 陣から除外**し現状維持。「使用実績が既定値だけ」を「不要」と読んだのが誤りで、**消費者の有無を見ずに D 基準を当てかけた**。M2 `time.deltaT.detectNaNInterval` の別名も無効ではない (`solverConfig.cpp:402` が読み `main.cpp:1900` が検査 step を変える。下位だけに 100 を書いた入力を無視すると間隔が 100 → 1 になり、「読まれていない」という警告文も事実と違う) → S1 では**読みを維持**し「トップレベルへ移せ」と警告するだけ、解決順序 (トップ → 旧別名 → 既定 1) と下限 1 補正を明記、拒否は S2。M3 S1 の作業順序は第 1 陣と**逆** (先に生成器を直すと旧パーサの必須キーが無くて起動できない) → 任意化 → `forge` と `convertGmshToForge` を両方リビルド確認 → 生成器から除去、の順に確定。runner は `solver_density_cuda/build/` を見るので `build-native` だけでは生産チェーンの確認にならない点も明記。`gpu`/`initial` は混ぜない。M4 §6 が第 1 陣用 (起動時拒否) のままで S1 の合格条件が無い + `check_solver_config.py` は既知の値キーを無条件に通すのでパーサ側の警告だけでは WARN 検出を保証できない → **§6.5' に S1 専用の受入表**を新設 (旧入力の受理・省略入力・不正値拒否・別名の解決順序・残置キーの巻き込み検出・検査ツール・生成〜変換〜起動の通し) と、短い固定 step の数値回帰で足りる旨。m5 S2 の条件「実運用ログから警告が消えた」では旧 run の再実行が救えない → **S2 は別マイルストーン**とし、移送確認か互換打ち切りの明示を条件に追加。**第 2 陣は 6 件 → 5 件**に確定 |
| plan | `2026-09-18` | [2026-09-18-config-key-pruning-plan-4.md](../../notes/reviews/2026-09-18-config-key-pruning-plan-4.md) | **GO-with-changes**, C0/M5/m1 | **全採用 (2026-09-18)**: M1 移送ツールが **anchor/alias を共有する入力で残すべき `time.deltaT.control` まで消して照合を通過**していた (`deepcopy` が共有参照を保つので期待値側も一緒に消える) → anchor / alias / merge key を含む入力は**書き換えず手動対応として報告**。M2 **重複キー**を含む入力で照合が solver の実効設定を表さない (PyYAML は後勝ち、yaml-cpp は**先勝ち**。実例 `forge-chem/.../run_0022_dbg_rk3` に `cfl` が 2 つ) → 全階層の重複キーを検出して拒否。M3 引数なしの実行が plan で禁じた**他 worktree まで書き換えて**いた + 手書きの引数処理で `--dry-run` の誤記が書き込みになる → 既定をこのツリーに限定し `argparse` で未知オプションを拒否。M4 `recommended-settings.md` §1 の**現行雛型に対象キーが残って**おり、正本に従う新規 config が S2 で起動できない → 雛型から除去し、これを切替条件にした。M5 §6.7' が実在せず、`--dry-run` が編集・照合の前に終了するので移送不能な入力を見つけられない → **§6.7' を実際に追加**し、dry-run も編集・照合まで同じことをするように直し、許容を §6.2' と同じ「ノイズ床 × 2」に統一。m6 非 `hdf5` の `meshFormat` (このツリーに 2 本) を意味を見ずに削除していた → 自動移送せず手動対応として報告。実施結果は §6.8' |

## 7. 影響範囲

solver の config 読込と分岐、`procedures/` の設定文書、skill `forge-config`、既存 run config (削除キーを書いているものは
再実行時に起動エラーになる → case README に注記)。

## 8. 完了条件

- [x] 分類表 (§5.1 #2) が live な値キー 167 パスを覆い、codex plan レビュー 3 回を通っている
- [x] 第 1 陣の削除 (2 パス + 定数化 5 件 + 撤去 1 件) と旧キーの起動時エラー化
- [x] §6.2' の事前確定基準で §6.4' の検証が全 PASS (4 経路 + 拒否 8 + 受理 9 + 到達 3 + opt-in 3)
- [x] `solver-settings.md` / `recommended-settings.md` / `methods/turbulence` / `procedures/verification` の同期
- [x] 第 2 陣 S1 (対象 5 件の任意化 + 1 回警告) — codex plan レビュー 4 回目を通し、§6.6' で受入 8 項目 + 数値回帰 2 本が全 PASS
- [x] 第 2 陣 S2 (起動時エラー化) — 移送 (このツリー 3491 本) と互換打ち切りの明示の両方で条件を満たし、§6.8' で全 PASS

## 9. 変更ログ

- `2026-09-18` — **第 2 陣 S2 を実装** (§4.2, §6.7', §6.8')。codex plan レビュー 5 回目 **GO-with-changes (C0/M5/m1)** を全採用してから着手。**移送ツールに実際に壊す欠陥が 2 つ**あった: (a) anchor/alias を共有する入力で、残すべき `time.deltaT.control` まで消えるのに照合を通過する (`deepcopy` が共有参照を保つため期待値側も一緒に消える)、(b) 重複キーを含む入力で PyYAML (後勝ち) の照合が solver の yaml-cpp (先勝ち) の実効設定を表さない。どちらも**書き換えず手動対応として報告**するようにし、反例を単体試験に入れた。併せて既定の対象範囲をこのツリーに限定 (引数なしで他 worktree を書き換えていた)、`argparse` で未知オプションを拒否、`--dry-run` も編集・照合まで同じ処理を行うよう修正、非 `hdf5` の `meshFormat` は自動移送しない扱いに。`recommended-settings.md` §1 の雛型からも対象キーを除去した (正本に従う新規 config が起動できなくなるため)。**移送はこのツリー 3497 本中 3491 本、手動対応 2 本 (旧 smac 入力)、失敗 0**、2 回目は変更ゼロ。受入 11 項目 + 数値回帰 (最大比 1.15) が全 PASS。

- `2026-09-18` — **第 2 陣 S1 を実装** (§4.1, §6.6')。codex plan レビュー 4 回目 **GO-with-changes (C0/M4/m1)** を全採用したうえで、対象を **5 件**に絞って着手: `physProp.isCompressible` / `physProp.ro` / `time.last.control` はメンバごと撤去 (読む場所がゼロ)、`mesh.meshFormat` は省略時 `hdf5` (不正値の拒否は維持)、`time.deltaT.detectNaNInterval` は**読みを維持**して「トップレベルへ移せ」と警告するだけ。**`speciesImplicitRelax` は除外** (M1: 化学種更新の増分に掛かる実効キーで、無効ではない)。作業順序は第 1 陣と逆にし、任意化 → `build/` の `forge` と `convertGmshToForge` を両方リビルド確認 → 生成器 3 本から除去、とした (先に生成器を直すと旧パーサの必須キーが無くて起動できない)。`check_solver_config.py` に `STAGED_KEYS` を持たせ (既知の値キーは無条件に通るのでパーサ側の警告だけでは足りない)、起動ログに `[config] detectNaNInterval effective: N` を足して解決順序の観測点にした。**受入 8 項目 + 数値回帰 2 本が全 PASS**。S2 (起動時エラー化) は別マイルストーンで、移送確認か互換打ち切りの明示が条件。

- `2026-09-18` — **167 パスの分類を完了** (§5.1 #2、§5.2「分類の結論」)。素データは [notes/investigations/config-key-inventory-2026-09-18.md](../../notes/investigations/config-key-inventory-2026-09-18.md) (4063 config)。内訳は 非既定の使用あり 119 / 必須 30 / 既定不明 (リスト・マップ値) 5 / 未記載 10 / 既定のみ 3。**「消してよい」と言えるのは第 2 陣の候補 6 件だけ**で、残りは残すか保留。候補は `physProp.isCompressible`・`physProp.ro`・`time.last.control` (いずれもパーサ外に消費者ゼロ。`endTimeControl`/`isCompressible`/`ro` を grep して確認、対照の `dtControl` は `main.cpp:510` ほかで読む)、`mesh.meshFormat` (合法値 1 つ)、`time.deltaT.detectNaNInterval` (0 run の 2 つ目の綴り)、`time.deltaT.speciesImplicitRelax` (10 run 全て既定, 要判断)。**実装は保留** (codex plan-3 の推奨)、着手前に 4 回目の plan レビュー。

- `2026-09-18` — codex plan レビュー **3 回目 GO-with-changes (C0/M5/m1)**。M1–M4・m6 を採用、M5 を棄却 (§6.1)。主なもの: (a) §5.3 m で `time.deltaT.detectNaN` を「0 run」と書いたのは**2 パスの取り違え**で実際は 2781 run が非既定 → 削除決定を撤回。(b) 判定器が **NaN 入りの候補を PASS にしていた** → 不備検査を数値判定の前に入れ、反例を単体試験に追加。(c) §6.2' が必須とした**壁量が比較から脱落**していた (境界出力ファイルにしか無い) → `--boundary` を足して 4 経路を再判定、全 PASS。(d) 受理 9 件は**パーサ受理**であり分岐到達ではない → 3 件に到達確認を追加 (診断印字・残差 14 桁変化・トレーサ場 14% 変化)。(e) この過程で、**ノイズ床を片側 3 反復で測るとカオス的 DDES で 15 倍の過小評価**になり偽の不合格を出すことが判明 → 床を両側の全ペアから測る形に直し、LES/DES は両側 3 本以上を要求することにした。

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

- **cell 側の共有変更が未検証のまま残る** (codex plan-3 M5 を棄却した代償): SST 点陰的緩和 (`update_d.cu:356`) と
  WMLES 壁モデルカーネルは cell/node 共有だが、2026-09-16 の決定により回帰は node のみ。cell を再び使うことに
  なったら、共有変更の cell 検証をその時点でやり直す。

- 「診断・A/B で要る」の線引き: 過去 1 回しか使っていないスイッチを残すか。**案**: 切り分けの再現手順が plan/notes に
  書かれているものだけ残し、それ以外は削除して必要になったら復活させる (git にある)。
- 既定値を変えるか: `lineKFreeze` のように「dual-time では常用推奨」なのに既定 0 のものは、削除ではなく既定変更が適切。
