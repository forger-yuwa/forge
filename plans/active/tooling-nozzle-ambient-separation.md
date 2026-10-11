# ノズルの外部領域と過膨張剥離の評価 (③④共通) と、DLR TIC test case 1A による検証

## メタ

- **area**: `tooling / optimization`
- **status**: `draft`  <!-- 2026-10-11 起票。同日 codex diagnose で NO-GO → 指摘を反映した改訂版。外気 BC の方式 (§8-1) はユーザ判断待ち。実装前に --stage plan のレビューを受ける -->
- **related_docs**:
  - [`methods/design/overview.md`](../../methods/design/overview.md) (実装着手前に「外部領域と剥離評価」節を足す — §5 ステップ 0)
  - [`notes/investigations/nozzle-flow-separation-validation-survey.md`](../../notes/investigations/nozzle-flow-separation-validation-survey.md) (検証データの調査。ケース選定の根拠)
  - [`design/CAPABILITIES.md`](../../design/CAPABILITIES.md) §2「高度・背圧」📋・§4「プリューム/チャンバー」🧩 (本計画で上げる)
- **related_plans**:
  - 親: [`tooling-nozzle-design-tool.md`](tooling-nozzle-design-tool.md) §3 確定事項 5〜8・10・12、§5.1 #1・#2
  - 兄弟: [`tooling-nozzle-bell-two-mode.md`](tooling-nozzle-bell-two-mode.md) (本計画の外部領域を海面作動の比較に使う。推力の会計は §4.4 を共通定義にする)
  - 資産元: [`../accepted/boundary-outlet-characteristic-backflow.md`](../accepted/boundary-outlet-characteristic-backflow.md) (出口の逆流は内部の等エントロピーで作る — §3)、[`../accepted/boundary-node-farfield-characteristic.md`](../accepted/boundary-node-farfield-characteristic.md) (`farfield` は**軸対称を起動時に拒否**)
- **created**: `2026-10-11`
- **owner**: `sano` — **実施は別セッション** (2026-10-11 ユーザ決定。設計側の plan はこの機能が使える前提で進める)。引き継ぎ: [`notes/sessions/2026-10-11-handoff-nozzle-ambient-separation.md`](../../notes/sessions/2026-10-11-handoff-nozzle-ambient-separation.md)

## 1. 目的

設計側の軸対称ノズル評価 (③ベル・④デュアルベル) に、ノズル出口より外の外部領域と外気圧 $p_a$ を入れ、過膨張で壁から流れが剥がれる位置を
計算できるようにする。forge の SST がこの剥離位置をどの程度に予測するかを、形状・壁圧・衝撃位置が表で公開された
DLR TIC test case 1A (以下 TC1A) で確かめる。forge 固有の誤差の有無は、同じ格子・同じ BC・同じモデル項の SU2 SST との比較で切り分ける。
完了時には、問題定義 YAML に $p_a$ (作動点のリスト) を書けば、外部領域つきで評価して剥離位置・推力・定常性を台帳に返す状態を得る。

## 2. スコープ

- **やる**
  - 外部領域つきの軸対称構造格子 (ノズル内部 + 出口下流のプルーム領域 + その外側) の決定的な生成と品質ゲート
  - 外気の境界条件 (§4.2。方式は §8-1 のユーザ判断)、段階起動、外気圧を変える継続計算 (NPR を上げる系列・下げる系列)
  - 剥離位置の抽出 (実験と同じ「MOC の全流れ壁圧」を主参照にする定義) と、その定常性判定
  - 閉じた検査体積での推力会計 (§4.4、兄弟 plan と共通)
  - 検証: TC1A (NPR 25.25 の 1 点。剥離位置・壁圧・マッハディスク位置)。forge と SU2 SST を同じ格子で比べる。2 番目の NASA MSFC TIC (NPR 依存) は表の転記が済んでから §6 に参照値を書いて回す
- **やらない**
  - 横荷重・非定常の剥離 (URANS)、3D。RSS とそのヒステリシスは定常 RANS の継続計算で扱える範囲だけ (MSFC PAR は本計画の外)
  - 外部流のある作動 (飛行中の外部流)。本計画は静止大気
  - 乱流モデルの改修 (偏りは記録し、改修は別 plan)
  - デュアルベル固有の部品 (壁圧指定の延長部 MOC)

## 3. 関連 docs と前提

- **現状の評価器**: ③の `evaluate/runner.py` は計算領域がノズル出口面で終わり、出口は `outlet_statPress` ($P_s=p_a$)。剥離域へ外気が巻き込まれる流れと剥離衝撃波を解けない。CPG・`viscMethod: 1` (空気の Sutherland)・SST (既定)・断熱壁 (`runner.py:90,103`)。
- **既存の境界条件の契約 (2026-10-11 にソースで確認)**:
  - `outlet_statPress` の逆流 ($U_n<0$) は**内部の等エントロピーと指定静圧**から状態を作り (`boundaryCond_d.cu:615-660`)、指定した $P_t,T_t$ を逆流の状態に使わない。旧 stagnation 逆流分岐は撤去済み (plan boundary-outlet-characteristic-backflow §2)。下流面から大気の温度・組成を規定することはできない。
  - `inlet_Pressure` は流入方向を境界法線の内向きに固定し、流出を検出すると速度を 0 にする (`boundaryCond_d.cu:1155-1158`)。流出に対しては壁に近い振る舞いになる。
  - `farfield` は軸対称を起動時に拒否する。
- **既存の外部領域の資産**: [`case/27.axi_nozzle_plume`](../../case/27.axi_nozzle_plume/README.md) (Rao M4 ノズル + プルーム、TP N2、cell・旧手法)。出口端面 `slip`、上面 `inlet_Pressure` ($P_t=p_a$)、下流面 `outlet_statPress`。README は「プルームの衝撃セルが非定常で残差 $2\times10^{-3}$ で頭打ち、ノズル内部は定常」と書くが、VERDICT による裏付けは無い (README の解釈)。[`case/23.axi_nozzle`](../../case/23.axi_nozzle/) の 6 ブロック `.geo`。
- **設計側メッシャ**: `meshing/mesh2d.py` は軸対称 1 ブロック (msh4.1 直書き、TFI)。⑤の `meshing/mesh_sern.py` に複数ブロックを節点共有で継ぐ実装がある。
- **推力の会計の手本**: ⑤の `metrics/sern_forces.py:88` は $F=\dot m u_{in}+(p_{in}-p_a)A_{in}+T_{wall}$ で入口の圧力項を含む。
- **壁処理**: SST の壁関数は使わない。低 Re で $y_1^+\lesssim1$ (`check_wall_resolution.py`)。AR ≤ 1000 (構造格子の境界層だけ ≤ 5000 の例外)。
- **ガス**: 全機種 TP 基準 (親 §3-8)。TC1A は乾燥 N2・$T_0$ 283 K なので TP (`thermalMethod: 2`, `species: [N2]`) で回す。
- **文献の基準 (TC1A、Stark & Hagemann EUCASS 2007 の表 1、原典で確認)**: 実験の剥離位置 $X_{sep}/R^*=6.6$ (定義: 最大壁圧勾配の接線と**MOC の全流れ壁圧**の交点)。同条件の RANS 19 解のうち標準 SST は −9.1〜−17.4 %、SA −1.1 % / +4.1 %、k-ω Wilcox +2.7 %・Sarkar −2.4 %。マッハディスク $X_{sh}/R^*=9.85$ (文献の SST の帯は調査ノートの記載で原典未確認)。実験の不確かさは原典に数値が無い。

## 4. 設計方針

### 4.1 計算領域とメッシュ

- **トポロジ**: 軸対称の 3 ブロック (ノズル内部 / 出口下流のプルーム (ノズル半径以内) / その外側)。ブロック間は節点共有。
- **リップ**: 有限厚 $t_{lip}$ の端面を断熱壁にする。リップの外側 (ノズル外壁) は、模型の外形が分かる範囲で壁として入れ、外形が分からない上流側の境界を壁で塞がない (§8-3)。
- **大きさ**: 既定は case/27 と同じ半径 $r_e+5D_e$・下流 $x_e+10D_e$。領域の感度は §6 V2 で測る。
- **格子**: ノズル壁は $y_1^+\lesssim1$。壁法線の層はリップでプルーム側へ滑らかに引き継ぐ (SERN の後縁で格子線の折れが冷点の経路になった)。剥離衝撃・マッハディスクの通る領域 (出口から下流 $2D_e$) は流れ方向に細かくする。
- **系統格子**: 3 水準 (各方向に一定比で細分) を同じ生成器で作れるようにする (§6 V2)。
- **生成**: `meshing/mesh_plume.py` (仮)。`mesh2d` の壁クラスタリング関数を流用し msh4.1 を直書き、node 構成・no-slip 壁 bcond で変換、`check_mesh_quality.py`。

### 4.2 境界条件と初期場

- **外気の境界**: 方式は §8-1 (ユーザ判断)。どちらの方式でも、外側の境界を流入・流出する質量流量と、流入する気体の全温を境界ごとに診断として出力し、台帳に残す。
- **ノズル入口**: `inlet_Pressure` (燃焼室 $P_t,T_t$。TC1A: 2.5 MPa / 283 K)。乱流の入口条件は §8-4。
- **小さい並行流**: 静止大気を外周の小さい並行流で近似する案は、安定化の手段として扱わず、結果を変えうる別の境界条件として扱う (採るなら §6 で感度を測る)。
- **初期場と段階起動**: [`procedures/divergence-and-startup.md`](../../procedures/divergence-and-startup.md) に従う。ノズル内部は高 NPR の全流れ解から始め、外部は $p_a$ の静止状態。NPR を下げる系列は外気圧を段階的に上げ、上げる系列は低い NPR の収束解から再始動する。同一格子の再始動は `restart_field.py`。
- 数値設定は [`procedures/recommended-settings.md`](../../procedures/recommended-settings.md) (skill `forge-config`)。SST のモデル項は §6 V3 の比較のため明示して固定する (§4.6)。

### 4.3 剥離位置の抽出 (事前登録する手順)

- **主参照の全流れ壁圧 = MOC**: 実験と同じく、TC1A の輪郭に沿った前進 MOC (`geometry/moc_kernel.py`、TP N2、スロート初期線は兄弟 plan と同じ CFD ピン、Hall は検算) の非粘性壁圧を使う。原典の図 3c の MOC 壁圧と突き合わせ、差を記録する。
- **副参照**: forge 自身の全流れ解 (NPR を十分高くした run) の壁圧。主参照との差を台帳に残す。
- **抽出手順**: 壁圧を実験のポート位置 (2.5 mm 間隔) に再標本化 → 最大勾配の点を中心差分で求め、その点の接線を引く → 主参照との交点を、最大勾配の点から上流 $2R^*$ の区間で探す。forge 格子の節点上の壁圧で同じ手順も行い、ポート再標本化との差を「抽出の幅」として位置の単位 ($R^*$) で記録する。
- **副定義**: $\tau_{w,x}$ が最初に負になる位置。マッハディスクは軸上静圧の最大勾配の位置。剥離域の最大壁圧 $p_{b,max}/p_a$ も記録。
- 実装は `metrics/separation.py` (仮)。合成壁圧で接線交点を解析値どおりに出す単体試験を置く。

### 4.4 推力の会計 (兄弟 plan と共通の定義)

ノズル内壁・入口面・出口面で閉じた検査体積 CV を 1 つ決め、同じ CV から 2 つの会計を作る。

- **壁会計**: 入口面の (運動量流束 $+\,(p-p_a)$ の面積積分 + 粘性応力) + 内壁の $(p-p_a)$ と摩擦の軸方向積分。
- **出口面会計**: 出口面の (運動量流束 $+\,(p-p_a)$) の符号付き積分 (逆流を含む)。
- 2 つの差を CV の運動量収支の残差とし、**付着流でも剥離流でも**受理基準 (§6 V1・V3) を置く。
- **総推力と流れ推力を分ける**: リップの端面・外壁にかかる力を含めたものを総推力、ノズル内部 CV だけのものを流れ推力とし、両方を記録する。

### 4.5 定常性

定常計算は局所擬似時間で進むので、`check_quasisteady.py` が `OSCILLATING` を返しても実験の平均±振幅に対応させない。剥離位置・マッハディスク位置・推力が `STEADY` でなければ、その run は**判定不能**とする (本計画では非定常の統計を扱わない)。`check_convergence.py` の VERDICT は区間付きで併記する。

### 4.6 SST のモデル項

TC1A の基準比較は**文献の「標準 SST」に近い項** (`dilatationCorrection: 0`、`katoLaunder: 0`。他の `turbulence.sst*` キーは既定) で行い、SU2 の SST をこれに合わせる。現行の推奨レシピ (`dilatationCorrection: 2`、`katoLaunder: 1`、recommended-settings §2) は感度として別に回す。SU2 側の対応する項 (圧縮性補正・生成項の形) は §5.1 #7 で照合表を作ってから決める。

## 5. 実装ステップ

0. `methods/design/overview.md` に「外部領域と剥離評価」節を書く (§4 が plan 段レビューを通った後、実装前)
1. TC1A の表 (原典の表 3: 輪郭約 150 点・壁圧・衝撃位置) と図 3c (MOC 壁圧) を転記し、図と突き合わせる (`case/67.nozzle_separation_validation/data/`、新設)
2. 外気 BC (§8-1 の決定に従う)
3. `meshing/mesh_plume.py`: 3 ブロック・系統格子 3 水準、品質ゲート、決定性の試験
4. runner: 外気 BC の生成、段階起動、NPR の継続計算、境界の流入出の診断
5. `metrics/separation.py`: 剥離位置 (主・副)、マッハディスク、$p_{b,max}$、推力の 2 会計と CV 残差
6. TC1A の MOC 全流れ壁圧 (§4.3)
7. TC1A の run (全流れ参照 → NPR 25.25、3 格子・領域感度) と SU2 SST (同じ格子・BC・モデル項)
8. 判定 (§6)、上位への諮問、codex result 段

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | ~~§4・§6 を上位に諮る~~ | **判断: 2026-10-11 codex (diagnose) [記録](../../notes/reviews/2026-10-11-thruster-nozzle-child-plans-diagnose.md) — NO-GO。C1・A1〜A4 を採用して §3・§4・§6 を改訂** (採否の詳細は §6.1) | F |
| 1b | **外気 BC の方式を決める** | §8-1。ユーザ判断 → 決まったら codex `--stage plan` | F |
| 2 | TC1A の表の転記と検査 | 原典 PDF (`papers/nozzle_separation/Stark_Hagemann_2007_*.pdf`、git 管理外) の表 3 と図 3c。調査役の指摘 (重複点 1 点、反射衝撃と滑り線が同値) を原典の図と突き合わせて扱いを決める。合格: 輪郭 CSV を図に重ねて一致、点数と範囲 ($X/R^*$ −4.5〜9.0) が原典と一致 | O |
| 3 | 外部領域のメッシャ | §4.1。合格: `check_mesh_quality.py` PASS、同一入力でビット同一、全流れ解で壁の $y_1^+\le1$ (`check_wall_resolution.py`) | O |
| 4 | runner と剥離抽出・推力会計 | §4.2〜§4.4。合格: 合成壁圧の単体試験、全流れ解で CV 残差 ≤ 0.1 % | O |
| 5 | TC1A の MOC 全流れ壁圧 | §4.3。合格: 原典図 3c との差を記録 (合否は §6 V1 の CPG/TP 差の確認後に決める) | O |
| 6 | TC1A の run と判定 | §6 V1〜V3。run は case/67 に置き README の run 一覧を維持 | O (実行) / F (解釈) |
| 7 | SU2 SST の照合 | 同じ格子・BC で SU2 SST を回す。先に forge と SU2 の SST の項の照合表 (圧縮性補正・生成項・壁の ω 境界・入口の k・ω) を作る。[`procedures/su2-cross-check.md`](../../procedures/su2-cross-check.md) | O |
| 8 | MSFC TIC の参照値・許容差の登録 | 付録 A の壁圧表 (スキャン) の転記と、MOC での輪郭の再構成 (スロート直後の段差の扱いを含む)。§6 に追記してから run | F |

## 6. 検証

**参照値と許容差は run の前にここで固定する** (エスカレーション条件 3)。codex diagnose (2026-10-11) の指摘を反映した案で、§8-1 の決定後に plan 段レビューを受けて確定する。

- **V0 (単体)**: メッシャの決定性・品質ゲート、剥離抽出の合成データ試験、推力会計の検算。
- **V1 (全流れ参照)**: TC1A の形状で NPR を十分高くした全流れ run。ノズル内で剥離しない (全壁で $\tau_{w,x}>0$)、CV 残差 ≤ 0.1 %、`STEADY`。実験の剥離より上流 ($X/R^*$ 1〜5) の壁圧が実験値と各ポートで相対 3 % 以内 (上流の場の整合確認であり、剥離位置の精度を保証するものではない)。MOC の主参照と forge 全流れ解の壁圧の差を、交点付近 ($X/R^*$ 5〜7) で記録する。
- **V2 (数値の幅)**: NPR 25.25 で次を測り、**保守的に足し合わせて** $U$ とする: 系統格子 3 水準の収束傾向から推定した格子の幅、領域 (半径・下流をそれぞれ 2 倍) の変化、反復の幅 (判定区間の前半と後半の差)、抽出の幅 (§4.3)。$U\le0.13R^*$ ($0.02\,X_{sep,exp}$) を目標とする (実験の不確かさではない)。3 水準が単調に収束しなければ判定不能。
- **V3 (forge 固有の誤差の切り分け)**: 同じ格子・同じ BC・照合したモデル項 (§4.6) で forge と SU2 の $X_{sep}$ を比べる。$|X_{sep,forge}-X_{sep,SU2}|\le U_{forge}+U_{SU2}$ なら **「forge 固有の誤差は $U$ の範囲で見えない」**、超えれば forge 固有の誤差の疑いとして、原因を書く前に上位に諮る。マッハディスク位置と $p_{b,max}/p_a$ も同じ比較を記録する。
  - 実験との差と文献の帯 (標準 SST −9.1〜−17.4 %) は**参考**として併記するだけで、判定には使わない。
  - 推奨レシピ (§4.6) との差は感度として記録する。
- **V4 (MSFC TIC)**: §5.1 #8 で参照値・許容差を書いてから回す (未登録)。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| diagnose | `2026-10-11` | [`2026-10-11-thruster-nozzle-child-plans-diagnose.md`](../../notes/reviews/2026-10-11-thruster-nozzle-child-plans-diagnose.md) (兄弟 plan と一括) | NO-GO, C1/M4 (本 plan 分: 共通 C1・A1〜A4) | **C1 採用** (入口の圧力項の欠落を確認 — `sern_forces.py:88` は含む。§4.4 を閉じた CV の 2 会計に)。**A1 採用 (契約の誤認を訂正)** — `boundaryCond_d.cu:615-660` と accepted plan で確認したので、codex 提案の境界カーネル A/B (逆流が $T_t$ に応答するか) はソースと accepted plan の記述で答えが出ており省略。対処の方式は §8-1 (ユーザ判断)。**A2 採用** (主参照を MOC に、抽出手順を事前登録)。**A3 採用** (SU2 SST を同一格子で必須、文献帯は参考に)。**A4 採用** (3 水準・保守的合成、定常計算の OSCILLATING は判定不能)。前提への異議 (case/27 の README は VERDICT の裏付けなし) も採用し §3 に明記 |

## 7. 影響範囲

- 新規: `design/forge_design/meshing/mesh_plume.py` (仮)、`design/forge_design/metrics/separation.py` (仮)、`case/67.nozzle_separation_validation/` (新設、README の run 一覧)
- 変更: `evaluate/runner.py` (外気 BC・作動点・継続計算)、`probdef.py` (作動点 $p_a$ のリスト)、`design/CAPABILITIES.md`、`methods/design/overview.md`
- ソルバ: §8-1 で (a) を選ぶ場合は `farfield` の軸対称化を別 plan で行う (数値カーネルの変更なので上位に諮る)。(b) なら変更しない

## 8. 未確定事項

1. **外気の境界の方式** (codex A1、ユーザ判断):
   - (a) 軸対称で使える双方向の開放境界 (流入は大気の温度・組成、流出は内部状態) を先に作る。`farfield` の軸対称化を別 plan で検証する。拒否フラグを外すだけにはしない。
   - (b) 現行の `inlet_Pressure` (上面) + `outlet_statPress` (下流面) のまま、境界の流入出の診断と領域感度 (§6 V2) を「境界が結果を変えていない」ことの受入れ条件にして先に進む。受入れ条件を満たさなければ (a) へ。
2. TC1A の周囲温度・外気組成 (乾燥空気と仮定するか) と、リップの厚み・外形 (原典で読めるか)
3. リップ外側の上流境界をどこで切るか (外形の分かる範囲で壁、その先は開放境界)
4. ノズル入口の乱流条件 ($k$, $\omega$ または乱流強度・粘性比)。剥離位置への感度を 1 組だけ測る
5. TC1A の MOC 全流れ壁圧を自前で作るか、原典図 3c の読み取りで済ませるか (自前を主、読み取りを照合とする案)

## 9. 変更ログ

- `2026-10-11` — codex diagnose (NO-GO) を反映して改訂: 推力会計を閉じた CV の 2 会計に (入口圧力項の追加、総推力と流れ推力の区別)、既存 BC の契約を訂正 (§3) し外気 BC の方式を §8-1 に、剥離抽出の主参照を MOC にして手順を事前登録、SU2 SST の同一格子比較を必須に (文献帯は参考)、系統格子 3 水準と保守的合成、定常計算の OSCILLATING は判定不能。
- `2026-10-11` — 起票 (親計画 §5.1 #1・#2)。§4・§6 は上位に諮る前の下書き。
