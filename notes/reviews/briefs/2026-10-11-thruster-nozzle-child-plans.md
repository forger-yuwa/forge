# 諮問ブリーフ: 推力ノズルの子 plan 2 本 (外部領域と剥離評価 / ベルの 2 モード比較) の §4・§6

- 日付: 2026-10-11
- 依頼元: 主セッション (Opus)。AGENTS.md「モデル分担とエスカレーション」条件 1 (§4・§6 を新規に書く) と 3 (参照値・許容差の事前登録)
- 対象: `plans/active/tooling-nozzle-ambient-separation.md` (以下 A)、`plans/active/tooling-nozzle-bell-two-mode.md` (以下 B)。親 `plans/active/tooling-nozzle-design-tool.md`
- 段階: **まだ実装も run もしていない**。§4 (設計方針) と §6 (検証・判定基準) の下書きの妥当性を問う

## 1. ユーザが決めたこと (前提。覆さない)

親 plan §3 確定事項 5〜12 (2026-10-11):

- ③ベル・④デュアルベルは**海面作動を含む**。③④共通で**外気圧を与えて過膨張剥離を評価**できるようにする。
- **全機種 thermally perfect (TP)** を基準 (MOC と forge 評価の両方)。
- ③は「**特性曲線で辿るモード (Rao を MOC で作る)**」と「**形状を直接振るモード**」の 2 つを作って比べる。
- ③の基準ケース: LOX/LH2・燃焼室圧 10 MPa・O/F 6・面積比 ε 50 前後、CEA の燃焼室平衡組成を凍結した TP、**壁は当面断熱**。
- ③の壁は**両モードとも 2 回微分連続 (C2)**。スロートの曲率は 1 つの値 1/R_t。直接モードは拡大部の 5 次 B-spline + 出口壁角 θe を明示の設計変数。特性曲線モードは曲率連続の kernel + Rao 輪郭の位置・壁角の同時当てはめ (`joint_fit_wall`)。初期線は CFD ピン既定・Hall は検算。
- 剥離の第一検証ケース = **DLR TIC test case 1A** (Stark & Hagemann, EUCASS 2007)、2 番目 = NASA MSFC TIC (Ruf ら 2010)。

## 2. 観測事実と出典 (既にあるもの)

- ③の現行評価器 (`design/forge_design/evaluate/runner.py:75-133`): CPG、粘性係数 1.8e-5 定数、断熱壁、領域は出口面で終わり `outlet_statPress` (Ps = p_a)。外部領域が無いので剥離域への巻き込みを解けない。
- Phase 2 (2026-08-13, `plans/accepted/tooling-nozzle-moo-loop.md` の変更ログ): TOP (θn, θe, L) の MOO 77 評価。Rao 80 % のチャート点は同じ長さの最適 η の 0.133 % 下 (判定帯 ±0.15 %)。固定長で θn 24.5〜28.5°・θe 8〜11° が平坦な台地。polish の頭打ち閾値 |Δη| < 3e-4。
- TOP の内部衝撃波 (`notes/investigations/nozzle-top-internal-shock-diagnosis.md`): forge と SU2 (RANS・Euler) で同じ構造。Rd 0.4/0.8/1.5 で内部の ΔM は −0.43/−0.40/−0.40 (接続点の曲率の跳びを小さくしても強さはほぼ不変)。
- TC1A の文献値 (原典 `papers/nozzle_separation/Stark_Hagemann_2007_*.pdf` の表 1、主セッションが本文テキストで確認): 実験 X_sep/R* = 6.6、標準 SST は Alziary −9.1/−9.1/−9.7 %、Nebbache −17.4 %、SA −1.1 % / +4.1 % (SA 改良版 −21.8 %)、SST 改良版 −15.5〜−34.4 %、Baldwin-Lomax −13.2 %、Östlund SST 改良 +3.3 %。マッハディスク X_sh/R* = 9.85 (標準 SST −6.8〜−8.7 %: 調査ノートの記載、主セッションは未確認)。
- `farfield` 境界は軸対称を起動時に拒否 (`plans/accepted/boundary-node-farfield-characteristic.md` §2)。case/27 (cell・旧手法) は上面 `inlet_Pressure` (Pt = p_a)、下流 `outlet_statPress` + 逆流 Pt/Tt、出口端面 `slip` で回り、プルームの非定常で残差 2e-3 で頭打ち (ノズル内部は定常)。
- 既存部品: `geometry/rao_planar.py` (平面の縁条件のみ、67 行)、`geometry/moc_kernel.py` (軸対称・semi-perfect 単位過程)、`geometry/wall_axismach.py:999` `joint_fit_wall`、`feedback/cfd_initial_line.py` (凍結源の形の照合は U→T 5 次 Hermite 前提)、`feedback/deltastar_integral.py` (CONTUR 積分法、粘性係数は空気の Sutherland)。

## 3. 問い

各 plan について、**§4 の設計と §6 の判定基準に Critical / Major があるか**を判定し、GO / GO-with-changes / NO-GO を付けてほしい。指摘は plan ごとに最大 5 件、各指摘の推奨は 1 つに絞る。特に次を見てほしい (これに縛られなくてよい):

- A-1 静止大気を `inlet_Pressure` 上面 + `outlet_statPress` 下流で表す構成 (farfield が使えない軸対称) の危険と、代わりの推奨。出口端面 (リップ) を壁にする扱い。
- A-2 剥離位置の主定義に「forge 自身の全流れ解の壁圧」を使うこと (実験は MOC の全流れ壁圧)。
- A-3 §6 V3 の分類: forge の X_sep を「文献の標準 SST の帯 [5.45, 6.00] R* を数値幅 U で広げた範囲」と比べて forge 固有の誤差の有無を判定する設計が妥当か。帯は別コード・別格子の 4 解から作っている。SU2 SST を同じ格子で**常に**回すべきか (現案は帯から外れたときだけ)。
- A-4 §6 V2 の U ≤ 0.02 X_sep,exp と、格子 1.5 倍の 2 水準で足りるか。V1 の壁圧 3 % 基準。
- B-1 比較の公平さ: (M) は設計変数 L の 1 個で掃引、(D) は 4 変数で MOO。(M) だけ δ* 補正 + 相似倍率の調整を通る。この比べ方で「どちらのモードが良いか」が言えるか。
- B-2 (M) の C2 化: kernel の終点曲率 κ_B を根探しで合わせて Rao 輪郭と C2 でつなぐ案は、問題として成立するか (Rao の輪郭の出だしの曲率は MOC の結果)。代案 (B をまたいで C2 の B-spline を当てはめ、残差を記録) の方が良いか。
- B-3 TP の Rao: 乗数条件が比熱比を陽に含まない形で書けて TP にそのまま使える、という見込みの当否と、確かめ方 (§6 V1)。
- B-4 §6 V2 (Euler と MOC の推力差 0.2 %) と V3 (前線差の閾値 1.5e-3 と U の作り方) の妥当性。
- B-5 着手順: 評価器の TP 化 → (D) → (M) の順でよいか。

## 4. 仮説 (主セッションの見立て。検証されていない)

- 推力系の目的は積分量なので、(M) の非粘性最適性は設計点 η でしか効かず、前線の差は小さい (Phase 2 の台地から 0.1〜0.2 % 程度) と見ている。差が出るのは海面作動 (剥離・内部衝撃波) の側。
- TC1A で forge の SST は文献の SST と同程度 (上流側に −10 % 前後) に出ると見ている。

## 5. 読んでよいもの (これ以外は読まない)

- `sed -n '1,144p' plans/active/tooling-nozzle-ambient-separation.md`
- `sed -n '1,144p' plans/active/tooling-nozzle-bell-two-mode.md`
- `sed -n '37,57p;117,175p;275,290p' plans/active/tooling-nozzle-design-tool.md`
- `sed -n '21,40p;57,62p;159,215p' notes/investigations/nozzle-flow-separation-validation-survey.md`
- `sed -n '1,67p' notes/investigations/nozzle-top-internal-shock-diagnosis.md`
- `sed -n '1,40p' case/27.axi_nozzle_plume/README.md`
- `sed -n '18,32p' plans/accepted/boundary-node-farfield-characteristic.md`
- `sed -n '75,135p' design/forge_design/evaluate/runner.py`
- `sed -n '1,67p' design/forge_design/geometry/rao_planar.py`
- `sed -n '999,1060p' design/forge_design/geometry/wall_axismach.py`
- `sed -n '1,40p' design/forge_design/feedback/deltastar_integral.py`
- `sed -n '312,345p' methods/design/overview.md` (CFD ピンの契約)
- grep は可

## 6. 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックス)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `*.pdf`, `plans/README.md` を読まない** (巨大)。
- 上に列挙した `sed -n` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か上の数値で示す。
