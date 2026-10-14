# 諮問ブリーフ: ③ベル (M) モードの S1 — 凍結 TP の軸対称 Rao 条件の導出と構成手順

- 日付: 2026-10-11
- 依頼元: 主セッション (Opus)。plan `plans/active/tooling-nozzle-bell-two-mode.md` §5.1 #4 (担当 F = 上位の判断) と §4.4・§6 V1。このセッションの諮問は diagnostician (Fable) — ユーザ指示
- 作業ツリー: `/home/sano/work/forge-thruster` (HEAD `6e1d7174`)

## 1. 依頼の内容

(M) モード = 特性曲線で辿る Rao の最大推力輪郭を、凍結 TP のガスで作る。その数学の部分 (S1) を固めたい。実装は後で implementer が行うので、ここでは**導出と構成手順と受理の試験**を、実装者が迷わない粒度で書いてほしい (ファイルは変更しない — 返答の本文に書く。主セッションが `methods/design/overview.md` と plan に写す)。

## 2. 前提 (決定済み、覆さない)

- 両モード提供 (親 plan §3-16)。(M) は S1 を通るまで比較に入れない。S1 が上位への諮問 2 回で閉じなければ TIC で代替 (§8-4)。
- ガス: 凍結 TP。CEA の燃焼室平衡組成を凍結した擬似種 EXH (`gas/rocket_exhaust.py` の `RocketExhaust`、MOC 用は `moc_gas` → `GasSemiPerfect`)。基準ケース LOX/LH2・10 MPa・O/F 6、$T_c$ 3523.96 K、γ* 1.198。
- 壁は C2: スロートの曲率 $1/R_t$ ($R_t$ = 1.0 $r_t$ 仮)、収縮部は U→T 5 次 Hermite。kernel (スロートから B までの膨張部) は曲率が連続に変わる曲線族で、B での曲率 $\kappa_B$ をパラメータに持つ。古典の円弧 kernel も回帰用に選べること。
- 与えるのは $(L,\varepsilon)$ (物理壁。δ* 補正後の拘束付き当てはめで合わせる — §4.2・§4.4)。初期線は Hall (n ≤ 321) または CFD ピン。前進 MOC は共用モジュール `geometry/moc_forward.py` (配列版 `interior_vec`、analytic・converge)。
- C2 接続: $g(\kappa_B)=\kappa_{Rao}(B^+;\kappa_B)-\kappa_B$ の根 (各評価で B・D を解き直す、枝の追跡、複数根の選択規則) — plan §4.4。

## 3. 問い (それぞれ具体的に)

1. **変分問題の定式化と最適性の条件**: 制御面 DE (出口リップ E から出る特性線) 上で、質量流量・長さ・出口半径 (固定 ε) の拘束の下で推力 $\int(p-p_a+\rho u^2)\,2\pi r\,dr$ (適切な形) を最大化する Lagrange 問題を、凍結 TP の $h(T)$・$s(T,p)$ (等エントロピー・等全エンタルピー) で書き、DE 上の条件と縁 E の条件を導く。Cain 2010 式 4.1–4.3 (軸対称、CPG で書かれている) と Rao 1958 の古典条件 (例: 制御面上で $V\cos(\theta-\mu)/\cos\mu$ と $r\rho V^2\sin^2\theta\tan\mu$ が一定) の TP での形はどうなるか。比熱比を陽に含まない形に書けるか (書けないなら何が残るか)。固定 ε が端点の変分をどう制限するか (E の位置が固定される)。
2. **構成手順**: 与えた $(L,\varepsilon)$ と kernel の族から、B・D を求め、領域 BDE を MOC で埋め、壁流線 BE を質量流量で取る、までのアルゴリズム。未知数と方程式の数、根探しの順序 (外側に $\kappa_B$、内側に B・D)、初期推定と枝の追跡、解が無い・複数ある場合の扱い。kernel の曲線族の具体形 (例: 壁角 θ(s) の多項式で、スロートで曲率 $1/R_t$、B で $\theta_B$ と $\kappa_B$)。
3. **受理の試験 (§6 V1 の具体化)**: CPG の極限で古典 Rao を再現する試験の参照値 (Rao チャートは γ 1.2 前後で読み取り誤差 ±1° — 他に式レベルで照合できる参照があれば)、TP の停留性の試験 (拘束を保つ壁の摂動族の作り方、一次変分が 0・振幅の掃引で $O(\epsilon^2)$・二次変分の符号)、保存の試験 (DE の流量 = スロート流量、推力の 2 会計)。
4. **成立の見込み**: この S1 は閉じる見込みが高いか。危ないところ (例: 固定 ε と長さ固定の同時拘束、TP での特性線の関係式、kernel 終点の曲率の根の存在) と、閉じない場合に TIC へ切り替える判断の目安。

推奨は 1 つに絞り、式は KaTeX ではなくプレーンな記法か Unicode でよい (主セッションが写す)。根拠は文献 (Cain 2010・Shyne & Keith 1990・Rao 1958・Zucrow & Hoffman) とリポジトリのファイル行で。

## 4. 読んでよいもの

- `plans/active/tooling-nozzle-bell-two-mode.md` (§4.2〜§4.4、§6 V1、§8)
- `design/forge_design/geometry/rao_planar.py` (平面の縁条件、67 行)、`design/forge_design/geometry/moc_forward.py` (`forward_moc` の入出力)、`design/forge_design/geometry/moc_kernel.py` (565 行付近 `interior_vec`、軸対称の適合式)、`design/forge_design/gas/semiperfect.py` (`GasSemiPerfect` の `nu`・`mach_of_nu`)、`design/forge_design/gas/rocket_exhaust.py` (`RocketExhaust`)
- `plans/active/tooling-nozzle-sern-chain.md` §4.3 (平面版の導出の経緯、Shyne の式 11–15 の扱い)
- `papers/nozzle_design/` の Cain 2010・Shyne & Keith 1990 の PDF は**主作業ツリー** `/home/sano/work/forge/papers/nozzle_design/` にある (git 管理外)。PDF を読むなら必要なページだけ (全体を読まない)
- 禁止: `*.log`, `res_*.h5`, `plans/README.md`。ファイルは変更しない。forge は起動しない
