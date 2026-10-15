# 諮問ブリーフ: M_knot 高速化 (逆 MOC の JIT と充填以外の削減) の result 段の解釈 (2026-10-15)

plan: `plans/active/tooling-nozzle-mknot-fast.md` (§4.4・§4.5・§4.6、§6 V3〜V5、§6.3・§6.4)。AGENTS.md のエスカレーション条件 7
(result 段の解釈を確定する前)。依頼者: forge-7d (Opus)。

## 問い (1 つ)

**「numba の版を既定にして、V5 (許容差の照合) と §6.4 の所要の記録で result 段に進んでよいか」**。進めない理由があれば、最大の残る危険を 1 つに
絞り、それを切り分ける試験を 1 つだけ設計してほしい (変えるのは 1 箇所、判定は下の指標)。下の (a)〜(c) は候補で、縛られなくてよい。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。`git show` / `git diff` は下に書いた範囲だけ。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- plan: `sed -n '106,170p' plans/active/tooling-nozzle-mknot-fast.md` (§4.4〜§4.6)、`sed -n '191,300p'` (§6、§6.2〜§6.4、§6.1)。
- 単位過程の numba の版: `sed -n '1,315p' design/forge_design/geometry/moc_fill_nb.py`
  (ν→M `_mach` 44–78、`_pass` 116–143、`_resid` 146–165、`pairs_nb` 168–315)。
- 網の検査の走査: `sed -n '316,416p' design/forge_design/geometry/moc_fill_nb.py`、C⁺ 流束閉包: `sed -n '417,500p'` 同ファイル。
- 呼び出し側: `sed -n '30,125p' design/forge_design/geometry/moc_inverse.py` (backend の選択・ガスの振り分け・差し替えの検出)、
  `sed -n '160,245p'` (集計 `MocFillDiag.add`)、`sed -n '295,460p'` (`_kernel`・`_pairs`・`_pairs_nb`・`fill_arrays`・`fill_levels`)、
  `sed -n '591,625p'` (`cplus_flux_wall` の振り分け)。
- numpy の版 (参照): `sed -n '400,676p' design/forge_design/geometry/moc_kernel.py` (`pm_mach_vec`・`_interior_pass`・`interior_vec`)。
- 網の検査: `sed -n '73,215p' design/forge_design/geometry/moc_diagnostics.py`。
- 試験: `sed -n '1,40p' design/tests/run_moc_jit_tests.py` (説明)、`sed -n '120,200p'` (照合の関数)、`sed -n '358,470p'` (部品・網・異常系)。
- 差分の全体を見るなら `git diff 9468e0d0 cc3b800f --stat` だけ (中身は上の範囲で)。

## 観測事実

- 変更の流れ (feature/mknot-fast): 9468e0d0 (引き継ぎ時点、numpy の版) → aeaf2216 (numba のビット同一の版) → 57f73940 (ユーザ決定でビット同一を
  やめて簡素化: numba は 1 段の単位過程 `pairs_nb` だけ、段のループと集計は Python で共通) → f7e185c4 (充填以外の削減: 集計・網の検査・
  `cplus_flux_wall`・単位過程の超越関数の削減) → cc3b800f (並列の見積もりの表)。
- **数値の照合 (V5、numba 対 `FORGE_MOC_KERNEL=numpy`)**:
  - 合成 (`run_moc_jit_tests.py`、ALL OK): 4 つの組 (legacy/analytic × fixed2/converge) × `fill_levels`・`fill_arrays` × 放射源流 n_ax 60・300、
    semi-perfect の表 (case/71 の燃焼ガス)、異常系 (欠損・反復中の非有限・平行・軸より下・上限 50 回・θ_r の欠如) で、NaN の位置・分類・分岐・
    `moc_gate` が同じ、位置の差 ≤ 2.4e−14 r_t、θ/ν ≤ 1.8e−14 rad。収束した対の修正子の回数が ±1 ずれた対が 8 対 (更新量が収束の許容差
    1e−12 のちょうど境目; 試験は収束した対だけ ±1 を許す)。最終残差の最大になる段・対は丸めで入れ替わるので比べていない (値は ≤ 1e−12 で照合)。
  - 実データ (`--real`、ALL OK、f7e185c4): M5 case/70 c2b (M_K 2.05)・M5.5 case/71 (r_t 0.1 m、M_K 2.80)・M6 case/45 prod (M_K 2.7、CFD の
    初期線) × N 600 (L_total_tol 0.1 r_t)・2400。差の最大は網の位置 2.0e−12 r_t・θ/ν 1.8e−14 rad・壁 1.6e−12 r_t・θ_max 1.2e−13°。
    `evaluate_mknot` の判定の欄 (ok・reasons・error)・集計の分類の数・`moc_gate` の pass は全部同じ (M6 は両方 ok=False)。
  - 受入 (事前登録、§6 V5): 位置 ≤ 1e−10 r_t、θ/ν ≤ 1e−10 rad、NaN の位置が同じ、壁 ≤ 1e−10 r_t、θ_max ≤ 1e−8°、判定の欄が同じ、分類の数が同じ。
- **出力が完全一致の変更** (参照と比べる試験つき): 集計 `MocFillDiag.add` (以前の書き方の写しを参照、42 通り)、網の検査 (git 57f73940 の版を
  参照、12 通り、反転・点の堆積を人工的に入れた網を含む; 点の堆積の数 100・79 対、seam 除外 9・89・78・62 対が一致。位置の一覧 (記録だけ、
  先頭 20 件) は並び順を段・添字の順に変えた — 20 件を超えると選ばれる 20 件が以前と変わる)、`cplus_flux_wall` (CPG γ 1.4・1.25・semi-perfect ×
  m* 3 通り、case/70 の網)。
- **既存の design の試験** (19 本) は全部通る (`run_moc_forward_tests.py` の「HEAD から無変更」の検査は commit 後に通る)。`run_moc_axis_limit_tests.py`
  の 8 (実装前 170b0d75 とのビット同一) は新しい側を `FORGE_MOC_KERNEL=numpy` で回すよう直した。`run_sern_moc_tests.py` は HEAD (変更前) でも
  `moc_sern.design_ramp` の ValueError で落ちる (moc_kernel だけを使う経路、今回の変更と無関係)。
- **所要 (AWS B、AMD EPYC 7R32、BLAS 1 スレッド、3 回、±0.3 %)**: N 2400 の 1 評価は引き継ぎ時点 (a) 21.4 s・(b) 8.52 s → (a) 9.17 s・(b) 3.58 s
  (2.34・2.38 倍)。V5 の当初の「3 倍」には届かず、**ユーザ決定 (2026-10-15)「今の速さで受入」で、結果を見た後に V5 の所要の条件を改訂した**。
  カード 1 枚 (V4、4 並列): case/70 の M_knot の所要 188 → 81 s、選ぶ値は同じ (2.05)。case/71 は M_knot ≈ 55 s、選ぶ値 2.80 (V3 と同じ) の後、
  格子の見積もり (節点数 58.2 万 > 上限 30 万) で止まる (新旧とも、既知)。
- **1 評価のメモリのピーク** (AWS B、新しいプロセス): N 600 0.15 → 0.26 GB、N 2400 1.08 → 1.30 GB。並列の見積もりの表を {600: 0.27, 2400: 1.31} に。
- 初回のコンパイル (キャッシュの無い機械で 1 回だけ、`cache=True`): 計 約 10 s。

## 期待値と出典

- 単位過程の式は numpy の版 (`moc_kernel._interior_pass`・`_interior_residual`) と同じ。違いは丸め: ν→M の Newton を対ごとに止めて前の反復の
  M から始める (numpy の版は配列全体の収束で止め、毎回 1 + max(ν, 1e−6) から)、sin μ = 1/M (numpy は sin(arcsin(1/M)))、tan を sin/cos から、
  表の補間は素直な二分探索と線形の式 (numpy の `np.interp` と傾きの式の丸めが違いうる)。
- 生産の既定は analytic + converge (`runner_axismach.py:238`)、`wall_mode: cplus` (`fill_levels`)。生産の逆 MOC に渡るガスは float の γ と
  `GasSemiPerfect` だけ (`runner_axismach.py:425`)。

## 実施済みの操作と結果 (潰した候補と証拠)

- Newton の停止を |ΔM| < 1e−8 に緩める案: M → 1 の近くで |ΔM| が 4.1e−9 に達した (ν ∝ (M−1)^{3/2} で 2 次収束の見積もりが崩れる) → 捨てた。
- `cplus_flux_wall` を numpy の配列演算で列にまとめる案: 完全一致でも速くならなかった → 捨てた (numba の版にした)。
- `joint_fit_wall` (有効制約法) の 2.2 s は、裏の試験と BLAS のスレッドを取り合った測定の誤り。負荷なし・BLAS 1 スレッドで 0.21 s → 対象外。
- 平行な特性線の合成の入力で位置の差が 0.25 出た件: 隣の対が θ+μ = 90° の向き (tan が発散) で丸めの差が増幅した入力の作り方の問題 →
  平行な対だけの前線に直した (生産の網では θ+μ ≈ 90° は起きないと考えているが、確かめてはいない)。

## 仮説・判断が要る点 (候補、縛られなくてよい)

- (a) 丸めの違いで判定が反転する危険: 修正子の回数 ±1・最終残差の最大の位置の入れ替わりは記録だけにした。`evaluate_mknot` の判定の欄が
  境目の近くで反転する入力 (例: 成立の境目の M_K、網の検査の `min_sin_angle` の閾値の近く) を V5 で試していない。
- (b) semi-perfect の表の補間を `np.interp` と違う丸めで書いた (単位過程の ν→M)。表の端 (ν が表の範囲外) の扱いは端の値で、numpy と同じ。
- (c) V3 (b) ③ (M5.5 の寸法ゲート) は、カードの生成が M_knot の後の格子の見積もりで止まるので判定できないまま。格子の見積もりは別 plan
  (`tooling-nozzle-card-local-mesh`、M6ノズルの設計、codex の NO-GO で書き直し待ち) の範囲。この plan を、③ をその plan に移して閉じてよいか、
  節点の上限を明示で許して (ユーザの了承が要る) ③ だけ判定してから閉じるべきか。
