forge (自作の圧縮性 FVM ソルバ。CUDA/float32、cell 中心と node 中心 median-dual の 2 離散化、現在は node 主体。
SLAU/Roe/KEEP、block-DPLUR 陰解法、SST、多成分 TP、凝縮、軸対称、ノズル設計ツール design/forge_design を含む) の
リポジトリに対する**外部レビュー**を依頼する。忖度なしで、主張はコードと実測 (run の数値) で検証すること。
結論が「この計画/結果は誤り」でも構わない。両論併記で逃げず、推奨は 1 つに絞ること。

ルール:
- **ファイルを変更しない** (read-only サンドボックスで動いている。読む・実行して確認するのは可)。
- 出力は日本語。識別子・ファイル名は原語のまま。
- 指摘は **Critical / Major / Minor** の重大度付きで、必ず根拠 (`ファイル:行` または `run_*` の数値) と対案をセットで書く。
- リポジトリのルールは `AGENTS.md`、現在仕様は `methods/`、運用手順は `procedures/`、設計判断は `plans/`。
  用語や設定の意味は推測せず `procedures/solver-settings.md` / `procedures/recommended-settings.md` を読むこと。
- 収束の判定は `solver_density_cuda/tools/check_convergence.py <run_dir>` (各 run の `CONVERGENCE_VERDICT.txt`)、
  派生量の定常性は `check_quasisteady.py` の VERDICT を根拠にする。`rms_ro` 単独やスナップショット 1 枚で判断しない。

## 依頼: 計画立案時レビュー (stage = plan)

対象の plan は下に全文を貼る (`plans/active/time_integration-line-viscous-jacobian-dt.md`)。これから実装に入る前の段階なので、次を順に評価せよ。

1. **目的とスコープ** (§1, §2): 解こうとしている課題は正しく同定されているか。既に解決済み/別 plan と重複していないか
   (`plans/README.md` と `plans/accepted/` を確認)。
2. **設計方針** (§4): 数学的・数値的に健全か。forge の既存構造 (`methods/architecture/overview.md`、該当 `methods/<area>/`)
   と整合するか。node/cell 両離散化、float32、陰解法 (block-DPLUR)、周期・軸対称などの既知の落とし穴に抵触しないか。
   代替案と比べて費用対効果は妥当か。
3. **実装ステップと残作業表** (§5, §5.1): 順序・粒度は妥当か。抜けている前提 (メッシュ品質、IC、段階起動) はないか。
4. **検証計画** (§6): 判定基準は定量的か。検証ケースの選択は `procedures/verification/README.md` と整合するか。
   「収束」「一致」を何で判定するかが書かれているか。
5. **見落としているリスク**: 我々が気づいていない構造的問題があれば挙げよ。

最後に「この計画で実装に進んでよいか」を **GO / GO-with-changes / NO-GO** の 1 語で判定し、
GO-with-changes なら実装前に直すべき点を優先順で列挙すること。

## 重点

§6.3 (2026-10-11 に追加した、上限 50 での全部入り L7・生産 L0・粘性入り L5 の 3 本を、速度 plan §6.15 の水準に届く step 数で比べる事前登録) と、台本 lvcmp.sh・事前のゲート lvcmp_pregate.py・本判定 lvcmp_judge.py (作業ツリーの未 commit の新規ファイル) だけを点検する。§6 (dt の A/B) の lvcdt_*.py と faceh の lvc75_*.py を元にし、見張りは速度 plan の m9_watch.py (--phase line) をそのまま使う。V(b) の運動量の薄層の項の独立計算、比較の閾値 15000 step・各腕 1 本という設計、分岐と言える範囲が妥当か、抜けている確認は何か。設計の諮問は diagnostician (notes/reviews/2026-10-11-lvcdt-result-diagnose.md の Q2)。

## plan 全文 (`plans/active/time_integration-line-viscous-jacobian-dt.md`)

```markdown
# 熱伝導の近傍 K まで入れたライン粘性 Jacobian (全部入り) を、生産の擬似 dt で成り立たせられるか

## メタ

- **area**: `time_integration`
- **status**: `in_progress`
- **related_docs**:
  - `methods/time_integration/implementation.md` の `lineViscCoupling: 2`・`3` の節 (診断用、本線不採用)
  - `procedures/solver-settings.md` の `lineDtDirectional`・`lineDtDirectionalCap`・`lineViscCoupling` の節
- **related_plans**:
  - [`time_integration-line-viscous-jacobian-faceh.md`](time_integration-line-viscous-jacobian-faceh.md) (前段。§6.2〜§6.13 の切り分けと、この plan へ移した §5.1 #8〜#10)
  - [`../accepted/time_integration-line-viscous-jacobian.md`](../accepted/time_integration-line-viscous-jacobian.md) (親。§4.1 の薄層 Jacobian、§6.0 の値 2 の上限 50、§6.17 の再開の条件)
  - [`time_integration-implicit-thermal-jacobian.md`](time_integration-implicit-thermal-jacobian.md) §4.4 (方向別 dt の熱拡散数の見積)
  - [`time_integration-line-implicit-speed.md`](time_integration-line-implicit-speed.md) (本線: 値 0 + 方向別 + キー 5 + 上限 50。粘性入り L5 = 値 3・マスク 5)
- **created**: `2026-10-11`
- **owner**: Claude (ユーザ 2026-10-10「速度上がる上がらないでなく、筋のいい手法がいい」、2026-10-11 方向別 dt の A/B で B を上限 50 にする案に「ほい」)

## 1. 目的

ユーザの方針は「速度でなく筋のいい手法」。ライン陰解法で壁法線の方向に粘性も熱伝導も結合する「全部入り」の Jacobian (値 3・マスク 7、最終的には値 2) を case/45 で成り立たせたい。
前段 (faceh plan) では、全部入りが**上限なしの方向別 dt** で 300 step 前後で非有限になることの切り分けをした。
面エンタルピーの精度・LHS の座標の差の精度は、単独では回避策にならなかった。熱伝導の近傍 K を外す (マスク 5) と、2000 step 以内の非有限化は避けられた (未収束)。
ただし、**生産で使う擬似 dt (方向別 + 上限 50) で全部入りが回るか**は試していない。この plan はまずそれを確かめ、回れば全部入りと回避策 (マスク 5) の収束を比べる。

## 2. スコープ

- **やる**: 既存のバイナリ (元のセッションの段 ③ の FP64) と設定・診断の切替だけで組む A/B。方向別 dt の上限の 3 水準 (上限なし / 50 / 方向別なし) の比較 (§6)。
  その結果に応じた次の比較 (上限 50 でのマスク 7/5)。faceh plan から移した §5.1 #8〜#10。
- **やらない**: ソルバの数値の変更、値 2・3 の既定化・本番化、マスク 5 の本番化、出発の場・精度・拘束の変更の混入、速度 (総時間) の評価。

## 3. 関連 docs と前提

- 共通関数 `block_dplur::accumulate_thinlayer_visc_jacobian` (`solver_density_cuda/cuda_forge/block_dplur_jacobian_d.cuh` 94〜175 行付近) は、薄層・凍結物性の前処理行列としての近似で、実残差の厳密な微分ではない。
  マスクのビット 2 は熱伝導の近傍 K だけを切り替え、熱伝導の D は常に入る。
- 時間項は D の全行の対角に `V/Δτ` で入る (`timeIntegration_d.cu` 868 行付近)。
- 方向別 dt (`setDT_d.cu` 135〜194 行付近) は、ライン上の節点でライン面の λ を CFL の max から除く。
  case/45 ではライン被覆が 100 % (`[lineImplicit] lines=4719 covered CVs=570999/570999`) なので、全節点が方向別の処理の対象になる。
  実際に Δτ が変わる節点は、残る面がもとから最大かどうかで決まるので実測する (書き出しの 605 節点の結果を全 570999 節点の実測とは扱わない。plan-dt レビュー m7)。
  `lineDtDirectionalCap: R` は Δτ ≤ R × (全面の λ で決めた Δτ) にする。`lineDtDirectional 0` (`cold_cfl.py --line only`) は全節点を point の dt にする。
- `dt_local` を読むのは、平均流の LHS・SST の点陰解法 (`scalarTransport_d.cu` 299 行付近)・化学種 (`speciesTransport_d.cu` 418・485 行付近)・k/ω の更新の対角 (`update_d.cu` 379 行付近)。
  dt の介入の効果は、平均流と SST・化学種で分けられない。
- thermal-jacobian plan §4.4 の見積: 上限なしの方向別 dt では、壁の第一層 (Δn ≈ 1e-6 m) の熱拡散数が約 23。熱拡散数 ≤ 0.5 にするには上限 R ≲ 80 が要り、R = 50 では約 0.2。
- 生産の値 0 は上限 50 (`run_0223`→`run_0224`、60000 step で NaN なし・NOT CONVERGED)。粘性入り L5 (値 3・マスク 5) は上限なしで水準まで到達した (速度 plan §6.14・§6.15)。
- 値 2 (スカラー対角なし) は上限 50 でも `run_0265` が step 99 で ρ の残差 98 倍・k 7000 倍に増えた (親 plan §6.0)。値 3・マスク 7 を上限 50 で回すのは初めて。
- faceh §6.13 の観測: 出発の 10 倍に達した step は、マスク 7 で `rms_roUy` 4 → `rms_roOmega` 14 → `rms_ro` 16 → `rms_roe` 34 → `rms_roUx` 37 → `rms_roK` 64、
  マスク 5 で `rms_roUy` 154 → `rms_ro` 169・174 → `rms_roOmega` 188 → `rms_roe` 217 → `rms_roUx` 220。どちらも `rms_roUy` が先。

## 4. 設計方針

- 1 因子 = 方向別 dt の上限。値 3・マスク 7 (`FORGE_LVC_TERMS=7`)・キー 5・cfl 4・緩和 0.7・sweep 5・ISP 0・出発の場 (`run_0183` の res_100000)・バイナリ (段 ③ の FP64) は全腕で同じにする。
- 水準は 3 つ。A = 上限なし (既知の破綻の再現)、B = 上限 50 (生産の dt)、C = 方向別なし (point の dt)。C は A・B とも 2 本ずつ非有限のときだけ回す。
  C は時間項が支配して何でも有限になりうるので、C の有限を仮説の支持には使わない。C の結果は次の調査の順を決めるだけに使う (§6 の分岐)。
- 介入の成立は、1 step の書き出しで測る。上限のキーは古いバイナリでは黙って無視されるので、`dt_local` が実際に変わったことを確かめる。D の差は時間項の変更だけで説明でき、K は変わらないこと。
- 仮説 (諮問 [記録](../../notes/reviews/2026-10-11-lvcdt-design-diagnose.md)。確度は諮問のもの):
  - **H1 (中)**: 上限なしの方向別 dt では、第一層の熱拡散数が ≫ 1 になり、全部入りの補正が Newton 的になる。
    連続の行はライン面の粘性の減衰を持たないので、K[4][0] を介した第一層の (ρ, T) の組が反復ごとに増える。上限 50 (熱拡散数 ≈ 0.2) なら消えるはず。
    諮問は帳簿で「第一層の ρ が先に動く」と読んだが、主セッションの照合では call 10 の相対変化は節点 4953〜4959 で ρ・T とも約 0.5 % と同程度で、先後は確認できていない。
  - **H2 (中)**: LHS (薄層・凍結物性) と実残差 (勾配から作る熱流束・TP の c_p(T)・node の弱形式) の不整合で、dt を下げても消えない成分がある。
  - **H3 (低)**: 全部入りのライン行列が対角優位を失う (対流の A± と薄層の K の符号の組み合わせ、faceh §6.8 の段 C の「対流の K の相殺」)。
  - H1〜H3 は排他的でない。

## 5. 実装ステップ

1. 台本 `case/45.isobutane_m6_d155/lvcdt.sh`、事前のゲート `lvcdt_pregate.py`、本判定 `lvcdt_judge.py` を、faceh の `lvc75.sh`・`lvc75_pregate.py`・`lvc75_judge.py` を元に作る。
2. codex plan 段のレビューを受け、採否を §6.1 に書いてから回す。
3. 結果を §6.2 以降に記録し、分岐に応じて次の登録をする。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 | 担当 |
| --- | --- | --- | --- |
| 1 | 方向別 dt の上限の A/B (§6。**完了 2026-10-11 §6.2: 回避 (上限なしは 2 本とも DIVERGED、上限 50 は 2 本とも 2000 step 有限・未収束、`rms_roOmega` は末尾で上昇)**。判断: 2026-10-11 diagnostician 諮問 — 回避は登録の文のまま採用、安定・成立とは書かない) | 書き出し 4 本 (A・B・A の再実行・C) → 事前のゲート → A・B 各 2 本 (最大 2000 step)。A・B とも 2 本ずつ非有限のときだけ C を 2 本。合格条件: 事前のゲートが PASS、本判定のゲートが全部通る、§6 の分岐で判定。run は case/45 の `run_0580`〜`run_0589` | F |
| 2 | 上限 50 での全部入り・粘性入り・生産の比較 (#1 が「回避」だった、ユーザ 2026-10-11「よち」、**§6.3 に事前登録**) | 本来の問い「全部入りが、回避策 (マスク 5) と生産 (値 0) より早く水準に届くか」。諮問の提案 (§6.2): 同じバイナリ・上限 50 で値 0・マスク 5・マスク 7 の 3 本を、速度 plan §6.15 の水準 (|欠損| ≤ 0.1 kg/s・θ_r のドリフト ≤ 0.05 %/2 万 step) に届く step 数で比べる (最大 15 万 step、差 ≥ 15000 で判別)。2000 step の残差では比べない。事前のゲートと codex plan 段の後に回す。GPU は 1 本約 1.2 h × FP64 の係数を逐次 | F |
| 3 | §6.8 の未解決の処置 (faceh §5.1 #8 から移した) | 段 C の不一致 (対流の K の相殺)・D の判別不能 3 件・対照 T の保留 (同じソースの監査のビルドと通常のビルドで D が 1 ulp 違う)。§6.9 は係数の変更の成立を示すだけで、これらを解消していない。H3 に進むときに併せて扱う | F |
| 4 | 実残差・更新の写像との整合 (faceh §5.1 #9 から移した、H2) | 実残差の微小応答と有限の補正の応答を、同じ評価の経路で測る (親 plan §6.16 の続き)。測る成分・凍結の範囲・拘束の向きを先に決める | F |

## 6. 検証 (事前登録、2026-10-11、回す前)

- **腕**: すべて段 ③ の FP64 (`~/forge-fgeom3-fp64/solver_density_cuda/build/forge`、sha256 129de3f4…、元のセッションのバイナリを読むだけで使う)。
  - A = `--line dir` (上限なし)。
  - B = `--line dir --cap 50` (`lineDtDirectionalCap: 50.0`)。
  - C (条件付き) = `--line only` (`lineDtDirectional` を書かない = 0)。
  - 共通: `--itj 5 --lvc 3 --cfl 4 --limiter-ref-from run_0183`、`FORGE_LVC_TERMS=7` (7 は表示なしで既定と同じ経路)、出発の場は `run_0183_ns_coldmesh_tw300_ext` の res_100000 (sha256 207d39f0…)。
  - 前提の設定: `axisTimestepBeta`・`lineViscousDtRelief` を書かない (既定 0)。下の用量の関係 Δτ_B = min(Δτ_A, 50 Δτ_C) の前提なので、事前のゲートと本判定で確かめる。
  - 最大 2000 step・100 ごと、序盤 200 step は faceh §6.9 と同じ 112 節点の帳簿。
- **run** (case/45、AWS `~/forge-faceh-audit/case/45.isobutane_m6_d155/`):
  - 書き出し: `run_0584_lvcdt_a_dump`・`run_0585_lvcdt_b_dump`・`run_0586_lvcdt_a_dump2`・`run_0587_lvcdt_c_dump`
  - 腕: `run_0580_lvcdt_a1`・`run_0581_lvcdt_b1`・`run_0582_lvcdt_a2`・`run_0583_lvcdt_b2`
  - 条件付き: `run_0588_lvcdt_c1`・`run_0589_lvcdt_c2`
  - 台本 `lvcdt.sh` (`dumps` / `arms` / `armsC` のモード)、事前のゲート `lvcdt_pregate.py`、本判定 `lvcdt_judge.py`。
- **事前のゲート** (腕の前、PASS のときだけ腕へ)。1 step の書き出し 4 本 (A・B・A の再実行・C、5 本のライン、残差の場は faceh §6.12 と同じ 10 場) から判定する。終了コード 0 = PASS、1 = INVALID、2 = 判別不能。
  - **I 入力**: faceh §6.12 と同じ項目に次を加える。
    - 設定の期待値は腕ごと: A は上限のキーなし、B は `lineDtDirectionalCap: 50.0`、C は `lineDtDirectional` なし。
    - 4 本の設定は、この 2 つのキー以外で全項目同じ。
  - **R 残差の不変**: 全節点の倍精度の残差の場に faceh §6.12 と同じ規則 (A の再実行どうしがビット一致なら B・C もビット一致、そうでなければ max|B − A|・max|C − A| ≤ 3 × max|A の再実行 − A|)。
    定常 (`unsteady 0`) なので残差は `dt_local` を読まないはずだが、測って確かめる。
  - **S 構造の照合**: 4 本それぞれで `rhs_s0` = 拘束の処理をした float(出力 step 1 の残差) が全要素で完全に一致 (faceh §6.12 と同じ)。
  - **V 介入の成立** (B と C のそれぞれについて、A に対して):
    - (a) `Kprev`・`Knext` が A とビット一致 (A の再実行どうしもビット一致のとき。そうでなければ R と同じ規則)。
    - (b) 書き出しの `dt_vol` から節点ごとに r = Δτ_A / Δτ_X を出す。
      - 全 605 節点で r ≥ 1 − 4 ulp (上限・方向別なしとも Δτ は縮む向きしかない)。
      - 3 本の壁法線のライン (節点 1572・4113・7985 の列、縮流部) の第一内部節点 (1571・4112・7984) で r ≥ 10。
      - 全節点で r = 1 なら INVALID (キーが無視された)。
      - C は全節点で r_C ≥ r_B − 4 ulp (point の dt は上限 50 の dt 以下)。
      - **上限 50 の用量** (plan-dt レビュー M2): 同じ状態の C (全面の λ = point の dt) を基準に、全 605 節点で Δτ_B = min(Δτ_A, 50 Δτ_C) が相対 1e-12 で成り立つ
        (`setDT_d.cu`: cfl_B = max(cfl_A, cfl_全面/50)、Δτ = cfl_pseudo / (cfl/dt)。FP64 のビルド)。上限が効いた節点の数と、第一内部節点で効いたかを記録する。
        人工入力で、codex の反例 (B = C = A/50) はこの項目で外れることを確かめた。
    - (c) D: 拘束の無い行の対角は、D_X − D_A を時間項の差 `V (1/Δτ_X − 1/Δτ_A)` (float の `v / dt` と同じ式で再計算) と要素ごとに 64 ulp 以内で再現する
      (尺度は |D_A|・|D_X|・2 つの時間項の大きさの最大)。対角以外の要素と拘束の行は A とビット一致。
    - (c′) 解像: 第一内部節点 1571・4112・7984 のそれぞれで、拘束の無い行のどれかの時間項の差が (c) の許容の 100 倍を超える。
      上限なしの方向別 dt では、第一内部節点の時間項 V/Δτ は D の対角の 1e-4 程度 (`run_0577` の書き出し: 節点 1571 で行 0 の 6.5e-5、行 4 の 4.2e-4) なので、
      行 0 だけでは r = 50 の人工入力でも許容の約 760 倍にとどまる。登録の前の試行で分かったので、行の最大と 100 倍にした。
  - 記録だけ: r の分布 (最小・中央・最大、第一層)、`rhs_s0` の差の ulp と要素の数。
- **分類** (faceh §6.2 と同じ): DIVERGED (detectNaN か残差の非有限) / FINITE (2000 step、全行が有限・非負、100 step ごとの場の保存量・P・T が有限で ρ・P・T が正) / INVALID。
- **分岐** (主判定は「上限 50 で、2000 step 以内の非有限化を避けたか」に限る。H1〜H3 の機構の支持・棄却には使わない。plan-dt レビュー M3):
  1. A が 2 本とも DIVERGED、B が 2 本とも FINITE → **回避**: この条件・2000 step では、上限 50 で全部入りの非有限化を避けた。次は上限 50 でマスク 7 と 5 を比べる (§5.1 #2)。
     全部入りが成り立った・収束した、H1 の機構が正しい、とは言わない (時間項の強化は H2・H3 のどの機構にも効きうる)。
  2. A・B とも 2 本ずつ DIVERGED → **回避せず**: 上限 50 でも避けられない。C を 2 本回す。C の起動は `lvcdt.sh armsC` で、起動の直前に今の証拠で A/B だけを判定し直し
     (`lvcdt_judge.py --ab-only`、ゲート合格・INVALID なし・4 本とも DIVERGED のときだけ終了コード 0)、それが通ったときだけ回す (plan-dt レビュー M5)。
     - C が 2 本とも FINITE → 次の調査の順は H3 から (C の有限は H1 の支持にしない)。行列を見るときは、単位に依る生の対角優位でなく、スケーリングと実際の反復の写像を見る。
     - C が 2 本とも DIVERGED → 次の調査の順は H2 から (§5.1 #4)。
     - C が分かれる → 判別不能。H2 と H3 はどちらも時間項の強化に応答しうるので、C だけでは分けられない。どちらに進むかは「順」であって判別ではない。
  3. A のどちらかが FINITE → 既知の破綻を再現しないので帰属不能。B が分かれる → 判別不能。事前のゲートが外れる → 腕を回さず判別不能 (INVALID なら証拠の不備)。
- **記録だけ (判定に使わない)**:
  - 帳簿の第一内部節点 (1571・4112・7984・4959・6169) の |Δρ/ρ| の最大。窓は 1〜200 呼び出し、破綻した run は 1〜min(200, 破綻の step − 2)。
    対象の (呼び出し, 節点, tag) の欠け・重複・非有限・基準の ρ ≤ 0 は INVALID (plan-dt レビュー M1。人工入力で、欠けた帳簿が INVALID になることを確かめた)。
    前の 3 節点は書き出しの r の節点と同じ列、後の 2 節点は faceh §6.10 で壊れた列 35〜52 をはさむ列 40・50。
    試行 (2026-10-11): faceh §6.13 の 7/5 の run に当てはめると、マスク 7 が 0.544・0.567、マスク 5 が 0.076・0.082。
    初期の変化が累積の擬似時間にほぼ比例するなら、dt を縮めるだけで小さくなるので、この量で「同じ過渡を遅く通っただけか」は決められない (レビュー M3)。
  - 各残差の列が出発の 10 倍に達した step (`first_step_over_10x`)。判定器が `run_0570`・`run_0572` の roUy 4・roOmega 14・ro 16・roe 34 を再現することを確かめた (レビュー m6)。
  - 変数別・列別の非有限の節点の集合は、本判定の後に要約のスクリプト (faceh の `notes/reviews/briefs/2026-10-10-lvc75-result-summary.py` と同じ方法) で出し、§6.2 に書く。
  - **累積の擬似時間は未測定**: 帳簿に `dt_local` は無く、腕の出力は `res_ro`・`volume` だけなので、反復ごとの Δτ の和は作れない。初回の Δτ × step で代用しない (レビュー M4)。
- **見る量 (記録)**: 全残差の列、最初の非有限の step、FINITE の腕の `check_convergence --segment` (区間 0〜1999)。上の「記録だけ」の量。
- **言えること・言えないこと**: 回避でも「この条件・2000 step・上限 50 で非有限化を避けた」まで。全部入りが成り立った (収束した)、H1 の機構が正しい、とは言わない。
  回避せずでも「熱拡散数 ≫ 1 の形の H1 を棄却」とは書かない (ゲートの r ≥ 10 は実際の熱拡散数 ≤ 0.5 を確かめていない。物性・状態も途中で変わり、dt は SST・化学種にも効く)。C の有限を H1 の支持とは書かない。
- **やらないこと**: 上限・cfl・マスクの同時の変更、事後の延長・閾値の変更、マスク 5 の本番化、出発の場・精度・拘束の変更の混入。
- **ディスクと GPU**: 腕は 1 本ずつ。2000 step の腕は 1 本約 3.3 GB。各 run の前に空き 4 GB 未満なら止める。調整板の「GPU の予定」に書いてから出す。判定が済んだら証拠を手元へ移して AWS から消す。

### 6.2 §6 の結果 (2026-10-11、段 ③ の FP64 129de3f4…、diagnostician 諮問 [記録](../../notes/reviews/2026-10-11-lvcdt-result-diagnose.md))

- **事前のゲート (`lvcdt_pregate.py`、`case/45.isobutane_m6_d155/_band_ab/cold_pair/lvcdt_pregate.json`)**: `PASS`。I 93・R 18・S 4・V 16 項目、外れ 0。
  - 上限 50 の第一内部節点の r = Δτ_A/Δτ_B は 1571: 183.7、4112: 170.8、7984: 171.6。上限が効いたのは書き出しの 605 節点のうち 231 (r の中央値 1)。
  - 用量 Δτ_B = min(Δτ_A, 50 Δτ_C) の相対誤差は最大 3.0e-16。
  - 方向別なし (C) の第一層の r は 8539〜9183 (中央値 12.4)。上限なしの方向別 dt は、壁の第一層で point の dt の約 9000 倍だった。
  - K は 3 本ともビット一致した。D の対角は時間項の差と許容の 0.04〜0.05 倍以内で一致した。
- **本判定 (`lvcdt_judge.py`、`_band_ab/cold_pair/lvcdt_judge.json`)**: `主 (値 3・マスク 7): 回避 — この条件・2000 step では、上限 50 で全部入りの非有限化を避けた`。
  ゲート 60 項目すべて合格、INVALID なし。分岐 1 なので C は回していない。

  | run (AWS-A `~/forge-faceh-audit/case/45.isobutane_m6_d155/`) | 方向別 dt | 結果 | `rms_ro`: 0 / 20 / 50 / 100 / 300 / 1000 / 1999 |
  | --- | --- | --- | --- |
  | `run_0580_lvcdt_a1` | 上限なし | step 367 で非有限 (`ro`、列 42〜52 の 1331 節点) | 7.97e-6 / 5.72e-5 / 3.83e-3 / 1.38e-2 / 1.21e-2 / — / — |
  | `run_0581_lvcdt_b1` | 上限 50 | 2000 step 有限・NOT CONVERGED | 7.97e-6 / 5.88e-6 / 5.04e-6 / 4.58e-6 / 4.80e-6 / 4.50e-6 / 4.03e-6 |
  | `run_0582_lvcdt_a2` | 上限なし | step 482 で非有限 (列 42〜52) | step 100 まで A1 と同じ |
  | `run_0583_lvcdt_b2` | 上限 50 | 2000 step 有限・NOT CONVERGED | B1 とほぼ同じ (step 1999 で 4.07e-6) |

- **主張の範囲 (諮問の採否を反映)**:
  - 書くこと: この条件 (段 ③ FP64・`run_0183` の res_100000・値 3・マスク 7・キー 5・cfl 4)・2000 step では、方向別 dt に上限 50 を付けると非有限化を避けた。
    上限なしは 2 本とも 366・481 step で非有限になった。上限 50 では、どの残差の列も出発の 10 倍に達しなかった (最大は `rms_roOmega` の 4.6 倍、`rms_roUy` の 2.4 倍)。
    ただし 2 本とも NOT CONVERGED (`check_convergence --segment` 区間 0〜1999、stalled/plateau。低下桁数 ro 0.4・roUx 0.1・roUy 0.6・roe 0.2・roK 0.0・roOmega 0.8) で、
    **`rms_roOmega` は末尾 500 step の傾きが +0.094・+0.093 桁で上昇している**。2000 step の外で非有限化する可能性は排除していない。
  - 書かないこと: 全部入りが成り立った・安定・収束した。H1 (熱拡散数 ≫ 1) を支持。上限 50 が必要十分 (用量は 1 水準だけ)。平均流の時間項の効果 (SST・化学種を分けていない)。同じ過渡を遅く通った。
- **観測 (記録だけ)**:
  - 帳簿の第一層の |Δρ/ρ| の最大は、上限なし 0.545・0.558、上限 50 は 2 本とも 0.0109。r (171〜184) に比例しない。ただし上限なしは非線形に飽和する過渡なので、解釈には使わない。
  - 上限なしの 10 倍到達の step (roUy 4・roOmega 14・ro 16・roe 34) は faceh §6.13 と同じだった。
  - `rms_roOmega` は全腕で step 1 に出発の 3 倍を超えた。restart で乱流の場が忠実に戻らない既知の過渡の型で、腕の判別には使わない。
- **諮問の採否** (全件採用): 「回避」は登録の文のまま採る。安定・成立とは書かない。Δρ/ρ は記録だけ。restart の共通の過渡は判別に使わない。
  別のバイナリの系列 (値 0・上限 50 の 60000 step、L5 の 11.5 万 step) は参考に限る。熱拡散数・SST と平均流の分離は書かない。
- **H1〜H3**: どれも消えていない (時間項の強化は 3 つの機構すべてに効く)。調査の順は §5.1 #2 の後に決める。
  諮問の提案は H2 (#4) → H3 (#3) → H1 の用量 (上限 20・100・200 の 2000 step)。用量は #2 でマスク 7 が DIVERGED になったときだけ前倒しする。
- **次の一手 (諮問の提案、§5.1 #2、着手はユーザの判断の後)**: 上限 50・同じバイナリで、値 0 (生産)・値 3・マスク 5・値 3・マスク 7 の 3 本を、速度 plan §6.15 の水準に届く step 数で比べる。
  - 水準は、|欠損| ≤ 0.1 kg/s かつ θ_r のドリフト ≤ 0.05 %/2 万 step。
  - 最大 15 万 step・2500 ごと。水準に入った最初の出力で止める。1 本約 1.2 h × FP64 の係数で、逐次に回す。
  - 判定は、到達 step の差が 15000 以上で「速い/遅い」、未満は判別不能。
- **証拠の場所**: AWS-A `~/forge-faceh-audit/case/45.isobutane_m6_d155/` の `run_0580`〜`run_0587`、台本のログ `lvcdt.log`、`lvcdt_pregate.stdout`・`lvcdt_judge.stdout`。result 段のレビューまで残す。

### 6.3 事前登録: 上限 50 での全部入り・粘性入り・生産の比較 (§5.1 #2、2026-10-11、回す前。ユーザ 2026-10-11「よち」「よい」、設計は diagnostician 諮問 [記録](../../notes/reviews/2026-10-11-lvcdt-result-diagnose.md) の Q2)

**問い**: 方向別 dt の上限 50 (生産の dt) で、全部入り (値 3・マスク 7) は、生産 (値 0) と粘性入り (値 3・マスク 5) より早く水準に届くか。
2000 step の窓ではどの腕も横ばいで比べられない (§6.2 の B の低下桁数 0.0〜0.8) ので、残差の桁数・傾きは主指標にしない。

- **腕** (各 1 本、1 本ずつ):
  - L7 = 値 3・`FORGE_LVC_TERMS=7` (全部入り)、L0 = 値 0 (生産)、L5 = 値 3・`FORGE_LVC_TERMS=5` (粘性入り)。
  - 共通: 段 ③ の FP64 (129de3f4…)、`run_0183` の res_100000、`--line dir --cap 50 --itj 5`、cfl 4・緩和 0.7・sweep 5・ISP 0、`axisTimestepBeta`・`lineViscousDtRelief` なし。
  - 値 0 の対照は同じバイナリで取り直す (速度 plan の L0・L5 の系列は別のバイナリなので比べない)。値 2 は足さない。
- **run** (case/45、AWS-A `~/forge-faceh-audit/case/45.isobutane_m6_d155/`):
  - 腕: `run_0553_lvcmp_l7`・`run_0554_lvcmp_l0`・`run_0555_lvcmp_l5` (この順)。
  - 書き出し: `run_0556_lvcmp_l7_dump`・`run_0557_lvcmp_l0_dump`・`run_0558_lvcmp_l5_dump`・`run_0559_lvcmp_l7_dump2`。
  - 番号は 05xx の空き (0590 以降は M6 ノズルのセッションが B で使っている)。台本 `lvcmp.sh`、事前のゲート `lvcmp_pregate.py`、本判定 `lvcmp_judge.py`。
- **長さと終わり**: 最大 150000 step・2500 ごと (水準のドリフトの窓 2 万 step に 8 点以上)。
  見張りは速度 plan の `m9_watch.py --phase line` を使う (水準 = |欠損| ≤ 0.1 kg/s かつ θ_r(40/70/94) のドリフト ≤ 0.05 %/2 万 step、窓は run の中、step 0 は使わない)。
  水準に入った最初の出力で止める (SIGTERM → 確かめてから後の出力を消す。5 万 step ごと・最新 5 つ・到達の出力は残す)。
  - 状態は REACHED / CENSORED (上限まで回って届かず) / DIVERGED / EXEC_ERROR / DATA_ERROR。
  - 起動後 120 s で実効のモードを確かめ、外れれば止める。確かめるのは、LAYOUT2 の分解・代入、値とマスクの表示、面エンタルピーの切替なし。
  - 1 本の目安は、単独で速度 plan の L0・L5 と同じ 1.2 h 前後。GPU を共有するとそれより延びる。
- **事前のゲート** (腕の前、PASS のときだけ腕へ)。1 step の書き出し 4 本 (L7・L0・L5・L7 の再実行):
  - I 入力: faceh §6.12 と同じ項目に次を加える。
    - 4 本の設定は `lineViscCoupling` 以外で同じ (上限 50)。
    - `dt_vol` は 4 本でビット一致し、§6.2 で用量を確かめた上限 50 の書き出し `run_0585_lvcdt_b_dump` ともビット一致する (値・マスクは dt を変えない)。
  - R 残差の不変: L7 の再実行どうしがビット一致なら L5・L0 もビット一致、そうでなければ ≤ 3 × 再実行の差。
  - S 構造の照合 (faceh §6.12 と同じ)。
  - V 介入の成立 (期待の接続の全部、端は 0):
    - (a) L7 − L5: D と K の行 0〜3 が不変。K の行 4 の変化が、独立に計算した熱伝導の近傍 K の項と許容内 (faceh §6.12 と同じ)。
    - (b) L5 − L0: K の行 0 が不変。K の行 1〜3 の変化が、独立に計算した運動量の薄層の項 β P ∂u_j/∂Q_j と許容内。
      ここで β = max(μ_f, 0)|S|²/|e·S|、P = I + n̂n̂ᵀ/3。許容は要素ごとに 1e-5 × 尺度 + 64 ulp。自節点の行の拘束・隣の速度固定では 0。
    - (a′)(b′) それぞれ、変化が解像される面が 100 面以上。
    - 試行 (2026-10-11): 運動量の項の式は、監査の記録と共通関数でマスク 5・4 を組んだ差 (行 1〜3) と、許容の 0.023 倍以内で一致した (1200 面、外れ 0)。
- **判定** (`lvcmp_judge.py`):
  - ゲート: 腕の入力 (バイナリ・出発の場・格子・入力ファイル・腕ごとの設定・ログ)、事前のゲートの `--verify`、腕の設定が同じ腕の書き出しと step 数・出力の間隔・`extraFields` 以外で同じ。
  - 腕ごと: REACHED なら、到達の step を系列から水準の規則で計算し直して見張りと一致すること (速度 plan の `run_0353` の系列で 115000 を再現した)、区間付きの収束 VERDICT が DIVERGED でないこと。
    EXEC_ERROR・DATA_ERROR・記録なし・CENSORED なのに系列が水準に入っている、は INVALID。
  - 比較: 2 本とも REACHED のとき、到達の step の差が 15000 以上なら速い/遅い、未満は判別不能。どちらかが届かなければ比較しない。閾値は事後に下げない。
- **分岐** (主は L7 と L0):
  1. L7 が L0 より 15000 step 以上早く届く → この条件で全部入りの利点を支持。用量 (上限 20・100・200) に進む前に result 段のレビューへ。
  2. 差が 15000 未満 → 「全部入りは生産と区別できない」で閉じる (速度 plan §6.15 の L5 と同じ型。長い run を足さない)。
  3. L7 が遅い、または L7 が届かない (DIVERGED・CENSORED) で L0 は届く → NO-GO。L0 が届かない → 判別不能。
  - L7 と L5 の差 (15000 以上) は「熱伝導の近傍 K が到達を変えるか」の記録に使う。L5 と L0 も記録する。
- **記録だけ**: 到達の窓 [N − 2 万, N] の θ_r・Q_w の `check_quasisteady --series-csv` (`--tail 1.0 --drift 0.0005 --min-snaps 5`)、各腕の区間付きの `check_convergence`。
  NOT CONVERGED は「水準への到達」の比較に限る (速度 plan §6.11 の規則)。
- **やらないこと**: 2000 step の残差で優劣を言う、既存の L0・L5 の系列 (別のバイナリ) と直接比べる、到達の閾値・終わりの条件の事後の変更、point の仕上げの追加。
- **ディスクと GPU**: 腕は 1 本ずつ。見張りが途中の場を消すので 1 本数 GB 以内。各腕の前に空き 6 GB 未満なら止める。調整板の「GPU の予定」に書いてから出す。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| 諮問 (§4・§6 の設計、diagnostician) | 2026-10-10 | [2026-10-11-lvcdt-design-diagnose.md](../../notes/reviews/2026-10-11-lvcdt-design-diagnose.md) | B を上限 50 に、`0` は条件付き C に。Major 4・Minor 2 | ユーザ 2026-10-11「ほい」で採用。介入のゲートに `dt_vol` の実測 (キーの無視を INVALID に)・K のビット一致・D の差 = 時間項を入れた。支持の条件に帳簿の第一層の Δρ/ρ を足した。この plan を新設し faceh §5.1 #8〜#10 を移した。帳簿の「ρ が先」の読みは主セッションの照合で裏付けられず、H1 の文言は確認できた範囲に留めた |
| plan | 2026-10-11 | [2026-10-11-time_integration-line-viscous-jacobian-dt-plan.md](../../notes/reviews/2026-10-11-time_integration-line-viscous-jacobian-dt-plan.md) | GO-with-changes, C0/M5/m2 | 全件採用 (人工入力で確認): M1 帳簿の欠け・重複・非有限・基準 ≤ 0 を INVALID に、窓を事前登録、M2 上限 50 の用量 Δτ_B = min(Δτ_A, 50 Δτ_C) を全 605 節点で相対 1e-12、前提の設定 (`axisTimestepBeta`・`lineViscousDtRelief` なし) を確かめる (反例 B = C = A/50 が外れる)、M3 主判定を「上限 50 で 2000 step 以内の非有限化を避けたか」に限り、Δρ/ρ は記録だけ、H1〜H3 は調査の順として扱う、M4 累積の擬似時間は未測定と明記、M5 C の起動の直前に `--ab-only` で A/B を判定し直す、m6 10 倍到達の step を出力 (既知の値を再現)、m7 「全節点の Δτ が変わる」を訂正 |
| 諮問 (§6 の結果、diagnostician) | 2026-10-11 | [2026-10-11-lvcdt-result-diagnose.md](../../notes/reviews/2026-10-11-lvcdt-result-diagnose.md) | 回避は登録の文のまま採用、Major 3・Minor 3 (安定・成立とは書かない・`rms_roOmega` の末尾の上昇を併記、Δρ/ρ は記録だけ、熱拡散数・SST の分離は書かない)、次は上限 50 での 3 本の水準までの比較 | 全件採用: §6.2 に限定した主張と観測を記録、§5.1 #1 を完了、#2 を 3 本の比較の提案に改めた (ユーザの判断待ち) |

## 7. 影響範囲

- ソルバの数値は変えない。case/45 の台本・判定と run を足す。結論は case/45・この出発の場・2000 step に限る。

## 変更ログ

- 2026-10-11: 起票 (draft)。faceh plan §6.13 の次の一手を、ユーザの判断 (B = 上限 50) で事前登録した。faceh §5.1 #8〜#10 をここへ移した。
- 2026-10-11: codex plan 段 (GO-with-changes, M5/m2) を全件採用し、判定・ゲート・台本を直して in_progress。
- 2026-10-11: §6 の A/B を回し、VERDICT 回避 (§6.2)。diagnostician 諮問を全件採用。次 (§5.1 #2、上限 50 での 3 本の比較) はユーザの判断待ち。
```

## 参考: `case/45.isobutane_m6_d155/lvcmp.sh`

```
#!/bin/bash
# plan time_integration-line-viscous-jacobian-dt §6.3 (2026-10-11): 方向別 dt の上限 50 で、全部入り (値 3・マスク 7)・粘性入り (値 3・マスク 5)・生産 (値 0) の 3 本を、
# 速度 plan §6.15 の水準 (|欠損| ≤ 0.1 kg/s かつ θ_r(40/70/94) のドリフト ≤ 0.05 %/2 万 step) に届く step 数で比べる。
# バイナリは全部元のセッションの段 ③ の FP64 (~/forge-fgeom3-fp64、読むだけ)。run_0183 の res_100000 から、キー 5・cfl 4・方向別・上限 50。
# 順序: ① 1 step の書き出しを L7・L0・L5・L7 の再実行で 4 本 (残差の場つき) → ② 事前のゲート lvcmp_pregate.py → 合格 (終了コード 0) のときだけ
# ③ L7 → L0 → L5 を 1 本ずつ、最大 150000 step・2500 ごと。見張り m9_watch.py (--phase line) が水準に入った最初の出力で止める。
# 引数: なし = ①〜③、dumps = ①② だけ、arms = ③ だけ (lvcmp_pregate.py --verify が通るときだけ)。本判定は lvcmp_judge.py (この台本は判定しない)。
# AWS-A の自分の作業ツリー (~/forge-faceh-audit/case/45.isobutane_m6_d155) で動かす。run_0183 は ~/forge-wallfit の run へのリンク。
set -uo pipefail
cd "$(dirname "$0")"
for v in $(compgen -e | grep '^FORGE_'); do case $v in FORGE_CUDA_BLOCKSIZE) ;; *) unset $v;; esac; done
export FORGE_CUDA_BLOCKSIZE=128
export REAL_CONVERTER=$HOME/forge-wallfit-bin-fp64/solver_density_cuda/build/convertGmshToForge FORGE_CONVERTER=$PWD/conv_tolerant.sh
BIN=$HOME/forge-fgeom3-fp64/solver_density_cuda/build/forge
KEY=fgeom3_fp64
TOOLS=$(cd ../../solver_density_cuda/tools && pwd)
SRC=run_0183_ns_coldmesh_tw300_ext
LOG=$PWD/lvcmp.log
NODES=1572,4113,7985,198560,264263
EXTRA=res_ro,res_roUx,res_roUy,res_roUz,res_roe,res_roK,res_roOmega,res_roY0,res_roY1,volume   # 書き出しの run だけ
BUDGET=150000

echo "== 開始 $(date -Is)" >> $LOG
[ "$(sha256sum $BIN | cut -c1-16)" = 129de3f4e7f67aa3 ] || { echo "バイナリの sha256 が登録と違う — 止める" >> $LOG; exit 1; }
[ -e $SRC ] || ln -s $HOME/forge-wallfit/case/45.isobutane_m6_d155/$SRC $SRC
last=$(ls $SRC | grep -E '^res_[0-9]+\.h5$' | sort -t_ -k2 -n | tail -1)
[ "$last" = res_100000.h5 ] && [ "$(sha256sum $SRC/res_100000.h5 | cut -c1-16)" = 207d39f0e7f4aa03 ] || { echo "出発の場が res_100000.h5 (207d39f0…) でない ($last) — 止める" >> $LOG; exit 1; }

meshsha() {  # meshsha <nozzle.h5> → stdout (MESH の全データセットのハッシュ、場は含めない)
  python3 -c '
import h5py, hashlib, sys
h = hashlib.sha256()
def visit(name, obj):
    if isinstance(obj, h5py.Dataset):
        h.update(name.encode()); h.update(obj[...].tobytes())
with h5py.File(sys.argv[1], "r") as f:
    f["MESH"].visititems(visit)
print(h.hexdigest())' "$1"
}

disk_ok() {  # 共有ディスクの空きが 6 GB 未満なら止める
  local free=$(df --output=avail -B1G ~ | tail -1 | tr -d ' ')
  [ "$free" -ge 6 ] || { echo "ディスクの空き ${free} GB < 6 GB — $1 の前で止める" >> $LOG; return 1; }
}
armargs() {  # armargs <l7|l5|l0>: cold_cfl.py prep に渡す値の引数
  case $1 in l7|l5) echo "--lvc 3";; l0) echo "";; *) return 1;; esac
}
masknum() {  # masknum <l7|l5|l0>: FORGE_LVC_TERMS ("-" = 渡さない)
  case $1 in l7) echo 7;; l5) echo 5;; l0) echo -;; esac
}
modes_ok() {  # modes_ok <run> <l7|l5|l0>: 起動ログの実効のモード
  local f=$1/forge_run.log
  grep -q "^\[line\] factor 1 回目: モード LAYOUT2$" $f && grep -q "^\[line\] solve 1 回目: モード LAYOUT2$" $f || return 1
  ! grep -q "FORGE_DIAG_FACE_H_DOUBLE" $f || return 1
  case $2 in
    l7) grep -q "'lineViscCoupling' in 'time.deltaT': 3" $f && ! grep -q "FORGE_LVC_TERMS=" $f ;;
    l5) grep -q "'lineViscCoupling' in 'time.deltaT': 3" $f && grep -q "診断のマスク FORGE_LVC_TERMS=5 " $f ;;
    l0) ! grep -q "'lineViscCoupling' in 'time.deltaT'" $f && ! grep -q "FORGE_LVC_TERMS=" $f ;;
  esac
}

dump1() {  # dump1 <run> <l7|l5|l0>: 1 step のライン行列の書き出し
  local r=$1 a=$2 mk=$(masknum $2)
  disk_ok $r || return 1
  ( export FORGE_BIN=$BIN COLD_ALT_BINARY=$KEY
    python3 cold_cfl.py prep $SRC $r --steps 1 --out 1 --cfl 4 --limiter-ref-from $SRC --line dir --cap 50 --itj 5 $(armargs $a) --extra $EXTRA >> $LOG 2>&1 ) \
    || { echo "prep $r 失敗 — 止める" >> $LOG; return 1; }
  rm -f $r/nozzle.msh; mkdir -p $r/linedump
  meshsha $r/nozzle.h5 > $r/MESH_SHA.txt || { echo "$r: 格子のハッシュに失敗 — 止める" >> $LOG; return 1; }
  echo $mk > $r/LVC_TERMS.txt
  ( export FORGE_BIN=$BIN COLD_ALT_BINARY=$KEY FORGE_LINE_DUMP_DIR=$PWD/$r/linedump FORGE_LINE_DUMP_CALL=1 FORGE_LINE_DUMP_NODES=$NODES
    [ "$mk" != - ] && export FORGE_LVC_TERMS=$mk
    python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1 ); local rc=$?
  echo "$r rc=$rc $(date -Is) $(grep -c '\[lineDump\] factor の直前を書いた' $r/forge_run.log)" >> $LOG
  [ $rc = 0 ] && grep -q '\[lineDump\] factor の直前を書いた' $r/forge_run.log || { echo "$r: 書き出しが無いか失敗 — 止める" >> $LOG; return 1; }
}

arm() {  # arm <run> <l7|l5|l0>: 水準に入るまで (見張りが止める)
  local r=$1 a=$2 mk=$(masknum $2)
  if [ -f $r/m9_watch.json ] && grep -q '"status": "\(REACHED\|CENSORED\|DIVERGED\|EXEC_ERROR\|DATA_ERROR\)"' $r/m9_watch.json; then echo "$r は既に終わっている" >> $LOG; return 0; fi
  [ -e $r ] && { echo "$r が途中まである — 止める (使い直さない)" >> $LOG; return 1; }
  disk_ok $r || return 1
  ( export FORGE_BIN=$BIN COLD_ALT_BINARY=$KEY
    python3 cold_cfl.py prep $SRC $r --steps $BUDGET --out 2500 --cfl 4 --limiter-ref-from $SRC --line dir --cap 50 --itj 5 $(armargs $a) --extra res_ro >> $LOG 2>&1 ) \
    || { echo "prep $r 失敗 — 止める" >> $LOG; return 1; }
  rm -f $r/nozzle.msh
  meshsha $r/nozzle.h5 > $r/MESH_SHA.txt || { echo "$r: 格子のハッシュに失敗 — 止める" >> $LOG; return 1; }
  echo $mk > $r/LVC_TERMS.txt
  ( export FORGE_BIN=$BIN COLD_ALT_BINARY=$KEY
    [ "$mk" != - ] && export FORGE_LVC_TERMS=$mk
    python3 cold_cfl.py run $r > $r/cold_pair_run_stdout.log 2>&1; echo $? > $r/forge_rc ) &
  local job=$!
  sleep 120
  if ! modes_ok $r $a; then
    echo "$r: 実効のモードが期待と違う — 止める" >> $LOG
    python3 m9_stop.py $r || echo "$r: forge を止められない" >> $LOG
    wait $job; return 1
  fi
  echo "$r: モード OK $(date -Is)" >> $LOG
  python3 m9_watch.py $r --phase line --budget $BUDGET >> $LOG 2>&1; local wrc=$?
  wait $job
  echo "$r 見張り rc=$wrc、forge rc=$(cat $r/forge_rc 2>/dev/null)、状態 $(python3 -c "import json;print(json.load(open('$r/m9_watch.json'))['status'])" 2>/dev/null) $(date -Is)" >> $LOG
  ( cd $r && python3 $TOOLS/check_convergence.py . --segment > CONVERGENCE_SEGMENT.txt 2>&1 ); echo "$r check_convergence --segment rc=$?" >> $LOG
}

MODE=${1:-all}
if [ $MODE = all ] || [ $MODE = dumps ]; then
dump1 run_0556_lvcmp_l7_dump l7 || exit 1
dump1 run_0557_lvcmp_l0_dump l0 || exit 1
dump1 run_0558_lvcmp_l5_dump l5 || exit 1
dump1 run_0559_lvcmp_l7_dump2 l7 || exit 1
python3 lvcmp_pregate.py > lvcmp_pregate.stdout 2>&1; prc=$?
echo "事前のゲート rc=$prc $(date -Is) $(grep '^VERDICT' lvcmp_pregate.stdout)" >> $LOG
[ $prc = 0 ] || { echo "事前のゲートが不合格 (1 = INVALID / 2 = 判別不能) — 腕を回さずに止める" >> $LOG; touch lvcmp.done; exit 1; }
fi
[ $MODE = dumps ] && { echo "== 書き出しと事前のゲートまで $(date -Is)" >> $LOG; touch lvcmp.done; exit 0; }
python3 lvcmp_pregate.py --verify >> $LOG 2>&1 \
  || { echo "事前のゲートの記録が PASS でないか、証拠が記録と違う — 腕を回さずに止める" >> $LOG; touch lvcmp.done; exit 1; }
arm run_0553_lvcmp_l7 l7 || exit 1
arm run_0554_lvcmp_l0 l0 || exit 1
arm run_0555_lvcmp_l5 l5 || exit 1
echo "== 終了 $(date -Is)" >> $LOG
touch lvcmp.done
```

## 参考: `case/45.isobutane_m6_d155/lvcmp_pregate.py`

```
#!/usr/bin/env python3
"""上限 50 での全部入り・粘性入り・生産の比較 (plan time_integration-line-viscous-jacobian-dt §6.3) の事前のゲート。本試験 (run_0553〜0555) の前に回し、合格したときだけ本試験へ進む。

同じ状態 (run_0183 の res_100000) から段 ③ の FP64・方向別 dt 上限 50・キー 5・cfl 4 で 1 step の書き出しを 4 本読む:
  m7  = run_0556_lvcmp_l7_dump  (値 3・FORGE_LVC_TERMS=7、全部入り)
  m7b = run_0559_lvcmp_l7_dump2 (m7 の再実行、atomicAdd による揺れを測る)
  m5  = run_0558_lvcmp_l5_dump  (値 3・FORGE_LVC_TERMS=5、粘性入り = 熱伝導の近傍 K なし)
  m0  = run_0557_lvcmp_l0_dump  (値 0、生産: ライン面の粘性はスカラーの対角だけ)
と、faceh §6.7 の監査の記録 (run_0550_lvcaudit、同じ状態・同じ 605 節点) の面ごとの生の入力を使い、
  I  入力のゲート (faceh §6.12 と同じ項目。設定は 4 本で lineViscCoupling 以外同じ、上限 50。dt_vol は 4 本と dt plan の run_0585 でビット一致)
  R  残差の不変: 出発の状態の残差 (全節点、倍精度) が、m7 の再実行どうしでビット一致なら m5・m0 もビット一致、そうでなければ ≤ 3 × 再実行の差
  S  構造の照合: 各書き出しで rhs_s0 = 拘束の処理をした float(出力 step 1 の残差) が全要素で完全に一致
  V  介入の成立 (期待の接続の全部、端は 0):
     (a) m7 − m5: D と K の行 0〜3 が不変、K の行 4 の変化が独立に計算した熱伝導の近傍 K の項と許容内 (faceh §6.12 と同じ)
     (b) m5 − m0: K の行 0 が不変、K の行 1〜3 の変化が独立に計算した運動量の薄層の項 β P ∂u_j/∂Q_j と許容内 (要素ごとに 1e-5 × |h| + 64 ulp)
     (c') それぞれ、変化が解像される面が 100 面以上
を判定する。出力は _band_ab/cold_pair/lvcmp_pregate.json。終了コード: 0 = 合格、1 = INVALID、2 = 判別不能。使い方: lvcmp_pregate.py [--verify]
"""
import hashlib
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
SRC, SRC_RES, SRC_SHA16 = "run_0183_ns_coldmesh_tw300_ext", "res_100000.h5", "207d39f0e7f4aa03"
SHA_NEW = "129de3f4e7f67aa3a80dbd30d5f5cb75df1582d3974e702d76b61c8998598cec"
DUMPS = {"m7": ("run_0556_lvcmp_l7_dump", SHA_NEW), "m7b": ("run_0559_lvcmp_l7_dump2", SHA_NEW),
         "m5": ("run_0558_lvcmp_l5_dump", SHA_NEW), "m0": ("run_0557_lvcmp_l0_dump", SHA_NEW)}
MASK = {"m7": 7, "m7b": 7, "m5": 5, "m0": None}          # None = 値 0 (FORGE_LVC_TERMS を渡さない、LVC_TERMS.txt は "-")
ARM_CFG = {"m7": {}, "m7b": {}, "m5": {}, "m0": {"time/deltaT/lineViscCoupling": None}}
CAP_REF = "run_0585_lvcdt_b_dump"                         # dt plan §6.2 で用量を確かめた上限 50 の書き出し (同じ状態)
AUDIT = "run_0550_lvcaudit"
WANT_NODES = [1572, 4113, 7985, 198560, 264263]
RES_FIELDS = ["res_ro", "res_roUx", "res_roUy", "res_roUz", "res_roe", "res_roK", "res_roOmega"]   # res_roUz は rhs_s0 の行 3 の照合に使う (§6.12)
RES_OPT = ["res_roY0", "res_roY1"]          # 化学種の残差 (3 本ともにあれば R に含める)
EXTRA = RES_FIELDS + RES_OPT + ["volume"]
SAME_FILES = ["bcondConfig.yaml", "probe.yaml", "species_meta.yaml", "wall_design.csv", "wall_physical.csv", "target_axis_M.csv", "wall_repr.json", "MESH_QUALITY.txt"]
EXPECT = {"time/deltaT/cfl": 4.0, "time/deltaT/cfl_pseudo": 4.0, "time/deltaT/implicitRelax": 0.7, "time/deltaT/implicitThermalJacobian": 5,
          "time/deltaT/lineDtDirectional": 1, "time/deltaT/lineImplicit": 1, "time/deltaT/lineViscCoupling": 3, "time/deltaT/blockDPLUR": 1,
          "time/deltaT/detectNaN": 1, "time/timeIntegration": 11, "time/nStepInner": 5, "time/last/nStepOuter": 1,
          "time/deltaT/lineDtDirectionalCap": 50.0, "time/deltaT/implicitSolvePrecision": None, "time/outStepInterval": 1,
          "output/extraFields": EXTRA, "time/deltaT/axisTimestepBeta": None, "time/deltaT/lineViscousDtRelief": None}
NF, FREC, NREC, HF = 12, 200, 160, 150
f32 = np.float32


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def meshsha(p):
    h = hashlib.sha256()

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            h.update(name.encode()); h.update(obj[...].tobytes())
    with h5py.File(p, "r") as f:
        f["MESH"].visititems(visit)
    return h.hexdigest()


def flat(d, p=""):
    out = {}
    for k, v in d.items():
        kk = f"{p}/{k}" if p else str(k)
        out.update(flat(v, kk) if isinstance(v, dict) else {kk: v})
    return out


def dump(run, name):
    d = HERE / run / "linedump"
    meta = {l.split()[0]: (int(l.split()[1]), int(l.split()[2])) for l in (d / "meta.txt").read_text().splitlines() if l.strip() and not l.startswith("#")}
    r, c = meta[name]
    a = np.fromfile(d / f"{name}.f64", dtype=np.float64)
    if a.size != r * c:
        raise ValueError(f"{run}/{name}: 大きさが違う")
    return a.reshape(r, c)


EVIDENCE_FILES = ["linedump/meta.txt", "linedump/D.f64", "linedump/Kprev.f64", "linedump/Knext.f64", "linedump/state_ro_roU_roe_cp_gamma.f64",
                  "linedump/dt_vol.f64", "linedump/flags_wall_iso_axis.f64", "linedump/node_line.f64", "linedump/rhs_s0.f64",
                  "res_0.h5", "res_1.h5", "solverConfig.yaml", "LVC_TERMS.txt", "forge_run.log", "RUN_PROVENANCE.txt"]
AUDIT_FILES = ["linedump/audit_face.f64", "linedump/audit_node.f64", "linedump/meta.txt", "linedump/node_line.f64",
               "linedump/state_ro_roU_roe_cp_gamma.f64", "linedump/dt_vol.f64", "linedump/flags_wall_iso_axis.f64"]


def evidence():
    """事前のゲートが読んだ証拠と、このスクリプト自身・書き出しの設定のハッシュ (本試験の起動と本判定で照合する)"""
    ev = {f"{run}/{fn}": sha(HERE / run / fn) for run, _ in DUMPS.values() for fn in EVIDENCE_FILES}
    ev.update({f"{AUDIT}/{fn}": sha(HERE / AUDIT / fn) for fn in AUDIT_FILES})
    ev[f"{CAP_REF}/linedump/dt_vol.f64"] = sha(HERE / CAP_REF / "linedump" / "dt_vol.f64")
    ev["lvcmp_pregate.py"] = sha(Path(__file__).resolve())
    return ev


def verify():
    """--verify: 記録が PASS で、記録したハッシュが今のファイルと同じなら 0、そうでなければ 1"""
    rec = json.loads((HERE / "_band_ab" / "cold_pair" / "lvcmp_pregate.json").read_text())
    if not str(rec.get("VERDICT", "")).startswith("PASS"):
        print("事前のゲートの記録が PASS でない:", rec.get("VERDICT")); return 1
    ev = rec.get("evidence_sha256") or {}
    if not ev:
        print("証拠のハッシュが記録に無い"); return 1
    now = evidence()
    bad = sorted(k for k in set(ev) | set(now) if ev.get(k) != now.get(k))
    if bad:
        print("記録と違う証拠:", bad[:10]); return 1
    print(f"照合 OK: {len(ev)} ファイル"); return 0


def ulp32(x):
    return np.spacing(np.abs(np.asarray(x, dtype=np.float32))).astype(np.float64)


def main():
    out = HERE / "_band_ab" / "cold_pair" / "lvcmp_pregate.json"
    rec = {"plan": "time_integration-line-viscous-jacobian-dt §6.3 (事前のゲート)"}
    checks = []

    def ok(stage, name, cond, val=None):
        checks.append((stage, name, bool(cond), val))
    verdict = None
    try:
        # ---- I: 入力 ----
        src_full = sha(HERE / SRC / SRC_RES)
        ok("I", "出発の場の sha256", src_full.startswith(SRC_SHA16), src_full[:16])
        mref = meshsha(HERE / SRC / "nozzle.h5")
        cfgs = {}
        for key, (run, bsha) in DUMPS.items():
            d = HERE / run
            prov = (d / "RUN_PROVENANCE.txt").read_text(errors="replace")
            ok("I", f"{run}: バイナリの sha256", bsha in prov, bsha[:16])
            cp = json.loads((d / "COLD_PAIR.json").read_text())
            ok("I", f"{run}: 出発の場", cp["parent"] == SRC and cp["parent_res"] == SRC_RES and cp["parent_res_sha256"] == src_full, cp["parent_res_sha256"][:16])
            m = meshsha(d / "nozzle.h5")
            ok("I", f"{run}: 格子の実体", m == mref and (d / "MESH_SHA.txt").read_text().strip() == m, m[:16])
            ok("I", f"{run}: 入力ファイル", all((d / fn).is_file() and sha(d / fn) == sha(HERE / SRC / fn) for fn in SAME_FILES), None)
            c = flat(yaml.safe_load((d / "solverConfig.yaml").read_text()))
            cfgs[key] = c
            exp_ = {**EXPECT, **ARM_CFG[key]}
            bad = {k: c.get(k) for k, v in exp_.items() if (k in c if v is None else c.get(k) != v)}
            ok("I", f"{run}: 設定 (腕ごとの期待どおり)", bad == {}, bad)
            log = (d / "forge_run.log").read_text(errors="replace")
            mk = MASK[key]
            shown = (f"診断のマスク FORGE_LVC_TERMS={mk} " in log) if mk == 5 else ("FORGE_LVC_TERMS=" not in log)
            vline = ("'lineViscCoupling' in 'time.deltaT': 3" in log) if mk is not None else ("'lineViscCoupling' in 'time.deltaT': 3" not in log)
            ok("I", f"{run}: 起動ログ (値・マスクの表示・LAYOUT2・面エンタルピーの切替なし・書き出しあり)",
               vline and shown and "Thomas の配列の並び: LAYOUT2" in log and "FORGE_DIAG_FACE_H_DOUBLE" not in log and "[lineDump] factor の直前を書いた" in log, None)
            lt = (d / "LVC_TERMS.txt").read_text().strip() if (d / "LVC_TERMS.txt").is_file() else None
            ok("I", f"{run}: 台本が渡した FORGE_LVC_TERMS = {mk if mk is not None else '-'}", lt == (str(mk) if mk is not None else "-"), lt)
            ok("I", f"{run}: RUN_RC = 0", (d / "RUN_RC").read_text().strip() == "0", None)
            ok("I", f"{run}: extraFields の無視の警告なし", "は確保されていない変数なので無視する" not in log, None)
        strip = lambda c: {k: v for k, v in c.items() if k != "time/deltaT/lineViscCoupling"}
        ok("I", "4 本の設定が lineViscCoupling 以外で完全に同じ", all(strip(cfgs[k]) == strip(cfgs["m7"]) for k in DUMPS), None)
        rec["dump_config"] = cfgs
        arrays = {}
        for key, (run, _) in DUMPS.items():
            d = HERE / run
            a = {n: dump(run, n) for n in ("node_line", "D", "Kprev", "Knext", "state_ro_roU_roe_cp_gamma", "dt_vol", "flags_wall_iso_axis", "rhs_s0")}
            nn = a["node_line"].shape[0]
            shapes = {"node_line": 2, "D": 25, "Kprev": 25, "Knext": 25, "state_ro_roU_roe_cp_gamma": 7, "dt_vol": 2, "flags_wall_iso_axis": 3, "rhs_s0": 5}
            ok("I", f"{run}: 配列の形", nn > 0 and all(a[n].shape == (nn, w) for n, w in shapes.items()), {n: a[n].shape for n in shapes})
            nodes_ = [int(x) for x in a["node_line"][:, 0]]
            ok("I", f"{run}: 要求した 5 節点・5 本のライン・重複なし", all(w in nodes_ for w in WANT_NODES) and len(set(a["node_line"][:, 1])) == 5 and len(set(nodes_)) == nn, None)
            ok("I", f"{run}: 有限", all(bool(np.all(np.isfinite(a[n]))) for n in shapes), None)
            st_ = a["state_ro_roU_roe_cp_gamma"]
            ok("I", f"{run}: ρ・c_p・dt・体積 > 0、γ > 1", bool(np.all(st_[:, 0] > 0) and np.all(st_[:, 5] > 0) and np.all(st_[:, 6] > 1) and np.all(a["dt_vol"] > 0)), None)
            with h5py.File(d / "res_1.h5", "r") as h:
                have = set(h["VALUE"].keys())
            ok("I", f"{run}: res_1.h5 に残差の場がそろう", all(f in have for f in RES_FIELDS), sorted(set(RES_FIELDS) - have))
            arrays[key] = a
            arrays[key]["_have"] = have
        for n in ("node_line", "state_ro_roU_roe_cp_gamma", "flags_wall_iso_axis", "dt_vol"):
            ok("I", f"{n}: 4 本でビット一致 (同じ入力・同じ上限 50 の dt)",
               all(arrays[k][n].shape == arrays["m7"][n].shape and np.array_equal(arrays[k][n].view(np.int64), arrays["m7"][n].view(np.int64)) for k in DUMPS), None)
        capref = dump(CAP_REF, "dt_vol")
        ok("I", f"dt_vol が {CAP_REF} (dt plan §6.2 で用量を確かめた上限 50) とビット一致", capref.shape == arrays["m7"]["dt_vol"].shape and np.array_equal(capref.view(np.int64), arrays["m7"]["dt_vol"].view(np.int64)), None)
        def h5vals(run):
            with h5py.File(HERE / run / "res_0.h5", "r") as h:
                return {k: h["VALUE"][k][...] for k in h["VALUE"].keys()}
        v0 = {k: h5vals(DUMPS[k][0]) for k in DUMPS}
        keys0 = sorted(v0["m7"].keys())
        ok("I", "res_0.h5 の全データセットが 4 本でビット一致",
           len(keys0) > 0 and all(sorted(v0[k].keys()) == keys0 and all(np.array_equal(v0[k][n], v0["m7"][n]) for n in keys0) for k in DUMPS), keys0)
        del v0
        aud = {n: dump(AUDIT, n) for n in ("node_line", "state_ro_roU_roe_cp_gamma", "flags_wall_iso_axis")}
        for key in DUMPS:
            for n in ("node_line", "state_ro_roU_roe_cp_gamma", "flags_wall_iso_axis"):
                ok("I", f"{DUMPS[key][0]}: {n} が監査の記録とビット一致", arrays[key][n].shape == aud[n].shape and np.array_equal(arrays[key][n].view(np.int64), aud[n].view(np.int64)), None)
        Nchk = np.fromfile(HERE / AUDIT / "linedump" / "audit_node.f64").reshape(aud["node_line"].shape[0], NREC)
        ok("I", "監査の記録: 記録の節点番号 = 書き出しの順、面の数が枠に入る、skipDiag なし",
           np.array_equal(Nchk[:, 0], aud["node_line"][:, 0]) and bool(np.all(Nchk[:, 5] == 0)) and bool(np.all(Nchk[:, 4] <= NF)) and bool(np.all(Nchk[:, 22] == 0)), None)
        if not all(c[2] for c in checks):
            raise ValueError("入力のゲートが不合格")

        # ---- R: 残差の不変 ----
        def res_field(run, f):
            with h5py.File(HERE / run / "res_1.h5", "r") as h:
                return h["VALUE"][f][...].astype(np.float64)
        resv = {}
        opt = [f for f in RES_OPT if all(f in arrays[k]["_have"] for k in DUMPS)]
        for f in RES_FIELDS + opt:
            fa, fa2 = res_field(DUMPS["m7"][0], f), res_field(DUMPS["m7b"][0], f)
            fx = {k: res_field(DUMPS[k][0], f) for k in ("m5", "m0")}
            if not (np.all(np.isfinite(fa)) and np.all(np.isfinite(fa2)) and all(np.all(np.isfinite(v)) for v in fx.values())):
                raise ValueError(f"{f} に非有限がある (出発の状態の残差)")
            rbit = bool(np.array_equal(fa.view(np.int64), fa2.view(np.int64)))
            d_oo = float(np.max(np.abs(fa2 - fa)))
            resv[f] = {"rerun_bit": rbit, "rerun": d_oo}
            for k, v in fx.items():
                d_no = float(np.max(np.abs(v - fa)))
                ok("R", f"{f}: {k} と m7 の差 ({'ビット一致を要求' if rbit else '≤ 3 × 再実行の差'})",
                   bool(np.array_equal(v.view(np.int64), fa.view(np.int64))) if rbit else d_no <= 3 * d_oo, [d_no, d_oo])
                resv[f][f"{k}_vs_m7"] = d_no
        rec["residual"] = resv
        # ---- S: 構造の照合 ----
        rows = [(0, "res_ro"), (1, "res_roUx"), (2, "res_roUy"), (3, "res_roUz"), (4, "res_roe")]
        struct = {}
        for key, (run, _) in DUMPS.items():
            a = arrays[key]
            idx = [int(x) for x in a["node_line"][:, 0]]
            fl_ = a["flags_wall_iso_axis"]
            with h5py.File(HERE / run / "res_1.h5", "r") as h:
                R = {f: h["VALUE"][f][...] for _, f in rows}
            mism = {}
            for row, f in rows:
                raw = R[f][idx]
                if not np.all(np.isfinite(raw)):
                    raise ValueError(f"{run}: {f} の書き出しの節点に非有限")
                exp = raw.astype(np.float32).astype(np.float64)
                if row in (1, 2, 3):
                    exp = np.where(fl_[:, 0] == 1, 0.0, exp)
                if row == 2:
                    exp = np.where(fl_[:, 2] == 1, 0.0, exp)
                if row == 4:
                    exp = np.where(fl_[:, 1] == 1, 0.0, exp)
                mism[row] = int(np.sum(a["rhs_s0"][:, row] != exp))
            struct[run] = mism
            ok("S", f"{run}: rhs_s0 = 拘束の処理をした float(残差) が全要素で完全に一致", sum(mism.values()) == 0, mism)
        rec["struct_rhs"] = struct

        # ---- V: 介入の成立 ----
        ad = HERE / AUDIT / "linedump"
        nn = aud["node_line"].shape[0]
        F = np.fromfile(ad / "audit_face.f64").reshape(nn, NF, FREC)
        N = np.fromfile(ad / "audit_node.f64").reshape(nn, NREC)
        st = arrays["m7"]["state_ro_roU_roe_cp_gamma"]
        nodes = [int(x) for x in arrays["m7"]["node_line"][:, 0]]
        pos = {n: k for k, n in enumerate(nodes)}
        flags = arrays["m7"]["flags_wall_iso_axis"]
        nl = arrays["m7"]["node_line"]

        def same_or_noise(name, x, ref, sl):
            a, a2, b = arrays[ref][name][:, sl], arrays["m7b"][name][:, sl], arrays[x][name][:, sl]
            a7, a7b = arrays["m7"][name][:, sl], arrays["m7b"][name][:, sl]
            rbit = bool(np.array_equal(a7.view(np.int64), a7b.view(np.int64)))
            d_ba, d_rr = float(np.max(np.abs(b - a))), float(np.max(np.abs(a7b - a7)))
            return (bool(np.array_equal(b.view(np.int64), a.view(np.int64))) if rbit else d_ba <= 3 * d_rr), {"rerun_bit": rbit, "diff": d_ba, "rerun": d_rr}
        same = {}
        okD, same["D_m5_vs_m7"] = same_or_noise("D", "m5", "m7", slice(0, 25))
        ok("V", "(a) D が m7 と m5 で不変 (熱伝導の D は常に入る)", okD, same["D_m5_vs_m7"])
        for nm in ("Kprev", "Knext"):
            okK, same[f"{nm}_rows0to3_m5_vs_m7"] = same_or_noise(nm, "m5", "m7", slice(0, 20))
            ok("V", f"(a) {nm} の行 0〜3 が m7 と m5 で不変", okK, same[f"{nm}_rows0to3_m5_vs_m7"])
            okK0, same[f"{nm}_row0_m0_vs_m5"] = same_or_noise(nm, "m0", "m5", slice(0, 5))
            ok("V", f"(b) {nm} の行 0 が m5 と m0 で不変 (連続の行に粘性の流束は無い)", okK0, same[f"{nm}_row0_m0_vs_m5"])
        rec["same"] = same
        exp_conn = []
        for k in range(nn):
            same_prev = k > 0 and nl[k - 1, 1] == nl[k, 1]
            same_next = k + 1 < nn and nl[k + 1, 1] == nl[k, 1]
            exp_conn.append((k, 0, nodes[k - 1] if same_prev else -1))
            exp_conn.append((k, 1, nodes[k + 1] if same_next else -1))
        conn_bad = []
        res = {"heat": [0.0, 0, 0], "mom": [0.0, 0, 0]}     # [比の最大, 比べた面, 解像された面]
        for k, side, nb in exp_conn:
            nd = N[k]
            if int(nd[0]) != nodes[k] or int(nd[1 + side]) != nb:
                conn_bad.append((nodes[k], side, "記録の line_prev/next が並びと違う")); continue
            wall, iso, ax = int(flags[k, 0]), int(flags[k, 1]), int(flags[k, 2])
            exp_dec = (14 if wall == 1 else 0) | (16 if iso == 1 else 0) | (4 if ax == 1 else 0)
            if int(nd[18]) != exp_dec:
                conn_bad.append((nodes[k], side, "拘束の行のビットがフラグと違う")); continue
            nm = "Kprev" if side == 0 else "Knext"
            K7, K5, K0 = (arrays[x][nm][k].reshape(5, 5) for x in ("m7", "m5", "m0"))
            if nb < 0:
                if not (np.all(K7 == 0) and np.all(K5 == 0) and np.all(K0 == 0)):
                    conn_bad.append((nodes[k], side, "ラインの端なのに K が 0 でない"))
                continue
            fs = [s_ for s_ in range(int(min(nd[4], NF))) if F[k, s_, 5] != 0 and int(F[k, s_, 3]) == nb]
            if len(fs) != 1:
                conn_bad.append((nodes[k], side, "期待の接続にライン面がちょうど 1 つでない")); continue
            fr = F[k, fs[0]]
            kb = pos.get(nb)
            back = [s_ for s_ in range(int(min(N[kb, 4], NF))) if F[kb, s_, 5] != 0 and int(F[kb, s_, 3]) == nodes[k]] if kb is not None else []
            if not (int(fr[6]) == side and int(fr[45]) == 1 and {int(fr[1]), int(fr[2])} == {nodes[k], nb} and kb is not None
                    and len(back) == 1 and int(F[kb, back[0], 0]) == int(fr[0]) and int(F[kb, back[0], 6]) == 1 - side):
                conn_bad.append((nodes[k], side, "面の向き・枝・両端・往復の面が合わない")); continue
            if bool(fr[43] != 0) != (int(flags[kb, 0]) == 1) or bool(fr[44] != 0) != (int(flags[kb, 1]) == 1):
                conn_bad.append((nodes[k], side, "隣の速度・温度固定のフラグが書き出しのフラグと違う")); continue
            f0 = fr[7]
            ic0_is_k = int(fr[1]) == nodes[k]
            e = np.array([fr[15] - fr[12], fr[16] - fr[13], fr[17] - fr[14]])
            dod = fr[11] ** 2 / abs(float(e @ fr[8:11]))
            rj, uj, vj, wj, rEj, gj = fr[52], fr[53], fr[54], fr[55], fr[56], fr[51]
            # 熱伝導の近傍 K (m7 − m5、行 4)
            kf = f0 * fr[186] + (1 - f0) * fr[187] + (f0 * fr[184] + (1 - f0) * fr[185]) * (f0 * fr[182] + (1 - f0) * fr[183]) / fr[50]
            kap = max(kf, 0.0) * dod
            cpj = max(st[kb, 5], 1e-30)
            q2 = uj * uj + vj * vj + wj * wj
            ej = rEj / rj - 0.5 * q2
            hh = kap * (gj / cpj) / rj * np.array([-(ej - 0.5 * q2), -uj, -vj, -wj, 1.0])
            if iso == 1 or int(flags[kb, 1]) == 1:
                hh = np.zeros(5)
            # 運動量の薄層 (m5 − m0、行 1〜3): K[1+a][0] = −β Σ_b P_ab u_j,b / ρ_j、K[1+a][1+b] = β P_ab / ρ_j
            mlam0, mlam1 = (fr[188], fr[189]) if ic0_is_k else (fr[189], fr[188])
            muf = f0 * (mlam0 + fr[182]) + (1 - f0) * (mlam1 + fr[183])
            beta = max(muf, 0.0) * dod
            nvec = np.array(fr[8:11]) / fr[11]
            P = np.eye(3) + np.outer(nvec, nvec) / 3.0
            ujv = np.array([uj, vj, wj])
            hm = np.zeros((3, 5)); hs = np.zeros((3, 5))
            for a_ in range(3):
                hm[a_, 0] = -beta * float(P[a_] @ ujv) / rj
                hs[a_, 0] = beta * float(np.abs(P[a_]) @ np.abs(ujv)) / rj
                hm[a_, 1:4] = beta * P[a_] / rj
                hs[a_, 1:4] = np.abs(hm[a_, 1:4])
            for a_ in range(3):
                if (wall == 1) or (ax == 1 and a_ == 1) or int(flags[kb, 0]) == 1:
                    hm[a_] = 0.0; hs[a_] = 0.0
            for tag, dK, h, hsc in (("heat", K7[4] - K5[4], hh, np.abs(hh)), ("mom", (K5[1:4] - K0[1:4]).ravel(), hm.ravel(), hs.ravel())):
                scale = np.maximum.reduce([np.abs(dK), hsc, np.abs(h)])
                if tag == "heat":
                    scale = np.maximum.reduce([scale, np.abs(K7[4]), np.abs(K5[4])])
                else:
                    scale = np.maximum.reduce([scale, np.abs(K5[1:4]).ravel(), np.abs(K0[1:4]).ravel()])
                floor = 64 * np.maximum(ulp32(scale), np.finfo(np.float32).tiny)
                tol = 1e-5 * hsc + floor
                ratio = np.abs(dK - h) / tol
                if not (np.isfinite(kap) and np.isfinite(beta) and np.all(np.isfinite(h)) and np.all(np.isfinite(ratio))):
                    raise ValueError(f"節点 {nodes[k]} 向き {side}: {tag} の比に非有限 (入力の不備)")
                res[tag][0] = max(res[tag][0], float(np.max(ratio))); res[tag][1] += 1
                if float(np.max(np.abs(h))) > 1e3 * float(np.max(floor)):
                    res[tag][2] += 1
        rec["intervention"] = {"heat_m7_minus_m5": res["heat"], "mom_m5_minus_m0": res["mom"], "connection_problems": conn_bad[:10], "expected_connections": len(exp_conn)}
        if conn_bad:
            raise ValueError(f"期待の接続と監査の記録が合わない ({len(conn_bad)} 件、最初 {conn_bad[0]})")
        ok("V", "(a) K の行 4 の変化 (m7 − m5) が独立に計算した熱伝導の項と許容内 (要素ごとに 1e-5 × |h| + 64 ulp)、端は 0", res["heat"][1] > 0 and res["heat"][0] <= 1.0, res["heat"])
        ok("V", "(a') 熱伝導の項が解像される面が 100 面以上", res["heat"][2] >= 100, res["heat"][2])
        ok("V", "(b) K の行 1〜3 の変化 (m5 − m0) が独立に計算した運動量の薄層の項と許容内 (要素ごとに 1e-5 × 尺度 + 64 ulp)、端は 0", res["mom"][1] > 0 and res["mom"][0] <= 1.0, res["mom"])
        ok("V", "(b') 運動量の項が解像される面が 100 面以上", res["mom"][2] >= 100, res["mom"][2])
        i_ok = all(c[2] for c in checks if c[0] == "I")
        r_ok = all(c[2] for c in checks if c[0] in ("R", "S"))
        v_ok = all(c[2] for c in checks if c[0] == "V")
        if not i_ok:
            raise ValueError("入力のゲートに外れがある")
        verdict = "PASS (本試験へ進む)" if (r_ok and v_ok) else "INDETERMINATE (残差の変化・構造の照合の外れ・介入の不成立のどれか: 本試験へ進まない)"
        if verdict.startswith("PASS"):
            rec["evidence_sha256"] = evidence()       # 本試験・本判定がこの記録と照合する (M3)
    except Exception as e:
        verdict = f"INVALID ({type(e).__name__}: {e})"
    rec["checks"] = [{"stage": a, "check": b, "ok": c, "value": d} for a, b, c, d in checks]
    rec["VERDICT"] = verdict
    out.write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str))
    for st in ("I", "R", "S", "V"):
        lst = [c for c in checks if c[0] == st]
        bad = [c for c in lst if not c[2]]
        print(f"{st}: {len(lst)} 項目、外れ {len(bad)}  {[(c[1], c[3]) for c in bad[:3]]}")
    for k in ("residual", "struct_rhs", "same", "intervention"):
        if k in rec:
            print(k, rec[k])
    print("VERDICT:", verdict)
    return 0 if verdict.startswith("PASS") else (1 if verdict.startswith("INVALID") else 2)


if __name__ == "__main__":
    sys.exit(verify() if sys.argv[1:] == ["--verify"] else main())
```

## 参考: `case/45.isobutane_m6_d155/lvcmp_judge.py`

```
#!/usr/bin/env python3
"""上限 50 での全部入り・粘性入り・生産の比較の判定 (plan time_integration-line-viscous-jacobian-dt §6.3、2026-10-11 事前登録)。

lvcmp.sh の run_0553 (L7 = 値 3・マスク 7)・run_0554 (L0 = 値 0、生産)・run_0555 (L5 = 値 3・マスク 5) と、事前のゲートの記録
(lvcmp_pregate.json、run_0556〜0559 の 1 step の書き出しから) を読み、
  (1) 入力のゲート: バイナリ・出発の場・格子の実体・入力ファイル・腕ごとの期待どおりの設定・起動ログ (値・マスク・LAYOUT2)、
      事前のゲートの --verify、腕の設定が同じ腕の書き出しの設定と step 数・出力の間隔・extraFields 以外で同じ
  (2) 腕ごと: 見張り (m9_watch.json) の状態。REACHED なら、到達の step を系列から水準の規則で計算し直して一致すること、
      区間付きの収束 VERDICT があって DIVERGED でないこと。EXEC_ERROR・DATA_ERROR・記録なしは INVALID
  (3) 比較: 2 本とも REACHED のとき、到達の step の差が 15000 以上なら速い/遅い、未満は判別不能。どちらかが届かなければ比較しない
  (4) 記録だけ: 到達の窓 [N − 2 万, N] の θ_r・Q_w の check_quasisteady (--tail 1.0 --drift 0.0005 --min-snaps 5)、各列の区間付きの収束
を行い、§6.3 の分岐どおりに判定する。出力は _band_ab/cold_pair/lvcmp_judge.json (異常でも必ず書く)。ゲートの不合格か INVALID があれば判定せず、終了コード 1。
"""
import csv
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
D = HERE / "_band_ab" / "cold_pair"
SHA_NEW = "129de3f4e7f67aa3a80dbd30d5f5cb75df1582d3974e702d76b61c8998598cec"
SRC, SRC_RES, SRC_SHA16 = "run_0183_ns_coldmesh_tw300_ext", "res_100000.h5", "207d39f0e7f4aa03"
ARMS = {"L7": ("run_0553_lvcmp_l7", "7", "m7"), "L0": ("run_0554_lvcmp_l0", "-", "m0"), "L5": ("run_0555_lvcmp_l5", "5", "m5")}
DUMPS = ["run_0556_lvcmp_l7_dump", "run_0557_lvcmp_l0_dump", "run_0558_lvcmp_l5_dump", "run_0559_lvcmp_l7_dump2"]
SAME_FILES = ["bcondConfig.yaml", "probe.yaml", "species_meta.yaml", "wall_design.csv", "wall_physical.csv", "target_axis_M.csv", "wall_repr.json", "MESH_QUALITY.txt"]
BUDGET, OUT, THR = 150000, 2500, 15000
EXPECT = {"time/deltaT/cfl": 4.0, "time/deltaT/cfl_pseudo": 4.0, "time/deltaT/implicitRelax": 0.7, "time/deltaT/implicitThermalJacobian": 5,
          "time/deltaT/lineDtDirectional": 1, "time/deltaT/lineImplicit": 1, "time/deltaT/lineDtDirectionalCap": 50.0, "time/deltaT/blockDPLUR": 1,
          "time/deltaT/detectNaN": 1, "time/timeIntegration": 11, "time/nStepInner": 5, "time/last/nStepOuter": BUDGET, "time/outStepInterval": OUT,
          "time/deltaT/implicitSolvePrecision": None, "time/deltaT/axisTimestepBeta": None, "time/deltaT/lineViscousDtRelief": None,
          "output/extraFields": ["res_ro"]}
ARM_CFG = {"L7": {"time/deltaT/lineViscCoupling": 3}, "L5": {"time/deltaT/lineViscCoupling": 3}, "L0": {"time/deltaT/lineViscCoupling": None}}
KEYS = ("theta_r_40", "theta_r_70", "theta_r_94", "Q_w")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def meshsha(p):
    h = hashlib.sha256()

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            h.update(name.encode()); h.update(obj[...].tobytes())
    with h5py.File(p, "r") as f:
        f["MESH"].visititems(visit)
    return h.hexdigest()


def flat(d, p=""):
    out = {}
    for k, v in d.items():
        kk = f"{p}/{k}" if p else str(k)
        out.update(flat(v, kk) if isinstance(v, dict) else {kk: v})
    return out


def cfg(run):
    return flat(yaml.safe_load((HERE / run / "solverConfig.yaml").read_text()))


def drift(rows, key):
    s_ = np.array([r["step"] for r in rows], float); v = np.array([r[key] for r in rows], float)
    if s_[-1] - s_[0] < 20000:
        return None
    j = int(np.searchsorted(s_, s_[-1] - 20000))
    return 100 * (v[-1] / v[j] - 1) * 20000 / (s_[-1] - s_[j])


def reach_level(rows):
    """水準 (速度 plan §6.15、窓は run の中・step 0 は使わない) に初めて入った出力。非有限なら diverged_at"""
    rows = sorted([x for x in rows if x["step"] > 0], key=lambda x: x["step"])
    for n in range(len(rows)):
        pre = rows[: n + 1]; r = pre[-1]
        if r.get("nonfinite") or r.get("nonfinite_all"):
            return {"diverged_at": r["step"]}
        dr = [drift(pre, k) for k in KEYS[:3]]
        if r.get("deficit") is not None and None not in dr and abs(r["deficit"]) <= 0.1 and all(abs(x) <= 0.05 for x in dr):
            return {"step": r["step"], "drift": dr, "deficit": r["deficit"], "window": [x for x in pre if x["step"] >= r["step"] - 20000]}
    return None


def qs(window, tag):
    p = D / f"lvcmp_qs_{tag}.csv"
    with open(p, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["step"] + list(KEYS))
        for x in window:
            w.writerow([x["step"]] + [x[k] for k in KEYS])
    o = subprocess.run([sys.executable, str(HERE / "../../solver_density_cuda/tools/check_quasisteady.py"), "--series-csv", str(p), "--series-cols", ",".join(KEYS),
                        "--tail", "1.0", "--drift", "0.0005", "--min-snaps", "5"], capture_output=True, text=True)
    txt = (o.stdout + o.stderr).strip()
    ov = [l for l in txt.splitlines() if "OVERALL" in l]
    return {"rc": o.returncode, "overall": ov[-1] if ov else None, "tail": txt.splitlines()[-6:]}


def segverdict(run):
    p = HERE / run / "CONVERGENCE_SEGMENT.txt"
    t = p.read_text(errors="replace") if p.is_file() else ""
    m = re.search(r"->\s*([A-Z][A-Z ]+[A-Z])", t)
    return (m.group(1) if (m and "判定区間" in t) else None), t


def input_gates():
    g = []
    src_full = sha(HERE / SRC / SRC_RES)
    g.append(("出発の場の sha256", src_full.startswith(SRC_SHA16), src_full[:16]))
    ref_mesh = meshsha(HERE / SRC / "nozzle.h5")
    for name, (run, mk, _) in ARMS.items():
        d = HERE / run
        try:
            prov = (d / "RUN_PROVENANCE.txt").read_text(errors="replace")
            g.append((f"{run}: forge の sha256", SHA_NEW in prov, SHA_NEW[:16]))
            cp = json.loads((d / "COLD_PAIR.json").read_text())
            g.append((f"{run}: 出発の場", cp["parent"] == SRC and cp["parent_res"] == SRC_RES and cp["parent_res_sha256"] == src_full, cp["parent_res_sha256"][:16]))
            m = meshsha(d / "nozzle.h5")
            g.append((f"{run}: 格子の実体", m == ref_mesh and (d / "MESH_SHA.txt").read_text().strip() == m, m[:16]))
            diff_in = [fn for fn in SAME_FILES if not (d / fn).is_file() or sha(HERE / SRC / fn) != sha(d / fn)]
            g.append((f"{run}: 入力ファイル", diff_in == [], diff_in))
            c = cfg(run)
            exp_ = {**EXPECT, **ARM_CFG[name]}
            bad = {k: c.get(k) for k, v in exp_.items() if (k in c if v is None else c.get(k) != v)}
            g.append((f"{run}: 設定が腕 {name} の期待どおり", bad == {}, bad))
            log = (d / "forge_run.log").read_text(errors="replace")
            g.append((f"{run}: LAYOUT2 の分解・代入", "[line] factor 1 回目: モード LAYOUT2" in log and "[line] solve 1 回目: モード LAYOUT2" in log, None))
            g.append((f"{run}: 面エンタルピーの切替なし", "FORGE_DIAG_FACE_H_DOUBLE" not in log, None))
            vshow = "'lineViscCoupling' in 'time.deltaT': 3" in log
            mshow = "診断のマスク FORGE_LVC_TERMS=5 " in log
            okv = {"L7": vshow and "FORGE_LVC_TERMS=" not in log, "L5": vshow and mshow, "L0": ("'lineViscCoupling' in 'time.deltaT'" not in log) and "FORGE_LVC_TERMS=" not in log}[name]
            g.append((f"{run}: 値とマスクの表示が腕 {name} と合う", okv, None))
            lt = (d / "LVC_TERMS.txt").read_text().strip() if (d / "LVC_TERMS.txt").is_file() else None
            g.append((f"{run}: 台本が渡した FORGE_LVC_TERMS = {mk}", lt == mk, lt))
        except Exception as e:
            g.append((f"{run}: 入力のゲートを評価できない", False, f"{type(e).__name__}: {e}"))
    try:
        vr = subprocess.run([sys.executable, str(HERE / "lvcmp_pregate.py"), "--verify"], capture_output=True, text=True)
        g.append(("事前のゲートの記録が PASS で、証拠のハッシュが今のファイルと同じ", vr.returncode == 0, (vr.stdout + vr.stderr).strip()[-300:]))
        pg = json.loads((D / "lvcmp_pregate.json").read_text())
        dcs = pg.get("dump_config") or {}
        allowed = {"time/last/nStepOuter", "time/outStepInterval", "output/extraFields"}
        for name, (run, _, dk) in ARMS.items():
            c = cfg(run); dc = dcs.get(dk) or {}
            dd = sorted(k for k in set(c) | set(dc) if k not in allowed and c.get(k) != dc.get(k))
            g.append((f"{run}: 設定が同じ腕の書き出しと同じ (差は step 数・出力の間隔・extraFields だけ)", bool(dc) and dd == [], dd[:10]))
        for run in DUMPS:
            g.append((f"{run}: 書き出しが残っている", (HERE / run / "linedump" / "D.f64").is_file(), None))
    except Exception as e:
        g.append(("事前のゲートの記録を読めない", False, f"{type(e).__name__}: {e}"))
    return g


def assess(name):
    run = ARMS[name][0]
    info = {"run": run}
    p = HERE / run / "m9_watch.json"
    if not p.is_file():
        return "INVALID", {**info, "why": "見張りの記録が無い"}
    w = json.loads(p.read_text())
    st = w.get("status"); info["watch_status"] = st
    rows = w.get("rows") or []
    info["n_rows"] = len(rows)
    seg, _ = segverdict(run)
    info["segment_verdict"] = seg
    if st in ("EXEC_ERROR", "DATA_ERROR") or st is None:
        return "INVALID", {**info, "why": f"見張りの状態 {st}", "note": w.get("note")}
    rl = reach_level(rows)
    if st == "REACHED":
        info["reach_step_watch"] = w.get("reach_step")
        if rl is None or "step" not in rl or int(rl["step"]) != int(w.get("reach_step", -1)):
            return "INVALID", {**info, "why": f"到達の step が系列からの再計算と合わない ({rl and rl.get('step')} ≠ {w.get('reach_step')})"}
        if seg is None or seg.startswith("DIVERGED"):
            return "INVALID", {**info, "why": f"区間付きの収束 VERDICT が無いか DIVERGED ({seg})"}
        info.update(reach_step=int(rl["step"]), drift=rl["drift"], deficit=rl["deficit"], quasisteady=qs(rl["window"], name))
        return "REACHED", info
    if st == "DIVERGED":
        info["diverged_at"] = (rl or {}).get("diverged_at")
        return "DIVERGED", info
    if st == "CENSORED":
        if rl is not None:
            return "INVALID", {**info, "why": "CENSORED なのに系列は水準に入っている"}
        return "CENSORED", info
    return "INVALID", {**info, "why": f"未知の状態 {st}"}


def cmp(cls, x, y):
    cx, cy = cls[x], cls[y]
    if cx[0] == "REACHED" and cy[0] == "REACHED":
        d = cx[1]["reach_step"] - cy[1]["reach_step"]
        v = (f"{x} が速い" if d < 0 else f"{x} が遅い") if abs(d) >= THR else "判別不能"
        return {"verdict": v, "diff_steps": d, "x": cx[1]["reach_step"], "y": cy[1]["reach_step"]}
    return {"verdict": "比較しない", "x": cx[0], "y": cy[0]}


def main():
    rec = {"plan": "time_integration-line-viscous-jacobian-dt §6.3"}
    out = D / "lvcmp_judge.json"
    try:
        gates = input_gates()
    except Exception as e:
        gates = [("入力のゲートを評価できない", False, f"{type(e).__name__}: {e}")]
    cls = {}
    for n in ARMS:
        try:
            cls[n] = assess(n)
        except Exception as e:
            cls[n] = ("INVALID", {"why": f"例外 {type(e).__name__}: {e}"})
    rec["gates"] = [{"check": g[0], "ok": bool(g[1]), "value": g[2]} for g in gates]
    rec["gates_ok"] = all(g[1] for g in gates)
    rec["arms"] = {n: {"class": c[0], **c[1]} for n, c in cls.items()}
    invalid = [n for n in ARMS if cls[n][0] == "INVALID"]
    if not rec["gates_ok"] or invalid:
        rec["VERDICT"] = f"INVALID (ゲート不合格 {sum(not g[1] for g in gates)} 件、INVALID の腕 {invalid}): 判定しない"
    else:
        c70, c75, c50 = cmp(cls, "L7", "L0"), cmp(cls, "L7", "L5"), cmp(cls, "L5", "L0")
        rec["comparisons"] = {"L7_vs_L0": c70, "L7_vs_L5": c75, "L5_vs_L0": c50}
        s7, s0 = cls["L7"][0], cls["L0"][0]
        if s0 != "REACHED":
            main_v = "判別不能 (生産 L0 が水準に届かない)"
        elif s7 != "REACHED":
            main_v = f"NO-GO (全部入り L7 は上限 50 で水準に届かない: {s7})"
        elif c70["verdict"] == "L7 が速い":
            main_v = "全部入りが生産より早く水準に届く (この条件)"
        elif c70["verdict"] == "L7 が遅い":
            main_v = "NO-GO (全部入りは生産より遅い)"
        else:
            main_v = "全部入りと生産は区別できない"
        rec["main"] = main_v
        rec["VERDICT"] = f"主 (L7 と L0): {main_v}; L7 と L5: {c75['verdict']}; L5 と L0: {c50['verdict']}"
    out.write_text(json.dumps(rec, ensure_ascii=False, indent=1, default=str))
    for g in gates:
        if not g[1]:
            print("ゲート不合格:", g[0], g[2])
    for n, c in cls.items():
        print(n, c[0], {k: c[1].get(k) for k in ("run", "reach_step", "diverged_at", "segment_verdict", "why")})
    print(rec["VERDICT"])
    return 0 if rec["gates_ok"] and not invalid else 1


if __name__ == "__main__":
    sys.exit(main())
```

## 参考: `case/45.isobutane_m6_d155/m9_watch.py`

```
"""§5.1 #9 の見張り v2 (plan time_integration-line-implicit-speed §6.11、2026-10-10、codex plan-2 の指摘の反映)。
usage: python3 m9_watch.py <run> --phase line|point --budget <最大 step> → <run>/m9_watch.json (状態と系列をまとめて原子的に保存)
- 書き終わった res_<N>.h5 (N > 0、60 s 以上前に書かれたもの、または forge が終わった後) を 1 つずつ `cold_series.snapshot` で 1 行にする。
  VALUE の全ての浮動小数点の量の非有限も数える (cold_series は ρ・u・T・k だけ)。読めないものは 3 回まで再試行、それでも読めなければ DATA_ERROR。
- 終わりの判定 (step 0 は使わない): 水準 = |欠損| ≤ 0.1 kg/s かつ θ_r(40/70/94) のドリフト ≤ 0.05 %/2 万 step (窓はこの run の中)。
  --phase point では加えて θ_r(40/70/94)・Q_w が参照 R (m9_ref.json) から相対 0.1 % 以内 (E2)。満たした最初の出力で REACHED。
- 状態: REACHED / CENSORED (forge が上限まで正常に回って届かず) / DIVERGED (非有限) / EXEC_ERROR (forge が上限より前に異常終了) / DATA_ERROR。
- REACHED・DIVERGED では cwd がこの run の forge に SIGTERM、180 s 待って残れば SIGKILL、終了を確かめてからそれより後の出力を消し、到達の出力の sha256 を記録する。
- 消す: 5 万 step ごと・最新 5 つ (準定常の判定の窓)・到達の出力以外の res_<N>.h5 と同じ N の res_* (状態を保存した後に消す)。
- 再開: 既存の m9_watch.json を読み、最終の状態なら何もしない。"""
import argparse, hashlib, json, os, re, signal, subprocess, sys, time
from pathlib import Path
import h5py
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cold_series as CS  # noqa: E402
import cold_xcheck as XC  # noqa: E402
def _emergency_stop(run_name, why):
    """初期化の失敗でも、cwd がその run の forge を止めて状態を残す (codex plan-4 M1)。"""
    import signal as _sg, subprocess as _sp, time as _t
    rd = os.path.realpath(str(HERE / run_name))
    def _p():
        r_ = []
        for q in _sp.run(["pgrep", "-x", "forge"], capture_output=True, text=True).stdout.split():
            try:
                if os.path.realpath(f"/proc/{q}/cwd") == rd: r_.append(int(q))
            except OSError: pass
        return r_
    for sig_, w_ in ((_sg.SIGTERM, 180), (_sg.SIGKILL, 60)):
        for q in _p():
            try: os.kill(q, sig_)
            except OSError: pass
        t0_ = _t.time()
        while _p() and _t.time() - t0_ < w_: _t.sleep(3)
    try:
        sf_ = HERE / run_name / "m9_watch.json"
        st_ = json.loads(sf_.read_text()) if sf_.exists() else {"run": run_name, "rows": []}
        st_.update(status="EXEC_ERROR", note=f"見張りの初期化の失敗: {why}")
        sf_.write_text(json.dumps(st_, indent=1, ensure_ascii=False))
    except Exception: pass
    print(f"[m9_watch] {run_name}: 初期化の失敗 ({why}) — forge を止めた", flush=True); sys.exit(1)

ap = argparse.ArgumentParser(); ap.add_argument("run")
ap.add_argument("--phase", choices=("line", "point", "line_e2"), required=True)   # line_e2 = ラインのまま E (水準 + E2) まで (§6.11 の追加の腕 L5L)
ap.add_argument("--budget", type=int, required=True)
ap.add_argument("--consec", type=int, default=1, help="水準 (point と line_e2 は E2 も) を何出力続けて満たしたら到達とするか (既定 1 = 従来。§6.18 は 2、diagnostician 2026-10-10)")
ap.add_argument("--inherit", default=None, help="line_e2: 同じ構成の分岐元の run。その到達の出力までの系列を step − 到達 (≤ 0) で引き継ぐ (codex plan-4 M2)")
try:
    a = ap.parse_args()
except SystemExit:
    if len(sys.argv) > 1 and not sys.argv[1].startswith("-"): _emergency_stop(sys.argv[1], "引数の解釈")
    raise
try:
    run = HERE / a.run
    SF = run / "m9_watch.json"
    KEYS = ("theta_r_40", "theta_r_70", "theta_r_94", "Q_w")
    REF = None
    if a.phase in ("point", "line_e2"):                    # E2 の参照は数値 4 つだけ (説明の文字列は使わない、codex plan-3 M2)
        _r = json.loads((HERE / "m9_ref.json").read_text())
        REF = {k: float(_r[k]) for k in KEYS}
        if not all(np.isfinite(v) and v > 0 for v in REF.values()): raise ValueError("参照が有限・正でない")
    if (a.phase == "line_e2") != (a.inherit is not None): raise ValueError("--inherit は line_e2 のときだけ、必ず付ける")
    FINAL = ("REACHED", "CENSORED", "DIVERGED", "EXEC_ERROR", "DATA_ERROR")
    if SF.exists():
        st = json.loads(SF.read_text())
    else:
        st = {"run": run.name, "phase": a.phase, "budget": a.budget, "consec": a.consec, "status": "running", "rows": [], "tries": {}}
        if a.inherit:
            src = json.loads((HERE / a.inherit / "m9_watch.json").read_text())
            if src.get("status") != "REACHED": raise ValueError(f"分岐元 {a.inherit} が REACHED でない")
            nb = int(src["reach_step"])
            st["inherit"] = {"run": a.inherit, "reach_step": nb, "reach_sha256": src.get("reach_sha256")}
            st["rows"] = [dict(r, step=r["step"] - nb, inherited=True) for r in src["rows"] if 0 < r["step"] <= nb]
    if a.consec < 1: raise ValueError("--consec は 1 以上")
    if st.get("phase") != a.phase or st.get("budget") != a.budget or st.get("consec", 1) != a.consec:   # 終端の状態でも先に照合する (codex plan-7 M3)
        raise ValueError("状態の phase/budget/consec が違う")
    if st["status"] in FINAL:
        print(f"[m9_watch] {run.name}: 既に {st['status']}"); sys.exit(0)
except SystemExit:
    raise
except Exception as e:
    _emergency_stop(a.run, repr(e))
yb_x, yb = XC.common_yb()

def save():
    tmp = SF.with_suffix(".tmp"); tmp.write_text(json.dumps(st, indent=1, ensure_ascii=False)); os.replace(tmp, SF)
def pids():
    out = subprocess.run(["pgrep", "-x", "forge"], capture_output=True, text=True).stdout.split(); r = []
    for p in out:
        try:
            if os.path.realpath(f"/proc/{p}/cwd") == str(run.resolve()): r.append(int(p))
        except OSError: pass
    return r
def stop():
    ps = pids()
    for p in ps:
        try: os.kill(p, signal.SIGTERM)
        except OSError: pass
    t0 = time.time()
    while pids() and time.time() - t0 < 180: time.sleep(5)
    for p in pids():
        try: os.kill(p, signal.SIGKILL)
        except OSError: pass
    t0 = time.time()
    while pids() and time.time() - t0 < 60: time.sleep(2)
    return not pids()
def steps():
    return sorted(int(m.group(1)) for f in run.glob("res_*.h5") if (m := re.fullmatch(r"res_(\d+)\.h5", f.name)))
def drift(rows, key):
    s_ = np.array([r["step"] for r in rows], float); v = np.array([r[key] for r in rows], float)
    if s_[-1] - s_[0] < 20000: return None
    j = int(np.searchsorted(s_, s_[-1] - 20000)); return float(100 * (v[-1] / v[j] - 1) * 20000 / (s_[-1] - s_[j]))
def nonfinite_all(n):
    c = 0
    with h5py.File(run / f"res_{n}.h5", "r") as h:
        for k, d in h["VALUE"].items():
            x = np.asarray(d[:])
            if x.dtype.kind == "f": c += int(np.count_nonzero(~np.isfinite(x)))
    return c
def finish(status, **kw):
    st.update(status=status, **kw); save(); print(f"[m9_watch] {run.name}: {status} {kw}", flush=True)
def cleanup(keep_extra=()):
    done = sorted(r["step"] for r in st["rows"] if r["step"] > 0 and not r.get("inherited"))
    keep = set(done[-1:]) | set(keep_extra)               # 最新と到達の出力だけ残す (準定常は系列から判定するので不要、ディスクのため。2026-10-10)
    for m in done:
        if m in keep: continue
        for f in list(run.glob(f"res_{m}.h5")) + list(run.glob(f"res_*_{m}.h5")) + list(run.glob(f"res*_{m}.xmf")): f.unlink()

def done_at(rows):
    """保存済みの系列の最初の終わりを ("REACH", step) か ("DIV", step) で返す (引き継いだ行の step は ≤ 0)。無ければ None。
    REACH は水準 (と E2) を --consec 個の連続した出力で満たした、その最後の出力 (§6.18)。"""
    run_ok = 0
    for n_ in range(len(rows)):
        r_ = rows[n_]
        if r_.get("nonfinite") or r_.get("nonfinite_all"): return ("DIV", r_["step"])
        dr_ = [drift(rows[: n_ + 1], k) for k in KEYS[:3]]
        e2_ = {k: 100 * (r_[k] / REF[k] - 1) for k in KEYS} if REF else None
        ok_ = r_.get("deficit") is not None and None not in dr_ and abs(r_["deficit"]) <= 0.1 and all(abs(x) <= 0.05 for x in dr_) \
              and (e2_ is None or all(abs(v) <= 0.1 for v in e2_.values()))
        run_ok = run_ok + 1 if ok_ else 0
        if run_ok >= a.consec:
            return ("REACH", r_["step"])
    return None
def conclude(n):
    """終わりの出力 n で止める (何度呼んでも同じ結果): forge を止め、後の出力を消し、sha256 を記録して REACHED。"""
    ok = stop()
    if not ok: finish("EXEC_ERROR", note="forge を止められない", step=n); sys.exit(1)
    for m in steps():
        if m > n:
            for f in list(run.glob(f"res_{m}.h5")) + list(run.glob(f"res_*_{m}.h5")) + list(run.glob(f"res*_{m}.xmf")): f.unlink()
    if n <= 0:                                           # 分岐点 (引き継いだ行) で既に E を満たす: 追加の step は 0
        finish("REACHED", reach_step=0, note="分岐点で既に E を満たす (追加 0 step)"); sys.exit(0)
    if not (run / f"res_{n}.h5").exists(): finish("DATA_ERROR", note=f"到達の出力 res_{n}.h5 が無い", step=n); sys.exit(1)
    sha = hashlib.sha256((run / f"res_{n}.h5").read_bytes()).hexdigest()
    finish("REACHED", reach_step=n, reach_sha256=sha); cleanup(keep_extra=(n,)); sys.exit(0)
# 再開: 保存済みの系列で先に判定し直す (codex plan-3 M4)
_d = done_at(st["rows"])
if _d is not None:
    if _d[0] == "DIV": stop(); finish("DIVERGED", fail_step=_d[1]); sys.exit(0)
    conclude(_d[1])

T0 = time.time()
def main_loop():
  while True:
      alive = bool(pids())
      have = {r["step"] for r in st["rows"]}
      new = [n for n in steps() if n > 0 and n not in have and (not alive or time.time() - (run / f"res_{n}.h5").stat().st_mtime > 60)]
      progressed = False
      for n in new:
          try:
              rec = CS.snapshot(run, n, yb_x, yb); rec["nonfinite_all"] = nonfinite_all(n)
          except Exception as e:
              st["tries"][str(n)] = st["tries"].get(str(n), 0) + 1; save()
              print(f"[m9_watch] {run.name} {n}: 読めない ({e}) — {st['tries'][str(n)]} 回目", flush=True)
              if st["tries"][str(n)] >= 3: finish("DATA_ERROR", step=n); stop(); sys.exit(1)
              break
          st["rows"].append(rec); st["rows"].sort(key=lambda r: r["step"]); progressed = True
          rows = st["rows"]
          dr = [drift(rows, k) for k in ("theta_r_40", "theta_r_70", "theta_r_94")]
          e2 = {k: 100 * (rec[k] / REF[k] - 1) for k in KEYS} if REF else None
          rec["drift"] = dr; rec["e2"] = e2
          print(f"[m9_watch] {run.name} {n}: 欠損 {rec['deficit']:.4f}、ドリフト {dr}、E2 {e2}、非有限 {rec['nonfinite']}/{rec['nonfinite_all']}", flush=True)
          if rec["nonfinite"] or rec["nonfinite_all"]:
              save(); ok = stop(); finish("DIVERGED", fail_step=n, stopped=ok); sys.exit(0)
          _dn = done_at(st["rows"])                       # 連続の数え方は再開時と同じ関数で (§6.18 の --consec)
          if _dn is not None and _dn[0] == "REACH":
              save(); conclude(_dn[1])
          save(); cleanup()
      started = (run / "RUN_RC").exists() or time.time() - T0 > 600     # forge の起動前に見張りが始まっても待つ (起動の猶予 10 分)
      if started and not progressed and not alive and not [n for n in steps() if n > 0 and n not in {r["step"] for r in st["rows"]}]:
          rcf = run / "RUN_RC"; rc = rcf.read_text().strip() if rcf.exists() else "?"
          last = max((r["step"] for r in st["rows"]), default=0)
          nan = sorted(run.glob("res_nan_*.h5"))
          if rc == "0" and last >= a.budget: finish("CENSORED", last_step=last)
          elif nan: finish("DIVERGED", fail_step=int(nan[0].stem.split("_")[-1]), note="forge の detectNaN が止めた")
          else: finish("EXEC_ERROR", rc=rc, last_step=last)
          sys.exit(0)
      time.sleep(30)

try:
    main_loop()
except SystemExit:
    raise
except Exception as e:                                   # 見張りの例外でも forge を止めて状態を残す (codex plan-3 M2)
    import traceback; traceback.print_exc()
    ok = stop(); finish("EXEC_ERROR", note=f"見張りの例外: {e!r}", stopped=ok); sys.exit(1)
```

## 参考: `notes/reviews/2026-10-11-lvcdt-result-diagnose.md`

```
# 諮問: 方向別 dt の上限の A/B の結果の解釈と次の設計 (diagnostician、Fable)

- **plan**: [`plans/active/time_integration-line-viscous-jacobian-dt.md`](../../plans/active/time_integration-line-viscous-jacobian-dt.md) §6
- **諮問先**: `diagnostician` サブエージェント (Fable)。設計の諮問 ([2026-10-11-lvcdt-design-diagnose.md](2026-10-11-lvcdt-design-diagnose.md)) と同じエージェントに、文脈を引き継いだまま結果を送った (ユーザの指定「diagnostitianにはかって」の続き)。
- **date**: 2026-10-11
- **主セッションの照合**: `rms_roOmega` の末尾 500 step の傾き (+0.094・+0.093 桁)、全腕で step 1 に 3 倍を超えること (`lvcdt_judge.json` の `tail_slope500`・`first_step_over_3x`) は記録と一致した。

本文はサブエージェントの最終報告をそのまま転記 (採否は plan §6.1・§6.2 に書く)。

---

## 結論

**結論**: 「回避」は事前登録 (dt plan §6 分岐 1) の文のまま採る。次は §5.1 #2 を「上限 50・同一バイナリで、値 0 (生産) / 値 3・マスク 5 / 値 3・マスク 7 の 3 腕を、速度 plan §6.15 の**水準** (|欠損| ≤ 0.1 kg/s かつ θ_r(40/70/94) のドリフト ≤ 0.05 %/2 万 step) に届く step 数で比べる」設計にする。2000 step の窓では比べない。機構 (H1〜H3) は今回の結果で 1 つも消えていないので、調査の順は #2 の結果の後に決める (下の Q3)。

## 採否表

| 重大度 | 対象 | 採否 | 根拠と対案 |
| --- | --- | --- | --- |
| Major | 主判定「回避」(`lvcdt_judge.json` VERDICT) | **採用 (登録の文のまま)** | A 2 本 DIVERGED (366・481 step、`ro`、列 42〜52)、B 2 本 FINITE、ゲート 60 合格・INVALID なし。事前のゲート PASS で、介入が dt だけ (K 3 本ビット一致、用量 Δτ_B = min(Δτ_A, 50Δτ_C) 相対 3e-16、第一内部節点 r 171〜184)。 |
| Major | 「マスク 7 は上限 50 で安定」「全部入りが成り立った」 | **却下 (書かない)** | B は 2 本とも `check_convergence --segment` NOT CONVERGED (低下桁数 ro 0.4・roK 0.0)。さらに **`rms_roOmega` の末尾 500 step の傾きが +0.094・+0.093 桁/500 step で上昇** (`lvcdt_judge.json` b1/b2 `tail_slope500`)。2000 step の窓の外で非有限化する可能性は排除されていない。plan には「2000 step 有限・2 本とも未収束・`rms_roOmega` は末尾で上昇」と書く。 |
| Major | 帳簿の第一層 |Δρ/ρ| (A 0.55、B 0.011) からの解釈 | **記録だけ (登録どおり)** | 比 50 は r (171〜184) に比例しないが、A は非線形に飽和する過渡なので「同じ過渡を遅く通っただけか」は決められない (plan §6 の M3 の留保のまま)。累積の擬似時間は未測定。 |
| Minor | `rms_roOmega` が全腕で step 1 に 3 倍超 (B は最大 4.6 倍) | **共通の restart 過渡として記録** | A・B とも `first_step_over_3x` = 1。restart で乱流場が忠実に戻らない既知の挙動 (vis_turb/壁 ω) の型で、腕の判別には使わない。 |
| Minor | 参考に挙げた値 0・上限 50 の 60000 step・L5 115000 step | **参考に限る** | 別バイナリ (古い・float)・別の終わりの条件。同一バイナリの対照は #2 で取り直す。 |
| Minor | 「熱拡散数 ≈ 0.2 になった」「SST と平均流を分けた」 | **書かない** | r ≥ 10 は Δτ の比だけで κ/(V/Δτ) は測っていない。dt は SST・化学種にも効く (登録 §6 の「言えないこと」)。 |

## Q1: plan に書いてよい文 / いけない文

- 書いてよい: 「この条件 (段 ③ FP64・`run_0183` res_100000・値 3・マスク 7・キー 5・cfl 4)・2000 step では、方向別 dt に上限 50 を付けると非有限化を避けた。上限なしは 2 本とも 366・481 step で非有限 (列 42〜52、`ro`)。上限 50 ではどの残差列も出発の 10 倍に達しなかった (最大 `rms_roOmega` 4.6 倍、`rms_roUy` 2.4 倍)。2 本とも NOT CONVERGED で `rms_roOmega` は末尾 500 step で +0.09 桁上昇。帳簿の第一層 |Δρ/ρ| は 0.55 → 0.011 (記録)。」
- 書いてはいけない: 「全部入りが成り立った/安定/収束した」「H1 (熱拡散数 ≫ 1) を支持」「上限 50 が必要十分」(用量は 1 水準のみ)「平均流の時間項の効果」(SST・化学種を分けていない)「同じ過渡を遅く通った」。

## Q2: §5.1 #2 の設計 (事前登録の案)

- **問い**: 上限 50 で、全部入り (値 3・マスク 7) は生産 (値 0) と回避策 (マスク 5) より**水準に早く届くか**。2000 step では全腕が横ばいで比べられない (B の低下桁数 0.0〜0.8)。残差の桁数・傾きは床と有界振動に支配されるので主指標にしない (速度 plan §6.14・§6.15 で L0/L5 とも NOT CONVERGED のまま水準に届いた実績)。
- **腕 (3 本、各 1 run)**: L0c = 値 0・キー 5・方向別・上限 50 (生産の構成、同一バイナリで取り直す)、L5c = 値 3・`FORGE_LVC_TERMS=5`・上限 50、L7c = 値 3・マスク 7・上限 50。共通: 段 ③ FP64 129de3f4、`run_0183` res_100000、cfl 4・緩和 0.7・sweep 5・ISP 0。値 0 の対照は必須 (全部入りの目的は生産を超えること。L0 の既存系列は別バイナリなので使わない)。値 2 は足さない。
- **長さ**: 最大 150000 step、出力 2500 ごと (水準のドリフト窓 2 万 step に ≥ 8 点)。水準に入った最初の出力で止める (速度 plan の `m9_watch.py`/`levelreach.py` の手順。SIGTERM → 確認 → 不要な出力を消す)。1 腕 ≈ 1.2 h × FP64 の係数 (単価は `m9_unit.py` で同一バイナリで測り直す)。GPU の予定に書いてから逐次。
- **事前のゲート** (1 step の書き出し 3 本 + L7c の再実行): I・R・S は faceh §6.12 と同じ。V: (a) `dt_vol` が 3 腕でビット一致 (`lineViscCoupling` は `setDT_d.cu` を触らない)、(b) L7c − L5c の K の差は faceh §6.12 の V(c) (行 4 の熱伝導の項 κ∂T_j/∂Q_j と一致、行 0〜3 不変) を流用、(c) L5c − L0c は K・D のライン面の項が薄層の式どおり (監査の記録で照合、無理なら「係数の変化の成立」だけ)、(d) 用量 Δτ = min(Δτ_dir, 50Δτ_point) は今回の記録と同じ節点でビット一致。
- **判定**: 到達 step N の差が **15000 step (出力 6 回ぶん) 以上**で「速い/遅い」、未満は判別不能。forge は同一入力でも再実行で数 % 揺れるので、判別の閾値を事後に下げない。DIVERGED の腕は「その構成は上限 50 で水準に届かない」と記録 (比較しない)。NOT CONVERGED は到達の比較に限る (速度 plan §6.11 の規則)。副指標 (記録): 到達の窓 [N − 2 万, N] の θ_r・Q_w の `check_quasisteady --series-csv`、`check_convergence --segment`、各列の末尾の傾き。
- **分岐**: L7c が L0c より速い → 全部入りの利点を支持 (この条件)、R の用量 (20/100/200) に進む前に result 段へ。L7c ≈ L0c (判別不能) → 「全部入りは生産と区別できない」で閉じる (速度 plan §6.15 の L5 と同じ結論の型。長い run を足さない)。L7c が遅い/DIVERGED → NO-GO。L5c は「熱伝導の K が効くか」の分離に使う (L7c と L5c の差 ≥ 15000 で判別)。
- **やらないこと**: 2000 step の残差で優劣を言う、既存の L0/L5 系列 (別バイナリ) との直接比較、到達の閾値や終わりの条件の事後変更、point 仕上げの追加。

## Q3: H1〜H3

- 言えること: 時間項の強化 (上限 50) が破綻を 2000 step 以内では止める。帳簿の先行変数 (roUy 4・roOmega 14・ro 16・roe 34 step) は faceh §6.13 と再現した。
- 言えないこと: H1・H2・H3 のどれも消えていない (時間項の強化は 3 つの機構すべてに効く。登録どおり主判定は機構に使わない)。`rms_roOmega` の上昇は SST の関与が残ることの観測で、機構の判別ではない。
- 順 (提案): まず #2 (実用の判断。L7c が長期で DIVERGED になればそれ自体が「上限 50 は崖の近く」の情報)。機構は #2 の後、**H2 (#4、実残差の方向微分と薄層 K の整合を上限 50 の状態で取り直す)** → H3 (#3、行列のスケール済み反復写像) → H1 の用量 (R = 20/100/200 の 2000 step、閾値 R* の有無)。用量は #2 で L7c が DIVERGED のときだけ前倒しする。

## 呼び出し側の前提への異議

無し (登録と判定の文は整合している)。ただし「上限 50 の 2 本とも 10 倍に達しなかった」は `rms_roOmega` の末尾上昇と併記しないと安定の印象を与える。

## 不足情報

無し (判定 JSON・事前のゲート・plan §6 で足りた)。plan 未反映 (§6.2・§5.1 #2 の事前登録は呼び出し側)。
```

## 出力形式

1. 冒頭に **判定 (GO / GO-with-changes / NO-GO)** と 3 行以内の要約。
2. 指摘一覧 (Critical → Major → Minor の順、番号付き。各項目に根拠と対案)。
3. 推奨 (1 つに絞る)。
4. 末尾に `指摘数: Critical N / Major N / Minor N` の 1 行。
