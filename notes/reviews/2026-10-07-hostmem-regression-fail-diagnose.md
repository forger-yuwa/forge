# codex 諮問 (diagnose): hostmem-regression-fail

- **brief**: [`notes/reviews/briefs/2026-10-07-hostmem-regression-fail.md`](../../notes/reviews/briefs/2026-10-07-hostmem-regression-fail.md)
- **plan**: [`plans/active/architecture-solver-host-memory.md`](../../plans/active/architecture-solver-host-memory.md)
- **date**: 2026-10-07
- **commit**: `cd725134` (feature/sern-design)
- **codex**: effort `xhigh`, 4.7 min, rc=0
- **結論**: ?
- **extra**: `case/66.hostmem_regression/README.md`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 諮問事項 | 判断 | 理由・対案 |
|---|---|---|
| 登録判定を FAIL のまま残す | **採用** | 原判定を保存し、比較器の欠陥を検証する追加診断を別記する。現時点で R1–R3 を合格にはしない。 |
| (i) 初期化後の全デバイス配列ハッシュで決着 | **却下〔Major〕** | 未初期化領域と、初期化中の非決定的演算を含む。一致しても、その後の挙動が同等である証明にはならない。 |
| (ii) 各10回に増やして順位和検定 | **今は却下** | 比較指標の修正が先。「有意差なし」は同等性の証明にならず、必要反復数も10回とは決められない。 |
| (iii) 組立直後・`atomicAdd` 前のハッシュ | **その定義では却下〔Major〕** | 組立にも、その前の初期化にも `atomicAdd` がある。RMS 集計前の採取は可能だが、非決定性より前の状態にはならない。 |

**結論:** 登録 FAIL を維持し、次は追加の forge 実行ではなく、既存6反復の比較を「ペアごとの相対差」から「絶対 L∞ 差」に替える A/B を行う。

**第1仮説:** (b) の巨大な超過には、比較器の非対称な正規化が主要因として入っている。**確度: 高**。比較器の欠陥自体は再現確認済みで、7量すべてを説明するかは未確認。

**根拠〔Major〕:**

[compare_runs.py:134](/home/sano/work/forge-sern-design/case/66.hostmem_regression/compare_runs.py:134) は `m(A,B)=max|A−B|/max|A|` を使う。一方、[同:383](/home/sano/work/forge-sern-design/case/66.hostmem_regression/compare_runs.py:383) は同ビルド内を `i<j` の片方向だけ、ビルド間を base→new の全組合せで評価する。**S と D の分母が揃っておらず、反復の並び順に判定が依存する。**

実際の `metric()` を読み込んだ最小再現は次の結果だった。これは CFD の実測値ではなく、比較器の検証用入力である。

| base の3値 | new の3値 | S | D | 判定 |
|---|---|---:|---:|---|
| `[100, 2, 1]` | `[100, 2, 1]` | 0.99 | 99 | FAIL |
| `[100, 2, 1]` | `[1, 2, 100]` | 99 | 99 | PASS |

**同一の値集合でも FAIL になり、new の順序変更だけで PASS に変わる。** 対案は、同じデータセットの全ペアに共通の尺度を使うこと。今回の `D≤2S` の診断なら、正規化を外した絶対差で十分である。

実データとの接続もある。`case/66.hostmem_regression/run_0034_c44dual_restart100_base_r1/` を B1、`run_0069_c44dual_restart100_base_r3/` を B3 とする[比較報告:104](/home/sano/work/forge-sern-design/case/66.hostmem_regression/results/2026-10-07_base9c9f623c_vs_new93e55957/c44dual_restart100.txt:104)では、`condClampCorrQ_0` の `max|B1|=1.433e10`、最悪ペア B3–N2 の絶対差 `7.938e6`、D=`4.386e3`。したがって、報告の丸め値から `max|B3|≈1.81e3` と逆算できる。

ここから逆向きの `m(B3,B1)` は少なくとも約 `7.92e6`。**base の順序を逆にするだけで S が増え、D を変えずにこの量は PASS へ反転する計算になる。** これは報告値からの推論であり、生の配列による再比較は未実施。

**反証条件:** 共通尺度で再比較しても `D_abs>2S_abs` が残る量については、「正規化だけで超過を説明できる」という仮説を棄却する。

**第2仮説:** 既存の演算順序による非決定性と、凝縮診断量自身の微小分母が差を増幅している。**確度: 中**。`condClampCorrQ_0` は絶対補正量ではなく、`|ΔQ|/max(|Q_before|,1e−30)` の最大値を保持する（[condensationRealizability_d.cuh:313](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:313)）。`1e23` という数値だけでは、保存量の大きな誤差を意味しない。ただし、変更起因の差が無いことも証明しない。

**判別 A/B:** 変更点は比較尺度の1点だけ。**追加計算は0 step**。現在の base/new 各3本、同じ保存時点・同じ配列・同じCSV行を使う。

- **A:** 現行の `m(A,B)` による登録済み結果。
- **B:** `d(A,B)=max|A−B|`。同ビルド内6対の最大を `S_abs`、ビルド間9対の最大を `D_abs` とする。FAIL 7量だけでなく、従来 PASS の量も再評価する。
- **事前登録する追加診断条件:** 全入力が有限で、shape・列・行キーが対応すること。`D_abs≤2S_abs`、`S_abs=0` なら `D_abs=0`。整数・構造情報は厳密一致。反復の並べ替えと base/new 交換で結果が変わらないこと。
- **Bで超過が消える場合:** その量の登録 FAIL は尺度依存で説明できる。「追加診断で反復内差の2倍以内」と記録する。
- **Bでも超過する場合:** 比較器だけでは説明できない。残った量について、変更起因の差を候補に戻す。

この追加診断の PASS は、旧判定の書換えや「非決定性だけだった」という証明には使わない。

**やらない方がよいこと:**

- run の順序を選んで PASS にすること。
- FAIL 量を事後的に除外すること、反復を足して S が増えるまで続けること。
- 全配列をゼロ初期化してハッシュを揃えること。比較対象の挙動を変更してしまう。
- ハッシュ一致、または順位和検定の非有意だけで R1–R3 を承認すること。

**呼び出し側の前提への異議:**

1. **〔Major〕比較器はゼロ基準で誤 PASS も返す。**  
   [compare_runs.py:146](/home/sano/work/forge-sern-design/case/66.hostmem_regression/compare_runs.py:146) はゼロでない差をゼロ基準で割る場合に `inf` を返し、[同:402](/home/sano/work/forge-sern-design/case/66.hostmem_regression/compare_runs.py:402) は `inf≤2×inf` を真として扱う。最小再現で base=`[0,1,1]`、new=`[100,100,100]` が PASS になった。今回の各量で発生したかは未確認だが、既存 PASS も再評価対象にする必要がある。対案は上記の絶対差と、非有限・比較不能を合格にしない明示的な処理。

2. **〔Major〕「初期化後は全配列が決定的」は成立しない。**  
   [variables.cpp:315](/home/sano/work/forge-sern-design/solver_density_cuda/variables.cpp:315) は一部の名前だけをゼロ初期化する。さらに初期化中の [main.cpp:1828](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:1828) から呼ぶ境界勾配には、実際に到達する [calcGradient_d.cu:251](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/calcGradient_d.cu:251) の `atomicAdd` がある。  
   ハッシュを使うなら、**比較時点・初期化済みの有効領域・対象配列を先に定義する**必要がある。`c_d` 外の `qacc_d` も存在する（[variables.hpp:41](/home/sano/work/forge-sern-design/solver_density_cuda/variables.hpp:41)）。提示された配列集合の一致だけでは、全実行状態を覆ったとは言えない。

3. **〔Major〕(a) は「初回組立」と「最初の物理step全体」を混同している。**  
   [compare_runs.py:198](/home/sano/work/forge-sern-design/case/66.hostmem_regression/compare_runs.py:198) は最初の step 番号の全行を拾う。dual-time は内反復ごとに組立・更新するため、22行を初期状態の検査として扱えない（[main.cpp:2679](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:2679)）。  
   今後は、初期化済みの決定的状態のビット一致と、最初の組立後の残差比較を別ゲートにする。残差は `(step, inner, phase)` を固定し、方程式別の絶対・相対許容差を事前登録する。許容差は独立した base の校正と必要な検出幅から決め、「観測した2値のどちらか」や全列共通の最大 ULP 幅を使わない。

4. **〔Major〕「差はすべて非決定性」はまだ解釈である。**  
   初期出力の一致とカーネルソースの不変更は支持材料だが、十分条件ではない。また `roUz`・`CHECKPOINT/roUzN` は名目ゼロでも保存量であり、FAIL を診断量だけとして除外できない（[c20cell_dual.txt:73](/home/sano/work/forge-sern-design/case/66.hostmem_regression/results/2026-10-07_base9c9f623c_vs_new93e55957/c20cell_dual.txt:73)）。対案は、ゼロ近傍を扱える尺度で全量を評価すること。

**不足情報:** ローカルには対象 `run_*` と `CONVERGENCE_VERDICT.txt` がなく、AWS の生の HDF5・残差CSVを独立再集計できていない。今回確認したのはソース、持ち帰り比較報告、比較器の最小再現である。統計的同等性を最終条件にする場合の量別許容幅も未定義。

ファイル変更・forge 実行はしていない。**plan 未反映**。依頼どおり、呼び出し側が `architecture-solver-host-memory.md` §6・§5.1 に反映するための診断として返す。
