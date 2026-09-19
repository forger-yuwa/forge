# codex レビュー: case-hypersonic-gap-heating-validation (plan)

- **plan**: [`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `ff68487e` (feature/sern-design)
- **codex**: effort `high`, 6.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M8/m0
- **focus**: 1 巡目 (NO-GO, C1/M8/m1) を全件採用した改訂稿。採用が実質的か (特に C1: T1b の降格と T2 = TN D-8233 の新設で case/49 への受け渡しが正当化できるか、M2 幾何、M4 開口面積平均と共通分母、M7 不確かさ伝播、M8 定常化、M9 量別格子と離散収支) を検証してほしい。新設の T2 (NASA TN D-8233: 単独すきま W=2.29mm D=45.7mm, M10.3, 乱流トンネル壁 BL, T_surf/T_gap 1.0-1.6) が 2D 再現可能か、入口 BL プロファイルの与え方と h_smooth の定義が閉じているかも見てほしい
- **extra**: `notes/investigations/hypersonic-gap-cavity-heating-survey.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: NO-GO**

T1 の幾何・開口面積平均・共通分母、T1b の探索への降格は実質的な改善です。  
ただし T2 の二次元近似・入口条件・正規化がまだ閉じておらず、case/49 のモデル誤差を定量化する根拠にはできません。  
ファイルは変更していません。以下の提案は **plan 未反映**です。

`plans/README.md` と `plans/accepted/` を確認した範囲では、本検証を置き換える完了計画はありません。平面 node・SLAU・block-DPLUR の採用も現行方針に整合します。今回はソルバ共通部の変更を計画していないため、cell・周期・軸対称の回帰試験追加は要求しません。

TH76 はリポジトリの `papers/gap_heating/` にはありませんでしたが、[作業用PDF](/tmp/claude-1000/-home-sano-work-forge/25329893-9712-4618-8f8d-d3ff183ad7d3/scratchpad/throckmorton_tnd8233.pdf)を読み、図も確認しました。単独すきま・乱流・壁温比を変える実験を選んだ方向自体は適切です。[NASA の原報情報](https://ntrs.nasa.gov/citations/19760019344)

1. **Major — 「掃引角 0°＝二次元」は成立せず、入口 BL の再構成も不足している**

   **根拠:** [plan:52](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:52)、[plan:137](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:137)。

   原報を確認すると、厚さを示す Fig.4 は **TH76 自身の図**です。TH75 の Fig.4 は模型写真で、速度分布は Fig.11–12、積分厚さは Fig.13 です。

   さらに [TH75](/home/sano/work/forge/papers/gap_heating/Throckmorton_1975_NASA-TN-D-7939_RSI_tile_array_gaps_M10.3.pdf) 本文 pp.8–11 は、横方向の BL 分布変化、壁圧と自由流静圧の差、トンネル圧縮領域の影響を明記しています。掃引角をゼロにしても、この流入場が二次元になるわけではありません。

   \(\delta^*,\theta\) の二つの積分値だけでは \(U,T,p,k,\omega\) の分布は決まりません。厚さの ±20% は、熱的履歴や乱流分布の不確かさを代替しません。また [`boundaryCond_d.cu:814`](/home/sano/work/forge/solver_density_cuda/cuda_forge/boundaryCond_d.cu:814) では、速度入口の亜音速部分は CPG なら特性関係、TP なら圧力外挿で状態を作り、CSV の全量固定にはなりません。

   **対案:** T2 をまず「中央断面の二次元近似」と定義する。TH75 の分布・熱的再構成と TH76 Fig.4 を使い、入口から開口直前までの \(U,T,p,k,\omega\)、加熱区間長、実効 \(\delta^*,\theta\) を固定する。**smooth 模型の圧力・熱伝達と開口直前 BL を再現できることを、T2 本体への投入条件にする。** この条件を満たせなければ、二次元 T2 は探索に降格する。

2. **Major — T2 の \(h_{smooth}\) は依然として定義が違い、温度条件も未確定**

   **根拠:** [plan:76](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:76)、[plan:122](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:122)、[plan:140](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:140)。

   TH76 本文 p.6 と Appendix p.10 の定義は、

   \[
   h=\frac{q}{T_{aw}-T_w},\qquad
   T_{aw}=T_\infty+0.89(T_t-T_\infty),\qquad
   h_{fp}=C(Re'_\infty)^{0.69}
   \]

   です。分母は **等温基準模型の測定から作った相関値**で、温度比ごとに加熱 smooth CFD から作り直す値ではありません。係数 \(C\) は Appendix 図から確定する必要があります。指数だけでも、Re を \(1.5\to7.8\times10^6/\mathrm m\) とすると分母は **3.119 倍**になります。G13 の向きはこの正規化に依存します。

   温度比だけでも条件は閉じません。\(T_t,T_\infty,T_{gap}\)、上流加熱板の長さ・温度、非加熱の下流板の扱いが必要です。Fig.5(b) の実際の温度比は **1.00／1.18／1.34／1.52**で、公称 1.0／1.2／1.4／1.6 とも異なります。

   **対案:** T2 専用の定義表を作り、各測定系列の実条件・\(T_{aw}\)・共通の実験由来 \(h_{fp}\) を固定する。CFD も原報と同じ演算で \(h/h_{fp}\) を作る。smooth CFD は分母診断に限定し、実験値との比較には使い回さない。

3. **Major — 壁温不連続を共有 node に直接置くと、境界条件の順序依存と熱流束特異性が生じる**

   **根拠:** [plan:153](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:153)。[`nodeWallDirichlet_d.cu:74`](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:74) は壁 node の温度を `Ts` に代入し、[154行](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:154) は各 `bcond` について順に実行します。異なる `Ts` の壁が角 node を共有すれば、後の境界条件が上書きします。

   これは実装上の問題だけではありません。異温度の壁が接する理想的な鋭角では、局所的な伝導解に
   \[
   |\partial T/\partial n|\sim |\Delta T|/(\alpha r)
   \]
   が現れます。壁温ジャンプを正則化せず角まで積分した個別壁熱量には、対数的な格子依存が生じ得ます。

   TH76 の Fig.1 は、加熱板と薄肉すきまの間に断熱・冷却構造を持っています。単純な二つの等温面への置換を無条件には認められません。

   **対案:** 加熱板終端、冷たい接続部分、リップ形状、温度遷移長を実験形状として定義する。共有 node に矛盾する温度を与えず、温度遷移長を固定して格子細分化する。寸法が不明なら、その感度を入力不確かさに含める。

4. **Major — T2 の合格を、case/49 のモデル誤差幅に移す論理はまだ成立しない**

   **根拠:** [plan:33](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:33)、[plan:210](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:210)、[G11–G13](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:262)。

   T1b の降格は正しいです。しかし \(D/W=20\) の一致は、深部熱伝達の検証範囲の一致を意味しません。TH76 の基本分布 Fig.5–7 は主に **\(z/W\lesssim6\)** を示し、本文は黒塗り記号を精度に疑義のある低加熱データとしています。底面 \(z/W\simeq20\) の熱量を拘束する根拠にはできません。

   また、加熱等温板と断熱板は、ある位置の温度比が似ていても熱的境界条件・履歴が異なります。外挿幅を列挙するだけでは、その外挿誤差の大きさは決まりません。

   **対案:** 受け渡す表を「**TH76 の測定範囲で観測した予測誤差**」と「**case/49 の適用不確かさ：未定量**」に明確に分ける。T2 合格で認めるのは測定された深さ・量・条件の予測実績までとし、case/49 の深部・底面・総熱量への誤差幅の移植は行わない。

5. **Major — M7 は部分採用に留まり、合否帯が依然として数値になっていない**

   **根拠:** [plan:199](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:199)、[plan:203](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:203)、[G4–G6](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:255)。

   ±25% の撤去は改善ですが、「施設間・配列間・相関のばらつき」を工学帯にするだけでは、その数値も条件差の除去方法も不明です。異なる幾何・流入条件による差を許容誤差に取り込むと、今回排除した相関混用が戻ります。

   また、W70 の ±0.68 kW/m² は報告された精度限界であり、標準偏差と確認した値ではありません。現在の線形加算式を統計的不確かさ伝播と呼ぶのは不適切です。積分熱量には測定点間の相関と前壁理論補完が必要で、T2 にはさらに \(T_w,T_{aw},h_{fp}\) の不確かさが必要です。TH76 の伝導補正に関する「10%未満」も総合測定不確かさそのものではありません。

   **対案:** CFD 結果を見る前に、点別・積分量別の判定表を確定する。精度限界なら保守的な区間として、統計量なら共分散を使って伝播する。TH76 の黒塗り点は精度検証から除外する。根拠のある工学帯を作れない量は、合否ではなく誤差報告に留める。

6. **Major — 「離散エネルギー収支」の取得方法がなく、ソルバ無改修というスコープと接続していない**

   **根拠:** [plan:179](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:179)、[plan:280](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:280)。

   [`nodeWallDirichlet_d.cu:85`](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:85) は壁のエネルギー残差をゼロ化します。さらに [`variables.hpp:223`](/home/sano/work/forge/solver_density_cuda/variables.hpp:223) の出力登録には、エネルギー数値流束や拘束前の `res_roe` がありません。[`output.cpp:55`](/home/sano/work/forge/solver_density_cuda/output/output.cpp:55) により、未登録量は `extraFields` に書いても出ません。

   保存場から開口の \(\rho u_nh_0\)、熱伝導、粘性仕事を積分することはできますが、それは **SLAU の数値散逸・再構成・node の拘束操作を含む離散収支そのものではありません**。

   **対案:** 非拘束 CV 群と、その境界を横切る実際の数値流束を対象に収支を定義する。取得に必要な診断出力を別 plan の先行依存として明示する。保存場からの物理収支は別の診断として残し、G8 の代用にしない。許容値も片道の大きな流入量ではなく、検証する壁熱量の誤差予算に対して設定する。

7. **Major — 量別格子検証が T1 にしかなく、受け渡しの核心である T2 が未検証になる**

   **根拠:** [残作業表:237](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:237)、[G7:258](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:258)、[G10–G13](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:261)。

   G7 は明示的に T1 の二幅だけです。T2 の \(y_1^+\le1\) は壁解像の必要条件であって、開口せん断層、温度不連続、前後壁差、反転深さの格子誤差を保証しません。GCI を計算するだけでも、合格にはなりません。

   **対案:** T2 の代表条件、少なくとも温度比の両端で三段階格子を追加する。実測点の \(h/h_{fp}\)、区間熱量、反転深さごとに数値誤差上限を決める。例えば許容比較誤差の一定割合を数値誤差予算に割り当て、その割合を事前固定する。非単調収束なら GCI の形式適用をせず、未解決として扱う。

8. **Major — 窓二倍化だけでは遅い緩和を排除できず、IC 判定の尺度も量ごとに未定義**

   **根拠:** [plan:164](/home/sano/work/forge/plans/active/case-hypersonic-gap-heating-validation.md:164)、[`check_quasisteady.py:278`](/home/sano/work/forge/solver_density_cuda/tools/check_quasisteady.py:278)。

   判定器を変更せず、合成系列
   \[
   R(s)=0.01(2-e^{-s/100000})
   \]
   で窓延長を確認しました。

   | 点数 | VERDICT | 末尾平均 | 漸近値 |
   |---|---|---:|---:|
   | 12 | `STEADY` | 0.01000090 | 0.020000 |
   | 24 | `STEADY` | 0.01000185 | 0.020000 |

   代表値差は約 **0.0095%**ですが、真の極限との差は依然ほぼ100%です。これは IC 独立性を含む全ゲートへの反例ではありませんが、**窓二倍化が未緩和を排除するという根拠にはならない**ことを示します。

   また「相対2%かつ絶対 \(0.005q_{fp}\)」は、圧力比・W/m の区間熱量・T2 の熱伝達係数に共通適用できません。ゼロ近傍では相対許容が消失します。

   **対案:** 各量の単位・固定尺度・絶対床・比較式を定義する。窓延長は計算継続と評価窓拡大を区別し、単調系列では推定残余変化もゲートに入れる。代表 URANS は物理時間で熱・渦度拡散時間との関係を示す。IC 二本は同一の最終方程式で、比較時点にも十分な初期差が与えられたことを記録する。

   なお、case/48 の扱いを慎重にした修正は妥当です。再実行した  
   `case/48.flat_plate_cooled_m4/run_0011_Bplain_tw300_y3/` の**本段単独**判定は、
   `NOT CONVERGED (stalled/plateau)`、`rms_roe` 最終 \(1.42\times10^{-1}\)、低下0.1桁でした。同一設定区間全体の判定とは区別します。

**推奨は、T1/T2 の構成を維持しつつ、T2 に「試験再構成の成立ゲート」を先行追加してから実装することです。**

実装前の優先順は、①T2 の実条件・共通分母・二次元近似の成立条件、②壁温接続部と離散収支の取得方法、③量別の不確かさ・格子・定常性の数値基準、④case/49 への受け渡しを検証実績に限定する文言、です。これらが閉じるまでは、現稿のまま本計算へ進むことを勧めません。

指摘数: Critical 0 / Major 8 / Minor 0
