# 凝縮物性の蒸気比熱 c_p,v をどこから取るか

## メタ

- **area**: `condensation`
- **status**: `in_progress`
- **related_docs**:
  - [`methods/condensation.md`](../../methods/condensation.md) §8b 物性の外挿規約
- **related_plans**: [`plans/accepted/condensation-air.md`](../accepted/condensation-air.md) (70 K 未満の Kirchhoff 外挿の出典)、
  [`plans/active/condensation-followups.md`](condensation-followups.md) (F-cf9 として起票、本計画へ分離)
- **created**: `2026-09-14`
- **owner**: `Claude`

## 1. 目的

N2 潜熱の 70 K 未満の外挿 $L(T)=L(70)+(c_{p,v}-c_l)(T-70)$ と、50 K 未満の飽和圧 (同じ $L$ の積分) に入る
$c_{p,v}$ を、内蔵定数 1038.8 J/(kg·K) にするか config の `physProp.cp` にするかを決める。
ユーザ指摘: 「CPG なんだから config のその値を使うのではないのか。謎の使い分けはやめてほしい」。

## 2. スコープ

- **やる**: `c_p,v` の出どころを 1 つの規則で決める。実装・単体・case/34 (node/cell) と case/16 (node/cell) の回帰。
- **やる**: そのとき壊してはいけない不変量 (Kantrowitz の $\gamma_v=c_p/c_v$、$c_p-c_v=R$) を明示し試験にする。
- **やらない**: CPG carrier の二相 EOS が液相の顕熱に混合気の熱容量 ($c_v+R_w$) を当てている件 (§5.1 #1)。
- **やらない**: 液相比熱 `condN2LiquidCp` の値そのもの (別の感度課題)。

## 3. 関連 docs と前提

- $L$ の 70 K 未満の外挿と $p_{sat}$ の 50 K 未満の外挿は [`plans/accepted/condensation-air.md`](../accepted/condensation-air.md) §4.2 で導入。
  そこでは $c_{p,v}=1038.8$ (N2 蒸気) 固定だった。
- `physProp.cp` は `thermalMethod: 0` (CPG) のときの**気相の**定圧比熱。純 N2 の run では 1038.8、
  空気キャリアの run では混合気の 1008.7。TP (`thermalMethod: 2`) では使われない。

## 4. 設計方針

**規則は 1 つ**: `thermalMethod: 0` (CPG) なら `physProp.cp`、それ以外 (TP) は内蔵の種固有値。

CPG は「気相の比熱はこの 1 つの定数」というモデルなので、物性相関だけ別の定数を持たせると同じ run の中に
2 つの $c_p$ が並ぶ。TP には「気相の $c_p$ 定数」が存在しないので内蔵値を使うほかない。

厳密には Kirchhoff の関係に入るべきは**凝縮する窒素蒸気**の $c_p$ で、空気キャリアの 1008.7 はその値ではない。
ただし CPG 二相 EOS は液相の顕熱にも混合気の熱容量 ($c_v+R_w=1014.7$) を当てており、**凝縮種固有の $c_p$ を
そもそも表現できない**。config の値に統一するのはその枠内で自己整合を取る選択で、差 30 J/(kg·K) は
液相比熱 $c_l$ 自体の感度幅 ($\pm500$ J/(kg·K)) の内側にある。

### 4.1 触ってはいけないもの (codex result M1)

$c_p$ を config 値で**上書きしてはならない**。`CondSpeciesProps` の `cp`/`cv`/`R` は種固有の整合した組で、
Kantrowitz の非等温補正が $\gamma_v=c_p/c_v$ を使う。空気の 1008.7 を `cp` に入れると
$c_p-c_v=266.7\neq R=296.8$、$\gamma_v=1.359$ となり、補正係数 $2(\gamma_v-1)/(\gamma_v+1)$ が
0.3333 → 0.3047 に動く。**Kirchhoff 専用の別フィールド `kirchhoffCpv` を持つ**。

### 4.2 途中で混入した不具合 (本計画で解消)

`580dd8be` で `CondPropOpts::gasCp` を入れたとき、TP にも 1038.8 を渡して `CondSpeciesProps::cp` を
上書きしていた。H2O の $c_p$ が 1855 → 1038.8 になり $\gamma_v=0.745<1$、Kantrowitz mode 1 の
$\theta=2(\gamma-1)/(\gamma+1)\,b(b-0.5)$ が負になって 0 にクリップされ、**非等温補正が無効化**されていた。
`580dd8be`–`e0ffca34` の間の H2O 凝縮 run はこの影響を受ける (対象は本セッションの回帰 run のみ)。
§4.1 の分離で解消する。

## 5. 実装ステップ

1. `CondSpeciesProps` に `kirchhoffCpv` を追加し、`n2_latent_ex` / `n2_psat_ex` はこれを使う。`cp`/`cv`/`R` は不変。
2. `cond_prop_opts`: `o.gasCp = (thermalMethod == 0 && cfg.cp > 0) ? cfg.cp : 0`。0 は内蔵値。
3. 単体試験に $\gamma_v$ と $c_p-c_v=R$ の不変量、H2O の $c_p$ 不変と $\theta>0$ を追加。
4. `methods/condensation.md` §8b の該当節を書き直す。

### 5.1 残作業 (優先順)

| # | 項目 | 内容 |
| --- | --- | --- |
| 1 | CPG carrier EOS の液相顕熱 | 二相 EOS が液相にも混合気の $c_v+R_w$ (空気 1014.7) を当てており、実効液比熱が意図した 2000 でなく 1976 になる (1.2 %)。定式化そのものの近似。直すには液相の顕熱を凝縮種の $c_p$ で持つ必要がある |
| 2 | cell のノイズ床の取り方 | case/16 cell は 3 反復では床を過小評価し、同じ比較が FAIL→PASS に変わった (§8.1)。凝縮 cell 系では反復数を増やすか床の取り方を規約化する |

## 6. 検証

- **単体**: 凝縮 8 本すべて PASS。新規に `kirchhoffCpv` の上書き、`cp`/`cv` 不変、$\gamma_v=1.4$ (N2) と
  1.331 (H2O)、$c_p-c_v=R$、H2O の Kantrowitz $\theta>0$、70 K 未満の傾きが $(c_{p,v}^{new}-c_{p,v}^{old})$ だけ動くこと。
- **case/34 (空気, node と cell)**: 内蔵 1038.8 と config 1008.7 の A/B。**onset (格子点判定) が動かないこと**が合否。
  場の差は変数別ノイズ床 + `diff_res.py --factor 2` で報告する (**PASS は期待しない。物性が変わるので場は動く**)。
- **case/34 (純 N2, node)**: config が 1038.8 なので**ノイズ床以内で不変**であること。
- **case/16 (H2O, node と cell)**: §4.2 の不具合が解消し、`483de03b` (gasCp 導入前) とノイズ床以内で一致すること。
- **判定基準**: (1) 単体全 PASS、(2) 空気 node/cell の onset が両判定基準とも不変、(3) 純 N2 と case/16 が床以内。

### 6.1 レビュー記録 (codex)

| 段階 | 日付 | 記録 | 判定 / 指摘 (C/M/m) | 対応 / 免除理由 |
| --- | --- | --- | --- | --- |
| plan 免除 | `2026-09-14` | — | — | F-cf9 として `condensation-followups.md` に起票済みの課題を、codex result レビューの指摘 (M3: 集約 plan に埋めるな) を受けて分離したもの。設計は下の result レビューで評価済み |
| result | `2026-09-14` | [`notes/reviews/2026-09-14-condensation-followups-result.md`](../../notes/reviews/2026-09-14-condensation-followups-result.md) | NO-GO, C0/M3/m2 | **全件採用** (§9) |

## 7. 影響範囲

- `cuda_forge/condensationProperties_d.cuh` (`CondSpeciesProps::kirchhoffCpv`, `condProps_make`, `cond_latent`, `cond_psat`)
- `cuda_forge/condensationTransport_d.cuh` (`cond_prop_opts`)
- `tests/unit/test_cond_air.cpp`、`methods/condensation.md`
- 場への影響: 空気キャリア (CPG carrier) のみ。純 N2・H2O・dry は不変。

## 8. 完了条件

- [x] `methods/condensation.md` を更新
- [x] 実装・検証完了 (§6、結果は §8.1)
- [x] codex レビューを §6.1 に記録し採否を反映
- [ ] `status` を `done` にし `accepted/` へ移動 (result 再レビュー後)

## 8.1 検証結果 (2026-09-14)

**単体**: 凝縮 8 本すべて PASS。$\gamma_v$ は N2 1.400 / H2O 1.331 で不変、$c_p-c_v=R$ を満たす。

**case/34 空気 (12000 step, `onset_analysis.py --series` は両者 STEADY)**:

| 量 | node 1038.8 | node 1008.7 | cell 1038.8 | cell 1008.7 |
| --- | --- | --- | --- | --- |
| onset ($\Delta p/p_{dry}>1\%$) | 2.20 in | **2.20 in** | 2.21 in | **2.21 in** |
| onset ($g>10^{-4}$) | 2.10 in | **2.10 in** | 2.13 in | **2.13 in** |
| $g_{exit}$ | 0.0583 | 0.0581 | 0.0590 | 0.0589 |
| 保存モーメント $\rho Q_0$ の最大差/基準場最大値 | — | 2.43 % | — | 2.92 % |

差の尺度は `max|Δfield|/max|base field|` (局所相対差でも出口値でもない)。`diff_res.py --factor 2` は
node/cell とも exit 1 (= 床超) で、これは**物性が変わったのだから当然**。合否は onset で見る。

**case/34 純 N2 (node)**: `diff_res.py --tolfile --factor 2` exit 0 (床以内)。config が 1038.8 で内蔵値と同じため。

**case/16 H2O (300 step, `483de03b` 直前のバイナリとの比較)**: node exit 0。
cell は 3 反復床では exit 1 だったが、**5 反復 (現行) + 3 反復 (旧) の合成床では exit 0**。
cell の atomicAdd 非決定性に対し 3 反復は床の推定が不十分だった (§5.1 #2)。

**run パス**: `case/34.arthur_n2_nozzle/run_0105_gascp_air_node_{base,new,new_r2,new_r3,dry}`,
`run_0106_gascp_n2_node_{base,new,new_r2,new_r3}`, `run_0107_gascp_air_cell_{base,new,new_r2,new_r3}`,
`case/16.nozzle_wys/run_0466_h2ocp_node2d_{pre,cur,cur_r2,cur_r3}`,
`run_0467_h2ocp_cell2d_{pre,pre_r2,pre_r3,cur,cur_r2..r5}`。
いずれも 12000 / 300 step の A/B で、`check_convergence.py` は NOT CONVERGED (plateau)。
case/34 は `onset_analysis.py --series` が STEADY なので、報告する onset・$g_{exit}$・壁圧比は定常値。

## 9. 変更ログ

- `2026-09-14` — 起票 (F-cf9 から分離)。codex result レビュー (C0/M3/m2) を全件採用:
  M1 (cp 上書きで $\gamma_v$ が壊れる) → `kirchhoffCpv` を分離し不変量を試験化。
  M2 (TP/H2O の検証漏れ・cell 回帰なし) → §4.2 の不具合を特定し case/16 node/cell を `483de03b` と比較、
  case/34 は空気 cell を追加。m4 (測定量の表記) → 「保存モーメント $\rho Q_0$ の最大差/基準場最大値」に統一。
  m5 (判定の適用範囲) → 未収束の準定常比較であること、onset は格子点判定であることを明記。
  M3 (集約 plan に埋めるな) → 本計画として分離。
