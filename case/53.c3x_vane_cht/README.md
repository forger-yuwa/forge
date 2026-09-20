# case/53 C3X タービン翼 (CHT の公知データ検証 V5 — 亜音速出口)

**壁温が結果として測られている**公知試験に対して forge の共役熱伝達を当てる。
主対象は超音速出口の Mark II ([case/54](../54.markii_vane_cht/)) で、**本 case はその足場** (亜音速出口)。
計画は [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md) §4.9 / §6 V5、
選定理由は [`notes/investigations/cht-validation-case-survey.md`](../../notes/investigations/cht-validation-case-survey.md)。

## 一次資料

**NASA CR-168015** (Hylton, Mihelc, Turner, Nealy, York, 1983; NTRS **19830020105**)。
PDF は `papers/cht/NASA-CR-168015_Hylton_1983.pdf` (**git 追跡外**。再取得は
`curl -L -o papers/cht/NASA-CR-168015_Hylton_1983.pdf https://ntrs.nasa.gov/api/citations/19830020105/downloads/19830020105.pdf`)。

抽出は [`tools/extract_vane_data.py`](tools/extract_vane_data.py)。**スキャン PDF の OCR なので検証してから使う**:

- 表 II / III は **cm と inch の両方**を載せているので、全点に `cm == 2.54 × in` (相対 0.5 %) を課す。
  通らない点は**捨てて**、ページ画像を目視して読んだ値 (`REPAIRS`、各行を同じ 2.54 則で検算済み) で埋める。
  点が 1 つでも欠けたら**書き出さずに落ちる**。
- 現状: C3X **73/78 が OCR 検証を通過、5 点を目視補修**。Mark II **54/60 通過、6 点を補修**。
  形状の自己検査: 連続点の最大間隔は C3X 0.618 cm / Mark II 1.626 cm (後者は**鈍頭前縁** $R_{LE}$=1.280 cm の
  円弧区間で、点列に無いのが正しい)。`ref/vane_shapes.png` で目視確認。

### 報告自身の不整合 (2 件、いずれも記録して回避)

1. **表 VIII / IX の SI 圧力列が psia 列と 51.7 倍ずれている** (例: C3X 4411 は "6177 (46.34 psia)" だが
   46.34 psia = 319.5 kPa)。**psia 側を正**とする (文献が引用する ~3.2 atm と一致する)。
2. **C3X 表 III の点 29 は inch 値が誤植**: "0.4115 cm (0.0162 in)" だが 0.0162 in = 0.0411 cm。
   前後の点間隔 (~0.6 cm) から **cm 値 0.4115 が幾何的に正しい** (0.0411 だと点 30 と 0.047 cm しか離れない)。

## 抽出済みデータ

| ファイル | 中身 |
| --- | --- |
| `ref/vane_c3x.csv` | C3X 翼型 78 点 (x, y [cm])。$R_{LE}$=1.168 cm, $R_{TE}$=0.173 cm |
| [`../54.markii_vane_cht/ref/vane_markii.csv`](../54.markii_vane_cht/ref/vane_markii.csv) | Mark II 翼型 60 点。$R_{LE}$=1.280 cm, $R_{TE}$=0 (blunt) |
| `ref/test_conditions.csv` | 表 VIII / IX の主要 run |
| `ref/vane_shapes.png` | 形状の目視確認図 |

主要 run (表 VIII / IX、ページ画像から読取):

| vane | code | run | $P_{T1}$ | $T_{T1}$ | $M_1$ | $Re_1$ | $M_2$ | $Re_2$ | Tu | $T_w/T_g$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **C3X** | 4411 | **108** | 46.34 psia (319.5 kPa) | 786 K | 0.17 | 0.52e6 | **0.90** | 1.99e6 | 6.5 % | 0.73 |
| **Mark II** | 5411 | **42** | 48.89 psia (337.1 kPa) | 788 K | 0.19 | 0.56e6 | **1.04** | 2.01e6 | 6.5 % | 0.68 |

4411 (亜音速出口) と 5411 (超音速出口) は入口条件がほぼ同じなので、**出口マッハだけを変えた対**として使える。

## まだ取っていないもの (次の作業)

1. **材料熱伝導率 $k_s(T)$** (ASTM 310 ステンレス) と**冷却孔の配置・径** (図 6/7 = 有限要素モデル図)。
2. **冷却孔ごとの熱伝達係数と冷却剤温度** (報告 p.21: 入口助走補正つき相関 + 入口/出口実測からの推定)。
3. **測定壁温・熱伝達率の分布** (各 run)。報告のデータ整理は**薄肉ではなく 2 次元伝導 FEM**
   ("a 2-D plane of the test vane as a fluxmeter"、内外の境界条件を測って定常熱伝導を解く) なので、
   forge 側も `fem2d` で合わせる。
4. カスケード幾何 (表 IV: ピッチ・翼弦・取付角) → メッシュ生成。

## 計算 run 一覧

| `run_*` | 目的・主要設定差分 | 主要結果・成果物 | 状態 |
| --- | --- | --- | --- |
| (まだ無し) | — | 一次資料の抽出のみ完了 | — |
