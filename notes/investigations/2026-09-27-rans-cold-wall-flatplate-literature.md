# 冷壁の極超音速乱流平板で RANS (SST) は熱流束を過大か過小か — 文献精査 (2026-09-27)

目的: D-7275 試験 26 (M 6.64、Tw/Tt ≈ 0.16、Tw/Taw ≈ 0.17) で forge SST の St が実験の 1.31 倍・Eckert の 0.99〜1.06 倍になった
([`case/60.flatplate_d7275_m7/README.md`](../../case/60.flatplate_d7275_m7/README.md)「Eckert 参照温度法との比較」) ことを、
公開文献の RANS/DNS/実験の比較と照らす。ユーザ指摘: AIAA 2025-1689 は「RANS は DNS より過小」と言っている。

調査はサブエージェント (Web) で行い、主要な引用 3 件 (NTRS 2020 の 14.0 %/15.2 %・21.8 %、Rumsey の 40−100 %・"too low") は原文で再確認した。

## 結論

- **比較の基準と対象量で符号が変わる。** AIAA 2025-1689 (Vaughan/Roy/Jordan、改訂版は VT 修論 2025) の「過小」は主に TKE・速度・外層温度のこと。
  壁熱流束 C_h は **Re_τ を揃えた比較**で「低 Re で過小、高 Re で過大、概ね ±10 %」(SST・SA とも圧縮性補正なし)。
- **x を揃えた冷壁の比較では素の SST は過大。** M11 冷壁 (Tw/Tr 0.20、CUBRC):
  SST は DNS より q_w +15.2 %・τ_w +14.0 % (x = 1.35 m)、DNS は実験より +21.8 % (x = 1.0 m) → SST は実験比で最大 +35 % (Gnoffo らの結果として PRF 2022 が引用)。
  HIFiRE-1 (M 6.6〜7.2、冷壁) でも 1〜2 方程式モデルは熱流束を過大。原因は τ_w の Re 依存が緩い (減り方が遅い) ことで、2C_h/C_f は RANS ≈ 1.19 で DNS と近い。
- Rumsey (NASA/TM-2009-215705): 素の SST の C_f は冷壁 M 10 で相関 (VD-II/Spalding–Chi) より 40〜100 % 過大。Sarkar/Zeman 型の散逸補正は C_f・熱伝達を下げすぎる傾向。
- Zeman 補正 (PRF 2022、M11) は q_w を 12.5〜13 % 下げる。
- DNS と相関: 冷壁 (Tw/Tr ≲ 0.3) で VD-II は DNS の C_f を 10〜20 % 過大 (JFM 2022)。**DNS を Eckert と直接比べた論文は見つからなかった。**

## forge の結果への含意

- forge の St ≈ Eckert ≈ 1.31 × 実験は、「冷壁・x 基準で素の SST も相関も実験より高い」という文献の構図と同じ向き・同じ大きさ (M11 で実験比 +30〜35 %)。
- **forge の `dilatationCorrection 2` は Sarkar/Zeman 型の散逸補正ではない** ([`methods/turbulence/implementation.md`](../../methods/turbulence/implementation.md) §5.3:
  生産項を偏差ひずみ + −⅔ρk∇·u にするだけ)。文献の「補正なし SST」に当たる。サブエージェントの「補正が入った状態で 1.31 倍」という推定は誤りとして採らない。
- 文献の内訳 (RANS > DNS が約 15 %、DNS > 実験が約 20 %) に照らすと、1.31 倍のうち「SST の過大」と「実験が DNS・相関より低い」の配分は D-7275 では未切り分け。
- 試験 26 の Re_θ・Re_τ (2026-09-27 計算、`run_0002_t26_ext/res_60000.h5` の壁法線分布、δ99 は u = 0.99 u_e、u_e は y 0.06〜0.1 m の中央値、Re_θ は外縁物性、Re_τ は壁物性):

  | x [m] | δ99 [mm] | θ [mm] | Re_θ | Re_τ | Re_x (外縁) | T_e [K] |
  |---|---|---|---|---|---|---|
  | 1.146 | 17.4 | 1.04 | 4137 | 1532 | 4.57e6 | 254 |
  | 1.879 (位置 II) | 27.6 | 1.54 | 6453 | 2346 | 7.86e6 | 247 |
  | 2.546 | 36.9 | 1.96 | 8496 | 3060 | 1.10e7 | 243 |

  比較域の Re_τ 1500〜3100 は DNS データベースの上限 (Re_τ ≤ 1244、M11 は 774/1172) を超える → VT (Re_τ 基準) の「高 Re で過大」側、
  x 基準の M11 比較 (SST が DNS より +15 %) とも同じ側。T_e が T∞ 228 K より高いのは前縁の弱い衝撃波 (P ≈ 1.08 p∞) の背後のため。

## 出典 (読んだもの)

- VT 修論 (Vaughan 2025、AIAA 2025-1689 を改訂・置換と明記、全文): https://vtechworks.lib.vt.edu/server/api/core/bitstreams/e9b98180-2cd2-4e23-a236-2ae6d804c7cc/content
- AIAA 2025-1689 (403、抄録スニペットのみ): https://arc.aiaa.org/doi/10.2514/6.2025-1689
- Huang, Nicholson, Duan, Choudhari, Bowersox, AIAA 2020-0571 (NTRS 20200002806、全文): https://ntrs.nasa.gov/api/citations/20200002806/downloads/20200002806.pdf
- Barone, Nicholson, Duan, Phys. Rev. Fluids 7, 084604 (2022、受理稿全文): https://link.aps.org/accepted/10.1103/PhysRevFluids.7.084604
- Huang, Bretzke, Duan, Fluids 4, 37 (2019、全文): https://doi.org/10.3390/fluids4010037
- Huang, Duan, Choudhari, JFM 937 A3 (2022、OA 全文): https://doi.org/10.1017/jfm.2022.80
- Rumsey, NASA/TM-2009-215705 (全文): https://ntrs.nasa.gov/api/citations/20090015399/downloads/20090015399.pdf
