# codex 諮問 (diagnose): bell-m-corrected-wall-result

- **brief**: [`notes/reviews/briefs/2026-10-13-bell-m-corrected-wall-result.md`](../../notes/reviews/briefs/2026-10-13-bell-m-corrected-wall-result.md)
- **plan**: [`plans/active/tooling-nozzle-bell-two-mode.md`](../../plans/active/tooling-nozzle-bell-two-mode.md)
- **date**: 2026-10-13
- **commit**: `09f644a0` (feature/thruster-nozzle-design)
- **codex**: effort `high`, 4.4 min, rc=0
- **結論**: **原因の記入を保留し、既存場の CPU 後処理で、壁圧差を「共通座標での圧力差」と「そこから各壁までの圧力増分の差」に分解する A/B を先に行う。**

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

| 判断対象 | 採否 | 理由 |
|---|---|---|
| 凍結した判定で「補正は効く」 | **採用** | 記録された VERDICT はそのとおり。壁圧 RMS の減少率は V-ext 65.7 %、V-zero 65.3 %。「今回の2通りとも登録線を満たす」という意味に限定する。[判定結果:118](/home/sano/work/forge-thruster/case/68.bell_two_mode/rao_m/corrected_wall_judgement.json:118) |
| 残る負の壁圧差の原因を確定する | **要再検証** | 壁圧差だけでは補正過大と、壁近傍の圧力分布の違いを分離できない。 |
| x < 1.17 の既定を V-ext にする | **今回の L100 条件の暫定既定として採用** | V-zero より `pw_t` が小さく、補正区間の最大壁圧誤差も 0.873 % 対 3.427 %。ただし、推定厚さの正しさや他条件への一般化を証明した結果ではない。[判定結果:56](/home/sano/work/forge-thruster/case/68.bell_two_mode/rao_m/corrected_wall_judgement.json:56) |

結論: **原因の記入を保留し、既存場の CPU 後処理で、壁圧差を「共通座標での圧力差」と「そこから各壁までの圧力増分の差」に分解する A/B を先に行う。**

第 1 仮説: **残る壁圧差の主成分は、共通の外縁位置から各壁までの圧力増分の違いである。** 確度: **中。ただし寄与の大きさは未確認。**

根拠: 現行の `pw` は、同じ x にある**異なる半径の壁面**の圧力を比較する。一方、`edM` は NS と Euler の**共通座標**で比較する。[judge_corrected_wall.py:98](/home/sano/work/forge-thruster/case/68.bell_two_mode/rao_m/judge_corrected_wall.py:98)、[同:109](/home/sano/work/forge-thruster/case/68.bell_two_mode/rao_m/judge_corrected_wall.py:109)。V-ext の `edM` は 0.00178 まで小さくなったが、`pw` は 0.739 % 残る。これは本仮説を調べる動機になるが、Mach の近さだけでは圧力の近さを保証しない。

**[Major] 最も重い誤読は、「δ_r の補正が正しければ壁圧差は格子差程度になるはず、したがって負の差は補正過大」という推論である。** 流量欠損で較正した δ_r と、壁圧一致との間を、今回の資料は検証していない。対案は、係数を動かす前に次の圧力差分解を行うこと。

反証条件: 壁までの圧力増分の差を除いても、共通座標での圧力差がほぼ残ること。下記では「70 %以上除去できる」を第1仮説の具体的な判別線とする。

第 2 仮説: **補正前の縁条件による一発補正が、補正後の必要厚さとずれている。** 確度: 低〜中。`CONTUR` は `run_0018` に固定されているが、補正後の厚さとのずれは未測定。[run_corrected_wall.py:48](/home/sano/work/forge-thruster/case/68.bell_two_mode/rao_m/run_corrected_wall.py:48)

第 3 仮説: **固定した物理スロートに残る排除効果が、一部を説明する。** 確度: 低〜中。スロートの補正量とその1・2階微分はゼロ。[offset_vext.json:28](/home/sano/work/forge-thruster/case/68.bell_two_mode/rao_m/rao_corr_parts/offset_vext.json:28)。ただし、提示された準1次元概算の −0.23〜−0.25 %を採っても、全差の説明には足りない。

判別 A/B: **変えるのは、後処理で壁近傍の圧力増分差を差し引くかどうかの1点。CFD は0ステップ。**

対象は `case/68.bell_two_mode/run_0026_tp_m_ns_l100_corr_vext/` と `case/68.bell_two_mode/run_0024_tp_m_euler_l100_inv/`。既存 `edM` と同じ一定 x 補間・凍結した `r_s = r_wall,NS − y_b` を使い、両場の圧力を抽出する。壁を w、共通座標を s として、

```text
e_w    = (p_Nw − p_Ew) / p_Ew
e_s    = (p_Ns − p_Es) / p_Ew
e_near = [(p_Nw − p_Ns) − (p_Ew − p_Es)] / p_Ew

e_w = e_s + e_near
```

比較するのは現行 `e_w` と、圧力増分差を除いた `e_w − e_near`。**分母、x = 1.17〜L、Δx 重み、RMS は現行壁圧指標のまま**にする。NS 側の壁−外縁差だけを引くと、Euler 側にもある圧力変化を原因に数えてしまうため、必ず両者の差を取る。

- **A:** 差し引き後の RMS が元の30 %以下、すなわち現集計値を基準に約 **0.222 %以下**で、`e_near` の平均が負。→ 壁近傍の圧力増分差が大きな寄与を持つという第1仮説を支持する。ただし、それを曲率など特定の物理機構に帰属する段階は残る。
- **B:** 満たさない。→ 「壁近傍の圧力増分差で70 %以上説明できる」を棄却する。共通座標でも残る不一致として扱い、第2・第3仮説のどちらかをこの結果だけで確定しない。

同じ比較を保存済み本段時系列と既存の Euler 細格子でも確認する。圧力指標の時系列には `check_quasisteady.py` を適用し、未定常なら物理原因の判定を留保する。この70 %は**新しい診断用の線**であり、凍結した補正効果の合否は変更しない。

やらない方がよいこと: **壁圧差の符号だけで `kδ` を下げること。** また、`p_NS` を流量比で一律に正規化して差を消すこと。どちらも、まだ分離していない寄与を較正係数や正規化に吸収する。

呼び出し側の前提への異議:

- 「採り方に依らない」は**合否についてのみ**採用する。`axM` の系列差 0.0309 は登録許容幅 0.0171 より大きく、`pw_t` の差も0.459ポイントある。推定区間の影響が小さいという意味には使えない。[判定結果:101](/home/sano/work/forge-thruster/case/68.bell_two_mode/rao_m/corrected_wall_judgement.json:101)
- **[Minor]「一様な −0.7 %」は範囲が広すぎる。** V-ext の区間平均は x = 1.17〜2 で −0.432 %、8〜13 で −0.797 %。対案は「下流 x ≥ 8 で約 −0.8 %、上流側には変化がある」と記録すること。[内訳:190](/home/sano/work/forge-thruster/case/68.bell_two_mode/rao_m/corrected_wall_breakdown.json:190)、[同:220](/home/sano/work/forge-thruster/case/68.bell_two_mode/rao_m/corrected_wall_breakdown.json:220)

不足情報: 共通座標の圧力、壁圧・Mach 指標自身の準定常 VERDICT、NS 側の格子感度。現行判定は最終場を読むため、推力系の `STEADY` だけで壁圧の定常性まで確認したとは扱えない。[judge_corrected_wall.py:54](/home/sano/work/forge-thruster/case/68.bell_two_mode/rao_m/judge_corrected_wall.py:54)

指定された読み取り範囲でレビューし、禁止された実場ファイルは未読、A/B は未実施、ファイル変更なし。**plan 未反映**。呼び出し側で `plans/active/tooling-nozzle-bell-two-mode.md` §5.1 #7 に採否と診断線を記録すること。
