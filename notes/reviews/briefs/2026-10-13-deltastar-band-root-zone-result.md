# 諮問ブリーフ: δ* 抽出の修正ほか 4 つのコード変更の result 段の解釈 (plan tooling-nozzle-deltastar-band-root-zone、main へ入れる前)

AGENTS.md 条件 7 (result 段の解釈を確定する前、`--stage result` に回す前)。判断役のスイッチが codex。
ユーザ決定 (2026-10-13): M5 デモ (V5 の出口 M) は判定保留のまま区切る (plan §9)。今回 main へ入れるのは**コード変更**で、その正しさと後方互換を判定する。
V5 (冷却壁の NS の遅いモードで出口 M が定常にならない) は本 result の合否の対象外とし、plan §5.1 #11 に残す。

## 変更と検証 (commit と試験)

| 変更 | commit | 検証 (事前登録の節) | 結果 |
| --- | --- | --- | --- |
| `band_local_deficit` の等価半径の探索を壁〜y_b に限る、非単調は `nonmono_zone` hard、符号分岐の前に探索域の q_ref の有限・正を検査 (`nonpos_ref` hard) | ae34f09e | §6 V1・V2 | V1: 旧の δ_r = r_w → 新は真値と相対 1.3e-7、codex の負例 (δ_r −4.98e13 soft) → hard NaN、F が単調な既存の合成断面 13 個はビット同一。**V2: case/45 run_0022 × run_0001 (欠陥に当たる呼び出し 0 件) で `delta_r_equiv.csv`・`delta_r_next.csv` がバイト同一**。case/70 の抽出は測定用の差し替え (`ab_root_zone.py` の zone の腕) とバイト同一 |
| `solve_rt_throat --prev`: 使用域 [0, x_e] の hard/held で停止、表の端を前の設計の入口端・出口に両側 \|Δ\| ≤ ε_f32 で揃える、区間で正規化した座標 s で今の設計へ写す; `wt_card --ns-correction`・`deltastar_loop --delta-r-csv` | ae34f09e | §6 V4 (a)〜(h)・カード (i)〜(iii) | 全合格。同じ設計への写しは恒等、`sizing_ab.py` の腕 B の r_t がビット同一 (0.02507572544991335)。NS 補正なしのカードは時刻を除き同一 |
| 設計壁の単調性の検査を、4000 点の隣接差分 (< −1e-9) から連続関数の複合ゲート (最小の壁角 ≥ −1e-4°・最大の累積減少 ≤ 1e-7 r_t) に | 35ab6821 | §5.1 #4b の事前登録の A/B | 実壁 33 本 (旧で落ちる 14 本を含む) が全部合格、最悪 min θ −1.18e-5°・D 2.9e-8 r_t; 負例 4 つは予定の理由だけで拒否 (旧検査は点の間の狭い収縮を見逃していた); 対照 5 つは合格; 既存の YAML 5 本 (case/70 の 4 本・case/45 生産) の設計はビット同一 |
| cplus の閉包で x_F が終端の壁点の丸め下に来たとき、その点を落とさず補間した交点を足す (deliverables plan #15) | de35cb0b | §5.1 #4c の事前登録の A/B | 落ちていた 4 設計 (case/70 c2b・c2・旧 ns2 +0.846 mm、case/45 生産 +5.353 mm) は移動 = r_t(x_F − x_e,旧) (差 0)、落ちていない 2 設計はビット同一。試験 3 本を「旧規則に戻すとビット同一」を確かめたうえで更新 |

- 回帰: design/tests の 43 本 (結果は呼び出し側が添付する; 実行中)。
- M5 デモ: 作り直した設計は寸法ゲート PASS (拡大部長 −0.478 mm)。最終 NS は 60k で収束 PASS・出口コア M 5.0022、延長で 4.9991 (cfl 1、通算 240k)・4.9979 (cfl 5) と下がり続ける (遅いモード)。
- M6 のセッションには 3 つの修正を伝え、取り込む前提として先方の plan に記録済み。

## 当方の解釈 (確定前)

4 つのコード変更は、それぞれ事前登録の判定に合格し、欠陥の無い入力では出力がビット同一 (抽出・単調性・出口の終点)。main へ入れてよい。
変わる振る舞い: (1) δ_r = r_w・巨大な負値を返していた断面が hard (NaN) になる、(2) 出口直前の微小な非単調で設計が落ちなくなり、点の間の収縮は落ちる、
(3) 終点が落ちていた設計の壁・メッシュ・STEP が 1 点間隔長くなる (case/45 生産 +5.35 mm)、(4) `solve_rt_throat --prev` は前の run の `wall_repr.json` が必須。

## 問い

**上の結果のうち、「4 つのコード変更を main へ入れてよい」という主張を支えていない (あるいは過大に言っている) 点を 1 つ挙げ、PR の前に足すべき確認を 1 つ推奨してほしい。**
候補 (縛られなくてよい): (a) (3) は case/45 の生産の壁を変えるので、その NS の結果への影響を見ていない、(b) `nonpos_ref` が case/45 の手元の対で各 ~6000 呼び出し (梯子の段だけ) あり、
生産の抽出 (AWS) で最終値に出ると NS 後の寸法が止まる、(c) V2 は 1 対だけ、(d) 単調性のゲートの閾値の根拠は実用許容で、雑音の上限の測定ではない。

## 禁止事項 (厳守)

- ファイルを変更しない (read-only サンドボックスで動いている)。
- **`*.log`, `residual_history*.csv`, `res_*.h5`, `*.vtu`, `*.msh`, `*.cgns`, `*.stl`, `plans/README.md`, `delta_r_equiv_diag.json` を読まない** (巨大)。
- 下に列挙した `sed -n 'A,Bp' <file>` 以外のファイル読みをしない。grep は可。
- 両論併記で逃げず、**推奨は 1 つに絞る**。根拠は `ファイル:行` か提示した数値で示す。

## 読んでよいもの

- `sed -n '1,260p' plans/active/tooling-nozzle-deltastar-band-root-zone.md`
- `sed -n '260,310p;580,600p' design/forge_design/metrics/deltastar.py`
- `sed -n '150,335p;500,545p;630,700p' design/forge_design/feedback/deltastar_loop.py`
- `sed -n '40,140p;225,245p' design/forge_design/geometry/wall_axismach.py`
- `sed -n '570,650p' design/forge_design/geometry/moc_inverse.py`
- `cat case/70.m5_air_wt_demo/_exit_cut_ab/exit_cut_ab.json | head -120`
- `grep -n "adopt_B\|verdict" case/70.m5_air_wt_demo/_monotone_gate_ab/monotone_gate_ab.json | head`
