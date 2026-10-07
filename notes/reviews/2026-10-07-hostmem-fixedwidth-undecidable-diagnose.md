# codex 諮問 (diagnose): hostmem-fixedwidth-undecidable

- **brief**: [`notes/reviews/briefs/2026-10-07-hostmem-fixedwidth-undecidable.md`](../../notes/reviews/briefs/2026-10-07-hostmem-fixedwidth-undecidable.md)
- **plan**: [`plans/active/architecture-solver-host-memory.md`](../../plans/active/architecture-solver-host-memory.md)
- **date**: 2026-10-07
- **commit**: `056e3576` (feature/sern-design)
- **codex**: effort `xhigh`, 6.3 min, rc=0
- **結論**: ?
- **extra**: `case/66.hostmem_regression/fixedwidth_c44dual_ckpt100/RESULT.txt`

本文は codex の最終メッセージをそのまま転記。結論は**仮説**であり、提案 A/B で確かめる。採否・反映は plan 側 (§5.1 担当列 F の判断欄) に書く。

---

採否表

| 諮問事項 | 判断 | 根拠と対案 |
|---|---|---|
| 固定幅試験を「判定不能」で終了する | **採用** | 凍結した SHA は再照合できた。保存結果では base も4/6本が超過し、[判定実装:274](/home/sano/work/forge-sern-design/case/66.hostmem_regression/fixedwidth_eval.py:274)どおり。今回の幅・本数は変更しない。 |
| 現状で限定事項を付けて `done` にする | **却下〔Major〕** | [plan:217](/home/sano/work/forge-sern-design/plans/active/architecture-solver-host-memory.md:217)には物理的な許容誤差が未定義と明記され、§8の完了条件も満たしていない。**互換性・メモリ削減の確認済み部分と、数値的受入れの未決着を分けて記録する。** |
| 両側2/6なので頻度差は解消した | **却下〔Major〕** | [保存結果:33](/home/sano/work/forge-sern-design/case/66.hostmem_regression/fixedwidth_c44dual_ckpt100/RESULT.txt:33)。独立二項標本を仮定しても、2/6の母比率の95%正確区間は約 **4.3–77.7%**。同数という観測では頻度差を十分に制限できない。対案は検出幅と本数を固定した別の頻度試験。 |
| 巨大値は診断量の分母の問題なので無害 | **却下〔Major〕** | [追加記録:25](/home/sano/work/forge-sern-design/case/66.hostmem_regression/fixedwidth_c44dual_ckpt100/notes.txt:25)では、`run_0106` の節点19954で `roQ1_0` が他5本の約2.1倍、`roQ2_0` が約4.8倍。診断比だけでなく保存量も変動している。補正量・保存量・履歴の許容差を定義して評価する。 |
| ばらつきの小さい別構成の合格で今回を代替する | **却下〔Major〕** | 凝縮 dual-time 固有の更新・FCT・射影は [main.cpp:2756](/home/sano/work/forge-sern-design/solver_density_cuda/main.cpp:2756)以降を通る。他構成の合格では、この経路の未解決差を除外できない。次の試験も同じ構成にする。 |

**結論:** R1–R3 の `done` は保留し、`c44dual_ckpt100` の**step 100 における閾値超過率**を比較する独立 A/B を、別試験として事前登録して一度だけ行う。

**第1仮説:** 既存の実行間変動によってモーメントの異なる状態が現れ、微小分母を持つ補正診断がその差を増幅している。**確度: 中**。

- **根拠:** 同じ入力の変更前バイナリでも、`case/66.hostmem_regression/run_0196_fw_c44dual_b2/` と `run_0200_fw_c44dual_b4/` に約1.45–1.46e23が記録されている。[run別記録:5](/home/sano/work/forge-sern-design/case/66.hostmem_regression/fixedwidth_c44dual_ckpt100/result/describe.txt:5)。実際の dual-time 終端では射影が `roQ1/roQ2` を書き換え、その補正を分母床1e−30の比で記録する。[condensationRealizability_d.cuh:377](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/condensationRealizability_d.cuh:377)。
- **反証条件:** 独立標本で事前に決めた大きさのビルド間頻度差が確認されたら、「既存の変動だけで説明できる」という説明を棄却する。ただし、微小分母による増幅機構そのものの反証にはならない。

**第2仮説:** ホスト確保・初期化順序の変更が、当該状態の出現確率を変えている。**確度: 低、未確認**。カーネルソース不変更と初期出力一致は支持材料だが、確率変化の除外にはならない。

**判別 A/B:** 変更因子は**実行バイナリだけ**。同じ `c44dual_ckpt100` 入力・環境・GPUで、base/new **各128本、各100 step**。128組の実行順を組内で無作為化し、投入前に固定する。既存標本は混ぜない。

- 主評価量は `E = 1[max(VALUE/condClampCorrQ_0 in res_100.h5) > 1e22]`。1 runを1標本とし、節点や内反復を独立標本に数えない。
- 診断用の検出幅を **頻度差±20パーセントポイント**、判定を **差 p_new−p_base の90% Newcombe–Wilson区間**で固定する。これは**物理的受入幅ではない**。
- **A：区間全体が −0.20〜+0.20 内**なら、この終点事象に関する「大きな頻度差」を退ける。
- **B：区間全体がその外側**なら、「既存変動だけ」という説明を退け、変更起因の確率変化を追う。
- 境界をまたぐ場合は判定不能で終了する。有限標本で必ず二択になる保証はできない。本数を追加して結論を取りに行かない。

この本数は、独立二項モデルで両側の真の率が1/3なら、上記Aになる確率が約92.8%となることを計算確認した。保存記録の12本83秒からは、実行約30分の規模となる。全出力・NaN・構造の検査は維持し、主評価量以外の異常も記録する。

**この試験のAだけで `done` にはしない。** 数値的受入れには、newだけが超過した `roQ0_0`・`roQ1_0_fctH` を含め、保存量・履歴・補正量について用途に基づく許容差が必要である。それが未定義のままなら、別構成や頻度試験の合格を追加しても受入れ判断は完成しない。

**やらない方がよいこと:** 今回のTの拡大、合格するまでの追加反復、newだけの超過量の除外、凝縮や射影を無効にした試験による代替、非有意を同等性とする判断。

**呼び出し側の前提への異議:** 「その状態に入った本数」は定義を修正する必要がある。抽出器は [fixedwidth_describe.py:46](/home/sano/work/forge-sern-design/case/66.hostmem_regression/fixedwidth_describe.py:46)で `res_100.h5` だけを読み、診断配列は [condensationUpdateLimiter_d.cuh:138](/home/sano/work/forge-sern-design/solver_density_cuda/cuda_forge/condensationUpdateLimiter_d.cuh:138)で上書きされる。したがって、**100 step中に一度でも発生した割合ではない**。

現時点の限定事項は、次の書き方が妥当である。

> ホストメモリ削減と試験した初期化・出力・checkpointの互換性を確認した。凝縮 dual-time の固定幅試験は、変更前対照も超過して判定不能だった。step 100の巨大診断値は両ビルドで観測されたが、出現確率、保存量への影響、および物理的に許容できる差は未確定である。

`procedures` には、**決定的状態の厳密比較／許容差を定義した数値比較／非決定的事象の頻度比較**を別判定として残すべきである。校正標本と評価標本を分離し、観測時点・検出幅・本数・判定不能時の停止を事前固定する。

**不足情報:** 対象runの生HDF5・残差CSV・個別VERDICTはローカルに存在せず、数値は保存された評価記録までの確認である。特に不足しているのは、保存量・checkpoint・補正量の**物理的許容差**。収束・準定常性は今回独立確認しておらず、それらの主張はしない。

ファイル変更・`forge`実行なし。**plan未反映**。呼び出し側で §5.1・§6.3・§8へ反映する診断として返す。
