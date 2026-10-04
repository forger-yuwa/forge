# 諮問ブリーフ: #1b 事前登録の最終確定 (再正規化の許容) (2026-10-02)

作業ツリー `/home/sano/work/forge-species` (HEAD は git log -1)。plan `plans/active/condensation-two-phase-transport.md` §5.1 #1b (前回の判断を反映済み) と #1b-pre (実装結果)。前回諮問 `notes/reviews/2026-10-02-twophase-1b-preregistration-diagnose.md`。設計メモ `notes/investigations/2026-10-02-twophase-diffusion-kernel-design.md` §16 (特に §16.5 許容の案)。

## 問い

1. 再正規化の許容を確定してほしい。案: 末尾 10 % の全更新 (monitorInterval をその窓に合わせる) で (i) 局所係数偏差 max_i|f_i−1| ≤ 2·n_s·ε₃₂ (case/16 n_s=2 で 4.8e-7)、(ii) ρY_w・ρg・ρQ2・ρQ1・ρQ0 の成分別相対補正 Σ|Δq|V/ΣqV ≤ 2·n_s·ε₃₂ (総量 0 は絶対 0、総液量比は使わない)。根拠: 定常の固定点では化学種の増分と δρ が 0 なので f−1 は格納丸め ~(n_s/2+1)ε₃₂、Δq = (f−1)q なので成分別相対補正は max|f−1| で上から抑えられる。偏差が下がらない run は未収束扱い。局所と積分のどちらで判定するか (両方か)、係数 2·n_s、窓。
2. #1b-pre の 5 項目 (メモ §16) で前回 NO-GO の不足 3 点 (A の旧作用素監査、再正規化の成分別計測、全更新集計) は満たされたか。満たされていれば、前回の判断 + 上の許容で #1b の事前登録を「確定」としてよいか。追加で欠けているものがあれば 1 つ。

## 実測 (2 step、ほぼ乾いた初期場、判定外)

- A 監査の自己検査 ≤0.601/0.630 ε·max A (全成分)、2 step で NOT CONVERGED (正しい)。既定キー OFF は G0 4 構成バイト一致、res_0 120 データセット一致。
- 再正規化: B max|f−1| 7.87e-5/6.15e-5、成分別相対補正 ρY_w 2.3e-7/1.9e-7・ρg 1.1e-7/1.4e-7・Q 1.2〜1.6e-7。A は液・Q に掛からない (0)。局所と積分で約 2.5 桁違う。既存の総液量比 renorm は乾いた場で 0.4 と無意味。
- θ 全更新集計: 更新 θ<1 は A 24・B 22 セル、θ_src<1 は 0。

## 禁止事項 (厳守)

- ファイルを変更しない。`*.log`, `residual_history.csv`, `res_*.h5`, `*.vtu`, `plans/README.md` を読まない。推奨は問いごとに 1 つ。

## 読んでよいもの

- plan 全体、前回諮問、設計メモ §14–§16、`solver_density_cuda/cuda_forge/speciesTransport_d.cu` の再正規化と監査 (grep `renorm`・`audit`)、`condensationTransport_d.cu` の新ログ行
