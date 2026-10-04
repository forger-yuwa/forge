# codex レビュー: tooling-energy-balance-diagnostics (plan)

- **plan**: [`plans/active/tooling-energy-balance-diagnostics.md`](../../plans/active/tooling-energy-balance-diagnostics.md)
- **stage**: `plan`
- **date**: 2026-09-19
- **commit**: `8e0523eb` (feature/sern-design)
- **codex**: effort `high`, 3.8 min, rc=0
- **判定**: **GO-with-changes**, 指摘 C0/M3/m0
- **focus**: 1 巡目 (NO-GO, C0/M7/m1) を全件採用した改訂稿。M1 収支式と拘束反力 (§4.1), M2 採取位相 (§4.2), M3 全加算経路 (§4.3), M4 CV/面スキーマ (§4.4), M5 opt-in 契約 (§4.5), M6 CV 選択・周期・単位 (§4.6), M7 解除試験 4 本立てと符号・丸め誤差限界 (§6), m1 現状説明の訂正 (§1) が実質的に閉じたかを検証してほしい。実装に進めない穴が残っていれば指摘し、無ければ GO 判定とその条件を明示してほしい
- **extra**: `plans/active/case-hypersonic-gap-heating-validation.md`

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**判定: GO-with-changes**  
目的と初版のスコープは妥当で、既存の出力関連 plan とも重複していません。  
旧 M2・M5・M6・m1 は解消と判断します。旧 M1・M3・M4・M7 は、以下の実装前修正が必要です。

1. **Major — 壁別の実効熱量が未定義で、共有角ノードの拘束反力を二重計上できる**

   **根拠:** [plan:58](/home/sano/work/forge/plans/active/tooling-energy-balance-diagnostics.md:58) は「物理境界流束 + 拘束反力」「別名で定義する」としていますが、集計式と帰属規則がありません。[gmshReader.hpp:2292](/home/sano/work/forge/solver_density_cuda/mesh/gmshReader.hpp:2292) では角ノードが複数 `bcond.iCells` に重複します。一方、[nodeWallDirichlet_d.cu:170](/home/sano/work/forge/solver_density_cuda/cuda_forge/nodeWallDirichlet_d.cu:170) がゼロ化するエネルギー残差は、その CV に一つだけです。

   したがって、各壁の `iCells` で `C_i=-R_i^{raw}` を積分すると、共有ノードの反力を重複計上します。異なる壁温の競合を避けても、**同温度の側壁と底面に反力をどう分けるか**は残ります。面の `physID` だけでは解決しません。

   **対案:** §4.1・§4.4 に次を確定してください。

   - 流体への供給を正とする実効熱量は、重複しない拘束 CV 集合 \(W\) と対応する物理壁面集合 \(B_W\) に対して  
     \[
     Q_{\mathrm{eff}}(W)=-\sum_{f\in B_W}F_f+\sum_{i\in W}C_i
     \]
     と明記する。外向き `F` と供給方向 `C` は、そのまま加算できません。
   - 第一内部列への供給は、壁 CV と内部 CV の間の**全数値流束**として別定義する。壁部分領域では接線方向輸送・ソースもあるため、上式との無条件な同一視を禁止する。
   - `C_i` は CV ごとに一度だけ保存する。複数壁が共有する反力は初版では**接合部の別勘定**とし、壁別には無理に配賦しない。壁∩出口などの非壁面流束も保持する。
   - 共有角を含む試験で、壁別勘定と接合部勘定の合計が領域全体に戻ることを検証する。

2. **Major — 「全加算経路」の表に欠落があり、受理する機能の集合も未確定**

   **根拠:** [plan:77](/home/sano/work/forge/plans/active/tooling-energy-balance-diagnostics.md:77) の表には化学反応熱がありません。しかし [chemistry_d.cu:108](/home/sano/work/forge/solver_density_cuda/cuda_forge/chemistry_d.cu:108) は `res_roe += V*Qdot` を実行し、[main.cpp:1402](/home/sano/work/forge/solver_density_cuda/main.cpp:1402) から残差組立て中に呼ばれます。「未対応なら拒否」は正しい方針ですが、何を未対応とするかがまだ決まっていません。

   また、[speciesTransport_d.cu:502](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:502) も `res_roe` を変更しますが、これは [main.cpp:1610](/home/sano/work/forge/solver_density_cuda/main.cpp:1610) の陰的連成補正です。空間残差の物理ソースとして採取すると誤ります。

   **対案:** §4.3 を、各経路について「有効条件／採取対象／拒否条件」を持つ表にしてください。初版は発注元に必要な**非反応 TP・層流／低 Re SST**を中核とし、化学・凝縮・壁モデルなど未検証の組合せは具体的な設定条件で拒否するのが妥当です。

   採取区間も「境界状態・物性・勾配の評価後から、壁残差射影直前まで」と固定し、陰的連成補正を除外してください。`sstEnergyIncludesK` を受理するなら、[ransTransport_d.cu:163](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransTransport_d.cu:163) の分割保持に従い、**残差は `roe+roK` の式、保存配列 `roe` は平均流エネルギー**という対応をスキーマに明記する必要があります。

3. **Major — 解除試験は整理されたが、非ゼロ経路の網羅と定量的な合格式がまだない**

   **根拠:** [plan:143](/home/sano/work/forge/plans/active/tooling-energy-balance-diagnostics.md:143) 以降は試験の分類であり、ケース条件・期待値・許容式は未定です。低 Re SST を実行するだけでは、[ransSource_d.cu:230](/home/sano/work/forge/solver_density_cuda/cuda_forge/ransSource_d.cu:230) の条件付きエネルギーソースを検証できません。一様組成の TP ケースでは、[speciesTransport_d.cu:263](/home/sano/work/forge/solver_density_cuda/cuda_forge/speciesTransport_d.cu:263) の種拡散もゼロになり得ます。

   また、丸め誤差限界を `Σ|F|` だけから作ると、[bodyForce_d.cu:43](/home/sano/work/forge/solver_density_cuda/cuda_forge/bodyForce_d.cu:43) などのソース加算を含む残差の誤差を覆えません。

   **対案:** §6 に以下を事前登録してください。

   - **受理する加算経路ごとに非ゼロになる試験**を割り当てる。組成勾配、`k` 勾配、非平衡の `Pk-Dk`、体積力仕事など、対象項が実際に発動したことも検査する。拒否対象には起動拒否試験を用意する。
   - 組立て試験の基準を、各 CV に実際に加えた項 \(a_{ij}\) に対して  
     \[
     \left|R_i^{raw}-\sum_j a_{ij}\right|\le B_i,\qquad
     B_i\sim\gamma_{n_i}\sum_j|a_{ij}|,\quad
     \gamma_n=\frac{nu}{1-nu},\ u=2^{-24}
     \]
     と具体化する。ソース・保存時の丸め・集計誤差も含め、領域和は FP64 で評価する。
   - **組立て誤差の範囲内**と、**壁熱量を必要精度で検算できる**を分ける。誤差限界が対象熱量の誤差予算を超える場合は、後者を「判定不能」とする。
   - 純伝導は上下壁それぞれの符号・期待熱量・離散化誤差許容を固定する。`C=-Rraw` から作った恒等式が閉じるだけでは、熱量の正しさの独立検証になりません。

旧 run の扱いの訂正は確認できました。`case/24.laminar_channel_bl/run_isoT_condN_node/` に判定器を再実行した結果は **`NOT CONVERGED (stalled/plateau — needs scheme change, not more steps)`**。`rms_roe=4.42e-4`、`rms_roUy=8.27e-7` はともに rising でした。最終壁ファイル `res_wall_bot_3_10000.h5`・`res_wall_top_4_10000.h5` の `qwall` も各17点すべてゼロです。[既存 run 索引](/home/sano/work/forge/case/24.laminar_channel_bl/README.md:35) の結果を解除試験に継承しない判断は正しいです。

**推奨は、初版スコープを維持し、上記 1 → 2 → 3 の順で仕様を確定してから実装することです。** 採取位相の分離、CV／面の別スキーマ、opt-in、周期・軸対称の拒否という骨格は採用できます。再設計は不要ですが、現稿の「確定する」という残作業を、そのまま実装者の判断に渡す段階ではありません。

ファイル変更なし。上記提案は **plan 未反映**です。

指摘数: Critical 0 / Major 3 / Minor 0
