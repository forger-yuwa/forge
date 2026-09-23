# codex レビュー: boundary-conjugate-heat-transfer (result)

- **plan**: [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md)
- **stage**: `result` (diff base `main`)
- **date**: 2026-09-23
- **commit**: `1ba87d56` (feature/cht-phase2-fem2d)
- **codex**: effort `high`, 9.4 min, rc=0
- **判定**: **NO-GO**, 指摘 C0/M5/m1
- **focus**: 1 巡目 (notes/reviews/2026-09-23-boundary-conjugate-heat-transfer-result.md) の Major 7 + Minor 2 への対処が正しいかを最優先で見る。特に (1) 残差補正形 delta = -A_old^{-1} r による更新が固定点を現在物性で保存しているか、(2) 積分済み荷重 iface_Qf_eff の受け渡しが角で保存するか、(3) G-if の内部残差必須化と 2N 待機、(4) 再開経路の同値 (run_0152/run_0153 で差 0.0084 K)、(5) 安全停止の誤爆対策が発散を見逃さないか。再判定は run_0151_fixed_avg42 (C3X) と run_0035_fixed_avg42 (Mark II)

本文は codex の最終メッセージをそのまま転記 (指摘の採否と対応は plan 側の「レビュー記録」表と残作業表に書く)。

---

**NO-GO**

C3XのG-if・準定常の合格は再現でき、固体内部残差も改善しています。  
ただし、温度の精度差による偽の固定点、ゲートの誤合格、再開・安全停止の欠陥が残っており、`accepted`への移動は認められません。

1. **Major — 残差に`D_f`が残り、量子化による停滞をPASSにできます。**

   [残差組立て](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:701)はdoubleの固体温度`u`と、float32の`Ts`を混用します。[更新後](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:778)も両者を異なる精度で保持するため、実際に判定する残差は
   `K(u)u − b − Qf + Df(Eu − Ts)`です。「`D_f`の項が厳密に相殺する」というコメントは成立しません。

   **再現:** 4節点・2三角形の伝導問題で、界面荷重各0.5 W、`Df=1e5 W/K`とすると、正解602 Kに対して約600.000005 Kで停滞しました。物理的な界面不釣合いは**99.99975%**ですが、現行ゲートは80更新を**`VERDICT: PASS`**と判定します。今回のC3X保存場では混入項は最大約0.66 W/m²で小さいものの、一般的な固定点保存は未達です。

   **対案:** 更新・判定に使う物理残差を`K(u)u−b−EᵀQf`として直接組み、界面温度の精度整合を明示してください。量子化で更新できない場合は、物理残差を残したまま合格させないこと。

2. **Major — G-ifの待機が不足し、最新の不良データを捨ててPASSします。**

   [check_cht_interface.py:115](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_interface.py:115)は`update >= 2N`を使います。仕様は「充填後さらに2N」なので、初回充填がN更新目なら開始は約3N更新目です。

   **再現:** `N=42、n_consec=80`で、開始から163更新しかない入力がPASSしました。仕様上、この時点では待機後の行は38行しかありません。さらに末尾へ「バッファ再初期化・残差1e9」の行を追加しても、`rows_ok`のフィルタで除外され、過去の80行からPASSしました。

   **対案:** 充填完了時点と再初期化を追跡し、**実際の末尾に連続する80更新**すべてが待機・残差条件を満たす場合だけ合格にしてください。

3. **Major — 再開の互換性と更新時刻を保証できていません。**

   [再開時の検査](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:476)は界面ハッシュと節点数だけです。しかし[ハッシュの生成](/home/sano/work/forge-cht/solver_density_cuda/tools/solid_mesh_to_h5.py:168)は順序を除いた界面座標のみで、内部節点順・接続・物性・冷却条件を識別しません。

   **再現:** C3Xの内部2節点を入れ替えると、節点数と界面ハッシュは同じまま、保存温度の対応が**214.33 K**ずれます。現行検査では拒否できません。

   また、保存した`step`を読まず、[更新判定](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:826)を再開後のローカルstepで行います。`interval=50`の倍数でない位置からの再開では更新位相が変わります。[状態出力](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:966)も最終stepの強制保存がなく、流体最終場と固体チェックポイントが異なる時刻になり得ます。

   **対案:** 節点順を含む固体入力ハッシュと連成設定、累積step・更新位相を保存／検査し、流体と固体を同一時刻で保存してください。`run_0152/0153`の差**+0.0084 K**は再現できましたが、`warmup=0`かつ50の倍数で分割した経路だけの根拠です。

4. **Major — 保存荷重の修正が外部ループとG-consに届いていません。**

   ソルバ内の`iface_Qf_eff`直接転送は正しい修正です。一方、[cht_loop.py:330](/home/sano/work/forge-cht/solver_density_cuda/tools/cht_loop.py:330)と[check_cht_balance.py:132](/home/sano/work/forge-cht/solver_density_cuda/tools/check_cht_balance.py:132)は、依然として`iface_q_eff × 固体集中辺長`を使います。

   **実測:** `case/54.markii_vane_cht/run_0035_fixed_avg42/res_wall_5_40000.h5`の後縁角では、直接荷重**1.65788 W/m**に対し、このPython経路は**2.15344 W/m（+29.89%）**を作ります。Phase 1とPhase 2が同じ界面契約になっていません。

   **対案:** 外部ループ・平均処理・G-consも積分済み荷重へ統一してください。G-consのソルバ内連成対応では、初期`wall_profile`から復元した固体ではなく、対応時刻の保存固体状態を使う必要があります。

5. **Major — 増加開始値が0だと、発散検知が永久に成立しません。**

   [安全停止条件](/home/sano/work/forge-cht/solver_density_cuda/conjugateWall.cpp:788)は開始値0を保存した後、`dTgrowStart > 0`を要求します。

   **再現:** `dTw = 0 → 0.001 → 0.001×1.1 → …`では、81回連続増加して**2.0484 K/更新**になっても、この停止条件は成立しません。正の開始値では、許容0.01 Kを超えたところで停止することも確認しました。

   **対案:** 増幅の基準を`max(開始値, 登録許容)`などで定義し、0からの増大も検知してください。許容内の微小振動を除外する方針自体は維持できます。

6. **Minor — 現在仕様と残作業表の同期が不十分です。**

   [methodsの設定例](/home/sano/work/forge-cht/methods/boundary.md:462)には、今回必須にした`tol_solid`がありません。同文書には旧来の直接求解・`D_f`凍結説明が残り、[plan §4.6a](/home/sano/work/forge-cht/plans/active/boundary-conjugate-heat-transfer.md:314)も荷重を`q_eff × 集中辺長`と定義しています。

   **対案:** 修正記録だけでなく仕様本文・設定例・索引を更新し、V2・V3・V6などの未完了事項とPhase 2の承認範囲を整理してください。

再実行した判定は以下です。準定常許容は壁温`0.0006`、熱量`0.002`です。

| run | G-if | 準定常 |
|---|---|---|
| `case/53.c3x_vane_cht/run_0151_fixed_avg42/` | **PASS**。内部残差最大`2.1783e-4 W/m` | **ALL STEADY** |
| `case/54.markii_vane_cht/run_0035_fixed_avg42/` | **NOT CONVERGED**。内部残差のみ許容内 | 壁温平均・熱量は**STEADY**、最大壁温は**OSCILLATING** |

両runの流体残差は**`NOT CONVERGED (stalled/plateau)`**です。保存された流体・固体場にNaN/Infはありませんでした。`test_solid_shell.py`と`test_solid_fem2d.py`はともに**`PASS (all)`**でした。runの対応は各case READMEの[ C3X台帳](/home/sano/work/forge-cht/case/53.c3x_vane_cht/README.md:778)・[Mark II台帳](/home/sano/work/forge-cht/case/54.markii_vane_cht/README.md:82)で確認しました。

推奨は、**`active`に留め、上記の順に修正・反例の回帰試験・V4b再判定を行うこと**です。ファイルは変更していません。判定器の書込みも抑止しており、提案は**plan未反映**です。

指摘数: Critical 0 / Major 5 / Minor 1
