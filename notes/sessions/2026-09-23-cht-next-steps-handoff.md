# 引き継ぎ — 冷却翼 CHT の検証が一段落した時点で、次にやること (2026-09-23)

**この文書は写しとポインタだけ**。残作業の正本は各 plan の §5.1
([AGENTS.md「方針の反映トリガは『実装』ではなく『決定』」](../../AGENTS.md))。ここが古くなっていたら plan を信じること。

## 0. いまの状態 (やり直さないために)

| 対象 | 状態 | 正本 |
| --- | --- | --- |
| 外部弱連成 CHT (Phase 1) | **動く**。C3X・Mark II とも 1 パスで収束。固体は報告の公開値だけ (合わせ込みなし、`solid_published.json`) | [`plans/active/boundary-conjugate-heat-transfer.md`](../../plans/active/boundary-conjugate-heat-transfer.md) |
| 公知データ検証 (NASA CR-168015) | **一段落**。C3X run 108 / Mark II run 42、固定壁温と連成の両方 | [`case/53…/README.md`](../../case/53.c3x_vane_cht/README.md), [`case/54…/README.md`](../../case/54.markii_vane_cht/README.md) |
| SU2 とのコード間照合 | 済 (同一メッシュ・同一節点)。乱流で 2.1 %、遷移モデルどうしで領域平均 −0.4〜+1.4 % | 同上 |
| 遷移モデル ($\gamma$–$Re_{\theta t}$) | 実装・平板検証・翼への適用まで済み、**2026-09-23 ユーザ判断で凍結** | [`plans/active/turbulence-transition-lm2009.md`](../../plans/active/turbulence-transition-lm2009.md) |
| 報告 | Artifact「Cooled Vane CHT Validation」**V42** (09 章が未解決一覧) | `notes/reports/vane-cht-validation/` |

**凍結の理由 (再開判断の材料)**: 遷移モデルは動くが、**翼の全体一致は入口の乱れの減衰の速さ (入口粘性比) に支配され**、
報告がその値を与えない。Mark II で粘性比 10 → 40 にすると全体 rms 32.0 → 27.3 %、正圧面 −27.7 → −9.2 % と動く一方、
文献が振っているモデル定数 ($Re_{\theta t,min}$) は 20 → 200 で 1 ポイントしか動かない (plan §6.5、報告 図 26)。
**つまりモデルを直しても「合った」と言えない構造**にある。

## 1. 次にやること (推奨順)

### (1) ソルバ内連成 (Phase 2) — 本命

外部ループは**流体計算を 30〜70 回**呼ぶ。1 回 4000 step なので C3X で 27 万 step、Mark II で 46 万 step。
すきま・キャビティ (下の (3)) に持っていくとこの回数が効く。設計は plan §4.6 に書いてある
(`wall_isothermal` + `ints: {conjugate: 1}`、node 限定、host 側で固体を double で解き $K$ step ごとに `bvar_d["Ts"]` を更新)。
界面契約 (`iface_q_eff` = 物理境界流束 + 拘束反力) と `fem2d` バックエンドは Phase 1 で実証済みなので、**そのまま使える**。

- 合格条件は plan §4.6/§6: **Phase 1 と同じ壁温に落ちること** (C3X で 587.09 K・Mark II で 565.09 K が現行値)、
  界面残差が熱負荷の 0.1 % 以内、速度向上は**主張でなく測定項目** (§4.5 の注意: 流体の熱場緩和は Phase 2 でも残る)。
- 副産物として、報告 09 章の「連成ループのコストとノイズ」が閉じる。

### (2) 引用している数字が定常解のものか、を 1 回だけ確かめる

**全 run が `NOT CONVERGED (stalled/plateau)`** で、引用量は `check_quasisteady` の STEADY に依っている。
C3X は衝撃足が動き続けている。報告 09 章の 2 番目。**生産設定を 1 つ選んで dual-time で回し、壁熱流束を時間平均**して、
定常擬似時間の値と比べる。ここが通れば「未収束だが量は定常」という但し書きを外せる。
(平板側では倍精度対照で床が丸めだと確かめた: `case/57…/run_0016_t3a_lm_double`。翼は未実施。)

### (3) 本来の適用先へ渡す — 極超音速すきま / 深いキャビティ

CHT を作った動機はここ ([`notes/investigations/cht-validation-case-survey.md`](../investigations/cht-validation-case-survey.md) §4)。
翼で得た誤差幅を付けて case/49 (環状深キャビティ) と case/50/51/55/56 (すきま加熱) に適用する。
**注意**: すきま側の文献 3 本は定量的な受け渡しに使えないと分かっている
([`plans/active/case-hypersonic-gap-heating-validation.md`](../../plans/active/case-hypersonic-gap-heating-validation.md))。
したがって (3) は「検証」でなく**予測に誤差幅を付ける**作業になる。誤差幅の出所は翼の領域別 bias (正圧面 ±10 %・遷移後 ±20 % 程度)。

### (4) 残っている小さいもの (どれも独立、手が空いたときに)

| 項目 | なぜ | 正本 |
| --- | --- | --- |
| **未知の設定キーを拒否する** | `turbulence.kInf` / `omegaInf` は**ソルバに存在しないのに黙って無視**され、4 run が同一結果になってから気付いた。同じ事故が起きる | CHT plan §5.1 #52、報告 09 章 |
| 固体単体の −7〜−10 K バイアス | 実測熱流束を課した固体だけで 16〜17 K rms。金属熱伝導率は外部出典 (報告に値が無い)、最終孔が計測範囲外。固体メッシュと $k_s$ を振っていない | 報告 09 章、CHT plan §5.1 |
| forge と SU2 の残り 2.6 % / 2.8 % | 壁面流束どうしで比べた差。原因未特定 | CHT plan §5.1 #41 |
| Mark II の壁方向メッシュ収束 | 遷移後の bias が壁節点数とともにまだ上がっている (340/480/680) | CHT plan §5.1 #47 |
| 壁半 CV の粘性仕事が相対 1 次精度 | MMS で確認済み。第一層 2 µm で報告熱流束の 2.8 %、1 µm で 1.4 %。実翼での符号・大きさは未評価 | [`notes/investigations/2026-09-22-mms-node-face-weight.md`](../investigations/2026-09-22-mms-node-face-weight.md) |

### (5) 遷移モデルを再開するなら (凍結解除の条件)

**合わせ込みを続けない**。効くのは次のどちらか。

- (a) **壁近傍の乱流量まで測っている試験に替える** — 入口の乱れの減衰が既知なら、モデルの誤差と入力の誤差を分けられる。
- (b) **遷移位置を実測に固定した計算** — 遷移の機構でなく、遷移後の熱伝達だけを見る。報告 09 章が挙げている検査。

未完の検証は plan §5.1 の #10 (T3B の実験照合、ERCOFTAC サーバに接続できず未入手)、#11 (両コードの報告量の時系列、SU2 側は済)、
#12 (1 µm・実測壁温での SU2 照合)、#13 (Mark II の粘性比 1)、#14 (CHT)、#16 (比較量を壁温へ)。

## 2. 触るときに知っておくこと (罠)

- **run は長い**。遷移つきの翼は 15 万〜30 万 step 回さないと正圧面の符号すら変わる (Mark II 粘性比 30 の正圧面は 6 万 step で +18 %、20 万 step で −22 %)。
  `h_series.py` + `check_quasisteady.py --series-csv` で止めること。
- **壁解像ゲート**は遷移モデルの有無で結果が変わる (Mark II 1 µm は遷移ありで PASS、なしで FAIL — 完全乱流の方が壁せん断が大きい)。
- **報告の $h$ は導出量** (実測壁温を境界値として固体を FEM で解いたもの)。文献は壁温で比較している。CHT が回るなら壁温で比べる方が筋が良い。
- メモリの関連項目: [[transition-lm2009-status]], [[vane-cht-transition-is-the-error]], [[vane-cht-hole-placement-and-loop]], [[model-tiering-workflow]]。
