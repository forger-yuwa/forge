# 諮問記録 (diagnostician、2026-10-01): 共役平板 C の結果の解釈

条件 7。diagnostician (Fable) に諮った内容と結論の写し。codex 代替に回せる形で残す。
- 入力: plan `plans/active/boundary-cht-conjugate-flat-plate.md` §4.6・§5.1 #1g・#2d、run_0020・run_0022〜0027 の判定出力、`lim0_chain.sh`・`lim0_regate.sh`、`eval_conj.py` の修正 `695e5de`。
- 結論: C2 n16 の判定不能は参照側の U による → n16 2 本に参照 4 水準の再判定を 1 回だけ回し、事前登録した 3 分岐 (PASS → 6 本全合格 / FAIL → n32・n64 の限定結果 / 判定不能 → 同) に従う。限定で閉じる枝は承認済み。C1 の軸方向熱量ゲートの修正は登録文 (発注元 §6) への整合として採用 (事後緩和ではない)。
- 報告で必ず書く限定: limiter 0 は事後改訂設定 (生産リミッタ経路の検証ではない)、cfl 0.5 は判別条件 4 で決まった (cfl 2 は未検証)、IC は limiter 2 の停滞場からの同一格子 restart (一様流 IC の対照なし)、④ の訂正、C1 ゲートの整合は FAIL を見た後、C2 n16 の扱い、流れ自体の正しさは未検証、収束は下降中のまま PASS。
- 実行結果: run_0025 L3 FAIL (2.354 % + 0.979 %)、run_0024 L3 PASS → C は n32・n64 の限定結果。
