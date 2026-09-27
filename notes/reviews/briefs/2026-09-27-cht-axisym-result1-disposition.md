# 諮問ブリーフ: 軸対称 fem2d の result 1 巡目 (M1 / m1) の反映を確認し、accepted へ移してよいか (2026-09-27)

AGENTS.md 条件 **5** (codex の Major の採否)。

## 読んでよいファイル
- レビュー: `notes/reviews/2026-09-27-boundary-cht-axisymmetric-fem2d-result.md`
- 差分: `git show d4801bd4` (conjugateWall.cpp の再開契約・check_cht_balance.py・methods・README・template)
- 負例の証拠: `case/61.conjugate_annulus/vax4_contract_evidence/README.md` (9 件、C++ と評価器の一致)
- plan `plans/active/boundary-cht-axisymmetric-fem2d.md` の §5.1 #9・#12 と §6.1 の result 行

## 採否 (全件採用・却下 0)
- M1: 旧形式として受理するのを「3 属性すべて無い かつ 平面」に限定、部分欠落は C++・評価器とも拒否。負例 9 件で一致 (平面の旧状態・3 属性は受理、各 1 属性欠落は拒否、軸対称 3 属性は受理、各 1 属性欠落は拒否)。
  平面の受理 2 件は復元後 step 0 で固体温度の安全停止 (流体場を引き継がない再開の既知の起動の罠) — 契約判定とは無関係と読んでいる。
- m1: methods/boundary.md (軸対称の面内伝導を未定義から外す)、plans/README、case/61 README (warmup 20000、感度の保留を決着)、両 case の template (warmup 20000)、case/62 README (非一様は純伝導 IC)。

## 諮りたいこと (推奨を 1 つに)
1. この反映で M1/m1 は閉じたか。不足があれば 1 つ。
2. `status: done` → `plans/accepted/` へ移してよいか (result 2 巡目が要るか)。
