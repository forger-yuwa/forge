#!/bin/bash
# lim0_chain.sh の後判定のやり直し (2026-10-01):
#  - check_cht_interface.py は run ディレクトリに CHT_INTERFACE_VERDICT.txt を自分で書くため、チェーンでは 2 回目 (旧 1e-9) が
#    1 回目 (登録式の tol_solid) を上書きしていた。別ファイルに分けて残す。
#  - eval_conj.py の C1 軸方向熱量ゲートを登録文 (発注元 §6「C1 は判定しない」) に合わせたので主判定を再実行。
T=~/forge-graetz/solver_density_cuda/tools; C=~/forge-graetz/case/65.conjugate_flat_plate
cd ~/forge-graetz && GIT_SSH_COMMAND="ssh -i ~/.ssh/forge-deploy -o StrictHostKeyChecking=no" git fetch -q origin feature/cht-conjugate-benchmarks && git checkout -q -f -B feature/cht-conjugate-benchmarks FETCH_HEAD; git log --oneline -1
for spec in "run_0022_c1_n64_lim0 5.322e-07" "run_0023_c2_n64_lim0 7.832e-07" "run_0024_c1_n16_lim0 2.311e-06" "run_0025_c2_n16_lim0 3.409e-06" "run_0026_c1_n32_lim0 1.094e-06" "run_0027_c2_n32_lim0 1.610e-06"; do
  set -- $spec; r=$C/$1
  python3 $T/check_cht_interface.py $r --phys-id 5 --eps-abs 1.0 --eps-rel 1e-3 --dt-k 5e-4 --tol-solid 1e-9 --n-consec 80 > /dev/null 2>&1
  cp $r/CHT_INTERFACE_VERDICT.txt $r/CHT_INTERFACE_VERDICT_old_tol1e-9.txt
  python3 $T/check_cht_interface.py $r --phys-id 5 --eps-abs 1.0 --eps-rel 1e-3 --dt-k 5e-4 --tol-solid $2 --n-consec 80 > /dev/null 2>&1
  sed -i "1i # 登録式の tol_solid $2 (plan §4.6、事後改訂)。旧転記値 1e-9 の結果は CHT_INTERFACE_VERDICT_old_tol1e-9.txt" $r/CHT_INTERFACE_VERDICT.txt
  (cd $C; python3 ../64.conjugate_pipe_wall/eval_conj.py C $r > $r/EVAL_CONJ.txt 2>&1; echo "rc=$?" >> $r/EVAL_CONJ.txt)
  echo "$1: G-if 登録 [$(grep VERDICT $r/CHT_INTERFACE_VERDICT.txt)] 旧 [$(grep VERDICT $r/CHT_INTERFACE_VERDICT_old_tol1e-9.txt)] | 主判定 [$(grep VERDICT $r/EVAL_CONJ.txt)] | conv [$(grep -o '\-> [A-Za-z ()]*' $r/CONVERGENCE_CHECK.txt | head -1)] | series [$(grep VERDICT $r/SERIES_CONJ.txt|tail -1)] | G-cons [$(grep -E 'VERDICT|REFUSED' $r/CHT_BALANCE_VERDICT.txt | tail -1)]"
done
