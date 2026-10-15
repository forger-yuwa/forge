#!/bin/bash
# p2_chain.sh (準備) と p2_follow.sh (NS) を続けて回す (間を空けない)。AWS A。
D=$(cd "$(dirname "$0")" && pwd)
bash "$D/p2_chain.sh"; rc=$?
[ $rc -eq 0 ] && bash "$D/p2_follow.sh"
