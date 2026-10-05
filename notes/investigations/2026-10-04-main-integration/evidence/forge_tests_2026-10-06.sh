#!/bin/bash
# forge 必須の Python 単体試験・二相 harness・潜熱 0 step A/B を AWS で回す (species-transport f01165a4、solver = 029dc631 と同一)
set +e
H=/home/ubuntu; O=/home/ubuntu/integ/vt; ST=/home/ubuntu/integ/vt/status.txt; mkdir -p /home/ubuntu/integ/vt; : > $ST
R=/home/ubuntu/forge-integ
cd $R && git fetch -q origin +refs/heads/feature/species-transport:refs/remotes/origin/feature/species-transport; git checkout -q --detach f01165a4
[ "$(git -C $R rev-parse --short=8 HEAD)" = "f01165a4" ] || { echo "WRONG HEAD" >> $ST; exit 1; }
export FORGE_BIN=/home/ubuntu/integ/v7e/forge_noarr FORGE_CUDA_BLOCKSIZE=128
cd $R/solver_density_cuda
G=/home/ubuntu/integ/vt/gen; mkdir -p $G
cmake -DIN=data/species/forge_species_v1.yaml -DOUT=$G/forge_species_data.hpp -P cmake/embed_species_data.cmake > /home/ubuntu/integ/vt/gen.log 2>&1
I="-I/usr/local/cuda/include -I/usr/local/cuda/include/cccl"
g++ -O2 -std=c++17 $I -I . -I $G -c input/speciesDB.cpp -o /home/ubuntu/integ/vt/speciesDB.o > /home/ubuntu/integ/vt/ab.build 2>&1 && g++ -O2 -std=c++17 $I -I . -I $G -c input/speciesTransportDB.cpp -o /home/ubuntu/integ/vt/speciesTransportDB.o >> /home/ubuntu/integ/vt/ab.build 2>&1 && nvcc -std=c++17 --expt-relaxed-constexpr -arch=sm_86 -I . -I $G -o /home/ubuntu/integ/vt/test_cond_latent_ab tests/unit/test_cond_latent_ab.cu /home/ubuntu/integ/vt/speciesDB.o /home/ubuntu/integ/vt/speciesTransportDB.o -lyaml-cpp >> /home/ubuntu/integ/vt/ab.build 2>&1 && /home/ubuntu/integ/vt/test_cond_latent_ab > /home/ubuntu/integ/vt/latent_ab.txt 2>&1
echo "latent A/B rc=$? $(tail -1 /home/ubuntu/integ/vt/latent_ab.txt)" >> $ST
for t in test_dmix_complement test_species_attrs_entry test_species_lump_solver test_species_record_solver test_thermo_intervals test_transport_gas_phase test_transport_gpu test_twophase_diffusion_harness; do
  timeout 3600 python3 tests/unit/$t.py --forge $FORGE_BIN > /home/ubuntu/integ/vt/$t.out 2>&1; rc=$?
  if grep -q "unrecognized arguments: --forge" /home/ubuntu/integ/vt/$t.out; then timeout 3600 python3 tests/unit/$t.py > /home/ubuntu/integ/vt/$t.out 2>&1; rc=$?; fi
  echo "$rc $t | $(tail -1 /home/ubuntu/integ/vt/$t.out | cut -c1-150)" >> $ST
done
echo DONE >> $ST
