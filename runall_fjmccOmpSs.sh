#!/bin/bash
#PJM -L elapse=07:00:00
#PJM -L node=1
#PJM -j
#PJM -o opt.out

module load fuji
export LANG=C
HOME_LOCAL=/ompss
export PATH=${HOME_LOCAL}/bin:$PATH
	
# compile
fjmcc -O3 --ompss -c -I${HOME_LOCAL}/include -DDESCPARAM=100 -DCOST_THRE=50000 -DUSE_NESTED_TASK -DUSE_MULTIDEPENDENCY cholmod_super_numeric.c -o "cholmod_super_numeric.o"
fjmcc -O3 --ompss -I${HOME_LOCAL}/include cholmod_demo.c cholmod_super_numeric.o -o cholmod_demo -L${HOME_LOCAL}/lib -lcholmod -lamd -lcolamd -lsuitesparseconfig -lccolamd -lcamd -lmetis -lm -lrt

export LD_LIBRARY_PATH=${HOME_LOCAL}/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=12
export XOS_MMM_L_HPAGE_TYPE=none
export NX_ARGS='--summary --stack-size=10000000'

iter=100
source groups.sh

for f in "${group1[@]}"
do
	for i in `seq 1 $iter`
	do
		./cholmod_demo "inputs/"$f".mtx"
        done
done
