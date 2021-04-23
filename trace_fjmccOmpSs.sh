#!/bin/bash
#PJM -L node=1
#PJM -L elapse=3600
#PJM -o testing.out
#PJM -e testing.err
#PJM -s

module load fuji
module load EXTRAE

export LANG=C
HOME_LOCAL=/fefs/scratch/bsc28/bsc28151/CHOLMOD-Tetsuzo/ompss/ompss
export PATH=${HOME_LOCAL}/bin:$PATH
# compile
#fjmcc -O3 --ompss --instrument -c -I${HOME_LOCAL}/include -DUSE_NESTED_TASK -DUSE_MULTIDEPENDENCY cholmod_super_numeric_nested.c
fjmcc -O3 --ompss --instrument -c -I${HOME_LOCAL}/include -DDESCPARAM=1 -DCOST_THRE=0 -DUSE_NESTED_TASK -DUSE_MULTIDEPENDENCY cholmod_super_numeric.c
fjmcc -O3 --ompss --instrument -I${HOME_LOCAL}/include cholmod_demo.c cholmod_super_numeric.o -o cholmod_demo-extrae -L${HOME_LOCAL}/lib -lcholmod -lamd -lcolamd -lsuitesparseconfig -lccolamd -lcamd -lmetis -lm -lrt

#export LD_LIBRARY_PATH=${HOME_LOCAL}/lib:home/local/gnu/gcc-7.3.0/lib64:/opt/intel/compilers_and_libraries_2020.2.254/linux/mkl/lib/intel64_lin:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${HOME_LOCAL}/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=12
export XOS_MMM_L_HPAGE_TYPE=none

#export EXTRAE_CONFIG_FILE=extrae.xml
export EXTRAE_HOME=/fefs/scratch/bsc28/bsc28151/CHOLMOD-Tetsuzo/ompss/ompss
export EXTRAE_OMP_COUNTERS_ON=1
export EXTRAE_COUNTERS=PAPI_TOT_INS,PAPI_TOT_CYC,PAPI_DP_OPS,PAPI_L2_DCM,PAPI_L1_DCM
export EXTRAE_SKIP_AUTO_LIBRARY_INITIALIZE=1
export EXTRAE_ON=1
export EXTRAE_PROGRAM_NAME=testing
export NX_ARGS="--summary --instrumentation=extrae --stack-size=10000000"

./cholmod_demo-extrae ../inputs/plbuckle.mtx
