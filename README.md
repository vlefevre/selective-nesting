# selective-nesting
Additional files (source code) for paper "A Selective Nesting Approach for the Sparse Multi-threaded Cholesky Factorization", Valentin Le Fèvre, Tetsuzo Usui and Marc Casas, submitted to CLUSTER 2021.

runall_fjmccOmpSs.sh is an example of a script file that compiles the cholmod_super_numeric.c and cholmod_demo.c files, linking with BLAS/SuiteSparse, for A64FX. You may need to adapt the directories to point towards your installation folders/librairies. Note that by default we use fjmcc, a mercurium compiler (for OmpSs) built using fcc (Fujitsu's compiler).  
To change the implementation file of the numerical factorization, lines 86,88 and 90 need to be changed in cholmod_super_numeric.c. When using the "general" implementation, threshold D can be set during compilation by using -DDESCPARAM=D (by default 100) and the threshold C on the cost of inner tasks can be set using -DTHRE_COST=C (by default 50,000).  
The script also runs the demo program on several matrices. By default, it runs on a group of matrices (Group 1 to 4 as defined in the paper) and the demo program is run iter times for each input matrix.

We also provide some raw data used in the paper in folder "raw_results". Figures 6 to 9 can be reproduced using this data and the python script plot_heuristics_opt_eval.py. The syntax is the following:  
python3 plot_heuristics_opt_eval.py group1/opt_12_g1  
which will generate Figure 1, for matrices of Group 1. Note that the opt_12_g1 file was generated using the script file raw_results/group1/out2table.sh.
The syntax for other groups is similar by changing group1 to groupX and g1 to gX, with X=2,3,4.

If you want to conduct similar experiments, please use the provided out2table.sh and adapt them to the number of iterations run and so on. The structure of the file opt_12_g1 is:  
1 line for matrix names (header)  
1 line for each heuristic that gives the average for each input matrix in the following order:  
Opt-D (or Opt-D-Cost)
Desc-100  
Desc-300  
Desc-500  
Desc-700  
SN-50k  
SN-70k  
SN-90k  
Non-nested  
Nested  
