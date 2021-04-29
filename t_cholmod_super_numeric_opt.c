/* ========================================================================== */
/* === Supernodal/t_cholmod_super_numeric =================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Supernodal Module.  Copyright (C) 2005-2012, Timothy A. Davis
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* Template routine for cholmod_super_numeric.  All xtypes supported, except
 * that a zomplex A and F result in a complex L (there is no supernodal
 * zomplex L).
 */

/* ========================================================================== */
/* === complex arithmetic =================================================== */
/* ========================================================================== */

#include "cholmod_template.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <sys/time.h>

#undef L_ENTRY
#undef L_CLEAR
#undef L_ASSIGN
#undef L_MULTADD
#undef L_ASSEMBLE
#undef L_ASSEMBLESUB

#ifdef REAL

/* -------------------------------------------------------------------------- */
/* A, F, and L are all real */
/* -------------------------------------------------------------------------- */

#define L_ENTRY 1
#define L_CLEAR(Lx,p)               Lx [p] = 0
#define L_ASSIGN(Lx,q, Ax,Az,p)     Lx [q] = Ax [p]
#define L_MULTADD(Lx,q, Ax,Az,p, f) Lx [q] += Ax [p] * f [0]
#define L_ASSEMBLE(Lx,q,b)          Lx [q] += b [0]
#define L_ASSEMBLESUB(Lx,q,C,p)     Lx [q] -= C [p]

#else

/* -------------------------------------------------------------------------- */
/* A and F are complex or zomplex, L and C are complex */
/* -------------------------------------------------------------------------- */

#define L_ENTRY 2
#define L_CLEAR(Lx,p)               Lx [2*(p)] = 0 ; Lx [2*(p)+1] = 0
#define L_ASSEMBLE(Lx,q,b)          Lx [2*(q)] += b [0] ;
#define L_ASSEMBLESUB(Lx,q,C,p)                 \
    Lx [2*(q)  ] -= C [2*(p)  ] ;               \
    Lx [2*(q)+1] -= C [2*(p)+1] ;

#ifdef COMPLEX

/* -------------------------------------------------------------------------- */
/* A, F, L, and C are all complex */
/* -------------------------------------------------------------------------- */

#define L_ASSIGN(Lx,q, Ax,Az,p)                 \
    Lx [2*(q)  ] = Ax [2*(p)  ] ;               \
    Lx [2*(q)+1] = Ax [2*(p)+1]

#define L_MULTADD(Lx,q, Ax,Az,p, f)                                     \
    Lx [2*(q)  ] += Ax [2*(p)  ] * f [0] - Ax [2*(p)+1] * f [1] ;       \
    Lx [2*(q)+1] += Ax [2*(p)+1] * f [0] + Ax [2*(p)  ] * f [1]

#else

/* -------------------------------------------------------------------------- */
/* A and F are zomplex, L and C is complex */
/* -------------------------------------------------------------------------- */

#define L_ASSIGN(Lx,q, Ax,Az,p)                 \
    Lx [2*(q)  ] = Ax [p] ;                     \
    Lx [2*(q)+1] = Az [p] ;

#define L_MULTADD(Lx,q, Ax,Az,p, f)                     \
    Lx [2*(q)  ] += Ax [p] * f [0] - Az [p] * f [1] ;   \
    Lx [2*(q)+1] += Az [p] * f [0] + Ax [p] * f [1]

#endif
#endif

/* ========================================================================== */
/* === t_cholmod_super_numeric ============================================== */
/* ========================================================================== */
/* This function returns FALSE only if integer overflow occurs in the BLAS.
 * It returns TRUE otherwise whether or not the matrix is positive definite. */

static int TEMPLATE (cholmod_super_numeric)
(
    /* ---- input ---- */
    cholmod_sparse *A,  /* matrix to factorize */
    cholmod_sparse *F,  /* F = A' or A(:,f)' */
    double beta [2],    /* beta*I is added to diagonal of matrix to factorize */
    /* ---- in/out --- */
    cholmod_factor *L,  /* factorization */
    /* -- workspace -- */
    cholmod_dense *Cwork,       /* size (L->maxcsize)-by-1 */
    /* --------------- */
    cholmod_common *Common
    )
{
    double one [2], zero [2], tstart ;
    double *Lx, *Ax, *Fx, *Az, *Fz, *C ;
    Int *Super, *Ls, *Lpi, *Lpx, *Map, *SuperMap, *RelativeMap,
        *Fp, *Fi, *Fnz, *Ap, *Ai, *Anz, *Iwork ;
    Int nsuper, n, j, i, k, s, p, pend, k1, k2, nscol, psi, psx, psend, nsrow,
        pj, d, kd1, kd2, info, ndcol, ndrow, pdi, pdx, pdend, pdi1, pdi2, pdx1,
        ndrow1, ndrow2, px, dnext, nsrow2, ndrow3, pk, pf,
        pfend, stype, Apacked, Fpacked, q, imap, nscol2, ss,
        tail;
    Int idxS, nnz_super, *Sp, *Si, *Sx, *STp, *STi, *STwork, *Sparent,
        *Liaison, csize;
    float cost;
    int maxChildren, optD;
    /* ---------------------------------------------------------------------- */
    /* guard against integer overflow in the BLAS */
    /* ---------------------------------------------------------------------- */

    /* If integer overflow occurs in the BLAS, Common->status is set to
     * CHOLMOD_TOO_LARGE, and the contents of Lx are undefined. */
    Common->blas_ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    nsuper = L->nsuper ;
    n = L->n ;

    C = Cwork->x ;      /* workspace of size L->maxcsize */

    one [0] =  1.0 ;    /* ALPHA for *syrk, *herk, *gemm, and *trsm */
    one [1] =  0. ;
    zero [0] = 0. ;     /* BETA for *syrk, *herk, and *gemm */
    zero [1] = 0. ;

    /* Iwork must be of size 2n + 4*nsuper +3, allocated in the caller,
     * cholmod_super_numeric.  The memory cannot be allocated here because the
     * cholmod_super_numeric initializes SuperMap, and cholmod_allocate_work
     * does not preserve existing workspace if the space needs to be increase
     * in size. */

    /* allocate integer workspace */
    Iwork = Common->Iwork ;
    SuperMap    = Iwork ;                                   /* size n (i/i/l) */
    RelativeMap = Iwork + n ;                               /* size n (i/i/l) */
    Sp       = Iwork + 2*((size_t) n) ;                        /* size nsuper+1*/
    STp      = Iwork + 2*((size_t) n) + nsuper+1 ;             /* size nsuper+1*/
    STwork   = Iwork + 2*((size_t) n) + 2*((size_t) nsuper+1) ;/* size nsuper+1*/
    Sparent  = Iwork + 2*((size_t) n) + 3*((size_t) nsuper+1) ;/* size nsuper*/

    Map  = Common->Flag ;   /* size n, use Flag as workspace for Map array */
    // Note: C, RelativeMap and Map above are not used, because thread-private
    //       areas are required.

    Ls = L->s ;
    Lpi = L->pi ;
    Lpx = L->px ;

    Super = L->super ;

    Lx = L->x ;

#ifndef NTIMER
    /* clear GPU / CPU statistics */
    Common->CHOLMOD_CPU_GEMM_CALLS  = 0 ;
    Common->CHOLMOD_CPU_SYRK_CALLS  = 0 ;
    Common->CHOLMOD_CPU_TRSM_CALLS  = 0 ;
    Common->CHOLMOD_CPU_POTRF_CALLS = 0 ;
    Common->CHOLMOD_GPU_GEMM_CALLS  = 0 ;
    Common->CHOLMOD_GPU_SYRK_CALLS  = 0 ;
    Common->CHOLMOD_GPU_TRSM_CALLS  = 0 ;
    Common->CHOLMOD_GPU_POTRF_CALLS = 0 ;
    Common->CHOLMOD_CPU_GEMM_TIME   = 0 ;
    Common->CHOLMOD_CPU_SYRK_TIME   = 0 ;
    Common->CHOLMOD_CPU_TRSM_TIME   = 0 ;
    Common->CHOLMOD_CPU_POTRF_TIME  = 0 ;
    Common->CHOLMOD_GPU_GEMM_TIME   = 0 ;
    Common->CHOLMOD_GPU_SYRK_TIME   = 0 ;
    Common->CHOLMOD_GPU_TRSM_TIME   = 0 ;
    Common->CHOLMOD_GPU_POTRF_TIME  = 0 ;
    Common->CHOLMOD_ASSEMBLE_TIME   = 0 ;
    Common->CHOLMOD_ASSEMBLE_TIME2  = 0 ;
#endif

    stype = A->stype ;

    if (stype != 0)
    {
        /* F not accessed */
        Fp = NULL ;
        Fi = NULL ;
        Fx = NULL ;
        Fz = NULL ;
        Fnz = NULL ;
        Fpacked = TRUE ;
    }
    else
    {
        Fp = F->p ;
        Fi = F->i ;
        Fx = F->x ;
        Fz = F->z ;
        Fnz = F->nz ;
        Fpacked = F->packed ;
    }

    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Az = A->z ;
    Anz = A->nz ;
    Apacked = A->packed ;

    /* clear the Map so that changes in the pattern of A can be detected */

    for (i = 0 ; i < n ; i++)
    {
        Map [i] = EMPTY ;
    }

    /* When the matrix is not positive definite, the supernode s containing
     * the first zero or negative diagonal entry of L can not be repeated
     * (because of asyncronous execution of tasks).
     * This can not provide MATLAB with [R,p]=chol(A) which require columns 1 to
     * p-1 of L=R', where L(p,p) is the problematic diagonal entry.
     * Thus the repeat_supernode flag is not used. */
//    repeat_supernode = FALSE ;
    L->minor = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* supernodal numerical factorization */
    /* ---------------------------------------------------------------------- */

// PROTOTYPE OF TASK-PARALLELIZATION written by Tetsuzo Usui
// Modified code is subject to CHOLMOD licence.
    // estimate nnz_super (table size)
    idxS = 0 ;
    for (s = 0 ; s < nsuper ; s++){
        Sp [s] = idxS ;
        psi = Lpi [s] ;
        psend = Lpi [s+1] ;
        d = SuperMap [Ls [psi]] ;
        idxS++ ;
        for (pdi2 = psi ; pdi2 < psend ; pdi2++ ){
            dnext = SuperMap [Ls [pdi2]] ;
            if ( dnext != d){
                d = dnext ;
                idxS++ ;
            }
        }
    }
    nnz_super = idxS ;
    Si = CHOLMOD(malloc) (nnz_super, sizeof (Int), Common) ;
    Sx = CHOLMOD(malloc) (nnz_super, sizeof (Int), Common) ;
    STi = CHOLMOD(malloc) (nnz_super, sizeof (Int), Common) ;
    // prepare offset table in Sp,Si,Sx
    idxS = 0 ;
    for (s = 0 ; s < nsuper ; s++){
        Sp [s] = idxS ;
        psi = Lpi [s] ;
        psend = Lpi [s+1] ;
        d = SuperMap [Ls [psi]] ;
        Si [idxS] = d ;
        Sx [idxS++] = 0 ;
        for (pdi2 = psi ; pdi2 < psend ; pdi2++ ){
            dnext = SuperMap [Ls [pdi2]] ;
            if ( dnext != d){
                d = dnext ;
                Si [idxS] = d ;
                Sx [idxS++] = pdi2 - psi ;
            }
        }
    }
    Sp [s] = idxS ;
    // relaxed supernodal etree (which is identical to Sparent in cholmod_super_symbolic.c)
    for (s = 0 ; s < nsuper ; s++) {
        if (Sp[s+1]-Sp[s] > 1 ) {
            Sparent[s] = Si[Sp[s]+1] ;
        }else{
            Sparent[s] = EMPTY ;
        }
    }

    // convert (Sp,Si) to (STp,STi), which will provide d-supernode-list for each task
    for (s = 0 ; s < nsuper ; s++) STwork[s] = 0;
    for (idxS = 0 ; idxS < Sp[nsuper] ; idxS++) STwork[Si[idxS]]++; // row counts
    idxS = 0 ;
    maxChildren = 0;
    for (s = 0 ; s < nsuper ; s++) {
        STp[s] = idxS ;
        idxS += STwork[s] ;
        STwork[s] = STp[s] ;
        if (s>0)
	    maxChildren = (maxChildren < STp[s]-STp[s-1] ? STp[s]-STp[s-1] : maxChildren);
    }
    STp[s] = idxS ;
    maxChildren = (maxChildren < STp[s]-STp[s-1] ? STp[s]-STp[s-1] : maxChildren);
    for (s = 0 ; s < nsuper ; s++) {
        for (p = Sp[s] ; p < Sp[s+1] ; p++) {
            STi[STwork[Si[p]]++] = s ;
        }
    }

   //We compute the optimal value of D
   int distribChildren[maxChildren];
   for (s = 0; s < maxChildren; s++)
	distribChildren[s] = 0;
   for (s = 0; s < nsuper; s++)
   {
	distribChildren[STp[s+1]-STp[s]-1]++;
   }
   int goalTasks=(int) ((double)n/14.0 < 1.1*(double)nsuper ? 1.1*(double)nsuper : (double)n/14.0);
   int currentTasks=nsuper;
   optD = maxChildren;
   int minOuter = 0.001*(double)nsuper;
   int nbOuter = 0;
   while ( (currentTasks < goalTasks || (double)optD/(double)maxChildren > 0.3 || nbOuter < minOuter)  && optD > 0)
   {
	optD--;
	currentTasks += (distribChildren[optD])*optD;
	nbOuter += distribChildren[optD];
   }
   optD++; //from index (0..maxChildren-1) to real value (1..maxChildren)
   printf("ALGO: goalTasks = %d, optD = %d\n",goalTasks,optD);

#ifdef USE_MULTIDEPENDENCY
    Int** dep_in[nsuper];
    Int* dep_out[nsuper];
    int num_in[nsuper];
    int size = 20;
    Liaison = CHOLMOD(malloc) (nsuper, sizeof (Int), Common); // dummy array for dependency-configuration
    for (s = 0 ; s < nsuper ; s++) {
        dep_in[s] = malloc(size * sizeof(Int*));
        int numIn = 0;
        for (int si = 0; si < nsuper; si++) {
            if (si == s) //Find the output dep
                dep_out[s] = &Liaison[s];
            else if (Sparent[si] == s) {
                if (numIn >= size) {
                    size += size;
                    dep_in[s] = realloc(dep_in[s], size * sizeof(Int*));
                    assert(dep_in[s] != NULL);
                }
                dep_in[s][numIn++] = &Liaison[si];
            }
        }
        num_in[s] = numIn;
    }
#else
    Int* dep_in[nsuper];
    Int* dep_out[nsuper];
    Liaison = CHOLMOD(malloc) (2*nsuper, sizeof (Int), Common); // dummy array for dependency-configuration
    for (s = 0 ; s < nsuper ; s++)
    {
      dep_in[s] = &Liaison[s];
      if (Sparent[s] != EMPTY) {
        dep_out[s] = &Liaison[Sparent[s]];
      }else{
        dep_out[s] = &Liaison[s+nsuper];
      }
    }
#endif
 
    /*int threshold = (int) nsuper * 0.9;
    int sum_sze = 0;
    for (int i = threshold; i < nsuper; ++i)
        sum_sze += STp[i+1] - STp[i];
    int thv = sum_sze / (nsuper - threshold);*/
    
    struct timeval  start, end;
    gettimeofday(&start, NULL);
#ifdef USE_OPENMP_TASK
#pragma omp parallel
  {
#pragma omp master
   {
#endif
    for (s = 0 ; s < nsuper ; s++)
    {

        k1 = Super [s] ;            /* s contains columns k1 to k2-1 of L */
        k2 = Super [s+1] ;
        nscol = k2 - k1 ;           /* # of columns in all of s */
        psi = Lpi [s] ;             /* pointer to first row of s in Ls */
        psx = Lpx [s] ;             /* pointer to first row of s in Lx */
        psend = Lpi [s+1] ;         /* pointer just past last row of s in Ls */
        nsrow = psend - psi ;       /* # of rows in all of s */
        int isFinal = 1;
        //printf("Stp: %d nsrow*nscol: %d\n", STp[s+1]-STp[s], nsrow*nscol);
        //printf("--------------------------------------------\n");
        if (STp[s+1] - STp[s] >= optD)
        //if ( nsrow * nscol >= 90000)
            isFinal = 0;
        /*printf("STAT Supernode %d with %d children: %d\n",s,STp[s+1]-STp[s],isFinal);
	if (STp[s+1]-STp[s] == 1)
	   printf("SN %d => %d\n",s,STi[STp[s]]);*/
        //isFinal = 1;
#ifdef USE_MULTIDEPENDENCY
#pragma omp task in({*(dep_in[s][ii]), ii=0;num_in[s]}) out(*dep_out[s]) \
                 default(none) shared(Ap,Ai,Ax,Az,Anz,Ls,Lx,Lpi,Lpx,stype,zero,one,Super, \
                 Sp,Si,Sx,STp,STi,Common,beta,Apacked,Fpacked,Fz,Fnz,Fi,Fx,Fp, \
                 Liaison, n, L, nsuper, Sparent, dep_in,num_in) firstprivate(s) \
                 private(C,Map, i,j,k,k1,k2, \
                 psi,psx,psend,pend,p,pk,px,q,nscol,nscol2,nsrow,nsrow2, \
                 d,kd1,kd2,ndcol,pdi1,pdi2,pdi,pdx,pdx1,pdend, \
                 ndrow,ndrow1,ndrow2,ndrow3,pf,imap,pfend,idxS,tstart,info, \
                 csize,cost) final(isFinal) label(outer)
#else
#ifdef USE_OPENMP_TASK
#pragma omp task depend(inout:dep_in[s][0]) depend(out:dep_out[s][0]) \
                 default(none) shared(Ap,Ai,Ax,Az,Anz,Ls,Lx,Lpi,Lpx,stype,zero,one,Super, \
                 Sp,Si,Sx,STp,STi,Common,beta,Apacked,Fpacked,Fz,Fnz,Fi,Fx,Fp, \
                 Liaison, n, L, nsuper, Sparent,dep_in,dep_out) firstprivate(s) \
                 private(C,Map, i,j,k,k1,k2, \
                 psi,psx,psend,pend,p,pk,px,q,nscol,nscol2,nsrow,nsrow2, \
                 d,kd1,kd2,ndcol,pdi1,pdi2,pdi,pdx,pdx1,pdend, \
                 ndrow,ndrow1,ndrow2,ndrow3,pf,imap,pfend,idxS,tstart,info, \
                 csize,cost) final(isFinal) 
#else
#pragma omp task in(*dep_in[s]) out(*dep_out[s]) \
                 default(none) shared(Ap,Ai,Ax,Az,Anz,Ls,Lx,Lpi,Lpx,stype,zero,one,Super, \
                 Sp,Si,Sx,STp,STi,Common,beta,Apacked,Fpacked,Fz,Fnz,Fi,Fx,Fp, \
                 Liaison, n, L, nsuper, Sparent) firstprivate(s) \
                 private(C,Map, i,j,k,k1,k2, \
                 psi,psx,psend,pend,p,pk,px,q,nscol,nscol2,nsrow,nsrow2, \
                 d,kd1,kd2,ndcol,pdi1,pdi2,pdi,pdx,pdx1,pdend, \
                 ndrow,ndrow1,ndrow2,ndrow3,pf,imap,pfend,idxS,tstart,info, \
                 csize,cost) label(outer)
#endif
#endif
      {
#ifdef USE_NESTED_TASK
        // prepare a private lock for each outer task, which will be shared for inner tasks
        omp_lock_t omp_lock;
        omp_init_lock(&omp_lock);
#endif
        // allocate private work array for each task in heap area
        Map = CHOLMOD(malloc) (n, sizeof (Int), Common) ;
        //for (i = 0 ; i < n ; i++) { Map [i] = EMPTY ; }

        /* ------------------------------------------------------------------ */
        /* get the size of supernode s */
        /* ------------------------------------------------------------------ */
        k1 = Super [s] ;            /* s contains columns k1 to k2-1 of L */
        k2 = Super [s+1] ;
        nscol = k2 - k1 ;           /* # of columns in all of s */
        psi = Lpi [s] ;             /* pointer to first row of s in Ls */
        psx = Lpx [s] ;             /* pointer to first row of s in Lx */
        psend = Lpi [s+1] ;         /* pointer just past last row of s in Ls */
        nsrow = psend - psi ;       /* # of rows in all of s */

        PRINT1 (("====================================================\n"
                 "S "ID" k1 "ID" k2 "ID" nsrow "ID" nscol "ID" psi "ID" psend "
                 ""ID" psx "ID"\n", s, k1, k2, nsrow, nscol, psi, psend, psx)) ;

        /* ------------------------------------------------------------------ */
        /* zero the supernode s */
        /* ------------------------------------------------------------------ */
        ASSERT ((size_t) (psx + nsrow*nscol) <= L->xsize) ;
        pend = psx + nsrow * nscol ;        /* s is nsrow-by-nscol */
        {
            for (p = psx ; p < pend ; p++) {
                L_CLEAR (Lx,p);
            }
        }

        /* ------------------------------------------------------------------ */
        /* construct the scattered Map for supernode s */
        /* ------------------------------------------------------------------ */

        /* If row i is the kth row in s, then Map [i] = k.  Similarly, if
         * column j is the kth column in s, then  Map [j] = k. */

        for (k = 0 ; k < nsrow ; k++)
        {
            PRINT1 (("  "ID" map "ID"\n", Ls [psi+k], k)) ;
            Map [Ls [psi + k]] = k ;
        }

        /* ------------------------------------------------------------------ */
        /* copy matrix into supernode s (lower triangular part only) */
        /* ------------------------------------------------------------------ */

        pk = psx ;

        for (k = k1 ; k < k2 ; k++)
        {
            if (stype != 0)
            {
                /* copy the kth column of A into the supernode */
                p = Ap [k] ;
                pend = (Apacked) ? (Ap [k+1]) : (p + Anz [k]) ;
                for ( ; p < pend ; p++)
                {
                    /* row i of L is located in row Map [i] of s */
                    i = Ai [p] ;
                    if (i >= k)
                    {
                        /* This test is here simply to avoid a segfault.  If
                         * the test is false, the numeric factorization of A
                         * is undefined.  It does not detect all invalid
                         * entries, only some of them (when debugging is
                         * enabled, and Map is cleared after each step, then
                         * all entries not in the pattern of L are detected). */
                        imap = Map [i] ;
                        if (imap >= 0 && imap < nsrow)
                        {
                            /* Lx [Map [i] + pk] = Ax [p] ; */
                            L_ASSIGN (Lx,(imap+(psx+(k-k1)*nsrow)), Ax,Az,p) ;
                        }
                    }
                }
            }
            else
            {
                double fjk[2];
                /* copy the kth column of A*F into the supernode */
                pf = Fp [k] ;
                pfend = (Fpacked) ? (Fp [k+1]) : (p + Fnz [k]) ;
                for ( ; pf < pfend ; pf++)
                {
                    j = Fi [pf] ;

                    /* fjk = Fx [pf] ; */
                    L_ASSIGN (fjk,0, Fx,Fz,pf) ;

                    p = Ap [j] ;
                    pend = (Apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
                    for ( ; p < pend ; p++)
                    {
                        i = Ai [p] ;
                        if (i >= k)
                        {
                            /* See the discussion of imap above. */
                            imap = Map [i] ;
                            if (imap >= 0 && imap < nsrow)
                            {
                                /* Lx [Map [i] + pk] += Ax [p] * fjk ; */
                                L_MULTADD (Lx,(imap+(psx+(k-k1)*nsrow)),
                                           Ax,Az,p, fjk) ;
                            }
                        }
                    }
                }
            }
        }

        /* add beta to the diagonal of the supernode, if nonzero */
        if (beta [0] != 0.0)
        {
            /* note that only the real part of beta is used */
            pk = psx ;
            for (k = k1 ; k < k2 ; k++)
            {
                /* Lx [pk] += beta [0] ; */
                L_ASSEMBLE (Lx,pk, beta) ;
                pk += nsrow + 1 ;       /* advance to the next diagonal entry */
            }
        }

        PRINT1 (("Supernode with just A: repeat: "ID"\n", repeat_supernode)) ;
        DEBUG (CHOLMOD(dump_super) (s, Super, Lpi, Ls, Lpx, Lx, L_ENTRY,
                                    Common)) ;
        PRINT1 (("\n\n")) ;

        /* ------------------------------------------------------------------ */
        /* update supernode s with each pending descendant d */
        /* ------------------------------------------------------------------ */

#ifndef NDEBUG
        for (idxS = STp [s] ; idxS < STp[s+1] ; idxS++)
        {
            d = STi[idxS] ;
            PRINT1 (("\nWill update "ID" with Child: "ID"\n", s, d)) ;
            DEBUG (CHOLMOD(dump_super) (d, Super, Lpi, Ls, Lpx, Lx, L_ENTRY,
                                        Common)) ;
        }
        PRINT1 (("\nNow factorizing supernode "ID":\n", s)) ;
#endif

        for (idxS = STp [s] ; idxS < STp[s+1] ; idxS++) {
            d = STi[idxS] ;
            if (d==s) continue ;

            /* -------------------------------------------------------------- */
            /* get the size of supernode d */
            /* -------------------------------------------------------------- */

            kd1 = Super [d] ;       /* d contains cols kd1 to kd2-1 of L */
            kd2 = Super [d+1] ;
            ndcol = kd2 - kd1 ;     /* # of columns in all of d */
            pdi = Lpi [d] ;         /* pointer to first row of d in Ls */
            pdx = Lpx [d] ;         /* pointer to first row of d in Lx */
            pdend = Lpi [d+1] ;     /* pointer just past last row of d in Ls */
            ndrow = pdend - pdi ;   /* # rows in all of d */

            PRINT1 (("Child: ")) ;
            DEBUG (CHOLMOD(dump_super) (d, Super, Lpi, Ls, Lpx, Lx, L_ENTRY,
                                        Common)) ;

            /* -------------------------------------------------------------- */
            /* find the range of rows of d that affect rows k1 to k2-1 of s */
            /* -------------------------------------------------------------- */

            for (p=Sp[d] ; Si[p]<s && p<Sp[d+1] ; p++);
            p = Sx[p] ;             /* offset of 1st row of d affecting s */
            pdi1 = pdi + p ;        /* ptr to 1st row of d affecting s in Ls */
            pdx1 = pdx + p ;        /* ptr to 1st row of d affecting s in Lx */

            /* there must be at least one row remaining in d to update s */
            ASSERT (pdi1 < pdend) ;
            ASSERT (Ls [pdi1] >= k1 && Ls [pdi1] < k2) ;

            for (pdi2 = pdi1 ; pdi2 < pdend && Ls [pdi2] < k2 ; pdi2++) ;
            // pdi2 is Sx[p+1], provided that p was preserved such as Si[p]=s and p+1<Sp[d+1]
            ndrow1 = pdi2 - pdi1 ;      /* # rows in first part of d */
            ndrow2 = pdend - pdi1 ;     /* # rows in remaining d */

            /* rows Ls [pdi1 ... pdi2-1] are in the range k1 to k2-1.  Since d
             * affects s, this set cannot be empty. */
            ASSERT (pdi1 < pdi2 && pdi2 <= pdend) ;
            PRINT1 (("ndrow1 "ID" ndrow2 "ID"\n", ndrow1, ndrow2)) ;
            DEBUG (for (p = pdi1 ; p < pdi2 ; p++)
                       PRINT1 (("Ls["ID"] "ID"\n", p, Ls[p]))) ;
	   
            ndrow3 = ndrow2 - ndrow1 ;  /* number of rows of C2 */
            ASSERT (ndrow3 >= 0) ;

            cost = (float)ndrow1*ndrow1*ndcol*2/2; // dsyrk
            cost += (float)ndrow3*ndrow1*ndcol*2;  // dgemm
            const int cost_thre = COST_THRE;
#ifdef USE_NESTED_TASK
#ifdef USE_OPENMP_TASK
#pragma omp task default(none) shared(Common,Super,Lpi,Lpx,zero,one,Lx, \
                 n,Sp,Si,Sx,Ls, omp_lock) \
                 firstprivate(d,s,k2,Map,psx,nsrow, \
                 kd1,kd2,ndcol,pdi,pdx,pdend,ndrow, pdi1,pdx1,ndrow1,ndrow2,ndrow3) \
		 private(C, tstart,i,j,px,q, csize)
#else
#pragma omp task default(none) shared(Common,Super,Lpi,Lpx,zero,one,Lx, \
                 n,Sp,Si,Sx,Ls, omp_lock) \
                 firstprivate(d,s,k2,Map,psx,nsrow, \
                 kd1,kd2,ndcol,pdi,pdx,pdend,ndrow,pdi1,pdx1,ndrow1,ndrow2,ndrow3) \
                 private(C, tstart,i,j,px,q, csize) label(inner)
#endif
#endif
          {
            /* -------------------------------------------------------------- */
            /* construct the update matrix C for this supernode d */
            /* -------------------------------------------------------------- */

            /* C = L (k1:n-1, kd1:kd2-1) * L (k1:k2-1, kd1:kd2-1)', except
             * that k1:n-1 refers to all of the rows in L, but many of the
             * rows are all zero.  Supernode d holds columns kd1 to kd2-1 of L.
             * Nonzero rows in the range k1:k2-1 are in the list
             * Ls [pdi1 ... pdi2-1], of size ndrow1.  Nonzero rows in the range
             * k2:n-1 are in the list Ls [pdi2 ... pdend], of size ndrow2.  Let
             * L1 = L (Ls [pdi1 ... pdi2-1], kd1:kd2-1), and let
             * L2 = L (Ls [pdi2 ... pdend],  kd1:kd2-1).  C is ndrow2-by-ndrow1.
             * Let C1 be the first ndrow1 rows of C and let C2 be the last
             * ndrow2-ndrow1 rows of C.  Only the lower triangular part of C1
             * needs to be computed since C1 is symmetric.
             */

            // allocate private work arrays for each task in heap area
            csize = ndrow2 * ndrow1 ;
            C = CHOLMOD(malloc) (csize, L_ENTRY*sizeof (double), Common) ;
            
            /* maxcsize is the largest size of C for all pairs (d,s) */
            ASSERT (ndrow2 * ndrow1 <= ((Int) L->maxcsize)) ;

            /* compute leading ndrow1-by-ndrow1 lower triangular block of C,
             * C1 = L1*L1' */

            {
#ifndef NTIMER

                Common->CHOLMOD_CPU_SYRK_CALLS++ ;
                tstart = SuiteSparse_time () ;
#endif
#ifdef REAL
                BLAS_dsyrk ("L", "N",
                    ndrow1, ndcol,              /* N, K: L1 is ndrow1-by-ndcol*/
                    one,                        /* ALPHA:  1 */
                    Lx + L_ENTRY*pdx1, ndrow,   /* A, LDA: L1, ndrow */
                    zero,                       /* BETA:   0 */
                    C, ndrow2) ;                /* C, LDC: C1 */
#else
                BLAS_zherk ("L", "N",
                    ndrow1, ndcol,              /* N, K: L1 is ndrow1-by-ndcol*/
                    one,                        /* ALPHA:  1 */
                    Lx + L_ENTRY*pdx1, ndrow,   /* A, LDA: L1, ndrow */
                    zero,                       /* BETA:   0 */
                    C, ndrow2) ;                /* C, LDC: C1 */
#endif
#ifndef NTIMER
                Common->CHOLMOD_CPU_SYRK_TIME += SuiteSparse_time () - tstart ;
#endif
                /* compute remaining (ndrow2-ndrow1)-by-ndrow1 block of C,
                 * C2 = L2*L1' */
                if (ndrow3 > 0)
                {
#ifndef NTIMER
                    Common->CHOLMOD_CPU_GEMM_CALLS++ ;
                    tstart = SuiteSparse_time () ;
#endif
#ifdef REAL
                    BLAS_dgemm ("N", "C",
                        ndrow3, ndrow1, ndcol,          /* M, N, K */
                        one,                            /* ALPHA:  1 */
                        Lx + L_ENTRY*(pdx1 + ndrow1),   /* A, LDA: L2 */
                        ndrow,                          /* ndrow */
                        Lx + L_ENTRY*pdx1,              /* B, LDB: L1 */
                        ndrow,                          /* ndrow */
                        zero,                           /* BETA:   0 */
                        C + L_ENTRY*ndrow1,             /* C, LDC: C2 */
                        ndrow2) ;
#else
                    BLAS_zgemm ("N", "C",
                        ndrow3, ndrow1, ndcol,          /* M, N, K */
                        one,                            /* ALPHA:  1 */
                        Lx + L_ENTRY*(pdx1 + ndrow1),   /* A, LDA: L2 */
                        ndrow,                          /* ndrow */
                        Lx + L_ENTRY*pdx1,              /* B, LDB: L1, ndrow */
                        ndrow,
                        zero,                           /* BETA:   0 */
                        C + L_ENTRY*ndrow1,             /* C, LDC: C2 */
                        ndrow2) ;
#endif
#ifndef NTIMER
                    Common->CHOLMOD_CPU_GEMM_TIME +=
                        SuiteSparse_time () - tstart ;
#endif
                }

                /* ---------------------------------------------------------- */
                /* assemble C into supernode s using the relative map */
                /* ---------------------------------------------------------- */

#ifdef USE_NESTED_TASK
#ifdef USE_CRITICAL_SECTION
#pragma omp critical
#else
               omp_set_lock(&omp_lock);
#endif
#endif
               {
                for (j = 0 ; j < ndrow1 ; j++)              /* cols k1:k2-1 */
                {
                    px = psx + Map [Ls [pdi1 + j]] * nsrow ;
                    for (i = j ; i < ndrow2 ; i++)          /* rows k1:n-1 */
                    {
                        //fprintf(stderr, "start: j %d i %d\n", j, i);
                        /* Lx [px + RelativeMap [i]] -= C [i + pj] ; */
                        q = px + Map [Ls [pdi1 + i]] ;
                        //fprintf(stderr, "assemble: q %d i+ndrow2*j %d\n", q, i+ndrow2*j);
                        //fprintf(stderr, "assemble: Lx[%d] %d C[%d] %d\n", q, Lx[q], i+ndrow2*j, C[i+ndrow2*j]);
                        L_ASSEMBLESUB (Lx,q, C, i+ndrow2*j) ;
                        //fprintf(stderr, "end: j %d i %d\n", j, i);
                    }
                }
               } // end critical
#ifdef USE_NESTED_TASK
#ifndef USE_CRITICAL_SECTION
               omp_unset_lock(&omp_lock);
#endif
#endif
            }
            C = CHOLMOD(free) (csize, L_ENTRY*sizeof (double), C, Common) ;
          } // end inner task
        }  /* end of descendant supernode loop */
#ifdef USE_NESTED_TASK
#pragma omp taskwait
#endif

        PRINT1 (("\nSupernode with contributions A: repeat: "ID"\n",
                 repeat_supernode)) ;
       DEBUG (CHOLMOD(dump_super) (s, Super, Lpi, Ls, Lpx, Lx, L_ENTRY,
                                    Common)) ;
        PRINT1 (("\n\n")) ;

        /* ------------------------------------------------------------------ */
        /* factorize diagonal block of supernode s in LL' */
        /* ------------------------------------------------------------------ */

        /* The current supernode s is ready to factorize.  It has been updated
         * by all descendant supernodes.  Let S = the current supernode, which
         * holds rows k1:n-1 and columns k1:k2-1 of the updated matrix.   It
         * splits into two parts:  the square diagonal block S1, and the
         * rectangular part S2.  Here, S1 is factorized into L1*L1' and
         * overwritten by L1.
         *
         * If supernode s is being repeated, only factorize it up to but not
         * including the column containing the problematic entry.
         */

        nscol2 = nscol ;

        {
#ifndef NTIMER
            Common->CHOLMOD_CPU_POTRF_CALLS++ ;
            tstart = SuiteSparse_time () ;
#endif
#ifdef REAL
            LAPACK_dpotrf ("L",
                nscol2,                     /* N: nscol2 */
                Lx + L_ENTRY*psx, nsrow,    /* A, LDA: S1, nsrow */
                info) ;                     /* INFO */
#else
            LAPACK_zpotrf ("L",
                nscol2,                     /* N: nscol2 */
                Lx + L_ENTRY*psx, nsrow,    /* A, LDA: S1, nsrow */
                info) ;                     /* INFO */
#endif
#ifndef NTIMER
            Common->CHOLMOD_CPU_POTRF_TIME += SuiteSparse_time ()- tstart ;
#endif
        }

        /* ------------------------------------------------------------------ */
        /* check if the matrix is not positive definite */
        /* ------------------------------------------------------------------ */

        /* info is set to one in LAPACK_*potrf if blas_ok is FALSE.  It is
         * set to zero in dpotrf/zpotrf if the factorization was successful. */
        if (CHECK_BLAS_INT && !Common->blas_ok)
        {
            ERROR (CHOLMOD_TOO_LARGE, "problem too large for the BLAS") ;
        }

        if (info != 0)
        {
            /* Matrix is not positive definite.  dpotrf/zpotrf do NOT report an
             * error if the diagonal of L has NaN's, only if it has a zero. */
            if (Common->status == CHOLMOD_OK)
            {
                ERROR (CHOLMOD_NOT_POSDEF, "matrix not positive definite") ;
            }

            /* L->minor is the column of L that contains a zero or negative
             * diagonal term. */
            if(L->minor == EMPTY) L->minor = k1 + info - 1 ;
        }

        /* ------------------------------------------------------------------ */
        /* compute the subdiagonal block and prepare supernode for its parent */
        /* ------------------------------------------------------------------ */

        nsrow2 = nsrow - nscol2 ;
        if (nsrow2 > 0)
        {
            /* The current supernode is columns k1 to k2-1 of L.  Let L1 be the
             * diagonal block (factorized by dpotrf/zpotrf above; rows/cols
             * k1:k2-1), and L2 be rows k2:n-1 and columns k1:k2-1 of L.  The
             * triangular system to solve is L2*L1' = S2, where S2 is
             * overwritten with L2.  More precisely, L2 = S2 / L1' in MATLAB
             * notation.
             */

            {
#ifndef NTIMER
                Common->CHOLMOD_CPU_TRSM_CALLS++ ;
                tstart = SuiteSparse_time () ;
#endif
#ifdef REAL
                BLAS_dtrsm ("R", "L", "C", "N",
                    nsrow2, nscol2,                 /* M, N */
                    one,                            /* ALPHA: 1 */
                    Lx + L_ENTRY*psx, nsrow,        /* A, LDA: L1, nsrow */
                    Lx + L_ENTRY*(psx + nscol2),    /* B, LDB, L2, nsrow */
                    nsrow) ;
#else
                BLAS_ztrsm ("R", "L", "C", "N",
                    nsrow2, nscol2,                 /* M, N */
                    one,                            /* ALPHA: 1 */
                    Lx + L_ENTRY*psx, nsrow,        /* A, LDA: L1, nsrow */
                    Lx + L_ENTRY*(psx + nscol2),    /* B, LDB, L2, nsrow */
                    nsrow) ;
#endif
#ifndef NTIMER
                Common->CHOLMOD_CPU_TRSM_TIME += SuiteSparse_time () - tstart ;
#endif
            }

            if (CHECK_BLAS_INT && !Common->blas_ok)
            {
                ERROR (CHOLMOD_TOO_LARGE, "problem too large for the BLAS") ;
            }

        }

        /* clear the Map (debugging only, to detect changes in pattern of A) */
        DEBUG (for (k = 0 ; k < nsrow ; k++) Map [Ls [psi + k]] = EMPTY) ;
        DEBUG (CHOLMOD(dump_super) (s, Super, Lpi, Ls, Lpx, Lx, L_ENTRY,
                                    Common)) ;
        Map = CHOLMOD(free) (n, sizeof (Int), Map, Common) ;
#ifdef USE_NESTED_TASK
        omp_destroy_lock(&omp_lock);
#endif
      } // end outer task
    } // end for(s)
#pragma omp taskwait
#ifdef USE_OPENMP_TASK
   } // end master
  } // end parallel
#endif
    gettimeofday(&end, NULL);
    double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    printf(">>>> Time::::: %lf ms\n", delta);

    Si = CHOLMOD(free) (nnz_super, sizeof (Int), Si, Common) ;
    Sx = CHOLMOD(free) (nnz_super, sizeof (Int), Sx, Common) ;
    STi = CHOLMOD(free) (nnz_super, sizeof (Int), STi, Common) ;

#ifdef USE_MULTIDEPENDENCY
    for (s = 0 ; s < nsuper ; s++) {
        free(dep_in[s]);
        dep_in[s] = NULL;
    }
    Liaison = CHOLMOD(free) (nsuper, sizeof (Int), Liaison, Common) ;
#else
    Liaison = CHOLMOD(free) (2*nsuper, sizeof (Int), Liaison, Common) ;
#endif

    if( L->minor == EMPTY) {
      /* success; matrix is positive definite */
      L->minor = n ;
    }
    return (Common->status >= CHOLMOD_OK) ;
}

#undef PATTERN
#undef REAL
#undef COMPLEX
#undef ZOMPLEX
