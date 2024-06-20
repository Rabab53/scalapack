/*
    j=j+1
 * Copyright (c) 2009-2010 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2010      University of Denver, Colorado.
 */

#ifndef MYSCALAPACK_H
#define MYSCALAPACK_H
#ifdef __cplusplus

void start_timing(int myrank, const char *name, double *time);
void stop_timing(int myrank, const char *name, double *time);

int checking_under_determined(int mpirank, int myrow, int mycol, float *Acpy, int *descAA, float *Bcpy, int *descBT, 
                              float *BT, int MF, int M, int MP, int NF, int RHSF, int nb, int procs, int prow, int pcol, 
                              float eps, double elapsed);

int l2_under_determined(int mpirank, int myrow, int mycol, float *A, int *descA, float *B, int *descB, float *AA, int *descAA, int myAArows,
                float *BT, int *descBT, float *X, int *descX, float *Acpy, float *Bcpy, int MF, int M, int MP, int NF, int RHSF,
                float lambda, const char* kernal, int prow, int pcol, int nb, int check);

int checking_over_determined(int mpirank, int myrow, int mycol, float *Acpy, int *descAA, float *Bcpy, int *descX,
                float *X, int MF, int M, int MP, int NF, int RHSF, int nb, int procs, int prow, int pcol,
                float eps, double elapsed);

int l2_over_determined(int mpirank, int myrow, int mycol, float *A, int *descA, float *B, int *descB, float *AA, int *descAA, int myAArows,
                float *BT, int *descBT, float *X, int *descX, float *Acpy, float *Bcpy, int MF, int M, int MP, int NF, int RHSF,
                float lambda, const char* kernal, int prow, int pcol, int nb, int check);

extern "C" {
#endif
#ifdef SUP_
#define pdpotrf_  pdpotrf
#define pdpotri_  pdpotri
#define pdpotrs_  pdpotrs
#define pdsymm_   pdsymm
#define pdsyr2k_  pdsyr2k
#define pdgeadd_  pdgeadd
#define pdgetrf_  pdgetrf 
#define pdlansy_  pdlansy 
#define pdmatgen_ pdmatgen
#define pdtrsm_   pdtrsm
#define psgesv_   psgesv
#define pdgesv_   pdgesv
#define psgemm_   psgemm
#define pdgemm_   pdgemm
#define numroc_   numroc
#define pslange_  pslange
#define pdlange_  pdlange
#define pslacpy_  pslacpy
#define pdlacpy_  pdlacpy
#define pdgeqrf_  pdgeqrf
#define pdormqr_  pdormqr
#define psgesvd_  psgesvd
#define pdgesvd_  pdgesvd
#define pslapiv_  pslapiv
#define pslaset_  pslaset
#define pdlaset_  pdlaset
#define pselset_  pselset
#define pdelset_  pdelset
#define pslamch_  pslamch
#define pdlamch_  pdlamch
#define indxg2p_  indxg2p
#define indxg2l_  indxg2l
#define descinit_ descinit
#define pslawrite_ pslawrite
#define pdlawrite_ pdlawrite
#define blacs_get_      blacs_get
#define blacs_pinfo_    blacs_pinfo
#define blacs_gridmap_  blacs_gridmap
#define blacs_gridinit_ blacs_gridinit
#define blacs_gridinfo_ blacs_gridinfo
#define blacs_gridexit_ blacs_gridexit
#define blacs_exit_     blacs_exit
#endif

extern void Cblacs_pinfo( int* mypnum, int* nprocs);
extern void Cblacs_get( int context, int request, int* value);
extern int  Cblacs_gridinit( int* context, const char* order, int np_row, int np_col);
extern void Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
extern void Cblacs_gridmap( int *context, int*  imap, int lda, int np_row, int np_col);
extern void Cblacs_gridexit( int context);
extern void Cblacs_exit( int error_code);
extern void blacs_barrier( int *context, const char* *all);

extern void blacs_pinfo_( int *mypnum, int *nprocs);
extern void blacs_get_( int *context, int *request, int* value);
extern void blacs_gridinit_( int* context, const char *order, int *np_row, int *np_col);
extern void blacs_gridinfo_( int *context, int *np_row, int *np_col, int *my_row, int *my_col);
extern void blacs_gridexit_( int *context);
extern void blacs_exit_( int *error_code);

extern void pdgeqrf_( int *m, int *n, double *a, int *ia, int *ja, int *desca, double *tau, double *work, int *lwork, int *info );
extern void pdormqr_( const char *side, const char *trans, int *m, int *n, int *k, double *a, int *ia,
		int *ja, int *desca, double *tau, double *c, int *ic, int *jc, int *descc, double *work, int *lwork, int *info );
extern void pdtrsm_ ( const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, double *alpha, double *a, int *ia,
		int *ja, int *desca, double *b, int *ib, int *jb, int *descb );

extern float  pslange_( const char* norm, int *m, int *n, float     *A, int *ia, int *ja, int *descA, float *work);
extern double pdlange_( const char* norm, int *m, int *n, double    *A, int *ia, int *ja, int *descA, double *work);


extern void pslacpy_( const char* uplo, int *m, int *n, float     *A, int *ia, int *ja, int *descA,
                                                    float     *B, int *ib, int *jb, int *descB);
extern void pdlacpy_( const char* uplo, int *m, int *n, double     *A, int *ia, int *ja, int *descA,
                                                    double     *B, int *ib, int *jb, int *descB);

extern void psgesv_( int *n, int *nrhs, float     *A, int *ia, int *ja, int *descA, int *ipiv,
                                          float     *B, int *ib, int *jb, int *descB, int *info);
extern void pdgesv_( int *n, int *nrhs, double    *A, int *ia, int *ja, int *descA, int *ipiv,
                                          double    *B, int *ib, int *jb, int *descB, int *info);

extern void psgemm_( const char *transa, const char *transb, int *M, int *N, int *K,
                                          float     *alpha,
                                          float     *A, int *ia, int *ja, int *descA,
                                          float     *B, int *ib, int *jb, int *descB,
                                          float     *beta,
                                          float     *C, int *ic, int *jc, int *descC );
extern void pdgemm_( const char *transa, const char *transb, int *M, int *N, int *K,
                                          double    *alpha,
                                          double    *A, int *ia, int *ja, int *descA,
                                          double    *B, int *ib, int *jb, int *descB,
                                          double    *beta,
                                          double    *C, int *ic, int *jc, int *descC );

extern void psgesvd_( const char *jobu, const char *jobvt, int *m, int *n,
                                  float     *A, int *ia, int *ja, int *descA,
                                  float     *s,
                                  float     *U, int *iu, int *ju, int *descU,
                                  float     *VT, int *ivt, int *jvt, int *descVT,
                                  float     *work, int *lwork, int *info);
extern void pdgesvd_( const char *jobu, const char *jobvt, int *m, int *n,
                                  double    *A, int *ia, int *ja, int *descA,
                                  double    *s,
                                  double    *U, int *iu, int *ju, int *descU,
                                  double    *VT, int *ivt, int *jvt, int *descVT,
                                  double    *work, int *lwork, int *info);

extern void pslaset_( const char *uplo, int *m, int *n, float     *alpha, float     *beta, float     *A, int *ia, int *ja, int *descA );
extern void pdlaset_( const char *uplo, int *m, int *n, double    *alpha, double    *beta, double    *A, int *ia, int *ja, int *descA );

extern void pselset_( float     *A, int *ia, int *ja, int *descA, float     *alpha);
extern void pdelset_( double    *A, int *ia, int *ja, int *descA, double    *alpha);

extern void pslawrite_( const char** filenam, int *m, int *n, float  *A, int *ia, int *ja, int *descA, int *irwrit, int *icwrit, float  *work);
extern void pdlawrite_( const char** filenam, int *m, int *n, double *A, int *ia, int *ja, int *descA, int *irwrit, int *icwrit, double *work);

extern float pslamch_( int *ictxt, const char *cmach);
extern double pdlamch_( int *ictxt, const char*cmach);

#ifdef __cplusplus
extern "C" int indxg2p_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern "C" int indxg2l_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern "C" int indxl2g_( int *indxloc, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern "C" void Cdgerv2d(int, int, int, double*, int, int, int);
extern "C" void Cdgesd2d(int, int, int, double*, int, int, int);
#else
extern int indxg2p_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern int indxg2l_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern int indxl2g_( int *indxloc, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern void Cdgerv2d(int, int, int, double*, int, int, int);
extern void Cdgesd2d(int, int, int, double*, int, int, int);
#endif
extern int numroc_( int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern void descinit_( int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc,
				int *ictxt, int *lld, int *info);

extern void   pdgetrf_ ( int* m, int *n, double *a, int *i1, int *i2, int *desca, int* ipiv, int *info );
extern void   pdgetrs_ ( const char* trans, int* n, int* nrhs, double* A, int* ia, int* ja, int* descA, int* ippiv, double* B, int* ib, int* jb, int* descB, int* info);
extern double pdlansy_ ( const char*norm, const char*uplo, int *n, double *a, int *ia, int *ja, int *desca, double *work );
extern void   pdmatgen_( int *ictxt, const char*aform, const char*diag, int *m, int *n, int *mb, int *nb, double *a, int *lda, int *iarow, int *iacol, int *iseed, int *iroff, int *irnum, int *icoff, int *icnum, int *myrow, int *mycol, int *nprow, int *npcol );

extern void  psorglq_(int *M, int *N, int *K, float *Q, int *ia, int *ja, int *descA, float *taup, float *work,  int *lwork, int *info);
extern void  pdorglq_(int *M, int *N, int *K, double *Q, int *ia, int *ja, int *descA, double *taup, double *work,  int *lwork, int *info);

extern void pdpotrf_( const char*uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info );
extern void pdpotri_( const char*uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info );
extern void pdpotrs_( const char*uplo, int *n, int *nrhs, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb, int *info );

extern void pdposv_(const char*uplo, int *M, int *RHS, double *A, int *ia, int *ja, int *descAA, double *B, int *ib, int *jb, int *descB, int *info);
extern void psposv_(const char*uplo, int *M, int *RHS, float *A, int *ia, int *ja, int *descAA, float *B, int *ib, int *jb, int *descB, int *info);

extern void pdsymm_(const char* side, const char* uplo, int *m, int *n, double *alpha, double *A, int *ia, int *ja, int *desca, double *B, int *ib, int *jb, int *descb, double *beta, double *C, int *ic, int *jc, int *descc);
extern void pdsyr2k_(const char* uplo, const char*trans, int *n, int *k, double *alpha, double *A, int *ia, int *ja, int *desca, double *B, int *ib, int *jb, int *descb, double *beta, double *C, int *ic, int *jc, int *descc);
extern void pdgeadd_(const char* trans, int *m, int *n, double *alpha, double *A, int *ia, int *ja, int *desca, double *beta, double *C, int *ic, int *jc, int *descc);

extern void psgels_( const char* trans, int *m, int *n, int *rhs, float *A, int *ia, int *ja, int *desca, float* B, int *ib, int *jb, int *descb, float *work, long int *lwork, int *info);

extern void pssyrk_(const char* uplo, const char* trans, int *N, int *M, float *alpha, float *A, int *ia, int *ja, int *descA, float *beta, float *B, int *ib, int *jb, int *descB);

extern void psgeadd_(const char* trans, int *m, int *n, float *alpha, float *A, int *ia, int *ja, int *desca, float *beta, float *C, int *ic, int *jc, int *descc);
extern void sscal(int *myUrows, float *ss, float *U, int *ia );

extern void infog2l_( int *i, int *j, int *desca, int *nprow, int *npcol, int *myrow, int *mycol, int *ii, int *jj, int *prow, int *pcol );

#ifdef __cplusplus
}
#endif
#endif
