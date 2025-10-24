#ifndef _F2C_H_
#define _F2C_H_

#include <complex.h>
#include "settings.h"
#include "cmumps_c.h"

#define INTEL 1

#ifdef INTEL
        #define FORTRAN(x) x##_
#else
        #define FORTRAN(x) x##_
#endif

void cgemv(const char *trans, const int *m, const int *n, const field *alpha,
           const field *a, const int *lda, const field *x, const int *incx,
           const field *beta, field *y, const int *incy);

void caxpy(const int *n, const field *alpha, const field *x, const int *incx, field *y, const int *incy);

void cgemm(const char *transa, const char *transb, const int *m, const int *n, const int *k,
           const field *alpha, const field *a, const int *lda,
           const field *b, const int *ldb, const field *beta,
           field *c, const int *ldc);

void    cscal(const int *n, const field *a, field *x, const int *incx);

void ctrmm(const char *side, const char *uplo, const char *transa, const char *diag,
           const int *m, const int *n, const field *alpha,
           const field *a, const int *lda, field *b, const int *ldb);

void casum_(const int *N, const field *DX, const int *INCX, real *csum);

void ctrtrs( const char* uplo, const char* trans, const char* diag,
             const int* n, const int* nrhs, const field* a,
             const int* lda, field* b, const int* ldb,
             int* info );

void cgeru(const int *m, const int *n, const field *alpha,
           const field *x, const int *incx, const field *y, const int *incy,
           field *a, const int *lda);

void ctrsm(const char *side, const char *uplo, const char *transa, const char *diag,
           const int *m, const int *n, const field *alpha,
           const field *a, const int *lda, field *b, const int *ldb);

void cgesvd( const char* jobu, const char* jobvt, const int* m,
             const int* n, field* a, const int* lda, real* s,
             field* u, const int* ldu, field* vt,
             const int* ldvt, field* work, const int* lwork,
             real* rwork, int* info );

void cgeqrf( const int* m, const int* n, field* a,
             const int* lda, field* tau, field* work,
             const int* lwork, int* info );

void cunmqr( const char* side, const char* trans, const int* m,
             const int* n, const int* k, const field* a,
             const int* lda, const field* tau, field* c,
             const int* ldc, field* work, const int* lwork,
             int* info );

void 
mycdot_(int *dim, const field *v, const int *i1, const field *v2, const int *i2, field *alpha);

void FORTRAN(opfile)(char *runname, int *clength, int *mxmat);

 void FORTRAN(read_fem_inp)(int *mvolnode, int *mvolele, int *mvoledge, int *msurnode, int *msurele, int *msuredge, int *msurstore,
                        int *maxnode, int *maxpatch, int *maxedge, int *kein, int *ipol, int *tyrcs, real *freq, real *ittol,
                        int (*thetai)[3], int (*phii)[3], int (*thetas)[3], int (*phis)[3], int *mxmat, field *ep, field *mu, 
                        field *xi, int *mdef1, int *mdef2);


 void FORTRAN(read_geom_fem)(int *mvolnode, int *mvolele, int *mvoledge,
                          int *msurnode, int *msurele, int *msuredge, int *kein, int *kc,
                         int *mvolstore, int *msurstore,
                         real *xyz, int *lv, int *le, int *lp, int *lc, int *ipin, int *kefemabc, int *kefemsrc);

 void FORTRAN(fem)(int *mvolnode, int *mvolele, int *mvoledge, int *msuredge, int *msurele, int *bsurele, 
                   int *hfemstore, int *hfemnum, int *kein, int *ipin, real *freq,
                   int *mxmat, field *ep, field *mu, field *xi, real *xyz, int *lv, int *le, 
                   int *lp, int *lc, int *bfface, int *blpt, int *edgesign2, int *matrixi, int *matrixj, field *matrix, int *kefemabc, int* kefemsrc,
                   int *nsurele, int *nfface, int *nlpt, real *et, field *kz);

void FORTRAN(sparse2csr)(int *mvoledge, int *hfemnum, int *matrixi, int *matrixj, field *matrix, int *iqvol, int *ipvol, field *volmat);

void FORTRAN(getfaceinfo)(int *mvolele, int *mvoledge, int *mvolface, int *lv, int *lf, int *le, int *fface, int *fedge);

void FORTRAN(matset)(int *mvoledge, int *hfemnum, int *iqvol, int *ipvol, field *volmat, int *irn, int *jcn, CMUMPS_COMPLEX *a);

void FORTRAN(find_face)(int *mvolnode, int *mvolele, real *xyz, int *lv, int *sfface, int *ssurele,
                        int *facetotal,int *mdef1, int *mdef2);

void FORTRAN(add_source)(real *wthetai, real *wphii, int *ipol, real *freq, int *mvolnode, int *mvolele,
                         int *mvoledge, real *xyz, int *lv, int *le, int *bfface, int *blpt, int *, field *crhs);

void FORTRAN(farfld)(real *wthetas,real *wphis, field *crhs, real *freq, int *mvolnode, int *mvolele, int *mvoledge,
                     real *xyz, int *lv, int *le, int *sfface, int *ssurele, real *totalp);


void FORTRAN(initialvector)(field *v, int *n);
void FORTRAN(vnorm)(field *v, int *iLen, double *norm );
void FORTRAN(csub)(int *iLen, field *v1, field *v2, field *v3);
void FORTRAN(zdotc)(int *iLen, field *v1, field *v2, double complex *value);

void FORTRAN(bicgstab1)(int *iLen, int *itr, int *iflg,
                                double complex *rho1, double complex *rho, double complex *alpha, double complex *wj, double complex *beta,
                                field *p, field *v, field *r );

void FORTRAN(bicgstab2)(int *iLen,
                        double complex *alpha, double complex *rho, double complex *cTmp,
                        field *x, field *r, field *v, field *PHAT);

void FORTRAN(bicgstab3)(int *iLen,
                        double complex *cTmp1, double complex *cTmp2, double complex *wj,
                        field *x, field *SHAT, field *r, field *T);

void FORTRAN(matrix_vector)(int *mvoledge, int *iqvol, int *ipvol, field *volmat, field *cwork, field *cwork2);

#endif // _PH_F2C_H_

