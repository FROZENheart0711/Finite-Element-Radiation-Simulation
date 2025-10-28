#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "basic.h"
#include "f2c.h"
#include "time.h"
#include <mpi.h>
#include "cmumps_c.h"

#include <stdlib.h>    // for system() and file IO
#include <stdbool.h>   // for bool type
#include <math.h>      // for fabs(), sqrt()
#include "adaptive.h"  // 自定义头文件，声明 compute_error/refine_mesh

#include <complex.h>
#include "compute_element_errors.h"


#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

// ==========================================================
// 封装原 main() 主体为一个函数（原求解流程）
// ==========================================================
void run_fem_solver(const char *mesh_file)
{
    int i,j,n;
    clock_t start, finish;
    double duration;
    CMUMPS_STRUC_C mumps_par;

    // ---- for FEM ---- //
    int mxmat;
    int bsurele, hfemstore, hfemnum;
    int nsurele;
    field *ep, *mu, *xi;
    int ipol, tyrcs;
    real freq, ittol;
    int thetai[3],phii[3],thetas[3],phis[3];
    real rad,pi;
    real wthetai,wphii,wthetas,wphis;
    real amp,phase;
    int dthetai,dphii,dthetas,dphis;
    int thetai_ang,phii_ang;
    int thetait,phiit,thetast,phist;
    int npol,ipl;
    int *sfface, ssurele, facetotal;
    double totaltime;
    int mdef1, mdef2;
    field *kz;

    // ---- geometry and FEM variables ---- //
    int max_pardepth;
    int mvolnode,mvoledge,mvolele;
    int msurnode, msurele, msuredge;
    int mvolstore,msurstore;
    int maxnode, maxpatch, maxedge;
    int *lv,*le;
    real *xyz;
    int jmem;
    int clength;
    char pchInNam[256],pchInFil[256],runname[256];
    FILE *fp;
    int *lp,*lc,*ipin,*ipout,*ipvol,*iqvol;
    field *volmat, *crhs;
    int kein,kc;
    int kefemabc,kefemsrc;
    int **edgesign;
    double ttime=0,ttime2=0;
    int mvolface;
    int *lf, *fface,*fedge,(*edgeco)[2];
    int *edgsign;
    int *edgsigntmp;
    double tmem=0;
    int *admismat,inum1=0,inum2=0;

    pi=3.1415926;
    rad=pi/180.0;
    omp_set_num_threads(10);

    strcpy(pchInFil, mesh_file);
    strcpy(runname, pchInFil);
    strcat(runname, " ");
    clength = strlen(runname);

    // =======================================================
    // 1. 读取输入、初始化 FEM 参数
    // =======================================================
    FORTRAN(opfile)(runname,&clength, &mxmat);

    ep = (field(*)) allocfield(3*3*mxmat);
    mu = (field(*)) allocfield(3*3*mxmat);
    xi = (field(*)) allocfield(mxmat);
    kz = (field(*)) allocfield(1);

    FORTRAN(read_fem_inp)(&mvolnode,&mvolele,&mvoledge,&msurnode, &msurele, &msuredge, &msurstore,
                          &maxnode, &maxpatch,&maxedge,&kein, &ipol, &tyrcs, &freq, &ittol,
                          &thetai, &phii, &thetas, &phis, &mxmat, ep, mu, xi, &mdef1, &mdef2);

    xyz = (real(*)) allocreal(3*mvolnode);
    lv = (int(*)) allocint(5*mvolele);
    le = (int(*)) allocint(6*mvolele);
    lp = (int(*)) allocint(3*msurele);
    lc = (int(*)) allocint(4*msurele);
    ipin = (int*) allocint(kein);
    edgesign = (int**)malloc(sizeof(int*) *mvoledge);
    for(i=0;i<mvoledge;i++) edgesign[i]=NULL;

    printf("Began to read_geom_fem\n");
    start = clock();

    FORTRAN(read_geom_fem)(&mvolnode, &mvolele, &mvoledge,
                           &msurnode, &msurele, &msuredge, &kein, &kc,
                           &mvolstore, &msurstore,
                           xyz, lv, le, lp, lc, ipin,
                           &kefemabc, &kefemsrc);

    sfface = (int *) allocint(16*mvolele);

    FORTRAN(find_face)(&mvolnode,&mvolele,xyz,lv,sfface,&ssurele,
                       &facetotal, &mdef1, &mdef2);

    // =======================================================
    // 2. 组装 FEM 系统
    // =======================================================
    printf("Began to assemble_fem\n");

    int *matrixi, *matrixj;
    field *matrix;
    int *edgesign2;
    int *bfface, *blpt;
    int *nfface, *nlpt;
    real *et;

    matrixi = (int(*)) allocint(36*mvolele+msurele*msurele*9);
    matrixj = (int(*)) allocint(36*mvolele+msurele*msurele*9);
    matrix  = (field(*)) allocfield(36*mvolele+msurele*msurele*9);
    edgesign2 = (int(*)) allocint(mvoledge);

    bfface = (int(*)) allocint(4*4*mvolele);
    blpt   = (int(*)) allocint(3*4*mvolele);
    nfface = (int(*)) allocint(4*4*mvolele);
    nlpt   = (int(*)) allocint(3*4*mvolele);
    et     = (real(*)) allocreal(mvoledge*2+msurele*2);

    FORTRAN(fem)(&mvolnode,&mvolele,&mvoledge,&msuredge,&msurele, &bsurele,&hfemstore,&hfemnum,&kein, ipin, &freq,
                 &mxmat, ep, mu, xi, xyz, lv, le, lp, lc, bfface, blpt, edgesign2,matrixi, matrixj, matrix, 
                 &kefemabc, &kefemsrc, &nsurele, nfface, nlpt, et, kz);

    iqvol= (int(*)) allocint(hfemnum);
    ipvol= (int(*)) allocint(hfemnum);
    volmat= (field(*)) allocfield(hfemnum);

    FORTRAN(sparse2csr)(&mvoledge, &hfemnum, matrixi, matrixj, matrix, iqvol, ipvol, volmat);

    // =======================================================
    // 3. 调用 MUMPS 求解器
    // =======================================================
    printf(" ================== Begin calling MUMPS ====================\n");

    mumps_par.job=JOB_INIT;
    mumps_par.par=1;
    mumps_par.sym=0;
    mumps_par.comm_fortran=USE_COMM_WORLD;
    cmumps_c(&mumps_par);

    mumps_par.n=mvoledge;
    mumps_par.nz=hfemnum;

    mumps_par.irn = (int *) allocint(hfemnum);
    mumps_par.jcn = (int *) allocint(hfemnum);
    mumps_par.a   = (CMUMPS_COMPLEX *) allocfield(hfemnum);
    mumps_par.rhs = (CMUMPS_COMPLEX *) allocfield(mvoledge);

    FORTRAN(matset)(&mvoledge, &hfemnum, iqvol, ipvol, volmat, mumps_par.irn, mumps_par.jcn, mumps_par.a);
    start = clock();
    mumps_par.job=4;
    cmumps_c(&mumps_par);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf("Time for MUMPS factorization: %f s\n", duration);

    field *cwork = (field *) allocfield(20*mvoledge);
    crhs = (field *) allocfield(mvoledge);
    for(i=0;i<20*mvoledge;i++) cwork[i]=0;

    // =======================================================
    // 4. 加载右端项 + 求解
    // =======================================================
    FORTRAN(rhside)(&mvolnode,&mvoledge,&mvolele,xyz,&nsurele,
                    nfface,nlpt,le,et,crhs,kz, ep, mu, &mxmat, lv, &amp, &phase);

    for(i=0;i<mvoledge;i++)
        if(edgesign2[i]==1)
            crhs[i]=0;

    totaltime=0;
    for(j=0;j<20;j++)
    {
        for(i=0;i<mvoledge;i++)
        {
            mumps_par.rhs[i].r=creal(crhs[i]);
            mumps_par.rhs[i].i=cimag(crhs[i]);
        }
        start = clock();
        mumps_par.job=3;
        cmumps_c(&mumps_par);
        finish = clock();
        duration = (double)(finish - start) / CLOCKS_PER_SEC;
        totaltime+=duration;
    }
    printf("Average MUMPS solve time: %f s\n", totaltime/20.0);

    for(i=0;i<mvoledge;i++)
        crhs[i] = (mumps_par.rhs[i].r + I*mumps_par.rhs[i].i);

    // =======================================================
    // 5. 输出结果 (Far field, S11, etc.)
    // =======================================================
    FORTRAN(s11_zinvol)(&mvolnode,&mvolele,&mvoledge,
                        xyz,lv,le,nfface,nlpt,&nsurele,
                        crhs,et,&mxmat,ep,mu,&amp,&phase);

    // =======================================================
    // 6. 导出解向量供 compute_element_errors 使用
    // =======================================================
    printf("Exporting solution vector to solution.dat for error estimation...\n");

    FILE *fsol = fopen("solution.dat", "w");
    if (!fsol)
    {
        fprintf(stderr, "Error: cannot open solution.dat for writing!\n");
    }
    else
    {
        for (i = 0; i < mvoledge; i++)
        {
            fprintf(fsol, "%d %.16e %.16e\n", i+1, creal(crhs[i]), cimag(crhs[i]));
        }
        fclose(fsol);
        printf("Solution exported: %d entries written.\n", mvoledge);
    }

    // 注意：crhs 是求解完成后的复数电场向量，对应每个边的自由度
    // Fortran 侧会从 solution.dat 中读取这些值并计算单元误差
}


// ==========================================================
// 主程序入口（外层自适应控制）
// ==========================================================
/* 外部函数声明 */

extern void compute_element_errors_(void *ptr_crhs, void *ptr_xyz, void *ptr_lv, void *ptr_le,
                                     int mvolnode, int mvolele, int mvoledge, void *ptr_errors);

int main(int argc, char **argv)
{
    int adapt_iter = 0;
    const int max_adapt_iter = 10;
    const double tol = 1e-3;
    double global_error = 1.0;
    char mesh_file[256];

    if (argc < 2) {
        printf("Usage: %s mesh_file\n", argv[0]);
        return -1;
    }

    strcpy(mesh_file, argv[1]);
    printf("=== Starting Adaptive FEM Solver ===\n");

    while (adapt_iter < max_adapt_iter && global_error > tol)
    {
        printf("\n=============================\n");
        printf(" Adaptive iteration %d\n", adapt_iter);
        printf("=============================\n");

        int mvolnode, mvolele, mvoledge;
        double *xyz = NULL;
        int *lv = NULL, *le = NULL;
        double complex *crhs = NULL;

        // (1) FEM 求解
        run_fem_solver(mesh_file, &mvolnode, &mvolele, &mvoledge, &xyz, &lv, &le, &crhs);

        // (2) Fortran 误差估计（直接传指针）
        double *errors = (double*) malloc(sizeof(double) * mvolele);
        FORTRAN(compute_element_errors)((void*)crhs, (void*)xyz, (void*)lv, (void*)le,
                                 mvolnode, mvolele, mvoledge, (void*)errors);

        // (3) 全局误差
        global_error = 0.0;
        for (int i = 0; i < mvolele; ++i)
            global_error += errors[i] * errors[i];
        global_error = sqrt(global_error);
        printf("Current global error estimate = %.6e\n", global_error);

        // (4) 收敛判断
        if (global_error < tol) {
            printf("Converged! Adaptive refinement stopped.\n");
            free(errors);
            free(xyz); free(lv); free(le); free(crhs);
            break;
        }

        // (5) 网格细化（临时仍用文件方式）
        FILE *f = fopen("errors.dat", "w");
        for (int i = 0; i < mvolele; i++)
            fprintf(f, "%d %.8e\n", i+1, errors[i]);
        fclose(f);

        printf("Refining mesh...\n");
        refine_mesh(mesh_file, "mesh_refined.msh", "errors.dat", 0.3);
        strcpy(mesh_file, "mesh_refined.msh");

        free(errors);
        free(xyz); free(lv); free(le); free(crhs);
        adapt_iter++;
    }

    printf("Adaptive loop finished after %d iterations.\n", adapt_iter);
    return 0;
}

