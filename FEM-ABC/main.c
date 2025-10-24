#include <stdio.h>
#include <string.h>
#include <omp.h>
#include "basic.h"
#include "f2c.h"
#include "time.h"

#include <mpi.h>
#include "cmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

//---- below added for the FEM -----//
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
// for nested dissection 
  int mvolface;
  int *lf, *fface,*fedge,(*edgeco)[2];
// help to define which one has been treated
  int *edgsign;
// temp memory
  int *edgsigntmp;
  double tmem=0;
  int *admismat,inum1=0,inum2=0;
int main(int argc,char **argv)
{
  int i,j,n;
  clock_t start, finish;
  double  duration;
    CMUMPS_STRUC_C mumps_par;
// for fem
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
 
  pi=3.1415926;
  rad=pi/180.0;

  omp_set_num_threads (10);

// local variables
  int *ip;
// read geometry ...
//-------------------------------
//ccccccccccc le and lv should minus one to make it started from zero.
  strcpy((pchInFil),*(argv+1));
  strcpy((runname),pchInFil);
  strcat(runname," ");
  clength=0;
  clength=strlen(runname);

// read FEM geom info ======================================================

  FORTRAN(opfile)(runname,&clength, &mxmat);

  ep = (field(*)) allocfield(3*3*mxmat);
  mu = (field(*)) allocfield(3*3*mxmat);
  xi = (field(*)) allocfield(mxmat);
  kz = (field(*)) allocfield(1);

  FORTRAN(read_fem_inp)(&mvolnode,&mvolele,&mvoledge,&msurnode, &msurele, &msuredge, &msurstore,
                        &maxnode, &maxpatch,&maxedge,&kein, &ipol, &tyrcs, &freq, &ittol,
                        &thetai, &phii, &thetas, &phis, &mxmat, ep, mu, xi, &mdef1, &mdef2);

  printf("***********************************************************\n");

//-------------------------------
  xyz = (real(*)) allocreal(3* mvolnode);
  lv = (int(*)) allocint(5*mvolele);
  le = (int(*)) allocint(6*mvolele);
  lp = (int(*)) allocint(3*msurele);
  lc = (int(*)) allocint(4*msurele);
  ipin = (int*) allocint(kein);
  edgesign = (int**)malloc(sizeof(int*) *mvoledge);

  for(i=0;i<mvoledge;i++)
  edgesign[i]=NULL;

//====read FEM geom  
      printf("Began to read_geom_fem\n");
   start = clock();

  FORTRAN(read_geom_fem)(&mvolnode, &mvolele, &mvoledge,
                         &msurnode, &msurele, &msuredge, &kein, &kc,
                         &mvolstore, &msurstore,
                         xyz, lv, le, lp, lc, ipin,
                         &kefemabc, &kefemsrc);

   sfface= (int *) allocint(16*mvolele);

   FORTRAN(find_face)(&mvolnode,&mvolele,xyz,lv,sfface,&ssurele,
                      &facetotal, &mdef1, &mdef2);

   printf("***********************************************************\n");
   printf("Began to assemble_fem\n");

  int *matrixi, *matrixj;
  field *matrix;
  int *edgesign2;
  int *bfface, *blpt;
  int *nfface, *nlpt;
  real *et;
 
  matrixi = (int(*)) allocint(36*mvolele+msurele*msurele*9);
  matrixj = (int(*)) allocint(36*mvolele+msurele*msurele*9);
  matrix = (field(*)) allocfield(36*mvolele+msurele*msurele*9);
  edgesign2 = (int(*)) allocint(mvoledge);

  bfface = (int(*)) allocint(4*4*mvolele);
  blpt = (int(*)) allocint(3*4*mvolele);
  nfface = (int(*)) allocint(4*4*mvolele);
  nlpt = (int(*)) allocint(3*4*mvolele);
  et = (real(*)) allocreal(mvoledge*2+msurele*2);


  FORTRAN(fem)(&mvolnode,&mvolele,&mvoledge,&msuredge,&msurele, &bsurele,&hfemstore,&hfemnum,&kein, ipin, &freq,
               &mxmat, ep, mu, xi, xyz, lv, le, lp, lc, bfface, blpt, edgesign2,matrixi, matrixj, matrix, 
               &kefemabc, &kefemsrc, &nsurele, nfface, nlpt, et, kz);
   printf(" FEM calling over!\n");
   iqvol= (int(*)) allocint(hfemnum);
   ipvol= (int(*)) allocint(hfemnum);
   volmat= (field(*)) allocfield(hfemnum);

   FORTRAN(sparse2csr)(&mvoledge, &hfemnum, matrixi, matrixj, matrix, iqvol, ipvol, volmat);

   finish = clock(); 
   duration = (double)(finish - start) / CLOCKS_PER_SEC; 
   printf( "Time for FEM is %f seconds\n\n", duration );
//=======================================================================
//   get surface information 
  mvolface=0;
  fface = (int(*)) allocint(5*4*mvolele); // this can be rewitten later
  fedge = (int(*)) allocint(3*4*mvolele);
  lf = (int(*)) allocint(4*mvolele);

  FORTRAN(getfaceinfo)(&mvolele,&mvoledge,&mvolface,lv,lf,le,fface,fedge);  
  
  edgsign= (int(*)) allocint(mvoledge);
  edgsigntmp= (int(*)) allocint(mvoledge);

  for(i=0;i<mvoledge;i++)
    edgsign[i]=0;

  for(i=0;i<mvoledge;i++)
    edgsign[i]=edgesign2[i];

//=======================================================================
  edgeco = (int(*)[2]) allocint(sizeof(int[2]) * mvoledge);
   for(i=0;i<mvoledge;i++)
    for(j=0;j<2;j++)
     edgeco[i][j]=0;

   for(i=0;i<mvolele;i++)
    {
     edgeco[le[i*6+0]-1][0]=lv[i*5+0]-1;
     edgeco[le[i*6+0]-1][1]=lv[i*5+1]-1;

     edgeco[le[i*6+1]-1][0]=lv[i*5+0]-1;
     edgeco[le[i*6+1]-1][1]=lv[i*5+2]-1;

     edgeco[le[i*6+2]-1][0]=lv[i*5+0]-1;
     edgeco[le[i*6+2]-1][1]=lv[i*5+3]-1;

     edgeco[le[i*6+3]-1][0]=lv[i*5+1]-1;
     edgeco[le[i*6+3]-1][1]=lv[i*5+2]-1;

     edgeco[le[i*6+4]-1][0]=lv[i*5+1]-1;
     edgeco[le[i*6+4]-1][1]=lv[i*5+3]-1;

     edgeco[le[i*6+5]-1][0]=lv[i*5+2]-1;
     edgeco[le[i*6+5]-1][1]=lv[i*5+3]-1;
    }
   
    printf(" ================== Begin calling MUMPS ====================\n");
    mumps_par.job=JOB_INIT;
    mumps_par.par=1;
    mumps_par.sym=0;
    mumps_par.comm_fortran=USE_COMM_WORLD;
    cmumps_c(&mumps_par);
    mumps_par.n=mvoledge;
    mumps_par.nz=hfemnum;

    mumps_par.irn= (int *) allocint(hfemnum);
    mumps_par.jcn= (int *) allocint(hfemnum);
    mumps_par.a= (CMUMPS_COMPLEX *) allocfield(hfemnum);
    mumps_par.rhs= (CMUMPS_COMPLEX *) allocfield(mvoledge);

      for(j=0;j<mumps_par.nz;j++)
      {
        mumps_par.irn[j]=0;mumps_par.jcn[j]=0;mumps_par.a[j].r=0;mumps_par.a[j].i=0;
                        }
      for(j=0;j<mumps_par.n;j++)
      {
        mumps_par.rhs[j].r=0;
        mumps_par.rhs[j].i=0;
                        }
    FORTRAN(matset)(&mvoledge, &hfemnum, iqvol, ipvol, volmat, mumps_par.irn, mumps_par.jcn, mumps_par.a);
   start = clock();
    mumps_par.job=4;
    cmumps_c(&mumps_par);
   finish = clock();
   duration = (double)(finish - start) / CLOCKS_PER_SEC;
   printf( "Time for mumps factorization is %f seconds\n\n", duration );

    mumps_par.icntl[3]=1;
    mumps_par.icntl[2]=0;
   

   field *cwork;
   cwork= (field *) allocfield(20*mvoledge);
   crhs= (field *) allocfield(mvoledge);

   for(i=0;i<20*mvoledge;i++)
      cwork[i]=0;

// addsource          
                    
      FORTRAN(rhside)(&mvolnode,&mvoledge,&mvolele,xyz,&nsurele,
             nfface,nlpt,le,et,crhs,kz, ep, mu, &mxmat, lv, &amp, &phase);

        for(i=0;i<mvoledge;i++)
          if(edgesign2[i]==1)
            crhs[i]=0;

//      for(i=0;i<mvoledge;i++){
//              printf("crhs: %f %f %d\n",creal(crhs[i]),cimag(crhs[i]),i);
//      }

      totaltime=0;
      for(j=0;j<40;j++)
     {
      for(i=0;i<mvoledge;i++)
     {mumps_par.rhs[i].r=creal(crhs[(i)]);
      mumps_par.rhs[i].i=cimag(crhs[(i)]);
     }
    start = clock();
    mumps_par.job=3;
    cmumps_c(&mumps_par);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    totaltime=totaltime+duration;
     }
     printf( "Average time for mumps Solve is %f seconds\n\n", totaltime/20.0 );

      for(i=0;i<mvoledge;i++)
     {crhs[i]=(mumps_par.rhs[i].r+I*mumps_par.rhs[i].i);}


        for(i=0;i<mvoledge;i++)
          if(edgesign2[i]==1)
            crhs[i]=0;

      if(tyrcs==1)
        {
         thetas[0]=thetai_ang;
         thetas[1]=thetai_ang;
         thetas[2]=1;
         phis[0]=phii_ang;
         phis[1]=phii_ang;
         phis[2]=1;
        }
        
    int nleg;
    real *fxgl,*fwgl;
    real *ftheta,*fphi;
    real cweigh;
    real totalp;

    totalp=0;

    nleg=90;   
    fxgl = (real(*)) allocreal(nleg);
    fwgl = (real(*)) allocreal(nleg);  
    ftheta = (real(*)) allocreal(2*nleg);
    fphi = (real(*)) allocreal(4*nleg);

    dthetas=(180-0)/nleg;
    dphis = (360-0)/(nleg*2);

    FORTRAN(numint)(&nleg,&nleg,fxgl,fwgl,ftheta,fphi);

        for(thetast=0; thetast<nleg;thetast++) {
          
          cweigh = rad*dphis*fwgl[thetast];

           for(phist=0; phist<(nleg*2);phist++)
            {  wthetas=rad*dthetas*thetast;
               wphis=rad*dphis*phist;
  
         FORTRAN(farfld1)(&wthetas,&wphis,crhs,&freq,&mvolnode,&mvolele,&mvoledge,
                         xyz,lv,le, sfface, &ssurele,&totalp,&cweigh);

         
            } }

      
      dthetas = 0;
      if (thetas[2]>1) 
          dthetas=(thetas[1]-thetas[0])/(thetas[2]-1.0);
      dphis = 0;
      if (phis[2]>1) 
         dphis = (phis[1]-phis[0])/(phis[2]-1.0);

      for(thetast=1; thetast<=thetas[2];thetast++)
         for(phist=1; phist<=phis[2];phist++)
          {  wthetas=rad*(thetas[0]+dthetas*(thetast-1));
             wphis=rad*(phis[0]+dphis*(phist-1));

       FORTRAN(farfld)(&wthetas,&wphis,crhs,&freq,&mvolnode,&mvolele,&mvoledge,
                       xyz,lv,le, sfface, &ssurele, &totalp);
          }

      FORTRAN(s11_zinvol)(&mvolnode,&mvolele,&mvoledge,
          xyz,lv,le,nfface,nlpt,&nsurele,
          crhs,et,&mxmat,ep,mu,&amp,&phase);

//  for(i=1;i<10;i++)
//  printf("solution: %f %f %d\n",creal(crhs[mvoledge-i]),cimag(crhs[mvoledge-i]),mvoledge-i);


  exit(0);
}

