/*
 *
 *
 */

#ifndef _GEOMINFO_H_
#define _GEOMINFO_H_

#include <memory.h>
#include <assert.h>

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#define CHECKERR(ierr) assert(!(ierr))
#define CHECKTURE(iture) assert((iture))
#define FALSE 0
#define INITMEM(m, base) {\
  long I_; for(I_=0;I_<m;I_++) (base)[I_]=0;\
}

#define CEMS_MALLOC1D(base, nmem, type,m) {\
  (m) = (nmem) * sizeof(type);\
  (base) = (type*)malloc( (m) ); \
  assert((base) != NULL); \
  INITMEM((m),(char*)(base));\
  (m) = (nmem)* sizeof(type);;\
}

#define CEMS_MALLOC2D(base, row, col, type,i, m) {\
  CHECKERR(((row>0)&&(col>0)) == FALSE); \
  (base) = (type **) malloc( (row)* sizeof(type *)); \
  i = (row)*(col);\
  CEMS_MALLOC1D((base)[0], i, type,m);\
  for(i=1; i<row; i++) (base)[i] = (base)[i-1] + col;\
}

#define CEMS_FREE1D(base) { \
  if((base) != NULL) free((base)); (base) = NULL; \
}

#define CEMS_FREE2D(base) {\
  if((base) != NULL) { \
    CEMS_FREE1D((base)[0]); CEMS_FREE1D((base));\
  }\
}

#define SWAP(v1,v2,tmp) {\
  tmp = (v1); (v1)=(v2); (v2)=tmp;\
}


typedef struct {
  int lev;
  int used;
  int *pos;
  float **arr;
} ArrFlt;
/*
typedef struct {
  int used;
  int *pos;
  int **arr;
} ArrInt;

typedef struct{
  int lev;
  ArrFlt *nodcrd;
  ArrFlt *edglen;
  ArrInt *edgnodpth;
  ArrInt *pthedg;
} PthInfo;
*/

typedef struct{
  int nnod,npth,nedg,nbrd,nsrc;
  int *edgsrc;
  float *edglen;
  int *pthedg;
  int *edgnodpth;
} InfoEdgPth;
	
/* Function */
//PthInfo *PthInfoCreate();
//void PthInfoDestroy(PthInfo **pth);
ArrFlt *ArrFltCreate();
void ArrFltDestroy(ArrFlt **arr);
InfoEdgPth *InfoEdgPthCreate();
void InfoEdgPthDestroy(InfoEdgPth **pth);

void PthEdgLen(int iTotEdg, int iTotPth, int (*ppiEdgPth)[4],
    float (*ppfNodCrd)[3], int (*)[3],int (*ppiPth3Edg)[3], float *pfLen);

double AvgEdgLen(int iTot, float *pfLen);

float *LocateNod(int i, ArrFlt *nodcrd);

int FineMesh(int iLev, InfoEdgPth **oldpth, InfoEdgPth **newpth, ArrFlt *nodcrd, double *dnewAvgLen,
    float fNedLen, double *dMem);

void GetPthNod(int oldpthtag, int *pthedg, InfoEdgPth *oldpth, int *nod);

void ShapeOfPatch(char *pchOut, InfoEdgPth *pth);

int TessCvt(InfoEdgPth **oldpth, InfoEdgPth **newpth, ArrFlt *nodcrd, 
    float **ppfNodCrd, double *dOrgAvgLen, float fNedLen, float fLambda, int *iTesLev, int *iLev, double *dMem);

int DelOpenBd(InfoEdgPth *edgpth);


#endif /* _GEOMINFO_H_ */
