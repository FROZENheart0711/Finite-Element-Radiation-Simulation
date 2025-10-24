#include "geominfo.h"
#include "num2char.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

#define FILE_MODE (S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH)
#define DIR_MODE (FILE_MODE |S_IXUSR |S_IXGRP |S_IXOTH)
#define _GNU_SOURCE
extern void fstrcat_(char *pchIn,char *tailname,char *pcFilNam);
void writenodcrd_(int *writeFlg, char *pchIn, int *nnode,float *nodcrd){

  char pcFilNam[256],tailname[256];
  FILE *fp;
  int i, k, j;
  strcpy(tailname,".bxyz");
  fstrcat_(pchIn,tailname,pcFilNam);
    if((fp=fopen(pcFilNam,"wb"))==NULL){
      printf("Cannot open the file %s !\n", pcFilNam);
      exit(0);
    }
    if(fwrite(nnode,sizeof(int),1,fp) != 1){
      printf("file(%s.bxyz) write error 1!\n", pcFilNam);
      exit(0);
             }
  i=3*(*nnode);
  if(fwrite(nodcrd, sizeof(float), i, fp) !=i){
    printf("file(%s) write error!\n",pcFilNam);
    exit(0);
  }
  fclose(fp);
  printf(" NodCtr: %s\n",pcFilNam);

}
void writepthnod_(int *writeFlg, char *pchIn, int *iMaxPth, int *ppiPthNod){

  char pcFilNam[256],tailname[256];
  FILE *fp;
  int i,k, nod[3], j;
  strcpy(tailname,".bipat");
  fstrcat_(pchIn,tailname,pcFilNam);
//  printf(" %11d%11d%11d\n",*(ppiPthNod+3*(*iMaxPth-1)),*(ppiPthNod+3*(*iMaxPth-1)+1),*(ppiPthNod+3*(*iMaxPth-1)+2));
//  strcpy(pcFilNam, pchIn);
  if((*writeFlg)==0){
//    strcat(pcFilNam, ".bipat");
    if((fp=fopen(pcFilNam,"wb"))==NULL){
      printf("Cannot open the file %s !\n", pcFilNam);
      exit(0);
    }
    if(fwrite(iMaxPth,sizeof(int),1,fp) !=1){
      printf("file(%s) write error!\n",pcFilNam);
      exit(0);
    }
    i = 3*(*iMaxPth);
    if(fwrite( ppiPthNod, sizeof(int),i,fp) !=i){
      printf("file(%s) write error!\n",pcFilNam);
      exit(0);
    }
  }
  else{
    strcat(pcFilNam, ".ipat");
    if((fp=fopen(pcFilNam,"w"))==NULL){
      printf("Cannot open the file %s !\n", pcFilNam);
      exit(0);
    }



    /* if(fwrite(&(iMaxPth),sizeof(int),1,fp) !=1){
       printf("file(%s) write error!\n",pcFilNam);
       exit(0);
       }*/
    fprintf(fp,"%11d\n",*iMaxPth);
    i = 3*(*iMaxPth); 
    /*  if(fwrite( ppiPthNod, sizeof(int),i,fp) !=i){
	printf("file(%s) write error!\n",pcFilNam);
	exit(0);
	}*/
    for(j=0;j<(*iMaxPth);j++){
      fprintf(fp,"%11d%11d%11d\n",*(ppiPthNod+3*j),*(ppiPthNod+3*j+1),*(ppiPthNod+3*j+2));
    }
  }

  fclose(fp);

  printf(" PthNod: %s\n",pcFilNam);
}
void writeedgctr_(int *writeFlg, char *pchIn, int *iMaxEdg, float *pfEdgCtr){

  char pcFilNam[256],tailname[256];
  FILE *fp;
  int i,j;
  strcpy(tailname,".bedgctr");
//  slength=strlen(*pchIn);
//  strcpy(pcFilNam, pchIn);
//  strcat(pcFilNam, ".mommax");
  fstrcat_(pchIn,tailname,pcFilNam);  
//  strcpy(pcFilNam, pchIn);
  if((*writeFlg)==0){
//  strcat(pcFilNam, ".bedgctr");
  if((fp=fopen(pcFilNam,"wb"))==NULL){
    printf("Cannot open the file %s !\n", pcFilNam);
    exit(0);
  }
  if(fwrite(iMaxEdg,sizeof(int),1,fp) !=1){
    printf("file(%s) write error!\n",pcFilNam);
    exit(0);
  }
  i=3*(*iMaxEdg);
  if(fwrite(pfEdgCtr, sizeof(float), i, fp) !=i){
    printf("file(%s) write error!\n",pcFilNam);
    exit(0);
  }
  }else{
    strcat(pcFilNam, ".edgctr");
    if((fp=fopen(pcFilNam,"w"))==NULL){
      printf("Cannot open the file %s !\n", pcFilNam);
      exit(0);
    }
  fprintf(fp,"%11d\n",*iMaxEdg);
  
  
  i=3*(*iMaxEdg);
/*  if(fwrite(pfEdgCtr, sizeof(float), i, fp) !=i){
    printf("file(%s) write error!\n",pcFilNam);
    exit(0);
  }*/
  for(j=0;j<(*iMaxEdg);j++){
    fprintf(fp,"%14.6e%14.6e%14.6e\n",pfEdgCtr[j*3],pfEdgCtr[j*3+1],pfEdgCtr[j*3+2]);
  }
  }
  fclose(fp);
  printf(" EdgCtr: %s\n",pcFilNam);
  
}
void writeedgpth_(int *writeFlg, char *pchIn, int *iMaxEdg, int *ppiEdgPth){

  char pcFilNam[256],tailname[256];
  FILE *fp;
  int i,j;
  strcpy(tailname,".bedge");
  fstrcat_(pchIn,tailname,pcFilNam);  
//  strcpy(pcFilNam, pchIn);
  if((*writeFlg)==0){
//  strcat(pcFilNam, ".bedge");
  if((fp=fopen(pcFilNam,"wb"))==NULL){
    printf("Cannot open the file %s !\n", pcFilNam);
    exit(0);
  }
  if(fwrite(iMaxEdg,sizeof(int),1,fp) !=1){
    printf("file(%s) write error!\n",pcFilNam);
    exit(0);
  }
  i=4*(*iMaxEdg);
  if(fwrite(ppiEdgPth, sizeof(int), i, fp) !=i){
    printf("file(%s) write error!\n",pcFilNam);
    exit(0);
  }
  }else{
    strcat(pcFilNam, ".edge");
    if((fp=fopen(pcFilNam,"w"))==NULL){
      printf("Cannot open the file %s !\n", pcFilNam);
      exit(0);}		  
      /*     if(fwrite(&iMaxEdg,sizeof(int),1,fp) !=1){
	     printf("file(%s) write error!\n",pcFilNam);
	     exit(0);
	     }*/
      fprintf(fp,"%11d\n",*iMaxEdg);


      i=4*(*iMaxEdg);
      /*if(fwrite(ppiEdgPth, sizeof(int), i, fp) !=i){
	printf("file(%s) write error!\n",pcFilNam);
	exit(0);
	}*/
      for(j=0;j<(*iMaxEdg);j++){
	fprintf(fp,"%11d%11d%11d%11d\n",ppiEdgPth[4*j],ppiEdgPth[4*j+1],ppiEdgPth[4*j+2],ppiEdgPth[4*j+3]);
      }
  }
  fclose(fp);
  printf(" EdgPth: %s\n",pcFilNam);

}
void writemommax_(int *iMaxNod, int *iMaxPth, int *iMaxEdg, int *iMaxSrc, char *pchIn){

  char pcFilNam[256],tailname[256];
  FILE *fp;
  int i;
  int slength;
  i = 1;
  strcpy(tailname,".mommax");
//  slength=strlen(*pchIn);
//  strcpy(pcFilNam, pchIn);
//  strcat(pcFilNam, ".mommax");
  fstrcat_(pchIn,tailname,pcFilNam);
  if((fp=fopen(pcFilNam,"w"))==NULL){
    printf("Cannot open the file %s !\n", pcFilNam);
    exit(0);
  }
  fprintf(fp,"%10d %10d %10d %10d\n",*iMaxNod, *iMaxPth,*iMaxEdg,i);
  fclose(fp);
/*
  printf(" Total Nodes   = %10d\n", *iMaxNod);
  printf(" Total Patches = %10d\n", *iMaxPth);
  printf(" Total Edges   = %10d\n", *iMaxEdg);
  printf(" Total SrcEdg  = %10d\n\n",*iMaxSrc);
  printf(" =======Write Geometry Information To Files=====\n");
*/
  printf(" Mommax: %s\n",pcFilNam);
}
