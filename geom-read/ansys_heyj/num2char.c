
#include <stdio.h>
#include <string.h>

//using namespace std;

int Num2Char(int iNum,
    char* pcNum,
    int *iLen){
  int i,j;

  if(iNum<0){
    printf("Wrong input integer %d !",iNum);
    return(-1);
  }
  else if(iNum<10){
    *iLen=1;
    pcNum[0]=iNum+48;
    strcpy((pcNum+(*iLen)),"\0");
  }
  else if(iNum<100){
    *iLen=2;
    i=iNum/10;
    pcNum[0]=i+48;
    pcNum[1]=iNum-i*10+48;
    strcpy((pcNum+(*iLen)),"\0");
  }
  else if(iNum<1000){
    *iLen=3;
    j=iNum/100;
    i=(iNum-j*100)/10;
    pcNum[0]=j+48;
    pcNum[1]=i+48;
    pcNum[2]=iNum-j*100-i*10+48;
    strcpy((pcNum+(*iLen)),"\0");
  }
  else{
    printf("Input number is great than 1000, iNum=%d\n",iNum);
    return(-2);
  }
  return(0);

}



