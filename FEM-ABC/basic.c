
/* ------------------------------------------------------------
   This is the file "basic.c" of the H2Lib package.
   All rights reserved, Steffen Boerm 2009
   ------------------------------------------------------------ */

#include "basic.h"
#include "string.h"
#include <stdio.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif


/* ------------------------------------------------------------
   Memory management
   ------------------------------------------------------------ */

void     *
_h2_allocmem(size_t sz)
{
  void     *ptr;
  ptr = malloc(sz);

  int i;

  memset(ptr,0,sz);

  if (ptr == NULL && sz > 0) {
    printf("Memory allocation of %lu bytes failed \n",
		   (unsigned long) sz);
    abort();
  }

  return ptr;
}

int     *
_h2_allocint(size_t sz)
{
  int     *ptr;
  size_t    dsz;
  int i;
  dsz = sizeof(int) * sz;
  if (dsz / sizeof(int) != sz) {
   printf("Integer overflow in vector allocation\n");
    abort();
  }

  ptr = (int *) malloc(dsz);

  for(i=0;i<sz;i++)
   ptr[i]=0;


  if (ptr == NULL && dsz > 0) {
   printf("Vector allocation of %lu entries failed \n",(unsigned long) sz);
    abort();
  }

  return ptr;
}


real     *
_h2_allocreal(size_t sz)
{
  real     *ptr;
  size_t    dsz;
  int i;

  dsz = sizeof(real) * sz;
  if (dsz / sizeof(real) != sz) {
    printf("Integer overflow in vector allocation\n");
    abort();
  }

  ptr = (real *) malloc(dsz);

  for(i=0;i<sz;i++)
   ptr[i]=0.0;  

  if (ptr == NULL && dsz > 0) {
  printf("Vector allocation of %lu entries failed\n",
		   (unsigned long) sz);
    abort();
  }

  return ptr;
}

field    *
_h2_allocfield(size_t sz)
{
  field    *ptr;
  size_t    dsz;
  int i;

  dsz = sizeof(field) * sz;
  if (dsz / sizeof(field) != sz) {
    printf("Integer overflow in vector allocation \n");
    abort();
  }

  ptr = (field *) malloc(dsz);

  for(i=0;i<sz;i++)
   ptr[i]=0.0+0.0*I;

  if (ptr == NULL && dsz > 0) {
   printf("Vector allocation of %lu entries failed \n",(unsigned long) sz);
    abort();
  }

  return ptr;
}

field    *
_h2_allocmatrix(size_t rows, size_t cols)
{
  field    *ptr;
  size_t    dsz;
  int i,sz;

  sz=rows*cols;
  dsz = sizeof(field) * rows * cols;
  if (dsz / sizeof(field) != rows * cols) {
    printf("Integer overflow in matrix allocation");
    abort();
  }

  ptr = (field *) malloc(dsz);

  for(i=0;i<sz;i++)
   ptr[i]=0.0+I*0.0;


  if (ptr == NULL && dsz > 0) {
    printf("Matrix allocation with %lu rows and %lu columns failed \n",
		   (unsigned long) rows, (unsigned long) cols);
    abort();
  }

  return ptr;
}

void
freemem(void *ptr)
{
  free(ptr);
}
