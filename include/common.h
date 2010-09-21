#ifndef __common_h_wanghan__
#define __common_h_wanghan__

#include <stdlib.h>
#include "systemDefines.h"
#include "error.h"

#include "compile_error_mixcode.h"

#ifdef DEVICE_CODE
static dim3 toGridDim (IndexType nob)
{
  dim3 tmp;
  // if (nob <= 65535){
    tmp.x = nob;
    tmp.y = 1;
    tmp.z = 1;
  // }
  // else{
  //   double nob2=nob;
  //   double root = pow (nob2, 1./3.);
  //   IndexType rootd = IndexType (root + 1e-8);
  //   if (rootd * rootd * rootd != nob){
  //     fprintf(stderr, "toGridDim: inconsistent nob %d %d\n", rootd, nob);
  //     exit (1);
  //   }
  //   tmp.x = rootd * rootd;
  //   tmp.y = rootd;
  //   tmp.z = 1;
  // }
  return tmp;
}
#endif


inline void freeAPointer (void ** ptr)
{
  if (*ptr != NULL) {
    free (*ptr);
    *ptr = NULL;
  }
}


#ifdef DEVICE_CODE
inline void cudaFreeAPointer (void ** ptr)
{
  if (*ptr != NULL) {
    cudaFree (*ptr);
    *ptr = NULL;
  }
}
#endif 


#endif 
