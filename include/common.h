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
  tmp.x = nob;
  tmp.y = 1;
  tmp.z = 1;
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
