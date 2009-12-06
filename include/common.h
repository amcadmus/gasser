#ifndef __common_h_wanghan__
#define __common_h_wanghan__

#include "systemDefines.h"
#include <stdlib.h>

#ifndef CPP_FILE
#include "error.h"
#endif

#include "MDException.h"

#ifndef CPP_FILE
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

#ifndef CPP_FILE
inline void cudaFreeAPointer (void ** ptr)
{
  if (*ptr != NULL) {
    cudaFree (*ptr);
    *ptr = NULL;
  }
}
#endif 


#endif 
