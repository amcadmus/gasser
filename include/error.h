#ifndef __error_h_wanghan___
#define __error_h_wanghan___

#include <stdio.h>
#include "MDException.h"

#ifdef DEVICE_CODE
inline void checkCUDAError (const char *msg)
{
  cudaError_t err = cudaGetLastError();
  if (cudaSuccess != err){
    fprintf (stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString(err));
    // exit (EXIT_FAILURE);
    throw MDExcptCuda();
  }
}
#endif

enum mdError{
  mdSuccess			= 0,
  mdErrorShortCellList		= 1,
  mdErrorShortNeighborList	= 2,
  mdErrorOverFlowCellIdx	= 3,
  mdErrorBreakFENEBond		= 4
};
typedef enum mdError mdError_t;

#endif
