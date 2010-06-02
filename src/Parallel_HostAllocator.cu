#define DEVICE_CODE

#include "common.h"
#include "Parallel_HostAllocator.h"
#include "MDException.h"
#include "error.h"

#include "compile_error_mixcode.h"

void Parallel::HostAllocator::
hostMalloc (void** ptr,
	    size_t size,
	    HostMallocType_t type)
{
  cudaError_t error;
  
  switch (type){
  case hostMallocCDefault:
      *ptr = malloc (size);
      if (*ptr == NULL){
	throw MDExcptFailedMallocOnHost ("Parallel::HostAllocator::",
					 "ptr", size);
      }
      break;
  case hostMallocPageLocked:
      error = cudaHostAlloc (ptr, size, cudaHostAllocDefault);
      if (error != cudaSuccess){
	throw MDExcptFailedMallocOnHost ("Parallel::HostAllocator::",
					 "ptr (page-locked)", size);
      }
      break;
  default:
      break;
  }
}

void Parallel::HostAllocator::
hostFree   (void** ptr,
	    HostMallocType_t type)
{
  // cudaError_t error;

  switch (type){
  case hostMallocCDefault:
      free (*ptr);
      *ptr = NULL;
      break;
  case hostMallocPageLocked:
      cudaFreeHost (*ptr);
      *ptr = NULL;
      checkCUDAError ("Parallel::HostAllocator::hostFree");
      break;
  default:
      break;
  }
}

