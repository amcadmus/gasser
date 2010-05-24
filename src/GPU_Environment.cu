#define DEVICE_CODE

#include "common.h"
#include "GPU_Environment.h"

#include "compile_error_mixcode.h"

unsigned GPU::Environment::numThreadsInCell = 0;
int GPU::Environment::my_deviceId = -1;
int GPU::Environment::memActiveDevice = 32;
int GPU::Environment::numActiveDevice = 0;
int * GPU::Environment::activeDeviceId = NULL;

void GPU::Environment::
initialize (const unsigned & numThreadsInCell_,
	    const char * activeDeviceName)
{
  numThreadsInCell = numThreadsInCell_;
  
  int numDevice = 0;
  cudaGetDeviceCount(&numDevice);

  if (activeDeviceId == NULL){
    size_t tmp_size = memActiveDevice * sizeof(int);
    activeDeviceId = (int * ) malloc (tmp_size);
    if (activeDeviceId == NULL){
      throw MDExcptFailedMallocOnHost ("GPU::Environment::Environment",
				       "activeDeviceId",
				       tmp_size);
    }
  }
  numActiveDevice = 0;
  for (int i = 0; i < numDevice; ++i){
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, i);
    // printf ("%s\n", deviceProp.name);
    if (strcmp(activeDeviceName, deviceProp.name) == 0){
      if (numActiveDevice == memActiveDevice){
	memActiveDevice *= 2;
	size_t tmp_size = memActiveDevice * sizeof(int);
	realloc (activeDeviceId, tmp_size);
	throw MDExcptFailedReallocOnHost("GPU::Environment::Environment",
					 "activeDeviceId",
					 tmp_size);
      }
      activeDeviceId[numActiveDevice ++] = i;
    }
  }
}

void GPU::Environment::
setDeviceId (int id)
{
  cudaThreadExit() ;
  my_deviceId = id;
  cudaSetDevice (my_deviceId);
  checkCUDAError ("GPU::Environment::set device");
}


void GPU::Environment::
finalize ()
{
  free (activeDeviceId);
  cudaThreadExit() ;
}
