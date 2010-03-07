#define DEVICE_CODE

#include "common.h"
#include "GPU_Environment.h"

GPU::Environment::
Environment()
{
}

void GPU::Environment::
setDeviceId (int id)
{
  cudaSetDevice (id);
  checkCUDAError ("set device");
}

