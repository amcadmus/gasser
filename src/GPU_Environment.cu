#define DEVICE_CODE

#include "GPU_Environment.h"
#include "error.h"

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

