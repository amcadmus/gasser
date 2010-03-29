#ifndef __Parallel_Auxiliary_h_wanghan__
#define __Parallel_Auxiliary_h_wanghan__


#ifdef DEVICE_CODE
namespace Parallel {
  namespace Auxiliary {
    template <typename T>
    __global__ void
    setValue (T *data, const T value)
    {
      data[threadIdx.x + (blockIdx.x + gridDim.x * blockIdx.y) * blockDim.x] = value;
    }
  }
}
#endif

#endif
