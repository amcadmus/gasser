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
    template <typename T>
    __global__ void
    setValue (T *data, IndexType num, const T value)
    {
      IndexType ii = threadIdx.x + (blockIdx.x + gridDim.x * blockIdx.y) * blockDim.x;
      if (ii < num){
	data[ii] = value;
      }
    }
  }
}
#endif

#endif
