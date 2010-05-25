#ifndef __Parallel_Auxiliary_h_wanghan__
#define __Parallel_Auxiliary_h_wanghan__

namespace Parallel{
  enum HostMallocType {
    hostMallocDefault			= 1,
    hostMallocPageLocked		= 2
  };
  typedef enum HostMallocType HostMallocType_t;
  
  namespace HostAllocator{
    void hostMalloc (void** ptr, size_t size, HostMallocType_t type);
    void hostFree   (void** ptr, HostMallocType_t type);
  }
}


#ifdef DEVICE_CODE
namespace Parallel {
  namespace Auxiliary {
    template <typename T>
    extern __global__ void setValue (T *data, const T value);
    template <typename T>
    extern __global__ void setValue (T *data, unsigned num, const T value);
    template <typename T>
    extern __global__ void copyMatrix (const unsigned  m,
				       const unsigned  n,
				       const T * matA,
				       const unsigned  strideA,
				       T * matB,
				       const unsigned  strideB);
    
  }
}

// template <typename T>
// static __global__ void Parallel::Auxiliary::
// setValue (T *data, const T value)
// {
//   data[threadIdx.x + (blockIdx.x + gridDim.x * blockIdx.y) * blockDim.x] = value;
// }


// template <typename T>
// static __global__ void Parallel::Auxiliary::
// setValue (T *data, unsigned num, const T value)
// {
//   unsigned ii = threadIdx.x + (blockIdx.x + gridDim.x * blockIdx.y) * blockDim.x;
//   if (ii < num){
//     data[ii] = value;
//   }
// }

// template <typename T>
// static __global__ void Parallel::Auxiliary::
// copyMatrix (const unsigned  m,
// 	    const unsigned  n,
// 	    const T * matA,
// 	    const unsigned  strideA,
// 	    T * matB,
// 	    const unsigned  strideB)
// {
//   unsigned ii = threadIdx.x + (blockIdx.x + gridDim.x * blockIdx.y) * blockDim.x;
//   unsigned shiftA = 0;
//   unsigned shiftB = 0;
//   if (ii < n){
//     for (unsigned jj = 0; jj < m; ++jj){
//       matB[shiftB + ii] = matA[shiftA + ii];
//       shiftA += strideA;
//       shiftB += strideB;
//     }
//   }
// }      


#endif

#endif
