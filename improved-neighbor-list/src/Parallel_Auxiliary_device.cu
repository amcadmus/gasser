#define DEVICE_CODE
#include "Parallel_Auxiliary.h"
#include "common.h"

#include "compile_error_mixcode.h"

// template __global__ void setValue <ScalorType> 
template <typename T>
__global__ void Parallel::Auxiliary::
setValue (T *data, const T value)
{
  data[threadIdx.x + (blockIdx.x + gridDim.x * blockIdx.y) * blockDim.x] = value;
}


template <typename T>
__global__ void Parallel::Auxiliary::
setValue (T *data, unsigned num, const T value)
{
  unsigned ii = threadIdx.x + (blockIdx.x + gridDim.x * blockIdx.y) * blockDim.x;
  if (ii < num){
    data[ii] = value;
  }
}

template <typename T>
__global__ void Parallel::Auxiliary::
copyMatrix (const unsigned  m,
	    const unsigned  n,
	    const T * matA,
	    const unsigned  strideA,
	    T * matB,
	    const unsigned  strideB)
{
  unsigned ii = threadIdx.x + (blockIdx.x + gridDim.x * blockIdx.y) * blockDim.x;
  unsigned shiftA = 0;
  unsigned shiftB = 0;
  if (ii < n){
    for (unsigned jj = 0; jj < m; ++jj){
      matB[shiftB + ii] = matA[shiftA + ii];
      shiftA += strideA;
      shiftB += strideB;
    }
  }
}      


template void Parallel::Auxiliary::setValue (int *, const int);
template void Parallel::Auxiliary::setValue (int3 *, const int3);
template void Parallel::Auxiliary::setValue (int4 *, const int4);
template void Parallel::Auxiliary::setValue (unsigned *, const unsigned);
template void Parallel::Auxiliary::setValue (float *, const float);
template void Parallel::Auxiliary::setValue (float3 *, const float3);
template void Parallel::Auxiliary::setValue (float4 *, const float4);

template void Parallel::Auxiliary::setValue (int *, unsigned,  const int);
template void Parallel::Auxiliary::setValue (int3 *, unsigned,  const int3);
template void Parallel::Auxiliary::setValue (int4 *, unsigned,  const int4);
template void Parallel::Auxiliary::setValue (unsigned *, unsigned,  const unsigned);
template void Parallel::Auxiliary::setValue (float *, unsigned,  const float);
template void Parallel::Auxiliary::setValue (float3 *, unsigned,  const float3);
template void Parallel::Auxiliary::setValue (float4 *, unsigned,  const float4);
