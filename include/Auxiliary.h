#ifndef __Auxiliary_h_wanghan__
#define __Auxiliary_h_wanghan__
#include <stdio.h>

#include "SumVector.h"
#include "common.h"


__device__ void      addKthBit (IndexType *a, IndexType k);
__device__ IndexType getKthBit (IndexType  a, IndexType k);

// __device__ void      sortBlockIndex (IndexType * data, IndexType deepth,
// 				     IndexType * result);

template <typename T>
__device__ void cpyGlobalDataToSharedBuff (T * globalData, T * sharedBuff,
					   IndexType length)
{
  IndexType start = 0;
  do {
    IndexType ii = threadIdx.x + start;
    if (ii < length) {
      sharedBuff[ii] = globalData[ii];
    }
    start += blockDim.x;
  } while (start < length);
  __syncthreads();
}

template <typename T>
__device__ void cpyGlobalDataToSharedBuff (T * globalData, volatile T * sharedBuff,
					   IndexType length)
{
  IndexType start = 0;
  do {
    IndexType ii = threadIdx.x + start;
    if (ii < length) {
      sharedBuff[ii] = globalData[ii];
    }
    start += blockDim.x;
  } while (start < length);
  __syncthreads();
}

template <typename T>
__device__ void setGlobalData (T * globalData, IndexType length, T value)
{
  IndexType start = 0;
  do {
    IndexType ii = threadIdx.x + start;
    if (ii < length) {
      globalData[ii] = value;
    }
    start += blockDim.x;
  } while (start < length);
  __syncthreads();
}

__device__ static IndexType roundUp4 (IndexType x)
{
  if (x & 3 == 0){
    return x;
  }
  else {
    return ((x >> 2) + 1) << 2;
  }
}

template<typename T>
__global__ void rescaleProperty (T * data, IndexType N,
				 const T scalor)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < N) {
    data[ii] *= scalor;
  }
}





__device__ void addKthBit (IndexType * a, IndexType k)
{
  // IndexType tmp = (1 << k);
  // *a += tmp;
  *a += (1 << k);
}

__device__ IndexType getKthBit (IndexType a, IndexType k)
{
  // return a & (1 << k);
  a <<= NUintBit - k - 1;
  return a >> NUintBit - 1;
}


// __device__ void  sortBlockIndex (IndexType * data, 
// 				 IndexType start,
// 				 IndexType deepth,
// 				 IndexType * result)
// {
//   IndexType tid = threadIdx.x;

//   __shared__ IndexType dataBuff  [MaxThreadsPerBlock];
//   __shared__ IndexType bitBuff   [MaxThreadsPerBlock * 2];
//   __shared__ IndexType resultBuff[MaxThreadsPerBlock];
//   __shared__ IndexType count1    [32];

//   dataBuff[tid] = data[tid + start];
//   bitBuff [tid + blockDim.x] = 0;
//   resultBuff[tid] = 0;

//   // printf ("0 tid: %d, \n", tid);
  
//   for (IndexType i = 0; i < deepth; ++i){
//     bitBuff[tid] = getKthBit (dataBuff[tid], i);
//     __syncthreads();
//     IndexType tmpsum = sumVectorBlockBuffer (bitBuff, blockDim.x);
//     if (tid == 0) count1[i] = tmpsum;
//   }
//   // printf ("1 tid: %d, \n", tid);
//   __syncthreads();
//   for (IndexType i = 0; i < deepth; ++i){
//     if (tid >= blockDim.x - count1[i])
//       addKthBit(&resultBuff[tid], i);
//     __syncthreads();
//   }
//   result[tid + start] = resultBuff[tid];
//   printf ("2 tid: %d, %d\n", tid, resultBuff[tid]);
//   __syncthreads();
// }



#endif
