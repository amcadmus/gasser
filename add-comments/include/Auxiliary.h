#ifndef __Auxiliary_h_wanghan__
#define __Auxiliary_h_wanghan__
#include <stdio.h>

#include "common.h"
#include "SumVector.h"

static __device__ void      addKthBit (IndexType *a, IndexType k);
static __device__ IndexType getKthBit (IndexType  a, IndexType k);

// __device__ void      sortBlockIndex (IndexType * data, IndexType deepth,
// 				     IndexType * result);

template <typename T>
__device__ void
cpyGlobalDataToSharedBuff (const T * globalData,
			   T * sharedBuff,
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

// template <typename T>
// __device__ void
// cpyGlobalDataToSharedBuff (const T * globalData,
// 			   volatile T * sharedBuff,
// 			   IndexType length)
// {
//   IndexType start = 0;
//   do {
//     IndexType ii = threadIdx.x + start;
//     if (ii < length) {
//       sharedBuff[ii] = globalData[ii];
//     }
//     start += blockDim.x;
//   } while (start < length);
//   __syncthreads();
// }


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
__global__ void
setValue (T * data,
	  const IndexType N,
	  const T value)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < N) {
    data[ii] = value;
  }
}


template<typename T>
__global__ void
rescaleProperty (T * data,
		 const IndexType N,
		 const T scalor)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < N) {
    data[ii] *= scalor;
  }
}

template<typename T>
__global__ void
rescaleProperty (T * data,
		 const IndexType start,
		 const IndexType N,
		 const T scalor)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < N) {
    data[start+ii] *= scalor;
  }
}
  

#ifdef COORD_IN_ONE_VEC
static __global__ void
rescaleCoord (CoordType * data,
	      const IndexType N,
	      const CoordType scalor)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < N) {
    data[ii].x *= scalor.x;
    data[ii].y *= scalor.y;
    data[ii].z *= scalor.z;
  }
}
#endif


template<typename T>
__global__ void
cpyProperty (T *		to,
	     const T *		from,
	     const IndexType	numAtom)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    to[ii] = from[ii];
  }
}


template<typename T>
__global__ void
addProperty (T *		to,
	     const T *		from,
	     const IndexType	numAtom)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    to[ii] += from[ii];
  }
}



__device__ void
addKthBit (IndexType * a, IndexType k)
{
  // IndexType tmp = (1 << k);
  // *a += tmp;
  *a += (1 << k);
}

__device__ IndexType
getKthBit (IndexType a, IndexType k)
{
  // return a & (1 << k);
  a <<= NUintBit - k - 1;
  return a >> NUintBit - 1;
}


#endif
