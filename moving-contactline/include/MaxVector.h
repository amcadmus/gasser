#ifndef __MaxVector_h_wanghan__
#define __MaxVector_h_wanghan__

#include <stdio.h>
#include "common.h"
#include "Auxiliary.h"

template <typename SCALORTYPE>
class MaxVector 
{
  bool deviceMalloced;
  IndexType buffLength;
  IndexType NMax;
  IndexType * NBlock;
  dim3 myBlockDim;
  IndexType bitBlockDim;
  IndexType sharedBuffSize;
  void clear();
  void init   (IndexType NumberOfMax,
	       IndexType NThread);
  IndexType roundUpDivide (IndexType v,
			   IndexType bit);
public:
  SCALORTYPE * buff;
public:
  MaxVector();
  ~MaxVector();
  void reinit (IndexType NumberOfMax,
	       IndexType NThread);
  // SCALORTYPE * getBuff () {return buff;}
  void maxBuff (SCALORTYPE * result,
		IndexType posi,
		cudaStream_t stream = 0);
  void maxBuffAdd (SCALORTYPE * result,
		   IndexType posi,
		   cudaStream_t stream = 0);
};



template<typename SCALORTYPE>
IndexType MaxVector<SCALORTYPE>::
roundUpDivide (IndexType v,
	       IndexType bit)
{
  IndexType result = v >> bit;
  if (v == result * (1 << bit)) return result;
  else return result + 1;
}

template <typename SCALORTYPE>
static __global__ void
maxVectorInitSetZero (SCALORTYPE * buff,
		      IndexType length)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  if (bid == 0){
    setGlobalData (buff, length, SCALORTYPE(0));
  }
}

template<typename SCALORTYPE>
MaxVector<SCALORTYPE>::MaxVector()
{
  deviceMalloced = false;
  NBlock = NULL;
  buffLength = 0;
  NMax = 0;
  bitBlockDim = 0;
  sharedBuffSize = 0;
}

template<typename SCALORTYPE>
void MaxVector<SCALORTYPE>::clear()
{
  freeAPointer ((void**)&NBlock);
  // free (NBlock);
  if (deviceMalloced){
    cudaFree (buff);
    deviceMalloced = false;
  }
  buffLength = 0;
  NMax = 0;
  bitBlockDim = 0;
  sharedBuffSize = 0;
}

template<typename SCALORTYPE>
MaxVector<SCALORTYPE>::~MaxVector()
{
  clear();
}

template<typename SCALORTYPE>
void MaxVector<SCALORTYPE>::
init (IndexType NumberOfMax,
      IndexType NThread)
{
  deviceMalloced = false;
  NBlock = NULL;
  buffLength = 0;
  NMax = 0;
  bitBlockDim = 0;
  sharedBuffSize = 0;

  NThread <<= 1;
  bitBlockDim = 1;
  while ((NThread >> bitBlockDim) != 0) bitBlockDim ++;
  bitBlockDim --;
  myBlockDim.x = (1 << bitBlockDim) >> 1;
  // printf ("blockdim is %d, bitblockdim is %d\n", myBlockDim.x, 1<<bitBlockDim);
  
  NMax = 0;
  IndexType tmpNumber = NumberOfMax;
  while ((tmpNumber = roundUpDivide (tmpNumber, bitBlockDim)) != 1) {
    NMax ++;
  }
  NMax ++;
  NBlock = (IndexType *) malloc (sizeof(IndexType) * (NMax));
  buffLength = 0;
  IndexType i = 0;
  tmpNumber = NumberOfMax;
  while ((tmpNumber = roundUpDivide (tmpNumber, bitBlockDim)) != 1) {
    NBlock[i] = tmpNumber;
    buffLength += NBlock[i];
    i++;
  }
  NBlock[i] = tmpNumber;
  buffLength += NBlock[i];
  buffLength <<= bitBlockDim;
  cudaMalloc ((void**)&buff, sizeof(IndexType) * buffLength);
  maxVectorInitSetZero <<<1, myBlockDim>>> (buff, buffLength);
  checkCUDAError ("MaxVector::init");
  deviceMalloced = true;
  sharedBuffSize = (1 << bitBlockDim) * sizeof(SCALORTYPE);
}

template<typename SCALORTYPE>
void MaxVector<SCALORTYPE>::
reinit (IndexType NumberOfMax,
	IndexType NThread)
{
  clear();
  init (NumberOfMax, NThread);
}


extern __shared__ volatile int sbuffpub [];
template<typename SCALORTYPE>
static __global__ void
maxVectorMaxUnit (SCALORTYPE *		buff,
		  IndexType		dataStart,
		  SCALORTYPE *		result,
		  IndexType		resultStart,
		  IndexType		bitBlockDim)
{
  volatile SCALORTYPE * sbuff = (volatile SCALORTYPE *) sbuffpub;

  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType mystart = dataStart + (bid << bitBlockDim) + threadIdx.x;
  sbuff[threadIdx.x] = buff[mystart];
  sbuff[threadIdx.x + blockDim.x] = buff[mystart + blockDim.x];

  IndexType skip = blockDim.x;
  while (skip != 0){
    __syncthreads();
    if (threadIdx.x < skip){
      if (sbuff[threadIdx.x + skip] > sbuff[threadIdx.x]){
	sbuff[threadIdx.x] = sbuff[threadIdx.x + skip];
      }
    }
    skip >>= 1;
  }  
  if (threadIdx.x == 0){
    result[resultStart + bid] = sbuff[0];
  }
}

template<typename SCALORTYPE>
static __global__ void
maxVectorMaxUnitMax (SCALORTYPE *	buff,
		     IndexType		dataStart,
		     SCALORTYPE *	result,
		     IndexType		resultStart,
		     IndexType		bitBlockDim)
{
  volatile SCALORTYPE * sbuff = (volatile SCALORTYPE *) sbuffpub;

  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  // IndexType mystart = dataStart + bid * (blockDim.x << 1) + threadIdx.x;
  IndexType mystart = dataStart + (bid << bitBlockDim) + threadIdx.x;
  sbuff[threadIdx.x] = buff[mystart];
  sbuff[threadIdx.x + blockDim.x] = buff[mystart + blockDim.x];

  IndexType skip = blockDim.x;
  while (skip != 0){
    __syncthreads();
    if (threadIdx.x < skip){
      if (sbuff[threadIdx.x + skip] > sbuff[threadIdx.x]){
	sbuff[threadIdx.x] = sbuff[threadIdx.x + skip];
      }
      // sbuff[threadIdx.x] += sbuff[threadIdx.x + skip];
    }
    skip >>= 1;
  }
  if (threadIdx.x == 0){
    if (sbuff[0] > result[resultStart + bid]){
      result[resultStart + bid] = sbuff[0];
      // result[resultStart + bid] += sbuff[0];
    }
  }
}

template<typename SCALORTYPE>
void MaxVector<SCALORTYPE>::maxBuff (SCALORTYPE * result,
				     IndexType posi,
				     cudaStream_t stream)
{
  IndexType startBuff = 0;
  IndexType startResult = NBlock[0] << bitBlockDim;
  for (IndexType i = 0; i < NMax-1; ++i){
    maxVectorMaxUnit
	<<<NBlock[i], myBlockDim, sharedBuffSize, stream>>> (
	    buff, startBuff,
	    buff, startResult, bitBlockDim);
    startBuff = startResult;
    startResult += NBlock[i+1] << bitBlockDim;
  }
  maxVectorMaxUnit
      <<<NBlock[NMax-1], myBlockDim, sharedBuffSize, stream>>> (
	  buff, startBuff,
	  result, posi, bitBlockDim);
  checkCUDAError ("maxVector::maxBuff");
}

template<typename SCALORTYPE>
void MaxVector<SCALORTYPE>::maxBuffAdd (SCALORTYPE * result,
					IndexType posi,
					cudaStream_t stream)
{
  IndexType startBuff = 0;
  IndexType startResult = NBlock[0] << bitBlockDim;
  for (IndexType i = 0; i < NMax-1; ++i){
    maxVectorMaxUnit
	<<<NBlock[i], myBlockDim, sharedBuffSize, stream>>> (
	    buff, startBuff,
	    buff, startResult, bitBlockDim);
    startBuff = startResult;
    startResult += NBlock[i+1] << bitBlockDim;
  }
  maxVectorMaxUnitMax
      <<<NBlock[NMax-1], myBlockDim, sharedBuffSize, stream>>> (
	  buff, startBuff,
	  result, posi, bitBlockDim);
  checkCUDAError ("maxVector::maxBuff");
}


#endif
