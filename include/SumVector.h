#ifndef __SumVector_h_wanghan__
#define __SumVector_h_wanghan__
#include <stdio.h>
#include "common.h"
#include "Auxiliary.h"

//#define MaxThreadsPerBlocktmptmp  64

template <typename SCALORTYPE>
class SumVector 
{
  bool deviceMalloced;
  IndexType buffLength;
  IndexType NSum;
  IndexType * NBlock;
  dim3 myBlockDim;
  IndexType bitBlockDim;
  IndexType sharedBuffSize;
  void clear();
  void init   (IndexType NumberOfSum, IndexType NThread);
public:
  SCALORTYPE * buff;
public:
  SumVector();
  ~SumVector();
  void reinit (IndexType NumberOfSum, IndexType NThread);
  // SCALORTYPE * getBuff () {return buff;}
  void sumBuff (SCALORTYPE * result, IndexType posi,
		cudaStream_t stream = 0);
  void sumBuffAdd (SCALORTYPE * result, IndexType posi,
		   cudaStream_t stream = 0);
};



static IndexType roundUpDivide (IndexType v, IndexType bit)
{
  IndexType result = v >> bit;
  if (v == result * (1 << bit)) return result;
  else return result + 1;
}

template <typename SCALORTYPE>
static __global__ void
sumVectorInitSetZero (SCALORTYPE * buff, IndexType length)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  if (bid == 0){
    setGlobalData (buff, length, SCALORTYPE(0));
  }
}

template<typename SCALORTYPE>
SumVector<SCALORTYPE>::SumVector()
{
  deviceMalloced = false;
  NBlock = NULL;
  buffLength = 0;
  NSum = 0;
  bitBlockDim = 0;
  sharedBuffSize = 0;
}

template<typename SCALORTYPE>
void SumVector<SCALORTYPE>::clear()
{
  freeAPointer ((void**)&NBlock);
  // free (NBlock);
  if (deviceMalloced){
    cudaFree (buff);
    deviceMalloced = false;
  }
  buffLength = 0;
  NSum = 0;
  bitBlockDim = 0;
  sharedBuffSize = 0;
}

template<typename SCALORTYPE>
SumVector<SCALORTYPE>::~SumVector()
{
  clear();
}

template<typename SCALORTYPE>
void SumVector<SCALORTYPE>::
init (IndexType NumberOfSum, IndexType NThread)
{
  deviceMalloced = false;
  NBlock = NULL;
  buffLength = 0;
  NSum = 0;
  bitBlockDim = 0;
  sharedBuffSize = 0;

  NThread <<= 1;
  bitBlockDim = 1;
  while ((NThread >> bitBlockDim) != 0) bitBlockDim ++;
  bitBlockDim --;
  myBlockDim.x = (1 << bitBlockDim) >> 1;
  // printf ("blockdim is %d, bitblockdim is %d\n", myBlockDim.x, 1<<bitBlockDim);
  
  NSum = 0;
  IndexType tmpNumber = NumberOfSum;
  while ((tmpNumber = roundUpDivide (tmpNumber, bitBlockDim)) != 1) {
    NSum ++;
  }
  NSum ++;
  NBlock = (IndexType *) malloc (sizeof(IndexType) * (NSum));
  buffLength = 0;
  IndexType i = 0;
  tmpNumber = NumberOfSum;
  while ((tmpNumber = roundUpDivide (tmpNumber, bitBlockDim)) != 1) {
    NBlock[i] = tmpNumber;
    buffLength += NBlock[i];
    i++;
  }
  NBlock[i] = tmpNumber;
  buffLength += NBlock[i];
  buffLength <<= bitBlockDim;
  cudaMalloc ((void**)&buff, sizeof(IndexType) * buffLength);
  sumVectorInitSetZero <<<1, myBlockDim>>> (buff, buffLength);
  checkCUDAError ("SumVector::init");
  deviceMalloced = true;
  sharedBuffSize = (1 << bitBlockDim) * sizeof(SCALORTYPE);
}

template<typename SCALORTYPE>
void SumVector<SCALORTYPE>::
reinit (IndexType NumberOfSum, IndexType NThread)
{
  clear();
  init (NumberOfSum, NThread);
}


extern __shared__ volatile int sbuffpub [];
template<typename SCALORTYPE>
static __global__ void sumVectorSumUnit (SCALORTYPE * buff, IndexType dataStart,
				  SCALORTYPE * result, IndexType resultStart,
				  IndexType bitBlockDim)
{
  volatile SCALORTYPE * sbuff = (volatile SCALORTYPE *) sbuffpub;

  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  // IndexType mystart = dataStart + bid * (blockDim.x << 1) + threadIdx.x;
  IndexType mystart = dataStart + (bid << bitBlockDim) + threadIdx.x;
  sbuff[threadIdx.x] = buff[mystart];
  sbuff[threadIdx.x + blockDim.x] = buff[mystart + blockDim.x];

  IndexType skip = blockDim.x;
#pragma unroll 2
  while (skip != 0){
    __syncthreads();
    if (threadIdx.x < skip){
      sbuff[threadIdx.x] += sbuff[threadIdx.x + skip];
    }
    skip >>= 1;
  }

  // if (skip == 32){
  //   if (threadIdx.x < 32){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 32];
  //   }
  //   if (threadIdx.x < 16){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 16];
  //   }
  //   if (threadIdx.x < 8){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 8];
  //   }
  //   if (threadIdx.x < 4){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 4];
  //   }
  //   if (threadIdx.x < 2){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 2];
  //   }
  //   if (threadIdx.x < 1){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 1];
  //   }
  // }
  // else if (skip == 64){
  //   __syncthreads();
  //   if (threadIdx.x < 64){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 64];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 32){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 32];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 16){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 16];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 8){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 8];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 4){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 4];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 2){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 2];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 1){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 1];
  //   }
  //   __syncthreads();
  // }
  // else if (skip == 128){
  //   __syncthreads();
  //   if (threadIdx.x < 128){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 128];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 64){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 64];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 32){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 32];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 16){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 16];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 8){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 8];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 4){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 4];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 2){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 2];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 1){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 1];
  //   }
  //   __syncthreads();
  // }
  // else {
  //   return ;
  // }  
  
  if (threadIdx.x == 0){
    result[resultStart + bid] = sbuff[0];
  }
}

template<typename SCALORTYPE>
static __global__ void sumVectorSumUnitAdd (SCALORTYPE * buff, IndexType dataStart,
				     SCALORTYPE * result, IndexType resultStart,
				     IndexType bitBlockDim)
{
  volatile SCALORTYPE * sbuff = (volatile SCALORTYPE *) sbuffpub;

  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  // IndexType mystart = dataStart + bid * (blockDim.x << 1) + threadIdx.x;
  IndexType mystart = dataStart + (bid << bitBlockDim) + threadIdx.x;
  sbuff[threadIdx.x] = buff[mystart];
  sbuff[threadIdx.x + blockDim.x] = buff[mystart + blockDim.x];

  IndexType skip = blockDim.x;
#pragma unroll 2
  while (skip != 0){
    __syncthreads();
    if (threadIdx.x < skip){
      sbuff[threadIdx.x] += sbuff[threadIdx.x + skip];
    }
    skip >>= 1;
  }

  // if (skip == 32){
  //   if (threadIdx.x < 32){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 32];
  //   }
  //   if (threadIdx.x < 16){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 16];
  //   }
  //   if (threadIdx.x < 8){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 8];
  //   }
  //   if (threadIdx.x < 4){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 4];
  //   }
  //   if (threadIdx.x < 2){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 2];
  //   }
  //   if (threadIdx.x < 1){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 1];
  //   }
  // }
  // else if (skip == 64){
  //   __syncthreads();
  //   if (threadIdx.x < 64){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 64];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 32){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 32];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 16){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 16];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 8){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 8];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 4){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 4];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 2){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 2];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 1){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 1];
  //   }
  //   __syncthreads();
  // }
  // else if (skip == 128){
  //   __syncthreads();
  //   if (threadIdx.x < 128){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 128];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 64){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 64];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 32){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 32];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 16){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 16];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 8){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 8];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 4){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 4];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 2){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 2];
  //   }
  //   __syncthreads();
  //   if (threadIdx.x < 1){
  //     sbuff[threadIdx.x] += sbuff[threadIdx.x + 1];
  //   }
  //   __syncthreads();
  // }
  // else {
  //   return ;
  // }

  if (threadIdx.x == 0){
    result[resultStart + bid] += sbuff[0];
  }
}

template<typename SCALORTYPE>
void SumVector<SCALORTYPE>::sumBuff (SCALORTYPE * result,
				     IndexType posi,
				     cudaStream_t stream)
{
  IndexType startBuff = 0;
  IndexType startResult = NBlock[0] << bitBlockDim;
  for (IndexType i = 0; i < NSum-1; ++i){
    sumVectorSumUnit
	<<<NBlock[i], myBlockDim, sharedBuffSize, stream>>> (
	    buff, startBuff,
	    buff, startResult, bitBlockDim);
    startBuff = startResult;
    startResult += NBlock[i+1] << bitBlockDim;
  }
  sumVectorSumUnit
      <<<NBlock[NSum-1], myBlockDim, sharedBuffSize, stream>>> (
	  buff, startBuff,
	  result, posi, bitBlockDim);
  checkCUDAError ("sumVector::sumBuff");
}

template<typename SCALORTYPE>
void SumVector<SCALORTYPE>::sumBuffAdd (SCALORTYPE * result,
					IndexType posi,
					cudaStream_t stream)
{
  IndexType startBuff = 0;
  IndexType startResult = NBlock[0] << bitBlockDim;
  for (IndexType i = 0; i < NSum-1; ++i){
    sumVectorSumUnit
	<<<NBlock[i], myBlockDim, sharedBuffSize, stream>>> (
	    buff, startBuff,
	    buff, startResult, bitBlockDim);
    startBuff = startResult;
    startResult += NBlock[i+1] << bitBlockDim;
  }
  sumVectorSumUnitAdd
      <<<NBlock[NSum-1], myBlockDim, sharedBuffSize, stream>>> (
	  buff, startBuff,
	  result, posi, bitBlockDim);
  checkCUDAError ("sumVector::sumBuff");
}


// class MaxVector 
// {
//   bool deviceMalloced;
//   ScalorType * buff;
//   IndexType buffLength;
//   IndexType NMax;
//   IndexType * NBlock;
//   dim3 myBlockDim;
//   IndexType bitBlockDim;
//   IndexType sharedBuffSize;
// public:
//   MaxVector();
//   ~MaxVector();
//   void init (IndexType NumberOfMax, IndexType NThread);
//   ScalorType * getBuff () {return buff;}
//   void maxBuff (ScalorType * result, IndexType posi,
// 		cudaStream_t stream = 0);
//   void maxBuffAdd (ScalorType * result, IndexType posi,
// 		   cudaStream_t stream = 0);
// };


/** 
 * Copy partialSum to the coresponding position in *buff and the last
 * block sums all values in *buff
 * 
 * @param partialSum Partial summation of a block
 * @param buff The buffer where each block stores their summation, the size should
 * be at least the same as number of blocks times sizeof(ScalorType).
 * @param counter A counter used by this program, the size should be at least 
 * sizeof(IndexType).
 * @param result Where the result is stored.
 * 
 * @return If this thread gives result, return 1. Otherwise return 0.
 */
static __device__ int sumPartialSum (ScalorType partialSum,
			      ScalorType * buff,
			      IndexType * counter, 
			      ScalorType * result);
/** 
 * Copy partialSum to the coresponding position in *buff and the last
 * block sums all values in *buff
 * 
 * @param partialSum Partial summation of a block
 * @param buff The buffer where each block stores their summation, the size should
 * be at least the same as number of blocks times sizeof(ScalorType).
 * @param counter A counter used by this program, the size should be at least 
 * sizeof(IndexType).
 * @param sharedBuff User provided shared buffer, the size should be at least
 * sizeof(ScalorType) * #threadPerBlock
 * @param result Where the result is stored.
 * 
 * @return If this thread gives result, return 1. Otherwise return 0.
 */
static __device__ int sumPartialSum (ScalorType partialSum,
			      ScalorType * buff,
			      IndexType * counter,
			      ScalorType * sharedBuff,
			      ScalorType * result);





// the number of data stored in pdata should be exactly equal to the
// number of threads in all blocks, the user garentee it.
//
// we provide the initSumVector to produce a suitable buff, the the
// size of which should be larger than the number of blocks and
// aligned to the number of threads in one block.
//
// the thread receive the reture value true will have the right sum in
// result.
static __device__ int sumVector (ScalorType * data,
			  ScalorType * buff, IndexType * counter,
			  ScalorType * result);
static __device__ int sumVector (ScalorType * data, IndexType * N,
			  ScalorType * buff, IndexType * counter,
			  ScalorType * result);



// sum the data stored in pdata, sum NThread items.
static __device__ ScalorType sumVectorBlock (ScalorType * pdata);
// sum the data stored in pdata, sum N items, N should be <= NThread
static __device__ ScalorType sumVectorBlock (ScalorType * pdata, IndexType N);
// the summation will start at start
static __device__ ScalorType sumVectorBlock (ScalorType * data, IndexType start,
				      IndexType N);
static __device__ ScalorType sumVectorBlock (ScalorType * data,
				      IndexType start,
				      IndexType N,
				      volatile ScalorType * buff);
static __device__ ScalorType maxVectorBlock (ScalorType * data, IndexType start,
				      IndexType N);
// summation in given buffer (shared memrory) the size of buffer
// should be larger than half number of threads
static __device__ ScalorType sumVectorBlockBuffer (ScalorType * sharedbuffer, IndexType N);




static __device__ int sumPartialSum (ScalorType partialSum,
			      ScalorType * buff, IndexType * counter, 
			      ScalorType * result)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  __shared__ volatile bool isLastBlockDone;

  if (threadIdx.x == 0) {
    buff[bid] = partialSum;
    IndexType value = atomicInc(counter, gridDim.x*gridDim.y);
    isLastBlockDone = (value == (gridDim.x*gridDim.y - 1));
  }
  __threadfence();
  __syncthreads();
  
  if (isLastBlockDone){
    IndexType p = 0;
    ScalorType tmpsum = 0.f;
    while (p < gridDim.x*gridDim.y){
      IndexType tmp = gridDim.x*gridDim.y - p;
      IndexType n = ((tmp < blockDim.x) ? tmp : blockDim.x);
      tmpsum += sumVectorBlock (buff, p, n);
      p += blockDim.x;
    }
    if (threadIdx.x == 0){
      *counter = 0;
      *result = tmpsum;
      return 1;
    }
  }
  return 0;
}

static __device__ int
sumPartialSum (ScalorType partialSum,
	       ScalorType * buff,
	       IndexType * counter,
	       volatile ScalorType * sharedBuff,
	       ScalorType * result)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  __shared__ volatile bool isLastBlockDone;

  if (threadIdx.x == 0) {
    buff[bid] = partialSum;
    IndexType value = atomicInc(counter, gridDim.x*gridDim.y);
    isLastBlockDone = (value == (gridDim.x*gridDim.y - 1));
  }
  __threadfence();
  __syncthreads();
  
  if (isLastBlockDone){
    IndexType p = 0;
    ScalorType tmpsum = 0.f;
    while (p < gridDim.x*gridDim.y){
      IndexType tmp = gridDim.x*gridDim.y - p;
      IndexType n = ((tmp < blockDim.x) ? tmp : blockDim.x);
      tmpsum += sumVectorBlock (buff, p, n, sharedBuff);
      p += blockDim.x;
    }
    if (threadIdx.x == 0){
      *counter = 0;
      *result = tmpsum;
      return 1;
    }
  }
  return 0;
}


static __device__ int
maxPartialMax (ScalorType partialMax,
	       ScalorType * buff, IndexType * counter, 
	       ScalorType * result)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  __shared__ volatile bool isLastBlockDone;

  if (threadIdx.x == 0) {
    buff[bid] = partialMax;
    IndexType value = atomicInc(counter, gridDim.x*gridDim.y);
    isLastBlockDone = (value == (gridDim.x*gridDim.y - 1));
  }
  __threadfence();
  __syncthreads();
  
  if (isLastBlockDone){
    IndexType p = 0;
    ScalorType tmpmax = 0;
    while (p < gridDim.x*gridDim.y){
      IndexType tmp = gridDim.x*gridDim.y - p;
      IndexType n = ((tmp < blockDim.x) ? tmp : blockDim.x);
      tmpmax = fmax(tmpmax, maxVectorBlock (buff, p, n));
      p += blockDim.x;
    }
    if (threadIdx.x == 0){
      *counter = 0;
      *result = tmpmax;
      return 1;
    }
  }
  return 0;
}
    
    
  

__device__ int sumVector (ScalorType * data,
			  ScalorType * buff, IndexType * counter,
			  ScalorType * result)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;  
  ScalorType partialSum = sumVectorBlock (data);
  __shared__ volatile bool isLastBlockDone_for_sumVector;

  if (tid == 0){
    buff[bid] = partialSum;
    IndexType value = atomicInc(counter, gridDim.x*gridDim.y);
    isLastBlockDone_for_sumVector = (value == (gridDim.x*gridDim.y - 1));
  }
  __threadfence();
  __syncthreads();

  if (isLastBlockDone_for_sumVector){
    IndexType p = 0;
    ScalorType tmpsum = 0.f;
    while (p < gridDim.x*gridDim.y){
      IndexType tmp = gridDim.x*gridDim.y - p;
      IndexType n = ((tmp < blockDim.x) ? tmp : blockDim.x);
      // printf ("%d\n", blockDim.x);
      tmpsum += sumVectorBlock (buff, p, n);
      p += blockDim.x;
    }
    if (tid == 0){
      *counter = 0;
      *result = tmpsum;
      return 1;
    }
  }
  return 0;
}

static __device__ int
sumVector (ScalorType * data, IndexType N, 
	   ScalorType * buff, IndexType * counter,
	   ScalorType * result)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;  
  ScalorType partialSum;
  IndexType value; 
  IntScalorType tmp  = IntScalorType((bid+1)*blockDim.x) - N ;
  if (tmp >= IntScalorType(blockDim.x)){
    partialSum = 0;
  }
  else if (tmp > 0){
    IndexType tmpN = blockDim.x - tmp;
    partialSum = sumVectorBlock (data, tmpN);
  }
  else {
    partialSum = sumVectorBlock (data);
  }

  __shared__ volatile bool isLastBlockDone_for_sumVector;

  // printf ("bid:%d, tid:%d, counter: %d\n", bid, tid, *counter);
  if (tid == 0){
    buff[bid] = partialSum;
    value = atomicInc(counter, gridDim.x*gridDim.y);
    isLastBlockDone_for_sumVector = (value == (gridDim.x*gridDim.y - 1));
  }
  __threadfence();
  __syncthreads();

  // if (tid == 0)
  //   printf ("bid:%d, tid:%d, value: %d\n", bid, tid, value);
  if (isLastBlockDone_for_sumVector){
    IndexType p = 0;
    ScalorType tmpsum = 0.f;
    while (p < gridDim.x*gridDim.y){
      IndexType tmp = gridDim.x*gridDim.y - p;
      IndexType n = ((tmp < blockDim.x) ? tmp : blockDim.x);
      // printf ("%d\n", blockDim.x);
      tmpsum += sumVectorBlock (buff, p, n);
      p += blockDim.x;
    }
    if (tid == 0){
      *counter = 0;
      *result = tmpsum;
      return 1;
    }
  }
  return 0;
}


__device__ ScalorType sumVectorBlock (ScalorType * data)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType num = blockDim.x ;
  IndexType skip = 1;
  __shared__ volatile ScalorType buff [MaxThreadsPerBlock * 2];
  buff[tid] = 0;
  buff[tid + blockDim.x] = 0;
  buff[tid] = data[tid + bid * blockDim.x];
  __syncthreads();
  
  while (num != 1){
    ScalorType tmp = buff[tid] + buff[tid+skip];
    __syncthreads();
    buff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return buff[0];
}

__device__ ScalorType sumVectorBlock (ScalorType * data, IndexType N)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType num = N;
  IndexType skip = 1;
  __shared__ volatile ScalorType buff [MaxThreadsPerBlock * 2];
  buff[tid] = 0;
  buff[tid + blockDim.x] = 0;
  if (tid < num)
    buff[tid] = data[tid + bid * blockDim.x];
  __syncthreads();
  
  while (num != 1){
    ScalorType tmp = buff[tid] + buff[tid+skip];
    __syncthreads();
    buff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return buff[0];
}


__device__ ScalorType sumVectorBlock (ScalorType * data, IndexType start,
				      IndexType N)
{
  IndexType tid = threadIdx.x;
  IndexType num = N;
  IndexType skip = 1;
  __shared__ volatile ScalorType buff [MaxThreadsPerBlock * 2];
  buff[tid] = 0;
  buff[tid + blockDim.x] = 0;
  if (tid < num)
    buff[tid] = data[tid + start];
  __syncthreads();
  
  while (num != 1){
    ScalorType tmp = buff[tid] + buff[tid+skip];
    __syncthreads();
    buff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return buff[0];
}

static __device__ ScalorType
sumVectorBlock (volatile ScalorType * data, IndexType start,
		IndexType N)
{
  IndexType tid = threadIdx.x;
  IndexType num = N;
  IndexType skip = 1;
  __shared__ volatile ScalorType buff [MaxThreadsPerBlock * 2];
  buff[tid] = 0;
  buff[tid + blockDim.x] = 0;
  if (tid < num)
    buff[tid] = data[tid + start];
  __syncthreads();
  
  while (num != 1){
    ScalorType tmp = buff[tid] + buff[tid+skip];
    __syncthreads();
    buff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return buff[0];
}

__device__ ScalorType sumVectorBlock (ScalorType * data,
				      IndexType start,
				      IndexType N,
				      volatile ScalorType * buff)
{
  IndexType tid = threadIdx.x;
  IndexType num = N;
  IndexType skip = 1;
  buff[tid] = 0;
  buff[tid + blockDim.x] = 0;
  if (tid < num)
    buff[tid] = data[tid + start];
  __syncthreads();
  
  while (num != 1){
    ScalorType tmp = buff[tid] + buff[tid+skip];
    __syncthreads();
    buff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return buff[0];
}


__device__ ScalorType maxVectorBlock (ScalorType * data, IndexType start,
				      IndexType N)
{
  IndexType tid = threadIdx.x;
  IndexType num = N;
  IndexType skip = 1;
  __shared__ volatile ScalorType buff [MaxThreadsPerBlock * 2];
  buff[tid] = 0;
  buff[tid + blockDim.x] = 0;
  if (tid < num)
    buff[tid] = data[tid + start];
  __syncthreads();
  
  while (num != 1){
    ScalorType tmp = fmaxf(buff[tid],  buff[tid+skip]);
    __syncthreads();
    buff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return buff[0];
}



static __device__ IndexType
maxVectorBlock (IndexType * data, IndexType start,
		IndexType N)
{
  IndexType tid = threadIdx.x;
  IndexType num = N;
  IndexType skip = 1;
  __shared__ volatile IndexType buff [MaxThreadsPerBlock * 2];
  buff[tid] = 0;
  buff[tid + blockDim.x] = 0;
  if (tid < num)
    buff[tid] = data[tid + start];
  __syncthreads();
  
  while (num != 1){
    IndexType tmp = (buff[tid] >  buff[tid+skip]) ? buff[tid] : buff[tid+skip];
    __syncthreads();
    buff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return buff[0];
}




__device__ ScalorType sumVectorBlockBuffer (ScalorType * sharedbuff, IndexType N)
{
  IndexType tid = threadIdx.x;
  IndexType num = N;
  IndexType skip = 1;
  
  while (num != 1){
    ScalorType tmp = sharedbuff[tid] + sharedbuff[tid+skip];
    __syncthreads();
    sharedbuff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return sharedbuff[0];
}

static __device__ ScalorType
sumVectorBlockBuffer (volatile ScalorType * sharedbuff, IndexType N)
{
  IndexType tid = threadIdx.x;
  IndexType num = N;
  IndexType skip = 1;
  
  while (num != 1){
    ScalorType tmp = sharedbuff[tid] + sharedbuff[tid+skip];
    __syncthreads();
    sharedbuff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return sharedbuff[0];
}

template <typename T>
__device__ void sumVectorBlockBuffer_2 (volatile T * sharedbuff)
{
  IndexType skip = blockDim.x;
  IndexType presentBond = skip;
  if (skip & 1) skip ++; 
  skip >>= 1;
  while (skip != 1) {
    __syncthreads();
    if (threadIdx.x < skip && threadIdx.x + skip < presentBond){
      sharedbuff[threadIdx.x] += sharedbuff[threadIdx.x + skip];
    }
    presentBond = skip;
    if (skip & 1) skip ++;
    skip >>= 1;
  }
  __syncthreads();
  if (threadIdx.x == 0) sharedbuff[0] += sharedbuff[1];
}

template <typename T>
__device__ void maxVectorBlockBuffer_2 (volatile T * sharedbuff)
{
  IndexType skip = blockDim.x;
  IndexType presentBond = skip;
  if (skip & 1) skip ++; 
  skip >>= 1;
  while (skip != 1) {
    __syncthreads();
    if (threadIdx.x < skip && threadIdx.x + skip < presentBond){
      if (sharedbuff[threadIdx.x + skip] > sharedbuff[threadIdx.x]){
	sharedbuff[threadIdx.x] = sharedbuff[threadIdx.x + skip];
      }
    }
    presentBond = skip;
    if (skip & 1) skip ++;
    skip >>= 1;
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    if (sharedbuff[1] > sharedbuff[0]){
      sharedbuff[0] = sharedbuff[1];
    }
  }
}


static __device__ ScalorType
maxVectorBlockBuffer (ScalorType * sharedbuff, IndexType N)
{
  IndexType tid = threadIdx.x;
  IndexType num = N;
  IndexType skip = 1;
  
  while (num != 1){
    ScalorType tmp = fmaxf(sharedbuff[tid], sharedbuff[tid+skip]);
    __syncthreads();
    sharedbuff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return sharedbuff[0];
}

static __device__ ScalorType
maxVectorBlockBuffer (volatile ScalorType * sharedbuff, IndexType N)
{
  IndexType tid = threadIdx.x;
  IndexType num = N;
  IndexType skip = 1;
  
  while (num != 1){
    ScalorType tmp = fmaxf(sharedbuff[tid], sharedbuff[tid+skip]);
    __syncthreads();
    sharedbuff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return sharedbuff[0];
}

static __device__ IndexType
sumVectorBlockBuffer (IndexType * sharedbuff, IndexType N)
{
  IndexType tid = threadIdx.x;
  IndexType num = N;
  IndexType skip = 1;
  
  while (num != 1){
    IndexType tmp = sharedbuff[tid] + sharedbuff[tid+skip];
    __syncthreads();
    sharedbuff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return sharedbuff[0];
}

static __device__ IndexType
sumVectorBlockBuffer (volatile IndexType * sharedbuff, IndexType N)
{
  IndexType tid = threadIdx.x;
  IndexType num = N;
  IndexType skip = 1;
  
  while (num != 1){
    IndexType tmp = sharedbuff[tid] + sharedbuff[tid+skip];
    __syncthreads();
    sharedbuff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return sharedbuff[0];
}

static __device__ IndexType
maxVectorBlockBuffer (IndexType * sharedbuff, IndexType N)
{
  IndexType tid = threadIdx.x;
  IndexType num = N;
  IndexType skip = 1;
  
  while (num != 1){
    IndexType tmp = (sharedbuff[tid] > sharedbuff[tid+skip]) ?
	sharedbuff[tid] : sharedbuff[tid + skip];
    __syncthreads();
    sharedbuff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return sharedbuff[0];
}

static __device__ IndexType
maxVectorBlockBuffer (volatile IndexType * sharedbuff, IndexType N)
{
  IndexType tid = threadIdx.x;
  IndexType num = N;
  IndexType skip = 1;
  
  while (num != 1){
    IndexType tmp = (sharedbuff[tid] > sharedbuff[tid+skip]) ?
	sharedbuff[tid] : sharedbuff[tid + skip];
    __syncthreads();
    sharedbuff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return sharedbuff[0];
}



static __device__ ScalorType
sumVectorBlockMomentum (ScalorType * mass,
			ScalorType * velo,
			IndexType N)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType num = N;
  IndexType skip = 1;
  __shared__ volatile ScalorType buff [MaxThreadsPerBlock * 2];
  buff[tid] = 0;
  buff[tid + blockDim.x] = 0;
  if (tid < num)
    buff[tid] = mass[tid + bid * blockDim.x] * velo[tid + bid * blockDim.x];
  __syncthreads();
  
  while (num != 1){
    ScalorType tmp = buff[tid] + buff[tid+skip];
    __syncthreads();
    buff[tid] = tmp;
    num += 1;
    num >>= 1;
    skip <<= 1;
    __syncthreads();
  }
  return buff[0];
}


static __device__ int
sumVectorMomentum (ScalorType * mass,
		   ScalorType * velo,
		   IndexType N, 
		   ScalorType * buff, IndexType * counter,
		   ScalorType * result)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;  
  ScalorType partialSum;
  IndexType value; 
  IntScalorType tmp  = IntScalorType((bid+1)*blockDim.x) - N ;
  if (tmp >= IntScalorType(blockDim.x)){
    partialSum = 0;
  }
  else if (tmp > 0){
    IndexType tmpN = blockDim.x - tmp;
    partialSum = sumVectorBlockMomentum (mass, velo, tmpN);
  }
  else {
    partialSum = sumVectorBlockMomentum (mass, velo, blockDim.x);
  }

  __shared__ volatile bool isLastBlockDone_for_sumVector;

  // printf ("bid:%d, tid:%d, counter: %d\n", bid, tid, *counter);
  if (tid == 0){
    buff[bid] = partialSum;
    value = atomicInc(counter, gridDim.x*gridDim.y);
    isLastBlockDone_for_sumVector = (value == (gridDim.x*gridDim.y - 1));
  }
  __threadfence();
  __syncthreads();

  // if (tid == 0)
  //   printf ("bid:%d, tid:%d, value: %d\n", bid, tid, value);
  if (isLastBlockDone_for_sumVector){
    IndexType p = 0;
    ScalorType tmpsum = 0.f;
    while (p < gridDim.x*gridDim.y){
      IndexType tmp = gridDim.x*gridDim.y - p;
      IndexType n = ((tmp < blockDim.x) ? tmp : blockDim.x);
      // printf ("%d\n", blockDim.x);
      tmpsum += sumVectorBlock (buff, p, n);
      p += blockDim.x;
    }
    if (tid == 0){
      *counter = 0;
      *result = tmpsum;
      return 1;
    }
  }
  return 0;
}
















#endif
