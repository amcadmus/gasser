// #include <assert.h>
// #include <ctype.h>
// #include <errno.h>
// #include <limits.h>
// #include <string.h>
// #include <stdarg.h>
// #include <stdlib.h>
// #include <stdio.h>

// #include "common.h"
// #include "SumBlock.h"
// #include "SumVector.h"

// #define N 128

// __global__ void compare1 (ScalorType * a, ScalorType * result)
// {
//   if (blockDim.x != N) return ;
//   __shared__ ScalorType buff[N*2];
//   buff[threadIdx.x] = a[threadIdx.x];
//   sumVectorBlockBuffer (buff, N);
//   if (threadIdx.x == 0 && blockIdx.x == 0) *result = buff[0];
// }

// __global__ void compare2 (ScalorType * a, ScalorType * result)
// {
//   if (blockDim.x != N) return ;
//   __shared__ ScalorType buff[N];
//   buff[threadIdx.x] = a[threadIdx.x];
//   sumVectorBlockBuffer_2 (buff);
//   if (threadIdx.x == 0 && blockIdx.x == 0) *result = buff[0];
// }


// __global__ void sum32 (ScalorType * a, ScalorType * result)
// {
//   if (blockDim.x != 32) return ;
//   __shared__ ScalorType buff[32];
//   buff[threadIdx.x] = a[threadIdx.x];
//   SumBlock::sum32_1bsize (buff);
//   if (threadIdx.x == 0 && blockIdx.x == 0) *result = buff[0];
// }

// __global__ void sum64 (ScalorType * a, ScalorType * result)
// {
//   if (blockDim.x != 64) return ;
//   __shared__ ScalorType buff[64];
//   buff[threadIdx.x] = a[threadIdx.x];
//   SumBlock::sum64_1bsize (buff);
//   if (threadIdx.x == 0 && blockIdx.x == 0) *result = buff[0];
// }

// __global__ void sum128 (ScalorType * a, ScalorType * result)
// {
//   if (blockDim.x != 128) return ;
//   __shared__ ScalorType buff[128];
//   buff[threadIdx.x] = a[threadIdx.x];
//   SumBlock::sum128_1bsize (buff);
//   if (threadIdx.x == 0 && blockIdx.x == 0) *result = buff[0];
// }


// __global__ void sum32_2b (ScalorType * a, ScalorType * result)
// {
//   if (blockDim.x != 32) return ;
//   __shared__ ScalorType buff[32*2];
//   buff[threadIdx.x] = a[threadIdx.x];
//   SumBlock::sum32_2bsize (buff);
//   if (threadIdx.x == 0 && blockIdx.x == 0) *result = buff[0];
// }


// #define NBlock 3000
// #define NTime 100

// int main(int argc, char * argv[])
// {
//   int deviceNum = 1;
//   printf ("# setting device to %d\n", deviceNum);
//   cudaSetDevice (deviceNum);
//   checkCUDAError ("set device");

  
//   ScalorType hdata [N];
//   for (IndexType i = 0; i < N; ++i){
//     hdata[i] = i;
//   }
//   ScalorType * ddata;
//   cudaMalloc ((void**)&ddata, sizeof(ScalorType) * N);
//   cudaMemcpy (ddata, hdata, sizeof(ScalorType)*N, cudaMemcpyHostToDevice);

//   ScalorType hresule;
//   ScalorType * dresult;
//   cudaMalloc ((void**)&dresult, sizeof(ScalorType));

//   cudaEvent_t start, stop;
//   float tmptime = 0;
//   cudaEventCreate (&start);
//   cudaEventCreate (&stop );
//   cudaEventRecord(start, 0);
//   for (unsigned i = 0; i < NTime; ++i){
//     sum128 <<<NBlock, N>>> (ddata, dresult);
//   }
//   cudaEventRecord(stop, 0);
//   cudaEventSynchronize (stop);
//   cudaEventElapsedTime (&tmptime, start, stop);

//   cudaMemcpy (&hresule, dresult, sizeof(ScalorType), cudaMemcpyDeviceToHost);
//   printf ("result: %f, time %.3e\n", hresule, tmptime);

    
//   cudaFree(dresult);
//   cudaEventDestroy (start);
//   cudaEventDestroy (stop);

//   return 0;
// }

int main(int argc, char * argv[])
{
    
}
