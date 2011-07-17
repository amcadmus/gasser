#define DEVICE_CODE

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>

#include "common.h"
#include "SumBlock.h"
#include "SumVector.h"
#include "MaxVector.h"
#include "RandomGenerator.h"
#include "Displacement_interface.h"

#define N 128
#define M 200

int main(int argc, char * argv[])
{
  // ScalorType * dresult, hresult;
  // cudaMalloc ((void**)&dresult, sizeof(ScalorType));
  // MaxVector<ScalorType> max;
  // max.reinit (M, N);

  RandomGenerator_MT19937::init_genrand((11));
  // for (unsigned i = 0; i < M; ++i){
  //   max.getBuff()[i] = RandomGenerator_MT19937::genrand_real2();
  // }

  // max.maxBuff(dresult, 0);
  // cudaMemcpy (&hresult, dresult, sizeof(ScalorType), cudaMemcpyDeviceToHost);

  // printf ("max is %f\n", hresult);
  // ScalorType maxValue = max.getBuff()[0];
  // for (unsigned i = 0; i < M; ++i){
  //   if (max.getBuff()[i] > maxValue){
  //     maxValue = max.getBuff()[i];
  //   }
  // }
  // printf ("max is %f\n", maxValue);
  
  // return 0;


  ////////////////////////////////////////////////////////////////////////////////
  // char * filename;
  
  // if (argc != 4){
  //   printf ("Usage:\n%s conf.gro nstep device\n", argv[0]);
  //   return 1;
  // }
  // if (argc != 1){
  //   filename = argv[1];
  // }
  // printf ("# setting device to %d\n", atoi(argv[3]));
  // cudaSetDevice (atoi(argv[3]));
  // checkCUDAError ("set device");

  // MDSystem sys;
  // sys.initConfig (filename);
  // Topology::System sysTop;
  // Topology::Molecule mol;
  // mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  // sysTop.addMolecules (mol, sys.hdata.numAtom);

  // sys.initTopology (sysTop);
  // sys.initDeviceData ();

  // Displacement_max disp(sys, N);
  // disp.recordCoord (sys);

  // ScalorType maxValue = 0;
  // for (unsigned i = 0; i < sys.hdata.numAtom; ++i){
  //   HostVectorType hcoord ;
  //   hcoord.x = sys.ddata.coord[i].x;
  //   hcoord.y = sys.ddata.coord[i].y;
  //   hcoord.z = sys.ddata.coord[i].z;
  //   sys.ddata.coord[i].x += RandomGenerator_MT19937::genrand_real2() - 0.5;
  //   sys.ddata.coord[i].y += RandomGenerator_MT19937::genrand_real2() - 0.5;
  //   sys.ddata.coord[i].z += RandomGenerator_MT19937::genrand_real2() - 0.5;
  //   ScalorType tmp = 0;
  //   tmp += (hcoord.x - sys.ddata.coord[i].x) * (hcoord.x - sys.ddata.coord[i].x);
  //   tmp += (hcoord.y - sys.ddata.coord[i].y) * (hcoord.y - sys.ddata.coord[i].y);
  //   tmp += (hcoord.z - sys.ddata.coord[i].z) * (hcoord.z - sys.ddata.coord[i].z);
  //   tmp = sqrtf(tmp);
  //   if (tmp > maxValue){
  //     maxValue = tmp;
  //   }
  // }

  // ScalorType result = disp.calMaxDisplacemant (sys);
  // printf ("max is %f, %f\n", result, maxValue);
  
  // return 0;


  char * filename;
  
  if (argc != 4){
    printf ("Usage:\n%s conf.gro nstep device\n", argv[0]);
    return 1;
  }
  if (argc != 1){
    filename = argv[1];
  }
  printf ("# setting device to %d\n", atoi(argv[3]));
  cudaSetDevice (atoi(argv[3]));
  checkCUDAError ("set device");

  MDSystem sys;
  sys.initConfig (filename);
  Topology::System sysTop;
  Topology::Molecule mol;
  mol.pushAtom (Topology::Atom (1.0, 0.0, 0));
  sysTop.addMolecules (mol, sys.hdata.numAtom);

  sys.initTopology (sysTop);
  sys.initDeviceData ();

  Displacement_mean disp(sys, N);
  disp.recordCoord (sys);

  double meanValue = 0;
  for (unsigned i = 0; i < sys.hdata.numAtom; ++i){
    HostVectorType hcoord ;
    hcoord.x = sys.ddata.coord[i].x;
    hcoord.y = sys.ddata.coord[i].y;
    hcoord.z = sys.ddata.coord[i].z;
    sys.ddata.coord[i].x += RandomGenerator_MT19937::genrand_real2() - 0.5;
    sys.ddata.coord[i].y += RandomGenerator_MT19937::genrand_real2() - 0.5;
    sys.ddata.coord[i].z += RandomGenerator_MT19937::genrand_real2() - 0.5;
    double tmp = 0;
    tmp += (hcoord.x - sys.ddata.coord[i].x) * (hcoord.x - sys.ddata.coord[i].x);
    tmp += (hcoord.y - sys.ddata.coord[i].y) * (hcoord.y - sys.ddata.coord[i].y);
    tmp += (hcoord.z - sys.ddata.coord[i].z) * (hcoord.z - sys.ddata.coord[i].z);
    tmp = sqrtf(tmp);
    meanValue += tmp;
  }

  ScalorType result = disp.calMeanDisplacemant (sys);
  printf ("mean is %.12f, %.12f\n", result, meanValue/sys.ddata.numAtom);
  
  return 0;

}





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

