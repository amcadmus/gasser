#include "Statistic.h"


// __global__ void clearStatisticData	(StatisticData sdata)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   if (threadIdx.x + bid == 0)
//     for (IndexType i = 0; i < Size_StatisticData; ++i){
//       sdata.data[i] = 0;
//     }
// }


// __host__ void initHostStatistic (HostStatistic * hst)
// {
//   hst->bondedP = 0;
//   hst->nonBondedP = 0;
//   hst->electrostaticP = 0;
//   hst->totalEnergy = 0;
//   hst->kinetic = 0;
//   hst->virial = 0;
//   hst->temperature = 0;
//   hst->pressure = 0;
// }

// __host__ void initDeviceStatistic (DeviceStatistic *dst)
// {
//   size_t size = sizeof(ScalorType);
//   cudaMalloc ((void**)&(dst->bondedP), size);
//   cudaMalloc ((void**)&(dst->nonBondedP), size);
//   cudaMalloc ((void**)&(dst->electrostaticP), size);
//   cudaMalloc ((void**)&(dst->totalEnergy), size);
//   cudaMalloc ((void**)&(dst->kinetic), size);
//   cudaMalloc ((void**)&(dst->virial), size);
//   cudaMalloc ((void**)&(dst->temperature), size);
//   cudaMalloc ((void**)&(dst->pressure), size);
// }

// __host__ void cpyDeviceStatisticToHost (const DeviceStatistic * dst,
// 					HostStatistic *hst)
// {
//   size_t size = sizeof(ScalorType);
//   cudaMemcpy (dst->bondedP, &(hst->bondedP), size, cudaMemcpyHostToDevice);
//   cudaMemcpy (dst->nonBondedP, &(hst->nonBondedP), size, cudaMemcpyHostToDevice);
//   cudaMemcpy (dst->electrostaticP, &(hst->electrostaticP), size, cudaMemcpyHostToDevice);
//   cudaMemcpy (dst->totalEnergy, &(hst->totalEnergy), size, cudaMemcpyHostToDevice);
//   cudaMemcpy (dst->kinetic, &(hst->kinetic), size, cudaMemcpyHostToDevice);
//   cudaMemcpy (dst->virial, &(hst->virial), size, cudaMemcpyHostToDevice);
//   cudaMemcpy (dst->temperature, &(hst->temperature), size, cudaMemcpyHostToDevice);
//   cudaMemcpy (dst->pressure, &(hst->pressure), size, cudaMemcpyHostToDevice);
// }

// __global__ void clearDeviceStatistic (DeviceStatistic * dst)
// {
//   *(dst->bondedP) = 0;
//   *(dst->nonBondedP) = 0;
//   *(dst->electrostaticP) = 0;
//   *(dst->totalEnergy) = 0;
//   *(dst->kinetic) = 0;
//   *(dst->virial) = 0;
//   *(dst->temperature) = 0;
//   *(dst->pressure) = 0;
// }

