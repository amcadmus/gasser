#ifndef __Statistic_h_wanghan__
#define __Statistic_h_wanghan__

#include "common.h"

#define Size_StatisticData 10



struct StatisticData
{
  ScalorType * data;
};

// needs 1 thread
__global__ void clearStatisticData	(StatisticData sdata);

static __device__ ScalorType * bondedP		(StatisticData sdata);
static __device__ ScalorType * nonBondedP	(StatisticData sdata);
static __device__ ScalorType * electrostaticP	(StatisticData sdata);
static __device__ ScalorType * kineticE	(StatisticData sdata);
static __device__ ScalorType * virialxx	(StatisticData sdata);
static __device__ ScalorType * virialyy	(StatisticData sdata);
static __device__ ScalorType * virialzz	(StatisticData sdata);
static __device__ ScalorType * virialxy	(StatisticData sdata);
static __device__ ScalorType * virialxz	(StatisticData sdata);
static __device__ ScalorType * virialyz	(StatisticData sdata);

inline ScalorType * ptr_bondedP	(StatisticData &sdata);
inline ScalorType * ptr_nonBondedP	(StatisticData &sdata);
inline ScalorType * ptr_electrostaticP (StatisticData &sdata);
inline ScalorType * ptr_kineticE	(StatisticData &sdata);
inline ScalorType * ptr_virialxx	(StatisticData &sdata);
inline ScalorType * ptr_virialyy	(StatisticData &sdata);
inline ScalorType * ptr_virialzz	(StatisticData &sdata);
inline ScalorType * ptr_virialxy	(StatisticData &sdata);
inline ScalorType * ptr_virialxz	(StatisticData &sdata);
inline ScalorType * ptr_virialyz	(StatisticData &sdata);


__device__ ScalorType * bondedP		(StatisticData sdata){return &sdata.data[0];}
__device__ ScalorType * nonBondedP	(StatisticData sdata){return &sdata.data[1];}
__device__ ScalorType * electrostaticP	(StatisticData sdata){return &sdata.data[2];}
__device__ ScalorType * kineticE	(StatisticData sdata){return &sdata.data[3];}
__device__ ScalorType * virialxx	(StatisticData sdata){return &sdata.data[4];}
__device__ ScalorType * virialyy	(StatisticData sdata){return &sdata.data[5];}
__device__ ScalorType * virialzz	(StatisticData sdata){return &sdata.data[6];}
__device__ ScalorType * virialxy	(StatisticData sdata){return &sdata.data[7];}
__device__ ScalorType * virialxz	(StatisticData sdata){return &sdata.data[8];}
__device__ ScalorType * virialyz	(StatisticData sdata){return &sdata.data[9];}

ScalorType * ptr_bondedP	(StatisticData &sdata){return &sdata.data[0];}
ScalorType * ptr_nonBondedP	(StatisticData &sdata){return &sdata.data[1];}
ScalorType * ptr_electrostaticP(StatisticData &sdata){return &sdata.data[2];}
ScalorType * ptr_kineticE	(StatisticData &sdata){return &sdata.data[3];}
ScalorType * ptr_virialxx	(StatisticData &sdata){return &sdata.data[4];}
ScalorType * ptr_virialyy	(StatisticData &sdata){return &sdata.data[5];}
ScalorType * ptr_virialzz	(StatisticData &sdata){return &sdata.data[6];}
ScalorType * ptr_virialxy	(StatisticData &sdata){return &sdata.data[7];}
ScalorType * ptr_virialxz	(StatisticData &sdata){return &sdata.data[8];}
ScalorType * ptr_virialyz	(StatisticData &sdata){return &sdata.data[9];}

// struct DeviceStatistic
// {
//   ScalorType * bondedP;
//   ScalorType * nonBondedP;
//   ScalorType * electrostaticP;
//   ScalorType * totalEnergy;
//   ScalorType * kinetic;
//   ScalorType * virial;
//   ScalorType * temperature;
//   ScalorType * pressure;
// };

// struct HostStatistic
// {
//   ScalorType bondedP;
//   ScalorType nonBondedP;
//   ScalorType electrostaticP;
//   ScalorType totalEnergy;
//   ScalorType kinetic;
//   ScalorType virial;
//   ScalorType temperature;
//   ScalorType pressure;
// };

// __host__ void initHostStatistic (HostStatistic * hst);
// __host__ void initDeviceStatistic (DeviceStatistic * dst);
// __host__ void cpyDeviceStatisticToHost (const DeviceStatistic * dst,
// 					HostStatistic * hst);
// // needs 1 thread.
// __global__ void clearDeviceStatistic (DeviceStatistic * dst);





#endif
