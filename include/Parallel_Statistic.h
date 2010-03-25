#ifndef __Parallel_Statistic_h_wanghan__
#define __Parallel_Statistic_h_wanghan__

#include "common.h"

#define NumberOfStatisticItems	32

namespace Parallel{
  enum mdStatisticItem {
    mdStatisticBondedPotential			= 0,
    mdStatisticNonBondedPotential		= 1,
    mdStatisticElectrostaticPotential		= 2,
    mdStatisticKineticEnergyXX			= 10,
    mdStatisticKineticEnergyYY			= 11,
    mdStatisticKineticEnergyZZ			= 12,
    mdStatisticVirialXX				= 4,
    mdStatisticVirialYY				= 5,
    mdStatisticVirialZZ				= 6,
    mdStatisticVirialXY				= 7,
    mdStatisticVirialXZ				= 8,
    mdStatisticVirialYZ				= 9
  };
  typedef enum mdStatisticItem mdStatisticItem_t;
}
  
#ifdef DEVICE_CODE

#include "Parallel_CellList.h"

namespace Parallel{ 
  class DeviceStatistic 
  {
private:
    bool dmalloced;
    ScalorType volume;
    void clear ();
public:
    ScalorType *ddata;
public:
    DeviceStatistic  () : dmalloced(false) {}
    DeviceStatistic  (const DeviceCellListedMDData & sys)
	: dmalloced(false) {reinit(sys);}
    ~DeviceStatistic ();
    void reinit (const DeviceCellListedMDData & sys);
public:
    ScalorType * dptr_statisticData () {return ddata;}
    const ScalorType * dptr_statisticData () const {return ddata;}
    void clearDevice ();
    void copy (const DeviceStatistic & st);
    void add  (const DeviceStatistic & st);
  };

  namespace CudaGlobal{
    __global__ void 
    clearStatisticData (ScalorType *ddata);
    __global__ void 
    addStatisticData (ScalorType * ddata,
		      const ScalorType * cddata);
  };

#endif
  
}

#endif
