#ifndef __Parallel_Statistic_h_wanghan__
#define __Parallel_Statistic_h_wanghan__

#include "common.h"
#include "Parallel_CellList.h"

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

  class HostStatistic
  {
    ScalorType * localData;
    ScalorType * globalData;
    ScalorType volume;
    ScalorType volumei;
    size_t size;
public:
    HostStatistic ();
    ~HostStatistic ();
    HostStatistic (const HostCellListedMDData & sys);
    void reinit (const HostCellListedMDData & sys);
public:
    ScalorType * cptr_localStatisticData  () {return localData;}
    ScalorType * cptr_globalStatisticData () {return globalData;}
    const ScalorType * cptr_localStatisticData  () const {return localData;}
    const ScalorType * cptr_globalStatisticData () const {return globalData;}
    void collectData    ();
    void collectDataAll ();
    void collectData    (const mdStatisticItem_t item);
    void collectDataAll (const mdStatisticItem_t item);
public:
    ScalorType kineticEnergy ();
    ScalorType pressureXX ();
    ScalorType pressureYY ();
    ScalorType pressureZZ ();
    ScalorType pressure   ();
    ScalorType virial ();
    ScalorType virialXX ();
    ScalorType virialYY ();
    ScalorType virialZZ ();
    ScalorType NonBondedEnergy ();
    ScalorType BondedEnergy ();
  };
}

#ifdef DEVICE_CODE

#include "Parallel_CellList.h"

namespace Parallel{ 
  class DeviceStatistic 
  {
private:
    bool dmalloced;
    ScalorType volume;
    size_t size;
    // void clear ();
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
    void copyToHost (HostStatistic & hst);
    void add  (const DeviceStatistic & st);
  };

  namespace CudaGlobal{
    __global__ void 
    clearStatisticData (ScalorType *ddata);
    __global__ void 
    addStatisticData (ScalorType * ddata,
		      const ScalorType * cddata);
  };  
}
#endif


inline ScalorType Parallel::HostStatistic::NonBondedEnergy()
{
  return globalData[mdStatisticNonBondedPotential];
}

inline ScalorType Parallel::HostStatistic::BondedEnergy ()
{
  return globalData[mdStatisticBondedPotential];
}

inline ScalorType Parallel::HostStatistic::kineticEnergy ()
{
  return globalData[mdStatisticKineticEnergyXX] +
      globalData[mdStatisticKineticEnergyYY] +
      globalData[mdStatisticKineticEnergyZZ];
}

inline ScalorType Parallel::HostStatistic::pressureXX ()
{
  return 2. * volumei * (globalData[mdStatisticKineticEnergyXX] -
			 globalData[mdStatisticVirialXX] * 0.5);
}

inline ScalorType Parallel::HostStatistic::pressureYY ()
{
  return 2. * volumei * (globalData[mdStatisticKineticEnergyYY] -
			 globalData[mdStatisticVirialYY] * 0.5);
}

inline ScalorType Parallel::HostStatistic::pressureZZ ()
{
  return 2. * volumei * (globalData[mdStatisticKineticEnergyZZ] -
			 globalData[mdStatisticVirialZZ] * 0.5);
}

inline ScalorType Parallel::HostStatistic::pressure ()
{
  return (pressureXX() + pressureYY() + pressureZZ()) * .33333333333333333333;
}

inline ScalorType Parallel::HostStatistic::virial()
{
  return virialXX () + virialYY() + virialZZ();
}

inline ScalorType Parallel::HostStatistic::virialXX ()
{
  return globalData[mdStatisticVirialXX];
}

inline ScalorType Parallel::HostStatistic::virialYY ()
{
  return globalData[mdStatisticVirialYY];
}

inline ScalorType Parallel::HostStatistic::virialZZ ()
{
  return globalData[mdStatisticVirialZZ];
}


#endif
