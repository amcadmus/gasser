#ifndef __Parallel_Statistic_h_wanghan__
#define __Parallel_Statistic_h_wanghan__

#include "common.h"
#include "Parallel_CellList.h"

#define NumberOfStatisticItems	32

namespace Parallel{
  enum mdStatisticItem {
    mdStatistic_BondedPotential			= 0,
    mdStatistic_NonBondedPotential		= 1,
    mdStatistic_ElectrostaticPotential		= 2,
    mdStatistic_KineticEnergyXX			= 10,
    mdStatistic_KineticEnergyYY			= 11,
    mdStatistic_KineticEnergyZZ			= 12,
    mdStatistic_VirialXX				= 4,
    mdStatistic_VirialYY				= 5,
    mdStatistic_VirialZZ				= 6,
    mdStatistic_VirialXY				= 7,
    mdStatistic_VirialXZ				= 8,
    mdStatistic_VirialYZ				= 9
  };
  typedef enum mdStatisticItem mdStatisticItem_t;

  class HostStatistic
  {
    ScalorType * localData;
    ScalorType * globalData;
    size_t size;
    void clear ();
public:
    HostStatistic ();
    ~HostStatistic ();
    HostStatistic (const HostCellListedMDData & sys);
    // void reinit (const HostCellListedMDData & sys);
public:
    ScalorType * cptr_localStatisticData  () {return localData;}
    ScalorType * cptr_globalStatisticData () {return globalData;}
    const ScalorType * cptr_localStatisticData  () const {return localData;}
    const ScalorType * cptr_globalStatisticData () const {return globalData;}
    void clearData	();
    void collectData    ();
    void collectDataAll ();
    void collectData    (const mdStatisticItem_t item);
    void collectDataAll (const mdStatisticItem_t item);
public:
    ScalorType kineticEnergy ();
    ScalorType pressureXX (const RectangularBox & box);
    ScalorType pressureYY (const RectangularBox & box);
    ScalorType pressureZZ (const RectangularBox & box);
    ScalorType pressure   (const RectangularBox & box);
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
    size_t size;
    void clear ();
private:
    ScalorType *ddata;
public:
    DeviceStatistic  ();
    ~DeviceStatistic ();
    // void reinit (const DeviceCellListedMDData & ddata);
public:
    ScalorType * dptr_statisticData () {return ddata;}
    const ScalorType * dptr_statisticData () const {return ddata;}
    void clearData ();
    void copy (const DeviceStatistic & st);
    void copyToHost (HostStatistic & hst);
    void add  (const DeviceStatistic & st);
    void rescale (const IndexType  & position,
		  const IndexType  & num,
		  const ScalorType & scale);
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
  return globalData[mdStatistic_NonBondedPotential];
}

inline ScalorType Parallel::HostStatistic::BondedEnergy ()
{
  return globalData[mdStatistic_BondedPotential];
}

inline ScalorType Parallel::HostStatistic::kineticEnergy ()
{
  return globalData[mdStatistic_KineticEnergyXX] +
      globalData[mdStatistic_KineticEnergyYY] +
      globalData[mdStatistic_KineticEnergyZZ];
}

inline ScalorType Parallel::HostStatistic::
pressureXX (const RectangularBox & box)
{
  double volumei = box.sizei.x * box.sizei.y * box.sizei.z;
  return 2. * volumei * (globalData[mdStatistic_KineticEnergyXX] -
			 globalData[mdStatistic_VirialXX] * 0.5);
}

inline ScalorType Parallel::HostStatistic::
pressureYY (const RectangularBox & box)
{
  double volumei = box.sizei.x * box.sizei.y * box.sizei.z;
  return 2. * volumei * (globalData[mdStatistic_KineticEnergyYY] -
			 globalData[mdStatistic_VirialYY] * 0.5);
}

inline ScalorType Parallel::HostStatistic::
pressureZZ (const RectangularBox & box)
{
  double volumei = box.sizei.x * box.sizei.y * box.sizei.z;
  return 2. * volumei * (globalData[mdStatistic_KineticEnergyZZ] -
			 globalData[mdStatistic_VirialZZ] * 0.5);
}

inline ScalorType Parallel::HostStatistic::
pressure (const RectangularBox & box)
{
  return (pressureXX(box) +
	  pressureYY(box) +
	  pressureZZ(box)) * .33333333333333333333;
}

inline ScalorType Parallel::HostStatistic::virial()
{
  return virialXX () + virialYY() + virialZZ();
}

inline ScalorType Parallel::HostStatistic::virialXX ()
{
  return globalData[mdStatistic_VirialXX];
}

inline ScalorType Parallel::HostStatistic::virialYY ()
{
  return globalData[mdStatistic_VirialYY];
}

inline ScalorType Parallel::HostStatistic::virialZZ ()
{
  return globalData[mdStatistic_VirialZZ];
}


#endif
