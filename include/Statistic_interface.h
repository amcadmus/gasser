#ifndef __Statistic_interface_h_wanghan__
#define __Statistic_interface_h_wanghan__

#include "common.h"
#include "Statistic.h"
#include "MDSystem_interface.h"

#define NumberOfStatisticItems	32

enum mdStatisticItem {
  mdStatisticBondedPotential			= 0,
  mdStatisticNonBondedPotential			= 1,
  mdStatisticElectrostaticPotential		= 2,
  mdStatisticKineticEnergyXX			= 10,
  mdStatisticKineticEnergyYY			= 11,
  mdStatisticKineticEnergyZZ			= 12,
  mdStatisticVirialXX				= 4,
  mdStatisticVirialYY				= 5,
  mdStatisticVirialZZ				= 6,
  mdStatisticVirialXY				= 7,
  mdStatisticVirialXZ				= 8,
  mdStatisticVirialYZ				= 9,
  mdStatisticEnergyCorrection			= 13,
  mdStatisticPressureCorrection			= 14,
  mdStatisticTwinRangeEnergyCorrection		= 15,
  mdStatisticTwinRangePressureCorrection	= 16
};
typedef enum mdStatisticItem mdStatisticItem_t;

class MDStatistic
{
  bool dmalloced;
  void clear ();
  // ScalorType volume;
public:
  ScalorType *hdata;
  ScalorType *ddata;
public:
  MDStatistic ();
  MDStatistic (const MDSystem & sys);
  ~MDStatistic ();
  void reinit (const MDSystem & sys);
public:
  void clearDevice ();
  void updateHost ();
//  ScalorType getStatistic (mdStatisticItem_t item) {return hdata[item];}
public:
  void deviceCopy (const MDStatistic & st);
  void deviceAdd  (const MDStatistic & st);
public:
  ScalorType nonBondedEnergy ();
  ScalorType bondedEnergy ();
  ScalorType kineticEnergy ();
  ScalorType pressureXX (const RectangularBox & box);
  ScalorType pressureYY (const RectangularBox & box);
  ScalorType pressureZZ (const RectangularBox & box);
  ScalorType pressure   (const RectangularBox & box);
  ScalorType virial ();
  ScalorType virialXX ();
  ScalorType virialYY ();
  ScalorType virialZZ ();
};



inline ScalorType MDStatistic::nonBondedEnergy()
{
  return hdata[mdStatisticNonBondedPotential] +
      hdata[mdStatisticEnergyCorrection] +
      hdata[mdStatisticTwinRangeEnergyCorrection];
}

inline ScalorType MDStatistic::bondedEnergy ()
{
  return hdata[mdStatisticBondedPotential];
}

inline ScalorType MDStatistic::kineticEnergy ()
{
  return hdata[mdStatisticKineticEnergyXX] +
      hdata[mdStatisticKineticEnergyYY] +
      hdata[mdStatisticKineticEnergyZZ];
}

inline ScalorType MDStatistic::pressureXX (const RectangularBox & box)
{
  return 2. * box.sizei.x * box.sizei.y * box.sizei.z * 
      (hdata[mdStatisticKineticEnergyXX] -
       hdata[mdStatisticVirialXX] * 0.5) +
      hdata[mdStatisticPressureCorrection] +
      hdata[mdStatisticTwinRangePressureCorrection];
}

inline ScalorType MDStatistic::pressureYY (const RectangularBox & box)
{
  return 2. * box.sizei.x * box.sizei.y * box.sizei.z *
      (hdata[mdStatisticKineticEnergyYY] -
       hdata[mdStatisticVirialYY] * 0.5) +
      hdata[mdStatisticPressureCorrection] +
      hdata[mdStatisticTwinRangePressureCorrection];
}

inline ScalorType MDStatistic::pressureZZ (const RectangularBox & box)
{
  return 2. * box.sizei.x * box.sizei.y * box.sizei.z *
      (hdata[mdStatisticKineticEnergyZZ] -
       hdata[mdStatisticVirialZZ] * 0.5) +
      hdata[mdStatisticPressureCorrection] +
      hdata[mdStatisticTwinRangePressureCorrection];
}

inline ScalorType MDStatistic::pressure (const RectangularBox & box)
{
  return (pressureXX(box) + pressureYY(box) + pressureZZ(box)) / 3.;
}

inline ScalorType MDStatistic::virial()
{
  return (virialXX () + virialYY() + virialZZ()) / 3.;
}

inline ScalorType MDStatistic::virialXX ()
{
  return hdata[mdStatisticVirialXX];
}

inline ScalorType MDStatistic::virialYY ()
{
  return hdata[mdStatisticVirialYY];
}

inline ScalorType MDStatistic::virialZZ ()
{
  return hdata[mdStatisticVirialZZ];
}


// class Statistic
// {
//   dim3 atomGridDim;
//   dim3 myBlockDim;
// public:
//   StatisticData deviceData;
//   StatisticData hostData;
//   ScalorType * statistic_buff;
// public:
//   ~Statistic();
//   void init (const MDSystem & sys, 
// 	     const IndexType & NThread);
//   void clearDevice ();
//   void updateHost ();
// public:
//   ScalorType get_bondedP	() {return *(ptr_bondedP(hostData));}
//   ScalorType get_nonBondedP	() {return *(ptr_nonBondedP(hostData));}
//   ScalorType get_electrostaticP () {return *(ptr_electrostaticP(hostData));}
//   ScalorType get_kineticEp	() {return *(ptr_kineticE(hostData));}
//   ScalorType get_virialxx	() {return *(ptr_virialxx(hostData));}
//   ScalorType get_virialyy	() {return *(ptr_virialyy(hostData));}
//   ScalorType get_virialzz	() {return *(ptr_virialzz(hostData));}
//   ScalorType get_virialxy	() {return *(ptr_virialxy(hostData));}
//   ScalorType get_virialxz	() {return *(ptr_virialxz(hostData));}
//   ScalorType get_virialyz	() {return *(ptr_virialyz(hostData));}
// };


// typedef ScalorType *	deviceStatisticPtr;
  


#endif
