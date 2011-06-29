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
  mutable ScalorType *ddata;
public:
  MDStatistic ();
  MDStatistic (const MDSystem & sys);
  ~MDStatistic ();
  void reinit (const MDSystem & sys);
public:
  void clearDevice ();
  void updateHost () const;
  // ScalorType getStatistic (mdStatisticItem_t item) {return hdata[item];}
public:
  void deviceCopy (const MDStatistic & st);
  void deviceAdd  (const MDStatistic & st);
  void copy (const MDStatistic & st,
	     const IndexType num,
	     const mdStatisticItem_t items[NumberOfStatisticItems]);
  void add  (const MDStatistic & st,
	     const IndexType num,
	     const mdStatisticItem_t items[NumberOfStatisticItems]);
public:
  ScalorType kineticEnergy () const;
  ScalorType pressureXX (const RectangularBox & box) const;
  ScalorType pressureYY (const RectangularBox & box) const;
  ScalorType pressureZZ (const RectangularBox & box) const;
  ScalorType pressure   (const RectangularBox & box) const;
  ScalorType virial () const;
  ScalorType virialXX () const;
  ScalorType virialYY () const;
  ScalorType virialZZ () const;
  ScalorType nonBondedEnergy () const;
  ScalorType bondedEnergy () const;
public:
  void setEnergyCorr (const ScalorType & energyCorr);
  void setPressureCorr (const ScalorType & pressureCorr);
};

inline ScalorType MDStatistic::nonBondedEnergy() const
{
  return hdata[mdStatisticNonBondedPotential] +
      hdata[mdStatisticEnergyCorrection] +
      hdata[mdStatisticTwinRangeEnergyCorrection];
}

inline ScalorType MDStatistic::bondedEnergy () const
{
  return hdata[mdStatisticBondedPotential];
}

inline ScalorType MDStatistic::kineticEnergy () const
{
  return hdata[mdStatisticKineticEnergyXX] +
      hdata[mdStatisticKineticEnergyYY] +
      hdata[mdStatisticKineticEnergyZZ];
}

inline ScalorType MDStatistic::pressureXX (const RectangularBox & box) const
{
  return 2. * box.sizei.x * box.sizei.y * box.sizei.z * 
      (hdata[mdStatisticKineticEnergyXX] -
       hdata[mdStatisticVirialXX] * 0.5) +
      hdata[mdStatisticPressureCorrection] +
      hdata[mdStatisticTwinRangePressureCorrection];
}

inline ScalorType MDStatistic::pressureYY (const RectangularBox & box) const
{
  return 2. * box.sizei.x * box.sizei.y * box.sizei.z *
      (hdata[mdStatisticKineticEnergyYY] -
       hdata[mdStatisticVirialYY] * 0.5) +
      hdata[mdStatisticPressureCorrection] +
      hdata[mdStatisticTwinRangePressureCorrection];
}

inline ScalorType MDStatistic::pressureZZ (const RectangularBox & box) const
{
  return 2. * box.sizei.x * box.sizei.y * box.sizei.z *
      (hdata[mdStatisticKineticEnergyZZ] -
       hdata[mdStatisticVirialZZ] * 0.5) +
      hdata[mdStatisticPressureCorrection] +
      hdata[mdStatisticTwinRangePressureCorrection];
}

inline ScalorType MDStatistic::pressure (const RectangularBox & box) const
{
  return (pressureXX(box) + pressureYY(box) + pressureZZ(box)) / 3.;
}

inline ScalorType MDStatistic::virial() const
{
  return (virialXX () + virialYY() + virialZZ()) / 3.;
}

inline ScalorType MDStatistic::virialXX () const
{
  return hdata[mdStatisticVirialXX];
}

inline ScalorType MDStatistic::virialYY () const
{
  return hdata[mdStatisticVirialYY];
}

inline ScalorType MDStatistic::virialZZ () const
{
  return hdata[mdStatisticVirialZZ];
}





#endif
