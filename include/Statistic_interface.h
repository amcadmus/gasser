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

/// The MD statistic, calculating physical properties of interest.
/**
 * Now can calculate a series of thermodynamic quantities: non-bonded
 * energy, bonded energy, kinetic energy, pressure (\f$ p_{xx}\f$, \f$
 * p_{yy}\f$ and \f$ p_{zz}\f$) and virial (\f$ v_{xx}\f$, \f$
 * v_{yy}\f$ and \f$ v_{zz}\f$).
 */

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
  /** 
   * Calculate the kinetic energy.
   * 
   * @return The kinetic energy.
   */
  ScalorType kineticEnergy () const;
  /** 
   * Calculate the xx component of the pressure tensor.
   * 
   * @return The xx component of the pressure tensor.
   */
  ScalorType pressureXX (const RectangularBox & box) const;
  /** 
   * Calculate the yy component of the pressure tensor.
   * 
   * @return The yy component of the pressure tensor.
   */
  ScalorType pressureYY (const RectangularBox & box) const;
  /** 
   * Calculate the zz component of the pressure tensor.
   * 
   * @return The zz component of the pressure tensor.
   */
  ScalorType pressureZZ (const RectangularBox & box) const;
  /** 
   * Calculate the pressure.
   * 
   * @return The pressure.
   */
  ScalorType pressure   (const RectangularBox & box) const;
  /** 
   * Calculate the virial.
   * 
   * @return The virial.
   */
  ScalorType virial () const;
  /** 
   * Calculate the xx component of the virial tensor.
   * 
   * @return The xx component of the virial tensor.
   */
  ScalorType virialXX () const;
  /** 
   * Calculate the yy component of the virial tensor.
   * 
   * @return The yy component of the virial tensor.
   */
  ScalorType virialYY () const;
  /** 
   * Calculate the zz component of the virial tensor.
   * 
   * @return The zz component of the virial tensor.
   */
  ScalorType virialZZ () const;
  /** 
   * Calculate the non-bonded energy.
   * 
   * @return The non-bonded energy.
   */
  ScalorType nonBondedEnergy () const;
  /** 
   * Calculate the bonded energy.
   * 
   * @return The bonded energy.
   */
  ScalorType bondedEnergy () const;
public:
  /** 
   * Set the energy correction according to the cutoff of short-range
   * interactions.
   * 
   * @param energyCorr The energy correction.
   */
  void setEnergyCorr (const ScalorType & energyCorr);
  /** 
   * Set the pressure correction according to the cutoff of short-range
   * interactions.
   * 
   * @param pressureCorr The pressure correction.
   */
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
