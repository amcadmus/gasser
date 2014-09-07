#ifndef __Thermostat_h_wanghan__
#define __Thermostat_h_wanghan__

#include "common.h"

class Thermostat 
{
};

class Thermostat_VRescale : public Thermostat
{
public:
  virtual ScalorType calScale (const ScalorType & Kientic_E) const = 0;
};

class Thermostat_VCouple : public Thermostat
{
public:
  // virtual ScalorType calCouple (const ScalorType & Kientic_E) const = 0;
  virtual ScalorType getCouple () const = 0;
};

class Thermostat_VelocityRescale : public Thermostat_VRescale
{
  ScalorType tau;
  ScalorType dt;
  IndexType  Nf;
  ScalorType refK;
  ScalorType scalor1, scalor2;
public:
  void reinit (const ScalorType & refT,
	       const ScalorType & dt,
	       const ScalorType & tau_T,
	       const IndexType & NKFreedom);
  ScalorType calScale (const ScalorType & kientic_E) const;
};

class Thermostat_Berendsen : public Thermostat_VRescale
{
  ScalorType tau;
  ScalorType dt;
  IndexType  Nf;
  ScalorType refK;
  ScalorType scalor;
public:
  void reinit (const ScalorType & refT,
	       const ScalorType & dt,
	       const ScalorType & tau_T,
	       const IndexType & NKFreedom);
  ScalorType calScale (const ScalorType & kientic_E) const;
};

class Thermostat_NoseHoover : public Thermostat_VCouple
{
  ScalorType dt;
  ScalorType tau;
  ScalorType refT;
  IndexType  Nf;
  ScalorType refK;
  ScalorType xi;
  ScalorType scalor;
public:
  void reinit (const ScalorType & refT,
	       const ScalorType & dt,
	       const ScalorType & tau_T,
	       const IndexType & NKFreedom);
  void integrate_LeapFrog (const ScalorType & dt,
			   const ScalorType kineticE);
  void integrate_VelocityVerlet (const ScalorType & dt,
				 const ScalorType kineticE);
  ScalorType calCouple () const {return xi;}
};


#endif
