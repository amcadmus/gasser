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
  virtual ScalorType calCouple (const ScalorType & Kientic_E) const = 0;
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
  mutable ScalorType xi;
  ScalorType scalor;
public:
  void reinit (const ScalorType & refT,
	       const ScalorType & dt,
	       const ScalorType & tau_T,
	       const IndexType & NKFreedom);
  ScalorType calCouple (const ScalorType & kientic_E) const;
};


#ifdef DEVICE_CODE

#include "SumVector.h"
#include "MDSystem_interface.h"
#include "Statistic_interface.h"
#include "MDTimer_interface.h"

class NoseHoover_Chains2
{
  ScalorType * dK;
  ScalorType hK;
private:
  dim3 atomGridDim;
  dim3 myBlockDim;
  IndexType sharedBuffSize;
public:
  ScalorType refT;
  ScalorType LL;
  ScalorType xi1, xi2;
  ScalorType vxi1, vxi2;
  ScalorType Q1, Q2;
  SumVector<ScalorType> sum_kxx;
  SumVector<ScalorType> sum_kyy;
  SumVector<ScalorType> sum_kzz;
  SumVector<ScalorType> sum_k;
public:
  NoseHoover_Chains2 ();
  ~NoseHoover_Chains2 ();
public:
  void reinit (const MDSystem &sys,
	       const IndexType & NThread,
	       const ScalorType & ref_T,
	       const ScalorType & tau_T);
  void operator_L_xi (const ScalorType & dt);
  void operator_L_Cv (const ScalorType & dt,
		      MDSystem & sys,
		      MDStatistic & st);
  void operator_L_Cv (const ScalorType & dt,
		      MDSystem & sys);
  void operator_L_G1 (const ScalorType & dt,
		      const MDSystem & sys);
  void operator_L_G2 (const ScalorType & dt);
  void operator_L_vxi1 (const ScalorType & dt);
  ScalorType HamiltonianContribution () const;

  void operator_L (const ScalorType & dt,
		   MDSystem & sys,
		   MDTimer * timer = NULL);
  void operator_L (const ScalorType & dt,
		   MDSystem & sys,
		   MDStatistic & st,
		   MDTimer * timer = NULL);
}
    ;
#endif


#endif
