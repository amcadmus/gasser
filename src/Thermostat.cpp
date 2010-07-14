#define MPI_CODE

#include "Thermostat.h"
#include "RandomGenerator.h"
#include <cmath>

#include "compile_error_mixcode.h"

void Thermostat_VelocityRescale::
reinit (const ScalorType & refT,
	const ScalorType & dt_,
	const ScalorType & tau_T,
	const IndexType & NKFreedom)
{
  dt  = dt_;
  tau = tau_T;
  Nf = NKFreedom;
  refK = Nf * 0.5 * refT;
  
  scalor1 = 1./tau;
  scalor2 = 1. / sqrt(Nf) / sqrt(tau) * 2; 
}
	
ScalorType Thermostat_VelocityRescale::
calScale (const ScalorType & kineticE) const
{
  double tmp;
  RandomGenerator_MT19937::genrand_Gaussian (0., sqrt(dt), &tmp);
  ScalorType tmp2 = ScalorType(tmp);
  tmp2 *= scalor2;

  ScalorType newK = kineticE + (refK - kineticE) * scalor1 * dt
      + sqrt(refK * kineticE) * tmp2;
  return sqrt(newK / kineticE);
}

void Thermostat_Berendsen::
reinit  (const ScalorType & refT,
	 const ScalorType & dt_,
	 const ScalorType & tau_T,
	 const IndexType & NKFreedom)
{
  tau = tau_T;
  scalor = 1./tau;
  dt = dt_;
  Nf = NKFreedom;
  refK = Nf * 0.5 * refT;
}

ScalorType Thermostat_Berendsen::
calScale (const ScalorType & kineticE) const
{
  return sqrtf (1.f + dt * scalor * (refK / kineticE - 1));
}

void Thermostat_NoseHoover::
reinit (const ScalorType & refT_,
	const ScalorType & dt_,
	const ScalorType & tau_T,
	const IndexType & NKFreedom)
{
  xi = 0;
  dt = dt_;
  tau = tau_T;
  refT = refT_;
  Nf = NKFreedom;
  refK = Nf * 0.5 * refT;
  
  scalor = 1./ (tau * tau * refT / (4. * M_PI * M_PI)) * 2. / Nf ;
}

ScalorType Thermostat_NoseHoover::
calCouple (const ScalorType & kineticE ) const
{
  // printf ("scalor is %f, xi is %f, refK is %f, nowK is %f \n",
  // 	  scalor, xi + dt * scalor * (kineticE - refK), refK, kineticE);
  return (xi += dt * scalor * (kineticE - refK));
}

