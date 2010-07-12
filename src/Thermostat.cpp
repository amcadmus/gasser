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

