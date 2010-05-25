#include "Parallel_TransferEngine.h"
#include "Parallel_Integrator.h"
#include "MDException.h"
#include <stdlib.h>

#include "compile_error_mixcode.h"

Parallel::HostSystemMomentum::
HostSystemMomentum ()
{
  p = (double *) malloc (3 * sizeof(double));
  if (p == NULL) {
    throw MDExcptFailedMallocOnHost ("HostSystemMomentum::HostSystemMomentum",
				     "p", 3 * sizeof(double));
  }
  sump = (double *) malloc (3 * sizeof(double));
  if (sump == NULL) {
    throw MDExcptFailedMallocOnHost ("HostSystemMomentum::HostSystemMomentum",
				     "sump", 3 * sizeof(double));
  }
  p_f = (ScalorType *) malloc (3 * sizeof(ScalorType));
  if (p_f == NULL) {
    throw MDExcptFailedMallocOnHost ("HostSystemMomentum::HostSystemMomentum",
				     "p_f", 3 * sizeof(ScalorType));
  }
  p[0] = p[1] = p[2] = 0.;
  p_f[0] = p_f[1] = p_f[2] = 0.;
}

Parallel::HostSystemMomentum::
~HostSystemMomentum()
{
  free (p);
  free (sump);
  free (p_f);
}

void Parallel::HostSystemMomentum::
sumAll ()
{
  p[0] = p_f[0];
  p[1] = p_f[1];
  p[2] = p_f[2];
  
  Parallel::SummationEngine sume;
  sume.sumAll (p, 3, sump);

  p_f[0] = sump[0];
  p_f[1] = sump[1];
  p_f[2] = sump[2];
}

void Parallel::HostSystemMass::
sumAll()
{
  mass = p_mass;
  Parallel::SummationEngine sume;
  sume.sumAll (&mass, 1, &sumMass);
  p_mass = sumMass;
}

