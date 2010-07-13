#define MPI_CODE

#include "Barostat.h"
#include "RandomGenerator.h"
#include <cmath>

#include "compile_error_mixcode.h"


void Barostat_Berendsen::
reinit (const ScalorType & dt_,
	const ScalorType & tau_P)
{
  dt  = dt_;
  tau = tau_P;
  scalor = dt / (tau * 3.);

  numPcoupleGroup = 0;
  numDirInGroup[0] = numDirInGroup[1] = numDirInGroup[2] = 0;
  refP[0] = refP[1] = refP[2] = 0.;
  beta[0] = beta[1] = beta[2] = 0.;
}

void Barostat_Berendsen::
assignGroup (const BoxDirection_t & dirs,
	     const ScalorType & refP_,
	     const ScalorType & beta_)
{
  PcoupleDirections[numPcoupleGroup] = dirs;
  if (dirs & mdRectBoxDirectionX){
    dirIndex[numPcoupleGroup][numDirInGroup[numPcoupleGroup]] = 0;
    numDirInGroup[numPcoupleGroup] ++;
  }
  if (dirs & mdRectBoxDirectionY){
    dirIndex[numPcoupleGroup][numDirInGroup[numPcoupleGroup]] = 1;
    numDirInGroup[numPcoupleGroup] ++;
  }
  if (dirs & mdRectBoxDirectionZ){
    dirIndex[numPcoupleGroup][numDirInGroup[numPcoupleGroup]] = 2;
    numDirInGroup[numPcoupleGroup] ++;
  }

  for (IndexType i = 0; i < numDirInGroup[numPcoupleGroup]; ++i){
    refP[dirIndex[numPcoupleGroup][i]] = refP_;
    beta[dirIndex[numPcoupleGroup][i]] = beta_;
  }
  numPcoupleGroup ++;
}

void Barostat_Berendsen::
clearGroups ()
{
  numPcoupleGroup = 0;
  numDirInGroup[0] = numDirInGroup[1] = numDirInGroup[2] = 0;
  refP[0] = refP[1] = refP[2] = 0.;
  beta[0] = beta[1] = beta[2] = 0.;
}

void Barostat_Berendsen::
calScale (const ScalorType * nowP,
	  ScalorType * scale) const
{
  ScalorType tmpP[3];
  tmpP[0] = nowP[0];
  tmpP[1] = nowP[1];
  tmpP[2] = nowP[2];
  
  for (IndexType i = 0; i < numPcoupleGroup; ++i) {
    if (numDirInGroup[i] > 1){
      ScalorType tmp = 0;
      for (IndexType j = 0; j < numDirInGroup[i]; ++j) {
	tmp += tmpP[dirIndex[i][j]];
      }
      tmp /= numDirInGroup[i];
      for (IndexType j = 0; j < numDirInGroup[i]; ++j){
	tmpP[dirIndex[i][j]] = tmp;
      }
    }
  }
  
  scale[0] = 1 + scalor * beta[0] * (tmpP[0] - refP[0]);
  scale[1] = 1 + scalor * beta[1] * (tmpP[1] - refP[1]);
  scale[2] = 1 + scalor * beta[2] * (tmpP[2] - refP[2]);
}


