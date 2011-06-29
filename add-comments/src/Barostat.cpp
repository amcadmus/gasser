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





Barostat_ParrinelloRahman::
Barostat_ParrinelloRahman()
{
  clearGroups();
}


void Barostat_ParrinelloRahman::
reinit (const ScalorType & dt_,
	const ScalorType & tau_P,
	const RectangularBox & box)
{
  dt  = dt_;
  tau = tau_P;
  scalor = 4 * M_PI * M_PI / (3. * tau * tau);
  b[0] = box.size.x;
  b[1] = box.size.y;
  b[2] = box.size.z;
  c[2] = c[1] = c[0] = 0;  
}

void Barostat_ParrinelloRahman::
assignGroup (const BoxDirection_t & dirs,
	     const ScalorType & refP_,
	     const ScalorType & beta_)
{
  if (numPcoupleGroup == 3) return;
  
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

void Barostat_ParrinelloRahman::
clearGroups ()
{
  numPcoupleGroup = 0;
  numDirInGroup[0] = numDirInGroup[1] = numDirInGroup[2] = 0;
  refP[0] = refP[1] = refP[2] = 0.;
  beta[0] = beta[1] = beta[2] = 0.;
}

void Barostat_ParrinelloRahman::
calCouple (const ScalorType * nowP,
	   ScalorType * lambda,
	   RectangularBox & box) const
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

  ScalorType cold[3];
  cold[0] = c[0];
  cold[1] = c[1];
  cold[2] = c[2];
  ScalorType L;
  L = (b[0] > b[1]) ? b[0] : b[1];
  L = (L > b[2]) ? L : b[2];
  
  ScalorType V = b[0] * b[1] * b[2];
  c[0] += dt * V * scalor * beta[0] / L * (tmpP[0] - refP[0]) / b[0];
  c[1] += dt * V * scalor * beta[1] / L * (tmpP[1] - refP[1]) / b[1];
  c[2] += dt * V * scalor * beta[2] / L * (tmpP[2] - refP[2]) / b[2];

  // if (dt * c[0] > 0.02) {
  //   c[0] = 0.02;
  // }
  // else if (dt * c[0] < -0.02){
  //   c[0] = -0.02;
  // }
  // if (dt * c[1] > 0.02) {
  //   c[1] = 0.02;
  // }
  // else if (dt * c[1] < -0.02){
  //   c[1] = -0.02;
  // }
  // if (dt * c[2] > 0.02) {
  //   c[2] = 0.02;
  // }
  // else if (dt * c[2] < -0.02){
  //   c[2] = -0.02;
  // }
  

  lambda[0] = (c[0] + cold[0]) / b[0];
  lambda[1] = (c[1] + cold[1]) / b[1];
  lambda[2] = (c[2] + cold[2]) / b[2];
  
  b[0] += dt * c[0];
  b[1] += dt * c[1];
  b[2] += dt * c[2];

  setBoxSize (b[0], b[1], b[2], box);
}



