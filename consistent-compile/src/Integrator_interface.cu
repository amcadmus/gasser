#define DEVICE_CODE

#include "Integrator_interface.h"
#include "Integrator.h"
#include <math.h>
#include "RandomGenerator.h"

__global__ void fillSums0 (ScalorType* sums)
{
  sums[0] = sums[1] = sums[2] = 0;
}

TranslationalFreedomRemover::~TranslationalFreedomRemover ()
{
  cudaFree (sums);
}

void TranslationalFreedomRemover::init (const MDSystem &sys,
					const IndexType & NThread)
{
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NThread;
  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);

  sharedBuffSize = NThread * 2 * sizeof(ScalorType);
  sum_x.reinit (nob, NThreadForSum);
  sum_y.reinit (nob, NThreadForSum);
  sum_z.reinit (nob, NThreadForSum);
  cudaMalloc ((void**)&sums,  3 * sizeof(ScalorType));
  fillSums0 <<<1,1>>> (sums);
}

void TranslationalFreedomRemover::remove (MDSystem & sys,
					  MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeRemoveTransFreedom);
  prepareRemoveTranslationalFreedom
      <<<atomGridDim, myBlockDim, sharedBuffSize>>>(
	  sys.ddata.numAtom,
	  sys.ddata.mass,
	  sys.ddata.velox,
	  sys.ddata.veloy,
	  sys.ddata.veloz,
	  sum_x.buff,
	  sum_y.buff,
	  sum_z.buff);
  checkCUDAError("TranslationalFreedomRemover::remove, prepare");
  sum_x.sumBuff (sums, 0);
  sum_y.sumBuff (sums, 1);
  sum_z.sumBuff (sums, 2);
  removeFreedom 
      <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.numAtom,
	  sys.ddata.velox,
	  sys.ddata.veloy,
	  sys.ddata.veloz,
	  sys.ddata.totalMassi,
	  sums);
  checkCUDAError("TranslationalFreedomRemover::remove, remove");
  if (timer != NULL) timer->toc(mdTimeRemoveTransFreedom);
}
  

LeapFrog::~LeapFrog()
{
}

void LeapFrog::
reinit (const MDSystem &sys,
	const IndexType & NThread)
{
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NThread;
  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);  

  sum_kxx.reinit (nob, NThreadForSum);
  sum_kyy.reinit (nob, NThreadForSum);
  sum_kzz.reinit (nob, NThreadForSum);  

  sharedBuffSize = NThread * sizeof(ScalorType);
}

void LeapFrog::step (MDSystem & sys,
		     const ScalorType & dt,
		     MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  leapFrog1Step
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, sys.ddata.massi,
#ifndef COORD_IN_ONE_VEC
	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
#else
	  sys.ddata.coord,
#endif
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz,
	  dt);
  checkCUDAError ("LeapFrog::step");
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}


void LeapFrog::step (MDSystem & sys,
		     const ScalorType & dt,
		     MDStatistic & st,
		     MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  leapFrog1Step
      <<<atomGridDim, myBlockDim, sharedBuffSize>>> (
	  sys.ddata.numAtom, sys.ddata.mass, sys.ddata.massi,
	  sys.ddata.coord,
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz,
	  dt,
	  sum_kxx.buff,
	  sum_kyy.buff,
	  sum_kzz.buff);
  checkCUDAError ("LeapFrog::step (with statistic)");
  sum_kxx.sumBuffAdd (st.ddata, mdStatisticKineticEnergyXX, 0);
  sum_kyy.sumBuffAdd (st.ddata, mdStatisticKineticEnergyYY, 0);
  sum_kzz.sumBuffAdd (st.ddata, mdStatisticKineticEnergyZZ, 0);
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}


void LeapFrog::stepX (MDSystem & sys,
		      const ScalorType & dt,
		      MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  leapFrogStepX
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, sys.ddata.massi,
#ifndef COORD_IN_ONE_VEC
	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
#else
	  sys.ddata.coord,
#endif
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  dt);
  checkCUDAError ("LeapFrog::stepX");
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}

void LeapFrog::stepV (MDSystem & sys,
		      const ScalorType & dt,
		      MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  leapFrogStepV
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, sys.ddata.massi,
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz,
	  dt);
  checkCUDAError ("LeapFrog::stepV");
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}

void LeapFrog::stepV (MDSystem & sys,
		      const ScalorType & dt,
		      MDStatistic & st,
		      MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  leapFrogStepV
      <<<atomGridDim, myBlockDim, sharedBuffSize>>> (
	  sys.ddata.numAtom, sys.ddata.mass, sys.ddata.massi,
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz,
	  dt,
	  sum_kxx.buff,
	  sum_kyy.buff,
	  sum_kzz.buff);
  checkCUDAError ("LeapFrog::stepV (with statistic)");
  sum_kxx.sumBuffAdd (st.ddata, mdStatisticKineticEnergyXX, 0);
  sum_kyy.sumBuffAdd (st.ddata, mdStatisticKineticEnergyYY, 0);
  sum_kzz.sumBuffAdd (st.ddata, mdStatisticKineticEnergyZZ, 0);
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}

void LeapFrog::
stepV_VCouple (MDSystem & sys,
	       const ScalorType & dt,
	       const ScalorType * lambda,
	       MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  leapFrogStepV_VCouple
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, sys.ddata.massi,
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz,
	  lambda[0], lambda[1], lambda[2],
	  dt);
  checkCUDAError ("LeapFrog::stepV");
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}

void LeapFrog::
stepV_VCouple (MDSystem & sys,
	       const ScalorType & dt,
	       const ScalorType * lambda,
	       MDStatistic & st,
	       MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  leapFrogStepV_VCouple
      <<<atomGridDim, myBlockDim, sharedBuffSize>>> (
	  sys.ddata.numAtom, sys.ddata.mass, sys.ddata.massi,
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz,
	  lambda[0], lambda[1], lambda[2],
	  dt,
	  sum_kxx.buff,
	  sum_kyy.buff,
	  sum_kzz.buff);
  checkCUDAError ("LeapFrog::stepV (with statistic)");
  sum_kxx.sumBuffAdd (st.ddata, mdStatisticKineticEnergyXX, 0);
  sum_kyy.sumBuffAdd (st.ddata, mdStatisticKineticEnergyYY, 0);
  sum_kzz.sumBuffAdd (st.ddata, mdStatisticKineticEnergyZZ, 0);
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}



VelocityVerlet::~VelocityVerlet()
{
}

void VelocityVerlet::init (const MDSystem &sys,
			   const IndexType & NThread)
{
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NThread;
  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);  

  sum_kxx.reinit (nob, NThreadForSum);
  sum_kyy.reinit (nob, NThreadForSum);
  sum_kzz.reinit (nob, NThreadForSum);  

  sharedBuffSize = NThread * sizeof(ScalorType);
}



void VelocityVerlet::step1 (MDSystem & sys,
			    const ScalorType & dt,
			    MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  velocityVerlet_part1
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, sys.ddata.massi,
#ifndef COORD_IN_ONE_VEC
	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
#else
	  sys.ddata.coord,
#endif
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz,
	  dt);
  checkCUDAError ("VelocityVerlet::step1");
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}

void VelocityVerlet::step2 (MDSystem & sys,
			    const ScalorType & dt,
			    MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  velocityVerlet_part2
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, sys.ddata.massi,
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz,
	  dt);
  checkCUDAError ("VelocityVerlet::step2");
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}

void VelocityVerlet::step2 (MDSystem & sys,
			    const ScalorType & dt,
			    MDStatistic & st,
			    MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  velocityVerlet_part2
      <<<atomGridDim, myBlockDim, sharedBuffSize>>> (
	  sys.ddata.numAtom, sys.ddata.mass, sys.ddata.massi,
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz,
	  dt,
	  sum_kxx.buff,
	  sum_kyy.buff,
	  sum_kzz.buff);
  checkCUDAError ("VelocityVerlet::step2 (with statistic)");
  sum_kxx.sumBuffAdd (st.ddata, mdStatisticKineticEnergyXX, 0);
  sum_kyy.sumBuffAdd (st.ddata, mdStatisticKineticEnergyYY, 0);
  sum_kzz.sumBuffAdd (st.ddata, mdStatisticKineticEnergyZZ, 0);
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}


VelocityRescale::~VelocityRescale()
{
  cudaFree (buff);
  cudaFree (kineticE);
}


void VelocityRescale::init (const MDSystem & sys, 
			    const IndexType & NThread,
			    const ScalorType & refT,
			    const ScalorType & tau_)
{
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NThread;
  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);

  tau = tau_;
  Nf = sys.ddata.NFreedom - 3;
  scalor1 = 1./tau;
  scalor2 = 1. / sqrt(Nf) / sqrt(tau) * 2; 

  refK = (Nf) * 0.5 * refT;
  printf ("# refK is %f\n", refK);
  
  cudaMalloc ((void**)&kineticE, sizeof(ScalorType) * 2);
  cudaMalloc ((void**)&buff, sizeof(ScalorType) * nob);
  checkCUDAError ("VelocityVerlet::init allocation");

  sum_kxx.reinit (nob, NThreadForSum);
  sum_kyy.reinit (nob, NThreadForSum);
  sum_kzz.reinit (nob, NThreadForSum);
  sum_k.reinit (nob, NThreadForSum);
  
  sharedBuffSize = NThread * sizeof(ScalorType);
}

void VelocityRescale::step1 (MDSystem & sys,
			     const ScalorType & dt,
			     MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  velocityVerlet_part1
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, sys.ddata.massi,
#ifndef COORD_IN_ONE_VEC
	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
#else
	  sys.ddata.coord,
#endif
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz,
	  dt);
  double tmp;
  RandomGenerator_MT19937::genrand_Gaussian (0., sqrt(dt), &tmp);
  tmp2 = ScalorType(tmp);
  tmp2 *= scalor2;
  checkCUDAError ("VelocityRescale::step1");
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}

void VelocityRescale::step2 (MDSystem & sys,
			     const ScalorType & dt,
			     MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  velocityVerlet_part2a
      <<<atomGridDim, myBlockDim, sharedBuffSize>>> (
	  sys.ddata.numAtom, sys.ddata.mass, sys.ddata.massi,
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz,
	  dt,
	  sum_k.buff);
  checkCUDAError ("VelocityRescale:: Veloccity Verlet step2");
  sum_k.sumBuff (kineticE, 0, 0);
  cudaMemcpy (&hkineticE, kineticE, sizeof(ScalorType), cudaMemcpyDeviceToHost);
  newK = hkineticE + (refK - hkineticE) * scalor1 * dt + sqrt(refK * hkineticE) * tmp2;
  alpha = sqrt(newK / hkineticE);
  velocityRescale_rescale 
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom,
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  alpha);
  checkCUDAError ("VelocityRescale::step2 rescale");
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}

void VelocityRescale::step2 (MDSystem & sys,
			     const ScalorType & dt,
			     MDStatistic &st,
			     MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  velocityVerlet_part2a
      <<<atomGridDim, myBlockDim, sharedBuffSize>>> (
	  sys.ddata.numAtom, sys.ddata.mass, sys.ddata.massi,
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz,
	  dt,
	  sum_k.buff);
  checkCUDAError ("VelocityRescale:: Veloccity Verlet step2");
  sum_k.sumBuff (kineticE, 0, 0);
  cudaMemcpy (&hkineticE, kineticE, sizeof(ScalorType), cudaMemcpyDeviceToHost);
  newK = hkineticE + (refK - hkineticE) * scalor1 * dt + sqrt(refK * hkineticE) * tmp2;
  alpha = sqrt(newK / hkineticE);
  velocityRescale_rescale 
      <<<atomGridDim, myBlockDim, sharedBuffSize>>> (
	  sys.ddata.numAtom,
	  sys.ddata.mass,
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  alpha,
	  sum_kxx.buff,
	  sum_kyy.buff,
	  sum_kzz.buff);
  checkCUDAError ("VelocityRescale::step2 rescale");
  sum_kxx.sumBuffAdd (st.ddata, mdStatisticKineticEnergyXX, 0);
  sum_kyy.sumBuffAdd (st.ddata, mdStatisticKineticEnergyYY, 0);
  sum_kzz.sumBuffAdd (st.ddata, mdStatisticKineticEnergyZZ, 0);
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}


BerendsenLeapFrog::BerendsenLeapFrog ()
{
  TCoupleOn = false;
  PCoupleOn = false;
  NPCoupleGroup = 0;
  ptr_inter = NULL;
  ptr_nlist = NULL;
  ptr_bdInterList = NULL;
  nstep = 0;
}


void
BerendsenLeapFrog::init (const MDSystem &sys,
			 const IndexType & NThread,
			 const ScalorType & dt_,
			 InteractionEngine &inter,
			 NeighborList & nlist,
			 const ScalorType & rebt,
			 BondedInteractionList * ptr_bdInterList_)
{
  dt = dt_;
  TCoupleOn = false;
  PCoupleOn = false;
  NPCoupleGroup = 0;
  nstep = 0;

  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NThread;
  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);  
  // cudaMalloc((void**)&(statistic_buff), sizeof(ScalorType) * nob);
  // checkCUDAError ("VelocityVerlet::init allocation");

  lpfrog.reinit (sys, NThread);
  ptr_inter = & inter;
  ptr_nlist = & nlist;
  rebuildThreshold = rebt;
  ptr_bdInterList = ptr_bdInterList_;
}


void
BerendsenLeapFrog::TCouple (const ScalorType & refT_,
			    const ScalorType & tauT_)
{
  TCoupleOn = true;
  refT = refT_;
  tauT = tauT_;
}

void
BerendsenLeapFrog::addPcoupleGroup (const PCoupleDirection_t & direction,
				    const ScalorType & refP_,
				    const ScalorType & tauP_,
				    const ScalorType & betaP_)
{
  if (direction == 0) return;
  PCoupleOn = true;

  if (NPCoupleGroup == 3){
    fprintf (stderr, "# too many P couple groups, add nothing" );
    return ;
  }
  refP[NPCoupleGroup] = refP_;
  tauP[NPCoupleGroup] = tauP_;
  betaP[NPCoupleGroup] = betaP_;
  PCoupleDirections[NPCoupleGroup] = direction;
  
  NPCoupleGroup ++;
}

// void
// BerendsenLeapFrog::firstStep (MDSystem & sys, MDStatistic &st, MDTimer * timer)
// {
//   myst.clearDevice ();
//   ptr_inter->clearInteraction (sys);
//   if (ptr_nlist != NULL){
//     ptr_inter->applyNonBondedInteraction (sys, *ptr_nlist, timer);
//   }
//   if (ptr_bdInterList != NULL){
//     ptr_inter->applyBondedInteraction (sys, *ptr_bdInterList, timer);
//   }
//   lpfrog.step (sys, dt, myst, timer);
//   if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
//     // printf("# rebuild at step %d\n", nstep);
//     // fflush(stdout);
//     ptr_nlist->build(sys, timer);
//   }
//   ptr_inter->clearInteraction (sys);
//   if (ptr_nlist != NULL){
//     ptr_inter->applyNonBondedInteraction (sys, *ptr_nlist, myst, timer);
//   }
//   if (ptr_bdInterList != NULL){
//     ptr_inter->applyBondedInteraction (sys, *ptr_bdInterList, myst, timer);
//   }
//   st.deviceAdd (myst);
//   nstep ++;
// }

// void
// BerendsenLeapFrog::firstStep (MDSystem & sys,  MDTimer * timer)
// {
//   myst.clearDevice ();
//   ptr_inter->clearInteraction (sys);
//   if (ptr_nlist != NULL){
//     ptr_inter->applyNonBondedInteraction (sys, *ptr_nlist, timer);
//   }
//   if (ptr_bdInterList != NULL){
//     ptr_inter->applyBondedInteraction (sys, *ptr_bdInterList, timer);
//   }
//   lpfrog.step (sys, dt, myst, timer);
//   if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
//     // printf("# rebuild at step %d\n", nstep);
//     // fflush(stdout);
//     ptr_nlist->build(sys, timer);
//   }
//   ptr_inter->clearInteraction (sys);
//   if (ptr_nlist != NULL){
//     ptr_inter->applyNonBondedInteraction (sys, *ptr_nlist, myst, timer);
//   }
//   if (ptr_bdInterList != NULL){
//     ptr_inter->applyBondedInteraction (sys, *ptr_bdInterList, myst, timer);
//   }
//   nstep ++;
// }


// void
// BerendsenLeapFrog::oneStep (MDSystem & sys, MDTimer * timer)
// {
//   ScalorType nowT, lambda;
//   ScalorType nowP[3], mu[3];
//   IndexType nDir[3];
  
//   if (timer != NULL) timer->tic (mdTimeIntegrator);
//   if (nstep != 0) {
//     myst.updateHost();
//     if (TCoupleOn){
//       nowT = myst.kineticEnergy();
//       nowT *= 2.f / (sys.ddata.NFreedom - 3);
//       lambda = sqrtf(1.f + dt / tauT * (refT / nowT - 1.f));
//     }
//     if (PCoupleOn){
//       for (IndexType i = 0; i < NPCoupleGroup; ++i){
// 	nowP[i] = 0;
// 	nDir[i] = 0;
// 	if ((PCoupleDirections[i] & PCoupleX) != 0){
// 	  nowP[i] += myst.pressureXX(sys.box);
// 	  nDir[i] ++;
// 	}
// 	if ((PCoupleDirections[i] & PCoupleY) != 0){
// 	  nowP[i] += myst.pressureYY(sys.box);
// 	  nDir[i] ++;
// 	}
// 	if ((PCoupleDirections[i] & PCoupleZ) != 0){
// 	  nowP[i] += myst.pressureZZ(sys.box);
// 	  nDir[i] ++;
// 	}
// 	nowP[i] /= ScalorType(nDir[i]);
// 	mu [i] = powf (1.f + dt / tauP[i] * betaP[i] * (nowP[i] - refP[i]), 1.f/3.f);
//       }
//     }
  
//     myst.clearDevice();
//     lpfrog.stepV (sys, dt, myst);
//     if (TCoupleOn){
//       rescaleProperty <<<atomGridDim, myBlockDim>>>(
// 	  sys.ddata.velox, sys.ddata.numAtom,
// 	  lambda);
//       rescaleProperty <<<atomGridDim, myBlockDim>>>(
// 	  sys.ddata.veloy, sys.ddata.numAtom,
// 	  lambda);
//       rescaleProperty <<<atomGridDim, myBlockDim>>>(
// 	  sys.ddata.veloz, sys.ddata.numAtom,
// 	  lambda);
//       rescaleProperty <<<1, 3>>>(
// 	  myst.ddata, mdStatisticKineticEnergyXX, 3,
// 	  lambda * lambda);
//     }
//     lpfrog.stepX (sys, dt);
//     if (PCoupleOn){
//       ScalorType newBoxX(sys.box.size.x);
//       ScalorType newBoxY(sys.box.size.y);
//       ScalorType newBoxZ(sys.box.size.z);
//       CoordType coordScalor ;
//       coordScalor.x = 1.f;
//       coordScalor.y = 1.f;
//       coordScalor.z = 1.f;
//       for (IndexType i = 0; i < NPCoupleGroup; ++i){
// 	if ((PCoupleDirections[i] & PCoupleX) != 0){
// 	  coordScalor.x *= mu[i];
// 	  newBoxX *= mu[i];
// 	}
// 	if ((PCoupleDirections[i] & PCoupleY) != 0){
// 	  coordScalor.y *= mu[i];
// 	  newBoxY *= mu[i];
// 	}
// 	if ((PCoupleDirections[i] & PCoupleZ) != 0){
// 	  coordScalor.z *= mu[i];
// 	  newBoxZ *= mu[i];
// 	}
//       }
//       rescaleCoord <<<atomGridDim, myBlockDim>>> (
// 	  sys.ddata.coord, sys.ddata.numAtom,
// 	  coordScalor);
//       sys.setBoxSize (newBoxX, newBoxY, newBoxZ);
//     }
//     nstep ++;
//     if (timer != NULL) timer->toc (mdTimeIntegrator);
//     if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
//       // printf("# rebuild at step %d\n", nstep);
//       // fflush(stdout);
//       ptr_nlist->build(sys, timer);
//     }
//     ptr_inter->clearInteraction (sys);
//     if (ptr_nlist != NULL){
//       ptr_inter->applyNonBondedInteraction (sys, *ptr_nlist, myst, timer);
//     }
//     if (ptr_bdInterList != NULL){
//       ptr_inter->applyBondedInteraction (sys, *ptr_bdInterList, myst, timer);
//     }
//   }
//   else {
//     firstStep (sys, timer);
//   }
// }

// void
// BerendsenLeapFrog::oneStep (MDSystem & sys, MDStatistic &st, MDTimer * timer)
// {
//   ScalorType nowT, lambda;
//   ScalorType nowP[3], mu[3];
//   IndexType nDir[3];
  
//   if (timer != NULL) timer->tic (mdTimeIntegrator);
//   if (nstep != 0) {
//     myst.updateHost();
//     if (TCoupleOn){
//       nowT = myst.kineticEnergy();
//       nowT *= 2.f / (sys.ddata.NFreedom - 3);
//       lambda = sqrtf(1.f + dt / tauT * (refT / nowT - 1.f));
//     }
//     if (PCoupleOn){
//       for (IndexType i = 0; i < NPCoupleGroup; ++i){
// 	nowP[i] = 0;
// 	nDir[i] = 0;
// 	if ((PCoupleDirections[i] & PCoupleX) != 0){
// 	  nowP[i] += myst.pressureXX(sys.box);
// 	  nDir[i] ++;
// 	}
// 	if ((PCoupleDirections[i] & PCoupleY) != 0){
// 	  nowP[i] += myst.pressureYY(sys.box);
// 	  nDir[i] ++;
// 	}
// 	if ((PCoupleDirections[i] & PCoupleZ) != 0){
// 	  nowP[i] += myst.pressureZZ(sys.box);
// 	  nDir[i] ++;
// 	}
// 	nowP[i] /= ScalorType(nDir[i]);
// 	mu [i] = powf (1.f + dt / tauP[i] * betaP[i] * (nowP[i] - refP[i]), 1.f/3.f);
//       }
//     }
  
//     myst.clearDevice();
//     lpfrog.stepV (sys, dt, myst);
//     if (TCoupleOn){
//       rescaleProperty <<<atomGridDim, myBlockDim>>>(
// 	  sys.ddata.velox, sys.ddata.numAtom,
// 	  lambda);
//       rescaleProperty <<<atomGridDim, myBlockDim>>>(
// 	  sys.ddata.veloy, sys.ddata.numAtom,
// 	  lambda);
//       rescaleProperty <<<atomGridDim, myBlockDim>>>(
// 	  sys.ddata.veloz, sys.ddata.numAtom,
// 	  lambda);
//       rescaleProperty <<<1, 3>>>(
// 	  myst.ddata, mdStatisticKineticEnergyXX, 3,
// 	  lambda * lambda);
//     }
//     lpfrog.stepX (sys, dt);
//     if (PCoupleOn){
//       ScalorType newBoxX(sys.box.size.x);
//       ScalorType newBoxY(sys.box.size.y);
//       ScalorType newBoxZ(sys.box.size.z);
//       CoordType coordScalor ;
//       coordScalor.x = 1.f;
//       coordScalor.y = 1.f;
//       coordScalor.z = 1.f;
//       for (IndexType i = 0; i < NPCoupleGroup; ++i){
// 	if ((PCoupleDirections[i] & PCoupleX) != 0){
// 	  coordScalor.x *= mu[i];
// 	  newBoxX *= mu[i];
// 	}
// 	if ((PCoupleDirections[i] & PCoupleY) != 0){
// 	  coordScalor.y *= mu[i];
// 	  newBoxY *= mu[i];
// 	}
// 	if ((PCoupleDirections[i] & PCoupleZ) != 0){
// 	  coordScalor.z *= mu[i];
// 	  newBoxZ *= mu[i];
// 	}
//       }
//       rescaleCoord <<<atomGridDim, myBlockDim>>> (
// 	  sys.ddata.coord, sys.ddata.numAtom,
// 	  coordScalor);
//       sys.setBoxSize (newBoxX, newBoxY, newBoxZ);
//     }
//     nstep ++;
//     if (timer != NULL) timer->toc (mdTimeIntegrator);
//     if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
//       // printf("# rebuild at step %d\n", nstep);
//       // fflush(stdout);
//       ptr_nlist->build(sys, timer);
//     }
//     ptr_inter->clearInteraction (sys);
//     if (ptr_nlist != NULL){
//       ptr_inter->applyNonBondedInteraction (sys, *ptr_nlist, myst, timer);
//     }
//     if (ptr_bdInterList != NULL){
//       ptr_inter->applyBondedInteraction (sys, *ptr_bdInterList, myst, timer);
//     }
//     st.deviceAdd (myst);
//   }
//   else {
//     firstStep (sys, st, timer);
//   }
// }














void LeapFrog_TPCouple_Rescale::
nullPointers ()
{
  ptr_thermostat = NULL;
  ptr_barostat = NULL;
  // ptr_barostat[0] = NULL;
  // ptr_barostat[1] = NULL;
  // ptr_barostat[2] = NULL;
  // NPCoupleGroup = 0;
  // PCoupleOn = false;
  // NPCoupleGroup = 0;
  ptr_inter = NULL;
  ptr_nlist = NULL;
  ptr_bdInterList = NULL;
  nstep = 0;
}
  

LeapFrog_TPCouple_Rescale::
LeapFrog_TPCouple_Rescale ()
{
  nullPointers();
}

void LeapFrog_TPCouple_Rescale::
init (const MDSystem &sys,
      const IndexType & NThread,
      const ScalorType & dt_,
      InteractionEngine &inter,
      NeighborList & nlist,
      const ScalorType & rebt,
      BondedInteractionList * ptr_bdInterList_)
{
  nullPointers();
  dt = dt_;
  nstep = 0;

  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NThread;
  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);  

  lpfrog.reinit (sys, NThread);
  ptr_inter = & inter;
  ptr_nlist = & nlist;
  rebuildThreshold = rebt;
  ptr_bdInterList = ptr_bdInterList_;
}


// void LeapFrog_TPCouple_Rescale::
// firstStep (MDSystem & sys,
// 	   MDStatistic &st,
// 	   MDTimer * timer)
// {
//   myst.clearDevice ();
//   ptr_inter->clearInteraction (sys);
//   if (ptr_nlist != NULL){
//     ptr_inter->applyNonBondedInteraction (sys, *ptr_nlist, timer);
//   }
//   if (ptr_bdInterList != NULL){
//     ptr_inter->applyBondedInteraction (sys, *ptr_bdInterList, timer);
//   }
//   lpfrog.step (sys, dt, myst, timer);
//   if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
//     // printf("# rebuild at step %d\n", nstep);
//     // fflush(stdout);
//     ptr_nlist->build(sys, timer);
//   }
//   ptr_inter->clearInteraction (sys);
//   if (ptr_nlist != NULL){
//     ptr_inter->applyNonBondedInteraction (sys, *ptr_nlist, myst, timer);
//   }
//   if (ptr_bdInterList != NULL){
//     ptr_inter->applyBondedInteraction (sys, *ptr_bdInterList, myst, timer);
//   }
//   st.deviceAdd (myst);
//   nstep ++;
// }

// void LeapFrog_TPCouple_Rescale::
// firstStep (MDSystem & sys,
// 	   MDTimer * timer)
// {
//   myst.clearDevice ();
//   ptr_inter->clearInteraction (sys);
//   if (ptr_nlist != NULL){
//     ptr_inter->applyNonBondedInteraction (sys, *ptr_nlist, timer);
//   }
//   if (ptr_bdInterList != NULL){
//     ptr_inter->applyBondedInteraction (sys, *ptr_bdInterList, timer);
//   }
//   lpfrog.step (sys, dt, myst, timer);
//   if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
//     // printf("# rebuild at step %d\n", nstep);
//     // fflush(stdout);
//     ptr_nlist->build(sys, timer);
//   }
//   ptr_inter->clearInteraction (sys);
//   if (ptr_nlist != NULL){
//     ptr_inter->applyNonBondedInteraction (sys, *ptr_nlist, myst, timer);
//   }
//   if (ptr_bdInterList != NULL){
//     ptr_inter->applyBondedInteraction (sys, *ptr_bdInterList, myst, timer);
//   }
//   nstep ++;
// }


// void LeapFrog_TPCouple_Rescale::
// oneStep (MDSystem & sys,
// 	 MDTimer * timer)
// {
//   ScalorType nowK, lambda;
//   ScalorType nowP[3], mu[3];
  
//   if (timer != NULL) timer->tic (mdTimeIntegrator);
//   if (nstep != 0) {
//     myst.updateHost();
//     if (ptr_thermostat != NULL){
//       nowK = myst.kineticEnergy();
//       lambda = ptr_thermostat->calScale (nowK);
//     }
//     if (ptr_barostat != 0){
//       nowP[0] = myst.pressureXX (sys.box);
//       nowP[1] = myst.pressureYY (sys.box);
//       nowP[2] = myst.pressureZZ (sys.box);
//       ptr_barostat->calScale (nowP, mu);
//       // for (IndexType i = 0; i < NPCoupleGroup; ++i){
//       // 	nowP[i] = 0;
//       // 	nDir[i] = 0;
//       // 	if ((PCoupleDirections[i] & PCoupleX) != 0){
//       // 	  nowP[i] += myst.pressureXX(sys.box);
//       // 	  nDir[i] ++;
//       // 	}
//       // 	if ((PCoupleDirections[i] & PCoupleY) != 0){
//       // 	  nowP[i] += myst.pressureYY(sys.box);
//       // 	  nDir[i] ++;
//       // 	}
//       // 	if ((PCoupleDirections[i] & PCoupleZ) != 0){
//       // 	  nowP[i] += myst.pressureZZ(sys.box);
//       // 	  nDir[i] ++;
//       // 	}
//       // 	nowP[i] /= ScalorType(nDir[i]);
//       // 	mu [i] = powf (1.f + dt / tauP[i] * betaP[i] * (nowP[i] - refP[i]), 1.f/3.f);
//       // }
//     }
  
//     myst.clearDevice();
//     lpfrog.stepV (sys, dt, myst);
//     if (ptr_thermostat != NULL){
//       rescaleProperty <<<atomGridDim, myBlockDim>>>(
// 	  sys.ddata.velox, sys.ddata.numAtom,
// 	  lambda);
//       rescaleProperty <<<atomGridDim, myBlockDim>>>(
// 	  sys.ddata.veloy, sys.ddata.numAtom,
// 	  lambda);
//       rescaleProperty <<<atomGridDim, myBlockDim>>>(
// 	  sys.ddata.veloz, sys.ddata.numAtom,
// 	  lambda);
//       rescaleProperty <<<1, 3>>>(
// 	  myst.ddata, mdStatisticKineticEnergyXX, 3,
// 	  lambda * lambda);
//     }
//     lpfrog.stepX (sys, dt);
//     if (ptr_barostat != NULL){
//       ScalorType newBoxX(sys.box.size.x);
//       ScalorType newBoxY(sys.box.size.y);
//       ScalorType newBoxZ(sys.box.size.z);
//       newBoxX *= mu[0];
//       newBoxY *= mu[1];
//       newBoxZ *= mu[2];
//       CoordType coordScalor ;
//       coordScalor.x = mu[0];
//       coordScalor.y = mu[1];
//       coordScalor.z = mu[2];      
//       // for (IndexType i = 0; i < NPCoupleGroup; ++i){
//       // 	if ((PCoupleDirections[i] & PCoupleX) != 0){
//       // 	  coordScalor.x *= mu[i];
//       // 	  newBoxX *= mu[i];
//       // 	}
//       // 	if ((PCoupleDirections[i] & PCoupleY) != 0){
//       // 	  coordScalor.y *= mu[i];
//       // 	  newBoxY *= mu[i];
//       // 	}
//       // 	if ((PCoupleDirections[i] & PCoupleZ) != 0){
//       // 	  coordScalor.z *= mu[i];
//       // 	  newBoxZ *= mu[i];
//       // 	}
//       // }
//       rescaleCoord <<<atomGridDim, myBlockDim>>> (
//     	  sys.ddata.coord, sys.ddata.numAtom,
//     	  coordScalor);
//       sys.setBoxSize (newBoxX, newBoxY, newBoxZ);
//     }
//     nstep ++;
//     if (timer != NULL) timer->toc (mdTimeIntegrator);
//     if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
//       ptr_nlist->build(sys, timer);
//     }
//     ptr_inter->clearInteraction (sys);
//     if (ptr_nlist != NULL){
//       ptr_inter->applyNonBondedInteraction (sys, *ptr_nlist, myst, timer);
//     }
//     if (ptr_bdInterList != NULL){
//       ptr_inter->applyBondedInteraction (sys, *ptr_bdInterList, myst, timer);
//     }
//   }
//   else {
//     firstStep (sys, timer);
//   }
// }

// void LeapFrog_TPCouple_Rescale::
// oneStep (MDSystem & sys,
// 	 MDStatistic &st,
// 	 MDTimer * timer)
// {
//   ScalorType nowK, lambda;
//   ScalorType nowP[3], mu[3];
//   // IndexType nDir[3];
  
//   if (timer != NULL) timer->tic (mdTimeIntegrator);
//   if (nstep != 0) {
//     myst.updateHost();
//     if (ptr_thermostat != NULL){
//       nowK = myst.kineticEnergy();
//       lambda = ptr_thermostat->calScale (nowK);
//     }
//     if (ptr_barostat != 0){
//       nowP[0] = myst.pressureXX (sys.box);
//       nowP[1] = myst.pressureYY (sys.box);
//       nowP[2] = myst.pressureZZ (sys.box);
//       ptr_barostat->calScale (nowP, mu);
//     }
  
//     myst.clearDevice();
//     lpfrog.stepV (sys, dt, myst);
//     if (ptr_thermostat != NULL){
//       rescaleProperty <<<atomGridDim, myBlockDim>>>(
// 	  sys.ddata.velox, sys.ddata.numAtom,
// 	  lambda);
//       rescaleProperty <<<atomGridDim, myBlockDim>>>(
// 	  sys.ddata.veloy, sys.ddata.numAtom,
// 	  lambda);
//       rescaleProperty <<<atomGridDim, myBlockDim>>>(
// 	  sys.ddata.veloz, sys.ddata.numAtom,
// 	  lambda);
//       rescaleProperty <<<1, 3>>>(
// 	  myst.ddata, mdStatisticKineticEnergyXX, 3,
// 	  lambda * lambda);
//     }
//     lpfrog.stepX (sys, dt);
//     if (ptr_barostat != NULL){
//       ScalorType newBoxX(sys.box.size.x);
//       ScalorType newBoxY(sys.box.size.y);
//       ScalorType newBoxZ(sys.box.size.z);
//       newBoxX *= mu[0];
//       newBoxY *= mu[1];
//       newBoxZ *= mu[2];
//       CoordType coordScalor ;
//       coordScalor.x = mu[0];
//       coordScalor.y = mu[1];
//       coordScalor.z = mu[2];      
//       rescaleCoord <<<atomGridDim, myBlockDim>>> (
//     	  sys.ddata.coord, sys.ddata.numAtom,
//     	  coordScalor);
//       sys.setBoxSize (newBoxX, newBoxY, newBoxZ);
//     }
//     nstep ++;
//     if (timer != NULL) timer->toc (mdTimeIntegrator);
//     if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
//       // printf("# rebuild at step %d\n", nstep);
//       // fflush(stdout);
//       ptr_nlist->build(sys, timer);
//     }
//     ptr_inter->clearInteraction (sys);
//     if (ptr_nlist != NULL){
//       ptr_inter->applyNonBondedInteraction (sys, *ptr_nlist, myst, timer);
//     }
//     if (ptr_bdInterList != NULL){
//       ptr_inter->applyBondedInteraction (sys, *ptr_bdInterList, myst, timer);
//     }
//     st.deviceAdd (myst);
//   }
//   else {
//     firstStep (sys, st, timer);
//   }
// }


void LeapFrog_TPCouple_Rescale::
addThermostat (const Thermostat_VRescale & thermostat)
{
  ptr_thermostat = &thermostat;
}

void LeapFrog_TPCouple_Rescale::
disableThermostat ()
{
  ptr_thermostat = NULL;
}

void LeapFrog_TPCouple_Rescale::
addBarostat (const Barostat_XRescale & barostat)
{
  ptr_barostat = &barostat;
}

void LeapFrog_TPCouple_Rescale::
disableBarostat ()
{
  ptr_barostat = NULL;
}










  
LeapFrog_TPCouple_VCouple::
LeapFrog_TPCouple_VCouple ()
{
  ptr_thermostat = NULL;
  ptr_barostat = NULL;
  isFirstStep = true;
}

LeapFrog_TPCouple_VCouple::
LeapFrog_TPCouple_VCouple (const MDSystem &sys,
			   const IndexType & NThread)
{
  reinit (sys, NThread);
}

void LeapFrog_TPCouple_VCouple::
reinit (const MDSystem &		sys,
	const IndexType &		NThread)
{
  ptr_thermostat = NULL;
  ptr_barostat = NULL;
  isFirstStep = true;

  lpfrog.reinit (sys, NThread);

  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NThread;
  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);
}

void LeapFrog_TPCouple_VCouple::
reset ()
{
  isFirstStep = true;
}


void LeapFrog_TPCouple_VCouple::
firstStep (MDSystem & sys,
	   const ScalorType & dt,
	   MDStatistic &thisStepSt,
	   MDTimer * timer)
{
  lpfrog.step (sys, dt, thisStepSt, timer);
  isFirstStep = false;
}

// void LeapFrog_TPCouple_VCouple::
// firstStep (MDSystem & sys,
// 	   MDTimer * timer)
// {
//   myst.clearDevice ();
//   ptr_inter->clearInteraction (sys);
//   if (ptr_nlist != NULL){
//     ptr_inter->applyNonBondedInteraction (sys, *ptr_nlist, timer);
//   }
//   if (ptr_bdInterList != NULL){
//     ptr_inter->applyBondedInteraction (sys, *ptr_bdInterList, timer);
//   }
//   lpfrog.step (sys, dt, myst, timer);
//   if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
//     // printf("# rebuild at step %d\n", nstep);
//     // fflush(stdout);
//     ptr_nlist->build(sys, timer);
//   }
//   ptr_inter->clearInteraction (sys);
//   if (ptr_nlist != NULL){
//     ptr_inter->applyNonBondedInteraction (sys, *ptr_nlist, myst, timer);
//   }
//   if (ptr_bdInterList != NULL){
//     ptr_inter->applyBondedInteraction (sys, *ptr_bdInterList, myst, timer);
//   }
//   nstep ++;
// }


// void LeapFrog_TPCouple_VCouple::
// oneStep (MDSystem & sys,
// 	 MDTimer * timer)
// {
//   // printf ("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
//   ScalorType nowK;
//   ScalorType lambda[3];
//   ScalorType nowP[3], mu[3];
  
//   lambda[2] = lambda[1] = lambda[0] = 0.f;

//   if (timer != NULL) timer->tic (mdTimeIntegrator);
//   if (nstep != 0) {
//     myst.updateHost();
//     if (ptr_thermostat != NULL){
//       nowK = myst.kineticEnergy();
//       lambda[0] = ptr_thermostat->calCouple (nowK);
//       lambda[2] = lambda[1] = lambda[0];
//     }
//     if (ptr_barostat != NULL){
//       RectangularBox tmpBox;
//       ScalorType tmpLambda[3];
//       nowP[0] = myst.pressureXX (sys.box);
//       nowP[1] = myst.pressureYY (sys.box);
//       nowP[2] = myst.pressureZZ (sys.box);
//       ptr_barostat->calCouple (nowP, tmpLambda, tmpBox);
//       lambda[0] += tmpLambda[0];
//       lambda[1] += tmpLambda[1];
//       lambda[2] += tmpLambda[2];
//       mu[0] = tmpBox.size.x / sys.box.size.x;
//       mu[1] = tmpBox.size.y / sys.box.size.y;
//       mu[2] = tmpBox.size.z / sys.box.size.z;
//       sys.box = (tmpBox);
//     }
  
//     myst.clearDevice();
//     lpfrog.stepV_VCouple (sys, dt, lambda, myst);
//     lpfrog.stepX (sys, dt);
//     if (ptr_barostat != NULL){
//       CoordType coordScalor ;
//       coordScalor.x = mu[0];
//       coordScalor.y = mu[1];
//       coordScalor.z = mu[2];      
//       rescaleCoord <<<atomGridDim, myBlockDim>>> (
//     	  sys.ddata.coord, sys.ddata.numAtom,
//     	  coordScalor);
//       rescaleProperty <<<atomGridDim, myBlockDim>>>(
// 	  sys.ddata.velox, sys.ddata.numAtom,
// 	  mu[0]);
//       rescaleProperty <<<atomGridDim, myBlockDim>>>(
// 	  sys.ddata.veloy, sys.ddata.numAtom,
// 	  mu[1]);
//       rescaleProperty <<<atomGridDim, myBlockDim>>>(
// 	  sys.ddata.veloz, sys.ddata.numAtom,
// 	  mu[2]);
//       rescaleProperty <<<1, 1>>>(
// 	  myst.ddata, mdStatisticKineticEnergyXX, 1,
// 	  mu[0] * mu[0]);
//       rescaleProperty <<<1, 1>>>(
// 	  myst.ddata, mdStatisticKineticEnergyYY, 1,
// 	  mu[1] * mu[1]);
//       rescaleProperty <<<1, 1>>>(
// 	  myst.ddata, mdStatisticKineticEnergyZZ, 1,
// 	  mu[2] * mu[2]);
//     }
//     nstep ++;
//     if (timer != NULL) timer->toc (mdTimeIntegrator);
//     if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
//       ptr_nlist->build(sys, timer);
//     }
//     ptr_inter->clearInteraction (sys);
//     if (ptr_nlist != NULL){
//       ptr_inter->applyNonBondedInteraction (sys, *ptr_nlist, myst, timer);
//     }
//     if (ptr_bdInterList != NULL){
//       ptr_inter->applyBondedInteraction (sys, *ptr_bdInterList, myst, timer);
//     }
//   }
//   else {
//     firstStep (sys, timer);
//   }
// }

void LeapFrog_TPCouple_VCouple::
oneStep (MDSystem &		sys,
	 const ScalorType &	dt,
	 MDStatistic &		lastStepSt,
	 MDStatistic &		thisStepSt,
	 MDTimer *		timer)
{
  ScalorType nowK;
  ScalorType lambda[3];
  ScalorType nowP[3], mu[3];
  
  lambda[2] = lambda[1] = lambda[0] = 0.f;
  
  if (! isFirstStep ) {
    if (timer != NULL) timer->tic (mdTimeIntegrator);
    lastStepSt.updateHost();
    if (ptr_thermostat != NULL){
      nowK = lastStepSt.kineticEnergy();
      lambda[2] = ptr_thermostat->calCouple (nowK);
      lambda[0] = lambda[1] = lambda[2];
      // printf ("lambda %f %f %f\n", lambda[0], lambda[1], lambda[2]);
    }
    if (ptr_barostat != NULL){
      RectangularBox tmpBox;
      ScalorType tmpLambda[3];
      nowP[0] = lastStepSt.pressureXX (sys.box);
      nowP[1] = lastStepSt.pressureYY (sys.box);
      nowP[2] = lastStepSt.pressureZZ (sys.box);
      ptr_barostat->calCouple (nowP, tmpLambda, tmpBox);
      lambda[0] += tmpLambda[0];
      lambda[1] += tmpLambda[1];
      lambda[2] += tmpLambda[2];
      mu[0] = tmpBox.size.x / sys.box.size.x;
      mu[1] = tmpBox.size.y / sys.box.size.y;
      mu[2] = tmpBox.size.z / sys.box.size.z;
      sys.box = (tmpBox);
    }
  
    lpfrog.stepV_VCouple (sys, dt, lambda, thisStepSt);
    lpfrog.stepX (sys, dt);
    if (ptr_barostat != NULL){
      CoordType coordScalor ;
      coordScalor.x = mu[0];
      coordScalor.y = mu[1];
      coordScalor.z = mu[2];      
      rescaleCoord <<<atomGridDim, myBlockDim>>> (
    	  sys.ddata.coord, sys.ddata.numAtom,
    	  coordScalor);
      rescaleProperty <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.velox, sys.ddata.numAtom,
	  mu[0]);
      rescaleProperty <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.veloy, sys.ddata.numAtom,
	  mu[1]);
      rescaleProperty <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.veloz, sys.ddata.numAtom,
	  mu[2]);
      rescaleProperty <<<1, 1>>>(
	  thisStepSt.ddata, mdStatisticKineticEnergyXX, 1,
	  mu[0] * mu[0]);
      rescaleProperty <<<1, 1>>>(
	  thisStepSt.ddata, mdStatisticKineticEnergyYY, 1,
	  mu[1] * mu[1]);
      rescaleProperty <<<1, 1>>>(
	  thisStepSt.ddata, mdStatisticKineticEnergyZZ, 1,
	  mu[2] * mu[2]);
    }
    if (timer != NULL) timer->toc (mdTimeIntegrator);
  }
  else {
    firstStep (sys, dt, thisStepSt, timer);
  }
}


void LeapFrog_TPCouple_VCouple::
addThermostat (const Thermostat_VCouple & thermostat)
{
  ptr_thermostat = &thermostat;
}

void LeapFrog_TPCouple_VCouple::
disableThermostat ()
{
  ptr_thermostat = NULL;
}

void LeapFrog_TPCouple_VCouple::
addBarostat (const Barostat_VCouple & barostat)
{
  ptr_barostat = &barostat;
}

void LeapFrog_TPCouple_VCouple::
disableBarostat ()
{
  ptr_barostat = NULL;
}
  

