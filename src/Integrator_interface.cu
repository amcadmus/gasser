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
  cudaFree (buffx);
  cudaFree (buffy);
  cudaFree (buffz);
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

  IndexType buffsize = nob * sizeof(ScalorType);
  cudaMalloc ((void**)&buffx, buffsize);
  cudaMalloc ((void**)&buffy, buffsize);
  cudaMalloc ((void**)&buffz, buffsize);
  cudaMalloc ((void**)&sums,  3 * sizeof(ScalorType));
  checkCUDAError ("TranslationalFreedomRemover::init allocation");

  fillSums0 <<<1,1>>> (sums);
  initRemoveTranslationalFreedom <<<1,1>>> ();
  checkCUDAError ("TranslationalFreedomRemover::init");
}

void TranslationalFreedomRemover::remove (MDSystem & sys,
					  MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeRemoveTransFreedom);
  prepareRemoveTranslationalFreedom
      <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.numAtom, sys.ddata.mass,
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  buffx, buffy, buffz, sums);
  checkCUDAError("TranslationalFreedomRemover::remove, prepare");
  removeFreedom 
      <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.numAtom,
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz, sys.ddata.totalMassi, sums);
  checkCUDAError("TranslationalFreedomRemover::remove, remove");
  if (timer != NULL) timer->toc(mdTimeRemoveTransFreedom);
}
  

LeapFrog::~LeapFrog()
{
}

void LeapFrog::init (const MDSystem &sys,
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

  sum_kxx.init (nob, NThreadForSum);
  sum_kyy.init (nob, NThreadForSum);
  sum_kzz.init (nob, NThreadForSum);  

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
#ifndef COORD_IN_ONE_VEC
	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
#else
	  sys.ddata.coord,
#endif
	  sys.ddata.velox, sys.ddata.veloy, sys.ddata.veloz,
	  sys.ddata.forcx, sys.ddata.forcy, sys.ddata.forcz,
	  dt,
	  sum_kxx.getBuff(),
	  sum_kyy.getBuff(),
	  sum_kzz.getBuff());
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
	  sum_kxx.getBuff(),
	  sum_kyy.getBuff(),
	  sum_kzz.getBuff());
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

  sum_kxx.init (nob, NThreadForSum);
  sum_kyy.init (nob, NThreadForSum);
  sum_kzz.init (nob, NThreadForSum);  

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
	  sum_kxx.getBuff(),
	  sum_kyy.getBuff(),
	  sum_kzz.getBuff());
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
  Nf = sys.ddata.NFreedom;
  scalor1 = 1./tau;
  scalor2 = 1. / sqrt(Nf) / sqrt(tau) * 2; 

  refK = (Nf-3) * 0.5 * refT;
  printf ("# refK is %f\n", refK);
  
  cudaMalloc ((void**)&kineticE, sizeof(ScalorType) * 2);
  cudaMalloc ((void**)&buff, sizeof(ScalorType) * nob);
  checkCUDAError ("VelocityVerlet::init allocation");

  sum_kxx.init (nob, NThreadForSum);
  sum_kyy.init (nob, NThreadForSum);
  sum_kzz.init (nob, NThreadForSum);
  sum_k.init (nob, NThreadForSum);
  
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
	  sum_k.getBuff());
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
	  sum_k.getBuff());
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
	  sum_kxx.getBuff(),
	  sum_kyy.getBuff(),
	  sum_kzz.getBuff());
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
  nstep = 0;
}


void
BerendsenLeapFrog::init (const MDSystem &sys,
			 const IndexType & NThread,
			 const ScalorType & dt_,
			 InteractionEngine_interface &inter,
			 NeighborList & nlist,
			 const ScalorType & rebt)
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

  lpfrog.init (sys, NThread);
  ptr_inter = & inter;
  ptr_nlist = & nlist;
  rebuildThreshold = rebt;
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

void
BerendsenLeapFrog::firstStep (MDSystem & sys, MDStatistic &st, MDTimer * timer)
{
  myst.clearDevice ();
  ptr_inter->applyInteraction (sys, *ptr_nlist, timer);
  lpfrog.step (sys, dt, myst, timer);
  if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
    // printf("# rebuild at step %d\n", nstep);
    // fflush(stdout);
    ptr_nlist->reBuild(sys, timer);
  }
  ptr_inter->applyInteraction (sys, *ptr_nlist, myst, timer);
  st.deviceAdd (myst);
  nstep ++;
}

void
BerendsenLeapFrog::firstStep (MDSystem & sys,  MDTimer * timer)
{
  myst.clearDevice ();
  ptr_inter->applyInteraction (sys, *ptr_nlist, timer);
  lpfrog.step (sys, dt, myst, timer);
  if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
    // printf("# rebuild at step %d\n", nstep);
    // fflush(stdout);
    ptr_nlist->reBuild(sys, timer);
  }
  ptr_inter->applyInteraction (sys, *ptr_nlist, myst, timer);
  nstep ++;
}


#ifndef COORD_IN_ONE_VEC
void
BerendsenLeapFrog::oneStep (MDSystem & sys, MDTimer * timer)
{
  ScalorType nowT, lambda;
  ScalorType nowP[3], mu[3];
  IndexType nDir[3];
  
  if (timer != NULL) timer->tic (mdTimeIntegrator);
  if (nstep != 0) {
    myst.updateHost();
    if (TCoupleOn){
      nowT = myst.kineticEnergy();
      nowT *= 2.f / (sys.ddata.NFreedom - 3);
      lambda = sqrtf(1.f + dt / tauT * (refT / nowT - 1.f));
    }
    if (PCoupleOn){
      for (IndexType i = 0; i < NPCoupleGroup; ++i){
	nowP[i] = 0;
	nDir[i] = 0;
	if ((PCoupleDirections[i] & PCoupleX) != 0){
	  nowP[i] += myst.pressureXX(sys.box);
	  nDir[i] ++;
	}
	if ((PCoupleDirections[i] & PCoupleY) != 0){
	  nowP[i] += myst.pressureYY(sys.box);
	  nDir[i] ++;
	}
	if ((PCoupleDirections[i] & PCoupleZ) != 0){
	  nowP[i] += myst.pressureZZ(sys.box);
	  nDir[i] ++;
	}
	nowP[i] /= ScalorType(nDir[i]);
	mu [i] = powf (1.f + dt / tauP[i] * betaP[i] * (nowP[i] - refP[i]), 1.f/3.f);
      }
    }
  
    myst.clearDevice();
    lpfrog.stepV (sys, dt, myst);
    if (TCoupleOn){
      rescaleProperty <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.velox, sys.ddata.numAtom,
	  lambda);
      rescaleProperty <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.veloy, sys.ddata.numAtom,
	  lambda);
      rescaleProperty <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.veloz, sys.ddata.numAtom,
	  lambda);
      rescaleProperty <<<1, 3>>>(
	  myst.ddata, mdStatisticKineticEnergyXX, 3,
	  lambda * lambda);
    }
    lpfrog.stepX (sys, dt);
    if (PCoupleOn){
      ScalorType newBoxX(sys.box.size.x);
      ScalorType newBoxY(sys.box.size.y);
      ScalorType newBoxZ(sys.box.size.z);
      for (IndexType i = 0; i < NPCoupleGroup; ++i){
	if ((PCoupleDirections[i] & PCoupleX) != 0){
	  rescaleProperty <<<atomGridDim, myBlockDim>>> (
	      sys.ddata.coordx, sys.ddata.numAtom,
	      mu[i]);
	  newBoxX *= mu[i];
	}
	if ((PCoupleDirections[i] & PCoupleY) != 0){
	  rescaleProperty <<<atomGridDim, myBlockDim>>> (
	      sys.ddata.coordy, sys.ddata.numAtom,
	      mu[i]);
	  newBoxY *= mu[i];
	}
	if ((PCoupleDirections[i] & PCoupleZ) != 0){
	  rescaleProperty <<<atomGridDim, myBlockDim>>> (
	      sys.ddata.coordz, sys.ddata.numAtom,
	      mu[i]);
	  newBoxZ *= mu[i];
	}
	sys.setBoxSize (newBoxX, newBoxY, newBoxZ);
      }
    }
    nstep ++;
    if (timer != NULL) timer->toc (mdTimeIntegrator);
    if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
      // printf("# rebuild at step %d\n", nstep);
      // fflush(stdout);
      ptr_nlist->reBuild(sys, timer);
    }
    ptr_inter->applyInteraction (sys, *ptr_nlist, myst, timer);
  }
  else {
    firstStep (sys, timer);
  }
}
void
BerendsenLeapFrog::oneStep (MDSystem & sys, MDStatistic &st, MDTimer * timer)
{
  ScalorType nowT, lambda;
  ScalorType nowP[3], mu[3];
  IndexType nDir[3];
  
  if (timer != NULL) timer->tic (mdTimeIntegrator);
  if (nstep != 0) {
    myst.updateHost();
    if (TCoupleOn){
      nowT = myst.kineticEnergy();
      nowT *= 2.f / (sys.ddata.NFreedom - 3);
      lambda = sqrtf(1.f + dt / tauT * (refT / nowT - 1.f));
    }
    if (PCoupleOn){
      for (IndexType i = 0; i < NPCoupleGroup; ++i){
	nowP[i] = 0;
	nDir[i] = 0;
	if ((PCoupleDirections[i] & PCoupleX) != 0){
	  nowP[i] += myst.pressureXX(sys.box);
	  nDir[i] ++;
	}
	if ((PCoupleDirections[i] & PCoupleY) != 0){
	  nowP[i] += myst.pressureYY(sys.box);
	  nDir[i] ++;
	}
	if ((PCoupleDirections[i] & PCoupleZ) != 0){
	  nowP[i] += myst.pressureZZ(sys.box);
	  nDir[i] ++;
	}
	nowP[i] /= ScalorType(nDir[i]);
	mu [i] = powf (1.f + dt / tauP[i] * betaP[i] * (nowP[i] - refP[i]), 1.f/3.f);
      }
    }
  
    myst.clearDevice();
    lpfrog.stepV (sys, dt, myst);
    if (TCoupleOn){
      rescaleProperty <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.velox, sys.ddata.numAtom,
	  lambda);
      rescaleProperty <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.veloy, sys.ddata.numAtom,
	  lambda);
      rescaleProperty <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.veloz, sys.ddata.numAtom,
	  lambda);
      rescaleProperty <<<1, 3>>>(
	  myst.ddata, mdStatisticKineticEnergyXX, 3,
	  lambda * lambda);
    }
    lpfrog.stepX (sys, dt);
    if (PCoupleOn){
      ScalorType newBoxX(sys.box.size.x);
      ScalorType newBoxY(sys.box.size.y);
      ScalorType newBoxZ(sys.box.size.z);
      for (IndexType i = 0; i < NPCoupleGroup; ++i){
	if ((PCoupleDirections[i] & PCoupleX) != 0){
	  rescaleProperty <<<atomGridDim, myBlockDim>>> (
	      sys.ddata.coordx, sys.ddata.numAtom,
	      mu[i]);
	  newBoxX *= mu[i];
	}
	if ((PCoupleDirections[i] & PCoupleY) != 0){
	  rescaleProperty <<<atomGridDim, myBlockDim>>> (
	      sys.ddata.coordy, sys.ddata.numAtom,
	      mu[i]);
	  newBoxY *= mu[i];
	}
	if ((PCoupleDirections[i] & PCoupleZ) != 0){
	  rescaleProperty <<<atomGridDim, myBlockDim>>> (
	      sys.ddata.coordz, sys.ddata.numAtom,
	      mu[i]);
	  newBoxZ *= mu[i];
	}
	sys.setBoxSize (newBoxX, newBoxY, newBoxZ);
      }
    }
    nstep ++;
    if (timer != NULL) timer->toc (mdTimeIntegrator);
    if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
      // printf("# rebuild at step %d\n", nstep);
      // fflush(stdout);
      ptr_nlist->reBuild(sys, timer);
    }
    ptr_inter->applyInteraction (sys, *ptr_nlist, myst, timer);
    st.deviceAdd (myst);
  }
  else {
    firstStep (sys, st, timer);
  }
}
#else
void
BerendsenLeapFrog::oneStep (MDSystem & sys, MDTimer * timer)
{
  ScalorType nowT, lambda;
  ScalorType nowP[3], mu[3];
  IndexType nDir[3];
  
  if (timer != NULL) timer->tic (mdTimeIntegrator);
  if (nstep != 0) {
    myst.updateHost();
    if (TCoupleOn){
      nowT = myst.kineticEnergy();
      nowT *= 2.f / (sys.ddata.NFreedom - 3);
      lambda = sqrtf(1.f + dt / tauT * (refT / nowT - 1.f));
    }
    if (PCoupleOn){
      for (IndexType i = 0; i < NPCoupleGroup; ++i){
	nowP[i] = 0;
	nDir[i] = 0;
	if ((PCoupleDirections[i] & PCoupleX) != 0){
	  nowP[i] += myst.pressureXX(sys.box);
	  nDir[i] ++;
	}
	if ((PCoupleDirections[i] & PCoupleY) != 0){
	  nowP[i] += myst.pressureYY(sys.box);
	  nDir[i] ++;
	}
	if ((PCoupleDirections[i] & PCoupleZ) != 0){
	  nowP[i] += myst.pressureZZ(sys.box);
	  nDir[i] ++;
	}
	nowP[i] /= ScalorType(nDir[i]);
	mu [i] = powf (1.f + dt / tauP[i] * betaP[i] * (nowP[i] - refP[i]), 1.f/3.f);
      }
    }
  
    myst.clearDevice();
    lpfrog.stepV (sys, dt, myst);
    if (TCoupleOn){
      rescaleProperty <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.velox, sys.ddata.numAtom,
	  lambda);
      rescaleProperty <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.veloy, sys.ddata.numAtom,
	  lambda);
      rescaleProperty <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.veloz, sys.ddata.numAtom,
	  lambda);
      rescaleProperty <<<1, 3>>>(
	  myst.ddata, mdStatisticKineticEnergyXX, 3,
	  lambda * lambda);
    }
    lpfrog.stepX (sys, dt);
    if (PCoupleOn){
      ScalorType newBoxX(sys.box.size.x);
      ScalorType newBoxY(sys.box.size.y);
      ScalorType newBoxZ(sys.box.size.z);
      CoordType coordScalor ;
      coordScalor.x = 1.f;
      coordScalor.y = 1.f;
      coordScalor.z = 1.f;
      for (IndexType i = 0; i < NPCoupleGroup; ++i){
	if ((PCoupleDirections[i] & PCoupleX) != 0){
	  coordScalor.x *= mu[i];
	  newBoxX *= mu[i];
	}
	if ((PCoupleDirections[i] & PCoupleY) != 0){
	  coordScalor.y *= mu[i];
	  newBoxY *= mu[i];
	}
	if ((PCoupleDirections[i] & PCoupleZ) != 0){
	  coordScalor.z *= mu[i];
	  newBoxZ *= mu[i];
	}
      }
      rescaleCoord <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.coord, sys.ddata.numAtom,
	  coordScalor);
      sys.setBoxSize (newBoxX, newBoxY, newBoxZ);
    }
    nstep ++;
    if (timer != NULL) timer->toc (mdTimeIntegrator);
    if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
      // printf("# rebuild at step %d\n", nstep);
      // fflush(stdout);
      ptr_nlist->reBuild(sys, timer);
    }
    ptr_inter->applyInteraction (sys, *ptr_nlist, myst, timer);
  }
  else {
    firstStep (sys, timer);
  }
}

void
BerendsenLeapFrog::oneStep (MDSystem & sys, MDStatistic &st, MDTimer * timer)
{
  ScalorType nowT, lambda;
  ScalorType nowP[3], mu[3];
  IndexType nDir[3];
  
  if (timer != NULL) timer->tic (mdTimeIntegrator);
  if (nstep != 0) {
    myst.updateHost();
    if (TCoupleOn){
      nowT = myst.kineticEnergy();
      nowT *= 2.f / (sys.ddata.NFreedom - 3);
      lambda = sqrtf(1.f + dt / tauT * (refT / nowT - 1.f));
    }
    if (PCoupleOn){
      for (IndexType i = 0; i < NPCoupleGroup; ++i){
	nowP[i] = 0;
	nDir[i] = 0;
	if ((PCoupleDirections[i] & PCoupleX) != 0){
	  nowP[i] += myst.pressureXX(sys.box);
	  nDir[i] ++;
	}
	if ((PCoupleDirections[i] & PCoupleY) != 0){
	  nowP[i] += myst.pressureYY(sys.box);
	  nDir[i] ++;
	}
	if ((PCoupleDirections[i] & PCoupleZ) != 0){
	  nowP[i] += myst.pressureZZ(sys.box);
	  nDir[i] ++;
	}
	nowP[i] /= ScalorType(nDir[i]);
	mu [i] = powf (1.f + dt / tauP[i] * betaP[i] * (nowP[i] - refP[i]), 1.f/3.f);
      }
    }
  
    myst.clearDevice();
    lpfrog.stepV (sys, dt, myst);
    if (TCoupleOn){
      rescaleProperty <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.velox, sys.ddata.numAtom,
	  lambda);
      rescaleProperty <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.veloy, sys.ddata.numAtom,
	  lambda);
      rescaleProperty <<<atomGridDim, myBlockDim>>>(
	  sys.ddata.veloz, sys.ddata.numAtom,
	  lambda);
      rescaleProperty <<<1, 3>>>(
	  myst.ddata, mdStatisticKineticEnergyXX, 3,
	  lambda * lambda);
    }
    lpfrog.stepX (sys, dt);
    if (PCoupleOn){
      ScalorType newBoxX(sys.box.size.x);
      ScalorType newBoxY(sys.box.size.y);
      ScalorType newBoxZ(sys.box.size.z);
      CoordType coordScalor ;
      coordScalor.x = 1.f;
      coordScalor.y = 1.f;
      coordScalor.z = 1.f;
      for (IndexType i = 0; i < NPCoupleGroup; ++i){
	if ((PCoupleDirections[i] & PCoupleX) != 0){
	  coordScalor.x *= mu[i];
	  newBoxX *= mu[i];
	}
	if ((PCoupleDirections[i] & PCoupleY) != 0){
	  coordScalor.y *= mu[i];
	  newBoxY *= mu[i];
	}
	if ((PCoupleDirections[i] & PCoupleZ) != 0){
	  coordScalor.z *= mu[i];
	  newBoxZ *= mu[i];
	}
      }
      rescaleCoord <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.coord, sys.ddata.numAtom,
	  coordScalor);
      sys.setBoxSize (newBoxX, newBoxY, newBoxZ);
    }
    nstep ++;
    if (timer != NULL) timer->toc (mdTimeIntegrator);
    if (ptr_nlist->judgeRebuild (sys, rebuildThreshold, timer)){
      // printf("# rebuild at step %d\n", nstep);
      // fflush(stdout);
      ptr_nlist->reBuild(sys, timer);
    }
    ptr_inter->applyInteraction (sys, *ptr_nlist, myst, timer);
    st.deviceAdd (myst);
  }
  else {
    firstStep (sys, st, timer);
  }
}
#endif


