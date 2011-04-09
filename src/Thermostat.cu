#define DEVICE_CODE

#include "Thermostat.h"
#include "MDSystem_interface.h"
#include <cmath>

#include "compile_error_mixcode.h"

NoseHoover_Chains2::
NoseHoover_Chains2 ()
{
  cudaMalloc ((void**)&dK, sizeof(ScalorType));
  checkCUDAError ("NoseHoover_Chains2::operator_L_G1, malloc dK");
}

NoseHoover_Chains2::
~NoseHoover_Chains2 ()
{
  cudaFree (dK);
  checkCUDAError ("NoseHoover_Chains2::operator_L_G1, free dK");
}

void NoseHoover_Chains2::
reinit (const MDSystem &sys,
	const IndexType & NThread,
	const ScalorType & ref_T,
	const ScalorType & tau_T)
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
  sum_k.reinit   (nob, NThreadForSum);
  
  sharedBuffSize = NThread * 2 * sizeof(ScalorType);

  LL = 3.f * sys.ddata.numAtom;
  Q2 = Q1 = tau_T * tau_T * ref_T / (4 * M_PI * M_PI) * LL;
  xi1 = xi2 = 1.f;
  vxi1 = vxi2 = 0.f;
  refT = ref_T;
}


void NoseHoover_Chains2::
operator_L_xi (const ScalorType & dt)
{
  xi1 += vxi1 * dt;
  xi2 += vxi2 * dt;
}

static void __global__
NHC_L_Cv (ScalorType a,
	  IndexType numAtom,
	  ScalorType * velox,
	  ScalorType * veloy,
	  ScalorType * veloz,
	  ScalorType * mass,
	  ScalorType * st_buffx,
	  ScalorType * st_buffy,
	  ScalorType * st_buffz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  ScalorType vx, vy, vz;
  
  if (ii < numAtom) {
    velox[ii] *= a;
    veloy[ii] *= a;
    veloz[ii] *= a;
    vx = velox[ii];
    vy = veloy[ii];
    vz = veloz[ii];
  }
  extern __shared__ ScalorType buff [];

  ScalorType scalor;
  if (ii < numAtom) scalor = 0.5f * mass[ii];
  else scalor = 0.f;

  if (ii < numAtom){
    buff[threadIdx.x] = scalor * vx * vx;
  }
  else {
    buff[threadIdx.x] = 0.f;
  }
  __syncthreads();
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) st_buffx[bid] = buff[0];
  __syncthreads();

  if (ii < numAtom){
    buff[threadIdx.x] = scalor * vy * vy;
  }
  __syncthreads();
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) st_buffy[bid] = buff[0];
  __syncthreads();

  if (ii < numAtom){
    buff[threadIdx.x] = scalor * vz * vz;
  }
  __syncthreads();
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) st_buffz[bid] = buff[0];
}


static void __global__
NHC_L_Cv (ScalorType a,
	  IndexType numAtom,
	  ScalorType * velox,
	  ScalorType * veloy,
	  ScalorType * veloz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  
  if (ii < numAtom) {
    velox[ii] *= a;
    veloy[ii] *= a;
    veloz[ii] *= a;
  }
}


void NoseHoover_Chains2::
operator_L_Cv (const ScalorType & dt,
	       MDSystem & sys,
	       MDStatistic & st)
{
  ScalorType a = expf (-vxi1 * dt);
  NHC_L_Cv
      <<<atomGridDim, myBlockDim, sharedBuffSize>>> (
	  a,
	  sys.ddata.numAtom,
	  sys.ddata.velox,
	  sys.ddata.veloy,
	  sys.ddata.veloz,
	  sys.ddata.mass,
	  sum_kxx.buff,
	  sum_kyy.buff,
	  sum_kzz.buff);
  checkCUDAError ("NoseHoover_Chains2::operator_L_Cv NHC_L_Cv");
  sum_kxx.sumBuffAdd (st.ddata, mdStatisticKineticEnergyXX, 0);
  sum_kyy.sumBuffAdd (st.ddata, mdStatisticKineticEnergyYY, 0);
  sum_kzz.sumBuffAdd (st.ddata, mdStatisticKineticEnergyZZ, 0);
  checkCUDAError ("NoseHoover_Chains2::operator_L_Cv sum");
}

void NoseHoover_Chains2::
operator_L_Cv (const ScalorType & dt,
	       MDSystem & sys)
{
  ScalorType a = expf (-vxi1 * dt);
  NHC_L_Cv
      <<<atomGridDim, myBlockDim>>> (
	  a,
	  sys.ddata.numAtom,
	  sys.ddata.velox,
	  sys.ddata.veloy,
	  sys.ddata.veloz);
}


static void __global__
calKinetic (IndexType numAtom,
	    ScalorType * velox,
	    ScalorType * veloy,
	    ScalorType * veloz,
	    ScalorType * mass,
	    ScalorType * st_buff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  ScalorType tmp;
  if (ii < numAtom) {
    tmp = 0.5f * mass[ii] * (
	velox[ii] * velox[ii] +
	veloy[ii] * veloy[ii] +
	veloz[ii] * veloz[ii]);
  }

  extern __shared__ ScalorType buff [];
  if (ii < numAtom){
    buff[threadIdx.x] = tmp;
  }
  else {
    buff[threadIdx.x] = 0.f;
  }
  __syncthreads();
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) st_buff[bid] = buff[0];
  __syncthreads();
}

void NoseHoover_Chains2::	
operator_L_G1 (const ScalorType & dt,
	       const MDSystem & sys)
{
  ScalorType G1;  
  calKinetic
      <<<atomGridDim, myBlockDim, sharedBuffSize>>> (
	  sys.ddata.numAtom,
	  sys.ddata.velox,
	  sys.ddata.veloy,
	  sys.ddata.veloz,
	  sys.ddata.mass,
	  sum_k.buff);
  sum_k.sumBuff (dK, 0, 0);
  checkCUDAError ("NoseHoover_Chains2::operator_L_G1, calKinetic");
  cudaMemcpy (&hK, dK, sizeof(ScalorType), cudaMemcpyDeviceToHost);
  checkCUDAError ("NoseHoover_Chains2::operator_L_G1, cpy dK");
  G1 = (2.f * hK - LL * refT) / Q1;

  vxi1 += G1 * dt;
}


void NoseHoover_Chains2::
operator_L_vxi1 (const ScalorType & dt)
{
  vxi1 *= expf (-vxi2 * dt);
}

void NoseHoover_Chains2::
operator_L_G2 (const ScalorType & dt)
{
  ScalorType G2 = (Q1 * vxi1 * vxi1 - refT) / Q2;
  vxi2 += G2 * dt;
}


ScalorType NoseHoover_Chains2::
HamiltonianContribution () const 
{
  return (0.5 * Q1 * vxi1 * vxi1 +
	  0.5 * Q2 * vxi2 * vxi2 +
	  LL * refT * xi1 +
	  refT * xi2 );
}

void NoseHoover_Chains2::
operator_L (const ScalorType & dt,
	    MDSystem & sys,
	    MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  operator_L_G2		(dt * 0.5);
  
  operator_L_vxi1	(dt * 0.25);
  operator_L_G1		(dt * 0.5, sys);
  operator_L_vxi1	(dt * 0.25);

  operator_L_Cv		(dt, sys);
  operator_L_xi		(dt);

  operator_L_vxi1	(dt * 0.25);
  operator_L_G1		(dt * 0.5, sys);
  operator_L_vxi1	(dt * 0.25);

  operator_L_G2		(dt * 0.5);
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}

void NoseHoover_Chains2::
operator_L (const ScalorType & dt,
	    MDSystem & sys,
	    MDStatistic & st,
	    MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  operator_L_G2		(dt * 0.5);
  
  operator_L_vxi1	(dt * 0.25);
  operator_L_G1		(dt * 0.5, sys);
  operator_L_vxi1	(dt * 0.25);

  operator_L_Cv		(dt, sys, st);
  operator_L_xi		(dt);

  operator_L_vxi1	(dt * 0.25);
  operator_L_G1		(dt * 0.5, sys);
  operator_L_vxi1	(dt * 0.25);

  operator_L_G2		(dt * 0.5);
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}




