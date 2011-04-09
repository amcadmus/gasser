#define DEVICE_CODE

#include "Barostat.h"
#include <cmath>

#include "compile_error_mixcode.h"


void NoseHoover_Chains2_Isobaric::
reinit (const MDSystem &sys,
	const IndexType & NThread,
	const ScalorType & ref_T,
	const ScalorType & tau_T,
	const ScalorType & ref_P,
	const ScalorType & tau_P,
	const ScalorType & pressureCorr)
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

  Nf = 3.f * sys.ddata.numAtom - 3.f;
  // LL = 3.f * sys.ddata.numAtom;
  // Q2 = Q1 = tau_T * tau_T * ref_T / (4 * M_PI * M_PI) * LL;
  Q2 = Q1 = tau_T * tau_T * ref_T * Nf;
  // Q2 /= 10;
  // xi1 = xi2 = 1.f;
  // vxi1 = vxi2 = 0.f;
  xi1 = 1.f;
  vxi1 = 0.f;
  refT = ref_T;

  // NN = sys.ddata.numAtom;
  // WW = tau_P * tau_P * ref_T / (4 * M_PI * M_PI) * 3.f * NN;
  WW = tau_P * tau_P * ref_T * Nf;
  ScalorType volume = sys.box.size.x * sys.box.size.y * sys.box.size.z;
  ep = log(volume) / 3.;
  vep = 0.;
  refP = ref_P;
  
  sum_kxx.reinit (nob, NThreadForSum);
  sum_kyy.reinit (nob, NThreadForSum);
  sum_kzz.reinit (nob, NThreadForSum);
  
  tmp_st.reinit (sys);
  tmp_st.setEnergyCorr (0.f);
  tmp_st.setPressureCorr (pressureCorr);
  virial_array[0] = mdStatisticVirialXX;
  virial_array[1] = mdStatisticVirialYY;
  virial_array[2] = mdStatisticVirialZZ;
  kinetic_array[0] = mdStatisticKineticEnergyXX;
  kinetic_array[1] = mdStatisticKineticEnergyYY;
  kinetic_array[2] = mdStatisticKineticEnergyZZ;
}


void NoseHoover_Chains2_Isobaric::
operator_L_ep (const ScalorType & dt)
{
  ep += vep * dt;
}

static void __global__
NHC_P_L_R (const IndexType numAtom,
	   const ScalorType scalor1,
	   const ScalorType scalor2,
	   CoordType * coord,
	   ScalorType * velox,
	   ScalorType * veloy,
	   ScalorType * veloz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (ii < numAtom) {
    coord[ii].x = scalor1 * coord[ii].x + scalor2 * velox[ii];
    coord[ii].y = scalor1 * coord[ii].y + scalor2 * veloy[ii];
    coord[ii].z = scalor1 * coord[ii].z + scalor2 * veloz[ii];
  }
}

void NoseHoover_Chains2_Isobaric::
operator_L_r (const ScalorType & dt,
	      MDSystem & sys,
	      MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  ScalorType xx = vep * dt * 0.5;
  ScalorType tmp = exp (xx);
  
  ScalorType scalor1 = exp(vep * dt);
  // ScalorType scalor2 = dt * tmp * 0.5 * (tmp - 1./tmp) / (vep * 0.5 * dt);
  // ScalorType scalor2 = dt;
  ScalorType x2 = xx * xx;
  ScalorType x4 = x2 * x2;
  ScalorType x6 = x4 * x2;
  ScalorType x8 = x6 * x2;
  ScalorType x10 = x8 * x2;
  ScalorType scalor2 = dt * tmp * (1. +
				   1./6.*x2 +
				   1./120.*x4 +
				   1./5040.*x6 + 
				   1./362880.*x8 +
				   1./39916800.*x10);
  
  NHC_P_L_R
      <<<atomGridDim, myBlockDim>>> (
	  sys.ddata.numAtom,
	  scalor1,
	  scalor2,
	  sys.ddata.coord,
	  sys.ddata.velox,
	  sys.ddata.veloy,
	  sys.ddata.veloz);
  checkCUDAError ("NoseHoover_Chains2_Isobaric::operator_L_r");
  if (timer != NULL) timer->tic(mdTimeIntegrator);
}


void NoseHoover_Chains2_Isobaric::
operator_L_xi (const ScalorType & dt)
{
  xi1 += vxi1 * dt;
  // xi2 += vxi2 * dt;
}

static void __global__
NHC_P_L_Cv (ScalorType a,
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

void NoseHoover_Chains2_Isobaric::
operator_L_Cv (const ScalorType & dt,
	       MDSystem & sys,
	       MDStatistic & output_st)
{
  // ScalorType a = expf (-(vxi1 + (1. + 1./NN) * vep) * dt);
  ScalorType a = expf (-(vxi1 + (1. + 3./Nf) * vep) * dt);
  NHC_P_L_Cv
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
  checkCUDAError ("NoseHoover_Chains2_Isobaric::operator_L_Cv NHC_P_L_Cv");
  sum_kxx.sumBuffAdd (output_st.ddata, mdStatisticKineticEnergyXX, 0);
  sum_kyy.sumBuffAdd (output_st.ddata, mdStatisticKineticEnergyYY, 0);
  sum_kzz.sumBuffAdd (output_st.ddata, mdStatisticKineticEnergyZZ, 0);
  checkCUDAError ("NoseHoover_Chains2_Isobaric::operator_L_Cv sum");
}

void NoseHoover_Chains2_Isobaric::
operator_L_Gep (const ScalorType & dt,
		const MDSystem & sys,
		const MDStatistic & input_st)
{
  ScalorType volume = sys.box.size.x * sys.box.size.y * sys.box.size.z;

  // ScalorType Gep = (2. / NN * input_st.kineticEnergy() +
  // 		    3. * volume * (input_st.pressure (sys.ddata.numAtom, sys.box) -
  // 				   refP)
  //     ) / WW;
  ScalorType Gep = (6. / Nf * input_st.kineticEnergy() +
  		    3. * volume * (input_st.pressure (sys.box) -
  				   refP)
      ) / WW;

  // ScalorType Gep1 = ((1. + 1./NN) * 2. * input_st.kineticEnergy() +
  // 		     - input_st.virial() -
  // 		     3. * volume * refP) / WW;
  // ScalorType pressure = input_st.pressure (sys.ddata.numAtom, sys.box);
  // printf ("%f %f %f \n", Gep, Gep1, input_st.pressure (sys.box));
  
  vep += Gep * dt;
}


void NoseHoover_Chains2_Isobaric::
operator_L_vep (const ScalorType & dt)
{
  vep *= expf (-vxi1 * dt);
}

void NoseHoover_Chains2_Isobaric::	
operator_L_G1 (const ScalorType & dt,
	       const MDSystem & sys,
	       const MDStatistic & input_st)
{
  ScalorType G1;
  ScalorType hK = input_st.kineticEnergy();
  
  // G1 = (2.f * hK + WW * vep * vep - (LL + 1.) * refT) / Q1;
  // G1 = (2.f * hK - (LL + 0.) * refT) / Q1;
  G1 = (2.f * hK + WW * vep * vep - (Nf + 1.) * refT) / Q1;

  vxi1 += G1 * dt;
}

// void NoseHoover_Chains2_Isobaric::
// operator_L_vxi1 (const ScalorType & dt)
// {
//   vxi1 *= expf (-vxi2 * dt);
// }

// void NoseHoover_Chains2_Isobaric::
// operator_L_G2 (const ScalorType & dt)
// {
//   ScalorType G2 = (Q1 * vxi1 * vxi1 - refT) / Q2;
//   vxi2 += G2 * dt;
// }


static void __global__
NHC_P_L_V (const ScalorType dt,
	   const IndexType numAtom,
	   const ScalorType * forcx,
	   const ScalorType * forcy,
	   const ScalorType * forcz,
	   const ScalorType * mass,
	   ScalorType * velox,
	   ScalorType * veloy,
	   ScalorType * veloz,
	   ScalorType * st_buffx,
	   ScalorType * st_buffy,
	   ScalorType * st_buffz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  ScalorType vx, vy, vz;
  
  if (ii < numAtom) {
    ScalorType dtmi = dt /mass[ii];
    vx = (velox[ii] += dtmi * forcx[ii]);
    vy = (veloy[ii] += dtmi * forcy[ii]);
    vz = (veloz[ii] += dtmi * forcz[ii]);
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
NHC_P_L_V (const ScalorType dt,
	   const IndexType numAtom,
	   const ScalorType * forcx,
	   const ScalorType * forcy,
	   const ScalorType * forcz,
	   const ScalorType * mass,
	   ScalorType * velox,
	   ScalorType * veloy,
	   ScalorType * veloz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  
  if (ii < numAtom) {
    ScalorType dtmi = dt /mass[ii];
    velox[ii] += dtmi * forcx[ii];
    veloy[ii] += dtmi * forcy[ii];
    veloz[ii] += dtmi * forcz[ii];
  }
}

  

void NoseHoover_Chains2_Isobaric::
operator_L_v (const ScalorType & dt,
	      MDSystem & sys,
	      MDStatistic & output_st,
	      MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  NHC_P_L_V
      <<<atomGridDim, myBlockDim, sharedBuffSize>>> (
	  dt,
	  sys.ddata.numAtom,
	  sys.ddata.forcx,
	  sys.ddata.forcy,
	  sys.ddata.forcz,
	  sys.ddata.mass,
	  sys.ddata.velox,
	  sys.ddata.veloy,
	  sys.ddata.veloz,
	  sum_kxx.buff,
	  sum_kyy.buff,
	  sum_kzz.buff);
  sum_kxx.sumBuffAdd (output_st.ddata, mdStatisticKineticEnergyXX, 0);
  sum_kyy.sumBuffAdd (output_st.ddata, mdStatisticKineticEnergyYY, 0);
  sum_kzz.sumBuffAdd (output_st.ddata, mdStatisticKineticEnergyZZ, 0);
  checkCUDAError ("NoseHoover_Chains2_Isobaric::operator_L_v");
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}

void NoseHoover_Chains2_Isobaric::
operator_L_v (const ScalorType & dt,
	      MDSystem & sys,
	      MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  NHC_P_L_V
      <<<atomGridDim, myBlockDim, sharedBuffSize>>> (
	  dt,
	  sys.ddata.numAtom,
	  sys.ddata.forcx,
	  sys.ddata.forcy,
	  sys.ddata.forcz,
	  sys.ddata.mass,
	  sys.ddata.velox,
	  sys.ddata.veloy,
	  sys.ddata.veloz);
  checkCUDAError ("NoseHoover_Chains2_Isobaric::operator_L_v");
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}


void NoseHoover_Chains2_Isobaric::
operator_L_box (const ScalorType & dt,
		RectangularBox & box,
		MDTimer * timer)
{
  ScalorType tmp = exp (vep * dt);
  ScalorType tmpx = box.size.x * tmp;
  ScalorType tmpy = box.size.y * tmp;
  ScalorType tmpz = box.size.z * tmp;

  ep += dt * vep;
  
  setBoxSize (tmpx, tmpy, tmpz, box);
}



void NoseHoover_Chains2_Isobaric::
operator_L_CP (const ScalorType & dt,
	       MDSystem & sys,
	       const MDStatistic & input_st,
	       MDStatistic & output_st,
	       MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  // operator_L_G2 (0.5 * dt);

  input_st.updateHost();
  // operator_L_vxi1 (0.25 * dt);
  operator_L_G1   (0.5 *  dt, sys, input_st);
  // operator_L_vxi1 (0.25 * dt);

  operator_L_vep  (0.25 * dt);
  operator_L_Gep  (0.5  * dt, sys, input_st);
  operator_L_vep  (0.25 * dt);

  tmp_st.clearDevice();
  operator_L_Cv  (dt, sys, tmp_st);
  tmp_st.add (input_st, 3, virial_array);
  operator_L_xi  (dt);

  operator_L_vep  (0.25 * dt);
  operator_L_Gep  (0.5  * dt, sys, tmp_st);
  operator_L_vep  (0.25 * dt);

  // operator_L_vxi1 (0.25 * dt);
  operator_L_G1   (0.5 *  dt, sys, tmp_st);
  // operator_L_vxi1 (0.25 * dt);

  output_st.add (tmp_st, 3, kinetic_array);
  
  // operator_L_G2 (0.5 * dt);
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}


void NoseHoover_Chains2_Isobaric::
operator_L_CP (const ScalorType & dt,
	       MDSystem & sys,
	       const MDStatistic & input_st,
	       MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeIntegrator);
  // operator_L_G2 (0.5 * dt);

  input_st.updateHost();
  // operator_L_vxi1 (0.25 * dt);
  operator_L_G1   (0.5 *  dt, sys, input_st);
  // operator_L_vxi1 (0.25 * dt);

  operator_L_vep  (0.25 * dt);
  operator_L_Gep  (0.5  * dt, sys, input_st);
  operator_L_vep  (0.25 * dt);

  tmp_st.clearDevice();
  operator_L_Cv  (dt, sys, tmp_st);
  tmp_st.add (input_st, 3, virial_array);
  operator_L_xi  (dt);

  operator_L_vep  (0.25 * dt);
  operator_L_Gep  (0.5  * dt, sys, tmp_st);
  operator_L_vep  (0.25 * dt);

  // operator_L_vxi1 (0.25 * dt);
  operator_L_G1   (0.5 *  dt, sys, tmp_st);
  // operator_L_vxi1 (0.25 * dt);

  // operator_L_G2 (0.5 * dt);
  if (timer != NULL) timer->toc(mdTimeIntegrator);
}


ScalorType NoseHoover_Chains2_Isobaric::
HamiltonianContribution (const RectangularBox & box) const
{
  ScalorType volume = box.size.x * box.size.y * box.size.z;
  return (0.5 * WW * vep * vep +
	  0.5 * Q1 * vxi1 * vxi1 +
	  // 0.5 * Q2 * vxi2 * vxi2 +
	  // (3.*NN + 1.) * refT * xi1 +
	  (Nf + 1.) * refT * xi1 +
	  // refT * xi2 +
	  refP * volume);
}




