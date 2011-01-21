#define DEVICE_CODE

#include "Parallel_Integrator.h"
#include "Parallel_Interface.h"
#include "Parallel_Auxiliary.h"
#include "compile_error_mixcode.h"

Parallel::TranslationalFreedomRemover::
~TranslationalFreedomRemover ()
{
  clear();
}

void Parallel::TranslationalFreedomRemover::
clear ()
{
  if (malloced){
    cudaFree (sums);
    cudaFree (sumM);
    free (hsums);
    malloced = false;
  }
}

void Parallel::TranslationalFreedomRemover::
reinit (const DeviceCellListedMDData & data)
{
  clear ();
  
  IndexType totalNumCell = data.getNumCell().x * data.getNumCell().y * data.getNumCell().z;
  gridDim = toGridDim (totalNumCell);
  numThreadsInCell = Parallel::Interface::numThreadsInCell();
  sharedBuffSize = numThreadsInCell * sizeof(ScalorType);

  sum_x.reinit (totalNumCell, NThreadForSum);
  sum_y.reinit (totalNumCell, NThreadForSum);
  sum_z.reinit (totalNumCell, NThreadForSum);

  if (! malloced){
    cudaMalloc ((void**)&sums, 3 * sizeof(ScalorType));
    Parallel::Auxiliary::setValue <<<1, 3>>> (sums, 3, ScalorType (0.f));
    cudaMalloc ((void**)&sumM, 1 * sizeof(ScalorType));
    checkCUDAError ("TranslationalFreedomRemover::reinit, malloc sums");
    hsums = (ScalorType *) malloc (3 * sizeof(ScalorType));
    if (hsums == NULL) {
      throw MDExcptFailedMallocOnHost ("TranslationalFreedomRemover::reinit",
				       "hsums", 3 * sizeof(ScalorType));
    }				       
    malloced = true;
  }

  SumVector<ScalorType > sum_mass;
  sum_mass.reinit (totalNumCell, NThreadForSum);
  Parallel::CudaGlobal::prepareCalTotalMass
      <<<gridDim, numThreadsInCell, sharedBuffSize>>>(
	  data.dptr_numAtomInCell(),
	  data.dptr_mass(),
	  sum_mass.buff);
  sum_mass.sumBuff (sumM, 0);
  cudaMemcpy (&totalMassi, sumM, sizeof(ScalorType), cudaMemcpyDeviceToHost);
  checkCUDAError ("TranslationalFreedomRemover::reinit, cpy sumM");
  Parallel::HostSystemMass sysMass;
  sysMass.setMass (totalMassi);
  sysMass.sumAll ();
  totalMassi = sysMass.getMass();
  totalMassi = 1.f/totalMassi;
}

void Parallel::TranslationalFreedomRemover::
remove (DeviceCellListedMDData & data)
{
  Parallel::CudaGlobal::prepareRemoveTranslationalFreedom
      <<<gridDim, numThreadsInCell, sharedBuffSize>>> (
	  data.dptr_numAtomInCell(),
	  data.dptr_mass(),
	  data.dptr_velocityX(),
	  data.dptr_velocityY(),
	  data.dptr_velocityZ(),
	  sum_x.buff,
	  sum_y.buff,
	  sum_z.buff);
  checkCUDAError ("TranslationalFreedomRemover::remove, prepare");
  sum_x.sumBuff (sums, 0);
  sum_y.sumBuff (sums, 1);
  sum_z.sumBuff (sums, 2);
  cudaMemcpy (hsums, sums, 3 * sizeof(ScalorType), cudaMemcpyDeviceToHost);
  checkCUDAError ("TranslationalFreedomRemover::remove, cpy p to host");
  hmomentum.setMomentunX (hsums[0]);
  hmomentum.setMomentunY (hsums[1]);
  hmomentum.setMomentunZ (hsums[2]);
  hmomentum.sumAll ();
  hsums[0] = hmomentum.getMomentumX();
  hsums[1] = hmomentum.getMomentumY();
  hsums[2] = hmomentum.getMomentumZ();
  cudaMemcpy (sums, hsums, 3 * sizeof(ScalorType), cudaMemcpyHostToDevice);
  checkCUDAError ("TranslationalFreedomRemover::remove, cpy p to device");
  Parallel::CudaGlobal::removeTranslationalFreedom
      <<<gridDim, numThreadsInCell>>> (
	  data.dptr_numAtomInCell(),
	  totalMassi,
	  sums,
	  data.dptr_velocityX(),
	  data.dptr_velocityY(),
	  data.dptr_velocityZ());
  checkCUDAError("TranslationalFreedomRemover::remove, remove");
}


__global__ void Parallel::CudaGlobal::
prepareCalTotalMass (const IndexType * numAtomInCell,
		     const ScalorType * mass,
		     ScalorType * mass_buff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  IndexType this_numAtomInCell = numAtomInCell[bid];
  if (this_numAtomInCell == 0){
    if (threadIdx.x == 0){
      mass_buff[bid] = 0;
    }
    return;
  }
  extern __shared__ ScalorType buff[];
  if (threadIdx.x < this_numAtomInCell){
    buff[threadIdx.x] = mass[ii];
  }
  else {
    buff[threadIdx.x] = 0.f;
  }
  __syncthreads ();
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) mass_buff[bid] = buff[0];
}

__global__ void Parallel::CudaGlobal::
prepareRemoveTranslationalFreedom (const IndexType * numAtomInCell,
				   const ScalorType * mass,
				   const ScalorType * velox,
				   const ScalorType * veloy,
				   const ScalorType * veloz,
				   ScalorType * st_buff_x,
				   ScalorType * st_buff_y,
				   ScalorType * st_buff_z)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;  
  IndexType this_numAtomInCell = numAtomInCell[bid];
  
  if (this_numAtomInCell == 0){
    if (threadIdx.x == 0){
      st_buff_x[bid] = 0;
      st_buff_y[bid] = 0;
      st_buff_z[bid] = 0;
    }
    return;
  }

  extern __shared__  ScalorType buff[];
  if (threadIdx.x < this_numAtomInCell){
    buff[threadIdx.x] = mass[ii] * velox[ii];
  }
  else {
    buff[threadIdx.x] = 0.f;
  }
  __syncthreads();
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) st_buff_x[bid] = buff[0];
  __syncthreads();

  if (threadIdx.x < this_numAtomInCell){
    buff[threadIdx.x] = mass[ii] * veloy[ii];
  }
  __syncthreads();
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) st_buff_y[bid] = buff[0];
  __syncthreads();

  if (threadIdx.x < this_numAtomInCell){
    buff[threadIdx.x] = mass[ii] * veloz[ii];
  }
  __syncthreads();
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) st_buff_z[bid] = buff[0];
}

__global__ void Parallel::CudaGlobal::
removeTranslationalFreedom (const IndexType * numAtomInCell,
			    const ScalorType totalMassi,
			    const ScalorType * sums,
			    ScalorType * velox,
			    ScalorType * veloy,
			    ScalorType * veloz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;  
  IndexType this_numAtomInCell = numAtomInCell[bid];
  
  if (this_numAtomInCell == 0){
    return;
  }

  __shared__ ScalorType buffSums[3];
  if (threadIdx.x < 3){
    buffSums[threadIdx.x] = sums[threadIdx.x];
  }
  __syncthreads();

  if (threadIdx.x < this_numAtomInCell){
    velox[ii] -= buffSums[0] * totalMassi;
    veloy[ii] -= buffSums[1] * totalMassi;
    veloz[ii] -= buffSums[2] * totalMassi;
  }
}


void Parallel::Integrator::VelocityVerlet::
reinit (const DeviceCellListedMDData & ddata)
{
  IndexType totalNumCell = ddata.getNumCell().x *
      ddata.getNumCell().y * ddata.getNumCell().z;
  gridDim = toGridDim (totalNumCell);

  sum_kxx.reinit (totalNumCell, NThreadForSum);
  sum_kyy.reinit (totalNumCell, NThreadForSum);
  sum_kzz.reinit (totalNumCell, NThreadForSum);
  
  sharedBuffSize = Parallel::Interface::numThreadsInCell() * sizeof(ScalorType);
}

void Parallel::Integrator::VelocityVerlet::
step1 (DeviceCellListedMDData & ddata,
       const ScalorType & dt)
{
  Parallel::CudaGlobal::velocityVerlet_step1
      <<<gridDim, Parallel::Interface::numThreadsInCell()>>>(
	  ddata.dptr_numAtomInCell(),
	  ddata.dptr_forceX(),
	  ddata.dptr_forceY(),
	  ddata.dptr_forceZ(),
	  ddata.dptr_mass(),
	  dt,
	  ddata.dptr_coordinate(),
	  ddata.dptr_velocityX(),
	  ddata.dptr_velocityY(),
	  ddata.dptr_velocityZ());
  checkCUDAError ("Integrator::VelocityVerlet::step1");	  
}

void Parallel::Integrator::VelocityVerlet::
step2 (DeviceCellListedMDData & data,
       const ScalorType & dt)
{
  Parallel::CudaGlobal::velocityVerlet_step2
      <<<gridDim, Parallel::Interface::numThreadsInCell()>>>(
	  data.dptr_numAtomInCell(),
	  data.dptr_forceX(),
	  data.dptr_forceY(),
	  data.dptr_forceZ(),
	  data.dptr_mass(),
	  dt,
	  data.dptr_velocityX(),
	  data.dptr_velocityY(),
	  data.dptr_velocityZ());
  checkCUDAError ("Integrator::VelocityVerlet::step2, no st");
}

void Parallel::Integrator::VelocityVerlet::
step2 (DeviceCellListedMDData & data,
       const ScalorType & dt,
       DeviceStatistic & st)
{
  Parallel::CudaGlobal::velocityVerlet_step2
      <<<gridDim, Parallel::Interface::numThreadsInCell(), sharedBuffSize>>> (
	  data.dptr_numAtomInCell(),
	  data.dptr_forceX(),
	  data.dptr_forceY(),
	  data.dptr_forceZ(),
	  data.dptr_mass(),
	  dt,
	  data.dptr_velocityX(),
	  data.dptr_velocityY(),
	  data.dptr_velocityZ(),
	  sum_kxx.buff,
	  sum_kyy.buff,
	  sum_kzz.buff);
  sum_kxx.sumBuffAdd (st.dptr_statisticData(), mdStatistic_KineticEnergyXX, 0);
  sum_kyy.sumBuffAdd (st.dptr_statisticData(), mdStatistic_KineticEnergyYY, 0);
  sum_kzz.sumBuffAdd (st.dptr_statisticData(), mdStatistic_KineticEnergyZZ, 0);

  checkCUDAError ("Integrator::VelocityVerlet::step2, with st");
}
	  

	  

__global__ void Parallel::CudaGlobal::
velocityVerlet_step1 (const IndexType * numAtomInCell,
		      const ScalorType * forcx,
		      const ScalorType * forcy,
		      const ScalorType * forcz,
		      const ScalorType * mass,
		      const ScalorType   dt,
		      CoordType * coord,
		      ScalorType * velox,
		      ScalorType * veloy,
		      ScalorType * veloz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  if (tid < numAtomInCell[bid]){
    // // velox[ii] = veloy[ii] = veloz[ii] = 1;
    // velox[ii] = 1;
    // veloy[ii] = veloz[ii] = 0;
    ScalorType hdtmi = 0.5 * dt / mass[ii];
    velox[ii]   += hdtmi * forcx[ii];
    coord[ii].x += dt * velox[ii];
    veloy[ii]   += hdtmi * forcy[ii];
    coord[ii].y += dt * veloy[ii];
    veloz[ii]   += hdtmi * forcz[ii];
    coord[ii].z += dt * veloz[ii];
  }
}

__global__ void Parallel::CudaGlobal::
velocityVerlet_step2 (const IndexType * numAtomInCell,
		      const ScalorType * forcx,
		      const ScalorType * forcy,
		      const ScalorType * forcz,
		      const ScalorType * mass,
		      const ScalorType   dt,
		      ScalorType * velox,
		      ScalorType * veloy,
		      ScalorType * veloz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  if (tid < numAtomInCell[bid]){
    ScalorType hdtmi = 0.5f * dt / mass[ii];
    velox[ii] += hdtmi * forcx[ii];
    veloy[ii] += hdtmi * forcy[ii];
    veloz[ii] += hdtmi * forcz[ii];
  }
}

__global__ void Parallel::CudaGlobal::
velocityVerlet_step2 (const IndexType * numAtomInCell,
		      const ScalorType * forcx,
		      const ScalorType * forcy,
		      const ScalorType * forcz,
		      const ScalorType * mass,
		      const ScalorType   dt,
		      ScalorType * velox,
		      ScalorType * veloy,
		      ScalorType * veloz,
		      ScalorType * statistic_buffxx,
		      ScalorType * statistic_buffyy,
		      ScalorType * statistic_buffzz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  
  extern __shared__  ScalorType buff [];
  ScalorType vx(0.f), vy(0.f), vz(0.f);
  IndexType this_numAtomInCell = numAtomInCell[bid];
  
  if (threadIdx.x < this_numAtomInCell) {
    ScalorType hdtmi = 0.5f * dt / mass[ii];
    vx = (velox[ii] += hdtmi * forcx[ii]);
    vy = (veloy[ii] += hdtmi * forcy[ii]);
    vz = (veloz[ii] += hdtmi * forcz[ii]);
  }

  ScalorType scalor = 0.5f * mass[ii];
  // if (threadIdx.x < this_numAtomInCell) scalor = 0.5f * mass[ii];
  // else scalor = 0.f;
  
  buff[threadIdx.x] = scalor * vx * vx;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffxx[bid] = buff[0];
  __syncthreads();
  buff[threadIdx.x] = scalor * vy * vy;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffyy[bid] = buff[0];
  __syncthreads();
  buff[threadIdx.x] = scalor * vz * vz;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffzz[bid] = buff[0];
}

  

void Parallel::Integrator::LeapFrog::
reinit (const DeviceCellListedMDData & ddata)
{
  IndexType totalNumCell = ddata.getNumCell().x *
      ddata.getNumCell().y * ddata.getNumCell().z;
  gridDim = toGridDim (totalNumCell);

  sum_kxx.reinit (totalNumCell, NThreadForSum);
  sum_kyy.reinit (totalNumCell, NThreadForSum);
  sum_kzz.reinit (totalNumCell, NThreadForSum);
  
  sharedBuffSize = Parallel::Interface::numThreadsInCell() * sizeof(ScalorType);
}
  
void Parallel::Integrator::LeapFrog::
stepX (DeviceCellListedMDData & ddata,
       const ScalorType & dt)
{
  Parallel::CudaGlobal::leapFrogStepX
      <<<gridDim, Parallel::Interface::numThreadsInCell()>>> (
	  ddata.dptr_numAtomInCell(),
	  ddata.dptr_velocityX(),
	  ddata.dptr_velocityY(),
	  ddata.dptr_velocityZ(),
	  dt,
	  ddata.dptr_coordinate());
  checkCUDAError ("Integrator::LeapFrog::stepX");	  
}

void Parallel::Integrator::LeapFrog::
stepV (DeviceCellListedMDData & ddata,
       const ScalorType & dt)
{
  Parallel::CudaGlobal::leapFrogStepV
      <<<gridDim, Parallel::Interface::numThreadsInCell()>>>(
	  ddata.dptr_numAtomInCell(),
	  ddata.dptr_forceX(),
	  ddata.dptr_forceY(),
	  ddata.dptr_forceZ(),
	  ddata.dptr_mass(),
	  dt,
	  ddata.dptr_velocityX(),
	  ddata.dptr_velocityY(),
	  ddata.dptr_velocityZ());
  checkCUDAError ("Integrator::LeapFrog::stepV, without st");
}

void Parallel::Integrator::LeapFrog::
stepV (DeviceCellListedMDData & ddata,
       const ScalorType & dt,
       DeviceStatistic & st)
{
  Parallel::CudaGlobal::leapFrogStepV
      <<<gridDim, Parallel::Interface::numThreadsInCell(), sharedBuffSize>>>(
	  ddata.dptr_numAtomInCell(),
	  ddata.dptr_forceX(),
	  ddata.dptr_forceY(),
	  ddata.dptr_forceZ(),
	  ddata.dptr_mass(),
	  dt,
	  ddata.dptr_velocityX(),
	  ddata.dptr_velocityY(),
	  ddata.dptr_velocityZ(),
	  sum_kxx.buff,
	  sum_kyy.buff,
	  sum_kzz.buff);
  sum_kxx.sumBuffAdd (st.dptr_statisticData(), mdStatistic_KineticEnergyXX, 0);
  sum_kyy.sumBuffAdd (st.dptr_statisticData(), mdStatistic_KineticEnergyYY, 0);
  sum_kzz.sumBuffAdd (st.dptr_statisticData(), mdStatistic_KineticEnergyZZ, 0);
  checkCUDAError ("Integrator::LeapFrog::stepV, with st");
}


__global__ void Parallel::CudaGlobal::
leapFrogStepX (const IndexType * numAtomInCell,
	       const ScalorType * velox,
	       const ScalorType * veloy,
	       const ScalorType * veloz,
	       const ScalorType dt,
	       CoordType * coord)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType this_numAtomInCell = numAtomInCell[bid];
  if (this_numAtomInCell == 0){
    return;
  }
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (threadIdx.x < this_numAtomInCell){
    coord[ii].x += dt * velox[ii];
    coord[ii].y += dt * veloy[ii];
    coord[ii].z += dt * veloz[ii];
  }
}

__global__ void Parallel::CudaGlobal::
leapFrogStepV (const IndexType * numAtomInCell,
	       const ScalorType * forcx,
	       const ScalorType * forcy,
	       const ScalorType * forcz,
	       const ScalorType * mass,
	       const ScalorType dt,
	       ScalorType * velox,
	       ScalorType * veloy,
	       ScalorType * veloz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType this_numAtomInCell = numAtomInCell[bid];
  if (this_numAtomInCell == 0){
    return;
  }
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (threadIdx.x < this_numAtomInCell){
    ScalorType dtmi = dt / mass[ii];
    velox[ii] += dtmi * forcx[ii];
    veloy[ii] += dtmi * forcy[ii];
    veloz[ii] += dtmi * forcz[ii];
  }
}


__global__ void Parallel::CudaGlobal::
leapFrogStepV (const IndexType * numAtomInCell,
	       const ScalorType * forcx,
	       const ScalorType * forcy,
	       const ScalorType * forcz,
	       const ScalorType * mass,
	       const ScalorType dt,
	       ScalorType * velox,
	       ScalorType * veloy,
	       ScalorType * veloz,
	       ScalorType * statistic_buffxx,
	       ScalorType * statistic_buffyy,
	       ScalorType * statistic_buffzz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType this_numAtomInCell = numAtomInCell[bid];
  IndexType ii = threadIdx.x + bid * blockDim.x;
  extern __shared__ ScalorType buff [];
  ScalorType vx(0.f), vy(0.f), vz(0.f);
  ScalorType scalor = 0.f;
  
  if (threadIdx.x < this_numAtomInCell){
    ScalorType dtmi = dt / mass[ii];
    scalor = 0.5f * mass[ii];
    vx = (velox[ii] += dtmi * forcx[ii]);
    vy = (veloy[ii] += dtmi * forcy[ii]);
    vz = (veloz[ii] += dtmi * forcz[ii]);
  }

  buff[threadIdx.x] = scalor * vx * vx;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffxx[bid] = buff[0];
  __syncthreads();
  buff[threadIdx.x] = scalor * vy * vy;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffyy[bid] = buff[0];
  __syncthreads();
  buff[threadIdx.x] = scalor * vz * vz;
  sumVectorBlockBuffer_2 (buff);
  if (threadIdx.x == 0) statistic_buffzz[bid] = buff[0];
}


using namespace RectangularBoxGeometry;


Parallel::Integrator::BerendsenLeapFrog::
BerendsenLeapFrog ()
{
  TCoupleOn = false;
  PCoupleOn = false;
  NPCoupleGroup = 0;
  ptr_inter = NULL;
  ptr_bdInterList = NULL;
  nstep = 0;
}


void Parallel::Integrator::BerendsenLeapFrog::
reinit (const MDSystem & sys,
	const ScalorType & dt_,
	InteractionEngine * ptr_inter_,
	TranslationalFreedomRemover * ptr_trRemover_,
	const IndexType removeFeq_,
	DeviceBondList * ptr_bdInterList_) 
{
  dt = dt_;
  TCoupleOn = false;
  PCoupleOn = false;
  NPCoupleGroup = 0;
  nstep = 0;

  IntVectorType numCell (sys.deviceData.getNumCell());
  IndexType nob = numCell.x * numCell.y * numCell.z;
  gridDim = toGridDim (nob);  

  lpfrog.reinit (sys.deviceData);
  ptr_inter = ptr_inter_;
  ptr_trRemover = ptr_trRemover_;
  removeFeq = removeFeq_;
  ptr_bdInterList = ptr_bdInterList_;

  relation.rebuild (sys.deviceData);
  if (ptr_bdInterList != NULL){
    Parallel::SubCellList ghost, innerShell;
    sys.deviceData.buildSubListGhostCell (ghost);
    sys.deviceData.buildSubListInnerShell (innerShell);
    ghost.add (innerShell);
    relation_buildBdList.rebuild (sys.deviceData, innerShell, ghost);
  }
}


void Parallel::Integrator::BerendsenLeapFrog::
TCouple (const ScalorType & refT_,
	 const ScalorType & tauT_)
{
  TCoupleOn = true;
  refT = refT_;
  tauT = tauT_;
}

void Parallel::Integrator::BerendsenLeapFrog::
addPcoupleGroup (const BoxDirection_t & direction,
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

using namespace Parallel::Timer;

void Parallel::Integrator::BerendsenLeapFrog::
firstStep (MDSystem & sys,
	   DeviceStatistic &st)
{
  myst.clearData ();
  DeviceTimer::tic (item_ClearInteraction);
  ptr_inter->clearInteraction (sys.deviceData);
  DeviceTimer::toc (item_ClearInteraction);
  HostTimer::tic (item_TransferGhost);
  sys.transferGhost();  
  HostTimer::toc (item_TransferGhost);
  
  DeviceTimer::tic (item_NonBondedInterStatistic);
  ptr_inter->applyNonBondedInteraction (sys.deviceData, relation);
  DeviceTimer::toc (item_NonBondedInterStatistic);
  if (ptr_bdInterList != NULL){
    DeviceTimer::tic (item_BuildBondList);
    buildDeviceBondList (sys.deviceData, relation_buildBdList, *ptr_bdInterList);
    DeviceTimer::toc (item_BuildBondList);
    DeviceTimer::tic (item_BondedInterStatistic);
    ptr_inter->applyBondedInteraction (sys.deviceData, *ptr_bdInterList);
    DeviceTimer::toc (item_BondedInterStatistic);
  }
  
  HostTimer::tic (item_TransferGhost);
  sys.clearGhost();
  DeviceTimer::toc (item_ClearInteraction);

  DeviceTimer::tic (item_Integrate);
  lpfrog.stepX (sys.deviceData, dt);
  lpfrog.stepV (sys.deviceData, dt, myst);
  DeviceTimer::toc (item_Integrate);
  
  DeviceTimer::tic (item_BuildCellList);
  if (ptr_bdInterList == NULL){
    sys.deviceData.rebuild ();
  }
  else{
    sys.deviceData.rebuild (*ptr_bdInterList);
  }
  DeviceTimer::toc (item_BuildCellList);
  HostTimer::tic (item_Redistribute);
  sys.redistribute ();
  HostTimer::toc (item_Redistribute);

  DeviceTimer::tic (item_ClearInteraction);
  ptr_inter->clearInteraction (sys.deviceData);
  DeviceTimer::toc (item_ClearInteraction);
  HostTimer::tic (item_TransferGhost);
  sys.transferGhost();  
  HostTimer::toc (item_TransferGhost);
  
  DeviceTimer::tic (item_NonBondedInterStatistic);
  ptr_inter->applyNonBondedInteraction (sys.deviceData, relation, myst);
  DeviceTimer::toc (item_NonBondedInterStatistic);
  if (ptr_bdInterList != NULL){
    DeviceTimer::tic (item_BuildBondList);
    buildDeviceBondList (sys.deviceData, relation_buildBdList, *ptr_bdInterList);
    DeviceTimer::toc (item_BuildBondList);
    DeviceTimer::tic (item_BondedInterStatistic);
    ptr_inter->applyBondedInteraction (sys.deviceData, *ptr_bdInterList, myst);
    DeviceTimer::toc (item_BondedInterStatistic);
  }
  
  HostTimer::tic (item_TransferGhost);
  sys.clearGhost();
  HostTimer::toc (item_TransferGhost);
  st.add (myst);
  nstep ++;
}



void Parallel::Integrator::BerendsenLeapFrog::
oneStep (MDSystem & sys,
	 DeviceStatistic &st)
{
  ScalorType nowT, lambda;
  ScalorType nowP[3], mu[3];
  IndexType nDir[3];

  DeviceTimer::tic (item_RemoveTransFreedom);
  if (nstep % removeFeq == 0 && ptr_trRemover != NULL){
    ptr_trRemover->remove (sys.deviceData);
  }
  DeviceTimer::toc (item_RemoveTransFreedom);
  
  if (nstep != 0) {
    DeviceTimer::tic (item_Integrate);
    myst.copyToHost (myhst);
    myhst.collectDataAll ();
    if (TCoupleOn){
      nowT = myhst.kineticEnergy();
      nowT *= 2.f / (sys.getNumFreedom() - 3);
      lambda = sqrtf(1.f + dt / tauT * (refT / nowT - 1.f));
    }
    if (PCoupleOn){
      for (IndexType i = 0; i < NPCoupleGroup; ++i){
	nowP[i] = 0;
	nDir[i] = 0;
	if ((PCoupleDirections[i] & mdRectBoxDirectionX) != 0){
	  nowP[i] += myhst.pressureXX(sys.deviceData.getGlobalBox());
	  nDir[i] ++;
	}
	if ((PCoupleDirections[i] & mdRectBoxDirectionY) != 0){
	  nowP[i] += myhst.pressureYY(sys.deviceData.getGlobalBox());
	  nDir[i] ++;
	}
	if ((PCoupleDirections[i] & mdRectBoxDirectionZ) != 0){
	  nowP[i] += myhst.pressureZZ(sys.deviceData.getGlobalBox());
	  nDir[i] ++;
	}
	nowP[i] /= ScalorType(nDir[i]);
	mu [i] = powf (1.f + dt / tauP[i] * betaP[i] * (nowP[i] - refP[i]), 1.f/3.f);
      }
    }
  
    myst.clearData();
    lpfrog.stepV (sys.deviceData, dt, myst);
    if (TCoupleOn){
      HostVectorType scale;
      scale.x = (scale.y = (scale.z = lambda));
      sys.deviceData.rescaleVelocity (scale);
      myst.rescale (mdStatistic_KineticEnergyXX, 3, lambda*lambda);
    }
    // st.add (myst);
    // myst.clearData();
    lpfrog.stepX (sys.deviceData, dt);
    HostVectorType coordScalor ;
    if (PCoupleOn){
      coordScalor.x = (coordScalor.y = (coordScalor.z = 1.f));
      for (IndexType i = 0; i < NPCoupleGroup; ++i){
	if ((PCoupleDirections[i] & mdRectBoxDirectionX) != 0){
	  coordScalor.x *= mu[i];
	}
	if ((PCoupleDirections[i] & mdRectBoxDirectionY) != 0){
	  coordScalor.y *= mu[i];
	}
	if ((PCoupleDirections[i] & mdRectBoxDirectionZ) != 0){
	  coordScalor.z *= mu[i];
	}
      }
    }
    nstep ++;
    DeviceTimer::tic (item_Integrate);
    DeviceTimer::tic (item_BuildCellList);
    if (ptr_bdInterList == NULL){
      sys.deviceData.rebuild ();
    }
    else {
      sys.deviceData.rebuild (*ptr_bdInterList);
    }
    DeviceTimer::toc (item_BuildCellList);
    HostTimer::tic (item_Redistribute);
    sys.redistribute();
    HostTimer::toc (item_Redistribute);
  
    if (PCoupleOn){
      DeviceTimer::tic (item_Integrate);
      sys.deviceData.rescaleCoordinate (coordScalor);
      DeviceTimer::toc (item_Integrate);
      DeviceTimer::tic (item_BuildCellList);
      // bool reinited = sys.reinitCellStructure (3.334);
      bool reinited = sys.reinitCellStructure (ptr_inter->getMaxRcut() + 0.1f);
      DeviceTimer::toc (item_BuildCellList);
      if (reinited){
	printf ("# change too much, reinit \n");
	fflush (stdout);
	HostTimer::tic (item_Redistribute);
	sys.redistribute();
	HostTimer::toc (item_Redistribute);
	DeviceTimer::tic (item_NonBondedInteraction);
	ptr_inter->reinit (sys.deviceData);
	DeviceTimer::toc (item_NonBondedInteraction);
	if (ptr_trRemover != NULL){
	  DeviceTimer::tic (item_RemoveTransFreedom);
	  ptr_trRemover->reinit (sys.deviceData);
	  DeviceTimer::toc (item_RemoveTransFreedom);
	}
	if (ptr_bdInterList != NULL){
	  DeviceTimer::tic (item_BondedInteraction);
	  ptr_bdInterList->reinit (sys.deviceData);
	  DeviceTimer::toc (item_BondedInteraction);
	}
	IntVectorType numCell (sys.deviceData.getNumCell());
	gridDim = toGridDim (numCell.x * numCell.y * numCell.z);  
	DeviceTimer::tic (item_Integrate);
	lpfrog.reinit (sys.deviceData);
	DeviceTimer::toc (item_Integrate);
	relation.rebuild (sys.deviceData);
	if (ptr_bdInterList != NULL){
	  buildDeviceBondList (sys.deviceData, relation, *ptr_bdInterList);
	  Parallel::SubCellList ghost, innerShell;
	  sys.deviceData.buildSubListGhostCell (ghost);
	  sys.deviceData.buildSubListInnerShell (innerShell);
	  ghost.add (innerShell);
	  relation_buildBdList.rebuild (sys.deviceData, innerShell, ghost);
	}
	DeviceTimer::toc (item_Integrate);
      }
    }    
    
    DeviceTimer::tic (item_ClearInteraction);
    ptr_inter->clearInteraction (sys.deviceData);
    DeviceTimer::toc (item_ClearInteraction);
    HostTimer::tic (item_TransferGhost);
    sys.transferGhost();
    HostTimer::toc (item_TransferGhost);
    
    DeviceTimer::tic (item_NonBondedInterStatistic);
    ptr_inter->applyNonBondedInteraction (sys.deviceData, relation, myst);
    DeviceTimer::toc (item_NonBondedInterStatistic);
    if (ptr_bdInterList != NULL){
      DeviceTimer::tic (item_BuildBondList);
      buildDeviceBondList (sys.deviceData, relation_buildBdList, *ptr_bdInterList);
      DeviceTimer::toc (item_BuildBondList);
      DeviceTimer::tic (item_BondedInterStatistic);
      ptr_inter->applyBondedInteraction (sys.deviceData, *ptr_bdInterList, myst);
      DeviceTimer::toc (item_BondedInterStatistic);
    }
    
    HostTimer::tic (item_TransferGhost);
    sys.clearGhost();
    HostTimer::toc (item_TransferGhost);
    st.add (myst);
  }
  else {
    firstStep (sys, st);
  }
}



void Parallel::Integrator::BerendsenLeapFrog::
firstStep (MDSystem & sys)
{
  myst.clearData ();
  DeviceTimer::tic (item_ClearInteraction);
  ptr_inter->clearInteraction (sys.deviceData);
  DeviceTimer::toc (item_ClearInteraction);
  HostTimer::tic (item_TransferGhost);
  sys.transferGhost();  
  HostTimer::toc (item_TransferGhost);
  DeviceTimer::tic (item_NonBondedInteraction);
  ptr_inter->applyNonBondedInteraction (sys.deviceData, relation);
  DeviceTimer::toc (item_NonBondedInteraction);
  if (ptr_bdInterList != NULL){
    DeviceTimer::tic (item_BuildBondList);
    buildDeviceBondList (sys.deviceData, relation_buildBdList, *ptr_bdInterList);
    DeviceTimer::toc (item_BuildBondList);
    DeviceTimer::tic (item_BondedInteraction);
    ptr_inter->applyBondedInteraction (sys.deviceData, *ptr_bdInterList);
    DeviceTimer::toc (item_BondedInteraction);
  }
  HostTimer::tic (item_TransferGhost);
  sys.clearGhost();
  DeviceTimer::toc (item_ClearInteraction);

  DeviceTimer::tic (item_Integrate);
  lpfrog.stepX (sys.deviceData, dt);
  lpfrog.stepV (sys.deviceData, dt, myst);
  DeviceTimer::toc (item_Integrate);
  
  DeviceTimer::tic (item_BuildCellList);
  if (ptr_bdInterList == NULL){
    sys.deviceData.rebuild ();
  }
  else {
    sys.deviceData.rebuild (*ptr_bdInterList);
  }
  DeviceTimer::toc (item_BuildCellList);
  HostTimer::tic (item_Redistribute);
  sys.redistribute ();
  HostTimer::toc (item_Redistribute);

  DeviceTimer::tic (item_ClearInteraction);
  ptr_inter->clearInteraction (sys.deviceData);
  DeviceTimer::toc (item_ClearInteraction);
  HostTimer::tic (item_TransferGhost);
  sys.transferGhost();  
  HostTimer::toc (item_TransferGhost);
  DeviceTimer::tic (item_NonBondedInterStatistic);
  ptr_inter->applyNonBondedInteraction (sys.deviceData, relation, myst);
  DeviceTimer::toc (item_NonBondedInterStatistic);
  if (ptr_bdInterList != NULL){
    DeviceTimer::tic (item_BuildBondList);
    buildDeviceBondList (sys.deviceData, relation_buildBdList, *ptr_bdInterList);
    DeviceTimer::toc (item_BuildBondList);
    DeviceTimer::tic (item_BondedInterStatistic);
    ptr_inter->applyBondedInteraction (sys.deviceData, *ptr_bdInterList, myst);
    DeviceTimer::toc (item_BondedInterStatistic);
  }
  HostTimer::tic (item_TransferGhost);
  sys.clearGhost();
  HostTimer::toc (item_TransferGhost);
  nstep ++;
}



void Parallel::Integrator::BerendsenLeapFrog::
oneStep (MDSystem & sys)
{
  ScalorType nowT, lambda;
  ScalorType nowP[3], mu[3];
  IndexType nDir[3];

  DeviceTimer::tic (item_RemoveTransFreedom);
  if (nstep % removeFeq == 0 && ptr_trRemover != NULL){
    ptr_trRemover->remove (sys.deviceData);
  }
  DeviceTimer::toc (item_RemoveTransFreedom);
  
  if (nstep != 0) {
    DeviceTimer::tic (item_Integrate);
    myst.copyToHost (myhst);
    myhst.collectDataAll ();
    if (TCoupleOn){
      nowT = myhst.kineticEnergy();
      nowT *= 2.f / (sys.getNumFreedom() - 3);
      lambda = sqrtf(1.f + dt / tauT * (refT / nowT - 1.f));
    }
    if (PCoupleOn){
      for (IndexType i = 0; i < NPCoupleGroup; ++i){
	nowP[i] = 0;
	nDir[i] = 0;
	if ((PCoupleDirections[i] & mdRectBoxDirectionX) != 0){
	  nowP[i] += myhst.pressureXX(sys.deviceData.getGlobalBox());
	  nDir[i] ++;
	}
	if ((PCoupleDirections[i] & mdRectBoxDirectionY) != 0){
	  nowP[i] += myhst.pressureYY(sys.deviceData.getGlobalBox());
	  nDir[i] ++;
	}
	if ((PCoupleDirections[i] & mdRectBoxDirectionZ) != 0){
	  nowP[i] += myhst.pressureZZ(sys.deviceData.getGlobalBox());
	  nDir[i] ++;
	}
	nowP[i] /= ScalorType(nDir[i]);
	mu [i] = powf (1.f + dt / tauP[i] * betaP[i] * (nowP[i] - refP[i]), 1.f/3.f);
      }
    }
  
    myst.clearData();
    if (TCoupleOn){
      lpfrog.stepV (sys.deviceData, dt, myst);
    }
    else {
      lpfrog.stepV (sys.deviceData, dt);
    }
    if (TCoupleOn){
      HostVectorType scale;
      scale.x = (scale.y = (scale.z = lambda));
      sys.deviceData.rescaleVelocity (scale);
      myst.rescale (mdStatistic_KineticEnergyXX, 3, lambda*lambda);
    }
    // st.add (myst);
    // myst.clearData();
    lpfrog.stepX (sys.deviceData, dt);
    HostVectorType coordScalor ;
    if (PCoupleOn){
      coordScalor.x = (coordScalor.y = (coordScalor.z = 1.f));
      for (IndexType i = 0; i < NPCoupleGroup; ++i){
	if ((PCoupleDirections[i] & mdRectBoxDirectionX) != 0){
	  coordScalor.x *= mu[i];
	}
	if ((PCoupleDirections[i] & mdRectBoxDirectionY) != 0){
	  coordScalor.y *= mu[i];
	}
	if ((PCoupleDirections[i] & mdRectBoxDirectionZ) != 0){
	  coordScalor.z *= mu[i];
	}
      }
    }
    nstep ++;
    DeviceTimer::tic (item_Integrate);
    DeviceTimer::tic (item_BuildCellList);
    if (ptr_bdInterList == NULL){
      sys.deviceData.rebuild ();
    }
    else {
      sys.deviceData.rebuild (*ptr_bdInterList);
    }
    DeviceTimer::toc (item_BuildCellList);
    HostTimer::tic (item_Redistribute);
    sys.redistribute();
    HostTimer::toc (item_Redistribute);
    
    if (PCoupleOn){
      DeviceTimer::tic (item_Integrate);
      sys.deviceData.rescaleCoordinate (coordScalor);
      DeviceTimer::toc (item_Integrate);
      DeviceTimer::tic (item_BuildCellList);
      bool reinited = sys.reinitCellStructure (ptr_inter->getMaxRcut() + 0.1f);
      DeviceTimer::toc (item_BuildCellList);
      if (reinited){
	printf ("# change too much, reinit 2 \n");
	HostTimer::tic (item_Redistribute);
	sys.redistribute();
	HostTimer::toc (item_Redistribute);
	DeviceTimer::tic (item_NonBondedInteraction);
	ptr_inter->reinit (sys.deviceData);
	DeviceTimer::toc (item_NonBondedInteraction);
	if (ptr_trRemover != NULL){
	  DeviceTimer::tic (item_RemoveTransFreedom);
	  ptr_trRemover->reinit (sys.deviceData);
	  DeviceTimer::toc (item_RemoveTransFreedom);
	}
	if (ptr_bdInterList != NULL){
	  DeviceTimer::tic (item_BondedInteraction);
	  ptr_bdInterList->reinit (sys.deviceData);
	  DeviceTimer::toc (item_BondedInteraction);
	}
	IntVectorType numCell (sys.deviceData.getNumCell());
	gridDim = toGridDim (numCell.x * numCell.y * numCell.z);  
	DeviceTimer::tic (item_Integrate);
	lpfrog.reinit (sys.deviceData);
	DeviceTimer::toc (item_Integrate);
	relation.rebuild (sys.deviceData);
	if (ptr_bdInterList != NULL){
	  buildDeviceBondList (sys.deviceData, relation, *ptr_bdInterList);
	  Parallel::SubCellList ghost, innerShell;
	  sys.deviceData.buildSubListGhostCell (ghost);
	  sys.deviceData.buildSubListInnerShell (innerShell);
	  ghost.add (innerShell);
	  relation_buildBdList.rebuild (sys.deviceData, innerShell, ghost);
	}
	DeviceTimer::toc (item_Integrate);
      }
    }    
    
    DeviceTimer::tic (item_ClearInteraction);
    ptr_inter->clearInteraction (sys.deviceData);
    DeviceTimer::toc (item_ClearInteraction);
    HostTimer::tic (item_TransferGhost);
    sys.transferGhost();  
    HostTimer::toc (item_TransferGhost);
    if (PCoupleOn){
      DeviceTimer::tic (item_NonBondedInterStatistic);
      ptr_inter->applyNonBondedInteraction (sys.deviceData, relation, myst);
      DeviceTimer::toc (item_NonBondedInterStatistic);
    }
    else {
      DeviceTimer::tic (item_NonBondedInteraction);
      ptr_inter->applyNonBondedInteraction (sys.deviceData, relation);
      DeviceTimer::toc (item_NonBondedInteraction);
    }
    if (ptr_bdInterList != NULL){
      DeviceTimer::tic (item_BuildBondList);
      buildDeviceBondList (sys.deviceData, relation_buildBdList, *ptr_bdInterList);
      DeviceTimer::toc (item_BuildBondList);
      if (PCoupleOn){
	DeviceTimer::tic (item_BondedInterStatistic);
	ptr_inter->applyBondedInteraction (sys.deviceData, *ptr_bdInterList, myst);
	DeviceTimer::toc (item_BondedInterStatistic);
      }
      else {
	DeviceTimer::tic (item_BondedInteraction);
	ptr_inter->applyBondedInteraction (sys.deviceData, *ptr_bdInterList);	
	DeviceTimer::toc (item_BondedInteraction);
      }
    }
    HostTimer::tic (item_TransferGhost);
    sys.clearGhost();
    HostTimer::toc (item_TransferGhost);
  }
  else {
    firstStep (sys);
  }
}

