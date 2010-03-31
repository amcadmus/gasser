#define DEVICE_CODE

#include "Parallel_Integrator.h"
#include "Parallel_Interface.h"
#include "Parallel_Auxiliary.h"
#include "compile_error_mixcode.h"

Parallel::TranslationalFreedomRemover::
~TranslationalFreedomRemover ()
{
  if (malloced){
    cudaFree (sums);
  }
}

void Parallel::TranslationalFreedomRemover::
reinit  (const DeviceCellListedMDData & data)
{
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
    malloced = true;
  }

  SumVector<ScalorType > sum_mass;
  sum_mass.reinit (totalNumCell, NThreadForSum);
  Parallel::CudaGlobal::prepareCalTotalMass
      <<<gridDim, numThreadsInCell, sharedBuffSize>>>(
	  data.dptr_numAtomInCell(),
	  data.dptr_mass(),
	  sum_mass.getBuff());
  sum_mass.sumBuff (sumM, 0);
  cudaMemcpy (&totalMassi, sumM, sizeof(ScalorType), cudaMemcpyDeviceToHost);
  checkCUDAError ("TranslationalFreedomRemover::reinit, cpy sumM");
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
	  sum_x.getBuff(),
	  sum_y.getBuff(),
	  sum_z.getBuff());
  checkCUDAError ("TranslationalFreedomRemover::remove, prepare");
  sum_x.sumBuff (sums, 0);
  sum_y.sumBuff (sums, 1);
  sum_z.sumBuff (sums, 2);
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
init (const DeviceCellListedMDData & ddata)
{
  IndexType totalNumCell = ddata.getNumCell().x *
      ddata.getNumCell().y * ddata.getNumCell().z;
  gridDim = toGridDim (totalNumCell);

  sum_kxx.init (totalNumCell, NThreadForSum);
  sum_kyy.init (totalNumCell, NThreadForSum);
  sum_kzz.init (totalNumCell, NThreadForSum);
  
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
  checkCUDAError ("interface::VelocityVerlet::step2, no st");
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
	  sum_kxx.getBuff(),
	  sum_kyy.getBuff(),
	  sum_kzz.getBuff());
  sum_kxx.sumBuffAdd (st.dptr_statisticData(), mdStatistic_KineticEnergyXX, 0);
  sum_kyy.sumBuffAdd (st.dptr_statisticData(), mdStatistic_KineticEnergyYY, 0);
  sum_kzz.sumBuffAdd (st.dptr_statisticData(), mdStatistic_KineticEnergyZZ, 0);

  checkCUDAError ("interface::VelocityVerlet::step2, with st");
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

  
