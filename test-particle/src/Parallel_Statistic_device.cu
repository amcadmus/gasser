#define DEVICE_CODE

#include "Parallel_Statistic.h"
#include "compile_error_mixcode.h"

Parallel::DeviceStatistic::
~DeviceStatistic ()
{
  clear ();
}

void Parallel::DeviceStatistic::
clear ()
{
  if (dmalloced){
    cudaFree (ddata);
    dmalloced = false;
  }
}

Parallel::DeviceStatistic::
DeviceStatistic ()
    : dmalloced (false)
{
  size = sizeof (ScalorType) * NumberOfStatisticItems;
  cudaMalloc ((void**)&ddata, size);
  checkCUDAError("DeviceStatistic::init, malloc");
  dmalloced = true;
  clearData ();
}

// void Parallel::DeviceStatistic::
// reinit (const DeviceCellListedMDData & ddata)
// {
//   clearData ();
// }

__global__ void Parallel::CudaGlobal::
clearStatisticData (ScalorType *ddata)
{
  if (threadIdx.x < NumberOfStatisticItems){
    ddata[threadIdx.x] = 0.f;
  }
}

__global__ void Parallel::CudaGlobal::
addStatisticData (ScalorType * ddata, const ScalorType * cddata)
{
  if (threadIdx.x < NumberOfStatisticItems){
    ddata[threadIdx.x] += cddata[threadIdx.x];
  }
}

void Parallel::DeviceStatistic::
clearData ()
{
  Parallel::CudaGlobal::clearStatisticData <<<1, NumberOfStatisticItems>>> (ddata);
  checkCUDAError("DeviceStatistic::clearDevice");
}

void Parallel::DeviceStatistic::
add (const DeviceStatistic & st)
{
  Parallel::CudaGlobal::addStatisticData <<<1, NumberOfStatisticItems>>> (ddata, st.ddata);
  checkCUDAError("DeviceStatistic::add");
}

void Parallel::DeviceStatistic::
copy (const DeviceStatistic & st)
{
  cudaMemcpy (ddata, st.ddata, sizeof(ScalorType) * NumberOfStatisticItems,
	      cudaMemcpyDeviceToDevice);
  checkCUDAError ("DeviceStatistic::copy");
}

void Parallel::DeviceStatistic::
copyToHost (HostStatistic & hst)
{
  cudaMemcpy (hst.cptr_localStatisticData(), ddata, size, cudaMemcpyDeviceToHost);
  checkCUDAError ("DeviceStatistic::copyToHost");
}

static __global__ void
localRescale (const IndexType  posistion,
	      const ScalorType scale,
	      ScalorType * ddata)
{
  ddata[posistion + threadIdx.x] *= scale;
}

void Parallel::DeviceStatistic::
rescale (const IndexType  & position,
	 const IndexType  & num,
	 const ScalorType & scale)
{
  localRescale <<<1, num>>> (position, scale, ddata);
}

