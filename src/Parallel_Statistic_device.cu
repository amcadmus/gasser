#define DEVICE_CODE

#include "Parallel_Statistic.h"
#include "compile_error_mixcode.h"

Parallel::DeviceStatistic::
~DeviceStatistic ()
{
  if (dmalloced){
    cudaFree (ddata);
  }
}


void Parallel::DeviceStatistic::
reinit (const DeviceCellListedMDData & data)
{
  volume = data.getGlobalBox().size.x * data.getGlobalBox().size.y *
      data.getGlobalBox().size.z;

  size = sizeof (ScalorType) * NumberOfStatisticItems;
  if (!dmalloced){
    cudaMalloc ((void**)&ddata, size);
    checkCUDAError("DeviceStatistic::init, malloc");
    dmalloced = true;
  }
  
  clearData ();
}

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

