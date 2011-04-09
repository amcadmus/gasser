#define DEVICE_CODE

#include "Statistic_interface.h"


void MDStatistic::deviceCopy (const MDStatistic & st)
{
  cudaMemcpy (ddata, st.ddata, sizeof(ScalorType) * NumberOfStatisticItems,
	      cudaMemcpyDeviceToDevice);
  checkCUDAError ("Statistic::deviceCopy");
}

MDStatistic::
MDStatistic ()
    : hdata (NULL), dmalloced (false)
{
}

MDStatistic::
MDStatistic (const MDSystem & sys)
    : hdata (NULL), dmalloced (false)
{
  reinit(sys);
}

void MDStatistic::
clear ()
{
  freeAPointer ((void**)&hdata);
  if (dmalloced){
    cudaFree (ddata);
    dmalloced = false;
  }
}


void MDStatistic::
reinit (const MDSystem & sys)
{
  clear ();
  
  // malloc and init system
  hdata = (ScalorType *) malloc (sizeof(ScalorType) * NumberOfStatisticItems);
  if (hdata == NULL){
    throw MDExcptFailedMallocOnHost ("MDStatistic::MDStatistic", "hdata",
				     sizeof(ScalorType) * NumberOfStatisticItems);
  }
  cudaMalloc ((void**)&ddata, sizeof(ScalorType) * NumberOfStatisticItems);
  checkCUDAError("MDStatistic::MDStatistic allocate for ddata");
  dmalloced = true;

  for (IndexType i = 0; i < NumberOfStatisticItems; ++i){
    hdata[i] = 0.f;
  }
  cudaMemcpy (ddata, hdata, sizeof(ScalorType) * NumberOfStatisticItems,
	      cudaMemcpyHostToDevice);
  checkCUDAError("MDStatistic::MDStatistic cpy form host to device");
}

MDStatistic::
~MDStatistic ()
{
  clear ();
}

__global__ void clearStatisticData (ScalorType *ddata)
{
  if (threadIdx.x < NumberOfStatisticItems){
    ddata[threadIdx.x] = 0.f;
  }
}

__global__ void addStatisticData (ScalorType * ddata, const ScalorType * cddata)
{
  if (threadIdx.x < NumberOfStatisticItems){
    ddata[threadIdx.x] += cddata[threadIdx.x];
  }
}

void MDStatistic::
clearDevice ()
{
  clearStatisticData <<<1, NumberOfStatisticItems>>> (ddata);
  checkCUDAError("Statistic::clearDevice");
}

void MDStatistic::
updateHost () const
{
  cudaMemcpy (hdata, ddata, sizeof(ScalorType) * NumberOfStatisticItems, 
	      cudaMemcpyDeviceToHost);
  checkCUDAError("Statistic::updateHost");
}

void MDStatistic::
deviceAdd (const MDStatistic & st)
{
  addStatisticData <<<1, NumberOfStatisticItems>>> (ddata, st.ddata);
}


void MDStatistic::
copy (const MDStatistic & st,
      const IndexType num,
      const mdStatisticItem_t items[NumberOfStatisticItems])
{
  ScalorType tmp[NumberOfStatisticItems];
  cudaMemcpy (tmp, st.ddata, sizeof(ScalorType) * NumberOfStatisticItems,
	      cudaMemcpyDeviceToHost);
  checkCUDAError ("MDStatistic::deviceCopy");
  updateHost();
  for (unsigned i = 0; i < num; ++i){
    hdata[items[i]] = tmp[items[i]];
  }
  cudaMemcpy (ddata, hdata, sizeof(ScalorType) * NumberOfStatisticItems,
	      cudaMemcpyHostToDevice);
}

void MDStatistic::
add  (const MDStatistic & st,
      const IndexType num,
      const mdStatisticItem_t items[NumberOfStatisticItems])
{
  ScalorType tmp[NumberOfStatisticItems];
  cudaMemcpy (tmp, st.ddata, sizeof(ScalorType) * NumberOfStatisticItems,
	      cudaMemcpyDeviceToHost);
  checkCUDAError ("MDStatistic::deviceCopy");
  updateHost();
  for (unsigned i = 0; i < num; ++i){
    hdata[items[i]] += tmp[items[i]];
  }
  cudaMemcpy (ddata, hdata, sizeof(ScalorType) * NumberOfStatisticItems,
	      cudaMemcpyHostToDevice);
}

  
__global__ static void
syncData (ScalorType * tmpddata,
	  ScalorType * ddata,
	  int * flag)
{
  if (threadIdx.x < NumberOfStatisticItems){
    int myflag = flag[threadIdx.x];
    if (myflag != -1){
      ddata[myflag] = tmpddata[myflag];
    }
  }
}

void MDStatistic::
setEnergyCorr (const ScalorType & energyCorr_)
{
  hdata[mdStatisticEnergyCorrection] = energyCorr_;
  int flag[NumberOfStatisticItems];
  for (unsigned i = 0; i < NumberOfStatisticItems; ++i){
    flag[i] = -1;
  }
  flag[0] = mdStatisticEnergyCorrection;
  ScalorType * tmpddata;
  int * dflag;
  cudaMalloc ((void**)&tmpddata, sizeof(ScalorType) * NumberOfStatisticItems);
  cudaMalloc ((void**)&dflag, sizeof(int) * NumberOfStatisticItems);
  checkCUDAError("MDStatistic::setEnergyCorr allocate for tmpddata");
  cudaMemcpy (tmpddata, hdata, sizeof(ScalorType) * NumberOfStatisticItems,
	      cudaMemcpyHostToDevice);
  cudaMemcpy (dflag, flag, sizeof(int) * NumberOfStatisticItems,
	      cudaMemcpyHostToDevice);
  syncData<<<1, NumberOfStatisticItems>>> (tmpddata, ddata, dflag);
  checkCUDAError("MDStatistic::setEnergyCorr sync for tmpddata");
  cudaFree (tmpddata);
  cudaFree (dflag);
}

void MDStatistic::
setPressureCorr (const ScalorType & pressureCorr_)
{
  hdata[mdStatisticPressureCorrection] = pressureCorr_;
  // printf ("# setting pressureCorr_ to %f\n", pressureCorr_);
  int flag[NumberOfStatisticItems];
  for (unsigned i = 0; i < NumberOfStatisticItems; ++i){
    flag[i] = -1;
  }
  flag[0] = mdStatisticPressureCorrection;
  ScalorType * tmpddata;
  int * dflag;
  cudaMalloc ((void**)&tmpddata, sizeof(ScalorType) * NumberOfStatisticItems);
  cudaMalloc ((void**)&dflag, sizeof(int) * NumberOfStatisticItems);
  checkCUDAError("MDStatistic::setPressureCorr allocate for tmpddata");
  cudaMemcpy (tmpddata, hdata, sizeof(ScalorType) * NumberOfStatisticItems,
	      cudaMemcpyHostToDevice);
  cudaMemcpy (dflag, flag, sizeof(int) * NumberOfStatisticItems,
	      cudaMemcpyHostToDevice);
  syncData<<<1, NumberOfStatisticItems>>> (tmpddata, ddata, dflag);
  checkCUDAError("MDStatistic::setPressureCorr sync for tmpddata");
  cudaFree (tmpddata);
  cudaFree (dflag);
}
