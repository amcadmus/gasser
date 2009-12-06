#include "Statistic_interface.h"


// __global__ void initBuff (ScalorType * buff, IndexType n);
// void Statistic::init(const MDSystem & sys, 
// 		     const IndexType & NThread)
// {
//   myBlockDim.y = 1;
//   myBlockDim.z = 1;
//   myBlockDim.x = NThread;
//   IndexType nob;
//   if (sys.ddata.numAtom % myBlockDim.x == 0){
//     nob = sys.ddata.numAtom / myBlockDim.x;
//   } else {
//     nob = sys.ddata.numAtom / myBlockDim.x + 1;
//   }
//   atomGridDim = toGridDim (nob);

//   hostData.data = (ScalorType *) malloc (sizeof(ScalorType) * Size_StatisticData);

//   cudaMalloc((void**)&(deviceData.data), sizeof(ScalorType) * Size_StatisticData);
//   clearStatisticData <<<1, 1>>> (deviceData);

//   updateHost();
  
//   cudaMalloc((void**)&(statistic_buff), sizeof(ScalorType) * nob);
//   initBuff <<<1, 1>>> (statistic_buff, nob);
//   checkCUDAError("Statistic::init");
// }

// void Statistic::clearDevice ()
// {
//   clearStatisticData <<<1, 1>>> (deviceData);
//   checkCUDAError("Statistic::clearDevice");
// }

// void Statistic::updateHost()
// {
//   cudaMemcpy (hostData.data, deviceData.data, 
// 	      sizeof(ScalorType) * Size_StatisticData, cudaMemcpyDeviceToHost);
//   checkCUDAError("Statistic::updateHost");
// }

// Statistic::~Statistic()
// {
//   free (hostData.data);
//   cudaFree (deviceData.data);
//   cudaFree (statistic_buff);
//   checkCUDAError("Statistic::~Statistic");
// }


// __global__ void initBuff (ScalorType * buff, IndexType n)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   if (bid + tid == 0)
//     for (IndexType i = 0; i < n; ++i){
//       buff[i] = 0;
//     }
// }


MDStatistic::
MDStatistic ()
{
  hdata = NULL;
  dmalloced = false;
  volume = 0;
}


void MDStatistic::
init (const MDSystem & sys)
{
  // recorde system infomation
  volume = sys.box.size.x * sys.box.size.y * sys.box.size.z;

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
  freeAPointer ((void**)&hdata);
  if (dmalloced){
    cudaFree (ddata);
  }
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
updateHost ()
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

  
