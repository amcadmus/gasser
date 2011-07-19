#include "MDSystem_interface.h"

#include "Reshuffle_interface.h"
#include "Reshuffle.h"
// #include "NonBondedInteraction.h"
// #include "HSFC.h"

void Reshuffle::
reinit (const MDSystem & sys)
{
  clear();
  cudaMalloc ((void**)&indexTable, sizeof(IndexType) * sys.hdata.numAtom);
  checkCUDAError ("Reshuffle::reinit");
}

Reshuffle::
Reshuffle (const MDSystem & sys)
    : malloced(false)
{
  reinit (sys);
}

Reshuffle::
~Reshuffle()
{
  clear();
}

void Reshuffle::
clear()
{
  if (malloced){
    cudaFree (indexTable);
    malloced = false;
  }
}


__global__ void
Reshuffle_calPosiList  (const IndexType * cellNumbers,
			const IndexType nob,
			IndexType * posiList)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  if (tid + bid == 0){
    posiList[0] = 0;
    for (IndexType ii = 1; ii < nob; ++ii){
      posiList[ii] = posiList[ii-1] + cellNumbers[ii-1];
    }
  }
}

__global__ void
Reshuffle_calIndexTable (const IndexType * clistData,
			 const IndexType * posiList,
			 IndexType * idxTable)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  IndexType fromid = clistData[bid * blockDim.x + tid];
  IndexType toid = posiList[bid] + tid;
  if (fromid != MaxIndexValue){
    idxTable[fromid] = toid;
  }
}

bool Reshuffle::
calIndexTable (const CellList & clist,
	       MDTimer * timer)
{
  if (clist.isempty()) return false;

  if (timer != NULL) timer->tic(mdTimeReshuffleSystem);
  dim3 cellGridDim = clist.getCellGrimDim();
  IndexType nob = cellGridDim.x * cellGridDim.y;
  IndexType * posiBuff;
  cudaMalloc ((void**)&posiBuff, sizeof(IndexType)*nob);  
  checkCUDAError ("Reshuffle::calIndexTable malloc posiBuff");
  
  Reshuffle_calPosiList
      <<<1, 1>>> (
   	  clist.dclist.numbers,
	  nob,
	  posiBuff);
  checkCUDAError ("Reshuffle::calIndexTable, cal posi");
  Reshuffle_calIndexTable 
      <<<clist.getCellGrimDim(), clist.getCellBlockDim()>>> (
   	  clist.dclist.data,
	  posiBuff,
	  indexTable);
  checkCUDAError ("Reshuffle::calIndexTable, cal idxTable");
  
  cudaFree(posiBuff);
  checkCUDAError ("Reshuffle::calIndexTable free posiBuff");
  if (timer != NULL) timer->toc(mdTimeReshuffleSystem);

  return true;
}
