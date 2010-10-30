#include "MDSystem_interface.h"

#include "Reshuffle_interface.h"
#include "Reshuffle.h"
// #include "NonBondedInteraction.h"
// #include "HSFC.h"

void Reshuffle::
init (const MDSystem & sys,
      const NeighborList & nlist, 
      const IndexType & NTread)
{
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NTread;
  cellGridDim = nlist.cellGridDim;
  
  nob = cellGridDim.x * cellGridDim.y;;
  cudaMalloc ((void**)&posiBuff, sizeof(IndexType) * nob);
  cudaMalloc ((void**)&indexTable, sizeof(IndexType) * sys.hdata.numAtom);
}

Reshuffle::
~Reshuffle()
{
  cudaFree (posiBuff);
  cudaFree (indexTable);
}


static __global__ void
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

static __global__ void
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

// bool Reshuffle::
// calIndexTable (const NeighborList & nlist,
// 	       MDTimer * timer)
// {
//   if (nlist.mode != CellListBuilt) return false;

//   if (timer != NULL) timer->tic(mdTimeReshuffleSystem);
//   cellGridDim = nlist.cellGridDim;
//   IndexType nob = cellGridDim.x * cellGridDim.y;
//   cudaFree(posiBuff);
//   cudaMalloc ((void**)&posiBuff, sizeof(IndexType)*nob);  

//   Reshuffle_calPosiList
//       <<<1, 1>>> (
//    	  nlist.dclist.numbers, nob, posiBuff);
//   checkCUDAError ("Reshuffle::calIndexTable, cal posi");
//   Reshuffle_calIndexTable 
//       <<<cellGridDim, myBlockDim>>> (
//    	  nlist.dclist.data, posiBuff, indexTable);
//   checkCUDAError ("Reshuffle::calIndexTable, cal idxTable");
//   if (timer != NULL) timer->toc(mdTimeReshuffleSystem);
//   return true;
// }


						       
