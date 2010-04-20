#define DEVICE_CODE

#include "Parallel_BondList.h"
#include "Parallel_Interface.h"

#include "compile_error_mixcode.h"

Parallel::HostBondList::
HostBondList ()
    : global_numBond (NULL),
      global_neighborIndex (NULL),
      global_bondIndex (NULL),
      neighborIndex (NULL)
      // bondActive (NULL)
{
}

Parallel::HostBondList::
~HostBondList()
{
  clear ();
}

void Parallel::HostBondList::
clear ()
{
  freeAPointer ((void**)&global_numBond);
  freeAPointer ((void**)&global_neighborIndex);
  freeAPointer ((void**)&global_bondIndex);
  freeAPointer ((void**)&neighborIndex);
  // freeAPointer ((void**)&bondActive);
}

void Parallel::HostBondList::
easyMalloc (const IndexType & stride,
	    const IndexType & length)
{
  myStride = stride;
  maxLength = length;  
  clear ();

  size_t size0 = sizeof(IndexType) * myStride * maxLength;
  size_t size1 = sizeof(IndexType) * myStride ;

  global_numBond = (IndexType*) malloc (size1);
  if (global_numBond == NULL){
    throw MDExcptFailedMallocOnHost ("HostBondList::easyMalloc",
				     "global_numBond", size1);
  }
  global_neighborIndex = (IndexType*) malloc (size0);
  if (global_neighborIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostBondList::easyMalloc",
				     "global_neighborIndex", size0);
  }
  global_bondIndex = (IndexType*) malloc (size0);
  if (global_bondIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostBondList::easyMalloc",
				     "global_bondIndex", size0);
  }
  neighborIndex = (IndexType*) malloc (size0);
  if (neighborIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostBondList::easyMalloc",
				     "neighborIndex", size0);
  }
  // bondActive = (IndexType*) malloc (size0);
  // if (bondActive == NULL){
  //   throw MDExcptFailedMallocOnHost ("HostBondList::easyMalloc",
  // 				     "bondActive", size0);
  // }
}

inline IndexType Parallel::HostBondList::
convertIndex (const IndexType & localIndex,
	      const IndexType & ithBond) const
{
  return ithBond * stride() + localIndex;
}

inline IndexType Parallel::DeviceBondList::
convertIndex (const IndexType & localIndex,
	      const IndexType & ithBond) const
{
  return ithBond * stride() + localIndex;
}

inline IndexType & Parallel::HostBondList::
item_global_numBond (const IndexType & localIndex)
{
  return global_numBond[localIndex];
}

inline IndexType & Parallel::HostBondList::
item_global_neighborIndex (const IndexType & localIndex,
			   const IndexType & ithBond)
{
  return global_neighborIndex[convertIndex(localIndex, ithBond)];
}

inline IndexType & Parallel::HostBondList::
item_global_bondIndex (const IndexType & localIndex,
		       const IndexType & ithBond)
{
  return global_bondIndex[convertIndex(localIndex, ithBond)];
}

inline IndexType & Parallel::HostBondList::
item_neighborIndex (const IndexType & localIndex,
		       const IndexType & ithBond)
{
  return neighborIndex[convertIndex(localIndex, ithBond)];
}

// inline IndexType & Parallel::HostBondList::
// item_bondActive (const IndexType & localIndex,
// 		       const IndexType & ithBond)
// {
//   return bondActive[convertIndex(localIndex, ithBond)];
// }


inline const IndexType & Parallel::HostBondList::
item_global_numBond (const IndexType & localIndex) const
{
  return global_numBond[localIndex];
}

inline const IndexType & Parallel::HostBondList::
item_global_neighborIndex (const IndexType & localIndex,
			   const IndexType & ithBond) const
{
  return global_neighborIndex[convertIndex(localIndex, ithBond)];
}

inline const IndexType & Parallel::HostBondList::
item_global_bondIndex (const IndexType & localIndex,
		       const IndexType & ithBond) const
{
  return global_bondIndex[convertIndex(localIndex, ithBond)];
}

inline const IndexType & Parallel::HostBondList::
item_neighborIndex (const IndexType & localIndex,
		    const IndexType & ithBond) const
{
  return neighborIndex[convertIndex(localIndex, ithBond)];
}

// inline const IndexType & Parallel::HostBondList::
// item_bondActive (const IndexType & localIndex,
// 		       const IndexType & ithBond) const
// {
//   return bondActive[convertIndex(localIndex, ithBond)];
// }


void Parallel::HostBondList::
initZero ()
{
  for (unsigned k = 0; k < maxLength; ++k){
    for (unsigned pIndex = 0; pIndex < myStride; ++pIndex){
      item_global_numBond (pIndex) = 0;
      item_global_neighborIndex (pIndex, k) = MaxIndexValue;
      item_global_bondIndex     (pIndex, k) = MaxIndexValue;
      item_neighborIndex	  (pIndex, k) = MaxIndexValue;
      // item_bondActive		  (pIndex, k) = MaxIndexValue;
    }
  }
}

void Parallel::HostBondList::
reinit (HostCellListedMDData & hcellData,
	Topology::System & sysTop,
	SystemBondedInteraction & sysBdInter)
{
  IndexType maxNumBond = 0;
  for (unsigned i = 0; i < sysBdInter.bondIndex.size(); ++i){
    for (unsigned j = 0; j < sysBdInter.bondIndex[i].size(); ++j){
      IndexType c ;
      if ((c=sysBdInter.bondIndex[i][j].size()) > maxNumBond){
	maxNumBond = c;
      }
    }
  }
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell ();
  HostIntVectorType numCell = hcellData.getNumCell();
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  myStride = totalNumCell * numThreadsInCell;
  maxLength = maxNumBond;

  easyMalloc (myStride, maxLength);
  initZero ();
  
  for (unsigned i = 0; i < totalNumCell; ++i){
    for (unsigned j = 0; j < hcellData.cptr_numAtomInCell()[i]; ++j){
      IndexType pIndex  = i * numThreadsInCell + j;
      IndexType pGlobal = hcellData.cptr_globalIndex()[pIndex];
      IndexType pTopMolIndex;
      IndexType pTopAtomIndex;
      sysTop.calMolTopPosition (pGlobal, pTopMolIndex, pTopAtomIndex);
      const std::vector<IndexType > & topNeighborIndex (
	  sysBdInter.getTopBondNeighborIndex (pTopMolIndex, pTopAtomIndex));
      const std::vector<IndexType > & topBondIndex (
	  sysBdInter.getTopBondIndex (pTopMolIndex, pTopAtomIndex));
      for (unsigned k = 0; k < topNeighborIndex.size(); ++k){
	item_global_neighborIndex(pIndex, k) = pGlobal + topNeighborIndex[k] - pTopAtomIndex;
	item_global_bondIndex(pIndex, k) = topBondIndex[k];
      }
      item_global_numBond(pIndex) = topNeighborIndex.size();
    }
  }
}


Parallel::DeviceBondList::
DeviceBondList ()
    : malloced (false), myStride(0), maxLength(0)
{
}

Parallel::DeviceBondList::
DeviceBondList (const HostBondList & hbdlist)
    : malloced (false), myStride(0), maxLength(0)
{
  reinit (hbdlist);
}

inline void Parallel::DeviceBondList::
reinit (const HostBondList & hbdlist)
{
  copyFromHost (hbdlist);
}

Parallel::DeviceBondList::
~DeviceBondList()
{
  clear ();
}

void Parallel::DeviceBondList::
clear ()
{
  if (malloced){
    cudaFree (global_numBond);
    cudaFree (global_neighborIndex);
    cudaFree (global_bondIndex);
    cudaFree (neighborIndex);
    checkCUDAError ("DeviceBondList::clear");
    malloced = false;
  }
}

void Parallel::DeviceBondList::
easyMalloc (const IndexType & stride,
	    const IndexType & length)
{
  clear ();
  myStride = stride;
  maxLength = length;
  
  size_t size0 = sizeof(IndexType) * myStride * maxLength;
  size_t size1 = sizeof(IndexType) * myStride ;
  
  cudaMalloc ((void**)&global_numBond, size1);
  cudaMalloc ((void**)&global_neighborIndex, size0);
  cudaMalloc ((void**)&global_bondIndex, size0);
  cudaMalloc ((void**)&neighborIndex, size0);

  checkCUDAError ("DeviceBondList::easyMalloc");
}

void Parallel::DeviceBondList::
copyFromHost (const HostBondList & hbdlist)
{
  if (myStride < hbdlist.myStride ||
      myStride * maxLength < hbdlist.myStride * hbdlist.maxLength){
    easyMalloc (hbdlist.myStride, hbdlist.maxLength);
  }  

  size_t size0 = sizeof(IndexType) * myStride * maxLength;
  size_t size1 = sizeof(IndexType) * myStride ;
  
  cudaMemcpy (global_numBond, hbdlist.global_numBond, size1,
	      cudaMemcpyHostToDevice);
  cudaMemcpy (global_neighborIndex, hbdlist.global_neighborIndex, size0,
	      cudaMemcpyHostToDevice);
  cudaMemcpy (global_bondIndex, hbdlist.global_bondIndex, size0,
	      cudaMemcpyHostToDevice);
  cudaMemcpy (neighborIndex, hbdlist.neighborIndex, size0,
	      cudaMemcpyHostToDevice);

  checkCUDAError ("DeviceBondList::copyFromHost");
}

void Parallel::DeviceBondList::
copyToHost (HostBondList & hbdlist) const 
{
  if (hbdlist.myStride < myStride ||
      hbdlist.myStride * hbdlist.maxLength < myStride * maxLength){
    hbdlist.easyMalloc (myStride, maxLength);
  }

  size_t size0 = sizeof(IndexType) * myStride * maxLength;
  size_t size1 = sizeof(IndexType) * myStride ;

  cudaMemcpy (hbdlist.global_numBond, global_numBond, size1,
	      cudaMemcpyDeviceToHost);
  cudaMemcpy (hbdlist.global_neighborIndex, global_neighborIndex, size0,
	      cudaMemcpyDeviceToHost);
  cudaMemcpy (hbdlist.global_bondIndex, global_bondIndex, size0,
	      cudaMemcpyDeviceToHost);
  cudaMemcpy (hbdlist.neighborIndex, neighborIndex, size0,
	      cudaMemcpyDeviceToHost);

  checkCUDAError ("DeviceBondList::copyToHost");
}


void Parallel::
buildDeviceBondList (const DeviceCellListedMDData & ddata,
		     const DeviceCellRelation & relation,
		     DeviceBondList & dbdlist)
{
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  IndexType totalNumCell = ddata.getNumCell().x * ddata.getNumCell().y * ddata.getNumCell().z;
  dim3 gridDim = toGridDim (totalNumCell);
  size_t sbuff_size = sizeof(IndexType) * numThreadsInCell;
  
  Parallel::CudaGlobal::buildDeviceBondList
      <<<gridDim, numThreadsInCell, sbuff_size>>> (
	  ddata.dptr_numAtomInCell(),
	  ddata.dptr_globalIndex(),
	  relation.dptr_numNeighborCell(),
	  relation.dptr_neighborCellIndex(),
	  relation.stride_neighborCellIndex(),
	  dbdlist.dptr_global_neighborIndex(),
	  dbdlist.dptr_global_numBond(),
	  dbdlist.stride (),
	  dbdlist.dptr_neighborIndex());  
}


__global__ void Parallel::CudaGlobal::
buildDeviceBondList (const IndexType * numAtomInCell,
		     const IndexType * globalIndex,
		     const IndexType * numNeighborCell,
		     const IndexType * neighborCellIndex,
		     const IndexType   cellRelationStride,
		     const IndexType * global_neighborIndex,
		     const IndexType * global_numBond,
		     const IndexType   bondListStride,
		     IndexType * neighborIndex)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  IndexType this_numAtom = numAtomInCell[bid];
  if (this_numAtom == 0) return;
  IndexType this_numNeighborCell = numNeighborCell[bid];
  if (this_numNeighborCell == 0) return;
  IndexType my_numBond;
  if (tid < this_numAtom) my_numBond = global_numBond[ii];

  extern __shared__ ScalorType buff_globalIndex[];
  

  for (IndexType kk = 0; kk < this_numNeighborCell; ++kk){
    __syncthreads();
    IndexType target_cellIndex = neighborCellIndex[bid * cellRelationStride + kk];
    IndexType indexShift = target_cellIndex * blockDim.x;
    IndexType target_numAtomInCell = numAtomInCell[target_cellIndex];
    IndexType jj = indexShift + tid;
    if (tid < target_numAtomInCell){
      buff_globalIndex[tid] = globalIndex[jj];
    }
    __syncthreads();
    
    if (tid < this_numAtom){
      for (IndexType ll = 0; ll < my_numBond; ++ll){
	IndexType tofind_globalIndex =
	    global_neighborIndex[
		Parallel::DeviceBondList_cudaDevice::
		indexConvert(bondListStride, ii, ll)
		];
	for (IndexType mm = 0; mm < target_numAtomInCell; ++mm){
	  if (tofind_globalIndex == buff_globalIndex[mm]){
	    neighborIndex[
		Parallel::DeviceBondList_cudaDevice::
		indexConvert(bondListStride, ii, ll)
		] = mm + indexShift;
	  }
	}
      }
    }
  }
}


  

  
