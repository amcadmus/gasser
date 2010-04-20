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
easyMalloc (const IndexType & totalNumCell_,
		     const IndexType & maxNumBond_)
{
  totalNumCell = totalNumCell_;
  maxNumBond = maxNumBond_;
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  clear ();

  size_t size0 = sizeof(IndexType) * totalNumCell * numThreadsInCell * maxNumBond;
  size_t size1 = sizeof(IndexType) * totalNumCell * numThreadsInCell ;

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
  return localIndex * stride() + ithBond;
}

inline IndexType Parallel::DeviceBondList::
convertIndex (const IndexType & localIndex,
	      const IndexType & ithBond) const
{
  return localIndex * stride() + ithBond;
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
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  for (unsigned i = 0; i < totalNumCell; ++i){
    for (unsigned j = 0; j < numThreadsInCell; ++j){
      IndexType pIndex = i * numThreadsInCell + j;
      item_global_numBond (pIndex) = 0;
      for (unsigned k = 0; k < maxNumBond; ++k){
	item_global_neighborIndex (pIndex, k) = MaxIndexValue;
	item_global_bondIndex     (pIndex, k) = MaxIndexValue;
	item_neighborIndex	  (pIndex, k) = MaxIndexValue;
	// item_bondActive		  (pIndex, k) = MaxIndexValue;
      }
    }
  }
}

void Parallel::HostBondList::
reinit (HostCellListedMDData & hcellData,
	Topology::System & sysTop,
	SystemBondedInteraction & sysBdInter)
{
  maxNumBond = 0;
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
  totalNumCell = numCell.x * numCell.y * numCell.z;

  easyMalloc (totalNumCell, maxNumBond);
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
    : malloced (false), totalNumCell(0), maxNumBond(0)
{
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
easyMalloc (const IndexType & totalNumCell_,
	    const IndexType & maxNumBond_)
{
  clear ();

  totalNumCell = totalNumCell_;
  maxNumBond = maxNumBond_;
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();

  size_t size0 = sizeof(IndexType) * totalNumCell * numThreadsInCell * maxNumBond;
  size_t size1 = sizeof(IndexType) * totalNumCell * numThreadsInCell ;
  
  cudaMalloc ((void**)&global_numBond, size1);
  cudaMalloc ((void**)&global_neighborIndex, size0);
  cudaMalloc ((void**)&global_bondIndex, size0);
  cudaMalloc ((void**)&neighborIndex, size0);

  checkCUDAError ("DeviceBondList::easyMalloc");
}

void Parallel::DeviceBondList::
copyFromHost (const HostBondList & hbdlist)
{
  if (totalNumCell < hbdlist.totalNumCell ||
      totalNumCell * maxNumBond < hbdlist.totalNumCell * hbdlist.maxNumBond){
    easyMalloc (hbdlist.totalNumCell, hbdlist.maxNumBond);
  }
  
  totalNumCell = hbdlist.totalNumCell;
  maxNumBond = hbdlist.maxNumBond;
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();

  size_t size0 = sizeof(IndexType) * totalNumCell * numThreadsInCell * maxNumBond;
  size_t size1 = sizeof(IndexType) * totalNumCell * numThreadsInCell ;
  
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
  if (hbdlist.totalNumCell < totalNumCell ||
      hbdlist.totalNumCell * hbdlist.maxNumBond < totalNumCell * maxNumBond){
    hbdlist.easyMalloc (totalNumCell, maxNumBond);
  }

  hbdlist.totalNumCell = totalNumCell;
  hbdlist.maxNumBond = maxNumBond;
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();

  size_t size0 = sizeof(IndexType) * totalNumCell * numThreadsInCell * maxNumBond;
  size_t size1 = sizeof(IndexType) * totalNumCell * numThreadsInCell ;

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


