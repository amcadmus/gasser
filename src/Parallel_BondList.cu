#define DEVICE_CODE

#include "Parallel_BondList.h"
#include "Parallel_Interface.h"
#include "Parallel_Auxiliary.h"

#include "compile_error_mixcode.h"

Parallel::HostBondList::
HostBondList ()
    : stride(0),
      maxNumBond(0),
      maxNumAngle(0),
      maxNumDihedral(0),
      bondNeighbor_localIndex(NULL),
      angleNeighbor_localIndex(NULL),
      dihedralNeighbor_localIndex(NULL)
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
  stride = 0;
  maxNumBond = 0;
  maxNumAngle = 0;
  maxNumDihedral = 0;
  
  freeAPointer ((void**)&bondNeighbor_localIndex);
  freeAPointer ((void**)&angleNeighbor_localIndex);
  freeAPointer ((void**)&dihedralNeighbor_localIndex);
}

void Parallel::HostBondList::
easyMalloc (const IndexType & stride_,
	    const IndexType & maxNumBond_,
	    const IndexType & maxNumAngle_,
	    const IndexType & maxNumDihedral_)
{
  clear ();

  stride = stride_;
  maxNumBond = maxNumBond_;
  maxNumAngle = maxNumAngle_;
  maxNumDihedral = maxNumDihedral_;

  size_t size;
  size = sizeof(IndexType) * stride * maxNumBond;
  bondNeighbor_localIndex = (IndexType*) malloc (size);
  if (bondNeighbor_localIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostBondList::easyMalloc",
				     "bondNeighbor_localIndex", size);
  }
  size = sizeof(IndexType) * stride * maxNumAngle * 2;
  angleNeighbor_localIndex = (IndexType*) malloc (size);
  if (angleNeighbor_localIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostBondList::easyMalloc",
				     "angleNeighbor_localIndex", size);
  }
  size = sizeof(IndexType) * stride * maxNumDihedral * 3;
  dihedralNeighbor_localIndex = (IndexType*) malloc (size);
  if (dihedralNeighbor_localIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostBondList::easyMalloc",
				     "dihedralNeighbor_localIndex", size);
  }
}


void Parallel::HostBondList::
fillzero ()
{
  for (unsigned i = 0; i < maxNumBond * stride; ++i){
      bondNeighbor_localIndex[i] = 0;
  }
  for (unsigned i = 0; i < maxNumAngle * stride; ++i){
      angleNeighbor_localIndex[i] = 0;
  }
  for (unsigned i = 0; i < maxNumDihedral * stride; ++i){
      dihedralNeighbor_localIndex[i] = 0;
  }  
}



Parallel::DeviceBondList::
DeviceBondList ()
    : malloced (false),
      stride(0),
      maxNumBond(0), 
      maxNumAngle(0), 
      maxNumDihedral(0)
{
}

Parallel::DeviceBondList::
DeviceBondList (const DeviceCellListedMDData & data)
    : malloced (false),
      stride(0),
      maxNumBond(0), 
      maxNumAngle(0), 
      maxNumDihedral(0)
{
  reinit (data);
}

void Parallel::DeviceBondList::
reinit (const DeviceCellListedMDData & data)
{
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  dim3 gridDim;

  easyMalloc (data.DeviceMDData::memSize(),
	      data.getMaxNumBond(),
	      data.getMaxNumAngle(),
	      data.getMaxNumDihedral());
  gridDim = toGridDim (maxNumBond * stride / numThreadsInCell + 1);
  Parallel::Auxiliary::setValue <<<gridDim, numThreadsInCell>>> (
      bondNeighbor_localIndex,
      maxNumBond * stride,
      IndexType (0));
  gridDim = toGridDim (maxNumAngle * stride * 2 / numThreadsInCell + 1);
  Parallel::Auxiliary::setValue <<<gridDim, numThreadsInCell>>> (
      angleNeighbor_localIndex,
      maxNumAngle * stride * 2,
      IndexType (0));
  gridDim = toGridDim (maxNumDihedral * stride * 3 / numThreadsInCell + 1);
  Parallel::Auxiliary::setValue <<<gridDim, numThreadsInCell>>> (
      dihedralNeighbor_localIndex,
      maxNumDihedral * stride * 3,
      IndexType (0));
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
    stride = 0;
    maxNumBond = 0;
    maxNumAngle = 0;
    maxNumDihedral = 0;
    cudaFree(bondNeighbor_localIndex);
    cudaFree(angleNeighbor_localIndex);
    cudaFree(dihedralNeighbor_localIndex);
    checkCUDAError ("DeviceBondList::clear");
    malloced = false;
  }
}

void Parallel::DeviceBondList::
easyMalloc (const IndexType & stride_,
	    const IndexType & maxNumBond_,
	    const IndexType & maxNumAngle_,
	    const IndexType & maxNumDihedral_)
{
  clear ();

  stride = stride_;
  maxNumBond = maxNumBond_;
  maxNumAngle = maxNumAngle_;
  maxNumDihedral = maxNumDihedral_;
  
  cudaMalloc ((void**)&bondNeighbor_localIndex,
	      sizeof(IndexType) * maxNumBond * stride);
  cudaMalloc ((void**)&angleNeighbor_localIndex,
	      sizeof(IndexType) * maxNumAngle * stride * 2);
  cudaMalloc ((void**)&dihedralNeighbor_localIndex,
	      sizeof(IndexType) * maxNumDihedral * stride * 3);
  
  checkCUDAError ("DeviceBondList::easyMalloc");
  malloced = true;
}

void Parallel::DeviceBondList::
copyFromHost (const HostBondList & hbdlist)
{
  if (getStride() != hbdlist.getStride() ||
      getMaxNumBond() != hbdlist.getMaxNumBond() ||
      getMaxNumAngle() != hbdlist.getMaxNumAngle() ||
      getMaxNumDihedral() != hbdlist.getMaxNumDihedral() ){
    easyMalloc (hbdlist.getStride(),
		hbdlist.getMaxNumBond(),
		hbdlist.getMaxNumAngle(),
		hbdlist.getMaxNumDihedral());
  }  

  size_t size;
  size = sizeof(IndexType) * stride * getMaxNumBond();
  cudaMemcpy (bondNeighbor_localIndex,
	      hbdlist.cptr_bondNeighbor_localIndex(),
	      size,
	      cudaMemcpyHostToDevice);
  checkCUDAError ("DeviceBondList::copyFromHost, bondNeighbor_localIndex");
  size = sizeof(IndexType) * stride * getMaxNumAngle() * 2;
  cudaMemcpy (angleNeighbor_localIndex,
	      hbdlist.cptr_angleNeighbor_localIndex(),
	      size,
	      cudaMemcpyHostToDevice);
  checkCUDAError ("DeviceBondList::copyFromHost, angleNeighbor_localIndex");
  size = sizeof(IndexType) * stride * getMaxNumDihedral() * 3;
  cudaMemcpy (dihedralNeighbor_localIndex,
	      hbdlist.cptr_dihedralNeighbor_localIndex(),
	      size,
	      cudaMemcpyHostToDevice);
  checkCUDAError ("DeviceBondList::copyFromHost, dihedralNeighbor_localIndex");
}

void Parallel::DeviceBondList::
copyToHost (HostBondList & hbdlist) const 
{
  if (getStride() != hbdlist.getStride() ||
      getMaxNumBond() != hbdlist.getMaxNumBond() ||
      getMaxNumAngle() != hbdlist.getMaxNumAngle() ||
      getMaxNumDihedral() != hbdlist.getMaxNumDihedral() ){
    hbdlist.easyMalloc (getStride(),
			getMaxNumBond(),
			getMaxNumAngle(),
			getMaxNumDihedral());
  }

  size_t size;
  size = sizeof(IndexType) * stride * getMaxNumBond();
  cudaMemcpy (hbdlist.cptr_bondNeighbor_localIndex(),
	      bondNeighbor_localIndex,
	      size,
	      cudaMemcpyDeviceToHost);
  checkCUDAError ("DeviceBondList::copyFromHost, bondNeighbor_localIndex");
  size = sizeof(IndexType) * stride * getMaxNumAngle() * 2;
  cudaMemcpy (hbdlist.cptr_angleNeighbor_localIndex(),
	      angleNeighbor_localIndex,
	      size,
	      cudaMemcpyDeviceToHost);
  checkCUDAError ("DeviceBondList::copyFromHost, angleNeighbor_localIndex");
  size = sizeof(IndexType) * stride * getMaxNumDihedral() * 3;
  cudaMemcpy (hbdlist.cptr_dihedralNeighbor_localIndex(),
	      dihedralNeighbor_localIndex,
	      size,
	      cudaMemcpyDeviceToHost);
  checkCUDAError ("DeviceBondList::copyFromHost, dihedralNeighbor_localIndex");
}



void Parallel::
buildDeviceBondList (const DeviceCellListedMDData & ddata,
		     const DeviceCellRelation & relation,
		     DeviceBondList & dbdlist)
{
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  IndexType totalNumCell =
      ddata.getNumCell().x * ddata.getNumCell().y * ddata.getNumCell().z;
  dim3 gridDim = toGridDim (totalNumCell);
  size_t sbuff_size = sizeof(IndexType) * numThreadsInCell;
  
  Parallel::CudaGlobal::buildDeviceBondList
      <<<gridDim, numThreadsInCell, sbuff_size>>> (
	  ddata.dptr_numAtomInCell(),
	  ddata.dptr_globalIndex(),
	  relation.dptr_numNeighborCell(),
	  relation.dptr_neighborCellIndex(),
	  relation.stride_neighborCellIndex(),
	  ddata.getMaxNumBond(),
	  ddata.dptr_numBond(),
	  ddata.dptr_bondNeighbor_globalIndex(),
	  dbdlist.dptr_bondNeighbor_localIndex(),
	  ddata.getMaxNumAngle(),
	  ddata.dptr_numAngle(),
	  ddata.dptr_angleNeighbor_globalIndex(),
	  dbdlist.dptr_angleNeighbor_localIndex(),
	  ddata.getMaxNumDihedral(),
	  ddata.dptr_numDihedral(),
	  ddata.dptr_dihedralNeighbor_globalIndex(),
	  dbdlist.dptr_dihedralNeighbor_localIndex(),
	  ddata.bondTopStride(),
	  dbdlist.getStride());  
}


__global__ void Parallel::CudaGlobal::
buildDeviceBondList (const IndexType * numAtomInCell,
		     const IndexType * globalIndex,
		     const IndexType * numNeighborCell,
		     const IndexType * neighborCellIndex,
		     const IndexType   cellRelationStride,
		     const IndexType   maxNumBond,
		     const IndexType * numBond,
		     const IndexType * bondNeighbor_globalIndex,
		     IndexType * bondNeighbor_localIndex,
		     const IndexType   maxNumAngle,
		     const IndexType * numAngle,
		     const IndexType * angleNeighbor_globalIndex,
		     IndexType * angleNeighbor_localIndex,
		     const IndexType   maxNumDihedral,
		     const IndexType * numDihedral,
		     const IndexType * dihedralNeighbor_globalIndex,
		     IndexType * dihedralNeighbor_localIndex,
		     const IndexType   bondTopStride,
		     const IndexType   listTopStride)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  extern __shared__ ScalorType buff_globalIndex[];  
  bool hasBond = (maxNumBond != 0);
  bool hasAngle = (maxNumAngle != 0);
  bool hasDihedral = (maxNumDihedral != 0);

  if (hasBond && hasAngle && hasDihedral) return;
  IndexType this_numAtom = numAtomInCell[bid];
  if (this_numAtom == 0) return;
  IndexType this_numNeighborCell = numNeighborCell[bid];
  if (this_numNeighborCell == 0) return;
  IndexType my_numBond(0), my_numAngle(0), my_numDihedral(0);
  if (hasBond && tid < this_numAtom) my_numBond = numBond[ii];
  if (hasAngle && tid < this_numAtom) my_numAngle = numAngle[ii];
  if (hasDihedral && tid < this_numAtom) my_numDihedral = numDihedral[ii];
  
  for (IndexType kk = 0; kk < this_numNeighborCell; ++kk){
    __syncthreads();
    IndexType target_cellIndex = neighborCellIndex[bid * cellRelationStride + kk];
    // if (bid == 31 && tid == 0) {
    //   printf ("targetcell id: %d\n", target_cellIndex);
    // }
    IndexType indexShift = target_cellIndex * blockDim.x;
    IndexType target_numAtomInCell = numAtomInCell[target_cellIndex];
    IndexType jj = indexShift + tid;
    if (tid < target_numAtomInCell){
      buff_globalIndex[tid] = globalIndex[jj];
    }
    __syncthreads();
    if (hasBond && tid < this_numAtom){
      IndexType topIndexShift = 0;
      IndexType listIndexShift = 0;
      for (IndexType ll = 0; ll < my_numBond; ++ll){
	IndexType tofind_globalIndex = bondNeighbor_globalIndex[topIndexShift + ii];
	for (IndexType mm = 0; mm < target_numAtomInCell; ++mm){
	  if (tofind_globalIndex == buff_globalIndex[mm]){
	    bondNeighbor_localIndex[listIndexShift + ii] = mm + indexShift;
	  }
	}
	topIndexShift += bondTopStride;
	listIndexShift += listTopStride;
      }
    }
    if (hasAngle && tid < this_numAtom){
      IndexType topIndexShift = 0;
      IndexType listIndexShift = 0;
      for (IndexType ll = 0; ll < my_numAngle * 2; ++ll){
	IndexType tofind_globalIndex = angleNeighbor_globalIndex[topIndexShift + ii];
	for (IndexType mm = 0; mm < target_numAtomInCell; ++mm){
	  if (tofind_globalIndex == buff_globalIndex[mm]){
	    angleNeighbor_localIndex[topIndexShift + ii] = mm + indexShift;
	  }
	}
	topIndexShift += bondTopStride;
	listIndexShift += listTopStride;
      }
    }
    if (hasDihedral && tid < this_numAtom){
      IndexType topIndexShift = 0;
      IndexType listIndexShift = 0;
      for (IndexType ll = 0; ll < my_numDihedral * 3; ++ll){
	IndexType tofind_globalIndex = dihedralNeighbor_globalIndex[topIndexShift + ii];
	for (IndexType mm = 0; mm < target_numAtomInCell; ++mm){
	  if (tofind_globalIndex == buff_globalIndex[mm]){
	    dihedralNeighbor_localIndex[topIndexShift + ii] = mm + indexShift;
	  }
	}
	topIndexShift += bondTopStride;
	listIndexShift += listTopStride;
      }
    }
  }
}




// buildDeviceBondList (const IndexType * numAtomInCell,
// 		     const IndexType * globalIndex,
// 		     const IndexType * numNeighborCell,
// 		     const IndexType * neighborCellIndex,
// 		     const IndexType   cellRelationStride,
// 		     const IndexType * global_neighborIndex,
// 		     const IndexType * global_numBond,
// 		     const IndexType   bondListStride,
// 		     IndexType * neighborIndex)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   IndexType ii = tid + bid * blockDim.x;

//   IndexType this_numAtom = numAtomInCell[bid];
//   if (this_numAtom == 0) return;
//   IndexType this_numNeighborCell = numNeighborCell[bid];
//   if (this_numNeighborCell == 0) return;
//   IndexType my_numBond;
//   if (tid < this_numAtom) my_numBond = global_numBond[ii];

//   extern __shared__ ScalorType buff_globalIndex[];  

//   for (IndexType kk = 0; kk < this_numNeighborCell; ++kk){
//     __syncthreads();
//     IndexType target_cellIndex = neighborCellIndex[bid * cellRelationStride + kk];
//     IndexType indexShift = target_cellIndex * blockDim.x;
//     IndexType target_numAtomInCell = numAtomInCell[target_cellIndex];
//     IndexType jj = indexShift + tid;
//     if (tid < target_numAtomInCell){
//       buff_globalIndex[tid] = globalIndex[jj];
//     }
//     __syncthreads();
    
//     if (tid < this_numAtom){
//       for (IndexType ll = 0; ll < my_numBond; ++ll){
// 	IndexType tofind_globalIndex =
// 	    global_neighborIndex[
// 		Parallel::DeviceBondList_cudaDevice::
// 		indexConvert(bondListStride, ii, ll)
// 		];
// 	for (IndexType mm = 0; mm < target_numAtomInCell; ++mm){
// 	  if (tofind_globalIndex == buff_globalIndex[mm]){
// 	    neighborIndex[
// 		Parallel::DeviceBondList_cudaDevice::
// 		indexConvert(bondListStride, ii, ll)
// 		] = mm + indexShift;
// 	  }
// 	  // __syncthreads();
// 	}
//       }
//     }
//   }
// }


  

  
