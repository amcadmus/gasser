#include "NeighborList.h"
#include<stdio.h>
#include "Auxiliary.h"
#include "NonBondedInteraction.h"

// extern texture<float, 1, cudaReadModeElementType> texRef;

__global__ void prepare_naivlyBuildDeviceCellList (DeviceCellList clist)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  clist.data[bid * clist.stride + threadIdx.x] = MaxIndexValue;
  if (threadIdx.x == 0) clist.numbers[bid] = 0;
}


#ifndef COORD_IN_ONE_VEC
__global__ void naivlyBuildDeviceCellList (IndexType numAtom,
					   ScalorType * coordx,
					   ScalorType * coordy,
					   ScalorType * coordz,
					   RectangularBox box,
					   DeviceCellList clist,
					   mdError_t * ptr_de,
					   IndexType * erridx,
					   ScalorType * errsrc)
{
  // normalizeSystem (box, numAtom, coordx, coordy, coordz, coordNoix, coordNoiy, coordNoiz);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  
  // calculate target cell id
  IndexType cellid ;
  if (ii < numAtom){
    IndexType targetCelli, targetCellj, targetCellk;
    targetCelli = IndexType(coordx[ii] * box.sizei.x * ScalorType (clist.NCell.x));
    targetCellj = IndexType(coordy[ii] * box.sizei.y * ScalorType (clist.NCell.y));
    targetCellk = IndexType(coordz[ii] * box.sizei.z * ScalorType (clist.NCell.z));
    cellid = D3toD1 (clist, targetCelli, targetCellj, targetCellk);
    if (ptr_de != NULL && 
	(targetCelli >= clist.NCell.x || 
	 targetCellj >= clist.NCell.y || 
	 targetCellk >= clist.NCell.z)){
      *ptr_de = mdErrorOverFlowCellIdx;
      if (targetCelli >= IndexType(clist.NCell.x)){
	*erridx = targetCelli;
	*errsrc = coordx[ii];
	return;
      }
      if (targetCellj >= IndexType(clist.NCell.y)){
	*erridx = targetCellj;
	*errsrc = coordy[ii];
	return;
      }
      if (targetCellk >= IndexType(clist.NCell.z)){
	*erridx = targetCellk;
	*errsrc = coordz[ii];
	return;
      }
    }
  }
  else {
    cellid = MaxIndexValue;
  }
  
  // write indexes to clist 
  if (cellid != MaxIndexValue){
    IndexType pid = atomicInc(&(clist.numbers[cellid]), blockDim.x);
    clist.data[cellid * clist.stride + pid]
	= ii;
    if (pid == blockDim.x && *ptr_de != NULL){
      *ptr_de = mdErrorShortCellList;
    }
  }
}


__global__ void naivlyBuildDeviceCellList2 (IndexType numAtom,
					    ScalorType * coordx,
					    ScalorType * coordy,
					    ScalorType * coordz,
					    RectangularBox box,
					    DeviceCellList clist,
					    mdError_t * ptr_de,
					    IndexType * erridx,
					    ScalorType * errsrc)
{
  // normalizeSystem (box, numAtom, coordx, coordy, coordz, coordNoix, coordNoiy, coordNoiz);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  
  // calculate target cell id
  __shared__ volatile IndexType targetCellid[MaxThreadsPerBlock];
  if (ii < numAtom){
    IndexType targetCelli, targetCellj, targetCellk;
    targetCelli = IndexType(coordx[ii] * box.sizei.x * ScalorType (clist.NCell.x));
    targetCellj = IndexType(coordy[ii] * box.sizei.y * ScalorType (clist.NCell.y));
    targetCellk = IndexType(coordz[ii] * box.sizei.z * ScalorType (clist.NCell.z));
    targetCellid[tid] = D3toD1 (clist, targetCelli, targetCellj, targetCellk);
    if (ptr_de != NULL && 
	(targetCelli >= clist.NCell.x || 
	 targetCellj >= clist.NCell.y || 
	 targetCellk >= clist.NCell.z)){
      *ptr_de = mdErrorOverFlowCellIdx;
      if (targetCelli >= IndexType(clist.NCell.x)){
	*erridx = targetCelli;
	*errsrc = coordx[ii];
	return;
      }
      if (targetCellj >= IndexType(clist.NCell.y)){
	*erridx = targetCellj;
	*errsrc = coordy[ii];
	return;
      }
      if (targetCellk >= IndexType(clist.NCell.z)){
	*erridx = targetCellk;
	*errsrc = coordz[ii];
	return;
      }      
    }
  }
  else {
    targetCellid[tid] = MaxIndexValue;
  }
  __syncthreads();
  // write indexes to clist only the first thread in the block to that
  if (threadIdx.x == 0){
    for (IndexType i = 0; i < blockDim.x; ++i){
      IndexType cellid = targetCellid[i];
      if (cellid != MaxIndexValue){
	IndexType pid = atomicInc(&clist.numbers[cellid], blockDim.x);
  	clist.data[cellid * clist.stride + pid] = i + bid * blockDim.x;
	if (pid == blockDim.x && ptr_de != NULL){
	  *ptr_de = mdErrorShortCellList;
	}
      }
      else break;
    }
  }
}



__global__ void buildDeviceNeighborList_AllPair  (IndexType numAtom,
						  ScalorType * coordx,
						  ScalorType * coordy, 
						  ScalorType * coordz,
						  TypeType * type,
						  RectangularBox box,
						  DeviceNeighborList nlist,
						  ForceIndexType * nbForceTable,
						  IndexType NatomType,
						  bool sharednbForceTable,
						  mdError_t * ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  
  IndexType Nneighbor = 0;
  IndexType numberAtom = numAtom;
  IndexType ii = tid + bid * blockDim.x;

  extern __shared__ volatile char pub_sbuff[];
  
  volatile IndexType * targetIndexes =
      (volatile IndexType *) pub_sbuff;
  volatile ScalorType * targetx =
      (volatile ScalorType *) &targetIndexes[roundUp4(blockDim.x)];
  volatile ScalorType * targety =
      (volatile ScalorType *) &targetx[roundUp4(blockDim.x)];
  volatile ScalorType * targetz =
      (volatile ScalorType *) &targety[roundUp4(blockDim.x)];
  volatile TypeType * targettype =
      (volatile TypeType *) &targetz[roundUp4(blockDim.x)];
  volatile ForceIndexType * nbForceTableBuff = NULL;

  IndexType nbForceTableLength = AtomNBForceTable::dCalDataLength(NatomType);
  if (sharednbForceTable){
    nbForceTableBuff = (volatile ForceIndexType *) &targettype[roundUp4(blockDim.x)];
    cpyGlobalDataToSharedBuff (nbForceTable, nbForceTableBuff, nbForceTableLength);
  }
  __syncthreads();

  // __shared__ volatile  ScalorType targetx    [MaxThreadsPerBlock];
  // __shared__ volatile  ScalorType targety    [MaxThreadsPerBlock];
  // __shared__ volatile  ScalorType targetz    [MaxThreadsPerBlock];
  // __shared__ volatile  TypeType   targettype [MaxThreadsPerBlock];
  // __shared__ volatile  ForceIndexType nbForceTableBuff [MaxNBForceTableBuffSize];

  ScalorType refx, refy, refz;
  TypeType reftype;
  if (ii < numberAtom){
    refx = coordx[ii];
    refy = coordy[ii];
    refz = coordz[ii];
    reftype = type[ii];
  }
  // int nbForceTag;
  // IndexType nbForceTableLength = AtomNBForceTable::dCalDataLength(NatomType);
  // if (nbForceTableLength < MaxNBForceTableBuffSize){
  //   cpyGlobalDataToSharedBuff (nbForceTable, nbForceTableBuff, nbForceTableLength);
  //   nbForceTag = 1;
  // }
  // else {
  //   nbForceTag = 0;
  // }

  for (IndexType targetBlockId = 0;
       targetBlockId * blockDim.x < numberAtom; ++targetBlockId){
    IndexType jj = tid + targetBlockId * blockDim.x;
    __syncthreads();
    if (jj < numberAtom){
      targetx[tid] = coordx[jj];
      targety[tid] = coordy[jj];
      targetz[tid] = coordz[jj];
      targettype[tid] = type[jj];
    }
    __syncthreads();
    if (ii < numberAtom){
      for (IndexType kk = 0; kk < blockDim.x; ++kk){
	if (kk + targetBlockId * blockDim.x >= numberAtom) break;
	ScalorType diffx = targetx[kk] - refx;
	ScalorType diffy = targety[kk] - refy;
	ScalorType diffz = targetz[kk] - refz;
	RectangularBoxGeometry::shortestImage (box, &diffx, &diffy, &diffz);
	if ((diffx*diffx+diffy*diffy+diffz*diffz) < nlist.rlist*nlist.rlist &&
	    kk + targetBlockId * blockDim.x != ii){
	  IndexType listIdx = Nneighbor * nlist.stride + ii;
	  nlist.data[listIdx] = kk + targetBlockId * blockDim.x;
	  if (sharednbForceTable){
	    nlist.forceIndex[listIdx] 
		= AtomNBForceTable::calForceIndex (
		    nbForceTableBuff, NatomType, reftype, targettype[kk]);
	  }
	  else {
	    nlist.forceIndex[listIdx] 
		= AtomNBForceTable::calForceIndex (
		    nbForceTable, NatomType, reftype, targettype[kk]);
	  }	
	  // if (nlist.forceIndex[listIdx] == 0){
	  //   printf ("%d  %d  reftype:%d targettype:%d\n",
	  // 	    ii, kk, reftype, targettype[kk]);
	  // }
	  Nneighbor ++;
	}
      }
    }
  }
  if (ii < numberAtom){
    nlist.Nneighbor[ii] = Nneighbor;
    if (Nneighbor > nlist.listLength && ptr_de != NULL){
      *ptr_de = mdErrorShortNeighborList;
    }
  }
}







//////////////////////////////////////////////////
// coord in one vec
//////////////////////////////////////////////////
#else
__global__ void naivlyBuildDeviceCellList (IndexType numAtom,
					   CoordType * coord,
					   RectangularBox box,
					   DeviceCellList clist,
					   mdError_t * ptr_de,
					   IndexType * erridx,
					   ScalorType * errsrc)
{
  // normalizeSystem (box, numAtom, coordx, coordy, coordz, coordNoix, coordNoiy, coordNoiz);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  
  // calculate target cell id
  IndexType cellid ;
  if (ii < numAtom){
    IndexType targetCelli, targetCellj, targetCellk;
    targetCelli = IndexType(coord[ii].x * box.sizei.x * ScalorType (clist.NCell.x));
    targetCellj = IndexType(coord[ii].y * box.sizei.y * ScalorType (clist.NCell.y));
    targetCellk = IndexType(coord[ii].z * box.sizei.z * ScalorType (clist.NCell.z));
    cellid = D3toD1 (clist, targetCelli, targetCellj, targetCellk);
    if (ptr_de != NULL && 
	(targetCelli >= clist.NCell.x || 
	 targetCellj >= clist.NCell.y || 
	 targetCellk >= clist.NCell.z)){
      *ptr_de = mdErrorOverFlowCellIdx;
      if (targetCelli >= IndexType(clist.NCell.x)){
	*erridx = targetCelli;
	*errsrc = coord[ii].x;
	return;
      }
      if (targetCellj >= IndexType(clist.NCell.y)){
	*erridx = targetCellj;
	*errsrc = coord[ii].y;
	return;
      }
      if (targetCellk >= IndexType(clist.NCell.z)){
	*erridx = targetCellk;
	*errsrc = coord[ii].z;
	return;
      }
    }
  }
  else {
    cellid = MaxIndexValue;
  }
  
  // write indexes to clist 
  if (cellid != MaxIndexValue){
    IndexType pid = atomicInc(&(clist.numbers[cellid]), blockDim.x);
    clist.data[cellid * clist.stride + pid]
	= ii;
    if (pid == blockDim.x && *ptr_de != NULL){
      *ptr_de = mdErrorShortCellList;
    }
  }
}


__global__ void naivlyBuildDeviceCellList2 (IndexType numAtom,
					    CoordType * coord,
					    RectangularBox box,
					    DeviceCellList clist,
					    mdError_t * ptr_de,
					    IndexType * erridx,
					    ScalorType * errsrc)
{
  // normalizeSystem (box, numAtom, coordx, coordy, coordz, coordNoix, coordNoiy, coordNoiz);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  
  // calculate target cell id
  __shared__ volatile IndexType targetCellid[MaxThreadsPerBlock];
  if (ii < numAtom){
    IndexType targetCelli, targetCellj, targetCellk;
    targetCelli = IndexType(coord[ii].x * box.sizei.x * ScalorType (clist.NCell.x));
    targetCellj = IndexType(coord[ii].y * box.sizei.y * ScalorType (clist.NCell.y));
    targetCellk = IndexType(coord[ii].z * box.sizei.z * ScalorType (clist.NCell.z));
    targetCellid[tid] = D3toD1 (clist, targetCelli, targetCellj, targetCellk);
    if (ptr_de != NULL && 
	(targetCelli >= clist.NCell.x || 
	 targetCellj >= clist.NCell.y || 
	 targetCellk >= clist.NCell.z)){
      *ptr_de = mdErrorOverFlowCellIdx;
      if (targetCelli >= IndexType(clist.NCell.x)){
	*erridx = targetCelli;
	*errsrc = coord[ii].x;
	return;
      }
      if (targetCellj >= IndexType(clist.NCell.y)){
	*erridx = targetCellj;
	*errsrc = coord[ii].y;
	return;
      }
      if (targetCellk >= IndexType(clist.NCell.z)){
	*erridx = targetCellk;
	*errsrc = coord[ii].z;
	return;
      }      
    }
  }
  else {
    targetCellid[tid] = MaxIndexValue;
  }
  __syncthreads();
  // write indexes to clist only the first thread in the block to that
  if (threadIdx.x == 0){
    for (IndexType i = 0; i < blockDim.x; ++i){
      IndexType cellid = targetCellid[i];
      if (cellid != MaxIndexValue){
	IndexType pid = atomicInc(&clist.numbers[cellid], blockDim.x);
  	clist.data[cellid * clist.stride + pid] = i + bid * blockDim.x;
	if (pid == blockDim.x && ptr_de != NULL){
	  *ptr_de = mdErrorShortCellList;
	}
      }
      else
  	break;
    }
  }
}



__global__ void buildDeviceNeighborList_AllPair  (IndexType numAtom,
						  CoordType * coord,
						  TypeType * type,
						  RectangularBox box,
						  DeviceNeighborList nlist,
						  ForceIndexType * nbForceTable,
						  IndexType NatomType,
						  bool sharednbForceTable,
						  mdError_t * ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  
  IndexType Nneighbor = 0;
  IndexType numberAtom = numAtom;
  IndexType ii = tid + bid * blockDim.x;

  extern __shared__ volatile char pub_sbuff[];

  volatile IndexType * targetIndexes =
      (volatile IndexType *) pub_sbuff;
  volatile CoordType * target =
      (volatile CoordType *) &targetIndexes[roundUp4(blockDim.x)];
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];
  volatile ForceIndexType * nbForceTableBuff = NULL;

  IndexType nbForceTableLength = AtomNBForceTable::dCalDataLength(NatomType);
  if (sharednbForceTable){
    nbForceTableBuff = (volatile ForceIndexType *) &targettype[roundUp4(blockDim.x)];
    cpyGlobalDataToSharedBuff (nbForceTable, nbForceTableBuff, nbForceTableLength);
  }
  __syncthreads();
  
  CoordType ref;
  TypeType reftype;
  if (ii < numberAtom){
    ref = coord[ii];
    reftype = type[ii];
  }

  for (IndexType targetBlockId = 0;
       targetBlockId * blockDim.x < numberAtom; ++targetBlockId){
    IndexType jj = tid + targetBlockId * blockDim.x;
    __syncthreads();
    if (jj < numberAtom){
      target[tid].x = coord[jj].x;
      target[tid].y = coord[jj].y;
      target[tid].z = coord[jj].z;
      targettype[tid] = type[jj];
    }
    __syncthreads();
    if (ii < numberAtom){
      for (IndexType kk = 0; kk < blockDim.x; ++kk){
	if (kk + targetBlockId * blockDim.x >= numberAtom) break;
	ScalorType diffx = target[kk].x - ref.x;
	ScalorType diffy = target[kk].y - ref.y;
	ScalorType diffz = target[kk].z - ref.z;
	RectangularBoxGeometry::shortestImage (box, &diffx, &diffy, &diffz);
	if ((diffx*diffx+diffy*diffy+diffz*diffz) < nlist.rlist*nlist.rlist &&
	    kk + targetBlockId * blockDim.x != ii){
	  IndexType listIdx = Nneighbor * nlist.stride + ii;
	  nlist.data[listIdx] = kk + targetBlockId * blockDim.x;
	  if (sharednbForceTable){
	    nlist.forceIndex[listIdx] 
		= AtomNBForceTable::calForceIndex (
		    nbForceTableBuff, NatomType, reftype, targettype[kk]);
	  }
	  else {
	    nlist.forceIndex[listIdx] 
		= AtomNBForceTable::calForceIndex (
		    nbForceTable, NatomType, reftype, targettype[kk]);
	  }	
	  // if (nlist.forceIndex[listIdx] == 0){
	  //   printf ("%d  %d  reftype:%d targettype:%d\n",
	  // 	    ii, kk, reftype, targettype[kk]);
	  // }
	  Nneighbor ++;
	}
      }
    }
  }
  if (ii < numberAtom){
    nlist.Nneighbor[ii] = Nneighbor;
    if (Nneighbor > nlist.listLength && ptr_de != NULL){
      *ptr_de = mdErrorShortNeighborList;
    }
  }
}

#endif




