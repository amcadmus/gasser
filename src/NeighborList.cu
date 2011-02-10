#define DEVICE_CODE

#include "NeighborList.h"
#include<stdio.h>
#include "Auxiliary.h"
#include "NonBondedInteraction.h"

// extern texture<float, 1, cudaReadModeElementType> texRef;

__global__ void
prepare_naivelyBuildDeviceCellList (DeviceCellList clist)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  clist.data[bid * clist.stride + threadIdx.x] = MaxIndexValue;
  if (threadIdx.x == 0) clist.numbers[bid] = 0;
}



//////////////////////////////////////////////////
// coord in one vec
//////////////////////////////////////////////////

__global__ void
naivelyBuildDeviceCellList (const IndexType		numAtom,
			    const CoordType *		coord,
			    const RectangularBox	box,
			    DeviceCellList		clist,
			    mdError_t *			ptr_de,
			    IndexType *			erridx,
			    ScalorType *		errsrc)
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
    cellid = D3toD1 (clist.NCell, targetCelli, targetCellj, targetCellk);
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


__global__ void
naivelyBuildDeviceCellList2 (const IndexType		numAtom,
			     CoordType *		coord,
			     IntScalorType *		coordNoix,
			     IntScalorType *		coordNoiy,
			     IntScalorType *		coordNoiz,
			     const RectangularBox	box,
			     DeviceCellList		clist,
			     mdError_t *		ptr_de,
			     IndexType *		erridx,
			     ScalorType *		errsrc)
{
  // normalizeSystem (box, numAtom, coordx, coordy, coordz, coordNoix, coordNoiy, coordNoiz);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  
  // calculate target cell id
  extern __shared__ volatile IndexType sbuff[];
  volatile IndexType * targetCellid = (volatile IndexType *) sbuff;
  
  if (ii < numAtom){
    IndexType targetCelli, targetCellj, targetCellk;
    targetCelli = IndexType(coord[ii].x * box.sizei.x * ScalorType (clist.NCell.x));
    targetCellj = IndexType(coord[ii].y * box.sizei.y * ScalorType (clist.NCell.y));
    targetCellk = IndexType(coord[ii].z * box.sizei.z * ScalorType (clist.NCell.z));
    if (targetCelli == clist.NCell.x){
      targetCelli -= clist.NCell.x;
      coord[ii].x -= box.size.x;
      coordNoix[ii] ++;
    }
    if (targetCellj == clist.NCell.y){
      targetCellj -= clist.NCell.y;
      coord[ii].y -= box.size.y;
      coordNoiy[ii] ++;
    }
    if (targetCellk == clist.NCell.z){
      targetCellk -= clist.NCell.z;
      coord[ii].z -= box.size.z;
      coordNoiz[ii] ++;
    }
    targetCellid[tid] = D3toD1 (clist.NCell, targetCelli, targetCellj, targetCellk);
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
	IndexType pid = atomicInc(&clist.numbers[cellid], clist.stride);
  	clist.data[cellid * clist.stride + pid] = i + bid * blockDim.x;
	if (pid == clist.stride && ptr_de != NULL){
	  *ptr_de = mdErrorShortCellList;
	}
      }
      else
  	break;
    }
  }
}


__global__ void
buildCellNeighborhood (DeviceCellList clist,
		       const IndexType divide,
		       const HostVectorType boxSize)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  bool oneCellX(false), oneCellY(false), oneCellZ(false);
  if (clist.NCell.x == 1) oneCellX = true;
  if (clist.NCell.y == 1) oneCellY = true;
  if (clist.NCell.z == 1) oneCellZ = true;

  IndexType upperX, upperY, upperZ;
  oneCellX ? upperX = 1 : upperX = 2*divide+1;
  oneCellY ? upperY = 1 : upperY = 2*divide+1;
  oneCellZ ? upperZ = 1 : upperZ = 2*divide+1;
  
  
  if (tid == 0) {
    ScalorType rlist2 = clist.rlist * clist.rlist;
    clist.numNeighborCell[bid] = 0;
    int centerx, centery, centerz;
    D1toD3 (clist.NCell, int(bid), centerx, centery, centerz);
    for (int ix = 0; ix < upperX; ++ix){
      for (int iy = 0; iy < upperY; ++iy){
	for (int iz = 0; iz < upperZ; ++iz){
	  int myx = ix - clist.divide + int(centerx) ;
	  int myy = iy - clist.divide + int(centery) ;
	  int myz = iz - clist.divide + int(centerz) ;
	  ScalorType scalorx = boxSize.x / clist.NCell.x ;
	  ScalorType scalory = boxSize.y / clist.NCell.y ;
	  ScalorType scalorz = boxSize.z / clist.NCell.z ;
      
	  ScalorType min = 1e9;
#pragma unroll 27
	  for (int dx = -1; dx <= 1; ++dx){
	    for (int dy = -1; dy <= 1; ++dy){
	      for (int dz = -1; dz <= 1; ++dz){
		ScalorType diffx ((-centerx + myx + dx) * scalorx);
		ScalorType diffy ((-centery + myy + dy) * scalory);
		ScalorType diffz ((-centerz + myz + dz) * scalorz);
		// shortestImage (box, &diffx, &diffy, &diffz);
		ScalorType diff2 (diffx * diffx + diffy * diffy + diffz * diffz);
		if (diff2 < min){
		  min = diff2;
		}
	      }
	    }
	  }
	  CoordNoiType shift;
	  shift.x = shift.y = shift.z = 0;
	  if (min < rlist2){
	    if (myx < 0) {
	      myx += clist.NCell.x;
	      shift.x += 1;
	    }
	    else if (myx >= clist.NCell.x){
	      myx -= clist.NCell.x;
	      shift.x -= 1;
	    }
	    if (myy < 0) {
	      myy += clist.NCell.y;
	      shift.y += 1;
	    }
	    else if (myy >= clist.NCell.y){
	      myy -= clist.NCell.y;
	      shift.y -= 1;
	    }
	    if (myz < 0) {
	      myz += clist.NCell.z;
	      shift.z += 1;
	    }
	    else if (myz >= clist.NCell.z){
	      myz -= clist.NCell.z;
	      shift.z -= 1;
	    }

	    pushNeighborCell (clist,
			      bid,
			      D3toD1 (clist.NCell, myx, myy, myz),
			      shift);
	  }
	}
      }
    }
  }
}

  
  // if (tid == 0) clist.numNeighborCell[bid] = 0;
  // __syncthreads();
  
  // int centerx, centery, centerz;
  // D1toD3 (clist.NCell, int(bid), centerx, centery, centerz);
  // IntVectorType thisCubic;
  // thisCubic.z = thisCubic.y = thisCubic.x = 2 * clist.divide + 1;
  // int myx, myy, myz;
  // D1toD3 (thisCubic, int(tid), myx, myy, myz);
  // myx += int(centerx) - clist.divide;
  // myy += int(centery) - clist.divide ;
  // myz += int(centerz) - clist.divide ;
  // ScalorType scalorx = boxSize.x / clist.NCell.x / divide;
  // ScalorType scalory = boxSize.y / clist.NCell.y / divide;
  // ScalorType scalorz = boxSize.z / clist.NCell.z / divide;

  // ScalorType rlist2 = clist.rlist * clist.rlist;
  // ScalorType min = 1e9;
  // for (int dx = -1; dx <= 1; ++dx){
  //   for (int dy = -1; dy <= 1; ++dy){
  //     for (int dz = -1; dz <= 1; ++dz){
  // 	ScalorType diffx ((centerx - myx + dx) * scalorx);
  // 	ScalorType diffy ((centery - myy + dy) * scalory);
  // 	ScalorType diffz ((centerz - myz + dz) * scalorz);
  // 	ScalorType diff2 (diffx * diffx + diffy * diffy + diffz * diffz);
  // 	if (diff2 < min){
  // 	  min = diff2;
  // 	}
  //     }
  //   }
  // }
  // CoordType shift;
  // shift.x = shift.y = shift.z = 0.f;
  // if (min < rlist2){
  //   IndexType index = atomicInc (&clist.numNeighborCell[bid], blockDim.x);
  //   if (myx < 0) {
  //     myx += clist.NCell.x;
  //     shift.x -= boxSize.x;
  //   }
  //   else if (myx > clist.NCell.x){
  //     myx -= clist.NCell.x;
  //     shift.x += boxSize.x;
  //   }
  //   if (myy < 0) {
  //     myy += clist.NCell.y;
  //     shift.y -= boxSize.y;
  //   }
  //   else if (myy > clist.NCell.y){
  //     myy -= clist.NCell.y;
  //     shift.y += boxSize.y;
  //   }
  //   if (myz < 0) {
  //     myz += clist.NCell.z;
  //     shift.z -= boxSize.z;
  //   }
  //   else if (myz > clist.NCell.z){
  //     myz -= clist.NCell.z;
  //     shift.z += boxSize.z;
  //   }

  //   clist.neighborCellIndex[index] = D3toD1 (clist.NCell, myx, myy, myz);
  //   clist.neighborCellShift[index] = shift;
  // }
// }  



__global__ void
buildDeviceNeighborList_AllPair  (const IndexType		numAtom,
				  const CoordType *		coord,
				  const TypeType *		type,
				  const RectangularBox		box,
				  DeviceNeighborList		nlist,
				  const IndexType *		nbForceTable,
				  const IndexType		NatomType,
				  const bool			sharednbForceTable,
				  mdError_t *			ptr_de)
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
  IndexType * nbForceTableBuff = NULL;

  IndexType nbForceTableLength = AtomNBForceTable::dCalDataLength(NatomType);
  if (sharednbForceTable){
    nbForceTableBuff = (IndexType *) &targettype[roundUp4(blockDim.x)];
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
		    nbForceTable, NatomType, reftype, TypeType(targettype[kk]));
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


__global__ void
Reshuffle_backupDeviceNeighborList (const IndexType numAtom,
				    const IndexType * nlistData1,
				    const IndexType * nbForceIndex1,
				    const IndexType stride,
				    const IndexType * Nneighbor1,
				    IndexType * nlistData2,
				    IndexType * nbForceIndex2,
				    IndexType * Nneighbor2)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  // __shared__ volatile IndexType myNumbers [MaxThreadsPerBlock * 2];
  extern __shared__ volatile IndexType myNumbers[];
  
  IndexType N;
  if ((bid + 1) * blockDim.x < numAtom) N = blockDim.x;
  else if (bid * blockDim.x >= numAtom) N = 0;
  else N = numAtom - bid * blockDim.x;

  myNumbers[tid] = 0;
  myNumbers[tid + blockDim.x] = 0;
  if (ii < numAtom){
    Nneighbor2[ii] = myNumbers[tid] = Nneighbor1[ii];
  }
  __syncthreads();
  IndexType maxNum = maxVectorBlockBuffer (myNumbers, N);
  __syncthreads();

  for (IndexType jj = 0; jj < maxNum; ++jj){
    if (jj < myNumbers[tid] && ii < numAtom ){
      nlistData2   [jj * stride + ii] = nlistData1   [jj * stride + ii];
      nbForceIndex2[jj * stride + ii] = nbForceIndex1[jj * stride + ii];
    }
  }
}


__global__ void
Reshuffle_reshuffleDeviceCellList (IndexType * clistData,
				   const IndexType * idxTable)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  IndexType fromid = clistData[bid * blockDim.x + tid];
  if (fromid != MaxIndexValue){
    IndexType toid = idxTable[fromid];
    clistData[bid * blockDim.x + tid] = toid;
  }
}


__global__ void
Reshuffle_reshuffleDeviceNeighborList (const IndexType numAtom,
				       const IndexType * nlistData1,
				       const IndexType* nbForceIndex1,
				       const IndexType stride,
				       const IndexType * Nneighbor1,
				       const IndexType * idxTable,
				       IndexType * nlistData2,
				       IndexType * nbForceIndex2,
				       IndexType * Nneighbor2)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  // fromid is  ii

  extern __shared__ volatile IndexType myNumbers [];
  // __shared__ volatile IndexType myNumbers [MaxThreadsPerBlock * 2];

  IndexType toid;
  if (ii < numAtom)
    toid = idxTable[ii];
  IndexType myNum = 0;
  
  myNumbers[tid] = 0;
  myNumbers[tid + blockDim.x] = 0;
  if (ii < numAtom){
    myNum = myNumbers[tid] = Nneighbor1[ii];
  }
  __syncthreads();
  IndexType maxNum = maxVectorBlockBuffer (myNumbers, blockDim.x);
  __syncthreads();

  if (ii < numAtom){
    Nneighbor2[toid] = Nneighbor1[ii];
  }
  for (unsigned jj = 0; jj < maxNum; ++jj){
    if (jj < myNum && ii < numAtom){
      nlistData2[jj * stride + toid] = idxTable[nlistData1[jj * stride + ii]];
      nbForceIndex2[jj * stride + toid] = nbForceIndex1[jj * stride + ii];
    }
  }
}




