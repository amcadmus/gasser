/**
 * @file   NeighborList_interface.cu
 * @author Han Wang <han_wang@math.pku.edu.cn>
 * @date   Thu Nov 19 12:53:42 2009
 * 
 * @brief  Implementation of neighbor list
 * 
 * 
 */

#define DEVICE_CODE

#include "NeighborList_interface.h"
#include "Auxiliary.h"
#include "NeighborList.h"
#include <stdio.h>
#include "NonBondedInteraction.h"
#include "Reshuffle_interface.h"

/** 
 * these are textures for a fast reference of particle position.
 * 
 */

texture<CoordType, 1, cudaReadModeElementType> global_texRef_neighbor_coord;
texture<TypeType, 1, cudaReadModeElementType> global_texRef_neighbor_type;


void NeighborList::
clearDeviceNeighborList()
{
  if ( mallocedDeviceNeighborList ){
    cudaFree (dnlist.data);
    cudaFree (dnlist.Nneighbor);
    cudaFree (dnlist.forceIndex);
    mallocedDeviceNeighborList = false;
    checkCUDAError ("NeighborList::clearDeviceNeighborList");
  }
}

void NeighborList::
clearNonBondedForce ()
{
  if (mallocedNonBondedForceTable == true){
    cudaFree (nbForceTable);
    mallocedNonBondedForceTable = false;
  }
}

void NeighborList::
clear()
{
  clearDeviceNeighborList();
  clearNonBondedForce();
  unbindGlobalTexture ();
}  

void NeighborList::
unbindGlobalTexture ()
{
  if ( initedGlobalTexture ){
    cudaUnbindTexture(global_texRef_neighbor_coord);
    cudaUnbindTexture(global_texRef_neighbor_type);
    initedGlobalTexture = false;
    checkCUDAError ("NeighborList::unbindGlobalTexture");
  }
}

void NeighborList::
bindGlobalTexture (const MDSystem & sys)
{
  size_t sizetype   = sizeof(TypeType)  *sys.ddata.numMem;
  size_t sizecoord  = sizeof(CoordType) *sys.ddata.numMem;
  cudaBindTexture(0, global_texRef_neighbor_coord, sys.ddata.coord, sizecoord);
  cudaBindTexture(0, global_texRef_neighbor_type, sys.ddata.type, sizetype);
  checkCUDAError ("NeighborList::init texture");
  initedGlobalTexture = true;
}


NeighborList::~NeighborList()
{
  clear();
}

static IndexType hroundUp4 (IndexType x)
{
  if (x & 3 == 0){
    return x;
  }
  else {
    return ((x >> 2) + 1) << 2;
  }
}


void NeighborList::
buildDeviceNeighborListCellList (const MDSystem & sys,
				 const CellList & clist)
{
  dim3 cellBlockDim = clist.getCellBlockDim();
  bool sharednbForceTable (true);
  size_t buildDeviceNeighborList_DeviceCellList_sbuffSize =
      sizeof(IndexType) *	hroundUp4(cellBlockDim.x) +
      sizeof(CoordType) *	hroundUp4(cellBlockDim.x) +
      sizeof(TypeType)  *	hroundUp4(cellBlockDim.x) +
      sizeof(IndexType) *	hroundUp4(nbForceTableLength);
  if (buildDeviceNeighborList_DeviceCellList_sbuffSize >=
      SystemSharedBuffSize - GlobalFunctionParamSizeLimit){
    sharednbForceTable = false;
    buildDeviceNeighborList_DeviceCellList_sbuffSize =
	sizeof(IndexType) *	hroundUp4(cellBlockDim.x) +
	sizeof(CoordType) *	hroundUp4(cellBlockDim.x) +
	sizeof(TypeType)  *	hroundUp4(cellBlockDim.x);
  }
  buildDeviceNeighborList_DeviceCellList 
      <<<clist.getCellGrimDim(), cellBlockDim, 
      buildDeviceNeighborList_DeviceCellList_sbuffSize>>> (
	  // <<<cellGridDim, myBlockDim>>> (
	  sys.ddata.numAtom, 
	  sys.ddata.coord,
	  sys.ddata.type,
	  sys.box,
	  clist.dclist,
	  dnlist,
	  nbForceTable,
	  NatomType,
	  sharednbForceTable,
	  err.ptr_de);
  err.check("NeighborList::buildDeviceNeighborListCellList");
  checkCUDAError ("NeighborList::buildDeviceNeighborListCellList");
}

void NeighborList::
buildDeviceNeighborListAllPair (const MDSystem & sys)
{
  bool sharednbForceTable (true);
  size_t buildDeviceNeighborList_AllPair_sbuffSize =
      sizeof(IndexType) *	hroundUp4(myBlockDim.x) +
      sizeof(CoordType) *	hroundUp4(myBlockDim.x) +
      sizeof(TypeType)  *	hroundUp4(myBlockDim.x) +
      sizeof(IndexType) *	hroundUp4(nbForceTableLength);
  if (buildDeviceNeighborList_AllPair_sbuffSize >=
      SystemSharedBuffSize - GlobalFunctionParamSizeLimit){
    sharednbForceTable = false;
    buildDeviceNeighborList_AllPair_sbuffSize =
	sizeof(IndexType) *	hroundUp4(myBlockDim.x) +
	sizeof(CoordType) *	hroundUp4(myBlockDim.x) +
	sizeof(TypeType)  *	hroundUp4(myBlockDim.x);
  }
  buildDeviceNeighborList_AllPair 
      <<<atomGridDim, myBlockDim,
      buildDeviceNeighborList_AllPair_sbuffSize>>>(
	  sys.ddata.numAtom,
	  sys.ddata.coord,
	  sys.ddata.type,
	  sys.box,
	  dnlist,
	  nbForceTable,
	  NatomType,
	  sharednbForceTable,
	  err.ptr_de);
  err.check("NeighborList::build, build neighbor list all pair");
  checkCUDAError ("NeighborList::build, build neighbor list all pair");
}


void NeighborList::
initNonBondedInteraction (const SystemNonBondedInteraction & sysNbInter)
{
  if (! sysNbInter.beBuilt()) {
    throw MDExcptUnbuiltNonBondedInteraction ("NeighborList");
  }
  NatomType = sysNbInter.numberOfAtomTypes();
  nbForceTableLength = sysNbInter.interactionTableSize();
  cudaMalloc ((void**)&nbForceTable,
	      nbForceTableLength * sizeof(IndexType));
  cudaMemcpy (nbForceTable,
	      sysNbInter.interactionTable(),
	      nbForceTableLength * sizeof(IndexType),
	      cudaMemcpyHostToDevice);
  checkCUDAError ("AtomNBForceTable::deviceInitTable");
  mallocedNonBondedForceTable = true;
}

void NeighborList::
mallocDeviceNeighborList (const MDSystem & sys,
			  const ScalorType & DeviceNeighborListExpansion)
{
  ScalorType density = sys.ddata.numAtom / (sys.box.size.x * sys.box.size.y * sys.box.size.z);
  ScalorType expectedNumberInList 
      = 4./3. * M_PI * myrlist * myrlist * myrlist * density;
  dnlist.listLength = IndexType(expectedNumberInList * DeviceNeighborListExpansion);
  if (dnlist.listLength < 30){
    dnlist.listLength = 30;
  }
  cudaMalloc ((void**)&(dnlist.data), sizeof(IndexType) * dnlist.stride * dnlist.listLength);
  cudaMalloc ((void**)&(dnlist.Nneighbor), sizeof(IndexType) * sys.ddata.numAtom);
  cudaMalloc ((void**)&(dnlist.forceIndex), sizeof(IndexType) *  dnlist.stride * dnlist.listLength);
  // reshuffle backup things
  cudaMalloc ((void**)&(bkdnlistData), sizeof(IndexType) * dnlist.stride * dnlist.listLength);
  cudaMalloc ((void**)&(bkdnlistNneighbor), sizeof(IndexType) * sys.ddata.numAtom);
  cudaMalloc ((void**)&(bkdnlistForceIndex), sizeof(IndexType) *  dnlist.stride * dnlist.listLength);
  checkCUDAError ("NeighborList::mallocDeviceNeighborList");
  mallocedDeviceNeighborList = true;
}



void NeighborList::
reinit (const SystemNonBondedInteraction & sysNbInter,
	const MDSystem & sys,
	const ScalorType & rlist,
	const IndexType & NTread,
	const ScalorType & DeviceNeighborListExpansion)
{
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NTread;

  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);
  
  myrlist = rlist;
  dnlist.rlist = myrlist;
  dnlist.stride = sys.ddata.numAtom;
  
  // init neighbor list
  clearDeviceNeighborList ();
  mallocDeviceNeighborList (sys, DeviceNeighborListExpansion);  

  clearNonBondedForce ();
  initNonBondedInteraction (sysNbInter);

  unbindGlobalTexture ();
  bindGlobalTexture (sys);
  
  //init shared memory size
}

void NeighborList::
rebuild (const MDSystem & sys,
	 const CellList & clist,
	 MDTimer * timer)
{
  if (clist.isempty()){
    if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
    // printf ("rlist is %f\n", dnlist.rlist);
    buildDeviceNeighborListAllPair (sys);
    if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
  }
  else {
    if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
    buildDeviceNeighborListCellList (sys, clist);
    if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
  }  
}


void NeighborList::
reshuffle (const IndexType * indexTable,
	   const IndexType & numAtom,
	   MDTimer *timer)
{
  if (timer != NULL) timer->tic(mdTimeReshuffleSystem);

  // Reshuffle_reshuffleDeviceCellList
  //     <<<cellGridDim, myBlockDim>>> (
  // 	  dclist.data, indexTable);
  // cudaMemcpy (bkbackupCoord, backupCoord,
  // 	      sizeof (CoordType) * numAtom,
  // 	      cudaMemcpyDeviceToDevice);
  // Reshuffle_reshuffleArray
  //     <<<atomGridDim, myBlockDim>>> 
  //     (bkbackupCoord, numAtom, indexTable, backupCoord);
  Reshuffle_backupDeviceNeighborList
      <<<atomGridDim, myBlockDim,
      2 * myBlockDim.x * sizeof(IndexType)>>> (
	  numAtom,
	  dnlist.data,
	  dnlist.forceIndex,
	  dnlist.stride,
	  dnlist.Nneighbor,
	  bkdnlistData,
	  bkdnlistForceIndex,
	  bkdnlistNneighbor);
  checkCUDAError ("NeighborList::reshuffle backup");
  Reshuffle_reshuffleDeviceNeighborList
      <<<atomGridDim, myBlockDim,
      2 * myBlockDim.x * sizeof(IndexType)>>> (
	  numAtom,
	  bkdnlistData,
	  bkdnlistForceIndex,
	  dnlist.stride,
	  bkdnlistNneighbor,
	  indexTable,
	  dnlist.data,
	  dnlist.forceIndex,
	  dnlist.Nneighbor);
  checkCUDAError ("NeighborList::reshuffle reshuffle");
  if (timer != NULL) timer->toc(mdTimeReshuffleSystem);
}

NeighborList::
NeighborList (const SystemNonBondedInteraction & sysNbInter,
	      const MDSystem & sys,
	      const ScalorType & rlist,
	      const IndexType & NTread,
	      const ScalorType & DeviceNeighborListExpansion)
    : mallocedDeviceNeighborList (false),
      mallocedNonBondedForceTable (false),
      initedGlobalTexture (false)
{
  reinit (sysNbInter, sys, rlist, NTread, DeviceNeighborListExpansion);
}





////////////////////////////////////////////////////////////
// for the reason of using texture, we place this function here. it
// should be placed in NeighborList.cu
////////////////////////////////////////////////////////////

using namespace RectangularBoxGeometry;

__device__ IndexType
shiftedD3toD1 (DeviceCellList clist,
	       RectangularBox box,
	       int ix,
	       int iy,
	       int iz,
	       ScalorType * shiftx ,
	       ScalorType * shifty,
	       ScalorType * shiftz)
{
  int tmp;
  ix += (tmp = -int(floorf(ix * clist.NCelli.x))) * clist.NCell.x;
  *shiftx = tmp * box.size.x;
  iy += (tmp = -int(floorf(iy * clist.NCelli.y))) * clist.NCell.y;
  *shifty = tmp * box.size.y;
  iz += (tmp = -int(floorf(iz * clist.NCelli.z))) * clist.NCell.z;
  *shiftz = tmp * box.size.z;
  return D3toD1 (clist.NCell, ix, iy, iz);
}

__global__ void
buildDeviceNeighborList_DeviceCellList (const IndexType		numAtom,
					const CoordType *	coord,
					const TypeType *	type,
					const RectangularBox	box,
					const DeviceCellList	clist,
					DeviceNeighborList	nlist,
					const IndexType *	nbForceTable,
					const IndexType		NatomType,
					const bool		sharednbForceTable,
					mdError_t *		ptr_de )
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType bidx, bidy, bidz;
  D1toD3 (clist.NCell, bid, bidx, bidy, bidz);
  
  // set number of neighbor to 0
  IndexType Nneighbor = 0;
  // load index
  IndexType ii = getDeviceCellListData (clist, bid, tid);
  // load iith coordinate // use texturefetch instead
  CoordType ref;
  TypeType reftype;
  if (ii != MaxIndexValue){
#ifdef COMPILE_NO_TEX
    ref = coord[ii];
    reftype = type[ii];
#else
    ref = tex1Dfetch (global_texRef_neighbor_coord, ii);
    reftype = tex1Dfetch(global_texRef_neighbor_type, ii);
#endif
  }
  ScalorType rlist = nlist.rlist;

  // the target index and coordinates are shared

  extern __shared__ volatile char pub_sbuff[];
  
  volatile IndexType * targetIndexes =
      (volatile IndexType *) pub_sbuff;
  CoordType * target =
      (CoordType *) &targetIndexes[roundUp4(blockDim.x)];
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];
  IndexType * nbForceTableBuff = NULL;

  IndexType nbForceTableLength = AtomNBForceTable::dCalDataLength(NatomType);
  if (sharednbForceTable){
    nbForceTableBuff = (IndexType *) &targettype[roundUp4(blockDim.x)];
    cpyGlobalDataToSharedBuff (nbForceTable, nbForceTableBuff, nbForceTableLength);
  }
  __syncthreads();

  // __shared__ volatile  IndexType  targetIndexes [MaxThreadsPerBlock];
  // __shared__ volatile CoordType target [MaxThreadsPerBlock];
  // __shared__ volatile  TypeType   targettype    [MaxThreadsPerBlock];
  // __shared__ volatile  IndexType nbForceTableBuff [MaxNBForceTableBuffSize];

  // IndexType nbForceTableLength = AtomNBForceTable::dCalDataLength(NatomType);
  // if (sharednbForceTable){
  //   cpyGlobalDataToSharedBuff (nbForceTable, nbForceTableBuff, nbForceTableLength);
  // }
  // __syncthreads();

  bool oneCellX(false), oneCellY(false), oneCellZ(false);
  if (clist.NCell.x == 1) oneCellX = true;
  if (clist.NCell.y == 1) oneCellY = true;
  if (clist.NCell.z == 1) oneCellZ = true;
  ScalorType rlist2 = rlist * rlist;
  
  // loop over 27 neighbor cells
  for (IndexType i = 0; i < clist.numNeighborCell[bid]; ++i){
    // if (threadIdx.x == 0){
    //   printf ("%d %d\n", bid, clist.numNeighborCell[bid]);
    // }
    __syncthreads();
    IndexType targetCellIdx = getNeighborCellIndex    (clist, bid, i);
    CoordNoiType shiftNoi   = getNeighborCellShiftNoi (clist, bid, i);
    CoordType shift;
    shift.x = shiftNoi.x * box.size.x;
    shift.y = shiftNoi.y * box.size.y;
    shift.z = shiftNoi.z * box.size.z;
    targetIndexes[tid] = getDeviceCellListData(clist, targetCellIdx, tid);  
    if (targetIndexes[tid] != MaxIndexValue){
      target[tid] = tex1Dfetch(global_texRef_neighbor_coord, targetIndexes[tid]);
      targettype[tid] = tex1Dfetch(global_texRef_neighbor_type, targetIndexes[tid]);
    }
    __syncthreads();
	
    // find neighbor
    if (ii != MaxIndexValue){
      for (IndexType jj = 0; jj < clist.numbers[targetCellIdx]; ++jj){
	ScalorType diffx = target[jj].x - shift.x - ref.x;
	ScalorType diffy = target[jj].y - shift.y - ref.y;
	ScalorType diffz = target[jj].z - shift.z - ref.z;
	if (oneCellX) shortestImage (box.size.x, box.sizei.x, &diffx);
	if (oneCellY) shortestImage (box.size.y, box.sizei.y, &diffy);
	if (oneCellZ) shortestImage (box.size.z, box.sizei.z, &diffz);
	if ((diffx*diffx+diffy*diffy+diffz*diffz) < rlist2 &&
	    targetIndexes[jj] != ii){
	  IndexType fidx;
	  if (sharednbForceTable){
	    fidx = AtomNBForceTable::calForceIndex (
		nbForceTableBuff, NatomType, reftype, targettype[jj]);
	  }
	  else {
	    fidx = AtomNBForceTable::calForceIndex (
		nbForceTable, NatomType, reftype, targettype[jj]);
	  }
	  // if (fidx != mdForceNULL) {
	  IndexType listIdx = Nneighbor * nlist.stride + ii;
	  nlist.data[listIdx] = targetIndexes[jj];
	  nlist.forceIndex[listIdx] = fidx;
	  Nneighbor ++;
	  // }
	}
      }
    }
  }

  if (ii != MaxIndexValue) {
    if (Nneighbor > nlist.listLength && ptr_de != NULL){
      *ptr_de = mdErrorShortNeighborList;
      return;
    }
    nlist.Nneighbor[ii] = Nneighbor;
    // printf ("%d %d\n", ii, Nneighbor);
  }
}



__global__ void
buildDeviceCellList_step1 (IndexType		numAtom,
			   CoordType *		coord,
			   IntScalorType *	coordNoix,
			   IntScalorType *	coordNoiy,
			   IntScalorType *	coordNoiz,
			   RectangularBox	box,
			   DeviceCellList	clist,
			   IndexType *		sendBuff,
			   IndexType *		targetBuff,
			   mdError_t *		ptr_de,
			   IndexType *		erridx,
			   ScalorType *		errsrc)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  extern __shared__ volatile IndexType sbuff[];
  volatile IndexType * originalData = (volatile IndexType *) sbuff;
  volatile IndexType * targetCellid = (volatile IndexType *) &originalData[blockDim.x];
  
  // __shared__ volatile  IndexType originalData[MaxThreadsPerBlock];
  // __shared__ volatile  IndexType targetCellid[MaxThreadsPerBlock];
  
  // copy data from cell list
  originalData[tid] = clist.data[bid*clist.stride + tid];
  IndexType originalNumber = clist.numbers[bid];

  // calculate the target cell
  if (originalData[tid] != MaxIndexValue){
    IndexType targetCelli, targetCellj, targetCellk;
    IndexType thisid = originalData[tid];
#ifdef COMPILE_NO_TEX
    ref = coord[thisid];
#else
    CoordType ref (tex1Dfetch(global_texRef_neighbor_coord, thisid));
#endif
    targetCelli = IndexType(ref.x * box.sizei.x * ScalorType (clist.NCell.x));
    targetCellj = IndexType(ref.y * box.sizei.y * ScalorType (clist.NCell.y));
    targetCellk = IndexType(ref.z * box.sizei.z * ScalorType (clist.NCell.z));
    if (targetCelli == clist.NCell.x){
      targetCelli -= clist.NCell.x;
      coord[thisid].x -= box.size.x;
      coordNoix[thisid] ++;
    }
    if (targetCellj == clist.NCell.y){
      targetCellj -= clist.NCell.y;
      coord[thisid].y -= box.size.y;
      coordNoiy[thisid] ++;
    }
    if (targetCellk == clist.NCell.z){
      targetCellk -= clist.NCell.z;
      coord[thisid].z -= box.size.z;
      coordNoiz[thisid] ++;
    }
    targetCellid[tid] = D3toD1 (clist.NCell, targetCelli, targetCellj, targetCellk);
    if (ptr_de != NULL && 
	(targetCelli >= clist.NCell.x || 
	 targetCellj >= clist.NCell.y || 
	 targetCellk >= clist.NCell.z)){
      *ptr_de = mdErrorOverFlowCellIdx;
      if (targetCelli >= IndexType(clist.NCell.x)){
	*erridx = targetCelli;
	*errsrc = ref.x;
	// return;
      }
      if (targetCellj >= IndexType(clist.NCell.y)){
	*erridx = targetCellj;
	*errsrc = ref.y;
	// return;
      }
      if (targetCellk >= IndexType(clist.NCell.z)){
	*erridx = targetCellk;
	*errsrc = ref.z;
	// return;
      }      
    }
  }
  else {
    targetCellid[tid] = MaxIndexValue;
  }

  // mark particles to be send 
  IndexType mark = MaxIndexValue - (MaxIndexValue >> 1);
  if (tid < originalNumber && targetCellid[tid] != bid){
    originalData[tid] += mark;
  }
  
  // head sort
  IndexType total1 = headSort (originalData, targetCellid);
  IndexType total0 = blockDim.x - total1;
  
  // unmark and copy to send buff
  if (tid < originalNumber && targetCellid[tid] != bid){
    sendBuff  [bid*clist.stride + tid - total0] = originalData[tid] - mark;
    targetBuff[bid*clist.stride + tid - total0] = targetCellid[tid];
    originalData[tid] = MaxIndexValue;
  }
  __syncthreads();
  
  // modify cell list
  clist.data[bid*clist.stride + tid] = originalData[tid];
  if (tid == 0) clist.numbers[bid] = total0;
}



__global__ void
buildDeviceCellList_initBuff (IndexType * sendBuff,
			      IndexType * targetBuff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  sendBuff[ii] = MaxIndexValue;
  targetBuff[ii] = 0;
}

__global__ void
buildDeviceCellList_clearBuff (IndexType * sendBuff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  sendBuff[ii] = MaxIndexValue;
}

  
  
__global__ void
buildDeviceCellList_step2 (RectangularBox	box,
			   DeviceCellList	clist,
			   IndexType *		sendBuff,
			   IndexType *		targetBuff,
			   IndexType		bitDeepth,
			   mdError_t *		ptr_de)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  IndexType thisid;
  IndexType ii = 0;
  IndexType buffPosi;
  if (tid == 0){
    while ((thisid = sendBuff[buffPosi = (bid*clist.stride + ii)]) != MaxIndexValue){
      IndexType cellid = targetBuff[buffPosi];
      IndexType tailIdx = atomicInc(&clist.numbers[cellid], blockDim.x);
      if (tailIdx >= blockDim.x&& ptr_de != NULL) {
	*ptr_de = mdErrorShortCellList;
	return;
      } 
      clist.data[cellid * clist.stride + tailIdx] = thisid;
      sendBuff[buffPosi] = MaxIndexValue;
      ii ++;
    }
  }
}










// void NeighborList::
// DecideNeighboringMethod (const MDSystem & sys,
// 			 const ScalorType & rlist,
// 			 const BoxDirection_t & bdir,
// 			 const IndexType & divide,
// 			 NeighborListBuiltMode & mode,
// 			 IntVectorType & NCell)
// {
//   bool CellOnX, CellOnY, CellOnZ;
//   CellOnX = bdir & RectangularBoxGeometry::mdRectBoxDirectionX;
//   CellOnY = bdir & RectangularBoxGeometry::mdRectBoxDirectionY;
//   CellOnZ = bdir & RectangularBoxGeometry::mdRectBoxDirectionZ;
//   double rlisti = 1./rlist;
//   if (CellOnX ) NCell.x = int ( floor(sys.box.size.x * rlisti) );
//   else NCell.x = 1;
//   if (CellOnY ) NCell.y = int ( floor(sys.box.size.y * rlisti) );
//   else NCell.y = 1;
//   if (CellOnZ ) NCell.z = int ( floor(sys.box.size.z * rlisti) );
//   else NCell.z = 1;
  
//   mode = CellListBuilt;
//   if (CellOnX && dclist.NCell.x < 4) mode = AllPairBuilt;
//   if (CellOnY && dclist.NCell.y < 4) mode = AllPairBuilt;
//   if (CellOnZ && dclist.NCell.z < 4) mode = AllPairBuilt;

//   if (mode == CellListBuilt){
//     if (CellOnX) NCell.x *= divide;
//     if (CellOnY) NCell.y *= divide;
//     if (CellOnZ) NCell.z *= divide;
//   } 
// }

// void NeighborList::
// mallocDeviceCellList (const IntVectorType & NCell,
// 		      const HostVectorType & boxSize)
// {
//   if (NCell.x == 1){
//     dclist.NCell.x = 1;
//     dclist.NCelli.x = 1./boxSize.x;
//   }
//   else{
//     dclist.NCelli.x = 1./ (dclist.NCell.x = NCell.x);
//   }
//   if (NCell.y == 1){
//     dclist.NCell.y = 1;
//     dclist.NCelli.y = 1./boxSize.y;
//   }
//   else{
//     dclist.NCelli.y = 1./ (dclist.NCell.y = NCell.y);
//   }
//   if (NCell.z == 1){
//     dclist.NCell.z = 1;
//     dclist.NCelli.z = 1./boxSize.z;
//   }
//   else{
//     dclist.NCelli.z = 1./ (dclist.NCell.z = NCell.z);
//   }

//   IndexType numCell = dclist.NCell.x * dclist.NCell.y * dclist.NCell.z;
//   dclist.rlist = myrlist;
//   dclist.divide = mydivide;
//   // suppose the number of atoms in any cell is smaller or equal
//   // to the number of threads in a block
//   dclist.stride = myBlockDim.x;
//   cellGridDim = toGridDim (numCell);

//   cudaMalloc ((void**)&(dclist.data), 
// 	      sizeof(ScalorType) * numCell * dclist.stride);
//   cudaMalloc ((void**)&(dclist.numbers), sizeof(unsigned) * numCell);
//   cudaMalloc ((void**)&mySendBuff, 
// 	      sizeof(ScalorType) * numCell * dclist.stride);
//   cudaMalloc ((void**)&myTargetBuff, 
// 	      sizeof(ScalorType) * numCell * dclist.stride);
//   checkCUDAError ("NeighborList::init cell list");

//   IndexType maxNumNeighborCell = (2*mydivide+1) * (2*mydivide+1) * (2*mydivide+1);
//   dclist.maxNumNeighborCell = maxNumNeighborCell;
//   cudaMalloc ((void**)&(dclist.numNeighborCell),
// 	      sizeof(IndexType) * numCell);
//   cudaMalloc ((void**)&(dclist.neighborCellIndex),
// 	      sizeof(IndexType) * maxNumNeighborCell * numCell);
//   cudaMalloc ((void**)&(dclist.neighborCellShiftNoi),
// 	      sizeof(CoordNoiType) * maxNumNeighborCell * numCell);
//   checkCUDAError ("NeighborList::maxNumNeighborCell cell list buff");

//   mallocedDeviceCellList = true;

//   buildDeviceCellList_initBuff<<<cellGridDim, myBlockDim>>> 
//       (mySendBuff, myTargetBuff);
//   checkCUDAError ("NeighborList::mallocedDeviceCellList cell list buff");
// }


// void NeighborList::
// clearDeviceCellList () 
// {
//   if ( mallocedDeviceCellList ){
//     cudaFree (dclist.data);
//     cudaFree (dclist.numbers);
//     cudaFree (dclist.numNeighborCell);
//     cudaFree (dclist.neighborCellIndex);
//     cudaFree (dclist.neighborCellShiftNoi);
//     cudaFree (mySendBuff);
//     cudaFree (myTargetBuff);
//     mallocedDeviceCellList = false;
//     checkCUDAError ("NeighborList::clearDeviceCellList");
//   }
// }

// void NeighborList::
// clearDeviceNeighborList()
// {
//   if ( mallocedDeviceNeighborList ){
//     cudaFree (dnlist.data);
//     cudaFree (dnlist.Nneighbor);
//     cudaFree (dnlist.forceIndex);
//     mallocedDeviceNeighborList = false;
//     checkCUDAError ("NeighborList::clearDeviceNeighborList");
//   }
// }

// void NeighborList::
// unbindGlobalTexture ()
// {
//   if ( initedGlobalTexture ){
// #ifndef COORD_IN_ONE_VEC
//     cudaUnbindTexture(global_texRef_neighbor_coordx);
//     cudaUnbindTexture(global_texRef_neighbor_coordy);
//     cudaUnbindTexture(global_texRef_neighbor_coordz);
// #else
//     cudaUnbindTexture(global_texRef_neighbor_coord);
// #endif
//     cudaUnbindTexture(global_texRef_neighbor_type);
//     initedGlobalTexture = false;
//     checkCUDAError ("NeighborList::unbindGlobalTexture");
//   }
// }

// void NeighborList::
// clearJudgeStuff ()
// {
//   if ( mallocedJudgeStuff ){
// #ifndef COORD_IN_ONE_VEC
//     cudaFree (backupCoordx);
//     cudaFree (backupCoordy);
//     cudaFree (backupCoordz);
// #else
//     cudaFree (backupCoord);
// #endif
//     cudaFree (judgeRebuild_tag);
//     mallocedJudgeStuff = false;
//     checkCUDAError ("NeighborList::clearJudgeStuff");
//   }
// }


// void NeighborList::
// clear()
// {
//   clearDeviceCellList();
//   clearDeviceNeighborList();
//   unbindGlobalTexture ();
//   clearJudgeStuff ();
// }  

// NeighborList::~NeighborList()
// {
//   clear();
//   if (mallocedNonBondedForceTable == true){
//     cudaFree (nbForceTable);
//   }
// }


// void NeighborList::
// naivelyBuildDeviceCellList (const MDSystem & sys)
// {  
//   // buildDeviceCellList_clearBuff
//   //     <<<cellGridDim, myBlockDim>>> (
//   // 	  mySendBuff);
//   // err.check ("NeighborList::rebuild, clear buff");
//   prepare_naivelyBuildDeviceCellList
//       <<<cellGridDim, myBlockDim>>> (dclist);
//   naivelyBuildDeviceCellList2        
//       <<<atomGridDim, myBlockDim,
//       buildDeviceCellList_step1_sbuffSize>>> (
// 	  sys.ddata.numAtom,
// #ifndef COORD_IN_ONE_VEC
// 	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz, 
// #else
// 	  sys.ddata.coord,
// #endif
// 	  sys.ddata.coordNoix,
// 	  sys.ddata.coordNoiy,
// 	  sys.ddata.coordNoiz,
// 	  sys.box, dclist,
// 	  err.ptr_de,
// 	  err.ptr_dindex, err.ptr_dscalor);
//   err.check ("NeighborList::naivelyBuildDeviceCellList");
//   checkCUDAError ("NeighborList::naivelyBuildDeviceCellList");

//   buildCellNeighborhood
//       <<<cellGridDim, 1>>> (
// 	  dclist,
// 	  mydivide,
// 	  sys.box.size);
// }

// void NeighborList::
// buildDeviceCellList (const MDSystem & sys)
// {
//   buildDeviceCellList_clearBuff
//       <<<cellGridDim, myBlockDim>>> (
// 	  mySendBuff);
//   err.check ("NeighborList::rebuild, clear buff");
//   buildDeviceCellList_step1 
//       <<<cellGridDim, myBlockDim,
//       buildDeviceCellList_step1_sbuffSize>>> (
// 	  sys.ddata.numAtom, 
// #ifndef COORD_IN_ONE_VEC
// 	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
// #else
// 	  sys.ddata.coord,
// #endif
// 	  sys.ddata.coordNoix, sys.ddata.coordNoiy, sys.ddata.coordNoiz,
// 	  sys.box, dclist, mySendBuff, myTargetBuff,
// 	  err.ptr_de, err.ptr_dindex, err.ptr_dscalor);
//   err.check ("NeighborList::buildDeviceCellList step 1");
//   checkCUDAError ("NeighborList::buildDeviceCellList step 1");
//   buildDeviceCellList_step2
//       <<<cellGridDim, myBlockDim>>> (
// 	  sys.box, dclist, mySendBuff, myTargetBuff, bitDeepth,
// 	  err.ptr_de);
//   err.check ("NeighborList::buildDeviceCellList step2");
//   checkCUDAError ("NeighborList::buildDeviceCellList step2");
// }

// void NeighborList::
// buildDeviceNeighborListCellList (const MDSystem & sys)
// {
//   buildDeviceNeighborList_DeviceCellList 
//       <<<cellGridDim, myBlockDim,
//       buildDeviceNeighborList_DeviceCellList_sbuffSize>>> (
// 	  // <<<cellGridDim, myBlockDim>>> (
// 	  sys.ddata.numAtom, 
// #ifndef COORD_IN_ONE_VEC
// 	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
// #else
// 	  sys.ddata.coord,
// #endif
// 	  sys.ddata.type,
// 	  sys.box, dclist, dnlist,
// 	  nbForceTable, NatomType,
// 	  sharednbForceTable,
// 	  err.ptr_de);
//   err.check("NeighborList::buildDeviceNeighborListCellList");
//   checkCUDAError ("NeighborList::buildDeviceNeighborListCellList");
// }

// void NeighborList::
// buildDeviceNeighborListAllPair (const MDSystem & sys)
// {
//   buildDeviceNeighborList_AllPair 
//       <<<atomGridDim, myBlockDim,
//       buildDeviceNeighborList_AllPair_sbuffSize>>>(
// 	  sys.ddata.numAtom,
// #ifndef COORD_IN_ONE_VEC
// 	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
// #else
// 	  sys.ddata.coord,
// #endif
// 	  sys.ddata.type,
// 	  sys.box, dnlist,
// 	  nbForceTable, NatomType,
// 	  sharednbForceTable,
// 	  err.ptr_de);
//   err.check("NeighborList::build, build neighbor list all pair");
//   checkCUDAError ("NeighborList::build, build neighbor list all pair");
// }

// static IndexType calDeepth (IndexType N)
// {
//   IndexType deepth = 0;
//   while (N != 0){
//     N >>= 1;
//     deepth ++;
//   }
//   return deepth;
// }

// static IndexType hroundUp4 (IndexType x)
// {
//   if (x & 3 == 0){
//     return x;
//   }
//   else {
//     return ((x >> 2) + 1) << 2;
//   }
// }

// void NeighborList::
// initNonBondedInteraction (const SystemNonBondedInteraction & sysNbInter)
// {
//   if (! sysNbInter.beBuilt()) {
//     throw MDExcptUnbuiltNonBondedInteraction ("NeighborList");
//   }
//   NatomType = sysNbInter.numberOfAtomTypes();
//   nbForceTableLength = sysNbInter.interactionTableSize();
//   cudaMalloc ((void**)&nbForceTable,
// 	      nbForceTableLength * sizeof(IndexType));
//   cudaMemcpy (nbForceTable,
// 	      sysNbInter.interactionTable(),
// 	      nbForceTableLength * sizeof(IndexType),
// 	      cudaMemcpyHostToDevice);
//   checkCUDAError ("AtomNBForceTable::deviceInitTable");
//   mallocedNonBondedForceTable = true;
// }

// void NeighborList::
// mallocDeviceNeighborList (const MDSystem & sys,
// 			  const ScalorType & DeviceNeighborListExpansion)
// {
//   dnlist.rlist = myrlist;
//   dnlist.stride = sys.ddata.numAtom;
//   ScalorType density = sys.ddata.numAtom / (sys.box.size.x * sys.box.size.y * sys.box.size.z);
//   IndexType expectedNumberInList 
//       = 4./3. * M_PI * myrlist * myrlist * myrlist * density;
//   dnlist.listLength = IndexType(expectedNumberInList * DeviceNeighborListExpansion);
//   if (dnlist.listLength < 20){
//     dnlist.listLength = 20;
//   }
//   cudaMalloc ((void**)&(dnlist.data), sizeof(IndexType) * dnlist.stride * dnlist.listLength);
//   cudaMalloc ((void**)&(dnlist.Nneighbor), sizeof(IndexType) * sys.ddata.numAtom);
//   cudaMalloc ((void**)&(dnlist.forceIndex), sizeof(IndexType) *  dnlist.stride * dnlist.listLength);
//   // reshuffle backup things
//   cudaMalloc ((void**)&(bkdnlistData), sizeof(IndexType) * dnlist.stride * dnlist.listLength);
//   cudaMalloc ((void**)&(bkdnlistNneighbor), sizeof(IndexType) * sys.ddata.numAtom);
//   cudaMalloc ((void**)&(bkdnlistForceIndex), sizeof(IndexType) *  dnlist.stride * dnlist.listLength);
//   checkCUDAError ("NeighborList::mallocDeviceNeighborList");
//   mallocedDeviceNeighborList = true;
// }

// void NeighborList::
// bindGlobalTexture (const MDSystem & sys)
// {
//   size_t sizetype   = sizeof(TypeType)  *sys.ddata.numMem;
// #ifndef COORD_IN_ONE_VEC
//   size_t sizescalor = sizeof(ScalorType)*sys.ddata.numMem;
//   cudaBindTexture(0, global_texRef_neighbor_coordx, sys.ddata.coordx, sizescalor);
//   cudaBindTexture(0, global_texRef_neighbor_coordy, sys.ddata.coordy, sizescalor);
//   cudaBindTexture(0, global_texRef_neighbor_coordz, sys.ddata.coordz, sizescalor);
// #else
//   size_t sizecoord  = sizeof(CoordType) *sys.ddata.numMem;
//   cudaBindTexture(0, global_texRef_neighbor_coord, sys.ddata.coord, sizecoord);
// #endif
//   cudaBindTexture(0, global_texRef_neighbor_type, sys.ddata.type, sizetype);
//   checkCUDAError ("NeighborList::init texture");
//   initedGlobalTexture = true;
// }

// void NeighborList::
// mallocJudgeStuff(const MDSystem & sys)
// {
//   IndexType nob;
//   if (sys.ddata.numAtom % myBlockDim.x == 0){
//     nob = sys.ddata.numAtom / myBlockDim.x;
//   } else {
//     nob = sys.ddata.numAtom / myBlockDim.x + 1;
//   }
//   cudaMalloc ((void **)& backupCoord,  sizeof(CoordType) *sys.ddata.numAtom);
//   // reshuffle backup
//   cudaMalloc ((void **)& bkbackupCoord,sizeof(CoordType) *sys.ddata.numAtom);
//   cudaMalloc ((void **)& judgeRebuild_tag,  sizeof(IndexType));
//   sum_judge.reinit (nob, NThreadForSum);
//   checkCUDAError ("NeighborList::init judge build allocations");
//   mallocedJudgeStuff = true;
// }


// void NeighborList::
// init (const SystemNonBondedInteraction & sysNbInter,
//       const MDSystem & sys,
//       const ScalorType & rlist,
//       const IndexType & NTread,
//       const ScalorType & DeviceNeighborListExpansion,
//       const BoxDirection_t & bdir,
//       const IndexType & divide)
// {
//   myBlockDim.y = 1;
//   myBlockDim.z = 1;
//   myBlockDim.x = NTread;
//   bitDeepth = calDeepth(sys.ddata.numAtom);
//   printf ("# the bit deepth is %d\n", bitDeepth);
  
//   IndexType nob;
//   if (sys.ddata.numAtom % myBlockDim.x == 0){
//     nob = sys.ddata.numAtom / myBlockDim.x;
//   } else {
//     nob = sys.ddata.numAtom / myBlockDim.x + 1;
//   }
//   atomGridDim = toGridDim (nob);
  
//   mybdir = bdir;
//   myrlist = rlist;
//   mydivide = divide;
//   DecideNeighboringMethod (sys, myrlist, mybdir, mydivide, mode, dclist.NCell);
//   if (mode == CellListBuilt){
//     mallocDeviceCellList (dclist.NCell, sys.box.size);
//   }
  
//   // init neighbor list
//   mallocDeviceNeighborList (sys, DeviceNeighborListExpansion);
  
//   // // nonbonded interaction
//   initNonBondedInteraction (sysNbInter);
  
//   // bind texture
//   bindGlobalTexture (sys);
  
//   // init backup Coords
//   mallocJudgeStuff (sys);

// // #ifndef COORD_IN_ONE_VEC
// //   judgeRebuild_backupCoord
// //       <<<atomGridDim, myBlockDim>>> (sys.ddata.numAtom, sys.ddata.coordx, backupCoordx);
// //   judgeRebuild_backupCoord
// //       <<<atomGridDim, myBlockDim>>> (sys.ddata.numAtom, sys.ddata.coordy, backupCoordy);
// //   judgeRebuild_backupCoord
// //       <<<atomGridDim, myBlockDim>>> (sys.ddata.numAtom, sys.ddata.coordz, backupCoordz);
// // #else
// //   cpyProperty <<<atomGridDim, myBlockDim>>> (
// //       backupCoord, sys.ddata.coord, sys.ddata.numAtom);
// // #endif
// //   checkCUDAError ("NeighborList::init backup coords");

//   //init shared memory size
//   sharednbForceTable = true;
//   buildDeviceNeighborList_DeviceCellList_sbuffSize =
//       sizeof(IndexType) *      hroundUp4(myBlockDim.x) +
// #ifndef COORD_IN_ONE_VEC
//       sizeof(ScalorType) * 3 * hroundUp4(myBlockDim.x) +
// #else
//       sizeof(CoordType) * hroundUp4(myBlockDim.x) +
// #endif
//       sizeof(TypeType)  *      hroundUp4(myBlockDim.x) +
//       sizeof(IndexType) * hroundUp4(nbForceTableLength);
//   if (buildDeviceNeighborList_DeviceCellList_sbuffSize >=
//       SystemSharedBuffSize - GlobalFunctionParamSizeLimit){
//     sharednbForceTable = false;
//     buildDeviceNeighborList_DeviceCellList_sbuffSize =
// 	sizeof(IndexType) *      hroundUp4(myBlockDim.x) +
// #ifndef COORD_IN_ONE_VEC
// 	sizeof(ScalorType) * 3 * hroundUp4(myBlockDim.x) +
// #else
// 	sizeof(CoordType) * hroundUp4(myBlockDim.x) +
// #endif
// 	sizeof(TypeType)  *      hroundUp4(myBlockDim.x);
//   }  
//   buildDeviceNeighborList_AllPair_sbuffSize =
//       buildDeviceNeighborList_DeviceCellList_sbuffSize;
//   buildDeviceCellList_step1_sbuffSize = 2 * myBlockDim.x * sizeof(IndexType);
//   judgeRebuild_judgeCoord_block_sbuffSize = myBlockDim.x * sizeof(IndexType);
// }




// void NeighborList::
// reinit (const MDSystem & sys,
// 	const ScalorType & rlist,
// 	const IndexType & NTread,
// 	const ScalorType & DeviceNeighborListExpansion,
// 	const BoxDirection_t & bdir,
// 	const IndexType & divide)
// {
//   myBlockDim.y = 1;
//   myBlockDim.z = 1;
//   myBlockDim.x = NTread;
//   bitDeepth = calDeepth(sys.ddata.numAtom);
//   printf ("# the bit deepth is %d\n", bitDeepth);

//   IndexType nob;
//   if (sys.ddata.numAtom % myBlockDim.x == 0){
//     nob = sys.ddata.numAtom / myBlockDim.x;
//   } else {
//     nob = sys.ddata.numAtom / myBlockDim.x + 1;
//   }
//   atomGridDim = toGridDim (nob);
  
//   // init cell list
//   mybdir = bdir;
//   myrlist = rlist;
//   mydivide = divide;
//   clearDeviceCellList ();
//   DecideNeighboringMethod (sys, myrlist, mybdir, mydivide, mode, dclist.NCell);
//   if (mode == CellListBuilt){
//     mallocDeviceCellList (dclist.NCell, sys.box.size);
//   }
  
//   // init neighbor list
//   clearDeviceNeighborList ();
//   mallocDeviceNeighborList (sys, DeviceNeighborListExpansion);  
  
//   // bind texture
//   unbindGlobalTexture ();
//   bindGlobalTexture (sys);  
  
//   // init backup Coords
//   clearJudgeStuff ();
//   mallocJudgeStuff (sys);
  
//   //init shared memory size
//   sharednbForceTable = true;
//   buildDeviceNeighborList_DeviceCellList_sbuffSize =
//       sizeof(IndexType) *      hroundUp4(myBlockDim.x) +
// #ifndef COORD_IN_ONE_VEC
//       sizeof(ScalorType) * 3 * hroundUp4(myBlockDim.x) +
// #else
//       sizeof(CoordType) * hroundUp4(myBlockDim.x) +
// #endif
//       sizeof(TypeType)  *      hroundUp4(myBlockDim.x) +
//       sizeof(IndexType) * hroundUp4(nbForceTableLength);
//   if (buildDeviceNeighborList_DeviceCellList_sbuffSize >=
//       SystemSharedBuffSize - GlobalFunctionParamSizeLimit){
//     sharednbForceTable = false;
//     buildDeviceNeighborList_DeviceCellList_sbuffSize =
// 	sizeof(IndexType) *      hroundUp4(myBlockDim.x) +
// #ifndef COORD_IN_ONE_VEC
// 	sizeof(ScalorType) * 3 * hroundUp4(myBlockDim.x) +
// #else
// 	sizeof(CoordType) * hroundUp4(myBlockDim.x) +
// #endif
// 	sizeof(TypeType)  *      hroundUp4(myBlockDim.x);
//   }  
//   buildDeviceNeighborList_AllPair_sbuffSize =
//       buildDeviceNeighborList_DeviceCellList_sbuffSize;
//   buildDeviceCellList_step1_sbuffSize = 2 * myBlockDim.x * sizeof(IndexType);
//   judgeRebuild_judgeCoord_block_sbuffSize = myBlockDim.x * sizeof(IndexType);
// }


// void NeighborList::
// normalizeMDSystem (const MDSystem & sys)
// {
//   normalizeSystem
//       <<<atomGridDim, myBlockDim>>> (
// 	  sys.box, sys.ddata.numAtom,
// #ifndef COORD_IN_ONE_VEC
// 	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
// #else
// 	  sys.ddata.coord,
// #endif
// 	  sys.ddata.coordNoix, sys.ddata.coordNoiy, sys.ddata.coordNoiz);
//   checkCUDAError ("NeighborList::rebuild, normalize System");
// }

// void NeighborList::
// backupJudgeCoord (const MDSystem & sys)
// {
// #ifndef COORD_IN_ONE_VEC
//   judgeRebuild_backupCoord
//       <<<atomGridDim, myBlockDim>>> (sys.ddata.numAtom, sys.ddata.coordx, backupCoordx);
//   judgeRebuild_backupCoord
//       <<<atomGridDim, myBlockDim>>> (sys.ddata.numAtom, sys.ddata.coordy, backupCoordy);
//   judgeRebuild_backupCoord
//       <<<atomGridDim, myBlockDim>>> (sys.ddata.numAtom, sys.ddata.coordz, backupCoordz);
// #else
//   cpyProperty <<<atomGridDim, myBlockDim>>> (
//       backupCoord, sys.ddata.coord, sys.ddata.numAtom);
// #endif
//   checkCUDAError ("NeighborList::init backup coords");
// }

// void NeighborList::build (const MDSystem & sys,
// 			  MDTimer * timer)
// {
//   backupJudgeCoord (sys);
//   if (timer != NULL) timer->tic(mdTimeNormalizeSys);
//   normalizeMDSystem (sys);
//   if (timer != NULL) timer->toc(mdTimeNormalizeSys);
//   if (mode == AllPairBuilt){
//     if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
//     buildDeviceNeighborListAllPair (sys);
//     if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
//   }
//   else if (mode == CellListBuilt){
//     if (timer != NULL) timer->tic(mdTimeBuildCellList);
//     naivelyBuildDeviceCellList (sys);
//     if (timer != NULL) timer->toc(mdTimeBuildCellList);
//     if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
//     buildDeviceNeighborListCellList (sys);
//     if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
//   }
//   else {
//     throw WrongBuildMethod();
//   }
// }



// void NeighborList::reBuild (const MDSystem & sys,
// 			    MDTimer * timer)
// {
//   if (timer != NULL) timer->tic(mdTimeNormalizeSys);
//   normalizeMDSystem (sys);
//   if (timer != NULL) timer->toc(mdTimeNormalizeSys);
//   NeighborListBuiltMode tmpMode;
//   IntVectorType tmpNCell;
//   if (timer != NULL) timer->tic(mdTimeBuildCellList);
//   DecideNeighboringMethod (sys, myrlist, mybdir, mydivide, tmpMode, tmpNCell);

//   // printf("# rebuild %d %d %d\n", tmpNCell.x, tmpNCell.y, tmpNCell.z);
  
//   if (tmpMode == mode) {
//     if (mode == AllPairBuilt){
//       if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
//       buildDeviceNeighborListAllPair (sys);
//       if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
//     }
//     else if (mode == CellListBuilt){
//       if (timer != NULL) timer->tic(mdTimeBuildCellList);
//       if (tmpNCell.x != dclist.NCell.x ||
// 	  tmpNCell.y != dclist.NCell.y ||
// 	  tmpNCell.z != dclist.NCell.z ){
// 	printf ("# box size change too much, rebuild cell list\n");
// 	clearDeviceCellList ();
// 	mallocDeviceCellList (tmpNCell, sys.box.size);
// 	naivelyBuildDeviceCellList (sys);
//       }
//       else {
// 	buildDeviceCellList (sys);
//       }
//       if (timer != NULL) timer->toc(mdTimeBuildCellList);
//       if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
//       buildDeviceNeighborListCellList (sys);
//       if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
//     }
//     else {
//       throw WrongBuildMethod();
//     }    
//   }
//   else if (mode == AllPairBuilt && tmpMode == CellListBuilt){
//     mode = tmpMode;
//     if (timer != NULL) timer->tic(mdTimeBuildCellList);
//     mallocDeviceCellList (tmpNCell, sys.box.size);
//     naivelyBuildDeviceCellList (sys);
//     if (timer != NULL) timer->toc(mdTimeBuildCellList);
//     if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
//     buildDeviceNeighborListCellList (sys);
//     if (timer != NULL) timer->toc(mdTimeBuildNeighborList);    
//   }
//   else if (mode == CellListBuilt && tmpMode == AllPairBuilt){
//     mode = tmpMode;
//     clearDeviceCellList();
//     if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
//     buildDeviceNeighborListAllPair (sys);
//     if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
//   }
//   else {
//     throw WrongBuildMethod();
//   }
// }


// void NeighborList::buildCellList (const MDSystem & sys,
// 				  MDTimer * timer)
// {
//   backupJudgeCoord (sys);
//   if (timer != NULL) timer->tic(mdTimeNormalizeSys);
//   normalizeMDSystem (sys);
//   if (timer != NULL) timer->toc(mdTimeNormalizeSys);
//   if (mode == AllPairBuilt){
//     if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
//     buildDeviceNeighborListAllPair (sys);
//     if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
//   }
//   else if (mode == CellListBuilt){
//     if (timer != NULL) timer->tic(mdTimeBuildCellList);
//     naivelyBuildDeviceCellList (sys);
//     if (timer != NULL) timer->toc(mdTimeBuildCellList);
//   }
//   else {
//     throw WrongBuildMethod();
//   }
// }



// void NeighborList::reBuildCellList (const MDSystem & sys,
// 				    MDTimer * timer)
// {
//   if (timer != NULL) timer->tic(mdTimeNormalizeSys);
//   normalizeMDSystem (sys);
//   if (timer != NULL) timer->toc(mdTimeNormalizeSys);
//   NeighborListBuiltMode tmpMode;
//   IntVectorType tmpNCell;
//   if (timer != NULL) timer->tic(mdTimeBuildCellList);
//   DecideNeighboringMethod (sys, myrlist, mybdir, mydivide, tmpMode, tmpNCell);

//   // printf("# rebuild %d %d %d\n", tmpNCell.x, tmpNCell.y, tmpNCell.z);
  
//   if (tmpMode == mode) {
//     if (mode == AllPairBuilt){
//       if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
//       buildDeviceNeighborListAllPair (sys);
//       if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
//     }
//     else if (mode == CellListBuilt){
//       if (timer != NULL) timer->tic(mdTimeBuildCellList);
//       if (tmpNCell.x != dclist.NCell.x ||
// 	  tmpNCell.y != dclist.NCell.y ||
// 	  tmpNCell.z != dclist.NCell.z ){
// 	printf ("# box size change too much, rebuild cell list\n");
// 	clearDeviceCellList ();
// 	mallocDeviceCellList (tmpNCell, sys.box.size);
// 	naivelyBuildDeviceCellList (sys);
//       }
//       else {
// 	buildDeviceCellList (sys);
//       }
//       if (timer != NULL) timer->toc(mdTimeBuildCellList);
//     }
//     else {
//       throw WrongBuildMethod();
//     }    
//   }
//   else if (mode == AllPairBuilt && tmpMode == CellListBuilt){
//     mode = tmpMode;
//     if (timer != NULL) timer->tic(mdTimeBuildCellList);
//     mallocDeviceCellList (tmpNCell, sys.box.size);
//     naivelyBuildDeviceCellList (sys);
//     if (timer != NULL) timer->toc(mdTimeBuildCellList);
//   }
//   else if (mode == CellListBuilt && tmpMode == AllPairBuilt){
//     mode = tmpMode;
//     clearDeviceCellList();
//     if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
//     buildDeviceNeighborListAllPair (sys);
//     if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
//   }
//   else {
//     throw WrongBuildMethod();
//   }
// }



// bool NeighborList::judgeRebuild (const MDSystem & sys,
// 				 const ScalorType &diffTol,
// 				 MDTimer *timer)
// {
//   if (timer != NULL) timer->tic(mdTimeJudgeRebuild);
//   IndexType mytag;
//   ScalorType diffTol2 = diffTol * diffTol;
  
//   // judgeRebuild_judgeCoord
//   //     <<<atomGridDim, myBlockDim>>> (
//   // 	  sys.ddata.numAtom, sys.box,
//   // 	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
//   // 	  backupCoordx, backupCoordy, backupCoordz,
//   // 	  diffTol2,
//   // 	  judgeRebuild_buff, judgeRebuild_tag);
//   judgeRebuild_judgeCoord_block 
//       <<<atomGridDim, myBlockDim,
//       judgeRebuild_judgeCoord_block_sbuffSize>>> (
//   	  sys.ddata.numAtom, sys.box,
// #ifndef COORD_IN_ONE_VEC
//   	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
//   	  backupCoordx, backupCoordy, backupCoordz,
// #else
// 	  sys.ddata.coord,
// 	  backupCoord,
// #endif
//   	  diffTol2,
//   	  sum_judge.getBuff());
//   sum_judge.sumBuff (judgeRebuild_tag, 0);	  
//   cudaMemcpy (&mytag, judgeRebuild_tag, sizeof(IndexType), cudaMemcpyDeviceToHost);
//   if (mytag != 0){
// #ifndef COORD_IN_ONE_VEC
//     judgeRebuild_backupCoord
// 	<<<atomGridDim, myBlockDim>>> (
// 	    sys.ddata.numAtom, sys.ddata.coordx, backupCoordx);
//     judgeRebuild_backupCoord
// 	<<<atomGridDim, myBlockDim>>> (
// 	    sys.ddata.numAtom, sys.ddata.coordy, backupCoordy);
//     judgeRebuild_backupCoord
// 	<<<atomGridDim, myBlockDim>>> (
// 	    sys.ddata.numAtom, sys.ddata.coordz, backupCoordz);
// #else
//     cpyProperty <<<atomGridDim, myBlockDim>>> (
// 	backupCoord, sys.ddata.coord, sys.ddata.numAtom);
// #endif    
//     if (timer != NULL) timer->toc(mdTimeJudgeRebuild);
//     // fprintf (stderr, "rebuild!\n");
//     return true;
//   }
//   else{
//     if (timer != NULL) timer->toc(mdTimeJudgeRebuild);
//     return false;
//   }
// }


// void NeighborList::
// reshuffleCell (const IndexType * indexTable,
// 	       const IndexType & numAtom,
// 	       MDTimer *timer)
// {
//   if (timer != NULL) timer->tic(mdTimeReshuffleSystem);

//   Reshuffle_reshuffleDeviceCellList
//       <<<cellGridDim, myBlockDim>>> (
// 	  dclist.data, indexTable);
//   if (timer != NULL) timer->toc(mdTimeReshuffleSystem);
// }


// void NeighborList::
// reshuffle (const IndexType * indexTable,
// 	   const IndexType & numAtom,
// 	   MDTimer *timer)
// {
//   if (timer != NULL) timer->tic(mdTimeReshuffleSystem);

//   Reshuffle_reshuffleDeviceCellList
//       <<<cellGridDim, myBlockDim>>> (
// 	  dclist.data, indexTable);
//   cudaMemcpy (bkbackupCoord, backupCoord,
// 	      sizeof (CoordType) * numAtom,
// 	      cudaMemcpyDeviceToDevice);
//   Reshuffle_reshuffleArray
//       <<<atomGridDim, myBlockDim>>> 
//       (bkbackupCoord, numAtom, indexTable, backupCoord);
//   Reshuffle_backupDeviceNeighborList
//       <<<atomGridDim, myBlockDim,
//       2 * myBlockDim.x * sizeof(IndexType)>>> (
// 	  numAtom,
// 	  dnlist.data,
// 	  dnlist.forceIndex,
// 	  dnlist.stride,
// 	  dnlist.Nneighbor,
// 	  bkdnlistData,
// 	  bkdnlistForceIndex,
// 	  bkdnlistNneighbor);
//   Reshuffle_reshuffleDeviceNeighborList
//       <<<atomGridDim, myBlockDim,
//       2 * myBlockDim.x * sizeof(IndexType)>>> (
// 	  numAtom,
// 	  bkdnlistData,
// 	  bkdnlistForceIndex,
// 	  dnlist.stride,
// 	  bkdnlistNneighbor,
// 	  indexTable,
// 	  dnlist.data,
// 	  dnlist.forceIndex,
// 	  dnlist.Nneighbor);
//   checkCUDAError ("NeighborList::reshuffle");
//   if (timer != NULL) timer->toc(mdTimeReshuffleSystem);
// }



// ////////////////////////////////////////////////////////////
// // for the reason of using texture, we place this function here. it
// // should be placed in NeighborList.cu
// ////////////////////////////////////////////////////////////


// __device__ IndexType shiftedD3toD1 (
//     DeviceCellList clist,
//     RectangularBoxGeometry::RectangularBox box,
//     int ix, int iy, int iz,
//     ScalorType * shiftx , ScalorType * shifty, ScalorType * shiftz)
// {
//   int tmp;
//   ix += (tmp = -int(floorf(ix * clist.NCelli.x))) * clist.NCell.x;
//   *shiftx = tmp * box.size.x;
//   iy += (tmp = -int(floorf(iy * clist.NCelli.y))) * clist.NCell.y;
//   *shifty = tmp * box.size.y;
//   iz += (tmp = -int(floorf(iz * clist.NCelli.z))) * clist.NCell.z;
//   *shiftz = tmp * box.size.z;
//   return D3toD1 (clist.NCell, ix, iy, iz);
// }

// __global__ void buildDeviceNeighborList_DeviceCellList (
//     IndexType numAtom,
//     CoordType * coord,
//     TypeType * type,
//     RectangularBox box,
//     DeviceCellList clist,
//     DeviceNeighborList nlist,
//     const IndexType *nbForceTable,
//     IndexType NatomType,
//     bool sharednbForceTable,
//     mdError_t * ptr_de)
// {
//   // RectangularBoxGeometry::normalizeSystem (box, &ddata);
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   IndexType bidx, bidy, bidz;
//   D1toD3 (clist.NCell, bid, bidx, bidy, bidz);
  
//   // set number of neighbor to 0
//   IndexType Nneighbor = 0;
//   // load index
//   IndexType ii = getDeviceCellListData (clist, bid, tid);
//   // load iith coordinate // use texturefetch instead
//   CoordType ref;
//   TypeType reftype;
//   if (ii != MaxIndexValue){
// #ifdef COMPILE_NO_TEX
//     ref = coord[ii];
//     reftype = type[ii];
// #else
//     ref = tex1Dfetch (global_texRef_neighbor_coord, ii);
//     reftype = tex1Dfetch(global_texRef_neighbor_type, ii);
// #endif
//   }
//   ScalorType rlist = clist.rlist;

//   // the target index and coordinates are shared

//   extern __shared__ volatile char pub_sbuff[];
  
//   volatile IndexType * targetIndexes =
//       (volatile IndexType *) pub_sbuff;
//   CoordType * target =
//       (CoordType *) &targetIndexes[roundUp4(blockDim.x)];
//   volatile TypeType * targettype =
//       (volatile TypeType *) &target[roundUp4(blockDim.x)];
//   IndexType * nbForceTableBuff = NULL;

//   IndexType nbForceTableLength = AtomNBForceTable::dCalDataLength(NatomType);
//   if (sharednbForceTable){
//     nbForceTableBuff = (IndexType *) &targettype[roundUp4(blockDim.x)];
//     cpyGlobalDataToSharedBuff (nbForceTable, nbForceTableBuff, nbForceTableLength);
//   }
//   __syncthreads();

//   // __shared__ volatile  IndexType  targetIndexes [MaxThreadsPerBlock];
//   // __shared__ volatile CoordType target [MaxThreadsPerBlock];
//   // __shared__ volatile  TypeType   targettype    [MaxThreadsPerBlock];
//   // __shared__ volatile  IndexType nbForceTableBuff [MaxNBForceTableBuffSize];

//   // IndexType nbForceTableLength = AtomNBForceTable::dCalDataLength(NatomType);
//   // if (sharednbForceTable){
//   //   cpyGlobalDataToSharedBuff (nbForceTable, nbForceTableBuff, nbForceTableLength);
//   // }
//   // __syncthreads();

//   bool oneCellX(false), oneCellY(false), oneCellZ(false);
//   if (clist.NCell.x == 1) oneCellX = true;
//   if (clist.NCell.y == 1) oneCellY = true;
//   if (clist.NCell.z == 1) oneCellZ = true;
//   ScalorType rlist2 = rlist * rlist;
  
//   // loop over 27 neighbor cells
//   for (IndexType i = 0; i < clist.numNeighborCell[bid]; ++i){
//     __syncthreads();
//     IndexType targetCellIdx = getNeighborCellIndex    (clist, bid, i);
//     CoordNoiType shiftNoi   = getNeighborCellShiftNoi (clist, bid, i);
//     CoordType shift;
//     shift.x = shiftNoi.x * box.size.x;
//     shift.y = shiftNoi.y * box.size.y;
//     shift.z = shiftNoi.z * box.size.z;
//     targetIndexes[tid] = getDeviceCellListData(clist, targetCellIdx, tid);  
//     if (targetIndexes[tid] != MaxIndexValue){
//       target[tid] = tex1Dfetch(global_texRef_neighbor_coord, targetIndexes[tid]);
//       targettype[tid] = tex1Dfetch(global_texRef_neighbor_type, targetIndexes[tid]);
//     }
//     __syncthreads();
	
//     // find neighbor
//     if (ii != MaxIndexValue){
//       for (IndexType jj = 0; jj < clist.numbers[targetCellIdx]; ++jj){
// 	ScalorType diffx = target[jj].x - shift.x - ref.x;
// 	ScalorType diffy = target[jj].y - shift.y - ref.y;
// 	ScalorType diffz = target[jj].z - shift.z - ref.z;
// 	if (oneCellX) shortestImage (box.size.x, box.sizei.x, &diffx);
// 	if (oneCellY) shortestImage (box.size.y, box.sizei.y, &diffy);
// 	if (oneCellZ) shortestImage (box.size.z, box.sizei.z, &diffz);
// 	if ((diffx*diffx+diffy*diffy+diffz*diffz) < rlist2 &&
// 	    targetIndexes[jj] != ii){
// 	  IndexType fidx;
// 	  if (sharednbForceTable){
// 	    fidx = AtomNBForceTable::calForceIndex (
// 		nbForceTableBuff, NatomType, reftype, targettype[jj]);
// 	  }
// 	  else {
// 	    fidx = AtomNBForceTable::calForceIndex (
// 		nbForceTable, NatomType, reftype, targettype[jj]);
// 	  }
// 	  // if (fidx != mdForceNULL) {
// 	  IndexType listIdx = Nneighbor * nlist.stride + ii;
// 	  nlist.data[listIdx] = targetIndexes[jj];
// 	  nlist.forceIndex[listIdx] = fidx;
// 	  Nneighbor ++;
// 	  // }
// 	}
//       }
//     }
//   }

//   if (ii != MaxIndexValue) {
//     if (Nneighbor > nlist.listLength && ptr_de != NULL){
//       *ptr_de = mdErrorShortNeighborList;
//       return;
//     }
//     nlist.Nneighbor[ii] = Nneighbor;
//   }
// }



// __global__ void buildDeviceCellList_step1 (IndexType numAtom,
// 					   CoordType * coord,
// 					   IntScalorType * coordNoix,
// 					   IntScalorType * coordNoiy,
// 					   IntScalorType * coordNoiz,
// 					   RectangularBox box,
// 					   DeviceCellList clist,
// 					   IndexType * sendBuff,
// 					   IndexType * targetBuff,
// 					   mdError_t * ptr_de,
// 					   IndexType * erridx,
// 					   ScalorType * errsrc)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;

//   extern __shared__ volatile IndexType sbuff[];
//   volatile IndexType * originalData = (volatile IndexType *) sbuff;
//   volatile IndexType * targetCellid = (volatile IndexType *) &originalData[blockDim.x];
  
//   // __shared__ volatile  IndexType originalData[MaxThreadsPerBlock];
//   // __shared__ volatile  IndexType targetCellid[MaxThreadsPerBlock];
  
//   // copy data from cell list
//   originalData[tid] = clist.data[bid*clist.stride + tid];
//   IndexType originalNumber = clist.numbers[bid];

//   // calculate the target cell
//   if (originalData[tid] != MaxIndexValue){
//     IndexType targetCelli, targetCellj, targetCellk;
//     IndexType thisid = originalData[tid];
// #ifdef COMPILE_NO_TEX
//     ref = coord[thisid];
// #else
//     CoordType ref (tex1Dfetch(global_texRef_neighbor_coord, thisid));
// #endif
//     targetCelli = IndexType(ref.x * box.sizei.x * ScalorType (clist.NCell.x));
//     targetCellj = IndexType(ref.y * box.sizei.y * ScalorType (clist.NCell.y));
//     targetCellk = IndexType(ref.z * box.sizei.z * ScalorType (clist.NCell.z));
//     if (targetCelli == clist.NCell.x){
//       targetCelli -= clist.NCell.x;
//       coord[thisid].x -= box.size.x;
//       coordNoix[thisid] ++;
//     }
//     if (targetCellj == clist.NCell.y){
//       targetCellj -= clist.NCell.y;
//       coord[thisid].y -= box.size.y;
//       coordNoiy[thisid] ++;
//     }
//     if (targetCellk == clist.NCell.z){
//       targetCellk -= clist.NCell.z;
//       coord[thisid].z -= box.size.z;
//       coordNoiz[thisid] ++;
//     }
//     targetCellid[tid] = D3toD1 (clist.NCell, targetCelli, targetCellj, targetCellk);
//     if (ptr_de != NULL && 
// 	(targetCelli >= clist.NCell.x || 
// 	 targetCellj >= clist.NCell.y || 
// 	 targetCellk >= clist.NCell.z)){
//       *ptr_de = mdErrorOverFlowCellIdx;
//       if (targetCelli >= IndexType(clist.NCell.x)){
// 	*erridx = targetCelli;
// 	*errsrc = ref.x;
// 	return;
//       }
//       if (targetCellj >= IndexType(clist.NCell.y)){
// 	*erridx = targetCellj;
// 	*errsrc = ref.y;
// 	return;
//       }
//       if (targetCellk >= IndexType(clist.NCell.z)){
// 	*erridx = targetCellk;
// 	*errsrc = ref.z;
// 	return;
//       }      
//     }
//   }
//   else {
//     targetCellid[tid] = MaxIndexValue;
//   }

//   // mark particles to be send 
//   IndexType mark = MaxIndexValue - (MaxIndexValue >> 1);
//   if (tid < originalNumber && targetCellid[tid] != bid){
//     originalData[tid] += mark;
//   }
  
//   // head sort
//   IndexType total1 = headSort (originalData, targetCellid);
//   IndexType total0 = blockDim.x - total1;
  
//   // unmark and copy to send buff
//   if (tid < originalNumber && targetCellid[tid] != bid){
//     sendBuff  [bid*clist.stride + tid - total0] = originalData[tid] - mark;
//     targetBuff[bid*clist.stride + tid - total0] = targetCellid[tid];
//     originalData[tid] = MaxIndexValue;
//   }
//   __syncthreads();
  
//   // modify cell list
//   clist.data[bid*clist.stride + tid] = originalData[tid];
//   if (tid == 0) clist.numbers[bid] = total0;
// }



// __global__ void buildDeviceCellList_initBuff (IndexType * sendBuff,
// 					      IndexType * targetBuff)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   IndexType ii = tid + bid * blockDim.x;
//   sendBuff[ii] = MaxIndexValue;
//   targetBuff[ii] = 0;
// }

// __global__ void buildDeviceCellList_clearBuff (IndexType * sendBuff)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType ii = threadIdx.x + bid * blockDim.x;
//   sendBuff[ii] = MaxIndexValue;
// }

  
  
// __global__ void buildDeviceCellList_step2 (RectangularBox box,
// 					   DeviceCellList clist,
// 					   IndexType * sendBuff,
// 					   IndexType * targetBuff,
// 					   IndexType bitDeepth,
// 					   mdError_t * ptr_de)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;

//   IndexType thisid;
//   IndexType ii = 0;
//   IndexType buffPosi;
//   if (tid == 0){
//     while ((thisid = sendBuff[buffPosi = (bid*clist.stride + ii)]) != MaxIndexValue){
//       IndexType cellid = targetBuff[buffPosi];
//       IndexType tailIdx = atomicInc(&clist.numbers[cellid], blockDim.x);
//       if (tailIdx >= blockDim.x&& ptr_de != NULL) {
// 	*ptr_de = mdErrorShortCellList;
// 	return;
//       } 
//       clist.data[cellid * clist.stride + tailIdx] = thisid;
//       sendBuff[buffPosi] = MaxIndexValue;
//       ii ++;
//     }
//   }
// }

// __global__ void judgeRebuild_backupCoord (const IndexType numAtom,
// 					  const ScalorType * from,
// 					  ScalorType * to)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType ii = threadIdx.x + bid * blockDim.x;
//   if (ii < numAtom){
//     to[ii] = from[ii];
//   }
// }


// __device__ IndexType judgeRebuild_counter = 0;
// __global__ void judgeRebuild_judgeCoord (const IndexType numAtom,
// 					 const RectangularBox box,
// 					 const CoordType * coord,
// 					 const CoordType * backupCoord,
// 					 const ScalorType diffTol2,
// 					 ScalorType * judgeRebuild_buff,
// 					 IndexType * tag)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType ii = threadIdx.x + bid * blockDim.x;

//   __shared__ volatile ScalorType mydiff2[MaxThreadsPerBlock * 2];

//   ScalorType dx=0.f, dy=0.f, dz=0.f;
//   if (ii < numAtom){
//     dx = coord[ii].x - backupCoord[ii].x;
//     dy = coord[ii].y - backupCoord[ii].y;
//     dz = coord[ii].z - backupCoord[ii].z;
//     shortestImage (box, &dx, &dy, &dz);
//   }
//   mydiff2[threadIdx.x] = dx*dx + dy*dy + dz*dz;
//   mydiff2[threadIdx.x + blockDim.x] = 0.f;
//   __syncthreads();
//   IndexType num = ((bid+1) * blockDim.x > numAtom) ? 
//       (blockDim.x - ((bid+1) * blockDim.x) + numAtom) : blockDim.x;
//   ScalorType partialMax = maxVectorBlockBuffer (mydiff2, num);
//   ScalorType maxDiff2;
  
//   if (maxPartialMax (partialMax, judgeRebuild_buff, &judgeRebuild_counter, &maxDiff2)){
//     if (maxDiff2 > diffTol2){
//       *tag = 1;
//     }
//     else {
//       *tag = 0;
//     }
//   }
// }


// __global__ void judgeRebuild_judgeCoord_block (const IndexType numAtom,
// 					       const RectangularBox box,
// 					       const CoordType * coord,
// 					       const CoordType * backupCoord,
// 					       const ScalorType diffTol2,
// 					       IndexType * judgeRebuild_buff)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType ii = threadIdx.x + bid * blockDim.x;

//   extern __shared__ volatile IndexType sbuff[];  
//   // __shared__ volatile IndexType sbuff [MaxThreadsPerBlock];

//   ScalorType dx=0.f, dy=0.f, dz=0.f;
//   if (ii < numAtom){
//     dx = coord[ii].x - backupCoord[ii].x;
//     dy = coord[ii].y - backupCoord[ii].y;
//     dz = coord[ii].z - backupCoord[ii].z;
//     shortestImage (box, &dx, &dy, &dz);
//   }

//   if (ii < numAtom){
//     sbuff[threadIdx.x] = ((dx*dx + dy*dy + dz*dz) > diffTol2);
//   }
//   else {
//     sbuff[threadIdx.x] = 0;
//   }
//   __syncthreads();
//   sumVectorBlockBuffer_2 (sbuff);
//   if (threadIdx.x == 0) judgeRebuild_buff[bid] = (sbuff[0] != 0);
// }

