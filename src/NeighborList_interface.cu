/**
 * @file   NeighborList_interface.cu
 * @author Han Wang <han_wang@math.pku.edu.cn>
 * @date   Thu Nov 19 12:53:42 2009
 * 
 * @brief  Implementation of neighbor list
 * 
 * 
 */
#include "NeighborList_interface.h"
#include "Auxiliary.h"
#include "NeighborList.h"
#include <stdio.h>

/** 
 * these are textures for a fast reference of particle position.
 * 
 */
#ifndef COORD_IN_ONE_VEC
texture<ScalorType, 1, cudaReadModeElementType> global_texRef_neighbor_coordx;
texture<ScalorType, 1, cudaReadModeElementType> global_texRef_neighbor_coordy;
texture<ScalorType, 1, cudaReadModeElementType> global_texRef_neighbor_coordz;
#else
texture<CoordType, 1, cudaReadModeElementType> global_texRef_neighbor_coord;
#endif
texture<TypeType, 1, cudaReadModeElementType> global_texRef_neighbor_type;

NeighborList::~NeighborList()
{
  if (malloced){
    if (mode == CellListBuilt){
      cudaFree (dclist.data);
      cudaFree (dclist.numbers);
    }
    cudaFree (dnlist.data);
    cudaFree (dnlist.Nneighbor);
    cudaFree (dnlist.forceIndex);

    cudaFree (nbForceTable);
#ifndef COORD_IN_ONE_VEC
    cudaUnbindTexture(global_texRef_neighbor_coordx);
    cudaUnbindTexture(global_texRef_neighbor_coordy);
    cudaUnbindTexture(global_texRef_neighbor_coordz);
#else
    cudaUnbindTexture(global_texRef_neighbor_coord);
#endif
    cudaUnbindTexture(global_texRef_neighbor_type);
#ifndef COORD_IN_ONE_VEC
    cudaFree (backupCoordx);
    cudaFree (backupCoordy);
    cudaFree (backupCoordz);
#else
    cudaFree (backupCoord);
#endif
    cudaFree (judgeRebuild_tag);
    checkCUDAError ("NeighborList::~NeighborList");
  }
}

static IndexType calDeepth (IndexType N)
{
  IndexType deepth = 0;
  while (N != 0){
    N >>= 1;
    deepth ++;
  }
  return deepth;
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

void NeighborList::initCellList (const MDSystem & sys,
				 const ScalorType & rlist,
				 const BoxDirection_t & bdir)
{
  bool CellOnX, CellOnY, CellOnZ;
  CellOnX = bdir & RectangularBoxGeometry::mdRectBoxDirectionX;
  CellOnY = bdir & RectangularBoxGeometry::mdRectBoxDirectionY;
  CellOnZ = bdir & RectangularBoxGeometry::mdRectBoxDirectionZ;
  IndexType numCell;
  double rlisti = 1./rlist;
  if (CellOnX ){
    dclist.NCelli.x = 1./ (dclist.NCell.x = floorf(sys.box.size.x * rlisti));
  }
  else {
    dclist.NCell.x = 1;
    dclist.NCelli.x = 1./sys.box.size.x;
  }
  if (CellOnY ){
    dclist.NCelli.y = 1./ (dclist.NCell.y = floorf(sys.box.size.y * rlisti));
  }
  else {
    dclist.NCell.y = 1;
    dclist.NCelli.y = 1./sys.box.size.y;
  }
  if (CellOnZ){
    dclist.NCelli.z = 1./ (dclist.NCell.z = floorf(sys.box.size.z * rlisti));
  }
  else{
    dclist.NCell.z = 1;
    dclist.NCelli.z = 1./sys.box.size.z;
  }
  numCell = dclist.NCell.x * dclist.NCell.y * dclist.NCell.z;
  dclist.rlist = rlist;
  // suppose the number of atoms in any cell is smaller or equal
  // to the number of threads in a block
  dclist.stride = myBlockDim.x;

  cellGridDim = toGridDim (numCell);
  
  mode = CellListBuilt;
  if (CellOnX && dclist.NCell.x < 4) mode = AllPairBuilt;
  if (CellOnY && dclist.NCell.y < 4) mode = AllPairBuilt;
  if (CellOnZ && dclist.NCell.z < 4) mode = AllPairBuilt;
  if (mode == CellListBuilt){
    cudaMalloc ((void**)&(dclist.data), 
		sizeof(ScalorType) * numCell * dclist.stride);
    cudaMalloc ((void**)&mySendBuff, 
		sizeof(ScalorType) * numCell * dclist.stride);
    cudaMalloc ((void**)&myTargetBuff, 
		sizeof(ScalorType) * numCell * dclist.stride);
    cudaMalloc ((void**)&(dclist.numbers), sizeof(unsigned) * numCell);
    checkCUDAError ("NeighborList::init cell list");
    buildDeviceCellList_initBuff<<<cellGridDim, myBlockDim>>> 
	(mySendBuff, myTargetBuff);
    checkCUDAError ("NeighborList::init cell list buff");
  }			  
}

void NeighborList::reinitCellList (const MDSystem & sys,
				   const ScalorType & rlist,
				   const BoxDirection_t & bdir)
{
  if (mode == CellListBuilt){
    cudaFree (dclist.data);
    cudaFree (dclist.numbers);
    cudaFree (mySendBuff);
    cudaFree (myTargetBuff);
  }
  initCellList (sys, rlist, bdir);
}


void NeighborList::initNonBondedInteraction (const MDSystem & sys)
{
  if (! sys.nbInter.isBuilt) {
    throw MDExcptUnbuiltNonBondedInteraction ("NeighborList");
  }
  NatomType = sys.nbInter.numAtomTypes;
  nbForceTableLength = sys.nbInter.interactionTableLength;
  cudaMalloc ((void**)&nbForceTable, nbForceTableLength * sizeof(ForceIndexType));
  cudaMemcpy (nbForceTable, sys.nbInter.interactionTable,
	      nbForceTableLength * sizeof(ForceIndexType),
	      cudaMemcpyHostToDevice);
  checkCUDAError ("AtomNBForceTable::deviceInitTable");
}


void NeighborList::init (const MDSystem & sys,
			 const ScalorType & rlist,
			 const IndexType & NTread,
			 const IndexType & DeviceNeighborListExpansion,
			 const BoxDirection_t & bdir)
{
  malloced = false;
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NTread;
  bitDeepth = calDeepth(sys.ddata.numAtom);
  printf ("# the bit deepth is %d\n", bitDeepth);

  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);
  
  // init cell list
  mybdir = bdir;
  initCellList (sys, rlist, mybdir);
  
  // init neighbor list
  dnlist.rlist = rlist;
  dnlist.stride = sys.ddata.numAtom;
  ScalorType density = sys.ddata.numAtom / (sys.box.size.x * sys.box.size.y * sys.box.size.z);
  IndexType expectedNumberInList 
      = 4./3. * M_PI * rlist * rlist * rlist * density;
  dnlist.listLength = expectedNumberInList * DeviceNeighborListExpansion;
  if (dnlist.listLength < 10){
    dnlist.listLength = 10;
  }
  cudaMalloc ((void**)&(dnlist.data), sizeof(IndexType) * dnlist.stride * dnlist.listLength);
  cudaMalloc ((void**)&(dnlist.Nneighbor), sizeof(IndexType) * sys.ddata.numAtom);
  cudaMalloc ((void**)&(dnlist.forceIndex), sizeof(ForceIndexType) *  dnlist.stride * dnlist.listLength);
  checkCUDAError ("NeighborList::init allocate dnlist");

  // // nonbonded interaction
  // NatomType = sys.nbForce.indexTable.NAtomType;
  // // AtomNBForceTable::deviceInitTable (sys.nbForce.indexTable.data,
  // // 				     &nbForceTable, NatomType);
  // nbForceTableLength = AtomNBForceTable::calDataLength(NatomType);
  // cudaMalloc ((void**)&nbForceTable, nbForceTableLength * sizeof(ForceIndexType));
  // cudaMemcpy (nbForceTable, sys.nbForce.indexTable.data,
  // 	      nbForceTableLength * sizeof(ForceIndexType),
  // 	      cudaMemcpyHostToDevice);
  // checkCUDAError ("AtomNBForceTable::deviceInitTable");
  initNonBondedInteraction (sys);
  
  // bind texture
  // printf ("numMem is %d\n", sys.ddata.numMem);
  size_t sizetype   = sizeof(TypeType)  *sys.ddata.numMem;
#ifndef COORD_IN_ONE_VEC
  size_t sizescalor = sizeof(ScalorType)*sys.ddata.numMem;
  cudaBindTexture(0, global_texRef_neighbor_coordx, sys.ddata.coordx, sizescalor);
  cudaBindTexture(0, global_texRef_neighbor_coordy, sys.ddata.coordy, sizescalor);
  cudaBindTexture(0, global_texRef_neighbor_coordz, sys.ddata.coordz, sizescalor);
#else
  size_t sizecoord  = sizeof(CoordType) *sys.ddata.numMem;
  cudaBindTexture(0, global_texRef_neighbor_coord, sys.ddata.coord, sizecoord);
#endif
  cudaBindTexture(0, global_texRef_neighbor_type, sys.ddata.type, sizetype);
  checkCUDAError ("NeighborList::init texture");

  // init backup Coords
#ifndef COORD_IN_ONE_VEC
  cudaMalloc ((void **)& backupCoordx, sizeof(ScalorType)*sys.ddata.numAtom);
  cudaMalloc ((void **)& backupCoordy, sizeof(ScalorType)*sys.ddata.numAtom);
  cudaMalloc ((void **)& backupCoordz, sizeof(ScalorType)*sys.ddata.numAtom);
#else
  cudaMalloc ((void **)& backupCoord,  sizeof(CoordType) *sys.ddata.numAtom);
#endif
  cudaMalloc ((void **)& judgeRebuild_tag,  sizeof(IndexType));
  sum_judge.init (nob, NThreadForSum);
  checkCUDAError ("NeighborList::init judge build allocations");
#ifndef COORD_IN_ONE_VEC
  judgeRebuild_backupCoord
      <<<atomGridDim, myBlockDim>>> (sys.ddata.numAtom, sys.ddata.coordx, backupCoordx);
  judgeRebuild_backupCoord
      <<<atomGridDim, myBlockDim>>> (sys.ddata.numAtom, sys.ddata.coordy, backupCoordy);
  judgeRebuild_backupCoord
      <<<atomGridDim, myBlockDim>>> (sys.ddata.numAtom, sys.ddata.coordz, backupCoordz);
#else
  cpyProperty <<<atomGridDim, myBlockDim>>> (
      backupCoord, sys.ddata.coord, sys.ddata.numAtom);
#endif
  checkCUDAError ("NeighborList::init backup coords");
  malloced = true;

  //init shared memory size
  sharednbForceTable = true;
  buildDeviceNeighborList_DeviceCellList_sbuffSize =
      sizeof(IndexType) *      hroundUp4(myBlockDim.x) +
#ifndef COORD_IN_ONE_VEC
      sizeof(ScalorType) * 3 * hroundUp4(myBlockDim.x) +
#else
      sizeof(CoordType) * hroundUp4(myBlockDim.x) +
#endif
      sizeof(TypeType)  *      hroundUp4(myBlockDim.x) +
      sizeof(ForceIndexType) * hroundUp4(nbForceTableLength);
  if (buildDeviceNeighborList_DeviceCellList_sbuffSize >=
      SystemSharedBuffSize - GlobalFunctionParamSizeLimit){
    sharednbForceTable = false;
    buildDeviceNeighborList_DeviceCellList_sbuffSize =
	sizeof(IndexType) *      hroundUp4(myBlockDim.x) +
#ifndef COORD_IN_ONE_VEC
	sizeof(ScalorType) * 3 * hroundUp4(myBlockDim.x) +
#else
	sizeof(CoordType) * hroundUp4(myBlockDim.x) +
#endif
	sizeof(TypeType)  *      hroundUp4(myBlockDim.x);
  }  
  buildDeviceNeighborList_AllPair_sbuffSize =
      buildDeviceNeighborList_DeviceCellList_sbuffSize;
  buildDeviceCellList_step1_sbuffSize = 2 * myBlockDim.x * sizeof(IndexType);
  judgeRebuild_judgeCoord_block_sbuffSize = myBlockDim.x * sizeof(IndexType);
}

void NeighborList::build (const MDSystem & sys,
			  MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeNormalizeSys);
  normalizeSystem
      <<<atomGridDim, myBlockDim>>> (
	  sys.box, sys.ddata.numAtom,
#ifndef COORD_IN_ONE_VEC
	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
#else
	  sys.ddata.coord,
#endif
	  sys.ddata.coordNoix, sys.ddata.coordNoiy, sys.ddata.coordNoiz);
  checkCUDAError ("NeighborList::build, normalize System");
  if (timer != NULL) timer->toc(mdTimeNormalizeSys);
  
  if (mode == CellListBuilt){
    // buildDeviceNeighborList_AllPair 
    // 	<<<atomGridDim, myBlockDim>>>(
    // 	    sys.ddata.numAtom,
    // 	    sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
    // 	    sys.box, dnlist);
    // checkCUDAError ("NeighborList::build, build neighbor list all pair");
    
    if (timer != NULL) timer->tic(mdTimeBuildCellList);
    prepare_naivlyBuildDeviceCellList
    	<<<cellGridDim, myBlockDim>>> (dclist);
    naivlyBuildDeviceCellList2        
    	<<<atomGridDim, myBlockDim,
	buildDeviceCellList_step1_sbuffSize>>> (
    	    sys.ddata.numAtom,
#ifndef COORD_IN_ONE_VEC
    	    sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz, 
#else
	    sys.ddata.coord,
#endif
	    sys.ddata.coordNoix,
	    sys.ddata.coordNoiy,
	    sys.ddata.coordNoiz,
    	    sys.box, dclist,
	    err.ptr_de,
	    err.ptr_dindex, err.ptr_dscalor);
    err.check ("NeighborList::build, naivly build cell list");
    checkCUDAError ("NeighborList::build, naivly build cell list");
    if (timer != NULL) timer->toc(mdTimeBuildCellList);
    if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
    buildDeviceNeighborList_DeviceCellList 
    	<<<cellGridDim, myBlockDim,
	buildDeviceNeighborList_DeviceCellList_sbuffSize>>> (
	    // <<<cellGridDim, myBlockDim>>> (
    	    sys.ddata.numAtom, 
#ifndef COORD_IN_ONE_VEC
    	    sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
#else
	    sys.ddata.coord,
#endif
	    sys.ddata.type,
	    sys.box, dclist, dnlist,
	    nbForceTable, NatomType,
	    sharednbForceTable,
	    err.ptr_de);
    err.check("NeighborList::build, build neighbor list from cell list");
    checkCUDAError ("NeighborList::build, build neighbor list from cell list");
    if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
  }
  else if (mode == AllPairBuilt){
    if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
    buildDeviceNeighborList_AllPair 
	<<<atomGridDim, myBlockDim,
	buildDeviceNeighborList_AllPair_sbuffSize>>>(
	    sys.ddata.numAtom,
#ifndef COORD_IN_ONE_VEC
	    sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
#else
	    sys.ddata.coord,
#endif
	    sys.ddata.type,
	    sys.box, dnlist,
	    nbForceTable, NatomType,
	    sharednbForceTable,
	    err.ptr_de);
    err.check("NeighborList::build, build neighbor list all pair");
    checkCUDAError ("NeighborList::build, build neighbor list all pair");
    if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
  }
}


void NeighborList::reBuild (const MDSystem & sys,
			    MDTimer * timer)
{
  if (timer != NULL) timer->tic(mdTimeNormalizeSys);
  normalizeSystem
      <<<atomGridDim, myBlockDim>>> (
	  sys.box, sys.ddata.numAtom,
#ifndef COORD_IN_ONE_VEC
	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
#else
	  sys.ddata.coord,
#endif
	  sys.ddata.coordNoix, sys.ddata.coordNoiy, sys.ddata.coordNoiz);
  checkCUDAError ("NeighborList::rebuild, normalize System");
  if (timer != NULL) timer->toc(mdTimeNormalizeSys);

  if (mode == CellListBuilt){
    if (timer != NULL) timer->tic(mdTimeBuildCellList);
    ScalorType rlisti = 1./dclist.rlist;
    bool CellOnX = mybdir & RectangularBoxGeometry::mdRectBoxDirectionX;
    bool CellOnY = mybdir & RectangularBoxGeometry::mdRectBoxDirectionY;
    bool CellOnZ = mybdir & RectangularBoxGeometry::mdRectBoxDirectionZ;
    if ((CellOnX && fabsf(floorf(sys.box.size.x * rlisti) - dclist.NCell.x) > 0.1) ||
	(CellOnY && fabsf(floorf(sys.box.size.y * rlisti) - dclist.NCell.y) > 0.1) ||
	(CellOnZ && fabsf(floorf(sys.box.size.z * rlisti) - dclist.NCell.z) > 0.1) ){
    // if (true){
      printf ("### box changes too much naively rebuild\n");
      reinitCellList (sys, dclist.rlist, mybdir);
      prepare_naivlyBuildDeviceCellList
	  <<<cellGridDim, myBlockDim>>> (dclist);
      naivlyBuildDeviceCellList2        
	  <<<atomGridDim, myBlockDim,
	  buildDeviceCellList_step1_sbuffSize>>> (
	      sys.ddata.numAtom,
#ifndef COORD_IN_ONE_VEC
	      sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz, 
#else
	      sys.ddata.coord,
#endif
	      sys.ddata.coordNoix,
	      sys.ddata.coordNoiy,
	      sys.ddata.coordNoiz,
	      sys.box, dclist,
	      err.ptr_de,
	      err.ptr_dindex, err.ptr_dscalor);
      err.check ("NeighborList::rebuild, naivly build cell list");
      checkCUDAError ("NeighborList::rebuild, naivly build cell list");
    }
    else {
      buildDeviceCellList_clearBuff
	  <<<cellGridDim, myBlockDim>>> (
	      mySendBuff);
      err.check ("NeighborList::rebuild, clear buff");
      buildDeviceCellList_step1 
	  <<<cellGridDim, myBlockDim,
	  buildDeviceCellList_step1_sbuffSize>>> (
	      sys.ddata.numAtom, 
#ifndef COORD_IN_ONE_VEC
	      sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
#else
	      sys.ddata.coord,
#endif
	      sys.ddata.coordNoix, sys.ddata.coordNoiy, sys.ddata.coordNoiz,
	      sys.box, dclist, mySendBuff, myTargetBuff,
	      err.ptr_de, err.ptr_dindex, err.ptr_dscalor);
      err.check ("NeighborList::rebuild, build cell list step1");
      checkCUDAError ("NeighborList::rebuild, build cell list step1");
      buildDeviceCellList_step2
	  <<<cellGridDim, myBlockDim>>> (
	      sys.box, dclist, mySendBuff, myTargetBuff, bitDeepth,
	      err.ptr_de);
      err.check ("NeighborList::rebuild, build cell list step2");
    }
    checkCUDAError ("NeighborList::rebuild, build cell list step2");
    if (timer != NULL) timer->toc(mdTimeBuildCellList);
    if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
    buildDeviceNeighborList_DeviceCellList 
    	<<<cellGridDim, myBlockDim,
    	buildDeviceNeighborList_DeviceCellList_sbuffSize>>> (
	    // <<<cellGridDim, myBlockDim>>>(
    	    sys.ddata.numAtom, 
#ifndef COORD_IN_ONE_VEC
	    sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
#else
	    sys.ddata.coord,
#endif
    	    sys.ddata.type,
    	    sys.box, dclist, dnlist,
    	    nbForceTable, NatomType,
    	    sharednbForceTable,
    	    err.ptr_de);
    err.check ("NeighborList::rebuild, build neighbor list from cell list");
    checkCUDAError ("NeighborList::rebuild, build neighbor list from cell list");
    if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
  }
  else if (mode == AllPairBuilt){
    if (timer != NULL) timer->tic(mdTimeBuildNeighborList);
    buildDeviceNeighborList_AllPair 
	<<<atomGridDim, myBlockDim,
	buildDeviceNeighborList_AllPair_sbuffSize>>>(
	    sys.ddata.numAtom,
#ifndef COORD_IN_ONE_VEC
	    sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
#else
	    sys.ddata.coord,
#endif
	    sys.ddata.type,
	    sys.box, dnlist,
	    nbForceTable, NatomType,
	    sharednbForceTable,
	    err.ptr_de);
    err.check ("NeighborList::rebuild, build neighbor list all pair");
    checkCUDAError ("NeighborList::rebuild, build neighbor list all pair");
    if (timer != NULL) timer->toc(mdTimeBuildNeighborList);
  }
}


bool NeighborList::judgeRebuild (const MDSystem & sys,
				 const ScalorType &diffTol,
				 MDTimer *timer)
{
  if (timer != NULL) timer->tic(mdTimeJudgeRebuild);
  IndexType mytag;
  ScalorType diffTol2 = diffTol * diffTol;
  
  // judgeRebuild_judgeCoord
  //     <<<atomGridDim, myBlockDim>>> (
  // 	  sys.ddata.numAtom, sys.box,
  // 	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
  // 	  backupCoordx, backupCoordy, backupCoordz,
  // 	  diffTol2,
  // 	  judgeRebuild_buff, judgeRebuild_tag);
  judgeRebuild_judgeCoord_block 
      <<<atomGridDim, myBlockDim,
      judgeRebuild_judgeCoord_block_sbuffSize>>> (
  	  sys.ddata.numAtom, sys.box,
#ifndef COORD_IN_ONE_VEC
  	  sys.ddata.coordx, sys.ddata.coordy, sys.ddata.coordz,
  	  backupCoordx, backupCoordy, backupCoordz,
#else
	  sys.ddata.coord,
	  backupCoord,
#endif
  	  diffTol2,
  	  sum_judge.getBuff());
  sum_judge.sumBuff (judgeRebuild_tag, 0);	  
  cudaMemcpy (&mytag, judgeRebuild_tag, sizeof(IndexType), cudaMemcpyDeviceToHost);
  if (mytag != 0){
#ifndef COORD_IN_ONE_VEC
    judgeRebuild_backupCoord
	<<<atomGridDim, myBlockDim>>> (
	    sys.ddata.numAtom, sys.ddata.coordx, backupCoordx);
    judgeRebuild_backupCoord
	<<<atomGridDim, myBlockDim>>> (
	    sys.ddata.numAtom, sys.ddata.coordy, backupCoordy);
    judgeRebuild_backupCoord
	<<<atomGridDim, myBlockDim>>> (
	    sys.ddata.numAtom, sys.ddata.coordz, backupCoordz);
#else
    cpyProperty <<<atomGridDim, myBlockDim>>> (
	backupCoord, sys.ddata.coord, sys.ddata.numAtom);
#endif    
    if (timer != NULL) timer->toc(mdTimeJudgeRebuild);
    // fprintf (stderr, "rebuild!\n");
    return true;
  }
  else{
    if (timer != NULL) timer->toc(mdTimeJudgeRebuild);
    return false;
  }
}

    


////////////////////////////////////////////////////////////
// for the reason of using texture, we place this function here. it
// should be placed in NeighborList.cu
////////////////////////////////////////////////////////////


__device__ IndexType shiftedD3toD1 (
    DeviceCellList clist,
    RectangularBoxGeometry::RectangularBox box,
    int ix, int iy, int iz,
    ScalorType * shiftx , ScalorType * shifty, ScalorType * shiftz)
{
  int tmp;
  ix += (tmp = -int(floorf(ix * clist.NCelli.x))) * clist.NCell.x;
  *shiftx = tmp * box.size.x;
  iy += (tmp = -int(floorf(iy * clist.NCelli.y))) * clist.NCell.y;
  *shifty = tmp * box.size.y;
  iz += (tmp = -int(floorf(iz * clist.NCelli.z))) * clist.NCell.z;
  *shiftz = tmp * box.size.z;
  return D3toD1 (clist, ix, iy, iz);
}


#ifndef COORD_IN_ONE_VEC
__global__ void buildDeviceNeighborList_DeviceCellList (
    IndexType numAtom,
    ScalorType * coordx, ScalorType * coordy, ScalorType * coordz,
    TypeType * type,
    RectangularBox box,
    DeviceCellList clist,
    DeviceNeighborList nlist,
    ForceIndexType *nbForceTable,
    IndexType NatomType,
    bool sharednbForceTable,
    mdError_t * ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType bidx, bidy, bidz;
  D1toD3 (clist, bid, bidx, bidy, bidz);
  
  // set number of neighbor to 0
  IndexType Nneighbor = 0;
  // load index
  IndexType ii = getDeviceCellListData (clist, bid, tid);
  // load iith coordinate // use texturefetch instead
  ScalorType refx, refy, refz;
  TypeType reftype;
  if (ii != MaxIndexValue){
#ifdef COMPILE_NO_TEX
    refx = coordx[ii];
    refy = coordy[ii];
    refz = coordz[ii];
    reftype = type[ii];
#else
    refx = tex1Dfetch(global_texRef_neighbor_coordx, ii);
    refy = tex1Dfetch(global_texRef_neighbor_coordy, ii);
    refz = tex1Dfetch(global_texRef_neighbor_coordz, ii);
    reftype = tex1Dfetch(global_texRef_neighbor_type, ii);
#endif
  }
  ScalorType rlist = clist.rlist;

  // the target index and coordinates are shared

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
  
  // __shared__ volatile  IndexType  targetIndexes [MaxThreadsPerBlock];
  // __shared__ volatile  ScalorType targetx       [MaxThreadsPerBlock];
  // __shared__ volatile  ScalorType targety       [MaxThreadsPerBlock];
  // __shared__ volatile  ScalorType targetz       [MaxThreadsPerBlock];
  // __shared__ volatile  TypeType   targettype    [MaxThreadsPerBlock];
  // __shared__ volatile  ForceIndexType nbForceTableBuff [MaxNBForceTableBuffSize];

  // IndexType nbForceTableLength = AtomNBForceTable::dCalDataLength(NatomType);
  // if (sharednbForceTable){
  //   cpyGlobalDataToSharedBuff (nbForceTable, nbForceTableBuff, nbForceTableLength);
  // }
  // __syncthreads();
  
  int upperx, lowerx;
  int uppery, lowery;
  int upperz, lowerz;
  (clist.NCell.x != 1) ? (upperx =  1) : (upperx = 0);
  (clist.NCell.x != 1) ? (lowerx = -1) : (lowerx = 0);
  (clist.NCell.y != 1) ? (uppery =  1) : (uppery = 0);
  (clist.NCell.y != 1) ? (lowery = -1) : (lowery = 0);
  (clist.NCell.z != 1) ? (upperz =  1) : (upperz = 0);
  (clist.NCell.z != 1) ? (lowerz = -1) : (lowerz = 0);
  
  // loop over 27 neighbor cells
  for (int di = lowerx; di <= upperx; ++di){
    for (int dj = lowery; dj <= uppery; ++dj){
      for (int dk = lowerz; dk <= upperz; ++dk){
  // for (int di = -1; di <= 1; ++di){
  //   for (int dj = -1; dj <= 1; ++dj){
  //     for (int dk = -1; dk <= 1; ++dk){
	__syncthreads();
	// the shift value of a cell is pre-computed
	ScalorType xshift, yshift, zshift;
	int nci = di + bidx;
	int ncj = dj + bidy;
	int nck = dk + bidz;
	IndexType targetCellIdx = shiftedD3toD1 (clist, box, 
						 nci, ncj, nck, 
						 &xshift, &yshift, &zshift);
	// load target index and coordinates form global memary
	IndexType tmp = (targetIndexes[tid] = 
			 getDeviceCellListData(clist, targetCellIdx, tid));
	if (tmp != MaxIndexValue){
#ifdef COMPILE_NO_TEX
	  targetx[tid] = coordx[tmp];
	  targety[tid] = coordy[tmp];
	  targetz[tid] = coordz[tmp];
	  targettype[tid] = type[tmp];
#else
	  targetx[tid] = tex1Dfetch(global_texRef_neighbor_coordx, tmp);
	  targety[tid] = tex1Dfetch(global_texRef_neighbor_coordy, tmp);
	  targetz[tid] = tex1Dfetch(global_texRef_neighbor_coordz, tmp);
	  targettype[tid] = tex1Dfetch(global_texRef_neighbor_type, tmp);
#endif
	}
	__syncthreads();
	// find neighbor
	if (ii != MaxIndexValue){
	  for (IndexType jj = 0; jj < blockDim.x; ++jj){
	    if (targetIndexes[jj] == MaxIndexValue) break;
	    ScalorType diffx = targetx[jj] - xshift - refx;
	    ScalorType diffy = targety[jj] - yshift - refy;
	    ScalorType diffz = targetz[jj] - zshift - refz;
	    if (clist.NCell.x == 1) shortestImage (box.size.x, box.sizei.x, &diffx);
	    if (clist.NCell.y == 1) shortestImage (box.size.y, box.sizei.y, &diffy);
	    if (clist.NCell.z == 1) shortestImage (box.size.z, box.sizei.z, &diffz);
	    if ((diffx*diffx+diffy*diffy+diffz*diffz) < rlist*rlist &&
		targetIndexes[jj] != ii){
	      ForceIndexType fidx;
	      if (sharednbForceTable){
		fidx = AtomNBForceTable::calForceIndex (
		    nbForceTableBuff, NatomType, reftype, targettype[jj]);
	      }
	      else {
		fidx = AtomNBForceTable::calForceIndex (
		    nbForceTable, NatomType, reftype, targettype[jj]);
	      }
	      if (fidx != mdForceNBNull) {
		IndexType listIdx = Nneighbor * nlist.stride + ii;
		nlist.data[listIdx] = targetIndexes[jj];
		nlist.forceIndex[listIdx] = fidx;
		Nneighbor ++;
	      }
	    }
	  }
	}
      }
    }
  }
  if (ii != MaxIndexValue) {
    nlist.Nneighbor[ii] = Nneighbor;
    if (Nneighbor > nlist.listLength && ptr_de != NULL){
      *ptr_de = mdErrorShortNeighborList;
    }
  }
}



__global__ void buildDeviceCellList_step1 (IndexType numAtom,
					   ScalorType * coordx,
					   ScalorType * coordy,
					   ScalorType * coordz,
					   IntScalorType * coordNoix,
					   IntScalorType * coordNoiy,
					   IntScalorType * coordNoiz,
					   RectangularBox box,
					   DeviceCellList clist,
					   IndexType * sendBuff,
					   IndexType * targetBuff,
					   mdError_t * ptr_de,
					   IndexType * erridx,
					   ScalorType * errsrc)
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
    ScalorType refx, refy, refz;
#ifdef COMPILE_NO_TEX
    refx = coordx[thisid];
    refy = coordy[thisid];
    refz = coordz[thisid];
#else
    refx = tex1Dfetch(global_texRef_neighbor_coordx, thisid);
    refy = tex1Dfetch(global_texRef_neighbor_coordy, thisid);
    refz = tex1Dfetch(global_texRef_neighbor_coordz, thisid);
#endif
    // if (box.size.x - refx < box.size.x * 1e-5){
    //   refx -= box.size.x * 1e-5;
    // }
    // if (box.size.y - refy < box.size.y * 1e-5){
    //   refy -= box.size.y * 1e-5;
    // }
    // if (box.size.z - refz < box.size.z * 1e-5){
    //   refz -= box.size.z * 1e-5;      
    // }
    targetCelli = IndexType(refx * box.sizei.x * ScalorType (clist.NCell.x));
    targetCellj = IndexType(refy * box.sizei.y * ScalorType (clist.NCell.y));
    targetCellk = IndexType(refz * box.sizei.z * ScalorType (clist.NCell.z));
    if (targetCelli == clist.NCell.x){
      targetCelli -= clist.NCell.x;
      coordx[thisid] -= box.size.x;
      coordNoix[thisid] ++;
    }
    if (targetCellj == clist.NCell.y){
      targetCellj -= clist.NCell.y;
      coordy[thisid] -= box.size.y;
      coordNoiy[thisid] ++;
    }
    if (targetCellk == clist.NCell.z){
      targetCellk -= clist.NCell.z;
      coordz[thisid] -= box.size.z;
      coordNoiz[thisid] ++;
    }
    targetCellid[tid] = D3toD1 (clist, targetCelli, targetCellj, targetCellk);
    if (ptr_de != NULL && 
	(targetCelli >= clist.NCell.x || 
	 targetCellj >= clist.NCell.y || 
	 targetCellk >= clist.NCell.z)){
      *ptr_de = mdErrorOverFlowCellIdx;
      if (targetCelli >= IndexType(clist.NCell.x)){
	*erridx = targetCelli;
	*errsrc = refx;
	return;
      }
      if (targetCellj >= IndexType(clist.NCell.y)){
	*erridx = targetCellj;
	*errsrc = refy;
	return;
      }
      if (targetCellk >= IndexType(clist.NCell.z)){
	*erridx = targetCellk;
	*errsrc = refz;
	return;
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
  // if (tid == 0 && originalNumber - total0 != 0){
  //   printf ("bid: %d, Nsend is %d, targetid is %d, targetCell is %d\n", 
  // 	    bid, originalNumber - total0, 
  // 	    originalData[total0] - mark,
  // 	    targetCellid[total0]);
  // }
  
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



__global__ void buildDeviceCellList_initBuff (IndexType * sendBuff,
					      IndexType * targetBuff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  sendBuff[ii] = MaxIndexValue;
  targetBuff[ii] = 0;
}
__global__ void buildDeviceCellList_clearBuff (IndexType * sendBuff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  sendBuff[ii] = MaxIndexValue;
}

  
  
__global__ void buildDeviceCellList_step2 (RectangularBox box,
					   DeviceCellList clist,
					   IndexType * sendBuff,
					   IndexType * targetBuff,
					   IndexType bitDeepth,
					   mdError_t * ptr_de)
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

// __global__ void buildDeviceCellList_step2 (RectangularBox box,
// 					   DeviceCellList clist,
// 					   IndexType * sendBuff,
// 					   IndexType * targetBuff,
// 					   IndexType bitDeepth,
// 					   mdError_t * ptr_de)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   IndexType bidx, bidy, bidz;
//   D1toD3 (clist, bid, bidx, bidy, bidz);

//   __shared__ IndexType myAtom [MaxThreadsPerBlock];
//   __shared__ IndexType neighborSendIdx[MaxThreadsPerBlock];
//   __shared__ IndexType neighborTarget [MaxThreadsPerBlock];
//   __shared__ IndexType myNum;
//   myAtom[tid] = clist.data[bid*clist.stride + tid];
//   if (tid == 0) myNum = clist.numbers[bid];
  
//   for (int di = -1; di <= 1; ++di){
//     for (int dj = -1; dj <=1; ++dj){
//       for (int dk = -1; dk <=1; ++dk){
// 	if (di == 0 && dj == 0 && dk == 0) continue;
// 	__syncthreads();
// 	ScalorType xshift, yshift, zshift;
// 	int nci = di + bidx;
// 	int ncj = dj + bidy;
// 	int nck = dk + bidz;
// 	IndexType neighborCellIdx = shiftedD3toD1 (clist, box, 
// 						   nci, ncj, nck, 
// 						   &xshift, &yshift, &zshift);
// 	neighborSendIdx[tid] = sendBuff  [neighborCellIdx*clist.stride + tid];
// 	neighborTarget [tid] = targetBuff[neighborCellIdx*clist.stride + tid];
// 	__syncthreads();
// 	if (tid == 0){
// 	  IndexType ii = 0;
// 	  while (neighborSendIdx[ii] != MaxIndexValue){
// 	    if (neighborTarget[ii] == bid){
// 	      myAtom[myNum++] = neighborSendIdx[ii];
// 	      // printf("copy from %d, to %d, atom %d\n", 
// 	      // 	     neighborCellIdx, bid, neighborSendIdx[ii]);
// 	    }
// 	    ++ii;
// 	  }
// 	}
//       }
//     }
//   }
//   __syncthreads();
// //  sortList (myAtom, bitDeepth);
//   clist.data[bid*clist.stride + tid] = myAtom[tid];
//   if (tid == 0) {
//     clist.numbers[bid] = myNum;
//     if (myNum > blockDim.x && ptr_de != NULL) {
//       *ptr_de = mdErrorShortCellList;
//     }
//   }
// }




__global__ void judgeRebuild_backupCoord (const IndexType numAtom,
					  const ScalorType * from,
					  ScalorType * to)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    to[ii] = from[ii];
  }
}


__device__ IndexType judgeRebuild_counter = 0;
__global__ void judgeRebuild_judgeCoord (const IndexType numAtom,
					 const RectangularBox box,
					 const ScalorType* coordx,
					 const ScalorType* coordy,
					 const ScalorType* coordz,
					 const ScalorType* backupCoordx,
					 const ScalorType* backupCoordy,
					 const ScalorType* backupCoordz,
					 const ScalorType diffTol2,
					 ScalorType * judgeRebuild_buff,
					 IndexType * tag)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  __shared__ volatile ScalorType mydiff2[MaxThreadsPerBlock * 2];

  ScalorType dx=0.f, dy=0.f, dz=0.f;
  if (ii < numAtom){
    dx = coordx[ii] - backupCoordx[ii];
    dy = coordy[ii] - backupCoordy[ii];
    dz = coordz[ii] - backupCoordz[ii];
    shortestImage (box, &dx, &dy, &dz);
  }
  mydiff2[threadIdx.x] = dx*dx + dy*dy + dz*dz;
  mydiff2[threadIdx.x + blockDim.x] = 0.f;
  __syncthreads();
  IndexType num = ((bid+1) * blockDim.x > numAtom) ? 
      (blockDim.x - ((bid+1) * blockDim.x) + numAtom) : blockDim.x;
  ScalorType partialMax = maxVectorBlockBuffer (mydiff2, num);
  ScalorType maxDiff2;
  
  if (maxPartialMax (partialMax, judgeRebuild_buff, &judgeRebuild_counter, &maxDiff2)){
    if (maxDiff2 > diffTol2){
      *tag = 1;
    }
    else {
      *tag = 0;
    }
  }
}


__global__ void judgeRebuild_judgeCoord_block (const IndexType numAtom,
					       const RectangularBox box,
					       const ScalorType* coordx,
					       const ScalorType* coordy,
					       const ScalorType* coordz,
					       const ScalorType* backupCoordx,
					       const ScalorType* backupCoordy,
					       const ScalorType* backupCoordz,
					       const ScalorType diffTol2,
					       IndexType * judgeRebuild_buff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  extern __shared__ volatile IndexType sbuff[];  
  // __shared__ volatile IndexType sbuff [MaxThreadsPerBlock];

  ScalorType dx=0.f, dy=0.f, dz=0.f;
  if (ii < numAtom){
    dx = coordx[ii] - backupCoordx[ii];
    dy = coordy[ii] - backupCoordy[ii];
    dz = coordz[ii] - backupCoordz[ii];
    shortestImage (box, &dx, &dy, &dz);
  }

  if (ii < numAtom){
    sbuff[threadIdx.x] = ((dx*dx + dy*dy + dz*dz) > diffTol2);
  }
  else {
    sbuff[threadIdx.x] = 0;
  }
  __syncthreads();
  sumVectorBlockBuffer_2 (sbuff);
  if (threadIdx.x == 0) judgeRebuild_buff[bid] = (sbuff[0] != 0);
}



//////////////////////////////////////////////////
// coord in one vec
//////////////////////////////////////////////////

#else
__global__ void buildDeviceNeighborList_DeviceCellList (
    IndexType numAtom,
    CoordType * coord,
    TypeType * type,
    RectangularBox box,
    DeviceCellList clist,
    DeviceNeighborList nlist,
    ForceIndexType *nbForceTable,
    IndexType NatomType,
    bool sharednbForceTable,
    mdError_t * ptr_de)
{
  // RectangularBoxGeometry::normalizeSystem (box, &ddata);
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType bidx, bidy, bidz;
  D1toD3 (clist, bid, bidx, bidy, bidz);
  
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
  ScalorType rlist = clist.rlist;

  // the target index and coordinates are shared

  extern __shared__ volatile char pub_sbuff[];
  
  volatile IndexType * targetIndexes =
      (volatile IndexType *) pub_sbuff;
  CoordType * target =
      (CoordType *) &targetIndexes[roundUp4(blockDim.x)];
  volatile TypeType * targettype =
      (volatile TypeType *) &target[roundUp4(blockDim.x)];
  volatile ForceIndexType * nbForceTableBuff = NULL;

  IndexType nbForceTableLength = AtomNBForceTable::dCalDataLength(NatomType);
  if (sharednbForceTable){
    nbForceTableBuff = (volatile ForceIndexType *) &targettype[roundUp4(blockDim.x)];
    cpyGlobalDataToSharedBuff (nbForceTable, nbForceTableBuff, nbForceTableLength);
  }
  __syncthreads();

  // __shared__ volatile  IndexType  targetIndexes [MaxThreadsPerBlock];
  // __shared__ volatile CoordType target [MaxThreadsPerBlock];
  // __shared__ volatile  TypeType   targettype    [MaxThreadsPerBlock];
  // __shared__ volatile  ForceIndexType nbForceTableBuff [MaxNBForceTableBuffSize];

  // IndexType nbForceTableLength = AtomNBForceTable::dCalDataLength(NatomType);
  // if (sharednbForceTable){
  //   cpyGlobalDataToSharedBuff (nbForceTable, nbForceTableBuff, nbForceTableLength);
  // }
  // __syncthreads();

  bool oneCellX(false), oneCellY(false), oneCellZ(false);
  if (clist.NCell.x == 1) oneCellX = true;
  if (clist.NCell.y == 1) oneCellY = true;
  if (clist.NCell.z == 1) oneCellZ = true;
  int upperx(1), lowerx(-1);
  int uppery(1), lowery(-1);
  int upperz(1), lowerz(-1);
  if (oneCellX) {lowerx =  0; upperx = 0;}
  if (oneCellY) {lowery =  0; uppery = 0;}
  if (oneCellZ) {lowerz =  0; upperz = 0;}
  ScalorType rlist2 = rlist * rlist;
  
  // loop over 27 neighbor cells
  #pragma unroll 3
  for (int di = lowerx; di <= upperx; ++di){
    for (int dj = lowery; dj <= uppery; ++dj){
      for (int dk = lowerz; dk <= upperz; ++dk){
	__syncthreads();
	// the shift value of a cell is pre-computed
	ScalorType xshift, yshift, zshift;
	int nci = di + bidx;
	int ncj = dj + bidy;
	int nck = dk + bidz;
	IndexType targetCellIdx = shiftedD3toD1 (clist, box, 
						 nci, ncj, nck, 
						 &xshift, &yshift, &zshift);
	// load target index and coordinates form global memary
	IndexType tmp = (targetIndexes[tid] = 
			 getDeviceCellListData(clist, targetCellIdx, tid));
	if (tmp != MaxIndexValue){
#ifdef COMPILE_NO_TEX
	  target[tid] = coord[tmp];
	  targettype[tid] = type[tmp];
#else
	  target[tid] = tex1Dfetch(global_texRef_neighbor_coord, tmp);
	  targettype[tid] = tex1Dfetch(global_texRef_neighbor_type, tmp);
#endif
	}
	__syncthreads();
	// find neighbor
	if (ii != MaxIndexValue){
	  for (IndexType jj = 0; jj < blockDim.x; ++jj){
	    if (targetIndexes[jj] == MaxIndexValue) break;
	    ScalorType diffx = target[jj].x - xshift - ref.x;
	    ScalorType diffy = target[jj].y - yshift - ref.y;
	    ScalorType diffz = target[jj].z - zshift - ref.z;
	    if (oneCellX) shortestImage (box.size.x, box.sizei.x, &diffx);
	    if (oneCellY) shortestImage (box.size.y, box.sizei.y, &diffy);
	    if (oneCellZ) shortestImage (box.size.z, box.sizei.z, &diffz);
	    if ((diffx*diffx+diffy*diffy+diffz*diffz) < rlist2 &&
		targetIndexes[jj] != ii){
	      ForceIndexType fidx;
	      if (sharednbForceTable){
		fidx = AtomNBForceTable::calForceIndex (
		    nbForceTableBuff, NatomType, reftype, targettype[jj]);
	      }
	      else {
		fidx = AtomNBForceTable::calForceIndex (
		    nbForceTable, NatomType, reftype, targettype[jj]);
	      }
	      if (fidx != mdForceNULL) {
		IndexType listIdx = Nneighbor * nlist.stride + ii;
		nlist.data[listIdx] = targetIndexes[jj];
		nlist.forceIndex[listIdx] = fidx;
		Nneighbor ++;
	      }
	    }
	  }
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
  }
}



__global__ void buildDeviceCellList_step1 (IndexType numAtom,
					   CoordType * coord,
					   IntScalorType * coordNoix,
					   IntScalorType * coordNoiy,
					   IntScalorType * coordNoiz,
					   RectangularBox box,
					   DeviceCellList clist,
					   IndexType * sendBuff,
					   IndexType * targetBuff,
					   mdError_t * ptr_de,
					   IndexType * erridx,
					   ScalorType * errsrc)
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
    targetCellid[tid] = D3toD1 (clist, targetCelli, targetCellj, targetCellk);
    if (ptr_de != NULL && 
	(targetCelli >= clist.NCell.x || 
	 targetCellj >= clist.NCell.y || 
	 targetCellk >= clist.NCell.z)){
      *ptr_de = mdErrorOverFlowCellIdx;
      if (targetCelli >= IndexType(clist.NCell.x)){
	*erridx = targetCelli;
	*errsrc = ref.x;
	return;
      }
      if (targetCellj >= IndexType(clist.NCell.y)){
	*erridx = targetCellj;
	*errsrc = ref.y;
	return;
      }
      if (targetCellk >= IndexType(clist.NCell.z)){
	*erridx = targetCellk;
	*errsrc = ref.z;
	return;
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



__global__ void buildDeviceCellList_initBuff (IndexType * sendBuff,
					      IndexType * targetBuff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  sendBuff[ii] = MaxIndexValue;
  targetBuff[ii] = 0;
}
__global__ void buildDeviceCellList_clearBuff (IndexType * sendBuff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  sendBuff[ii] = MaxIndexValue;
}

  
  
__global__ void buildDeviceCellList_step2 (RectangularBox box,
					   DeviceCellList clist,
					   IndexType * sendBuff,
					   IndexType * targetBuff,
					   IndexType bitDeepth,
					   mdError_t * ptr_de)
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

__global__ void judgeRebuild_backupCoord (const IndexType numAtom,
					  const ScalorType * from,
					  ScalorType * to)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;
  if (ii < numAtom){
    to[ii] = from[ii];
  }
}


__device__ IndexType judgeRebuild_counter = 0;
__global__ void judgeRebuild_judgeCoord (const IndexType numAtom,
					 const RectangularBox box,
					 const CoordType * coord,
					 const CoordType * backupCoord,
					 const ScalorType diffTol2,
					 ScalorType * judgeRebuild_buff,
					 IndexType * tag)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  __shared__ volatile ScalorType mydiff2[MaxThreadsPerBlock * 2];

  ScalorType dx=0.f, dy=0.f, dz=0.f;
  if (ii < numAtom){
    dx = coord[ii].x - backupCoord[ii].x;
    dy = coord[ii].y - backupCoord[ii].y;
    dz = coord[ii].z - backupCoord[ii].z;
    shortestImage (box, &dx, &dy, &dz);
  }
  mydiff2[threadIdx.x] = dx*dx + dy*dy + dz*dz;
  mydiff2[threadIdx.x + blockDim.x] = 0.f;
  __syncthreads();
  IndexType num = ((bid+1) * blockDim.x > numAtom) ? 
      (blockDim.x - ((bid+1) * blockDim.x) + numAtom) : blockDim.x;
  ScalorType partialMax = maxVectorBlockBuffer (mydiff2, num);
  ScalorType maxDiff2;
  
  if (maxPartialMax (partialMax, judgeRebuild_buff, &judgeRebuild_counter, &maxDiff2)){
    if (maxDiff2 > diffTol2){
      *tag = 1;
    }
    else {
      *tag = 0;
    }
  }
}


__global__ void judgeRebuild_judgeCoord_block (const IndexType numAtom,
					       const RectangularBox box,
					       const CoordType * coord,
					       const CoordType * backupCoord,
					       const ScalorType diffTol2,
					       IndexType * judgeRebuild_buff)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType ii = threadIdx.x + bid * blockDim.x;

  extern __shared__ volatile IndexType sbuff[];  
  // __shared__ volatile IndexType sbuff [MaxThreadsPerBlock];

  ScalorType dx=0.f, dy=0.f, dz=0.f;
  if (ii < numAtom){
    dx = coord[ii].x - backupCoord[ii].x;
    dy = coord[ii].y - backupCoord[ii].y;
    dz = coord[ii].z - backupCoord[ii].z;
    shortestImage (box, &dx, &dy, &dz);
  }

  if (ii < numAtom){
    sbuff[threadIdx.x] = ((dx*dx + dy*dy + dz*dz) > diffTol2);
  }
  else {
    sbuff[threadIdx.x] = 0;
  }
  __syncthreads();
  sumVectorBlockBuffer_2 (sbuff);
  if (threadIdx.x == 0) judgeRebuild_buff[bid] = (sbuff[0] != 0);
}

#endif
