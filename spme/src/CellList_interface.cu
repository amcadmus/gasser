#define DEVICE_CODE

#include "NeighborList.h"
#include "CellList_interface.h"

void CellList::
mallocDeviceCellList (const IntVectorType & NCell,
		      const HostVectorType & boxSize)
{
  if (NCell.x == 1){
    dclist.NCell.x = 1;
    dclist.NCelli.x = 1./boxSize.x;
  }
  else{
    dclist.NCelli.x = 1./ (dclist.NCell.x = NCell.x);
  }
  if (NCell.y == 1){
    dclist.NCell.y = 1;
    dclist.NCelli.y = 1./boxSize.y;
  }
  else{
    dclist.NCelli.y = 1./ (dclist.NCell.y = NCell.y);
  }
  if (NCell.z == 1){
    dclist.NCell.z = 1;
    dclist.NCelli.z = 1./boxSize.z;
  }
  else{
    dclist.NCelli.z = 1./ (dclist.NCell.z = NCell.z);
  }

  IndexType numCell = dclist.NCell.x * dclist.NCell.y * dclist.NCell.z;
  dclist.rlist = mycellSize;
  dclist.divide = mydivide;
  // suppose the number of atoms in any cell is smaller or equal
  // to the number of threads in a block
  dclist.stride = cellBlockDim.x;
  cellGridDim = toGridDim (numCell);

  cudaMalloc ((void**)&(dclist.data), 
	      sizeof(IndexType) * numCell * dclist.stride);
  cudaMalloc ((void**)&(dclist.numbers), sizeof(IndexType) * numCell);
  cudaMalloc ((void**)&mySendBuff, 
	      sizeof(IndexType) * numCell * dclist.stride);
  cudaMalloc ((void**)&myTargetBuff, 
	      sizeof(IndexType) * numCell * dclist.stride);
  checkCUDAError ("CellList::init cell list");

  IndexType maxNumNeighborCell = (2*mydivide+1) * (2*mydivide+1) * (2*mydivide+1);
  dclist.maxNumNeighborCell = maxNumNeighborCell;
  cudaMalloc ((void**)&(dclist.numNeighborCell),
	      sizeof(IndexType) * numCell);
  cudaMalloc ((void**)&(dclist.neighborCellIndex),
	      sizeof(IndexType) * maxNumNeighborCell * numCell);
  cudaMalloc ((void**)&(dclist.neighborCellShiftNoi),
	      sizeof(CoordNoiType) * maxNumNeighborCell * numCell);
  checkCUDAError ("CellList::maxNumNeighborCell cell list buff");

  mallocedDeviceCellList = true;

  buildDeviceCellList_initBuff<<<cellGridDim, cellBlockDim>>> 
      (mySendBuff, myTargetBuff);
  checkCUDAError ("CellList::mallocedDeviceCellList cell list buff");
}

void CellList::
clearDeviceCellList () 
{
  if ( mallocedDeviceCellList ){
    cudaFree (dclist.data);
    cudaFree (dclist.numbers);
    cudaFree (dclist.numNeighborCell);
    cudaFree (dclist.neighborCellIndex);
    cudaFree (dclist.neighborCellShiftNoi);
    cudaFree (mySendBuff);
    cudaFree (myTargetBuff);
    mallocedDeviceCellList = false;
    checkCUDAError ("CellList::clearDeviceCellList");
  }
}

void CellList::
naivelyBuildDeviceCellList (const MDSystem & sys)
{  
  // buildDeviceCellList_clearBuff
  //     <<<cellGridDim, myBlockDim>>> (
  // 	  mySendBuff);
  // err.check ("CellList::rebuild, clear buff");
  prepare_naivelyBuildDeviceCellList
      <<<cellGridDim, cellBlockDim>>> (dclist);
  naivelyBuildDeviceCellList2        
      <<<atomGridDim, atomBlockDim,
      buildDeviceCellList_step1_sbuffSize>>> (
	  sys.ddata.numAtom,
	  sys.ddata.coord,
	  sys.ddata.coordNoix,
	  sys.ddata.coordNoiy,
	  sys.ddata.coordNoiz,
	  sys.box, dclist,
	  err.ptr_de,
	  err.ptr_dindex, err.ptr_dscalor);
  err.check ("CellList::naivelyBuildDeviceCellList");
  checkCUDAError ("CellList::naivelyBuildDeviceCellList");

  buildCellNeighborhood
      <<<cellGridDim, 1>>> (
	  dclist,
	  mydivide,
	  sys.box.size);
}

void CellList::
buildDeviceCellList (const MDSystem & sys)
{
  buildDeviceCellList_clearBuff
      <<<cellGridDim, cellBlockDim>>> (
	  mySendBuff);
  err.check ("CellList::rebuild, clear buff");
  size_t sbuffSize = sizeof (IndexType) * 2 * cellBlockDim.x;
  buildDeviceCellList_step1 
      <<<cellGridDim, cellBlockDim, sbuffSize>>> (
	  sys.ddata.numAtom, 
	  sys.ddata.coord,
	  sys.ddata.coordNoix,
	  sys.ddata.coordNoiy,
	  sys.ddata.coordNoiz,
	  sys.box,
	  dclist,
	  mySendBuff,
	  myTargetBuff,
	  err.ptr_de, err.ptr_dindex, err.ptr_dscalor);
  err.check ("CellList::buildDeviceCellList step 1");
  checkCUDAError ("CellList::buildDeviceCellList step 1");
  buildDeviceCellList_step2
      <<<cellGridDim, cellBlockDim>>> (
	  sys.box,
	  dclist,
	  mySendBuff,
	  myTargetBuff,
	  bitDeepth,
	  err.ptr_de);
  err.check ("CellList::buildDeviceCellList step2");
  checkCUDAError ("CellList::buildDeviceCellList step2");
}


bool CellList::
DecideNumCell (const MDSystem &		sys,
	       const ScalorType &	cellSize,
	       const BoxDirection_t &	bdir,
	       const IndexType &	divide,
	       IntVectorType &		NCell)
{
  bool CellOnX, CellOnY, CellOnZ;
  CellOnX = bdir & RectangularBoxGeometry::mdRectBoxDirectionX;
  CellOnY = bdir & RectangularBoxGeometry::mdRectBoxDirectionY;
  CellOnZ = bdir & RectangularBoxGeometry::mdRectBoxDirectionZ;
  double cellSizei = 1./cellSize;
  if (CellOnX ) NCell.x = int ( floor(sys.box.size.x * cellSizei) );
  else NCell.x = 1;
  if (CellOnY ) NCell.y = int ( floor(sys.box.size.y * cellSizei) );
  else NCell.y = 1;
  if (CellOnZ ) NCell.z = int ( floor(sys.box.size.z * cellSizei) );
  else NCell.z = 1;
  
  bool mode = true;
  if (CellOnX && NCell.x < 4) mode = false;
  if (CellOnY && NCell.y < 4) mode = false;
  if (CellOnZ && NCell.z < 4) mode = false;

  if (mode == true){
    if (CellOnX) NCell.x *= divide;
    if (CellOnY) NCell.y *= divide;
    if (CellOnZ) NCell.z *= divide;
  }
  return mode;
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

void CellList::
reinit (const MDSystem &		sys,
	const ScalorType &		cellSize,
	const IndexType &		NTreadCell,
	const IndexType &		NTreadAtom,
	const BoxDirection_t &		bdir,
	const IndexType &		divide)
{
  cellBlockDim.x = NTreadCell;
  atomBlockDim.x = NTreadAtom;
  bitDeepth = calDeepth(sys.ddata.numAtom);

  IndexType nob;
  nob = (sys.ddata.numAtom + atomBlockDim.x - 1) / atomBlockDim.x;
  atomGridDim = toGridDim (nob);

  mybdir = bdir;
  mycellSize = cellSize;
  mydivide = divide;

  buildDeviceCellList_step1_sbuffSize = 2 * atomBlockDim.x * sizeof(IndexType);
}

void CellList::
rebuild (const MDSystem &	sys,
	 MDTimer *		timer)
{
  
  if (timer != NULL) timer->tic(mdTimeBuildCellList);
  
  IntVectorType tmpNCell;
  bool buildNew = DecideNumCell (sys, mycellSize, mybdir, mydivide, tmpNCell);
  
  if (isempty()){
    if (buildNew){
      printf ("# box size change too much, build cell list\n");
      mallocDeviceCellList (tmpNCell, sys.box.size);
      naivelyBuildDeviceCellList (sys);
    }
    else {
      return ;
    }
  }
  else {
    if (buildNew){
      if (tmpNCell.x != dclist.NCell.x ||
	  tmpNCell.y != dclist.NCell.y ||
	  tmpNCell.z != dclist.NCell.z ){
	printf ("# box size change too much, rebuild cell list\n");
	clearDeviceCellList ();
	mallocDeviceCellList (tmpNCell, sys.box.size);
	naivelyBuildDeviceCellList (sys);
      }
      else{
	buildDeviceCellList (sys);
      }
    }
    else {
      clearDeviceCellList ();
    }
  }
  if (timer != NULL) timer->toc(mdTimeBuildCellList);
}

void CellList::
reshuffle (const IndexType * indexTable,
	   const IndexType & numAtom,
	   MDTimer *timer)
{
  if (! isempty() ){
    if (timer != NULL) timer->tic(mdTimeReshuffleSystem);
    Reshuffle_reshuffleDeviceCellList
	<<<cellGridDim, cellBlockDim>>> (
	    dclist.data,
	    indexTable);
    checkCUDAError ("CellList::reshuffle");
    if (timer != NULL) timer->toc(mdTimeReshuffleSystem);
  }
}


CellList::
CellList (const MDSystem &		sys,
	  const ScalorType &		cellSize,
	  const IndexType &		NThreadCell,
	  const IndexType &		NThreadAtom,
	  const BoxDirection_t &	bdir,
	  const IndexType &		divide)
    : mallocedDeviceCellList(false)
{
  reinit (sys, cellSize, NThreadCell, NThreadAtom, bdir, divide);
}

CellList::
~CellList ()
{
  clearDeviceCellList();
}
