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
  dclist.stride = myBlockDim.x;
  cellGridDim = toGridDim (numCell);

  cudaMalloc ((void**)&(dclist.data), 
	      sizeof(ScalorType) * numCell * dclist.stride);
  cudaMalloc ((void**)&(dclist.numbers), sizeof(unsigned) * numCell);
  cudaMalloc ((void**)&mySendBuff, 
	      sizeof(ScalorType) * numCell * dclist.stride);
  cudaMalloc ((void**)&myTargetBuff, 
	      sizeof(ScalorType) * numCell * dclist.stride);
  checkCUDAError ("NeighborList::init cell list");

  IndexType maxNumNeighborCell = (2*mydivide+1) * (2*mydivide+1) * (2*mydivide+1);
  dclist.maxNumNeighborCell = maxNumNeighborCell;
  cudaMalloc ((void**)&(dclist.numNeighborCell),
	      sizeof(IndexType) * numCell);
  cudaMalloc ((void**)&(dclist.neighborCellIndex),
	      sizeof(IndexType) * maxNumNeighborCell * numCell);
  cudaMalloc ((void**)&(dclist.neighborCellShiftNoi),
	      sizeof(CoordNoiType) * maxNumNeighborCell * numCell);
  checkCUDAError ("NeighborList::maxNumNeighborCell cell list buff");

  mallocedDeviceCellList = true;

  buildDeviceCellList_initBuff<<<cellGridDim, myBlockDim>>> 
      (mySendBuff, myTargetBuff);
  checkCUDAError ("NeighborList::mallocedDeviceCellList cell list buff");
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
    checkCUDAError ("NeighborList::clearDeviceCellList");
  }
}

void CellList::
naivelyBuildDeviceCellList (const MDSystem & sys)
{  
  // buildDeviceCellList_clearBuff
  //     <<<cellGridDim, myBlockDim>>> (
  // 	  mySendBuff);
  // err.check ("NeighborList::rebuild, clear buff");
  prepare_naivelyBuildDeviceCellList
      <<<cellGridDim, myBlockDim>>> (dclist);
  naivelyBuildDeviceCellList2        
      <<<atomGridDim, myBlockDim,
      buildDeviceCellList_step1_sbuffSize>>> (
	  sys.ddata.numAtom,
	  sys.ddata.coord,
	  sys.ddata.coordNoix,
	  sys.ddata.coordNoiy,
	  sys.ddata.coordNoiz,
	  sys.box, dclist,
	  err.ptr_de,
	  err.ptr_dindex, err.ptr_dscalor);
  err.check ("NeighborList::naivelyBuildDeviceCellList");
  checkCUDAError ("NeighborList::naivelyBuildDeviceCellList");

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
      <<<cellGridDim, myBlockDim>>> (
	  mySendBuff);
  err.check ("NeighborList::rebuild, clear buff");
  buildDeviceCellList_step1 
      <<<cellGridDim, myBlockDim,
      buildDeviceCellList_step1_sbuffSize>>> (
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
  err.check ("NeighborList::buildDeviceCellList step 1");
  checkCUDAError ("NeighborList::buildDeviceCellList step 1");
  buildDeviceCellList_step2
      <<<cellGridDim, myBlockDim>>> (
	  sys.box,
	  dclist,
	  mySendBuff,
	  myTargetBuff,
	  bitDeepth,
	  err.ptr_de);
  err.check ("NeighborList::buildDeviceCellList step2");
  checkCUDAError ("NeighborList::buildDeviceCellList step2");
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
  if (CellOnX && dclist.NCell.x < 4) mode = false;
  if (CellOnY && dclist.NCell.y < 4) mode = false;
  if (CellOnZ && dclist.NCell.z < 4) mode = false;

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
	const IndexType &		NTread,
	const BoxDirection_t &		bdir,
	const IndexType &		divide)
{
  myBlockDim.y = 1;
  myBlockDim.z = 1;
  myBlockDim.x = NTread;
  bitDeepth = calDeepth(sys.ddata.numAtom);

  IndexType nob;
  if (sys.ddata.numAtom % myBlockDim.x == 0){
    nob = sys.ddata.numAtom / myBlockDim.x;
  } else {
    nob = sys.ddata.numAtom / myBlockDim.x + 1;
  }
  atomGridDim = toGridDim (nob);

  mybdir = bdir;
  mycellSize = cellSize;
  mydivide = divide;

  buildDeviceCellList_step1_sbuffSize = 2 * myBlockDim.x * sizeof(IndexType);
}

void CellList::
rebuild (const MDSystem &	sys,
	 MDTimer *		timer)
{
  IntVectorType tmpNCell;
  bool buildNew = DecideNumCell (sys, mycellSize, mybdir, mydivide, tmpNCell);

  if (timer != NULL) timer->tic(mdTimeBuildCellList);
  if (isempty()){
    if (buildNew){
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

CellList::
CellList (const MDSystem &		sys,
	  const ScalorType &		cellSize,
	  const IndexType &		NTread,
	  const BoxDirection_t &	bdir,
	  const IndexType &		divide)
    : mallocedDeviceCellList(false)
{
  reinit (sys, cellSize, NTread, bdir, divide);
}
