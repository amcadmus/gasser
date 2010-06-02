#define DEVICE_CODE

#include "Parallel_CellList.h"
#include "Parallel_Interface.h"
#include "Parallel_CellList_device.h"
#include "Parallel_Algorithm.h"
#include "Parallel_BondList.h"
#include "Auxiliary.h"
#include "Parallel_Timer.h"
#include "Parallel_Interface.h"
#include "Parallel_Auxiliary.h"

#include "compile_error_mixcode.h"

void Parallel::DeviceCellListedMDData::
initZeroCell ()
{
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  dim3 gridDim = toGridDim(numCell.x*numCell.y*numCell.z);
  Parallel::CudaGlobal::initZeroCell
      <<<gridDim, numThreadsInCell>>>(
	  numCell, 
	  numAtomInCell);
  checkCUDAError ("DeviceCellListedMDData::initZeroCell");
}

void Parallel::DeviceCellListedMDData::
calNumCell (const ScalorType & input_rlist,
	    const IndexType & input_dividelevel,
	    const BoxDirection_t & bdir,
	    IntVectorType & result_numCell,
	    VectorType & result_frameLow,
	    VectorType & result_frameUp)
{
  int Nx, Ny, Nz;
  Parallel::Interface::numProcDim (Nx, Ny, Nz);
  int ix, iy, iz;
  Parallel::Interface::rankToCartCoord (Parallel::Interface::myRank(), ix, iy, iz);
  double dx, dy, dz;
  dx = getGlobalBoxSize().x / double(Nx);
  dy = getGlobalBoxSize().y / double(Ny);
  dz = getGlobalBoxSize().z / double(Nz);
  result_frameLow.x = dx * ix;
  result_frameLow.y = dy * iy;
  result_frameLow.z = dz * iz;
  result_frameUp.x = result_frameLow.x + dx;
  result_frameUp.y = result_frameLow.y + dy;
  result_frameUp.z = result_frameLow.z + dz;
  
  bool CellOnX, CellOnY, CellOnZ;
  CellOnX = bdir & RectangularBoxGeometry::mdRectBoxDirectionX;
  CellOnY = bdir & RectangularBoxGeometry::mdRectBoxDirectionY;
  CellOnZ = bdir & RectangularBoxGeometry::mdRectBoxDirectionZ;
  double input_rlisti = 1./input_rlist;

  if (CellOnX ) result_numCell.x = int ( floor(dx * input_rlisti) );
  else result_numCell.x = 1;
  if (CellOnY ) result_numCell.y = int ( floor(dy * input_rlisti) );
  else result_numCell.y = 1;
  if (CellOnZ ) result_numCell.z = int ( floor(dz * input_rlisti) );
  else result_numCell.z = 1;

  if ((CellOnX && result_numCell.x < 3) ||
      (CellOnY && result_numCell.y < 3) ||
      (CellOnZ && result_numCell.z < 3) ){
    printf("input_rlist is %f, nx %d ny %d nz %d\n", input_rlist, result_numCell.x, result_numCell.y, result_numCell.z);
    throw MDExcptCellList ("Number of cell on one direction is less than 3");
  }

  // add ghost cell
  VectorType dcell;
  dcell.x = (result_frameUp.x - result_frameLow.x) / result_numCell.x;
  dcell.y = (result_frameUp.y - result_frameLow.y) / result_numCell.y;
  dcell.z = (result_frameUp.z - result_frameLow.z) / result_numCell.z;
  result_frameUp.x += dcell.x;
  result_frameUp.y += dcell.y;
  result_frameUp.z += dcell.z;
  result_frameLow.x -= dcell.x;
  result_frameLow.y -= dcell.y;
  result_frameLow.z -= dcell.z;
  result_numCell.x += 2;
  result_numCell.y += 2;
  result_numCell.z += 2;
  
  if (CellOnX) result_numCell.x *= input_dividelevel;
  if (CellOnY) result_numCell.y *= input_dividelevel;
  if (CellOnZ) result_numCell.z *= input_dividelevel;
}

void Parallel::DeviceCellListedMDData::
initCellStructure (const ScalorType & rlist_,
		   const IndexType & devideLevel_,
		   const BoxDirection_t & bdir)
{
  rlist = rlist_;
  devideLevel = devideLevel_;
  calNumCell (rlist, devideLevel, bdir, numCell, frameLow, frameUp);
  
  DeviceMDData bkData (*this);
  
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell ();
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  // maxNumNeighborCell = 1;
  // if (CellOnX) maxNumNeighborCell *= devideLevel * 2 + 1;
  // if (CellOnY) maxNumNeighborCell *= devideLevel * 2 + 1;
  // if (CellOnZ) maxNumNeighborCell *= devideLevel * 2 + 1;
  
  if (numThreadsInCell * totalNumCell > DeviceMDData::memSize()){
    DeviceMDData::easyMalloc (numThreadsInCell * totalNumCell);
  }
  numData() = totalNumCell * numThreadsInCell;

  // printf ("rank %d, numcell %d\n", Parallel::Interface::myRank(), totalNumCell);
  // getchar ();
  // mallocCell (totalNumCell, maxNumNeighborCell);
  easyMallocCell (totalNumCell);
  initZeroCell ();

  IndexType numThreadBlock = numThreadsInCell;
  dim3 gridDim = toGridDim(numCell.x*numCell.y*numCell.z);

  Parallel::CudaGlobal::formCellStructure
      <<<gridDim, numThreadBlock >>>(
	  frameLow,
	  frameUp,
	  numCell,
	  numAtomInCell,
	  bkData.numData(),
	  bkData.dptr_coordinate(),
	  bkData.dptr_coordinateNoi(),
	  bkData.dptr_velocityX(),
	  bkData.dptr_velocityY(),
	  bkData.dptr_velocityZ(),
	  bkData.dptr_globalIndex(),
	  bkData.dptr_type(),
	  bkData.dptr_mass(),
	  bkData.dptr_charge(),
	  coord,
	  coordNoi,
	  velox,
	  veloy,
	  veloz,
	  globalIndex,
	  type,
	  mass,
	  charge,
	  err.ptr_de);
  checkCUDAError ("DeviceCellListedMDData::formCellStructure");
  err.updateHost();
  err.check ("DeviceCellListedMDData::formCellSturcture");
}

bool Parallel::DeviceCellListedMDData::
reinitCellStructure (const ScalorType & rlist_,
		     const IndexType & devideLevel_,
		     const BoxDirection_t & bdir)
{
  IntVectorType tmpNumCell;
  VectorType tmpFrameLow, tmpFrameUp;
  calNumCell (rlist_, devideLevel_, bdir, tmpNumCell, tmpFrameLow, tmpFrameUp);
  if (tmpNumCell.x == numCell.x &&
      tmpNumCell.y == numCell.y &&
      tmpNumCell.z == numCell.z ){
    rlist = rlist_;
    devideLevel = devideLevel_;
    return false;
  }
  
  DeviceCellListedMDData bkData ;
  bkData.copyFromDevice (*this, MDDataItemMask_All);
  rlist = rlist_;
  devideLevel = devideLevel_;
  numCell = tmpNumCell;
  frameLow = tmpFrameLow;
  frameUp = tmpFrameUp;

  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell ();
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  
  if (numThreadsInCell * totalNumCell > DeviceMDData::memSize()){
    DeviceMDData::easyMalloc (numThreadsInCell * totalNumCell,
			      getMaxNumBond(), getMaxNumAngle(), getMaxNumDihedral());
  }
  numData() = totalNumCell * numThreadsInCell;

  easyMallocCell (totalNumCell);
  initZeroCell ();

  IndexType numThreadBlock = numThreadsInCell;
  IntVectorType numCell_ = bkData.getNumCell();
  dim3 gridDim = toGridDim(numCell_.x*numCell_.y*numCell_.z);
  IndexType * forwardMap;
  cudaMalloc ((void**)&forwardMap,
	      sizeof(IndexType) * numCell_.x * numCell_.y * numCell_.z * numThreadsInCell);
  checkCUDAError ("DeviceCellListedMDData::reinitCellStructure malloc forwardMap");

  Parallel::CudaGlobal::reinitCellStructure_calForwardMap
      <<<gridDim, numThreadBlock >>>(
	  bkData.dptr_numAtomInCell(),
	  bkData.dptr_coordinate(),
	  frameLow,
	  frameUp,
	  numCell,
	  numAtomInCell,
	  forwardMap,
	  err.ptr_de);
  checkCUDAError ("DeviceCellListedMDData::reinitCellStructure_calForwardMap");
  err.updateHost();
  err.check ("DeviceCellListedMDData::reinitCellStructure_calForwardMap");

  Parallel::CudaGlobal::reinitCellStructure_step1
      <<<gridDim, numThreadBlock >>>(
	  bkData.dptr_numAtomInCell(),
	  bkData.dptr_coordinate(),
	  bkData.dptr_coordinateNoi(),
	  bkData.dptr_velocityX(),
	  bkData.dptr_velocityY(),
	  bkData.dptr_velocityZ(),
	  bkData.dptr_forceX(),
	  bkData.dptr_forceY(),
	  bkData.dptr_forceZ(),
	  bkData.dptr_globalIndex(),
	  bkData.dptr_type(),
	  bkData.dptr_mass(),
	  bkData.dptr_charge(),
	  forwardMap,
	  coord,
	  coordNoi,
	  velox,
	  veloy,
	  veloz,
	  forcx,
	  forcy,
	  forcz,
	  globalIndex,
	  type,
	  mass,
	  charge);
  checkCUDAError ("DeviceCellListedMDData::reinitCellStructure_step1");
  Parallel::CudaGlobal::reinitCellStructure_step2
      <<<gridDim, numThreadBlock>>>(
	  bkData.dptr_numAtomInCell(),
	  bkData.getMaxNumBond(),
	  bkData.dptr_numBond(),
	  bkData.dptr_bondIndex(),
	  bkData.dptr_bondNeighbor_globalIndex(),
	  bkData.getMaxNumAngle(),
	  bkData.dptr_numAngle(),
	  bkData.dptr_angleIndex(),
	  bkData.dptr_anglePosi(),
	  bkData.dptr_angleNeighbor_globalIndex(),
	  bkData.getMaxNumDihedral(),
	  bkData.dptr_numDihedral(),
	  bkData.dptr_dihedralIndex(),
	  bkData.dptr_dihedralPosi(),
	  bkData.dptr_dihedralNeighbor_globalIndex(),
	  bkData.getBondTopStride(),
	  forwardMap,
	  dptr_numBond(),
	  dptr_bondIndex(),
	  dptr_bondNeighbor_globalIndex(),
	  dptr_numAngle(),
	  dptr_angleIndex(),
	  dptr_anglePosi(),
	  dptr_angleNeighbor_globalIndex(),
	  dptr_numDihedral(),
	  dptr_dihedralIndex(),
	  dptr_dihedralPosi(),
	  dptr_dihedralNeighbor_globalIndex(),
	  getBondTopStride());      
  cudaFree (forwardMap);

  return true;
}
  

void Parallel::DeviceCellListedMDData::
rebuild ()
{
  IndexType numThreadBlock = Parallel::Interface::numThreadsInCell();
  IndexType totalNumCell = numCell.x*numCell.y*numCell.z;
  dim3 gridDim = toGridDim(totalNumCell);

  IndexType * bk_numAtomInCell;
  cudaMalloc ((void**)&bk_numAtomInCell, totalNumCell * sizeof(IndexType));
  cudaMemcpy (bk_numAtomInCell, numAtomInCell, totalNumCell * sizeof(IndexType),
	      cudaMemcpyDeviceToDevice);
  checkCUDAError ("DeviceCellListedMDData::rebuild malloc backup");

  if (getMaxNumBond() == 0 &&
      getMaxNumAngle() == 0 &&
      getMaxNumDihedral() == 0){
    Parallel::CudaGlobal::rebuildCellList_step1
	<<<gridDim, numThreadBlock>>> (
	    frameLow,
	    frameUp,
	    numCell,
	    bk_numAtomInCell,
	    numAtomInCell,
	    coord,
	    coordNoi,
	    velox,
	    veloy,
	    veloz,
	    forcx,
	    forcy,
	    forcz,
	    globalIndex,
	    type,
	    mass,
	    charge,
	    err.ptr_de);
    checkCUDAError ("DeviceCellListedMDData::rebuild step1");
    err.updateHost();
    err.check ("DeviceCellListedMDData::rebuild step1");
    Parallel::CudaGlobal::rebuildCellList_step2
	<<<gridDim, numThreadBlock, numThreadBlock*sizeof(IndexType)*3>>> (
	    numAtomInCell,
	    coord,
	    coordNoi,
	    velox,
	    veloy,
	    veloz,
	    forcx,
	    forcy,
	    forcz,
	    globalIndex,
	    type,
	    mass,
	    charge,
	    err.ptr_de);
    checkCUDAError ("DeviceCellListedMDData::rebuild step2");
    err.updateHost();
    err.check ("DeviceCellListedMDData::rebuild step2");  
  }
  else{
    Parallel::CudaGlobal::rebuildCellList_step1
	<<< gridDim, numThreadBlock>>> (
	    frameLow,
	    frameUp,
	    numCell,
	    bk_numAtomInCell,
	    numAtomInCell,
	    coord,
	    coordNoi,
	    velox,
	    veloy,
	    veloz,
	    forcx,
	    forcy,
	    forcz,
	    globalIndex,
	    type,
	    mass,
	    charge,
	    forwardMap_step1,
	    err.ptr_de,
	    err.ptr_dindex,
	    err.ptr_dscalor);
    checkCUDAError ("DeviceCellListedMDData::rebuild step1");
    err.updateHost();
    err.check ("DeviceCellListedMDData::rebuild step1");
    Parallel::CudaGlobal::rebuildCellList_step1_mapBondTop
	<<<gridDim, numThreadBlock>>> (
	    forwardMap_step1,
	    bk_numAtomInCell,
	    getMaxNumBond(),
	    dptr_numBond(),
	    dptr_bondIndex(),
	    dptr_bondNeighbor_globalIndex(),
	    getMaxNumAngle(),
	    dptr_numAngle(),
	    dptr_angleIndex(),
	    dptr_anglePosi(),
	    dptr_angleNeighbor_globalIndex(),
	    getMaxNumDihedral(),
	    dptr_numDihedral(),
	    dptr_dihedralIndex(),
	    dptr_dihedralPosi(),
	    dptr_angleNeighbor_globalIndex(),
	    bondTopStride());
    checkCUDAError ("DeviceCellListedMDData::rebuild map top step1");

    cudaMemcpy (bk_numAtomInCell, dptr_numAtomInCell(),
		totalNumCell * sizeof(IndexType), cudaMemcpyDeviceToDevice);
    
    Parallel::CudaGlobal::rebuildCellList_step2
	<<<gridDim, numThreadBlock, numThreadBlock*sizeof(IndexType)*3>>> (
	    numAtomInCell,
	    coord,
	    coordNoi,
	    velox,
	    veloy,
	    veloz,
	    forcx,
	    forcy,
	    forcz,
	    globalIndex,
	    type,
	    mass,
	    charge,
	    forwardMap_step2,
	    err.ptr_de);
    checkCUDAError ("DeviceCellListedMDData::rebuild step2");
    err.updateHost();
    err.check ("DeviceCellListedMDData::rebuild step2");
    Parallel::CudaGlobal::rebuildCellList_step2_mapBondTop
	<<<gridDim, numThreadBlock>>> (
	    forwardMap_step2,
	    bk_numAtomInCell,
	    getMaxNumBond(),
	    dptr_numBond(),
	    dptr_bondIndex(),
	    dptr_bondNeighbor_globalIndex(),
	    getMaxNumAngle(),
	    dptr_numAngle(),
	    dptr_angleIndex(),
	    dptr_anglePosi(),
	    dptr_angleNeighbor_globalIndex(),
	    getMaxNumDihedral(),
	    dptr_numDihedral(),
	    dptr_dihedralIndex(),
	    dptr_dihedralPosi(),
	    dptr_angleNeighbor_globalIndex(),
	    bondTopStride());
    checkCUDAError ("DeviceCellListedMDData::rebuild map top step2");
  }

  cudaFree (bk_numAtomInCell);
  checkCUDAError ("DeviceCellListedMDData::rebuild free backup");  
}


void Parallel::DeviceCellListedMDData::
rebuild (DeviceBondList & dbdlist)
{
  IndexType numThreadBlock = Parallel::Interface::numThreadsInCell();
  IndexType totalNumCell = numCell.x*numCell.y*numCell.z;
  dim3 gridDim = toGridDim(totalNumCell);

  IndexType * bk_numAtomInCell;
  cudaMalloc ((void**)&bk_numAtomInCell, totalNumCell * sizeof(IndexType));
  cudaMemcpy (bk_numAtomInCell, numAtomInCell, totalNumCell * sizeof(IndexType),
	      cudaMemcpyDeviceToDevice);
  checkCUDAError ("DeviceCellListedMDData::rebuild malloc backup");

  if (getMaxNumBond() == 0 &&
      getMaxNumAngle() == 0 &&
      getMaxNumDihedral() == 0){
    Parallel::CudaGlobal::rebuildCellList_step1
	<<<gridDim, numThreadBlock>>> (
	    frameLow,
	    frameUp,
	    numCell,
	    bk_numAtomInCell,
	    numAtomInCell,
	    coord,
	    coordNoi,
	    velox,
	    veloy,
	    veloz,
	    forcx,
	    forcy,
	    forcz,
	    globalIndex,
	    type,
	    mass,
	    charge,
	    err.ptr_de);
    checkCUDAError ("DeviceCellListedMDData::rebuild step1");
    err.updateHost();
    err.check ("DeviceCellListedMDData::rebuild step1");
    Parallel::CudaGlobal::rebuildCellList_step2
	<<<gridDim, numThreadBlock, numThreadBlock*sizeof(IndexType)*3>>> (
	    numAtomInCell,
	    coord,
	    coordNoi,
	    velox,
	    veloy,
	    veloz,
	    forcx,
	    forcy,
	    forcz,
	    globalIndex,
	    type,
	    mass,
	    charge,
	    err.ptr_de);
    checkCUDAError ("DeviceCellListedMDData::rebuild step2");
    err.updateHost();
    err.check ("DeviceCellListedMDData::rebuild step2");
  }
  else {
    Parallel::CudaGlobal::rebuildCellList_step1
	<<< gridDim, numThreadBlock>>> (
	    frameLow,
	    frameUp,
	    numCell,
	    bk_numAtomInCell,
	    numAtomInCell,
	    coord,
	    coordNoi,
	    velox,
	    veloy,
	    veloz,
	    forcx,
	    forcy,
	    forcz,
	    globalIndex,
	    type,
	    mass,
	    charge,
	    forwardMap_step1,
	    err.ptr_de,
	    err.ptr_dindex,
	    err.ptr_dscalor);
    checkCUDAError ("DeviceCellListedMDData::rebuild step1");
    err.updateHost();
    err.check ("DeviceCellListedMDData::rebuild step1");
    Parallel::CudaGlobal::rebuildCellList_step1_mapBondTop
	<<<gridDim, numThreadBlock>>> (
	    forwardMap_step1,
	    bk_numAtomInCell,
	    getMaxNumBond(),
	    dptr_numBond(),
	    dptr_bondIndex(),
	    dptr_bondNeighbor_globalIndex(),
	    dbdlist.dptr_bondNeighbor_localIndex(),
	    getMaxNumAngle(),
	    dptr_numAngle(),
	    dptr_angleIndex(),
	    dptr_anglePosi(),
	    dptr_angleNeighbor_globalIndex(),
	    dbdlist.dptr_angleNeighbor_localIndex(),
	    getMaxNumDihedral(),
	    dptr_numDihedral(),
	    dptr_dihedralIndex(),
	    dptr_dihedralPosi(),
	    dptr_angleNeighbor_globalIndex(),
	    dbdlist.dptr_dihedralNeighbor_localIndex(),
	    bondTopStride());
    checkCUDAError ("DeviceCellListedMDData::rebuild map top step1");

    cudaMemcpy (bk_numAtomInCell, dptr_numAtomInCell(),
		totalNumCell * sizeof(IndexType), cudaMemcpyDeviceToDevice);
    
    Parallel::CudaGlobal::rebuildCellList_step2
	<<<gridDim, numThreadBlock, numThreadBlock*sizeof(IndexType)*3>>> (
	    numAtomInCell,
	    coord,
	    coordNoi,
	    velox,
	    veloy,
	    veloz,
	    forcx,
	    forcy,
	    forcz,
	    globalIndex,
	    type,
	    mass,
	    charge,
	    forwardMap_step2,
	    err.ptr_de);
    checkCUDAError ("DeviceCellListedMDData::rebuild step2");
    err.updateHost();
    err.check ("DeviceCellListedMDData::rebuild step2");
    Parallel::CudaGlobal::rebuildCellList_step2_mapBondTop
	<<<gridDim, numThreadBlock>>> (
	    forwardMap_step2,
	    bk_numAtomInCell,
	    getMaxNumBond(),
	    dptr_numBond(),
	    dptr_bondIndex(),
	    dptr_bondNeighbor_globalIndex(),
	    dbdlist.dptr_bondNeighbor_localIndex(),
	    getMaxNumAngle(),
	    dptr_numAngle(),
	    dptr_angleIndex(),
	    dptr_anglePosi(),
	    dptr_angleNeighbor_globalIndex(),
	    dbdlist.dptr_angleNeighbor_localIndex(),
	    getMaxNumDihedral(),
	    dptr_numDihedral(),
	    dptr_dihedralIndex(),
	    dptr_dihedralPosi(),
	    dptr_angleNeighbor_globalIndex(),
	    dbdlist.dptr_dihedralNeighbor_localIndex(),
	    bondTopStride());
    checkCUDAError ("DeviceCellListedMDData::rebuild map top step2");
  }

  cudaFree (bk_numAtomInCell);
  checkCUDAError ("DeviceCellListedMDData::rebuild free backup");  
}

  

__global__ void Parallel::CudaGlobal::
initZeroCell (const IntVectorType numCell,
	      IndexType * numAtomInCell)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;

  if (ii < totalNumCell){
    numAtomInCell[ii] = 0;
    // numNeighborCell[ii] = 0;
  }
}


__global__ void Parallel::CudaGlobal::
formCellStructure (const VectorType frameLow,
		   const VectorType frameUp,
		   const IntVectorType numCell,
		   IndexType * numAtomInCell,
		   const IndexType numAtom,
		   const CoordType  * bk_coord,
		   const CoordNoiType * bk_coordNoi,
		   const ScalorType * bk_velox,
		   const ScalorType * bk_veloy,
		   const ScalorType * bk_veloz,
		   const IndexType  * bk_globalIndex,
		   const TypeType   * bk_type,
		   const ScalorType * bk_mass,
		   const ScalorType * bk_charge,
		   CoordType  * coord,
		   CoordNoiType * coordNoi,
		   ScalorType * velox,
		   ScalorType * veloy,
		   ScalorType * veloz,
		   IndexType  * globalIndex,
		   TypeType   * type,
		   ScalorType * mass,
		   ScalorType * charge,
		   mdError_t * ptr_de)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  __shared__ ScalorType dcellxi;
  __shared__ ScalorType dcellyi;
  __shared__ ScalorType dcellzi;
  if (tid == 0) {
    dcellxi = ScalorType(numCell.x) / (frameUp.x - frameLow.x);
  }
  if (tid == 1){
    dcellyi = ScalorType(numCell.y) / (frameUp.y - frameLow.y);
  }
  if (tid == 2){
    dcellzi = ScalorType(numCell.z) / (frameUp.z - frameLow.z);
  }
  __syncthreads();
  
  IndexType targetIndex;
  if (ii < numAtom){
    IndexType targetCellx, targetCelly, targetCellz;
    targetCellx = IndexType((bk_coord[ii].x - frameLow.x) * dcellxi);
    targetCelly = IndexType((bk_coord[ii].y - frameLow.y) * dcellyi);
    targetCellz = IndexType((bk_coord[ii].z - frameLow.z) * dcellzi);
    // if (targetCellx == numCell.x){
    //   targetCellx = numCell.x - 1;
    // }
    // if (targetCelly == numCell.y){
    //   targetCelly = numCell.y - 1;
    // }
    // if (targetCellz == numCell.z){
    //   targetCellz = numCell.z - 1;
    // }
    if (ptr_de != NULL && 
	(targetCellx >= numCell.x || 
	 targetCelly >= numCell.y || 
	 targetCellz >= numCell.z)){
      *ptr_de = mdErrorOverFlowCellIdx;
      return;
    }

    IndexType cellid = CudaDevice::D3toD1
	(numCell, targetCellx, targetCelly, targetCellz);

    IndexType pid = atomicAdd (&numAtomInCell[cellid], 1);
    if (pid >= blockDim.x){
      *ptr_de = mdErrorShortCellList;
      pid = 0;
    }
    targetIndex = pid + cellid * blockDim.x;
    coord[targetIndex] = bk_coord[ii];
    coordNoi[targetIndex].x = bk_coordNoi[ii].x;
    coordNoi[targetIndex].y = bk_coordNoi[ii].y;
    coordNoi[targetIndex].z = bk_coordNoi[ii].z;
    velox[targetIndex] = bk_velox[ii];
    veloy[targetIndex] = bk_veloy[ii];
    veloz[targetIndex] = bk_veloz[ii];
    globalIndex[targetIndex] = bk_globalIndex[ii];
    type[targetIndex] = bk_type[ii];
    mass[targetIndex] = bk_mass[ii];
    charge[targetIndex] = bk_charge[ii];
  }
}


__global__ void Parallel::CudaGlobal::
reinitCellStructure_calForwardMap (const IndexType * bk_numAtomInCell,
				   const CoordType * bk_coord,
				   const VectorType frameLow,
				   const VectorType frameUp,
				   const IntVectorType numCell,
				   IndexType * numAtomInCell,
				   IndexType * forwardMap,
				   mdError_t * ptr_de)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  IndexType this_numAtom = bk_numAtomInCell[bid];
  if (this_numAtom == 0) return;

  ScalorType dcellxi;
  ScalorType dcellyi;
  ScalorType dcellzi;
  dcellxi = ScalorType(numCell.x) / (frameUp.x - frameLow.x);
  dcellyi = ScalorType(numCell.y) / (frameUp.y - frameLow.y);
  dcellzi = ScalorType(numCell.z) / (frameUp.z - frameLow.z);

  if (tid < this_numAtom){
    IndexType targetCellx, targetCelly, targetCellz;
    targetCellx = IndexType((bk_coord[ii].x - frameLow.x) * dcellxi);
    targetCelly = IndexType((bk_coord[ii].y - frameLow.y) * dcellyi);
    targetCellz = IndexType((bk_coord[ii].z - frameLow.z) * dcellzi);
    if (ptr_de != NULL && 
	(targetCellx >= numCell.x || 
	 targetCelly >= numCell.y || 
	 targetCellz >= numCell.z)){
      *ptr_de = mdErrorOverFlowCellIdx;
      return;
    }
    IndexType cellid = CudaDevice::D3toD1
	(numCell, targetCellx, targetCelly, targetCellz);
    IndexType pid = atomicAdd (&numAtomInCell[cellid], 1);
    if (pid >= blockDim.x){
      *ptr_de = mdErrorShortCellList;
      pid = 0;
    }
    IndexType targetIndex = pid + cellid * blockDim.x;
    forwardMap[ii] = targetIndex;
  }
}


__global__ void Parallel::CudaGlobal::
reinitCellStructure_step1 (const IndexType * bk_numAtomInCell,
			   const CoordType  * bk_coord,
			   const CoordNoiType * bk_coordNoi,
			   const ScalorType * bk_velox,
			   const ScalorType * bk_veloy,
			   const ScalorType * bk_veloz,
			   const ScalorType * bk_forcx,
			   const ScalorType * bk_forcy,
			   const ScalorType * bk_forcz,
			   const IndexType  * bk_globalIndex,
			   const TypeType   * bk_type,
			   const ScalorType * bk_mass,
			   const ScalorType * bk_charge,
			   const IndexType * forwardMap,
			   CoordType  * coord,
			   CoordNoiType * coordNoi,
			   ScalorType * velox,
			   ScalorType * veloy,
			   ScalorType * veloz,
			   ScalorType * forcx,
			   ScalorType * forcy,
			   ScalorType * forcz,
			   IndexType  * globalIndex,
			   TypeType   * type,
			   ScalorType * mass,
			   ScalorType * charge)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  IndexType this_numAtom = bk_numAtomInCell[bid];
  if (this_numAtom == 0) return;

  if (tid < this_numAtom){
    IndexType targetIndex = forwardMap[ii];
    coord[targetIndex] = bk_coord[ii];
    coordNoi[targetIndex].x = bk_coordNoi[ii].x;
    coordNoi[targetIndex].y = bk_coordNoi[ii].y;
    coordNoi[targetIndex].z = bk_coordNoi[ii].z;
    velox[targetIndex] = bk_velox[ii];
    veloy[targetIndex] = bk_veloy[ii];
    veloz[targetIndex] = bk_veloz[ii];
    forcx[targetIndex] = bk_forcx[ii];
    forcy[targetIndex] = bk_forcy[ii];
    forcz[targetIndex] = bk_forcz[ii];
    globalIndex[targetIndex] = bk_globalIndex[ii];
    type[targetIndex] = bk_type[ii];
    mass[targetIndex] = bk_mass[ii];
    charge[targetIndex] = bk_charge[ii];
  }
}

__global__ void Parallel::CudaGlobal::
reinitCellStructure_step2 (const IndexType * bk_numAtomInCell,
			   const IndexType bk_maxNumBond,
			   const IndexType * bk_numBond,
			   const IndexType * bk_bondIndex,
			   const IndexType * bk_bondNeighbor_globalIndex,
			   const IndexType bk_maxNumAngle,
			   const IndexType * bk_numAngle,
			   const IndexType * bk_angleIndex,
			   const IndexType * bk_anglePosi,			   
			   const IndexType * bk_angleNeighbor_globalIndex,
			   const IndexType bk_maxNumDihedral,
			   const IndexType * bk_numDihedral,
			   const IndexType * bk_dihedralIndex,
			   const IndexType * bk_dihedralPosi,
			   const IndexType * bk_dihedralNeighbor_globalIndex,
			   const IndexType bk_bondTopStride,
			   const IndexType * forwardMap,
			   IndexType * numBond,
			   IndexType * bondIndex,
			   IndexType * bondNeighbor_globalIndex,
			   IndexType * numAngle,
			   IndexType * angleIndex,
			   IndexType * anglePosi,
			   IndexType * angleNeighbor_globalIndex,
			   IndexType * numDihedral,
			   IndexType * dihedralIndex,
			   IndexType * dihedralPosi,
			   IndexType * dihedralNeighbor_globalIndex,
			   const IndexType bondTopStride)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  IndexType this_numAtom = bk_numAtomInCell[bid];
  if (this_numAtom == 0) return;

  if (tid < this_numAtom){
    IndexType my_targetIndex = forwardMap[ii];
    if (bk_maxNumBond != 0){
      IndexType tmpNum = numBond[my_targetIndex] = bk_numBond[ii];
      IndexType tmp_target = my_targetIndex;
      IndexType tmp_src = ii;
      for (IndexType jj = 0; jj < tmpNum; ++jj){
	bondIndex[tmp_target] = bk_bondIndex[tmp_src];
	bondNeighbor_globalIndex[tmp_target] = bk_bondNeighbor_globalIndex[tmp_src];
	tmp_target += bondTopStride;
	tmp_src    += bk_bondTopStride;
      }
    }
    if (bk_maxNumAngle != 0){
      IndexType tmpNum = numAngle[my_targetIndex] = bk_numAngle[ii];
      IndexType tmp_target = my_targetIndex;
      IndexType tmp_src = ii;
      for (IndexType jj = 0; jj < tmpNum; ++jj){
	angleIndex[tmp_target] = bk_angleIndex[tmp_src];
	anglePosi [tmp_target] = bk_anglePosi [tmp_src];
	tmp_target += bondTopStride;
	tmp_src    += bk_bondTopStride;
      }
      tmp_target = my_targetIndex;
      tmp_src = ii;
      for (IndexType jj = 0; jj < tmpNum * 2; ++jj){
	bondNeighbor_globalIndex[tmp_target] = bk_bondNeighbor_globalIndex[tmp_src];
	tmp_target += bondTopStride;
	tmp_src    += bk_bondTopStride;
      }  
    }
    if (bk_maxNumDihedral != 0){
      IndexType tmpNum = numDihedral[my_targetIndex] = bk_numDihedral[ii];
      IndexType tmp_target = my_targetIndex;
      IndexType tmp_src = ii;
      for (IndexType jj = 0; jj < tmpNum; ++jj){
	dihedralIndex[tmp_target] = bk_dihedralIndex[tmp_src];
	dihedralPosi [tmp_target] = bk_dihedralPosi [tmp_src];
	tmp_target += bondTopStride;
	tmp_src    += bk_bondTopStride;
      }
      tmp_target = my_targetIndex;
      tmp_src = ii;
      for (IndexType jj = 0; jj < tmpNum * 2; ++jj){
	bondNeighbor_globalIndex[tmp_target] = bk_bondNeighbor_globalIndex[tmp_src];
	tmp_target += bondTopStride;
	tmp_src    += bk_bondTopStride;
      }
    }
  }
}


__global__ void Parallel::CudaGlobal::
rebuildCellList_step1 (const VectorType frameLow,
		       const VectorType frameUp,
		       const IntVectorType numCell,
		       const IndexType * bk_numAtomInCell,
		       IndexType * numAtomInCell,
		       CoordType * coord,
		       CoordNoiType * coordNoi,
		       ScalorType * velox,
		       ScalorType * veloy,
		       ScalorType * veloz,
		       ScalorType * forcx,
		       ScalorType * forcy,
		       ScalorType * forcz,
		       IndexType  * globalIndex,
		       TypeType   * type,
		       ScalorType * mass,
		       ScalorType * charge,
		       mdError_t * ptr_de)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  __shared__ ScalorType dcellxi;
  __shared__ ScalorType dcellyi;
  __shared__ ScalorType dcellzi;
  if (tid == 0) {
    dcellxi = ScalorType(numCell.x) / (frameUp.x - frameLow.x);
  }
  if (tid == 1){
    dcellyi = ScalorType(numCell.y) / (frameUp.y - frameLow.y);
  }
  if (tid == 2){
    dcellzi = ScalorType(numCell.z) / (frameUp.z - frameLow.z);
  }
  __syncthreads();

  // IndexType mark = MaxIndexValue - (MaxIndexValue >> 1);
  // IndexType mystat ;
  
  if (tid < bk_numAtomInCell[bid]){
    // mystat = globalIndex[ii];
    IndexType targetCellx, targetCelly, targetCellz;
    targetCellx = IndexType((coord[ii].x - frameLow.x) * dcellxi);
    targetCelly = IndexType((coord[ii].y - frameLow.y) * dcellyi);
    targetCellz = IndexType((coord[ii].z - frameLow.z) * dcellzi);
    // printf ("%d %d %d %d %f %f %f\n", ii, targetCellx, targetCelly, targetCellz,
    // 	    coord[ii].x, coord[ii].y, coord[ii].z);
    // if (targetCellx == numCell.x){
    //   targetCellx = numCell.x - 1;
    // }
    // if (targetCelly == numCell.y){
    //   targetCelly = numCell.y - 1;
    // }
    // if (targetCellz == numCell.z){
    //   targetCellz = numCell.z - 1;
    // }
    if (ptr_de != NULL && 
	(targetCellx >= numCell.x || 
	 targetCelly >= numCell.y || 
	 targetCellz >= numCell.z)){
      *ptr_de = mdErrorOverFlowCellIdx;
      return;
    }
    IndexType cellid = CudaDevice::D3toD1
	(numCell, targetCellx, targetCelly, targetCellz);
    if (cellid != bid){
      // IndexType pid = atomicAdd (&numAtomInCell[cellid], 1);
      // if (pid >= blockDim.x){
      // 	*ptr_de = mdErrorShortCellList;
      // 	pid = 0;
      // }
      IndexType pid = atomicAdd (&numAtomInCell[cellid], 1);
      if (pid >= blockDim.x){
	*ptr_de = mdErrorShortCellList;
	pid = 0;
      }
      IndexType targetIndex = pid + cellid * blockDim.x;
      coord[targetIndex] = coord[ii];
      coordNoi[targetIndex].x = coordNoi[ii].x;
      coordNoi[targetIndex].y = coordNoi[ii].y;
      coordNoi[targetIndex].z = coordNoi[ii].z;
      velox[targetIndex] = velox[ii];
      veloy[targetIndex] = veloy[ii];
      veloz[targetIndex] = veloz[ii];
      forcx[targetIndex] = forcx[ii];
      forcy[targetIndex] = forcy[ii];
      forcz[targetIndex] = forcz[ii];
      globalIndex[targetIndex] = globalIndex[ii];
      globalIndex[ii] = MaxIndexValue;
      type[targetIndex] = type[ii];
      mass[targetIndex] = mass[ii];
      charge[targetIndex] = charge[ii];
    }
  }
  // globalIndex[ii] = mystat;
  return;
}


static __device__ IndexType
headSort (volatile IndexType * index,
	  volatile IndexType * sbuff)
{
  IndexType k = NUintBit - 1;
  IndexType tid = threadIdx.x;
  sbuff[tid] = getKthBit(index[tid], k);
  sbuff[tid+blockDim.x] = 0;
  
  __syncthreads();
  IndexType total1 = sumVectorBlockBuffer (sbuff, blockDim.x);
  IndexType target, mydata = index[tid];
  __syncthreads();
  
  if (getKthBit(index[tid], k)) {
    target = blockDim.x - sbuff[tid];
  }
  else {
    target = tid + sbuff[tid] - total1;
  }
  __syncthreads();
  index[target] = mydata;
  __syncthreads();

  return total1;
}


__global__ void Parallel::CudaGlobal::
rebuildCellList_step2 (IndexType * numAtomInCell,
		       CoordType  * coord,
		       CoordNoiType * coordNoi,
		       ScalorType * velox,
		       ScalorType * veloy,
		       ScalorType * veloz,
		       ScalorType * forcx,
		       ScalorType * forcy,
		       ScalorType * forcz,
		       IndexType  * globalIndex,
		       TypeType   * type,
		       ScalorType * mass,
		       ScalorType * charge,
		       mdError_t * ptr_de)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  // IndexType k = NUintBit - 1;
  
  extern __shared__ volatile IndexType sbuff[];
  volatile IndexType * myIndex = (volatile IndexType * )sbuff;

  if (tid >= numAtomInCell[bid] || globalIndex[ii] == MaxIndexValue){
    myIndex[tid] = MaxIndexValue;
  }
  else {
    myIndex[tid] = tid;
  }
  __syncthreads();
  IndexType total = headSort (myIndex, &sbuff[blockDim.x]);
  total = blockDim.x - total;

  IndexType fromId;
  CoordType bk_coord;
  CoordNoiType bk_coordNoi;
  ScalorType bk_velox, bk_veloy, bk_veloz;
  ScalorType bk_forcx, bk_forcy, bk_forcz;
  IndexType bk_globalIndex;
  TypeType bk_type;
  ScalorType bk_mass;
  ScalorType bk_charge;
  
  if (tid < total){
    fromId = myIndex[tid] + bid * blockDim.x;
    if (ii != fromId){
      bk_coord = coord[fromId];
      bk_coordNoi.x = coordNoi[fromId].x;
      bk_coordNoi.y = coordNoi[fromId].y;
      bk_coordNoi.z = coordNoi[fromId].z;
      bk_velox = velox[fromId];
      bk_veloy = veloy[fromId];
      bk_veloz = veloz[fromId];
      bk_forcx = forcx[fromId];
      bk_forcy = forcy[fromId];
      bk_forcz = forcz[fromId];
      bk_globalIndex = globalIndex[fromId];
      bk_type = type[fromId];
      bk_mass = mass[fromId];
      bk_charge = charge[fromId];
    }
  }
  __syncthreads();

  if (tid < total && ii != fromId){
    coord[ii] = bk_coord;
    coordNoi[ii].x = bk_coordNoi.x;
    coordNoi[ii].y = bk_coordNoi.y;
    coordNoi[ii].z = bk_coordNoi.z;
    velox[ii] = bk_velox;
    veloy[ii] = bk_veloy;
    veloz[ii] = bk_veloz;
    forcx[ii] = bk_forcx;
    forcy[ii] = bk_forcy;
    forcz[ii] = bk_forcz;
    globalIndex[ii] = bk_globalIndex;
    type[ii] = bk_type;
    mass[ii] = bk_mass;
    charge[ii] = bk_charge;
  }

    
  // else {
  //   globalIndex[ii] = MaxIndexValue;
  // }  

  if (tid == 0){
    numAtomInCell[bid] = total;
  }
}

__global__ void Parallel::CudaGlobal::
rebuildCellList_step1 (const VectorType frameLow,
		       const VectorType frameUp,
		       const IntVectorType numCell,
		       const IndexType * bk_numAtomInCell,
		       IndexType * numAtomInCell,
		       CoordType * coord,
		       CoordNoiType * coordNoi,
		       ScalorType * velox,
		       ScalorType * veloy,
		       ScalorType * veloz,
		       ScalorType * forcx,
		       ScalorType * forcy,
		       ScalorType * forcz,
		       IndexType  * globalIndex,
		       TypeType   * type,
		       ScalorType * mass,
		       ScalorType * charge,
		       IndexType * forwardMap,
		       mdError_t * ptr_de,
		       IndexType * erridx,
		       ScalorType * errsrc)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  __shared__ ScalorType dcellxi;
  __shared__ ScalorType dcellyi;
  __shared__ ScalorType dcellzi;
  if (tid == 0) {
    dcellxi = ScalorType(numCell.x) / (frameUp.x - frameLow.x);
    dcellyi = ScalorType(numCell.y) / (frameUp.y - frameLow.y);
    dcellzi = ScalorType(numCell.z) / (frameUp.z - frameLow.z);
  }
  forwardMap[ii] = MaxIndexValue;
  __syncthreads();
  
  bool inRange = tid < bk_numAtomInCell[bid];
  if (inRange){
    IndexType targetCellx, targetCelly, targetCellz;
    targetCellx = IndexType((coord[ii].x - frameLow.x) * dcellxi);
    targetCelly = IndexType((coord[ii].y - frameLow.y) * dcellyi);
    targetCellz = IndexType((coord[ii].z - frameLow.z) * dcellzi);
    if (ptr_de != NULL && 
	(targetCellx >= numCell.x || 
	 targetCelly >= numCell.y || 
	 targetCellz >= numCell.z)){
      *ptr_de = mdErrorOverFlowCellIdx;
      if (targetCellx >= numCell.x){
	*erridx = targetCellx;
	*errsrc = coord[ii].x;
	errsrc[1] = 0;
	errsrc[2] = frameLow.x;
	errsrc[3] = dcellxi;
	errsrc[4] = numCell.x;
      }
      if (targetCelly >= numCell.y) {
	*erridx = targetCelly;
	*errsrc = coord[ii].y;
	errsrc[1] = 1;
	errsrc[2] = frameLow.y;
	errsrc[3] = dcellyi;
	errsrc[4] = numCell.y;
      }
      if (targetCellz >= numCell.z) {
	*erridx = targetCellz;      
	*errsrc = coord[ii].z;
	errsrc[1] = 2;
	errsrc[2] = frameLow.z;
	errsrc[3] = dcellzi;
	errsrc[4] = numCell.z;
      }
      return;
    }
    IndexType cellid = CudaDevice::D3toD1
	(numCell, targetCellx, targetCelly, targetCellz);
    forwardMap[ii] = ii;
    if (cellid != bid){
      IndexType pid = atomicAdd (&numAtomInCell[cellid], 1);
      if (pid >= blockDim.x){
	*ptr_de = mdErrorShortCellList;
	pid = 0;
      }
      IndexType targetIndex = pid + cellid * blockDim.x;
      coord[targetIndex] = coord[ii];
      coordNoi[targetIndex].x = coordNoi[ii].x;
      coordNoi[targetIndex].y = coordNoi[ii].y;
      coordNoi[targetIndex].z = coordNoi[ii].z;
      velox[targetIndex] = velox[ii];
      veloy[targetIndex] = veloy[ii];
      veloz[targetIndex] = veloz[ii];
      forcx[targetIndex] = forcx[ii];
      forcy[targetIndex] = forcy[ii];
      forcz[targetIndex] = forcz[ii];
      globalIndex[targetIndex] = globalIndex[ii];
      globalIndex[ii] = MaxIndexValue;
      type[targetIndex] = type[ii];
      mass[targetIndex] = mass[ii];
      charge[targetIndex] = charge[ii];
      forwardMap[ii] = targetIndex;
    }
  }
  return;
}


__global__ void Parallel::CudaGlobal::
rebuildCellList_step2 (IndexType * numAtomInCell,
		       CoordType  * coord,
		       CoordNoiType * coordNoi,
		       ScalorType * velox,
		       ScalorType * veloy,
		       ScalorType * veloz,
		       ScalorType * forcx,
		       ScalorType * forcy,
		       ScalorType * forcz,
		       IndexType  * globalIndex,
		       TypeType   * type,
		       ScalorType * mass,
		       ScalorType * charge,
		       IndexType * forwardMap,
		       mdError_t * ptr_de)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  // IndexType k = NUintBit - 1;
  IndexType this_numAtomInCell = numAtomInCell[bid];
  if (tid < this_numAtomInCell) {
    forwardMap[ii] = ii;
  }
  else {
    forwardMap[ii] = MaxIndexValue;
  }
  __syncthreads();
  if (this_numAtomInCell == 0) return;
  
  extern __shared__ volatile IndexType sbuff[];
  volatile IndexType * myIndex = (volatile IndexType * )sbuff;

  if (tid >= this_numAtomInCell || globalIndex[ii] == MaxIndexValue){
    myIndex[tid] = MaxIndexValue;
  }
  else {
    myIndex[tid] = tid;
  }
  __syncthreads();
  IndexType total = headSort (myIndex, &sbuff[blockDim.x]);
  total = blockDim.x - total;

  IndexType fromId;
  CoordType bk_coord;
  CoordNoiType bk_coordNoi;
  ScalorType bk_velox, bk_veloy, bk_veloz;
  ScalorType bk_forcx, bk_forcy, bk_forcz;
  IndexType bk_globalIndex;
  TypeType bk_type;
  ScalorType bk_mass;
  ScalorType bk_charge;

  __syncthreads();
  
  if (tid < total){
    fromId = myIndex[tid] + bid * blockDim.x;
    if (ii != fromId){
      bk_coord = coord[fromId];
      bk_coordNoi.x = coordNoi[fromId].x;
      bk_coordNoi.y = coordNoi[fromId].y;
      bk_coordNoi.z = coordNoi[fromId].z;
      bk_velox = velox[fromId];
      bk_veloy = veloy[fromId];
      bk_veloz = veloz[fromId];
      bk_forcx = forcx[fromId];
      bk_forcy = forcy[fromId];
      bk_forcz = forcz[fromId];
      bk_globalIndex = globalIndex[fromId];
      bk_type = type[fromId];
      bk_mass = mass[fromId];
      bk_charge = charge[fromId];
      forwardMap[fromId] = ii;
    }    
  }
  __syncthreads();

  if (tid < total && ii != fromId){
    coord[ii] = bk_coord;
    coordNoi[ii].x = bk_coordNoi.x;
    coordNoi[ii].y = bk_coordNoi.y;
    coordNoi[ii].z = bk_coordNoi.z;
    velox[ii] = bk_velox;
    veloy[ii] = bk_veloy;
    veloz[ii] = bk_veloz;
    forcx[ii] = bk_forcx;
    forcy[ii] = bk_forcy;
    forcz[ii] = bk_forcz;
    globalIndex[ii] = bk_globalIndex;
    type[ii] = bk_type;
    mass[ii] = bk_mass;
    charge[ii] = bk_charge;
  }

  if (tid == 0){
    numAtomInCell[bid] = total;
  }
}


__global__ void Parallel::CudaGlobal::
rebuildCellList_step1_mapBondTop (const IndexType * forwardMap,
				  const IndexType * bk_numAtomInCell,
				  IndexType maxNumBond,
				  IndexType * numBond,
				  IndexType * bondIndex,
				  IndexType * bondNeighbor_globalIndex,
				  IndexType * bondNeighbor_localIndex,
				  IndexType maxNumAngle,
				  IndexType * numAngle,
				  IndexType * angleIndex,
				  IndexType * anglePosi,
				  IndexType * angleNeighbor_globalIndex,
				  IndexType * angleNeighbor_localIndex,
				  IndexType maxNumDihedral,
				  IndexType * numDihedral,
				  IndexType * dihedralIndex,
				  IndexType * dihedralPosi,
				  IndexType * dihedralNeighbor_globalIndex,
				  IndexType * dihedralNeighbor_localIndex,
				  IndexType bondTopStride)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType this_numAtomInCell = bk_numAtomInCell[bid];
  if (this_numAtomInCell == 0) return;
  IndexType fromIndex = tid + bid * blockDim.x;
  IndexType toIndex   = forwardMap[fromIndex];
  
  if (maxNumBond != 0){
    IndexType my_numBond;
    if (tid < this_numAtomInCell){
      my_numBond = (numBond[fromIndex]);
      IndexType my_index = fromIndex;
      for (IndexType i = 0; i < my_numBond; ++i){
      	IndexType tmp0 = bondNeighbor_localIndex[my_index];
      	IndexType tmp1 = forwardMap[tmp0];
      	if (tmp1 != MaxIndexValue){
      	  bondNeighbor_localIndex[my_index] = tmp1;
      	}    
      	my_index += bondTopStride;
      }
    }
    __syncthreads();
    if (toIndex != MaxIndexValue && toIndex != fromIndex){
      numBond[toIndex] = my_numBond;
      numBond[fromIndex] = 0;
      IndexType my_fromIndex = fromIndex;
      IndexType my_toIndex   = toIndex;
      for (IndexType i = 0; i < my_numBond; ++i){
	bondIndex[my_toIndex] = bondIndex[my_fromIndex];
	bondNeighbor_globalIndex[my_toIndex] =
	    bondNeighbor_globalIndex[my_fromIndex];
	bondNeighbor_localIndex [my_toIndex] =
	    bondNeighbor_localIndex [my_fromIndex];
	my_fromIndex += bondTopStride;
	my_toIndex   += bondTopStride;
      }
    }
  }

  if (maxNumAngle != 0){
    IndexType my_numAngle;
    if (tid < this_numAtomInCell){
      my_numAngle = (numAngle[fromIndex]);
      IndexType my_index = fromIndex;
      for (IndexType i = 0; i < my_numAngle * 2; ++i){
      	IndexType tmp0 = angleNeighbor_localIndex[my_index];
      	IndexType tmp1 = forwardMap[tmp0];
      	if (tmp1 != MaxIndexValue){
      	  angleNeighbor_localIndex[my_index] = tmp1;
      	}    
      	my_index += bondTopStride;
      }
    }
    __syncthreads();
    if (toIndex != MaxIndexValue && toIndex != fromIndex){
      numAngle[toIndex] = my_numAngle;
      numAngle[fromIndex] = 0;
      IndexType my_fromIndex = fromIndex;
      IndexType my_toIndex   = toIndex;
      for (IndexType i = 0; i < my_numAngle; ++i){
	angleIndex[my_toIndex] = angleIndex[my_fromIndex];
	anglePosi [my_toIndex] = anglePosi [my_fromIndex];
	my_fromIndex += bondTopStride;
	my_toIndex   += bondTopStride;
      }
      my_fromIndex = fromIndex;
      my_toIndex   = toIndex;
      for (IndexType i = 0; i < my_numAngle * 2; ++i){
	angleNeighbor_globalIndex[my_toIndex] =
	    angleNeighbor_globalIndex[my_fromIndex];
	angleNeighbor_localIndex [my_toIndex] =
	    angleNeighbor_localIndex [my_fromIndex];
	my_fromIndex += bondTopStride;
	my_toIndex   += bondTopStride;
      }
    }
  }

  if (maxNumDihedral != 0){
    IndexType my_numDihedral;
    if (tid < this_numAtomInCell){
      my_numDihedral = (numDihedral[fromIndex]);
      IndexType my_index = fromIndex;
      for (IndexType i = 0; i < my_numDihedral * 3; ++i){
      	IndexType tmp0 = dihedralNeighbor_localIndex[my_index];
      	IndexType tmp1 = forwardMap[tmp0];
      	if (tmp1 != MaxIndexValue){
      	  dihedralNeighbor_localIndex[my_index] = tmp1;
      	}    
      	my_index += bondTopStride;
      }
    }
    __syncthreads();
    if (toIndex != MaxIndexValue && toIndex != fromIndex){
      numDihedral[toIndex] = my_numDihedral;
      numDihedral[fromIndex] = 0;
      IndexType my_fromIndex = fromIndex;
      IndexType my_toIndex   = toIndex;
      for (IndexType i = 0; i < my_numDihedral; ++i){
	dihedralIndex[my_toIndex] = dihedralIndex[my_fromIndex];
	dihedralPosi [my_toIndex] = dihedralPosi [my_fromIndex];
	my_fromIndex += bondTopStride;
	my_toIndex   += bondTopStride;
      }
      my_fromIndex = fromIndex;
      my_toIndex   = toIndex;
      for (IndexType i = 0; i < my_numDihedral * 3; ++i){
	dihedralNeighbor_globalIndex[my_toIndex] =
	    dihedralNeighbor_globalIndex[my_fromIndex];
	dihedralNeighbor_localIndex [my_toIndex] =
	    dihedralNeighbor_localIndex [my_fromIndex];
	my_fromIndex += bondTopStride;
	my_toIndex   += bondTopStride;
      }
    }
  }

      
  // if (toIndex != MaxIndexValue){
  //   if (maxNumBond != 0){
  //     IndexType my_numBond = (numBond[toIndex] = numBond[fromIndex]);
  //     IndexType my_fromIndex = fromIndex;
  //     IndexType my_toIndex   = toIndex;
  //     for (IndexType i = 0; i < my_numBond; ++i){
  // 	bondIndex[my_toIndex] = bondIndex[my_fromIndex];
  // 	bondNeighbor_globalIndex[my_toIndex] =
  // 	    bondNeighbor_globalIndex[my_fromIndex];
  // 	IndexType tmpIndex0 = bondNeighbor_localIndex[my_fromIndex];
  // 	IndexType tmpIndex1 = forwardMap[tmpIndex0];
  // 	if (tmpIndex1 == MaxIndexValue){
  // 	  bondNeighbor_localIndex [my_toIndex] = tmpIndex0;
  // 	}
  // 	else {
  // 	  bondNeighbor_localIndex [my_toIndex] = tmpIndex1;
  // 	}
  // 	my_fromIndex += bondTopStride;
  // 	my_toIndex   += bondTopStride;
  //     }
  //   }
  //   if (maxNumAngle != 0){
  //     IndexType my_numAngle = (numAngle[toIndex] = numAngle[fromIndex]);
  //     IndexType my_fromIndex = fromIndex;
  //     IndexType my_toIndex   = toIndex;
  //     for (IndexType i = 0; i < my_numAngle; ++i){
  // 	angleIndex[my_toIndex] = angleIndex[my_fromIndex];
  // 	anglePosi [my_toIndex] = anglePosi [my_fromIndex];
  // 	my_fromIndex += bondTopStride;
  // 	my_toIndex   += bondTopStride;
  //     }
  //     my_fromIndex = fromIndex;
  //     my_toIndex   = toIndex;
  //     for (IndexType i = 0; i < my_numAngle * 2; ++i){      
  // 	angleNeighbor_globalIndex[my_toIndex] =
  // 	    angleNeighbor_globalIndex[my_fromIndex];
  // 	IndexType tmpIndex0 = angleNeighbor_localIndex[my_fromIndex];
  // 	IndexType tmpIndex1 = forwardMap[tmpIndex0];
  // 	if (tmpIndex1 == MaxIndexValue){
  // 	  angleNeighbor_localIndex [my_toIndex] = tmpIndex0;
  // 	}
  // 	else {
  // 	  angleNeighbor_localIndex [my_toIndex] = tmpIndex1;
  // 	}
  // 	my_fromIndex += bondTopStride;
  // 	my_toIndex   += bondTopStride;
  //     }
  //   }
  //   if (maxNumDihedral != 0){
  //     IndexType my_numDihedral = (numDihedral[toIndex] = numDihedral[fromIndex]);
  //     IndexType my_fromIndex = fromIndex;
  //     IndexType my_toIndex   = toIndex;
  //     for (IndexType i = 0; i < my_numDihedral; ++i){
  // 	dihedralIndex[my_toIndex] = dihedralIndex[my_fromIndex];
  // 	dihedralPosi [my_toIndex] = dihedralPosi [my_fromIndex];
  // 	my_fromIndex += bondTopStride;
  // 	my_toIndex   += bondTopStride;
  //     }
  //     my_fromIndex = fromIndex;
  //     my_toIndex   = toIndex;
  //     for (IndexType i = 0; i < my_numDihedral * 3; ++i){      
  // 	dihedralNeighbor_globalIndex[my_toIndex] =
  // 	    dihedralNeighbor_globalIndex[my_fromIndex];
  // 	IndexType tmpIndex0 = dihedralNeighbor_localIndex[my_fromIndex];
  // 	IndexType tmpIndex1 = forwardMap[tmpIndex0];
  // 	if (tmpIndex1 == MaxIndexValue){
  // 	  dihedralNeighbor_localIndex [my_toIndex] = tmpIndex0;
  // 	}
  // 	else {
  // 	  dihedralNeighbor_localIndex [my_toIndex] = tmpIndex1;
  // 	}
  // 	my_fromIndex += bondTopStride;
  // 	my_toIndex   += bondTopStride;
  //     }
  //   }
  // }
}

__global__ void Parallel::CudaGlobal::
rebuildCellList_step2_mapBondTop (const IndexType * forwardMap,
				  const IndexType * bk_numAtomInCell,
				  IndexType maxNumBond,
				  IndexType * numBond,
				  IndexType * bondIndex,
				  IndexType * bondNeighbor_globalIndex,
				  IndexType * bondNeighbor_localIndex,
				  IndexType maxNumAngle,
				  IndexType * numAngle,
				  IndexType * angleIndex,
				  IndexType * anglePosi,
				  IndexType * angleNeighbor_globalIndex,
				  IndexType * angleNeighbor_localIndex,
				  IndexType maxNumDihedral,
				  IndexType * numDihedral,
				  IndexType * dihedralIndex,
				  IndexType * dihedralPosi,
				  IndexType * dihedralNeighbor_globalIndex,
				  IndexType * dihedralNeighbor_localIndex,
				  IndexType bondTopStride)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType this_numAtomInCell = bk_numAtomInCell[bid];
  if (this_numAtomInCell == 0) return;
  IndexType fromIndex = tid + bid * blockDim.x;
  IndexType toIndex = forwardMap[fromIndex];
  bool docopy = (toIndex != MaxIndexValue && fromIndex != toIndex);

  IndexType bk_num, bk_Index, bk_Posi;
  IndexType bk_Neighbor_globalIndex, bk_Neighbor_localIndex;

  if (maxNumBond != 0){
    IndexType my_numBond;
    if (tid < this_numAtomInCell){
      my_numBond = (numBond[fromIndex]);
      IndexType my_index = fromIndex;
      for (IndexType i = 0; i < my_numBond; ++i){
      	IndexType tmp0 = bondNeighbor_localIndex[my_index];
      	IndexType tmp1 = forwardMap[tmp0];
      	if (tmp1 != MaxIndexValue){
      	  bondNeighbor_localIndex[my_index] = tmp1;
      	}
      	my_index += bondTopStride;
      }
    }
    __syncthreads();
    if (maxNumBond != 0){
      if (docopy){
	bk_num = my_numBond;
      }
      __syncthreads();
      if (docopy){
	numBond[toIndex] = bk_num;
      }
      __syncthreads();
      IndexType my_fromIndex = fromIndex;
      IndexType my_toIndex = toIndex;
      for (IndexType kk = 0; kk < maxNumBond; ++kk){
	if (docopy){
	  bk_Index = bondIndex[my_fromIndex];
	  bk_Neighbor_globalIndex = bondNeighbor_globalIndex[my_fromIndex];
	  bk_Neighbor_localIndex  = bondNeighbor_localIndex [my_fromIndex];
	}
	__syncthreads();
	if (docopy){
	  bondIndex[my_toIndex] = bk_Index;
	  bondNeighbor_globalIndex[my_toIndex] = bk_Neighbor_globalIndex;
	  bondNeighbor_localIndex [my_toIndex] = bk_Neighbor_localIndex;
	}
	__syncthreads();
	my_fromIndex += bondTopStride;
	my_toIndex   += bondTopStride;
      }
    }
  }

  if (maxNumAngle != 0){
    IndexType my_numAngle;
    if (tid < this_numAtomInCell){
      my_numAngle = (numAngle[fromIndex]);
      IndexType my_index = fromIndex;
      for (IndexType i = 0; i < my_numAngle * 2; ++i){
      	IndexType tmp0 = angleNeighbor_localIndex[my_index];
      	IndexType tmp1 = forwardMap[tmp0];
      	if (tmp1 != MaxIndexValue){
      	  angleNeighbor_localIndex[my_index] = tmp1;
      	}
      	my_index += bondTopStride;
      }
    }
    __syncthreads();
    if (maxNumAngle != 0){
      if (docopy){
	bk_num = my_numAngle;
      }
      __syncthreads();
      if (docopy){
	numAngle[toIndex] = bk_num;
      }
      __syncthreads();
      IndexType my_fromIndex = fromIndex;
      IndexType my_toIndex = toIndex;
      for (IndexType kk = 0; kk < maxNumAngle; ++kk){
	if (docopy){
	  bk_Index = angleIndex[my_fromIndex];
	  bk_Posi  = anglePosi [my_fromIndex];
	}
	__syncthreads();
	if (docopy){
	  angleIndex[my_toIndex] = bk_Index;
	  anglePosi [my_toIndex] = bk_Posi ;
	}
	__syncthreads();
	my_fromIndex += bondTopStride;
	my_toIndex   += bondTopStride;
      }
      my_fromIndex = fromIndex;
      my_toIndex = toIndex;
      for (IndexType kk = 0; kk < maxNumAngle * 2; ++kk){
	if (docopy){
	  bk_Neighbor_globalIndex = angleNeighbor_globalIndex[my_fromIndex];
	  bk_Neighbor_localIndex  = angleNeighbor_localIndex [my_fromIndex];
	}
	__syncthreads();
	if (docopy){
	  angleNeighbor_globalIndex[my_toIndex] = bk_Neighbor_globalIndex;
	  angleNeighbor_localIndex [my_toIndex] = bk_Neighbor_localIndex;
	}
	__syncthreads();
	my_fromIndex += bondTopStride;
	my_toIndex   += bondTopStride;
      }
    }
  }

  
  if (maxNumDihedral != 0){
    IndexType my_numDihedral;
    if (tid < this_numAtomInCell){
      my_numDihedral = (numDihedral[fromIndex]);
      IndexType my_index = fromIndex;
      for (IndexType i = 0; i < my_numDihedral * 3; ++i){
      	IndexType tmp0 = dihedralNeighbor_localIndex[my_index];
      	IndexType tmp1 = forwardMap[tmp0];
      	if (tmp1 != MaxIndexValue){
      	  dihedralNeighbor_localIndex[my_index] = tmp1;
      	}
      	my_index += bondTopStride;
      }
    }
    __syncthreads();
    if (maxNumDihedral != 0){
      if (docopy){
	bk_num = my_numDihedral;
      }
      __syncthreads();
      if (docopy){
	numDihedral[toIndex] = bk_num;
      }
      __syncthreads();
      IndexType my_fromIndex = fromIndex;
      IndexType my_toIndex = toIndex;
      for (IndexType kk = 0; kk < maxNumDihedral; ++kk){
	if (docopy){
	  bk_Index = dihedralIndex[my_fromIndex];
	  bk_Posi  = dihedralPosi [my_fromIndex];
	}
	__syncthreads();
	if (docopy){
	  dihedralIndex[my_toIndex] = bk_Index;
	  dihedralPosi [my_toIndex] = bk_Posi ;
	}
	__syncthreads();
	my_fromIndex += bondTopStride;
	my_toIndex   += bondTopStride;
      }
      my_fromIndex = fromIndex;
      my_toIndex = toIndex;
      for (IndexType kk = 0; kk < maxNumDihedral * 3; ++kk){
	if (docopy){
	  bk_Neighbor_globalIndex = dihedralNeighbor_globalIndex[my_fromIndex];
	  bk_Neighbor_localIndex  = dihedralNeighbor_localIndex [my_fromIndex];
	}
	__syncthreads();
	if (docopy){
	  dihedralNeighbor_globalIndex[my_toIndex] = bk_Neighbor_globalIndex;
	  dihedralNeighbor_localIndex [my_toIndex] = bk_Neighbor_localIndex;
	}
	__syncthreads();
	my_fromIndex += bondTopStride;
	my_toIndex   += bondTopStride;
      }
    }
  }


  // if (maxNumBond != 0){
  //   if (docopy){
  //     bk_num = numBond[fromIndex];
  //   }
  //   __syncthreads();
  //   if (docopy){
  //     numBond[toIndex] = bk_num;
  //   }
  //   __syncthreads();
  //   IndexType my_fromIndex = fromIndex;
  //   IndexType my_toIndex = toIndex;
  //   for (IndexType kk = 0; kk < maxNumBond; ++kk){
  //     if (docopy){
  // 	bk_Index = bondIndex[my_fromIndex];
  // 	bk_Neighbor_globalIndex = bondNeighbor_globalIndex[my_fromIndex];
  // 	IndexType tmpIndex0 = bondNeighbor_localIndex[my_fromIndex];
  // 	IndexType tmpIndex1 = forwardMap[tmpIndex0];
  // 	if (tmpIndex1 == MaxIndexValue){
  // 	  bk_Neighbor_localIndex = tmpIndex0;
  // 	}
  // 	else {
  // 	  bk_Neighbor_localIndex = tmpIndex1;
  // 	}
  //     }
  //     __syncthreads();
  //     if (docopy){
  // 	bondIndex[my_toIndex] = bk_Index;
  // 	bondNeighbor_globalIndex[my_toIndex] = bk_Neighbor_globalIndex;
  // 	bondNeighbor_localIndex [my_toIndex] = bk_Neighbor_localIndex;
  //     }
  //     __syncthreads();
  //     my_fromIndex += bondTopStride;
  //     my_toIndex   += bondTopStride;
  //   }
  // }
  
  // if (maxNumAngle != 0){
  //   if (docopy){
  //     bk_num = numAngle[fromIndex];
  //   }
  //   __syncthreads();
  //   if (docopy){
  //     numAngle[toIndex] = bk_num;
  //   }
  //   __syncthreads();
  //   IndexType my_fromIndex = fromIndex;
  //   IndexType my_toIndex = toIndex;
  //   for (IndexType kk = 0; kk < maxNumAngle; ++kk){
  //     if (docopy){
  // 	bk_Index = angleIndex[my_fromIndex];
  // 	bk_Posi  = anglePosi [my_fromIndex];
  //     }
  //     __syncthreads();
  //     if (docopy){
  // 	angleIndex[my_toIndex] = bk_Index;
  // 	anglePosi [my_toIndex] = bk_Posi ;
  //     }
  //     __syncthreads();
  //     my_fromIndex += bondTopStride;
  //     my_toIndex   += bondTopStride;
  //   }
  //   my_fromIndex = fromIndex;
  //   my_toIndex = toIndex;
  //   for (IndexType kk = 0; kk < maxNumAngle * 2; ++kk){
  //     if (docopy){
  // 	bk_Neighbor_globalIndex = angleNeighbor_globalIndex[my_fromIndex];
  // 	IndexType tmpIndex0 = angleNeighbor_localIndex[my_fromIndex];
  // 	IndexType tmpIndex1 = forwardMap[tmpIndex0];
  // 	if (tmpIndex1 == MaxIndexValue){
  // 	  bk_Neighbor_localIndex = tmpIndex0;
  // 	}
  // 	else {
  // 	  bk_Neighbor_localIndex = tmpIndex1;
  // 	}
  //     }
  //     __syncthreads();
  //     if (docopy){
  // 	angleNeighbor_globalIndex[my_toIndex] = bk_Neighbor_globalIndex;
  // 	angleNeighbor_localIndex [my_toIndex] = bk_Neighbor_localIndex;
  //     }
  //     __syncthreads();
  //     my_fromIndex += bondTopStride;
  //     my_toIndex   += bondTopStride;
  //   }
  // }

  // if (maxNumDihedral != 0){
  //   if (docopy){
  //     bk_num = numDihedral[fromIndex];
  //   }
  //   __syncthreads();
  //   if (docopy){
  //     numDihedral[toIndex] = bk_num;
  //   }
  //   __syncthreads();
  //   IndexType my_fromIndex = fromIndex;
  //   IndexType my_toIndex = toIndex;
  //   for (IndexType kk = 0; kk < maxNumDihedral; ++kk){
  //     if (docopy){
  // 	bk_Index = dihedralIndex[my_fromIndex];
  // 	bk_Posi  = dihedralPosi [my_fromIndex];
  //     }
  //     __syncthreads();
  //     if (docopy){
  // 	dihedralIndex[my_toIndex] = bk_Index;
  // 	dihedralPosi [my_toIndex] = bk_Posi ;
  //     }
  //     __syncthreads();
  //     my_fromIndex += bondTopStride;
  //     my_toIndex   += bondTopStride;
  //   }
  //   my_fromIndex = fromIndex;
  //   my_toIndex = toIndex;
  //   for (IndexType kk = 0; kk < maxNumDihedral * 3; ++kk){
  //     if (docopy){
  // 	bk_Neighbor_globalIndex = dihedralNeighbor_globalIndex[my_fromIndex];
  // 	IndexType tmpIndex0 = dihedralNeighbor_localIndex[my_fromIndex];
  // 	IndexType tmpIndex1 = forwardMap[tmpIndex0];
  // 	if (tmpIndex1 == MaxIndexValue){
  // 	  bk_Neighbor_localIndex = tmpIndex0;
  // 	}
  // 	else {
  // 	  bk_Neighbor_localIndex = tmpIndex1;
  // 	}
  //     }
  //     __syncthreads();
  //     if (docopy){
  // 	dihedralNeighbor_globalIndex[my_toIndex] = bk_Neighbor_globalIndex;
  // 	dihedralNeighbor_localIndex [my_toIndex] = bk_Neighbor_localIndex;
  //     }
  //     __syncthreads();
  //     my_fromIndex += bondTopStride;
  //     my_toIndex   += bondTopStride;
  //   }
  // }
}




__global__ void Parallel::CudaGlobal::
rebuildCellList_step1_mapBondTop (const IndexType * forwardMap,
				  const IndexType * bk_numAtomInCell,
				  IndexType maxNumBond,
				  IndexType * numBond,
				  IndexType * bondIndex,
				  IndexType * bondNeighbor_globalIndex,
				  IndexType maxNumAngle,
				  IndexType * numAngle,
				  IndexType * angleIndex,
				  IndexType * anglePosi,
				  IndexType * angleNeighbor_globalIndex,
				  IndexType maxNumDihedral,
				  IndexType * numDihedral,
				  IndexType * dihedralIndex,
				  IndexType * dihedralPosi,
				  IndexType * dihedralNeighbor_globalIndex,
				  IndexType bondTopStride)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType this_numAtomInCell = bk_numAtomInCell[bid];
  if (this_numAtomInCell == 0) return;
  IndexType fromIndex = tid + bid * blockDim.x;
  IndexType toIndex   = forwardMap[fromIndex];
  
  if (maxNumBond != 0){
    if (toIndex != MaxIndexValue && toIndex != fromIndex){
      IndexType my_numBond = numBond[toIndex] = numBond[fromIndex];
      numBond[fromIndex] = 0;
      IndexType my_fromIndex = fromIndex;
      IndexType my_toIndex   = toIndex;
      for (IndexType i = 0; i < my_numBond; ++i){
	bondIndex[my_toIndex] = bondIndex[my_fromIndex];
	bondNeighbor_globalIndex[my_toIndex] =
	    bondNeighbor_globalIndex[my_fromIndex];
	my_fromIndex += bondTopStride;
	my_toIndex   += bondTopStride;
      }
    }
  }

  if (maxNumAngle != 0){
    if (toIndex != MaxIndexValue && toIndex != fromIndex){
      IndexType my_numAngle = numAngle[toIndex] = numAngle[fromIndex];
      numAngle[fromIndex] = 0;
      IndexType my_fromIndex = fromIndex;
      IndexType my_toIndex   = toIndex;
      for (IndexType i = 0; i < my_numAngle; ++i){
	angleIndex[my_toIndex] = angleIndex[my_fromIndex];
	anglePosi [my_toIndex] = anglePosi [my_fromIndex];
	my_fromIndex += bondTopStride;
	my_toIndex   += bondTopStride;
      }
      my_fromIndex = fromIndex;
      my_toIndex   = toIndex;
      for (IndexType i = 0; i < my_numAngle * 2; ++i){
	angleNeighbor_globalIndex[my_toIndex] =
	    angleNeighbor_globalIndex[my_fromIndex];
	my_fromIndex += bondTopStride;
	my_toIndex   += bondTopStride;
      }
    }
  }

  if (maxNumDihedral != 0){
    if (toIndex != MaxIndexValue && toIndex != fromIndex){
      IndexType my_numDihedral = numDihedral[toIndex] = numDihedral[fromIndex];
      numDihedral[fromIndex] = 0;
      IndexType my_fromIndex = fromIndex;
      IndexType my_toIndex   = toIndex;
      for (IndexType i = 0; i < my_numDihedral; ++i){
	dihedralIndex[my_toIndex] = dihedralIndex[my_fromIndex];
	dihedralPosi [my_toIndex] = dihedralPosi [my_fromIndex];
	my_fromIndex += bondTopStride;
	my_toIndex   += bondTopStride;
      }
      my_fromIndex = fromIndex;
      my_toIndex   = toIndex;
      for (IndexType i = 0; i < my_numDihedral * 3; ++i){
	dihedralNeighbor_globalIndex[my_toIndex] =
	    dihedralNeighbor_globalIndex[my_fromIndex];
	my_fromIndex += bondTopStride;
	my_toIndex   += bondTopStride;
      }
    }
  }
}


__global__ void Parallel::CudaGlobal::
rebuildCellList_step2_mapBondTop (const IndexType * forwardMap,
				  const IndexType * bk_numAtomInCell,
				  IndexType maxNumBond,
				  IndexType * numBond,
				  IndexType * bondIndex,
				  IndexType * bondNeighbor_globalIndex,
				  IndexType maxNumAngle,
				  IndexType * numAngle,
				  IndexType * angleIndex,
				  IndexType * anglePosi,
				  IndexType * angleNeighbor_globalIndex,
				  IndexType maxNumDihedral,
				  IndexType * numDihedral,
				  IndexType * dihedralIndex,
				  IndexType * dihedralPosi,
				  IndexType * dihedralNeighbor_globalIndex,
				  IndexType bondTopStride)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType this_numAtomInCell = bk_numAtomInCell[bid];
  if (this_numAtomInCell == 0) return;
  IndexType fromIndex = tid + bid * blockDim.x;
  IndexType toIndex = forwardMap[fromIndex];
  bool docopy = (toIndex != MaxIndexValue && fromIndex != toIndex);

  IndexType bk_num, bk_Index, bk_Posi;
  IndexType bk_Neighbor_globalIndex;

  if (maxNumBond != 0){
    if (docopy){
      bk_num = (numBond[fromIndex]);
    }
    __syncthreads();
    if (docopy){
      numBond[toIndex] = bk_num;
    }
    __syncthreads();
    IndexType my_fromIndex = fromIndex;
    IndexType my_toIndex = toIndex;
    for (IndexType kk = 0; kk < maxNumBond; ++kk){
      if (docopy){
	bk_Index = bondIndex[my_fromIndex];
	bk_Neighbor_globalIndex = bondNeighbor_globalIndex[my_fromIndex];
      }
      __syncthreads();
      if (docopy){
	bondIndex[my_toIndex] = bk_Index;
	bondNeighbor_globalIndex[my_toIndex] = bk_Neighbor_globalIndex;
      }
      __syncthreads();
      my_fromIndex += bondTopStride;
      my_toIndex   += bondTopStride;
    }
  }

  if (maxNumAngle != 0){
    if (docopy){
      bk_num = (numAngle[fromIndex]);
    }
    __syncthreads();
    if (docopy){
      numAngle[toIndex] = bk_num;
    }
    __syncthreads();
    IndexType my_fromIndex = fromIndex;
    IndexType my_toIndex = toIndex;
    for (IndexType kk = 0; kk < maxNumAngle; ++kk){
      if (docopy){
	bk_Index = angleIndex[my_fromIndex];
	bk_Posi  = anglePosi [my_fromIndex];
      }
      __syncthreads();
      if (docopy){
	angleIndex[my_toIndex] = bk_Index;
	anglePosi [my_toIndex] = bk_Posi ;
      }
      __syncthreads();
      my_fromIndex += bondTopStride;
      my_toIndex   += bondTopStride;
    }
    my_fromIndex = fromIndex;
    my_toIndex = toIndex;
    for (IndexType kk = 0; kk < maxNumAngle * 2; ++kk){
      if (docopy){
	bk_Neighbor_globalIndex = angleNeighbor_globalIndex[my_fromIndex];
      }
      __syncthreads();
      if (docopy){
	angleNeighbor_globalIndex[my_toIndex] = bk_Neighbor_globalIndex;
      }
      __syncthreads();
      my_fromIndex += bondTopStride;
      my_toIndex   += bondTopStride;
    }
  }

  if (maxNumDihedral != 0){
    if (docopy){
      bk_num = (numDihedral[fromIndex]);
    }
    __syncthreads();
    if (docopy){
      numDihedral[toIndex] = bk_num;
    }
    __syncthreads();
    IndexType my_fromIndex = fromIndex;
    IndexType my_toIndex = toIndex;
    for (IndexType kk = 0; kk < maxNumDihedral; ++kk){
      if (docopy){
	bk_Index = dihedralIndex[my_fromIndex];
	bk_Posi  = dihedralPosi [my_fromIndex];
      }
      __syncthreads();
      if (docopy){
	dihedralIndex[my_toIndex] = bk_Index;
	dihedralPosi [my_toIndex] = bk_Posi ;
      }
      __syncthreads();
      my_fromIndex += bondTopStride;
      my_toIndex   += bondTopStride;
    }
    my_fromIndex = fromIndex;
    my_toIndex = toIndex;
    for (IndexType kk = 0; kk < maxNumDihedral * 3; ++kk){
      if (docopy){
	bk_Neighbor_globalIndex = dihedralNeighbor_globalIndex[my_fromIndex];
      }
      __syncthreads();
      if (docopy){
	dihedralNeighbor_globalIndex[my_toIndex] = bk_Neighbor_globalIndex;
      }
      __syncthreads();
      my_fromIndex += bondTopStride;
      my_toIndex   += bondTopStride;
    }
  }
}



// __global__ void Parallel::CudaGlobal::
// rebuildCellList_step1_mapBondTop (const IndexType * forwardMap,
// 				  IndexType maxNumBond,
// 				  IndexType * numBond,
// 				  IndexType * bondIndex,
// 				  IndexType * bondNeighbor_globalIndex,
// 				  IndexType maxNumAngle,
// 				  IndexType * numAngle,
// 				  IndexType * angleIndex,
// 				  IndexType * anglePosi,
// 				  IndexType * angleNeighbor_globalIndex,
// 				  IndexType maxNumDihedral,
// 				  IndexType * numDihedral,
// 				  IndexType * dihedralIndex,
// 				  IndexType * dihedralPosi,
// 				  IndexType * dihedralNeighbor_globalIndex,
// 				  IndexType bondTopStride)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   IndexType fromIndex = tid + bid * blockDim.x;
//   IndexType toIndex   = forwardMap[fromIndex];
  
//   if (toIndex != MaxIndexValue){
//     if (maxNumBond != 0){
//       IndexType my_numBond = (numBond[toIndex] = numBond[fromIndex]);
//       IndexType my_fromIndex = fromIndex;
//       IndexType my_toIndex   = toIndex;
//       for (IndexType i = 0; i < my_numBond; ++i){
// 	bondIndex[my_toIndex] = bondIndex[my_fromIndex];
// 	bondNeighbor_globalIndex[my_toIndex] =
// 	    bondNeighbor_globalIndex[my_fromIndex];
// 	my_fromIndex += bondTopStride;
// 	my_toIndex   += bondTopStride;
//       }
//     }
//     if (maxNumAngle != 0){
//       IndexType my_numAngle = (numAngle[toIndex] = numAngle[fromIndex]);
//       IndexType my_fromIndex = fromIndex;
//       IndexType my_toIndex   = toIndex;
//       for (IndexType i = 0; i < my_numAngle; ++i){
// 	angleIndex[my_toIndex] = angleIndex[my_fromIndex];
// 	anglePosi [my_toIndex] = anglePosi [my_fromIndex];
// 	my_fromIndex += bondTopStride;
// 	my_toIndex   += bondTopStride;
//       }
//       my_fromIndex = fromIndex;
//       my_toIndex   = toIndex;
//       for (IndexType i = 0; i < my_numAngle * 2; ++i){      
// 	angleNeighbor_globalIndex[my_toIndex] =
// 	    angleNeighbor_globalIndex[my_fromIndex];
// 	my_fromIndex += bondTopStride;
// 	my_toIndex   += bondTopStride;
//       }
//     }
//     if (maxNumDihedral != 0){
//       IndexType my_numDihedral = (numDihedral[toIndex] = numDihedral[fromIndex]);
//       IndexType my_fromIndex = fromIndex;
//       IndexType my_toIndex   = toIndex;
//       for (IndexType i = 0; i < my_numDihedral; ++i){
// 	dihedralIndex[my_toIndex] = dihedralIndex[my_fromIndex];
// 	dihedralPosi [my_toIndex] = dihedralPosi [my_fromIndex];
// 	my_fromIndex += bondTopStride;
// 	my_toIndex   += bondTopStride;
//       }
//       my_fromIndex = fromIndex;
//       my_toIndex   = toIndex;
//       for (IndexType i = 0; i < my_numDihedral * 3; ++i){      
// 	dihedralNeighbor_globalIndex[my_toIndex] =
// 	    dihedralNeighbor_globalIndex[my_fromIndex];
// 	my_fromIndex += bondTopStride;
// 	my_toIndex   += bondTopStride;
//       }
//     }
//   }
// }

// __global__ void Parallel::CudaGlobal::
// rebuildCellList_step2_mapBondTop (const IndexType * forwardMap,
// 				  IndexType maxNumBond,
// 				  IndexType * numBond,
// 				  IndexType * bondIndex,
// 				  IndexType * bondNeighbor_globalIndex,
// 				  IndexType maxNumAngle,
// 				  IndexType * numAngle,
// 				  IndexType * angleIndex,
// 				  IndexType * anglePosi,
// 				  IndexType * angleNeighbor_globalIndex,
// 				  IndexType maxNumDihedral,
// 				  IndexType * numDihedral,
// 				  IndexType * dihedralIndex,
// 				  IndexType * dihedralPosi,
// 				  IndexType * dihedralNeighbor_globalIndex,
// 				  IndexType bondTopStride)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   IndexType fromIndex = tid + bid * blockDim.x;
//   IndexType toIndex = forwardMap[fromIndex];
//   bool docopy = (toIndex != MaxIndexValue);


//   IndexType bk_num, bk_Index, bk_Posi;
//   IndexType bk_Neighbor_globalIndex;
//   if (maxNumBond != 0){
//     if (docopy){
//       bk_num = numBond[fromIndex];
//     }
//     __syncthreads();
//     if (docopy){
//       numBond[toIndex] = bk_num;
//     }
//     __syncthreads();
//     IndexType my_fromIndex = fromIndex;
//     IndexType my_toIndex = toIndex;
//     for (IndexType kk = 0; kk < maxNumBond; ++kk){
//       if (docopy){
// 	bk_Index = bondIndex[my_fromIndex];
// 	bk_Neighbor_globalIndex = bondNeighbor_globalIndex[my_fromIndex];
//       }
//       __syncthreads();
//       if (docopy){
// 	bondIndex[my_toIndex] = bk_Index;
// 	bondNeighbor_globalIndex[my_toIndex] = bk_Neighbor_globalIndex;
//       }
//       __syncthreads();
//       my_fromIndex += bondTopStride;
//       my_toIndex   += bondTopStride;
//     }
//   }
  
//   if (maxNumAngle != 0){
//     if (docopy){
//       bk_num = numAngle[fromIndex];
//     }
//     __syncthreads();
//     if (docopy){
//       numAngle[toIndex] = bk_num;
//     }
//     __syncthreads();
//     IndexType my_fromIndex = fromIndex;
//     IndexType my_toIndex = toIndex;
//     for (IndexType kk = 0; kk < maxNumAngle; ++kk){
//       if (docopy){
// 	bk_Index = angleIndex[my_fromIndex];
// 	bk_Posi  = anglePosi [my_fromIndex];
//       }
//       __syncthreads();
//       if (docopy){
// 	angleIndex[my_toIndex] = bk_Index;
// 	anglePosi [my_toIndex] = bk_Posi ;
//       }
//       __syncthreads();
//       my_fromIndex += bondTopStride;
//       my_toIndex   += bondTopStride;
//     }
//     my_fromIndex = fromIndex;
//     my_toIndex = toIndex;
//     for (IndexType kk = 0; kk < maxNumAngle * 2; ++kk){
//       if (docopy){
// 	bk_Neighbor_globalIndex = angleNeighbor_globalIndex[my_fromIndex];
//       }
//       __syncthreads();
//       if (docopy){
// 	angleNeighbor_globalIndex[my_toIndex] = bk_Neighbor_globalIndex;
//       }
//       __syncthreads();
//       my_fromIndex += bondTopStride;
//       my_toIndex   += bondTopStride;
//     }
//   }

//   if (maxNumDihedral != 0){
//     if (docopy){
//       bk_num = numDihedral[fromIndex];
//     }
//     __syncthreads();
//     if (docopy){
//       numDihedral[toIndex] = bk_num;
//     }
//     __syncthreads();
//     IndexType my_fromIndex = fromIndex;
//     IndexType my_toIndex = toIndex;
//     for (IndexType kk = 0; kk < maxNumDihedral; ++kk){
//       if (docopy){
// 	bk_Index = dihedralIndex[my_fromIndex];
// 	bk_Posi  = dihedralPosi [my_fromIndex];
//       }
//       __syncthreads();
//       if (docopy){
// 	dihedralIndex[my_toIndex] = bk_Index;
// 	dihedralPosi [my_toIndex] = bk_Posi ;
//       }
//       __syncthreads();
//       my_fromIndex += bondTopStride;
//       my_toIndex   += bondTopStride;
//     }
//     my_fromIndex = fromIndex;
//     my_toIndex = toIndex;
//     for (IndexType kk = 0; kk < maxNumDihedral * 3; ++kk){
//       if (docopy){
// 	bk_Neighbor_globalIndex = dihedralNeighbor_globalIndex[my_fromIndex];
//       }
//       __syncthreads();
//       if (docopy){
// 	dihedralNeighbor_globalIndex[my_toIndex] = bk_Neighbor_globalIndex;
//       }
//       __syncthreads();
//       my_fromIndex += bondTopStride;
//       my_toIndex   += bondTopStride;
//     }
//   }
// }


// __global__ void Parallel::CudaGlobal::
// rebuildCellList_step1 (const VectorType frameLow,
// 		       const VectorType frameUp,
// 		       const IntVectorType numCell,
// 		       const IndexType * bk_numAtomInCell,
// 		       IndexType * numAtomInCell,
// 		       CoordType * coord,
// 		       CoordNoiType * coordNoi,
// 		       ScalorType * velox,
// 		       ScalorType * veloy,
// 		       ScalorType * veloz,
// 		       ScalorType * forcx,
// 		       ScalorType * forcy,
// 		       ScalorType * forcz,
// 		       IndexType  * globalIndex,
// 		       TypeType   * type,
// 		       ScalorType * mass,
// 		       ScalorType * charge,
// 		       mdError_t * ptr_de)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   IndexType ii = tid + bid * blockDim.x;

//   __shared__ ScalorType dcellxi;
//   __shared__ ScalorType dcellyi;
//   __shared__ ScalorType dcellzi;
//   if (tid == 0) {
//     dcellxi = ScalorType(numCell.x) / (frameUp.x - frameLow.x);
//     dcellyi = ScalorType(numCell.y) / (frameUp.y - frameLow.y);
//     dcellzi = ScalorType(numCell.z) / (frameUp.z - frameLow.z);
//   }
//   __syncthreads();
  
//   bool inRange = tid < bk_numAtomInCell[bid];
//   if (inRange){
//     IndexType targetCellx, targetCelly, targetCellz;
//     targetCellx = IndexType((coord[ii].x - frameLow.x) * dcellxi);
//     targetCelly = IndexType((coord[ii].y - frameLow.y) * dcellyi);
//     targetCellz = IndexType((coord[ii].z - frameLow.z) * dcellzi);
//     if (ptr_de != NULL && 
// 	(targetCellx >= numCell.x || 
// 	 targetCelly >= numCell.y || 
// 	 targetCellz >= numCell.z)){
//       *ptr_de = mdErrorOverFlowCellIdx;
//       return;
//     }
//     IndexType cellid = CudaDevice::D3toD1
// 	(numCell, targetCellx, targetCelly, targetCellz);
//     if (cellid != bid){
//       IndexType pid = atomicAdd (&numAtomInCell[cellid], 1);
//       if (pid >= blockDim.x){
// 	*ptr_de = mdErrorShortCellList;
// 	pid = 0;
//       }
//       IndexType targetIndex = pid + cellid * blockDim.x;
//       coord[targetIndex] = coord[ii];
//       coordNoi[targetIndex].x = coordNoi[ii].x;
//       coordNoi[targetIndex].y = coordNoi[ii].y;
//       coordNoi[targetIndex].z = coordNoi[ii].z;
//       velox[targetIndex] = velox[ii];
//       veloy[targetIndex] = veloy[ii];
//       veloz[targetIndex] = veloz[ii];
//       forcx[targetIndex] = forcx[ii];
//       forcy[targetIndex] = forcy[ii];
//       forcz[targetIndex] = forcz[ii];
//       globalIndex[targetIndex] = globalIndex[ii];
//       globalIndex[ii] = MaxIndexValue;
//       type[targetIndex] = type[ii];
//       mass[targetIndex] = mass[ii];
//       charge[targetIndex] = charge[ii];
//       IndexType my_numBond = (numBond[targetIndex] = numBond[ii]);
//       for (IndexType kk = 0; kk < my_numBond; ++kk){
// 	IndexType fromListIndex = Parallel::DeviceBondList_cudaDevice::indexConvert
// 	    (bondTopStride, ii, kk);
// 	IndexType toListIndex   = Parallel::DeviceBondList_cudaDevice::indexConvert
// 	    (bondTopStride, targetIndex, kk);
// 	bondNeighbor_globalIndex[toListIndex] = bondNeighbor_globalIndex[fromListIndex];
// 	bondIndex[toListIndex] = bondIndex[fromListIndex];
// 	bondNeighbor_localIndex[toListIndex] = bondNeighbor_localIndex[fromListIndex];
//       }
//     }
//   }
//   return;
// }


// __global__ void Parallel::CudaGlobal::
// rebuildCellList_step2 (IndexType * numAtomInCell,
// 		       CoordType  * coord,
// 		       CoordNoiType * coordNoi,
// 		       ScalorType * velox,
// 		       ScalorType * veloy,
// 		       ScalorType * veloz,
// 		       ScalorType * forcx,
// 		       ScalorType * forcy,
// 		       ScalorType * forcz,
// 		       IndexType  * globalIndex,
// 		       TypeType   * type,
// 		       ScalorType * mass,
// 		       ScalorType * charge,
// 		       IndexType * numBond,
// 		       IndexType * bondNeighbor_globalIndex,
// 		       IndexType * bondIndex,
// 		       IndexType * bondNeighbor_localIndex,
// 		       IndexType bondTopStride,
// 		       IndexType maxNumBond,
// 		       mdError_t * ptr_de)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;
//   IndexType ii = tid + bid * blockDim.x;
//   // IndexType k = NUintBit - 1;
  
//   extern __shared__ volatile IndexType sbuff[];
//   volatile IndexType * myIndex = (volatile IndexType * )sbuff;

//   if (tid >= numAtomInCell[bid] || globalIndex[ii] == MaxIndexValue){
//     myIndex[tid] = MaxIndexValue;
//   }
//   else {
//     myIndex[tid] = tid;
//   }
//   __syncthreads();
//   IndexType total = headSort (myIndex, &sbuff[blockDim.x]);
//   total = blockDim.x - total;

//   IndexType fromId;
//   CoordType bk_coord;
//   CoordNoiType bk_coordNoi;
//   ScalorType bk_velox, bk_veloy, bk_veloz;
//   ScalorType bk_forcx, bk_forcy, bk_forcz;
//   IndexType bk_globalIndex;
//   TypeType bk_type;
//   ScalorType bk_mass;
//   ScalorType bk_charge;
//   IndexType bk_numBond;
  
//   if (tid < total){
//     fromId = myIndex[tid] + bid * blockDim.x;
//     if (ii != fromId){
//       bk_coord = coord[fromId];
//       bk_coordNoi.x = coordNoi[fromId].x;
//       bk_coordNoi.y = coordNoi[fromId].y;
//       bk_coordNoi.z = coordNoi[fromId].z;
//       bk_velox = velox[fromId];
//       bk_veloy = veloy[fromId];
//       bk_veloz = veloz[fromId];
//       bk_forcx = forcx[fromId];
//       bk_forcy = forcy[fromId];
//       bk_forcz = forcz[fromId];
//       bk_globalIndex = globalIndex[fromId];
//       bk_type = type[fromId];
//       bk_mass = mass[fromId];
//       bk_charge = charge[fromId];
//       bk_numBond = numBond[fromId];
//     }
//   }
//   __syncthreads();

//   bool docopy = (tid < total && ii != fromId);
//   if (docopy){
//     coord[ii] = bk_coord;
//     coordNoi[ii].x = bk_coordNoi.x;
//     coordNoi[ii].y = bk_coordNoi.y;
//     coordNoi[ii].z = bk_coordNoi.z;
//     velox[ii] = bk_velox;
//     veloy[ii] = bk_veloy;
//     veloz[ii] = bk_veloz;
//     forcx[ii] = bk_forcx;
//     forcy[ii] = bk_forcy;
//     forcz[ii] = bk_forcz;
//     globalIndex[ii] = bk_globalIndex;
//     type[ii] = bk_type;
//     mass[ii] = bk_mass;
//     charge[ii] = bk_charge;
//     numBond[ii] = bk_numBond;
//   }

//   IndexType bk_bondNeighbor_globalIndex;
//   IndexType bk_bondIndex;
//   IndexType bk_bondNeighbor_localIndex;
  
//   for (IndexType kk = 0; kk < maxNumBond; ++kk){
//     __syncthreads();
//     if (docopy){
//       bk_bondNeighbor_globalIndex = bondNeighbor_globalIndex[fromId];
//       bk_bondIndex = bondIndex[fromId];
//       bk_bondNeighbor_localIndex = bondNeighbor_localIndex[fromId];
//     }
//     __syncthreads();
//     if (docopy){
//       bondNeighbor_globalIndex [ii] = bk_bondNeighbor_globalIndex;
//       bondIndex[ii] = bk_bondIndex;
//       bondNeighbor_localIndex[ii] = bk_bondNeighbor_localIndex;
//     }
//   }
	  
//   if (tid == 0){
//     numAtomInCell[bid] = total;
//   }
// }





void Parallel::DeviceCellListedMDData::
easyMallocCell (const IndexType & totalNumCell)
{
  if (totalNumCell == 0) return;
  // if (totalNumCell == numCell.x * numCell.y * numCell.z) return;
  // maxNumNeighborCell = maxNumNeighborCell_;
  clearCell ();
  cudaMalloc ((void**)&numAtomInCell, sizeof(IndexType) * totalNumCell);
  cudaMalloc ((void**)&forwardMap_step1,
	      sizeof(IndexType) * totalNumCell * Parallel::Interface::numThreadsInCell());
  cudaMalloc ((void**)&forwardMap_step2,
	      sizeof(IndexType) * totalNumCell * Parallel::Interface::numThreadsInCell());
  // cudaMalloc ((void**)&numNeighborCell, sizeof(IndexType) * totalNumCell);
  // cudaMalloc ((void**)&neighborCellIndex,
  // 	      sizeof(IndexType) * totalNumCell * maxNumNeighborCell);
  memSize = totalNumCell;
  checkCUDAError ("malloc Cell");
  malloced = true;
}

void Parallel::DeviceCellListedMDData::
clearCell()
{
  if (malloced){
    cudaFree (numAtomInCell);
    cudaFree (forwardMap_step2); 
    cudaFree (forwardMap_step1);
    memSize = 0;
    // cudaFree (numNeighborCell);
    // cudaFree (neighborCellIndex);
    malloced = false;
  }
}

Parallel::DeviceCellListedMDData::
DeviceCellListedMDData ()
{
  rlist = 0;
  devideLevel = 0;
  frameLow.x = frameLow.y = frameLow.z = 0;
  frameUp.x  = frameUp.y  = frameUp.z  = 0;
  numCell.x  = numCell.y  = numCell.z  = 0;
  memSize = 0;
  // maxNumNeighborCell = 0;
  malloced = false;
}

Parallel::DeviceCellListedMDData::
~DeviceCellListedMDData()
{
  clearCell();
}

Parallel::SubCellList::
SubCellList ()
{
}


void Parallel::SubCellList::
build ()
{
  Parallel::Interface::sort (this->begin(), this->end());
}

bool Parallel::SubCellList::
isBuilt ()
{
  return (Parallel::Interface::is_sorted (this->begin(), this->end()));
}

void Parallel::SubCellList::
add (const SubCellList & a)
{
  for (std::vector<IndexType>::const_iterator it = a.begin();
       it != a.end(); ++it){
    push_back (*it);
  }
  Parallel::Interface::unique (this->begin(), this->end());
  Parallel::Interface::sort   (this->begin(), this->end());
}

void Parallel::SubCellList::
sub (const SubCellList & a)
{
  std::vector<IndexType > result (this->size());
  std::vector<IndexType >::iterator newend =
      Parallel::Interface::set_difference (this->begin(), this->end(),
					   a.begin(), a.end(),
					   result.begin());
  std::vector<IndexType >::iterator newend2 =
      Parallel::Interface::copy (result.begin(), newend, this->begin());
  this->erase (newend2, this->end());
}



void Parallel::DeviceCellListedMDData::
buildSubList (const IndexType & xIdLo,
	      const IndexType & xIdUp,
	      const IndexType & yIdLo,
	      const IndexType & yIdUp,
	      const IndexType & zIdLo,
	      const IndexType & zIdUp,
	      SubCellList & subList) const
{
  if (xIdUp > numCell.x){
    throw MDExcptCellList ("x up index exceeds number of cells on x");
  }
  if (yIdUp > numCell.y){
    throw MDExcptCellList ("y up index exceeds number of cells on y");
  }
  if (zIdUp > numCell.z){
    throw MDExcptCellList ("z up index exceeds number of cells on z");
  }

  subList.clear();
  
  for (IndexType i = xIdLo; i < xIdUp; ++i){
    for (IndexType j = yIdLo; j < yIdUp; ++j){
      for (IndexType k = zIdLo; k < zIdUp; ++k){
	subList.push_back ( D3toD1 (i, j, k));
      }
    }
  }
}


Parallel::DeviceTransferPackage::
DeviceTransferPackage ()
    : numCell (0), memSize(0), hcellIndex(NULL), hcellStartIndex(NULL),
      myMask (MDDataItemMask_All)
{
}

void Parallel::DeviceTransferPackage::
clearMe ()
{
  if (memSize != 0){
    cudaFree (cellIndex);
    cudaFree (cellStartIndex);
    freeAPointer ((void**)&hcellIndex);
    freeAPointer ((void**)&hcellStartIndex);
    memSize = 0;
    numCell = 0;
  }
}

Parallel::DeviceTransferPackage::
~DeviceTransferPackage ()
{
  clearMe();
}

void Parallel::DeviceTransferPackage::
easyMallocMe (IndexType memSize_)
{
  if (memSize_ == 0) return;
  // if (memSize == memSize_) return;
  clearMe ();
  memSize = memSize_;
  size_t size = memSize * sizeof(IndexType);
  size_t size1 = (memSize+1) * sizeof(IndexType);
  cudaMalloc ((void**)&cellIndex, size);
  cudaMalloc ((void**)&cellStartIndex, size1);
  checkCUDAError ("DeviceTransferPackage::mallocMe failed malloc");
  hcellIndex = (IndexType *) malloc (size);
  if (hcellIndex == NULL){
    throw MDExcptFailedMallocOnHost ("DeviceTransferPackage::reinit",
				     "hcellIndex", size);
  }
  hcellStartIndex = (IndexType *) malloc (size1);
  if (hcellStartIndex == NULL){
    throw MDExcptFailedMallocOnHost ("DeviceTransferPackage::reinit",
				     "hcellStartIndex", size1);
  }
}

void Parallel::DeviceTransferPackage::
reinit (const SubCellList & subCellList)
{
  if (memSize < subCellList.size()){
    easyMallocMe (subCellList.size()*MemAllocExtension);
  }
  numCell = subCellList.size();
  for (IndexType i = 0; i < numCell; ++i){
    hcellIndex[i] = subCellList[i];
  }
  size_t size = memSize * sizeof(IndexType);
  cudaMemcpy (cellIndex, hcellIndex, size, cudaMemcpyHostToDevice);
  checkCUDAError ("DeviceTransferPackage::reinit memcpy");
}

void Parallel::DeviceTransferPackage::
pack (const DeviceCellListedMDData & ddata,
      const MDDataItemMask_t mask)
{
  if (numCell == 0) return;
  myMask = mask;
  
  IndexType totalNumCell = ddata.numCell.x * ddata.numCell.y * ddata.numCell.z;
  IndexType * numAtomInCell ;
  size_t size = totalNumCell * sizeof(IndexType);
  
  numAtomInCell = (IndexType *) malloc (size);
  if (numAtomInCell == NULL){
    throw MDExcptFailedMallocOnHost ("DeviceTransferPackage::reinit",
				     "numAtomInCell", size);
  }
  cudaMemcpy (numAtomInCell, ddata.numAtomInCell, size,
	      cudaMemcpyDeviceToHost);
  checkCUDAError ("DeviceTransferPackage::pack cpy numAtomInCell to host");
  
  hcellStartIndex[0] = 0;
  for (IndexType i = 1; i < numCell+1; ++i){
    hcellStartIndex[i] = hcellStartIndex[i-1] + numAtomInCell[hcellIndex[i-1]];
  }
  IndexType & expectedNumData (hcellStartIndex[numCell]);
  cudaMemcpy (cellStartIndex, hcellStartIndex, (numCell+1) * sizeof(IndexType),
	      cudaMemcpyHostToDevice);
  checkCUDAError ("DeviceTransferPackage::pack cpy cellStartIndex to device");
  
  free (numAtomInCell);

  IndexType expectedNumBond(0), expectedNumAngle(0), expectedNumDihedral(0);
  bool copyBond = (mask & MDDataItemMask_Bond);
  bool copyAngle = (mask & MDDataItemMask_Angle);
  bool copyDihedral = (mask & MDDataItemMask_Dihedral);  
  if (copyBond) expectedNumBond = ddata.getMaxNumBond();
  if (copyAngle) expectedNumAngle = ddata.getMaxNumAngle();
  if (copyDihedral) expectedNumDihedral = ddata.getMaxNumDihedral();  

  DeviceMDData::setGlobalBox (ddata.getGlobalBox());
  if (expectedNumData > DeviceMDData::memSize() ){
    // printf ("# DeviceTransferPackage::pack, realloc\n");
    DeviceMDData::easyMalloc(
	expectedNumData * MemAllocExtension,
	expectedNumBond, expectedNumAngle, expectedNumDihedral);
  }
  else if ((copyBond && (getMaxNumBond() != ddata.getMaxNumBond())) ||
	   (copyAngle && (getMaxNumAngle() != ddata.getMaxNumAngle())) ||
	   (copyDihedral && (getMaxNumDihedral() != ddata.getMaxNumDihedral())) ){
    // printf ("# DeviceTransferPackage::pack, realloc\n");
    DeviceMDData::easyMalloc(
	DeviceMDData::memSize(),
	expectedNumBond, expectedNumAngle, expectedNumDihedral);
  }
  DeviceMDData::numData() = expectedNumData;

  checkCUDAError ("DeviceTransferPackage::pack, packDeviceMDData, before");
  Parallel::CudaGlobal::packDeviceMDData
      <<<numCell, Parallel::Interface::numThreadsInCell()>>> (
	  cellIndex,
	  ddata.dptr_numAtomInCell(),
	  cellStartIndex,
	  mask,
	  ddata.dptr_coordinate(),
	  ddata.dptr_coordinateNoi(),
	  ddata.dptr_velocityX(),
	  ddata.dptr_velocityY(),
	  ddata.dptr_velocityZ(),
	  ddata.dptr_forceX(),
	  ddata.dptr_forceY(),
	  ddata.dptr_forceZ(),
	  ddata.dptr_globalIndex(),
	  ddata.dptr_type(),
	  ddata.dptr_mass(),
	  ddata.dptr_charge(),
	  this->dptr_coordinate(),
	  this->dptr_coordinateNoi(),
	  this->dptr_velocityX(),
	  this->dptr_velocityY(),
	  this->dptr_velocityZ(),
	  this->dptr_forceX(),
	  this->dptr_forceY(),
	  this->dptr_forceZ(),
	  this->dptr_globalIndex(),
	  this->dptr_type(),
	  this->dptr_mass(),
	  this->dptr_charge());
  Parallel::CudaGlobal::packDeviceMDData_bondTop
      <<<numCell, Parallel::Interface::numThreadsInCell()>>> (
	  cellIndex,
	  ddata.dptr_numAtomInCell(),
	  cellStartIndex,
	  mask,
	  ddata.dptr_numBond(),
	  ddata.dptr_bondIndex(),
	  ddata.dptr_bondNeighbor_globalIndex(),
	  ddata.dptr_numAngle(),
	  ddata.dptr_angleIndex(),
	  ddata.dptr_anglePosi(),
	  ddata.dptr_angleNeighbor_globalIndex(),
	  ddata.dptr_numDihedral(),
	  ddata.dptr_dihedralIndex(),
	  ddata.dptr_dihedralPosi(),
	  ddata.dptr_dihedralNeighbor_globalIndex(),
	  ddata.bondTopStride(),
	  ddata.getMaxNumBond(),
	  ddata.getMaxNumAngle(),
	  ddata.getMaxNumDihedral(),
	  this->dptr_numBond(),
	  this->dptr_bondIndex(),
	  this->dptr_bondNeighbor_globalIndex(),
	  this->dptr_numAngle(),
	  this->dptr_angleIndex(),
	  this->dptr_anglePosi(),
	  this->dptr_angleNeighbor_globalIndex(),
	  this->dptr_numDihedral(),
	  this->dptr_dihedralIndex(),
	  this->dptr_dihedralPosi(),
	  this->dptr_dihedralNeighbor_globalIndex(),
	  this->bondTopStride()
	  );
  checkCUDAError ("DeviceTransferPackage::pack, packDeviceMDData");
}


__global__ void Parallel::CudaGlobal::
packDeviceMDData (const IndexType * cellIndex,
		  const IndexType * numAtomInCell,
		  const IndexType * cellStartIndex,
		  const MDDataItemMask_t mask,
		  const CoordType  * source_coord,
		  const CoordNoiType * source_coordNoi,
		  const ScalorType * source_velox,
		  const ScalorType * source_veloy,
		  const ScalorType * source_veloz,
		  const ScalorType * source_forcx,
		  const ScalorType * source_forcy,
		  const ScalorType * source_forcz,
		  const IndexType  * source_globalIndex,
		  const TypeType   * source_type,
		  const ScalorType * source_mass,
		  const ScalorType * source_charge,
		  CoordType  * coord,
		  CoordNoiType * coordNoi,
		  ScalorType * velox,
		  ScalorType * veloy,
		  ScalorType * veloz,
		  ScalorType * forcx,
		  ScalorType * forcy,
		  ScalorType * forcz,
		  IndexType  * globalIndex,
		  TypeType   * type,
		  ScalorType * mass,
		  ScalorType * charge)		  
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  IndexType cellIdx = cellIndex[bid];
  IndexType fromid = tid + cellIdx * blockDim.x;
  IndexType toid = tid + cellStartIndex[bid];
  
  if (tid < numAtomInCell[cellIdx]){
    if (mask & MDDataItemMask_Coordinate){
      coord[toid] = source_coord[fromid];
    }
    if (mask & MDDataItemMask_CoordinateNoi){
      coordNoi[toid].x = source_coordNoi[fromid].x;
    }
    if (mask & MDDataItemMask_Velocity){
      velox[toid] = source_velox[fromid];
      veloy[toid] = source_veloy[fromid];
      veloz[toid] = source_veloz[fromid];
    }
    if (mask & MDDataItemMask_Force){
      forcx[toid] = source_forcx[fromid];
      forcy[toid] = source_forcy[fromid];
      forcz[toid] = source_forcz[fromid];
    }
    if (mask & MDDataItemMask_GlobalIndex){
      globalIndex[toid] = source_globalIndex[fromid];
    }
    if (mask & MDDataItemMask_Type){
      type[toid] = source_type[fromid];
    }
    if (mask & MDDataItemMask_Mass){
      mass[toid] = source_mass[fromid];
    }
    if (mask & MDDataItemMask_Charge){
      charge[toid] = source_charge[fromid];
    }
  }
}


__global__ void Parallel::CudaGlobal::
packDeviceMDData_bondTop (const IndexType * cellIndex,
			  const IndexType * numAtomInCell,
			  const IndexType * cellStartIndex,
			  const MDDataItemMask_t mask,
			  const IndexType * source_numBond,
			  const IndexType * source_bondIndex,
			  const IndexType * source_bondNeighbor_globalIndex,
			  const IndexType * source_numAngle,
			  const IndexType * source_angleIndex,
			  const IndexType * source_anglePosi ,
			  const IndexType * source_angleNeighbor_globalIndex,
			  const IndexType * source_numDihedral,
			  const IndexType * source_dihedralIndex,
			  const IndexType * source_dihedralPosi ,
			  const IndexType * source_dihedralNeighbor_globalIndex,
			  const IndexType   source_bondTopStride,
			  const IndexType   maxNumBond,
			  const IndexType   maxNumAngle,
			  const IndexType   maxNumDihedral,
			  IndexType * numBond,
			  IndexType * bondIndex,
			  IndexType * bondNeighbor_globalIndex,
			  IndexType * numAngle,
			  IndexType * angleIndex,
			  IndexType * anglePosi ,
			  IndexType * angleNeighbor_globalIndex,
			  IndexType * numDihedral,
			  IndexType * dihedralIndex,
			  IndexType * dihedralPosi ,
			  IndexType * dihedralNeighbor_globalIndex,
			  const IndexType bondTopStride)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  IndexType cellIdx = cellIndex[bid];
  IndexType fromid = tid + cellIdx * blockDim.x;
  IndexType toid = tid + cellStartIndex[bid];

  if (tid < numAtomInCell[cellIdx]){
    if ((mask & MDDataItemMask_Bond) && maxNumBond != 0) {
      IndexType my_fromid = fromid;
      IndexType my_toid   = toid;
      numBond[toid] = source_numBond[fromid];
      for (IndexType k = 0; k < maxNumBond; ++k){
	bondIndex[my_toid] = source_bondIndex[my_fromid];
	bondNeighbor_globalIndex[my_toid] = source_bondNeighbor_globalIndex[my_fromid];
	my_fromid += source_bondTopStride;
	my_toid   += bondTopStride;
      }
    }
    if ((mask & MDDataItemMask_Angle) && maxNumAngle != 0) {
      IndexType my_fromid = fromid;
      IndexType my_toid   = toid;
      numAngle[toid] = source_numAngle[fromid];
      for (IndexType k = 0; k < maxNumAngle; ++k){
	angleIndex[my_toid] = source_angleIndex[my_fromid];
	anglePosi [my_toid] = source_anglePosi [my_fromid];
	my_fromid += source_bondTopStride;
	my_toid   += bondTopStride;
      }
      my_fromid = fromid;
      my_toid   = toid;
      for (IndexType k = 0; k < maxNumAngle * 2; ++k){
	angleNeighbor_globalIndex[my_toid] = source_angleNeighbor_globalIndex[my_fromid];
	my_fromid += source_bondTopStride;
	my_toid   += bondTopStride;
      }
    }
    if ((mask & MDDataItemMask_Dihedral) && maxNumDihedral != 0) {
      IndexType my_fromid = fromid;
      IndexType my_toid   = toid;
      numDihedral[toid] = source_numDihedral[fromid];
      for (IndexType k = 0; k < maxNumDihedral; ++k){
	dihedralIndex[my_toid] = source_dihedralIndex[my_fromid];
	dihedralPosi [my_toid] = source_dihedralPosi [my_fromid];
	my_fromid += source_bondTopStride;
	my_toid   += bondTopStride;
      }
      my_fromid = fromid;
      my_toid   = toid;
      for (IndexType k = 0; k < maxNumDihedral * 3; ++k){
	dihedralNeighbor_globalIndex[my_toid] = source_dihedralNeighbor_globalIndex[my_fromid];
	my_fromid += source_bondTopStride;
	my_toid   += bondTopStride;
      }
    }
  }
}

	


void Parallel::DeviceTransferPackage::
copyToHost (HostTransferPackage & hpkg) const
{
  HostMDData & hdata(hpkg);
  const DeviceMDData & ddata(*this);  
  ddata.copyToHost (hdata, myMask);
  // alloc memory
  if (hpkg.getMemSize() < numCell){
    hpkg.easyMallocMe (numCell * MemAllocExtension);
  }
  hpkg.getTotalNumCell() = numCell;
  hpkg.getMask() = myMask;
  for (IndexType i = 0; i < numCell; ++i){
    hpkg.getCellIndex()[i] = hcellIndex[i];
    hpkg.getCellStartIndex()[i] = hcellStartIndex[i];
  }
  hpkg.getCellStartIndex()[numCell] = hcellStartIndex[numCell];  
}

void Parallel::DeviceTransferPackage::
copyFromHost (const HostTransferPackage & hpkg)
{
  const HostMDData & hdata(hpkg);
  DeviceMDData & ddata(*this);
  myMask = hpkg.getMask();
  ddata.copyFromHost (hdata, myMask);

  if (memSize < hpkg.getTotalNumCell()){
    easyMallocMe (hpkg.getTotalNumCell() * MemAllocExtension);
  }
  numCell = hpkg.getTotalNumCell();
  for (IndexType i = 0; i < numCell; ++i){
    hcellIndex[i] = hpkg.getCellIndex()[i];
    hcellStartIndex[i] = hpkg.getCellStartIndex()[i];
  }
  hcellStartIndex[numCell] = hpkg.getCellStartIndex()[numCell];

  cudaMemcpy (cellIndex, hcellIndex, sizeof(IndexType)*numCell,
	      cudaMemcpyHostToDevice);
  cudaMemcpy (cellStartIndex, hcellStartIndex, sizeof(IndexType)*(numCell+1),
	      cudaMemcpyHostToDevice);
  checkCUDAError ("DeviceTransferPackage::copyFromHost cpy from host");
}


void Parallel::DeviceTransferPackage::
unpack_replace (DeviceCellListedMDData & ddata) const
{
  Parallel::CudaGlobal::unpackDeviceMDData_replace
      <<<numCell, Parallel::Interface::numThreadsInCell()>>> (
	  cellIndex,
	  cellStartIndex,
	  myMask,
	  this->dptr_coordinate(),
	  this->dptr_coordinateNoi(),
	  this->dptr_velocityX(),
	  this->dptr_velocityY(),
	  this->dptr_velocityZ(),
	  this->dptr_forceX(),
	  this->dptr_forceY(),
	  this->dptr_forceZ(),
	  this->dptr_globalIndex(),
	  this->dptr_type(),
	  this->dptr_mass(),
	  this->dptr_charge(),
	  ddata.numAtomInCell,
	  ddata.dptr_coordinate(),
	  ddata.dptr_coordinateNoi(),
	  ddata.dptr_velocityX(),
	  ddata.dptr_velocityY(),
	  ddata.dptr_velocityZ(),
	  ddata.dptr_forceX(),
	  ddata.dptr_forceY(),
	  ddata.dptr_forceZ(),
	  ddata.dptr_globalIndex(),
	  ddata.dptr_type(),
	  ddata.dptr_mass(),
	  ddata.dptr_charge());
  Parallel::CudaGlobal::unpackDeviceMDData_bondTop_replace
      <<<numCell, Parallel::Interface::numThreadsInCell()>>>(
	  cellIndex,
	  cellStartIndex,
	  myMask,
	  this->dptr_numBond(),
	  this->dptr_bondIndex(),
	  this->dptr_bondNeighbor_globalIndex(),
	  this->dptr_numAngle(),
	  this->dptr_angleIndex(),
	  this->dptr_anglePosi(),
	  this->dptr_angleNeighbor_globalIndex(),
	  this->dptr_numDihedral(),
	  this->dptr_dihedralIndex(),
	  this->dptr_dihedralPosi(),
	  this->dptr_dihedralNeighbor_globalIndex(),
	  this->bondTopStride(),
	  this->getMaxNumBond(),
	  this->getMaxNumAngle(),
	  this->getMaxNumDihedral(),
	  ddata.dptr_numBond(),
	  ddata.dptr_bondIndex(),
	  ddata.dptr_bondNeighbor_globalIndex(),
	  ddata.dptr_numAngle(),
	  ddata.dptr_angleIndex(),
	  ddata.dptr_anglePosi(),
	  ddata.dptr_angleNeighbor_globalIndex(),
	  ddata.dptr_numDihedral(),
	  ddata.dptr_dihedralIndex(),
	  ddata.dptr_dihedralPosi(),
	  ddata.dptr_dihedralNeighbor_globalIndex(),
	  ddata.bondTopStride());
  checkCUDAError ("DeviceTransferPackage::unpack_replace");
}



void Parallel::DeviceTransferPackage::
unpack_add (DeviceCellListedMDData & ddata) const
{
  IndexType totalNumInCell = ddata.getNumCell().x * ddata.getNumCell().y * ddata.getNumCell().z;
  IndexType * oldNumAtomInCell;
  size_t size = totalNumInCell * sizeof(IndexType);
  cudaMalloc ((void**)&oldNumAtomInCell, size);
  checkCUDAError ("DeviceTransferPackage::unpack_add, malloc oldNumAtomInCell");
  cudaMemcpy (oldNumAtomInCell, ddata.dptr_numAtomInCell(), size,
	      cudaMemcpyDeviceToDevice);
  checkCUDAError ("DeviceTransferPackage::unpack_add, copy oldNumAtomInCell");
  
  Parallel::CudaGlobal::unpackDeviceMDData_add
      <<<numCell, Parallel::Interface::numThreadsInCell()>>> (
	  cellIndex,
	  cellStartIndex,
	  myMask,
	  this->dptr_coordinate(),
	  this->dptr_coordinateNoi(),
	  this->dptr_velocityX(),
	  this->dptr_velocityY(),
	  this->dptr_velocityZ(),
	  this->dptr_forceX(),
	  this->dptr_forceY(),
	  this->dptr_forceZ(),
	  this->dptr_globalIndex(),
	  this->dptr_type(),
	  this->dptr_mass(),
	  this->dptr_charge(),
	  ddata.dptr_numAtomInCell(),
	  ddata.dptr_coordinate(),
	  ddata.dptr_coordinateNoi(),
	  ddata.dptr_velocityX(),
	  ddata.dptr_velocityY(),
	  ddata.dptr_velocityZ(),
	  ddata.dptr_forceX(),
	  ddata.dptr_forceY(),
	  ddata.dptr_forceZ(),
	  ddata.dptr_globalIndex(),
	  ddata.dptr_type(),
	  ddata.dptr_mass(),
	  ddata.dptr_charge(),
	  err.ptr_de);
  Parallel::CudaGlobal::unpackDeviceMDData_bondTop_add
      <<<numCell, Parallel::Interface::numThreadsInCell()>>>(
	  cellIndex,
	  cellStartIndex,
	  oldNumAtomInCell,
	  myMask,
	  this->dptr_numBond(),
	  this->dptr_bondIndex(),
	  this->dptr_bondNeighbor_globalIndex(),
	  this->dptr_numAngle(),
	  this->dptr_angleIndex(),
	  this->dptr_anglePosi(),
	  this->dptr_angleNeighbor_globalIndex(),
	  this->dptr_numDihedral(),
	  this->dptr_dihedralIndex(),
	  this->dptr_dihedralPosi(),
	  this->dptr_dihedralNeighbor_globalIndex(),
	  this->bondTopStride(),
	  this->getMaxNumBond(),
	  this->getMaxNumAngle(),
	  this->getMaxNumDihedral(),
	  ddata.dptr_numBond(),
	  ddata.dptr_bondIndex(),
	  ddata.dptr_bondNeighbor_globalIndex(),
	  ddata.dptr_numAngle(),
	  ddata.dptr_angleIndex(),
	  ddata.dptr_anglePosi(),
	  ddata.dptr_angleNeighbor_globalIndex(),
	  ddata.dptr_numDihedral(),
	  ddata.dptr_dihedralIndex(),
	  ddata.dptr_dihedralPosi(),
	  ddata.dptr_dihedralNeighbor_globalIndex(),
	  ddata.bondTopStride(),
	  err.ptr_de);	  
  checkCUDAError ("DeviceTransferPackage::unpack_add");
  err.updateHost ();
  cudaFree (oldNumAtomInCell);
  err.check ("DeviceTransferPackage::unpack_add");
}



__global__ void Parallel::CudaGlobal::
unpackDeviceMDData_replace (const IndexType * cellIndex,
			  const IndexType * cellStartIndex,
			  const MDDataItemMask_t mask,
			  const CoordType  * source_coord,
			  const CoordNoiType * source_coordNoi,
			  const ScalorType * source_velox,
			  const ScalorType * source_veloy,
			  const ScalorType * source_veloz,
			  const ScalorType * source_forcx,
			  const ScalorType * source_forcy,
			  const ScalorType * source_forcz,
			  const IndexType  * source_globalIndex,
			  const TypeType   * source_type,
			  const ScalorType * source_mass,
			  const ScalorType * source_charge,
			  IndexType * numAtomInCell,
			  CoordType  * coord,
			  CoordNoiType * coordNoi,
			  ScalorType * velox,
			  ScalorType * veloy,
			  ScalorType * veloz,
			  ScalorType * forcx,
			  ScalorType * forcy,
			  ScalorType * forcz,
			  IndexType  * globalIndex,
			  TypeType   * type,
			  ScalorType * mass,
			  ScalorType * charge)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  IndexType cellIdx = cellIndex[bid];
  IndexType startIdx = cellStartIndex[bid];
  IndexType numAtomInThisCell = cellStartIndex[bid+1] - startIdx;
  IndexType toid   = tid + cellIdx * blockDim.x;
  IndexType fromid = tid + startIdx;
  if (tid == blockDim.x-1){
    numAtomInCell[cellIdx] = numAtomInThisCell;
  }
  
  if (tid < numAtomInThisCell){
    if (mask & MDDataItemMask_Coordinate){
      coord[toid] = source_coord[fromid];
    }
    if (mask & MDDataItemMask_CoordinateNoi){
      coordNoi[toid].x = source_coordNoi[fromid].x;
      coordNoi[toid].y = source_coordNoi[fromid].y;
      coordNoi[toid].z = source_coordNoi[fromid].z;
    }
    if (mask & MDDataItemMask_Velocity){
      velox[toid] = source_velox[fromid];
      veloy[toid] = source_veloy[fromid];
      veloz[toid] = source_veloz[fromid];
    }
    if (mask & MDDataItemMask_Force){
      forcx[toid] = source_forcx[fromid];
      forcy[toid] = source_forcy[fromid];
      forcz[toid] = source_forcz[fromid];
    }
    if (mask & MDDataItemMask_GlobalIndex){
      globalIndex[toid] = source_globalIndex[fromid];
    }
    if (mask & MDDataItemMask_Type){
      type[toid] = source_type[fromid];
    }
    if (mask & MDDataItemMask_Mass){
      mass[toid] = source_mass[fromid];
    }
    if (mask & MDDataItemMask_Charge){
      charge[toid] = source_charge[fromid];
    }
  }
}

__global__ void Parallel::CudaGlobal::
unpackDeviceMDData_bondTop_replace (const IndexType * cellIndex,
				    const IndexType * cellStartIndex,
				    const MDDataItemMask_t mask,
				    const IndexType * source_numBond,
				    const IndexType * source_bondIndex,
				    const IndexType * source_bondNeighbor_globalIndex,
				    const IndexType * source_numAngle,
				    const IndexType * source_angleIndex,
				    const IndexType * source_anglePosi ,
				    const IndexType * source_angleNeighbor_globalIndex,
				    const IndexType * source_numDihedral,
				    const IndexType * source_dihedralIndex,
				    const IndexType * source_dihedralPosi ,
				    const IndexType * source_dihedralNeighbor_globalIndex,
				    const IndexType source_bondTopStride,
				    const IndexType maxNumBond,
				    const IndexType maxNumAngle,
				    const IndexType maxNumDihedral,
				    IndexType * numBond,
				    IndexType * bondIndex,
				    IndexType * bondNeighbor_globalIndex,
				    IndexType * numAngle,
				    IndexType * angleIndex,
				    IndexType * anglePosi ,
				    IndexType * angleNeighbor_globalIndex,
				    IndexType * numDihedral,
				    IndexType * dihedralIndex,
				    IndexType * dihedralPosi ,
				    IndexType * dihedralNeighbor_globalIndex,
				    const IndexType bondTopStride)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  IndexType cellIdx = cellIndex[bid];
  IndexType startIdx = cellStartIndex[bid];
  IndexType numAtomInThisCell = cellStartIndex[bid+1] - startIdx;
  IndexType toid   = tid + cellIdx * blockDim.x;
  IndexType fromid = tid + startIdx;
  
  if (tid < numAtomInThisCell){
    if ((mask & MDDataItemMask_Bond) && maxNumBond != 0) {
      IndexType my_fromid = fromid;
      IndexType my_toid   = toid;
      numBond[toid] = source_numBond[fromid];
      for (IndexType k = 0; k < maxNumBond; ++k){
	bondIndex[my_toid] = source_bondIndex[my_fromid];
	bondNeighbor_globalIndex[my_toid] = source_bondNeighbor_globalIndex[my_fromid];
	my_fromid += source_bondTopStride;
	my_toid   += bondTopStride;
      }
    }
    if ((mask & MDDataItemMask_Angle) && maxNumAngle != 0) {
      IndexType my_fromid = fromid;
      IndexType my_toid   = toid;
      numAngle[toid] = source_numAngle[fromid];
      for (IndexType k = 0; k < maxNumAngle; ++k){
	angleIndex[my_toid] = source_angleIndex[my_fromid];
	anglePosi [my_toid] = source_anglePosi [my_fromid];
	my_fromid += source_bondTopStride;
	my_toid   += bondTopStride;
      }
      my_fromid = fromid;
      my_toid   = toid;
      for (IndexType k = 0; k < maxNumAngle * 2; ++k){
	angleNeighbor_globalIndex[my_toid] = source_angleNeighbor_globalIndex[my_fromid];
	my_fromid += source_bondTopStride;
	my_toid   += bondTopStride;
      }
    }
    if ((mask & MDDataItemMask_Dihedral) && maxNumDihedral != 0) {
      IndexType my_fromid = fromid;
      IndexType my_toid   = toid;
      numDihedral[toid] = source_numDihedral[fromid];
      for (IndexType k = 0; k < maxNumDihedral; ++k){
	dihedralIndex[my_toid] = source_dihedralIndex[my_fromid];
	dihedralPosi [my_toid] = source_dihedralPosi [my_fromid];
	my_fromid += source_bondTopStride;
	my_toid   += bondTopStride;
      }
      my_fromid = fromid;
      my_toid   = toid;
      for (IndexType k = 0; k < maxNumDihedral * 3; ++k){
	dihedralNeighbor_globalIndex[my_toid] = source_dihedralNeighbor_globalIndex[my_fromid];
	my_fromid += source_bondTopStride;
	my_toid   += bondTopStride;
      }
    }
  }
}



__global__ void Parallel::CudaGlobal::
unpackDeviceMDData_add (const IndexType * cellIndex,
			const IndexType * cellStartIndex,
			const MDDataItemMask_t mask,
			const CoordType  * source_coord,
			const CoordNoiType * source_coordNoi,
			const ScalorType * source_velox,
			const ScalorType * source_veloy,
			const ScalorType * source_veloz,
			const ScalorType * source_forcx,
			const ScalorType * source_forcy,
			const ScalorType * source_forcz,
			const IndexType  * source_globalIndex,
			const TypeType   * source_type,
			const ScalorType * source_mass,
			const ScalorType * source_charge,
			IndexType * numAtomInCell,
			CoordType  * coord,
			CoordNoiType * coordNoi,
			ScalorType * velox,
			ScalorType * veloy,
			ScalorType * veloz,
			ScalorType * forcx,
			ScalorType * forcy,
			ScalorType * forcz,
			IndexType  * globalIndex,
			TypeType   * type,
			ScalorType * mass,
			ScalorType * charge,
			mdError_t * ptr_de)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  __shared__ volatile IndexType tmpbuff[2];
  if (tid < 2){
    tmpbuff[tid] = cellStartIndex[bid+tid];
  }
  __syncthreads();
  IndexType startIdx = tmpbuff[0];
  IndexType numAdded = tmpbuff[1] - startIdx;
  // IndexType startIdx = cellStartIndex[bid];
  // IndexType numAdded = cellStartIndex[bid+1] - startIdx;
  if (numAdded == 0) return;
  
  IndexType cellIdx = cellIndex[bid];
  IndexType alreadyInCell = numAtomInCell[cellIdx];
  IndexType toid   = tid + cellIdx * blockDim.x + alreadyInCell;
  IndexType fromid = tid + startIdx;
  // __shared__ IndexType numAtomInThisCell;
  
  __syncthreads();
  __shared__ volatile bool failed ;
  if (tid == 0){
    if (((numAtomInCell[cellIdx] += numAdded)) > blockDim.x &&
	ptr_de != NULL){
      *ptr_de = mdErrorShortCellList;
      failed = true;
    }
    else {
      failed = false;
    }
  }
  __syncthreads();
  if (failed) return;
  
  if (tid < numAdded){
    if (mask & MDDataItemMask_Coordinate){
      coord[toid] = source_coord[fromid];
    }
    if (mask & MDDataItemMask_CoordinateNoi){
      coordNoi[toid].x = source_coordNoi[fromid].x;
      coordNoi[toid].y = source_coordNoi[fromid].y;
      coordNoi[toid].z = source_coordNoi[fromid].z;
    }
    if (mask & MDDataItemMask_Velocity){
      velox[toid] = source_velox[fromid];
      veloy[toid] = source_veloy[fromid];
      veloz[toid] = source_veloz[fromid];
    }
    if (mask & MDDataItemMask_Force){
      forcx[toid] = source_forcx[fromid];
      forcy[toid] = source_forcy[fromid];
      forcz[toid] = source_forcz[fromid];
    }
    if (mask & MDDataItemMask_GlobalIndex){
      globalIndex[toid] = source_globalIndex[fromid];
    }
    if (mask & MDDataItemMask_Type){
      type[toid] = source_type[fromid];
    }
    if (mask & MDDataItemMask_Mass){
      mass[toid] = source_mass[fromid];
    }
    if (mask & MDDataItemMask_Charge){
      charge[toid] = source_charge[fromid];
    }
  }
}

__global__ void Parallel::CudaGlobal::
unpackDeviceMDData_bondTop_add (const IndexType * cellIndex,
				const IndexType * cellStartIndex,
				const IndexType * oldNumAtomInCell,
				const MDDataItemMask_t mask,
				const IndexType * source_numBond,
				const IndexType * source_bondIndex,
				const IndexType * source_bondNeighbor_globalIndex,
				const IndexType * source_numAngle,
				const IndexType * source_angleIndex,
				const IndexType * source_anglePosi ,
				const IndexType * source_angleNeighbor_globalIndex,
				const IndexType * source_numDihedral,
				const IndexType * source_dihedralIndex,
				const IndexType * source_dihedralPosi ,
				const IndexType * source_dihedralNeighbor_globalIndex,
				const IndexType source_bondTopStride,
				const IndexType maxNumBond,
				const IndexType maxNumAngle,
				const IndexType maxNumDihedral,
				IndexType * numBond,
				IndexType * bondIndex,
				IndexType * bondNeighbor_globalIndex,
				IndexType * numAngle,
				IndexType * angleIndex,
				IndexType * anglePosi ,
				IndexType * angleNeighbor_globalIndex,
				IndexType * numDihedral,
				IndexType * dihedralIndex,
				IndexType * dihedralPosi ,
				IndexType * dihedralNeighbor_globalIndex,
				const IndexType bondTopStride,
				mdError_t * ptr_de)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  __shared__ volatile IndexType tmpbuff[2];
  if (tid < 2){
    tmpbuff[tid] = cellStartIndex[bid+tid];
  }
  __syncthreads();
  IndexType startIdx = tmpbuff[0];
  IndexType numAdded = tmpbuff[1] - startIdx;
  if (numAdded == 0) return;
  
  IndexType cellIdx = cellIndex[bid];
  IndexType toid   = tid + cellIdx * blockDim.x + oldNumAtomInCell[cellIdx];
  IndexType fromid = tid + startIdx;
  
  if (tid < numAdded){
    if ((mask & MDDataItemMask_Bond) && maxNumBond != 0) {
      IndexType my_fromid = fromid;
      IndexType my_toid   = toid;
      numBond[toid] = source_numBond[fromid];
      for (IndexType k = 0; k < maxNumBond; ++k){
	bondIndex[my_toid] = source_bondIndex[my_fromid];
	bondNeighbor_globalIndex[my_toid] = source_bondNeighbor_globalIndex[my_fromid];
	my_fromid += source_bondTopStride;
	my_toid   += bondTopStride;
      }
    }
    if ((mask & MDDataItemMask_Angle) && maxNumAngle != 0) {
      IndexType my_fromid = fromid;
      IndexType my_toid   = toid;
      numAngle[toid] = source_numAngle[fromid];
      for (IndexType k = 0; k < maxNumAngle; ++k){
	angleIndex[my_toid] = source_angleIndex[my_fromid];
	anglePosi [my_toid] = source_anglePosi [my_fromid];
	my_fromid += source_bondTopStride;
	my_toid   += bondTopStride;
      }
      my_fromid = fromid;
      my_toid   = toid;
      for (IndexType k = 0; k < maxNumAngle * 2; ++k){
	angleNeighbor_globalIndex[my_toid] = source_angleNeighbor_globalIndex[my_fromid];
	my_fromid += source_bondTopStride;
	my_toid   += bondTopStride;
      }
    }
    if ((mask & MDDataItemMask_Dihedral) && maxNumDihedral != 0) {
      IndexType my_fromid = fromid;
      IndexType my_toid   = toid;
      numDihedral[toid] = source_numDihedral[fromid];
      for (IndexType k = 0; k < maxNumDihedral; ++k){
	dihedralIndex[my_toid] = source_dihedralIndex[my_fromid];
	dihedralPosi [my_toid] = source_dihedralPosi [my_fromid];
	my_fromid += source_bondTopStride;
	my_toid   += bondTopStride;
      }
      my_fromid = fromid;
      my_toid   = toid;
      for (IndexType k = 0; k < maxNumDihedral * 3; ++k){
	dihedralNeighbor_globalIndex[my_toid] = source_dihedralNeighbor_globalIndex[my_fromid];
	my_fromid += source_bondTopStride;
	my_toid   += bondTopStride;
      }
    }
  }
}



void Parallel::DeviceCellListedMDData::
copyToHost (HostCellListedMDData & hdata,
	    const MDDataItemMask_t mask) const
{
  const DeviceMDData & ddata (*this);
  ddata.copyToHost (hdata, mask);

  hdata.rlist = rlist;
  hdata.devideLevel = devideLevel;
  hdata.numCell.x = numCell.x;
  hdata.numCell.y = numCell.y;
  hdata.numCell.z = numCell.z;
  hdata.frameUp.x = frameUp.x;
  hdata.frameUp.y = frameUp.y;
  hdata.frameUp.z = frameUp.z;
  hdata.frameLow.x = frameLow.x;
  hdata.frameLow.y = frameLow.y;
  hdata.frameLow.z = frameLow.z;

  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  size_t size = totalNumCell * sizeof (IndexType);

  if (hdata.HostCellListedMDData::memSize < totalNumCell){
    hdata.easyReallocCell (totalNumCell * MemAllocExtension);
  }
  cudaMemcpy (hdata.numAtomInCell, numAtomInCell, size,
	      cudaMemcpyDeviceToHost);
  checkCUDAError ("DeviceCellListedMDData::copyToHost copy numAtomInCell");
}

  
void Parallel::DeviceCellListedMDData::
copyFromHost (const HostCellListedMDData & hdata,
	      const MDDataItemMask_t mask)
{
  DeviceMDData & ddata(*this);
  ddata.copyFromHost (hdata, mask);

  rlist = hdata.rlist;
  devideLevel = hdata.devideLevel;
  numCell.x = hdata.numCell.x;
  numCell.y = hdata.numCell.y;
  numCell.z = hdata.numCell.z;
  frameUp.x = hdata.frameUp.x;
  frameUp.y = hdata.frameUp.y;
  frameUp.z = hdata.frameUp.z;
  frameLow.x = hdata.frameLow.x;
  frameLow.y = hdata.frameLow.y;
  frameLow.z = hdata.frameLow.z;
  
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  size_t size = totalNumCell * sizeof (IndexType);

  if (memSize < totalNumCell){
    easyMallocCell (totalNumCell * MemAllocExtension);
  }
  cudaMemcpy (numAtomInCell, hdata.numAtomInCell, size,
	      cudaMemcpyHostToDevice);
  checkCUDAError ("DeviceCellListedMDData::copyFromHost copy numAtomInCell");
}

void Parallel::DeviceCellListedMDData::
copyFromDevice (const DeviceCellListedMDData & ddata,
		const MDDataItemMask_t mask)
{  
  DeviceMDData & me(*this);
  me.copyFromDevice (ddata, mask);

  rlist = ddata.rlist;
  devideLevel = ddata.devideLevel;
  numCell.x = ddata.numCell.x;
  numCell.y = ddata.numCell.y;
  numCell.z = ddata.numCell.z;
  frameUp.x = ddata.frameUp.x;
  frameUp.y = ddata.frameUp.y;
  frameUp.z = ddata.frameUp.z;
  frameLow.x = ddata.frameLow.x;
  frameLow.y = ddata.frameLow.y;
  frameLow.z = ddata.frameLow.z;
  
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  size_t size = totalNumCell * sizeof (IndexType);

  if (memSize < totalNumCell){
    easyMallocCell (totalNumCell * MemAllocExtension);
  }

  cudaMemcpy (numAtomInCell, ddata.numAtomInCell, size,
	      cudaMemcpyDeviceToDevice);
  checkCUDAError ("DeviceCellListedMDData::copyFromDevice copy numAtomInCell");
}


void Parallel::DeviceCellListedMDData::
clearData (const SubCellList & subList)
{
  IndexType * tmpList;
  IndexType * deviceList;
  IndexType num = subList.size();
  size_t size = sizeof(IndexType) * num;
  
  tmpList = (IndexType * )malloc (size);
  if (tmpList == NULL){
    throw MDExcptFailedMallocOnHost ("DeviceCellListedMDData::clearData",
				     "tmpList", size);
  }
  cudaMalloc ((void**)&deviceList, size);
  checkCUDAError ("DeviceCellListedMDData::clearData, malloc deviceList");

  for (IndexType i = 0; i < num; ++i){
    tmpList[i] = subList[i];
  }
  cudaMemcpy (deviceList, tmpList, size, cudaMemcpyHostToDevice);
  freeAPointer ((void**)&tmpList);
  
  Parallel::CudaGlobal::clearCellListData
      <<<(num + DefaultNThreadPerBlock -1) / DefaultNThreadPerBlock,
      DefaultNThreadPerBlock>>> (
	  deviceList,
	  num,
	  numAtomInCell);
  checkCUDAError ("DeviceCellListedMDData::clearData, clearCellListData");
  cudaFree (deviceList);
  checkCUDAError ("DeviceCellListedMDData::clearData, free deviceList");  
}


void __global__
Parallel::CudaGlobal::
clearCellListData (const IndexType * deviceList,
		   IndexType num,
		   IndexType * numAtomInCell)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;

  if (ii < num){
    numAtomInCell[deviceList[ii]] = 0;
  }
}



__global__ void Parallel::CudaGlobal::
normalizeSystem_CellListed (RectangularBox box,
			    const IndexType * numAtomInCell,
			    CoordType * coord,
			    CoordNoiType * coordNoi)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
 
  if (tid < numAtomInCell[bid]) {
    RectangularBoxGeometry::moveParticleToBox_1image (
	box.size.x, &(coord[ii].x), &(coordNoi[ii].x));
    RectangularBoxGeometry::moveParticleToBox_1image (
	box.size.y, &(coord[ii].y), &(coordNoi[ii].y));
    RectangularBoxGeometry::moveParticleToBox_1image (
	box.size.z, &(coord[ii].z), &(coordNoi[ii].z));
  }
}


void Parallel::DeviceCellListedMDData::
applyPeriodicBondaryCondition ()
{
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell ();
  
  Parallel::CudaGlobal::normalizeSystem_CellListed
      <<<totalNumCell, numThreadsInCell>>> (
	  globalBox,
	  numAtomInCell,
	  coord,
	  coordNoi);
  checkCUDAError ("DeviceCellListedMDData::applyPeriodicBondaryCondition");
}

__global__ void Parallel::CudaGlobal::
normalizeSystemOnGhost_CellListed (RectangularBox box,
				   const IndexType * numAtomInCell,
				   const IntVectorType numCell,
				   const IndexType divideLevel,
				   const int nx,
				   const int ny,
				   const int nz,
				   const int rankx,
				   const int ranky,
				   const int rankz,
				   CoordType * coord,
				   CoordNoiType * coordNoi)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  IndexType this_numAtomInCell = numAtomInCell[bid];
  
  IndexType ix, iy, iz;
  Parallel::CudaDevice::D1toD3 (numCell, bid, ix, iy, iz);

  if (rankx == 0 && ix < divideLevel){
    if (tid < this_numAtomInCell){
      coord[ii].x += box.size.x;
      coordNoi[ii].x -= 1;
    }
  }
  if (rankx == nx - 1 && ix >= numCell.x - divideLevel){
    if (tid < this_numAtomInCell){
      coord[ii].x -= box.size.x;
      coordNoi[ii].x += 1;
    }
  }
  if (ranky == 0 && iy < divideLevel){
    if (tid < this_numAtomInCell){
      coord[ii].y += box.size.y;
      coordNoi[ii].y -= 1;
    }
  }
  if (ranky == ny - 1 && iy >= numCell.y - divideLevel){
    if (tid < this_numAtomInCell){
      coord[ii].y -= box.size.y;
      coordNoi[ii].y += 1;
    }
  }
  if (rankz == 0 && iz < divideLevel){
    if (tid < this_numAtomInCell){
      coord[ii].z += box.size.z;
      coordNoi[ii].z -= 1;
    }
  }
  if (rankz == nz - 1 && iz >= numCell.z - divideLevel){
    if (tid < this_numAtomInCell){
      coord[ii].z -= box.size.z;
      coordNoi[ii].z += 1;
    }
  }
  // if (tid < numAtomInCell[bid]) {
  //   RectangularBoxGeometry::moveParticleToBox_1image (
  // 	box.size.x, &(coord[ii].x), &(coordNoi[ii].x));
  //   RectangularBoxGeometry::moveParticleToBox_1image (
  // 	box.size.y, &(coord[ii].y), &(coordNoi[ii].y));
  //   RectangularBoxGeometry::moveParticleToBox_1image (
  // 	box.size.z, &(coord[ii].z), &(coordNoi[ii].z));
  // }
}



void Parallel::DeviceCellListedMDData::
applyPeriodicBondaryConditionOnGhostCells ()
{
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell ();
  int Nx, Ny, Nz;
  Parallel::Interface::numProcDim (Nx, Ny, Nz);
  int rankx, ranky, rankz;
  Parallel::Interface::rankToCartCoord (Parallel::Interface::myRank(),
					rankx, ranky, rankz);

  Parallel::CudaGlobal::normalizeSystemOnGhost_CellListed
      <<<totalNumCell, numThreadsInCell>>> (
	  globalBox,
	  numAtomInCell,
	  numCell,
	  devideLevel,
	  Nx, Ny, Nz,
	  rankx, ranky, rankz,
	  coord,
	  coordNoi);
  checkCUDAError ("DeviceCellListedMDData::applyPeriodicBondaryCondition");
}


void Parallel::DeviceCellRelation::
easyMalloc (const IndexType & totalNumCell_,
	    const IndexType & MaxNeiPerCell_)
{
  if (malloced){
    clear ();
  }
  totalNumCell = totalNumCell_;
  MaxNeiPerCell = MaxNeiPerCell_;
  
  cudaMalloc ((void**)&numNeighbor, totalNumCell * sizeof(IndexType));
  cudaMalloc ((void**)&neighborCellIndex,
	      totalNumCell * MaxNeiPerCell * sizeof (IndexType));
  cudaMalloc ((void**)&neighborShiftNoi,
	      totalNumCell * MaxNeiPerCell * sizeof (CoordNoiType));
  checkCUDAError ("DeviceCellRelation::easyMalloc");
  malloced = true;
}

void Parallel::DeviceCellRelation::
clear ()
{
  if (malloced){
    cudaFree (numNeighbor);
    cudaFree (neighborCellIndex);
    cudaFree (neighborShiftNoi);
    malloced = false;
  }
}

Parallel::DeviceCellRelation::
~DeviceCellRelation ()
{
  clear();
}

void Parallel::DeviceCellRelation::
rebuild (const DeviceCellListedMDData & list)
{
  ptr_list = &list;

  IntVectorType numCell = list.getNumCell ();
  totalNumCell = numCell.x * numCell.y * numCell.z;
  MaxNeiPerCell = 2 * list.getDevideLevel() + 1;
  MaxNeiPerCell = MaxNeiPerCell * MaxNeiPerCell * MaxNeiPerCell;
  int Nx, Ny, Nz;
  Parallel::Interface::numProcDim (Nx, Ny, Nz);
  
  easyMalloc (totalNumCell, MaxNeiPerCell);
  int rankx, ranky, rankz;
  Parallel::Interface::rankToCartCoord (Parallel::Interface::myRank(), rankx, ranky, rankz);
  
  Parallel::CudaGlobal::buildCellNeighborhood
      <<<toGridDim(totalNumCell), 1>>> (
	  numCell,
	  list.getDevideLevel(),
	  list.getRlist(),
	  list.getGlobalBoxSize(),
	  rankx, ranky, rankz,
	  Nx, Ny, Nz,
	  numNeighbor,
	  neighborCellIndex,
	  neighborShiftNoi,
	  MaxNeiPerCell);
  checkCUDAError ("DeviceCellRelation::build, buildCellNeighborhood");
}


__global__ void Parallel::CudaGlobal::
buildCellNeighborhood (const IntVectorType numCell,
		       const IndexType devideLevel,
		       const ScalorType rlist,
		       const HostVectorType boxSize,
		       IndexType * numNeighbor,
		       IndexType * neighborCellIndex,
		       const IndexType stride)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  __shared__ bool stop;

  numNeighbor[bid] = 0;
  int centerx, centery, centerz;
  bool oneCellX(false), oneCellY(false), oneCellZ(false);
  if (numCell.x == 3) oneCellX = true;
  if (numCell.y == 3) oneCellY = true;
  if (numCell.z == 3) oneCellZ = true;

  if (tid == 0){
    stop = false;
    Parallel::CudaDevice::D1toD3 (numCell, int(bid), centerx, centery, centerz);
    if (oneCellX){
      if (centerx == 0 || centerx == 2){
	stop = true;
      }
    }
    else {
      if (centerx < devideLevel || centerx >= numCell.x - devideLevel){
	stop = true;
      }
    }
    if (oneCellY){
      if (centery == 0 || centery == 2){
	stop = true;
      }
    }
    else {
      if (centery < devideLevel || centery >= numCell.y - devideLevel){
	stop = true;
      }
    }
    if (oneCellZ){
      if (centerz == 0 || centerz == 2){
	stop = true;
      }
    }
    else {
      if (centerz < devideLevel || centerz >= numCell.z - devideLevel){
	stop = true;
      }
    }
  }
  __syncthreads();
  if (stop) return;

  if (tid == 0){
    int lowerX (-devideLevel);
    int lowerY (-devideLevel);
    int lowerZ (-devideLevel);
    if (oneCellX) lowerX = -1;
    if (oneCellY) lowerY = -1;
    if (oneCellZ) lowerZ = -1;
    int upperX (devideLevel+1);
    int upperY (devideLevel+1);
    int upperZ (devideLevel+1);
    if (oneCellX) upperX = 2;
    if (oneCellY) upperY = 2;
    if (oneCellZ) upperZ = 2;
    ScalorType scalorx, scalory, scalorz;
    oneCellX ? scalorx = boxSize.x :
	scalorx = boxSize.x / ScalorType(numCell.x - (devideLevel << 1));
    oneCellY ? scalory = boxSize.y :
	scalory = boxSize.y / ScalorType(numCell.y - (devideLevel << 1));
    oneCellZ ? scalorz = boxSize.z :
	scalorz = boxSize.z / ScalorType(numCell.z - (devideLevel << 1));
    
    ScalorType rlist2 = rlist * rlist;
    for (int ix = lowerX; ix < upperX; ++ix){
      for (int iy = lowerY; iy < upperY; ++iy){
	for (int iz = lowerZ; iz < upperZ; ++iz){
	  int myx = ix + int(centerx);
	  int myy = iy + int(centery);
	  int myz = iz + int(centerz);
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
	  if (min < rlist2){
	    IndexType tmp = Parallel::CudaDevice::D3toD1 (numCell, myx, myy, myz);
	    neighborCellIndex[(numNeighbor[bid]++) + bid * stride] = tmp;
	  }
	}
      }
    }
  }
}


__global__ void Parallel::CudaGlobal::
buildCellNeighborhood (const IntVectorType numCell,
		       const IndexType devideLevel,
		       const ScalorType rlist,
		       const HostVectorType globalBoxSize,
		       const int rankx,
		       const int ranky,
		       const int rankz,
		       const int nProcDimx,
		       const int nProcDimy,
		       const int nProcDimz,
		       IndexType * numNeighbor,
		       IndexType * neighborCellIndex,
		       CoordNoiType * neighborShiftNoi,
		       const IndexType stride)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  __shared__ bool stop;

  numNeighbor[bid] = 0;
  int centerx, centery, centerz;
  bool oneCellX(false), oneCellY(false), oneCellZ(false);
  if (numCell.x == 3) oneCellX = true;
  if (numCell.y == 3) oneCellY = true;
  if (numCell.z == 3) oneCellZ = true;

  HostVectorType boxSize;
  if (tid == 0){
    boxSize.x = globalBoxSize.x / nProcDimx;
    boxSize.y = globalBoxSize.y / nProcDimy;
    boxSize.z = globalBoxSize.z / nProcDimz;
    stop = false;
    Parallel::CudaDevice::D1toD3 (numCell, int(bid), centerx, centery, centerz);
    if (oneCellX){
      if (centerx == 0 || centerx == 2){
	stop = true;
      }
    }
    else {
      if (centerx < devideLevel || centerx >= numCell.x - devideLevel){
	stop = true;
      }
    }
    if (oneCellY){
      if (centery == 0 || centery == 2){
	stop = true;
      }
    }
    else {
      if (centery < devideLevel || centery >= numCell.y - devideLevel){
	stop = true;
      }
    }
    if (oneCellZ){
      if (centerz == 0 || centerz == 2){
	stop = true;
      }
    }
    else {
      if (centerz < devideLevel || centerz >= numCell.z - devideLevel){
	stop = true;
      }
    }
  }
  __syncthreads();
  if (stop) return;

  if (tid == 0){
    int lowerX (-devideLevel);
    int lowerY (-devideLevel);
    int lowerZ (-devideLevel);
    if (oneCellX) lowerX = -1;
    if (oneCellY) lowerY = -1;
    if (oneCellZ) lowerZ = -1;
    int upperX (devideLevel+1);
    int upperY (devideLevel+1);
    int upperZ (devideLevel+1);
    if (oneCellX) upperX = 2;
    if (oneCellY) upperY = 2;
    if (oneCellZ) upperZ = 2;
    ScalorType scalorx, scalory, scalorz;
    oneCellX ? scalorx = boxSize.x :
	scalorx = boxSize.x / ScalorType(numCell.x - (devideLevel << 1));
    oneCellY ? scalory = boxSize.y :
	scalory = boxSize.y / ScalorType(numCell.y - (devideLevel << 1));
    oneCellZ ? scalorz = boxSize.z :
	scalorz = boxSize.z / ScalorType(numCell.z - (devideLevel << 1));
    
    ScalorType rlist2 = rlist * rlist;
    CoordNoiType myshift ;
    for (int ix = lowerX; ix < upperX; ++ix){
      for (int iy = lowerY; iy < upperY; ++iy){
	for (int iz = lowerZ; iz < upperZ; ++iz){
	  int myx = ix + int(centerx);
	  int myy = iy + int(centery);
	  int myz = iz + int(centerz);
	  myshift.x = myshift.y = myshift.z = 0;
	  if (rankx == 0 && myx < devideLevel) {
	    myshift.x = - 1;
	  }
	  else if (rankx == nProcDimx - 1 && myx >= int(numCell.x - devideLevel)){
	    myshift.x = 1;
	  }
	  if (ranky == 0 && myy < devideLevel) {
	    myshift.y = - 1;
	  }
	  else if (ranky == nProcDimy - 1 && myy >= int(numCell.y - devideLevel)){
	    myshift.y = 1;
	  }
	  if (rankz == 0 && myz < devideLevel) {
	    myshift.z = - 1;
	  }
	  else if (rankz == nProcDimz - 1 && myz >= int(numCell.z - devideLevel)){
	    myshift.z = 1;
	  }
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
	  if (min < rlist2){
	    IndexType tmp = Parallel::CudaDevice::D3toD1 (numCell, myx, myy, myz);
	    neighborShiftNoi [(numNeighbor[bid]  ) + bid * stride] = myshift;
	    neighborCellIndex[(numNeighbor[bid]++) + bid * stride] = tmp;
	  }
	}
      }
    }
  }
}


void Parallel::DeviceCellRelation::
copyToHost (HostCellRelation & hrelation)
{
  hrelation.easyMalloc (totalNumCell, MaxNeiPerCell);
  cudaMemcpy (hrelation.cptr_numNeighborCell(),
	      numNeighbor,
	      sizeof(IndexType) * totalNumCell,
	      cudaMemcpyDeviceToHost);
  cudaMemcpy (hrelation.cptr_neighborCellIndex(),
	      neighborCellIndex,
	      sizeof(IndexType) * totalNumCell * MaxNeiPerCell,
	      cudaMemcpyDeviceToHost);
  cudaMemcpy (hrelation.cptr_neighborShiftNoi(),
	      neighborShiftNoi,
	      sizeof(CoordNoiType) * totalNumCell * MaxNeiPerCell,
	      cudaMemcpyDeviceToHost);
  checkCUDAError ("DeviceCellRelation::copyToHost");
}


void Parallel::DeviceCellRelation::
rebuild (const DeviceCellListedMDData & list,
	 const SubCellList & sub0,
	 const SubCellList & sub1)
{
  IndexType * tmpHost;
  IndexType * deviceSub0;
  IndexType * deviceSub1;
  size_t size0 = sizeof (IndexType) * sub0.size();
  size_t size1 = sizeof (IndexType) * sub1.size();
      
  tmpHost = (IndexType *) malloc (size0);
  if (tmpHost == NULL){
    throw MDExcptFailedMallocOnHost ("DeviceCellListedMDData::build",
				     "tmpHost", size0);
  }
  cudaMalloc ((void**)&deviceSub0, size0);
  checkCUDAError ("DeviceCellListedMDData::build, malloc deviceSub0");
  for (IndexType i = 0; i < sub0.size(); ++i){
    tmpHost[i] = sub0[i];
  }
  cudaMemcpy (deviceSub0, tmpHost, size0, cudaMemcpyHostToDevice);
  checkCUDAError ("DeviceCellListedMDData::build, copy to deviceSub0");
  free (tmpHost);

  tmpHost = (IndexType *) malloc (size1);
  if (tmpHost == NULL){
    throw MDExcptFailedMallocOnHost ("DeviceCellListedMDData::build",
				     "tmpHost", size1);
  }
  cudaMalloc ((void**)&deviceSub1, size1);
  checkCUDAError ("DeviceCellListedMDData::build, malloc deviceSub1");
  for (IndexType i = 0; i < sub1.size(); ++i){
    tmpHost[i] = sub1[i];
  }
  cudaMemcpy (deviceSub1, tmpHost, size1, cudaMemcpyHostToDevice);
  checkCUDAError ("DeviceCellListedMDData::build, copy to deviceSub1");
  free (tmpHost);

  ptr_list = & list;
  IntVectorType numCell = list.getNumCell ();
  totalNumCell = numCell.x * numCell.y * numCell.z;
  MaxNeiPerCell = 2 * list.getDevideLevel() + 1;
  MaxNeiPerCell = MaxNeiPerCell * MaxNeiPerCell * MaxNeiPerCell;
  int Nx, Ny, Nz;
  Parallel::Interface::numProcDim (Nx, Ny, Nz);  
  easyMalloc (totalNumCell, MaxNeiPerCell);
  int rankx, ranky, rankz;
  Parallel::Interface::rankToCartCoord
      (Parallel::Interface::myRank(), rankx, ranky, rankz);
  
  checkCUDAError ("DeviceCellRelation::build, buildCellNeighborhood");

  Parallel::Auxiliary::setValue
      <<<toGridDim(totalNumCell / DefaultNThreadPerBlock + 1),
      DefaultNThreadPerBlock>>> (
	  numNeighbor, totalNumCell, IndexType(0));
  Parallel::CudaGlobal::buildCellNeighborhood
      <<<sub0.size(), 1>>> (
	  numCell,
	  list.getDevideLevel(),
	  list.getRlist(),
	  list.getGlobalBoxSize(),
	  rankx, ranky, rankz,
	  Nx, Ny, Nz,
	  deviceSub0,
	  sub0.size(),
	  deviceSub1,
	  sub1.size(),
	  numNeighbor,
	  neighborCellIndex,
	  neighborShiftNoi,
	  MaxNeiPerCell);

  cudaFree(deviceSub0);
  cudaFree(deviceSub1);
  checkCUDAError ("DeviceCellListedMDData::build, cuda free");
}


__global__ void Parallel::CudaGlobal::
buildCellNeighborhood (const IntVectorType numCell,
		       const IndexType devideLevel,
		       const ScalorType rlist,
		       const HostVectorType globalBoxSize,
		       const int rankx,
		       const int ranky,
		       const int rankz,
		       const int nProcDimx,
		       const int nProcDimy,
		       const int nProcDimz,
		       const IndexType * subList0,
		       const IndexType length0,
		       const IndexType * subList1,
		       const IndexType length1,
		       IndexType * numNeighbor,
		       IndexType * neighborCellIndex,
		       CoordNoiType * neighborShiftNoi,
		       const IndexType stride)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;

  if (bid >= length0) return ;
  
  bool oneCellX(false), oneCellY(false), oneCellZ(false);
  if (numCell.x == 3) oneCellX = true;
  if (numCell.y == 3) oneCellY = true;
  if (numCell.z == 3) oneCellZ = true;

  HostVectorType boxSize;
  ScalorType scalorx, scalory, scalorz;
  if (tid == 0){
    boxSize.x = globalBoxSize.x / nProcDimx;
    boxSize.y = globalBoxSize.y / nProcDimy;
    boxSize.z = globalBoxSize.z / nProcDimz;
    oneCellX ? scalorx = boxSize.x :
	scalorx = boxSize.x / ScalorType(numCell.x - (devideLevel << 1));
    oneCellY ? scalory = boxSize.y :
	scalory = boxSize.y / ScalorType(numCell.y - (devideLevel << 1));
    oneCellZ ? scalorz = boxSize.z :
	scalorz = boxSize.z / ScalorType(numCell.z - (devideLevel << 1));
  }

  IndexType my_cellIndex = subList0[bid];
  IntScalorType my_cellIndexx, my_cellIndexy, my_cellIndexz;
  if (tid == 0){
    Parallel::CudaDevice::D1toD3 (numCell,
				  int(my_cellIndex),
				  my_cellIndexx, my_cellIndexy, my_cellIndexz);
    ScalorType rlist2 = rlist * rlist;
    for (IndexType ii = 0; ii < length1; ++ii){
      IndexType target_cellIndex = subList1[ii];
      IntScalorType target_cellIndexx, target_cellIndexy, target_cellIndexz;
      Parallel::CudaDevice::D1toD3 (numCell,
				    int(target_cellIndex),
				    target_cellIndexx,
				    target_cellIndexy,
				    target_cellIndexz);
      if (abs(target_cellIndexx - my_cellIndexx) > devideLevel ||
	  abs(target_cellIndexy - my_cellIndexy) > devideLevel ||
	  abs(target_cellIndexz - my_cellIndexz) > devideLevel ){
	continue;
      }
      CoordNoiType myshift ;
      myshift.x = myshift.y = myshift.z = 0;
      if (rankx == 0 &&
	  target_cellIndexx < devideLevel) {
	myshift.x = - 1;
      }
      else if (rankx == nProcDimx - 1 &&
	       target_cellIndexx >= int(numCell.x - devideLevel)){
	myshift.x = 1;
      }
      if (ranky == 0 &&
	  target_cellIndexy < devideLevel) {
	myshift.y = - 1;
      }
      else if (ranky == nProcDimy - 1 &&
	       target_cellIndexy >= int(numCell.y - devideLevel)){
	myshift.y = 1;
      }
      if (rankz == 0 &&
	  target_cellIndexz < devideLevel) {
	myshift.z = - 1;
      }
      else if (rankz == nProcDimz - 1 &&
	       target_cellIndexz >= int(numCell.z - devideLevel)){
	myshift.z = 1;
      }
      ScalorType min = 1e9;
      for (int dx = -1; dx <= 1; ++dx){
	for (int dy = -1; dy <= 1; ++dy){
	  for (int dz = -1; dz <= 1; ++dz){
	    ScalorType diffx ((-my_cellIndexx + target_cellIndexx + dx) * scalorx);
	    ScalorType diffy ((-my_cellIndexy + target_cellIndexy + dy) * scalory);
	    ScalorType diffz ((-my_cellIndexz + target_cellIndexz + dz) * scalorz);
	    ScalorType diff2 (diffx * diffx + diffy * diffy + diffz * diffz);
	    if (diff2 < min){
	      min = diff2;
	    }
	  }
	}
      }
      if (min < rlist2){
	neighborShiftNoi [(numNeighbor[my_cellIndex]  ) + my_cellIndex * stride]
	    = myshift;
	neighborCellIndex[(numNeighbor[my_cellIndex]++) + my_cellIndex * stride]
	    = target_cellIndex;
      }
    }
  }
}
  

// __global__ void Parallel::CudaGlobal::
// buildCellNeighborhood (const IntVectorType numCell,
// 		       const IndexType devideLevel,
// 		       const ScalorType rlist,
// 		       const HostVectorType globalBoxSize,
// 		       const int rankx,
// 		       const int ranky,
// 		       const int rankz,
// 		       const int nProcDimx,
// 		       const int nProcDimy,
// 		       const int nProcDimz,
// 		       const IndexType * subList0,
// 		       const IndexType length0,
// 		       const IndexType * subList1,
// 		       const IndexType length1,
// 		       IndexType * numNeighbor,
// 		       IndexType * neighborCellIndex,
// 		       const IndexType stride)
// {
//   IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
//   IndexType tid = threadIdx.x;

//   if (bid >= length0) return ;
  
//   int centerx, centery, centerz;
  
//   if (tid == 0){
//     HostVectorType boxSize;
//     boxSize.x = globalBoxSize.x / nProcDimx;
//     boxSize.y = globalBoxSize.y / nProcDimy;
//     boxSize.z = globalBoxSize.z / nProcDimz;
//     ScalorType rlist2 = rlist;

//     numNeighbor[bid] = 0;
//     Parallel::CudaDevice::D1toD3 (numCell, int(bid), centerx, centery, centerz);

//     bool oneCellX(false), oneCellY(false), oneCellZ(false);
//     if (numCell.x == 3) oneCellX = true;
//     if (numCell.y == 3) oneCellY = true;
//     if (numCell.z == 3) oneCellZ = true;
//     ScalorType scalorx, scalory, scalorz;
//     oneCellX ? scalorx = boxSize.x :
// 	scalorx = boxSize.x / ScalorType(numCell.x - (devideLevel << 1));
//     oneCellY ? scalory = boxSize.y :
// 	scalory = boxSize.y / ScalorType(numCell.y - (devideLevel << 1));
//     oneCellZ ? scalorz = boxSize.z :
// 	scalorz = boxSize.z / ScalorType(numCell.z - (devideLevel << 1));

//     int targetx, targety, targetz;
//     IndexType targetIndex;
//     for (IndexType i = 0; i < length1; ++i){
//       targetIndex = subList1[i];
//       Parallel::CudaDevice::D1toD3 (numCell, int(targetIndex),
// 				    targetx, targety, targetz);

//       ScalorType min = 1e9;
// #pragma unroll 27
//       for (int dx = -1; dx <= 1; ++dx){
// 	for (int dy = -1; dy <= 1; ++dy){
// 	  for (int dz = -1; dz <= 1; ++dz){
// 	    ScalorType diffx ((-centerx + targetx + dx) * scalorx);
// 	    ScalorType diffy ((-centery + targety + dy) * scalory);
// 	    ScalorType diffz ((-centerz + targetz + dz) * scalorz);
// 	    // shortestImage (box, &diffx, &diffy, &diffz);
// 	    ScalorType diff2 (diffx * diffx + diffy * diffy + diffz * diffz);
// 	    if (diff2 < min){
// 	      min = diff2;
// 	    }
// 	  }
// 	}
//       }
//       if (min < rlist2){
// 	IndexType tmp = Parallel::CudaDevice::D3toD1 (numCell,
// 						      targetx, targety, targetz);
// 	neighborCellIndex[(numNeighbor[bid]++) + bid * stride] = tmp;
//       }
//     }
//   }
// }
  

__global__ void Parallel::CudaGlobal::
rescaleCoordinate (const IndexType * numAtomInCell,
		   const HostVectorType scale,
		   CoordType * coord)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType this_numAtomInCell = numAtomInCell[bid];
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (threadIdx.x < this_numAtomInCell){
    coord[ii].x *= scale.x;
    coord[ii].y *= scale.y;
    coord[ii].z *= scale.z;
  }
}  

__global__ void Parallel::CudaGlobal::
rescaleVelocity (const IndexType * numAtomInCell,
		 const HostVectorType scale,
		 ScalorType * velox,
		 ScalorType * veloy,
		 ScalorType * veloz)
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType this_numAtomInCell = numAtomInCell[bid];
  IndexType ii = threadIdx.x + bid * blockDim.x;

  if (threadIdx.x < this_numAtomInCell){
    velox[ii] *= scale.x;
    veloy[ii] *= scale.y;
    veloz[ii] *= scale.z;
  }
}  


void Parallel::DeviceCellListedMDData::
rescaleCoordinate (const HostVectorType & scale)
{
  HostVectorType globalBoxSize = getGlobalBoxSize();
  IntVectorType numCell = getNumCell();
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  dim3 gridDim = toGridDim (totalNumCell);
  
  Parallel::CudaGlobal::rescaleCoordinate
      <<<gridDim, Parallel::Interface::numThreadsInCell()>>>(
	  numAtomInCell,
	  scale,
	  coord);
  checkCUDAError ("DeviceCellListedMDData::rescaleCoordinate");
  
  globalBoxSize.x *= scale.x;
  globalBoxSize.y *= scale.y;
  globalBoxSize.z *= scale.z;
  setGlobalBox (globalBoxSize.x, globalBoxSize.y, globalBoxSize.z);

  frameLow.x *= scale.x;
  frameLow.y *= scale.y;
  frameLow.z *= scale.z;
  frameUp.x *= scale.x;
  frameUp.y *= scale.y;
  frameUp.z *= scale.z;
}

void Parallel::DeviceCellListedMDData::
rescaleVelocity  (const HostVectorType & scale)
{
  HostVectorType globalBoxSize = getGlobalBoxSize();
  IntVectorType numCell = getNumCell();
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  dim3 gridDim = toGridDim (totalNumCell);

  Parallel::CudaGlobal::rescaleVelocity
      <<<gridDim, Parallel::Interface::numThreadsInCell()>>>(
	  numAtomInCell,
	  scale,
	  velox, veloy, veloz);
  checkCUDAError ("DeviceCellListedMDData::rescaleVelocity");
}


void Parallel::DeviceCellListedMDData::
mallocFromHost (const HostCellListedMDData & hdata)
{
  DeviceMDData::mallocFromHost (hdata);
  
  rlist = hdata.getRlist();
  devideLevel = hdata.getDevideLevel();
  numCell.x = hdata.getNumCell().x;
  numCell.y = hdata.getNumCell().y;
  numCell.z = hdata.getNumCell().z;
  frameLow.x = hdata.getFrameLow().x;
  frameLow.y = hdata.getFrameLow().y;
  frameLow.z = hdata.getFrameLow().z;
  frameUp.x  = hdata.getFrameUp().x;
  frameUp.y  = hdata.getFrameUp().y;
  frameUp.z  = hdata.getFrameUp().z;
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  easyMallocCell (totalNumCell);
  setValue <<<totalNumCell / DefaultNThreadPerBlock + 1, DefaultNThreadPerBlock>>>
      (numAtomInCell, totalNumCell, IndexType(0));
}

void Parallel::DeviceCellListedMDData::
mallocToHost (HostCellListedMDData & hdata) const
{
  DeviceMDData::mallocToHost (hdata);

  hdata.rlist = rlist;
  hdata.devideLevel = devideLevel;
  hdata.frameLow.x = frameLow.x;
  hdata.frameLow.y = frameLow.y;
  hdata.frameLow.z = frameLow.z;
  hdata.frameUp.x  = frameUp.x;
  hdata.frameUp.y  = frameUp.y;
  hdata.frameUp.z  = frameUp.z;
  hdata.numCell.x  = numCell.x;
  hdata.numCell.y  = numCell.y;
  hdata.numCell.z  = numCell.z;

  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  hdata.easyReallocCell (totalNumCell);
  for (IndexType i = 0; i < totalNumCell; ++i){
    hdata.numAtomInCell[i] = 0;
  }
}

void Parallel::DeviceCellListedMDData::
mallocFromDevice (const DeviceCellListedMDData & ddata)
{
  DeviceMDData::mallocFromDevice (ddata);

  rlist = ddata.getRlist();
  devideLevel = ddata.getDevideLevel();
  numCell.x = ddata.getNumCell().x;
  numCell.y = ddata.getNumCell().y;
  numCell.z = ddata.getNumCell().z;
  frameLow.x = ddata.getFrameLow().x;
  frameLow.y = ddata.getFrameLow().y;
  frameLow.z = ddata.getFrameLow().z;
  frameUp.x  = ddata.getFrameUp().x;
  frameUp.y  = ddata.getFrameUp().y;
  frameUp.z  = ddata.getFrameUp().z;

  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  easyMallocCell (totalNumCell);
  setValue <<<totalNumCell / DefaultNThreadPerBlock + 1, DefaultNThreadPerBlock>>>
      (numAtomInCell, totalNumCell, IndexType(0));
}



