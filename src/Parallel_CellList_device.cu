#define DEVICE_CODE

#include "Parallel_CellList.h"
#include "Parallel_Interface.h"
#include "Parallel_CellList_device.h"
#include "Parallel_Algorithm.h"
#include "Auxiliary.h"

#include "compile_error_mixcode.h"

void Parallel::DeviceCellListedMDData::
initZeroCell ()
{
  IndexType numThreadBlock = Parallel::Interface::numThreadsInCell();
  dim3 gridDim = toGridDim(numCell.x*numCell.y*numCell.z);
  Parallel::CudaGlobal::initZeroCell
      <<<gridDim, numThreadBlock>>>(
	  numCell, 
	  numAtomInCell,
	  numNeighborCell);
}

void Parallel::DeviceCellListedMDData::
initCellStructure (const ScalorType & rlist,
		   const IndexType & devideLevel_,
		   const BoxDirection_t & bdir)
{
  int Nx, Ny, Nz;
  Parallel::Interface::numProcDim (Nx, Ny, Nz);
  int ix, iy, iz;
  Parallel::Interface::rankToCartCoord (Parallel::Interface::myRank(), ix, iy, iz);
  double dx, dy, dz;
  dx = getGlobalBoxSize().x / double(Nx);
  dy = getGlobalBoxSize().x / double(Nx);
  dz = getGlobalBoxSize().x / double(Nx);
  frameLow.x = dx * ix;
  frameLow.y = dy * iy;
  frameLow.z = dz * iz;
  frameUp.x = frameLow.x + dx;
  frameUp.y = frameLow.y + dy;
  frameUp.z = frameLow.z + dz;
  
  bool CellOnX, CellOnY, CellOnZ;
  CellOnX = bdir & RectangularBoxGeometry::mdRectBoxDirectionX;
  CellOnY = bdir & RectangularBoxGeometry::mdRectBoxDirectionY;
  CellOnZ = bdir & RectangularBoxGeometry::mdRectBoxDirectionZ;
  double rlisti = 1./rlist;

  if (CellOnX ) numCell.x = int ( floor(getGlobalBoxSize().x * rlisti) );
  else numCell.x = 1;
  if (CellOnY ) numCell.y = int ( floor(getGlobalBoxSize().y * rlisti) );
  else numCell.y = 1;
  if (CellOnZ ) numCell.z = int ( floor(getGlobalBoxSize().z * rlisti) );
  else numCell.z = 1;

  if ((CellOnX && numCell.x < 3) ||
      (CellOnY && numCell.y < 3) ||
      (CellOnZ && numCell.z < 3) ){
    throw MDExcptCellList ("Number of cell on one direction is less than 3");
  }

  // add ghost cell
  VectorType dcell;
  dcell.x = (frameUp.x - frameLow.x) / numCell.x;
  dcell.y = (frameUp.y - frameLow.y) / numCell.y;
  dcell.z = (frameUp.z - frameLow.z) / numCell.z;
  frameUp.x += dcell.x;
  frameUp.y += dcell.y;
  frameUp.z += dcell.z;
  frameLow.x -= dcell.x;
  frameLow.y -= dcell.y;
  frameLow.z -= dcell.z;
  numCell.x += 2;
  numCell.y += 2;
  numCell.z += 2;
  
  devideLevel = devideLevel_;
  if (CellOnX) numCell.x *= devideLevel;
  if (CellOnY) numCell.y *= devideLevel;
  if (CellOnZ) numCell.z *= devideLevel;

  DeviceMDData bkData (*this);
  
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell ();
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  maxNumNeighborCell = 1;
  if (CellOnX) maxNumNeighborCell *= devideLevel * 2 + 1;
  if (CellOnY) maxNumNeighborCell *= devideLevel * 2 + 1;
  if (CellOnZ) maxNumNeighborCell *= devideLevel * 2 + 1;
  
  if (numThreadsInCell * totalNumCell > memSize_){
    mallocAll (numThreadsInCell * totalNumCell);
    initZero();
  }
  mallocCell (totalNumCell, maxNumNeighborCell);
  initZeroCell ();

  IndexType numThreadBlock = numThreadsInCell;
  dim3 gridDim = toGridDim(numCell.x*numCell.y*numCell.z);

  Parallel::CudaGlobal::formCellStructure
      <<<gridDim, numThreadBlock >>>(
	  frameLow,
	  frameUp,
	  numCell,
	  numAtomInCell,
	  numAtom_,
	  bkData.dptr_coordinate(),
	  bkData.dptr_coordinateNoiX(),
	  bkData.dptr_coordinateNoiY(),
	  bkData.dptr_coordinateNoiZ(),
	  bkData.dptr_velocityX(),
	  bkData.dptr_velocityY(),
	  bkData.dptr_velocityZ(),
	  bkData.dptr_globalIndex(),
	  bkData.dptr_type(),
	  bkData.dptr_mass(),
	  bkData.dptr_charge(),
	  coord,
	  coordNoix,
	  coordNoiy,
	  coordNoiz,
	  velox,
	  veloy,
	  veloz,
	  globalIndex,
	  type,
	  mass,
	  charge,
	  err.ptr_de);
  checkCUDAError ("Parallel::formCellStructure");
  err.updateHost();
  err.check ("Parallel::formCellSturcture");
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
  checkCUDAError ("Parallel::rebuild malloc backup");
  
  Parallel::CudaGlobal::rebuildCellList_step1
      <<<gridDim, numThreadBlock>>> (
	  frameLow,
	  frameUp,
	  numCell,
	  bk_numAtomInCell,
	  numAtomInCell,
	  coord,
	  coordNoix,
	  coordNoiy,
	  coordNoiz,
	  velox,
	  veloy,
	  veloz,
	  globalIndex,
	  type,
	  mass,
	  charge,
	  err.ptr_de);
  checkCUDAError ("Parallel::rebuild step1");
  err.updateHost();
  err.check ("Parallel::rebuild step1");
  Parallel::CudaGlobal::rebuildCellList_step2
      <<<gridDim, numThreadBlock, numThreadBlock*sizeof(IndexType)*3>>> (
	  numAtomInCell,
	  coord,
	  coordNoix,
	  coordNoiy,
	  coordNoiz,
	  velox,
	  veloy,
	  veloz,
	  globalIndex,
	  type,
	  mass,
	  charge,
	  err.ptr_de);
  checkCUDAError ("Parallel::rebuild step2");
  err.updateHost();
  err.check ("Parallel::rebuild step2");
  cudaFree (bk_numAtomInCell);
  checkCUDAError ("Parallel::rebuild free backup");  
}
	  

__global__ void Parallel::CudaGlobal::
initZeroCell (const IntVectorType numCell,
	      IndexType * numAtomInCell,
	      IndexType * numNeighborCell )
{
  IndexType bid = blockIdx.x + gridDim.x * blockIdx.y;
  IndexType tid = threadIdx.x;
  IndexType ii = tid + bid * blockDim.x;
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;

  if (ii < totalNumCell){
    numAtomInCell[ii] = 0;
    numNeighborCell[ii] = 0;
  }
}


__global__ void Parallel::CudaGlobal::
formCellStructure (const VectorType frameLow,
		   const VectorType frameUp,
		   const IntVectorType numCell,
		   IndexType * numAtomInCell,
		   const IndexType numAtom,
		   const CoordType  * bk_coord,
		   const IntScalorType * bk_coordNoix,
		   const IntScalorType * bk_coordNoiy,
		   const IntScalorType * bk_coordNoiz,
		   const ScalorType * bk_velox,
		   const ScalorType * bk_veloy,
		   const ScalorType * bk_veloz,
		   const IndexType  * bk_globalIndex,
		   const TypeType   * bk_type,
		   const ScalorType * bk_mass,
		   const ScalorType * bk_charge,
		   CoordType  * coord,
		   IntScalorType * coordNoix,
		   IntScalorType * coordNoiy,
		   IntScalorType * coordNoiz,
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
    if (targetCellx == numCell.x){
      targetCellx -= numCell.x;
    }
    if (targetCelly == numCell.y){
      targetCelly -= numCell.y;
    }
    if (targetCellz == numCell.z){
      targetCellz -= numCell.z;
    }
    if (ptr_de != NULL && 
	(targetCellx >= numCell.x || 
	 targetCelly >= numCell.y || 
	 targetCellz >= numCell.z)){
      *ptr_de = mdErrorOverFlowCellIdx;
      return;
    }

    IndexType cellid = CudaDevice::D3toD1
	(numCell, targetCellx, targetCelly, targetCellz);

    IndexType pid = atomicInc (&numAtomInCell[cellid], blockDim.x);
    targetIndex = pid + cellid * blockDim.x;
    coord[targetIndex] = bk_coord[ii];
    coordNoix[targetIndex] = bk_coordNoix[ii];
    coordNoiy[targetIndex] = bk_coordNoiy[ii];
    coordNoiz[targetIndex] = bk_coordNoiz[ii];
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
rebuildCellList_step1 (const VectorType frameLow,
		       const VectorType frameUp,
		       const IntVectorType numCell,
		       const IndexType * bk_numAtomInCell,
		       IndexType * numAtomInCell,
		       CoordType  * coord,
		       IntScalorType * coordNoix,
		       IntScalorType * coordNoiy,
		       IntScalorType * coordNoiz,
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

  // IndexType mark = MaxIndexValue - (MaxIndexValue >> 1);
  // IndexType mystat ;
  
  if (tid < bk_numAtomInCell[bid]){
    // mystat = globalIndex[ii];
    IndexType targetCellx, targetCelly, targetCellz;
    targetCellx = IndexType((coord[ii].x - frameLow.x) * dcellxi);
    targetCelly = IndexType((coord[ii].y - frameLow.y) * dcellyi);
    targetCellz = IndexType((coord[ii].z - frameLow.z) * dcellzi);
    printf ("%d %d %d %d %f %f %f\n", ii, targetCellx, targetCelly, targetCellz,
	    coord[ii].x, coord[ii].y, coord[ii].z);
    if (targetCellx == numCell.x){
      targetCellx -= numCell.x;
    }
    if (targetCelly == numCell.y){
      targetCelly -= numCell.y;
    }
    if (targetCellz == numCell.z){
      targetCellz -= numCell.z;
    }
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
      IndexType pid = atomicAdd (&numAtomInCell[cellid], 1);
      if (pid == blockDim.x){
	*ptr_de = mdErrorShortCellList;
	pid -= blockDim.x;
      }
      IndexType targetIndex = pid + cellid * blockDim.x;
      coord[targetIndex] = coord[ii];
      coordNoix[targetIndex] = coordNoix[ii];
      coordNoiy[targetIndex] = coordNoiy[ii];
      coordNoiz[targetIndex] = coordNoiz[ii];
      velox[targetIndex] = velox[ii];
      veloy[targetIndex] = veloy[ii];
      veloz[targetIndex] = veloz[ii];
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
		       IntScalorType * coordNoix,
		       IntScalorType * coordNoiy,
		       IntScalorType * coordNoiz,
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
  
  if (tid < total){
    IndexType fromId = myIndex[tid] + bid * blockDim.x;
    if (ii != fromId){
      coord[ii] = coord[fromId];
      coordNoix[ii] = coordNoix[fromId];
      coordNoiy[ii] = coordNoiy[fromId];
      coordNoiz[ii] = coordNoiz[fromId];
      velox[ii] = velox[fromId];
      veloy[ii] = veloy[fromId];
      veloz[ii] = veloz[fromId];
      globalIndex[ii] = globalIndex[fromId];
      type[ii] = type[fromId];
      mass[ii] = mass[fromId];
      charge[ii] = charge[fromId];
    }
  }
  // else {
  //   globalIndex[ii] = MaxIndexValue;
  // }  

  if (tid == 0){
    numAtomInCell[bid] = total;
  }
}


void Parallel::DeviceCellListedMDData::
mallocCell (const IndexType & totalNumCell,
	    const IndexType & maxNumNeighborCell_)
{
  maxNumNeighborCell = maxNumNeighborCell_;
  if (malloced){
    clearCell ();
  }
  cudaMalloc ((void**)&numAtomInCell, sizeof(IndexType) * totalNumCell);
  cudaMalloc ((void**)&numNeighborCell, sizeof(IndexType) * totalNumCell);
  cudaMalloc ((void**)&neighborCellIndex,
	      sizeof(IndexType) * totalNumCell * maxNumNeighborCell);
  checkCUDAError ("malloc Cell");
  malloced = true;
}

void Parallel::DeviceCellListedMDData::
clearCell()
{
  if (malloced){
    cudaFree (numAtomInCell);
    cudaFree (numNeighborCell);
    cudaFree (neighborCellIndex);
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
  maxNumNeighborCell = 0;
  malloced = false;
}

Parallel::DeviceCellListedMDData::
~DeviceCellListedMDData()
{
  clearCell();
}

Parallel::DeviceSubCellList::
DeviceSubCellList ()
    : ptr_data (NULL)
{
}


void Parallel::DeviceSubCellList::
bond (DeviceCellListedMDData & ddata)
{
  ptr_data = &ddata;
}

void Parallel::DeviceSubCellList::
build ()
{
  Parallel::Interface::sort (this->begin(), this->end());
}

bool Parallel::DeviceSubCellList::
isBuilt ()
{
  return (ptr_data != NULL) &&
      (Parallel::Interface::is_sorted (this->begin(), this->end()));
}

void Parallel::DeviceSubCellList::
add (const DeviceSubCellList & a)
{
  for (std::vector<IndexType>::const_iterator it = a.begin();
       it != a.end(); ++it){
    push_back (*it);
  }
  Parallel::Interface::unique (this->begin(), this->end());
  Parallel::Interface::sort   (this->begin(), this->end());
}

void Parallel::DeviceSubCellList::
sub (const DeviceSubCellList & a)
{
  std::vector<IndexType > result;
  Parallel::Interface::set_difference (this->begin(), this->end(),
				       a.begin(), a.end(),
				       result.begin());
  this->std::vector<IndexType>::operator= (result);
}



void Parallel::DeviceCellListedMDData::
buildSubList (const IndexType & xIdLo,
	      const IndexType & xIdUp,
	      const IndexType & yIdLo,
	      const IndexType & yIdUp,
	      const IndexType & zIdLo,
	      const IndexType & zIdUp,
	      DeviceSubCellList & subList)
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
  subList.bond(*this);
  
  for (IndexType i = xIdLo; i < xIdUp; ++i){
    for (IndexType j = yIdLo; j < yIdUp; ++j){
      for (IndexType k = zIdLo; k < zIdUp; ++k){
	subList.push_back ( D3toD1 (i, j, k));
      }
    }
  }
}




