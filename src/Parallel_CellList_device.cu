#define DEVICE_CODE

#include "Parallel_CellList.h"
#include "Parallel_Interface.h"
#include "Parallel_CellList_device.h"
#include "Parallel_Algorithm.h"
#include "Auxiliary.h"
#include "Parallel_Timer.h"

#include "compile_error_mixcode.h"

void Parallel::DeviceCellListedMDData::
initZeroCell ()
{
  IndexType numThreadBlock = Parallel::Interface::numThreadsInCell();
  dim3 gridDim = toGridDim(numCell.x*numCell.y*numCell.z);
  Parallel::CudaGlobal::initZeroCell
      <<<gridDim, numThreadBlock>>>(
	  numCell, 
	  numAtomInCell);
  checkCUDAError ("DeviceCellListedMDData::initZeroCell");
}

void Parallel::DeviceCellListedMDData::
initCellStructure (const ScalorType & rlist_,
		   const IndexType & devideLevel_,
		   const BoxDirection_t & bdir)
{
  int Nx, Ny, Nz;
  Parallel::Interface::numProcDim (Nx, Ny, Nz);
  int ix, iy, iz;
  Parallel::Interface::rankToCartCoord (Parallel::Interface::myRank(), ix, iy, iz);
  double dx, dy, dz;
  dx = getGlobalBoxSize().x / double(Nx);
  dy = getGlobalBoxSize().y / double(Ny);
  dz = getGlobalBoxSize().z / double(Nz);
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
  rlist = rlist_;
  double rlisti = 1./rlist;

  if (CellOnX ) numCell.x = int ( floor(dx * rlisti) );
  else numCell.x = 1;
  if (CellOnY ) numCell.y = int ( floor(dy * rlisti) );
  else numCell.y = 1;
  if (CellOnZ ) numCell.z = int ( floor(dz * rlisti) );
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
  // maxNumNeighborCell = 1;
  // if (CellOnX) maxNumNeighborCell *= devideLevel * 2 + 1;
  // if (CellOnY) maxNumNeighborCell *= devideLevel * 2 + 1;
  // if (CellOnZ) maxNumNeighborCell *= devideLevel * 2 + 1;
  
  if (numThreadsInCell * totalNumCell > DeviceMDData::memSize_){
    DeviceMDData::easyMalloc (numThreadsInCell * totalNumCell);
    DeviceMDData::initZero();
  }
  numData_ = totalNumCell * numThreadsInCell;

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
  checkCUDAError ("Parallel::rebuild step1");
  err.updateHost();
  err.check ("Parallel::rebuild step1");
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
  checkCUDAError ("Parallel::rebuild step2");
  err.updateHost();
  err.check ("Parallel::rebuild step2");
  cudaFree (bk_numAtomInCell);
  checkCUDAError ("Parallel::rebuild free backup");  
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


void Parallel::DeviceCellListedMDData::
easyMallocCell (const IndexType & totalNumCell)
{
  if (totalNumCell == 0) return;
  // if (totalNumCell == numCell.x * numCell.y * numCell.z) return;
  // maxNumNeighborCell = maxNumNeighborCell_;
  clearCell ();
  cudaMalloc ((void**)&numAtomInCell, sizeof(IndexType) * totalNumCell);
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
	      SubCellList & subList)
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
  this->numData() = hcellStartIndex[numCell];
  cudaMemcpy (cellStartIndex, hcellStartIndex, (numCell+1) * sizeof(IndexType),
	      cudaMemcpyHostToDevice);
  checkCUDAError ("DeviceTransferPackage::pack cpy cellStartIndex to device");
  
  free (numAtomInCell);

  this->DeviceMDData::setGlobalBox (ddata.getGlobalBox());
  if (this->DeviceMDData::numData() > this->DeviceMDData::memSize()){
    printf ("# DeviceTransferPackage::pack, realloc\n");
    
    this->DeviceMDData::easyMalloc (this->DeviceMDData::numData() * MemAllocExtension);
  }

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
  checkCUDAError ("DeviceTransferPackage::unpack_replace");
}


void Parallel::DeviceTransferPackage::
unpack_add (DeviceCellListedMDData & ddata) const
{
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
	  ddata.dptr_charge(),
	  err.ptr_de);
  checkCUDAError ("DeviceTransferPackage::unpack_add");
  err.updateHost ();
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

void Parallel::DeviceCellRelation::
easyMalloc (const IndexType & totalNumCell,
	    const IndexType & MaxNeiPerCell)
{
  if (malloced){
    clear ();
  }
  cudaMalloc ((void**)&numNeighbor, totalNumCell * sizeof(IndexType));
  cudaMalloc ((void**)&neighborCellIndex,
	      totalNumCell * MaxNeiPerCell * sizeof (IndexType));
  cudaMalloc ((void**)&neighborShift,
	      totalNumCell * MaxNeiPerCell * sizeof (CoordType));
  checkCUDAError ("DeviceCellRelation::easyMalloc");
  malloced = true;
}

void Parallel::DeviceCellRelation::
clear ()
{
  cudaFree (numNeighbor);
  cudaFree (neighborCellIndex);
  cudaFree (neighborShift);
  malloced = false;
}

Parallel::DeviceCellRelation::
~DeviceCellRelation ()
{
  clear();
}

void Parallel::DeviceCellRelation::
build (DeviceCellListedMDData & list)
{
  ptr_list = &list;

  IntVectorType numCell = list.getNumCell ();
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  MaxNeiPerCell = 2 * list.getDevideLevel() + 1;
  MaxNeiPerCell = MaxNeiPerCell * MaxNeiPerCell * MaxNeiPerCell;
  int Nx, Ny, Nz;
  Parallel::Interface::numProcDim (Nx, Ny, Nz);
  
  easyMalloc (totalNumCell, MaxNeiPerCell);
  // setValue <<<
  //     (totalNumCell + DefaultNThreadPerBlock - 1) / DefaultNThreadPerBlock,
  //     DefaultNThreadPerBlock >>> (
  // 	  numNeighbor,
  // 	  totalNumCell,
  // 	  0);
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
	  neighborShift,
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
		       CoordType * neighborShift,
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
    CoordType myshift ;
    myshift.x = myshift.y = myshift.z = 0;
    for (int ix = lowerX; ix < upperX; ++ix){
      for (int iy = lowerY; iy < upperY; ++iy){
	for (int iz = lowerZ; iz < upperZ; ++iz){
	  int myx = ix + int(centerx);
	  int myy = iy + int(centery);
	  int myz = iz + int(centerz);
	  myshift.x = myshift.y = myshift.z = 0.f;
	  if (rankx == 0 && myx < devideLevel) {
	    myshift.x = - globalBoxSize.x;
	  }
	  else if (rankx == nProcDimx - 1 && myx >= int((numCell.x-1) * devideLevel)){
	    myshift.x = globalBoxSize.x;
	  }
	  if (ranky == 0 && myy < devideLevel) {
	    myshift.y = - globalBoxSize.y;
	  }
	  else if (ranky == nProcDimy - 1 && myy >= int((numCell.y-1) * devideLevel)){
	    myshift.y = globalBoxSize.y;
	  }
	  if (rankz == 0 && myz < devideLevel) {
	    myshift.z = - globalBoxSize.z;
	  }
	  else if (rankz == nProcDimz - 1 && myz >= int((numCell.z-1) * devideLevel)){
	    myshift.z = globalBoxSize.z;
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
	    neighborShift    [(numNeighbor[bid]  ) + bid * stride] = myshift;
	    neighborCellIndex[(numNeighbor[bid]++) + bid * stride] = tmp;
	  }
	}
      }
    }
  }
}

  

  






