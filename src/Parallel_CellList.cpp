#define MPI_CODE

#include "Parallel_CellList.h"
#include <algorithm>
#include "Parallel_Interface.h"

#include "compile_error_mixcode.h"

Parallel::HostCellListedMDData::
HostCellListedMDData ()
{
  rlist = 0;
  devideLevel = 0;
  frameLow.x = frameLow.y = frameLow.z = 0;
  frameUp.x  = frameUp.y  = frameUp.z  = 0;
  numCell.x  = numCell.y  = numCell.z  = 0;
  // maxNumNeighborCell = 0;

  numAtomInCell = NULL;
  // numNeighborCell = NULL;
  // maxNumNeighborCell = NULL;
}

Parallel::HostCellListedMDData::
~HostCellListedMDData ()
{
  freeAPointer ((void**)&numAtomInCell);
  // freeAPointer ((void**)&numNeighborCell);
  // freeAPointer ((void**)&maxNumNeighborCell);
}

void Parallel::HostCellListedMDData::
reallocAll (const IndexType & totalNumCell,
	    const IndexType & maxNumNeighborCell)
{
  numAtomInCell = (IndexType *) realloc (
      numAtomInCell, totalNumCell * sizeof(IndexType));
  if (numAtomInCell == NULL){
    throw MDExcptFailedReallocOnHost ("HostCellListedMDData::realloc",
				      "numAtomInCell",
				      totalNumCell * sizeof(IndexType));
  }
  // numNeighborCell = (IndexType *) realloc (
  //     numNeighborCell, totalNumCell * sizeof(IndexType));
  // if (numNeighborCell == NULL){
  //   throw MDExcptFailedReallocOnHost ("HostCellListedMDData::realloc",
  // 				      "numNeighborCell",
  // 				      totalNumCell * sizeof(IndexType));
  // }
  // neighborCellIndex = (IndexType *) realloc (
  //     neighborCellIndex, totalNumCell * maxNumNeighborCell * sizeof(IndexType));
  // if (neighborCellIndex == NULL){
  //   throw MDExcptFailedReallocOnHost ("HostCellListedMDData::realloc",
  // 				      "neighborCellIndex",
  // 				      totalNumCell * sizeof(IndexType));
  // }
}

  

// void Parallel::HostCellList::
// formCellStructure (const ScalorType & rlist,
// 		   const IndexType & devideLevel_,
// 		   const BoxDirection_t & bdir)
// {
//   if (myData == NULL) return;
  
//   bool CellOnX, CellOnY, CellOnZ;
//   CellOnX = bdir & RectangularBoxGeometry::mdRectBoxDirectionX;
//   CellOnY = bdir & RectangularBoxGeometry::mdRectBoxDirectionY;
//   CellOnZ = bdir & RectangularBoxGeometry::mdRectBoxDirectionZ;
//   double rlisti = 1./rlist;

//   if (CellOnX ) numCell.x = int ( floor(myData.getGlobalBox().size.x * rlisti) );
//   else numCell.x = 1;
//   if (CellOnY ) numCell.y = int ( floor(myData.getGlobalBox().size.y * rlisti) );
//   else numCell.y = 1;
//   if (CellOnZ ) numCell.z = int ( floor(myData.getGlobalBox().size.z * rlisti) );
//   else numCell.z = 1;

//   if ((CellOnX && numCell.x < 3) ||
//       (CellOnY && numCell.y < 3) ||
//       (CellOnZ && numCell.z < 3) ){
//     throw MDExcptTooFewCells ();
//   }

//   devideLevel = devideLevel_;
//   if (CellOnX) numCell.x *= devideLevel;
//   if (CellOnY) numCell.y *= devideLevel;
//   if (CellOnZ) numCell.z *= devideLevel;

//   HostMDData bkData (myData);

//   int Nx, Ny, Nz;
//   Parallel::Interface::numProcDim (Nx, Ny, Nz);
//   int ix, iy, iz;
//   Parallel::Interface::randToCartCoord (Parallel::Interface::myRank(), ix, iy, iz);

//   dcell.x = myData.getGlobalBox().size.x / ScalorType(numCell.x);
//   dcell.y = myData.getGlobalBox().size.y / ScalorType(numCell.y);
//   dcell.z = myData.getGlobalBox().size.z / ScalorType(numCell.z);  
//   frameLow.x = ix * dcell.x;
//   frameLow.y = iy * dcell.y;
//   frameLow.z = iz * dcell.z;
//   frameUp.x = frameLow.x + dcell.x;
//   frameUp.y = frameLow.y + dcell.y;
//   frameUp.z = frameLow.z + dcell.z;
  
//   IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell ();

//   IndexType totalNumCell;
//   maxNumNeighborCell = 1;
//   if (CellOnX) maxNumNeighborCell *= devideLevel * 2 + 1;
//   if (CellOnY) maxNumNeighborCell *= devideLevel * 2 + 1;
//   if (CellOnZ) maxNumNeighborCell *= devideLevel * 2 + 1;
  
//   reallocAll (totalNumCell, numThreadsInCell, maxNumNeighborCell);

//   for (IndexType i = 0; i < totalNumCell; ++i){
//     cellStartIndex[i] = i * numThreadsInCell;
//     numAtomInCell[i] = 0;
//     numNeighborCell[i] = 0;
//   }

//   for (IndexType ii = 0; ii < bkData.numAtom(); ++ii){
//     IndexType cellIdx, cellIdy, cellIdz;
//     if (bkData.cptr_coordinate()[ii].x <  frameLow.x ||
// 	bkData.cptr_coordinate()[ii].x >= frameUp.x){
//       throw MDExcptWrongAtomInProc ();
//     }
//     cellIdx = IndexType (double(bkData.cptr_coordinate()[ii] - frameLow.x) /
// 			 double(dcell.x) * numCell.x );
//     if (cellIdx >= numCell.x) {
//       throw MDExcptWrongCellIndex ();
//     }
//     if (bkData.cptr_coordinate()[ii].y <  frameLow.y ||
// 	bkData.cptr_coordinate()[ii].y >= frameUp.y){
//       throw MDExcptWrongAtomInProc ();
//     }
//     cellIdy = IndexType (double(bkData.cptr_coordinate()[ii] - frameLow.y) /
// 			 double(dcell.y) * numCell.y );
//     if (cellIdy >= numCell.y) {
//       throw MDExcptWrongCellIndex ();
//     }
//     if (bkData.cptr_coordinate()[ii].z <  frameLow.z ||
// 	bkData.cptr_coordinate()[ii].z >= frameUp.z){
//       throw MDExcptWrongAtomInProc ();
//     }
//     cellIdz = IndexType (double(bkData.cptr_coordinate()[ii] - frameLow.z) /
// 			 double(dcell.z) * numCell.z );
//     if (cellIdz >= numCell.z) {
//       throw MDExcptWrongCellIndex ();
//     }
    
//     IndexType cellIndex;
//     D3ToD1 (cellIdx, cellIdy, cellIdz, cellIndex);
//     IndexType atomIndex = cellIndex * numThreadsInCell + (numAtomInCell[cellIndex]++);
//     myData.cptr_coordinate()[atomIndex] = bkData.cptr_coordinate()[ii];
//     myData.cptr_coordinateNoiX()[atomIndex] = bkData.cptr_coordinateNoiX()[ii];
//     myData.cptr_coordinateNoiY()[atomIndex] = bkData.cptr_coordinateNoiY()[ii];
//     myData.cptr_coordinateNoiZ()[atomIndex] = bkData.cptr_coordinateNoiZ()[ii];
//     myData.cptr_velocityX()[atomIndex] = bkData.cptr_velocityX()[ii];
//     myData.cptr_velocityY()[atomIndex] = bkData.cptr_velocityY()[ii];
//     myData.cptr_velocityZ()[atomIndex] = bkData.cptr_velocityZ()[ii];
//     myData.cptr_globalIndex()[atomIndex] = bkData.cptr_globalIndex()[ii];
//     myData.cptr_type()[atomIndex] = bkData.cptr_type()[ii];
//     myData.cptr_mass()[atomIndex] = bkData.cptr_mass()[ii];
//     myData.cptr_charge()[atomIndex] = bkData.cptr_charge()[ii];
//   }
// }




Parallel::HostTransferPackage::
HostTransferPackage ()
    : numCell (0), memSize(0), cellIndex(NULL), cellStartIndex(NULL),
      myMask (MDDataItemMask_All)
{
}

void Parallel::HostTransferPackage::
clearMe ()
{
  if (memSize != 0){
    freeAPointer ((void**)&cellIndex);
    freeAPointer ((void**)&cellStartIndex);
    memSize = 0;
    numCell = 0;
  }
}

Parallel::HostTransferPackage::
~HostTransferPackage ()
{
  clearMe();
}

void Parallel::HostTransferPackage::
easyMallocMe (IndexType memSize_)
{
  clearMe ();
  memSize = memSize_;
  size_t size = memSize * sizeof(IndexType);
  size_t size1 = (memSize+1) * sizeof(IndexType);
  cellIndex = (IndexType *) malloc (size);
  if (cellIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostTransferPackage::reinit",
				     "cellIndex", size);
  }
  cellStartIndex = (IndexType *) malloc (size1);
  if (cellStartIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostTransferPackage::reinit",
				     "cellStartIndex", size1);
  }
}

void Parallel::HostTransferPackage::
reinit (const SubCellList & subCellList)
{
  if (memSize < subCellList.size()){
    easyMallocMe (subCellList.size());
  }
  numCell = subCellList.size();
  for (IndexType i = 0; i < numCell; ++i){
    cellIndex[i] = subCellList[i];
  }
}

void Parallel::HostTransferPackage::
pack (const HostCellListedMDData & hdata,
      const MDDataItemMask_t mask)
{
  if (numCell == 0) return;
  myMask = mask;
  
  IndexType totalNumCell = hdata.numCell.x * hdata.numCell.y * hdata.numCell.z;  
  
  cellStartIndex[0] = 0;
  for (IndexType i = 1; i < numCell+1; ++i){
    cellStartIndex[i] = cellStartIndex[i-1] + hdata.numAtomInCell[cellIndex[i-1]];
  }
  this->numData() = cellStartIndex[numCell];  

  this->HostMDData::setGlobalBox (hdata.getGlobalBox());
  if (this->HostMDData::numData() > this->HostMDData::memSize()){
    this->HostMDData::reallocAll (this->HostMDData::numData() * 2);
  }

  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  IndexType toid = 0;
  for (IndexType i = 0; i < numCell; ++i){
    IndexType cellid = cellIndex[i];
    IndexType fromid = cellid * numThreadsInCell;
    for (IndexType j = 0; j < hdata.numAtomInCell[cellid]; ++j){
      this->cptr_coordinate()[toid++] = hdata.cptr_coordinate()[fromid++];
    }
  }  
}




