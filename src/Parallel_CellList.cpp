#define MPI_CODE

#include "Parallel_CellList.h"
#include <algorithm>
#include "Parallel_Interface.h"

#include "compile_error_mixcode.h"

Parallel::HostCellListedMDData::
HostCellListedMDData ()
    : memSize (0)
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
easyReallocCell (const IndexType & totalNumCell)
{
  numAtomInCell = (IndexType *) realloc (
      numAtomInCell, totalNumCell * sizeof(IndexType));
  if (numAtomInCell == NULL){
    throw MDExcptFailedReallocOnHost ("HostCellListedMDData::realloc",
				      "numAtomInCell",
				      totalNumCell * sizeof(IndexType));
  }
  memSize = totalNumCell;
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


void Parallel::HostCellListedMDData::
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


void Parallel::HostCellListedMDData::
buildSubList (const IndexType & xIdLo,
	      const IndexType & xIdUp,
	      const IndexType & yIdLo,
	      const IndexType & yIdUp,
	      const IndexType & zIdLo,
	      const IndexType & zIdUp,
	      HostSubCellList & subList)
{
  subList.setHostData (*this);
  buildSubList (xIdLo, xIdUp, yIdLo, yIdUp, zIdLo, zIdUp, (SubCellList &)(subList));
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
  if (memSize == memSize_) return;
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
    this->HostMDData::easyRealloc (this->HostMDData::numData() * MemAllocExtension);
  }

  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  // IndexType toid = 0;
  for (IndexType i = 0; i < numCell; ++i){
    IndexType toid = cellStartIndex[i];
    IndexType cellid = cellIndex[i];
    IndexType fromid = cellid * numThreadsInCell;
    IndexType numInThisCell = hdata.numAtomInCell[cellid];
    if (mask & MDDataItemMask_Coordinate){
      for (IndexType j = 0; j < numInThisCell; ++j){
	this->cptr_coordinate()[toid+j] = hdata.cptr_coordinate()[fromid+j];
      }
    }
    if (mask & MDDataItemMask_CoordinateNoi){
      for (IndexType j = 0; j < numInThisCell; ++j){
	this->cptr_coordinateNoiX()[toid+j] = hdata.cptr_coordinateNoiX()[fromid+j];
	this->cptr_coordinateNoiY()[toid+j] = hdata.cptr_coordinateNoiY()[fromid+j];
	this->cptr_coordinateNoiZ()[toid+j] = hdata.cptr_coordinateNoiZ()[fromid+j];
      }
    }
    if (mask & MDDataItemMask_Velocity){
      for (IndexType j = 0; j < numInThisCell; ++j){
	this->cptr_velocityX()[toid+j] = hdata.cptr_velocityX()[fromid+j];
	this->cptr_velocityY()[toid+j] = hdata.cptr_velocityY()[fromid+j];
	this->cptr_velocityZ()[toid+j] = hdata.cptr_velocityZ()[fromid+j];
      }
    }    
    if (mask & MDDataItemMask_Force){
      for (IndexType j = 0; j < numInThisCell; ++j){
	this->cptr_forceX()[toid+j] = hdata.cptr_forceX()[fromid+j];
	this->cptr_forceY()[toid+j] = hdata.cptr_forceY()[fromid+j];
	this->cptr_forceZ()[toid+j] = hdata.cptr_forceZ()[fromid+j];
      }
    }
    if (mask & MDDataItemMask_GlobalIndex){
      for (IndexType j = 0; j < numInThisCell; ++j){
	this->cptr_globalIndex()[toid+j] = hdata.cptr_globalIndex()[fromid+j];
      }
    }
    if (mask & MDDataItemMask_Mass){
      for (IndexType j = 0; j < numInThisCell; ++j){
	this->cptr_mass()[toid+j] = hdata.cptr_mass()[fromid+j];
      }
    }
    if (mask & MDDataItemMask_Type){
      for (IndexType j = 0; j < numInThisCell; ++j){
	this->cptr_type()[toid+j] = hdata.cptr_type()[fromid+j];
      }
    }
    if (mask & MDDataItemMask_Charge){
      for (IndexType j = 0; j < numInThisCell; ++j){
	this->cptr_charge()[toid+j] = hdata.cptr_charge()[fromid+j];
      }
    }
    
  }  
}


void Parallel::HostTransferPackage::
unpack_replace (HostCellListedMDData & hdata) const 
{
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();

  for (IndexType i = 0; i < numCell; ++i){
    IndexType fromid = cellStartIndex[i];
    IndexType numInThisCell = cellStartIndex[i+1] - fromid;
    if (numInThisCell == 0) continue;
    IndexType targetCellId = cellIndex[i];
    IndexType toid = targetCellId * numThreadsInCell;
    hdata.numAtomInCell[targetCellId] = numInThisCell;
    
    if (myMask & MDDataItemMask_Coordinate){
      for (IndexType j = 0; j < numInThisCell; ++j){
	hdata.cptr_coordinate()[toid+j] = this->cptr_coordinate()[fromid+j];
      }
    }
    if (myMask & MDDataItemMask_CoordinateNoi){
      for (IndexType j = 0; j < numInThisCell; ++j){
	hdata.cptr_coordinateNoiX()[toid+j] = this->cptr_coordinateNoiX()[fromid+j];
	hdata.cptr_coordinateNoiY()[toid+j] = this->cptr_coordinateNoiY()[fromid+j];
	hdata.cptr_coordinateNoiZ()[toid+j] = this->cptr_coordinateNoiZ()[fromid+j];
      }
    }
    if (myMask & MDDataItemMask_Velocity){
      for (IndexType j = 0; j < numInThisCell; ++j){
	hdata.cptr_velocityX()[toid+j] = this->cptr_velocityX()[fromid+j];
	hdata.cptr_velocityY()[toid+j] = this->cptr_velocityY()[fromid+j];
	hdata.cptr_velocityZ()[toid+j] = this->cptr_velocityZ()[fromid+j];
      }
    }    
    if (myMask & MDDataItemMask_Force){
      for (IndexType j = 0; j < numInThisCell; ++j){
	hdata.cptr_forceX()[toid+j] = this->cptr_forceX()[fromid+j];
	hdata.cptr_forceY()[toid+j] = this->cptr_forceY()[fromid+j];
	hdata.cptr_forceZ()[toid+j] = this->cptr_forceZ()[fromid+j];
      }
    }
    if (myMask & MDDataItemMask_GlobalIndex){
      for (IndexType j = 0; j < numInThisCell; ++j){
	hdata.cptr_globalIndex()[toid+j] = this->cptr_globalIndex()[fromid+j];
      }
    }
    if (myMask & MDDataItemMask_Mass){
      for (IndexType j = 0; j < numInThisCell; ++j){
	hdata.cptr_mass()[toid+j] = this->cptr_mass()[fromid+j];
      }
    }
    if (myMask & MDDataItemMask_Type){
      for (IndexType j = 0; j < numInThisCell; ++j){
	hdata.cptr_type()[toid+j] = this->cptr_type()[fromid+j];
      }
    }
    if (myMask & MDDataItemMask_Charge){
      for (IndexType j = 0; j < numInThisCell; ++j){
	hdata.cptr_charge()[toid+j] = this->cptr_charge()[fromid+j];
      }
    }
    
  }
}



void Parallel::HostTransferPackage::
unpack_add (HostCellListedMDData & hdata) const 
{
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();

  for (IndexType i = 0; i < numCell; ++i){
    IndexType fromid = cellStartIndex[i];
    IndexType numAdded = cellStartIndex[i+1] - fromid;
    if (numAdded == 0) continue;
    IndexType targetCellId = cellIndex[i];
    IndexType toid = targetCellId * numThreadsInCell + hdata.numAtomInCell[targetCellId];
    hdata.numAtomInCell[targetCellId] += numAdded;
    
    if (myMask & MDDataItemMask_Coordinate){
      for (IndexType j = 0; j < numAdded; ++j){
	hdata.cptr_coordinate()[toid+j] = this->cptr_coordinate()[fromid+j];
      }
    }
    if (myMask & MDDataItemMask_CoordinateNoi){
      for (IndexType j = 0; j < numAdded; ++j){
	hdata.cptr_coordinateNoiX()[toid+j] = this->cptr_coordinateNoiX()[fromid+j];
	hdata.cptr_coordinateNoiY()[toid+j] = this->cptr_coordinateNoiY()[fromid+j];
	hdata.cptr_coordinateNoiZ()[toid+j] = this->cptr_coordinateNoiZ()[fromid+j];
      }
    }
    if (myMask & MDDataItemMask_Velocity){
      for (IndexType j = 0; j < numAdded; ++j){
	hdata.cptr_velocityX()[toid+j] = this->cptr_velocityX()[fromid+j];
	hdata.cptr_velocityY()[toid+j] = this->cptr_velocityY()[fromid+j];
	hdata.cptr_velocityZ()[toid+j] = this->cptr_velocityZ()[fromid+j];
      }
    }    
    if (myMask & MDDataItemMask_Force){
      for (IndexType j = 0; j < numAdded; ++j){
	hdata.cptr_forceX()[toid+j] = this->cptr_forceX()[fromid+j];
	hdata.cptr_forceY()[toid+j] = this->cptr_forceY()[fromid+j];
	hdata.cptr_forceZ()[toid+j] = this->cptr_forceZ()[fromid+j];
      }
    }
    if (myMask & MDDataItemMask_GlobalIndex){
      for (IndexType j = 0; j < numAdded; ++j){
	hdata.cptr_globalIndex()[toid+j] = this->cptr_globalIndex()[fromid+j];
      }
    }
    if (myMask & MDDataItemMask_Mass){
      for (IndexType j = 0; j < numAdded; ++j){
	hdata.cptr_mass()[toid+j] = this->cptr_mass()[fromid+j];
      }
    }
    if (myMask & MDDataItemMask_Type){
      for (IndexType j = 0; j < numAdded; ++j){
	hdata.cptr_type()[toid+j] = this->cptr_type()[fromid+j];
      }
    }
    if (myMask & MDDataItemMask_Charge){
      for (IndexType j = 0; j < numAdded; ++j){
	hdata.cptr_charge()[toid+j] = this->cptr_charge()[fromid+j];
      }
    }
    
  }
}

inline static void
reallocBuffAndSize (IndexType memSize,
		    void *** buffs,
		    size_t ** sizes)
{
  *buffs = (void **) realloc (*buffs, memSize * sizeof(void*));
  if (*buffs == NULL){
    throw MDExcptFailedReallocOnHost ("HostSubCellList::collectBuffInfo",
				      "buffs",
				      memSize * sizeof(void*));
  }
  *sizes = (size_t*) realloc (*sizes, memSize * sizeof(size_t));
  if (*sizes == NULL){
    throw MDExcptFailedReallocOnHost ("HostSubCellList::collectBuffInfo",
				      "sizes",
				      memSize * sizeof(size_t));
  }
}


IndexType Parallel::HostSubCellList::
collectBuffInfo (const MDDataItemMask_t mask,
		 IndexType * num,
		 void *** buffs,
		 size_t ** sizes)
{
  // *num = this->size();
  *num = 0;
  IndexType numItem = 0;
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  
  IndexType memSize = 1;
  reallocBuffAndSize (memSize, buffs, sizes);
  
  if (mask & MDDataItemMask_Coordinate){
    numItem ++;
    for (IndexType i = 0; i < this->size(); ++i){
      if (memSize == *num){
	memSize ++;
	memSize <<= 1;
	reallocBuffAndSize (memSize, buffs, sizes);
      }
      (*buffs)[*num] = (void*) &(ptr_hdata->cptr_coordinate()
			      [this->operator[](i) * numThreadsInCell]);
      (*sizes)[*num] = ptr_hdata->getNumAtomInCell()
	  [this->operator[](i)] * sizeof (HostCoordType);
      (*num) ++;
    }
  }
    
  if (mask & MDDataItemMask_CoordinateNoi){
    numItem += 3;
    for (IndexType i = 0; i < this->size(); ++i){
      if (memSize < *num + 3){
	memSize += 3;
	memSize <<= 1;
	reallocBuffAndSize (memSize, buffs, sizes);
      }
      (*buffs)[*num+0] = (void*) &(ptr_hdata->cptr_coordinateNoiX()
				[this->operator[](i) * numThreadsInCell]);
      (*buffs)[*num+1] = (void*) &(ptr_hdata->cptr_coordinateNoiY()
				[this->operator[](i) * numThreadsInCell]);
      (*buffs)[*num+2] = (void*) &(ptr_hdata->cptr_coordinateNoiZ()
				[this->operator[](i) * numThreadsInCell]);
      size_t size =  ptr_hdata->getNumAtomInCell()[this->operator[](i)] *
	  sizeof (IntScalorType);
      (*sizes)[*num+0] = size;
      (*sizes)[*num+1] = size;
      (*sizes)[*num+2] = size;
      (*num) += 3;
    }
  }

  if (mask & MDDataItemMask_Velocity){
    numItem += 3;
    for (IndexType i = 0; i < this->size(); ++i){
      if (memSize < *num + 3){
	memSize += 3;
	memSize <<= 1;
	reallocBuffAndSize (memSize, buffs, sizes);
      }
      (*buffs)[*num+0] = (void*) &(ptr_hdata->cptr_velocityX()
				[this->operator[](i) * numThreadsInCell]);
      (*buffs)[*num+1] = (void*) &(ptr_hdata->cptr_velocityY()
				[this->operator[](i) * numThreadsInCell]);
      (*buffs)[*num+2] = (void*) &(ptr_hdata->cptr_velocityZ()
				[this->operator[](i) * numThreadsInCell]);
      size_t size =  ptr_hdata->getNumAtomInCell()[this->operator[](i)] *
	  sizeof (ScalorType);
      (*sizes)[*num+0] = size;
      (*sizes)[*num+1] = size;
      (*sizes)[*num+2] = size;
      (*num) += 3;
    }
  }

  if (mask & MDDataItemMask_Force){
    numItem += 3;
    for (IndexType i = 0; i < this->size(); ++i){
      if (memSize < *num + 3){
	memSize += 3;
	memSize <<= 1;
	reallocBuffAndSize (memSize, buffs, sizes);
      }
      (*buffs)[*num+0] = (void*) &(ptr_hdata->cptr_forceX()
				[this->operator[](i) * numThreadsInCell]);
      (*buffs)[*num+1] = (void*) &(ptr_hdata->cptr_forceY()
				[this->operator[](i) * numThreadsInCell]);
      (*buffs)[*num+2] = (void*) &(ptr_hdata->cptr_forceZ()
				[this->operator[](i) * numThreadsInCell]);
      size_t size =  ptr_hdata->getNumAtomInCell()[this->operator[](i)] *
	  sizeof (ScalorType);
      (*sizes)[*num+0] = size;
      (*sizes)[*num+1] = size;
      (*sizes)[*num+2] = size;
      (*num) += 3;
    }
  }

  if (mask & MDDataItemMask_GlobalIndex){
    numItem ++;
    for (IndexType i = 0; i < this->size(); ++i){
      if (memSize == *num){
	memSize ++;
	memSize <<= 1;
	reallocBuffAndSize (memSize, buffs, sizes);
      }
      (*buffs)[*num] = (void*) &(ptr_hdata->cptr_globalIndex()
			      [this->operator[](i) * numThreadsInCell]);
      (*sizes)[*num] =  ptr_hdata->getNumAtomInCell()
	  [this->operator[](i)] * sizeof (IndexType);
      (*num) ++;
    }
  }

  if (mask & MDDataItemMask_Type){
    numItem ++;
    for (IndexType i = 0; i < this->size(); ++i){
      if (memSize == *num){
	memSize ++;
	memSize <<= 1;
	reallocBuffAndSize (memSize, buffs, sizes);
      }
      (*buffs)[*num] = (void*) &(ptr_hdata->cptr_type()
			      [this->operator[](i) * numThreadsInCell]);
      (*sizes)[*num] =  ptr_hdata->getNumAtomInCell()
	  [this->operator[](i)] * sizeof (TypeType);
      (*num) ++;
    }
  }

  if (mask & MDDataItemMask_Mass){
    numItem ++;
    for (IndexType i = 0; i < this->size(); ++i){
      if (memSize == *num){
	memSize ++;
	memSize <<= 1;
	reallocBuffAndSize (memSize, buffs, sizes);
      }
      (*buffs)[*num] = (void*) &(ptr_hdata->cptr_mass()
			      [this->operator[](i) * numThreadsInCell]);
      (*sizes)[*num] =  ptr_hdata->getNumAtomInCell()
	  [this->operator[](i)] * sizeof (ScalorType);
      (*num) ++;
    }
  }

  if (mask & MDDataItemMask_Charge){
    numItem ++;
    for (IndexType i = 0; i < this->size(); ++i){
      if (memSize == *num){
	memSize ++;
	memSize <<= 1;
	reallocBuffAndSize (memSize, buffs, sizes);
      }
      (*buffs)[*num] = (void*) &(ptr_hdata->cptr_charge()
			      [this->operator[](i) * numThreadsInCell]);
      (*sizes)[*num] =  ptr_hdata->getNumAtomInCell()
	  [this->operator[](i)] * sizeof (ScalorType);
      (*num) ++;
    }
  }

  return memSize;
}

void Parallel::HostCellListedMDData::
clearData (const SubCellList & subList)
{
  for (IndexType i = 0; i < subList.size(); ++i){
    numAtomInCell[i] = 0;
  }
}



