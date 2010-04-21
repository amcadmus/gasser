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
HostCellListedMDData (const HostCellListedMDData & hdata)
    : memSize (0)
{
  rlist = 0;
  devideLevel = 0;
  frameLow.x = frameLow.y = frameLow.z = 0;
  frameUp.x  = frameUp.y  = frameUp.z  = 0;
  numCell.x  = numCell.y  = numCell.z  = 0;
  numAtomInCell = NULL;
  this->copy (hdata);
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
  //   throw MDExcptFailedRealloOcnHost ("HostCellListedMDData::realloc",
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
	this->cptr_coordinateNoi()[toid+j].x = hdata.cptr_coordinateNoi()[fromid+j].x;
	this->cptr_coordinateNoi()[toid+j].y = hdata.cptr_coordinateNoi()[fromid+j].y;
	this->cptr_coordinateNoi()[toid+j].z = hdata.cptr_coordinateNoi()[fromid+j].z;
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
	hdata.cptr_coordinateNoi()[toid+j].x = this->cptr_coordinateNoi()[fromid+j].x;
	hdata.cptr_coordinateNoi()[toid+j].y = this->cptr_coordinateNoi()[fromid+j].y;
	hdata.cptr_coordinateNoi()[toid+j].z = this->cptr_coordinateNoi()[fromid+j].z;
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
    if ((hdata.numAtomInCell[targetCellId] += numAdded) > numThreadsInCell){
      throw MDExcptCellList ("HostTransferPackage::unpack_add: num Atom exceed cell limit");
    }
    
    if (myMask & MDDataItemMask_Coordinate){
      for (IndexType j = 0; j < numAdded; ++j){
	hdata.cptr_coordinate()[toid+j] = this->cptr_coordinate()[fromid+j];
      }
    }
    if (myMask & MDDataItemMask_CoordinateNoi){
      for (IndexType j = 0; j < numAdded; ++j){
	hdata.cptr_coordinateNoi()[toid+j].x = this->cptr_coordinateNoi()[fromid+j].x;
	hdata.cptr_coordinateNoi()[toid+j].y = this->cptr_coordinateNoi()[fromid+j].y;
	hdata.cptr_coordinateNoi()[toid+j].z = this->cptr_coordinateNoi()[fromid+j].z;
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


// void Parallel::HostSubCellList::
// collectBuffInfo (const MDDataItemMask_t mask,
// 		 IndexType * num,
// 		 void *** buffs,
// 		 size_t ** sizes)
// {
//   // *num = this->size();
//   *num = 0;
//   IndexType numItem = 0;
//   IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  
//   IndexType memSize = 1;
//   reallocBuffAndSize (memSize, buffs, sizes);
  
//   if (mask & MDDataItemMask_Coordinate){
//     numItem ++;
//     for (IndexType i = 0; i < this->size(); ++i){
//       if (memSize == *num){
// 	memSize ++;
// 	memSize <<= 1;
// 	reallocBuffAndSize (memSize, buffs, sizes);
//       }
//       (*buffs)[*num] = (void*) &(ptr_hdata->cptr_coordinate()
// 			      [this->operator[](i) * numThreadsInCell]);
//       (*sizes)[*num] = ptr_hdata->getNumAtomInCell()
// 	  [this->operator[](i)] * sizeof (HostCoordType);
//       (*num) ++;
//     }
//   }
    
//   if (mask & MDDataItemMask_CoordinateNoi){
//     numItem += 3;
//     for (IndexType i = 0; i < this->size(); ++i){
//       if (memSize < *num + 3){
// 	memSize += 3;
// 	memSize <<= 1;
// 	reallocBuffAndSize (memSize, buffs, sizes);
//       }
//       (*buffs)[*num+0] = (void*) &(ptr_hdata->cptr_coordinateNoiX()
// 				[this->operator[](i) * numThreadsInCell]);
//       (*buffs)[*num+1] = (void*) &(ptr_hdata->cptr_coordinateNoiY()
// 				[this->operator[](i) * numThreadsInCell]);
//       (*buffs)[*num+2] = (void*) &(ptr_hdata->cptr_coordinateNoiZ()
// 				[this->operator[](i) * numThreadsInCell]);
//       size_t size =  ptr_hdata->getNumAtomInCell()[this->operator[](i)] *
// 	  sizeof (IntScalorType);
//       (*sizes)[*num+0] = size;
//       (*sizes)[*num+1] = size;
//       (*sizes)[*num+2] = size;
//       (*num) += 3;
//     }
//   }

//   if (mask & MDDataItemMask_Velocity){
//     numItem += 3;
//     for (IndexType i = 0; i < this->size(); ++i){
//       if (memSize < *num + 3){
// 	memSize += 3;
// 	memSize <<= 1;
// 	reallocBuffAndSize (memSize, buffs, sizes);
//       }
//       (*buffs)[*num+0] = (void*) &(ptr_hdata->cptr_velocityX()
// 				[this->operator[](i) * numThreadsInCell]);
//       (*buffs)[*num+1] = (void*) &(ptr_hdata->cptr_velocityY()
// 				[this->operator[](i) * numThreadsInCell]);
//       (*buffs)[*num+2] = (void*) &(ptr_hdata->cptr_velocityZ()
// 				[this->operator[](i) * numThreadsInCell]);
//       size_t size =  ptr_hdata->getNumAtomInCell()[this->operator[](i)] *
// 	  sizeof (ScalorType);
//       (*sizes)[*num+0] = size;
//       (*sizes)[*num+1] = size;
//       (*sizes)[*num+2] = size;
//       (*num) += 3;
//     }
//   }

//   if (mask & MDDataItemMask_Force){
//     numItem += 3;
//     for (IndexType i = 0; i < this->size(); ++i){
//       if (memSize < *num + 3){
// 	memSize += 3;
// 	memSize <<= 1;
// 	reallocBuffAndSize (memSize, buffs, sizes);
//       }
//       (*buffs)[*num+0] = (void*) &(ptr_hdata->cptr_forceX()
// 				[this->operator[](i) * numThreadsInCell]);
//       (*buffs)[*num+1] = (void*) &(ptr_hdata->cptr_forceY()
// 				[this->operator[](i) * numThreadsInCell]);
//       (*buffs)[*num+2] = (void*) &(ptr_hdata->cptr_forceZ()
// 				[this->operator[](i) * numThreadsInCell]);
//       size_t size =  ptr_hdata->getNumAtomInCell()[this->operator[](i)] *
// 	  sizeof (ScalorType);
//       (*sizes)[*num+0] = size;
//       (*sizes)[*num+1] = size;
//       (*sizes)[*num+2] = size;
//       (*num) += 3;
//     }
//   }

//   if (mask & MDDataItemMask_GlobalIndex){
//     numItem ++;
//     for (IndexType i = 0; i < this->size(); ++i){
//       if (memSize == *num){
// 	memSize ++;
// 	memSize <<= 1;
// 	reallocBuffAndSize (memSize, buffs, sizes);
//       }
//       (*buffs)[*num] = (void*) &(ptr_hdata->cptr_globalIndex()
// 			      [this->operator[](i) * numThreadsInCell]);
//       (*sizes)[*num] =  ptr_hdata->getNumAtomInCell()
// 	  [this->operator[](i)] * sizeof (IndexType);
//       (*num) ++;
//     }
//   }

//   if (mask & MDDataItemMask_Type){
//     numItem ++;
//     for (IndexType i = 0; i < this->size(); ++i){
//       if (memSize == *num){
// 	memSize ++;
// 	memSize <<= 1;
// 	reallocBuffAndSize (memSize, buffs, sizes);
//       }
//       (*buffs)[*num] = (void*) &(ptr_hdata->cptr_type()
// 			      [this->operator[](i) * numThreadsInCell]);
//       (*sizes)[*num] =  ptr_hdata->getNumAtomInCell()
// 	  [this->operator[](i)] * sizeof (TypeType);
//       (*num) ++;
//     }
//   }

//   if (mask & MDDataItemMask_Mass){
//     numItem ++;
//     for (IndexType i = 0; i < this->size(); ++i){
//       if (memSize == *num){
// 	memSize ++;
// 	memSize <<= 1;
// 	reallocBuffAndSize (memSize, buffs, sizes);
//       }
//       (*buffs)[*num] = (void*) &(ptr_hdata->cptr_mass()
// 			      [this->operator[](i) * numThreadsInCell]);
//       (*sizes)[*num] =  ptr_hdata->getNumAtomInCell()
// 	  [this->operator[](i)] * sizeof (ScalorType);
//       (*num) ++;
//     }
//   }

//   if (mask & MDDataItemMask_Charge){
//     numItem ++;
//     for (IndexType i = 0; i < this->size(); ++i){
//       if (memSize == *num){
// 	memSize ++;
// 	memSize <<= 1;
// 	reallocBuffAndSize (memSize, buffs, sizes);
//       }
//       (*buffs)[*num] = (void*) &(ptr_hdata->cptr_charge()
// 			      [this->operator[](i) * numThreadsInCell]);
//       (*sizes)[*num] =  ptr_hdata->getNumAtomInCell()
// 	  [this->operator[](i)] * sizeof (ScalorType);
//       (*num) ++;
//     }
//   }

//   return memSize;
// }

void Parallel::HostCellListedMDData::
clearData (const SubCellList & subList)
{
  for (IndexType i = 0; i < subList.size(); ++i){
    numAtomInCell[subList[i]] = 0;
  }
}


// IndexType Parallel::HostSubCellList::
// calculateNumAtomSent (IndexType * numSent)
// {
//   for (IndexType i = 0; i < this->size(); ++i){
//     numSent[i] = ptr_hdata->getNumAtomInCell()[this->operator[](i)];
//   }
//   return this->size();
// }


// void Parallel::HostSubCellList::
// mallocNumAtomSend (IndexType ** numSent, size_t * size)
// {
//   *size = sizeof (IndexType) * this->size();
  
//   *numSent = (IndexType *)malloc(*size);
//   if (numSent == NULL){
//     throw MDExcptFailedMallocOnHost ("HostSubCellList::calculateNumAtomSent",
// 				     "numSent", *size);
//   }
// }


Parallel::TransNumAtomInSubList::
TransNumAtomInSubList ()
    : buffs(NULL), sizes(NULL), num(0)
{
}

void Parallel::TransNumAtomInSubList::
clear ()
{
  freeAPointer ((void**)&buffs);
  freeAPointer ((void**)&sizes);
}

Parallel::TransNumAtomInSubList::
~TransNumAtomInSubList ()
{
  clear ();
}

void Parallel::TransNumAtomInSubList::
reinit (HostSubCellList & list)
{
  clear ();
  IndexType number = list.size();
  num = number;
  buffs = (void**) malloc (sizeof(void*) * number);
  if (buffs == NULL){
    throw MDExcptFailedMallocOnHost ("TransNumAtomInSubList::reinit",
				     "buffs", sizeof(void*)*number);
  }
  sizes = (size_t*) malloc (sizeof(size_t) * number);
  if (sizes == NULL){
    throw MDExcptFailedMallocOnHost ("TransNumAtomInSubList::reinit",
				     "sizes", sizeof(size_t)*number);
  }
  
  for (unsigned i = 0; i < number; ++i){
    buffs[i] = (void*) & (list.host_ptr()->cptr_numAtomInCell()[list.operator[](i)]);
    sizes[i] = sizeof(IndexType);
  }
}

void Parallel::TransNumAtomInSubList::
getTransBuffs  (IndexType * num_,
		void *** buffs_,
		size_t ** sizes_)
{
  *num_ = num;
  *buffs_ = buffs;
  *sizes_ = sizes;
}



Parallel::TransSubListData::
TransSubListData ()
    : buffs(NULL), sizes(NULL), num(0), myMask(MDDataItemMask_All)
{
}

void Parallel::TransSubListData::
clear ()
{
  freeAPointer ((void**)&buffs);
  freeAPointer ((void**)&sizes);
}

Parallel::TransSubListData::
~TransSubListData ()
{
  clear ();
}

void Parallel::TransSubListData::
reinit (HostSubCellList & list,
	const MDDataItemMask_t mask)
{
  clear ();
  IndexType number = 0;
  if (mask & MDDataItemMask_Coordinate) number+=1;
  if (mask & MDDataItemMask_CoordinateNoi) number+=1;
  if (mask & MDDataItemMask_Velocity) number+=3;
  if (mask & MDDataItemMask_Force) number+=3;
  if (mask & MDDataItemMask_GlobalIndex) number+=1;
  if (mask & MDDataItemMask_Type) number+=1;
  if (mask & MDDataItemMask_Mass) number+=1;
  if (mask & MDDataItemMask_Charge) number+=1;
  number *= list.size();
  
  buffs = (void**) malloc (sizeof(void*) * number);
  if (buffs == NULL){
    throw MDExcptFailedMallocOnHost ("TransSubListData::reinit",
				     "buffs", sizeof(void*)*number);
  }
  sizes = (size_t*) malloc (sizeof(size_t) * number);
  if (sizes == NULL){
    throw MDExcptFailedMallocOnHost ("TransSubListData::reinit",
				     "sizes", sizeof(size_t)*number);
  }

  myMask = mask;
  ptr_list = & list;
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  
  if (myMask & MDDataItemMask_Coordinate){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      (buffs)[num++] = (void*) &(ptr_list->host_ptr()->cptr_coordinate()
				 [ptr_list->operator[](i) * numThreadsInCell]);
    }
  }

  if (myMask & MDDataItemMask_CoordinateNoi){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      (buffs)[num++] = (void*) &(ptr_list->host_ptr()->cptr_coordinateNoi()
				 [ptr_list->operator[](i) * numThreadsInCell]);
    }
  }

  if (myMask & MDDataItemMask_Velocity){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      (buffs)[num++] = (void*) &(ptr_list->host_ptr()->cptr_velocityX()
				 [ptr_list->operator[](i) * numThreadsInCell]);
      (buffs)[num++] = (void*) &(ptr_list->host_ptr()->cptr_velocityY()
				 [ptr_list->operator[](i) * numThreadsInCell]);  
      (buffs)[num++] = (void*) &(ptr_list->host_ptr()->cptr_velocityZ()
				 [ptr_list->operator[](i) * numThreadsInCell]);
    }
  }

  if (myMask & MDDataItemMask_Force){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      (buffs)[num++] = (void*) &(ptr_list->host_ptr()->cptr_forceX()
				 [ptr_list->operator[](i) * numThreadsInCell]);
      (buffs)[num++] = (void*) &(ptr_list->host_ptr()->cptr_forceY()
				 [ptr_list->operator[](i) * numThreadsInCell]);  
      (buffs)[num++] = (void*) &(ptr_list->host_ptr()->cptr_forceZ()
				 [ptr_list->operator[](i) * numThreadsInCell]);
    }
  }

  if (myMask & MDDataItemMask_GlobalIndex){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      (buffs)[num++] = (void*) &(ptr_list->host_ptr()->cptr_globalIndex()
				 [ptr_list->operator[](i) * numThreadsInCell]);
    }
  }

  if (myMask & MDDataItemMask_Type){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      (buffs)[num++] = (void*) &(ptr_list->host_ptr()->cptr_type()
				 [ptr_list->operator[](i) * numThreadsInCell]);
    }
  }

  if (myMask & MDDataItemMask_Mass){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      (buffs)[num++] = (void*) &(ptr_list->host_ptr()->cptr_mass()
				 [ptr_list->operator[](i) * numThreadsInCell]);
    }
  }

  if (myMask & MDDataItemMask_Charge){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      (buffs)[num++] = (void*) &(ptr_list->host_ptr()->cptr_charge()
				 [ptr_list->operator[](i) * numThreadsInCell]);
    }
  }

  if (num != number){
    throw MDExcptCellList ("TransSubListData::, inconsistent num and number");
  }
}

void Parallel::TransSubListData::
build ()
{
  IndexType count  = 0;
  if (myMask & MDDataItemMask_Coordinate){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      (sizes)[count++] = ptr_list->host_ptr()->cptr_numAtomInCell()
	  [ptr_list->operator[](i)] * sizeof(HostCoordType);
    }
  }

  if (myMask & MDDataItemMask_CoordinateNoi){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      size_t tmp = ptr_list->host_ptr()->cptr_numAtomInCell()
	  [ptr_list->operator[](i)] * sizeof(HostCoordNoiType);
      (sizes)[count++] = tmp;
    }
  }

  if (myMask & MDDataItemMask_Velocity){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      size_t tmp = ptr_list->host_ptr()->cptr_numAtomInCell()
	  [ptr_list->operator[](i)] * sizeof(ScalorType);
      (sizes)[count++] = tmp;
      (sizes)[count++] = tmp;
      (sizes)[count++] = tmp;
    }
  }

  if (myMask & MDDataItemMask_Force){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      size_t tmp = ptr_list->host_ptr()->cptr_numAtomInCell()
	  [ptr_list->operator[](i)] * sizeof(ScalorType);
      (sizes)[count++] = tmp;
      (sizes)[count++] = tmp;
      (sizes)[count++] = tmp;
    }
  }

  if (myMask & MDDataItemMask_GlobalIndex){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      (sizes)[count++] = ptr_list->host_ptr()->cptr_numAtomInCell()
	  [ptr_list->operator[](i)] * sizeof(IndexType);
    }
  }

  if (myMask & MDDataItemMask_Type){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      (sizes)[count++] = ptr_list->host_ptr()->cptr_numAtomInCell()
	  [ptr_list->operator[](i)] * sizeof(TypeType);
    }
  }

  if (myMask & MDDataItemMask_Mass){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      (sizes)[count++] = ptr_list->host_ptr()->cptr_numAtomInCell()
	  [ptr_list->operator[](i)] * sizeof(ScalorType);
    }
  }

  if (myMask & MDDataItemMask_Charge){
    for (IndexType i = 0; i < ptr_list->size(); ++i){
      (sizes)[count++] = ptr_list->host_ptr()->cptr_numAtomInCell()
	  [ptr_list->operator[](i)] * sizeof(ScalorType);
    }
  } 
}



void Parallel::TransSubListData::
getTransBuffs  (IndexType * num_,
		void *** buffs_,
		size_t ** sizes_)
{
  *num_ = num;
  *buffs_ = buffs;
  *sizes_ = sizes;
}


void Parallel::HostCellListedMDData::
copy (const HostCellListedMDData & cellListData,
      const MDDataItemMask_t mask)
{
  HostMDData & hdata (*this);
  hdata.copy (cellListData, mask);

  rlist = cellListData.getRlist();
  devideLevel = cellListData.getDevideLevel();
  frameLow = cellListData.getFrameLow();
  frameUp  = cellListData.getFrameUp ();
  numCell  = cellListData.getNumCell ();

  if (numCell.x * numCell.y * numCell.z > memSize){
    easyReallocCell (numCell.x * numCell.y * numCell.z * MemAllocExtension);
  }
  
  for (IndexType i = 0; i < numCell.x * numCell.y * numCell.z; ++i){
    numAtomInCell[i] = cellListData.cptr_numAtomInCell()[i];
  }
}

  
void Parallel::HostCellListedMDData::
add (const HostCellListedMDData & hdata,
     const MDDataItemMask_t mask)
{
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();
  
  for (IndexType i = 0; i < totalNumCell; ++i){
    if (numAtomInCell[i] + hdata.cptr_numAtomInCell()[i] > numThreadsInCell){
      throw MDExcptCellList ("HostCellListedMDData::add: num Atom exceed cell limit");
    }

    IndexType tmp0 = numAtomInCell[i] + i * numThreadsInCell;
    IndexType tmp2 = numAtomInCell[i];
    numAtomInCell[i] += hdata.cptr_numAtomInCell()[i];
    IndexType tmp1 = numAtomInCell[i] + i * numThreadsInCell;

    if (mask & MDDataItemMask_Coordinate){
      for (IndexType j = tmp0; j < tmp1; ++j){
	this->cptr_coordinate()[j] = hdata.cptr_coordinate()[j-tmp2];
      }
    }
    if (mask & MDDataItemMask_CoordinateNoi){
      for (IndexType j = tmp0; j < tmp1; ++j){
	this->cptr_coordinateNoi()[j].x = hdata.cptr_coordinateNoi()[j-tmp2].x;
	this->cptr_coordinateNoi()[j].y = hdata.cptr_coordinateNoi()[j-tmp2].y;
	this->cptr_coordinateNoi()[j].z = hdata.cptr_coordinateNoi()[j-tmp2].z;
      }
    }
    if (mask & MDDataItemMask_Velocity){
      for (IndexType j = tmp0; j < tmp1; ++j){
	this->cptr_velocityX()[j] = hdata.cptr_velocityX()[j-tmp2];
	this->cptr_velocityY()[j] = hdata.cptr_velocityY()[j-tmp2];
	this->cptr_velocityZ()[j] = hdata.cptr_velocityZ()[j-tmp2];
      }
    }
    if (mask & MDDataItemMask_Force){
      for (IndexType j = tmp0; j < tmp1; ++j){
	this->cptr_forceX()[j] = hdata.cptr_forceX()[j-tmp2];
	this->cptr_forceY()[j] = hdata.cptr_forceY()[j-tmp2];
	this->cptr_forceZ()[j] = hdata.cptr_forceZ()[j-tmp2];
      }
    }
    if (mask & MDDataItemMask_GlobalIndex){
      for (IndexType j = tmp0; j < tmp1; ++j){
	this->cptr_globalIndex()[j] = hdata.cptr_globalIndex()[j-tmp2];
      }
    }
    if (mask & MDDataItemMask_Type){
      for (IndexType j = tmp0; j < tmp1; ++j){
	this->cptr_type()[j] = hdata.cptr_type()[j-tmp2];
      }
    }
    if (mask & MDDataItemMask_Mass){
      for (IndexType j = tmp0; j < tmp1; ++j){
	this->cptr_mass()[j] = hdata.cptr_mass()[j-tmp2];
      }
    }
    if (mask & MDDataItemMask_Charge){
      for (IndexType j = tmp0; j < tmp1; ++j){
	this->cptr_charge()[j] = hdata.cptr_charge()[j-tmp2];
      }
    }
  }
}



void Parallel::HostSubCellList::
add (const HostSubCellList & clist,
     const MDDataItemMask_t mask)
{
  IndexType numThreadsInCell = Parallel::Interface::numThreadsInCell();

  for (IndexType i = 0; i < this->size(); ++i){
    IndexType cellid = this->operator[](i);
    IndexType tmp0 = ptr_hdata->cptr_numAtomInCell()[cellid];
    IndexType tmp1 = clist.host_ptr()->cptr_numAtomInCell()[cellid];
    IndexType tmp2 = tmp0 + tmp1;
    if ((ptr_hdata->cptr_numAtomInCell()[cellid] += tmp1) > numThreadsInCell){
      throw MDExcptCellList ("HostSubCellList::add: num Atom exceed cell limit");
    }
    tmp1  = tmp0;
    tmp0 += cellid * numThreadsInCell;
    tmp2 += cellid * numThreadsInCell;
    
    if (mask & MDDataItemMask_Coordinate){
      for (IndexType j = tmp0; j < tmp2; ++j){
	ptr_hdata->cptr_coordinate()[j] = clist.ptr_hdata->cptr_coordinate()[j-tmp1];
      }
    }

    if (mask & MDDataItemMask_CoordinateNoi){
      for (IndexType j = tmp0; j < tmp2; ++j){
	ptr_hdata->cptr_coordinateNoi()[j].x =
	    clist.ptr_hdata->cptr_coordinateNoi()[j-tmp1].x;
	ptr_hdata->cptr_coordinateNoi()[j].y =
	    clist.ptr_hdata->cptr_coordinateNoi()[j-tmp1].y;
	ptr_hdata->cptr_coordinateNoi()[j].z =
	    clist.ptr_hdata->cptr_coordinateNoi()[j-tmp1].z;
      }
    }

    if (mask & MDDataItemMask_Velocity){
      for (IndexType j = tmp0; j < tmp2; ++j){
	ptr_hdata->cptr_velocityX()[j] = clist.ptr_hdata->cptr_velocityX()[j-tmp1];
	ptr_hdata->cptr_velocityY()[j] = clist.ptr_hdata->cptr_velocityY()[j-tmp1];
	ptr_hdata->cptr_velocityZ()[j] = clist.ptr_hdata->cptr_velocityZ()[j-tmp1];
      }
    }

    if (mask & MDDataItemMask_Force){
      for (IndexType j = tmp0; j < tmp2; ++j){
	ptr_hdata->cptr_forceX()[j] = clist.ptr_hdata->cptr_forceX()[j-tmp1];
	ptr_hdata->cptr_forceY()[j] = clist.ptr_hdata->cptr_forceY()[j-tmp1];
	ptr_hdata->cptr_forceZ()[j] = clist.ptr_hdata->cptr_forceZ()[j-tmp1];
      }
    }

    if (mask & MDDataItemMask_GlobalIndex){
      for (IndexType j = tmp0; j < tmp2; ++j){
	ptr_hdata->cptr_globalIndex()[j] =
	    clist.ptr_hdata->cptr_globalIndex()[j-tmp1];
      }
    }

    if (mask & MDDataItemMask_Type){
      for (IndexType j = tmp0; j < tmp2; ++j){
	ptr_hdata->cptr_type()[j] = clist.ptr_hdata->cptr_type()[j-tmp1];
      }
    }

    if (mask & MDDataItemMask_Mass){
      for (IndexType j = tmp0; j < tmp2; ++j){
	ptr_hdata->cptr_mass()[j] = clist.ptr_hdata->cptr_mass()[j-tmp1];
      }
    }

    if (mask & MDDataItemMask_Charge){
      for (IndexType j = tmp0; j < tmp2; ++j){
	ptr_hdata->cptr_charge()[j] = clist.ptr_hdata->cptr_charge()[j-tmp1];
      }
    }
  }
}


void Parallel::HostCellListedMDData::
clearData()
{
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;  
  for (IndexType i = 0; i < totalNumCell; ++i){
    numAtomInCell[i] = 0;
  }
}

void Parallel::HostCellListedMDData::
writeData_SimpleFile (const char * filename)
{
  FILE * fp = fopen (filename, "w");
  if (fp == NULL){
    throw MDExcptCannotOpenFile(filename);
  }
  IndexType totalNumCell = numCell.x * numCell.y * numCell.z;
  IndexType totalAtom = 0;
  for (IndexType i = 0; i < totalNumCell; ++i){
    totalAtom += numAtomInCell[i];
  }
  fprintf (fp, "# %d\n", totalAtom);

  for (IndexType i = 0; i < totalNumCell; ++i){
    IndexType cellShift = i * Parallel::Interface::numThreadsInCell();
    for (IndexType j = cellShift; j < cellShift + numAtomInCell[i]; ++j){
      fprintf (fp, "%8.3f %8.3f %8.3f\n", coord[j].x, coord[j].y, coord[j].z);
    }
  }
  fclose (fp);
}


Parallel::HostCellRelation::
HostCellRelation ()
    : numNeighbor (NULL),
      neighborCellIndex (NULL),
      neighborShift (NULL)
{
}

Parallel::HostCellRelation::
~HostCellRelation ()
{
  clear ();
}

void Parallel::HostCellRelation::
clear ()
{
  freeAPointer ((void**)&numNeighbor);
  freeAPointer ((void**)&neighborCellIndex);
  freeAPointer ((void**)&neighborShift);
}

void Parallel::HostCellRelation::
easyMalloc (const IndexType & totalNumCell_,
	    const IndexType & MaxNeiPerCell_)
{
  clear ();
  totalNumCell = totalNumCell_;
  MaxNeiPerCell = MaxNeiPerCell_;
  
  numNeighbor = (IndexType *) malloc (sizeof (IndexType) * totalNumCell);
  if (numNeighbor == NULL){
    throw MDExcptFailedMallocOnHost ("HostCellRelation::easyMalloc",
				     "numNeighbor", sizeof(IndexType)*totalNumCell);
  }
  neighborCellIndex = (IndexType *) malloc (sizeof(IndexType) * 
					    totalNumCell * MaxNeiPerCell);
  if (neighborCellIndex == NULL){
    throw MDExcptFailedMallocOnHost ("HostCellRelation::easyMalloc",
				     "neighborCellIndex",
				     sizeof(IndexType)*totalNumCell*MaxNeiPerCell);
  }
  neighborShift = (HostCoordType *) malloc (sizeof(HostCoordType) *
					    totalNumCell * MaxNeiPerCell);
  if (neighborShift == NULL){
    throw MDExcptFailedMallocOnHost ("HostCellRelation::easyMalloc",
				     "neighborShift",
				     sizeof(HostCoordType)*totalNumCell*MaxNeiPerCell);
  }
}

