#define HOST_CODE

#include "Parallel_MDData.h"
#include "Parallel_TransferEngine.h"
#include "Parallel_Interface.h"

#include "GromacsFileManager.h"

#include "compile_error_mixcode.h"


Parallel::HostMDData::
HostMDData()
    : numAtom_(0), memSize_(0),
      coord (NULL),
      coordNoix (NULL), coordNoiy(NULL), coordNoiz(NULL),
      velox (NULL), veloy(NULL), veloz(NULL),
      forcx (NULL), forcy(NULL), forcz(NULL),
      globalIndex (NULL),
      type (NULL), mass (NULL), charge (NULL)
{
}

Parallel::HostMDData::
~HostMDData ()
{
  clearAll ();
}

void Parallel::HostMDData::
clearAll ()
{
  numAtom_ = 0;
  memSize_ = 0;
  
  freeAPointer ((void**)&coord);
  freeAPointer ((void**)&coordNoix);
  freeAPointer ((void**)&coordNoiy);
  freeAPointer ((void**)&coordNoiz);
  freeAPointer ((void**)&velox);
  freeAPointer ((void**)&veloy);
  freeAPointer ((void**)&veloz);
  freeAPointer ((void**)&forcx);
  freeAPointer ((void**)&forcy);
  freeAPointer ((void**)&forcz);
  freeAPointer ((void**)&globalIndex);
  freeAPointer ((void**)&type);
  freeAPointer ((void**)&mass);
  freeAPointer ((void**)&charge);
}


void Parallel::HostMDData::
reallocAll (const IndexType & memSize__)
{
  reallocCoordinate	(memSize__);
  reallocVelocity	(memSize__);
  reallocForce		(memSize__);
  reallocGlobalIndex	(memSize__);
  reallocTopProperty	(memSize__);
  // reallocGroProperty	(memSize__);
}


void Parallel::HostMDData::
reallocCoordinate (const IndexType & memSize__)
{
  if (memSize_ != 0 && memSize_ != memSize__){
    throw MDExcptInconsistentMemorySizeOnHostMDData ();
  }
  if (memSize__ == 0) return ;
  memSize_ = memSize__;

  size_t sizecoord =memSize_ * sizeof(HostCoordType);
  size_t sizei = memSize_ * sizeof(IntScalorType);

  coord = (HostCoordType *) realloc (coord, sizecoord);
  if (coord == NULL) throw (MDExcptFailedReallocOnHost("coord", sizecoord));  
  coordNoix = (IntScalorType *) realloc (coordNoix, sizei);
  if (coordNoix == NULL) throw MDExcptFailedReallocOnHost ("coordNoix", sizei);
  coordNoiy = (IntScalorType *) realloc (coordNoiy, sizei);
  if (coordNoiy == NULL) throw MDExcptFailedReallocOnHost ("coordNoiy", sizei);
  coordNoiz = (IntScalorType *) realloc (coordNoiz, sizei);
  if (coordNoiz == NULL) throw MDExcptFailedReallocOnHost ("coordNoiz", sizei);
}

void Parallel::HostMDData::
reallocVelocity (const IndexType & memSize__)
{
  if (memSize_ != 0 && memSize_ != memSize__){
    throw MDExcptInconsistentMemorySizeOnHostMDData ();
  }
  if (memSize__ == 0) return ;
  memSize_ = memSize__;
  
  size_t sizef = memSize_ * sizeof(ScalorType);

  velox = (ScalorType *) realloc (velox, sizef);
  if (velox == NULL) throw MDExcptFailedReallocOnHost("velox", sizef);
  veloy = (ScalorType *) realloc (veloy, sizef);
  if (veloy == NULL) throw MDExcptFailedReallocOnHost("veloy", sizef);
  veloz = (ScalorType *) realloc (veloz, sizef);
  if (veloz == NULL) throw MDExcptFailedReallocOnHost("veloz", sizef);
}

void Parallel::HostMDData::
reallocForce (const IndexType & memSize__)
{
  if (memSize_ != 0 && memSize_ != memSize__){
    throw MDExcptInconsistentMemorySizeOnHostMDData ();
  }
  if (memSize__ == 0) return ;
  memSize_ = memSize__;
  
  size_t sizef = memSize_ * sizeof(ScalorType);
  
  forcx = (ScalorType *) realloc (forcx, sizef);
  if (forcx == NULL) throw MDExcptFailedReallocOnHost("forcx", sizef);
  forcy = (ScalorType *) realloc (forcy, sizef);
  if (forcy == NULL) throw MDExcptFailedReallocOnHost("forcy", sizef);
  forcz = (ScalorType *) realloc (forcz, sizef);
  if (forcz == NULL) throw MDExcptFailedReallocOnHost("forcz", sizef);
}

void Parallel::HostMDData::
reallocGlobalIndex (const IndexType & memSize__)
{
  if (memSize_ != 0 && memSize_ != memSize__){
    throw MDExcptInconsistentMemorySizeOnHostMDData ();
  }
  if (memSize__ == 0) return ;
  memSize_ = memSize__;

  size_t sizeidx = memSize_ * sizeof(IndexType);

  globalIndex = (IndexType *) realloc (globalIndex, sizeidx);
  if (globalIndex == NULL) throw MDExcptFailedReallocOnHost("globalIndex", sizeidx);
}

void Parallel::HostMDData::
reallocTopProperty (const IndexType & memSize__)
{
  if (memSize_ != 0 && memSize_ != memSize__){
    throw MDExcptInconsistentMemorySizeOnHostMDData ();
  }
  if (memSize__ == 0) return ;
  memSize_ = memSize__;

  size_t sizef = memSize_ * sizeof(ScalorType);

  type = (TypeType *) realloc (type, memSize_ * sizeof(TypeType));
  if (type == NULL) throw MDExcptFailedReallocOnHost("type", memSize_ * sizeof(TypeType));
  mass = (ScalorType *) realloc (mass, sizef);
  if (mass == NULL) throw MDExcptFailedReallocOnHost("mass", sizef);
  charge = (ScalorType *) realloc(charge, sizef);
  if (charge == NULL) throw MDExcptFailedReallocOnHost("charge", sizef);
}


IndexType Parallel::GlobalHostMDData::
numAtomInGroFile (const char * filename)
{
  FILE * fpc = fopen (filename, "r");
  if (fpc == NULL) {
    throw MDExcptCannotOpenFile ("HostMDData::numAtomInGroFile", filename);
  }
  while (fgetc(fpc) != '\n');

  IndexType numAtom;
  if (fscanf (fpc, "%d", &(numAtom)) != 1){
    throw MDExcptWrongFileFormat ("HostMDData::numAtomInGroFile", filename);
  }

  return numAtom;
}

void Parallel::HostMDData::
pushBackAtom  (const HostCoordType & coord_,
	       const IntScalorType & coordNoix_,
	       const IntScalorType & coordNoiy_,
	       const IntScalorType & coordNoiz_,
	       const ScalorType & velox_,
	       const ScalorType & veloy_,
	       const ScalorType & veloz_,
	       const IndexType & globalIndex_,
	       const TypeType & type_,
	       const ScalorType & mass_,
	       const ScalorType & charge_)
{
  if (numAtom_ == memSize_){
    memSize_ ++;
    memSize_ <<= 1;
    reallocAll (memSize_);
  }
  coord[numAtom_] = coord_;
  coordNoix[numAtom_] = coordNoix_;
  coordNoiy[numAtom_] = coordNoiy_;
  coordNoiz[numAtom_] = coordNoiz_;
  velox[numAtom_] = velox_;
  veloy[numAtom_] = veloy_;
  veloz[numAtom_] = veloz_;
  globalIndex[numAtom_] = globalIndex_;
  type[numAtom_] = type_;
  mass[numAtom_] = mass_;
  charge[numAtom_] = charge_;
  numAtom_ ++;
}



void Parallel::GlobalHostMDData::
initConf_GroFile (const char * filename,
		  char * atomName, IndexType * atomIndex,
		  char * resdName, IndexType * resdIndex)
{
  FILE * fpc = fopen (filename, "r");
  if (fpc == NULL) {
    throw MDExcptCannotOpenFile ("HostMDData::initConf_GroFile", filename);
  }
  while (fgetc(fpc) != '\n');

  if (fscanf (fpc, "%d", &(numAtom_)) != 1){
    throw MDExcptWrongFileFormat ("HostMDData::initConf_GroFile", filename);
  }

  if (numAtom_ > memSize_) {
    throw MDExcptNumAtomMoreThanMemSize ();
  }
  

  ScalorType bx, by, bz;
  ScalorType * tmpx, * tmpy, * tmpz;
  tmpx = (ScalorType *)malloc (sizeof(ScalorType) * memSize_);
  if (tmpx == NULL){
    throw MDExcptFailedMallocOnHost ("HostMDData::initConf_GroFile",
				     "tmpx",
				     sizeof(ScalorType) * memSize_);
  }
  tmpy = (ScalorType *)malloc (sizeof(ScalorType) * memSize_);
  if (tmpy == NULL){
    throw MDExcptFailedMallocOnHost ("HostMDData::initConf_GroFile",
				     "tmpy",
				     sizeof(ScalorType) * memSize_);
  }
  tmpz = (ScalorType *)malloc (sizeof(ScalorType) * memSize_);
  if (tmpz == NULL){
    throw MDExcptFailedMallocOnHost ("HostMDData::initConf_GroFile",
				     "tmpz",
				     sizeof(ScalorType) * memSize_);
  }
  GromacsFileManager::readGroFile (filename,
				   resdIndex, resdName, 
				   atomName, atomIndex,
				   tmpx, tmpy, tmpz,
				   velox,  veloy,  veloz,
				   &bx, &by, &bz) ;
  for (IndexType i = 0; i < numAtom_; ++i){
    globalIndex[i] = i;
    coord[i].x = tmpx[i];
    coord[i].y = tmpy[i];
    coord[i].z = tmpz[i];
    coordNoix[i] = 0;
    coordNoiy[i] = 0;
    coordNoiz[i] = 0;
  }
  free (tmpx);
  free (tmpy);
  free (tmpz);

  RectangularBoxGeometry::setBoxSize (bx, by, bz, &globalBox);

  fclose (fpc);  
}


IndexType Parallel::GlobalHostMDData::
findMolIndex (const Topology::System & sysTop,
	      const IndexType & globalIndex)
{
  const std::vector<IndexType > & indexShift (sysTop.indexShift);
  if (indexShift.size() < 2){
    throw MDExcptInvalidTopology ();
  }

  IndexType min = 0;
  IndexType max = indexShift.size() - 1;
  while (max - min > 1){
    IndexType mid = (max + min) >> 1;
    if (globalIndex >= indexShift[mid]){
      min = mid;
    }
    else {
      max = mid;
    }
  }
  return min;
}


void Parallel::GlobalHostMDData::
initTopology (const Topology::System & sysTop)
{
  for (IndexType i = 0; i < numAtom_; ++i){
    IndexType molIndex = findMolIndex (sysTop, globalIndex[i]);
    IndexType atomIndex = (globalIndex[i] - sysTop.indexShift[molIndex]) %
	(sysTop.molecules[molIndex].size());
    mass[i] = sysTop.molecules[molIndex].atoms[atomIndex].mass;
    charge[i] = sysTop.molecules[molIndex].atoms[atomIndex].charge;
    type[i] = sysTop.molecules[molIndex].atoms[atomIndex].type;
  }
}


void Parallel::HostMDData::
writeData_SimpleFile (const char * filename)
{
  FILE * fp = fopen (filename, "w");
  if (fp == NULL){
    throw MDExcptCannotOpenFile(filename);
  }
  fprintf (fp, "# %d\n", numAtom_);
  for (IndexType i = 0; i < numAtom_; ++i){
    fprintf (fp, "%8.3f %8.3f %8.3f\n", coord[i].x, coord[i].y, coord[i].z);
  }
  fclose (fp);
}


void Parallel::distributeGlobalMDData (const GlobalHostMDData & gdata,
				       HostMDData & ldata)
{
  Parallel::TransferEngine transSend ;
  Parallel::TransferEngine transRecv ;
  
  int Nx, Ny, Nz;
  Parallel::Interface::numProcDim (Nx, Ny, Nz);
  ScalorType hx(0), hy(0), hz(0);

  HostMDData sdata;
  int myRank = Parallel::Interface::myRank();

  int nx, ny, nz;
  Parallel::Interface::numProcDim(nx, ny, nz);
  
  if (myRank == 0){
    IndexType naiveLocalMemSize  = (gdata.numAtom() /
				    Parallel::Interface::numProc() * 2);
    if (naiveLocalMemSize < 100)  naiveLocalMemSize = 100 ;
    sdata.reallocAll (naiveLocalMemSize);
    sdata.setGlobalBox (gdata.getGlobalBox());
    
    hx = (gdata.getGlobalBox().size.x) / Nx;
    hy = (gdata.getGlobalBox().size.y) / Ny;
    hz = (gdata.getGlobalBox().size.z) / Nz;
  }
  
  for (int ix = 0; ix < Nx; ++ix){
    for (int iy = 0; iy < Ny; ++iy){
      for (int iz = 0; iz < Nz; ++iz){

	if (myRank == 0){
	  sdata.clearData();
	  ScalorType lx = ix * hx;
	  ScalorType ux = lx + hx;
	  ScalorType ly = iy * hy;
	  ScalorType uy = ly + hy;
	  ScalorType lz = iz * hz;
	  ScalorType uz = lz + hz;
	  for (IndexType i = 0; i < gdata.numAtom(); ++i){
	    if (gdata.cptr_coordinate()[i].x >= lx &&
		gdata.cptr_coordinate()[i].x <  ux &&
		gdata.cptr_coordinate()[i].y >= ly &&
		gdata.cptr_coordinate()[i].y <  uy &&
		gdata.cptr_coordinate()[i].z >= lz &&
		gdata.cptr_coordinate()[i].z <  uz){
	      sdata.pushBackAtom (
		  gdata.cptr_coordinate()[i],
		  gdata.cptr_coordinateNoiX()[i],
		  gdata.cptr_coordinateNoiY()[i],
		  gdata.cptr_coordinateNoiZ()[i],
		  gdata.cptr_velocityX()[i],
		  gdata.cptr_velocityY()[i],
		  gdata.cptr_velocityZ()[i],
		  gdata.cptr_globalIndex()[i],
		  gdata.cptr_type()[i],
		  gdata.cptr_mass()[i],
		  gdata.cptr_charge()[i]);
	    }
	  }
	}

	int workingRank;
	Parallel::Interface::cartCoordToRank (ix, iy, iz, workingRank);
	int recvNumAtom, recvMemSize;
	int sendNumAtom, sendMemSize;
	
	if (myRank == 0){
	  sendNumAtom = sdata.numAtom();
	  sendMemSize = sdata.memSize();
	  transSend.clearRegistered ();
	  transSend.registerBuff (&sendNumAtom, sizeof(IndexType));
	  transSend.registerBuff (&sendMemSize, sizeof(IndexType));
	  transSend.build ();
	  transSend.Isend (workingRank, 0);
	  // transSend.wait ();
	}
	if (myRank == workingRank){
	  transRecv.clearRegistered();
	  transRecv.registerBuff (&recvNumAtom, sizeof(IndexType));
	  transRecv.registerBuff (&recvMemSize, sizeof(IndexType));
	  transRecv.build ();
	  transRecv.Irecv (0, 0);
	  transRecv.wait ();
	  ldata.reallocAll (recvMemSize);
	}
	if (myRank == 0){
	  transSend.wait ();
	}
	if (myRank == 0){
	  transSend.clearRegistered ();
	  transSend.registerBuff (sdata.cptr_coordinate(),
				  sizeof(HostCoordType) * sdata.numAtom());
	  transSend.registerBuff (sdata.cptr_coordinateNoiX(),
				  sizeof(IntScalorType) * sdata.numAtom());
	  transSend.registerBuff (sdata.cptr_coordinateNoiY(),
				  sizeof(IntScalorType) * sdata.numAtom());
	  transSend.registerBuff (sdata.cptr_coordinateNoiZ(),
				  sizeof(IntScalorType) * sdata.numAtom());
	  transSend.registerBuff (sdata.cptr_velocityX(),
				  sizeof(ScalorType) * sdata.numAtom());
	  transSend.registerBuff (sdata.cptr_velocityY(),
				  sizeof(ScalorType) * sdata.numAtom());
	  transSend.registerBuff (sdata.cptr_velocityZ(),
				  sizeof(ScalorType) * sdata.numAtom());
	  transSend.registerBuff (sdata.cptr_globalIndex(),
				  sizeof(IndexType) * sdata.numAtom());
	  transSend.registerBuff (sdata.cptr_type(),
				  sizeof(TypeType) * sdata.numAtom());
	  transSend.registerBuff (sdata.cptr_mass(),
				  sizeof(ScalorType) * sdata.numAtom());
	  transSend.registerBuff (sdata.cptr_charge(),
				  sizeof(ScalorType) * sdata.numAtom());
	  ScalorType tmpx, tmpy, tmpz;
	  tmpx = sdata.getGlobalBox().size.x;
	  tmpy = sdata.getGlobalBox().size.y;
	  tmpz = sdata.getGlobalBox().size.z;
	  transSend.registerBuff (&tmpx,
				  sizeof(ScalorType) * 1);
	  transSend.registerBuff (&tmpy,
				  sizeof(ScalorType) * 1);
	  transSend.registerBuff (&tmpz,
				  sizeof(ScalorType) * 1);
	  transSend.build ();
	  transSend.Isend (workingRank, 0);
	}
	if (myRank == workingRank){
	  transRecv.clearRegistered();
	  transRecv.registerBuff (ldata.cptr_coordinate(),
				  sizeof(HostCoordType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_coordinateNoiX(),
				  sizeof(IntScalorType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_coordinateNoiY(),
				  sizeof(IntScalorType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_coordinateNoiZ(),
				  sizeof(IntScalorType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_velocityX(),
				  sizeof(ScalorType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_velocityY(),
				  sizeof(ScalorType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_velocityZ(),
				  sizeof(ScalorType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_globalIndex(),
				  sizeof(IndexType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_type(),
				  sizeof(TypeType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_mass(),
				  sizeof(ScalorType) * recvNumAtom);
	  transRecv.registerBuff (ldata.cptr_charge(),
				  sizeof(ScalorType) * recvNumAtom);
	  ScalorType tmpx, tmpy, tmpz;
	  transRecv.registerBuff (&tmpx,
				  sizeof(ScalorType) * 1);
	  transRecv.registerBuff (&tmpy,
				  sizeof(ScalorType) * 1);
	  transRecv.registerBuff (&tmpz,
				  sizeof(ScalorType) * 1);
	  transRecv.build ();
	  transRecv.Irecv (0, 0);
	  transRecv.wait ();
	  ldata.setGlobalBox (tmpx, tmpy, tmpz);
	  ldata.numAtom () = recvNumAtom;
	}
	if (myRank == 0){
	  transSend.wait ();
	}
      }
    }
  }

  (Parallel::Interface::barrier());
}




void Parallel::HostMDData::
formDataTransferBlock (const IndexType & startIndex,
		       const IndexType & num,
		       DataTransferBlock & block)
{
  block.pointer[0] = &coord[startIndex];
  block.pointer[1] = &coordNoix[startIndex];
  block.pointer[2] = &coordNoiy[startIndex];
  block.pointer[3] = &coordNoiz[startIndex];
  block.pointer[4] = &velox[startIndex];
  block.pointer[5] = &veloy[startIndex];
  block.pointer[6] = &veloz[startIndex];
  // block.pointer[7] = &forcx[startIndex];
  // block.pointer[8] = &forcy[startIndex];
  // block.pointer[9] = &forcz[startIndex];
  block.pointer[10] = &globalIndex[startIndex];
  block.pointer[11] = &type[startIndex];
  block.pointer[12] = &mass[startIndex];
  block.pointer[13] = &charge[startIndex];
  block.size[0] = num * sizeof(HostCoordType);
  block.size[1] = block.size[2] = block.size[3]
      = num * sizeof(IntScalorType);
  block.size[4] = block.size[5] = block.size[6]
      // = block.size[7] = block.size[8] = block.size[9]
      = block.size[12] = block.size[13]
      = num * sizeof(ScalorType);
  block.size[7] = block.size[8] = block.size[9] = 0;
  block.size[10] = num * sizeof(IndexType);
  block.size[11] = num * sizeof(TypeType);
}


Parallel::HostMDData::
HostMDData (const HostMDData & hdata)
    : numAtom_(0), memSize_(0)
{
  this->copy (hdata);
}

void Parallel::HostMDData::
copy (const HostMDData & hdata)
{
  if (memSize_ < hdata.numAtom_){
    reallocAll (hdata.numAtom_);
  }
  for (IndexType i = 0 ; i < hdata.numAtom_; ++i){
    coord[i] = hdata.coord[i];
    coordNoix[i] = hdata.coordNoix[i];
    coordNoiy[i] = hdata.coordNoiy[i];
    coordNoiz[i] = hdata.coordNoiz[i];
    velox[i] = hdata.velox[i];
    veloy[i] = hdata.veloy[i];
    veloz[i] = hdata.veloz[i];
    globalIndex[i] = hdata.globalIndex[i];
    type[i] = hdata.type[i];
    mass[i] = hdata.mass[i];
    charge[i] = hdata.charge[i];
  }
  // return *this;
}

// void Parallel::MDSystem::
// redistribute ()
// {
//   Parallel::TransferEngine tr_xsend0;
//   Parallel::TransferEngine tr_xrecv0;
//   Parallel::TransferEngine tr_xsend1;
//   Parallel::TransferEngine tr_xrecv1;
//   Parallel::TransferEngine tr_ysend0;
//   Parallel::TransferEngine tr_yrecv0;
//   Parallel::TransferEngine tr_ysend1;
//   Parallel::TransferEngine tr_yrecv1;
//   Parallel::TransferEngine tr_zsend0;
//   Parallel::TransferEngine tr_zrecv0;
//   Parallel::TransferEngine tr_zsend1;
//   Parallel::TransferEngine tr_zrecv1;

//   int xsend0_nei, xrecv0_nei;
//   int xsend1_nei, xrecv1_nei;
//   int ysend0_nei, yrecv0_nei;
//   int ysend1_nei, yrecv1_nei;
//   int zsend0_nei, zrecv0_nei;
//   int zsend1_nei, zrecv1_nei;
  
//   Parallel::Environment env;
//   env.neighborProcIndex (Parallel::CoordXIndex,  1, xsend0_nei, xrecv0_nei);
//   env.neighborProcIndex (Parallel::CoordXIndex, -1, xsend1_nei, xrecv1_nei);
//   env.neighborProcIndex (Parallel::CoordYIndex,  1, ysend0_nei, yrecv0_nei);
//   env.neighborProcIndex (Parallel::CoordYIndex, -1, ysend1_nei, yrecv1_nei);
//   env.neighborProcIndex (Parallel::CoordZIndex,  1, zsend0_nei, zrecv0_nei);
//   env.neighborProcIndex (Parallel::CoordZIndex, -1, zsend1_nei, zrecv1_nei);
  
//   Parallel::HostCellList::iterator it;
//   for (it = xsend0.begin(); it != xsend0.end(); ++it){
//     tr_xsend0.registerBuff (&it.numAtomInCell(), sizeof(IndexType));
//   }
//   tr_xsend0.build();
//   tr_xsend0.Isend(xsend0_nei, 0);
//   IndexType * numAtomInCell =
//       (IndexType *) malloc (sizeof(IndexType) * xrecv0.numCell());
//   if (numAtomInCell == 0){
//     throw MDExcptFailedMallocOnHost ("HostMDData::redistribute",
// 				     "numAtomInCell",
// 				     sizeof(IndexType) * xrecv0.numCell());
//   }
//   tr_xrecv0.registerBuff (numAtomInCell, sizeof(IndexType) * xrecv0.numCell());
//   tr_xrecv0.build();
//   tr_xrecv0.Irecv(xrecv0_nei, 0);
//   tr_xrecv0.wait();
//   tr_xsend0.wait();
//   tr_xrecv0.clearRegistered();
//   tr_xsend0.clearRegistered();
//   // for (it = xsend0.begin(); it != xsend.end(); ++it){
//   //   tr_xsend0.registerBuff (
//   // 	it.atomIndexInCell(), sizeof(IndexType)*it.numAtomInCell());
// }

  
  
  
  // for (itsend = xsend1.begin(); itsend != xsend1.end(); ++itsend){
  //   tr_xsend1.registerBuff (&it.numAtomInCell(), sizeof(IndexType));
  // }
  // tr_xrecv0.build();
  // tr_xsend0.Isend(xsend1_nei, 0);
  // for (itsend = ysend0.begin(); itsend != ysend0.end(); ++itsend){
  //   tr_ysend0.registerBuff (&it.numAtomInCell(), sizeof(IndexType));
  // }
  // tr_ysend0.build();
  // tr_ysend0.Isend(ysend0_nei, 0);
  // for (itsend = ysend1.begin(); itsend != ysend1.end(); ++itsend){
  //   tr_ysend1.registerBuff (&it.numAtomInCell(), sizeof(IndexType));
  // }
  // tr_yrecv0.build();
  // tr_ysend1.Isend(ysend1_nei, 0);
  // for (itsend = zsend0.begin(); itsend != zsend0.end(); ++itsend){
  //   tr_zsend0.registerBuff (&it.numAtomInCell(), sizeof(IndexType));
  // }
  // tr_zsend0.build();
  // tr_zsend0.Isend(zsend0_nei, 0);
  // for (itsend = zsend1.begin(); itsend != zsend1.end(); ++itsend){
  //   tr_zsend1.registerBuff (&it.numAtomInCell(), sizeof(IndexType));
  // }
  // tr_zrecv0.build();
  // tr_zsend1.Isend(zsend1_nei, 0);
  
  
